/*   sequin2.c
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
* File Name:  sequin2.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   1/22/95
*
* $Revision: 6.201 $
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

#include "sequin.h"
#include <document.h>
#include <sequtil.h>
#include <biosrc.h>
#include <cdrgn.h>
#include <seqsub.h>
#include <tofasta.h>
#include <gather.h>
#include <subutil.h>
#include <suggslp.h>
#include <toasn3.h>
#include <toporg.h>
#include <salfiles.h>
#include <salsap.h>
#include <salign.h>
#include <edutil.h>
#include <vsm.h>
#include <accentr.h>
#include <accutils.h>
#include <pmfapi.h>
#include <explore.h>
#include <aliparse.h>
#include <algo/blast/api/twoseq_api.h>
#ifdef WIN_MOTIF
#include <netscape.h>
#endif

extern EnumFieldAssoc  orgmod_subtype_alist [];
extern EnumFieldAssoc  subsource_subtype_alist [];
extern EnumFieldAssoc  biosource_genome_simple_alist [];

static void FindFirstTitle (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  CharPtr PNTR  ttlptr;

  if (mydata == NULL) return;
  ttlptr = (CharPtr PNTR) mydata;
  if (*ttlptr != NULL) return;
  *ttlptr = SeqEntryGetTitle (sep);
}

static void FindFirstSeqDescrTitle (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  ValNodePtr PNTR  vnpptr;

  if (mydata == NULL) return;
  vnpptr = (ValNodePtr PNTR) mydata;
  if (*vnpptr != NULL) return;
  *vnpptr = SeqEntryGetSeqDescr (sep, Seq_descr_title, NULL);
}

static void FindFirstSeqEntryTitle (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  SeqEntryPtr PNTR  sepptr;

  if (mydata == NULL) return;
  sepptr = (SeqEntryPtr PNTR) mydata;
  if (*sepptr != NULL) return;
  if (SeqEntryGetSeqDescr (sep, Seq_descr_title, NULL) != NULL) {
   *sepptr = sep;
  }
}

typedef struct fastapage {
  DIALOG_MESSAGE_BLOCK
  Char         path [PATH_MAX];
  SeqEntryPtr  list;
  ValNodePtr   errmsgs;
  DoC          doc;
  GrouP        instructions;
  Boolean      is_na;
  Boolean      is_mrna;
  Boolean      parseSeqId;
  Boolean      single;
  ButtoN       import_btn;
  ButtoN       create_btn;
  ButtoN       clear_btn;
} FastaPage, PNTR FastaPagePtr;

static ParData faParFmt = {FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0};
static ColData faColFmt = {0, 0, 80, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, TRUE};

static void ResetFastaPage (FastaPagePtr fpp)

{
  SeqEntryPtr  next;
  SeqEntryPtr  sep;

  if (fpp != NULL) {
    sep = fpp->list;
    while (sep != NULL) {
      next = sep->next;
      sep->next = NULL;
      SeqEntryFree (sep);
      sep = next;
    }
    fpp->list = NULL;
    fpp->errmsgs = ValNodeFreeData (fpp->errmsgs);
  }
}

extern void MakeSearchStringFromAlist (CharPtr str, CharPtr name)

{
  Char     ch;
  CharPtr  ptr;

  StringCpy (str, "[");
  StringCat (str, name);
  StringCat (str, "=");
  ptr = str;
  ch = *ptr;
  while (*ptr != '\0') {
    *ptr = TO_LOWER (ch);
    ptr++;
    ch = *ptr;
  }
}

static Boolean LookForSearchString (CharPtr title, CharPtr str, CharPtr tmp, size_t maxsize)

{
  CharPtr  ptr;

  ptr = StringISearch (title, str);
  if (ptr != NULL) {
    StringNCpy_0 (tmp, ptr + StringLen (str), maxsize);
     ptr = StringChr (tmp, ']');
     if (ptr != NULL) {
       *ptr = '\0';
     }
    return TRUE;
  }
  return FALSE;
}

static void AddReportLine (CharPtr str, CharPtr name, CharPtr tmp)

{
  StringCat (str, name);
  StringCat (str, ": ");
  StringCat (str, tmp);
  StringCat (str, "\n");
}

static CharPtr singlewarn = "\
ERROR - You may not enter multiple segments for a single sequence submission.\
You should either clear the nucleotide and import a single FASTA record, or\
return to the Sequence Format form and choose the proper submission type.\n\n";

#define FastaFormatBufLen 1000

static void FormatFastaDoc (FastaPagePtr fpp)

{
  EnumFieldAssocPtr  ap;
  BioseqPtr          bsp;
  Boolean            hasErrors;
  CharPtr            label;
  Int4               len;
  Char               lookfor [256];
  CharPtr            measure;
  SeqEntryPtr        nsep;
  Int2               num;
  CharPtr            plural;
  CharPtr            ptr;
  SeqIdPtr           sip;
  SeqEntryPtr        sep;
  CharPtr            str;
  CharPtr            title;
  CharPtr            ttl;
  CharPtr            tmp;
  ValNodePtr         vnp;

  if (fpp != NULL) {
    str = MemNew (sizeof (char) * FastaFormatBufLen);
    tmp = MemNew (sizeof (char) * FastaFormatBufLen);
    if (str == NULL || tmp == NULL) return;
    num = 0;
    len = 0;
    hasErrors = FALSE;
    for (sep = fpp->list; sep != NULL; sep = sep->next) {
      num++;
      if (IS_Bioseq (sep)) {
        bsp = (BioseqPtr) sep->data.ptrvalue;
        if (bsp != NULL) {
          len += bsp->length;
        }
      } else if (IS_Bioseq_set (sep)) {
        nsep = FindNucSeqEntry (sep);
        if (nsep != NULL && IS_Bioseq (nsep)) {
          bsp = (BioseqPtr) nsep->data.ptrvalue;
          if (bsp != NULL) {
            len += bsp->length;
          }
        }
      }
    }
    if (num > 1) {
      plural = "s";
    } else {
      plural = "";
    }
    if (fpp->single && num > 1) {
      AppendText (fpp->doc, singlewarn, &faParFmt, &faColFmt, programFont);
      hasErrors = TRUE;
    }
    if (fpp->is_mrna) {
      label = "Message";
      measure = "nucleotides";
    } else if (fpp->is_na) {
      label = "Segment";
      measure = "nucleotides";
    } else {
      label = "Sequence";
      measure = "amino acids";
    }
    if (fpp->is_mrna) {
      sprintf (str, "%d transcript sequence%s, total length %ld %s\n",
               (int) num, plural, (long) len, measure);
    } else if (fpp->is_na) {
      sprintf (str, "%d nucleotide segment%s, total length %ld %s\n",
               (int) num, plural, (long) len, measure);
    } else {
      sprintf (str, "%d protein sequence%s, total length %ld %s\n",
               (int) num, plural, (long) len, measure);
    }
    AppendText (fpp->doc, str, &faParFmt, &faColFmt, programFont);
    AppendText (fpp->doc, "\nChoose Clear from the Edit menu to clear these sequences",
                &faParFmt, &faColFmt, programFont);
    vnp = fpp->errmsgs;
    num = 0;
    for (sep = fpp->list; sep != NULL; sep = sep->next) {
      num++;
      len = 0;
      sip = NULL;
      tmp [0] = '\0';
      if (IS_Bioseq (sep)) {
        bsp = (BioseqPtr) sep->data.ptrvalue;
        if (bsp != NULL) {
          len = bsp->length;
          sip = SeqIdFindWorst (bsp->id);
          SeqIdWrite (sip, tmp, PRINTID_REPORT, FastaFormatBufLen);
        }
      } else if (IS_Bioseq_set (sep)) {
        nsep = FindNucSeqEntry (sep);
        if (nsep != NULL && IS_Bioseq (nsep)) {
          bsp = (BioseqPtr) nsep->data.ptrvalue;
          if (bsp != NULL) {
            len = bsp->length;
            sip = SeqIdFindWorst (bsp->id);
            SeqIdWrite (sip, tmp, PRINTID_REPORT, FastaFormatBufLen);
          }
        }
      }
      sprintf (str, "\n%s %d\nLength: %ld %s\nSequence ID: %s\n", label,
               (int) num, (long) len, measure, tmp);
      ttl = NULL;
      SeqEntryExplore (sep, (Pointer) (&ttl), FindFirstTitle);
      title = StringSaveNoNull (ttl);
      if (title != NULL && (! fpp->is_na)) {
        if (LookForSearchString (title, "[gene=", tmp, FastaFormatBufLen - 1)) {
          AddReportLine (str, "Gene", tmp);
        } else {
          StringCat (str, "No gene name detected\n");
        }
        if (LookForSearchString (title, "[prot=", tmp, FastaFormatBufLen - 1)) {
          AddReportLine (str, "Protein", tmp);
        } else if (LookForSearchString (title, "[protein=", tmp, FastaFormatBufLen - 1)) {
          AddReportLine (str, "Protein", tmp);
        } else {
          StringCat (str, "No protein name detected\n");
        }
        if (LookForSearchString (title, "[gene_syn=", tmp, FastaFormatBufLen - 1)) {
          AddReportLine (str, "Gene Syn", tmp);
        } 
        if ((LookForSearchString (title, "[protein_desc=", tmp, FastaFormatBufLen - 1)) ||
	    (LookForSearchString (title, "[prot_desc=", tmp, FastaFormatBufLen - 1))) {
          AddReportLine (str, "Protein Desc", tmp);
        } 
        ptr = StringISearch (title, "[orf]");
        if (ptr != NULL) {
        StringCat (str, "ORF indicated\n");
        }
        if (LookForSearchString (title, "[comment=", tmp, FastaFormatBufLen - 1)) {
          AddReportLine (str, "Comment", tmp);
        }
      }
      if (title != NULL && fpp->is_na && (! fpp->is_mrna)) {
        if (LookForSearchString (title, "[org=", tmp, FastaFormatBufLen - 1)) {
          AddReportLine (str, "Organism", tmp);
        }
        if (LookForSearchString (title, "[organism=", tmp, FastaFormatBufLen - 1)) {
          AddReportLine (str, "Organism", tmp);
        }
        if (LookForSearchString (title, "[lineage=", tmp, FastaFormatBufLen - 1)) {
          AddReportLine (str, "Lineage", tmp);
        }
        for (ap = orgmod_subtype_alist; ap->name != NULL; ap++) {
          MakeSearchStringFromAlist (lookfor, ap->name);
          if (LookForSearchString (title, lookfor, tmp, FastaFormatBufLen - 1)) {
            AddReportLine (str, ap->name, tmp);
          }
        }
        for (ap = subsource_subtype_alist; ap->name != NULL; ap++) {
          MakeSearchStringFromAlist (lookfor, ap->name);
          if (LookForSearchString (title, lookfor, tmp, FastaFormatBufLen - 1)) {
            AddReportLine (str, ap->name, tmp);
          }
        }
        if (LookForSearchString (title, "[note=", tmp, FastaFormatBufLen - 1)) {
          AddReportLine (str, "Note", tmp);
        }
        if (LookForSearchString (title, "[comment=", tmp, FastaFormatBufLen - 1)) {
          AddReportLine (str, "Comment", tmp);
        }
        if (LookForSearchString (title, "[subsource=", tmp, FastaFormatBufLen - 1)) {
          AddReportLine (str, "Note", tmp);
        }
        if (LookForSearchString (title, "[molecule=", tmp, FastaFormatBufLen - 1)) {
          AddReportLine (str, "Molecule", tmp);
        }
        if (LookForSearchString (title, "[moltype=", tmp, FastaFormatBufLen - 1)) {
          AddReportLine (str, "MolType", tmp);
        }
        if (LookForSearchString (title, "[location=", tmp, FastaFormatBufLen - 1)) {
          AddReportLine (str, "Location", tmp);
        }
        /*
        for (ap = biosource_genome_simple_alist; ap->name != NULL; ap++) {
          MakeSearchStringFromAlist (lookfor, ap->name);
          ptr = StringISearch (lookfor, "=");
          if (ptr != NULL) {
            *ptr = '\0';
          }
          if (LookForSearchString (title, lookfor, tmp, FastaFormatBufLen - 1)) {
            AddReportLine (str, ap->name, tmp);
          }
        }
        ptr = StringISearch (title, "[dna]");
        if (ptr != NULL) {
          AddReportLine (str, "DNA", "");
        }
        ptr = StringISearch (title, "[rna]");
        if (ptr != NULL) {
          AddReportLine (str, "RNA", "");
        }
        */
      }
      if (title != NULL && fpp->is_na && fpp->is_mrna) {
        if (LookForSearchString (title, "[gene=", tmp, FastaFormatBufLen - 1)) {
          AddReportLine (str, "Gene", tmp);
        } else {
          StringCat (str, "No gene name detected\n");
        }
        if (LookForSearchString (title, "[mrna=", tmp, FastaFormatBufLen - 1)) {
          AddReportLine (str, "mRNA", tmp);
        } else if (LookForSearchString (title, "[cdna=", tmp, FastaFormatBufLen - 1)) {
          AddReportLine (str, "cDNA", tmp);
        } else {
          StringCat (str, "No mRNA name detected\n");
        }
        if (LookForSearchString (title, "[comment=", tmp, FastaFormatBufLen - 1)) {
          AddReportLine (str, "Comment", tmp);
        }
      }
      MemFree (title);
      ttl = NULL;
      SeqEntryExplore (sep, (Pointer) (&ttl), FindFirstTitle);
      title = StringSaveNoNull (ttl);
      if (title != NULL) {
        if (! fpp->is_na) {
          ExciseString (title, "[gene=", "]");
          ExciseString (title, "[gene_syn=", "]");
          ExciseString (title, "[prot=", "]");
          ExciseString (title, "[protein=", "]");
          ExciseString (title, "[prot_desc=", "]");
          ExciseString (title, "[protein_desc=", "]");
          ExciseString (title, "[orf", "]");
          ExciseString (title, "[comment", "]");
        } else if (fpp->is_mrna) {
          ExciseString (title, "[gene=", "]");
          ExciseString (title, "[mrna=", "]");
          ExciseString (title, "[cdna=", "]");
          ExciseString (title, "[comment=", "]");
        } else {
          ExciseString (title, "[org=", "]");
          ExciseString (title, "[organism=", "]");
          ExciseString (title, "[lineage=", "]");
          for (ap = orgmod_subtype_alist; ap->name != NULL; ap++) {
            MakeSearchStringFromAlist (lookfor, ap->name);
            ExciseString (title, lookfor, "]");
          }
          for (ap = subsource_subtype_alist; ap->name != NULL; ap++) {
            MakeSearchStringFromAlist (lookfor, ap->name);
            ExciseString (title, lookfor, "]");
          }
          ExciseString (title, "[note=", "]");
          ExciseString (title, "[comment=", "]");
          ExciseString (title, "[subsource=", "]");
          /*
          ExciseString (title, "[clone=", "]");
          ExciseString (title, "[strain=", "]");
          ExciseString (title, "[isolate=", "]");
          */
          ExciseString (title, "[molecule=", "]");
          ExciseString (title, "[moltype=", "]");
          ExciseString (title, "[location=", "]");
          /*
          for (ap = biosource_genome_simple_alist; ap->name != NULL; ap++) {
            MakeSearchStringFromAlist (lookfor, ap->name);
            ptr = StringISearch (lookfor, "=");
            if (ptr != NULL) {
              *ptr = '\0';
            }
            ExciseString (title, lookfor, "]");
          }
          ExciseString (title, "[dna", "]");
          ExciseString (title, "[rna", "]");
          */
        }
        TrimSpacesAroundString (title);
        if (! StringHasNoText (title)) {
          StringCat (str, "Title: ");
          StringNCat (str, title, 128);
          StringCat (str, "\n");
        } else {
          StringCat (str, "No title detected\n");
        }
      }
      MemFree (title);
      if (vnp != NULL && vnp->data.ptrvalue != NULL) {
        hasErrors = TRUE;
        StringCat (str, (CharPtr) vnp->data.ptrvalue);
        StringCat (str, "\n");
      }
      AppendText (fpp->doc, str, &faParFmt, &faColFmt, programFont);
      if (vnp != NULL) {
        vnp = vnp->next;
      }
    }
    MemFree (str);
    MemFree (tmp);
    UpdateDocument (fpp->doc, 0, 0);
    if (hasErrors) {
      Beep ();
      Beep ();
      Beep ();
    }
  }
}
   
static Boolean ImportFastaDialog (DialoG d, CharPtr filename)

{
  BioseqPtr     bsp;
  Int4          count;
  CharPtr       errormsg;
  CharPtr       extension;
  FILE          *f;
  FastaPagePtr  fpp;
  ValNodePtr    head;
  Boolean       insegset;
  Boolean       isLocalUnknownID;
  SeqEntryPtr   last;
  Char          lastchar;
  SeqEntryPtr   lastsep;
  SeqEntryPtr   nextsep;
  ObjectIdPtr   oid;
  Char          path [PATH_MAX];
  RecT          r;
  SeqEntryPtr   sep;
  SeqIdPtr      sip;
  Char          str [32];
  ValNodePtr    vnp;

  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  fpp = (FastaPagePtr) GetObjectExtra (d);
  if (fpp != NULL) {
    extension = NULL;
    if (fpp->is_mrna) {
      extension = GetAppProperty ("FastaNucExtension");
    } else if (fpp->is_na) {
      extension = GetAppProperty ("FastaNucExtension");
    } else {
      extension = GetAppProperty ("FastaProtExtension");
    }
    if (path [0] != '\0' || GetInputFileName (path, sizeof (path), extension, "TEXT")) {
      WatchCursor ();
      StringCpy (fpp->path, path);
      ObjectRect (fpp->doc, &r);
      InsetRect (&r, 4, 4);
      faColFmt.pixWidth = r.right - r.left;
      /*
      ResetFastaPage (fpp);
      */
      Reset (fpp->doc);
      Update ();
      sep = fpp->list;
      head = fpp->errmsgs;
      errormsg = NULL;
      f = FileOpen (fpp->path, "r");
      if (f != NULL) {
        count = 0;
        last = sep;
        while (last != NULL) {
          count++;
          last = last->next;
        }
        last = sep;
        while (last != NULL && last->next != NULL) {
          last = last->next;
        }
        lastsep = NULL;
        insegset = FALSE;
        nextsep = SequinFastaToSeqEntryEx (f, fpp->is_na, &errormsg, fpp->parseSeqId, &lastchar);
        while (nextsep != NULL ||
               (lastchar != (Char) EOF && lastchar != NULLB && lastchar != 255)) {
          if (nextsep != NULL) {
            count++;
            if (IS_Bioseq (nextsep) && nextsep->data.ptrvalue != NULL) {
              bsp = (BioseqPtr) nextsep->data.ptrvalue;
              isLocalUnknownID = FALSE;
              sip = bsp->id;
              if (sip != NULL && sip->choice == SEQID_LOCAL) {
                oid = (ObjectIdPtr) sip->data.ptrvalue;
                if (oid != NULL && oid->str != NULL) {
                  isLocalUnknownID = (Boolean) (StringICmp (oid->str, "unknown") == 0);
                }
              }
              if ((! fpp->parseSeqId) || isLocalUnknownID) {
                oid = ObjectIdNew ();
                if (oid != NULL) {
                  if (fpp->is_na) {
                    sprintf (str, "nuc %ld", (long) count);
                  } else {
                    sprintf (str, "prot %ld", (long) count);
                  }
                  oid->str = StringSave (str);
                  sip = ValNodeNew (NULL);
                  if (sip != NULL) {
                    sip->choice = SEQID_LOCAL;
                    sip->data.ptrvalue = (Pointer) oid;
                    bsp->id = SeqIdFree (bsp->id);
                    bsp->id = sip;
                    SeqMgrReplaceInBioseqIndex (bsp);
                  } else {
                    ObjectIdFree (oid);
                  }
                }
              }
            }
            SeqEntryPack (nextsep);
            if (sep != NULL) {
              if (insegset) {
                if (lastsep != NULL) {
                  AddSeqEntryToSeqEntry (lastsep, nextsep, TRUE);
                } else {
                  lastsep = nextsep;
                  last->next = nextsep;
                  last = nextsep;
                }
              } else {
                last->next = nextsep;
                last = nextsep;
              }
            } else {
              if (insegset && lastsep == NULL) {
                lastsep = nextsep;
                sep = nextsep;
                last = sep;
              } else {
                sep = nextsep;
                last = sep;
              }
            }
            vnp = ValNodeNew (head);
            if (head == NULL) {
              head = vnp;
            }
            if (vnp != NULL) {
              vnp->data.ptrvalue = errormsg;
              errormsg = NULL;
            }
          } else if (lastchar == '[') {
            insegset = TRUE;
            lastsep = NULL;
          } else if (lastchar == ']') {
            insegset = FALSE;
          }
          nextsep = SequinFastaToSeqEntryEx (f, fpp->is_na, &errormsg, fpp->parseSeqId, &lastchar);
        }
        FileClose (f);
        MemFree (errormsg);
      }
      fpp->list = sep;
      fpp->errmsgs = head;
      SafeHide (fpp->instructions);
      Update ();
      Disable (fpp->import_btn);
      Disable (fpp->create_btn);
      Enable (fpp->clear_btn);
      FormatFastaDoc (fpp);
      SafeShow (fpp->doc);
      ArrowCursor ();
      Update ();
      return TRUE;
    }
  }
  return FALSE;
}

static void CleanupFastaDialog (GraphiC g, VoidPtr data)

{
  FastaPagePtr  fpp;

  fpp = (FastaPagePtr) data;
  if (fpp != NULL) {
    ResetFastaPage (fpp);
  }
  MemFree (data);
}

static CharPtr  fastaNucMsg = "\
\nPlease enter information about the nucleotide \
sequence in the spaces above.  Then click on either \
'Create My FASTA File' to create a new FASTA file or \
'Import Nucleotide FASTA' to read a previously generated \
FASTA file that contains the sequence (which can be in segments).\
  The FASTA definition lines may be of the following form:\n\n\
>ID [organism=scientific name] [strain=name] [clone=name] title\n\n\
where the [ and ] brackets are actually in the text.";

static CharPtr  fastaGenMsg = "\
\nPlease enter information about the genomic \
sequence in the spaces above.  Then click on either \
'Create My FASTA file' to create a new FASTA file or \
'Import Genomic FASTA' to read a previously generated FASTA file that \
contains the sequence (which can be in segments).  The \
FASTA definition lines may be of the following form:\n\n\
>ID [organism=scientific name] [strain=name] [clone=name] title\n\n\
where the [ and ] brackets are actually in the text.";

static CharPtr  fastaMrnaMsg  = "\
\nPlease enter information about the transcript \
sequences in the spaces above.  Then click on \
'Import Transcript FASTA' to read a FASTA file that \
contains the sequence (which can be in segments).  The \
FASTA definition lines may be of the following form:\n\n\
>ID [gene=symbol] [mrna=name] title\n\n\
where the [ and ] brackets are actually in the text.";

static CharPtr  fastaProtMsg = "\
\nPlease enter information about the protein \
sequences in the spaces above.  Then click on \
'Import Protein FASTA' to read a FASTA file that \
contains the sequences.  The FASTA definition lines should \
be of the following form:\n\n\
>ID [gene=symbol] [protein=name] title\n\n\
where the [ and ] brackets are actually in the text.";

extern DialoG CreateFastaDialog (GrouP h, CharPtr title,
                                 Boolean is_na, Boolean is_mrna, CharPtr text,
                                 Boolean parseSeqId, Boolean single)

{
  FastaPagePtr  fpp;
  GrouP         g;
  GrouP         m;
  GrouP         p;
  GrouP         s;

  p = HiddenGroup (h, 1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  fpp = (FastaPagePtr) MemNew (sizeof (FastaPage));
  if (fpp != NULL) {

    SetObjectExtra (p, fpp, CleanupFastaDialog);
    fpp->dialog = (DialoG) p;
    fpp->todialog = NULL;
    fpp->fromdialog = NULL;
    fpp->importdialog = ImportFastaDialog;

    if (title != NULL && title [0] != '\0') {
      s = NormalGroup (p, 0, -2, title, systemFont, NULL);
    } else {
      s = HiddenGroup (p, 0, -2, NULL);
    }
    m = HiddenGroup (s, -1, 0, NULL);

    fpp->path [0] = '\0';
    fpp->is_na = is_na;
    fpp->is_mrna = is_mrna;
    fpp->parseSeqId = parseSeqId;
    fpp->single = single;

    g = HiddenGroup (m, 0, 0, NULL);
    fpp->instructions = MultiLinePrompt (g, text, 27 * stdCharWidth, programFont);
    fpp->doc = DocumentPanel (g, stdCharWidth * 27, stdLineHeight * 8);
    SetDocAutoAdjust (fpp->doc, FALSE);
    Hide (fpp->doc);
    AlignObjects (ALIGN_CENTER, (HANDLE) fpp->instructions,
                  (HANDLE) fpp->doc, NULL);
  }

  return (DialoG) p;
}

typedef struct phylippage {
  DIALOG_MESSAGE_BLOCK
  Uint1        format;
  Char         path [PATH_MAX];
  SeqEntryPtr  sep;
  ValNodePtr   errmsgs;
  DoC          doc;
  GrouP        instructions;
  Char         extension [10];
  /* new for alignment */
  TexT  missing;
  TexT  beginning_gap;
  TexT  middle_gap;
  TexT  end_gap;
  TexT  match;
  Int4  type;
  DialoG genbio;
} PhylipPage, PNTR PhylipPagePtr;

static ParData phParFmt = {FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0};
static ColData phColFmt = {0, 0, 80, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, TRUE};

#define PhylipFormatBufLen 1000

static void FormatPhylipDoc (PhylipPagePtr ppp)

{
  EnumFieldAssocPtr  ap;
  BioseqPtr          bsp;
  BioseqSetPtr       bssp;
  CharPtr            label;
  Int4               len;
  Char               lookfor [128];
  CharPtr            measure;
  SeqEntryPtr        nsep;
  Int2               num;
  CharPtr            plural;
  SeqIdPtr           sip;
  SeqEntryPtr        sep;
  CharPtr            str;
  CharPtr            title;
  CharPtr            ttl;
  CharPtr            tmp;
  ValNodePtr         vnp;

  if (ppp != NULL) {
    str = MemNew (sizeof (char) * PhylipFormatBufLen);
    tmp = MemNew (sizeof (char) * PhylipFormatBufLen);
    if (str == NULL || tmp == NULL) return;
    num = 0;
    len = 0;
    sep = ppp->sep;
    if (sep != NULL && IS_Bioseq_set (sep)) {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      if (bssp != NULL && (bssp->_class == 7 ||
                           (IsPopPhyEtcSet (bssp->_class)))) {
        for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
          num++;
          if (IS_Bioseq (sep)) {
            bsp = (BioseqPtr) sep->data.ptrvalue;
            if (bsp != NULL) {
              len += bsp->length;
            }
          } else if (IS_Bioseq_set (sep)) {
            nsep = FindNucSeqEntry (sep);
            if (nsep != NULL && IS_Bioseq (nsep)) {
              bsp = (BioseqPtr) nsep->data.ptrvalue;
              if (bsp != NULL) {
                len += bsp->length;
              }
            }
          }
        }
      }
    }
    if (num > 1) {
      plural = "s";
    } else {
      plural = "";
    }
    label = "Sequence";
    measure = "nucleotides";
    sprintf (str, "%d nucleotide sequence%s, total length %ld %s\n",
             (int) num, plural, (long) len, measure);
    AppendText (ppp->doc, str, &faParFmt, &faColFmt, programFont);
    AppendText (ppp->doc, "\nChoose Clear from the Edit menu to clear these sequences",
                &faParFmt, &faColFmt, programFont);
    vnp = ppp->errmsgs;
    num = 0;
    sep = ppp->sep;
    if (sep != NULL && IS_Bioseq_set (sep)) {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      if (bssp != NULL && (bssp->_class == 7 ||
                           (IsPopPhyEtcSet (bssp->_class)))) {
        for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
          num++;
          len = 0;
          sip = NULL;
          tmp [0] = '\0';
          if (IS_Bioseq (sep)) {
            bsp = (BioseqPtr) sep->data.ptrvalue;
            if (bsp != NULL) {
              len = bsp->length;
              sip = SeqIdFindWorst (bsp->id);
              SeqIdWrite (sip, tmp, PRINTID_REPORT, FastaFormatBufLen);
            }
          } else if (IS_Bioseq_set (sep)) {
            nsep = FindNucSeqEntry (sep);
            if (nsep != NULL && IS_Bioseq (nsep)) {
              bsp = (BioseqPtr) nsep->data.ptrvalue;
              if (bsp != NULL) {
                len = bsp->length;
                sip = SeqIdFindWorst (bsp->id);
                SeqIdWrite (sip, tmp, PRINTID_REPORT, FastaFormatBufLen);
              }
            }
          }
          sprintf (str, "\n%s %d\nLength: %ld %s\nSequence ID: %s\n", label,
                   (int) num, (long) len, measure, tmp);
          ttl = NULL;
          SeqEntryExplore (sep, (Pointer) (&ttl), FindFirstTitle);
          title = StringSaveNoNull (ttl);
          if (title != NULL) {
            if (LookForSearchString (title, "[org=", tmp, PhylipFormatBufLen - 1)) {
              AddReportLine (str, "Organism", tmp);
            }
            if (LookForSearchString (title, "[organism=", tmp, PhylipFormatBufLen - 1)) {
              AddReportLine (str, "Organism", tmp);
            }
            if (LookForSearchString (title, "[lineage=", tmp, PhylipFormatBufLen - 1)) {
              AddReportLine (str, "Lineage", tmp);
            }
            for (ap = orgmod_subtype_alist; ap->name != NULL; ap++) {
              MakeSearchStringFromAlist (lookfor, ap->name);
              if (LookForSearchString (title, lookfor, tmp, PhylipFormatBufLen - 1)) {
                AddReportLine (str, ap->name, tmp);
              }
            }
            for (ap = subsource_subtype_alist; ap->name != NULL; ap++) {
              MakeSearchStringFromAlist (lookfor, ap->name);
              if (LookForSearchString (title, lookfor, tmp, PhylipFormatBufLen - 1)) {
                AddReportLine (str, ap->name, tmp);
              }
            }
            if (LookForSearchString (title, "[note=", tmp, PhylipFormatBufLen - 1)) {
              AddReportLine (str, "Note", tmp);
            }
            if (LookForSearchString (title, "[subsource=", tmp, PhylipFormatBufLen - 1)) {
              AddReportLine (str, "Note", tmp);
            }
            if (LookForSearchString (title, "[molecule=", tmp, FastaFormatBufLen - 1)) {
              AddReportLine (str, "Molecule", tmp);
            }
            if (LookForSearchString (title, "[moltype=", tmp, FastaFormatBufLen - 1)) {
              AddReportLine (str, "MolType", tmp);
            }
            if (LookForSearchString (title, "[location=", tmp, FastaFormatBufLen - 1)) {
              AddReportLine (str, "Location", tmp);
            }
            /*
            for (ap = biosource_genome_simple_alist; ap->name != NULL; ap++) {
              MakeSearchStringFromAlist (lookfor, ap->name);
              ptr = StringISearch (lookfor, "=");
              if (ptr != NULL) {
                *ptr = '\0';
              }
              if (LookForSearchString (title, lookfor, tmp, PhylipFormatBufLen - 1)) {
                AddReportLine (str, ap->name, tmp);
              }
            }
            ptr = StringISearch (title, "[dna]");
            if (ptr != NULL) {
              AddReportLine (str, "DNA", "");
            }
            ptr = StringISearch (title, "[rna]");
            if (ptr != NULL) {
              AddReportLine (str, "RNA", "");
            }
            */
          }
          MemFree (title);
          ttl = NULL;
          SeqEntryExplore (sep, (Pointer) (&ttl), FindFirstTitle);
          title = StringSaveNoNull (ttl);
          if (title != NULL) {
            ExciseString (title, "[org=", "]");
            ExciseString (title, "[organism=", "]");
            ExciseString (title, "[lineage=", "]");
            for (ap = orgmod_subtype_alist; ap->name != NULL; ap++) {
              MakeSearchStringFromAlist (lookfor, ap->name);
              ExciseString (title, lookfor, "]");
            }
            for (ap = subsource_subtype_alist; ap->name != NULL; ap++) {
              MakeSearchStringFromAlist (lookfor, ap->name);
              ExciseString (title, lookfor, "]");
            }
            ExciseString (title, "[note=", "]");
            ExciseString (title, "[subsource=", "]");
            /*
            ExciseString (title, "[clone=", "]");
            ExciseString (title, "[strain=", "]");
            ExciseString (title, "[isolate=", "]");
            */
            ExciseString (title, "[molecule=", "]");
            ExciseString (title, "[moltype=", "]");
            ExciseString (title, "[location=", "]");
            /*
            for (ap = biosource_genome_simple_alist; ap->name != NULL; ap++) {
              MakeSearchStringFromAlist (lookfor, ap->name);
              ptr = StringISearch (lookfor, "=");
              if (ptr != NULL) {
                *ptr = '\0';
              }
              ExciseString (title, lookfor, "]");
            }
            ExciseString (title, "[dna", "]");
            ExciseString (title, "[rna", "]");
            */
            TrimSpacesAroundString (title);
            if (! StringHasNoText (title)) {
              StringCat (str, "Title: ");
              StringNCat (str, title, 128);
              StringCat (str, "\n");
            } else {
              StringCat (str, "No title detected\n");
            }
          }
          MemFree (title);
          if (vnp != NULL && vnp->data.ptrvalue != NULL) {
            StringCat (str, (CharPtr) vnp->data.ptrvalue);
            StringCat (str, "\n");
          }
          AppendText (ppp->doc, str, &faParFmt, &faColFmt, programFont);
          if (vnp != NULL) {
            vnp = vnp->next;
          }
        }
      }
    }
    MemFree (str);
    MemFree (tmp);
    UpdateDocument (ppp->doc, 0, 0);
  }
}

static void ResetPhylipPage (PhylipPagePtr ppp)

{
  if (ppp != NULL) {
    ppp->sep = SeqEntryFree (ppp->sep);
    ppp->errmsgs = ValNodeFreeData (ppp->errmsgs);
  }
}

static void LookForATitle (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BoolPtr  rsult;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  rsult = (BoolPtr) mydata;
  if (SeqEntryGetSeqDescr (sep, Seq_descr_title, NULL) != NULL) {
    *rsult = FALSE;
  }
}

static Boolean SeqEntryHasNoTitles (SeqEntryPtr sep)

{
  Boolean hasNoTitles = TRUE;

  SeqEntryExplore (sep, (Pointer) &hasNoTitles, LookForATitle);
  return hasNoTitles;
}

static CharPtr noOrgInTitleWarning =
"sequences have organism information in titles. " \
"It is critical to annotate the data file with organism and source information. " \
"Please quit Sequin and read the Sequin Quick Guide section on preparing the data files before proceeding.";

static CharPtr noSrcInTitleWarning =
"sequences have source information in titles. " \
"It is critical to annotate the data file with organism and source information. " \
"Please quit Sequin and read the Sequin Quick Guide section on preparing the data files before proceeding.";

static void CountTitlesWithoutOrganisms (SeqEntryPtr sep)
{
  BioseqSetPtr bssp;
  SeqEntryPtr  tmp;
  Int2         seqtitles = 0;
  Int2         seqtotals = 0;
  CharPtr      ttl;
  Char         str [256];

  /* count titles without organisms */
  if (sep != NULL && IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && (bssp->_class == 7 || (IsPopPhyEtcSet (bssp->_class)))) {
      for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
        ttl = NULL;
        SeqEntryExplore (tmp, (Pointer) (&ttl), FindFirstTitle);
        if (ttl != NULL) {
          if (bssp->_class == BioseqseqSet_class_phy_set) {
            if (StringISearch (ttl, "[org=") != NULL || StringISearch (ttl, "[organism=") != NULL) {
              seqtitles++;
            }
          } else if (StringISearch (ttl, "[") != NULL) {
            seqtitles++;
          }
        }
        seqtotals++;
      }
      if (seqtotals != seqtitles) {
        sprintf (str, "None");
        if (seqtitles > 0) {
          sprintf (str, "Only %d", (int) seqtitles);
        }
        ArrowCursor ();
        Update ();
        Beep ();
        if (bssp->_class == BioseqseqSet_class_phy_set) {
          Message (MSG_OK, "%s of %d %s", str, (int) seqtotals, noOrgInTitleWarning);
        } else {
          Message (MSG_OK, "%s of %d %s", str, (int) seqtotals, noSrcInTitleWarning);
        }
      }
    }
  }
}

static CharPtr  phylipNucMsg = "\
\nPlease enter information about the nucleotide \
sequence in the spaces above.  Then click on \
'Import Nucleotide Alignment' to read a file that \
contains the sequences.\n\n\
Beginning Gap: When some of the sequences in an \
alignment are shorter or longer than others, beginning \
gap characters are added to the beginning of the sequence \
to maintain the correct spacing.  These will not appear \
in your sequence file.\n\
Middle Gap: These characters are used to maintain the spacing \
inside an alignment.  These are not nucleotides and will \
not appear as part of your sequence file.\n\
End Gap: When some of the sequences in an alignment are shorter \
or longer than others, end gap characters are added to the end \
of the sequence to maintain the correct spacing.  These will \
not appear in your sequence file.\n\
Ambiguous/Unknown: These characters are used to represent \
indeterminate/ambiguous nucleotides.  These will appear in your \
sequence file as 'n'.\n\
Match: These characters are used to indicate positions where \
sequences are identical to the first sequence.  These will be \
replaced by the actual characters from the first sequence.";

static void SetPhylipDocInstructions (PhylipPagePtr ppp)
{
  if (ppp == NULL || ppp->doc == NULL) return;
  Reset (ppp->doc);
  AppendText (ppp->doc, phylipNucMsg, &faParFmt, &faColFmt, programFont);
  UpdateDocument (ppp->doc, 0, 0);
  Update ();
}

static Boolean ImportPhylipDialog (DialoG d, CharPtr filename)
{
  Char           path [PATH_MAX];
  PhylipPagePtr  ppp;
  SeqEntryPtr    sep;
  RecT           r;
  FILE           *fp;
  ObjMgrDataPtr  omdptop;
  ObjMgrData     omdata;
  Uint2          parenttype;
  Pointer        parentptr;
  Char           errStr [PATH_MAX + 64];
  Char           missing [15];
  Char           beginning_gap [15];
  Char           middle_gap [15];
  Char           end_gap [15];
  Char           match [15];
  BioSourcePtr   biop;
  CharPtr        no_org_err_msg;
  const CharPtr  pop_no_org_err_msg = "You must supply organism information either "
                  "in the alignment file or on the organism page.";
  const CharPtr  phy_no_org_err_msg = "You must supply organism information in the "
                  "alignment file.";

  if (d == NULL || filename == NULL) return FALSE;

  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  ppp = (PhylipPagePtr) GetObjectExtra (d);
  if (ppp == NULL) {
    return FALSE;
  }

  if (path [0] != '\0' || GetInputFileName (path, sizeof (path), ppp->extension, "TEXT")) {
    WatchCursor ();
    StringCpy (ppp->path, path);
    ObjectRect (ppp->doc, &r);
    InsetRect (&r, 4, 4);
    faColFmt.pixWidth = r.right - r.left;
    Reset (ppp->doc);
    Update ();
    ppp->sep = SeqEntryFree (ppp->sep);
    fp = FileOpen (path, "r");
    if (fp != NULL) {
      GetTitle (ppp->missing, missing, sizeof (missing) -1);
      GetTitle (ppp->beginning_gap, beginning_gap, sizeof (beginning_gap) - 1);
      GetTitle (ppp->middle_gap, middle_gap, sizeof (middle_gap) - 1);
      GetTitle (ppp->end_gap, end_gap, sizeof (end_gap) - 1);
      GetTitle (ppp->match, match, sizeof (match) - 1);
      if (ppp->type == SEQ_PKG_POPULATION) {
        biop = (BioSourcePtr) DialogToPointer (ppp->genbio);
        if (biop == NULL) {
          no_org_err_msg = pop_no_org_err_msg;
        } else {
          no_org_err_msg = NULL;
          BioSourceFree (biop);
        }
      } else {
        no_org_err_msg = phy_no_org_err_msg;
      }
      ppp->sep = SeqEntryFromAlignmentFile (fp, missing, match, 
                                            beginning_gap, middle_gap, end_gap,
                                            "ABCDGHKMRSTUVWXYabcdghkmrstuvwxy",
                                            Seq_mol_na, no_org_err_msg);
      sep = ppp->sep;
      if (sep != NULL) {
        SaveSeqEntryObjMgrData (ppp->sep, &omdptop, &omdata);
        GetSeqEntryParent (ppp->sep, &parentptr, &parenttype);
        SeqMgrLinkSeqEntry (sep, parenttype, parentptr);
        RestoreSeqEntryObjMgrData (sep, omdptop, &omdata);

        FormatPhylipDoc (ppp);
        SafeShow (ppp->doc);

        CountTitlesWithoutOrganisms (sep);
      } else {
        SendHelpScrollMessage (helpForm, "Organism and Sequences Form", "Nucleotide Page");
        SetPhylipDocInstructions (ppp);
      }
    } else {
      SetPhylipDocInstructions (ppp);
    }
  } else {
	sprintf (errStr, "ERROR: Unable to open file %s\n\n", path);
	AppendText (ppp->doc, errStr, &faParFmt, &faColFmt, programFont);
	AppendText (ppp->doc, strerror(errno), &faParFmt, &faColFmt, programFont);
	SafeShow (ppp->doc);
    Update ();
  }
  ArrowCursor ();
  Update ();
  return TRUE;
}

static void CleanupPhylipDialog (GraphiC g, VoidPtr data)

{
  PhylipPagePtr  ppp;

  ppp = (PhylipPagePtr) data;
  if (ppp != NULL) {
    ResetPhylipPage (ppp);
  }
  MemFree (data);
}


static DialoG CreatePhylipDialog (GrouP h, CharPtr title, CharPtr text,
                                  Uint1 format, CharPtr extension,
                                  Int4 type, DialoG genbio)

{
  PhylipPagePtr  ppp;
  GrouP          g;
  GrouP          m;
  GrouP          p;
  GrouP          s;
  GrouP          a;
  RecT          r;

  p = HiddenGroup (h, 1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  ppp = (PhylipPagePtr) MemNew (sizeof (PhylipPage));
  if (ppp != NULL) {

    SetObjectExtra (p, ppp, CleanupPhylipDialog);
    ppp->dialog = (DialoG) p;
    ppp->todialog = NULL;
    ppp->fromdialog = NULL;
    ppp->importdialog = ImportPhylipDialog;
    ppp->type = type;
    ppp->genbio = genbio;

    if (title != NULL && title [0] != '\0') {
      s = NormalGroup (p, 0, -2, title, systemFont, NULL);
    } else {
      s = HiddenGroup (p, 0, -2, NULL);
    }
    m = HiddenGroup (s, -1, 0, NULL);

    ppp->format = format;
    ppp->path [0] = '\0';
    StringNCpy_0 (ppp->extension, extension, sizeof (ppp->extension));

    /* new for alignment */
    a = NormalGroup (m, 4, 0, "Sequence Characters", systemFont, NULL);
    StaticPrompt (a, "Beginning Gap", 0, dialogTextHeight, systemFont, 'c');
    ppp->beginning_gap = DialogText (a, "-.Nn?", 5, NULL);
    StaticPrompt (a, "Ambiguous/Unknown", 0, dialogTextHeight, systemFont, 'c');
    ppp->missing = DialogText (a, "?Nn", 5, NULL);
    StaticPrompt (a, "Middle Gap", 0, dialogTextHeight, systemFont, 'c');
    ppp->middle_gap = DialogText (a, "-.", 5, NULL);
    StaticPrompt (a, "Match", 0, dialogTextHeight, systemFont, 'c');
    ppp->match = DialogText (a, ":", 5, NULL);
    StaticPrompt (a, "End Gap", 0, dialogTextHeight, systemFont, 'c');
    ppp->end_gap = DialogText (a, "-.Nn?", 5, NULL);
  
    g = HiddenGroup (m, 0, 0, NULL);
    ppp->doc = DocumentPanel (g, stdCharWidth * 27, stdLineHeight * 8);
    ObjectRect (ppp->doc, &r);
    InsetRect (&r, 4, 4);
    faColFmt.pixWidth = r.right - r.left;
    SetPhylipDocInstructions (ppp);
  }

  return (DialoG) p;
}

#define ORGANISM_PAGE     0
#define NUCLEOTIDE_PAGE   1
#define MRNA_PAGE         2
#define PROTEIN_PAGE      3
#define ANNOTATE_PAGE     4

static ENUM_ALIST(biomol_nucX_alist)
  {"Genomic DNA",            253},
  {"Genomic RNA",            254},
  {"Precursor RNA",            2},
  {"mRNA [cDNA]",              3},
  {"Ribosomal RNA",            4},
  {"Transfer RNA",             5},
  {"Small nuclear RNA",        6},
  {"Small cytoplasmic RNA",    7},
  {"Other-Genetic",            9},
  {"cRNA",                    11},
  {"Small nucleolar RNA",     12},
END_ENUM_ALIST

static ENUM_ALIST(biomol_nucGen_alist)
  {"Genomic DNA",            253},
  {"Genomic RNA",            254},
END_ENUM_ALIST

static ENUM_ALIST(topology_nuc_alist)
{"Linear",          TOPOLOGY_LINEAR},
{"Circular",        TOPOLOGY_CIRCULAR},
END_ENUM_ALIST

/*---------------------------------------------------------------------*/
/*                                                                     */
/* HasZeroLengthSequence () -- Checks to see if any of a submission's  */
/*                             sequences are missing (ie -- zero       */
/*                             length).                                */
/*                                                                     */
/*---------------------------------------------------------------------*/

extern Boolean HasZeroLengthSequence (ForM newForm)
{
  SequencesFormPtr  sqfp;
  FastaPagePtr      fpp;
  SeqEntryPtr       sep;
  BioseqPtr         bsp;

  /* Get the list of Bioseqs to check */

  sqfp = (SequencesFormPtr) GetObjectExtra (newForm);
  if (NULL == sqfp)
    return TRUE;

  fpp = GetObjectExtra (sqfp->dnaseq);
  sep = fpp->list;

  /* Check the list */

  while (NULL != sep) {
    if (sep->choice == 1) { 
      bsp = (BioseqPtr) sep->data.ptrvalue;
      if (bsp->length <= 0)
	return TRUE;
    }
    sep = sep->next;
  }

  /* If we made it to here, then */
  /* there were none found.      */

  return FALSE;
}

extern Boolean SequencesFormHasProteins (ForM f)

{
  FastaPagePtr      fpp;
  SequencesFormPtr  sqfp;

  sqfp = (SequencesFormPtr) GetObjectExtra (f);
  if (sqfp != NULL) {
    if (sqfp->seqPackage >= SEQ_PKG_POPULATION) return TRUE;
    fpp = GetObjectExtra (sqfp->protseq);
    if (fpp != NULL) {
      if (fpp->path [0] != '\0') {
        return TRUE;
      }
    }
  }
  return FALSE;
}

extern Boolean SequencesFormHasTooManyNucleotides (ForM f)

{
  FastaPagePtr      fpp;
  SequencesFormPtr  sqfp;

  sqfp = (SequencesFormPtr) GetObjectExtra (f);
  if (sqfp != NULL) {
    if (sqfp->seqPackage > SEQ_PKG_SINGLE) return FALSE;
    fpp = GetObjectExtra (sqfp->dnaseq);
    if (fpp != NULL) {
      if (fpp->list != NULL && fpp->list->next != NULL) {
        return TRUE;
      }
    }
  }
  return FALSE;
}

static Boolean UpdateBspMolGatherFunc (GatherContextPtr gcp)

{
  BioseqPtr         bsp;
  SequencesFormPtr  sqfp;
  UIEnum            val;
  
  if (gcp == NULL) return TRUE;
  sqfp = (SequencesFormPtr) gcp->userdata;
  if (sqfp == NULL) return TRUE;
  if (gcp->thistype == OBJ_BIOSEQ && gcp->thisitem != NULL) {
    bsp = (BioseqPtr) gcp->thisitem;
    if (bsp->mol != Seq_mol_dna && bsp->mol != Seq_mol_rna) {
      bsp->mol = (Uint1) sqfp->dnamolfrommolinfo;
    }
    if (GetEnumPopup (sqfp->topologyPopup, topology_nuc_alist, &val)) {
      bsp->topology = (Uint1) val;
    }
  }
  return TRUE;
}

extern DialoG CreateTagListDialogEx (GrouP h, Uint2 rows, Uint2 cols,
                                     Int2 spacing, Uint2Ptr types,
                                     Uint2Ptr textWidths, EnumFieldAssocPtr PNTR alists,
                                     Boolean useBar, Boolean noExtend,
                                     ToDialogFunc tofunc, FromDialogFunc fromfunc);

typedef struct fixnucorm {
  FORM_MESSAGE_BLOCK
  PopuP              modType;
  GrouP              modGrp;
  GrouP              pptGrp;
  DialoG             modifiers;
  Int2               oldModVal;
  PrompT             modifier_intro;
  PrompT             modifier_list;
  
  SequencesFormPtr   sqfp;
} FixNucForm, PNTR FixNucFormPtr;

static ENUM_ALIST(combined_subtype_alist)
  {" ",                    0},
  {"Authority",           31},
  {"Anamorph",            36},
  {"Breed",               38},
  {"Cell-line",            2},
  {"Cell-type",            3},
  {"Chromosome",           4},
  {"Clone",                5},
  {"Clone-lib",            6},
  {"Country",             29},
  {"Cultivar",             7},
  {"Dev-stage",            8},
  {"Ecotype",             34},
  {"Forma",               32},
  {"Forma-specialis",     33},
  {"Haplotype",            9},
  {"Isolate",             10},
  {"Lab-host",            11},
  {"Lineage",             26},
  {"Location",            25},
  {"Map",                 12},
  {"Molecule",            24},
  {"Old Name",            28},
  {"Organism",             1},
  {"Plasmid-name",        14},
  {"Plastid-name",        15},
  {"Segment",             30},
  {"Sex",                 16},
  {"Specific-host",       13},
  {"Specimen-voucher",    27},
  {"Strain",              17},
  {"Sub-species",         18},
  {"Synonym",             35},
  {"Teleomorph",          37},
  {"Tissue-lib",          19},
  {"Tissue-type",         20},
  {"Variety",             22},
  {"Note",                23},
END_ENUM_ALIST

static void ModifierDialogToSeqEntryPtr (DialoG d, PopuP t, Pointer data, FixNucFormPtr fnfp)

{
  EnumFieldAssocPtr  ap;
  BioseqPtr          bsp;
  size_t             len;
  SeqEntryPtr        list;
  Char               lookfor [128];
  SeqEntryPtr        nsep;
  CharPtr            str;
  CharPtr            title;
  TagListPtr         tlp;
  CharPtr            tmp;
  ValNodePtr         ttlvnp;
  UIEnum             val;
  ValNodePtr         vnp;

  tlp = (TagListPtr) GetObjectExtra (d);
  list = (SeqEntryPtr) data;
  if (tlp != NULL && tlp->vnp != NULL) {
    val = (UIEnum) fnfp->oldModVal;
    for (ap = combined_subtype_alist; ap->name != NULL && ap->value != val; ap++) continue;
    if (ap->name == NULL) return;
    if (StringICmp (ap->name, "Organism") == 0) {
      MakeSearchStringFromAlist (lookfor, "Org");
    } else {
      MakeSearchStringFromAlist (lookfor, ap->name);
    }
    for (vnp = tlp->vnp; vnp != NULL && list != NULL; vnp = vnp->next, list = list->next) {
      bsp = NULL;
      if (IS_Bioseq (list)) {
        bsp = (BioseqPtr) list->data.ptrvalue;
      } else if (IS_Bioseq_set (list)) {
        nsep = FindNucSeqEntry (list);
        if (nsep != NULL && IS_Bioseq (nsep)) {
          bsp = (BioseqPtr) nsep->data.ptrvalue;
        }
      }
      if (bsp != NULL) {
        ttlvnp = NULL;
        SeqEntryExplore (list, (Pointer) (&ttlvnp), FindFirstSeqDescrTitle);
        if (ttlvnp == NULL) {
          ttlvnp = CreateNewDescriptor (list, Seq_descr_title);
        }
        if (ttlvnp != NULL) {
          title = (CharPtr) ttlvnp->data.ptrvalue;
          ExciseString (title, lookfor, "]");
          str = ExtractTagListColumn ((CharPtr) vnp->data.ptrvalue, 1);
          TrimSpacesAroundString (str);
          if (! StringHasNoText (str)) {
            len = StringLen (title) + StringLen (str) + StringLen (lookfor) + 10;
            tmp = MemNew (len);
            if (tmp != NULL) {
              StringCpy (tmp, title);
              StringCat (tmp, " ");
              StringCat (tmp, lookfor);
              StringCat (tmp, str);
              StringCat (tmp, "]");
              ttlvnp->data.ptrvalue = MemFree (ttlvnp->data.ptrvalue);
              ttlvnp->data.ptrvalue = StringSave (tmp);
            }
            MemFree (tmp);
          }
        }
      }
    }
  }
}

typedef struct modifierinfo 
{
  CharPtr name;
  Uint1   subtype;
  CharPtr value;
  Boolean is_org;
} ModifierInfoData, PNTR ModifierInfoPtr;

static ModifierInfoPtr ModifierInfoNew (void)
{
  ModifierInfoPtr mip;
  mip = (ModifierInfoPtr) MemNew (sizeof (ModifierInfoData));
  if (mip == NULL) return NULL;
  mip->name = NULL;
  mip->value = NULL;
  mip->is_org = FALSE;
  return mip;
}

static ModifierInfoPtr ModifierInfoFree (ModifierInfoPtr mip)
{
  if (mip == NULL) return NULL;
  mip->name = MemFree (mip->name);
  mip->value = MemFree (mip->value);
  mip = MemFree (mip);
  return mip;
}

static ValNodePtr ModifierInfoListFree (ValNodePtr list)
{
  if (list == NULL) return NULL;
  ModifierInfoListFree (list->next);
  list->next = NULL;
  list->data.ptrvalue = ModifierInfoFree (list->data.ptrvalue);
  ValNodeFree (list);
  return NULL;
}

static ModifierInfoPtr ParseOneBracketedModifier (CharPtr str)
{
  CharPtr         start, stop, eq_loc;
  ModifierInfoPtr mip;
  Int4            value_len, name_len;
  
  start = StringChr (str, '[');
  stop = StringChr (str, ']');
  eq_loc = StringChr (start + 1, '=');
  
  if (start == NULL || stop == NULL || eq_loc == NULL) return NULL;
  
  mip = ModifierInfoNew();
  if (mip == NULL) return NULL;
  name_len = eq_loc - start + 1;
  value_len = stop - eq_loc + 1;
  mip->value = (CharPtr) MemNew (value_len * sizeof (Char));
  mip->name = (CharPtr) MemNew (name_len * sizeof (Char));
  if (mip->value == NULL || mip->name == NULL)
  {
  	ModifierInfoFree (mip);
  	return NULL;
  }
  StringNCpy (mip->value, eq_loc + 1, value_len - 2);
  mip->value [value_len - 1] = 0;
  StringNCpy (mip->name, start + 1, name_len - 2);
  mip->name [name_len - 1] = 0;
  
  if (StringICmp (mip->name, "organism") == 0)
  {
	mip->is_org = TRUE;
  }
  else
  {
  	mip->subtype = FindTypeForModNameText (mip->name);
  }
  return mip;
}

static ValNodePtr ParseAllBracketedModifiers (CharPtr str)
{
  CharPtr         stop, cp;
  ValNodePtr      list = NULL;
  ValNodePtr      vnp;
  ModifierInfoPtr mip;
  
  cp = str;
  stop = StringChr (cp, ']');
  while (stop != NULL)
  {
  	mip = ParseOneBracketedModifier (cp);
  	if (mip == NULL)
  	{
  	  stop = NULL;
  	}
  	else
  	{
  	  vnp = ValNodeAdd (&list);
  	  if (vnp != NULL)
  	  {
  	  	vnp->data.ptrvalue = mip;
  	  }
  	  cp = stop + 1;
  	  stop = StringChr (cp, ']');
  	}
  }
  return list;
}

static ValNodePtr BuildModifierTypeList (ValNodePtr type_list, CharPtr new_title)
{
  ValNodePtr      modifier_info_list;
  ValNodePtr      info_vnp, type_vnp;
  ModifierInfoPtr mip;
  
  modifier_info_list = ParseAllBracketedModifiers (new_title);
  for (info_vnp = modifier_info_list; info_vnp != NULL; info_vnp = info_vnp->next)
  {
    mip = (ModifierInfoPtr)info_vnp->data.ptrvalue;
    if (mip == NULL) continue;
    if (mip->is_org) continue;
  	for (type_vnp = type_list; type_vnp != NULL && type_vnp->choice != mip->subtype; type_vnp = type_vnp->next)
  	{
  	}
  	if (type_vnp == NULL)
  	{
  	  type_vnp = ValNodeNew (type_list);
  	  if (type_list == NULL) type_list = type_vnp;
  	  if (type_vnp != NULL)
  	  {
  	  	type_vnp->choice = mip->subtype;
  	  	type_vnp->data.ptrvalue = StringSave (mip->name);
  	  }
  	}
  }
  ModifierInfoListFree (modifier_info_list);
  return type_list;
}

static void SetModifierListPrompt (FixNucFormPtr fnfp, ValNodePtr mod_list)
{
  ValNodePtr vnp;
  Int4       num_modifiers = 0;
  Int4       text_len = 0;
  CharPtr    text;
  CharPtr    text_fmt = "Already have values for:";
    
  if (mod_list == NULL)
  {
  	SetTitle (fnfp->modifier_intro, "No modifiers are present.");
  	SetTitle (fnfp->modifier_list, "");
  }
  else
  {
    for (vnp = mod_list; vnp != NULL; vnp = vnp->next)
    {
      num_modifiers ++;
      text_len += StringLen (vnp->data.ptrvalue);
    }
    text_len += num_modifiers * 6;
    text = (CharPtr) MemNew (text_len * sizeof (Char));
    if (text != NULL)
    {
      text[0] = 0;
      for (vnp = mod_list; vnp != NULL; vnp = vnp->next)
      {
        if (vnp->next == NULL && num_modifiers > 1)
        {
          StringCat (text, "and ");
        }
      	StringCat (text, vnp->data.ptrvalue);
      	if (vnp->next != NULL)
      	{
      	  if (num_modifiers > 2)
      	  {
      	  	StringCat (text, ", ");
      	  }
      	  else
      	  {
      	  	StringCat (text, " ");
      	  }
      	}
      }
      SetTitle (fnfp->modifier_intro, text_fmt);
      SetTitle (fnfp->modifier_list, text);
    }
  }	
}

static void SeqEntryPtrToModifierDialog (DialoG d, PopuP p, Pointer data, FixNucFormPtr fnfp)

{
  EnumFieldAssocPtr  ap;
  BioseqPtr          bsp;
  ValNodePtr         head;
  Int2               j;
  size_t             len;
  SeqEntryPtr        list;
  Char               lookfor [128];
  SeqEntryPtr        nsep;
  CharPtr            ptr;
  SeqIdPtr           sip;
  CharPtr            str;
  Char               text [128];
  CharPtr            title;
  TagListPtr         tlp;
  CharPtr            ttl;
  Char               tmp [128];
  UIEnum             val;
  ValNodePtr         vnp;
  ValNodePtr         found_modifiers = NULL;


  tlp = (TagListPtr) GetObjectExtra (d);
  list = (SeqEntryPtr) data;

  if (tlp != NULL) {
    if (! GetEnumPopup (p, combined_subtype_alist, &val)) return;
    fnfp->oldModVal = (Int2) val;
    for (ap = combined_subtype_alist; ap->name != NULL && ap->value != val; ap++) continue;
    if (ap->name == NULL) return;
    if (StringICmp (ap->name, "Organism") == 0) {
      MakeSearchStringFromAlist (lookfor, "Org");
    } else {
      MakeSearchStringFromAlist (lookfor, ap->name);
    }
    head = NULL;
    while (list != NULL) {
      bsp = NULL;
      if (IS_Bioseq (list)) {
        bsp = (BioseqPtr) list->data.ptrvalue;
      } else if (IS_Bioseq_set (list)) {
        nsep = FindNucSeqEntry (list);
        if (nsep != NULL && IS_Bioseq (nsep)) {
          bsp = (BioseqPtr) nsep->data.ptrvalue;
        }
      }
      if (bsp != NULL) {
        vnp = ValNodeNew (head);
        if (head == NULL) {
          head = vnp;
        }
        if (vnp != NULL) {
          ttl = NULL;
          SeqEntryExplore (list, (Pointer) (&ttl), FindFirstTitle);
          title = StringSaveNoNull (ttl);
          found_modifiers = BuildModifierTypeList (found_modifiers, title);

          sip = SeqIdFindWorst (bsp->id);
          SeqIdWrite (sip, tmp, PRINTID_REPORT, sizeof (tmp));
          ptr = StringChr (tmp, '|');
          if (ptr == NULL) {
            ptr = tmp;
          } else {
            ptr++;
          }
          text [0] = '\0';
          if (! LookForSearchString (title, lookfor, text, sizeof (text) - 1)) {
            StringCpy (text, " ");
          }
          len = StringLen (ptr) + StringLen (text);
          str = MemNew (len + 4);
          if (str != NULL) {
            StringCpy (str, ptr);
            StringCat (str, "\t");
            StringCat (str, text);
            StringCat (str, "\n");
          }
          vnp->data.ptrvalue = str;
        }
      }
      list = list->next;
    }
    SetModifierListPrompt (fnfp, found_modifiers);
    found_modifiers = ValNodeFreeData (found_modifiers);
  
    SendMessageToDialog (tlp->dialog, VIB_MSG_RESET);
    tlp->vnp = head;
    SendMessageToDialog (tlp->dialog, VIB_MSG_REDRAW);
    for (j = 0, vnp = tlp->vnp; vnp != NULL; j++, vnp = vnp->next) {
    }
    tlp->max = MAX ((Int2) 0, (Int2) (j - tlp->rows));
    CorrectBarMax (tlp->bar, tlp->max);
    CorrectBarPage (tlp->bar, tlp->rows - 1, tlp->rows - 1);
    if (tlp->max > 0) {
      SafeShow (tlp->bar);
    } else {
      SafeHide (tlp->bar);
    }
    SendMessageToDialog (tlp->dialog, VIB_MSG_ENTER);
  }
}

static CharPtr noOrgInTitleReject =
"sequences have organism information in titles. " \
"It is critical to annotate the data file with organism and source information. " \
"Sequin will not continue until you supply this information.";

static Boolean NotEnoughOrgTitles (SequencesFormPtr sqfp, SeqEntryPtr sep)

{
  MsgAnswer  ans;
  Int2       seqtitles = 0;
  Int2       seqtotals = 0;
  Char       str [256];
  CharPtr    title;

  if (sqfp == NULL || sep == NULL) return FALSE;

  while (sep != NULL) {
    title = NULL;
    SeqEntryExplore (sep, (Pointer) (&title), FindFirstTitle);
    if (title != NULL) {
      if (sqfp->seqPackage == SEQ_PKG_PHYLOGENETIC) {
        if (StringISearch (title, "[org=") != NULL ||
            StringISearch (title, "[organism=") != NULL) {
          seqtitles++;
        }
      } else if (StringISearch (title, "[") != NULL) {
        seqtitles++;
      }
    }
    sep = sep->next;
    seqtotals++;
  }

  if (seqtotals != seqtitles) {
    sprintf (str, "None");
    if (seqtitles > 0) {
      sprintf (str, "Only %d", (int) seqtitles);
    }
    ArrowCursor ();
    Update ();
    Beep ();
    if (! indexerVersion) {
      if (sqfp->seqPackage == SEQ_PKG_PHYLOGENETIC) {
        Message (MSG_OK, "%s of %d %s", str, (int) seqtotals, noOrgInTitleReject);
        return TRUE;
      }
    } else {
      if (sqfp->seqPackage == SEQ_PKG_PHYLOGENETIC) {
        ans = Message (MSG_YN, "%s of %d %s (Regular version will block here - continue?)", str, (int) seqtotals, noOrgInTitleReject);
        if (ans == ANS_NO) return TRUE;
      }
    }
  }

  return FALSE;
}

static void AcceptNucleotideFixup (ButtoN b)

{
  FixNucFormPtr     fnfp;
  SeqEntryPtr       sep;
  SequencesFormPtr  sqfp;

  fnfp = (FixNucFormPtr) GetObjectExtra (b);
  if (fnfp == NULL) return;
  if (fnfp->sqfp != NULL) {
    Hide (fnfp->form);
    Update ();
    sqfp = fnfp->sqfp;
    if (fnfp->oldModVal != 0) {
      sep = sqfp->currConfirmSeq;
      ModifierDialogToSeqEntryPtr (fnfp->modifiers, fnfp->modType, sep, fnfp);
    }
    sep = sqfp->currConfirmSeq;
    if (NotEnoughOrgTitles (sqfp, sep)) {
      Show (fnfp->form);
      Update ();
      return;
    }
    Remove (fnfp->form);
    Update ();
    if (sqfp->putItAllTogether != NULL) {
      sqfp->putItAllTogether (sqfp->form);
    }
  }
}

static void ChangeModifier (PopuP p)

{
  FixNucFormPtr     fnfp;
  SeqEntryPtr       sep;
  SequencesFormPtr  sqfp;
  UIEnum            val;

  fnfp = (FixNucFormPtr) GetObjectExtra (p);
  if (fnfp == NULL) return;
  sqfp = fnfp->sqfp;
  if (sqfp == NULL) return;
  sep = sqfp->currConfirmSeq;
  if (! GetEnumPopup (p, combined_subtype_alist, &val)) {
    val = 0;
  }
  if (val > 0) {
    SafeHide (fnfp->pptGrp);
    SafeShow (fnfp->modGrp);
  } else {
    SafeHide (fnfp->modGrp);
    SafeShow (fnfp->pptGrp);
  }
  if (fnfp->oldModVal != 0) {
    ModifierDialogToSeqEntryPtr (fnfp->modifiers, fnfp->modType, sep, fnfp);
  }
  SeqEntryPtrToModifierDialog (fnfp->modifiers, fnfp->modType, sep, fnfp);
}

static void FixNucMessage (ForM f, Int2 mssg)

{
  FixNucFormPtr  fnfp;

  fnfp = (FixNucFormPtr) GetObjectExtra (f);
  if (fnfp != NULL) {
    switch (mssg) {
      case VIB_MSG_QUIT :
        break;
      case VIB_MSG_CLOSE :
        break;
      case VIB_MSG_CUT :
        StdCutTextProc (NULL);
        break;
      case VIB_MSG_COPY :
        StdCopyTextProc (NULL);
        break;
      case VIB_MSG_PASTE :
        StdPasteTextProc (NULL);
        break;
      case VIB_MSG_DELETE :
        StdDeleteTextProc (NULL);
        break;
      default :
        if (fnfp->appmessage != NULL) {
          fnfp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static void FixUpNucActivate (WindoW w)

{
  IteM           closeItm;
  IteM           exportItm;
  FixNucFormPtr  fnfp;
  IteM           importItm;

  fnfp = (FixNucFormPtr) GetObjectExtra (w);
  if (fnfp != NULL) {
    if (fnfp->activate != NULL) {
      fnfp->activate (w);
    }
    importItm = FindFormMenuItem ((BaseFormPtr) fnfp, VIB_MSG_IMPORT);
    exportItm = FindFormMenuItem ((BaseFormPtr) fnfp, VIB_MSG_EXPORT);
    closeItm = FindFormMenuItem ((BaseFormPtr) fnfp, VIB_MSG_CLOSE);
    SafeDisable (importItm);
    SafeDisable (exportItm);
    SafeDisable (closeItm);
  }
}

static Uint2 modedit_types [] = {
  TAGLIST_PROMPT, TAGLIST_TEXT
};

static Uint2 modedit_widths [] = {
  0, 0,
};

static CharPtr  sourceModMsg = "\
\nThe Modifier popup lets you select between organism name, \
strain, isolate, and other source modifiers.\n\n\
The resulting spreadsheet lets you enter or edit source \
information for all sequence components.\n\n\
You may enter data for multiple modifiers. \
The data for each modifier is saved when you use the menu \
to change to a different modifier.\n\n\
Scientific names should not be abbreviated (use 'Drosophila \
melanogaster' instead of 'D. melanogaster').";

static void LetUserFixNucleotideInfo (SequencesFormPtr sqfp)

{
  ButtoN             b;
  GrouP              c;
  FixNucFormPtr      fnfp;
  GrouP              g;
  GrouP              h;
  GrouP              k;
  SeqEntryPtr        sep;
  StdEditorProcsPtr  sepp;
  GrouP              t;
  WindoW             w;
  GrouP              x;

  if (sqfp == NULL) return;

  sep = sqfp->currConfirmSeq;

  w = NULL;
  fnfp = (FixNucFormPtr) MemNew (sizeof (FixNucForm));
  if (fnfp != NULL) {
    fnfp->sqfp = sqfp;
    w = FixedWindow (-50, -33, -10, -10, "Source Modifiers", NULL);
    SetObjectExtra (w, fnfp, StdCleanupFormProc);
    fnfp->form = (ForM) w;
    fnfp->formmessage = FixNucMessage;

#ifndef WIN_MAC
    CreateStdEditorFormMenus (w);
#endif

    fnfp->activate = NULL;
    sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
    if (sepp != NULL) {
      fnfp->activate = sepp->activateForm;
      fnfp->appmessage = sepp->handleMessages;
    }
    SetActivate (w, FixUpNucActivate);

    h = HiddenGroup (w, -1, 0, NULL);
    SetGroupSpacing (h, 10, 10);

    g = HiddenGroup (h, 2, 0, NULL);
    StaticPrompt (g, "Modifier", 0, popupMenuHeight, programFont, 'l');
    fnfp->modType = PopupList (g, TRUE, ChangeModifier);
    InitEnumPopup (fnfp->modType, combined_subtype_alist, NULL);
    SetEnumPopup (fnfp->modType, combined_subtype_alist, 0);
    SetObjectExtra (fnfp->modType, fnfp, NULL);

    k = HiddenGroup (h, 0, 0, NULL);

    fnfp->modGrp = HiddenGroup (k, 1, 0, NULL);
    t = HiddenGroup (fnfp->modGrp, 2, 0, NULL);
    StaticPrompt (t, "SeqID", 7 * stdCharWidth, 0, programFont, 'l');
    StaticPrompt (t, "Value", 18 * stdCharWidth, 0, programFont, 'l');

    modedit_widths [0] = 7;
    modedit_widths [1] = 18;
    fnfp->modifiers = CreateTagListDialogEx (fnfp->modGrp, 5, 2, 2,
                                             modedit_types, modedit_widths,
                                             NULL, TRUE, TRUE, NULL, NULL);
    fnfp->oldModVal = 0;
    Hide (fnfp->modGrp);

    fnfp->pptGrp = HiddenGroup (k, 1, 0, NULL);
    MultiLinePrompt (fnfp->pptGrp, sourceModMsg, 25 * stdCharWidth, programFont);
    
    x = HiddenGroup (h, 0, 2, NULL);
    fnfp->modifier_intro = StaticPrompt (x, "", 30 * stdCharWidth, popupMenuHeight, programFont, 'l');
    fnfp->modifier_list = StaticPrompt (x, "", 30 * stdCharWidth, popupMenuHeight, programFont, 'l');

    c = HiddenGroup (h, 2, 0, NULL);
    b = DefaultButton (c, "OK", AcceptNucleotideFixup);
    SetObjectExtra (b, fnfp, NULL);

    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) x, (HANDLE) fnfp->modGrp,
                  (HANDLE) fnfp->pptGrp, (HANDLE) c, NULL);

    RealizeWindow (w);

    SeqEntryPtrToModifierDialog (fnfp->modifiers, fnfp->modType, sep, fnfp);
    Show (w);
    Select (w);
    Update ();
    SendHelpScrollMessage (helpForm, "Source Modifiers Form", NULL);
  }
}

typedef struct fixprotform {
  FORM_MESSAGE_BLOCK
  TexT               geneDesc;
  TexT               geneSymbol;
  ButtoN             orf;
  TexT               protDesc;
  TexT               protName;
  TexT               seqID;
  TexT               comment;
  TexT               title;

  SequencesFormPtr   sqfp;
  SeqEntryPtr        esep;
  SeqEntryPtr        nsep;
  SeqEntryPtr        psep;
  BioseqPtr          pbsp;
} FixProtForm, PNTR FixProtFormPtr;

static void LetUserFixProteinInfo (SequencesFormPtr sqfp);

static void AcceptThisFixup (ButtoN b)

{
  FixProtFormPtr    fpfp;
  Boolean           needsOpener;
  SeqEntryPtr       old;
  BioseqPtr         pbsp;
  SeqEntryPtr       psep;
  CharPtr           ptr;
  SeqLocPtr         slp;
  SequencesFormPtr  sqfp;
  Char              str [256];
  CharPtr           title;
  Char              tmp [256];
  ValNodePtr        vnp;

  fpfp = (FixProtFormPtr) GetObjectExtra (b);
  if (fpfp == NULL) return;
  if (fpfp->sqfp != NULL) {
    Hide (fpfp->form);
    Update ();
    sqfp = fpfp->sqfp;
    psep = fpfp->psep;
    pbsp = fpfp->pbsp;
    GetTitle (fpfp->seqID, str, sizeof (str));
    if (! StringHasNoText (str)) {
      pbsp->id = SeqIdFree (pbsp->id);
      pbsp->id = MakeSeqID (str);
      SeqMgrReplaceInBioseqIndex (pbsp);
    } else {
      pbsp->id = SeqIdFree (pbsp->id);
      slp = CreateWholeInterval (fpfp->nsep);
      old = SeqEntrySetScope (fpfp->esep);
      pbsp->id = MakeNewProteinSeqId (slp, NULL);
      SeqMgrReplaceInBioseqIndex (pbsp);
      SeqEntrySetScope (old);
      SeqLocFree (slp);
    }
    title = MemNew (1000 * sizeof (Char));
    if (title != NULL) {
      tmp [0] = '\0';
      needsOpener = TRUE;

      GetTitle (fpfp->geneSymbol, str, sizeof (str));
      TrimSpacesAroundString (str);
      if (! StringHasNoText (str)) {
        if (needsOpener) {
          StringCat (title, "[gene=");
          needsOpener = FALSE;
        }
        StringCat (title, str);
      }
      GetTitle (fpfp->geneDesc, str, sizeof (str));
      TrimSpacesAroundString (str);
      if (! StringHasNoText (str)) {
        if (needsOpener) {
          StringCat (title, "[gene=");
          needsOpener = FALSE;
        } else {
          StringCat (title, ";");
        }
        StringCat (title, str);
      }
      if (! needsOpener) {
        StringCat (title, "] ");
      }

      needsOpener = TRUE;
      GetTitle (fpfp->protName, str, sizeof (str));
      TrimSpacesAroundString (str);
      if (! StringHasNoText (str)) {
        if (needsOpener) {
          StringCat (title, "[prot=");
          needsOpener = FALSE;
        }
        StringCat (title, str);
      }
      GetTitle (fpfp->protDesc, str, sizeof (str));
      TrimSpacesAroundString (str);
      if (! StringHasNoText (str)) {
        if (needsOpener) {
          StringCat (title, "[prot=");
          needsOpener = FALSE;
        } else {
          StringCat (title, ";");
        }
        StringCat (title, str);
      }
      if (! needsOpener) {
        StringCat (title, "] ");
      }

      needsOpener = TRUE;
      GetTitle (fpfp->comment, str, sizeof (str));
      TrimSpacesAroundString (str);
      if (! StringHasNoText (str)) {
        if (needsOpener) {
          StringCat (title, "[comment=");
          needsOpener = FALSE;
        }
        StringCat (title, str);
      }
      if (! needsOpener) {
        StringCat (title, "] ");
      }

      if (GetStatus (fpfp->orf)) {
        StringCat (title, "[orf] ");
      }

      ptr = SaveStringFromText (fpfp->title);
      StringCat (title, ptr);
      MemFree (ptr);

      if (pbsp->descr != NULL) {
        vnp = ValNodeFindNext (pbsp->descr, NULL, Seq_descr_title);
        if (vnp != NULL && vnp->data.ptrvalue != NULL) {
          vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
          ValNodeExtract (&(pbsp->descr), Seq_descr_title);
        }
      }
      vnp = CreateNewDescriptor (psep, Seq_descr_title);
      if (vnp != NULL) {
        vnp->data.ptrvalue = StringSave (title);
      }
    }
    MemFree (title);

    Remove (fpfp->form);
    Update ();
    if (sqfp->currConfirmSeq != NULL) {
      sqfp->currConfirmSeq = sqfp->currConfirmSeq->next;
    }
    LetUserFixProteinInfo (sqfp);
  }
}

static void FixProtMessage (ForM f, Int2 mssg)

{
  FixProtFormPtr  fpfp;

  fpfp = (FixProtFormPtr) GetObjectExtra (f);
  if (fpfp != NULL) {
    switch (mssg) {
      case VIB_MSG_QUIT :
        break;
      case VIB_MSG_CLOSE :
        break;
      case VIB_MSG_CUT :
        StdCutTextProc (NULL);
        break;
      case VIB_MSG_COPY :
        StdCopyTextProc (NULL);
        break;
      case VIB_MSG_PASTE :
        StdPasteTextProc (NULL);
        break;
      case VIB_MSG_DELETE :
        StdDeleteTextProc (NULL);
        break;
      default :
        if (fpfp->appmessage != NULL) {
          fpfp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static void FixUpProtActivate (WindoW w)

{
  IteM            closeItm;
  IteM            exportItm;
  FixProtFormPtr  fpfp;
  IteM            importItm;

  fpfp = (FixProtFormPtr) GetObjectExtra (w);
  if (fpfp != NULL) {
    if (fpfp->activate != NULL) {
      fpfp->activate (w);
    }
    importItm = FindFormMenuItem ((BaseFormPtr) fpfp, VIB_MSG_IMPORT);
    exportItm = FindFormMenuItem ((BaseFormPtr) fpfp, VIB_MSG_EXPORT);
    closeItm = FindFormMenuItem ((BaseFormPtr) fpfp, VIB_MSG_CLOSE);
    SafeDisable (importItm);
    SafeDisable (exportItm);
    SafeDisable (closeItm);
  }
}

static Boolean HasMinimalInformation (BioseqPtr pbsp)

{
  CharPtr     ptr;
  SeqIdPtr    sip;
  Char        str [128];
  CharPtr     title;
  ValNodePtr  vnp;

  if (pbsp == NULL) return FALSE;
  sip = SeqIdFindWorst (pbsp->id);
  SeqIdWrite (sip, str, PRINTID_REPORT, sizeof (str));
  if (StringHasNoText (str)) return FALSE;
  if (pbsp->descr == NULL) return FALSE;
  vnp = ValNodeFindNext (pbsp->descr, NULL, Seq_descr_title);
  if (vnp == NULL || vnp->data.ptrvalue == NULL) return FALSE;
  title = (CharPtr) vnp->data.ptrvalue;
  ptr = StringISearch (title, "[gene=");
  if (ptr == NULL) {
    ptr = StringISearch (title, "[gene_syn=");
  }
  if (ptr == NULL) return FALSE;
  StringNCpy_0 (str, ptr + 6, sizeof (str));
  ptr = StringChr (str, ']');
  if (ptr == NULL) return FALSE;
  *ptr = '\0';
  ptr = StringChr (str, ';');
  if (ptr != NULL) {
    *ptr = '\0';
    ptr++;
  }
  if (StringHasNoText (str)) return FALSE;
  ptr = StringISearch (title, "[prot=");
  if (ptr != NULL) {
    StringNCpy_0 (str, ptr + 6, sizeof (str));
    ptr = StringChr (str, ']');
  } else {
    ptr = StringISearch (title, "[protein=");
    if (ptr != NULL) {
      StringNCpy_0 (str, ptr + 9, sizeof (str));
      ptr = StringChr (str, ']');
    }
  }
  if (ptr == NULL) return FALSE;
  *ptr = '\0';
  ptr = StringChr (str, ';');
  if (ptr != NULL) {
    *ptr = '\0';
    ptr++;
  }
  if (StringHasNoText (str)) return FALSE;
  return TRUE;
}

static void LetUserFixProteinInfo (SequencesFormPtr sqfp)

{
  ButtoN             acceptBtn;
  GrouP              c;
  SeqEntryPtr        esep;
  FixProtFormPtr     fpfp;
  GrouP              g;
  GrouP              k;
  Int2               len;
  SeqEntryPtr        nsep;
  BioseqPtr          pbsp;
  SeqEntryPtr        psep;
  CharPtr            ptr;
  Boolean            quickmode;
  StdEditorProcsPtr  sepp;
  SeqIdPtr           sip;
  Char               str [256];
  CharPtr            title;
  ValNodePtr         vnp;
  WindoW             w;

  if (sqfp == NULL) return;

  quickmode = TRUE;
  if (GetAppParam ("SEQUIN", "PREFERENCES", "QUICKMODE", NULL, str, sizeof (str))) {
    if (StringICmp (str, "FALSE") == 0) {
      quickmode = FALSE;
    }
  }

  while (quickmode) {
    psep = sqfp->currConfirmSeq;
    if (psep == NULL) {
      if (sqfp->putItAllTogether != NULL) {
        sqfp->putItAllTogether (sqfp->form);
        return;
      }
    }
    pbsp = (BioseqPtr) psep->data.ptrvalue;
    if (pbsp != NULL) {
      if (HasMinimalInformation (pbsp)) {
        (sqfp->currConfirmCount)++;
        sqfp->currConfirmSeq = psep->next;
      } else {
        quickmode = FALSE;
      }
    } else {
      quickmode = FALSE;
    }
  }

  psep = sqfp->currConfirmSeq;
  if (psep == NULL) {
    if (sqfp->putItAllTogether != NULL) {
      sqfp->putItAllTogether (sqfp->form);
      return;
    }
  }

  esep = sqfp->topSeqForConfirm;
  if (esep == NULL || psep == NULL) return;
  nsep = FindNucSeqEntry (esep);
  if (nsep == NULL) return;
  pbsp = (BioseqPtr) psep->data.ptrvalue;
  if (pbsp == NULL) return;

  (sqfp->currConfirmCount)++;

  w = NULL;
  fpfp = (FixProtFormPtr) MemNew (sizeof (FixProtForm));
  if (fpfp != NULL) {

    fpfp->sqfp = sqfp;
    fpfp->esep = esep;
    fpfp->nsep = nsep;
    fpfp->psep = psep;
    fpfp->pbsp = pbsp;

    sprintf (str, "Product %d, length %ld amino acids",
             (int) sqfp->currConfirmCount, (long) pbsp->length);
    w = FixedWindow (-50, -33, -10, -10, str, NULL);
    SetObjectExtra (w, fpfp, StdCleanupFormProc);
    fpfp->form = (ForM) w;
    fpfp->formmessage = FixProtMessage;

#ifndef WIN_MAC
    CreateStdEditorFormMenus (w);
#endif

    fpfp->activate = NULL;
    sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
    if (sepp != NULL) {
      fpfp->activate = sepp->activateForm;
      fpfp->appmessage = sepp->handleMessages;
    }
    SetActivate (w, FixUpProtActivate);

    SelectFont (programFont);
    len = StringWidth ("Protein Name") + 2;
    SelectFont (systemFont);

    g = HiddenGroup (w, 1, 0, NULL);

    k = HiddenGroup (g, 2, 0, NULL);
    StaticPrompt (k, "Sequence ID", len, dialogTextHeight, programFont, 'l');
    fpfp->seqID = DialogText (k, "", 20, NULL);
    StaticPrompt (k, "Title", len, dialogTextHeight, programFont, 'l');
    fpfp->title = DialogText (k, "", 20, NULL);

    k = HiddenGroup (g, 2, 0, NULL);
    StaticPrompt (k, "Gene Symbol", len, dialogTextHeight, programFont, 'l');
    fpfp->geneSymbol = DialogText (k, "", 20, NULL);
    StaticPrompt (k, "Description", len, dialogTextHeight, programFont, 'l');
    fpfp->geneDesc = DialogText (k, "", 20, NULL);

    k = HiddenGroup (g, 2, 0, NULL);
    StaticPrompt (k, "Protein Name", len, dialogTextHeight, programFont, 'l');
    fpfp->protName = DialogText (k, "", 20, NULL);
    StaticPrompt (k, "Description", len, dialogTextHeight, programFont, 'l');
    fpfp->protDesc = DialogText (k, "", 20, NULL);

    k = HiddenGroup (g, 2, 0, NULL);
    StaticPrompt (k, "Comment", len, 3 * Nlm_stdLineHeight, programFont, 'l');
    fpfp->comment = ScrollText (k, 20, 3, programFont, TRUE, NULL);

    AlignObjects (ALIGN_RIGHT, (HANDLE) fpfp->seqID, (HANDLE) fpfp->title,
                  (HANDLE) fpfp->geneSymbol, (HANDLE) fpfp->geneDesc,
                  (HANDLE) fpfp->protName, (HANDLE) fpfp->protDesc,
                  (HANDLE) fpfp->comment, NULL);

    fpfp->orf = NULL;
    if (pbsp->descr != NULL) {
      vnp = ValNodeFindNext (pbsp->descr, NULL, Seq_descr_title);
      if (vnp != NULL && vnp->data.ptrvalue != NULL) {
        title = (CharPtr) vnp->data.ptrvalue;
        if (StringISearch (title, "[orf]") != NULL) {
          fpfp->orf = CheckBox (w, "Open Reading Frame", NULL);
        }
      }
    }

    c = HiddenGroup (w, 2, 0, NULL);
    acceptBtn = DefaultButton (c, "OK", AcceptThisFixup);
    SetObjectExtra (acceptBtn, fpfp, NULL);

    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, (HANDLE) fpfp->orf, NULL);

    RealizeWindow (w);

    str [0] = '\0';
    sip = SeqIdFindWorst (pbsp->id);
    if (sip != NULL) {
      if (sip->choice == SEQID_LOCAL) {
        SeqIdWrite (sip, str, PRINTID_REPORT, sizeof (str));
      } else {
        SeqIdWrite (sip, str, PRINTID_FASTA_SHORT, sizeof (str));
      }
    }
    SetTitle (fpfp->seqID, str);
    if (pbsp->descr != NULL) {
      vnp = ValNodeFindNext (pbsp->descr, NULL, Seq_descr_title);
      if (vnp != NULL && vnp->data.ptrvalue != NULL) {
        title = (CharPtr) vnp->data.ptrvalue;
        ptr = StringISearch (title, "[gene=");
        if (ptr != NULL) {
          StringNCpy_0 (str, ptr + 6, sizeof (str));
          ptr = StringChr (str, ']');
          if (ptr != NULL) {
            *ptr = '\0';
            ptr = StringChr (str, ';');
            if (ptr != NULL) {
              *ptr = '\0';
              ptr++;
            }
            SetTitle (fpfp->geneSymbol, str);
            SetTitle (fpfp->geneDesc, ptr);
          }
        }
        ptr = StringISearch (title, "[prot=");
        if (ptr != NULL) {
          StringNCpy_0 (str, ptr + 6, sizeof (str));
          ptr = StringChr (str, ']');
        } else {
          ptr = StringISearch (title, "[protein=");
          if (ptr != NULL) {
            StringNCpy_0 (str, ptr + 9, sizeof (str));
            ptr = StringChr (str, ']');
          }
        }
        if (ptr != NULL) {
          *ptr = '\0';
          ptr = StringChr (str, ';');
          if (ptr != NULL) {
            *ptr = '\0';
            ptr++;
          }
          SetTitle (fpfp->protName, str);
          SetTitle (fpfp->protDesc, ptr);
        }
        ptr = StringISearch (title, "[comment=");
        if (ptr != NULL) {
          StringNCpy_0 (str, ptr + 9, sizeof (str));
          ptr = StringChr (str, ']');
          if (ptr != NULL) {
            *ptr = '\0';
            SetTitle (fpfp->comment, str);
          }
        }
        ptr = StringISearch (title, "[orf]");
        if (ptr != NULL) {
          SetStatus (fpfp->orf, TRUE);
        } else {
          Hide (fpfp->orf);
        }
        ExciseString (title, "[gene=", "]");
        ExciseString (title, "[prot=", "]");
        ExciseString (title, "[protein=", "]");
        ExciseString (title, "[orf", "]");
        ExciseString (title, "[comment", "]");
        TrimSpacesAroundString (title);
        SetTitle (fpfp->title, title);
      }
    }

    Select (fpfp->seqID);
    Show (w);
    Select (w);
    Update ();
  }
}

extern void ConfirmSequencesFormParsing (ForM f, FormActnFunc putItAllTogether)

{
  BioseqSetPtr      bssp;
  FastaPagePtr      fpp;
  PhylipPagePtr     ppp;
  SeqEntryPtr       sep;
  SequencesFormPtr  sqfp;

  sqfp = (SequencesFormPtr) GetObjectExtra (f);
  if (sqfp != NULL) {
    sqfp->putItAllTogether = putItAllTogether;
    sqfp->topSeqForConfirm = NULL;
    sqfp->currConfirmSeq = NULL;
    if (sqfp->seqPackage < SEQ_PKG_POPULATION) {
      fpp = (FastaPagePtr) GetObjectExtra (sqfp->dnaseq);
      if (fpp != NULL) {
        sqfp->topSeqForConfirm = fpp->list;
      }
      sqfp->currConfirmSeq = NULL;
      fpp = (FastaPagePtr) GetObjectExtra (sqfp->protseq);
      if (fpp != NULL) {
        sqfp->currConfirmSeq = fpp->list;
      }
      sqfp->currConfirmCount = 0;
      LetUserFixProteinInfo (sqfp);
    } else {
      if (sqfp->seqFormat == SEQ_FMT_FASTA) {
        fpp = (FastaPagePtr) GetObjectExtra (sqfp->dnaseq);
        if (fpp != NULL) {
          sqfp->currConfirmSeq = fpp->list;
        }
      } else if (sqfp->seqFormat == SEQ_FMT_ALIGNMENT) {
        ppp = (PhylipPagePtr) GetObjectExtra (sqfp->dnaseq);
        if (ppp != NULL) {
          sep = ppp->sep;
          if (sep != NULL && IS_Bioseq_set (sep)) {
            bssp = (BioseqSetPtr) sep->data.ptrvalue;
            if (bssp != NULL) {
              sqfp->currConfirmSeq = bssp->seq_set;
            }
          }
        }
      }
      LetUserFixNucleotideInfo (sqfp);
    }
  }
}

extern void AddToSubSource (BioSourcePtr biop, CharPtr title, CharPtr label, Uint1 subtype)

{
  CharPtr       ptr;
  SubSourcePtr  ssp;
  CharPtr       str;
  SubSourcePtr  tmpssp;

  if (biop == NULL || title == NULL || label == NULL) return;
  str = MemNew (StringLen (title));
  if (str == NULL) return;
  ptr = StringISearch (title, label);
  if (ptr != NULL) {
    StringCpy (str, ptr + StringLen (label));
    ptr = StringChr (str, ']');
    if (ptr != NULL) {
      *ptr = '\0';
      TrimSpacesAroundString (str);
      ssp = SubSourceNew ();
      if (biop->subtype == NULL) {
        biop->subtype = ssp;
      } else {
        tmpssp = biop->subtype;
        while (tmpssp->next != NULL) {
          tmpssp = tmpssp->next;
        }
        tmpssp->next = ssp;
      }
      if (ssp != NULL) {
        ssp->subtype = subtype;
        ssp->name = StringSave (str);
      }
    }
  }
  MemFree (str);
}

extern void AddToOrgMod (BioSourcePtr biop, CharPtr title, CharPtr label, Uint1 subtype)

{
  OrgModPtr   mod;
  OrgNamePtr  onp;
  OrgRefPtr   orp;
  CharPtr     ptr;
  CharPtr     str;
  OrgModPtr   tmpmod;

  if (biop == NULL || title == NULL || label == NULL) return;
  str = MemNew (StringLen (title));
  if (str == NULL) return;
  ptr = StringISearch (title, label);
  if (ptr != NULL) {
    StringCpy (str, ptr + StringLen (label));
    ptr = StringChr (str, ']');
    if (ptr != NULL) {
      *ptr = '\0';
      TrimSpacesAroundString (str);
      orp = biop->org;
      if (orp == NULL) {
        orp = OrgRefNew ();
        biop->org = orp;
      }
      if (orp != NULL) {
        onp = orp->orgname;
        if (onp == NULL) {
          onp = OrgNameNew ();
          orp->orgname = onp;
        }
        if (onp != NULL) {
          mod = OrgModNew ();
          if (onp->mod == NULL) {
            onp->mod = mod;
          } else {
            tmpmod = onp->mod;
            while (tmpmod->next != NULL) {
              tmpmod = tmpmod->next;
            }
            tmpmod->next = mod;
          }
          if (mod != NULL) {
            mod->subtype = subtype;
            mod->subname = StringSave (str);
          }
        }
      }
    }
  }
  MemFree (str);
}

#define PROC_NUC_STR_SIZE 4096

extern Boolean ProcessOneNucleotideTitle (Int2 seqPackage, DialoG genbio, PopuP genome,
                                          PopuP gencode, SeqEntryPtr nsep, SeqEntryPtr top,
                                          BioSourcePtr masterbiop);
extern Boolean ProcessOneNucleotideTitle (Int2 seqPackage, DialoG genbio, PopuP genome,
                                          PopuP gencode, SeqEntryPtr nsep, SeqEntryPtr top,
                                          BioSourcePtr masterbiop)

{
  EnumFieldAssocPtr  ap;
  Uint1              biomol;
  BioSourcePtr       biop;
  BioseqSetPtr       bssp;
  Int2               code;
  CharPtr            lin;
  OrgNamePtr         masteronp = NULL;
  OrgRefPtr          masterorp = NULL;
  MolInfoPtr         mip;
  BioseqPtr          nbsp;
  Boolean            needbiop;
  OrgNamePtr         onp;
  OrgRefPtr          orp;
  CharPtr            ptr;
  SeqEntryPtr        sep;
  CharPtr            str;
  CharPtr            title;
  UIEnum             val;
  ValNodePtr         vnp;

  if (nsep == NULL || top == NULL) return FALSE;
  nbsp = (BioseqPtr) nsep->data.ptrvalue;
  if (nbsp == NULL) return FALSE;
  if (! ISA_na (nbsp->mol)) return FALSE;
  str = MemNew (PROC_NUC_STR_SIZE * sizeof (Char));
  if (str == NULL) return FALSE;
  sep = NULL;
  SeqEntryExplore (top, (Pointer) &sep, FindFirstSeqEntryTitle);
  if (sep != NULL) {
    vnp = SeqEntryGetSeqDescr (sep, Seq_descr_title, NULL);
    if (vnp != NULL && vnp->data.ptrvalue != NULL) {
      title = (CharPtr) vnp->data.ptrvalue;
      needbiop = FALSE;
      if (seqPackage >= SEQ_PKG_POPULATION && seqPackage <= SEQ_PKG_GENBANK) {
        needbiop = TRUE;
        if (GetAppParam ("SEQUIN", "PREFERENCES", "BIOSRCONALL", NULL, str, PROC_NUC_STR_SIZE)) {
          if (StringICmp (str, "FALSE") == 0) {
            needbiop = FALSE;
          }
        }
      }
      if ((! needbiop) && StringISearch (title, "[") != NULL) {
        for (ap = orgmod_subtype_alist; ap->name != NULL; ap++) {
          MakeSearchStringFromAlist (str, ap->name);
          if (StringISearch (title, str) != NULL) {
            needbiop = TRUE;
          }
        }
        for (ap = subsource_subtype_alist; ap->name != NULL; ap++) {
          MakeSearchStringFromAlist (str, ap->name);
          if (StringISearch (title, str) != NULL) {
            needbiop = TRUE;
          }
        }
        if (StringISearch (title, "[note=") != NULL) {
          needbiop = TRUE;
        }
        if (StringISearch (title, "[comment=") != NULL) {
          needbiop = TRUE;
        }
        if (StringISearch (title, "[subsource=") != NULL) {
          needbiop = TRUE;
        }
        /*
        if (StringISearch (title, "[strain=") != NULL ||
            StringISearch (title, "[isolate=") != NULL ||
            StringISearch (title, "[clone=") != NULL) {
          needbiop = TRUE;
        }
        */
      }
      biop = NULL;
      ptr = StringISearch (title, "[org=");
      if (ptr != NULL) {
        StringNCpy_0 (str, ptr + 5, PROC_NUC_STR_SIZE);
      } else {
        ptr = StringISearch (title, "[organism=");
        if (ptr != NULL) {
          StringNCpy_0 (str, ptr + 10, PROC_NUC_STR_SIZE);
        }
      }
      if (ptr != NULL) {
        ptr = StringChr (str, ']');
        if (ptr != NULL) {
          *ptr = '\0';
          if (SetBioSourceDialogTaxName (genbio, str)) {
            biop = (BioSourcePtr) DialogToPointer (genbio);
          } else {
            biop = BioSourceNew ();
            if (biop != NULL) {
              orp = OrgRefNew ();
              biop->org = orp;
              if (orp != NULL) {
                TrimSpacesAroundString (str);
                SetTaxNameAndRemoveTaxRef (orp, StringSave (str));
                if (masterbiop != NULL) {
                  masterorp = masterbiop->org;
                  if (masterorp != NULL) {
                    masteronp = masterorp->orgname;
                    if (masteronp != NULL) {
                      onp = OrgNameNew ();
                      orp->orgname = onp;
                      if (onp != NULL) {
                        onp->gcode = masteronp->gcode;
                        onp->mgcode = masteronp->mgcode;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      } else if (needbiop && masterbiop != NULL) {
        masterorp = masterbiop->org;
        if (masterorp != NULL) {
          if (SetBioSourceDialogTaxName (genbio, masterorp->taxname)) {
            biop = (BioSourcePtr) DialogToPointer (genbio);
          } else {
            biop = BioSourceNew ();
            if (biop != NULL) {
              orp = OrgRefNew ();
              biop->org = orp;
              if (orp != NULL) {
                SetTaxNameAndRemoveTaxRef (orp, StringSave (masterorp->taxname));
                orp->common = StringSave (masterorp->common);
                masteronp = masterorp->orgname;
                if (masteronp != NULL) {
                  onp = OrgNameNew ();
                  orp->orgname = onp;
                  if (onp != NULL) {
                    onp->gcode = masteronp->gcode;
                    onp->mgcode = masteronp->mgcode;
                  }
                }
              }
            }
          }
        }
      }
      if (biop != NULL) {
        ptr = StringISearch (title, "[lineage=");
        if (ptr != NULL) {
          lin = StringSave (ptr + 9);
          if (lin != NULL) {
            ptr = StringChr (lin, ']');
            if (ptr != NULL) {
              *ptr = '\0';
              orp = biop->org;
              if (orp != NULL) {
                onp = orp->orgname;
                if (onp == NULL) {
                  onp = OrgNameNew ();
                  orp->orgname = onp;
                }
                if (onp != NULL) {
                  onp->lineage = MemFree (onp->lineage);
                  onp->lineage = StringSave (lin);
                }
              }
            }
            MemFree (lin);
          }
        }
        if (seqPackage == SEQ_PKG_PHYLOGENETIC) {
          orp = biop->org;
          if (orp != NULL) {
            onp = orp->orgname;
            if (onp == NULL) {
              onp = OrgNameNew ();
              orp->orgname = onp;
            }
            if (onp != NULL) {
              if (onp->gcode == 0 && onp->mgcode == 0) {
                code = gcIndexToId [GetValue (gencode)];
                if (GetEnumPopup (genome, biosource_genome_simple_alist, &val)) {
                  if (val == 4 || val == 5) {
                    onp->mgcode = gcIndexToId [GetValue (gencode)];
                  } else {
                    onp->gcode = gcIndexToId [GetValue (gencode)];
                  }
                  biop->genome = (Uint1) val;
                }
              } else if (GetEnumPopup (genome, biosource_genome_simple_alist, &val)) {
                biop->genome = (Uint1) val;
              }
            }
          }
        }
        vnp = CreateNewDescriptor (top, Seq_descr_source);
        if (vnp != NULL) {
          vnp->data.ptrvalue = (Pointer) biop;
        }
        for (ap = orgmod_subtype_alist; ap->name != NULL; ap++) {
          MakeSearchStringFromAlist (str, ap->name);
          AddToOrgMod (biop, title, str, (Int2) ap->value);
          ExciseString (title, str, "]");
        }
        for (ap = subsource_subtype_alist; ap->name != NULL; ap++) {
          MakeSearchStringFromAlist (str, ap->name);
          AddToSubSource (biop, title, str, (Int2) ap->value);
          ExciseString (title, str, "]");
        }
        AddToOrgMod (biop, title, "[note=", 255);
        ExciseString (title, "[note=", "]");
        AddToOrgMod (biop, title, "[comment=", 255);
        ExciseString (title, "[comment=", "]");
        AddToSubSource (biop, title, "[subsource=", 255);
        ExciseString (title, "[subsource=", "]");
        /*
        AddToOrgMod (biop, title, "[strain=", 2);
        AddToOrgMod (biop, title, "[isolate=", 17);
        AddToSubSource (biop, title, "[clone=", 3);
        ExciseString (title, "[strain=", "]");
        ExciseString (title, "[isolate=", "]");
        ExciseString (title, "[clone=", "]");
        */
        ptr = StringISearch (title, "[molecule=");
        if (ptr != NULL) {
          StringNCpy_0 (str, ptr + 10, PROC_NUC_STR_SIZE);
          ptr = StringChr (str, ']');
          if (ptr != NULL) {
            *ptr = '\0';
            if (StringICmp (str, "dna") == 0) {
              nbsp->mol = Seq_mol_dna;
            } else if (StringICmp (str, "rna") == 0) {
              nbsp->mol = Seq_mol_rna;
            }
          }
        }
        ptr = StringISearch (title, "[location=");
        if (ptr != NULL) {
          StringNCpy_0 (str, ptr + 10, PROC_NUC_STR_SIZE);
          ptr = StringChr (str, ']');
          if (ptr != NULL) {
            *ptr = '\0';
            if (StringICmp (str, "Mitochondrial") == 0) { /* alternative spelling */
              biop->genome = 5;
            }
            for (ap = biosource_genome_simple_alist; ap->name != NULL; ap++) {
              if (StringICmp (str, ap->name) == 0) {
                biop->genome = (Uint1) ap->value;
              }
            }
          }
        }
        ptr = StringISearch (title, "[moltype=");
        if (ptr != NULL) {
          biomol = 0;
          StringNCpy_0 (str, ptr + 9, PROC_NUC_STR_SIZE);
          ptr = StringChr (str, ']');
          if (ptr != NULL) {
            *ptr = '\0';
            if (StringICmp (str, "genomic") == 0) {
              biomol = MOLECULE_TYPE_GENOMIC;
            } else if (StringICmp (str, "mRNA") == 0) {
              biomol = MOLECULE_TYPE_MRNA;
            }
            if (biomol != 0) {
              mip = MolInfoNew ();
              if (mip != NULL) {
                mip->biomol = biomol;
                vnp = CreateNewDescriptor (nsep, Seq_descr_molinfo);
                if (vnp != NULL) {
                  vnp->data.ptrvalue = (Pointer) mip;
                }
              }
            }
          }
        }
        /*
        for (ap = biosource_genome_simple_alist; ap->name != NULL; ap++) {
          MakeSearchStringFromAlist (str, ap->name);
          ptr = StringISearch (str, "=");
          if (ptr != NULL) {
            *ptr = '\0';
          }
          ptr = StringISearch (title, str);
          if (ptr != NULL) {
            biop->genome = (Uint1) ap->value;
          }
          ExciseString (title, str, "]");
        }
        ptr = StringISearch (title, "[dna]");
        if (ptr != NULL) {
          nbsp->mol = Seq_mol_dna;
          ExciseString (title, "[dna", "]");
        }
        ptr = StringISearch (title, "[rna]");
        if (ptr != NULL) {
          nbsp->mol = Seq_mol_rna;
          ExciseString (title, "[rna", "]");
        }
        */
      }
      ExciseString (title, "[org=", "]");
      ExciseString (title, "[organism=", "]");
      ExciseString (title, "[lineage=", "]");
      ExciseString (title, "[molecule=", "]");
      ExciseString (title, "[moltype=", "]");
      ExciseString (title, "[location=", "]");
      TrimSpacesAroundString (title);
      if (StringHasNoText (title) || sep != top) {
        vnp = NULL;
        if (IS_Bioseq (sep)) {
          nbsp = (BioseqPtr) sep->data.ptrvalue;
          vnp = ValNodeExtract (&(nbsp->descr), Seq_descr_title);
        } else if (IS_Bioseq_set (sep)) {
          bssp = (BioseqSetPtr) sep->data.ptrvalue;
          vnp = ValNodeExtract (&(bssp->descr), Seq_descr_title);
        }
        if (vnp != NULL && StringHasNoText ((CharPtr) vnp->data.ptrvalue)) {
          vnp = ValNodeFreeData (vnp);
        }
        if (sep != top && vnp != NULL) {
          if (IS_Bioseq (top)) {
            nbsp = (BioseqPtr) top->data.ptrvalue;
            ValNodeLink (&(nbsp->descr), vnp);
          } else if (IS_Bioseq_set (top)) {
            bssp = (BioseqSetPtr) top->data.ptrvalue;
            ValNodeLink (&(bssp->descr), vnp);
          }
        }
      }
    }
  } else {
    needbiop = FALSE;
    if (seqPackage >= SEQ_PKG_POPULATION && seqPackage <= SEQ_PKG_GENBANK) {
      needbiop = TRUE;
      if (GetAppParam ("SEQUIN", "PREFERENCES", "BIOSRCONALL", NULL, str, PROC_NUC_STR_SIZE)) {
        if (StringICmp (str, "FALSE") == 0) {
          needbiop = FALSE;
        }
      }
    }
    if (needbiop && masterbiop != NULL) {
      masterorp = masterbiop->org;
      if (masterorp != NULL) {
        biop = NULL;
        if (SetBioSourceDialogTaxName (genbio, masterorp->taxname)) {
          biop = (BioSourcePtr) DialogToPointer (genbio);
        } else {
          biop = BioSourceNew ();
          if (biop != NULL) {
            orp = OrgRefNew ();
            biop->org = orp;
            if (orp != NULL) {
              SetTaxNameAndRemoveTaxRef (orp, StringSave (masterorp->taxname));
              orp->common = StringSave (masterorp->common);
              onp = orp->orgname;
              if (onp == NULL) {
                orp->orgname = OrgNameNew ();
                onp = orp->orgname;
              }
              if (onp != NULL) {
                masteronp = masterorp->orgname;
                if (masteronp != NULL) {
                  if (onp->gcode == 0 && onp->mgcode == 0) {
                    onp->gcode = masteronp->gcode;
                    onp->mgcode = masteronp->mgcode;
                  }
                }
                biop->genome = masterbiop->genome;
                if (onp->gcode == 0 && onp->mgcode == 0) {
                  code = gcIndexToId [GetValue (gencode)];
                  if (GetEnumPopup (genome, biosource_genome_simple_alist, &val)) {
                    if (val == 4 || val == 5) {
                      onp->mgcode = gcIndexToId [GetValue (gencode)];
                    } else {
                      onp->gcode = gcIndexToId [GetValue (gencode)];
                    }
                    biop->genome = (Uint1) val;
                  }
                }
              }
            }
          }
        }
        if (biop != NULL) {
          vnp = CreateNewDescriptor (nsep, Seq_descr_source);
          if (vnp != NULL) {
            vnp->data.ptrvalue = (Pointer) biop;
          }
        }
      }
    }
  }
  MemFree (str);
  return TRUE;
}

static Boolean AutomaticNucleotideProcess (SequencesFormPtr sqfp, SeqEntryPtr nsep,
                                           SeqEntryPtr top, BioSourcePtr masterbiop)

{
  BioseqSetPtr  bssp;
  Boolean       rsult;
  SeqEntryPtr   tmp;

  if (sqfp == NULL || nsep == NULL || top == NULL) return FALSE;
  if (IS_Bioseq_set (nsep)) {
    bssp = (BioseqSetPtr) nsep->data.ptrvalue;
    rsult = FALSE;
    if (bssp != NULL) {
      for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
        if (AutomaticNucleotideProcess (sqfp, tmp, top, masterbiop)) {
          rsult = TRUE;
        }
      }
    }
    return rsult;
  }
  return ProcessOneNucleotideTitle (sqfp->seqPackage, sqfp->genbio, sqfp->genome,
                                    sqfp->gencode, nsep, top, masterbiop);
}

typedef struct idlist {
  BioseqPtr  bsp;
  CharPtr    key;
  struct idlist PNTR left;
  struct idlist PNTR right;
} IdList, PNTR IdListPtr;

static void BuildTree (IdListPtr PNTR head, BioseqPtr bsp, CharPtr x)

{
  Int2       comp;
  IdListPtr  idlist;
  SeqIdPtr   sip;
  Char       str [64];

  if (*head != NULL) {
    idlist = *head;
    comp = StringICmp (idlist->key, x);
    if (comp < 0) {
      BuildTree (&(idlist->right), bsp, x);
    } else if (comp > 0) {
      BuildTree (&(idlist->left), bsp, x);
    } else {
      sip = MakeNewProteinSeqId (NULL, NULL);
      if (sip != NULL) {
        bsp->id = SeqIdFree (bsp->id);
        bsp->id = sip;
        SeqMgrReplaceInBioseqIndex (bsp);
        sip = SeqIdFindWorst (bsp->id);
        SeqIdWrite (sip, str, PRINTID_REPORT, sizeof (str));
        BuildTree (head, bsp, str);
      }
    }
  } else {
    idlist = MemNew (sizeof (IdList));
    if (idlist != NULL) {
      *head = idlist;
      idlist->bsp = bsp;
      sip = SeqIdFindWorst (bsp->id);
      SeqIdWrite (sip, str, PRINTID_REPORT, sizeof (str));
      idlist->key = StringSave (str);
      idlist->left = NULL;
      idlist->right = NULL;
    }
  }
}

static void FreeTree (IdListPtr PNTR head)

{
  IdListPtr  idlist;

  if (head != NULL && *head != NULL) {
    idlist = *head;
    FreeTree (&(idlist->left));
    FreeTree (&(idlist->right));
    MemFree (idlist->key);
    MemFree (idlist);
  }
}

static void ResolveCollidingIDs (IdListPtr PNTR head, SeqEntryPtr list)

{
  BioseqPtr  bsp;
  SeqIdPtr   sip;
  Char       str [64];

  if (head == NULL) return;
  while (list != NULL) {
    if (IS_Bioseq (list)) {
      bsp = (BioseqPtr) list->data.ptrvalue;
      if (bsp != NULL) {
        sip = SeqIdFindWorst (bsp->id);
        SeqIdWrite (sip, str, PRINTID_REPORT, sizeof (str));
        BuildTree (head, bsp, str);
      }
    }
    list = list->next;
  }
}


static void PutMolInfoOnSeqEntry (SequencesFormPtr sqfp, SeqEntryPtr sep, Uint1 biomol)

{
  BioseqPtr    bsp;
  BioseqSetPtr bssp;
  MolInfoPtr   mip;
  Boolean      partial5;
  Boolean      partial3;
  UIEnum       val;
  ValNodePtr   vnp;


  if (sqfp != NULL && sep != NULL) {
    if (IS_Bioseq_set (sep))
    {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) 
      {
      	PutMolInfoOnSeqEntry (sqfp, sep, biomol);
      }
      return;
    }

    mip = MolInfoNew ();
    if (mip != NULL) {
      mip->biomol = biomol;
      partial5 = GetStatus (sqfp->partial5);
      partial3 = GetStatus (sqfp->partial3);
      if (partial5 && partial3) {
        mip->completeness = 5;
      } else if (partial5) {
        mip->completeness = 3;
      } else if (partial3) {
        mip->completeness = 4;
      }
      vnp = CreateNewDescriptor (sep, Seq_descr_molinfo);
      if (vnp != NULL) {
        vnp->data.ptrvalue = (Pointer) mip;
        if (mip->biomol == 1) {
          sqfp->dnamolfrommolinfo = Seq_mol_na;
        } else if (mip->biomol >= 2  && mip->biomol <= 7) {
          sqfp->dnamolfrommolinfo = Seq_mol_rna;
        } else if (mip->biomol == 9) {
          sqfp->dnamolfrommolinfo = Seq_mol_dna;
        } else if (mip->biomol == 253) {
          sqfp->dnamolfrommolinfo = Seq_mol_dna;
          mip->biomol = 1;
        } else if (mip->biomol == 254) {
          sqfp->dnamolfrommolinfo = Seq_mol_rna;
          mip->biomol = 1;
        } else if (mip->biomol == 255) {
          sqfp->dnamolfrommolinfo = Seq_mol_other;
        }
      }
      if (sqfp->seqPackage > SEQ_PKG_GENOMICCDNA && IS_Bioseq (sep))
      {
        bsp = (BioseqPtr) sep->data.ptrvalue;
        if (bsp != NULL) {
          if (bsp->mol != Seq_mol_dna && bsp->mol != Seq_mol_rna) {
            bsp->mol = (Uint1) sqfp->dnamolfrommolinfo;
          }
          if (GetEnumPopup (sqfp->topologyPopup, topology_nuc_alist, &val)) {
            bsp->topology = (Uint1) val;
          }
        }		
      }
    }
  }
}

static void PrefixOrgToDefline (SeqEntryPtr sep)

{
  BioSourcePtr  biop;
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  CharPtr       def;
  OrgRefPtr     orp;
  CharPtr       ptr;
  CharPtr       str;
  Char          taxname [64];
  ValNodePtr    ttl;
  ValNodePtr    vnp;

  if (sep == NULL) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && (bssp->_class == 7 ||
                         (IsPopPhyEtcSet (bssp->_class)))) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        PrefixOrgToDefline (sep);
      }
      return;
    }
  }

  if (! IS_Bioseq (sep)) return;
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL) return;

  taxname [0] = '\0';
  orp = NULL;
  biop = NULL;
  ttl = NULL;
  vnp = bsp->descr;
  for (vnp = bsp->descr; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == Seq_descr_source) {
      biop = (BioSourcePtr) vnp->data.ptrvalue;
    } else if (vnp->choice == Seq_descr_org) {
      orp = (OrgRefPtr) vnp->data.ptrvalue;
    } else if (vnp->choice == Seq_descr_title) {
      ttl = vnp;
    }
  }
  if (orp == NULL && biop != NULL) {
    orp = biop->org;
  }
  if (orp == NULL) return;
  if (ttl == NULL) return;
  StringNCpy_0 (taxname, orp->taxname, sizeof (taxname));
  ptr = StringSearch (taxname, "(");
  if (ptr != NULL) {
    *ptr = '\0';
  }
  TrimSpacesAroundString (taxname);
  if ((StringICmp (taxname, "Human immunodeficiency virus type 1") == 0) ||
      (StringICmp (taxname, "Human immunodeficiency virus 1") == 0)) {
    StringCpy (taxname, "HIV-1");
  } else if ((StringICmp (taxname,"Human immunodeficiency virus type 2")==0) ||
	     (StringICmp (taxname,"Human immunodeficiency virus 2")==0)) {
    StringCpy (taxname, "HIV-2");
  }

  def = (CharPtr) ttl->data.ptrvalue;
  if (StringHasNoText (def)) return;

  ptr = StringISearch (def, taxname);
  if (ptr != NULL && ptr == def) return;
  str = MemNew ((StringLen (taxname) + StringLen (def) + 4) * sizeof (Char));
  if (str == NULL) return;
  StringCpy (str, taxname);
  StringCat (str, " ");
  StringCat (str, def);
  ttl->data.ptrvalue = MemFree (ttl->data.ptrvalue);
  ttl->data.ptrvalue = str;
}

static CharPtr onecomponent = "\
Multiple sequence components are expected in this submission.\n\
They should all be read in at the same time from the same file.";

static void OnlyOneComponentWarning (SequencesFormPtr sqfp)

{
  CharPtr  type;

  if (sqfp != NULL) {
    if (sqfp->seqPackage == SEQ_PKG_GENOMICCDNA) return;
    switch (sqfp->seqPackage) {
      case SEQ_PKG_SEGMENTED :
        type = "segmented sequence";
        break;
      case SEQ_PKG_POPULATION :
        type = "population set";
        break;
      case SEQ_PKG_PHYLOGENETIC :
        type = "phylogenetic set";
        break;
      case SEQ_PKG_MUTATION :
        type = "mutation set";
        break;
      case SEQ_PKG_ENVIRONMENT :
        type = "environmental samples";
        break;
      case SEQ_PKG_GENBANK :
        type = "batch submission";
        break;
      default :
        type = "unknown set";
        break;
    }
    Message (MSG_OK, "WARNING - There is only one component in this %s.\n%s",
             type, onecomponent);
  }
}

static void AutomaticMrnaProcess (SeqEntryPtr nucsep, SeqEntryPtr mrnasep, Boolean partial5, Boolean partial3)

{
  CharPtr     allele = NULL;
  BioseqPtr   bsp;
  GeneRefPtr  grp;
  SeqLocPtr   gslp;
  Boolean     hasNulls;
  MolInfoPtr  mip;
  BioseqPtr   mrnabsp;
  BioseqPtr   nucbsp;
  SeqLocPtr   oldslp;
  CharPtr     ptr;
  RnaRefPtr   rrp;
  SeqFeatPtr  sfp;
  SeqIdPtr    sip;
  SeqLocPtr   slp;
  Char        str [128];
  CharPtr     ttl;
  ValNodePtr  vnp;

  if (nucsep == NULL || mrnasep == NULL) return;
  if (IS_Bioseq (nucsep) && IS_Bioseq (mrnasep)) {
    nucbsp = (BioseqPtr) nucsep->data.ptrvalue;
    mrnabsp = (BioseqPtr) mrnasep->data.ptrvalue;
    if (nucbsp == NULL || mrnabsp == NULL) return;
    slp = AlignmRNA2genomic (nucbsp, mrnabsp);
    if (slp == NULL) return;
    sip = SeqLocId (slp);
    if (sip != NULL) {
      bsp = BioseqFind (sip);
      if (bsp != NULL) {
        if (bsp->repr == Seq_repr_seg) {
          oldslp = slp;
          slp = SegLocToParts (bsp, oldslp);
          FreeAllFuzz (slp);
          SeqLocFree (oldslp);
        }
      }
    }
    StripLocusFromSeqLoc (slp);
    str [0] = '\0';
    ttl = NULL;
    vnp = ValNodeFindNext (mrnabsp->descr, NULL, Seq_descr_title);
    if (vnp != NULL) {
      ttl = (CharPtr) vnp->data.ptrvalue;
    }
    if (ttl != NULL) {
      ptr = StringISearch (ttl, "[gene=");
      if (ptr != NULL) {
        StringNCpy_0 (str, ptr + 6, sizeof (str));
        ptr = StringChr (str, ']');
        if (ptr != NULL) {
          *ptr = '\0';
          ptr = StringChr (str, ';');
          if (ptr != NULL) {
            *ptr = '\0';
            ptr++;
            allele = StringChr (ptr, ';');
            if (allele != NULL) {
              *allele = '\0';
              allele++;
            }
          }
          grp = CreateNewGeneRef (str, allele, ptr, FALSE);
          if (grp != NULL) {
            if (ExtendGene (grp, nucsep, slp)) {
              grp = GeneRefFree (grp);
            } else {
              sfp = CreateNewFeature (nucsep, NULL, SEQFEAT_GENE, NULL);
              if (sfp != NULL) {
                sfp->data.value.ptrvalue = (Pointer) grp;
                sfp->location = SeqLocFree (sfp->location);
                sfp->location = AsnIoMemCopy ((Pointer) slp,
                                              (AsnReadFunc) SeqLocAsnRead,
                                              (AsnWriteFunc) SeqLocAsnWrite);
                sip = SeqLocId (sfp->location);
                if (sip != NULL) {
                  bsp = BioseqFind (sip);
                } else {
                  bsp = nucbsp;
                }
                if (bsp != NULL) {
                  gslp = SeqLocMerge (bsp, sfp->location, NULL, TRUE, FALSE, FALSE);
                  if (gslp != NULL) {
                    sfp->location = SeqLocFree (sfp->location);
                    sfp->location = gslp;
                    if (bsp->repr == Seq_repr_seg) {
                      gslp = SegLocToPartsEx (bsp, sfp->location, TRUE);
                      sfp->location = SeqLocFree (sfp->location);
                      sfp->location = gslp;
                      hasNulls = LocationHasNullsBetween (sfp->location);
                      sfp->partial = (sfp->partial || hasNulls);
                    }
                    FreeAllFuzz (gslp);
                  }
                }
              }
            }
          }
        }
      }
      str [0] = '\0';
      ptr = StringISearch (ttl, "[mrna=");
      if (ptr != NULL) {
        StringNCpy_0 (str, ptr + 6, sizeof (str));
        ptr = StringChr (str, ']');
        if (ptr != NULL) {
          *ptr = '\0';
        }
      } else {
        ptr = StringISearch (ttl, "[cdna=");
        if (ptr != NULL) {
          StringNCpy_0 (str, ptr + 6, sizeof (str));
          ptr = StringChr (str, ']');
          if (ptr != NULL) {
            *ptr = '\0';
          }
        }
      }
    }
    rrp = RnaRefNew ();
    if (rrp != NULL) {
      rrp->type = 2;
      if (! StringHasNoText (str)) {
        rrp->ext.choice = 1;
        rrp->ext.value.ptrvalue = StringSave (str);
      }
      sfp = CreateNewFeature (nucsep, NULL, SEQFEAT_RNA, NULL);
      if (sfp != NULL) {
        sfp->data.value.ptrvalue = (Pointer) rrp;
        sfp->location = SeqLocFree (sfp->location);
        sfp->location = AsnIoMemCopy ((Pointer) slp,
                                      (AsnReadFunc) SeqLocAsnRead,
                                      (AsnWriteFunc) SeqLocAsnWrite);
        SetSeqFeatProduct (sfp, mrnabsp);
        SetSeqLocPartial (sfp->location, partial5, partial3);
        sfp->partial = (sfp->partial || partial5 || partial3);
        if (ttl != NULL) {
          ptr = StringISearch (ttl, "[comment=");
          if (ptr != NULL) {
            StringNCpy_0 (str, ptr + 9, sizeof (str));
            ptr = StringChr (str, ']');
            if (ptr != NULL) {
              *ptr = '\0';
            }
            if (! StringHasNoText (str)) {
              sfp->comment = StringSave (str);
            }
          }
        }
      }
    }
    SeqLocFree (slp);
    ExciseString (ttl, "[mrna", "]");
    ExciseString (ttl, "[cdna", "]");
    ExciseString (ttl, "[gene", "]");
    ExciseString (ttl, "[comment", "]");
    TrimSpacesAroundString (ttl);
    if (StringHasNoText (ttl)) {
      ValNodeExtract (&(mrnabsp->descr), Seq_descr_title);
    }
    mip = MolInfoNew ();
    if (mip != NULL) {
      mip->biomol = 3;
      if (partial5 && partial3) {
        mip->completeness = 5;
      } else if (partial5) {
        mip->completeness = 3;
      } else if (partial3) {
        mip->completeness = 4;
      }
      vnp = CreateNewDescriptor (mrnasep, Seq_descr_molinfo);
      if (vnp != NULL) {
        vnp->data.ptrvalue = (Pointer) mip;
      }
    }
    mrnabsp->mol = Seq_mol_rna;
  }
}

static void ExciseStringFromBioseq (SeqEntryPtr sep, CharPtr from, CharPtr to)

{
  BioseqPtr   bsp;
  CharPtr     title;
  ValNodePtr  vnp;

  if (sep == NULL) return;
  if (! IS_Bioseq (sep)) return;
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL || bsp->descr == NULL) return;
  vnp = SeqEntryGetSeqDescr (sep, Seq_descr_title, NULL);
  if (vnp == NULL) return;
  title = (CharPtr) vnp->data.ptrvalue;
  if (title == NULL) return;
  ExciseString (title, from, to);
  TrimSpacesAroundString (title);
  if (StringHasNoText (title)) {
    ValNodeExtract (&(bsp->descr), Seq_descr_title);
  }
}

static Boolean LookForStringInBioseq (SeqEntryPtr sep, Uint1 mol, CharPtr str, CharPtr tmp, size_t maxsize)

{
  BioseqPtr   bsp;
  CharPtr     title;
  ValNodePtr  vnp;

  if (sep == NULL || tmp == NULL) return FALSE;
  if (! IS_Bioseq (sep)) return FALSE;
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL || bsp->mol != mol || bsp->descr == NULL) return FALSE;
  vnp = ValNodeFindNext (bsp->descr, NULL, Seq_descr_title);
  if (vnp == NULL || vnp->data.ptrvalue == NULL) return FALSE;
  title = (CharPtr) vnp->data.ptrvalue;
  return LookForSearchString (title, str, tmp, maxsize);
}

static void FindBioseqWithString (SeqEntryPtr sep, Uint1 mol, CharPtr tag, CharPtr str, SeqEntryPtr PNTR rsult)

{
  BioseqPtr     bsp = NULL;
  BioseqSetPtr  bssp = NULL;
  Char          tmp [256];

  if (sep == NULL || sep->data.ptrvalue == NULL || rsult == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (LookForStringInBioseq (sep, mol, tag, tmp, sizeof (tmp) - 1)) {
      if (StringICmp (str, tmp) == 0) {
        *rsult = sep;
        return;
      }
    }
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
      FindBioseqWithString (sep, mol, tag, str, rsult);
    }
  }
}

static SeqEntryPtr FindRnaByRefOnRna (SeqEntryPtr sep, SeqEntryPtr psep)

{
  SeqEntryPtr  msep;
  Char         tmp [256];

  msep = NULL;
  if (sep == NULL || psep == NULL) return NULL;
  if (LookForStringInBioseq (psep, Seq_mol_aa, "[prot=", tmp, sizeof (tmp) - 1)) {
    FindBioseqWithString (sep, Seq_mol_rna, "[prot=", tmp, &msep);
    if (msep != NULL) {
      ExciseStringFromBioseq (msep, "[prot=", "]");
      return msep;
    }
  } else if (LookForStringInBioseq (psep, Seq_mol_aa, "[protein=", tmp, sizeof (tmp) - 1)) {
    FindBioseqWithString (sep, Seq_mol_rna, "[protein=", tmp, &msep);
    if (msep != NULL) {
      ExciseStringFromBioseq (msep, "[protein=", "]");
      return msep;
    }
  }
  return msep;
}

static void FindRnaByName (SeqEntryPtr sep, CharPtr str, SeqEntryPtr PNTR msep)

{
  BioseqPtr     bsp = NULL;
  BioseqSetPtr  bssp = NULL;
  RnaRefPtr     rrp;
  SeqAnnotPtr   sap;
  SeqFeatPtr    sfp;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  if (str == NULL || msep == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
  } else return;
  while (sap != NULL) {
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      while (sfp != NULL) {
        if (sfp->data.choice == SEQFEAT_RNA) {
          rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
          if (rrp != NULL && rrp->type == 2 && rrp->ext.choice == 1 && sfp->product != NULL) {
            if (StringICmp (rrp->ext.value.ptrvalue, str) == 0) {
              bsp = BioseqFind (SeqLocId (sfp->product));
              if (bsp != NULL) {
                *msep = SeqMgrGetSeqEntryForData (bsp);
              }
            }
          }
        }
        sfp = sfp->next;
      }
    }
    sap = sap->next;
  }
  if (bssp != NULL) {
    for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
      FindRnaByName (sep, str, msep);
    }
  }
}

static SeqEntryPtr FindRnaByRefOnProtein (SeqEntryPtr sep, SeqEntryPtr psep)

{
  SeqEntryPtr  msep;
  Char         tmp [256];

  msep = NULL;
  if (sep == NULL || psep == NULL) return NULL;
  if (LookForStringInBioseq (psep, Seq_mol_aa, "[mrna=", tmp, sizeof (tmp) - 1)) {
    FindRnaByName (sep, tmp, &msep);
    if (msep != NULL) {
      ExciseStringFromBioseq (psep, "[mrna=", "]");
      return msep;
    }
  }
  return msep;
}

static void FindRnaByLocationOverlap (SeqEntryPtr sep, SeqLocPtr slp,
                                      Int4Ptr mindiff, SeqEntryPtr PNTR msep)

{
  BioseqPtr     bsp = NULL;
  BioseqSetPtr  bssp = NULL;
  Int4          diff;
  RnaRefPtr     rrp;
  SeqAnnotPtr   sap;
  SeqFeatPtr    sfp;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  if (slp == NULL || mindiff == NULL || msep == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
  } else return;
  while (sap != NULL) {
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      while (sfp != NULL) {
        if (sfp->data.choice == SEQFEAT_RNA) {
          rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
          if (rrp != NULL && rrp->type == 2 && sfp->product != NULL) {
            diff = SeqLocAinB (slp, sfp->location);
            if (diff >= 0) {
              if (diff < *mindiff) {
                bsp = BioseqFind (SeqLocId (sfp->product));
                if (bsp != NULL) {
                  *mindiff = diff;
                  *msep = SeqMgrGetSeqEntryForData (bsp);
                }
              }
            }
          }
        }
        sfp = sfp->next;
      }
    }
    sap = sap->next;
  }
  if (bssp != NULL) {
    for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
      FindRnaByLocationOverlap (sep, slp, mindiff, msep);
    }
  }
}

static void FuseNucProtBiosources (SeqEntryPtr sep)

{
  BioSourcePtr  biop1, biop2;
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  ValNodePtr    PNTR prev;
  ValNodePtr    sdp1, sdp2;
  SeqEntryPtr   tmp;

  if (sep == NULL) return;
  if (! IS_Bioseq_set (sep)) return;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp == NULL || bssp->_class != BioseqseqSet_class_nuc_prot) return;
  tmp = FindNucSeqEntry (sep);
  if (tmp == NULL) return;
  if (! IS_Bioseq (tmp)) return;
  bsp = (BioseqPtr) tmp->data.ptrvalue;
  if (bsp == NULL) return;
  prev = &(bssp->descr);
  sdp1 = bssp->descr;
  while (sdp1 != NULL && sdp1->choice != Seq_descr_source) {
    prev = &(sdp1->next);
    sdp1 = sdp1->next;
  }
  if (sdp1 == NULL) return;
  sdp2 = SeqEntryGetSeqDescr (tmp, Seq_descr_source, NULL);
  if (sdp2 == NULL) return;
  biop1 = (BioSourcePtr) sdp1->data.ptrvalue;
  biop2 = (BioSourcePtr) sdp2->data.ptrvalue;
  if (CmpOrgById (biop1, biop2)) {
    *prev = sdp1->next;
    sdp1->next = NULL;
    SeqDescrFree (sdp1);
  }
}

static Pointer FastaSequencesFormToSeqEntryPtr (ForM f)

{
  Int2              ambig;
  Int2              annotType;
  Uint1             biomol;
  BioSourcePtr      biop;
  BioseqPtr         bsp;
  BioseqSetPtr      bssp;
  Int2              code;
  Int2              count;
  SeqAnnotPtr       curr;
  DatePtr           dp;
  FastaPagePtr      fpp;
  GatherScope       gs;
  IdListPtr         head;
  SeqEntryPtr       list;
  Int4              mindiff;
  MolInfoPtr        mip;
  MonitorPtr        mon;
  SeqEntryPtr       msep;
  SeqEntryPtr       next;
  SeqEntryPtr       nsep;
  BioseqPtr         nucbsp;
  SeqEntryPtr       nucsep;
  ValNodePtr        nulvnp;
  ValNodePtr        nxtvnp;
  Boolean           partialN;
  Boolean           partialC;
  Boolean           partialmRNA5;
  Boolean           partialmRNA3;
  BioseqPtr         protbsp;
  CharPtr           plural;
  SeqAnnotPtr       sap;
  SeqAnnotPtr PNTR  sapp;
  BioseqPtr         segseq;
  SeqEntryPtr       sep;
  SeqIdPtr          sip;
  SeqLocPtr         slp;
  SequencesFormPtr  sqfp;
  Char              str [128];
  Char              tmp [128];
  UIEnum            val;
  ValNodePtr        vnp;
  ObjectIdPtr       oip;
  UserFieldPtr      ufp;
  UserObjectPtr     uop;

  sep = NULL;
  sqfp = (SequencesFormPtr) GetObjectExtra (f);
  if (sqfp != NULL) {
    WatchCursor ();
    Update ();
    head = NULL;
    fpp = (FastaPagePtr) GetObjectExtra (sqfp->dnaseq);
    if (fpp != NULL) {
      ResolveCollidingIDs (&head, fpp->list);
    }
    fpp = (FastaPagePtr) GetObjectExtra (sqfp->protseq);
    if (fpp != NULL) {
      ResolveCollidingIDs (&head, fpp->list);
    }
    fpp = (FastaPagePtr) GetObjectExtra (sqfp->mrnaseq);
    if (fpp != NULL) {
      ResolveCollidingIDs (&head, fpp->list);
    }
    FreeTree (&head);
    biop = (BioSourcePtr) DialogToPointer (sqfp->genbio);
    if (biop != NULL) {
      code = BioSourceToGeneticCode (biop);
      if (code == 0) {
        code = gcIndexToId [GetValue (sqfp->gencode)];
      }
    } else if (sqfp->seqPackage == SEQ_PKG_PHYLOGENETIC) {
      code = gcIndexToId [GetValue (sqfp->gencode)];
    } else {
      code = 1;
    }
    list = NULL;
    fpp = (FastaPagePtr) GetObjectExtra (sqfp->dnaseq);
    if (fpp != NULL) {
      list = fpp->list;
      fpp->list = NULL;
    }
    if (sqfp->seqPackage >= SEQ_PKG_GENOMICCDNA) {
      bssp = BioseqSetNew ();
      if (bssp != NULL) {
        switch (sqfp->seqPackage) {
          case SEQ_PKG_GENOMICCDNA :
            bssp->_class = BioseqseqSet_class_gen_prod_set;
            break;
          case SEQ_PKG_POPULATION :
            bssp->_class = 14;
            break;
          case SEQ_PKG_PHYLOGENETIC :
            bssp->_class = 15;
            break;
          case SEQ_PKG_MUTATION :
            bssp->_class = 13;
            break;
          case SEQ_PKG_ENVIRONMENT :
            bssp->_class = 16;
            break;
          case SEQ_PKG_GENBANK :
            bssp->_class = 7;
            break;
          default :
            bssp->_class = 7;
            break;
        }
        sep = SeqEntryNew ();
        if (sep != NULL) {
          sep->choice = 2;
          sep->data.ptrvalue = (Pointer) bssp;
        }
      }
    }
    biomol = 0;
    if (GetEnumPopup (sqfp->moltypePopup, sqfp->moltypeAlist, &val)) {
      biomol = (Uint1) val;
    }
    if (list != NULL) {
      if (sqfp->seqPackage >= SEQ_PKG_SEGMENTED) {
        if (list->next == NULL) {
          OnlyOneComponentWarning (sqfp);
        }
      }
      while (list != NULL) {
        next = list->next;
        list->next = NULL;
        if (sep != NULL) {
          AddSeqEntryToSeqEntry (sep, list, TRUE);
          AutomaticNucleotideProcess (sqfp, list, list, biop);
        } else {
          sep = list;
          AutomaticNucleotideProcess (sqfp, list, list, biop);
        }
        if (sqfp->seqPackage > SEQ_PKG_SEGMENTED) {
          PutMolInfoOnSeqEntry (sqfp, list, biomol);
        }
        list = next;
      }
    }
    if (sep != NULL) {
      sqfp->dnamolfrommolinfo = 0;
      if (sqfp->seqPackage < SEQ_PKG_GENOMICCDNA) {
        PutMolInfoOnSeqEntry (sqfp, sep, biomol);
        if (sqfp->dnamolfrommolinfo > 0 ||
            (GetEnumPopup (sqfp->topologyPopup, topology_nuc_alist, &val) && val > 0)) {
          MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
          gs.seglevels = 1;
          gs.scope = sep;
          MemSet ((Pointer) (gs.ignore), (int) (TRUE), (size_t) (OBJ_MAX * sizeof (Boolean)));
          gs.ignore[OBJ_BIOSEQ] = FALSE;
          gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
          gs.ignore[OBJ_SEQANNOT] = FALSE;
          gs.ignore[OBJ_SEQDESC] = FALSE;
          GatherSeqEntry (sep, (Pointer) sqfp, UpdateBspMolGatherFunc, &gs);
        }
      } else {
      }
      dp = DateCurr ();
      if (dp != NULL) {
        vnp = CreateNewDescriptor (sep, Seq_descr_create_date);
        if (vnp != NULL) {
          vnp->data.ptrvalue = (Pointer) dp;
        }
      }
    }
    if (sep != NULL && sqfp->seqPackage == SEQ_PKG_SEGMENTED) {
      nsep = FindNucSeqEntry (sep);
      if (nsep != NULL && nsep->choice == 1) {
        segseq = (BioseqPtr) nsep->data.ptrvalue;
        if (segseq != NULL && segseq->repr == Seq_repr_seg && segseq->seq_ext_type == 1) {
          vnp = (ValNodePtr) segseq->seq_ext;
          while (vnp != NULL) {
            nxtvnp = vnp->next;
            if (nxtvnp != NULL && vnp->choice != SEQLOC_NULL) {
              nulvnp = ValNodeNew (NULL);
              if (nulvnp != NULL) {
                nulvnp->choice = SEQLOC_NULL;
                nulvnp->next = nxtvnp;
                vnp->next = nulvnp;
              }
            }
            vnp = nxtvnp;
          }
        }
      }
    }
    if (sqfp->seqPackage >= SEQ_PKG_POPULATION &&
        sqfp->seqPackage <= SEQ_PKG_GENBANK) {
      if (! TextHasNoText (sqfp->defline)) {
        ApplyAnnotationToAll (ADD_TITLE, sep, sqfp->partialLft, sqfp->partialRgt,
                              NULL, NULL, NULL, NULL, NULL, sqfp->defline);
      }
      if (GetStatus (sqfp->orgPrefix)) {
        PrefixOrgToDefline (sep);
      }
    }

    if (sep != NULL && sqfp->seqPackage == SEQ_PKG_GENOMICCDNA) {
      list = NULL;
      fpp = (FastaPagePtr) GetObjectExtra (sqfp->mrnaseq);
      if (fpp != NULL) {
        list = fpp->list;
        /* now we will keep instantiated mrna bioseqs */
        fpp->list = NULL;
      }
      /*
      if (list != NULL) {
        nucsep = FindNucSeqEntry (sep);
        if (nucsep != NULL) {
          while (list != NULL) {
            next = list->next;
            AutomaticMrnaProcess (nucsep, list, FALSE, FALSE);
            list = next;
          }
        }
      }
      */
      if (list != NULL) {
        nucsep = FindNucSeqEntry (sep);
        if (nucsep != NULL) {
          partialmRNA5 = GetStatus (sqfp->partialmRNA5);
          partialmRNA3 = GetStatus (sqfp->partialmRNA3);
          while (list != NULL) {
            next = list->next;
            list->next = NULL;
            AddSeqEntryToSeqEntry (sep, list, TRUE);
            AutomaticMrnaProcess (nucsep, list, partialmRNA5, partialmRNA3);
            list = next;
          }
        }
      }
    }

    list = NULL;
    fpp = (FastaPagePtr) GetObjectExtra (sqfp->protseq);
    if (fpp != NULL) {
      list = fpp->list;
      fpp->list = NULL;
    }
    if (list != NULL) {
      nucbsp = NULL;
      nucsep = FindNucSeqEntry (sep);
      if (nucsep != NULL && IS_Bioseq (nucsep)) {
        nucbsp = (BioseqPtr) nucsep->data.ptrvalue;
      }
      if (nucbsp != NULL) {
        SetBatchSuggestNucleotide (nucbsp, code);
      }
      mon = MonitorStrNewEx ("Predicting Coding Region", 20, FALSE);
      count = 0;
      while (list != NULL) {
        next = list->next;
        list->next = NULL;
        count++;
        if (mon != NULL) {
          str [0] = '\0';
          tmp [0] = '\0';
          bsp = (BioseqPtr) list->data.ptrvalue;
          if (bsp != NULL) {
            sip = SeqIdFindWorst (bsp->id);
            SeqIdWrite (sip, tmp, PRINTID_REPORT, sizeof (tmp));
          }
          sprintf (str, "Processing sequence %d [%s]", (int) count, tmp);
          MonitorStrValue (mon, str);
          Update ();
        }
        mip = MolInfoNew ();
        if (mip != NULL) {
          mip->biomol = 8;
          if (GetStatus (sqfp->protTechBoth)) {
            mip->tech = 10;
          } else {
            mip->tech = 13;
          }
          partialN = GetStatus (sqfp->partialN);
          partialC = GetStatus (sqfp->partialC);
          if (partialN && partialC) {
            mip->completeness = 5;
          } else if (partialN) {
            mip->completeness = 3;
          } else if (partialC) {
            mip->completeness = 4;
          }
          vnp = CreateNewDescriptor (list, Seq_descr_molinfo);
          if (vnp != NULL) {
            vnp->data.ptrvalue = (Pointer) mip;
          }
        }
        if (sep != NULL) {
          msep = NULL;
          if (sqfp->seqPackage == SEQ_PKG_GENOMICCDNA) {
            ClearBatchSuggestNucleotide ();
            msep = FindRnaByRefOnProtein (sep, list);
            if (msep == NULL) {
              msep = FindRnaByRefOnRna (sep, list);
            }
            if (msep == NULL && nucbsp != NULL && IS_Bioseq (list)) {
              protbsp = (BioseqPtr) list->data.ptrvalue;
              if (protbsp != NULL) {
                slp = PredictCodingRegion (nucbsp, protbsp, code);
                if (slp != NULL) {
                  mindiff = INT4_MAX;
                  FindRnaByLocationOverlap (sep, slp, &mindiff, &msep);
                }
                SeqLocFree (slp);
              }
            }
          }
          if (msep != NULL) {
            msep = GetBestTopParentForDataEx (ObjMgrGetEntityIDForChoice (msep),
                                              (BioseqPtr) msep->data.ptrvalue, TRUE);
          }
          if (msep == NULL) {
            msep = sep;
          }
          AddSeqEntryToSeqEntry (msep, list, TRUE);
          AutomaticProteinProcess (msep, list, code, sqfp->makeMRNA, NULL);
        } else {
          sep = list;
          AutomaticProteinProcess (sep, list, code, sqfp->makeMRNA, NULL);
        }
        list = next;
      }
      if (nucbsp != NULL) {
        ClearBatchSuggestNucleotide ();
      }
      mon = MonitorFree (mon);
      Update ();
    }
    if (biop != NULL) {
      vnp = CreateNewDescriptor (sep, Seq_descr_source);
      if (vnp != NULL) {
        vnp->data.ptrvalue = (Pointer) biop;
      }
    }
    if (sqfp->seqPackage >= SEQ_PKG_POPULATION &&
        sqfp->seqPackage <= SEQ_PKG_GENBANK) {
        if (GetStatus (sqfp->makeAlign)) {
        sap = SeqEntryToSeqAlign (sep, Seq_mol_na);
        if (sap != NULL && sap->type == 2) {

          oip = ObjectIdNew ();
          oip->str = StringSave ("Hist Seqalign");
          uop = UserObjectNew ();
          uop->type = oip;

          oip = ObjectIdNew();
          oip->str = StringSave ("Hist Seqalign");
          ufp = UserFieldNew ();
          ufp->choice = 4;
          ufp->data.boolvalue = TRUE;
          ufp->label = oip;

          uop->data = ufp;
	      ValNodeAddPointer (&(sap->desc), Annot_descr_user, (Pointer) uop);

          sapp = NULL;
          if (IS_Bioseq (sep)) {
            bsp = (BioseqPtr) sep->data.ptrvalue;
            sapp = &(bsp->annot);
          } else if (IS_Bioseq_set (sep)) {
            bssp = (BioseqSetPtr) sep->data.ptrvalue;
            sapp = &(bssp->annot);
          }
          if (sapp != NULL) {
            if (*sapp != NULL) {
              curr = *sapp;
              while (curr->next != NULL) {
                curr = curr->next;
              }
              curr->next = sap;
            } else {
              *sapp = sap;
            }
          }
        }
      }
      annotType = GetValue (sqfp->annotType);
      if (annotType > 0) {
        switch (annotType) {
          case 1 :
            ApplyAnnotationToAll (ADD_IMP, sep, sqfp->partialLft, sqfp->partialRgt,
                                  sqfp->geneName, NULL, NULL, NULL,
                                  sqfp->featcomment, NULL);
            break;
          case 2 :
            ApplyAnnotationToAll (ADD_RRNA, sep, sqfp->partialLft, sqfp->partialRgt,
                                  sqfp->geneName, NULL, NULL, sqfp->protOrRnaName,
                                  sqfp->featcomment, NULL);
            break;
          case 3 :
            ambig = ApplyAnnotationToAll (ADD_CDS, sep, sqfp->partialLft, sqfp->partialRgt,
                                          sqfp->geneName, sqfp->protOrRnaName, sqfp->protDesc, NULL,
                                          sqfp->featcomment, NULL);
            if (ambig > 0) {
              if (ambig > 1) {
                 plural = "records";
               } else {
                 plural = "record";
               }
               Message (MSG_OK, "Possible ambiguous frames detected in %d %s",
                        (int) ambig, plural);
            }
            break;
          default :
            break;
        }
      }
    }
    ArrowCursor ();
    Update ();
  }
  FuseNucProtBiosources (sep);
  return (Pointer) sep;
}

static ValNodePtr FastaTestSequencesForm (ForM f)

{
  BioSourcePtr      biop;
  BioseqPtr         bsp;
  Int2              code;
  FastaPagePtr      fpp;
  ValNodePtr        head;
  ObjectIdPtr       oip;
  CharPtr           ptr;
  SeqEntryPtr       sep;
  SeqIdPtr          sip;
  SequencesFormPtr  sqfp;
  Char              str [128];
  CharPtr           title;

  head = NULL;

  sqfp = (SequencesFormPtr) GetObjectExtra (f);
  if (sqfp != NULL) {

    if (sqfp->seqPackage != SEQ_PKG_PHYLOGENETIC) {
      head = TestDialog (sqfp->genbio);
      if (head != NULL) {
        head->choice = 1;
      } else {
        if (sqfp->seqPackage <= SEQ_PKG_GENOMICCDNA) {
          biop = (BioSourcePtr) DialogToPointer (sqfp->genbio);
          code = BioSourceToGeneticCode (biop);
          BioSourceFree (biop);
          if (code == 0) {
            head = AddStringToValNodeChain (head, "genetic code", 1);
          }
        }
      }
    } else {
      /*
      head = AddStringToValNodeChain (head, "organism name (in FASTA definition line)", 1);
      */
    }

    fpp = (FastaPagePtr) GetObjectExtra (sqfp->dnaseq);
    if (fpp != NULL) {
      if (fpp->list == NULL) {
        head = AddStringToValNodeChain (head, "nucleotide sequence", 1);
      } else if (head != NULL) {
        sep = fpp->list;
        if (sep != NULL) {
          title = NULL;
          SeqEntryExplore (sep, (Pointer) (&title), FindFirstTitle);
          if (title != NULL) {
            ptr = StringISearch (title, "[org=");
            if (ptr != NULL) {
              StringNCpy_0 (str, ptr + 5, sizeof (str));
              ptr = StringChr (str, ']');
              if (ptr != NULL) {
                *ptr = '\0';
                head = ValNodeFreeData (head);
              }
            } else {
              ptr = StringISearch (title, "[organism=");
              if (ptr != NULL) {
                StringNCpy_0 (str, ptr + 10, sizeof (str));
                ptr = StringChr (str, ']');
                if (ptr != NULL) {
                  *ptr = '\0';
                  head = ValNodeFreeData (head);
                }
              }
            }
          }
        }
      }
      sep = fpp->list;
      if (sep != NULL && sqfp->seqPackage == SEQ_PKG_SINGLE && IS_Bioseq (sep)) {
        bsp = (BioseqPtr) sep->data.ptrvalue;
        if (bsp != NULL) {
          sip = bsp->id;
          if (sip != NULL && sip->choice == SEQID_LOCAL) {
            oip = (ObjectIdPtr) sip->data.ptrvalue;
            if (oip != NULL && oip->str != NULL) {
              if (StringICmp (oip->str, "nuc 1") == 0) {
                GetTitle (sqfp->singleSeqID, str, sizeof (str));
                if (! StringHasNoText (str)) {
                  oip->str = MemFree (oip->str);
                  oip->str = StringSaveNoNull (str);
                  SeqMgrReplaceInBioseqIndex (bsp);
                }
                if (StringICmp (oip->str, "nuc 1") == 0) {
                  head = AddStringToValNodeChain (head, "unique identifier for your nucleotide sequence", 1);
                  SafeShow (sqfp->singleIdGrp);
                }
              }
            }
          }
        }
      }
    }

  }
  return head;
}

static void LaunchSequinQuickGuide (void)

{
  Char       str [256];
#ifdef WIN_MOTIF
  NS_Window  window = NULL;
#endif

  sprintf (str,
           "http://www.ncbi.nlm.nih.gov/Sequin/QuickGuide/sequin.htm#before");
#ifdef WIN_MAC
  Nlm_SendURLAppleEvent (str, "MOSS", NULL);
#endif
#ifdef WIN_MSWIN
  Nlm_MSWin_OpenDocument (str);
#endif
#ifdef WIN_MOTIF
  NS_OpenURL (&window, str, NULL, TRUE);
  NS_WindowFree (window);
#endif
}

extern Boolean allowUnableToProcessMessage;

static CharPtr noOrgInTitleAbort =
"sequences have organism information in titles. " \
"It is critical to annotate the data file with organism and source information. " \
"Sequin will not continue processing this submission. " \
"Please read the Sequin Quick Guide section on preparing the data files before proceeding. " \
"Do you wish to launch your browser on the Sequin Quick Guide automatically?";

static CharPtr pleaseReadLocalGuide =
"Please read your local copy of the Sequin Quick Guide before annotating your data file.";

static CharPtr noSrcInTitleAbort =
"sequences have source information in titles. " \
"It is critical to annotate the data file with organism and source information. " \
"Sequin will continue processing this submission. " \
"However, please consider reading the Sequin Quick Guide section on preparing the data files before proceeding.";

static Pointer PhylipSequencesFormToSeqEntryPtr (ForM f)

{
  Int2              ambig;
  Int2              annotType;
  MsgAnswer         ans;
  Uint2             biomol;
  BioSourcePtr      biop;
  BioseqSetPtr      bssp;
  Int2              code;
  DatePtr           dp;
  CharPtr           plural;
  PhylipPagePtr     ppp;
  SeqEntryPtr       sep;
  Int2              seqtitles;
  Int2              seqtotals;
  Char              str [256];
  SequencesFormPtr  sqfp;
  SeqEntryPtr       tmp;
  CharPtr           ttl;
  UIEnum            val;
  ValNodePtr        vnp;

  sep = NULL;
  sqfp = (SequencesFormPtr) GetObjectExtra (f);
  if (sqfp != NULL) {
    biop = (BioSourcePtr) DialogToPointer (sqfp->genbio);
    if (biop != NULL) {
      code = BioSourceToGeneticCode (biop);
      if (code == 0) {
        code = gcIndexToId [GetValue (sqfp->gencode)];
      }
    } else if (sqfp->seqPackage == SEQ_PKG_PHYLOGENETIC) {
      code = gcIndexToId [GetValue (sqfp->gencode)];
    } else {
      code = 1;
    }
    ppp = (PhylipPagePtr) GetObjectExtra (sqfp->dnaseq);
    if (ppp != NULL) {
      sep = ppp->sep;
      ppp->sep = NULL;
    }
    if (sep != NULL) {
      if (biop != NULL) {
        vnp = CreateNewDescriptor (sep, Seq_descr_source);
        if (vnp != NULL) {
          vnp->data.ptrvalue = (Pointer) biop;
        }
      }
      biomol = 0;
      if (GetEnumPopup (sqfp->moltypePopup, sqfp->moltypeAlist, &val)) {
        biomol = (Uint1) val;
      }
      if (IS_Bioseq_set (sep)) {
        bssp = (BioseqSetPtr) sep->data.ptrvalue;
        if (bssp != NULL && (bssp->_class == 7 ||
                             (IsPopPhyEtcSet (bssp->_class)))) {
          seqtitles = 0;
          seqtotals = 0;
          for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
            /*
            ttl = SeqEntryGetTitle (tmp);
            */
            ttl = NULL;
            SeqEntryExplore (sep, (Pointer) (&ttl), FindFirstTitle);
            if (ttl != NULL) {
              if (sqfp->seqPackage == SEQ_PKG_PHYLOGENETIC) {
                if (StringISearch (ttl, "[org=") != NULL ||
                    StringISearch (ttl, "[organism=") != NULL) {
                  seqtitles++;
                }
              } else if (StringISearch (ttl, "[") != NULL) {
                seqtitles++;
              }
            }
            seqtotals++;
          }
          if (seqtotals != seqtitles) {
            sprintf (str, "None");
            if (seqtitles > 0) {
              sprintf (str, "Only %d", (int) seqtitles);
            }
            ArrowCursor ();
            Update ();
            Beep ();
            if (! indexerVersion) {
              if (sqfp->seqPackage == SEQ_PKG_PHYLOGENETIC) {
                ans = Message (MSG_YN, "%s of %d %s", str, (int) seqtotals, noOrgInTitleAbort);
                if (ans == ANS_YES) {
                  LaunchSequinQuickGuide ();
                } else {
                  Message (MSG_OK, "%s", pleaseReadLocalGuide);
                }
                allowUnableToProcessMessage = FALSE;
                QuitProgram ();
                return NULL; /* aborting */
              } else {
                Message (MSG_OK, "%s of %d %s", str, (int) seqtotals, noSrcInTitleAbort);
              }
            } else {
              if (sqfp->seqPackage == SEQ_PKG_PHYLOGENETIC) {
                Message (MSG_OK, "%s of %d %s (Regular version will abort here.)", str, (int) seqtotals, noOrgInTitleAbort);
              } else {
                Message (MSG_OK, "%s of %d %s (Regular version will continue here.)", str, (int) seqtotals, noSrcInTitleAbort);
              }
            }
          }
        }
        if (bssp != NULL) {
          switch (sqfp->seqPackage) {
            case SEQ_PKG_POPULATION :
              bssp->_class = 14;
              break;
            case SEQ_PKG_PHYLOGENETIC :
              bssp->_class = 15;
              break;
            case SEQ_PKG_MUTATION :
              bssp->_class = 13;
              break;
            case SEQ_PKG_ENVIRONMENT :
              bssp->_class = 16;
              break;
            case SEQ_PKG_GENBANK :
              bssp->_class = 7;
              break;
            default :
              bssp->_class = 7;
              break;
          }
          tmp = bssp->seq_set;
          if (tmp == NULL || tmp->next == NULL) {
            OnlyOneComponentWarning (sqfp);
          }
          for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
            AutomaticNucleotideProcess (sqfp, tmp, tmp, biop);
            PutMolInfoOnSeqEntry (sqfp, tmp, biomol);
          }
        }
      } else {
        OnlyOneComponentWarning (sqfp);
        PutMolInfoOnSeqEntry (sqfp, sep, biomol);
      }
      dp = DateCurr ();
      if (dp != NULL) {
        vnp = CreateNewDescriptor (sep, Seq_descr_create_date);
        if (vnp != NULL) {
          vnp->data.ptrvalue = (Pointer) dp;
        }
      }
    }
    if (sqfp->seqPackage >= SEQ_PKG_POPULATION &&
        sqfp->seqPackage <= SEQ_PKG_GENBANK) {
      if (! TextHasNoText (sqfp->defline)) {
        ApplyAnnotationToAll (ADD_TITLE, sep, sqfp->partialLft, sqfp->partialRgt,
                              NULL, NULL, NULL, NULL, NULL, sqfp->defline);
      }
      if (GetStatus (sqfp->orgPrefix)) {
        PrefixOrgToDefline (sep);
      }
    }
    annotType = GetValue (sqfp->annotType);
    if (annotType > 0) {
      switch (annotType) {
        case 1 :
          ApplyAnnotationToAll (ADD_IMP, sep, sqfp->partialLft, sqfp->partialRgt,
                                sqfp->geneName, NULL, NULL, NULL,
                                sqfp->featcomment, NULL);
          break;
        case 2 :
          ApplyAnnotationToAll (ADD_RRNA, sep, sqfp->partialLft, sqfp->partialRgt,
                                sqfp->geneName, NULL, NULL, sqfp->protOrRnaName,
                                sqfp->featcomment, NULL);
          break;
        case 3 :
          ambig = ApplyAnnotationToAll (ADD_CDS, sep, sqfp->partialLft, sqfp->partialRgt,
                                        sqfp->geneName, sqfp->protOrRnaName, sqfp->protDesc, NULL,
                                        sqfp->featcomment, NULL);
          if (ambig > 0) {
            if (ambig > 1) {
               plural = "records";
             } else {
               plural = "record";
             }
             Message (MSG_OK, "Possible ambiguous frames detected in %d %s",
                      (int) ambig, plural);
          }
          break;
        default :
          break;
      }
    }
  }
  FuseNucProtBiosources (sep);
  return (Pointer) sep;
}

static ValNodePtr PhylipTestSequencesForm (ForM f)

{
  ValNodePtr        head;
  PhylipPagePtr     ppp;
  CharPtr           ptr;
  SeqEntryPtr       sep;
  SequencesFormPtr  sqfp;
  Char              str [128];
  CharPtr           title;

  head = NULL;

  sqfp = (SequencesFormPtr) GetObjectExtra (f);
  if (sqfp != NULL) {

    if (sqfp->seqPackage != SEQ_PKG_PHYLOGENETIC) {
      head = TestDialog (sqfp->genbio);
      if (head != NULL) {
        head->choice = 1;
      }
    } else {
      /*
      head = AddStringToValNodeChain (head, "organism name (in FASTA definition line)", 1);
      */
    }

    ppp = (PhylipPagePtr) GetObjectExtra (sqfp->dnaseq);
    if (ppp != NULL) {
      if (ppp->sep == NULL) {
        head = AddStringToValNodeChain (head, "nucleotide sequence", 1);
      } else if (head != NULL) {
        sep = ppp->sep;
        if (sep != NULL) {
          title = NULL;
          SeqEntryExplore (sep, (Pointer) (&title), FindFirstTitle);
          if (title != NULL) {
            ptr = StringISearch (title, "[org=");
            if (ptr != NULL) {
              StringNCpy_0 (str, ptr + 5, sizeof (str));
              ptr = StringChr (str, ']');
              if (ptr != NULL) {
                *ptr = '\0';
                head = ValNodeFreeData (head);
              }
            } else {
              ptr = StringISearch (title, "[organism=");
              if (ptr != NULL) {
                StringNCpy_0 (str, ptr + 10, sizeof (str));
                ptr = StringChr (str, ']');
                if (ptr != NULL) {
                  *ptr = '\0';
                  head = ValNodeFreeData (head);
                }
              }
            }
          }
        }
      }
    }

  }
  return head;
}

static Boolean ImportSequencesForm (ForM f, CharPtr filename)

{
  SequencesFormPtr  sqfp;

  sqfp = (SequencesFormPtr) GetObjectExtra (f);
  if (sqfp != NULL) {
    switch (sqfp->tagFromPage [sqfp->currentPage]) {
      case NUCLEOTIDE_PAGE :
        return ImportDialog (sqfp->dnaseq, "");
      case MRNA_PAGE :
        return ImportDialog (sqfp->mrnaseq, "");
      case PROTEIN_PAGE :
        return ImportDialog (sqfp->protseq, "");
      default :
        break;
    }
  }
  return FALSE;
}

static void ImportBtnProc (ButtoN b)

{
  SequencesFormPtr  sqfp;

  sqfp = (SequencesFormPtr) GetObjectExtra (b);
  if (sqfp) {
    ImportSequencesForm (sqfp->form, "");
  }
}

static void SetOrgNucProtImportExportItems (SequencesFormPtr sqfp)

{
  IteM  exportItm;
  IteM  importItm;

  if (sqfp != NULL) {
    importItm = FindFormMenuItem ((BaseFormPtr) sqfp, VIB_MSG_IMPORT);
    exportItm = FindFormMenuItem ((BaseFormPtr) sqfp, VIB_MSG_EXPORT);
    switch (sqfp->tagFromPage [sqfp->currentPage]) {
      case ORGANISM_PAGE :
        SafeSetTitle (importItm, "Import...");
        SafeSetTitle (exportItm, "Export...");
        SafeDisable (importItm);
        SafeDisable (exportItm);
        break;
      case NUCLEOTIDE_PAGE :
        switch (sqfp->seqFormat) {
          case SEQ_FMT_FASTA :
            if (sqfp->seqPackage == SEQ_PKG_GENOMICCDNA) {
              SafeSetTitle (importItm, "Import Genomic FASTA...");
            } else {
              SafeSetTitle (importItm, "Import Nucleotide FASTA...");
            }
            break;
          case SEQ_FMT_ALIGNMENT :
            SafeSetTitle (importItm, "Import Nucleotide Alignment...");
            break;
          default :
            SafeSetTitle (importItm, "Import Nucleotide FASTA...");
            break;
        }
        SafeSetTitle (exportItm, "Export...");
        SafeEnable (importItm);
        SafeDisable (exportItm);
        break;
      case MRNA_PAGE :
        SafeSetTitle (importItm, "Import Transcript FASTA...");
        SafeSetTitle (exportItm, "Export...");
        SafeEnable (importItm);
        SafeDisable (exportItm);
        break;
      case PROTEIN_PAGE :
        SafeSetTitle (importItm, "Import Protein FASTA...");
        SafeSetTitle (exportItm, "Export...");
        SafeEnable (importItm);
        SafeDisable (exportItm);
        break;
      case ANNOTATE_PAGE :
        SafeSetTitle (importItm, "Import...");
        SafeSetTitle (exportItm, "Export...");
        SafeDisable (importItm);
        SafeDisable (exportItm);
        break;
      default :
        break;
    }
  }
}

static void ChangeSequencesPage (VoidPtr data, Int2 newval, Int2 oldval)

{
  SequencesFormPtr  sqfp;

  sqfp = (SequencesFormPtr) data;
  if (sqfp != NULL) {
    sqfp->currentPage = newval;
    SafeHide (sqfp->pages [oldval]);
    Update ();
    switch (sqfp->tagFromPage [newval]) {
      case ORGANISM_PAGE :
        SendMessageToDialog (sqfp->genbio, VIB_MSG_ENTER);
        break;
      case NUCLEOTIDE_PAGE :
        SendMessageToDialog (sqfp->dnaseq, VIB_MSG_ENTER);
        break;
      case MRNA_PAGE :
        SendMessageToDialog (sqfp->mrnaseq, VIB_MSG_ENTER);
        break;
      case PROTEIN_PAGE :
        SendMessageToDialog (sqfp->protseq, VIB_MSG_ENTER);
        break;
      default :
        break;
    }
    if (newval == 0) {
      SafeSetTitle (sqfp->prevBtn, "<< Prev Form");
    } else {
      SafeSetTitle (sqfp->prevBtn, "<< Prev Page");
    }
    if (newval == sqfp->numPages - 1) {
      SafeSetTitle (sqfp->nextBtn, "Next Form >>");
    } else {
      SafeSetTitle (sqfp->nextBtn, "Next Page >>");
    }
    SetOrgNucProtImportExportItems (sqfp);
    SafeShow (sqfp->pages [newval]);
    Update ();
    switch (sqfp->tagFromPage [newval]) {
      case ORGANISM_PAGE :
        SendHelpScrollMessage (helpForm, "Organism and Sequences Form", "Organism Page");
        break;
      case NUCLEOTIDE_PAGE :
        if (sqfp->seqPackage == SEQ_PKG_GENOMICCDNA) {
          SendHelpScrollMessage (helpForm, "Organism and Sequences Form", "Genomic Page");
        } else {
          SendHelpScrollMessage (helpForm, "Organism and Sequences Form", "Nucleotide Page");
        }
        break;
      case MRNA_PAGE :
        SendHelpScrollMessage (helpForm, "Organism and Sequences Form", "Transcripts Page");
        break;
      case PROTEIN_PAGE :
        SendHelpScrollMessage (helpForm, "Organism and Sequences Form", "Proteins Page");
        break;
      case ANNOTATE_PAGE :
        SendHelpScrollMessage (helpForm, "Organism and Sequences Form", "Annotation Page");
        break;
      default :
        break;
    }
  }
}

static void NextSequencesFormBtn (ButtoN b)

{
  SequencesFormPtr  sqfp;

  sqfp = (SequencesFormPtr) GetObjectExtra (b);
  if (sqfp != NULL) {
    if (sqfp->currentPage < sqfp->numPages - 1) {
      SetValue (sqfp->tbs, sqfp->currentPage + 1);
    } else if (sqfp->goToNext != NULL) {
      (sqfp->goToNext) (b);
    }
  }
}

static void PrevSequencesFormBtn (ButtoN b)

{
  SequencesFormPtr  sqfp;

  sqfp = (SequencesFormPtr) GetObjectExtra (b);
  if (sqfp != NULL) {
    if (sqfp->currentPage > 0) {
      SetValue (sqfp->tbs, sqfp->currentPage - 1);
    } else if (sqfp->goToPrev != NULL) {
      (sqfp->goToPrev) (b);
    }
  }
}

static void SequencesFormDeleteProc (Pointer formDataPtr)

{
  FastaPagePtr      fpp;
  PhylipPagePtr     ppp;
  SequencesFormPtr  sqfp;

  sqfp = (SequencesFormPtr) formDataPtr;
  if (sqfp != NULL) {
    switch (sqfp->tagFromPage [sqfp->currentPage]) {
      case ORGANISM_PAGE :
        ClearText (CurrentVisibleText ());
        break;
      case NUCLEOTIDE_PAGE :
        if (sqfp->seqFormat == SEQ_FMT_FASTA) {
          fpp = (FastaPagePtr) GetObjectExtra (sqfp->dnaseq);
          if (fpp != NULL) {
            ResetFastaPage (fpp);
            fpp->path [0] = '\0';
            SafeHide (fpp->doc);
            Reset (fpp->doc);
            SafeShow (fpp->instructions);
            Update ();
          }
          Enable (fpp->import_btn);
          Enable (fpp->create_btn);
          Disable (fpp->clear_btn);
        } else if (sqfp->seqFormat == SEQ_FMT_ALIGNMENT) {
          ppp = (PhylipPagePtr) GetObjectExtra (sqfp->dnaseq);
          if (ppp != NULL) {
            ResetPhylipPage (ppp);
            ppp->path [0] = '\0';
            SetPhylipDocInstructions (ppp);
          }
        }
        break;
      case MRNA_PAGE :
        if (sqfp->seqFormat == SEQ_FMT_FASTA) {
          fpp = (FastaPagePtr) GetObjectExtra (sqfp->mrnaseq);
          if (fpp != NULL) {
            ResetFastaPage (fpp);
            fpp->path [0] = '\0';
            SafeHide (fpp->doc);
            Reset (fpp->doc);
            SafeShow (fpp->instructions);
            Update ();
          }
        }
        break;
      case PROTEIN_PAGE :
        if (sqfp->seqPackage <= SEQ_PKG_GENOMICCDNA) {
          if (sqfp->seqFormat == SEQ_FMT_FASTA) {
            fpp = (FastaPagePtr) GetObjectExtra (sqfp->protseq);
            if (fpp != NULL) {
              ResetFastaPage (fpp);
              fpp->path [0] = '\0';
              SafeHide (fpp->doc);
              Reset (fpp->doc);
              SafeShow (fpp->instructions);
              Update ();
            }
          }
        } else {
          ClearText (CurrentVisibleText ());
        }
        break;
      default :
        break;
    }
  }
}

static CharPtr  seqSegFormTabs [] = {
  "Organism", "Nucleotide", "Proteins", NULL
};

static CharPtr  cdnaGenFormTabs [] = {
  "Organism", "Genomic", "Transcripts", "Proteins", NULL
};

static CharPtr  popPhyMutFormTabs [] = {
  "Organism", "Nucleotide", "Annotation", NULL
};

static void PasteIntoDialog (DialoG seq)

{
  Char     ch;
  FILE     *fp;
  Char     path [PATH_MAX];
  CharPtr  ptr;
  CharPtr  str;

  if (Nlm_ClipboardHasString ()) {
    TmpNam (path);
    fp = FileOpen (path, "w");
    if (fp == NULL) return;
    str = ClipboardToString ();
    if (str != NULL) {
      ptr = str;
      ch = *ptr;
      while (ch != '\0') {
        if (ch == '\r') {
          *ptr = '\n';
        }
        ptr++;
        ch = *ptr;
      }
      FilePuts (str, fp);
      MemFree (str);
    }
    FileClose (fp);
    ImportDialog (seq, path);
    FileRemove (path);
  }
}

static void SequencesFormMessage (ForM f, Int2 mssg)

{
  SequencesFormPtr  sqfp;

  sqfp = (SequencesFormPtr) GetObjectExtra (f);
  if (sqfp != NULL) {
    switch (mssg) {
      case VIB_MSG_IMPORT :
        ImportSequencesForm (f, NULL);
        break;
      case VIB_MSG_CUT :
        StdCutTextProc (NULL);
        break;
      case VIB_MSG_COPY :
        StdCopyTextProc (NULL);
        break;
      case VIB_MSG_PASTE :
        switch (sqfp->tagFromPage [sqfp->currentPage]) {
          case ORGANISM_PAGE :
            StdPasteTextProc (NULL);
            break;
          case NUCLEOTIDE_PAGE :
            PasteIntoDialog (sqfp->dnaseq);
            break;
          case MRNA_PAGE :
            PasteIntoDialog (sqfp->mrnaseq);
            break;
          case PROTEIN_PAGE :
            PasteIntoDialog (sqfp->protseq);
            break;
          default :
            StdPasteTextProc (NULL);
            break;
        }
        break;
      case VIB_MSG_DELETE :
        SequencesFormDeleteProc (sqfp);
        break;
      default :
        if (sqfp->appmessage != NULL) {
          sqfp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static void InitOrgNucProtFormActivate (WindoW w)

{
  SequencesFormPtr  sqfp;

  sqfp = (SequencesFormPtr) GetObjectExtra (w);
  if (sqfp != NULL) {
    if (sqfp->activate != NULL) {
      sqfp->activate (w);
    }
    SetOrgNucProtImportExportItems (sqfp);
  }
}

static void ChangeDnaParse (ButtoN b)

{
  FastaPagePtr      fpp;
  SequencesFormPtr  sqfp;

  sqfp = (SequencesFormPtr) GetObjectExtra (b);
  if (sqfp != NULL) {
    fpp = (FastaPagePtr) GetObjectExtra (sqfp->dnaseq);
    if (fpp != NULL) {
      fpp->parseSeqId = GetStatus (b);
      if (fpp->parseSeqId) {
        WriteSequinAppParam ("PREFERENCES", "PARSENUCSEQID", "TRUE");
        SafeHide (sqfp->singleIdGrp);
      } else {
        WriteSequinAppParam ("PREFERENCES", "PARSENUCSEQID", "FALSE");
        SafeShow (sqfp->singleIdGrp);
      }
    }
  }
}

static void ChangeMrnaParse (ButtoN b)

{
  FastaPagePtr      fpp;
  SequencesFormPtr  sqfp;

  sqfp = (SequencesFormPtr) GetObjectExtra (b);
  if (sqfp != NULL) {
    fpp = (FastaPagePtr) GetObjectExtra (sqfp->mrnaseq);
    if (fpp != NULL) {
      fpp->parseSeqId = GetStatus (b);
      if (fpp->parseSeqId) {
        WriteSequinAppParam ("PREFERENCES", "PARSEMRNASEQID", "TRUE");
      } else {
        WriteSequinAppParam ("PREFERENCES", "PARSEMRNASEQID", "FALSE");
      }
    }
  }
}

static void ChangeProtParse (ButtoN b)

{
  FastaPagePtr      fpp;
  SequencesFormPtr  sqfp;

  sqfp = (SequencesFormPtr) GetObjectExtra (b);
  if (sqfp != NULL) {
    fpp = (FastaPagePtr) GetObjectExtra (sqfp->protseq);
    if (fpp != NULL) {
      fpp->parseSeqId = GetStatus (b);
      if (fpp->parseSeqId) {
        WriteSequinAppParam ("PREFERENCES", "PARSEPROTSEQID", "TRUE");
      } else {
        WriteSequinAppParam ("PREFERENCES", "PARSEPROTSEQID", "FALSE");
      }
    }
  }
}

static void ChangeMrnaFlag (ButtoN b)

{
  SequencesFormPtr  sqfp;

  sqfp = (SequencesFormPtr) GetObjectExtra (b);
  if (sqfp != NULL) {
    sqfp->makeMRNA = GetStatus (b);
    if (sqfp->makeMRNA) {
      WriteSequinAppParam ("PREFERENCES", "CREATEMRNA", "TRUE");
    } else {
      WriteSequinAppParam ("PREFERENCES", "CREATEMRNA", "FALSE");
    }
  }
}

static void ChangeAnnotType (GrouP g)

{
  SequencesFormPtr  sqfp;
  Int2              val;

  sqfp = (SequencesFormPtr) GetObjectExtra (g);
  if (sqfp == NULL) return;
  val = GetValue (g);
  switch (val) {
    case 1 :
      SafeHide (sqfp->protOrRnaPpt);
      SafeHide (sqfp->protOrRnaName);
      SafeHide (sqfp->protDescPpt);
      SafeHide (sqfp->protDesc);
      SafeShow (sqfp->annotGrp);
      Select (sqfp->geneName);
      break;
    case 2 :
      SafeSetTitle (sqfp->protOrRnaPpt, "rRNA Name");
      SafeShow (sqfp->protOrRnaPpt);
      SafeShow (sqfp->protOrRnaName);
      SafeHide (sqfp->protDescPpt);
      SafeHide (sqfp->protDesc);
      SafeShow (sqfp->annotGrp);
      Select (sqfp->protOrRnaName);
      break;
    case 3 :
      SafeSetTitle (sqfp->protOrRnaPpt, "Protein Name");
      SafeShow (sqfp->protOrRnaPpt);
      SafeShow (sqfp->protOrRnaName);
      SafeShow (sqfp->protDescPpt);
      SafeShow (sqfp->protDesc);
      SafeShow (sqfp->annotGrp);
      Select (sqfp->protOrRnaName);
      break;
    default :
      SafeHide (sqfp->annotGrp);
      break;
  }
  Update ();
}

static CharPtr  phyloOrgFastaMsg = "\
\nFor phylogenetic studies, you should encode the \
organism names in the individual nucleotide sequence \
FASTA definition lines. These should be of the \
following form:\n\n\
>[organism=scientific name]\n\n\
Additional information, e.g., [strain=name], \
can also be included. See the help documentation for \
full details";

static CharPtr  phyloOrgPhylipMsg = "\
\nFor phylogenetic studies, you should encode the \
organism names FASTA-like definition lines after \
the normal PHYLIP, NEXUS or MACAW file. These should \
be of the following form:\n\n\
>[organism=scientific name]\n\n\
Additional information, e.g., [strain=name], \
can also be included. See the help documentation for \
full details";

static void PopulateGeneticCodePopup (PopuP gc)

{
  Int2  i;

   if (gc != NULL) {
    PopupItem (gc, " ");
    for (i = 1; i <= numGeneticCodes; i++) {
      PopupItem (gc, gcNames [i]);
    }
  }
}

typedef struct createfastadata 
{
  FORM_MESSAGE_BLOCK
  
  ValNodePtr   sequence_list;
  FastaPagePtr fpp;
  ButtoN       use_id_from_fasta_defline;
  ButtoN       import_btn;
  ButtoN       create_btn;
  ButtoN       clear_btn;
  PopuP        sequence_selector;
  ValNodePtr   selected_qual_list;
  PopuP        qual_list;
  CharPtr      default_tax_name;
  ButtoN       prev_btn;
  ButtoN       next_btn;
  Boolean      auto_id;
  Boolean      is_segmented;
  Boolean      is_popset;
} CreateFastaData, PNTR CreateFastaPtr;

typedef struct createfastaseqdata
{
  TexT           local_id_txt;
  TexT           sequence_txt;
  CharPtr        local_id;
  CharPtr        sequence;
  CharPtr        title;
  GrouP          grp;
  
  /* organism information*/
  DoC            orglist;
  TexT           taxName;
  Int2           orgrow;
  CharPtr        tax_name;
  
  DialoG         mods;
  ListPairData   mod_values;
  CreateFastaPtr master_set; /* This points to the master set */
} CreateFastaSeqData, PNTR CreateFastaSeqPtr;

static void AddFastaSequence (CreateFastaPtr cfp);
static void RemoveFastaSequence (CreateFastaPtr cfp);
static WindoW CreateFastaWindow (CreateFastaPtr cfp);
static void CollectAllSequenceInformation (ValNodePtr list);

static CharPtr valid_iupac_characters = "atgcbdhkmnrsuvwy";

static void FreeSequenceGroup (ValNodePtr list)
{
  CreateFastaSeqPtr cfsp;
  
  if (list == NULL) return;
  FreeSequenceGroup (list->next);
  list->next = NULL;
  cfsp = (CreateFastaSeqPtr) list->data.ptrvalue;
  if (cfsp != NULL)
  {
    cfsp->local_id = MemFree (cfsp->local_id);
  	cfsp->sequence = MemFree (cfsp->sequence);
  	cfsp->tax_name = MemFree (cfsp->tax_name);
  	cfsp->mod_values.selected_names_list = ValNodeFree (cfsp->mod_values.selected_names_list);
  	cfsp->mod_values.selected_values_list = ValNodeFreeData (cfsp->mod_values.selected_values_list);
  	MemFree (cfsp);
  }
  ValNodeFree (list);
}

static void FreeCreateFasta (CreateFastaPtr cfp)
{
  if (cfp == NULL) return;
  
  cfp->selected_qual_list = ValNodeFree (cfp->selected_qual_list);
  FreeSequenceGroup (cfp->sequence_list);
  cfp->sequence_list = NULL;
  cfp->default_tax_name = MemFree (cfp->default_tax_name);
  MemFree (cfp);
}

static CharPtr CreateNewLocalID (ValNodePtr list, CharPtr auto_id_fmt, Int4 to_remove)
{
  CharPtr           new_id;
  Int4              id_len = 15;
  long int          val;
  ValNodePtr        vnp;
  CreateFastaSeqPtr cfsp;
  Int4              new_seq_id = 1;
  Int4              seq_num;
  
  if (auto_id_fmt != NULL)
  {
    id_len += StringLen (auto_id_fmt);
  }
  
  new_id = (CharPtr) MemNew (id_len * sizeof (Char));
  if (new_id == NULL) return NULL;
  
  if (list == NULL)
  {
    sprintf (new_id, auto_id_fmt, 1);	
  }
  else
  {
    for (seq_num = 1, vnp = list; vnp != NULL; vnp = vnp->next, seq_num++)
    {
      if (seq_num == to_remove) continue;
      cfsp = (CreateFastaSeqPtr) vnp->data.ptrvalue;
      if (cfsp != NULL && !StringHasNoText (cfsp->local_id))
      {
        if (sscanf (cfsp->local_id, auto_id_fmt, &val) == 1) 
        {
          if (new_seq_id <= val)
          {
          	new_seq_id = val + 1;
          }
        }
      }
    }
    sprintf (new_id, auto_id_fmt, new_seq_id);
  }
  return new_id;
}

static ValNodePtr AddSequenceGroup 
(ValNodePtr     list, 
 ValNodePtr     qual_list,
 CharPtr        def_taxname, 
 Boolean        auto_id, 
 Int4           to_remove,
 CreateFastaPtr cfp)
{
  CreateFastaSeqPtr cfsp;
  ValNodePtr        vnp, name_vnp, val_vnp;
  CharPtr           auto_id_fmt = "seq_%d";
  
  cfsp = (CreateFastaSeqPtr) MemNew (sizeof (CreateFastaSeqData));
  if (cfsp == NULL) return NULL;
  
  cfsp->local_id_txt = NULL;
  cfsp->local_id = NULL;
  cfsp->sequence_txt = NULL;
  cfsp->sequence = NULL;
  cfsp->master_set = cfp;
  
  /* use default tax name if we have one */
  cfsp->orgrow = -1;
  if (def_taxname != NULL)
  {
  	cfsp->tax_name = MemNew (sizeof (Char) * (StringLen (def_taxname) + 1));
  	StringCpy (cfsp->tax_name, def_taxname);
  }
  
  /* add source modifiers to match existing lists */  
  cfsp->mod_values.selected_names_list = NULL;
  cfsp->mod_values.selected_values_list = NULL;
  for (vnp = qual_list; vnp != NULL; vnp = vnp->next)
  {
    name_vnp = ValNodeNew (cfsp->mod_values.selected_names_list);
    if (cfsp->mod_values.selected_names_list == NULL)
    {
      cfsp->mod_values.selected_names_list = name_vnp;
    }
    if (name_vnp != NULL)
    {
      name_vnp->choice = vnp->choice;
      name_vnp->data.ptrvalue = vnp->data.ptrvalue;
    }
    val_vnp = ValNodeNew (cfsp->mod_values.selected_values_list); 
    if (cfsp->mod_values.selected_values_list == NULL)	
    {
      cfsp->mod_values.selected_values_list = val_vnp;
    }
    if (val_vnp != NULL)
    {
      val_vnp->data.ptrvalue = NULL;
    }
  }
  
  if (auto_id)
  {
    cfsp->local_id = CreateNewLocalID (list, auto_id_fmt, to_remove);
  }
  vnp = ValNodeNew (list);
  if (vnp == NULL) return NULL;
  vnp->data.ptrvalue = cfsp;
  if (list == NULL) list = vnp;
  return list;
}

static void ShowSequenceGroup (PopuP p)
{
  CreateFastaPtr cfp;
  Int4           val;
  Int4           seq_num = 1;
  ValNodePtr     vnp;
  CreateFastaSeqPtr cfsp;
  
  cfp = (CreateFastaPtr) GetObjectExtra (p);
  if (cfp == NULL) return;
  
  val = GetValue (p);
  for (vnp = cfp->sequence_list; vnp != NULL; vnp = vnp->next, seq_num++)
  {
    cfsp = (CreateFastaSeqPtr) vnp->data.ptrvalue;
    if (cfsp == NULL) continue;
    if (seq_num == val)
    {
      Show (cfsp->grp);
      if (val > 1)
      {
      	Enable (cfp->prev_btn);
      }
      else
      {
      	Disable (cfp->prev_btn);
      }
      if (vnp->next == NULL)
      {
      	Disable (cfp->next_btn);
      }
      else
      {
      	Enable (cfp->next_btn);
      }
    }
    else
    {
      Hide (cfsp->grp);
    }
  }
  
}

static void GetMasterQualList (CreateFastaPtr cfp)
{
  ValNodePtr        vnp, name_vnp, master_vnp;
  CreateFastaSeqPtr cfsp;
  if (cfp == NULL) return;
  
  cfp->selected_qual_list = ValNodeFree (cfp->selected_qual_list);
  /* add modifiers listed in any sequence to selected list */
  for (vnp = cfp->sequence_list; vnp != NULL; vnp = vnp->next)
  {
  	cfsp = (CreateFastaSeqPtr) vnp->data.ptrvalue;
  	if (cfsp != NULL)
  	{
  	  for (name_vnp = cfsp->mod_values.selected_names_list; name_vnp != NULL; name_vnp = name_vnp->next)
  	  {
        master_vnp = FindExactStringInStrings (cfp->selected_qual_list, name_vnp->data.ptrvalue);
  	  	if (master_vnp == NULL)
  	  	{
  	  	  master_vnp = ValNodeNew (cfp->selected_qual_list);
  	  	  if (cfp->selected_qual_list == NULL)
  	  	  {
  	  	  	cfp->selected_qual_list = master_vnp;
  	  	  }
  	  	  master_vnp->choice = name_vnp->choice;
  	  	  master_vnp->data.ptrvalue = name_vnp->data.ptrvalue;
  	  	}
  	  }
  	}
  }
}

static void 
ReportMissingQualsOrValues 
(ValNodePtr missing_list,
 CharPtr    sequence_str_fmt,
 CharPtr    local_id,
 CharPtr    sequence_num_fmt,
 Int4       seq_num)
{
  ValNodePtr vnp;
  Int4       msg_len = 0;
  CharPtr    msg;
  
  if (missing_list == NULL || sequence_str_fmt == NULL || sequence_num_fmt == NULL)
  {
  	return;
  }
  
  if (StringHasNoText (local_id))
  {
  	msg_len = StringLen (sequence_str_fmt) + 15;
  }
  else
  {
  	msg_len = StringLen (sequence_num_fmt) + StringLen (local_id);
  }
  
  for (vnp = missing_list; vnp != NULL; vnp = vnp->next)
  {
  	msg_len += StringLen (vnp->data.ptrvalue) + 3;
  }

  msg = (CharPtr) MemNew (sizeof (Char) * msg_len);
  if (msg != NULL)
  {
    if (StringHasNoText (local_id))
    {
      sprintf (msg, sequence_num_fmt, seq_num);
    }
    else
    {
      sprintf (msg, sequence_str_fmt, local_id);
    }
    for (vnp = missing_list; vnp != NULL; vnp = vnp->next)
    {
      StringCat (msg, vnp->data.ptrvalue);
      if (vnp->next != NULL)
      {
        StringCat (msg, ", ");
        
      }
    }
    Message (MSG_ERROR, msg);
    MemFree (msg);
  }
}

static Boolean 
CheckSeqQualListForConsistency 
(CreateFastaSeqPtr cfsp,
 ValNodePtr        master_list,
 Int4              seq_num)
{
  ValNodePtr name_vnp, value_vnp, master_vnp, missing_vnp;
  ValNodePtr missing_list = NULL;
  Int4       num_missing = 0;
  CharPtr    sequence_str_fmt = "Sequence %s is missing the following quals: ";
  CharPtr    sequence_num_fmt = "Sequence %d is missing the following quals: ";
  Boolean    rval = TRUE;
  
  if (cfsp == NULL) return FALSE;
  
  for (master_vnp = master_list;
       master_vnp != NULL;
       master_vnp = master_vnp->next)
  {
   	if (FindExactStringInStrings (cfsp->mod_values.selected_names_list,
   	                              master_vnp->data.ptrvalue) == NULL)
   	{
   	  missing_vnp = ValNodeNew (missing_list);
   	  if (missing_list == NULL)
   	  {
   	  	missing_list = missing_vnp;
   	  }
   	  if (missing_vnp != NULL)
   	  {
   	  	missing_vnp->data.ptrvalue = master_vnp->data.ptrvalue;
   	  	missing_vnp->choice = master_vnp->choice;
   	  }
   	  num_missing ++;
   	}
  }
  
  for (name_vnp = cfsp->mod_values.selected_names_list,
       value_vnp = cfsp->mod_values.selected_values_list;
       name_vnp != NULL && value_vnp != NULL;
       name_vnp = name_vnp->next, value_vnp = value_vnp->next)
  {
  	if (StringHasNoText (value_vnp->data.ptrvalue))
  	{
  	  num_missing ++;
  	  missing_vnp = ValNodeNew (missing_list);
  	  if (missing_list == NULL)
  	  {
  	  	missing_list = missing_vnp;
  	  }
  	  if (missing_vnp != NULL)
  	  {
  	  	missing_vnp->data.ptrvalue = name_vnp->data.ptrvalue;
  	  	missing_vnp->choice = name_vnp->choice;
  	  }
  	}
  }
  
  if (num_missing > 0)
  {
    rval = FALSE;
    ReportMissingQualsOrValues (missing_list, sequence_str_fmt, cfsp->local_id,
                                sequence_num_fmt, seq_num);
    ValNodeFree (missing_list);
  }
  return rval;
}

static Boolean CheckQualListForConsistency (CreateFastaPtr cfp)
{
  ValNodePtr        vnp;
  CreateFastaSeqPtr cfsp;
  Boolean           rval = TRUE;
  Int4              seq_num = 1;
  
  if (cfp == NULL) return FALSE;

  GetMasterQualList (cfp);

  for (vnp = cfp->sequence_list; vnp != NULL; vnp = vnp->next, seq_num++)
  {
  	cfsp = (CreateFastaSeqPtr) vnp->data.ptrvalue;
  	if (cfsp != NULL)
  	{
  	  rval &= CheckSeqQualListForConsistency (cfsp, cfp->selected_qual_list, seq_num); 	  
    }
  }
  return rval;
}
  
/* code for scientific name selection controls */
static ValNodePtr orglist = NULL;

static void LoadOrganismList (void)
{
  Char     str [PATH_MAX];
  ReadBufferData    rbd;
  CharPtr  line;
  ValNodePtr vnp;
  CharPtr           ptr;
  FILE              *f;

  if (orglist != NULL) return;
  ProgramPath (str, sizeof (str));
  ptr = StringRChr (str, DIRDELIMCHR);
  if (ptr != NULL) {
    *ptr = '\0';
    FileBuildPath (str, NULL, "taxlist.txt");
    f = FileOpen (str, "r");
    if (f == NULL) {
      if (GetAppParam ("NCBI", "NCBI", "DATA", "", str, sizeof (str))) {
        FileBuildPath (str, NULL, "taxlist.txt");
        f = FileOpen (str, "r");
      }
    }
    if (f != NULL) {
      rbd.fp = f;
      rbd.current_data = NULL;
      line = AbstractReadFunction (&rbd);
      line = AbstractReadFunction (&rbd);
      while (line != NULL)
      {
        ptr = StringChr (line, '\t');
        if (ptr != NULL)
        {
          *ptr = 0;
        }
      	vnp = ValNodeNew (orglist);
      	if (vnp != NULL)
      	{
      	  vnp->data.ptrvalue = line;
      	}
      	if (orglist == NULL)
      	{
      	  orglist = vnp;
      	}
      	line = AbstractReadFunction (&rbd);
      }
      FileClose (f);
    }
  }
}

static CharPtr GetTextForOrgPos (Int4 pos)
{
  ValNodePtr vnp;
  Int4       val;
  
  for (vnp = orglist, val = 1; vnp != NULL && val < pos; vnp = vnp->next, val++)
  {
  }
  if (vnp != NULL)
  {
  	return (CharPtr)vnp->data.ptrvalue;
  }
  else
  {
  	return NULL;
  }
}

static void GetOrgPosForText (CharPtr cp, Int4Ptr pos, Boolean PNTR match)
{
  ValNodePtr vnp;
  Int4       val = 1;
  CharPtr    dat;
  Int4       res;
  
  if (cp == NULL || pos == NULL || match == NULL) return;
  for (vnp = orglist; vnp != NULL; vnp = vnp->next)
  {
  	dat = (CharPtr) vnp->data.ptrvalue;
  	res = StringCmp (cp, dat);
  	if (res < 0)
  	{
  	  *pos = val;
  	  *match = FALSE;
  	  return;
  	}
  	else if (res == 0)
  	{
  	  *pos = val;
  	  *match = TRUE;
  	  return;
  	}
  	val++;
  }
  *pos = val - 1;
  *match = FALSE;
}

static void SetAllOrganismNames (CreateFastaSeqPtr cfsp)
{
  ValNodePtr        vnp;
  CreateFastaSeqPtr cfsp_set;
  if (cfsp == NULL || cfsp->master_set == NULL) return;
  
  for (vnp = cfsp->master_set->sequence_list;
       vnp != NULL;
       vnp = vnp->next)
  {
    cfsp_set = vnp->data.ptrvalue;
    if (cfsp_set == NULL || cfsp_set == cfsp)
    {
      continue; /* don't need to set the values in this one */
    }
    SetTitle (cfsp_set->taxName, cfsp->tax_name);
    cfsp_set->orgrow = cfsp->orgrow;
  }
}

static void SetOrganismDoc (DoC d, PoinT pt)
{
  Int2 item, row, prevrow;

  CreateFastaSeqPtr cfsp;
  
  cfsp = (CreateFastaSeqPtr) GetObjectExtra (d);
  if (cfsp == NULL) return;
  
  MapDocPoint (d, pt, &item, &row, NULL, NULL);
  if (item > 0 && row > 0) {
    prevrow = cfsp->orgrow;
    cfsp->orgrow = item;
    if (item != prevrow)
    {
      if (prevrow != -1)
      {
        InvalDocRows (d, prevrow, 1, 1);
      }
      InvalDocRows (d, item, 1, 1);
      SetTitle (cfsp->taxName, GetTextForOrgPos (item));
      if (cfsp->master_set != NULL && cfsp->master_set->is_popset)
      {
        cfsp->tax_name = MemFree (cfsp->tax_name);
        cfsp->tax_name = SaveStringFromText (cfsp->taxName);
        SetAllOrganismNames (cfsp);
      }
    }  	
  }
}

static void SetOrganismText (TexT t)
{
  CreateFastaSeqPtr cfsp;
  Int4              pos, prevpos;
  Boolean           match;
  
  cfsp = (CreateFastaSeqPtr) GetObjectExtra (t);
  if (cfsp == NULL) return;
  if (cfsp->tax_name != NULL)
  {
  	MemFree (cfsp->tax_name);
  }
  cfsp->tax_name = SaveStringFromText (cfsp->taxName);
  if (cfsp->tax_name != NULL)
  {
  	cfsp->tax_name [0] = TO_UPPER (cfsp->tax_name[0]);
  }
  GetOrgPosForText (cfsp->tax_name, &pos, &match);
  SetOffset (cfsp->orglist, 0, pos - 1);
  if (pos != cfsp->orgrow)
  {
    prevpos = cfsp->orgrow;
    if (match)
    { 
      cfsp->orgrow = pos;
      SetTitle (cfsp->taxName, cfsp->tax_name);
    }
    else
    {
      cfsp->orgrow = -1;
    }
  	if (prevpos != -1)
    {
  	  InvalDocRows (cfsp->orglist, prevpos, 1, 1);
    }
    if (match)
    {
      InvalDocRows (cfsp->orglist, cfsp->orgrow, 1, 1);
    }
  }
  else if (!match)
  {
  	cfsp->orgrow = -1;
    InvalDocRows (cfsp->orglist, pos, 1, 1);	
  }
  if (cfsp->master_set != NULL && cfsp->master_set->is_popset)
  {
    cfsp->tax_name = MemFree (cfsp->tax_name);
    cfsp->tax_name = SaveStringFromText (cfsp->taxName);
    SetAllOrganismNames (cfsp);
  }
}

static Boolean OrgNameHighlight (DoC doc, Int2 item, Int2 row, Int2 col)
{
  CreateFastaSeqPtr cfsp;
  
  cfsp = (CreateFastaSeqPtr) GetObjectExtra (doc);
  if (cfsp == NULL) return FALSE;
  
  if (item == cfsp->orgrow) return TRUE;
  return FALSE;
}

static ParData orgListPar = {FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0};
static ColData orgListCol [] = {
  {0, 0, 80, 0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, FALSE},
  {0, 0,  0, 0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, FALSE},
  {0, 0,  0, 0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, FALSE},
  {0, 0,  0, 0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, FALSE},
  {0, 0,  0, 0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, FALSE},
  {0, 0,  0, 0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, FALSE},
  {0, 0,  0, 0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, FALSE},
  {0, 0,  0, 0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, TRUE}};

static void DrawOrganismSelection (GrouP g, CreateFastaSeqPtr cfsp, CreateFastaPtr cfp)
{
  GrouP           f, h;
  Int2            height;
  ValNodePtr      vnp;

  LoadOrganismList ();	
  h = NormalGroup (g, 2, 0, "Organism Information", programFont, NULL);
  f = NormalGroup (h, 0, 3, "Scientific Name", programFont, NULL);
  cfsp->taxName = DialogText (f, "", 20, SetOrganismText);
  SetObjectExtra (cfsp->taxName, cfsp, NULL);
  if (cfsp->tax_name != NULL)
  {
  	SetTitle (cfsp->taxName, cfsp->tax_name);
  	SetOrganismText (cfsp->taxName);
  }
  if (cfp->is_popset)
  {
    StaticPrompt (f, "All organism names must be identical for a population study",
                   0, dialogTextHeight, programFont, 'l');
  }
/*  f = HiddenGroup (h, 1, 0, NULL); */
  SelectFont (programFont);
  height = LineHeight ();
  SelectFont (systemFont);
  cfsp->orglist = DocumentPanel (f, stdCharWidth * 25, height * 6);
  SetObjectExtra (cfsp->orglist, cfsp, NULL);
/*  SetDocAutoAdjust (cfsp->orglist, FALSE); */
  for (vnp = orglist; vnp != NULL; vnp = vnp->next)
  {
  	AppendText (cfsp->orglist, vnp->data.ptrvalue, &orgListPar, orgListCol, programFont);
  }
  SetDocAutoAdjust (cfsp->orglist, FALSE);
  SetDocProcs (cfsp->orglist, SetOrganismDoc, NULL, NULL, NULL);
  SetDocShade (cfsp->orglist, NULL, NULL, OrgNameHighlight, NULL);
  
  f = NormalGroup (h, 1, 0, "Organism Qualifiers", programFont, NULL);
  cfsp->mods = CreateModifierTagList (f, &(cfsp->mod_values));
}

static void CreateFastaPrev (ButtoN b)
{
  CreateFastaPtr cfp;
  Int4           val;
  
  cfp = (CreateFastaPtr) GetObjectExtra (b);
  if (cfp == NULL) return;
  val = GetValue (cfp->sequence_selector);
  if (val <= 1) return;
  SetValue (cfp->sequence_selector, val - 1);
  ShowSequenceGroup (cfp->sequence_selector);
}

static void CreateFastaNext (ButtoN b)
{
  CreateFastaPtr cfp;
  Int4           val;
  Int4           num_sequences;
  
  cfp = (CreateFastaPtr) GetObjectExtra (b);
  if (cfp == NULL) return;
  num_sequences = ValNodeLen (cfp->sequence_list);
  val = GetValue (cfp->sequence_selector);
  if (val >= num_sequences) return;
  SetValue (cfp->sequence_selector, val + 1);
  ShowSequenceGroup (cfp->sequence_selector);
}

static void DrawSequenceSelector (CreateFastaPtr cfp, GrouP g)
{
  GrouP k;
  Char              str [200];
  ValNodePtr        vnp; 
  Int4              seq_num;
  CreateFastaSeqPtr cfsp;
  Int4              num_sequences;

  if (cfp == NULL) return;
  
  num_sequences = ValNodeLen (cfp->sequence_list);
  
  k = NormalGroup (g, 4, 0, "Select Sequence", programFont, NULL);

  cfp->prev_btn = PushButton (k, "<<", CreateFastaPrev);
  SetObjectExtra (cfp->prev_btn, cfp, NULL);
  Disable (cfp->prev_btn);
  cfp->sequence_selector = PopupList (k, TRUE, ShowSequenceGroup);
  SetObjectExtra (cfp->sequence_selector, cfp, NULL);
  for (vnp = cfp->sequence_list, seq_num = 1; vnp != NULL; vnp = vnp->next, seq_num++)
  {
    cfsp = (CreateFastaSeqPtr) vnp->data.ptrvalue;
    if (cfsp != NULL)
    {
      if (!StringHasNoText (cfsp->local_id)) 
      {
        PopupItem (cfp->sequence_selector, cfsp->local_id);
      }
      else
      {
        sprintf (str, "Unlabeled sequence %d", seq_num);
        PopupItem (cfp->sequence_selector, str);
      }
    }
  }
  cfp->next_btn = PushButton (k, ">>", CreateFastaNext);
  SetObjectExtra (cfp->next_btn, cfp, NULL);
  if (num_sequences < 2)
  {
  	Disable (cfp->next_btn);
  }
}
 
static void DrawSequenceGroups (CreateFastaPtr cfp, GrouP g)
{
  ValNodePtr        vnp;
  CreateFastaSeqPtr cfsp;
  GrouP             k, n;
  Char              str [200];
  Int4              seq_num = 1;
  Int4              num_sequences;
  
  if (cfp == NULL || cfp->sequence_list == NULL || g == NULL) return;
  cfp->sequence_selector = NULL;
  n = HiddenGroup (g, 0, 0, NULL);
  num_sequences = ValNodeLen (cfp->sequence_list);
  for (vnp = cfp->sequence_list; vnp != NULL; vnp = vnp->next)
  {
    cfsp = (CreateFastaSeqPtr) vnp->data.ptrvalue;
    if (cfsp != NULL)
    {
      cfsp->grp = HiddenGroup (n, 0, 3, NULL);
      k = HiddenGroup (cfsp->grp, 4, 0, NULL);
      StaticPrompt (k, "Sequence ID", 0, dialogTextHeight, programFont, 'l');
      if (cfsp->local_id == NULL)
      {
      	cfsp->local_id_txt = DialogText (k, "", 20, NULL);
      }
      else
      {
      	cfsp->local_id_txt = DialogText (k, cfsp->local_id, 20, NULL);
      }
      sprintf (str, "%d of %d", seq_num, num_sequences);
      StaticPrompt (k, str, 0, dialogTextHeight, programFont, 'l');
      if (!StringHasNoText (cfsp->title))
      {
      	StaticPrompt (k, cfsp->title, 0, dialogTextHeight, programFont, 'l');
      }
      
      DrawOrganismSelection (cfsp->grp, cfsp, cfp);
      
      k = NormalGroup (cfsp->grp, 1, 0, "Sequence Characters", programFont, NULL);
      sprintf (str, "Please enter the nucleotides for your sequence.  "
                    "You may only use the valid IUPAC characters (%s).", 
                    valid_iupac_characters);
      MultiLinePrompt (k, str, 60 * stdCharWidth, programFont);
      cfsp->sequence_txt = ScrollText (k, 60, 10, programFont, FALSE, NULL);
      if (cfsp->sequence != NULL)
      {
        SetTitle (cfsp->sequence_txt, cfsp->sequence);
      }
      seq_num++;      
    }
  }
}

static Boolean SeqCharsOk( CharPtr seq_chars, Int4 seq_num)
{
  CharPtr cp;
  Char    ch;
  Boolean at_least_one = FALSE;
  Int4    len;
  CharPtr badchars;
  
  if (seq_chars == NULL)
  {
  	Message (MSG_ERROR, "There are no sequence characters for sequence %d.  Please enter some.", seq_num);
  	return FALSE;
  }
  len = StringLen (seq_chars);
  if (len == 0)
  {
  	Message (MSG_ERROR, "There are no sequence characters for sequence %d.  Please enter some.", seq_num);
  	return FALSE;
  }
  badchars = (CharPtr) MemNew (sizeof (Char) * (len + 1));
  if (badchars == NULL) return FALSE;
  badchars[0] = 0;
  len = 0;
  for (cp = seq_chars; *cp != 0; cp++)
  {
    ch = TO_LOWER (*cp);
  	if (isspace (ch))
  	{
  	  /* space allowed */
  	}
  	else if (StringChr (valid_iupac_characters, ch) == NULL)
  	{
  	  if (StringChr (badchars, *cp) == NULL)
  	  {
  	  	badchars [len] = ch;
  	  	len++; 
  	  	badchars [len] = 0;
  	  }
  	}
  	else 
  	{
  	  at_least_one = TRUE;
  	}
  }
  if (len > 0)
  {
  	Message (MSG_ERROR, 
  	         "There were %d illegal characters were found in sequence %d: %s."
  	         "  Repeated characters are listed only once. "
  	         "  You may only have IUPAC characters in your sequence "
  	         "(%s).\n", len, seq_num, badchars, valid_iupac_characters);
  	return FALSE;
  }
  if (!at_least_one)
  {
  	Message (MSG_ERROR, "There are no sequence characters in sequence %d.  Please enter some.", seq_num);
  	return FALSE;
  }
  return TRUE;
}

static CharPtr ReformatSequenceText (CharPtr seq_text)
{
  CharPtr src, dst;
  CharPtr new_text;
  Int4    num_lines;
  Int4    len;
  Int4    line_len = 80;
  Int4    counter;

  if (StringHasNoText (seq_text))
  {
  	MemFree (seq_text);
  	return NULL;
  }
  len = StringLen (seq_text);
  num_lines = len / line_len;
  len += num_lines + 2;
  new_text = (CharPtr) MemNew (len * sizeof (Char));
  if (new_text == NULL)
  {
  	return seq_text;
  }
  dst = new_text;
  counter = 0;
  for (src = seq_text; *src != 0; src++)
  {
  	if (!isspace (*src))
  	{
  	  *dst = *src;
  	  dst++;
  	  counter++;
  	  if (counter == line_len)
  	  {
  	  	*dst = '\n';
  	  	dst++;
  	  	counter = 0;
  	  }
  	}
  }
  *dst = 0;
  MemFree (seq_text);
  return new_text;
}

static CharPtr ReformatLocalId (CharPtr local_id)
{
  CharPtr cp, new_local_id;
  
  cp = local_id;
  while (*cp == '>')
  {
    cp ++;
  }
  while (isspace (*cp))
  {
  	cp++;
  }
  new_local_id = StringSave (cp);
  cp = new_local_id;
  while (*cp != 0)
  {
    if (isspace (*cp))
    {
      *cp = '_';
    }
    cp++;
  }
  MemFree (local_id);
  return new_local_id;
}

static void CollectOneSequenceInformation (CreateFastaSeqPtr cfsp)
{
  ListPairPtr     lpp;
  
  if (cfsp == NULL) return;

  /* get sequence ID */
  if (cfsp->local_id != NULL)
  {
  	MemFree (cfsp->local_id);
  }
  cfsp->local_id = SaveStringFromText (cfsp->local_id_txt);
  cfsp->local_id = ReformatLocalId (cfsp->local_id);

  /* get sequence text */
  if (cfsp->sequence != NULL)
  {
  	MemFree (cfsp->sequence);
  }
  cfsp->sequence = SaveStringFromText (cfsp->sequence_txt);
  cfsp->sequence = ReformatSequenceText (cfsp->sequence);

  /* get organism name */
  cfsp->tax_name = SaveStringFromText (cfsp->taxName);
  
  
  /* copy in qual information */
  cfsp->mod_values.selected_names_list = ValNodeFree (cfsp->mod_values.selected_names_list);
  cfsp->mod_values.selected_values_list = ValNodeFreeData (cfsp->mod_values.selected_values_list);
  lpp = GetModifierList (cfsp->mods);
  if (lpp != NULL)
  {
  	cfsp->mod_values.selected_names_list = lpp->selected_names_list;
  	cfsp->mod_values.selected_values_list = lpp->selected_values_list;
  	MemFree (lpp);
  }

}

static void CollectAllSequenceInformation (ValNodePtr list)
{
  ValNodePtr vnp;
  
  for (vnp = list; vnp != NULL; vnp = vnp->next)
  {
  	CollectOneSequenceInformation (vnp->data.ptrvalue);
  }
}

static Boolean ValidateOneSequenceInformation (CreateFastaSeqPtr cfsp, Int4 seq_num)
{
  Boolean         rval = TRUE;
  CharPtr        localid_msg = "You must include a sequence ID for sequence %d."
        "  A sequence ID is a temporary ID which will be replaced with a unique"
        " GenBank accession number by the GenBank curators.  The sequence ID should"
        " be unique for each sequence in a record.  It could represent a clone,"
        " isolate, a laboratory designation for your specimen, or some other useful"
        " information, but this is not required."
        "  A sequence ID may not begin with a '>' character or contain spaces.";
  
  if (cfsp == NULL) return FALSE;
  
  CollectOneSequenceInformation (cfsp);
  
  /* complain if sequence ID is empty */
  if (StringHasNoText (cfsp->local_id)) {
    Message (MSG_ERROR, localid_msg, seq_num);
    rval = FALSE;
  }
   
  /* get sequence text, complain if empty or bad chars */
  if (!SeqCharsOk (cfsp->sequence, seq_num))
  {
    rval = FALSE;
  }
  /* complain if no organism name */
  if (StringHasNoText (cfsp->tax_name))
  {
  	Message (MSG_ERROR, "Sequence %d has no organism name.  Please supply one.", seq_num);
  }
    
  return rval;	
}

static Boolean CheckLocalIDs (ValNodePtr list, Int4 skip)
{
  ValNodePtr        vnp, check_vnp, ignore_list, ignore_vnp;
  Int4              seq_num, check_seq_num;
  CreateFastaSeqPtr cfsp;
  CharPtr           local_id; 
  Boolean           rval = TRUE;
  Int4              num_appearances;
  ValNodePtr        appearance_positions, pos_vnp;

  ignore_list = NULL;  
  for (vnp = list, seq_num = 1; vnp != NULL; vnp = vnp->next, seq_num++)
  {
    if (seq_num == skip) continue;
    cfsp = (CreateFastaSeqPtr) vnp->data.ptrvalue;
    local_id = cfsp->local_id;
    if (StringHasNoText (local_id)) continue;
    for (ignore_vnp = ignore_list; 
         ignore_vnp != NULL && StringCmp (local_id, ignore_vnp->data.ptrvalue) != 0;
         ignore_vnp = ignore_vnp->next)
    {
    }
    if (ignore_vnp != NULL) continue;
    num_appearances = 1;
    appearance_positions = NULL;
    for (check_vnp = vnp->next, check_seq_num = seq_num + 1;
         check_vnp != NULL;
         check_vnp = check_vnp->next)
    {
      cfsp = (CreateFastaSeqPtr) check_vnp->data.ptrvalue;
      if (StringCmp (local_id, cfsp->local_id) == 0)
      {
        rval = FALSE;
        /* add to list of positions for report */
      	pos_vnp = ValNodeNew (appearance_positions);
      	if (appearance_positions == NULL)
      	{
      	  appearance_positions = pos_vnp;
          /* add to list of strings to ignore in further checking */
      	  ignore_vnp = ValNodeNew (ignore_list);
      	  if (ignore_list == NULL)
      	  {
      	    ignore_list = ignore_vnp;
      	  }
      	  if (ignore_vnp != NULL)
      	  {
      	    ignore_vnp->data.ptrvalue = local_id;
      	  }
      	}
      	if (pos_vnp != NULL)
      	{
      	  pos_vnp->data.intvalue = check_seq_num;
      	}
      }
    }
    if (appearance_positions != NULL)
    {
      Message (MSG_ERROR, "Sequence ID %s is not unique.  "
               "Please change your sequence IDs to be unique within this record.",
               local_id);
      ValNodeFree (appearance_positions);
    }
  	
  }
  
  ValNodeFree (ignore_list);
  return rval;
}

static CharPtr GetNthSeqId (ValNodePtr sequence_list, Int4 n)
{
  Int4              id_pos;
  ValNodePtr        vnp;
  CreateFastaSeqPtr cfsp;
  CharPtr           txt;
  CharPtr           fmt = "Unlabeled sequence %d";
  
  if (sequence_list == NULL) return NULL;
  for (id_pos = 0, vnp = sequence_list; id_pos < n - 1 && vnp != NULL; id_pos++, vnp = vnp->next)
  {	
  }
  if (vnp == NULL) return NULL;
  cfsp = (CreateFastaSeqPtr) vnp->data.ptrvalue;
  if (cfsp == NULL) return NULL;
  if (cfsp->local_id == NULL)
  {
  	txt = (CharPtr) MemNew (StringLen (fmt) + 25);
  	if (txt != NULL)
  	{
  	  sprintf (txt, fmt, n);
  	}
  }
  else
  {
   txt = StringSave (cfsp->local_id);
  }
  return txt;
}

static CharPtr ReportMultipleAppearances (ValNodePtr app_pos, ValNodePtr sequence_list)
{
  Int4       num_pos, pos_ctr;
  ValNodePtr vnp, seq_id_vnp;
  CharPtr    fmt = "Identical sequences appear for ";
  CharPtr    msg_txt;
  ValNodePtr seq_id_list = NULL;
  Int4       txt_len = StringLen (fmt) + 1;
  
  if (app_pos == NULL || sequence_list == NULL) return NULL;
  num_pos = ValNodeLen (app_pos);
  if (num_pos == 1) return NULL;
  for (vnp = app_pos; vnp != NULL; vnp = vnp->next)
  {
    msg_txt = GetNthSeqId (sequence_list, vnp->data.intvalue);
    if (msg_txt == NULL) continue;
  	seq_id_vnp = ValNodeNew (seq_id_list);
  	if (seq_id_list == NULL) seq_id_list = seq_id_vnp;
  	if (seq_id_vnp != NULL)
  	{
  	  seq_id_vnp->data.ptrvalue = msg_txt;
  	  txt_len += StringLen (msg_txt) + 5;
  	}
  }
  num_pos = ValNodeLen (seq_id_list);
  msg_txt = (CharPtr) MemNew (txt_len * sizeof (Char));
  if (msg_txt != NULL)
  {
    sprintf (msg_txt, fmt);
  	for (seq_id_vnp = seq_id_list, pos_ctr = 0;
  	     seq_id_vnp != NULL;
  	     seq_id_vnp = seq_id_vnp->next, pos_ctr ++)
  	{
  	  if (seq_id_vnp->next == NULL)
  	  {
  	  	StringCat (msg_txt, "and ");
  	  }
  	  StringCat (msg_txt, seq_id_vnp->data.ptrvalue);
  	  if (seq_id_vnp->next != NULL)
  	  {
  	    if (num_pos > 2)
  	    {
	      StringCat (msg_txt, ", ");
  	    }
  	    else
  	    {
  	  	  StringCat (msg_txt, " ");
  	    }
  	  }
  	}
  }
  ValNodeFreeData (seq_id_list);
  return msg_txt;
}


static Boolean CheckSequenceUniqueness (ValNodePtr list, Int4 skip)
{
  ValNodePtr        vnp, check_vnp, ignore_list, ignore_vnp;
  Int4              seq_num, check_seq_num;
  CreateFastaSeqPtr cfsp;
  Boolean           rval = TRUE;
  Int4              num_appearances;
  ValNodePtr        appearance_positions, pos_vnp;
  CharPtr           sequence;

  ignore_list = NULL;  
  for (vnp = list, seq_num = 1; vnp != NULL; vnp = vnp->next, seq_num++)
  {
    if (seq_num == skip) continue;
    cfsp = (CreateFastaSeqPtr) vnp->data.ptrvalue;
    sequence = cfsp->sequence;
    if (StringHasNoText (sequence)) continue;
    for (ignore_vnp = ignore_list; 
         ignore_vnp != NULL && StringCmp (sequence, ignore_vnp->data.ptrvalue) != 0;
         ignore_vnp = ignore_vnp->next)
    {
    }
    if (ignore_vnp != NULL) continue;
    num_appearances = 1;
    appearance_positions = NULL;
    for (check_vnp = vnp->next, check_seq_num = seq_num + 1;
         check_vnp != NULL;
         check_vnp = check_vnp->next, check_seq_num++)
    {
      cfsp = (CreateFastaSeqPtr) check_vnp->data.ptrvalue;
      if (StringCmp (sequence, cfsp->sequence) == 0)
      {
        rval = FALSE;
        /* add to list of positions for report */
      	pos_vnp = ValNodeNew (appearance_positions);
      	if (appearance_positions == NULL)
      	{
      	  appearance_positions = pos_vnp;
      	  /* add in first location */
      	  pos_vnp->data.intvalue = seq_num;
      	  pos_vnp = ValNodeNew (appearance_positions);
          /* add to list of strings to ignore in further checking */
      	  ignore_vnp = ValNodeNew (ignore_list);
      	  if (ignore_list == NULL)
      	  {
      	    ignore_list = ignore_vnp;
      	  }
      	  if (ignore_vnp != NULL)
      	  {
      	    ignore_vnp->data.ptrvalue = sequence;
      	  }
      	}
      	if (pos_vnp != NULL)
      	{
      	  pos_vnp->data.intvalue = check_seq_num;
      	}
      }
    }
    if (appearance_positions != NULL)
    {
      Message (MSG_ERROR, ReportMultipleAppearances (appearance_positions, list));
      ValNodeFree (appearance_positions);
    }
  	
  }
  
  ValNodeFree (ignore_list);
  return rval;
}

static Boolean IsValuePairInListPair (CharPtr qname, CharPtr qvalue, ListPairPtr lpp)
{
  ValNodePtr name_vnp, val_vnp;
  if (qname == NULL || qvalue == NULL || lpp == NULL) return FALSE;
  
  for (name_vnp = lpp->selected_names_list, val_vnp = lpp->selected_values_list;
       name_vnp != NULL && val_vnp != NULL;
       name_vnp = name_vnp->next, val_vnp = val_vnp->next)
  {
  	if (StringCmp (name_vnp->data.ptrvalue, qname) == 0)
  	{
  	  if (StringCmp (val_vnp->data.ptrvalue, qvalue) == 0)
  	  {
  	  	return TRUE;
  	  }
  	  else
  	  {
  	  	return FALSE;
  	  }
  	}
  }
  return FALSE;
}

static Boolean AreOrgsIdentical (CreateFastaSeqPtr seq1, CreateFastaSeqPtr seq2)
{
  ValNodePtr name_vnp, val_vnp;

  if (seq1 == NULL && seq2 == NULL) return TRUE;
  if (seq1 == NULL || seq2 == NULL) return FALSE;
  
  if (StringHasNoText (seq1->tax_name) && StringHasNoText (seq2->tax_name))
  {
  	return TRUE;
  }
  if (StringHasNoText (seq1->tax_name) || StringHasNoText (seq2->tax_name))
  {
  	return FALSE;
  }
  if (StringCmp (seq1->tax_name, seq2->tax_name) != 0)
  {
  	return FALSE;
  }
  /* check modifiers */
  for (name_vnp = seq1->mod_values.selected_names_list, val_vnp = seq1->mod_values.selected_values_list;
       name_vnp != NULL && val_vnp != NULL;
       name_vnp = name_vnp->next, val_vnp = val_vnp->next)
  {
    if (! StringHasNoText (name_vnp->data.ptrvalue)
        && !StringHasNoText (val_vnp->data.ptrvalue)
        && ! IsValuePairInListPair (name_vnp->data.ptrvalue, val_vnp->data.ptrvalue, &(seq2->mod_values)))
    {
      return FALSE;
    }
  }
  
  for (name_vnp = seq2->mod_values.selected_names_list, val_vnp = seq2->mod_values.selected_values_list;
       name_vnp != NULL && val_vnp != NULL;
       name_vnp = name_vnp->next, val_vnp = val_vnp->next)
  {
    if (! StringHasNoText (name_vnp->data.ptrvalue)
        && !StringHasNoText (val_vnp->data.ptrvalue)
        && ! IsValuePairInListPair (name_vnp->data.ptrvalue, val_vnp->data.ptrvalue, &(seq1->mod_values)))
    {
      return FALSE;
    }
  }
  
  return TRUE;
}

static Boolean ValidateOrganismList (ValNodePtr list, Int4 skip, Boolean is_popset)
{
  ValNodePtr        vnp, check_vnp, ignore_list, ignore_vnp;
  Int4              seq_num, check_seq_num;
  CreateFastaSeqPtr cfsp, check_cfsp;
  CharPtr           local_id; 
  Boolean           rval = TRUE;
  Int4              num_appearances;
  ValNodePtr        appearance_positions, pos_vnp;
  Char              str [255];
  Char              tmp [20];
  Boolean           tax_names_need_identical = FALSE;

  ignore_list = NULL;  
  for (vnp = list, seq_num = 1; vnp != NULL; vnp = vnp->next, seq_num++)
  {
    if (seq_num == skip) continue;
    cfsp = (CreateFastaSeqPtr) vnp->data.ptrvalue;
    local_id = cfsp->local_id;
    for (ignore_vnp = ignore_list; 
         ignore_vnp != NULL && ! AreOrgsIdentical (cfsp, ignore_vnp->data.ptrvalue);
         ignore_vnp = ignore_vnp->next)
    {
    }
    if (ignore_vnp != NULL) continue;
    num_appearances = 1;
    appearance_positions = NULL;
    for (check_vnp = vnp->next, check_seq_num = seq_num + 1;
         check_vnp != NULL;
         check_vnp = check_vnp->next, check_seq_num++)
    {
      check_cfsp = (CreateFastaSeqPtr) check_vnp->data.ptrvalue;
      if (is_popset)
      {
        if (StringCmp (cfsp->tax_name, check_cfsp->tax_name) != 0)
        {
          tax_names_need_identical = TRUE;
        }
      }
      if (AreOrgsIdentical (cfsp, check_cfsp))
      {
        rval = FALSE;
        num_appearances ++;
        /* add to list of positions for report */
      	pos_vnp = ValNodeNew (appearance_positions);
      	if (appearance_positions == NULL)
      	{
      	  appearance_positions = pos_vnp;
      	  if (pos_vnp != NULL)
      	  {
      	  	pos_vnp->data.intvalue = seq_num;
      	  }
      	  pos_vnp = ValNodeNew (appearance_positions);
      	  
          /* add to list of strings to ignore in further checking */
      	  ignore_vnp = ValNodeNew (ignore_list);
      	  if (ignore_list == NULL)
      	  {
      	    ignore_list = ignore_vnp;
      	  }
      	  if (ignore_vnp != NULL)
      	  {
      	    ignore_vnp->data.ptrvalue = cfsp;
      	  }
      	}
      	if (pos_vnp != NULL)
      	{
      	  pos_vnp->data.intvalue = check_seq_num;
      	}
      }
    }
    if (appearance_positions != NULL)
    {
      if (num_appearances == 2)
      {
      	sprintf (str, "sequences %d and %d", seq_num, appearance_positions->data.intvalue);
      }
      else if (num_appearances < 25)
      {
        str[0] = 0;
        for (pos_vnp = appearance_positions; pos_vnp->next != NULL; pos_vnp = pos_vnp->next)
        {
      	  sprintf (tmp, "%d, ", pos_vnp->data.intvalue);
          StringCat (str, tmp);
        }
        sprintf (tmp, " and %d", pos_vnp->data.intvalue);
        StringCat (str, tmp);
      }
      else 
      {
      	sprintf (str, "%d sequences", num_appearances);
      }
      
      if (is_popset)
      {
        Message (MSG_ERROR, "The qualifiers in %s are identical.  "
                 "Please change or add qualifiers to make your organisms unique within this record.",
                 str);
      }
      else
      {
        Message (MSG_ERROR, "The organisms in %s are identical.  "
                 "Please change your organisms or add qualifiers to make your organisms unique within this record.",
                 str);
       
      }
      ValNodeFree (appearance_positions);
    }
  }
  if (tax_names_need_identical && is_popset)
  {
    Message (MSG_ERROR, "This is a population study.  All of the organism names should be identical, but they are not.");
  }

  
  ValNodeFree (ignore_list);
  return rval;
}

static Boolean ValidateSequenceInformation (ValNodePtr list, Int4 skip)
{
  Boolean           rval = TRUE;
  ValNodePtr        vnp;
  CreateFastaSeqPtr cfsp;
  Int4              seq_num = 1;
  
  if (list == NULL) 
  {
    Message (MSG_ERROR, "You have no sequences.  Please add some.");
  	return FALSE;
  }
  for (vnp = list; vnp != NULL; vnp = vnp->next)
  {
  	cfsp = (CreateFastaSeqPtr) vnp->data.ptrvalue;
  	if (seq_num != skip)
  	{
      rval &= ValidateOneSequenceInformation (cfsp, seq_num);
  	}
    seq_num ++;
  }
  rval &= CheckLocalIDs (list, skip);
  return rval;
}

static void WriteSequenceInformationToFile (ValNodePtr list, FILE *fp)
{
  ValNodePtr        vnp, name_vnp, val_vnp;
  CreateFastaSeqPtr cfsp;
  CharPtr           cp;
  Int4              counter;
  
  for (vnp = list; vnp != NULL; vnp = vnp->next)
  {
  	cfsp = (CreateFastaSeqPtr) vnp->data.ptrvalue;
  	if (cfsp != NULL)
  	{
  	  if (cfsp->local_id == NULL)
  	  {
  	  	fprintf (fp, ">");
  	  }
  	  else
  	  {
        fprintf (fp, ">%s", cfsp->local_id);
  	  }
      if (!StringHasNoText (cfsp->tax_name))
      {
      	fprintf (fp, " [organism=%s]", cfsp->tax_name);
      }
      
      for (name_vnp = cfsp->mod_values.selected_names_list,
            val_vnp = cfsp->mod_values.selected_values_list;
           name_vnp != NULL && val_vnp != NULL;
           name_vnp = name_vnp->next, val_vnp = val_vnp->next)
      {
      	if (!StringHasNoText(name_vnp->data.ptrvalue)
      	    && ! StringHasNoText (val_vnp->data.ptrvalue))
      	{
      	  fprintf (fp, "[%s=%s]", name_vnp->data.ptrvalue, val_vnp->data.ptrvalue);
      	}
      }
      if (cfsp->title != NULL)
      {
      	fprintf (fp, " %s", cfsp->title);
      }
      fprintf (fp, "\n");
      counter = 0;
      if (cfsp->sequence != NULL)
      {
        for (cp = cfsp->sequence; *cp != 0; cp++)
        {
       	  if (isalpha (*cp))
      	  {
      	    fprintf (fp, "%c", *cp);
      	    counter++;
      	    if (counter == 80)
      	    {
      	  	  fprintf (fp, "\n");
      	  	  counter = 0;
      	    }
      	  }
        }    	
      }
      if (counter != 0)
      {
      	fprintf (fp, "\n");
      }
  	}    
  }
}

/* This function checks the values of the inputs for the seq_num sequence to determine
 * whether the user has made any changes from the defaults - i.e., did he enter
 * any qualifier values?  any sequence data?  a different organism name from the
 * default?  a sequence ID?
 * If so, the sequence is not "empty".
 */
static Boolean IsThisSequenceEmpty (CreateFastaPtr cfp, Int4 seq_num)
{
  Int4              val;
  ValNodePtr        vnp;
  CreateFastaSeqPtr cfsp;
  
  if (cfp == NULL) return TRUE;
  
  for (val = 1, vnp = cfp->sequence_list; vnp != NULL && val < seq_num; vnp = vnp->next, val++)
  {
  }
  if (val != seq_num || vnp == NULL)
  {
  	return TRUE;
  }
  cfsp = (CreateFastaSeqPtr) vnp->data.ptrvalue;
  if (cfsp == NULL) return TRUE;
  CollectOneSequenceInformation (cfsp);

  if (StringHasNoText (cfp->default_tax_name))
  {
  	if (!StringHasNoText (cfsp->tax_name))
  	{
  	  return FALSE;
  	}
  }
  else if (StringCmp (cfp->default_tax_name, cfsp->tax_name) != 0 && !StringHasNoText (cfsp->tax_name))
  {
  	return FALSE;
  }
  
  if (!StringHasNoText (cfsp->sequence))
  {
  	return FALSE;
  }
  
  if (! cfp->auto_id && ! StringHasNoText (cfsp->local_id))
  {
  	return FALSE;
  }

  for (vnp = cfsp->mod_values.selected_values_list; vnp != NULL; vnp = vnp->next)
  {
  	if (!StringHasNoText (vnp->data.ptrvalue))
  	{
  	  return FALSE;
  	}
  }
  	
  return TRUE;
}


static void DoCreateFASTAFile (CreateFastaPtr cfp)
{
  Char           path [PATH_MAX];
  FILE           *fp;
  Int4           to_remove = -1;
  Boolean        remove_current = FALSE;
  Int4           num_seq;
  ValNodePtr     vnp, prev;
  Int4           val;

  if (cfp == NULL) return;
  
  /* if we are creating a multiple sequence set, the user may have hit
   * Add Sequence to "finish" his last sequence when he didn't need to.
   * If the sequence dialog is "empty" and we have enough sequences for
   * a set, we can remove the empty sequence.
   */
  if (!cfp->fpp->single)
  {
    to_remove = GetValue (cfp->sequence_selector);
    num_seq = ValNodeLen (cfp->sequence_list);
    if (num_seq > 2)
    {
      remove_current = IsThisSequenceEmpty (cfp, to_remove);	
    }
    if (! remove_current)
    {
      to_remove = -1;
    }
  }
  
  if (! ValidateSequenceInformation (cfp->sequence_list, to_remove)) 
  {
    return; 
  }
  
  if (cfp->sequence_list->next == NULL && ! cfp->fpp->single)
  {
  	Message (MSG_ERROR, "You must have at least two sequences for a set.");
  	return;
  }
  if (!ValidateOrganismList (cfp->sequence_list, to_remove, cfp->is_popset))
  {
    if (Message (MSG_YN, "Do you wish to edit your organisms?") == ANS_YES) return;
  }
  if (!CheckQualListForConsistency (cfp)) 
  {
    if (Message (MSG_YN, "Do you wish to edit your modifiers?") == ANS_YES) return;
  }
  if (! CheckSequenceUniqueness (cfp->sequence_list, to_remove))
  {
  	if (Message (MSG_YN, "Do you wish to edit your sequences?") == ANS_YES) return;
  }
  
  if (cfp->fpp->single)
  {
    if (Message (MSG_YN, "This will complete the creation of your FASTA file.  "
                 "The new file will automatically be imported into Sequin.  "
                 "You will not be able to change the sequence in this file if you continue."
                 "Are you sure that you are done editing this sequence?  "
                 "(Click Yes to continue and save your file)") == ANS_NO) return;
  }
  else
  {
    if (Message (MSG_YN, "This will complete the creation of your FASTA file.  "
                 "The new file will automatically be imported into Sequin.  "
                 "You have %d sequences.  "
                 "You will not be able to add more sequences to this file if you continue.  "
                 "Are you sure that you are done adding sequences?  "
                 "(Click Yes to continue and save your file)",
                 ValNodeLen (cfp->sequence_list)) == ANS_NO) return;
  }

  if (! GetOutputFileName (path, sizeof (path), NULL)) 
  {
    return;
  }
  fp = FileOpen (path, "w");
  if (fp == NULL) 
  {
    Message (MSG_ERROR, "Unable to write to %s", path);
    return;
  }
  
  /* remove current sequence if empty */
  if (remove_current)
  {
    prev = NULL;
    for (val = 1, vnp = cfp->sequence_list; vnp != NULL && val < to_remove; vnp = vnp->next, val++)
    {
  	  prev = vnp;
    }
    if (vnp == NULL) return;
    if (prev == NULL)
    {
  	  cfp->sequence_list = vnp->next;
    }
    else 
    {
  	  prev->next = vnp->next;
    }
    vnp->next = NULL;
    FreeSequenceGroup (vnp);  	
  }  
  
  WriteSequenceInformationToFile (cfp->sequence_list, fp);
  FileClose (fp);
  SetObjectExtra (cfp->form, cfp->fpp, NULL);
  ImportFastaDialog ((DialoG) cfp->form, path);
  SetObjectExtra (cfp->form, cfp, NULL);
  SetStatus (cfp->use_id_from_fasta_defline, TRUE);
  cfp->fpp->parseSeqId = TRUE;
  Remove (cfp->form);
}

static void DoCreateFASTAFileButton (ButtoN b)
{
  CreateFastaPtr cfp;
  
  cfp = (CreateFastaPtr) GetObjectExtra (b);
  DoCreateFASTAFile (cfp);
}

static void DoCreateFASTAFileItem (IteM i)
{
  CreateFastaPtr cfp;
  
  cfp = (CreateFastaPtr) GetObjectExtra (i);
  DoCreateFASTAFile (cfp);
}

static void PreviewFastaFile (CreateFastaPtr cfp)
{
  Char         path [PATH_MAX];
  FILE         *fp;
  
  if (cfp == NULL) return;

  TmpNam (path);
  fp = FileOpen (path, "wb");
  if (fp == NULL) 
  {
    Message (MSG_ERROR, "Unable to open temporary file");
    return;
  }

  ValidateSequenceInformation (cfp->sequence_list, -1);
  WriteSequenceInformationToFile (cfp->sequence_list, fp);

  FileClose (fp);
  LaunchGeneralTextViewer (path, "FASTA File");
  FileRemove (path);
}

static void PreviewFastaFileButton (ButtoN b)
{
  CreateFastaPtr cfp;
  
  cfp = (CreateFastaPtr) GetObjectExtra (b);
  PreviewFastaFile (cfp);  
}

static void PreviewFastaFileItem (IteM i)
{
  CreateFastaPtr cfp;
  
  cfp = (CreateFastaPtr) GetObjectExtra (i);
  PreviewFastaFile (cfp);  
}

static void CleanupCreateFastaForm (
  GraphiC g,
  VoidPtr data
)

{
  CreateFastaPtr  cfp;

  cfp = (CreateFastaPtr) data;
  FreeCreateFasta (cfp);
}

static void CreateFastaCancel (CreateFastaPtr cfp)
{
  if (cfp == NULL) return;
  
  if (cfp->fpp->single)
  {
    if (Message (MSG_YN, "Your sequence will not be saved if you cancel.  Do you wish to cancel?") == ANS_NO)
    {
      return;
    }
  }
  else
  {
    if (Message (MSG_YN, "None of your sequences will be saved if you cancel.  Do you wish to cancel?") == ANS_NO)
    {
      return;
    }
  }
  Remove (cfp->form);
}

static void CreateFastaCancelButton (ButtoN b)
{
  CreateFastaPtr  cfp;
  cfp = GetObjectExtra (b);
  CreateFastaCancel (cfp);
}

static void CreateFastaCancelItem (IteM i)
{
  CreateFastaPtr  cfp;
  cfp = GetObjectExtra (i);
  CreateFastaCancel (cfp);
}

static void AddFastaSequenceButton (ButtoN b)
{
  CreateFastaPtr   cfp;
  cfp = (CreateFastaPtr) GetObjectExtra (b);
  AddFastaSequence (cfp);
}

static void AddFastaSequenceItem (IteM i)
{
  CreateFastaPtr   cfp;
  cfp = (CreateFastaPtr) GetObjectExtra (i);
  AddFastaSequence (cfp);
}

static void RemoveFastaSequenceButton (ButtoN b)
{
  CreateFastaPtr   cfp;
  cfp = (CreateFastaPtr) GetObjectExtra (b);
  RemoveFastaSequence (cfp);
}

static void RemoveFastaSequenceItem (IteM i)
{
  CreateFastaPtr   cfp;
  cfp = (CreateFastaPtr) GetObjectExtra (i);
  RemoveFastaSequence (cfp);
}

static void SequenceIdHelp (IteM i)
{
  Message (MSG_OK,
        "  A sequence ID is a temporary ID which will be replaced with a unique"
        " GenBank accession number by the GenBank curators.  The sequence ID should"
        " be unique for each sequence in a record.  It could represent a clone,"
        " isolate, a laboratory designation for your specimen, or some other useful"
        " information, but this is not required."
        "  A sequence ID may not begin with a '>' character or contain spaces.");
}

static void OrganismNameHelp (IteM i)
{
  Message (MSG_OK, "You must enter a scientific name for your organism. "
           "The name does not need to be present in the list in the dialog.");
}

static void OrganismQualifiersHelp (IteM i)
{
  Message (MSG_OK, "If you have multiple organisms with the same scientific name, "
           "please use modifiers to distinguish the organisms from one another.  "
           "Strain, clone, isolate, and specimen voucher are modifiers frequently "
           "used for this purpose, but you may select any applicable modifiers.");
}

static void SequenceCharactersHelp (IteM i)
{
  Message (MSG_OK, "Please enter the nucleotides for your sequence into the  "
                   "sequence characters area.  You may only use the valid "
                   "IUPAC characters (%s).  You may not use *, -, ., or any other "
                   "alignment characters or punctuation in your sequence.  If you "
                   "are trying to import an alignment, you should use the Cancel "
                   "button to exit this dialog, then hit the Prev Page button and then "
                   "the Prev Form button to get to the Sequence Format dialog, and "
                   "select 'Alignment' instead of FASTA for the sequence data format.  "
                   "You will need to have a file prepared for import.", 
                    valid_iupac_characters);
}

const CharPtr bracket_mismatch_msg = "You have mismatched brackets at line %d.";
const CharPtr missing_equals_msg = "Your bracketed pairs must be in the form [qualifier=value].  "
  	                    "You are missing an equals sign in line %d.";
  	                    
/* later, add handling for equals signs without brackets? */
static CharPtr SuggestCorrectBracketing (CharPtr str)
{
  CharPtr    start, stop, next_start, next_stop, eq_loc;
  ValNodePtr pieces = NULL;
  CharPtr    cp = str;
  Boolean    done = FALSE;
  ValNodePtr vnp;
  Int4       len;
  CharPtr    txt;

  if (str == NULL) return NULL;
  
  while (!done && *cp != 0)
  {  	
    start = StringChr (cp, '[');
    stop = StringChr (cp, ']');
    eq_loc = StringChr (cp, '=');
    if (start == NULL)
    {
      next_start = NULL;
    }
    else
    {
      next_start = StringChr (start + 1, '[');    	
    }

    if (start == NULL && stop == NULL) 
    {
      len = StringLen (cp) + 1;
      if (len > 1)
      {
        txt = (CharPtr) MemNew (len * sizeof (Char));
        if (txt != NULL)
        {
      	  StringNCpy (txt, cp, len - 1);
      	  txt [len - 1] = 0;
          vnp = ValNodeAdd (&pieces);
  		  if (vnp != NULL)
  		  {
            vnp->data.ptrvalue = txt;
  		  }
        }
      }
      done = TRUE;
    }
    else if (start != NULL && stop != NULL && eq_loc != NULL
             && start < eq_loc && eq_loc < stop 
             && (next_start == NULL || next_start > stop))
    {
      len = stop - cp + 2;
      txt = (CharPtr) MemNew (len * sizeof (Char));
      if (txt != NULL)
      {
      	StringNCpy (txt, cp, len - 1);
      	txt [len - 1] = 0;
        vnp = ValNodeAdd (&pieces);
  		if (vnp != NULL)
  		{
          vnp->data.ptrvalue = txt;
  		}
      }
      cp = stop + 1;
    }
    else if (start == NULL || (stop != NULL && start > stop))
    {
  	  eq_loc = StringChr (cp, '=');
  	  if (eq_loc == NULL || eq_loc == cp || eq_loc > stop)
  	  {
  	    /* if there is no equals sign, remove the offending bracket */
  		len = stop - cp + 1;
  		txt = (CharPtr) MemNew (len * sizeof (Char));
  		if (txt != NULL)
  		{
  		  StringNCpy (txt, cp, len - 1);
  		  txt [len - 1] = 0;
  		  vnp = ValNodeAdd (&pieces);
  		  if (vnp != NULL)
  		  {
  		  	vnp->data.ptrvalue = txt;
  		  }
  		}
  		cp = stop + 1;
  	  }
  	  else 
  	  {
  	  	/* find the first non-alphabet character before the equals sign and put in a bracket */
  	  	start = eq_loc - 1;
  	  	/* skip over whitespace before equals sign */
  	  	while (cp != start && isspace (*start))
  	  	{
  	  	  start --;
  	  	}
  	  	/* back up past qualifier name */
  	  	while (cp != start && isalpha (*start))
  	  	{
  	  	  start --;
  	  	}
  	  	/* now insert left bracket */
  	  	len = stop - cp + 1;
  	  	txt = (CharPtr) MemNew (len * sizeof (Char));
  		if (txt != NULL)
  		{
  		  if (start > cp)
  		  {
  		    StringNCpy (txt, cp, start - cp - 1);
  		  }
  		  StringCat (txt, "[");
  		  StringNCat (txt, start, stop - start);
  		  txt [len - 1] = 0;
  		  vnp = ValNodeAdd (&pieces);
  		  if (vnp != NULL)
  		  {
  		  	vnp->data.ptrvalue = txt;
  		  }
  		}
  		cp = stop + 1;
  	  }
    }
    else if (stop != NULL && eq_loc != NULL && eq_loc > stop)
    {
      next_stop = StringChr (stop + 1, ']');
      if (next_stop != NULL && next_stop < next_start && eq_loc < next_stop)
      {
      	/* remove the intermediate stop */
      	len = next_stop - cp;
      	txt = (CharPtr) MemNew (len * sizeof (Char));
  		if (txt != NULL)
  		{
  		  StringNCpy (txt, cp, stop - cp);
  		  StringNCat (txt, stop + 1, next_stop - stop);
  		  vnp = ValNodeAdd (&pieces);
  		  if (vnp != NULL)
  		  {
  		  	vnp->data.ptrvalue = txt;
  		  }
  		}
  		cp = next_stop + 1;
      }
      else
      {
      	/* remove both the start and stop */
      	len = stop - cp - 1;
      	txt = (CharPtr) MemNew (len * sizeof (Char));
  		if (txt != NULL)
  		{
          if (start > cp)
          {
          	StringNCpy (txt, cp, start - cp);
          }
          StringNCat (txt, start + 1, stop - start - 1);
  		  vnp = ValNodeAdd (&pieces);
  		  if (vnp != NULL)
  		  {
  		  	vnp->data.ptrvalue = txt;
  		  }
  		}
  		cp = stop + 1;
      }
    }
    else
    {
      /* we have a start without a stop */
      eq_loc = StringChr (start, '=');
      next_start = StringChr (start + 1, '[');
      if (eq_loc == NULL || (next_start != NULL && eq_loc > next_start))
      {
      	/* if we have no equals sign, remove the offending bracket */
      	if (next_start == NULL)
      	{
      	  /* if there are no more starts, copy the rest of the string and finish */
      	  len = StringLen (cp);
      	  done = TRUE;
      	}
      	else
      	{
      	  /* copy up to the next start */
      	  len = next_start - cp;      		
      	}
      	if (len > 1)
      	{
      	  txt = (CharPtr) MemNew (len * sizeof (Char));
  		  if (txt != NULL)
  		  {
  		    if (cp < start)
  		    {
 		      StringNCpy (txt, cp, start - cp);
 		      if (next_start - start > 1)
 		      {
 		      	StringNCat (txt, start + 1, next_start - start - 1);
 		      }
  		    }
  		    else
  		    {
  		      StringNCpy (txt, start + 1, len);  		    	
  		    }
  		    vnp = ValNodeAdd (&pieces);
  		    if (vnp != NULL)
  		    {
  		      vnp->data.ptrvalue = txt;
  		    }
      	  }
      	}
      	cp += len;
      }
      else
      {
      	/* put everything before the next start inside the bracket */
      	if (next_start == NULL)
      	{
      	  len = StringLen (cp) + 2;
      	  done = TRUE;
      	}
      	else
      	{
      	  len = next_start - cp + 2;
      	}
      	txt = (CharPtr) MemNew (len * sizeof (Char));
  		if (txt != NULL)
  		{
  		  StringNCpy (txt, cp, len - 2);
  		  StringCat (txt, "]");
  		  vnp = ValNodeAdd (&pieces);
  		  if (vnp != NULL)
  		  {
  		    vnp->data.ptrvalue = txt;
  		  }
      	}
      	cp += len - 2;
      }
    }
  }
  
  txt = MergeValNodeStrings (pieces, FALSE);
  ValNodeFreeData (pieces);
  return txt;
}

static CharPtr DetectBadBracketing (CharPtr str)
{
  CharPtr start, stop, next_start, next_stop, eq_loc;
  
  
  if (str == NULL) return NULL;
  
  start = StringChr (str, '[');
  stop = StringChr (str, ']');
  if (start == NULL && stop == NULL) return NULL;
  if ((start != NULL && stop == NULL)
   || (start == NULL && stop != NULL)
   || (start > stop))
  {
  	return bracket_mismatch_msg;
  }
  eq_loc = StringChr (start + 1, '=');
  if (eq_loc == NULL || eq_loc > stop)
  {
    return missing_equals_msg;
  }
  while (start != NULL && stop != NULL)
  {
    next_start = StringChr (start + 1, '[');
    next_stop = StringChr (stop + 1, ']');
    if (next_start == NULL && next_stop == NULL) return FALSE;

    if ((next_start != NULL && next_stop == NULL)
     || (next_start == NULL && next_stop != NULL)
     || (next_start > next_stop)
     || (next_start < stop))
    {
  	  return bracket_mismatch_msg;
    }
    eq_loc = StringChr (next_start + 1, '=');
    if (eq_loc == NULL || eq_loc > next_stop)
    {
      return missing_equals_msg;
    }
    
    start = next_start;
    stop = next_stop;    
  }
  return NULL;
}

typedef struct fixbadlineform 
{
  WindoW  w;
  PrompT  prompt;
  TexT    new_line;	
  CharPtr new_text;
  Int4    line_num;
  Boolean done;
  Boolean cancelled;
  CharPtr orig_txt;
} FixBadLineFormData, PNTR FixBadLineFormPtr;

static void SetBadLineFormPrompt (FixBadLineFormPtr fp, CharPtr msg)
{
  CharPtr str;
  if (fp == NULL || msg == NULL) return;
  
  str = (CharPtr) MemNew (sizeof (Char) * (StringLen (msg) + 15));
  if (str != NULL)
  {
    sprintf (str, msg, fp->line_num);
  }
  SetTitle (fp->prompt, str);
}

static void FixBadLineOk (ButtoN b)
{
  FixBadLineFormPtr fp;
  CharPtr           msg;
  
  fp = (FixBadLineFormPtr) GetObjectExtra (b);
  if (fp == NULL) return;
  fp->new_text = MemFree (fp->new_text);
  fp->new_text = SaveStringFromText (fp->new_line);
  msg = DetectBadBracketing (fp->new_text);
  if (msg != NULL)
  {
  	SetBadLineFormPrompt (fp, msg);
  	return;
  }
  Remove (fp->w);
  fp->cancelled = FALSE;
  fp->done = TRUE;
}

static void FixBadLineCancel (ButtoN b)
{
  FixBadLineFormPtr fp;
  
  fp = (FixBadLineFormPtr) GetObjectExtra (b);
  if (fp == NULL) return;
  fp->new_text = MemFree (fp->new_text);
  fp->cancelled = TRUE;
  
  Remove (fp->w);
  fp->done = TRUE;
}

static void SuggestBracketFix (ButtoN b)
{
  FixBadLineFormPtr fp;
  
  fp = (FixBadLineFormPtr) GetObjectExtra (b);
  if (fp == NULL) return;
  fp->new_text = MemFree (fp->new_text);
  fp->new_text = SaveStringFromText (fp->new_line);
  SetTitle (fp->new_line, SuggestCorrectBracketing(fp->new_text));
  return;
}

static void ResetBracketFixText (ButtoN b)
{
  FixBadLineFormPtr fp;
  
  fp = (FixBadLineFormPtr) GetObjectExtra (b);
  if (fp == NULL) return;
  SetTitle (fp->new_line, fp->orig_txt);
  return;
}


static CharPtr FixBadBracketing (CharPtr bad_line, CharPtr msg, Int4 line_num, BoolPtr cancelled)
{
  GrouP  g, c;
  ButtoN b;
  FixBadLineFormData fd;
  Int4               len;
  
  fd.w = ModalWindow(-20, -13, -10, -10, NULL);
  g = HiddenGroup(fd.w, 0, 4, NULL);
  len = StringLen (bad_line) + 5;
  fd.orig_txt = bad_line;
  fd.line_num = line_num;
  fd.new_text = NULL;
  fd.prompt = StaticPrompt (g, "", 0, popupMenuHeight, programFont, 'l');
  SetBadLineFormPrompt (&fd, msg);
  fd.new_line = DialogText (g, bad_line, len, NULL);
  c = HiddenGroup (g, 4, 0, NULL);
  b = PushButton(c, "OK", FixBadLineOk);
  SetObjectExtra (b, &fd, NULL);
  b = PushButton(c, "Cancel", FixBadLineCancel);
  SetObjectExtra (b, &fd, NULL);
  b = PushButton (c, "Suggest Correction", SuggestBracketFix);
  SetObjectExtra (b, &fd, NULL);
  b = PushButton (c, "Reset to original text", ResetBracketFixText);
  SetObjectExtra (b, &fd, NULL);
  
  Show(fd.w); 
  Select (fd.w);
  fd.done = FALSE;
  while (!fd.done)
  {
    ProcessExternalEvent ();
    Update ();
  }
  ProcessAnEvent ();
  if (fd.cancelled)
  {
  	if (cancelled != NULL)
  	{
      *cancelled = TRUE;
  	}
  	return NULL;
  }
  return fd.new_text;
}

static CharPtr GetUnbracketedText (CharPtr str)
{
  CharPtr unbr_text, src, dst;
  
  if (StringHasNoText (str))
  {
  	return NULL;
  }
  unbr_text = StringSave (str);
  dst = StringChr (unbr_text, '[');
  src = StringChr (unbr_text, ']');
  while (dst != NULL && src != NULL && *src != 0)
  {
    src++;
    while (*src != '[' && *src != 0)
    {
      *dst = *src;
      dst++;
      src++;
    }
  	src = StringChr (src, ']');
  }
  if (dst != NULL)
  {
    *dst = 0;    
  }
  return unbr_text;
}

typedef struct stringpair
{
  CharPtr findstr;
  CharPtr replstr;
  Int4    replint;
  Boolean is_org;
} StringPairData, PNTR StringPairPtr;

static ValNodePtr FreeStringPairList (ValNodePtr list)
{
  StringPairPtr spp;
  
  if (list == NULL) return NULL;
  list->next = FreeStringPairList (list->next);
  spp = (StringPairPtr) list->data.ptrvalue;
  if (spp != NULL)
  {
  	MemFree (spp->findstr);
  	MemFree (spp->replstr);
  	MemFree (spp);
  	list->data.ptrvalue = NULL;
  }
  list = ValNodeFree (list);
  return list;
}

typedef struct fixmodifierform
{
  WindoW     w;
  PopuP      qualtype_selector;
  ValNodePtr qual_list;
  ButtoN     add_to_fixes;
  Boolean    done;
  Boolean    cancelled;
} FixModifierFormData, PNTR FixModifierFormPtr;

static void FixQualTypeOk (ButtoN b)
{
  FixModifierFormPtr fp;
  Int4               val;
  
  fp = (FixModifierFormPtr) GetObjectExtra (b);
  if (fp == NULL) return;
  val = GetValue (fp->qualtype_selector);
  if (val == 0) return;
  
  Remove (fp->w);
  fp->cancelled = FALSE;
  fp->done = TRUE;	
}

static void FixQualTypeCancel (ButtoN b)
{
  FixModifierFormPtr fp;
  
  fp = (FixModifierFormPtr) GetObjectExtra (b);
  if (fp == NULL) return;
  
  Remove (fp->w);
  fp->cancelled = TRUE;
  fp->done = TRUE;	
}

static StringPairPtr 
GetModifierFix (ModifierInfoPtr mip, Int4 line_num, BoolPtr cancelled, BoolPtr add_to_fixes)
{
  GrouP  g, c;
  ButtoN b;
  FixModifierFormData fd;
  CharPtr            str = NULL;
  CharPtr            prompt_fmt = "In line %d, %s is not a valid qualifier type.  "
                                  "Please choose a valid qualifier type:";
  ValNodePtr         vnp_new, vnp;
  StringPairPtr      spp = NULL;
  Int4               val;
  Int4               pos;
  CharPtr            old_name;
  
  fd.w = ModalWindow(-20, -13, -10, -10, NULL);
  g = HiddenGroup(fd.w, 0, 4, NULL);
  str = MemNew ((StringLen (mip->name) + StringLen (prompt_fmt) + 15) * sizeof (Char));
  sprintf (str, prompt_fmt, line_num, mip->name);
  MultiLinePrompt (g, str, 27 * stdCharWidth, programFont);
  fd.qual_list = GetQualList ();
  vnp_new = ValNodeNew (NULL);
  if (vnp_new != NULL)
  {
  	vnp_new->data.ptrvalue = "organism";
  	vnp_new->choice = 0;
  	vnp_new->next = fd.qual_list;
  	fd.qual_list = vnp_new;
  }
  fd.qualtype_selector = PopupList (g, TRUE, NULL); 
  InitValNodePopup (fd.qual_list, fd.qualtype_selector);
  fd.add_to_fixes = CheckBox (g, "Replace all instances", NULL);
  c = HiddenGroup (g, 4, 0, NULL);
  b = PushButton(c, "OK", FixQualTypeOk);
  SetObjectExtra (b, &fd, NULL);
  b = PushButton(c, "Cancel", FixQualTypeCancel);
  SetObjectExtra (b, &fd, NULL);
  
  Show(fd.w); 
  Select (fd.w);
  fd.done = FALSE;
  while (!fd.done)
  {
    ProcessExternalEvent ();
    Update ();
  }
  ProcessAnEvent ();
  
  if (fd.cancelled)
  {
  	if (cancelled != NULL)
  	{
  	  *cancelled = TRUE;
  	}
  }
  else
  {
  	val = GetValue (fd.qualtype_selector);
  	for (vnp = fd.qual_list, pos = 1; vnp != NULL && pos < val; vnp = vnp->next, pos++)
  	{
  	}
  	if (vnp != NULL)
  	{
  	  old_name = mip->name;
  	  mip->name = StringSave (vnp->data.ptrvalue);
  	  if (StringCmp (vnp->data.ptrvalue, "organism") == 0)
  	  {
  	  	mip->is_org = TRUE;
  	  }
      else
      {
      	mip->subtype = vnp->choice;
      }  	  
      if (add_to_fixes != NULL)
      {
        *add_to_fixes = GetStatus (fd.add_to_fixes);
      }
      if (GetStatus (fd.add_to_fixes))
      {
        spp = (StringPairPtr) MemNew (sizeof(StringPairData));
  	    if (spp != NULL)
  	    {
    	  spp->findstr = old_name;
    	  spp->replstr = StringSave (vnp->data.ptrvalue);
    	  if (mip->is_org)
    	  {
    	  	spp->is_org = TRUE;
    	  }
    	  else
    	  {
    	  	spp->replint = mip->subtype;
    	  }
  	    }
  	    else 
  	    {
  	      MemFree (old_name);
  	    }
      }
      else
      {
      	MemFree (old_name);
      }
  	}
  }
  ValNodeFree (fd.qual_list);
  return spp;
}

static ValNodePtr FixUnrecognizedModifier 
(ModifierInfoPtr mip, ValNodePtr fixes, Int4 line_num, BoolPtr cancelled)
{
  ValNodePtr    vnp;
  StringPairPtr spp;
  Boolean       found = FALSE;
  Boolean       add_to_fixes = FALSE;
  
  if (mip == NULL) return fixes;
  
  /* try to find mip->name in list of automatic fixes */
  for (vnp = fixes; vnp != NULL && !found; vnp = vnp->next)
  {
  	spp = (StringPairPtr) vnp->data.ptrvalue;
  	if (spp != NULL && StringCmp (spp->findstr, mip->name) == 0)
  	{
  	  found = TRUE;
  	  mip->name = MemFree (mip->name);
  	  mip->name = StringSave (spp->replstr);
  	  if (spp->is_org)
  	  {
  	  	mip->is_org = TRUE;
  	  }
  	  else
  	  {
  	  	mip->subtype = spp->replint;
  	  }
  	}
  }
  if (!found)
  {
  	/* get new fix and add to list*/
    spp = GetModifierFix (mip, line_num, cancelled, &add_to_fixes);
    if (add_to_fixes && spp != NULL)
    {
      vnp = ValNodeNew (fixes);
      if (fixes == NULL)
      {
      	fixes = vnp;
      }
      if (vnp != NULL)
      {
      	vnp->data.ptrvalue = spp;
      }
    }
  }
  return fixes;
}

/* remove bracketed pairs, return remaining title characters */
static CharPtr ParseImportFastaOrgAndQuals 
(CreateFastaSeqPtr cfsp, CharPtr str, Int4 line_num, BoolPtr cancelled, ValNodePtr PNTR fixes)
{
  CharPtr    cp;
  CharPtr    msg;
  ValNodePtr mod_list;
  CharPtr    title;
  ValNodePtr vnp, name_vnp, val_vnp;
  ModifierInfoPtr mip;
  
  if (cfsp == NULL || str == NULL) return NULL;
 
  cp = StringSave (str);
  msg = DetectBadBracketing (cp);
  if (msg != NULL)
  {
  	str = FixBadBracketing (cp, msg, line_num, cancelled);
  	MemFree (cp);
  	if (*cancelled)
  	{
  	  return NULL;
  	}
  	cp = str;
  }

  mod_list = ParseAllBracketedModifiers (cp);
  title = GetUnbracketedText (cp);
  
  for (vnp = mod_list; vnp != NULL && !*cancelled; vnp = vnp->next)
  {
    mip = (ModifierInfoPtr)vnp->data.ptrvalue;
    if (mip == NULL) continue;
    if (!mip->is_org && mip->subtype == 255)
    {
      /* fix qualifier name */
      *fixes =  FixUnrecognizedModifier (mip, *fixes, line_num, cancelled);
    }
    if (! *cancelled)
    {
      if (mip->is_org)
      {
        cfsp->tax_name = MemFree (cfsp->tax_name);
        cfsp->tax_name = StringSave (mip->value);
      }
      else
      {
        name_vnp = ValNodeNew (cfsp->mod_values.selected_names_list);
        if (cfsp->mod_values.selected_names_list == NULL)
        {
          cfsp->mod_values.selected_names_list = name_vnp;
        }
        if (name_vnp != NULL)
        {
          name_vnp->choice = mip->subtype;
          name_vnp->data.ptrvalue = StringSave (mip->name);
        }
        val_vnp = ValNodeNew (cfsp->mod_values.selected_values_list); 
        if (cfsp->mod_values.selected_values_list == NULL)	
        {
          cfsp->mod_values.selected_values_list = val_vnp;
        }
        if (val_vnp != NULL)
        {
          val_vnp->data.ptrvalue = StringSave (mip->value);
        }      	
      }	
    }
  }  

  mod_list = ModifierInfoListFree (mod_list);
  MemFree (cp);
  
  return title;
}

static CreateFastaSeqPtr ParseImportFastaDefLine 
(CharPtr         line,
 Int4            line_num,
 CreateFastaPtr  cfp,
 BoolPtr         cancelled, 
 ValNodePtr PNTR qualfixes,
 Int4            to_remove)
{
  ValNodePtr        vnp;
  CreateFastaSeqPtr cfsp;
  CharPtr           cp, cpend, cp_return;
  
  if (line == NULL) return NULL;
  
  cfp->sequence_list = AddSequenceGroup (cfp->sequence_list, 
                                         cfp->selected_qual_list,
                                         cfp->default_tax_name,
                                         cfp->auto_id,
                                         to_remove,
                                         cfp);
  vnp = cfp->sequence_list;
  while (vnp != NULL && vnp->next != NULL)
  {
    vnp = vnp->next;
  }
  if (vnp != NULL)
  {
    cfsp = (CreateFastaSeqPtr) vnp->data.ptrvalue;
  }
  if (cfsp == NULL) return NULL;

  cp = line;
  if (*cp == '>')
  {
  	cp++;
  }
  while (isspace (*cp))
  {
  	cp++;
  }
  
  /* find bracketed pairs and parse to quals or organism name, remove from line */
  cp_return = ParseImportFastaOrgAndQuals (cfsp, cp, line_num, cancelled, qualfixes);
  if (*cancelled)
  {
  	return NULL;
  }
  
  
  /* take remaining characters and use as sequence ID or title */
  if (cfp->auto_id)
  {
  	cfsp->title = StringSave (cp_return);
  }
  else
  {
  	cpend = cp_return + 1;
  	while (*cpend != 0 && !isspace (*cpend))
  	{
  	  cpend ++;
  	}
  	if (*cpend == 0)
  	{
  	  cfsp->local_id = StringSave (cp_return);
  	}
  	else
  	{
  	  *cpend = 0;
  	  cfsp->local_id = StringSave (cp_return);
  	  cfsp->title = StringSave (cpend + 1);
  	}
  }
  if (cp_return != cp)
  {
  	MemFree (cp_return);
  }
  return cfsp;	
}

static CreateFastaSeqPtr ParseImportFastaSeqLine 
(CreateFastaSeqPtr cfsp,
 CharPtr           line,     
 Int4              line_num,
 CreateFastaPtr    cfp,
 BoolPtr           cancelled,
 ValNodePtr PNTR   qualfixes,
 Int4              to_remove)
{
  CharPtr    cp;
  Char       ch;
  MsgAnswer  ans;
  ValNodePtr vnp;
  Int4       new_len = 0;
  CharPtr    new_seq;
  
  if (line == NULL) return NULL;	
  
  /* check for non-IUPAC characters */
  for (cp = line; *cp != 0; cp++)
  {
    ch = TO_LOWER (*cp);
  	if (isspace (ch))
  	{
  	  /* space allowed */
  	}
  	else if (StringChr (valid_iupac_characters, ch) == NULL)
  	{
  	  ans = Message (MSG_YN, "Line %d (%s) contains invalid IUPAC characters "
  	                 "(should contain only %s).  Is this a definition line?  "
  	                 "If not, you will need to remove or replace the invalid "
  	                 "characters in the sequence later.",
  	                 line_num, line, valid_iupac_characters);
  	  if (ans == ANS_YES)
  	  {
  	  	return ParseImportFastaDefLine (line, line_num, cfp, cancelled, qualfixes, to_remove);
  	  }
  	  else
  	  {
  	  	break;
  	  }
  	}
  }

  if (cfsp == NULL)
  {	
    cfp->sequence_list = AddSequenceGroup (cfp->sequence_list, 
                                           cfp->selected_qual_list,
                                           cfp->default_tax_name,
                                           TRUE,
                                           to_remove,
                                           cfp);
    vnp = cfp->sequence_list;
    while (vnp != NULL && vnp->next != NULL)
    {
      vnp = vnp->next;
    }
    if (vnp != NULL)
    {
      cfsp = (CreateFastaSeqPtr) vnp->data.ptrvalue;
    }
  }
  if (cfsp == NULL) return NULL;
  
  new_len = StringLen (line) + StringLen (cfsp->sequence) + 1;
  new_seq = (CharPtr) MemNew (new_len * sizeof (Char));
  if (new_seq != NULL)
  {
  	StringCpy (new_seq, cfsp->sequence);
  	StringCat (new_seq, line);
  	MemFree (cfsp->sequence);
  	cfsp->sequence = new_seq;
  	cfsp->sequence = ReformatSequenceText (cfsp->sequence);
  }
  return cfsp;
}

static void ImportFastaCreationFile (CreateFastaPtr cfp, BoolPtr cancelled, BoolPtr segmented)
{
  CharPtr           extension;
  Char              path [PATH_MAX];
  ReadBufferData    rbd;
  CreateFastaSeqPtr cfsp = NULL;
  Int4              line_num = 1;
  Int4              last;
  CharPtr           line;
  WindoW            w;
  Int4              current_selection = 1;
  Int4              to_remove = -1;
  Boolean           remove_current = FALSE;
  ValNodePtr        vnp, prev;
  Int4              val;
  ValNodePtr        qualfixes = NULL;

  if (cfp == NULL) return;
  
  if (cfp->fpp->single)
  {
    to_remove = 1;
  	remove_current = TRUE;
  }
  else if (cfp->sequence_list != NULL)
  {
    current_selection = GetValue (cfp->sequence_selector);
    to_remove = current_selection;
    remove_current = IsThisSequenceEmpty (cfp, to_remove);	
  }
  if (!remove_current)
  {
  	to_remove = -1;
  }

  if (cfp->sequence_list != NULL && ! ValidateSequenceInformation (cfp->sequence_list, to_remove)) 
  {
    return; 
  }
  GetMasterQualList (cfp);
  last = ValNodeLen (cfp->sequence_list);

  extension = GetAppProperty ("FastaNucExtension");
  if (! GetInputFileName (path, sizeof (path), extension, "TEXT")) return;
  rbd.fp = FileOpen (path, "r");
  if (rbd.fp == NULL) return;
  
  rbd.current_data = NULL;
  line = AbstractReadFunction (&rbd);
  while (line != NULL && ! *cancelled)
  {
    if (line [0] == '[')
    {
      *cancelled = TRUE;
      *segmented = TRUE;
      if (cfp->form != NULL)
      {
      	Message (MSG_ERROR, "Cannot import segmented sequences");
      }
    }
    else if (line [0] == '#')
    {
      cfsp = NULL;
    }
    else if (StringHasNoText (line))
    {
      /* finish current sequence */
      cfsp = NULL;
    }
    else if (line [0] == '>')
    {
      cfsp = ParseImportFastaDefLine (line, line_num, cfp, cancelled, &qualfixes, to_remove);
    }
    else
    {
      /* parse as sequence line */
      cfsp = ParseImportFastaSeqLine (cfsp, line, line_num, cfp, cancelled, &qualfixes, to_remove);
    }
    line = MemFree (line);
    line_num ++;
    line = AbstractReadFunction (&rbd);
  }
  /* if sequence in progress, finish sequence */
  FileClose (rbd.fp);
  qualfixes = FreeStringPairList (qualfixes);  
  
  if (*cancelled)
  {
  	/* free newly created sequences */
    prev = NULL;
    for (val = 1, vnp = cfp->sequence_list; vnp != NULL && val <= last; vnp = vnp->next, val++)
    {
  	  prev = vnp;
    }
    if (vnp != NULL) 
    {    	
      if (prev == NULL || last == 0)
      {
  	    cfp->sequence_list = NULL;
      }
      else 
      {
  	    prev->next = NULL;
      }
      FreeSequenceGroup (vnp);  	
    }  	
    last = current_selection;
  }
  else
  { 	
    /* remove current sequence if empty */
    if (remove_current)
    {
      prev = NULL;
      for (val = 1, vnp = cfp->sequence_list; vnp != NULL && val < to_remove; vnp = vnp->next, val++)
      {
  	    prev = vnp;
      }
      if (vnp == NULL) return;
      if (prev == NULL)
      {
  	    cfp->sequence_list = vnp->next;
      }
      else 
      {
  	    prev->next = vnp->next;
      }
      vnp->next = NULL;
      FreeSequenceGroup (vnp);  	
    }

    if (remove_current && to_remove < last)
    {
      last --;
    }
    else if (!remove_current)
    {
  	  last++;
    }
  
    if (cfp->fpp->single && cfp->sequence_list->next != NULL)
    {
  	  Message (MSG_ERROR, "Only one sequence can be imported!");
  	  vnp = cfp->sequence_list->next;
  	  cfp->sequence_list->next = NULL;
  	  FreeSequenceGroup (vnp);
    }
  }

  if (cfp->form == NULL) return;
  
  SetObjectExtra (cfp->form, NULL, NULL);
  Remove (cfp->form);
  w = CreateFastaWindow (cfp);
  cfp->form = (ForM) w;
  SetValue (cfp->sequence_selector, last);
  ShowSequenceGroup (cfp->sequence_selector);
  RealizeWindow (w);
  Show (w);
  Select (w);
}

static void ImportFastaCreationFileItem (IteM i)
{
  CreateFastaPtr cfp;
  Boolean        cancelled = FALSE;
  Boolean        segmented = FALSE;
  
  cfp = (CreateFastaPtr) GetObjectExtra (i);
  ImportFastaCreationFile (cfp, &cancelled, &segmented);
}

static WindoW CreateFastaWindow (CreateFastaPtr cfp)
{
  WindoW w;
  GrouP  h, g, c1, c2, c3;
  ButtoN a;
  MenU   m;
  IteM   i;
  
  if (cfp == NULL) return NULL;
  
  GetMasterQualList (cfp);
  
  w = FixedWindow (-50, -33, -10, -10, "FASTA File", NULL);
  SetObjectExtra (w, cfp, CleanupCreateFastaForm);
  
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 20);
  
  g = HiddenGroup (h, 1, 0, NULL);
  DrawSequenceGroups (cfp, g);
  if (!cfp->fpp->single)
  {
    c1 = HiddenGroup (h, 2, 0, NULL);
    a = PushButton (c1, "Add Sequence", AddFastaSequenceButton);
    SetObjectExtra (a, cfp, NULL);
    a = PushButton (c1, "Remove Sequence", RemoveFastaSequenceButton);
    SetObjectExtra (a, cfp, NULL);
  
    c2 = HiddenGroup (h, 4, 0, NULL);
    DrawSequenceSelector (cfp, c2);
    if (cfp->sequence_selector != NULL)
    {
      SetValue (cfp->sequence_selector, 1);
  	  ShowSequenceGroup (cfp->sequence_selector);
    }
  }
  c3 = HiddenGroup (h, 7, 0, NULL);
  SetGroupSpacing (c3, 10, 3);
  a = DefaultButton (c3, "Finish FASTA File", DoCreateFASTAFileButton);
  SetObjectExtra (a, cfp, NULL);
  if (! cfp->fpp->single && ValNodeLen (cfp->sequence_list) < 2)
  {
  	Disable (a);
  }
  a = PushButton (c3, "Cancel", CreateFastaCancelButton); 
  SetObjectExtra (a, cfp, NULL);
  a = PushButton (c3, "Preview File", PreviewFastaFileButton);
  SetObjectExtra (a, cfp, NULL);

  if (cfp->fpp->single)
  {
    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c3, NULL);
  }
  else
  {
    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c1, (HANDLE) c2, (HANDLE) c3, NULL);
  }


  /* add menus */  
  m = PulldownMenu (w, "File");
  i = CommandItem (m, "Import sequence data from file", ImportFastaCreationFileItem);
  SetObjectExtra (i, cfp, NULL);
  i = CommandItem (m, "Finish FASTA File", DoCreateFASTAFileItem);
  SetObjectExtra (i, cfp, NULL);
  i = CommandItem (m, "Cancel", CreateFastaCancelItem);
  SetObjectExtra (i, cfp, NULL);
  i = CommandItem (m, "Preview File", PreviewFastaFileItem);
  SetObjectExtra (i, cfp, NULL);
  if (!cfp->fpp->single)
  {
    m = PulldownMenu (w, "Edit");
    i = CommandItem (m, "Add Sequence", AddFastaSequenceItem);
    SetObjectExtra (i, cfp, NULL);
    i = CommandItem (m, "Remove Sequence", RemoveFastaSequenceItem);
    SetObjectExtra (i, cfp, NULL);
  }
  m = PulldownMenu (w, "Help");
  i = CommandItem (m, "Sequence ID", SequenceIdHelp);
  i = CommandItem (m, "Organism Name", OrganismNameHelp);
  i = CommandItem (m, "Organism Qualifiers", OrganismQualifiersHelp);
  i = CommandItem (m, "Sequence Characters", SequenceCharactersHelp);
  return w;
}

static void AddFastaSequence (CreateFastaPtr   cfp)
{
  WindoW           w;
  Int4             last;
  
  if (cfp == NULL) return;
  
  if (cfp->sequence_list != NULL && ! ValidateSequenceInformation (cfp->sequence_list, -1))
  {
  	return;
  }
  GetMasterQualList (cfp);

  SetObjectExtra (cfp->form, NULL, NULL);
  Remove (cfp->form);
  cfp->sequence_list = AddSequenceGroup (cfp->sequence_list, 
                                         cfp->selected_qual_list,
                                         cfp->default_tax_name,
                                         cfp->auto_id,
                                         -1,
                                         cfp);

  w = CreateFastaWindow (cfp);
  cfp->form = (ForM) w;
  last = ValNodeLen (cfp->sequence_list);
  SetValue (cfp->sequence_selector, last);
  ShowSequenceGroup (cfp->sequence_selector);
  RealizeWindow (w);
  Show (w);
  Select (w);
}

static void RemoveFastaSequence (CreateFastaPtr   cfp)
{
  WindoW           w;
  Int4             to_remove, val;
  ValNodePtr       prev, vnp;
  
  if (cfp == NULL) return;

  to_remove = GetValue (cfp->sequence_selector);
  
  ValidateSequenceInformation (cfp->sequence_list, to_remove);
  SetObjectExtra (cfp->form, NULL, NULL);
  
  prev = NULL;
  for (val = 1, vnp = cfp->sequence_list; vnp != NULL && val < to_remove; vnp = vnp->next, val++)
  {
  	prev = vnp;
  }
  if (vnp == NULL) return;
  if (prev == NULL)
  {
  	cfp->sequence_list = vnp->next;
  }
  else 
  {
  	prev->next = vnp->next;
  }
  vnp->next = NULL;
  FreeSequenceGroup (vnp);
  
  Remove (cfp->form);

  w = CreateFastaWindow (cfp);
  cfp->form = (ForM) w;
  
  if (to_remove != 1)
  {
  	to_remove--;
  }
  SetValue (cfp->sequence_selector, to_remove);
  ShowSequenceGroup (cfp->sequence_selector);
  RealizeWindow (w);
  Show (w);
  Select (w);
}

static void CreateFASTAFile (ButtoN b)
{
  CreateFastaPtr   cfp;
  WindoW           w;
  SequencesFormPtr sqfp;
  FastaPagePtr     fpp;
  BioSourcePtr     biop;
  
  sqfp = (SequencesFormPtr) GetObjectExtra (b);
  if (sqfp == NULL || sqfp->dnaseq == NULL) return;
  fpp = (FastaPagePtr) GetObjectExtra (sqfp->dnaseq);
  if (fpp == NULL) return;
  
  cfp = (CreateFastaPtr) MemNew (sizeof (CreateFastaData));
  if (cfp == NULL) return;
  cfp->fpp = fpp;
  cfp->use_id_from_fasta_defline = sqfp->use_id_from_fasta_defline;
  cfp->auto_id = TRUE;
  
  if (sqfp->seqPackage == SEQ_PKG_SEGMENTED)
  {
    cfp->is_segmented = TRUE;
    cfp->is_popset = FALSE;
  }
  else if (sqfp->seqPackage == SEQ_PKG_POPULATION)
  {
    cfp->is_popset = TRUE;
    cfp->is_segmented = FALSE;
  }
  else
  {
    cfp->is_popset = FALSE;
    cfp->is_segmented = FALSE;
  }
  
  biop = (BioSourcePtr) DialogToPointer (sqfp->genbio);
  if (biop != NULL && biop->org != NULL && !StringHasNoText (biop->org->taxname))
  {
  	cfp->default_tax_name = MemNew (sizeof (Char) * (StringLen (biop->org->taxname) + 1));
  	StringCpy (cfp->default_tax_name, biop->org->taxname);
  }
  
  cfp->sequence_list = AddSequenceGroup (NULL, cfp->selected_qual_list, 
                                         cfp->default_tax_name, cfp->auto_id,
                                         -1, cfp);  
  
  w = CreateFastaWindow (cfp);
  cfp->form = (ForM) w;

  SendHelpScrollMessage (helpForm, "Organism and Sequences Form", "FASTA Format for Nucleotide Sequences");
  
  RealizeWindow (w);
  Show (w);
  Select (w);
}

static void ClearSequencesButton (ButtoN b)
{
  SequencesFormPtr   sqfp;
  
  sqfp = (SequencesFormPtr) GetObjectExtra (b);
  if (sqfp == NULL) return;
  SequencesFormDeleteProc (sqfp);
}

extern ForM CreateInitOrgNucProtForm (Int2 left, Int2 top, CharPtr title,
                                      FormatBlockPtr format,
                                      BtnActnProc goToNext,
                                      BtnActnProc goBack,
                                      WndActnProc activateForm)

{
  ButtoN             b, b2 = NULL;
  GrouP              c;
  GrouP              f1, f2, f3;
  GrouP              g;
  GrouP              h;
  Handle             h1, h2, h3, h4;
  GrouP              j;
  GrouP              k;
  ButtoN             mrna;
  GrouP              mult;
  GrouP              p;
  Int2               page;
  Boolean            parseSeqId;
  PrompT             ppt1, ppt2;
  ButtoN             prs;
  GrouP              q;
  StdEditorProcsPtr  sepp;
  Boolean            single;
  SequencesFormPtr   sqfp;
  Char               str [32];
  WindoW             w;
  GrouP              x;
  GrouP              y;
  GrouP              z;
  GrouP              import_btn_grp;
  FastaPagePtr       fpp;

  w = NULL;
  sqfp = MemNew (sizeof (SequencesForm));
  if (sqfp != NULL) {

    if (format != NULL) {
      sqfp->seqPackage = format->seqPackage;
      sqfp->seqFormat = format->seqFormat;
      sqfp->numSeqs = format->numSeqs;
      sqfp->submType = format->submType;
    } else {
      sqfp->seqPackage = SEQ_PKG_SINGLE;
      sqfp->seqFormat = SEQ_FMT_FASTA;
      sqfp->numSeqs = 0;
      sqfp->submType = SEQ_ORIG_SUBMISSION;
    }

    w = FixedWindow (left, top, -10, -10, title, NULL);
    SetObjectExtra (w, sqfp, StdCleanupFormProc);
    sqfp->form = (ForM) w;
    sqfp->toform = NULL;
    if (sqfp->seqFormat == SEQ_FMT_FASTA) {
      sqfp->fromform = FastaSequencesFormToSeqEntryPtr;
      sqfp->testform = FastaTestSequencesForm;
    } else if (sqfp->seqFormat == SEQ_FMT_ALIGNMENT) {
      sqfp->fromform = PhylipSequencesFormToSeqEntryPtr;
      sqfp->testform = PhylipTestSequencesForm;
    }
    sqfp->importform = ImportSequencesForm;
    sqfp->formmessage = SequencesFormMessage;

#ifndef WIN_MAC
    CreateSqnInitialFormMenus (w);
#endif

    sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
    if (sepp != NULL) {
      sqfp->appmessage = sepp->handleMessages;
    }

    SetGroupSpacing (w, 10, 10);

    j = HiddenGroup (w, 10, 0, NULL);

    if (sqfp->seqPackage <= SEQ_PKG_SEGMENTED) {
      sqfp->tbs = CreateFolderTabs (j, seqSegFormTabs, ORGANISM_PAGE,
                                    0, 0, SYSTEM_FOLDER_TAB,
                                    ChangeSequencesPage, (Pointer) sqfp);
    } else if (sqfp->seqPackage == SEQ_PKG_GENOMICCDNA) {
      sqfp->tbs = CreateFolderTabs (j, cdnaGenFormTabs, ORGANISM_PAGE,
                                    0, 0, SYSTEM_FOLDER_TAB,
                                    ChangeSequencesPage, (Pointer) sqfp);
    } else {
      sqfp->tbs = CreateFolderTabs (j, popPhyMutFormTabs, ORGANISM_PAGE,
                                    0, 0, SYSTEM_FOLDER_TAB,
                                    ChangeSequencesPage, (Pointer) sqfp);
    }
    sqfp->currentPage = 0;
    page = 0;

    h = HiddenGroup (w, 0, 0, NULL);

    q = HiddenGroup (h, 0, 0, NULL);
    SetGroupSpacing (q, 10, 20);
    sqfp->genbio = CreateSimpleBioSourceDialog (q, "");
    if (sqfp->seqPackage == SEQ_PKG_PHYLOGENETIC) {
      p = HiddenGroup (q, -1, 0, NULL);
      SetGroupSpacing (p, 10, 20);
      if (sqfp->seqFormat == SEQ_FMT_FASTA) {
        mult = MultiLinePrompt (p, phyloOrgFastaMsg, 27 * stdCharWidth, programFont);
      } else if (sqfp->seqFormat == SEQ_FMT_ALIGNMENT) {
        mult = MultiLinePrompt (p, phyloOrgPhylipMsg, 27 * stdCharWidth, programFont);
      } else {
        mult = MultiLinePrompt (p, phyloOrgFastaMsg, 27 * stdCharWidth, programFont);
      }
      f3 = HiddenGroup (p, -1, 0, NULL);
      SetGroupSpacing (f3, 10, 10);

      f1 = HiddenGroup (f3, 3, 0, NULL);
      StaticPrompt (f1, "Location of Sequence",
                    0, popupMenuHeight, programFont, 'l');
      sqfp->genome = PopupList (f1, TRUE, SetGenome);
      SetObjectExtra (sqfp->genome, GetObjectExtra (sqfp->genbio), NULL);
      InitEnumPopup (sqfp->genome, biosource_genome_simple_alist, NULL);
      SetValue (sqfp->genome, 2);
      ReplaceBioSourceGenomePopup (sqfp->genbio, sqfp->genome);

      f2 = HiddenGroup (f3, 3, 0, NULL);
      StaticPrompt (f2, "Genetic Code for Translation", 0, popupMenuHeight, programFont, 'l');
      sqfp->gencode = PopupList (f2, TRUE, NULL);
      PopulateGeneticCodePopup (sqfp->gencode);
      SetValue (sqfp->gencode, 1);

      AlignObjects (ALIGN_CENTER, (HANDLE) mult, (HANDLE) f3, NULL);
      Hide (sqfp->genbio);
    }
    sqfp->pages [page] = q;
    Hide (sqfp->pages [page]);
    sqfp->tagFromPage [page] = ORGANISM_PAGE;
    page++;

    q = HiddenGroup (h, -1, 0, NULL);
    SetGroupSpacing (q, 10, 20);
    g = HiddenGroup (q, -1, 0, NULL);
    SetGroupSpacing (g, 10, 10);
    x = HiddenGroup (g, -4, 0, NULL);
    StaticPrompt (x, "Molecule", 0, popupMenuHeight, programFont, 'l');
    sqfp->moltypePopup = PopupList (x, TRUE, NULL);
    if (sqfp->seqPackage == SEQ_PKG_GENOMICCDNA) {
      sqfp->moltypeAlist = biomol_nucGen_alist;
    } else {
      sqfp->moltypeAlist = biomol_nucX_alist;
    }
    InitEnumPopup (sqfp->moltypePopup, sqfp->moltypeAlist, NULL);
    if (sqfp->seqPackage >= SEQ_PKG_SEGMENTED) {
      SetEnumPopup (sqfp->moltypePopup, sqfp->moltypeAlist, 253);
    } else {
      SetEnumPopup (sqfp->moltypePopup, sqfp->moltypeAlist, 3);
    }
    StaticPrompt (x, "Topology", 0, popupMenuHeight, programFont, 'c');
    sqfp->topologyPopup = PopupList (x, TRUE, NULL);
    InitEnumPopup (sqfp->topologyPopup, topology_nuc_alist, NULL);
    SetEnumPopup (sqfp->topologyPopup, topology_nuc_alist, TOPOLOGY_LINEAR);
    y = HiddenGroup (g, -2, 0, NULL);
    SetGroupSpacing (y, 10, 2);
    sqfp->partial5 = CheckBox (y, "Incomplete at 5' end", NULL);
    sqfp->partial3 = CheckBox (y, "Incomplete at 3' end", NULL);
    prs = NULL;
    if (sqfp->seqFormat == SEQ_FMT_FASTA) {
      prs = CheckBox (g, "Fasta definition line starts with sequence ID", ChangeDnaParse);
      sqfp->use_id_from_fasta_defline = prs;
      SetObjectExtra (prs, sqfp, NULL);
      parseSeqId = FALSE;
      if (GetAppParam ("SEQUIN", "PREFERENCES", "PARSENUCSEQID", NULL, str, sizeof (str))) {
        if (StringICmp (str, "TRUE") == 0) {
          parseSeqId = TRUE;
        }
      }
      SetStatus (prs, parseSeqId);
      if (sqfp->seqPackage == SEQ_PKG_SINGLE) {
        sqfp->singleIdGrp = HiddenGroup (g, 2, 0, NULL);
        StaticPrompt (sqfp->singleIdGrp, "Enter unique identifier for this sequence", 0, dialogTextHeight, programFont, 'l');
        sqfp->singleSeqID = DialogText (sqfp->singleIdGrp, "", 6, NULL);
        if (parseSeqId) {
          Hide (sqfp->singleIdGrp);
        }
      }
    }
    if (((sqfp->seqPackage == SEQ_PKG_POPULATION) ||
	 (sqfp->seqPackage == SEQ_PKG_PHYLOGENETIC) ||
	 (sqfp->seqPackage == SEQ_PKG_MUTATION) ||
	 (sqfp->seqPackage == SEQ_PKG_ENVIRONMENT)) &&
	(sqfp->seqFormat == SEQ_FMT_FASTA)) {
      sqfp->makeAlign = CheckBox (g, "Create Alignment", NULL);
      /*
      if (sqfp->seqPackage < SEQ_PKG_GENBANK) {
        SetStatus (sqfp->makeAlign, TRUE);
      }
      */
    }

    b = NULL;
    k = HiddenGroup (g, 0, 2, NULL);
    if (sqfp->seqFormat == SEQ_FMT_FASTA) {
      single = (Boolean) (sqfp->seqPackage == SEQ_PKG_SINGLE);
      if (sqfp->seqPackage == SEQ_PKG_GENOMICCDNA) {
        sqfp->dnaseq = CreateFastaDialog (k, "", TRUE, FALSE, fastaGenMsg, parseSeqId, single);
        fpp = (FastaPagePtr) GetObjectExtra (sqfp->dnaseq);
        import_btn_grp = HiddenGroup (g, 4, 0, NULL);
        fpp->import_btn = PushButton (import_btn_grp, "Import Genomic FASTA", ImportBtnProc);
      } else {
        sqfp->dnaseq = CreateFastaDialog (k, "", TRUE, FALSE, fastaNucMsg, parseSeqId, single);
        fpp = (FastaPagePtr) GetObjectExtra (sqfp->dnaseq);
        import_btn_grp = HiddenGroup (g, 4, 0, NULL);
        fpp->import_btn = PushButton (import_btn_grp, "Import Nucleotide FASTA", ImportBtnProc);
      }
      SetObjectExtra (fpp->import_btn, sqfp, NULL);
#ifdef USE_CREATE_MY_FASTA
      fpp->create_btn = PushButton (import_btn_grp, "Create My FASTA File", CreateFASTAFile);
      SetObjectExtra (fpp->create_btn, sqfp, NULL);
#endif
      fpp->clear_btn = PushButton (import_btn_grp, "Clear sequences", ClearSequencesButton);
      SetObjectExtra (fpp->clear_btn, sqfp, NULL);
      Disable (fpp->clear_btn);
    } else if (sqfp->seqFormat == SEQ_FMT_ALIGNMENT) {
      sqfp->dnaseq = CreatePhylipDialog (k, "", phylipNucMsg, sqfp->seqFormat, "",
                                         sqfp->seqPackage, sqfp->genbio);
      import_btn_grp = HiddenGroup (g, 4, 0, NULL);
      import_btn_grp = HiddenGroup (g, 4, 0, NULL);
      b = PushButton (import_btn_grp, "Import Nucleotide Alignment", ImportBtnProc);
      SetObjectExtra (b, sqfp, NULL);
    }
    if (sqfp->makeAlign != NULL) {
      h1 = (Handle) sqfp->makeAlign;
      h2 = (Handle) import_btn_grp;
      h3 = (Handle) prs;
      h4 = (Handle) sqfp->singleIdGrp;
    } else {
      h1 = import_btn_grp;
      h2 = (Handle) prs;
      h3 = (Handle) sqfp->singleIdGrp;
      h4 = NULL;
    }
    AlignObjects (ALIGN_CENTER, (HANDLE) x, (HANDLE) y, (HANDLE) k,
                  (HANDLE) h1, (HANDLE) h2, (HANDLE) h3, (HANDLE) h4, NULL);
    sqfp->pages [page] = q;
    Hide (sqfp->pages [page]);
    sqfp->tagFromPage [page] = NUCLEOTIDE_PAGE;
    page++;

    if (sqfp->seqPackage == SEQ_PKG_GENOMICCDNA) {
      q = HiddenGroup (h, -1, 0, NULL);
      SetGroupSpacing (q, 10, 20);
      g = HiddenGroup (q, -1, 0, NULL);
      SetGroupSpacing (g, 10, 10);
      y = HiddenGroup (g, -2, 0, NULL);
      SetGroupSpacing (y, 10, 2);
      sqfp->partialmRNA5 = CheckBox (y, "Incomplete at 5' end", NULL);
      sqfp->partialmRNA3 = CheckBox (y, "Incomplete at 3' end", NULL);
      prs = CheckBox (g, "Fasta definition line starts with sequence ID", ChangeMrnaParse);
      SetObjectExtra (prs, sqfp, NULL);
      parseSeqId = FALSE;
      if (GetAppParam ("SEQUIN", "PREFERENCES", "PARSEMRNASEQID", NULL, str, sizeof (str))) {
        if (StringICmp (str, "TRUE") == 0) {
          parseSeqId = TRUE;
        }
      }
      SetStatus (prs, parseSeqId);

      k = HiddenGroup (g, 0, 2, NULL);
      sqfp->mrnaseq = CreateFastaDialog (k, "", TRUE, TRUE, fastaMrnaMsg, parseSeqId, FALSE);
      b = PushButton (g, "Import Transcript FASTA", ImportBtnProc);
      SetObjectExtra (b, sqfp, NULL);

      AlignObjects (ALIGN_CENTER, (HANDLE) y, (HANDLE) k, (HANDLE) prs, (HANDLE) b, NULL);
      sqfp->pages [page] = q;
      Hide (sqfp->pages [page]);
      sqfp->tagFromPage [page] = MRNA_PAGE;
      page++;
    }

    if (sqfp->seqPackage <= SEQ_PKG_GENOMICCDNA) {
      q = HiddenGroup (h, -1, 0, NULL);
      SetGroupSpacing (q, 10, 20);
      g = HiddenGroup (q, -1, 0, NULL);
      SetGroupSpacing (g, 10, 10);
      x = HiddenGroup (g, -1, 0, NULL);
      sqfp->protTechBoth = CheckBox (x,
             "Conceptual translation confirmed by peptide sequencing", NULL);
      y = HiddenGroup (g, -2, 0, NULL);
      SetGroupSpacing (y, 10, 2);
      sqfp->partialN = CheckBox (y, "Incomplete at NH2 end", NULL);
      sqfp->partialC = CheckBox (y, "Incomplete at CO2H end", NULL);
      prs = CheckBox (g, "Fasta definition line starts with sequence ID", ChangeProtParse);
      SetObjectExtra (prs, sqfp, NULL);
      parseSeqId = FALSE;
      if (GetAppParam ("SEQUIN", "PREFERENCES", "PARSEPROTSEQID", NULL, str, sizeof (str))) {
        if (StringICmp (str, "TRUE") == 0) {
          parseSeqId = TRUE;
        }
      }
      SetStatus (prs, parseSeqId);
      sqfp->makeMRNA = FALSE;
      mrna = NULL;
      if (sqfp->seqPackage != SEQ_PKG_GENOMICCDNA) {
        mrna = CheckBox (g, "Create initial mRNA with CDS intervals", ChangeMrnaFlag);
        SetObjectExtra (mrna, sqfp, NULL);
        if (GetAppParam ("SEQUIN", "PREFERENCES", "CREATEMRNA", NULL, str, sizeof (str))) {
          if (StringICmp (str, "TRUE") == 0) {
            sqfp->makeMRNA = TRUE;
          }
        }
      }
      SetStatus (mrna, sqfp->makeMRNA);
      k = HiddenGroup (g, 0, 2, NULL);
      sqfp->protseq = CreateFastaDialog (k, "", FALSE, FALSE, fastaProtMsg, parseSeqId, FALSE);
      b = PushButton (g, "Import Protein FASTA", ImportBtnProc);
      SetObjectExtra (b, sqfp, NULL);
      AlignObjects (ALIGN_CENTER, (HANDLE) x, (HANDLE) y, (HANDLE) k,
                    (HANDLE) prs, (HANDLE) b, (HANDLE) mrna, NULL);
      sqfp->pages [page] = q;
      Hide (sqfp->pages [page]);
      sqfp->tagFromPage [page] = PROTEIN_PAGE;
      page++;
    } else if (sqfp->seqPackage >= SEQ_PKG_POPULATION) {
      q = HiddenGroup (h, -1, 0, NULL);
      SetGroupSpacing (q, 10, 10);
      ppt1 = StaticPrompt (q, "Add feature across full length of all sequences",
                          0, 0, programFont, 'l');
      sqfp->annotType = HiddenGroup (q, 5, 0, ChangeAnnotType);
      SetObjectExtra (sqfp->annotType, sqfp, NULL);
      RadioButton (sqfp->annotType, "Gene");
      RadioButton (sqfp->annotType, "rRNA");
      RadioButton (sqfp->annotType, "CDS");
      RadioButton (sqfp->annotType, "None");
      SetValue (sqfp->annotType, 1);
      sqfp->annotGrp = HiddenGroup (q, -1, 0, NULL);
      SetGroupSpacing (sqfp->annotGrp, 10, 10);
      x = HiddenGroup (sqfp->annotGrp, 2, 0, NULL);
      sqfp->partialLft = CheckBox (x, "Incomplete at 5' end", NULL);
      sqfp->partialRgt = CheckBox (x, "Incomplete at 3' end", NULL);
      y = HiddenGroup (sqfp->annotGrp, 2, 0, NULL);
      sqfp->protOrRnaPpt = StaticPrompt (y, "Protein Name", 0, dialogTextHeight, programFont, 'l');
      sqfp->protOrRnaName = DialogText (y, "", 20, NULL);
      sqfp->protDescPpt = StaticPrompt (y, "Protein Description", 0, dialogTextHeight, programFont, 'l');
      sqfp->protDesc = DialogText (y, "", 20, NULL);
      StaticPrompt (y, "Gene Symbol", 0, dialogTextHeight, programFont, 'l');
      sqfp->geneName = DialogText (y, "", 20, NULL);
      StaticPrompt (y, "Comment", 0, 3 * Nlm_stdLineHeight, programFont, 'l');
      sqfp->featcomment = ScrollText (y, 20, 3, programFont, TRUE, NULL);
      ppt2 = StaticPrompt (q, "Add title to all sequences if not in definition line",
                           0, 0, programFont, 'c');
      z = HiddenGroup (q, 2, 0, NULL);
      StaticPrompt (z, "Title       ", 0, 3 * Nlm_stdLineHeight, programFont, 'c');
      sqfp->defline = ScrollText (z, 20, 3, programFont, TRUE, NULL);
      sqfp->orgPrefix = CheckBox (q, "Prefix title with organism name", NULL);
      AlignObjects (ALIGN_CENTER, (HANDLE) ppt1, (HANDLE) sqfp->annotType,
                    (HANDLE) x, (HANDLE) y, (HANDLE) ppt2, (HANDLE) z,
                    (HANDLE) sqfp->orgPrefix, NULL);
      Hide (sqfp->protOrRnaPpt);
      Hide (sqfp->protOrRnaName);
      Hide (sqfp->protDescPpt);
      Hide (sqfp->protDesc);
      /* Hide (sqfp->annotGrp); */
      sqfp->pages [page] = q;
      Hide (sqfp->pages [page]);
      sqfp->tagFromPage [page] = ANNOTATE_PAGE;
      page++;
    }
    sqfp->numPages = page;

    c = HiddenGroup (w, 3, 0, NULL);
    SetGroupSpacing (c, 10, 2);
    sqfp->goToPrev = goBack;
    sqfp->prevBtn = PushButton (c, " << Prev Form ", PrevSequencesFormBtn);
    SetObjectExtra (sqfp->prevBtn, sqfp, NULL);
    sqfp->goToNext = goToNext;
    sqfp->nextBtn = PushButton (c, " Next Page >> ", NextSequencesFormBtn);
    SetObjectExtra (sqfp->nextBtn, sqfp, NULL);

    AlignObjects (ALIGN_CENTER, (HANDLE) j, (HANDLE) c,
                  (HANDLE) sqfp->pages [0], (HANDLE) sqfp->pages [1],
                  (HANDLE) sqfp->pages [2], (HANDLE) sqfp->pages [3], NULL);

    RealizeWindow (w);

    SafeSetTitle (sqfp->prevBtn, "<< Prev Form");
    SafeSetTitle (sqfp->nextBtn, "Next Page >>");

    sqfp->activate = activateForm;
    SetActivate (w, InitOrgNucProtFormActivate);

    SendMessageToDialog (sqfp->tbs, VIB_MSG_INIT);
    SendMessageToDialog (sqfp->genbio, VIB_MSG_INIT);
    SendMessageToDialog (sqfp->dnaseq, VIB_MSG_INIT);
    SendMessageToDialog (sqfp->protseq, VIB_MSG_INIT);

    Show (sqfp->pages [sqfp->currentPage]);
    SendMessageToDialog (sqfp->genbio, VIB_MSG_ENTER);
  }
  return (ForM) w;
}

static void MakePubAndDefLine (SequinBlockPtr sbp, SeqEntryPtr sep)

{
  AffilPtr     affil;
  AuthListPtr  alp;
  CitGenPtr    cgp;
  PubdescPtr   pdp;
  ValNodePtr   pep;
  ValNodePtr   vnp;
  /*
  BioseqSetPtr  bssp;
  Char          str [256];
  CharPtr       ttl;
  */

  if (sep == NULL) return;
  /*
  if (SeqEntryGetTitle (sep) != NULL) return;
  ttl = NULL;
  SeqEntryExplore (sep, (Pointer) (&ttl), FindFirstTitle);
  if (ttl != NULL) {
    vnp = CreateNewDescriptor (sep, Seq_descr_title);
    if (vnp != NULL) {
      StringNCpy_0 (str, ttl, sizeof (str) - 32);
      if (IS_Bioseq_set (sep)) {
        bssp = (BioseqSetPtr) sep->data.ptrvalue;
        if (bssp != NULL && bssp->_class == 1) {
          StringCat (str, ", and translated products");
        }
      }
      vnp->data.ptrvalue = StringSave (str);
    }
  }
  */
  if (sbp == NULL || sbp->citsubauthors == NULL) return;
  pdp = PubdescNew ();
  if (pdp != NULL) {
    vnp = CreateNewDescriptor (sep, Seq_descr_pub);
    if (vnp != NULL) {
      vnp->data.ptrvalue = (Pointer) pdp;
      pdp->reftype = 0;
      pep = ValNodeNew (NULL);
      pdp->pub = pep;
      if (pep != NULL) {
        cgp = CitGenNew ();
        if (cgp != NULL) {
          pep->choice = PUB_Gen;
          pep->data.ptrvalue = cgp;
          cgp->cit = StringSave ("unpublished");
          alp = AsnIoMemCopy ((Pointer) sbp->citsubauthors,
                              (AsnReadFunc) AuthListAsnRead,
                              (AsnWriteFunc) AuthListAsnWrite);
          cgp->authors = alp;
          if (alp != NULL) {
            affil = AsnIoMemCopy ((Pointer) sbp->citsubaffil,
                                  (AsnReadFunc) AffilAsnRead,
                                  (AsnWriteFunc) AffilAsnWrite);
            alp->affil = affil;
            if (affil != NULL) {
              affil->phone = MemFree (affil->phone);
              affil->fax = MemFree (affil->fax);
              affil->email = MemFree (affil->email);
            }
          }
          cgp->title = sbp->citsubtitle;
          sbp->citsubtitle = NULL;
        }
      }
    }
  }
}

extern SubmitBlockPtr ConvertSequinBlockToSubmitBlock (SequinBlockPtr sqp);

extern SubmitBlockPtr ConvertSequinBlockToSubmitBlock (SequinBlockPtr sqp)

{
  AffilPtr        affil;
  AuthorPtr       ap;
  AuthListPtr     authors;
  ContactInfoPtr  cip;
  CitSubPtr       csp;
  DatePtr         dp;
  CharPtr         os;
  SubmitBlockPtr  sbp;
  Char            str [64];

  sbp = NULL;
  if (sqp != NULL) {
    sbp = SubmitBlockNew ();
    if (sbp != NULL) {
      sbp->subtype = 1;
      os = GetOpSysString ();
      if (os != NULL) {
        sprintf (str, "Sequin %s - %s", SEQUIN_APPLICATION, os);
      } else {
        sprintf (str, "Sequin %s", SEQUIN_APPLICATION);
      }
      sbp->tool = StringSave (str);
      MemFree (os);
      sbp->reldate = sqp->releasedate;
      dp = sbp->reldate;
      if (dp != NULL && dp->data [0] == 1 && dp->data [1] > 0) {
        if (dp->data [2] == 0) {
          dp->data [2] = 1;
        }
        if (dp->data [3] == 0) {
          switch (dp->data [2]) {
            case 4 :
            case 6 :
            case 9 :
            case 11 :
              dp->data [3] = 30;
              break;
            case 2 :
              dp->data [3] = 28;
              break;
            default :
              dp->data [3] = 31;
              break;
          }
        }
      }
      cip = ContactInfoNew ();
      if (cip != NULL) {
        ap = sqp->contactperson;
        cip->contact = ap;
        if (ap != NULL) {
          affil = sqp->citsubaffil;
          if (affil != NULL) {
            if (ap->affil != NULL) {
              affil->phone = MemFree (affil->phone);
              affil->fax = MemFree (affil->fax);
              affil->email = MemFree (affil->email);
              affil->phone = StringSave (ap->affil->phone);
              affil->fax = StringSave (ap->affil->fax);
              affil->email = StringSave (ap->affil->email);
              ap->affil = AffilFree (ap->affil);
            }
            ap->affil = affil;
          }
        }
      }
      sbp->contact = cip;
      csp = CitSubFromContactInfo (cip);
      sbp->cit = csp;
      if (csp != NULL) {
        authors = csp->authors;
        if (authors != NULL) {
          affil = authors->affil;
          authors->affil = NULL;
          csp->authors = AuthListFree (csp->authors);
          csp->authors = sqp->citsubauthors;
          authors = csp->authors;
          if (authors != NULL) {
            authors->affil = affil;
            if (affil != NULL) {
              affil->phone = MemFree (affil->phone);
              affil->fax = MemFree (affil->fax);
              affil->email = MemFree (affil->email);
            }
          }
        }
      }
      sbp->hup = sqp->holduntilpublished;
    }
    MemFree (sqp);
  }
  return sbp;
}

extern Uint2 PackageFormResults (SequinBlockPtr sbp, SeqEntryPtr sep, Boolean makePubAndDefLine)

{
  Uint2         entityID;
  SeqSubmitPtr  ssp;

  entityID = 0;
  if (sep != NULL) {
    if (sbp != NULL) {
      ssp = SeqSubmitNew ();
      if (ssp != NULL) {
        ssp->datatype = 1;
        ssp->data = (Pointer) sep;
        if (makePubAndDefLine) {
          MakePubAndDefLine (sbp, sep);
        }
        sbp->citsubtitle = MemFree (sbp->citsubtitle);
        ssp->sub = ConvertSequinBlockToSubmitBlock (sbp);
        ObjMgrConnect (OBJ_SEQENTRY, sep->data.ptrvalue, OBJ_SEQSUB, (Pointer) ssp);
        if (! ObjMgrRegister (OBJ_SEQSUB, (Pointer) ssp)) {
          ErrPostEx (SEV_ERROR, 0, 0, "ObjMgrRegister failed.");
        }
      } else {
        if (! ObjMgrRegister (OBJ_SEQENTRY, (Pointer) sep)) {
          ErrPostEx (SEV_ERROR, 0, 0, "ObjMgrRegister failed.");
        }
      }
    } else {
      if (! ObjMgrRegister (OBJ_SEQENTRY, (Pointer) sep)) {
        ErrPostEx (SEV_ERROR, 0, 0, "ObjMgrRegister failed.");
      }
    }
    if (EntrezASN1Detected (sep)) {
      ErrPostEx (SEV_WARNING, 0, 0, "This record was retrieved from Entrez.");
    }
    entityID = ObjMgrGetEntityIDForChoice (sep);
  }
  return entityID;
}

static void GetRawBsps (BioseqPtr bsp, Pointer userdata)

{
  ValNodePtr PNTR  head;

  if (bsp->repr != Seq_repr_raw) return;
  head = (ValNodePtr PNTR) userdata;
  ValNodeAddPointer (head, 0, (Pointer) bsp);
}

static void ParseInMoreProteinsCommon (IteM i, Boolean doSuggest)

{
  SeqEntryPtr  addhere;
  MsgAnswer    ans;
  BaseFormPtr  bfp;
  BioseqPtr    bsp;
  BioseqSetPtr  bssp;
  Int2         code;
  Int4         count;
  ValNodePtr   descr;
  CharPtr      errormsg;
  CharPtr      extension;
  FILE         *fp;
  ValNodePtr   head;
  Boolean      isLocalUnknownID;
  SeqEntryPtr  last;
  SeqEntryPtr  list;
  Boolean      makeMRNA;
  MolInfoPtr   mip;
  MonitorPtr   mon;
  SeqEntryPtr  nextsep;
  BioseqPtr    nucbsp;
  SeqEntryPtr  nucsep;
  ObjectIdPtr  oid;
  Boolean      parseSeqId;
  Char         path [PATH_MAX];
  ValNodePtr   rawBsps = NULL;
  ValNodePtr   rawvnp;
  SeqEntryPtr  sep;
  SeqIdPtr     sip;
  SeqLocPtr    slp;
  Char         str [64];
  BioseqPtr    target = NULL;
  Char         tmp [128];
  ValNodePtr   vnp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  extension = GetAppProperty ("FastaProtExtension");
  if (! GetInputFileName (path, sizeof (path), extension, "TEXT")) return;
  fp = FileOpen (path, "r");
  if (fp == NULL) return;

  ans = Message (MSG_YN, "Do FASTA definition lines start with seqID?");
  parseSeqId = (Boolean) (ans == ANS_YES);

  WatchCursor ();
  Update ();

  count = 0;
  list = NULL;
  last = NULL;
  head = NULL;
  errormsg = NULL;
  code = SeqEntryToGeneticCode (sep, NULL, NULL, 0);

  nucsep = FindNucSeqEntry (sep);
  slp = CreateWholeInterval (sep);

  nextsep = SequinFastaToSeqEntryEx (fp, FALSE, &errormsg, parseSeqId, NULL);
  while (nextsep != NULL) {
    count++;
    if (IS_Bioseq (nextsep) && nextsep->data.ptrvalue != NULL) {
      bsp = (BioseqPtr) nextsep->data.ptrvalue;
      isLocalUnknownID = FALSE;
      sip = bsp->id;
      if (sip != NULL && sip->choice == SEQID_LOCAL) {
        oid = (ObjectIdPtr) sip->data.ptrvalue;
        if (oid != NULL && oid->str != NULL) {
          isLocalUnknownID = (Boolean) (StringICmp (oid->str, "unknown") == 0);
        }
      }
      if ((! parseSeqId) || isLocalUnknownID) {
        sip = MakeNewProteinSeqId (slp, NULL);
        if (sip != NULL) {
          bsp->id = SeqIdFree (bsp->id);
          bsp->id = sip;
          SeqMgrReplaceInBioseqIndex (bsp);
        }
      }
    }
    SeqEntryPack (nextsep);
    if (last != NULL) {
      last->next = nextsep;
      last = nextsep;
    } else {
      list = nextsep;
      last = list;
    }
    if (! StringHasNoText (errormsg)) {
      vnp = ValNodeNew (head);
      if (head == NULL) {
        head = vnp;
      }
      if (vnp != NULL) {
        vnp->data.ptrvalue = errormsg;
        errormsg = NULL;
      }
    }
    nextsep = SequinFastaToSeqEntryEx (fp, FALSE, &errormsg, parseSeqId, NULL);
  }

  SeqLocFree (slp);
  FileClose (fp);

  ArrowCursor ();
  Update ();

  if (head != NULL) {
    for (vnp = head; vnp != NULL; vnp = vnp->next) {
      Message (MSG_POSTERR, "%s\n", (CharPtr) vnp->data.ptrvalue);
    }
    ValNodeFreeData (head);
    ans = Message (MSG_YN, "Errors detected - do you wish to proceed?");
    if (ans == ANS_NO) {
      sep = list;
      while (sep != NULL) {
        nextsep = sep->next;
        sep->next = NULL;
        SeqEntryFree (sep);
        sep = nextsep;
      }
      return;
    }
  }

  if (list == NULL) return;
  
  ans = Message (MSG_YN, "Do you wish to make default mRNAs?");
  makeMRNA = (Boolean) (ans == ANS_YES);

  WatchCursor ();
  Update ();

  nucbsp = NULL;
  nucsep = FindNucSeqEntry (sep);
  if (nucsep != NULL && IS_Bioseq (nucsep)) {
    nucbsp = (BioseqPtr) nucsep->data.ptrvalue;
  }
  if (nucbsp != NULL) {
    SetBatchSuggestNucleotide (nucbsp, code);
  }
  descr = ExtractBioSourceAndPubs (sep);
  mon = MonitorStrNewEx ("Predicting Coding Region", 20, FALSE);

  rawBsps = NULL;
  if (! doSuggest) {
    VisitBioseqsInSep (sep, (Pointer) &rawBsps, GetRawBsps);
  }
  rawvnp = rawBsps;

  count = 0;
  while (list != NULL) {
    nextsep = list->next;
    list->next = NULL;
    count++;
    if (mon != NULL) {
      str [0] = '\0';
      tmp [0] = '\0';
      bsp = (BioseqPtr) list->data.ptrvalue;
      if (bsp != NULL) {
        sip = SeqIdFindWorst (bsp->id);
        SeqIdWrite (sip, tmp, PRINTID_REPORT, sizeof (tmp));
      }
      sprintf (str, "Processing sequence %d [%s]", (int) count, tmp);
      MonitorStrValue (mon, str);
      Update ();
    }
    mip = MolInfoNew ();
    if (mip != NULL) {
      mip->biomol = 8;
      mip->tech = 13;
      vnp = CreateNewDescriptor (list, Seq_descr_molinfo);
      if (vnp != NULL) {
        vnp->data.ptrvalue = (Pointer) mip;
      }
    }
    addhere = sep;
    if (IS_Bioseq_set (addhere)) {
      bssp = (BioseqSetPtr) addhere->data.ptrvalue;
      if (bssp != NULL && bssp->_class == BioseqseqSet_class_genbank) {
        addhere = bssp->seq_set;
      }
    }
    target = NULL;
    if (! doSuggest) {
      if (rawvnp != NULL) {
        target = (BioseqPtr) rawvnp->data.ptrvalue;
        if (SeqMgrGetParentOfPart (target, NULL) == NULL) {
          addhere = SeqMgrGetSeqEntryForData (target);
        }
        rawvnp = rawvnp->next;
      }
    }
    AddSeqEntryToSeqEntry (addhere, list, TRUE);
    AutomaticProteinProcess (addhere, list, code, makeMRNA, target);
    list = nextsep;
  }

  ValNodeFree (rawBsps);

  mon = MonitorFree (mon);
  if (nucbsp != NULL) {
    ClearBatchSuggestNucleotide ();
  }
  ReplaceBioSourceAndPubs (sep, descr);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update ();
}

extern void ParseInMoreProteins (IteM i)

{
  ParseInMoreProteinsCommon (i, TRUE);
}

extern void ParseInProteinsInOrder (IteM i)

{
  ParseInMoreProteinsCommon (i, FALSE);
}

extern void ParseInMoreMRNAs (IteM i)

{
  MsgAnswer    ans;
  BaseFormPtr  bfp;
  BioseqPtr    bsp;
  Int4         count;
  CharPtr      errormsg;
  CharPtr      extension;
  FILE         *fp;
  ValNodePtr   head;
  Boolean      isLocalUnknownID;
  SeqEntryPtr  last;
  SeqEntryPtr  list;
  MonitorPtr   mon;
  SeqEntryPtr  nextsep;
  SeqEntryPtr  nucsep;
  ObjectIdPtr  oid;
  Boolean      parseSeqId;
  Char         path [PATH_MAX];
  SeqEntryPtr  sep;
  SeqIdPtr     sip;
  SeqLocPtr    slp;
  Char         str [32];
  Char         tmp [128];
  ValNodePtr   vnp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  extension = GetAppProperty ("FastaNucExtension");
  if (! GetInputFileName (path, sizeof (path), extension, "TEXT")) return;
  fp = FileOpen (path, "r");
  if (fp == NULL) return;

  ans = Message (MSG_YN, "Do FASTA definition lines start with seqID?");
  parseSeqId = (Boolean) (ans == ANS_YES);

  WatchCursor ();
  Update ();

  count = 0;
  list = NULL;
  last = NULL;
  head = NULL;
  errormsg = NULL;

  nucsep = FindNucSeqEntry (sep);
  slp = CreateWholeInterval (sep);

  nextsep = SequinFastaToSeqEntryEx (fp, TRUE, &errormsg, parseSeqId, NULL);
  while (nextsep != NULL) {
    count++;
    if (IS_Bioseq (nextsep) && nextsep->data.ptrvalue != NULL) {
      bsp = (BioseqPtr) nextsep->data.ptrvalue;
      isLocalUnknownID = FALSE;
      sip = bsp->id;
      if (sip != NULL && sip->choice == SEQID_LOCAL) {
        oid = (ObjectIdPtr) sip->data.ptrvalue;
        if (oid != NULL && oid->str != NULL) {
          isLocalUnknownID = (Boolean) (StringICmp (oid->str, "unknown") == 0);
        }
      }
      if ((! parseSeqId) || isLocalUnknownID) {
        sip = MakeNewProteinSeqId (slp, NULL);
        if (sip != NULL) {
          bsp->id = SeqIdFree (bsp->id);
          bsp->id = sip;
          SeqMgrReplaceInBioseqIndex (bsp);
        }
      }
    }
    SeqEntryPack (nextsep);
    if (last != NULL) {
      last->next = nextsep;
      last = nextsep;
    } else {
      list = nextsep;
      last = list;
    }
    if (! StringHasNoText (errormsg)) {
      vnp = ValNodeNew (head);
      if (head == NULL) {
        head = vnp;
      }
      if (vnp != NULL) {
        vnp->data.ptrvalue = errormsg;
        errormsg = NULL;
      }
    }
    nextsep = SequinFastaToSeqEntryEx (fp, TRUE, &errormsg, parseSeqId, NULL);
  }

  SeqLocFree (slp);
  FileClose (fp);

  ArrowCursor ();
  Update ();

  if (head != NULL) {
    for (vnp = head; vnp != NULL; vnp = vnp->next) {
      Message (MSG_POSTERR, "%s\n", (CharPtr) vnp->data.ptrvalue);
    }
    ValNodeFreeData (head);
    ans = Message (MSG_YN, "Errors detected - do you wish to proceed?");
    if (ans == ANS_NO) {
      sep = list;
      while (sep != NULL) {
        nextsep = sep->next;
        sep->next = NULL;
        SeqEntryFree (sep);
        sep = nextsep;
      }
      return;
    }
  }

  if (list == NULL) return;
  
  WatchCursor ();
  Update ();

  nucsep = FindNucSeqEntry (sep);
  if (nucsep == NULL) return;

  mon = MonitorStrNewEx ("Reading mRNA sequences", 20, FALSE);
  count = 0;
  while (list != NULL) {
    nextsep = list->next;
    list->next = NULL;
    count++;
    if (mon != NULL) {
      str [0] = '\0';
      tmp [0] = '\0';
      bsp = (BioseqPtr) list->data.ptrvalue;
      if (bsp != NULL) {
        sip = SeqIdFindWorst (bsp->id);
        SeqIdWrite (sip, tmp, PRINTID_REPORT, sizeof (tmp));
      }
      sprintf (str, "Processing sequence %d [%s]", (int) count, tmp);
      MonitorStrValue (mon, str);
      Update ();
    }
    AutomaticMrnaProcess (nucsep, list, FALSE, FALSE);
    SeqEntryFree (list);
    list = nextsep;
  }
  mon = MonitorFree (mon);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update ();
}

/*#ifdef ALLOW_DOWNLOAD*/
typedef struct fetchform {
  FORM_MESSAGE_BLOCK
  GrouP           accntype;
  TexT            accession;
  ButtoN          accept;
} FetchForm, PNTR FetchFormPtr;

static void FetchFormMessage (ForM f, Int2 mssg)

{
  FetchFormPtr  ffp;

  ffp = (FetchFormPtr) GetObjectExtra (f);
  if (ffp != NULL) {
    switch (mssg) {
      case VIB_MSG_CUT :
        StdCutTextProc (NULL);
        break;
      case VIB_MSG_COPY :
        StdCopyTextProc (NULL);
        break;
      case VIB_MSG_PASTE :
        StdPasteTextProc (NULL);
        break;
      case VIB_MSG_DELETE :
        StdDeleteTextProc (NULL);
        break;
      default :
        if (ffp->appmessage != NULL) {
          ffp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static void ExamineIdProc (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr  bsp;
  Int2       i;
  BoolPtr    idTypes;
  SeqIdPtr   sip;

  if (sep == NULL || sep->data.ptrvalue == NULL || mydata == NULL) return;
  idTypes = (BoolPtr) mydata;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sip = bsp->id;
    while (sip != NULL) {
      i = (Int2) sip->choice;
      if (i >= 0 && i < NUM_SEQID) {
        (idTypes [i])++;
      }
      sip = sip->next;
    }
  }
}

static Boolean OwnedByOtherDatabase (SeqEntryPtr sep, BoolPtr idTypes)

{
  Int2  i;

  if (sep == NULL || idTypes == NULL) return FALSE;
  for (i = 0; i < NUM_SEQID; i++) {
    idTypes [i] = FALSE;
  }
  BioseqExplore (sep, (Pointer) idTypes, ExamineIdProc);
  if (! (idTypes [SEQID_GENBANK])) return TRUE;
  if (idTypes [SEQID_EMBL] || idTypes [SEQID_DDBJ]) return TRUE;
  if (! FindNucSeqEntry (sep)) return TRUE;
  return FALSE;
}

static CharPtr cantUpdateMsg = "Sequin updating is\n\
currently only being tested for GenBank records.\n\
Please send e-mail to info@ncbi.nlm.nih.gov if you\n\
have any questions.";

static Int4 AccessionToGi (CharPtr string)
{
   /*
   CharPtr str;
   LinkSetPtr lsp;
   Int4 gi;

   str = MemNew (StringLen (string) + 10);
   sprintf (str, "\"%s\" [ACCN]", string);
   lsp = EntrezTLEvalString (str, TYP_NT, -1, NULL, NULL);
   MemFree (str);
   if (lsp == NULL) return 0;
   if (lsp->num <= 0) {
       LinkSetFree (lsp);
       return 0;
   }
   gi = lsp->uids [0];
   LinkSetFree (lsp);
   return gi;
   */
   Int4      gi;
   SeqIdPtr  sip;

   sip = SeqIdFromAccessionDotVersion (string);
   if (sip == NULL) return 0;
   gi = GetGIForSeqId (sip);
   SeqIdFree (sip);
   return gi;
}

static void LookForReplacedByCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr   bsp;
  SeqHistPtr  hist;
  BoolPtr     rsult;

  if (! IS_Bioseq (sep)) return;
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL) return;
  hist = bsp->hist;
  if (hist == NULL) return;
  if (hist->replaced_by_ids != NULL) {
    rsult = (BoolPtr) mydata;
    if (rsult == NULL) return;
    *rsult = TRUE;
  }
}

#ifdef USE_SMARTNET
extern Pointer ReadFromDirSub (CharPtr accn, Uint2Ptr datatype, Uint2Ptr entityID);
extern Pointer ReadFromSmart (CharPtr accn, Uint2Ptr datatype, Uint2Ptr entityID);
extern Pointer ReadFromTPASmart (CharPtr accn, Uint2Ptr datatype, Uint2Ptr entityID);
#endif

static void DownloadProc (ButtoN b)

{
  CharPtr       accn = NULL;
  MsgAnswer     ans;
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Pointer       dataptr = NULL;
  Uint2         datatype;
  CharPtr       dbname;
  Uint2         entityID;
  FetchFormPtr  ffp;
  Int2          handled;
  Boolean       idTypes [NUM_SEQID];
  Boolean       isReplaced = FALSE;
  SeqEntryPtr   sep;
  Char          str [32];
  Int4          uid;
  ForM          w;

  ffp = (FetchFormPtr) GetObjectExtra (b);
  if (ffp == NULL) return;
  w = ffp->form;
  Hide (w);
  WatchCursor ();
  Update ();
  GetTitle (ffp->accession, str, sizeof (str));
  TrimSpacesAroundString (str);
  if (StringHasNoText (str)) {
    Message (MSG_OK, "Please enter an accession number or gi");
    Show (w);
    Select (w);
    Select (ffp->accession);
    return;
  }
  sep = NULL;
  uid = 0;
  /*
  if (! EntrezIsInited ()) {
    if (! SequinEntrezInit ("Sequin", FALSE, NULL)) {
      Remove (w);
      Show (startupForm);
      Select (startupForm);
      ArrowCursor ();
      return;
    }
  }
  */
  if (GetValue (ffp->accntype) == 1) {
    /*
    sip = ValNodeNew (NULL);
    if (sip != NULL) {
      tsip = TextSeqIdNew ();
      if (tsip != NULL) {
        tsip->accession = StringSave (str);
        sip->choice = SEQID_GENBANK;
        sip->data.ptrvalue = (Pointer) tsip;
        uid = EntrezFindSeqId (sip);
        if (uid == 0) {
          sip->choice = SEQID_EMBL;
          uid = EntrezFindSeqId (sip);
        }
        if (uid == 0) {
          sip->choice = SEQID_DDBJ;
          uid = EntrezFindSeqId (sip);
        }
      }
    }
    SeqIdFree (sip);
    */
    uid = AccessionToGi (str);
    accn = str;
  } else {
    if (! StrToLong (str, &uid)) {
     uid = 0;
    }
  }
  if (uid > 0) {
    sep = PubSeqSynchronousQuery (uid, 0, -1);
    /* EntrezFini (); */
    if (sep == NULL) {
      ArrowCursor ();
      Message (MSG_OK, "Unable to find this record in the database.");
      Show (w);
      Select (w);
      Select (ffp->accession);
      return;
    }
    if (IS_Bioseq (sep)) {
      datatype = OBJ_BIOSEQ;
    } else if (IS_Bioseq_set (sep)) {
      datatype = OBJ_BIOSEQSET;
    } else {
      ArrowCursor ();
      Message (MSG_OK, "Unable to find this record in the database.");
      Show (w);
      Select (w);
      Select (ffp->accession);
      return;
    }
    Remove (w);
    SeqEntryExplore (sep, (Pointer) (&isReplaced), LookForReplacedByCallback);
    if (isReplaced) {
      ans = Message (MSG_YN, "This record has been replaced.  Are you sure you want to edit it?");
      if (ans == ANS_NO) {
        SeqEntryFree (sep);
        Show (startupForm);
        Select (startupForm);
        ArrowCursor ();
        return;
      }
    }
    dataptr = (Pointer) sep->data.ptrvalue;
  } else if (! StringHasNoText (accn)) {
#ifdef USE_SMARTNET
    if (accn != NULL) {
      dataptr = ReadFromTPASmart (accn, &datatype, NULL);
      if (dataptr == NULL) {
        dataptr = ReadFromSmart (accn, &datatype, NULL);
        if (dataptr == NULL) {
          dataptr = ReadFromDirSub (accn, &datatype, NULL);
        }
      }
    }
#endif
  }
  if (dataptr != NULL) {
    entityID = ObjMgrRegister (datatype, dataptr);
    if (dataptr != NULL && entityID > 0) {
      if (datatype == OBJ_SEQSUB || datatype == OBJ_SEQENTRY ||
          datatype == OBJ_BIOSEQ || datatype == OBJ_BIOSEQSET) {
        WatchCursor ();
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
        if (sep != NULL && OwnedByOtherDatabase (sep, idTypes)) {
          dbname = NULL;
          if (idTypes [SEQID_EMBL]) {
            dbname = "EMBL";
          } else if (idTypes [SEQID_DDBJ]) {
            dbname = "DDBJ";
          }
          /*
          if (dbname != NULL) {
            Message (MSG_OK, "This record is owned by %s.  %s", dbname, cantUpdateMsg);
          } else {
            Message (MSG_OK, "This record is not owned by GenBank.  %s", cantUpdateMsg);
          }
          */
        }
        if (datatype != OBJ_SEQSUB && uid > 0) {
          ArrowCursor ();
          Update ();
          if (Message (MSG_YN, repackageMsg) == ANS_YES) {
            globalEntityID = entityID;
            globalsep = sep;
            StringNCpy_0 (globalPath, str, sizeof (globalPath));
            WatchCursor ();
            Update ();
            w = CreateSubmitBlockForm (-50, -33, "Submitting Authors",
                                       FALSE, TRUE, NULL, JustRegisterSeqEntryBtn,
                                       AddSubmitBlockToSeqEntry);
            ArrowCursor ();
            if (w != NULL) {
              Show (w);
              Select (w);
              SendHelpScrollMessage (helpForm, "Submitting Authors Form", NULL);
              return;
            } else {
              Message (MSG_FATAL, "Unable to create window.");
              SeqEntryFree (sep);
              Show (startupForm);
              Select (startupForm);
              return;
            }
          }
        }
        seqviewprocs.filepath = str;
        seqviewprocs.forceSeparateViewer = TRUE;
        handled = GatherProcLaunch (OMPROC_VIEW, FALSE, entityID, 1,
                                    OBJ_BIOSEQ, 0, 0, OBJ_BIOSEQ, 0);
        seqviewprocs.filepath = NULL;
        ArrowCursor ();
        if (handled != OM_MSG_RET_DONE || handled == OM_MSG_RET_NOPROC) {
          Message (MSG_FATAL, "Unable to launch viewer.");
          SeqEntryFree (sep);
          Show (startupForm);
          Select (startupForm);
          return;
        } else {
          SendHelpScrollMessage (helpForm, "Editing the Record", NULL);
        }
        ObjMgrSetOptions (OM_OPT_FREE_IF_NO_VIEW, entityID);
        ObjMgrSetDirtyFlag (entityID, TRUE);
      } else {
        Message (MSG_ERROR, "Unable to process object type %d.", (int) datatype);
        ObjMgrDelete (datatype, dataptr);
        Show (startupForm);
        Select (startupForm);
        ArrowCursor ();
      }
    } else {
      Show (startupForm);
      Select (startupForm);
      ArrowCursor ();
    }
  } else {
    /* EntrezFini (); */
    ArrowCursor ();
    Message (MSG_OK, "Unable to find this record in the database.");
    Show (w);
    Select (w);
    Select (ffp->accession);
  }
}

extern void DownloadAndUpdateProc (ButtoN b)

{
  FetchFormPtr  ffp;
  SeqEntryPtr   sep;
  Char          str [32];
  Int4          uid;

  ffp = (FetchFormPtr) GetObjectExtra (b);
  if (ffp == NULL) return;
  Hide (ParentWindow (b));
  WatchCursor ();
  Update ();
  GetTitle (ffp->accession, str, sizeof (str));
  if (StringHasNoText (str)) {
    Remove (ParentWindow (b));
    ArrowCursor ();
    return;
  }
  sep = NULL;
  uid = 0;
  /*
  if (! EntrezIsInited ()) {
    if (! SequinEntrezInit ("Sequin", FALSE, NULL)) {
      Remove (ParentWindow (b));
      ArrowCursor ();
      return;
    }
  }
  */
  if (GetValue (ffp->accntype) == 1) {
    uid = AccessionToGi (str);
  } else {
    if (! StrToLong (str, &uid)) {
     uid = 0;
    }
  }
  Remove (ParentWindow (b));
  ArrowCursor ();
  Update ();
  if (uid > 0) {
    sep = PubSeqSynchronousQuery (uid, 0, -1);
    /* EntrezFini (); */
    if (sep == NULL) {
      Message (MSG_OK, "Unable to find this record in the database.");
      return;
    }
    if (IS_Bioseq (sep)) {
    } else if (IS_Bioseq_set (sep)) {
    } else {
      Message (MSG_OK, "Unable to find this record in the database.");
      return;
    }
  }

  SqnReadAlignView ((BaseFormPtr) ffp, updateTargetBspKludge, sep, TRUE);
}

extern void DownloadAndExtendProc (ButtoN b)

{
  FetchFormPtr  ffp;
  SeqEntryPtr   sep;
  Char          str [32];
  Int4          uid;

  ffp = (FetchFormPtr) GetObjectExtra (b);
  if (ffp == NULL) return;
  Hide (ParentWindow (b));
  WatchCursor ();
  Update ();
  GetTitle (ffp->accession, str, sizeof (str));
  if (StringHasNoText (str)) {
    Remove (ParentWindow (b));
    ArrowCursor ();
    return;
  }
  sep = NULL;
  uid = 0;
  /*
  if (! EntrezIsInited ()) {
    if (! SequinEntrezInit ("Sequin", FALSE, NULL)) {
      Remove (ParentWindow (b));
      ArrowCursor ();
      return;
    }
  }
  */
  if (GetValue (ffp->accntype) == 1) {
    uid = AccessionToGi (str);
  } else {
    if (! StrToLong (str, &uid)) {
     uid = 0;
    }
  }
  Remove (ParentWindow (b));
  ArrowCursor ();
  Update ();
  if (uid > 0) {
    sep = PubSeqSynchronousQuery (uid, 0, -1);
    /* EntrezFini (); */
    if (sep == NULL) {
      Message (MSG_OK, "Unable to find this record in the database.");
      return;
    }
    if (IS_Bioseq (sep)) {
    } else if (IS_Bioseq_set (sep)) {
    } else {
      Message (MSG_OK, "Unable to find this record in the database.");
      return;
    }
  }

  SqnReadAlignView ((BaseFormPtr) ffp, updateTargetBspKludge, sep, FALSE);
}

static void CancelFetchProc (ButtoN b)

{
  StdCancelButtonProc (b);
  Show (startupForm);
  Select (startupForm);
}

static void FetchTextProc (TexT t)

{
  Boolean       alldigits;
  Char          ch;
  FetchFormPtr  ffp;
  CharPtr       ptr;
  Char          str [32];

  ffp = (FetchFormPtr) GetObjectExtra (t);
  if (ffp == NULL) return;
  GetTitle (t, str, sizeof (str));
  if (StringHasNoText (str)) {
    SafeDisable (ffp->accept);
  } else {
    SafeEnable (ffp->accept);
    TrimSpacesAroundString (str);
    alldigits = TRUE;
    ptr = str;
    ch = *ptr;
    while (ch != '\0') {
      if (! IS_DIGIT (ch)) {
        alldigits = FALSE;
      }
      ptr++;
      ch = *ptr;
    }
    if (alldigits) {
      SafeSetValue (ffp->accntype, 2);
    } else {
      SafeSetValue (ffp->accntype, 1);
    }
  }
}

extern void CommonFetchFromNet (BtnActnProc actn, BtnActnProc cancel)

{
  GrouP              c;
  FetchFormPtr       ffp;
  GrouP              g;
  StdEditorProcsPtr  sepp;
  WindoW             w;

  Hide (startupForm);
  Update ();
  w = NULL;
  ffp = MemNew (sizeof (FetchForm));
  if (ffp != NULL) {
    w = FixedWindow (-50, -33, -10, -10, "Download From Entrez", NULL);
    SetObjectExtra (w, ffp, StdCleanupFormProc);
    ffp->form = (ForM) w;
    ffp->formmessage = FetchFormMessage;

#ifndef WIN_MAC
    CreateSqnInitialFormMenus (w);
#endif

    sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
    if (sepp != NULL) {
      ffp->appmessage = sepp->handleMessages;
    }
    SetGroupSpacing (w, 10, 10);

    g = HiddenGroup (w, -3, 0, NULL);
    StaticPrompt (g, "Type", 0, stdLineHeight, programFont, 'l');
    ffp->accntype = HiddenGroup (g, 4, 0, NULL);
    RadioButton (ffp->accntype, "Accession");
    RadioButton (ffp->accntype, "GI");
    SetValue (ffp->accntype, 1);
    ffp->accession = DialogText (g, "", 6, FetchTextProc);
    SetObjectExtra (ffp->accession, ffp, NULL);

    c = HiddenGroup (w, 4, 0, NULL);
    SetGroupSpacing (c, 10, 2);
    ffp->accept = DefaultButton (c, "Retrieve", actn);
    SetObjectExtra (ffp->accept, ffp, NULL);
    Disable (ffp->accept);
    PushButton (c, "Cancel", cancel);

    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
    RealizeWindow (w);

    if (sepp != NULL) {
      SetActivate (w, sepp->activateForm);
    }
    Select (ffp->accession);
    Show (w);
    Select (w);
    Update ();
  } else {
    Show (startupForm);
    Select (startupForm);
  }
}

extern void FetchFromNet (ButtoN b)

{
  CommonFetchFromNet (DownloadProc, CancelFetchProc);
}

/*#else
#define FetchFromNet NULL
#endif*/

/*
static Boolean FindPerfectSubMatch (CharPtr prot, CharPtr trans, Int4 start,
                                    Int4 len, Uint1 frame, Int2 strand,
                                    Int4Ptr fromPtr, Int4Ptr toPtr)

{
  int      ch;
  Int2     d [256];
  Int4     from;
  int      i;
  int      j;
  int      k;
  size_t   protLen;
  Boolean  rsult;
  Int4     to;
  size_t   transLen;

  rsult = FALSE;
  from = 0;
  to = 0;
  if (prot != NULL && trans != NULL) {
    protLen = StringLen (prot);
    transLen = StringLen (trans);
    if (protLen <= transLen) {
      for (ch = 0; ch < 256; ch++) {
        d [ch] = protLen;
      }
      for (j = 0; j < protLen - 1; j++) {
        d [(int) prot [j]] = protLen - j - 1;
      }
      i = protLen;
      do {
        j = protLen;
        k = i;
        do {
          k--;
          j--;
        } while (j >= 0 && prot [j] == trans [k]);
        if (j >= 0) {
          i += d [(int) trans [i - 1]];
        }
      } while (j >= 0 && i <= transLen);
      if (j < 0) {
        i -= protLen;
        from = (long) (i * 3 + (frame - 1));
        to = from + 3 * protLen;
        if (trans [i + protLen] == '*') {
          to += 3;
        }
        if (strand == Seq_strand_plus) {
          from += 1;
        } else if (strand == Seq_strand_minus) {
          from = len - from;
          to = len - to + 1;
        }
        rsult = TRUE;
      }
    }
  }
  if (fromPtr != NULL) {
    *fromPtr = from + start;
  }
  if (toPtr != NULL) {
    *toPtr = to + start;
  }
  return rsult;
}

static Boolean CheckOneFrame (BioseqPtr bsp, Int4 start, Int4 len,
                              CharPtr prot, Int2 gencode,
                              Uint1 frame, Int2 strand,
                              Int4Ptr fromPtr, Int4Ptr toPtr)

{
  ByteStorePtr  bs;
  Char          ch;
  ValNodePtr    code;
  CdRegionPtr   crp;
  CharPtr       ptr;
  Boolean       rsult;
  SeqFeatPtr    sfp;
  CharPtr       trans;
  ValNodePtr    vnp;

  rsult = FALSE;
  if (bsp != NULL && gencode > 0) {
    sfp = SeqFeatNew ();
    if (sfp != NULL) {
      sfp->data.choice = SEQFEAT_CDREGION;
      crp = CdRegionNew ();
      sfp->data.value.ptrvalue = (Pointer) crp;
      if (crp != NULL) {
        crp->orf = FALSE;
        crp->conflict = FALSE;
        crp->frame = frame;
        crp->gaps = 0;
        crp->mismatch = 0;
        crp->stops = 0;
        code = ValNodeNew (NULL);
        if (code != NULL) {
          code->choice = 254;
          vnp = ValNodeNew (NULL);
          code->data.ptrvalue = vnp;
          if (vnp != NULL) {
            vnp->choice = 2;
            vnp->data.intvalue = (Int4) gencode;
          }
        }
        crp->genetic_code = code;
        crp->code_break = NULL;
        AddIntToSeqFeat (sfp, start, start + len - 1, bsp, -1, -1, strand);
        trans = NULL;
        bs = ProteinFromCdRegion (sfp, TRUE);
        if (bs != NULL) {
          trans = BSMerge (bs, NULL);
          BSFree (bs);
        }
        if (trans != NULL) {
          ptr = trans;
          ch = *ptr;
          while (ch != '\0') {
            *ptr = TO_UPPER (ch);
            ptr++;
            ch = *ptr;
          }
          if (trans [0] == '-') {
            trans [0] = prot [0];
          }
          rsult = FindPerfectSubMatch (prot, trans, start, len,
                                       frame, strand, fromPtr, toPtr);
          MemFree (trans);
        }
      }
      SeqFeatFree (sfp);
    }
  }
  return rsult;
}

#define PREDICT_BLOCK_SIZE 30000L

static SeqLocPtr FindSingleCodingInterval (BioseqPtr nuc, BioseqPtr prot, Int2 genCode)

{
  Int4        cdsFrom;
  Int4        cdsTo;
  Char        ch;
  Int4        cntr;
  Uint1       frame;
  Int4        from;
  Int4        incr;
  Int4        len;
  Boolean     matched;
  size_t      protLen;
  CharPtr     protstr;
  CharPtr     ptr;
  SeqFeatPtr  sfp;
  SeqLocPtr   slp;
  Int4        start;
  Int2        strand;
  Int4        tmp;
  Int4        to;

  slp = NULL;
  if (nuc != NULL && prot != NULL) {
    cdsFrom = 0;
    cdsTo = 0;
    strand = Seq_strand_unknown;
    protstr = NULL;
    if (prot->length > 0) {
      protstr = BSMerge (prot->seq_data, NULL);
      if (protstr != NULL) {
        ptr = protstr;
        ch = *ptr;
        while (ch != '\0') {
          *ptr = TO_UPPER (ch);
          ptr++;
          ch = *ptr;
        }
        protLen = StringLen (protstr);
        matched = FALSE;
        for (frame = 1; frame <= 3 && (! matched); frame++) {
          strand = Seq_strand_plus;
          start = 0;
          cntr = nuc->length;
          len = MIN (cntr, (Int4) (PREDICT_BLOCK_SIZE + (Int4) protLen * 3L));
          while (len > 0 && (! matched)) {
            incr = MIN (cntr, PREDICT_BLOCK_SIZE);
            matched = CheckOneFrame (nuc, start, len, protstr, genCode, frame,
                                     strand, &cdsFrom, &cdsTo);
            start += incr;
            cntr -= incr;
            len = MIN (cntr, (Int4) (PREDICT_BLOCK_SIZE + (Int4) protLen * 3L));
          }
        }
        for (frame = 1; frame <= 3 && (! matched); frame++) {
          strand = Seq_strand_minus;
          start = 0;
          cntr = nuc->length;
          len = MIN (cntr, (Int4) (PREDICT_BLOCK_SIZE + (Int4) protLen * 3L));
          while (len > 0 && (! matched)) {
            incr = MIN (cntr, PREDICT_BLOCK_SIZE);
            matched = CheckOneFrame (nuc, start, len, protstr, genCode, frame,
                                     strand, &cdsFrom, &cdsTo);
            start += incr;
            cntr -= incr;
            len = MIN (cntr, (Int4) (PREDICT_BLOCK_SIZE + (Int4) protLen * 3L));
          }
        }
        if (matched) {
          sfp = SeqFeatNew ();
          if (sfp != NULL) {
            from = cdsFrom - 1;
            to = cdsTo - 1;
            if (from > to) {
              tmp = from;
              from = to;
              to = tmp;
            }
            AddIntToSeqFeat (sfp, from, to, nuc, -1, -1, strand);
            slp = sfp->location;
            sfp->location = NULL;
          }
          SeqFeatFree (sfp);
        }
      }
      MemFree (protstr);
    }
  }
  return slp;
}
*/

static Boolean ReplaceTPAAccessionNumbers (
  CharPtr    seqid,
  ValNodePtr acc_list,
  SeqEntryPtr sep
)
{
  BioseqSetPtr      bssp;
  BioseqPtr         bsp;
  SeqDescrPtr       sdp;
  Char              str [128];
  SeqMgrDescContext context;
  UserObjectPtr     uop;
  ValNodePtr        vnp;

  if (sep == NULL || sep->data.ptrvalue == NULL) return FALSE;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    /* this also delves into nuc-prot sets */
    if (bssp != NULL && (bssp->_class == 7 ||
                         (IsPopPhyEtcSet (bssp->_class)) ||
                         bssp->_class == 1)) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next)
      {
        if (ReplaceTPAAccessionNumbers (seqid, acc_list, sep))
        {
          return TRUE;
        }
      }
      return FALSE;
    }
  }
  if (!IS_Bioseq (sep)) return FALSE;

  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL) return FALSE;
  SeqIdWrite (bsp->id, str, PRINTID_REPORT, sizeof (str));
  if (StringCmp (str, seqid) != 0) return FALSE;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &context);
  while (sdp != NULL
    && ((uop = (UserObjectPtr)sdp->data.ptrvalue) == NULL
      || StringICmp (uop->type->str, "TpaAssembly") != 0))
  {
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &context);
  }
  if (sdp == NULL)
  {
    sdp = CreateNewDescriptor (sep, Seq_descr_user);
    if (sdp == NULL) return FALSE;
    uop = CreateTpaAssemblyUserObject ();
    if (uop == NULL) return FALSE;
    sdp->data.ptrvalue = uop;
  }

  for (vnp = acc_list; vnp != NULL; vnp = vnp->next)
  {
    AddAccessionToTpaAssemblyUserObject (uop, vnp->data.ptrvalue, 0, 0);
  }
  ValNodeFreeData (acc_list);
  
  return TRUE;
}

extern void LoadTPAAccessionNumbersFromFile (
  IteM i
)
{
  BaseFormPtr   bfp;
  SeqEntryPtr   sep;
  Char          path [PATH_MAX];
  FILE          *fp;
  Char          str [8192];
  size_t        len = 8192;
  Boolean       need_seqid;
  Char          seqid[100];
  Int4          seqid_len;
  CharPtr       cp;
  CharPtr       acc_end;
  Boolean       found_end;
  ValNodePtr    acc_list;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  path [0] = '\0';
  if (! GetInputFileName (path, sizeof (path), NULL, "TEXT")) return;
  
  fp = FileOpen (path, "r");
  if (fp == NULL) return;

  need_seqid = TRUE;
  acc_list = NULL;
  ReadLine (fp, str, len);
  while (Nlm_fileDone) 
  {
    cp = str;
    if (strlen (str) == 0)
    {
      ReadLine (fp, str, len);
      continue;
    }
    if (need_seqid)
    {
      seqid_len = StringCSpn (str, " \t");
      if (seqid_len > 0)
      {
        StringNCpy (seqid, str, seqid_len);
        seqid [seqid_len] = 0;
        need_seqid = FALSE;
      }
      cp = str + seqid_len + 1;
    }
    if (need_seqid)
    {
      ReadLine (fp, str, len);
      continue;
    }
    if (str [strlen (str) - 1] != ',')
    {
      need_seqid = TRUE;
    }
    
    found_end = FALSE;
    while (*cp != 0)
    {
      if (*cp == ' ' || *cp == ',' || *cp == '\t')
      {
        cp++;
      }
      else
      {
        acc_end = cp + 1;
        while (*acc_end != 0 && *acc_end != ',')
        {
          acc_end++;
        }
        if (*acc_end == 0)
        {
          found_end = TRUE;
        }
        else
        {
          *acc_end = 0;
        }
        ValNodeAddStr (&acc_list, 0, StringSave (cp));
        if (found_end)
        {
          cp = acc_end;
        }
        else
        {
          cp = acc_end + 1;
        }
      }
    }

    if (need_seqid == TRUE)
    {
      /* do something with accession list */
      if ( ! ReplaceTPAAccessionNumbers (seqid, acc_list, sep))
      {
        Message (MSG_ERROR,
                 "Unable to update accession numbers for %s (not found)",
                 seqid);
      }
      acc_list = NULL;
    }
      
    ReadLine (fp, str, len);
  }
  if (acc_list != NULL
    && ! ReplaceTPAAccessionNumbers (seqid, acc_list, sep))
  {
    Message (MSG_ERROR,
             "Unable to update accession numbers for %s (not found)",
             seqid);
  }

  FileClose (fp);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update ();
  return;
}

static void AddHistory (
  BioseqPtr  bsp,
  ValNodePtr acc_list
)
{
  SeqHistPtr      hist;
  ValNodePtr      vnp;
  SeqIdPtr        sip;
  TextSeqIdPtr    tsip;
  Uint4           whichdb;
  Char            prefix [20];

  if (bsp == NULL || acc_list == NULL) return;
  hist = bsp->hist;
  if (hist == NULL)
  {
    hist = SeqHistNew ();
    if (hist == NULL) return;
    bsp->hist = hist;
  }
  for (vnp = acc_list; vnp != NULL; vnp = vnp->next) {
    tsip = TextSeqIdNew ();
    if (tsip == NULL) return;
    tsip->accession = StringSave (vnp->data.ptrvalue);

    sip = ValNodeNew (hist->replace_ids);
    if (hist->replace_ids == NULL) {
      hist->replace_ids = sip;
    }
    if (sip == NULL) return;

    sip->data.ptrvalue = (Pointer) tsip;

    StringNCpy_0 (prefix, (CharPtr) vnp->data.ptrvalue, sizeof (prefix));
    whichdb = WHICH_db_accession (prefix);
    if (ACCN_IS_EMBL (whichdb)) {
      sip->choice = SEQID_EMBL;
    } else if (ACCN_IS_DDBJ (whichdb)) {
      sip->choice = SEQID_DDBJ;
    } else {
      sip->choice = SEQID_GENBANK;
    }
  }
  if (hist != NULL
    && hist->assembly == NULL
    && hist->replace_date == NULL
    && hist->replace_ids == NULL
    && hist->replaced_by_date == NULL
    && hist->replaced_by_ids == NULL
    && hist->deleted_date == NULL
    && ! hist->deleted)
  {
      bsp->hist = SeqHistFree (bsp->hist);
  }
}

static Boolean DoIDsMatch (CharPtr seqid, BioseqPtr bsp, Boolean AllowLocal)
{
  Char         str [128];
  Int4         seqid_len;
  SeqIdPtr     sip;

  if (bsp == NULL) return FALSE;

  SeqIdWrite (bsp->id, str, PRINTID_REPORT, sizeof (str));
  seqid_len = StringCSpn (str, ".");
  if (seqid_len > 0)
  {
    str [ seqid_len ] = 0;
  }
  if (StringCmp (str, seqid) == 0) return TRUE;

  for (sip = bsp->id; sip != NULL; sip = sip->next)
  {
    if (sip->choice == SEQID_LOCAL && AllowLocal)
    {
      SeqIdWrite (bsp->id, str, PRINTID_REPORT, sizeof (str));
      seqid_len = StringCSpn (str, ".");
      if (seqid_len > 0)
      {
        str [ seqid_len ] = 0;
      }
      if (StringCmp (str, seqid) == 0) return TRUE;
    }
  }
  return FALSE;
}

static Boolean AddAccessionToGenbankBlock (
  CharPtr     seqid,
  ValNodePtr  acc_list,
  SeqEntryPtr sep,
  Boolean     add_hist
)
{
  BioseqPtr    bsp;
  BioseqSetPtr bssp;
  GBBlockPtr   gbp;
  ValNodePtr   last_one;
  SeqDescrPtr       sdp;

  if (seqid == NULL || acc_list == NULL
    || sep == NULL || sep->data.ptrvalue == NULL) return FALSE;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    /* this also delves into nuc-prot sets */
    if (bssp != NULL && (bssp->_class == 7 ||
                         (IsPopPhyEtcSet (bssp->_class)) ||
                         bssp->_class == 1)) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next)
      {
        if (AddAccessionToGenbankBlock (seqid, acc_list, sep, add_hist))
        {
          return TRUE;
        }
      }
      return FALSE;
    }
  }
  if (!IS_Bioseq (sep)) return FALSE;

  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL) return FALSE;
  if (! DoIDsMatch (seqid, bsp, TRUE)) return FALSE;

  sdp = BioseqGetSeqDescr (bsp, Seq_descr_genbank, NULL);

  if (sdp == NULL)
  {
    sdp = CreateNewDescriptor (sep, Seq_descr_genbank);
    if (sdp == NULL) return FALSE;
  }
 
  if (sdp->data.ptrvalue == NULL)
  {
    sdp->data.ptrvalue = GBBlockNew ();
    if (sdp->data.ptrvalue == NULL) return FALSE;
  }
 
  gbp = (GBBlockPtr) sdp->data.ptrvalue;
  
  for (last_one = gbp->extra_accessions;
       last_one != NULL && last_one->next != NULL;
       last_one = last_one->next)
  {}
  if (last_one == NULL)
  {
    gbp->extra_accessions = acc_list;
  }
  else
  {
    last_one->next = acc_list;
  }
  if (add_hist)
  {
    AddHistory (bsp, acc_list);
  }
  return TRUE;
}
 
static void LoadSecondaryAccessionNumbersPlusHistFromFile (
  IteM    i,
  Boolean add_hist
)
{
  BaseFormPtr   bfp;
  SeqEntryPtr   sep;
  Char          path [PATH_MAX];
  FILE          *fp;
  Char          str [8192];
  size_t        len = 8192;
  Boolean       need_seqid;
  Char          seqid[100];
  Int4          seqid_len;
  CharPtr       cp;
  CharPtr       acc_end;
  Boolean       found_end;
  ValNodePtr    acc_list;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  path [0] = '\0';
  if (! GetInputFileName (path, sizeof (path), NULL, "TEXT")) return;
  
  fp = FileOpen (path, "r");
  if (fp == NULL) return;

  need_seqid = TRUE;
  acc_list = NULL;
  ReadLine (fp, str, len);
  while (Nlm_fileDone || str[0] != 0) 
  {
    cp = str;
    if (strlen (str) == 0)
    {
      ReadLine (fp, str, len);
      continue;
    }
    seqid_len = StringCSpn (str, " \t");
    if (seqid_len > 0)
    {
      StringNCpy (seqid, str, seqid_len);
      seqid [seqid_len] = 0;
      cp = str + seqid_len + 1;
    }
    else
    {
      ReadLine (fp, str, len);
      continue;
    }
    
    found_end = FALSE;
    while (*cp != 0)
    {
      if (*cp == ' ' || *cp == ' ')
      {
        cp++;
      }
      else
      {
        acc_end = cp + 1;
        while (*acc_end != 0 && *acc_end != ' ')
        {
          acc_end++;
        }
        if (*acc_end == 0)
        {
          found_end = TRUE;
        }
        else
        {
          *acc_end = 0;
        }
        ValNodeAddStr (&acc_list, 0, StringSave (cp));
        if (found_end)
        {
          cp = acc_end;
        }
        else
        {
          cp = acc_end + 1;
        }
      }
    }

    /* do something with accession list */
    if ( ! AddAccessionToGenbankBlock (seqid, acc_list, sep, add_hist))
    {
      Message (MSG_ERROR,
               "Unable to update accession numbers for %s (not found)",
               seqid);
    }
    acc_list = NULL;
      
    ReadLine (fp, str, len);
  }
  if (acc_list != NULL
    && ! AddAccessionToGenbankBlock (seqid, acc_list, sep, add_hist))
  {
    Message (MSG_ERROR,
             "Unable to update accession numbers for %s (not found)",
             seqid);
  }

  FileClose (fp);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update ();
  return;
}

extern void LoadSecondaryAccessionNumbersFromFile (
  IteM i
)
{
  LoadSecondaryAccessionNumbersPlusHistFromFile (i, FALSE);
}

extern void LoadHistoryAccessionNumbersFromFile (
  IteM i
)
{
  LoadSecondaryAccessionNumbersPlusHistFromFile (i, TRUE);
}

/* This section of code is used for managing lists of features.
 * Sometimes the features will be displayed alphabetically, sometimes
 * they will be displayed alphabetically with a list of the most used features
 * also appearing at the top of the list.
 */

/* This is used to compare feature names with the special alphabetical order */
static int CompareFeatureNames (CharPtr cp1, CharPtr cp2)
{
  /* NULL name goes at the end */
  if (cp1 == NULL && cp2 == NULL) return 0;
  if (cp1 == NULL) return 1;
  if (cp2 == NULL) return -1;

  /* starts with a space goes at the beginning */
  if (cp1 [0] == ' ' && cp2 [0] == ' ') return 0;
  if (cp1 [0] == ' ') return -1;
  if (cp2 [0] == ' ') return 1;

  /* Is "All" or [ALL FEATURES] goes at the beginning */
  if ((StringCmp (cp1, "All") == 0
    || StringCmp (cp1, "[ALL FEATURES]") == 0)
    && (StringCmp (cp2, "All") == 0
    || StringCmp (cp2, "[ALL FEATURES]") == 0))
  {
    return 0;
  }
  if (StringCmp (cp1, "All") == 0
    || StringCmp (cp1, "[ALL FEATURES]") == 0)
  {
    return -1;
  }
  if (StringCmp (cp2, "All") == 0
    || StringCmp (cp2, "[ALL FEATURES]") == 0)
  {
    return 1;
  }

  /* starts with a number -> goes at the end */
  if (cp1 [0] >= '0' && cp1 [0] <= '9'
   && cp2 [0] >= '0' && cp2 [0] <= '9')
  {
    return StringICmp (cp1, cp2);
  }
  if (cp1 [0] >= '0' && cp1 [0] <= '9')
  {
    return 1;
  }
  if (cp2 [0] >= '0' && cp2 [0] <= '9')
  {
    return -1;
  }

  /* starts with a tilde or dash - sort with other tildes, put before numbers after alphas */
  if (cp1 [0] == '~' && cp2 [0] == '~') 
  {
    return StringICmp (cp1 + 1, cp2 + 1);
  }
  if (cp1 [0] == '~') return 1;
  if (cp2 [0] == '~') return -1;

  if (cp1 [0] == '-' && cp2 [0] == '-') 
  {
    return StringICmp (cp1 + 1, cp2 + 1);
  }
  if (cp1 [0] == '-') return 1;
  if (cp2 [0] == '-') return -1;

  return StringICmp (cp1, cp2);
}

extern int LIBCALLBACK CompareFeatureValNodeStrings (VoidPtr ptr1, VoidPtr ptr2)
{
  ValNodePtr vnp1, vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;

  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);

  if (vnp1 == NULL || vnp2 == NULL) return 0;

  return CompareFeatureNames (vnp1->data.ptrvalue, vnp2->data.ptrvalue);
}

extern int LIBCALLBACK CompareImpFeatEnumFieldAssoc (VoidPtr ptr1, VoidPtr ptr2)
{
  ValNodePtr        vnp1, vnp2;
  EnumFieldAssocPtr ap1, ap2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;

  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;

  ap1 = (EnumFieldAssocPtr) vnp1->data.ptrvalue;
  ap2 = (EnumFieldAssocPtr) vnp2->data.ptrvalue;
  if (ap1 == NULL || ap2 == NULL) return 0;

  return CompareFeatureNames (ap1->name, ap2->name);
}

CharPtr MostUsedFeatureList[] = { 
  "CDS",
  "exon",
  "Gene",
  "intron",
  "mRNA",
  "rRNA",
  "RNA"
};

extern ValNodePtr InsertMostUsedFeatureValNodes (ValNodePtr old_list)
{
  ValNodePtr new_list, new_item, old_item;
  Int4       index;

  new_list = NULL;
  for (index = 0;
       index < sizeof (MostUsedFeatureList) / sizeof (CharPtr);
       index ++)
  {
    old_item = FindExactStringInStrings ( old_list, MostUsedFeatureList [index])
;
    if (old_item == NULL) continue;
    new_item = ValNodeNew ( new_list);
    if (new_item == NULL) return old_list;
    new_item->choice = old_item->choice;
    new_item->data.ptrvalue = StringSave (MostUsedFeatureList [index]);
    if (new_list == NULL) new_list = new_item;
  }
  if (new_item != NULL)
  {
    if (old_list != NULL &&
      ( StringCmp (old_list->data.ptrvalue, "All") == 0
       || StringCmp (old_list->data.ptrvalue, "[ALL FEATURES]") == 0))
    {
      new_item->next = old_list->next;
      old_list->next = new_list;
      new_list = old_list;
    }
    else
    {
      new_item->next = old_list;
    }
  }
  else
  {
    new_list = old_list;
  }
  return new_list;
}

static EnumFieldAssocPtr FindEnumFieldAssoc (
  EnumFieldAssocPtr alist,
  CharPtr findStr
)
{
  EnumFieldAssocPtr ap;
  
  for (ap = alist; ap != NULL && ap->name != NULL; ap++)
  {
    if (StringCmp (ap->name, findStr) == 0) return ap;
  }
  return NULL;
}

static void CopyEnumFieldAssoc (EnumFieldAssocPtr ap1, EnumFieldAssocPtr ap2)
{
  if (ap1 == NULL || ap2 == NULL) return;

  ap1->name = StringSave (ap2->name);
  ap1->value = ap2->value;
}

extern EnumFieldAssocPtr InsertMostUsedFeatureEnumFieldAssoc (
  EnumFieldAssocPtr alist
)
{
  Int4              num_total_fields, index, new_index;
  EnumFieldAssocPtr ap, new_alist, old_ap;

  num_total_fields = sizeof (MostUsedFeatureList) / sizeof (CharPtr);

  for (ap = alist; ap != NULL && ap->name != NULL; ap++)
  {
    num_total_fields ++;
  }
  /* need the last null field */
  num_total_fields ++;

  new_alist = MemNew (num_total_fields * sizeof (EnumFieldAssoc));
  if (new_alist == NULL) return alist;

  /* copy the first item if wildcard */
  if (StringCmp (alist->name, "[ALL FEATURES]") == 0)
  {
    CopyEnumFieldAssoc (new_alist, alist);
    new_index = 1;
  }
  else
  {
    new_index = 0;
  }

  for (index = 0;
       index < sizeof (MostUsedFeatureList) / sizeof (CharPtr);
       index ++)
  {
    old_ap = FindEnumFieldAssoc (alist, MostUsedFeatureList [index]);
    if (old_ap == NULL) continue;
    CopyEnumFieldAssoc (new_alist + new_index++, old_ap);
  }

  for (ap = alist; ap != NULL && ap->name != NULL; ap++)
  {
    CopyEnumFieldAssoc (new_alist + new_index ++, ap);
  }
  /* copy over the last null field */
  if (ap != NULL)
  {
    CopyEnumFieldAssoc (new_alist + new_index ++, ap);
  }
  return new_alist;
  
}

extern void SortEnumFieldAssocPtrArray (EnumFieldAssocPtr alist, CompareFunc compar)
{
  ValNodePtr        head, vnp;
  EnumFieldAssocPtr ap;
  Int4              index;

  /* first, create ValNode list so we can sort the data */
  head = NULL;
  for (ap = alist; ap != NULL && ap->name != NULL; ap++)
  {
    vnp = ValNodeNew (head);
    if (vnp == NULL) return;
    vnp->data.ptrvalue = MemNew (sizeof (EnumFieldAssoc));
    if (vnp->data.ptrvalue == NULL) return;
    MemCpy (vnp->data.ptrvalue, ap, sizeof (EnumFieldAssoc));
    if (head == NULL) head = vnp;
  }

  /* Now sort the ValNode list */
  head = SortValNode (head, compar);

  /* Now repopulate the EnumFieldAssoc list */
  index = 0;
  for (vnp = head; vnp != NULL; vnp = vnp->next)
  {
    MemCpy (alist + index++, vnp->data.ptrvalue, sizeof (EnumFieldAssoc));
  }

  /* And free the ValNode list */
  ValNodeFreeData (head);
}

static Uint2 UnusualFeatureTypes [] = {
  FEATDEF_ORG,
  FEATDEF_mutation,
  FEATDEF_site_ref,
  FEATDEF_gap,
  FEATDEF_NON_STD_RESIDUE,
  FEATDEF_NUM
};
 
extern ValNodePtr BuildFeatureValNodeList (
  Boolean prefer_most_used,
  CharPtr wild_card_name,
  Int4    wild_card_value,
  Boolean skip_unusual,
  Boolean skip_import
)
{
  FeatDefPtr  curr;
  ValNodePtr  head, vnp;
  Uint1       key;
  CharPtr     label = NULL;
  Uint2       subtype;
  Int4        index;
  Boolean     skip;
  Char        str [256];

  head = NULL;
  curr = FeatDefFindNext (NULL, &key, &label, FEATDEF_ANY, TRUE);
  while (curr != NULL) {
    skip = FALSE;
    if (skip_unusual)
    {
      for (index = 0;
           ! skip && index < sizeof ( UnusualFeatureTypes ) / sizeof (Uint2);
           index ++)
      {
        if (curr->featdef_key == UnusualFeatureTypes [ index ]) skip = TRUE;
      }
    }
    if (key != FEATDEF_BAD && ! skip) {
      
      subtype = curr->featdef_key;
	  if (subtype == FEATDEF_PUB)
	  {
        StringNCpy_0 (str, curr->typelabel, sizeof (str) - 15);
        StringCat (str, " (Publication)");
	  }
	  else if (subtype != FEATDEF_misc_RNA &&
          subtype != FEATDEF_precursor_RNA &&
          subtype != FEATDEF_mat_peptide &&
          subtype != FEATDEF_sig_peptide &&
          subtype != FEATDEF_transit_peptide &&
          subtype != FEATDEF_Imp_CDS)
      {
        StringNCpy_0 (str, curr->typelabel, sizeof (str) - 1);
      }
      else if (! skip_import)
      {
        StringNCpy_0 (str, curr->typelabel, sizeof (str) - 10);
        StringCat (str, "_imp");
      }
      else
      {
        skip = TRUE;
      }
      if (! skip)
      {
        vnp = ValNodeNew (head);
        if (head == NULL) {
          head = vnp;
        }
        if (vnp != NULL) {
          vnp->choice = subtype;
          vnp->data.ptrvalue = StringSave (str);
        }
      }
    }
    curr = FeatDefFindNext (curr, &key, &label, FEATDEF_ANY, TRUE);
  }
  if (head != NULL) {
    head = SortValNode (head, CompareFeatureValNodeStrings);
    head = InsertMostUsedFeatureValNodes (head);
    if (wild_card_name != NULL)
    {
      vnp = ValNodeNew (NULL);
      if (vnp != NULL) {
        vnp->choice = wild_card_value;
        vnp->data.ptrvalue = StringSave (wild_card_name);
        vnp->next = head;
        head = vnp;
      }
    }
  }
  return head;
}

extern void SetTaxNameAndRemoveTaxRef (OrgRefPtr orp, CharPtr taxname)
{
  ValNodePtr      vnp, next;
  ValNodePtr PNTR prev;
  DbtagPtr        dbt;
  Boolean         remove_taxrefs = FALSE;

  if (orp == NULL) return;

  if ( taxname == NULL || orp->taxname == NULL
    || StringCmp (taxname, orp->taxname) != 0)
  {
    remove_taxrefs = TRUE;
  }
  MemFree (orp->taxname);
  orp->taxname = taxname;

  if (! remove_taxrefs) return;

  orp->common = MemFree (orp->common);

  vnp = orp->db;
  if (vnp == NULL) return;
  prev = (ValNodePtr PNTR) &(orp->db);
  while (vnp != NULL) {
    next = vnp->next;
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt != NULL && StringICmp ((CharPtr) dbt->db, "taxon") == 0) {
      *prev = vnp->next;
      vnp->next = NULL;
      DbtagFree (dbt);
      ValNodeFree (vnp);
    } else {
      prev = (ValNodePtr PNTR) &(vnp->next);
    }
    vnp = next;
  }
}

static Boolean
FindMatchingProprotein 
(SeqFeatPtr sfp,
 SeqMgrFeatContextPtr fcontext,
 BioseqPtr prot_bsp)
{
  SeqFeatPtr        prot_sfp;
  SeqMgrFeatContext pcontext;
  CharPtr           start;

  if (prot_bsp == NULL || fcontext == NULL) return FALSE;
  if (StringNICmp (fcontext->label, "encodes ", 8) == 0) {
    start = fcontext->label + 8;
  } else {
    start = fcontext->label;
  }
  prot_sfp = NULL;
  while ((prot_sfp = SeqMgrGetNextFeature (prot_bsp, prot_sfp, 
                                           0, 0, &pcontext)) != NULL) {
    if (StringCmp (pcontext.label, start) == 0) {
      return TRUE;
    } 
  }
  return FALSE;
}


static void 
RemoveRedundantProproteinMiscFeatsOnBioseq
(BioseqPtr bsp,
 Pointer userdata)
{
  SeqFeatPtr        sfp, cds;
  SeqMgrFeatContext fcontext, cds_context;
  BioseqPtr         bsp_prot;

  sfp = NULL;

  /* list misc feats */
  while ((sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext)) != NULL) {
    if (fcontext.featdeftype == FEATDEF_misc_feature
        &&  StringStr(fcontext.label, "proprotein") != NULL) {
      cds = NULL;
      while ((cds = SeqMgrGetNextFeature (bsp, cds, SEQFEAT_CDREGION, 0, &cds_context)) != NULL) {
        if (cds_context.left <= fcontext.left
            &&  cds_context.right >= fcontext.right) {
          /* Get Protein sequence, look for matching proprotein feat */
          bsp_prot = BioseqFind (SeqLocId(cds->product));
          if (FindMatchingProprotein (sfp, &fcontext, bsp_prot)) {
            sfp->idx.deleteme = TRUE;
          }
        }
      }
    }
  }

}


extern void RemoveRedundantProproteinMiscFeats (IteM i)
{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  /* Visit each bioseq to remove redundant proprotein misc feats */
  VisitBioseqsInSep (sep, NULL, RemoveRedundantProproteinMiscFeatsOnBioseq);

  DeleteMarkedObjects (bfp->input_entityID, 0, NULL);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update ();
}

typedef struct typestraindata
{
  FORM_MESSAGE_BLOCK
  
  TexT   find_this_txt;
  ButtoN when_string_not_found_btn;
  ButtoN case_insensitive_btn;
  GrouP  string_loc_grp;
  PopuP  field_choice_popup;
  ButtoN remove_found_text_btn;
  
  Boolean when_string_not_found;
  Boolean case_insensitive;
  Int4    string_loc;
  Int4    field_choice;
  CharPtr find_this;
  Boolean remove_found_text;
} TypeStrainData, PNTR TypeStrainPtr;

static Boolean MeetsTypeStrainConstraint (BioSourcePtr biop, TypeStrainPtr tsp)
{
  CharPtr      string_found = NULL;
  CharPtr      search_text = NULL;
  OrgModPtr    mod = NULL, prev_mod = NULL;
  SubSourcePtr ssp = NULL, prev_ssp = NULL;
  Boolean      rval;
  CharPtr      cp, destp;
  
  if (biop == NULL) return FALSE;
  if (tsp == NULL) return TRUE;
  if (StringHasNoText (tsp->find_this)) return TRUE;
  if (biop->org == NULL)
  {
  	if (tsp->when_string_not_found)
  	{
  	  return TRUE;
  	}
  	else
  	{
  	  return FALSE;
  	}
  }
  
  if (tsp->field_choice == 1)
  {
  	/* look for strain field */
  	if (biop->org->orgname == NULL)
  	{
  	  if (tsp->when_string_not_found)
  	  {
  	  	return TRUE;
  	  }
  	  else
  	  {
  	  	return FALSE;
  	  }
  	}
  	for (mod = biop->org->orgname->mod;
  	     mod != NULL && mod->subtype != ORGMOD_strain;
  	     mod = mod->next)
  	{
  	  prev_mod = mod;
  	}
  	if (mod != NULL)
  	{
  	  search_text = mod->subname;
  	}
  }
  else if (tsp->field_choice == 2)
  {
    /* look for biosource comment */
  	for (mod = biop->org->orgname->mod;
  	     mod != NULL && mod->subtype != 255;
  	     mod = mod->next)
  	{
  	  prev_mod = mod;
  	}
  	if (mod != NULL)
  	{
  	  search_text = mod->subname;
  	}
  	else
  	{
  	  for (ssp = biop->subtype; ssp != NULL && ssp->subtype != 255; ssp = ssp->next)
  	  {
  	  	prev_ssp = ssp;
  	  }
  	  if (ssp != NULL)
  	  {
  	  	search_text = ssp->name;
  	  }
  	}
  }
  else
  {
  	return FALSE;
  }
  if (search_text != NULL)
  {
  	if (tsp->case_insensitive)
  	{
  	  string_found = StringISearch (search_text, tsp->find_this);
  	}
  	else
  	{
  	  string_found = StringSearch (search_text, tsp->find_this);
  	}
  	if (string_found != NULL)
  	{
  	  if (tsp->string_loc == 2 && string_found != search_text)
  	  {
  	  	string_found = NULL;
  	  }
  	  else if (tsp->string_loc == 3)
  	  {
  	  	while (string_found != NULL && string_found[StringLen (tsp->find_this)] != 0)
  	  	{
  	      if (tsp->case_insensitive)
  	      {
  	        string_found = StringISearch (string_found + 1, tsp->find_this);
          }
          else
          {
            string_found = StringSearch (string_found + 1, tsp->find_this);
          }
  	  	}
  	  }
  	}
  }
  
  if (string_found == NULL) 
  {
    if (tsp->when_string_not_found)
  	{
  	  rval = TRUE;
  	}
  	else
  	{
  	  rval = FALSE;
  	}
  }
  else
  {
    if (tsp->when_string_not_found)
  	{
  	  rval = FALSE;
  	}
  	else
  	{
  	  rval = TRUE;
  	  if (tsp->remove_found_text)
  	  {
  	  	if (string_found == search_text)
  	  	{
  	  	  if (StringLen (string_found) == StringLen (tsp->find_this))
  	  	  {
  	  	  	/* remove entire mod or ssp */
  	  	  	if (mod != NULL)
  	  	  	{
  	  	  	  if (prev_mod == NULL)
  	  	  	  {
  	  	  	  	biop->org->orgname->mod = mod->next;
  	  	  	  }
  	  	  	  else
  	  	  	  {
  	  	  	  	prev_mod->next = mod->next;
  	  	  	  }
  	  	  	  mod->next = NULL;
  	  	  	  OrgModFree (mod);
  	  	  	}
  	  	  	else if (ssp != NULL)
  	  	  	{
  	  	  	  if (prev_ssp == NULL)
  	  	  	  {
  	  	  	  	biop->subtype = ssp->next;
  	  	  	  }
  	  	  	  else
  	  	  	  {
  	  	  	  	prev_ssp->next = ssp->next;
  	  	  	  	ssp->next = NULL;
  	  	  	  	SubSourceFree (ssp);
  	  	  	  }
  	  	  	  ssp->next = NULL;
  	  	  	  SubSourceFree (ssp);
  	  	  	}
  	  	  }
  	  	  else
  	  	  {
  	  	  	/* remove first part of string and shift remainder */
  	  	  	destp = search_text;
  	  	  	for (cp = search_text + StringLen (tsp->find_this); *cp != 0; cp++)
  	  	  	{
  	  	  	  *destp++ = *cp;
  	  	  	}
  	  	  	*destp = 0;
  	  	  }
  	  	}
  	  	else
  	  	{
  	  	  /* keep first part of string, skip match, keep remainder */
  	  	  destp = string_found;
  	  	  for (cp = string_found + StringLen (tsp->find_this); *cp != 0; cp++)
  	  	  {
  	  	  	*destp++ = *cp;
  	  	  }
  	  	  *destp = 0;
  	  	}
  	  }
  	}
  }
  return rval;
}

static void AddTypeStrainCommentsProc (BioSourcePtr biop, Pointer userdata)
{
  SubSourcePtr       ssp, last_ssp;
  TypeStrainPtr      tsp;
  CharPtr            tmp;
  CharPtr            short_format = "type strain of %s";
  CharPtr            long_format = "%s; type strain of %s";

  if (biop == NULL || biop->org == NULL || biop->org->taxname == NULL) return;

  tsp = (TypeStrainPtr) userdata;
  
  if (! MeetsTypeStrainConstraint (biop, tsp)) return;
  
  ssp = biop->subtype;
  last_ssp = NULL;
  while (ssp != NULL && ssp->subtype != 255) {
    last_ssp = ssp;
    ssp = ssp->next;
  }
  if (ssp != NULL) {
    if (StringStr (ssp->name, "type strain of") != NULL) return;
    tmp = (CharPtr) MemNew (StringLen (long_format) + StringLen (ssp->name) 
                            + StringLen (biop->org->taxname) + 1);
    if (tmp != NULL) {
      sprintf (tmp, long_format, ssp->name, biop->org->taxname);
      MemFree (ssp->name);
      ssp->name = tmp;
    }
  } else {
    ssp = SubSourceNew ();
    if (ssp != NULL) {
      ssp->subtype = 255;
      tmp = (CharPtr) MemNew (StringLen (short_format)
                            + StringLen (biop->org->taxname) + 1);
      if (tmp != NULL) {
        sprintf (tmp, short_format, biop->org->taxname);
        ssp->name = tmp;
      }
      if (last_ssp == NULL) {
        biop->subtype = ssp;
      } else {
        last_ssp->next = ssp;
      }
    }
  }
}

static void CleanupTypeStrainForm (GraphiC g, VoidPtr data)

{
  TypeStrainPtr tsp;

  tsp = (TypeStrainPtr) data;
  if (tsp != NULL)
  {
  	tsp->find_this = MemFree (tsp->find_this);
  }
  MemFree (tsp);
  StdCleanupFormProc (g, data);
}

static void AddTypeStrainCommentsWithConstraintProc (ButtoN b)
{
  TypeStrainPtr tsp;
  SeqEntryPtr   sep;
  
  tsp = (TypeStrainPtr) GetObjectExtra (b);
  if (tsp == NULL) return;
  sep = GetTopSeqEntryForEntityID (tsp->input_entityID);
  if (sep == NULL) return;

  tsp->find_this = SaveStringFromText (tsp->find_this_txt);  
  tsp->when_string_not_found = GetStatus (tsp->when_string_not_found_btn);
  tsp->case_insensitive = GetStatus (tsp->case_insensitive_btn);
  tsp->string_loc = GetValue (tsp->string_loc_grp);
  tsp->field_choice = GetValue (tsp->field_choice_popup);
  tsp->remove_found_text = GetStatus (tsp->remove_found_text_btn);
  
  /* Visit each bioseq to remove redundant proprotein misc feats */
  VisitBioSourcesInSep (sep, tsp, AddTypeStrainCommentsProc);

  ObjMgrSetDirtyFlag (tsp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, tsp->input_entityID, 0, 0);
  Remove (tsp->form);
  ArrowCursor ();
  Update ();
}

extern void AddTypeStrainCommentsWithConstraint (IteM i)
{
  BaseFormPtr    bfp;
  TypeStrainPtr  tsp;
  WindoW         w;
  GrouP          h, k, l, m, c;
  ButtoN         b;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
	
  tsp = (TypeStrainPtr) MemNew (sizeof (TypeStrainData));
  if (tsp == NULL) return;
  tsp->input_entityID = bfp->input_entityID;

  w = FixedWindow (-50, -33, -10, -10, "Remove Sequences From Alignment", StdCloseWindowProc);
  if (w == NULL) {
	MemFree (tsp);
	return;
  }
  tsp->form = (ForM) w;
  SetObjectExtra (w, tsp, CleanupTypeStrainForm);
  
  h = HiddenGroup (w, 1, 0, NULL);
  k = HiddenGroup (h, 2, 0, NULL);

  StaticPrompt (k, "When this text is present", 0, dialogTextHeight, systemFont, 'c');
  tsp->find_this_txt = DialogText (k, "", 15, NULL);
  l = HiddenGroup (h, 2, 0, NULL);
  StaticPrompt (l, "In ", 0, dialogTextHeight, systemFont, 'c');
  tsp->field_choice_popup = PopupList (l, TRUE, NULL);
  PopupItem (tsp->field_choice_popup, "Strain");
  PopupItem (tsp->field_choice_popup, "Comment");
  SetValue (tsp->field_choice_popup, 1);
  tsp->string_loc_grp = HiddenGroup (h, 3, 0, NULL);
  RadioButton (tsp->string_loc_grp, "Anywhere in field");
  RadioButton (tsp->string_loc_grp, "At beginning of field");
  RadioButton (tsp->string_loc_grp, "At end of field");
  SetValue (tsp->string_loc_grp, 3);
  m = HiddenGroup (h, 2, 0, NULL);
  tsp->case_insensitive_btn = CheckBox (m, "Case Insensitive", NULL);
  tsp->when_string_not_found_btn = CheckBox (m, "When string is not found", NULL);  
  tsp->remove_found_text_btn = CheckBox (m, "Remove found text", NULL);

  c = HiddenGroup (h, 4, 0, NULL);
  b = DefaultButton (c, "Accept", AddTypeStrainCommentsWithConstraintProc);
  SetObjectExtra (b, tsp, NULL);
  b = PushButton (c, "Cancel", StdCancelButtonProc); 
  SetObjectExtra (b, tsp, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) k, (HANDLE) l, (HANDLE) tsp->string_loc_grp, 
                (HANDLE) m, (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Update ();
}

extern void AddTypeStrainCommentsToAll (IteM i)
{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  /* Visit each bioseq to remove redundant proprotein misc feats */
  VisitBioSourcesInSep (sep, NULL, AddTypeStrainCommentsProc);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ArrowCursor ();
  Update ();
}

extern void SqnNewAlign (BioseqPtr bsp1, BioseqPtr bsp2, SeqAlignPtr PNTR salp)
{
  BLAST_SummaryOptions *options = NULL;
  Uint1 mol_was;

  if (bsp1 == NULL || bsp2 == NULL || salp == NULL) return;

  *salp = NULL;
  if (ISA_na (bsp1->mol) != ISA_na (bsp2->mol)) return;

  mol_was = bsp2->mol;
  bsp2->mol = bsp1->mol;
  BLAST_SummaryOptionsInit(&options);

  options->filter_string = StringSave ("F");
  BLAST_TwoSequencesSearch(options, bsp1, bsp2, salp);
  bsp2->mol = mol_was;
  BLAST_SummaryOptionsFree(options);
  
}

/* This section of code is for the Remove Sequences From Alignments function. */

typedef struct alignmentsequencelist {
  SeqIdPtr sip;
  Char     descr[255];
} AlignmentSequenceListData, PNTR AlignmentSequenceListPtr;

typedef struct removeseqfromaligndata {
  FORM_MESSAGE_BLOCK
  LisT        sequence_list_ctrl;
  ValNodePtr  sequence_list;
  SeqEntryPtr sep;
  Boolean     remove_all_from_alignments;
  Boolean     no_remove_all_from_alignments;
  Boolean     remove_all_products;
  Boolean     no_remove_all_products;
} RemoveSeqFromAlignData, PNTR RemoveSeqFromAlignPtr;

/* This function will remove DenDiag and pairwise alignments if they contain
 * the sequence identified by sip, otherwise it will remove the sequence from
 * the alignment.
 */
static SeqAlignPtr RemoveOneSequenceFromAlignment (SeqIdPtr sip, SeqAlignPtr salphead)
{
  Uint4       seqid_order;
  SeqIdPtr    tmpsip;
  SeqAlignPtr salp, salp_next, prev_salp, remove_salp, last_remove;
  
  if (!FindSeqIdinSeqAlign (salphead, sip)) return NULL;
  
  salp = salphead;
  prev_salp = NULL;
  remove_salp = NULL;
  last_remove = NULL;
  while (salp != NULL)
  {
    salp_next = salp->next;
    tmpsip = SeqIdPtrFromSeqAlign (salp);
    seqid_order = SeqIdOrderInBioseqIdList(sip, tmpsip);
    if (seqid_order == 0)
    {
      /* do nothing for this subalignment */
      prev_salp = salp;
    }
    else if (salp->dim == 2 || salphead->segtype ==1)
    {
      /* This is for a pairwise alignment or a DENDIAG alignment */
      if (prev_salp == NULL)
      {
      	salphead = salp->next;
      }
      else
      {
      	prev_salp->next = salp->next;
      }
      /* save the alignments that we want to free in a list and get rid of them
       * at the end - freeing them beforehand causes problems with listing the
       * IDs in the alignment.
       */
      salp->next = NULL;
      if (remove_salp == NULL)
      {
      	remove_salp = salp;
      }
      else
      {
      	last_remove->next = salp;
      }
      last_remove = salp;
    }
    else 
    {
      SeqAlignBioseqDeleteById (salphead, sip);  
      prev_salp = salp;
    }
    salp = salp_next;
  }
  /* Now we can free the alignment */
  SeqAlignFree (remove_salp);
  return salphead;
}

static void RemoveSequenceFromAlignmentsCallback (SeqAnnotPtr sap, Pointer userdata)
{
  SeqAlignPtr salp;
  SeqIdPtr    sip;

  if (sap == NULL || sap->type != 2 || userdata == NULL) return;
  salp = (SeqAlignPtr) sap->data;
  if (salp == NULL) return;
  sip = (SeqIdPtr) userdata;
  sap->data = RemoveOneSequenceFromAlignment (sip, salp);
  /* if we've deleted all of the alignments, get rid of the annotation as well */
  if (sap->data == NULL)
  {
  	sap->idx.deleteme = TRUE;
  }
}

typedef struct checkforremovesequencefromalignments
{
  Boolean  found_problem;
  SeqIdPtr sip;
} CheckForRemoveSequenceFromAlignmentsData, PNTR CheckForRemoveSequenceFromAlignmentsPtr;

/* This is the callback function for looking for pairwise alignments.
/* If we delete the first sequence in a pairwise alignment, we end up deleting
 * the entire alignment because that sequence is paired with every other sequence.
 */
static void CheckForRemoveSequenceFromAlignmentsProblemsCallback (SeqAnnotPtr sap, Pointer userdata)
{
  CheckForRemoveSequenceFromAlignmentsPtr p;
  SeqAlignPtr salphead, salp;
  Uint4       seqid_order;
  SeqIdPtr    tmpsip;
  
  if (sap == NULL || sap->type != 2
      || (p = (CheckForRemoveSequenceFromAlignmentsPtr)userdata) == NULL
      || p->found_problem)
  {
  	return;
  }
  salphead = (SeqAlignPtr) sap->data;
  if (salphead == NULL) return;
  
  if (!FindSeqIdinSeqAlign (salphead, p->sip))
  {
  	return;
  }
  for (salp = salphead; salp != NULL; salp = salp->next)
  {
    tmpsip = SeqIdPtrFromSeqAlign (salp);
    seqid_order = SeqIdOrderInBioseqIdList(p->sip, tmpsip);
    if (seqid_order == 0)
    {
      continue;
    }
    else if (seqid_order == 1 && salp->dim == 2)
    {
      p->found_problem = TRUE;      
    }
  }
}

static void DoRemoveSequencesFromAlignment (ButtoN b)
{
  RemoveSeqFromAlignPtr    rp;
  WindoW                   w;
  ValNodePtr               vnp;
  Int2                     val;
  AlignmentSequenceListPtr aslp;
  CheckForRemoveSequenceFromAlignmentsData data;
  
  if (b == NULL) return;
  rp = (RemoveSeqFromAlignPtr) GetObjectExtra (b);
  if (rp == NULL) return;
  
  w = (WindoW) rp->form;
  Hide (w);
  /* first, check for pairwise alignments */
  val = 1;
  for (vnp = rp->sequence_list; vnp != NULL; vnp = vnp->next) {
    aslp = vnp->data.ptrvalue;
	if (aslp == NULL) continue;
	if (GetItemStatus (rp->sequence_list_ctrl, val)) {
	  data.sip = aslp->sip;
	  data.found_problem = FALSE;
	  VisitAnnotsInSep (rp->sep, (Pointer) &data, CheckForRemoveSequenceFromAlignmentsProblemsCallback);
	  if (data.found_problem)
	  {
	  	Message (MSG_ERROR, "One of the selected sequences is the first in a pairwise alignment."
	  	"  You must convert the alignment to a multiple alignment before trying to remove this sequence.");
        Remove (rp->form);  
        return;
	  }
	}
	val++;
  }

  val = 1;
  for (vnp = rp->sequence_list; vnp != NULL; vnp = vnp->next) {
    aslp = vnp->data.ptrvalue;
	if (aslp == NULL) continue;
	if (GetItemStatus (rp->sequence_list_ctrl, val)) {
	  VisitAnnotsInSep (rp->sep, (Pointer) aslp->sip, RemoveSequenceFromAlignmentsCallback);
	}
	val++;
  }
 
  ValNodeFree (rp->sequence_list);
  rp->sequence_list = NULL; 
  DeleteMarkedObjects (rp->input_entityID, 0, NULL);
  ObjMgrSetDirtyFlag (rp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, rp->input_entityID, 0, 0);
  Remove (rp->form);  
}

/* This function is used so that a sequence ID will only appear once in the list,
 * even if it appears in more than one alignment or subalignment.
 */
static Boolean IsIDAlreadyInList (SeqIdPtr sip, ValNodePtr list)
{
  ValNodePtr vnp;
  AlignmentSequenceListPtr aslp;
  
  if (sip == NULL) return FALSE;
  
  for (vnp = list; vnp != NULL; vnp = vnp->next)
  {
    aslp = (AlignmentSequenceListPtr) vnp->data.ptrvalue;
    if (aslp != NULL && SeqIdComp (aslp->sip, sip) == SIC_YES)
    {
      return TRUE;
    }
  }
  return FALSE;
}

/* This function creates the list of sequence IDs and descriptions to use in 
 * the Remove Sequences From Alignments dialog.
 */
static void ListSequencesInAlignmentsCallback (SeqAnnotPtr sap, Pointer userdata)
{
  SeqAlignPtr salp;
  SeqIdPtr    sip_list, sip, bsp_sip;
  ValNodePtr PNTR list;
  ValNodePtr  vnp; 
  AlignmentSequenceListPtr aslp;
  BioseqPtr                bsp;
  Int4                     offset;

  if (sap == NULL || sap->type != 2 || userdata == NULL) return;
  salp = (SeqAlignPtr) sap->data;
  while (salp != NULL) 
  {
    list = (ValNodePtr PNTR)userdata;
    sip_list = SeqAlignIDList (salp);
    if (sip_list == NULL) return;
    for (sip = sip_list; sip != NULL; sip = sip->next) {
      if (IsIDAlreadyInList (sip, *list)) continue;
      aslp = (AlignmentSequenceListPtr) MemNew (sizeof (AlignmentSequenceListData));
	  if (aslp == NULL) return;
	  aslp->sip = sip;
	  bsp = BioseqFindCore (sip);
	  if (bsp != NULL) {
		  aslp->descr[0] = 0;
		  aslp->descr[253] = 0;
		  offset = 0;
		  for (bsp_sip = bsp->id; bsp_sip != NULL && offset < 250; bsp_sip = bsp_sip->next) {
			if (aslp->descr[0] != 0) {
			  aslp->descr[offset] = '\t';
			  offset ++;
			}
		    SeqIdWrite (bsp_sip, aslp->descr + offset, PRINTID_TEXTID_ACCESSION, 254 - offset);
			offset = StringLen (aslp->descr);
		  }
	  } else {
        SeqIdWrite (sip, aslp->descr, PRINTID_TEXTID_ACCESSION, 254);	    
	  }
	  vnp = ValNodeNew (*list);
	  vnp->data.ptrvalue = aslp;
	  if (*list == NULL) {
		  *list = vnp;
	  }	  
    }
    salp = salp->next;
  }
}

static ValNodePtr ListSequencesInAlignments (SeqEntryPtr sep)
{
	ValNodePtr list = NULL;
    VisitAnnotsInSep (sep, (Pointer) &list, ListSequencesInAlignmentsCallback);
    return list;
}

static void CleanupRemoveSequencesFromAlignmentForm (
  GraphiC g,
  VoidPtr data
)

{
  RemoveSeqFromAlignPtr rp;

  rp = (RemoveSeqFromAlignPtr) data;
  if (rp != NULL) {
    if (rp->sequence_list != NULL)
    {
	  ValNodeFree (rp->sequence_list);
	  rp->sequence_list = NULL;
    }
  }
  StdCleanupFormProc (g, data);
}

extern void RemoveSequencesFromAlignment (IteM i)
{
  BaseFormPtr              bfp;
  WindoW                   w;
  RemoveSeqFromAlignPtr    rp;
  GrouP                    h, k, c;
  ButtoN                   b;
  ValNodePtr               vnp;
  AlignmentSequenceListPtr aslp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif

  if (bfp == NULL) return;

  rp = (RemoveSeqFromAlignPtr) MemNew (sizeof (RemoveSeqFromAlignData));
  if (rp == NULL) return;
  rp->input_entityID = bfp->input_entityID;
  rp->sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (rp->sep == NULL) {
	MemFree (rp);
	return;
  }

  rp->sequence_list = ListSequencesInAlignments (rp->sep);
  if (rp->sequence_list == NULL) {
    Message (MSG_ERROR, "There are no sequences in alignments");
	MemFree (rp);
	return;
  }

  w = FixedWindow (-50, -33, -10, -10, "Remove Sequences From Alignment", StdCloseWindowProc);
  if (w == NULL) {
	MemFree (rp);
	return;
  }
  rp->form = (ForM) w;
  SetObjectExtra (w, rp, CleanupRemoveSequencesFromAlignmentForm);
  
  h = HiddenGroup (w, -1, 0, NULL);
  k = HiddenGroup (h, 2, 0, NULL);

  rp->sequence_list_ctrl = MultiList (k, 16, 16, NULL);
  for (vnp = rp->sequence_list; vnp != NULL; vnp = vnp->next) {
    aslp = vnp->data.ptrvalue;
	if (aslp != NULL) {
      ListItem (rp->sequence_list_ctrl, aslp->descr);
	}
  }

  c = HiddenGroup (h, 4, 0, NULL);
  b = DefaultButton (c, "Accept", DoRemoveSequencesFromAlignment);
  SetObjectExtra (b, rp, NULL);
  b = PushButton (c, "Cancel", StdCancelButtonProc); 
  SetObjectExtra (b, rp, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) k, (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Update ();
}

/* End of Remove Sequences From Alignments function code. */

/* This section of code is used for removing sequences from the record. */

static void ListSequencesInSeqEntry (SeqEntryPtr sep, ValNodePtr PNTR list)
{
  BioseqPtr                bsp;
  BioseqSetPtr             bssp;
  ValNodePtr               vnp;
  AlignmentSequenceListPtr aslp;
  Int4                     offset;
  SeqIdPtr                 bsp_sip;
  
  if (sep == NULL) return;
  
  if (IS_Bioseq (sep))
  {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp == NULL) return;
    aslp = (AlignmentSequenceListPtr) MemNew (sizeof (AlignmentSequenceListData));
    if (aslp == NULL) return;
    aslp->sip = bsp->id;
    aslp->descr[0] = 0;
	aslp->descr[253] = 0;
    offset = 0;
    for (bsp_sip = bsp->id; bsp_sip != NULL && offset < 250; bsp_sip = bsp_sip->next) {
	  if (aslp->descr[0] != 0) {
	    aslp->descr[offset] = '\t';
	    offset ++;
	  }
      SeqIdWrite (bsp_sip, aslp->descr + offset, PRINTID_TEXTID_ACCESSION, 254 - offset);
      offset = StringLen (aslp->descr);
	}
    vnp = ValNodeNew (*list);
    if (vnp != NULL)
    {
      vnp->data.ptrvalue = aslp;
    }
    if (*list == NULL)
    {
      *list = vnp;
    }
  }
  else
  {
  	bssp = (BioseqSetPtr) sep->data.ptrvalue;
    for (sep = bssp->seq_set; sep != NULL; sep = sep->next) 
    {
      ListSequencesInSeqEntry (sep, list);
    }
  }
}

typedef struct bioseqinalignmentdata {
	Boolean   found;
	BioseqPtr lookingfor;
} BioseqInAlignmentData, PNTR BioseqInAlignmentPtr;

static Boolean IsBioseqInThisAlignment (SeqAlignPtr salp, BioseqPtr bsp)
{
  SeqIdPtr sip;
  Boolean found = FALSE;

  for (sip = bsp->id; sip != NULL && ! found; sip = sip->next) 
  {
    found = SeqAlignFindSeqId (salp, sip);
  }
  return found;
}

static void FindAlignmentCallback (SeqAnnotPtr sap, Pointer userdata)
{
  BioseqInAlignmentPtr biap;
  SeqAlignPtr          salp;

  if (sap == NULL || sap->type != 2 || userdata == NULL) 
  {
    return;
  }
  biap = (BioseqInAlignmentPtr) userdata;
  if (biap->found) return;
  salp = (SeqAlignPtr) sap->data;
  if (salp == NULL) return;
  biap->found = IsBioseqInThisAlignment (salp, biap->lookingfor);

}

static Boolean IsBioseqInAnyAlignment (BioseqPtr bsp, Uint2 input_entityID)
{
  SeqEntryPtr           topsep;
  BioseqInAlignmentData biad;

  topsep = GetTopSeqEntryForEntityID (input_entityID);
  biad.found = FALSE;
  biad.lookingfor = bsp;

  VisitAnnotsInSep (topsep, &biad, FindAlignmentCallback);
  return biad.found;
}

static void DoesBioseqHaveFeaturesWithProductsCallback (SeqFeatPtr sfp, Pointer userdata)
{
  ValNodePtr PNTR list;
  ValNodePtr vnp;
  
  if (sfp == NULL || userdata == NULL) return;
  list = (ValNodePtr PNTR) userdata;
  
  if (sfp->product != NULL)
  {
  	vnp = ValNodeNew (*list);
  	if (vnp != NULL)
  	{
  	  vnp->data.ptrvalue = sfp;
  	}
  	if (*list == NULL)
  	{
  	  *list = vnp;
  	}
  }
}

static void RemoveBioseq (BioseqPtr bsp, RemoveSeqFromAlignPtr rp);

static void RemoveBioseqProducts (ValNodePtr product_feature_list, RemoveSeqFromAlignPtr rp)
{
  ValNodePtr vnp;
  SeqFeatPtr sfp;
  BioseqPtr  bsp;
  
  for (vnp = product_feature_list; vnp != NULL; vnp = vnp->next)
  {
    sfp = (SeqFeatPtr) vnp->data.ptrvalue;
    if (sfp != NULL)
    {
  	  bsp = BioseqFindFromSeqLoc (sfp->product);
  	  sfp->product = SeqLocFree (sfp->product);
  	  RemoveBioseq (bsp, rp);
    }
  }
}

static void RemoveEmptyNucProtSet (SeqEntryPtr sep)
{
  BioseqSetPtr bssp;
  BioseqPtr    bsp;

  if (sep == NULL || !IS_Bioseq_set (sep)) return;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp->_class != BioseqseqSet_class_nuc_prot) return;

  for (sep = bssp->seq_set; sep != NULL; sep = sep->next)
  {
  	if (!IS_Bioseq (sep)) return;
  	bsp = sep->data.ptrvalue;
  	if (bsp != NULL && !bsp->idx.deleteme) return;
  }
  bssp->idx.deleteme = TRUE;
}

typedef struct removealnorproductans 
{
  WindoW  w;
  Boolean ans;
  Boolean do_all;
  Boolean done;
} RemoveAlnOrProductAnsData, PNTR RemoveAlnOrProductAnsPtr;

static void RemoveAlnOrProductYes (ButtoN b)
{
  RemoveAlnOrProductAnsPtr rp;
  
  rp = (RemoveAlnOrProductAnsPtr) GetObjectExtra (b);
  if (rp == NULL) return;
  rp->ans = TRUE;
  rp->do_all = FALSE;
  Remove (rp->w);
  rp->done = TRUE;
}

static void RemoveAlnOrProductYesAll (ButtoN b)
{
  RemoveAlnOrProductAnsPtr rp;
  
  rp = (RemoveAlnOrProductAnsPtr) GetObjectExtra (b);
  if (rp == NULL) return;
  rp->ans = TRUE;
  rp->do_all = TRUE;
  Remove (rp->w);
  rp->done = TRUE;
}

static void RemoveAlnOrProductNo (ButtoN b)
{
  RemoveAlnOrProductAnsPtr rp;
  
  rp = (RemoveAlnOrProductAnsPtr) GetObjectExtra (b);
  if (rp == NULL) return;
  rp->ans = FALSE;
  rp->do_all = FALSE;
  Remove (rp->w);
  rp->done = TRUE;
}

static void RemoveAlnOrProductNoAll (ButtoN b)
{
  RemoveAlnOrProductAnsPtr rp;
  
  rp = (RemoveAlnOrProductAnsPtr) GetObjectExtra (b);
  if (rp == NULL) return;
  rp->ans = FALSE;
  rp->do_all = TRUE;
  Remove (rp->w);
  rp->done = TRUE;
}

static Boolean GetRemoveAlignments (RemoveSeqFromAlignPtr rp, CharPtr idstr)
{
  RemoveAlnOrProductAnsData rd;

  GrouP                    g, h, c;
  ButtoN                   b;
  CharPtr                  prompt_fmt = "%s is part of an alignment - would you like to remove it from the alignment before deleting it?";
  CharPtr                  prompt_str = NULL;
  
  if (rp == NULL || idstr == NULL) return FALSE;
  if (rp->remove_all_from_alignments) return TRUE;
  if (rp->no_remove_all_from_alignments) return FALSE;

  prompt_str = (CharPtr) MemNew (sizeof (Char) * (StringLen (prompt_fmt) + StringLen (idstr)));
  if (prompt_str == NULL) return FALSE;
  sprintf (prompt_str, prompt_fmt, idstr);
  rd.w = ModalWindow(-20, -13, -10, -10, NULL);
  h = HiddenGroup(rd.w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  rd.done = FALSE;
  g = HiddenGroup (h, 1, 0, NULL);
  StaticPrompt (g, prompt_str, 0, popupMenuHeight, programFont, 'l');
  c = HiddenGroup (h, 4, 0, NULL);
  b = PushButton(c, "Yes", RemoveAlnOrProductYes);
  SetObjectExtra (b, &rd, NULL);
  b = PushButton(c, "Remove All", RemoveAlnOrProductYesAll);
  SetObjectExtra (b, &rd, NULL);
  b = DefaultButton(c, "No", RemoveAlnOrProductNo);
  SetObjectExtra (b, &rd, NULL);
  b = DefaultButton(c, "Remove None", RemoveAlnOrProductNoAll);
  SetObjectExtra (b, &rd, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
  prompt_str = MemFree (prompt_str);
  
  Show(rd.w); 
  Select (rd.w);
  rd.done = FALSE;
  while (!rd.done)
  {
    ProcessExternalEvent ();
    Update ();
  }
  ProcessAnEvent ();
  if (rd.do_all)
  {
    if (rd.ans)
    {
  	  rp->remove_all_from_alignments = TRUE;
  	  rp->no_remove_all_from_alignments = FALSE;
    }
    else
    {
  	  rp->remove_all_from_alignments = FALSE;
  	  rp->no_remove_all_from_alignments = TRUE;
    }
  }
  return rd.ans;
}


static Boolean GetRemoveProducts (RemoveSeqFromAlignPtr rp, CharPtr idstr)
{
  RemoveAlnOrProductAnsData rd;

  GrouP                    g, h, c;
  ButtoN                   b;
  CharPtr                  prompt_fmt = "%s contains features that have products (proteins, etc.).  Would you like to remove the product sequences?";
  CharPtr                  prompt_str = NULL;
  
  if (rp == NULL || idstr == NULL) return FALSE;
  if (rp->remove_all_products) return TRUE;
  if (rp->no_remove_all_products) return FALSE;

  prompt_str = (CharPtr) MemNew (sizeof (Char) * (StringLen (prompt_fmt) + StringLen (idstr)));
  if (prompt_str == NULL) return FALSE;
  sprintf (prompt_str, prompt_fmt, idstr);
  rd.w = ModalWindow(-20, -13, -10, -10, NULL);
  h = HiddenGroup(rd.w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  rd.done = FALSE;
  g = HiddenGroup (h, 1, 0, NULL);
  StaticPrompt (g, prompt_str, 0, popupMenuHeight, programFont, 'l');
  c = HiddenGroup (h, 4, 0, NULL);
  b = PushButton(c, "Yes", RemoveAlnOrProductYes);
  SetObjectExtra (b, &rd, NULL);
  b = PushButton(c, "Remove All", RemoveAlnOrProductYesAll);
  SetObjectExtra (b, &rd, NULL);
  b = DefaultButton(c, "No", RemoveAlnOrProductNo);
  SetObjectExtra (b, &rd, NULL);
  b = DefaultButton(c, "Remove None", RemoveAlnOrProductNoAll);
  SetObjectExtra (b, &rd, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
  prompt_str = MemFree (prompt_str);
  
  Show(rd.w); 
  Select (rd.w);
  rd.done = FALSE;
  while (!rd.done)
  {
    ProcessExternalEvent ();
    Update ();
  }
  ProcessAnEvent ();
  if (rd.do_all)
  {
    if (rd.ans)
    {
  	  rp->remove_all_products = TRUE;
  	  rp->no_remove_all_products = FALSE;
    }
    else
    {
  	  rp->remove_all_products = FALSE;
  	  rp->no_remove_all_products = TRUE;
    }
  }
  return rd.ans;
}


static void RemoveBioseq (BioseqPtr bsp, RemoveSeqFromAlignPtr rp)
{
  ValNodePtr   product_feature_list = NULL;
  Char         str [128];
  SeqEntryPtr  sep;
  
  if (bsp == NULL || rp == NULL) return;	
  
  SeqIdWrite (bsp->id, str, PRINTID_REPORT, sizeof (str));

  if (IsBioseqInAnyAlignment (bsp, rp->input_entityID))
  {
    if (GetRemoveAlignments (rp, str))
    {
 	  VisitAnnotsInSep (rp->sep, (Pointer) bsp->id, RemoveSequenceFromAlignmentsCallback);
    }
  }
  VisitFeaturesOnBsp (bsp, &product_feature_list, DoesBioseqHaveFeaturesWithProductsCallback);
  if (product_feature_list != NULL)
  {
    if (GetRemoveProducts (rp, str))
    {
      RemoveBioseqProducts (product_feature_list, rp);
    }
  }
        
  bsp->idx.deleteme = TRUE;
  /* remove nuc-prot set if we are deleting the nucleotide and its proteins */
  sep = GetBestTopParentForData (rp->input_entityID, bsp);
  RemoveEmptyNucProtSet (sep);

  ValNodeFree (product_feature_list);
  
}


static void DoRemoveSequencesFromRecord (ButtoN b)
{
  RemoveSeqFromAlignPtr    rp;
  WindoW                   w;
  ValNodePtr               vnp;
  Int2                     val;
  AlignmentSequenceListPtr aslp;
  BioseqPtr                bsp;
  
  if (b == NULL) return;
  rp = (RemoveSeqFromAlignPtr) GetObjectExtra (b);
  if (rp == NULL) return;
  
  w = (WindoW) rp->form;
  Hide (w);

  val = 1;
  for (vnp = rp->sequence_list; vnp != NULL; vnp = vnp->next) {
    aslp = vnp->data.ptrvalue;
	if (aslp == NULL) continue;
	if (GetItemStatus (rp->sequence_list_ctrl, val)) {
	  bsp = BioseqFind (aslp->sip);
	  if (bsp != NULL)
	  {
	    RemoveBioseq (bsp, rp);
	  }
	}
	val++;
  }
 
  ValNodeFree (rp->sequence_list);
  rp->sequence_list = NULL; 
  DeleteMarkedObjects (rp->input_entityID, 0, NULL);
  ObjMgrSetDirtyFlag (rp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, rp->input_entityID, 0, 0);
  Remove (rp->form);  
}

extern void RemoveSequencesFromRecord (IteM i)
{
  BaseFormPtr              bfp;
  WindoW                   w;
  RemoveSeqFromAlignPtr    rp;
  GrouP                    h, k, c;
  ButtoN                   b;
  ValNodePtr               vnp;
  AlignmentSequenceListPtr aslp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif

  if (bfp == NULL) return;

  rp = (RemoveSeqFromAlignPtr) MemNew (sizeof (RemoveSeqFromAlignData));
  if (rp == NULL) return;
  rp->input_entityID = bfp->input_entityID;
  rp->sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (rp->sep == NULL) {
	MemFree (rp);
	return;
  }
  ListSequencesInSeqEntry (rp->sep, &rp->sequence_list);
  if (rp->sequence_list == NULL) {
    Message (MSG_ERROR, "There are no sequences in alignments");
	MemFree (rp);
	return;
  }
  
  rp->remove_all_from_alignments = FALSE;
  rp->remove_all_products = FALSE;
  rp->no_remove_all_from_alignments = FALSE;
  rp->no_remove_all_products = FALSE;
  
  w = FixedWindow (-50, -33, -10, -10, "Remove Sequences From Record", StdCloseWindowProc);
  if (w == NULL) {
	MemFree (rp);
	return;
  }
  rp->form = (ForM) w;
  SetObjectExtra (w, rp, CleanupRemoveSequencesFromAlignmentForm);
  
  h = HiddenGroup (w, -1, 0, NULL);
  k = HiddenGroup (h, 2, 0, NULL);

  rp->sequence_list_ctrl = MultiList (k, 16, 16, NULL);
  for (vnp = rp->sequence_list; vnp != NULL; vnp = vnp->next) {
    aslp = vnp->data.ptrvalue;
	if (aslp != NULL) {
      ListItem (rp->sequence_list_ctrl, aslp->descr);
	}
  }

  c = HiddenGroup (h, 4, 0, NULL);
  b = DefaultButton (c, "Accept", DoRemoveSequencesFromRecord);
  SetObjectExtra (b, rp, NULL);
  b = PushButton (c, "Cancel", StdCancelButtonProc); 
  SetObjectExtra (b, rp, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) k, (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Update ();  
}
