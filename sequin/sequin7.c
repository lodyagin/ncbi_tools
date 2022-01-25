/*   sequin7.c
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
* File Name:  sequin7.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   1/3/98
*
* $Revision: 6.44 $
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

#ifndef CODECENTER
static char *date_of_compilation = __DATE__;
#else
static char *date_of_compilation = "today";
#endif

#include "sequin.h"
#include <gather.h>
#include <edutil.h>
#include <cdrgn.h>
#include <subutil.h>
#include <tofasta.h>
#include <vsm.h>
#include <document.h>
#include <maputil.h>
#include <asn2ff.h>
#include <ffprint.h>
#include <bspview.h>
#include <findrepl.h>
#include <toasn3.h>
#include <taxutil.h>
#include <utilpub.h>
#include <salsap.h>

NLM_EXTERN SeqEntryPtr FastaToSeqEntryInternal
(
 VoidPtr input,          /* input pointer (file or memory) */
 Int4 type,              /* type of inquiry FASTA_MEM_IO or FASTA_FILE_IO */
 CharPtr PNTR last_char, /* returned pointer to next FASTA sequence */
 Boolean is_na,          /* type of sequence */
 CharPtr PNTR errormsg,  /* error messge for debugging */
 Boolean parseSeqId,     /* Parse SeqID from def line */
 CharPtr special_symbol     /* Returns special symbol if no SeqEntry */
 );

NLM_EXTERN SeqEntryPtr SequinFastaToSeqEntryEx 
  (
    FILE *fp,               /* file to get sequence from */ 
    Boolean is_na,          /* type of sequence */
    CharPtr PNTR errormsg,  /* error message for debugginq */
    Boolean parseSeqId,     /* Parse SeqID from def line */
    CharPtr special_symbol  /* Returns special symbol if no SeqEntry */
  )
{
  return FastaToSeqEntryInternal((void *)fp, 2, NULL,is_na, errormsg, parseSeqId, special_symbol);	 
}

static FonT  titleFont = NULL;

#ifndef WIN_MAC
void CreateSqnInitialFormMenus (WindoW w)

{
  BaseFormPtr   bfp;
  MenU          m;

  bfp = (BaseFormPtr) GetObjectExtra (w);
  if (bfp != NULL) {
    m = PulldownMenu (w, "File");
    AddAboutAndHelpMenuItems (m);
    if (bfp->importform != NULL || bfp->exportform != NULL) {
      if (bfp->importform != NULL) {
        FormCommandItem (m, "Import...", bfp, VIB_MSG_IMPORT);
      }
      if (bfp->exportform != NULL) {
        FormCommandItem (m, "Export...", bfp, VIB_MSG_EXPORT);
      }
      SeparatorItem (m);
    }
    FormCommandItem (m, "Quit", bfp, VIB_MSG_QUIT);
    m = PulldownMenu (w, "Edit");
    FormCommandItem (m, CUT_MENU_ITEM, bfp, VIB_MSG_CUT);
    FormCommandItem (m, COPY_MENU_ITEM, bfp, VIB_MSG_COPY);
    FormCommandItem (m, PASTE_MENU_ITEM, bfp, VIB_MSG_PASTE);
    FormCommandItem (m, CLEAR_MENU_ITEM, bfp, VIB_MSG_DELETE);
  }
}
#endif

static void DefaultMessageProc (ForM f, Int2 mssg)

{
  StdEditorProcsPtr  sepp;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    if (sepp->handleMessages != NULL) {
      sepp->handleMessages (f, mssg);
    }
  }
}

typedef struct startupform {
  FORM_MESSAGE_BLOCK
} StartupForm, PNTR StartupFormPtr;

static void ChangeDestination (GrouP g)

{
  Char  str [64];
  Int2  val;

  val = GetValue (g);
  switch (val) {
    case 1 :
      RemoveAppProperty ("SequinUseEMBLStyle");
      RemoveAppProperty ("SequinUseDDBJStyle");
      if (GetAppParam ("SEQUIN", "PREFERENCES", "DATABASE", NULL, str, sizeof (str))) {
        if (! StringHasNoText (str)) {
          if (StringICmp (str, "GenBank") != 0) {
            SetAppParam ("SEQUIN", "PREFERENCES", "DATABASE", "GenBank");
          }
        }
      }
      break;
    case 2 :
      SetAppProperty ("SequinUseEMBLStyle", (void *) 1024);
      RemoveAppProperty ("SequinUseDDBJStyle");
      SetAppParam ("SEQUIN", "PREFERENCES", "DATABASE", "EMBL");
      break;
    case 3 :
      RemoveAppProperty ("SequinUseEMBLStyle");
      SetAppProperty ("SequinUseDDBJStyle", (void *) 1024);
      SetAppParam ("SEQUIN", "PREFERENCES", "DATABASE", "DDBJ");
      break;
    default :
      break;
  }
  SetupBioseqPageList ();
}

static void CenterString (RectPtr rptr, CharPtr text, FonT fnt, Int2 inc)

{
  if (fnt != NULL) {
    SelectFont (fnt);
  }
  rptr->bottom = rptr->top + LineHeight ();
  DrawString (rptr, text, 'c', FALSE);
  rptr->top = rptr->bottom + inc;
}

extern void DrawAbout (PaneL p)

{
  RecT     r;
  CharPtr  ptr;
  Char     sequinServices [60];
  Char     sequinVersion [60];


  if (titleFont == NULL) {
#ifdef WIN_MAC
    titleFont = GetFont ("Geneva", 18, TRUE, TRUE, FALSE, "");
#endif
#ifdef WIN_MSWIN
    titleFont = GetFont ("Arial", 24, TRUE, TRUE, FALSE, "");
#endif
#ifdef WIN_MOTIF
    titleFont = GetFont ("Courier", 24, TRUE, TRUE, FALSE, "");
#endif
  }

  sprintf (sequinVersion, "Sequin Application Version %s", SEQUIN_APP_VERSION);
  ptr = "Standard Release";
/*#ifdef USE_ENTREZ*/
  if (useEntrez || useBlast) {
    ptr = "Network Aware";
  }
/*#endif*/
/*#ifdef INTERNAL_NCBI_SEQUIN*/
  if (indexerVersion) {
    ptr = "Indexer Services";
  }
/*#endif*/
  if (genomeCenter != NULL) {
    ptr = "Genome Center";
  }
  sprintf (sequinServices, "%s [%s]", ptr, date_of_compilation);

  ObjectRect (p, &r);
  InsetRect (&r, 4, 4);
  r.top += 10;
  Blue ();
  CenterString (&r, "Sequin", titleFont, 5);
  CenterString (&r, sequinVersion, programFont, 5);
  CenterString (&r, sequinServices, programFont, 10);
  CenterString (&r, "National Center for Biotechnology Information", systemFont, 5);
  CenterString (&r, "National Library of Medicine", systemFont, 5);
  CenterString (&r, "National Institutes of Health", systemFont, 10);
  CenterString (&r, "(301) 496-2475", systemFont, 5);
  CenterString (&r, "info@ncbi.nlm.nih.gov", systemFont, 0);
}

extern Int2 AboutBoxWidth (void)

{
  Int2     max;
  CharPtr  ptr;
  Char     sequinServices [60];
  Char     sequinVersion [60];
  Int2     wid;


  if (titleFont == NULL) {
#ifdef WIN_MAC
    titleFont = GetFont ("Geneva", 18, TRUE, TRUE, FALSE, "");
#endif
#ifdef WIN_MSWIN
    titleFont = GetFont ("Arial", 24, TRUE, TRUE, FALSE, "");
#endif
#ifdef WIN_MOTIF
    titleFont = GetFont ("Courier", 24, TRUE, TRUE, FALSE, "");
#endif
  }

  sprintf (sequinVersion, "Sequin Application Version %s", SEQUIN_APP_VERSION);
  ptr = "Standard Release";
/*#ifdef USE_ENTREZ*/
  if (useEntrez || useBlast) {
    ptr = "Network Aware";
  }
/*#endif*/
/*#ifdef INTERNAL_NCBI_SEQUIN*/
  if (indexerVersion) {
    ptr = "Indexer Services";
  }
/*#endif*/
  if (genomeCenter != NULL) {
    ptr = "Genome Center";
  }
  sprintf (sequinServices, "%s [%s]", ptr, date_of_compilation);

  SelectFont (titleFont);
  max = StringWidth ("Sequin");
  SelectFont (programFont);
  wid = StringWidth (sequinVersion);
  if (wid > max) {
    max = wid;
  }
  wid = StringWidth (sequinServices);
  if (wid > max) {
    max = wid;
  }
  SelectFont (systemFont);
  wid = StringWidth ("National Center for Biotechnology Information");
  if (wid > max) {
    max = wid;
  }
  max += 2 * stdCharWidth + 2;
  return max;
}

extern ForM CreateStartupForm (Int2 left, Int2 top, CharPtr title,
                               BtnActnProc startFa2htgs,
                               BtnActnProc startPhrap,
                               BtnActnProc buildContig,
                               BtnActnProc startNew,
                               BtnActnProc readExisting,
                               BtnActnProc fetchFromNet,
                               BtnActnProc showHelp,
                               BtnActnProc quitProgram,
                               WndActnProc activateForm)

{
  ButtoN          b;
  GrouP           c;
  GrouP           d;
  GrouP           k;
  PaneL           p;
  StartupFormPtr  sfp;
  Char            str [32];
  WindoW          w;
#ifndef WIN_MAC
  MenU            m;
#endif

  w = NULL;
  sfp = MemNew (sizeof (StartupForm));
  if (sfp != NULL) {
    w = FixedWindow (left, top, -10, -10, title, NULL);
    SetObjectExtra (w, sfp, StdCleanupFormProc);
    sfp->form = (ForM) w;
    sfp->formmessage = DefaultMessageProc;

#ifndef WIN_MAC
    m = PulldownMenu (w, "Misc");
    CommandItem (m, "Net Configure...", NetConfigureProc);
    if (useEntrez) {
      SeparatorItem (m);
      CommandItem (m, "Entrez Query...", EntrezQueryProc);
      if (extraServices) {
        SeparatorItem (m);
        CommandItem (m, "Process FASTA Nucleotide Updates", ParseInNucUpdates);
      }
    }
    if (useDesktop) {
      SeparatorItem (m);
      VSMAddToMenu (m, VSM_DESKTOP);
    }
#endif

    p = SimplePanel (w, AboutBoxWidth (), 14 * stdLineHeight, DrawAbout);

    k = HiddenGroup (w, 4, 0, NULL);
    SetGroupSpacing (k, 3, 10);
    StaticPrompt (k, "Database for submission", 0, stdLineHeight, programFont, 'l');
    d = HiddenGroup (k, 4, 0, ChangeDestination);
    RadioButton (d, "GenBank");
    RadioButton (d, "EMBL");
    RadioButton (d, "DDBJ");
    if (GetAppParam ("SEQUIN", "PREFERENCES", "DATABASE", NULL, str, sizeof (str))) {
      if (StringICmp (str, "GenBank") == 0) {
        SetValue (d, 1);
      } else if (StringICmp (str, "EMBL") == 0) {
        SetValue (d, 2);
      } else if (StringICmp (str, "DDBJ") == 0) {
        SetValue (d, 3);
      } else {
        SetValue (d, 1);
      }
    } else {
      SetValue (d, 1);
    }
    ChangeDestination (d);

    c = HiddenGroup (w, 1, 0, NULL);
    SetGroupSpacing (c, 10, 5);

    if (startFa2htgs != NULL) {
      b = PushButton (c, "New FA2HTGS Submission", startFa2htgs);
      SetObjectExtra (b, sfp, NULL);
    }
    if (startPhrap != NULL) {
      b = PushButton (c, "New PHRAP Submission", startPhrap);
      SetObjectExtra (b, sfp, NULL);
    }
    if (buildContig != NULL) {
      b = PushButton (c, "Read CONTIG Instructions", buildContig);
      SetObjectExtra (b, sfp, NULL);
    }
    b = PushButton (c, "Start New Submission", startNew);
    SetObjectExtra (b, sfp, NULL);
    b = PushButton (c, "Read Existing Record", readExisting);
    SetObjectExtra (b, sfp, NULL);
    if (fetchFromNet != NULL) {
      b = PushButton (c, "Download From Entrez", fetchFromNet);
      SetObjectExtra (b, sfp, NULL);
    }
    b = PushButton (c, "Show Help", showHelp);
    SetObjectExtra (b, sfp, NULL);
    b = PushButton (c, "Quit Program", quitProgram);
    SetObjectExtra (b, sfp, NULL);

    AlignObjects (ALIGN_CENTER, (HANDLE) p, (HANDLE) c, (HANDLE) k, NULL);

    RealizeWindow (w);

    if (activateForm != NULL) {
      SetActivate (w, activateForm);
    }
  }
  return (ForM) w;
}

typedef struct formatform {
  FORM_MESSAGE_BLOCK

  GrouP           package;
  GrouP           format;
  ButtoN          nexusButton;
  ButtoN          paupButton;
  TexT            numseqs;

  Int2            restoreFormatTo;
} FormatForm, PNTR FormatFormPtr;

static Boolean allowGenomicPlusCDNA = FALSE;

static void FormatBlockPtrToFormatForm (ForM f, Pointer data)

{
  FormatBlockPtr  fbp;
  FormatFormPtr   ffp;
  Char            str [32];

  ffp = (FormatFormPtr) GetObjectExtra (f);
  fbp = (FormatBlockPtr) data;
  if (ffp == NULL) return;
  if (fbp != NULL) {
    if (fbp->seqPackage > 0 && fbp->seqPackage <= NUM_SEQ_PKG) {
      if ((! allowGenomicPlusCDNA) && fbp->seqPackage >= SEQ_PKG_GENOMICCDNA) {
        SafeSetValue (ffp->package, fbp->seqPackage - 1);
      } else {
        SafeSetValue (ffp->package, fbp->seqPackage);
      }
      if (fbp->seqPackage <= SEQ_PKG_GENOMICCDNA || fbp->seqPackage == SEQ_PKG_GENBANK) {
        SafeDisable (ffp->nexusButton);
        SafeDisable (ffp->paupButton);
      } else {
        SafeEnable (ffp->nexusButton);
        SafeEnable (ffp->paupButton);
      }
    } else {
      SafeSetValue (ffp->package, SEQ_PKG_SINGLE);
      SafeDisable (ffp->nexusButton);
      SafeDisable (ffp->paupButton);
    }
    if (fbp->seqFormat > 0 && fbp->seqFormat <= NUM_SEQ_FMT) {
      SafeSetValue (ffp->format, fbp->seqFormat);
    } else {
      SafeSetValue (ffp->format, SEQ_FMT_FASTA);
    }
    if (fbp->numSeqs > 0) {
      IntToStr (fbp->numSeqs, str, 0, sizeof (str));
      SafeSetTitle (ffp->numseqs, str);
    } else {
      SafeSetTitle (ffp->numseqs, "");
    }
  } else {
    SafeSetValue (ffp->package, SEQ_PKG_SINGLE);
    SafeDisable (ffp->nexusButton);
    SafeDisable (ffp->paupButton);
    SafeSetValue (ffp->format, SEQ_FMT_FASTA);
    SafeSetTitle (ffp->numseqs, "");
    ffp->restoreFormatTo = SEQ_FMT_FASTA;
  }
}

static Pointer FormatFormToFormatBlockPtr (ForM f)

{
  FormatBlockPtr  fbp;
  FormatFormPtr   ffp;
  Char            str [32];
  Int2            val;

  fbp = NULL;
  ffp = (FormatFormPtr) GetObjectExtra (f);
  if (ffp == NULL) return NULL;
  fbp = (FormatBlockPtr) MemNew (sizeof (FormatBlock));
  if (fbp == NULL) return NULL;
  fbp->seqPackage = GetValue (ffp->package);
  if ((! allowGenomicPlusCDNA) && fbp->seqPackage >= SEQ_PKG_GENOMICCDNA) {
    (fbp->seqPackage)++;
  }
  fbp->seqFormat = GetValue (ffp->format);
  GetTitle (ffp->numseqs, str, sizeof (str));
  if (StrToInt (str, &val) && val > 0) {
    fbp->numSeqs = val;
  } else {
    fbp->numSeqs = 0;
  }
  return (Pointer) fbp;
}

static void EnableOrDisableFormats (GrouP g)

{
  FormatFormPtr  ffp;
  Int2           val;

  ffp = (FormatFormPtr) GetObjectExtra (g);
  if (ffp == NULL) return;
  val = GetValue (g);
  if ((! allowGenomicPlusCDNA) && val >= SEQ_PKG_GENOMICCDNA) {
    val++;
  }
  if (val <= SEQ_PKG_GENOMICCDNA || val == SEQ_PKG_GENBANK) {
    if (Enabled (ffp->nexusButton)) {
      ffp->restoreFormatTo = GetValue (ffp->format);
    }
    SafeSetValue (ffp->format, SEQ_FMT_FASTA);
    SafeDisable (ffp->nexusButton);
    SafeDisable (ffp->paupButton);
  } else {
    if (! Enabled (ffp->nexusButton)) {
      SafeSetValue (ffp->format, ffp->restoreFormatTo);
    }
    SafeEnable (ffp->nexusButton);
    SafeEnable (ffp->paupButton);
  }
}

extern ForM CreateFormatForm (Int2 left, Int2 top, CharPtr title,
                              BtnActnProc goToNext,
                              BtnActnProc goBack,
                              WndActnProc activateForm)

{
  ButtoN         b;
  GrouP          c;
  FormatFormPtr  ffp;
  GrouP          g1, g2;
  GrouP          h;
  PrompT         ppt;
  Char           str [32];
  WindoW         w;

  w = NULL;
  ffp = MemNew (sizeof (FormatForm));
  if (ffp != NULL) {
    w = FixedWindow (left, top, -10, -10, title, NULL);
    SetObjectExtra (w, ffp, StdCleanupFormProc);
    ffp->form = (ForM) w;
    ffp->toform = FormatBlockPtrToFormatForm;
    ffp->fromform = FormatFormToFormatBlockPtr;
    ffp->formmessage = DefaultMessageProc;

    SetGroupSpacing (w, 10, 10);

    h = HiddenGroup (w, -1, 0, NULL);
    SetGroupSpacing (h, 3, 10);

    g1 = HiddenGroup (h, 2, 0, NULL);

    allowGenomicPlusCDNA = FALSE;
    if (GetAppParam ("SEQUIN", "SETTINGS", "GENOMICPLUSTRANSCRIPTS", NULL, str, sizeof (str))) {
      if (StringICmp (str, "TRUE") == 0) {
        allowGenomicPlusCDNA = TRUE;
      }
    }

    ppt = StaticPrompt (g1, "Submission type", 0, 0, programFont, 'l');
    ffp->package = HiddenGroup (g1, 2, 0, EnableOrDisableFormats);
    SetObjectExtra (ffp->package, ffp, NULL);
    RadioButton (ffp->package, "Single sequence");
    RadioButton (ffp->package, "Segmented sequence");
    if (allowGenomicPlusCDNA) {
      RadioButton (ffp->package, "Genomic + Transcripts");
    }
    RadioButton (ffp->package, "Population study");
    RadioButton (ffp->package, "Phylogenetic study");
    RadioButton (ffp->package, "Mutation study");
/*#ifdef INTERNAL_NCBI_SEQUIN*/
    /* if (indexerVersion) { */
      RadioButton (ffp->package, "Batch submission");
    /* } */
/*#endif*/
    SetValue (ffp->package, SEQ_PKG_SINGLE);
    AlignObjects (ALIGN_MIDDLE, (HANDLE) ppt, (HANDLE) ffp->package, NULL);

    g2 = HiddenGroup (h, 2, 0, NULL);

    ppt = StaticPrompt (g2, "Sequence data format", 0, 0, programFont, 'l');
    ffp->format = HiddenGroup (g2, -1, 0, NULL);
    SetObjectExtra (ffp->format, ffp, NULL);
    RadioButton (ffp->format, "FASTA");
    ffp->paupButton = RadioButton (ffp->format, "Contiguous (FASTA+GAP, NEXUS, MACAW)");
    ffp->nexusButton = RadioButton (ffp->format, "Interleaved (PHYLIP, NEXUS)");
    Disable (ffp->nexusButton);
    Disable (ffp->paupButton);
    SetValue (ffp->format, SEQ_FMT_FASTA);
    ffp->restoreFormatTo = SEQ_FMT_FASTA;
    AlignObjects (ALIGN_MIDDLE, (HANDLE) ppt, (HANDLE) ffp->format, NULL);

    c = HiddenGroup (w, 4, 0, NULL);
    SetGroupSpacing (c, 10, 2);
    b = PushButton (c, " << Prev Form ", goBack);
    SetObjectExtra (b, ffp, NULL);
    b = PushButton (c, " Next Form >> ", goToNext);
    SetObjectExtra (b, ffp, NULL);

    AlignObjects (ALIGN_CENTER, (HANDLE) g1, (HANDLE) g2, (HANDLE) c, NULL);

    RealizeWindow (w);

    if (activateForm != NULL) {
      SetActivate (w, activateForm);
    }
  }
  return (ForM) w;
}

extern SequinBlockPtr SequinBlockFree (SequinBlockPtr sbp)

{
  if (sbp != NULL) {
    AuthorFree (sbp->contactperson);
    AuthListFree (sbp->citsubauthors);
    AffilFree (sbp->citsubaffil);
    MemFree (sbp->citsubtitle);
    DateFree (sbp->releasedate);
  }
  return NULL;
}

extern void ExciseString (CharPtr str, CharPtr from, CharPtr to)

{
  Char     ch;
  CharPtr  ptrf;
  CharPtr  ptrt;

  if (str == NULL || from == NULL || to == NULL) return;
  ptrf = StringStr (str, from);
  if (ptrf == NULL) return;
  ptrt = StringStr (ptrf, to);
  if (ptrt == NULL) return;
  ptrt += StringLen (to);
  ch = *ptrt;
  while (ch != '\0') {
    *ptrf = ch;
    ptrf++;
    ptrt++;
    ch = *ptrt;
  }
  *ptrf = '\0';
}

typedef struct geneextendlist {
  GeneRefPtr  grp;
  SeqLocPtr   slp;
  ObjMgrPtr   omp;
  Boolean     rsult;
  Char        label [41];
} GeneExtendList, PNTR GeneExtendPtr;

static Boolean GeneExtendFunc (GatherContextPtr gcp)

{
  BioseqPtr      bsp;
  GeneExtendPtr  gep;
  GeneRefPtr     grp;
  Boolean        hasNulls;
  ObjMgrTypePtr  omtp;
  SeqFeatPtr     sfp;
  SeqLocPtr      slp;
  Char           thislabel [41];

  if (gcp == NULL) return TRUE;

  gep = (GeneExtendPtr) gcp->userdata;
  if (gep == NULL ) return TRUE;

  thislabel [0] = '\0';

  if (gcp->thistype == OBJ_SEQFEAT) {
    sfp = (SeqFeatPtr) gcp->thisitem;
    if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE && sfp->data.value.ptrvalue != NULL) {
      grp = (GeneRefPtr) sfp->data.value.ptrvalue;
      omtp = ObjMgrTypeFind (gep->omp, gcp->thistype, NULL, NULL);
      if (omtp == NULL) {
        return TRUE;
      }
      if (omtp->labelfunc != NULL) {
        (*(omtp->labelfunc)) (gcp->thisitem, thislabel, 40, OM_LABEL_CONTENT);
      }
      if (thislabel [0] != '\0') {
        if (StringICmp (thislabel, gep->label) == 0) {
          if (SeqLocCompare (gep->slp, sfp->location) != SLC_NO_MATCH) {
            bsp = GetBioseqGivenSeqLoc (sfp->location, gcp->entityID);
            if (bsp != NULL) {
              slp = SeqLocMerge (bsp, sfp->location, gep->slp, TRUE, FALSE, FALSE);
              if (slp != NULL) {
                sfp->location = SeqLocFree (sfp->location);
                sfp->location = slp;
                if (bsp->repr == Seq_repr_seg) {
                  slp = SegLocToPartsEx (bsp, sfp->location, TRUE);
                  sfp->location = SeqLocFree (sfp->location);
                  sfp->location = slp;
                  hasNulls = LocationHasNullsBetween (sfp->location);
                  sfp->partial = (sfp->partial || hasNulls);
                }
                FreeAllFuzz (slp);
                gep->rsult = TRUE;
              }
            }
          }
          return FALSE;
        }
      }
    }
  }
  return TRUE;
}

extern Boolean ExtendGene (GeneRefPtr grp, SeqEntryPtr nsep, SeqLocPtr slp)

{
  GeneExtendList  gel;
  GatherScope     gs;
  ObjMgrTypePtr   omtp;
  SeqFeatPtr      sfp;

  if (grp == NULL || nsep == NULL || slp == NULL) return FALSE;
  gel.grp = grp;
  gel.slp = slp;
  gel.omp = ObjMgrGet ();
  gel.label [0] = '\0';
  gel.rsult = FALSE;
  omtp = ObjMgrTypeFind (gel.omp, OBJ_SEQFEAT, NULL, NULL);
  if (omtp != NULL && omtp->labelfunc != NULL) {
    sfp = SeqFeatNew ();
    if (sfp != NULL) {
      sfp->data.choice = SEQFEAT_GENE;
      sfp->data.value.ptrvalue = (Pointer) grp;
      (*(omtp->labelfunc)) ((Pointer) sfp, gel.label, 40, OM_LABEL_CONTENT);
      sfp->data.value.ptrvalue = NULL;
      SeqFeatFree (sfp);
    }
  }
  MemSet ((Pointer)(&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = TRUE;
  MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  GatherSeqEntry (nsep, (Pointer) &gel, GeneExtendFunc, &gs);
  return gel.rsult;
}

static void CreateGeneAndProtFeats (SeqEntryPtr nsep, SeqEntryPtr psep,
                                    SeqLocPtr slp, CdRegionPtr crp, CharPtr title,
                                    CharPtr best, size_t maxsize, CharPtr PNTR ttl)

{
  Char        activity [256];
  BioseqPtr   bsp;
  Char        ec [32];
  GeneRefPtr  grp;
  SeqLocPtr   gslp;
  Boolean     hasNulls;
  BioseqPtr   nbsp;
  BioseqPtr   pbsp;
  ProtRefPtr  prp;
  CharPtr     ptr;
  SeqFeatPtr  sfp;
  SeqIdPtr    sip;
  Char        str [256];

  if (nsep != NULL && psep != NULL && slp != NULL && crp != NULL && title != NULL) {
    if (best != NULL) {
      best [0] = '\0';
    }
    if (IS_Bioseq (nsep) && IS_Bioseq (psep)) {
      nbsp = (BioseqPtr) nsep->data.ptrvalue;
      pbsp = (BioseqPtr) psep->data.ptrvalue;
      if (nbsp != NULL && pbsp != NULL) {
        ptr = StringStr (title, "[gene=");
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
            StringNCpy_0 (best, str, maxsize);
            if (StringHasNoText (best)) {
              StringNCpy_0 (best, ptr, maxsize);
            }
            grp = CreateNewGeneRef (str, NULL, ptr, FALSE);
            if (grp != NULL) {
              if (ExtendGene (grp, nsep, slp)) {
                grp = GeneRefFree (grp);
              } else {
                sfp = CreateNewFeature (nsep, NULL, SEQFEAT_GENE, NULL);
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
                    bsp = nbsp;
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
        activity [0] = '\0';
        ec [0] = '\0';
        ptr = StringStr (title, "[function=");
        if (ptr != NULL) {
          StringNCpy_0 (activity, ptr + 10, sizeof (str));
          ptr = StringChr (activity, ']');
          if (ptr != NULL) {
            *ptr = '\0';
          }
        }
        ptr = StringStr (title, "[EC_number=");
        if (ptr != NULL) {
          StringNCpy_0 (ec, ptr + 11, sizeof (str));
          ptr = StringChr (ec, ']');
          if (ptr != NULL) {
            *ptr = '\0';
          }
        }
        ptr = StringStr (title, "[prot=");
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
            StringNCpy_0 (best, str, maxsize);
            if (StringHasNoText (best)) {
              StringNCpy_0 (best, ptr, maxsize);
            }
            prp = CreateNewProtRef (str, ptr, ec, activity);
            if (prp != NULL) {
              sfp = CreateNewFeature (psep, NULL, SEQFEAT_PROT, NULL);
              if (sfp != NULL) {
                sfp->data.value.ptrvalue = (Pointer) prp;
              }
            }
          }
        }
        ptr = StringStr (title, "[orf]");
        if (ptr != NULL) {
          crp->orf = TRUE;
        }
        ptr = StringStr (title, "[comment=");
        if (ptr != NULL) {
          StringNCpy_0 (str, ptr + 9, sizeof (str));
          ptr = StringChr (str, ']');
          if (ptr != NULL) {
            *ptr = '\0';
            if (ttl != NULL) {
              *ttl = StringSave (str);
            }
          }
        }
        ptr = StringStr (title, "[note=");
        if (ptr != NULL) {
          StringNCpy_0 (str, ptr + 6, sizeof (str));
          ptr = StringChr (str, ']');
          if (ptr != NULL) {
            *ptr = '\0';
            if (ttl != NULL && *ttl == NULL) {
              *ttl = StringSave (str);
            }
          }
        }
      }
    }
  }
}

static Boolean  intBoxUp;
static Boolean  intBoxRsult;

static void AcceptAskProc (ButtoN b)

{
  intBoxRsult = TRUE;
  intBoxUp = FALSE;
}

static void CancelAskProc (ButtoN b)

{
  intBoxRsult = FALSE;
  intBoxUp = FALSE;
}

static SeqLocPtr AskForInterval (SeqEntryPtr sep, BioseqPtr nuc, BioseqPtr prot)

{
  GrouP       c;
  DialoG      d;
  GrouP       g;
  GrouP       m;
  SeqIdPtr    sip;
  SeqLocPtr   slp;
  Char        str [128];
  ValNodePtr  vnp;
  WindoW      w;

  slp = NULL;
  if (sep == NULL || nuc == NULL || prot == NULL) return NULL;

  if (GetAppParam ("SEQUIN", "PREFERENCES", "ASKIFSUGGESTFAILED", NULL, str, sizeof (str))) {
    if (StringICmp (str, "FALSE") == 0) {
      sip = SeqIdFindWorst (prot->id);
      SeqIdWrite (sip, str, PRINTID_REPORT, sizeof (str));
      Message (MSG_POSTERR, "Suggest failure for %s", str);
      return NULL;
    }
  }

  w = MovableModalWindow (-50, -33, -10, -10, "Enter coding region interval", NULL);
  g = HiddenGroup (w, -1, 0, NULL);
  m = NULL;
  SetGroupSpacing (g, 3, 10);
  if (prot->descr != NULL) {
    vnp = ValNodeFindNext (prot->descr, NULL, Seq_descr_title);
    if (vnp != NULL && vnp->data.ptrvalue != NULL) {
      m = MultiLinePrompt (g, (CharPtr) vnp->data.ptrvalue, stdCharWidth * 28, programFont);
    }
  }
  d = CreateIntervalEditorDialog (g, NULL, 4, 2, sep, TRUE, FALSE);
  c = HiddenGroup (g, 2, 0, NULL);
  SetGroupSpacing (c, 10, 2);
  PushButton (c, "Accept", AcceptAskProc);
  PushButton (c, "Cancel", CancelAskProc);
  AlignObjects (ALIGN_CENTER, (HANDLE) d, (HANDLE) c, (HANDLE) m, NULL);
  Show (w);
  Select (w);
  intBoxUp = TRUE;
  intBoxRsult = FALSE;
  while (intBoxUp) {
    ProcessEventOrIdle ();
  }
  ProcessAnEvent ();
  if (intBoxRsult) {
    slp = (SeqLocPtr) DialogToPointer (d);
  }
  Remove (w);
  return slp;
}

extern Boolean AutomaticProteinProcess (SeqEntryPtr esep, SeqEntryPtr psep,
                                        Int2 code, Boolean makeMRNA)

{
  SeqFeatPtr   cds;
  CdRegionPtr  crp;
  Char         mRnaName [128];
  BioseqPtr    nbsp;
  SeqEntryPtr  nsep;
  BioseqPtr    pbsp;
  SeqFeatPtr   rna;
  RnaRefPtr    rrp;
  SeqLocPtr    slp;
  CharPtr      ttl;
  ValNodePtr   vnp;
  CharPtr      vnpstr;

  if (esep == NULL || psep == NULL) return FALSE;

  nsep = FindNucSeqEntry (esep);
  if (nsep == NULL || (! IS_Bioseq (nsep)) || (! IS_Bioseq (psep))) return FALSE;

  nbsp = (BioseqPtr) nsep->data.ptrvalue;
  pbsp = (BioseqPtr) psep->data.ptrvalue;
  if (nbsp == NULL || pbsp == NULL) return FALSE;

  cds = NULL;
  WatchCursor ();
  Update ();
  slp = PredictCodingRegion (nbsp, pbsp, code);
  if (slp == NULL) {
    ArrowCursor ();
    Update ();
    slp = AskForInterval (nsep, nbsp, pbsp);
  }
  if (slp == NULL) return FALSE;

  mRnaName [0] = '\0';
  ttl = NULL;
  crp = CreateNewCdRgn (0, FALSE, code);
  if (crp != NULL) {
    if (pbsp->descr != NULL) {
      vnp = ValNodeFindNext (pbsp->descr, NULL, Seq_descr_title);
      if (vnp != NULL && vnp->data.ptrvalue != NULL) {
        vnpstr = (CharPtr) vnp->data.ptrvalue;
        CreateGeneAndProtFeats (nsep, psep, slp, crp, vnpstr, mRnaName, sizeof (mRnaName), &ttl);
        ExciseString (vnpstr, "[gene=", "]");
        ExciseString (vnpstr, "[prot=", "]");
        ExciseString (vnpstr, "[function=", "]");
        ExciseString (vnpstr, "[EC_number=", "]");
        ExciseString (vnpstr, "[orf", "]");
        ExciseString (vnpstr, "[comment", "]");
        ExciseString (vnpstr, "[note", "]");
        TrimSpacesAroundString (vnpstr);
        if (StringHasNoText (vnpstr)) {
          ValNodeExtract (&(pbsp->descr), Seq_descr_title);
        }
      }
    }
    if (makeMRNA) {
      rrp = RnaRefNew ();
      if (rrp != NULL) {
        rrp->type = 2;
        if (! StringHasNoText (mRnaName)) {
          rrp->ext.choice = 1;
          rrp->ext.value.ptrvalue = StringSave (mRnaName);
        }
        rna = CreateNewFeature (nsep, NULL, SEQFEAT_RNA, NULL);
        if (rna != NULL) {
          rna->data.value.ptrvalue = (Pointer) rrp;
          rna->location = SeqLocFree (rna->location);
          rna->location = AsnIoMemCopy ((Pointer) slp,
                                        (AsnReadFunc) SeqLocAsnRead,
                                        (AsnWriteFunc) SeqLocAsnWrite);
        }
      }
    }
    cds = CreateNewFeature (nsep, NULL, SEQFEAT_CDREGION, NULL);
    if (cds != NULL) {
      cds->data.value.ptrvalue = (Pointer) crp;
      cds->location = SeqLocFree (cds->location);
      cds->location = slp;
      slp = NULL;
      SetSeqFeatProduct (cds, pbsp);
      if (! StringHasNoText (ttl)) {
        cds->comment = ttl;
      }
    }
  }

  SeqLocFree (slp);
  return TRUE;
}

typedef struct fa2htgsform {
  FORM_MESSAGE_BLOCK

  SeqSubmitPtr       ssp;
  SeqEntryPtr        sep;

  GrouP              templateblock;
  GrouP              fastablock;
  GrouP              orderblock;
  GrouP              controlblock;
  GrouP              contigtype;

  DialoG             contigorder;
  EnumFieldAssocPtr  alists [1];
  GrouP              htgsphase;
  TexT               orgname;
  TexT               seqname;
  ButtoN             update;
  TexT               accession;
  TexT               knownlength;
  TexT               remark;
  TexT               clone;
  TexT               chromosome;
  TexT               title;
  DialoG             secondaries;

  SeqEntryPtr        seplist;

  ButtoN             okBtn;
  BtnActnProc        finish;
  BtnActnProc        cancel;
  Boolean            readPhrap;
  Boolean            buildContig;

} Fa2htgsForm, PNTR Fa2htgsFormPtr;

/*------------- MakeAc2GBSeqId() -----------------------*/
/***************************************************************
*   MakeAc2GBSeqId:
*   -- return NULL if acnum == null
*                                             Hsiu-Chuan 4-18-97
****************************************************************/
static SeqIdPtr  SqnMakeAc2GBSeqId(CharPtr accession)
{
   TextSeqIdPtr tsip;
   SeqIdPtr sip;

   if (accession == NULL || *accession == '\0')
      return NULL;

   sip = ValNodeNew(NULL);
   sip->choice = SEQID_GENBANK;
   tsip = TextSeqIdNew();
   sip->data.ptrvalue = tsip;
   tsip->accession = StringSave(accession);

   return sip;

} /* MakeAc2GBSeqId */

/*----------- AddExtraAc2Entry() ----------------------------*/
/***************************************************************
*   AddExtraAc2Entry:
*                                             Hsiu-Chuan 4-11-97, modified by JK
****************************************************************/
static Boolean SqnAddExtraAc2Entry (SeqEntryPtr entry , ValNodePtr extra_accs )
{
   BioseqPtr  bsp;
   ValNodePtr vnp;
   GBBlockPtr gbp;
   Char       acnum[17];
   CharPtr    p;
   Int4       i, j;
   SeqHistPtr shp;
   SeqIdPtr   sip;
   ValNodePtr tmp;

   if ((entry == NULL) || (extra_accs == NULL))
      return FALSE;

   bsp = (BioseqPtr)(entry->data.ptrvalue);

   for (gbp= NULL, vnp = bsp->descr; vnp != NULL; vnp = vnp->next)
   {
       if (vnp->choice == Seq_descr_genbank)
       {
          gbp = vnp->data.ptrvalue;
          break;
       }
   }

   shp = bsp->hist; 

   if (gbp == NULL)
   {
      vnp = (ValNodePtr) NewDescrOnSeqEntry (entry, Seq_descr_genbank);
      gbp = GBBlockNew();
      vnp->data.ptrvalue = (Pointer)gbp;
   }
   
   for (tmp = extra_accs; tmp != NULL; tmp = tmp->next)
   {
       p = (CharPtr) tmp->data.ptrvalue;
       if (p == NULL) continue;
       for (i = 0; isalnum(*p) && *p != '\0'; ++p, ++i)
           acnum[i] = *p;
       acnum[i] = '\0'; 
               /* check one_letter+5digits or two_letter+6digits */
       if (i == 6 || i == 8)
       {
          if (!isalpha(acnum[0]) || (!(isdigit(acnum[1]) && i == 6) &&
              !(isalpha(acnum[1]) && i == 8)))
          {
             ErrPostEx(SEV_ERROR,0,0,
 "Invalid accession (one_letter+5digits or two_letter+6digits): %s",
                                                           acnum);
             return FALSE;
          }

          for (j = 2; j < i; ++j)
          {
              if (!(isdigit(acnum[j])))
              {
                 ErrPostEx(SEV_ERROR,0,0,
 "Invalid accession (one_letter+5digits or two_letter+6digits): %s",
                                                           acnum);
                 return FALSE;
              }
          }

          ValNodeCopyStr(&gbp->extra_accessions, 0, acnum);
          sip = SqnMakeAc2GBSeqId (acnum);
          if (shp == NULL)
          {
             shp = SeqHistNew();
             bsp->hist = shp;
          }
          ValNodeLink(&shp->replace_ids, sip);
       }
       else
       {
          ErrPostEx(SEV_ERROR,0,0,
 "Invalid accession (one_letter+5digits or two_letter+6digits): %s",
                                                           acnum);
          return FALSE;
       }

       while (!isalnum(*p) && *p != '\0')
           ++p;
   }

   return TRUE;

} /* AddExtraAc2Entry */

static void RescueSeqGraphs (BioseqPtr bsp, Int2 index, ValNodePtr PNTR vnpp)

{
  SeqAnnotPtr   nextsap;
  SeqGraphPtr   nextsgp;
  Pointer PNTR  prevsap;
  Pointer PNTR  prevsgp;
  SeqAnnotPtr   sap;
  SeqGraphPtr   sgp;

  if (bsp == NULL || vnpp == NULL) return;
  sap = bsp->annot;
  prevsap = (Pointer PNTR) &(bsp->annot);
  while (sap != NULL) {
    nextsap = sap->next;
    if (sap->type == 3) {
      sgp = (SeqGraphPtr) sap->data;
      prevsgp = (Pointer PNTR) &(sap->data);
      while (sgp != NULL) {
        nextsgp = sgp->next;
        *(prevsgp) = sgp->next;
        sgp->next = NULL;
        ValNodeAddPointer (vnpp, index, (Pointer) sgp);
        sgp = nextsgp;
      }
    }
    if (sap->data == NULL) {
      *(prevsap) = sap->next;
      sap->next = NULL;
      SeqAnnotFree (sap);
    } else {
      prevsap = (Pointer PNTR) &(sap->next);
    }
    sap = nextsap;
  }
}

static SeqAnnotPtr NewSeqAnnotType3 (CharPtr name, SeqGraphPtr sgp)

{
  SeqAnnotPtr  sap = NULL;

  if (sgp == NULL) return NULL;
  sap = SeqAnnotNew ();
  if (sap == NULL) return NULL;

  if (! StringHasNoText (name)) {
    ValNodeAddPointer (&(sap->desc), Annot_descr_name, StringSave (name));
  }
  sap->type = 3;
  sap->data = (Pointer) sgp;

  return sap;
}

static void OffsetAndLinkSeqGraph (BioseqPtr bsp, SeqGraphPtr sgp, Int2 index)

{
  DeltaSeqPtr  dsp;
  SeqGraphPtr  lastsgp;
  Int4         len;
  SeqLitPtr    litp;
  SeqAnnotPtr  sap;
  SeqIntPtr    sintp;
  SeqLocPtr    slp;

  if (bsp == NULL || sgp == NULL || index < 1) return;
  len = 0;
  if (bsp->repr == Seq_repr_delta && bsp->seq_ext_type == 4) {
    for (dsp = (DeltaSeqPtr) (bsp->seq_ext);
         dsp != NULL && index > 1; dsp = dsp->next, index--) {
      if (dsp->choice == 1) {
        len += SeqLocLen ((SeqLocPtr) dsp->data.ptrvalue);
      } else if (dsp->choice == 2) {
        litp = (SeqLitPtr) dsp->data.ptrvalue;
        if (litp != NULL) {
          len += litp->length;
        }
      }
    }
  }
  slp = sgp->loc;
  if (slp != NULL && slp->choice == SEQLOC_INT) {
    sintp = (SeqIntPtr) slp->data.ptrvalue;
    if (sintp != NULL) {
      sintp->from += len;
      sintp->to += len;
      sintp->id = SeqIdFree (sintp->id);
      sintp->id = SeqIdDup (bsp->id);
    }
  }
  for (sap = bsp->annot; sap != NULL; sap = sap->next) {
    if (sap->type == 3) {
      for (lastsgp = sap->data; lastsgp->next != NULL; lastsgp = lastsgp->next) {
        continue;
      }
      lastsgp->next = sgp;
      break;
    }
  }
  if (sap == NULL) {
    if (bsp->annot != NULL) {
      for (sap = bsp->annot; sap->next != NULL; sap = sap->next) {
        continue;
      }
      sap->next = NewSeqAnnotType3 ("Graphs", sgp);
    } else {
      bsp->annot = NewSeqAnnotType3 ("Graphs", sgp);
    }
  }
}

static CharPtr phrapBoilerPlate = "Sequence Quality Assessment:~ \
This entry has been annotated with sequence quality~ \
estimates computed by the Phrap assembly program.~ \
All manually edited bases have been reduced to quality zero.~ \
Quality levels above 40 are expected to have less than ~ \
1 error in 10,000 bp.~ \
Base-by-base quality values are not generally visible from the~ \
GenBank flat file format but are available as part~ \
of this entry's ASN.1 file.~----------------------~";

static void ProcessFa2htgs (Fa2htgsFormPtr ffp, SeqSubmitPtr ssp)

{
  SeqEntryPtr  sep, oldsep, the_entry, nextsep;
  NCBISubPtr nsp;
  Uint1 htgs_phase;
  CharPtr seqname = NULL, accession = NULL, orgname = NULL;
  CharPtr clone = NULL, chromosome = NULL;
  CharPtr remark = NULL, title = NULL, seqbuf = NULL;
  Int4 length = 0, cumlength = 0;
  BioseqPtr bsp;
  BioseqSetPtr bssp;
  SeqLitPtr slp;
  ValNodePtr vnp, PNTR prevpnt, next, extra_accs;
  Boolean lastwasraw;
  Char str [32];
  long int val;
  Int2 index = 0;
  ValNodePtr rescuedsgps = NULL;
  ObjectIdPtr    oip;
  UserFieldPtr   ufp;
  UserObjectPtr  uop;
  DatePtr dp;

  if (ffp == NULL || ssp == NULL) return;

  htgs_phase = (Uint1) GetValue (ffp->htgsphase);
  orgname = SaveStringFromText (ffp->orgname);
  seqname = SaveStringFromText (ffp->seqname);
  if (GetStatus (ffp->update)) {
    accession = SaveStringFromText (ffp->accession);
  }
  clone = SaveStringFromText (ffp->clone);
  chromosome = SaveStringFromText (ffp->chromosome);
  remark = SaveStringFromText (ffp->remark);
  title = SaveStringFromText (ffp->title);
  extra_accs = DialogToPointer (ffp->secondaries);

  length = 0;
/* may need to really calculate length */
  GetTitle (ffp->knownlength, str, sizeof (str));
  if (! StringHasNoText (str)) {
    if (sscanf (str, "%ld", &val) == 1 && val > 0) {
      length = (Int4) val;
    }
  }

/* modified from fa2htgs */
   oldsep = (SeqEntryPtr)(ssp->data);  /* clear out template */
   ssp->data = NULL;
   MemFree(ssp->sub->tool);
   ssp->sub->tool = StringSave(SEQUIN_APP_VERSION);
   nsp = MemNew(sizeof(NCBISub));
   nsp->ssp = ssp;
   nsp->submittor_key = StringSave (genomeCenter);
   /*
   MemFree(ssp->sub->cit->descr);
   ssp->sub->cit->descr = remark;
   */

   cumlength = 0;
   index = 0;

   if (ffp->buildContig) {
     ssp->data = (Pointer) ffp->seplist;
     the_entry = ffp->seplist;
     sep = the_entry;

     oip = ObjectIdNew ();
     oip->str = StringSave ("info");
     uop = UserObjectNew ();
     uop->type = oip;
     uop->_class = StringSave ("Genomes");

     oip = ObjectIdNew ();
     oip->id = 0;
     ufp = UserFieldNew ();
     ufp->choice = 2;
     ufp->data.intvalue = 0;
     ufp->label = oip;

     uop->data = ufp;

     if (IS_Bioseq (sep)) {
       bsp = (BioseqPtr) sep->data.ptrvalue;
       vnp = ValNodeNew (NULL);
       vnp->choice = Seq_descr_user;
       vnp->data.ptrvalue = (Pointer) uop;
       vnp->next = bsp->descr;
       bsp->descr = vnp;
       cumlength = bsp->length;
     }

   }
   else if (htgs_phase < 3)
   {
      the_entry = AddDeltaSeqOnlyToSubmission (
                        nsp,
                        seqname,
                        NULL,
                        accession,
                        0,
                        MOLECULE_CLASS_DNA,
                        MOLECULE_TYPE_GENOMIC,
                        length,
                        TOPOLOGY_LINEAR,
                        STRANDEDNESS_DOUBLE);

      sep = ffp->seplist;
      lastwasraw = FALSE;
      while (sep != NULL)
      {
         nextsep = sep->next;
         sep->next = NULL;
         bsp = (BioseqPtr)(sep->data.ptrvalue);
         if (bsp->repr == Seq_repr_raw)
         {
            if (lastwasraw) {
               AddGapToDeltaSeq(nsp, the_entry, 0);
               index++;
            }
            BioseqRawConvert(bsp, Seq_code_iupacna);
            seqbuf = BSMerge((ByteStorePtr)(bsp->seq_data), NULL);
            slp = AddLiteralToDeltaSeq(nsp, the_entry,
               bsp->length);
            AddBasesToLiteral(nsp, slp, seqbuf);
            MemFree(seqbuf);
            lastwasraw = TRUE;
            index++;
         }
         else
         {
            if (bsp->length < 0)
               bsp->length = 0;  /* -1 may be set */
            AddGapToDeltaSeq(nsp, the_entry,
               bsp->length);
            lastwasraw = FALSE;
            index++;
         }
         cumlength += bsp->length;
         RescueSeqGraphs (bsp, index, &rescuedsgps);
         SeqEntryFree(sep);
         sep = nextsep;
      }
   }
   else
   {
    the_entry = AddSeqOnlyToSubmission (
                        nsp,
                        seqname,
                        NULL,
                        accession,
                        0,
                        MOLECULE_CLASS_DNA,
                        MOLECULE_TYPE_GENOMIC,
                        length,
                        TOPOLOGY_LINEAR,
                        STRANDEDNESS_DOUBLE);

      sep = ffp->seplist;
      nextsep = sep->next;
      sep->next = NULL;
      bsp = (BioseqPtr)(sep->data.ptrvalue);
      if (bsp->repr == Seq_repr_raw)
      {
         BioseqRawConvert(bsp, Seq_code_iupacna);
         seqbuf = BSMerge((ByteStorePtr)(bsp->seq_data), NULL);
         AddBasesToBioseq(nsp, the_entry, seqbuf);
         MemFree(seqbuf);
         index++;
      }
      cumlength += bsp->length;
      RescueSeqGraphs (bsp, index, &rescuedsgps);
      SeqEntryFree(sep);
      if (nextsep != NULL) {
        ErrPostEx (SEV_ERROR ,0, 0, "Only the first contig was used for HTGS 3");
      }
      if (length > 0 && length != cumlength) {
        ErrPostEx (SEV_ERROR ,0, 0, "Length is exactly %ld, not %ld",
                   (long) cumlength, (long) length);
        length = cumlength;
      }
   }

    /* get data from template: pub, organism, and comment */
   if (IS_Bioseq(oldsep))
   {
      bsp = (BioseqPtr)(oldsep->data.ptrvalue);
      prevpnt = &(bsp->descr);
   }
   else
   {
      bssp = (BioseqSetPtr)(oldsep->data.ptrvalue);
      prevpnt = &(bssp->descr);
   }

   bsp = (BioseqPtr)(the_entry->data.ptrvalue);
   if (bsp != NULL) {
     bsp->length = MAX (cumlength, length);
   }

   for (vnp = *prevpnt; vnp != NULL; vnp = next)
   {
      next = vnp->next;
      if (vnp->choice == Seq_descr_pub)
      {
         *prevpnt = next;
         vnp->next = NULL;
         ValNodeLink(&(bsp->descr), vnp);
      }
      else
         prevpnt = &(vnp->next);
   }
   if (remark != NULL) {
     vnp = ValNodeNew (NULL);
     if (vnp != NULL) {
       vnp->choice = Seq_descr_comment;
       vnp->data.ptrvalue = remark;
       ValNodeLink(&(bsp->descr), vnp);
     }
   }

   SeqEntryFree(oldsep);

   AddOrganismToEntryNew(nsp, the_entry, orgname, NULL, NULL, NULL,
                         NULL, NULL, NULL, NULL);

   AddGenomeToEntry(nsp, the_entry, 1);
   if (clone != NULL)
      AddSubSourceToEntry(nsp, the_entry, 3, clone);
   if (chromosome != NULL)
       AddSubSourceToEntry(nsp, the_entry, 1, chromosome);
   if (title != NULL)
      AddTitleToEntry(nsp, the_entry, title);
   if (ffp->readPhrap) {
      AddCommentToEntry(nsp, the_entry, phrapBoilerPlate);
   }

   if (extra_accs != NULL) {
      SqnAddExtraAc2Entry(the_entry, extra_accs);
   }

   AddBiomolToEntry(nsp, the_entry, 1);
   if (ffp->buildContig) {
   } else {
     AddTechToEntry(nsp, the_entry, (Uint1)(MI_TECH_htgs_1 + htgs_phase - 1));
   }
   vnp = NewDescrOnSeqEntry (the_entry, Seq_descr_create_date);
   if (vnp != NULL) {
     dp = DateCurr ();
     vnp->data.ptrvalue = (Pointer) dp;
   }

   if (bsp != NULL) {
     for (vnp = rescuedsgps; vnp != NULL; vnp = vnp->next) {
       OffsetAndLinkSeqGraph (bsp, (SeqGraphPtr) vnp->data.ptrvalue, (Int2) vnp->choice);
       vnp->data.ptrvalue = NULL;
     }
   }
   rescuedsgps = ValNodeFreeData (rescuedsgps);

   MemFree (nsp);

}

static Pointer Fa2htgsToSeqSubmitPtr (ForM f)

{
  Fa2htgsFormPtr  ffp;
  SeqSubmitPtr    ssp;

  ffp = (Fa2htgsFormPtr) GetObjectExtra (f);
  if (ffp == NULL) return NULL;
  ProcessFa2htgs (ffp, ffp->ssp);
  ssp = ffp->ssp;
  ffp->ssp = NULL;
  ffp->sep = NULL;
  ffp->seplist = NULL;
  return (Pointer) ssp;
}

static void Fa2htgsFormMessage (ForM f, Int2 mssg)

{
  Fa2htgsFormPtr  ffp;

  ffp = (Fa2htgsFormPtr) GetObjectExtra (f);
  if (ffp) {
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

#ifndef WIN_MAC
void CreateFa2htgsFormMenus (WindoW w)

{
  BaseFormPtr  bfp;
  MenU         m;

  bfp = (BaseFormPtr) GetObjectExtra (w);
  if (bfp != NULL) {
    m = PulldownMenu (w, "File");
    AddAboutAndHelpMenuItems (m);
    FormCommandItem (m, "Quit", bfp, VIB_MSG_QUIT);
    m = PulldownMenu (w, "Edit");
    FormCommandItem (m, CUT_MENU_ITEM, bfp, VIB_MSG_CUT);
    FormCommandItem (m, COPY_MENU_ITEM, bfp, VIB_MSG_COPY);
    FormCommandItem (m, PASTE_MENU_ITEM, bfp, VIB_MSG_PASTE);
    FormCommandItem (m, CLEAR_MENU_ITEM, bfp, VIB_MSG_DELETE);
  }
}
#endif

static void FillInFa2htgsFields (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr       bsp;
  BioseqSetPtr    bssp;
  Fa2htgsFormPtr  ffp;
  ValNodePtr      sdp = NULL;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  ffp = (Fa2htgsFormPtr) mydata;
  if (ffp == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sdp = bsp->descr;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sdp = bssp->descr;
  } else return;
  while (sdp != NULL) {
    if (sdp->choice == Seq_descr_comment) {
      SetTitle (ffp->remark, (CharPtr) sdp->data.ptrvalue);
    }
    sdp = sdp->next;
  }
}

static void ReadTemplate (ButtoN b)

{
  AsnIoPtr        aip;
  Fa2htgsFormPtr  ffp;
  Char            path [PATH_MAX];
  SeqEntryPtr     sep;
  SeqSubmitPtr    ssp;
  Char            str [128];

  ffp = (Fa2htgsFormPtr) GetObjectExtra (b);
  if (ffp == NULL) return;
  if (GetInputFileName (path, sizeof (path), "", "TEXT")) {
    ffp->ssp = SeqSubmitFree (ffp->ssp);
    aip = AsnIoOpen (path, "r");
    if (aip != NULL) {
      ffp->ssp = SeqSubmitAsnRead (aip, NULL);
    }
    AsnIoClose (aip);
  }
  if (ffp->ssp != NULL) {
    SafeShow (ffp->fastablock);
    SafeDisable (b);
    ssp = ffp->ssp;
    if (ssp->datatype == 1) {
      sep = (SeqEntryPtr) ssp->data;
      if (sep != NULL) {
        SeqEntryToGeneticCode (sep, NULL, str, sizeof (str));
        SetTitle (ffp->orgname, str);
        SeqEntryExplore (sep, (Pointer) ffp, FillInFa2htgsFields);
      }
    }
  }
}

static Boolean Fa2htgsFormOkay (Fa2htgsFormPtr ffp)

{
  if (ffp == NULL) return FALSE;
  if (ffp->ssp == NULL) return FALSE;
  if (ffp->seplist == NULL) return FALSE;
  if (GetValue (ffp->htgsphase) == 0 && (! ffp->buildContig)) return FALSE;
  if (TextHasNoText (ffp->orgname)) return FALSE;
  if (TextHasNoText (ffp->seqname)) return FALSE;
  return TRUE;
}

static void ReadFastaHtgsFile (ButtoN b)

{
  Fa2htgsFormPtr  ffp;
  FILE            *fp;
  SeqEntryPtr     head;
  Char            path [PATH_MAX];
  SeqEntryPtr     sep;
  CharPtr         ttl;

  ffp = (Fa2htgsFormPtr) GetObjectExtra (b);
  if (ffp == NULL) return;
  if (! GetInputFileName (path, sizeof (path), "", "TEXT")) return;
  if (! StringHasNoText (path)) {
    fp = FileOpen (path, "r");
    if (fp != NULL) {
      head = NULL;
      sep = FastaToSeqEntry (fp, TRUE);
      while (sep != NULL) {
        ValNodeLink (&head, sep);
        sep = FastaToSeqEntry (fp, TRUE);
      }
      if (head != NULL) {
        ttl = SeqEntryGetTitle (head);
        SetTitle (ffp->title, ttl);
      }
      ffp->seplist = head;
    }
    FileClose (fp);
    if (ffp->seplist != NULL) {
      SafeShow (ffp->controlblock);
      SafeDisable (b);
      if (TextHasNoText (ffp->orgname)) {
        Select (ffp->orgname);
      } else {
        Select (ffp->seqname);
      }
      if (Fa2htgsFormOkay (ffp)) {
        SafeEnable (ffp->okBtn);
      } else {
        SafeDisable (ffp->okBtn);
      }
    }
  }
}

static EnumFieldAssocPtr MakePhrapAlists (SeqEntryPtr head)

{
  EnumFieldAssocPtr  alist = NULL;
  BioseqPtr          bsp;
  Int2               count;
  Int2               j;
  SeqEntryPtr        sep;
  SeqIdPtr           sip;
  Char               str [128];

  if (head == NULL) return NULL;
  count = ValNodeLen (head);
  alist = MemNew (sizeof (EnumFieldAssoc) * (size_t) (count + 4));
  if (alist == NULL) return NULL;

  j = 0;
  alist [j].name = StringSave ("                              ");
  alist [j].value = (UIEnum) 0;
  for (j = 1, sep = head; j <= count && sep != NULL; j++, sep = sep->next) {
    if (sep != NULL && sep->choice == 1 && sep->data.ptrvalue != NULL) {
      bsp = (BioseqPtr) sep->data.ptrvalue;
      sip = SeqIdFindWorst (bsp->id);
      SeqIdWrite (sip, str, PRINTID_REPORT, sizeof (str));
      str [30] = '\0';
      alist [j].name = StringSave (str);
      alist [j].value = (UIEnum) j;
    }
  }
  j = count + 1;
  alist [j].name = NULL;
  alist [j].value = (UIEnum) 0;


  return alist;
}

static void ReadAPhrapFile (ButtoN b)

{
  PopuP           control;
  Int2            count;
  Fa2htgsFormPtr  ffp;
  FILE            *fp;
  SeqEntryPtr     head;
  Int2            i;
  Char            path [PATH_MAX];
  RecT            r;
  TagListPtr      tlp;

  ffp = (Fa2htgsFormPtr) GetObjectExtra (b);
  if (ffp == NULL) return;
  if (! GetInputFileName (path, sizeof (path), "", "TEXT")) return;
  if (! StringHasNoText (path)) {
    WatchCursor ();
    fp = FileOpen (path, "r");
    if (fp != NULL) {
      head = ReadPhrapFile (fp);
      ffp->seplist = head;
    }
    FileClose (fp);
    tlp = (TagListPtr) GetObjectExtra (ffp->contigorder);
    if (tlp != NULL) {
      ffp->alists [0] = MakePhrapAlists (ffp->seplist);
      tlp->alists = ffp->alists;
      for (i = 0; i < tlp->rows; i++) {
        control = (PopuP) tlp->control [i * MAX_TAGLIST_COLS + 0];
        if (control != NULL) {
          GetPosition (control, &r);
          Reset (control);
          InitEnumPopup (control, ffp->alists [0], NULL);
          SetEnumPopup (control, ffp->alists [0], 0);
          SetPosition (control, &r);
        }
      }
      count = ValNodeLen (ffp->seplist);
      for (i = 0; i < count; i++) {
        ValNodeCopyStr (&(tlp->vnp), 0, "0");
      }
      tlp->max = MAX ((Int2) 0, (Int2) (count - tlp->rows));
      CorrectBarMax (tlp->bar, tlp->max);
      CorrectBarPage (tlp->bar, tlp->rows - 1, tlp->rows - 1);
    }
    SafeShow (ffp->orderblock);
    SafeDisable (b);
    ArrowCursor ();
  }
}

static void PhrapOrderChosen (ButtoN b)

{
  SeqEntryPtr     PNTR collision;
  Int2            count;
  Fa2htgsFormPtr  ffp;
  Int2            i;
  Int2            j;
  SeqEntryPtr     lastsep;
  SeqEntryPtr     nextsep;
  Boolean         okay;
  SeqEntryPtr     PNTR order;
  SeqEntryPtr     sep;
  CharPtr         str;
  TagListPtr      tlp;
  int             val;
  ValNodePtr      vnp;
  /*
  Char            contigs [256];
  Char            tmp [256];
  BioseqPtr       bsp;
  SeqIdPtr        sip;
  Boolean         notfirst;
  */

  ffp = (Fa2htgsFormPtr) GetObjectExtra (b);
  if (ffp == NULL) return;
  tlp = (TagListPtr) GetObjectExtra (ffp->contigorder);
  if (tlp == NULL) return;
  count = ValNodeLen (tlp->vnp);
  order = MemNew (sizeof (SeqEntryPtr) * (count + 1));
  if (order == NULL) return;
  collision = MemNew (sizeof (SeqEntryPtr) * (count + 1));
  if (collision == NULL) return;

  okay = TRUE;
  for (i = 0, vnp = tlp->vnp; i < count && vnp != NULL; i++, vnp = vnp->next) {
    str = ExtractTagListColumn ((CharPtr) vnp->data.ptrvalue, 0);
    if (str != NULL && sscanf (str, "%d", &val) == 1 && val > 0) {
      val--;
      for (j = 0, sep = ffp->seplist; j < (Int2) val && sep != NULL; j++, sep = sep->next) continue;
      if (sep != NULL) {
        if (collision [j] != NULL) {
          okay = FALSE;
        }
        collision [j] = sep;
        order [i] = sep;
      }
    }
    MemFree (str);
  }
  if (! okay) {
    Message (MSG_ERROR, "You must not select a contig more than once");
    MemFree (order);
    MemFree (collision);
    return;
  }

  okay = FALSE;
  for (i = 0; i < count; i++) {
    if (order [i] != NULL) {
      okay = TRUE;
    }
  }
  if (! okay) {
    Message (MSG_ERROR, "You must select at least one contig");
    MemFree (order);
    MemFree (collision);
    return;
  }

  /* use spreadsheet to reorder ffp->seplist, delete unwanted items */

  /*
  contigs [0] = '\0';
  notfirst = FALSE;
  for (i = 0; i < count; i++) {
    if (order [i] != NULL) {
      sep = order [i];
      bsp = (BioseqPtr) sep->data.ptrvalue;
      sip = SeqIdFindWorst (bsp->id);
      SeqIdWrite (sip, tmp, PRINTID_REPORT, sizeof (tmp));
      if (notfirst) {
        StringCat (contigs, " ");
      }
      StringCat (contigs, tmp);
      notfirst = TRUE;
    }
  }
  ffp->seplist = SetPhrapContigOrder (ffp->seplist, contigs);
  */

  for (i = 0; i < count; i++) {
    if (order [i] != NULL) {
      sep = ffp->seplist;
      lastsep = NULL;
      while (sep != NULL && sep != order [i]) {
        lastsep = sep;
        sep = sep->next;
      }
      if (sep != NULL) {
        if (lastsep != NULL) {
          lastsep->next = sep->next;
          sep->next = NULL;
        } else {
          ffp->seplist = sep->next;
          sep->next = NULL;
        }
      }
    }
  }
  sep = ffp->seplist;
  while (sep != NULL) {
    nextsep = sep->next;
    sep->next = NULL;
    SeqEntryFree (sep);
    sep = nextsep;
  }
  ffp->seplist = NULL;

  for (i = 0; i < count; i++) {
    if (order [i] != NULL) {
      ValNodeLink (&(ffp->seplist), order [i]);
    }
  }

  MemFree (order);
  MemFree (collision);
  SafeDisable (b);
  SafeHide (ffp->orderblock);
  SafeShow (ffp->controlblock);
  if (TextHasNoText (ffp->orgname)) {
    Select (ffp->orgname);
  } else {
    Select (ffp->seqname);
  }
  if (Fa2htgsFormOkay (ffp)) {
    SafeEnable (ffp->okBtn);
  } else {
    SafeDisable (ffp->okBtn);
  }
}

static void ReadContigFile (ButtoN b)

{
  BioseqPtr       bsp;
  Fa2htgsFormPtr  ffp;
  FILE            *fp;
  Boolean         onMaster;
  Char            path [PATH_MAX];
  SeqEntryPtr     sep = NULL;

  ffp = (Fa2htgsFormPtr) GetObjectExtra (b);
  if (ffp == NULL) return;
  if (! GetInputFileName (path, sizeof (path), "", "TEXT")) return;
  if (! StringHasNoText (path)) {
    WatchCursor ();
    fp = FileOpen (path, "r");
    if (fp != NULL) {
      onMaster = (Boolean) (GetValue (ffp->contigtype) == 2);
      sep = ReadContigList (fp, onMaster);
      if (sep != NULL && sep->choice == 1) {
        ffp->seplist = sep;
        bsp = (BioseqPtr) sep->data.ptrvalue;
        SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, sep);
      }
    }
    FileClose (fp);
    ArrowCursor ();
    if (ffp->seplist != NULL) {
      SafeShow (ffp->controlblock);
      SafeDisable (b);
      SafeDisable (ffp->contigtype);
      if (TextHasNoText (ffp->orgname)) {
        Select (ffp->orgname);
      } else {
        Select (ffp->seqname);
      }
      if (Fa2htgsFormOkay (ffp)) {
        SafeEnable (ffp->okBtn);
      } else {
        SafeDisable (ffp->okBtn);
      }
    }
  }
}

static void AcceptFa2htgs (ButtoN b)

{
  Fa2htgsFormPtr  ffp;

  ffp = (Fa2htgsFormPtr) GetObjectExtra (b);
  if (ffp == NULL) return;
  Hide (ffp->form);
  Update ();
  if (ffp->finish != NULL) {
    ffp->finish (b);
  }
  SeqSubmitFree (ffp->ssp);
  SeqEntryFree (ffp->sep);
  Remove (ffp->form);
}

static void CancelFa2htgs (ButtoN b)

{
  Fa2htgsFormPtr  ffp;

  ffp = (Fa2htgsFormPtr) GetObjectExtra (b);
  if (ffp == NULL) return;
  Hide (ffp->form);
  Update ();
  if (ffp->cancel != NULL) {
    ffp->cancel (b);
  }
  SeqSubmitFree (ffp->ssp);
  SeqEntryFree (ffp->sep);
  Remove (ffp->form);
}

static void SetFa2htgsAcceptBtn (Handle control)

{
  Fa2htgsFormPtr  ffp;

  ffp = (Fa2htgsFormPtr) GetObjectExtra (control);
  if (ffp == NULL) return;
  if (Fa2htgsFormOkay (ffp)) {
    SafeEnable (ffp->okBtn);
  } else {
    SafeDisable (ffp->okBtn);
  }
}

static void SetFa2htgsUpdate (ButtoN b)

{
  Fa2htgsFormPtr  ffp;

  ffp = (Fa2htgsFormPtr) GetObjectExtra (b);
  if (ffp == NULL) return;
  if (GetStatus (b)) {
    SafeEnable (ffp->accession);
  } else {
    SafeDisable (ffp->accession);
  }
  SetFa2htgsAcceptBtn ((Handle) b);
}

static void CleanupGenomeCenterForm (GraphiC g, VoidPtr data)

{
  Fa2htgsFormPtr  ffp;
  SeqEntryPtr     next;
  SeqEntryPtr     sep;

  ffp = (Fa2htgsFormPtr) data;
  if (ffp != NULL) {
    sep = ffp->seplist;
    while (sep != NULL) {
      next = sep->next;
      sep->next = NULL;
      SeqEntryFree (sep);
      sep = next;
    }
    if (ffp->alists [0] != NULL) {
      FreeEnumFieldAlist (ffp->alists [0]);
    }
  }
  StdCleanupFormProc (g, data);
}

static CharPtr secaccstrings [] = {"Secondary", "Accessions", NULL};

static Uint2 contigorder_types [] = {
  TAGLIST_POPUP
};

static ENUM_ALIST(contigdefault_alist)
  {"                              ", 0},
END_ENUM_ALIST

static EnumFieldAssocPtr contigdefault_alists [] = {
  contigdefault_alist
};

extern DialoG CreateTagListDialogEx (GrouP h, Uint2 rows, Uint2 cols,
                                     Int2 spacing, Uint2Ptr types,
                                     Uint2Ptr textWidths, EnumFieldAssocPtr PNTR alists,
                                     Boolean useBar, Boolean noExtend,
                                     ToDialogFunc tofunc, FromDialogFunc fromfunc);

extern ForM CreateGenomeCenterForm (Int2 left, Int2 top, CharPtr title,
                                    BtnActnProc finish,
                                    BtnActnProc cancel,
                                    Boolean readPhrap,
                                    Boolean buildContig,
                                    WndActnProc activateForm)

{
  ButtoN             b;
  GrouP              c;
  Fa2htgsFormPtr     ffp;
  GrouP              g;
  GrouP              h;
  PrompT             ppt;
  GrouP              q;
  GrouP              sa;
  StdEditorProcsPtr  sepp;
  WindoW             w;
  Int2               wid;
  GrouP              x;
  GrouP              y;

  ffp = (Fa2htgsFormPtr) MemNew (sizeof (Fa2htgsForm));
  if (ffp == NULL) return NULL;
  ffp->finish = finish;
  ffp->cancel = cancel;
  ffp->readPhrap = readPhrap;
  ffp->buildContig = buildContig;

  w = FixedWindow (left, top, -10, -10, title, NULL);
  if (w == NULL) return NULL;
  SetObjectExtra (w, ffp, CleanupGenomeCenterForm);

  ffp->form = (ForM) w;
  ffp->toform = NULL;
  ffp->fromform = Fa2htgsToSeqSubmitPtr;
  ffp->formmessage = Fa2htgsFormMessage;

  ffp->ssp = NULL;
  ffp->sep = NULL;
  ffp->seplist = NULL;
  ffp->alists [0] = NULL;

#ifndef WIN_MAC
  CreateFa2htgsFormMenus (w);
#endif

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    ffp->appmessage = sepp->handleMessages;
  }

  SetGroupSpacing (w, 10, 10);

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  ffp->templateblock = HiddenGroup (h, -1, 0, NULL);
  b = PushButton (ffp->templateblock, "Read Seq-submit Template", ReadTemplate);
  SetObjectExtra (b, ffp, NULL);

  ffp->fastablock = HiddenGroup (h, -1, 0, NULL);
  if (readPhrap) {
    b = PushButton (ffp->fastablock, "Read PHRAP File", ReadAPhrapFile);
  } else if (buildContig) {
    SetGroupSpacing (ffp->fastablock, 10, 10);
    y = HiddenGroup (ffp->fastablock, 1, 0, NULL);
    StaticPrompt (y, "Coordinates are on", 0, stdLineHeight, programFont, 'c');
    ffp->contigtype = HiddenGroup (y, 3, 0, NULL);
    RadioButton (ffp->contigtype, "Individual Accessions");
    RadioButton (ffp->contigtype, "Master Sequence");
    SetValue (ffp->contigtype, 1);
    b = PushButton (ffp->fastablock, "Read CONTIG Instructions", ReadContigFile);
    AlignObjects (ALIGN_CENTER, (HANDLE) y, (HANDLE) b, NULL);
  } else {
    b = PushButton (ffp->fastablock, "Read FASTA File", ReadFastaHtgsFile);
  }
  SetObjectExtra (b, ffp, NULL);
  Hide (ffp->fastablock);

  x = HiddenGroup (h, 0, 0, NULL);

  ffp->orderblock = HiddenGroup (x, -1, 0, NULL);
  SetGroupSpacing (ffp->orderblock, 5, 5);
  ppt = StaticPrompt (ffp->orderblock, "Enter contigs in desired order", 0, 0, programFont, 'c');
  ffp->contigorder = CreateTagListDialogEx (ffp->orderblock, 8, 1, 2,
                                            contigorder_types, NULL,
                                            contigdefault_alists,
                                            TRUE, TRUE, NULL, NULL);
  b = PushButton (ffp->orderblock, "Proceed", PhrapOrderChosen);
  SetObjectExtra (b, ffp, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) ppt, (HANDLE) ffp->contigorder,
                (HANDLE) b, NULL);
  Hide (ffp->orderblock);

  ffp->controlblock = HiddenGroup (x, -1, 0, NULL);
  g = HiddenGroup (ffp->controlblock, 2, 0, NULL);

  ffp->htgsphase = NULL;
  if (! buildContig) {
    StaticPrompt (g, "Phase", 0, stdLineHeight, programFont, 'l');
    ffp->htgsphase = HiddenGroup (g, 3, 0, (GrpActnProc) SetFa2htgsAcceptBtn);
    SetObjectExtra (ffp->htgsphase, ffp, NULL);
    RadioButton (ffp->htgsphase, "HTGS-1");
    RadioButton (ffp->htgsphase, "HTGS-2");
    RadioButton (ffp->htgsphase, "HTGS-3");
  }

  StaticPrompt (g, "Organism", 0, dialogTextHeight, programFont, 'l');
  ffp->orgname = DialogText (g, "", 10, (TxtActnProc) SetFa2htgsAcceptBtn);
  SetObjectExtra (ffp->orgname, ffp, NULL);

  StaticPrompt (g, "Sequence name", 0, dialogTextHeight, programFont, 'l');
  ffp->seqname = DialogText (g, "", 10, (TxtActnProc) SetFa2htgsAcceptBtn);
  SetObjectExtra (ffp->seqname, ffp, NULL);

  StaticPrompt (g, "Length", 0, dialogTextHeight, programFont, 'l');
  ffp->knownlength = DialogText (g, "", 10, NULL);
  SetObjectExtra (ffp->knownlength, ffp, NULL);

  ffp->update = CheckBox (g, "Update", SetFa2htgsUpdate);
  SetObjectExtra (ffp->update, ffp, NULL);
  ffp->accession = DialogText (g, "", 10, (TxtActnProc) SetFa2htgsAcceptBtn);
  SetObjectExtra (ffp->accession, ffp, NULL);
  Disable (ffp->accession);

  StaticPrompt (g, "Chromosome", 0, dialogTextHeight, programFont, 'l');
  ffp->chromosome = DialogText (g, "", 10, (TxtActnProc) SetFa2htgsAcceptBtn);
  SetObjectExtra (ffp->chromosome, ffp, NULL);

  StaticPrompt (g, "Clone", 0, dialogTextHeight, programFont, 'l');
  ffp->clone = DialogText (g, "", 10, (TxtActnProc) SetFa2htgsAcceptBtn);
  SetObjectExtra (ffp->clone, ffp, NULL);

  wid = MaxStringWidths (secaccstrings) + 2;
  sa = MultiLinePrompt (g, "Secondary Accessions", wid, programFont);
  ffp->secondaries = CreateVisibleStringDialog (g, 3, -1, 15);
  AlignObjects (ALIGN_MIDDLE, (HANDLE) sa, (HANDLE) ffp->secondaries, NULL);

  q = HiddenGroup (ffp->controlblock, 1, 0, NULL);

  StaticPrompt (q, "Title", 0, 0, programFont, 'c');
  ffp->title = ScrollText (q, 20, 4, programFont, TRUE, (TxtActnProc) SetFa2htgsAcceptBtn);
  SetObjectExtra (ffp->title, ffp, NULL);

  StaticPrompt (q, "Remark", 0, 0, programFont, 'c');
  ffp->remark = ScrollText (q, 20, 4, programFont, TRUE, (TxtActnProc) SetFa2htgsAcceptBtn);
  SetObjectExtra (ffp->remark, ffp, NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) q, NULL);
  Hide (ffp->controlblock);

  c = HiddenGroup (h, 2, 0, NULL);
  ffp->okBtn = PushButton (c, "Accept", AcceptFa2htgs);
  SetObjectExtra (ffp->okBtn, ffp, NULL);
  Disable (ffp->okBtn);
  b = PushButton (c, "Cancel", CancelFa2htgs);
  SetObjectExtra (b, ffp, NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) ffp->templateblock, (HANDLE) ffp->fastablock,
                (HANDLE) ffp->orderblock, (HANDLE) ffp->controlblock, (HANDLE) c, NULL);

  RealizeWindow (w);

  if (activateForm != NULL) {
    SetActivate (w, activateForm);
  }

  Show (w);
  Select (w);
  Select (ffp->orgname);

  return (ForM) w;
}

typedef struct convcdsdata {
  FEATURE_FORM_BLOCK

  SeqEntryPtr    sep;
  TexT           geneName;
  TexT           protName;
  TexT           featcomment;
  ButtoN         retain;
  Uint2          subtype;
  Int2           errcount;
  Char           findThis [128];
  SeqLocPtr      slp;
  SeqEntryPtr    nsep;
  ObjMgrPtr      omp;
  ObjMgrTypePtr  omtp;
} ConvCdsData, PNTR ConvCdsPtr;

static Boolean CollectCDSAncestorGatherFunc (GatherContextPtr gcp)

{
  BioseqPtr      bsp;
  ConvCdsPtr     ccp;
  Boolean        noLeft;
  Boolean        noRight;
  ObjMgrTypePtr  omtp;
  SeqEntryPtr    sep;
  SeqFeatPtr     sfp;
  SeqLocPtr      slp;
  Uint2          subtype;

  if (gcp == NULL) return TRUE;

  ccp = (ConvCdsPtr) gcp->userdata;
  if (ccp == NULL ) return TRUE;

  if (gcp->thistype != OBJ_SEQFEAT) return TRUE;
  omtp = ccp->omtp;
  if (omtp == NULL || omtp->subtypefunc == NULL) return TRUE;

  sfp = (SeqFeatPtr) gcp->thisitem;
  subtype = (*(omtp->subtypefunc)) ((Pointer) sfp);
  if (subtype != ccp->subtype) return TRUE;

  bsp = GetBioseqGivenSeqLoc (sfp->location, gcp->entityID);
  if (bsp == NULL) return TRUE;
  if (ISA_aa (bsp->mol)) return TRUE;
  if (bsp->repr != Seq_repr_seg) {
    sep = ccp->nsep;
    if (sep == NULL || sep->choice != 1) return TRUE;
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp == NULL) return TRUE;
  }
  CheckSeqLocForPartial (sfp->location, &noLeft, &noRight);
  slp = SeqLocMerge (bsp, sfp->location, ccp->slp, FALSE, TRUE, FALSE);
  if (slp == NULL) return TRUE;
  SetSeqLocPartial (slp, noLeft, noRight);

  ccp->slp = SeqLocFree (ccp->slp);
  ccp->slp = slp;

  return TRUE;
}

static void FinishConvertingToCDS (SeqEntryPtr sep, Uint2 entityID, ConvCdsPtr ccp)

{
  ByteStorePtr  bs;
  BioseqPtr     bsp;
  Char          ch;
  CdRegionPtr   crp;
  ValNodePtr    descr;
  Uint1         frame;
  Int2          genCode;
  Int2          i;
  Int4          len;
  Int4          lens [4];
  Int4          max;
  Boolean       noLeft;
  Boolean       noRight;
  MolInfoPtr    mip;
  SeqEntryPtr   nsep;
  SeqEntryPtr   old;
  CharPtr       prot;
  ProtRefPtr    prp;
  SeqEntryPtr   psep;
  CharPtr       ptr;
  SeqFeatPtr    sfp;
  Char          str [128];
  ValNodePtr    vnp;

  if (sep == NULL || ccp == NULL) return;
  genCode = SeqEntryToGeneticCode (sep, NULL, NULL, 0);
  crp = CreateNewCdRgn (1, FALSE, genCode);
  if (crp == NULL) return;
  sfp = CreateNewFeature (ccp->nsep, NULL, SEQFEAT_CDREGION, NULL);
  if (sfp == NULL) {
    CdRegionFree (crp);
    return;
  }
  sfp->data.value.ptrvalue = (Pointer) crp;
  sfp->location = SeqLocFree (sfp->location);
  sfp->location = AsnIoMemCopy ((Pointer) ccp->slp,
                                (AsnReadFunc) SeqLocAsnRead,
                                (AsnWriteFunc) SeqLocAsnWrite);
  CheckSeqLocForPartial (sfp->location, &noLeft, &noRight);
  sfp->partial = (sfp->partial || noLeft || noRight);
  if (! TextHasNoText (ccp->featcomment)) {
    sfp->comment = SaveStringFromTextAndStripNewlines (ccp->featcomment);
  }
  max = 0;
  frame = 0;
  for (i = 1; i <= 3; i++) {
    crp->frame = (Uint1) i;
    bs = ProteinFromCdRegionEx (sfp, FALSE, FALSE);
    len = BSLen (bs);
    BSFree (bs);
    lens [i] = len;
    if (len > max) {
      max = len;
      frame = (Uint1) i;
    }
  }
  for (i = 1; i <= 3; i++) {
    if (lens [i] == max && i != frame) {
      (ccp->errcount)++;
    }
  }
  crp->frame = frame;
  bs = ProteinFromCdRegionEx (sfp, TRUE, FALSE);
  if (bs == NULL) return;
  prot = BSMerge (bs, NULL);
  bs = BSFree (bs);
  if (prot == NULL) return;
  ptr = prot;
  ch = *ptr;
  while (ch != '\0') {
    *ptr = TO_UPPER (ch);
    ptr++;
    ch = *ptr;
  }
  i = (Int2) StringLen (prot);
  if (i > 0 && prot [i - 1] == '*') {
    prot [i - 1] = '\0';
  }
  bs = BSNew (1000);
  if (bs != NULL) {
    ptr = prot;
    /*
    if (prot [0] == '-') {
      ptr++;
    }
    */
    BSWrite (bs, (VoidPtr) ptr, (Int4) StringLen (ptr));
  }
  MemFree (prot);
  if (bs == NULL) return;
  bsp = BioseqNew ();
  if (bsp == NULL) return;
  bsp->repr = Seq_repr_raw;
  bsp->mol = Seq_mol_aa;
  bsp->seq_data_type = Seq_code_ncbieaa;
  bsp->seq_data = bs;
  bsp->length = BSLen (bs);
  bs = NULL;
  old = SeqEntrySetScope (sep);
  bsp->id = MakeNewProteinSeqId (sfp->location, NULL);
  SeqMgrAddToBioseqIndex (bsp);
  SeqEntrySetScope (old);
  psep = SeqEntryNew ();
  if (psep != NULL) {
    psep->choice = 1;
    psep->data.ptrvalue = (Pointer) bsp;
    SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, psep);
    mip = MolInfoNew ();
    if (mip != NULL) {
      mip->biomol = 8;
      mip->tech = 8;
      if (noLeft && noRight) {
        mip->completeness = 5;
      } else if (noLeft) {
        mip->completeness = 3;
      } else if (noRight) {
        mip->completeness = 4;
      }
      vnp = CreateNewDescriptor (psep, Seq_descr_molinfo);
      if (vnp != NULL) {
        vnp->data.ptrvalue = (Pointer) mip;
      }
    }
    descr = ExtractBioSourceAndPubs (sep);
    /*
    AddSeqEntryToSeqEntry (sep, psep, FALSE);
    */
    AddSeqEntryToSeqEntry (sep, psep, TRUE);
    nsep = FindNucSeqEntry (sep);
    ReplaceBioSourceAndPubs (sep, descr);
    SetSeqFeatProduct (sfp, bsp);
    GetTitle (ccp->protName, str, sizeof (str));
    if (! StringHasNoText (str)) {
      prp = CreateNewProtRef (str, NULL, NULL, NULL);
      if (prp != NULL) {
        sfp = CreateNewFeature (psep, NULL, SEQFEAT_PROT, NULL);
        if (sfp != NULL) {
          sfp->data.value.ptrvalue = (Pointer) prp;
          SetSeqLocPartial (sfp->location, noLeft, noRight);
          sfp->partial = (sfp->partial || noLeft || noRight);
        }
      }
    }
  }
}

static void RemoveCDSAncestorCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  ConvCdsPtr     ccp;
  SeqAnnotPtr    nextsap;
  SeqFeatPtr     nextsfp;
  ObjMgrTypePtr  omtp;
  Pointer PNTR   prevsap;
  Pointer PNTR   prevsfp;
  SeqAnnotPtr    sap;
  SeqFeatPtr     sfp;
  Uint2          subtype;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  ccp = (ConvCdsPtr) mydata;
  if (ccp == NULL) return;
  omtp = ccp->omtp;
  if (omtp == NULL || omtp->subtypefunc == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
    prevsap = (Pointer PNTR) &(bsp->annot);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
    prevsap = (Pointer PNTR) &(bssp->annot);
  } else return;
  while (sap != NULL) {
    nextsap = sap->next;
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      prevsfp = (Pointer PNTR) &(sap->data);
      while (sfp != NULL) {
        nextsfp = sfp->next;
        subtype = (*(omtp->subtypefunc)) ((Pointer) sfp);
        if (subtype == ccp->subtype) {
          *(prevsfp) = sfp->next;
          sfp->next = NULL;
          SeqFeatFree (sfp);
        } else {
          prevsfp = (Pointer PNTR) &(sfp->next);
        }
        sfp = nextsfp;
      }
    }
    if (sap->data == NULL) {
      *(prevsap) = sap->next;
      sap->next = NULL;
      SeqAnnotFree (sap);
    } else {
      prevsap = (Pointer PNTR) &(sap->next);
    }
    sap = nextsap;
  }
}

static void ConvertToCDSCallback (SeqEntryPtr sep, Uint2 entityID, ConvCdsPtr ccp)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  GeneRefPtr    grp;
  GatherScope   gs;
  SeqLocPtr     gslp;
  Boolean       hasNulls;
  Boolean       noLeft;
  Boolean       noRight;
  SeqEntryPtr   nsep;
  SeqFeatPtr    sfp;
  SeqIdPtr      sip;
  SeqLocPtr     slp;
  Char          str [128];

  if (sep == NULL || ccp == NULL) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && (bssp->_class == 7 || bssp->_class == 13 ||
                         bssp->_class == 14 || bssp->_class == 15)) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        ConvertToCDSCallback (sep, entityID, ccp);
      }
      return;
    }
  }
  ccp->nsep = FindNucSeqEntry (sep);
  if (ccp->nsep == NULL) return;
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = FALSE;
  MemSet ((Pointer) (gs.ignore), (int)(TRUE), (size_t) (OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  gs.scope = sep;
  ccp->slp = NULL;
  GatherEntity (entityID, (Pointer) ccp, CollectCDSAncestorGatherFunc, &gs);
  if (ccp->slp != NULL) {
    CheckSeqLocForPartial (ccp->slp, &noLeft, &noRight);
    sip = SeqLocId (ccp->slp);
    if (sip != NULL) {
      bsp = BioseqFind (sip);
      if (bsp != NULL && ISA_na (bsp->mol)) {
        slp = SegLocToParts (bsp, ccp->slp);
        if (slp != NULL) {
          ccp->slp = SeqLocFree (ccp->slp);
          ccp->slp = slp;
          FreeAllFuzz (ccp->slp);
          SetSeqLocPartial (ccp->slp, noLeft, noRight);
        }
      }
    }
    FinishConvertingToCDS (sep, entityID, ccp);
    nsep = FindNucSeqEntry (sep);
    GetTitle (ccp->geneName, str, sizeof (str));
    if (! StringHasNoText (str)) {
      grp = CreateNewGeneRef (str, NULL, NULL, FALSE);
      if (grp != NULL) {
        if (ExtendGene (grp, nsep, ccp->slp)) {
          grp = GeneRefFree (grp);
        } else {
          sfp = CreateNewFeature (nsep, NULL, SEQFEAT_GENE, NULL);
          if (sfp != NULL) {
            sfp->data.value.ptrvalue = (Pointer) grp;
            sfp->location = SeqLocFree (sfp->location);
            sfp->location = AsnIoMemCopy ((Pointer) ccp->slp,
                                          (AsnReadFunc) SeqLocAsnRead,
                                          (AsnWriteFunc) SeqLocAsnWrite);
            sip = SeqLocId (sfp->location);
            if (sip != NULL) {
              bsp = BioseqFind (sip);
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
                  SetSeqLocPartial (sfp->location, noLeft, noRight);
                  sfp->partial = (sfp->partial || noLeft || noRight);
                }
              }
            }
          }
        }
      }
    }
  }
  ccp->slp = SeqLocFree (ccp->slp);
}

static void DoConvertToCDS (ButtoN b)

{
  ConvCdsPtr  ccp;
  CharPtr     plural;
  Char        str [128];

  ccp = GetObjectExtra (b);
  if (ccp == NULL) {
    Remove (ParentWindow (b));
    return;
  }
  GetTitle (ccp->protName, str, sizeof (str));
  if (StringHasNoText (str)) {
    Message (MSG_OK, "Protein name is required");
    return;
  }
  Hide (ccp->form);
  WatchCursor ();

  ccp->omp = ObjMgrGet ();
  if (ccp->omp != NULL) {
    ccp->omtp = ObjMgrTypeFind (ccp->omp, OBJ_SEQFEAT, NULL, NULL);
    if (ccp->omtp != NULL && ccp->omtp->subtypefunc != NULL) {
      ConvertToCDSCallback (ccp->sep, ccp->input_entityID, ccp);
      if (! GetStatus (ccp->retain)) {
        SeqEntryExplore (ccp->sep, (Pointer) ccp, RemoveCDSAncestorCallback);
      }
    }
  }

  ArrowCursor ();
  Update ();

  if (ccp->errcount > 0) {
    if (ccp->errcount > 1) {
      plural = "records";
    } else {
      plural = "record";
    }
    Message (MSG_ERROR, "Possible ambiguous frames detected in %d %s",
             (int) ccp->errcount, plural);
  }

  ObjMgrSetDirtyFlag (ccp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, ccp->input_entityID, 0, 0);
  Remove (ccp->form);
  Update ();
}

extern void PrepareToConvertToCDS (SeqEntryPtr sep, Uint2 entityID,
                                   Uint2 subtype, CharPtr findthis)

{
  ButtoN             b;
  GrouP              c;
  ConvCdsPtr         ccp;
  GrouP              g;
  GrouP              h;
  StdEditorProcsPtr  sepp;
  WindoW             w;

  if (sep == NULL || entityID == 0 || subtype == 0) return;

  ccp = (ConvCdsPtr) MemNew (sizeof (ConvCdsData));
  if (ccp == NULL) return;
  w = FixedWindow (-50, -33, -10, -10, "Convert to CDS", StdCloseWindowProc);
  SetObjectExtra (w, ccp, StdCleanupFormProc);
  ccp->form = (ForM) w;
  ccp->formmessage = DefaultMessageProc;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    ccp->appmessage = sepp->handleMessages;
  }

  ccp->input_entityID = entityID;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  ccp->sep = sep;
  ccp->subtype = subtype;
  StringNCpy_0 (ccp->findThis, findthis, sizeof (ccp->findThis));
  ccp->errcount = 0;

  g = HiddenGroup (h, 2, 0, NULL);
  StaticPrompt (g, "Gene Symbol", 0, dialogTextHeight, programFont, 'l');
  ccp->geneName = DialogText (g, "", 20, NULL);
  StaticPrompt (g, "Protein Name", 0, dialogTextHeight, programFont, 'l');
  ccp->protName = DialogText (g, "", 20, NULL);
  StaticPrompt (g, "Comment", 0, 4 * Nlm_stdLineHeight, programFont, 'l');
  ccp->featcomment = ScrollText (g, 20, 4, programFont, TRUE, NULL);

  ccp->retain = CheckBox (h, "Retain original features", NULL);

  c = HiddenGroup (h, 4, 0, NULL);
  b = PushButton (c, "Accept", DoConvertToCDS);
  SetObjectExtra (b, ccp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, (HANDLE) ccp->retain, NULL);
  RealizeWindow (w);
  Show (w);
  Select (ccp->geneName);
  Update ();
}

typedef struct alignform {
  FORM_MESSAGE_BLOCK
  DoC            doc;
} AlignForm, PNTR AlignFormPtr;

static ParData bioseqParFmt = {FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0};
static ColData bioseqColFmt [] = {
  {0, 5, 0, 0, NULL, 'l', TRUE, FALSE, FALSE, FALSE, TRUE}
};

static ParData annotParFmt = {FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0};
static ColData annotColFmt [] = {
  {0, 10, 0, 0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, TRUE}
};

static ParData alignParFmt = {FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0};
static ColData alignColFmt [] = {
  {0, 15, 0, 0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, FALSE},
  {0, 5, 0, 0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, FALSE},
  {0, 0, 0, 0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, TRUE}
};

static Boolean AlignFormPopulateProc (GatherContextPtr gcp)

{
  AlignFormPtr      afp;
  SeqAlignPtr       align;
  Char              annotDB [32];
  Uint1             annot_type;
  BioseqContextPtr  bcp;
  BioseqPtr         bsp;
  DenseDiagPtr      ddp;
  DenseSegPtr       dsp;
  ProtRefPtr        prp;
  CharPtr           ptr;
  SeqAnnotPtr       sap;
  SeqEntryPtr       sep;
  SeqFeatPtr        sfp;
  SeqIdPtr          sip;
  StdSegPtr         ssp;
  Char              str [256];
  Char              tmp [128];

  if (gcp == NULL) return TRUE;

  afp = (AlignFormPtr) gcp->userdata;
  if (afp == NULL ) return TRUE;

  switch (gcp->thistype) {
    case OBJ_BIOSEQ :
      bsp = (BioseqPtr) gcp->thisitem;
      if (bsp != NULL) {
        SeqIdWrite (bsp->id, str, PRINTID_REPORT, sizeof (str));
        bcp = BioseqContextNew (bsp);
        sfp = BioseqContextGetSeqFeat (bcp, SEQFEAT_PROT, NULL, NULL, 0);
        BioseqContextFree (bcp);
        if (sfp != NULL) {
          prp = (ProtRefPtr) sfp->data.value.ptrvalue;
          if (prp != NULL) {
            if (prp->name != NULL && (! StringHasNoText (prp->name->data.ptrvalue))) {
              StringCat (str, " (");
              StringCat (str, (CharPtr) prp->name->data.ptrvalue);
              StringCat (str, ")");
            } else if (! StringHasNoText (prp->desc)) {
              StringCat (str, " (");
              StringCat (str, (CharPtr) prp->desc);
              StringCat (str, ")");
            }
          }
        }
        StringCat (str, ":\n");
        AppendText (afp->doc, str, &bioseqParFmt, bioseqColFmt, programFont);
      }
      break;
    case OBJ_SEQANNOT :
      sap = (SeqAnnotPtr) gcp->thisitem;
      if (sap != NULL && sap->type == 2) {
        get_align_annot_qual (sap, annotDB, sizeof (annotDB), &annot_type);
        if (annot_type == ANNOT_BLAST) {
          StringCpy (str, annotDB);
          StringCat (str, ":\n");
          AppendText (afp->doc, str, &annotParFmt, annotColFmt, programFont);
        }
      }
      break;
    case OBJ_SEQALIGN :
    case OBJ_SEQHIST_ALIGN :
      align = (SeqAlignPtr) gcp->thisitem;
      sip = NULL;
      if (align->segtype == 1) {
        ddp = (DenseDiagPtr) align->segs;
        if (ddp != NULL) {
          for (sip = ddp->id; sip != NULL && sip->next != NULL; sip = sip->next) {
          }
        }
      } else if (align->segtype == 2) {
        dsp = (DenseSegPtr) align->segs;
        if (dsp != NULL) {
          for (sip = dsp->ids; sip != NULL && sip->next != NULL; sip = sip->next) {
          }
        }
      } else if (align->segtype == 3) {
        ssp = (StdSegPtr) align->segs;
        if (ssp != NULL) {
          for (sip = ssp->ids; sip != NULL && sip->next != NULL; sip = sip->next) {
          }
        }
      }
      if (sip != NULL) {
        bsp = BioseqLockById (sip);
        if (bsp != NULL) {
          /* SeqIdWrite (bsp->id, str, PRINTID_FASTA_LONG, sizeof (str)); */
          SeqIdWrite (bsp->id, str, PRINTID_REPORT, sizeof (str));
          StringNCpy_0 (tmp, BioseqGetTitle (bsp), sizeof (tmp));
          StringCat (str, "\t");
          if (! StringHasNoText (tmp)) {
            StringCat (str, tmp);
          }
          tmp [0] = '\0';
          sep = SeqMgrGetSeqEntryForData (bsp);
          if (sep != NULL) {
            SeqEntryToBioSource (sep, NULL, tmp, sizeof (tmp), NULL);
            if (! StringHasNoText (tmp)) {
              ptr = StringStr (tmp, "(");
              if (ptr != NULL) {
                *ptr = '\0';
              }
              StringCat (str, " [");
              StringCat (str, tmp);
              StringCat (str, "]");
            }
          }
          StringCat (str, "\t");
          sprintf (tmp, "%d %d %d", (int) gcp->entityID,
                   (int) gcp->itemID, (int) gcp->thistype);
          StringCat (str, tmp);
          StringCat (str, "\n");
          AppendText (afp->doc, str, &alignParFmt, alignColFmt, programFont);
        }
        BioseqUnlock (bsp);
      }
      return TRUE;
    default :
      break;
  }
  return TRUE;
}

static void PopulateAlignForm (ForM f, SeqEntryPtr sep)

{
  AlignFormPtr  afp;
  GatherScope   gs;
  RecT          r;
  Int2          width;

  if (f == NULL || sep == NULL) return;
  afp = (AlignFormPtr) GetObjectExtra (f);
  if (afp == NULL) return;
  Reset (afp->doc);
  SetDocAutoAdjust (afp->doc, FALSE);
  ObjectRect (afp->doc, &r);
  InsetRect (&r, 4, 4);
  width = r.right - r.left;
  alignColFmt [0].pixWidth = width / 5;
  alignColFmt [1].pixWidth = width - alignColFmt [0].pixWidth;
  bioseqColFmt [0].pixWidth = width;
  annotColFmt [0].pixWidth = width;
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  MemSet((Pointer) (gs.ignore), (int) (TRUE), (size_t) (OBJ_MAX * sizeof (Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQALIGN] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  gs.ignore[OBJ_SEQHIST] = FALSE;
  gs.ignore[OBJ_SEQHIST_ALIGN] = FALSE;
  gs.scope = sep;
  GatherEntity (afp->input_entityID, (Pointer) afp, AlignFormPopulateProc, &gs);
  AdjustDocScroll (afp->doc);
}

static void AlignMessageProc (ForM f, Int2 mssg)

{
  FILE          *fp;
  AlignFormPtr  afp;
  Char          path [PATH_MAX];

  afp = (AlignFormPtr) GetObjectExtra (f);
  if (afp != NULL) {
    switch (mssg) {
      case VIB_MSG_CLOSE :
        Remove (f);
        break;
      case VIB_MSG_PRINT :
        PrintDocument (afp->doc);
        break;
      case VIB_MSG_EXPORT :
        if (GetOutputFileName (path, sizeof (path), "align.txt")) {
          WatchCursor ();
#ifdef WIN_MAC
          fp = FileOpen (path, "r");
          if (fp != NULL) {
            FileClose (fp);
          } else {
            FileCreate (path, "TEXT", "ttxt");
          }
#endif
          fp = FileOpen (path, "w");
          if (fp != NULL) {
            SaveDocument (afp->doc, fp);
            FileClose (fp);
          }
          ArrowCursor ();
        }
        break;
      default :
        if (afp->appmessage != NULL) {
          afp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static void AlignFormActivate (WindoW w)

{
  AlignFormPtr  afp;
  IteM          exportItm;
  IteM          printItm;

  afp = (AlignFormPtr) GetObjectExtra (w);
#ifdef WIN_MAC
  currentFormDataPtr = (VoidPtr) afp;
#endif
  if (afp != NULL) {
    if (afp->activate != NULL) {
      afp->activate (w);
    }
    exportItm = FindFormMenuItem ((BaseFormPtr) afp, VIB_MSG_EXPORT);
    SafeSetTitle (exportItm, "Export Align Summary...");
    SafeEnable (exportItm);
    printItm = FindFormMenuItem ((BaseFormPtr) afp, VIB_MSG_PRINT);
    SafeEnable (printItm);
  }
}

static void AlignNotifyProc (DoC d, Int2 item, Int2 row, Int2 col, Boolean dblClick)

{
  unsigned int  entityID;
  unsigned int  itemID;
  unsigned int  itemtype;
  CharPtr       txt;

  if (d == NULL || item < 1 || row < 1 || col < 1) return;
  txt = GetDocText (d, item, row, 3);
  if (! StringHasNoText (txt)) {
    if (sscanf (txt, "%u %u %u", &entityID, &itemID, &itemtype) == 3) {
      if (shftKey) {
        ObjMgrAlsoSelect (entityID, itemID, itemtype, 0, NULL);
      } else {
        ObjMgrSelect (entityID, itemID, itemtype, 0, NULL);
      }
    }
  }
  MemFree (txt);
}

static ForM CreateAlignForm (Uint2 entityID)

{
  AlignFormPtr       afp;
  GrouP              c;
  GrouP              h;
  StdEditorProcsPtr  sepp;
  WindoW             w;
#ifndef WIN_MAC
  MenU               m;
#endif

  w = NULL;
  afp = MemNew (sizeof (AlignForm));
  if (afp != NULL) {
    w = FixedWindow (-50, -33, -10, -10, "Alignment Summary", StdCloseWindowProc);
    SetObjectExtra (w, afp, StdCleanupFormProc);
    afp->form = (ForM) w;
    afp->formmessage = AlignMessageProc;
    afp->input_entityID = entityID;

    sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
    if (sepp != NULL) {
      afp->appmessage = sepp->handleMessages;
    }

#ifndef WIN_MAC
    m = PulldownMenu (w, "File");
    FormCommandItem (m, "Close", (BaseFormPtr) afp, VIB_MSG_CLOSE);
    SeparatorItem (m);
    FormCommandItem (m, "Export...", (BaseFormPtr) afp, VIB_MSG_EXPORT);
    SeparatorItem (m);
    FormCommandItem (m, "Print...", (BaseFormPtr) afp, VIB_MSG_PRINT);
#endif

    h = HiddenGroup (w, -1, 0, NULL);
    SetGroupSpacing (h, 10, 10);

    afp->doc = DocumentPanel (h, 35 * stdCharWidth, 20 * stdLineHeight);
    SetObjectExtra (afp->doc, afp, NULL);
    SetDocNotify (afp->doc, AlignNotifyProc);

    c = HiddenGroup (w, 4, 0, NULL);
    PushButton (c, "Cancel", StdCancelButtonProc);

    AlignObjects (ALIGN_CENTER, (HANDLE) afp->doc, (HANDLE) c, NULL);

    RealizeWindow (w);
    SetActivate (w, AlignFormActivate);
  }
  return (ForM) w;
}

extern void ViewAlignmentSummary (IteM i)

{
  BaseFormPtr  bfp;
  ForM         f;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL || bfp->input_entityID == 0) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  WatchCursor ();
  Update ();
  f = CreateAlignForm (bfp->input_entityID);
  PopulateAlignForm (f, sep);
  ArrowCursor ();
  Update ();
  Show (f);
  Select (f);
}

#define REMOVE_CDSET   1
#define CONVERT_CDSET  2
#define EDIT_CDSET     3
#define ADD_CDSET      4

static ENUM_ALIST(cds_gene_prot_field_alist)
  {" ",                    0},
  {"CDS comment",          1},
  {"Gene locus",           2},
  {"Gene description",     3},
  {"Gene allele",          4},
  {"Gene maploc",          5},
  {"Gene synonym",         6},
  {"Gene comment",         7},
  {"Protein name",         8},
  {"Protein description",  9},
  {"Protein E.C. number", 10},
  {"Protein activity",    11},
  {"Protein comment",     12},
END_ENUM_ALIST

static ENUM_ALIST(prot_subtype_alist)
  {" ",                    5},
  {"Product",              0},
  {"Pre-protein",          1},
  {"Mature",               2},
  {"Signal peptide",       3},
  {"Transit peptide",      4},
END_ENUM_ALIST

typedef struct cdsetformdata {
  FEATURE_FORM_BLOCK

  Int2           type;
  PopuP          fromfield;
  PopuP          tofield;
  PrompT         example;
  ButtoN         leaveOn;

  Int2           fromval;
  Int2           toval;

  Boolean        onlyGene;
  Boolean        onlyProt;

  Boolean        replaceOldAsked;
  Boolean        doReplaceAll;

  TexT           findthis;
  TexT           replacewith;
  CharPtr        findStr;
  CharPtr        replaceStr;
  Char           temp [128];

  PopuP          protSubType;
  GrouP          protSubGrp;
} CdsetFormData, PNTR CdsetFormPtr;

typedef struct geneprotfind {
  SeqFeatPtr     cds;
  SeqFeatPtr     gene;
  SeqFeatPtr     prot;
  SeqLocPtr      location;
  SeqLocPtr      product;
  Int4           mingene;
  Int4           minprot;
} GeneProtFind, PNTR GeneProtPtr;

static Boolean GeneProtFindFunc (GatherContextPtr gcp)

{
  Int4         diff;
  GeneProtPtr  gpp;
  SeqFeatPtr   sfp;

  if (gcp == NULL) return TRUE;
  gpp = (GeneProtPtr) gcp->userdata;
  if (gpp == NULL) return TRUE;
  if (gcp->thistype != OBJ_SEQFEAT) return TRUE;
  sfp = (SeqFeatPtr) gcp->thisitem;
  if (sfp == NULL || sfp->data.value.ptrvalue == NULL) return TRUE;
  if (sfp->data.choice == SEQFEAT_GENE && gpp->location != NULL) {
    diff = SeqLocAinB (gpp->location, sfp->location);
    if (diff >= 0) {
      if (diff < gpp->mingene) {
        gpp->gene = sfp;
        gpp->mingene = diff;
      }
    }
  } else if (sfp->data.choice == SEQFEAT_PROT && gpp->product != NULL) {
    diff = SeqLocAinB (gpp->product, sfp->location);
    if (diff >= 0) {
      if (diff < gpp->minprot) {
        gpp->prot = sfp;
        gpp->minprot = diff;
      }
    }
  }
  return TRUE;
}

void FindGeneAndProtForCDS (Uint2 entityID, SeqFeatPtr cds,
                            SeqFeatPtr PNTR gene, SeqFeatPtr PNTR prot)

{
  GeneProtFind  gpf;
  GatherScope   gs;

  if (gene != NULL) {
    *gene = NULL;
  }
  if (prot != NULL) {
    *prot = NULL;
  }
  if (entityID == 0 || cds == NULL) return;
  MemSet ((Pointer) (&gpf), 0, sizeof (GeneProtFind));
  gpf.cds = cds;
  gpf.gene = NULL;
  gpf.prot = NULL;
  gpf.location = cds->location;
  gpf.product = cds->product;
  gpf.mingene = INT4_MAX;
  gpf.minprot = INT4_MAX;
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = FALSE;
  MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  GatherEntity (entityID, (Pointer) &gpf, GeneProtFindFunc, &gs);
  if (gene != NULL) {
    *gene = gpf.gene;
  }
  if (prot != NULL) {
    *prot = gpf.prot;
  }
}

static Boolean ShouldReplaceString (CdsetFormPtr cfp)

{
  MsgAnswer  ans;

  if (cfp == NULL) return FALSE;
  if (cfp->type == EDIT_CDSET) {
    cfp->doReplaceAll = TRUE;
    cfp->replaceOldAsked = TRUE;
  }
  if (! cfp->replaceOldAsked) {
    cfp->replaceOldAsked = TRUE;
    ArrowCursor ();
    ans = Message (MSG_YN, "Do you wish to overwrite existing contents?");
    WatchCursor ();
    Update ();
    cfp->doReplaceAll = (Boolean) (ans == ANS_YES);
  }
  return cfp->doReplaceAll;
}

static CharPtr SaveOrReplaceString (CdsetFormPtr cfp, CharPtr str, CharPtr current)

{
  size_t   len;
  CharPtr  tmp;

  if (cfp == NULL) return str;
  if (StringHasNoText (str) && cfp->type != EDIT_CDSET) return current;
  if (current != NULL) {
    if (ShouldReplaceString (cfp)) {
      MemFree (current);
    } else {
      len = StringLen (str) + StringLen (current) + 4;
      tmp = MemNew (len);
      if (tmp != NULL) {
        StringCpy (tmp, current);
        StringCat (tmp, "; ");
        StringCat (tmp, str);
        MemFree (str);
        MemFree (current);
        return tmp;
      }
    }
  }
  return str;
}

static void EditCdsetString (CharPtr PNTR strptr, CdsetFormPtr cfp, CharPtr foundit)

{
  size_t   diff;
  size_t   foundlen;
  size_t   replen;
  CharPtr  newstring;
  CharPtr  tmp;
  CharPtr  tmp2;

  if (cfp == NULL || strptr == NULL || foundit == NULL) return;
  foundlen = StringLen (cfp->findStr);
  replen = StringLen (cfp->replaceStr);
  if (replen > foundlen) {
    diff = replen - foundlen;
  } else {
    diff = foundlen - replen;
  }
  newstring = MemNew (StringLen (*strptr) + diff + 1);
  tmp = *strptr;
  tmp2 = newstring;
  while (tmp != foundit) {
    *tmp2 = *tmp;
    tmp++;
    tmp2++;
  }
  if (cfp->replaceStr != NULL) {
    tmp2 = MemCopy (tmp2, cfp->replaceStr, replen);
  }
  tmp2 += replen;
  tmp += foundlen;
  tmp2 = StringMove (tmp2, tmp);
  *strptr = MemFree (*strptr);
  *strptr = newstring;
}

static Boolean ProcessEachCDSFunc (GatherContextPtr gcp)

{
  SeqFeatPtr    cds;
  CdsetFormPtr  cfp;
  CharPtr       foundit;
  SeqFeatPtr    gene;
  GeneRefPtr    grp;
  SeqFeatPtr    prot;
  ProtRefPtr    prp;
  SeqFeatPtr    sfp;
  CharPtr       str;
  UIEnum        val;
  ValNodePtr    vnp;

  if (gcp == NULL) return TRUE;
  cfp = (CdsetFormPtr) gcp->userdata;
  if (cfp == NULL) return TRUE;
  if (gcp->thistype != OBJ_SEQFEAT) return TRUE;
  sfp = (SeqFeatPtr) gcp->thisitem;
  if (sfp == NULL || sfp->data.value.ptrvalue == NULL) return TRUE;
  cds = NULL;
  gene = NULL;
  grp = NULL;
  prot = NULL;
  prp = NULL;
  if (cfp->onlyGene) {
    if (sfp->data.choice == SEQFEAT_GENE) {
      gene = sfp;
      grp = (GeneRefPtr) gene->data.value.ptrvalue;
    } else return TRUE;
  } else if (cfp->onlyProt) {
    if (sfp->data.choice == SEQFEAT_PROT) {
      prot = sfp;
      prp = (ProtRefPtr) prot->data.value.ptrvalue;
    } else return TRUE;
  } else if (sfp->data.choice == SEQFEAT_CDREGION) {
    cds = sfp;
    FindGeneAndProtForCDS (cfp->input_entityID, cds, &gene, &prot);
    if (gene != NULL) {
      grp = (GeneRefPtr) gene->data.value.ptrvalue;
    }
    if (prot != NULL) {
      prp = (ProtRefPtr) prot->data.value.ptrvalue;
    }
  } else return TRUE;
  if (cfp->protSubType != NULL && prot != NULL && prp != NULL) {
    if (GetEnumPopup (cfp->protSubType, prot_subtype_alist, &val)) {
      if (val < 5 && (Uint1) val != prp->processed) {
        prot = NULL;
        prp = NULL;
      }
    }
  }
  str = NULL;
  switch (cfp->type) {
    case REMOVE_CDSET :
      break;
    case ADD_CDSET :
      str = StringSave (cfp->findStr);
      break;
    case EDIT_CDSET :
    case CONVERT_CDSET :
      switch (cfp->fromval) {
        case 1 :
          if (cds != NULL) {
            str = StringSave (cds->comment);
          }
          break;
        case 2 :
          if (grp != NULL) {
            str = StringSave (grp->locus);
          }
          break;
        case 3 :
          if (grp != NULL) {
            str = StringSave (grp->desc);
          }
          break;
        case 4 :
          if (grp != NULL) {
            str = StringSave (grp->allele);
          }
          break;
        case 5 :
          if (grp != NULL) {
            str = StringSave (grp->maploc);
          }
          break;
        case 6 :
          if (grp != NULL) {
            vnp = grp->syn;
            if (vnp != NULL) {
              str = StringSave (vnp->data.ptrvalue);
            }
          }
          break;
        case 7 :
          if (gene != NULL) {
            str = StringSave (gene->comment);
          }
          break;
        case 8 :
          if (prp != NULL) {
            vnp = prp->name;
            if (vnp != NULL) {
              str = StringSave (vnp->data.ptrvalue);
            }
          }
          break;
        case 9 :
          if (prp != NULL) {
            str = StringSave (prp->desc);
          }
          break;
        case 10 :
          if (prp != NULL) {
            vnp = prp->ec;
            if (vnp != NULL) {
              str = StringSave (vnp->data.ptrvalue);
            }
          }
          break;
        case 11 :
          if (prp != NULL) {
            vnp = prp->activity;
            if (vnp != NULL) {
              str = StringSave (vnp->data.ptrvalue);
            }
          }
          break;
        case 12 :
          if (prot != NULL) {
            str = StringSave (prot->comment);
          }
          break;
        default :
          break;
      }
      break;
    default :
      break;
  }
  if (cfp->type == EDIT_CDSET) {
    foundit = StringISearch (str, cfp->findStr);
    if (foundit != NULL) {
      EditCdsetString (&str, cfp, foundit);
    }
  }
  switch (cfp->type) {
    case REMOVE_CDSET :
      break;
    case CONVERT_CDSET :
    case ADD_CDSET :
    case EDIT_CDSET :
      switch (cfp->toval) {
        case 1 :
          if (cds != NULL) {
            cds->comment = SaveOrReplaceString (cfp, str, cds->comment);
            str = NULL;
          }
          break;
        case 2 :
          if (grp != NULL) {
            grp->locus = SaveOrReplaceString (cfp, str, grp->locus);
            str = NULL;
          }
          break;
        case 3 :
          if (grp != NULL) {
            grp->desc = SaveOrReplaceString (cfp, str, grp->desc);
            str = NULL;
          }
          break;
        case 4 :
          if (grp != NULL) {
            grp->allele = SaveOrReplaceString (cfp, str, grp->allele);
            str = NULL;
          }
          break;
        case 5 :
          if (grp != NULL) {
            grp->maploc = SaveOrReplaceString (cfp, str, grp->maploc);
            str = NULL;
          }
          break;
        case 6 :
          if (grp != NULL) {
            if (! StringHasNoText (str)) {
              if (grp->syn != NULL && ShouldReplaceString (cfp)) {
                vnp = grp->syn;
                vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
                vnp->data.ptrvalue = str;
                str = NULL;
              } else {
                vnp = ValNodeNew (NULL);
                if (vnp != NULL) {
                  vnp->next = grp->syn;
                  grp->syn = vnp;
                  vnp->data.ptrvalue = str;
                  str = NULL;
                }
              }
            }
          }
          break;
        case 7 :
          if (gene != NULL) {
            gene->comment = SaveOrReplaceString (cfp, str, gene->comment);
            str = NULL;
          }
          break;
        case 8 :
          if (prp != NULL) {
            if (! StringHasNoText (str)) {
              if (prp->name != NULL && ShouldReplaceString (cfp)) {
                vnp = prp->name;
                vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
                vnp->data.ptrvalue = str;
                str = NULL;
              } else {
                vnp = ValNodeNew (NULL);
                if (vnp != NULL) {
                  vnp->next = prp->name;
                  prp->name = vnp;
                  vnp->data.ptrvalue = str;
                  str = NULL;
                }
              }
            }
          }
          break;
        case 9 :
          if (prp != NULL) {
            prp->desc = SaveOrReplaceString (cfp, str, prp->desc);
            str = NULL;
          }
          break;
        case 10 :
          if (prp != NULL) {
            if (! StringHasNoText (str)) {
              if (prp->ec != NULL && ShouldReplaceString (cfp)) {
                vnp = prp->ec;
                vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
                vnp->data.ptrvalue = str;
                str = NULL;
              } else {
                vnp = ValNodeNew (NULL);
                if (vnp != NULL) {
                  vnp->next = prp->ec;
                  prp->ec = vnp;
                  vnp->data.ptrvalue = str;
                  str = NULL;
                }
              }
            }
          }
          break;
        case 11 :
          if (prp != NULL) {
            if (! StringHasNoText (str)) {
              if (prp->activity != NULL && ShouldReplaceString (cfp)) {
                vnp = prp->activity;
                vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
                vnp->data.ptrvalue = str;
                str = NULL;
              } else {
                vnp = ValNodeNew (NULL);
                if (vnp != NULL) {
                  vnp->next = prp->activity;
                  prp->activity = vnp;
                  vnp->data.ptrvalue = str;
                  str = NULL;
                }
              }
            }
          }
          break;
        case 12 :
          if (prot != NULL) {
            prot->comment = SaveOrReplaceString (cfp, str, prot->comment);
            str = NULL;
          }
          break;
        default :
          break;
      }
      break;
    default :
      break;
  }
  MemFree (str);
  return TRUE;
}

static Boolean CleanupEachCDSFunc (GatherContextPtr gcp)

{
  SeqFeatPtr    cds;
  CdsetFormPtr  cfp;
  SeqFeatPtr    gene;
  GeneRefPtr    grp;
  SeqFeatPtr    prot;
  ProtRefPtr    prp;
  SeqFeatPtr    sfp;

  if (gcp == NULL) return TRUE;
  cfp = (CdsetFormPtr) gcp->userdata;
  if (cfp == NULL) return TRUE;
  if (gcp->thistype != OBJ_SEQFEAT) return TRUE;
  sfp = (SeqFeatPtr) gcp->thisitem;
  if (sfp == NULL || sfp->data.value.ptrvalue == NULL) return TRUE;
  cds = NULL;
  gene = NULL;
  grp = NULL;
  prot = NULL;
  prp = NULL;
  if (sfp->data.choice == SEQFEAT_GENE) {
    gene = sfp;
    grp = (GeneRefPtr) gene->data.value.ptrvalue;
  } else if (sfp->data.choice == SEQFEAT_PROT) {
    prot = sfp;
    prp = (ProtRefPtr) prot->data.value.ptrvalue;
  } else if (sfp->data.choice == SEQFEAT_CDREGION) {
    cds = sfp;
    /*
    FindGeneAndProtForCDS (cfp->input_entityID, cds, &gene, &prot);
    if (gene != NULL) {
      grp = (GeneRefPtr) gene->data.value.ptrvalue;
    }
    if (prot != NULL) {
      prp = (ProtRefPtr) prot->data.value.ptrvalue;
    }
    */
  } else return TRUE;
  switch (cfp->type) {
    case REMOVE_CDSET :
    case CONVERT_CDSET :
      switch (cfp->fromval) {
        case 1 :
          if (cds != NULL) {
            cds->comment = MemFree (cds->comment);
          }
          break;
        case 2 :
          if (grp != NULL) {
            grp->locus = MemFree (grp->locus);
          }
          break;
        case 3 :
          if (grp != NULL) {
            grp->desc = MemFree (grp->desc);
          }
          break;
        case 4 :
          if (grp != NULL) {
            grp->allele = MemFree (grp->allele);
          }
          break;
        case 5 :
          if (grp != NULL) {
            grp->maploc = MemFree (grp->maploc);
          }
          break;
        case 6 :
          if (grp != NULL) {
            grp->syn = ValNodeFreeData (grp->syn);
          }
          break;
        case 7 :
          if (gene != NULL) {
            gene->comment = MemFree (gene->comment);
          }
          break;
        case 8 :
          if (prp != NULL) {
            prp->name = ValNodeFreeData (prp->name);
          }
          break;
        case 9 :
          if (prp != NULL) {
            prp->desc = MemFree (prp->desc);
          }
          break;
        case 10 :
          if (prp != NULL) {
            prp->ec = ValNodeFreeData (prp->ec);
          }
          break;
        case 11 :
          if (prp != NULL) {
            prp->activity = ValNodeFreeData (prp->activity);
          }
          break;
        case 12 :
          if (prot != NULL) {
            prot->comment = MemFree (prot->comment);
          }
          break;
        default :
          break;
      }
      break;
    case ADD_CDSET :
      break;
    case EDIT_CDSET :
      break;
    default :
      break;
  }
  return TRUE;
}

static void AddProtFeatIfMissing (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr    bsp;
  ProtRefPtr   prp;
  SeqAnnotPtr  sap;
  SeqFeatPtr   sfp;

  if (sep == NULL) return;
  if (! IS_Bioseq (sep)) return;
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL) return;
  if (! ISA_aa (bsp->mol)) return;
  sap = bsp->annot;
  while (sap != NULL) {
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      while (sfp != NULL) {
        if (sfp->data.choice == SEQFEAT_PROT) {
          return;
        }
      }
      sfp = sfp->next;
    }
    sap = sap->next;
  }
  sfp = CreateNewFeature (sep, NULL, SEQFEAT_PROT, NULL);
  if (sfp == NULL) return;
  prp = ProtRefNew ();
  sfp->data.value.ptrvalue = (Pointer) prp;
}

CharPtr JustSaveStringFromText (TexT t)

{
  size_t   len;
  CharPtr  str;

  len = TextLength (t);
  if (len > 0) {
    str = (CharPtr) MemNew(len + 1);
    if (str != NULL) {
      GetTitle (t, str, len + 1);
      /* TrimSpacesAroundString (str); */
      if (StringHasNoText (str)) {
        str = (CharPtr) MemFree(str);
      }
      return str;
    } else {
      return NULL;
    }
  } else {
    return NULL;
  }
}

extern void CleanupEmptyFeatCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent);
static void DoProcessCDSet (ButtoN b)

{
  CdsetFormPtr  cfp;
  GatherScope   gs;
  SeqEntryPtr   sep;
  UIEnum        val;

  cfp = (CdsetFormPtr) GetObjectExtra (b);
  if (cfp == NULL || cfp->input_entityID == 0) return;
  Hide (cfp->form);
  WatchCursor ();
  Update ();
  cfp->fromval = 0;
  cfp->toval = 0;
  if (GetEnumPopup (cfp->fromfield, cds_gene_prot_field_alist, &val)) {
    cfp->fromval = (Int2) val;
  }
  if (GetEnumPopup (cfp->tofield, cds_gene_prot_field_alist, &val)) {
    cfp->toval = (Int2) val;
  }
  if (cfp->type == ADD_CDSET || cfp->type == EDIT_CDSET) {
    cfp->toval = cfp->fromval;
  }
  cfp->onlyGene = (Boolean) ((cfp->toval >= 2 && cfp->toval <= 7) &&
                             (cfp->fromval >= 2 && cfp->fromval <= 7));
  cfp->onlyProt = (Boolean) ((cfp->toval >= 8 && cfp->toval <= 12) &&
                             (cfp->fromval >= 8 && cfp->fromval <= 12));
  if (cfp->type == ADD_CDSET || cfp->type == CONVERT_CDSET) {
    if (cfp->toval >= 8 && cfp->toval <= 12) {
      sep = GetTopSeqEntryForEntityID (cfp->input_entityID);
      SeqEntryExplore (sep, NULL, AddProtFeatIfMissing);
    }
  }
  cfp->replaceOldAsked = FALSE;
  cfp->doReplaceAll = FALSE;
  cfp->findStr = JustSaveStringFromText (cfp->findthis);
  cfp->replaceStr = JustSaveStringFromText (cfp->replacewith);
  if (cfp->type != EDIT_CDSET) {
    TrimSpacesAroundString (cfp->findStr);
    TrimSpacesAroundString (cfp->replaceStr);
  }
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = FALSE;
  MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  if (cfp->type != REMOVE_CDSET) {
    GatherEntity (cfp->input_entityID, (Pointer) cfp, ProcessEachCDSFunc, &gs);
  }
  if (cfp->type == REMOVE_CDSET ||
      (cfp->type == CONVERT_CDSET && (! GetStatus (cfp->leaveOn)))) {
    GatherEntity (cfp->input_entityID, (Pointer) cfp, CleanupEachCDSFunc, &gs);
  }
  /* SeqEntryExplore (sep, (Pointer) cfp, CleanupEmptyFeatCallback); */
  ArrowCursor ();
  Update ();
  if (indexerVersion && cfp->type == ADD_CDSET) {
    sep = GetTopSeqEntryForEntityID (cfp->input_entityID);
    if (CountSeqEntryComponents (sep) == 1) {
      switch (cfp->fromval) {
        case 1 :
          Message (MSG_OK, "When only one record present, edit the CdRgn feature directly");
          break;
        case 2 :
        case 3 :
        case 4 :
        case 5 :
        case 6 :
        case 7 :
          Message (MSG_OK, "When only one record present, edit the gene feature directly");
          break;
        case 8 :
        case 9 :
        case 10 :
        case 11 :
        case 12 :
          Message (MSG_OK, "When only one record present, edit the protein feature directly");
          break;
        default :
          break;
      }
    }
  }
  ObjMgrSetDirtyFlag (cfp->input_entityID, TRUE);
  sep = GetTopSeqEntryForEntityID (cfp->input_entityID);
  SeqEntryExplore (sep, NULL, CleanupEmptyFeatCallback);
  ObjMgrSendMsg (OM_MSG_UPDATE, cfp->input_entityID, 0, 0);
  Remove (cfp->form);
}

static Boolean SetCDSExample (GatherContextPtr gcp)

{
  SeqFeatPtr    cds;
  CdsetFormPtr  cfp;
  SeqFeatPtr    gene;
  GeneRefPtr    grp;
  size_t        len;
  SeqFeatPtr    prot;
  ProtRefPtr    prp;
  SeqFeatPtr    sfp;
  CharPtr       str;
  CharPtr       tmp;
  ValNodePtr    vnp;

  if (gcp == NULL) return TRUE;
  cfp = (CdsetFormPtr) gcp->userdata;
  if (cfp == NULL) return TRUE;
  if (gcp->thistype != OBJ_SEQFEAT) return TRUE;
  sfp = (SeqFeatPtr) gcp->thisitem;
  if (sfp == NULL || sfp->data.value.ptrvalue == NULL) return TRUE;
  cds = NULL;
  gene = NULL;
  grp = NULL;
  prot = NULL;
  prp = NULL;
  if (sfp->data.choice == SEQFEAT_GENE) {
    gene = sfp;
    grp = (GeneRefPtr) gene->data.value.ptrvalue;
  } else if (sfp->data.choice == SEQFEAT_PROT) {
    prot = sfp;
    prp = (ProtRefPtr) prot->data.value.ptrvalue;
  } else if (sfp->data.choice == SEQFEAT_CDREGION) {
    cds = sfp;
  } else return TRUE;
  str = NULL;
  switch (cfp->fromval) {
    case 1 :
      if (cds != NULL) {
        str = cds->comment;
      }
      break;
    case 2 :
      if (grp != NULL) {
        str = grp->locus;
      }
      break;
    case 3 :
      if (grp != NULL) {
        str = grp->desc;
      }
      break;
    case 4 :
      if (grp != NULL) {
        str = grp->allele;
      }
      break;
    case 5 :
      if (grp != NULL) {
        str = grp->maploc;
      }
      break;
    case 6 :
      if (grp != NULL) {
        vnp = grp->syn;
        if (vnp != NULL) {
          str = vnp->data.ptrvalue;
        }
      }
      break;
    case 7 :
      if (gene != NULL) {
        str = gene->comment;
      }
      break;
    case 8 :
      if (prp != NULL) {
        vnp = prp->name;
        if (vnp != NULL) {
          str = vnp->data.ptrvalue;
        }
      }
      break;
    case 9 :
      if (prp != NULL) {
        str = prp->desc;
      }
      break;
    case 10 :
      if (prp != NULL) {
        vnp = prp->ec;
        if (vnp != NULL) {
          str = vnp->data.ptrvalue;
        }
      }
      break;
    case 11 :
      if (prp != NULL) {
        vnp = prp->activity;
        if (vnp != NULL) {
          str = vnp->data.ptrvalue;
        }
      }
      break;
    case 12 :
      if (prot != NULL) {
        str = prot->comment;
      }
      break;
    default :
      break;
  }
  if (str != NULL) {
    if (! StringHasNoText (str)) {
      if (cfp->type == EDIT_CDSET) {
        SafeSetTitle (cfp->findthis, str);
      }
      StringNCpy_0 (cfp->temp, str, sizeof (cfp->temp));
      len = StringLen (str) + 12;
      tmp = MemNew (len);
      if (tmp != NULL) {
        StringCpy (tmp, "Example (");
        StringCat (tmp, str);
        StringCat (tmp, ")");
        SafeSetTitle (cfp->example, tmp);
        MemFree (tmp);
        return FALSE;
      } else {
        SafeSetTitle (cfp->example, str);
        return FALSE;
      }
    }
  }
  return TRUE;
}

static void ChangeCDSetExample (PopuP p)

{
  CdsetFormPtr  cfp;
  GatherScope   gs;
  UIEnum        val;

  cfp = (CdsetFormPtr) GetObjectExtra (p);
  if (cfp == NULL || cfp->input_entityID == 0) return;
  SafeSetTitle (cfp->example, "");
  if (GetEnumPopup (cfp->fromfield, cds_gene_prot_field_alist, &val)) {
    cfp->fromval = (Int2) val;
    MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
    gs.seglevels = 1;
    gs.get_feats_location = FALSE;
    MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
    gs.ignore[OBJ_BIOSEQ] = FALSE;
    gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
    gs.ignore[OBJ_SEQFEAT] = FALSE;
    gs.ignore[OBJ_SEQANNOT] = FALSE;
    GatherEntity (cfp->input_entityID, (Pointer) cfp, SetCDSExample, &gs);
    if (cfp->fromval >= 8 && cfp->fromval <= 12) {
      SafeShow (cfp->protSubGrp);
    } else {
      SafeHide (cfp->protSubGrp);
    }
  } else {
    SafeSetTitle (cfp->findthis, "");
  }
  Update ();
}

static void CDSetMessageProc (ForM f, Int2 mssg)

{
  CdsetFormPtr  cfp;

  cfp = (CdsetFormPtr) GetObjectExtra (f);
  if (cfp != NULL) {
    if (cfp->appmessage != NULL) {
      cfp->appmessage (f, mssg);
    }
  }
}

static CharPtr  cdsetLegend = "\
CDS\n\
        /gene=\"gene locus\"\n\
        /map=\"gene maploc\" (moves to source with taxon fixup)\n\
        /EC_number=\"protein E.C. number\"\n\
        /function=\"protein activity\"\n\
        /note=\"gene description; gene allele; gene synonym;\n\
        cds comment; protein description; protein comment\"\n\
        /product=\"protein name\"\n\
gene\n\
        /gene=\"gene locus\"\n\
        /allele=\"gene allele\"\n\
        /note=\"gene comment; gene description; gene synonym\"\n\
        /map=\"gene maploc\" (moves to source with taxon fixup)\n";

static void CDSetLegend (ButtoN b)

{
  GrouP   legend;
  WindoW  w;

  w = FixedWindow (-90, -33, -10, -10, "CDS-Gene-Prot Legend", StdCloseWindowProc);
  legend = MultiLinePromptEx (w, cdsetLegend, stdCharWidth * 35, programFont, FALSE);
  b = PushButton (w, "Cancel", StdCancelButtonProc);
  AlignObjects (ALIGN_CENTER, (HANDLE) legend, (HANDLE) b, NULL);
  RealizeWindow (w);
  Show (w);
  Select (w);
}

static void CleanupCdsetForm (GraphiC g, VoidPtr data)

{
  CdsetFormPtr  cfp;

  cfp = (CdsetFormPtr) data;
  if (cfp != NULL) {
    MemFree (cfp->findStr);
    MemFree (cfp->replaceStr);
  }
  StdCleanupFormProc (g, data);
}

static void ProcessCDSet (IteM i, Int2 type)

{
  EnumFieldAssocPtr  ap;
  ButtoN             b, b1;
  BaseFormPtr        bfp;
  GrouP              c;
  CdsetFormPtr       cfp;
  GrouP              g;
  GatherScope        gs;
  GrouP              h;
  GrouP              samples;
  SeqEntryPtr        sep;
  StdEditorProcsPtr  sepp;
  CharPtr            title;
  WindoW             w;
  GrouP              x;
  GrouP              y;
  GrouP              z;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  cfp = (CdsetFormPtr) MemNew (sizeof (CdsetFormData));
  if (cfp == NULL) return;
  cfp->type = type;
  switch (type) {
    case REMOVE_CDSET :
      title = "CDS-Gene-Prot Removal";
      break;
    case CONVERT_CDSET :
      title = "CDS-Gene-Prot Conversion";
      break;
    case EDIT_CDSET :
      title = "Edit CDS-Gene-Prot";
      break;
    case ADD_CDSET :
      title = "Add CDS-Gene-Prot";
      break;
    default :
      title = "?";
      break;
  }
  w = FixedWindow (-50, -33, -10, -10, title, StdCloseWindowProc);
  SetObjectExtra (w, cfp, CleanupCdsetForm);
  cfp->form = (ForM) w;
  cfp->formmessage = CDSetMessageProc;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    cfp->appmessage = sepp->handleMessages;
  }

  cfp->input_entityID = bfp->input_entityID;
  cfp->input_itemID = bfp->input_itemID;
  cfp->input_itemtype = bfp->input_itemtype;

  z = HiddenGroup (w, -2, 0, NULL);
  if (indexerVersion) {
    samples = HiddenGroup (z, 2, 0, NULL);
    for (ap = cds_gene_prot_field_alist; ap->name != NULL; ap++) {
      if (ap->value != 0) {
        cfp->fromval = (Int2) ap->value;
        cfp->temp [0] = '\0';
        MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
        gs.seglevels = 1;
        gs.get_feats_location = FALSE;
        MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
        gs.ignore[OBJ_BIOSEQ] = FALSE;
        gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
        gs.ignore[OBJ_SEQFEAT] = FALSE;
        gs.ignore[OBJ_SEQANNOT] = FALSE;
        GatherEntity (cfp->input_entityID, (Pointer) cfp, SetCDSExample, &gs);
        StaticPrompt (samples, ap->name, 0, 0, programFont, 'l');
        StaticPrompt (samples, cfp->temp, 14 * stdCharWidth, 0, programFont, 'l');
      }
    }
  }

  h = HiddenGroup (z, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  g = HiddenGroup (h, 4, 0, NULL);
  switch (type) {
    case ADD_CDSET :
      StaticPrompt (g, "Type", 0, popupMenuHeight, programFont, 'l');
      cfp->fromfield = PopupList (g, TRUE, ChangeCDSetExample);
      SetObjectExtra (cfp->fromfield, cfp, NULL);
      InitEnumPopup (cfp->fromfield, cds_gene_prot_field_alist, NULL);
      SetEnumPopup (cfp->fromfield, cds_gene_prot_field_alist, 0);
      break;
    case REMOVE_CDSET :
      StaticPrompt (g, "Remove", 0, popupMenuHeight, programFont, 'l');
      cfp->fromfield = PopupList (g, TRUE, ChangeCDSetExample);
      SetObjectExtra (cfp->fromfield, cfp, NULL);
      InitEnumPopup (cfp->fromfield, cds_gene_prot_field_alist, NULL);
      SetEnumPopup (cfp->fromfield, cds_gene_prot_field_alist, 0);
      break;
    case CONVERT_CDSET :
      StaticPrompt (g, "From", 0, popupMenuHeight, programFont, 'l');
      cfp->fromfield = PopupList (g, TRUE, ChangeCDSetExample);
      SetObjectExtra (cfp->fromfield, cfp, NULL);
      InitEnumPopup (cfp->fromfield, cds_gene_prot_field_alist, NULL);
      SetEnumPopup (cfp->fromfield, cds_gene_prot_field_alist, 0);

      StaticPrompt (g, "To", 0, popupMenuHeight, programFont, 'l');
      cfp->tofield = PopupList (g, TRUE, NULL);
      InitEnumPopup (cfp->tofield, cds_gene_prot_field_alist, NULL);
      SetEnumPopup (cfp->tofield, cds_gene_prot_field_alist, 0);
      break;
    case EDIT_CDSET :
      StaticPrompt (g, "Type", 0, popupMenuHeight, programFont, 'l');
      cfp->fromfield = PopupList (g, TRUE, ChangeCDSetExample);
      SetObjectExtra (cfp->fromfield, cfp, NULL);
      InitEnumPopup (cfp->fromfield, cds_gene_prot_field_alist, NULL);
      SetEnumPopup (cfp->fromfield, cds_gene_prot_field_alist, 0);
      break;
    default :
      break;
  }

  x = HiddenGroup (h, 4, 0, NULL);
  /*StaticPrompt (x, "Example", 0, 0, programFont, 'l');*/
  cfp->example = StaticPrompt (x, "", 14 * stdCharWidth, 0, programFont, 'l');

  y = NULL;
  cfp->protSubGrp = NULL;
  switch (type) {
    case REMOVE_CDSET :
      /*
      y = HiddenGroup (h, 2, 0, NULL);
      StaticPrompt (y, "Optional string constraint", 0, dialogTextHeight, programFont, 'c');
      cfp->findthis = DialogText (y, "", 14, NULL);
      */
      break;
    case CONVERT_CDSET :
      /*
      y = HiddenGroup (h, 2, 0, NULL);
      StaticPrompt (y, "Optional string constraint", 0, dialogTextHeight, programFont, 'c');
      cfp->findthis = DialogText (y, "", 14, NULL);
      */
      y = HiddenGroup (h, 2, 0, NULL);
      cfp->leaveOn = CheckBox (y, "Leave on original", NULL);
      break;
    case ADD_CDSET :
      y = HiddenGroup (h, 2, 0, NULL);
      StaticPrompt (y, "Text", 0, dialogTextHeight, programFont, 'c');
      cfp->findthis = DialogText (y, "", 14, NULL);
      break;
    case EDIT_CDSET :
      cfp->protSubGrp = HiddenGroup (h, 2, 0, NULL);
      StaticPrompt (cfp->protSubGrp, "Protein subtype", 0, dialogTextHeight, programFont, 'l');
      cfp->protSubType = PopupList (cfp->protSubGrp, TRUE, NULL);
      InitEnumPopup (cfp->protSubType, prot_subtype_alist, NULL);
      SetEnumPopup (cfp->protSubType, prot_subtype_alist, (UIEnum) 5);
      Hide (cfp->protSubGrp);

      y = HiddenGroup (h, 2, 0, NULL);
      StaticPrompt (y, "Find", 0, dialogTextHeight, programFont, 'c');
      cfp->findthis = DialogText (y, "", 14, NULL);
      StaticPrompt (y, "Replace", 0, dialogTextHeight, programFont, 'c');
      cfp->replacewith = DialogText (y, "", 14, NULL);
      break;
    default :
      break;
  }

  b1 = PushButton (h, "Legend", CDSetLegend);

  c = HiddenGroup (w, 4, 0, NULL);
  b = PushButton (c, "Accept", DoProcessCDSet);
  SetObjectExtra (b, cfp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) b1, (HANDLE) y,
                (HANDLE) cfp->protSubGrp, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) h, (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Select (cfp->findthis);
  Update ();
}

extern void AddCDSet (IteM i)

{
  ProcessCDSet (i, ADD_CDSET);
}

extern void RemoveCDSet (IteM i)

{
  ProcessCDSet (i, REMOVE_CDSET);
}

extern void EditCDSet (IteM i)

{
  ProcessCDSet (i, EDIT_CDSET);
}

extern void ConvertCDSet (IteM i)

{
  ProcessCDSet (i, CONVERT_CDSET);
}

typedef struct findform {
  FORM_MESSAGE_BLOCK
  TexT            findTxt;
  TexT            replaceTxt;
  ButtoN          caseCounts;
  ButtoN          wholeWord;
  ButtoN          findAllBtn;
  ButtoN          replaceAllBtn;
} FindForm, PNTR FindFormPtr;

extern CharPtr CompressSpaces (CharPtr str)

{
  Char     ch;
  CharPtr  dst;
  Char     last;
  CharPtr  ptr;

  if (str != NULL && str [0] != '\0') {
    dst = str;
    ptr = str;
    ch = *ptr;
    while (ch != '\0' && ch <= ' ') {
      ptr++;
      ch = *ptr;
    }
    while (ch != '\0') {
      *dst = ch;
      dst++;
      ptr++;
      last = ch;
      ch = *ptr;
      if (ch != '\0' && ch < ' ') {
        *ptr = ' ';
        ch = *ptr;
      }
      while (ch != '\0' && last <= ' ' && ch <= ' ') {
        ptr++;
        ch = *ptr;
      }
    }
    *dst = '\0';
    dst = NULL;
    ptr = str;
    ch = *ptr;
    while (ch != '\0') {
      if (ch != ' ') {
        dst = NULL;
      } else if (dst == NULL) {
        dst = ptr;
      }
      ptr++;
      ch = *ptr;
    }
    if (dst != NULL) {
      *dst = '\0';
    }
  }
  return str;
}

extern CharPtr SearchForString (CharPtr str, CharPtr sub, Boolean case_counts, Boolean whole_word)

{
  CharPtr  ptr = NULL;
  CharPtr  tmp;

  if (case_counts) {
    ptr = StringSearch (str, sub);
  } else {
    ptr = StringISearch (str, sub);
  }
  if (ptr == NULL) return NULL;
  if (whole_word) {
    if (ptr > str) {
      tmp = ptr - 1;
      if (! IS_WHITESP (*tmp)) return NULL;
    }
    tmp = ptr + StringLen (sub);
    if (*tmp != '\0' && (! IS_WHITESP (*tmp))) return NULL;
  }
  return ptr;
}

static void FindInFlatFile (Uint2 entityID, Uint1 format, Uint1 mode,
                            CharPtr sub, Boolean case_counts, Boolean whole_word)

{
  Asn2ffJobPtr     ajp;
  Boolean          already;
  DescrStructPtr   descr;
  unsigned int     entID;
  Int4             index;
  unsigned int     itemID;
  unsigned int     itemtype;
  ErrSev           level;
  FFPrintArrayPtr  pap = NULL;
  Int4             pap_size;
  SeqEntryPtr      sep;
  SelStructPtr     sel;
  CharPtr          string;

  ajp = (Asn2ffJobPtr) MemNew (sizeof (Asn2ffJob));
  if (ajp == NULL) return;
  level = ErrSetMessageLevel (SEV_MAX);
  WatchCursor ();
  Update ();

  sep = GetTopSeqEntryForEntityID (entityID);
  ajp->sep = sep;
  ajp->mode = mode;
  ajp->format = format;
  ajp->gb_style = TRUE;
  ajp->show_seq = TRUE;
  ajp->show_gi = TRUE;
  ajp->error_msgs = FALSE;
  ajp->non_strict = TRUE;
  ajp->Spop = spop;
  ajp->show_gene = TRUE;
  if (IsAGenomeRecord (sep)) {
    ajp->only_one = TRUE;
    ajp->genome_view = TRUE;
  }

  pap_size = asn2ff_setup (ajp, &pap);
  if (pap_size > 0) {
    asn2ff_set_output (NULL, "\n");
    for (index = 0; index < pap_size; index++) {
      string = FFPrint (pap, index, pap_size);
      if (string != NULL && *string != '\0') {
        CompressSpaces (string);
        if (SearchForString (string, sub, case_counts, whole_word) != NULL) {
          if (pap [index].descr) {
            descr = pap [index].descr;
            entID = descr->entityID;
            itemID = descr->itemID;
            itemtype = descr->itemtype;
            pap [index].descr = MemFree (pap [index].descr);
            already = FALSE;
            for (sel = ObjMgrGetSelected (); sel != NULL; sel = sel->next) {
              if (sel->entityID == entID && sel->itemID == itemID && sel->itemtype == itemtype) {
                already = TRUE;
              }
            }
            if (! already) {
              ObjMgrAlsoSelect (entID, itemID, itemtype, 0, NULL);
            }
          }
        }
      }
      MemFree (string);
    }
    MemFree (pap);
  }
  asn2ff_cleanup (ajp);

  MemFree (ajp);
  ErrSetMessageLevel (level);
  ArrowCursor ();
  Update ();
}

static void FindFlatProc (ButtoN b)

{
  Boolean      caseCounts;
  FindFormPtr  ffp;
  Char         findme [256];
  Boolean      wholeWord;

  ffp = (FindFormPtr) GetObjectExtra (b);
  if (ffp == NULL) return;
  GetTitle (ffp->findTxt, findme, sizeof (findme) - 1);
  ObjMgrDeSelect (0, 0, 0, 0, NULL);
  caseCounts = GetStatus (ffp->caseCounts);
  wholeWord = GetStatus (ffp->wholeWord);
  CompressSpaces (findme);
  FindInFlatFile (ffp->input_entityID, GENBANK_FMT, SEQUIN_MODE,
                  findme, caseCounts, wholeWord);
}

static void FindTextProc (TexT t)

{
  FindFormPtr  ffp;

  ffp = (FindFormPtr) GetObjectExtra (t);
  if (ffp != NULL) {
    if (TextLength (t) > 0) {
      SafeEnable (ffp->findAllBtn);
      SafeEnable (ffp->replaceAllBtn);
    } else {
      SafeDisable (ffp->findAllBtn);
      SafeDisable (ffp->replaceAllBtn);
    }
  }
}

static void CommonFindReplaceProc (ButtoN b, Boolean replace, Boolean replaceAll)

{
  Boolean      caseCounts;
  Char         changeme [256];
  FindFormPtr  ffp;
  Char         findme [256];
  Boolean      wholeWord;

  ffp = (FindFormPtr) GetObjectExtra (b);
  if (ffp != NULL) {
    GetTitle (ffp->findTxt, findme, sizeof (findme) - 1);
    ObjMgrDeSelect (0, 0, 0, 0, NULL);
    caseCounts = GetStatus (ffp->caseCounts);
    wholeWord = GetStatus (ffp->wholeWord);
    if (replace) {
      GetTitle (ffp->replaceTxt, changeme, sizeof (changeme) - 1);
      FindInEntity (ffp->input_entityID, findme, changeme, TRUE, UPDATE_NEVER, caseCounts, wholeWord, replaceAll);
      GetRidOfEmptyFeatsDescStrings (ffp->input_entityID, NULL);
      ObjMgrSetDirtyFlag (ffp->input_entityID, TRUE);
      ObjMgrSendMsg (OM_MSG_UPDATE, ffp->input_entityID, 0, 0);
    } else {
      FindInEntity (ffp->input_entityID, findme, NULL, TRUE, UPDATE_ONCE, caseCounts, wholeWord, FALSE);
    }
    Update ();
  }
}

static void FindAllProc (ButtoN b)

{
  CommonFindReplaceProc (b, FALSE, FALSE);
}

static void ReplaceAllProc (ButtoN b)

{
  CommonFindReplaceProc (b, TRUE, TRUE);
}

static void ClearFindTextProc (ButtoN b)

{
  FindFormPtr  ffp;

  ffp = (FindFormPtr) GetObjectExtra (b);
  if (ffp == NULL) return;
  SafeSetTitle (ffp->findTxt, "");
  SafeSetTitle (ffp->replaceTxt, "");
  Select (ffp->findTxt);
}

static void FindFormMessage (ForM f, Int2 mssg)

{
  FindFormPtr  ffp;

  ffp = (FindFormPtr) GetObjectExtra (f);
  if (ffp != NULL) {
    switch (mssg) {
      case VIB_MSG_CLOSE :
        Remove (f);
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
        if (ffp->appmessage != NULL) {
          ffp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static ForM CreateFindForm (Int2 left, Int2 top, CharPtr title,
                            Uint2 entityID, Boolean flatfile)

{
  ButtoN             b;
  GrouP              c;
  GrouP              g;
  FindFormPtr        ffp;
  GrouP              j;
  GrouP              q;
  StdEditorProcsPtr  sepp;
  WindoW             w;

  w = NULL;
  ffp = MemNew (sizeof (FindForm));
  if (ffp != NULL) {
    w = FixedWindow (left, top, -10, -10, title, StdCloseWindowProc);
    SetObjectExtra (w, ffp, StdCleanupFormProc);
    ffp->form = (ForM) w;
    ffp->formmessage = FindFormMessage;
    ffp->input_entityID = entityID;

#ifndef WIN_MAC
    CreateStdEditorFormMenus (w);
#endif

    sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
    if (sepp != NULL) {
      SetActivate (w, sepp->activateForm);
      ffp->appmessage = sepp->handleMessages;
    }

    j = HiddenGroup (w, -1, 0, NULL);
    SetGroupSpacing (j, 10, 10);

    g = HiddenGroup (j, 2, 0, NULL);
    StaticPrompt (g, "Find", 0, dialogTextHeight, programFont, 'l');
    ffp->findTxt = DialogText (g, "", 25, FindTextProc);
    SetObjectExtra (ffp->findTxt, ffp, NULL);
    if (! flatfile) {
      StaticPrompt (g, "Replace", 0, dialogTextHeight, programFont, 'l');
      ffp->replaceTxt = DialogText (g, "", 25, NULL);
      SetObjectExtra (ffp->replaceTxt, ffp, NULL);
    }

    q = HiddenGroup (w, 2, 0, NULL);
    ffp->caseCounts = CheckBox (q, "Case Sensitive", NULL);
    ffp->wholeWord = CheckBox (q, "Entire Word", NULL);

    c = HiddenGroup (w, 4, 0, NULL);
    SetGroupSpacing (c, 10, 2);
    if (! flatfile) {
      ffp->findAllBtn = DefaultButton (c, "Find All", FindAllProc);
      SetObjectExtra (ffp->findAllBtn, ffp, NULL);
      Disable (ffp->findAllBtn);
      ffp->replaceAllBtn = PushButton (c, "Replace All", ReplaceAllProc);
      SetObjectExtra (ffp->replaceAllBtn, ffp, NULL);
      Disable (ffp->replaceAllBtn);
    } else {
      ffp->findAllBtn = DefaultButton (c, "Find All", FindFlatProc);
      SetObjectExtra (ffp->findAllBtn, ffp, NULL);
      Disable (ffp->findAllBtn);
    }
    b = PushButton (c, "Clear", ClearFindTextProc);
    SetObjectExtra (b, ffp, NULL);
    PushButton (c, "Cancel", StdCancelButtonProc);

    AlignObjects (ALIGN_CENTER, (HANDLE) j, (HANDLE) q, (HANDLE) c, NULL);

    RealizeWindow (w);
    Select (ffp->findTxt);
  }
  return (ForM) w;
}

extern void FindStringProc (IteM i)

{
  BaseFormPtr  bfp;
  ForM         w;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp != NULL) {
    w = CreateFindForm (-90, -66, "Find", bfp->input_entityID, FALSE);
    Show (w);
    Select (w);
  }
}

extern void FindFlatfileProc (IteM i)

{
  BaseFormPtr  bfp;
  ForM         w;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp != NULL) {
    w = CreateFindForm (-90, -66, "Flat File Find", bfp->input_entityID, TRUE);
    Show (w);
    Select (w);
  }
}

extern EnumFieldAssoc  orgmod_subtype_alist [];
extern EnumFieldAssoc  subsource_subtype_alist [];

static void RemoveRedundantSourceNote (CharPtr str, CharPtr keyword, CharPtr name)

{
  Char     ch;
  Boolean  extract;
  CharPtr  tmp1, tmp2;

  if (str == NULL || keyword == NULL || name == NULL) return;
  while (str != NULL && *str != '\0') {
    extract = TRUE;
    tmp1 = str;
    ch = *tmp1;
    while (ch == ' ' || ch == ';') {
      tmp1++;
      ch = *tmp1;
    }
    tmp2 = keyword;
    ch = TO_UPPER (*tmp1);
    while (ch != '\0' && ch == TO_UPPER (*tmp2)) {
      tmp1++;
      tmp2++;
      ch = TO_UPPER (*tmp1);
    }
    if (*tmp2 == '\0' && ch == ' ') {
      while (ch != '\0' && ch == ' ') {
        tmp1++;
        ch = *tmp1;
      }
      tmp2 = name;
      while (ch != '\0' && ch == *tmp2) {
        tmp1++;
        tmp2++;
        ch = *tmp1;
      }
      if (*tmp2 == '\0' && (ch == '\0' || ch == ' ' || ch == ';')) {
        while (ch != '\0' && ch == ' ') {
          tmp1++;
          ch = *tmp1;
        }
        /* now ready to extract */
      } else {
        extract = FALSE;
      }
    } else {
      extract = FALSE;
    }
    if (extract) {
      while (ch == ' ' || ch == ';') {
        tmp1++;
        ch = *tmp1;
      }
      tmp2 = str;
      while (ch != '\0') {
        *tmp2 = ch;
        tmp1++;
        tmp2++;
        ch = *tmp1;
      }
      *tmp2 = '\0';
    } else {
      ch = *tmp1;
      while (ch != '\0' && ch != ';') {
        tmp1++;
        ch = *tmp1;
      }
      if (ch == ';') {
        tmp1++;
        ch = *tmp1;
      }
      while (ch != '\0' && (ch == ' ' || ch == ';')) {
        tmp1++;
        ch = *tmp1;
      }
      str = tmp1;
    }
  }
}

static CharPtr TrimSpacesAndSemicolonsAroundString (CharPtr str)

{
  Uchar    ch;	/* to use 8bit characters in multibyte languages */
  CharPtr  dst;
  CharPtr  ptr;

  if (str != NULL && str [0] != '\0') {
    dst = str;
    ptr = str;
    ch = *ptr;
    while (ch != '\0' && (ch <= ' ' || ch == ';')) {
      ptr++;
      ch = *ptr;
    }
    while (ch != '\0') {
      *dst = ch;
      dst++;
      ptr++;
      ch = *ptr;
    }
    *dst = '\0';
    dst = NULL;
    ptr = str;
    ch = *ptr;
    while (ch != '\0') {
      if (ch != ' ' && ch != ';') {
        dst = NULL;
      } else if (dst == NULL) {
        dst = ptr;
      }
      ptr++;
      ch = *ptr;
    }
    if (dst != NULL) {
      *dst = '\0';
    }
  }
  return str;
}

extern void GetRidOfRedundantSourceNotes (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent);
extern void GetRidOfRedundantSourceNotes (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  EnumFieldAssocPtr  ap;
  BioSourcePtr       biop;
  BioseqPtr          bsp;
  BioseqSetPtr       bssp;
  OrgModPtr          mod;
  OrgNamePtr         onp;
  OrgRefPtr          orp;
  ValNodePtr         sdp;
  SubSourcePtr       ssp;
  CharPtr            str1, str2;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sdp = bsp->descr;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sdp = bssp->descr;
  } else return;
  while (sdp != NULL) {
    if (sdp->choice == Seq_descr_source) {
      str1 = NULL;
      str2 = NULL;
      onp = NULL;
      biop = (BioSourcePtr) sdp->data.ptrvalue;
      if (biop != NULL) {
        orp = biop->org;
        if (orp != NULL) {
          onp = orp->orgname;
          if (onp != NULL) {
            mod = onp->mod;
            while (mod != NULL) {
              if (mod->subtype == 255) {
                str1 = mod->subname;
              }
              mod = mod->next;
            }
          }
        }
        ssp = biop->subtype;
        while (ssp != NULL) {
          if (ssp->subtype == 255) {
            str2 = ssp->name;
          }
          ssp = ssp->next;
        }
        if (str1 != NULL || str2 != NULL) {
          if (onp != NULL) {
            mod = onp->mod;
            while (mod != NULL) {
              if (mod->subtype != 255) {
                for (ap = orgmod_subtype_alist; ap->name != NULL; ap++) {
                  if (ap->value == mod->subtype) {
                    RemoveRedundantSourceNote (str1, ap->name, mod->subname);
                    RemoveRedundantSourceNote (str2, ap->name, mod->subname);
                  }
                }
              }
              mod = mod->next;
            }
          }
          ssp = biop->subtype;
          while (ssp != NULL) {
            if (ssp->subtype != 255) {
              for (ap = subsource_subtype_alist; ap->name != NULL; ap++) {
                if (ap->value == ssp->subtype) {
                  RemoveRedundantSourceNote (str1, ap->name, ssp->name);
                  RemoveRedundantSourceNote (str2, ap->name, ssp->name);
                }
              }
            }
            ssp = ssp->next;
          }
        }
      }
      TrimSpacesAndSemicolonsAroundString (str1);
      TrimSpacesAndSemicolonsAroundString (str2);
    }
    sdp = sdp->next;
  }
}

typedef struct convaccdata {
  CharPtr        currID;
  CharPtr        newID;
} ConvAccData, PNTR ConvAccPtr;

static void ChangeSeqIdListToAccLocalID (SeqIdPtr sip, CharPtr currID, CharPtr newID)

{
  ObjectIdPtr   oip;

  while (sip != NULL) {
    switch (sip->choice) {
      case SEQID_LOCAL :
        oip = (ObjectIdPtr) sip->data.ptrvalue;
        if (oip != NULL) {
          if (StringCmp (oip->str, currID) == 0) {
            MemFree (oip->str);
            oip->str = StringSave (newID);
          }
        }
        break;
      default :
        break;
    }
    sip = sip->next;
  }
}

static void ChangeSeqLocListToAccLocalID (SeqLocPtr slp, CharPtr currID, CharPtr newID)

{
  SeqLocPtr      loc;
  PackSeqPntPtr  psp;
  SeqBondPtr     sbp;
  SeqIntPtr      sinp;
  SeqIdPtr       sip;
  SeqPntPtr      spp;

  while (slp != NULL) {
    switch (slp->choice) {
      case SEQLOC_NULL :
        break;
      case SEQLOC_EMPTY :
      case SEQLOC_WHOLE :
        sip = (SeqIdPtr) slp->data.ptrvalue;
        ChangeSeqIdListToAccLocalID (sip, currID, newID);
        break;
      case SEQLOC_INT :
        sinp = (SeqIntPtr) slp->data.ptrvalue;
        if (sinp != NULL) {
          sip = sinp->id;
          ChangeSeqIdListToAccLocalID (sip, currID, newID);
        }
        break;
      case SEQLOC_PNT :
        spp = (SeqPntPtr) slp->data.ptrvalue;
        if (spp != NULL) {
          sip = spp->id;
          ChangeSeqIdListToAccLocalID (sip, currID, newID);
        }
        break;
      case SEQLOC_PACKED_PNT :
        psp = (PackSeqPntPtr) slp->data.ptrvalue;
        if (psp != NULL) {
          sip = psp->id;
          ChangeSeqIdListToAccLocalID (sip, currID, newID);
        }
        break;
      case SEQLOC_PACKED_INT :
      case SEQLOC_MIX :
      case SEQLOC_EQUIV :
        loc = (SeqLocPtr) slp->data.ptrvalue;
        while (loc != NULL) {
          ChangeSeqIdListToAccLocalID (loc, currID, newID);
          loc = loc->next;
        }
        break;
      case SEQLOC_BOND :
        sbp = (SeqBondPtr) slp->data.ptrvalue;
        if (sbp != NULL) {
          spp = (SeqPntPtr) sbp->a;
          if (spp != NULL) {
            sip = spp->id;
            ChangeSeqIdListToAccLocalID (sip, currID, newID);
          }
          spp = (SeqPntPtr) sbp->b;
          if (spp != NULL) {
            sip = spp->id;
            ChangeSeqIdListToAccLocalID (sip, currID, newID);
          }
        }
        break;
      case SEQLOC_FEAT :
        break;
      default :
        break;
    }
    slp = slp->next;
  }
}

static void ChangeAlignListToAccLocalID (SeqAlignPtr align, CharPtr currID, CharPtr newID)

{
  DenseDiagPtr  ddp;
  DenseSegPtr   dsp;
  StdSegPtr     ssp;

  if (align == NULL) return;
  if (align->segtype == 1) {
    ddp = (DenseDiagPtr) align->segs;
    if (ddp != NULL) {
      ChangeSeqIdListToAccLocalID (ddp->id, currID, newID);
    }
  } else if (align->segtype == 2) {
    dsp = (DenseSegPtr) align->segs;
    if (dsp != NULL) {
      ChangeSeqIdListToAccLocalID (dsp->ids, currID, newID);
    }
  } else if (align->segtype == 3) {
    ssp = (StdSegPtr) align->segs;
    if (ssp != NULL) {
       ChangeSeqLocListToAccLocalID (ssp->loc, currID, newID);
    }
  }
}

static void CopyAccToLocalCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  ConvAccPtr    cap;
  SeqAlignPtr   sal;
  SeqAnnotPtr   sap;
  SeqFeatPtr    sfp;
  SeqIdPtr      sip;
  SeqLocPtr     slp;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  cap = (ConvAccPtr) mydata;
  if (cap == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sip = bsp->id;
    ChangeSeqIdListToAccLocalID (sip, cap->currID, cap->newID);
    SeqMgrReplaceInBioseqIndex (bsp);
    sap = bsp->annot;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
  } else return;
  while (sap != NULL) {
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      while (sfp != NULL) {
        slp = sfp->location;
        ChangeSeqLocListToAccLocalID (slp, cap->currID, cap->newID);
        slp = sfp->product;
        ChangeSeqLocListToAccLocalID (slp, cap->currID, cap->newID);
        sfp = sfp->next;
      }
    } else if (sap->type == 2) {
      sal = (SeqAlignPtr) sap->data;
      while (sal != NULL) {
        ChangeAlignListToAccLocalID (sal, cap->currID, cap->newID);
        sal = sal->next;
      }
    }
    sap = sap->next;
  }
}

static void RecordAccToConvert (ValNodePtr PNTR vnpp, BioseqPtr bsp, CharPtr newID)

{
  ConvAccPtr   cap;
  ObjectIdPtr  oip;
  SeqIdPtr     sip;

  if (vnpp == NULL || bsp == NULL || StringHasNoText (newID)) return;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_LOCAL) {
      oip = (ObjectIdPtr) sip->data.ptrvalue;
      if (oip != NULL) {
        if (! StringHasNoText (oip->str)) {
          cap = (ConvAccPtr) MemNew (sizeof (ConvAccData));
          if (cap == NULL) return;
          cap->currID = StringSave (oip->str);
          cap->newID = StringSave (newID);
          ValNodeAddPointer (vnpp, 0, (Pointer) cap);
          return;
        }
      }
    }
  }
}

static Boolean IsLegalAccession (CharPtr acnum, Int2 i)

{
  Int2  j;

  if (! isalpha (acnum [0])) return FALSE;
  if (!(isdigit(acnum[1]) && i == 6) && !(isalpha(acnum[1]) && i == 8)) return FALSE;
  for (j = 2; j < i; j++) {
    if (!(isdigit(acnum[j]))) return FALSE;
  }
  return TRUE;
}

static void FindAccInDefCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  Char         acnum [17];
  BioseqPtr    bsp;
  Int2         i;
  CharPtr      p;
  Char         str [128];
  ValNodePtr   ttl;
  ValNodePtr   PNTR vnpp;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  if (! IS_Bioseq (sep)) return;
  vnpp = (ValNodePtr PNTR) mydata;
  if (vnpp == NULL) return;
  ttl = SeqEntryGetSeqDescr (sep, Seq_descr_title, NULL);
  while (ttl != NULL) {
    StringNCpy_0 (str, (CharPtr) ttl->data.ptrvalue, sizeof (str));
    TrimSpacesAroundString (str);
    if (! StringHasNoText (str)) {
      if (str [0] == 'a' && str [1] == 'c' && str [2] == 'c') {
        p = str + 3;
        for (i = 0; isalnum (*p) && *p != '\0'; p++, i++) {
          acnum [i] = *p;
        }
        acnum [i] = '\0';
        if (i == 6 || i == 8) {
          if (IsLegalAccession (acnum, i)) {
            bsp = (BioseqPtr) sep->data.ptrvalue;
            RecordAccToConvert (vnpp, bsp, str);
          }
        }
      }
    }
    ttl = SeqEntryGetSeqDescr (sep, Seq_descr_title, ttl);
  }
}

static void ConvertAccInDefToLocalIdProc (SeqEntryPtr sep)

{
  ConvAccPtr  cap;
  ValNodePtr  head;
  ValNodePtr  vnp;

  if (sep == NULL) return;
  head = NULL;
  SeqEntryExplore (sep, (Pointer) &head, FindAccInDefCallback);
  if (head == NULL) return;
  if (Message (MSG_YN, "Convert accLNNNNN or accLLNNNNNN in title to local ID?") == ANS_NO) return;
  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    cap = (ConvAccPtr) vnp->data.ptrvalue;
    if (cap != NULL) {
      SeqEntryExplore (sep, (Pointer) cap, CopyAccToLocalCallback);
    }
  }
  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    cap = (ConvAccPtr) vnp->data.ptrvalue;
    if (cap != NULL) {
      MemFree (cap->currID);
      MemFree (cap->newID);
      vnp->data.ptrvalue = MemFree (cap);
    }
  }
  ValNodeFree (head);
}

static Int2  taxonCount;

static Int4 DoSeqEntryToAsn3 (SeqEntryPtr sep, Boolean strip, Boolean correct,
                              Boolean force, Boolean dotaxon, MonitorPtr mon)

{
  BioseqSetPtr  bssp;
  SeqEntryPtr   oldscope;
  Int4          rsult;
  Char          str [32];

  rsult = 0;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && (bssp->_class == 7 || bssp->_class == 13 ||
                         bssp->_class == 14 || bssp->_class == 15)) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        rsult += DoSeqEntryToAsn3 (sep, strip, correct, force, dotaxon, mon);
      }
      return rsult;
    }
  }
/*#ifdef USE_TAXON*/
  if (dotaxon && mon != NULL) {
    taxonCount++;
    sprintf (str, "Processing Component %d", (int) taxonCount);
    MonitorStrValue (mon, str);
  }
/*#endif*/
/*#ifdef INTERNAL_NCBI_SEQUIN*/
  if (indexerVersion) {
    if ((! force) && (! NoBiosourceOrTaxonId (sep))) return 0;
  } else {
/*#else*/
    if ((! force) && SeqEntryGetSeqDescr (sep, Seq_descr_source, NULL) != NULL) return 0;
  }
/*#endif*/
  oldscope = SeqEntrySetScope (sep);
  if (dotaxon) {
    rsult = SeqEntryToAsn3Ex (sep, strip, correct, TRUE, GetTaxserverOrg, TaxMergeBSinDescr);
  } else {
    rsult = SeqEntryToAsn3Ex (sep, strip, correct, FALSE, NULL, NULL);
  }
  SeqEntrySetScope (oldscope);
  return rsult;
}

/* CheckSeqAlignCallback copied from salsa.c */
static void CheckSeqAlignCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)
{
  BioseqPtr          bsp;
  BioseqSetPtr       bssp;
  SeqAnnotPtr        sap,
                     pre;

  if (sep != NULL && sep->data.ptrvalue) {
     if (IS_Bioseq(sep)) {
        bsp = (BioseqPtr) sep->data.ptrvalue;
        if (bsp!=NULL) {
           pre=NULL;
           sap=bsp->annot;
           while (sap!= NULL)
           {
              if (sap->type == 2) {
                 if (is_dim1seqalign ((SeqAlignPtr) sap->data)) {
                    if (pre==NULL) {
                       bsp->annot=sap->next;
                       sap->next=NULL;
                       SeqAnnotFree (sap);
                       sap=bsp->annot; 
                    }
                    else {
                       pre->next=sap->next;
                       pre=sap->next;
                       sap->next=NULL; 
                       SeqAnnotFree (sap);
                       sap=pre;
                    }
                 }
                 else {
                    pre=sap;
                    sap=sap->next;
                 }
              }
              else {
                 pre=sap;
                 sap=sap->next;
              }
           }
        }
     }
     else if(IS_Bioseq_set(sep)) {
        bssp = (BioseqSetPtr)sep->data.ptrvalue;
        if (bssp!=NULL) {
           pre=NULL;
           sap=bssp->annot;
           while (sap!= NULL)
           {
              if (sap->type == 2) {
                 if (is_dim1seqalign ((SeqAlignPtr) sap->data)) {
                    if (pre==NULL) {
                       bssp->annot=sap->next;
                       sap->next=NULL;
                       SeqAnnotFree (sap);
                       sap=bssp->annot; 
                    }
                    else {
                       pre=sap->next;
                       sap->next=NULL; 
                       SeqAnnotFree (sap);
                       sap=sap->next;
                    }
                 }
                 else {
                    pre=sap;
                    sap=sap->next;
                 }
              }
              else {
                 pre=sap;
                 sap=sap->next;
              }
           }
        }
     }
  }
}

extern Int4 MySeqEntryToAsn3Ex (SeqEntryPtr sep, Boolean strip, Boolean correct, Boolean force, Boolean dotaxon);
extern Int4 MySeqEntryToAsn3Ex (SeqEntryPtr sep, Boolean strip, Boolean correct, Boolean force, Boolean dotaxon)

{
  MonitorPtr  mon;
  Int4        rsult;
  ErrSev      sev;

  rsult = 0;
  sev = ErrSetMessageLevel (SEV_FATAL);
  BasicSeqEntryCleanup (sep);
  EntryChangeImpFeat(sep);     /* change any CDS ImpFeat to real CdRegion */
  /* NormalizePeriodsOnInitials (sep); */ /* put periods on author initials */
  /* MoveRnaGBQualProductToName (sep); */ /* move rna gbqual product to rna-ref.ext.name */
  /* MoveProtGBQualProductToName (sep); */ /* move prot gbqual product to prot-ref.name */
  /* MoveCdsGBQualProductToName (sep); */ /* move cds gbqual product to prot-ref.name */
  /* MoveFeatGBQualsToFields (sep); */ /* move feature partial, exception to fields */
  if (indexerVersion) {
    StripTitleFromProtsInNucProts (sep);
    move_cds (sep); /* move CDS features to nuc-prot set */
  }
  ExtendGeneFeatIfOnMRNA (0, sep); /* gene on mRNA is full length */
  RemoveBioSourceOnPopSet (sep, NULL);
  /* SeqEntryExplore (sep, NULL, CleanupEmptyFeatCallback); */
  SeqEntryExplore (sep, NULL, DeleteMultipleTitles);
  SeqEntryPack (sep);
  if (indexerVersion) {
    ConvertAccInDefToLocalIdProc (sep); /* if title is accXNNNNN, convert to seqID */
  }
  if (useEntrez) {
    ValidateSeqAlignInSeqEntry (sep, FALSE, TRUE); /* remove components with seqID accXNNNNN */
  }
  SeqEntryExplore (sep, NULL, CheckSeqAlignCallback); /* remove alignments with single dimension */
  if ((! force) && (! NoBiosourceOrTaxonId (sep))) {
    /* EntryStripSerialNumber(sep); */ /* strip citation serial numbers */
    EntryChangeGBSource (sep);   /* at least remove redundant information in GBBlocks */
    EntryCheckGBBlock (sep);
    /* SeqEntryMoveDbxrefs (sep); */ /* db_xref gbqual to sfp->dbxref */
    EntryMergeDupBioSources (sep);
    GetRidOfEmptyFeatsDescStrings (0, sep);
    GetRidOfLocusInSeqIds (0, sep);
    ErrSetMessageLevel (sev);
    return rsult;
  }
  if (force && useTaxon) {
    dotaxon = TRUE;
  }
  mon = NULL;
  taxonCount = 0;
/*#ifdef USE_TAXON*/
  if (dotaxon) {
    WatchCursor ();
    mon = MonitorStrNewEx ("Taxonomy Lookup", 40, FALSE);
    MonitorStrValue (mon, "Connecting to Taxon");
    Update ();
    if (! TaxArchInit ()) {
      ArrowCursor ();
      Update ();
      Message (MSG_ERROR, "Couldn't connect to Taxon");
      MonitorFree (mon);
      Update ();
      ErrSetMessageLevel (sev);
      return 0;
    }
    MonitorStrValue (mon, "Processing Organism Info");
    Update ();
  }
/*#endif*/
  rsult = DoSeqEntryToAsn3 (sep, strip, correct, force, dotaxon, mon);
/*#ifdef USE_TAXON*/
  if (dotaxon) {
    MonitorStrValue (mon, "Closing Taxon");
    Update ();
    TaxArchFini ();
    MonitorFree (mon);
    ArrowCursor ();
    Update ();
  }
/*#endif*/
  /* EntryStripSerialNumber(sep); */ /* strip citation serial numbers */
  EntryChangeGBSource (sep);   /* remove redundant information in GBBlocks again */
  EntryCheckGBBlock (sep);
  /* SeqEntryMoveDbxrefs (sep); */ /* db_xref gbqual to sfp->dbxref */
  EntryMergeDupBioSources (sep);
  GetRidOfEmptyFeatsDescStrings (0, sep);
  GetRidOfLocusInSeqIds (0, sep);
  ErrSetMessageLevel (sev);
  ErrClear ();
  ErrShow ();
  return rsult;
}

typedef struct partialformdata {
  FEATURE_FORM_BLOCK

  ValNodePtr         featlist;
  LisT               feature;
  PopuP              change5;
  PopuP              change3;
  GrouP              nucprotchoice;
  Int2               subtype;
  Int2               leftpolicy;
  Int2               rightpolicy;
  Int2               nucprotpolicy;
  ObjMgrPtr          omp;
  ObjMgrTypePtr      omtp;
} PartialFormData, PNTR PartialFormPtr;

static void PartialCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  Boolean         atEnd;
  BioseqPtr       bsp;
  BioseqSetPtr    bssp;
  SeqLocPtr       firstSlp;
  SeqLocPtr       lastSlp;
  ObjMgrTypePtr   omtp;
  Boolean         partial5;
  Boolean         partial3;
  PartialFormPtr  pfp;
  SeqAnnotPtr     sap;
  SeqFeatPtr      sfp;
  SeqLocPtr       slp;
  Uint1           strand;
  Uint2           subtype;

  if (mydata == NULL) return;
  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  pfp = (PartialFormPtr) mydata;
  if (pfp == NULL) return;
  omtp = pfp->omtp;
  if (omtp == NULL || omtp->subtypefunc == NULL) return;
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
        subtype = (*(omtp->subtypefunc)) ((Pointer) sfp);
        if (pfp->subtype == 0 || subtype == pfp->subtype) {
          CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
          firstSlp = NULL;
          lastSlp = NULL;
          slp = SeqLocFindNext (sfp->location, NULL);
          while (slp != NULL) {
            if (firstSlp == NULL) {
              firstSlp = slp;
            }
            lastSlp = slp;
            slp = SeqLocFindNext (sfp->location, slp);
          }
          if (firstSlp != NULL) {
            bsp = BioseqFind (SeqLocId (firstSlp));
            if (bsp != NULL) {
              if (((pfp->nucprotpolicy == 1 || pfp->nucprotpolicy == 2) && ISA_na (bsp->mol)) ||
                  ((pfp->nucprotpolicy == 1 || pfp->nucprotpolicy == 3) && ISA_aa (bsp->mol))) {
                strand = SeqLocStrand (firstSlp);
                atEnd = FALSE;
                if (strand == Seq_strand_minus) {
                  if (GetOffsetInBioseq (firstSlp, bsp, SEQLOC_START) == bsp->length - 1) {
                    atEnd = TRUE;
                  }
                } else {
                  if (GetOffsetInBioseq (firstSlp, bsp, SEQLOC_START) == 0) {
                    atEnd = TRUE;
                  }
                }
                switch (pfp->leftpolicy) {
                  case 1 :
                    partial5 = TRUE;
                    break;
                  case 2 :
                    if (atEnd) {
                      partial5 = TRUE;
                    }
                    break;
                  case 3 :
                    partial5 = FALSE;
                   break;
                  case 4 :
                    if (! atEnd) {
                      partial5 = FALSE;
                    }
                    break;
                  default :
                    break;
                }
              }
            }
          }
          if (lastSlp != NULL) {
            bsp = BioseqFind (SeqLocId (firstSlp));
            if (bsp != NULL) {
              if (((pfp->nucprotpolicy == 1 || pfp->nucprotpolicy == 2) && ISA_na (bsp->mol)) ||
                  ((pfp->nucprotpolicy == 1 || pfp->nucprotpolicy == 3) && ISA_aa (bsp->mol))) {
                strand = SeqLocStrand (firstSlp);
                atEnd = FALSE;
                if (strand == Seq_strand_minus) {
                  if (GetOffsetInBioseq (lastSlp, bsp, SEQLOC_STOP) == 0) {
                    atEnd = TRUE;
                  }
                } else {
                  if (GetOffsetInBioseq (lastSlp, bsp, SEQLOC_STOP) == bsp->length - 1) {
                    atEnd = TRUE;
                  }
                }
                switch (pfp->rightpolicy) {
                  case 1 :
                    partial3 = TRUE;
                    break;
                  case 2 :
                    if (atEnd) {
                      partial3 = TRUE;
                    }
                    break;
                  case 3 :
                    partial3 = FALSE;
                    break;
                  case 4 :
                    if (! atEnd) {
                      partial3 = FALSE;
                    }
                    break;
                  default :
                    break;
                }
              }
            }
          }
          SetSeqLocPartial (sfp->location, partial5, partial3);
          sfp->partial = (sfp->partial || partial5 || partial3 || LocationHasNullsBetween (sfp->location));
        }
        sfp = sfp->next;
      }
    }
    sap = sap->next;
  }
}

static void DoEditPartials (ButtoN b)

{
  PartialFormPtr  pfp;
  SeqEntryPtr     sep;
  Int2            val;
  ValNodePtr      vnp;

  pfp = (PartialFormPtr) GetObjectExtra (b);
  if (pfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (pfp->input_entityID);
  if (sep == NULL) return;
  Hide (pfp->form);
  WatchCursor ();
  Update ();

  pfp->leftpolicy = GetValue (pfp->change5);
  pfp->rightpolicy = GetValue (pfp->change3);
  pfp->nucprotpolicy = GetValue (pfp->nucprotchoice);

  val = GetValue (pfp->feature);
  vnp = NULL;
  if (val > 0) {
    vnp = pfp->featlist;
    while (vnp != NULL && val > 1) {
      val--;
      vnp = vnp->next;
    }
  }
  if (vnp != NULL) {
    pfp->subtype = vnp->choice;
  } else {
    pfp->subtype = 0;
  }

  pfp->omp = ObjMgrGet ();
  pfp->omtp = NULL;
  if (pfp->omp != NULL) {
    pfp->omtp = ObjMgrTypeFind (pfp->omp, OBJ_SEQFEAT, NULL, NULL);
  }
  SeqEntryExplore (sep, (Pointer) pfp, PartialCallback);

  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (pfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, pfp->input_entityID, 0, 0);
  Remove (pfp->form);
}

static void PartialMessageProc (ForM f, Int2 mssg)

{
  PartialFormPtr  pfp;

  pfp = (PartialFormPtr) GetObjectExtra (f);
  if (pfp != NULL) {
    if (pfp->appmessage != NULL) {
      pfp->appmessage (f, mssg);
    }
  }
}

static void CleanupPartialPage (GraphiC g, VoidPtr data)

{
  PartialFormPtr  pfp;

  pfp = (PartialFormPtr) data;
  if (pfp != NULL) {
    ValNodeFreeData (pfp->featlist);
  }
  StdCleanupFormProc (g, data);
}

extern void EditFeaturePartials (IteM i)

{
  BaseFormPtr        bfp;
  ButtoN             b;
  GrouP              c;
  FeatDefPtr         curr;
  GrouP              g;
  GrouP              h;
  ValNodePtr         head;
  GrouP              k;
  Uint1              key;
  CharPtr            label = NULL;
  Int2               listHeight;
  PartialFormPtr     pfp;
  SeqEntryPtr        sep;
  StdEditorProcsPtr  sepp;
  Uint2              subtype;
  ValNodePtr         vnp;
  WindoW             w;
  GrouP              x;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  pfp = (PartialFormPtr) MemNew (sizeof (PartialFormData));
  if (pfp == NULL) return;
  w = FixedWindow (-50, -33, -10, -10, "Feature Partial Editor", StdCloseWindowProc);
  SetObjectExtra (w, pfp, CleanupPartialPage);
  pfp->form = (ForM) w;
  pfp->formmessage = PartialMessageProc;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    pfp->appmessage = sepp->handleMessages;
  }

  pfp->input_entityID = bfp->input_entityID;
  pfp->input_itemID = bfp->input_itemID;
  pfp->input_itemtype = bfp->input_itemtype;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  g = HiddenGroup (h, 0, 2, NULL);
  StaticPrompt (g, "Feature", 0, 0, programFont, 'c');
  if (indexerVersion) {
    listHeight = 16;
  } else {
    listHeight = 8;
  }
  pfp->feature = SingleList (g, 16, listHeight, NULL);

  head = ValNodeNew (NULL);
  if (head != NULL) {
    head->choice = 0;
    head->data.ptrvalue = StringSave ("[ALL FEATURES]");
  }
  curr = FeatDefFindNext (NULL, &key, &label, FEATDEF_ANY, TRUE);
  while (curr != NULL) {
    if (key != FEATDEF_BAD) {
      subtype = curr->featdef_key;
      if (subtype != FEATDEF_misc_RNA &&
          subtype != FEATDEF_precursor_RNA &&
          subtype != FEATDEF_mat_peptide &&
          subtype != FEATDEF_sig_peptide &&
          subtype != FEATDEF_transit_peptide) {
        vnp = ValNodeNew (head);
        if (head == NULL) {
          head = vnp;
        }
        if (vnp != NULL) {
          vnp->choice = subtype;
          vnp->data.ptrvalue = StringSave (curr->typelabel);
        }
      }
    }
    curr = FeatDefFindNext (curr, &key, &label, FEATDEF_ANY, TRUE);
  }
  if (head != NULL) {
    head = SortValNode (head, SortByVnpChoice);
    for (vnp = head; vnp != NULL; vnp = vnp->next) {
      ListItem (pfp->feature, (CharPtr) vnp->data.ptrvalue);
    }
  }
  pfp->featlist = head;

  k = HiddenGroup (h, 2, 0, NULL);

  StaticPrompt (k, "5' policy", 0, stdLineHeight, programFont, 'l');
  pfp->change5 = PopupList (k, TRUE, NULL);
  PopupItem (pfp->change5, "Set");
  PopupItem (pfp->change5, "Set only if at 5' end");
  PopupItem (pfp->change5, "Clear");
  PopupItem (pfp->change5, "Clear if not at 5' end");
  PopupItem (pfp->change5, "Do not change");
  SetValue (pfp->change5, 5);

  StaticPrompt (k, "3' policy", 0, stdLineHeight, programFont, 'l');
  pfp->change3 = PopupList (k, TRUE, NULL);
  PopupItem (pfp->change3, "Set");
  PopupItem (pfp->change3, "Set only if at 3' end");
  PopupItem (pfp->change3, "Clear");
  PopupItem (pfp->change3, "Clear if not at 3' end");
  PopupItem (pfp->change3, "Do not change");
  SetValue (pfp->change3, 5);

  x = HiddenGroup (h, 2, 0, NULL);
  /* StaticPrompt (x, "Constraints", 0, stdLineHeight, programFont, 'l'); */
  pfp->nucprotchoice = HiddenGroup (x, 4, 0, NULL);
  RadioButton (pfp->nucprotchoice, "All Sequences");
  RadioButton (pfp->nucprotchoice, "Nucleotides");
  RadioButton (pfp->nucprotchoice, "Proteins");
  SetValue (pfp->nucprotchoice, 1);

  c = HiddenGroup (h, 4, 0, NULL);
  b = PushButton (c, "Accept", DoEditPartials);
  SetObjectExtra (b, pfp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) k, (HANDLE) x, (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Update ();
}

static Boolean AddOrgToDefGatherFunc (GatherContextPtr gcp)

{
  CharPtr     def;
  CharPtr     ptr;
  ValNodePtr  sdp;
  CharPtr     str;
  CharPtr     text;

  if (gcp == NULL || gcp->thisitem == NULL) return TRUE;
  if (gcp->thistype != OBJ_SEQDESC) return TRUE;
  text = (CharPtr) gcp->userdata;
  if (text == NULL || StringHasNoText (text)) return TRUE;
  sdp = (ValNodePtr) gcp->thisitem;
  if (sdp->choice != Seq_descr_title) return TRUE;
  def = (CharPtr) sdp->data.ptrvalue;
  if (StringHasNoText (def)) return TRUE;

  ptr = StringISearch (def, text);
  if (ptr != NULL && ptr == def) return TRUE;
  str = MemNew ((StringLen (text) + StringLen (def) + 4) * sizeof (Char));
  if (str != NULL) {
    StringCpy (str, text);
    StringCat (str, " ");
    StringCat (str, def);
    sdp->data.ptrvalue = MemFree (sdp->data.ptrvalue);
    sdp->data.ptrvalue = str;
    ObjMgrSetDirtyFlag (gcp->entityID, TRUE);
  }
  return TRUE;
}

static void AddOrgToDef (Uint2 entityID, SeqEntryPtr sep, Int2 orgmod, Int2 subsource)

{
  EnumFieldAssocPtr  ap;
  BioSourcePtr       biop;
  BioseqSetPtr       bssp;
  Char               ch;
  GatherScope        gs;
  OrgModPtr          mod;
  OrgNamePtr         onp;
  OrgRefPtr          orp;
  CharPtr            ptr;
  SubSourcePtr       ssp;
  Char               str [96];
  Char               text [64];

  if (sep == NULL) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && (bssp->_class == 7 || bssp->_class == 13 ||
                         bssp->_class == 14 || bssp->_class == 15)) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        AddOrgToDef (entityID, sep, orgmod, subsource);
      }
      return;
    }
  }
  biop = NULL;
  text [0] = '\0';
  str [0] = '\0';
  SeqEntryToBioSource (sep, NULL, str, sizeof (str) - 1, &biop);
  if (orgmod == 0 && subsource == 0) {
    /*
    ptr = StringSearch (str, "(");
    if (ptr != NULL) {
      *ptr = '\0';
    }
    */
    TrimSpacesAroundString (str);
    if (StringICmp (str, "Human immunodeficiency virus type 1") == 0) {
      StringCpy (str, "HIV-1");
    } else if (StringICmp (str, "Human immunodeficiency virus type 2") == 0) {
      StringCpy (str, "HIV-2");
    }
    str [0] = TO_UPPER (str [0]);
  } else if (biop != NULL && biop->org != NULL) {
    text [0] = '\0';
    str [0] = '\0';
    orp = biop->org;
    if (orgmod > 0) {
      onp = orp->orgname;
      if (onp != NULL) {
        mod = onp->mod;
        while (mod != NULL) {
          if (mod->subtype == orgmod) {
            StringNCpy_0 (text, mod->subname, sizeof (text));
            for (ap = orgmod_subtype_alist; ap->name != NULL; ap++) {
              if (ap->value == orgmod) {
                StringNCpy_0 (str, ap->name, sizeof (str));
              }
            }
          }
          mod = mod->next;
        }
      }
    } else if (subsource > 0) {
      ssp = biop->subtype;
      while (ssp != NULL) {
        if (ssp->subtype == subsource) {
          StringNCpy_0 (text, ssp->name, sizeof (text));
          for (ap = subsource_subtype_alist; ap->name != NULL; ap++) {
            if (ap->value == subsource) {
              StringNCpy_0 (str, ap->name, sizeof (str));
            }
          }
        }
        ssp = ssp->next;
      }
    }
    if (StringHasNoText (text)) {
      str [0] = '\0';
      text [0] = '\0';
    } else {
      StringCat (str, " ");
      ptr = str;
      while (*ptr != '\0') {
        ch = *ptr;
        *ptr = TO_LOWER (ch);
        ptr++;
      }
      StringCat (str, text);
    }
  }
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = FALSE;
  MemSet ((Pointer) (gs.ignore), (int)(TRUE), (size_t) (OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  gs.ignore[OBJ_SEQDESC] = FALSE;
  gs.scope = sep;
  GatherSeqEntry (sep, (Pointer) str, AddOrgToDefGatherFunc, &gs);
}

extern void CommonAddOrgOrModsToDefLines (IteM i, Int2 orgmod, Int2 subsource, ButtoN b)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

  if (b != NULL) {
    bfp = GetObjectExtra (b);
  } else {
#ifdef WIN_MAC
    bfp = currentFormDataPtr;
#else
    bfp = GetObjectExtra (i);
#endif
  }
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  WatchCursor ();
  Update ();
  AddOrgToDef (bfp->input_entityID, sep, orgmod, subsource);
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static void RemoveAlignmentCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  SeqAlignPtr    nextsalp;
  SeqAnnotPtr    nextsap;
  Pointer PNTR   prevsalp;
  Pointer PNTR   prevsap;
  SeqAlignPtr    salp;
  SeqAnnotPtr    sap;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
    prevsap = (Pointer PNTR) &(bsp->annot);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
    prevsap = (Pointer PNTR) &(bssp->annot);
  } else return;
  while (sap != NULL) {
    nextsap = sap->next;
    if (sap->type == 2) {
      salp = (SeqAlignPtr) sap->data;
      prevsalp = (Pointer PNTR) &(sap->data);
      while (salp != NULL) {
        nextsalp = salp->next;
        if (salp->type >= 0 && salp->type <= 255) {
          *(prevsalp) = salp->next;
          salp->next = NULL;
          SeqAlignFree (salp);
        } else {
          prevsalp = (Pointer PNTR) &(salp->next);
        }
        salp = nextsalp;
      }
    }
    if (sap->data == NULL) {
      *(prevsap) = sap->next;
      sap->next = NULL;
      SeqAnnotFree (sap);
    } else {
      prevsap = (Pointer PNTR) &(sap->next);
    }
    sap = nextsap;
  }
}

extern void RemoveAlignment (IteM i)

{
  BaseFormPtr bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  SeqEntryExplore (sep, NULL, RemoveAlignmentCallback);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ObjMgrDeSelect (0, 0, 0, 0, NULL);
  Update ();
}

static void RemoveGraphCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  SeqAnnotPtr    nextsap;
  SeqGraphPtr    nextsgp;
  Pointer PNTR   prevsap;
  Pointer PNTR   prevsgp;
  SeqAnnotPtr    sap;
  SeqGraphPtr    sgp;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
    prevsap = (Pointer PNTR) &(bsp->annot);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
    prevsap = (Pointer PNTR) &(bssp->annot);
  } else return;
  while (sap != NULL) {
    nextsap = sap->next;
    if (sap->type == 3) {
      sgp = (SeqGraphPtr) sap->data;
      prevsgp = (Pointer PNTR) &(sap->data);
      while (sgp != NULL) {
        nextsgp = sgp->next;
        if (sgp->flags [2] >= 1 && sgp->flags [2] <= 3) {
          *(prevsgp) = sgp->next;
          sgp->next = NULL;
          SeqGraphFree (sgp);
        } else {
          prevsgp = (Pointer PNTR) &(sgp->next);
        }
        sgp = nextsgp;
      }
    }
    if (sap->data == NULL) {
      *(prevsap) = sap->next;
      sap->next = NULL;
      SeqAnnotFree (sap);
    } else {
      prevsap = (Pointer PNTR) &(sap->next);
    }
    sap = nextsap;
  }
}

extern void RemoveGraph (IteM i)

{
  BaseFormPtr bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  SeqEntryExplore (sep, NULL, RemoveGraphCallback);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ObjMgrDeSelect (0, 0, 0, 0, NULL);
  Update ();
}

static void MarkProteinCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr        bsp;
  ValNodePtr PNTR  vnpp;

  if (mydata == NULL) return;
  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  vnpp = (ValNodePtr PNTR) mydata;
  if (vnpp == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (ISA_aa (bsp->mol)) {
      ValNodeAddPointer (vnpp, 0, (Pointer) bsp);
    }
  }
}

extern void RemoveProteins (IteM i)

{
  BaseFormPtr    bfp;
  BioseqPtr      bsp;
  Uint2          itemID;
  ObjMgrDataPtr  omdptop;
  ObjMgrData     omdata;
  OMProcControl  ompc;
  Uint2          parenttype;
  Pointer        parentptr;
  SeqEntryPtr    sep;
  ValNodePtr     tmp;
  ValNodePtr     vnp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  SaveSeqEntryObjMgrData (sep, &omdptop, &omdata);
  GetSeqEntryParent (sep, &parentptr, &parenttype);
  vnp = NULL;
  SeqEntryExplore (sep, (Pointer) &vnp, MarkProteinCallback);
  for (tmp = vnp; tmp != NULL; tmp = tmp->next) {
    bsp = (BioseqPtr) tmp->data.ptrvalue;
    itemID = GetItemIDGivenPointer (bfp->input_entityID, OBJ_BIOSEQ, (Pointer) bsp);
    if (itemID > 0) {
      MemSet ((Pointer) (&ompc), 0, sizeof (OMProcControl));
      ompc.do_not_reload_from_cache = TRUE;
      ompc.input_entityID = bfp->input_entityID;
      ompc.input_itemID = itemID;
      ompc.input_itemtype = OBJ_BIOSEQ;
      if (! DetachDataForProc (&ompc, FALSE)) {
        Message (MSG_POSTERR, "DetachDataForProc failed");
      }
    }
  }
  ValNodeFree (vnp);
  SeqMgrLinkSeqEntry (sep, parenttype, parentptr);
  RestoreSeqEntryObjMgrData (sep, omdptop, &omdata);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  ObjMgrDeSelect (0, 0, 0, 0, NULL);
  Update ();
}
