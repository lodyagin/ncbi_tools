/*   cdrgn.c
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
* File Name:  cdrgn.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   1/22/95
*
* $Revision: 6.19 $
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

#include <cdrgn.h>
#include <biosrc.h>
#include <seqport.h>
#include <gather.h>
#include <tofasta.h>
#include <document.h>
#include <suggslp.h>
#include <salutil.h>
#include <edutil.h>
#include <bspview.h>
#include <toasn3.h>

#define CODING_REGION_PAGE    0
#define GENE_PAGE             0
#define PROT_PAGE             0
#define RNA_PAGE              0
#define COMMON_PAGE           1
#define LOCATION_PAGE         2

#define NUM_PAGES  8

typedef struct cdrgnform {
  FEATURE_FORM_BLOCK
  SeqEntryPtr   sep;
  GrouP         pages [NUM_PAGES];
  DialoG        foldertabs;
  Int2          currentPage;
  ButtoN        conceptTransA;
  GrouP         protTextGrp;
  TexT          protNameText;
  TexT          protDescText;
  GrouP         protPromptGrp;
  DialoG        ecProcs;
  Uint2         protEntityID;
  Uint2         protItemID;
  Uint2         protItemtype;
  Boolean       protFound;
  DialoG        protseq;
  GrouP         protSeqIdGrp;
  TexT          protSeqIdTxt;
  ButtoN        launchBtn;
  ButtoN        edProtBtn;
  ButtoN        edSeqBtn;
  BioseqPtr     usethisbioseq;
  Boolean       autoUpdatePartials;
  Boolean       autoRetranslate;
  Boolean       locvisited;
  ButtoN        convertToMiscFeat;
  Boolean       saveAsMiscFeat;
} CdRgnForm, PNTR CdRgnFormPtr;

typedef struct cdrgnpage {
  DIALOG_MESSAGE_BLOCK
  ButtoN             orf;
  PopuP              frame;
  ButtoN             conflict;
  PopuP              geneticCode;
  DialoG             cdBrk;
  Int2               count;
  EnumFieldAssoc     PNTR alist;
  EnumFieldAssocPtr  alists [2];
  GrouP              cdRgnGrp [4];
  GrouP              commonGrp;
  Int2               currentPage;
  CdRgnFormPtr       cfp;
} CdRgnPage, PNTR CdRgnPagePtr;

typedef struct rawbioseqpage {
  DIALOG_MESSAGE_BLOCK
  ByteStorePtr  bs;
  Uint1         mol;
  Uint1         repr;
  Uint1         seq_data_type;
  DoC           doc;
  CdRgnFormPtr  cfp;
} RawBioseqPage, PNTR RawBioseqPagePtr;

static ParData brParFmt = {FALSE, FALSE, FALSE, FALSE, FALSE, 0, 0};
static ColData brColFmt = {0, 0, 80, 0, NULL, 'l', FALSE, FALSE, FALSE, FALSE, TRUE};

static Boolean RawFastaFileFunc (BioseqPtr bsp, Int2 key, CharPtr buf, Uint4 buflen, Pointer data)

{
  FILE  *fp;

  fp = (FILE *) data;
  switch (key) {
    case FASTA_SEQLINE :
      fprintf (fp, "%s\n", buf);
      break;
    default :
      break;
  }
  return TRUE;
}

static void BioseqPtrToRawBioseqPage (DialoG d, Pointer data)

{
  BioseqPtr         bsp;
  Char              buf [256];
  FILE              *fp;
  MyFsa             mfa;
  ErrSev            oldErrSev;
  Char              path [PATH_MAX];
  RecT              r;
  RawBioseqPagePtr  rbpp;

  rbpp = (RawBioseqPagePtr) GetObjectExtra (d);
  bsp = (BioseqPtr) data;
  if (rbpp != NULL) {
    rbpp->bs = BSFree (rbpp->bs);
    Reset (rbpp->doc);
    if (bsp != NULL) {
      if (bsp->repr == Seq_repr_raw || bsp->repr == Seq_repr_const) {
        rbpp->bs = BSDup (bsp->seq_data);
        rbpp->seq_data_type = bsp->seq_data_type;
        rbpp->mol = bsp->mol;
        rbpp->repr = bsp->repr;
        TmpNam (path);
        fp = FileOpen (path, "w");
        if (fp != NULL) {
          MemSet ((Pointer) (&mfa), 0, sizeof (MyFsa));
          mfa.buf = buf;
          mfa.buflen = 255;
          mfa.seqlen = 50;
          mfa.mydata = (Pointer) fp;
          mfa.myfunc = RawFastaFileFunc;
          mfa.bad_asn1 = FALSE;
          mfa.order = 0;
          mfa.accession = NULL;
          mfa.organism = NULL;
          mfa.tech = 0;
          mfa.do_virtual = FALSE;
	  mfa.formatdb	= TRUE;

	  if (ISA_na(bsp->mol))
	      mfa.code = Seq_code_iupacna;
	  else
	      mfa.code = Seq_code_ncbieaa;

          oldErrSev = ErrSetMessageLevel (SEV_MAX);
          BioseqToFastaX (bsp, &mfa, ISA_na (bsp->mol));
          ErrSetMessageLevel (oldErrSev);
        }
        FileClose (fp);
        ObjectRect (rbpp->doc, &r);
        InsetRect (&r, 4, 4);
        brColFmt.pixWidth = r.right - r.left;
        DisplayFancy (rbpp->doc, path, &brParFmt, &brColFmt, programFont, 4);
        FileRemove (path);
      }
    }
    Update ();
  }
}

static Pointer RawBioseqPageToBioseqPtr (DialoG d)

{
  BioseqPtr         bsp;
  Char              ch;
  CdRgnFormPtr      cfp;
  Int4              cnt;
  Int4              cntr;
  Int2              i;
  Char              ids [32];
  Int4              nxt;
  Char              str [102];
  RawBioseqPagePtr  rbpp;

  bsp = NULL;
  rbpp = (RawBioseqPagePtr) GetObjectExtra (d);
  if (rbpp != NULL) {
    bsp = BioseqNew ();
    if (bsp != NULL) {
      cfp = rbpp->cfp;
      ids [0] = '\0';
      if (cfp != NULL) {
        if (Visible (cfp->protSeqIdGrp)) {
          GetTitle (cfp->protSeqIdTxt, ids, sizeof (ids) - 1);
        }
      }
      if (StringHasNoText (ids)) {
        if (ISA_na (rbpp->mol)) {
          bsp->id = MakeSeqID ("dna_seq");
        } else if (ISA_aa (rbpp->mol)) {
          bsp->id = MakeSeqID ("prt_seq");
        } else {
          bsp->id = MakeSeqID ("oth_seq");
        }
      } else {
        bsp->id = MakeSeqID (ids);
      }
      bsp->mol = rbpp->mol;
      bsp->repr = rbpp->repr;
      /*
      bsp->seq_data_type = ISA_na (bsp->mol) ? Seq_code_iupacna : Seq_code_ncbieaa;
      */
      bsp->seq_data_type = rbpp->seq_data_type;
      bsp->length = BSLen (rbpp->bs);
      if (bsp->length == 0) {
        bsp = BioseqFree (bsp);
        return NULL;
      }
      bsp->seq_data = BSNew (bsp->length);
      if (bsp->seq_data != NULL) {
        BSSeek (rbpp->bs, 0, SEEK_SET);
        cntr = BSLen (rbpp->bs);
        cnt = MIN (cntr, 100L);
        while (cnt > 0) {
          cntr -= cnt;
          nxt = MIN (cntr, 100L);
          BSRead (rbpp->bs, str, cnt * sizeof (Char));
          for (i = 0; i < cnt; i++) {
            ch = str [i];
            str [i] = TO_UPPER (ch);
          }
          if (nxt == 0 && cnt > 0 && bsp->seq_data_type == Seq_code_ncbieaa) {
            if (str [cnt - 1] == '*') {
              cnt--;
              if (bsp->length > 0) {
                (bsp->length)--;
              }
            }
          }
          BSWrite (bsp->seq_data, str, cnt * sizeof (Char));
          cnt = nxt;
        }
      }
      /*
      if (bsp->seq_data_type == Seq_code_ncbieaa) {
        BioseqRawConvert (bsp, Seq_code_iupacaa);
      }
      */
      SeqMgrAddToBioseqIndex (bsp);
    }
  }
  return (Pointer) bsp;
}

static void CleanupRawBioseqPage (GraphiC g, VoidPtr data)

{
  RawBioseqPagePtr  rbpp;

  rbpp = (RawBioseqPagePtr) data;
  if (rbpp != NULL) {
    BSFree (rbpp->bs);
  }
  MemFree (data);
}

static DialoG CreateBioseqRawDialog (GrouP h, CharPtr title, Uint1 mol, CdRgnFormPtr cfp)

{
  Int2              charwid;
  RawBioseqPagePtr  rbpp;
  GrouP             m;
  GrouP             p;
  GrouP             s;

  p = HiddenGroup (h, 1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  rbpp = (RawBioseqPagePtr) MemNew (sizeof (RawBioseqPage));
  if (rbpp != NULL) {

    SetObjectExtra (p, rbpp, CleanupRawBioseqPage);
    rbpp->dialog = (DialoG) p;
    rbpp->todialog = BioseqPtrToRawBioseqPage;
    rbpp->fromdialog = RawBioseqPageToBioseqPtr;
    rbpp->testdialog = NULL;

    if (title != NULL && title [0] != '\0') {
      s = NormalGroup (p, 0, -2, title, systemFont, NULL);
    } else {
      s = HiddenGroup (p, 0, -2, NULL);
    }
    m = HiddenGroup (s, -1, 0, NULL);
    SetGroupSpacing (m, 10, 10);

    rbpp->cfp = cfp;
    rbpp->mol = mol;
    SelectFont (programFont);
    charwid = CharWidth ('0');
    SelectFont (systemFont);
    rbpp->doc = DocumentPanel (m, charwid * 50 + 2, stdLineHeight * 6);
  }

  return (DialoG) p;
}

static void ConvertToAACodeBreak (CdRgnFormPtr cfp, CdRegionPtr crp, CodeBreakPtr cb)

{
  SeqFeatPtr  sfp;
  SeqLocPtr   slp;

  if (cfp == NULL || crp == NULL || cb == NULL) return;
  sfp = SeqFeatNew ();
  if (sfp == NULL) return;
  sfp->data.choice = SEQFEAT_CDREGION;
  sfp->data.value.ptrvalue = crp;
  sfp->location = DialogToPointer (cfp->location);
  sfp->product = DialogToPointer (cfp->product);
  while (cb != NULL) {
    slp = dnaLoc_to_aaLoc (sfp, cb->loc, TRUE, NULL);
    cb->loc = SeqLocFree (cb->loc);
    cb->loc = slp;
    cb = cb->next;
  }
  SeqLocFree (sfp->product);
  SeqLocFree (sfp->location);
  MemFree (sfp);
}

static void ConvertToNTCodeBreak (CdRgnFormPtr cfp, CdRegionPtr crp, CodeBreakPtr cb)

{
  SeqFeatPtr  sfp;
  SeqLocPtr   slp;

  if (cfp == NULL || crp == NULL || cb == NULL) return;
  sfp = SeqFeatNew ();
  if (sfp == NULL) return;
  sfp->data.choice = SEQFEAT_CDREGION;
  sfp->data.value.ptrvalue = crp;
  sfp->location = DialogToPointer (cfp->location);
  sfp->product = DialogToPointer (cfp->product);
  while (cb != NULL) {
    slp = aaLoc_to_dnaLoc (sfp, cb->loc);
    cb->loc = SeqLocFree (cb->loc);
    cb->loc = slp;
    cb = cb->next;
  }
  SeqLocFree (sfp->product);
  SeqLocFree (sfp->location);
  MemFree (sfp);
}

static void CdRegionPtrToCdRgnPage (DialoG d, Pointer data)

{
  CodeBreakPtr    cbhead;
  CodeBreakPtr    cbp;
  CdRgnPagePtr    cpp;
  CdRegionPtr     crp;
  CodeBreakPtr    curr;
  GeneticCodePtr  gcp;
  Int4            id;
  CodeBreakPtr    last;
  Char            name [256];
  ValNodePtr      tmp;

  cpp = (CdRgnPagePtr) GetObjectExtra (d);
  crp = (CdRegionPtr) data;
  if (cpp != NULL) {
    if (crp != NULL) {
      SetStatus (cpp->orf, crp->orf);
      if (crp->orf) {
        SafeShow (cpp->orf);
      }
      SetStatus (cpp->conflict, crp->conflict);
      if (crp->conflict) {
        SafeShow (cpp->conflict);
      }
      SetValue (cpp->frame, crp->frame + 1);
      id = 0;
      name [0] = '\0';
      gcp = crp->genetic_code;
      if (gcp != NULL) {
        tmp = (ValNodePtr) gcp->data.ptrvalue;
        for (tmp = (ValNodePtr) gcp->data.ptrvalue; tmp != NULL; tmp = tmp->next) {
          switch (tmp->choice) {
            case 1 :
              if (name [0] == '\0') {
                StringNCpy_0 (name, (CharPtr) tmp->data.ptrvalue, sizeof (name));
              }
              break;
            case 2 :
              id = tmp->data.intvalue;
              break;
            default :
              break;
          }
        }
        if (id == 0) {
          gcp = GeneticCodeFind (id, name);
          if (gcp != NULL) {
            for (tmp = (ValNodePtr) gcp->data.ptrvalue; tmp != NULL; tmp = tmp->next) {
              switch (tmp->choice) {
                case 2 :
                  id = tmp->data.intvalue;
                  break;
                default :
                  break;
              }
            }
          }
        }
      }
      if (id > 0 && id < NUM_GENETIC_CODES) {
        SetValue (cpp->geneticCode, gcIdToIndex [(Int2) id]);
      } else {
        SetValue (cpp->geneticCode, 1);
      }
      cbhead = NULL;
      last = NULL;
      curr = crp->code_break;
      while (curr != NULL) {
        cbp = AsnIoMemCopy ((Pointer) curr,
                            (AsnReadFunc) CodeBreakAsnRead,
                            (AsnWriteFunc) CodeBreakAsnWrite);
        if (cbhead == NULL) {
          cbhead = cbp;
        }
        if (last != NULL) {
          last->next = cbp;
        }
        last = cbp;
        curr = curr->next;
      }
      ConvertToAACodeBreak (cpp->cfp, crp, cbhead);
      PointerToDialog (cpp->cdBrk, cbhead);
      CodeBreakFree (cbhead);
    } else {
      SafeSetStatus (cpp->orf, FALSE);
      SafeSetStatus (cpp->conflict, FALSE);
      SafeSetValue (cpp->frame, 1);
      SafeSetValue (cpp->geneticCode, 1);
    }
  }
}

static SeqIdPtr  globalid;

static Pointer CdRgnPageToCdRegionPtr (DialoG d)

{
  BioseqPtr     bsp;
  CodeBreakPtr  cbp;
  CdRgnFormPtr  cfp;
  ValNodePtr    code;
  CdRgnPagePtr  cpp;
  CdRegionPtr   crp;
  ImpFeatPtr    ifp;
  SeqIdPtr      sip;
  SeqLocPtr     slp;
  Int2          val;
  ValNodePtr    vnp;

  crp = NULL;
  cpp = (CdRgnPagePtr) GetObjectExtra (d);
  if (cpp != NULL) {
    cfp = cpp->cfp;
    if (cfp != NULL && cfp->saveAsMiscFeat) {
      ifp = ImpFeatNew ();
      if (ifp != NULL) {
        ifp->key = StringSave ("misc_feature");
      }
      return (Pointer) ifp;
    }
    crp = CdRegionNew ();
    if (crp != NULL) {
      crp->orf = GetStatus (cpp->orf);
      crp->conflict = GetStatus (cpp->conflict);
      crp->frame = (Uint1) GetValue (cpp->frame) - 1;
      crp->gaps = 0;
      crp->mismatch = 0;
      crp->stops = 0;
      val = gcIndexToId [GetValue (cpp->geneticCode)];
      if (val == 0) {
        val = 1;  /* default to standard genetic code if not set */
      }
      if (val > 0 && val < NUM_GENETIC_CODES) {
        code = ValNodeNew (NULL);
        if (code != NULL) {
          code->choice = 254;
          vnp = ValNodeNew (NULL);
          code->data.ptrvalue = vnp;
          if (vnp != NULL) {
            vnp->choice = 2;
            vnp->data.intvalue = (Int4) val;
          }
        }
        crp->genetic_code = code;
        globalid = NULL;
        cfp = cpp->cfp;
        if (cfp != NULL) {
          slp = DialogToPointer (cfp->product);
          if (slp != NULL) {
            sip = SeqLocId (slp);
            if (sip != NULL) {
              bsp = BioseqFind (sip);
              if (bsp != NULL) {
                globalid = SeqIdFindBest (bsp->id, 0);
              }
            }
            SeqLocFree (slp);
          }
        }
        if (globalid != NULL) {
          cbp = DialogToPointer (cpp->cdBrk);
          ConvertToNTCodeBreak (cpp->cfp, crp, cbp);
          crp->code_break = cbp;
        }
      }
    }
  }
  return (Pointer) crp;
}

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

static void CleanupCdRgnPage (GraphiC g, VoidPtr data)

{
  CdRgnPagePtr  cpp;
  Int2          j;

  cpp = (CdRgnPagePtr) data;
  if (cpp != NULL) {
    if (cpp->alist != NULL) {
      for (j = 0; j < cpp->count; j++) {
        MemFree (cpp->alist [j].name);
      }
    }
    MemFree (cpp->alist);
  }
  MemFree (data);
}

static void SetCdRgnImportExportItems (CdRgnFormPtr cfp)

{
  CdRgnPagePtr  cpp;
  IteM          exportItm;
  IteM          importItm;

  if (cfp != NULL) {
    importItm = FindFormMenuItem ((BaseFormPtr) cfp, VIB_MSG_IMPORT);
    exportItm = FindFormMenuItem ((BaseFormPtr) cfp, VIB_MSG_EXPORT);
    switch (cfp->currentPage) {
      case CODING_REGION_PAGE :
        cpp = (CdRgnPagePtr) GetObjectExtra (cfp->data);
        if (cpp != NULL) {
          if (cpp->currentPage == 3) {
            SafeSetTitle (importItm, "Import Protein FASTA...");
            SafeSetTitle (exportItm, "Export Protein FASTA...");
            SafeEnable (importItm);
            SafeEnable (exportItm);
          } else {
            SafeSetTitle (importItm, "Import...");
            SafeSetTitle (exportItm, "Export...");
            SafeDisable (importItm);
            SafeDisable (exportItm);
          }
        } else {
          SafeSetTitle (importItm, "Import...");
          SafeSetTitle (exportItm, "Export...");
          SafeDisable (importItm);
          SafeDisable (exportItm);
        }
        break;
      case COMMON_PAGE :
        SafeSetTitle (importItm, "Import...");
        SafeSetTitle (exportItm, "Export...");
        SafeDisable (importItm);
        SafeDisable (exportItm);
        break;
      case LOCATION_PAGE :
        SafeSetTitle (importItm, "Import SeqLoc...");
        SafeSetTitle (exportItm, "Export SeqLoc...");
        SafeEnable (importItm);
        SafeEnable (exportItm);
        break;
      default :
        break;
    }
  }
}

static CharPtr cdRgnTabs [] = {
  "General", "Exceptions", "Protein", "Product", NULL
};

static void ChangeCdRgnSubPage (VoidPtr data, Int2 newval, Int2 oldval)

{
  CdRgnPagePtr  cpp;

  cpp = (CdRgnPagePtr) data;
  if (cpp != NULL) {
    cpp->currentPage = newval;
    if (oldval >= 0 && oldval <= 3) {
      SafeHide (cpp->cdRgnGrp [oldval]);
    }
    /*
    if (newval == 3) {
      SafeHide (cpp->commonGrp);
    }
    */
    if (newval >= 0 && newval <= 3) {
      SafeShow (cpp->cdRgnGrp [newval]);
    }
    /*
    if (newval >= 0 && newval < 3) {
      SafeShow (cpp->commonGrp);
    }
    */
    SetCdRgnImportExportItems (cpp->cfp);
    Update ();
  }
}

static Boolean ProtViewLaunchFunc (GatherContextPtr gcp)

{
  BioseqPtr        bsp;
  Int2             handled;
  SeqViewProcsPtr  svpp;

  if (gcp == NULL) return TRUE;
  bsp = (BioseqPtr) gcp->userdata;
  if (bsp == NULL) return TRUE;
  if (gcp->thistype == OBJ_BIOSEQ) {
    if (bsp == (BioseqPtr) gcp->thisitem) {
      WatchCursor ();
      svpp = (SeqViewProcsPtr) GetAppProperty ("SeqDisplayForm");
      if (svpp != NULL) {
        svpp->forceSeparateViewer = TRUE;
      }
      handled = GatherProcLaunch (OMPROC_VIEW, FALSE, gcp->entityID, gcp->itemID,
                                  OBJ_BIOSEQ, 0, 0, OBJ_BIOSEQ, 0);
      ArrowCursor ();
      if (handled != OM_MSG_RET_DONE || handled == OM_MSG_RET_NOPROC) {
        Message (MSG_ERROR, "Unable to launch viewer on protein product.");
      }
      return FALSE;
    }
  }
  return TRUE;
}

static void LaunchProteinViewer (ButtoN b)

{
  BioseqPtr     bsp;
  CdRgnFormPtr  cfp;
  GatherScope   gs;
  SeqIdPtr      sip;
  SeqLocPtr     slp;

  cfp = (CdRgnFormPtr) GetObjectExtra (b);
  if (cfp != NULL) {
    slp = DialogToPointer (cfp->product);
    if (slp != NULL) {
      sip = SeqLocId (slp);
      if (sip != NULL) {
        bsp = BioseqFind (sip);
        if (bsp != NULL) {
          MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
          gs.seglevels = 1;
          gs.get_feats_location = TRUE;
          MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
          gs.ignore[OBJ_BIOSEQ] = FALSE;
          gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
          GatherEntity (cfp->input_entityID, (Pointer) bsp, ProtViewLaunchFunc, &gs);
        }
      }
      SeqLocFree (slp);
    }
  }
}

static Boolean EdProtSeqLaunchFunc (GatherContextPtr gcp)

{
  BioseqPtr  bsp;
  Int2       handled;

  if (gcp == NULL) return TRUE;
  bsp = (BioseqPtr) gcp->userdata;
  if (bsp == NULL) return TRUE;
  if (gcp->thistype == OBJ_BIOSEQ) {
    if (bsp == (BioseqPtr) gcp->thisitem) {
      WatchCursor ();
      handled = GatherProcLaunch (OMPROC_EDIT, FALSE, gcp->entityID, gcp->itemID,
                                  OBJ_BIOSEQ, OBJ_BIOSEQ, Seq_repr_raw, OBJ_BIOSEQ, Seq_repr_raw);
      ArrowCursor ();
      if (handled != OM_MSG_RET_DONE || handled == OM_MSG_RET_NOPROC) {
        Message (MSG_ERROR, "Unable to launch editor on protein product.");
      }
      return FALSE;
    }
  }
  return TRUE;
}

static void LaunchProtSeqEd (ButtoN b)

{
  BioseqPtr     bsp;
  CdRgnFormPtr  cfp;
  GatherScope   gs;
  SeqIdPtr      sip;
  SeqLocPtr     slp;

  cfp = (CdRgnFormPtr) GetObjectExtra (b);
  if (cfp != NULL) {
    slp = DialogToPointer (cfp->product);
    if (slp != NULL) {
      sip = SeqLocId (slp);
      if (sip != NULL) {
        bsp = BioseqFind (sip);
        if (bsp != NULL) {
Message (MSG_OK, "Please close and reopen the coding region window after editing the sequence");
          MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
          gs.seglevels = 1;
          gs.get_feats_location = TRUE;
          MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
          gs.ignore[OBJ_BIOSEQ] = FALSE;
          gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
          GatherEntity (cfp->input_entityID, (Pointer) bsp, EdProtSeqLaunchFunc, &gs);
        }
      }
      SeqLocFree (slp);
    }
  }
}

static void LaunchProtFeatEd (ButtoN b)

{
  CdRgnFormPtr  cfp;
  Int2          handled;

  cfp = (CdRgnFormPtr) GetObjectExtra (b);
  if (cfp != NULL && cfp->protFound) {
    WatchCursor ();
    SafeHide (cfp->protTextGrp);
    SafeShow (cfp->protPromptGrp);
    Update ();
    handled = GatherProcLaunch (OMPROC_EDIT, FALSE, cfp->protEntityID, cfp->protItemID,
                                cfp->protItemtype, 0, 0, cfp->protItemtype, 0);
    ArrowCursor ();
    Update ();
    if (handled != OM_MSG_RET_DONE || handled == OM_MSG_RET_NOPROC) {
      Message (MSG_ERROR, "Unable to launch editor on protein feature.");
    }
  }
}

static void ChangeCdRgnPage (VoidPtr data, Int2 newval, Int2 oldval)

{
  CdRgnFormPtr  cfp;

  cfp = (CdRgnFormPtr) data;
  if (cfp != NULL) {
    cfp->currentPage = newval;
    SafeHide (cfp->pages [oldval]);
    SafeShow (cfp->pages [newval]);
    switch (newval) {
      case CODING_REGION_PAGE :
        SendMessageToDialog (cfp->data, VIB_MSG_ENTER);
        break;
      case COMMON_PAGE :
        break;
      case LOCATION_PAGE :
        SendMessageToDialog (cfp->location, VIB_MSG_ENTER);
        cfp->locvisited = TRUE;
        break;
      default :
        break;
    }
    SetCdRgnImportExportItems (cfp);
    Update ();
  }
}

typedef struct protgatherlist {
  ObjMgrPtr  omp;
  SeqLocPtr  slp;
  Int4       min;
  Char       protName [128];
  Char       protDesc [128];
  Uint2      entityID;
  Uint2      itemID;
  Uint2      itemtype;
  Boolean    protFound;
} ProtGatherList, PNTR ProtGatherPtr;

static Boolean ProtMatchFunc (GatherContextPtr gcp)

{
  Int4           diff;
  ObjMgrTypePtr  omtp;
  ProtGatherPtr  pgp;
  ProtRefPtr     prp;
  SeqFeatPtr     sfp;
  Char           thislabel [41];
  ValNodePtr     vnp;

  if (gcp == NULL) return TRUE;

  pgp = (ProtGatherPtr) gcp->userdata;
  if (pgp == NULL) return TRUE;

  thislabel [0] = '\0';

  if (gcp->thistype == OBJ_SEQFEAT) {
    sfp = (SeqFeatPtr) gcp->thisitem;
    if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT && sfp->data.value.ptrvalue != NULL) {
      omtp = ObjMgrTypeFind (pgp->omp, gcp->thistype, NULL, NULL);
      if (omtp == NULL) {
        return TRUE;
      }
      if (omtp->labelfunc != NULL) {
        (*(omtp->labelfunc)) (gcp->thisitem, thislabel, 40, OM_LABEL_CONTENT);
      }
      if (thislabel [0] != '\0') {
        diff = SeqLocAinB (sfp->location, pgp->slp);
        if (diff >= 0) {
          if (diff < pgp->min) {
            pgp->min = diff;
            prp = (ProtRefPtr) sfp->data.value.ptrvalue;
            pgp->protName [0] = '\0';
            pgp->protDesc [0] = '\0';
            vnp = prp->name;
            if (vnp != NULL) {
              StringNCpy_0 (pgp->protName, (CharPtr) vnp->data.ptrvalue, sizeof (pgp->protName));
            }
            StringNCpy_0 (pgp->protDesc, prp->desc, sizeof (pgp->protDesc));
            pgp->entityID = gcp->entityID;
            pgp->itemID = gcp->itemID;
            pgp->itemtype = gcp->thistype;
            pgp->protFound = TRUE;
          }
        }
      }
    }
  }

  return TRUE;
}

static void SetBestProteinFeature (CdRgnFormPtr cfp, SeqLocPtr product)

{
  GatherScope     gs;
  ProtGatherList  pgl;

  if (cfp == NULL) return;
  pgl.protName [0] = '\0';
  pgl.protDesc [0] = '\0';
  pgl.entityID = 0;
  pgl.itemID = 0;
  pgl.itemtype = 0;
  pgl.protFound = FALSE;
  if (product != NULL) {
    pgl.omp = ObjMgrGet ();
    pgl.slp = product;
    pgl.min = INT4_MAX;
    MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
    gs.seglevels = 1;
    gs.get_feats_location = TRUE;
    MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
    gs.ignore[OBJ_BIOSEQ] = FALSE;
    gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
    gs.ignore[OBJ_SEQFEAT] = FALSE;
    gs.ignore[OBJ_SEQANNOT] = FALSE;
    GatherEntity (cfp->input_entityID, (Pointer) &pgl, ProtMatchFunc, &gs);
  }
  SafeSetTitle (cfp->protNameText, pgl.protName);
  SafeSetTitle (cfp->protDescText, pgl.protDesc);
  cfp->protEntityID = pgl.entityID;
  cfp->protItemID = pgl.itemID;
  cfp->protItemtype = pgl.itemtype;
  cfp->protFound = pgl.protFound;
}

static void CdRegionPtrToForm (ForM f, Pointer data)

{
  BioseqPtr       bsp;
  CdRgnFormPtr    cfp;
  MolInfoPtr      mip;
  ProtRefPtr      prp;
  ValNodePtr      sdp;
  SeqEntryPtr     sep;
  SeqFeatPtr      sfp;
  SeqIdPtr        sip;
  ValNodePtr      vnp;
  SeqFeatXrefPtr  xref;

  cfp = (CdRgnFormPtr) GetObjectExtra (f);
  if (cfp != NULL) {
    cfp->locvisited = TRUE;  /* even passing NULL will clear this protection */
    sfp = (SeqFeatPtr) data;
    if (sfp != NULL) {
      PointerToDialog (cfp->location, sfp->location);
      PointerToDialog (cfp->product, sfp->product);
      PointerToDialog (cfp->data, sfp->data.value.ptrvalue);
      SeqFeatPtrToCommon ((FeatureFormPtr) cfp, sfp);
      SetBestProteinFeature (cfp, sfp->product);
      if (! cfp->protFound) {
        xref = sfp->xref;
        while (xref != NULL && xref->data.choice != SEQFEAT_PROT) {
          xref = xref->next;
        }
        if (xref != NULL) {
          prp = (ProtRefPtr) xref->data.value.ptrvalue;
          if (prp != NULL) {
            vnp = prp->name;
            if (vnp != NULL) {
              SafeSetTitle (cfp->protNameText, (CharPtr) vnp->data.ptrvalue);
            }
            SafeSetTitle (cfp->protDescText, prp->desc);
          }
        }
      }
      if (sfp->product != NULL) {
        SafeHide (cfp->protSeqIdGrp);
        sip = SeqLocId (sfp->product);
        if (sip != NULL) {
          bsp = BioseqFind (sip);
          if (bsp != NULL) {
            PointerToDialog (cfp->protseq, bsp);
            SafeEnable (cfp->launchBtn);
            SafeEnable (cfp->edProtBtn);
            SafeEnable (cfp->edSeqBtn);
            SafeSetStatus (cfp->conceptTransA, FALSE);
            sep = SeqMgrGetSeqEntryForData (bsp);
            if (sep == NULL) return;
            sdp = SeqEntryGetSeqDescr (sep, Seq_descr_molinfo, NULL);
            if (sdp == NULL) return;
            mip = (MolInfoPtr) sdp->data.ptrvalue;
            if (mip == NULL) return;
            SafeSetStatus (cfp->conceptTransA, (Boolean) (mip->tech == 13));
            return;
          }
        }
      }
    }
    SafeSetTitle (cfp->protSeqIdTxt, "");
    SafeShow (cfp->protSeqIdGrp);
    SafeDisable (cfp->launchBtn);
    SafeDisable (cfp->edProtBtn);
    SafeDisable (cfp->edSeqBtn);
    SafeSetStatus (cfp->conceptTransA, FALSE);
  }
}

typedef struct productpage {
  DIALOG_MESSAGE_BLOCK
  Int2          count;
  SeqEntryPtr   PNTR bsptr;
  Handle        product;
  Boolean       usePopupForProduct;
  Int2          oldval;
  Boolean       mRNAnotCDS;
  CdRgnFormPtr  cfp;
} ProductPage, PNTR ProductPagePtr;

static void FillInProducts (SeqEntryPtr sep, Pointer mydata,
                            Int4 index, Int2 indent)

{
  ProductPagePtr  ppp;

  if (sep != NULL && mydata != NULL) {
    ppp = (ProductPagePtr) mydata;
    ppp->bsptr [index] = sep;
  }
}

static void PopulateProductControl (Handle pr, ProductPagePtr ppp)

{
  BioseqPtr    bsp;
  Int2         i;
  size_t       len;
  Boolean      okay;
  SeqEntryPtr  sep;
  Char         str [128];
  /*
  Char         ch;
  CharPtr      ttl;
  */

  if (pr != NULL && ppp != NULL && ppp->bsptr != NULL) {
    if (ppp->usePopupForProduct) {
      PopupItem (pr, "     ");
    } else {
      ListItem (pr, "     ");
    }
    for (i = 0; i < ppp->count; i++) {
      sep = ppp->bsptr [i];
      if (sep != NULL && sep->choice == 1 && sep->data.ptrvalue != NULL) {
        bsp = (BioseqPtr) sep->data.ptrvalue;
        okay = FALSE;
        if (ppp->mRNAnotCDS) {
          if (ISA_na (bsp->mol)) {
            okay = TRUE;
          }
        } else {
          if (ISA_aa (bsp->mol)) {
            okay = TRUE;
          }
        }
        if (okay) {
          SeqIdWrite (bsp->id, str, PRINTID_REPORT, sizeof (str));
          /*
          StringCat (str, " ");
          len = StringLen (str);
          ttl = BioseqGetTitle (bsp);
          StringNCat (str, ttl, sizeof (str) - len - 2);
          ttl = str;
          ch = *ttl;
          while (ch != '\0') {
            if (ch == ';') {
              *ttl = ',';
            }
            ttl++;
            ch = *ttl;
          }
          */
          len = StringLen (str);
          if (len > 60) {
            StringCpy (str + 60, "...");
          }
          if (ppp->usePopupForProduct) {
            PopupItem (pr, str);
          } else {
            ListItem (pr, str);
          }
        }
      }
    }
    SetValue (pr, 1);
    ppp->oldval = 1;
  }
}

static void SeqLocPtrToProduct (DialoG d, Pointer data)

{
  BioseqPtr       bsp;
  Int2            found;
  Int2            i;
  Boolean         okay;
  ProductPagePtr  ppp;
  Int2            prots;
  SeqEntryPtr     sep;
  SeqIdPtr        sip;
  SeqLocPtr       slp;
  BioseqPtr       theBsp;

  ppp = (ProductPagePtr) GetObjectExtra (d);
  slp = (SeqLocPtr) data;
  if (ppp != NULL) {
    found = 0;
    if (slp != NULL) {
      sip = SeqLocId (slp);
      if (sip != NULL) {
        theBsp = BioseqFind (sip);
        if (theBsp != NULL) {
          prots = 0;
          for (i = 0; i < ppp->count; i++) {
            sep = ppp->bsptr [i];
            if (sep != NULL && sep->choice == 1 && sep->data.ptrvalue != NULL) {
              bsp = (BioseqPtr) sep->data.ptrvalue;
              okay = FALSE;
              if (ppp->mRNAnotCDS) {
                if (ISA_na (bsp->mol)) {
                  okay = TRUE;
                }
              } else {
                if (ISA_aa (bsp->mol)) {
                  okay = TRUE;
                }
              }
              if (okay) {
                prots++;
                if (bsp == theBsp) {
                  found = prots;
                }
              }
            }
          }
        }
      }
    }
    SetValue (ppp->product, found + 1);
    ppp->oldval = found + 1;
    if (! ppp->usePopupForProduct) {
      SetOffset (ppp->product, 0, found);
    }
  }
}

static Pointer ProductToSeqLocPtr (DialoG d)

{
  BioseqPtr       bsp;
  CdRgnFormPtr    cfp;
  Int2            i;
  Boolean         okay;
  ProductPagePtr  ppp;
  Int2            prots;
  SeqEntryPtr     sep;
  SeqLocPtr       slp;
  BioseqPtr       theBsp;
  Int2            val;

  slp = NULL;
  ppp = (ProductPagePtr) GetObjectExtra (d);
  if (ppp != NULL && ppp->cfp != NULL) {
    cfp = ppp->cfp;
    if (cfp->saveAsMiscFeat) {
      return NULL;
    }
  }
  if (ppp != NULL && ppp->bsptr != NULL) {
    if (ppp->mRNAnotCDS) {
      if (! Visible (d)) return NULL; /* not an mRNA, do not make product */
    }
    val = GetValue (ppp->product) - 1;
    prots = 0;
    theBsp = NULL;
    for (i = 0; i < ppp->count; i++) {
      sep = ppp->bsptr [i];
      if (sep != NULL && sep->choice == 1 && sep->data.ptrvalue != NULL) {
        bsp = (BioseqPtr) sep->data.ptrvalue;
        okay = FALSE;
        if (ppp->mRNAnotCDS) {
          if (ISA_na (bsp->mol)) {
            okay = TRUE;
          }
        } else {
          if (ISA_aa (bsp->mol)) {
            okay = TRUE;
          }
        }
        if (okay) {
          prots++;
          if (prots == val) {
            theBsp = bsp;
          }
        }
      }
    }
    if (ppp->cfp != NULL && ppp->cfp->usethisbioseq != NULL) {
      theBsp = ppp->cfp->usethisbioseq;
    }
    if (theBsp != NULL && theBsp->id != NULL) {
      slp = ValNodeNew (NULL);
      if (slp != NULL) {
        slp->choice = 3;
        slp->data.ptrvalue = SeqIdStripLocus (SeqIdDup (SeqIdFindBest (theBsp->id, 0)));
      }
    }
  }
  return (Pointer) slp;
}

static void CleanupProductPage (GraphiC g, VoidPtr data)

{
  ProductPagePtr  ppp;

  ppp = (ProductPagePtr) data;
  if (ppp != NULL) {
    MemFree (ppp->bsptr);
  }
  MemFree (data);
}

static CharPtr protChangeMsg =
"This will change the protein product assigned to this CDS.\n\
It should only be done if the assignment was incorrect.\n\
Are you really sure you want to do this?";

static void ChangeProteinView (Handle obj)

{
  MsgAnswer       ans;
  BioseqPtr       bsp;
  CdRgnFormPtr    cfp;
  ProductPagePtr  ppp;
  SeqIdPtr        sip;
  SeqLocPtr       slp;

  ppp = GetObjectExtra (obj);
  if (ppp != NULL && ppp->cfp != NULL) {
    if (ppp->mRNAnotCDS) return; /* using for mRNA product */
    ans = Message (MSG_OKC, "%s", protChangeMsg);
    if (ans == ANS_CANCEL) {
      SetValue (obj, ppp->oldval);
      return;
    }
    cfp = ppp->cfp;
    if (GetValue (obj) == 1) {
      SafeSetTitle (cfp->protSeqIdTxt, "");
      SafeShow (cfp->protSeqIdGrp);
      SafeDisable (cfp->launchBtn);
      SafeDisable (cfp->edProtBtn);
      SafeDisable (cfp->edSeqBtn);
      SetBestProteinFeature (cfp, NULL);
      PointerToDialog (cfp->protseq, NULL);
      Update ();
      return;
    }
    SafeHide (cfp->protSeqIdGrp);
    slp = (SeqLocPtr) ProductToSeqLocPtr (ppp->dialog);
    SetBestProteinFeature (cfp, slp);
    if (slp != NULL) {
      sip = SeqLocId (slp);
      if (sip != NULL) {
        bsp = BioseqFind (sip);
        if (bsp != NULL) {
          PointerToDialog (cfp->protseq, bsp);
          SafeEnable (cfp->launchBtn);
          SafeEnable (cfp->edProtBtn);
          SafeEnable (cfp->edSeqBtn);
        }
      }
      SeqLocFree (slp);
    }
    Update ();
  }
}


static DialoG CreateProteinOrMRNAProductDialog (GrouP h, CharPtr title, Boolean mRNAnotCDS,
                                                SeqEntryPtr sep, CdRgnFormPtr cfp)

{
  GrouP           m;
  GrouP           p;
  ProductPagePtr  ppp;
  GrouP           s;

  p = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  ppp = (ProductPagePtr) MemNew (sizeof (ProductPage));
  if (ppp != NULL) {

    SetObjectExtra (p, ppp, CleanupProductPage);
    ppp->dialog = (DialoG) p;
    ppp->todialog = SeqLocPtrToProduct;
    ppp->fromdialog = ProductToSeqLocPtr;
    ppp->testdialog = NULL;
    ppp->mRNAnotCDS = mRNAnotCDS;

    if (title != NULL && title [0] != '\0') {
      s = NormalGroup (p, 0, -2, title, systemFont, NULL);
    } else {
      s = HiddenGroup (p, 0, -2, NULL);
    }
    m = HiddenGroup (s, -2, 0, NULL);
    /*
    SetGroupSpacing (m, 10, 10);
    */

    StaticPrompt (m, "Product", 0, popupMenuHeight, programFont, 'l');
    ppp->count = (Int2) BioseqCount (sep);
    ppp->bsptr = MemNew (sizeof (BioseqPtr) * ppp->count);
    if (ppp->bsptr != NULL) {
      BioseqExplore (sep, (Pointer) ppp, FillInProducts);
    }
    ppp->cfp = cfp;
    if (ppp->count < 32) {
      ppp->usePopupForProduct = TRUE;
      ppp->product = (Handle) PopupList (m, TRUE, (PupActnProc) ChangeProteinView);
    } else {
      ppp->usePopupForProduct = FALSE;
      ppp->product = (Handle) SingleList (m, 10, 3, (LstActnProc) ChangeProteinView);
    }
    SetObjectExtra (ppp->product, ppp, NULL);
    PopulateProductControl (ppp->product, ppp);
  }

  return (DialoG) p;
}

extern SeqLocPtr PredictCodingRegion (BioseqPtr nuc, BioseqPtr prot, Int2 genCode)

{
  BioseqPtr    bsp;
  SeqLocPtr    oldslp;
  SeqAnnotPtr  sap;
  SeqFeatPtr   sfp;
  SeqIdPtr     sip;
  SeqLocPtr    slp;

  if (nuc == NULL || prot == NULL) return NULL;
  if (nuc->length + 5 >= MAXALLOC) {
    Message (MSG_OK, "Nucleotide sequence is too large for coding region prediction.");
    return NULL;
  }
  slp = NULL;
  /*
  slp = FindSingleCodingInterval (nuc, prot, genCode);
  */
  if (slp == NULL) {
    sap = SuggestCodingRegion (nuc, prot, genCode);
    if (sap != NULL && sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      if (sfp != NULL && sfp->data.choice == SEQFEAT_CDREGION) {
        slp = sfp->location;
        sfp->location = NULL;
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
      }
    }
    sap = SeqAnnotFree (sap);
  }
  StripLocusFromSeqLoc (slp);
  return slp;
}

static void DoPredictCdRegion (CdRgnFormPtr cfp)

{
  Int2          code;
  CdRgnPagePtr  cpp;
  BioseqPtr     nucbsp;
  BioseqPtr     protbsp;
  SeqEntryPtr   psep;
  SeqEntryPtr   sep;
  SeqLocPtr     slp;
  Int2          val;

  if (cfp != NULL) {
    WatchCursor ();
    Update ();
    code = 0;
    cpp = (CdRgnPagePtr) GetObjectExtra (cfp->data);
    if (cpp != NULL) {
      val = gcIndexToId [GetValue (cpp->geneticCode)];
      if (val > 0 && val < NUM_GENETIC_CODES) {
        code = val;
      }
    }
    psep = SeqEntryNew ();
    if (psep != NULL) {
      protbsp = (BioseqPtr) DialogToPointer (cfp->protseq);
      psep->choice = 1;
      psep->data.ptrvalue = protbsp;
      SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) protbsp, psep);
      sep = GetBestTopParentForItemID (cfp->input_entityID,
                                       cfp->input_itemID,
                                       cfp->input_itemtype);
      /*
      sep = GetTopSeqEntryForEntityID (cfp->input_entityID);
      */
      if (sep != NULL) {
        sep = FindNucSeqEntry (sep);
      }
      if (sep != NULL && sep->data.ptrvalue != NULL &&
          protbsp != NULL && protbsp->length > 0) {
        nucbsp = (BioseqPtr) sep->data.ptrvalue;
        slp = PredictCodingRegion (nucbsp, protbsp, code);
        PointerToDialog (cfp->location, slp);
        if (slp != NULL) {
          cfp->locvisited = TRUE;
        }
        SeqLocFree (slp);
      }
      SeqEntryFree (psep);
    }
    ArrowCursor ();
    Update ();
  }
}

static void PredictCdRegion (ButtoN b)

{
  CdRgnFormPtr  cfp;

  cfp = (CdRgnFormPtr) GetObjectExtra (b);
  DoPredictCdRegion (cfp);
}

static void DoTranslateProtein (CdRgnFormPtr cfp)

{
  ByteStorePtr  bs;
  BioseqPtr     bsp;
  Char          ch;
  CharPtr       prot;
  CharPtr       ptr;
  SeqFeatPtr    sfp;

  if (cfp != NULL) {
    sfp = SeqFeatNew ();
    if (sfp != NULL) {
      sfp->data.choice = SEQFEAT_CDREGION;
      sfp->data.value.ptrvalue = DialogToPointer (cfp->data);
      sfp->location = DialogToPointer (cfp->location);
      if (sfp->location != NULL) {
        bs = ProteinFromCdRegionEx (sfp, TRUE, FALSE);
        if (bs != NULL) {
          prot = BSMerge (bs, NULL);
          bs = BSFree (bs);
          if (prot != NULL) {
            ptr = prot;
            ch = *ptr;
            while (ch != '\0') {
              *ptr = TO_UPPER (ch);
              ptr++;
              ch = *ptr;
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
          }
          bsp = BioseqNew ();
          if (bsp != NULL) {
            bsp->repr = Seq_repr_raw;
            bsp->mol = Seq_mol_aa;
            bsp->seq_data_type = Seq_code_ncbieaa;
            bsp->seq_data = bs;
            bsp->length = BSLen (bs);
            PointerToDialog (cfp->protseq, bsp);
            bsp->seq_data = NULL;
            BioseqFree (bsp);
          }
          MemFree (prot);
        }
        BSFree (bs);
      }
    }
    SeqFeatFree (sfp);
  }
}

extern void CdRgnTranslateWithFrame (ForM f, Uint1 frame)

{
  CdRgnFormPtr  cfp;
  CdRgnPagePtr  cpp;

  cfp = (CdRgnFormPtr) GetObjectExtra (f);
  if (cfp != NULL) {
    cpp = (CdRgnPagePtr) GetObjectExtra (cfp->data);
    if (cpp != NULL) {
      if (frame < 2) {
        SafeSetValue (cpp->frame, frame + 1);
      } else {
        SafeSetValue (cpp->frame, 0);
      }
    }
    DoTranslateProtein (cfp);
  }
}

static void TranslateProtein (ButtoN b)

{
  CdRgnFormPtr  cfp;

  cfp = (CdRgnFormPtr) GetObjectExtra (b);
  if (cfp != NULL) {
    DoTranslateProtein (cfp);
  }
}

static void SetSynchProc (ButtoN b)

{
  CdRgnFormPtr  cfp;

  cfp = (CdRgnFormPtr) GetObjectExtra (b);
  if (cfp != NULL) {
    cfp->autoUpdatePartials = GetStatus (b);
  }
}

static void SetRetransProc (ButtoN b)

{
  CdRgnFormPtr  cfp;

  cfp = (CdRgnFormPtr) GetObjectExtra (b);
  if (cfp != NULL) {
    cfp->autoRetranslate = GetStatus (b);
  }
}

EnumFieldAssocPtr codebreak_alists [] = {
  NULL, NULL
};

Uint2 codebreak_types [] = {
  TAGLIST_TEXT, TAGLIST_POPUP
};

Uint2 codebreak_widths [] = {
  5, 0
};

static void CreateCodeBreakAlist (CdRgnPagePtr cpp)

{
  Char             ch;
  Int2             count;
  Uint1            first;
  Uint1            i;
  Int2             j;
  Char             item [77];
  Uint1            last;
  SeqCodeTablePtr  sctp;
  CharPtr          str;

  if (cpp != NULL) {
    sctp = SeqCodeTableFind (Seq_code_ncbieaa);
    first = FirstResidueInCode (sctp);
    last = LastResidueInCode (sctp);
    count = last - 65 + 4;
    cpp->alist = MemNew (sizeof (EnumFieldAssoc) * (size_t) count);
    cpp->count = count;
    if (cpp->alist != NULL) {
      j = 0;
      cpp->alist [j].name = StringSave ("     ");
      cpp->alist [j].value = (UIEnum) 0;
      j++;
      for (i = 65; i <= last; i++) {
        if (i == 74 || i == 79) {
          continue;
        }
        ch = GetSymbolForResidue (sctp, i);
        str = (CharPtr) GetNameForResidue (sctp, i);
        cpp->alist [j].name = StringSave (str);
        cpp->alist [j].value = (UIEnum) ch;
        sprintf (item, "%c    %s", ch, str);
        j++;
      }
      i = 42;
      ch = GetSymbolForResidue (sctp, i);
      str = (CharPtr) GetNameForResidue (sctp, i);
      cpp->alist [j].name = StringSave (str);
      cpp->alist [j].value = (UIEnum) ch;
      sprintf (item, "%c    %s", ch, str);
      j++;
      cpp->alist [j].name = NULL;
      cpp->alist [j].value = (UIEnum) 0;
    }
    cpp->alists [0] = NULL;
    cpp->alists [1] = cpp->alist;
  }
}

static void CodeBreakPtrToCodeBreakDialog (DialoG d, Pointer data)

{
  Uint1           aa;
  BioseqPtr       bsp;
  ValNodePtr      head;
  Int2            j;
  CodeBreakPtr    list;
  Int4            pos;
  SeqIdPtr        sip;
  SeqMapTablePtr  smtp;
  Char            str [64];
  TagListPtr      tlp;
  Char            tmp [32];
  ValNodePtr      vnp;

  tlp = (TagListPtr) GetObjectExtra (d);
  list = (CodeBreakPtr) data;
  if (tlp != NULL) {
    head = NULL;
    while (list != NULL) {
      if (list->loc != NULL) {
        vnp = ValNodeNew (head);
        if (head == NULL) {
          head = vnp;
        }
        if (vnp != NULL) {
          pos = 0;
          sip = SeqLocId (list->loc);
          if (sip != NULL) {
            bsp = BioseqFind (sip);
            if (bsp != NULL) {
              pos = GetOffsetInBioseq (list->loc, bsp, SEQLOC_START);
              aa = 0;
              switch (list->aa.choice) {
                case 1 :
                  aa = (Uint1) list->aa.value.intvalue;
                  break;
                case 2 :
                  smtp = SeqMapTableFind (Seq_code_ncbieaa, Seq_code_ncbi8aa);
                  if (smtp != NULL) {
                    aa = SeqMapTableConvert (smtp, (Uint1) list->aa.value.intvalue);
                  }
                  break;
                case 3 :
                  smtp = SeqMapTableFind (Seq_code_ncbieaa, Seq_code_ncbistdaa);
                  if (smtp != NULL) {
                    aa = SeqMapTableConvert (smtp, (Uint1) list->aa.value.intvalue);
                  }
                  break;
                default :
                  aa = (Uint1) list->aa.value.intvalue;
                  break;
              }
              sprintf (tmp, "%ld", (long) (pos + 1));
              StringCpy (str, tmp);
              StringCat (str, "\t");
              sprintf (tmp, "%d", (int) list->aa.value.intvalue);
              StringCat (str, tmp);
              StringCat (str, "\n");
              vnp->data.ptrvalue = StringSave (str);
            }
          }
        }
      }
      list = list->next;
    }
    SendMessageToDialog (tlp->dialog, VIB_MSG_RESET);
    tlp->vnp = head;
    SendMessageToDialog (tlp->dialog, VIB_MSG_REDRAW);
    for (j = 0, vnp = tlp->vnp; vnp != NULL; j++, vnp = vnp->next) {
    }
    tlp->max = MAX ((Int2) 0, (Int2) (j - tlp->rows + 1));
    CorrectBarMax (tlp->bar, tlp->max);
    CorrectBarPage (tlp->bar, tlp->rows - 1, tlp->rows - 1);
  }
}

static Pointer CodeBreakDialogToCodeBreakPtr (DialoG d)

{
  Uint1         aa;
  CodeBreakPtr  cbp;
  CodeBreakPtr  cbplast;
  Char          ch;
  CodeBreakPtr  head;
  Int2          j;
  Int2          len;
  Boolean       okay;
  SeqLocPtr     slp;
  SeqPntPtr     spntp;
  CharPtr       str;
  TagListPtr    tlp;
  CharPtr       tmp;
  int           val;
  ValNodePtr    vnp;

  head = NULL;
  tlp = (TagListPtr) GetObjectExtra (d);
  if (tlp != NULL && tlp->vnp != NULL) {
    cbp = NULL;
    cbplast = NULL;
    for (vnp = tlp->vnp; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      okay = FALSE;
      len = StringLen (str);
      for (j = 0; j < len; j++) {
        ch = str [j];
        if (ch != ' ' && ch != '\t' && ch != '\n') {
          okay = TRUE;
        }
      }
      if (okay) {
        tmp = ExtractTagListColumn ((CharPtr) vnp->data.ptrvalue, 1);
        if (tmp != NULL && sscanf (tmp, "%d", &val) == 1 && val != 0) {
          aa = (Uint1) val;
          MemFree (tmp);
          tmp = ExtractTagListColumn ((CharPtr) vnp->data.ptrvalue, 0);
          if (tmp != NULL && sscanf (tmp, "%d", &val) == 1 && val != 0) {
            cbp = CodeBreakNew ();
            if (cbplast == NULL) {
              head = cbp;
            } else {
              cbplast->next = cbp;
            }
            cbplast = cbp;
            if (cbp != NULL) {
              if (globalid != NULL) {
                spntp = SeqPntNew ();
                slp = ValNodeNew (NULL);
                slp->choice = SEQLOC_PNT;
                slp->data.ptrvalue = (Pointer) spntp;
                spntp->point = (Int4) (val - 1);
                spntp->id = SeqIdStripLocus (SeqIdDup (globalid));
                cbp->loc = slp;
              }
              cbp->aa.value.intvalue = aa;
              cbp->aa.choice = 1;
            }
          } else {
            MemFree (tmp);
          }
        } else {
          MemFree (tmp);
        }
      }
    }
  }
  return (Pointer) head;
}

static CharPtr cdsToMiscFeatWarning =
"You are indicating that this CDS is non-functional.\n\
Database policy is to convert this into a misc_feature,\n\
and to remove any protein product attached to the CDS.\n\
If this is the case, you should add some explanatory\n\
text as a comment before pressing the Accept button.\n\
Are you sure you want to check this box?";

static void ConvMiscFeatWarn (ButtoN b)

{
  MsgAnswer     ans;
  CdRgnFormPtr  cfp;

  cfp = (CdRgnFormPtr) GetObjectExtra (b);
  if (cfp == NULL) return;
  if (GetStatus (b)) {
    ans = Message (MSG_YN, "%s", cdsToMiscFeatWarning);
    if (ans == ANS_NO) {
      SetStatus (b, FALSE);
    }
  }
}

static DialoG CreateCdRgnDialog (GrouP h, CharPtr title, Int2 genCode,
                                 SeqFeatPtr sfp, CdRgnFormPtr cfp)

{
  ButtoN        b;
  GrouP         c;
  CdRgnPagePtr  cpp;
  GrouP         f;
  GrouP         g;
  Boolean       internalVersion;
  GrouP         j;
  GrouP         k;
  GrouP         m;
  GrouP         n;
  GrouP         p;
  PrompT        p2;
  GrouP         q;
  GrouP         s;
  DialoG        tbs;
  TagListPtr    tlp;
  GrouP         x;
  GrouP         y;
  GrouP         z;

  p = HiddenGroup (h, 1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  cpp = (CdRgnPagePtr) MemNew (sizeof (CdRgnPage));
  if (cpp != NULL) {

    SetObjectExtra (p, cpp, CleanupCdRgnPage);
    cpp->dialog = (DialoG) p;
    cpp->todialog = CdRegionPtrToCdRgnPage;
    cpp->fromdialog = CdRgnPageToCdRegionPtr;
    cpp->testdialog = NULL;

    if (title != NULL && title [0] != '\0') {
      s = NormalGroup (p, 0, -2, title, systemFont, NULL);
    } else {
      s = HiddenGroup (p, 0, -2, NULL);
    }
    m = HiddenGroup (s, -1, 0, NULL);
    SetGroupSpacing (m, 10, 10);

    cpp->cfp = cfp;
    tbs = CreateFolderTabs (m, cdRgnTabs, 0, 0, 0,
                            PROGRAM_FOLDER_TAB,
                            ChangeCdRgnSubPage, (Pointer) cpp);
    cpp->currentPage = 0;
    k = HiddenGroup (m, 0, 0, NULL);

    j = HiddenGroup (k, -1, 0, NULL);
    SetGroupSpacing (j, 10, 10);
    n = HiddenGroup (j, 0, 0, NULL);

    cpp->cdRgnGrp [0] = HiddenGroup (n, -1, 0, NULL);
    SetGroupSpacing (cpp->cdRgnGrp [0], 10, 20);
    q = HiddenGroup (cpp->cdRgnGrp [0], -1, 0, NULL);
    f = HiddenGroup (q, 0, -2, NULL);

    internalVersion = (Boolean) (GetAppProperty ("InternalNcbiSequin") != NULL);

    /*
    if (! internalVersion) {
      StaticPrompt (f, "Genetic Code", 0, popupMenuHeight, programFont, 'l');
      StaticPrompt (f, "Reading Frame", 0, popupMenuHeight, programFont, 'l');

      cpp->geneticCode = PopupList (f, TRUE, NULL);
      PopulateGeneticCodePopup (cpp->geneticCode);
      SetValue (cpp->geneticCode, gcIdToIndex [genCode]);
    }
    */

    g = HiddenGroup (f, 5, 0, NULL);

    /*
    if (! internalVersion) {
      cpp->frame = PopupList (g, TRUE, NULL);
      PopupItem (cpp->frame, " ");
      PopupItem (cpp->frame, "One");
      PopupItem (cpp->frame, "Two");
      PopupItem (cpp->frame, "Three");
      SetValue (cpp->frame, 1);
    }
    */

    p2 = StaticPrompt (g, "Flags", 0, popupMenuHeight, programFont, 'l');
    if (! internalVersion) {
      Hide (p2);
    }
    cpp->conflict = CheckBox (g, "Conflict", NULL);
    if (! internalVersion) {
      Hide (cpp->conflict);
    }
    cpp->orf = CheckBox (g, "ORF", NULL);
    Hide (cpp->orf);

    AlignObjects (ALIGN_MIDDLE, (HANDLE) cpp->orf,
                  (HANDLE) cpp->conflict, (HANDLE) cpp->frame, NULL);

    cfp->conceptTransA = NULL;
    if (internalVersion) {
      cfp->conceptTransA = CheckBox (q, "Author-supplied conceptual translation", NULL);
    }

    StaticPrompt (q, "", 10, stdLineHeight, programFont, 'l');
    cfp->convertToMiscFeat = CheckBox (q, "Non-functional CDS, convert to misc_feature", ConvMiscFeatWarn);
    SetObjectExtra (cfp->convertToMiscFeat, cfp, NULL);
    cfp->saveAsMiscFeat = FALSE;

    cpp->cdRgnGrp [1] = HiddenGroup (n, -1, 0, NULL);
    CreateCodeBreakAlist (cpp);
    codebreak_widths [0] = 5;
    f = HiddenGroup (cpp->cdRgnGrp [1], 4, 0, NULL);
    StaticPrompt (f, "Position", 5 * stdCharWidth, 0, programFont, 'c');
    p2 = StaticPrompt (f, "Amino Acid", 0, 0, programFont, 'c');
    cpp->cdBrk = CreateTagListDialog (cpp->cdRgnGrp [1],
                                      3, 2, 2, codebreak_types,
                                      codebreak_widths, cpp->alists,
                                      CodeBreakPtrToCodeBreakDialog,
                                      CodeBreakDialogToCodeBreakPtr);
    tlp = (TagListPtr) GetObjectExtra (cpp->cdBrk);
    if (tlp != NULL) {
      AlignObjects (ALIGN_JUSTIFY, (HANDLE) tlp->control [1], (HANDLE) p2, NULL);
    }
    Hide (cpp->cdRgnGrp [1]);

    cpp->cdRgnGrp [2] = HiddenGroup (n, -1, 0, NULL);
    f = HiddenGroup (cpp->cdRgnGrp [2], 0, 0, NULL);
    cfp->protTextGrp = HiddenGroup (f, -1, 0, NULL);
    x = HiddenGroup (cfp->protTextGrp, 2, 0, NULL);
    SetGroupSpacing (x, 3, 5);
    StaticPrompt (x, "Name", 0, dialogTextHeight, programFont, 'l');
    cfp->protNameText = DialogText (x, "", 15, NULL);
    StaticPrompt (x, "Description", 0, dialogTextHeight, programFont, 'l');
    cfp->protDescText = DialogText (x, "", 15, NULL);
    cfp->protPromptGrp = HiddenGroup (f, -1, 0, NULL);
    StaticPrompt (cfp->protPromptGrp,
                  "Press Edit Protein Feature to change protein name",
                  0, 0, programFont, 'l');
    Hide (cfp->protPromptGrp);
    z = HiddenGroup (cpp->cdRgnGrp [2], 3, 0, NULL);
    cfp->edProtBtn = PushButton (z, "Edit Protein Feature", LaunchProtFeatEd);
    SetObjectExtra (cfp->edProtBtn, cfp, NULL);
    Disable (cfp->edProtBtn);
    cfp->launchBtn = PushButton (z, "Launch Product Viewer", LaunchProteinViewer);
    SetObjectExtra (cfp->launchBtn, cfp, NULL);
    Disable (cfp->launchBtn);
    AlignObjects (ALIGN_CENTER, (HANDLE) cfp->protTextGrp, (HANDLE) cfp->protPromptGrp,
                  (HANDLE) z, NULL);
    Hide (cpp->cdRgnGrp [2]);

    /*
    cpp->commonGrp = HiddenGroup (j, -1, 0, NULL);
    SetGroupSpacing (cpp->commonGrp, 3, 5);
    y = HiddenGroup (cpp->commonGrp, 3, 0, NULL);
    x = HiddenGroup (cpp->commonGrp, 2, 0, NULL);
    b = CheckBox (x, "Retranslate on Accept", SetRetransProc);
    SetObjectExtra (b, cfp, NULL);
    cfp->autoRetranslate = FALSE;
    b = CheckBox (x, "Synchronize Partials", SetSynchProc);
    SetObjectExtra (b, cfp, NULL);
    SetStatus (b, TRUE);
    cfp->autoUpdatePartials = TRUE;
    AlignObjects (ALIGN_CENTER, (HANDLE) x, (HANDLE) y, NULL);
    */

    cpp->cdRgnGrp [3] = HiddenGroup (k, -1, 0, NULL);
    g = HiddenGroup (cpp->cdRgnGrp [3], -1, 0, NULL);

    y = NULL;
    /* if (internalVersion) { */
      y = HiddenGroup (g, 0, -2, NULL);
      StaticPrompt (y, "Genetic Code", 0, popupMenuHeight, programFont, 'l');
      StaticPrompt (y, "Reading Frame", 0, popupMenuHeight, programFont, 'l');
      cpp->geneticCode = PopupList (y, TRUE, NULL);
      PopulateGeneticCodePopup (cpp->geneticCode);
      SetValue (cpp->geneticCode, gcIdToIndex [genCode]);
      cpp->frame = PopupList (y, TRUE, NULL);
      PopupItem (cpp->frame, " ");
      PopupItem (cpp->frame, "One");
      PopupItem (cpp->frame, "Two");
      PopupItem (cpp->frame, "Three");
      SetValue (cpp->frame, 1);
    /* } */

    x = HiddenGroup (g, -1, 0, NULL);
    f = HiddenGroup (x, 4, 0, NULL);
    cfp->product = CreateProteinOrMRNAProductDialog (f, NULL, FALSE, cfp->sep, cfp);
    cfp->protSeqIdGrp = HiddenGroup (f, 2, 0, NULL);
    StaticPrompt (cfp->protSeqIdGrp, "SeqID", 0, dialogTextHeight, programFont, 'l');
    cfp->protSeqIdTxt = DialogText (cfp->protSeqIdGrp, "", 6, NULL);

    cfp->protseq = CreateBioseqRawDialog (x, NULL, Seq_mol_aa, cfp);
    cfp->usethisbioseq = NULL;

    c = HiddenGroup (g, 4, 0, NULL);
    b = PushButton (c, "Predict Interval", PredictCdRegion);
    SetObjectExtra (b, cfp, NULL);
    b = PushButton (c, "Translate Product", TranslateProtein);
    SetObjectExtra (b, cfp, NULL);
    cfp->edSeqBtn = PushButton (c, "Edit Protein Sequence", LaunchProtSeqEd);
    SetObjectExtra (cfp->edSeqBtn, cfp, NULL);
    Disable (cfp->edSeqBtn);
    AlignObjects (ALIGN_LEFT, (HANDLE) f, (HANDLE) cfp->protseq, NULL);
    AlignObjects (ALIGN_CENTER, (HANDLE) x, (HANDLE) c, (HANDLE) y, NULL);
    Hide (cpp->cdRgnGrp [3]);

    AlignObjects (ALIGN_CENTER, (HANDLE) cpp->cdRgnGrp [0],
                  (HANDLE) cpp->cdRgnGrp [1], (HANDLE) cpp->cdRgnGrp [2],
                  (HANDLE) cpp->cdRgnGrp [3], NULL);
    AlignObjects (ALIGN_CENTER, (HANDLE) tbs, (HANDLE) k, NULL);
  }

  return (DialoG) p;
}

typedef struct genefinddata {
  ObjMgrPtr     omp;
  CdRgnFormPtr  cfp;
  CharPtr       title;
  Int2          idx;
  Int2          val;
  Boolean       geneFound;
} GeneFindData, PNTR GeneFindPtr;

static Boolean GeneFindFunc (GatherContextPtr gcp)

{
  GeneFindPtr    gfp;
  ObjMgrTypePtr  omtp;
  SeqFeatPtr     sfp;
  Char           thislabel [41];

  if (gcp == NULL) return TRUE;

  gfp = (GeneFindPtr) gcp->userdata;
  if (gfp == NULL ) return TRUE;

  thislabel [0] = '\0';

  if (gcp->thistype == OBJ_SEQFEAT) {
    sfp = (SeqFeatPtr) gcp->thisitem;
    if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE) {
      omtp = ObjMgrTypeFind (gfp->omp, gcp->thistype, NULL, NULL);
      if (omtp == NULL) {
        return TRUE;
      }
      if (omtp->labelfunc != NULL) {
        (*(omtp->labelfunc)) (gcp->thisitem, thislabel, 40, OM_LABEL_CONTENT);
      }
      if (thislabel [0] != '\0') {
        gfp->idx++;
        if (StringICmp (thislabel, gfp->title) == 0) {
          gfp->val = gfp->idx;
          gfp->geneFound = TRUE;
          return FALSE;
        }
      }
    }
  }
  return TRUE;
}

static void SetGeneAndProtFields (CdRgnFormPtr cfp, CharPtr title)

{
  GeneFindData  gfd;
  GatherScope   gs;
  CharPtr       ptr;
  Char          str [128];

  if (cfp != NULL && title != NULL) {
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
        if (! StringHasNoText (str)) {
          SetTitle (cfp->protNameText, str);
        }
        if (! StringHasNoText (ptr)) {
          SetTitle (cfp->protDescText, ptr);
        }
      }
    }
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
        if (! StringHasNoText (str)) {
          SetTitle (cfp->geneSymbol, str);
        }
        if (! StringHasNoText (ptr)) {
          SetTitle (cfp->geneDesc, ptr);
        }
      }
      gfd.omp = ObjMgrGet ();
      gfd.cfp = cfp;
      gfd.title = str;
      gfd.idx = 2;
      gfd.val = 1;
      gfd.geneFound = FALSE;
      MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
      gs.seglevels = 1;
      gs.get_feats_location = TRUE;
	  MemSet ((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
	  gs.ignore[OBJ_BIOSEQ] = FALSE;
	  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
	  gs.ignore[OBJ_SEQFEAT] = FALSE;
	  gs.ignore[OBJ_SEQANNOT] = FALSE;
      gs.scope = GetBestTopParentForItemID (cfp->input_entityID,
                                            cfp->input_itemID,
                                            cfp->input_itemtype);
      GatherEntity (cfp->input_entityID, (Pointer) &gfd, GeneFindFunc, &gs);
      if (gfd.geneFound) {
        SetValue (cfp->gene, gfd.val);
        SafeHide (cfp->newGeneGrp);
        SafeShow (cfp->editGeneBtn);
      } else if (TextHasNoText (cfp->geneSymbol)) {
        SetValue (cfp->gene, 1);
        SafeHide (cfp->newGeneGrp);
        SafeShow (cfp->editGeneBtn);
      } else {
        SetValue (cfp->gene, 2);
        SafeHide (cfp->editGeneBtn);
        SafeShow (cfp->newGeneGrp);
      }
    }
  }
}

static Boolean ImportCdRgnForm (ForM f, CharPtr filename)

{
  BioseqPtr     bsp;
  CdRgnFormPtr  cfp;
  CharPtr       extension;
  FILE          *fp;
  Char          path [PATH_MAX];
  SeqEntryPtr   sep;
  SeqLocPtr     slp;
  Char          str [128];
  ValNodePtr    vnp;
  CharPtr       vnpstr;

  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  cfp = (CdRgnFormPtr) GetObjectExtra (f);
  if (cfp != NULL) {
    switch (cfp->currentPage) {
      case CODING_REGION_PAGE :
        extension = GetAppProperty ("FastaProtExtension");
        if (path [0] != '\0' || GetInputFileName (path, sizeof (path), extension, "TEXT")) {
          fp = FileOpen (path, "r");
          if (fp != NULL) {
            sep = FastaToSeqEntry (fp, FALSE);
            if (sep != NULL && sep->choice == 1 && sep->data.ptrvalue != NULL) {
              bsp = (BioseqPtr) sep->data.ptrvalue;
              PointerToDialog (cfp->protseq, bsp);
              SeqIdWrite (bsp->id, str, PRINTID_REPORT, sizeof (str));
              SafeSetTitle (cfp->protSeqIdTxt, str);
              if (bsp->descr != NULL) {
                vnp = ValNodeFindNext (bsp->descr, NULL, Seq_descr_title);
                if (vnp != NULL && vnp->data.ptrvalue != NULL) {
                  vnpstr = (CharPtr) vnp->data.ptrvalue;
                  SetGeneAndProtFields (cfp, vnpstr);
                }
              }
              slp = DialogToPointer (cfp->location);
              if (slp == NULL || (! cfp->locvisited)) {
                DoPredictCdRegion (cfp);
              }
              SeqLocFree (slp);
            }
            SeqEntryFree (sep);
          }
          FileClose (fp);
        }
        break;
      case COMMON_PAGE :
        break;
      case LOCATION_PAGE :
        return ImportDialog (cfp->location, filename);
      default :
        break;
    }
  }
  return FALSE;
}

static Boolean ExportCdRgnForm (ForM f, CharPtr filename)

{
  BioseqPtr     bsp;
  CdRgnFormPtr  cfp;
  FILE          *fp;
  ErrSev        oldErrSev;
  Char          path [PATH_MAX];
  SeqEntryPtr   sep;

  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  cfp = (CdRgnFormPtr) GetObjectExtra (f);
  if (cfp != NULL) {
    switch (cfp->currentPage) {
      case CODING_REGION_PAGE :
        if (path [0] != '\0' || GetOutputFileName (path, sizeof (path), NULL)) {
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
            sep = SeqEntryNew ();
            if (sep != NULL) {
              bsp = DialogToPointer (cfp->protseq);
              sep->choice = 1;
              sep->data.ptrvalue = (Pointer) bsp;
              oldErrSev = ErrSetMessageLevel (SEV_MAX);
              SeqEntryToFasta (sep, fp, FALSE);
              ErrSetMessageLevel (oldErrSev);
            }
            FileClose (fp);
            SeqEntryFree (sep);
            return TRUE;
          }
        }
        break;
      case COMMON_PAGE :
        break;
      case LOCATION_PAGE :
        return ExportDialog (cfp->location, filename);
      default :
        break;
    }
  }
  return FALSE;
}

static CharPtr  cdRgnFormTabs [] = {
  "Coding Region", "Properties", "Location", NULL
};

static void CdRgnFormMessage (ForM f, Int2 mssg)

{
  CdRgnFormPtr  cfp;

  cfp = (CdRgnFormPtr) GetObjectExtra (f);
  if (cfp != NULL) {
    switch (mssg) {
      case VIB_MSG_INIT :
        StdInitFeatFormProc (f);
        break;
      case VIB_MSG_IMPORT :
        ImportCdRgnForm (f, NULL);
        break;
      case VIB_MSG_EXPORT :
        ExportCdRgnForm (f, NULL);
        break;
      case VIB_MSG_PRINT :
        break;
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
        if (cfp->currentPage == LOCATION_PAGE) {
          PointerToDialog (cfp->location, NULL);
        } else {
          StdDeleteTextProc (NULL);
        }
        break;
      default :
        if (cfp->appmessage != NULL) {
          cfp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static CharPtr noLocMessage =
"Feature must have a location!  You can import a\n\
FASTA protein file from the Product subpage, which\n\
will generate a location from the protein sequence.";

static CharPtr cdsToMiscFeatMsg =
"You have indicated that this CDS is non-functional.\n\
Database policy is to convert this into a misc_feature,\n\
and to remove any protein product attached to the CDS.\n\
Are you sure you want to trigger this conversion now?";

static void CdRgnFormAcceptButtonProc (ButtoN b)

{
  MsgAnswer      ans;
  BioseqPtr      bsp;
  CdRgnFormPtr   cfp;
  GBQualPtr      gbq;
  SeqEntryPtr    old;
  OMProcControl  ompc;
  Boolean        pseudo;
  GBQualPtr      qual;
  SeqEntryPtr    sep;
  SeqIdPtr       sip;
  SeqLocPtr      slp;
  Char           str [128];
  WindoW         w;

  if (b != NULL) {
    w = ParentWindow (b);
    cfp = (CdRgnFormPtr) GetObjectExtra (b);
    if (cfp != NULL && cfp->form != NULL && cfp->actproc != NULL) {
      if (! cfp->locvisited) {
        ErrPostEx (SEV_ERROR, 0, 0, "%s", noLocMessage);
        return;
      }
      slp = DialogToPointer (cfp->location);
      if (slp == NULL) {
        ErrPostEx (SEV_ERROR, 0, 0, "%s", noLocMessage);
        return;
      }

      cfp->saveAsMiscFeat = FALSE;
      if (GetStatus (cfp->convertToMiscFeat)) {
        ans = Message (MSG_OKC, "%s", cdsToMiscFeatMsg);
        if (ans == ANS_CANCEL) return;
        SeqLocFree (slp);
        Hide (w);
        Update ();
        bsp = NULL;
        slp = DialogToPointer (cfp->product); /* need to get current product reference */
        if (slp != NULL) {
          sip = SeqLocId (slp);
          if (sip != NULL) {
            bsp = BioseqFind (sip);
          }
        }
        SeqLocFree (slp);
        cfp->this_itemtype = OBJ_SEQFEAT;
        cfp->this_subtype = FEATDEF_IMP;
        cfp->saveAsMiscFeat = TRUE; /* now cfp->product dialog will return NULL */
        FeatFormReplaceWithoutUpdateProc (cfp->form);
        GetRidOfEmptyFeatsDescStrings (cfp->input_entityID, NULL);
        if (GetAppProperty ("InternalNcbiSequin") != NULL) {
          ExtendGeneFeatIfOnMRNA (cfp->input_entityID, NULL);
        }
        if (bsp != NULL) {
          MemSet ((Pointer) (&ompc), 0, sizeof (OMProcControl));
          ompc.do_not_reload_from_cache = TRUE;
          ompc.input_entityID = cfp->input_entityID;
          ompc.input_itemID = GetItemIDGivenPointer (cfp->input_entityID, OBJ_BIOSEQ, (Pointer) bsp);
          ompc.input_itemtype = OBJ_BIOSEQ;
          if (! DetachDataForProc (&ompc, FALSE)) {
            Message (MSG_ERROR, "DetachDataForProc failed");
          }
          sep = GetBestTopParentForItemID (cfp->input_entityID,
                                           cfp->input_itemID,
                                           cfp->input_itemtype);
          RenormalizeNucProtSets (sep, TRUE);
        }
        ObjMgrSendMsg (OM_MSG_UPDATE, cfp->input_entityID,
                       cfp->input_itemID, cfp->input_itemtype);
        ObjMgrDeSelect (0, 0, 0, 0, NULL);
        Update ();
        Remove (w);
        return;
      }

      ans = ANS_OK;
      if (cfp->autoRetranslate) {
        DoTranslateProtein (cfp);
      }
      pseudo = FALSE;
      qual = (GBQualPtr) DialogToPointer (cfp->gbquals);
      gbq = qual;
      while (gbq != NULL) {
        if (StringICmp (gbq->qual, "pseudo") == 0) {
          pseudo = TRUE;
        }
        gbq = gbq->next;
      }
      GBQualFree (qual);
      bsp = (BioseqPtr) DialogToPointer (cfp->protseq);
      if (bsp == NULL) {
        if (! pseudo) {
          ans = Message (MSG_OKC, "The coding region has no peptide sequence.  Are you sure?");
        }
      } else if (pseudo) {
        /* ask about deleting the product bioseq? */
      }
      BioseqFree (bsp);
      if (ans == ANS_OK) {
        Hide (w);
        Update ();
        if (Visible (cfp->protSeqIdGrp)) {
          GetTitle (cfp->protSeqIdTxt, str, sizeof (str) - 1);
          if (StringHasNoText (str)) {
            sep = GetBestTopParentForItemID (cfp->input_entityID,
                                             cfp->input_itemID,
                                             cfp->input_itemtype);
            old = SeqEntrySetScope (sep);
            sip = MakeNewProteinSeqId (slp, NULL);
            SeqEntrySetScope (old);
            if (sip != NULL) {
              SeqIdWrite (sip, str, PRINTID_REPORT, sizeof (str));
              SafeSetTitle (cfp->protSeqIdTxt, str);
            }
            SeqIdFree (sip);
          }
        }
        (cfp->actproc) (cfp->form);
        SeqLocFree (slp);
      } else {
        SeqLocFree (slp);
        return;
      }
    }
    Update ();
    Remove (w);
  }
}

static void CdRgnFormActivate (WindoW w)

{
  CdRgnFormPtr  cfp;

  cfp = (CdRgnFormPtr) GetObjectExtra (w);
  if (cfp != NULL) {
    if (cfp->activate != NULL) {
      cfp->activate (w);
    }
    SetCdRgnImportExportItems (cfp);
  }
}

extern ForM CreateCdRgnForm (Int2 left, Int2 top, CharPtr title,
                             SeqFeatPtr sfp, SeqEntryPtr sep,
                             FormActnFunc actproc)

{
  ButtoN             b;
  GrouP              c;
  CdRgnFormPtr       cfp;
  SeqFeatPtr         copysfp;
  CdRgnPagePtr       cpp;
  GrouP              g;
  GBQualPtr          gbq;
  Int2               genCode;
  GrouP              h;
  GBQualPtr          qual;
  GrouP              s;
  StdEditorProcsPtr  sepp;
  WindoW             w;
  GrouP              x;
  GrouP              y;

  w = NULL;
  cfp = (CdRgnFormPtr) MemNew (sizeof (CdRgnForm));
  if (cfp != NULL) {
    w = FixedWindow (left, top, -10, -10, title, StdCloseWindowProc);
    SetObjectExtra (w, cfp, StdFeatFormCleanupProc);
    cfp->form = (ForM) w;
    cfp->actproc = actproc;
    cfp->toform = CdRegionPtrToForm;
    cfp->fromform = NULL;
    cfp->formmessage = CdRgnFormMessage;
    cfp->testform = NULL;
    cfp->importform = ImportCdRgnForm;
    cfp->exportform = ExportCdRgnForm;

#ifndef WIN_MAC
    CreateStdEditorFormMenus (w);
#endif

    cfp->activate = NULL;
    sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
    if (sepp != NULL) {
      cfp->activate = sepp->activateForm;
      cfp->appmessage = sepp->handleMessages;
    }
    SetActivate (w, CdRgnFormActivate);

    g = HiddenGroup (w, -1, 0, NULL);
    SetGroupSpacing (g, 3, 10);

    cfp->sep = sep;
    cfp->foldertabs = CreateFolderTabs (g, cdRgnFormTabs, CODING_REGION_PAGE,
                                        0, 0, SYSTEM_FOLDER_TAB,
                                        ChangeCdRgnPage, (Pointer) cfp);
    cfp->currentPage = CODING_REGION_PAGE;

    h = HiddenGroup (g, 0, 0, NULL);

    s = HiddenGroup (h, -1, 0, NULL);
    genCode = SeqEntryToGeneticCode (sep, NULL, NULL, 0);
    cfp->data = CreateCdRgnDialog (s, NULL, genCode, sfp, cfp);
    cfp->pages [CODING_REGION_PAGE] = s;
    Hide (cfp->pages [CODING_REGION_PAGE]);

    s = HiddenGroup (h, -1, 0, NULL);
    copysfp = NULL;
    if (sfp != NULL) {
      copysfp = AsnIoMemCopy ((Pointer) sfp,
                              (AsnReadFunc) SeqFeatAsnRead,
                              (AsnWriteFunc) SeqFeatAsnWrite);
    } else {
      copysfp = SeqFeatNew ();
      if (copysfp != NULL) {
        copysfp->data.choice = SEQFEAT_CDREGION;
      }
    }
    if (copysfp != NULL) {
      if (GetAppProperty ("InternalNcbiSequin") != NULL) {
        qual = GBQualNew ();
        if (qual != NULL) {
          qual->qual = StringSave ("pseudo");
          qual->val = StringSave ("");
          gbq = copysfp->qual;
          if (gbq != NULL) {
            while (gbq->next != NULL) {
              gbq = gbq->next;
            }
            gbq->next = qual;
          } else {
            copysfp->qual = qual;
          }
        }
        qual = GBQualNew ();
        if (qual != NULL) {
          qual->qual = StringSave ("exception");
          qual->val = StringSave ("");
          gbq = copysfp->qual;
          if (gbq != NULL) {
            while (gbq->next != NULL) {
              gbq = gbq->next;
            }
            gbq->next = qual;
          } else {
            copysfp->qual = qual;
          }
        }
      }
    }
    CreateCommonFeatureGroup (s, (FeatureFormPtr) cfp, copysfp, TRUE, TRUE);
    SeqFeatFree (copysfp);
    cfp->pages [COMMON_PAGE] = s;
    Hide (cfp->pages [COMMON_PAGE]);

    s = HiddenGroup (h, -1, 0, NULL);
    cfp->location = CreateIntervalEditorDialogEx (s, NULL, 4, 2, sep, TRUE, FALSE,
                                                  TRUE, TRUE, FALSE,
                                                  (FeatureFormPtr) cfp,
                                                  StdFeatIntEdPartialCallback);
    cfp->pages [LOCATION_PAGE] = s;
    Hide (cfp->pages [LOCATION_PAGE]);
    cfp->locvisited = FALSE;

    AlignObjects (ALIGN_CENTER, (HANDLE) cfp->pages [CODING_REGION_PAGE],
                  (HANDLE) cfp->pages [COMMON_PAGE],
                  (HANDLE) cfp->pages [LOCATION_PAGE], NULL);
    AlignObjects (ALIGN_CENTER, (HANDLE) cfp->foldertabs, (HANDLE) h, NULL);

    cpp = (CdRgnPagePtr) GetObjectExtra (cfp->data);
    cpp->commonGrp = HiddenGroup (w, -1, 0, NULL);
    SetGroupSpacing (cpp->commonGrp, 3, 5);
    y = HiddenGroup (cpp->commonGrp, 3, 0, NULL);
    x = HiddenGroup (cpp->commonGrp, 2, 0, NULL);
    b = CheckBox (x, "Retranslate on Accept", SetRetransProc);
    SetObjectExtra (b, cfp, NULL);
    cfp->autoRetranslate = FALSE;
    b = CheckBox (x, "Synchronize Partials", SetSynchProc);
    SetObjectExtra (b, cfp, NULL);
    SetStatus (b, TRUE);
    cfp->autoUpdatePartials = TRUE;
    AlignObjects (ALIGN_CENTER, (HANDLE) x, (HANDLE) y, NULL);

    c = HiddenGroup (w, 2, 0, NULL);
    b = PushButton (c, "Accept", CdRgnFormAcceptButtonProc);
    SetObjectExtra (b, cfp, NULL);
    PushButton (c, "Cancel", StdCancelButtonProc);
    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, (HANDLE) cpp->commonGrp, NULL);
    RealizeWindow (w);

    SendMessageToDialog (cfp->data, VIB_MSG_INIT);
    SendMessageToDialog (cfp->location, VIB_MSG_INIT);
    SendMessageToDialog (cfp->product, VIB_MSG_INIT);
    Show (cfp->pages [cfp->currentPage]);
    SendMessageToDialog (cfp->data, VIB_MSG_ENTER);
    Update ();
  }
  return (ForM) w;
}

typedef struct partialtrio {
  Boolean       partial;
  Boolean       partial5;
  Boolean       partial3;
  SeqEntryPtr   sep;
} PartialTrio, PNTR PartialTrioPtr;

static Boolean SetProteinPartial (GatherContextPtr gcp)

{
  SeqFeatPtr      sfp;
  PartialTrioPtr  trioPtr;

  trioPtr = (PartialTrioPtr) gcp->userdata;
  sfp = (SeqFeatPtr) gcp->thisitem;
  if (sfp != NULL && trioPtr != NULL) {
    sfp->partial = trioPtr->partial;
    if (trioPtr->sep != NULL) {
      sfp->location = SeqLocFree (sfp->location);
      sfp->location =  CreateWholeInterval (trioPtr->sep);
    }
    SetSeqLocPartial (sfp->location, trioPtr->partial5, trioPtr->partial3);
  }
  return TRUE;
}

static Boolean UpdateProteinName (GatherContextPtr gcp)

{
  CdRgnFormPtr  cfp;
  ProtRefPtr    prp;
  SeqFeatPtr    sfp;
  ValNodePtr    vnp;

  cfp = (CdRgnFormPtr) gcp->userdata;
  sfp = (SeqFeatPtr) gcp->thisitem;
  if (sfp == NULL || cfp == NULL) return TRUE;
  if (sfp->data.choice != SEQFEAT_PROT) return TRUE;
  prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  if (prp == NULL) return TRUE;
  if (prp->name != NULL) {
    if (! TextHasNoText (cfp->protNameText)) {
      vnp = prp->name;
      vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
      vnp->data.ptrvalue = SaveStringFromText (cfp->protNameText);
    } else {
      vnp = prp->name;
      vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
      prp->name = vnp->next;
      MemFree (vnp);
    }
  } else if (! TextHasNoText (cfp->protNameText)) {
    vnp = ValNodeNew (NULL);
    prp->name = vnp;
    if (vnp != NULL) {
      vnp->data.ptrvalue = SaveStringFromText (cfp->protNameText);
    }
  }
  prp->desc = MemFree (prp->desc);
  prp->desc = SaveStringFromText (cfp->protDescText);
  return TRUE;
}

extern void CdRgnFeatFormActnProc (ForM f)

{
  BioseqPtr     bsp;
  CdRgnFormPtr  cfp;
  Char          desc [128];
  ValNodePtr    descr;
  GBQualPtr     gbq;
  MolInfoPtr    mip;
  Char          name [128];
  Boolean       partial;
  Boolean       partial5;
  Boolean       partial3;
  SeqEntryPtr   prot;
  ProtRefPtr    prp;
  GBQualPtr     quals;
  SeqEntryPtr   sep;
  Uint1         seq_data_type;
  SeqFeatPtr    sfp;
  SeqIdPtr      sip;
  SeqLocPtr     slp;
  BioseqPtr     target;
  PartialTrio   trio;
  ValNodePtr    vnp;

  cfp = (CdRgnFormPtr) GetObjectExtra (f);
  sep = NULL;
  if (cfp != NULL) {
    if (! cfp->locvisited) {
      ErrPostEx (SEV_ERROR, 0, 0, "Feature must have a location!");
      return;
    }
    slp = DialogToPointer (cfp->location);
    if (slp == NULL) {
      ErrPostEx (SEV_ERROR, 0, 0, "Feature must have a location!");
      return;
    }
	/** this was causing a coredump ***
    SeqLocFree (slp);
	**********************************/
    /*
    if (cfp->autoRetranslate) {
      DoTranslateProtein (cfp);
    }
    */
    CheckSeqLocForPartial (slp, &partial5, &partial3);
    SeqLocFree (slp);
    partial = GetStatus (cfp->partial);
    slp = DialogToPointer (cfp->product);
    if (slp != NULL) {
      quals = DialogToPointer (cfp->gbquals);
      cfp->usethisbioseq = NULL;
      gbq = quals;
      while (gbq != NULL) {
        if (StringICmp (gbq->qual, "pseudo") == 0) {
          GBQualFree (quals);
          SeqLocFree (slp);
          FeatFormReplaceWithoutUpdateProc (f);
          GetRidOfEmptyFeatsDescStrings (cfp->input_entityID, NULL);
          if (GetAppProperty ("InternalNcbiSequin") != NULL) {
            ExtendGeneFeatIfOnMRNA (cfp->input_entityID, NULL);
          }
          ObjMgrSendMsg (OM_MSG_UPDATE, cfp->input_entityID,
                         cfp->input_itemID, cfp->input_itemtype);
          return;
        }
        gbq = gbq->next;
      }
      GBQualFree (quals);
      if (cfp->protFound) {
        sip = SeqLocId (slp);
        if (sip != NULL) {
          target = BioseqFind (sip);
          if (target != NULL) {
            trio.sep = NULL;
            bsp = (BioseqPtr) DialogToPointer (cfp->protseq);
            if (bsp != NULL) {
              seq_data_type = target->seq_data_type;
              target->length = bsp->length;
              target->seq_data = BSFree (target->seq_data);
              target->seq_data = bsp->seq_data;
              target->seq_data_type = bsp->seq_data_type;
              bsp->seq_data = NULL;
              BioseqConvert (target, seq_data_type);
              /*BioseqConvert (target, Seq_code_ncbistdaa);*/
            }
            BioseqFree (bsp);
            trio.sep = SeqMgrGetSeqEntryForData (target);
            if (cfp->autoUpdatePartials) {
              partial = (Boolean) (partial || partial5 || partial3);
              trio.partial = partial;
              trio.partial5 = partial5;
              trio.partial3 = partial3;
              GatherItem (cfp->protEntityID, cfp->protItemID,
                          cfp->protItemtype, (Pointer) &trio, SetProteinPartial);
            }
            if (Visible (cfp->protTextGrp)) {
              GatherItem (cfp->protEntityID, cfp->protItemID,
                          cfp->protItemtype, (Pointer) cfp, UpdateProteinName);
            }
          }
        }
      }
      FeatFormReplaceWithoutUpdateProc (f);
      sip = SeqLocId (slp);
      if (sip != NULL) {
        target = BioseqFind (sip);
        if (target != NULL) {
          bsp = (BioseqPtr) DialogToPointer (cfp->protseq);
          if (bsp != NULL) {
            seq_data_type = target->seq_data_type;
            target->length = bsp->length;
            target->seq_data = BSFree (target->seq_data);
            target->seq_data = bsp->seq_data;
            target->seq_data_type = bsp->seq_data_type;
            bsp->seq_data = NULL;
            BioseqConvert (target, seq_data_type);
            /*BioseqConvert (target, Seq_code_ncbistdaa);*/
          }
          BioseqFree (bsp);
          if ((! cfp->protFound) && Visible (cfp->protTextGrp)) {
            GetTitle (cfp->protNameText, name, sizeof (name));
            GetTitle (cfp->protDescText, desc, sizeof (desc));
            prp = CreateNewProtRef (name, desc, NULL, NULL);
            if (prp != NULL) {
              sep = SeqMgrGetSeqEntryForData (target);
              if (sep != NULL) {
                sfp = CreateNewFeature (sep, NULL, SEQFEAT_PROT, NULL);
                if (sfp != NULL) {
                  sfp->data.value.ptrvalue = (Pointer) prp;
                  if (cfp->autoUpdatePartials) {
                    sfp->partial = (Boolean) (partial || partial5 || partial3);
                    if (sfp->partial) {
                      SetSeqLocPartial (sfp->location, partial5, partial3);
                    }
                  }
                }
              }
            }
          }
          sep = SeqMgrGetSeqEntryForData (target);
          if (sep != NULL && cfp->autoUpdatePartials) {
            vnp = SeqEntryGetSeqDescr (sep, Seq_descr_molinfo, NULL);
            if (vnp == NULL) {
              vnp = CreateNewDescriptor (sep, Seq_descr_molinfo);
              if (vnp != NULL) {
                mip = MolInfoNew ();
                vnp->data.ptrvalue = (Pointer) mip;
                if (mip != NULL) {
                  mip->biomol = 8;
                  mip->tech = 13;
                }
              }
            }
            if (vnp != NULL) {
              mip = (MolInfoPtr) vnp->data.ptrvalue;
              if (mip != NULL) {
                if (cfp->conceptTransA != NULL) {
                  if (GetStatus (cfp->conceptTransA)) {
                    mip->tech = 13;
                  } else {
                    if (mip->tech == 13) {
                      mip->tech = 8;
                    }
                  }
                }
                if (partial5 && partial3) {
                  mip->completeness = 5;
                } else if (partial5) {
                  mip->completeness = 3;
                } else if (partial3) {
                  mip->completeness = 4;
                } else if (partial) {
                  mip->completeness = 2;
                } else {
                  mip->completeness = 0;
                }
              }
            }
          }
        }
      }
    } else {
      bsp = (BioseqPtr) DialogToPointer (cfp->protseq);
      cfp->usethisbioseq = bsp;
      FeatFormReplaceWithoutUpdateProc (f);
      if (bsp != NULL) {
        prot = SeqEntryNew ();
        if (prot != NULL) {
          prot->choice = 1;
          prot->data.ptrvalue = (Pointer) bsp;
          SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, prot);
          sep = GetBestTopParentForItemID (cfp->input_entityID,
                                           cfp->input_itemID,
                                           cfp->input_itemtype);
          /*
          sep = GetTopSeqEntryForEntityID (cfp->input_entityID);
          */
          descr = ExtractBioSourceAndPubs (sep);
          AddSeqEntryToSeqEntry (sep, prot, TRUE);
          ReplaceBioSourceAndPubs (sep, descr);
          cfp->input_entityID = ObjMgrGetEntityIDForChoice (sep);
          mip = MolInfoNew ();
          if (mip != NULL) {
            mip->biomol = 8;
            mip->tech = 13;
            if (cfp->conceptTransA != NULL) {
              if (! GetStatus (cfp->conceptTransA)) {
                mip->tech = 8;
              }
            }
            if (partial5 && partial3) {
              mip->completeness = 5;
            } else if (partial5) {
              mip->completeness = 3;
            } else if (partial3) {
              mip->completeness = 4;
            } else if (partial) {
              mip->completeness = 2;
            } else {
              mip->completeness = 0;
            }
            vnp = CreateNewDescriptor (prot, Seq_descr_molinfo);
            if (vnp != NULL) {
              vnp->data.ptrvalue = (Pointer) mip;
            }
          }
        }
        if (bsp != NULL && Visible (cfp->protTextGrp)) {
          GetTitle (cfp->protNameText, name, sizeof (name));
          GetTitle (cfp->protDescText, desc, sizeof (desc));
          prp = CreateNewProtRef (name, desc, NULL, NULL);
          if (prp != NULL) {
            sep = SeqMgrGetSeqEntryForData (bsp);
            if (sep != NULL) {
              sfp = CreateNewFeature (sep, NULL, SEQFEAT_PROT, NULL);
              if (sfp != NULL) {
                sfp->data.value.ptrvalue = (Pointer) prp;
                sfp->partial = (Boolean) (partial || partial5 || partial3);
                if (sfp->partial) {
                  SetSeqLocPartial (sfp->location, partial5, partial3);
                }
              }
            }
          }
        }
      }
    }
    SeqLocFree (slp);
    GetRidOfEmptyFeatsDescStrings (cfp->input_entityID, NULL);
    if (GetAppProperty ("InternalNcbiSequin") != NULL) {
      ExtendGeneFeatIfOnMRNA (cfp->input_entityID, NULL);
    }
    ObjMgrSendMsg (OM_MSG_UPDATE, cfp->input_entityID,
                   cfp->input_itemID, cfp->input_itemtype);
  }
}

extern Int2 LIBCALLBACK CdRgnGenFunc (Pointer data)

{
  BioseqPtr         bsp;
  CdRgnFormPtr      cfp;
  HelpMessageFunc   helpfunc;
  OMProcControlPtr  ompcp;
  OMUserDataPtr     omudp;
  SelStructPtr      sel;
  SeqEntryPtr       sep;
  SeqFeatPtr        sfp;
  WindoW            w;

  ompcp = (OMProcControlPtr) data;
  sfp = NULL;
  sep = NULL;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  switch (ompcp->input_itemtype) {
    case OBJ_SEQFEAT :
      sfp = (SeqFeatPtr) ompcp->input_data;
      if (sfp != NULL && sfp->data.choice != SEQFEAT_CDREGION) {
        return OM_MSG_RET_ERROR;
      }
      break;
    case OBJ_BIOSEQ :
      break;
    case OBJ_BIOSEQSET :
      break;
    case 0 :
      break;
    default :
      return OM_MSG_RET_ERROR;
  }
  omudp = ItemAlreadyHasEditor (ompcp->input_entityID, ompcp->input_itemID,
                                ompcp->input_itemtype, ompcp->proc->procid);
  if (omudp != NULL) {
    cfp = (CdRgnFormPtr) omudp->userdata.ptrvalue;
    if (cfp != NULL) {
      Select (cfp->form);
    }
    return OM_MSG_RET_DONE;
  }
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  w = (WindoW) CreateCdRgnForm (-50, -33, "Coding Region", sfp, sep,
                                CdRgnFeatFormActnProc);
  cfp = (CdRgnFormPtr) GetObjectExtra (w);
  if (cfp != NULL) {
    cfp->input_entityID = ompcp->input_entityID;
    cfp->input_itemID = ompcp->input_itemID;
    cfp->input_itemtype = ompcp->input_itemtype;
    cfp->this_itemtype = OBJ_SEQFEAT;
    cfp->this_subtype = FEATDEF_CDS;
    cfp->procid = ompcp->proc->procid;
    cfp->proctype = ompcp->proc->proctype;
    cfp->userkey = OMGetNextUserKey ();
    omudp = ObjMgrAddUserData (ompcp->input_entityID, ompcp->proc->procid,
	                           OMPROC_EDIT, cfp->userkey);
    if (omudp != NULL) {
      omudp->userdata.ptrvalue = (Pointer) cfp;
      omudp->messagefunc = StdVibrantEditorMsgFunc;
    }
    SendMessageToForm (cfp->form, VIB_MSG_INIT);
    if (sfp != NULL) {
      PointerToForm (cfp->form, (Pointer) sfp);
      if (SetClosestParentIfDuplicating ((BaseFormPtr) cfp)) {
        SetBestProteinFeature (cfp, NULL);
        PointerToDialog (cfp->product, NULL);
        PointerToDialog (cfp->protseq, NULL);
      }
    } else {
      /* restored, but locvisited flag protects against accidental CDS */
      SetNewFeatureDefaultInterval ((FeatureFormPtr) cfp);
      sel = ObjMgrGetSelected ();
      if (sel != NULL && sel->next == NULL && sel->entityID == cfp->input_entityID &&
          sel->itemID == cfp->input_itemID && sel->itemtype == cfp->input_itemtype) {
        if (sel->regiontype == 1 && sel->region != NULL) {
          bsp = GetBioseqGivenIDs (cfp->input_entityID, cfp->input_itemID, cfp->input_itemtype);
          if (bsp != NULL && GetBioseqGivenSeqLoc ((SeqLocPtr) sel->region, cfp->input_entityID) == bsp) {
            DoTranslateProtein (cfp);
          }
        }
      }
    }
  }
  Show (w);
  Select (w);
  helpfunc = (HelpMessageFunc) GetAppProperty ("HelpMessageProc");
  if (helpfunc != NULL) {
    helpfunc ("Features", "CDS");
  }
  return OM_MSG_RET_DONE;
}

typedef struct genepage {
  DIALOG_MESSAGE_BLOCK
  TexT          locus;
  TexT          allele;
  TexT          desc;
  TexT          maploc;
  ButtoN        pseudo;
  DialoG        db;
  DialoG        syn;
  GrouP         geneGrp [3];
} GenePage, PNTR GenePagePtr;

typedef struct geneform {
  FEATURE_FORM_BLOCK
  SeqEntryPtr   sep;
  GrouP         pages [NUM_PAGES];
  DialoG        foldertabs;
  Int2          currentPage;
} GeneForm, PNTR GeneFormPtr;

static void GeneRefPtrToGenePage (DialoG d, Pointer data)

{
  GenePagePtr  gpp;
  GeneRefPtr   grp;

  gpp = (GenePagePtr) GetObjectExtra (d);
  grp = (GeneRefPtr) data;
  if (gpp != NULL) {
    if (grp != NULL) {
      SafeSetTitle (gpp->locus, grp->locus);
      SafeSetTitle (gpp->allele, grp->allele);
      SafeSetTitle (gpp->desc, grp->desc);
      SafeSetTitle (gpp->maploc, grp->maploc);
      SafeSetStatus (gpp->pseudo, grp->pseudo);
      PointerToDialog (gpp->db, grp->db);
      PointerToDialog (gpp->syn, grp->syn);
    } else {
      SafeSetTitle (gpp->locus, "");
      SafeSetTitle (gpp->allele, "");
      SafeSetTitle (gpp->desc, "");
      SafeSetTitle (gpp->maploc, "");
      SafeSetStatus (gpp->pseudo, FALSE);
      PointerToDialog (gpp->db, NULL);
      PointerToDialog (gpp->syn, NULL);
    }
  }
}

static Pointer GenePageToGeneRefPtr (DialoG d)

{
  GenePagePtr  gpp;
  GeneRefPtr   grp;

  grp = NULL;
  gpp = (GenePagePtr) GetObjectExtra (d);
  if (gpp != NULL) {
    grp = GeneRefNew ();
    if (grp != NULL) {
      grp->locus = SaveStringFromText (gpp->locus);
      grp->allele = SaveStringFromText (gpp->allele);
      grp->desc = SaveStringFromText (gpp->desc);
      grp->maploc = SaveStringFromText (gpp->maploc);
      grp->pseudo = GetStatus (gpp->pseudo);
      grp->db = DialogToPointer (gpp->db);
      grp->syn = DialogToPointer (gpp->syn);
    }
  }
  return (Pointer) grp;
}

static CharPtr geneTabs [] = {
  "General", "Synonyms", NULL
};

static CharPtr geneTabsXref [] = {
  "General", "Synonyms", "Cross-Refs", NULL
};

static void ChangeGeneSubPage (VoidPtr data, Int2 newval, Int2 oldval)

{
  GenePagePtr  gpp;

  gpp = (GenePagePtr) data;
  if (gpp != NULL) {
    if (oldval >= 0 && oldval <= 2) {
      SafeHide (gpp->geneGrp [oldval]);
    }
    if (newval >= 0 && newval <= 2) {
      SafeShow (gpp->geneGrp [newval]);
    }
    Update ();
  }
}

static DialoG CreateGeneDialog (GrouP h, CharPtr title, GeneRefPtr grp)

{
  GrouP        f;
  GrouP        g;
  GenePagePtr  gpp;
  Char         just;
  GrouP        k;
  GrouP        m;
  GrouP        p;
  GrouP        q;
  GrouP        s;
  Boolean      showXrefs;
  GrouP        t;
  DialoG       tbs;

  p = HiddenGroup (h, 1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  gpp = (GenePagePtr) MemNew (sizeof (GenePage));
  if (gpp != NULL) {

    SetObjectExtra (p, gpp, StdCleanupExtraProc);
    gpp->dialog = (DialoG) p;
    gpp->todialog = GeneRefPtrToGenePage;
    gpp->fromdialog = GenePageToGeneRefPtr;
    gpp->testdialog = NULL;

    if (title != NULL && title [0] != '\0') {
      s = NormalGroup (p, 0, -2, title, systemFont, NULL);
    } else {
      s = HiddenGroup (p, 0, -2, NULL);
    }
    m = HiddenGroup (s, -1, 0, NULL);
    SetGroupSpacing (m, 10, 10);

    showXrefs = FALSE;
    if (grp != NULL && grp->db != NULL) {
      showXrefs = TRUE;
    }

    if (showXrefs) {
      tbs = CreateFolderTabs (m, geneTabsXref, 0, 0, 0,
                              PROGRAM_FOLDER_TAB,
                              ChangeGeneSubPage, (Pointer) gpp);
    } else {
      tbs = CreateFolderTabs (m, geneTabs, 0, 0, 0,
                              PROGRAM_FOLDER_TAB,
                              ChangeGeneSubPage, (Pointer) gpp);
    }
    k = HiddenGroup (m, 0, 0, NULL);

    gpp->geneGrp [0] = HiddenGroup (k, -1, 0, NULL);
    f = HiddenGroup (gpp->geneGrp [0], 2, 0, NULL);
    StaticPrompt (f, "Locus", 0, dialogTextHeight, programFont, 'l');
    gpp->locus = DialogText (f, "", 15, NULL);
    StaticPrompt (f, "Allele", 0, dialogTextHeight, programFont, 'l');
    gpp->allele = DialogText (f, "", 15, NULL);
    StaticPrompt (f, "Description", 0, dialogTextHeight, programFont, 'l');
    gpp->desc = DialogText (f, "", 15, NULL);
    StaticPrompt (f, "Map Location", 0, dialogTextHeight, programFont, 'l');
    gpp->maploc = DialogText (f, "", 15, NULL);

    gpp->pseudo = CheckBox (gpp->geneGrp [0], "PseudoGene", NULL);
    AlignObjects (ALIGN_CENTER, (HANDLE) f, (HANDLE) gpp->pseudo, NULL);

    gpp->geneGrp [1] = HiddenGroup (k, -1, 0, NULL);
    g = HiddenGroup (gpp->geneGrp [1], 0, 2, NULL);
    StaticPrompt (g, "Synonyms", 0, 0, programFont, 'c');
    gpp->syn = CreateVisibleStringDialog (g, 3, -1, 15);
    Hide (gpp->geneGrp [1]);

    gpp->geneGrp [2] = HiddenGroup (k, -1, 0, NULL);
    if (showXrefs) {
      q = HiddenGroup (gpp->geneGrp [2], -1, 0, NULL);
      if (GetAppProperty ("ReadOnlyDbTags") == NULL) {
        just = 'c';
      } else {
        just = 'l';
        StaticPrompt (q, "This page is read-only", 15 * stdCharWidth, 0, programFont, 'c');
      }
      t = HiddenGroup (q, 2, 0, NULL);
      StaticPrompt (t, "Database", 7 * stdCharWidth, 0, programFont, just);
      StaticPrompt (t, "Object ID", 8 * stdCharWidth, 0, programFont, just);
      gpp->db = CreateDbtagDialog (q, 3, -1, 7, 8);
    }
    Hide (gpp->geneGrp [2]);

    AlignObjects (ALIGN_CENTER, (HANDLE) tbs,
                  (HANDLE) gpp->geneGrp [0], (HANDLE) gpp->geneGrp [1],
                  (HANDLE) gpp->geneGrp [2], NULL);
  }

  return (DialoG) p;
}

static void SetGeneImportExportItems (GeneFormPtr gfp)

{
  IteM  exportItm;
  IteM  importItm;

  if (gfp != NULL) {
    importItm = FindFormMenuItem ((BaseFormPtr) gfp, VIB_MSG_IMPORT);
    exportItm = FindFormMenuItem ((BaseFormPtr) gfp, VIB_MSG_EXPORT);
    switch (gfp->currentPage) {
      case GENE_PAGE :
        SafeSetTitle (importItm, "Import...");
        SafeSetTitle (exportItm, "Export...");
        SafeDisable (importItm);
        SafeDisable (exportItm);
        break;
      case COMMON_PAGE :
        SafeSetTitle (importItm, "Import...");
        SafeSetTitle (exportItm, "Export...");
        SafeDisable (importItm);
        SafeDisable (exportItm);
        break;
      case LOCATION_PAGE :
        SafeSetTitle (importItm, "Import SeqLoc...");
        SafeSetTitle (exportItm, "Export SeqLoc...");
        SafeEnable (importItm);
        SafeEnable (exportItm);
        break;
      default :
        break;
    }
  }
}

static void ChangeGenePage (VoidPtr data, Int2 newval, Int2 oldval)

{
  GeneFormPtr   gfp;

  gfp = (GeneFormPtr) data;
  if (gfp != NULL) {
    gfp->currentPage = newval;
    SafeHide (gfp->pages [oldval]);
    SafeShow (gfp->pages [newval]);
    switch (newval) {
      case GENE_PAGE :
        SendMessageToDialog (gfp->data, VIB_MSG_ENTER);
        break;
      case COMMON_PAGE :
        break;
      case LOCATION_PAGE :
        SendMessageToDialog (gfp->location, VIB_MSG_ENTER);
        break;
      default :
        break;
    }
    SetGeneImportExportItems (gfp);
    Update ();
  }
}

static Boolean ImportGeneForm (ForM f, CharPtr filename)

{
  GeneFormPtr  gfp;

  gfp = (GeneFormPtr) GetObjectExtra (f);
  if (gfp != NULL) {
    switch (gfp->currentPage) {
      case LOCATION_PAGE :
        return ImportDialog (gfp->location, filename);
      default :
        break;
    }
  }
  return FALSE;
}

static Boolean ExportGeneForm (ForM f, CharPtr filename)

{
  GeneFormPtr  gfp;

  gfp = (GeneFormPtr) GetObjectExtra (f);
  if (gfp != NULL) {
    switch (gfp->currentPage) {
      case LOCATION_PAGE :
        return ExportDialog (gfp->location, filename);
      default :
        break;
    }
  }
  return FALSE;
}

static CharPtr  geneFormTabs [] = {
  "Gene", "Properties", "Location", NULL
};

static void GeneFormMessage (ForM f, Int2 mssg)

{
  GeneFormPtr  gfp;

  gfp = (GeneFormPtr) GetObjectExtra (f);
  if (gfp != NULL) {
    switch (mssg) {
      case VIB_MSG_INIT :
        StdInitFeatFormProc (f);
        break;
      case VIB_MSG_IMPORT :
        ImportGeneForm (f, NULL);
        break;
      case VIB_MSG_EXPORT :
        ExportGeneForm (f, NULL);
        break;
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
        if (gfp->currentPage == LOCATION_PAGE) {
          PointerToDialog (gfp->location, NULL);
        } else {
          StdDeleteTextProc (NULL);
        }
        break;
      default :
        if (gfp->appmessage != NULL) {
          gfp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static void GeneFormActivate (WindoW w)

{
  GeneFormPtr  gfp;

  gfp = (GeneFormPtr) GetObjectExtra (w);
  if (gfp != NULL) {
    if (gfp->activate != NULL) {
      gfp->activate (w);
    }
    SetGeneImportExportItems (gfp);
  }
}

extern ForM CreateGeneForm (Int2 left, Int2 top, CharPtr title,
                            SeqFeatPtr sfp, SeqEntryPtr sep,
                            FormActnFunc actproc)

{
  ButtoN             b;
  GrouP              c;
  GrouP              g;
  GeneFormPtr        gfp;
  GeneRefPtr         grp;
  GrouP              h;
  GrouP              s;
  StdEditorProcsPtr  sepp;
  WindoW             w;

  w = NULL;
  gfp = (GeneFormPtr) MemNew (sizeof (GeneForm));
  if (gfp != NULL) {
    w = FixedWindow (left, top, -10, -10, title, StdCloseWindowProc);
    SetObjectExtra (w, gfp, StdFeatFormCleanupProc);
    gfp->form = (ForM) w;
    gfp->actproc = actproc;
    gfp->toform = StdSeqFeatPtrToFeatFormProc;
    gfp->fromform = NULL;
    gfp->formmessage = GeneFormMessage;
    gfp->testform = NULL;
    gfp->importform = ImportGeneForm;
    gfp->exportform = ExportGeneForm;

#ifndef WIN_MAC
    CreateStdEditorFormMenus (w);
#endif

    gfp->activate = NULL;
    sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
    if (sepp != NULL) {
      gfp->activate = sepp->activateForm;
      gfp->appmessage = sepp->handleMessages;
    }
    SetActivate (w, GeneFormActivate);

    g = HiddenGroup (w, -1, 0, NULL);
    SetGroupSpacing (g, 3, 10);

    gfp->sep = sep;
    gfp->foldertabs = CreateFolderTabs (g, geneFormTabs, GENE_PAGE,
                                        0, 0, SYSTEM_FOLDER_TAB,
                                        ChangeGenePage, (Pointer) gfp);
    gfp->currentPage = GENE_PAGE;

    h = HiddenGroup (g, 0, 0, NULL);

    s = HiddenGroup (h, -1, 0, NULL);
    grp = NULL;
    if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE) {
      grp = sfp->data.value.ptrvalue;
    }
    gfp->data = CreateGeneDialog (s, NULL, grp);
    gfp->pages [GENE_PAGE] = s;
    Hide (gfp->pages [GENE_PAGE]);

    s = HiddenGroup (h, -1, 0, NULL);
    CreateCommonFeatureGroup (s, (FeatureFormPtr) gfp, sfp, FALSE, TRUE);
    gfp->pages [COMMON_PAGE] = s;
    Hide (gfp->pages [COMMON_PAGE]);

    s = HiddenGroup (h, -1, 0, NULL);
    gfp->location = CreateIntervalEditorDialogEx (s, NULL, 4, 2, sep, TRUE, FALSE,
                                                  TRUE, TRUE, FALSE,
                                                  (FeatureFormPtr) gfp,
                                                  StdFeatIntEdPartialCallback);
    gfp->pages [LOCATION_PAGE] = s;
    Hide (gfp->pages [LOCATION_PAGE]);

    AlignObjects (ALIGN_CENTER, (HANDLE) gfp->pages [GENE_PAGE],
                  (HANDLE) gfp->pages [COMMON_PAGE], (HANDLE) gfp->pages [LOCATION_PAGE],
                  NULL);
    AlignObjects (ALIGN_CENTER, (HANDLE) gfp->foldertabs, (HANDLE) h, NULL);

    c = HiddenGroup (w, 2, 0, NULL);
    b = PushButton (c, "Accept", StdFeatFormAcceptButtonProc);
    SetObjectExtra (b, gfp, NULL);
    PushButton (c, "Cancel", StdCancelButtonProc);
    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
    RealizeWindow (w);

    SendMessageToDialog (gfp->data, VIB_MSG_INIT);
    SendMessageToDialog (gfp->location, VIB_MSG_INIT);
    Show (gfp->pages [gfp->currentPage]);
    SendMessageToDialog (gfp->data, VIB_MSG_ENTER);
    Update ();
  }
  return (ForM) w;
}

static Boolean CDSPseudoByGeneFunc (GatherContextPtr gcp, Boolean pseudo)

{
  GBQualPtr   gbq;
  GBQualPtr   next;
  GBQualPtr PNTR  prev;
  GBQualPtr   qual;
  SeqFeatPtr  sfp;
  SeqLocPtr   slp;

  if (gcp == NULL) return TRUE;
  slp = (SeqLocPtr) gcp->userdata;
  if (slp == NULL) return TRUE;
  if (gcp->thistype == OBJ_SEQFEAT) {
    sfp = (SeqFeatPtr) gcp->thisitem;
    if (sfp != NULL && sfp->data.choice == SEQFEAT_CDREGION) {
      if (SeqLocAinB (sfp->location, slp) >= 0) {
        gbq = sfp->qual;
        prev = (GBQualPtr PNTR) &(sfp->qual);
        while (gbq != NULL) {
          next = gbq->next;
          if (StringICmp (gbq->qual, "pseudo") == 0) {
            if (pseudo) return TRUE;
            *(prev) = gbq->next;
            gbq->next = NULL;
            GBQualFree (gbq);
          } else {
            prev = (GBQualPtr PNTR) &(gbq->next);
          }
          gbq = next;
        }
        if (! pseudo) return TRUE;
        qual = GBQualNew ();
        if (qual != NULL) {
          qual->qual = StringSave ("pseudo");
          qual->val = StringSave ("");
          gbq = sfp->qual;
          if (gbq != NULL) {
            while (gbq->next != NULL) {
              gbq = gbq->next;
            }
            gbq->next = qual;
          } else {
            sfp->qual = qual;
          }
        }
      }
    }
  }
  return TRUE;
}

static Boolean CDSPseudoOnFunc (GatherContextPtr gcp)

{
  return CDSPseudoByGeneFunc (gcp, TRUE);
}

static Boolean CDSPseudoOffFunc (GatherContextPtr gcp)

{
  return CDSPseudoByGeneFunc (gcp, FALSE);
}

static void GeneFeatFormActnProc (ForM f)

{
  GeneFormPtr  gfp;
  GenePagePtr  gpp;
  GatherScope  gs;
  Boolean      pseudo = FALSE;
  SeqEntryPtr  sep;
  SeqLocPtr    slp = NULL;

  gfp = (GeneFormPtr) GetObjectExtra (f);
  if (gfp != NULL) {
    slp = DialogToPointer (gfp->location);
    gpp = GetObjectExtra (gfp->data);
    if (gpp != NULL) {
      pseudo = GetStatus (gpp->pseudo);
    }
  }
  if (FeatFormReplaceWithoutUpdateProc (f)) {
    if (gfp != NULL) {
      sep = GetBestTopParentForItemID (gfp->input_entityID,
                                       gfp->input_itemID,
                                       gfp->input_itemtype);
      if (slp != NULL) {
        MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
        gs.seglevels = 1;
        gs.get_feats_location = FALSE;
	    MemSet ((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
	    gs.ignore[OBJ_BIOSEQ] = FALSE;
	    gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
	    gs.ignore[OBJ_SEQFEAT] = FALSE;
	    gs.ignore[OBJ_SEQANNOT] = FALSE;
        gs.scope = sep;
        if (pseudo) {
          GatherEntity (gfp->input_entityID, (Pointer) slp, CDSPseudoOnFunc, &gs);
        } else {
          GatherEntity (gfp->input_entityID, (Pointer) slp, CDSPseudoOffFunc, &gs);
        }
      }
      GetRidOfEmptyFeatsDescStrings (gfp->input_entityID, NULL);
      if (GetAppProperty ("InternalNcbiSequin") != NULL) {
        ExtendGeneFeatIfOnMRNA (gfp->input_entityID, NULL);
      }
      ObjMgrSendMsg (OM_MSG_UPDATE, gfp->input_entityID, 0, 0);
    }
  }
  SeqLocFree (slp);
}

extern Int2 LIBCALLBACK GeneGenFunc (Pointer data)

{
  GeneFormPtr       gfp;
  HelpMessageFunc   helpfunc;
  OMProcControlPtr  ompcp;
  OMUserDataPtr     omudp;
  SeqEntryPtr       sep;
  SeqFeatPtr        sfp;
  WindoW            w;

  ompcp = (OMProcControlPtr) data;
  sfp = NULL;
  sep = NULL;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  switch (ompcp->input_itemtype) {
    case OBJ_SEQFEAT :
      sfp = (SeqFeatPtr) ompcp->input_data;
      if (sfp != NULL && sfp->data.choice != SEQFEAT_GENE) {
        return OM_MSG_RET_ERROR;
      }
      break;
    case OBJ_BIOSEQ :
      break;
    case OBJ_BIOSEQSET :
      break;
    case 0 :
      break;
    default :
      return OM_MSG_RET_ERROR;
  }
  omudp = ItemAlreadyHasEditor (ompcp->input_entityID, ompcp->input_itemID,
                                ompcp->input_itemtype, ompcp->proc->procid);
  if (omudp != NULL) {
    gfp = (GeneFormPtr) omudp->userdata.ptrvalue;
    if (gfp != NULL) {
      Select (gfp->form);
    }
    return OM_MSG_RET_DONE;
  }
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  w = (WindoW) CreateGeneForm (-50, -33, "Gene", sfp, sep,
                               GeneFeatFormActnProc);
  gfp = (GeneFormPtr) GetObjectExtra (w);
  if (gfp != NULL) {
    gfp->input_entityID = ompcp->input_entityID;
    gfp->input_itemID = ompcp->input_itemID;
    gfp->input_itemtype = ompcp->input_itemtype;
    gfp->this_itemtype = OBJ_SEQFEAT;
    gfp->this_subtype = FEATDEF_GENE;
    gfp->procid = ompcp->proc->procid;
    gfp->proctype = ompcp->proc->proctype;
    gfp->userkey = OMGetNextUserKey ();
    omudp = ObjMgrAddUserData (ompcp->input_entityID, ompcp->proc->procid,
	                           OMPROC_EDIT, gfp->userkey);
    if (omudp != NULL) {
      omudp->userdata.ptrvalue = (Pointer) gfp;
      omudp->messagefunc = StdVibrantEditorMsgFunc;
    }
    SendMessageToForm (gfp->form, VIB_MSG_INIT);
    if (sfp != NULL) {
      PointerToForm (gfp->form, (Pointer) sfp);
      SetClosestParentIfDuplicating ((BaseFormPtr) gfp);
    } else {
      SetNewFeatureDefaultInterval ((FeatureFormPtr) gfp);
    }
  }
  Show (w);
  Select (w);
  helpfunc = (HelpMessageFunc) GetAppProperty ("HelpMessageProc");
  if (helpfunc != NULL) {
    helpfunc ("Features", "Gene");
  }
  return OM_MSG_RET_DONE;
}

typedef struct protpage {
  DIALOG_MESSAGE_BLOCK
  DialoG        name;
  TexT          desc;
  DialoG        ec;
  DialoG        activity;
  DialoG        db;
  PopuP         processed;
  GrouP         protGrp [4];
} ProtPage, PNTR ProtPagePtr;

typedef struct protform {
  FEATURE_FORM_BLOCK
  SeqEntryPtr   sep;
  GrouP         pages [NUM_PAGES];
  DialoG        foldertabs;
  Int2          currentPage;
} ProtForm, PNTR ProtFormPtr;

static void ProtRefPtrToProtPage (DialoG d, Pointer data)

{
  ProtPagePtr  ppp;
  ProtRefPtr   prp;

  ppp = (ProtPagePtr) GetObjectExtra (d);
  prp = (ProtRefPtr) data;
  if (ppp != NULL) {
    if (prp != NULL) {
      PointerToDialog (ppp->name, prp->name);
      SafeSetTitle (ppp->desc, prp->desc);
      PointerToDialog (ppp->ec, prp->ec);
      PointerToDialog (ppp->activity, prp->activity);
      PointerToDialog (ppp->db, prp->db);
      SafeSetValue (ppp->processed, prp->processed + 1);
    } else {
      PointerToDialog (ppp->name, NULL);
      SafeSetTitle (ppp->desc, "");
      PointerToDialog (ppp->ec, NULL);
      PointerToDialog (ppp->activity, NULL);
      PointerToDialog (ppp->db, NULL);
      SafeSetValue (ppp->processed, 1);
    }
  }
}

static Pointer ProtPageToProtRefPtr (DialoG d)

{
  ProtPagePtr  ppp;
  ProtRefPtr   prp;

  prp = NULL;
  ppp = (ProtPagePtr) GetObjectExtra (d);
  if (ppp != NULL) {
    prp = ProtRefNew ();
    if (prp != NULL) {
      prp->name = DialogToPointer (ppp->name);
      prp->desc = SaveStringFromText (ppp->desc);
      prp->ec = DialogToPointer (ppp->ec);
      prp->activity = DialogToPointer (ppp->activity);
      prp->db = DialogToPointer (ppp->db);
      prp->processed = GetValue (ppp->processed) - 1;
    }
  }
  return (Pointer) prp;
}

static CharPtr protTabs [] = {
  "General", "E.C. Number", "Activity", NULL
};

static CharPtr protTabsXref [] = {
  "General", "E.C. Number", "Activity", "Cross-Refs", NULL
};

static void ChangeProtSubPage (VoidPtr data, Int2 newval, Int2 oldval)

{
  ProtPagePtr  ppp;

  ppp = (ProtPagePtr) data;
  if (ppp != NULL) {
    if (oldval >= 0 && oldval <= 3) {
      SafeHide (ppp->protGrp [oldval]);
    }
    if (newval >= 0 && newval <= 3) {
      SafeShow (ppp->protGrp [newval]);
    }
    Update ();
  }
}

static DialoG CreateProtDialog (GrouP h, CharPtr title, ProtRefPtr prp)

{
  GrouP        f;
  GrouP        g;
  Char         just;
  GrouP        k;
  GrouP        m;
  GrouP        p;
  ProtPagePtr  ppp;
  GrouP        q;
  GrouP        r;
  GrouP        s;
  Boolean      showXrefs;
  GrouP        t;
  DialoG       tbs;

  p = HiddenGroup (h, 1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  ppp = (ProtPagePtr) MemNew (sizeof (ProtPage));
  if (ppp != NULL) {

    SetObjectExtra (p, ppp, StdCleanupExtraProc);
    ppp->dialog = (DialoG) p;
    ppp->todialog = ProtRefPtrToProtPage;
    ppp->fromdialog = ProtPageToProtRefPtr;
    ppp->testdialog = NULL;

    if (title != NULL && title [0] != '\0') {
      s = NormalGroup (p, 0, -2, title, systemFont, NULL);
    } else {
      s = HiddenGroup (p, 0, -2, NULL);
    }
    m = HiddenGroup (s, -1, 0, NULL);
    SetGroupSpacing (m, 10, 10);

    showXrefs = FALSE;
    if (prp != NULL && prp->db != NULL) {
      showXrefs = TRUE;
    }

    if (showXrefs) {
      tbs = CreateFolderTabs (m, protTabsXref, 0, 0, 0,
                              PROGRAM_FOLDER_TAB,
                              ChangeProtSubPage, (Pointer) ppp);
    } else {
      tbs = CreateFolderTabs (m, protTabs, 0, 0, 0,
                              PROGRAM_FOLDER_TAB,
                              ChangeProtSubPage, (Pointer) ppp);
    }
    k = HiddenGroup (m, 0, 0, NULL);

    ppp->protGrp [0] = HiddenGroup (k, -1, 0, NULL);
    g = HiddenGroup (ppp->protGrp [0], 0, 10, NULL);
    StaticPrompt (g, "Protein Names", 0, 0, programFont, 'c');
    ppp->name = CreateVisibleStringDialog (g, 3, -1, 15);

    f = HiddenGroup (ppp->protGrp [0], 0, 4, NULL);
    StaticPrompt (f, "Description", 0, dialogTextHeight, programFont, 'c');
    ppp->desc = DialogText (f, "", 20, NULL);

    r = HiddenGroup (ppp->protGrp [0], 2, 0, NULL);
    StaticPrompt (r, "Processing", 0, dialogTextHeight, programFont, 'l');
    ppp->processed = PopupList (r, TRUE, NULL);
    PopupItem (ppp->processed, " ");
    PopupItem (ppp->processed, "Pre-protein");
    PopupItem (ppp->processed, "Mature");
    PopupItem (ppp->processed, "Signal peptide");
    PopupItem (ppp->processed, "Transit peptide");
    SetValue (ppp->processed, 1);

    ppp->protGrp [1] = HiddenGroup (k, 0, 10, NULL);
    StaticPrompt (ppp->protGrp [1], "Enzyme Commission Number", 0, 0, programFont, 'c');
    ppp->ec = CreateVisibleStringDialog (ppp->protGrp [1], 3, -1, 15);
    Hide (ppp->protGrp [1]);

    ppp->protGrp [2] = HiddenGroup (k, 0, 10, NULL);
    StaticPrompt (ppp->protGrp [2], "Activity", 0, 0, programFont, 'c');
    ppp->activity = CreateVisibleStringDialog (ppp->protGrp [2], 3, -1, 15);
    Hide (ppp->protGrp [2]);

    ppp->protGrp [3] = HiddenGroup (k, -1, 0, NULL);
    if (showXrefs) {
      q = HiddenGroup (ppp->protGrp [3], -1, 0, NULL);
      if (GetAppProperty ("ReadOnlyDbTags") == NULL) {
        just = 'c';
      } else {
        just = 'l';
        StaticPrompt (q, "This page is read-only", 15 * stdCharWidth, 0, programFont, 'c');
      }
      t = HiddenGroup (q, 2, 0, NULL);
      StaticPrompt (t, "Database", 7 * stdCharWidth, 0, programFont, just);
      StaticPrompt (t, "Object ID", 8 * stdCharWidth, 0, programFont, just);
      ppp->db = CreateDbtagDialog (q, 3, -1, 7, 8);
    }
    Hide (ppp->protGrp [3]);

    AlignObjects (ALIGN_CENTER, (HANDLE) ppp->protGrp [0],
                  (HANDLE) ppp->protGrp [1], (HANDLE) ppp->protGrp [2],
                  (HANDLE) ppp->protGrp [3], (HANDLE) tbs,
                  (HANDLE) g, (HANDLE) r, NULL);
  }

  return (DialoG) p;
}

static void SetProtImportExportItems (ProtFormPtr pfp)

{
  IteM  exportItm;
  IteM  importItm;

  if (pfp != NULL) {
    importItm = FindFormMenuItem ((BaseFormPtr) pfp, VIB_MSG_IMPORT);
    exportItm = FindFormMenuItem ((BaseFormPtr) pfp, VIB_MSG_EXPORT);
    switch (pfp->currentPage) {
      case PROT_PAGE :
        SafeSetTitle (importItm, "Import...");
        SafeSetTitle (exportItm, "Export...");
        SafeDisable (importItm);
        SafeDisable (exportItm);
        break;
      case COMMON_PAGE :
        SafeSetTitle (importItm, "Import...");
        SafeSetTitle (exportItm, "Export...");
        SafeDisable (importItm);
        SafeDisable (exportItm);
        break;
      case LOCATION_PAGE :
        SafeSetTitle (importItm, "Import SeqLoc...");
        SafeSetTitle (exportItm, "Export SeqLoc...");
        SafeEnable (importItm);
        SafeEnable (exportItm);
        break;
      default :
        break;
    }
  }
}

static void ChangeProtPage (VoidPtr data, Int2 newval, Int2 oldval)

{
  ProtFormPtr   pfp;

  pfp = (ProtFormPtr) data;
  if (pfp != NULL) {
    pfp->currentPage = newval;
    SafeHide (pfp->pages [oldval]);
    SafeShow (pfp->pages [newval]);
    switch (newval) {
      case PROT_PAGE :
        SendMessageToDialog (pfp->data, VIB_MSG_ENTER);
        break;
      case COMMON_PAGE :
        break;
      case LOCATION_PAGE :
        SendMessageToDialog (pfp->location, VIB_MSG_ENTER);
        break;
      default :
        break;
    }
    SetProtImportExportItems (pfp);
    Update ();
  }
}

static Boolean ImportProtForm (ForM f, CharPtr filename)

{
  ProtFormPtr  pfp;

  pfp = (ProtFormPtr) GetObjectExtra (f);
  if (pfp != NULL) {
    switch (pfp->currentPage) {
      case LOCATION_PAGE :
        return ImportDialog (pfp->location, filename);
      default :
        break;
    }
  }
  return FALSE;
}

static Boolean ExportProtForm (ForM f, CharPtr filename)

{
  ProtFormPtr  pfp;

  pfp = (ProtFormPtr) GetObjectExtra (f);
  if (pfp != NULL) {
    switch (pfp->currentPage) {
      case LOCATION_PAGE :
        return ExportDialog (pfp->location, filename);
      default :
        break;
    }
  }
  return FALSE;
}

static CharPtr  protFormTabs [] = {
  "Protein", "Properties", "Location", NULL
};

static void ProtFormMessage (ForM f, Int2 mssg)

{
  ProtFormPtr  pfp;

  pfp = (ProtFormPtr) GetObjectExtra (f);
  if (pfp != NULL) {
    switch (mssg) {
      case VIB_MSG_INIT :
        StdInitFeatFormProc (f);
        break;
      case VIB_MSG_IMPORT :
        ImportProtForm (f, NULL);
        break;
      case VIB_MSG_EXPORT :
        ExportProtForm (f, NULL);
        break;
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
        if (pfp->currentPage == LOCATION_PAGE) {
          PointerToDialog (pfp->location, NULL);
        } else {
          StdDeleteTextProc (NULL);
        }
        break;
      default :
        if (pfp->appmessage != NULL) {
          pfp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static void ProtFormActivate (WindoW w)

{
  ProtFormPtr  pfp;

  pfp = (ProtFormPtr) GetObjectExtra (w);
  if (pfp != NULL) {
    if (pfp->activate != NULL) {
      pfp->activate (w);
    }
    SetProtImportExportItems (pfp);
  }
}

extern ForM CreateProtForm (Int2 left, Int2 top, CharPtr title,
                            SeqFeatPtr sfp, SeqEntryPtr sep,
                            FormActnFunc actproc)

{
  ButtoN             b;
  GrouP              c;
  GrouP              g;
  GrouP              h;
  ProtFormPtr        pfp;
  ProtRefPtr         prp;
  GrouP              s;
  StdEditorProcsPtr  sepp;
  WindoW             w;

  w = NULL;
  pfp = (ProtFormPtr) MemNew (sizeof (ProtForm));
  if (pfp != NULL) {
    w = FixedWindow (left, top, -10, -10, title, StdCloseWindowProc);
    SetObjectExtra (w, pfp, StdFeatFormCleanupProc);
    pfp->form = (ForM) w;
    pfp->actproc = actproc;
    pfp->toform = StdSeqFeatPtrToFeatFormProc;
    pfp->fromform = NULL;
    pfp->formmessage = ProtFormMessage;
    pfp->testform = NULL;
    pfp->importform = ImportProtForm;
    pfp->exportform = ExportProtForm;

#ifndef WIN_MAC
    CreateStdEditorFormMenus (w);
#endif

    pfp->activate = NULL;
    sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
    if (sepp != NULL) {
      pfp->activate = sepp->activateForm;
      pfp->appmessage = sepp->handleMessages;
    }
    SetActivate (w, ProtFormActivate);

    g = HiddenGroup (w, -1, 0, NULL);
    SetGroupSpacing (g, 3, 10);

    pfp->sep = sep;
    pfp->foldertabs = CreateFolderTabs (g, protFormTabs, PROT_PAGE,
                                        0, 0, SYSTEM_FOLDER_TAB,
                                        ChangeProtPage, (Pointer) pfp);
    pfp->currentPage = PROT_PAGE;

    h = HiddenGroup (g, 0, 0, NULL);

    s = HiddenGroup (h, -1, 0, NULL);
    prp = NULL;
    if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT) {
      prp = sfp->data.value.ptrvalue;
    }
    pfp->data = CreateProtDialog (s, NULL, prp);
    pfp->pages [PROT_PAGE] = s;
    Hide (pfp->pages [PROT_PAGE]);

    s = HiddenGroup (h, -1, 0, NULL);
    CreateCommonFeatureGroup (s, (FeatureFormPtr) pfp, sfp, FALSE, TRUE);
    pfp->pages [COMMON_PAGE] = s;
    Hide (pfp->pages [COMMON_PAGE]);

    s = HiddenGroup (h, -1, 0, NULL);
    pfp->location = CreateIntervalEditorDialogEx (s, NULL, 4, 2, sep, FALSE, TRUE,
                                                  TRUE, TRUE, FALSE,
                                                  (FeatureFormPtr) pfp,
                                                  StdFeatIntEdPartialCallback);
    pfp->pages [LOCATION_PAGE] = s;
    Hide (pfp->pages [LOCATION_PAGE]);

    AlignObjects (ALIGN_CENTER, (HANDLE) pfp->pages [PROT_PAGE],
                  (HANDLE) pfp->pages [COMMON_PAGE], (HANDLE) pfp->pages [LOCATION_PAGE],
                  NULL);
    AlignObjects (ALIGN_CENTER, (HANDLE) pfp->foldertabs, (HANDLE) h, NULL);

    c = HiddenGroup (w, 2, 0, NULL);
    b = PushButton (c, "Accept", StdFeatFormAcceptButtonProc);
    SetObjectExtra (b, pfp, NULL);
    PushButton (c, "Cancel", StdCancelButtonProc);
    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
    RealizeWindow (w);

    SendMessageToDialog (pfp->data, VIB_MSG_INIT);
    SendMessageToDialog (pfp->location, VIB_MSG_INIT);
    Show (pfp->pages [pfp->currentPage]);
    SendMessageToDialog (pfp->data, VIB_MSG_ENTER);
    Update ();
  }
  return (ForM) w;
}

extern Int2 LIBCALLBACK ProtGenFunc (Pointer data)

{
  ProtFormPtr       pfp;
  HelpMessageFunc   helpfunc;
  ObjMgrPtr         omp;
  OMProcControlPtr  ompcp;
  ObjMgrTypePtr     omtp;
  OMUserDataPtr     omudp;
  ProtPagePtr       ppp;
  ObjMgrProcPtr     proc;
  SeqEntryPtr       sep;
  SeqFeatPtr        sfp;
  Uint2             subtype;
  WindoW            w;

  ompcp = (OMProcControlPtr) data;
  sfp = NULL;
  sep = NULL;
  subtype = FEATDEF_PROT;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  proc = ompcp->proc;
  switch (ompcp->input_itemtype) {
    case OBJ_SEQFEAT :
      sfp = (SeqFeatPtr) ompcp->input_data;
      if (sfp != NULL && sfp->data.choice != SEQFEAT_PROT) {
        return OM_MSG_RET_ERROR;
      }
      break;
    case OBJ_BIOSEQ :
      break;
    case OBJ_BIOSEQSET :
      break;
    case 0 :
      break;
    default :
      return OM_MSG_RET_ERROR;
  }
  omudp = ItemAlreadyHasEditor (ompcp->input_entityID, ompcp->input_itemID,
                                ompcp->input_itemtype, ompcp->proc->procid);
  if (omudp != NULL) {
    pfp = (ProtFormPtr) omudp->userdata.ptrvalue;
    if (pfp != NULL) {
      Select (pfp->form);
    }
    return OM_MSG_RET_DONE;
  }
  if (sfp == NULL && proc->subinputtype > 0) {
    subtype = proc->subinputtype;
  }
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  w = (WindoW) CreateProtForm (-50, -33, "Protein", sfp, sep,
                               StdFeatFormActnProc);
  pfp = (ProtFormPtr) GetObjectExtra (w);
  if (pfp != NULL) {
    pfp->input_entityID = ompcp->input_entityID;
    pfp->input_itemID = ompcp->input_itemID;
    pfp->input_itemtype = ompcp->input_itemtype;
    pfp->this_itemtype = OBJ_SEQFEAT;
    pfp->this_subtype = subtype;
    pfp->procid = ompcp->proc->procid;
    pfp->proctype = ompcp->proc->proctype;
    pfp->userkey = OMGetNextUserKey ();
    omudp = ObjMgrAddUserData (ompcp->input_entityID, ompcp->proc->procid,
	                           OMPROC_EDIT, pfp->userkey);
    if (omudp != NULL) {
      omudp->userdata.ptrvalue = (Pointer) pfp;
      omudp->messagefunc = StdVibrantEditorMsgFunc;
    }
    if (sfp != NULL) {
      omp = ObjMgrGet ();
      if (omp != NULL) {
        omtp = ObjMgrTypeFind (omp, OBJ_SEQFEAT, NULL, NULL);
        if (omtp != NULL && omtp->subtypefunc != NULL) {
          pfp->this_subtype = (*(omtp->subtypefunc)) (sfp);
        }
      }
    }
    SendMessageToForm (pfp->form, VIB_MSG_INIT);
    if (sfp != NULL) {
      PointerToForm (pfp->form, (Pointer) sfp);
      SetClosestParentIfDuplicating ((BaseFormPtr) pfp);
    } else {
      ppp = (ProtPagePtr) GetObjectExtra (pfp->data);
      if (ppp != NULL) {
        if (subtype == FEATDEF_PROT) {
          SetValue (ppp->processed, 1);
        } else if (subtype >= FEATDEF_preprotein && subtype <= FEATDEF_transit_peptide_aa) {
          SetValue (ppp->processed, subtype - FEATDEF_preprotein + 2);
        } else {
          SetValue (ppp->processed, 1);
        }
      }
      SetNewFeatureDefaultInterval ((FeatureFormPtr) pfp);
    }
  }
  Show (w);
  Select (w);
  helpfunc = (HelpMessageFunc) GetAppProperty ("HelpMessageProc");
  if (helpfunc != NULL) {
    helpfunc ("Features", "Protein");
  }
  return OM_MSG_RET_DONE;
}

typedef struct rnaform {
  FEATURE_FORM_BLOCK
  SeqEntryPtr   sep;
  GrouP         pages [NUM_PAGES];
  DialoG        foldertabs;
  Int2          currentPage;
} RnaForm, PNTR RnaFormPtr;

typedef struct rnapage {
  DIALOG_MESSAGE_BLOCK
  PopuP         type;
  ButtoN        pseudo;
  GrouP         nameGrp;
  PrompT        rrnaPrompt;
  TexT          name;
  GrouP         trnaGrp;
  PopuP         AAitem;
  DialoG        codons;
  DialoG        anticodon;
  GrouP         trnaPages [3];
  RnaFormPtr    rfp;
} RnaPage, PNTR RnaPagePtr;

static ENUM_ALIST(rna_type_alist)
  {" ",            0},
  {"preRna",       1},
  {"mRNA",         2},
  {"tRNA",         3},
  {"rRNA",         4},
  {"snRNA",        5},
  {"scRNA",        6},
  {"Other RNA",  255},
END_ENUM_ALIST

static Uint1 check_rna_type (Uint1 type)

{
  if (type > 6 && type != 255) return 0;
  return type;
}

/* trna aa conversion modified from Tatiana Tatusov's code */

static void RnaRefPtrToRnaPage (DialoG d, Pointer data)

{
  Uint1           aa;
  Uint1           codon [4];
  Uint1           from;
  ValNodePtr      head;
  Uint1           i;
  Int2            j;
  RnaPagePtr      rpp;
  RnaRefPtr       rrp;
  Uint1           shift;
  SeqMapTablePtr  smtp;
  Char            str [8];
  tRNAPtr         trna;
  ValNodePtr      vnp;

  rpp = (RnaPagePtr) GetObjectExtra (d);
  rrp = (RnaRefPtr) data;
  if (rpp != NULL) {
    SetEnumPopup (rpp->type, rna_type_alist,
                  (UIEnum) check_rna_type (rrp->type));
    SafeSetStatus (rpp->pseudo, rrp->pseudo);
    switch (rrp->type) {
      case 0 :
        SafeHide (rpp->trnaGrp);
        SafeHide (rpp->nameGrp);
        SafeHide (rpp->rrnaPrompt);
        break;
      case 3 :
        SafeSetTitle (rpp->name, "");
        SafeHide (rpp->nameGrp);
        SafeHide (rpp->rrnaPrompt);
        SafeSetValue (rpp->AAitem, 1);
        if (rrp->ext.choice == 2) {
          trna = rrp->ext.value.ptrvalue;
          if (trna != NULL) {
            aa = 0;
            if (trna->aatype == 2) {
              aa = trna->aa;
            } else {
              from = 0;
              switch (trna->aatype) {
                case 0 :
                  from = 0;
                  break;
                case 1 :
                  from = Seq_code_iupacaa;
                  break;
                case 2 :
                  from = Seq_code_ncbieaa;
                  break;
                case 3 :
                  from = Seq_code_ncbi8aa;
                  break;
                case 4 :
                  from = Seq_code_ncbistdaa;
                  break;
                default:
                  break;
              }
              smtp = SeqMapTableFind (Seq_code_ncbieaa, from);
              if (smtp != NULL) {
                aa = SeqMapTableConvert (smtp, trna->aa);
              }
            }
            if (aa > 0 && aa != 255) {
              if (aa <= 74) {
                shift = 0;
              } else if (aa > 79) {
                shift = 2;
              } else {
                shift = 1;
              } 
              i = aa - (64 + shift);
              SetValue (rpp->AAitem, (Int2) i + 1);
            }
            head = NULL;
            for (j = 0; j < 6; j++) {
              if (trna->codon [j] < 64) {
                if (CodonForIndex (trna->codon [j], Seq_code_iupacna, codon)) {
                  StringCpy (str, (CharPtr) codon);
                  str [3] = '\0';
                  vnp = ValNodeNew (head);
                  if (head == NULL) {
                    head = vnp;
                  }
                  if (vnp != NULL) {
                    vnp->data.ptrvalue = StringSave (str);
                  }
                }
              }
            }
            if (head != NULL) {
              PointerToDialog (rpp->codons, head);
              ValNodeFreeData (head);
            }
            PointerToDialog (rpp->anticodon, trna->anticodon);
          }
        }
        SafeShow (rpp->trnaGrp);
        break;
      case 4 :
        SafeHide (rpp->trnaGrp);
        if (rrp->ext.choice == 1 && rrp->ext.value.ptrvalue != NULL) {
          SafeSetTitle (rpp->name, (CharPtr) rrp->ext.value.ptrvalue);
        } else {
          SafeSetTitle (rpp->name, "");
        }
        SafeShow (rpp->rrnaPrompt);
        SafeShow (rpp->nameGrp);
        break;
      default :
        SafeHide (rpp->trnaGrp);
        if (rrp->ext.choice == 1 && rrp->ext.value.ptrvalue != NULL) {
          SafeSetTitle (rpp->name, (CharPtr) rrp->ext.value.ptrvalue);
        } else {
          SafeSetTitle (rpp->name, "");
        }
        SafeHide (rpp->rrnaPrompt);
        SafeShow (rpp->nameGrp);
        break;
    }
  } else {
    SetEnumPopup (rpp->type, rna_type_alist, (UIEnum) 0);
    SafeSetStatus (rpp->pseudo, FALSE);
    SafeSetTitle (rpp->name, "");
    SafeHide (rpp->nameGrp);
    SafeHide (rpp->rrnaPrompt);
    SafeSetValue (rpp->AAitem, 0);
    SafeHide (rpp->trnaGrp);
  }
}

static Pointer RnaPageToRnaRefPtr (DialoG d)

{
  Char        ch;
  Uint1       code;
  Uint1       codon [4];
  ValNodePtr  head;
  Int2        i;
  Int2        j;
  Int2        k;
  Int2        q;
  RnaPagePtr  rpp;
  RnaRefPtr   rrp;
  Uint1       shift;
  Char        str [8];
  tRNAPtr     trna;
  UIEnum      val;
  ValNodePtr  vnp;

  rrp = NULL;
  rpp = (RnaPagePtr) GetObjectExtra (d);
  if (rpp != NULL) {
    rrp = RnaRefNew ();
    if (rrp != NULL) {
      if (GetEnumPopup (rpp->type, rna_type_alist, &val)) {
        rrp->type = (Uint1) val;
      }
      rrp->pseudo = GetStatus (rpp->pseudo);
      switch (rrp->type) {
        case 3 :
          rrp->ext.choice = 2;
          trna = (tRNAPtr) MemNew (sizeof (tRNA));
          rrp->ext.value.ptrvalue = (Pointer) trna;
          if (trna != NULL) {
            i = GetValue (rpp->AAitem) - 1;
            if (i > 0) {
              trna->aatype = 2;
              if (i < 10) {
                shift = 0;
              } else if (i > 13) {
                shift = 2;
              } else {
                shift = 1;
              } 
              trna->aa = (Uint1) i + 64 + shift;
            }
            for (j = 0; j < 6; j++) {
              trna->codon [j] = 255;
            }
            head = (ValNodePtr) DialogToPointer (rpp->codons);
            if (head != NULL) {
              for (vnp = head, j = 0; vnp != NULL && j < 6; vnp = vnp->next) {
                str [0] = '\0';
                StringNCpy_0 (str, (CharPtr) vnp->data.ptrvalue, sizeof (str));
                if (str [0] != '\0') {
                  k = 0;
                  q = 0;
                  ch = str [k];
                  while (ch != '\0' && q < 3) {
                    ch = TO_UPPER (ch);
                    if (StringChr ("ACGTU", ch) != NULL) {
                      if (ch == 'U') {
                        ch = 'T';
                      }
                      codon [q] = (Uint1) ch;
                      q++;
                    }
                    k++;
                    ch = str [k];
                  }
                  codon [q] = 0;
                  if (q == 3) {
                    code = IndexForCodon (codon, Seq_code_iupacna);
                    if (code != INVALID_RESIDUE) {
                      trna->codon [j] = code;
                    }
                  }
                  j++;
                }
              }
              ValNodeFreeData (head);
            }
            trna->anticodon = DialogToPointer (rpp->anticodon);
          }
          break;
        default :
          if (! TextHasNoText (rpp->name)) {
            rrp->ext.choice = 1;
            rrp->ext.value.ptrvalue = SaveStringFromText (rpp->name);
          }
          break;
      }
    }
  }
  return (Pointer) rrp;
}

static void SetRnaImportExportItems (RnaFormPtr rfp)

{
  IteM        exportItm;
  IteM        importItm;
  RnaPagePtr  rpp;
  UIEnum      val;

  if (rfp != NULL) {
    importItm = FindFormMenuItem ((BaseFormPtr) rfp, VIB_MSG_IMPORT);
    exportItm = FindFormMenuItem ((BaseFormPtr) rfp, VIB_MSG_EXPORT);
    switch (rfp->currentPage) {
      case RNA_PAGE :
        rpp = (RnaPagePtr) GetObjectExtra (rfp->data);
        if (rpp != NULL && GetEnumPopup (rpp->type, rna_type_alist, &val) && val == 2) {
          SafeSetTitle (importItm, "Import cDNA FASTA...");
          SafeSetTitle (exportItm, "Export...");
          SafeEnable (importItm);
          SafeDisable (exportItm);
        } else {
          SafeSetTitle (importItm, "Import...");
          SafeSetTitle (exportItm, "Export...");
          SafeDisable (importItm);
          SafeDisable (exportItm);
        }
        break;
      case COMMON_PAGE :
        SafeSetTitle (importItm, "Import...");
        SafeSetTitle (exportItm, "Export...");
        SafeDisable (importItm);
        SafeDisable (exportItm);
        break;
      case LOCATION_PAGE :
        SafeSetTitle (importItm, "Import SeqLoc...");
        SafeSetTitle (exportItm, "Export SeqLoc...");
        SafeEnable (importItm);
        SafeEnable (exportItm);
        break;
      default :
        break;
    }
  }
}

static void ChangeRnaProc (PopuP p)

{
  RnaFormPtr  rfp;
  RnaPagePtr  rpp;
  UIEnum      val;

  rpp = (RnaPagePtr) GetObjectExtra (p);
  if (rpp != NULL) {
    rfp = rpp->rfp;
    if (GetEnumPopup (rpp->type, rna_type_alist, &val)) {
      switch (val) {
        case 0 :
          SafeHide (rpp->nameGrp);
          SafeHide (rpp->rrnaPrompt);
          SafeHide (rpp->trnaGrp);
          SafeHide (rfp->product);
          break;
        case 2 :
          SafeHide (rpp->trnaGrp);
          SafeHide (rpp->rrnaPrompt);
          SafeShow (rpp->nameGrp);
          SafeShow (rfp->product);
          break;
        case 3 :
          SafeHide (rpp->nameGrp);
          SafeHide (rpp->rrnaPrompt);
          SafeShow (rpp->trnaGrp);
          SafeHide (rfp->product);
          break;
        case 4 :
          SafeHide (rpp->trnaGrp);
          SafeShow (rpp->rrnaPrompt);
          SafeShow (rpp->nameGrp);
          SafeHide (rfp->product);
          break;
        default :
          SafeHide (rpp->trnaGrp);
          SafeHide (rpp->rrnaPrompt);
          SafeShow (rpp->nameGrp);
          SafeHide (rfp->product);
          break;
      }
    }
    SetRnaImportExportItems (rpp->rfp);
    Update ();
  }
}

static void PopulateAAPopup (PopuP AAitem)

{
  Char             ch;
  Uint1            first;
  Uint1            i;
  Char             item [77];
  Uint1            last;
  SeqCodeTablePtr  sctp;
  CharPtr          str;

  sctp = SeqCodeTableFind (Seq_code_ncbieaa);
  first = FirstResidueInCode (sctp);
  last = LastResidueInCode (sctp);
  PopupItem (AAitem, " ");
  for (i = 65; i <= last; i++) {
    if (i == 74 || i == 79) {
      continue;
    }
    ch = GetSymbolForResidue (sctp, i);
    str = (CharPtr) GetNameForResidue (sctp, i);
    sprintf (item, "%c    %s", ch, str);
    PopupItem (AAitem, item);
  } 
  SetValue (AAitem, 1); 
}

static void ChangetRNASubPage (VoidPtr data, Int2 newval, Int2 oldval)

{
  RnaPagePtr  rpp;

  rpp = (RnaPagePtr) data;
  if (rpp != NULL) {
    if (oldval >= 0 && oldval <= 2) {
      SafeHide (rpp->trnaPages [oldval]);
    }
    if (newval >= 0 && newval <= 2) {
      SafeShow (rpp->trnaPages [newval]);
    }
    SetRnaImportExportItems (rpp->rfp);
    Update ();
  }
}

static CharPtr trnaTabs [] = {
  "Amino Acid", "Codons", "Anticodon", NULL
};

static DialoG CreateRnaDialog (GrouP h, CharPtr title,
                               Uint2 subtype, SeqEntryPtr sep,
                               SeqFeatPtr sfp, RnaFormPtr rfp)

{
  GrouP       f;
  GrouP       g;
  Boolean     indexerVersion;
  GrouP       k;
  GrouP       m;
  GrouP       p;
  GrouP       q;
  RnaPagePtr  rpp;
  RnaRefPtr   rrp;
  GrouP       s;
  Boolean     showpseudo;
  PrompT      t;
  DialoG      tbs;

  p = HiddenGroup (h, 1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  rpp = (RnaPagePtr) MemNew (sizeof (RnaPage));
  if (rpp != NULL) {

    SetObjectExtra (p, rpp, StdCleanupExtraProc);
    rpp->dialog = (DialoG) p;
    rpp->todialog = RnaRefPtrToRnaPage;
    rpp->fromdialog = RnaPageToRnaRefPtr;
    rpp->testdialog = NULL;

    if (title != NULL && title [0] != '\0') {
      s = NormalGroup (p, 0, -2, title, systemFont, NULL);
    } else {
      s = HiddenGroup (p, 0, -2, NULL);
    }
    m = HiddenGroup (s, -1, 0, NULL);
    SetGroupSpacing (m, 10, 10);

    rpp->rfp = rfp;

    f = HiddenGroup (m, -2, 0, NULL);
    StaticPrompt (f, "Type", 0, popupMenuHeight, programFont, 'l');
    rpp->type = PopupList (f, TRUE, ChangeRnaProc);
    SetObjectExtra (rpp->type, rpp, NULL);
    InitEnumPopup (rpp->type, rna_type_alist, NULL);
    SetEnumPopup (rpp->type, rna_type_alist, (UIEnum) 0);

    showpseudo = FALSE;
    if (sfp != NULL && sfp->data.choice == SEQFEAT_RNA) {
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp != NULL && rrp->pseudo) {
        showpseudo = TRUE;
      }
    }
    indexerVersion = (Boolean) (GetAppProperty ("InternalNcbiSequin") != NULL);
    if (indexerVersion) {
      showpseudo = TRUE;
    }
    if (showpseudo) {
      rpp->pseudo = CheckBox (m, "PseudoRNA", NULL);
    }

    g = HiddenGroup (m, 0, 0, NULL);

    rpp->nameGrp = HiddenGroup (g, -1, 0, NULL);
    q = HiddenGroup (rpp->nameGrp, -2, 0, NULL);
    StaticPrompt (q, "Name", 0, dialogTextHeight, programFont, 'l');
    rpp->name = DialogText (q, "", 20, NULL);
    rpp->rrnaPrompt = StaticPrompt (rpp->nameGrp, "E.g., 16S ribosomal RNA", 0, 0, programFont, 'c');
    if ((sfp != NULL && sfp->product != NULL) || indexerVersion) {
      rfp->product = CreateProteinOrMRNAProductDialog (rpp->nameGrp, NULL, TRUE, sep, NULL);
    }
    AlignObjects (ALIGN_CENTER, (HANDLE) q, (HANDLE) rpp->rrnaPrompt, NULL);
    Hide (rpp->nameGrp);

    rpp->trnaGrp = HiddenGroup (g, -1, 0, NULL);
    SetGroupSpacing (rpp->trnaGrp, 10, 10);
    tbs = CreateFolderTabs (rpp->trnaGrp, trnaTabs, 0, 0, 0,
                            PROGRAM_FOLDER_TAB,
                            ChangetRNASubPage, (Pointer) rpp);
    k = HiddenGroup (rpp->trnaGrp, 0, 0, NULL);
    rpp->trnaPages [0] = HiddenGroup (k, -2, 0, NULL);
    StaticPrompt (rpp->trnaPages [0], "Amino Acid", 0, popupMenuHeight, programFont, 'l');
    rpp->AAitem = PopupList (rpp->trnaPages [0], TRUE, NULL);
    PopulateAAPopup (rpp->AAitem);
    rpp->trnaPages [1] = HiddenGroup (k, -1, 0, NULL);
    t = StaticPrompt (rpp->trnaPages [1], "Recognized Codons", 0, 0, programFont, 'c');
    rpp->codons = CreateVisibleStringDialog (rpp->trnaPages [1], 3, -1, 4);
    AlignObjects (ALIGN_CENTER, (HANDLE) t, (HANDLE) rpp->codons, NULL);
    Hide (rpp->trnaPages [1]);
    rpp->trnaPages [2] = HiddenGroup (k, -1, 0, NULL);
    rpp->anticodon = CreateIntervalEditorDialogEx (rpp->trnaPages [2], NULL,
                                                   3, 2, sep, TRUE, FALSE,
                                                   FALSE, FALSE, FALSE, NULL, NULL);
    Hide (rpp->trnaPages [2]);
    Hide (rpp->trnaGrp);

    AlignObjects (ALIGN_CENTER, (HANDLE) tbs,
                  (HANDLE) rpp->trnaPages [0], (HANDLE) rpp->trnaPages [1],
                  (HANDLE) rpp->trnaPages [2], NULL);
    AlignObjects (ALIGN_CENTER, (HANDLE) f, (HANDLE) rpp->nameGrp,
                  (HANDLE) rpp->trnaGrp,  (HANDLE) rpp->pseudo,NULL);

    switch (subtype) {
      case 0 :
        SafeHide (rpp->nameGrp);
        SafeHide (rpp->rrnaPrompt);
        SafeHide (rpp->trnaGrp);
        SafeHide (rfp->product);
        break;
      case FEATDEF_mRNA :
        SafeHide (rpp->trnaGrp);
        SafeHide (rpp->rrnaPrompt);
        SafeShow (rpp->nameGrp);
        SafeShow (rfp->product);
        break;
      case FEATDEF_tRNA :
        SafeHide (rpp->nameGrp);
        SafeHide (rpp->rrnaPrompt);
        SafeShow (rpp->trnaGrp);
        SafeHide (rfp->product);
        break;
      case FEATDEF_rRNA :
        SafeHide (rpp->trnaGrp);
        SafeShow (rpp->rrnaPrompt);
        SafeShow (rpp->nameGrp);
        SafeHide (rfp->product);
        break;
      default :
        SafeHide (rpp->trnaGrp);
        SafeHide (rpp->rrnaPrompt);
        SafeShow (rpp->nameGrp);
        SafeHide (rfp->product);
        break;
    }
  }

  return (DialoG) p;
}

static void ChangeRnaPage (VoidPtr data, Int2 newval, Int2 oldval)

{
  RnaFormPtr    rfp;

  rfp = (RnaFormPtr) data;
  if (rfp != NULL) {
    rfp->currentPage = newval;
    SafeHide (rfp->pages [oldval]);
    SafeShow (rfp->pages [newval]);
    switch (newval) {
      case RNA_PAGE :
        SendMessageToDialog (rfp->data, VIB_MSG_ENTER);
        break;
      case COMMON_PAGE :
        break;
      case LOCATION_PAGE :
        SendMessageToDialog (rfp->location, VIB_MSG_ENTER);
        break;
      default :
        break;
    }
    SetRnaImportExportItems (rfp);
    Update ();
  }
}

static Boolean ImportRnaForm (ForM f, CharPtr filename)

{
  BioseqPtr    bsp;
  SeqLocPtr    currslp;
  CharPtr      extension;
  FILE         *fp;
  BioseqPtr    mbsp;
  BioseqPtr    nbsp;
  SeqEntryPtr  nsep;
  SeqLocPtr    oldslp;
  Char         path [PATH_MAX];
  RnaFormPtr   rfp;
  SeqIdPtr     sip;
  SeqEntryPtr  sep;
  SeqLocPtr    slp;

  path [0] = '\0';
  StringNCpy_0 (path, filename, sizeof (path));
  rfp = (RnaFormPtr) GetObjectExtra (f);
  if (rfp != NULL) {
    switch (rfp->currentPage) {
      case RNA_PAGE :
        extension = GetAppProperty ("FastaNucExtension");
        if (path [0] != '\0' || GetInputFileName (path, sizeof (path), extension, "TEXT")) {
          fp = FileOpen (path, "r");
          if (fp != NULL) {
            sep = FastaToSeqEntry (fp, TRUE);
            if (sep != NULL && sep->choice == 1 && sep->data.ptrvalue != NULL) {
              mbsp = (BioseqPtr) sep->data.ptrvalue;
              nsep = GetBestTopParentForItemID (rfp->input_entityID,
                                                rfp->input_itemID,
                                                rfp->input_itemtype);
              if (nsep != NULL) {
                nsep = FindNucSeqEntry (nsep);
              }
              if (nsep != NULL && nsep->data.ptrvalue != NULL) {
                nbsp = (BioseqPtr) nsep->data.ptrvalue;
                slp = AlignmRNA2genomic (nbsp, mbsp);
                if (slp != NULL) {
                  currslp = DialogToPointer (rfp->location);
                  if (currslp != NULL) {
                    if (Message (MSG_OKC, "Replace existing location?") == ANS_CANCEL) {
                      slp = SeqLocFree (slp);
                    }
                    Select (rfp->form);
                    Update ();
                  }
                  SeqLocFree (currslp);
                  if (slp != NULL) {
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
                    PointerToDialog (rfp->location, slp);
                    SeqLocFree (slp);
                  }
                }
              }
            }
            SeqEntryFree (sep);
          }
          FileClose (fp);
        }
        break;
      case LOCATION_PAGE :
        return ImportDialog (rfp->location, filename);
      default :
        break;
    }
  }
  return FALSE;
}

static Boolean ExportRnaForm (ForM f, CharPtr filename)

{
  RnaFormPtr  rfp;

  rfp = (RnaFormPtr) GetObjectExtra (f);
  if (rfp != NULL) {
    switch (rfp->currentPage) {
      case LOCATION_PAGE :
        return ExportDialog (rfp->location, filename);
      default :
        break;
    }
  }
  return FALSE;
}

static CharPtr  rnaFormTabs [] = {
  NULL, "Properties", "Location", NULL
};

static void RnaFormMessage (ForM f, Int2 mssg)

{
  RnaFormPtr  rfp;

  rfp = (RnaFormPtr) GetObjectExtra (f);
  if (rfp != NULL) {
    switch (mssg) {
      case VIB_MSG_INIT :
        StdInitFeatFormProc (f);
        break;
      case VIB_MSG_IMPORT :
        ImportRnaForm (f, NULL);
        break;
      case VIB_MSG_EXPORT :
        ExportRnaForm (f, NULL);
        break;
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
        if (rfp->currentPage == LOCATION_PAGE) {
          PointerToDialog (rfp->location, NULL);
        } else {
          StdDeleteTextProc (NULL);
        }
        break;
      default :
        if (rfp->appmessage != NULL) {
          rfp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static void RnaFormActivate (WindoW w)

{
  RnaFormPtr  rfp;

  rfp = (RnaFormPtr) GetObjectExtra (w);
  if (rfp != NULL) {
    if (rfp->activate != NULL) {
      rfp->activate (w);
    }
    SetRnaImportExportItems (rfp);
  }
}

static void RnaRefPtrToForm (ForM f, Pointer data)

{
  FeatureFormPtr  ffp;
  SeqFeatPtr      sfp;
  Int4            val;

  ffp = (FeatureFormPtr) GetObjectExtra (f);
  if (ffp != NULL) {
    sfp = (SeqFeatPtr) data;
    if (sfp != NULL) {
      switch (sfp->data.choice) {
        case SEQFEAT_BOND :
        case SEQFEAT_SITE :
        case SEQFEAT_PSEC_STR :
          val = (Int4) sfp->data.value.intvalue;
          PointerToDialog (ffp->data, (Pointer) &(val));
          break;
        case SEQFEAT_COMMENT:
          break;
        default :
          PointerToDialog (ffp->data, sfp->data.value.ptrvalue);
          break;
      }
      SeqFeatPtrToCommon (ffp, sfp);
      PointerToDialog (ffp->location, sfp->location);
      PointerToDialog (ffp->product, sfp->product);
    }
  }
}

extern ForM CreateRnaForm (Int2 left, Int2 top, CharPtr title,
                           SeqFeatPtr sfp, SeqEntryPtr sep,
                           Uint2 subtype, FormActnFunc actproc)

{
  ButtoN             b;
  GrouP              c;
  GrouP              g;
  GrouP              h;
  RnaFormPtr         rfp;
  GrouP              s;
  StdEditorProcsPtr  sepp;
  WindoW             w;

  w = NULL;
  rfp = (RnaFormPtr) MemNew (sizeof (RnaForm));
  if (rfp != NULL) {
    w = FixedWindow (left, top, -10, -10, title, StdCloseWindowProc);
    SetObjectExtra (w, rfp, StdFeatFormCleanupProc);
    rfp->form = (ForM) w;
    rfp->actproc = actproc;
    rfp->toform = RnaRefPtrToForm;
    rfp->fromform = NULL;
    rfp->formmessage = RnaFormMessage;
    rfp->testform = NULL;
    rfp->importform = ImportRnaForm;
    rfp->exportform = ExportRnaForm;

#ifndef WIN_MAC
    CreateStdEditorFormMenus (w);
#endif

    rfp->activate = NULL;
    sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
    if (sepp != NULL) {
      rfp->activate = sepp->activateForm;
      rfp->appmessage = sepp->handleMessages;
    }
    SetActivate (w, RnaFormActivate);

    g = HiddenGroup (w, -1, 0, NULL);
    SetGroupSpacing (g, 3, 10);

    rfp->sep = sep;
    rnaFormTabs [0] = NULL;
    if (title != NULL && *title != '\0') {
      rnaFormTabs [0] = title;
    } else {
      rnaFormTabs [0] = "RNA";
    }
    rfp->foldertabs = CreateFolderTabs (g, rnaFormTabs, RNA_PAGE,
                                        0, 0, SYSTEM_FOLDER_TAB,
                                        ChangeRnaPage, (Pointer) rfp);
    rfp->currentPage = RNA_PAGE;

    h = HiddenGroup (g, 0, 0, NULL);

    s = HiddenGroup (h, -1, 0, NULL);
    SetGroupSpacing (s, 3, 10);
    rfp->data = CreateRnaDialog (s, NULL, subtype, sep, sfp, rfp);
    rfp->pages [RNA_PAGE] = s;
    Hide (rfp->pages [RNA_PAGE]);
    rnaFormTabs [0] = NULL;

    s = HiddenGroup (h, -1, 0, NULL);
    CreateCommonFeatureGroup (s, (FeatureFormPtr) rfp, sfp, TRUE, TRUE);
    rfp->pages [COMMON_PAGE] = s;
    Hide (rfp->pages [COMMON_PAGE]);

    s = HiddenGroup (h, -1, 0, NULL);
    rfp->location = CreateIntervalEditorDialogEx (s, NULL, 4, 2, sep, TRUE, FALSE,
                                                  TRUE, TRUE, FALSE,
                                                  (FeatureFormPtr) rfp,
                                                  StdFeatIntEdPartialCallback);
    rfp->pages [LOCATION_PAGE] = s;
    Hide (rfp->pages [LOCATION_PAGE]);

    AlignObjects (ALIGN_CENTER, (HANDLE) rfp->pages [RNA_PAGE],
                  (HANDLE) rfp->pages [COMMON_PAGE], (HANDLE) rfp->pages [LOCATION_PAGE],
                  NULL);
    AlignObjects (ALIGN_CENTER, (HANDLE) rfp->foldertabs, (HANDLE) h, NULL);

    c = HiddenGroup (w, 2, 0, NULL);
    b = PushButton (c, "Accept", StdFeatFormAcceptButtonProc);
    SetObjectExtra (b, rfp, NULL);
    PushButton (c, "Cancel", StdCancelButtonProc);
    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
    RealizeWindow (w);

    SendMessageToDialog (rfp->data, VIB_MSG_INIT);
    SendMessageToDialog (rfp->location, VIB_MSG_INIT);
    SendMessageToDialog (rfp->product, VIB_MSG_INIT);
    Show (rfp->pages [rfp->currentPage]);
    SendMessageToDialog (rfp->data, VIB_MSG_ENTER);
    Update ();
  }
  return (ForM) w;
}

extern Int2 LIBCALLBACK RnaGenFunc (Pointer data)

{
  FeatDefPtr        curr;
  HelpMessageFunc   helpfunc;
  Uint1             key;
  CharPtr           label = NULL;
  ObjMgrPtr         omp;
  OMProcControlPtr  ompcp;
  ObjMgrTypePtr     omtp;
  OMUserDataPtr     omudp;
  ObjMgrProcPtr     proc;
  RnaFormPtr        rfp;
  RnaPagePtr        rpp;
  RnaRefPtr         rrp;
  SeqEntryPtr       sep;
  SeqFeatPtr        sfp;
  Uint2             subtype;
  Char              title [64];
  WindoW            w;

  ompcp = (OMProcControlPtr) data;
  sfp = NULL;
  sep = NULL;
  subtype = FEATDEF_otherRNA;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  proc = ompcp->proc;
  switch (ompcp->input_itemtype) {
    case OBJ_SEQFEAT :
      sfp = (SeqFeatPtr) ompcp->input_data;
      if (sfp != NULL && sfp->data.choice != SEQFEAT_RNA) {
        return OM_MSG_RET_ERROR;
      }
      break;
    case OBJ_BIOSEQ :
      break;
    case OBJ_BIOSEQSET :
      break;
    case 0 :
      break;
    default :
      return OM_MSG_RET_ERROR;
  }
  omudp = ItemAlreadyHasEditor (ompcp->input_entityID, ompcp->input_itemID,
                                ompcp->input_itemtype, ompcp->proc->procid);
  if (omudp != NULL) {
    rfp = (RnaFormPtr) omudp->userdata.ptrvalue;
    if (rfp != NULL) {
      Select (rfp->form);
    }
    return OM_MSG_RET_DONE;
  }
  StringCpy (title, "Unknown RNA");
  if (sfp != NULL && sfp->data.value.ptrvalue != NULL) {
    rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
    StringNCpy_0 (title,
                GetEnumName ((UIEnum) rrp->type, rna_type_alist),
                sizeof (title));
    if (StringHasNoText (title)) {
      StringCpy (title, "Unknown RNA");
    } else if (proc->subinputtype > 0) {
      subtype = proc->subinputtype;
    }
  } else if (proc->subinputtype > 0) {
    subtype = proc->subinputtype;
    curr = FeatDefFindNext (NULL, &key, &label, FEATDEF_ANY, TRUE);
    while (curr != NULL) {
      if (key != FEATDEF_BAD && curr->seqfeat_key == SEQFEAT_RNA) {
        if (subtype == curr->featdef_key) {
          StringNCpy_0 (title, curr->typelabel, sizeof (title));
          break;
        }
      }
      curr = FeatDefFindNext (curr, &key, &label, FEATDEF_ANY, TRUE);
    }
  }
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  w = (WindoW) CreateRnaForm (-50, -33, title, sfp, sep,
                              subtype, StdFeatFormActnProc);
  rfp = (RnaFormPtr) GetObjectExtra (w);
  if (rfp != NULL) {
    rfp->input_entityID = ompcp->input_entityID;
    rfp->input_itemID = ompcp->input_itemID;
    rfp->input_itemtype = ompcp->input_itemtype;
    rfp->this_itemtype = OBJ_SEQFEAT;
    rfp->this_subtype = subtype;
    rfp->procid = ompcp->proc->procid;
    rfp->proctype = ompcp->proc->proctype;
    rfp->userkey = OMGetNextUserKey ();
    omudp = ObjMgrAddUserData (ompcp->input_entityID, ompcp->proc->procid,
	                           OMPROC_EDIT, rfp->userkey);
    if (omudp != NULL) {
      omudp->userdata.ptrvalue = (Pointer) rfp;
      omudp->messagefunc = StdVibrantEditorMsgFunc;
    }
    if (sfp != NULL) {
      omp = ObjMgrGet ();
      if (omp != NULL) {
        omtp = ObjMgrTypeFind (omp, OBJ_SEQFEAT, NULL, NULL);
        if (omtp != NULL && omtp->subtypefunc != NULL) {
          rfp->this_subtype = (*(omtp->subtypefunc)) (sfp);
        }
      }
    }
    SendMessageToForm (rfp->form, VIB_MSG_INIT);
    if (sfp != NULL) {
      PointerToForm (rfp->form, (Pointer) sfp);
      SetClosestParentIfDuplicating ((BaseFormPtr) rfp);
    } else {
      rpp = (RnaPagePtr) GetObjectExtra (rfp->data);
      if (rpp != NULL) {
        if (subtype >= FEATDEF_preRNA && subtype <= FEATDEF_scRNA) {
          SetEnumPopup (rpp->type, rna_type_alist,
                        (UIEnum) check_rna_type (subtype - FEATDEF_preRNA + 1));
        } else if (subtype == FEATDEF_otherRNA) {
          SetEnumPopup (rpp->type, rna_type_alist,
                        (UIEnum) check_rna_type (255));
        } else {
          SetEnumPopup (rpp->type, rna_type_alist, (UIEnum) 0);
        }
      }
      SetNewFeatureDefaultInterval ((FeatureFormPtr) rfp);
    }
  }
  Show (w);
  Select (w);
  helpfunc = (HelpMessageFunc) GetAppProperty ("HelpMessageProc");
  if (helpfunc != NULL) {
    helpfunc ("Features", title);
  }
  return OM_MSG_RET_DONE;
}

