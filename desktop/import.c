/*   import.c
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
* File Name:  import.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   6/18/95
*
* $Revision: 6.12 $
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

#include <import.h>
#include <objfdef.h>
#include <gbfeat.h>
#include <gbftdef.h>
#include <gather.h>
#include <subutil.h>    /* TOPOLOGY_xxx definitions */

#define IMPORT_PAGE      0
#define ENUM_PAGE        0
#define REGION_PAGE      0
#define COMMON_PAGE      1
#define LOCATION_PAGE    2

#define NUM_PAGES  8

typedef struct imprtform {
  FEATURE_FORM_BLOCK
  SeqEntryPtr   sep;
  GrouP         pages [NUM_PAGES];
  DialoG        foldertabs;
  Int2          currentPage;
} ImprtForm, PNTR ImprtFormPtr;

typedef struct importpage {
  DIALOG_MESSAGE_BLOCK
  ImprtFormPtr       ifp;
  Handle             key;
  TexT               loc;
  EnumFieldAssoc     PNTR alist;
} ImportPage, PNTR ImportPagePtr;

extern EnumFieldAssocPtr import_featdef_alist (Boolean notJustImpFeats, Boolean allowPeptideFeats, Boolean allowRnaImpFeats)

{
  EnumFieldAssocPtr  alist;
  EnumFieldAssocPtr  ap;
  FeatDefPtr         curr;
  Int2               fdn;
  Int2               i;
  Uint1              key;
  CharPtr            label = NULL;
  CharPtr            str;
  Uint2              subtype;

  fdn = FeatDefNum ();
  alist = MemNew (sizeof (EnumFieldAssoc) * (fdn + 2));
  if (alist == NULL) {
    Message (MSG_ERROR, "in import_featdef_alist: no room");
  } else {
    i = 0;
    ap = alist;
    ap->name = StringSave ("     ");
    ap->value = (UIEnum) 0;
    ap++;
    i++;
    curr = FeatDefFindNext (NULL, &key, &label, FEATDEF_ANY, TRUE);
    while (curr != NULL) {
      if (key != FEATDEF_BAD) {
        if (notJustImpFeats || curr->seqfeat_key == SEQFEAT_IMP) {
          subtype = curr->featdef_key;
          if (subtype != FEATDEF_Imp_CDS &&
               subtype != FEATDEF_source &&
              subtype != FEATDEF_site_ref) {
            if (allowPeptideFeats ||
                (subtype != FEATDEF_mat_peptide &&
                 subtype != FEATDEF_sig_peptide &&
                 subtype != FEATDEF_transit_peptide)) {
              if (allowRnaImpFeats ||
                (subtype != FEATDEF_misc_RNA &&
                 subtype != FEATDEF_precursor_RNA)) {
                str = StringSave (curr->typelabel);
                if (*str == '-') {
                  *str = '~';
                }
                ap->name = str;
                ap->value = key;
                ap++;
                i++;
              }
            }
          }
        }
      }
      curr = FeatDefFindNext (curr, &key, &label, FEATDEF_ANY, TRUE);
    }
    ap->name = NULL;
  }
  return alist;
}

static void ImpFeatPtrToImportPage (DialoG d, Pointer data)

{
  EnumFieldAssocPtr  ap;
  Int2               i;
  ImportPagePtr      ipp;
  ImpFeatPtr         ifp;

  ipp = (ImportPagePtr) GetObjectExtra (d);
  ifp = (ImpFeatPtr) data;
  if (ipp != NULL) {
    if (ifp != NULL) {
      for (i = 1, ap = ipp->alist; ap->name != NULL; i++, ap++) {
        if (StringCmp (ap->name, ifp->key) == 0) {
          SetValue (ipp->key, i);
          SetTitle (ipp->loc, ifp->loc);
          return;
        }
      }
    }
    SetValue (ipp->key, 1);
    SetTitle (ipp->loc, "");
  }
}

static Pointer ImportPageToImpFeatPtr (DialoG d)

{
  EnumFieldAssocPtr  ap;
  Int2               i;
  ImportPagePtr      ipp;
  ImpFeatPtr         ifp;
  UIEnum             val;

  ifp = NULL;
  ipp = (ImportPagePtr) GetObjectExtra (d);
  if (ipp != NULL) {
    ifp = ImpFeatNew ();
    if (ifp != NULL) {
      ifp->loc = SaveStringFromText (ipp->loc);
      val = GetValue (ipp->key);
      if (val > 1) {
        for (i = 1, ap = ipp->alist; ap->name != NULL; i++, ap++) {
          if (i == val) {
            ifp->key = StringSave (ap->name);
            return (Pointer) ifp;
          }
        }
      } else {
        ifp->key = StringSave ("misc_feature");
      }
    }
  }
  return (Pointer) ifp;
}

static void CleanupImportPage (GraphiC g, VoidPtr data)

{
  ImportPagePtr  ipp;
  Int2           j;

  ipp = (ImportPagePtr) data;
  if (ipp != NULL) {
    if (ipp->alist != NULL) {
      for (j = 0; ipp->alist [j].name != NULL; j++) {
        MemFree (ipp->alist [j].name);
      }
    }
    MemFree (ipp->alist);
  }
  MemFree (data);
}

static void ChangeKey (Handle obj)

{
  Char             ch;
  Int2             expev;
  Int2             geneval;
  HelpMessageFunc  helpfunc;
  ImprtFormPtr     ifp;
  ImpFeatPtr       imp;
  ImprtFormPtr     newifp;
  ObjMgrPtr        omp;
  ObjMgrTypePtr    omtp;
  CharPtr          ptr;
  SeqEntryPtr      sep;
  SeqFeatPtr       sfp;
  Char             title [128];
  WindoW           w;

  ifp = (ImprtFormPtr) GetObjectExtra (obj);
  if (ifp != NULL) {
    sfp = SeqFeatNew ();
    if (sfp != NULL) {
      sfp->data.choice = FindFeatFromFeatDefType (ifp->this_subtype);
      sfp->data.value.ptrvalue = DialogToPointer (ifp->data);
      sfp->comment = SaveStringFromText (ifp->comment);
      ptr = sfp->comment;
      if (ptr != NULL) {
        ch = *ptr;
        while (ch != '\0') {
          if (ch < ' ' || ch > '~') {
            *ptr = '~';
          }
          ptr++;
          ch = *ptr;
        }
      }
      expev = GetValue (ifp->evidence);
      if (expev > 0 && expev <= 3) {
        sfp->exp_ev = expev - 1;
      } else {
        sfp->exp_ev = 0;
      }
      sfp->partial = GetStatus (ifp->partial);
      sfp->excpt = GetStatus (ifp->exception);
      sfp->title = NULL;
      sfp->product = DialogToPointer (ifp->product);
      sfp->location = DialogToPointer (ifp->location);
      sfp->cit = DialogToPointer (ifp->featcits);
      sfp->dbxref = DialogToPointer (ifp->dbxrefs);
      sfp->qual = DialogToPointer (ifp->gbquals);
      geneval = GetValue (ifp->gene);
      sep = GetTopSeqEntryForEntityID (ifp->input_entityID);
      StringCpy (title, "Feature");
      if (sfp != NULL && sfp->data.value.ptrvalue != NULL) {
        imp = (ImpFeatPtr) sfp->data.value.ptrvalue;
        StringNCpy_0 (title, imp->key, sizeof (title));
        if (StringHasNoText (title)) {
          StringCpy (title, "Feature");
        }
      }
      w = (WindoW) CreateImportForm (-50, -33, title, sfp, sep,
                                     StdFeatFormActnProc);
      newifp = (ImprtFormPtr) GetObjectExtra (w);
      if (newifp != NULL) {
        newifp->input_entityID = ifp->input_entityID;
        newifp->input_itemID = ifp->input_itemID;
        newifp->input_itemtype = ifp->input_itemtype;
        newifp->this_itemtype = ifp->this_itemtype;
        newifp->this_subtype = ifp->this_subtype;
        if (sfp != NULL) {
          omp = ObjMgrGet ();
          if (omp != NULL) {
            omtp = ObjMgrTypeFind (omp, OBJ_SEQFEAT, NULL, NULL);
            if (omtp != NULL && omtp->subtypefunc != NULL) {
              newifp->this_subtype = (*(omtp->subtypefunc)) (sfp);
            }
          }
        }
        SendMessageToForm (newifp->form, VIB_MSG_INIT);
        SetValue (newifp->gene, geneval);
        if (sfp != NULL) {
          PointerToForm (newifp->form, (Pointer) sfp);
        }
      }
      Remove (ifp->form);
      Show (w);
      Select (w);
      helpfunc = (HelpMessageFunc) GetAppProperty ("HelpMessageProc");
      if (helpfunc != NULL) {
        helpfunc ("Features", title);
      }
    }
    SeqFeatFree (sfp);
    Update ();
  }
}

static DialoG CreateImportDialog (GrouP h, CharPtr title, ImprtFormPtr ifp, Boolean allowPeptideFeats)

{
  EnumFieldAssocPtr  ap;
  GrouP              f;
  ImportPagePtr      ipp;
  GrouP              m;
  GrouP              p;
  PrompT             ppt;
  GrouP              s;
  CharPtr            str;
  Boolean            usePopup;
  PrompT             x;

  p = HiddenGroup (h, 1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  ipp = (ImportPagePtr) MemNew (sizeof (ImportPage));
  if (ipp != NULL) {

    SetObjectExtra (p, ipp, CleanupImportPage);
    ipp->dialog = (DialoG) p;
    ipp->todialog = ImpFeatPtrToImportPage;
    ipp->fromdialog = ImportPageToImpFeatPtr;
    ipp->testdialog = NULL;

    if (title != NULL && title [0] != '\0') {
      s = NormalGroup (p, 0, -2, title, systemFont, NULL);
    } else {
      s = HiddenGroup (p, 0, -2, NULL);
    }
    m = HiddenGroup (s, -1, 0, NULL);
    /*
    SetGroupSpacing (m, 10, 10);
    */

    ipp->ifp = ifp;
    ipp->alist = import_featdef_alist (FALSE, allowPeptideFeats, FALSE);

    ppt = StaticPrompt (m, "Changing feature key will recreate the window.",
                        0, 0, programFont, 'c');
    f = HiddenGroup (m, -2, 0, NULL);

#ifdef WIN_MOTIF
    usePopup = FALSE;
#else
    usePopup = TRUE;
#endif

    x = NULL;
#ifndef WIN_MOTIF
    x = StaticPrompt (f, "Key", 0, 0, programFont, 'l');
#endif
    if (usePopup) {
      ipp->key = (Handle) PopupList (f, TRUE, (PupActnProc) ChangeKey);
      SetObjectExtra (ipp->key, ifp, NULL);
      InitEnumPopup ((PopuP) ipp->key, ipp->alist, NULL);
      SetEnumPopup ((PopuP) ipp->key, ipp->alist, (UIEnum) 0);
    } else {
      ipp->key = (Handle) SingleList (f, 12, 3, (LstActnProc) ChangeKey);
      SetObjectExtra (ipp->key, ifp, NULL);
      for (ap = ipp->alist; ap->name != NULL; ap++) {
        ListItem ((LisT) ipp->key, ap->name);
      }
      SetValue (ipp->key, 0);
    }
    if (Nlm_GetAppProperty ("SequinUseImpFeatLocField") != NULL) {
      StaticPrompt (f, "Loc", 0, dialogTextHeight, programFont, 'l');
      ipp->loc = DialogText (f, "", 15, NULL);
    }
    AlignObjects (ALIGN_VERTICAL, (HANDLE) ipp->key, (HANDLE) x, NULL);
    for (ap = ipp->alist; ap->name != NULL; ap++) {
      str = ap->name;
      if (*str == '~') {
        *str = '-';
      }
    }
    AlignObjects (ALIGN_CENTER, (HANDLE) ppt, (HANDLE) f, NULL);
  }

  return (DialoG) p;
}

static void SetImportImportExportItems (ImprtFormPtr ifp)

{
  IteM  exportItm;
  IteM  importItm;

  if (ifp != NULL) {
    importItm = FindFormMenuItem ((BaseFormPtr) ifp, VIB_MSG_IMPORT);
    exportItm = FindFormMenuItem ((BaseFormPtr) ifp, VIB_MSG_EXPORT);
    switch (ifp->currentPage) {
      case IMPORT_PAGE :
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

static void ChangeImportPage (VoidPtr data, Int2 newval, Int2 oldval)

{
  ImprtFormPtr  ifp;

  ifp = (ImprtFormPtr) data;
  if (ifp != NULL) {
    ifp->currentPage = newval;
    SafeHide (ifp->pages [oldval]);
    SafeShow (ifp->pages [newval]);
    switch (newval) {
      case IMPORT_PAGE :
        SendMessageToDialog (ifp->data, VIB_MSG_ENTER);
        break;
      case COMMON_PAGE :
        break;
      case LOCATION_PAGE :
        SendMessageToDialog (ifp->location, VIB_MSG_ENTER);
        break;
      default :
        break;
    }
    SetImportImportExportItems (ifp);
    Update ();
  }
}

static Boolean ImportImportForm (ForM f, CharPtr filename)

{
  ImprtFormPtr  ifp;

  ifp = (ImprtFormPtr) GetObjectExtra (f);
  if (ifp != NULL) {
    switch (ifp->currentPage) {
      case LOCATION_PAGE :
        return ImportDialog (ifp->location, filename);
      default :
        break;
    }
  }
  return FALSE;
}

static Boolean ExportImportForm (ForM f, CharPtr filename)

{
  ImprtFormPtr  ifp;

  ifp = (ImprtFormPtr) GetObjectExtra (f);
  if (ifp != NULL) {
    switch (ifp->currentPage) {
      case LOCATION_PAGE :
        return ExportDialog (ifp->location, filename);
      default :
        break;
    }
  }
  return FALSE;
}

static CharPtr  importFormTabs [] = {
  NULL, "Properties", "Location", NULL
};

extern DialoG CreateQualsDialog (GrouP h, Uint2 rows, Int2 spacing,
                                 Int2 width1, Int2 width2);

static void ImportFormMessage (ForM f, Int2 mssg)

{
  ImprtFormPtr  ifp;

  ifp = (ImprtFormPtr) GetObjectExtra (f);
  if (ifp != NULL) {
    switch (mssg) {
      case VIB_MSG_INIT :
        StdInitFeatFormProc (f);
        break;
      case VIB_MSG_IMPORT :
        ImportImportForm (f, NULL);
        break;
      case VIB_MSG_EXPORT :
        ExportImportForm (f, NULL);
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
        if (ifp->currentPage == LOCATION_PAGE) {
          PointerToDialog (ifp->location, NULL);
        } else {
          StdDeleteTextProc (NULL);
        }
        break;
      default :
        if (ifp->appmessage != NULL) {
          ifp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static void ImportFormActivate (WindoW w)

{
  ImprtFormPtr  ifp;

  ifp = (ImprtFormPtr) GetObjectExtra (w);
  if (ifp != NULL) {
    if (ifp->activate != NULL) {
      ifp->activate (w);
    }
    SetImportImportExportItems (ifp);
  }
}

extern ForM CreateImportForm (Int2 left, Int2 top, CharPtr title,
                              SeqFeatPtr sfp, SeqEntryPtr sep,
                              FormActnFunc actproc)

{
  Boolean            allowPeptideFeats;
  Boolean            allowProductGBQual;
  ButtoN             b;
  GrouP              c;
  GrouP              g;
  GrouP              h;
  ImprtFormPtr       ifp;
  ImpFeatPtr         imp;
  GrouP              s;
  StdEditorProcsPtr  sepp;
  WindoW             w;
  GrouP              x;
  GrouP              z;

  w = NULL;
  ifp = (ImprtFormPtr) MemNew (sizeof (ImprtForm));
  if (ifp != NULL) {
    w = FixedWindow (left, top, -10, -10, title, StdCloseWindowProc);
    SetObjectExtra (w, ifp, StdFeatFormCleanupProc);
    ifp->form = (ForM) w;
    ifp->actproc = actproc;
    ifp->toform = StdSeqFeatPtrToFeatFormProc;
    ifp->fromform = NULL;
    ifp->formmessage = ImportFormMessage;
    ifp->testform = NULL;
    ifp->importform = ImportImportForm;
    ifp->exportform = ExportImportForm;

#ifndef WIN_MAC
    CreateStdEditorFormMenus (w);
#endif

    ifp->activate = NULL;
    sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
    if (sepp != NULL) {
      ifp->activate = sepp->activateForm;
      ifp->appmessage = sepp->handleMessages;
    }
    SetActivate (w, ImportFormActivate);

    g = HiddenGroup (w, -1, 0, NULL);
    SetGroupSpacing (g, 3, 10);

    ifp->sep = sep;
    importFormTabs [0] = NULL;
    if (title != NULL && *title != '\0') {
      importFormTabs [0] = title;
    } else {
      importFormTabs [0] = "misc_feature";
    }
    ifp->foldertabs = CreateFolderTabs (g, importFormTabs, IMPORT_PAGE,
                                        0, 0, SYSTEM_FOLDER_TAB,
                                        ChangeImportPage, (Pointer) ifp);
    ifp->currentPage = IMPORT_PAGE;

    h = HiddenGroup (g, 0, 0, NULL);

    s = HiddenGroup (h, -1, 0, NULL);
    SetGroupSpacing (s, 3, 10);
    allowPeptideFeats = FALSE;
    allowProductGBQual = FALSE;
    if (sfp != NULL && sfp->data.choice == SEQFEAT_IMP) {
      imp = (ImpFeatPtr) sfp->data.value.ptrvalue;
      if (imp != NULL) {
        if (StringICmp (imp->key, "mat_peptide") == 0 ||
            StringICmp (imp->key, "sig_peptide") == 0 ||
            StringICmp (imp->key, "transit_peptide") == 0) {
          allowPeptideFeats = TRUE;
          allowProductGBQual = TRUE;
        }
      }
    }
    ifp->data = CreateImportDialog (s, NULL, ifp, allowPeptideFeats);
    x = HiddenGroup (s, -1, 0, NULL);
    if (StringICmp (importFormTabs [0], "source") == 0) {
      z = HiddenGroup (x, -4, 0, NULL);
      SetGroupSpacing (z, -1, 0);
      StaticPrompt (z, "Qualifier", 7 * stdCharWidth, 0, programFont, 'c');
      StaticPrompt (z, "Value", 10 * stdCharWidth, 0, programFont, 'c');
      ifp->gbquals = CreateQualsDialog (x, 5, -1, 7, 10);
    } else {
      ifp->gbquals = CreateImportFields (x, importFormTabs [0], sfp, allowProductGBQual);
    }
    AlignObjects (ALIGN_CENTER, (HANDLE) ifp->data, (HANDLE) x, NULL);
    ifp->pages [IMPORT_PAGE] = s;
    Hide (ifp->pages [IMPORT_PAGE]);
    importFormTabs [0] = NULL;

    s = HiddenGroup (h, -1, 0, NULL);
    CreateCommonFeatureGroup (s, (FeatureFormPtr) ifp, sfp, TRUE, TRUE);
    ifp->pages [COMMON_PAGE] = s;
    Hide (ifp->pages [COMMON_PAGE]);

    s = HiddenGroup (h, -1, 0, NULL);
    ifp->location = CreateIntervalEditorDialogEx (s, NULL, 4, 2, sep, TRUE, FALSE,
                                                  TRUE, TRUE, FALSE,
                                                  (FeatureFormPtr) ifp,
                                                  StdFeatIntEdPartialCallback);
    ifp->pages [LOCATION_PAGE] = s;
    Hide (ifp->pages [LOCATION_PAGE]);

    AlignObjects (ALIGN_CENTER, (HANDLE) ifp->pages [IMPORT_PAGE],
                  (HANDLE) ifp->pages [COMMON_PAGE], (HANDLE) ifp->pages [LOCATION_PAGE],
                  NULL);
    AlignObjects (ALIGN_CENTER, (HANDLE) ifp->foldertabs, (HANDLE) h, NULL);

    c = HiddenGroup (w, 2, 0, NULL);
    b = PushButton (c, "Accept", StdFeatFormAcceptButtonProc);
    SetObjectExtra (b, ifp, NULL);
    PushButton (c, "Cancel", StdCancelButtonProc);
    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
    RealizeWindow (w);

    SendMessageToDialog (ifp->data, VIB_MSG_INIT);
    SendMessageToDialog (ifp->location, VIB_MSG_INIT);
    Show (ifp->pages [ifp->currentPage]);
    SendMessageToDialog (ifp->data, VIB_MSG_ENTER);
    Update ();
  }
  return (ForM) w;
}

extern Int2 LIBCALLBACK ImportGenFunc (Pointer data)

{
  EnumFieldAssocPtr  ap;
  FeatDefPtr         curr;
  HelpMessageFunc    helpfunc;
  Int2               i;
  ImprtFormPtr       ifp;
  ImpFeatPtr         imp;
  ImportPagePtr      ipp;
  Uint1              key;
  CharPtr            label = NULL;
  ObjMgrPtr          omp;
  OMProcControlPtr   ompcp;
  ObjMgrTypePtr      omtp;
  OMUserDataPtr      omudp;
  ObjMgrProcPtr      proc;
  SeqEntryPtr        sep;
  SeqFeatPtr         sfp;
  Uint2              subtype;
  Char               title [64];
  WindoW             w;

  ompcp = (OMProcControlPtr) data;
  sfp = NULL;
  sep = NULL;
  imp = NULL;
  subtype = FEATDEF_IMP;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  proc = ompcp->proc;
  switch (ompcp->input_itemtype) {
    case OBJ_SEQFEAT :
      sfp = (SeqFeatPtr) ompcp->input_data;
      if (sfp != NULL && sfp->data.choice != SEQFEAT_IMP) {
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
    ifp = (ImprtFormPtr) omudp->userdata.ptrvalue;
    if (ifp != NULL) {
      Select (ifp->form);
    }
    return OM_MSG_RET_DONE;
  }
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  StringCpy (title, "Feature");
  if (sfp != NULL && sfp->data.value.ptrvalue != NULL) {
    imp = (ImpFeatPtr) sfp->data.value.ptrvalue;
    StringNCpy_0 (title, imp->key, sizeof (title));
    if (StringHasNoText (title)) {
      StringCpy (title, "Feature");
    }
  } else if (proc->subinputtype > 0) {
    subtype = proc->subinputtype;
    curr = FeatDefFindNext (NULL, &key, &label, FEATDEF_ANY, TRUE);
    while (curr != NULL) {
      if (key != FEATDEF_BAD && curr->seqfeat_key == SEQFEAT_IMP) {
        if (subtype == curr->featdef_key) {
          StringNCpy_0 (title, curr->typelabel, sizeof (title));
          break;
        }
      }
      curr = FeatDefFindNext (curr, &key, &label, FEATDEF_ANY, TRUE);
    }
  }
  w = (WindoW) CreateImportForm (-50, -33, title, sfp, sep,
                                 StdFeatFormActnProc);
  ifp = (ImprtFormPtr) GetObjectExtra (w);
  if (ifp != NULL) {
    ifp->input_entityID = ompcp->input_entityID;
    ifp->input_itemID = ompcp->input_itemID;
    ifp->input_itemtype = ompcp->input_itemtype;
    ifp->this_itemtype = OBJ_SEQFEAT;
    ifp->this_subtype = subtype;
    ifp->procid = ompcp->proc->procid;
    ifp->proctype = ompcp->proc->proctype;
    ifp->userkey = OMGetNextUserKey ();
    omudp = ObjMgrAddUserData (ompcp->input_entityID, ompcp->proc->procid,
	                           OMPROC_EDIT, ifp->userkey);
    if (omudp != NULL) {
      omudp->userdata.ptrvalue = (Pointer) ifp;
      omudp->messagefunc = StdVibrantEditorMsgFunc;
    }
    if (sfp != NULL) {
      omp = ObjMgrGet ();
      if (omp != NULL) {
        omtp = ObjMgrTypeFind (omp, OBJ_SEQFEAT, NULL, NULL);
        if (omtp != NULL && omtp->subtypefunc != NULL) {
          ifp->this_subtype = (*(omtp->subtypefunc)) (sfp);
        }
      }
    }
    SendMessageToForm (ifp->form, VIB_MSG_INIT);
    if (sfp != NULL) {
      PointerToForm (ifp->form, (Pointer) sfp);
      SetClosestParentIfDuplicating ((BaseFormPtr) ifp);
    } else {
      ipp = (ImportPagePtr) GetObjectExtra (ifp->data);
      if (ipp != NULL) {
        for (i = 1, ap = ipp->alist; ap->name != NULL; i++, ap++) {
          if (StringCmp (ap->name, title) == 0) {
            SetValue (ipp->key, i);
            SetTitle (ipp->loc, "");
            break;
          }
        }
      }
      SetNewFeatureDefaultInterval ((FeatureFormPtr) ifp);
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


typedef struct enumpage {
  DIALOG_MESSAGE_BLOCK
  PopuP              key;
  EnumFieldAssocPtr  alist;
} EnumPage, PNTR EnumPagePtr;

typedef struct enumform {
  FEATURE_FORM_BLOCK
  SeqEntryPtr   sep;
  GrouP         pages [NUM_PAGES];
  DialoG        foldertabs;
  Int2          currentPage;
} EnumForm, PNTR EnumFormPtr;

extern EnumFieldAssoc  enum_bond_alist [];
ENUM_ALIST(enum_bond_alist)
  {"Disulfide",              1},
  {"Thioester",              2},
  {"Crosslink",              3},
  {"Thioether",              4},
  {"Other",                255},
END_ENUM_ALIST

extern EnumFieldAssoc  enum_site_alist [];
ENUM_ALIST(enum_site_alist)
  {"Active",                       1},
  {"Binding",                      2},
  {"Cleavage",                     3},
  {"Inhibit",                      4},
  {"Modified",                     5},
  {"Glycosylation",                6},
  {"Myristoylation",               7},
  {"Mutagenized",                  8},
  {"Metal-binding",                9},
  {"Phosphorylation",             10},
  {"Acetylation",                 11},
  {"Amidation",                   12},
  {"Methylation",                 13},
  {"Hydroxylation",               14},
  {"Sulfatation",                 15},
  {"Oxidative-deamination",       16},
  {"Pyrrolidone-carboxylic-acid", 17},
  {"Gamma-carboxyglutamic-acid",  18},
  {"Blocked",                     19},
  {"Lipid-binding",               20},
  {"np-binding",                  21},
  {"DNA-binding",                 22},
  {"Signal-peptide",              23},
  {"Transit-peptide",             24},
  {"Transmembrane-region",        25},
  {"Other",                      255},
END_ENUM_ALIST

static ENUM_ALIST(enum_psec_alist)
  {"Helix",                  1},
  {"Sheet",                  2},
  {"Turn",                   3},
 END_ENUM_ALIST

static void EnumFeatPtrToEnumPage (DialoG d, Pointer data)

{
  EnumPagePtr  epp;
  Int4Ptr      intptr;
  Int4         val;

  epp = (EnumPagePtr) GetObjectExtra (d);
  intptr = (Int4Ptr) data;
  if (epp != NULL) {
    val = *intptr;
    SetEnumPopup (epp->key, epp->alist, (UIEnum) val);
  }
}

static Pointer EnumPageToEnumFeatPtr (DialoG d)

{
  EnumPagePtr  epp;
  Int4Ptr      intptr;
  UIEnum       val;

  intptr = NULL;
  epp = (EnumPagePtr) GetObjectExtra (d);
  if (epp != NULL) {
    if (GetEnumPopup (epp->key, epp->alist, &val)) {
      epp->intvalue = (Int4) val;
      intptr = (&epp->intvalue);
    }
  }
  return (Pointer) intptr;
}

static DialoG CreateEnumDialog (GrouP h, CharPtr title, Uint2 subtype)

{
  EnumPagePtr  epp;
  GrouP        f;
  GrouP        m;
  GrouP        p;
  GrouP        s;

  p = HiddenGroup (h, 1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  epp = (EnumPagePtr) MemNew (sizeof (EnumPage));
  if (epp != NULL) {

    SetObjectExtra (p, epp, StdCleanupExtraProc);
    epp->dialog = (DialoG) p;
    epp->todialog = EnumFeatPtrToEnumPage;
    epp->fromdialog = EnumPageToEnumFeatPtr;
    epp->testdialog = NULL;

    if (title != NULL && title [0] != '\0') {
      s = NormalGroup (p, 0, -2, title, systemFont, NULL);
    } else {
      s = HiddenGroup (p, 0, -2, NULL);
    }
    m = HiddenGroup (s, -1, 0, NULL);
    /*
    SetGroupSpacing (m, 10, 10);
    */

    switch (subtype) {
      case FEATDEF_BOND :
        epp->alist = enum_bond_alist;
        break;
      case FEATDEF_SITE :
        epp->alist = enum_site_alist;
        break;
      case FEATDEF_PSEC_STR :
        epp->alist = enum_psec_alist;
        break;
      default :
        Message (MSG_FATAL, "Unknown enumerated feature type");
        return NULL;
    }

    f = HiddenGroup (m, -2, 0, NULL);
    StaticPrompt (f, "Type", 0, popupMenuHeight, programFont, 'l');
    epp->key = PopupList (f, TRUE, NULL);
    SetObjectExtra (epp->key, epp, NULL);
    InitEnumPopup (epp->key, epp->alist, NULL);
    /*
    SetEnumPopup (epp->key, epp->alist, (UIEnum) 0);
    */
  }

  return (DialoG) p;
}

static void SetEnumImportExportItems (EnumFormPtr efp)

{
  IteM  exportItm;
  IteM  importItm;

  if (efp != NULL) {
    importItm = FindFormMenuItem ((BaseFormPtr) efp, VIB_MSG_IMPORT);
    exportItm = FindFormMenuItem ((BaseFormPtr) efp, VIB_MSG_EXPORT);
    switch (efp->currentPage) {
      case ENUM_PAGE :
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

static void ChangeEnumPage (VoidPtr data, Int2 newval, Int2 oldval)

{
  EnumFormPtr  efp;

  efp = (EnumFormPtr) data;
  if (efp != NULL) {
    efp->currentPage = newval;
    SafeHide (efp->pages [oldval]);
    SafeShow (efp->pages [newval]);
    switch (newval) {
      case ENUM_PAGE :
        SendMessageToDialog (efp->data, VIB_MSG_ENTER);
        break;
      case COMMON_PAGE :
        break;
      case LOCATION_PAGE :
        SendMessageToDialog (efp->location, VIB_MSG_ENTER);
        break;
      default :
        break;
    }
    SetEnumImportExportItems (efp);
    Update ();
  }
}

static Boolean ImportEnumForm (ForM f, CharPtr filename)

{
  EnumFormPtr  efp;

  efp = (EnumFormPtr) GetObjectExtra (f);
  if (efp != NULL) {
    switch (efp->currentPage) {
      case LOCATION_PAGE :
        return ImportDialog (efp->location, filename);
      default :
        break;
    }
  }
  return FALSE;
}

static Boolean ExportEnumForm (ForM f, CharPtr filename)

{
  EnumFormPtr  efp;

  efp = (EnumFormPtr) GetObjectExtra (f);
  if (efp != NULL) {
    switch (efp->currentPage) {
      case LOCATION_PAGE :
        return ExportDialog (efp->location, filename);
      default :
        break;
    }
  }
  return FALSE;
}

static CharPtr  enumFormTabs [] = {
  NULL, "Properties", "Location", NULL
};

static void EnumFormMessage (ForM f, Int2 mssg)

{
  EnumFormPtr  efp;

  efp = (EnumFormPtr) GetObjectExtra (f);
  if (efp != NULL) {
    switch (mssg) {
      case VIB_MSG_INIT :
        StdInitFeatFormProc (f);
        break;
      case VIB_MSG_IMPORT :
        ImportEnumForm (f, NULL);
        break;
      case VIB_MSG_EXPORT :
        ExportEnumForm (f, NULL);
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
        if (efp->currentPage == LOCATION_PAGE) {
          PointerToDialog (efp->location, NULL);
        } else {
          StdDeleteTextProc (NULL);
        }
        break;
      default :
        if (efp->appmessage != NULL) {
          efp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static void EnumFormActivate (WindoW w)

{
  EnumFormPtr  efp;

  efp = (EnumFormPtr) GetObjectExtra (w);
  if (efp != NULL) {
    if (efp->activate != NULL) {
      efp->activate (w);
    }
    SetEnumImportExportItems (efp);
  }
}

extern DialoG CreateBondEditorDialog (GrouP h, CharPtr title, SeqEntryPtr sep);

extern ForM CreateEnumForm (Int2 left, Int2 top, CharPtr title,
                            SeqFeatPtr sfp, SeqEntryPtr sep,
                            Uint2 subtype, FormActnFunc actproc)

{
  ButtoN             b;
  GrouP              c;
  EnumFormPtr        efp;
  GrouP              g;
  GrouP              h;
  GrouP              s;
  StdEditorProcsPtr  sepp;
  WindoW             w;

  w = NULL;
  efp = (EnumFormPtr) MemNew (sizeof (EnumForm));
  if (efp != NULL) {
    w = FixedWindow (left, top, -10, -10, title, StdCloseWindowProc);
    SetObjectExtra (w, efp, StdFeatFormCleanupProc);
    efp->form = (ForM) w;
    efp->actproc = actproc;
    efp->toform = StdSeqFeatPtrToFeatFormProc;
    efp->fromform = NULL;
    efp->formmessage = EnumFormMessage;
    efp->testform = NULL;
    efp->importform = ImportEnumForm;
    efp->exportform = ExportEnumForm;

#ifndef WIN_MAC
    CreateStdEditorFormMenus (w);
#endif

    efp->activate = NULL;
    sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
    if (sepp != NULL) {
      efp->activate = sepp->activateForm;
      efp->appmessage = sepp->handleMessages;
    }
    SetActivate (w, EnumFormActivate);

    g = HiddenGroup (w, -1, 0, NULL);
    SetGroupSpacing (g, 3, 10);

    efp->sep = sep;
    enumFormTabs [0] = NULL;
    if (title != NULL && *title != '\0') {
      enumFormTabs [0] = title;
    } else if (subtype == FEATDEF_BOND) {
      enumFormTabs [0] = "Bond";
    } else if (subtype == FEATDEF_SITE) {
      enumFormTabs [0] = "Site";
    } else if (subtype == FEATDEF_PSEC_STR) {
      enumFormTabs [0] = "Secondary Structure";
    }
    efp->foldertabs = CreateFolderTabs (g, enumFormTabs, ENUM_PAGE,
                                        0, 0, SYSTEM_FOLDER_TAB,
                                        ChangeEnumPage, (Pointer) efp);
    efp->currentPage = ENUM_PAGE;

    h = HiddenGroup (g, 0, 0, NULL);

    s = HiddenGroup (h, -1, 0, NULL);
    SetGroupSpacing (s, 3, 10);
    efp->data = CreateEnumDialog (s, NULL, subtype);
    efp->pages [ENUM_PAGE] = s;
    Hide (efp->pages [ENUM_PAGE]);
    enumFormTabs [0] = NULL;

    s = HiddenGroup (h, -1, 0, NULL);
    CreateCommonFeatureGroup (s, (FeatureFormPtr) efp, sfp, TRUE, TRUE);
    efp->pages [COMMON_PAGE] = s;
    Hide (efp->pages [COMMON_PAGE]);

    s = HiddenGroup (h, -1, 0, NULL);
    if (subtype == FEATDEF_BOND) {
      efp->location = CreateBondEditorDialog (s, NULL, sep);
    } else {
      efp->location = CreateIntervalEditorDialogEx (s, NULL, 4, 2, sep, TRUE, TRUE,
                                                    TRUE, TRUE, FALSE,
                                                    (FeatureFormPtr) efp,
                                                    StdFeatIntEdPartialCallback);
    }
    efp->pages [LOCATION_PAGE] = s;
    Hide (efp->pages [LOCATION_PAGE]);

    AlignObjects (ALIGN_CENTER, (HANDLE) efp->pages [ENUM_PAGE],
                  (HANDLE) efp->pages [COMMON_PAGE], (HANDLE) efp->pages [LOCATION_PAGE],
                  NULL);
    AlignObjects (ALIGN_CENTER, (HANDLE) efp->foldertabs, (HANDLE) h, NULL);

    c = HiddenGroup (w, 2, 0, NULL);
    b = PushButton (c, "Accept", StdFeatFormAcceptButtonProc);
    SetObjectExtra (b, efp, NULL);
    PushButton (c, "Cancel", StdCancelButtonProc);
    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
    RealizeWindow (w);

    SendMessageToDialog (efp->data, VIB_MSG_INIT);
    SendMessageToDialog (efp->location, VIB_MSG_INIT);
    Show (efp->pages [efp->currentPage]);
    SendMessageToDialog (efp->data, VIB_MSG_ENTER);
    Update ();
  }
  return (ForM) w;
}

extern Int2 LIBCALLBACK EnumGenFunc (Pointer data)

{
  EnumFormPtr       efp;
  HelpMessageFunc   helpfunc;
  OMProcControlPtr  ompcp;
  OMUserDataPtr     omudp;
  ObjMgrProcPtr     proc;
  SeqEntryPtr       sep;
  SeqFeatPtr        sfp;
  Uint2             subtype;
  WindoW            w;

  ompcp = (OMProcControlPtr) data;
  w = NULL;
  sfp = NULL;
  sep = NULL;
  subtype = 0;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  proc = ompcp->proc;
  switch (ompcp->input_itemtype)
  {
    case OBJ_SEQFEAT :
      sfp = (SeqFeatPtr) ompcp->input_data;
      if (sfp != NULL &&
         (sfp->data.choice != SEQFEAT_BOND &&
          sfp->data.choice != SEQFEAT_SITE &&
          sfp->data.choice != SEQFEAT_PSEC_STR)) {
        return OM_MSG_RET_ERROR;
      }
      if (sfp->data.choice == SEQFEAT_BOND) {
        subtype = FEATDEF_BOND;
      } else if (sfp->data.choice == SEQFEAT_SITE) {
        subtype = FEATDEF_SITE;
      } else if (sfp->data.choice == SEQFEAT_PSEC_STR) {
        subtype = FEATDEF_PSEC_STR;
      } else {
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
    efp = (EnumFormPtr) omudp->userdata.ptrvalue;
    if (efp != NULL) {
      Select (efp->form);
    }
    return OM_MSG_RET_DONE;
  }
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  if (sfp == NULL) {
    subtype = proc->subinputtype;
  }
  if (subtype == FEATDEF_BOND) {
    w = (WindoW) CreateEnumForm (-50, -33, "Bond",
                                 sfp, sep, subtype,
                                 StdFeatFormActnProc);
  } else if (subtype == FEATDEF_SITE) {
    w = (WindoW) CreateEnumForm (-50, -33, "Site",
                                 sfp, sep, subtype,
                                 StdFeatFormActnProc);
  } else if (subtype == FEATDEF_PSEC_STR) {
    w = (WindoW) CreateEnumForm (-50, -33, "Secondary Structure",
                                 sfp, sep, subtype,
                                 StdFeatFormActnProc);
  } else {
    return OM_MSG_RET_ERROR;
  }
  efp = (EnumFormPtr) GetObjectExtra (w);
  if (efp != NULL) {
    efp->input_entityID = ompcp->input_entityID;
    efp->input_itemID = ompcp->input_itemID;
    efp->input_itemtype = ompcp->input_itemtype;
    efp->this_itemtype = OBJ_SEQFEAT;
    efp->this_subtype = subtype;
    efp->procid = ompcp->proc->procid;
    efp->proctype = ompcp->proc->proctype;
    efp->userkey = OMGetNextUserKey ();
    omudp = ObjMgrAddUserData (ompcp->input_entityID, ompcp->proc->procid,
	                           OMPROC_EDIT, efp->userkey);
    if (omudp != NULL) {
      omudp->userdata.ptrvalue = (Pointer) efp;
      omudp->messagefunc = StdVibrantEditorMsgFunc;
    }

    SendMessageToForm (efp->form, VIB_MSG_INIT);
    if (sfp != NULL) {
      PointerToForm (efp->form, (Pointer) sfp);
      SetClosestParentIfDuplicating ((BaseFormPtr) efp);
    } else {
      SetNewFeatureDefaultInterval ((FeatureFormPtr) efp);
    }
  }
  Show (w);
  Select (w);
  helpfunc = (HelpMessageFunc) GetAppProperty ("HelpMessageProc");
  if (helpfunc != NULL) {
    if (subtype == FEATDEF_BOND) {
      helpfunc ("Features", "Bond");
    } else if (subtype == FEATDEF_SITE) {
      helpfunc ("Features", "Site");
    } else if (subtype == FEATDEF_PSEC_STR) {
      helpfunc ("Features", "Secondary Structure");
    }
  }
  return OM_MSG_RET_DONE;
}


typedef struct regionform {
  FEATURE_FORM_BLOCK
  SeqEntryPtr   sep;
  GrouP         pages [NUM_PAGES];
  DialoG        foldertabs;
  Int2          currentPage;
} RegionForm, PNTR RegionFormPtr;

typedef struct regionpage {
  DIALOG_MESSAGE_BLOCK
  TexT           region;
  ButtoN         convertToMiscFeat;
  RegionFormPtr  rfp;
} RegionPage, PNTR RegionPagePtr;

static void CharPtrToRegionPage (DialoG d, Pointer data)

{
  RegionFormPtr  rfp;
  RegionPagePtr  rpp;
  CharPtr        str;

  rpp = (RegionPagePtr) GetObjectExtra (d);
  str = (CharPtr) data;
  if (rpp != NULL) {
    rfp = rpp->rfp;
    if (rfp->this_subtype == FEATDEF_COMMENT) return;
    if (rfp->this_subtype == FEATDEF_misc_feature) return;
    if (str != NULL) {
      SafeSetTitle (rpp->region, str);
    } else {
      SafeSetTitle (rpp->region, "");
    }
  }
}

static Pointer RegionPageToCharPtr (DialoG d)

{
  ImpFeatPtr     ifp;
  RegionFormPtr  rfp;
  RegionPagePtr  rpp;
  CharPtr        str;

  str = NULL;
  rpp = (RegionPagePtr) GetObjectExtra (d);
  if (rpp != NULL) {
    rfp = rpp->rfp;
    if (rfp->this_subtype == FEATDEF_misc_feature) {
      if (GetStatus (rpp->convertToMiscFeat)) {
        ifp = ImpFeatNew ();
        if (ifp != NULL) {
          ifp->key = StringSave ("misc_feature");
        }
        return (Pointer) ifp;
      } else {
        return NULL;
      }
    }
    if (rfp->this_subtype == FEATDEF_COMMENT) return NULL;
    str = SaveStringFromText (rpp->region);
    if (str == NULL) {
      str = StringSave ("");
    }
  }
  return (Pointer) str;
}

static DialoG CreateRegionDialog (GrouP h, CharPtr title, CharPtr str,
                                  Uint2 subtype, RegionFormPtr rfp)

{
  GrouP          m;
  GrouP          p;
  RegionPagePtr  rpp;
  GrouP          s;

  p = HiddenGroup (h, 1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  rpp = (RegionPagePtr) MemNew (sizeof (RegionPage));
  if (rpp != NULL) {

    SetObjectExtra (p, rpp, StdCleanupExtraProc);
    rpp->dialog = (DialoG) p;
    rpp->todialog = CharPtrToRegionPage;
    rpp->fromdialog = RegionPageToCharPtr;
    rpp->testdialog = NULL;

    rpp->rfp = rfp;

    if (title != NULL && title [0] != '\0') {
      s = NormalGroup (p, 0, -2, title, systemFont, NULL);
    } else {
      s = HiddenGroup (p, 0, -2, NULL);
    }
    m = HiddenGroup (s, 2, 0, NULL);

    if (subtype == FEATDEF_REGION) {
      StaticPrompt (m, "Name", 0, dialogTextHeight, programFont, 'c');
      rpp->region = DialogText (m, "", 20, NULL);
    } else if (subtype == FEATDEF_COMMENT) {
      StaticPrompt (m, "Enter comment in Properties page", 0, 0, programFont, 'c');
    }
    if (GetAppProperty ("InternalNcbiSequin") != NULL) {
      rpp->convertToMiscFeat = CheckBox (p, "Convert to misc_feat", NULL);
      AlignObjects (ALIGN_CENTER, (HANDLE) s, (HANDLE) rpp->convertToMiscFeat, NULL);
    }
  }

  return (DialoG) p;
}

static void SetRegionImportExportItems (RegionFormPtr rfp)

{
  IteM  exportItm;
  IteM  importItm;

  if (rfp != NULL) {
    importItm = FindFormMenuItem ((BaseFormPtr) rfp, VIB_MSG_IMPORT);
    exportItm = FindFormMenuItem ((BaseFormPtr) rfp, VIB_MSG_EXPORT);
    switch (rfp->currentPage) {
      case REGION_PAGE :
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

static void ChangeRegionPage (VoidPtr data, Int2 newval, Int2 oldval)

{
  RegionFormPtr  rfp;

  rfp = (RegionFormPtr) data;
  if (rfp != NULL) {
    rfp->currentPage = newval;
    SafeHide (rfp->pages [oldval]);
    SafeShow (rfp->pages [newval]);
    switch (newval) {
      case REGION_PAGE :
        break;
      case COMMON_PAGE :
        break;
      case LOCATION_PAGE :
        SendMessageToDialog (rfp->location, VIB_MSG_ENTER);
        break;
      default :
        break;
    }
    SetRegionImportExportItems (rfp);
    Update ();
  }
}

static Boolean ImportRegionForm (ForM f, CharPtr filename)

{
  RegionFormPtr  rfp;

  rfp = (RegionFormPtr) GetObjectExtra (f);
  if (rfp != NULL) {
    switch (rfp->currentPage) {
      case LOCATION_PAGE :
        return ImportDialog (rfp->location, filename);
      default :
        break;
    }
  }
  return FALSE;
}

static Boolean ExportRegionForm (ForM f, CharPtr filename)

{
  RegionFormPtr  rfp;

  rfp = (RegionFormPtr) GetObjectExtra (f);
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

static CharPtr  regionFormTabs [] = {
  "Region", "Properties", "Location", NULL
};

static CharPtr  commentFormTabs [] = {
  "Comment", "Properties", "Location", NULL
};

static void RegionFormMessage (ForM f, Int2 mssg)

{
  RegionFormPtr  rfp;

  rfp = (RegionFormPtr) GetObjectExtra (f);
  if (rfp != NULL) {
    switch (mssg) {
      case VIB_MSG_INIT :
        StdInitFeatFormProc (f);
        break;
      case VIB_MSG_IMPORT :
        ImportRegionForm (f, NULL);
        break;
      case VIB_MSG_EXPORT :
        ExportRegionForm (f, NULL);
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

static void RegionFormActivate (WindoW w)

{
  RegionFormPtr  rfp;

  rfp = (RegionFormPtr) GetObjectExtra (w);
  if (rfp != NULL) {
    if (rfp->activate != NULL) {
      rfp->activate (w);
    }
    SetRegionImportExportItems (rfp);
  }
}

extern ForM CreateRegionOrCommentForm (Int2 left, Int2 top, CharPtr title,
                                       SeqFeatPtr sfp, SeqEntryPtr sep,
                                       Uint2 subtype, FormActnFunc actproc)

{
  ButtoN             b;
  GrouP              c;
  GrouP              g;
  GrouP              h;
  RegionFormPtr      rfp;
  GrouP              s;
  StdEditorProcsPtr  sepp;
  CharPtr            str;
  WindoW             w;

  w = NULL;
  rfp = (RegionFormPtr) MemNew (sizeof (RegionForm));
  if (rfp != NULL) {
    w = FixedWindow (left, top, -10, -10, title, StdCloseWindowProc);
    SetObjectExtra (w, rfp, StdFeatFormCleanupProc);
    rfp->form = (ForM) w;
    rfp->actproc = actproc;
    rfp->toform = StdSeqFeatPtrToFeatFormProc;
    rfp->fromform = NULL;
    rfp->formmessage = RegionFormMessage;
    rfp->testform = NULL;
    rfp->importform = ImportRegionForm;
    rfp->exportform = ExportRegionForm;

#ifndef WIN_MAC
    CreateStdEditorFormMenus (w);
#endif

    rfp->activate = NULL;
    sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
    if (sepp != NULL) {
      rfp->activate = sepp->activateForm;
      rfp->appmessage = sepp->handleMessages;
    }
    SetActivate (w, RegionFormActivate);

    g = HiddenGroup (w, -1, 0, NULL);
    SetGroupSpacing (g, 3, 10);

    rfp->sep = sep;
    if (subtype == FEATDEF_REGION) {
      rfp->foldertabs = CreateFolderTabs (g, regionFormTabs, REGION_PAGE,
                                          0, 0, SYSTEM_FOLDER_TAB,
                                          ChangeRegionPage, (Pointer) rfp);
    } else if (subtype == FEATDEF_COMMENT) {
      rfp->foldertabs = CreateFolderTabs (g, commentFormTabs, REGION_PAGE,
                                          0, 0, SYSTEM_FOLDER_TAB,
                                          ChangeRegionPage, (Pointer) rfp);
    } else {
      Message (MSG_FATAL, "Region or Comment Feature error");
    }
    rfp->currentPage = REGION_PAGE;

    h = HiddenGroup (g, 0, 0, NULL);

    s = HiddenGroup (h, -1, 0, NULL);
    str = NULL;
    if (sfp != NULL && sfp->data.choice == SEQFEAT_REGION) {
      str = sfp->data.value.ptrvalue;
    }
    rfp->data = CreateRegionDialog (s, NULL, str, subtype, rfp);
    /*
    if (subtype == FEATDEF_COMMENT) {
      rfp->data = NULL;
    }
    */
    rfp->pages [REGION_PAGE] = s;
    Hide (rfp->pages [REGION_PAGE]);

    s = HiddenGroup (h, -1, 0, NULL);
    CreateCommonFeatureGroup (s, (FeatureFormPtr) rfp, sfp, TRUE, TRUE);
    rfp->pages [COMMON_PAGE] = s;
    Hide (rfp->pages [COMMON_PAGE]);

    s = HiddenGroup (h, -1, 0, NULL);
    rfp->location = CreateIntervalEditorDialogEx (s, NULL, 4, 2, sep, TRUE, TRUE,
                                                  TRUE, TRUE, FALSE,
                                                  (FeatureFormPtr) rfp,
                                                  StdFeatIntEdPartialCallback);
    rfp->pages [LOCATION_PAGE] = s;
    Hide (rfp->pages [LOCATION_PAGE]);

    AlignObjects (ALIGN_CENTER, (HANDLE) rfp->pages [REGION_PAGE],
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
    Show (rfp->pages [rfp->currentPage]);
    SendMessageToDialog (rfp->data, VIB_MSG_ENTER);
    Update ();
  }
  return (ForM) w;
}

static void RegionOrCommentFeatFormActnProc (ForM f)

{
  RegionFormPtr  rfp;
  RegionPagePtr  rpp;

  rfp = (RegionFormPtr) GetObjectExtra (f);
  if (rfp != NULL) {
    rpp = (RegionPagePtr) GetObjectExtra (rfp->data);
    if (rpp != NULL && GetStatus (rpp->convertToMiscFeat)) {
      rfp->this_subtype = FEATDEF_misc_feature;
    }
  }
  if (FeatFormReplaceWithoutUpdateProc (f)) {
    if (rfp != NULL) {
      GetRidOfEmptyFeatsDescStrings (rfp->input_entityID, NULL);
      if (GetAppProperty ("InternalNcbiSequin") != NULL) {
        ExtendGeneFeatIfOnMRNA (rfp->input_entityID, NULL);
      }
      ObjMgrSendMsg (OM_MSG_UPDATE, rfp->input_entityID,
                     rfp->input_itemID, rfp->input_itemtype);
    }
  }
}

extern Int2 LIBCALLBACK RegionOrCommentGenFunc (Pointer data)

{
  OMProcControlPtr  ompcp;
  OMUserDataPtr     omudp;
  ObjMgrProcPtr     proc;
  RegionFormPtr     rfp;
  SeqEntryPtr       sep;
  SeqFeatPtr        sfp;
  Uint2             subtype;
  CharPtr           title;
  WindoW            w;

  ompcp = (OMProcControlPtr) data;
  sfp = NULL;
  sep = NULL;
  subtype = 0;
  title = NULL;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  proc = ompcp->proc;
  switch (ompcp->input_itemtype) {
    case OBJ_SEQFEAT :
      sfp = (SeqFeatPtr) ompcp->input_data;
      if (sfp != NULL &&
          (sfp->data.choice != SEQFEAT_REGION &&
           sfp->data.choice != SEQFEAT_COMMENT)) {
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
    rfp = (RegionFormPtr) omudp->userdata.ptrvalue;
    if (rfp != NULL) {
      Select (rfp->form);
    }
    return OM_MSG_RET_DONE;
  }
  if (sfp != NULL) {
    if (sfp->data.choice == SEQFEAT_REGION) {
      subtype = FEATDEF_REGION;
    } else if (sfp->data.choice == SEQFEAT_COMMENT) {
      subtype = FEATDEF_COMMENT;
    }
  } else if (proc->subinputtype > 0) {
    subtype = proc->subinputtype;
  }
  if (subtype == FEATDEF_REGION) {
    title = "Region";
  } else if (subtype == FEATDEF_COMMENT) {
    title = "Comment";
  } else {
    return OM_MSG_RET_ERROR;
  }
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  w = (WindoW) CreateRegionOrCommentForm (-50, -33, title, sfp, sep,
                                          subtype, RegionOrCommentFeatFormActnProc);
  rfp = (RegionFormPtr) GetObjectExtra (w);
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
    SendMessageToForm (rfp->form, VIB_MSG_INIT);
    if (sfp != NULL) {
      PointerToForm (rfp->form, (Pointer) sfp);
      SetClosestParentIfDuplicating ((BaseFormPtr) rfp);
    } else {
      SetNewFeatureDefaultInterval ((FeatureFormPtr) rfp);
    }
  }
  Show (w);
  Select (w);
  return OM_MSG_RET_DONE;
}


typedef struct molinfopage {
  DIALOG_MESSAGE_BLOCK
  PopuP           moltype;
  PopuP           technique;
  PopuP           complete;
  GrouP           expGrp;
  TexT            explain;
  EnumFieldAssoc  PNTR moltypeAlist;
  EnumFieldAssoc  PNTR techniqueAlist;
  EnumFieldAssoc  PNTR completeAlist;
} MolInfoPage, PNTR MolInfoPagePtr;

typedef struct molinfoform {
  DESCRIPTOR_FORM_BLOCK
  PopuP           molPopup;
  PopuP           topologyPopup;
  PopuP           strandPopup;
} MolInfoForm, PNTR MolInfoFormPtr;

static ENUM_ALIST(molinfo_biomol_alist)
  {" ",                      0},
  {"Genomic DNA or RNA",     1},
  {"Precursor RNA",          2},
  {"mRNA [cDNA]",            3},
  {"Ribosomal RNA",          4},
  {"Transfer RNA",           5},
  {"Small nuclear RNA",      6},
  {"Small cytoplasmic RNA",  7},
  {"Peptide",                8},
  {"Other-Genetic",          9},
  {"Genomic-mRNA",          10},
  {"Other",                255},
END_ENUM_ALIST

static ENUM_ALIST(molinfo_biomol_nuc_alist)
  {" ",                      0},
  {"Genomic DNA or RNA",     1},
  {"Precursor RNA",          2},
  {"mRNA [cDNA]",            3},
  {"Ribosomal RNA",          4},
  {"Transfer RNA",           5},
  {"Small nuclear RNA",      6},
  {"Small cytoplasmic RNA",  7},
  {"Other-Genetic",          9},
  {"Genomic-mRNA",          10},
  {"Other",                255},
END_ENUM_ALIST

static ENUM_ALIST(molinfo_biomol_nucX_alist)
  {" ",                      0},
  {"Genomic DNA",          253},
  {"Genomic RNA",          254},
  {"Precursor RNA",          2},
  {"mRNA [cDNA]",            3},
  {"Ribosomal RNA",          4},
  {"Transfer RNA",           5},
  {"Small nuclear RNA",      6},
  {"Small cytoplasmic RNA",  7},
  {"Other-Genetic",          9},
  {"Genomic-mRNA",          10},
  {"Other",                255},
END_ENUM_ALIST

static ENUM_ALIST(molinfo_biomol_prot_alist)
  {" ",                      0},
  {"Peptide",                8},
  {"Other",                255},
END_ENUM_ALIST

static ENUM_ALIST(molinfo_tech_alist)
  {" ",                  0},
  {"Standard",           1},
  {"EST",                2},
  {"STS",                3},
  {"Survey",             4},
  {"Genetic Map",        5},
  {"Physical Map",       6},
  {"Derived",            7},
  {"Concept-Trans",      8},
  {"Seq-Pept",           9},
  {"Both",              10},
  {"Seq-Pept-Overlap",  11},
  {"Seq-Pept-Homol",    12},
  {"Concept-Trans-A",   13},
  {"HTGS 1",            14},
  {"HTGS 2",            15},
  {"HTGS 3",            16},
  {"FLI_cDNA",          17},
  {"Other:",            255},
END_ENUM_ALIST

static ENUM_ALIST(molinfo_tech_nuc_alist)
  {" ",                  0},
  {"Standard",           1},
  {"EST",                2},
  {"STS",                3},
  {"Survey",             4},
  {"Genetic Map",        5},
  {"Physical Map",       6},
  {"Derived",            7},
  {"HTGS 1",            14},
  {"HTGS 2",            15},
  {"HTGS 3",            16},
  {"FLI_cDNA",          17},
  {"Other:",            255},
END_ENUM_ALIST

static ENUM_ALIST(molinfo_tech_prot_alist)
  {" ",                  0},
  {"Concept-Trans",      8},
  {"Seq-Pept",           9},
  {"Both",              10},
  {"Seq-Pept-Overlap",  11},
  {"Seq-Pept-Homol",    12},
  {"Concept-Trans-A",   13},
  {"Other:",            255},
END_ENUM_ALIST

static ENUM_ALIST(molinfo_complete_alist)
  {" ",         0},
  {"Complete",  1},
  {"Partial",   2},
  {"No Left",   3},
  {"No Right",  4},
  {"No Ends",   5},
  {"Has Left",  6},
  {"Has Right", 7},
  {"Other",   255},
END_ENUM_ALIST

static ENUM_ALIST(molinfo_complete_nuc_alist)
  {" ",         0},
  {"Complete",  1},
  {"Partial",   2},
  {"No 5'",     3},
  {"No 3'",     4},
  {"No Ends",   5},
  {"Has 5'",    6},
  {"Has 3'",    7},
  {"Other",   255},
END_ENUM_ALIST

static ENUM_ALIST(molinfo_complete_prot_alist)
  {" ",         0},
  {"Complete",  1},
  {"Partial",   2},
  {"No NH3",    3},
  {"No CO2",    4},
  {"No Ends",   5},
  {"Has NH3",   6},
  {"Has CO2",   7},
  {"Other",   255},
END_ENUM_ALIST

static Uint1 check_biomol (Uint1 biomol)

{
  if (biomol > 10 && biomol < 253) return 0;
  return biomol;
}

static Uint1 check_technique (Uint1 tech)

{
  if (tech > 17 && tech != 255) return 0;
  return tech;
}

static Uint1 check_complete (Uint1 completeness)

{
  if (completeness > 7 && completeness != 255) return 0;
  return completeness;
}

static void MolInfoPtrToMolInfoPage (DialoG d, Pointer data)

{
  MolInfoPtr      mip;
  MolInfoPagePtr  mpp;

  mpp = (MolInfoPagePtr) GetObjectExtra (d);
  mip = (MolInfoPtr) data;
  if (mpp != NULL) {
    SetEnumPopup (mpp->moltype, mpp->moltypeAlist,
                  (UIEnum) check_biomol (mip->biomol));
    SetEnumPopup (mpp->technique, mpp->techniqueAlist,
                  (UIEnum) check_technique (mip->tech));
    SetEnumPopup (mpp->complete, mpp->completeAlist,
                  (UIEnum) check_complete (mip->completeness));
    SetTitle (mpp->explain, mip->techexp);
    if (mip->tech == 255) {
      SafeShow (mpp->expGrp);
    } else {
      SafeHide (mpp->expGrp);
    }
  }
}

static Pointer MolInfoPageToMolInfoPtr (DialoG d)

{
  MolInfoPtr      mip;
  MolInfoPagePtr  mpp;
  UIEnum          val;

  mip = NULL;
  mpp = (MolInfoPagePtr) GetObjectExtra (d);
  if (mpp != NULL) {
    mip = MolInfoNew ();
    if (mip != NULL) {
      if (GetEnumPopup (mpp->moltype, mpp->moltypeAlist, &val)) {
        mip->biomol = (Uint1) val;
      }
      if (GetEnumPopup (mpp->technique, mpp->techniqueAlist, &val)) {
        mip->tech = (Uint1) val;
      }
      if (GetEnumPopup (mpp->complete, mpp->completeAlist, &val)) {
        mip->completeness = (Uint1) val;
      }
      if (mip->tech == 255) {
        mip->techexp = SaveStringFromText (mpp->explain);
      }
    }
  }
  return (Pointer) mip;
}

static void ChangeTech (PopuP p)

{
  MolInfoPagePtr  mpp;
  UIEnum          val;

  mpp = (MolInfoPagePtr) GetObjectExtra (p);
  if (mpp != NULL) {
    if (GetEnumPopup (mpp->technique, mpp->techniqueAlist, &val)) {
      if (val == 255) {
        SafeShow (mpp->expGrp);
      } else {
        SafeHide (mpp->expGrp);
      }
    }
  }
}

static CharPtr  labels [] = {
  "Molecule", "Technique", NULL
};

extern DialoG CreateMolInfoDialog (GrouP h, CharPtr title, Uint1 biomol, Uint1 tech,
                                   Boolean showComplete, Boolean nucsOK, Boolean protsOK,
                                   Boolean separateGenomicDNAandRNA)

{
  GrouP           f;
  GrouP           m;
  MolInfoPagePtr  mpp;
  GrouP           p;
  GrouP           s;
  Int2            wid;
  GrouP           x;

  p = HiddenGroup (h, 1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  mpp = (MolInfoPagePtr) MemNew (sizeof (MolInfoPage));
  if (mpp != NULL) {

    SetObjectExtra (p, mpp, StdCleanupExtraProc);
    mpp->dialog = (DialoG) p;
    mpp->todialog = MolInfoPtrToMolInfoPage;
    mpp->fromdialog = MolInfoPageToMolInfoPtr;
    mpp->testdialog = NULL;

    if (title != NULL && title [0] != '\0') {
      s = NormalGroup (p, 0, -2, title, systemFont, NULL);
    } else {
      s = HiddenGroup (p, 0, -2, NULL);
    }
    m = HiddenGroup (s, -1, 0, NULL);
    /*
    SetGroupSpacing (m, 10, 10);
    */

    if (showComplete) {
      f = HiddenGroup (m, -3, 0, NULL);
    } else {
      f = HiddenGroup (m, -2, 0, NULL);
    }

    SelectFont (programFont);
    wid = MaxStringWidths (labels);

    if (nucsOK && protsOK) {
      mpp->moltypeAlist = molinfo_biomol_alist;
      mpp->techniqueAlist = molinfo_tech_alist;
      mpp->completeAlist = molinfo_complete_alist;
    } else if (nucsOK) {
      if (separateGenomicDNAandRNA) {
        mpp->moltypeAlist = molinfo_biomol_nucX_alist;
        mpp->techniqueAlist = molinfo_tech_nuc_alist;
        mpp->completeAlist = molinfo_complete_nuc_alist;
      } else {
        mpp->moltypeAlist = molinfo_biomol_nuc_alist;
        mpp->techniqueAlist = molinfo_tech_nuc_alist;
        mpp->completeAlist = molinfo_complete_nuc_alist;
      }
    } else if (protsOK) {
      mpp->moltypeAlist = molinfo_biomol_prot_alist;
      mpp->techniqueAlist = molinfo_tech_prot_alist;
      mpp->completeAlist = molinfo_complete_prot_alist;
    } else {
      mpp->moltypeAlist = molinfo_biomol_alist;
      mpp->techniqueAlist = molinfo_tech_alist;
      mpp->completeAlist = molinfo_complete_alist;
    }

    StaticPrompt (f, "Molecule",
                  wid, popupMenuHeight, programFont, 'l');
    mpp->moltype = PopupList (f, TRUE, NULL);
    SetObjectExtra (mpp->moltype, mpp, NULL);
    InitEnumPopup (mpp->moltype, mpp->moltypeAlist, NULL);
    SetEnumPopup (mpp->moltype, mpp->moltypeAlist, biomol);

    x = NULL;
    if (showComplete) {
      x = HiddenGroup (f, 2, 0, NULL);
      StaticPrompt (x, "Completedness",
                    0, popupMenuHeight, programFont, 'l');
      mpp->complete = PopupList (x, TRUE, NULL);
      SetObjectExtra (mpp->complete, mpp, NULL);
      InitEnumPopup (mpp->complete, mpp->completeAlist, NULL);
      SetEnumPopup (mpp->complete, mpp->completeAlist, 0);
    }

    StaticPrompt (f, "Technique",
                  wid, popupMenuHeight, programFont, 'l');
    mpp->technique = PopupList (f, TRUE, ChangeTech);
    SetObjectExtra (mpp->technique, mpp, NULL);
    InitEnumPopup (mpp->technique, mpp->techniqueAlist, NULL);
    SetEnumPopup (mpp->technique, mpp->techniqueAlist, 0);
    SetEnumPopup (mpp->technique, mpp->techniqueAlist, tech);

    mpp->expGrp = HiddenGroup (f, -2, 0, NULL);
    SetGroupSpacing (mpp->expGrp, 0, 0);
    /*
    StaticPrompt (mpp->expGrp, "Explanation",
                  0, dialogTextHeight, programFont, 'l');
    */
    mpp->explain = DialogText (mpp->expGrp, "", 10, NULL);
    SetObjectExtra (mpp->explain, mpp, NULL);
    if (showComplete) {
      AlignObjects (ALIGN_LEFT, (HANDLE) mpp->explain,
                    (HANDLE) x, NULL);
      AlignObjects (ALIGN_RIGHT, (HANDLE) mpp->explain,
                    (HANDLE) mpp->complete, NULL);
      AlignObjects (ALIGN_RIGHT, (HANDLE) mpp->explain,
                    (HANDLE) x, NULL);
    } else {
      AlignObjects (ALIGN_RIGHT, (HANDLE) mpp->technique,
                    (HANDLE) mpp->explain, NULL);
      AlignObjects (ALIGN_LEFT, (HANDLE) mpp->technique,
                    (HANDLE) mpp->explain, NULL);
    }
    Hide (mpp->expGrp);
  }

  return (DialoG) p;
}

/************ Bioseq field editors ************/

static ENUM_ALIST(mol_alist)
{" ",               0},             /* Unknown? */
{"DNA",             Seq_mol_dna},   /* 1 */
{"RNA",             Seq_mol_rna},   /* 2 */
{"Protein",         Seq_mol_aa},    /* 3 */
{"Nucleotide",      Seq_mol_na},    /* 4 */
{"Other",           Seq_mol_other}, /* 255 */
END_ENUM_ALIST

static ENUM_ALIST(topology_alist)
{" ",               0},                 /* Unknown? */
{"Linear",          TOPOLOGY_LINEAR},   /* 1 */
{"Circular",        TOPOLOGY_CIRCULAR}, /* 2 */
{"Tandem",          TOPOLOGY_TANDEM},   /* 3 */
{"Other",           255},               /* Other? */
END_ENUM_ALIST

static ENUM_ALIST(strand_alist)
{" ",               Seq_strand_unknown},  /* 0 */
{"Single",          Seq_strand_plus},     /* 1 */
{"Double",          Seq_strand_minus},    /* 2 */
{"Mixed",           Seq_strand_both},     /* 3 */
{"Mixed Rev",       Seq_strand_both_rev}, /* 4 */
{"Other",           Seq_strand_other},    /* 255 */
END_ENUM_ALIST

static Uint1 check_mol (Uint1 mol)
{
  if (mol > Seq_mol_na && mol != Seq_mol_other) return 0;
  return mol;
}

static Uint1 check_topology (Uint1 topology)
{
  if (topology > TOPOLOGY_TANDEM && topology != 255) return 0;
  return topology;
}

static Uint1 check_strand (Uint1 strand)
{
  /* if (strand == Seq_strand_both_rev) return Seq_strand_other; ??? */
  if (strand > Seq_strand_both_rev && strand != Seq_strand_other)
    return Seq_strand_unknown;
  return strand;
}

static Boolean UpdateSeqInstGatherFunc (GatherContextPtr gcp)

{
  BioseqPtr       bsp;
  MolInfoFormPtr  mfp;
  UIEnum          val;
  
  if (gcp == NULL) return TRUE;
  mfp = (MolInfoFormPtr) gcp->userdata;
  if (mfp == NULL) return TRUE;
  if (gcp->thistype == OBJ_BIOSEQ && gcp->thisitem != NULL) {
    bsp = (BioseqPtr) gcp->thisitem;
    if (GetEnumPopup (mfp->molPopup, mol_alist, &val)) {
      bsp->mol = (Uint1) val;
    }
    if (GetEnumPopup (mfp->strandPopup, strand_alist, &val)) {
      bsp->strand = (Uint1) val;
    }
    if (GetEnumPopup (mfp->topologyPopup, topology_alist, &val)) {
      bsp->topology = (Uint1) val;
    }
  }
  return TRUE;
}

static Boolean UpdateSeqInstFlags (GatherContextPtr gcp)

{
  GatherScope     gs;
  MolInfoFormPtr  mfp;
  SeqEntryPtr     sep;

  if (gcp == NULL || gcp->userdata == NULL) return TRUE;
  mfp = (MolInfoFormPtr) gcp->userdata;
  if (gcp->parenttype == OBJ_BIOSEQ || gcp->parenttype == OBJ_BIOSEQSET) {
    sep = SeqMgrGetSeqEntryForData (gcp->parentitem);
    if (sep != NULL) {
      MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
      gs.seglevels = 1;
      gs.scope = sep;
      MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
      gs.ignore[OBJ_BIOSEQ] = FALSE;
      gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
      gs.ignore[OBJ_SEQANNOT] = FALSE;
      gs.ignore[OBJ_SEQDESC] = FALSE;
      GatherEntity (mfp->input_entityID, (Pointer) mfp, UpdateSeqInstGatherFunc, &gs);
    }
  }
  return TRUE;
}

static Boolean CollectSeqInstGatherFunc (GatherContextPtr gcp)

{
  BioseqPtr       bsp;
  MolInfoFormPtr  mfp;
  
  if (gcp == NULL) return TRUE;
  mfp = (MolInfoFormPtr) gcp->userdata;
  if (mfp == NULL) return TRUE;
  if (gcp->thistype == OBJ_BIOSEQ && gcp->thisitem != NULL) {
    bsp = (BioseqPtr) gcp->thisitem;
    SetEnumPopup (mfp->molPopup, mol_alist,
                  (UIEnum) check_mol (bsp->mol));
    SetEnumPopup (mfp->strandPopup, strand_alist,
                  (UIEnum) check_strand (bsp->strand));
    SetEnumPopup (mfp->topologyPopup, topology_alist,
                  (UIEnum) check_topology (bsp->topology));
  }
  return TRUE;
}

static Boolean CollectSeqInstFlags (GatherContextPtr gcp)

{
  GatherScope     gs;
  MolInfoFormPtr  mfp;
  SeqEntryPtr     sep;

  if (gcp == NULL || gcp->userdata == NULL) return TRUE;
  mfp = (MolInfoFormPtr) gcp->userdata;
  if (gcp->parenttype == OBJ_BIOSEQ || gcp->parenttype == OBJ_BIOSEQSET) {
    sep = SeqMgrGetSeqEntryForData (gcp->parentitem);
    if (sep != NULL) {
      MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
      gs.seglevels = 1;
      gs.scope = sep;
      MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
      gs.ignore[OBJ_BIOSEQ] = FALSE;
      gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
      gs.ignore[OBJ_SEQANNOT] = FALSE;
      gs.ignore[OBJ_SEQDESC] = FALSE;
      GatherEntity (mfp->input_entityID, (Pointer) mfp, CollectSeqInstGatherFunc, &gs);
    }
  }
  return TRUE;
}

static void MolInfoDescFormActnProc (ForM f)

{
  MolInfoFormPtr  mfp;

  mfp = (MolInfoFormPtr) GetObjectExtra (f);
  if (mfp != NULL) {
    DescFormReplaceWithoutUpdateProc (f);
    if (mfp->input_itemtype == OBJ_SEQDESC) {
      GatherItem (mfp->input_entityID, mfp->input_itemID, mfp->input_itemtype,
                  (Pointer) mfp, UpdateSeqInstFlags);
    }
    GetRidOfEmptyFeatsDescStrings (mfp->input_entityID, NULL);
    if (GetAppProperty ("InternalNcbiSequin") != NULL) {
      ExtendGeneFeatIfOnMRNA (mfp->input_entityID, NULL);
    }
    ObjMgrSendMsg (OM_MSG_UPDATE, mfp->input_entityID,
                   mfp->input_itemID, mfp->input_itemtype);
  }
}

static void MolInfoFormMessage (ForM f, Int2 mssg)

{
  MolInfoFormPtr  mfp;

  mfp = (MolInfoFormPtr) GetObjectExtra (f);
  if (mfp != NULL) {
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
        if (mfp->appmessage != NULL) {
          mfp->appmessage (f, mssg);
        }
        break;
    }
  }
}

extern ForM CreateMolInfoForm (Int2 left, Int2 top, Int2 width, Int2 height,
                               CharPtr title, Uint1 biomol, Uint1 tech,
                               Boolean nucsOK, Boolean protsOK, FormActnFunc actproc)

{
  ButtoN             b;
  GrouP              c;
  GrouP              g;
  MolInfoFormPtr     mfp;
  PrompT             p;
  GrouP              q;
  StdEditorProcsPtr  sepp;
  WindoW             w;

  w = NULL;
  mfp = (MolInfoFormPtr) MemNew (sizeof (MolInfoForm));
  if (mfp != NULL) {
    w = FixedWindow (left, top, width, height, title, StdCloseWindowProc);
    SetObjectExtra (w, mfp, StdDescFormCleanupProc);
    mfp->form = (ForM) w;
    mfp->actproc = actproc;
    mfp->formmessage = MolInfoFormMessage;

#ifndef WIN_MAC
    CreateStdEditorFormMenus (w);
#endif
    sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
    if (sepp != NULL) {
      SetActivate (w, sepp->activateForm);
      mfp->appmessage = sepp->handleMessages;
    }

    g = HiddenGroup (w, -1, 0, NULL);
    mfp->data = CreateMolInfoDialog (g, NULL, biomol, tech, TRUE, nucsOK, protsOK, FALSE);

    p = StaticPrompt (w, "Biological characteristics of sequences:",
                      0, stdLineHeight, programFont, 'c');

    q = HiddenGroup (w, 6, 0, NULL);
    StaticPrompt (q, "Class", 0, popupMenuHeight, programFont, 'c');
    mfp->molPopup = PopupList (q, TRUE, NULL);
    InitEnumPopup (mfp->molPopup, mol_alist, NULL);
    StaticPrompt (q, "Topology", 0, popupMenuHeight, programFont, 'c');
    mfp->topologyPopup = PopupList (q, TRUE, NULL);
    InitEnumPopup (mfp->topologyPopup, topology_alist, NULL);
    StaticPrompt (q, "Strand", 0, popupMenuHeight, programFont, 'c');
    mfp->strandPopup = PopupList (q, TRUE, NULL);
    InitEnumPopup (mfp->strandPopup, strand_alist, NULL);

    SetGroupSpacing (q, 15, 2);
    c = HiddenGroup (w, 2, 0, NULL);
    b = PushButton (c, "Accept", StdAcceptFormButtonProc);
    SetObjectExtra (b, mfp, NULL);
    PushButton (c, "Cancel", StdCancelButtonProc);
    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) p, (HANDLE) q, (HANDLE) c, NULL);
    RealizeWindow (w);
  }
  return (ForM) w;
}

extern Int2 LIBCALLBACK MolInfoGenFunc (Pointer data)

{
  MolInfoFormPtr    mfp;
  OMProcControlPtr  ompcp;
  OMUserDataPtr     omudp;
  ValNodePtr        sdp;
  WindoW            w;

  ompcp = (OMProcControlPtr) data;
  sdp = NULL;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  switch (ompcp->input_itemtype) {
    case OBJ_SEQDESC :
      sdp = (ValNodePtr) ompcp->input_data;
      if (sdp != NULL && sdp->choice != Seq_descr_molinfo) {
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
    mfp = (MolInfoFormPtr) omudp->userdata.ptrvalue;
    if (mfp != NULL) {
      Select (mfp->form);
    }
    return OM_MSG_RET_DONE;
  }
  w = (WindoW) CreateMolInfoForm (-50, -33, -10, -10,
                                  "Molecule Information", 0, 0, TRUE, TRUE,
                                  MolInfoDescFormActnProc);
  mfp = (MolInfoFormPtr) GetObjectExtra (w);
  if (mfp != NULL) {
    mfp->input_entityID = ompcp->input_entityID;
    mfp->input_itemID = ompcp->input_itemID;
    mfp->input_itemtype = ompcp->input_itemtype;
    mfp->this_itemtype = OBJ_SEQDESC;
    mfp->this_subtype = Seq_descr_molinfo;
    mfp->procid = ompcp->proc->procid;
    mfp->proctype = ompcp->proc->proctype;
    mfp->userkey = OMGetNextUserKey ();
    omudp = ObjMgrAddUserData (ompcp->input_entityID, ompcp->proc->procid,
	                           OMPROC_EDIT, mfp->userkey);
    if (omudp != NULL) {
      omudp->userdata.ptrvalue = (Pointer) mfp;
      omudp->messagefunc = StdVibrantEditorMsgFunc;
    }
    SendMessageToDialog (mfp->data, VIB_MSG_INIT);
    if (sdp != NULL) {
      PointerToDialog (mfp->data, (Pointer) sdp->data.ptrvalue);
      if (mfp->input_itemtype == OBJ_SEQDESC) {
        GatherItem (mfp->input_entityID, mfp->input_itemID, mfp->input_itemtype,
                    (Pointer) mfp, CollectSeqInstFlags);
      }
      SetClosestParentIfDuplicating ((BaseFormPtr) mfp);
    }
  }
  Show (w);
  Select (w);
  return OM_MSG_RET_DONE;
}


typedef struct gbblockpage {
  DIALOG_MESSAGE_BLOCK
  TexT          source;
  TexT          origin;
  TexT          date;
  TexT          div;
  TexT          taxonomy;
  DialoG        kywds;
  DialoG        xaccns;
  DialoG        entryDate;
} GenBankPage, PNTR GenBankPagePtr;

typedef struct gbblockform {
  DESCRIPTOR_FORM_BLOCK
  DialoG        gbppxaccns;
  ButtoN        xaccnstohistory;
} GenBankForm, PNTR GenBankFormPtr;

static void GBBlockPtrToGenBankPage (DialoG d, Pointer data)

{
  GBBlockPtr      gbp;
  GenBankPagePtr  gpp;

  gpp = (GenBankPagePtr) GetObjectExtra (d);
  gbp = (GBBlockPtr) data;
  if (gpp != NULL) {
    SetTitle (gpp->source, gbp->source);
    SetTitle (gpp->origin, gbp->origin);
    SetTitle (gpp->date, gbp->date);
    SetTitle (gpp->div, gbp->div);
    SetTitle (gpp->taxonomy, gbp->taxonomy);
    PointerToDialog (gpp->kywds, gbp->keywords);
    PointerToDialog (gpp->xaccns, gbp->extra_accessions);
    PointerToDialog (gpp->entryDate, gbp->entry_date);
  }
}

static Pointer GenBankPageToGBBlockPtr (DialoG d)

{
  GBBlockPtr      gbp;
  GenBankPagePtr  gpp;

  gbp = NULL;
  gpp = (GenBankPagePtr) GetObjectExtra (d);
  if (gpp != NULL) {
    gbp = GBBlockNew ();
    if (gbp != NULL) {
      gbp->source = SaveStringFromText (gpp->source);
      gbp->origin = SaveStringFromText (gpp->origin);
      gbp->date = SaveStringFromText (gpp->date);
      gbp->div = SaveStringFromText (gpp->div);
      gbp->taxonomy = SaveStringFromTextAndStripNewlines (gpp->taxonomy);
      gbp->keywords = DialogToPointer (gpp->kywds);
      gbp->extra_accessions = DialogToPointer (gpp->xaccns);
      gbp->entry_date = DialogToPointer (gpp->entryDate);
    }
  }
  return (Pointer) gbp;
}

static void DeleteKeywordProc (ButtoN b)

{
  GenBankPagePtr  gpp;

  gpp = (GenBankPagePtr) GetObjectExtra (b);
  if (gpp != NULL) {
    PointerToDialog (gpp->kywds, NULL);
  }
}

static DialoG CreateGenBankDialog (GrouP h, CharPtr title, ValNodePtr sdp, GenBankFormPtr gfp)

{
  ButtoN          b;
  GrouP           f1, f2, f3, f4;
  GBBlockPtr      gbp;
  Boolean         genome;
  GenBankPagePtr  gpp;
  Boolean         internal;
  GrouP           m;
  GrouP           p;
  GrouP           s;

  p = HiddenGroup (h, 1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  gpp = (GenBankPagePtr) MemNew (sizeof (GenBankPage));
  if (gpp) {

    SetObjectExtra (p, gpp, StdCleanupExtraProc);
    gpp->dialog = (DialoG) p;
    gpp->todialog = GBBlockPtrToGenBankPage;
    gpp->fromdialog = GenBankPageToGBBlockPtr;
    gpp->testdialog = NULL;

    if (title != NULL && title [0] != '\0') {
      s = NormalGroup (p, 0, -2, title, systemFont, NULL);
    } else {
      s = HiddenGroup (p, 0, -2, NULL);
    }
    m = HiddenGroup (s, -1, 0, NULL);
    /*
    SetGroupSpacing (m, 10, 10);
    */

    gbp = NULL;
    if (sdp != NULL && sdp->choice == Seq_descr_genbank && sdp->data.ptrvalue != NULL) {
      gbp = (GBBlockPtr) sdp->data.ptrvalue;
    }
    internal = (Boolean) (GetAppProperty ("InternalNcbiSequin") != NULL);
    genome = (Boolean) (GetAppProperty ("GenomeCenterSequin") != NULL);

    f1 = HiddenGroup (m, 1, 0, NULL);

    if ((! genome) || (gbp != NULL && gbp->div != NULL)) {
      StaticPrompt (f1, "Division", 0, 0, programFont, 'c');
      gpp->div = DialogText (f1, "", 15, NULL);
    }

    if (internal || (gbp != NULL && gbp->origin != NULL)) {
      StaticPrompt (f1, "Origin", 0, 0, programFont, 'c');
      gpp->origin = DialogText (f1, "", 15, NULL);
    }

    if (internal || (gbp != NULL && gbp->date != NULL)) {
      StaticPrompt (f1, "(Old Date)", 0, 0, programFont, 'c');
      gpp->date = DialogText (f1, "", 15, NULL);
    }

    if (internal || (gbp != NULL && gbp->source != NULL)) {
      StaticPrompt (f1, "Source", 0, 0, programFont, 'c');
      gpp->source = DialogText (f1, "", 15, NULL);
    }

    if (internal || (gbp != NULL && gbp->taxonomy != NULL)) {
      StaticPrompt (f1, "Taxonomy", 0, 0, programFont, 'c');
      gpp->taxonomy = ScrollText (f1, 15, 5, programFont, TRUE, NULL);
    }

    f3 = HiddenGroup (m, 0, 3, NULL);
    if (internal || genome || (gbp != NULL && gbp->extra_accessions != NULL)) {
      StaticPrompt (f3, "Secondary Accessions", 0, 0, programFont, 'c');
      gpp->xaccns = CreateVisibleStringDialog (f3, 3, -1, 15);
      if (internal || genome) {
        gfp->gbppxaccns = gpp->xaccns;
        gfp->xaccnstohistory = CheckBox (f3, "Copy to Bioseq-history.replaces", NULL);
        if (genome) {
          SetStatus (gfp->xaccnstohistory, TRUE);
          Hide (gfp->xaccnstohistory);
        }
      }
    }

    b = NULL;
    f2 = HiddenGroup (m, 0, 2, NULL);
    if (internal || (gbp != NULL && gbp->keywords != NULL)) {
      StaticPrompt (f2, "Keywords", 0, 0, programFont, 'c');
      gpp->kywds = CreateVisibleStringDialog (f2, 3, -1, 15);
      if (internal) {
        b = PushButton (m, "Delete All Keywords", DeleteKeywordProc);
        SetObjectExtra (b, gpp, NULL);
      }
    }

    f4 = HiddenGroup (m, 2, 0, NULL);
    if (internal || (gbp != NULL && gbp->entry_date != NULL)) {
      gpp->entryDate = CreateDateDialog (f4, "");
    }

    AlignObjects (ALIGN_CENTER, (HANDLE) f1, (HANDLE) f2, (HANDLE) f3,
                  (HANDLE) f4, (HANDLE) b, NULL);
  }

  return (DialoG) p;
}

static void GenBankFormMessage (ForM f, Int2 mssg)

{
  GenBankFormPtr  gfp;

  gfp = (GenBankFormPtr) GetObjectExtra (f);
  if (gfp != NULL) {
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
        if (gfp->appmessage != NULL) {
          gfp->appmessage (f, mssg);
        }
        break;
    }
  }
}

extern ForM CreateGenBankForm (Int2 left, Int2 top, CharPtr title,
                               ValNodePtr sdp, SeqEntryPtr sep,
                               FormActnFunc actproc)

{
  ButtoN             b;
  GrouP              c;
  GrouP              g;
  GenBankFormPtr     gfp;
  StdEditorProcsPtr  sepp;
  WindoW             w;

  w = NULL;
  gfp = (GenBankFormPtr) MemNew (sizeof (GenBankForm));
  if (gfp != NULL) {
    w = FixedWindow (left, top, -10, -10, title, StdCloseWindowProc);
    SetObjectExtra (w, gfp, StdDescFormCleanupProc);
    gfp->form = (ForM) w;
    gfp->actproc = actproc;
    gfp->formmessage = GenBankFormMessage;

#ifndef WIN_MAC
    CreateStdEditorFormMenus (w);
#endif

    sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
    if (sepp != NULL) {
      SetActivate (w, sepp->activateForm);
      gfp->appmessage = sepp->handleMessages;
    }

    g = HiddenGroup (w, -1, 0, NULL);
    gfp->xaccnstohistory = NULL;
    gfp->gbppxaccns = NULL;
    gfp->data = CreateGenBankDialog (g, NULL, sdp, gfp);

    c = HiddenGroup (w, 2, 0, NULL);
    b = PushButton (c, "Accept", StdAcceptFormButtonProc);
    SetObjectExtra (b, gfp, NULL);
    PushButton (c, "Cancel", StdCancelButtonProc);
    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
    RealizeWindow (w);
  }
  return (ForM) w;
}

static Boolean GetLowestStackSeqEntryForGBB (GatherContextPtr gcp)

{
  BaseFormPtr  bfp;
  Int2         i;

  if (gcp == NULL) return TRUE;
  bfp = (BaseFormPtr) gcp->userdata;
  if (bfp == NULL) return TRUE;
  if (gcp->gatherstack != NULL && gcp->numstack > 0) {
    for (i = 0; i < gcp->numstack; i++) {
      if (gcp->gatherstack [i].itemtype == OBJ_BIOSEQ ||
          gcp->gatherstack [i].itemtype == OBJ_BIOSEQSET) {
        bfp->input_itemID = gcp->gatherstack [i].itemID;
        bfp->input_itemtype = gcp->gatherstack [i].itemtype;
      }
    }
  }
  return FALSE;
}

static void GenBankFormActnProc (ForM f)

{
  BioseqPtr       bsp;
  Boolean         changehist;
  GenBankFormPtr  gfp;
  SeqHistPtr      hist;
  SeqIdPtr        sip;
  TextSeqIdPtr    tsip;
  ValNodePtr      tmp;
  ValNodePtr      vnp;

  gfp = (GenBankFormPtr) GetObjectExtra (f);
  if (gfp != NULL) {
    vnp = NULL;
    changehist = (Boolean) (gfp->xaccnstohistory != NULL &&
                            GetStatus (gfp->xaccnstohistory));
    if (changehist) {
      vnp = DialogToPointer (gfp->gbppxaccns);
    }
    if (DescFormReplaceWithoutUpdateProc (f)) {
      GetRidOfEmptyFeatsDescStrings (gfp->input_entityID, NULL);
      if (GetAppProperty ("InternalNcbiSequin") != NULL) {
        ExtendGeneFeatIfOnMRNA (gfp->input_entityID, NULL);
      }
      if (changehist) {
        GatherItem (gfp->input_entityID, gfp->input_itemID, gfp->input_itemtype,
                    (Pointer) gfp, GetLowestStackSeqEntryForGBB);
        bsp = GetBioseqGivenIDs (gfp->input_entityID, gfp->input_itemID, gfp->input_itemtype);
        if (bsp != NULL) {
          hist = bsp->hist;
          if (hist != NULL) {
            hist->replace_ids = SeqIdSetFree (hist->replace_ids);
          } else {
            hist = SeqHistNew ();
            bsp->hist = hist;
          }
          if (hist != NULL && vnp != NULL) {
            for (tmp = vnp; tmp != NULL; tmp = tmp->next) {
              sip = ValNodeNew (hist->replace_ids);
              if (hist->replace_ids == NULL) {
                hist->replace_ids = sip;
              }
              if (sip != NULL) {
                tsip = TextSeqIdNew ();
                sip->choice = SEQID_GENBANK;
                sip->data.ptrvalue = (Pointer) tsip;
                if (tsip != NULL) {
                  tsip->accession = StringSave (tmp->data.ptrvalue);
                }
              }
            }
          }
          if (hist != NULL) {
            if (hist->assembly == NULL && hist->replace_date == NULL &&
                hist->replace_ids == NULL && hist->replaced_by_date == NULL &&
                hist->replaced_by_ids == NULL && hist->deleted_date == NULL &&
                (! hist->deleted)) {
              bsp->hist = SeqHistFree (bsp->hist);
            }
          }
        }
        ValNodeFreeData (vnp);
      }
      ObjMgrSendMsg (OM_MSG_UPDATE, gfp->input_entityID,
                     gfp->input_itemID, gfp->input_itemtype);
    }
  }
}

extern Int2 LIBCALLBACK GenBankGenFunc (Pointer data)

{
  GenBankFormPtr    gfp;
  OMProcControlPtr  ompcp;
  OMUserDataPtr     omudp;
  ValNodePtr        sdp;
  SeqEntryPtr       sep;
  WindoW            w;

  ompcp = (OMProcControlPtr) data;
  sdp = NULL;
  sep = NULL;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  switch (ompcp->input_itemtype) {
    case OBJ_SEQDESC :
      sdp = (ValNodePtr) ompcp->input_data;
      if (sdp != NULL && sdp->choice != Seq_descr_genbank) {
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
    gfp = (GenBankFormPtr) omudp->userdata.ptrvalue;
    if (gfp != NULL) {
      Select (gfp->form);
    }
    return OM_MSG_RET_DONE;
  }
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  w = (WindoW) CreateGenBankForm (-50, -33, "GenBank Block", sdp, sep,
                                  GenBankFormActnProc);
  gfp = (GenBankFormPtr) GetObjectExtra (w);
  if (gfp != NULL) {
    gfp->input_entityID = ompcp->input_entityID;
    gfp->input_itemID = ompcp->input_itemID;
    gfp->input_itemtype = ompcp->input_itemtype;
    gfp->this_itemtype = OBJ_SEQDESC;
    gfp->this_subtype = Seq_descr_genbank;
    gfp->procid = ompcp->proc->procid;
    gfp->proctype = ompcp->proc->proctype;
    gfp->userkey = OMGetNextUserKey ();
    omudp = ObjMgrAddUserData (ompcp->input_entityID, ompcp->proc->procid,
	                           OMPROC_EDIT, gfp->userkey);
    if (omudp != NULL) {
      omudp->userdata.ptrvalue = (Pointer) gfp;
      omudp->messagefunc = StdVibrantEditorMsgFunc;
    }
    SendMessageToDialog (gfp->data, VIB_MSG_INIT);
    if (sdp != NULL) {
      PointerToDialog (gfp->data, (Pointer) sdp->data.ptrvalue);
      SetClosestParentIfDuplicating ((BaseFormPtr) gfp);
    }
  }
  Show (w);
  Select (w);
  return OM_MSG_RET_DONE;
}


typedef struct visstrpage {
  DIALOG_MESSAGE_BLOCK
  TexT          title;
} VisStrPage, PNTR VisStrPagePtr;

typedef struct visstrform {
  DESCRIPTOR_FORM_BLOCK
} VisStrForm, PNTR VisStrFormPtr;

static void CharPtrToVisStrPage (DialoG d, Pointer data)

{
  CharPtr        title;
  VisStrPagePtr  vpp;

  vpp = (VisStrPagePtr) GetObjectExtra (d);
  title = (CharPtr) data;
  if (vpp != NULL) {
    SetTitle (vpp->title, title);
  }
}

static Pointer VisStrPageToCharPtr (DialoG d)

{
  CharPtr        title;
  VisStrPagePtr  vpp;

  title = NULL;
  vpp = (VisStrPagePtr) GetObjectExtra (d);
  if (vpp != NULL) {
    title = SaveStringFromTextAndStripNewlines (vpp->title);
  }
  return (Pointer) title;
}

static void VisStrDialogMessage (DialoG d, Int2 mssg)

{
  VisStrPagePtr  vpp;

  vpp = (VisStrPagePtr) GetObjectExtra (d);
  if (vpp != NULL) {
    switch (mssg) {
      case VIB_MSG_ENTER :
        Select (vpp->title);
        break;
      default :
        break;
    }
  }
}

static DialoG CreateVisStrDialog (GrouP h, CharPtr title, Uint2 subtype)

{
  GrouP          f;
  CharPtr        label;
  GrouP          m;
  GrouP          p;
  GrouP          s;
  VisStrPagePtr  vpp;

  p = HiddenGroup (h, 1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  vpp = (VisStrPagePtr) MemNew (sizeof (VisStrPage));
  if (vpp) {

    SetObjectExtra (p, vpp, StdCleanupExtraProc);
    vpp->dialog = (DialoG) p;
    vpp->todialog = CharPtrToVisStrPage;
    vpp->fromdialog = VisStrPageToCharPtr;
    vpp->testdialog = NULL;
    vpp->dialogmessage = VisStrDialogMessage;

    if (title != NULL && title [0] != '\0') {
      s = NormalGroup (p, 0, -2, title, systemFont, NULL);
    } else {
      s = HiddenGroup (p, 0, -2, NULL);
    }
    m = HiddenGroup (s, -1, 0, NULL);
    /*
    SetGroupSpacing (m, 10, 10);
    */

    f = HiddenGroup (m, 0, 2, NULL);
    if (subtype == Seq_descr_title) {
      label = "Title";
    } else if (subtype == Seq_descr_comment) {
      label = "Comment";
    } else if (subtype == Seq_descr_name) {
      label = "Name";
    } else if (subtype == Seq_descr_region) {
      label = "Region";
    } else {
      label = "Title";
    }
    StaticPrompt (f, label, 0, 0, programFont, 'c');
    vpp->title = ScrollText (f, 25, 15, programFont, TRUE, NULL);
  }

  return (DialoG) p;
}

static void VisStrFormMessage (ForM f, Int2 mssg)

{
  VisStrFormPtr  vfp;

  vfp = (VisStrFormPtr) GetObjectExtra (f);
  if (vfp != NULL) {
    switch (mssg) {
      case VIB_MSG_ENTER :
        SendMessageToDialog (vfp->data, VIB_MSG_ENTER);
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
        StdDeleteTextProc (NULL);
        break;
      default :
        if (vfp->appmessage != NULL) {
          vfp->appmessage (f, mssg);
        }
        break;
    }
  }
}

extern ForM CreateVisStrForm (Int2 left, Int2 top, CharPtr title,
                              Uint2 subtype, FormActnFunc actproc)

{
  ButtoN             b;
  GrouP              c;
  GrouP              g;
  StdEditorProcsPtr  sepp;
  VisStrFormPtr      vfp;
  WindoW             w;

  w = NULL;
  vfp = (VisStrFormPtr) MemNew (sizeof (VisStrForm));
  if (vfp != NULL) {
    w = FixedWindow (left, top, -10, -10, title, StdCloseWindowProc);
    SetObjectExtra (w, vfp, StdDescFormCleanupProc);
    vfp->form = (ForM) w;
    vfp->actproc = actproc;
    vfp->formmessage = VisStrFormMessage;

#ifndef WIN_MAC
    CreateStdEditorFormMenus (w);
#endif

    sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
    if (sepp != NULL) {
      SetActivate (w, sepp->activateForm);
      vfp->appmessage = sepp->handleMessages;
    }

    g = HiddenGroup (w, -1, 0, NULL);
    vfp->data = CreateVisStrDialog (g, NULL, subtype);

    c = HiddenGroup (w, 2, 0, NULL);
    b = PushButton (c, "Accept", StdAcceptFormButtonProc);
    SetObjectExtra (b, vfp, NULL);
    PushButton (c, "Cancel", StdCancelButtonProc);
    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
    RealizeWindow (w);
  }
  return (ForM) w;
}

static void VisStrDescFormActnProc (ForM f)

{
  DescriptorFormPtr  dfp;
  OMProcControl      ompc;
  CharPtr            str;

  dfp = (DescriptorFormPtr) GetObjectExtra (f);
  if (dfp != NULL) {
    str = DialogToPointer (dfp->data);
    if (str == NULL || StringHasNoText (str)) {
      MemSet ((Pointer) &ompc, 0, sizeof (OMProcControl));
      ompc.do_not_reload_from_cache = TRUE;
      ompc.input_entityID = dfp->input_entityID;
      ompc.input_itemID = dfp->input_itemID;
      ompc.input_itemtype = dfp->input_itemtype;
      if (! DetachDataForProc (&ompc, FALSE)) {
        Message (MSG_ERROR, "DetachDataForProc failed");
      }
      GetRidOfEmptyFeatsDescStrings (dfp->input_entityID, NULL);
      ObjMgrSetDirtyFlag (dfp->input_entityID, TRUE);
      ObjMgrSendMsg (OM_MSG_UPDATE, dfp->input_entityID,
                     dfp->input_itemID, dfp->input_itemtype);
      MemFree (str);
      return;
    }
    MemFree (str);
    if (DescFormReplaceWithoutUpdateProc (f)) {
      GetRidOfEmptyFeatsDescStrings (dfp->input_entityID, NULL);
      if (GetAppProperty ("InternalNcbiSequin") != NULL) {
        ExtendGeneFeatIfOnMRNA (dfp->input_entityID, NULL);
      }
      ObjMgrSetDirtyFlag (dfp->input_entityID, TRUE);
      ObjMgrSendMsg (OM_MSG_UPDATE, dfp->input_entityID,
                     dfp->input_itemID, dfp->input_itemtype);
    }
  }
}

extern Int2 LIBCALLBACK VisStrGenFunc (Pointer data)

{
  OMProcControlPtr  ompcp;
  OMUserDataPtr     omudp;
  ObjMgrProcPtr     proc;
  ValNodePtr        sdp;
  Uint2             subtype;
  VisStrFormPtr     vfp;
  WindoW            w;

  ompcp = (OMProcControlPtr) data;
  sdp = NULL;
  subtype = 0;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  proc = ompcp->proc;
  switch (ompcp->input_itemtype) {
    case OBJ_SEQDESC :
      sdp = (ValNodePtr) ompcp->input_data;
      if (sdp != NULL &&
         (sdp->choice != Seq_descr_title &&
          sdp->choice != Seq_descr_comment &&
          sdp->choice != Seq_descr_name &&
          sdp->choice != Seq_descr_region)) {
        return OM_MSG_RET_ERROR;
      }
      subtype = sdp->choice;
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
    vfp = (VisStrFormPtr) omudp->userdata.ptrvalue;
    if (vfp != NULL) {
      Select (vfp->form);
    }
    return OM_MSG_RET_DONE;
  }
  if (sdp == NULL) {
    subtype = proc->subinputtype;
  }
  if (subtype == Seq_descr_title) {
      w = (WindoW) CreateVisStrForm (-50, -33, "Title",
                                     subtype, VisStrDescFormActnProc);
  } else if (subtype == Seq_descr_comment) {
      w = (WindoW) CreateVisStrForm (-50, -33, "Comment",
                                     subtype, VisStrDescFormActnProc);
  } else if (subtype == Seq_descr_name) {
      w = (WindoW) CreateVisStrForm (-50, -33, "Name",
                                     subtype, VisStrDescFormActnProc);
  } else if (subtype == Seq_descr_region) {
      w = (WindoW) CreateVisStrForm (-50, -33, "Region",
                                     subtype, VisStrDescFormActnProc);
  } else {
    return OM_MSG_RET_ERROR;
  }
  vfp = (VisStrFormPtr) GetObjectExtra (w);
  if (vfp != NULL) {
    vfp->input_entityID = ompcp->input_entityID;
    vfp->input_itemID = ompcp->input_itemID;
    vfp->input_itemtype = ompcp->input_itemtype;
    vfp->this_itemtype = OBJ_SEQDESC;
    vfp->this_subtype = subtype;
    vfp->procid = ompcp->proc->procid;
    vfp->proctype = ompcp->proc->proctype;
    vfp->userkey = OMGetNextUserKey ();
    omudp = ObjMgrAddUserData (ompcp->input_entityID, ompcp->proc->procid,
	                           OMPROC_EDIT, vfp->userkey);
    if (omudp != NULL) {
      omudp->userdata.ptrvalue = (Pointer) vfp;
      omudp->messagefunc = StdVibrantEditorMsgFunc;
    }
    SendMessageToDialog (vfp->data, VIB_MSG_INIT);
    if (sdp != NULL) {
      PointerToDialog (vfp->data, (Pointer) sdp->data.ptrvalue);
      SetClosestParentIfDuplicating ((BaseFormPtr) vfp);
    }
  }
  Show (w);
  Select (w);
  return OM_MSG_RET_DONE;
}


typedef struct dateform {
  DESCRIPTOR_FORM_BLOCK
} DateForm, PNTR DateFormPtr;

static void DateFormMessage (ForM f, Int2 mssg)

{
  DateFormPtr  dfp;

  dfp = (DateFormPtr) GetObjectExtra (f);
  if (dfp != NULL) {
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
        if (dfp->appmessage != NULL) {
          dfp->appmessage (f, mssg);
        }
        break;
    }
  }
}

extern ForM CreateDateForm (Int2 left, Int2 top, CharPtr title,
                            FormActnFunc actproc)

{
  ButtoN             b;
  GrouP              c;
  DateFormPtr        dfp;
  GrouP              g;
  StdEditorProcsPtr  sepp;
  WindoW             w;

  w = NULL;
  dfp = (DateFormPtr) MemNew (sizeof (DateForm));
  if (dfp != NULL) {
    w = FixedWindow (left, top, -10, -10, title, StdCloseWindowProc);
    SetObjectExtra (w, dfp, StdDescFormCleanupProc);
    dfp->form = (ForM) w;
    dfp->actproc = actproc;
    dfp->formmessage = DateFormMessage;

#ifndef WIN_MAC
    CreateStdEditorFormMenus (w);
#endif

    sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
    if (sepp != NULL) {
      SetActivate (w, sepp->activateForm);
      dfp->appmessage = sepp->handleMessages;
    }

    g = HiddenGroup (w, -1, 0, NULL);
    dfp->data = CreateDateDialog (g, NULL);

    c = HiddenGroup (w, 2, 0, NULL);
    b = PushButton (c, "Accept", StdAcceptFormButtonProc);
    SetObjectExtra (b, dfp, NULL);
    PushButton (c, "Cancel", StdCancelButtonProc);
    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
    RealizeWindow (w);
  }
  return (ForM) w;
}

extern Int2 LIBCALLBACK DateGenFunc (Pointer data)

{
  DateFormPtr       dfp;
  Uint2             itemtype;
  OMProcControlPtr  ompcp;
  OMUserDataPtr     omudp;
  ObjMgrProcPtr     proc;
  ValNodePtr        sdp;
  CharPtr           title;
  Uint2             subtype;
  WindoW            w;

  ompcp = (OMProcControlPtr) data;
  sdp = NULL;
  itemtype = 0;
  subtype = 0;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  title = "Date";
  proc = ompcp->proc;
  switch (ompcp->input_itemtype) {
    case OBJ_SEQDESC :
      sdp = (ValNodePtr) ompcp->input_data;
      if (sdp != NULL &&
          (sdp->choice != Seq_descr_create_date &&
          sdp->choice != Seq_descr_update_date)) {
        return OM_MSG_RET_ERROR;
      }
      itemtype = OBJ_SEQDESC;
      subtype = sdp->choice;
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
    dfp = (DateFormPtr) omudp->userdata.ptrvalue;
    if (dfp != NULL) {
      Select (dfp->form);
    }
    return OM_MSG_RET_DONE;
  }
  if (sdp != NULL) {
    if (sdp->choice == Seq_descr_create_date) {
      title = "Create Date";
    } else {
      title = "Update Date";
    }
  } else {
    itemtype = proc->inputtype;
    subtype = proc->subinputtype;
    if (itemtype == OBJ_SEQDESC && subtype == Seq_descr_create_date) {
      title = "Create Date";
    } else if (itemtype == OBJ_SEQDESC && subtype == Seq_descr_update_date) {
      title = "Update Date";
    } else {
      return OM_MSG_RET_ERROR;
    }
  }
  w = (WindoW) CreateDateForm (-50, -33, title, StdDescFormActnProc);
  dfp = (DateFormPtr) GetObjectExtra (w);
  if (dfp != NULL) {
    dfp->input_entityID = ompcp->input_entityID;
    dfp->input_itemID = ompcp->input_itemID;
    dfp->input_itemtype = ompcp->input_itemtype;
    dfp->this_itemtype = itemtype;
    dfp->this_subtype = subtype;
    dfp->procid = ompcp->proc->procid;
    dfp->proctype = ompcp->proc->proctype;
    dfp->userkey = OMGetNextUserKey ();
    omudp = ObjMgrAddUserData (ompcp->input_entityID, ompcp->proc->procid,
	                           OMPROC_EDIT, dfp->userkey);
    if (omudp != NULL) {
      omudp->userdata.ptrvalue = (Pointer) dfp;
      omudp->messagefunc = StdVibrantEditorMsgFunc;
    }
    SendMessageToDialog (dfp->data, VIB_MSG_INIT);
    if (sdp != NULL) {
      PointerToDialog (dfp->data, (Pointer) sdp->data.ptrvalue);
      SetClosestParentIfDuplicating ((BaseFormPtr) dfp);
    }
  }
  Show (w);
  Select (w);
  return OM_MSG_RET_DONE;
}











