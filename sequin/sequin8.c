/*   sequin8.c
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
* File Name:  sequin8.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   2/3/98
*
* $Revision: 6.128 $
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
#include <objsub.h>
#include <valid.h>
#include <cdrgn.h>
#include <suggslp.h>
#include <urkptpf.h>
#include <entrez.h>
#include <accentr.h>
#include <urlquery.h>
#include <toasn3.h>
#include <subutil.h>
#include <explore.h>
#include <medarch.h>
#include <medutil.h>
#include <vecscnapi.h>
#include <qblastapi.h>

typedef struct evidenceformdata {
  FEATURE_FORM_BLOCK

  LisT           objlist;
  TexT           findthis;
  Uint2          itemtype;
  Uint2          subtype;
  PopuP          evdence;
  Uint2          exp_ev;
  ObjMgrPtr      omp;
  ObjMgrTypePtr  omtp;
  ValNodePtr     head;
  Boolean        stringfound;
  Char           findStr [128];
} EvidenceFormData, PNTR EvidenceFormPtr;

static void LIBCALLBACK AsnWriteEvidenceCallBack (AsnExpOptStructPtr pAEOS)

{
  EvidenceFormPtr  efp;
  CharPtr          pchFind;
  CharPtr          pchSource;

  efp = (EvidenceFormPtr) pAEOS->data;
  if (ISA_STRINGTYPE (AsnFindBaseIsa (pAEOS->atp))) {
	pchSource = (CharPtr) pAEOS->dvp->ptrvalue;
	pchFind = efp->findStr;
	if (StringSearch (pchSource, pchFind) != NULL) {
	  efp->stringfound = TRUE;
	}
  }
}

static Boolean EvidenceHasSubstring (ObjMgrTypePtr omtp, AsnIoPtr aip, Pointer ptr, EvidenceFormPtr efp)

{
  efp->stringfound = FALSE;
  (omtp->asnwrite) (ptr, aip, NULL);
  return efp->stringfound;
}

static void EvidenceCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  AsnExpOptPtr     aeop;
  AsnIoPtr         aip;
  BioseqPtr        bsp;
  BioseqSetPtr     bssp;
  EvidenceFormPtr  efp;
  Boolean          notext;
  ObjMgrTypePtr    omtp;
  SeqAnnotPtr      sap;
  SeqFeatPtr       sfp;
  Uint2            subtype;

  if (mydata == NULL) return;
  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  efp = (EvidenceFormPtr) mydata;
  if (efp == NULL) return;
  omtp = efp->omtp;
  if (omtp == NULL || omtp->subtypefunc == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
  } else return;
  GetTitle (efp->findthis, efp->findStr, sizeof (efp->findStr) - 1);
  notext = StringHasNoText (efp->findStr);
  aip = AsnIoNullOpen ();
  aeop = AsnExpOptNew (aip, NULL, NULL, AsnWriteEvidenceCallBack);
  if (aeop != NULL) {
    aeop->user_data = (Pointer) efp;
  }
  while (sap != NULL) {
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      while (sfp != NULL) {
        subtype = (*(omtp->subtypefunc)) ((Pointer) sfp);
        if (efp->subtype == 0 || subtype == efp->subtype ||
           (efp->subtype == FEATDEF_IMP &&
            subtype >= FEATDEF_allele && subtype <= FEATDEF_site_ref)) {
          if (notext || EvidenceHasSubstring (omtp, aip, (Pointer) sfp, efp)) {
            sfp->exp_ev = (Uint1) efp->exp_ev;
          }
        }
        sfp = sfp->next;
      }
    }
    sap = sap->next;
  }
  AsnIoClose (aip);
}

static void DoEvidence (ButtoN b)

{
  EvidenceFormPtr  efp;
  SeqEntryPtr      sep;
  Uint2            subtype;
  Int2             val;
  ValNodePtr       vnp;

  efp = GetObjectExtra (b);
  if (efp == NULL) return;
  sep = GetTopSeqEntryForEntityID (efp->input_entityID);
  if (sep == NULL) return;
  Hide (efp->form);
  WatchCursor ();
  Update ();
  efp->itemtype = OBJ_SEQFEAT;
  subtype = 0;
  vnp = NULL;
  val = GetValue (efp->objlist);
  if (val > 0) {
    vnp = efp->head;
    while (vnp != NULL && val > 1) {
      val--;
      vnp = vnp->next;
    }
  }
  if (vnp != NULL) {
    subtype = vnp->choice;
  }
  efp->omp = ObjMgrGet ();
  efp->omtp = NULL;
  if (efp->omp != NULL) {
    efp->omtp = ObjMgrTypeFind (efp->omp, efp->itemtype, NULL, NULL);
  }
  efp->subtype = subtype;
  efp->exp_ev = GetValue (efp->evdence) - 1;
  if (efp->itemtype != 0 && efp->omtp != NULL) {
    SeqEntryExplore (sep, (Pointer) efp, EvidenceCallback);
  }
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (efp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, efp->input_entityID, 0, 0);
  Remove (efp->form);
}

static void EvidenceMessageProc (ForM f, Int2 mssg)

{
  EvidenceFormPtr  efp;

  efp = (EvidenceFormPtr) GetObjectExtra (f);
  if (efp != NULL) {
    if (efp->appmessage != NULL) {
      efp->appmessage (f, mssg);
    }
  }
}

static void CleanupEvidencePage (GraphiC g, VoidPtr data)

{
  EvidenceFormPtr  efp;

  efp = (EvidenceFormPtr) data;
  if (efp != NULL) {
    ValNodeFreeData (efp->head);
  }
  StdCleanupFormProc (g, data);
}

extern void EditEvidenceFlag (IteM i)

{
  BaseFormPtr        bfp;
  ButtoN             b;
  GrouP              c;
  FeatDefPtr         curr;
  EvidenceFormPtr    efp;
  GrouP              g;
  GrouP              h;
  ValNodePtr         head;
  GrouP              k;
  Uint1              key;
  CharPtr            label = NULL;
  Int2               listHeight;
  GrouP              q;
  SeqEntryPtr        sep;
  StdEditorProcsPtr  sepp;
  Uint2              subtype;
  ValNodePtr         vnp;
  WindoW             w;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  efp = (EvidenceFormPtr) MemNew (sizeof (EvidenceFormData));
  if (efp == NULL) return;
  w = FixedWindow (-50, -33, -10, -10, "Feature Evidence", StdCloseWindowProc);
  SetObjectExtra (w, efp, CleanupEvidencePage);
  efp->form = (ForM) w;
  efp->formmessage = EvidenceMessageProc;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    efp->appmessage = sepp->handleMessages;
  }

  efp->input_entityID = bfp->input_entityID;
  efp->input_itemID = bfp->input_itemID;
  efp->input_itemtype = bfp->input_itemtype;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  g = HiddenGroup (h, 0, 2, NULL);
  StaticPrompt (g, "Feature", 0, 0, programFont, 'c');
  if (indexerVersion) {
    listHeight = 16;
  } else {
    listHeight = 8;
  }
  efp->objlist = SingleList (g, 16, listHeight, NULL);
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
        /* if (head == NULL) {
          head = vnp;
        } */
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
      ListItem (efp->objlist, (CharPtr) vnp->data.ptrvalue);
    }
  }
  efp->head = head;

  q = HiddenGroup (h, 2, 0, NULL);
  StaticPrompt (q, "Evidence", 0, popupMenuHeight, programFont, 'c');
  efp->evdence = PopupList (q, TRUE, NULL);
  PopupItem (efp->evdence, " ");
  PopupItem (efp->evdence, "Experimental");
  PopupItem (efp->evdence, "Non-Experimental");
  SetValue (efp->evdence, 1);

  k = HiddenGroup (h, 0, 2, NULL);
  StaticPrompt (k, "Optional string constraint", 0, dialogTextHeight, programFont, 'c');
  efp->findthis = DialogText (k, "", 14, NULL);

  c = HiddenGroup (h, 4, 0, NULL);
  b = PushButton (c, "Accept", DoEvidence);
  SetObjectExtra (b, efp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) q, (HANDLE) k, (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Update ();
}

typedef struct exceptionformdata {
  FEATURE_FORM_BLOCK

  LisT           objlist;
  TexT           findthis;
  Uint2          itemtype;
  Uint2          subtype;
  GrouP          xception;
  TexT           xceptText;
  ButtoN         rescueExpl;
  Boolean        excpt;
  CharPtr        except_text;
  Boolean        rescue;
  ObjMgrPtr      omp;
  ObjMgrTypePtr  omtp;
  ValNodePtr     head;
  Boolean        stringfound;
  Char           findStr [128];
} ExceptionFormData, PNTR ExceptionFormPtr;

static void LIBCALLBACK AsnWriteExceptionCallBack (AsnExpOptStructPtr pAEOS)

{
  ExceptionFormPtr  efp;
  CharPtr          pchFind;
  CharPtr          pchSource;

  efp = (ExceptionFormPtr) pAEOS->data;
  if (ISA_STRINGTYPE (AsnFindBaseIsa (pAEOS->atp))) {
	pchSource = (CharPtr) pAEOS->dvp->ptrvalue;
	pchFind = efp->findStr;
	if (StringSearch (pchSource, pchFind) != NULL) {
	  efp->stringfound = TRUE;
	}
  }
}

static Boolean ExceptionHasSubstring (ObjMgrTypePtr omtp, AsnIoPtr aip, Pointer ptr, ExceptionFormPtr efp)

{
  efp->stringfound = FALSE;
  (omtp->asnwrite) (ptr, aip, NULL);
  return efp->stringfound;
}

static void ExceptionCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  AsnExpOptPtr      aeop;
  AsnIoPtr          aip;
  BioseqPtr         bsp;
  BioseqSetPtr      bssp;
  ExceptionFormPtr  efp;
  Boolean           notext;
  ObjMgrTypePtr     omtp;
  SeqAnnotPtr       sap;
  SeqFeatPtr        sfp;
  CharPtr           str;
  Uint2             subtype;

  if (mydata == NULL) return;
  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  efp = (ExceptionFormPtr) mydata;
  if (efp == NULL) return;
  omtp = efp->omtp;
  if (omtp == NULL || omtp->subtypefunc == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
  } else return;
  GetTitle (efp->findthis, efp->findStr, sizeof (efp->findStr) - 1);
  notext = StringHasNoText (efp->findStr);
  aip = AsnIoNullOpen ();
  aeop = AsnExpOptNew (aip, NULL, NULL, AsnWriteExceptionCallBack);
  if (aeop != NULL) {
    aeop->user_data = (Pointer) efp;
  }
  while (sap != NULL) {
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      while (sfp != NULL) {
        subtype = (*(omtp->subtypefunc)) ((Pointer) sfp);
        if (efp->subtype == 0 || subtype == efp->subtype ||
           (efp->subtype == FEATDEF_IMP &&
            subtype >= FEATDEF_allele && subtype <= FEATDEF_site_ref)) {
          if (notext || ExceptionHasSubstring (omtp, aip, (Pointer) sfp, efp)) {
            sfp->excpt = efp->excpt;
            if (sfp->excpt) {
              if (! StringHasNoText (efp->except_text)) {
                sfp->except_text = MemFree (sfp->except_text);
                sfp->except_text = StringSave (efp->except_text);
              }
            } else {
              if (efp->rescue) {
                if (sfp->comment == NULL) {
                  sfp->comment = sfp->except_text;
                } else {
                  str = MemNew (StringLen (sfp->comment) + StringLen (sfp->except_text) + 5);
                  if (str != NULL) {
                    StringCpy (str, sfp->comment);
                    StringCat (str, "; ");
                    StringCat (str, sfp->except_text);
                    sfp->comment = MemFree (sfp->comment);
                    sfp->comment = str;
                  }
                }
                sfp->except_text = NULL;
              } else {
                sfp->except_text = MemFree (sfp->except_text);
              }
            }
          }
        }
        sfp = sfp->next;
      }
    }
    sap = sap->next;
  }
  AsnIoClose (aip);
}

static void DoException (ButtoN b)

{
  ExceptionFormPtr  efp;
  SeqEntryPtr      sep;
  Uint2            subtype;
  Int2             val;
  ValNodePtr       vnp;

  efp = GetObjectExtra (b);
  if (efp == NULL) return;
  sep = GetTopSeqEntryForEntityID (efp->input_entityID);
  if (sep == NULL) return;
  Hide (efp->form);
  WatchCursor ();
  Update ();
  efp->itemtype = OBJ_SEQFEAT;
  subtype = 0;
  vnp = NULL;
  val = GetValue (efp->objlist);
  if (val > 0) {
    vnp = efp->head;
    while (vnp != NULL && val > 1) {
      val--;
      vnp = vnp->next;
    }
  }
  if (vnp != NULL) {
    subtype = vnp->choice;
  }
  efp->omp = ObjMgrGet ();
  efp->omtp = NULL;
  if (efp->omp != NULL) {
    efp->omtp = ObjMgrTypeFind (efp->omp, efp->itemtype, NULL, NULL);
  }
  efp->subtype = subtype;
  efp->excpt = (Boolean) (GetValue (efp->xception) == 2);
  if (efp->excpt) {
    efp->except_text = SaveStringFromText (efp->xceptText);
  } else {
    efp->except_text = NULL;
  }
  efp->rescue = GetStatus (efp->rescueExpl);
  if (efp->itemtype != 0 && efp->omtp != NULL) {
    SeqEntryExplore (sep, (Pointer) efp, ExceptionCallback);
  }
  MemFree (efp->except_text);
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (efp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, efp->input_entityID, 0, 0);
  Remove (efp->form);
}

static void ExceptionMessageProc (ForM f, Int2 mssg)

{
  ExceptionFormPtr  efp;

  efp = (ExceptionFormPtr) GetObjectExtra (f);
  if (efp != NULL) {
    if (efp->appmessage != NULL) {
      efp->appmessage (f, mssg);
    }
  }
}

static void CleanupExceptionPage (GraphiC g, VoidPtr data)

{
  ExceptionFormPtr  efp;

  efp = (ExceptionFormPtr) data;
  if (efp != NULL) {
    ValNodeFreeData (efp->head);
  }
  StdCleanupFormProc (g, data);
}

extern void EditExceptionFlag (IteM i)

{
  BaseFormPtr        bfp;
  ButtoN             b;
  GrouP              c;
  FeatDefPtr         curr;
  ExceptionFormPtr   efp;
  GrouP              g;
  GrouP              h;
  ValNodePtr         head;
  GrouP              k;
  Uint1              key;
  CharPtr            label = NULL;
  Int2               listHeight;
  PrompT             ppt;
  GrouP              q;
  SeqEntryPtr        sep;
  StdEditorProcsPtr  sepp;
  Uint2              subtype;
  ValNodePtr         vnp;
  WindoW             w;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  efp = (ExceptionFormPtr) MemNew (sizeof (ExceptionFormData));
  if (efp == NULL) return;
  w = FixedWindow (-50, -33, -10, -10, "Feature Exception", StdCloseWindowProc);
  SetObjectExtra (w, efp, CleanupExceptionPage);
  efp->form = (ForM) w;
  efp->formmessage = ExceptionMessageProc;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    efp->appmessage = sepp->handleMessages;
  }

  efp->input_entityID = bfp->input_entityID;
  efp->input_itemID = bfp->input_itemID;
  efp->input_itemtype = bfp->input_itemtype;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  g = HiddenGroup (h, 0, 2, NULL);
  StaticPrompt (g, "Feature", 0, 0, programFont, 'c');
  if (indexerVersion) {
    listHeight = 16;
  } else {
    listHeight = 8;
  }
  efp->objlist = SingleList (g, 16, listHeight, NULL);
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
        /* if (head == NULL) {
          head = vnp;
        } */
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
      ListItem (efp->objlist, (CharPtr) vnp->data.ptrvalue);
    }
  }
  efp->head = head;

  q = HiddenGroup (h, -2, 0, NULL);
  ppt = StaticPrompt (q, "Exception  ", 0, stdLineHeight, programFont, 'l');
  efp->xception = HiddenGroup (q, -2, 0, NULL);
  SetObjectExtra (efp->xception, efp, NULL);
  RadioButton (efp->xception, "Clear");
  RadioButton (efp->xception, "Set");
  SetValue (efp->xception, 1);
  StaticPrompt (q, "Explanation", 0, dialogTextHeight, programFont, 'l');
  efp->xceptText = DialogText (q, "", 12, NULL);
  AlignObjects (ALIGN_MIDDLE, (HANDLE) ppt, (HANDLE) efp->xception, NULL);
  efp->rescueExpl = CheckBox (q, "Move explanation to comment", NULL);

  k = HiddenGroup (h, 0, 2, NULL);
  StaticPrompt (k, "Optional string constraint", 0, dialogTextHeight, programFont, 'c');
  efp->findthis = DialogText (k, "", 14, NULL);

  c = HiddenGroup (h, 4, 0, NULL);
  b = PushButton (c, "Accept", DoException);
  SetObjectExtra (b, efp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) q, (HANDLE) k, (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Update ();
}

static void BreakIntoAGroup (BioseqSetPtr parent, Uint1 _class, SeqEntryPtr list)

{
  BioseqSetPtr  bssp;
  Int2          count;
  SeqEntryPtr   sep;
  SeqEntryPtr   tmp;

  while (list != NULL) {
    bssp = BioseqSetNew ();
    if (bssp == NULL) return;
    bssp->_class = _class;
    sep = SeqEntryNew ();
    if (sep == NULL) return;
    sep->choice = 2;
    sep->data.ptrvalue = (Pointer) bssp;
    if (parent->seq_set != NULL) {
      tmp = parent->seq_set;
      while (tmp->next != NULL) {
        tmp = tmp->next;
      }
      tmp->next = sep;
    } else {
      parent->seq_set = sep;
    }
    bssp->seq_set = list;
    for (tmp = list, count = 0; tmp != NULL && count < 99; tmp = tmp->next, count++) continue;
    if (tmp != NULL) {
      list = tmp->next;
      tmp->next = NULL;
    } else {
      list = NULL;
    }
  }
}

extern void MakeGroupsOf200 (IteM i)

{
  AsnIoPtr       aip;
  BaseFormPtr    bfp;
  BioseqSetPtr   bssp;
  Uint1          _class;
  Int2           count;
  Char           file [FILENAME_MAX];
  SeqEntryPtr    list;
  SeqEntryPtr    next;
  ObjMgrDataPtr  omdptop;
  ObjMgrData     omdata;
  Char           output [PATH_MAX];
  Uint2          parenttype;
  Pointer        parentptr;
  Char           path [PATH_MAX];
  CharPtr        ptr;
  SeqEntryPtr    sep;
  SeqSubmitPtr   ssp;
  Char           str [FILENAME_MAX];
  SeqEntryPtr    tmp;
#ifdef WIN_MAC
  FILE           *f;
#endif

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  if (! IS_Bioseq_set (sep)) return;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp == NULL) return;
  _class = bssp->_class;
  if (_class != 7 && _class != 13 && _class != 14 && _class != 15) return;

  SaveSeqEntryObjMgrData (sep, &omdptop, &omdata);
  GetSeqEntryParent (sep, &parentptr, &parenttype);

  list = bssp->seq_set;
  bssp->seq_set = NULL;
  bssp->_class = 7;
  BreakIntoAGroup (bssp, _class, list);

  SeqMgrLinkSeqEntry (sep, parenttype, parentptr);
  RestoreSeqEntryObjMgrData (sep, omdptop, &omdata);
  PropagateFromGenBankBioseqSet (sep, TRUE);

  if (parenttype == OBJ_SEQSUB) {
    if (GetOutputFileName (path, sizeof (path), "")) {
      ssp = (SeqSubmitPtr) parentptr;
      if (ssp != NULL && ssp->datatype == 1) {
        sep = (SeqEntryPtr) ssp->data;
        ptr = StringRChr (path, DIRDELIMCHR);
        if (ptr != NULL) {
          ptr++;
          StringNCpy_0 (file, ptr, sizeof (file));
          *ptr = '\0';
          tmp = bssp->seq_set;
          count = 0;
          while (tmp != NULL) {
            next = tmp->next;
            tmp->next = NULL;
            ssp->data = (Pointer) tmp;
            StringCpy (output, path);
            count++;
            if (count < 10) {
              sprintf (str, "%s0%1d", file, (int) count);
            } else {
              sprintf (str, "%s%2d", file, (int) count);
            }
            FileBuildPath (output, NULL, str);
#ifdef WIN_MAC
            f = FileOpen (output, "r");
            if (f != NULL) {
              FileClose (f);
            } else {
              FileCreate (output, "TEXT", "ttxt");
            }
#endif
            aip = AsnIoOpen (output, "w");
            if (aip != NULL) {
              SeqSubmitAsnWrite (ssp, aip, NULL);
            }
            AsnIoClose (aip);
            tmp->next = next;
            tmp = next;
          }
          ssp->data = (Pointer) sep;
        }
      }
    }
  }

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update ();
}

extern void ParseInNucUpdates (IteM i)

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
  Message (MSG_OK, "Not yet implemented");
}

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

NLM_EXTERN void CdTransCheck(ValidStructPtr vsp, SeqFeatPtr sfp);
static Boolean CdsTranslatesProperly (SeqFeatPtr sfp, GatherContextPtr gcp)

{
  Int2            errors;
  Int2            j;
  ErrSev          oldErrSev;
  ValidStructPtr  vsp;

  if (sfp == NULL) return TRUE;
  vsp = ValidStructNew ();
  if (vsp == NULL) return TRUE;
  vsp->gcp = gcp;
  oldErrSev = ErrSetMessageLevel (SEV_MAX);
  CdTransCheck (vsp, sfp);
  ErrSetMessageLevel (oldErrSev);
  ErrClear ();
  ErrShow ();
  errors = 0;
  for (j = 0; j < 6; j++) {
    errors += vsp->errors [j];
  }
  ValidStructFree (vsp);
  return (Boolean) (errors == 0);
}

typedef struct recompdata {
  Int4        count;
  MonitorPtr  mon;
  BioseqPtr   batchbsp;
} RecompData, PNTR RecompDataPtr;

static Int2 GeneticCodeFromCrp (CdRegionPtr crp)

{
  Int2            code;
  GeneticCodePtr  gcp;
  Char            name [256];
  ValNodePtr      tmp;

  code = 0;
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
          code = tmp->data.intvalue;
          break;
        default :
          break;
      }
    }
    if (code == 0) {
      gcp = GeneticCodeFind (code, name);
      if (gcp != NULL) {
        for (tmp = (ValNodePtr) gcp->data.ptrvalue; tmp != NULL; tmp = tmp->next) {
          switch (tmp->choice) {
            case 2 :
              code = tmp->data.intvalue;
              break;
            default :
              break;
          }
        }
      }
    }
  }
  return code;
}

static Boolean RecomputeSuggCallback (GatherContextPtr gcp)

{
  Int2           code;
  CdRegionPtr    crp;
  BioseqPtr      nucbsp;
  BioseqPtr      protbsp;
  RecompDataPtr  rdp;
  SeqFeatPtr     sfp;
  SeqIdPtr       sip;
  SeqLocPtr      slp;
  Char           str [256];
  Char           tmp [256];

  if (gcp == NULL) return TRUE;
  if (gcp->thistype != OBJ_SEQFEAT) return TRUE;
  rdp = (RecompDataPtr) gcp->userdata;
  if (rdp == NULL) return TRUE;
  sfp = (SeqFeatPtr) gcp->thisitem;
  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return TRUE;
  crp = (CdRegionPtr) sfp->data.value.ptrvalue;
  if (crp == NULL) return TRUE;
  /* if (CdsTranslatesProperly (sfp, gcp)) return TRUE; */

  code = GeneticCodeFromCrp (crp);

  nucbsp = GetBioseqGivenSeqLoc (sfp->location, gcp->entityID);
  if (nucbsp != NULL && rdp->batchbsp != NULL && nucbsp != rdp->batchbsp) {
    ClearBatchSuggestNucleotide ();
    rdp->batchbsp = NULL;
    Message (MSG_POSTERR, "Recompute Suggest is reverting to slower processing");
  }
  sip = SeqLocId (sfp->product);
  if (sip != NULL) {
    protbsp = BioseqFind (sip);
    if (nucbsp != NULL && protbsp != NULL &&
        ISA_na (nucbsp->mol) && ISA_aa (protbsp->mol) &&
        nucbsp->length > 0 && protbsp->length > 0) {
      str [0] = '\0';
      tmp [0] = '\0';
      sip = SeqIdFindWorst (protbsp->id);
      SeqIdWrite (sip, tmp, PRINTID_REPORT, sizeof (tmp));
      (rdp->count)++;
      sprintf (str, "Processing sequence %d [%s]", (int) rdp->count, tmp);
      if (rdp->mon != NULL) {
        MonitorStrValue (rdp->mon, str);
        Update ();
      }
      slp = PredictCodingRegion (nucbsp, protbsp, code);
      if (slp == NULL) return TRUE;
      sfp->location = SeqLocFree (sfp->location);
      sfp->location = slp;
      sfp->partial = LocationHasNullsBetween (sfp->location);
    }
  }
  return TRUE;
}

extern void RecomputeSuggest (IteM i)

{
  BaseFormPtr  bfp;
  Int2         code;
  GatherScope  gs;
  SeqEntryPtr  nucsep;
  RecompData   rd;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  WatchCursor ();
  Update ();
  rd.count = 0;
  rd.mon = MonitorStrNewEx ("Correcting Coding Regions", 20, FALSE);
  rd.batchbsp = NULL;
  nucsep = FindNucSeqEntry (sep);
  if (nucsep != NULL && IS_Bioseq (nucsep)) {
    rd.batchbsp = (BioseqPtr) nucsep->data.ptrvalue;
  }
  if (rd.batchbsp != NULL) {
    code = SeqEntryToGeneticCode (sep, NULL, NULL, 0);
    SetBatchSuggestNucleotide (rd.batchbsp, code);
  }
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = FALSE;
  MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  GatherEntity (bfp->input_entityID, (Pointer) (&rd), RecomputeSuggCallback, &gs);
  MonitorFree (rd.mon);
  if (rd.batchbsp != NULL) {
    ClearBatchSuggestNucleotide ();
  }
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

extern Boolean RetranslateOneCDS (SeqFeatPtr sfp, Uint2 entityID)

{
  SeqFeatPtr    bestprot;
  ByteStorePtr  bs;
  BioseqPtr     bsp;
  Char          ch;
  SeqEntryPtr   master;
  MolInfoPtr    mip;
  SeqEntryPtr   old;
  Boolean       partial5;
  Boolean       partial3;
  CharPtr       prot;
  CharPtr       ptr;
  SeqEntryPtr   sep;
  SeqIdPtr      sip;
  SeqLocPtr     slp;
  ValNodePtr    vnp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return TRUE;
  slp = SeqLocFindNext (sfp->location, NULL);
  if (slp == NULL) return TRUE;
  CheckSeqLocForPartial (slp, &partial5, &partial3);
  /* if (CdsTranslatesProperly (sfp, gcp)) return TRUE; */

  if (sfp->product == NULL) {
    master = NULL;
    old = NULL;
    bsp = GetBioseqGivenSeqLoc (sfp->location, entityID);
    if (bsp != NULL) {
      master = GetBestTopParentForData (entityID, bsp);
    }
    bsp = BioseqNew ();
    if (bsp != NULL) {
      bsp->mol = Seq_mol_aa;
      bsp->repr = Seq_repr_raw;
      bsp->seq_data_type = Seq_code_ncbieaa;
      bsp->length = 0;
      bsp->seq_data = BSNew (0);
      if (master != NULL) {
        old = SeqEntrySetScope (master);
      }
      bsp->id = MakeNewProteinSeqId (sfp->location, NULL);
      SeqMgrAddToBioseqIndex (bsp);
      if (master != NULL) {
        SeqEntrySetScope (old);
      }
      sep = SeqEntryNew ();
      if (sep != NULL) {
        sep->choice = 1;
        sep->data.ptrvalue = (Pointer) bsp;
        SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, sep);
      }
      SetSeqFeatProduct (sfp, bsp);
      if (master != NULL && sep != NULL) {
        AddSeqEntryToSeqEntry (master, sep, TRUE);
      }
    }
  }

  sip = SeqLocId (sfp->product);
  if (sip != NULL) {
    bsp = BioseqFind (sip);
    if (bsp != NULL && ISA_aa (bsp->mol) && bsp->repr == Seq_repr_raw) {
      bestprot = FindBestProtein (entityID, sfp->product);
      bs = ProteinFromCdRegionEx (sfp, FALSE, FALSE);
      if (bs == NULL) return TRUE;
      prot = BSMerge (bs, NULL);
      bs = BSFree (bs);
      if (prot == NULL) return TRUE;
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
      bsp->repr = Seq_repr_raw;
      bsp->mol = Seq_mol_aa;
      bsp->seq_data_type = Seq_code_ncbieaa;
      bsp->seq_data = BSFree (bsp->seq_data);
      bsp->seq_data = bs;
      bsp->length = BSLen (bs);
      sep = SeqMgrGetSeqEntryForData (bsp);
      if (sep == NULL) return TRUE;
      if (bestprot != NULL) {
        bestprot->location = SeqLocFree (bestprot->location);
        bestprot->location = CreateWholeInterval (sep);
        SetSeqLocPartial (bestprot->location, partial5, partial3);
        bestprot->partial = (partial5 || partial3);
      }
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
          if (partial5 && partial3) {
            mip->completeness = 5;
          } else if (partial5) {
            mip->completeness = 3;
          } else if (partial3) {
            mip->completeness = 4;
          /*
          } else if (partial) {
            mip->completeness = 2;
          */
          } else {
            mip->completeness = 0;
          }
        }
      }
    }
  }
  return TRUE;
}

static Boolean RetranslateCDSCallback (GatherContextPtr gcp)

{
  RecompDataPtr  rdp;
  SeqFeatPtr     sfp;

  if (gcp == NULL) return TRUE;
  if (gcp->thistype != OBJ_SEQFEAT) return TRUE;
  rdp = (RecompDataPtr) gcp->userdata;
  if (rdp == NULL) return TRUE;
  sfp = (SeqFeatPtr) gcp->thisitem;
  return RetranslateOneCDS (sfp, gcp->entityID);
}

extern void RetranslateCdRegions (IteM i)

{
  BaseFormPtr  bfp;
  GatherScope  gs;
  RecompData   rd;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  WatchCursor ();
  Update ();
  rd.count = 0;
  rd.mon = MonitorStrNewEx ("Correcting Coding Regions", 20, FALSE);
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = FALSE;
  MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  GatherEntity (bfp->input_entityID, (Pointer) (&rd), RetranslateCDSCallback, &gs);
  MonitorFree (rd.mon);
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static Boolean CorrectGenCodeCallback (GatherContextPtr gcp)

{
  CdRegionPtr     crp;
  GeneticCodePtr  gc;
  Int2            genCode;
  SeqEntryPtr     sep;
  SeqFeatPtr      sfp;
  ValNodePtr      vnp;

  if (gcp == NULL) return TRUE;
  sep = (SeqEntryPtr) gcp->userdata;
  if (sep == NULL ) return TRUE;
  if (gcp->thistype != OBJ_SEQFEAT) return TRUE;
  sfp = (SeqFeatPtr) gcp->thisitem;
  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return TRUE;
  crp = (CdRegionPtr) sfp->data.value.ptrvalue;
  if (crp == NULL) return TRUE;
  genCode = SeqEntryToGeneticCode (sep, NULL, NULL, 0);
  gc = GeneticCodeNew ();
  if (gc == NULL) return TRUE;
  crp->genetic_code = GeneticCodeFree (crp->genetic_code);
  vnp = ValNodeNew (NULL);
  gc->data.ptrvalue = vnp;
  if (vnp != NULL) {
    vnp->choice = 2;
    vnp->data.intvalue = (Int4) genCode;
  }
  crp->genetic_code = gc;
  return TRUE;
}

extern void CorrectGenCodes (SeqEntryPtr sep, Uint2 entityID)

{
  BioseqSetPtr  bssp;
  GatherScope   gs;

  if (sep == NULL) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && (bssp->_class == 7 ||
                         (bssp->_class >= 13 && bssp->_class <= 16))) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        CorrectGenCodes (sep, entityID);
      }
      return;
    }
  }
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = FALSE;
  MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  gs.scope = sep;
  GatherEntity (entityID, (Pointer) sep, CorrectGenCodeCallback, &gs);
}

extern void CorrectCDSGenCodes (IteM i)

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
  CorrectGenCodes (sep, bfp->input_entityID);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

static BioseqPtr SqnGetBioseqGivenSeqLoc (SeqLocPtr slp, Uint2 entityID)

{
  BioseqPtr    bsp;
  SeqEntryPtr  sep;
  SeqIdPtr     sip;
  SeqLocPtr    tmp;

  if (slp == NULL) return NULL;
  bsp = NULL;
  sip = SeqLocId (slp);
  if (sip != NULL) {
    bsp = BioseqFind (sip);
  } else {
    tmp = SeqLocFindNext (slp, NULL);
    if (tmp != NULL) {
      sip = SeqLocId (tmp);
      if (sip != NULL) {
        bsp = BioseqFind (sip);
        if (bsp != NULL) {
          sep = SeqMgrGetSeqEntryForData (bsp);
          entityID = ObjMgrGetEntityIDForChoice (sep);
          bsp = GetBioseqGivenSeqLoc (slp, entityID);
        }
      }
    }
  }
  return bsp;
}

static BioseqPtr GetBioseqReferencedByAnnot (SeqAnnotPtr sap, Uint2 entityID)

{
  SeqAlignPtr   align;
  BioseqPtr     bsp;
  DenseDiagPtr  ddp;
  DenseSegPtr   dsp;
  SeqFeatPtr    feat;
  SeqGraphPtr   graph;
  SeqIdPtr      sip;
  SeqLocPtr     slp;
  StdSegPtr     ssp;
  SeqLocPtr     tloc;

  if (sap == NULL) return NULL;
  switch (sap->type) {
    case 1 :
      feat = (SeqFeatPtr) sap->data;
      while (feat != NULL) {
        slp = feat->location;
        if (slp != NULL) {
          bsp = SqnGetBioseqGivenSeqLoc (slp, entityID);
          if (bsp != NULL) return bsp;
        }
        feat = feat->next;
      }
      break;
    case 2 :
      align = (SeqAlignPtr) sap->data;
      while (align != NULL) {
        if (align->segtype == 1) {
          ddp = (DenseDiagPtr) align->segs;
          if (ddp != NULL) {
            for (sip = ddp->id; sip != NULL; sip = sip->next) {
              bsp = BioseqFind (sip);
              if (bsp != NULL) return bsp;
            }
          }
        } else if (align->segtype == 2) {
          dsp = (DenseSegPtr) align->segs;
          if (dsp != NULL) {
            for (sip = dsp->ids; sip != NULL; sip = sip->next) {
              bsp = BioseqFind (sip);
              if (bsp != NULL) return bsp;
            }
          }
        } else if (align->segtype == 3) {
          ssp = (StdSegPtr) align->segs;
          if (ssp != NULL && ssp->loc != NULL) {
            for (tloc = ssp->loc; tloc != NULL; tloc = tloc->next) {
              bsp = BioseqFind (SeqLocId (tloc));
              if (bsp != NULL) return bsp;
            }
          }
        }
        align = align->next;
      }
      break;
    case 3 :
      graph = (SeqGraphPtr) sap->data;
      while (graph != NULL) {
        slp = graph->loc;
        if (slp != NULL) {
          bsp = SqnGetBioseqGivenSeqLoc (slp, entityID);
          if (bsp != NULL) return bsp;
        }
        graph = graph->next;
      }
      break;
    default :
      break;
  }
  return NULL;
}

static Int4 GetScore (ScorePtr score)

{
  ObjectIdPtr  id;

  while (score != NULL) {
    id = score->id;
    if (id != NULL) {
      if (StringICmp (id->str, "score") == 0) {
        if (score->choice == 1) {
          return (score->value.intvalue);
        }
      }
    }
    score = score->next;
  }
  return 0;
}

static Int4 FindScore (SeqAlignPtr align)

{
  if (align == NULL) return 0;
  if (align->score != NULL) {
    return GetScore (align->score);
  }
  return 0;
}

static int LIBCALLBACK SortByScoreCallback (VoidPtr ptr1, VoidPtr ptr2)

{
  SeqAlignPtr   sap1;
  SeqAlignPtr   sap2;
  Int4          score1;
  Int4          score2;

  if (ptr1 != NULL && ptr2 != NULL) {
    sap1 = *((SeqAlignPtr PNTR) ptr1);
    sap2 = *((SeqAlignPtr PNTR) ptr2);
    if (sap1 != NULL && sap2 != NULL) {
      score1 = FindScore (sap1);
      score2 = FindScore (sap2);
      if (score1 < score2) {
        return 1;
      } else if (score1 > score2) {
        return -1;
      } else {
        return 0;
      }
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}

static SeqAlignPtr SortBySeqAlignScore (SeqAlignPtr list)

{
  SeqAlignPtr  align;
  Int4         count, i;
  SeqAlignPtr  PNTR head;

  if (list == NULL) return 0;
  count = 0;
  for (align = list; align != NULL; align = align->next) {
    count++;
  }
  head = MemNew (sizeof (SeqAlignPtr) * (size_t) (count + 1));
  if (head == NULL) return 0;
  for (align = list, i = 0; align != NULL && i < count; i++) {
    head [i] = align;
    align = align->next;
  }
  HeapSort (head, (size_t) count, sizeof (SeqAlignPtr), SortByScoreCallback);
  for (i = 0; i < count; i++) {
    align = head [i];
    align->next = head [i + 1];
  }
  list = head [0];
  MemFree (head);
  return list;
}

static void TakeTop10Alignments (SeqAnnotPtr sap)

{
  SeqAlignPtr  align;
  MsgAnswer    ans;
  Int2         count;
  SeqAlignPtr  next;

  if (sap == NULL || sap->type != 2 || sap->data == NULL) return;
  count = 0;
  for (align = (SeqAlignPtr) sap->data; align != NULL; align = align->next) {
    count++;
  }
  if (count <= 10) return;
  ans = Message (MSG_YN, "Do you want to take only the top 10 (out of %d) alignments?", (int) count);
  if (ans == ANS_NO) return;
  sap->data = SortBySeqAlignScore ((SeqAlignPtr) sap->data);
  for (align = (SeqAlignPtr) sap->data, count = 0; align != NULL && count < 10; align = align->next) {
    count++;
  }
  next = align->next;
  align->next = NULL;
  align = next;
  while (align != NULL) {
    next = align->next;
    align->next = NULL;
    SeqAlignFree (align);
    align = next;
  }
}

static void DoOnePub (PubdescPtr pdp)

{
  ValNodePtr    citartptr = NULL;
  Int4          muid = 0;
  Int4          pmid = 0;
  ValNodePtr    tmp = NULL;
  ValNodePtr    vnp;

  if (pdp != NULL) {
    for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
      if (vnp->choice == PUB_Muid) {
        muid = vnp->data.intvalue;
      } else if (vnp->choice == PUB_PMid) {
        pmid = vnp->data.intvalue;
      } else if (vnp->choice == PUB_Article) {
        citartptr = vnp;
      }
    }
    if (pmid != 0) {
      tmp = MedArchGetPubPmId (pmid);
      muid = MedArchPm2Mu (pmid);
    } else if (muid != 0) {
      tmp = MedArchGetPub (muid);
      pmid = MedArchMu2Pm (muid);
    } else if (citartptr != NULL) {
      muid = MedArchCitMatch (citartptr);
      if (muid != 0) {
        tmp = MedArchGetPub (muid);
        pmid = MedArchMu2Pm (muid);
      }
    }
    if (tmp != NULL) {
      MedlineToISO (tmp);
      if (pmid != 0) {
        ValNodeAddInt (&tmp, PUB_PMid, pmid);
      }
      if (muid != 0) {
        ValNodeAddInt (&tmp, PUB_Muid, muid);
      }
      pdp->pub = PubEquivFree (pdp->pub);
      pdp->pub = tmp;
    }
  }
}

static void DoLookupPub (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  PubdescPtr    pdp;
  SeqAnnotPtr   sap;
  ValNodePtr    sdp;
  SeqFeatPtr    sfp;

  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
    sdp = bsp->descr;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
    sdp = bssp->descr;
  } else return;
  while (sap != NULL) {
    if (sap->type == 1) {
      for (sfp = (SeqFeatPtr) sap->data; sfp != NULL; sfp = sfp->next) {
        if (sfp->data.choice == SEQFEAT_PUB) {
          pdp = (PubdescPtr) sfp->data.value.ptrvalue;
          DoOnePub (pdp);
        }
      }
    }
    sap = sap->next;
  }
  while (sdp != NULL) {
    if (sdp->choice == Seq_descr_pub) {
      pdp = (PubdescPtr) sdp->data.ptrvalue;
      DoOnePub (pdp);
    }
    sdp = sdp->next;
  }
}

extern void LookupAllPubs (IteM i);
extern void LookupAllPubs (IteM i)

{
  BaseFormPtr  bfp;
  MonitorPtr   mon = NULL;
  SeqEntryPtr  sep;
  ErrSev       sev;


  if (! useMedarch) return;
#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  sev = ErrSetMessageLevel (SEV_FATAL);
  WatchCursor ();
  mon = MonitorStrNewEx ("Processing Publications", 40, FALSE);
  MonitorStrValue (mon, "Connecting to MedArch");
  Update ();
  if (! MedArchInit ()) {
    MonitorFree (mon);
    ArrowCursor ();
    Update ();
    Message (MSG_POST, "Unable to connect to MedArch");
    return;
  }
  SeqEntryExplore (sep, NULL, DoLookupPub);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  MonitorStrValue (mon, "Closing MedArch");
  Update ();
  MedArchFini ();
  MonitorFree (mon);
  ArrowCursor ();
  Update ();
  ErrSetMessageLevel (sev);
  ErrClear ();
  ErrShow ();
}

static void LookupPublications (SeqAnnotPtr sap)

{
  MonitorPtr  mon = NULL;
  PubdescPtr  pdp;
  SeqFeatPtr  sfp;
  ValNodePtr  tmp;
  Int4        uid;
  Boolean     usingMedarch = FALSE;
  ValNodePtr  vnp;

  if (! useMedarch) return;
  if (sap == NULL || sap->type != 1) return;
  for (sfp = (SeqFeatPtr) sap->data; sfp != NULL; sfp = sfp->next) {
    if (sfp->data.choice == SEQFEAT_PUB) {
      pdp = (PubdescPtr) sfp->data.value.ptrvalue;
      if (pdp != NULL) {
        vnp = pdp->pub;
        if (vnp != NULL && vnp->next == NULL) {
          if (vnp->choice == PUB_Muid || vnp->choice == PUB_PMid) {
            if (! usingMedarch) {
              WatchCursor ();
              mon = MonitorStrNewEx ("Processing Publications", 40, FALSE);
              MonitorStrValue (mon, "Connecting to MedArch");
              Update ();
              if (MedArchInit ()) {
                usingMedarch = TRUE;
              } else {
                MonitorFree (mon);
                ArrowCursor ();
                Update ();
                Message (MSG_POST, "Unable to connect to MedArch");
                return;
              }
            }
          }
          tmp = NULL;
          if (vnp->choice == PUB_Muid) {
            uid = vnp->data.intvalue;
            tmp = MedArchGetPub (uid);
          } else if (vnp->choice == PUB_PMid) {
            uid = vnp->data.intvalue;
            tmp = MedArchGetPubPmId (uid);
          }
          if (tmp != NULL) {
            MedlineToISO (tmp);
            tmp->next = vnp;
            pdp->pub = tmp;
          }
        }
      }
    }
  }
  if (usingMedarch) {
    MonitorStrValue (mon, "Closing MedArch");
    Update ();
    MedArchFini ();
    MonitorFree (mon);
    ArrowCursor ();
    Update ();
  }
}

static void PromotePubs (SeqFeatPtr first, BioseqPtr bsp, Uint2 entityID)

{
  MsgAnswer    ans;
  Boolean      asked = FALSE;
  PubdescPtr   pdp;
  SeqDescrPtr  sdp;
  SeqEntryPtr  sep;
  SeqFeatPtr   sfp;
  ValNode      vn;

  MemSet ((Pointer) &vn, 0, sizeof (ValNode));
  vn.choice = SEQLOC_WHOLE;
  vn.data.ptrvalue = (Pointer) SeqIdFindBest (bsp->id, 0);
  vn.next = NULL;

  for (sfp = first; sfp != NULL; sfp = sfp->next) {
    if (sfp->data.choice == SEQFEAT_PUB) {
      if (SeqLocCompare (sfp->location, &vn) == SLC_A_EQ_B) {
        if (! asked) {
          ans = Message (MSG_YN, "Do you wish to convert full-length publication features to descriptors?");
          if (ans == ANS_NO) return;
          asked = TRUE;
        }
      }
    }
  }

  sep = GetBestTopParentForData (entityID, bsp);
  for (sfp = first; sfp != NULL; sfp = sfp->next) {
    if (sfp->data.choice == SEQFEAT_PUB) {
      if (SeqLocCompare (sfp->location, &vn) == SLC_A_EQ_B) {
        sfp->idx.deleteme = TRUE;
        sfp->data.choice = SEQFEAT_COMMENT;
        pdp = (PubdescPtr) sfp->data.value.ptrvalue;
        sfp->data.value.ptrvalue = NULL;
        sdp = CreateNewDescriptor (sep, Seq_descr_pub);
        if (sdp != NULL) {
          sdp->data.ptrvalue = (Pointer) pdp;
        }
      }
    }
  }

  DeleteMarkedObjects (entityID, 0, NULL);
}

extern Uint2 SmartAttachSeqAnnotToSeqEntry (Uint2 entityID, SeqAnnotPtr sap)

{
  BioseqPtr      bsp;
  Int2           genCode;
  SeqEntryPtr    oldscope;
  OMProcControl  ompc;
  SeqEntryPtr    sep;
  SeqFeatPtr     sfp = NULL;

  if (sap == NULL) return entityID;
  bsp = GetBioseqReferencedByAnnot (sap, entityID);
  if (bsp == NULL) {
    oldscope = SeqEntrySetScope (NULL);
    if (oldscope != NULL) {
      bsp = GetBioseqReferencedByAnnot (sap, entityID);
      SeqEntrySetScope (oldscope);
    }
  }
  if (bsp != NULL) {
    sep = SeqMgrGetSeqEntryForData (bsp);
    entityID = ObjMgrGetEntityIDForChoice (sep);
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      sep = GetBestTopParentForData (entityID, bsp);
      genCode = SeqEntryToGeneticCode (sep, NULL, NULL, 0);
      SetEmptyGeneticCodes (sap, genCode);
      LookupPublications (sap);
    } else if (sap->type == 2) {
      TakeTop10Alignments (sap);
    }
    MemSet ((Pointer) &ompc, 0, sizeof (OMProcControl));
    ompc.input_entityID = entityID;
    ompc.input_itemID = GetItemIDGivenPointer (entityID, OBJ_BIOSEQ, (Pointer) bsp);
    ompc.input_itemtype = OBJ_BIOSEQ;
    ompc.output_itemtype = OBJ_SEQANNOT;
    ompc.output_data = (Pointer) sap;
    if (! AttachDataForProc (&ompc, FALSE)) {
      Message (MSG_ERROR, "SmartAttachSeqAnnotToSeqEntry failed");
    } else if (sfp != NULL) {
      PromoteXrefs (sfp, bsp, entityID);
      PromotePubs (sfp, bsp, entityID);
    }
  } else {
    Message (MSG_ERROR, "Feature table identifiers do not match record");
  }
  return entityID;
}

typedef struct removeformdata {
  FEATURE_FORM_BLOCK

  Boolean        is_feature;
  LisT           objlist;
  TexT           findthis;
  Uint2          itemtype;
  Uint2          subtype;
  ObjMgrPtr      omp;
  ObjMgrTypePtr  omtp;
  ValNodePtr     head;
  Boolean        stringfound;
  Char           findStr [128];
  ValNodePtr     bsplist;
} RemoveFormData, PNTR RemoveFormPtr;

static void LIBCALLBACK AsnWriteRemoveForDCallBack (AsnExpOptStructPtr pAEOS)

{
  CharPtr        pchFind;
  CharPtr        pchSource;
  RemoveFormPtr  rfp;

  rfp = (RemoveFormPtr) pAEOS->data;
  if (ISA_STRINGTYPE (AsnFindBaseIsa (pAEOS->atp))) {
	pchSource = (CharPtr) pAEOS->dvp->ptrvalue;
	pchFind = rfp->findStr;
	if (StringSearch (pchSource, pchFind) != NULL) {
	  rfp->stringfound = TRUE;
	}
  }
}

static Boolean ObjectHasSubstring (ObjMgrTypePtr omtp, AsnIoPtr aip, Pointer ptr, RemoveFormPtr rfp)

{
  rfp->stringfound = FALSE;
  (omtp->asnwrite) (ptr, aip, NULL);
  return rfp->stringfound;
}

static void RemoveFeatureCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  AsnExpOptPtr   aeop;
  AsnIoPtr       aip;
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  SeqAnnotPtr    nextsap;
  SeqFeatPtr     nextsfp;
  Boolean        notext;
  ObjMgrTypePtr  omtp;
  Pointer PNTR   prevsap;
  Pointer PNTR   prevsfp;
  BioseqPtr      productbsp;
  RemoveFormPtr  rfp;
  SeqAnnotPtr    sap;
  SeqFeatPtr     sfp;
  SeqIdPtr       sip;
  Uint2          subtype;

  if (mydata == NULL) return;
  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  rfp = (RemoveFormPtr) mydata;
  if (rfp == NULL) return;
  omtp = rfp->omtp;
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
  GetTitle (rfp->findthis, rfp->findStr, sizeof (rfp->findStr) - 1);
  notext = StringHasNoText (rfp->findStr);
  aip = AsnIoNullOpen ();
  aeop = AsnExpOptNew (aip, NULL, NULL, AsnWriteRemoveForDCallBack);
  if (aeop != NULL) {
    aeop->user_data = (Pointer) rfp;
  }
  while (sap != NULL) {
    nextsap = sap->next;
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      prevsfp = (Pointer PNTR) &(sap->data);
      while (sfp != NULL) {
        nextsfp = sfp->next;
        subtype = (*(omtp->subtypefunc)) ((Pointer) sfp);
        if (subtype == rfp->subtype ||
           (rfp->subtype == FEATDEF_IMP &&
            subtype >= FEATDEF_allele && subtype <= FEATDEF_site_ref)) {
          if (notext || ObjectHasSubstring (omtp, aip, (Pointer) sfp, rfp)) {
            if (sfp->data.choice == SEQFEAT_CDREGION) {
              if (sfp->product != NULL) {
                sip = SeqLocId (sfp->product);
                if (sip != NULL) {
                  productbsp = BioseqFind (sip);
                  if (productbsp != NULL) {
                    ValNodeAddPointer (&(rfp->bsplist), 0, (Pointer) productbsp);
                  }
                }
              }
            }
            *(prevsfp) = sfp->next;
            sfp->next = NULL;
            SeqFeatFree (sfp);
          } else {
            prevsfp = (Pointer PNTR) &(sfp->next);
          }
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
  AsnIoClose (aip);
}

static void RemoveDescriptorCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  AsnExpOptPtr   aeop;
  AsnIoPtr       aip;
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  ValNodePtr     nextsdp;
  Boolean        notext;
  ObjMgrTypePtr  omtp;
  Pointer PNTR   prevsdp;
  RemoveFormPtr  rfp;
  ValNodePtr     sdp;
  Uint2          subtype;

  if (mydata == NULL) return;
  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  rfp = (RemoveFormPtr) mydata;
  if (rfp == NULL) return;
  omtp = rfp->omtp;
  if (omtp == NULL || omtp->subtypefunc == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sdp = bsp->descr;
    prevsdp = (Pointer PNTR) &(bsp->descr);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sdp = bssp->descr;
    prevsdp = (Pointer PNTR) &(bssp->descr);
  } else return;
  GetTitle (rfp->findthis, rfp->findStr, sizeof (rfp->findStr) - 1);
  notext = StringHasNoText (rfp->findStr);
  aip = AsnIoNullOpen ();
  aeop = AsnExpOptNew (aip, NULL, NULL, AsnWriteRemoveForDCallBack);
  if (aeop != NULL) {
    aeop->user_data = (Pointer) rfp;
  }
  while (sdp != NULL) {
    nextsdp = sdp->next;
    subtype = (*(omtp->subtypefunc)) ((Pointer) sdp);
    if (subtype == rfp->subtype) {
      if (notext || ObjectHasSubstring (omtp, aip, (Pointer) sdp, rfp)) {
        *(prevsdp) = sdp->next;
        sdp->next = NULL;
        SeqDescFree (sdp);
      } else {
        prevsdp = (Pointer PNTR) &(sdp->next);
      }
    } else {
      prevsdp = (Pointer PNTR) &(sdp->next);
    }
    sdp = nextsdp;
  }
  AsnIoClose (aip);
}

static void DoRemoveAsnObject (ButtoN b)

{
  MsgAnswer      ans;
  BioseqPtr      bsp;
  Uint2          itemID;
  OMProcControl  ompc;
  RemoveFormPtr  rfp;
  SeqEntryPtr    sep;
  ValNodePtr     tmp;
  Int2           val;
  ValNodePtr     vnp;

  rfp = GetObjectExtra (b);
  if (rfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (rfp->input_entityID);
  if (sep == NULL) return;
  Hide (rfp->form);
  WatchCursor ();
  Update ();
  if (rfp->is_feature) {
    rfp->itemtype = OBJ_SEQFEAT;
  } else {
    rfp->itemtype = OBJ_SEQDESC;
  }
  vnp = NULL;
  val = GetValue (rfp->objlist);
  if (val > 0) {
    vnp = rfp->head;
    while (vnp != NULL && val > 1) {
      val--;
      vnp = vnp->next;
    }
  }
  if (vnp != NULL) {
    rfp->omp = ObjMgrGet ();
    rfp->omtp = NULL;
    if (rfp->omp != NULL) {
      rfp->omtp = ObjMgrTypeFind (rfp->omp, rfp->itemtype, NULL, NULL);
    }
    rfp->subtype = vnp->choice;
    if (rfp->itemtype != 0 && rfp->subtype != 0 && rfp->omtp != NULL) {
      if (rfp->is_feature) {
        SeqEntryExplore (sep, (Pointer) rfp, RemoveFeatureCallback);
        if (rfp->bsplist != NULL) {
          ans = Message (MSG_YN, "Remove protein products?");
          if (ans == ANS_YES) {
            for (tmp = rfp->bsplist; tmp != NULL; tmp = tmp->next) {
              bsp = (BioseqPtr) tmp->data.ptrvalue;
              itemID = GetItemIDGivenPointer (rfp->input_entityID, OBJ_BIOSEQ, (Pointer) bsp);
              if (itemID > 0) {
                MemSet ((Pointer) (&ompc), 0, sizeof (OMProcControl));
                ompc.do_not_reload_from_cache = TRUE;
                ompc.input_entityID = rfp->input_entityID;
                ompc.input_itemID = itemID;
                ompc.input_itemtype = OBJ_BIOSEQ;
                if (! DetachDataForProc (&ompc, FALSE)) {
                  Message (MSG_POSTERR, "DetachDataForProc failed");
                }
              }
            }
          }
        }
      } else {
        SeqEntryExplore (sep, (Pointer) rfp, RemoveDescriptorCallback);
      }
    }
  }
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (rfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, rfp->input_entityID, 0, 0);
  ObjMgrDeSelect (0, 0, 0, 0, NULL);
  Remove (rfp->form);
}

int LIBCALLBACK SortByVnpChoice (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr   vnp1;
  ValNodePtr   vnp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      if (vnp1->choice > vnp2->choice) {
        return 1;
      } else if (vnp1->choice < vnp2->choice) {
        return -1;
      } else {
        return 0;
      }
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}

static void RemoveMessageProc (ForM f, Int2 mssg)

{
  RemoveFormPtr  rfp;

  rfp = (RemoveFormPtr) GetObjectExtra (f);
  if (rfp != NULL) {
    if (rfp->appmessage != NULL) {
      rfp->appmessage (f, mssg);
    }
  }
}

static void CleanupRemovePage (GraphiC g, VoidPtr data)

{
  RemoveFormPtr  rfp;

  rfp = (RemoveFormPtr) data;
  if (rfp != NULL) {
    ValNodeFreeData (rfp->head);
    ValNodeFree (rfp->bsplist);
  }
  StdCleanupFormProc (g, data);
}

static CharPtr descNames [] = {
  " ", "MolType", "Modif", "Method", "Name",
  "Title", "Organism", "Comment", "Numbering",
  "MapLoc", "PIR", "GenBank", "Publication",
  "Region", "User", "SWISS-PROT", "dbXREF",
  "EMBL", "Create Date", "Update Date", "PRF",
  "PDB", "Heterogen", "BioSource", "MolInfo", NULL
};

/*
#ifdef INTERNAL_NCBI_SEQUIN
#define LISTHEIGHT 16
#else
#define LISTHEIGHT 8
#endif
*/

static void RemoveAsnObject (IteM i, Boolean feature)

{
  BaseFormPtr        bfp;
  ButtoN             b;
  GrouP              c;
  FeatDefPtr         curr;
  GrouP              g;
  GrouP              h;
  ValNodePtr         head;
  Int2               j;
  GrouP              k;
  Uint1              key;
  CharPtr            label = NULL;
  Int2               listHeight;
  RemoveFormPtr      rfp;
  SeqEntryPtr        sep;
  StdEditorProcsPtr  sepp;
  Char               str [256];
  Uint2              subtype;
  CharPtr            title;
  ValNodePtr         vnp;
  WindoW             w;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  rfp = (RemoveFormPtr) MemNew (sizeof (RemoveFormData));
  if (rfp == NULL) return;
  if (feature) {
    title = "Feature Removal";
  } else {
    title = "Descriptor Removal";
  }
  w = FixedWindow (-50, -33, -10, -10, title, StdCloseWindowProc);
  SetObjectExtra (w, rfp, CleanupRemovePage);
  rfp->form = (ForM) w;
  rfp->formmessage = RemoveMessageProc;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    rfp->appmessage = sepp->handleMessages;
  }

  rfp->input_entityID = bfp->input_entityID;
  rfp->input_itemID = bfp->input_itemID;
  rfp->input_itemtype = bfp->input_itemtype;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  g = HiddenGroup (h, 0, 2, NULL);
  rfp->is_feature = feature;
  if (feature) {
    StaticPrompt (g, "Feature", 0, 0, programFont, 'c');
  } else {
    StaticPrompt (g, "Descriptor", 0, 0, programFont, 'c');
  }
  if (indexerVersion) {
    listHeight = 16;
  } else {
    listHeight = 8;
  }
  rfp->objlist = SingleList (g, 16, listHeight, NULL);
  head = NULL;
  if (feature) {
    curr = FeatDefFindNext (NULL, &key, &label, FEATDEF_ANY, TRUE);
    while (curr != NULL) {
      if (key != FEATDEF_BAD) {
        subtype = curr->featdef_key;
        vnp = ValNodeNew (head);
        if (head == NULL) {
          head = vnp;
        }
        if (vnp != NULL) {
          vnp->choice = subtype;
          if (subtype != FEATDEF_misc_RNA &&
              subtype != FEATDEF_precursor_RNA &&
              subtype != FEATDEF_mat_peptide &&
              subtype != FEATDEF_sig_peptide &&
              subtype != FEATDEF_transit_peptide) {
            vnp->data.ptrvalue = StringSave (curr->typelabel);
          } else {
            StringNCpy_0 (str, curr->typelabel, sizeof (str) - 10);
            StringCat (str, "_imp");
            vnp->data.ptrvalue = StringSave (str);
          }
        }
      }
      curr = FeatDefFindNext (curr, &key, &label, FEATDEF_ANY, TRUE);
    }
  } else {
    for (j = 1; descNames [j] != NULL; j++) {
      vnp = ValNodeNew (head);
      if (head == NULL) {
        head = vnp;
      }
      if (vnp != NULL) {
        vnp->choice = j;
        vnp->data.ptrvalue = StringSave (descNames [j]);
      }
    }
  }
  if (head != NULL) {
    head = SortValNode (head, SortByVnpChoice);
    for (vnp = head; vnp != NULL; vnp = vnp->next) {
      ListItem (rfp->objlist, (CharPtr) vnp->data.ptrvalue);
    }
  }
  rfp->head = head;
  rfp->bsplist = NULL;

  k = HiddenGroup (h, 0, 2, NULL);
  StaticPrompt (k, "Optional string constraint", 0, dialogTextHeight, programFont, 'c');
  rfp->findthis = DialogText (k, "", 14, NULL);

  c = HiddenGroup (h, 4, 0, NULL);
  b = PushButton (c, "Accept", DoRemoveAsnObject);
  SetObjectExtra (b, rfp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) k, (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Update ();
}

extern void RemoveFeature (IteM i)

{
  RemoveAsnObject (i, TRUE);
}

extern void RemoveDescriptor (IteM i)

{
  RemoveAsnObject (i, FALSE);
}

#define SLCT_FEAT    1
#define SLCT_DESC    2
#define SLCT_BIOSEQ  3

typedef struct selectformdata {
  FEATURE_FORM_BLOCK

  Int2           type;
  LisT           objlist;
  TexT           findthis;
  Uint2          itemtype;
  Uint2          subtype;
  ObjMgrPtr      omp;
  ObjMgrTypePtr  omtp;
  ValNodePtr     head;
  Boolean        stringfound;
  Char           findStr [128];
} SelectFormData, PNTR SelectFormPtr;

static Boolean SelectObjCallback (GatherContextPtr gcp)

{
  ObjMgrTypePtr  omtp;
  SelectFormPtr  selfp;
  Uint2          subtype;

  if (gcp == NULL) return TRUE;
  selfp = (SelectFormPtr) gcp->userdata;
  if (selfp == NULL) return TRUE;
  if (gcp->thistype != selfp->itemtype) return TRUE;
  omtp = selfp->omtp;
  if (omtp == NULL) return TRUE;
  subtype = (*(omtp->subtypefunc)) ((Pointer) gcp->thisitem);
  if (subtype != selfp->subtype) return TRUE;
  ObjMgrAlsoSelect (gcp->entityID, gcp->itemID, gcp->thistype, 0, NULL);
  return TRUE;
}

static void DoSelectAsnObject (ButtoN b)

{
  GatherScope    gs;
  SelectFormPtr  selfp;
  SeqEntryPtr    sep;
  Int2           val;
  ValNodePtr     vnp;

  selfp = GetObjectExtra (b);
  if (selfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (selfp->input_entityID);
  if (sep == NULL) return;
  Hide (selfp->form);
  switch (selfp->type) {
    case SLCT_FEAT :
      selfp->itemtype = OBJ_SEQFEAT;
      break;
    case SLCT_DESC :
      selfp->itemtype = OBJ_SEQDESC;
      break;
    case SLCT_BIOSEQ :
      selfp->itemtype = OBJ_BIOSEQ;
      break;
    default :
      Remove (selfp->form);
      Update ();
      return;
  }
  WatchCursor ();
  Update ();
  vnp = NULL;
  val = GetValue (selfp->objlist);
  if (val > 0) {
    vnp = selfp->head;
    while (vnp != NULL && val > 1) {
      val--;
      vnp = vnp->next;
    }
  }
  if (selfp->type == SLCT_BIOSEQ || vnp != NULL) {
    selfp->omp = ObjMgrGet ();
    selfp->omtp = NULL;
    if (selfp->omp != NULL) {
      selfp->omtp = ObjMgrTypeFind (selfp->omp, selfp->itemtype, NULL, NULL);
    }
    if (selfp->type == SLCT_BIOSEQ) {
      selfp->subtype = Seq_repr_raw;
    } else {
      selfp->subtype = vnp->choice;
    }
    if (selfp->itemtype != 0 && selfp->subtype != 0 && selfp->omtp != NULL) {
      MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
      gs.seglevels = 1;
      MemSet((Pointer)(gs.ignore), (int)(TRUE), (size_t)(OBJ_MAX * sizeof(Boolean)));
      gs.ignore[OBJ_BIOSEQ] = FALSE;
      gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
      gs.ignore[OBJ_SEQANNOT] = FALSE;
      gs.ignore[OBJ_SEQFEAT] = FALSE;
      gs.ignore[OBJ_SEQDESC] = FALSE;
      GatherEntity (selfp->input_entityID, (Pointer) selfp, SelectObjCallback, &gs);
    }
  }
  ArrowCursor ();
  Update ();
  /* ObjMgrSendMsg (OM_MSG_UPDATE, selfp->input_entityID, 0, 0); */
  Remove (selfp->form);
}

static void SelectMessageProc (ForM f, Int2 mssg)

{
  SelectFormPtr  selfp;

  selfp = (SelectFormPtr) GetObjectExtra (f);
  if (selfp != NULL) {
    if (selfp->appmessage != NULL) {
      selfp->appmessage (f, mssg);
    }
  }
}

static void CleanupSelectPage (GraphiC g, VoidPtr data)

{
  SelectFormPtr  selfp;

  selfp = (SelectFormPtr) data;
  if (selfp != NULL) {
    ValNodeFreeData (selfp->head);
  }
  StdCleanupFormProc (g, data);
}

static void SelectAsnObject (IteM i, Int2 type)

{
  BaseFormPtr        bfp;
  ButtoN             b;
  GrouP              c;
  FeatDefPtr         curr;
  GrouP              g;
  GrouP              h;
  ValNodePtr         head;
  Int2               j;
  Uint1              key;
  CharPtr            label = NULL;
  Int2               listHeight;
  SelectFormPtr      selfp;
  SeqEntryPtr        sep;
  StdEditorProcsPtr  sepp;
  Uint2              subtype;
  CharPtr            title;
  ValNodePtr         vnp;
  WindoW             w;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  selfp = (SelectFormPtr) MemNew (sizeof (SelectFormData));
  if (selfp == NULL) return;
  switch (type) {
    case SLCT_FEAT :
      title = "Feature Selection";
      break;
    case SLCT_DESC :
      title = "Descriptor Selection";
      break;
    case SLCT_BIOSEQ :
      title = "Sequence Selection";
      break;
    default :
      title = "? Selection";
      break;
  }
  w = FixedWindow (-50, -33, -10, -10, title, StdCloseWindowProc);
  SetObjectExtra (w, selfp, CleanupSelectPage);
  selfp->form = (ForM) w;
  selfp->formmessage = SelectMessageProc;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    selfp->appmessage = sepp->handleMessages;
  }

  selfp->input_entityID = bfp->input_entityID;
  selfp->input_itemID = bfp->input_itemID;
  selfp->input_itemtype = bfp->input_itemtype;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  g = HiddenGroup (h, 0, 2, NULL);
  selfp->type = type;
  switch (type) {
    case SLCT_FEAT :
      StaticPrompt (g, "Feature", 0, 0, programFont, 'c');
      break;
    case SLCT_DESC :
      StaticPrompt (g, "Descriptor", 0, 0, programFont, 'c');
      break;
    case SLCT_BIOSEQ :
      StaticPrompt (g, "Sequence", 0, 0, programFont, 'c');
      break;
    default :
      break;
  }
  if (indexerVersion) {
    listHeight = 16;
  } else {
    listHeight = 8;
  }
  selfp->objlist = SingleList (g, 16, listHeight, NULL);
  head = NULL;
  if (type == SLCT_FEAT) {
    curr = FeatDefFindNext (NULL, &key, &label, FEATDEF_ANY, TRUE);
    while (curr != NULL) {
      if (key != FEATDEF_BAD) {
        subtype = curr->featdef_key;
        vnp = ValNodeNew (head);
        if (head == NULL) {
          head = vnp;
        }
        if (vnp != NULL) {
          vnp->choice = subtype;
          vnp->data.ptrvalue = StringSave (curr->typelabel);
        }
      }
      curr = FeatDefFindNext (curr, &key, &label, FEATDEF_ANY, TRUE);
    }
  } else if (type == SLCT_DESC) {
    for (j = 1; descNames [j] != NULL; j++) {
      vnp = ValNodeNew (head);
      if (head == NULL) {
        head = vnp;
      }
      if (vnp != NULL) {
        vnp->choice = j;
        vnp->data.ptrvalue = StringSave (descNames [j]);
      }
    }
  }
  if (head != NULL) {
    head = SortValNode (head, SortByVnpChoice);
    for (vnp = head; vnp != NULL; vnp = vnp->next) {
      ListItem (selfp->objlist, (CharPtr) vnp->data.ptrvalue);
    }
  }
  selfp->head = head;

  c = HiddenGroup (h, 4, 0, NULL);
  b = PushButton (c, "Accept", DoSelectAsnObject);
  SetObjectExtra (b, selfp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
  RealizeWindow (w);
  if (type == SLCT_BIOSEQ) {
    DoSelectAsnObject (b);
    Update ();
    return;
  }
  Show (w);
  Update ();
}

extern void SelectFeature (IteM i)

{
  SelectAsnObject (i, SLCT_FEAT);
}

extern void SelectDescriptor (IteM i)

{
  SelectAsnObject (i, SLCT_DESC);
}

extern void SelectBioseq (IteM i)

{
  SelectAsnObject (i, SLCT_BIOSEQ);
}

typedef struct fuseformdata {
  FEATURE_FORM_BLOCK

  LisT           objlist;
  Uint2          subtype;
  ObjMgrPtr      omp;
  ObjMgrTypePtr  omtp;
  ValNodePtr     head;
} FuseFormData, PNTR FuseFormPtr;

static void FuseFeatureCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr      bsp;
  FuseFormPtr    ffp;
  SeqFeatPtr     first;
  SeqFeatPtr     nextsfp;
  ObjMgrTypePtr  omtp;
  Pointer PNTR   prevsfp;
  SeqAnnotPtr    sap;
  SeqFeatPtr     sfp;
  SeqLocPtr      slp;
  Uint2          subtype;
  BioseqPtr      target;

  if (mydata == NULL) return;
  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  ffp = (FuseFormPtr) mydata;
  if (ffp == NULL) return;
  omtp = ffp->omtp;
  if (omtp == NULL || omtp->subtypefunc == NULL) return;
  if (! IS_Bioseq (sep)) return;
  bsp = (BioseqPtr) sep->data.ptrvalue;
  sap = bsp->annot;
  first = NULL;
  while (sap != NULL) {
     if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      prevsfp = (Pointer PNTR) &(sap->data);
      while (sfp != NULL) {
        nextsfp = sfp->next;
        subtype = (*(omtp->subtypefunc)) ((Pointer) sfp);
        if (subtype == ffp->subtype ||
           (ffp->subtype == FEATDEF_IMP &&
            subtype >= FEATDEF_allele && subtype <= FEATDEF_site_ref)) {
          if (first != NULL) {
            target = GetBioseqGivenSeqLoc (sfp->location, ffp->input_entityID);
            if (target != NULL) {
              slp = SeqLocMerge (target, sfp->location, first->location, FALSE, TRUE, FALSE);
              first->location = SeqLocFree (first->location);
              first->location = slp;
              first->partial = CheckSeqLocForPartial (slp, NULL, NULL);
            }
            *(prevsfp) = sfp->next;
            sfp->next = NULL;
            SeqFeatFree (sfp);
          } else {
            first = sfp;
            prevsfp = (Pointer PNTR) &(sfp->next);
          }
        }
        sfp = nextsfp;
      }
    }
    sap = sap->next;
  }
}

static void DoFuseFeature (ButtoN b)

{
  FuseFormPtr  ffp;
  SeqEntryPtr  sep;
  Int2         val;
  ValNodePtr   vnp;

  ffp = (FuseFormPtr) GetObjectExtra (b);
  if (ffp == NULL) return;
  sep = GetTopSeqEntryForEntityID (ffp->input_entityID);
  if (sep == NULL) return;
  Hide (ffp->form);
  WatchCursor ();
  Update ();

  vnp = NULL;
  val = GetValue (ffp->objlist);
  if (val > 0) {
    vnp = ffp->head;
    while (vnp != NULL && val > 1) {
      val--;
      vnp = vnp->next;
    }
  }
  if (vnp != NULL) {
    ffp->omp = ObjMgrGet ();
    ffp->omtp = NULL;
    if (ffp->omp != NULL) {
      ffp->omtp = ObjMgrTypeFind (ffp->omp, OBJ_SEQFEAT, NULL, NULL);
    }
    ffp->subtype = vnp->choice;
    if (ffp->subtype != 0 && ffp->omtp != NULL) {
      SeqEntryExplore (sep, (Pointer) ffp, FuseFeatureCallback);
    }
  }

  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (ffp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, ffp->input_entityID, 0, 0);
  ObjMgrDeSelect (0, 0, 0, 0, NULL);
  Remove (ffp->form);
}

static void FuseMessageProc (ForM f, Int2 mssg)

{
  FuseFormPtr  ffp;

  ffp = (FuseFormPtr) GetObjectExtra (f);
  if (ffp != NULL) {
    if (ffp->appmessage != NULL) {
      ffp->appmessage (f, mssg);
    }
  }
}

static void CleanupFusePage (GraphiC g, VoidPtr data)

{
  FuseFormPtr  ffp;

  ffp = (FuseFormPtr) data;
  if (ffp != NULL) {
    ValNodeFreeData (ffp->head);
  }
  StdCleanupFormProc (g, data);
}

extern void FuseFeature (IteM i)

{
  BaseFormPtr        bfp;
  ButtoN             b;
  GrouP              c;
  FeatDefPtr         curr;
  FuseFormPtr        ffp;
  GrouP              g;
  GrouP              h;
  ValNodePtr         head;
  Uint1              key;
  CharPtr            label = NULL;
  Int2               listHeight;
  SeqEntryPtr        sep;
  StdEditorProcsPtr  sepp;
  Uint2              subtype;
  ValNodePtr         vnp;
  WindoW             w;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  ffp = (FuseFormPtr) MemNew (sizeof (FuseFormData));
  if (ffp == NULL) return;
  w = FixedWindow (-50, -33, -10, -10, "Fuse Feature", StdCloseWindowProc);
  SetObjectExtra (w, ffp, CleanupFusePage);
  ffp->form = (ForM) w;
  ffp->formmessage = FuseMessageProc;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    ffp->appmessage = sepp->handleMessages;
  }

  ffp->input_entityID = bfp->input_entityID;
  ffp->input_itemID = bfp->input_itemID;
  ffp->input_itemtype = bfp->input_itemtype;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  g = HiddenGroup (h, 0, 2, NULL);
  StaticPrompt (g, "Feature", 0, 0, programFont, 'c');
  if (indexerVersion) {
    listHeight = 16;
  } else {
    listHeight = 8;
  }
  ffp->objlist = SingleList (g, 16, listHeight, NULL);
  head = NULL;
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
      ListItem (ffp->objlist, (CharPtr) vnp->data.ptrvalue);
    }
  }
  ffp->head = head;

  c = HiddenGroup (h, 4, 0, NULL);
  b = PushButton (c, "Accept", DoFuseFeature);
  SetObjectExtra (b, ffp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Update ();
}


extern Int2 LIBCALLBACK RefGeneUserGenFunc (Pointer data);

#define REFGENE_ASSEMBLY   1
#define REFGENE_RELATED    2
#define REFGENE_SPLICEVAR  3
#define REFGENE_RELDREK    4
#define REFGENE_REJECT     5
#define REFGENE_UNKNOWN    6

typedef struct refgeneuserdialog {
  DIALOG_MESSAGE_BLOCK
  GrouP         status;
  DialoG        fields;
} RefgeneUserDialog, PNTR RefgeneUserDialogPtr;

typedef struct refgeneuserform {
  FEATURE_FORM_BLOCK
  SeqEntryPtr   sep;
} RefgeneUserForm, PNTR RefgeneUserFormPtr;

static ENUM_ALIST(changeflags_alist)
  {" ",           0},
  {"Sequence",    1},
  {"Annotation",  2},
  {"Both",        3},
END_ENUM_ALIST

static ENUM_ALIST(refgene_alist)
  {" ",          0},
  {"Assembly",    REFGENE_ASSEMBLY},
  {"Related",     REFGENE_RELATED},
  {"SpliceVar",   REFGENE_SPLICEVAR},
  {"RelatedDrek", REFGENE_RELDREK},
  {"Reject",      REFGENE_REJECT},
  {"Unknown",     REFGENE_UNKNOWN},
END_ENUM_ALIST

static Uint2 refgene_types [] = {
  TAGLIST_TEXT, TAGLIST_TEXT, TAGLIST_TEXT, TAGLIST_POPUP, TAGLIST_POPUP
};

static Uint2 refgene_widths [] = {
  9, 7, 15, 0, 0
};

static EnumFieldAssocPtr refgene_popups [] = {
  NULL, NULL, NULL, changeflags_alist, refgene_alist
};

static CharPtr refgene_labels [] = {
  "", "Assembly", "Related", "SpliceVar", "RelatedDrek", "Reject", "Unknown", NULL
};

static CharPtr refgene_fields [] = {
  "Accession", "GI", "Comment", "Change", "Type", NULL
};

static void AccessionUserFieldPtrToVisStringDialog (DialoG d, Pointer data)

{
  CharPtr       accession;
  Boolean       annotChange;
  CharPtr       comment;
  UserFieldPtr  curr;
  UserFieldPtr  entry;
  Int2          field;
  Int2          flags;
  Int4          from;
  Int4          gi;
  ValNodePtr    head;
  Int2          i;
  Int2          j;
  ObjectIdPtr   oip;
  Boolean       seqChange;
  CharPtr       str;
  TagListPtr    tlp;
  Int4          to;
  UserFieldPtr  ufp;
  ValNodePtr    vnp;

  tlp = (TagListPtr) GetObjectExtra (d);
  if (tlp == NULL) return;
  str = MemNew (sizeof (Char) * 1024);
  head = NULL;
  curr = (UserFieldPtr) data;
  while (curr != NULL) {
    oip = curr->label;
    if (oip != NULL) {
      field = 0;
      for (i = REFGENE_ASSEMBLY; i <= REFGENE_UNKNOWN; i++) {
        if (StringICmp (oip->str, refgene_labels [i]) == 0 && curr->choice == 11) {
          field = i;
        }
      }
      if (field > 0) {
        entry = (UserFieldPtr) curr->data.ptrvalue;
        while (entry != NULL && entry->choice == 11) {
          accession = NULL;
          comment = NULL;
          gi = 0;
          from = 0;
          to = 0;
          annotChange = FALSE;
          seqChange = FALSE;
          ufp = (UserFieldPtr) entry->data.ptrvalue;
          while (ufp != NULL) {
            oip = ufp->label;
            if (oip != NULL && oip->str != NULL) {
              if (StringICmp (oip->str, "accession") == 0 && ufp->choice == 1) {
                accession = (CharPtr) ufp->data.ptrvalue;
              } else if (StringICmp (oip->str, "gi") == 0 && ufp->choice == 2) {
                gi = ufp->data.intvalue;
              } else if (StringICmp (oip->str, "from") == 0 && ufp->choice == 2) {
                from = ufp->data.intvalue;
              } else if (StringICmp (oip->str, "to") == 0 && ufp->choice == 2) {
                to = ufp->data.intvalue;
              } else if (StringICmp (oip->str, "sequenceChange") == 0 && ufp->choice == 4) {
                seqChange = ufp->data.boolvalue;
              } else if (StringICmp (oip->str, "annotationChange") == 0 && ufp->choice == 4) {
                annotChange = ufp->data.boolvalue;
              } else if (StringICmp (oip->str, "comment") == 0 && ufp->choice == 1) {
                comment = (CharPtr) ufp->data.ptrvalue;
              }
            }
            ufp = ufp->next;
          }
          if (accession != NULL) {
            if (comment == NULL) {
              comment = "";
            }
            flags = 0;
            if (seqChange) {
              flags++;
            }
            if (annotChange) {
              flags += 2;
            }
            if (gi != 0) {
              sprintf (str, "%s\t%ld\t%s\t%d\t%d\n", accession, (long) gi, comment, (int) flags, (int) field);
            } else {
              sprintf (str, "%s\t\t%s\t%d\t%d\n", accession, comment, (int) flags, (int) field);
            }
            vnp = ValNodeNew (head);
            if (head == NULL) {
              head = vnp;
            }
            if (vnp != NULL) {
              vnp->data.ptrvalue = StringSave (str);
            }
          }
          entry = entry->next;
        }
      }
    }
    curr = curr->next;
  }
  MemFree (str);
  SendMessageToDialog (tlp->dialog, VIB_MSG_RESET);
  tlp->vnp = head;
  SendMessageToDialog (tlp->dialog, VIB_MSG_REDRAW);
  for (j = 0, vnp = tlp->vnp; vnp != NULL; j++, vnp = vnp->next) {
  }
  tlp->max = MAX ((Int2) 0, (Int2) (j - tlp->rows + 1));
  CorrectBarMax (tlp->bar, tlp->max);
  CorrectBarPage (tlp->bar, tlp->rows - 1, tlp->rows - 1);
}

static Pointer VisStringDialogToUserFieldPtr (DialoG d)

{
  return NULL;
}

static void UserObjectPtrToRefGeneDialog (DialoG d, Pointer data)

{
  UserFieldPtr          curr;
  ObjectIdPtr           oip;
  RefgeneUserDialogPtr  rdp;
  Int2                  status = 0;
  CharPtr               str;
  UserObjectPtr         uop;

  rdp = (RefgeneUserDialogPtr) GetObjectExtra (d);
  if (rdp == NULL) return;
  uop = (UserObjectPtr) data;
  if (uop == NULL || uop->type == NULL || StringICmp (uop->type->str, "RefGeneTracking") != 0) {
    SetValue (rdp->status, 0);
    PointerToDialog (rdp->fields, NULL);
    return;
  }
  PointerToDialog (rdp->fields, uop->data);
  for (curr = uop->data; curr != NULL; curr = curr->next) {
    oip = curr->label;
    if (oip != NULL && StringICmp (oip->str, "Status") == 0) {
      break;
    }
  }
  if (curr != NULL && curr->choice == 1) {
    str = (CharPtr) curr->data.ptrvalue;
    if (StringICmp (str, "Predicted") == 0) {
      status = 1;
    } else if (StringICmp (str, "Provisional") == 0) {
      status = 2;
    } else if (StringICmp (str, "Reviewed") == 0) {
      status = 3;
    }
  }
  SetValue (rdp->status, status);
}

static Pointer RefGeneDialogToUserObjectPtr (DialoG d)

{
  Boolean               annotChange;
  Char                  ch;
  Int2                  i;
  Int2                  j;
  size_t                len;
  Int4                  num [5];
  Boolean               okay;
  CharPtr               organism = NULL;
  RefgeneUserDialogPtr  rdp;
  Boolean               seqChange;
  Int2                  status;
  CharPtr               str;
  TagListPtr            tlp;
  CharPtr               txt [5];
  UserObjectPtr         uop;
  long int              val;
  ValNodePtr            vnp;

  rdp = (RefgeneUserDialogPtr) GetObjectExtra (d);
  if (rdp == NULL) return NULL;

  uop = CreateRefGeneTrackUserObject ();
  if (uop == NULL) return NULL;

  status = GetValue (rdp->status);
  if (status == 1) {
    AddStatusToRefGeneTrackUserObject (uop, "Predicted");
  } else if (status == 2) {
    AddStatusToRefGeneTrackUserObject (uop, "Provisional");
  } else if (status == 3) {
    AddStatusToRefGeneTrackUserObject (uop, "Reviewed");
  }

  tlp = (TagListPtr) GetObjectExtra (rdp->fields);
  if (tlp != NULL && tlp->vnp != NULL) {
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
        for (j = 0; j < 5; j++) {
          txt [j] = ExtractTagListColumn ((CharPtr) vnp->data.ptrvalue, j);
          num [j] = 0;
        }
        for (j = 1; j < 5; j++) {
          if (j != 2) {
            num [j] = 0;
            if (txt [j] != NULL && sscanf (txt [j], "%ld", &val) == 1) {
              num [j] = val;
            }
          }
        }
        annotChange = FALSE;
        seqChange = FALSE;
        if (num [3] >= 2) {
          annotChange = TRUE;
          (num [3]) -= 2;
        }
        if (num [3] > 0) {
          seqChange = TRUE;
        }
        i = num [4];
        if (i >= REFGENE_ASSEMBLY && i <= REFGENE_UNKNOWN) {
          if (! StringHasNoText (txt [0])) {
            AddAccessionToRefGeneTrackUserObject (uop, refgene_labels [i],
                                                  txt [0], num [1],
                                                  seqChange, annotChange, txt [2]);
          }
          for (j = 0; j < 4; j++) {
            txt [j] = MemFree (txt [j]);
          }
        }
      }
    }
  }

  return uop;
}

static DialoG CreateRefGeneDialog (GrouP g)

{
  Int2                  i;
  PrompT                lastppt;
  GrouP                 p;
  PrompT                ppt;
  GrouP                 q;
  RefgeneUserDialogPtr  rdp;
  TagListPtr            tlp;
  GrouP                 x;

  p = HiddenGroup (g, -1, 0, NULL);
  SetGroupSpacing (p, 10, 10);

  rdp = (RefgeneUserDialogPtr) MemNew (sizeof (RefgeneUserDialog));
  if (rdp == NULL) return NULL;

  SetObjectExtra (p, rdp, NULL);
  rdp->dialog = (DialoG) p;
  rdp->todialog = UserObjectPtrToRefGeneDialog;
  rdp->fromdialog = RefGeneDialogToUserObjectPtr;

  x = HiddenGroup (p, 4, 0, NULL);
  /* StaticPrompt (x, "Status", 0, stdLineHeight, programFont, 'l'); */
  rdp->status = HiddenGroup (x, 4, 0, NULL);
  SetObjectExtra (rdp->status, rdp, NULL);
  RadioButton (rdp->status, "Predicted");
  RadioButton (rdp->status, "Provisional");
  RadioButton (rdp->status, "Reviewed");

  q = HiddenGroup (p, -6, 0, NULL);
  lastppt = NULL;
  ppt = NULL;
  for (i = 0; i < 5; i++) {
    lastppt = ppt;
    ppt = StaticPrompt (q, refgene_fields [i], refgene_widths [i] * stdCharWidth, 0, systemFont, 'c');
  }
  rdp->fields = CreateTagListDialog (p, 6, 5, STD_TAG_SPACING,
                                     refgene_types, refgene_widths, refgene_popups,
                                     AccessionUserFieldPtrToVisStringDialog,
                                     VisStringDialogToUserFieldPtr);

  tlp = (TagListPtr) GetObjectExtra (rdp->fields);
  if (tlp != NULL) {
    AlignObjects (ALIGN_JUSTIFY, (HANDLE) tlp->control [3], (HANDLE) lastppt, NULL);
    AlignObjects (ALIGN_JUSTIFY, (HANDLE) tlp->control [4], (HANDLE) ppt, NULL);
  }

  AlignObjects (ALIGN_CENTER, (HANDLE) x, (HANDLE) q, (HANDLE) rdp->fields, NULL);
  return (DialoG) p;
}

static void RefgeneUserFormMessage (ForM f, Int2 mssg)

{
  RefgeneUserFormPtr  rfp;

  rfp = (RefgeneUserFormPtr) GetObjectExtra (f);
  if (rfp != NULL) {
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
        if (rfp->appmessage != NULL) {
          rfp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static ForM CreateRefGeneDescForm (Int2 left, Int2 top, Int2 width,
                                   Int2 height, CharPtr title, ValNodePtr sdp,
                                   SeqEntryPtr sep, FormActnFunc actproc)

{
  ButtoN              b;
  GrouP               c;
  GrouP               g;
  RefgeneUserFormPtr  rfp;
  StdEditorProcsPtr   sepp;
  WindoW              w;

  w = NULL;
  rfp = (RefgeneUserFormPtr) MemNew (sizeof (RefgeneUserForm));
  if (rfp != NULL) {
    w = FixedWindow (left, top, width, height, title, StdCloseWindowProc);
    SetObjectExtra (w, rfp, StdDescFormCleanupProc);
    rfp->form = (ForM) w;
    rfp->actproc = actproc;
    rfp->formmessage = RefgeneUserFormMessage;

    rfp->sep = sep;

#ifndef WIN_MAC
    CreateStdEditorFormMenus (w);
#endif
    sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
    if (sepp != NULL) {
      SetActivate (w, sepp->activateForm);
      rfp->appmessage = sepp->handleMessages;
    }

    g = HiddenGroup (w, -1, 0, NULL);
    rfp->data = CreateRefGeneDialog (g);

    c = HiddenGroup (w, 2, 0, NULL);
    b = PushButton (c, "Accept", StdAcceptFormButtonProc);
    SetObjectExtra (b, rfp, NULL);
    PushButton (c, "Cancel", StdCancelButtonProc);
    AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
    RealizeWindow (w);
  }
  return (ForM) w;
}

extern Int2 LIBCALLBACK RefGeneUserGenFunc (Pointer data)

{
  ObjectIdPtr         oip;
  OMProcControlPtr    ompcp;
  OMUserDataPtr       omudp;
  ObjMgrProcPtr       proc;
  RefgeneUserFormPtr  rfp;
  ValNodePtr          sdp;
  SeqEntryPtr         sep;
  UserObjectPtr       uop;
  WindoW              w;

  ompcp = (OMProcControlPtr) data;
  w = NULL;
  sdp = NULL;
  sep = NULL;
  uop = NULL;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  proc = ompcp->proc;
  switch (ompcp->input_itemtype) {
    case OBJ_SEQDESC :
      sdp = (ValNodePtr) ompcp->input_data;
      if (sdp != NULL && sdp->choice != Seq_descr_user) {
        return OM_MSG_RET_ERROR;
      }
      uop = (UserObjectPtr) sdp->data.ptrvalue;
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
    if (StringCmp (proc->procname, "Edit RefGene UserTrack Desc") == 0) {
      rfp = (RefgeneUserFormPtr) omudp->userdata.ptrvalue;
      if (rfp != NULL) {
        Select (rfp->form);
      }
      return OM_MSG_RET_DONE;
    } else {
      return OM_MSG_RET_OK; /* not this type, check next registered user object editor */
    }
  }
  if (uop != NULL) {
    oip = uop->type;
    if (oip == NULL || oip->str == NULL) return OM_MSG_RET_OK;
    if (StringCmp (oip->str, "RefGeneTracking") != 0) return OM_MSG_RET_OK;
  }
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  w = (WindoW) CreateRefGeneDescForm (-50, -33, -10, -10,
                                      "Reference Gene Tracking", sdp, sep,
                                      StdDescFormActnProc);
  rfp = (RefgeneUserFormPtr) GetObjectExtra (w);
  if (rfp != NULL) {
    rfp->input_entityID = ompcp->input_entityID;
    rfp->input_itemID = ompcp->input_itemID;
    rfp->input_itemtype = ompcp->input_itemtype;
    rfp->this_itemtype = OBJ_SEQDESC;
    rfp->this_subtype = Seq_descr_user;
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
    if (sdp != NULL) {
      PointerToDialog (rfp->data, (Pointer) sdp->data.ptrvalue);
      SetClosestParentIfDuplicating ((BaseFormPtr) rfp);
    }
  }
  Show (w);
  Select (w);
  return OM_MSG_RET_DONE;
}

/*
static void TestGeneRefStuff (void)

{
  UserObjectPtr uop;
  ValNodePtr    sdp;

  uop = CreateRefGeneTrackUserObject ();
  AddAccessionToRefGeneTrackUserObject (uop, "Assembly", "U12345", 57, 29, 1995);
  AddAccessionToRefGeneTrackUserObject (uop, "Assembly", "L97531", 142, 66, 963);
  AddAccessionToRefGeneTrackUserObject (uop, "Assembly", "M66778", 823, 7677, 343);
  AddAccessionToRefGeneTrackUserObject (uop, "Related", "P34345", 445, 0, 0);
  AddAccessionToRefGeneTrackUserObject (uop, "Reject", "S19635", 1765, 0, 0);
  AddAccessionToRefGeneTrackUserObject (uop, "Related", "Q14884", 664, 35, 97);
  sdp = ValNodeNew (NULL);
  sdp->choice = Seq_descr_user;
  sdp->data.ptrvalue = (Pointer) uop;
  if (! ObjMgrRegister (OBJ_SEQDESC, (Pointer) sdp)) {
     ErrPostEx (SEV_ERROR, 0, 0, "ObjMgrRegister failed.");
  }
}
*/





typedef struct historyformdata {
  FEATURE_FORM_BLOCK

  BioseqPtr      bsp;
  DialoG         replace_date;
  DialoG         replace_ids;
  ButtoN         secondary_on_part;
  DialoG         replaced_by_date;
  DialoG         replaced_by_ids;
  ButtoN         deleted;
  DialoG         deleted_date;
} HistoryFormData, PNTR HistoryFormPtr;

static SeqIdPtr VisStrDialogToSeqIdSet (DialoG d)

{
  long          gi;
  SeqIdPtr      head = NULL;
  ValNodePtr    list;
  CharPtr       str;
  TextSeqIdPtr  tsip;
  ValNodePtr    vnp;

  if (d == NULL) return NULL;
  list = DialogToPointer (d);
  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (str != NULL) {
      if (sscanf (str, "%ld", &gi)) {
        /*
        ValNodeAddInt (&head, SEQID_GI, (Int4) gi);
        */
      } else {
        tsip = TextSeqIdNew ();
        if (tsip != NULL) {
          tsip->accession = StringSaveNoNull (str);
          ValNodeAddPointer (&head, SEQID_GENBANK, (Pointer) tsip);
        }
      }
    }
  }
  ValNodeFreeData (list);
  return head;
}

static int LIBCALLBACK SortByName (VoidPtr ptr1, VoidPtr ptr2)

{
  CharPtr     str1;
  CharPtr     str2;
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      str1 = (CharPtr) vnp1->data.ptrvalue;
      str2 = (CharPtr) vnp2->data.ptrvalue;
      if (str1 != NULL && str2 != NULL) {
        return StringICmp (str1, str2);
      } else {
        return 0;
      }
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}

static ValNodePtr GetStringsForSeqIDs (SeqIdPtr sip)

{
  Char          buf [40];
  ValNodePtr    head = NULL;
  TextSeqIdPtr  tsip;

  if (sip == NULL) return NULL;
  while (sip != NULL) {
    buf [0] = '\0';
    switch (sip->choice) {
      case SEQID_GENBANK :
      case SEQID_EMBL :
      case SEQID_DDBJ :
        tsip = (TextSeqIdPtr) sip->data.ptrvalue;
        if (tsip != NULL && (! StringHasNoText (tsip->accession))) {
          StringNCpy_0 (buf, tsip->accession, sizeof (buf));
        }
        break;
      case SEQID_GI :
        /*
        gi = sip->data.intvalue;
        if (gi > 0) {
          sprintf (buf, "%ld", (long) gi);
        }
        */
        break;
      default :
        break;
    }
    if (! StringHasNoText (buf)) {
      ValNodeCopyStr (&head, 0, buf);
    }
    sip = sip->next;
  }
  return head;
}

static void AddGenBankBlockToBioseq (BioseqPtr bsp, ValNodePtr head1, ValNodePtr head2)

{
  GBBlockPtr       gbp = NULL;
  CharPtr          last = NULL;
  ValNodePtr       next;
  ValNodePtr PNTR  prev;
  ValNodePtr       sdp;
  SeqEntryPtr      sep;
  CharPtr          str1, str2;
  ValNodePtr       vnp, vnp1, vnp2;

  sdp = BioseqGetSeqDescr (bsp, Seq_descr_genbank, NULL);
  if (sdp != NULL) {
    gbp = (GBBlockPtr) sdp->data.ptrvalue;
    if (gbp != NULL) {
      for (vnp1 = head1; vnp1 != NULL; vnp1 = vnp1->next) {
        str1 = (CharPtr) vnp1->data.ptrvalue;
        if (str1 != NULL) {
          for (vnp2 = gbp->extra_accessions; vnp2 != NULL; vnp2 = vnp2->next) {
            str2 = (CharPtr) vnp2->data.ptrvalue;
            if (StringICmp (str1, str2) == 0) {
              vnp2->data.ptrvalue = MemFree (vnp2->data.ptrvalue);
            }
          }
        }
      }
    }
  }
  if (sdp == NULL) {
    sep = SeqMgrGetSeqEntryForData (bsp);
    sdp = CreateNewDescriptor (sep, Seq_descr_genbank);
    if (sdp != NULL) {
      sdp->data.ptrvalue = GBBlockNew ();
    }
  }
  if (sdp != NULL) {
    gbp = (GBBlockPtr) sdp->data.ptrvalue;
    if (gbp != NULL) {
      while (head2 != NULL) {
        ValNodeCopyStr (&(gbp->extra_accessions), 0, (CharPtr) head2->data.ptrvalue);
        head2 = head2->next;
      }
      /*
      ValNodeLink (&(gbp->extra_accessions), head2);
      head2 = NULL;
      */
      gbp->extra_accessions = SortValNode (gbp->extra_accessions, SortByName);
      prev = &(gbp->extra_accessions);
      vnp = gbp->extra_accessions;
      last = NULL;
      while (vnp != NULL) {
        next = vnp->next;
        str2 = (CharPtr) vnp->data.ptrvalue;
        if (str2 == NULL || StringHasNoText (str2) || StringICmp (last, str2) == 0) {
          *prev = next;
          vnp->next = NULL;
          MemFree (vnp);
          vnp = next;
        } else {
          last = str2;
          prev = &(vnp->next);
          vnp = next;
        }
      }
    }
  }

}

static void DoChangeHistory (ButtoN b)

{
  MsgAnswer       ans;
  BioseqPtr       bsp;
  ValNodePtr      head1 = NULL, head2 = NULL;
  HistoryFormPtr  hfp;
  SeqHistPtr      hist;
  BioseqPtr       pbsp;
  SeqEntryPtr     sep;
  SeqLocPtr       slp;
  CharPtr         str1, str2;
  ValNodePtr      vnp1, vnp2;

  hfp = (HistoryFormPtr) GetObjectExtra (b);
  if (hfp == NULL) return;
  ans = Message (MSG_OKC, "Are you sure you want to edit the history?");
  if (ans == ANS_CANCEL) {
    return;
  }
  Hide (hfp->form);
  Update ();
  bsp = hfp->bsp;
  hist = bsp->hist;
  if (hist == NULL) {
    hist = SeqHistNew ();
    bsp->hist = hist;
  }
  if (hist != NULL) {

    hist->replace_date = DateFree (hist->replace_date);
    hist->replace_date = DialogToPointer (hfp->replace_date);
    head1 = GetStringsForSeqIDs (hist->replace_ids);
    hist->replace_ids = SeqIdSetFree (hist->replace_ids);
    hist->replace_ids = VisStrDialogToSeqIdSet (hfp->replace_ids);
    head2 = GetStringsForSeqIDs (hist->replace_ids);

    hist->replaced_by_date = DateFree (hist->replaced_by_date);
    hist->replaced_by_date = DialogToPointer (hfp->replaced_by_date);
    hist->replaced_by_ids = SeqIdSetFree (hist->replaced_by_ids);
    hist->replaced_by_ids = VisStrDialogToSeqIdSet (hfp->replaced_by_ids);

    hist->deleted = GetStatus (hfp->deleted);
    hist->deleted_date = DateFree (hist->deleted_date);
    hist->deleted_date = DialogToPointer (hfp->deleted_date);
  }

  if (hist->assembly == NULL &&
      hist->replace_date == NULL && hist->replace_ids == NULL &&
      hist->replaced_by_date == NULL && hist->replaced_by_ids == NULL &&
      (! hist->deleted) && hist->deleted_date == NULL) {
    bsp->hist = SeqHistFree (bsp->hist);
  }

  head1 = SortValNode (head1, SortByName);
  head2 = SortValNode (head2, SortByName);
  for (vnp1 = head1; vnp1 != NULL; vnp1 = vnp1->next) {
    str1 = (CharPtr) vnp1->data.ptrvalue;
    for (vnp2 = head2; vnp2 != NULL; vnp2 = vnp2->next) {
      str2 = (CharPtr) vnp2->data.ptrvalue;
      if (StringICmp (str1, str2) == 0) {
        vnp1->data.ptrvalue = MemFree (vnp1->data.ptrvalue);
      }
    }
  }

  AddGenBankBlockToBioseq (bsp, head1, head2);

  if (GetStatus (hfp->secondary_on_part)) {
    if (bsp->repr == Seq_repr_seg) {
      for (slp = (SeqLocPtr) bsp->seq_ext; slp != NULL; slp = slp->next) {
        pbsp = BioseqFind (SeqLocId (slp));
        if (pbsp != NULL) {
          AddGenBankBlockToBioseq (pbsp, head1, head2);
        }
      }
    }
  }

  ValNodeFreeData (head1);
  ValNodeFreeData (head2);

  sep = GetTopSeqEntryForEntityID (hfp->input_entityID);
  EntryCheckGBBlock (sep);

  Update ();
  ObjMgrSetDirtyFlag (hfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, hfp->input_entityID, 0, 0);
  ObjMgrDeSelect (0, 0, 0, 0, NULL);
  Remove (hfp->form);
}

static void SeqIdSetToVisStrDialog (DialoG d, SeqIdPtr sip)

{
  ValNodePtr  head = NULL;

  if (d == NULL || sip == NULL) return;
  head = GetStringsForSeqIDs (sip);
  PointerToDialog (d, head);
  ValNodeFreeData (head);
}

static void HistoryMessageProc (ForM f, Int2 mssg)

{
  HistoryFormPtr  hfp;

  hfp = (HistoryFormPtr) GetObjectExtra (f);
  if (hfp != NULL) {
    if (hfp->appmessage != NULL) {
      hfp->appmessage (f, mssg);
    }
  }
}

static void CleanupHistoryPage (GraphiC g, VoidPtr data)

{
  HistoryFormPtr  hfp;

  hfp = (HistoryFormPtr) data;
  if (hfp != NULL) {
  }
  StdCleanupFormProc (g, data);
}

extern void EditSequenceHistory (IteM i)

{
  ButtoN             b;
  BaseFormPtr        bfp;
  BioseqPtr          bsp;
  GrouP              c;
  GrouP              g;
  GrouP              h;
  HistoryFormPtr     hfp;
  SeqHistPtr         hist;
  GrouP              j;
  GrouP              k;
  PrompT             ppt1, ppt2, ppt3;
  SeqEntryPtr        sep;
  StdEditorProcsPtr  sepp;
  WindoW             w;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  bsp = GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID, bfp->input_itemtype);
  if (bsp == NULL) return;
  hist = bsp->hist;

  hfp = (HistoryFormPtr) MemNew (sizeof (HistoryFormData));
  if (hfp == NULL) return;
  w = FixedWindow (-50, -33, -10, -10, "Sequence History", StdCloseWindowProc);
  SetObjectExtra (w, hfp, CleanupHistoryPage);
  hfp->form = (ForM) w;
  hfp->formmessage = HistoryMessageProc;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    hfp->appmessage = sepp->handleMessages;
  }

  hfp->input_entityID = bfp->input_entityID;
  hfp->input_itemID = bfp->input_itemID;
  hfp->input_itemtype = bfp->input_itemtype;

  hfp->bsp = bsp;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  ppt1 = StaticPrompt (h, "Replaces", 0, 0, systemFont, 'c');

  g = HiddenGroup (h, -1, 0, NULL);
  hfp->replace_ids = CreateVisibleStringDialog (g, 4, -1, 10);
  hfp->secondary_on_part = CheckBox (g, "Secondary on Parts", NULL);
  hfp->replace_date = CreateDateDialog (g, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) hfp->replace_ids, (HANDLE) hfp->secondary_on_part, (HANDLE) hfp->replace_date, NULL);

  ppt2 = StaticPrompt (h, "Replaced By", 0, 0, systemFont, 'c');

  j = HiddenGroup (h, -1, 0, NULL);
  hfp->replaced_by_ids = CreateVisibleStringDialog (j, 4, -1, 10);
  hfp->replaced_by_date = CreateDateDialog (j, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) hfp->replaced_by_ids, (HANDLE) hfp->replaced_by_date, NULL);

  ppt3 = StaticPrompt (h, "Status", 0, 0, systemFont, 'c');

  k = HiddenGroup (h, -1, 0, NULL);
  hfp->deleted = CheckBox (k, "Deleted", NULL);
  hfp->deleted_date = CreateDateDialog (k, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) hfp->deleted, (HANDLE) hfp->deleted_date, NULL);

  c = HiddenGroup (h, 4, 0, NULL);
  b = PushButton (c, "Accept", DoChangeHistory);
  SetObjectExtra (b, hfp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) j, (HANDLE) k, (HANDLE) c,
                (HANDLE) ppt1, (HANDLE) ppt2, (HANDLE) ppt3, NULL);
  RealizeWindow (w);

  if (bsp->repr != Seq_repr_seg) {
    Disable (hfp->secondary_on_part);
  } else {
    SetStatus (hfp->secondary_on_part, TRUE);
  }

  if (hist != NULL) {
    PointerToDialog (hfp->replace_date, hist->replace_date);
    SeqIdSetToVisStrDialog (hfp->replace_ids, hist->replace_ids);
    PointerToDialog (hfp->replaced_by_date, hist->replaced_by_date);
    SeqIdSetToVisStrDialog (hfp->replaced_by_ids, hist->replaced_by_ids);
    PointerToDialog (hfp->deleted_date, hist->deleted_date);
    SetStatus (hfp->deleted, hist->deleted);
  }

  Show (w);
  Update ();
}

extern void HandleProjectAsn (ProjectPtr proj, Uint2 entityID)

{
  Int2              db = -1;
  EntrezGlobalsPtr  egp;
  Int4              i;
  ValNodePtr        list;
  Int4              num;
  ValNodePtr        pip;
  Int4Ptr           uids;
  ValNodePtr        vnp;

  if (proj == NULL) return;
  if (! useEntrez) return;
  egp = (EntrezGlobalsPtr) GetAppProperty ("EntrezGlobals");
  if (egp == NULL) return;
  pip = proj->data;
  if (pip == NULL) return;
  list = (ValNodePtr) pip->data.ptrvalue;
  if (list == NULL) return;
  if (pip->choice >= ProjectItem_protent && pip->choice <= ProjectItem_genomeent) {
    if (egp->retrieveProjectProc == NULL) return;
    if (! EntrezIsInited ()) {
      SequinEntrezInit ("Sequin", FALSE, NULL);
    }
    egp->retrieveProjectProc (NULL, (Pointer) proj);
    Update ();
    return;
  }
  if (egp->retrieveDocsProc == NULL) return;
  switch (pip->choice) {
    case ProjectItem_pmuid :
      db = 0;
      break;
    case ProjectItem_protuid :
      db = 1;
      break;
    case ProjectItem_nucuid :
      db = 2;
      break;
    case ProjectItem_genomeuid :
      db = 4;
      break;
    case ProjectItem_structuid :
      db = 3;
      break;
    default :
      break;
  }
  if (db == -1) return;
  if (! EntrezIsInited ()) {
    SequinEntrezInit ("Sequin", FALSE, NULL);
  }
  num = ValNodeLen (list);
  uids = MemNew ((size_t) (num * sizeof (Int4)));
  if (uids == NULL) return;
  for (i = 0, vnp = list; i < 32766 && vnp != NULL; i++, vnp = vnp->next) {
    uids [i] = vnp->data.intvalue;
  }
  if (egp->retrieveDocsProc != NULL) {
    egp->retrieveDocsProc (NULL, (Int2) num, 0, uids, db);
  }
  MemFree (uids);
  Update ();
}

/* BioseqViewOrDocSumChoice allows a single callback for each analysis item */
static Int2 BioseqViewOrDocSumChoice (NewObjectPtr nop)

{
  Int2  which = 0;   /* 1 = bioseq viewer, 2 = docsum window */

  if (nop == NULL) return 0;
#ifdef WIN_MAC
  if (bioseqViewUp) {
    which = 1;
  } else if (docSumUp) {
    which = 2;
  }
#else
  if (nop->bspOK) {
    which = 1;
  } else if (nop->dsmOK) {
    which = 2;
  }
#endif
  return which;
}

/*
#define TEST_APPLE_EVENT_MESSAGING
*/

#ifndef TEST_APPLE_EVENT_MESSAGING
static void AddRestrictionSite (SeqAnnotPtr annot, PackSeqPntPtr pspp, CharPtr name)

{
  DbtagPtr     dbt;
  ObjectIdPtr  oip;
  RsiteRefPtr  rrp;
  SeqFeatPtr   sfp, lastsfp;
  SeqLocPtr    slp;

  if (annot == NULL || pspp == NULL || name == NULL) return;
  slp = ValNodeNew (NULL);
  if (slp == NULL) return;
  slp->choice = SEQLOC_PACKED_PNT;
  slp->data.ptrvalue = (Pointer) pspp;
  sfp = SeqFeatNew ();
  if (sfp == NULL) return;

  sfp->data.choice = SEQFEAT_RSITE;
  sfp->location = slp;
  dbt = DbtagNew ();
  if (dbt != NULL) {
    dbt->db = StringSave ("REBASE");
    oip = ObjectIdNew ();
    if (oip != NULL) {
      oip->str = StringSave (name);
    }
    dbt->tag = oip;
  }
  rrp = ValNodeNew (NULL);
  if (rrp != NULL) {
    rrp->choice = 2;
    rrp->data.ptrvalue = dbt;
  }
  sfp->data.value.ptrvalue = (Pointer) rrp;

  if (annot->data == NULL) {
    annot->data = (Pointer) sfp;
  } else {
    lastsfp = (SeqFeatPtr) annot->data;
    while (lastsfp->next != NULL) {
      lastsfp = lastsfp->next;
    }
    lastsfp->next = sfp;
  }
}

static void RestrictionSearchProc (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  SeqAnnotPtr    annot;
  BioseqPtr      bsp;
  ComPatPtr      cpp, cpph;
  ValNodePtr     desc;
  SeqAnnotPtr    lastannot;
  PackSeqPntPtr  pspp;
  Int4           pt;
  SeqAlignPtr    sap;
  SeqLocPtr      slp, slpn;

  if (sep == NULL) return;
  if (! IS_Bioseq (sep)) return;
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL) return;
  if (! ISA_na (bsp->mol)) return;

  desc = ValNodeNew (NULL);
  desc->choice = Annot_descr_name;
  desc->data.ptrvalue = StringSave ("cutsites");

  annot = SeqAnnotNew ();
  annot->type = 1;
  annot->desc = desc;
  annot->data = NULL;

  cpph = (ComPatPtr) mydata;
  cpp = cpph;
  while (cpp != NULL) {
    sap = PatternMatchBioseq (bsp, cpp, 0);
    slp = MatchSa2Sl (&sap);
    if (slp != NULL) {
      pspp = PackSeqPntNew ();
      pspp->id = SeqIdDup (SeqIdFindBest (bsp->id, 0));
      while (slp != NULL) {
        pt = SeqLocStart (slp);
        PackSeqPntPut (pspp, pt);
        slpn = slp->next;
        SeqLocFree (slp);
        slp = slpn;
      }
      AddRestrictionSite (annot, pspp, cpp->name);
    }
    cpp = cpp->nextpattern;
  }

  if (annot->data == NULL) {
    SeqAnnotFree (annot);
    return;
  }
  if (bsp->annot == NULL) {
    bsp->annot = annot;
  } else {
    lastannot = bsp->annot;
    while (lastannot->next != NULL) {
      lastannot = lastannot->next;
    }
    lastannot->next = annot;
  }
}
#endif

static void SimpleRsiteProc (IteM i)

{
  BaseFormPtr   bfp;
  BioseqPtr     bsp;
  NewObjectPtr  nop;
  SeqEntryPtr   sep = NULL;
  Int2          which;
#ifdef TEST_APPLE_EVENT_MESSAGING
  AsnIoPtr      aip;
  Char          tmp [PATH_MAX];
#else
  ComPatPtr     cpph;
  Char          enzyme [64];
  Int2          j;
  Char          temp [32];
  ValNodePtr    enzymes;
#endif

  nop = (NewObjectPtr) GetObjectExtra (i);
  if (nop == NULL) return;
#ifdef WIN_MAC
  bfp = (BaseFormPtr) currentFormDataPtr;
#else
  bfp = nop->bfp;
#endif
  if (bfp == NULL) return;
  which = BioseqViewOrDocSumChoice (nop);
  if (which != 1) return;
  bsp =  GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID, bfp->input_itemtype);
  if (bsp != NULL) {
    sep = SeqMgrGetSeqEntryForData (bsp);
  } else {
    sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  }
  if (sep == NULL) return;

#ifdef TEST_APPLE_EVENT_MESSAGING
  TmpNam (tmp);
  aip = AsnIoOpen (tmp, "w");
  if (aip != NULL) {
    SeqEntryAsnWrite (sep, aip, NULL);
    AsnIoClose (aip);
    /* Nlm_SendOpenDocAppleEventEx (tmp, "REST", NULL, TRUE); */
    Nlm_SendOpenDocAppleEventEx (tmp, NULL, "RsiteFind", TRUE);
  }
#else
  enzymes = NULL;
  j = 1;
  sprintf (temp, "ENZ_%d", (int) j);
  while (GetAppParam ("SEQNCGIS", "ENZYMES", temp, NULL, enzyme, sizeof (enzyme) - 1)) {
    ValNodeCopyStr (&enzymes, 0, enzyme);
    j++;
    sprintf (temp, "ENZ_%d", (int) j);
  }
  if (enzymes == NULL) {
    ValNodeCopyStr (&enzymes, 0, "BamHI");
    ValNodeCopyStr (&enzymes, 0, "EcoRI");
    ValNodeCopyStr (&enzymes, 0, "HindIII");
  }
  cpph = ReadRENPattern ("ncbiren.dat", FALSE, enzymes);
  /* PalindromeCheck (cpph); */
  SeqEntryExplore (sep, (Pointer) cpph, RestrictionSearchProc);
  ComPatFree (cpph);
  ValNodeFreeData (enzymes);

  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  Update ();
#endif
}

static VQUEUE  vsquerylist = NULL;

static void LIBCALLBACK AnnounceCallback (CharPtr requestID, CharPtr seqID, Int2 estimatedSeconds)

{
  if (StringHasNoText (requestID)) {
    requestID = "?";
  }
  if (StringHasNoText (seqID)) {
    seqID = "?";
  }
  Message (MSG_POST, "Queued rID %s, seqID %s, estimated seconds = %d",
           requestID, seqID, (int) estimatedSeconds);
}

static Boolean LIBCALLBACK VecScreenCallback (
  CharPtr filename,
  VoidPtr userdata,
  CharPtr requestID,
  CharPtr seqID,
  Boolean success
)

{
  if (StringHasNoText (requestID)) {
    requestID = "?";
  }
  if (StringHasNoText (seqID)) {
    seqID = "?";
  }
  if (success) {
    if (! SequinHandleNetResults (filename)) {
      /* LaunchGeneralTextViewer (path, "QueueFastaQueryToURL failed"); */
    }
  } else {
    Message (MSG_POST, "Failure of rID '%s', seqID %s", requestID, seqID);
  }
  return TRUE;
}

static void SimpleUniVecScreenProc (IteM i)

{
  BaseFormPtr   bfp;
  BioseqPtr     bsp;
  NewObjectPtr  nop;
  SeqEntryPtr   sep = NULL;
  Int2          which;

  nop = (NewObjectPtr) GetObjectExtra (i);
  if (nop == NULL) return;
#ifdef WIN_MAC
  bfp = (BaseFormPtr) currentFormDataPtr;
#else
  bfp = nop->bfp;
#endif
  if (bfp == NULL) return;
  which = BioseqViewOrDocSumChoice (nop);
  if (which != 1) return;
  bsp =  GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID, bfp->input_itemtype);
  if (bsp != NULL) {
    sep = SeqMgrGetSeqEntryForData (bsp);
  } else {
    sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  }
  if (sep == NULL) return;

  VecScreenAsynchronousRequest ("UniVec", bsp, &vsquerylist, VecScreenCallback, AnnounceCallback, NULL);
}

static void SimpleUniVecCoreScreenProc (IteM i)

{
  BaseFormPtr   bfp;
  BioseqPtr     bsp;
  NewObjectPtr  nop;
  SeqEntryPtr   sep = NULL;
  Int2          which;

  nop = (NewObjectPtr) GetObjectExtra (i);
  if (nop == NULL) return;
#ifdef WIN_MAC
  bfp = (BaseFormPtr) currentFormDataPtr;
#else
  bfp = nop->bfp;
#endif
  if (bfp == NULL) return;
  which = BioseqViewOrDocSumChoice (nop);
  if (which != 1) return;
  bsp =  GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID, bfp->input_itemtype);
  if (bsp != NULL) {
    sep = SeqMgrGetSeqEntryForData (bsp);
  } else {
    sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  }
  if (sep == NULL) return;

  VecScreenAsynchronousRequest ("UniVec_Core", bsp, &vsquerylist, VecScreenCallback, AnnounceCallback, NULL);
}

static QBQUEUE  qbquerylist = NULL;

static void LIBCALLBACK QBAnnounceCallback (CharPtr requestID, CharPtr seqID, Int2 estimatedSeconds)

{
  if (StringHasNoText (requestID)) {
    requestID = "?";
  }
  if (StringHasNoText (seqID)) {
    seqID = "?";
  }
  Message (MSG_POST, "Queued rID %s, seqID %s, estimated seconds = %d",
           requestID, seqID, (int) estimatedSeconds);
}

static Boolean LIBCALLBACK QBCallback (
  CharPtr filename,
  VoidPtr userdata,
  CharPtr requestID,
  CharPtr seqID,
  Boolean success
)

{
  if (StringHasNoText (requestID)) {
    requestID = "?";
  }
  if (StringHasNoText (seqID)) {
    seqID = "?";
  }
  if (success) {
    if (! SequinHandleNetResults (filename)) {
      /* LaunchGeneralTextViewer (path, "QueueFastaQueryToURL failed"); */
    }
  } else {
    Message (MSG_POST, "Failure of rID '%s', seqID %s", requestID, seqID);
  }
  return TRUE;
}

static void SimpleQBlastProc (IteM i)

{
  BaseFormPtr   bfp;
  BioseqPtr     bsp;
  NewObjectPtr  nop;
  SeqEntryPtr   sep = NULL;
  Int2          which;

  nop = (NewObjectPtr) GetObjectExtra (i);
  if (nop == NULL) return;
#ifdef WIN_MAC
  bfp = (BaseFormPtr) currentFormDataPtr;
#else
  bfp = nop->bfp;
#endif
  if (bfp == NULL) return;
  which = BioseqViewOrDocSumChoice (nop);
  if (which != 1) return;
  bsp =  GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID, bfp->input_itemtype);
  if (bsp != NULL) {
    sep = SeqMgrGetSeqEntryForData (bsp);
  } else {
    sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  }
  if (sep == NULL) return;

  QBlastAsynchronousRequest ("nr", "blastn", bsp, &qbquerylist, QBCallback, QBAnnounceCallback, NULL);
}

/* Analysis menu can launch external programs or use Web services */

static QUEUE  urlquerylist = NULL;

static Int4 pendingqueries = 0;

static Boolean LIBCALLBACK SubmitToNCBIResultProc (CONN conn, VoidPtr userdata, EIO_Status status)

{
  AsnIoPtr     aop;
  FILE         *fp;
  Char         path [PATH_MAX];
  SeqEntryPtr  sep;

  TmpNam (path);
  fp = FileOpen (path, "wb");
  QUERY_CopyResultsToFile (conn, fp);
  FileClose (fp);
  aop = AsnIoOpen (path, "rb");
  sep = SeqEntryAsnRead (aop, NULL);
  AsnIoClose (aop);
  aop = AsnIoOpen (path, "w");
  SeqEntryAsnWrite (sep, aop, NULL);
  AsnIoFlush (aop);
  AsnIoClose (aop);
  LaunchGeneralTextViewer (path, "Echo binary transformation of Seq-entry");
  FileRemove (path);
  return TRUE;
}

extern void SubmitToNCBI (IteM i);
extern void SubmitToNCBI (IteM i)

{
  AsnIoPtr     aop;
  BaseFormPtr  bfp;
  CONN         conn;
  FILE         *fp;
  Char         path [PATH_MAX];
  Char         progname [64];
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  sprintf (progname, "Sequin/%s", SEQUIN_APPLICATION);

  TmpNam (path);

  aop = AsnIoOpen (path, "wb");
  SeqEntryAsnWrite (sep, aop, NULL);
  AsnIoFlush (aop);
  AsnIoClose (aop);

  conn = QUERY_OpenUrlQuery ("cruncher.nlm.nih.gov", 80,
                             "/cgi-bin/Sequin/testcgi.cgi", "request=echo",
                             progname, 30, eMIME_T_NcbiData, eMIME_AsnBinary,
                             eENCOD_Url,
                             fHCC_UrlDecodeInput | fHCC_UrlEncodeOutput);

  fp = FileOpen (path, "rb");
  QUERY_CopyFileToQuery (conn, fp);
  FileClose (fp);

  QUERY_SendQuery (conn);

  QUERY_AddToQueue (&urlquerylist, conn, SubmitToNCBIResultProc, NULL, TRUE);

  pendingqueries++;

  FileRemove (path);
}

extern void SequinCheckSocketsProc (void)

{
  Int4  remaining;

  remaining = QUERY_CheckQueue (&urlquerylist);
  if (remaining < pendingqueries) {
    Beep ();
    pendingqueries--;
  }
  remaining = VecScreenCheckQueue (&vsquerylist);
  remaining = QBlastCheckQueue (&qbquerylist);
}

static Boolean LIBCALLBACK DemoModeResultProc (CONN conn, VoidPtr userdata, EIO_Status status)

{
  FILE  *fp;
  Char  path [PATH_MAX];

  TmpNam (path);
  fp = FileOpen (path, "w");
  QUERY_CopyResultsToFile (conn, fp);
  FileClose (fp);
  LaunchGeneralTextViewer (path, "QueueFastaQueryToURL results");
  FileRemove (path);
  return TRUE;
}

static Boolean LIBCALLBACK SequinHandleURLResults (CONN conn, VoidPtr userdata, EIO_Status status)

{
  FILE  *fp;
  Char  path [PATH_MAX];

  TmpNam (path);
  fp = FileOpen (path, "w");
  QUERY_CopyResultsToFile (conn, fp);
  FileClose (fp);
  if (! SequinHandleNetResults (path)) {
    /* LaunchGeneralTextViewer (path, "QueueFastaQueryToURL failed"); */
  }
  FileRemove (path);
  return TRUE;
}

static void FinishURLProc (NewObjectPtr nop, CharPtr arguments, CharPtr path)

{
  CONN             conn;
  FILE             *fp;
  Char             progname [64];
  QueryResultProc  resultproc;
  EMIME_SubType    subtype;

  sprintf (progname, "Sequin/%s", SEQUIN_APPLICATION);

  if (nop->demomode) {
    resultproc = DemoModeResultProc;
  } else {
    resultproc = nop->resultproc;
  }
  if (nop->format == 1) {
    subtype = eMIME_Fasta;
  } else if (nop->format == 2) {
    subtype = eMIME_AsnText;
  } else {
    subtype = eMIME_Unknown;
  }

  conn = QUERY_OpenUrlQuery (nop->host_machine, nop->host_port,
                             nop->host_path, arguments,
                             progname, nop->timeoutsec,
                             eMIME_T_NcbiData, subtype, eENCOD_Url,
                             fHCC_UrlDecodeInput | fHCC_UrlEncodeOutput);
  if (conn == NULL) return;

  fp = FileOpen (path, "r");
  QUERY_CopyFileToQuery (conn, fp);
  FileClose (fp);

  QUERY_SendQuery (conn);

  QUERY_AddToQueue (&urlquerylist, conn, resultproc, NULL, TRUE);

  pendingqueries++;
}

static void DoAnalysisProc (NewObjectPtr nop, BaseFormPtr bfp, Int2 which, CharPtr arguments, ResultProc dotheanalysis)

{
  AsnIoPtr     aop;
  BioseqPtr    bsp;
  Char         path1 [PATH_MAX];
  SeqEntryPtr  sep;

  if (nop == NULL || bfp == NULL) return;
  switch (which) {
    case 1 :
      if (BioseqViewCanSaveFasta (bfp->form, nop->fastaNucOK, nop->fastaProtOK, nop->onlyBspTarget)) {
        TmpNam (path1);
        switch (nop->format) {
          case 1 :
            ExportBioseqViewFasta (bfp->form, path1, nop->fastaNucOK, nop->fastaProtOK, nop->onlyBspTarget);
            break;
          case 2 :
            sep = NULL;
            if (nop->onlyBspTarget) {
              bsp =  GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID, bfp->input_itemtype);
              sep = SeqMgrGetSeqEntryForData (bsp);
            } else {
              sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
            }
            if (sep != NULL) {
              aop = AsnIoOpen (path1, "w");
              SeqEntryAsnWrite (sep, aop, NULL);
              AsnIoFlush (aop);
              AsnIoClose (aop);
            }
            break;
          default :
            break;
        }
        if (dotheanalysis != NULL) {
          dotheanalysis (path1);
        } else {
          FinishURLProc (nop, arguments, path1);
        }
        FileRemove (path1);
      } else {
        ErrPostEx (SEV_ERROR, 0, 0, "BioseqView cannot save fasta format");
      }
      break;
    case 2 :
      if (DocSumCanSaveFasta (bfp->form, nop->fastaNucOK, nop->fastaProtOK)) {
        TmpNam (path1);
        ExportDocSumFasta (bfp->form, path1, nop->fastaNucOK, nop->fastaProtOK);
        if (dotheanalysis != NULL) {
          dotheanalysis (path1);
        } else {
          FinishURLProc (nop, arguments, path1);
        }
        FileRemove (path1);
      } else {
        ErrPostEx (SEV_ERROR, 0, 0, "DocSum cannot save fasta format");
      }
      break;
    default :
      break;
  }
}

/* encodes spaces as %20 in URLs */
static CharPtr StrSaveNoNullEncodeSpaces (CharPtr from)

{
  Char     ch;
  size_t   len = 0;
  CharPtr  p;
  CharPtr  q;
  CharPtr  to;

  if (StringHasNoText (from)) return NULL;
  p = from;
  ch = *p;
  while (ch != '\0') {
    if (ch == ' ') {
      len += 3;
    } else {
      len++;
    }
    p++;
    ch = *p;
  }
  to = MemNew (len + 1);
  if (to == NULL) return NULL;

  q = to;
  p = from;
  ch = *p;
  while (ch != '\0') {
    if (ch == ' ') {
      *q = '%';
      q++;
      *q = '2';
      q++;
      *q = '0';
      q++;
    } else {
      *q = ch;
      q++;
    }
    p++;
    ch = *p;
  }
  *q = '\0';
  return to;
}

typedef struct urlargform {
  FORM_MESSAGE_BLOCK

  NewObjectPtr nop;
  BaseFormPtr  bfp;
  ValNodePtr   controls;
  ValNodePtr   helps;
  Int2         which;
} UrlArgForm, PNTR UrlArgFormPtr;

static void AcceptArgumentFormProc (ButtoN b)

{
  CharPtr        args = NULL;
  CharPtr        arguments = NULL;
  ButtoN         btn;
  Char           ch;
  Int2           choice;
  Char           cpy [256];
  GrouP          grp;
  ValNodePtr     head = NULL;
  Int2           i;
  CharPtr        itms;
  CharPtr        last;
  size_t         len;
  LisT           lst;
  NewObjectPtr   nop;
  Boolean        notFirst = FALSE;
  PopuP          pop;
  ValNodePtr     ppt;
  CharPtr        ptr;
  CharPtr        str;
  Char           tmp [256];
  TexT           txt;
  UrlArgFormPtr  ufp;
  UrlParamPtr    upp;
  Int2           val;
  ValNodePtr     vnp;

  ufp = (UrlArgFormPtr) GetObjectExtra (b);
  if (ufp == NULL) return;
  Hide (ufp->form);
  Update ();
  nop = ufp->nop;
  if (nop != NULL) {
    if (! StringHasNoText (nop->prefix)) {
      ValNodeCopyStr (&head, 0, nop->prefix);
    }
    for (vnp = ufp->controls, ppt = nop->paramlist;
         vnp != NULL && ppt != NULL;
         vnp = vnp->next, ppt = ppt->next) {
      upp = (UrlParamPtr) ppt->data.ptrvalue;
      if (upp == NULL) continue;
      choice = vnp->choice;
      switch (upp->type) {
        case 1 :
          txt = (TexT) vnp->data.ptrvalue;
          str = SaveStringFromText (txt);
          if (str != NULL) {
            sprintf (tmp, "%s=%s", upp->param, str);
            ValNodeCopyStr (&head, ppt->choice, tmp);
            MemFree (str);
          }
          break;
        case 2 :
          btn = (ButtoN) vnp->data.ptrvalue;
          if (GetStatus (btn)) {
            sprintf (tmp, "%s=TRUE", upp->param);
          } else {
            sprintf (tmp, "%s=FALSE", upp->param);
          }
          ValNodeCopyStr (&head, ppt->choice, tmp);
          break;
        case 3 :
          pop = (PopuP) vnp->data.ptrvalue;
          val = GetValue (pop);
          if (val > 0) {
            i = 0;
            itms = upp->choices;
            StringNCpy_0 (tmp, itms, sizeof (tmp));
            last = tmp;
            ptr = last;
            ch = *ptr;
            while (ch != '\0') {
              if (ch == ',') {
                *ptr = '\0';
                i++;
                if (val == i) {
                  sprintf (cpy, "%s=%s", upp->param, last);
                  ValNodeCopyStr (&head, ppt->choice, cpy);
                }
                ptr++;
                last = ptr;
                ch = *ptr;
              } else {
                ptr++;
                ch = *ptr;
              }
            }
            if (! StringHasNoText (last)) {
              i++;
              if (val == i) {
                sprintf (cpy, "%s=%s", upp->param, last);
                ValNodeCopyStr (&head, ppt->choice, cpy);
              }
            }
          }
          break;
        case 4 :
          grp = (GrouP) vnp->data.ptrvalue;
          val = GetValue (grp);
          if (val > 0) {
            i = 0;
            itms = upp->choices;
            StringNCpy_0 (tmp, itms, sizeof (tmp));
            last = tmp;
            ptr = last;
            ch = *ptr;
            while (ch != '\0') {
              if (ch == ',') {
                *ptr = '\0';
                i++;
                if (val == i) {
                  sprintf (cpy, "%s=%s", upp->param, last);
                  ValNodeCopyStr (&head, ppt->choice, cpy);
                }
                ptr++;
                last = ptr;
                ch = *ptr;
              } else {
                ptr++;
                ch = *ptr;
              }
            }
            if (! StringHasNoText (last)) {
              i++;
              if (val == i) {
                sprintf (cpy, "%s=%s", upp->param, last);
                ValNodeCopyStr (&head, ppt->choice, cpy);
              }
            }
          }
          break;
        case 5 :
          lst = (LisT) vnp->data.ptrvalue;
          val = GetValue (lst);
          if (val > 0) {
            i = 0;
            itms = upp->choices;
            StringNCpy_0 (tmp, itms, sizeof (tmp));
            last = tmp;
            ptr = last;
            ch = *ptr;
            while (ch != '\0') {
              if (ch == ',') {
                *ptr = '\0';
                i++;
                if (val == i) {
                  sprintf (cpy, "%s=%s", upp->param, last);
                  ValNodeCopyStr (&head, ppt->choice, cpy);
                }
                ptr++;
                last = ptr;
                ch = *ptr;
              } else {
                ptr++;
                ch = *ptr;
              }
            }
            if (! StringHasNoText (last)) {
              i++;
              if (val == i) {
                sprintf (cpy, "%s=%s", upp->param, last);
                ValNodeCopyStr (&head, ppt->choice, cpy);
              }
            }
          }
          break;
        default :
          break;
      }
    }
    head = SortValNode (head, SortByVnpChoice);
    if (! StringHasNoText (nop->suffix)) {
      ValNodeCopyStr (&head, 0, nop->suffix);
    }
    for (len = 0, vnp = head; vnp != NULL; vnp = vnp->next) {
      len += StringLen ((CharPtr) vnp->data.ptrvalue) + 1;
    }
    if (len > 0) {
      arguments = MemNew (len + 5);
      if (arguments != NULL) {
        vnp = head;
        notFirst = FALSE;
        while (vnp != NULL) {
          if (notFirst) {
            StringCat (arguments, "&");
          }
          StringCat (arguments, (CharPtr) vnp->data.ptrvalue);
          notFirst = TRUE;
          vnp = vnp->next;
        }
      }
    }
    args = /* StrSaveNoNullEncodeSpaces */ StringSave (arguments);
    MemFree (arguments);
    DoAnalysisProc (nop, ufp->bfp, ufp->which, args, NULL);
    MemFree (args);
  }
  Remove (ufp->form);
}

static void ShowArgumentHelp (ButtoN b)

{
  NewObjectPtr   nop;
  ValNodePtr     ppt;
  CharPtr        str;
  UrlArgFormPtr  ufp;
  UrlParamPtr    upp;
  ValNodePtr     vnp;

  ufp = (UrlArgFormPtr) GetObjectExtra (b);
  if (ufp == NULL) return;
  nop = ufp->nop;
  if (nop == NULL) return;
  for (vnp = ufp->helps, ppt = nop->paramlist;
         vnp != NULL && ppt != NULL;
         vnp = vnp->next, ppt = ppt->next) {
    upp = (UrlParamPtr) ppt->data.ptrvalue;
    if (upp == NULL) continue;
    if ((Pointer) b == (Pointer) vnp->data.ptrvalue) {
      str = upp->help;
      Message (MSG_OK, "%s", str);
      return;
    }
  }
  Beep ();
}

static void ArgumentFormMessage (ForM f, Int2 mssg)

{
  UrlArgFormPtr  ufp;

  ufp = (UrlArgFormPtr) GetObjectExtra (f);
  if (ufp != NULL) {
    switch (mssg) {
      case VIB_MSG_CLOSE :
        Remove (f);
        break;
      default :
        if (ufp->appmessage != NULL) {
          ufp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static void CleanupArgumentForm (GraphiC g, VoidPtr data)

{
  UrlArgFormPtr  ufp;

  ufp = (UrlArgFormPtr) data;
  if (ufp != NULL) {
    ValNodeFree (ufp->controls);
    ValNodeFree (ufp->helps);
  }
  StdCleanupFormProc (g, data);
}

static ValNodePtr RearrangeParamList (ValNodePtr paramlist)

{
  ValNodePtr       curr;
  CharPtr          group;
  ValNodePtr       head = NULL;
  ValNodePtr       list;
  ValNodePtr       next;
  ValNodePtr PNTR  prev;
  ValNodePtr       ppt;
  UrlParamPtr      upp;

  ppt = paramlist;
  while (ppt != NULL) {
    list = ppt->next;
    ppt->next = NULL;
    ValNodeLink (&head, ppt);
    upp = (UrlParamPtr) ppt->data.ptrvalue;
    if (upp == NULL) {
      ppt = list;
      continue;
    }
    group = upp->group;
    curr = list;
    prev = &list;
    while (curr != NULL) {
      next = curr->next;
      upp = (UrlParamPtr) curr->data.ptrvalue;
      if (upp == NULL) {
        prev = &(curr->next);
        curr = next;
        continue;
      }
      if (StringICmp (upp->group, group) == 0) {
        *prev = next;
        curr->next = NULL;
        ValNodeLink (&head, curr);
      } else {
        prev = &(curr->next);
      }
      curr = next;
    }
    ppt = list;
  }
  return head;
}

static void BuildArgumentForm (NewObjectPtr nop, BaseFormPtr bfp, Int2 which)

{
  ButtoN             b;
  ButtoN             btn;
  GrouP              c;
  Char               ch;
  CharPtr            def;
  Int2               delta;
  TexT               first = NULL;
  GrouP              g;
  GrouP              grp;
  GrouP              h;
  ValNodePtr         hlp;
  Int2               i;
  CharPtr            itms;
  CharPtr            last;
  CharPtr            lastGroup = " ";
  LisT               lst;
  GrouP              m;
  Int2               max;
  Int2               min;
  ValNodePtr         moveMe = NULL;
  Nlm_Handle         obj1, obj2;
  PopuP              pop;
  PrompT             prmpt;
  ValNodePtr         ppt;
  CharPtr            ptr;
  RecT               r;
  StdEditorProcsPtr  sepp;
  CharPtr            str;
  Char               tmp [128];
  TexT               txt;
  UrlArgFormPtr      ufp;
  UrlParamPtr        upp;
  Int2               val;
  ValNodePtr         vnp;
  WindoW             w;

  if (nop == NULL || bfp == NULL) return;
  ufp = (UrlArgFormPtr) MemNew (sizeof (UrlArgForm));
  if (ufp == NULL) return;

  nop->paramlist = RearrangeParamList (nop->paramlist);

  w = FixedWindow (-50, -33, -10, -10, "Arguments", NULL);
  SetObjectExtra (w, ufp, CleanupArgumentForm);
  ufp->form = (ForM) w;
  ufp->formmessage = ArgumentFormMessage;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    ufp->appmessage = sepp->handleMessages;
  }

  ufp->bfp = bfp;
  ufp->nop = nop;
  ufp->which = which;

  m = HiddenGroup (w, 1, 0, NULL);

  g = NULL;
  for (ppt = nop->paramlist;
       ppt != NULL;
       ppt = ppt->next) {
    upp = (UrlParamPtr) ppt->data.ptrvalue;
    if (upp == NULL) continue;
    if (StringICmp (upp->group, lastGroup) != 0) {
      if (StringHasNoText (upp->group)) {
        if (StringHasNoText (lastGroup)) {
          g = HiddenGroup (m, 3, 0, NULL);
        } else {
          g = NormalGroup (m, 3, 0, "", programFont, NULL);
        }
      } else {
        g = NormalGroup (m, 3, 0, upp->group, programFont, NULL);
      }
      lastGroup = upp->group;
    }
    if (g == NULL) {
      g = HiddenGroup (m, 3, 0, NULL);
    }
    switch (upp->type) {
      case 1 :
        str = upp->prompt;
        StaticPrompt (g, str, 0, dialogTextHeight, programFont, 'l');
        def = upp->dfault;
        if (StringHasNoText (def)) {
          def = "";
        }
        txt = DialogText (g, def, 10, NULL);
        if (first == NULL) {
          first = txt;
        }
        ValNodeAddPointer (&(ufp->controls), 1, (Pointer) txt);
        ValNodeAddPointer (&moveMe, 0, (Pointer) txt);
        b = PushButton (g, "?", ShowArgumentHelp);
        SetObjectExtra (b, ufp, NULL);
        ValNodeAddPointer (&(ufp->helps), 0, (Pointer) b);
        break;
      case 2 :
        str = upp->prompt;
        btn = CheckBox (g, str, NULL);
        def = upp->dfault;
        if (StringICmp (def, "TRUE") == 0) {
          SetStatus (btn, TRUE);
        }
        prmpt = StaticPrompt (g, "", 0, 0, programFont, 'l');
        ValNodeAddPointer (&moveMe, 0, (Pointer) prmpt);
        ValNodeAddPointer (&(ufp->controls), 2, (Pointer) btn);
        b = PushButton (g, "?", ShowArgumentHelp);
        SetObjectExtra (b, ufp, NULL);
        ValNodeAddPointer (&(ufp->helps), 0, (Pointer) b);
        break;
      case 3 :
        str = upp->prompt;
        StaticPrompt (g, str, 0, dialogTextHeight, programFont, 'l');
        h = HiddenGroup (g, 1, 0, NULL);
        pop = PopupList (h, TRUE, NULL);
        def = upp->dfault;
        val = 0;
        i = 0;
        itms = upp->choices;
        StringNCpy_0 (tmp, itms, sizeof (tmp));
        last = tmp;
        ptr = last;
        ch = *ptr;
        while (ch != '\0') {
          if (ch == ',') {
            *ptr = '\0';
            PopupItem (pop, last);
            i++;
            if (StringICmp (def, last) == 0) {
              val = i;
            }
            ptr++;
            last = ptr;
            ch = *ptr;
          } else {
            ptr++;
            ch = *ptr;
          }
        }
        if (! StringHasNoText (last)) {
          PopupItem (pop, last);
          i++;
          if (StringICmp (def, last) == 0) {
            val = i;
          }
        }
        if (val > 0) {
          SetValue (pop, val);
        }
        ValNodeAddPointer (&(ufp->controls), 3, (Pointer) pop);
        ValNodeAddPointer (&moveMe, 0, (Pointer) pop);
        b = PushButton (g, "?", ShowArgumentHelp);
        SetObjectExtra (b, ufp, NULL);
        ValNodeAddPointer (&(ufp->helps), 0, (Pointer) b);
        break;
      case 4 :
        str = upp->prompt;
        StaticPrompt (g, str, 0, dialogTextHeight, programFont, 'l');
        h = HiddenGroup (g, 1, 0, NULL);
        grp = HiddenGroup (h, -3, 0, NULL);
        def = upp->dfault;
        val = 0;
        i = 0;
        itms = upp->choices;
        StringNCpy_0 (tmp, itms, sizeof (tmp));
        last = tmp;
        ptr = last;
        ch = *ptr;
        while (ch != '\0') {
          if (ch == ',') {
            *ptr = '\0';
            RadioButton (grp, last);
            i++;
            if (StringICmp (def, last) == 0) {
              val = i;
            }
            ptr++;
            last = ptr;
            ch = *ptr;
          } else {
            ptr++;
            ch = *ptr;
          }
        }
        if (! StringHasNoText (last)) {
          RadioButton (grp, last);
          i++;
          if (StringICmp (def, last) == 0) {
            val = i;
          }
        }
        if (val > 0) {
          SetValue (grp, val);
        }
        ValNodeAddPointer (&(ufp->controls), 4, (Pointer) grp);
        ValNodeAddPointer (&moveMe, 0, (Pointer) grp);
        b = PushButton (g, "?", ShowArgumentHelp);
        SetObjectExtra (b, ufp, NULL);
        ValNodeAddPointer (&(ufp->helps), 0, (Pointer) b);
        break;
      case 5 :
        str = upp->prompt;
        StaticPrompt (g, str, 0, dialogTextHeight, programFont, 'l');
        h = HiddenGroup (g, 1, 0, NULL);
        lst = SingleList (h, 10, 3, NULL);
        def = upp->dfault;
        val = 0;
        i = 0;
        itms = upp->choices;
        StringNCpy_0 (tmp, itms, sizeof (tmp));
        last = tmp;
        ptr = last;
        ch = *ptr;
        while (ch != '\0') {
          if (ch == ',') {
            *ptr = '\0';
            ListItem (lst, last);
            i++;
            if (StringICmp (def, last) == 0) {
              val = i;
            }
            ptr++;
            last = ptr;
            ch = *ptr;
          } else {
            ptr++;
            ch = *ptr;
          }
        }
        if (! StringHasNoText (last)) {
          ListItem (lst, last);
          i++;
          if (StringICmp (def, last) == 0) {
            val = i;
          }
        }
        if (val > 0) {
          SetValue (lst, val);
        }
        ValNodeAddPointer (&(ufp->controls), 5, (Pointer) lst);
        ValNodeAddPointer (&moveMe, 0, (Pointer) lst);
        b = PushButton (g, "?", ShowArgumentHelp);
        SetObjectExtra (b, ufp, NULL);
        ValNodeAddPointer (&(ufp->helps), 0, (Pointer) b);
        break;
      default :
        break;
    }
  }

  min = 0;
  max = 0;
  for (vnp = moveMe; vnp != NULL; vnp = vnp->next) {
    obj1 = (Nlm_Handle) vnp->data.ptrvalue;
    GetPosition (obj1, &r);
    min = MAX (min, r.left);
  }
  for (vnp = moveMe; vnp != NULL; vnp = vnp->next) {
    obj1 = (Nlm_Handle) vnp->data.ptrvalue;
    GetPosition (obj1, &r);
    delta = min - r.left;
    OffsetRect (&r, delta, 0);
    SetPosition (obj1, &r);
    AdjustPrnt (obj1, &r, FALSE);
    max = MAX (max, r.right);
  }
  max += 3;
  for (vnp = moveMe, hlp = ufp->helps;
       vnp != NULL && hlp != NULL;
       vnp = vnp->next, hlp = hlp->next) {
    obj2 = (Nlm_Handle) hlp->data.ptrvalue;
    GetPosition (obj2, &r);
    delta = max - r.left;
    OffsetRect (&r, delta, 0);
    SetPosition (obj2, &r);
    AdjustPrnt (obj2, &r, TRUE);
  }

  c = HiddenGroup (w, 2, 0, NULL);
  SetGroupSpacing (c, 10, 3);
  b = PushButton (c, "Accept", AcceptArgumentFormProc);
  SetObjectExtra (b, ufp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) m, (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Select (w);
  if (first != NULL) {
    Select (first);
  }
}

static void DoURLProc (IteM i)

{
  CharPtr       args = NULL;
  BaseFormPtr   bfp;
  size_t        len;
  NewObjectPtr  nop;
  Int2          which;

  nop = (NewObjectPtr) GetObjectExtra (i);
  if (nop == NULL) return;
#ifdef WIN_MAC
  bfp = (BaseFormPtr) currentFormDataPtr;
#else
  bfp = nop->bfp;
#endif
  if (bfp == NULL) return;
  which = BioseqViewOrDocSumChoice (nop);
  if (nop->paramlist == NULL) {
    len = StringLen (nop->prefix) + StringLen (nop->suffix);
    if (len > 0) {
      args = MemNew (sizeof (Char) * (len + 2));
      StringCpy (args, nop->prefix);
      if (! StringHasNoText (nop->suffix)) {
        StringCat (args, "&");
        StringCat (args, nop->suffix);
      }
    }
    DoAnalysisProc (nop, bfp, which, args, NULL);
  } else {
    BuildArgumentForm (nop, bfp, which);
  }
}

extern void EnableAnalysisItems (BaseFormPtr bfp, Boolean isDocSum)

{
  Boolean       hasFastaNuc;
  Boolean       hasFastaProt;
  NewObjectPtr  nop;

  if (bfp == NULL) return;
#ifdef WIN_MAC
  nop = (NewObjectPtr) macUserDataPtr;
#else
  nop = (NewObjectPtr) bfp->userDataPtr;
#endif
  if (isDocSum) {
  } else {
  }
  while (nop != NULL) {
    if (nop->kind == 1) {
      /* annotate menu item, ignore it */
    } else if (isDocSum) {
      if (nop->dsmOK) {
        hasFastaNuc = DocSumCanSaveFasta (bfp->form, TRUE, FALSE);
        hasFastaProt = DocSumCanSaveFasta (bfp->form, FALSE, TRUE);
        if (nop->fastaNucOK && hasFastaNuc) {
          SafeEnable (nop->item);
        } else if (nop->fastaProtOK && hasFastaProt) {
          SafeEnable (nop->item);
        } else {
          SafeDisable (nop->item);
        }
      } else {
        SafeDisable (nop->item);
      }
    } else {
      if (nop->bspOK) {
        hasFastaNuc = BioseqViewCanSaveFasta (bfp->form, TRUE, FALSE, nop->onlyBspTarget);
        hasFastaProt = BioseqViewCanSaveFasta (bfp->form, FALSE, TRUE, nop->onlyBspTarget);
        if (nop->fastaNucOK && hasFastaNuc) {
          SafeEnable (nop->item);
        } else if (nop->fastaProtOK && hasFastaProt) {
          SafeEnable (nop->item);
        } else {
          SafeDisable (nop->item);
        }
      } else {
        SafeDisable (nop->item);
      }
    }
    nop = nop->next;
  }
}

static VoidPtr LinkNewObjectLists (NewObjectPtr list1, NewObjectPtr list2)

{
  NewObjectPtr  nop;

  if (list1 == NULL) return list2;
  nop = list1;
  while (nop->next != NULL) {
    nop = nop->next;
  }
  nop->next = list2;
  return list1;
}

static void CleanupAnalysisExtraProc (GraphiC g, VoidPtr data)

{
  NewObjectPtr  nop;
  ValNodePtr    ppt;
  UrlParamPtr   upp;

  nop = (NewObjectPtr) data;
  if (nop != NULL) {
    MemFree (nop->host_machine);
    MemFree (nop->host_path);
    for (ppt = nop->paramlist; ppt != NULL; ppt = ppt->next) {
      upp = (UrlParamPtr) ppt->data.ptrvalue;
      if (upp == NULL) continue;
      MemFree (upp->param);
      MemFree (upp->prompt);
      MemFree (upp->dfault);
      MemFree (upp->choices);
      MemFree (upp->group);
      MemFree (upp->help);
    }
    ValNodeFreeData (nop->paramlist);
    MemFree (nop->prefix);
    MemFree (nop->suffix);
  }
  MemFree (data);
}

typedef struct sbstruc {
  CharPtr    name;
  MenU       menu;
} Sbstruc, PNTR SbstrucPtr;

static ValNodePtr  analysissubmenulist = NULL;

static void AddAnalysisItem (MenU m, BaseFormPtr bfp,
                             Boolean bspviewOK, Boolean docsumOK,
                             Boolean nucOK, Boolean protOK, Boolean onlyBspTarget,
                             CharPtr host_machine, Uint2 host_port,
                             CharPtr host_path, CharPtr program,
                             Uint2 timeoutsec, Int2 format, Boolean demomode,
                             QueryResultProc resultproc, ValNodePtr paramlist,
                             CharPtr prefix, CharPtr suffix,
                             CharPtr title, CharPtr submenu,
                             ItmActnProc actn, NewObjectPtr PNTR head)

{
  IteM          i;
  NewObjectPtr  last;
  size_t        len;
  NewObjectPtr  nop;
  SbstrucPtr    sbp;
  CharPtr       tmp;
  ValNodePtr    vnp;
  MenU          x;

  if (m == NULL || actn == NULL) return;
  x = NULL;
  if (! StringHasNoText (submenu)) {
    vnp = analysissubmenulist;
    while (vnp != NULL && x == NULL) {
      sbp = (SbstrucPtr) vnp->data.ptrvalue;
      if (sbp != NULL && StringICmp (sbp->name, submenu) == 0) {
        x = sbp->menu;
      }
      vnp = vnp->next;
    }
    if (x == NULL) {
      sbp = (SbstrucPtr) MemNew (sizeof (Sbstruc));
      if (sbp != NULL) {
        sbp->name = StringSave (submenu);
        sbp->menu = SubMenu (m, sbp->name);
        x = sbp->menu;
        ValNodeAddPointer (&analysissubmenulist, 0, (VoidPtr) sbp);
      }
    }
  }
  if (x == NULL) {
    x = m;
  }
  i = CommandItem (x, title, actn);
  nop = (NewObjectPtr) MemNew (sizeof (NewObjectData));
  if (nop != NULL) {
    nop->kind = 2; /* analysis menu item */
    nop->bfp = bfp;
    nop->item = i;
    nop->bspOK = bspviewOK;
    nop->dsmOK = docsumOK;
    nop->fastaNucOK = nucOK;
    nop->fastaProtOK = protOK;
    nop->onlyBspTarget = onlyBspTarget;
    nop->host_machine = /* StrSaveNoNullEncodeSpaces */ StringSave (host_machine);
    nop->host_port = host_port;
    len = StringLen (host_path);
    tmp = MemNew (len + StringLen (program) + 5);
    if (tmp != NULL) {
      StringCpy (tmp, host_path);
      if (len > 1 && tmp [len - 1] != '/') {
        StringCat (tmp, "/");
      }
      StringCat (tmp, program);
    }
    nop->host_path = /* StrSaveNoNullEncodeSpaces */ StringSave (tmp);
    MemFree (tmp);
    nop->query = NULL;
    /*
    nop->host_path = StrSaveNoNullEncodeSpaces (host_path);
    nop->query = StrSaveNoNullEncodeSpaces (program);
    */
    nop->timeoutsec = timeoutsec;
    nop->format = format;
    nop->demomode = demomode;
    nop->resultproc = resultproc;
    nop->paramlist = paramlist;
    nop->prefix = StringSaveNoNull (prefix);
    nop->suffix = StringSaveNoNull (suffix);
  }
  SetObjectExtra (i, (Pointer) nop, CleanupAnalysisExtraProc);
  if (head == NULL) return;
  last = *head;
  if (last != NULL) {
    while (last->next != NULL) {
     last = last->next;
    }
    last->next = nop;
  } else {
    *head = nop;
  }
}

/* Sample seqncgis.cnf/seqncgis.ini/.seqncgisrc/sequincgi.cfg config file.
   PATH can contain query (separated by ? symbol), or separate QUERY item can
   be used, or multiple QUERY and TITLE items can also be used.

[SERVICES]
PATH=mydisk:Common Files:services:

[ORDER]
ORDER_1=tRNAscan
ORDER_2=Seg

[tRNAscan]
PROGRAM=testcgi.cgi?request=trnascan
HOST=www.myserver.myschool.edu
PORT=80
PATH=/MyServices/cgi-bin/testcgi.cgi
SUBMENU=Search
FORMATIN=FASTA
FLAGS=SEQ,NUC,TRG
TIMEOUT=30

[Seg]
PROGRAM=segify
HOST=www.myserver.myschool.edu
PORT=80
PATH=/MyServices/cgi-bin/testcgi.cgi
FORMATIN=fasta
FLAGS=SEQ,DOC,PRT,TRG
SUBMENU=Secondary structure prediction
PROMPT_1=Window Size
PARAM_1=window
DESCRIPTION_1=window size for determining low-complexity segments
TYPE_1=text
DEFAULT_1=12
REQUIRED_1=FALSE
IMPORTANCE_1=
GROUP_1=Algorithm
HELP_1=window size for determining low-complexity segments
PROMPT_2=Trigger Complexity
PARAM_2=trigger
DESCRIPTION_2=trigger complexity for determining low-complexity segments
TYPE_2=text
DEFAULT_2=2.2
REQUIRED_2=FALSE
IMPORTANCE_2=
GROUP_2=Algorithm
HELP_2=trigger complexity for determining low-complexity segments
...

[ENZYMES]
ENZ_1=BamHI
ENZ_2=EcoRI
ENZ_3=HindIII

*/

static Int2 GetServiceParam (ValNodePtr head, CharPtr type, CharPtr buf, Int2 buflen)

{
  size_t      len;
  Boolean     seenBracket = FALSE;
  CharPtr     str;
  ValNodePtr  vnp;

  if (buf == NULL || buflen <= 0) return 0;
  *buf = '\0';
  len = StringLen (type);
  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (str != NULL) {
      if (str [0] == '[') {
        if (seenBracket) return 0;
        seenBracket = TRUE;
      } else if (StringNICmp (str, type, len) == 0) {
        str += len;
        StringNCpy_0 (buf, str, buflen);
        return (Int2) StringLen (buf);
      }
    }
  }
  return 0;
}

static ValNodePtr GetConfigParamAndPromptLists (CharPtr sect)

{
  Int2         i;
  ValNodePtr   paramlist = NULL;
  Uint1        paramtype;
  Char         title [512];
  Char         tmp [32];
  UrlParamPtr  upp;

  if (sect == NULL) return NULL;
  i = 1;
  sprintf (tmp, "PARAM_%d", (int) i);
  while (GetAppParam ("SEQNCGIS", sect, tmp, NULL, title, sizeof (title) - 1)) {
    upp = (UrlParamPtr) MemNew (sizeof (UrlParamData));
    if (upp == NULL) continue;
    upp->param = StringSave (title);
    sprintf (tmp, "TYPE_%d", (int) i);
    paramtype = 1;
    if (GetAppParam ("SEQNCGIS", sect, tmp, NULL, title, sizeof (title) - 1)) {
      if (StringICmp (title, "text") == 0) {
        paramtype = 1;
      } else if (StringICmp (title, "checkbox") == 0) {
        paramtype = 2;
      } else if (StringICmp (title, "popup") == 0) {
        paramtype = 3;
      } else if (StringICmp (title, "radio") == 0) {
        paramtype = 4;
      } else if (StringICmp (title, "list") == 0) {
        paramtype = 5;
      }
    }
    upp->type = paramtype;
    sprintf (tmp, "PROMPT_%d", (int) i);
    if (GetAppParam ("SEQNCGIS", sect, tmp, NULL, title, sizeof (title) - 1)) {
      upp->prompt = StringSave (title);
    } else {
      upp->prompt = StringSave (upp->param);
    }
    sprintf (tmp, "DEFAULT_%d", (int) i);
    if (GetAppParam ("SEQNCGIS", sect, tmp, NULL, title, sizeof (title) - 1)) {
      upp->dfault = StringSave (title);
    } else {
      upp->dfault = StringSave (" ");
    }
    sprintf (tmp, "CHOICES_%d", (int) i);
    if (GetAppParam ("SEQNCGIS", sect, tmp, NULL, title, sizeof (title) - 1)) {
      upp->choices = StringSave (title);
    } else {
      upp->choices = StringSave (" ");
    }
    sprintf (tmp, "GROUP_%d", (int) i);
    if (GetAppParam ("SEQNCGIS", sect, tmp, NULL, title, sizeof (title) - 1)) {
      upp->group = StringSave (title);
    } else {
      upp->group = StringSave (" ");
    }
    sprintf (tmp, "HELP_%d", (int) i);
    if (GetAppParam ("SEQNCGIS", sect, tmp, NULL, title, sizeof (title) - 1)) {
      upp->help = StringSave (title);
    } else {
      upp->help = StringSave (" ");
    }
    ValNodeAddPointer (&paramlist, i, (Pointer) upp);
    i++;
    sprintf (tmp, "PARAM_%d", (int) i);
  }
  return paramlist;
}

static ValNodePtr GetServiceParamAndPromptLists (ValNodePtr list)

{
  Int2         i;
  ValNodePtr   paramlist = NULL;
  Uint1        paramtype;
  Char         title [512];
  Char         tmp [32];
  UrlParamPtr  upp;

  if (list == NULL) return NULL;
  i = 1;
  sprintf (tmp, "PARAM_%d=", (int) i);
  while (GetServiceParam (list, tmp, title, sizeof (title) - 1)) {
    upp = (UrlParamPtr) MemNew (sizeof (UrlParamData));
    if (upp == NULL) continue;
    upp->param = StringSave (title);
    sprintf (tmp, "TYPE_%d", (int) i);
    paramtype = 1;
    if (GetServiceParam (list, tmp, title, sizeof (title) - 1)) {
      if (StringICmp (title, "text") == 0) {
        paramtype = 1;
      } else if (StringICmp (title, "checkbox") == 0) {
        paramtype = 2;
      } else if (StringICmp (title, "popup") == 0) {
        paramtype = 3;
      } else if (StringICmp (title, "radio") == 0) {
        paramtype = 4;
      } else if (StringICmp (title, "list") == 0) {
        paramtype = 5;
      }
    }
    upp->type = paramtype;
    sprintf (tmp, "PROMPT_%d=", (int) i);
    if (GetServiceParam (list, tmp, title, sizeof (title) - 1)) {
      upp->prompt = StringSave (title);
    } else {
      upp->prompt = StringSave (upp->param);
    }
    sprintf (tmp, "DEFAULT_%d=", (int) i);
    if (GetServiceParam (list, tmp, title, sizeof (title) - 1)) {
      upp->dfault = StringSave (title);
    } else {
      upp->dfault = StringSave (" ");
    }
    sprintf (tmp, "CHOICES_%d=", (int) i);
    if (GetServiceParam (list, tmp, title, sizeof (title) - 1)) {
      upp->choices = StringSave (title);
    } else {
      upp->choices = StringSave (" ");
    }
    sprintf (tmp, "GROUP_%d=", (int) i);
    if (GetServiceParam (list, tmp, title, sizeof (title) - 1)) {
      upp->group = StringSave (title);
    } else {
      upp->group = StringSave (" ");
    }
    sprintf (tmp, "HELP_%d=", (int) i);
    if (GetServiceParam (list, tmp, title, sizeof (title) - 1)) {
      upp->help = StringSave (title);
    } else {
      upp->help = StringSave (" ");
    }
    ValNodeAddPointer (&paramlist, i, (Pointer) upp);
    i++;
    sprintf (tmp, "PARAM_%d=", (int) i);
  }
  return paramlist;
}

static void ReadAnalysisConfigFile (CharPtr sect, MenU m, BaseFormPtr bfp,
                                    Boolean bspviewOK, Boolean docsumOK,
                                    NewObjectPtr PNTR head)

{
  Boolean     demomode = FALSE;
  Int2        format = 1;
  Char        host [128];
  Boolean     nucOK = FALSE;
  Boolean     onlyBspTarget = FALSE;
  ValNodePtr  paramlist = NULL;
  Char        program [128];
  Char        path [256];
  Uint2       port = 80;
  Char        prefix [128];
  Boolean     protOK = FALSE;
  Char        submenu [128];
  Char        suffix [128];
  Uint2       timeoutsec = 30;
  Char        title [128];
  Char        tmp [32];
  unsigned    int  val;

  if (! GetAppParam ("SEQNCGIS", sect, "TITLE", NULL, title, sizeof (title) - 1)) {
    StringNCpy_0 (title, sect, sizeof (title));
  }
  if (GetAppParam ("SEQNCGIS", sect, "HOST", NULL, host, sizeof (host) - 1)) {
    if (GetAppParam ("SEQNCGIS", sect, "FLAGS", NULL, tmp, sizeof (tmp) - 1)) {
      if (StringStr (tmp, "SEQ") == NULL) {
        bspviewOK = FALSE;
      }
      if (StringStr (tmp, "DOC") == NULL) {
        docsumOK = FALSE;
      }
      if (StringStr (tmp, "NUC") != NULL) {
        nucOK = TRUE;
      }
      if (StringStr (tmp, "PRT") != NULL) {
        protOK = TRUE;
      }
      if (StringStr (tmp, "TRG") != NULL) {
        onlyBspTarget = TRUE;
      }
    }

    if ((! bspviewOK) && (! docsumOK)) return;

    if (GetAppParam ("SEQNCGIS", sect, "PORT", NULL, tmp, sizeof (tmp) - 1) && 
        sscanf (tmp, "%u", &val) == 1) {
      port = (Uint2) val;
    } else {
      port = 80;
    }
    if (GetAppParam ("SEQNCGIS", sect, "FORMATIN", NULL, tmp, sizeof (tmp) - 1)) {
      if (StringICmp (tmp, "FASTA") == 0) {
        format = 1;
      } else if (StringICmp (tmp, "ASN.1") == 0) {
        format = 2;
      }
    }
    if (GetAppParam ("SEQNCGIS", sect, "TIMEOUT", NULL, tmp, sizeof (tmp) - 1) && 
        sscanf (tmp, "%u", &val) == 1) {
      timeoutsec = (Uint2) val;
    } else {
      timeoutsec = 30;
    }
    submenu [0] = '\0';
    GetAppParam ("SEQNCGIS", sect, "SUBMENU", NULL, submenu, sizeof (submenu) - 1);
    if (GetAppParam ("SEQNCGIS", sect, "DEMO", NULL, tmp, sizeof (tmp) - 1)) {
      if (StringICmp (tmp, "TRUE") == 0) {
        demomode = TRUE;
      }
    }

    if (GetAppParam ("SEQNCGIS", sect, "PATH", NULL, path, sizeof (path) - 1)) {
      if (GetAppParam ("SEQNCGIS", sect, "PROGRAM", NULL, program, sizeof (program) - 1)) {
        paramlist = GetConfigParamAndPromptLists (sect);
        prefix [0] = '\0';
        GetAppParam ("SEQNCGIS", sect, "PREFIX", NULL, prefix, sizeof (prefix) - 1);
        suffix [0] = '\0';
        GetAppParam ("SEQNCGIS", sect, "SUFFIX", NULL, suffix, sizeof (suffix) - 1);
        AddAnalysisItem (m, bfp, bspviewOK, docsumOK,
                         nucOK, protOK, onlyBspTarget,
                         host, port, path, program, timeoutsec, format, demomode,
                         SequinHandleURLResults, paramlist, prefix, suffix,
                         title, submenu, DoURLProc, head);
      }
    }
  }
}

static void ReadServiceConfigFile (CharPtr pathbase, ValNodePtr config,
                                   MenU m, BaseFormPtr bfp,
                                   Boolean bspviewOK, Boolean docsumOK,
                                   NewObjectPtr PNTR head)

{
  Char          ch;
  ValNodePtr    choicelist = NULL;
  Boolean       demomode = FALSE;
  ValNodePtr    dfaultlist = NULL;
  Int2          format = 1;
  FILE          *fp;
  Boolean       goOn = TRUE;
  ValNodePtr    helplist = NULL;
  Char          host [128];
  Boolean       keepGoing;
  ValNodePtr    list = NULL;
  Boolean       nucOK = FALSE;
  Boolean       onlyBspTarget = FALSE;
  ValNodePtr    paramlist = NULL;
  Char          program [128];
  ValNodePtr    promptlist = NULL;
  Char          path [PATH_MAX];
  Uint2         port = 80;
  Char          prefix [128];
  Boolean       protOK = FALSE;
  CharPtr       ptr;
  Boolean       seenBracket;
  Char          str [256];
  Char          submenu [128];
  Char          suffix [128];
  Uint2         timeoutsec = 30;
  Char          title [128];
  Char          tmp [32];
  unsigned int  val;
  ValNodePtr    vnp;

  if (path == NULL || config == NULL || config->data.ptrvalue == NULL) return;
  StringNCpy_0 (path, pathbase, sizeof (path));
  FileBuildPath (path, NULL, (CharPtr) config->data.ptrvalue);
  fp = FileOpen (path, "r");
  if (fp == NULL) return;
  while (fgets (str, sizeof (str), fp) != NULL) {
    ptr = str;
    ch = *ptr;
    while (ch != '\0' && ch != '\n' && ch != '\r') {
      ptr++;
      ch = *ptr;
    }
    *ptr = '\0';
    ValNodeCopyStr (&list, 1, str);
  }
  FileClose (fp);
  while (goOn) {
    goOn = FALSE;
    title [0] = '\0';
    if (GetServiceParam (list, "TITLE=", tmp, sizeof (tmp) - 1)) {
      StringNCpy_0 (title, tmp, sizeof (title));
    }
    if (StringHasNoText (title)) {
      if (GetServiceParam (list, "[", title, sizeof (title) - 1)) {
        ptr = StringChr (title, ']');
        if (ptr != NULL) {
          *ptr = '\0';
        }
      }
    }
    if (title [0] != '\0' && GetServiceParam (list, "HOST=", host, sizeof (host) - 1)) {
      if (GetServiceParam (list, "FLAGS=", tmp, sizeof (tmp) - 1)) {
        if (StringStr (tmp, "SEQ") == NULL) {
          bspviewOK= FALSE;
        }
        if (StringStr (tmp, "DOC") == NULL) {
          docsumOK= FALSE;
        }
        if (StringStr (tmp, "NUC") != NULL) {
          nucOK= TRUE;
        }
        if (StringStr (tmp, "PRT") != NULL) {
          protOK= TRUE;
        }
        if (StringStr (tmp, "TRG") != NULL) {
          onlyBspTarget= TRUE;
        }
      }

      if (bspviewOK || docsumOK) {

        if (GetServiceParam (list, "PORT=", tmp, sizeof (tmp) - 1) && 
            sscanf (tmp, "%u", &val) == 1) {
          port = (Uint2) val;
        } else {
          port = 80;
        }
        if (GetServiceParam (list, "FORMATIN=", tmp, sizeof (tmp) - 1)) {
          if (StringICmp (tmp, "FASTA") == 0) {
            format = 1;
          } else if (StringICmp (tmp, "ASN.1") == 0) {
            format = 2;
          }
        }
       if (GetServiceParam (list, "TIMEOUT=", tmp, sizeof (tmp) - 1) && 
            sscanf (tmp, "%u", &val) == 1) {
          timeoutsec = (Uint2) val;
        } else {
          timeoutsec = 30;
        }
        submenu [0] = '\0';
        GetServiceParam (list, "SUBMENU=", submenu, sizeof (submenu) - 1);
        if (GetServiceParam (list, "DEMO=", tmp, sizeof (tmp) - 1)) {
          if (StringICmp (tmp, "TRUE") == 0) {
            demomode = TRUE;
          }
        }

        if (GetServiceParam (list, "PATH=", path, sizeof (path) - 1)) {
          if (GetServiceParam (list, "PROGRAM=", program, sizeof (program) - 1)) {
            paramlist = GetServiceParamAndPromptLists (list);
            prefix [0] = '\0';
            GetServiceParam (list, "PREFIX=", prefix, sizeof (prefix) - 1);
            suffix [0] = '\0';
            GetServiceParam (list, "SUFFIX=", suffix, sizeof (suffix) - 1);
            AddAnalysisItem (m, bfp, bspviewOK, docsumOK,
                             nucOK, protOK, onlyBspTarget,
                             host, port, path, program, timeoutsec, format, demomode,
                             SequinHandleURLResults, paramlist, prefix, suffix,
                             title, submenu, DoURLProc, head);
          }
        }

      }
    }

    seenBracket = FALSE;
    keepGoing = TRUE;
    for (vnp = list; vnp != NULL && keepGoing; vnp = vnp->next) {
      ptr = (CharPtr) vnp->data.ptrvalue;
      if (ptr != NULL) {
        if (ptr [0] == '[') {
          if (seenBracket) {
            keepGoing = FALSE;
          } else {
            seenBracket = TRUE;
          }
        }
        if (keepGoing) {
          vnp->data.ptrvalue = MemFree (vnp->data.ptrvalue);
        }
      }
    }

  }

  ValNodeFreeData (list);
}

extern MenU CreateAnalysisMenu (WindoW w, BaseFormPtr bfp, Boolean bspviewOK, Boolean docsumOK)

{
  NewObjectPtr  first;
  ValNodePtr    head1 = NULL, head2 = NULL;
  Int2          i;
  size_t        len;
  MenU          m;
  Char          path1 [PATH_MAX];
  Char          path2 [PATH_MAX];
  CharPtr       ptr;
  SbstrucPtr    sbp;
  Char          sect [256];
  Char          temp [32];
  ValNodePtr    vnp;

  ProgramPath (path1, sizeof (path1));
  ptr = StringRChr (path1, DIRDELIMCHR);
  if (ptr != NULL) {
    ptr++;
    *ptr = '\0';
  }
  FileBuildPath (path1, "services", NULL);
  head1 = DirCatalog (path1);

  if (GetAppParam ("SEQNCGIS", "SERVICES", "PATH", NULL, path2, sizeof (path2) - 1)) {
    len = StringLen (path2);
    if (path2 [len - 1] != DIRDELIMCHR) {
      StringCat (path2, DIRDELIMSTR);
    }
    if (StringCmp (path1, path2) != 0) {
      head2 = DirCatalog (path2);
    }
  }

  if ((! extraServices) && (! indexerVersion) && (! genomeCenter) &&
      head1 == NULL && head2 == NULL) {
    if (! GetAppParam ("SEQNCGIS", "ORDER", NULL, NULL, sect, sizeof (sect) - 1)) {
      return NULL;
    }
  }
  m = PulldownMenu (w, "Analysis");
  if (m == NULL) return NULL;
  analysissubmenulist = NULL;
  first = NULL;
  if (bspviewOK) {
    AddAnalysisItem (m, bfp, bspviewOK, FALSE, TRUE, FALSE, TRUE,
                     NULL, 0, NULL, NULL, 0, 0, FALSE, NULL, NULL, NULL, NULL,
                     "Restriction Search", "Search",
                     SimpleRsiteProc, &first);
    if (indexerVersion) {
      AddAnalysisItem (m, bfp, bspviewOK, FALSE, TRUE, FALSE, TRUE,
                       NULL, 0, NULL, NULL, 0, 0, FALSE, NULL, NULL, NULL, NULL,
                       "Vector Screen - UniVec", "Search",
                       SimpleUniVecScreenProc, &first);
      AddAnalysisItem (m, bfp, bspviewOK, FALSE, TRUE, FALSE, TRUE,
                       NULL, 0, NULL, NULL, 0, 0, FALSE, NULL, NULL, NULL, NULL,
                       "Vector Screen - UniVec Core", "Search",
                       SimpleUniVecCoreScreenProc, &first);
      AddAnalysisItem (m, bfp, bspviewOK, FALSE, TRUE, FALSE, TRUE,
                       NULL, 0, NULL, NULL, 0, 0, FALSE, NULL, NULL, NULL, NULL,
                       "QBlast Test", "Search",
                       SimpleQBlastProc, &first);
    }
  }
  if (bspviewOK || docsumOK) {
    if (useEntrez) {
      i = 1;
      sprintf (temp, "ORDER_%d", (int) i);
      while (GetAppParam ("SEQNCGIS", "ORDER", temp, NULL, sect, sizeof (sect) - 1)) {
        ReadAnalysisConfigFile (sect, m, bfp, bspviewOK, docsumOK, &first);
        i++;
        sprintf (temp, "ORDER_%d", (int) i);
      }
      for (vnp = head1; vnp != NULL; vnp = vnp->next) {
        if (vnp->choice == 0) {
          ReadServiceConfigFile (path1, vnp, m, bfp, bspviewOK, docsumOK, &first);
        }
      }
      for (vnp = head2; vnp != NULL; vnp = vnp->next) {
        if (vnp->choice == 0) {
          ReadServiceConfigFile (path2, vnp, m, bfp, bspviewOK, docsumOK, &first);
        }
      }
    }
  }
  if (bspviewOK) {
  }
  if (docsumOK) {
  }
#ifdef WIN_MAC
  macUserDataPtr = LinkNewObjectLists (macUserDataPtr, first);
#else
  bfp->userDataPtr = LinkNewObjectLists (bfp->userDataPtr, first);
#endif
  for (vnp = analysissubmenulist; vnp != NULL; vnp = vnp->next) {
    sbp = (SbstrucPtr) vnp->data.ptrvalue;
    if (sbp != NULL) {
      sbp->name = MemFree (sbp->name);
    }
  }
  analysissubmenulist = ValNodeFreeData (analysissubmenulist);
  ValNodeFreeData (head1);
  ValNodeFreeData (head2);
  return m;
}

