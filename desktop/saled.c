/*   saled.c
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
* File Name:  saled.c
*
* Author:  Colombe Chappey
*
* Version Creation Date:   1/27/96
*
* $Revision: 6.58 $
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
#include <saledit.h>
#include <saled.h>
#include <salutil.h>
#include <salsap.h>
#include <salfiles.h>
#include <salpanel.h>
#include <salstruc.h>
#include <salparam.h>
#include <edutil.h>
#include <dlogutil.h>
#include <aliread.h>
#include <alignmgr2.h>

#define SIDLAND          1
#define SEQLAND          2
#define BADLAND          3
#define HOLDDNLAND       4
#define HOLDUPLAND       5

#define NOCLICK          0
#define CLICKID_BIOSEQ   1
#define CLICKID_SEQFEAT  2
#define CLICKID_FEAT     3
#define CLICKSEQ_BIOSEQ  4
#define CLICKSEQ_SEQFEAT 5
#define CLICKSEQ_FEAT    6

#define MOUSE_NOACTION   0
#define MOUSE_CLICK      1
#define MOUSE_DRAG       2
#define MOUSE_RELEASE    3

/**********************************************************/

static Int4    timerCount = 0;
static Boolean typing_timerInUse   = FALSE;
static Boolean deleting_timerInUse = FALSE;
static Boolean selecting_timerInUse= FALSE;

static void click_caret (RecT *rp, Int4 position, SelStructPtr csp, SelStructPtr caret, EditAlignDataPtr adp, Int2 line);

static Int4 BioseqLength (SeqIdPtr sip)
{ 
  BioseqPtr bsp;
  Int4      position = APPEND_RESIDUE;

  bsp=BioseqLockById(sip);
  if (bsp) {
     position=bsp->length-1;
     BioseqUnlock(bsp);
  }
  return position;
}
/*********************************************************
***
***  ??????
***
**********************************************************/
extern Uint2 OBJ_ (Uint2 feattype)
{
  if ( feattype == FEATDEF_BAD ) return OBJ_BIOSEQ;
  return OBJ_SEQFEAT;
}

/******************************************************************/
static void get_client_rect (PaneL p, RectPtr prc)
{
  ObjectRect (p, prc);
  InsetRect (prc, HRZ_BORDER_WIDTH, VER_BORDER_WIDTH);
}

/******************************************************************/
extern SelStructPtr locate_region (SelStructPtr ssp, Int4 from, Int4 to, SeqIdPtr sip, Uint1 strand, Boolean is_fuzz)
{
  SeqLocPtr slp=NULL;
  SeqIdPtr  siptmp=NULL;
 
  if (sip)
     siptmp = SeqIdDup(sip);
  if (ssp->region != NULL) 
         ssp->region = SeqLocFree ((SeqLocPtr) ssp->region);
  if (is_fuzz)
         slp = fuzz_loc (from, to, strand, siptmp, TRUE, TRUE);
  else { 
         slp = SeqLocIntNew (from, to, strand, siptmp);
  }
  if ( slp != NULL ) {
         ssp->regiontype = OM_REGION_SEQLOC;
         ssp->region = (Pointer) slp;
  }
  else 
         ssp->regiontype = 0;
  if (siptmp)
     SeqIdFree (siptmp);
  return ssp;
}

/******************************************************************/
static SelStructPtr change_to_seqpos (SelStructPtr ssp, SeqAlignPtr salp, ValNodePtr sqloc_list)
{
  SelStructPtr tmp = ssp;
  SeqLocPtr slp;
  SeqIntPtr sit;
  SeqIdPtr  sip;

  slp = (SeqLocPtr) tmp->region;
  if (slp) {
     sit = (SeqIntPtr) slp->data.ptrvalue;
     sip = SeqLocId (slp);
     sit->from = AlignCoordToSeqCoord (SeqLocStart(slp), sip, salp, sqloc_list, 0);
     if (sit->from < 0)
        sit->from = 0;
     sit->to = AlignCoordToSeqCoord (SeqLocStop(slp), sip, salp, sqloc_list, 0);
     if (sit->to < 0)
        sit->to = BioseqLength (sip);
  }
  return tmp;
}

/******************************************************************/
extern SelStructPtr replace_region (SelStructPtr ssp, Uint2 ssp_ed, Uint2 ssp_id, Uint2 ssp_it, Int4 from, Int4 to, SeqIdPtr sip, Uint1 strand, Boolean is_fuzz)
{
  if (from < 0 || sip == NULL) 
         return NULL;
  ssp->entityID = ssp_ed;
  ssp->itemID = ssp_id;
  ssp->itemtype = ssp_it;
  if (ssp->regiontype == OM_REGION_SEQLOC) 
         ssp->regiontype = 0;
  if (ssp->region != NULL) 
         ssp->region = SeqLocFree ((SeqLocPtr) ssp->region);
  ssp = locate_region (ssp, from, to, sip, strand, is_fuzz);
  ssp->next = NULL;
  ssp->prev = NULL;
  return ssp;
}

/******************************************************************/
static Int2 nextlineup (Int2 line, Uint2Ptr itemtype, Uint2Ptr itemsubtype)
{
  Int2 j;
  for (j = line-1; j >= 0; j--) 
  {
         if (itemtype[j] == itemtype[line] 
         && itemsubtype[j] == itemsubtype[line]) break;
  }
  return j;
}
static Int2 nextlinedown (Int2 line, Int2 to, Uint2Ptr itemtype, Uint2Ptr itemsubtype)
{
  Int2 j;
  for (j = line+1; j <= to; j++) {
         if (itemtype[j] == itemtype[line] 
         && itemsubtype[j] == itemsubtype[line]) break;
  }
  if (j > to) return -1;
  return j;
}

static Int2 moveup_scrollbar (PaneL pnl, Int2 adpvoffset, Int2 offset)
{
  if (adpvoffset > offset) 
  {
     adpvoffset -= offset;
     SeqEdSetValueScrollBar (pnl, (Int2)adpvoffset);
  }
  return adpvoffset;
}

static Int2 movedown_scrollbar (PaneL pnl, Int2 adpvoffset, Int2 offset)
{
  adpvoffset += offset;
  SeqEdSetValueScrollBar (pnl, (Int2)(adpvoffset));
  return adpvoffset;
}

/******************************************************************
***
***  GoToButton 
***  LookAtButton
***
*******************************************************************/
static void GoToFunc (SeqEditViewFormPtr wdp)
{
  WindoW             temport;
  PaneL              pnl;
  EditAlignDataPtr   adp;
  RecT               rp;
  Char               str [16];
  Int4               line, column;
  Int4               val; 
  Boolean            goOn;

  if (wdp==NULL)
     return;
  pnl = wdp->pnl;
  if ( ( adp = GetAlignDataPanel (pnl) ) == NULL ) return;
  temport = SavePort (pnl);
  Select (pnl);
 
  GetTitle (wdp->gototxt, str, 15);
  goOn = CCStrToLong (str, &val);
  if (!goOn) {
     RestorePort (temport);
     CaptureSlateFocus ((SlatE) pnl);
     return;
  }
  ResetClip ();
  if (val < 0) {
     val = 0;
     sprintf (str, "%ld", (long) val);
     SetTitle (wdp->gototxt, str);
  }
  else if (val > adp->length) {
     val = adp->length;
     sprintf (str, "%ld", (long) val);
     SetTitle (wdp->gototxt, str);
  }
  if (val < adp->hoffset || val >=  adp->hoffset + adp->visibleLength) 
  {
     adp->voffset = hoffset2voffset(adp, adp->anp_list, adp->visibleWidth, 0, adp->length-1, val);
     data_collect_arrange (adp, TRUE);
     SeqEdSetValueScrollBar (pnl, (Int2) adp->voffset);
  }
  get_client_rect (pnl, &rp);
  SeqPosToLineColumn (adp->caret.itemID, adp->caret.entityID, adp->caret.itemtype, val, &line, &column, adp->hoffset, adp);
  click_caret (&rp, val, &(adp->caret), &(adp->caret), adp, line);
  if (!adp->display_panel)
     to_update_prompt (pnl, &(adp->caret), (SeqAlignPtr) adp->sap_align->data, adp->sqloc_list, FALSE, adp->printid);
  RestorePort (temport);
  CaptureSlateFocus ((SlatE) pnl);
}

extern void GoToButton (ButtoN b)
{
  SeqEditViewFormPtr wdp;

  wdp = (SeqEditViewFormPtr) GetObjectExtra ((WindoW)ParentWindow (b));
  GoToFunc (wdp);
}

extern void LookAtButton (ButtoN b)
{
  WindoW             temport;
  PaneL              pnl;
  SeqEditViewFormPtr wdp;
  EditAlignDataPtr   adp;
  Char               str [16];
  Int4               val;
  Boolean            goOn;

  wdp = (SeqEditViewFormPtr) GetObjectExtra ((WindoW)ParentWindow (b));
  pnl = wdp->pnl;
  if ( ( adp = GetAlignDataPanel (pnl) ) == NULL ) return;
  temport = SavePort (pnl);
  Select (pnl);

  GetTitle (wdp->lookattxt, str, 15);
  goOn = CCStrToLong (str, &val);
  if (!goOn) return;

  ResetClip ();
  if (val < 0) {
     val = 0;
     sprintf (str, "%ld", (long) val);
     SetTitle (wdp->lookattxt, str);
  }
  else if (val > adp->length) {
     val = adp->length;
     sprintf (str, "%ld", (long) val);
     SetTitle (wdp->lookattxt, str);
  }
  if (val < adp->hoffset || val >=  adp->hoffset + adp->visibleLength) 
  {
     adp->voffset = hoffset2voffset(adp, adp->anp_list, adp->visibleWidth, 0, adp->length-1, val);
     data_collect_arrange (adp, TRUE);
     SeqEdSetValueScrollBar (pnl, (Int2) adp->voffset);
  }
  inval_panel (pnl, -1, -1);
  RestorePort (temport);
  CaptureSlateFocus ((SlatE) pnl);
}

static Boolean checkfeatdel (SeqLocPtr slp, ValNodePtr feathead)
{
  ValNodePtr     vnpfeat;
  SelEdStructPtr feat;
  SeqLocPtr      slpfeat;
  Int4           lg;
  Boolean        include;
  Boolean        overlap;

  for  (vnpfeat = feathead; vnpfeat != NULL; vnpfeat = vnpfeat->next) 
  {
     feat = (SelEdStructPtr) vnpfeat->data.ptrvalue;
     while ( feat != NULL ) 
     {
        slpfeat = (SeqLocPtr) feat->region;
        include = include_ssp (slp, slpfeat); 
        if (include) 
           return FALSE;
        overlap = overlapp_ssp (slp, slpfeat);
        if (overlap) {
           lg = 0;
           if (SeqLocStart(slp) > SeqLocStart(slpfeat)) 
              lg += SeqLocStart(slp) - SeqLocStart(slpfeat);
           if (SeqLocStop(slpfeat) > SeqLocStop(slp)) 
              lg += SeqLocStop(slpfeat) - SeqLocStop(slp);
           if (lg < 2) {
              return FALSE;
           }
        }
        feat = feat->next;
     }
  }
  return TRUE;
}

static void CutBtn (ButtoN b)
{
  WindoW             wdialog;
  DialogBoxDataPtr   dbdp;
  EditAlignDataPtr   adp;
  SelStructPtr       ssp; 
  SeqLocPtr          slp;
  BioseqPtr          bsp;

  wdialog = ParentWindow (b);
  Hide (wdialog);
  Update();
  dbdp = (DialogBoxDataPtr) GetObjectExtra (wdialog);
  if ( ( adp = GetAlignEditData (dbdp->w) ) == NULL ) return;
  WatchCursor ();
  ssp = (SelStructPtr) adp->extra_data;
  slp = (SeqLocPtr) ssp->region;
  if (slp!=NULL) {
     bsp = BioseqCopy (NULL, SeqLocId(slp), SeqLocStart(slp), SeqLocStop(slp), Seq_strand_plus, TRUE);
     if (bsp == NULL) {
         ErrPostEx (SEV_ERROR, 0, 0, "fail at BioseqCopy [9]");
     }
     else {
        ObjMgrAddToClipBoard (0, (Pointer) bsp);
        SeqDeleteByLoc (slp, TRUE, adp->spliteditmode);
        adp->dirty = TRUE;
        ObjMgrSetDirtyFlag (ssp->entityID, TRUE);
        ObjMgrSendMsg(OM_MSG_UPDATE, ssp->entityID, ssp->itemID, ssp->itemtype);
        ObjMgrDeSelect(ssp->entityID, ssp->itemID, ssp->itemtype, ssp->regiontype, ssp->region);
     }
     adp->extra_data = NULL;
  }
  ArrowCursor ();
  Remove (wdialog);
  Update ();
}

static void do_cut_dialog (WindoW w)
{
  WindoW           wdialog;
  DialogBoxDataPtr dbdp;
  GrouP            g1;
 
  wdialog = FixedWindow (-50, -33, -10, -10, "Cut", StdCloseWindowProc);
  dbdp = (DialogBoxDataPtr) MemNew (sizeof (DialogBoxData));
  SetObjectExtra (wdialog, (Pointer) dbdp, StdCleanupExtraProc);
  dbdp->w = w;
  g1 = HiddenGroup (wdialog, 2, 0, NULL);
  StaticPrompt (g1, "Cut is going to delete a feature",0,popupMenuHeight,systemFont, 'l');
  g1 = HiddenGroup (wdialog, 2, 0, NULL);
  PushButton (g1, "Process", CutBtn);
  PushButton (g1, "Dismiss", StdCancelButtonProc);
  RealizeWindow (wdialog);
  Show (wdialog);
  return;
}
 
extern Boolean do_cut (PaneL pnl, EditAlignDataPtr adp, SelStructPtr ssp, Boolean cut)
{
  BioseqPtr        bsp;
  SeqLocPtr        slp;
  SeqIntPtr        sit;
  Int4             to;

  StringToClipboard (NULL);
  if (ssp == NULL) {
         ErrPostEx (SEV_ERROR, 0, 0, "Can not cut [1]");
         return FALSE;
  }
  if (ssp->itemtype != OBJ_BIOSEQ) {
         ErrPostEx (SEV_ERROR, 0, 0, "Can not cut [2]");
         return FALSE;
  }
  slp = (SeqLocPtr) ssp->region;
  if (slp == NULL) {
         ErrPostEx (SEV_ERROR, 0, 0, "Can not cut [3]");
         return FALSE;
  }
  if (SeqLocStop(slp) == APPEND_RESIDUE) {
     to = BioseqLength (SeqLocId(slp));
     if (to > 0) {
        sit = (SeqIntPtr) slp->data.ptrvalue;
        sit->to = to; 
     }
     else {
        ErrPostEx (SEV_ERROR, 0, 0, "Can not cut [6]");
        return FALSE;
     }
  }
  if (SeqLocLen(slp) <= 0 || SeqLocLen(slp) >= adp->length) {
         ErrPostEx (SEV_ERROR, 0, 0, "Can not cut [5]");
         return FALSE;
  }

  if ( !checkfeatdel (slp, adp->seqfeat) ) {
         adp->extra_data = (Pointer) ssp;
         do_cut_dialog ((WindoW)ParentWindow(pnl)); 
         return TRUE;
  }
  if (cut) {
     WatchCursor ();
     bsp = BioseqCopy (NULL, SeqLocId(slp), SeqLocStart(slp), SeqLocStop(slp), Seq_strand_plus, TRUE);
     ArrowCursor ();
     if (bsp == NULL) {
         ErrPostEx (SEV_ERROR, 0, 0, "Can not cut [4]");
         return FALSE;
     }
  }
  WatchCursor ();
  if (cut)
     ObjMgrAddToClipBoard (0, (Pointer) bsp);
  SeqDeleteByLoc (slp, TRUE, adp->spliteditmode);
  adp->dirty = TRUE;
  ObjMgrSetDirtyFlag (ssp->entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, ssp->entityID, ssp->itemID, ssp->itemtype);
  ObjMgrDeSelect (ssp->entityID, ssp->itemID, ssp->itemtype, ssp->regiontype, ssp->region);
  ArrowCursor ();
  Update ();
  return TRUE;
}

/******************************************************************/
extern Boolean do_paste (PaneL pnl, EditAlignDataPtr adp, SeqIdPtr sourceid)
{
  SelStructPtr     caret;
  SeqLocPtr        slpcaret;
  Int4             caretpos;
  Boolean          ok;

  caret = &(adp->caret);
  slpcaret = (SeqLocPtr) caret->region;
  if ( slpcaret == NULL ) {
         ErrPostEx (SEV_ERROR, 0, 0, "Fail in do_paste [1]");
         return FALSE;
  }
  caretpos = SeqLocStart (slpcaret);
  if ( caretpos < 0 ) {
         ErrPostEx (SEV_ERROR, 0, 0, "Fail in do_paste [2]");
         return FALSE;
  }
  if ( caretpos == adp->length ) caretpos = -2;

  WatchCursor ();
  ok = BioseqInsert (sourceid, FIRST_RESIDUE, LAST_RESIDUE, Seq_strand_plus, SeqLocId (slpcaret), caretpos, TRUE, TRUE, adp->spliteditmode);
  ArrowCursor ();
  if (!ok) {
         ErrPostEx (SEV_ERROR, 0, 0, "Fail in do_paste [3]");
         return FALSE;
  }
  Update ();
  return TRUE;
}

/******************************************************************/
extern void do_copy (IteM i)
{
  SelStructPtr     ssp = NULL;
  SeqLocPtr        slp;
  BioseqPtr        bsp;
  Int4             from, to;

  WatchCursor ();
  StringToClipboard (NULL);
  if ( checkOMss_for_itemtype (OBJ_BIOSEQ) ) {
     ssp = ObjMgrGetSelected ();
     for (; ssp != NULL; ssp = ssp->next)
        if ( ssp->itemtype == OBJ_BIOSEQ && checkssp_for_editor (ssp) ) 
           break;
     if ( ssp != NULL ) {
        slp = (SeqLocPtr) ssp->region;
        from = SeqLocStart(slp);
        to = SeqLocStop(slp);
        if (to < 0) 
           to = BioseqLength (SeqLocId(slp));
        if (to > 0) {
           bsp=BioseqCopy(NULL, SeqLocId(slp), from, to, Seq_strand_plus, TRUE);
           if (bsp != NULL) 
              ObjMgrAddToClipBoard (0, (Pointer) bsp);
        }
     }
  }
  ArrowCursor ();
  return;
}

static Int4 CleanBufferProc (ValNodePtr bufvnp, Int2 nbseq, Int4 buflength)
{
  TextAlignBufPtr tdp;
  CharPtr PNTR    str;
  ValNodePtr      vnp;
  Boolean         goOn = TRUE;
  Int2            j;

  str = (CharPtr PNTR) MemNew ((size_t) ((nbseq +1) * sizeof(CharPtr)));
  for (vnp =bufvnp, j =0; vnp !=NULL && j <nbseq; vnp =vnp->next, j++) {
         tdp = (TextAlignBufPtr) vnp->data.ptrvalue;
         str[j] = tdp->buf;
  }
  while (goOn) {
         for (j = 0; j < nbseq; j++) {
                if ( *(str[j] +buflength -1) != '-') break;
         }
         goOn = (Boolean) (j == nbseq);
         if (goOn) {
                for (j=0; j<nbseq; j++)
                {
                       *(str[j] +buflength -1) = '\0';
                }
                buflength--;
         }
  }
  MemFree (str);
  return buflength;
}

static Boolean AlignDataGapAddProc (Char ch, Int4 cursorx, ValNodePtr linebuff, Int4 bufferlength, Uint2 entityID, Uint2 itemID)
{
  TextAlignBufPtr  tap;
  ValNodePtr       vnp;
  CharPtr          buf;
  Boolean          insert = FALSE;
  Int4             k;
 
  for (vnp = linebuff; vnp != NULL; vnp = vnp->next)
  {
         tap = (TextAlignBufPtr) vnp->data.ptrvalue;
         if (tap != NULL) {
            if (OBJ_(tap->feattype) == OBJ_BIOSEQ) {
                buf = (CharPtr) tap->buf;
                if((tap->seqEntityID == entityID && tap->bsp_itemID == itemID)
                || is_selectedbyID(tap->seqEntityID,tap->bsp_itemID,OBJ_BIOSEQ))
                {
                      for (k = bufferlength; k > cursorx; k--)
                         buf [k] = buf [k-1];
                      buf [cursorx] = ch;
                      buf [bufferlength +1] = '\0';
                      insert = TRUE;
                }
                else {
                      buf [bufferlength] = '-';
                      buf [bufferlength +1] = '\0';
                }
            }
         }       
  } 
  return insert;
}

/******************************************************************
***   
***   CLICK - DRAG - RELEASE procedure
***   
******************************************************************/
static void update_prompt (PaneL pnl, CharPtr label, Int4 position)
{
  WindoW             w;
  SeqEditViewFormPtr wdp;
  SelStructPtr       ssp;
  Char               str[128];
  Int4               from, to;
  Int2               k;

  w = getwindow_frompanel (pnl);
  wdp = (SeqEditViewFormPtr) GetObjectExtra (w);
  ResetClip ();
  if (position >= 0) {
     sprintf (str, "%s Position %ld", label, (long)position);
     SetTitle (wdp->pos, str);
  }
  else if (position == LAST_RESIDUE) {
     sprintf (str, "%s", label);
     SetTitle (wdp->pos, str);
  }
  else {
     k = checkOMss_for_itemtype (OBJ_BIOSEQ);
     if (k == 1) {
        ssp = getOMselect_for_itemtype (OBJ_BIOSEQ);
        from = SeqLocStart ((SeqLocPtr) ssp->region);
        to = SeqLocStop ((SeqLocPtr) ssp->region);
        sprintf (str, "%s Selection %ld - %ld", label, (long)(from+1), (long)(to+1));
        SetTitle (wdp->pos, str);
     }
     else {
        sprintf (str, " ");
        SetTitle (wdp->pos, str);
     }
  }
}

static void update_promptsel (PaneL pnl, CharPtr label)
{
  WindoW             w;
  SeqEditViewFormPtr wdp;
 
  w = getwindow_frompanel (pnl);
  wdp = (SeqEditViewFormPtr) GetObjectExtra (w);
  ResetClip ();
  SetTitle (wdp->pos, label);
}

static CharPtr seqid_tolabel (SeqIdPtr sip, Uint1 choice)
{
  BioseqPtr bsp;
  SeqIdPtr  tmp = NULL;
  Char      str[120];
  CharPtr   strp;

  str[0] = '\0';
  if (choice && choice <= PRINTID_GIcc && sip)
  {
     bsp = BioseqLockById(sip);
     if (bsp!=NULL)
     {
        if (choice==PRINTID_GIcc)
           tmp = SeqIdFindBest (bsp->id, 0);
        else
           tmp = SeqIdFindWorst (bsp->id);
        if (tmp)
        {
           SeqIdWrite (tmp, str, choice, 120);
           BioseqUnlock(bsp);
           strp = StringSave (str);
           return strp;
        }
        BioseqUnlock(bsp);
     }
  }
  if (choice!=0 && sip) {
     SeqIdWrite (sip, str, choice, 120);
     strp = StringSave (str);
     return strp;
  }
  strp = StringSave (str);
  return strp;
}

extern void to_update_prompt (PaneL pnl, SelStructPtr ssp, SeqAlignPtr salp, ValNodePtr sqloc_list, Boolean sel, Uint1 choice)
{
  SeqIdPtr  sip;
  Int4      pos, posto;
  CharPtr   seqid=NULL;
  Char      str[255];
  Char      tmp[252];
  CharPtr   strp;
  SeqLocPtr slp;
  SelStructPtr selp;

  if (ssp==NULL)
     return;
  if (ssp->region == NULL)
     return;

  str[0] = '0';
  slp = (SeqLocPtr) ssp->region;
  sip = SeqLocId(slp);
  if (choice == 0)
     choice = PRINTID_FASTA_LONG;  /*PRINTID_TEXTID_ACCESSION;*/
  seqid = seqid_tolabel (sip, choice);
  selp = is_selectedbyID (ssp->entityID, ssp->itemID, ssp->itemtype);
  if (selp != NULL) 
  {
     strp = str;
     pos = SeqLocStart((SeqLocPtr)selp->region);
     posto = SeqLocStop((SeqLocPtr)selp->region);
     if (posto == APPEND_RESIDUE)
        posto = BioseqLength (sip);
     if (selp->next == NULL)
        sprintf (tmp, "%s Selection %ld - %ld", seqid, (long)(pos+1), (long)(posto+1));
     else
        sprintf (tmp, "%s Selections %ld - %ld", seqid, (long)(pos+1), (long)(posto+1));
     strp = StringMove (strp, tmp);
     selp = selp->next;
     while (selp != NULL && (StringLen(str) < 200)) {
        if (checkssp_for_editor (selp) && is_sameId (selp->entityID, selp->itemID, selp->itemtype, 255, ssp->entityID, ssp->itemID, ssp->itemtype, 255))
        {
           pos = SeqLocStart((SeqLocPtr)selp->region);
           posto = SeqLocStop((SeqLocPtr)selp->region);
           sprintf (tmp, ", %ld - %ld", (long)(pos+1), (long)(posto+1));
           strp = StringMove (strp, tmp);
        }
        selp = selp->next;
     }
     update_promptsel (pnl, str); 
  }
  else {
     if (salp != NULL) {
        pos = AlignCoordToSeqCoord (SeqLocStart(slp), sip, salp, sqloc_list, 1);
        if (pos == APPEND_RESIDUE)
           pos = BioseqLength (sip);
     }
     else pos = SeqLocStart(slp);
     if (sel) {
        if (salp != NULL) {
           posto = AlignCoordToSeqCoord (SeqLocStop(slp), sip, salp, sqloc_list, 0);
        }
        else posto = SeqLocStop(slp);
        if (posto == APPEND_RESIDUE) 
           posto = BioseqLength (sip);
        sprintf (str, "%s Selection %ld - %ld", seqid, (long)(pos+1), (long)(posto+1));
        update_promptsel (pnl, str);
     }
     else { 
        update_prompt (pnl, seqid, pos);
     }
  }
  if (seqid)
     MemFree (seqid);
  return;
}

extern void update_edititem (PaneL pnl)
{
  WindoW             w;
  EditAlignDataPtr   adp;
  SeqEditViewFormPtr wdp;
  Int2               k;

  w = getwindow_frompanel (pnl);
  if ((adp = GetAlignEditData (w)) == NULL) 
     return;
  if (adp->edit_mode == SEQ_EDIT || adp->edit_mode == ALIGN_EDIT) 
  {
     wdp = (SeqEditViewFormPtr) GetObjectExtra (w);
     ResetClip ();
     k = checkOMss_for_itemtype (OBJ_BIOSEQ);
     if (k > 0) {
        Enable (wdp->cutitem);
        Enable (wdp->copyitem);
     }
     else {
        Disable (wdp->cutitem);
        Disable (wdp->copyitem);
     }
  }
}

static void update_savefeatitem (PaneL pnl)
{
  WindoW             w;
  EditAlignDataPtr   adp;
  SeqEditViewFormPtr wdp;
  Int2               k;

  w = getwindow_frompanel (pnl);
  if ((adp = GetAlignEditData (w)) == NULL) 
     return;
  if (adp->edit_mode == SEQ_EDIT || adp->edit_mode == ALIGN_EDIT) {
     wdp = (SeqEditViewFormPtr) GetObjectExtra (w);
     ResetClip ();
     k = checkOMss_for_itemtype(OBJ_SEQFEAT);
     if (k > 0) {
        Enable (wdp->savefeatitem);
        Enable (wdp->savefeatbt);
     }
     else {
        Disable (wdp->savefeatitem);
        Disable (wdp->savefeatbt);
     }
  }
}

extern void update_translateitem (PaneL pnl, ValNodePtr seqfeathead, ValNodePtr feathead)
{
  WindoW             w;
  SeqEditViewFormPtr wdp;
  Int2               k;

  w = getwindow_frompanel (pnl);
  wdp = (SeqEditViewFormPtr) GetObjectExtra (w);
  ResetClip ();
  k = checkCDSselect_forprotein (seqfeathead, feathead, FALSE);
  if (k == 0) {
     Disable (wdp->translateitem);
     Disable (wdp->translatebt);
  }
  else {
     Enable  (wdp->translateitem);
     Enable  (wdp->translatebt);
  }
}

extern void update_codonstartbt (PaneL pnl, ValNodePtr seqfeathead, ValNodePtr feathead)
{
  WindoW             w;
  SeqEditViewFormPtr wdp;
  SelEdStructPtr     feat;
  Char               str[128];
  Int2               k;

  w = getwindow_frompanel (pnl);
  wdp = (SeqEditViewFormPtr) GetObjectExtra (w);
  ResetClip ();
  k = checkCDSselect_forprotein (seqfeathead, feathead, TRUE);
  if (k == 1) {
     feat = getCDSselect (seqfeathead, feathead);
     sprintf (str, "Codon start %d", (int)(feat->codonstart));
     SetTitle (wdp->codonstitem, str);
     Enable  (wdp->codonstitem);
  }
  else {
     sprintf (str, "Codon start");
     SetTitle (wdp->codonstitem, str);
     Disable (wdp->codonstitem);
  }
}

extern void on_time (WindoW w)
{
  EditAlignDataPtr adp;

  if ((adp = GetAlignEditData (w)) == NULL) return; 
  timerCount++;
  if (typing_timerInUse || deleting_timerInUse) {
     if (adp->edit_item.entityID != 0 && timerCount > 40) {
        typing_timerInUse = FALSE;
        deleting_timerInUse = FALSE;
        timerCount = 0;
        ObjMgrSendMsg (OM_MSG_UPDATE, adp->caret.entityID, adp->caret.itemID, adp->caret.itemtype);
     }
  }
  else if (selecting_timerInUse) {
     if (timerCount > 80) {
        selecting_timerInUse = FALSE;  
        timerCount = 0;
        ObjMgrSendMsg (OM_MSG_SELECT, adp->caret.entityID, adp->caret.itemID, adp->caret.itemtype);
     }
  }
  return;
}

static void update_select_when_release (EditAlignDataPtr adp)
{
   if (selecting_timerInUse) {
     selecting_timerInUse = FALSE;  
     timerCount = 0;
     ObjMgrSendMsg (OM_MSG_SELECT, adp->caret.entityID, adp->caret.itemID, adp->caret.itemtype);
  }
  return;
}

static void click_id (PaneL pnl, SelStructPtr csp, EditAlignDataPtr adp)
{
  SelStructPtr     spp;
  SeqIdPtr         sip;
  Int4             from, to;
  Int4             start, stop;
  ValNodePtr       vnp;
  AlignNodePtr     anp;
  SeqLocPtr        slp;
  SeqAlignPtr      salp = (SeqAlignPtr) adp->sap_align->data;

  if ( csp == NULL ) {
     ErrPostEx (SEV_ERROR, 0, 0, "fail in click_id [1]");
     return;
  }
  Select (pnl);
  spp = is_selectedbyID (csp->entityID, csp->itemID, csp->itemtype); 
  if (spp != NULL) {
     ObjMgrDeSelectAll ();
     if (! adp->display_panel)
        to_update_prompt (pnl, &(adp->caret), salp, adp->sqloc_list, FALSE, adp->printid); 
  } 
  else if (shftKey && ctrlKey) 
  {
     spp=ObjMgrGetSelected();
     if(spp!=NULL) {
          for(; spp!=NULL; spp=spp->next){
             if(spp->entityID==csp->entityID && spp->itemID==csp->itemID)
                break;
          }
          if (spp == NULL) {
             spp=ObjMgrGetSelected();
             slp=(SeqLocPtr)spp->region;
             start=SeqCoordToAlignCoord(SeqLocStart(slp), SeqLocId(slp), salp, 0, 0); 
             stop=SeqCoordToAlignCoord(SeqLocStop(slp), SeqLocId(slp), salp, 0, 0); 
             sip = SeqLocId ((SeqLocPtr)csp->region);
             from=AlignCoordToSeqCoord (start, sip, salp, adp->sqloc_list, 0);
             to=AlignCoordToSeqCoord (stop, sip, salp, adp->sqloc_list, 0);
             slp = SeqLocIntNew(from, to, SeqLocStrand((SeqLocPtr)csp->region), sip);  
             ObjMgrAlsoSelect (csp->entityID, csp->itemID, csp->itemtype, OM_REGION_SEQLOC, (Pointer)slp);
          }
     }
  } 
  else if (shftKey) 
  {
     csp = change_to_seqpos (csp, salp, adp->sqloc_list);
     for(spp=ObjMgrGetSelected(); spp!=NULL; spp=spp->next) 
           if (spp->entityID!=csp->entityID || spp->itemID!=csp->itemID)
              break;
     if (spp == NULL) 
     {
           ObjMgrAlsoSelect (csp->entityID, csp->itemID, csp->itemtype, csp->regiontype, csp->region); 
     } 
     else {
           for (vnp=adp->anp_list; vnp != NULL; vnp = vnp->next) {
              anp = (AlignNodePtr) vnp->data.ptrvalue; 
              if (anp != NULL)
                 if (anp->seq_entityID==spp->entityID && anp->bsp_itemID == spp->itemID )
                   break;
           }
           if (vnp!=NULL) {
              vnp = vnp->next;
              for (; vnp!=NULL; vnp = vnp->next) {
                 anp = (AlignNodePtr) vnp->data.ptrvalue;
                 if (anp!=NULL) {
                    spp = is_selectedbyID (anp->seq_entityID, anp->bsp_itemID, OBJ_BIOSEQ);
                    if (spp==NULL) {
                       slp = CollectSeqLocFromAlignNode(anp);
                       if (slp!=NULL) 
                       {
                          ObjMgrAlsoSelect (anp->seq_entityID, anp->bsp_itemID, OBJ_BIOSEQ, OM_REGION_SEQLOC, (Pointer)slp);
                       }
                       if (anp->seq_entityID==csp->entityID && anp->bsp_itemID == csp->itemID )             
                          break;
                    }
                 }
              }
           }
     }
  } 
  else if (ctrlKey) 
  {
     csp = change_to_seqpos (csp, salp, adp->sqloc_list);
     ObjMgrAlsoSelect (csp->entityID, csp->itemID, csp->itemtype, csp->regiontype, csp->region); 
  } 
  else {
     csp = change_to_seqpos (csp, salp, adp->sqloc_list);
     ObjMgrSelect (csp->entityID, csp->itemID, csp->itemtype, csp->regiontype, csp->region);
  } 
  Update ();
  return;
}

/******************************************************************/
static void click_caret (RecT *rp, Int4 position, SelStructPtr csp, SelStructPtr caret, EditAlignDataPtr adp, Int2 line)
{
  SeqIdPtr sip;

  if (csp != NULL) {
     if ( caret->regiontype != 0 )
        inval_selstructpos(adp,caret->entityID, caret->itemID, caret->itemtype, rp, SeqLocStart ((SeqLocPtr)caret->region));
     sip = SeqIdDup (SeqLocId((SeqLocPtr)csp->region));
     replace_region (caret, csp->entityID, csp->itemID, csp->itemtype, position, position, sip, SeqLocStrand((SeqLocPtr)csp->region), FALSE);
     SeqIdFree (sip);
     adp->caret_orig = position;
     adp->cur_pat = NULL;
     inval_selstructpos (adp, caret->entityID, caret->itemID, caret->itemtype, rp, position);
  }
}

/******************************************************************/
static void click_feat (PaneL pnl, Int4 position, SelStructPtr csp, EditAlignDataPtr adp, Int2 line)
{
  RecT             rp;
  SeqLocPtr        slp;
  SeqAlignPtr      salp;
  Int4             start, stop;
  Int2             width;

  if ( csp == NULL ) {
     ErrPostEx (SEV_ERROR, 0, 0, "fail in click_feat [1]");
     return;
  }
  get_client_rect (pnl, &rp);
  if ( ! shftKey ) {
     width = adp->visibleWidth;
     if (adp->columnpcell > 0) 
        width += (Int2) adp->visibleWidth / (Int2) adp->columnpcell;
  }
  salp = (SeqAlignPtr) adp->sap_align->data;
  slp = (SeqLocPtr) csp->region;
  start = SeqLocStart(slp);
  stop = SeqLocStop(slp);
  start = SeqCoordToAlignCoord (start, SeqLocId(slp), salp, 0, 0);
  stop = SeqCoordToAlignCoord (stop, SeqLocId(slp), salp, 0, 0);
  if (position > start && position < stop) 
     adp->click_feat = 3;
  else if (position <= start && position > start- 2) 
     adp->click_feat = 1;
  else if (position >= stop  && position < stop + 2) 
     adp->click_feat = 2;
  else adp->click_feat = 0;
  adp->feat_line = line;
  inval_selstructloc_forfeat (adp, csp->entityID, csp->itemID, csp->itemtype, 255, &rp, start, stop);
  csp = change_to_seqpos (csp, salp, adp->sqloc_list);
  if ( ! shftKey ) 
  {
     ObjMgrSelect (csp->entityID,csp->itemID, csp->itemtype, csp->regiontype, csp->region);
  } 
  else 
  { 
     ObjMgrAlsoSelect (csp->entityID, csp->itemID, csp->itemtype, csp->regiontype, csp->region);
  }
}

/******************************************************************/
static void invalminmax (RecT *rp, SelStructPtr spp, Int4 position, EditAlignDataPtr adp)
{
  SeqLocPtr  slp;
  Int4       pos, 
             posmin, posmax;

  slp = (SeqLocPtr) spp->region;
  if ( position > adp->caret_orig ) 
  {
     pos = MIN((Int4)SeqLocStart(slp),(Int4)SeqLocStop(slp));
     posmin = MIN(pos,(Int4) adp->edit_pos);
     posmax = MAX(pos,(Int4) adp->edit_pos);
     inval_selstructloc (adp, spp->entityID, spp->itemID, spp->itemtype, 255, rp, posmin, posmax);
  } 
  else if ( position < adp->caret_orig ) 
  {
     pos = MAX((Int4)SeqLocStart(slp),(Int4)SeqLocStop(slp));
     posmin = MIN(pos,(Int4) adp->edit_pos);
     posmax = MAX(pos,(Int4) adp->edit_pos);
     inval_selstructloc (adp, spp->entityID, spp->itemID, spp->itemtype, 255, rp, posmin, posmax);
  }
}


/******************************************************************/
static SelStructPtr get_closest_selection (SelStructPtr caret, SeqAlignPtr salp, SelStructPtr current_point, ValNodePtr sqloc_list)
{
  SelStructPtr ssp,
               car,
               this_selection=NULL;
  Int4         position;
  Int4         val , valmin = INT4_MAX;

  car = SelStructDup (caret);
  car = change_to_seqpos (car, salp, sqloc_list);
  if (car)
  {
     ssp = ObjMgrGetSelected();
     if (ssp==NULL)
        this_selection = current_point;
     else 
     {
        position = SeqLocStart((SeqLocPtr)car->region);
        for (; ssp != NULL; ssp = ssp->next)
        {
           if ( checkssp_for_editor (ssp) && is_sameId (ssp->entityID, ssp->itemID, ssp->itemtype, 255, car->entityID, car->itemID, car->itemtype, 255) )
           {
              val = (Int4)(ABS(SeqLocStart((SeqLocPtr)ssp->region)-position)+ABS(SeqLocStop((SeqLocPtr)ssp->region)-position));

              if ( val < valmin ) 
              {
                 this_selection = ssp;
                 valmin = val;
              }
           }
        }
     }
     SelStructDel (car);
  }
  return this_selection;
}


/******************************************************************
*** 
*** Open and Extend this_selection by AlsoSeqlect
***
**************** *******************************************/
static Boolean extend_selectregion (EditAlignDataPtr adp, SelStructPtr this_selection, Int4 position,Int4 old_caret_pos, Boolean IsOtherSel, PaneL pnl)
{
  SeqIdPtr         sip;
  SeqLocPtr        slp;
  SeqAlignPtr      salp;
  Int4             from, to, itmp;  
  Uint2            sspei, sspii, sspit;
  Boolean          retval=FALSE;

  sspei = (Uint2) adp->caret.entityID;
  sspii = (Uint2) adp->caret.itemID;
  sspit = (Uint2) adp->caret.itemtype;
  if (is_sameId (this_selection->entityID, this_selection->itemID, this_selection->itemtype, 255, sspei, sspii, sspit, 255))
  {  
     salp = (SeqAlignPtr) adp->sap_align->data;
	 sip = SeqLocId((SeqLocPtr)adp->caret.region);	
	 if (!IsOtherSel){/*no previous close selection*/
	    		   /*Create a SelStruct between old_caret_pos and position*/
		from = AlignCoordToSeqCoord (old_caret_pos, sip, salp, adp->sqloc_list, 0);
		to = AlignCoordToSeqCoord (position, sip, salp, adp->sqloc_list, 0);
		if (from > to){
			itmp=to; to=from; from=itmp;
		}
		slp = SeqLocIntNew (from, to, Seq_strand_plus, sip);
		retval=ObjMgrAlsoSelect (this_selection->entityID, this_selection->itemID, 
                   this_selection->itemtype, OM_REGION_SEQLOC, slp);
	 }
 } 
  return retval;
}

/******************************************************************/
static void dragcaret_selectrgn (PaneL pnl, Int4 position, SelStructPtr spp, SelStructPtr caret, EditAlignDataPtr adp, Int2 line)
{
  SeqLocPtr        slp;
  SeqAlignPtr      salp;
  SeqIdPtr         sip;
  RecT             rp;
  float hratio;
  
  if ( spp == NULL)  {
/*********
         ErrPostEx (SEV_ERROR, 0, 0, "fail in dragcaret_selectrgn [2]");
*********/
         return;
  }
  get_client_rect (pnl, &rp);
  hratio = (float)adp->hoffset / (float)adp->length;
  salp = (SeqAlignPtr) adp->sap_align->data;
  selecting_timerInUse = TRUE;
  if ( position == adp->caret_orig ) 
  {
     invalminmax (&rp, spp, position, adp);
     inval_selstructloc (adp, spp->entityID, spp->itemID, spp->itemtype, 255, &rp, SeqLocStart((SeqLocPtr) spp->region), SeqLocStop((SeqLocPtr) spp->region));
     adp->caret_line = line;
     data_collect_arrange (adp, TRUE);
     SeqEdSetCorrectBarMax (pnl, adp, hratio);
     if (!adp->display_panel)
         to_update_prompt (pnl, spp, NULL, adp->sqloc_list, FALSE, adp->printid);
     ObjMgrDeSelect (spp->entityID, spp->itemID, spp->itemtype, spp->regiontype, spp->region);
  } 
  else {
     slp = (SeqLocPtr) spp->region;
     sip = SeqLocId (slp);
     invalminmax (&rp, spp, position, adp);
     position = AlignCoordToSeqCoord (position, sip, salp, adp->sqloc_list, 0);
     if ( position > adp->caret_orig ) 
     {
        if (position < 0) position =BioseqLength (sip);
        slp = SeqLocIntNew (SeqLocStart(slp), position, Seq_strand_plus, sip);
     } 
     else 
     {
        if (position < 0) position = 0;
        slp = SeqLocIntNew (position, SeqLocStop(slp), Seq_strand_plus, sip);
     }
     ObjMgrAlsoSelect (spp->entityID, spp->itemID, spp->itemtype, OM_REGION_SEQLOC, slp);
     adp->caret_line = line;
     if (!adp->display_panel)
         to_update_prompt (pnl, spp, salp, adp->sqloc_list, TRUE, adp->printid);
  }
  return;
}



/******************************************************************/
static void drag_feat (PaneL pnl, Int4 position, SelEdStructPtr feat, EditAlignDataPtr adp, Int2 line, Int2 feattype, SelEdStructPtr prec, SelEdStructPtr next, SelStructPtr selssp)
{
  SeqLocPtr        slp;
  SeqLocPtr        slpprec = NULL, slpnext = NULL;
  SeqIntPtr        sitfeat;
  SeqAlignPtr      salp;
  RecT             rp;
  Int4             linebefmv, lineaftmv, line2befmv, line2aftmv;
  Int4             column;
  Int2             width;
  Int4             oldposfrom, oldposto;
  Int4             oldseqfrom, oldseqto;
  Int4             minpos, maxpos;
  Boolean          updateline = FALSE;
  SeqLocPtr selslp;
  SeqIntPtr selsit;
  float hratio;
  
  if ( adp->click_feat == 0 ) return;
  if ( feat == NULL)  {
         ErrPostEx (SEV_ERROR, 0, 0, "fail in drag_feat [0]");
         return;
  }
  get_client_rect (pnl, &rp);
  width = adp->visibleWidth;
  if (adp->columnpcell > 0) 
         width += (Int2) adp->visibleWidth / (Int2) adp->columnpcell;
  slp = (SeqLocPtr) feat->region;
  if (slp == NULL) {
         ErrPostEx (SEV_ERROR, 0, 0, "fail in drag_feat [2.1]");
         return;
  }
  sitfeat = (SeqIntPtr) slp->data.ptrvalue;
  if (sitfeat == NULL) {
         ErrPostEx (SEV_ERROR, 0, 0, "fail in drag_feat [2.2]");
         return;
  }
  hratio = (float)adp->hoffset / (float)adp->length;
  salp = (SeqAlignPtr) adp->sap_align->data;
  oldseqfrom = SeqLocStart (slp);
  oldseqto = SeqLocStop (slp);
  oldposfrom = SeqCoordToAlignCoord(oldseqfrom, SeqLocId(slp), salp, 0, 0);
  oldposto = SeqCoordToAlignCoord(oldseqto, SeqLocId(slp), salp, 0, 0); 
  if (prec != NULL) 
     slpprec = (SeqLocPtr) prec->region;
  if (next != NULL) 
     slpnext = (SeqLocPtr) next->region;
  if ( adp->click_feat == 1 ) 
  {
      if ( position >= 0 && position <= oldposto) 
      {
         sitfeat->from=AlignCoordToSeqCoord(position,SeqLocId(slp), salp, adp->sqloc_list, 0);
         if (!overlapp_ssp (slp, slpprec) && !overlapp_ssp (slp, slpnext)) 
         {
            if ((feattype == SEQFEAT_CDREGION || feattype == FEATDEF_CDS) 
            && feat->data != NULL ) {
               sesp_to_pept (feat, salp, adp->sqloc_list, FALSE);
               inval_selstruct (adp, feat->entityID, feat->itemID, feat->itemtype, FEATDEF_TRSL, &rp, adp->margin.left, (Int2) (width *adp->charw));
            }
            minpos = MIN (oldposfrom, position);
            maxpos = MAX (oldposfrom, position) + 1;
            inval_selstructloc_forfeat (adp, feat->entityID, feat->itemID, feat->itemtype, feattype, &rp,  minpos, maxpos);
            updateline = TRUE;
         }
         else sitfeat->from = oldseqfrom;
      }
  }
  else if ( adp->click_feat == 2 ) 
  {
      if ( position < adp->length && position >= oldposfrom) 
      {
         sitfeat->to =AlignCoordToSeqCoord(position,SeqLocId(slp), salp, adp->sqloc_list, 0);
         if (!overlapp_ssp (slp, slpprec) && !overlapp_ssp (slp, slpnext)) 
         {
            if ((feattype == SEQFEAT_CDREGION || feattype == FEATDEF_CDS) 
            && feat->data != NULL ) {
               sesp_to_pept (feat, salp, adp->sqloc_list, FALSE);
               inval_selstruct (adp, feat->entityID, feat->itemID, feat->itemtype, FEATDEF_TRSL, &rp, adp->margin.left, (Int2) (width *adp->charw));
            }
            minpos = MIN (oldposto, position) - 1;
            maxpos = MAX (oldposto, position);
            inval_selstructloc_forfeat (adp, feat->entityID, feat->itemID, feat->itemtype, feattype, &rp, minpos, maxpos);
            updateline = TRUE;
         }
         else sitfeat->to = oldseqto;
      }
  }
  else if ( adp->click_feat == 3 ) 
  {
      if (oldposfrom >= adp->feat_pos - position
      && oldposto + position - adp->feat_pos < adp->length) 
      {
         if(position < adp->feat_pos && oldposfrom > adp->feat_pos -position) 
         {
                sitfeat->from -= (adp->feat_pos - position);
                sitfeat->to -= (adp->feat_pos - position);
         } 
         else if (oldposto + adp->feat_pos - position < adp->length) 
         {
                sitfeat->from += (position - adp->feat_pos);
                sitfeat->to += (position - adp->feat_pos);
         }
         if (!overlapp_ssp (slp, slpprec) && !overlapp_ssp (slp, slpnext)) 
         {
            if ((feattype == SEQFEAT_CDREGION || feattype == FEATDEF_CDS) 
            && feat->data != NULL ) {
               sesp_to_pept (feat, salp, adp->sqloc_list, FALSE);
               inval_selstruct (adp, feat->entityID, feat->itemID, feat->itemtype, FEATDEF_TRSL, &rp, adp->margin.left, (Int2) (width *adp->charw));
            }
            maxpos = SeqCoordToAlignCoord(sitfeat->from, SeqLocId(slp), salp, 0, 0);
            minpos = MIN (oldposfrom, maxpos);
            maxpos = MAX (oldposfrom, sitfeat->from) + 1;
            inval_selstructloc_forfeat (adp, feat->entityID, feat->itemID, feat->itemtype, feattype, &rp, minpos, maxpos);
            maxpos = SeqCoordToAlignCoord(sitfeat->to, SeqLocId(slp), salp, 0, 0);
            minpos = MIN (oldposto, maxpos) - 1;
            maxpos = MAX (oldposto, maxpos);
            inval_selstructloc_forfeat (adp, feat->entityID, feat->itemID, feat->itemtype, feattype, &rp, minpos, maxpos);
            updateline = TRUE;
         }
         else {
            sitfeat->from = oldseqfrom;
            sitfeat->to = oldseqto;
         }
      }
  }
  selslp = (SeqLocPtr) selssp->region;
  selsit = (SeqIntPtr) selslp->data.ptrvalue;
  selsit->from = sitfeat->from;
  selsit->to = sitfeat->to;

  if (updateline)
  {
     SeqPosToLineColumn (feat->itemID, feat->entityID, feat->itemtype,  
         oldposfrom, &linebefmv, &column, adp->hoffset, adp);
     SeqPosToLineColumn (feat->itemID, feat->entityID, feat->itemtype, 
         oldposto, &line2befmv, &column, adp->hoffset, adp);
     SeqPosToLineColumn (feat->itemID, feat->entityID, feat->itemtype, 
         SeqLocStart(slp), &lineaftmv, &column, adp->hoffset, adp);
     SeqPosToLineColumn (feat->itemID, feat->entityID, feat->itemtype, 
         SeqLocStop(slp), &line2aftmv,&column, adp->hoffset, adp);
     if ( lineaftmv > linebefmv ) 
     {
         ResetClip ();
         inval_rect (rp.left, (Int2)(rp.top+ linebefmv *adp->lineheight),  
                     rp.right, rp.bottom);
         data_collect_arrange (adp, TRUE);
         SeqEdSetCorrectBarMax (pnl, adp, hratio);
         adp->voffset = moveup_scrollbar (pnl, adp->voffset, 1);
         Update ();

     }
     else if ( lineaftmv < linebefmv ) 
     {
         ResetClip ();
         inval_rect (rp.left, (Int2)(rp.top+ lineaftmv *adp->lineheight),  
                          rp.right, rp.bottom);
         data_collect_arrange (adp, TRUE);
         SeqEdSetCorrectBarMax (pnl, adp, hratio);
         adp->voffset = movedown_scrollbar (pnl, adp->voffset, (Int2)(1));
         Update ();
     }
     if ( line2aftmv < line2befmv ) {
         inval_rect (rp.left, (Int2)(rp.top+ line2befmv *adp->lineheight),  
                          rp.right, rp.bottom); 
         data_collect_arrange (adp, TRUE);
         SeqEdSetCorrectBarMax (pnl, adp, hratio);
     }
     else if ( line2aftmv > line2befmv ) {
         inval_rect (rp.left, (Int2)(rp.top+ line2aftmv *adp->lineheight),  
                          rp.right, rp.bottom);
         data_collect_arrange (adp, TRUE);
         SeqEdSetCorrectBarMax (pnl, adp, hratio);
     } 
     adp->feat_line = line;
  } 
  adp->dirty = TRUE;
  return;
}

/******************************************************************/
static void release_caret (PaneL pnl, Int4 position, SelStructPtr csp, SelStructPtr caret, EditAlignDataPtr adp, Int2 line)
{
  RecT             rp;

  if ( csp == NULL ) {
         return;
  }
  if  ( adp->caret.regiontype == 0 || adp->caret.region == NULL)  {
/*****
         ErrPostEx (SEV_ERROR, 0, 0, "fail in release_caret [2]");
*****/
         return;
  }
  if ( !is_samessp (caret, csp) ) {
         get_client_rect (pnl, &rp);
         click_caret (&rp, position, csp, caret, adp, line);
         if (!adp->display_panel)
            to_update_prompt (pnl, &(adp->caret), (SeqAlignPtr) adp->sap_align->data, adp->sqloc_list, FALSE, adp->printid);
         ObjMgrDeSelectAll ();
  }
  update_select_when_release (adp);
}

/******************************************************************/
static SelStructPtr make_bspssp (SelStructPtr csp, ValNodePtr anp_list, Int4 from, Int4 to)
{
  SeqIdPtr         sip;
  Uint1            strand;

  sip = (SeqIdPtr) SeqIdFromAlignNode (anp_list, csp->entityID, csp->itemID, 
         csp->itemtype);
  if (sip == NULL) return NULL;
  strand = StrandFromAlignNode (anp_list, csp->entityID, csp->itemID,
         csp->itemtype);
  locate_region (csp, from, to, sip, strand, FALSE);

  return csp;
}

/******************************************************************/
static SeqIdPtr get_featId_fromid (SelStructPtr csp, ValNodePtr feathead, Int2 itemsubtype)
{
  SelEdStructPtr   feat;
  SeqLocPtr        slp;

  feat = get_feat_fromid (feathead, itemsubtype, csp->entityID, csp->itemID, -1, NULL);
  if (feat == NULL) 
     return NULL;
  slp = (SeqLocPtr) feat->region;
  if (slp == NULL) 
     return NULL;
  return SeqLocId(slp);
}
/******************************************************************/
static SelStructPtr make_featssp (SelStructPtr csp, ValNodePtr feathead, Int2 itemsubtype, Int4 position, Int4 from, Int4 to)
{
  SelEdStructPtr   feat;
  SeqLocPtr        slp;

  feat = get_feat_fromid (feathead, itemsubtype, csp->entityID, csp->itemID, position, NULL);
  if (feat == NULL ) return NULL;
  slp = (SeqLocPtr) feat->region;
  if (slp == NULL) return NULL;
  if (from < 0) { 
      locate_region (csp, SeqLocStart(slp), SeqLocStop(slp), SeqLocId (slp),
                    SeqLocStrand(slp), FALSE);
  }
  else 
      locate_region (csp, from, to, SeqLocId(slp), SeqLocStrand(slp), FALSE);
  return csp;
}

/******************************************************************/
static SelStructPtr make_bufssp (SelStructPtr csp, SelStructPtr bufhead, Int2 itemsubtype)
{
  SelStructPtr     buf;
  SelStructPtr     ssp;
  SeqLocPtr        slp = NULL;
  Uint2            subtype, type;
  Boolean          locate = FALSE;

  if (itemsubtype==SEQFEAT_CDREGION) {
         subtype = SEQFEAT_CDREGION; type = OBJ_SEQFEAT; 
  }
  else if (itemsubtype==FEATDEF_CDS) {
         subtype = FEATDEF_CDS; type = OBJ_SEQFEAT; 
  }
  else {
         subtype = FEATDEF_BAD; type = OBJ_BIOSEQ; 
  }
  for (buf = bufhead; buf != NULL; buf = buf->next) {
         ssp = (SelStructPtr) buf->region;
         if ( buf->itemtype == subtype && ssp->itemtype == type) {
                slp = (SeqLocPtr) ssp->region;
         }
         else slp = NULL;
         if ( buf->itemtype == itemsubtype && slp != NULL ) {
                ssp = (SelStructPtr) buf->region;
                if ( is_samessp(ssp, csp) ) {
                   locate_region(csp, SeqLocStart(slp), SeqLocStop(slp),  
                       SeqLocId(slp), SeqLocStrand(slp), FALSE);
                   locate = TRUE;
                   break;
                }
         }
  }
  if (!locate) return NULL;
  return csp;
}

/*********************************************************
***
***  locate_point 
***
**********************************************************/
extern Uint1 locate_point (PoinT pt, RecT rp, Uint4 *item_id, Uint2 *the_entity_id, Uint2 *item_type, Uint2 *item_subtype, Int4 *position, Int2 *line, EditAlignDataPtr adp)
{
  Uint4 itemid = 0;
  Uint2 seqEntityid = 0;
  Uint2 itemtype = 0;
  Uint2 itemsubtype = 0;
  Int4  pos = -1;
  Int4  x = -1;
  Int4  dx;
  Int2  ligne = -1;
  Uint1 what = BADLAND;

  *item_id = *the_entity_id = (Uint4) 0;
  *item_type= *item_subtype = (Uint2) 0;
  *position = (Int4) -1;
  *line = (Int2) -1;
  if ( pt.x < rp.left || pt.x > rp.right ) return 0;
  ligne = (pt.y - rp.top) / adp->lineheight; 
  if ( ligne >= 0 && ligne < adp->pnlLine ) {
	 itemid = adp->item_id [ligne];
	 seqEntityid = adp->seqEntity_id [ligne];
	 itemtype = adp->itemtype [ligne];
	 itemsubtype = adp->itemsubtype [ligne];
  }
  if (ligne > adp->pnlLine) {
         return HOLDDNLAND;
  }
  if (ligne < 0) {
         return HOLDUPLAND;
  }
  pos = (pt.x - rp.left - adp->margin.left);
/**
  WriteLog ("locate ii %d ei %d it %d  ist %d\n", adp->item_id [ligne], adp->seqEntity_id [ligne], adp->itemtype [ligne], adp->itemsubtype [ligne]);
**/
  if ( adp->seqEntity_id [ligne] > 0 && adp->item_id [ligne] > 0 ) 
  {
         if ( pos >= 0 ) 
         {
                x = (Int4) (pos + (adp->charw/6)) / adp->charw; 
                if (adp->columnpcell > 0)
                {
                        dx = (Int4) x / (Int4) adp->columnpcell;
                        dx = (Int4) (x -dx) / (Int4) adp->columnpcell;
                        x -= dx;
                }
                x += adp->hoffset - adp->bufferstart
                        + adp->alignline[ligne] * adp->visibleWidth;
                if (adp->colonne[x]==-1 && adp->seqnumber==1 && x>adp->length) {
                   x = adp->length;
                } else {
                   x = adp->colonne[x];
                }
         } 
         else if ( pos <  0 && pos > -adp->charw ) 
         {
                x = adp->hoffset - adp->bufferstart
                        + adp->alignline[ligne] * adp->visibleWidth;
                x = adp->colonne[x];
         }
         /* else WriteLog ("locate_point in margin\n"); */ 
         if ( x >= 0 && x <= adp->length ) 
               what = SEQLAND;
         else if ( x < 0 ) 
               what = SIDLAND;
  }  
  if ( what == SIDLAND ) 
  {
         *item_id   = itemid;
         *the_entity_id = (Uint2) seqEntityid;
         *item_type = (Uint2) itemtype;
         *item_subtype = (Uint2) itemsubtype;
         *position = (Int4) -1;
         *line = (Int2) ligne;
  } 
  else if ( what == SEQLAND ) 
  {
         *item_id   = itemid;
         *the_entity_id = (Uint2) seqEntityid;
         *item_type = (Uint2) itemtype;
         *item_subtype = (Uint2) itemsubtype;
         *position = (Int4) x;
         *line = (Int2) ligne;
  }
  return what;
}

extern void on_click (PaneL pnl, PoinT pt)
{
  EditAlignDataPtr adp;
  SeqAlignPtr      salp;
  SeqIdPtr         the_sip;
  SelStruct        sp;
  SelStructPtr     ssp = NULL;
  SelStructPtr     this_selection = NULL;
  RecT             rp;  
  Int4             position;
  Int4             pos_inseq;
  Int2             line;
  Uint2            itemsubtype;
  Uint1            what;

  if ( ( adp = GetAlignDataPanel (pnl)) == NULL ) return;
  if ( adp->seqnumber == 0 ) return;
  get_client_rect (pnl, &rp);
  sp.regiontype = 0;
  sp.region = NULL;
  what = locate_point (pt, rp, &sp.itemID, &sp.entityID, &sp.itemtype, 
          &itemsubtype, &position, &line, adp);
  if ( what >= BADLAND ) 
          return;
  if ( position > adp->length +1) {
          ErrPostEx (SEV_ERROR, 0, 0, "position> length");
          return;
  }
  adp->select_block = NULL;
  adp->mouse_mode = MOUSE_CLICK;
/***
  WriteLog("click what %d  %d %d %d %d pos %d lin %d  \n", (int) what, 
    (int) sp.itemID, (int) sp.entityID, (int)sp.itemtype, (int)itemsubtype,
    (int) position, (int)line);
****/
  salp = (SeqAlignPtr) adp->sap_align->data;
  switch (what) 
  {
    case SIDLAND: 
          if ( itemsubtype == FEATDEF_BAD ) 
          {
             ssp = make_bspssp (&sp, adp->anp_list, 0, -2);
             adp->clickwhat = CLICKID_BIOSEQ;
             if (dblClick) {
                if (adp->edit_mode==SEQ_VIEW || adp->edit_mode==SEQ_EDIT || adp->edit_mode==ALIGN_EDIT) {
                   WatchCursor ();
                   Update ();
                   GatherProcLaunch (OMPROC_EDIT, FALSE, sp.entityID, sp.itemID, OBJ_BIOSEQ, 0, 0, OBJ_BIOSEQ, 0);
                   ArrowCursor ();
                   Update ();
                   ssp= NULL;
                }
             }
          }
          else if (sp.itemtype == OBJ_SEQFEAT) 
          {
             if (dblClick) 
             {
                WatchCursor ();
                Update ();
                GatherProcLaunch (OMPROC_EDIT, FALSE, sp.entityID, sp.itemID, sp.itemtype, 0, 0, sp.itemtype, 0);
                ArrowCursor ();
                Update ();
                ssp = NULL;
             }
             else {
                the_sip = get_featId_fromid(&sp, adp->seqfeat, itemsubtype);
                pos_inseq = AlignCoordToSeqCoord (position, the_sip, salp, adp->sqloc_list, 0);
                ssp = make_featssp (&sp, adp->seqfeat, itemsubtype, pos_inseq, -1, -1);
                adp->clickwhat = CLICKID_SEQFEAT;
             }
          }
          else if( itemsubtype >=EDITDEF_RF1 && itemsubtype<=EDITDEF_RF6 )
          {
             ssp = make_bufssp (&sp, adp->buffer, itemsubtype);
             adp->clickwhat = CLICKID_FEAT;
          }
          if (ssp != NULL) {
             click_id (pnl, ssp, adp);
          }
          break;
    case SEQLAND: 
          if ( position >=0 ) {
             adp->clickwhat = CLICKSEQ_BIOSEQ;
             if(adp->edit_mode==SEQ_VIEW || adp->edit_mode==SEQ_EDIT || adp->edit_mode == ALIGN_EDIT) 
             {
                if ( itemsubtype == FEATDEF_BAD ) 
                {
                   if (!shftKey) {
		      ssp=make_bspssp (&sp, adp->anp_list, position, position);
                      click_caret (&rp, position, ssp, &(adp->caret), adp, line);
                      if (!ctrlKey)
                         ObjMgrDeSelectAll ();
                   }
		   else 
                   {
                      this_selection = get_closest_selection(&(adp->caret),salp, &sp, adp->sqloc_list);
                      if (this_selection!=NULL) {
			 extend_selectregion (adp, this_selection, position, adp->caret_orig, FALSE, pnl);
                      }
                   }
                   if (!adp->display_panel)
                      to_update_prompt(pnl,&(adp->caret), salp, adp->sqloc_list, FALSE, adp->printid);
                   adp->edit_pos = position;
                   adp->caret_line = line;
                }
                else if (sp.itemtype == OBJ_SEQFEAT)
                {
                   if (dblClick) 
                   {
                      WatchCursor ();
                      Update ();
                      GatherProcLaunch (OMPROC_EDIT, FALSE, sp.entityID, sp.itemID, sp.itemtype, 0, 0, sp.itemtype, 0);
                      ArrowCursor ();
                      Update ();
                   }
                   else {
                      the_sip = get_featId_fromid(&sp, adp->seqfeat, itemsubtype);
                      pos_inseq = AlignCoordToSeqCoord (position, the_sip, salp, adp->sqloc_list, 0);
                      ssp = make_featssp (&sp, adp->seqfeat, itemsubtype, pos_inseq, -1, -1);
                      if (ssp != NULL) 
                      {
                         click_feat (pnl, position, ssp, adp, line);
                         adp->feat_pos = position;
                         adp->clickwhat = CLICKSEQ_SEQFEAT;
                      }
                   }
                }
             }
          } 
          break;
    default:
          adp->clickwhat = NOCLICK;
          break;
  }
  Update ();
}

/**************************************
*** 
***  Drag Procedure
***        2: drag on sequence
***
***************************************/
/*if CtrlKey is pushed and the user clicks within a selected region,
don't move the caret*/
static Int4 getcurrentcaretpos (Int4 position, Uint2 eID, Uint2 iID, Uint2 itype, Boolean direction)
{
   SelStructPtr sspsel=NULL;

 if (ctrlKey){
	sspsel=ObjMgrGetSelected();	
    	if (sspsel != NULL) 
    	{
       		for (; sspsel != NULL; sspsel = sspsel->next)
       		{
       		   if ( checkssp_for_editor (sspsel) && is_sameId (sspsel->entityID, sspsel->itemID, sspsel->itemtype, 255, eID, iID, itype, 255) )
		 {
            		  if ( (position>= SeqLocStart((SeqLocPtr)sspsel->region)) && 
				(position<= SeqLocStop((SeqLocPtr)sspsel->region))){
				if (direction)
					return SeqLocStart((SeqLocPtr)sspsel->region);
				else 
					return SeqLocStop((SeqLocPtr)sspsel->region);
			   }
        	   }
       		}
    	}
  }
  return -1;
}

extern void on_drag (PaneL pnl, PoinT pt)
{
  EditAlignDataPtr adp;
  SeqAlignPtr      salp;
  SeqIdPtr         the_sip;
  SelStruct        sp;
  SelStructPtr     this_selection = NULL;
  RecT             rp;  
  Int4             position;
  Int2             line;
  Uint2            itemsubtype;
  Int4             from, to,
                   caret_origseq, new_caret_orig;
  SeqLocPtr        slp;
  Uint1            direction = Seq_strand_plus;
  
  
  if ( ( adp = GetAlignDataPanel (pnl)) == NULL ) return;
  if ( adp->seqnumber == 0 ) return;
  get_client_rect (pnl, &rp);
  sp.regiontype = 0;
  sp.region = NULL;
  locate_point (pt, rp, &sp.itemID, &sp.entityID, &sp.itemtype, &itemsubtype, &position, &line, adp);
  if ( position > adp->length) {
          ErrPostEx (SEV_ERROR, 0, 0, "position> length");
          return;
  }
  salp = (SeqAlignPtr) adp->sap_align->data;
  switch ( adp->clickwhat ) 
  {
    case CLICKID_BIOSEQ: 
          break;
    case CLICKID_SEQFEAT: 
          break;
    case CLICKID_FEAT: 
          break;
    case CLICKSEQ_BIOSEQ: 
          if (position>=0) {
           if (position != adp->edit_pos || adp->mouse_mode == MOUSE_DRAG)
           {
            this_selection = &sp;
            if (is_sameId (this_selection->entityID, this_selection->itemID, 
               this_selection->itemtype, 255, adp->caret.entityID, adp->caret.itemID, 
               adp->caret.itemtype, 255))
            {  
                if (position < adp->edit_pos) {
                   direction = Seq_strand_minus;
                }
                the_sip = SeqLocId((SeqLocPtr)adp->caret.region);
                if (the_sip) {
                 caret_origseq = AlignCoordToSeqCoord (adp->caret_orig, the_sip, salp, adp->sqloc_list, 0);
                 if (adp->mouse_mode == MOUSE_CLICK)
                 {
                   new_caret_orig = getcurrentcaretpos (caret_origseq, adp->caret.entityID, adp->caret.itemID, adp->caret.itemtype, (Boolean)(direction==Seq_strand_plus));
                   if (new_caret_orig<0) 
                     new_caret_orig=caret_origseq;
                   else if (caret_origseq != new_caret_orig) {
                      adp->caret_orig = SeqCoordToAlignCoord (new_caret_orig, the_sip, salp, 0, 0);
                   }
                 }
                 else {
                   new_caret_orig=caret_origseq;
                 }
                 from = new_caret_orig;
                 to = AlignCoordToSeqCoord (position, the_sip, salp, adp->sqloc_list, 0);
                 if (from > to) {
                   from = to;
                   to = new_caret_orig;
                 }
                 slp = SeqLocIntNew (from, to, direction, the_sip);
                 ObjMgrAlsoSelect (this_selection->entityID, this_selection->itemID, 
                   this_selection->itemtype, OM_REGION_SEQLOC, slp);
                 adp->edit_pos = position;
                 adp->mouse_mode = MOUSE_DRAG;
                }
             }
           }
           else if (adp->mouse_mode == MOUSE_CLICK)
              adp->mouse_mode = MOUSE_DRAG;
          }
          break;
          
    case CLICKSEQ_SEQFEAT: 
/********************
          if ( select != NULL ) 
          {
             sspei = (Uint2) select->entityID;
             sspii = (Uint2) select->itemID;
             sspit = (Uint2) select->itemtype;
             sspist= (Uint2) select->regiontype;
             if (position >= 0 && position < adp->length && sspit == OBJ_SEQFEAT)
             {
                ssp = is_selectedbyID (sspei, sspii, sspit); 
                if ( ssp == NULL )  
                {
                   click_feat (pnl, position, adp->select, adp, line);
                   ssp = is_selectedbyID (sspei, sspii, sspit); 
                }
                if ( ssp != NULL && position != adp->feat_pos )  
                {
                   prec = NULL; 
                   pos_inseq = AlignCoordToSeqCoord (adp->feat_pos, SeqLocId((SeqLocPtr)ssp->region), salp, adp->sqloc_list, 0);
                   feat = get_feat_fromid (adp->seqfeat, sspist, sspei, sspii, pos_inseq, &prec);
                   if ( feat != NULL ) 
                   {
                       drag_feat (pnl, position, feat, adp, line, sspist, prec, feat->next, ssp);
                   }
                   adp->feat_pos = position;
                }
             }
          }
**********************/
          break;
          
    case NOCLICK: 
          break;
    default:
          break;
  }
  Update ();
}

/**************************************
*** 
***  Hold Procedure
***        4: hold on sequence
***        5: hold on sequence
***
***************************************/
extern void on_hold (PaneL pnl, PoinT pt)
{
  EditAlignDataPtr adp;
  SeqAlignPtr      salp;
  RecT             rp;  
  SelStruct        sp;
  SelStructPtr     this_selection = NULL;
  Int4             position;
  Int2             line;
  Uint2            itemsubtype;
  Uint1            what;

  SelEdStructPtr   feat, prec;

  if ( ( adp = GetAlignDataPanel (pnl)) == NULL ) return;
  if ( adp->seqnumber == 0 ) return;
  get_client_rect (pnl, &rp);
  sp.regiontype = 0;
  sp.region = NULL;
  what = locate_point (pt, rp, &sp.itemID, &sp.entityID, &sp.itemtype, &itemsubtype, &position, &line, adp);
  ResetClip ();
  salp = (SeqAlignPtr) adp->sap_align->data;
  if ( what == HOLDDNLAND ) {
          if (adp->voffset >= adp->nlines-1) return;
          adp->voffset = movedown_scrollbar (pnl, adp->voffset, (Int2)1);
          Update ();
          if ( adp->clickwhat == CLICKSEQ_BIOSEQ) 
          {
             position = adp->edit_pos + adp->visibleWidth;
             if (position > adp->length) position = adp->length;
             adp->caret_line++;
             dragcaret_selectrgn (pnl, position, get_closest_selection(&(adp->caret), salp, NULL, adp->sqloc_list), &(adp->caret), adp, adp->caret_line);
             setposition_tossp (&(adp->caret), position, position);
             adp->edit_pos = position;
          }
          else if (adp->clickwhat== CLICKSEQ_SEQFEAT || adp->clickwhat == CLICKSEQ_FEAT)
          {
             prec = NULL; 
             this_selection = get_closest_selection(&(adp->caret), salp, NULL, adp->sqloc_list);
             feat = get_feat_fromid (adp->feat, this_selection->regiontype, this_selection->entityID, this_selection->itemID, adp->feat_pos, &prec);
             if ( feat != NULL ) 
             {
                   adp->feat_line++;
                   position = adp->feat_pos + adp->visibleWidth;
                   drag_feat (pnl, position, feat, adp, adp->feat_line, this_selection->regiontype, prec, feat->next, NULL);
             }
             adp->feat_pos = position;
          }
  }
  else if ( what == HOLDUPLAND ) 
  {
          if (adp->voffset == 0 || adp->hoffset == 0) return;
          adp->voffset = moveup_scrollbar (pnl, adp->voffset, (Int2)1);
          Update ();
          if ( adp->clickwhat == CLICKSEQ_BIOSEQ) 
          {
             position = adp->edit_pos - adp->visibleWidth;
             if (position < 0) position = 0;
             adp->caret_line--;
             dragcaret_selectrgn (pnl,position,get_closest_selection(&(adp->caret), salp, NULL, adp->sqloc_list),&(adp->caret), adp, adp->caret_line);
             setposition_tossp (&(adp->caret), position, position);
             adp->edit_pos = position;
          }
          else if (adp->clickwhat == CLICKSEQ_SEQFEAT || adp->clickwhat == CLICKSEQ_FEAT)
          {
             this_selection = get_closest_selection(&(adp->caret), salp, NULL, adp->sqloc_list);
             feat = get_feat_fromid (adp->feat, this_selection->regiontype, this_selection->entityID, this_selection->itemID, adp->feat_pos, &prec);
             if ( feat != NULL ) 
             {
                   adp->feat_line--;
                   position = adp->feat_pos - adp->visibleWidth;
                   drag_feat (pnl, position, feat, adp, adp->feat_line, this_selection->regiontype, prec, feat->next, NULL);
             }
             adp->feat_pos = position;
          }
  }
  Update ();
}

/**************************************
*** 
***  Release Procedure
***
***************************************/
extern void on_release (PaneL pnl, PoinT pt)
{
  EditAlignDataPtr adp;
  SelStruct        sp;
  SelStructPtr     ssp = NULL;
  RecT             rp;  
  Int2             line;
  Int4             position;
  Uint2            itemsubtype;

  if ( ( adp = GetAlignDataPanel (pnl)) == NULL ) return;
  if ( adp->seqnumber == 0 ) return;
  get_client_rect (pnl, &rp);
  sp.regiontype = 0;
  sp.region = NULL;
  locate_point (pt, rp, &sp.itemID, &sp.entityID, &sp.itemtype, &itemsubtype, &position, &line, adp);
  if ( position > adp->length +1) {
          ErrPostEx (SEV_ERROR, 0, 0, "position> length");
          return;
  }
  switch ( adp->clickwhat ) 
  {
    case CLICKID_BIOSEQ: 
          break;      
    case CLICKSEQ_BIOSEQ: 
          if ( position >= 0 && adp->mouse_mode == MOUSE_DRAG) {
                ssp = make_bspssp (&sp, adp->anp_list, position, position);
                if (ssp != NULL) {
                   release_caret (pnl, position, ssp, &adp->caret, adp, line);
                   adp->edit_pos = position;
                }
          }
          break;
    case NOCLICK: 
          ObjMgrDeSelectAll ();
          break;
    default:
          break;
  }
  adp->clickwhat = NOCLICK;
  adp->mouse_mode = MOUSE_RELEASE;
  if (!adp->display_panel) {
     update_edititem (pnl);
     update_savefeatitem (pnl);  
     update_translateitem (pnl, adp->seqfeat, adp->feat);  
     update_codonstartbt (pnl, adp->seqfeat, adp->feat);
  }
  Update ();
}

#define FINDNEXT 14
#define FINDPREV 16
#define SNLM_DEL   8
#define NLM_RETURN 13

/*-----------------------*/
static Boolean gap_insert (EditAlignDataPtr adp)
{
  Int4    position;
  Int4    oldbufferlength;
  Int2    newvisibleLine,
          oldvisibleLine;
  Boolean insert = FALSE;

  if ( (adp->bufferlength+1) >= adp->minbufferlength + adp->editbuffer) {
     ErrPostEx (SEV_ERROR, 0, 0, "Save alignment before");
     return FALSE;
  }
  position = SeqLocStart((SeqLocPtr)adp->caret.region);
  insert = AlignDataGapAddProc ('-', position, adp->linebuff, adp->bufferlength, adp->caret.entityID, adp->caret.itemID); 
  if (insert) 
  {
     oldbufferlength = adp->bufferlength;
     adp->bufferlength = CleanBufferProc(adp->linebuff, adp->seqnumber, adp->bufferlength +1); 
     adp->edit_pos++;
     adp->colonne [adp->bufferlength] =adp->colonne [adp->bufferlength -1] +1;
     oldvisibleLine =1 +MIN ((Int2) ((oldbufferlength-1) / adp->visibleWidth), (Int2) ((adp->pnlLine -3)/ (adp->seqnumber +1 +2)));
     newvisibleLine =1 +MIN ((Int2) ((adp->bufferlength-1) / adp->visibleWidth), (Int2) ((adp->pnlLine -3)/ (adp->seqnumber +1 +2)));
     if ( newvisibleLine < 1)  
     newvisibleLine = 1;
     if (newvisibleLine != oldvisibleLine)
     {
        adp->visibleLength =adp->visibleWidth * newvisibleLine;
     }
     if(oldbufferlength<adp->bufferlength && adp->length<adp->minbufferlength)
     {
        adp->length++;
     }
     return TRUE;
  }
  return FALSE;
}

extern void on_key (SlatE s, Char ch)
{
  PaneL            pnl;
  EditAlignDataPtr adp;
  SeqAlignPtr      salp;
  RecT             rp;
  SeqLocPtr        slp;
  SelStructPtr     ssp;
  Int4             caretstart;
  Int4             new_caretstart;
  Int4             k;
  Int4             line, oldline, 
                   column, oldcolumn;
  Int4             j;
  Boolean          is_caretvisible;
  Boolean          minusstrand;
  CharPtr          str;
  Int4             width;

  if ( (int) ch == 0 ) return;

  pnl = (PaneL) s;
  Select (pnl);
  if ( (adp = GetAlignDataPanel (pnl)) == NULL ) return;
  if ( adp->seqnumber == 0 ) return;
  get_client_rect (pnl, &rp);
  if ( adp->caret.regiontype == 0 || adp->caret.region == NULL ) { 
           return;
  }
  slp = (SeqLocPtr) adp->caret.region;
  caretstart = SeqLocStart(slp);
  if (caretstart < 0) {
           ErrPostEx (SEV_ERROR, 0, 0, "CARET position negative");
           return;
  }
  is_caretvisible = (Boolean)(caretstart >= adp->hoffset 
                    && caretstart <= adp->hoffset + adp->visibleLength);
  if (!is_caretvisible) 
  {
     adp->voffset=(Int2)(caretstart /(adp->visibleWidth));
     ResetClip ();
     SeqEdSetValueScrollBar (pnl, (Int2) adp->voffset);
     data_collect_arrange (adp, TRUE);
     Select (pnl);
     inval_panel (pnl, -1, -1);
     Update ();
  }
  width = adp->visibleWidth;
  if (adp->columnpcell > 0) 
     width += (Int2) adp->visibleWidth / (Int2) adp->columnpcell;
  salp = (SeqAlignPtr) adp->sap_align->data;
  if (!ctrlKey && adp->edit_mode == SEQ_EDIT ) 
  {
     if ( (str = char_to_insert (&ch, 1, adp->mol_type)) != NULL) 
     {
        if (adp->input_format == OBJ_BIOSEQ) 
        {
           ssp = ObjMgrGetSelected ();
           if ( checkssp_for_editor(ssp) && ssp->itemtype == OBJ_BIOSEQ) 
           {
              do_cut (pnl, adp, ssp, FALSE);
           }
           else ssp = NULL;
           if (insertchar_atcaret (str, adp)) 
           {
              adp->dirty = TRUE;
              ObjMgrSetDirtyFlag (adp->caret.entityID, TRUE);

              ObjMgrSendMsg (OM_MSG_UPDATE, adp->caret.entityID, adp->caret.itemID, adp->caret.itemtype);

              if (adp->edit_item.entityID == 0) 
              {
                 adp->edit_item.entityID = adp->caret.entityID;
                 adp->edit_item.itemID = adp->caret.itemID;
                 adp->edit_item.itemtype = adp->caret.itemtype;
              }
              else if (adp->edit_item.entityID != adp->caret.entityID && adp->edit_item.itemID != adp->caret.itemID) {
                 ErrPostEx (SEV_ERROR, 0, 0, "Warning in SetWindowTimer");
                 adp->edit_item.entityID = 0; 
                 adp->edit_item.itemID = 0; 
                 adp->edit_item.itemtype = 0; 
              }
              if (! typing_timerInUse) 
              {
                 typing_timerInUse = TRUE;
                 timerCount = 0;
              }
           }
           MemFree (str);
        }
        else if (adp->input_format == OBJ_SEQALIGN)
        {
           Beep ();
        }
        return;
     }
  }
  if (!ctrlKey &&  adp->edit_mode == ALIGN_EDIT ) 
  {
     if (ch=='-') {
        if (adp->length < adp->minbufferlength) {
           if (gap_insert (adp)) {
              setposition_tossp (&adp->caret, caretstart+1, caretstart+1);
              adp->dirty = TRUE;
              inval_panel (pnl, -1, -1);
           }
        }
        return;
     }
  } 

  switch ((int) TO_UPPER(ch)) 
  {
     case SNLM_DEL:
     case NLM_DEL:
          if ( adp->edit_mode == SEQ_EDIT ) 
          {
               k = checkOMss_for_itemtype (OBJ_BIOSEQ);
               if ( k > 0 ) {
                  ssp = getOMselect_for_itemtype (OBJ_BIOSEQ);
               }
               else {
                  ssp = (SelStructPtr) MemNew (sizeof (SelStruct));
                  ssp->entityID = adp->caret.entityID;
                  ssp->itemID = adp->caret.itemID;
                  ssp->itemtype = adp->caret.itemtype;
                  locate_region (ssp, caretstart - 1, caretstart - 1, SeqLocId(slp), Seq_strand_plus, FALSE);
                  ObjMgrSelect (ssp->entityID, ssp->itemID, ssp->itemtype,
ssp->regiontype, ssp->region);
               }
               if ( ssp != NULL ) {
                  if (adp->input_format == OBJ_BIOSEQ)
                  {
                     if (do_cut (pnl, adp, ssp, FALSE)) 
                     { 
                        if (adp->edit_item.entityID == 0) 
                        {
                           adp->edit_item.entityID = adp->caret.entityID;
                           adp->edit_item.itemID = adp->caret.itemID;
                           adp->edit_item.itemtype = adp->caret.itemtype;
                        }  
                        else if (adp->edit_item.entityID != adp->caret.entityID && adp->edit_item.itemID != adp->caret.itemID) {
                           ErrPostEx (SEV_ERROR, 0, 0, "Warning in SetWindowTimer");
                           adp->edit_item.entityID = 0;
                           adp->edit_item.itemID = 0;
                           adp->edit_item.itemtype = 0;
                        }  
                        if (! deleting_timerInUse) {
                           deleting_timerInUse = TRUE;
                           timerCount = 0;
                        }  
                     }
                  }
                  else if (adp->input_format == OBJ_SEQALIGN)
                  {
                     Beep ();
                  }
               }
          }
          break;

     case NLM_LEFT: 
          if (ctrlKey) {
             adp->cur_pat = ShowPrecPattern(adp->match_pat, adp->cur_pat, adp->edit_pos);
          }
          else if ( SeqLocStart(slp) > 0 ) 
          {
               minusstrand=(Boolean)(SeqLocStrand((SeqLocPtr)adp->caret.region) == Seq_strand_minus);
               if (minusstrand)
                  new_caretstart = caretstart+1;
               else 
                  new_caretstart = caretstart-1;
               SeqPosToLineColumn (adp->caret.itemID, adp->caret.entityID,  adp->caret.itemtype, caretstart, &oldline, &oldcolumn, adp->hoffset, adp);
               SeqPosToLineColumn (adp->caret.itemID, adp->caret.entityID, adp->caret.itemtype, new_caretstart, &line, &column, adp->hoffset, adp);
               setposition_tossp (&adp->caret, new_caretstart, new_caretstart);
               if ( oldline>=0 && ( abs(oldcolumn - column) >2 || column <0)) 
                   inval_rect((Int2)(rp.left+adp->margin.left+(oldcolumn-2)*adp->charw),
                       (Int2)(rp.top + adp->lineheight * oldline),
                       (Int2)(rp.left+ adp->margin.left + (oldcolumn+2)*adp->charw),
                       (Int2)(rp.top + adp->lineheight * (oldline+1))); 
               if ( line>=0 && column>=0 ) {
                   inval_rect((Int2)(rp.left +adp->margin.left +(column-2)*adp->charw),
                       (Int2)(rp.top + adp->lineheight * line),
                       (Int2)(rp.left + adp->margin.left + (column+2)*adp->charw),
                       (Int2)(rp.top + adp->lineheight * (line+1))); 
               }
               caretstart = SeqLocStart(slp);
               if ( ! shftKey ) {
                   slp = (SeqLocPtr) adp->caret.region;
                   if (!adp->display_panel) {
                      to_update_prompt (pnl, &(adp->caret), (SeqAlignPtr) adp->sap_align->data, adp->sqloc_list, FALSE, adp->printid);
                   }
                   adp->caret_orig = caretstart;
                   if (minusstrand) 
                      adp->edit_pos--;
                   else 
                      adp->edit_pos++;
                   ObjMgrDeSelectAll ();
               } 
               else {
                   dragcaret_selectrgn (pnl, caretstart, get_closest_selection(&(adp->caret), salp, NULL, adp->sqloc_list), &(adp->caret), adp, line);
                   adp->edit_pos --;
               }
          }
          break;

     case NLM_RIGHT: 
          if (ctrlKey) {
             adp->cur_pat = ShowNextPattern(adp->match_pat, adp->cur_pat, adp->edit_pos);
          }
          else if ( SeqLocStart(slp) < adp->length ) 
          {
               minusstrand=(Boolean)(SeqLocStrand((SeqLocPtr)adp->caret.region) == Seq_strand_minus);
               if (minusstrand)
                  new_caretstart = caretstart-1;
               else 
                  new_caretstart = caretstart+1;
               SeqPosToLineColumn (adp->caret.itemID, adp->caret.entityID,  adp->caret.itemtype, caretstart, &oldline, &oldcolumn, adp->hoffset, adp);
               SeqPosToLineColumn (adp->caret.itemID, adp->caret.entityID, adp->caret.itemtype, new_caretstart, &line, &column, adp->hoffset, adp);
               setposition_tossp (&adp->caret, new_caretstart, new_caretstart);
               if ( oldline>=0 && ( abs(oldcolumn - column) >2 || column <0)) 
               {      
                  inval_rect ((Int2) (rp.left +adp->margin.left +(oldcolumn-2) *adp->charw), (Int2)(rp.top + adp->lineheight * oldline), (Int2)(rp.left+ adp->margin.left + (oldcolumn+2)*adp->charw), (Int2)(rp.top + adp->lineheight * (oldline+1))); 
               }
               if ( line>=0 && column>=0 ) {
                  inval_rect ((Int2) (rp.left +adp->margin.left +(column-3) *adp->charw),  (Int2)(rp.top + adp->lineheight * line), (Int2)(rp.left + adp->margin.left + (column+2)*adp->charw), (Int2)(rp.top + adp->lineheight * (line+1))); 
               }
               caretstart = SeqLocStart(slp);
               if ( ! shftKey ) {
                   slp = (SeqLocPtr) adp->caret.region;
                   if (!adp->display_panel)
                      to_update_prompt(pnl, &(adp->caret), (SeqAlignPtr) adp->sap_align->data, adp->sqloc_list, FALSE, adp->printid);
                   adp->caret_orig = caretstart;
                   if (minusstrand) 
                      adp->edit_pos += 1;
                   else 
                      adp->edit_pos -= 1;
                   ObjMgrDeSelectAll ();
               } 
               else {
                   dragcaret_selectrgn (pnl, caretstart, get_closest_selection(&(adp->caret), salp, NULL, adp->sqloc_list), &(adp->caret), adp, line);
                   adp->edit_pos ++;
               }
          }
          break;

     case NLM_UP: 
          SeqPosToLineColumn (adp->caret.itemID, adp->caret.entityID, 
                       adp->caret.itemtype, caretstart, &line, &column,
                       adp->hoffset, adp);
          adp->caret_line  = line;
          j = nextlineup (adp->caret_line, adp->itemtype, adp->itemsubtype); 
          if ( j < 0 ) 
          {
             Int4      pos;
             if (adp->voffset != 0 && adp->hoffset != 0) 
             {
                ResetClip ();
                adp->voffset = moveup_scrollbar (pnl, adp->voffset, 1);
                Update ();
                if (is_samess_ses (&(adp->caret), (SelEdStructPtr) adp->firstssp->region) )
                {
                   pos = SeqLocStart(slp) - adp->visibleWidth;
                   if (pos < 0) pos = 0;
                   adp->caret_line = 0;
                   if ( ! shftKey ) {
                      click_caret (&rp, pos, &(adp->caret), &(adp->caret), adp, adp->caret_line);
                      if (!adp->display_panel)
                         to_update_prompt (pnl, &(adp->caret), (SeqAlignPtr) adp->sap_align->data, adp->sqloc_list, FALSE, adp->printid);
                      ObjMgrDeSelectAll ();
                   } 
                   else {
                      dragcaret_selectrgn (pnl, pos, get_closest_selection(&(adp->caret), salp, NULL, adp->sqloc_list), &(adp->caret), adp, adp->caret_line);
                      setposition_tossp (&(adp->caret), pos, pos);
                      adp->edit_pos = pos;
                   }
                }
             }
          }
          else if ( j != adp->caret_line ) 
          {
               SeqIdPtr sip;
               Uint2 eid, it;
               Uint4 iid;
               eid = adp->caret.entityID;
               iid = adp->caret.itemID;
               it = adp->caret.itemtype;
               SeqPosToLineColumn (iid, eid, it, caretstart, &oldline, 
                                   &oldcolumn, adp->hoffset, adp);
               if (adp->alignline[adp->caret_line] != adp->alignline[j]) 
               {
                   line=abs(adp->alignline[adp->caret_line] -adp->alignline[j]);
                   caretstart = caretstart -(line *adp->visibleWidth);
                   if (caretstart < 0) caretstart = 0;
                   setposition_tossp (&(adp->caret), caretstart, caretstart);
               }
               else
                   caretstart = SeqLocStart(slp);
               if (adp->seqEntity_id[adp->caret_line] != adp->seqEntity_id[j]
                || adp->item_id[adp->caret_line] != adp->item_id[j]) 
               {
                   eid = adp->seqEntity_id[j];
                   iid = adp->item_id[j];
                   it = adp->itemtype[j];
                   adp->caret.entityID = eid; 
                   adp->caret.itemID = iid; 
                   adp->caret.itemtype = it; 
                   sip=(SeqIdPtr)SeqIdFromAlignNode(adp->anp_list,eid, iid, it);
                   replace_region (&(adp->caret), eid, iid, it, caretstart, caretstart, sip, Seq_strand_plus, FALSE);
               }
               adp->caret_line = j;
               SeqPosToLineColumn (iid, eid, it, caretstart, &line, &column,
                       adp->hoffset, adp);
               inval_rect((Int2)(rp.left+adp->margin.left+(column-2)* adp->charw), (Int2)(rp.top + adp->lineheight * oldline), (Int2)(rp.left+ adp->margin.left + (column+2)*adp->charw), (Int2)(rp.top + adp->lineheight * (oldline+1))); 
               inval_rect((Int2)(rp.left+adp->margin.left+(column-2)* adp->charw), (Int2)(rp.top + adp->lineheight * line), (Int2)(rp.left + adp->margin.left + (column+2)*adp->charw), (Int2)(rp.top + adp->lineheight * (line+1))); 
               if ( ! shftKey ) 
               {
                   slp = (SeqLocPtr) adp->caret.region;
                   if (!adp->display_panel)
                      to_update_prompt (pnl, &(adp->caret), (SeqAlignPtr) adp->sap_align->data, adp->sqloc_list, FALSE, adp->printid);
                   adp->caret_orig = caretstart;
                   adp->edit_pos = caretstart;
                   ObjMgrDeSelectAll ();
               } 
               else {
                   dragcaret_selectrgn (pnl, caretstart, get_closest_selection(&(adp->caret), salp, NULL, adp->sqloc_list), &(adp->caret), adp, line);
                   adp->edit_pos -= adp->visibleWidth;
               }
          } 
          break;

     case NLM_DOWN: 
          SeqPosToLineColumn (adp->caret.itemID, adp->caret.entityID, 
                       adp->caret.itemtype, caretstart, &line, &column,
                       adp->hoffset, adp);
          adp->caret_line  = line;
          j = nextlinedown (adp->caret_line, adp->pnlLine, adp->itemtype, adp->itemsubtype);
          if ( j < 0 ) 
          {
             Int4      pos;
             if (adp->voffset < adp->nlines-1) 
             {
                ResetClip ();
                adp->voffset = movedown_scrollbar (pnl, adp->voffset, (Int2)1);
                Update ();
                if (is_samess_ses (&(adp->caret), adp->lastses) ) 
                {
                   pos = SeqLocStart(slp) + adp->visibleWidth;
                   if (pos > adp->length) pos = adp->length;
                   adp->caret_line = adp->pnlLine -1;
                   if ( ! shftKey ) {
                      click_caret (&rp, pos, &(adp->caret), &(adp->caret), adp, adp->caret_line);
                      if (!adp->display_panel)
                         to_update_prompt (pnl, &(adp->caret), (SeqAlignPtr) adp->sap_align->data, adp->sqloc_list, FALSE, adp->printid);
                      ObjMgrDeSelectAll ();
                   }
                   else {
                      dragcaret_selectrgn (pnl, pos, get_closest_selection(&(adp->caret), salp, NULL, adp->sqloc_list), &(adp->caret), adp, adp->caret_line);
                      setposition_tossp (&(adp->caret), pos, pos);
                      adp->edit_pos = pos;
                   }
                }
             }
          }
          else if (j != adp->caret_line ) 
          {
               SeqIdPtr sip;
               Uint2 eid, it;
               Uint4 iid;
               eid = adp->caret.entityID;
               iid = adp->caret.itemID;
               it = adp->caret.itemtype;
               SeqPosToLineColumn (adp->caret.itemID, adp->caret.entityID, 
                       adp->caret.itemtype, caretstart, &oldline, 
                       &oldcolumn, adp->hoffset, adp);
               if (adp->alignline[adp->caret_line] != adp->alignline[j]) 
               {
                   line=abs(adp->alignline[adp->caret_line] -adp->alignline[j]);
                   caretstart = caretstart + (line * adp->visibleWidth);
                   if (caretstart > adp->length) 
                      caretstart = adp->length;
                   setposition_tossp (&(adp->caret), caretstart, caretstart);
               }
               else
                   caretstart = SeqLocStart(slp);
               if (adp->seqEntity_id[adp->caret_line] != adp->seqEntity_id[j]
                || adp->item_id[adp->caret_line] != adp->item_id[j]) 
               {
                   eid = adp->seqEntity_id[j];
                   iid = adp->item_id[j];
                   it = adp->itemtype[j];
                   adp->caret.entityID = eid; 
                   adp->caret.itemID = iid; 
                   adp->caret.itemtype = it; 
                   sip=(SeqIdPtr)SeqIdFromAlignNode(adp->anp_list,eid, iid, it);
                   replace_region (&(adp->caret), eid, iid, it, caretstart, caretstart, sip, Seq_strand_plus, FALSE);
               }
               adp->caret_line = j;

               SeqPosToLineColumn (iid, eid, it, caretstart, &line, &column, 
                       adp->hoffset, adp);
               inval_rect ((Int2) (rp.left+adp->margin.left +(column-2) * adp->charw), (Int2)(rp.top + adp->lineheight * oldline), (Int2)(rp.left+ adp->margin.left + (column+2)*adp->charw), (Int2)(rp.top + adp->lineheight * (oldline+1))); 
               inval_rect((Int2)(rp.left+adp->margin.left+(column-2) * adp->charw), (Int2)(rp.top + adp->lineheight * line), (Int2)(rp.left + adp->margin.left + (column+2)*adp->charw), (Int2)(rp.top + adp->lineheight * (line+1))); 
               if ( ! shftKey ) {
                   slp = (SeqLocPtr) adp->caret.region;
                   if (!adp->display_panel)
                      to_update_prompt (pnl, &(adp->caret), (SeqAlignPtr) adp->sap_align->data, adp->sqloc_list, FALSE, adp->printid);
                   adp->caret_orig = caretstart;
                   adp->edit_pos = caretstart;
                   ObjMgrDeSelectAll ();
               } 
               else {
                   dragcaret_selectrgn(pnl, caretstart, get_closest_selection(&(adp->caret), salp, NULL, adp->sqloc_list), &(adp->caret),  adp, line);
                   adp->edit_pos += adp->visibleWidth;
               }
          } 
          break;

     case FINDNEXT:
          adp->cur_pat = ShowNextPattern(adp->match_pat, adp->cur_pat, adp->edit_pos);
          break;

     case FINDPREV:
          adp->cur_pat = ShowPrecPattern(adp->match_pat, adp->cur_pat, adp->edit_pos);
          break;

     case NLM_RETURN:
{
  WindoW w;
  SeqEditViewFormPtr wdp;

          w = getwindow_frompanel (pnl);
          wdp = (SeqEditViewFormPtr) GetObjectExtra (w);
          GoToFunc (wdp);
}
          break;

     case NLM_ESC: 
          break;

     default:
          Beep ();
          break;
  }
  if (!adp->display_panel)
     update_edititem (pnl);
  Update ();
}


const char *nucleotide_alphabet = "ABCDGHKMRSTUVWYabcdghkmrstuvwy";
const char *protein_alphabet = "ABCDEFGHIKLMPQRSTUVWXYZabcdefghiklmpqrstuvwxyz";

typedef struct alnsettingsdlg {
  DIALOG_MESSAGE_BLOCK
  TexT missing;
  TexT match;
  TexT beginning_gap;
  TexT middle_gap;
  TexT end_gap;

  PopuP sequence_type;
} AlnSettingsDlgData, PNTR AlnSettingsDlgPtr;

static Pointer AlnSettingsDlgToData (DialoG d)
{
  AlnSettingsDlgPtr dlg;
  TSequenceInfoPtr  sequence_info;

  dlg = (AlnSettingsDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  sequence_info = SequenceInfoNew ();
  if (sequence_info == NULL) return NULL;

  sequence_info->missing = MemFree (sequence_info->missing);
  sequence_info->missing = SaveStringFromText (dlg->missing);

  sequence_info->beginning_gap = MemFree (sequence_info->beginning_gap);
  sequence_info->beginning_gap = SaveStringFromText (dlg->beginning_gap);

  sequence_info->middle_gap = MemFree (sequence_info->middle_gap);
  sequence_info->middle_gap = SaveStringFromText (dlg->middle_gap);

  sequence_info->end_gap = MemFree (sequence_info->end_gap);
  sequence_info->end_gap = SaveStringFromText (dlg->end_gap);

  sequence_info->match = MemFree (sequence_info->match);
  sequence_info->match = SaveStringFromText (dlg->match);

  if (dlg->sequence_type != NULL) 
  {
    if (GetValue (dlg->sequence_type) == 1) {
      sequence_info->alphabet = nucleotide_alphabet;
    } else {
      sequence_info->alphabet = protein_alphabet;
    }
  }
  else
  {
    sequence_info->alphabet = nucleotide_alphabet;
  }

  return sequence_info;
}


static void DataToAlnSettingsDlg (DialoG d, Pointer data)
{
  AlnSettingsDlgPtr dlg;
  TSequenceInfoPtr  sequence_info;

  dlg = (AlnSettingsDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  sequence_info = (TSequenceInfoPtr) data;

  if (sequence_info == NULL)
  {
    SetTitle (dlg->missing, "?Nn");
    SetTitle (dlg->beginning_gap, "-.Nn?");
    SetTitle (dlg->middle_gap, "-.");
    SetTitle (dlg->end_gap, "-.Nn?");
    SetTitle (dlg->match, ":");
    if (dlg->sequence_type != NULL)
    {
      SetValue (dlg->sequence_type, 1);
    }
  }
  else
  {
    SetTitle (dlg->missing, sequence_info->missing);
    SetTitle (dlg->beginning_gap, sequence_info->beginning_gap);
    SetTitle (dlg->middle_gap, sequence_info->middle_gap);
    SetTitle (dlg->end_gap, sequence_info->end_gap);
    SetTitle (dlg->match, sequence_info->match);

    if (dlg->sequence_type != NULL) 
    {
      if (StringCmp (sequence_info->alphabet, protein_alphabet) == 0) 
      {
        SetValue (dlg->sequence_type, 2);
      }
      else
      {
        SetValue (dlg->sequence_type, 1);
      }
    }
  }
}


static ValNodePtr TestAlnSettingsDlg (DialoG d)
{
  ValNodePtr        err_list = NULL;
  TSequenceInfoPtr  sequence_info;
  CharPtr           cp;
  CharPtr           fmt = "Character %c cannot appear in both %s and %s.";
  CharPtr           err_str;
  CharPtr           missing_name = "Ambiguous/Unknown";
  CharPtr           middle_gap_name = "Middle Gap";
  CharPtr           match_name = "Match";

  sequence_info = DialogToPointer (d);
  if (sequence_info == NULL) return NULL;

  /* missing and match cannot appear in middle gap list */
  cp = sequence_info->missing;
  while (cp != NULL && *cp != 0)
  {
    if (StringChr (sequence_info->middle_gap, *cp)) 
    {
      err_str = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) 
                                                   + StringLen (missing_name)
                                                   + StringLen (middle_gap_name)));
      sprintf (err_str, fmt, *cp, missing_name, middle_gap_name);
      ValNodeAddPointer (&err_list, 0, err_str);
    }
    cp++;
  }

  cp = sequence_info->match;
  while (cp != NULL && *cp != 0)
  {
    if (StringChr (sequence_info->middle_gap, *cp)) 
    {
      err_str = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) 
                                                   + StringLen (match_name)
                                                   + StringLen (middle_gap_name)));
      sprintf (err_str, fmt, *cp, match_name, middle_gap_name);
      ValNodeAddPointer (&err_list, 0, err_str);
    }
    cp++;
  }

  /* missing and match cannot share characters */
  cp = sequence_info->missing;
  while (cp != NULL && *cp != 0)
  {
    if (StringChr (sequence_info->match, *cp)) 
    {
      err_str = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) 
                                                   + StringLen (missing_name)
                                                   + StringLen (match_name)));
      sprintf (err_str, fmt, *cp, missing_name, match_name);
      ValNodeAddPointer (&err_list, 0, err_str);
    }
    cp++;
  }

  return err_list;
}


static CharPtr aln_settings_help = "\
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


NLM_EXTERN DialoG AlnSettingsDlg (GrouP h, Boolean allow_sequence_type)
{
  AlnSettingsDlgPtr dlg;
  GrouP             p, g, p_msg;
  
  dlg = (AlnSettingsDlgPtr) MemNew (sizeof (AlnSettingsDlgData));
  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);

  dlg->dialog = (DialoG) p;
  dlg->todialog = DataToAlnSettingsDlg;
  dlg->fromdialog = AlnSettingsDlgToData;
  dlg->testdialog = TestAlnSettingsDlg;

  g = HiddenGroup (p, 2, 0, NULL);
  StaticPrompt (g, "Ambiguous/Unknown", 0, dialogTextHeight, programFont, 'c');
  dlg->missing = DialogText (g, "?Nn", 5, NULL);
  StaticPrompt (g, "Match", 0, dialogTextHeight, programFont, 'c');
  dlg->match = DialogText (g, ".", 5, NULL);
  StaticPrompt (g, "Beginning Gap", 0, dialogTextHeight, programFont, 'c');
  dlg->beginning_gap = DialogText (g, "-.?nN", 5, NULL);
  StaticPrompt (g, "Middle Gap", 0, dialogTextHeight, programFont, 'c');
  dlg->middle_gap = DialogText (g, "-", 5, NULL);
  StaticPrompt (g, "End Gap", 0, dialogTextHeight, programFont, 'c');
  dlg->end_gap = DialogText (g, "-.?nN", 5, NULL);
  if (allow_sequence_type) {
    StaticPrompt (g, "Sequence Type", 0, dialogTextHeight, programFont, 'c');
    dlg->sequence_type = PopupList (g, TRUE, NULL);
    PopupItem (dlg->sequence_type, "Nucleotide");
    PopupItem (dlg->sequence_type, "Protein");
    SetValue (dlg->sequence_type, 1);
  }
  
  p_msg = MultiLinePrompt (p, aln_settings_help, 750, systemFont);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) p_msg, NULL);
  
  return (DialoG) p;
}


typedef struct alignmentoptionsform
{
  Boolean accepted;
  Boolean done;	
} AlignmentOptionsFormData, PNTR AlignmentOptionsFormPtr;

static void AcceptAlignmentOptions (ButtoN b)
{
  AlignmentOptionsFormPtr aofp;
  
  aofp = (AlignmentOptionsFormPtr) GetObjectExtra (b);
  if (aofp == NULL) return;
  aofp->accepted = TRUE;
  aofp->done = TRUE;
}

static void CancelAlignmentOptions (ButtoN b)
{
  AlignmentOptionsFormPtr aofp;
  
  aofp = (AlignmentOptionsFormPtr) GetObjectExtra (b);
  if (aofp == NULL) return;
  aofp->accepted = FALSE;
  aofp->done = TRUE;
}


NLM_EXTERN TSequenceInfoPtr GetAlignmentOptions (Uint1Ptr moltype, TSequenceInfoPtr sequence_info)
{
  ButtoN                   b;
  GrouP                    c, h;
  WindoW                   w;
  AlignmentOptionsFormData aofd;
  DialoG                   d;
  ValNodePtr               err_list;

  aofd.accepted = FALSE;
  aofd.done = FALSE;
  w = ModalWindow (-50, -33, -10, -10, NULL);

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  d = AlnSettingsDlg (h, moltype == NULL ? FALSE : TRUE);  
  PointerToDialog (d, sequence_info);

  c = HiddenGroup (h, 4, 0, NULL);
  b = DefaultButton (c, "Accept", AcceptAlignmentOptions);
  SetObjectExtra (b, &aofd, NULL);
  b = PushButton (c, "Cancel", CancelAlignmentOptions);
  SetObjectExtra (b, &aofd, NULL);

  AlignObjects (ALIGN_CENTER, (HANDLE) d, (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Update ();
 
  while (!aofd.done) {
    while (!aofd.done)
    {
      ProcessExternalEvent ();
      Update ();
    }
    ProcessAnEvent ();
    if (!aofd.accepted)
    {
      Remove (w);
  	  return NULL;
    }

    err_list = TestDialog (d);
    if (err_list != NULL) {
      Message (MSG_ERROR, err_list->data.ptrvalue);
      err_list = ValNodeFreeData (err_list);
      aofd.done = FALSE;
      aofd.accepted = FALSE;
    }
  }
  sequence_info = DialogToPointer (d);
  if (sequence_info == NULL) return NULL;
  if (StringCmp (sequence_info->alphabet, protein_alphabet) == 0) {
    if (moltype != NULL) {
      *moltype = Seq_mol_aa;
    }
  } else {
    if (moltype != NULL) {
      *moltype = Seq_mol_na;
    }
  }

  Remove (w);
  return sequence_info;
}


static void RemoveSequenceFromAlignmentFile (TAlignmentFilePtr afp, CharPtr str)
{
  int pos, k;

  if (afp == NULL) {
    return;
  }
  for (pos = 0; pos < afp->num_sequences; pos++) {
    if (StringCmp (str, afp->ids[pos]) == 0) {
      free (afp->ids[pos]);
      for (k = pos + 1; k < afp->num_sequences; k++) {
        afp->ids[k - 1] = afp->ids[k];
      }
      free (afp->sequences[pos]);
      for (k = pos + 1; k < afp->num_sequences; k++) {
        afp->sequences[k - 1] = afp->sequences[k];
      }
      afp->num_sequences--;
      if (pos < afp->num_organisms) {
        free (afp->organisms[pos]);
        for (k = pos + 1; k < afp->num_organisms; k++) {
          afp->organisms[k - 1] = afp->organisms[k];
        }
        afp->num_organisms--;
      }
      if (pos < afp->num_deflines) {
        free (afp->deflines[pos]);
        for (k = pos + 1; k < afp->num_deflines; k++) {
          afp->deflines[k - 1] = afp->deflines[k];
        }
        afp->num_deflines--;
      }
    }
  }
}


static void FixToFarPointer (TAlignmentFilePtr afp, Int4 index)
{
  CharPtr tmp_id_str;

  if (afp == NULL || index < -1 || index >= afp->num_sequences)
  {
    return;
  }
  tmp_id_str = (CharPtr) MemNew (sizeof (Char) * (StringLen (afp->ids [index]) + 4));
  if (tmp_id_str == NULL) 
  {
    return;
  }
  sprintf (tmp_id_str, "acc%s", afp->ids [index]);
  MemFree (afp->ids [index]);
  afp->ids[index] = tmp_id_str;
}


static void RemoveNthSequenceFromAlignment (TAlignmentFilePtr afp, Int4 n)
{
  Int4 i;
  
  if (afp == NULL || n < 0) {
    return;
  }
  
  if (afp->deflines != NULL && n < afp->num_deflines) {
    afp->deflines[n] = MemFree (afp->deflines[n]);
    for (i = n + 1; i < afp->num_deflines; i++) {
      afp->deflines[i - 1] = afp->deflines[i];
    }
    afp->deflines[afp->num_deflines - 1] = NULL;
    afp->num_deflines--;
  }
  
  if (afp->organisms != NULL && n < afp->num_organisms) {
    afp->organisms[n] = MemFree (afp->organisms[n]);
    for (i = n + 1; i < afp->num_organisms; i++) {
      afp->organisms[i - 1] = afp->organisms[i];
    }
    afp->organisms[afp->num_organisms - 1] = NULL;
    afp->num_organisms--;
  }
  
  if (afp->sequences != NULL && n < afp->num_sequences) {
    afp->sequences[n] = MemFree (afp->sequences[n]);
    afp->ids[n] = MemFree (afp->ids[n]);
    for (i = n + 1; i < afp->num_sequences; i++) {
      afp->sequences[i - 1] = afp->sequences[i];
      afp->ids[i - 1] = afp->ids[i];
    }
    afp->sequences[afp->num_sequences - 1] = NULL;
    afp->ids[afp->num_sequences - 1] = NULL;
    afp->num_sequences--;
  }
}


static void FixAlignmentIdsOkCancel (ButtoN b)
{
  BoolPtr bp;
  
  bp = (BoolPtr) GetObjectExtra (b);
  if (bp != NULL)
  {
    *bp = TRUE;
  }
}

static void EnableTextID (GrouP g)
{
  TexT id_text;
  
  id_text = (TexT) GetObjectExtra (g);
  if (id_text != NULL)
  {
    if (GetValue (g) > 5)
    {
      Enable (id_text);
    }
    else
    {
      Disable (id_text);
    }
  }
}


static Boolean ReplaceAlignmentIDsFromFile (TAlignmentFilePtr afp, Int4 index)
{
  ReadBufferData    rbd;
  Char              path [PATH_MAX];
  CharPtr           line, cp, first_id, second_id;
  Int4              k;
  Boolean           found_id;
  ValNodePtr        err_list = NULL;
  CharPtr           err_msg = NULL;
  CharPtr           err_msg_prefix = "Unable to find ";
  CharPtr           err_msg_suffix = " from file in alignment";
  Int4              err_msg_len = 0;
  ValNodePtr        vnp;

  if (afp == NULL || index < -1 || index >= afp->num_sequences)
  {
    return FALSE;
  }
  
  rbd.fp = NULL;
  while (rbd.fp == NULL)
  {
    if (!GetInputFileName (path, sizeof (path), NULL, NULL))
    {
      return FALSE;
    }
    rbd.fp = FileOpen (path, "r");
    if (rbd.fp == NULL)
    {
      Message (MSG_ERROR, "Unable to open %s", path);
    }
  }
    
  rbd.current_data = NULL;
  
  line = AbstractReadFunction (&rbd);
  while (line != NULL)
  {
    cp = line;
    while (isspace ((Int4)(*cp)))
    {
      cp++;
    }
    if (*cp != 0)
    {
      first_id = cp;
      while (!isspace ((Int4)(*cp)) && *cp != 0)
      {
        cp++;
      }
      while (isspace ((Int4)(*cp)))
      {
        *cp = 0;
        cp++;
      }
      second_id = cp;
      TrimSpacesAroundString (second_id); 
      if (*second_id != 0)
      {
        found_id = FALSE;
        for (k = index; k < afp->num_sequences && ! found_id; k++)
        {
          if (StringCmp (afp->ids[k], first_id) == 0)
          {
            found_id = TRUE;
            MemFree (afp->ids[k]);
            afp->ids[k] = StringSave (second_id);
          }
          else if (StringCmp (afp->ids[k], second_id) == 0)
          {
            found_id = TRUE;
            MemFree (afp->ids[k]);
            afp->ids[k] = StringSave (first_id);
          }
        }
        if (!found_id)
        {
          ValNodeAddPointer (&err_list, 0, StringSave (first_id));
          ValNodeAddPointer (&err_list, 0, StringSave (second_id));
          err_msg_len += StringLen (first_id) + StringLen (second_id) + 8;
        }
      }
    }
    
    line = AbstractReadFunction (&rbd);
  }
  FileClose (rbd.fp);
  if (err_list != NULL)
  {
    err_msg_len += StringLen (err_msg_prefix) + StringLen (err_msg_suffix) + 6;
    err_msg = (CharPtr) MemNew (err_msg_len * sizeof (Char));
    if (err_msg != NULL)
    {
      sprintf (err_msg, "%s", err_msg_prefix);
      vnp = err_list;
      while (vnp != NULL && vnp->next != NULL)
      {
        StringCat (err_msg, (CharPtr) vnp->data.ptrvalue);
        StringCat (err_msg, " or ");
        StringCat (err_msg, (CharPtr) vnp->next->data.ptrvalue);
        
        if (vnp->next->next != NULL)
        {
          StringCat (err_msg, ", ");
        }
        if (vnp->next->next != NULL && vnp->next->next->next != NULL && vnp->next->next->next->next == NULL)
        {
          StringCat (err_msg, " and ");
        }
        vnp = vnp->next->next;
      }
      StringCat (err_msg, err_msg_suffix);
      Message (MSG_ERROR, err_msg);
      MemFree (err_msg);
    }
  }
  return TRUE; 
}


static Boolean FixAlignmentIDs (TAlignmentFilePtr afp, Int4 index, BoolPtr all_far, BoolPtr all_skip, BoolPtr removed)
{
  WindoW  w;
  GrouP   h, choice_grp, c;
  ButtoN  b;
  Boolean done = FALSE;
  Boolean cancelled = FALSE;
  PrompT  p;
  CharPtr prompt_str;
  CharPtr prompt_str_fmt = "Unable to find sequence %s from alignment in set.";
  TexT    id_text;
  CharPtr id_str;
  Int2    fix_choice;
  
  if (afp == NULL || index < -1 || index >= afp->num_sequences)
  {
    return FALSE;
  }
  
  id_str = afp->ids[index];
  if (StringHasNoText (id_str))
  {
    return FALSE;
  }
  
  w = MovableModalWindow (-20, -13, -10, -10, "Source Assistant", NULL);
  h = HiddenGroup(w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  prompt_str = (CharPtr) MemNew ((StringLen (prompt_str_fmt) + StringLen (id_str)) * sizeof (Char));
  if (prompt_str != NULL)
  {
    sprintf (prompt_str, prompt_str_fmt, id_str);
  }
  p = StaticPrompt (h, prompt_str, 0, dialogTextHeight, systemFont, 'c');
  choice_grp = HiddenGroup (h, 0, 7, EnableTextID);
  RadioButton (choice_grp, "This is a far pointer");
  RadioButton (choice_grp, "All unmatched sequences are far pointers");
  RadioButton (choice_grp, "Read in a file that maps alignment IDs to sequence IDs");
  RadioButton (choice_grp, "Remove this sequence from the alignment");
  RadioButton (choice_grp, "Remove all unmatched sequences from the alignment");
  RadioButton (choice_grp, "Use this ID for this sequence");
  id_text = DialogText (choice_grp, "", 20, NULL);
  Disable (id_text);
  SetValue (choice_grp, 1);
  SetObjectExtra (choice_grp, id_text, NULL);
  c = HiddenGroup (h, 2, 0, NULL);
  b = PushButton(c, "OK", FixAlignmentIdsOkCancel);
  SetObjectExtra (b, &done, NULL);
  b = PushButton(c, "Cancel", FixAlignmentIdsOkCancel);
  SetObjectExtra (b, &cancelled, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) p,  (HANDLE) choice_grp, (HANDLE) c, NULL);
  
  Show(w); 
  Select (w);
  while (!done && !cancelled)
  {
    while (!done && !cancelled)
    {
      ProcessExternalEvent ();
      Update ();
    }
    ProcessAnEvent ();
    if (!cancelled)
    {
      fix_choice = GetValue (choice_grp);
      switch (fix_choice)
      {
        case 1:
          /* far pointer */
          FixToFarPointer (afp, index);
          break;
        case 2:
          /* all far pointers */
          if (all_far != NULL)
          {
            *all_far = TRUE;
          }
          FixToFarPointer (afp, index);
          break;
        case 3:
          /* read in file with replacements */
          if (!ReplaceAlignmentIDsFromFile (afp, index))
          {
            done = FALSE;
          }
          break;
        case 4:
          /* skip */
          RemoveNthSequenceFromAlignment (afp, index);
          if (removed != NULL) {
            *removed = TRUE;
          }     
          break;
        case 5:
          /* skip all */
          RemoveNthSequenceFromAlignment (afp, index);          
          if (all_skip != NULL)
          {
            *all_skip = TRUE;
          }
          if (removed != NULL) {
            *removed = TRUE;
          }     
          break;
        case 6:
          /* use single replacement */
          id_str = SaveStringFromText (id_text);
          if (StringHasNoText (id_str))
          {
            MemFree (id_str);
            Message (MSG_ERROR, "You did not specify text for the ID!");
            done = FALSE;
          }
          else
          {
            MemFree (afp->ids [index]);
            afp->ids [index] = id_str;
          }
          break;
      }
    }
  }
  Remove (w);
  if (cancelled)
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }
}


NLM_EXTERN Boolean CorrectAlignmentIDs (TAlignmentFilePtr afp, Uint1 moltype)
{
  Int4      index;
  CharPtr   seq_data;
  SeqIdPtr  sip;
  BioseqPtr bsp;
  CharPtr   tmp_id_str;
  Char      prot_str[200];
  Boolean   all_far = FALSE;
  Boolean   all_skip = FALSE;
  Boolean   removed;
  SeqEntryPtr nucprot_sep;
  BioseqSetPtr nucprot_bssp;
  
  for (index = 0; index < afp->num_sequences; index++) {
    seq_data = AlignmentStringToSequenceString (afp->sequences [index], moltype);
    if (! StringHasNoText (seq_data))
    {
      sip = MakeSeqID (afp->ids [index]);
      sip->next = SeqIdFree (sip->next);
      if (StringNCmp (afp->ids[index], "acc", 3) != 0)
      {
        bsp = BioseqFind (sip);
        if (bsp == NULL && StringChr (afp->ids[index], '|') == NULL)
        {
          sip = SeqIdFree (sip);
          tmp_id_str = (CharPtr) MemNew (sizeof (Char) * (StringLen (afp->ids [index]) + 4));
          sprintf (tmp_id_str, "gb|%s", afp->ids [index]);
          sip = MakeSeqID (tmp_id_str);
          MemFree (tmp_id_str);
          bsp = BioseqFind (sip);
        }

        if (bsp != NULL && moltype == Seq_mol_aa && !ISA_aa (bsp->mol)) 
        {
          /* IDs in alignment are for nucleotide sequences but this is 
           * protein sequence - if single prot in nuc-prot set, replace
           * with protein sequence ID.  Otherwise set to NULL - don't
           * want this ID in our alignment */
          nucprot_sep = GetBestTopParentForData (bsp->idx.entityID, bsp);
          bsp = NULL;
          if (nucprot_sep != NULL && IS_Bioseq_set (nucprot_sep) && nucprot_sep->data.ptrvalue != NULL) 
          {
            nucprot_bssp = nucprot_sep->data.ptrvalue;
            if (nucprot_bssp->seq_set != NULL 
                && nucprot_bssp->seq_set->next != NULL 
                && nucprot_bssp->seq_set->next->next == NULL
                && IS_Bioseq (nucprot_bssp->seq_set->next)) 
            {
              sip = SeqIdFree (sip);
              bsp = (BioseqPtr) nucprot_bssp->seq_set->next->data.ptrvalue;
              SeqIdWrite (SeqIdFindBest (bsp->id, 0), prot_str, PRINTID_FASTA_LONG, sizeof (prot_str) - 1);
              afp->ids[index] = MemFree (afp->ids[index]);
              afp->ids[index] = StringSave (prot_str);
            }
          }
        }

        if (bsp == NULL)
        {
          if (all_far)
          {
            FixToFarPointer (afp, index);
          }
          else if (all_skip) 
          {
            RemoveNthSequenceFromAlignment(afp, index);
            index--;
          }
          else
          {
            removed = FALSE;
            if (!FixAlignmentIDs (afp, index, &all_far, &all_skip, &removed))
            {
              /* bail - user does not want to fix IDs */
              return FALSE;
            }
            if (removed) {
              index--;
            }
          }
        }

      }
      sip = SeqIdFree (sip);
    }
    MemFree (seq_data);
  }
  return TRUE;
}


NLM_EXTERN SeqAlignPtr ReadAlignmentForSeqEntry (SeqEntryPtr sep, Boolean is_nuc, Boolean allow_options, Boolean from_clipboard)
{
  Char              path [PATH_MAX];
  FILE              *fp;
  SeqAlignPtr       salp=NULL,
                    salpnew, salp_return = NULL;
  SeqEntryPtr       sepnew=NULL;
  Boolean           ok = TRUE, 
                    dirty = FALSE;
  TSequenceInfoPtr  default_info, sequence_info;
  ReadBufferData    rbd;
  TErrorInfoPtr     error_list;
  TAlignmentFilePtr afp;
  Uint1             moltype;
  ErrSev            sev;
  SeqEntryPtr       scope;
  CharPtr           str;

  if (sep == NULL) return NULL;

  if (from_clipboard) {
    if (!Nlm_ClipboardHasString()) {
      Message (MSG_ERROR, "Clipboard is empty!");
      return NULL;
    }
    TmpNam (path);
    fp = FileOpen (path, "w");
    str = ClipboardToString();
    fprintf (fp, "%s", str);
    FileClose (fp);
    str = MemFree (str);
  } else if (! GetInputFileName (path, sizeof (path),"","TEXT")) {
    return NULL;
  }
  fp = FileOpen (path, "r");
  if (fp == NULL)
  {
    Message (MSG_ERROR, "Unable to open %s", path);
    if (from_clipboard) {
      FileRemove (path);
    }
    return NULL;
  }

  default_info = SequenceInfoNew ();
  if (is_nuc) {
    default_info->alphabet = nucleotide_alphabet;
  } else {
    default_info->alphabet = protein_alphabet;
  }
  default_info->beginning_gap = MemFree (default_info->beginning_gap);
  default_info->beginning_gap = StringSave ("-.Nn?");
  default_info->middle_gap = MemFree (default_info->middle_gap);
  default_info->middle_gap = StringSave ("-.#");
  default_info->end_gap = MemFree (default_info->end_gap);
  default_info->end_gap = StringSave ("-.Nn?");
  default_info->match = MemFree (default_info->match);
  default_info->match = StringSave (":");
  default_info->missing = MemFree (default_info->missing);
  default_info->missing = StringSave ("?Nn");

  if (allow_options) {
    sequence_info = GetAlignmentOptions (&moltype, default_info);
    SequenceInfoFree (default_info);
    if (sequence_info == NULL) return NULL;
    default_info = sequence_info;
  }

  WatchCursor();
  error_list = NULL;

  rbd.fp = fp;
  rbd.current_data = NULL;
  afp = ReadAlignmentFile ( AbstractReadFunction,
                            (Pointer) &rbd,
                            AbstractReportError,
                            (Pointer) &error_list,
                            default_info);
  FileClose (fp);
  if (from_clipboard) {
    FileRemove (path);
  }
  SequenceInfoFree (default_info);
  default_info = NULL;
  if (afp != NULL) 
  {
    RemoveSequenceFromAlignmentFile (afp, "Consensus");
    SeqEntrySetScope (sep);
    if (CorrectAlignmentIDs (afp, moltype))
    {
      sepnew = MakeSequinDataFromAlignmentEx (afp, moltype, TRUE); 
    }
  }
  if (sepnew == NULL) {
    ProduceAlignmentNotes (afp, error_list);
  }
  ErrorInfoFree (error_list);
  AlignmentFileFree (afp);
  if (sepnew) 
  {
    salpnew = (SeqAlignPtr) FindSeqAlignInSeqEntry (sepnew, OBJ_SEQALIGN);
    if (salpnew) {
      sev = ErrSetMessageLevel (SEV_FATAL);

      scope = SeqEntrySetScope (NULL);
      /* adjust the start positions for the sequences read in from the alignments. */
      CalculateAlignmentOffsets (sepnew, sep);
      /* ValidateSeqAlignandACCInSeqEntry will readjust the start positions for
       * the alignments for far pointer sequences.
       */
      SeqEntrySetScope (scope);
      ok = ValidateSeqAlignandACCInSeqEntry (sepnew, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE);
      
      ErrSetMessageLevel (sev);
      
      if (ok) {
        AlnMgr2IndexSeqAlignEx(salpnew, FALSE);
        ok = CheckAlignmentSequenceLengths (salpnew);
      }
      
      if (ok) {
        /* make copy, otherwise it will be removed when we delete sep_new */
        salp_return = (SeqAlignPtr) AsnIoMemCopy (salpnew, (AsnReadFunc) SeqAlignAsnRead,
                                                   (AsnWriteFunc) SeqAlignAsnWrite);
      }
    }
    /* this statement will free salpnew as part of sepnew */
    ObjMgrFree (OBJ_SEQENTRY, (Pointer)sepnew);
    sepnew=NULL;
  }
  ArrowCursor();
  Update ();
  return salp_return;
}


static void ReportPotentialDupIDs (TAlignmentFilePtr afp, FILE *fp)
{
  int     seq_index, k;
  int     curr_seg;
  int     num_sequences;
  Int4Ptr seq_len;
  BoolPtr may_be_dup;
  Int4    a, b;

  if (afp == NULL || afp->sequences == NULL || afp->num_sequences == 0) {
    return;
  }

  num_sequences = afp->num_sequences / afp->num_segments;

  may_be_dup = (BoolPtr) MemNew (sizeof (Boolean) * num_sequences);
  for (seq_index = 0; seq_index < num_sequences; seq_index++) {
    may_be_dup[seq_index] = FALSE;
  }

  seq_len = (Int4Ptr) MemNew (sizeof (Int4) * afp->num_sequences);
  for (seq_index = 0; seq_index < afp->num_sequences; seq_index++) {
    seq_len[seq_index] = StringLen (afp->sequences[seq_index]);
  }

  for (curr_seg = 0; curr_seg < afp->num_segments; curr_seg++) {
    for (seq_index = 0; seq_index < num_sequences - 1; seq_index++) {
      for (k = seq_index + 1; k < num_sequences; k++) {
        a = curr_seg * num_sequences + seq_index;
        b = curr_seg * num_sequences + k;
        if (seq_len[a] != seq_len[b]) {
          if (seq_len[a] % seq_len [b] == 0) {
            may_be_dup[seq_index] = TRUE;
          } else if (seq_len[b] % seq_len[a] == 0) {
            may_be_dup[k] = TRUE;
          }
        }
      }
    }
  }
  seq_len = MemFree (seq_len);

  for (seq_index = 0; seq_index < num_sequences; seq_index ++)
  {
    if (may_be_dup[seq_index]) {
      fprintf (fp, "Please check your file - %s may have been used as an ID for multiple sequences.\n", afp->ids[seq_index]);
    }
  }

  may_be_dup = MemFree (may_be_dup);
}


static void PrintExtraErrorInstructions (FILE *fp, CharPtr message)
{
  CharPtr explanation, end;
  Char    tmp = '\0';
  if (fp == NULL || message == NULL) return;

  if (StringStr (message, "bad characters") == NULL
      && StringStr (message, " found at position ") == NULL) {
    return;
  }

  explanation = StringRChr (message, '(');
  if (explanation == NULL) return;
  if (StringNCmp (explanation, "(expect only ", 13) == 0) {
    end = StringStr (explanation + 13, " here)");
    if (end != NULL) {
      tmp = *end;
      *end = 0;
    }
    fprintf (fp, 
             "Try changing the sequence character specifications for %s.\n",
             explanation + 13);
    if (StringNCmp (explanation + 13, "beginning", 9) == 0) {
      fprintf (fp, 
"\nWhen some of the sequences in an alignment are shorter or "
"longer than others, beginning gap characters are added to "
"the beginning of the sequence to maintain the correct spacing."
"  These will not appear in your sequence file.\n");
    } else if (StringNCmp (explanation + 13, "end", 3) == 0) {
      fprintf (fp,
"\nWhen some of the sequences in an alignment are shorter or "
"longer than others, end gap characters are added to "
"the end of the sequence to maintain the correct spacing."
"  These will not appear in your sequence file.\n");
    } else {
      fprintf (fp,
"\nMiddle gap characters are used to maintain the spacing "
"inside an alignment.  These are not nucleotides and will "
"not appear as part of your sequence file.\n"
"Ambiguous/unknown characters are used to represent indeterminate/ambiguous "
"nucleotides.  These will appear in your sequence file as 'n'.\n"
"Match characters are used to indicate positions where "
"sequences are identical to the first sequence.  These will be "
"replaced by the actual characters from the first sequence.\n");
    }
    if (end != NULL) {
      *end = tmp;
    }
  } else if (StringCmp (explanation, 
                        "(can't specify match chars in first sequence).") == 0) {
    fprintf (fp, "Try changing the match character specification.\n");
    fprintf (fp,
"\nMatch characters are used to indicate positions where "
"sequences are identical to the first sequence.  These will be "
"replaced by the actual characters from the first sequence.\n");
  }
}


static void PrintError (FILE *fp, TErrorInfoPtr eip)
{
  if (eip == NULL || fp == NULL) return;

  fprintf (fp, "*****\nError category %d\n", eip->category);
  if (eip->line_num > -1) {
    fprintf (fp, "Line number %d\n", eip->line_num);
  }
  if (eip->id != NULL) {
    fprintf (fp, "Sequence ID %s\n", eip->id);
  }
  if (eip->message != NULL) {
    fprintf (fp, "%s\n", eip->message);
    PrintExtraErrorInstructions (fp, eip->message);
  }
}
  
static Int4 CountNucleotides (CharPtr sequence)
{
  Int4    num = 0;
  CharPtr cp;
  
  if (sequence == NULL) return 0;	
  for (cp = sequence; *cp != 0; cp++)
  {
  	if (*cp != '-') 
  	{
  	  num++;
  	}
  }
  return num;
}

static void WalkErrorList (TErrorInfoPtr list, FILE *fp)
{
  TErrorInfoPtr eip;
  
  if (list == NULL || fp == NULL) return;

  for (eip = list; eip != NULL; eip = eip->next) {
    PrintError (fp, eip);
  }

}


static void PrintAlignmentSummary (TAlignmentFilePtr afp, FILE *fp)
{
  Int4 index;

  if (fp == NULL) return;

  if (afp == NULL) {
    fprintf (fp, "Catastrophic failure during reading\n");
  } else {
    fprintf (fp, "Found %d sequences\n", afp->num_sequences);
    fprintf (fp, "Found %d organisms\n", afp->num_organisms);
    if (afp->num_sequences == afp->num_segments * afp->num_organisms)
    {
      for (index = 0; index < afp->num_sequences; index++)
      {
        fprintf (fp, "\t%s\t%d nucleotides\t", afp->ids [index],
                 CountNucleotides (afp->sequences[index]));
        if (index / afp->num_segments < afp->num_organisms) {
          fprintf (fp, "%s\n", afp->organisms [index / afp->num_segments]);
        } else {
          fprintf (fp, "No organism information\n");
        }
      }    	
    }
    else
    {
      for (index = 0; index < afp->num_sequences; index++)
      {
        fprintf (fp, "\t%s\t%d nucleotides\t", afp->ids [index], 
                 CountNucleotides (afp->sequences[index]));
        if (index < afp->num_organisms) {
          fprintf (fp, "%s\n", afp->organisms [index]);
        } else {
          fprintf (fp, "No organism information\n");
        }
      }
      while (index < afp->num_organisms) {
        fprintf (fp, "Unclaimed organism: %s\n", afp->organisms [index]);
        index++;
      }	
    }
  }
}


NLM_EXTERN void 
ProduceAlignmentNotes 
(TAlignmentFilePtr afp,
 TErrorInfoPtr error_list)
{
  Char         path [PATH_MAX];
  FILE         *fp;
  Boolean      ok_to_import = FALSE;

  TmpNam (path);
  fp = FileOpen (path, "wb");
  if (fp == NULL) return;


  if (afp != NULL && DoSequenceLengthsMatch (afp)) {
    ok_to_import = TRUE;
  }

  
  if (ok_to_import && error_list != NULL) {
    fprintf (fp, "Congratulations, you have successfully created a sequin file;\nhowever, I had trouble reading part of your file.\nPlease check your data carefully before submitting to be sure that all of your sequences\nwere included correctly.\nIf your file is incomplete, or contains incorrect sequences, please use the error report below\nto find the problem.\n");
  }
  ReportPotentialDupIDs (afp, fp);

  WalkErrorList (error_list, fp);
  PrintAlignmentSummary (afp, fp);

  FileClose (fp);
  LaunchGeneralTextViewer (path, "Alignment reading summary");
  FileRemove (path);
}

