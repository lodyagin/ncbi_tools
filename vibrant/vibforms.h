/*   vibforms.h
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
* File Name:  vibforms.h
*
* Author:  Jonathan Kans, Sergei Egorov (EnumPopup code)
*
* Version Creation Date:   1/22/95
*
* $Revision: 6.22 $
*
* File Description: 
*
* Modifications:  
* --------------------------------------------------------------------------
*
* ==========================================================================
*/

#ifndef _VIBFORMS_
#define _VIBFORMS_

#ifndef _VIBRANT_
#include <vibrant.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif


typedef Nlm_Int4 Nlm_UIEnum, PNTR Nlm_UIEnumPtr;

typedef struct Nlm_enum_field_assoc {
   Nlm_CharPtr name; 
   Nlm_UIEnum value;
} Nlm_EnumFieldAssoc, PNTR Nlm_EnumFieldAssocPtr, Nlm_EnumFieldAlist[];

#define ENUM_ALIST(name) Nlm_EnumFieldAssoc name[] = {
#define END_ENUM_ALIST { NULL, 0 } };

#ifdef WIN_MAC
extern Pointer Nlm_currentFormDataPtr;
extern Nlm_IteM PNTR Nlm_globalMenuItemList;
extern Nlm_Int2 Nlm_globalMenuListSize;
#endif

/* enum value -> name */
extern Nlm_CharPtr Nlm_GetEnumName (Nlm_UIEnum val, Nlm_EnumFieldAssocPtr al);

/* enum field <-> popup list UI */
/* returns TRUE and enum in pdef if default is selected */
extern Nlm_Boolean Nlm_InitEnumPopup (Nlm_PopuP lst, Nlm_EnumFieldAssocPtr al, Nlm_UIEnumPtr pdef);
extern void Nlm_SetEnumPopup (Nlm_PopuP lst, Nlm_EnumFieldAssocPtr al, Nlm_UIEnum val);
/* returns TRUE and enum in pdef if selection exists */
extern Nlm_Boolean Nlm_GetEnumPopup (Nlm_PopuP lst, Nlm_EnumFieldAssocPtr al, Nlm_UIEnumPtr pval);

extern Nlm_CharPtr Nlm_GetEnumPopupByName (Nlm_PopuP lst, Nlm_EnumFieldAssocPtr al);
extern void Nlm_SetEnumPopupByName (Nlm_PopuP lst, Nlm_EnumFieldAssocPtr al, Nlm_CharPtr name);

extern Nlm_Boolean Nlm_WhereInEnumPopup (Nlm_EnumFieldAssocPtr al, Nlm_CharPtr name, Nlm_UIEnumPtr pval);
extern void Nlm_SortEnumFieldAlist (Nlm_EnumFieldAssocPtr alist);
extern Nlm_EnumFieldAssocPtr Nlm_DuplicateEnumFieldAlist (Nlm_EnumFieldAssocPtr alist);
extern Nlm_EnumFieldAssocPtr Nlm_FreeEnumFieldAlist (Nlm_EnumFieldAssocPtr alist);
extern Nlm_EnumFieldAssocPtr Nlm_MakeEnumFieldAlistFromValNodeList (ValNodePtr vlist);

/* convenience functions that create the popup, set extra object data, init the popup, and set the value */
extern Nlm_PopuP Nlm_CreateEnumPopupListInitVal (Nlm_GrouP prnt, Nlm_Boolean macLike, Nlm_PupActnProc actn,
                                                 Nlm_VoidPtr data, Nlm_EnumFieldAssocPtr al, Nlm_UIEnum val);
extern Nlm_PopuP Nlm_CreateEnumPopupListInitName (Nlm_GrouP prnt, Nlm_Boolean macLike, Nlm_PupActnProc actn,
                                                  Nlm_VoidPtr data, Nlm_EnumFieldAssocPtr al, Nlm_CharPtr name);

typedef void (*Nlm_HandleActnProc) PROTO((Nlm_Handle));
void CDECL Nlm_RepeatProcOnHandles VPROTO((Nlm_HandleActnProc proc, ...));

extern void Nlm_SafeShow (Nlm_Handle a);
extern void Nlm_SafeHide (Nlm_Handle a);
extern void Nlm_SafeEnable (Nlm_Handle a);
extern void Nlm_SafeDisable (Nlm_Handle a);
extern void Nlm_SafeSetValue (Nlm_Handle a, Nlm_Int2 value);
extern void Nlm_SafeSetStatus (Nlm_Handle a, Nlm_Boolean status);
extern void Nlm_SafeSetTitle (Nlm_Handle a, Nlm_CharPtr title);

extern Nlm_CharPtr Nlm_SaveStringFromText (Nlm_TexT t);
extern Nlm_CharPtr Nlm_JustSaveStringFromText (Nlm_TexT t);
extern ValNodePtr Nlm_SetTextFromVnp (Nlm_TexT t, ValNodePtr vnp);
extern ValNodePtr Nlm_GetVnpFromText (Nlm_TexT t, ValNodePtr vnp, Nlm_Boolean last);
extern Nlm_Int2 Nlm_MaxStringWidths (Nlm_CharPtr PNTR strs);
extern Nlm_CharPtr Nlm_WidestString (Nlm_CharPtr PNTR strs);
extern Nlm_Int2 Nlm_MaxAlistWidths (Nlm_EnumFieldAssocPtr al);
extern Nlm_CharPtr Nlm_WidestAlist (Nlm_EnumFieldAssocPtr al);
extern Nlm_GrouP Nlm_MultiLinePrompt (Nlm_GrouP prnt, Nlm_CharPtr text, Nlm_Int2 maxWidth, Nlm_FonT font);
extern Nlm_GrouP Nlm_MultiLinePromptEx (Nlm_GrouP prnt, Nlm_CharPtr text, Nlm_Int2 maxWidth, Nlm_FonT font, Nlm_Boolean stripSpaces);

extern void Nlm_JustInvalObject (Nlm_Handle a);

/* extern Nlm_CharPtr Nlm_TrimSpacesAroundString (Nlm_CharPtr str);  - now in ncbistr.h */
/* extern Nlm_Boolean Nlm_StringHasNoText (Nlm_CharPtr str); - now in ncbistr.h */
extern Nlm_Boolean Nlm_TextHasNoText (Nlm_TexT t);

/*****************************************************************************
*
*   DialoG is a stand-alone Vibrant complex object.
*
*****************************************************************************/

typedef  struct  Nlm_dialog {
  Nlm_VoidPtr  dummy;
} HNDL Nlm_DialoG;

typedef void (*Nlm_ToDialogFunc) (Nlm_DialoG, Pointer);
typedef Pointer (*Nlm_FromDialogFunc) (Nlm_DialoG);
typedef void (*Nlm_DialogActnFunc) (Nlm_DialoG);
typedef ValNodePtr (*Nlm_DialogTestFunc) (Nlm_DialoG);
typedef void (*Nlm_DialogMessageFunc) (Nlm_DialoG, Nlm_Int2);
typedef Nlm_Boolean (*Nlm_DialogInOutFunc) (Nlm_DialoG, Nlm_CharPtr filename);

/*
  Place DIALOG_MESSAGE_BLOCK at the head of the
  dialog-specific extra object data structure
*/

#define DIALOG_MESSAGE_BLOCK       \
  Nlm_DialoG             dialog;          \
  Nlm_DialogActnFunc     actproc;         \
  Nlm_ToDialogFunc       todialog;        \
  Nlm_FromDialogFunc     fromdialog;      \
  Nlm_DialogTestFunc     testdialog;      \
  Nlm_DialogMessageFunc  dialogmessage;   \
  Nlm_DialogInOutFunc    importdialog;    \
  Nlm_DialogInOutFunc    exportdialog;    \
  Nlm_Int4               intvalue;        \
  Nlm_VoidPtr            userdata;

typedef struct Nlm_basedialog {
  DIALOG_MESSAGE_BLOCK
} PNTR Nlm_BaseDialogPtr;

extern void Nlm_PointerToDialog (Nlm_DialoG d, Pointer data);
extern Pointer Nlm_DialogToPointer (Nlm_DialoG d);
extern ValNodePtr Nlm_TestDialog (Nlm_DialoG d);
extern void Nlm_SendMessageToDialog (Nlm_DialoG d, Nlm_Int2 mssg);
extern Nlm_Boolean Nlm_ImportDialog (Nlm_DialoG d, Nlm_CharPtr filename);
extern Nlm_Boolean Nlm_ExportDialog (Nlm_DialoG d, Nlm_CharPtr filename);

extern void Nlm_SetDialogActnProc (Nlm_DialoG d, Nlm_DialogActnFunc actproc);

/* popup list autonomous dialogs - copies alist, frees on cleanup */

typedef struct alistdialogdata {
  DIALOG_MESSAGE_BLOCK
  Nlm_PopuP                pop;
  Nlm_EnumFieldAssoc PNTR  alist;
} Nlm_AlistDialogData, PNTR Nlm_AlistDialogPtr;

/* convenience function that creates enumerated popup as a dialog, copies alist, cleans up on freeing */
extern Nlm_PopuP Nlm_CreateEnumPopupDialog (Nlm_GrouP prnt, Nlm_Boolean macLike, Nlm_PupActnProc actn,
                                            Nlm_EnumFieldAssocPtr al, Nlm_UIEnum val, Pointer userdata);

/* convenience function that creates enumerated single list as a dialog, copies alist, cleans up on freeing */
extern Nlm_LisT Nlm_CreateEnumListDialog (Nlm_GrouP prnt, Nlm_Int2 width, Nlm_Int2 height, Nlm_LstActnProc actn,
                                          Nlm_EnumFieldAssocPtr al, Nlm_UIEnum val, Pointer userdata);
/*****************************************************************************
*
*   A form is a general collection of associated Vibrant objects.
*   Like the DialoG, a handle to which ObjectExtra data is attached,
*   the ForM is itself a window to which data is attached.
*
*****************************************************************************/

typedef  struct  Nlm_form {
  Nlm_VoidPtr  dummy;
} HNDL Nlm_ForM;

typedef void (*Nlm_ToFormFunc) (Nlm_ForM, Pointer);
typedef Pointer (*Nlm_FromFormFunc) (Nlm_ForM);
typedef void (*Nlm_FormActnFunc) (Nlm_ForM);
typedef ValNodePtr (*Nlm_FormTestFunc) (Nlm_ForM);
typedef void (*Nlm_FormMessageFunc) (Nlm_ForM, Nlm_Int2);
typedef Nlm_Boolean (*Nlm_FormInOutFunc) (Nlm_ForM, Nlm_CharPtr filename);

#define FORM_MESSAGE_BLOCK         \
  Nlm_ForM             form;            \
  Nlm_FormActnFunc     actproc;         \
  Nlm_ToFormFunc       toform;          \
  Nlm_FromFormFunc     fromform;        \
  Nlm_FormTestFunc     testform;        \
  Nlm_FormMessageFunc  formmessage;     \
  Nlm_FormInOutFunc    importform;      \
  Nlm_FormInOutFunc    exportform;      \
  Nlm_IteM PNTR        menuitemlist;    \
  Nlm_Int2             menulistsize;    \
  Nlm_FormMessageFunc  appmessage;      \
  Nlm_Uint2            input_entityID;  \
  Nlm_Uint4            input_itemID;    \
  Nlm_Uint2            input_itemtype;  \
  Nlm_Uint2            this_itemtype;   \
  Nlm_Uint2            this_subtype;    \
  Nlm_Uint2            procid;          \
  Nlm_Uint2            proctype;        \
  BIG_ID               docuid;          \
  Nlm_Int2             doctype;         \
  Nlm_VoidPtr          objectDataPtr;   \
  Nlm_VoidPtr          userDataPtr;     \
  Nlm_FreeProc         cleanupuser;     \
  Nlm_Uint2            userkey;         \
  Nlm_CharPtr          filepath;        \
  WndActnProc          activate;

typedef struct Nlm_baseform {
  FORM_MESSAGE_BLOCK
} PNTR Nlm_BaseFormPtr;

extern void Nlm_PointerToForm (Nlm_ForM f, Pointer data);
extern Pointer Nlm_FormToPointer (Nlm_ForM f);
extern ValNodePtr Nlm_TestForm (Nlm_ForM f);
extern void Nlm_SendMessageToForm (Nlm_ForM f, Nlm_Int2 mssg);
extern Nlm_Boolean Nlm_ImportForm (Nlm_ForM f, Nlm_CharPtr filename);
extern Nlm_Boolean Nlm_ExportForm (Nlm_ForM f, Nlm_CharPtr filename);


/* The following functions remove a window */
extern void Nlm_StdCloseWindowProc (Nlm_WindoW w);
extern void Nlm_StdCancelButtonProc (Nlm_ButtoN b);

/* The following functions instead send a VIB_MSG_CLOSE message to a form */
extern void Nlm_StdSendCloseWindowMessageProc (Nlm_WindoW w);
extern void Nlm_StdSendCancelButtonMessageProc (Nlm_ButtoN b);

extern void Nlm_SetFormActnProc (Nlm_ForM f, Nlm_FormActnFunc actproc);

/* StdAcceptFormButtonProc hides the form, calls the actproc, and removes the form */
/* StdSendAcceptButtonMessageProc instead sends a VIB_MSG_ACCEPT message to a form */
extern void Nlm_StdAcceptFormButtonProc (Nlm_ButtoN b);
extern void Nlm_StdSendAcceptButtonMessageProc (Nlm_ButtoN b);

/* StdCleanupExtraProc MemFree's the extra data */
/* StdCleanupFormProc MemFree's the menuitemlist and then the extra data */
extern void Nlm_StdCleanupExtraProc (Nlm_GraphiC g, Nlm_VoidPtr data);
extern void Nlm_StdCleanupFormProc (Nlm_GraphiC g, Nlm_VoidPtr data);

/*****************************************************************************
*
*   Messages are meant to be passed to dialogs and forms.  Some messages are
*   sent by the application (e.g., init, enter, reset, redraw), while most are
*   from menu items.  Numbers greater than NUM_VIB_MSG may be used by individual
*   applications, but note that NUM_VIB_MSG may change over time, so application
*   messages should be defined as NUM_VIB_MSG + ###.
*
*   Forms can handle some messages, and could be designed to pass others back to
*   the application (e.g., communicated via an app property pointer, and set in
*   the form's appmessage slot) for special servicing.
*
*****************************************************************************/

#define VIB_MSG_INIT     1
#define VIB_MSG_ENTER    2
#define VIB_MSG_RESET    3
#define VIB_MSG_REDRAW   4

#define VIB_MSG_NEW      6
#define VIB_MSG_OPEN     7
#define VIB_MSG_SAVE     8
#define VIB_MSG_SAVE_AS  9
#define VIB_MSG_CLOSE   10
#define VIB_MSG_QUIT    11
#define VIB_MSG_ACCEPT  12

#define VIB_MSG_IMPORT  15
#define VIB_MSG_EXPORT  16
#define VIB_MSG_PRINT   17
#define VIB_MSG_ABOUT   18
#define VIB_MSG_HELP    19


#define VIB_MSG_CUT     21
#define VIB_MSG_COPY    22
#define VIB_MSG_PASTE   23
#define VIB_MSG_DELETE  24
#define VIB_MSG_UNDO    25

#define VIB_MSG_CHANGE  27
#define VIB_MSG_SELECT  28

#define NUM_VIB_MSG     30

/*****************************************************************************
*
*   FormCommandItem creates a menu item whose action is to send the designated
*   message to the specified form.  It will register the item handle in the
*   form's menuitemlist, or in the Mac's globalMenuItemList.  These handles can
*   then be accessed for disable/enable or for changing the item title.
*
*****************************************************************************/

extern Nlm_IteM Nlm_FormCommandItem (Nlm_MenU m, Nlm_CharPtr title, Nlm_BaseFormPtr bfp, Nlm_Int2 mssg);
extern Nlm_IteM Nlm_FindFormMenuItem (Nlm_BaseFormPtr bfp, Nlm_Int2 mssg);
extern void Nlm_SetFormMenuItem (Nlm_BaseFormPtr bfp, Nlm_Int2 mssg, Nlm_IteM itm);
extern Nlm_Int2 Nlm_RegisterFormMenuItemName (Nlm_CharPtr title);

/*****************************************************************************
*
*   Folder tabs object is implemented as an opaque dialog.  Tool bar has a
*   slightly different appearance.  Folder tabs should be used in forms where
*   every page should be visited.  Tool bars are for viewers where different
*   formats are shown in different pages.
*
*****************************************************************************/

typedef void (*Nlm_TabActnProc) (Nlm_VoidPtr, Nlm_Int2, Nlm_Int2);

extern Nlm_DialoG Nlm_CreateFolderTabs (Nlm_GrouP h, Nlm_CharPtr PNTR titles, Nlm_Int2 initPage,
                                        Nlm_Int2 maxPerLine, Nlm_Int2 indentNextLine,
                                        Nlm_FonT font, Nlm_Int2 horizMargin, Nlm_Int2 vertMargin,
                                        Nlm_Int2 spaceBtwn, Nlm_Int2 cornerTaper, Nlm_Int2 endExt,
                                        Nlm_TabActnProc changeView, Pointer userdata);
#ifdef WIN_MSWIN
#define SYSTEM_FOLDER_TAB programFont,20,5,2,6,10
#define PROGRAM_FOLDER_TAB systemFont,10,5,2,4,5
#else
#define SYSTEM_FOLDER_TAB systemFont,20,5,2,6,10
#define PROGRAM_FOLDER_TAB programFont,10,5,2,4,5
#endif

extern void Nlm_ReleaseFolderTabButtons (Nlm_PaneL p, Nlm_PoinT pt);
extern void Nlm_SetFolderTabButton (Nlm_GraphiC a, Nlm_Int2 value, Nlm_Boolean savePort);
extern Nlm_DialoG Nlm_CreateFolderTabButtons (Nlm_GrouP h, Nlm_CharPtr PNTR titles, Nlm_Int2 initPage,
					  Nlm_Int2 maxPerLine, Nlm_Int2 indentNextLine,
					  Nlm_FonT font, Nlm_Int2 horizMargin, Nlm_Int2 vertMargin,
					  Nlm_Int2 spaceBtwn, Nlm_Int2 cornerTaper, Nlm_Int2 endExt,
					  Nlm_TabActnProc changeView, Pointer userdata);
extern void Nlm_SetFolderTabTitle (Nlm_DialoG dialog, Nlm_Int2 iTab, Nlm_CharPtr title);

extern Nlm_DialoG Nlm_CreateTextTabs (Nlm_GrouP h, Nlm_CharPtr PNTR titles, Nlm_Int2 initPage,
                                      Nlm_Int2 maxPerLine, Nlm_Int2 indentNextLine,
                                      Nlm_FonT font, Nlm_Int2 horizMargin, Nlm_Int2 vertMargin,
                                      Nlm_Int2 spaceBtwn, Nlm_TabActnProc changeView,
                                      Pointer userdata);

#define SYSTEM_TEXT_TAB systemFont,16,16,5
#define PROGRAM_TEXT_TAB programFont,14,14,5

/*****************************************************************************
*
*   A tag list is an arbitrary opaque spreadsheet dialog.
*
*****************************************************************************/

#define MAX_TAGLIST_ROWS   8
#define MAX_TAGLIST_COLS   12
#define STD_TAG_SPACING    2

#define TAGLIST_TEXT     1
#define TAGLIST_POPUP    2
#define TAGLIST_LIST     3
#define TAGLIST_PROMPT   4

typedef void (*TaglistCallback) PROTO ((Pointer userdata));

typedef struct Nlm_taglist {
  DIALOG_MESSAGE_BLOCK
  Nlm_Int2               rows;
  Nlm_Int2               cols;
  Nlm_Uint2Ptr           types;
  Nlm_Handle             control [MAX_TAGLIST_ROWS * MAX_TAGLIST_COLS];
  Nlm_EnumFieldAssocPtr  PNTR alists;
  Nlm_BaR                bar;
  Nlm_BaR                left_bar;
  Nlm_Int2               max;
  Nlm_Boolean            noExtend;
  ValNodePtr             vnp;
  TaglistCallback        PNTR callbacks;
  Pointer                callback_data;
  ButtoN                 clear_btns [MAX_TAGLIST_ROWS];
  Boolean                ask_before_clear;
} Nlm_TagList, PNTR Nlm_TagListPtr;

extern Nlm_DialoG Nlm_CreateTagListDialog (Nlm_GrouP h, Nlm_Uint2 rows, Nlm_Uint2 cols,
                                           Nlm_Int2 spacing, Nlm_Uint2Ptr types,
                                           Nlm_Uint2Ptr textWidths,
                                           Nlm_EnumFieldAssocPtr PNTR alists,
                                           Nlm_ToDialogFunc tofunc, Nlm_FromDialogFunc fromfunc);

extern Nlm_DialoG CreateTagListDialogEx (Nlm_GrouP h, Nlm_Uint2 rows, Nlm_Uint2 cols,
                                         Nlm_Int2 spacing, Nlm_Uint2Ptr types,
                                         Nlm_Uint2Ptr textWidths, Nlm_EnumFieldAssocPtr PNTR alists,
                                         Nlm_Boolean useBar, Nlm_Boolean noExtend,
                                         Nlm_ToDialogFunc tofunc, Nlm_FromDialogFunc fromfunc);

extern Nlm_DialoG CreateTagListDialogExEx (Nlm_GrouP h, Nlm_Uint2 rows, Nlm_Uint2 cols,
                                       Nlm_Int2 spacing, Nlm_Uint2Ptr types,
                                       Nlm_Uint2Ptr textWidths, Nlm_EnumFieldAssocPtr PNTR alists,
                                       Nlm_Boolean useBar, Nlm_Boolean noExtend,
                                       Nlm_ToDialogFunc tofunc, Nlm_FromDialogFunc fromfunc,
                                       TaglistCallback PNTR callbacks, Pointer callback_data,
                                       Boolean useLeftBar);

extern Nlm_DialoG CreateTagListDialogEx3 (Nlm_GrouP h, Nlm_Uint2 rows, Nlm_Uint2 cols,
                                       Nlm_Int2 spacing, Nlm_Uint2Ptr types,
                                       Nlm_Uint2Ptr textWidths, Nlm_EnumFieldAssocPtr PNTR alists,
                                       Nlm_Boolean useBar, Nlm_Boolean noExtend,
                                       Nlm_ToDialogFunc tofunc, Nlm_FromDialogFunc fromfunc,
                                       TaglistCallback PNTR callbacks, Pointer callback_data,
                                       Boolean useLeftBar, Boolean useClearBtns);

extern Nlm_CharPtr Nlm_ExtractTagListColumn (Nlm_CharPtr source, Nlm_Int2 col);
extern Nlm_CharPtr ReplaceTagListColumn (CharPtr source, CharPtr new_value, Int2 col);

#ifdef WIN_MAC
#define currentFormDataPtr Nlm_currentFormDataPtr
#define globalMenuItemList Nlm_globalMenuItemList
#define globalMenuListSize Nlm_globalMenuListSize
#endif
#define UIEnum Nlm_UIEnum
#define UIEnumPtr Nlm_UIEnumPtr
#define EnumFieldAssoc Nlm_EnumFieldAssoc
#define EnumFieldAssocPtr Nlm_EnumFieldAssocPtr
#define EnumFieldAlist Nlm_EnumFieldAlist
#define GetEnumName Nlm_GetEnumName
#define InitEnumPopup Nlm_InitEnumPopup
#define SetEnumPopup Nlm_SetEnumPopup
#define GetEnumPopup Nlm_GetEnumPopup
#define GetEnumPopupByName Nlm_GetEnumPopupByName
#define SetEnumPopupByName Nlm_SetEnumPopupByName
#define WhereInEnumPopup Nlm_WhereInEnumPopup
#define SortEnumFieldAlist Nlm_SortEnumFieldAlist
#define DuplicateEnumFieldAlist Nlm_DuplicateEnumFieldAlist
#define FreeEnumFieldAlist Nlm_FreeEnumFieldAlist
#define MakeEnumFieldAlistFromValNodeList Nlm_MakeEnumFieldAlistFromValNodeList
#define CreateEnumPopupListInitVal Nlm_CreateEnumPopupListInitVal
#define CreateEnumPopupListInitName Nlm_CreateEnumPopupListInitName
#define HandleActnProc Nlm_HandleActnProc
#define RepeatProcOnHandles Nlm_RepeatProcOnHandles
#define SafeShow Nlm_SafeShow
#define SafeHide Nlm_SafeHide
#define SafeEnable Nlm_SafeEnable
#define SafeDisable Nlm_SafeDisable
#define SafeSetValue Nlm_SafeSetValue
#define SafeSetStatus Nlm_SafeSetStatus
#define SafeSetTitle Nlm_SafeSetTitle
#define TrimSpacesAroundString Nlm_TrimSpacesAroundString
#define SaveStringFromText Nlm_SaveStringFromText
#define JustSaveStringFromText Nlm_JustSaveStringFromText
#define SetTextFromVnp Nlm_SetTextFromVnp
#define GetVnpFromText Nlm_GetVnpFromText
#define MaxStringWidths Nlm_MaxStringWidths
#define WidestString Nlm_WidestString
#define MaxAlistWidths Nlm_MaxAlistWidths
#define WidestAlist Nlm_WidestAlist
#define MultiLinePrompt Nlm_MultiLinePrompt
#define MultiLinePromptEx Nlm_MultiLinePromptEx
#define JustInvalObject Nlm_JustInvalObject
#define StringHasNoText Nlm_StringHasNoText
#define TextHasNoText Nlm_TextHasNoText
#define DialoG Nlm_DialoG
#define ToDialogFunc Nlm_ToDialogFunc
#define FromDialogFunc Nlm_FromDialogFunc
#define DialogActnFunc Nlm_DialogActnFunc
#define DialogTestFunc Nlm_DialogTestFunc
#define DialogMessageFunc Nlm_DialogMessageFunc
#define DialogInOutFunc Nlm_DialogInOutFunc
#define BaseDialogPtr Nlm_BaseDialogPtr
#define PointerToDialog Nlm_PointerToDialog
#define DialogToPointer Nlm_DialogToPointer
#define TestDialog Nlm_TestDialog
#define SendMessageToDialog Nlm_SendMessageToDialog
#define ImportDialog Nlm_ImportDialog
#define ExportDialog Nlm_ExportDialog
#define SetDialogActnProc Nlm_SetDialogActnProc
#define AlistDialogData Nlm_AlistDialogData
#define AlistDialogPtr Nlm_AlistDialogPtr
#define CreateEnumPopupDialog Nlm_CreateEnumPopupDialog
#define CreateEnumListDialog Nlm_CreateEnumListDialog
#define ForM Nlm_ForM
#define ToFormFunc Nlm_ToFormFunc
#define FromFormFunc Nlm_FromFormFunc
#define FormActnFunc Nlm_FormActnFunc
#define FormTestFunc Nlm_FormTestFunc
#define FormMessageFunc Nlm_FormMessageFunc
#define FormInOutFunc Nlm_FormInOutFunc
#define BaseFormPtr Nlm_BaseFormPtr
#define PointerToForm Nlm_PointerToForm
#define FormToPointer Nlm_FormToPointer
#define TestForm Nlm_TestForm
#define SendMessageToForm Nlm_SendMessageToForm
#define ImportForm Nlm_ImportForm
#define ExportForm Nlm_ExportForm
#define StdCloseWindowProc Nlm_StdCloseWindowProc
#define StdCancelButtonProc Nlm_StdCancelButtonProc
#define StdSendCloseWindowMessageProc Nlm_StdSendCloseWindowMessageProc
#define StdSendCancelButtonMessageProc Nlm_StdSendCancelButtonMessageProc
#define SetFormActnProc Nlm_SetFormActnProc
#define StdAcceptFormButtonProc Nlm_StdAcceptFormButtonProc
#define StdSendAcceptButtonMessageProc Nlm_StdSendAcceptButtonMessageProc
#define StdCleanupExtraProc Nlm_StdCleanupExtraProc
#define StdCleanupFormProc Nlm_StdCleanupFormProc
#define TabActnProc Nlm_TabActnProc
#define FormCommandItem Nlm_FormCommandItem
#define FindFormMenuItem Nlm_FindFormMenuItem
#define SetFormMenuItem Nlm_SetFormMenuItem
#define RegisterFormMenuItemName Nlm_RegisterFormMenuItemName
#define CreateFolderTabs Nlm_CreateFolderTabs
#define ReleaseFolderTabButtons Nlm_ReleaseFolderTabButtons
#define SetFolderTabButton Nlm_SetFolderTabButton
#define CreateFolderTabButtons Nlm_CreateFolderTabButtons
#define SetFolderTabTitle Nlm_SetFolderTabTitle
#define CreateTextTabs Nlm_CreateTextTabs
#define TagList Nlm_TagList
#define TagListPtr Nlm_TagListPtr
#define CreateTagListDialog Nlm_CreateTagListDialog
#define ExtractTagListColumn Nlm_ExtractTagListColumn

#ifdef __cplusplus
}
#endif

#endif /* ndef _DLOGUTIL_ */

