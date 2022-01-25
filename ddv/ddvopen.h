/*  $Id: ddvopen.h,v 1.12 2000/01/12 21:52:17 durand Exp $
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
* File Name:  ddvopen.h
*
* Author:  Patrick Durand
*
* Version Creation Date:   06/19/99
*
* $Revision: 1.12 $
*
* File Description: 
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: ddvopen.h,v $
* Revision 1.12  2000/01/12 21:52:17  durand
* add import function; update menus when DDV is loaded from Cn3D
*
* Revision 1.11  1999/12/30 21:08:45  lewisg
* bioseq import dialog
*
* Revision 1.10  1999/12/29 22:55:03  lewisg
* get rid of seqalign id
*
* Revision 1.9  1999/12/20 20:20:41  lewisg
* allow cn3d to do color and ddv to do case when both are running
*
* Revision 1.8  1999/12/03 23:17:23  lewisg
* Patrick's new global update msg, argument passing when launching ddv, experimental editing
*
* Revision 1.7  1999/11/29 15:26:26  durand
* designed a new GUI to fix problems under MacOS, Linux and SGI
*
* Revision 1.6  1999/11/03 21:29:47  durand
* add CTRL and SHFT keys for mouse selection. redesign the loader functions of DDV to properly register the message callbacks
*
* Revision 1.5  1999/10/29 19:04:21  durand
* move DDVUpdateMSG in objmgr.h
*
* Revision 1.4  1999/10/29 14:15:39  durand
* add simple mouse selection functions
*
* Revision 1.3  1999/10/20 13:17:18  durand
* add display for disc. SeqAlign tails
*
* Revision 1.2  1999/10/16 15:02:25  durand
* fixes due to toolkit build failed
*
* Revision 1.1  1999/09/30 14:10:28  durand
* add ddv to toolkit
*
* Revision 1.8  1999/09/28 17:29:41  durand
* add include of udviewer.h
*
* Revision 1.7  1999/09/16 13:07:53  durand
* add File|Close and File|Open|Network commands
*
* Revision 1.6  1999/09/09 21:55:06  durand
* instantiate the Fle|Close command of DDV
*
* Revision 1.5  1999/07/29 12:43:07  durand
* update DDV_GetAndCheckSeqAlign
*
* Revision 1.4  1999/07/01 15:28:29  durand
* validate function loaders of DDV
*
* Revision 1.2  1999/06/28 22:07:21  durand
* add loader functions and clean the code with Lint and Purify
*
* Revision 1.1  1999/06/19 17:21:07  durand
* add Vibrant DDV code
*
*
*
* ==========================================================================
*/

#ifndef _DDVOPEN_
#define _DDVOPEN_

#ifdef __cplusplus
extern "C" {
#endif


#include <ncbi.h>
#include <udviewer.h>
#include <samutil.h>
#include <accentr.h>

/******************************************************************************

	ERROR / Information messages from DDV_OPEN module 

******************************************************************************/
#define DVV_MSG_O_I_NOMAINDATA   2000
#define DVV_MSG_O_E_READFILE     2001
#define DVV_MSG_O_E_READGI       2002
#define DVV_MSG_O_E_NOFETCHFUNC  2003
#define DVV_MSG_O_E_BADTYPE      2004
#define DVV_MSG_O_E_NOTHINGTODO  2005
#define DVV_MSG_O_E_OPENFILEFAIL 2006

/******************************************************************************

	other defines

******************************************************************************/
	/*use only by the standalone version of DDV */
#define REG_DDV_AUTO_EDIT ObjMgrProcLoad(OMPROC_EDIT, \
		"DDV", "MSA_Editor", OBJ_SEQALIGN, 0, OBJ_SEQALIGN, 0, \
		NULL, DDV_ObjRegMasterDDV, 0)	

	/*slave 1 : the editor in a separate window*/
#define REG_DDV_SLA_EDIT ObjMgrProcLoad(OMPROC_EDIT, \
		"DDV", "MSA_Editor", OBJ_SEQALIGN, 0, OBJ_SEQALIGN, 0, \
		NULL, DDV_ObjRegSlaveEditDDV, 0)	

	/*slave 2 : the viewer in a separate window*/
#define REG_DDV_SLA_VIEW ObjMgrProcLoad(OMPROC_VIEW, \
		"DDV", "MSA_Viewer", OBJ_SEQALIGN, 0, OBJ_SEQALIGN, 0, \
		NULL, DDV_ObjRegSlaveViewDDV, 0)	

/*the following are used to delete data when closing DDV*/
#define DDV_OPENTYPE_NOTRESP ((Uint1)1)/*if DDV is not responsible to delete
                                         data when closing*/
#define DDV_OPENTYPE_FILE    ((Uint1)2)/*if ReadAsnFastaOrFlatFile was used*/
#define DDV_OPENTYPE_SEP     ((Uint1)3)/*if get a SEP*/
#define DDV_OPENTYPE_GI      ((Uint1)4)/*if fetch GI from DB was used*/

/******************************************************************************

	Data structures

******************************************************************************/
typedef struct dlgfileopendata {/*use to manage the FileOpen dialog box*/
	WindoW 		parent;			/*main window of the application*/
	TexT		FNameEditCtrl;	/*handle of the file text control*/
	ButtoN		Ok;				/*handle of the Ok button*/
	GrouP		ReadMode;		/*handle of the file type control*/
} DlgFileOpenData, PNTR DlgFileOpenDataPtr;

typedef struct ddvopendata {
	Uint1           choice;
	SeqEntryPtr     sep;
	ValNodePtr      vnp;
	} DdvOpenData, PNTR DdvOpenDataPtr;

typedef struct ddvupdatelayoutdata{
	/*color display ?*/
	Boolean bUseColors;
	/*styles for a disc. seqalign*/
	Boolean ShowLeftTail;
	Boolean ShowRightTail;
	Uint1   DispDiscStyle;
	Uint1   SpacerSize;
    Uint1   DiscJustification;
	/*style for sequences*/
	Int4    nSeq;	/*number of rows (sequences) to update*/
	Int4Ptr SeqList; /*list of row number; one-based array*/
	Uint1   RulerStyle;
	/*DDV panel handle*/
	PaneL   ddv_panel;
	}DDVUpdateLayoutData, PNTR DDVUpdateLayoutDataPtr;

/* for the BLAST import dialog */
typedef struct _DDV_ImportDialog {
    DocType AAorNN; /* is this of TYP_AA or TYP_NA */
    WindoW DDV_wImport; /*the import dialog*/
    GrouP DDV_gAccType; /*the type of accession*/
    ButtoN DDV_bImportAccept; /*accept button*/
    TexT DDV_tAccession; /* the accession */
    SeqId *sip;  /* the master sequence */
    SeqAlign *sap; /* the seqalign to add to */
    void *userdata; /* for the update message */
} DDV_ImportDialog;


/******************************************************************************

	Extern functions

******************************************************************************/
extern Int2 LIBCALLBACK DDV_ObjRegMasterDDV (Pointer data);
extern Int2 LIBCALLBACK DDV_ObjRegSlaveEditDDV (Pointer data);
extern Int2 LIBCALLBACK DDV_ObjRegSlaveViewDDV (Pointer data);
extern void DDV_OpenFile(IteM i);
extern void DDV_OpenNetwork(IteM i);
extern ValNodePtr DDV_GetAndCheckSeqAlign(FILE *fp,Int4 gi,SeqEntryPtr sep2,
	UdvFetchSeqEntryProc fetchSepProc,DdvOpenDataPtr dodp,Uint2Ptr entityID);
extern void DDV_LaunchAlignEditor (Uint2 entityID,SeqAlignPtr sap);
extern void DDV_LaunchAlignViewer (Uint2 entityID,SeqAlignPtr sap);
extern WindoW DDV_StartMainWin_Slave(SAM_ViewGlobal *vgp);

extern ValNodePtr DDV_GetSelectedRegions(SelStructPtr om_slp, Uint2 bsp_eID,
	Uint2 bsp_iID);
extern Boolean DDV_IsLetterSelected(ValNodePtr vnp_bsp, Int4 bsp_pos);

NLM_EXTERN Int4 DDV_Accession2Gi (CharPtr string, DocType type);
NLM_EXTERN void DDV_ImportBioseqDlg(DDV_ImportDialog *idp);
NLM_EXTERN void DDV_ImportBioseq(IteM i);
NLM_EXTERN SeqAlign *DDV_Blast2Seqs(Bioseq *bsp1, Bioseq *bsp2, Boolean gapped,
                         Char *progname);
extern void DDV_ImportNucSeqAlign(IteM i);
extern void DDV_ImportProtSeqAlign(IteM i);


#ifdef __cplusplus
}
#endif

#endif /* ndef _DDVOPEN_ */

