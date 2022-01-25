/*   cn3dopen.c
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
* File Name:  cn3dopen.c
*
* Author:  Christopher Hogue, Yanli Wang, Lewis Geer
*
* First Version Creation Date:   1/31/96
*
* $Revision: 6.70 $
*
* File Description: Cn3d file opening routines 
*                   
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
* $Log: cn3dopen.c,v $
* Revision 6.70  2000/01/21 15:59:05  lewisg
* add check for binary/ascii files
*
* Revision 6.69  2000/01/14 21:40:41  lewisg
* add translucent spheres, ion labels, new cpk, fix misc bugs
*
* Revision 6.68  2000/01/04 15:55:51  lewisg
* don't hang on disconnected network and fix memory leak/hang at exit
*
* Revision 6.67  1999/12/28 23:06:35  lewisg
* udv/cn3d communication
*
* Revision 6.66  1999/12/28 15:55:21  lewisg
* remove remaining mediainfo code
*
* Revision 6.65  1999/12/15 23:17:47  lewisg
* make cn3d launch udv, fix launching of non-autonomous udv, add mouseup message
*
* Revision 6.64  1999/12/13 23:20:44  lewisg
* bug fixes: duplicate color structures eliminated, clear color doesn't clear case
*
* Revision 6.63  1999/12/11 01:30:35  lewisg
* fix bugs with sharing colors between ddv and cn3d
*
* Revision 6.62  1999/12/01 16:15:54  lewisg
* interim checkin to fix blocking memory leak
*
* Revision 6.61  1999/11/24 15:23:19  lewisg
* added color selection dialogs for SS
*
* Revision 6.60  1999/11/15 18:30:08  lewisg
* get rid of extra redraws when selecting
*
* Revision 6.59  1999/11/10 23:19:41  lewisg
* rewrite of selection code for ddv
*
* Revision 6.58  1999/11/02 23:06:07  lewisg
* fix cn3d to launch correctly if there is no seqentry associated with bioseq
*
* Revision 6.57  1999/10/29 14:15:30  thiessen
* ran all Cn3D source through GNU Indent to prettify
*
* Revision 6.56  1999/10/18 16:01:09  lewisg
* fix bssp to bsssp
*
* Revision 6.55  1999/10/18 15:32:50  lewisg
* move ClearSequences() to cn3dshim.c
*
* Revision 6.54  1999/10/15 20:56:39  lewisg
* append DDV_ColorGlobal as userdata.  free memory when cn3d terminates.
*
* Revision 6.53  1999/10/05 23:18:24  lewisg
* add ddv and udv to cn3d with memory management
*
* Revision 6.52  1999/09/16 17:16:21  ywang
* open multiple salsa window for data with multiple seq-annot data
*
* Revision 6.51  1999/08/04 21:18:01  lewisg
* modularized open operations to allow sequin to launch cn3d
*
* Revision 6.50  1999/07/12 19:37:44  ywang
* fix bug
*
* Revision 6.49  1999/07/07 20:45:37  ywang
* clear domaindata, mediadata, special feature before reading in new data in cn3d
*
* Revision 6.48  1999/07/07 16:41:30  ywang
* reset Num_Bioseq as 0 when sequences are freed
*
* Revision 6.47  1999/07/06 18:30:45  chappey
* Send FLUSH message on each Bioseq to close windows
*
* Revision 6.46  1999/07/06 17:18:38  ywang
* send OM_MSG_FLUSH message to salsa to close all its windows before new sequence window launched
*
* Revision 6.45  1999/07/02 20:58:37  ywang
*  return NULL sep upon NULL pmsdThis in Cn3DFetchSeqEntry
*
* Revision 6.44  1999/07/02 19:57:43  ywang
* go through every node of SeqEntryPtr for free
*
* Revision 6.43  1999/07/02 14:39:14  ywang
* clean seq-entry fetching and OM register for netopen and local read(non mimetype) functioning
*
* Revision 6.42  1999/07/01 22:14:58  ywang
* back to previous version to use EntrezSeqEntryGet to get sequences for additional loading, before this free existant sequences
*
* Revision 6.41  1999/07/01 22:08:08  lewisg
* add resize before all redraws
*
* Revision 6.40  1999/07/01 21:44:30  ywang
* work around EntrezSeqEntryGet core dump on redunant sequence loading and free existant sequences before additional loading
*
* Revision 6.39  1999/07/01 14:07:25  ywang
* cn3dwin.c
*
* Revision 6.38  1999/06/15 19:08:17  kans
* instantiate Cn3D_useEntrez in library
*
* Revision 6.37  1999/06/15 17:57:54  ywang
* rename useEntrez as Cn3D_useEntrez
*
* Revision 6.36  1999/06/14 17:40:26  ywang
* assign mimetype to structure data read in locally or via network with respect to availability of sequences
*
* Revision 6.35  1999/04/06 20:14:09  lewisg
* more opengl
*
* Revision 6.34  1999/03/30 22:36:19  ywang
* add functions to color salsa for NcbiMimeAsn1_strucseqs & code reorganization
*
* Revision 6.33  1999/03/18 22:28:58  ywang
* add functions for saveout+readin+index user defined features
*
* Revision 6.32  1999/03/03 23:17:22  lewisg
* one master struct at a time, list slaves in structure info, bug fixes
*
* Revision 6.31  1999/02/25 23:14:31  ywang
* move around menu item and callback function, change menu item names
*
* Revision 6.30  1999/02/24 23:00:49  ywang
* record mime type at MSD node
*
* Revision 6.29  1999/01/20 18:21:20  ywang
* include salmedia.h due to the move around of MediaInfo from cn3dmsg.h to the new created salmedia.h
*
* Revision 6.28  1999/01/19 23:42:48  ywang
* fix bugs over improving color msg
*
* Revision 6.27  1999/01/14 19:07:17  kans
* network availability is configurable
*
* Revision 6.26  1998/12/22 15:40:32  ywang
* restore sequences pointer for strucseq
*
* Revision 6.25  1998/12/21  18:47:37  addess
* fixed strucseq bug found by Aron
*
* Revision 6.24  1998/12/16  17:02:57  ywang
* pick up strucseqs mime type
*
* Revision 6.23  1998/10/21  15:49:21  ywang
* attach the whole vast alignment data to master structure
*
* Revision 6.22  1998/09/23  18:38:50  ywang
* add functions to control display on domain level
*
* Revision 6.21  1998/06/30  23:29:22  ywang
* fix bugs regarding to read in more structures
*
* Revision 6.20  1998/06/30  14:48:04  ywang
* launch salsa automatically and accordingly when more structures are readin
*
* Revision 6.19  1998/06/29  19:28:00  lewisg
* on the fly update of conservation color
*
* Revision 6.17  1998/06/15 14:26:06  ywang
* automatic launch salsa when mime data got in either through command line or via local reading in
*
* Revision 6.16  1998/06/12  21:21:56  ywang
* change sCompare to make cn3d accept mime strucseq from local readin
*
* Revision 6.15  1998/06/04  16:48:36  ywang
* fix bug triggered by automatic salsa launch
*
* Revision 6.14  1998/05/27  22:15:16  ywang
* add Cn3dObjRegister() to several places
*
* Revision 6.13  1998/05/18  22:09:14  ywang
* move codes around
*
* Revision 6.12  1998/05/12  21:47:05  lewisg
* stricter conservation coloring
*
* Revision 6.11  1998/05/12 16:34:01  ywang
* take strucseq mime data
*
* Revision 6.10  1998/05/06  23:50:24  lewisg
* fixed launching problem with sequin
*
* Revision 6.9  1998/05/06 14:20:09  lewisg
* get rid of NULL structure bugs
*
* Revision 6.8  1998/04/28 22:47:21  lewisg
* master/slave color in sync
*
* Revision 6.7  1998/04/28 19:40:01  lewisg
* codewarrior fixes
*
* Revision 6.6  1998/04/28 19:38:39  lewisg
* codewarrior fixes
*
* Revision 6.5  1998/04/28 18:53:25  ywang
* take NcbiMimeAsn1_alignseq to view alignment only
*
* Revision 6.4  1998/04/28  15:20:37  kans
* cast needed on FileGets for CodeWarrior
*
* Revision 6.3  1998/04/28 15:14:29  lewisg
* moved OpenMimeFileWithDeletion to cn3dopen
*
* Revision 6.2  1998/04/27 23:23:06  lewisg
* added ability to open mime files
*
* Revision 6.1  1998/04/16 00:32:25  lewisg
* corrected neighbor mode bugs
*
* Revision 6.0  1997/08/25 18:13:38  madden
* Revision changed to 6.0
*
* Revision 5.0  1996/05/28 14:05:44  ostell
* Set to revision 5.0
*
* Revision 1.5  1996/04/18  17:01:10  hogue
* Fixed calls to FetchBS for later integration into Entrez, added feature to turn off redundant backbone model for faster spinning...
*
* Revision 1.4  1996/03/30  23:41:01  hogue
* Redraw now saves camera
*
* Revision 1.3  1996/03/29  20:01:24  hogue
* Added call to reset active structure upon opening.
*
* Revision 1.2  1996/02/02  19:40:19  hogue
* Initial Revision
*
*
* ==========================================================================
*/

#include <vibrant.h>
#include <math.h>
#include <mmdbapi.h>
#include <accentr.h>
#include <objalign.h>
#include <objseq.h>
#include <objmgr.h>
#include <sequtil.h>
#include <saledit.h>
#include <lsqfetch.h>
#include <cn3dopen.h>
#include <cn3dmain.h>
#include <algorend.h>
#include <cn3dmsg.h>
#include <salmedia.h>
#include <cn3dmodl.h>
#include <cn3dshim.h>



static Boolean Cn3D_Open_InUse = FALSE;
static WindoW Cn3D_wNetOpen;
static TexT Cn3D_tOpen;
static GrouP Cn3D_gAccType;
static ButtoN Cn3D_bOpenAccept;
static GrouP Cn3D_gMdlLvl;
static GrouP Cn3D_gMdlNo;

static WindoW Cn3D_wOpen;
static ButtoN Cn3D_bOpenBrowse;

Boolean Mime_ReadIn = FALSE;

Int2 Cn3DMime = 0;

/*------------------------------------------------------*/
SeqEntryPtr Cn3DFetchSeqEntry(PMSD pmsdThis)
{
    PDNMM pdnmmHead = NULL;
    PMMD pmmdThis = NULL;
    SeqId si;
    Bioseq *bsp;

    if (pmsdThis == NULL)
        return NULL;

    pdnmmHead = pmsdThis->pdnmmHead;
    while (pdnmmHead) {
        pmmdThis = pdnmmHead->data.ptrvalue;
        if (pmmdThis == NULL)
            goto errot;
        if (pmmdThis->bWhat != (Byte) AM_PROT
            && pmmdThis->bWhat != (Byte) AM_RNA
            && pmmdThis->bWhat != (Byte) AM_DNA)
            goto errot;
        if (pmmdThis->iGi != 0) {
            si.choice = SEQID_GI;
            si.data.intvalue = pmmdThis->iGi;
            bsp = BioseqLockById(&si);
            BioseqUnlockById(&si);

            if (bsp != NULL)
                return SeqEntryFind(&si);
        }

      errot:
        pdnmmHead = pdnmmHead->next;
    }

    return NULL;
}

/*------------------------------------------------------*/
void fnMarkAlignedResiduesForStrucSeqs(PDNMS pdnmsThis)
{
    PMSD pmsdThis = NULL;
    PMMD pmmdThis = NULL;
    PMGD pmgdThis = NULL;

    SeqIdPtr sip = NULL;
    SeqAnnotPtr sap = NULL;
    SeqAlignPtr salp = NULL, salp_curr = NULL;
    DenseSegPtr dssp = NULL;
    Int4Ptr starts = NULL;
    Int4Ptr lens = NULL;

    Int4 from = 0, to = 0, nres = 0;
    Int4 numseg = 0;


    pmsdThis = pdnmsThis->data.ptrvalue;
    if (pmsdThis == NULL)
        return;

    sap = Cn3D_ColorData.sap;
    if (sap == NULL)
        return;

    while (sap) {
        if (sap->type == 2) {
            salp = sap->data;
            break;
        }
        sap = sap->next;
    }

    if (salp == NULL)
        return;
    salp_curr = salp;
    if (salp_curr) {
        dssp = salp->segs;
        sip = dssp->ids;
        pmmdThis = GetMMFromMSDBySeqId(pmsdThis, sip);
        /* assume the first sip is for the sequence with known structure */
    }

    if (pmmdThis == NULL)
        return;

    salp_curr = salp;
    while (salp_curr) {
        pmsdThis->bAligned++;
        salp_curr = salp_curr->next;
    }

    salp_curr = salp;
    while (salp) {
        dssp = salp->segs;
        starts = dssp->starts;
        lens = dssp->lens;

        for (numseg = 0; numseg < dssp->numseg; numseg++, lens++) {
            from = *starts;
            to = from + *lens;
            if (from == -1) {
                starts++;
                starts++;
                continue;
            }
            starts++;
            if (*starts != -1) {
                for (nres = from; nres < to; nres++) {
                    pmgdThis = GetMGFromMM(pmmdThis, nres + 1);
                    if (pmgdThis)
                        pmgdThis->bReserved++;
                }
            }
            starts++;
        }

        salp = salp->next;
    }

}

/*----------------------------------------------------*/
/**********
fnMarkAlignedResidue()
Given a master and model, mark all of the residues in the slave and master that are aligned.
**********/

ValNodePtr fnMarkAlignedResidues(PDNMS pdnmsMaster, PDNMS pdnmsSlave,
                                 BiostrucFeaturePtr pbsfThis)
{

    ValNodePtr pvnAlignment;
    ValNodePtr pvnThis = NULL;
    ChemGraphPntrsPtr pcgpThis;
    ValNodePtr pvnListMaster = NULL, pvnListSlave = NULL;
    PFB pfbMaster = NULL, pfbSlave = NULL;
    ChemGraphAlignmentPtr pcgaSlave;


    /* find the corresponding alignment */
    pvnAlignment =
        ValNodeFindNext(pbsfThis->Location_location, NULL,
                        Location_location_alignment);
    if (pvnAlignment == NULL)
        return NULL;


    pcgaSlave = pvnAlignment->data.ptrvalue;
    pvnThis = pcgaSlave->alignment;
    if (pvnThis) {
        pcgpThis = (ChemGraphPntrsPtr) pvnThis;
        pvnListMaster = MakeChemGraphNodeList(pdnmsMaster, pcgpThis);
        pcgpThis = (ChemGraphPntrsPtr) pvnThis->next;
        pvnListSlave = MakeChemGraphNodeList(pdnmsSlave, pcgpThis);
        if (!pvnListMaster || !pvnListSlave)
            return NULL;
        while (pvnListMaster && pvnListSlave) {
            pfbMaster = (PFB) pvnListMaster->data.ptrvalue;
            pfbMaster->bReserved++;
            pfbSlave = (PFB) pvnListSlave->data.ptrvalue;
            pfbSlave->bReserved++;
            ((PMGD) pfbSlave)->pbMasterReserved = &(pfbMaster->bReserved); /* use unused feature pointer to point slave back at master */
            pvnListSlave = pvnListSlave->next;
            pvnListMaster = pvnListMaster->next;
        }                       /* while pvnListMaster */
        ValNodeFree(pvnListMaster);
        ValNodeFree(pvnListSlave);
    }                           /* while pvnThis */
    return pvnAlignment;
}

/************
fnClearMarkedResidues
callback function used to clear all of the alignment pointers and counters used for protein conservation calculations
************/

void LIBCALLBACK fnClearMarkedResidues(PFB pfbThis, Int4 iModel,
                                       Int4 iIndex, Pointer ptr)
{
    PMGD pmgdThis;

    pfbThis->bReserved = 0;
    pmgdThis = (PMGD) pfbThis;
    pmgdThis->pbMasterReserved = NULL;
}

/******************************************************************************
*
* Traverses the given MSD and sets it up for displaying one model or animation
* Turns off backbone model if all atoms is loaded in?
*
******************************************************************************/

static void MMDB_OpenTraverse(PMSD pmsd)
{
    PDNML pdnmlThis = NULL;
    PMLD pmldThis = NULL;
    PMLD pmldOne = NULL;
    PMLD pmldAll = NULL;

    pdnmlThis = pmsd->pdnmlModels;
    /* set up for doing one model or animation */
    while (pdnmlThis) {
        pmldThis = (PMLD) pdnmlThis->data.ptrvalue;
        if (pmldThis->iType == Model_type_ncbi_backbone)
            pmldOne = pmldThis;
        if (pmldThis->iType == Model_type_ncbi_all_atom)
            pmldAll = pmldThis;
        pdnmlThis = pdnmlThis->next;
    }

    if (pmldOne && pmldAll)
        pmldOne->bSelected &= (Byte) 0xFE;
}

/******************************************************************************
*
* Sets up Cn3D for reading in a new structure and any potential alignments
* and bioseqs.
*
******************************************************************************/

NLM_EXTERN void Cn3D_OpenStart()
{
    Int4 iBioseq = 0;

#ifndef _OPENGL
    Cn3D_SaveActiveCam();
#endif

    ObjMgrDeSelectAll();
    Mime_ReadIn = FALSE;

    ClearSequences();
    ClearStructures();
    ClearRest();
}

static void Cn3D_StartNet()
{
    if (Cn3D_ColorData.EntrezOn) return;
    if (Cn3D_ColorData.UseEntrez) {
        if(EntrezBioseqFetchEnable("Cn3D", FALSE)) {
            if(EntrezInit("Cn3D", TRUE, NULL)) Cn3D_ColorData.EntrezOn = TRUE;
        }
    }
}

/******************************************************************************
*
* Finishes setting up Cn3D. Called after reading in a new structure and any
* alignments and bioseqs.  If necessary, will load in a sequence via
* netentrez if a network connection is open.
*
******************************************************************************/

/* call back to replace gi's with full seqid's */
static void LIBCALLBACK fnRewireSips(PFB pfbThis, Int4 iModel,
                                       Int4 iIndex, Pointer ptr)
{
    PMMD pmmdThis;
    Bioseq *bsp;
    SeqId *sip;
    
    pmmdThis = (PMMD) pfbThis;
    if(pmmdThis->pSeqId == NULL) return;
    bsp = BioseqLockById(pmmdThis->pSeqId);
    if(bsp == NULL) return;
    sip = SeqIdDupList(bsp->id);
    BioseqUnlockById(pmmdThis->pSeqId);
    if(sip == NULL) {
        BioseqUnlockById(pmmdThis->pSeqId);
        return;
    }
    SeqIdFree(pmmdThis->pSeqId);
    pmmdThis->pSeqId = sip;
    return;
}


NLM_EXTERN void Cn3D_OpenEnd()
{
    PMSD pmsdMaster = NULL;
    PDNMS pdnmsMaster = NULL, pdnmsSlave = NULL;
    SeqEntryPtr sep = NULL;

    Cn3DIndexUserDefinedFeature();

#ifdef _OPENGL
    OGL_Reset(Cn3D_ColorData.OGL_Data);
#else
    ZoomAll3D(Cn3D_v3d);
#endif


        
    pdnmsMaster = GetSelectedModelstruc();
    if(pdnmsMaster == NULL) return;
    pmsdMaster = pdnmsMaster->data.ptrvalue;
    if(pmsdMaster == NULL) return;

    if (!Mime_ReadIn || pmsdMaster->iMimeType == NcbiMimeAsn1_entrez) {
        Cn3D_StartNet();
        if (Cn3D_ColorData.EntrezOn) {
            sep = Cn3DFetchSeqEntry(pmsdMaster);
            if (sep != NULL) {
                SeqAnnot *sap;
                pmsdMaster->iMimeType = NcbiMimeAsn1_strucseq;
                Cn3D_RegisterSeqEntry(sep);
            }
            else pmsdMaster->iMimeType = NcbiMimeAsn1_entrez;
        }
        else pmsdMaster->iMimeType = NcbiMimeAsn1_entrez;
    }

    Cn3D_ColorData.IsUserData = TRUE;
    if(Cn3D_ColorData.pDDVColorGlobal == NULL) {
        Cn3D_ColorData.pDDVColorGlobal = DDV_CreateColorGlobal(FALSE);
        DDV_LoadSSColor(Cn3D_ColorData.pDDVColorGlobal, "CN3D");
        Cn3D_ColorData.IsUserData = FALSE;
        Cn3D_RegisterColor();
    }
    
    if (pmsdMaster->iMimeType != NcbiMimeAsn1_entrez) {    
        TraverseMolecules(pdnmsMaster, 0, 0, NULL,
            (pNodeFunc) fnRewireSips);
        for(pdnmsSlave = pmsdMaster->pdnmsSlaves; pdnmsSlave != NULL;
        pdnmsSlave = pdnmsSlave->next)
            TraverseMolecules(pdnmsSlave, 0, 0, NULL,
            (pNodeFunc) fnRewireSips);
    }

    Cn3D_ResetActiveStrucProc();
    if (pmsdMaster->iMimeType != NcbiMimeAsn1_entrez) {    
        
        LaunchSequenceWindow();
    }
    
}

/******************************************************************************
*
* Opens the file with filename and tries to read in NcbiMimeAsn1.  Doesn't 
* care if it binary or text encoded.  If removeIt is set, deletes the file
* afterwards.  The NcbiMimeAsn1 is read into the object manager and MMDBapi
*
******************************************************************************/

Boolean OpenMimeFileWithDeletion(CharPtr filename, Boolean removeIt)
{
    Boolean retval = FALSE;
    AsnIoPtr aip;
    Char buf[50];
    NcbiMimeAsn1Ptr mime;
    FILE *fp = FileOpen(filename, "r");

    WatchCursor();
    if (!fp)
        return FALSE;

    FileRead(buf, 1, StrLen(PRINT_FORM_MIME_NAME), fp);
    FileClose(fp);
    if (StrNCmp(buf, PRINT_FORM_MIME_NAME, StrLen(PRINT_FORM_MIME_NAME)) ==
        0) aip = AsnIoOpen(filename, "r");
    else
        aip = AsnIoOpen(filename, "rb");
    if (!aip)
        return FALSE;
    mime = NcbiMimeAsn1AsnRead(aip, NULL);
    AsnIoClose(aip);
    if (!mime)
        return FALSE;
    Cn3DMime = mime->choice;

    Cn3D_OpenStart();
    Mime_ReadIn = TRUE;
    retval = MMDB_ReadMime(mime);
    Cn3D_OpenEnd();

    if (removeIt) {
        FileRemove(filename);
    }
    ArrowCursor();

    return retval;
}

/******************************************************************************
*
* Meant to be an ncbiobj and mmdb dependent function to process mime input
* into objects and MSD's.  Still has a few dependencies to be eliminated
*
******************************************************************************/

Boolean MMDB_ReadMime(NcbiMimeAsn1Ptr mime)
{
    Boolean retval = FALSE;
    BiostrucFeaturePtr pbsfThis;
    PMSD pmsdSlave = NULL, pmsdMaster = NULL;
    PDNMS pdnmsMaster = NULL, pdnmsSlave = NULL;
    BiostrucAlignSeqPtr pbsasThis = NULL;
    SeqAnnotPtr sap = NULL;
    SeqAlignPtr salp = NULL;


    do { {                      /* TRY */
            EntrezGeneralPtr egp;
            BiostrucAlignPtr pbsaThis;
            BiostrucSeqPtr bssp;
            BiostrucSeqsPtr bsssp;
            PDNMS pdnms = NULL;
            PMSD pmsdThis = NULL;
            PDNML pdnmlThis = NULL;
            ValNodePtr pvnAlignment;
            PDNTRN pdnTransform = NULL;
            PARS parsThis;
            BiostrucPtr pbsThis;
            SeqAnnot *sap;

            switch (mime->choice) {

            case NcbiMimeAsn1_entrez:
                egp = (EntrezGeneralPtr) mime->data.ptrvalue;
                if (!egp || !egp->Data_data)
                    break;      /* THROW */

                retval = TRUE;
                switch (egp->Data_data->choice) {
                case Data_data_structure:
                    pdnms =
                        MakeAModelstruc((BiostrucPtr) egp->Data_data->data.
                                        ptrvalue);
                    if (!pdnms)
                        break;

                    pmsdThis = (PMSD) pdnms->data.ptrvalue;

                    pmsdThis->iMimeType = NcbiMimeAsn1_entrez;

                    MMDB_OpenTraverse(pmsdThis);
                    break;

                case Data_data_nuc:
                case Data_data_prot:
                case Data_data_ml:
                case Data_data_genome:
                default:
                    break;
                }
                break;
            case NcbiMimeAsn1_strucseq:
                bssp = (BiostrucSeqPtr) mime->data.ptrvalue;
                pdnms = MakeAModelstruc((BiostrucPtr) bssp->structure);
                if (!pdnms)
                    break;

                pmsdThis = (PMSD) pdnms->data.ptrvalue;
                MMDB_OpenTraverse(pmsdThis);

                Cn3D_RegisterSeqEntry(bssp->sequences);

                break;

            case NcbiMimeAsn1_strucseqs:
                bsssp = (BiostrucSeqsPtr) mime->data.ptrvalue;
                pdnms = MakeAModelstruc((BiostrucPtr) bsssp->structure);
                if (!pdnms)
                    break;

                pmsdThis = (PMSD) pdnms->data.ptrvalue;

                pmsdThis->iMimeType = NcbiMimeAsn1_strucseqs;

                MMDB_OpenTraverse(pmsdThis);
                fnMarkAlignedResiduesForStrucSeqs(pdnms);

                Cn3D_RegisterSeqEntry(bsssp->sequences);
                Cn3D_RegisterSeqAnnot(bsssp->seqalign);


                parsThis = NewStrucSeqsRenderSet();
                pmsdThis->pExtra = parsThis;

                break;

            case NcbiMimeAsn1_alignseq:
                pbsasThis = (BiostrucAlignSeqPtr) mime->data.ptrvalue;
                sap = pbsasThis->seqalign;
                Cn3D_RegisterSeqAnnot(pbsasThis->seqalign);
                Cn3D_RegisterSeqEntry(pbsasThis->sequences);

                return retval;
            case NcbiMimeAsn1_alignstruc: /* this is the code that received alignments */
                pbsaThis = (BiostrucAlignPtr) mime->data.ptrvalue;
                if (!pbsaThis)
                    break;      /* THROW */
                retval = TRUE;
                /* load in the master */
                pdnmsMaster = MakeAModelstruc((BiostrucPtr) pbsaThis->master); /* grab the master struc */
                if (!pdnmsMaster)
                    break;

                pmsdMaster = (PMSD) pdnmsMaster->data.ptrvalue;

                pmsdMaster->iMimeType = NcbiMimeAsn1_alignstruc;

                parsThis = NewAlignRenderSet();
                pmsdMaster->pExtra = parsThis;
                pmsdMaster->bVisible = TRUE;
                MMDB_OpenTraverse(pmsdMaster);

#ifdef CN3DDEBUG
                {
                    AsnIoPtr aip;
                    
                    aip = AsnIoOpen("align.out", "w");
                    SeqAnnotAsnWrite(pbsaThis->seqalign, aip, NULL);
                    AsnIoClose(aip);
                }
#endif /* CN3DDEBUG */
                

                /* add the alignment seq annot ptr, etc. */
                pmsdMaster->psaStrucAlignment = pbsaThis->alignments;
                Cn3D_RegisterSeqAnnot(pbsaThis->seqalign);
                Cn3D_RegisterSeqEntry(pbsaThis->sequences);


                SetNeighborOn(); /* turn on neighbor mode */

                SetMasterModelstruc(pdnmsMaster);

                /* do a slave */

                pbsThis = (BiostrucPtr) pbsaThis->slaves;
                pbsfThis =
                    (BiostrucFeaturePtr) pbsaThis->alignments->features->
                    features;
                pvnAlignment = NULL;
/* pmsdMaster->pdnsfFeatures = DValNodeAddPointer(NULL, 0, pbsfThis); *//* tack the alignments onto the master model struct. */
                TraverseGraphs(pdnmsMaster, 0, 0, NULL,
                               (pNodeFunc) fnClearMarkedResidues);

                while (pbsThis) {
                    pdnmsSlave = MakeAModelstruc(pbsThis);
                    if (!pdnmsSlave)
                        break;
                    TraverseGraphs(pdnmsSlave, 0, 0, NULL,
                                   (pNodeFunc) fnClearMarkedResidues);
                    pmsdSlave = (PMSD) pdnmsSlave->data.ptrvalue;
                    pmsdSlave->bMaster = FALSE; /* this is not a master struct */
                    pmsdMaster->bAligned++;
                    pmsdSlave->pbAligned = &(pmsdMaster->bAligned);
                    MMDB_OpenTraverse(pmsdSlave);
                    pmsdSlave->pExtra = parsThis;
                    pmsdSlave->bVisible = TRUE; /* turn them all on by default */

                    pvnAlignment =
                        fnMarkAlignedResidues(pdnmsMaster, pdnmsSlave,
                                              pbsfThis);
                    if (!pvnAlignment)
                        break;

                    /* use this to clear out the bReserved value. TraverseGraphs( DValNodePtr pdnModel, Int4 iModel, Int4 iIndex,
                       Pointer ptr, pNodeFunc pfnCallMe)
                       (*pfnCallMe)((PFB) pxxxxThis,  Int4 iMode, Int4 iIndex, Pointer ptr)
                     */

                    pdnTransform = NULL;
                    /* create the spatial transformation */
                    TransformToDNTRN(&pdnTransform,
                                     ((ChemGraphAlignmentPtr)
                                      pvnAlignment->data.ptrvalue)->
                                     transform);
                    if (pdnTransform == NULL)
                        break;
                    /* loop over the slave's models with the transformation */
                    pdnmlThis = pmsdSlave->pdnmlModels;

                    while (pdnmlThis) {
                        TraverseAtoms(pdnmsSlave, pdnmlThis->choice, 0,
                                      pdnTransform, DoApplyTransform);
                        TraverseSolids(pdnmsSlave, pdnmlThis->choice, 0,
                                       pdnTransform, DoApplyTransform);
                        pdnmlThis = pdnmlThis->next;
                    }

                    FreeDNTRN(pdnTransform);
                    pbsThis = pbsThis->next;
                    pbsfThis = pbsfThis->next;
                    /* turn the transform into a dvalnode for the traverse */
                }               /*while pbsThis */

                break;
            default:
                break;
            }

    }
    } while (0);                /* End-of-TRY-block */

    return retval;
}


static void Cn3D_OpenEnableProc(TexT t)
{
    Char str[32];
    GetTitle(Cn3D_tOpen, str, sizeof(str));
    if (StringLen(str) == 0) {
        Disable(Cn3D_bOpenAccept);
    } else {
        Enable(Cn3D_bOpenAccept);
    }
    return;
}

static void Cn3D_NetOpenAcceptProc(ButtoN b)
{
    Char str[32];
    BiostrucPtr pbsBiostruc = NULL;
    PDNMS pdnmsModelstruc = NULL;
    PMSD pmsdThis = NULL;
    Int4 MdlNo, MdlLvl;

    WatchCursor();
    Cn3D_OpenStart();

    GetTitle(Cn3D_tOpen, str, sizeof(str));

    switch (GetValue(Cn3D_gMdlLvl)) {
    case 1:
        MdlLvl = ONECOORDATOM;
        break;
    case 2:
        MdlLvl = ONECOORDRES;
        break;
    case 3:
        MdlLvl = ALLMDL;
        break;
    case 4:
        MdlLvl = VECMODEL;
        break;
    case 5:
        MdlLvl = ALLSIMPLEMDL;
        break;
    case 6:
        MdlLvl = BSEVERYTHING;
        break;
    default:
        MdlLvl = ONECOORDATOM;  /* set from config-file default */
    }
    switch (GetValue(Cn3D_gMdlNo)) {
    case 1:
        MdlNo = 1;
        break;
    case 2:
        MdlNo = 2;
        break;
    case 3:
        MdlNo = 5;
        break;
    case 4:
        MdlNo = 10;
        break;
    case 5:
        MdlNo = 15;
        break;
    case 6:
        MdlNo = 20;
        break;
    case 7:
        MdlNo = MAX_MDLNO;
        break;
    case 8:
    default:
        MdlNo = 1;              /* set from config-file default */
    }
    if (GetValue(Cn3D_gMdlLvl) == 6)
        MdlNo = MAX_MDLNO;      /* get all */
    if (GetValue(Cn3D_gAccType) == 1) { /* PDB */
        pbsBiostruc = FetchBiostrucPDB(str, MdlLvl, MdlNo);
    } else {                    /* MMDB */

        pbsBiostruc =
            FetchBS(str, INP_GI, MdlLvl, MdlNo, GetMMDBAPIbExtent());
    }
    ArrowCursor();
    if (pbsBiostruc != NULL) {
        WatchCursor();
        pdnmsModelstruc = MakeAModelstruc(pbsBiostruc);
        ArrowCursor();
    }
    if (pdnmsModelstruc == NULL) {
        /* return a not found error here */
        Remove(Cn3D_wNetOpen);
        Cn3D_EnableFileOps();
        Cn3D_Open_InUse = FALSE;
        return;
    }
    if (GetValue(Cn3D_gMdlLvl) == 5) { /* turn off backbone model if "All" models present */
        pmsdThis = (PMSD) pdnmsModelstruc->data.ptrvalue;
        MMDB_OpenTraverse(pmsdThis);
    }

    Remove(Cn3D_wNetOpen);
    Cn3D_EnableFileOps();

    /*CALL TO initialize the view */

    Cn3D_Open_InUse = FALSE;

    Cn3D_OpenEnd();

    return;
}


static void Cn3D_NetOpenCancelProc(GraphiC g)
{
    Remove(Cn3D_wNetOpen);
    Cn3D_EnableFileOps();
    Cn3D_Open_InUse = FALSE;
    return;
}

static void Cn3D_NetOpenBiostruc(IteM i)
{
    GrouP g, hg;
    ButtoN b;

    if (Cn3D_Open_InUse)
        return;
    else
        Cn3D_Open_InUse = TRUE;

    Cn3D_StartNet();
    if (!Cn3D_ColorData.EntrezOn) return;

    Cn3D_wNetOpen =
        FixedWindow(-30, -20, -10, -10, " Internet retrieve from MMDB ",
                    NULL);
    hg = HiddenGroup(Cn3D_wNetOpen, 3, 0, NULL);
    SetGroupSpacing(hg, 30, 30);
    g = NormalGroup(hg, 1, 0, " Enter accession code:", systemFont, NULL);
    SetGroupMargins(g, 10, 15);
    Cn3D_tOpen = DialogText(g, "", 10, (TxtActnProc) Cn3D_OpenEnableProc);
    Cn3D_gAccType =
        NormalGroup(hg, 1, 2, " accession type", systemFont, NULL);
    SetGroupMargins(Cn3D_gAccType, 10, 10);
    RadioButton(Cn3D_gAccType, "PDB Code");
    RadioButton(Cn3D_gAccType, "MMDB ID");

    g = HiddenGroup(hg, 1, 2, NULL);
    SetGroupSpacing(g, 15, 15);
    Cn3D_bOpenAccept =
        DefaultButton(g, "OK", (BtnActnProc) Cn3D_NetOpenAcceptProc);
    b = PushButton(g, "Cancel", (BtnActnProc) Cn3D_NetOpenCancelProc);

    Cn3D_gMdlLvl =
        NormalGroup(Cn3D_wNetOpen, 2, 3, " model complexity", systemFont,
                    NULL);
    SetGroupMargins(Cn3D_gMdlLvl, 10, 10);
    SetGroupSpacing(Cn3D_gMdlLvl, 10, 5);
    RadioButton(Cn3D_gMdlLvl, "a) NCBI one XYZ per atom model");
    RadioButton(Cn3D_gMdlLvl, "c) NCBI backbone model");
    RadioButton(Cn3D_gMdlLvl, "b) original PDB models 1-n");
    RadioButton(Cn3D_gMdlLvl, "d) NCBI vector model");
    RadioButton(Cn3D_gMdlLvl, "Viewing Subset (a, c and d)");
    RadioButton(Cn3D_gMdlLvl, "Everything");

    Cn3D_gMdlNo =
        NormalGroup(Cn3D_wNetOpen, 7, 0, " n = ", systemFont, NULL);
    RadioButton(Cn3D_gMdlNo, "1");
    RadioButton(Cn3D_gMdlNo, "2");
    RadioButton(Cn3D_gMdlNo, "5");
    RadioButton(Cn3D_gMdlNo, "10");
    RadioButton(Cn3D_gMdlNo, "15");
    RadioButton(Cn3D_gMdlNo, "20");
    RadioButton(Cn3D_gMdlNo, "maximum");


    SetValue(Cn3D_gMdlNo, 7);
    SetValue(Cn3D_gAccType, 1);
    SetValue(Cn3D_gMdlLvl, 5);
    Disable(Cn3D_bOpenAccept);
    Cn3D_DisableFileOps();
    Select(Cn3D_tOpen);
    Show(Cn3D_wNetOpen);
    return;
}

/*********************************************/
/* below this are the file i/o open-er procs */
/*********************************************/

static void Cn3D_OpenBrowseProc(GraphiC g)
{
    Char path[PATH_MAX];

    path[0] = '\0';

    if (GetInputFileName(path, sizeof(path), NULL, NULL)) {
        SetTitle(Cn3D_tOpen, path);
        Cn3D_OpenEnableProc(NULL);
    }

    return;
}


static void Cn3D_OpenCancelProc(ButtoN b)
{
    Remove(Cn3D_wOpen);
    Cn3D_EnableFileOps();
    Cn3D_Open_InUse = FALSE;
    return;
}


static void Cn3D_OpenAcceptProc(ButtoN b)
{
    Char str[PATH_MAX];
    unsigned char szBegin[10];
    BiostrucPtr pbsBiostruc;
    PDNMS pdnmsModelstruc;
    PMSD pmsdThis = NULL;
    Int4 MdlNo = MAX_MDLNO;
    FILE *hFile;
    Char buf[50];



    WatchCursor();
    GetTitle(Cn3D_tOpen, str, sizeof(str));
    hFile = FileOpen(str, "rb");
    FileGets((CharPtr) szBegin, 2, hFile);
    if (hFile == NULL)
        szBegin[0] = (Char) 0;
    FileClose(hFile);
    /* to make cn3d to take strucseq for which szBegin[0] is 78 */
    /* for single biostruc szBegin[0] is 48 */
    if (szBegin[0] > 70) { /* mime */
        if (!OpenMimeFileWithDeletion(str, FALSE)) {
            Remove(Cn3D_wOpen);
            Cn3D_EnableFileOps();
            Cn3D_Open_InUse = FALSE;
            return;
        }
    } else {                    /* not mime */
        
        Cn3D_OpenStart();
        hFile = FileOpen(str, "r");
        
        FileRead(buf, 1, StrLen(PRINT_FORM_BIOSTRUC), hFile);
        FileClose(hFile);
        if (StrNCmp(buf, PRINT_FORM_BIOSTRUC, StrLen(PRINT_FORM_BIOSTRUC)) ==
            0) { /* ASCII */
            /* these get everything in the file , ignoring modellevel */
            pbsBiostruc =
                FetchBS(str, INP_ASCII_FILE, BSEVERYTHING, MdlNo,
                CONVERT_ALL);
        } else {                /* Binary */
            pbsBiostruc =
                FetchBS(str, INP_BINARY_FILE, BSEVERYTHING, MdlNo,
                CONVERT_ALL);
        }
        ArrowCursor();
        if (pbsBiostruc != NULL) {
            WatchCursor();
            pdnmsModelstruc = MakeAModelstruc(pbsBiostruc);
            pmsdThis = (PMSD) pdnmsModelstruc->data.ptrvalue;
            MMDB_OpenTraverse(pmsdThis);

            ArrowCursor();
        } else {
            /* return a not found error here */

            Remove(Cn3D_wOpen);
            Cn3D_EnableFileOps();
            Cn3D_Open_InUse = FALSE;
            return;
        }

        Cn3D_OpenEnd();
    }                           /* switch between mime and non-mime */
    Remove(Cn3D_wOpen);
    Cn3D_EnableFileOps();
    Cn3D_Open_InUse = FALSE;

    return;
}

static void Cn3D_OpenBiostruc(IteM i)
{
    GrouP g;
    ButtoN b;

    if (Cn3D_Open_InUse)
        return;
    else
        Cn3D_Open_InUse = TRUE;

    Cn3D_wOpen =
        FixedWindow(-30, -20, -10, -10, " Open a local Biostruc ", NULL);
    g =
        NormalGroup(Cn3D_wOpen, 2, 1, " Enter Biostruc file name:",
                    systemFont, NULL);
    SetGroupMargins(g, 10, 10);
    SetGroupSpacing(g, 10, 20);
    Cn3D_tOpen = DialogText(g, "", 25, (TxtActnProc) Cn3D_OpenEnableProc);
    Cn3D_bOpenBrowse =
        PushButton(g, " browse...", (BtnActnProc) Cn3D_OpenBrowseProc);

    g = HiddenGroup(Cn3D_wOpen, 3, 1, NULL);
    SetGroupMargins(g, 10, 10);
    SetGroupSpacing(g, 30, 30);

/*    Cn3D_gBinAscii = NormalGroup(g, 2, 1, "file mode", systemFont, NULL);
    SetGroupMargins(Cn3D_gBinAscii, 10, 10);
    RadioButton(Cn3D_gBinAscii, "Ascii");
    RadioButton(Cn3D_gBinAscii, "Binary");
    SetValue(Cn3D_gBinAscii, 2);*/

    b = PushButton(g, "Cancel", (BtnActnProc) Cn3D_OpenCancelProc);
    Cn3D_bOpenAccept =
        DefaultButton(g, "OK", (BtnActnProc) Cn3D_OpenAcceptProc);

    Disable(Cn3D_bOpenAccept);
    Cn3D_DisableFileOps();
    Select(Cn3D_tOpen);
    Show(Cn3D_wOpen);
    return;
}


MenU LIBCALL Cn3D_OpenSub(MenU m)
{
    IteM i;
    MenU s;

    s = SubMenu(m, "Open");
    i = CommandItem(s, "Network/N", Cn3D_NetOpenBiostruc);
    if (!Cn3D_ColorData.UseEntrez) {
        Disable(i);
    }
    i = CommandItem(s, "Local File/B", Cn3D_OpenBiostruc);

    return s;
}
