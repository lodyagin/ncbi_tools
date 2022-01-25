/*   cn3dsave.c
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
* File Name:  cn3dsave.c
*
* Author:  Christopher Hogue, Yanli Wang, Lewis Geer
*
* First Version Creation Date:   1/31/96
*
* $Revision: 6.17 $
*
* File Description: Cn3d file saving routines 
*                   
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
* $Log: cn3dsave.c,v $
* Revision 6.17  1999/11/22 14:46:44  thiessen
* moved _OPENGL code blocks to only vibrant and ncbicn3d libraries
*
* Revision 6.16  1999/11/11 01:18:40  kans
* cannot valnodeaddpointer to a seqannot - reverting to old instructions for now
*
* Revision 6.15  1999/11/10 23:19:42  lewisg
* rewrite of selection code for ddv
*
* Revision 6.14  1999/10/29 14:15:30  thiessen
* ran all Cn3D source through GNU Indent to prettify
*
* Revision 6.13  1999/10/18 15:32:50  lewisg
* move ClearSequences() to cn3dshim.c
*
* Revision 6.12  1999/10/05 23:18:24  lewisg
* add ddv and udv to cn3d with memory management
*
* Revision 6.11  1999/08/04 21:18:01  lewisg
* modularized open operations to allow sequin to launch cn3d
*
* Revision 6.10  1999/07/01 22:14:47  ywang
* back to previous version to use EntrezSeqEntryGet to get sequences for additional loading, before this free existant sequences
*
* Revision 6.9  1999/07/01 21:44:31  ywang
* work around EntrezSeqEntryGet core dump on redunant sequence loading and free existant sequences before additional loading
*
* Revision 6.8  1999/07/01 14:09:05  ywang
* *** empty log message ***
*
* Revision 6.7  1999/07/01 14:07:25  ywang
* cn3dwin.c
*
* Revision 6.6  1999/06/03 21:31:48  ywang
* fix bug for saving structures--get MS coordinates back for the current run after creating ASN.1 for saving
*
* Revision 6.5  1999/03/19 19:20:11  kans
* prototype needed for Cn3DAddUserDefinedFeature
*
* Revision 6.4  1999/03/18 22:28:58  ywang
* add functions for saveout+readin+index user defined features
*
* Revision 6.3  1999/03/18 21:04:28  lewisg
* reverse transform for saving
*
* Revision 6.2  1999/02/26 15:44:07  kans
* boolean initialized to false, not null
*
* Revision 6.1  1999/02/25 23:13:41  ywang
* make cn3d export the entire input mime data
*
 * Revision 6.0  1997/08/25  18:13:42  madden
 * Revision changed to 6.0
 *
* Revision 5.1  1996/07/22 00:28:11  hogue
* Fixed binary ASN.1 saving.
*
 * Revision 5.0  1996/05/28  14:05:44  ostell
 * Set to revision 5.0
 *
 * Revision 1.3  1996/04/18  17:02:20  hogue
 * removed temporary reference to mmdbapi4.h
 *
 * Revision 1.2  1996/02/02  19:40:41  hogue
 * Initial Revision
 *
*
* ==========================================================================
*/

#include <vibrant.h>
#include <mmdbapi.h>
#include <mmdbapi1.h>
#include <cn3dmain.h>
#include <algorend.h>
#include <cn3dmodl.h>
#include <cn3dsave.h>
#include <cn3dmsel.h>
#include <cn3dmsg.h>
#include <asnmime.h>
#include <objmime.h>
#include <cn3dopen.h>
#include <cn3dshim.h>

static Boolean Cn3D_Save_InUse = FALSE;

static WindoW Cn3D_wAsnSave;
static TexT Cn3D_tAsnSave;
static ButtoN Cn3D_bAsnBrowse;
static ButtoN Cn3D_bAsnOk;
static GrouP Cn3D_gBinAscii;
static ButtoN Cn3D_bFeatOn;


static void Cn3D_AsnEnableProc(TexT t)
{
    Char str[32];
    GetTitle(Cn3D_tAsnSave, str, sizeof(str));
    if (StringLen(str) == 0) {
        Disable(Cn3D_bAsnOk);
    } else {
        Enable(Cn3D_bAsnOk);
    }
    return;
}

static void Cn3D_AsnBrowseProc(ButtoN b)
{
    Char dfault[32];
    Char path[256];

    path[0] = '\0';
    dfault[0] = '\0';
    if (GetOutputFileName(path, sizeof(path), dfault)) {
        SetTitle(Cn3D_tAsnSave, path);
        Cn3D_AsnEnableProc(NULL);
    }
    return;
}

static void Cn3D_ExportAsnNow(ButtoN b)
{

    Char path[256];
    Int2 iTest;
    Int4 iCount = 0;
    CharPtr pcSave = NULL;
    Byte bSave = 0;

    Int2Ptr i2Vec = NULL;

    Boolean iCn3d = TRUE, FirstSlave = FALSE;
    CharPtr pcSaveMode = NULL;

    PDNMS pdnmsMaster = NULL, pdnmsSlave = NULL;
    PMSD pmsdSlave = NULL, pmsdMaster = NULL;
    NcbiMimeAsn1Ptr mime = NULL;
    BiostrucAlignPtr bsap = NULL;
    BiostrucSeqPtr bssp = NULL;
    BiostrucSeqsPtr bsssp = NULL;
    BiostrucPtr bsSlaveHead = NULL, bsSlave = NULL, bsSlaveTail = NULL;
    BiostrucFeaturePtr pbsfThis;
    PDNTRN pdnTransform = NULL;
    AsnIoPtr aip = NULL;
    ValNodePtr pvnAlignment;
    PDNML pdnmlThis = NULL;

/*   if (GetStatus(Cn3D_bFeatOn) == FALSE)
       bSave = (Byte) (bSave | (Byte) NOT_FEATURES);  */
    /* save feature always for now -yanli */

    if (GetValue(Cn3D_gBinAscii) == 2)
        bSave = (Byte) (bSave | (Byte) SAVE_BINARY);

    if (bSave & (Byte) SAVE_BINARY)
        pcSaveMode = SAVE_BINARY_STRING;
    else
        pcSaveMode = SAVE_ASCII_STRING;

    i2Vec = PickedModels(&iCount);
    GetTitle(Cn3D_tAsnSave, path, sizeof(path));
    pdnmsMaster = GetSelectedModelstruc();
    iTest =
        WriteAsnModelList(pdnmsMaster, iCount, i2Vec, path, bSave, iCn3d);
    pdnmsMaster = (PDNMS) Cn3DAddUserDefinedFeature(pdnmsMaster);
    if (!iTest) {
        ErrClear();
        ErrPostEx(SEV_FATAL, 0, 0,
                  "Unable to Export\nPossibly Corrupt Data in Memory!\n");
        ErrShow();

    }

    pmsdMaster = pdnmsMaster->data.ptrvalue;

    aip = AsnIoOpen(path, pcSaveMode);

    mime = (NcbiMimeAsn1Ptr) ValNodeNew(NULL);
    switch (pmsdMaster->iMimeType) {
    case NcbiMimeAsn1_entrez:
        BiostrucAsnWrite(pmsdMaster->pbsBS, aip, NULL);
        break;
    case NcbiMimeAsn1_strucseq:
        mime->choice = NcbiMimeAsn1_strucseq;
        bssp = BiostrucSeqNew();
        bssp->structure = pmsdMaster->pbsBS;
        bssp->sequences = Cn3D_ColorData.pvnsep;
        mime->data.ptrvalue = bssp;
        NcbiMimeAsn1AsnWrite(mime, aip, NULL);
        break;
    case NcbiMimeAsn1_strucseqs:
        mime->choice = NcbiMimeAsn1_strucseqs;
        bsssp = BiostrucSeqsNew();
        bsssp->structure = pmsdMaster->pbsBS;
        bsssp->sequences = Cn3D_ColorData.pvnsep;
        /* ValNodeAddPointer(bsssp->seqalign, 0, Cn3D_ColorData.sap); */
        bsssp->seqalign = Cn3D_ColorData.sap;
        mime->data.ptrvalue = bsssp;
        NcbiMimeAsn1AsnWrite(mime, aip, NULL);
        break;
    case NcbiMimeAsn1_alignstruc:
        mime->choice = NcbiMimeAsn1_alignstruc;
        bsap = BiostrucAlignNew();
        bsap->master = pmsdMaster->pbsBS;
        pdnmsSlave = pmsdMaster->pdnmsSlaves;
        FirstSlave = TRUE;

        pbsfThis =
            (BiostrucFeaturePtr) pmsdMaster->psaStrucAlignment->features->
            features;


        while (pdnmsSlave) {

            /* begin reverse transform */

            pmsdSlave = pdnmsSlave->data.ptrvalue;
            pdnTransform = NULL;
            pvnAlignment =
                ValNodeFindNext(pbsfThis->Location_location, NULL,
                                Location_location_alignment);
            if (pvnAlignment == NULL)
                break;
            /* create the spatial transformation */
            TransformToDNTRN(&pdnTransform,
                             ((ChemGraphAlignmentPtr) pvnAlignment->data.
                              ptrvalue)->transform);
            if (pdnTransform == NULL)
                break;
            /* loop over the slave's models with the transformation */
            pdnmlThis = pmsdSlave->pdnmlModels;

            while (pdnmlThis) {
                TraverseAtoms(pdnmsSlave, pdnmlThis->choice, 0,
                              pdnTransform, DoReverseTransform);
                TraverseSolids(pdnmsSlave, pdnmlThis->choice, 0,
                               pdnTransform, DoReverseTransform);
                pdnmlThis = pdnmlThis->next;
            }
            pbsfThis = pbsfThis->next;

            /* end reverse transform */

            iTest =
                WriteAsnModelList(pdnmsSlave, iCount, i2Vec, path, bSave,
                                  iCn3d);
            pdnmlThis = pmsdSlave->pdnmlModels;
            while (pdnmlThis) {
                TraverseAtoms(pdnmsSlave, pdnmlThis->choice, 0,
                              pdnTransform, DoApplyTransform);
                TraverseSolids(pdnmsSlave, pdnmlThis->choice, 0,
                               pdnTransform, DoApplyTransform);
                pdnmlThis = pdnmlThis->next;
            }
            /* after creating ASN.1, get the coordinates back for this run, yanli */
            FreeDNTRN(pdnTransform);

            if (!iTest) {
                ErrClear();
                ErrPostEx(SEV_FATAL, 0, 0,
                          "Unable to Export\nPossibly Corrupt Data in Memory!\n");
                ErrShow();

            }

            if (!bsSlaveHead) {
                bsSlaveHead =
                    (BiostrucPtr) ((PMSD) pdnmsSlave->data.ptrvalue)->
                    pbsBS;
                bsSlaveTail = bsSlaveHead;
            } else {
                bsSlave =
                    (BiostrucPtr) ((PMSD) pdnmsSlave->data.ptrvalue)->
                    pbsBS;
                bsSlaveTail->next = bsSlave;
                bsSlaveTail = bsSlaveTail->next;
                bsSlaveTail->next = NULL;
            }
            pdnmsSlave = (PDNMS) Cn3DAddUserDefinedFeature(pdnmsSlave);
            pdnmsSlave = pdnmsSlave->next;
        }
        bsap->slaves = bsSlaveHead;
        bsap->sequences = Cn3D_ColorData.pvnsep;
        bsap->alignments = pmsdMaster->psaStrucAlignment;
        /* ValNodeAddPointer(bsap->seqalign, 0, Cn3D_ColorData.sap); */
        bsap->seqalign = Cn3D_ColorData.sap;
        mime->data.ptrvalue = bsap;
        NcbiMimeAsn1AsnWrite(mime, aip, NULL);
        break;
    }

    aip = AsnIoClose(aip);

    FreeRedundantAsn(pdnmsMaster);
    pdnmsSlave = pmsdMaster->pdnmsSlaves;
    while (pdnmsSlave) {
        FreeRedundantAsn(pdnmsSlave);
        pdnmsSlave = pdnmsSlave->next;
    }

    if (i2Vec)
        I2VectorFree(i2Vec, 0);
    Remove(Cn3D_wAsnSave);
    Cn3D_EnableFileOps();
    Cn3D_Save_InUse = FALSE;
    ArrowCursor();
    return;
}

static void Cn3D_CancelAsn(ButtoN b)
{
    Remove(Cn3D_wAsnSave);
    Cn3D_EnableFileOps();
    Cn3D_Save_InUse = FALSE;
    return;
}



static void Cn3D_SaveBiostruc(IteM i)
{
    PDNMS pdnmsThis = NULL;
    PMSD pmsdThis = NULL;
    Char pcSavestr[60];
    Char pcSavename[32];
    CharPtr Cn3D_pcAsnName;
    GrouP g, g2, g3, g4, g5, g6;
    GrouP gMS;
    ButtoN b;

    if (Cn3D_Save_InUse)
        return;
    else
        Cn3D_Save_InUse = TRUE;

    pdnmsThis = GetSelectedModelstruc();
    if (!pdnmsThis) {
        Cn3D_Save_InUse = FALSE;
        return;
    }
    pmsdThis = (PMSD) pdnmsThis->data.ptrvalue;

    Cn3D_wAsnSave = ModalWindow(-20, -13, -10, -10, NULL);

    /* set up a group encolosing structures - models selection lists and - "info strings" */
    Cn3D_pcAsnName = StringSave(GetStrucStrings(pdnmsThis, PDB_ACC));
    sprintf(pcSavestr, "Save Biostruc %s as an Asn.1 File...",
            Cn3D_pcAsnName);
    g = HiddenGroup(Cn3D_wAsnSave, 0, 4, NULL);
    SetGroupMargins(g, 10, 10);
    SetGroupSpacing(g, 10, 5);
    StaticPrompt(g, pcSavestr, 0, 0, systemFont, 'l');
    g2 = HiddenGroup(g, 2, 0, NULL);
    SetGroupMargins(g2, 10, 10);
    SetGroupSpacing(g2, 10, 5);
    StringNCpy(pcSavename, Cn3D_pcAsnName, 8);
    StringCat(pcSavename, ".prt");
    Cn3D_tAsnSave =
        DialogText(g2, pcSavename, 18, (TxtActnProc) Cn3D_AsnEnableProc);
    MemFree(Cn3D_pcAsnName);

    Cn3D_bAsnBrowse =
        PushButton(g2, " browse...", (BtnActnProc) Cn3D_AsnBrowseProc);
    g3 = HiddenGroup(g, 2, 0, NULL);
    gMS = Cn3D_ModelSelect(g3, TRUE); /*  vector models OK for Asn file saves */
    g4 = HiddenGroup(g3, 0, 2, NULL);
    SetGroupMargins(g4, 10, 10);
    SetGroupSpacing(g4, 10, 5);
    Cn3D_bAsnOk = PushButton(g4, "OK", Cn3D_ExportAsnNow);
    b = PushButton(g4, "Cancel", Cn3D_CancelAsn);
    g5 = HiddenGroup(g, 2, 0, NULL); /* for bin/ascii and features on/off */
    Cn3D_gBinAscii = NormalGroup(g5, 2, 0, "file mode", systemFont, NULL);
    SetGroupMargins(Cn3D_gBinAscii, 10, 0);
    SetGroupSpacing(Cn3D_gBinAscii, 10, 5);
    RadioButton(Cn3D_gBinAscii, "Ascii");
    RadioButton(Cn3D_gBinAscii, "Binary");
    SetValue(Cn3D_gBinAscii, 1);
    g6 = HiddenGroup(g5, 0, 2, NULL);
    SetGroupMargins(g6, 10, 10);
    StaticPrompt(g6, " ", 0, 0, systemFont, 'l');
/*  Cn3D_bFeatOn =  CheckBox(g6, "Include Features", NULL);
    SetStatus(Cn3D_bFeatOn,TRUE);  */
    /* save features always for now - yanli */
    Cn3D_AsnEnableProc(NULL);
    Select(Cn3D_bAsnOk);
    /* disable appropriate stuff here */
    Cn3D_DisableFileOps();
    Show(Cn3D_wAsnSave);

    return;
}




static void Cn3D_ExportDictNow(ButtoN b)
{

    Char path[256];
    Int2 iTest;
    Int4 iCount = 0;
    PDNMS pdnmsMain = NULL;
    CharPtr pcSave = NULL;
    Byte bSave = 0;

    GetTitle(Cn3D_tAsnSave, path, sizeof(path));
    pdnmsMain = GetSelectedModelstruc();
    if (GetValue(Cn3D_gBinAscii) == 2)
        bSave = (Byte) SAVE_BINARY;
    iTest = WriteAsnLocalDict(pdnmsMain, path, bSave, TRUE);
    if (!iTest) {
        ErrClear();
        ErrPostEx(SEV_FATAL, 0, 0,
                  "Unable to Export\nPossibly Corrupt Data in Memory!\n");
        ErrShow();
    }
    Remove(Cn3D_wAsnSave);
    Cn3D_EnableFileOps();
    Cn3D_Save_InUse = FALSE;
    ArrowCursor();
    return;
}



static void Cn3D_SaveDictionary(IteM i)
{
    PDNMS pdnmsThis = NULL;
    PMSD pmsdThis = NULL;
    Char pcSavestr[60];
    Char pcSavename[32];
    CharPtr Cn3D_pcAsnName;
    GrouP g, g2, g3, g4, g5;
    ButtoN b;

    if (Cn3D_Save_InUse)
        return;
    else
        Cn3D_Save_InUse = TRUE;

    pdnmsThis = GetSelectedModelstruc();
    if (!pdnmsThis) {
        Cn3D_Save_InUse = FALSE;
        return;
    }
    pmsdThis = (PMSD) pdnmsThis->data.ptrvalue;

    Cn3D_wAsnSave = ModalWindow(-20, -13, -10, -10, NULL);

    Cn3D_pcAsnName = StringSave(GetStrucStrings(pdnmsThis, PDB_ACC));
    sprintf(pcSavestr, "Save Local Dictionary from %s ...",
            Cn3D_pcAsnName);
    g = HiddenGroup(Cn3D_wAsnSave, 0, 3, NULL);
    SetGroupMargins(g, 10, 10);
    SetGroupSpacing(g, 10, 5);
    StaticPrompt(g, pcSavestr, 0, 0, systemFont, 'l');
    g2 = HiddenGroup(g, 2, 0, NULL);
    SetGroupMargins(g2, 10, 10);
    SetGroupSpacing(g2, 10, 5);
    StringNCpy(pcSavename, Cn3D_pcAsnName, 8);
    StringCat(pcSavename, "dict.prt");
    Cn3D_tAsnSave =
        DialogText(g2, pcSavename, 18, (TxtActnProc) Cn3D_AsnEnableProc);
    MemFree(Cn3D_pcAsnName);
    Cn3D_bAsnBrowse =
        PushButton(g2, " browse...", (BtnActnProc) Cn3D_AsnBrowseProc);

    g3 = HiddenGroup(g, 2, 0, NULL);
    g4 = HiddenGroup(g3, 2, 0, NULL);
    SetGroupMargins(g4, 10, 0);
    Cn3D_gBinAscii = NormalGroup(g4, 2, 1, "file mode", systemFont, NULL);
    SetGroupMargins(Cn3D_gBinAscii, 10, 10);
    RadioButton(Cn3D_gBinAscii, "Ascii");
    RadioButton(Cn3D_gBinAscii, "Binary");
    SetValue(Cn3D_gBinAscii, 1);

    g5 = HiddenGroup(g3, 2, 0, NULL);
    SetGroupMargins(g5, 10, 10);
    SetGroupSpacing(g5, 10, 5);
    Cn3D_bAsnOk = PushButton(g5, "OK", Cn3D_ExportDictNow);
    b = PushButton(g5, "Cancel", Cn3D_CancelAsn);

    Cn3D_AsnEnableProc(NULL);
    Select(Cn3D_bAsnOk);
    Cn3D_DisableFileOps();
    Show(Cn3D_wAsnSave);

    return;
}




MenU LIBCALL Cn3D_SaveSub(MenU m)
{
    IteM i;
    MenU s;

    s = SubMenu(m, "Save");
    i = CommandItem(s, "All/B", Cn3D_SaveBiostruc);
    /* i = CommandItem (s, "Feature-set/F", NULL); */
    i = CommandItem(s, "Dictionary/D", Cn3D_SaveDictionary);
    return s;
}
