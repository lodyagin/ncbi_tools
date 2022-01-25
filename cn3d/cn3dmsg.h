/*   cn3dmsg.h
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
* File Name:  cn3dmsg.h
*
* Author: Yanli Wang
*
* Version Creation Date:   3/26/98
*
* File Description: Main functions for building up cn3d/salsa communication
*
* Modifications:
* $Log: cn3dmsg.h,v $
* Revision 6.26  1998/09/23 22:09:42  ywang
* to record checkin log
*
* ===========================================================================  */

#ifndef _CN3DMSG_
#define _CN3DMSG_ 1

#include <ncbi.h>
#include <objseq.h>
#include <objmgr.h>
#include <objfdef.h>
#include <gather.h>
#include <vibrant.h>
#include <salsa.h>
#include <viewer3d.h>
#include <mmdbapi1.h>

#ifdef __cplusplus
extern "C" {
#endif

#define MaxObj 16

extern Boolean Cn3D_ObjMgrOpen;
extern Boolean Cn3D_SalsaOpen;
extern Boolean Salsa_BioseqUpdate;
extern Boolean Cn3D_ReColor;
extern Int4 Num_Bioseq, Num_Biostruc;
extern Int4 Num_ActiveSlave;
extern Boolean Mime_ReadIn;

extern Uint2 sap_entityID, sap_itemID;

typedef struct mediainfo{
SeqIdPtr sip;
Int4 Gi, length;
Uint2 entityID, itemID;
Byte bVisible;
} MediaInfo, PNTR MediaInfoPtr;

/* extern MediaInfo mediadata[MaxObj];  */
extern MediaInfo **mediadata;  /* put no limitation on mediadata size in principle */

typedef struct mediaview{
FORM_MESSAGE_BLOCK
} MediaView, PNTR MediaViewPtr;

#define REGISTER_BIOSEQ_BIOSTRUC_MEDIA ObjMgrProcLoad(OMPROC_VIEW, "Seq-Struc Communication", "Media", OBJ_BIOSEQ, 0, OBJ_BIOSEQ, 0, NULL, SeqStrucMediaFunc, 0)

extern Int2 LIBCALLBACK SeqStrucMediaFunc PROTO((Pointer data));
extern void MediaObjSelect(PDNMG pdnmgThis, PALD paldThis, Boolean highlight);
extern void MediaRegister(void);
extern void MediaLaunch(void);
extern void MediaDataLoad(SeqAlignPtr salp);
extern void MediaDataLoad2(PDNMM pdnmmHead);
extern void MediaHL(SelStructPtr sel, Boolean highlight);
extern void fnCHLresidueRedraw(void);
extern void fnPreCHLresidue(PDNMG pdnmgThis, Boolean highlight);
/* extern PMMD FindMM(Int4 Gi, Int4 iCount); */
extern PMMD FindMM(SeqIdPtr sip, Int4 iCount);   
extern void DoCHighlightSeg(PFB pfbThis, Int4 iModel, Int4 iIndex, Pointer ptr, Boolean highlight);
extern void fnCHLresidue(PDNMG pdnmgThis, Viewer3D  vvv, Boolean highlight);
extern void LaunchMediaViewer(BioseqPtr bsp);
extern void SalsaRegister(void);
extern void LaunchSalsa(SeqAlignPtr salp);
extern void ColorSalsa(Uint2 entityID, Uint2 itemID, SeqIdPtr sip, Int4 from, Int4 to, Uint1Ptr rgb);
extern void ColorSalsa_BYMG(PMGD pmgdThis, Uint1Ptr rgb);
extern void ResetSalsaColor(void);
extern Uint1Ptr GetRGB(Int2 iColor);
extern void Cn3dObjRegister(void);
extern void PrepareColorMsg(Int4 iCount, Int4 from, Int4 to, Uint1Ptr rgb);
extern void Cn3dObjMgrGetSelected(void);
extern void DoMediaHL(PMMD  pmmdThis, Int4 from, Int4 to, Boolean highlight);
extern void LaunchSequenceWindow(void);
extern void Cn3DLaunchAnnotAlignEditor (SeqAnnotPtr sap);
#ifdef __cplusplus
}
#endif


#endif /* _CN3DMSG_ */

