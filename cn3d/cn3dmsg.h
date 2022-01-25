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
* Revision 6.41  1999/04/06 20:12:57  lewisg
* more opengl
*
* Revision 6.40  1999/03/30 22:36:19  ywang
* add functions to color salsa for NcbiMimeAsn1_strucseqs & code reorganization
*
* Revision 6.39  1999/03/22 22:41:14  ywang
* remove argument in MediaObjSelect
*
* Revision 6.38  1999/02/12 15:34:01  ywang
* include Cn3DSendColorMsg
*
* Revision 6.37  1999/02/11 22:40:15  ywang
* rename functions
*
* Revision 6.36  1999/02/11 18:48:15  lewisg
* delete color index functions
*
* Revision 6.35  1999/02/10 23:49:43  lewisg
* use RGB values instead of indexed palette
*
* Revision 6.34  1999/01/20 22:57:25  ywang
* customize color for secondary structure & rearrange Option menu
*
* Revision 6.33  1999/01/20 16:06:48  durand
* move mediainfo to salmedia.h
*
* Revision 6.32  1998/12/16 22:49:39  ywang
* fix compiling warnings on Win32
*
 * Revision 6.31  1998/12/16  19:32:56  ywang
 * improve highlight residues function when rerendering
 *
 * Revision 6.30  1998/10/27  15:55:52  ywang
 * add functions for testing color by sequence conservation
 *
 * Revision 6.29  1998/10/19  20:16:06  ywang
 * add function FillSeqinfoForSeqEditViewProcs so that salsa can get color array
 *
 * Revision 6.28  1998/10/19  17:43:03  kans
 * prototype needed for RealColorSalsa
 *
* Revision 6.27  1998/10/16 22:06:08  ywang
* make global color array for sequence display
*
 * Revision 6.26  1998/09/23  22:09:42  ywang
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
#include <salmedia.h>

#ifdef _OPENGL
#include <shim3d.h>
#else
#include <viewer3d.h>
#endif

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


/* extern MediaInfo mediadata[MaxObj];  */
extern MediaInfo **mediadata;  /* put no limitation on mediadata size in principle */

typedef struct mediaview{
FORM_MESSAGE_BLOCK
} MediaView, PNTR MediaViewPtr;

#define REGISTER_BIOSEQ_BIOSTRUC_MEDIA ObjMgrProcLoad(OMPROC_VIEW, "Seq-Struc Communication", "Media", OBJ_BIOSEQ, 0, OBJ_BIOSEQ, 0, NULL, SeqStrucMediaFunc, 0)

extern Int2 LIBCALLBACK SeqStrucMediaFunc PROTO((Pointer data));
extern void MediaObjSelect(PDNMG pdnmgThis, Boolean highlight);
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
extern void fnCHLresidue(PDNMG pdnmgThis,
#ifndef _OPENGL
                         Viewer3D  vvv,
#endif
                         Boolean highlight);
extern void LaunchMediaViewer(BioseqPtr bsp);
extern void SalsaRegister(void);
extern void LaunchSalsa(SeqAlignPtr salp);
extern void ColorSalsa_BYMG(PMGD pmgdThis, Uint1Ptr rgb);
extern void ResetSalsaColor(void);
extern void Cn3DSendColorMsg(void);
extern void ColorSalsa(void);
extern void Cn3dObjRegister(void);
extern void Cn3dObjMgrGetSelected(void);
extern void DoMediaHL(PMMD  pmmdThis, Int4 from, Int4 to, Boolean highlight);
extern void LaunchSequenceWindow(void);
extern void Cn3DLaunchAnnotAlignEditor (SeqAnnotPtr sap);
extern ResidueColorCell * GetColorIndexForMG(PMGD pmgdThis);
extern void LIBCALLBACK Cn3DCheckAndDoHighlight PROTO((PFB pfbThis,Int4 iModel, Int4 iIndex, Pointer ptr));
extern void Cn3DSendColorMsgForBioseq(Int4 iCount);
extern void Cn3DColorSalsaForStrucSeqs(void);
extern void Cn3DCheckAlignmentStatusForStrucSeqs(void);
extern void Cn3DCheckAlignmentStatusForStrucSeqsForMasterSeq(void);
#ifdef __cplusplus
}
#endif


#endif /* _CN3DMSG_ */

