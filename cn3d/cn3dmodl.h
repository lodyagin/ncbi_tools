/*   cn3dmodl.h
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
* $Log: cn3dmodl.h,v $
* Revision 6.19  1999/07/07 20:45:36  ywang
* clear domaindata, mediadata, special feature before reading in new data in cn3d
*
* Revision 6.18  1999/06/07 19:39:51  ywang
* remove obsolete user defined features, an obsolete feature means a feature whose region has been completely overwirtten by a later defined feature
*
* Revision 6.17  1999/03/22 23:44:22  kans
* added Cn3DAddUserDefinedFeature to header
*
* Revision 6.16  1999/03/22 22:41:52  ywang
* redesign feature page, fix bugs
*
* Revision 6.14  1999/03/18 22:28:57  ywang
* add functions for saveout+readin+index user defined features
*
* Revision 6.13  1999/02/11 22:38:51  ywang
* fix bug on display highlight residues only--if no res are highlighted, cn3d sets that button status as FALSE and draw whole structurescn3dwin.c
*
* Revision 6.12  1999/02/04 16:16:50  ywang
* support delete added features
*
* Revision 6.10  1999/02/01 20:43:26  ywang
* improve 'Model' menu
*
* Revision 6.9  1999/01/26 17:14:35  ywang
* redesign Display menu and add 'display highlight residues only' function
*
 * Revision 6.8  1998/12/16  22:49:38  ywang
 * fix compiling warnings on Win32
 *
 * Revision 6.7  1998/11/04  00:06:20  ywang
 * add function for modeling: change render/color for special residue(s)
 *
 * Revision 6.6  1998/10/21  15:48:24  ywang
 * rearrange View Control menu
 *
 * Revision 6.5  1998/10/20  17:50:38  ywang
 * fix MMDB/VAST domain number inconsistency problem
 *
 * Revision 6.4  1998/10/07  23:10:45  ywang
 * merge align control with general display control
 *
 * Revision 6.3  1998/09/30  22:10:46  ywang
 * control display on three levels: structure, chain, domain
 *
 * Revision 6.2  1998/09/23  18:38:49  ywang
 * add functions to control display on domain level
 *
 * Revision 6.1  1998/09/22  18:02:54  ywang
 * panels and functions for display control
 *
*/


#ifndef _CN3DMODEL_
#define _CN3DMODEL_ 1
  
#ifdef __cplusplus
extern "C" {
#endif

typedef struct domain_info{
Char pcPDBName[20];      /* PDB code */
Char pcMolName[2];     /* Chain Id   */
Int2 iDomain;       /* domain number */
Int4 iStrucIndex, iChainIndex;       /* domain number */
Boolean bVisible, bVisibleParent, bAligned;
} DomainInfo, PNTR DomainInfoPtr;

typedef struct special_feature_info {
PARS parsSpecial;
/* PRK  prkSpecial;      */
Boolean On;
CharPtr title, description;
Int4 iRes;
} SpecialFeatureInfo, PNTR SpecialFeatureInfoPtr;

typedef ValNodePtr SpecialFeaturePtr;

extern DomainInfo **domaindata; 
extern Boolean Cn3D_DisplayHighlight;
extern Boolean Cn3D_NoSingleHL;

extern GrouP LIBCALL DisplayControls PROTO((Nlm_GrouP prnt));  
extern GrouP LIBCALL ModelControls PROTO((Nlm_GrouP prnt));  
extern void LIBCALL ResetDisplayCtrls(void);
extern void LIBCALL ResetModelCtrls(void);
extern void LIBCALL Cn3D_CountDomainProc(void);
extern void LIBCALLBACK DoLinkSpecialFeatureWithMGD  PROTO((PFB pfbThis,Int4 iModel, Int4 iIndex, Pointer ptr));
extern void  LIBCALLBACK DoCleanSpecialFeatureWithMGD(PFB pfbThis,Int4 iModel, Int4 iIndex, Pointer ptr);
extern void LIBCALLBACK DoCleanJustHLStatusWithMGD(PFB pfbThis,Int4 iModel, Int4 iIndex, Pointer ptr);
extern void LIBCALLBACK DoTurnOnSpecialFeatureWithMGD(PFB pfbThis,Int4 iModel, Int4 iIndex, Pointer ptr);
extern void LIBCALLBACK DoDeHighlightWithMGD(PFB pfbThis,Int4 iModel, Int4 iIndex, Pointer ptr);
extern SpecialFeaturePtr LIBCALL SpecialFeatureFree(SpecialFeaturePtr sfpThis);
extern void LIBCALLBACK DoUnLinkSpecialFeatureWithMGD(PFB pfbThis,Int4 iModel, Int4 iIndex, Pointer ptr);
extern void Cn3D_bDisplayHighlightStatusSet(Boolean Yes);
extern void Cn3DCheckHighlighted(void);
extern void LIBCALLBACK Cn3DCheckNoSingleHighlight(PFB pfbThis,Int4 iModel, Int4 iIndex, Pointer ptr);
extern PDNMS Cn3DAddUserDefinedFeature(PDNMS pdnmsThis);
extern void Cn3DIndexUserDefinedFeature(void);
extern void ClearRest(void);
#ifdef __cplusplus
}
#endif

#endif



