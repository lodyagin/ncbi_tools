/*   cn3dpane.h
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
* Revision 6.8  1998/12/16 22:49:38  ywang
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
PRK  prkSpecial; 
CharPtr title;
Boolean On;
} SpecialFeatureInfo, PNTR SpecialFeatureInfoPtr;

typedef ValNodePtr SpecialFeaturePtr;

extern DomainInfo **domaindata; 

extern GrouP LIBCALL DisplayControls PROTO((Nlm_GrouP prnt));  
extern GrouP LIBCALL ModelControls PROTO((Nlm_GrouP prnt));  
extern void LIBCALL ResetDisplayCtrls(void);
extern void LIBCALL ResetModelCtrls(void);
extern void LIBCALL Cn3D_CountDomainProc(void);
extern void LIBCALLBACK DoLinkSpecialFeatureWithMGD  PROTO((PFB pfbThis,Int4 iModel, Int4 iIndex, Pointer ptr));

#ifdef __cplusplus
}
#endif

#endif



