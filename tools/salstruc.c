/* ===========================================================================
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
* File Name:  salstruc.c
*
* Author:  Colombe Chappey
*
* Version Creation Date:   1/27/96
*
* $Log: salstruc.c,v $
* Revision 6.56  1999/05/03 22:15:37  chappey
* cs_option.slabel_format = PRINTID_FASTA_LONG
*
* Revision 6.55  1999/01/25 17:21:49  durand
* .
*
* Revision 6.54  1999/01/24 19:32:41  chappey
* MergeFunc
*
* Revision 6.53  1999/01/23 22:14:26  chappey
* bug in MergeFunc and related functions
*
* Revision 6.52  1999/01/19 03:26:18  chappey
* Merge function allows to specify the intervalle to merge
*
* Revision 6.51  1999/01/18 22:47:55  chappey
* merge update functions Merge3Func and Merge5Func
*
* Revision 6.50  1998/12/26 20:09:10  chappey
* .
*
* Revision 6.49  1998/12/02 07:00:38  chappey
* Copy Feature copies now all the names from the SeqFeat of the protein Bioseq
*
* Revision 6.48  1998/11/09 05:13:57  chappey
* renamed SeqAlign functions
*
* Revision 6.47  1998/10/29 19:07:06  chappey
* Feature propagation now extend CDS until the 1st stop codon
*
* Revision 6.46  1998/09/29 03:00:57  chappey
* Propagation of SEQFEAT_BOND, SEQFEAT_SITE
*
* Revision 6.45  1998/09/18 11:48:23  chappey
* SeqAlignGapCount and fix in BioseqTrimN
*
* Revision 6.44  1998/08/24 22:53:32  chappey
* .
*
* Revision 6.43  1998/07/29 01:16:20  chappey
* .
*
* Revision 6.42  1998/07/23 00:15:12  chappey
* bug in showfastagap_fromalign
*
* Revision 6.41  1998/07/22 23:51:03  chappey
* bug in export PHYLIP format
*
* Revision 6.40  1998/07/22 22:54:02  chappey
* .
*
* Revision 6.39  1998/07/17 19:07:57  chappey
* .
*
* Revision 6.38  1998/07/15 20:54:24  chappey
* seq_info
*
* Revision 6.37  1998/07/14 16:49:08  chappey
* .
*
* Revision 6.36  1998/06/12 03:08:34  chappey
* PropagateFeatureByApply only propagates new features
*
* Revision 6.35  1998/06/09 15:11:36  chappey
* .
*
* Revision 6.34  1998/06/02 05:16:49  chappey
* bug in multseqalign_from_pairseqalign
*
* Revision 6.33  1998/06/01 01:09:29  chappey
* Message -> ErrPostEx
*
* Revision 6.32  1998/05/23 04:27:40  chappey
* Blocks are visible only on sequences in list of Ids
*
* Revision 6.31  1998/05/21 17:38:23  chappey
* changed SEQFEAT to FEATDEF in PropagateFeatureBySeqLock
*
* Revision 6.30  1998/05/20 19:02:08  chappey
* fixes in PairSeqAlign2MultiSeqAlign
*
* Revision 6.29  1998/05/19 22:46:10  chappey
* fixes in data_collect_arrange for LOWER/UPPER cases
*
* Revision 6.28  1998/05/19 16:08:21  chappey
* Fixes in BioseqTrimN
*
* Revision 6.27  1998/05/19 02:36:03  chappey
* minor fixes
*
* Revision 6.26  1998/05/18 19:39:48  chappey
* fixes
*
* Revision 6.25  1998/05/12 21:23:30  chappey
* minor fixes
*
* Revision 6.24  1998/05/07 15:25:54  chappey
* ADD_REGION in ApplyBioFeatToSeqEntry
*
* Revision 6.23  1998/05/05 23:37:11  chappey
* Comments for CopySeqLocFromSeqAlign
*
* Revision 6.22  1998/05/05 04:30:42  chappey
* Changes in ApplyBioFeatToSeqEntry: do not re-translate the CDS because can not find the right Genetic code. Keep the 1rst translation (prot)
*
* Revision 6.21  1998/05/04 18:46:53  chappey
* Allows 1 base long feature
*
* Revision 6.20  1998/04/28 17:55:33  chappey
* BioseqTrimN new function to truncate nnn.s at sequence end
*
* Revision 6.19  1998/04/15 19:17:55  chappey
* not_aligned segments are written in LOWER cases in arrange_buffer
*
* Revision 6.18  1998/04/06 23:28:26  chappey
* minor change
*
* Revision 6.16  1998/04/02 23:42:58  chappey
* minor changes
*
* Revision 6.15  1998/04/01 18:26:09  chappey
* minor changes
*
* Revision 6.14  1998/03/23 22:36:36  chappey
* printf removed from showfastagap_fromalign
*
* Revision 6.13  1998/03/22 21:56:33  chappey
* new function for display alignment as fasta+gap showfastagap_fromalign
*
* Revision 6.12  1998/01/03 19:53:03  chappey
* bugs fixed in ShowAlignmentText
*
* Revision 6.11  1997/12/05 18:37:44  chappey
* casting variables
*
* Revision 6.10  1997/12/03 22:31:01  chappey
* no global variable
*
* Revision 6.9  1997/11/06 02:48:05  chappey
* Bugs fixed
*
* Revision 6.8  1997/10/21 23:02:47  chappey
* bugs fixed
*
* Revision 6.7  1997/10/21 22:49:01  chappey
* bugs fixed
*
* Revision 6.6  1997/10/02 05:58:07  chappey
* Bug fixes
*
* Revision 6.5  1997/09/18 14:36:43  kans
* final fixes to update sequence, removed vibrant dependencies (CC)
*
* Revision 6.4  1997/09/16 19:30:42  kans
* many fixes to speed up sequence merge/replace/copy features (CC)
*
* Revision 6.3  1997/09/05 15:04:09  kans
* ReplaceBioseq, Merge5Func, Merge3Func, and CpyFeatFunc (CC)
*
* Revision 6.2  1997/09/04 14:13:56  chappey
* bug fixes
*
* Revision 6.1  1997/08/26 21:00:32  kans
* various changes (CC)
*
* Revision 6.0  1997/08/25 18:54:13  madden
* Revision changed to 6.0
*
* Revision 5.67  1997/08/08 17:25:51  chappey
* bug fixes
*
* Revision 5.66  1997/07/14 12:38:25  kans
* needed to call MyGetTopSeqEntryForEntityID
*
* Revision 5.65  1997/07/14 04:30:14  chappey
* function FindSeqEntryForSeqId is now Callback fct in SeqEntryExplore
*
* Revision 5.64  1997/07/11 17:11:59  kans
* copied CheckSeqLocForPartial and SetSeqLocPartial temporarily
*
* Revision 5.63  1997/07/11 16:13:41  chappey
* move CopySeqLocFromSeqAlign from salfiles.c
*
* Revision 5.62  1997/07/08 05:46:33  chappey
* bug fixes
*
* Revision 5.61  1997/07/01 06:28:40  chappey
* fixes in FindSeqEntryFromId
*
* Revision 5.60  1997/06/26 20:22:47  chappey
* bug fixed in SeqAlignTranslate
*
* Revision 5.59  1997/06/24 01:13:49  chappey
* bugs fixed
*
* Revision 5.58  1997/06/02 19:20:41  kans
* remove vibrant dependency and unused variables
*
* Revision 5.57  1997/06/02 19:08:36  kans
* another round of changes
*
* Revision 5.56  1997/05/22 17:41:07  kans
* various fixes, supports paup format
*
* Revision 5.54  1997/05/12 17:47:43  kans
* various improvements
*
 * Revision 5.53  1997/05/04  21:39:40  kans
 * numerous fixes (CC)
 *
 * Revision 5.49  1997/04/03  06:32:18  chappey
 * *** empty log message ***
 *
 * Revision 5.48  1997/03/18  23:53:27  kans
 * *** empty log message ***
 *
 * Revision 5.47  1997/03/13  17:38:24  kans
 * *** empty log message ***
 *
 * Revision 5.35  1997/01/04  00:05:51  kans
 * *** empty log message ***
 *
 * Revision 5.26  1996/11/04  04:23:07  kans
 * in collect feature now gets interval on proper part
 *
 * Revision 5.23  1996/10/27  01:14:23  kans
 * *** empty log message ***
 *
 * Revision 5.22  1996/10/22  23:22:39  chappey
 * *** empty log message ***
 *
 * Revision 5.21  1996/10/04  19:13:20  kans
 * *** empty log message ***
 *
 * Revision 5.18  1996/09/12  02:38:07  chappey
 * New PropagateFeature function to propagate the features from one SeqEntry
 * change the SeqLocs of the features using given SeqAlign
 *
 * Revision 5.8  1996/07/15  15:20:19  kans
 * support propagation (preliminary)
 *
 * Revision 5.1  1996/05/30  17:29:01  kans
 * integrate single and multiple alignment code
 *
 * Revision 1.24  1996/05/14  16:46:15  kans
 * added calls to ExtractBioSourceAndPubs and ReplaceBioSourceAndPubs
 *
 * Revision 1.23  1996/05/06  19:16:36  kans
 * *** empty log message ***
 *
 * Revision 1.22  1996/05/02  21:33:07  kans
 * *** empty log message ***
 *
 * Revision 1.21  1996/04/29  21:59:57  vakatov
 * Sequince of length 1 now treated properly (infinitive loop eliminated)
 *
 * Revision 1.20  1996/04/25  18:59:21  vakatov
 * *** empty log message ***
 *
* ==========================================================================
*/
#include <salstruc.h>
#include <salsa.h>
#include <salutil.h>
#include <salsap.h>
#include <subutil.h>
#include <gather.h>
#include <satutil.h>
#include <sqnutils.h>

static Boolean stringhasnotext (CharPtr str)

{
  Char  ch;

  if (str != NULL) {
    ch = *str;
    while (ch != '\0') {
      if (ch > ' ' && ch <= '~') {
        return FALSE;
      }
      str++;
      ch = *str;
    }
  }
  return TRUE;
}

extern SelStructPtr BufferFree (SelStructPtr ssp)
{
  SelStructPtr   ssptmp, ssptmp2, next;
  SelEdStructPtr sep;
  ValNodePtr     vnp;

  for (ssptmp = ssp; ssptmp != NULL; ssptmp = ssptmp->next)
            ssptmp->region = NULL;
  ssptmp = ssp; 
  while (ssptmp != NULL)
  {
         next = ssptmp->next;
       ssptmp2 = ssptmp;
         if (ssptmp2->region != NULL)
         {
            sep = (SelEdStructPtr) ssptmp2->region;
            if (sep->region != NULL) 
               sep->region = SeqLocFree ((SeqLocPtr) sep->region);
            if (sep->data != NULL) {
               vnp = (ValNodePtr) sep->data;
               vnp->data.ptrvalue = NULL;
               sep->data = ValNodeFree (vnp);
            }
            ssptmp2->region = MemFree (sep);
         }
         ssptmp2 = MemFree (ssptmp2);
         ssptmp = next;
  }
  return NULL;
}

/*******************************************************************
***  
***    SetupDataBuffer
***    minbufferlength = (Int4) ( WINPERBUF * MAXCharLine * j );
*******************************************************************/
extern EditAlignDataPtr SetupDataBuffer (EditAlignDataPtr adp)
{
  if (adp==NULL) 
     return NULL;
  if ( adp->seqnumber == 0 ) return NULL;
  adp->minbufferlength = TMP_BUFFERLENGTH;
  adp->bufferlength = (Int4) MIN (adp->length, adp->minbufferlength);
  return adp;
}

/*******************************************************************
***  
***    SetupDataPanel
***
*******************************************************************/
extern EditAlignDataPtr SetupDataPanel (EditAlignDataPtr adp)
{
  Int4  j;
  Int4  lg;

  if (adp==NULL) 
     return NULL;
  if ( adp->seqnumber == 0 ) return NULL;
  lg = MAXLineWindow + 4;

  if ( adp->item_id != NULL ) adp->item_id = MemFree (adp->item_id);
  adp->item_id = NULL;
  adp->item_id = (Uint2Ptr)MemNew ((size_t) (lg * sizeof(Uint2)));
  for (j=0; j<MAXLineWindow; j++) adp->item_id[j] = 0;

  if ( adp->seqEntity_id != NULL ) adp->seqEntity_id = MemFree (adp->seqEntity_id);
  adp->seqEntity_id = NULL;
  adp->seqEntity_id =(Uint2Ptr)MemNew((size_t) (lg * sizeof(Uint2)) );
  for (j=0; j<MAXLineWindow; j++) adp->seqEntity_id[j] = 0;

  if ( adp->itemtype != NULL ) adp->itemtype = MemFree (adp->itemtype);
  adp->itemtype = NULL;
  adp->itemtype =(Uint2Ptr)MemNew ((size_t) (lg * sizeof(Uint2)));
  for (j=0; j<MAXLineWindow; j++) adp->itemtype[j] = 0;

  if ( adp->itemsubtype != NULL ) adp->itemsubtype = MemFree (adp->itemsubtype);
  adp->itemsubtype = NULL;
  adp->itemsubtype =(Uint2Ptr)MemNew((size_t) (lg * sizeof(Uint2)) );
  for (j=0; j<MAXLineWindow; j++) adp->itemsubtype[j] = 0;

  if ( adp->alignline != NULL ) adp->alignline = MemFree (adp->alignline);
  adp->alignline = NULL;
  adp->alignline= (Uint2Ptr)MemNew((size_t) (lg * sizeof(Uint2)));
  for (j=0; j<MAXLineWindow; j++) adp->alignline[j] = 0;

  if ( adp->colonne != NULL ) adp->colonne = MemFree (adp->colonne);
  adp->colonne = NULL;
  lg = adp->minbufferlength + adp->editbuffer + 4;
  adp->colonne = (Int4Ptr) MemNew ((size_t) (lg * sizeof(Int4)));
  for (j=0; j<adp->minbufferlength +adp->editbuffer; j++) adp->colonne[j] = -1;

  if (adp->item_id==NULL || adp->seqEntity_id==NULL || adp->itemsubtype==NULL 
  || adp->alignline==NULL || adp->colonne==NULL)  {
         adp->seqnumber = 0;
         return NULL;
  }
  return adp;
}

static Uint2 OBJ_ (Uint2 feattype)
{
  if ( feattype == FEATDEF_BAD ) return OBJ_BIOSEQ;
  return OBJ_SEQFEAT;
}

/***********************************************************
***
************************************************************/
typedef struct ccid {
  SeqIdPtr    sip;
  SeqEntryPtr sep;
  BioseqPtr   bsp;
} CcId, PNTR CcIdPtr;
 
typedef struct orgscan {
  ObjMgrPtr  omp;
  Int2       nuclCode;
  Int2       mitoCode;
  Boolean    mito;
  Char       taxname [64];
} OrgScan, PNTR OrgScanPtr;

static Boolean CC_OrgScanGatherFunc (GatherContextPtr gcp)

{
  BioSourcePtr   biop;
  ObjMgrTypePtr  omtp;
  OrgNamePtr     onp;
  OrgRefPtr      orp;
  OrgScanPtr     osp;
  ValNodePtr     sdp;
  SeqFeatPtr     sfp;
  Uint2          subtype;
  Int2           val;
  ValNodePtr     vnp;

  if (gcp == NULL || gcp->thisitem == NULL) return TRUE;
  if (gcp->thistype != OBJ_SEQFEAT  && gcp->thistype != OBJ_SEQDESC) return TRUE;

  osp = (OrgScanPtr) gcp->userdata;
  if (osp == NULL) return TRUE;

  subtype = 0;   
  if (gcp->thistype == OBJ_SEQFEAT  || gcp->thistype == OBJ_SEQDESC) {
    omtp = ObjMgrTypeFind (osp->omp, gcp->thistype, NULL, NULL);
    if (omtp == NULL) {
      return TRUE;
    }
    if (omtp->subtypefunc != NULL) {
      subtype = (*(omtp->subtypefunc)) (gcp->thisitem);
    }
  }  

  orp = NULL;
  biop = NULL;
  switch (gcp->thistype) {
    case OBJ_SEQFEAT :
      sfp = (SeqFeatPtr) gcp->thisitem;
      switch (subtype) {
        case FEATDEF_ORG :
          orp = (OrgRefPtr) sfp->data.value.ptrvalue;
          break;
        case FEATDEF_BIOSRC :
          biop = (BioSourcePtr) sfp->data.value.ptrvalue;
          break;
        default :
          break;
      }
      break;
    case OBJ_SEQDESC :
      sdp = (ValNodePtr) gcp->thisitem;
      switch (subtype) {
        case Seq_descr_modif :
          vnp = (ValNodePtr) sdp->data.ptrvalue;
          while (vnp != NULL) {
            val = (Int2) vnp->data.intvalue;
            if (val == MODIF_mitochondrial || val == MODIF_kinetoplast) {
              osp->mito = TRUE;
            }  
            vnp = vnp->next;
          }
          break;
        case Seq_descr_org :
          orp = (OrgRefPtr) sdp->data.ptrvalue;
          break;
        case Seq_descr_source :
          biop = (BioSourcePtr) sdp->data.ptrvalue;
          break;
        default :
          break;
      }
      break;
    default :
      break;
  }
 
  if (orp == NULL && biop != NULL) {
    orp = biop->org;
    osp->mito = (Boolean) (biop->genome == 4 || biop->genome == 5);
  }
  if (orp != NULL) {
    StringNCpy_0 (osp->taxname, orp->taxname, sizeof (osp->taxname));
    onp = orp->orgname;
    if (onp != NULL) {
      osp->nuclCode = onp->gcode;
      osp->mitoCode = onp->mgcode;
    }
  }
    
  return TRUE;
}

static Int2 CC_SeqEntryOrEntityIDToGeneticCode (SeqEntryPtr sep, Uint2 entityID, BoolPtr mito, CharPtr taxname, size_t maxsize)
{
  GatherScope  gs;
  OrgScan      osp;

  if (mito != NULL) {
    *mito = FALSE;
  }
  if (taxname != NULL && maxsize > 0) {
    *taxname = '\0';
  }
  osp.mito = FALSE;
  osp.nuclCode = 0;
  osp.mitoCode = 0;
  osp.omp = ObjMgrGet ();
  osp.taxname [0] = '\0';
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = TRUE;
  MemSet ((Pointer) (gs.ignore), (int)(TRUE), (size_t) (OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  gs.ignore[OBJ_SEQDESC] = FALSE;
  if (sep != NULL) {
    gs.scope = sep;
    GatherSeqEntry (sep, (Pointer) &osp, CC_OrgScanGatherFunc, &gs);
  } else if (entityID > 0) {
    GatherEntity (entityID, (Pointer) &osp, CC_OrgScanGatherFunc, &gs);
  }
  if (mito != NULL) {
    *mito = osp.mito;
  }
  if (taxname != NULL && maxsize > 0) {
    StringNCpy_0 (taxname, osp.taxname, maxsize);
  }
  if (osp.mito) {
    return osp.mitoCode;
  } else {
    return osp.nuclCode;
  }
}


static void FindSeqEntryForSeqIdCallback (SeqEntryPtr sep, Pointer mydata,
                                          Int4 index, Int2 indent)
{
  BioseqPtr          bsp;
  SeqIdPtr           sip;
  CcIdPtr            cip;
 
  if (sep != NULL && sep->data.ptrvalue && mydata != NULL) {
     cip = (CcIdPtr)mydata;
     if (cip->sep==NULL && IS_Bioseq(sep)) {
        bsp = (BioseqPtr) sep->data.ptrvalue;
        if (bsp!=NULL) {
           sip = SeqIdFindBest(bsp->id, 0);
           if (SeqIdForSameBioseq(cip->sip, sip)) {
              cip->sep = sep;
              cip->bsp = bsp;
           }
        }
     }   
  }
}
 
 
static Int2 CC_SeqEntryToGeneticCode (Uint2 entityID, SeqIdPtr sip)
{
  SeqEntryPtr sep_head,
              sep;
  CcId        ci;
  Int2        genCode = 0;

  sep_head  = GetTopSeqEntryForEntityID (entityID);
  ci.sip = SeqIdDup (sip);
  ci.sep = NULL;
  SeqEntryExplore(sep_head,(Pointer)&ci, FindSeqEntryForSeqIdCallback);
  sep = ci.sep;
  SeqIdFree (ci.sip);
  if (sep!=NULL) {
/*
     genCode = SeqEntryToGeneticCode (sep, NULL, NULL, 0);*/
     genCode = CC_SeqEntryOrEntityIDToGeneticCode (sep, 0, NULL, NULL, 0);
  }
  return genCode;
}

/************************************************** 
***  next_notemptyline: search for not empty line (not_empty_string) 
***      if not exists, go to next alignment line  
***  get_substring
***
***************************************************/
static CharPtr get_substring (CharPtr str, Int4 drw_start, Int4 drw_width)
{
  Int4            width;
  Int4            stringlens;
  CharPtr         strp;
  if (str == NULL ) return NULL; 
  stringlens = StringLen (str);
  if ( drw_start >= stringlens ) { return NULL; } 
  strp = str + drw_start;
  stringlens = StringLen (strp);
  if (stringlens == 0) return NULL; 
  width = MIN ((Int2) drw_width, (Int2) stringlens);
  if ( !not_empty_string (strp, width) ) return NULL;
  return strp;
}

extern CharPtr next_notemptyline (ValNodePtr anp_list, ValNodePtr linebuff, Int2 numberalignline, Int2 *index, Int4 start, Int4 *drw_width, TextAlignBufPtr *tdp, AlignNodePtr *anp)
{
  TextAlignBufPtr tdptmp;
  ValNodePtr      vnp;
  ValNodePtr      vnpanp;
  CharPtr         str;
  Int4            width = *drw_width;
  Int2            j = 1;

  if (*index > numberalignline) return NULL;
  vnp = linebuff;
  tdptmp = (TextAlignBufPtr) vnp->data.ptrvalue;
  vnpanp = anp_list;
  if (j < *index ) {
         vnp = vnp->next;
         while (j < *index) {
           tdptmp = (TextAlignBufPtr) vnp->data.ptrvalue;
           if (OBJ_(tdptmp->feattype) == OBJ_BIOSEQ) vnpanp = vnpanp->next;
           j++;
           if (j == *index ) break;
           vnp = vnp->next;
         }
  }
  while ( vnp != NULL && j <= numberalignline ) 
  {
         str = get_substring (tdptmp->buf, start, width);
         if ( str != NULL ) break;
         vnp = vnp->next;
         if (vnp == NULL) break;
         tdptmp = (TextAlignBufPtr) vnp->data.ptrvalue;
         if (OBJ_(tdptmp->feattype) == OBJ_BIOSEQ) vnpanp = vnpanp->next;
         j++;
  } 
  if ( j > numberalignline && vnp == NULL ) return NULL;
  *index = j;
  *drw_width = width;
  *tdp = tdptmp;
  if (OBJ_(tdptmp->feattype) == OBJ_BIOSEQ) 
     *anp= (AlignNodePtr) vnpanp->data.ptrvalue;
  else 
     *anp = NULL;
  return str;
}

/******************************************************
***  GetPerCol
***    .seqline    : index of the sequences  (0, 1..)
***    .alignline  : index of alignment line (0, 1..)  
***    
*******************************************************/
static void GetPerCol (EditAlignDataPtr adp, Int4 hoffset)
{
  SeqAlignPtr  salp;
  CompSegPtr   dsp;
  Int4Ptr      lenp;
  AlignNodePtr anp;
  AlignSegPtr  segs;
  Int4         j = 0; 
  Int4         m;
  Int4         k, l;

  if (adp->minbufferlength == 0) return;
  if (adp->sap_align != NULL && (salp = (SeqAlignPtr)adp->sap_align->data) != NULL) {
         j = 0;
         while ( j < adp->bufferlength && salp != NULL )
         {
                dsp = (CompSegPtr) salp->segs;
                lenp = dsp->lens;
                k = m = 0; 
                while (k < dsp->numseg && m < hoffset) 
                {
                    for (l = 0; l < *lenp && m < hoffset; l++, m++) continue;
                    if (m == hoffset) break;
                    lenp++;
                    k++;
                }
                for ( ; k <= dsp->numseg && j < adp->bufferlength; k++, lenp++)
                {
                   for (l=0; l < *lenp && j < adp->bufferlength; l++, m++, j++)
                   {
                      adp->colonne[j] = m;
                   }
                }
                if ( j < adp->bufferlength && salp->next != NULL ) 
                {
                   for (l=0; l < adp->intersalpwidth && j < adp->bufferlength; l++, j++)
                      adp->colonne[j] = -1;
                   salp = salp->next;
                }
                else break;
         }
  }
  else if ( adp->anp_list != NULL ) 
  {
         anp = (AlignNodePtr) adp->anp_list->data.ptrvalue;
         segs = anp->segs;
         j=0;
         for (k=0, m=0; segs != NULL && m < hoffset; k++, segs=segs->next)
            for (l=0; l<(segs->gr.right-segs->gr.left+1) && m<hoffset; l++, m++)
                continue;
         for (; segs != NULL && j < adp->bufferlength; k++, segs=segs->next)
            for(l=0; l<(segs->gr.right-segs->gr.left+1) && j< adp->bufferlength; l++, m++, j
++)
                adp->colonne[j] = m;
  }
  else {
         return;
  }
  if ( j < adp->minbufferlength + adp->editbuffer && adp->colonne[j-1] > -1) {
         adp->colonne[j] = adp->colonne[j-1]+1;
         j++;
  }
  if ( j < adp->minbufferlength + adp->editbuffer )
         for (; j < adp->minbufferlength + adp->editbuffer; j++)
                adp->colonne[j] = -1;
  return ;
}

/*********************************************************************
***  is_feature_to_buffer
*********************************************************************/
extern SelEdStructPtr is_feature_to_buffer (ValNodePtr vnphead, Uint2 bspitemID, Uint2 entityID, Int4 from, Int4 drw_width, SeqAlignPtr salp, Uint2 seqedit, ValNodePtr sqloc_list)
{
  SelEdStructPtr  cds;
  SeqIdPtr        sip;
  SeqLocPtr       slp;
  Int4            start, stop;
  Int4            start2, stop2;
  Int2            chklocp;
  Boolean         draw = FALSE;

  if (vnphead == NULL) return NULL;
  cds = (SelEdStructPtr) vnphead->data.ptrvalue;
  if ( cds == NULL ) 
     return NULL;
  sip = SeqLocId ((SeqLocPtr) cds->region);
  for (; cds != NULL && !draw; cds = cds->next) 
  {
     if (entityID == cds->entityID && bspitemID == cds->bsp_itemID)
     {
        if ( cds->regiontype == OM_REGION_SEQLOC && cds->region != NULL ) 
        {
           slp = (SeqLocPtr) cds->region;
           start2 = SeqLocStart(slp);
           chklocp =chkloc(sip, start2, sqloc_list, &start2);
           start= SeqCoordToAlignCoord(start2, sip, salp, 0, chklocp);
           stop2 = SeqLocStop(slp);
           chklocp =chkloc(sip, stop2, sqloc_list, &stop2);
           stop = SeqCoordToAlignCoord(stop2, sip, salp, 0, chklocp);
           if (start <= stop && start < from + drw_width && stop >= from)  {
                 draw = TRUE;
                 break;
           }
        }
     }
  }
  if (draw) 
     return cds;
  return NULL;
}

/******************************************************************/
extern ByteStorePtr cds_to_pept (SeqLocPtr slp, Uint1 frame, Int2 gencode, Boolean include_stop)
{
  ByteStorePtr  bs;
  ValNodePtr    code;
  CdRegionPtr   crp;
  SeqFeatPtr    sfp;
  ValNodePtr    vnp;

  if (slp == NULL) return NULL;
  sfp = SeqFeatNew ();
  if (sfp == NULL) return NULL;
  sfp->data.choice = SEQFEAT_CDREGION;
  crp = CdRegionNew ();
  sfp->data.value.ptrvalue = (Pointer) crp;
  if (crp == NULL) {
         SeqFeatFree (sfp);
         return NULL;
  }
  crp->orf = FALSE;
  crp->conflict = FALSE;
  crp->frame = frame;
  crp->gaps = 0;
  crp->mismatch = 0;
  crp->stops = 0;
  code = ValNodeNew (NULL);
  if (code != NULL) {
         code->choice = 254;
         vnp = ValNodeNew (NULL);
         code->data.ptrvalue = vnp;
         if (vnp != NULL) {
            vnp->choice = 2;
            vnp->data.intvalue = (Int4) gencode;
         }
  }
  crp->genetic_code = code;
  crp->code_break = NULL;
  sfp->location = slp;
  bs = ProteinFromCdRegion (sfp, include_stop);
  if (bs == NULL) {
     sfp->location = NULL;
     SeqFeatFree (sfp);
     return NULL;
  }
  return bs;
}

/******************************************************************/
static CharPtr SeqAlignTranslate (SeqAlignPtr salp, Uint2 entityID, Int4 from, Int4 to, Uint1 codonbase, SeqIdPtr sip, Int2 seqnumber, Uint1 strand, ValNodePtr sqlocs)
{
  SeqLocPtr    slp;
  SeqIntPtr    sit;
  ByteStorePtr bsp;
  CharPtr      pep = NULL, 
               pepPtr = NULL;
  CharPtr      str = NULL, 
               strPtr = NULL;
  CharPtr      buffer = NULL, 
               bufferPtr = NULL;
  Int4         k, 
               slplens,
               strlens;
  Int4         start, stop,
               len;
  Int2         cb;  
  Int2         genCode;
  BioseqPtr    bsq;

  if ( salp == NULL ) return NULL; 
  if ( salp->segtype != COMPSEG ) {
     return NULL; 
  }
  start= (Int4) AlignCoordToSeqCoord (from, sip, salp, sqlocs, 0);
  stop = (Int4) AlignCoordToSeqCoord (to, sip, salp, sqlocs, 0);

  bsq = BioseqLockById (sip);
  if (bsq != NULL) {
     if (start + stop >= bsq->length) 
        stop = bsq->length - 1;
     BioseqUnlock (bsq);
  }
  else return NULL;

  slp = fuzz_loc (start, stop, strand, sip, TRUE, TRUE);
  if (slp == NULL) {
     return NULL; 
  }
  slplens = SeqLocLen (slp);
  if ( slplens < (Int4) (3 + codonbase) ) {
     return NULL; 
  }
  genCode = CC_SeqEntryToGeneticCode (entityID, sip);
  if (genCode == 0)
     genCode = Seq_code_ncbieaa; 

  sit = (SeqIntPtr) slp->data.ptrvalue;
  if (strand == Seq_strand_plus) 
  {
     sit->from = sit->from + codonbase;
     slplens = SeqLocLen (slp);
     cb = (Int2)(slplens % (Int4) 3); 
     if (cb == 1) {
          sit->to --;
     } 
  }
  else if (strand == Seq_strand_minus) {
     sit->to = sit->to - codonbase;
     slplens = SeqLocLen (slp);
     cb = (Int2)(slplens % (Int4) 3); 
     if (cb == 1 && sit->from >0) {
          sit->from --;
     } else if (cb == 2) {
          sit->from ++;
     } 
     if (cb == 0) codonbase = 0;
     else if (cb == 1) codonbase = 1;
     else if (cb == 2) codonbase = 2;
  }
  slplens = SeqLocLen (slp);
  if ( slplens >= 3 ) {
     bsp = cds_to_pept (slp, 1, genCode, TRUE);
     str = (CharPtr) BSMerge (bsp, NULL);
     BSFree (bsp);
     pep = MemNew ((size_t) ((slplens +5) *sizeof(Char)));
     pep = emptystring (pep, (Int4)(slplens + 5));
     pep [slplens +3] = '\0';
     pepPtr = pep;
     *pepPtr = ' ';
     pepPtr += codonbase +1;
     strlens = 3*StringLen(str);
     if (slplens < strlens) {
        strlens=(Int4)(slplens/(Int4)3);
        str [strlens] ='\0';
     }
     if (strand == Seq_strand_minus)
          reverse_string (str);
     strlens = StringLen(str);
     strPtr = str;
     for (k = 0; k < strlens; k++, pepPtr += 3, strPtr++) {
          *pepPtr = *strPtr; 
     }
     MemFree (str);
     buffer = MemNew ((size_t) ((to -from +5) *sizeof(Char)));
     buffer = emptystring (buffer, (Int4)(to -from +5));
     buffer [to -from +3] = '\0';
     bufferPtr = buffer;
     *bufferPtr = ' ';
     buffer = ReadBufferFromSap (pep, buffer, salp, sip, from, to, &len);
     MemFree (pep);
  }
  SeqLocFree (slp);
  return buffer; 
}

static CharPtr prot_to_putprot (CharPtr trans)
{
  CharPtr strp;
  Boolean in = TRUE;

  for (strp = trans; *strp != '\0'; strp++)
  {
         if (*strp == 'M') in = TRUE;
         else if (*strp == '*') {
            if (!in) *strp = ' ';
            else in = FALSE;
         }
         else if (in) *strp = '~';
         else *strp = ' ';
  }
  return trans;
}

static CharPtr prot_to_rputprot (CharPtr trans)
{
  CharPtr strp, strptmp;
  Boolean in = TRUE;
  Int4    j, k;

  strptmp = strp = trans; 
  while (*strptmp != '\0') {
     for (j=0; *strptmp != '\0'; strptmp++, j++) {
         if (*strptmp == '*') {
            for (k=0; k<=j && *strp != '\0'; strp++, k++) {
               if (*strp != '*' && *strp != 'M' ) *strp = ' ';
            }
            in = TRUE; 
            break;
         }
         if ( *strptmp == 'M' ) {
            if (in) {
               for (k=0; k<=j && *strp != '\0'; strp++, k++) {
                  if (*strp != 'M' && *strp != '*') *strp = '~';
               }
               in = FALSE; 
            }
            else {
               for (k=0; k<j && *strp != '\0'; strp++, k++) {
                  if (*strp != 'M' && *strp != '*') *strp = ' ';
               }
            }
            break;
         }
     }
     if (strptmp == '\0') break;
     strptmp++;
     j++;
  }
  if (in) {
     for (k=0; k<j && *strp != '\0'; strp++, k++) 
         if (*strp != 'M' && *strp != '*') *strp = '~';
  }
  else {
     for (k=0; k<j && *strp != '\0'; strp++, k++) 
         if (*strp != 'M' && *strp != '*') *strp = ' ';
  }
  strptmp = strp = trans; 
/*
  strp++;
  for(; *strp != '\0'; strp++, strptmp++)
     if (*strptmp == '*' && *strp != '~' ) *strptmp = ' ';
*/
  return trans;
}

static SelStructPtr MakeRf (SelStructPtr buffer, Uint2 entityID, Uint2 itemID, Uint2 feattype, SeqIdPtr bspsip, Uint1 strand, Uint1 j, SeqAlignPtr salp, Int4 fromseq, Int4 toseq, EditAlignDataPtr adp, Int4 from_bufferstart, Int4 to_bufferstart, Uint2 typerf)
{
  SelEdStructPtr  rf;
  CharPtr         trans;

  trans = (CharPtr) SeqAlignTranslate (salp, entityID, fromseq, toseq, (Uint1)(j), bspsip, adp->seqnumber, strand, adp->sqloc_list);
  if (trans != NULL) {
     if ( adp->prot_mode == PUTPROT ) {
        if (strand == Seq_strand_minus) trans = prot_to_rputprot (trans);
        else trans = prot_to_putprot (trans); 
     }
     rf=new_seledstruct(entityID, itemID, feattype, itemID, from_bufferstart, to_bufferstart,  bspsip, strand, TRUE, NULL, (Pointer)trans, 0, 1);
     addssp (&(buffer), typerf, (Pointer) rf, itemID);
  }
  return buffer;
}


static SelStructPtr get_firstline (SelEdStructPtr sesp1, SelStructPtr buffer)
{
  SelEdStructPtr sesp;
  SelStructPtr   buf;

  if (buffer == NULL) {
         return NULL;
  }
  if (sesp1 == NULL) {
         return buffer;
  }
  for (buf=buffer; buf!=NULL; buf=buf->next) 
  {
         sesp = (SelEdStructPtr) buf->region;
         if (is_sameId (sesp1->entityID, sesp1->itemID, sesp1->itemtype, 255, sesp->entityID, sesp->itemID, sesp->itemtype, 255) ) 
            break;
  }
  if (buf == NULL) {
         return buffer;
  }
  return buf;
}

static Boolean has_complement (ValNodePtr params, Uint2 entityID, Uint2 itemID)
{
  ValNodePtr  vnp;
  SeqParamPtr prm;

  for (vnp = params; vnp != NULL; vnp = vnp->next)
  {
     prm = (SeqParamPtr) vnp->data.ptrvalue;
     if (prm->entityID == entityID && prm->itemID == itemID) {
        if ( prm->complement ) return TRUE;
     }
  }
  return FALSE;
}

static Boolean rf_on (ValNodePtr params, Uint2 entityID, Uint2 itemID, Uint2 rf)
{
  ValNodePtr  vnp;
  SeqParamPtr prm;

  for (vnp = params; vnp != NULL; vnp = vnp->next)
  {
     prm = (SeqParamPtr) vnp->data.ptrvalue;
     if (prm->entityID == entityID && prm->itemID == itemID) {
        return prm->rf[rf];
     }
  }
  return FALSE;
}

/*********************************************************************
***  arrange_buffer
*********************************************************************/
static void arrange_buffer (EditAlignDataPtr adp)
{
  TextAlignBufPtr tdp;
  AlignNodePtr    anp;
  SelEdStructPtr  ssp = NULL, 
                  bspssp = NULL;
  SelEdStructPtr  firstsesp = NULL;
  SeqIdPtr        bspsip;
  CharPtr         bufstr = NULL;
  Int4            fromseq = adp->bufferstart;
  Int4            toseq = adp->bufferstart + adp->bufferlength -1;
  Int4            from = 0;
  Int4            to = adp->bufferlength-1;
  Int4            length = adp->bufferlength;
  Int2            index = 1;
  Uint2           itemID, entityID;
  Int2            j, k;

  SeqAlignPtr     salp = (SeqAlignPtr) adp->sap_align->data;
  SelEdStructPtr  drawfeat = NULL;
  ValNodePtr      vnpfeat;

  if (adp->firstssp != NULL)
     firstsesp = SelEdStructDup ((SelEdStructPtr) adp->firstssp->region);
  BufferFree (adp->buffer);
  adp->buffer = NULL;
  if (adp->draw_scale) {
         ssp =new_seledstruct(LINE0,LINE0,LINE0,LINE0, from, to, NULL, 0, FALSE, NULL, NULL, 0, 1);
         addssp (&(adp->buffer), EDITDEF_SCA, (Pointer) ssp, EDITDEF_SCA);
  }
  if (adp->draw_bars) {
         ssp =new_seledstruct(LINE0,LINE0,LINE0,LINE0, from, to, NULL, 0, FALSE, NULL, NULL, 0, 1);
         addssp (&(adp->buffer), EDITDEF_SCB, (Pointer) ssp, EDITDEF_SCB);
  }
  bufstr =next_notemptyline (adp->anp_list, adp->linebuff, adp->numberalignline, &index, from, &length, &tdp, &anp);
  if ( index > adp->numberalignline ) {
         return;
  }
  while ( index <= adp->numberalignline && bufstr != NULL) 
  {
         if ( OBJ_(tdp->feattype) == OBJ_BIOSEQ ) {
            itemID = anp->bsp_itemID;
            if ( adp->input_format == OBJ_BIOSEQ ) {
               entityID = tdp->seqEntityID;
            }
            else if ( adp->input_format == OBJ_SEQALIGN ) {
               entityID = anp->seq_entityID;
            }
            bspsip = anp->sip;
/**
WriteLog ("ARRANGE %d %d   %d %d %d %d  %d %d %d %d  %d  \n", entityID, itemID, anp->entityID, anp->itemID, anp->seq_entityID, anp->bsp_itemID, tdp->entityID, tdp->itemID, tdp->seqEntityID, tdp->bsp_itemID, OBJ_(tdp->feattype));
**/
         } 
         else {
            itemID = tdp->itemID;
            entityID = tdp->entityID;
/**
WriteLog ("ARRANGFEATE %d %d   %d %d %d %d  %d   %d  \n", entityID, itemID, tdp->entityID, tdp->itemID, tdp->seqEntityID, tdp->bsp_itemID, tdp->feattype, OBJ_(tdp->feattype));
**/
         }
         if (OBJ_(tdp->feattype) == OBJ_BIOSEQ && (is_seqvisible(entityID, itemID, adp->seq_info)) )
         {
            bspssp = new_seledstruct (entityID, itemID, OBJ_(tdp->feattype), 
                 itemID, from + adp->bufferstart, to + adp->bufferstart,  
                 bspsip, Seq_strand_plus, TRUE, NULL, (Pointer) tdp, 0, 1);
            addssp (&(adp->buffer), tdp->feattype, (Pointer) bspssp, itemID);

            if (has_complement (adp->params, bspssp->entityID, itemID)) {
                addssp (&(adp->buffer), EDITDEF_CPL, (Pointer) bspssp, itemID);
            }
            for (j=0; j<3;j++) {
               if (rf_on (adp->params, bspssp->entityID, itemID, j))  
                  MakeRf (adp->buffer, entityID, itemID, OBJ_(tdp->feattype), bspsip, Seq_strand_plus, (Uint1)j, salp, fromseq, toseq, adp, from + adp->bufferstart, to + adp->bufferstart, (Uint2)(EDITDEF_RF1 + j));
            }
            for (j=3, k=2; j<6;j++, k--) {
               if (rf_on (adp->params, bspssp->entityID, itemID, j))  
                  MakeRf (adp->buffer, entityID, itemID, OBJ_(tdp->feattype), bspsip, Seq_strand_minus, (Uint1)k, salp, fromseq, toseq, adp, from + adp->bufferstart, to + adp->bufferstart, (Uint2)(EDITDEF_RF4 + k));
            }
            if ( adp->seqfeat != NULL ) 
            {
                vnpfeat = adp->seqfeat;
                for (; vnpfeat != NULL; vnpfeat = vnpfeat->next) 
                {
                   drawfeat = is_feature_to_buffer (vnpfeat, itemID, entityID, fromseq, length, salp, adp->input_format, adp->sqloc_list);
                   if (drawfeat != NULL) 
                   {
                      ssp = drawfeat;
                      addssp (&(adp->buffer), (Uint1)vnpfeat->choice, (Pointer) ssp, ssp->entityID);
                      if ( vnpfeat->choice == FEATDEF_CDS ) 
                      {
                         if ( drawfeat->data != NULL ) {
                            ssp = drawfeat;
                            addssp (&(adp->buffer), FEATDEF_TRSL,(Pointer) ssp, ssp->entityID);
                         }
                         tdp = TextAlignBufFind (adp->linebuff, drawfeat->entityID, drawfeat->itemID, drawfeat->itemtype);
                         if ( tdp != NULL )
                         {
                            ssp = new_seledstruct(tdp->entityID, tdp->itemID,  OBJ_SEQFEAT, itemID, from + adp->bufferstart, to + adp->bufferstart, bspsip, Seq_strand_plus, TRUE, NULL, (Pointer) tdp, 0, 1);
                            addssp(&(adp->buffer), FEATDEF_PROT, (Pointer) ssp, tdp->itemID);
                         }
                      }
                   }
                }
            }
            if ( adp->feat != NULL ) {
                vnpfeat = adp->feat;
                for (; vnpfeat != NULL; vnpfeat = vnpfeat->next)
                {
                   drawfeat = is_feature_to_buffer(vnpfeat, itemID, entityID, fromseq, length, salp, adp->input_format, adp->sqloc_list);
                   if (drawfeat != NULL) 
                   {
                      ssp = drawfeat;
                      addssp (&(adp->buffer), (Uint1)vnpfeat->choice, (Pointer) ssp, ssp->entityID);
                      if (vnpfeat->choice == SEQFEAT_CDREGION 
                      && drawfeat->data != NULL) {
                         ssp = drawfeat;
                         addssp (&(adp->buffer), FEATDEF_TRSL, (Pointer) ssp, ssp->entityID);
                      }
                   }
                }
            }
         }
         index++;
         length = adp->bufferlength;
         bufstr = next_notemptyline (adp->anp_list, adp->linebuff, 
                     adp->numberalignline, &index, from, &length, &tdp, &anp);
  }
  adp->buffertail = adp->buffer;
  if (adp->buffertail !=NULL)
     for (; adp->buffertail->next !=NULL; adp->buffertail =adp->buffertail->next)
        continue;
  adp->firstssp = get_firstline (firstsesp, adp->buffer);
  SelEdStructDel (firstsesp);
  return;
}


static Boolean update_fromalignnode (EditAlignDataPtr adp)
{
  AlignNodePtr  anp;
  ValNodePtr    curr;      /*for the list of AlignNodePtr*/
  ValNodePtr    list;      /*list of DrawText*/
  ValNodePtr    tdp_list;  /*list of DrawText*/
  TextAlignBufPtr tdp;
  Int4          start;
  Int4          p_stop=0;
  Uint2         entityID;
  SeqIdPtr      sip;
  ValNodePtr    vnp;
  SeqParamPtr   prm;
  SelEdStructPtr seq_info;
  Boolean       first_draw;
  Int2 j;

  if (adp->minbufferlength == 0) return FALSE;
  if (adp->linebuff != NULL) 
         adp->linebuff = (ValNodePtr) FreeTextAlignList(adp->linebuff);
  adp->linebuff = NULL;
  adp->numberalignline = 0;
  if ( adp->anp_list == NULL ) {
     ErrPostEx (SEV_ERROR, 0, 0, "fail in update_fromalignnode [1]");
  }
{
/***!!!!!!!!!!!!!!**/
SelEdStructPtr tmp;
  for (tmp=adp->seq_info; tmp!=NULL; tmp=tmp->next)
  {
     if (tmp->entityID==0)
        break;
  }
  first_draw=(Boolean)(tmp!=NULL);
}
  vnp = adp->params;
  for (curr = adp->anp_list, j=0; curr != NULL; curr = curr->next, j++)
  {
         anp = (AlignNodePtr) curr->data.ptrvalue;
         if ( anp == NULL ) {
                break;
         }
         start = (Int4) (adp->gr.left +adp->bufferstart);
         list = (ValNodePtr) ProcessTextAlignNode (anp, start, 
                start + adp->bufferlength + adp->editbuffer, &p_stop, NULL, 
                adp->visibleWidth, (Int1)0, (Uint4)0, NULL);  
         if ( list == NULL ) {
                break;
         }
         tdp_list = list;

         while (tdp_list != NULL)
	 {
                tdp = (TextAlignBufPtr) tdp_list->data.ptrvalue;
                if (adp->input_format == OBJ_BIOSEQ) {
                   entityID = tdp->seqEntityID;
                } else if ( adp->input_format == OBJ_SEQALIGN ) {
                   entityID = anp->seq_entityID;
                }
                if (adp->master.entityID == 0)
                {
                   if (adp->master.region != NULL) {
                      sip = SeqLocId((SeqLocPtr)adp->master.region);
                      if (SeqIdForSameBioseq(sip, anp->sip)) 
                      {
                         adp->master.entityID = entityID;
                         adp->master.itemID = anp->bsp_itemID;
                         adp->master.itemtype = OBJ_BIOSEQ;
                         adp->caret.entityID = entityID;
                         adp->caret.itemID = anp->bsp_itemID;
                         adp->caret.itemtype = OBJ_BIOSEQ;
                      }
                   }
                }
                if (vnp != NULL) {
                   prm = (SeqParamPtr) vnp->data.ptrvalue;
                   if (prm->entityID == 0) {
                      prm->entityID = entityID;
                      prm->itemID = anp->bsp_itemID;
                   }
                }
                if(tdp->buf != NULL)
                {
                       if (tdp->label == NULL) {
                          tdp->label=(CharPtr)MemNew((size_t)(64*sizeof(Char)));
                          tdp->label = StringSave ("unknown");
                       }
                       adp->numberalignline++;
                       ValNodeAddPointer (&adp->linebuff, 0, (Pointer) tdp);
                }
                if (first_draw)
                {
                   for (seq_info=adp->seq_info;seq_info!=NULL;seq_info=seq_info->next)
                   {
                      sip=SeqLocId((SeqLocPtr)seq_info->region);
                      if (SeqIdForSameBioseq(sip, anp->sip))
                      {
                         seq_info->entityID=entityID;
                         seq_info->itemID=anp->bsp_itemID;
                         seq_info->itemtype=OBJ_BIOSEQ;
                         break;
                      }
                   }
                }
                tdp_list = tdp_list->next;
         }
         if (vnp != NULL) vnp = vnp->next;
  }
  if (adp->numberalignline == 0)
  { 
         return FALSE;
  } 
  GetPerCol (adp, adp->bufferstart);
  return TRUE;
}

static Int4 get_bufferstart (EditAlignDataPtr adp)
{
  Int4 start, modstart = 0;

  if ( adp->minbufferlength == 0 ) return 0;
  if ( adp->hoffset < adp->minbufferlength/3 ) { 
         start = 0;
  }
  else if ( adp->length < adp->minbufferlength ) {
         start = 0;
  }
  else if ( adp->hoffset > adp->length- adp->minbufferlength) {
         start = adp->length - adp->minbufferlength;
  }
  else {
         start = adp->hoffset-(adp->minbufferlength/3);
  }
  if (start > 0) {
         modstart = start % adp->visibleWidth; 
  }
  if (modstart > 0) {
         start -= modstart;
  }
  return start;
}

static void count_feature_buf_line (ValNodePtr fnp_list, Int4 g_left, Int4 g_right, ValNodePtr PNTR feature_line)
{
        FeatNodePtr fnp;
        Int4 c_left, c_right;
        ValNodePtr vnp;
        Boolean found;

        if(fnp_list == NULL)
                return;
        
        while(fnp_list)
        {
           fnp = fnp_list->data.ptrvalue;
           c_left = fnp->extremes.left;
           c_right = fnp->extremes.right;
           if(!(c_left > g_right || c_right < g_left))
           {
                found = FALSE;
                for(vnp = *feature_line; vnp != NULL && !found; vnp = vnp->next)
                {
                        if(vnp->data.intvalue == (Int4)(fnp->itemID))
                                found = TRUE;
                }
                if(!found)
                        ValNodeAddInt(feature_line, 0, (Int4)(fnp->itemID));
           }
           fnp_list = fnp_list->next;
        }
}

static Int4 CountTextAlignNodeNum(AlignNodePtr anp, Int4 m_left, Int4 m_right)
{
        Int4 num_line = 0;
        Int4 g_left, g_right;

        AlignSegPtr asp;
        ValNodePtr feature_line, curr;  /*the number of lines for a feature*/

        g_left = anp->extremes.left;
        g_right = anp->extremes.right;
        if(m_left > g_right || m_right < g_left)
                return 0;

        num_line = 1;
        feature_line = NULL;

        /*process  the GAPs and the DIAGs segs*/
        for(asp = anp->segs; asp !=NULL; asp = asp->next)
        {
           g_left = asp->gr.left;
           g_right = asp->gr.right;
           if(!(g_left > m_right || g_right < m_left))
           {
              switch(asp->type)
              {  
                case GAP_SEG:
                   break;

                case REG_SEG:
                case DIAG_SEG:
                   g_left = MAX(m_left, g_left);
                   g_right = MIN(m_right, g_right);
                   count_feature_buf_line (asp->cnp, g_left, g_right, &feature_line);
                   break;
                default:
                   break;
              }
           }
           if(g_left > m_right)
                break;
        }
        if(feature_line != NULL)
        {
           for(curr = feature_line; curr != NULL; curr = curr->next)
                ++num_line;
           ValNodeFree(feature_line);
        }
 
        return num_line;
}

static Int4 count_total_line(ValNodePtr anp_list, Int4 line_len, Int4 left, Int4 right, Int4Ptr h_blocks)
{
        AlignNodePtr anp;
        Int4 c_start, c_stop;
        Int4 line_num = 0;
        ValNodePtr curr;
         
        if(anp_list == NULL)
                return line_num;
        if(h_blocks != NULL)
                *h_blocks = 0;
        anp = anp_list->data.ptrvalue;
        if(left == -1)
                left = anp->extremes.left;
        if(right == -1)
                right = anp->extremes.right;
        if(left > anp->extremes.right || right < anp->extremes.left)
                return line_num;
        left = MAX(left, anp->extremes.left);
        right = MIN(right, anp->extremes.right);
        if (left >= right)
                return line_num;
        c_start = left;
        while(c_start <=right)
        {
           c_stop = MIN(right, (c_start+line_len-1));
           for(curr = anp_list; curr != NULL; curr = curr->next)
           {  
              anp = curr->data.ptrvalue;
              line_num += CountTextAlignNodeNum(anp, c_start, c_stop);
           }
           if(h_blocks != NULL)
              ++(*h_blocks);
           c_start = c_stop+1;
        }
  return line_num;
}
 
static Int4 addline_perblock (EditAlignDataPtr adp, Int4 diffs)
{
  Int4        line = 0;
  ValNodePtr  vnp;
  SeqParamPtr prm;
  Int1        j;

  if (adp->draw_scale) line += (Int4) diffs;
  if (adp->draw_bars)  line += (Int4) diffs;
  for (vnp = adp->params; vnp != NULL; vnp = vnp->next)
  {
     prm = (SeqParamPtr) vnp->data.ptrvalue;
     if ( prm->complement ) line += (Int4) diffs;
     for (j=0; j<=6; j++) 
        if (prm->rf[j]) line += (Int4) diffs;
  }
  return line;
}  
 
static Int2 feat_linenum (Int4 slp_start, Int4 slp_stop, Int4 line_len, Int4 left,
Int4 right)
{
  Int4         modstart;
  Int4         modstop;
 
  slp_start = MAX (slp_start, left);
  modstart = slp_start % line_len;
  if ( modstart > 0) slp_start -= modstart;
 
  slp_stop = MIN (slp_stop, right);
  modstop = slp_stop % line_len;
  if ( modstop > 0) slp_stop += line_len;
 
  return (Int2)((slp_stop - slp_start) / line_len);
}
 
static Int4 CountFeatNum (ValNodePtr adpfeat, Int4 line_len, Int4 left, Int4 right){
  ValNodePtr   vnp;
  SelEdStructPtr sesp;
  SeqLocPtr    slp;
  Int4         line = 0;
 
  for (vnp = adpfeat; vnp != NULL; vnp = vnp->next)
  {
     sesp = (SelEdStructPtr) vnp->data.ptrvalue;
     for (; sesp != NULL; sesp = sesp->next) 
     {
        if (vnp->choice ==FEATDEF_CDS && sesp->regiontype ==OM_REGION_SEQLOC 
        && sesp->region !=NULL)
        {
           slp = (SeqLocPtr) sesp->region;;
           if (SeqLocStart(slp) > right || SeqLocStop(slp) < left) { }
           else {
              line+= feat_linenum (SeqLocStart(slp), SeqLocStop(slp), line_len, left, right);
           }
        }
     }
  }
  return line;
}
 
static Int4 get_tot_line (EditAlignDataPtr adp, Int4 line_len, Int4 left, Int4 right)
{
  Int4         diffs, line = 0;

  if (left >= right )           
      return line;
  line=count_total_line (adp->anp_list, line_len, left, right, NULL);
  diffs = right - left;
  if ( (diffs % line_len) > 0 ) diffs += line_len;
  diffs = (Int4) (diffs / line_len);
  line += (Int4) addline_perblock (adp, diffs);
  line += (Int4) CountFeatNum (adp->feat, line_len, left, right);
  line += (Int4) CountFeatNum (adp->seqfeat, line_len, left, right);
  return line;
}
 
extern void data_collect_arrange (EditAlignDataPtr adp, Boolean recollect)
{
  ValNodePtr         list = NULL, 
                     vnp = NULL;
  Int4               maxscroll;
  Int2               x;
  Boolean            goOn = FALSE;

  x = adp->seqnumber;
  if (adp->draw_scale) x++;
  if (adp->draw_bars) x++;
  maxscroll = adp->length * x / adp->visibleWidth;
  adp->vscrollunit = MAX ((FloatLo)1.00, ((FloatLo)maxscroll/(FloatLo)32767));
  adp->nlines = get_tot_line (adp, adp->visibleWidth, 0, adp->length-1);
  adp->nlines = adp->nlines - adp->pnlLine +3;
  adp->nlines = MAX ((Int4)0, adp->nlines);
  if (recollect) 
  {
     adp->bufferstart = get_bufferstart (adp);
     adp->bufferlength= MIN ((Int4) (adp->length - adp->bufferstart), 
                             (Int4) adp->minbufferlength);
     if (abs((adp->bufferstart +adp->bufferlength ) -adp->length) < 100) 
     {
        adp->bufferlength = adp->length - adp->bufferstart;
     }
     goOn = update_fromalignnode (adp);

     if (adp->blocks!=NULL) 
{
  TextAlignBufPtr tdp, 
                  tdp1;
  CharPtr         tmp;
  DenseDiagPtr    ddp;
  Int4            start, len, j;
  Int4            start1, stop1;
  SeqIdPtr        sip,
                  seqsip;
  SeqAlignPtr     salp = (SeqAlignPtr)(adp->sap_align->data);
  Int2        chklocp;

  if (adp->blocks!=NULL) {
     for (vnp=adp->linebuff; vnp!=NULL; vnp=vnp->next) 
     {
        tdp = (TextAlignBufPtr) vnp->data.ptrvalue; 
        for (tmp=tdp->buf; *tmp!='\0'; ++tmp) 
        {
           *tmp = TO_LOWER (*tmp);
        }
     }
     tdp1 = (TextAlignBufPtr) adp->linebuff->data.ptrvalue;
     ddp =(DenseDiagPtr)adp->blocks->segs;
     if (ddp != NULL) 
     {
        for (; ddp!=NULL; ddp=ddp->next) 
        {
           sip = ddp->id;
           start = *(ddp->starts); 
           len = ddp->len;
           if (start+len > adp->gr.left +adp->bufferstart) 
           {
              start = MAX(start, adp->gr.left +adp->bufferstart);
              start1=SeqCoordToAlignCoord (start, sip, salp, 0, 0);
              chklocp =chkloc (sip, start+len-1, adp->sqloc_list, &stop1);
              stop1 =SeqCoordToAlignCoord (stop1, sip, salp, 0, chklocp);
              for (vnp=adp->linebuff; vnp!=NULL; vnp=vnp->next) 
              {
                 sip = ddp->id;
                 tdp = (TextAlignBufPtr) vnp->data.ptrvalue;
                 seqsip=(SeqIdPtr)SeqIdFromAlignNode(adp->anp_list, tdp->seqEntityID, tdp->bsp_itemID, OBJ_BIOSEQ);
                 while (sip!=NULL)
                 {
                    if (SeqIdForSameBioseq(seqsip, sip) )
                    {
                       for (j=start1; j<=stop1; j++) 
                          if (tdp1->buf[j] != '-') 
                             tdp->buf[j] = TO_UPPER(tdp->buf[j]);
                       break;
                    }
                    sip=sip->next;
                 }
              }
           }
        }      
     }
  }
}

     if ( !goOn ) {
        adp->firstssp = NULL;
        return;
     }
  }
  if (goOn || adp->firstssp == NULL) {
     arrange_buffer (adp);
  }
}

/***************************************************************
***
****************************************************************/
static void print_scale (FILE *fout, EditAlignDataPtr adp, Int4 hoffset, Int2 leftmargin, Int2 scalelength)
{
  Int4   scal;
  Char   str[128];
  Int4   j;

  for (j = 0; j < leftmargin; j++)
         fprintf (fout, " ");
  for ( j = hoffset; j < hoffset + scalelength; j++) {
         scal = (Int4)(adp->gr.left + 1 + j);
         if (scal % 10 == 0) 
         {
            sprintf (str, "%d", (int)scal);
            fprintf (fout, "     ");
            fprintf (fout, "%5s", str);
         }
         if (adp->columnpcell > 0 && j > hoffset)
            if ((Int2) j % (Int2) adp->columnpcell == 0) 
               fprintf (fout, " ");
  }
  fprintf (fout, "\n");
}

static void print_bars (FILE *fout, EditAlignDataPtr adp, Int4 hoffset, Int2 leftmargin, Int2 scalelength)
{
  Int4   scal;
  Int4   j;

  for (j = 0; j < leftmargin; j++)
         fprintf (fout, " ");
  for ( j = hoffset; j < hoffset + scalelength; j++)  {
         scal = (Int4)(adp->gr.left + 1 + j);
         if (scal % 10 == 0)  {
            fprintf (fout, "|");
         }
         else if (scal % 5 == 0)  {
            fprintf (fout, ".");
         }
         else { 
            fprintf (fout, " ");
         }
         if (adp->columnpcell > 0 && j > hoffset)
            if ((Int2) j % (Int2) adp->columnpcell == 0) fprintf (fout, " ");
  }
  fprintf (fout, "\n");
}

static CharPtr gi_query1="<A HREF=\"http://www3.ncbi.nlm.nih.gov:80/htbin-post/Entrez/query?uid=";
static CharPtr gi_query2="&form=6&db=n&Dopt=g\">";

static void print_line (FILE *fout, CharPtr label, CharPtr txt, Int2 leftmargin, Uint1 columnpcell, SeqIdPtr gi_sip, Boolean  html_option, Boolean is_emptyline)
{
  Int4     strlens;
  Int4     j;
  CharPtr  txtp;
  Boolean  goOn=FALSE;

  txtp = txt;
  if (is_emptyline) {
     for (j=0; j< StringLen(txt); j++, txtp++) 
        if (*txtp != '-') { 
           goOn=TRUE; break; }
     if (!goOn) return;
  }
  if (html_option && gi_sip != NULL) 
  {
/*
        j = (Int4) EntrezFindSeqId (gi_sip); 
        j = 0;
        if (j > 0 && j < 10000000) { 
           sprintf (gip, "%ld", (long)j);
           fprintf (fout, "%s", gi_query1);
           fprintf (fout, "%s", gip);
           fprintf (fout, "%s", gi_query2);
        }
        else {
           gi_sip = NULL;
        }
*/
  }
  strlens = 0;
  if (label!=NULL)
     strlens = (Int2)StringLen(label);
  if (strlens > 0) {
     if (strlens > leftmargin)  
        label[leftmargin] = '\0';
     fprintf (fout, label);       
     if (html_option && gi_sip != NULL) 
        fprintf (fout, "</A>");
     for (j = strlens; j < leftmargin; j++)
        fprintf (fout, " ");
  }
  else {
     for (j = 0; j < leftmargin; j++)
        fprintf (fout, " ");
  }
  txtp = txt;
  strlens = StringLen(txt);
  for (j=0; j< strlens; j++, txtp++) {
         fprintf (fout, "%c", (char)(*txtp));
         if (columnpcell > 0 && j > 0)
            if ((Int4) (j+1) % (Int4) columnpcell == 0) {
               fprintf (fout, " ");
            }
  }
  fprintf (fout, "\n");
}

static Uint1 mycode (Char c)
{
  if (c == 't') return 0;
  if (c == 'g') return 1;
  if (c == 'a') return 2;
  if (c == 'c') return 3;
  return 4;
}

static Char mydecode (Uint1 c)
{
  if (c == 0) return 't'; 
  if (c == 1) return 'g';
  if (c == 2) return 'a';
  if (c == 3) return 'c';
  return '?';
}

static CharPtr MakeConsensus (ValNodePtr anp_list, ValNodePtr linebuff, Int4 width, Int2 numberalignline)
{
  TextAlignBufPtr tdp;
  AlignNodePtr    anp;
  SelStructPtr    ssp;
  CharPtr         strPtr;
  CharPtr         bufstr;
  CharPtr         cons;
  Int4            widthtmp;
  Int2            index;
  Int4            maxvalue;
  Uint2           itemID, entityID, itemtype;
  Int4            j;
  Uint1           k, kmaxvalue;

  Int4Ptr tab;

  if (GetNumberObjMgrSelect() == 0) return NULL;
  index = 1;
  widthtmp = width;
  bufstr =next_notemptyline (anp_list, linebuff, numberalignline, &index, 0, &widthtmp, &tdp, &anp);

  tab = (Int4Ptr) MemNew((size_t) ((4*widthtmp + 1) * sizeof(Int4)));
  for(j=0; j<4*widthtmp; j++) tab[j]=0;

  while ( index <= numberalignline && bufstr != NULL) 
  {
         if ( OBJ_(tdp->feattype) == OBJ_BIOSEQ ) {
              itemID = anp->bsp_itemID;
              entityID = tdp->seqEntityID;
         } else {
              itemID = tdp->itemID;
              entityID = tdp->entityID;
         }
         itemtype = OBJ_(tdp->feattype);
         ssp = is_selectedbyID (entityID, itemID, itemtype);
         if (ssp != NULL) {
              strPtr = tdp->buf;
              for(j=0; j<widthtmp; j++, strPtr++) 
              {
                   if ((k = mycode(*strPtr)) < 4)
                      tab[4 * j + k]++;
              }
         }
         index++;
         bufstr =next_notemptyline (anp_list, linebuff, numberalignline, &index, 0, &widthtmp, &tdp, &anp);
  }
  cons = (CharPtr) MemNew((size_t) ((widthtmp + 1) * sizeof(Char)));
  for (j = 0; j < widthtmp; j++) { 
         maxvalue = 0;
         for (k=0; k<4; k++) {
              if (tab[4*j + k] > maxvalue) {
                 maxvalue = tab[4*j + k];
                 kmaxvalue = k;
              }
         }
         cons[j] = (Char) mydecode(kmaxvalue);
  }
  MemFree (tab);
  return cons;
}

static CharPtr restrict_todiff (CharPtr str1, CharPtr str0)
{
  size_t   j;

  if (str1 != NULL && str0 != NULL) {
     for (j=0; j<MIN(StringLen(str1), StringLen(str0)); j++) {
        if (str1[j] == str0[j]) str1[j] = '.';
     }
  }
  return str1;
}



/************************************************
****  get_master sequence  
************************************************/
extern CharPtr get_master (ValNodePtr linebuff, Uint2 entityID, Uint2 itemID, Uint2 itemtype)
{
  ValNodePtr      vnp;
  TextAlignBufPtr tap;
  Uint2           tentityID;

  if ( linebuff == NULL ) return NULL;
  for (vnp = linebuff; vnp != NULL; vnp = vnp->next)
  {
         tap = (TextAlignBufPtr) vnp->data.ptrvalue;
         if ( tap != NULL)
         {
            if (OBJ_(tap->feattype) == OBJ_BIOSEQ) 
                 tentityID = tap->seqEntityID;
            else tentityID = tap->entityID;
            if (tentityID == entityID && tap->bsp_itemID == itemID && OBJ_(tap->feattype) == itemtype)
               break;
         }
  }
  if (vnp==NULL || tap == NULL) return NULL;
  if (tap->buf == NULL) return NULL; 
  return (tap->buf);
}


/************************************************************************
***  read_buffer_fromalignnode
*************************************************************************/
extern Boolean read_buffer_fromalignnode (EditAlignDataPtr adp, ValNodePtr *linebuff, Int4 bufferstart, Int4 minbufferlength, Int2 *numberalignline)
{
  AlignNodePtr  anp;
  ValNodePtr    curr;      /*for the list of AlignNodePtr*/
  ValNodePtr    list;      /*list of DrawText*/
  ValNodePtr    tdp_list;  /*list of DrawText*/
  ValNodePtr    linebufftmp = NULL; 
  TextAlignBufPtr tdp;
  Int4          p_stop=0;

  *linebuff = NULL;
  *numberalignline = 0;
  if ( adp->anp_list == NULL ) {
         ErrPostEx (SEV_ERROR, 0, 0, "fail in read_buffer_fromalignnode [1]");
  }
  for(curr = adp->anp_list; curr !=NULL; curr = curr->next)
  {
         anp = (AlignNodePtr) curr->data.ptrvalue;
         if ( anp == NULL ) {
                ErrPostEx (SEV_ERROR, 0, 0, "fail in read_buffer_fromalignnode [2]");
                break;
         }
                /*generate the DrawText buffer*/
         list = (ValNodePtr) ProcessTextAlignNode (anp, 
                (Int4) (adp->gr.left +bufferstart), 
                (Int4) (adp->gr.left + bufferstart + minbufferlength), 
                &p_stop, NULL, adp->visibleWidth, (Int1)0, (Uint4)0, NULL);
         if ( list == NULL ) {
                ErrPostEx (SEV_ERROR, 0, 0, "fail in read_buffer_fromalignnode [3]");
                break;
         }
	 tdp_list = list;
         while ( tdp_list != NULL )
	 {
                tdp = (TextAlignBufPtr) tdp_list->data.ptrvalue;
                if(tdp->buf != NULL)
                {
                       if (tdp->label == NULL) {
                          tdp->label=(CharPtr)MemNew((size_t)(64*sizeof(Char)));
                          emptystring (tdp->label, 64);
                          StringNCpy_0 (tdp->label, "unknown", 7);
                       }
                       (*numberalignline)++;
                       ValNodeAddPointer (&linebufftmp, 0, (Pointer) tdp);
                }
                tdp_list->data.ptrvalue = NULL;  /*<<<<<<<<<NEW */
                tdp_list = tdp_list->next;
         }
         ValNodeFree (list);                     /*<<<<<<<<<NEW */
  }
  if (linebufftmp == NULL)
     return FALSE;
  *linebuff = linebufftmp;
  return TRUE;
}

/**********************************
***
*** BioseqTrimN
***   truncates the nnn's at the extremities of a bioseq bsp
***   providing TopSeqEntry sep allows to modify the SeqAlign if any
***
***********************************/
static void SeqAlignDeleteByLocCallback (SeqEntryPtr sep, Pointer mydata,
                                          Int4 index, Int2 indent)
{
  BioseqPtr          bsp;
  BioseqSetPtr       bssp;
  SeqAlignPtr        salp,
                     salptmp;
  SeqLocPtr          slp;

  slp = (SeqLocPtr)(mydata);
  if (slp!=NULL && sep != NULL && sep->data.ptrvalue) {
     if (IS_Bioseq(sep)) {
        bsp = (BioseqPtr) sep->data.ptrvalue;
        if (bsp!=NULL) {
           salp = is_salp_in_sap (bsp->annot, 2);
           if (salp!=NULL) {
              for (salptmp=salp; salptmp!=NULL; salptmp=salptmp->next)
                 salptmp = SeqAlignDeleteByLoc  (slp, salptmp);
           }
        }
     }
     else if(IS_Bioseq_set(sep)) {
        bssp = (BioseqSetPtr)sep->data.ptrvalue;
        if (bssp!=NULL) {
           salp = is_salp_in_sap (bssp->annot, 2);
           if (salp!=NULL) {
              for (salptmp=salp; salptmp!=NULL; salptmp=salptmp->next)
                 salptmp = SeqAlignDeleteByLoc  (slp, salptmp);
           }
        }
     }
  }
}

extern Boolean BioseqTrimN (BioseqPtr bsp, SeqEntryPtr sep)
{
  SeqIdPtr      sip;
  SeqLocPtr     slp1 = NULL,
                slp2 = NULL;
  CharPtr       str;
  Int4          j,
                lens;
  Boolean       truncate = FALSE;
  
  if (bsp==NULL)
     return FALSE;
  sip = bsp->id;
  str = load_seq_data (sip, -1, -1, FALSE, &lens);
  if (str != NULL) 
  {
     j = lens-1;
     while (j>0) {
        if (str[j] != 'n' && str[j] != 'N') 
           break;
        j--;
     }
     if (j<lens-1) 
     {
        slp1 = SeqLocIntNew (j+1, lens-1, Seq_strand_plus, sip);
        SeqDeleteByLoc (slp1, TRUE, FALSE);
        truncate = TRUE;
     }
     j=0;
     while (j<lens) {
        if (str[j] != 'n' && str[j] != 'N') 
           break;
        j++;
     }
     if (j>0) {
        slp2 = SeqLocIntNew (0, j-1, Seq_strand_plus, sip);
        SeqDeleteByLoc (slp2, TRUE, FALSE);
        truncate = TRUE;
     }
     if (slp1!=NULL) {
        if (sep!=NULL)
           SeqEntryExplore (sep, (Pointer)slp1, SeqAlignDeleteByLocCallback);
        ValNodeFree (slp1);
     }
     if (slp2!=NULL) {
        if (sep!=NULL)
           SeqEntryExplore (sep, (Pointer)slp2, SeqAlignDeleteByLocCallback);
        ValNodeFree (slp2);
     }
  }
  return truncate;
}


/**********************************************************
***  GetFeatureForEditor
***
***********************************************************/
extern Boolean is_newfeat (ValNodePtr feathead, Uint2 eID, SeqLocPtr slp)
{
  SelEdStructPtr   psp;
  ValNodePtr       vnp = NULL;

  if (feathead == NULL)  { 
         return TRUE;
  }
  for (vnp = feathead; vnp != NULL; vnp = vnp->next)
  {
        psp = (SelEdStructPtr) vnp->data.ptrvalue;
        if (psp->entityID == eID)
        {
           if (SeqLocCompare(slp, (SeqLocPtr) psp->region) == SLC_A_EQ_B) {
              return FALSE;
           }
        }
  }
  return TRUE;
}

extern ValNodePtr AddFeatFunc (SelEdStructPtr feat, ValNodePtr *featlist, Uint2 itemsubtype)
{
  SelEdStructPtr   psp, 
                   prepsp, tmp;
  ValNodePtr       feathead;
  ValNodePtr       vnp = NULL;
  Int4             featstart, pspstart, pspnextstart;
  Int1             insert;
  Uint2            itemID;

  if (feat == NULL) 
         return *featlist;
  itemID = feat->itemID;

  feathead = *featlist;
  if (feathead == NULL)
  {
         feathead = ValNodeNew (NULL);
         if (feathead == NULL) {
                return feathead;
         }
         feathead->choice = (Uint1)itemsubtype;
         feathead->data.ptrvalue = (Pointer) feat;
         feat->prev = NULL;
         return feathead;
  }
  if (feathead != NULL) 
  {
         for (vnp = feathead; vnp != NULL; vnp = vnp->next)
         {
            if (vnp->choice == itemsubtype) 
            {
               psp = (SelEdStructPtr) vnp->data.ptrvalue;
               if (is_sameses (psp, feat))
               {
                 insert = 0;
                 featstart = SeqLocStart((SeqLocPtr)feat->region);
                 prepsp = NULL;
                 while (psp!= NULL) 
                 {
                   if (psp->next == NULL) {
                      if (itemID == psp->itemID) {
                         if (overlapp_ssp ((SeqLocPtr)feat->region, (SeqLocPtr) psp->region)
) { 
/*{
CharPtr str1, str2;
str1 = SeqLocPrint((SeqLocPtr) feat->region);
str2 = SeqLocPrint((SeqLocPtr) psp->region);
WriteLog("%d %d   %d %d>>> %s  >> %s\n", feat->entityID, feat->itemID, psp->entityID, psp->itemID, str1, str2);
}*/
                            insert = 99; /*break; */
                         }
                         pspstart = SeqLocStart((SeqLocPtr)psp->region);
                         if ( featstart < pspstart) { 
                            insert=-1; break; }
                         else { 
                            insert = +1; break; }
                      }
                      if (itemID < psp->itemID) { insert=-1; break; }
                      if (itemID > psp->itemID) { insert=+1; break; }
                   }
                   else {
                      pspstart = SeqLocStart((SeqLocPtr)psp->region);
                      pspnextstart = SeqLocStart((SeqLocPtr)psp->next->region);
                      if (itemID == psp->itemID ) { 
                         if (overlapp_ssp ((SeqLocPtr)feat->region, (SeqLocPtr) psp->region)) { 
/*{
CharPtr str1, str2;
str1 = SeqLocPrint((SeqLocPtr) feat->region);
str2 = SeqLocPrint((SeqLocPtr) psp->region);
WriteLog("%d %d   %d %d>>> %s  >> %s\n", feat->entityID, feat->itemID, psp->entityID, psp->itemID, str1, str2);
}*/
                            insert = 99; break; }
                         if (itemID == psp->itemID && featstart < pspstart) { 
                            insert=-1; break; 
                         }
                         if (itemID == psp->next->itemID) { 
                            if (overlapp_ssp ((SeqLocPtr)feat->region, (SeqLocPtr) psp->next->region)) { 
/*{
CharPtr str1, str2;
str1 = SeqLocPrint((SeqLocPtr) feat->region);
str2 = SeqLocPrint((SeqLocPtr) psp->region);
WriteLog("%d %d   %d %d>>> %s  >> %s\n", feat->entityID, feat->itemID, psp->entityID, psp->itemID, str1, str2);
}*/
                               insert = 99; break; }
                            if(featstart >pspstart && featstart <pspnextstart) {
                               insert=+1; break; 
                            }
                         }
                         if (itemID < psp->next->itemID) { 
                            if (featstart > pspstart ) {
                               insert=+1; break; 
                            }
                         }
                      }
                      if (itemID >psp->itemID && itemID ==psp->next->itemID) { 
                         if (overlapp_ssp((SeqLocPtr)feat->region, (SeqLocPtr) psp->next->region)) {
/*{
CharPtr str1, str2;
str1 = SeqLocPrint((SeqLocPtr) feat->region);
str2 = SeqLocPrint((SeqLocPtr) psp->region);
WriteLog("%d %d   %d %d>>> %s  >> %s\n", feat->entityID, feat->itemID, psp->entityID, psp->itemID, str1, str2);
}*/
                            insert=99; break; }
                         if (featstart < pspnextstart ) {
                            insert=+1; break; 
                         }
                      }
                      if (itemID <psp->itemID) { 
                         insert=-1; break; 
                      }
                   }
                   prepsp = psp;
                   psp = psp->next;
                 }
                 if (insert == 99) {
/* exons overlapp */
                    feat = MemFree (feat);
                    return feathead;
                 }
                 if (insert < 0 && prepsp == NULL) {
                    feat->next = psp;
                    vnp->data.ptrvalue = (Pointer) feat;
                    psp->prev = feat;
                    feat->prev = NULL;
                 }
                 else if (insert < 0 && prepsp != NULL) {
                    tmp = prepsp->next;
                    prepsp->next = feat;
                    feat->next = tmp;
                    psp->prev = feat;
                    feat->prev = prepsp;
                 }
                 else  if (insert > 0 && psp->next == NULL) {
                    psp->next = feat;
                    feat->prev = psp;
                 }
                 else  if (insert > 0 && psp->next != NULL) {
                    tmp = psp->next;
                    psp->next = feat;
                    feat->next = tmp;
                    tmp->prev = feat;
                    feat->prev = prepsp;
                 }
                 else 
                    ErrPostEx (SEV_ERROR, 0, 0, "fail in MadeFeatProc [99]"); 
                 break;
               }
            }
         }
  }
  if (vnp == NULL) {
         vnp = ValNodeAddPointer (&feathead, 0, (Pointer) feat);
         vnp->choice = (Uint1)itemsubtype;
         return feathead;
  }
  return feathead;
}

/****************************************************************
*
*   satcollfunc()
*   callback function for collecting features on Sequence 
*   alignment. It recalculates the feature intervals based on 
*   the intervals in the aligned segments
*
****************************************************************/
typedef struct alignfeat
{
   ObjMgrPtr   omp;
   ValNodePtr  data;
   Int2        filter_level;
   Uint2       entityID;
   CollectSeqOptionPtr csop;
   BioseqPtr   bsp;

}AlignFeat, PNTR AlignFeatPtr;

static Boolean slpfeatcollfunc(GatherContextPtr gcp)
{
  CollectSeqOptionPtr csop;
  ObjMgrTypePtr  omtp;
  SeqFeatPtr     sfp;
  AlignFeatPtr   afp;
  SeqLocPtr      slp = NULL;
  SeqLocPtr      curr;
  SelEdStructPtr feat;
  SeqIdPtr       sip;
  BioseqPtr      bsp;
  CdRegionPtr    crp;
  Char           label[101];
  Int4           start, stop;
  Uint2          eID, iID, bspID;
  Uint2          feat_subtype;   /*types defined by objfdef.h*/
  Int2           label_size;
  Uint1          strand;
  Uint1          codonstart;

  if(gcp->thistype != OBJ_SEQFEAT)
      return TRUE;
  afp = (AlignFeatPtr)(gcp->userdata);
  if(afp == NULL || afp->csop == NULL)
      return FALSE;
  if(afp->filter_level == gcp->seglevel+1)
      return TRUE;
  csop = afp->csop;
  if(csop->features == NULL)
      return FALSE;

  omtp = ObjMgrTypeFind (afp->omp, OBJ_SEQFEAT, NULL, NULL);
  if(omtp == NULL)
      return TRUE;

  feat_subtype = 0;
  if(omtp->subtypefunc !=NULL)
      feat_subtype = (*(omtp->subtypefunc)) (gcp->thisitem); 
   /*do not collect the current feature*/
  if(csop->features[feat_subtype] == FALSE)
      return TRUE;
  sfp = gcp->thisitem;
/*
  label_size = MIN(100, csop->label_size);
*/
  label_size = 50;
  label [0] = '\0';
  if (omtp->labelfunc != NULL)
  {
        (*(omtp->labelfunc))(sfp, label, label_size, csop->flabel_format [feat_subtype] );
  }
  if(gcp->product)
     slp = sfp->product;
  else {
     slp = sfp->location;
     eID = gcp->entityID;
     iID = gcp->itemID;
     bspID = afp->entityID;
     if (sfp->data.choice == SEQFEAT_CDREGION) {
        crp = (CdRegionPtr) sfp->data.value.ptrvalue;
        codonstart = crp->frame;
     }
     else codonstart = 0;
     if (slp->choice == SEQLOC_PACKED_INT || slp->choice == SEQLOC_MIX)
     {
        strand = SeqLocStrand (slp);
        curr = NULL;
        while ((curr = SeqLocFindNext(slp, curr)) != NULL)
        {
           bsp = afp->bsp;
           if (bsp != NULL) {
                 start = GetOffsetInBioseq (curr, bsp, SEQLOC_LEFT_END);
                 stop = GetOffsetInBioseq (curr, bsp, SEQLOC_RIGHT_END);
                 if (start != -1 || stop != -1) 
                 {
                    start = MAX (SeqLocStart (curr), 0);
                    stop = MIN ((bsp->length - 1), SeqLocStop (curr));
                    sip = SeqIdFindBest (bsp->id, 0);
                    feat = new_seledstruct (eID, iID, OBJ_SEQFEAT, bspID, start, stop,  sip, strand, FALSE, label, NULL, 0 , codonstart);
                    afp->data=AddFeatFunc(feat, &(afp->data), feat_subtype);
                 }
           }
           codonstart = 1;
        }
     }
     else  {
        bsp = afp->bsp;
        if (bsp != NULL) {
           start = GetOffsetInBioseq (slp, bsp, SEQLOC_LEFT_END);
           stop = GetOffsetInBioseq (slp, bsp, SEQLOC_RIGHT_END);
           if (start != -1 || stop != -1) 
           {
                 start = MAX (SeqLocStart (slp), 0);
                 stop = MIN ((bsp->length - 1), SeqLocStop (slp));
                 sip = SeqIdFindBest (bsp->id, 0);
                 feat = new_seledstruct (eID, iID, OBJ_SEQFEAT, bspID, start, stop,  sip, SeqLocStrand(slp), FALSE, label, NULL, 0, codonstart);
                 afp->data = AddFeatFunc (feat, &(afp->data), feat_subtype);
           }
        }
     }
  }
  return TRUE;
}

/******************************************************************
*
*    CollectFeatureForEditor (slp, anp, csop)
*   collect feature for the alignment
*   slp: the target Seq-loc
*   anp: the AlignNode belong to the target Seq-loc
*   csop: the option for gathering the features
*   
******************************************************************/
extern ValNodePtr CollectFeatureForEditor (SeqLocPtr slp, ValNodePtr seqfeat, Uint2 seq_entityID, Uint2 bsp_itemID, Uint1 *featOrder, Boolean all_feat)
{
  CollectSeqOption cs_option;
  GatherScope      gs;
  AlignFeat        af;
  BioseqPtr        bsp;
  ValNodePtr       vnp;
  ValNodePtr       vnp2, next;
  SelEdStructPtr   psp;
  Int2             j;
  Uint1Ptr         flabel_format = NULL;
  Boolean          show_feature, collect = FALSE;
   
  if(slp == NULL)
      return FALSE;
  if(seq_entityID == 0)
      return FALSE;
  cs_option.nointerval = FALSE;
  cs_option.slabel_format = PRINTID_FASTA_LONG;   /*PRINTID_TEXTID_ACCESSION;*/
  cs_option.seglevels = 0;
  for( j = 0; j < FEATDEF_ANY; ++j)   
  {
     if (all_feat)
        show_feature = TRUE;
     else
        show_feature = (Boolean)(featOrder[j] != 0);
     cs_option.features[j] = show_feature;
     if(show_feature) 
        collect = TRUE;
  }
  if(collect)
  {
     MemSet ((Pointer) (cs_option.flabel_format), OM_LABEL_CONTENT, (size_t) FEATDEF_ANY*sizeof(Uint1));
     bsp = BioseqLockById(SeqLocId(slp));
        
     MemSet((Pointer)&gs, 0, sizeof (GatherScope));
     gs.get_feats_location = TRUE;
     gs.get_feats_product =( bsp->mol == Seq_mol_aa);
     MemSet ((Pointer)(gs.ignore), (int)TRUE, (size_t) (OBJ_MAX) * sizeof (Boolean));
     gs.ignore[OBJ_SEQANNOT] = FALSE;
     gs.ignore[OBJ_SEQFEAT] = FALSE;
     gs.nointervals = FALSE;   
     gs.ignore_top = FALSE;
     gs.currlevel = 0;
     gs.offset = 0;               /*anp->extremes.left;*/
     gs.target = slp;

     af.data = NULL;
     af.csop = &cs_option;
     af.omp = ObjMgrGet();
     af.filter_level = 0;
     af.entityID = bsp_itemID;
     af.bsp = bsp;
     GatherEntity (seq_entityID, (Pointer)(&af), slpfeatcollfunc, &gs);
     BioseqUnlock (bsp);

     if (seqfeat != NULL) {
        for (vnp=seqfeat; vnp->next!=NULL; vnp=vnp->next)
           continue;
        vnp2 = af.data;
        while (vnp2!=NULL) {
           next = vnp2->next;
           vnp2->next = NULL;
           psp = (SelEdStructPtr) vnp2->data.ptrvalue;
           if (is_newfeat (seqfeat, psp->entityID, (SeqLocPtr)psp->region))
           {
              vnp->next = vnp2;
              vnp = vnp2;
           }
           else {
              SelEdStructDel (psp); 
              vnp2->data.ptrvalue = NULL;
              ValNodeFree (vnp2);
           }
           vnp2 = next;
        }
     }
     else seqfeat = af.data;
     return seqfeat;
  }
  return NULL;
}

/************************************************************
***
*** CopySeqLocFromSeqAlign
***  map source_loc to the target sequence defined by target_sip.
***  Calls map_one_location, that maps one Seq-loc
***    with gap_choice == TAKE_GAP_CHOICE
***    because one SEQLOCINT seqloc has to be stoppped at gaps
***    when mapped, and returns SEQLOCMIX
*** 
#define DEFAULT_GAP_CHOICE 0    * will split seqloc 
                                * if gaps length > MAX_GAP_LEN_BY_DEFAULT *
#define IGNORE_GAP_CHOICE 1     * extends over gaps *
#define TAKE_GAP_CHOICE 2       * split at gaps *
#define MAX_GAP_LEN_BY_DEFAULT 10       *the maximum allowd length of gap*
***
*************************************************************/
extern SeqLocPtr CopySeqLocFromSeqAlign (SeqFeatPtr source_sfp, SeqIdPtr target_id, SeqIdPtr s_id, SeqAlignPtr align, Uint1 gap_choice, Uint1 *frame)
{
   GatherRange gr;
   SeqLocPtr slp, source_slp, c_slp;
   SeqLocPtr process_slp, new_slp,
             new_location;
   Boolean check_gap;
   
   SeqIdPtr sip, source_id;
   Boolean map_to_source;
   BioseqPtr source_bsp;
   CdRegionPtr crp;
   Int4 cds_len;
   Int4 aa_start, aa_stop;
   Int4 a_start, a_stop;
   Int4 r_start, r_stop, e_start, e_stop;
   Int4 frame_offset, c_frame_offset;
   Int4 offset;
   Uint1     new_frame;
   Boolean   stop_here;
   Boolean   had_first_seg;
   Boolean   p5, p3;

   if(source_sfp == NULL || target_id == NULL || align == NULL)
      return NULL;

   if(source_sfp->location == NULL)
      return NULL;

   sip = SeqLocId(source_sfp->location);
   if(sip == NULL)
      return NULL;

   if(s_id == NULL)
      source_id = sip;
   else
      source_id = s_id;

   source_bsp = BioseqLockById(source_id);

   if(source_bsp == NULL)
   {
      return NULL;
   }

   map_to_source = FALSE;   /*segmented sequence may have locations on different segment sequences*/

   if(!BioseqMatch(source_bsp, sip))
   {
      ErrPostEx (SEV_ERROR, 0, 0, "Source Bioseq does not match the Seq-id of the Source Seq-feat");
      BioseqUnlock(source_bsp);
      return NULL;
   }
   else
      map_to_source = TRUE;

   source_slp = NULL;
   if(!map_to_source)
      source_slp = SeqLocIntNew(0, source_bsp->length-1, Seq_strand_plus, sip);
   if(source_sfp->data.choice == 3)
      crp = source_sfp->data.value.ptrvalue;
   else
      crp = NULL;

   if(gap_choice == DEFAULT_GAP_CHOICE)
   {
      if(source_sfp->data.choice == SEQFEAT_CDREGION || source_sfp->data.choice == SEQFEAT_RNA)
         check_gap = TRUE;
      else
         check_gap = FALSE;
   }
   else
      check_gap = (Boolean) (gap_choice == TAKE_GAP_CHOICE);

   slp = NULL;
   process_slp = NULL;
   cds_len = 0;
   new_frame = 0;
   a_start = -1;
   a_stop = -1;
   e_start = -1;
   e_stop = -1;
   stop_here = FALSE;
   had_first_seg = FALSE;
   while((slp = SeqLocFindNext(source_sfp->location, slp)) != NULL && !stop_here)
   {
      c_slp = NULL;
      new_slp = NULL;

      if(map_to_source == FALSE)
      {   /* map the location to the coordinates on the source */
         if(SeqLocOffset(source_slp, slp, &gr, 0))
         {
            if(gr.l_trunc == FALSE && gr.r_trunc == FALSE)
               c_slp = SeqLocIntNew(gr.left, gr.right, gr.strand, source_id);
         }
      }
      else
         c_slp = slp;

      if(c_slp != NULL)
      {

         new_slp = map_one_location(c_slp, align, target_id, gap_choice, &r_start, &r_stop);
         if(new_slp != NULL)
         {
            ValNodeLink(&process_slp, new_slp);
            if(e_start == -1)
               e_start = r_start;
            else
               e_start = MIN(e_start, r_start);
            if(e_stop == -1)
               e_stop = r_stop;
            else
               e_stop = MAX(e_stop, r_stop);

               if(crp != NULL)   /*for coding region features*/
               {
               /*calculate the frame for the first exon*/
               if(!had_first_seg)
               { 
                  if(crp->frame > 1)
                     frame_offset = (Int4)crp->frame -1L;
                  else
                     frame_offset = 0L;

                  c_frame_offset = frame_offset;
                  if(cds_len > 0)
                  {
                     c_frame_offset = (cds_len - frame_offset)%3;
                     if(c_frame_offset > 0)   /*left over*/
                        c_frame_offset = 3 - c_frame_offset;
                  }
                  if(SeqLocStrand(c_slp) == Seq_strand_minus)
                     frame_offset = (SeqLocStop(c_slp) - r_stop)%3;
                  else
                     frame_offset = (r_start - SeqLocStart(c_slp))%3;
                  if(frame_offset > 0)
                     c_frame_offset = (frame_offset + c_frame_offset)%3;
                  if(c_frame_offset > 1)
                     new_frame = (Uint1) (c_frame_offset + 1);
                  else
                     new_frame = 1;
               }

               /*calculate the position in the amino acid*/
               if(SeqLocStrand(c_slp) == Seq_strand_minus)
                  offset = SeqLocStop(c_slp) - r_stop;
               else
                  offset = r_start - SeqLocStart(c_slp);
                  
               if(crp->frame > 1)
                  frame_offset = (Int4)crp->frame -1L;
               else
                  frame_offset = 0L;
               aa_start = (cds_len + offset - frame_offset)/3;
               if(aa_start < 0)
                  aa_start = 0 ;

               if(SeqLocStrand(c_slp) == Seq_strand_minus)
                  offset = SeqLocStop(c_slp) - r_start;
               else
                  offset = r_stop- SeqLocStart(c_slp);
               aa_stop= (cds_len + offset - frame_offset)/3;
               if(aa_stop < 0 )
                  aa_stop = 0;
               if(a_start == -1)
                  a_start = aa_start;
               else
                  a_start = MIN(a_start, aa_start);

               if(a_stop == -1)
                  a_stop = aa_stop;
               else
                  a_stop = MAX(a_stop, aa_stop);

              }/*finishing processing the CDS region*/
            had_first_seg = TRUE;

         }
         else
            stop_here = TRUE;
      }
      else
         stop_here = TRUE;
      cds_len += SeqLocLen(c_slp);
      if(c_slp != NULL && !map_to_source)
         SeqLocFree(c_slp);
   }
  if(process_slp == NULL)
  {
      if(source_slp != NULL)
         SeqLocFree(source_slp);
      BioseqUnlock(source_bsp);
      return NULL;
  }
  new_location = SeqLocPackage(process_slp);
  CheckSeqLocForPartial (source_sfp->location, &p5, &p3);
  if (!p5 && e_start != -1) {
     p5 = (Boolean) (SeqLocStart(source_sfp->location) < e_start);
  }
  else p5 = TRUE;
  if (!p3 && e_stop != -1) {
     p3 = (Boolean) (SeqLocStop(source_sfp->location) > e_stop);
  }
  else p3 = TRUE;
  
  SetSeqLocPartial (new_location, p5, p3);
  BioseqUnlock(source_bsp);
  *frame = new_frame;
  return new_location;
}
 

/**********************************************************************
*** ApplyNewSeqFeat
***
**********************************************************************/
#define ADD_TITLE 1
#define ADD_RRNA  2
#define ADD_CDS   3
#define ADD_IMP   4
#define ADD_PUB   5
#define ADD_GENE  6
#define ADD_REGION 7
#define ADD_BOND  8
#define ADD_SITE  9

#define first_GBFeat  15
#define number_GBFeat 58
static CharPtr GBFeat[number_GBFeat] = {
"allele", "attenuator", "C_region", "CAAT_signal", "CDS", 
"conflict", "D-loop", "D_segment", "enhancer",  "exon",  
"GC_signal", "gene",  "intron",  "J_segment",  "LTR",  
"mat_peptide", "misc_binding",   "misc_difference",  
"misc_feature", "misc_recomb",  "misc_RNA",   "misc_signal",   
"misc_structure", "modified_base",  "mutation", "N_region", 
"old_sequence", "polyA_signal",  "polyA_site", "precursor_RNA",   
"prim_transcript", "primer_bind",  "promoter",   "protein_bind",  "RBS",  
"repeat_region", "repeat_unit",  "rep_origin",  "S_region",  "satellite",  
"sig_peptide", "source",  "stem_loop",  "STS",   "TATA_signal", 
"terminator", "transit_peptide",  "unsure",   "V_region",   "V_segment",   
"variation", "virion",   "3'clip",   "3'UTR",   "5'clip",  "5'UTR",  
"-10_signal", "-35_signal"};

typedef struct applyformdata {
  Int2           type;
  Int2           errcount;
  Boolean        noLeft;
  Boolean        noRight;
  CharPtr        geneName;
  ValNodePtr     protName;
  CharPtr        rnaName;
  CharPtr        featcomment;
  CharPtr        defline;
  Uint2          key;
  Uint2          input_entityID;
  CharPtr PNTR   alist;
} ApplyFormData, PNTR ApplyFormPtr;


static SeqLocPtr find_termination_before (CharPtr str, SeqIdPtr sip)
{
  Int4 strlens;
  Int4 j;

  if (str != NULL) {
     strlens = StringLen (str);
     for (j=0; j<strlens; j++) {
        if (str[j] == '*')
           break;
     }
     if (j<strlens) {
        return (SeqLocPtr)SeqLocIntNew(j, j, Seq_strand_plus, sip); 
     }
  }
  return NULL;
}

static SeqFeatPtr find_termination_after (SeqFeatPtr sfp, Int4 stop, BioseqPtr bsp)
{
  ByteStorePtr bs;
  SeqLocPtr    slp,
               slp_tmp;
  Boolean      include_stop = FALSE;
 
  slp_tmp = sfp->location;
  sfp->location = NULL;
  slp = SeqLocIntNew (SeqLocStart(slp_tmp), stop, SeqLocStrand (slp_tmp), SeqLocId(slp_tmp)); 
  sfp->location = slp;
  bs = ProteinFromCdRegion (sfp, include_stop);
  if (bs && bsp) {
     if (bsp->seq_data)
        bsp->seq_data = BSFree (bsp->seq_data);
     bsp->seq_data = bs;
     bsp->length = BSLen (bs);
     bs = NULL;
     ValNodeFree (slp_tmp);
     slp = sfp->location;
     stop = SeqLocStart(slp) + (3*bsp->length) + 2;
     slp_tmp=SeqLocIntNew (SeqLocStart(slp), stop, SeqLocStrand (slp), SeqLocId(slp));  
     sfp->location = slp_tmp;
     ValNodeFree (slp);
  } 
  else  {
     slp = sfp->location;
     sfp->location = slp_tmp;
     ValNodeFree (slp);
  }
  return sfp;
}

static ValNodePtr get_names_from_prot (ValNodePtr sfp_product)
{
  BioseqPtr bsp;
  SeqAnnotPtr annot;
  SeqFeatPtr sfpp;
  ProtRefPtr prp;
  ValNodePtr vnp,
           new_vnp=NULL;
  CharPtr    name;

  bsp=BioseqLockById (SeqLocId(sfp_product));
  if (bsp){
     if (bsp->annot){
        for (annot=bsp->annot;annot!=NULL;annot=annot->next){
           if (annot->type==1){
              sfpp=annot->data;
              if (sfpp->data.choice==SEQFEAT_PROT) {
                 prp=(ProtRefPtr)sfpp->data.value.ptrvalue;
                 for (vnp = prp->name;vnp!=NULL;vnp=vnp->next){
                    name = StringSave((CharPtr)vnp->data.ptrvalue); 
                    ValNodeAddPointer(&new_vnp, 0, name);
                 }
              }
           }
        } 
    }
    BioseqUnlock(bsp);
  } 
  return new_vnp;
}

static void ApplyBioFeatToSeqEntry (SeqEntryPtr sep, ApplyFormPtr afp, SeqLocPtr slp, Uint1 frame, Boolean include_stop, Boolean stoptransl, SeqFeatPtr sfp_source)

{
  ByteStorePtr  bs;
  BioseqPtr     bsp;
  BioseqPtr     nbsp;
  CdRegionPtr   crp;
  ValNodePtr    descr;
  GeneRefPtr    grp;
  MolInfoPtr    mip;
  SeqEntryPtr   nsep;
  SeqEntryPtr   old;
  CharPtr       prot;
  ProtRefPtr    prp;
  SeqEntryPtr   psep;
  ValNodePtr    sdp;
  SeqFeatPtr    sfp;
  SeqFeatPtr    sfp2;
  ValNodePtr    vnp;

  SeqIdPtr sfp_sip;
  ByteStorePtr bs2;
  Int4 strlens;
  Int4 start, stop, j;
  SeqLocPtr slp_tmp;
  Boolean   changed;

  if (sep == NULL || afp == NULL) return;
  nsep = FindNucSeqEntry (sep);
  if (nsep == NULL) 
     return;
  if (afp->type == ADD_TITLE) {
    if (! stringhasnotext(afp->defline)) {
      sdp = CreateNewDescriptor (sep, Seq_descr_title);
      if (sdp != NULL) {
        sdp->data.ptrvalue = (Pointer) StringSave (afp->defline);
      }
    }
  } else if (afp->type == ADD_GENE) {  
     sfp = CreateNewFeature (nsep, NULL, SEQFEAT_GENE, sfp_source);
  } else if (afp->type == ADD_CDS) {  
     sfp = CreateNewFeature (nsep, NULL, SEQFEAT_CDREGION, sfp_source);
  } else if (afp->type == ADD_RRNA) { 
     sfp = CreateNewFeature (nsep, NULL, SEQFEAT_RNA, sfp_source);
  } else if (afp->type == ADD_PUB) { 
     sfp = CreateNewFeature (nsep, NULL, SEQFEAT_PUB, sfp_source);
  } else if (afp->type == ADD_IMP) {
     sfp = CreateNewFeature (nsep, NULL, SEQFEAT_IMP, sfp_source);
  } else if (afp->type == ADD_REGION) {
     sfp = CreateNewFeature (nsep, NULL, SEQFEAT_REGION, sfp_source);
  } else if (afp->type == ADD_BOND) {
     sfp = CreateNewFeature (nsep, NULL, SEQFEAT_BOND, sfp_source);
  } else if (afp->type == ADD_SITE) {
     sfp = CreateNewFeature (nsep, NULL, SEQFEAT_SITE, sfp_source);
  }
  if (sfp != NULL)
  {
     sfp->location = SeqLocFree (sfp->location);
     sfp->location = (SeqLocPtr) AsnIoMemCopy ((Pointer) slp, (AsnReadFunc) SeqLocAsnRead, (AsnWriteFunc) SeqLocAsnWrite);
     SetSeqLocPartial (sfp->location, afp->noLeft, afp->noRight);
     sfp->partial = (afp->noLeft || afp->noRight);
     if (! stringhasnotext(afp->featcomment)) {
        sfp->comment = (Pointer) StringSave (afp->featcomment);
     }
     if (afp->type == ADD_CDS) 
     { 
        crp = (CdRegionPtr) (sfp->data.value.ptrvalue);
        crp->frame = frame;
        bs = ProteinFromCdRegion (sfp, include_stop);
        if (bs != NULL) 
        {
{
          if (stoptransl && include_stop && sfp->product!=NULL)
          {
            prot = BSMerge (bs, NULL);
            strlens = StringLen (prot);
            if (prot)
            {
             if (prot[(Int4)(strlens-1)] == '*')
             {
               prot [(Int4)(strlens-1)] = '\0';
               bs = BSFree (bs);
               bs = BSNew (StringLen (prot)+5);
               BSWrite (bs, (VoidPtr) prot, (Int4) StringLen (prot));
             }
             else {
               for (j=0; j<strlens; j++) {
                 if (prot[j] == '*')
                   break;
               }
               if (j<strlens) 
               {
                 prot [(Int4)(j)] = '\0';
                 bs = BSFree (bs);
                 bs = BSNew (StringLen (prot)+5);
                 BSWrite (bs, (VoidPtr) prot, (Int4) StringLen (prot));
               }
               else {
                 nbsp=(BioseqPtr)nsep->data.ptrvalue;
                 slp_tmp = sfp->location;
                 stop = (Int4)(nbsp->length-1);
                 j = stop-SeqLocStop(slp_tmp);
                 if (j>0) 
                 {
                  sfp->location = NULL;
                  sfp_sip = SeqLocId(slp_tmp);
                  slp=AsnIoMemCopy ((Pointer) slp_tmp, (AsnReadFunc) SeqLocAsnRead, (AsnWriteFunc) SeqLocAsnWrite);
                  slp = SeqLocInsert (slp, sfp_sip, SeqLocStop(slp_tmp), j, FALSE, NULL); 
                  sfp->location = slp;
                  bs2 = ProteinFromCdRegion (sfp, FALSE);
                  if (bs2)
                  {
                    bs = BSFree (bs);
                    bs = bs2;
                    ValNodeFree (slp_tmp);
                    slp = sfp->location;
                    j=(3*BSLen(bs))+3;
                    if (j<SeqLocLen(slp)) {
                       start=stop-(SeqLocLen(slp)-j)+1;
                       slp=SeqLocDelete(slp, sfp_sip, start, stop, FALSE, &changed);
                    }
                  }
                  else {
                    ValNodeFree (sfp->location);
                    sfp->location=slp_tmp;
                  }
                 }
                 bs2=NULL;
               }
             }
            }
            MemFree (prot);
          }

}
          bsp = BioseqNew ();
          if (bsp != NULL) 
          {
            bsp->repr = Seq_repr_raw;
            bsp->mol = Seq_mol_aa;
            bsp->seq_data_type = Seq_code_ncbieaa;
            bsp->seq_data = bs;
            bsp->length = BSLen (bs);
            bs = NULL;
            old = SeqEntrySetScope (sep);
            bsp->id = MakeNewProteinSeqId (sfp->location, NULL);
            SeqMgrAddToBioseqIndex (bsp);
            SeqEntrySetScope (old);
            psep = SeqEntryNew ();
            if (psep != NULL) {
              psep->choice = 1;
              psep->data.ptrvalue = (Pointer) bsp;
              SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, psep);
              mip = MolInfoNew ();
              if (mip != NULL) {
                mip->biomol = 8;
                mip->tech = 8;
                if (afp->noLeft && afp->noRight) {
                  mip->completeness = 5;
                } else if (afp->noLeft) {
                  mip->completeness = 3;
                } else if (afp->noRight) {
                  mip->completeness = 4;
                }
                vnp = CreateNewDescriptor (psep, Seq_descr_molinfo);
                if (vnp != NULL) {
                  vnp->data.ptrvalue = (Pointer) mip;
                }
              }
              descr = ExtractBioSourceAndPubs (sep);
              /*
              AddSeqEntryToSeqEntry (sep, psep, FALSE);
              */
              AddSeqEntryToSeqEntry (sep, psep, TRUE);
              ReplaceBioSourceAndPubs (sep, descr);
              SetSeqFeatProduct (sfp, bsp);

              if (afp->protName)
              {
               vnp=afp->protName;
               prot = (CharPtr)vnp->data.ptrvalue;
               if (! stringhasnotext(prot))  
               {
                prp = CreateNewProtRef (prot, NULL, NULL, NULL);
                if (prp != NULL) {
                  vnp=vnp->next;
                  while (vnp) {
                    ValNodeAddPointer (&(prp->name), 0, (CharPtr)vnp->data.ptrvalue);
                    vnp=vnp->next;
                  }
                  sfp2 = CreateNewFeature (psep, NULL, SEQFEAT_PROT, NULL);
                  if (sfp2 != NULL) {
                    sfp2->data.value.ptrvalue = (Pointer) prp;
                    SetSeqLocPartial (sfp2->location, afp->noLeft, afp->noRight);
                    sfp2->partial = (afp->noLeft || afp->noRight);
                  }
                }
               }
              }
            }
          }
        }
     }
     if (afp->type==ADD_CDS || afp->type==ADD_RRNA || afp->type==ADD_IMP) {
      if (! stringhasnotext (afp->geneName)) {
        grp = CreateNewGeneRef (afp->geneName, NULL, NULL, FALSE);
        if (grp != NULL) {
           sfp = CreateNewFeature (nsep, NULL, SEQFEAT_GENE, NULL);
           if (sfp != NULL) {
              sfp->data.value.ptrvalue = (Pointer) grp;
              sfp->location = SeqLocFree (sfp->location);
              sfp->location = AsnIoMemCopy ((Pointer) slp, (AsnReadFunc) SeqLocAsnRead, (AsnWriteFunc) SeqLocAsnWrite);
              SetSeqLocPartial (sfp->location, afp->noLeft, afp->noRight);
              sfp->partial = (afp->noLeft || afp->noRight);
           }
        }
      }
     }
  }
}

extern Boolean ApplyNewSeqFeat (ValNodePtr vnpfeat, ValNodePtr vnpsfp, Boolean stoptransl)
{
  ValNodePtr       vnps,
                   vnpf;
  ApplyFormData    af;
  SelEdStructPtr   sesp,
                   sesp1;
  SeqFeatPtr       sfp;
  SeqEntryPtr      sep;
  SeqEntryPtr      sep2;
  SeqLocPtr        new_slp;
  SeqEntryPtr      sep_head;
  Boolean          val = FALSE,
                   p5, p3;
  CdRegionPtr      crp;
  CodeBreakPtr     cbp;
  CcId             ci;

  if (vnpfeat!=NULL && vnpsfp!=NULL)
  {
     MemSet ((Pointer)(&af), 0, sizeof(ApplyFormData));
     af.errcount = 0;
     af.geneName = NULL;
     af.rnaName = NULL;
     af.protName = NULL;
     af.featcomment = NULL;
     af.defline = NULL;
     af.key = 0;
     af.alist = GBFeat;
     val = FALSE;
     vnps= vnpsfp;
     for (vnpf=vnpfeat; vnpf!=NULL; vnpf=vnpf->next, vnps=vnps->next) 
     {
        sesp = (SelEdStructPtr) vnpf->data.ptrvalue;
        new_slp = (SeqLocPtr) sesp->region;
        sep_head  = GetTopSeqEntryForEntityID (sesp->entityID);

        ci.sip = SeqIdDup (SeqIdFindBest(SeqLocId(new_slp), 0));
        ci.sep = NULL;
        ci.bsp = NULL;
        SeqEntryExplore(sep_head,(Pointer)&ci, FindSeqEntryForSeqIdCallback);
        af.input_entityID = SeqMgrGetEntityIDForSeqEntry(ci.sep);
        sep2 = GetBestTopParentForData (af.input_entityID, ci.bsp);
        SeqIdFree (ci.sip);

        if (!val) 
           sesp1 = sesp;
        if (sep2 != NULL && vnpf!=NULL && sesp!=NULL && new_slp!=NULL) 
        {
           af.input_entityID = SeqMgrGetEntityIDForSeqEntry(sep2);
           if (af.input_entityID != 0) {
              sep = GetTopSeqEntryForEntityID (af.input_entityID);
              if (sep != NULL) {
                 if (sesp->itemtype == FEATDEF_GENE) {
                    af.type = ADD_GENE;
                    af.geneName = sesp->label; 
                 }
                 else if (sesp->itemtype == FEATDEF_CDS) {
                    af.type = ADD_CDS;
                    sfp=(SeqFeatPtr)(vnps->data.ptrvalue);
/*
                    if (sfp->data.value.ptrvalue!=NULL){
                       grp = (GeneRefPtr) sfp->data.value.ptrvalue;
                       af.geneName= StringSave (grp->locus);  
                    } 
*/
                    if (sfp->data.choice == SEQFEAT_CDREGION) {
                       crp = (CdRegionPtr) sfp->data.value.ptrvalue;
                       if (crp->code_break !=NULL) {
                          for (cbp=crp->code_break; cbp!=NULL; cbp=cbp->next)
                             cbp->loc = SeqLocReplaceID (cbp->loc, SeqLocId(new_slp));
                       }
                       if (sfp->product!=NULL)
                          af.protName = get_names_from_prot (sfp->product);
                       
                    }
                    if (af.protName==NULL && sesp->label)
                       ValNodeAddPointer (&(af.protName), 0, sesp->label);
                 }
                 else if (sesp->itemtype == FEATDEF_preRNA ) {
                    af.type = ADD_RRNA;
                    af.key = RNA_TYPE_premsg;
                    af.rnaName = sesp->label;
                 }
                 else if (sesp->itemtype == FEATDEF_mRNA ) {
                    af.type = ADD_RRNA;
                    af.key = RNA_TYPE_mRNA;
                    af.rnaName = sesp->label;
                 }
                 else if (sesp->itemtype == FEATDEF_tRNA ) {
                    af.type = ADD_RRNA;
                    af.key = RNA_TYPE_tRNA;
                    af.rnaName = sesp->label;
                 }
                 else if (sesp->itemtype == FEATDEF_rRNA ) {
                    af.type = ADD_RRNA;
                    af.key = RNA_TYPE_rRNA;
                    af.rnaName = sesp->label;
                 }
                 else if (sesp->itemtype == FEATDEF_snRNA ) {
                    af.type = ADD_RRNA;
                    af.key = RNA_TYPE_snRNA;
                    af.rnaName = sesp->label;
                 }
                 else if (sesp->itemtype == FEATDEF_scRNA ) {
                    af.type = ADD_RRNA;
                    af.key = RNA_TYPE_scRNA;
                    af.rnaName = sesp->label;
                 }
                 else if (sesp->itemtype == FEATDEF_otherRNA ) {
                    af.type = ADD_RRNA;
                    af.key = RNA_TYPE_other;
                    af.rnaName = sesp->label;
                 }
                 else if (sesp->itemtype == FEATDEF_PUB ) {
                    af.type = ADD_PUB;
                 }
                 else if (sesp->itemtype == FEATDEF_REGION ) {
                    af.type = ADD_REGION;
                 }
                 else if (sesp->itemtype == FEATDEF_BOND ) {
                    af.type = ADD_BOND;
                 }
                 else if (sesp->itemtype == FEATDEF_SITE ) {
                    af.type = ADD_SITE;
                 }
                 else if (sesp->itemtype >= first_GBFeat && sesp->itemtype < number_GBFeat) { 
                    af.type = ADD_IMP;
                    af.key = sesp->itemtype;
                 }
                 else af.type = 0;
                 if (af.type > 0) {
                    CheckSeqLocForPartial (new_slp, &p5, &p3);
                    af.noLeft = p5;
                    af.noRight = p3;
                    ApplyBioFeatToSeqEntry(sep2, &af, new_slp, sesp->codonstart, TRUE, stoptransl, (SeqFeatPtr)(vnps->data.ptrvalue));
                    val = TRUE;
                 }
              } 
           } 
        } 
     }
  } 
  return val;
}

/**********************************************************************
*** PropagateFeatureBySeqLock
***   propagates features (seqfeat) from a bioseq (source_bspitemID)
***   to another (target_sep)
***
**********************************************************************/
static Boolean FindSqFeatItem (GatherContextPtr gcp)
{
  SeqFeatPtr PNTR sfpp;

  sfpp = (SeqFeatPtr PNTR) gcp->userdata;
  if (sfpp != NULL && gcp->thistype == OBJ_SEQFEAT) {
    *sfpp = (SeqFeatPtr) gcp->thisitem;
  }
  return TRUE;
}

extern void PropagateFeatureBySeqLock (SeqAnnotPtr sap, Uint2 source_bspitemID, Uint2 target_entityID, SeqEntryPtr target_sep, ValNodePtr seqfeat, Uint1 gap_choice)
{
  BioseqPtr        target_bsp;
  SeqFeatPtr       source_sfp;
  SeqAlignPtr      salp;
  SelEdStructPtr   feat;
  SeqIdPtr         featsip;
  SeqIdPtr         target_sip;
  SeqLocPtr        featslp;
  SeqLocPtr        new_slp;
  ValNodePtr       vnpf;
  ValNodePtr       vnpfeat = NULL;
  ValNodePtr       vnpsfp = NULL;
  Uint2            subtype;
  Uint1            frame;
  Uint1            gap_choice_subtype;
  Boolean          val;

  if (sap != NULL && target_sep != NULL && seqfeat != NULL) {
     if (sap != NULL)
     {   
        target_bsp = (BioseqPtr) target_sep->data.ptrvalue;
        salp = (SeqAlignPtr) sap->data;
        for (vnpf= seqfeat; vnpf!=NULL; vnpf=vnpf->next)
        {
           feat  = (SelEdStructPtr) vnpf->data.ptrvalue;
           val = (Boolean)(feat->bsp_itemID == source_bspitemID);
           if (val) {
              featslp = (SeqLocPtr) feat->region;
              featsip = SeqLocId (featslp);
              subtype = vnpf->choice;
              GatherItem (feat->entityID, feat->itemID, OBJ_SEQFEAT, (Pointer) (&source_sfp), FindSqFeatItem);
              if (source_sfp != NULL) {
                 target_sip=target_bsp->id;
                 if (subtype == FEATDEF_GENE)
                    gap_choice_subtype = IGNORE_GAP_CHOICE;
                 else if (subtype == FEATDEF_PUB)
                    gap_choice_subtype = IGNORE_GAP_CHOICE;
                 else
                    gap_choice_subtype = gap_choice;
                 new_slp = CopySeqLocFromSeqAlign (source_sfp, target_sip, featsip, salp, gap_choice_subtype, &frame);
                 if (new_slp != NULL) {
                    if (is_newfeat (seqfeat, target_entityID, new_slp))
                    {
                       new_slp = SeqLocReplaceID(new_slp, SeqLocId(source_sfp->location));
                       SeqLocFree (source_sfp->location);
                       source_sfp->location = new_slp;
                    }
                 }
              }
           }
        }
     }
  }
  return;
}

static void PropagateFeatureByApply (SeqAnnotPtr sap, Uint2 source_entityID, Uint2 source_bspitemID, Uint2 target_entityID, Uint2 target_bsp_itemID, SeqEntryPtr target_sep, ValNodePtr source_seqfeat, ValNodePtr target_seqfeat, Uint1 gap_choice, Boolean stoptransl)
{
  BioseqPtr        target_bsp;
  SeqFeatPtr       source_sfp;
  SeqAlignPtr      salp;
  SelEdStructPtr   feat;
  SeqIdPtr         featsip;
  SeqIdPtr         target_sip;
  SeqLocPtr        featslp;
  SeqLocPtr        new_slp;
  ValNodePtr       vnpf;
  ValNodePtr       vnpfeat = NULL;
  ValNodePtr       vnpsfp = NULL;
  Uint2            subtype;
  Uint1            frame;
  Uint1            gap_choice_subtype;
  Boolean          val;
  SeqFeatPtr       source_dup;
  SelEdStructPtr   sesp;

  if (sap != NULL && target_sep != NULL && source_seqfeat != NULL) {
     if (sap != NULL)
     {   
        target_bsp = (BioseqPtr) target_sep->data.ptrvalue;
        salp = (SeqAlignPtr) sap->data;
        for (vnpf= source_seqfeat; vnpf!=NULL; vnpf=vnpf->next)
        {
           feat  = (SelEdStructPtr) vnpf->data.ptrvalue;
           val = (Boolean)(feat->bsp_itemID == source_bspitemID);
           if (val) {
              featslp = (SeqLocPtr) feat->region;
              featsip = SeqLocId (featslp);
              subtype = vnpf->choice;
              GatherItem (feat->entityID, feat->itemID, OBJ_SEQFEAT, (Pointer) (&source_sfp), FindSqFeatItem);
              if (source_sfp != NULL) {
                 target_sip=target_bsp->id;
                 if (subtype == FEATDEF_GENE)
                    gap_choice_subtype = IGNORE_GAP_CHOICE;
                 else if (subtype == FEATDEF_PUB)
                    gap_choice_subtype = IGNORE_GAP_CHOICE;
                 else 
                    gap_choice_subtype = gap_choice;
                 new_slp = CopySeqLocFromSeqAlign (source_sfp, target_sip, featsip, salp, gap_choice_subtype, &frame);
                 if (new_slp != NULL) 
                 {
                    new_slp = SeqLocReplaceID (new_slp, SeqIdFindBest(target_sip, 0));
                    if (is_newfeat (target_seqfeat, target_entityID, new_slp))
                    {
                       source_dup = (SeqFeatPtr) AsnIoMemCopy((Pointer) source_sfp, (AsnReadFunc)SeqFeatAsnRead, (AsnWriteFunc)SeqFeatAsnWrite);
                       sesp = new_seledstruct_fromseqloc (target_entityID, target_bsp_itemID, subtype, target_bsp_itemID, new_slp, feat->label, NULL, 0, frame);
                       if (sesp != NULL) {
                          ValNodeAddPointer(&vnpfeat, 0, (Pointer) sesp);
                          ValNodeAddPointer(&vnpsfp,0, (Pointer)source_dup);
                       }
                    }
                 }
              }
           }
        }
        val = ApplyNewSeqFeat (vnpfeat, vnpsfp, stoptransl);
     }
  }
  return;
}


/**********************************************************************
*** ReplaceBioseq
***   replaces the target bioseq (target_id) by the source
***   bioseq (source_id)
***   The location of the features of the source bioseq are 
***   modified using SeqAlign (salp) as template (PropagateFeatureBySeqLock).
***
**********************************************************************/

extern void ReplaceBioseq (SeqIdPtr target_id, SeqIdPtr source_id, SeqAlignPtr salp, Uint1 gap_choice, Boolean stoptransl)
{
  SeqAnnotPtr sap;
  ValNodePtr  source_allseqfeat = NULL,
              target_allseqfeat = NULL;
  SeqLocPtr   source_slp,
              target_slp;
  BioseqPtr   source_bsp = NULL;
  BioseqPtr   target_bsp = NULL;
  SeqEntryPtr source_sep = NULL;
  SeqEntryPtr target_sep = NULL;
  Uint2       source_entityID;   
  Uint2       source_bsp_itemID;
  Uint2       target_entityID;
  Uint2       target_bsp_itemID;
  Boolean     goOn = TRUE;
 
  if (target_id==NULL || source_id == NULL)
     return;
  source_bsp = BioseqLockById (source_id);
  target_bsp = BioseqLockById (target_id);
  if (source_bsp!=NULL && target_bsp!=NULL) 
  {
     if (goOn)
     {  
           sap=SeqAnnotForSeqAlign(salp);
           source_sep = SeqMgrGetSeqEntryForData ((Pointer)source_bsp); 
           source_entityID = SeqMgrGetEntityIDForSeqEntry (source_sep);
           source_bsp_itemID = GetItemIDGivenPointer (source_entityID, OBJ_BIOSEQ, source_bsp);
           target_sep = SeqMgrGetSeqEntryForData ((Pointer)target_bsp);
           target_entityID = SeqMgrGetEntityIDForSeqEntry (target_sep);
           target_bsp_itemID = GetItemIDGivenPointer (target_entityID, OBJ_BIOSEQ, target_bsp);

           target_slp = SeqLocIntNew (0, target_bsp->length-1, Seq_strand_plus, target_id);
           target_allseqfeat = CollectFeatureForEditor (target_slp, target_allseqfeat, target_entityID, target_bsp_itemID, NULL, TRUE);
           PropagateFeatureBySeqLock (sap, target_bsp_itemID, source_entityID, source_sep, target_allseqfeat, gap_choice);

           source_slp = SeqLocIntNew (0, source_bsp->length-1, Seq_strand_plus,
source_id);
           source_allseqfeat = CollectFeatureForEditor (source_slp, source_allseqfeat, source_entityID, source_bsp_itemID, NULL, TRUE);
           PropagateFeatureByApply (sap, source_entityID, source_bsp_itemID, target_entityID, target_bsp_itemID, target_sep, source_allseqfeat, target_allseqfeat, gap_choice, stoptransl);

           target_bsp->seq_data = BSFree (target_bsp->seq_data);
           target_bsp->seq_data = BSDup(source_bsp->seq_data);
           target_bsp->seq_data_type = source_bsp->seq_data_type;
           target_bsp->length = source_bsp->length;
           BioseqRawPack (target_bsp);

           sap->data = NULL;
           SeqAnnotFree (sap);
     }
     BioseqUnlock (source_bsp);
     BioseqUnlock (target_bsp);
  }
  return;
}

static Int4 map_position_seqalign (SeqAlignPtr salp, Int4 pos, ValNodePtr sqlocs)
{
  SeqIdPtr sip1, sip2;
  Int4     tmp_pos;
  
  sip1 = SeqAlignId (salp, 0);
  sip2 = SeqAlignId (salp, 1);
  tmp_pos = SeqCoordToAlignCoord (pos, sip2, salp, 0, 0);
  tmp_pos = (Int4) AlignCoordToSeqCoord2 (tmp_pos, sip1, salp, sqlocs, 0);
  return tmp_pos;
}

extern Boolean MergeFunc (SeqIdPtr target_id, SeqIdPtr source_id, SeqAlignPtr salp, Int4 fromseq2, Int4 toseq2, ValNodePtr sqlocs, Boolean spliteditmode)
{
  BioseqPtr     bsp,
                bsp_target;
  SeqLocPtr     target_slp;
  Int4          fromseq1, toseq1,
                from1seq1, to1seq1,
                from_overlapp, to_overlapp;
  Int4          caret = -1, 
                overlapp=0;
  Boolean       ok = TRUE;

  if (target_id==NULL || source_id == NULL || toseq2 < fromseq2)
     return FALSE;
  
  bsp = BioseqCopy (NULL, source_id, fromseq2, toseq2, Seq_strand_plus, TRUE);
  if (bsp!=NULL) 
  {
     fromseq1= map_position_seqalign(salp,fromseq2, sqlocs);
     toseq1  = map_position_seqalign(salp,toseq2, sqlocs);
     from1seq1= map_position_seqalign(salp,fromseq2-1, sqlocs);
     to1seq1  = map_position_seqalign(salp,toseq2+1, sqlocs);
     overlapp = 0;
     if (toseq1>-1) {
        if (fromseq1>-1) {
           overlapp=toseq1-fromseq1 + 1;
           from_overlapp = fromseq1; 
           to_overlapp = toseq1;
        } else {
           overlapp=toseq1 + 1;
           from_overlapp = 0; 
           to_overlapp = toseq1;
        }
     }
     else if (to1seq1 > -1) {
        if (fromseq1>-1) {
           overlapp = to1seq1 - fromseq1;
           from_overlapp = fromseq1; 
           to_overlapp = to1seq1-1;
        } else if (from1seq1>-1) {
           overlapp = to1seq1 - from1seq1 -1;
           from_overlapp = from1seq1+1; 
           to_overlapp = to1seq1-1;
        }
     }
     else if (toseq1==-2 && to1seq1 ==-2)
     {
        if (fromseq1 > -1 || from1seq1 > -1) {
           if (fromseq1 < 0)
              fromseq1 = from1seq1+1;
           bsp_target = BioseqLockById (target_id);
           if (bsp_target) {
              if (bsp_target->length > fromseq1+1) {
                 overlapp = bsp_target->length - fromseq1;
                 from_overlapp = fromseq1;
                 to_overlapp = bsp_target->length-1;
              }
              BioseqUnlock (bsp_target);
           }
           else ok = FALSE;
        }
        else ok = FALSE;
     }
     if (to1seq1 == -2) 
        caret = -2;
     else if (to1seq1 == 0)
        caret = 0;
     else if (overlapp <= to1seq1) {
        if (overlapp == 0)
           caret = to1seq1;
        else 
           caret = to1seq1 - overlapp;
     }
     if (ok) {
        if (overlapp > 0) 
        { 
           target_slp=SeqLocIntNew (from_overlapp, to_overlapp, Seq_strand_plus, target_id);
           SeqDeleteByLoc (target_slp, TRUE, spliteditmode);
           SeqLocFree (target_slp);
        }
        ok = BioseqInsert (bsp->id, FIRST_RESIDUE, LAST_RESIDUE, Seq_strand_plus, target_id, caret, TRUE, TRUE, spliteditmode);
     }
  }
  return ok;
}

extern void CopyFeatFunc (SeqIdPtr target_id, SeqIdPtr source_id, SeqAlignPtr salp, Uint1 gap_choice, Boolean stoptransl)
{
  SeqAnnotPtr sap;
  ValNodePtr  source_allseqfeat = NULL,
              target_allseqfeat = NULL;
  SeqLocPtr   source_slp,
              target_slp;
  BioseqPtr   source_bsp = NULL;
  BioseqPtr   target_bsp = NULL;
  SeqEntryPtr source_sep = NULL;
  SeqEntryPtr target_sep = NULL;
  Uint2       source_entityID;
  Uint2       source_bsp_itemID;
  Uint2       target_entityID;
  Uint2       target_bsp_itemID;
  Boolean     goOn = TRUE;

  if (target_id==NULL || source_id == NULL)
     return;
  source_bsp = BioseqLockById (source_id);
  target_bsp = BioseqLockById (target_id);
  if (source_bsp!=NULL && target_bsp!=NULL)
  {
     if (goOn)
     {
           sap=SeqAnnotForSeqAlign(salp);
           source_sep = SeqMgrGetSeqEntryForData ((Pointer)source_bsp);
           source_entityID = SeqMgrGetEntityIDForSeqEntry (source_sep);
           source_bsp_itemID = GetItemIDGivenPointer (source_entityID, OBJ_BIOSEQ, source_bsp);
           target_sep = SeqMgrGetSeqEntryForData ((Pointer)target_bsp);
           target_entityID = SeqMgrGetEntityIDForSeqEntry (target_sep);
           target_bsp_itemID = GetItemIDGivenPointer (target_entityID, OBJ_BIOSEQ, target_bsp);

           source_slp = SeqLocIntNew (0, source_bsp->length-1, Seq_strand_plus, source_id);
           source_allseqfeat = CollectFeatureForEditor (source_slp, source_allseqfeat, source_entityID, source_bsp_itemID, NULL, TRUE);
           target_slp = SeqLocIntNew (0, target_bsp->length-1, Seq_strand_plus, target_id);
           target_allseqfeat = CollectFeatureForEditor (target_slp, target_allseqfeat, target_entityID, target_bsp_itemID, NULL, TRUE);
           PropagateFeatureByApply (sap, source_entityID, source_bsp_itemID, target_entityID, target_bsp_itemID, target_sep, source_allseqfeat, target_allseqfeat, gap_choice, stoptransl);

           sap->data = NULL;
           SeqAnnotFree (sap);
     }
     BioseqUnlock (source_bsp); 
     BioseqUnlock (target_bsp); 
  }
  return;
}

static void print_firstlinePHYLIP (FILE *fout, Int2 seqnumber, Int4 length)
{
  fprintf (fout, "%7d%5d\n", seqnumber, length);
}

extern void ShowAlignmentText (FILE *fout, EditAlignDataPtr adp, SelStructPtr ssp, Int2 leftmargin, Int4 printfrom, Int4 printto, Boolean html_option)
{
  ValNodePtr      linebuff=NULL;
  TextAlignBufPtr tdp;
  AlignNodePtr    anp;
  SeqIdPtr        bspsip;
  SeqAlignPtr     salp = (SeqAlignPtr) adp->sap_align->data;
  CharPtr         bufstr,
                  trans;
  CharPtr         masterbuf =NULL;
  Int4            width;
  Int4            widthtmp;
  Int4            from; 
  Uint2           itemID, entityID, itemtype;
  Int2            numberalignline = 0;
  Int2            index = 0;
  Int2            j, k;
  Boolean         firstline = TRUE;

  if (printfrom > adp->length || printto > adp->length) {
     Message(MSG_OK, "fail in ShowAlignment: %ld > %ld or %ld > %ld", (long)printfrom, (long)adp->length, (long)printto, (long)adp->length);
     return;
  }  
  if (adp->align_format == SALSA_PHYLIP)
     print_firstlinePHYLIP (fout, adp->seqnumber, adp->length);
  from = printfrom;
  width = (Int4)MIN ((Int4) adp->visibleWidth, (Int4)(printto - from +1));

  while ( read_buffer_fromalignnode (adp, &linebuff, from, width-1, &numberalignline) )
  {
     if (adp->charmode)
        masterbuf = get_master (linebuff, adp->master.entityID, adp->master.itemID, adp->master.itemtype);
     else masterbuf = NULL;

     if (adp->draw_scale && (adp->align_format != SALSA_PHYLIP)) {
         print_scale (fout, adp, from, leftmargin, (Int2)width);
     }
     if (adp->draw_bars && (adp->align_format != SALSA_PHYLIP)) {
         print_bars (fout, adp, from, leftmargin, (Int2)width);
     }
     index = 1;
     widthtmp = (Int4)width;
     bufstr = next_notemptyline (adp->anp_list, linebuff, numberalignline, 
                                &index, 0, &widthtmp, &tdp, &anp);
     while ( index <= numberalignline && bufstr != NULL) 
     {
         if ( OBJ_(tdp->feattype) == OBJ_BIOSEQ ) {
              itemID = anp->bsp_itemID;
              if (adp->input_format == OBJ_BIOSEQ) {
                 entityID = tdp->seqEntityID;
              }
              else {
                 entityID = anp->seq_entityID;
              }
              bspsip = anp->sip;
         } 
         else {
              itemID = tdp->itemID;
              entityID = tdp->entityID;
         }
         itemtype = OBJ_(tdp->feattype);
         if (not_empty_string (bufstr, width)) 
         {
            if (masterbuf != NULL && (adp->master.entityID != entityID || adp->master.itemID != itemID))
               bufstr = restrict_todiff (bufstr, masterbuf);
            if ((adp->align_format == SALSA_PHYLIP) && !firstline)
               print_line (fout, NULL, bufstr, leftmargin, adp->columnpcell, bspsip, html_option, FALSE);
            else
               print_line (fout, tdp->label, bufstr, leftmargin, adp->columnpcell, bspsip, html_option, FALSE);
         }
         if (OBJ_(tdp->feattype) == OBJ_BIOSEQ && not_empty_string (bufstr, width)) 
         {
            if (has_complement (adp->params, entityID, itemID)) {
               print_line (fout, "complement", complement_string(bufstr), leftmargin, adp->columnpcell, NULL, html_option, FALSE);
            }
            for (j=0; j<3;j++) {
               if (rf_on (adp->params, entityID, itemID, j))  {
                   trans = (CharPtr) SeqAlignTranslate (salp, entityID, from, from + widthtmp - 1,  (Uint1)(j), bspsip, adp->seqnumber, Seq_strand_plus, adp->sqloc_list);
                   if ( adp->prot_mode == PUTPROT )
                      trans = prot_to_putprot (trans);
                   print_line (fout, "RF >", trans, leftmargin, adp->columnpcell, NULL, html_option, FALSE);
               }
            }
            for (j=3, k=2; j<6;j++, k--) {
               if (rf_on (adp->params, entityID, itemID, j))  {
                   trans = (CharPtr) SeqAlignTranslate (salp, entityID, from, from + widthtmp - 1,  (Uint1)k, bspsip, adp->seqnumber, Seq_strand_minus, adp->sqloc_list);
                   if ( adp->prot_mode == PUTPROT )
                      trans = prot_to_rputprot (trans);
                   print_line (fout, "RF <", trans, leftmargin, adp->columnpcell, NULL, html_option, FALSE);
               }
            }
            if ( adp->feat != NULL ) 
            {
/*
                vnpfeat = adp->feat;
                for (; vnpfeat != NULL; vnpfeat = vnpfeat->next)
                {
                   drawfeat = is_feature_to_buffer (vnpfeat, entityID, itemID, from, widthtmp, salp, adp->input_format, adp->sqloc_list);
                   if (drawfeat != NULL) 
                   {
                      if (vnpfeat->choice == SEQFEAT_CDREGION && drawfeat->data != NULL) 
                      {
                         slp = (SeqLocPtr) drawfeat->region;
                         vnptmp = (ValNodePtr) drawfeat->data;
                         str = (CharPtr) vnptmp->data.ptrvalue;
                         if (from <= SeqLocStop(slp)) 
                         {
                          left = (Int4) MAX ((Int4) from, (Int4) SeqLocStart(slp));
                          right = (Int4) MIN ((Int4) from + widthtmp, (Int4) SeqLocStop(slp)+1);
                          if (right > left) 
                          {
                            for (j=0; j< leftmargin; j++) fprintf(fout," ");
                            for (j=0; j< left-from; j++) 
                            {
                              fprintf(fout," ");
                              if (adp->columnpcell > 0 && j > 0)
                                 if ((Int2) j % (Int2) adp->columnpcell == 0) 
                                    fprintf (fout, " ");
                            }
                            if (SeqLocStart(slp) < from) 
                              jj = (Int4)(from - SeqLocStart(slp));  
                            else jj=0;
                            for (k=0; k< right - left -1; j++, jj++, k++) 
                            {
                              fprintf (fout, "%c", (char)str[jj]);
                              if (adp->columnpcell > 0 && j > 0)
                                 if ((Int2) j % (Int2) adp->columnpcell == 0) 
                                    fprintf (fout, " ");
                            }
                            fprintf (fout, "\n");
                          }
                         }
                      }
                   }
                }
*/
            }
         }
         index++;
         widthtmp = width;
         bufstr = next_notemptyline (adp->anp_list, linebuff, numberalignline, 
                                     &index, 0, &widthtmp, &tdp, &anp);
     }
     firstline = FALSE;
     from += width;
     if (from >= printto) break;
     width = (Int4)MIN ((Int4) adp->visibleWidth, (Int4)(printto - from +1));
     fprintf (fout, "\n");
  }
  return;
}

extern void showfastagap_fromalign (SeqAlignPtr salp, Int4 line, FILE *f)
{
  BioseqPtr bsp;
  CharPtr   str,
            strp;
  Char      buffer [255];
  CharPtr   bufp;
  SeqIdPtr  sip;
  Int4      len,
            al_len;
  Int4      from, to, offset,
            ncar = 0;
  Int2      j;
  Boolean   goOn = TRUE,
            is_prot;

  Char      strLog[120];

  if (salp == NULL)
     return;
  for (j=0; j<salp->dim; j++) {
     sip = SeqAlignId (salp, j);
     if (sip!=NULL) {
      bsp = BioseqLockById (sip);
      if (bsp!=NULL) {
        is_prot = (Boolean) (ISA_aa (bsp->mol));
        SeqIdWrite (bsp->id, strLog, PRINTID_FASTA_LONG, 120);
        BioseqUnlock (bsp);
        fprintf(f, "> %s\n", strLog);
        str = load_seq_data (sip, -1, -1, is_prot, &len);
        strp = str;
        al_len = SeqAlignLength (salp);
        from = 0;
        offset = (Int4) MIN (line, (al_len-from));
        to = offset -1;
        goOn = (Boolean) (to > from);
        while (goOn) {
           fprintf(f, "%5d  ", from);
           bufp = ReadBufferFromSap (strp, buffer, salp, sip, from, to, &ncar);
           fprintf(f, "%s\n", bufp);
           strp+= ncar;
           from = to+1;
           offset = (Int4) MIN (line, (al_len-from));
           to = from +offset -1;
           goOn = (Boolean) (to > from);
        }
      }
     }
  }
}


extern SeqAnnotPtr multseqalign_from_pairseqalign (SeqAlignPtr salp)
{
  SeqAlignPtr      salptmp,
                   tmp;
  DenseSegPtr      dsp = NULL;
  SeqPortPtr       spp;
  SeqLocPtr        slp;
  SeqIdPtr         sip, 
                   siphead, siptmp, sipcur;
  BioseqPtr        bsp;
  CharPtr          bufferin[500];
  CharPtr          bufferout[500];
  CharPtr          buffertmp[2];
  Int4Ptr          starttmp;
  Int4Ptr          lenp;
  Uint1Ptr         strandp;

  Int4             algfrom, algto, alglens;
  Int4             seqfrom, seqto;
  Int4             cur;
  Int4             dsp_len, bsp_len, sum_lens;
  Int4             max_lens = 0;
  Int4             seqoffset;
  Int4             gapoffset;
  Int4             gaptotal;
  Int4             step = SALSALENLIMIT;
  Int2             curnseq;
  Int2             j;
  Int4             k, k1, k2, k3;
  Char             c1, c2;

  ValNodePtr       vnp;
  ValNodePtr       vnpfrom;
  ValNodePtr       vnpstrand;
  SeqAnnotPtr      sap;
 
  Int4             letter;
  Uint1            strandtmp = Seq_strand_unknown;
  Boolean          strand_nonull;
  Boolean          ok;
  CharPtr          str;

  if (salp==NULL)
     return NULL;

  for (salptmp=salp; salptmp!=NULL; salptmp=salptmp->next) {
     if (salptmp->type > 0 && salptmp->segtype==2) {
        dsp = (DenseSegPtr) salptmp->segs;
        break;
     }
  }  
  if (dsp==NULL)
     return NULL;
  bsp = BioseqLockById (dsp->ids);
  if (bsp == NULL)
     return NULL;
  bsp_len = bsp->length;
  sip = SeqIdFindBest(bsp->id, 0);
  BioseqUnlock (bsp);

  /*----- check if not all sequences start with gaps ----*/
  for (salptmp=salp; salptmp!=NULL; salptmp=salptmp->next) {
     if (salptmp->type > 0 && salptmp->segtype==2) {
	dsp = (DenseSegPtr) salptmp->segs;
        for (j=0, starttmp=dsp->starts; j<dsp->dim; j++, starttmp++)
           if ( *starttmp > -1) 
              break;
        ok = (Boolean) (j < dsp->dim);
        if (!ok)
           break;
     }
  }
  if (!ok)
     return NULL;
  /*---- find longest seqalign --------*/
  sum_lens = 0;
  for (salptmp=salp; salptmp!=NULL; salptmp=salptmp->next) {
     if (salptmp->type > 0 && salptmp->segtype==2) {
        dsp_len = SeqAlignStart(salptmp, 0) + SeqAlignLength (salptmp);
        if (dsp_len > sum_lens)
           sum_lens = dsp_len; 
     }
  }  
  if (sum_lens > SALSALENLIMIT) {
     ErrPostEx (SEV_ERROR, 0, 0, "Too long alignment.\n Wait for next Sequin version");
     return NULL;
  }
  step = 2*sum_lens;
  step = MIN ((Int4)step, (Int4)SALSALENLIMIT);
  for (j=0; j<2; j++) 
  {
     str = MemNew ((size_t) ((step+10) * sizeof(CharPtr)));
     buffertmp[j] = str;
     for (k=0; k<step; k++) buffertmp[j][k] = '-';
     str = MemNew ((size_t) ((step+10) * sizeof(CharPtr)));
     bufferout[j] = str;
     for (k=0; k<step; k++) bufferout[j][k] = '-';
     str = MemNew ((size_t) ((step+10) * sizeof(CharPtr)));
     bufferin[j] = str;
     for (k=0; k<step; k++) bufferin[j][k] = '-';
  }  
  siphead = SeqIdDup (sip);
  sipcur = siphead;
  vnpstrand = NULL;
  vnpfrom = NULL;
  ValNodeAddInt (&vnpfrom, 1, (Int4) 0);

  for (salptmp=salp; salptmp!=NULL; salptmp=salptmp->next) 
     if (salptmp->type > 0 && salptmp->segtype==2) 
        break;

  /**************** if 1rst sequence start with gaps in one alignment ***/
  gapoffset = gaptotal = 0;
  for (tmp=salptmp, j=0; tmp!=NULL; tmp=tmp->next, j++) {
     if (tmp->type > 0 && tmp->segtype==2) {
        dsp = (DenseSegPtr) tmp->segs;
        if (*(dsp->starts) < 0) { 
           if (j==0)
              gapoffset = *(dsp->lens);
           if (*(dsp->lens) > gaptotal) {
              gaptotal = *(dsp->lens);
           }
        }
     }
  }
  if (gaptotal > 0) {
     if (gaptotal <= gapoffset)
        gapoffset = gaptotal - gapoffset; 
  }
  algfrom = 0;
  algto = 1;
  for (cur = algfrom; cur < algto; cur += 2*step)
  {
     curnseq = 2;
     dsp = (DenseSegPtr) salptmp->segs;
     seqoffset = *(dsp->starts);
     strand_nonull = (Boolean) (dsp->strands != NULL);
     if (strand_nonull)
        strandp = dsp->strands;
     else
        strandtmp = Seq_strand_unknown;
     for (sip = dsp->ids, j=0; sip!=NULL; sip=sip->next, j++)
     {
        starttmp = dsp->starts;
        if (j==0) {
           seqfrom = 0;
        } else {
           starttmp += j;
           seqfrom = *starttmp;
        }
        lenp = dsp->lens;
        dsp_len = 0;
        sum_lens = 0;
        for(k=0; k < dsp->numseg; k++, lenp++, starttmp += dsp->dim) {
           if ( *starttmp >= 0 ) dsp_len += *lenp;
           sum_lens += *lenp;
        }
        bsp = BioseqLockById (sip);
        if (bsp == NULL) {
           return NULL;
        }
        seqto = MIN ((Int4) (seqfrom + step), (Int4) (bsp->length -1));
        if ( j != 0 ) {
           siptmp = SeqIdDup(SeqIdFindBest(bsp->id, 0));
           sipcur->next = siptmp;
           sipcur = siptmp;
           starttmp = dsp->starts;
           starttmp += j;
           lenp = dsp->lens;
           for(k=0; k < dsp->numseg; k++, lenp++, starttmp += dsp->dim) 
              if ( *starttmp >= 0 ) break;
           if (strand_nonull && *strandp==Seq_strand_minus)
              ValNodeAddInt (&vnpfrom, 1, (Int4) (*starttmp+*lenp));
           else 
              ValNodeAddInt (&vnpfrom, 1, (Int4) *starttmp);
        }
        BioseqUnlock (bsp);
        if (strand_nonull)
           strandtmp = *strandp;
        slp = SeqLocIntNew (0, seqto, strandtmp, sip);
        if ( slp == NULL) {
           return NULL;
        }
        spp = SeqPortNewByLoc (slp, Seq_code_iupacna);
/**
        SeqIdWrite (sip, strLog, PRINTID_FASTA_LONG, 120);
        if (spp!=NULL) WriteLog("mergesalp2 %s\n", strLog);
        else WriteLog("PRINTFspp NULL\n");
        WriteLog("1>  %ld  %ld    %ld %ld  > %ld   \n", dsp_len, sum_lens, seqto, seqoffset, seqfrom);
**/
        if (seqoffset < 0) seqoffset = 0;  /*!!!!!!!!!!!!!!!!!!!!*/

        if ( j == 0 && seqoffset > seqfrom) {
           seqfrom = ReadBufferFromSep (spp, bufferin[j], (Int4)seqfrom,
                                         (Int4)(seqfrom + seqoffset), 0);
        }
/**
WriteLog ("3> %d  %d   %ld  %ld \n", j, seqfrom, sum_lens, seqfrom + sum_lens);
**/
        alglens = readbuff_fromseqalign (spp, salptmp, j, bufferin[j], seqfrom, seqfrom + sum_lens, seqoffset+gapoffset, strandtmp);
        if (alglens == 0) {
           return NULL;
        }
/**
        WriteLog("4>  %d  %ld   %ld %ld   %ld < %ld   %d\n", j, alglens, sum_lens, seqfrom+sum_lens, seqfrom + dsp_len, bsp_len, strandtmp);
**/
        if ( j == 0 && (seqfrom + dsp_len) < bsp_len) {
           alglens = ReadBufferFromSep (spp, bufferin[j], 
                           (Int4)(seqfrom + dsp_len), (Int4)bsp_len, alglens);
           bufferin[j][alglens] = '-';
        }
        else { 
           for (k = alglens; k < step; k++) 
              if (bufferin[j][k] != '-') bufferin[j][k] = '-';
        }
        SeqPortFree(spp);
        ValNodeAddInt (&vnpstrand, 1, (Int4)(strandtmp));
        if (strand_nonull)
           strandp++;
     }
     for (salptmp = salptmp->next; salptmp != NULL; salptmp = salptmp->next)
     {
      if (salptmp->type > 0 && salptmp->segtype==2)
      {
/**
WriteLog ("\n\nCURRENT ALIGN %d\n", curnseq);
**/
       dsp = (DenseSegPtr) salptmp->segs;
       seqoffset = *(dsp->starts);
       if (seqoffset < 0) {
          gapoffset = gaptotal - *(dsp->lens);
       }
       else 
          gapoffset = 0;
       strand_nonull = (Boolean) (dsp->strands != NULL);
       if (strand_nonull)
           strandp = dsp->strands;
       else
           strandtmp = Seq_strand_unknown;
       for (sip = dsp->ids, j=0; sip!=NULL; sip=sip->next, j++)
       {
           dsp = (DenseSegPtr) salptmp->segs;
           for (k=0; k<step; k++) 
              buffertmp[j][k] = '-';
           starttmp = dsp->starts;
           starttmp += j;
           if (j == 0) seqfrom = 0;
           else seqfrom = *starttmp;
           lenp = dsp->lens;
           dsp_len = 0;
           sum_lens = 0;
           for(k=0; k < dsp->numseg; k++, lenp++, starttmp += dsp->dim) {
              if ( *starttmp >= 0 ) 
                 dsp_len += *lenp;
              sum_lens += *lenp;
           }
           bsp = BioseqLockById (sip);
           if (bsp == NULL){
              return NULL;
           }
           seqto = MIN ((Int4) (seqfrom + step), (Int4) (bsp->length -1));
           if ( j != 0 ) {
              siptmp = SeqIdDup(SeqIdFindBest(bsp->id, 0));
              sipcur->next = siptmp;
              sipcur = siptmp;
              starttmp = dsp->starts;
              starttmp += j;
              lenp = dsp->lens;
              for(k=0; k < dsp->numseg; k++, lenp++, starttmp += dsp->dim) 
                 if ( *starttmp >= 0 ) break;
              if (strand_nonull && *strandp==Seq_strand_minus)
                 ValNodeAddInt (&vnpfrom, 1, (Int4) (*starttmp+*lenp));
              else 
                 ValNodeAddInt (&vnpfrom, 1, (Int4) *starttmp);
           }
           BioseqUnlock (bsp);
           if (strand_nonull)
              strandtmp = *strandp;
           slp = SeqLocIntNew (0, seqto, strandtmp, sip);
           if ( slp == NULL) {
              return NULL;
           }
           spp = SeqPortNewByLoc (slp, Seq_code_iupacna);
/**
           SeqIdWrite (sip, strLog, PRINTID_FASTA_LONG, 120);
           if (spp!=NULL) WriteLog ("mergesalp2 %s\n", strLog);
           else WriteLog ("spp NULL\n");
           WriteLog ("1>  %ld  %ld  from %ld  to %ld  offset %ld  \n", dsp_len, sum_lens, seqfrom, seqto, seqoffset);
**/
           if (seqoffset < 0) seqoffset = 0;  /*!!!!!!!!!!!!!!!!!!!!*/
           if ( j == 0 && seqoffset > seqfrom) {
              seqfrom = ReadBufferFromSep (spp, buffertmp[j], (Int4)seqfrom,
                                         (Int4)(seqfrom + seqoffset), 0);
           }
/**
           if ( j == 0)
             WriteLog ("2> %d   %ld  %ld \n", seqfrom, sum_lens, seqfrom + sum_lens);
           else
             WriteLog ("2> %d   %ld  %ld \n", seqfrom, sum_lens, seqfrom + sum_lens);
**/
           alglens = readbuff_fromseqalign (spp, salptmp, j, buffertmp[j], seqfrom, seqfrom + sum_lens, (seqoffset + gapoffset), strandtmp);
           if (alglens == 0) {
              return NULL;
           }
/**
           WriteLog ("3> %d  %c%c%c \n", j, buffertmp[j][0], buffertmp[j][1], buffertmp[j][2]);
           WriteLog ("4>  %ld     %ld %ld  %ld < %ld\n", alglens, sum_lens, sum_lens + seqfrom, seqfrom + dsp_len, bsp_len);
**/
           if ( j == 0 && (seqfrom + dsp_len) < bsp_len) {
              alglens = ReadBufferFromSep (spp, buffertmp[j], (Int4)(seqfrom + dsp_len), (Int4)bsp_len, alglens);
           }
           else { 
              for (k = alglens; k < step; k++) 
                 if (buffertmp[j][k] != '-') buffertmp[j][k] = '-';
           }
           SeqPortFree(spp);
           if (j>0)
              ValNodeAddInt (&vnpstrand, 1, (Int4)strandtmp);
           if (strand_nonull)
              strandp++;
       }
       str = MemNew ((size_t) ((step+10) * sizeof(CharPtr)));
       bufferout [curnseq] = str;
       for (j=0; j < curnseq + 1; j++) {
           for (k=0; k<step; k++) 
              bufferout [j][k] = '-';
       }  
       k1=k2=k3=letter=0;
       while ( k1 < step && k2 < step && k3 < step)
       {
           c1 = bufferin[0][k1];
           c2 = buffertmp[0][k2];
           if ((c1 != '-' && c2 != '-') || (c1 == '-' && c2 == '-')) {
              if (c1 != '-' && c2 != '-' && c1 != c2) {
/**
                 WriteLog ("ERROR %d [%c]  %d [%c] %c%c%c%c%c  %c%c%c%c%c", k1, c1, k2, c2, c1,bufferin[0][k1+1], bufferin[0][k1+2], bufferin[0][k1+3], bufferin[0][k1+4], c2,buffertmp[0][k2+1], buffertmp[0][k2+2], buffertmp[0][k2+3], buffertmp[0][k2+4]); 
**/
                 break;
              }
              for (j = 0; j < curnseq; ++j) {
                 bufferout[j][k3] = bufferin[j][k1];
              }
              for (j = 1; j < 2; ++j) {
                 bufferout[curnseq + j -1][k3] = buffertmp[j][k2];
              }
              if (bufferout[0][k3] !='-') letter++;
/***
              WriteLog ("%d %d>%c ", (int)k3, (int)letter, c2);
              for (j=0; j<curnseq + 1; j++) WriteLog ("%c", bufferout[j][k3]);
              WriteLog ("\n");
***/
              k1++;
              k2++;
              k3++;
           }
           else if (c1 == '-' && c2 != '-') {
              for (j = 0; j < curnseq; ++j) {
                 bufferout[j][k3] = bufferin[j][k1];
              }
              if (bufferout[0][k3] !='-') letter++;
/***
              WriteLog ("%d %d*%c ", (int)k3, (int)letter, c2);
              for (j=0; j<curnseq + 1; j++) WriteLog ("%c", bufferout[j][k3]);
              WriteLog ("\n");
***/
              k1++;
              k3++;
           }
           else if (c1 != '-' && c2 == '-') {
              for (j = 1; j < 2; ++j) {
                 bufferout[curnseq + j -1][k3] = buffertmp[j][k2];
              }
              if (bufferout[0][k3] !='-') letter++;
/***
              WriteLog ("%d %d<%c ", (int)k3, (int)letter,c2);
              for (j=0; j<curnseq + 1; j++) WriteLog ("%c", bufferout[j][k3]);
              WriteLog ("\n");
***/
              k2++;
              k3++;
           }
       }
       if (k3 > 0) {
           if (k3 > max_lens) 
              max_lens = k3;
           for (j=0; j < curnseq ; j++) {
              for (k=0; k<step; k++) {
                 bufferin[j][k] = bufferout[j][k];
              }
           }
           str = MemNew ((size_t) ((step+10) * sizeof(CharPtr)));
           bufferin[j] = str;
           for (k=0; k<step; k++) {
              bufferin[j][k] = bufferout[j][k] ;
           }
           curnseq++;
       }
      }
     }
  }
  MemFree (buffertmp[0]);
  MemFree (buffertmp[1]);
  for (j=0; j < curnseq ; j++) 
     MemFree (bufferout[j]);
  
  vnp = NULL;
  for (j=0; j < curnseq; j++) {
     ValNodeAddPointer (&vnp, 0, (Pointer)(bufferin[j])); 
  }
  sap = LocalAlignToSeqAnnotDimn (vnp, siphead, vnpfrom, curnseq, max_lens, vnpstrand, TRUE);
  for (j=0; j < curnseq ; j++) 
     MemFree (bufferin[j]);
  return sap;
}

extern SeqAlignPtr PairSeqAlign2MultiSeqAlign (SeqAlignPtr salp)
{
  SeqAnnotPtr sap;
  SeqAlignPtr salptmp = NULL;

  if (salp!=NULL) 
  {
   if (is_dim2seqalign (salp)) 
   {
     sap = multseqalign_from_pairseqalign (salp);
     if (sap!=NULL && sap->type==2) 
     {
        salptmp = (SeqAlignPtr) sap->data;
        sap->data = NULL;
        SeqAnnotFree (sap);
        salp = SeqAlignFree (salp);
     }
   }
   else salptmp = salp;
  }
  return salptmp;   
}


extern Int2 LIBCALLBACK MultSeqAlignFromPairSeqAlign (Pointer data)

{
  OMProcControlPtr  ompcp;
  SeqAnnotPtr       sap;
  SeqAlignPtr       salp;
  Uint2             entityID;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  switch (ompcp->input_itemtype) {
    case OBJ_SEQALIGN :
      break;
    case 0 :
      return OM_MSG_RET_ERROR;
    default :
      return OM_MSG_RET_ERROR;
  }
  if (ompcp->input_data == NULL) return OM_MSG_RET_ERROR;
  salp = ompcp->input_data;
  if (salp == NULL) return OM_MSG_RET_ERROR;
  if (is_dim2seqalign (salp))  {
     SortSeqAlign (&salp);
     sap = multseqalign_from_pairseqalign (salp);
  } else {
     sap = SeqAnnotNew ();
     sap->type = 2;
     sap->data = (Pointer) salp;
  }
  if (sap != NULL) {
    entityID = ObjMgrRegister (OBJ_SEQANNOT, (Pointer) sap);
    if (entityID > 0) {
      ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
    }
  }
  return OM_MSG_RET_DONE;
}


/********
*** TEMP FUNCTION
*** if the SeqAlign arrives as DenseSeg 
*** -> DenseDiag
************/



