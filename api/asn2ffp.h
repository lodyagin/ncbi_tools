/*    asn2ffp.h
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
* File Name:  asn2ffp.h
*
* Author:  Karl Sirotkin, Tom Madden, Tatiana Tatusov
*
* Version Creation Date:   7/15/95
*
* $Revision: 6.37 $
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

#ifndef _ASN2FFP_
#define _ASN2FFP_

#include <asn2ffg.h>
#include <asn2ff.h>
#include <asn2ff6.h>


#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

NLM_EXTERN  Boolean asn2ff_flags[13];

#define ASN2FF_LOCAL_ID                 asn2ff_flags[0]
#define ASN2FF_LOOK_FOR_SEQ             asn2ff_flags[1]
#define ASN2FF_VALIDATE_FEATURES        asn2ff_flags[2]
#define ASN2FF_IGNORE_PATENT_PUBS       asn2ff_flags[3]
#define ASN2FF_DROP_SHORT_AA            asn2ff_flags[4]
#define ASN2FF_AVOID_LOCUS_COLL         asn2ff_flags[5]
#define ASN2FF_DATE_ERROR_MSG           asn2ff_flags[6]
#define ASN2FF_IUPACAA_ONLY             asn2ff_flags[7]
#define ASN2FF_TRANSL_TABLE             asn2ff_flags[8]
#define ASN2FF_REPORT_LOCUS_COLL        asn2ff_flags[9]
#define ASN2FF_SHOW_ALL_PUBS	        asn2ff_flags[10]
#define ASN2FF_SHOW_ERROR_MSG	        asn2ff_flags[11]
#define ASN2FF_SHOW_GB_STYLE	        asn2ff_flags[12]

#define DBNUM 58
NLM_EXTERN CharPtr dbtag[DBNUM];

NLM_EXTERN void FlatSpliceOff PROTO((SeqEntryPtr the_set, ValNodePtr desc));
NLM_EXTERN void FlatSpliceOn PROTO((SeqEntryPtr the_set, ValNodePtr desc));

NLM_EXTERN void PrintLocusLine PROTO((Asn2ffJobPtr ajp, GBEntryPtr gbp));
NLM_EXTERN void PrintAccessLine PROTO((Asn2ffJobPtr ajp, GBEntryPtr gbp));
NLM_EXTERN void PrintVersionLine PROTO((Asn2ffJobPtr ajp, GBEntryPtr gbp));
NLM_EXTERN void PrintNCBI_GI PROTO((Asn2ffJobPtr ajp, GBEntryPtr gbp));
NLM_EXTERN void PrintNID PROTO((Asn2ffJobPtr ajp, GBEntryPtr gbp));
NLM_EXTERN void GetDefinitionLine PROTO((Asn2ffJobPtr ajp, GBEntryPtr gbp));
NLM_EXTERN void PrintDefinitionLine PROTO((Asn2ffJobPtr ajp, GBEntryPtr gbp));
NLM_EXTERN void PrintKeywordLine PROTO((Asn2ffJobPtr ajp, GBEntryPtr gbp));
NLM_EXTERN void PrintOriginLine PROTO ((Asn2ffJobPtr ajp, GBEntryPtr gbp));
NLM_EXTERN void PrintOrganismLine PROTO ((Asn2ffJobPtr ajp, GBEntryPtr gbp));

NLM_EXTERN void PrintEPLocusLine PROTO((Asn2ffJobPtr ajp, GBEntryPtr gbp));
NLM_EXTERN void PrintSegmentLine PROTO((Asn2ffJobPtr ajp, GBEntryPtr gbp));

NLM_EXTERN void PrintGBSourceLine PROTO ((Asn2ffJobPtr ajp, GBEntryPtr gbp));
NLM_EXTERN void PrintGBOrganismLine PROTO ((Asn2ffJobPtr ajp, GBEntryPtr gbp));

NLM_EXTERN void PrintPubsByNumber PROTO ((Asn2ffJobPtr ajp, GBEntryPtr gbp));
NLM_EXTERN void PrintFeatHeader PROTO ((Asn2ffJobPtr ajp, GBEntryPtr gbp));
NLM_EXTERN void PrintSequence PROTO ((Asn2ffJobPtr ajp, GBEntryPtr gbp, Int4 start, Int4 stop));
NLM_EXTERN void PrintEPSequence PROTO ((Asn2ffJobPtr ajp, GBEntryPtr gbp, Int4 start, Int4 stop));
NLM_EXTERN void PrintBaseCount PROTO ((Asn2ffJobPtr ajp, GBEntryPtr gbp));
NLM_EXTERN ValNodePtr tie_next PROTO((ValNodePtr head, ValNodePtr vnp));
NLM_EXTERN ValNodePtr GatherDescrByChoice PROTO ((Asn2ffJobPtr ajp, GBEntryPtr gbp, Uint1 choice));
NLM_EXTERN ValNodePtr GatherDescrListByChoice PROTO ((Asn2ffJobPtr ajp, GBEntryPtr gbp, Uint1 choice));
NLM_EXTERN ValNodePtr GetOrgRef PROTO ((Asn2ffJobPtr ajp, GBEntryPtr gbp));
NLM_EXTERN ValNodePtr GetBiosource PROTO ((Asn2ffJobPtr ajp, GBEntryPtr gbp));
NLM_EXTERN void EMBL_PrintPubs PROTO ((Asn2ffJobPtr ajp, GBEntryPtr gbp, PubStructPtr psp));
NLM_EXTERN void GB_PrintPubs PROTO ((Asn2ffJobPtr ajp, GBEntryPtr gbp, PubStructPtr psp));
NLM_EXTERN void GR_PrintPubs PROTO ((Asn2ffJobPtr ajp, GBEntryPtr gbp, PubStructPtr psp));
NLM_EXTERN Boolean FlatIgnoreThisPatentPub PROTO ((BioseqPtr bsp, ValNodePtr best, Int4Ptr seqidPt));
NLM_EXTERN CharPtr FlatCleanEquals PROTO ((CharPtr retval));
NLM_EXTERN ValNodePtr GetAuthors PROTO((Asn2ffJobPtr ajp, ValNodePtr the_pub));
NLM_EXTERN CharPtr FlatJournal PROTO ((Asn2ffJobPtr ajp, GBEntryPtr gbp, ValNodePtr the_pub, Int4 pat_seqid, Boolean PNTR submit, Boolean make_index));
NLM_EXTERN ValNodePtr GetKeywordLine PROTO((Asn2ffJobPtr ajp, GBEntryPtr gbp));
NLM_EXTERN void PrintSourceFeat PROTO ((Asn2ffJobPtr ajp, GBEntryPtr gbp));
NLM_EXTERN Int2 PrintImpFeat PROTO ((Asn2ffJobPtr ajp, BioseqPtr bsp, SeqFeatPtr sfp));
NLM_EXTERN Int2 PrintImpFeatEx PROTO ((Asn2ffJobPtr ajp, BioseqPtr bsp, SeqFeatPtr sfp, BIG_ID gi, Int2 entityID, Uint4 itemID));
NLM_EXTERN void PrintNAFeatAwp PROTO ((Asn2ffJobPtr ajp, GBEntryPtr gbp));
NLM_EXTERN void PrintNAFeatByNumber PROTO ((Asn2ffJobPtr ajp, GBEntryPtr gbp));
NLM_EXTERN void PrintAAFeatByNumber PROTO ((Asn2ffJobPtr ajp, GBEntryPtr gbp));
NLM_EXTERN CharPtr FlatLoc PROTO ((BioseqPtr bsp, ValNodePtr location));
NLM_EXTERN Boolean FlatAnnotPartial PROTO ((SeqFeatPtr sfp, Boolean use_product));
NLM_EXTERN Boolean FlatIgnoreThisPatentPub PROTO ((BioseqPtr bsp, ValNodePtr best, Int4Ptr seqidPt));

NLM_EXTERN void PrintCommentByNumber PROTO((Asn2ffJobPtr aip, GBEntryPtr gbp));
NLM_EXTERN void PrintFirstComment PROTO((Asn2ffJobPtr aip, GBEntryPtr gbp));
NLM_EXTERN void GBDescrComFeat PROTO((Asn2ffJobPtr ajp, GBEntryPtr gbp));

NLM_EXTERN Int2 GB_GetSeqDescrComms  PROTO((Asn2ffJobPtr ajp, GBEntryPtr gbp));
NLM_EXTERN Int2 GP_GetSeqDescrComms  PROTO((Asn2ffJobPtr ajp, GBEntryPtr gbp));

NLM_EXTERN Int4 GetGibbsqNumber PROTO ((BioseqPtr bsp));
NLM_EXTERN Int4 GetGibbsqCommentLength PROTO ((GBEntryPtr gbp));
NLM_EXTERN CharPtr GetGibbsqComment PROTO ((GBEntryPtr gbp));
NLM_EXTERN Int4 GetGibbsqStatement PROTO ((GBEntryPtr gbp, CharPtr ptr));
NLM_EXTERN Int2 StorePubInfo PROTO ((Asn2ffJobPtr ajp, BioseqContextPtr bcp, BioseqPtr bsp, ValNodePtr PNTR vnpp));
NLM_EXTERN Int2 StoreNAPubInfo PROTO ((Asn2ffJobPtr ajp, GBEntryPtr gbp, ValNodePtr PNTR vnp, Boolean error_msgs));
NLM_EXTERN void GetGBDate PROTO ((Asn2ffJobPtr ajp, GBEntryPtr gbp));
NLM_EXTERN void GetGPDate PROTO ((Asn2ffJobPtr ajp, GBEntryPtr gbp));
NLM_EXTERN void GetEMBLDate PROTO ((Asn2ffJobPtr ajp, GBEntryPtr gbp));
NLM_EXTERN void GetEntryVersion PROTO ((GBEntryPtr gbp));
NLM_EXTERN Boolean GetGeneQuals PROTO ((SeqFeatPtr sfp_in, GeneStructPtr gsp));
NLM_EXTERN Boolean GetCdregionGeneXrefInfo PROTO ((Asn2ffJobPtr ajp, SeqFeatPtr sfp, GBEntryPtr gbp, Int2 index));
NLM_EXTERN void GetGeneRefInfo PROTO ((GeneStructPtr gsp, NoteStructPtr nsp, GeneRefPtr grp));
NLM_EXTERN void GetDBXrefFromGene PROTO ((GeneRefPtr grp, SeqFeatPtr sfp));
NLM_EXTERN Int2 CompareStringWithGsp PROTO ((GeneStructPtr gsp, CharPtr string));
NLM_EXTERN Boolean CheckNAFeat PROTO ((Boolean is_new, BioseqPtr bsp, SeqFeatPtr sfp));
NLM_EXTERN Boolean CheckAndGetNAFeatLoc PROTO ((BioseqPtr bsp, CharPtr PNTR buffer, SeqFeatPtr sfp, Boolean loc_return));
NLM_EXTERN void GetAAFeatLoc PROTO ((BioseqPtr bsp, CharPtr PNTR buffer, SeqFeatPtr sfp, Boolean use_product));
NLM_EXTERN CharPtr GetGBSourceLine PROTO ((GBBlockPtr gb));

NLM_EXTERN Int2 CheckPubs PROTO ((Asn2ffJobPtr ajp, BioseqPtr bsp, ValNodePtr PNTR vnpp));
NLM_EXTERN CharPtr FlatAuthor PROTO ((Asn2ffJobPtr ajp, ValNodePtr the_pub));
NLM_EXTERN CharPtr FlatPubTitle PROTO ((ValNodePtr the_pub));
NLM_EXTERN void PrintDBSourceLine PROTO ((Asn2ffJobPtr ajp, GBEntryPtr gbp));

NLM_EXTERN void PostARefErrMessage PROTO ((Asn2ffJobPtr ajp, BioseqPtr bsp, PubStructPtr psp, ValNodePtr ext_pub, Int2 status, CharPtr string));

NLM_EXTERN void SeparatePartSuppl PROTO((CharPtr vol_issue, CharPtr part_sub));
NLM_EXTERN void AddExtraAccessions PROTO ((Asn2ffJobPtr ajp, GBEntryPtr gbp));
NLM_EXTERN void PrintTerminator PROTO ((void));
NLM_EXTERN Boolean get_pubs PROTO ((GatherContextPtr gcp));
NLM_EXTERN void GatherItemWithLock PROTO((Uint2 entityID, Uint4 itemID, Uint2 itemtype,
                                   Pointer userdata, GatherItemProc userfunc));
NLM_EXTERN CharPtr format_article PROTO ((Asn2ffJobPtr ajp, BioseqPtr bsp, ValNodePtr the_pub, Boolean make_index));
NLM_EXTERN CharPtr format_bookarticle PROTO ((Asn2ffJobPtr ajp, BioseqPtr bsp, ValNodePtr the_pub, Boolean make_index));
NLM_EXTERN CharPtr format_jourarticle PROTO ((Asn2ffJobPtr ajp, BioseqPtr bsp, ValNodePtr the_pub, Boolean make_index));

NLM_EXTERN void GetProtRefInfo PROTO ((Uint1 format, GeneStructPtr gsp, NoteStructPtr nsp, ProtRefPtr prp));
NLM_EXTERN Int2 GetMapFeats PROTO((Asn2ffJobPtr ajp, GBEntryPtr gbp));
NLM_EXTERN Boolean find_item PROTO ((GatherContextPtr gcp));
NLM_EXTERN Boolean get_prot_feats PROTO((GatherContextPtr gcp));
NLM_EXTERN void AddSiteNoteQual PROTO((SeqFeatPtr sfp_in, SeqFeatPtr sfp));
NLM_EXTERN void MatchAAGeneToFeat PROTO((OrganizeFeatPtr ofp, SortStructPtr p));
NLM_EXTERN void MatchNAGeneToFeat PROTO ((Boolean non_strict, OrganizeFeatPtr ofp, SortStructPtr p));
NLM_EXTERN void SortOrganizeFeat PROTO((OrganizeFeatPtr ofp));
NLM_EXTERN void OrganizeSeqFeat PROTO ((Asn2ffJobPtr ajp, GBEntryPtr gbp));
NLM_EXTERN void GetSeqFeat PROTO ((Asn2ffJobPtr ajp));

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif
