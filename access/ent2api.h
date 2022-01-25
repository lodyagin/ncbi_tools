/*   ent2api.h
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
* File Name:  ent2api.h
*
* Author:  Jonathan Kans
*
* Version Creation Date:   7/29/99
*
* $Revision: 1.3 $
*
* File Description: 
*
* Modifications:  
* --------------------------------------------------------------------------
*
* ==========================================================================
*/

#ifndef _ENT2API_
#define _ENT2API_

#include <ncbi.h>
#include <asn.h>
#include <objent2.h>

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif


#ifdef __cplusplus
extern "C" {
#endif


NLM_EXTERN void EntrezSetProgramName (
  CharPtr progname
);

NLM_EXTERN Entrez2IdListPtr EntrezCreateEntrezIdList (
  CharPtr db,
  Int4 uid,
  Int4 num,
  Int4Ptr uids,
  ByteStorePtr bs
);

NLM_EXTERN Entrez2LimitsPtr EntrezCreateEntrezLimits (
  Int4 begin_date,
  Int4 end_date,
  CharPtr type_date,
  Int4 max_uids,
  Int4 offset_uids
);

NLM_EXTERN Entrez2RequestPtr EntrezCreateGetInfoRequest (
  void
);

NLM_EXTERN Entrez2RequestPtr EntrezCreateBooleanRequest (
  Boolean return_uids,
  Boolean return_parsed,
  CharPtr db,
  CharPtr query_string,
  Boolean do_not_explode,
  Boolean do_not_translate,
  Int4 begin_date,
  Int4 end_date,
  CharPtr type_date,
  Int4 max_uids,
  Int4 offset_uids
);

#define ENTREZ_OP_AND         1
#define ENTREZ_OP_OR          2
#define ENTREZ_OP_BUTNOT      3
#define ENTREZ_OP_RANGE       4
#define ENTREZ_OP_LEFT_PAREN  5
#define ENTREZ_OP_RIGHT_PAREN 6

NLM_EXTERN void EntrezAddToBooleanRequest (
  Entrez2RequestPtr e2rp,
  CharPtr query_string,
  Int4 op,
  CharPtr field,
  CharPtr term,
  Int4 uid,
  Int4 num,
  Int4Ptr uids,
  ByteStorePtr bs,
  Boolean do_not_explode,
  Boolean do_not_translate
);

NLM_EXTERN Entrez2RequestPtr EntrezCreateDocSumRequest (
  CharPtr db,
  Int4 uid,
  Int4 num,
  Int4Ptr uids,
  ByteStorePtr bs
);

NLM_EXTERN Entrez2RequestPtr EntrezCreateGetTermPositionRequest (
  CharPtr db,
  CharPtr field,
  CharPtr term
);

NLM_EXTERN Entrez2RequestPtr EntrezCreateGetTermListRequest (
  CharPtr db,
  CharPtr field,
  Int4 first_term_pos,
  Int4 num_terms
);

NLM_EXTERN Entrez2RequestPtr EntrezCreateGetTermHierarchyRequest (
  CharPtr db,
  CharPtr field,
  CharPtr term
);

NLM_EXTERN Entrez2RequestPtr EntrezCreateGetLinksRequest (
  CharPtr db,
  Int4 uid,
  Int4 num,
  Int4Ptr uids,
  ByteStorePtr bs,
  CharPtr linktype,
  Int4 max_uids,
  Boolean count_only,
  Boolean parents_persist
);

NLM_EXTERN Entrez2RequestPtr EntrezCreateGetLinkedRequest (
  CharPtr db,
  Int4 uid,
  Int4 num,
  Int4Ptr uids,
  ByteStorePtr bs,
  CharPtr linktype,
  Int4 max_uids,
  Boolean count_only,
  Boolean parents_persist
);

NLM_EXTERN Entrez2RequestPtr EntrezCreateGetLinkCountsRequest (
  CharPtr db,
  Int4 uid
);


#ifdef __cplusplus
}
#endif

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif /* _ENT2API_ */