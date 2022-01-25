/*   ent2api.c
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
* File Name:  ent2api.c
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

#include <ent2api.h>

#define ENTREZ_TOOL_PROPERTY "Entrez2Tool"
#define ENTREZ_TOOL_VERSION 1

NLM_EXTERN void EntrezSetProgramName (
  CharPtr progname
)

{
  CharPtr  str;

  str = (CharPtr) GetAppProperty (ENTREZ_TOOL_PROPERTY);
  if (str != NULL) {
    MemFree (str);
  }
  if (! StringHasNoText (progname)) {
    SetAppProperty (ENTREZ_TOOL_PROPERTY, StringSave (progname));
  } else {
    SetAppProperty (ENTREZ_TOOL_PROPERTY, NULL);
  }
}

static CharPtr EntrezGetProgramName (
  void
)

{
  Char     path [PATH_MAX];
  CharPtr  ptr;

  ptr = (CharPtr) GetAppProperty (ENTREZ_TOOL_PROPERTY);
  if (StringHasNoText (ptr)) {
    Nlm_ProgramPath (path, sizeof (path));
    ptr = StringRChr (path, DIRDELIMCHR);
    if (ptr != NULL) {
      ptr++;
      EntrezSetProgramName (ptr);
    }
  }
  return ptr;
}

static Entrez2RequestPtr CreateRequest (
  Uint1 choice, Pointer data
)

{
  Entrez2RequestPtr  e2rp;
  ValNodePtr         vnp;

  e2rp = Entrez2RequestNew ();
  if (e2rp == NULL) return NULL;

  e2rp->version = ENTREZ_TOOL_VERSION;
  e2rp->tool = StringSaveNoNull (EntrezGetProgramName ());

  vnp = ValNodeNew (NULL);
  if (vnp == NULL) return NULL;
  vnp->choice = choice;
  vnp->data.ptrvalue = data;
  vnp->next = NULL;

  e2rp->request = vnp;

  return e2rp;
}

NLM_EXTERN Entrez2IdListPtr EntrezCreateEntrezIdList (
  CharPtr db,
  Int4 uid,
  Int4 num,
  Int4Ptr uids,
  ByteStorePtr bs
)

{
  Entrez2IdListPtr  e2il;

  e2il = Entrez2IdListNew ();
  if (e2il == NULL) return NULL;

  e2il->db = StringSaveNoNull (db);

  if (uid != 0 && uids == NULL) {
    uids = &uid;
    num = 1;
  }

  if (uids != NULL && num > 0 && bs == NULL) {
    bs = BSNew (4 * num);
    if (bs == NULL) return NULL;
    BSWrite (bs, (Uint4Ptr) uids, num * sizeof (Uint4));
  }

  e2il->uids = (Pointer) bs;
  e2il->num = BSLen (bs) / sizeof (Uint4);

  return e2il;
}

NLM_EXTERN Entrez2RequestPtr EntrezCreateGetInfoRequest (
  void
)

{
  return CreateRequest (E2Request_get_info, NULL);
}

NLM_EXTERN Entrez2LimitsPtr EntrezCreateEntrezLimits (
  Int4 begin_date,
  Int4 end_date,
  CharPtr type_date,
  Int4 max_uids,
  Int4 offset_uids
)

{
  Entrez2DtFilterPtr  e2df;
  Entrez2LimitsPtr    e2lm;

  if (begin_date == 0 && end_date == 0 &&
      StringHasNoText (type_date) &&
      max_uids == 0 && offset_uids == 0) return NULL;

  e2lm = Entrez2LimitsNew ();
  if (e2lm == NULL) return NULL;

  e2lm->max_UIDs = max_uids;
  e2lm->offset_UIDs = offset_uids;

  if (begin_date == 0 && end_date == 0 &&
      StringHasNoText (type_date)) return e2lm;

  e2df = Entrez2DtFilterNew ();
  if (e2df == NULL) return NULL;

  e2df->begin_date = begin_date;
  e2df->end_date = end_date;
  e2df->type_date = StringSaveNoNull (type_date);

  e2lm->filter_date = e2df;

  return e2lm;
}

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
)

{
  Entrez2BooleanExpPtr   e2be;
  Entrez2EvalBooleanPtr  e2eb;
  Entrez2RequestPtr      e2rp;

  e2be = Entrez2BooleanExpNew ();
  if (e2be == NULL) return NULL;

  e2be->db = StringSaveNoNull (db);
  e2be->limits = EntrezCreateEntrezLimits (begin_date, end_date,
                                           type_date, max_uids, offset_uids);

  e2eb = Entrez2EvalBooleanNew ();
  if (e2eb == NULL) return NULL;

  e2eb->return_UIDs = return_uids;
  e2eb->return_parse = return_parsed;
  e2eb->query = e2be;

  e2rp = CreateRequest (E2Request_eval_boolean, (Pointer) e2eb);
  if (e2rp == NULL) return NULL;

  if (! StringHasNoText (query_string)) {
    EntrezAddToBooleanRequest (e2rp, query_string, 0, NULL, NULL, 0, 0,
                               NULL, NULL, do_not_explode, do_not_translate);
  }

  return e2rp;
}

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
)

{
  Entrez2BooleanExpPtr   e2be;
  Entrez2BooleanTermPtr  e2bt;
  Entrez2EvalBooleanPtr  e2eb;
  Entrez2IdListPtr       e2il;
  ValNodePtr             vnp;

  if (e2rp == NULL) return;
  vnp = e2rp->request;
  if (vnp == NULL || vnp->choice != E2Request_eval_boolean) return;

  e2eb = (Entrez2EvalBooleanPtr) vnp->data.ptrvalue;
  if (e2eb == NULL) return;

  e2be = e2eb->query;
  if (e2be == NULL) return;

  if (! StringHasNoText (query_string)) {
    ValNodeCopyStr (&(e2be->exp), Entrez2BooleanElement_str, query_string);

  } else if (op > 0) {
    ValNodeAddInt (&(e2be->exp), Entrez2BooleanElement_op, op);

  } else if ((! StringHasNoText (field)) && (! StringHasNoText (term))) {
    e2bt = Entrez2BooleanTermNew ();
    if (e2bt == NULL) return;

    e2bt->field = StringSaveNoNull (field);
    e2bt->term = StringSaveNoNull (term);
    e2bt->do_not_explode = do_not_explode;
    e2bt->do_not_translate = do_not_translate;

    ValNodeAddPointer (&(e2be->exp), Entrez2BooleanElement_term, (Pointer) e2bt);

  } else {

    e2il = EntrezCreateEntrezIdList (e2be->db, uid, num, uids, bs);
    if (e2il == NULL) return;

    ValNodeAddPointer (&(e2be->exp), Entrez2BooleanElement_ids, (Pointer) e2il);
  }
}

NLM_EXTERN Entrez2RequestPtr EntrezCreateDocSumRequest (
  CharPtr db,
  Int4 uid,
  Int4 num,
  Int4Ptr uids,
  ByteStorePtr bs
)

{
  Entrez2IdListPtr  e2il;

  e2il = EntrezCreateEntrezIdList (db, uid, num, uids, bs);
  if (e2il == NULL) return NULL;

  return CreateRequest (E2Request_get_docsum, (Pointer) e2il);
}

NLM_EXTERN Entrez2RequestPtr EntrezCreateGetTermPositionRequest (
  CharPtr db,
  CharPtr field,
  CharPtr term
)

{
  Entrez2TermQueryPtr  e2tq;

  e2tq = Entrez2TermQueryNew ();
  if (e2tq == NULL) return NULL;
  e2tq->db = StringSaveNoNull (db);
  e2tq->field = StringSaveNoNull (field);
  e2tq->term = StringSaveNoNull (term);

  return CreateRequest (E2Request_get_term_pos, (Pointer) e2tq);
}

NLM_EXTERN Entrez2RequestPtr EntrezCreateGetTermListRequest (
  CharPtr db,
  CharPtr field,
  Int4 first_term_pos,
  Int4 num_terms
)

{
  Entrez2TermPosPtr  e2tp;

  e2tp = Entrez2TermPosNew ();
  if (e2tp == NULL) return NULL;
  e2tp->db = StringSaveNoNull (db);
  e2tp->field = StringSaveNoNull (field);
  e2tp->first_term_pos = first_term_pos;
  e2tp->number_of_terms = num_terms;

  return CreateRequest (E2Request_get_term_list, (Pointer) e2tp);
}

NLM_EXTERN Entrez2RequestPtr EntrezCreateGetTermHierarchyRequest (
  CharPtr db,
  CharPtr field,
  CharPtr term
)

{
  Entrez2TermQueryPtr  e2tq;

  e2tq = Entrez2TermQueryNew ();
  if (e2tq == NULL) return NULL;
  e2tq->db = StringSaveNoNull (db);
  e2tq->field = StringSaveNoNull (field);
  e2tq->term = StringSaveNoNull (term);

  return CreateRequest (E2Request_get_term_hierarchy, (Pointer) e2tq);
}

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
)

{
  Entrez2GetLinksPtr  e2gl;
  Entrez2IdListPtr    e2il;

  e2il = EntrezCreateEntrezIdList (db, uid, num, uids, bs);
  if (e2il == NULL) return NULL;

  e2gl = Entrez2GetLinksNew ();
  if (e2gl == NULL) return NULL;

  e2gl->uids = e2il;
  e2gl->linktype = StringSaveNoNull (linktype);
  e2gl->max_UIDS = max_uids;
  e2gl->count_only = count_only;
  e2gl->parents_persist = parents_persist;

  return CreateRequest (E2Request_get_links, (Pointer) e2gl);
}

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
)

{
  Entrez2GetLinksPtr  e2gl;
  Entrez2IdListPtr    e2il;

  e2il = EntrezCreateEntrezIdList (db, uid, num, uids, bs);
  if (e2il == NULL) return NULL;

  e2gl = Entrez2GetLinksNew ();
  if (e2gl == NULL) return NULL;

  e2gl->uids = e2il;
  e2gl->linktype = StringSaveNoNull (linktype);
  e2gl->max_UIDS = max_uids;
  e2gl->count_only = count_only;
  e2gl->parents_persist = parents_persist;

  return CreateRequest (E2Request_get_linked, (Pointer) e2gl);
}

NLM_EXTERN Entrez2RequestPtr EntrezCreateGetLinkCountsRequest (
  CharPtr db,
  Int4 uid
)

{
  Entrez2IdPtr  e2id;

  e2id = Entrez2IdNew ();
  if (e2id == NULL) return NULL;

  e2id->db = StringSaveNoNull (db);
  e2id->uid = uid;

  return CreateRequest (E2Request_get_link_counts, (Pointer) e2id);
}








