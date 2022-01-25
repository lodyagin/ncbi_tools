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
* $Revision: 1.24 $
*
* File Description: 
*
* Modifications:  
* --------------------------------------------------------------------------
*
* ==========================================================================
*/

#include <ent2api.h>
#include <urlquery.h>

#define ENTREZ_TOOL_PROPERTY "Entrez2Tool"
#define ENTREZ_TOOL_VERSION 1

/* utility functions */

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
      ptr = (CharPtr) GetAppProperty (ENTREZ_TOOL_PROPERTY);
    }
  }
  return ptr;
}

/* override URL paths or service name */

static CharPtr  e2_host_machine = NULL;
static Uint2    e2_host_port = 0;
static CharPtr  e2_host_path = NULL;

static CharPtr  e2_service = NULL;

/* to be phased out, currently overrides ncbi named service, goes directly to URL */

NLM_EXTERN void EntrezSetServer (
  CharPtr host_machine,
  Uint2 host_port,
  CharPtr host_path
)

{
  if (! StringHasNoText (host_machine)) {
    e2_host_machine = MemFree (e2_host_machine);
    e2_host_machine = StringSaveNoNull (host_machine);
  }
  if (host_port != 0) {
    e2_host_port = host_port;
  }
  if (! StringHasNoText (host_path)) {
    e2_host_path = MemFree (e2_host_path);
    e2_host_path = StringSaveNoNull (host_path);
  }
}

/* use EntrezTest to override default Entrez ncbi named service */

NLM_EXTERN void EntrezSetService (
  CharPtr service
)

{
  if (! StringHasNoText (service)) {
    e2_service = MemFree (e2_service);
    e2_service = StringSaveNoNull (service);
  }
}

/* low-level connection functions */

NLM_EXTERN CONN EntrezOpenConnection (
  void
)

{
  CharPtr  host_machine = e2_host_machine;
  Uint2    host_port = e2_host_port;
  CharPtr  host_path = e2_host_path;
  CharPtr  host_service = e2_service;

  if (StringHasNoText (host_machine)) {
    host_machine = "www.ncbi.nlm.nih.gov";
  }
  if (host_port == 0) {
    host_port = 80;
  }
  if (StringHasNoText (host_path)) {
    host_path = "/entrez/utils/entrez2server.fcgi";
  }
  if (StringHasNoText (host_service)) {
    host_service = "Entrez2";
  }

  /* temporarily for direct access to URL, to be phased out */

  if (e2_host_machine != NULL || e2_host_port != 0 || e2_host_path != NULL) {
    return QUERY_OpenUrlQuery (host_machine, host_port, host_path,
                               NULL, EntrezGetProgramName (),
                               30, eMIME_T_NcbiData, eMIME_AsnBinary,
                               eENCOD_None, 0);
  }

  /* use new named service, Entrez or EntrezTest */

  return QUERY_OpenServiceQuery (host_service, NULL, 30);
}

#ifdef OS_MAC
#include <Events.h>
#endif

NLM_EXTERN Entrez2ReplyPtr EntrezWaitForReply (
  CONN conn
)

{
  AsnIoConnPtr     aicp;
  time_t           currtime, starttime;
  Entrez2ReplyPtr  e2ry = NULL;
  Int2             max = 0;
  EIO_Status       status;
  STimeout         timeout;
#ifdef OS_MAC
  EventRecord      currEvent;
#endif

  if (conn == NULL) return NULL;

#ifdef OS_MAC
  timeout.sec = 0;
  timeout.usec = 0;
#else
  timeout.sec = 100;
  timeout.usec = 0;
#endif

  starttime = GetSecs ();
  while ((status = CONN_Wait (conn, eIO_Read, &timeout)) != eIO_Success && max < 300) {
    currtime = GetSecs ();
    max = currtime - starttime;
#ifdef OS_MAC
    WaitNextEvent (0, &currEvent, 0, NULL);
#endif
  }
  if (status == eIO_Success) {
    aicp = QUERY_AsnIoConnOpen ("rb", conn);
    e2ry = Entrez2ReplyAsnRead (aicp->aip, NULL);
    QUERY_AsnIoConnClose (aicp);
  }
  CONN_Close (conn);

  return e2ry;
}

/* ent2api silently maintains entrez2 server session cookie */

static CharPtr e2cookie = NULL;

/* high-level connection functions */

NLM_EXTERN Entrez2ReplyPtr EntrezSynchronousQuery (
  Entrez2RequestPtr e2rq
)

{
  AsnIoConnPtr     aicp;
  CONN             conn;
  Entrez2ReplyPtr  e2ry;
  CharPtr          tempcookie = NULL;

  if (e2rq == NULL) return NULL;

  conn = EntrezOpenConnection ();

  aicp = QUERY_AsnIoConnOpen ("wb", conn);

  tempcookie = e2rq->cookie;
  if (e2rq->cookie == NULL && e2cookie != NULL) {
    e2rq->cookie = e2cookie;
  }

  Entrez2RequestAsnWrite (e2rq, aicp->aip, NULL);

  e2rq->cookie = tempcookie;

  AsnIoFlush (aicp->aip);
  QUERY_AsnIoConnClose (aicp);

  QUERY_SendQuery (conn);

  e2ry = EntrezWaitForReply (conn);

  if (e2ry != NULL && e2cookie == NULL && e2ry->cookie != NULL) {
    e2cookie = StringSave (e2ry->cookie);
  }

  return e2ry;
}

NLM_EXTERN Boolean EntrezAsynchronousQuery (
  Entrez2RequestPtr e2rq,
  QUEUE* queue,
  QueryResultProc resultproc,
  VoidPtr userdata
)

{
  AsnIoConnPtr  aicp;
  CONN          conn;
  CharPtr       tempcookie = NULL;

  if (e2rq == NULL) return FALSE;

  conn = EntrezOpenConnection ();

  aicp = QUERY_AsnIoConnOpen ("wb", conn);

  tempcookie = e2rq->cookie;
  if (e2rq->cookie == NULL && e2cookie != NULL) {
    e2rq->cookie = e2cookie;
  }

  Entrez2RequestAsnWrite (e2rq, aicp->aip, NULL);

  e2rq->cookie = tempcookie;

  AsnIoFlush (aicp->aip);
  QUERY_AsnIoConnClose (aicp);

  QUERY_SendQuery (conn);

  QUERY_AddToQueue (queue, conn, resultproc, userdata, TRUE);

  return TRUE;
}

NLM_EXTERN Int4 EntrezCheckQueue (QUEUE* queue)

{
  return QUERY_CheckQueue (queue);
}

NLM_EXTERN Entrez2ReplyPtr EntrezReadReply (
  CONN conn,
  EIO_Status status
)

{
  AsnIoConnPtr     aicp;
  Entrez2ReplyPtr  e2ry = NULL;

  if (conn != NULL && status == eIO_Success) {
    aicp = QUERY_AsnIoConnOpen ("rb", conn);
    e2ry = Entrez2ReplyAsnRead (aicp->aip, NULL);
    QUERY_AsnIoConnClose (aicp);
  }

  if (e2ry != NULL && e2cookie == NULL && e2ry->cookie != NULL) {
    e2cookie = StringSave (e2ry->cookie);
  }

  return e2ry;
}

/* request creation functions */

static Entrez2RequestPtr CreateRequest (
  Uint1 choice, Pointer data
)

{
  Entrez2RequestPtr  e2rq;
  ValNodePtr         vnp;

  e2rq = Entrez2RequestNew ();
  if (e2rq == NULL) return NULL;

  e2rq->version = ENTREZ_TOOL_VERSION;
  e2rq->tool = StringSaveNoNull (EntrezGetProgramName ());

  vnp = ValNodeNew (NULL);
  if (vnp == NULL) return NULL;
  vnp->choice = choice;
  vnp->data.ptrvalue = data;
  vnp->next = NULL;

  e2rq->request = vnp;

  return e2rq;
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
  Int4 begin_date,
  Int4 end_date,
  CharPtr type_date,
  Int4 max_uids,
  Int4 offset_uids
)

{
  Entrez2BooleanExpPtr   e2be;
  Entrez2EvalBooleanPtr  e2eb;
  Entrez2RequestPtr      e2rq;

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

  e2rq = CreateRequest (E2Request_eval_boolean, (Pointer) e2eb);
  if (e2rq == NULL) return NULL;

  if (! StringHasNoText (query_string)) {
    EntrezAddToBooleanRequest (e2rq, query_string, 0, NULL, NULL, 0, 0,
                               NULL, NULL, TRUE, TRUE);
  }

  return e2rq;
}

NLM_EXTERN void EntrezAddToBooleanRequest (
  Entrez2RequestPtr e2rq,
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

  if (e2rq == NULL) return;
  vnp = e2rq->request;
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
  CharPtr term,
  Int4 txid
)

{
  Entrez2HierQueryPtr  e2hq;

  e2hq = Entrez2HierQueryNew ();
  if (e2hq == NULL) return NULL;
  e2hq->db = StringSaveNoNull (db);
  e2hq->field = StringSaveNoNull (field);
  e2hq->term = StringSaveNoNull (term);
  e2hq->txid = txid;

  return CreateRequest (E2Request_get_term_hierarchy, (Pointer) e2hq);
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

/* reply extraction functions */

static Pointer GeneralEntrezExtractReply (
  Entrez2ReplyPtr e2ry,
  Uint1 choice,
  Int4Ptr termpos
)

{
  E2ReplyPtr  reply;
  Pointer     result = NULL;

  if (e2ry == NULL) return NULL;
  reply = e2ry->reply;
  if (reply == NULL) return NULL;

  if (reply->choice == choice) {
    if (termpos != NULL) {
      *termpos = reply->data.intvalue;
    } else {
      result = (Pointer) reply->data.ptrvalue;
      reply->data.ptrvalue = NULL;
    }
  }
  Entrez2ReplyFree (e2ry);

  return result;
}

NLM_EXTERN Entrez2InfoPtr EntrezExtractInfoReply (
  Entrez2ReplyPtr e2ry
)

{
  return (Entrez2InfoPtr) GeneralEntrezExtractReply (e2ry, E2Reply_get_info, NULL);
}

NLM_EXTERN Entrez2BooleanReplyPtr EntrezExtractBooleanReply (
  Entrez2ReplyPtr e2ry
)

{
  return (Entrez2BooleanReplyPtr) GeneralEntrezExtractReply (e2ry, E2Reply_eval_boolean, NULL);
}

NLM_EXTERN Entrez2DocsumListPtr EntrezExtractDocsumReply (
  Entrez2ReplyPtr e2ry
)

{
  return (Entrez2DocsumListPtr) GeneralEntrezExtractReply (e2ry, E2Reply_get_docsum, NULL);
}

NLM_EXTERN Int4 EntrezExtractTermPosReply (
  Entrez2ReplyPtr e2ry
)

{
  Int4  termpos = 0;

  GeneralEntrezExtractReply (e2ry, E2Reply_get_term_pos, &termpos);
  return termpos;
}

NLM_EXTERN Entrez2TermListPtr EntrezExtractTermListReply (
  Entrez2ReplyPtr e2ry
)

{
  return (Entrez2TermListPtr) GeneralEntrezExtractReply (e2ry, E2Reply_get_term_list, NULL);
}

NLM_EXTERN Entrez2HierNodePtr EntrezExtractHierNodeReply (
  Entrez2ReplyPtr e2ry
)

{
  return (Entrez2HierNodePtr) GeneralEntrezExtractReply (e2ry, E2Reply_get_term_hierarchy, NULL);
}

NLM_EXTERN Entrez2LinkSetPtr EntrezExtractLinksReply (
  Entrez2ReplyPtr e2ry
)

{
  return (Entrez2LinkSetPtr) GeneralEntrezExtractReply (e2ry, E2Reply_get_links, NULL);
}

NLM_EXTERN Entrez2IdListPtr EntrezExtractLinkedReply (
  Entrez2ReplyPtr e2ry
)

{
  return (Entrez2IdListPtr) GeneralEntrezExtractReply (e2ry, E2Reply_get_linked, NULL);
}
NLM_EXTERN Entrez2LinkCountListPtr EntrezExtractLinkCountReply (
  Entrez2ReplyPtr e2ry
)

{
  return (Entrez2LinkCountListPtr) GeneralEntrezExtractReply (e2ry, E2Reply_get_link_counts, NULL);
}

/* special SeqIdString to UID convenience function */

NLM_EXTERN Uint4 EntrezGetUIDforSeqIdString (
  CharPtr db,
  CharPtr seq_id_string
)

{
  Char                    ch;
  Entrez2BooleanReplyPtr  e2br;
  Entrez2IdListPtr        e2id;
  Entrez2RequestPtr       e2rq;
  Entrez2ReplyPtr         e2ry;
  CharPtr                 ptr;
  Char                    str [61];
  Uint4                   uid = 0;

  if (StringHasNoText (db) || StringHasNoText (seq_id_string)) return 0;

  StringNCpy_0 (str, seq_id_string, sizeof (str));
  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    if (ch == '|' || ch == '.') {
      *ptr = ' ';
    }
    ptr++;
    ch = *ptr;
  }
  TrimSpacesAroundString (str);
  if (StringStr (str, "[SQID]") == NULL) {
    StringCat (str, " [SQID]");
  }

  e2rq = EntrezCreateBooleanRequest (TRUE, FALSE, db, str,
                                     0, 0, NULL, 1, 0);
  if (e2rq == NULL) return 0;
  e2ry = EntrezSynchronousQuery (e2rq);
  e2rq = Entrez2RequestFree (e2rq);
  if (e2ry == NULL) return 0;
  e2br = EntrezExtractBooleanReply (e2ry);
  if (e2br == NULL) return 0;

  if (e2br->count > 0) {
    e2id = e2br->uids;
    if (e2id != NULL && e2id->num > 0 && e2id->uids != NULL) {
      BSSeek (e2id->uids, 0, SEEK_SET);
      uid = Nlm_BSGetUint4 (e2id->uids);
    }
  }

  Entrez2BooleanReplyFree (e2br);

  return uid;
}

