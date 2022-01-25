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
* $Revision: 1.154 $
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
#include <ncbithr.h>
#include <sqnutils.h>

#ifdef OS_UNIX
#include <sys/times.h>
#include <limits.h>
#endif

/* EUtilities replacement functions */

static CharPtr  eutils_host = "eutils.ncbi.nlm.nih.gov";

static CharPtr  einfo_url = "/entrez/eutils/einfo.fcgi";
static CharPtr  elink_url = "/entrez/eutils/elink.fcgi";
static CharPtr  esearch_url = "/entrez/eutils/esearch.fcgi";
static CharPtr  esummary_url = "/entrez/eutils/esummary.fcgi";

static Char _ToKey[256] = {
    0x00, 0x01, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20,

    0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20,

    0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20,

    0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x1F,
  /*  sp     !     "     #     $     %     &     ' */
    0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x27,
  /*   (     )     *     +     ,     -     .     / */
    0x20, 0x20, 0x20, 0x20, 0x2C, 0x20, 0x20, 0x2F,
  /*   0     1     2     3     4     5     6     7 */
    0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37,
  /*   8     9     :     ;     <     =     >     ? */
    0x38, 0x39, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20,
  /*   @     A     B     C     D     E     F     G */
    0x40, 0x61, 0x62, 0x63, 0x64, 0x65, 0x66, 0x67,
  /*   H     I     J     K     L     M     N     O */
    0x68, 0x69, 0x6A, 0x6B, 0x6C, 0x6D, 0x6E, 0x6F,
  /*   P     Q     R     S     T     U     V     W */
    0x70, 0x71, 0x72, 0x73, 0x74, 0x75, 0x76, 0x77,
  /*   X     Y     Z     [     \     ]     ^     _ */
    0x78, 0x79, 0x7A, 0x20, 0x20, 0x20, 0x20, 0x20,
  /*   `     a     b     c     d     e     f     g */
    0x60, 0x61, 0x62, 0x63, 0x64, 0x65, 0x66, 0x67,
  /*   h     i     j     k     l     m     n     o */
    0x68, 0x69, 0x6A, 0x6B, 0x6C, 0x6D, 0x6E, 0x6F,
  /*   p     q     r     s     t     u     v     w */
    0x70, 0x71, 0x72, 0x73, 0x74, 0x75, 0x76, 0x77,
  /*   x     y     z     {     |     }     ~   DEL */
    0x78, 0x79, 0x7A, 0x20, 0x20, 0x20, 0x20, 0x20,

    0x80, 0x81, 0x82, 0x83, 0x84, 0x85, 0x86, 0x87,

    0x88, 0x89, 0x8A, 0x8B, 0x8C, 0x8D, 0x8E, 0x8F,

    0x90, 0x91, 0x92, 0x93, 0x94, 0x95, 0x96, 0x97,

    0x98, 0x99, 0x9A, 0x9B, 0x9C, 0x9D, 0x9E, 0x9F,

    0xA0, 0xA1, 0xA2, 0xA3, 0xA4, 0xA5, 0xA6, 0xA7,

    0xA8, 0xA9, 0xAA, 0xAB, 0xAC, 0xAD, 0xAE, 0xAF,

    0xB0, 0xB1, 0xB2, 0xB3, 0xB4, 0xB5, 0xB6, 0xB7,

    0xB8, 0xB9, 0xBA, 0xBB, 0xBC, 0xBD, 0xBE, 0xBF,

    0xC0, 0xC1, 0xC2, 0xC3, 0xC4, 0xC5, 0xC6, 0xC7,

    0xC8, 0xC9, 0xCA, 0xCB, 0xCC, 0xCD, 0xCE, 0xCF,

    0xD0, 0xD1, 0xD2, 0xD3, 0xD4, 0xD5, 0xD6, 0xD7,

    0xD8, 0xD9, 0xDA, 0xDB, 0xDC, 0xDD, 0xDE, 0xDF,

    0xE0, 0xE1, 0xE2, 0xE3, 0xE4, 0xE5, 0xE6, 0xE7,

    0xE8, 0xE9, 0xEA, 0xEB, 0xEC, 0xED, 0xEE, 0xEF,

    0xF0, 0xF1, 0xF2, 0xF3, 0xF4, 0xF5, 0xF6, 0xF7,

    0xF8, 0xF9, 0xFA, 0xFB, 0xFC, 0xFD, 0xFE, 0xFF
};

NLM_EXTERN void ConvertToTermListForm (
  CharPtr str
)

{
  Char     ch;
  CharPtr  ptr;

  if (str == NULL) return;

  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    *ptr = _ToKey [(int) ch];
    ptr++;
    ch = *ptr;
  }

  TrimSpacesAroundString (str);
  CompressSpaces (str);
}

NLM_EXTERN Boolean LooksLikeISSN (
  CharPtr str
)

{
  Char  ch;
  Int2  i;

  if (StringHasNoText (str)) return FALSE;

  if (StringLen (str) != 9) return FALSE;
  ch = str [4];
  if (ch != '-' && ch != ' ' && ch != '+') return FALSE;

  for (i = 0; i < 9; i++) {
    ch = str [i];
    if (IS_DIGIT (ch)) continue;
    if (i == 4) {
      if (ch == '-' || ch == '+' || ch == ' ') continue;
    }
    if (i == 8) {
      if (ch == 'X' || ch == 'x') continue;
    }
    return FALSE;
  }

  return TRUE;
}

NLM_EXTERN CharPtr UrlEncodeString (
  CharPtr str
)

{
  size_t   arg_len, res_len, rd_len, wr_len;
  CharPtr  res;

  if (StringHasNoText (str)) return NULL;

  arg_len = StringLen (str);
  res_len = 3 * arg_len;

  res = (CharPtr) MemNew ((res_len + 2) * sizeof (Char));
  if (res == NULL) return NULL;

  URL_Encode (str, arg_len, &rd_len, res, res_len, &wr_len);

  if (arg_len != rd_len) {
    res = MemFree (res);
  }

  return res;
}

NLM_EXTERN CharPtr DoEinfoQuery (
  CharPtr database
)

{
  CharPtr  db, ptr = NULL, str;
  size_t   len;

  if (StringDoesHaveText (database)) {
    db = UrlEncodeString (database);

    len = StringLen (db);
    if (len < 1) return NULL;

    ptr = (CharPtr) MemNew (len + 50);
    if (ptr == NULL) return NULL;

    StringCpy (ptr, "db=");
    StringCat (ptr, db);

    MemFree (db);
  }

  str = QUERY_UrlSynchronousQuery (eutils_host, 80, einfo_url, ptr, NULL, NULL, NULL);

  MemFree (ptr);

  return str;
}

NLM_EXTERN CharPtr DoEsearchQuery (
  CharPtr database,
  CharPtr query,
  CharPtr suffix,
  CharPtr webenv,
  Int4 retstart,
  Int4 retmax,
  Boolean return_UIDs
)

{
  Char     buf [64];
  CharPtr  db, qu, sf, ptr, str;
  size_t   len;

  if (StringHasNoText (database) || StringHasNoText (query)) return NULL;

  db = UrlEncodeString (database);
  qu = UrlEncodeString (query);
  sf = UrlEncodeString (suffix);

  len = StringLen (db) + StringLen (qu) + StringLen (sf) + StringLen (webenv);
  if (len < 1) return NULL;

  ptr = (CharPtr) MemNew (len + 200);
  if (ptr == NULL) return NULL;

  StringCpy (ptr, "db=");
  StringCat (ptr, db);
  StringCat (ptr, "&term=");
  StringCat (ptr, qu);
  StringCat (ptr, sf);
  if (! return_UIDs) {
    StringCat (ptr, "&rettype=count");
  }
  StringCat (ptr, "&usehistory=y");
  if (StringDoesHaveText (webenv)) {
    StringCat (ptr, "&WebEnv=");
    StringCat (ptr, webenv);
  }
  if (retstart > 1) {
    sprintf (buf, "&retstart=%ld", (long) retstart);
    StringCat (ptr, buf);
  }
  if (retmax > 1) {
    sprintf (buf, "&retmax=%ld", (long) retmax);
    StringCat (ptr, buf);
  }

  MemFree (db);
  MemFree (qu);
  MemFree (sf);

  str = QUERY_UrlSynchronousQuery (eutils_host, 80, esearch_url, NULL, NULL, NULL, ptr);

  MemFree (ptr);

  return str;
}

NLM_EXTERN CharPtr DoElinkQuery (
  CharPtr dbfrom,
  CharPtr dbto,
  CharPtr linkname,
  CharPtr cmd,
  CharPtr ids,
  CharPtr querykey,
  CharPtr webenv
)

{
  CharPtr  cd, df, dt, id, lk, pfx, ptr, str;
  size_t   len;

  if (StringHasNoText (dbfrom) && StringHasNoText (dbto) && StringHasNoText (linkname)) return NULL;

  df = UrlEncodeString (dbfrom);
  dt = UrlEncodeString (dbto);
  lk = UrlEncodeString (linkname);
  cd = UrlEncodeString (cmd);
  id = UrlEncodeString (ids);

  len = StringLen (df) + StringLen (dt) + StringLen (lk) + StringLen (cd) + StringLen (id);
  if (len < 1) return NULL;

  ptr = (CharPtr) MemNew (len + 200);
  if (ptr == NULL) return NULL;

  pfx = NULL;
  if (StringDoesHaveText (dt)) {
    StringCat (ptr, "db=");
    StringCat (ptr, dt);
    pfx = "&";
  }
  if (StringDoesHaveText (df)) {
    StringCat (ptr, pfx);
    StringCat (ptr, "dbfrom=");
    StringCat (ptr, df);
    pfx = "&";
  }
  if (StringDoesHaveText (lk)) {
    StringCat (ptr, pfx);
    StringCat (ptr, "linkname=");
    StringCat (ptr, lk);
    pfx = "&";
  }
  if (StringDoesHaveText (cd)) {
    StringCat (ptr, pfx);
    StringCat (ptr, "cmd=");
    StringCat (ptr, cd);
    pfx = "&";
  }
  if (StringDoesHaveText (id)) {
    StringCat (ptr, "&id=");
    StringCat (ptr, id);
  } else if (StringDoesHaveText (querykey) && StringDoesHaveText (webenv)) {
    StringCat (ptr, "&query_key=");
    StringCat (ptr, querykey);
    StringCat (ptr, "&WebEnv=");
    StringCat (ptr, webenv);
  }

  MemFree (df);
  MemFree (dt);
  MemFree (lk);
  MemFree (cd);
  MemFree (id);

  str = QUERY_UrlSynchronousQuery (eutils_host, 80, elink_url, NULL, NULL, NULL, ptr);

  MemFree (ptr);

  return str;
}

NLM_EXTERN CharPtr DoEsummaryQuery (
  CharPtr database,
  CharPtr ids,
  CharPtr querykey,
  CharPtr webenv,
  Int4 retstart,
  Int4 retmax
)

{
  Char     buf [64];
  CharPtr  db, id, ptr, str;
  size_t   len;

  if (StringHasNoText (database)) return NULL;

  db = UrlEncodeString (database);
  id = UrlEncodeString (ids);

  len = StringLen (db) + StringLen (id) + StringLen (querykey) + StringLen (webenv);
  if (len < 1) return NULL;

  ptr = (CharPtr) MemNew (len + 200);
  if (ptr == NULL) return NULL;

  StringCpy (ptr, "db=");
  StringCat (ptr, db);
  if (StringDoesHaveText (id)) {
    StringCat (ptr, "&id=");
    StringCat (ptr, id);
  } else if (StringDoesHaveText (querykey) && StringDoesHaveText (webenv)) {
    StringCat (ptr, "&query_key=");
    StringCat (ptr, querykey);
    StringCat (ptr, "&WebEnv=");
    StringCat (ptr, webenv);
  }
  if (retstart > 1) {
    sprintf (buf, "&retstart=%ld", (long) retstart);
    StringCat (ptr, buf);
  }
  if (retmax > 1) {
    sprintf (buf, "&retmax=%ld", (long) retmax);
    StringCat (ptr, buf);
  }

  MemFree (db);
  MemFree (id);

  str = QUERY_UrlSynchronousQuery (eutils_host, 80, esummary_url, NULL, NULL, NULL, ptr);

  MemFree (ptr);

  return str;
}

static Int4 GetNumericValue (
  CharPtr str
)

{
  long int  val;

  if (StringHasNoText (str)) return 0;
  if (! StringIsAllDigits (str)) return 0;

  if (sscanf (str, "%ld", &val) != 1) return 0;

  return val;
}

static Boolean GetBooleanValue (
  CharPtr str
)

{
  Char  ch;

  if (StringHasNoText (str)) return FALSE;

  ch = *str;
  if (ch == 'Y' || ch == 'y') return TRUE;

  return FALSE;
}

static Entrez2DbInfoPtr GetDbInfo (
  CharPtr database
)

{
  Entrez2DbInfoPtr     e2db = NULL;
  Entrez2FieldInfoPtr  e2fp, lastfp = NULL;
  Entrez2LinkInfoPtr   e2lp, lastlp = NULL;
  XmlObjPtr            nxt, sub, tmp, xop;
  CharPtr              str;

  str = DoEinfoQuery (database);
  if (str == NULL) return NULL;

  xop = ParseXmlString (str);
  MemFree (str);
  if (xop == NULL) return NULL;

  e2db = Entrez2DbInfoNew ();
  if (e2db == NULL) return NULL;

  for (tmp = xop; tmp != NULL; tmp = nxt) {
    nxt = tmp->successor;
    if (XmlPathSuffixIs (tmp, "/DbInfo/DbName")) {
      e2db->db_name = StringSaveNoNull (tmp->contents);
    } else if (XmlPathSuffixIs (tmp, "/DbInfo/MenuName")) {
      e2db->db_menu = StringSaveNoNull (tmp->contents);
    } else if (XmlPathSuffixIs (tmp, "/DbInfo/Description")) {
      e2db->db_descr = StringSaveNoNull (tmp->contents);
    } else if (XmlPathSuffixIs (tmp, "/DbInfo/Count")) {
      e2db->doc_count = GetNumericValue (tmp->contents);
    } else if (XmlPathSuffixIs (tmp, "/FieldList/Field")) {
      e2fp = Entrez2FieldInfoNew ();
      if (e2fp != NULL) {
        for (sub = tmp->successor; sub != NULL && sub->level > tmp->level; sub = sub->successor) {
          if (XmlPathSuffixIs (sub, "/Field/Name")) {
            e2fp->field_name = StringSaveNoNull (sub->contents);
          } else if (XmlPathSuffixIs (sub, "/Field/FullName")) {
            e2fp->field_menu = StringSaveNoNull (sub->contents);
          } else if (XmlPathSuffixIs (sub, "/Field/Description")) {
            e2fp->field_descr = StringSaveNoNull (sub->contents);
          } else if (XmlPathSuffixIs (sub, "/Field/TermCount")) {
            e2fp->term_count = GetNumericValue (sub->contents);
          } else if (XmlPathSuffixIs (sub, "/Field/IsDate")) {
            e2fp->is_date = GetBooleanValue (sub->contents);
          } else if (XmlPathSuffixIs (sub, "/Field/IsNumerical")) {
            e2fp->is_numerical = GetBooleanValue (sub->contents);
          } else if (XmlPathSuffixIs (sub, "/Field/SingleToken")) {
            e2fp->single_token = GetBooleanValue (sub->contents);
          } else if (XmlPathSuffixIs (sub, "/Field/Hierarchy")) {
            e2fp->hierarchy_avail = GetBooleanValue (sub->contents);
          }
/* is_rangable and is_truncatable are not present in XML */
        }
        if (lastfp != NULL) {
          lastfp->next = e2fp;
        } else {
          e2db->fields = e2fp;
        }
        lastfp = e2fp;
        (e2db->field_count)++;
        nxt = sub;
      }
    } else if (XmlPathSuffixIs (tmp, "/LinkList/Link")) {
      e2lp = Entrez2LinkInfoNew ();
      if (e2lp != NULL) {
        e2lp->data_size = 4;
        for (sub = tmp->successor; sub != NULL && sub->level > tmp->level; sub = sub->successor) {
          if (XmlPathSuffixIs (sub, "/Link/Name")) {
            e2lp->link_name = StringSaveNoNull (sub->contents);
          } else if (XmlPathSuffixIs (sub, "/Link/Menu")) {
            e2lp->link_menu = StringSaveNoNull (sub->contents);
          } else if (XmlPathSuffixIs (sub, "/Link/Description")) {
            e2lp->link_descr = StringSaveNoNull (sub->contents);
          } else if (XmlPathSuffixIs (sub, "/Link/DbTo")) {
            e2lp->db_to = StringSaveNoNull (sub->contents);
          }
        }
        if (lastlp != NULL) {
          lastlp->next = e2lp;
        } else {
          e2db->links = e2lp;
        }
        lastlp = e2lp;
        (e2db->link_count)++;
        nxt = sub;
      }
    }
/* docsum_fields is not present in XML */
  }
  FreeXmlObject (xop);

  return e2db;
}

static Entrez2ReplyPtr GetE2info (
  Entrez2RequestPtr e2rq
)

{
  Entrez2DbInfoPtr  e2db, last;
  Entrez2ReplyPtr   e2ry;
  Entrez2InfoPtr    eip;
  ValNodePtr        head, tail, reply, vnp;
  XmlObjPtr         nxt, tmp, xop;
  CharPtr           str;

  if (e2rq == NULL) return NULL;

  vnp = e2rq->request;
  if (vnp == NULL) return NULL;
  if (vnp->choice != E2Request_get_info) return NULL;

  str = DoEinfoQuery (NULL);
  if (str == NULL) return NULL;

  xop = ParseXmlString (str);
  MemFree (str);
  if (xop == NULL) return NULL;

  eip = Entrez2InfoNew ();
  if (eip == NULL) return NULL;

  head = NULL;
  tail = NULL;
  for (tmp = xop; tmp != NULL; tmp = nxt) {
    nxt = tmp->successor;
    if (XmlPathSuffixIs (tmp, "/DbList/DbName")) {
      if (StringHasNoText (tmp->contents)) continue;
      ValNodeCopyStrEx (&head, &tail, 0, tmp->contents);
    }
  }
  FreeXmlObject (xop);
  if (head == NULL) {
    Entrez2InfoFree (eip);
    return NULL;
  }

  head = ValNodeSort (head, SortVnpByString);

  last = NULL;
  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    e2db = GetDbInfo (str);
    if (last != NULL) {
      last->next = e2db;
    } else {
      eip->db_info = e2db;
    }
    last = e2db;
    (eip->db_count)++;
  }

  ValNodeFreeData (head);

  e2ry = Entrez2ReplyNew ();
  if (e2ry == NULL) return NULL;
  reply = ValNodeNew (NULL);
  if (reply == NULL) return NULL;
  e2ry->reply = reply;
  reply->choice = E2Reply_get_info;
  reply->data.ptrvalue = (Pointer) eip;

  return e2ry;
}

static CharPtr ParseE2bool (
  ValNodePtr exp
)

{
  Entrez2BooleanTermPtr  e2bt;
  ValNodePtr             head = NULL, tail = NULL;
  CharPtr                str;

  if (exp == NULL) return NULL;

  while (exp != NULL) {
    switch (exp->choice) {
      case Entrez2BooleanElement_str :
        str = (CharPtr) exp->data.ptrvalue;
        if (StringHasNoText (str)) break;
        ValNodeCopyStrEx (&head, &tail, 0, str);
        break;
      case Entrez2BooleanElement_op :
        switch (exp->data.intvalue) {
          case 1 :
            ValNodeCopyStrEx (&head, &tail, 0, " AND ");
            break;
          case 2 :
            ValNodeCopyStrEx (&head, &tail, 0, " OR ");
            break;
          case 3 :
            ValNodeCopyStrEx (&head, &tail, 0, " NOT ");
            break;
          case 4 :
            ValNodeCopyStrEx (&head, &tail, 0, ":");
            break;
          case 5 :
            ValNodeCopyStrEx (&head, &tail, 0, "(");
            break;
          case 6 :
            ValNodeCopyStrEx (&head, &tail, 0, ")");
            break;
          default :
            break;
        }
        break;
      case Entrez2BooleanElement_term :
        e2bt = (Entrez2BooleanTermPtr) exp->data.ptrvalue;
        if (e2bt == NULL) break;
        if (StringHasNoText (e2bt->term)) break;
        ValNodeCopyStrExEx (&head, &tail, 0, e2bt->term, "\"", "\"");
        if (StringHasNoText (e2bt->field)) break;
        ValNodeCopyStrExEx (&head, &tail, 0, e2bt->field, "[", "]");
        break;
      case Entrez2BooleanElement_ids :
        break;
      case Entrez2BooleanElement_key :
        break;
      default :
        break;
    }
    exp = exp->next;
  }

  if (head == NULL) return NULL;

  str = ValNodeMergeStrs (head);
  ValNodeFreeData (head);

  return str;
}

static Entrez2ReplyPtr GetE2bool (
  Entrez2RequestPtr e2rq
)

{
  ByteStorePtr            bs;
  CharPtr                 cookie, count, database, key, query, str;
  Entrez2BooleanExpPtr    e2be;
  Entrez2BooleanReplyPtr  e2br;
  Entrez2EvalBooleanPtr   e2eb;
  Entrez2IdListPtr        e2id;
  Entrez2ReplyPtr         e2ry;
  ValNodePtr              head, tail, reply, vnp;
  Int4                    num, uid;
  XmlObjPtr               tmp, xop;

  if (e2rq == NULL) return NULL;

  vnp = e2rq->request;
  if (vnp == NULL) return NULL;
  if (vnp->choice != E2Request_eval_boolean) return NULL;

  e2eb = (Entrez2EvalBooleanPtr) vnp->data.ptrvalue;
  if (e2eb == NULL) return NULL;

  e2be = e2eb->query;
  if (e2be == NULL) return NULL;
  if (StringHasNoText (e2be->db)) return NULL;
  database = e2be->db;
  if (StringHasNoText (database)) return NULL;

  query = ParseE2bool (e2be->exp);
  if (query == NULL) return NULL;

  if (StringHasNoText (query)) {
    MemFree (query);
    return NULL;
  }

  str = DoEsearchQuery (database, query, NULL, /* e2rq->cookie */ NULL, 0, 100000, (Boolean) e2eb->return_UIDs);

  MemFree (query);
  if (str == NULL) return NULL;

  xop = ParseXmlString (str);
  MemFree (str);
  if (xop == NULL) return NULL;

  head = NULL;
  tail = NULL;
  count= NULL;
  key = NULL;
  cookie = NULL;
  for (tmp = xop; tmp != NULL; tmp = tmp->successor) {
    if (XmlPathSuffixIs (tmp, "/Id")) {
      ValNodeCopyStrEx (&head, &tail, 0, tmp->contents);
    } else if (XmlPathSuffixIs (tmp, "/eSearchResult/Count")) {
      count = StringSave (tmp->contents);
    } else if (XmlPathSuffixIs (tmp, "/eSearchResult/QueryKey")) {
      key = StringSave (tmp->contents);
    } else if (XmlPathSuffixIs (tmp, "/eSearchResult/WebEnv")) {
      cookie = StringSave (tmp->contents);
    }
  }

  FreeXmlObject (xop);

  num = GetNumericValue (count);
  MemFree (count);

  e2br = Entrez2BooleanReplyNew ();
  if (e2br == NULL) {
    ValNodeFreeData (head);
    return NULL;
  }

  e2br->count = num;
  if (num >= 0 && head != NULL) {
    e2id = Entrez2IdListNew ();
    if (e2id != NULL) {
      e2br->uids = e2id;
      e2id->db = StringSave (database);
      e2id->num = num;
      bs = BSNew (num * sizeof (Int4));
      if (bs != NULL) {
        e2id->uids = bs;
        for (vnp = head; vnp != NULL; vnp = vnp->next) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (StringHasNoText (str)) continue;
          uid = GetNumericValue (str);
          BSWrite (bs, (Pointer) &uid, sizeof (Int4));
        }
      }
    }
  }

  ValNodeFreeData (head);

  e2ry = Entrez2ReplyNew ();
  if (e2ry == NULL) return NULL;
  reply = ValNodeNew (NULL);
  if (reply == NULL) return NULL;
  e2ry->reply = reply;
  reply->choice = E2Reply_eval_boolean;
  reply->data.ptrvalue = (Pointer) e2br;
  e2ry-> key = key;
  e2ry->cookie = cookie;

  return e2ry;
}

static Entrez2ReplyPtr GetE2link (
  Entrez2RequestPtr e2rq
)

{
  ByteStorePtr        bs;
  Char                buf [32];
  Entrez2GetLinksPtr  e2gl;
  Entrez2IdListPtr    e2id;
  Entrez2LinkSetPtr   e2ls;
  Entrez2ReplyPtr     e2ry;
  ValNodePtr          head, tail, reply, vnp;
  CharPtr             ids, str, target;
  Int4                j, num, uid;
  XmlObjPtr           tmp, xop;

  if (e2rq == NULL) return NULL;

  vnp = e2rq->request;
  if (vnp == NULL) return NULL;
  if (vnp->choice != E2Request_get_links) return NULL;

  e2gl = (Entrez2GetLinksPtr) vnp->data.ptrvalue;
  if (e2gl == NULL) return NULL;

  e2id = e2gl->uids;
  if (e2id == NULL) return NULL;

  num = e2id->num;
  bs = e2id->uids;
  if (bs == NULL) return NULL;

  if (BSLen (bs) / sizeof (Int4) != num) {
    return NULL;
  }

  head = NULL;
  tail = NULL;
  BSSeek (bs, 0L, 0);

  for (j = 0; j < num; j++) {
    BSRead (bs, &uid, sizeof (Int4));
    sprintf (buf, "%ld", (long) uid);
    ValNodeCopyStrEx (&head, &tail, 0, buf);
  }

  if (head == NULL) return NULL;

  ids = ValNodeMergeStrsEx (head, ",");
  ValNodeFreeData (head);

  str = DoElinkQuery (e2id->db, NULL, e2gl->linktype, NULL, ids, NULL, NULL);

  MemFree (ids);
  if (str == NULL) return NULL;

  xop = ParseXmlString (str);
  MemFree (str);
  if (xop == NULL) return NULL;

  head = NULL;
  tail = NULL;
  target = NULL;
  for (tmp = xop; tmp != NULL; tmp = tmp->successor) {
    if (XmlPathSuffixIs (tmp, "/Link/Id")) {
      ValNodeCopyStrEx (&head, &tail, 0, tmp->contents);
    } else if (XmlPathSuffixIs (tmp, "/LinkSetDb/DbTo")) {
      target = StringSave (tmp->contents);
    }
  }

  FreeXmlObject (xop);

  if (head == NULL) return NULL;

  e2ls = Entrez2LinkSetNew ();
  if (e2ls == NULL) return NULL;

  num = ValNodeLen (head);
  if (num >= 0 && head != NULL) {
    e2id = Entrez2IdListNew ();
    if (e2id != NULL) {
      e2ls->ids = e2id;
      e2id->db = StringSave (target);
      e2id->num = num;
      bs = BSNew (num * sizeof (Int4));
      if (bs != NULL) {
        e2id->uids = bs;
        for (vnp = head; vnp != NULL; vnp = vnp->next) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (StringHasNoText (str)) continue;
          uid = GetNumericValue (str);
          BSWrite (bs, (Pointer) &uid, sizeof (Int4));
        }
      }
    }
  }

  MemFree (target);
  ValNodeFreeData (head);

  e2ry = Entrez2ReplyNew ();
  if (e2ry == NULL) return NULL;
  reply = ValNodeNew (NULL);
  if (reply == NULL) return NULL;
  e2ry->reply = reply;
  reply->choice = E2Reply_get_links;
  reply->data.ptrvalue = (Pointer) e2ls;

  return e2ry;
}

static Entrez2ReplyPtr GetE2summary (
  Entrez2RequestPtr e2rq
)

{
  XmlObjPtr             attr, lst, nxt, scc, sub, tmp, xop;
  ByteStorePtr          bs;
  Char                  buf [32];
  Entrez2DocsumDataPtr  e2dd, laste2dd;
  Entrez2DocsumListPtr  e2dl;
  Entrez2DocsumPtr      e2ds, laste2ds;
  Entrez2IdListPtr      e2id;
  Entrez2ReplyPtr       e2ry;
  ValNodePtr            head, tail, reply, vnp;
  CharPtr               ids, name, sep, str, type;
  Int4                  j, num, uid;

  if (e2rq == NULL) return NULL;

  vnp = e2rq->request;
  if (vnp == NULL) return NULL;
  if (vnp->choice != E2Request_get_docsum) return NULL;

  e2id = (Entrez2IdListPtr) vnp->data.ptrvalue;
  if (e2id == NULL) return NULL;

  num = e2id->num;
  bs = e2id->uids;
  if (bs == NULL) return NULL;

  if (BSLen (bs) / sizeof (Int4) != num) {
    return NULL;
  }

  head = NULL;
  tail = NULL;
  BSSeek (bs, 0L, 0);

  for (j = 0; j < num; j++) {
    BSRead (bs, &uid, sizeof (Int4));
    sprintf (buf, "%ld", (long) uid);
    ValNodeCopyStrEx (&head, &tail, 0, buf);
  }

  if (head == NULL) return NULL;

  ids = ValNodeMergeStrsEx (head, ",");
  ValNodeFreeData (head);

  str = DoEsummaryQuery (e2id->db, ids, NULL, NULL, 0, 100000);

  MemFree (ids);
  if (str == NULL) return NULL;

  xop = ParseXmlString (str);
  MemFree (str);
  if (xop == NULL) return NULL;

  e2dl = Entrez2DocsumListNew ();
  if (e2dl == NULL) return NULL;

  laste2ds = NULL;
  for (tmp = xop; tmp != NULL; tmp = nxt) {
    nxt = tmp->successor;
    if (XmlPathSuffixIs (tmp, "/DocSum")) {
      e2ds = Entrez2DocsumNew ();
      if (e2ds == NULL) continue;
      laste2dd = NULL;
      for (sub = tmp->successor; sub != NULL && sub->level > tmp->level; sub = scc) {
        scc = sub->successor;
        if (XmlPathSuffixIs (sub, "/DocSum/Id")) {
          e2ds->uid = GetNumericValue (sub->contents);
        } else if (XmlPathSuffixIs (sub, "/DocSum/Item")) {
          name = NULL;
          type = NULL;
          e2dd = NULL;
          for (attr = sub->attributes; attr != NULL; attr = attr->next) {
            if (StringICmp (attr->name, "Name") == 0) {
              name = attr->contents;
            } else if (StringICmp (attr->name, "Type") == 0) {
              type = attr->contents;
            }
          }
          if (StringHasNoText (name) || StringHasNoText (type)) {
            /* skip */
          } else if (StringICmp (type, "List") == 0) {
            head = NULL;
            tail = NULL;
            for (lst = sub->successor; lst != NULL && lst->level > sub->level; lst = lst->successor) {
              ValNodeCopyStrEx (&head, &tail, 0, lst->contents);
            }
            scc = lst;
            if (head != NULL) {
              sep = ",";
              if (StringICmp (name, "AuthorList") == 0) {
                sep = ", ";
                name = "Authors";
              }
              str = ValNodeMergeStrsEx (head, sep);
              e2dd = Entrez2DocsumDataNew ();
              if (e2dd != NULL) {
                e2dd->field_name = StringSave (name);
                e2dd->field_value = StringSave (str);
              }
              MemFree (str);
            }
            ValNodeFreeData (head);
          } else if (StringHasNoText (sub->contents)) {
            /* skip */
          } else {
            if (StringICmp (name, "HasAbstract") == 0) {
              name = "Attributes";
            }
            e2dd = Entrez2DocsumDataNew ();
            if (e2dd != NULL) {
              e2dd->field_name = StringSave (name);
              e2dd->field_value = StringSave (sub->contents);
            }
          }
          if (e2dd != NULL) {
            if (laste2dd != NULL) {
              laste2dd->next = e2dd;
            } else {
              e2ds->docsum_data = e2dd;
            }
          laste2dd = e2dd;
          }
        }
      }
      if (laste2ds != NULL) {
        laste2ds->next = e2ds;
      } else {
        e2dl->list = e2ds;
      }
      laste2ds = e2ds;
      (e2dl->count)++;
      nxt = sub;
    }
  }

  FreeXmlObject (xop);

  e2ry = Entrez2ReplyNew ();
  if (e2ry == NULL) return NULL;
  reply = ValNodeNew (NULL);
  if (reply == NULL) return NULL;
  e2ry->reply = reply;
  reply->choice = E2Reply_get_docsum;
  reply->data.ptrvalue = (Pointer) e2dl;

  return e2ry;
}

#define ENTREZ_TOOL_PROPERTY "Entrez2Tool"
#define ENTREZ_TOOL_VERSION 1

/* utility functions */

NLM_EXTERN void EntrezSetProgramName (
  const char* progname
)

{
  MemFree (GetAppProperty (ENTREZ_TOOL_PROPERTY));
  SetAppProperty (ENTREZ_TOOL_PROPERTY, (StringHasNoText (progname)
                                         ? NULL
                                         : StringSave (progname)));
}

static const char* EntrezGetProgramName (
  void
)

{
  char         path [PATH_MAX];
  const char*  ptr;

  ptr = (const char*) GetAppProperty (ENTREZ_TOOL_PROPERTY);
  if (StringHasNoText (ptr)) {
    Nlm_ProgramPath (path, sizeof (path));
    ptr = StringRChr (path, DIRDELIMCHR);
    if (ptr != NULL) {
      ptr++;
      EntrezSetProgramName (ptr);
      ptr = (const char*) GetAppProperty (ENTREZ_TOOL_PROPERTY);
    }
  }
  return ptr;
}

/* override service name */
static const char*  e2_service = NULL;

/* use EntrezTest to override default Entrez ncbi named service */

NLM_EXTERN void EntrezSetService (
  const char* service
)

{
  MemFree ((void*) e2_service);
  e2_service = StringSaveNoNull (service);
}

/* low-level connection functions */

static const char* GetDbFromE2Request (Entrez2RequestPtr e2rq)

{
  Entrez2BooleanExpPtr   e2be;
  Entrez2EvalBooleanPtr  e2eb;
  Entrez2GetLinksPtr     e2gl;
  Entrez2HierQueryPtr    e2hq;
  Entrez2IdPtr           e2id;
  Entrez2IdListPtr       e2il;
  Entrez2TermPosPtr      e2tp;
  Entrez2TermQueryPtr    e2tq;
  ValNodePtr             vnp;

  if (e2rq == NULL) return NULL;

  vnp = e2rq->request;
  if (vnp == NULL) return NULL;

  switch (vnp->choice) {
    case E2Request_get_info :
      break;
    case E2Request_eval_boolean :
      e2eb = (Entrez2EvalBooleanPtr) vnp->data.ptrvalue;
      if (e2eb == NULL) return NULL;
      e2be = e2eb->query;
      if (e2be == NULL) return NULL;
      if (StringDoesHaveText (e2be->db)) return e2be->db;
      break;
    case E2Request_get_docsum :
      e2il = (Entrez2IdListPtr) vnp->data.ptrvalue;
      if (e2il == NULL) return NULL;
      if (StringDoesHaveText (e2il->db)) return e2il->db;
      break;
    case E2Request_get_term_pos :
      e2tq = (Entrez2TermQueryPtr) vnp->data.ptrvalue;
      if (e2tq == NULL) return NULL;
      if (StringDoesHaveText (e2tq->db)) return e2tq->db;
      break;
    case E2Request_get_term_list :
      e2tp = (Entrez2TermPosPtr) vnp->data.ptrvalue;
      if (e2tp == NULL) return NULL;
      if (StringDoesHaveText (e2tp->db)) return e2tp->db;
      break;
    case E2Request_get_term_hierarchy :
      e2hq = (Entrez2HierQueryPtr) vnp->data.ptrvalue;
      if (e2hq == NULL) return NULL;
      if (StringDoesHaveText (e2hq->db)) return e2hq->db;
      break;
    case E2Request_get_links :
      e2gl = (Entrez2GetLinksPtr) vnp->data.ptrvalue;
      if (e2gl == NULL) return NULL;
      e2il = (Entrez2IdListPtr) e2gl->uids;
      if (e2il == NULL) return NULL;
      if (StringDoesHaveText (e2il->db)) return e2il->db;
      break;
    case E2Request_get_linked :
      e2gl = (Entrez2GetLinksPtr) vnp->data.ptrvalue;
      if (e2gl == NULL) return NULL;
      e2il = (Entrez2IdListPtr) e2gl->uids;
      if (e2il == NULL) return NULL;
      if (StringDoesHaveText (e2il->db)) return e2il->db;
      break;
    case E2Request_get_link_counts :
      e2id = (Entrez2IdPtr) vnp->data.ptrvalue;
      if (e2id == NULL) return NULL;
      if (StringDoesHaveText (e2id->db)) return e2id->db;
      break;
    default :
      break;
  }

  return NULL;
}

NLM_EXTERN CONN EntrezOpenConnection (
  Entrez2RequestPtr e2rq
)

{
  char         arg [128];
  const char*  db;

  db = GetDbFromE2Request (e2rq);
  if (StringDoesHaveText (db) && StringLen (db) < 100) {
    StrCpy (arg,    "DB=");
    StrCpy (arg + 3, db);
  } else
    *arg = '\0';

  return QUERY_OpenServiceQueryEx
    (StringHasNoText (e2_service) ? "Entrez2" : e2_service, NULL, 30, arg);
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
  time_t           max = 0;
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
  while ((status = CONN_Wait (conn, eIO_Read, &timeout)) == eIO_Timeout && max < 300) {
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

static TNlmTls e2cookie_tls = NULL;

/* high-level connection functions */

static Entrez2ReplyPtr EntrezViaEUtils (
  Entrez2RequestPtr e2rq
)

{
  ValNodePtr  vnp;

  if (e2rq == NULL) return NULL;

  vnp = e2rq->request;
  if (vnp == NULL) return NULL;

  switch (vnp->choice) {
    case E2Request_get_info :
      return GetE2info (e2rq);
    case E2Request_eval_boolean :
      return GetE2bool (e2rq);
    case E2Request_get_docsum :
      return GetE2summary (e2rq);
    case E2Request_get_links :
      return GetE2link (e2rq);
    default :
      return NULL;
  }

  return NULL;
}

NLM_EXTERN Entrez2ReplyPtr EntrezSynchronousQuery (
  Entrez2RequestPtr e2rq
)

{
  AsnIoConnPtr     aicp;
  CONN             conn;
  char*            e2cookie = NULL;
  Entrez2ReplyPtr  e2ry;
  char*            tempcookie = NULL;
#ifdef OS_UNIX
  Boolean          logtimes;
  clock_t          starttime;
  clock_t          stoptime;
  struct tms       timebuf;
#endif

  if (e2rq == NULL) return NULL;

/* conditionally try EUtils replacement first */
#ifdef OS_UNIX
  if (getenv ("NCBI_EUTILS_REPLACES_ENTREZ2") != NULL) {
    e2ry = EntrezViaEUtils (e2rq);
    if (e2ry != NULL) return e2ry;
  }
#endif

#ifdef OS_UNIX
  logtimes = (Boolean) ((getenv ("NCBI_LOG_SYNC_QUERY_TIMES")) != NULL);
#endif

  conn = EntrezOpenConnection (e2rq);

  if (conn == NULL) return NULL;

  aicp = QUERY_AsnIoConnOpen ("wb", conn);

  tempcookie = e2rq->cookie;
  if (NlmTlsGetValue (e2cookie_tls, (VoidPtr PNTR) &e2cookie)) {
    if (e2rq->cookie == NULL && e2cookie != NULL) {
      e2rq->cookie = e2cookie;
    }
  }

  Entrez2RequestAsnWrite (e2rq, aicp->aip, NULL);

  e2rq->cookie = tempcookie;

  AsnIoFlush (aicp->aip);
  QUERY_AsnIoConnClose (aicp);

  QUERY_SendQuery (conn);

#ifdef OS_UNIX
  if (logtimes) {
    starttime = times (&timebuf);
  }
#endif

  e2ry = EntrezWaitForReply (conn);

#ifdef OS_UNIX
  if (logtimes) {
    stoptime = times (&timebuf);
    printf ("EntrezWaitForReply %ld\n", (long) (stoptime - starttime));
  }
#endif

  if (e2ry != NULL && e2ry->cookie != NULL) {
    if (NlmTlsGetValue (e2cookie_tls, (VoidPtr PNTR) &e2cookie)) {
      e2cookie = MemFree (e2cookie);
      e2cookie = StringSave (e2ry->cookie);
      NlmTlsSetValue (&e2cookie_tls, (VoidPtr PNTR) e2cookie, NULL);
    }
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
  char*         e2cookie = NULL;
  char*         tempcookie = NULL;

  if (e2rq == NULL) return FALSE;

  conn = EntrezOpenConnection (e2rq);

  if (conn == NULL) return FALSE;

  aicp = QUERY_AsnIoConnOpen ("wb", conn);

  tempcookie = e2rq->cookie;
  if (NlmTlsGetValue (e2cookie_tls, (VoidPtr PNTR) &e2cookie)) {
    if (e2rq->cookie == NULL && e2cookie != NULL) {
      e2rq->cookie = e2cookie;
    }
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
  char*            e2cookie = NULL;
  Entrez2ReplyPtr  e2ry = NULL;

  if (conn != NULL && status == eIO_Success) {
    aicp = QUERY_AsnIoConnOpen ("rb", conn);
    e2ry = Entrez2ReplyAsnRead (aicp->aip, NULL);
    QUERY_AsnIoConnClose (aicp);
  }

  if (e2ry != NULL && e2ry->cookie != NULL) {
    if (NlmTlsGetValue (e2cookie_tls, (VoidPtr PNTR) &e2cookie)) {
      e2cookie = MemFree (e2cookie);
      e2cookie = StringSave (e2ry->cookie);
      NlmTlsSetValue (&e2cookie_tls, (VoidPtr PNTR) e2cookie, NULL);
    }
  }

  return e2ry;
}

/* request creation functions */

static Entrez2RequestPtr CreateRequest (
  Uint1 choice, Pointer data
)

{
  char*              e2cookie = NULL;
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

  if (NlmTlsGetValue (e2cookie_tls, (VoidPtr PNTR) &e2cookie)) {
    e2rq->cookie = StringSaveNoNull (e2cookie);
  }

  return e2rq;
}

/* history needs to be used for Boolean ids and key queries */

NLM_EXTERN void EntrezSetUseHistoryFlag (
  Entrez2RequestPtr e2rq
)

{
  if (e2rq == NULL) return;
  e2rq->use_history = TRUE;
}

NLM_EXTERN Entrez2IdListPtr EntrezCreateEntrezIdList (
  const char* db,
  Int4 uid,
  Int4 num,
  const Int4 uids[],
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
  const char* type_date,
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
  const char* db,
  const char* query_string,
  Int4 begin_date,
  Int4 end_date,
  const char* type_date,
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
    EntrezAddToBooleanRequest (e2rq, query_string, 0, NULL, NULL, NULL,
                               0, 0, NULL, NULL, TRUE, TRUE);
  }

  return e2rq;
}

NLM_EXTERN void EntrezAddToBooleanRequest (
  Entrez2RequestPtr e2rq,
  const char* query_string,
  Int4 op,
  const char* field,
  const char* term,
  const char* key,
  Int4 uid,
  Int4 num,
  const Int4 uids[],
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

  } else if (! StringHasNoText (key)) {
    ValNodeCopyStr (&(e2be->exp), Entrez2BooleanElement_key, key);

  } else {

    e2il = EntrezCreateEntrezIdList (e2be->db, uid, num, uids, bs);
    if (e2il == NULL) return;

    ValNodeAddPointer (&(e2be->exp), Entrez2BooleanElement_ids, (Pointer) e2il);
  }
}

NLM_EXTERN Entrez2RequestPtr EntrezCreateDocSumRequest (
  const char* db,
  Int4 uid,
  Int4 num,
  const Int4 uids[],
  ByteStorePtr bs
)

{
  Entrez2IdListPtr  e2il;

  e2il = EntrezCreateEntrezIdList (db, uid, num, uids, bs);
  if (e2il == NULL) return NULL;

  return CreateRequest (E2Request_get_docsum, (Pointer) e2il);
}

NLM_EXTERN Entrez2RequestPtr EntrezCreateGetTermPositionRequest (
  const char* db,
  const char* field,
  const char* term
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
  const char* db,
  const char* field,
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
  const char* db,
  const char* field,
  const char* term,
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
  const char* db,
  Int4 uid,
  Int4 num,
  const Int4 uids[],
  ByteStorePtr bs,
  const char* linktype,
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
  const char* db,
  Int4 uid,
  Int4 num,
  const Int4 uids[],
  ByteStorePtr bs,
  const char* linktype,
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
  const char* db,
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

NLM_EXTERN char* EntrezExtractErrorReply (
  Entrez2ReplyPtr e2ry
)

{
  return (char*) GeneralEntrezExtractReply (e2ry, E2Reply_error, NULL);
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
  const char* db,
  const char* seq_id_string
)

{
  char                    ch;
  Entrez2BooleanReplyPtr  e2br;
  Entrez2IdListPtr        e2id;
  Entrez2RequestPtr       e2rq;
  Entrez2ReplyPtr         e2ry;
  char*                   ptr;
  char                    str [61];
  Uint4                   uid = 0;

  if (StringHasNoText (db) || StringHasNoText (seq_id_string)) return 0;

  StringNCpy_0 (str, seq_id_string, sizeof (str) - 1);
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

/* result validation function */

static int LIBCALLBACK SortVnpByStr (VoidPtr ptr1, VoidPtr ptr2)

{
  const char*  str1;
  const char*  str2;
  ValNodePtr   vnp1;
  ValNodePtr   vnp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      str1 = (const char*) vnp1->data.ptrvalue;
      str2 = (const char*) vnp2->data.ptrvalue;
      if (str1 != NULL && str2 != NULL) {
        return StringCmp (str1, str2);
      }
    }
  }
  return 0;
}

static ValNodePtr SplitAtSpaces (CharPtr str)

{
  Char        ch;
  CharPtr     ptr;
  CharPtr     tmp;
  ValNodePtr  head = NULL;
  ValNodePtr  tail = NULL;

  if (StringHasNoText (str)) return NULL;

  tmp = StringSave (str);
  if (tmp == NULL) return NULL;

  str = tmp;
  while (str != NULL) {
    ptr = str;
    ch = *ptr;
    while (ch == ' ') {
      ptr++;
      ch = *ptr;
    }
    while (ch != '\0' && ch != ' ' && ch != '/' && ch != '-') {
      ptr++;
      ch = *ptr;
    }
    if (ch != '\0') {
      *ptr = '\0';
      ptr++;
    } else {
      ptr = NULL;
    }
    if (StringDoesHaveText (str)) {
      TrimSpacesAroundString (str);
      ValNodeCopyStrEx (&head, &tail, 0, str);
    }
    str = ptr;
  }

  MemFree (tmp);

  return head;
}

static Boolean LowerCaseWords (CharPtr str)

{
  Char        ch;
  ValNodePtr  head, vnp;
  Boolean     rsult = FALSE;

  if (StringHasNoText (str)) return FALSE;
  head = SplitAtSpaces (str);
  if (head == NULL) return FALSE;

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    if (StringCmp (str, "eISSN") == 0) continue;
    if (StringCmp (str, "pISSN") == 0) continue;
    if (StringCmp (str, "mRNA") == 0) continue;
    if (vnp != head && vnp->next != NULL) {
      if (StringCmp (str, "a") == 0) continue;
      if (StringCmp (str, "as") == 0) continue;
      if (StringCmp (str, "by") == 0) continue;
      if (StringCmp (str, "for") == 0) continue;
      if (StringCmp (str, "in") == 0) continue;
      if (StringCmp (str, "of") == 0) continue;
    }
    ch = *str;
    if (IS_ALPHA (ch)) {
      if (IS_LOWER (ch)) {
        rsult = TRUE;
      }
    }
  }

  ValNodeFreeData (head);
  return rsult;
}

static Boolean MixedCaseWords (CharPtr str)

{
  Char        ch;
  ValNodePtr  head, vnp;
  Boolean     have_seen_lower;
  Boolean     rsult = FALSE;

  if (StringHasNoText (str)) return FALSE;
  head = SplitAtSpaces (str);
  if (head == NULL) return FALSE;

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    if (StringCmp (str, "eISSN") == 0) continue;
    if (StringCmp (str, "pISSN") == 0) continue;
    if (StringCmp (str, "mRNA") == 0) continue;
    if (StringCmp (str, "PubMed") == 0) continue;
    if (StringCmp (str, "MeSH") == 0) continue;
    if (StringCmp (str, "LocusLink") == 0) continue;
    if (StringCmp (str, "UniGene") == 0) continue;
    if (StringCmp (str, "UniSTS") == 0) continue;
    have_seen_lower = FALSE;
    ch = *str;
    while (ch != '\0') {
      if (IS_ALPHA (ch)) {
        if (IS_UPPER (ch)) {
          if (have_seen_lower) {
            rsult = TRUE;
          }
        } else if (IS_LOWER (ch)) {
          have_seen_lower = TRUE;
        }
      }
      str++;
      ch = *str;
    }
  }

  ValNodeFreeData (head);
  return rsult;
}

NLM_EXTERN Boolean ValidateEntrez2InfoPtrExExEx (
  Entrez2InfoPtr e2ip,
  ValNodePtr PNTR head,
  Boolean checkMenuNameVariants,
  Boolean checkMenuNameFormat,
  Boolean fullMenuVariantReport
)

{
  Char                       buf [512];
  Char                       ch;
  CharPtr                    db;
  Int2                       dbcount;
  CharPtr                    dbnames [256];
  CharPtr                    dsf;
  Int2                       dsfcount;
  Entrez2DbInfoPtr           e2db;
  Entrez2DocsumFieldInfoPtr  e2dsp;
  Entrez2FieldInfoPtr        e2fip;
  Entrez2LinkInfoPtr         e2lip;
  CharPtr                    fld;
  Int2                       fldcount;
  Boolean                    hasLowCase;
  Int2                       i;
  CharPtr                    last;
  ValNodePtr                 lastvnp;
  size_t                     len1;
  size_t                     len2;
  CharPtr                    lnk;
  Int2                       lnkcount;
  ValNodePtr                 menuhead = NULL;
  Boolean                    notAlphNum;
  Boolean                    rsult = TRUE;
  CharPtr                    str;
  Char                       tmpdb [32];
  Char                       tmpdsf [32];
  Char                       tmpfld [32];
  Char                       tmplnk [32];
  ValNodePtr                 vnp;

  if (head != NULL) {
    *head = NULL;
  }
  if (e2ip == NULL) return FALSE;

  if (e2ip->db_count < 1 || e2ip->db_info == NULL) {
    sprintf (buf, "Entrez2 has no databases");
    ValNodeCopyStr (head, 0, buf);
    return FALSE;
  }

  for (i = 0; i < sizeof (dbnames) / sizeof (CharPtr); i++) {
    dbnames [i] = "?";
  }
  i = 0;
  for (e2db = e2ip->db_info; e2db != NULL; e2db = e2db->next) {
    i++;
    if (! StringHasNoText (e2db->db_name)) {
      dbnames [i] = e2db->db_name;
    } else if (! StringHasNoText (e2db->db_menu)) {
      dbnames [i] = e2db->db_menu;
    }
  }

  dbcount = 0;
  for (e2db = e2ip->db_info; e2db != NULL; e2db = e2db->next) {
    dbcount++;

    db = e2db->db_name;
    if (StringICmp (db, "gtr") == 0) continue;
    if (StringICmp (db, "genomeprj") == 0) continue;

    if (StringHasNoText (db)) {
      rsult = FALSE;
      if (StringHasNoText (e2db->db_menu)) {
        sprintf (tmpdb, "%d", (int) dbcount);
        db = tmpdb;
        sprintf (buf, "Database %d has no name", (int) dbcount);
        ValNodeCopyStr (head, 0, buf);
      } else {
        db = e2db->db_menu;
        sprintf (buf, "Database %s (%d) has no name", db, (int) dbcount);
        ValNodeCopyStr (head, 0, buf);
      }
    }

    if (StringHasNoText (e2db->db_menu)) {
      sprintf (buf, "Database %s has no menu name", db);
      ValNodeCopyStr (head, 0, buf);
      rsult = FALSE;
    }
    if (StringHasNoText (e2db->db_descr)) {
      sprintf (buf, "Database %s has no description", db);
      ValNodeCopyStr (head, 0, buf);
      rsult = FALSE;
    }

    if (e2db->doc_count < 1) {
      if (StringICmp (db, "Nucleotide") == 0) {
        /* now a virtual database consolidating NucCore, NucEst, NucGss */
      } else {
        sprintf (buf, "Database %s has no documents", db);
        ValNodeCopyStr (head, 0, buf);
        rsult = FALSE;
      }
    }
    if (e2db->field_count < 1 || e2db->fields == NULL) {
      sprintf (buf, "Database %s has no fields", db);
      ValNodeCopyStr (head, 0, buf);
      rsult = FALSE;
    }
    if (e2db->link_count < 1 || e2db->links == NULL) {
      if (StringICmp (db, "books") != 0 &&
          StringICmp (db, "gap") != 0 &&
          StringICmp (db, "gensat") != 0 &&
          StringICmp (db, "mesh") != 0 &&
          StringICmp (db, "ncbisearch") != 0 &&
          StringICmp (db, "nlmcatalog") != 0 &&
          StringICmp (db, "nucleotide") != 0 &&
          StringICmp (db, "nuccore") != 0 &&
          StringICmp (db, "nucgss") != 0 &&
          StringICmp (db, "nucest") != 0 &&
          StringICmp (db, "seqannot") != 0 &&
          StringICmp (db, "toolkit") != 0 &&
          StringICmp (db, "blastdbinfo") != 0 &&
          StringICmp (db, "virus") != 0 &&
          StringICmp (db, "toolkitall") != 0 &&
          StringICmp (db, "gencoll") != 0 &&
          StringICmp (db, "images") != 0 &&
          StringICmp (db, "geo") != 0 &&
          StringICmp (db, "journals") != 0 &&
          StringICmp (db, "genomeprj") != 0 &&
          StringICmp (db, "gtr") != 0 &&
          StringICmp (db, "gcassembly") != 0 &&
          StringICmp (db, "clone") != 0 &&
          StringICmp (db, "pubmedhealth") != 0 &&
          StringICmp (db, "assembly") != 0 &&
          StringICmp (db, "gapplus") != 0) {
        sprintf (buf, "Database %s has no links", db);
        ValNodeCopyStr (head, 0, buf);
        rsult = FALSE;
      }
    }
    if (e2db->docsum_field_count < 1 || e2db->docsum_fields == NULL) {
      sprintf (buf, "Database %s has no docsum fields", db);
      ValNodeCopyStr (head, 0, buf);
      rsult = FALSE;
    }

    fldcount = 0;
    for (e2fip = e2db->fields; e2fip != NULL; e2fip = e2fip->next) {
      fldcount++;

      fld = e2fip->field_name;
      if (StringHasNoText (fld)) {
        rsult = FALSE;
        if (StringHasNoText (e2fip->field_menu)) {
          sprintf (tmpfld, "%d", (int) dbcount);
          fld = tmpfld;
          sprintf (buf, "Database %s field %d has no name", db, (int) fldcount);
          ValNodeCopyStr (head, 0, buf);
        } else {
          fld = e2fip->field_menu;
          sprintf (buf, "Database %s field %s (%d) has no name", db, fld, (int) fldcount);
          ValNodeCopyStr (head, 0, buf);
        }
      } else if (StringCmp (fld, "SLEN") == 0 ||
          StringCmp (fld, "MLWT") == 0 ||
          StringCmp (fld, "PMID") == 0 ||
          StringCmp (fld, "LLID") == 0 ||
          StringCmp (fld, "UID") == 0) {
        if (! e2fip->is_numerical) {
          sprintf (buf, "Database %s field %s does not have is_numerical set", db, fld);
          ValNodeCopyStr (head, 0, buf);
          rsult = FALSE;
        }
      } else if (StringCmp (fld, "TEXT") == 0) {
        sprintf (buf, "Database %s field %s should be WORD", db, fld);
        ValNodeCopyStr (head, 0, buf);
        rsult = FALSE;
      } else if (StringCmp (fld, "ORGN") == 0) {
        if (e2fip->term_count == 0) {
          if (StringICmp (db, "Nucleotide") == 0) {
            /* now a virtual database consolidating NucCore, NucEst, NucGss */
          } else {
            sprintf (buf, "Database %s field %s term count is 0", db, fld);
            ValNodeCopyStr (head, 0, buf);
            rsult = FALSE;
           }
       }
      }
      if (StringLen (fld) > 4) {
        sprintf (buf, "Database %s field %s name is > 4 characters long", db, fld);
        ValNodeCopyStr (head, 0, buf);
        rsult = FALSE;
      }

      hasLowCase = FALSE;
      notAlphNum = FALSE;
      str = fld;
      ch = *str;
      while (ch != '\0') {
        if (IS_LOWER (ch)) {
          hasLowCase = TRUE;
        } else if (! (IS_ALPHANUM (ch))) {
          notAlphNum = TRUE;
        }
        str++;
        ch = *str;
      }
      if (hasLowCase) {
        sprintf (buf, "Database %s field %s has lower case letters", db, fld);
        ValNodeCopyStr (head, 0, buf);
        rsult = FALSE;
      }
      if (notAlphNum) {
        sprintf (buf, "Database %s field %s has non-alphanumeric characters", db, fld);
        ValNodeCopyStr (head, 0, buf);
        rsult = FALSE;
      }

      if (StringHasNoText (e2fip->field_menu)) {
        sprintf (buf, "Database %s field %s has no menu name", db, fld);
        ValNodeCopyStr (head, 0, buf);
        rsult = FALSE;
      } else {
        ValNodeCopyStr (&menuhead, (Int2) dbcount, e2fip->field_menu);
        if (StringStr (e2fip->field_menu, "Date") != NULL) {
          if (! e2fip->is_date) {
            sprintf (buf, "Database %s field %s does not have is_date set", db, fld);
            ValNodeCopyStr (head, 0, buf);
            rsult = FALSE;
          }
        } else if (StringICmp (e2fip->field_menu, "Mesh") == 0) {
          sprintf (buf, "Database %s field-menu %s should be MeSH Terms", db, e2fip->field_menu);
          ValNodeCopyStr (head, 0, buf);
          rsult = FALSE;
        }
      }
      if (StringHasNoText (e2fip->field_descr)) {
        sprintf (buf, "Database %s field %s has no description", db, fld);
        ValNodeCopyStr (head, 0, buf);
        rsult = FALSE;
      }
    }
    if (e2db->field_count != fldcount) {
      sprintf (buf, "Database %s field count %ld does not match fldcount %d", db, (long) e2db->field_count, (int) fldcount);
      ValNodeCopyStr (head, 0, buf);
      rsult = FALSE;
    }

    lnkcount = 0;
    for (e2lip = e2db->links; e2lip != NULL; e2lip = e2lip->next) {
      lnkcount++;

      lnk = e2lip->link_name;
      if (StringHasNoText (lnk)) {
        rsult = FALSE;
        if (StringHasNoText (e2lip->link_menu)) {
          sprintf (tmplnk, "%d", (int) lnkcount);
          lnk = tmplnk;
          sprintf (buf, "Database %s link %d has no name", db, (int) lnkcount);
          ValNodeCopyStr (head, 0, buf);
        } else {
          lnk = e2lip->link_menu;
          sprintf (buf, "Database %s link %s (%d) has no name", db, lnk, (int) lnkcount);
          ValNodeCopyStr (head, 0, buf);
        }
      }

      /*
      if (StringHasNoText (e2lip->link_menu)) {
        sprintf (buf, "Database %s link %s has no menu name", db, lnk);
        ValNodeCopyStr (head, 0, buf);
        rsult = FALSE;
      }
      */
      if (StringHasNoText (e2lip->link_descr)) {
        if (StringICmp (db, "nucest") == 0 && StringICmp (lnk, "nucest_gene_clust") == 0) {
        } else if (StringICmp (db, "gene") == 0 && StringICmp (lnk, "gene_nucest_clust") == 0) {
        } else {
          sprintf (buf, "Database %s link %s has no description", db, lnk);
          ValNodeCopyStr (head, 0, buf);
          rsult = FALSE;
        }
      }
      if (StringHasNoText (e2lip->db_to)) {
        sprintf (buf, "Database %s link %s has no target database", db, lnk);
        ValNodeCopyStr (head, 0, buf);
        rsult = FALSE;
      }
    }
    if (e2db->link_count != lnkcount) {
      sprintf (buf, "Database %s link count %ld does not match lnkcount %d", db, (long) e2db->link_count, (int) lnkcount);
      ValNodeCopyStr (head, 0, buf);
      rsult = FALSE;
    }

    dsfcount = 0;
    for (e2dsp = e2db->docsum_fields; e2dsp != NULL; e2dsp = e2dsp->next) {
      dsfcount++;

      dsf = e2dsp->field_name;
      if (StringHasNoText (dsf)) {
        rsult = FALSE;
        if (StringHasNoText (e2dsp->field_description)) {
          sprintf (tmpdsf, "%d", (int) dsfcount);
          dsf = tmpdsf;
          sprintf (buf, "Database %s link %d has no name", db, (int) dsfcount);
          ValNodeCopyStr (head, 0, buf);
        } else {
          dsf = e2dsp->field_description;
          sprintf (buf, "Database %s link %s (%d) has no name", db, dsf, (int) dsfcount);
          ValNodeCopyStr (head, 0, buf);
        }
      }

      if (StringHasNoText (e2dsp->field_description)) {
        sprintf (buf, "Database %s docsum %s has no description", db, dsf);
        ValNodeCopyStr (head, 0, buf);
        rsult = FALSE;
      }
      if (e2dsp->field_type < 0) {
        sprintf (buf, "Database %s docsum %s field type not indicated", db, dsf);
        ValNodeCopyStr (head, 0, buf);
        rsult = FALSE;
      }

    }
    if (e2db->docsum_field_count != dsfcount) {
      sprintf (buf, "Database %s docsum count %ld does not match lnkcount %d", db, (long) e2db->docsum_field_count, (int) dsfcount);
      ValNodeCopyStr (head, 0, buf);
      rsult = FALSE;
    }
  }

  if (e2ip->db_count != dbcount) {
    sprintf (buf, "Database count %ld does not match dbcount %d", (long) e2ip->db_count, (int) dbcount);
    ValNodeCopyStr (head, 0, buf);
    rsult = FALSE;
  }

  menuhead = ValNodeSort (menuhead, SortVnpByStr);

  last = NULL;
  lastvnp = NULL;
  for (vnp = menuhead; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    if (last != NULL && lastvnp != NULL) {
      if (StringICmp (last, str) == 0 && StringCmp (last, str) != 0) {
        if (StringICmp (last, "PmId") == 0 && StringICmp (str, "PMID") == 0) {
        } else if (StringICmp (last, "Object Type") == 0 && StringICmp (str, "Object type") == 0) {
          /* suppress for now */
        } else {
          sprintf (buf, "Menu names %s [%s] and %s [%s] differ in capitalization", last, dbnames [lastvnp->choice], str, dbnames [vnp->choice]);
          ValNodeCopyStr (head, 0, buf);
          rsult = FALSE;
        }
      } else if (checkMenuNameVariants) {
        len1 = StringLen (last);
        len2 = StringLen (str);
        if (len1 < len2) {
          if (StringNICmp (last, str, len1) == 0) {
            if (fullMenuVariantReport) {
              sprintf (buf, "Menu names %s [%s] and %s [%s] may be unintended variants", last, dbnames [lastvnp->choice], str, dbnames [vnp->choice]);
              ValNodeCopyStr (head, 0, buf);
              rsult = FALSE;
            } else if (StringICmp (last, "Gene Map") == 0 && StringICmp (str, "Gene Map Disorder") == 0) {
            } else if (StringICmp (last, "Organism") == 0 && StringICmp (str, "Organism unsynonymized") == 0) {
            } else if (StringICmp (last, "Reference") == 0 && StringICmp (str, "Reference Author") == 0) {
            } else if (StringICmp (last, "Reference") == 0 && StringICmp (str, "Reference SNP ID") == 0) {
            } else if (StringICmp (last, "Title") == 0 && StringICmp (str, "Title/Abstract") == 0) {
            } else if (StringICmp (last, "Rank") == 0 && StringICmp (str, "Ranked standard deviation") == 0) {
            } else if (StringICmp (last, "Book") == 0 && StringICmp (str, "Book's Topic") == 0) {
            } else if (StringICmp (last, "Gene Name") == 0 && StringICmp (str, "Gene Name or Description") == 0) {
            } else if (StringICmp (last, "Gene Name") == 0 && StringICmp (str, "Gene Name or Alias") == 0) {
            } else if (StringICmp (last, "Submitter") == 0 && StringICmp (str, "Submitter Handle") == 0) {
            } else if (StringICmp (last, "Abstract") == 0 && StringICmp (str, "Abstract/Index Tags") == 0) {
            } else if (StringICmp (last, "Author") == 0 && StringICmp (str, "Author Cluster ID") == 0) {
            } else if (StringICmp (last, "Author") == 0 && StringICmp (str, "Author Full Name") == 0) {
            } else if (StringICmp (last, "Expression") == 0 && StringICmp (str, "Expression Level") == 0) {
            } else if (StringICmp (last, "Chromosome") == 0 && StringICmp (str, "Chromosome GI") == 0) {
            } else if (StringICmp (last, "Disease") == 0 && StringICmp (str, "Disease or phenotype") == 0) {
            } else if (StringICmp (last, "GC") == 0 && StringICmp (str, "GC Content") == 0) {
            } else if (StringICmp (last, "Organism") == 0 && StringICmp (str, "Organism Motility") == 0) {
            } else if (StringICmp (last, "Publisher") == 0 && StringICmp (str, "Publisher ID") == 0) {
            } else if (StringICmp (last, "Disease") == 0 && StringICmp (str, "Disease-Stage") == 0) {
            } else if (StringICmp (last, "ActiveAid") == 0 && StringICmp (str, "ActiveAidCount") == 0) {
            } else if (StringICmp (last, "InactiveAid") == 0 && StringICmp (str, "InactiveAidCount") == 0) {
            } else if (StringICmp (last, "Phenotype") == 0 && StringICmp (str, "Phenotype Ontology ID") == 0) {
            } else if (StringICmp (last, "Title") == 0 && StringICmp (str, "Title Abbreviation") == 0) {
            } else if (StringICmp (last, "Library") == 0 && StringICmp (str, "Library Class") == 0) {
            } else if (StringICmp (last, "Sequence") == 0 && StringICmp (str, "Sequence Count") == 0) {
            } else if (StringICmp (last, "Journal") == 0 && StringICmp (str, "Journal List Identifier") == 0) {
            } else if (StringICmp (last, "CompoundID") == 0 && StringICmp (str, "CompoundIDActive") == 0) {
            } else if (StringICmp (last, "MeSHDescription") == 0 && StringICmp (str, "MeSHDescriptionActive") == 0) {
            } else if (StringICmp (last, "MeSHTerm") == 0 && StringICmp (str, "MeSHTermActive") == 0) {
            } else if (StringICmp (last, "PharmAction") == 0 && StringICmp (str, "PharmActionActive") == 0) {
            } else if (StringICmp (last, "SubstanceID") == 0 && StringICmp (str, "SubstanceIDActive") == 0) {
            } else if (StringICmp (last, "Synonym") == 0 && StringICmp (str, "SynonymActive") == 0) {
            } else if (StringICmp (last, "Definition") == 0 && StringICmp (str, "Definition Type") == 0) {
            } else if (StringICmp (last, "Reference") == 0 && StringICmp (str, "Reference Amino Acid") == 0) {
            } else if (StringICmp (last, "Reference SNP") == 0 && StringICmp (str, "Reference SNP ID") == 0) {
            } else if (StringICmp (last, "Analysis") == 0 && StringICmp (str, "Analysis ID") == 0) {
            } else if (StringICmp (last, "Document") == 0 && StringICmp (str, "Document ID") == 0) {
            } else if (StringICmp (last, "Study") == 0 && StringICmp (str, "Study Accession") == 0) {
            } else if (StringICmp (last, "Study") == 0 && StringICmp (str, "Study ID") == 0) {
            } else if (StringICmp (last, "Variable") == 0 && StringICmp (str, "Variable ID") == 0) {
            } else if (StringICmp (last, "COG") == 0 && StringICmp (str, "COG group") == 0) {
            } else if (StringICmp (last, "Locus Tag") == 0 && StringICmp (str, "Locus Tag Prefix") == 0) {
            } else if (StringICmp (last, "Attribute") == 0 && StringICmp (str, "Attributes") == 0) {
            } else if (StringICmp (last, "Genotype") == 0 && StringICmp (str, "Genotype Platform") == 0) {
            } else if (StringICmp (last, "Group") == 0 && StringICmp (str, "Group ID") == 0) {
            } else if (StringICmp (last, "Clinical Synopsis") == 0 && StringICmp (str, "Clinical Synopsis Date") == 0) {
            } else if (StringICmp (last, "Volume") == 0 && StringICmp (str, "Volume3D") == 0) {
            } else if (StringICmp (last, "InChI") == 0 && StringICmp (str, "InChIKey") == 0) {
            } else if (StringICmp (last, "Dataset") == 0 && StringICmp (str, "Dataset ID") == 0) {
            } else if (StringICmp (last, "Comment") == 0 && StringICmp (str, "Comments") == 0) {
            } else if (StringICmp (last, "SID") == 0 && StringICmp (str, "SidExternalID") == 0) {
            } else if (StringICmp (last, "Platform") == 0 && StringICmp (str, "Platform Reporter Type") == 0) {
            } else if (StringICmp (last, "Database") == 0 && StringICmp (str, "Database Name") == 0) {
            } else if (StringICmp (last, "Date") == 0 && StringICmp (str, "Date Discontinued") == 0) {
            } else if (StringICmp (last, "Accession") == 0 && StringICmp (str, "AccessionID") == 0) {
            } else if (StringICmp (last, "CategorizedComment") == 0 && StringICmp (str, "CategorizedCommentTitle") == 0) {
            } else if (StringICmp (last, "Journal") == 0 && StringICmp (str, "JournalName") == 0) {
            } else if (StringICmp (last, "Project") == 0 && StringICmp (str, "Project ID") == 0) {
            } else if (StringICmp (last, "Study") == 0 && StringICmp (str, "Study Has SRA Components") == 0) {
            } else if (StringICmp (last, "Subject") == 0 && StringICmp (str, "Subject Terms") == 0) {
            } else if (StringICmp (last, "Accession") == 0 && StringICmp (str, "Accession Version") == 0) {
            } else if (StringICmp (last, "Allele") == 0 && StringICmp (str, "Allele Origin") == 0) {
            } else if (StringICmp (last, "Allele") == 0 && StringICmp (str, "Allele Type") == 0) {
            } else if (StringICmp (last, "Attribute") == 0 && StringICmp (str, "Attribute Name") == 0) {
            } else if (StringICmp (last, "Chromosome") == 0 && StringICmp (str, "Chromosome Accession") == 0) {
            } else if (StringICmp (last, "Figure Caption") == 0 && StringICmp (str, "Figure Caption Title") == 0) {
            } else if (StringICmp (last, "Genome Project") == 0 && StringICmp (str, "Genome Projects ID") == 0) {
            } else if (StringICmp (last, "Journal") == 0 && StringICmp (str, "Journal Author") == 0) {
            } else if (StringICmp (last, "Library") == 0 && StringICmp (str, "Library Abbreviation") == 0) {
            } else if (StringICmp (last, "Method Type") == 0 && StringICmp (str, "Method Type Category") == 0) {
            } else if (StringICmp (last, "Method Type") == 0 && StringICmp (str, "Method Type Weight") == 0) {
            } else if (StringICmp (last, "Sample") == 0 && StringICmp (str, "Sample Count") == 0) {
            } else if (StringICmp (last, "Study") == 0 && StringICmp (str, "Study Sccession") == 0) {
            } else if (StringICmp (last, "Subject") == 0 && StringICmp (str, "Subject Phenotype Status") == 0) {
            } else if (StringICmp (last, "Validation Result") == 0 && StringICmp (str, "Validation Result Weight") == 0) {
            } else if (StringICmp (last, "Alias") == 0 && StringICmp (str, "Alias for a variant") == 0) {
            } else if (StringICmp (last, "CID") == 0 && StringICmp (str, "CID Count") == 0) {
            } else if (StringICmp (last, "SID") == 0 && StringICmp (str, "SID Count") == 0) {
            } else if (StringICmp (last, "BioProject") == 0 && StringICmp (str, "BioProject ID") == 0) {
            } else if (StringICmp (last, "Project") == 0 && StringICmp (str, "Project Accession") == 0) {
            } else if (StringICmp (last, "Last Update") == 0 && StringICmp (str, "Last Update Date") == 0) {
            } else if (StringICmp (last, "Accession") == 0 && StringICmp (str, "Accession ID") == 0) {
            } else if (StringICmp (last, "Book") == 0 && StringICmp (str, "Book Accession ID") == 0) {
            } else if (StringICmp (last, "Author") == 0 && StringICmp (str, "Author - Corporate") == 0) {
            } else if (StringICmp (last, "Investigator") == 0 && StringICmp (str, "Investigator - Full") == 0) {
            } else if (StringICmp (last, "Submission") == 0 && StringICmp (str, "Submission Date") == 0) {
            } else if (StringICmp (last, "Categorized Comment") == 0 && StringICmp (str, "Categorized Comment Title") == 0) {
            } else if (StringICmp (last, "Journal Publication date") == 0 && StringICmp (str, "Journal Publication Date") == 0) {
            } else if (StringICmp (last, "PharmAction") == 0 && StringICmp (str, "PharmActionID") == 0) {
            } else if (StringICmp (last, "Synonym") == 0 && StringICmp (str, "Synonym Active") == 0) {
            } else if (StringICmp (last, "Categorized Comment") == 0 && StringICmp (str, "Categorized Comment Title") == 0) {
            } else if (StringICmp (last, "Journal") == 0 && StringICmp (str, "Journal Name") == 0) {
            } else if (StringICmp (last, "Submitter") == 0 && StringICmp (str, "Submitter Affiliation") == 0) {
            } else if (StringICmp (last, "Clone Name") == 0 && StringICmp (str, "Clone Name Alias") == 0) {
            } else if (StringICmp (last, "Cultivar") == 0 && StringICmp (str, "Cultivar Accession") == 0) {
            } else if (StringICmp (last, "Organ") == 0 && StringICmp (str, "Organism") == 0) {
            } else if (StringICmp (last, "Placed") == 0 && StringICmp (str, "Placed Scaffolds Count") == 0) {
            } else if (StringICmp (last, "Population") == 0 && StringICmp (str, "Population Class") == 0) {
            } else if (StringICmp (last, "Vector") == 0 && StringICmp (str, "Vector Type") == 0) {
            } else if (StringICmp (last, "Disease") == 0 && StringICmp (str, "Disease/Phenotype") == 0) {
            } else if (StringICmp (last, "Variable") == 0 && StringICmp (str, "Variable Description") == 0) {
            } else {
              sprintf (buf, "Menu names %s [%s] and %s [%s] may be unintended variants", last, dbnames [lastvnp->choice], str, dbnames [vnp->choice]);
              ValNodeCopyStr (head, 0, buf);
              rsult = FALSE;
            }
          }
        }
      }
    } else if (StringICmp (str, "Title Word") == 0) {
      sprintf (buf, "Menu name Title Word should be replaced by Title");
      ValNodeCopyStr (head, 0, buf);
      rsult = FALSE;
    }
    last = str;
    lastvnp = vnp;
  }

  if (checkMenuNameFormat) {
    for (vnp = menuhead; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      if (LowerCaseWords (str)) {
        sprintf (buf, "Lower-case words in field %s [%s]", str, dbnames [vnp->choice]);
        ValNodeCopyStr (head, 0, buf);
        rsult = FALSE;
      }
    }
    for (vnp = menuhead; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      if (MixedCaseWords (str)) {
        sprintf (buf, "Mixed-case words in field %s [%s]", str, dbnames [vnp->choice]);
        ValNodeCopyStr (head, 0, buf);
        rsult = FALSE;
      }
    }
  }

  ValNodeFreeData (menuhead);

  return rsult;
}

NLM_EXTERN Boolean ValidateEntrez2InfoPtrExEx (
  Entrez2InfoPtr e2ip,
  ValNodePtr PNTR head,
  Boolean checkMenuNameVariants,
  Boolean checkMenuNameFormat
)

{
  return ValidateEntrez2InfoPtrExExEx (e2ip, head, checkMenuNameVariants, checkMenuNameFormat, FALSE);
}

NLM_EXTERN Boolean ValidateEntrez2InfoPtrEx (
  Entrez2InfoPtr e2ip,
  ValNodePtr PNTR head,
  Boolean checkMenuNameVariants
)

{
  return ValidateEntrez2InfoPtrExExEx (e2ip, head, checkMenuNameVariants, FALSE, FALSE);
}

NLM_EXTERN Boolean ValidateEntrez2InfoPtr (
  Entrez2InfoPtr e2ip,
  ValNodePtr PNTR head
)

{
  return ValidateEntrez2InfoPtrExExEx (e2ip, head, FALSE, FALSE, FALSE);
}

/* network connection test functions */

static CONN NetTestOpenConnection (void)

{
  char        buffer [64];
  CONN        conn;
  size_t      n_written;
  EIO_Status  status;

  conn = QUERY_OpenUrlQuery ("www.ncbi.nlm.nih.gov", 80, "/Service/bounce.cgi",
                             NULL, "Entrez2Tool", 0, eMIME_T_Text,
                             eMIME_Plain, eENCOD_None, 0);
  if (conn == NULL) return NULL;

  sprintf (buffer, "test\n");
  buffer [4] = '\012';
  status = CONN_Write (conn, (const void *) buffer, StringLen (buffer),
                       &n_written, eIO_WritePersist);
  if (status != eIO_Success) {
    CONN_Close (conn);
    return NULL;
  }

  return conn;
}

NLM_EXTERN Boolean NetTestAsynchronousQuery (
  QUEUE* queue,
  QueryResultProc resultproc,
  VoidPtr userdata
)

{
  CONN  conn;

  conn = NetTestOpenConnection ();

  if (conn == NULL) return FALSE;

  QUERY_SendQuery (conn);

  QUERY_AddToQueue (queue, conn, resultproc, userdata, TRUE);

  return TRUE;
}

NLM_EXTERN Boolean NetTestReadReply (
  CONN conn,
  EIO_Status status
)

{
  char         buffer [64];
  size_t       n_read;
  ErrSev       oldsev;
  Boolean      res = FALSE;

  if (conn != NULL && status == eIO_Success) {
    oldsev = ErrSetMessageLevel (SEV_MAX);
    status = CONN_Read (conn, buffer, sizeof (buffer), &n_read, eIO_ReadPlain);
    if (status == eIO_Success) {
      if (StringNCmp (buffer, "test", 4) == 0) {
        res = TRUE;
      }
    }
    ErrSetMessageLevel (oldsev);
  }
  return res;
}

NLM_EXTERN Int4 NetTestCheckQueue (
  QUEUE* queue
)

{
  return QUERY_CheckQueue (queue);
}

