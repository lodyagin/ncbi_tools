/*   urlquery.c
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
* File Name:  urlquery.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   4/16/98
*
* $Revision: 6.6 $
*
* File Description: 
*
* Modifications:  
* --------------------------------------------------------------------------
* $Log: urlquery.c,v $
* Revision 6.6  1999/07/28 21:09:15  vakatov
* Multiple fixes in QUERY_OpenUrlQuery() to make it work with a generic
* URL server;  also, pass arguments in the cmd.-line
*
* ==========================================================================
*/

#include "asnbuild.h"
#include <urlquery.h>


NLM_EXTERN CONN QUERY_OpenUrlQuery (
  Nlm_CharPtr host_machine, Nlm_Uint2 host_port,
  Nlm_CharPtr host_path, Nlm_CharPtr arguments,
  Nlm_CharPtr appName, Nlm_Uint4 timeoutsec,
  EMIME_SubType subtype, URLC_Flags flags)
{
  CONN            conn = NULL;
  CONNECTOR       connector;
  Nlm_Char        contentType[MAX_CONTENT_TYPE_LEN];
  EMIME_Encoding  encoding = eENCOD_None;
  SNetConnInfo*   info;
  Nlm_Char        path [PATH_MAX];
  EConnStatus     status;
  Nlm_CharPtr     userAgentName = NULL;
  Nlm_Char        user_header[sizeof(contentType) + 256];

  if (StringHasNoText (host_machine) || StringHasNoText (host_path))
    return NULL;

  /* allow the user to specify a prog. name, otherwise get from ProgramPath */
  if (! StringHasNoText (appName)) {
    userAgentName = appName;
  } else {
    Nlm_ProgramPath (path, sizeof (path));
    userAgentName = StringRChr (path, DIRDELIMCHR);
    if (userAgentName != NULL) {
      userAgentName++;
    }
  }
  if (StringHasNoText (userAgentName)) {
    userAgentName = "?";
  }

  /* set content type from parameters */
  if ((flags & URLC_URL_ENCODE_OUT) != 0) {
    encoding = eENCOD_Url;
  }
  VERIFY( MIME_ComposeContentType(subtype, encoding,
                                  contentType, sizeof(contentType)) );

  /* set HTML header with program name as user agent */
  sprintf (user_header, "%sUser-Agent: %s\r\n", contentType, userAgentName);

  /* fill in connection info fields and create the connection */
  info = NetConnInfo_Create(0, 0);
  ASSERT( info );

  if ( host_machine ) {
    MemFree(info->host);  info->host = StringSave(host_machine);
  }
  info->port = host_port;
  if ( !StringHasNoText(arguments) ) {
    MemFree(info->args);  info->args = StringSave(arguments);
  }
  MemFree(info->path);    info->path = StringSave(host_path);

  info->timeout.sec  = timeoutsec;
  info->timeout.usec = 0;

  connector = URL_CreateConnector(info, user_header, flags);
  status = CONN_Create(connector, &conn);
  if (status != eCONN_Success) {
    ErrPostEx(SEV_ERROR, 0, 0, "QUERY_OpenUrlQuery failed in CONN_Create");
  }

  /* cleanup & return */
  NetConnInfo_Destroy(&info);
  return conn;
}


NLM_EXTERN void QUERY_SendQuery (
  CONN conn
)

{
  STimeout  timeout;

  if (conn == NULL) return;

  /* flush buffer, sending query, without waiting for response */

  timeout.sec  = 0;
  timeout.usec = 0;
  CONN_Wait (conn, eCONN_Read, &timeout);
}

#define URL_QUERY_BUFLEN  4096

NLM_EXTERN void QUERY_CopyFileToQuery (
  CONN conn, FILE *fp
)

{
  Nlm_CharPtr  buffer;
  size_t       ct;
  Nlm_Uint4    n_written;
  EConnStatus  status;

  if (conn == NULL || fp == NULL) return;

  buffer = (Nlm_CharPtr) MemNew(URL_QUERY_BUFLEN + 1);
  if (buffer != NULL) {
    while ((ct = FileRead (buffer, 1, URL_QUERY_BUFLEN, fp)) > 0) {
      status = CONN_Write (conn, (const void *) buffer, ct, &n_written);
    }
  }
  MemFree (buffer);
}

NLM_EXTERN void QUERY_CopyResultsToFile (
  CONN conn, FILE *fp
)

{
  Nlm_CharPtr  buffer;
  Nlm_Uint4    n_read;
  EConnStatus  status;

  if (conn == NULL || fp == NULL) return;

  buffer = (Nlm_CharPtr) MemNew(URL_QUERY_BUFLEN + 1);
  if (buffer != NULL) {
    while ((status = CONN_Read (conn, buffer, URL_QUERY_BUFLEN, &n_read, eCR_Read)) == eCONN_Success) {
      FileWrite (buffer, 1, n_read, fp);
    }
  }
  MemFree (buffer);
}

static Nlm_Int2 LIBCALL AsnIoConnWrite (Pointer ptr, Nlm_CharPtr buf, Nlm_Uint2 count)

{
	Nlm_Uint4     bytes;
	AsnIoConnPtr  aicp;

	aicp = (AsnIoConnPtr) ptr;
	CONN_Write (aicp->conn, (const void *) buf, (size_t) count, &bytes);
	return (Nlm_Int2) bytes;
}

static Nlm_Int2 LIBCALL AsnIoConnRead (Pointer ptr, CharPtr buf, Nlm_Uint2 count)

{
	Nlm_Uint4     bytes;
	AsnIoConnPtr  aicp;

	aicp = (AsnIoConnPtr) ptr;
	CONN_Read (aicp->conn, (Pointer) buf, (Int4) count, &bytes, eCR_Read);
	return (Nlm_Int2) bytes;
}

NLM_EXTERN AsnIoConnPtr QUERY_AsnIoConnOpen (Nlm_CharPtr mode, CONN conn)

{
  Int1          type;
  AsnIoConnPtr  aicp;

  if (! StringCmp(mode, "r"))
    type = (ASNIO_IN | ASNIO_TEXT);
  else if (! StringCmp(mode, "rb"))
    type = (ASNIO_IN | ASNIO_BIN);
  else if (! StringCmp(mode, "w"))
    type = (ASNIO_OUT | ASNIO_TEXT);
  else if (! StringCmp(mode, "wb"))
    type = (ASNIO_OUT | ASNIO_BIN);
  else
  {
    AsnIoErrorMsg (NULL, 81, mode);
    return NULL;
  }

  aicp = (AsnIoConnPtr) MemNew (sizeof (AsnIoConn));
  aicp->aip = AsnIoNew (type, NULL, (Pointer) aicp, AsnIoConnRead, AsnIoConnWrite);
  aicp->conn = conn;
  return aicp;
}

NLM_EXTERN AsnIoConnPtr QUERY_AsnIoConnClose (AsnIoConnPtr aicp)

{
  if (aicp == NULL) return NULL;
  AsnIoClose (aicp->aip);
  return (AsnIoConnPtr) MemFree (aicp);
}

typedef struct SQueueTag {
  CONN                conn;
  QueryResultProc     resultproc;
  Nlm_VoidPtr         userdata;
  struct SQueueTag*   next;
} SConnQueue, PNTR QueuePtr;

NLM_EXTERN void QUERY_AddToQueue (
  QUEUE* queue, CONN conn, QueryResultProc resultproc, Nlm_VoidPtr userdata
)

{
  QueuePtr       cqp;
  QueuePtr PNTR  qptr;
  QueuePtr       tmp;

  if (conn == NULL || resultproc == NULL) return;

  /* allocate queue element */

  cqp = (QueuePtr) MemNew (sizeof (SConnQueue));
  if (cqp == NULL) return;

  cqp->conn = conn;
  cqp->resultproc = resultproc;
  cqp->userdata = userdata;

  /* add to polling queue */

  qptr = (QueuePtr PNTR) queue;
  if (qptr != NULL) {
    if (*qptr != NULL) {
      tmp = *qptr;
      if (tmp != NULL) {
        while (tmp->next != NULL) {
          tmp = tmp->next;
        }
        tmp->next = cqp;
      }
    } else {
      *qptr = cqp;
    }
  }
}

static void QUERY_RemoveFromQueue (
  QUEUE* queue, CONN conn
)

{
  QueuePtr       curr;
  QueuePtr       next;
  QueuePtr PNTR  prev;
  QueuePtr PNTR  qptr;

  qptr = (QueuePtr PNTR) queue;
  if (qptr == NULL || *qptr == NULL || conn == NULL) return;

  prev = qptr;
  curr = *qptr;

  while (curr != NULL) {
    next = curr->next;
    if (curr->conn == conn) {
      *(prev) = next;
      curr->next = NULL;
      MemFree (curr);
    } else {
      prev = &(curr->next);
    }
    curr = next;
  }
}

NLM_EXTERN Nlm_Int4 QUERY_CheckQueue (
  QUEUE* queue
)

{
  Nlm_Int4       count = 0;
  QueuePtr       curr;
  QueuePtr       next;
  QueuePtr PNTR  qptr;
  EConnStatus    status;
  STimeout       timeout;

  qptr = (QueuePtr PNTR) queue;
  if (qptr == NULL || *qptr == NULL) return 0;

  curr = *qptr;

  while (curr != NULL) {
    next = curr->next;

    if (curr->conn != NULL) {
      timeout.sec  = 0;
      timeout.usec = 0;
      status = CONN_Wait (curr->conn, eCONN_Read, &timeout);

      if (status == eCONN_Success) {
        if (curr->resultproc != NULL) {
          /* result could eventually be used to reconnect on timeout */
          curr->resultproc (curr->conn, curr->userdata, status);
        }
        CONN_Close (curr->conn);
        QUERY_RemoveFromQueue (queue, curr->conn);

      } else {
        count++;
      }
    }

    curr = next;
  }

  return count;
}

