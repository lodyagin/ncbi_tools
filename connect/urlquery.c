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
* $Revision: 6.1 $
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

#include <urlquery.h>

NLM_EXTERN CONN QUERY_OpenUrlQuery (
  Nlm_CharPtr host_machine, Nlm_Uint2 host_port,
  Nlm_CharPtr host_path, Nlm_CharPtr arguments,
  Nlm_CharPtr appName, Nlm_Uint4 timeoutsec,
  UrlContentType ctype, URLC_Flags flags
)

{
  CONN           conn = NULL;
  CONNECTOR      connector;
  Nlm_CharPtr    contentType = "application/x-www-form-urlencoded";
  SNetConnInfo*  info;
  Nlm_Uint4      n_written;
  Nlm_Char       path [PATH_MAX];
  EConnStatus    status;
  Nlm_CharPtr    supp_input = NULL;
  Nlm_CharPtr    userAgentName = NULL;
  Nlm_Char       user_header [256];

  if (StringHasNoText (host_machine) || StringHasNoText (host_path)) return NULL;

  /* allow the user to specify a program name, otherwise get from ProgramPath */

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

  /* set content type from parameter */

  if (ctype == wwwencoded) {
    contentType = "application/x-www-form-urlencoded";
  } else if (ctype == multipart) {
    contentType = "multipart/form-data";
  }

  /* set HTML header with program name as user agent */

  sprintf (user_header, "Content-type: %s\n"
           "Content-Disposition: form-data; name=supplemental_input\n"
           "User-Agent: %s\n",
           contentType, userAgentName);

  /* fill in connection info fields */

  info = NetConnInfo_Create (0, 0);
  if (info == NULL) return NULL;

  info->host = StringSave (host_machine);
  info->port = host_port;

  /* no longer using info->args */

  info->path = MemFree (info->path);
  info->path = StringSave (host_path);

  info->timeout.sec  = timeoutsec;
  info->timeout.usec = 0;

  connector = URL_CreateConnector (info, user_header, flags);

  status = CONN_Create (connector, &conn);

  NetConnInfo_Destroy (&info);

  if (status != eCONN_Success) {
    ErrPostEx (SEV_ERROR, 0, 0, "QUERY_OpenUrlQuery failed in CONN_Create");
  }

  /* POST arguments in data buffer followed by supplemental_input=\n, but do not URL encode this */

  if (! StringHasNoText (arguments)) {
    status = CONN_Write (conn, (const void *) arguments, StringLen (arguments), &n_written);
    supp_input = "&supplemental_input=\n";
  } else {
    supp_input = "supplemental_input=\n";
  }
   status = CONN_Write (conn, (const void *) supp_input, StringLen (supp_input), &n_written);

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

  buffer = MemNew (URL_QUERY_BUFLEN + 1);
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

  buffer = MemNew (URL_QUERY_BUFLEN + 1);
  if (buffer != NULL) {
    while ((status = CONN_Read (conn, buffer, URL_QUERY_BUFLEN, &n_read, eCR_Read)) == eCONN_Success) {
      FileWrite (buffer, 1, n_read, fp);
    }
  }
  MemFree (buffer);
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

