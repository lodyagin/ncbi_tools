/*  $Id: connectn.c,v 6.3 1999/04/05 15:32:53 vakatov Exp $
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================
*
* Author:  Denis Vakatov
*
* File Description:
*   Generic API to open and handle connection to an abstract service.
*   For more detail, see in "connectn.h".
*
* --------------------------------------------------------------------------
* $Log: connectn.c,v $
* Revision 6.3  1999/04/05 15:32:53  vakatov
* CONN_Wait():  be more mild and discrete about the posted error severity
*
* Revision 6.2  1999/04/02 20:41:53  vakatov
* Got rid of occasional comment inside comment(in CVS Log)
*
* Revision 6.1  1999/04/01 21:48:09  vakatov
* Fixed for the change in spec:  "n_written/n_read" args in
* CONN_Write/Read to be non-NULL and "*n_written / *n_read" := 0
*
* Revision 6.0  1999/03/25 23:04:57  vakatov
* Initial revision
*
* ==========================================================================
*/

#include <ncbi.h>
#include <connectn.h>
#include <connectr.h>
#include <ncbibuf.h>


/***********************************************************************
 *  INTERNAL
 ***********************************************************************/


/* Connection internal data
 */
typedef struct SConnectionTag {
  CONNECTOR              connector;  /* current connector */
  BUF                    buf;        /* storage for the Peek'd data */
#ifdef IMPLEMENTED__CONN_WaitAsync
  SConnectorAsyncHandler async_data; /* info on the curr async event handler */
#endif
  Nlm_Boolean            connected;  /* if it is already connected */
  /* "[c|r|w]_timeout" is either NULL or points to "[cc|rr|ww]_timeout" */
  STimeout* c_timeout;  /* connection timeout on connect */
  STimeout* r_timeout;  /* connection timeout on reading */
  STimeout* w_timeout;  /* connection timeout on writing */
  STimeout* l_timeout;  /* connection timeout on close */
  STimeout  cc_timeout; /* storage for "c_timeout" */
  STimeout  rr_timeout; /* storage for "r_timeout" */
  STimeout  ww_timeout; /* storage for "w_timeout" */
  STimeout  ll_timeout; /* storage for "l_timeout" */
} SConnection;


/* Standard error report
 */

#define s_ErrPost(severity, sub_code, descr, connector, status) \
  (void)((Nlm_ErrSetContext("",__FILE__,__LINE__,DBFLAG,0,0,0) != 0) ? \
  0 : sf_ErrPost(severity, sub_code, descr, connector, status))

static int sf_ErrPost(ErrSev severity, int sub_code, const Nlm_Char* descr,
                       CONNECTOR connector, EConnStatus status)
{
  const Nlm_Char* x_descr = descr ? descr : "undescribed error";
  const Nlm_Char* x_type  = (connector && connector->vtable.get_type) ?
    (*connector->vtable.get_type)(connector->handle) : "Unknown";

  Nlm_ErrPostEx(severity, CONN_ERRCODE, sub_code,
                "%s (connector \"%s\", error \"%s\")",
                x_descr, x_type, CONN_StatusString(status));
  return 0;
}



/***********************************************************************
 *  EXTERNAL
 ***********************************************************************/


NLM_EXTERN const Nlm_Char* CONN_StatusString
(EConnStatus status)
{
  static const Nlm_Char* s_StatusStr[] = {
    "Success", "Timeout", "Closed", "InvalidArg", "NotSupported", "Unknown"
  };
  return s_StatusStr[status];
}


NLM_EXTERN EConnStatus CONN_Create
(CONNECTOR connector,
 CONN* conn)
{
  *conn = (SConnection*)Nlm_MemNew(sizeof(SConnection));

  (*conn)->connector = connector;
  /* all other fields are filled with zero */

  return eCONN_Success;
}


NLM_EXTERN EConnStatus CONN_Reconnect
(CONN conn,
 CONNECTOR connector)
{
  EConnStatus status;
  CONNECTOR   x_connector = conn->connector;

  /* check arg */
  if (!connector  &&  !x_connector) {
    status = eCONN_Unknown;
    s_ErrPost(SEV_ERROR, 20,
              "[CONN_Connect]  Both current and new connectors are NULLs",
              x_connector, status);
    return status;    
  }

  /* reset and close current connector, if any */
  if ( x_connector ) {
#ifdef IMPLEMENTED__CONN_WaitAsync
    /* cancel async. i/o event handler */
    CONN_WaitAsync(conn, eCONN_ReadWrite, 0, 0, 0);
#endif
    /* flush unwritten data */
    CONN_Flush(conn);
    {{ /* erase unread data */
      Nlm_Uint4 buf_size = BUF_Size(conn->buf);
      VERIFY( BUF_Read(conn->buf, 0, buf_size) == buf_size );
    }}
    /* call current connector's "CLOSE" method */
    if (x_connector != connector  &&  x_connector->vtable.close) {
      status = (*x_connector->vtable.close)(x_connector, conn->l_timeout);
      if (status != eCONN_Success) {
        s_ErrPost(SEV_ERROR, 21, "[CONN_Connect]  Cannot close connection",
                  x_connector, status);
        return status;
      }
      conn->connector = 0;
    }
  }

  /* success -- setup the new connector, and schedule for the future connect */
  if ( connector )
    x_connector = conn->connector = connector;
  conn->connected = FALSE;

  return eCONN_Success;
}


/* Perform a "real" connect
 */
static EConnStatus s_Connect(CONN conn)
{
  CONNECTOR   x_connector = conn->connector;
  EConnStatus status;

  ASSERT( !conn->connected );
  if ( !x_connector ) {
    status = eCONN_Unknown;
    s_ErrPost(SEV_ERROR, 25, "[CONN_Connect]  Cannot connect, NULL connector",
              x_connector, status);
    return status;
  }

  /* call current connector's "CONNECT" method */
  status = x_connector->vtable.connect ?
    (*x_connector->vtable.connect)(x_connector->handle, conn->c_timeout) :
    eCONN_NotSupported;
  if (status != eCONN_Success) {
    s_ErrPost(SEV_ERROR, 26, "[CONN_Connect]  Cannot connect",
              x_connector, status);
    return status;
  }

  /* success */
  conn->connected = TRUE;
  return eCONN_Success;
}


NLM_EXTERN void CONN_SetTimeout
(CONN            conn,
 EConnDirection  direction,
 const STimeout* new_timeout)
{
  if (direction == eCONN_Connect) {
    if ( new_timeout ) {
      conn->cc_timeout = *new_timeout;
      conn->c_timeout  = &conn->cc_timeout;
    }
    else
      conn->c_timeout = 0;

    return;
  }

  if (direction == eCONN_Close) {
    if ( new_timeout ) {
      conn->ll_timeout = *new_timeout;
      conn->l_timeout  = &conn->ll_timeout;
    }
    else
      conn->l_timeout = 0;

    return;
  }

  if (direction == eCONN_ReadWrite  ||  direction == eCONN_Read) {
    if ( new_timeout ) {
      conn->rr_timeout = *new_timeout;
      conn->r_timeout  = &conn->rr_timeout;
    }
    else
      conn->r_timeout = 0;
  }

  if (direction == eCONN_ReadWrite  ||  direction == eCONN_Write) {
    if ( new_timeout ) {
      conn->ww_timeout = *new_timeout;
      conn->w_timeout  = &conn->ww_timeout;
    }
    else
      conn->w_timeout = 0;
  }
}


NLM_EXTERN const STimeout* CONN_GetTimeout
(CONN           conn,
 EConnDirection direction)
{
  switch ( direction ) {
    case eCONN_Connect:
      return conn->c_timeout;
    case eCONN_ReadWrite:
    case eCONN_Read:
      return conn->r_timeout;
    case eCONN_Write:
      return conn->w_timeout;
    case eCONN_Close:
      return conn->l_timeout;
  }
  ASSERT(0);  return 0;
}


NLM_EXTERN EConnStatus CONN_Wait
(CONN            conn,
 EConnDirection  direction,
 const STimeout* timeout)
{
  CONNECTOR x_connector = conn->connector;
  EConnStatus status;

  /* check args */
  if (direction != eCONN_Read  &&  direction != eCONN_Write)
    return eCONN_InvalidArg;

  /* check if there is a PEEK'ed data in the input */
  if (direction == eCONN_Read  &&  BUF_Size(conn->buf))
    return eCONN_Success;

  /* perform connect, if not connected yet */
  if (!conn->connected  &&  (status = s_Connect(conn)) != eCONN_Success)
    return status;

  /* call current connector's "WAIT" method */
  status = x_connector->vtable.wait ?
    (*x_connector->vtable.wait)(x_connector->handle, direction, timeout) :
    eCONN_NotSupported;

  if (status != eCONN_Success) {
    ErrSev severity = SEV_ERROR;
    if (status == eCONN_Timeout) {
      severity = (timeout  &&  !timeout->sec  &&  !timeout->usec) ?
        SEV_INFO : SEV_WARNING;
    }
    s_ErrPost(severity, 40, "[CONN_Wait]  Error blocking on i/o",
              x_connector, status);
  }

  return status;
}


#ifdef IMPLEMENTED__CONN_WaitAsync
/* Internal handler(wrapper for the user-provided handler) for CONN_WaitAsync()
 */
static void s_ConnectorAsyncHandler
(SConnectorAsyncHandler* data,
 EConnDirection          direction,
 EConnStatus             status)
{
  /* handle the async. event */
  (*data->handler)(data->conn, direction, status, data->data);

  /* reset */
  VERIFY( CONN_WaitAsync(data->conn, eCONN_ReadWrite, 0, 0, 0) ==
          eCONN_Success );
}


NLM_EXTERN EConnStatus CONN_WaitAsync
(CONN              conn,
 EConnDirection    direction,
 FConnAsyncHandler handler,
 void*             data,
 FConnAsyncCleanup cleanup)
{
  EConnStatus status;
  CONNECTOR x_connector = conn->connector;
  SConnectorAsyncHandler* x_data = &conn->async_data;

  /* perform connect, if not connected yet */
  if (!conn->connected  &&  (status = s_Connect(conn)) != eCONN_Success)
    return status;

  /* reset previous handler, cleanup its data */
  /* (call current connector's "WAIT_ASYNC" method with NULLs) */
  status = x_connector->vtable.wait_async ?
    (*x_connector->vtable.wait_async)(x_connector->handle, 0, 0) :
    eCONN_NotSupported;
  if (status != eCONN_Success) {
    s_ErrPost(SEV_ERROR, 50, "[CONN_WaitAsync]  Cannot reset the handler",
              x_connector, status);
    return status;
  }
  if ( x_data->cleanup )
    (*x_data->cleanup)(x_data->data);
  Nlm_MemSet(x_data, '\0', sizeof(*x_data));

  /* set new handler, if specified */
  /* (call current connector's "WAIT_ASYNC" method with new handler/data) */
  if ( !handler )
    return eCONN_Success;

  x_data->conn           = conn;
  x_data->wait_direction = direction;
  x_data->handler        = handler;
  x_data->data           = data;
  x_data->cleanup        = cleanup;

  status = (*x_connector->vtable.wait_async)(x_connector->handle,
                                             s_ConnectorAsyncHandler, x_data);
  if (status != eCONN_Success)
    s_ErrPost(SEV_ERROR, 51, "[CONN_WaitAsync]  Cannot set new handler",
              x_connector, status);
  return status;
}
#endif /* IMPLEMENTED__CONN_WaitAsync */


NLM_EXTERN EConnStatus CONN_Write
(CONN        conn,
 const void* buf,
 Nlm_Uint4   size,
 Nlm_Uint4*  n_written)
{
  CONNECTOR x_connector = conn->connector;
  EConnStatus status;

  *n_written = 0;

  /* perform connect, if not connected yet */
  if (!conn->connected  &&  (status = s_Connect(conn)) != eCONN_Success)
    return status;

  /* call current connector's "WRITE" method */
  status = x_connector->vtable.write ?
    (*x_connector->vtable.write)(x_connector->handle, buf, size, n_written,
                                 conn->w_timeout) :
    eCONN_NotSupported;

  if (status != eCONN_Success)
    s_ErrPost(SEV_ERROR, 60, "[CONN_Write]  Writing error",
              x_connector, status);
  return status;
}


NLM_EXTERN EConnStatus CONN_Flush
(CONN conn)
{
  CONNECTOR x_connector = conn->connector;
  EConnStatus status;

  /* do nothing, if not connected yet */
  if ( !conn->connected )
    return eCONN_Success;

  /* call current connector's "FLUSH" method */
  status = x_connector->vtable.flush ?
    (*x_connector->vtable.flush)(x_connector->handle, conn->w_timeout) :
    eCONN_NotSupported;

  if (status != eCONN_Success)
    s_ErrPost(SEV_WARNING, 70, "[CONN_Flush]  Cannot flush data",
              x_connector, status);
  return status;
}


/* Read or peek data from the input queue
 * See CONN_Read()
 */
static EConnStatus s_CONN_Read
(CONN        conn,
 void*       buf,
 Nlm_Uint4   size,
 Nlm_Uint4*  n_read,
 Nlm_Boolean peek)
{
  CONNECTOR x_connector = conn->connector;
  EConnStatus status;

  /* perform connect, if not connected yet */
  if (!conn->connected  &&  (status = s_Connect(conn)) != eCONN_Success)
    return status;

  /* flush the unwritten output data, if any */
  CONN_Flush(conn);

  /* check if the read method is specified at all */
  if ( !x_connector->vtable.read ) {
    status = eCONN_NotSupported;
    s_ErrPost(SEV_ERROR, 80, "[CONN_Read]  Cannot read data",
              x_connector, status);
    return status;
  }

  /* read data from the internal "peek-buffer", if any */
  *n_read = peek ?
    BUF_Peek(conn->buf, buf, size) : BUF_Read(conn->buf, buf, size);
  if (*n_read == size)
    return eCONN_Success;
  buf = (char*)buf + *n_read;
  size -= *n_read;

  /* read data from the connection */
  {{
    Nlm_Uint4 x_read = 0;
    /* call current connector's "READ" method */
    status = (*x_connector->vtable.read)(x_connector->handle,
                                         buf, size, &x_read, conn->r_timeout);
    *n_read += x_read;
    if (peek  &&  x_read)  /* save the newly read data in the "peek-buffer" */
      VERIFY ( BUF_Write(&conn->buf, buf, x_read) );
  }}

  if (status != eCONN_Success) {
    if ( *n_read ) {
      s_ErrPost(SEV_INFO,  81, "[CONN_Read]  Error in reading data",
                x_connector, status);
      return eCONN_Success;
    } else {
      s_ErrPost((status == eCONN_Closed ? SEV_WARNING : SEV_ERROR), 82,
                "[CONN_Read]  Cannot read data", x_connector, status);
      return status;  /* error or EOF */
    }
  }

  /* success */
  return eCONN_Success;
}


/* Persistently read data from the input queue
 * See CONN_Read()
 */
static EConnStatus s_CONN_ReadPersist
(CONN       conn,
 void*      buf,
 Nlm_Uint4  size,
 Nlm_Uint4* n_read)
{
  ASSERT( *n_read == 0 );
  while (*n_read != size) {
    Nlm_Uint4 x_read;
    EConnStatus status = CONN_Read(conn, (char*)buf + *n_read,
                                   size - *n_read, &x_read, eCR_Read);
    *n_read += x_read;
    if (status != eCONN_Success)
      return status;
  }

  return eCONN_Success;
}


NLM_EXTERN EConnStatus CONN_Read
(CONN       conn,
 void*      buf,
 Nlm_Uint4  size,
 Nlm_Uint4* n_read,
 EConnRead  method)
{
  *n_read = 0;

  switch ( method ) {
  case eCR_Read:
    return s_CONN_Read(conn, buf, size, n_read, FALSE);
  case eCR_Peek:
    return s_CONN_Read(conn, buf, size, n_read, TRUE);
  case eCR_Persist:
    return s_CONN_ReadPersist(conn, buf, size, n_read);
  }
  return eCONN_Unknown;
}


NLM_EXTERN EConnStatus CONN_Close
(CONN conn)
{
  CONN_Reconnect(conn, 0);
  BUF_Destroy(conn->buf);
  Nlm_MemFree(conn);
  return eCONN_Success;
}
