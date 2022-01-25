#ifndef CONNECTN__H
#define CONNECTN__H

/*  $Id: connectn.h,v 6.2 1999/04/02 20:41:53 vakatov Exp $
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
*   Several methods can be used to establish the connection, and each of them
*   yields in a simple handle(of type "CONN") that contains a handle(of type
*   "CONNECTOR") to a data and methods implementing generic connection i/o
*   operations. E.g. this API can be used to:
*     1) connect using HTTPD-based dispatcher(to NCBI services only!);
*     2) hit a CGI script;
*     3) connect to a bare socket at some "host:port";
*     4) whatever else can be fit into this paradigm -- see the
*        SConnectorTag-related structures;  e.g. it could be a plain file i/o
*        or even a memory area.
*
*  See in "connectr.h" for the detailed specification of the underlying
*  connector("CONNECTOR", "SConnectorTag") methods and structures.
*
* --------------------------------------------------------------------------
* $Log: connectn.h,v $
* Revision 6.2  1999/04/02 20:41:53  vakatov
* Got rid of occasional comment inside comment(in CVS Log)
*
* Revision 6.1  1999/04/01 21:48:10  vakatov
* Fixed for the change in spec:  "n_written/n_read" args in
* CONN_Write/Read to be non-NULL and "*n_written / *n_read" := 0
*
* Revision 6.0  1999/03/25 23:04:57  vakatov
* Initial revision
* ==========================================================================
*/


/* for ErrPost */
#define CONN_ERRCODE  778

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

#ifdef __cplusplus
extern "C" {
#endif

struct SConnectionTag;
typedef struct SConnectionTag* CONN;      /* connection handle */

struct SConnectorTag;
typedef struct SConnectorTag* CONNECTOR;  /* connector handle */


/* Status of the last operation on the connection
 */
typedef enum {
  eCONN_Success = 0,  /* everything is fine, no errors occured          */
  eCONN_Timeout,      /* timeout expired before the data could be i/o'd */
  eCONN_Closed,       /* peer has closed the connection                 */
  eCONN_InvalidArg,   /* bad argument value                             */
  eCONN_NotSupported, /* the requested operation is not supported       */

  eCONN_Unknown       /* unknown(most probably -- fatal) error          */
} EConnStatus;

/* Return verbal description for the passed status code
 */
NLM_EXTERN const Nlm_Char* CONN_StatusString(EConnStatus status);



/* Connection i/o direction
 */
typedef enum {
  eCONN_Connect,   /* on connect */
  eCONN_Read,      /* on read */
  eCONN_Write,     /* on write */
  eCONN_ReadWrite, /* on both read and write */
  eCONN_Close      /* on close */
} EConnDirection;



/* Compose all data necessary to establish a new connection
 * (merely bind it to the specified connector).
 * NOTE:  The real connection will not be established right away. Instead,
 *        it will be established at the moment of the first call to one of
 *        "Wait", "WaitAsync", "Write" or "Read".
 */
NLM_EXTERN EConnStatus CONN_Create
(CONNECTOR connector, /* [in]  connector */
 CONN*     conn       /* [out] handle of the created connection */
 );


/* Reconnect, using "connector".
 * If "conn" is already opened then close the current connection first,
 * even if "connector" is just the same as the current connector.
 * If "connector" is NULL then use the current connector to (re)connect.
 * NOTE:  Alhough it closes the previos connection immediately, however it
 *        does not establish the new connection right away -- it postpones
 *        until the first call to "Wait", "WaitAsync", "Write" or "Read".
 */
NLM_EXTERN EConnStatus CONN_Reconnect
(CONN      conn,     /* [in] connection handle */
 CONNECTOR connector /* [in] new connector */
 );


/* Specify timeout for the connection i/o (including "CONNECT" and "CLOSE").
 * This function can be called at any time during the connection lifetime.
 * NOTE: if "new_timeout" is NULL then set the timeout to the maximum.
 * NOTE: the default timeout is the maximum possible (wait "ad infinitum").
 */
NLM_EXTERN void CONN_SetTimeout
(CONN            conn,        /* [in]  connection handle */
 EConnDirection  direction,   /* [in]  i/o direction */
 const STimeout* new_timeout  /* [in]  new timeout */
 );


/* Retrieve current timeout(return NULL if indefinite).
 * The returned pointer is guaranteed to point to the valid timeout structure
 * (or NULL) until the next CONN_SetTimeout or CONN_Close function call.
 */
NLM_EXTERN const STimeout* CONN_GetTimeout
(CONN           conn,      /* [in]  connection handle */
 EConnDirection direction  /* [in]  i/o direction, not "eCONN_ReadWrite" */
 );


/* Block on the connection until it becomes available for either read or
 * write(dep. on "direction"), or until the timeout expires(if "timeout"
 * is NULL then assume it infinite), or until any error.
 */
NLM_EXTERN EConnStatus CONN_Wait
(CONN            conn,      /* [in] connection handle */
 EConnDirection  direction, /* [in] i/o, can only be eCONN_Read|Write" */
 const STimeout* timeout    /* [in] the blocking timeout */
 );


#ifdef IMPLEMENTED__CONN_WaitAsync
/* Wait for an asynchronous i/o event, then call the specified handler.
 * In the "handler" function:
 *   "direction" - is the i/o direction where the async. event happened
 *   "status"    - must be "eCONN_Success" if it is ready for i/o
 *   "data"      - callback data(passed as "data" in CONN_WaitAsync)
 * If "handler" is NULL then discard the current handler, if any.
 * The "cleanup" function to be called right after the call to "handler" or
 * by CONN_Close(), or if the handler is reset by calling CONN_WaitAsync()
 * again -- whichever happens first.
 */
typedef void (*FConnAsyncHandler)
(CONN           conn,
 EConnDirection direction,
 EConnStatus    status,
 void*          data
);
typedef void (*FConnAsyncCleanup)(void* data);

NLM_EXTERN EConnStatus CONN_WaitAsync
(CONN              conn,      /* [in] connection handle */
 EConnDirection    direction, /* [in] i/o direction */
 FConnAsyncHandler handler,   /* [in] callback function */
 void*             data,      /* [in] callback data */
 FConnAsyncCleanup cleanup    /* [in] cleanup procedure */
 );
#endif /* IMPLEMENTED__CONN_WaitAsync */


/* Write "size" bytes from the mem.buffer "buf" to "conn".
 * In "*n_written", return the number of successfully written bytes.
 * If cannot write all data and the timeout(see CONN_Timeout) is expired
 * then return eCONN_ETimeout.
 */
NLM_EXTERN EConnStatus CONN_Write
(CONN        conn,      /* [in]  connection handle */ 
 const void* buf,       /* [in]  pointer to the data buffer to write */ 
 Nlm_Uint4   size,      /* [in]  # of bytes to write */ 
 Nlm_Uint4*  n_written  /* [out] # of actually written bytes(not NULL!) */
 );


/* Explicitly flush(send to the connection) the data written by "CONN_Write()".
 * NOTE:  "CONN_Close|Peek|Read()" always call "CONN_Flush()" before reading.
 */
NLM_EXTERN EConnStatus CONN_Flush
(CONN conn  /* [in] connection handle */ 
 );


/* Read up to "size" bytes from "conn" to the mem.buffer pointed by "buf".
 * In "*n_read", return the number of succesfully read bytes.
 * If there is no data available to read and the timeout(see
 * CONN_Timeout()) is expired then return eCONN_ETimeout(and "*n_read" := 0).
 * NOTE:  CONN_Read(eCONN_Peek...) uses read(not peek!) function provided
 *        by the connector, and the "peeked" data gets stored in the
 *        connection internals.
 */
typedef enum {
  eCR_Read,   /* read presently available data */
  eCR_Peek,   /* eCONN_Read, but dont discard the read data from inp. queue */
  eCR_Persist /* try to read/peek exactly "size" bytes -- read again */
} EConnRead;

NLM_EXTERN EConnStatus CONN_Read
(CONN        conn,   /* [in]  connection handle */
 void*       buf,    /* [out] memory buffer to read to */
 Nlm_Uint4   size,   /* [in]  max. # of bytes to read */
 Nlm_Uint4*  n_read, /* [out] # of actually read bytes(not NULL!) */
 EConnRead   method  /* [in]  read/peek | persist */
 );


/* Close the connection, destroy relevant internal data.
 * NOTE: whatever error code is returned, this function cannot be
 *       called more than once for the same connection.
 */
NLM_EXTERN EConnStatus CONN_Close
(CONN conn  /* [in] connection handle */
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

#endif /* CONNECTN__H */
