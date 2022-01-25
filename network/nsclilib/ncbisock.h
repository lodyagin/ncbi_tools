#ifndef NCBISOCK__H
#define NCBISOCK__H

/*  $RCSfile: ncbisock.h,v $  $Revision: 4.5 $  $Date: 1998/08/12 13:12:43 $
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
*   Plain portable socket API
*
* --------------------------------------------------------------------------
* $Log: ncbisock.h,v $
* Revision 4.5  1998/08/12 13:12:43  kans
* moved high level query functions to ncbiurl.[ch]
*
* Revision 4.4  1998/08/10 23:24:25  kans
* added SOCK_SendURLQuery and SOCK_CheckURLQuery from Sequin
*
* Revision 4.3  1998/03/30 17:50:13  vakatov
* Ingrafted to the main NCBI CVS tree
*
* ==========================================================================
*/

#include <ncbistd.h>

/* for ErrPost */
#define SOCK_ERRCODE  777

#define LSOCK              Nlm_LSOCK
#define SOCK               Nlm_SOCK

#define ESOCK_ErrCode      Nlm_ESOCK_ErrCode
#define ESOCK_Mode         Nlm_ESOCK_Mode
#define STimeout           Nlm_STimeout

#define SOCK_ErrCodeStr    Nlm_SOCK_ErrCodeStr

#define SOCK_Initialize    Nlm_SOCK_Initialize
#define SOCK_Destroy       Nlm_SOCK_Destroy

#define LSOCK_Create       Nlm_LSOCK_Create
#define LSOCK_Accept       Nlm_LSOCK_Accept
#define LSOCK_Close        Nlm_LSOCK_Close

#define SOCK_Create        Nlm_SOCK_Create
#define SOCK_SetTimeout    Nlm_SOCK_SetTimeout
#define SOCK_Select        Nlm_SOCK_Select
#define SOCK_Read          Nlm_SOCK_Read
#define SOCK_Write         Nlm_SOCK_Write
#define SOCK_Close         Nlm_SOCK_Close


#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

#ifdef __cplusplus
extern "C" {
#endif

struct Nlm_LSOCKtag;                /* listening socket:  internal storage  */
typedef struct Nlm_LSOCKtag *LSOCK; /* listening socket:  handle */

struct Nlm_SOCKtag;               /* socket:  internal storage  */
typedef struct Nlm_SOCKtag *SOCK; /* socket:  handle */

typedef enum
{
  eSOCK_ESuccess = 0, /* everything is fine, no errors occured          */
  eSOCK_ETimeout,     /* timeout expired before the data could be i/o'd */
  eSOCK_EClosed,      /* peer has closed the connection                 */

  eSOCK_EUnknown      /* unknown(most probably -- fatal) error          */
} ESOCK_ErrCode;


typedef enum {
  eSOCK_OnRead,
  eSOCK_OnWrite,
  eSOCK_OnReadWrite
} ESOCK_Mode;

typedef struct {
 Nlm_Uint4 sec;   /* seconds(get truncated to platform-dep. max. limit) */
 Nlm_Uint4 usec;  /* microseconds(always get truncated to 10**6)        */
} STimeout;


/* Return (const) verbal description for the passed error code
 */
NLM_EXTERN const Nlm_Char *SOCK_ErrCodeStr
(ESOCK_ErrCode err_code
 );


/* [SERVER-side]  Create and initialize the server-side(listening) socket
 * (socket() + bind() + listen())
 */
NLM_EXTERN ESOCK_ErrCode LSOCK_Create
(Nlm_Uint2  port,     /* [in] the port to listen at            */
 Nlm_Uint2  n_listen, /* [in] maximal # of pending connections */
 LSOCK     *lsock     /* [out]  handle of the created listening socket  */
 );


/* [SERVER-side]  Accept connection from a client
 * NOTE: the "*timeout" is for this accept() only;  to set i/o timeout,
 *       use SOCK_Timeout(); (by default -- infinite)
 */
NLM_EXTERN ESOCK_ErrCode LSOCK_Accept
(LSOCK           lsock,    /* [in] handle of a listening socket   */
 const STimeout *timeout,  /* [in] timeout(infinite if NULL)      */
 SOCK           *sock      /* [out]  handle of the created socket */
 );


/* [SERVER-side]  Close the listening socket, destroy relevant internal data
 */
NLM_EXTERN ESOCK_ErrCode LSOCK_Close
(LSOCK lsock
 );


/* [CLIENT-side]  Connect client to another(server-side, listening) socket
 * (socket() + connect())
 */
NLM_EXTERN ESOCK_ErrCode SOCK_Create
(const Nlm_Char *host, /* [in] server host */
 Nlm_Uint2       port, /* [in] server port */
 SOCK           *sock  /* [out]  handle of the created socket */
 );


/* [CLIENT-side]  Close the connection, destroy relevant internal data
 * NOTE: if write timeout is specified then it blocks until either all
 *       unsent data are sent or until the timeout expires
 * NOTE: whatever error code is returned, this function cannot be
 *       called more than once for the same socket
 */
NLM_EXTERN ESOCK_ErrCode SOCK_Close
(SOCK sock
 );


/* Block on the socket until either read/write(dep. on "mode") is
 * available or timeout expires(if "timeout" is NULL then assume it infinite)
 */
NLM_EXTERN ESOCK_ErrCode SOCK_Select
(SOCK            sock,
 ESOCK_Mode      mode,
 const STimeout *timeout
 );


/* Specify timeout for the connection i/o(see SOCK_[Read|Write|Close] funcs).
 * NOTE: set the timeout to the maximum if "new_timeout" is NULL
 * NOTE: the default timeout is the maximum possible(wait "ad infinitum")
 */
#define SOCK_GET_TIMEOUT ((const STimeout *)~0)
NLM_EXTERN ESOCK_ErrCode SOCK_SetTimeout
(SOCK           sock,
 ESOCK_Mode     mode,
 const STimeout *new_timeout, /* (dont set if equal to SOCK_GET_TIMEOUT) */
 STimeout       *r_timeout,   /* if non-NULL, return previous read */
 STimeout       *w_timeout    /* and(or) write timeout values      */
 );


/* Read up to "size" bytes from "sock" to the mem.buffer pointed by "buf".
 * In "*n_read", return the number of succesfully read bytes.
 * If there is no data available to read and the timeout(see
 * SOCK_Timeout()) is expired then return eSOCK_ETimeout.
 * NOTE: eSOCK_Closed may indicate an empty message rather than a
 *       a real closure of connection
 */
NLM_EXTERN ESOCK_ErrCode SOCK_Read
(SOCK        sock,
 Nlm_VoidPtr buf,
 Nlm_Uint4   size,
 Nlm_Uint4  *n_read
 );

/* Operate just like SOCK_Read() but dont remove the read data from the
 * input queue
 */
NLM_EXTERN ESOCK_ErrCode SOCK_Peek
(SOCK        sock,
 Nlm_VoidPtr buf,
 Nlm_Uint4   size,
 Nlm_Uint4  *n_read
 );


/* Write "size" bytes from the mem.buffer "buf" to "sock".
 * In "*n_written", return the number of successfully written bytes.
 * If cannot write all data and the timeout(see SOCK_Timeout()) is expired
 * then return eSOCK_ETimeout.
 */
NLM_EXTERN ESOCK_ErrCode SOCK_Write
(SOCK        sock,
 const void *buf,
 Nlm_Uint4   size,
 Nlm_Uint4  *n_written
 );

/* Destroy internal data used by this module
 * NOTE: no function from this API can be used after the call to SOCK_Destroy
 */
NLM_EXTERN ESOCK_ErrCode SOCK_Destroy(void);


/* To be moved from this module to... somewhere else; later... maybe ;-)
 */
#define GetHostName Nlm_GetHostName
NLM_EXTERN Nlm_Boolean GetHostName
(Nlm_Char *name,
 Nlm_Uint4 namelen
 );

#define Uint4toInaddr Nlm_Uint4toInaddr
NLM_EXTERN Nlm_Boolean Uint4toInaddr
(Nlm_Uint4   ui4_addr,  /* NOTE: must be in the local host byte order  */
 Nlm_CharPtr buf,       /* to be filled by smth. like "123.45.67.89\0" */
 Nlm_Uint4   buf_len
 );

/* KLUDGE(dont use this beast, please) */
extern Nlm_Uint4 Nlm_htonl(Nlm_Uint4 value);


#ifdef __cplusplus
}
#endif

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif /* NCBISOCK__H */
