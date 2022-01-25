/* $Id: ncbisock.c,v 6.1 1999/10/18 15:39:05 vakatov Exp $
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
 *   Plain portable TCP/IP socket API for:  UNIX, MS-Win, MacOS
 *
 *   This is just a back-compatibility interface("Nlm_*") to the real
 *   SOCK API located in "ncbi_socket.[ch]".
 *   Unlike the "real" SOCK API, this API uses:
 *    a) "Nlm_" name prefix for structures, types and functions;
 *    b) "Nlm_*" fixed size integer types like "Nlm_Uint4";
 *    c) "Nlm_Boolean" rather than a native "int" for the boolean type;
 *    d) [MSWIN] "NLM_EXTERN" rather than just "extern" to ease the compilation
 *       for MSWIN DLL.
 *
 * ---------------------------------------------------------------------------
 * $Log: ncbisock.c,v $
 * Revision 6.1  1999/10/18 15:39:05  vakatov
 * Initial revision
 * This is actually just an interface for the back compatibility with the
 * former "ncbisock.[ch]"; the real code is in "ncbi_socket.[ch]"
 *
 * ===========================================================================
 */

#include <ncbithr.h>
#include <ncbisock.h>

/* undefine all "Nlm_SOCK_*" to clear the access to <ncbi_socket.h> API */
#undef LSOCK
#undef SOCK

#undef ESOCK_ErrCode
#undef eSOCK_ESuccess
#undef eSOCK_ETimeout
#undef eSOCK_EClosed
#undef eSOCK_EUnknown


#undef ESOCK_Mode
#undef eSOCK_OnRead
#undef eSOCK_OnWrite
#undef eSOCK_OnReadWrite

#undef SOCK_ErrCodeStr

#undef SOCK_Initialize
#undef SOCK_Destroy

#undef LSOCK_Create
#undef LSOCK_Accept
#undef LSOCK_Close
#undef LSOCK_GetOSHandle

#undef SOCK_Create
#undef SOCK_SetTimeout
#undef SOCK_Select
#undef SOCK_Read
#undef SOCK_ReadPersist
#undef SOCK_Peek
#undef SOCK_PushBack
#undef SOCK_Eof
#undef SOCK_Write
#undef SOCK_Reconnect
#undef SOCK_Close
#undef SOCK_GetOSHandle

#undef GetHostName
#undef Uint4toInaddr

/* this is the only place where both <ncbibuf.h> and <ncbi_buffer.h> can
 * be #include'd in one source module! */
#undef NCBISOCK__H

#include <ncbi_socket.h>


/* ESOCK_Status <--> ESOCK_ErrCode
 */
#define S2E(status)   ((Nlm_ESOCK_ErrCode) status)
#define E2S(err_code) ((ESOCK_Status)    err_code)

/* ESOCK_Mode -> ESOCK_Direction
 */
#define M2D(mode) ((ESOCK_Direction) mode)

/* STimeout <--> SSOCK_Timeout
 */
static STimeout* s_ss2tt(const SSOCK_Timeout *ss, STimeout *tt)
{
  if ( !ss )
    return 0;

  tt->sec  = ss->sec;
  tt->usec = ss->usec % 1000000;
  return tt;
}

static SSOCK_Timeout* s_tt2ss(const STimeout *tt, SSOCK_Timeout *ss)
{
  if ( !tt )
    return 0;

  ss->sec  = tt->sec;
  ss->usec = tt->usec % 1000000;
  return ss;
}


/* MT-protection callback
 */
#if defined(__cplusplus)
extern "C" {
  static void s_MT_CSection(void* data, ESOCK_MT_Locking what);
}
#endif /* __cplusplus */
static void s_MT_CSection(void* data, ESOCK_MT_Locking what)
{
  TNlmMutex* mtx = (TNlmMutex*) data;

  switch ( what ) {
  case eSOCK_MT_Lock:
    NlmMutexLockEx(mtx);
    break;
  case eSOCK_MT_Unlock:
    NlmMutexUnlock(*mtx);
    break;
  case eSOCK_MT_Cleanup:
    NlmMutexDestroy(*mtx);
    *mtx = (TNlmMutex) ~0;
    break;
  }
}


/* Setup MT-protection callback
 */
static void s_MT_SetCriticalSection(void)
{
  if ( !NlmThreadsAvailable() ) {
    static TNlmMutex s_Mutex;
    NlmMutexLockEx(&s_Mutex);
    SOCK_MT_SetCriticalSection(s_MT_CSection, &s_Mutex);
    NlmMutexUnlock(s_Mutex);
  }
}



/***********************************************************************
 *  EXTERNAL
 ***********************************************************************/

NLM_EXTERN const char* Nlm_SOCK_ErrCodeStr
(Nlm_ESOCK_ErrCode err_code)
{
  return SOCK_StatusStr( E2S(err_code) );
}


NLM_EXTERN Nlm_ESOCK_ErrCode Nlm_SOCK_Destroy(void)
{
  return S2E( SOCK_ShutdownAPI() );
}


NLM_EXTERN Nlm_ESOCK_ErrCode Nlm_LSOCK_Create
(Nlm_Uint2  port,
 Nlm_Uint2  n_listen,
 LSOCK     *lsock)
{
  s_MT_SetCriticalSection();
  return S2E( LSOCK_Create(port, n_listen, lsock) );
}


NLM_EXTERN Nlm_ESOCK_ErrCode Nlm_LSOCK_Accept
(LSOCK           lsock,
 const STimeout *timeout,
 SOCK           *sock)
{
  SSOCK_Timeout x_timeout;
  return S2E( LSOCK_Accept(lsock, s_tt2ss(timeout, &x_timeout), sock) );
}


NLM_EXTERN Nlm_ESOCK_ErrCode Nlm_LSOCK_Close
(LSOCK lsock)
{
  return S2E( LSOCK_Close(lsock) );
}


NLM_EXTERN Nlm_ESOCK_ErrCode Nlm_SOCK_Create
(const char*     host,
 Nlm_Uint2       port,
 const STimeout* timeout,
 SOCK*           sock)
{
  SSOCK_Timeout x_timeout;
  s_MT_SetCriticalSection();
  return S2E( SOCK_Create(host, port, s_tt2ss(timeout, &x_timeout), sock) );
}


NLM_EXTERN Nlm_ESOCK_ErrCode Nlm_SOCK_Reconnect
(SOCK            sock,
 const char*     host,
 Nlm_Uint2       port,
 const STimeout* timeout)
{
  SSOCK_Timeout x_timeout;
  return S2E( SOCK_Reconnect(sock, host, port, s_tt2ss(timeout, &x_timeout)) );
}


NLM_EXTERN Nlm_ESOCK_ErrCode Nlm_SOCK_Close
(SOCK sock)
{
  return S2E( SOCK_Close(sock) );
}


NLM_EXTERN Nlm_ESOCK_ErrCode Nlm_SOCK_Select
(SOCK            sock,
 Nlm_ESOCK_Mode  mode,
 const STimeout* timeout)
{
  SSOCK_Timeout x_timeout;
  return S2E( SOCK_Wait(sock, M2D(mode), s_tt2ss(timeout, &x_timeout)) );
}


NLM_EXTERN Nlm_ESOCK_ErrCode Nlm_SOCK_SetTimeout
(SOCK            sock,
 Nlm_ESOCK_Mode  mode,
 const STimeout* new_timeout,
 STimeout*       r_timeout,
 STimeout*       w_timeout)
{
  /* retrieve R and/or W timeouts, if requested */
  static const STimeout s_Infinite = { 99999999, 999999 };
  if ( r_timeout ) {
    if ( !s_ss2tt(SOCK_GetReadTimeout(sock), r_timeout) )
      *r_timeout = s_Infinite;
  }
  if ( w_timeout ) {
    if ( !s_ss2tt(SOCK_GetWriteTimeout(sock), w_timeout) )
      *w_timeout = s_Infinite;
  }

  /* special case -- do not change the timeout(s) */
  if (new_timeout == SOCK_GET_TIMEOUT)
    return Nlm_eSOCK_Success;

  /* change R and/or W timeouts */
  {{
    SSOCK_Timeout x_timeout;
    return S2E( SOCK_SetTimeout(sock, M2D(mode),
                                s_tt2ss(new_timeout, &x_timeout)) );
  }}
}


NLM_EXTERN Nlm_ESOCK_ErrCode Nlm_SOCK_Read
(SOCK       sock,
 void*      buf,
 Nlm_Uint4  size,
 Nlm_Uint4* n_read)
{
  size_t x_read;
  Nlm_ESOCK_ErrCode err_code =
    S2E( SOCK_Read(sock, buf, (size_t) size, &x_read, eSOCK_Read) );
  *n_read = (Nlm_Uint4) x_read;
  return err_code;
}


NLM_EXTERN Nlm_ESOCK_ErrCode Nlm_SOCK_ReadPersist
(SOCK       sock,
 void*      buf,
 Nlm_Uint4  size,
 Nlm_Uint4* n_read)
{
  size_t x_read;
  Nlm_ESOCK_ErrCode err_code =
    S2E( SOCK_Read(sock, buf, (size_t) size, &x_read, eSOCK_Persist) );
  *n_read = (Nlm_Uint4) x_read;
  return err_code;
}


NLM_EXTERN Nlm_ESOCK_ErrCode Nlm_SOCK_Peek
(SOCK       sock,
 void*      buf,
 Nlm_Uint4  size,
 Nlm_Uint4* n_read)
{
  size_t x_read;
  Nlm_ESOCK_ErrCode err_code =
    S2E( SOCK_Read(sock, buf, (size_t) size, &x_read, eSOCK_Peek) );
  *n_read = (Nlm_Uint4) x_read;
  return err_code;
}


NLM_EXTERN Nlm_ESOCK_ErrCode Nlm_SOCK_PushBack
(SOCK        sock,
 const void* buf,
 Nlm_Uint4   size)
{
  return S2E( SOCK_PushBack(sock, buf, (size_t) size) );
}


NLM_EXTERN Nlm_Boolean Nlm_SOCK_Eof(SOCK sock)
{
    return SOCK_Eof(sock) ? TRUE : FALSE;
}


NLM_EXTERN Nlm_ESOCK_ErrCode Nlm_SOCK_Write
(SOCK        sock,
 const void* buf,
 Nlm_Uint4   size,
 Nlm_Uint4*  n_written)
{
  size_t x_written;
  Nlm_ESOCK_ErrCode err_code =
    S2E( SOCK_Write(sock, buf, (size_t) size, &x_written) );
  if ( n_written )
    *n_written = (Nlm_Uint4) x_written;
  return err_code;
}


NLM_EXTERN void Nlm_SOCK_Address
(SOCK         sock,
 Nlm_Uint4*   host,
 Nlm_Uint2*   port,
 Nlm_Boolean  network_byte_order)
{
  if ( host ) {
    unsigned int x_host;
    SOCK_GetAddress(sock, &x_host, 0, network_byte_order ? 1 : 0);
    *host = x_host;
  }
  if ( port ) {
    unsigned short x_port;
    SOCK_GetAddress(sock, 0, &x_port, network_byte_order ? 1 : 0);
    *port = x_port;
  }
}


NLM_EXTERN Nlm_Boolean Nlm_GetHostName
(char*     name,
 Nlm_Uint4 namelen)
{
  return SOCK_gethostname(name, (size_t) namelen) ? FALSE : TRUE;
}


NLM_EXTERN Nlm_Boolean Nlm_Uint4toInaddr
(Nlm_Uint4 ui4_addr,
 char*     buf,
 Nlm_Uint4 buf_len)
{
  s_MT_SetCriticalSection();
  return SOCK_host2inaddr((unsigned int) ui4_addr, buf, (size_t) buf_len) ?
    TRUE : FALSE;
}

NLM_EXTERN Nlm_Uint4 Nlm_htonl(Nlm_Uint4 value)
{
  return (Nlm_Uint4) SOCK_htonl((Nlm_Uint4) value);
}
