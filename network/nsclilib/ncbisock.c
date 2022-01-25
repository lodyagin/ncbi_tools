/*  $RCSfile: ncbisock.c,v $  $Revision: 4.12 $  $Date: 1999/01/22 22:04:59 $
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
* $Log: ncbisock.c,v $
* Revision 4.12  1999/01/22 22:04:59  vakatov
* Uint4toInaddr() to take address in the network byte order
*
* Revision 4.11  1998/12/15 17:22:22  vakatov
* + SOCK_Address() -- to get the socket peer's host and port
*
* Revision 4.10  1998/08/12 13:12:42  kans
* moved high level query functions to ncbiurl.[ch]
*
* Revision 4.9  1998/08/10 23:24:24  kans
* added SOCK_SendURLQuery and SOCK_CheckURLQuery from Sequin
*
* Revision 4.8  1998/05/01 19:13:55  vakatov
* [OS_UNIX]  Ignore SIGPIPE signal
*
* Revision 4.7  1998/04/14 20:04:25  vakatov
* Added a "hand-made" data buffering on the socket read(#SOCK_NCBI_PEEK)
* in order to implement MSG_PEEK functionality
* [OS_MAC]  #define SOCK_NCBI_PEEK as MSG_PEEK was not implemented on Mac
* [OS_UNIX] LSOCK_Create():  use SO_REUSEADDR(see the comments)
*
* Revision 4.6  1998/04/07 17:39:13  kans
* comment out mac definition of TEST_MODULE__NCBISOCK and DO_CLIENT
*
* Revision 4.5  1998/03/30 17:50:11  vakatov
* Ingrafted to the main NCBI CVS tree
*
* ==========================================================================
*/

#include <ncbi.h>
#include <ncbithr.h>

#if defined(OS_UNIX)
#include <ncbiwin.h>
#include <sys/socket.h>
#include <unistd.h>
#include <netdb.h>
#include <netinet/in.h>
#include <fcntl.h>
#include <arpa/inet.h>
#include <signal.h>

#elif defined(OS_MSWIN)
#include <ncbiwin.h>
#include <winsock.h>

#elif defined(OS_MAC)
#include "ncbinet.h"
#include "ni_net.h"
#include <ncbiwin.h>
/* cannot write more than SLICE_SIZE at once on Mac, and we have to split
 * big output buffers into smaller(<SLICE_SIZE) slices before writing it to
 * socket */
#define SOCK_WRITE_SLICED
#define SLICE_SIZE 2048
/* MSG_PEEK is not implemented(ignored) on Mac, and we have to implement
 * this feature ourselves */
#define SOCK_NCBI_PEEK

#endif /* platform-specific #include's (UNIX, MSWIN, MAC) */


#include <ncbisock.h>

#ifdef SOCK_NCBI_PEEK
#include <ncbibuf.h>
#endif


/* patch to old SunOS/Solaris proto(cut&paste from Solaris 2.6 "unistd.h") */
#ifdef __cplusplus
extern "C" {
#endif

#if defined(OS_UNIX_SOL) || defined(OS_UNIX_SUN)
#if defined(_XPG4_2)
extern int gethostname(char *, size_t);
#elif  defined(__EXTENSIONS__) || \
        (!defined(_POSIX_C_SOURCE) && !defined(_XOPEN_SOURCE))
extern int gethostname(char *, int);
#endif
#endif

#ifdef __cplusplus
}
#endif



/***********************************************************************
 *  INTERNAL
 ***********************************************************************/

#define SOCK_MAX_TIMEOUT 99999999

#ifdef SOCK_NCBI_PEEK
#define SOCK_RECV(s,b,l,f) s_NCBI_Recv(s, (char *)b, l, f)
#else
#define SOCK_RECV(s,b,l,f) recv(s->sock, (char *)b, l, f)
#endif /* SOCK_NCBI_PEEK */

#define SOCK_WRITE(s,b,l)  send(s, (char *)b, l, 0)


#if defined(OS_MSWIN)

typedef SOCKET Nlm_Socket;
#define SOCK_INVALID     INVALID_SOCKET
#define SOCK_ERRNO       WSAGetLastError()
#define SOCK_EINTR       WSAEINTR
#define SOCK_EWOULDBLOCK WSAEWOULDBLOCK
#define SOCK_ECONNRESET  WSAECONNRESET
#define SOCK_EPIPE       WSAESHUTDOWN
#define SOCK_EAGAIN      WSAEINPROGRESS
#define SOCK_NFDS(s)     0
#define SOCK_CLOSE(s)    closesocket(s)

#ifdef WIN16
#ifndef MAKEWORD
#define MAKEWORD(a,b)  ((WORD)(((BYTE)(a)) | ((WORD)((BYTE)(b))) <<8)) 
#endif 
#endif
/* OS_MSWIN */

#elif defined(OS_UNIX)

typedef int Nlm_Socket;
#define SOCK_INVALID     (-1)
#define SOCK_ERRNO       errno
#define SOCK_EINTR       EINTR
#define SOCK_EWOULDBLOCK EWOULDBLOCK
#define SOCK_ECONNRESET  ECONNRESET
#define SOCK_EPIPE       EPIPE
#define SOCK_EAGAIN      EAGAIN
#define SOCK_NFDS(s)     (s + 1)
#define SOCK_CLOSE(s)    close(s)

#ifndef INADDR_NONE
#define INADDR_NONE      (Nlm_Uint4)(-1)
#endif
/* OS_UNIX */

#elif defined(OS_MAC)

typedef int Nlm_Socket;
#define SOCK_INVALID     (-1)
#ifndef SOCK_ERRNO
#define SOCK_ERRNO       errno
#endif
#define SOCK_EINTR       EINTR
#define SOCK_EWOULDBLOCK EWOULDBLOCK
#define SOCK_ECONNRESET  ECONNRESET
#define SOCK_EPIPE       EPIPE
#define SOCK_EAGAIN      EAGAIN
#define SOCK_NFDS(s)     (s + 1)
#define SOCK_CLOSE(s)    close(s)

extern int gethostname(char *machname, long buflen);
/* but see ni_lib.c line 2508 for gethostname substitute for Mac */

#endif /* !OS_MSWIN!OS_UNIX!OS_MAC */


/* Listening socket */
typedef struct Nlm_LSOCKtag
{
  Nlm_Socket sock;
} Nlm_LSOCKstruct;

/* Socket */
typedef struct Nlm_SOCKtag
{
  Nlm_Socket sock;
  struct timeval *r_timeout;
  struct timeval *w_timeout;
  struct timeval rr_timeout;
  struct timeval ww_timeout;
  Nlm_Uint4 host;  /* peer host (in the network byte order) */
  Nlm_Uint2 port;  /* peer port (in the network byte order) */
#ifdef SOCK_NCBI_PEEK
  BUF buf;
#endif
} Nlm_SOCKstruct;


/* STimeout <--> struct timeval  conversions
 */
static void s_tv2to(const struct timeval *tv, STimeout *to) {
  if ( tv ) {
    to->sec  = tv->tv_sec;
    to->usec = tv->tv_usec;
  }
  else {
    to->sec  = SOCK_MAX_TIMEOUT;
    to->usec = 999999;
  }
}

static void s_to2tv(const STimeout *to, struct timeval *tv) {
  if ( to ) {
    tv->tv_sec  = (to->sec < SOCK_MAX_TIMEOUT) ? to->sec : SOCK_MAX_TIMEOUT;
    tv->tv_usec = to->usec % 1000000;
  } else {
    tv->tv_sec  = SOCK_MAX_TIMEOUT;
    tv->tv_usec = 999999;
  }
}


/* Initialize the API internal data and/or underlying network code
 */
#if defined(OS_MSWIN)
static TNlmMutex   s_InitMutex;
static Nlm_Boolean s_Initialized = FALSE;
#endif

static ESOCK_ErrCode s_Initialize(void)
{
  ESOCK_ErrCode err_code = eSOCK_ESuccess;

#if defined(OS_MSWIN)
  if ( s_Initialized )
    return eSOCK_ESuccess;

  NlmMutexLockEx(&s_InitMutex);
  if ( !s_Initialized ) {
    WSADATA wsadata;
    int x_errno = WSAStartup(MAKEWORD(1,1), &wsadata);
    if (x_errno != 0) {
      ErrPostEx(SEV_ERROR, SOCK_ERRCODE, 0,
                "[s_Initialize] WSAStartup() failed, errno = %d",(int)x_errno);
      err_code = eSOCK_EUnknown;
    }
    s_Initialized = TRUE;
  }
  NlmMutexUnlock(s_InitMutex);

#elif defined(OS_UNIX)
  signal(SIGPIPE, SIG_IGN);
#endif

  return err_code;
}


/* Switch the specified socket i/o between blocking and non-blocking mode
 */
static Nlm_Boolean s_SetNonblock(Nlm_Socket sock, Nlm_Boolean nonblock)
{
#if defined(OS_MSWIN)
  unsigned long argp = nonblock ? 1 : 0;
  return (Nlm_Boolean)(ioctlsocket(sock, FIONBIO, &argp) == 0);
#elif defined(OS_UNIX)
  return (Nlm_Boolean)(fcntl(sock, F_SETFL,
                             nonblock ?
                             fcntl(sock, F_GETFL, 0) | O_NONBLOCK :
                             fcntl(sock, F_GETFL, 0) & (int)~O_NONBLOCK)
                       != -1);
#elif defined(OS_MAC)
  return (Nlm_Boolean)(fcntl(sock, F_SETFL,
                             nonblock ?
                             fcntl(sock, F_GETFL, 0) | O_NDELAY :
                             fcntl(sock, F_GETFL, 0) & (int)~O_NDELAY)
                       != -1);
#else
  return FALSE;
#endif /* else!OS_MSWIN!OS_UNIX */
}


/* Select on the socket i/o
 */
static ESOCK_ErrCode s_Select(Nlm_Socket            sock,
                              ESOCK_Mode            mode,
                              const struct timeval *timeout)
{
  int n_dfs;
  fd_set fds, *r_fds, *w_fds;

  /* setup i/o descriptor to select */
  r_fds = (mode == eSOCK_OnRead  ||  mode == eSOCK_OnReadWrite) ? &fds : 0;
  w_fds = (mode == eSOCK_OnWrite ||  mode == eSOCK_OnReadWrite) ? &fds : 0;

  do { /* auto-resume if interrupted by a signal */
    fd_set e_fds;
    struct timeval tmout;
    if ( timeout )
      tmout = *timeout;
    FD_ZERO(&fds);       FD_ZERO(&e_fds);
    FD_SET(sock, &fds);  FD_SET(sock, &e_fds);
    n_dfs = select(SOCK_NFDS(sock), r_fds, w_fds, &e_fds,
                   timeout ? &tmout : 0);
    ASSERT ( -1 <= n_dfs  &&  n_dfs <= 2 );
    if ((n_dfs < 0  &&  SOCK_ERRNO != SOCK_EINTR)  ||  FD_ISSET(sock, &e_fds))
      return eSOCK_EUnknown;
  } while (n_dfs < 0);

  /* timeout has expired */
  if (n_dfs == 0)
    return eSOCK_ETimeout;

  /* success;  can i/o now */
  ASSERT ( FD_ISSET(sock, &fds) );
  return eSOCK_ESuccess;
}


#ifdef SOCK_NCBI_PEEK
/* Emulate "peek" using the NCBI data buffering
 */
static int s_NCBI_Recv(SOCK        sock,
                       Nlm_CharPtr buffer,
                       Nlm_Uint4   size,
                       int         flags)
{
  Nlm_Uint4 n_readbuf;
  int       n_readsock;

  /* read(or peek)from the internal buffer */
  n_readbuf = (flags & MSG_PEEK) ?
    BUF_Peek(sock->buf, buffer, size) : BUF_Read(sock->buf, buffer, size);
  if (n_readbuf == size)
    return (int)n_readbuf;
  buffer += n_readbuf;
  size   -= n_readbuf;

  /* read(dont peek) from the socket */
  n_readsock = recv(sock->sock, buffer, size, 0);
  if (n_readsock <= 0)
    return (int)(n_readbuf ? n_readbuf : n_readsock);

  /* if "peek" -- store the new read data in the internal buffer */
  if (flags & MSG_PEEK)
    VERIFY ( BUF_Write(&sock->buf, buffer, n_readsock) );

  return (int)(n_readbuf + n_readsock);
}
#endif /* SOCK_NCBI_PEEK */


/* Read/Peek data from the socket
 */
static ESOCK_ErrCode s_Recv(SOCK         sock,
                            Nlm_VoidPtr  buf,
                            Nlm_Uint4    size,
                            Nlm_Uint4   *n_read,
                            Nlm_Boolean  peek)
{
  int x_errno;

  if ( n_read )
    *n_read = 0;
  
  for (;;) {
    /* try to read */
    int buf_read = SOCK_RECV(sock, buf, (int)size, peek ? MSG_PEEK : 0);
    if (buf_read > 0) {
      if ( n_read )
        *n_read = buf_read;
      return eSOCK_ESuccess; /* success */
    }

    x_errno = SOCK_ERRNO;

    if (buf_read == 0) {
      /* NOTE:  empty message may cause the same effect, and it does */
      /* not (!) get discarded from the input queue; therefore, the  */
      /* subsequent attemts to read will cause just the same effect) */  
      return eSOCK_EClosed;  /* the connection must have been closed by peer */
    }

    /* wait for a data to come(or exit on timeout/error) */
    if (x_errno == SOCK_EWOULDBLOCK  ||  x_errno == SOCK_EAGAIN) {
      ESOCK_ErrCode err_code = s_Select(sock->sock,
                                        eSOCK_OnRead, sock->r_timeout);
      if (err_code != eSOCK_ESuccess)
        return err_code;
      continue;
    }

    /* retry if interrupted by a signal */
    if (x_errno == SOCK_EINTR)
      continue;

    /* forcibly closed by peer, or shut down */
    if (x_errno == SOCK_ECONNRESET  ||  x_errno == SOCK_EPIPE)
      return eSOCK_EClosed;

    /* dont want to handle all possible errors... let it be "unknown" */  
    break;
  }
  return eSOCK_EUnknown;
}


/* Write data to the socket "as is"(the whole buffer at once)
 */
static ESOCK_ErrCode s_WriteWhole(SOCK        sock,
                                  const void *buf,
                                  Nlm_Uint4   size,
                                  Nlm_Uint4  *n_written)
{
  const Nlm_Char *x_buf  = (const Nlm_Char *)buf;
  int             x_size = (int)size;
  int             x_errno;

  if ( n_written )
    *n_written = 0;
  if ( !size )
    return eSOCK_ESuccess;

  for (;;) {
    /* try to write */
    int buf_written = SOCK_WRITE(sock->sock, x_buf, x_size);
    if (buf_written >= 0) {
      if ( n_written )
        *n_written += buf_written;
      if (buf_written == x_size)
        return eSOCK_ESuccess; /* all data has been successfully sent */
      x_buf  += buf_written;
      x_size -= buf_written;
      ASSERT ( x_size > 0 );
      continue; /* there is unsent data */
    }
    x_errno = SOCK_ERRNO;

    /* blocked;  retry if unblocked before the timeout is expired */
    if (x_errno == SOCK_EWOULDBLOCK  ||  x_errno == SOCK_EAGAIN) {
      ESOCK_ErrCode err_code = s_Select(sock->sock,
                                        eSOCK_OnWrite, sock->w_timeout);
      if (err_code != eSOCK_ESuccess)
        return err_code;
      continue;
    }

    /* retry if interrupted by a signal */
    if (x_errno == SOCK_EINTR)
      continue;

    /* forcibly closed by peer, or shut down */
    if (x_errno == SOCK_ECONNRESET  ||  x_errno == SOCK_EPIPE)
      return eSOCK_EClosed;

    /* dont want to handle all possible errors... let it be "unknown" */  
    break;
  }
  return eSOCK_EUnknown;
}


#ifdef SOCK_WRITE_SLICED
/* Split output buffer by slices(of size <= SLICE_SIZE) before writing
 * to the socket
 */
static ESOCK_ErrCode s_WriteSliced(SOCK        sock,
                                   const void *buf,
                                   Nlm_Uint4   size,
                                   Nlm_Uint4  *n_written)
{
  ESOCK_ErrCode err_code  = eSOCK_ESuccess;
  Nlm_Uint4     x_written = 0;

  while (size  &&  err_code == eSOCK_ESuccess) {
    Nlm_Uint4 n_io = (size > SLICE_SIZE) ? SLICE_SIZE : size;
    Nlm_Uint4 n_io_done;
    err_code = s_WriteWhole(sock, (char*)buf + x_written, n_io, &n_io_done);
    if ( n_io_done ) {
      x_written += n_io_done;
      size      -= n_io_done;
    }
  }

  if ( n_written )
    *n_written = x_written;
  return err_code;
}
#endif


/***********************************************************************
 *  EXTERNAL
 ***********************************************************************/

NLM_EXTERN const Nlm_Char *SOCK_ErrCodeStr(ESOCK_ErrCode err_code)
{
  static const Nlm_Char *s_ErrCodeStr[eSOCK_EUnknown+1] = {
    "Success",
    "Timeout",
    "Closed",
    "Unknown"
  };

  return s_ErrCodeStr[err_code];
}


NLM_EXTERN ESOCK_ErrCode SOCK_Destroy(void)
{
#if defined(OS_MSWIN)
  if ( s_Initialized ) {
    if (WSACleanup() != 0) {
      ErrPostEx(SEV_WARNING, SOCK_ERRCODE, 0,
                "[SOCK_Destroy] WSACleanup() failed, errno = %d",
                (int)SOCK_ERRNO);
      return eSOCK_EUnknown;
    }

    NlmMutexDestroy(s_InitMutex);
    s_InitMutex = 0;
    s_Initialized = FALSE;
  }
#endif

  return eSOCK_ESuccess;
}


NLM_EXTERN ESOCK_ErrCode LSOCK_Create(Nlm_Uint2  port,
                                      Nlm_Uint2  n_listen,
                                      LSOCK     *lsock)
{
  Nlm_Socket sock;
  struct sockaddr_in addr;

  /* Initialize internals */
  VERIFY ( s_Initialize() == eSOCK_ESuccess );

  /* Create new(listening) socket */
  if ((sock = socket(AF_INET, SOCK_STREAM, 0)) == SOCK_INVALID) {
    ErrPostEx(SEV_ERROR, SOCK_ERRCODE, 1,
              "[LSOCK_Create]  Cannot create listening socket, errno = %d",
              (int)SOCK_ERRNO);
    return eSOCK_EUnknown;
  }

#ifdef OS_UNIX
  {{
    /* Let more than one "bind()" to use the same address.
     * It was affirmed(?) that under Solaris 2.5 this precaution:
     * 1) makes the address to be released immediately after the process
     *    termination
     * 2) still issue EADDINUSE error on the attempt to bind() to the
     *    same address being in-use by a living process(if SOCK_STREAM)
     */
    int reuse_addr = 1;
    if (setsockopt(sock, SOL_SOCKET, SO_REUSEADDR, 
                   (const char *)&reuse_addr, sizeof(reuse_addr)) != 0) {
      ErrPostEx(SEV_ERROR, SOCK_ERRCODE, 2,
                "[LSOCK_Create]  setsockopt():  errno = %d", (int)SOCK_ERRNO);
      SOCK_CLOSE(sock);
      return eSOCK_EUnknown;
    }
  }}
#endif

  /* Bind */
  Nlm_MemSet(&addr, '\0', sizeof(addr));
  addr.sin_family = AF_INET;
  addr.sin_addr.s_addr = htonl(INADDR_ANY);
  addr.sin_port = htons(port);
  if (bind(sock, (struct sockaddr *)&addr, sizeof(struct sockaddr)) != 0) {
    ErrPostEx(SEV_ERROR, SOCK_ERRCODE, 3,
              "[LSOCK_Create]  bind():  errno = %d", (int)SOCK_ERRNO);
    SOCK_CLOSE(sock);
    return eSOCK_EUnknown;
  }

  /* Listen */
  if (listen(sock, n_listen) != 0) {
    ErrPostEx(SEV_ERROR, SOCK_ERRCODE, 4,
              "[LSOCK_Create]  listen():  errno = %d", (int)SOCK_ERRNO);
    SOCK_CLOSE(sock);
    return eSOCK_EUnknown;
  }

  /* Set to non-blocking mode */
  if ( !s_SetNonblock(sock, TRUE) ) {
    ErrPostEx(SEV_ERROR, SOCK_ERRCODE, 4,
              "[LSOCK_Create]  Cannot set to non-blocking mode");
    SOCK_CLOSE(sock);
    return eSOCK_EUnknown;
  }

  /* Success... */
  *lsock = (LSOCK)Nlm_MemGet(sizeof(Nlm_LSOCKstruct), 0);
  (*lsock)->sock = sock;
  return eSOCK_ESuccess;
}


NLM_EXTERN ESOCK_ErrCode LSOCK_Accept(LSOCK           lsock,
                                      const STimeout *timeout,
                                      SOCK           *sock)
{
  Nlm_Socket x_sock;
  Nlm_Uint4  x_host;
  Nlm_Uint2  x_port;

  {{ /* wait for the connection(up to timeout) */
    ESOCK_ErrCode code;
    struct timeval tv;
    if ( timeout )
      s_to2tv(timeout, &tv);

    code = s_Select(lsock->sock, eSOCK_OnRead, timeout ? &tv : 0);
    if (code != eSOCK_ESuccess)
      return code;
  }}

  {{ /* accept next connection */
    struct sockaddr addr;
    struct sockaddr_in *x_addr = (struct sockaddr_in *)&addr;
#ifdef OS_MAC
    long addrlen = sizeof(addr);
#else
    int addrlen = sizeof(addr);
#endif
    if ((x_sock = accept(lsock->sock, &addr, &addrlen))== SOCK_INVALID) {
      ErrPostEx(SEV_ERROR, SOCK_ERRCODE, 1,
                "[LSOCK_Accept]  accept():  errno = %d", (int)SOCK_ERRNO);
      return eSOCK_EUnknown;
    }
    x_host = x_addr->sin_addr.s_addr;
    x_port = x_addr->sin_port;
  }}

  /* success:  create new SOCK structure */
  *sock = (SOCK)Nlm_MemNew(sizeof(Nlm_SOCKstruct));
  (*sock)->sock = x_sock;
  VERIFY ( SOCK_SetTimeout(*sock, eSOCK_OnReadWrite, 0, 0, 0)
           == eSOCK_ESuccess );
  (*sock)->host = x_host;
  (*sock)->port = x_port;
  return eSOCK_ESuccess;
}


NLM_EXTERN ESOCK_ErrCode LSOCK_Close(LSOCK lsock)
{
  int code;

  /* Set the socket back to blocking mode */
  if ( !s_SetNonblock(lsock->sock, FALSE) ) {
    ErrPostEx(SEV_WARNING, SOCK_ERRCODE, 1,
              "[LSOCK_Close]  Cannot set socket back to blocking mode");
  }

  do {  /* auto-resume if interrupted by a signal */
    code = SOCK_CLOSE(lsock->sock);
    if (code != 0  &&  SOCK_ERRNO != SOCK_EINTR) {
      ErrPostEx(SEV_ERROR, SOCK_ERRCODE, 2,
                "[LSOCK_Close] close():  errno = %d", (int)SOCK_ERRNO);
      return eSOCK_EUnknown;
    }
  } while (code != 0);

  Nlm_MemFree(lsock);
  return eSOCK_ESuccess;
}


NLM_EXTERN ESOCK_ErrCode SOCK_Create(const Nlm_Char *host,
                                     Nlm_Uint2       port,
                                     SOCK           *sock)
{
  Nlm_Socket x_sock;
  Nlm_Uint4  x_host;
  Nlm_Uint2  x_port;

  struct sockaddr_in server;
  Nlm_MemSet(&server, '\0', sizeof(server));

  /* Initialize internals */
  VERIFY ( s_Initialize() == eSOCK_ESuccess );

  /* Get address of the remote host */
  x_host = inet_addr(host);
  if (x_host == INADDR_NONE) {
    struct hostent *hp;
#ifdef OS_UNIX_SOL
    struct hostent x_hp;
    char x_buf[1024];
    int  x_err;
    hp = gethostbyname_r(host, &x_hp, x_buf, sizeof(x_buf), &x_err);
#else
    hp = gethostbyname(host);
#endif
    if ( !hp ) {
      ErrPostEx(SEV_ERROR, SOCK_ERRCODE, 1,
                "[SOCK_Create]  Cannot resolve network address(\"%s\")",
                host);
      return eSOCK_EUnknown;
    }
    Nlm_MemCpy((void *)&x_host, (void *)hp->h_addr, sizeof(x_host));
  }
  x_port = htons(port);

  /* Fill out the "server" struct */
  Nlm_MemCpy((void *)&server.sin_addr, (void *)&x_host, sizeof(x_host));
  server.sin_family = AF_INET;
  server.sin_port = x_port;

  /* Create new socket */
  if ((x_sock = socket(AF_INET, SOCK_STREAM, 0)) == SOCK_INVALID) {
    ErrPostEx(SEV_ERROR, SOCK_ERRCODE, 2,
              "[SOCK_Create]  Cannot create socket, errno = %d",
              (int)SOCK_ERRNO);
    return eSOCK_EUnknown;
  }

  /* Establish connection to the server */  
  if (connect(x_sock, (struct sockaddr *)&server, sizeof(server)) != 0) {
    ErrPostEx(SEV_ERROR, SOCK_ERRCODE, 3,
              "[SOCK_Create]  Cannot connect to host \"%s\", port %d",
              host, (int)port);
    SOCK_CLOSE(x_sock);
    return eSOCK_EUnknown;
  }

  /* Set the socket i/o to non-blocking mode */
  if ( !s_SetNonblock(x_sock, TRUE) ) {
    ErrPostEx(SEV_ERROR, SOCK_ERRCODE, 4,
              "[SOCK_Create]  Cannot set socket i/o to non-blocking mode");
    SOCK_CLOSE(x_sock);
    return eSOCK_EUnknown;
  }

  /* Success */
  *sock = (SOCK)Nlm_MemNew(sizeof(Nlm_SOCKstruct));
  (*sock)->sock = x_sock;
  VERIFY ( SOCK_SetTimeout(*sock, eSOCK_OnReadWrite, 0, 0, 0)
           == eSOCK_ESuccess );
  (*sock)->host = x_host;
  (*sock)->port = x_port;

  return eSOCK_ESuccess;
}


NLM_EXTERN ESOCK_ErrCode SOCK_Close(SOCK sock)
{
  int code;

  /* Set the socket back to blocking mode */
  if ( !s_SetNonblock(sock->sock, FALSE) ) {
    ErrPostEx(SEV_ERROR, SOCK_ERRCODE, 1,
              "[SOCK_Close]  Cannot set socket back to blocking mode");
  }

  /* Set the linger period according to the write timeout */
  if ( sock->w_timeout ) {
    struct linger lgr;
    lgr.l_onoff  = 1;
    lgr.l_linger = sock->w_timeout->tv_sec ? sock->w_timeout->tv_sec : 1;
    if (setsockopt(sock->sock, SOL_SOCKET, SO_LINGER, 
                   (char *)&lgr, sizeof(lgr)) != 0) {
      ErrPostEx(SEV_WARNING, SOCK_ERRCODE, 1,
                "[SOCK_Close] setsockopt():  errno = %d", (int)SOCK_ERRNO);
    }
  }   

  do {
    code = SOCK_CLOSE(sock->sock);
    if (code != 0  &&  SOCK_ERRNO != SOCK_EINTR) {
      ErrPostEx(SEV_WARNING, SOCK_ERRCODE, 2,
                "[SOCK_Close] close():  errno = %d", (int)SOCK_ERRNO);
      return (SOCK_ERRNO == SOCK_ECONNRESET || SOCK_ERRNO == SOCK_EPIPE) ?
        eSOCK_EClosed : eSOCK_EUnknown;
    }
    /* auto-resume if interrupted by a signal */
  } while (code != 0);

#ifdef SOCK_NCBI_PEEK
  BUF_Destroy(sock->buf);
#endif
  Nlm_MemFree(sock);
  return eSOCK_ESuccess;
}


NLM_EXTERN ESOCK_ErrCode SOCK_Select(SOCK            sock,
                                     ESOCK_Mode      mode,
                                     const STimeout *timeout)
{
  struct timeval tv;
  if ( timeout )
    s_to2tv(timeout, &tv);

  return s_Select(sock->sock, mode, timeout ? &tv : 0);
}


NLM_EXTERN ESOCK_ErrCode SOCK_SetTimeout(SOCK           sock,
                                         ESOCK_Mode     mode,
                                         const STimeout *new_timeout,
                                         STimeout       *r_timeout,
                                         STimeout       *w_timeout)
{
  if ( r_timeout )
    s_tv2to(sock->r_timeout, r_timeout);
  if ( w_timeout )
    s_tv2to(sock->w_timeout, w_timeout);

  if (new_timeout == SOCK_GET_TIMEOUT)
    return eSOCK_ESuccess;

  switch ( mode )
    {
    case eSOCK_OnRead:
      if ( new_timeout ) {
        sock->r_timeout = &sock->rr_timeout;
        s_to2tv(new_timeout, sock->r_timeout);
      } else {
        sock->r_timeout = 0;
      }
      break;

    case eSOCK_OnWrite:
      if ( new_timeout ) {
        sock->w_timeout = &sock->ww_timeout;
        s_to2tv(new_timeout, sock->w_timeout);
      } else {
        sock->w_timeout = 0;
      }
      break;

    case eSOCK_OnReadWrite:
      if ( new_timeout ) {
        sock->r_timeout = &sock->rr_timeout;
        s_to2tv(new_timeout, sock->r_timeout);
        sock->w_timeout = &sock->ww_timeout;
        s_to2tv(new_timeout, sock->w_timeout);
      } else {
        sock->r_timeout = 0;
        sock->w_timeout = 0;
      }
      break;

    default:
      ASSERT ( 0 );  return eSOCK_EUnknown;
    }

  return eSOCK_ESuccess;
}


NLM_EXTERN ESOCK_ErrCode SOCK_Read(SOCK        sock,
                                   Nlm_VoidPtr buf,
                                   Nlm_Uint4   size,
                                   Nlm_Uint4  *n_read)
{
  return s_Recv(sock, buf, size, n_read, FALSE);
}


NLM_EXTERN ESOCK_ErrCode SOCK_Peek(SOCK        sock,
                                   Nlm_VoidPtr buf,
                                   Nlm_Uint4   size,
                                   Nlm_Uint4  *n_read)
{
  return s_Recv(sock, buf, size, n_read, TRUE);
}


NLM_EXTERN ESOCK_ErrCode SOCK_Write(SOCK        sock,
                                    const void *buf,
                                    Nlm_Uint4   size,
                                    Nlm_Uint4  *n_written)
{
#ifdef SOCK_WRITE_SLICED
  return s_WriteSliced(sock, buf, size, n_written);
#else
  return s_WriteWhole (sock, buf, size, n_written);
#endif
}


NLM_EXTERN void SOCK_Address(SOCK         sock,
                             Nlm_Uint4   *host,
                             Nlm_Uint2   *port,
                             Nlm_Boolean  network_byte_order) {
  if ( host )
    *host = network_byte_order ? sock->host : (Nlm_Uint4)ntohl(sock->host);
  if ( port )
    *port = network_byte_order ? sock->port : (Nlm_Uint2)ntohs(sock->port);
}


NLM_EXTERN Nlm_Boolean GetHostName(Nlm_Char *name,
                                   Nlm_Uint4 namelen)
{
  int x_errno;
  ASSERT ( namelen > 0 );
  VERIFY ( s_Initialize() == eSOCK_ESuccess );
  name[0] = name[namelen-1] = '\0';
  x_errno = gethostname(name, (int)namelen);
  if (x_errno  ||  name[namelen-1]) {
    ErrPostEx(SEV_ERROR, SOCK_ERRCODE, 0,
              "[GetHostName]  Cannot get local host name;  errno = %d",
              (int)x_errno);
    name[0] = '\0';
    return FALSE;
  }
  return TRUE;
}


NLM_EXTERN Nlm_Boolean Uint4toInaddr(Nlm_Uint4 ui4_addr,
                                     Nlm_CharPtr buf, Nlm_Uint4 buf_len)
{
  struct in_addr addr_struct;
  char          *addr_string;

  if ( !buf )
    return FALSE;

  addr_struct.s_addr = ui4_addr;
  addr_string = inet_ntoa(addr_struct);
  if (!addr_string  ||  StrLen(addr_string) >= buf_len) {
    ErrPostEx(SEV_ERROR, SOCK_ERRCODE, 0,
              "[Uint4toInaddr]  Cannot convert %lx", (unsigned long)ui4_addr);
    buf[0] = '\0';
    return FALSE;
  }

  StrCpy(buf, addr_string);
  return TRUE;
}

extern Nlm_Uint4 Nlm_htonl(Nlm_Uint4 value)
{
  return (Nlm_Uint4)htonl((Nlm_Uint4)value);
}







/*
#ifdef OS_MAC
#define TEST_MODULE__NCBISOCK
#define DO_CLIENT
#endif
*/

#ifdef TEST_MODULE__NCBISOCK
/***********************************************************************
 *  TEST
 ***********************************************************************/

/* [WIN16] NOTE: make sure that STACKSIZE in your *.def file is set, and
 * it is big enough(I would recommend 32K; see also N_FIELD * N_FIELD)
 */
#define TEST_BUFSIZE 8192

/* exit server before sending the data expected by client(Test 1) */
/*#define TEST_SRV1_SHUTDOWN*/

#ifndef TEST_SRV1_SHUTDOWN
/* exit server immediately after its first client is served(Test 1) */
/*#define TEST_SRV1_ONCE*/
#endif


#define ASS_RET(expr,retcode) \
if ( !(expr) ) { ASSERT ( 0 );  return retcode; } else {;}


/* The simplest randezvous(short request-reply) test functions
 *      "TEST__client_1(SOCK sock)"
 *      "TEST__server_1(SOCK sock)"
 */

static const Nlm_Char s_C1[] = "C1";
static const Nlm_Char s_S1[] = "S1";

static Nlm_Int2 TEST__client_1(SOCK sock)
{ /* reserved ret.codes [110-119] */
  ESOCK_ErrCode err_code;
  Nlm_Uint4 n_io, n_io_done;
  Nlm_Char  buf[TEST_BUFSIZE];

  ErrPostEx(SEV_INFO, SOCK_ERRCODE, 110, "TC1()");

  n_io = Nlm_StrLen(s_C1) + 1;
  err_code = SOCK_Write(sock, s_C1, n_io, &n_io_done);
  ASS_RET((err_code == eSOCK_ESuccess  &&  n_io == n_io_done), 101);

  n_io = Nlm_StrLen(s_S1) + 1;
  err_code = SOCK_Read(sock, buf, n_io, &n_io_done);
  if (err_code == eSOCK_EClosed) {
    ErrPostEx(SEV_WARNING, SOCK_ERRCODE, 103, "TC1(): connection closed");
    return 103;
  }
  ASS_RET((err_code == eSOCK_ESuccess  &&  n_io == n_io_done), 104);
  ASS_RET((Nlm_StrCmp(buf, s_S1) == 0), 105);

  return 0;
}


static Nlm_Int2 TEST__server_1(SOCK sock)
{ /* reserved ret.codes [210-219] */
  ESOCK_ErrCode err_code;
  Nlm_Uint4     n_io, n_io_done;
  Nlm_Char      buf[TEST_BUFSIZE];

  ErrPostEx(SEV_INFO, SOCK_ERRCODE, 210, "TS1()");

  n_io = Nlm_StrLen(s_C1) + 1;
  err_code = SOCK_Read(sock, buf, n_io, &n_io_done);
  ASS_RET((err_code == eSOCK_ESuccess  &&  n_io == n_io_done), 210);
  ASS_RET((Nlm_StrCmp(buf, s_C1) == 0), 211);

#ifdef TEST_SRV1_SHUTDOWN
  return 212;
#endif

  n_io = Nlm_StrLen(s_S1) + 1;
  err_code = SOCK_Write(sock, s_S1, n_io, &n_io_done);
  ASS_RET((err_code == eSOCK_ESuccess  &&  n_io == n_io_done), 213);

  return 0;
}


/* More complicated randezvous test functions
 *      "TEST__client_2(SOCK sock)"
 *      "TEST__server_2(SOCK sock)"
 */

static void s_DoubleTimeout(STimeout *to) {
  if (!to->sec  &&  !to->usec) {
    to->usec = 1;
  } else {
    to->sec   = 2 * to->sec + (2 * to->usec) / 1000000;
    to->usec = (2 * to->usec) % 1000000;
  }
}

static Nlm_Int2 TEST__client_2(SOCK sock)
{ /* reserved ret.codes [120-139] */
  ESOCK_ErrCode err_code;
  Nlm_Uint4     n_io, n_io_done, i;
#define W_FIELD  10
#define N_FIELD  1000
#define N_REPEAT 10
  Nlm_Char      buf[W_FIELD * N_FIELD + 1];

  ErrPostEx(SEV_INFO, SOCK_ERRCODE, 110, "TC2()");

  /* fill out a buffer to send to server */
  Nlm_MemSet(buf, '\0', sizeof(buf));
  for (i = 0;  i < N_FIELD;  i++) {
    sprintf(buf + i * W_FIELD, "%10lu", (unsigned long)i);
  }

  /* send the buffer to server, then get it back */
  for (i = 0;  i < N_REPEAT;  i++)
    {
      Nlm_Char    buf1[sizeof(buf)];
      STimeout    w_to, r_to;
      Nlm_Boolean w_timeout_on = (Nlm_Boolean)(i%2); /* if to start from     */
      Nlm_Boolean r_timeout_on = (Nlm_Boolean)(i%3); /* zero or inf. timeout */
      Nlm_CharPtr x_buf;

      /* send */
      w_to.sec  = 0;
      w_to.usec = 0;
      err_code = SOCK_SetTimeout(sock, eSOCK_OnWrite,
                                 (w_timeout_on ? &w_to : 0), 0, 0);
      ASS_RET((err_code == eSOCK_ESuccess), 111);

      x_buf = buf;
      n_io = sizeof(buf);
      do {
#if defined(OS_UNIX)
        sleep(1);
#elif defined(WIN32)
        Sleep(1000);
#endif
        err_code = SOCK_Write(sock, x_buf, n_io, &n_io_done);
        if (err_code == eSOCK_EClosed) {
          ErrPostEx(SEV_ERROR, SOCK_ERRCODE, 112,
                    "TC2:write: connection closed");
          return 112;
        }
        ErrPostEx(SEV_INFO, SOCK_ERRCODE, 113,
                  "TC2:write: i=%d, err_code=%d, n_io=%5lu, n_io_done=%5lu"
                  "\ntimeout(%d): %5lu sec, %6lu msec",
                  (int)i, (int)err_code,
                  (unsigned long)n_io, (unsigned long)n_io_done,
                  (int)w_timeout_on,
                  (unsigned long)w_to.sec, (unsigned long)w_to.usec);
        if ( !w_timeout_on ) {
          ASS_RET((err_code == eSOCK_ESuccess  &&  n_io_done == n_io), 113);
        } else {
          STimeout x_to;
          ASS_RET((err_code == eSOCK_ESuccess  ||  err_code == eSOCK_ETimeout),
                  114);
          ASS_RET((SOCK_SetTimeout(sock, eSOCK_OnWrite, SOCK_GET_TIMEOUT,
                                   0, &x_to) == eSOCK_ESuccess  &&
                   w_to.sec == x_to.sec  &&  w_to.usec == x_to.usec), 115);
        }
        n_io  -= n_io_done;
        x_buf += n_io_done;
        if (err_code == eSOCK_ETimeout)
          s_DoubleTimeout(&w_to);
        err_code = SOCK_SetTimeout(sock, eSOCK_OnWrite, &w_to, 0, 0);
        ASS_RET((err_code == eSOCK_ESuccess), 116);
        w_timeout_on = TRUE;
      } while ( n_io );

      /* get back the just sent data */
      r_to.sec  = 0;
      r_to.usec = 0;
      err_code = SOCK_SetTimeout(sock, eSOCK_OnRead,
                                 (r_timeout_on ? &r_to : 0), 0, 0);
      ASS_RET((err_code == eSOCK_ESuccess), 121);

      x_buf = buf1;
      n_io = sizeof(buf1);
      do {
        if (i%2 == 0)
          { /* peek a little piece twice and compare */
            Nlm_Char  xx_buf1[128], xx_buf2[128];
            Nlm_Uint4 xx_io_done1, xx_io_done2;
            if (SOCK_Peek(sock, xx_buf1, sizeof(xx_buf1), &xx_io_done1)
                == eSOCK_ESuccess  &&
                SOCK_Peek(sock, xx_buf2, xx_io_done1, &xx_io_done2)
                == eSOCK_ESuccess) {
              ASSERT ( xx_io_done1 >= xx_io_done2 );
              VERIFY ( !Nlm_MemCmp(xx_buf1, xx_buf2, xx_io_done2) );
            }
          }
        err_code = SOCK_Read(sock, x_buf, n_io, &n_io_done);
        if (err_code == eSOCK_EClosed) {
          ErrPostEx(SEV_ERROR, SOCK_ERRCODE, 122,
                    "TC2:read: connection closed");
          return 122;
        }
        ErrPostEx(SEV_INFO, SOCK_ERRCODE, 123,
                  "TC2:read: i=%d, err_code=%d, n_io=%5lu, n_io_done=%5lu"
                  "\ntimeout(%d): %5lu sec, %6lu msec",
                  (int)i, (int)err_code,
                  (unsigned long)n_io, (unsigned long)n_io_done,
                  (int)r_timeout_on,
                  (unsigned long)r_to.sec, (unsigned long)r_to.usec);
        if ( !r_timeout_on ) {
          ASS_RET((err_code == eSOCK_ESuccess  &&  n_io_done > 0), 124);
        } else {
          STimeout x_to;
          ASS_RET((err_code == eSOCK_ESuccess  ||  err_code == eSOCK_ETimeout),
                  125);
          ASS_RET((SOCK_SetTimeout(sock, eSOCK_OnRead, SOCK_GET_TIMEOUT,
                                   &x_to, 0) == eSOCK_ESuccess  &&
                   r_to.sec == x_to.sec  &&  r_to.usec == x_to.usec), 126);
        }

        n_io  -= n_io_done;
        x_buf += n_io_done;
        if (err_code == eSOCK_ETimeout)
          s_DoubleTimeout(&r_to);
        err_code = SOCK_SetTimeout(sock, eSOCK_OnRead, &r_to, 0, 0);
        ASS_RET((err_code == eSOCK_ESuccess), 127);
        r_timeout_on = TRUE;
      } while ( n_io );

      ASS_RET((!Nlm_MemCmp(buf, buf1, sizeof(buf))), 120); 
    }

  return 0;
}


static Nlm_Int2 TEST__server_2(SOCK sock)
{ /* reserved ret.codes [220-229] */
  ESOCK_ErrCode err_code;
  Nlm_Uint4     n_io, n_io_done;
  Nlm_Char      buf[TEST_BUFSIZE];
  STimeout      r_to, w_to;
  Nlm_Uint4     i;

  ErrPostEx(SEV_INFO, SOCK_ERRCODE, 220, "TS2()");

  r_to.sec  = 0;
  r_to.usec = 0;
  w_to = r_to;
  err_code = SOCK_SetTimeout(sock, eSOCK_OnReadWrite, &r_to, 0, 0);
  ASS_RET((err_code == eSOCK_ESuccess), 221);

  for (i = 0;  ;  i++) {
    Nlm_CharPtr x_buf;

    /* read data from socket */
    n_io = sizeof(buf);
    err_code = SOCK_Read(sock, buf, n_io, &n_io_done);
    switch ( err_code )
      {
      case eSOCK_ESuccess:
        ErrPostEx(SEV_INFO, SOCK_ERRCODE, 222,
                  "TS2:read:[%lu], err_code=%d, n_io=%5lu, n_io_done=%5lu",
                  (unsigned long)i, (int)err_code,
                  (unsigned long)n_io, (unsigned long)n_io_done);
        ASS_RET((n_io_done > 0), 222);
        break;
      case eSOCK_EClosed:
        ErrPostEx(SEV_INFO, SOCK_ERRCODE, 223,
                  "TS2:read: connection closed");
        return 0;
      case eSOCK_ETimeout:
        ErrPostEx(SEV_INFO, SOCK_ERRCODE, 224,
                  "TS2:read:[%lu] timeout expired: %5lu sec, %6lu msec",
                  (unsigned long)i,
                  (unsigned long)r_to.sec, (unsigned long)r_to.usec);
        ASS_RET((n_io_done == 0), 224);
        s_DoubleTimeout(&r_to);
        err_code = SOCK_SetTimeout(sock, eSOCK_OnRead, &r_to, 0, 0);
        ASS_RET((err_code == eSOCK_ESuccess), 225);
        break;
      default:
        ASS_RET(0, 226);
      }

    /* write(just the same) data back to client */
    n_io  = n_io_done;
    x_buf = buf;
    while ( n_io ) {
      err_code = SOCK_Write(sock, buf, n_io, &n_io_done);
      switch ( err_code )
        {
        case eSOCK_ESuccess:
          ErrPostEx(SEV_INFO, SOCK_ERRCODE, 231,
                    "TS2:write:[%lu], err_code=%d, n_io=%5lu, n_io_done=%5lu",
                    (unsigned long)i, (int)err_code,
                    (unsigned long)n_io, (unsigned long)n_io_done);
          ASS_RET((n_io_done > 0), 231);
          break;
        case eSOCK_EClosed:
          ErrPostEx(SEV_ERROR, SOCK_ERRCODE, 232,
                    "TS2:write: connection closed");
          return 230;
        case eSOCK_ETimeout:
          ErrPostEx(SEV_INFO, SOCK_ERRCODE, 233,
                    "TS2:write:[%lu] timeout expired: %5lu sec, %6lu msec",
                    (unsigned long)i,
                    (unsigned long)w_to.sec, (unsigned long)w_to.usec);
          ASS_RET((n_io_done == 0), 233);
          s_DoubleTimeout(&w_to);
          err_code = SOCK_SetTimeout(sock, eSOCK_OnWrite, &w_to, 0, 0);
          ASS_RET((err_code == eSOCK_ESuccess), 234);
          break;
        default:
          ASS_RET(0, 235);
        }

      n_io  -= n_io_done;
      x_buf += n_io_done;
    }
  }

  return 0;
}


/* Skeletons for the socket i/o test:
 *   TEST__client(...)
 *   TEST__server(...)
 *   establish and close connection;  call test i/o functions like
 *     TEST__[client|server]_[1|2|...] (...)
 */
static Nlm_Int2 TEST__client(const Nlm_Char *server_host,
                             Nlm_Uint2       server_port)
{ /* reserved ret.codes [100-109] */
  SOCK          sock;
  ESOCK_ErrCode err_code;
  Nlm_Int2      ret_code;

  ErrPostEx(SEV_INFO, SOCK_ERRCODE, 101,
            "TEST__client(host = \"%s\", port = %u)",
            server_host, (unsigned)server_port);

  /* Connect to server */
  err_code = SOCK_Create(server_host, server_port, &sock);
  ASS_RET((err_code == eSOCK_ESuccess), 100);

  /* Test the simplest randezvous(short request-reply)
   * The two peer functions are:
   *      "TEST__[client|server]_1(SOCK sock)"
   */
  ret_code = TEST__client_1(sock);
  ASS_RET((ret_code == 0), 101);

  /* Test a more complex case
   * The two peer functions are:
   *      "TEST__[client|server]_2(SOCK sock)"
   */
  ret_code = TEST__client_2(sock);
  ASS_RET((ret_code == 0), 102);

  /* Close connection and exit */
  err_code = SOCK_Close(sock);
  ASS_RET((err_code == eSOCK_ESuccess  ||  err_code == eSOCK_EClosed), 109);
  return 0;
}


static Nlm_Int2 TEST__server(Nlm_Uint2 port)
{ /* reserved ret.codes [200-209] */
  LSOCK lsock;
  ESOCK_ErrCode err_code;

  ErrPostEx(SEV_INFO, SOCK_ERRCODE, 201,
            "TEST__server(port = %u)", (unsigned)port);

  /* Create listening socket */
  err_code = LSOCK_Create(port, 1, &lsock);
  ASS_RET((err_code == eSOCK_ESuccess), 200);

  /* Accept connections from clients and run test sessions */
  for (;;)
    {
      Nlm_Int2 ret_code;

      /* Accept connection */
      SOCK sock;
      err_code = LSOCK_Accept(lsock, NULL, &sock);
      ASS_RET((err_code == eSOCK_ESuccess), 208);

      /* Test the simplest randezvous(short request-reply)
       * The two peer functions are:
       *      "TEST__[client|server]_1(SOCK sock)"
       */
      ret_code = TEST__server_1(sock);
      ASS_RET((ret_code == 0), 201);

      /* Test a more complex case
       * The two peer functions are:
       *      "TEST__[client|server]_2(SOCK sock)"
       */
      ret_code = TEST__server_2(sock);
      ASS_RET((ret_code == 0), 202);

      /* Close connection */
      err_code = SOCK_Close(sock);
      ASS_RET((err_code == eSOCK_ESuccess || err_code == eSOCK_EClosed), 209);

#ifdef TEST_SRV1_ONCE
      /* finish after the first session */
      break;
#endif
    }

  /* Close listening socket */
  err_code = LSOCK_Close(lsock);
  ASS_RET((err_code == eSOCK_ESuccess), 204);
  return 0;
}


/* Main function
 * Parse command-line options, initialize and cleanup API internals;
 * run client or server test
 */
extern Nlm_Int2 Nlm_Main(void)
{
#define MIN_PORT 5001
  Nlm_Int4     argc = Nlm_GetArgc();
  Nlm_CharPtr *argv = Nlm_GetArgv();

#ifdef WIN16
  {{ /* a kludge to make sure the "vibwndws.c"(MainProc) get linked */
    extern void Nlm_Metronome(Nlm_VoidPtr actn);  Nlm_Metronome(0);
  }}
#endif

  Nlm_ErrSetOpts(ERR_TEE, ERR_LOG_ON);
  ErrSetLogLevel(SEV_INFO);
  ErrSetMessageLevel(SEV_INFO);
  VERIFY ( Nlm_ErrSetLog("ncbisock.log") );

  {{
    Nlm_Char local_host[64];
    VERIFY ( GetHostName(local_host, sizeof(local_host)) );
    ErrPostEx(SEV_INFO, SOCK_ERRCODE, 200,
              "\nRunning NCBISOCK test on host \"%s\"", local_host);
  }}

#ifdef DO_SERVER
  argc = 2;
#endif
#ifdef DO_CLIENT
  argc = 3;
#endif

  switch ( argc )
    {
    case 2:
      { /* Server */
#ifdef DO_SERVER
        short port = 5555;
#else
        short port;
        if (sscanf(argv[1], "%hd", &port) != 1  ||
            port < MIN_PORT)
          break;
#endif
        ErrPostEx(SEV_INFO, SOCK_ERRCODE, 200,
                  "Starting NCBISOCK server test...");
        VERIFY ( Nlm_ErrSetLog("ncbisock.srv") );

        {{
          Nlm_Int2 ret_code = TEST__server((Nlm_Uint2)port);
          VERIFY ( SOCK_Destroy() == eSOCK_ESuccess );
          return ret_code;
        }}
      }

    case 3:
      { /* Client */
#ifdef DO_CLIENT
        Nlm_CharPtr server_host = "peony";
        short       server_port = 5555;
#else
        Nlm_CharPtr server_host = argv[1];
        short       server_port;
        if (sscanf(argv[2], "%hd", &server_port) != 1  ||
            server_port < MIN_PORT)
          break;
#endif
        ErrPostEx(SEV_INFO, SOCK_ERRCODE, 100,
                  "Starting NCBISOCK client test... ");
        VERIFY ( Nlm_ErrSetLog("ncbisock.cli") );

        {{
          Nlm_Int2 ret_code= TEST__client(server_host, (Nlm_Uint2)server_port);
          VERIFY ( SOCK_Destroy() == eSOCK_ESuccess );
          return ret_code;
        }}
      }
    }

  /* Bad cmd-line arguments;  Usage */
  ErrPostEx(SEV_ERROR, SOCK_ERRCODE, 666,
            "Usage:\n"
            "  Client: %s <srv_host> <port>\n"
            "  Server: %s <port>\n"
            " where <port> not less than %hd",
            argv[0], argv[0], (short)MIN_PORT);
  return 1;
}

#endif /* TEST_MODULE__NCBISOCK */

/* EOF */

