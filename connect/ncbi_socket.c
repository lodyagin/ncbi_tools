/*  $Id: ncbi_socket.c,v 6.6 2000/02/25 20:31:52 vakatov Exp $
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
 *     [UNIX ]   -DNCBI_OS_UNIX     -lresolv -lsocket -lnsl
 *     [MSWIN]   -DNCBI_OS_MSWIN    wsock32.lib
 *     [MacOS]   -DNCBI_OS_MAC      NCSASOCK -- BSD-style socket emulation lib
 *
 * ---------------------------------------------------------------------------
 * $Log: ncbi_socket.c,v $
 * Revision 6.6  2000/02/25 20:31:52  vakatov
 * Rollback to R6.4
 *
 * Revision 6.4  2000/02/23 22:34:35  vakatov
 * Can work both "standalone" and as a part of NCBI C++ or C toolkits
 *
 * Revision 6.3  1999/10/19 22:21:48  vakatov
 * Put the call to "gethostbyname()" into a CRITICAL_SECTION
 *
 * Revision 6.2  1999/10/19 16:16:01  vakatov
 * Try the NCBI C and C++ headers only if NCBI_OS_{UNIX, MSWIN, MAC} is
 * not #define'd
 *
 * Revision 6.1  1999/10/18 15:36:38  vakatov
 * Initial revision (derived from the former "ncbisock.[ch]")
 *
 * ===========================================================================
 */


/* OS must be specified in the command-line ("-D....") or in the conf. header
 */
#if !defined(NCBI_OS_UNIX) && !defined(NCBI_OS_MSWIN) && !defined(NCBI_OS_MAC)
#  error "Unknown OS, must be one of NCBI_OS_UNIX, NCBI_OS_MSWIN, NCBI_OS_MAC!"
#endif


/* Uncomment this (or specify "-DHAVE_GETHOSTBYNAME_R=1") only if:
 * 0) you are compiling this outside of the NCBI C or C++ Toolkits
 *    (USE_NCBICONF is not #define'd), and
 * 1) your platform has "gethostbyname_r()", and
 * 2) you are going to use this API code in multi-thread application, and
 * 3) "gethostbyname()" gets called somewhere else in your code
 */
/* #define HAVE_GETHOSTBYNAME_R 1 */


/* Platform-specific system headers
 */
#if defined(NCBI_OS_UNIX)
#  include <errno.h>
#  include <sys/socket.h>
#  include <sys/time.h>
#  include <unistd.h>
#  include <netdb.h>
#  include <netinet/in.h>
#  include <fcntl.h>
#  include <arpa/inet.h>
#  include <signal.h>

#elif defined(NCBI_OS_MSWIN)
#  include <winsock.h>

#elif defined(NCBI_OS_MAC)
#  include <macsockd.h>
#  define __TYPES__       /* avoid Mac <Types.h> */
#  define __MEMORY__      /* avoid Mac <Memory.h> */
#  define ipBadLapErr     /* avoid Mac <MacTCPCommonTypes.h> */
#  define APPL_SOCK_DEF
#  define SOCK_DEFS_ONLY
#  include <sock_ext.h>
extern void bzero(char* target, long numbytes);
#  include <netdb.h>
#  include <s_types.h>
#  include <s_socket.h>
#  include <s_ioctl.h>
#  include <neti_in.h>
#  include <a_inet.h>
#  include <s_time.h>
#  include <s_fcntl.h>
#  include <neterrno.h> /* missing error numbers on Mac */
/* cannot write more than SOCK_WRITE_SLICE at once on Mac; so we have to split
 * big output buffers into smaller(<SOCK_WRITE_SLICE) slices before writing
 * them to the socket
 */
#  define SOCK_WRITE_SLICE 2048

#else
#  error "Unsupported platform, must be one of NCBI_OS_UNIX, NCBI_OS_MSWIN, NCBI_OS_MAC !!!"
#endif /* platform-specific headers (for UNIX, MSWIN, MAC) */


/* NCBI core headers
 */
#include <connect/ncbi_util.h>
#include <connect/ncbi_buffer.h>
#include <connect/ncbi_socket.h>


/* Portable standard C headers
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>




/******************************************************************************
 *  TYPEDEF & MACRO
 */


/* Minimal size of the data buffer chunk in the socket internal buffer(s) */
#define SOCK_BUF_CHUNK_SIZE 4096


/* Macro #def for the platform-dependent constants, error codes and functions
 */
#if defined(NCBI_OS_MSWIN)

typedef SOCKET TSOCK_Handle;
#  define SOCK_INVALID     INVALID_SOCKET
#  define SOCK_ERRNO       WSAGetLastError()
#  define SOCK_EINTR       WSAEINTR
#  define SOCK_EWOULDBLOCK WSAEWOULDBLOCK
#  define SOCK_ECONNRESET  WSAECONNRESET
#  define SOCK_EPIPE       WSAESHUTDOWN
#  define SOCK_EAGAIN      WSAEINPROGRESS
#  define SOCK_EINPROGRESS WSAEINPROGRESS
#  define SOCK_NFDS(s)     0
#  define SOCK_CLOSE(s)    closesocket(s)
/* NCBI_OS_MSWIN */

#elif defined(NCBI_OS_UNIX)

typedef int TSOCK_Handle;
#  define SOCK_INVALID     (-1)
#  define SOCK_ERRNO       errno
#  define SOCK_EINTR       EINTR
#  define SOCK_EWOULDBLOCK EWOULDBLOCK
#  define SOCK_ECONNRESET  ECONNRESET
#  define SOCK_EPIPE       EPIPE
#  define SOCK_EAGAIN      EAGAIN
#  define SOCK_EINPROGRESS EINPROGRESS
#  define SOCK_NFDS(s)     (s + 1)
#  define SOCK_CLOSE(s)    close(s)
#  ifndef INADDR_NONE
#    define INADDR_NONE    (unsigned int)(-1)
#  endif
/* NCBI_OS_UNIX */

#elif defined(NCBI_OS_MAC)

typedef int TSOCK_Handle;
#  define SOCK_INVALID     (-1)
#  ifndef SOCK_ERRNO
#    define SOCK_ERRNO     errno
#  endif
#  define SOCK_EINTR       EINTR
#  define SOCK_EWOULDBLOCK EWOULDBLOCK
#  define SOCK_ECONNRESET  ECONNRESET
#  define SOCK_EPIPE       EPIPE
#  define SOCK_EAGAIN      EAGAIN
#  define SOCK_EINPROGRESS EINPROGRESS
#  define SOCK_NFDS(s)     (s + 1)
#  define SOCK_CLOSE(s)    close(s)

/* (but see ni_lib.c line 2508 for gethostname substitute for Mac) */
extern int gethostname(char* machname, long buflen);

#endif /* NCBI_OS_MSWIN, NCBI_OS_UNIX, NCBI_OS_MAC */


/* Listening socket
 */
typedef struct LSOCK_tag {
    TSOCK_Handle sock;  /* OS-specific socket handle */
} LSOCK_struct;


/* Socket
 */
typedef struct SOCK_tag {
    TSOCK_Handle    sock;       /* OS-specific socket handle */
    struct timeval* r_timeout;  /* zero (if infinite), or points to "r_tv" */
    struct timeval  r_tv;       /* finite read timeout value */
    STimeout        r_to;       /* finite read timeout value (aux., temp.) */
    struct timeval* w_timeout;  /* zero (if infinite), or points to "w_tv" */
    struct timeval  w_tv;       /* finite write timeout value */
    STimeout        w_to;       /* finite write timeout value (aux., temp.) */
    struct timeval* c_timeout;  /* zero (if infinite), or points to "c_tv" */
    struct timeval  c_tv;       /* finite close timeout value */
    STimeout        c_to;       /* finite close timeout value (aux., temp.) */
    unsigned int    host;       /* peer host (in the network byte order) */
    unsigned short  port;       /* peer port (in the network byte order) */
    BUF             buf;        /* read buffer */
    int/*bool*/     is_eof;     /* if EOF has been detected (on read) */
} SOCK_struct;



/******************************************************************************
 *  Multi-Thread SAFETY
 *
 *  To protect static data, such like "s_Log" and "s_Initialized".
 */

#define SOCK_LOCK_WRITE  verify( MT_LOCK_Do(s_MT_Lock, eMT_Lock    ) )
#define SOCK_LOCK_READ   verify( MT_LOCK_Do(s_MT_Lock, eMT_LockRead) )
#define SOCK_UNLOCK      verify( MT_LOCK_Do(s_MT_Lock, eMT_Unlock  ) )

static MT_LOCK s_MT_Lock;

extern void SOCK_SetLOCK(MT_LOCK lk)
{
    if (s_MT_Lock  &&  lk != s_MT_Lock) {
        MT_LOCK_Delete(s_MT_Lock);
    }
    s_MT_Lock = lk;
}



/******************************************************************************
 *  ERROR HANDLING
 */

#define SOCK_LOG(level, message)  do { \
    if ( s_Log ) { \
        SOCK_LOCK_READ; \
        LOG_WRITE(s_Log, level, message); \
        SOCK_UNLOCK; \
    } \
} while(0)

#define SOCK_LOG_ERRNO(x_errno, level, message)  do { \
    if ( s_Log ) { \
        char buf[2048]; \
        int xx_errno = x_errno; \
        SOCK_LOCK_READ; \
        LOG_WRITE(s_Log, level, \
                  MessagePlusErrno(message, xx_errno, "", buf, sizeof(buf))); \
        SOCK_UNLOCK; \
    } \
} while(0)


static LOG s_Log;

extern void SOCK_SetLOG(LOG lg)
{
    SOCK_LOCK_WRITE;
    if (s_Log  &&  lg != s_Log) {
        LOG_Delete(s_Log);
    }
    s_Log = lg;
    SOCK_UNLOCK;
}



/******************************************************************************
 *  API Initialization and Shutdown/Cleanup
 */


static int/*bool*/ s_Initialized = 0/*false*/;


extern EIO_Status SOCK_InitializeAPI(void)
{
    SOCK_LOCK_WRITE;

    if ( s_Initialized ) {
        SOCK_UNLOCK;
        return eIO_Success;
    }

#if defined(NCBI_OS_MSWIN)
    {{
        WSADATA wsadata;
        int x_errno = WSAStartup(MAKEWORD(1,1), &wsadata);
        if (x_errno != 0) {
            SOCK_UNLOCK;
            SOCK_LOG_ERRNO(x_errno, eLOG_Error,
                           "[SOCK_InitializeAPI]  Failed WSAStartup()");
            return eIO_Unknown;
        }
    }}
#elif defined(NCBI_OS_UNIX)
    signal(SIGPIPE, SIG_IGN);
#endif

    s_Initialized = 1/*true*/;
    SOCK_UNLOCK;
    return eIO_Success;
}


extern EIO_Status SOCK_ShutdownAPI(void)
{
    EIO_Status status = eIO_Success;

    SOCK_LOCK_WRITE;
    if ( s_Initialized ) {
        s_Initialized = 0/*false*/;
#if defined(NCBI_OS_MSWIN)
        {{
            int x_errno = WSACleanup() ? SOCK_ERRNO : 0;
            SOCK_UNLOCK;
            if ( x_errno ) {
                SOCK_LOG_ERRNO(x_errno, eLOG_Warning,
                               "[SOCK_ShutdownAPI]  Failed WSACleanup()");
                status = eIO_Unknown;
            }
        }}
#else
        SOCK_UNLOCK;
#endif
    } else {
        SOCK_UNLOCK;
    }

    SOCK_SetLOG(0);
    SOCK_SetLOCK(0);
    return status;
}



/******************************************************************************
 *  LSOCK & SOCK AUXILIARES
 */


/* STimeout <--> struct timeval  conversions
 */
static STimeout *s_tv2to(const struct timeval* tv, STimeout* to)
{
    if ( !tv )
        return 0;

    to->sec  = (unsigned int) tv->tv_sec;
    to->usec = (unsigned int) tv->tv_usec;
    return to;
}

static struct timeval* s_to2tv(const STimeout* to, struct timeval* tv)
{
    if ( !to )
        return 0;

    tv->tv_sec  = to->sec;
    tv->tv_usec = to->usec % 1000000;
    return tv;
}


/* Switch the specified socket i/o between blocking and non-blocking mode
 */
static int/*bool*/ s_SetNonblock(TSOCK_Handle sock, int/*bool*/ nonblock)
{
#if defined(NCBI_OS_MSWIN)
    unsigned long argp = nonblock ? 1 : 0;
    return (ioctlsocket(sock, FIONBIO, &argp) == 0);
#elif defined(NCBI_OS_UNIX)
    return (fcntl(sock, F_SETFL,
                  nonblock ?
                  fcntl(sock, F_GETFL, 0) | O_NONBLOCK :
                  fcntl(sock, F_GETFL, 0) & (int)~O_NONBLOCK) != -1);
#elif defined(NCBI_OS_MAC)
    return (fcntl(sock, F_SETFL,
                  nonblock ?
                  fcntl(sock, F_GETFL, 0) | O_NDELAY :
                  fcntl(sock, F_GETFL, 0) & (int)~O_NDELAY) != -1);
#else
    assert(0);
    return 0/*false*/;
#endif
}


/* Select on the socket i/o
 */
static EIO_Status s_Select(TSOCK_Handle          sock,
                           EIO_Event             event,
                           const struct timeval* timeout)
{
    int n_dfs;
    fd_set fds, *r_fds, *w_fds;

    /* just checking */
    if (sock == SOCK_INVALID) {
        SOCK_LOG(eLOG_Error,
                 "[SOCK::s_Select]  Attempted to wait on an invalid socket");
        assert(0);
        return eIO_Unknown;
    }

    /* setup i/o descriptor to select */
    r_fds = (event == eIO_Read  ||  event == eIO_ReadWrite) ? &fds : 0;
    w_fds = (event == eIO_Write ||  event == eIO_ReadWrite) ? &fds : 0;

    do { /* auto-resume if interrupted by a signal */
        fd_set e_fds;
        struct timeval tmout;
        if ( timeout )
            tmout = *timeout;
        FD_ZERO(&fds);       FD_ZERO(&e_fds);
        FD_SET(sock, &fds);  FD_SET(sock, &e_fds);
        n_dfs = select(SOCK_NFDS(sock), r_fds, w_fds, &e_fds,
                       timeout ? &tmout : 0);
        assert(-1 <= n_dfs  &&  n_dfs <= 2);
        if ((n_dfs < 0  &&  SOCK_ERRNO != SOCK_EINTR)  ||
            FD_ISSET(sock, &e_fds))
            return eIO_Unknown;
    } while (n_dfs < 0);

    /* timeout has expired */
    if (n_dfs == 0)
        return eIO_Timeout;

    /* success;  can i/o now */
    assert(FD_ISSET(sock, &fds));
    return eIO_Success;
}



/******************************************************************************
 *  LISTENING SOCKET
 */


extern EIO_Status LSOCK_Create(unsigned short port,
                               unsigned short backlog,
                               LSOCK*         lsock)
{
    TSOCK_Handle       x_lsock;
    struct sockaddr_in addr;

    /* Initialize internals */
    verify(s_Initialized  ||  SOCK_InitializeAPI() == eIO_Success);

    /* Create new(listening) socket */
    if ((x_lsock = socket(AF_INET, SOCK_STREAM, 0)) == SOCK_INVALID) {
        SOCK_LOG_ERRNO(SOCK_ERRNO, eLOG_Error,
                       "[LSOCK_Create]  Cannot create listening socket");
        return eIO_Unknown;
    }

#ifdef NCBI_OS_UNIX
    {{
        /* Let more than one "bind()" to use the same address.
         * It was affirmed(?) that at least under Solaris 2.5 this precaution:
         * 1) makes the address to be released immediately after the process
         *    termination
         * 2) still issue EADDINUSE error on the attempt to bind() to the
         *    same address being in-use by a living process(if SOCK_STREAM)
         */
        int reuse_addr = 1;
        if (setsockopt(x_lsock, SOL_SOCKET, SO_REUSEADDR, 
                       (const char*) &reuse_addr, sizeof(reuse_addr)) != 0) {
            SOCK_LOG_ERRNO(SOCK_ERRNO, eLOG_Error,
                           "[LSOCK_Create]  Failed setsockopt(REUSEADDR)");
            SOCK_CLOSE(x_lsock);
            return eIO_Unknown;
        }
    }}
#endif

    /* Bind */
    memset(&addr, 0, sizeof(addr));
    addr.sin_family      = AF_INET;
    addr.sin_addr.s_addr = htonl(INADDR_ANY);
    addr.sin_port        = htons(port);
    if (bind(x_lsock, (struct sockaddr*)&addr, sizeof(struct sockaddr)) != 0) {
        SOCK_LOG_ERRNO(SOCK_ERRNO, eLOG_Error,
                       "[LSOCK_Create]  Failed bind()");
        SOCK_CLOSE(x_lsock);
        return eIO_Unknown;
    }

    /* Listen */
    if (listen(x_lsock, backlog) != 0) {
        SOCK_LOG_ERRNO(SOCK_ERRNO, eLOG_Error,
                       "[LSOCK_Create]  Failed listen()");
        SOCK_CLOSE(x_lsock);
        return eIO_Unknown;
    }

    /* Set to non-blocking mode */
    if ( !s_SetNonblock(x_lsock, 1/*true*/) ) {
        SOCK_LOG(eLOG_Error,
                 "[LSOCK_Create]  Cannot set socket to non-blocking mode");
        SOCK_CLOSE(x_lsock);
        return eIO_Unknown;
    }

    /* Success... */
    *lsock = (LSOCK) calloc(1, sizeof(LSOCK_struct));
    (*lsock)->sock = x_lsock;
    return eIO_Success;
}


extern EIO_Status LSOCK_Accept(LSOCK           lsock,
                               const STimeout* timeout,
                               SOCK*           sock)
{
    TSOCK_Handle   x_sock;
    unsigned int   x_host;
    unsigned short x_port;

    {{ /* wait for the connection request to come (up to timeout) */
        struct timeval tv;
        EIO_Status status =
            s_Select(lsock->sock, eIO_Read, s_to2tv(timeout, &tv));
        if (status != eIO_Success)
            return status;
    }}

    {{ /* accept next connection */
        struct sockaddr_in addr;
#ifdef NCBI_OS_MAC
        long addrlen = sizeof(struct sockaddr);
#else
        int addrlen = sizeof(struct sockaddr);
#endif
        if ((x_sock = accept(lsock->sock, (struct sockaddr *)&addr, &addrlen))
            == SOCK_INVALID) {
            SOCK_LOG_ERRNO(SOCK_ERRNO, eLOG_Error,
                           "[LSOCK_Accept]  Failed accept()");
            return eIO_Unknown;
        }
        x_host = addr.sin_addr.s_addr;
        x_port = addr.sin_port;
    }}

    /* success:  create new SOCK structure */
    *sock = (SOCK) calloc(1, sizeof(SOCK_struct));
    (*sock)->sock = x_sock;
    verify(SOCK_SetTimeout(*sock, eIO_ReadWrite, 0) == eIO_Success);
    (*sock)->host = x_host;
    (*sock)->port = x_port;
    BUF_SetChunkSize(&(*sock)->buf, SOCK_BUF_CHUNK_SIZE);
    return eIO_Success;
}


extern EIO_Status LSOCK_Close(LSOCK lsock)
{
    /* Set the socket back to blocking mode */
    if ( !s_SetNonblock(lsock->sock, 0/*false*/) ) {
        SOCK_LOG(eLOG_Warning,
                 "[LSOCK_Close]  Cannot set socket back to blocking mode");
    }

    /* Close (persistently -- retry if interrupted by a signal) */
    for (;;) {
        /* success */
        if (SOCK_CLOSE(lsock->sock) == 0)
            break;

        /* error */
        if (SOCK_ERRNO != SOCK_EINTR) {
            SOCK_LOG_ERRNO(SOCK_ERRNO, eLOG_Error,
                           "[LSOCK_Close]  Failed close()");
            return eIO_Unknown;
        }
    }

    /* success:  cleanup & return */
    free(lsock);
    return eIO_Success;
}


extern EIO_Status LSOCK_GetOSHandle(LSOCK  lsock,
                                    void*  handle_buf,
                                    size_t handle_size)
{
    if (handle_size != sizeof(lsock->sock)) {
        SOCK_LOG(eLOG_Error, "[LSOCK_GetOSHandle]  Invalid handle size");
        assert(0);
        return eIO_Unknown;
    }

    memcpy(handle_buf, &lsock->sock, handle_size);
    return eIO_Success;
}



/******************************************************************************
 *  SOCKET
 */


/* Connect the (pre-allocated) socket to the specified "host:port" peer.
 * HINT: if "host" is NULL then assume(!) that the "sock" already exists,
 *       and connect to the same host;  the same is for zero "port".
 */
static EIO_Status s_Connect(SOCK            sock,
                            const char*     host,
                            unsigned short  port,
                            const STimeout* timeout)
{
    TSOCK_Handle   x_sock;
    unsigned int   x_host;
    unsigned short x_port;

    struct sockaddr_in server;
    memset(&server, 0, sizeof(server));

    /* Initialize internals */
    verify(s_Initialized  ||  SOCK_InitializeAPI() == eIO_Success);

    /* Get address of the remote host (assume the same host if it is NULL) */
    assert(host  ||  sock->host);
    x_host = host ? inet_addr(host) : sock->host;
    if (x_host == INADDR_NONE) {
        struct hostent *hp;
#if defined(HAVE_GETHOSTBYNAME_R)
        struct hostent x_hp;
        char x_buf[1024];
        int  x_err;
        hp = gethostbyname_r(host, &x_hp, x_buf, sizeof(x_buf), &x_err);
        if ( hp )
            memcpy(&x_host, hp->h_addr, sizeof(x_host));
#else
        SOCK_LOCK_WRITE;
        hp = gethostbyname(host);
        if ( hp )
            memcpy(&x_host, hp->h_addr, sizeof(x_host));
        SOCK_UNLOCK;
#endif
        if ( !hp ) {
            if ( s_Log ) {
                char str[128];
                sprintf(str, "[SOCK::s_Connect]  Failed gethostbyname(%.64s)",
                        host ? host : "");
                SOCK_LOG(eLOG_Error, str);
            }
            return eIO_Unknown;
        }
    }

    /* Set the port to connect to(assume the same port if "port" is zero) */
    assert(port  ||  sock->port);
    x_port = (unsigned short) (port ? htons(port) : sock->port);

    /* Fill out the "server" struct */
    memcpy(&server.sin_addr, &x_host, sizeof(x_host));
    server.sin_family = AF_INET;
    server.sin_port   = x_port;

    /* Create new socket */
    if ((x_sock = socket(AF_INET, SOCK_STREAM, 0)) == SOCK_INVALID) {
        SOCK_LOG_ERRNO(SOCK_ERRNO, eLOG_Error,
                       "[SOCK::s_Connect]  Cannot create socket");
        return eIO_Unknown;
    }

    /* Set the socket i/o to non-blocking mode */
    if ( !s_SetNonblock(x_sock, 1/*true*/) ) {
        SOCK_LOG(eLOG_Error,
                 "[SOCK::s_Connect]  Cannot set socket to non-blocking mode");
        SOCK_CLOSE(x_sock);
        return eIO_Unknown;
    }

    /* Establish connection to the peer */
    if (connect(x_sock, (struct sockaddr*) &server, sizeof(server)) != 0) {
        if (SOCK_ERRNO != SOCK_EINTR  &&  SOCK_ERRNO != SOCK_EINPROGRESS  &&
            SOCK_ERRNO != SOCK_EWOULDBLOCK) {
            if ( s_Log ) {
                char str[256];
                sprintf(str,
                        "[SOCK::s_Connect]  Failed connect() to %.64s:%d",
                        host ? host : "???", (int)ntohs(x_port));
                SOCK_LOG_ERRNO(SOCK_ERRNO, eLOG_Error, str);
            }
            SOCK_CLOSE(x_sock);
            return eIO_Unknown; /* unrecoverable error */
        }

        /* The connect could be interrupted by a signal or just cannot be
         * established immediately;  yet, the connect must have been in progess
         * (asynchroneous), so wait here for it to succeed (become writeable).
         */
        {{
            struct timeval tv;
            EIO_Status status =
                s_Select(x_sock, eIO_Write,s_to2tv(timeout, &tv));
            if (status != eIO_Success  &&  s_Log) {
                char str[256];
                sprintf(str,
                        "[SOCK::s_Connect]  Failed pending connect to "
                        "%.64s:%d (%.32s)",
                        host ? host : "???",
                        (int) ntohs(x_port),
                        IO_StatusStr(status));
                SOCK_LOG(eLOG_Error, str);
                SOCK_CLOSE(x_sock);
                return status;
            }
        }}
    }

    /* Success */
    /* NOTE:  it does not change the timeouts and the data buffer content */
    sock->sock   = x_sock;
    sock->host   = x_host;
    sock->port   = x_port;
    sock->is_eof = 0/*false*/;

    return eIO_Success;
}


/* Shutdown the socket (close its system file descriptor)
 */
static EIO_Status s_Shutdown(SOCK sock)
{
    /* set EOF flag */
    sock->is_eof = 1/*true*/;

  /* Just checking */
    if (sock->sock == SOCK_INVALID) {
        SOCK_LOG(eLOG_Warning,
                 "[SOCK::s_Shutdown]  Attempt to shutdown an invalid socket");
        return eIO_Unknown;
    }

    {{ /* Reset the auxiliary data buffer */
        size_t buf_size = BUF_Size(sock->buf);
        if (BUF_Read(sock->buf, 0, buf_size) != buf_size) {
            SOCK_LOG(eLOG_Error,
                     "[SOCK::s_Shutdown]  Cannot reset aux. data buffer");
            return eIO_Unknown;
        }
    }}

    /* Set the socket back to blocking mode */
    if ( !s_SetNonblock(sock->sock, 0/*false*/) ) {
        SOCK_LOG(eLOG_Error,
                 "[SOCK::s_Shutdown]  Cannot set socket to blocking mode");
    }

    /* Set the close()'s linger period be equal to the close timeout */
    if ( sock->c_timeout ) {
        struct linger lgr;
        lgr.l_onoff  = 1;
        lgr.l_linger = sock->c_timeout->tv_sec ? sock->c_timeout->tv_sec : 1;
        if (setsockopt(sock->sock, SOL_SOCKET, SO_LINGER, 
                       (char*) &lgr, sizeof(lgr)) != 0) {
            SOCK_LOG_ERRNO(SOCK_ERRNO, eLOG_Warning,
                           "[SOCK::s_Shutdown]  Failed setsockopt()");
        }
    }   

    /* Close (persistently -- retry if interrupted by a signal) */
    for (;;) {
        /* success */
        if (SOCK_CLOSE(sock->sock) == 0)
            break;

        /* error */
        if (SOCK_ERRNO != SOCK_EINTR) {
            SOCK_LOG_ERRNO(SOCK_ERRNO, eLOG_Warning,
                           "[SOCK::s_Shutdown]  Failed close()");
            sock->sock = SOCK_INVALID;
            return (SOCK_ERRNO == SOCK_ECONNRESET || SOCK_ERRNO == SOCK_EPIPE)
                ? eIO_Closed : eIO_Unknown;
        }
    }

    /* Success */
    sock->sock = SOCK_INVALID;
    return eIO_Success;
}


/* Emulate "peek" using the NCBI data buffering.
 * (MSG_PEEK is not implemented on Mac, and it is poorly implemented
 * on Win32, so we had to implement this feature by ourselves)
 */
static int s_NCBI_Recv(SOCK          sock,
                       void*         buffer,
                       size_t        size,
                       int  /*bool*/ peek,
                       int* /*bool*/ read_error)
{
    char*  x_buffer = (char*) buffer;
    size_t n_readbuf;
    int    n_readsock;

    *read_error = 0/*false*/;

    /* read(or peek) from the internal buffer */
    n_readbuf = peek ?
        BUF_Peek(sock->buf, x_buffer, size) :
        BUF_Read(sock->buf, x_buffer, size);
    if (n_readbuf == size)
        return (int)n_readbuf;
    x_buffer += n_readbuf;
    size     -= n_readbuf;

    /* read(dont peek) from the socket */
    n_readsock = recv(sock->sock, x_buffer, size, 0);
    if (n_readsock <= 0) {
        *read_error = 1/*true*/;
        return (int)(n_readbuf ? n_readbuf : n_readsock);
    }

    /* if "peek" -- store the new read data in the internal buffer */
    if ( peek )
        verify(BUF_Write(&sock->buf, x_buffer, n_readsock));

    return (int)(n_readbuf + n_readsock);
}


/* Read/Peek data from the socket
 */
static EIO_Status s_Recv(SOCK        sock,
                         void*       buf,
                         size_t      size,
                         size_t*     n_read,
                         int/*bool*/ peek)
{
    int/*bool*/ read_error;
    int         x_errno;

    /* just checking */
    if (sock->sock == SOCK_INVALID) {
        SOCK_LOG(eLOG_Error,
                 "[SOCK::s_Recv]  Attempted to read from an invalid socket");
        assert(0);
        return eIO_Unknown;
    }

    *n_read = 0;
    for (;;) {
        /* try to read */
        int buf_read = s_NCBI_Recv(sock, buf, size, peek, &read_error);
        if (buf_read > 0) {
            assert(buf_read <= (int)size);
            *n_read = buf_read;
            if ( read_error ) {
                x_errno = SOCK_ERRNO;
                if (x_errno != SOCK_EWOULDBLOCK  &&  x_errno != SOCK_EAGAIN  &&
                    x_errno != SOCK_EINTR)
                    sock->is_eof = 1/*true*/;
            }
            return eIO_Success; /* success */
        }

        x_errno = SOCK_ERRNO;

        if (buf_read == 0) {
            /* NOTE: an empty message may cause the same effect, and it does */
            /* not (!) get discarded from the input queue; therefore, the    */
            /* subsequent attempts to read will cause just the same effect)  */
            sock->is_eof = 1/*true*/;
            return eIO_Closed;  /* the connection has been closed by peer */
        }

        /* blocked -- wait for a data to come;  exit if timeout/error */
        if (x_errno == SOCK_EWOULDBLOCK  ||  x_errno == SOCK_EAGAIN) {
            EIO_Status status =
                s_Select(sock->sock, eIO_Read, sock->r_timeout);
            if (status != eIO_Success)
                return status;
            continue;
        }

        /* retry if interrupted by a signal */
        if (x_errno == SOCK_EINTR)
            continue;

        /* forcibly closed by peer, or shut down */
        if (x_errno == SOCK_ECONNRESET  ||  x_errno == SOCK_EPIPE) {
            sock->is_eof = 1/*true*/;
            return eIO_Closed;
        }

#if defined(NCBI_OS_MAC)
        if (buf_read == -1) {
            sock->is_eof = 1/*true*/;
            return eIO_Closed;
        }
#endif

        /* dont want to handle all possible errors... let them be "unknown" */
        /* and assume EOF */
        sock->is_eof = 1/*true*/;
        break;
    }

    return eIO_Unknown;
}


/* Write data to the socket "as is" (the whole buffer at once)
 */
static EIO_Status s_WriteWhole(SOCK        sock,
                               const void* buf,
                               size_t      size,
                               size_t*     n_written)
{
    const char* x_buf  = (const char*) buf;
    int         x_size = (int) size;
    int         x_errno;

    if ( n_written )
        *n_written = 0;
    if ( !size )
        return eIO_Success;

    for (;;) {
        /* try to write */
        int buf_written = send(sock->sock, (char*)x_buf, x_size, 0);
        if (buf_written >= 0) {
            if ( n_written )
                *n_written += buf_written;
            if (buf_written == x_size)
                return eIO_Success; /* all data has been successfully sent */
            x_buf  += buf_written;
            x_size -= buf_written;
            assert(x_size > 0);
            continue; /* there is unsent data */
        }
        x_errno = SOCK_ERRNO;

        /* blocked -- retry if unblocked before the timeout is expired */
        if (x_errno == SOCK_EWOULDBLOCK  ||  x_errno == SOCK_EAGAIN) {
            EIO_Status status =
                s_Select(sock->sock, eIO_Write, sock->w_timeout);
            if (status != eIO_Success)
                return status;
            continue;
        }

        /* retry if interrupted by a signal */
        if (x_errno == SOCK_EINTR)
            continue;

        /* forcibly closed by peer, or shut down */
        if (x_errno == SOCK_ECONNRESET  ||  x_errno == SOCK_EPIPE)
            return eIO_Closed;

        /* dont want to handle all possible errors... let it be "unknown" */  
        break;
    }
    return eIO_Unknown;
}


#if defined(SOCK_WRITE_SLICE)
/* Split output buffer by slices (of size <= SOCK_WRITE_SLICE) before writing
 * to the socket
 */
static EIO_Status s_WriteSliced(SOCK        sock,
                                const void* buf,
                                size_t      size,
                                size_t*     n_written)
{
    EIO_Status status    = eIO_Success;
    size_t       x_written = 0;

    while (size  &&  status == eIO_Success) {
        size_t n_io = (size > SOCK_WRITE_SLICE) ? SOCK_WRITE_SLICE : size;
        size_t n_io_done;
        status = s_WriteWhole(sock, (char*)buf + x_written, n_io, &n_io_done);
        if ( n_io_done ) {
            x_written += n_io_done;
            size      -= n_io_done;
        }
    }

    if ( n_written )
        *n_written = x_written;
    return status;
}
#endif /* SOCK_WRITE_SLICE */



extern EIO_Status SOCK_Create(const char*     host,
                              unsigned short  port,
                              const STimeout* timeout,
                              SOCK*           sock)
{
    /* Allocate memory for the internal socket structure */
    SOCK x_sock = (SOCK_struct*) calloc(1, sizeof(SOCK_struct));

  /* Connect */
    EIO_Status status = s_Connect(x_sock, host, port, timeout);
    if (status != eIO_Success) {
        free(x_sock);
        return status;
    }

    /* Setup the input data buffer properties */
    BUF_SetChunkSize(&x_sock->buf, SOCK_BUF_CHUNK_SIZE);

    /* Success */
    *sock = x_sock;
    return eIO_Success;
}


extern EIO_Status SOCK_Reconnect(SOCK            sock,
                                 const char*     host,
                                 unsigned short  port,
                                 const STimeout* timeout)
{
    /* Close the socket, if necessary */
    if (sock->sock != SOCK_INVALID)
        s_Shutdown(sock);

  /* Connect */
    return s_Connect(sock, host, port, timeout);
}


extern EIO_Status SOCK_Close(SOCK sock)
{
    EIO_Status status = s_Shutdown(sock);
    BUF_Destroy(sock->buf);
    free(sock);
    return status;
}


extern EIO_Status SOCK_Wait(SOCK            sock,
                            EIO_Event       event,
                            const STimeout* timeout)
{
    struct timeval tv;
    return s_Select(sock->sock, event, s_to2tv(timeout, &tv));
}


extern EIO_Status SOCK_SetTimeout(SOCK            sock,
                                  EIO_Event       event,
                                  const STimeout* timeout)
{
    switch ( event ) {
    case eIO_Read:
        sock->r_timeout = s_to2tv(timeout, &sock->r_tv);
        break;
    case eIO_Write:
        sock->w_timeout = s_to2tv(timeout, &sock->w_tv);
        break;
    case eIO_ReadWrite:
        sock->r_timeout = s_to2tv(timeout, &sock->r_tv);
        sock->w_timeout = s_to2tv(timeout, &sock->w_tv);
        break;
    case eIO_Close:
        sock->c_timeout = s_to2tv(timeout, &sock->c_tv);
        break;
    default:
        assert(0);  return eIO_Unknown;
    }

    return eIO_Success;
}


extern const STimeout* SOCK_GetTimeout(SOCK      sock,
                                       EIO_Event event)
{
    switch ( event ) {
    case eIO_Read:
        return s_tv2to(sock->r_timeout, &sock->r_to);
    case eIO_Write:
        return s_tv2to(sock->w_timeout, &sock->w_to);
    case eIO_Close:
        return s_tv2to(sock->c_timeout, &sock->c_to);
    default:
        assert(0);
    }
    return 0;
}


extern EIO_Status SOCK_Read(SOCK           sock,
                            void*          buf,
                            size_t         size,
                            size_t*        n_read,
                            EIO_ReadMethod how)
{
    switch ( how ) {
    case eIO_Plain: {
        return s_Recv(sock, buf, size, n_read, 0/*false*/);
    }

    case eIO_Peek: {
        return s_Recv(sock, buf, size, n_read, 1/*true*/);
    }

    case eIO_Persist: {
        EIO_Status status = eIO_Success;
        for (*n_read = 0;  size; ) {
            size_t x_read;
            status = SOCK_Read(sock, (char*)buf + *n_read, size, &x_read,
                               eIO_Plain);
            if (status != eIO_Success)
                return status;

            *n_read += x_read;
            size    -= x_read;
        }
        return eIO_Success;
    }
    }

    assert(0);
    return eIO_Unknown;
}


extern EIO_Status SOCK_PushBack(SOCK        sock,
                                const void* buf,
                                size_t      size)
{
    return BUF_PushBack(&sock->buf, buf, size) ? eIO_Success : eIO_Unknown;
}


extern int/*bool*/ SOCK_Eof(SOCK sock)
{
    return sock->is_eof;
}


extern EIO_Status SOCK_Write(SOCK        sock,
                             const void* buf,
                             size_t      size,
                             size_t*     n_written)
{
    /* just checking */
    if (sock->sock == SOCK_INVALID) {
        SOCK_LOG(eLOG_Error,
                 "[SOCK_Write]  Attempted to write to an invalid socket");
        assert(0);
        return eIO_Unknown;
    }

    /* write */
#if defined(SOCK_WRITE_SLICE)
    return s_WriteSliced(sock, buf, size, n_written);
#else
    return s_WriteWhole (sock, buf, size, n_written);
#endif
}


extern void SOCK_GetAddress(SOCK            sock,
                            unsigned int*   host,
                            unsigned short* port,
                            int/*bool*/     network_byte_order) {
    if ( host ) {
        *host = (unsigned int)
            (network_byte_order ? sock->host : ntohl(sock->host));
    }

    if ( port ) {
        *port = (unsigned short)
            (network_byte_order ? sock->port : ntohs(sock->port));
    }
}


extern EIO_Status SOCK_GetOSHandle
(SOCK   sock,
 void*  handle_buf,
 size_t handle_size)
{
    if (handle_size != sizeof(sock->sock)) {
        SOCK_LOG(eLOG_Error, "[SOCK_GetOSHandle]  Invalid handle size!");
        assert(0);
        return eIO_Unknown;
    }

    memcpy(handle_buf, &sock->sock, handle_size);
    return eIO_Success;
}



extern int SOCK_gethostname(char*  name,
                            size_t namelen)
{
    int x_errno;
    assert(namelen > 0);
    verify(s_Initialized  ||  SOCK_InitializeAPI() == eIO_Success);
    name[0] = name[namelen-1] = '\0';
    x_errno = gethostname(name, (int)namelen);
    if (x_errno  ||  name[namelen-1]) {
        SOCK_LOG_ERRNO(x_errno, eLOG_Error,
                       "[SOCK_gethostname]  Cannot get local hostname");
        name[0] = '\0';
        return 1/*failed*/;
    }
    return 0/*success*/;
}


extern int SOCK_host2inaddr(unsigned int host,
                            char*        buf,
                            size_t       buflen)
{
    struct in_addr addr_struct;
    char*          addr_string;

    if ( !buf )
        return 1/*failed*/;

    SOCK_LOCK_WRITE;
    addr_struct.s_addr = host;
    addr_string = inet_ntoa(addr_struct);
    if (!addr_string  ||  strlen(addr_string) >= buflen) {
        char str[96];
        sprintf(str, "[SOCK_host2inaddr]  Cannot convert %x", host);
        SOCK_LOG(eLOG_Error, str);
        buf[0] = '\0';
        return 1/*failed*/;
    }
    strcpy(buf, addr_string);
    SOCK_UNLOCK;

    return 0/*success*/;
}


extern unsigned int SOCK_htonl(unsigned int value)
{
    return htonl(value);
}
