#ifndef CONNECT___NCBI_SOCKETP__H
#define CONNECT___NCBI_SOCKETP__H

/* $Id: ncbi_socketp.h,v 1.10 2009/01/22 17:59:35 kazimird Exp $
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
 * Author:  Anton Lavrentiev
 *
 * File Description:
 *   Private API to define socket structure
 *
 */

#include "ncbi_config.h"
/* OS must be specified in the command-line ("-D....") or in the conf. header
 */
#if !defined(NCBI_OS_UNIX) && !defined(NCBI_OS_MSWIN)
#  error "Unknown OS, must be one of NCBI_OS_UNIX, NCBI_OS_MSWIN!"
#endif /*supported platforms*/

#include <connect/ncbi_socket.h>
#include <connect/ncbi_buffer.h>


/* Pull in a minial set of platform-specific system headers here.
 */

#if defined(NCBI_OS_MSWIN)
#  include <winsock2.h>
#elif defined(NCBI_OS_UNIX)
#  include <sys/socket.h>
#  include <sys/time.h>
#endif

/* Portable error codes.
 */
#include <errno.h>

#if   defined(NCBI_OS_MSWIN)

typedef SOCKET TSOCK_Handle;
typedef HANDLE TRIGGER_Handle;

#  define SOCK_EINTR            WSAEINTR
#  define SOCK_EWOULDBLOCK      WSAEWOULDBLOCK/*EAGAIN*/
#  define SOCK_EADDRINUSE       WSAEADDRINUSE
#  define SOCK_ECONNRESET       WSAECONNRESET
#  define SOCK_EPIPE            WSAESHUTDOWN
#  define SOCK_EAGAIN           WSAEINPROGRESS/*special-case missing in WSA*/
#  define SOCK_EINPROGRESS      WSAEINPROGRESS
#  define SOCK_EALREADY         WSAEALREADY
#  define SOCK_ENOTCONN         WSAENOTCONN
#  define SOCK_ECONNABORTED     WSAECONNABORTED
#  define SOCK_ECONNREFUSED     WSAECONNREFUSED
#  define SOCK_ENETRESET        WSAENETRESET
#  define SOCK_ETIMEDOUT        WSAETIMEDOUT
#  define SOCK_SHUTDOWN_RD      SD_RECEIVE
#  define SOCK_SHUTDOWN_WR      SD_SEND
#  define SOCK_SHUTDOWN_RDWR    SD_BOTH

#else

typedef int TSOCK_Handle;
typedef int TRIGGER_Handle;

#  define SOCK_EINTR            EINTR
#  define SOCK_EWOULDBLOCK      EWOULDBLOCK
#  define SOCK_EADDRINUSE       EADDRINUSE
#  define SOCK_ECONNRESET       ECONNRESET
#  define SOCK_EPIPE            EPIPE
#  define SOCK_EAGAIN           EAGAIN
#  define SOCK_EINPROGRESS      EINPROGRESS
#  define SOCK_EALREADY         EALREADY
#  define SOCK_ENOTCONN         ENOTCONN
#  define SOCK_ECONNABORTED     ECONNABORTED
#  define SOCK_ECONNREFUSED     ECONNREFUSED
#  define SOCK_ENETRESET        ENETRESET
#  define SOCK_ETIMEDOUT        ETIMEDOUT

#  ifndef SHUT_RD
#    define SHUT_RD           0
#  endif /*SHUT_RD*/
#  define SOCK_SHUTDOWN_RD    SHUT_RD
#  ifndef SHUT_WR
#    define SHUT_WR           1
#  endif /*SHUT_WR*/
#  define SOCK_SHUTDOWN_WR    SHUT_WR
#  ifndef SHUT_RDWR
#    define SHUT_RDWR         2
#  endif /*SHUT_RDWR*/
#  define SOCK_SHUTDOWN_RDWR  SHUT_RDWR

#endif

#if   defined(ENFILE)
#  define SOCK_ETOOMANY         ENFILE
#elif defined(EMFILE)
#  define SOCK_ETOOMANY         EMFILE
#elif defined(WSAEMFILE)
#  define SOCK_ETOOMANY         WSAEMFILE
#elif defined(EINVAL)
#  define SOCK_ETOOMANY         EINVAL
#else
#  define SOCK_ETOOMANY         0
#endif


#if 0/*defined(__GNUC__)*/
typedef ESwitch    EBSwitch;
typedef EIO_Status EBIO_Status;
typedef ESOCK_Side EBSOCK_Side;
#else
typedef unsigned   EBSwitch;
typedef unsigned   EBIO_Status;
typedef unsigned   EBSOCK_Side;
#endif


typedef enum {
    eInvalid   = 0,
    eTrigger   = 1,
    eSocket    = 2,
    eDatagram  = 3/*2|1*/,
    eListening = 4
} ESOCK_Type;
typedef unsigned char TSOCK_Type;


/* Event trigger
 */
typedef struct TRIGGER_tag {
    TRIGGER_Handle   fd;        /* OS-specific trigger handle                */
    unsigned int     id;        /* the internal ID (cf. "s_ID_Counter")      */

    volatile int     isset;     /* trigger state (UNIX only, otherwise MBZ)  */
    volatile int     isset_;    /* CAUTION: "isset" pointer protrusion area! */

    /* type, status, EOF, log, read-on-write etc bit-field indicators */
    TSOCK_Type          type;   /* eTrigger                                  */
    EBSwitch             log:2; /* how to log events                         */
    EBSOCK_Side         side:1; /* MBZ                                       */
    unsigned/*bool*/    keep:1; /* MBZ                                       */
    EBSwitch          r_on_w:2; /* MBZ                                       */
    EBSwitch        i_on_sig:2; /* eDefault                                  */
    EBIO_Status     r_status:3; /* MBZ (NB: eIO_Success)                     */
    unsigned/*bool*/     eof:1; /* MBZ                                       */
    EBIO_Status     w_status:3; /* MBZ (NB: eIO_Success)                     */
    unsigned/*bool*/ pending:1; /* MBZ                                       */

    unsigned        reserved:8; /* MBZ                                       */

#ifdef NCBI_OS_UNIX
    int              out;       /* write end of the pipe                     */
#endif /*NCBI_OS_UNIX*/
} TRIGGER_struct;


/* Listening socket [must be in one-2-one binary correspondene with TRIGGER]
 */
typedef struct LSOCK_tag {
    TSOCK_Handle     sock;      /* OS-specific socket handle                 */
    unsigned int     id;        /* the internal ID (see also "s_ID_Counter") */

    unsigned int     n_accept;  /* total number of accepted clients          */
    unsigned short   backlog;   /* (unused)                                  */
    unsigned short   port;      /* port on which listening (host byte order) */

    /* type, status, EOF, log, read-on-write etc bit-field indicators */
    TSOCK_Type          type;   /* eListening                                */
    EBSwitch             log:2; /* how to log events and data for this socket*/
    EBSOCK_Side         side:1; /* MBZ                                       */
    unsigned/*bool*/    keep:1; /* MBZ                                       */
    EBSwitch          r_on_w:2; /* MBZ                                       */
    EBSwitch        i_on_sig:2; /* eDefault                                  */
    EBIO_Status     r_status:3; /* MBZ (NB: eIO_Success)                     */
    unsigned/*bool*/     eof:1; /* MBZ                                       */
    EBIO_Status     w_status:3; /* MBZ (NB: eIO_Success)                     */
    unsigned/*bool*/ pending:1; /* MBZ                                       */

    unsigned          unused:1; /* MBZ                                       */
#ifdef NCBI_OS_MSWIN
    unsigned        readable:1; /* =1 if known to have a pending accept      */
    unsigned        reserved:6; /* MBZ                                       */
#else
    unsigned        reserved:7; /* MBZ                                       */
#endif /*NCBI_OS_MSWIN*/

#ifdef NCBI_OS_MSWIN
	WSAEVENT         event;     /* event bound to I/O                        */
#endif /*NCBI_OS_MSWIN*/

    void*            context;   /* per-server credentials                    */

#ifdef NCBI_OS_UNIX
    char             path[1];   /* must go last                              */
#endif /*NCBI_OS_UNIX*/
} LSOCK_struct;


/* Sides of connecting socket
 */
typedef enum {
    eSOCK_Client = 0,
    eSOCK_Server = 1
} ESOCK_Side;


/* Socket [it must be in 1-2-1 binary correspondence with LSOCK above]
 */
typedef struct SOCK_tag {
    TSOCK_Handle     sock;      /* OS-specific socket handle                 */
    unsigned int     id;        /* the internal ID (see also "s_ID_Counter") */

    /* connection point */
    unsigned int     host;      /* peer host (network byte order)            */
    unsigned short   port;      /* peer port (host byte order)               */
    unsigned short   myport;    /* this socket's port number, host byte order*/

    /* type, status, EOF, log, read-on-write etc bit-field indicators */
    TSOCK_Type          type;   /* |= eSocket ({ eSocket | eDatagram })      */
    EBSwitch             log:2; /* how to log events and data for this socket*/
    EBSOCK_Side         side:1; /* socket side: client- or server-side       */
    unsigned/*bool*/    keep:1; /* whether to keep OS handle upon close      */
    EBSwitch          r_on_w:2; /* enable/disable automatic read-on-write    */
    EBSwitch        i_on_sig:2; /* enable/disable I/O restart on signals     */

    EBIO_Status     r_status:3; /* read  status:  eIO_Closed if was shut down*/
    unsigned/*bool*/     eof:1; /* Stream sockets: 'End of file' seen on read
                                   Datagram socks: 'End of message' written  */
    EBIO_Status     w_status:3; /* write status:  eIO_Closed if was shut down*/
    unsigned/*bool*/ pending:1; /* =1 if connection is still initing         */

    unsigned       connected:1; /* =1 if remote end-point is fully connected */
#ifdef NCBI_OS_MSWIN
    unsigned        readable:1; /* =1 if known to be readable                */
    unsigned        closeing:1; /* =1 if FD_CLOSE posted (as ugly as spelled)*/
    unsigned        writable:1; /* =1 if known to be writeable               */
    unsigned        reserved:4; /* MBZ                                       */
#else
    unsigned        reserved:6; /* MBZ                                       */
    unsigned       crossexec:1; /* =1 if close-on-exec must NOT be set       */
#endif /*NCBI_OS_MSWIN*/

#ifdef NCBI_OS_MSWIN
	WSAEVENT         event;     /* event bound to I/O                        */
#endif /*NCBI_OS_MSWIN*/

    void*            session;   /* secure session id if secure, else 0       */

	/* timeouts */
    const struct timeval* r_timeout;/* NULL if infinite, or points to "r_tv" */
    struct timeval   r_tv;      /* finite read  timeout value                */
    STimeout         r_to;      /* finite read  timeout value (aux., temp.)  */
    const struct timeval* w_timeout;/* NULL if infinite, or points to "w_tv" */
    struct timeval   w_tv;      /* finite write timeout value                */
    STimeout         w_to;      /* finite write timeout value (aux., temp.)  */
    const struct timeval* c_timeout;/* NULL if infinite, or points to "c_tv" */
    struct timeval   c_tv;      /* finite close timeout value                */
    STimeout         c_to;      /* finite close timeout value (aux., temp.)  */

    /* aux I/O data */
    BUF              r_buf;     /* read  buffer                              */
    BUF              w_buf;     /* write buffer                              */
    size_t           w_len;     /* SOCK: how much data is pending for output */

    /* statistics */
    size_t           n_read;    /* DSOCK: total #; SOCK: last connect/ only  */
    size_t           n_written; /* DSOCK: total #; SOCK: last /session only  */
    size_t           n_in;      /* DSOCK: msg #; SOCK: total # of bytes read */
    size_t           n_out;     /* DSOCK: msg #; SOCK: total # of bytes sent */

#ifdef NCBI_OS_UNIX
    /* pathname for UNIX socket */
    char             path[1];   /* must go last                              */
#endif /*NCBI_OS_UNIX*/
} SOCK_struct;

/*
 * The following implementation details are worth noting:
 *
 * 1. w_buf is used for stream sockets to keep initial data segment
 *    that has to be sent upon connection establishment.
 *
 * 2. eof is used differently for stream and datagram sockets:
 *    =1 for stream sockets means that read has hit EOF;
 *    =1 for datagram sockets means that message in w_buf has been completed.
 *
 * 3. r_status keeps completion code of the last low-level read call;
 *    however, eIO_Closed is there when the socket is shut down for reading;
 *    see the table below for full details on stream sockets.
 *
 * 4. w_status keeps completion code of the last low-level write call;
 *    however, eIO_Closed is there when the socket is shut down for writing.
 *
 * 5. The following table depicts r_status and eof combinations and their
 *    meanings for stream sockets:
 * -------------------------------+--------------------------------------------
 *              Field             |
 * ---------------+---------------+                  Meaning
 * sock->r_status |   sock->eof   |           (stream sockets only)
 * ---------------+---------------+--------------------------------------------
 * eIO_Closed     |       0       |  Socket shut down for reading
 * eIO_Closed     |       1       |  Read severely failed
 * not eIO_Closed |       0       |  Read completed with r_status error
 * not eIO_Closed |       1       |  Read hit EOF (and [maybe later] r_status)
 * ---------------+---------------+--------------------------------------------
 */


#endif /* CONNECT___NCBI_SOCKETP__H */
