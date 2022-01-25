/*  $Id: ncbi_http_connector.c,v 6.36 2002/11/19 19:20:37 lavr Exp $
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
 * Author:  Denis Vakatov, Anton Lavrentiev
 *
 * File Description:
 *   Implement CONNECTOR for the HTTP-based network connection
 *
 *   See in "ncbi_connector.h" for the detailed specification of the underlying
 *   connector("CONNECTOR", "SConnectorTag") methods and structures.
 *
 */

#include "ncbi_ansi_ext.h"
#include "ncbi_priv.h"
#include <connect/ncbi_http_connector.h>
#include <connect/ncbi_socket.h>
#include <connect/ncbi_buffer.h>
#include <ctype.h>
#include <stdlib.h>


/***********************************************************************
 *  INTERNAL -- Auxiliary types and static functions
 ***********************************************************************/

/* If the connector is allowed to connect
 */
typedef enum {
    eCC_None,
    eCC_Once,
    eCC_Unlimited
} ECanConnect;


/* All internal data necessary to perform the (re)connect and I/O
 *
 * The following states are defined:
 *  "sock"  | "first_read" | State description
 * ---------+--------------+--------------------------------------------------
 *   NULL   |  <whatever>  | User "WRITE" mode: accumulate data in buffer
 * non-NULL |   non-zero   | HTTP header is being read
 * non-NULL |     zero     | HTTP body is being read (user "READ" mode)
 * ---------+--------------+--------------------------------------------------
 */
typedef struct {
    THCC_Flags           flags;           /* as passed to constructor        */
    SConnNetInfo*        net_info;        /* network configuration parameters*/
    FHttpParseHTTPHeader parse_http_hdr;  /* callback to parse HTTP reply hdr*/
    FHttpAdjustNetInfo   adjust_net_info; /* for on-the-fly net_info adjust  */
    FHttpAdjustCleanup   adjust_cleanup;  /* supplemental user data...       */
    void*                adjust_data;     /* ...and cleanup routine          */

    SOCK            sock;         /* socket;  NULL if not in the "READ" mode */
    BUF             obuf;         /* storage to accumulate output data       */
    BUF             ibuf;         /* storage to accumulate input data        */
    BUF             http;         /* storage for HTTP reply header           */
    const STimeout* c_timeout;    /* NULL(infinite), default or ptr to next  */
    STimeout        cc_timeout;   /* storage for a (finite) connect timeout  */
    const STimeout* w_timeout;    /* NULL(infinite), default or ptr to next  */
    STimeout        ww_timeout;   /* storage for a (finite) write timeout    */
    int/*bool*/     first_read;   /* if the 1st attempt to read after connect*/
    ECanConnect     can_connect;  /* if still permitted to make a connection */
    unsigned int    failure_count;/* incremented each failure since open     */
} SHttpConnector;


/* Try to fix connection parameters (called for an unconnected connector) */
static int/*bool*/ s_Adjust(SHttpConnector* uuu,
                            char**          redirect,
                            int/*bool*/     drop_unread)
{
    assert(!uuu->sock);
    /* we're here because something is going wrong */
    if (++uuu->failure_count >= uuu->net_info->max_try) {
        if (*redirect) {
            free(*redirect);
            *redirect = 0;
        }
        if (!drop_unread)
            CORE_LOG(eLOG_Error,"[HTTP]  Too many failed attempts, giving up");
        uuu->can_connect = eCC_None;
        return 0/*failure*/;
    }
    /* adjust info before another connection attempt */
    if (*redirect) {
        int status = ConnNetInfo_ParseURL(uuu->net_info, *redirect);
        free(*redirect);
        *redirect = 0;
        if (!status) {
            CORE_LOG(eLOG_Error, "[HTTP]  Unable to parse redirect");
            uuu->can_connect = eCC_None;
            return 0/*failure*/;
        }
    } else if (!uuu->adjust_net_info ||
               !uuu->adjust_net_info(uuu->net_info,
                                     uuu->adjust_data,
                                     uuu->failure_count)) {
        if (!drop_unread)
            CORE_LOG(eLOG_Error,"[HTTP]  Retry attempts exhausted, giving up");
        uuu->can_connect = eCC_None;
        return 0/*failure*/;
    }

    ConnNetInfo_AdjustForHttpProxy(uuu->net_info);

    if (uuu->net_info->debug_printout)
        ConnNetInfo_Log(uuu->net_info, CORE_GetLOG());
    return 1/*success*/;
}


/* Unconditionally drop the connection; timeout may specify time allowance */
static void s_DropConnection(SHttpConnector* uuu, const STimeout* timeout)
{
    size_t http_size = BUF_Size(uuu->http);
    assert(uuu->sock);
    if (BUF_Read(uuu->http, 0, http_size) != http_size) {
        CORE_LOG(eLOG_Error, "[HTTP]  Cannot drop HTTP header buffer");
        assert(0);
    }
    SOCK_SetTimeout(uuu->sock, eIO_Close, timeout);
    SOCK_Close(uuu->sock);
    uuu->sock = 0;
}


/* Connect to the HTTP server, specified by uuu->net_info's "port:host".
 * Return eIO_Success only if socket connection has succeeded and uuu->sock
 * is non-zero. If unsuccessful, try to adjust uuu->net_info by s_Adjust(),
 * and then re-try the connection attempt.
 */
static EIO_Status s_Connect(SHttpConnector* uuu, int/*bool*/ drop_unread)
{
    assert(!uuu->sock);
    if (uuu->can_connect == eCC_None) {
        CORE_LOG(eLOG_Error, "[HTTP]  Connector is no longer usable");
        return eIO_Closed;
    }

    /* the re-try loop... */
    for (;;) {
        char* null = 0;

        /* connect & send HTTP header */
        uuu->sock = URL_Connect
            (uuu->net_info->host, uuu->net_info->port,
             uuu->net_info->path, uuu->net_info->args,
             uuu->net_info->req_method, BUF_Size(uuu->obuf),
             uuu->c_timeout, uuu->w_timeout,
             uuu->net_info->http_user_header,
             (int/*bool*/) (uuu->flags & fHCC_UrlEncodeArgs),
             uuu->net_info->debug_printout == eDebugPrintout_Data ? eOn :
             uuu->net_info->debug_printout == eDebugPrintout_None ? eOff :
             eDefault);

        if (uuu->sock) {
            if (!(uuu->flags & fHCC_NoUpread))
                SOCK_SetReadOnWrite(uuu->sock, eOn);
            return eIO_Success;
        }

        /* connection failed, no socket was created */
        if (!s_Adjust(uuu, &null, drop_unread))
            break;
    }

    return eIO_Unknown;
}


/* Connect to the server specified by uuu->net_info, then compose and form
 * relevant HTTP header, and flush the accumulated output data(uuu->obuf)
 * after the HTTP header. If connection/write unsiccessful, retry to reconnect
 * and send the data again until permitted by s_Adjust().
 */
static EIO_Status s_ConnectAndSend(SHttpConnector* uuu,int/*bool*/ drop_unread)
{
    EIO_Status status;

    assert(!uuu->sock);
    uuu->first_read = 1/*true*/;
    for (;;) {
        char* null = 0;
        size_t size;

        if ((status = s_Connect(uuu, drop_unread)) != eIO_Success)
            break;
        assert(uuu->sock);

        if ((size = BUF_Size(uuu->obuf)) != 0) {
            /* flush the accumulated output data */
            size_t off = 0;

            SOCK_SetTimeout(uuu->sock, eIO_Write, uuu->w_timeout);
            do {
                char   buf[4096];
                size_t n_written;
                size_t len = BUF_PeekAt(uuu->obuf, off, buf, sizeof(buf));

                status = SOCK_Write(uuu->sock, buf, len,
                                    &n_written, eIO_WritePersist);
                if (status != eIO_Success) {
                    CORE_LOGF(eLOG_Error,  ("[HTTP]  Error writing body at "
                                            "offset %lu", (unsigned long)off));
                    break;
                }
                off += n_written;
            } while (off < size);
        }

        if (status == eIO_Success) {
            /* shutdown the socket for writing */
#if !defined(NCBI_OS_MSWIN)
            /* (on MS-Win, socket shutdown on writing apparently messes up (?!)
             *  with the later reading, esp. when reading a lot of data.) */
            SOCK_Shutdown(uuu->sock, eIO_Write);
#endif
            break;
        }

        /* write failed; close socket and try to use another server */
        s_DropConnection(uuu, 0/*no wait*/);
        if (!s_Adjust(uuu, &null, drop_unread))
            break;
    }

    return status;
}


/* Parse HTTP header */
static EIO_Status s_ReadHeader(SHttpConnector* uuu, char** redirect)
{
    int/*bool*/ moved = 0/*false*/;
    int         server_error = 0;
    int         http_status = 0;
    EIO_Status  status;
    char*       header;
    size_t      size;

    assert(uuu->sock && uuu->first_read);
    *redirect = 0;
    if (uuu->flags & fHCC_KeepHeader) {
        uuu->first_read = 0;
        return eIO_Success;
    }

    /* line by line input */
    for (;;) {
        status = SOCK_StripToPattern(uuu->sock, "\r\n", 2, &uuu->http, 0);
        if (status != eIO_Success) {
            const STimeout* tmo = SOCK_GetTimeout(uuu->sock, eIO_Read);
            if (tmo->sec || tmo->usec)
                CORE_LOG(eLOG_Error, "[HTTP]  Error reading header");
            return status;
        }

        size = BUF_Size(uuu->http);
        if (!(header = (char*) malloc(size + 1))) {
            CORE_LOGF(eLOG_Error, ("[HTTP]  Cannot allocate header, %lu bytes",
                                   (unsigned long) size));
            return eIO_Unknown;
        }
        verify(BUF_Peek(uuu->http, header, size) == size);
        header[size] = '\0';
        if (strcmp(&header[size - 4], "\r\n\r\n") == 0)
            break;
        free(header);
    }
    uuu->first_read = 0/*false*/; /* the entire header has been read */
    if (BUF_Read(uuu->http, 0, size) != size)
        CORE_LOG(eLOG_Warning, "[HTTP]  Cannot drop HTTP header buffer");

    /* HTTP status must come on the first line of the reply */
    if (sscanf(header, " HTTP/%*d.%*d %d ", &http_status) != 1  ||
        http_status < 200  ||  299 < http_status) {
        server_error = http_status;
        if (http_status == 301  ||  http_status == 302)
            moved = 1;
    }

    if (uuu->net_info->debug_printout) {
        CORE_DATA(header, size,
                  server_error ? "HTTP header (server error)" : "HTTP header");
    }

    if (uuu->parse_http_hdr) {
        if (!(*uuu->parse_http_hdr)
            (header, uuu->adjust_data, server_error))
            server_error = 1;
    }

    if (moved) {
        const char k_LocationTag[] = "\nLocation: ";
        char*      location = strstr(header, k_LocationTag);

        if (location) {
            char* s;

            location += sizeof(k_LocationTag) - 1;
            if (!(s = strchr(location, '\r')))
                *strchr(location, '\n') = 0;
            else
                *s = 0;
            while (*location) {
                if (isspace((unsigned char)(*location)))
                    location++;
                else
                    break;
            }
            for (s = location; *s; s++)
                if (isspace((unsigned char)(*s)))
                    break;
            *s = 0;
            if (*location)
                *redirect = strdup(location);
        }
    }

    if (header)
        free(header);

    /* skip & printout the content, if server error was flagged */
    if (server_error  &&  uuu->net_info->debug_printout) {
        size_t size;
        char*  body;
        BUF    buf = 0;
        SOCK_SetTimeout(uuu->sock, eIO_Read, 0);
        status = SOCK_StripToPattern(uuu->sock, 0, 0, &buf, 0);
        assert(status != eIO_Success); /* cause reading until EOF */
        if (!(size = BUF_Size(buf)))
            CORE_LOG(eLOG_Trace, "[HTTP]  No body received along with error");
        else if ((body = (char*) malloc(size)) != 0) {
            if (BUF_Read(buf, body, size) != size) {
                CORE_LOG(eLOG_Error, "[HTTP]  Cannot read HTTP body buffer");
                assert(0);
            }
            CORE_DATA(body, size, "HTTP server error");
            free(body);
        } else
            CORE_LOGF(eLOG_Error, ("[HTTP]  Cannot allocate body, %lu bytes",
                                   (unsigned long) size));
        BUF_Destroy(buf);
    }

    if (server_error)
        return eIO_Unknown;
    return status;
}


/* Prepare connector for reading. Open socket if necessary and
 * make initial connect and send, re-trying if possible until success.
 * Return codes:
 *   eIO_Success = success, connector is ready for reading (uuu->sock != NULL);
 *   eIO_Timeout = maybe (check uuu->sock) connected and no data yet available;
 *   other code  = error, not connected (uuu->sock == NULL).
 */
static EIO_Status s_PreRead(SHttpConnector* uuu,
                            const STimeout* timeout,
                            int/*bool*/     drop_unread)
{
    static const STimeout zero_timeout = {0, 0};
    char* redirect = 0;
    EIO_Status status;

    for (;;) {
        if (!uuu->sock) {
            if ((status = s_ConnectAndSend(uuu, drop_unread)) != eIO_Success)
                break;
            assert(uuu->sock);
        } else
            status = eIO_Success;

        /* set timeout */
        SOCK_SetTimeout(uuu->sock, eIO_Read, timeout);

        if (!uuu->first_read)
            break;

        if ((status = s_ReadHeader(uuu, &redirect)) == eIO_Success) {
            size_t obuf_size = BUF_Size(uuu->obuf);
            assert(!uuu->first_read);
            /* pending output data no longer needed */
            if (BUF_Read(uuu->obuf, 0, obuf_size) != obuf_size)
                CORE_LOG(eLOG_Warning, "[HTTP]  Cannot drop output buffer");
            break;
        }
        /* if polling then bail out with eIO_Timeout */
        if (status == eIO_Timeout && timeout &&
            memcmp(timeout, &zero_timeout, sizeof(STimeout)) == 0)
            break;

        /* HTTP header read error; disconnect and try to use another server */
        s_DropConnection(uuu, 0/*no wait*/);
        if (!s_Adjust(uuu, &redirect, drop_unread))
            break;
        assert(redirect == 0);
    }
    assert(redirect == 0);
    return status;
}


/* Read non-header data from connection */
static EIO_Status s_Read(SHttpConnector* uuu, void* buf,
                         size_t size, size_t* n_read)
{
    assert(uuu->sock);
    /* just read, with no URL-decoding */
    if (!(uuu->flags & fHCC_UrlDecodeInput))
        return SOCK_Read(uuu->sock, buf, size, n_read, eIO_ReadPlain);

    /* read and URL-decode */
    {{
        EIO_Status status;
        size_t     n_peeked, n_decoded;
        size_t     peek_size = 3 * size;
        void*      peek_buf  = malloc(peek_size);

        /* peek the data */
        status= SOCK_Read(uuu->sock,peek_buf,peek_size,&n_peeked,eIO_ReadPeek);
        if (status != eIO_Success) {
            *n_read = 0;
            free(peek_buf);
            return status;
        }

        /* decode, then discard the successfully decoded data from the input */
        if (URL_Decode(peek_buf, n_peeked, &n_decoded, buf, size, n_read)) {
            if (n_decoded) {
                size_t x_read;
                SOCK_Read(uuu->sock,peek_buf,n_decoded,&x_read,eIO_ReadPlain);
                assert(x_read == n_decoded);
                status = eIO_Success;
            } else if (SOCK_Status(uuu->sock, eIO_Read) == eIO_Closed)
                /* we are at EOF, and the remaining data cannot be decoded */
                status = eIO_Unknown;
        } else
            status = eIO_Unknown;
        
        if (status != eIO_Success)
            CORE_LOG(eLOG_Error, "[HTTP]  Cannot URL-decode data");

        free(peek_buf);
        return status;
    }}
}


/* Reset/readout input data and close socket */
static EIO_Status s_Disconnect(SHttpConnector* uuu,
                               int/*bool*/     drop_unread,
                               const STimeout* timeout)
{
    EIO_Status status = eIO_Success;

    if (drop_unread) {
        size_t ibuf_size = BUF_Size(uuu->ibuf);
        if (BUF_Read(uuu->ibuf, 0, ibuf_size) != ibuf_size) {
            CORE_LOG(eLOG_Error, "[HTTP]  Cannot drop input buffer");
            assert(0);
        }
    } else if ((status = s_PreRead(uuu, timeout, 0/*nodrop*/)) == eIO_Success){
        do {
            char     buf[4096];
            size_t   x_read;
            status = s_Read(uuu, buf, sizeof(buf), &x_read);
            if (!BUF_Write(&uuu->ibuf, buf, x_read))
                status = eIO_Unknown;
        } while (status == eIO_Success);
        if (status == eIO_Closed)
            status = eIO_Success;
    }

    if (uuu->sock) /* s_PreRead() might have dropped the connection already */
        s_DropConnection(uuu, timeout);
    if (uuu->can_connect == eCC_Once)
        uuu->can_connect = eCC_None;

    return status;
}


/* Send the accumulated output data(if any) to server, then close socket.
 * Regardless of the flush, clear both input and output buffer.
 * This function is only called to either re-open or close the connector.
 */
static void s_FlushAndDisconnect(SHttpConnector* uuu,
                                 const STimeout* timeout,
                                 int/*bool*/     close)
{
    size_t obuf_size = BUF_Size(uuu->obuf);

    /* store timeouts for later use */
    if (timeout) {
        uuu->cc_timeout = *timeout;
        uuu->c_timeout  = &uuu->cc_timeout;
        uuu->ww_timeout = *timeout;
        uuu->w_timeout  = &uuu->ww_timeout;
    } else {
        uuu->c_timeout  = timeout;
        uuu->w_timeout  = timeout;
    }

    if (close && !uuu->sock && uuu->can_connect != eCC_None &&
        ((uuu->flags & fHCC_SureFlush) || BUF_Size(uuu->obuf))) {
        /* "WRITE" mode and data (or just flag) pending */
        s_PreRead(uuu, timeout, 1/*drop*/);
    }
    s_Disconnect(uuu, 1/*drop_unread*/, timeout);
    assert(!uuu->sock);

    /* clear pending output data, if any */
    if (BUF_Read(uuu->obuf, 0, obuf_size) != obuf_size) {
        CORE_LOG(eLOG_Error, "[HTTP]  Cannot drop output buffer");
        assert(0);
    }
}


/***********************************************************************
 *  INTERNAL -- "s_VT_*" functions for the "virt. table" of connector methods
 ***********************************************************************/

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    static const char* s_VT_GetType(CONNECTOR       connector);
    static EIO_Status  s_VT_Open   (CONNECTOR       connector,
                                    const STimeout* timeout);
    static EIO_Status  s_VT_Wait   (CONNECTOR       connector,
                                    EIO_Event       event,
                                    const STimeout* timeout);
    static EIO_Status  s_VT_Write  (CONNECTOR       connector,
                                    const void*     buf,
                                    size_t          size,
                                    size_t*         n_written,
                                    const STimeout* timeout);
    static EIO_Status  s_VT_Flush  (CONNECTOR       connector,
                                    const STimeout* timeout);
    static EIO_Status  s_VT_Read   (CONNECTOR       connector,
                                    void*           buf,
                                    size_t          size,
                                    size_t*         n_read,
                                    const STimeout* timeout);
    static EIO_Status  s_VT_Status (CONNECTOR       connector,
                                    EIO_Event       dir);
    static EIO_Status  s_VT_Close  (CONNECTOR       connector,
                                    const STimeout* timeout);

    static void        s_Setup     (SMetaConnector *meta,
                                    CONNECTOR connector);
    static void        s_Destroy   (CONNECTOR connector);
#  ifdef IMPLEMENTED__CONN_WaitAsync
    static EIO_Status  s_VT_WaitAsync(CONNECTOR     connector,
                                      FConnectorAsyncHandler  func,
                                      SConnectorAsyncHandler* data);
#  endif
#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */


static const char* s_VT_GetType
(CONNECTOR connector)
{
    return "HTTP";
}


static EIO_Status s_VT_Open
(CONNECTOR       connector,
 const STimeout* timeout)
{
    SHttpConnector* uuu = (SHttpConnector*) connector->handle;

    /* NOTE: the real connect will be performed on the first "READ", or
     *"CLOSE", or on the "WAIT" on read -- see in "s_ConnectAndSend()";
     * we just close underlying socket and prepare to open it later */
    s_FlushAndDisconnect(uuu, timeout, 0/*open*/);

    /* reset the auto-reconnect feature */
    uuu->can_connect = uuu->flags & fHCC_AutoReconnect
        ? eCC_Unlimited : eCC_Once;
    uuu->failure_count = 0;

    return eIO_Success;
}


static EIO_Status s_VT_Wait
(CONNECTOR       connector,
 EIO_Event       event,
 const STimeout* timeout)
{
    SHttpConnector* uuu = (SHttpConnector*) connector->handle;

    switch (event) {
    case eIO_Read:
        if (uuu->can_connect == eCC_None)
            return eIO_Closed;
        if (!uuu->sock || uuu->first_read) {
            EIO_Status status = s_PreRead(uuu, timeout, 0/*no drop unread*/);
            if (status != eIO_Success || BUF_Size(uuu->ibuf))
                return status;
            assert(uuu->sock);
        }
        return SOCK_Wait(uuu->sock, eIO_Read, timeout);
    case eIO_Write:
        /* Return 'Closed' if no more writes are allowed (and now - reading) */
        return uuu->can_connect == eCC_None ||
            (uuu->sock && uuu->can_connect == eCC_Once)
            ? eIO_Closed : eIO_Success;
    default:
        assert(0);
        return eIO_InvalidArg;
    }
}


static EIO_Status s_VT_Write
(CONNECTOR       connector,
 const void*     buf,
 size_t          size,
 size_t*         n_written,
 const STimeout* timeout)
{
    SHttpConnector* uuu = (SHttpConnector*) connector->handle;

    /* if trying to "WRITE" after "READ" then close the socket,
     * and so switch to "WRITE" mode */
    if (uuu->sock) {
        EIO_Status status = s_Disconnect(uuu,
                                         uuu->flags & fHCC_DropUnread,
                                         timeout);
        if (status != eIO_Success)
            return status;
    }
    if (uuu->can_connect == eCC_None)
        return eIO_Closed; /* no more connects permitted */

    /* accumulate all output in the memory buffer */
    if (uuu->flags & fHCC_UrlEncodeOutput) {
        /* with URL-encoding */
        size_t dst_size = 3 * size;
        void* dst = malloc(dst_size);
        size_t dst_written;
        URL_Encode(buf, size, n_written, dst, dst_size, &dst_written);
        assert(*n_written == size);
        if (!BUF_Write(&uuu->obuf, dst, dst_written)) {
            free(dst);
            return eIO_Unknown;
        }
        free(dst);
    }
    else {
        /* "as is" (without URL-encoding) */
        if (!BUF_Write(&uuu->obuf, buf, size))
            return eIO_Unknown;
        *n_written = size;
    }

    /* store the write timeout */
    if (timeout) {
        uuu->ww_timeout = *timeout;
        uuu->w_timeout  = &uuu->ww_timeout;
    } else
        uuu->w_timeout  = timeout;

    return eIO_Success;
}


static EIO_Status s_VT_Flush
(CONNECTOR       connector,
 const STimeout* timeout)
{
    SHttpConnector* uuu = (SHttpConnector*) connector->handle;

    /* The real flush will be performed on the first "READ" (or "CLOSE"),
     * or on "WAIT". Here, we just store the write timeout, that's all...
     */
    if (timeout) {
        uuu->ww_timeout = *timeout;
        uuu->w_timeout  = &uuu->ww_timeout;
    } else
        uuu->w_timeout  = timeout;

    return eIO_Success;
}


static EIO_Status s_VT_Read
(CONNECTOR       connector,
 void*           buf,
 size_t          size,
 size_t*         n_read,
 const STimeout* timeout)
{
    SHttpConnector* uuu = (SHttpConnector*) connector->handle;
    EIO_Status status = s_PreRead(uuu, timeout, 0/*no drop unread*/);
    size_t x_read = BUF_Read(uuu->ibuf, buf, size);

    *n_read = x_read;
    if (x_read < size) {
        if (status != eIO_Success)
            return status;
        status = s_Read(uuu, (char*) buf + x_read, size - x_read, n_read);
        *n_read += x_read;
    } else
        status = eIO_Success;

    return status;
}


static EIO_Status s_VT_Status
(CONNECTOR connector,
 EIO_Event dir)
{
    SHttpConnector* uuu = (SHttpConnector*) connector->handle;
    return uuu->sock ? SOCK_Status(uuu->sock, dir) : eIO_Success;
}


static EIO_Status s_VT_Close
(CONNECTOR       connector,
 const STimeout* timeout)
{
    s_FlushAndDisconnect((SHttpConnector*) connector->handle,
                         timeout, 1/*close*/);
    return eIO_Success;
}


#ifdef IMPLEMENTED__CONN_WaitAsync
static EIO_Status s_VT_WaitAsync
(void*                   connector,
 FConnectorAsyncHandler  func,
 SConnectorAsyncHandler* data)
{
    return eIO_NotSupported;
}
#endif


static void s_Setup(SMetaConnector *meta, CONNECTOR connector)
{
    SHttpConnector* uuu = (SHttpConnector*) connector->handle;

    /* initialize virtual table */
    CONN_SET_METHOD(meta, get_type,   s_VT_GetType,   connector);
    CONN_SET_METHOD(meta, open,       s_VT_Open,      connector);
    CONN_SET_METHOD(meta, wait,       s_VT_Wait,      connector);
    CONN_SET_METHOD(meta, write,      s_VT_Write,     connector);
    CONN_SET_METHOD(meta, flush,      s_VT_Flush,     connector);
    CONN_SET_METHOD(meta, read,       s_VT_Read,      connector);
    CONN_SET_METHOD(meta, status,     s_VT_Status,    connector);
    CONN_SET_METHOD(meta, close,      s_VT_Close,     connector);
#ifdef IMPLEMENTED__CONN_WaitAsync
    CONN_SET_METHOD(meta, wait_async, s_VT_WaitAsync, connector);
#endif
    CONN_SET_DEFAULT_TIMEOUT(meta, uuu->net_info->timeout);
}


static void s_Destroy(CONNECTOR connector)
{
    SHttpConnector* uuu = (SHttpConnector*) connector->handle;

    ConnNetInfo_Destroy(uuu->net_info);
    if (uuu->adjust_cleanup)
        uuu->adjust_cleanup(uuu->adjust_data);
    BUF_Destroy(uuu->http);
    BUF_Destroy(uuu->obuf);
    BUF_Destroy(uuu->ibuf);
    free(uuu);
    connector->handle = 0;
    free(connector);
}


/***********************************************************************
 *  EXTERNAL -- the connector's "constructor"
 ***********************************************************************/

extern CONNECTOR HTTP_CreateConnector
(const SConnNetInfo* info,
 const char*         user_header,
 THCC_Flags          flags)
{
    SConnNetInfo *net_info;
    CONNECTOR connector;

    net_info = info ? ConnNetInfo_Clone(info) : ConnNetInfo_Create(0);
    if (user_header)
        ConnNetInfo_SetUserHeader(net_info, user_header);
    connector = HTTP_CreateConnectorEx(net_info, flags, 0, 0, 0, 0);
    ConnNetInfo_Destroy(net_info);
    return connector;
}


extern CONNECTOR HTTP_CreateConnectorEx
(const SConnNetInfo*  net_info,
 THCC_Flags           flags,
 FHttpParseHTTPHeader parse_http_hdr,
 FHttpAdjustNetInfo   adjust_net_info,
 void*                adjust_data,
 FHttpAdjustCleanup   adjust_cleanup)
{
    CONNECTOR       ccc = (SConnector    *) malloc(sizeof(SConnector    ));
    SHttpConnector* uuu = (SHttpConnector*) malloc(sizeof(SHttpConnector));

    /* initialize internal data structures */
    uuu->net_info = net_info
        ? ConnNetInfo_Clone(net_info) : ConnNetInfo_Create(0);
    ConnNetInfo_AdjustForHttpProxy(uuu->net_info);
    if (uuu->net_info->debug_printout)
        ConnNetInfo_Log(uuu->net_info, CORE_GetLOG());

    uuu->flags           = flags;
    if (flags & fHCC_UrlDecodeInput)
        uuu->flags      &= ~fHCC_KeepHeader;

    uuu->parse_http_hdr  = parse_http_hdr;

    uuu->adjust_net_info = adjust_net_info;
    uuu->adjust_data     = adjust_data;
    uuu->adjust_cleanup  = adjust_cleanup;

    uuu->sock = 0;
    uuu->obuf = 0;
    uuu->ibuf = 0;
    uuu->http = 0;
    uuu->c_timeout = CONN_DEFAULT_TIMEOUT; /* intentionally bad values --    */
    uuu->w_timeout = CONN_DEFAULT_TIMEOUT; /* should be reset prior to usage */
    uuu->can_connect = eCC_Once;           /* will be properly reset at open */

    /* initialize connector structure */
    ccc->handle  = uuu;
    ccc->next    = 0;
    ccc->meta    = 0;
    ccc->setup   = s_Setup;
    ccc->destroy = s_Destroy;

    return ccc;
}


/*
 * --------------------------------------------------------------------------
 * $Log: ncbi_http_connector.c,v $
 * Revision 6.36  2002/11/19 19:20:37  lavr
 * Server error parsing slightly changed
 *
 * Revision 6.35  2002/10/28 15:46:20  lavr
 * Use "ncbi_ansi_ext.h" privately
 *
 * Revision 6.34  2002/10/22 15:11:24  lavr
 * Zero connector's handle to crash if revisited
 *
 * Revision 6.33  2002/09/06 15:46:45  lavr
 * Few bug fixes about corner cases of I/O waiting, flushing and logging
 *
 * Revision 6.32  2002/08/12 15:12:23  lavr
 * Use persistent SOCK_Write()
 *
 * Revision 6.31  2002/08/07 16:33:04  lavr
 * Changed EIO_ReadMethod enums accordingly; log moved to end
 *
 * Revision 6.30  2002/07/25 20:20:43  lavr
 * Do not report header read error on {0, 0} timeout
 *
 * Revision 6.29  2002/07/25 13:59:35  lavr
 * More diagnostic messages added
 *
 * Revision 6.28  2002/06/19 18:08:02  lavr
 * Fixed some wrong assumptions on use of s_PreRead(); more comments added
 *
 * Revision 6.27  2002/06/10 19:51:20  lavr
 * Small prettifying
 *
 * Revision 6.26  2002/05/06 19:13:48  lavr
 * Output to stderr replaced with calls to logger
 *
 * Revision 6.25  2002/04/26 16:36:56  lavr
 * Added setting of default timeout in meta-connector's setup routine
 * Remove all checks for CONN_DEFAULT_TIMEOUT: now supplied good from CONN
 *
 * Revision 6.24  2002/04/22 19:31:33  lavr
 * Reading/waiting redesigned to be more robust in case of network errors
 *
 * Revision 6.23  2002/03/22 19:52:16  lavr
 * Do not include <stdio.h>: included from ncbi_util.h or ncbi_priv.h
 *
 * Revision 6.22  2002/02/11 20:36:44  lavr
 * Use "ncbi_config.h"
 *
 * Revision 6.21  2002/02/05 22:04:12  lavr
 * Included header files rearranged
 *
 * Revision 6.20  2001/12/31 14:53:46  lavr
 * #include <connect/ncbi_ansi_ext.h> for Mac to be happy (noted by J.Kans)
 *
 * Revision 6.19  2001/12/30 20:00:00  lavr
 * Redirect on non-empty location only
 *
 * Revision 6.18  2001/12/30 19:41:07  lavr
 * Process error codes 301 and 302 (document moved) and reissue HTTP request
 *
 * Revision 6.17  2001/09/28 20:48:23  lavr
 * Comments revised; parameter (and SHttpConnector's) names adjusted
 * Retry logic moved entirely into s_Adjust()
 *
 * Revision 6.16  2001/09/10 21:15:56  lavr
 * Readability issue: FParseHTTPHdr -> FParseHTTPHeader
 *
 * Revision 6.15  2001/05/23 21:52:44  lavr
 * +fHCC_NoUpread
 *
 * Revision 6.14  2001/05/17 15:02:51  lavr
 * Typos corrected
 *
 * Revision 6.13  2001/05/08 20:26:27  lavr
 * Patches in re-try code
 *
 * Revision 6.12  2001/04/26 20:21:34  lavr
 * Reorganized and and made slightly more effective
 *
 * Revision 6.11  2001/04/24 21:41:47  lavr
 * Reorganized code to allow reconnects in case of broken connections at
 * stage of connection, sending data and receiving HTTP header. Added
 * code to pull incoming data from connection while sending - stall protection.
 *
 * Revision 6.10  2001/03/02 20:08:47  lavr
 * Typo fixed
 *
 * Revision 6.9  2001/01/25 16:53:24  lavr
 * New flag for HTTP_CreateConnectorEx: fHCC_DropUnread
 *
 * Revision 6.8  2001/01/23 23:11:20  lavr
 * Status virtual method implemented
 *
 * Revision 6.7  2001/01/11 16:38:17  lavr
 * free(connector) removed from s_Destroy function
 * (now always called from outside, in METACONN_Remove)
 *
 * Revision 6.6  2001/01/03 22:32:43  lavr
 * Redundant calls to 'adjust_info' removed
 *
 * Revision 6.5  2000/12/29 17:57:16  lavr
 * Adapted for use of new connector structure;
 * parse header callback added; some internal functions renamed.
 *
 * Revision 6.4  2000/11/15 18:52:02  vakatov
 * Call SOCK_Shutdown() after the HTTP request is sent.
 * Use SOCK_Status() instead of SOCK_Eof().
 *
 * Revision 6.3  2000/10/12 21:43:14  vakatov
 * Minor cosmetic fix...
 *
 * Revision 6.2  2000/09/26 22:02:57  lavr
 * HTTP request method added
 *
 * Revision 6.1  2000/04/21 19:41:01  vakatov
 * Initial revision
 *
 * ==========================================================================
 */
