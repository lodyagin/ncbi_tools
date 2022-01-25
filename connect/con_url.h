#ifndef CON_URL__H
#define CON_URL__H

/*  $Id: con_url.h,v 6.4 1999/04/09 22:27:25 vakatov Exp $
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
*  Implement CONNECTOR to hit an URL.
*
*  NOTE: See in "connectr.h" for the detailed specification of the underlying
*        connector("CONNECTOR", "SConnectorTag") methods and structures.
*
* --------------------------------------------------------------------------
* $Log: con_url.h,v $
* Revision 6.4  1999/04/09 22:27:25  vakatov
* Added flag "URLC_URL_ENCODE_ARGS";  thus do not URL-encode CGI args by
* default
*
* Revision 6.3  1999/04/09 21:36:01  vakatov
* Split former CGI "args" into "path?args":
* - added "path" arg. to Ncbi_ConnectURL(), have it URL-encoded automagically;
* - added "path" field to SNetConnInfo
* - split ..._ENGINE_URL defs and config. names to ..._ENGINE_PATH / ARG
*
* Revision 6.2  1999/04/08 18:20:54  vakatov
* Added the on-the-fly URL-encoding/decoding functionality
* (see flags "URLC_URL_*" for URL_CreateConnector())
*
* Revision 6.1  1999/04/01 21:51:22  vakatov
* In comment (2):  must be no '/' before <info->args>!
* ==========================================================================
*/

#include <connectr.h>
#include <connutil.h>


#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

#ifdef __cplusplus
extern "C" {
#endif


/* Create new CONNECTOR structure to hit the specified URL.
 * Use the configuration values recorded in "info", and if "info" is
 * passed NULL then use default info(created by "NetConnInfo_Create(0,0)").
 *
 * In order to workaround some HTTP communication features, this code will:
 *  1) Accumulate all output data in an internal memory buffer until
 *     the first "Read"(or "Peek", or "Close", or "Wait" on read) is performed.
 *  2) On the first "Read"(or "Peek", or "Close", or "Wait" on read), compose
 *     and send the whole HTTP request as:
 *        POST <info->path>?<info->args> HTTP/1.0\n
 *        <user_header\n>
 *        Content-Length: <accumulated_data_length>\n\n
 *        <accumulated_data>
 *     NOTE:
 *       if <user->header> is not NULL/empty string then:
 *       - it must be terminated by a single '\n';
 *       - it gets inserted to the HTTP header "as is", without any
 *         automatic checking or encoding.
 *  3) Now you can "Read" the reply data sent to you by the peer CGI program.
 *  4) Then, if you "Write" again then the connection to the peer CGI program
 *     will be forcibly closed, and you cannot communicate with it anymore
 *     (this CGI process will die).
 *
 *     But if "URLC_AUTO_RECONNECT" is set in "flags" then the connector will
 *     make an automatic reconnect to the same HTTP server with just the
 *     same parameters, and you can repeat the (1,2,3) micro-session with
 *     another instance of your peer CGI program.
 *
 *     If "URLC_AUTO_RECONNECT" is not set then only one
 *     "Write ... Write Read ... Read" micro-session is allowed, and the next
 *     try to "Write" will fail with error status "eCONN_Closed".
 *
 * Setting the "URLC_SURE_FLUSH" flag guarantees that the connector
 * would try to send the HTTP header on "CLOSE" and re-"CONNECT" -- even
 * if no data was written.
 *
 * If the "URLC_HTTP_HEADER" flag is set then dont strip HTTP header
 * (i.e. everything up to the first "\r\n\r\n", including the "\r\n\r\n")
 * from the CGI script response.
 *
 * "URLC_URL_DECODE_INP" -- if set, then:
 *    strip the HTTP header from the input data;  assume the input
 *    data are single-part, URL-encoded; perform the URL-decoding on read.
 *    NOTE:  this flag discards the "URLC_HTTP_HEADER" flag.
 *
 * "URLC_URL_ENCODE_OUT" -- if set, then:
 *    URL-encode all output data(except of the HTTP header).
 *
 * "URLC_URL_ENCODE_ARGS" -- if set, then:
 *    URL-encode "info->args".
 *
 * NOTE: the URL encoding/decoding (in the "URLC_URL_*" cases and "info->args")
 *       is performed by "connutil.[ch]":  URL_Encode() and URL_Decode().
 */

#define URLC_AUTO_RECONNECT  0x1
#define URLC_SURE_FLUSH      0x2
#define URLC_HTTP_HEADER     0x4
#define URLC_URL_DECODE_INP  0x8
#define URLC_URL_ENCODE_OUT  0x10
#define URLC_URL_CODEC       (URLC_URL_DECODE_INP | URLC_URL_ENCODE_OUT)
#define URLC_URL_ENCODE_ARGS 0x20

typedef Nlm_Uint4 URLC_Flags;  /* see the #define URLC_* above */

NLM_EXTERN CONNECTOR URL_CreateConnector
(const SNetConnInfo* info,
 const Nlm_Char*     user_header,
 URLC_Flags          flags
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

#endif /* CON_URL__H */

