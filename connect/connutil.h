#ifndef CONNUTIL__H
#define CONNUTIL__H

/*  $Id: connutil.h,v 6.13 2000/06/29 17:32:21 vakatov Exp $
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
*   Auxiliary and miscellaneous functions to:
*    - get and handle "standard" connection-related configuration
*      from the program config.-file(s) and environment
*    - convert error codes to EConnStatus
*    - seek a pattern in the CONN or SOCK stream
*    - perform an URL encoding/decoding
*    - etc.
*
* --------------------------------------------------------------------------
* $Log: connutil.h,v $
* Revision 6.13  2000/06/29 17:32:21  vakatov
* Added more MIME sub-types
* Allow non-NCBI ("standard") MIME types
* Added MIME_*ContentTypeEx(), and tests for them
*
* Revision 6.12  2000/02/25 16:45:53  vakatov
* Redesigned to really share "ncbi_*.[ch]" etc. between the C and
* the C++ toolkits, and even to use them in a "standalone" fashion
*
* Revision 6.11  1999/09/13 15:54:32  vakatov
* Added URL_DecodeEx() -- for "relaxed" URL decoding: let the user to
* allow some of the symbols prohibited by the standard
*
* Revision 6.10  1999/07/26 17:58:36  vakatov
* Use "\r\n" rather than just "\n" as the HTTP header line terminator.
* Made the hard-coded names and defautl values of config.parameters
* (#define CFG_CONN_, DEF_CONN_) be public.
* +eMIME_Fasta.
*
* Revision 6.9  1999/07/23 13:16:10  vakatov
* MIME_ParseContentType() - dont crash if passed NULL string
*
* Revision 6.8  1999/07/23 00:37:21  vakatov
* Final version of NCBI MIME functions
*
* Revision 6.7  1999/07/21 21:32:32  vakatov
* Get rid of the "nph-" prefixes for DISPD/NCBID -- use new DISPD/NCBID
*
* Revision 6.6  1999/07/16 21:13:42  vakatov
* + NetConnInfo_Print()
* + NCBI MIME-type handling routines and typedefs for a future use (MIME_*)
*
* Revision 6.5  1999/07/09 22:53:38  vakatov
* Added more conf. parameters:  SRV_NCBID_PORT, SRV_NCBID_PATH, SRV_LB_DISABLE
*
* Revision 6.4  1999/06/28 16:28:02  vakatov
* SNetConnInfo:: separated the HTTP and the CERL-like(non-transparent)
* firewall proxies;  renamed config./env. parameters accordingly:
*   SRV_HTTP_PROXY_HOST/PORT replaced the former SRV_PROXY_HOST/PORT
*   SRV_PROXY_HOST now specifies the host name of non-transparent CERN proxy
*   SRV_PROXY_PORT is obsolete
*   SRV_PROXY_TRANSPARENT is obsolete
* Also:  NetConnInfo_AdjustForCernProxy --> ...HttpProxy
*
* Revision 6.3  1999/04/09 22:27:01  vakatov
* Ncbi_ConnectURL():  added "encode_args"
*
* Revision 6.2  1999/04/09 21:36:02  vakatov
* Split former CGI "args" into "path?args":
* - added "path" arg. to Ncbi_ConnectURL(), have it URL-encoded automagically;
* - added "path" field to SNetConnInfo
* - split ..._ENGINE_URL defs and config. names to ..._ENGINE_PATH / ARG
*
* Revision 6.1  1999/04/07 20:43:53  vakatov
* Added URL_Encode() and URL_Decode()
* ==========================================================================
*/

#include <connectn.h>
#include <ncbisock.h>
#include <ncbibuf.h>


#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

#ifdef __cplusplus
extern "C" {
#endif


/* Network connection-related configurable info struct
 */
typedef struct {
  Nlm_Char*    client_host;      /* effective client hostname */
  Nlm_Char*    host;             /* server:  host */
  Nlm_Uint2    port;             /* server:  service port */
  Nlm_Char*    path;             /* server:  path(e.g. path to a CGI script) */
  Nlm_Char*    args;             /* server:  args(e.g. args for CGI script)  */
  STimeout     timeout;          /* i/o (connection and read/write) timeout  */
  Nlm_Uint4    conn_try;         /* max.number of attempts to establish conn */
  Nlm_Char*    http_proxy_host;  /* host name of HTTP proxy (can be NULL) */
  Nlm_Uint2    http_proxy_port;  /* port #    of HTTP proxy server */
  Nlm_Char*    proxy_host;       /* host of CERN-like firewall proxy server */
  Nlm_Boolean  debug_printout;   /* printout some debug info */
  Nlm_Boolean  firewall;         /* if to work in the firewall mode */
  Nlm_Uint2    ncbid_port;       /* port of the NCBID CGI script(used by LB) */
  Nlm_Char*    ncbid_path;       /* path to the NCBID CGI script(used by LB) */
  Nlm_Boolean  lb_disable;       /* if to disable the local load-balancing */

  /* the following field(s) are for the internal use only! */
  Nlm_Boolean  http_proxy_adjusted;
} SNetConnInfo;


/* Defaults and conf.parameter names for the SNetConnInfo fields
 */

#define DEF_CONN_CONF_FILE        "ncbi"
#define DEF_CONN_CONF_SECTION     "NET_SERV"
                                  
#define CFG_CONN_ENGINE_HOST      "SRV_ENGINE_HOST"
#define DEF_CONN_ENGINE_HOST      "www.ncbi.nlm.nih.gov"
                                  
#define CFG_CONN_ENGINE_PORT      "SRV_ENGINE_PORT"
#define DEF_CONN_ENGINE_PORT      80
                                  
#define CFG_CONN_ENGINE_PATH      "SRV_ENGINE_PATH"
#define DEF_CONN_ENGINE_PATH       "/Service/dispd.cgi"
                                  
#define CFG_CONN_ENGINE_ARGS      "SRV_ENGINE_ARGS"
#define DEF_CONN_ENGINE_ARGS      ""
                                  
#define CFG_CONN_TIMEOUT          "SRV_CONN_TIMEOUT"
#define DEF_CONN_TIMEOUT          30.0
                                  
#define CFG_CONN_TRY              "SRV_CONN_TRY"
#define DEF_CONN_TRY              3
                                  
#define CFG_CONN_PROXY_HOST       "SRV_PROXY_HOST"
#define DEF_CONN_PROXY_HOST       ""
                                  
#define CFG_CONN_HTTP_PROXY_HOST  "SRV_HTTP_PROXY_HOST"
#define DEF_CONN_HTTP_PROXY_HOST  ""
                                  
#define CFG_CONN_HTTP_PROXY_PORT  "SRV_HTTP_PROXY_PORT"
#define DEF_CONN_HTTP_PROXY_PORT  80
                                  
#define CFG_CONN_DEBUG_PRINTOUT   "SRV_DEBUG_PRINTOUT"
#define DEF_CONN_DEBUG_PRINTOUT   ""
                                  
#define CFG_CONN_FIREWALL         "SRV_FIREWALL"
#define DEF_CONN_FIREWALL         ""
                                  
#define CFG_CONN_NCBID_PORT       "SRV_NCBID_PORT"
#define DEF_CONN_NCBID_PORT       80
                                  
#define CFG_CONN_NCBID_PATH       "SRV_NCBID_PATH"
#define DEF_CONN_NCBID_PATH       "/Service/ncbid.cgi"
                                  
#define CFG_CONN_LB_DISABLE       "SRV_LB_DISABLE"
#define DEF_CONN_LB_DISABLE       ""


/* Limits for some SNetConnInfo fields
 */
#define MAX_CONN_HOST_LEN  64
#define MAX_CONN_PATH_LEN  1024
#define MAX_CONN_ARGS_LEN  1024


/* This function to fill out the "*info" structure;  it looks for the field
 * values in the application transient parameters(<conf_file>[<conf_section>]),
 * and/or environment variables, and/or configuration
 * file (<conf_file> [<conf_section>]).
 * For the configurable parameter lookup rules see in
 * "ncbienv.h":Nlm_GetEnvParam[Ex]().
 *
 * All fields will be resolved using the following env.var. / app.param.
 * keys.
 * For the field conf.parameter name see CNF_CONN_<KEY> (right above).
 * For the field default value see DEF_CONN_<KEY> (right above).
 *
 *  -- ARGUMENT ----------- KEY --------------- REMARKS/EXAMPLES ---------
 *   client_host          <will be assigned to the local host name>
 *   host                 ENGINE_HOST
 *   port                 ENGINE_PORT
 *   path                 ENGINE_PATH
 *   args                 ENGINE_ARGS
 *   timeout              CONN_TIMEOUT        "<sec>.<usec>": "30.0", "0.005"
 *   conn_try             CONN_TRY  
 *   http_proxy_host      HTTP_PROXY_HOST     no HTTP proxy if empty/NULL
 *   http_proxy_port      HTTP_PROXY_PORT
 *   proxy_host           PROXY_HOST
 *   debug_printout       DEBUG_PRINTOUT
 *   firewall             FIREWALL
 *   ncbid_port           NCBID_PORT
 *   ncbid_path           NCBID_PATH
 *   lb_disable           LB_DISABLE
 */
NLM_EXTERN SNetConnInfo* NetConnInfo_Create
(const Nlm_Char* conf_file,    /* if NULL/empty then DEF_CONN_CONF_FILE    */
 const Nlm_Char* conf_section  /* if NULL/empty then DEF_CONN_CONF_SECTION */
 );


/* Adjust the "host:port" to "proxy_host:proxy_port", and
 * "path" to "http://host:port/path" to connect through a HTTP proxy.
 * Return FALSE if already adjusted(see the NOTE).
 * NOTE:  it does nothing if applied more then once to the same "info"(or its
 *        clone), or when "http_proxy_host" is NULL.
 */
NLM_EXTERN Nlm_Boolean NetConnInfo_AdjustForHttpProxy
(SNetConnInfo* info
 );


/* Make an exact and independent copy of "*info".
 */
NLM_EXTERN SNetConnInfo* NetConnInfo_Clone
(const SNetConnInfo* info
 );


/* Printout the "*info" to file "fp".
 */
NLM_EXTERN void NetConnInfo_Print
(const SNetConnInfo* info,
 FILE*               fp
 );


/* Destroy "**info" and deallocate "*info" if "*info" is not NULL.
 * Assign "*info" to NULL.
 */
NLM_EXTERN void NetConnInfo_Destroy
(SNetConnInfo** info
 );


/* Convert SOCK-specific error status to the CONN's one
 */
NLM_EXTERN EConnStatus ESOCK2ECONN
(ESOCK_ErrCode err_code
 );


/* Hit URL "http://host:port/path?args":
 *    POST <path>?<args> HTTP/1.0\r\n
 *    <user_header\r\n>
 *    Content-Length: <content_length>\r\n\r\n
 * If "encode_args" is TRUE then URL-encode the "args".
 * "args" can be NULL/empty -- then the '?' symbol does not get added.
 * Use "c_timeout" to specify timeout for the "CONNECT" stage;
 * "rw_timeout" specifies timeout for the "READ"/"WRITE"
 * (see also "ncbisock.h:SOCK_SetTimeout()" for more info).
 * "content_length" is mandatory, and it specifies an exact(!) amount
 * of data to be sent to the resultant socket.
 * String "user_header"(can be NULL) will be written in the header end, and
 * it must be terminated by a single '\r\n' (if it is not NULL/empty).
 * 
 * On success, return non-NULL handle to a readable&writable
 * "NCBI socket"(see "ncbisock.h:SOCK_[Peek|Read|Write]").
 * ATTENTION:  due to the very essence of the HTTP connection, you can
 *             perform only one { WRITE, ..., WRITE, READ, ..., READ } cycle.
 * On error, return NULL.
 *
 * If "encode_args" is TRUE then URL-encode the "args", if any.
 *
 * NOTE: the socket must be closed by "ncbisock.h:SOCK_Close()" when not
 *       needed anymore.
 */

NLM_EXTERN SOCK Ncbi_ConnectURL
(const Nlm_Char* host,
 Nlm_Uint2       port,
 const Nlm_Char* path,
 const Nlm_Char* args,
 size_t          content_length,
 const STimeout* c_timeout,
 const STimeout* rw_timeout,
 const Nlm_Char* user_header,
 Nlm_Boolean     encode_args
 );


/* Discard all input data before(and including) the first occurence of
 * "pattern". If "buf" is not NULL then add the discarded data(including
 * the "pattern") to it. If "n_discarded" is not NULL then "*n_discarded"
 * will return # of discarded bytes.
 */
NLM_EXTERN EConnStatus CONN_StripToPattern
(CONN        conn,
 const void* pattern,
 Nlm_Uint4   pattern_size,
 BUF*        buf,
 Nlm_Uint4*  n_discarded
 );

NLM_EXTERN EConnStatus SOCK_StripToPattern
(SOCK        sock,
 const void* pattern,
 Nlm_Uint4   pattern_size,
 BUF*        buf,
 Nlm_Uint4*  n_discarded
 );


/* URL-decode up to "src_size" symbols(bytes) from buffer "src_buf".
 * Write the decoded data to buffer "dst_buf", but no more than "dst_size"
 * bytes.
 * Assign "*src_read" to the # of bytes succesfully decoded from "src_buf".
 * Assign "*dst_written" to the # of bytes written to buffer "dst_buf".
 * Return FALSE only if cannot decode nothing, and an unrecoverable
 * URL-encoding error (such as an invalid symbol or a bad "%.." sequence)
 * has occured.
 * NOTE:  the unfinished "%.." sequence is fine -- return TRUE, but dont
 *        "read" it.
 */
NLM_EXTERN Nlm_Boolean URL_Decode
(const void* src_buf,    /* [in]     non-NULL */
 Nlm_Uint4   src_size,   /* [in]              */
 Nlm_Uint4*  src_read,   /* [out]    non-NULL */
 void*       dst_buf,    /* [in/out] non-NULL */
 Nlm_Uint4   dst_size,   /* [in]              */
 Nlm_Uint4*  dst_written /* [out]    non-NULL */
 );


/* Act just like URL_Decode (see above) but caller can allow the specified
 * non-standard URL symbols in the input buffer to be decoded "as is".
 * The extra allowed symbols are passed in a '\0'-terminated string
 * "allow_symbols" (it can be NULL or empty -- then it will be an exact
 * equivalent of URL_Decode).
 */
NLM_EXTERN Nlm_Boolean URL_DecodeEx
(const void* src_buf,      /* [in]     non-NULL  */
 Nlm_Uint4   src_size,     /* [in]               */
 Nlm_Uint4*  src_read,     /* [out]    non-NULL  */
 void*       dst_buf,      /* [in/out] non-NULL  */
 Nlm_Uint4   dst_size,     /* [in]               */
 Nlm_Uint4*  dst_written,  /* [out]    non-NULL  */
 Char*       allow_symbols /* [in]     '\0'-term */
 );


/* URL-encode up to "src_size" symbols(bytes) from buffer "src_buf".
 * Write the encoded data to buffer "dst_buf", but no more than "dst_size"
 * bytes.
 * Assign "*src_read" to the # of bytes succesfully encoded from "src_buf".
 * Assign "*dst_written" to the # of bytes written to buffer "dst_buf".
 */
NLM_EXTERN void URL_Encode
(const void* src_buf,    /* [in]     non-NULL */
 Nlm_Uint4   src_size,   /* [in]              */
 Nlm_Uint4*  src_read,   /* [out]    non-NULL */
 void*       dst_buf,    /* [in/out] non-NULL */
 Nlm_Uint4   dst_size,   /* [in]              */
 Nlm_Uint4*  dst_written /* [out]    non-NULL */
 );


/****************************************************************************
 * NCBI-specific MIME content type and sub-types
 * (the API to compose and parse them)
 *    Content-Type: <type>/<MIME_ComposeSubType()>\r\n
 *
 *    Content-Type: <type>/<subtype>-<encoding>\r\n
 *
 * where  MIME_ComposeSubType(EMIME_SubType subtype, EMIME_Encoding encoding):
 *   "x-<subtype>-<encoding>":
 *     "x-<subtype>",   "x-<subtype>-urlencoded",   "x-<subtype>-<encoding>",
 *     "x-asn-text",    "x-asn-text-urlencoded",    "x-asn-text-<encoding>
 *     "x-asn-binary",  "x-asn-binary-urlencoded",  "x-asn-binary-<encoding>"
 *     "x-www-form",    "x-www-form-urlencoded",    "x-www-form-<encoding>"
 *     "html",          "html-urlencoded",          "html-<encoding>"
 *     "x-unknown",     "x-unknown-urlencoded",     "x-unknown-<encoding>"
 *
 *  Note:  <subtype> and <encoding> are expected to contain only
 *         alphanumeric symbols, '-' and '_'. They are case-insensitive.
 ****************************************************************************/

/* Type
 */
typedef enum {
  eMIME_T_NcbiData = 0,  /* "x-ncbi-data"  (NCBI specific data) */
  eMIME_T_Text,          /* "text" */
  eMIME_T_Application,   /* "application" */
  /* eMIME_T_???, "<type>"  here go other types */
  eMIME_T_Unknown        /* "unknown" */
} EMIME_Type;

/* SubType
 */
typedef enum {
  eMIME_AsnText = 0,  /* "x-asn-text"    (text ASN.1 data) */
  eMIME_AsnBinary,    /* "x-asn-binary"  (binary ASN.1 data) */
  eMIME_Fasta,        /* "x-fasta"       (data in FASTA format) */
  eMIME_WwwForm,      /* "x-www-form" */
  /* standard MIMEs */
  eMIME_Html,         /* "html" */
  /* eMIME_???,          "<subtype>"   here go other NCBI subtypes */
  eMIME_Unknown       /* "x-unknown"     (an arbitrary binary data) */
} EMIME_SubType;

/* Encoding
 */
typedef enum {
  eENCOD_Url = 0,  /* "-urlencoded"  (the content is URL-encoded) */
  /* eENCOD_???,      "-<encoding>"   here go other NCBI encodings */
  eENCOD_None      /* ""              (the content is passed "as is") */
} EMIME_Encoding;


/* Write up to "buflen" bytes to "buf":
 *   Content-Type: <type>/[x-]<subtype>-<encoding>\r\n
 * Return pointer to the "buf".
 */
#define MAX_CONTENT_TYPE_LEN 64
NLM_EXTERN char* MIME_ComposeContentTypeEx
(EMIME_Type     type,
 EMIME_SubType  subtype,
 EMIME_Encoding encoding,
 char*          buf,
 size_t         buflen    /* must be at least MAX_CONTENT_TYPE_LEN */
 );

/* Exactly equivalent to MIME_ComposeContentTypeEx(eMIME_T_NcbiData, ...)
 */
NLM_EXTERN char* MIME_ComposeContentType
(EMIME_SubType  subtype,
 EMIME_Encoding encoding,
 char*          buf,
 size_t         buflen
 );

/* Parse the NCBI-specific content-type; the (case-insensitive) "str"
 * can be in the following two formats:
 *   Content-Type: <type>/x-<subtype>-<encoding>
 *   <type>/x-<subtype>-<encoding>
 * If it does not match any of NCBI MIME type/subtypes/encodings, then
 * return TRUE, eMIME_T_Unknown, eMIME_Unknown or eENCOD_None, respectively.
 * If the passed "str" has an invalid (non-HTTP ContentType) format
 * (or if it is NULL/empty), then
 * return FALSE, eMIME_T_Unknown, eMIME_Unknown, and eENCOD_None.
 */
NLM_EXTERN Nlm_Boolean MIME_ParseContentTypeEx
(const char*     str,      /* the HTTP "Content-Type:" header to parse */
 EMIME_Type*     type,     /* can be NULL */
 EMIME_SubType*  subtype,  /* can be NULL */
 EMIME_Encoding* encoding  /* can be NULL */
 );

/* Requires the MIME type be "x-ncbi-data"
 */
NLM_EXTERN Nlm_Boolean MIME_ParseContentType
(const char*     str,      /* the HTTP "Content-Type:" header to parse */
 EMIME_SubType*  subtype,  /* can be NULL */
 EMIME_Encoding* encoding  /* can be NULL */
 );


#ifdef __cplusplus
} /* extern "C" */
#endif

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif /* CONNUTIL__H */
