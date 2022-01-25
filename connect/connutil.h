#ifndef CONNUTIL__H
#define CONNUTIL__H

/*  $Id: connutil.h,v 6.3 1999/04/09 22:27:01 vakatov Exp $
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
*
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
  Nlm_STimeout timeout;          /* i/o (connection and read/write) timeout  */
  Nlm_Uint4    conn_try;         /* max.number of attempts to establish conn */
  Nlm_Char*    proxy_host;
  Nlm_Uint2    proxy_port;
  Nlm_Boolean  proxy_transparent;
  Nlm_Boolean  debug_printout;
  Nlm_Boolean  firewall;
  /* the following field(s) are for the internal use only! */
  Nlm_Boolean  proxy_adjusted;
} SNetConnInfo;


/* This function to fill out the "*info" structure; it looks for the field
 * values in the application transient parameters(<conf_file>[<conf_section>]),
 * and/or environment variables, and/or configuration
 * file (<conf_file> [<conf_section>]);  for the configurable parameter lookup
 * rules see in "ncbienv.h":Nlm_GetEnvParam[Ex]().
 *
 * All fields will be resolved using the following env.var. / app.param.
 * keys and defaults:
 *
 *  -- ARGUMENT ------------ DEFAULT ----------------------- KEY -----
 *   client_host          <will be assigned to the local host name>
 *   host                 "www.ncbi.nlm.nih.gov"        SRV_ENGINE_HOST
 *   port                 80                            SRV_ENGINE_PORT
 *   path                 "Service/nph-dispd.cgi"       SRV_ENGINE_PATH
 *   args                 ""                            SRV_ENGINE_ARGS
 *   timeout              30.0                          SRV_CONN_TIMEOUT
 *   conn_try             3                             SRV_CONN_TRY
 *   proxy_host           NULL(means no proxy stuff)    SRV_PROXY_HOST
 *   proxy_port           80                            SRV_PROXY_PORT
 *   proxy_transparent    TRUE                          SRV_PROXY_TRANSPARENT
 *   debug_printout       FALSE                         SRV_DEBUG_PRINTOUT
 *   firewall             FALSE                         SRV_FIREWALL
 */
NLM_EXTERN SNetConnInfo* NetConnInfo_Create
(const Nlm_Char* conf_file,    /* by default(if NULL/empty) -- "ncbi" */
 const Nlm_Char* conf_section  /* by default(if NULL/empty) -- "NET_SERV" */
 );

/* Adjust the "host:port" to "proxy_host:proxy_port", and
 * "path" to "http://host:port/path" to connect through a CERN-style proxy.
 * Return FALSE if already adjusted(see the NOTE).
 * NOTE:  it does nothing if applied more then once to the same "info"(or its
 *        clone), or when "proxy_host" is NULL.
 */
NLM_EXTERN Nlm_Boolean NetConnInfo_AdjustForCernProxy
(SNetConnInfo* info
 );

/* Make an exact and independent copy of "*info".
 */
NLM_EXTERN SNetConnInfo* NetConnInfo_Clone
(const SNetConnInfo* info
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


/* Hit URL "http://host:port/path?args".
 * If "encode_args" is TRUE then URL-encode the "args".
 * "args" can be NULL/empty -- then the '?' symbol does not get added.
 * Use "c_timeout" to specify timeout for the "CONNECT" stage;
 * "rw_timeout" specifies timeout for the "READ"/"WRITE"
 * (see also "ncbisock.h:SOCK_SetTimeout()" for more info).
 * "content_length" is mandatory, and it specifies an exact(!) amount
 * of data to be sent to the resultant socket.
 * "user_header"(can be NULL) will be written in the header end, and
 * it must be terminated by a single '\n' (if it is not NULL/empty).
 * 
 * On success, return non-NULL handle to a readable&writable
 * "NCBI socket"(see "ncbisock.h:SOCK_[Peek|Read|Write]").
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
(SOCK        conn,
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
 * Return FALSE on unrecoverable URL-encoding error, such as an invalid symbol
 * or a bad "%.." sequence.
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


#ifdef __cplusplus
}
#endif

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif /* CONNUTIL__H */
