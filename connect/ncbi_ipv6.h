#ifndef CONNECT___NCBI_IPV6__H
#define CONNECT___NCBI_IPV6__H

/* $Id: ncbi_ipv6.h,v 1.1 2016/12/31 22:35:10 ucko Exp $
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
 *   IPv6 addressing support
 *
 */

#include <connect/connect_export.h>
#include <stddef.h>


#ifdef __cplusplus
extern "C" {
#endif /*__cplusplus*/


typedef struct {
    unsigned char octet[16];  /* no alignment assumed */
} TNCBI_IPv6Addr;


/** Return non-zero if the address is empty (either as IPv6 or IPv4); return
 *  zero otherwise.
 */
extern NCBI_XCONNECT_EXPORT
int/*bool*/ NcbiIsEmptyIPv6(const TNCBI_IPv6Addr* addr);


/** Return non-zero if the address is either IPv4 compatible or a mapped IPv4
 *  address; return zero otherwise.
 */
extern NCBI_XCONNECT_EXPORT
int/*bool*/  NcbiIsIPv4    (const TNCBI_IPv6Addr* addr);


extern NCBI_XCONNECT_EXPORT
unsigned int NcbiIPv6ToIPv4(const TNCBI_IPv6Addr* addr, size_t pfxlen);


extern NCBI_XCONNECT_EXPORT
int/*bool*/  NcbiIPv4ToIPv6(TNCBI_IPv6Addr* addr,
                            unsigned int ipv4, size_t pfxlen);

/** Convert into an IPv4 address, the first "len" (or "strlen(str)" if "len" is
 *  0) bytes of "str" from a full-quad decimal notation; return a non-zero
 *  string pointer to the first non-converted character (which is neither a
 *  digit nor a dot); return 0 if conversion failed and no IPv4 address had
 *  been found.
 */
extern NCBI_XCONNECT_EXPORT
const char*  NcbiStringToIPv4(unsigned int* addr,
                              const char* str, size_t len);


/** Convert into an IPv6 address, the first "len" (or "strlen(str)" if "len" is
 *  0) bytes of "str" from a hexadecimal colon-separated notation (including
 *  full-quad trailing IPv4); return a non-zero string pointer to the first
 *  non-converted character (which is neither a hex-digit, nor a colon, nor a
 *  dot); return 0 if conversion failed and no IPv6 address had been found.
 */
extern NCBI_XCONNECT_EXPORT
const char*  NcbiStringToIPv6(TNCBI_IPv6Addr* addr,
                             const char* str, size_t len);


/** Convert into an IPv6 address, the first "len" (or "strlen(str)" if "len" is
 *  0) bytes of "str" from either a full-quad decimal IPv4 or a hexadecimal
 *  colon-separated IPv6; return a non-zero string pointer to the first
 *  non-converted character (which is neither a [hex-]digit, nor a colon, nor a
 *  dot); return 0 if  no conversion can be made.
 */
extern NCBI_XCONNECT_EXPORT
const char*  NcbiIPToAddr(TNCBI_IPv6Addr* addr,
                          const char* str, size_t len);


/** Convert into an IPv6 address, the first "len" (or "strlen(str)" if "len" is
 *  0) bytes of "str", which can be either of a full-quad decimal IPv4, a
 *  hexadecimal colon-separated IPv6, an .in-addr.arpa- or an .in6.arpa-domain
 *  names; return a non-zero string pointer to the first non-converted
 *  character (which is neither a [hex-]digit, nor a colon, nor a dot); return
 *  0 if no conversion can be made.
 */
extern NCBI_XCONNECT_EXPORT
const char*  NcbiStringToAddr(TNCBI_IPv6Addr* addr,
                              const char* str, size_t len);


/** Convert network byte order IPv4 into a full-quad text form and store the
 *  result in the "buf" of size "bufsize".  Return non-zero string address
 *  past the stored result, or 0 when the conversion failed for buffer being
 *  too small.
 */
extern NCBI_XCONNECT_EXPORT
char*        NcbiIPv4ToString(char* buf, size_t bufsize,
                              unsigned int addr);


/** Convert IPv6 address into a hex colon-separated text form and store the
 *  result in the "buf" of size "bufsize".  Return non-zero string address
 *  past the stored result, or 0 when the conversion failed for buffer being
 *  too small.
 */
extern NCBI_XCONNECT_EXPORT
char*        NcbiIPv6ToString(char* buf, size_t bufsize,
                              const TNCBI_IPv6Addr* addr);


/** Convert IPv6 address into either a full-quad text IPv4 (for IPv4-compatible
 *  IPv6 addresses) or a hex colon-separated text form(for all other) and store
 *  the result in the "buf" of size "bufsize".  Return non-zero string address
 *  past the stored result, or 0 when the conversion failed for buffer being
 *  too small.
 */
extern NCBI_XCONNECT_EXPORT
char*        NcbiAddrToString(char* buf, size_t bufsize,
                              const TNCBI_IPv6Addr* addr);


/** Convert IPv6 address into either .in-addr.arpa domain (for IPv4-compatible
 *  IPv6 addresses) or .ip6.arpa domain (for all other) and store the result in
 *  the "buf" of size "bufsize".  Return non-zero string address past the
 *  stored result, or 0 when the conversion failed for buffer being too small.
 */
extern NCBI_XCONNECT_EXPORT
const char*  NcbiAddrToDNS(char* buf, size_t bufsize,
                           const TNCBI_IPv6Addr* addr);


/** Return non-zero if "addr" belongs to the network specified as CIDR
 *  "base/bits"; return zero otherwise.
 */
extern NCBI_XCONNECT_EXPORT
int/*bool*/  NcbiIsInIPv6Network(const TNCBI_IPv6Addr* base,
                                 unsigned int          bits,
                                 const TNCBI_IPv6Addr* addr);


#ifdef __cplusplus
}
#endif /*__cplusplus*/


#endif  /* CONNECT___NCBI_IPV6__H */
