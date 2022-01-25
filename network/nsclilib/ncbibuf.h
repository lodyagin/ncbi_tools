#ifndef NCBIBUF__H
#define NCBIBUF__H

/*  $RCSfile: ncbibuf.h,v $  $Revision: 6.1 $  $Date: 1998/04/14 15:34:23 $
* ==========================================================================
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
* ==========================================================================
*
* Author:  Denis Vakatov
*
* File Description:
*   Memory-resided FIFO storage area(to be used e.g. in I/O buffering)
*
* --------------------------------------------------------------------------
* $Log: ncbibuf.h,v $
* Revision 6.1  1998/04/14 15:34:23  vakatov
* Initial revision
*
* ==========================================================================
*/

#include <ncbistd.h>


#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

#ifdef __cplusplus
extern "C" {
#endif


#define BUF         Nlm_BUF
#define BUF_Size    Nlm_BUF_Size
#define BUF_Write   Nlm_BUF_Write
#define BUF_Peek    Nlm_BUF_Peek
#define BUF_Read    Nlm_BUF_Read
#define BUF_Destroy Nlm_BUF_Destroy


struct Nlm_BufferTag;
typedef struct Nlm_BufferTag* BUF;  /* handle of the NCBI buffer */


/* Return the number of bytes stored in "buf".
 * NOTE: return 0 if "buf" == NULL
 */
NLM_EXTERN Nlm_Uint4 BUF_Size(BUF buf);


/* Add new data to "*pBuf".
 * NOTE:  if "*pBuf" == NULL then create it
 */
NLM_EXTERN Nlm_Boolean BUF_Write(BUF* pBuf, const void* data, Nlm_Uint4 size);


/* Copy up to "size" bytes stored in "buf" to "data".
 * Return the # of copied bytes(can be less than "size").
 * NOTE:  "buf" and "data" can be NULL; in both cases, do nothing
 *        and return 0.
 */
NLM_EXTERN Nlm_Uint4 BUF_Peek(BUF buf, void* data, Nlm_Uint4 size);


/* Copy up to "size" bytes stored in "buf" to "data" and remove
 * copied data from the "buf".
 * Return the # of copied-and/or-removed bytes(can be less than "size")
 * NOTE: if "buf"  == NULL then do nothing and return 0
 *       if "data" == NULL then do not copy data anywhere(still, remove it)
 */
NLM_EXTERN Nlm_Uint4 BUF_Read(BUF buf, void* data, Nlm_Uint4 size);


/* Destroy all internal data;  return NULL
 * NOTE: do nothing if "buf" == NULL
 */
NLM_EXTERN BUF BUF_Destroy(BUF buf);


#ifdef __cplusplus
}
#endif

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif /* NCBIBUF__H */
