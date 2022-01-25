/*  ncbibs.h
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
* File Name:  ncbibs.h
*
* Author:  Jim Ostell
*
* Version Creation Date:  1/1/91
*
* $Revision: 6.0 $
*
* File Description:
*   ByteStore typedefs, prototypes, and defines
*
* Modifications:
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
* 04-15-93 Schuler     Changed _cdecl to LIBCALL
*
*
* $Log: ncbibs.h,v $
* Revision 6.0  1997/08/25 18:15:16  madden
* Revision changed to 6.0
*
* Revision 5.1  1996/12/03 21:48:33  vakatov
* Adopted for 32-bit MS-Windows DLLs
*
 * Revision 5.0  1996/05/28  13:18:57  ostell
 * Set to revision 5.0
 *
 * Revision 4.1  1995/12/21  02:33:19  ostell
 * added use_min_size to BSAll
 *
 * Revision 4.0  1995/07/26  13:46:50  ostell
 * force revision to 4.0
 *
 * Revision 2.4  1995/05/15  18:45:58  ostell
 * added Log line
 *
*
*
* ==========================================================================
*/

#ifndef _NCBIBS_
#define _NCBIBS_

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif


#ifdef __cplusplus
extern "C" {
#endif

typedef struct bsunit {             /* for building multiline strings */
	Nlm_Handle str;            /* the string piece */
	Nlm_Int2 len_avail,
		 len;
	struct bsunit PNTR next; }       /* the next one */
Nlm_BSUnit, PNTR Nlm_BSUnitPtr;

typedef struct bytestore {
	Nlm_Int4 seekptr,       /* current position */
		totlen,             /* total stored data length in bytes */
		chain_offset;       /* offset in ByteStore of first byte in curchain */
	Nlm_BSUnitPtr chain,       /* chain of elements */
		curchain;           /* the BSUnit containing seekptr */
} Nlm_ByteStore, PNTR Nlm_ByteStorePtr;

NLM_EXTERN Nlm_VoidPtr LIBCALL Nlm_BSMerge PROTO((Nlm_ByteStorePtr ssp, Nlm_VoidPtr dest));
NLM_EXTERN Nlm_ByteStorePtr LIBCALL Nlm_BSNew PROTO((Nlm_Int4 len));
NLM_EXTERN Nlm_Int2 LIBCALL Nlm_BSSeek PROTO((Nlm_ByteStorePtr bsp, Nlm_Int4 offset, Nlm_Int2 origin));
NLM_EXTERN Nlm_Int4 LIBCALL Nlm_BSTell PROTO((Nlm_ByteStorePtr bsp));
NLM_EXTERN Nlm_Int4 LIBCALL Nlm_BSDelete PROTO((Nlm_ByteStorePtr bsp, Nlm_Int4 len));
NLM_EXTERN Nlm_Int4 LIBCALL Nlm_BSWrite PROTO((Nlm_ByteStorePtr bsp, Nlm_VoidPtr ptr, Nlm_Int4 len));
NLM_EXTERN Nlm_Int4 LIBCALL Nlm_BSInsert PROTO((Nlm_ByteStorePtr bsp, Nlm_VoidPtr ptr, Nlm_Int4 len));
NLM_EXTERN Nlm_Int4 LIBCALL Nlm_BSInsertFromBS PROTO((Nlm_ByteStorePtr bsp, Nlm_ByteStorePtr bsp2, Nlm_Int4 len));
NLM_EXTERN Nlm_Int4 LIBCALL Nlm_BSRead PROTO((Nlm_ByteStorePtr bsp, Nlm_VoidPtr ptr, Nlm_Int4 len));
NLM_EXTERN Nlm_Int2 LIBCALL Nlm_BSGetByte PROTO((Nlm_ByteStorePtr bsp));
NLM_EXTERN Nlm_Int2 LIBCALL Nlm_BSPutByte PROTO((Nlm_ByteStorePtr bsp, Nlm_Int2 value));
NLM_EXTERN Nlm_Int4 LIBCALL Nlm_BSLen PROTO((Nlm_ByteStorePtr ssp));
NLM_EXTERN Nlm_ByteStorePtr LIBCALL Nlm_BSFree PROTO((Nlm_ByteStorePtr ssp));
NLM_EXTERN Nlm_ByteStorePtr LIBCALL Nlm_BSDup PROTO((Nlm_ByteStorePtr source));
/*****************************************************************************
*
*   Int4 BSAdd(bsp, len, use_min_size)
*   	adds len bytes BEFORE current bsp->seekptr
*       bsp->seekptr returned pointing at first added byte
*   	returns bytes added
*       if (use_min_size) then does not add anything smaller than MIN_BSALLOC
*
*****************************************************************************/
NLM_EXTERN Nlm_Int4 LIBCALL Nlm_BSAdd PROTO((Nlm_ByteStorePtr bsp, Nlm_Int4 len, Nlm_Boolean use_min_size));


#define ByteStore Nlm_ByteStore
#define ByteStorePtr Nlm_ByteStorePtr
#define BSMerge Nlm_BSMerge
#define BSNew Nlm_BSNew
#define BSLen Nlm_BSLen
#define BSFree Nlm_BSFree
#define BSSeek Nlm_BSSeek
#define BSTell Nlm_BSTell
#define BSDelete Nlm_BSDelete
#define BSWrite Nlm_BSWrite
#define BSInsert Nlm_BSInsert
#define BSInsertFromBS Nlm_BSInsertFromBS
#define BSDup Nlm_BSDup
#define BSRead Nlm_BSRead
#define BSGetByte Nlm_BSGetByte
#define BSPutByte Nlm_BSPutByte

#ifdef __cplusplus
}
#endif


#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif

