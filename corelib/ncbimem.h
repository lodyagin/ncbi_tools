/*   ncbimem.h
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
* File Name:  ncbimem.h
*
* Author:  Gish, Kans, Ostell, Schuler
*
* Version Creation Date:   1/1/91
*
* $Revision: 6.4 $
*
* File Description:
*   	prototypes for ncbi memory functions
*
* Modifications:
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
* 9-19-91  Schuler     Modified existing prototypes for ANSI-resemblance
* 9-19-91  Schuler     Added new prototypes for Windows ANSI-like functions
* 9-19-91  Schuler     Changed all functions to _cdecl calling convention
* 04-15-93 Schuler     Changed _cdecl to LIBCALL
* 05-21-93 Schuler     Nlm_MemFreeTrace added for debugging MemFree
* 06-14-93 Schuler     Added dll_Malloc and dll_Free
*
* $Log: ncbimem.h,v $
* Revision 6.4  1998/06/23 16:30:03  vakatov
* MemMapInit(): return zero-filled structure(not a NULL) if file is empty
*
* Revision 6.3  1998/03/26 17:54:12  vakatov
* Added Nlm_CallocViaMalloc to replace native "calloc" on Solaris
*
* Revision 6.2  1998/03/20 17:09:39  vakatov
* [OS_UNIX] Added Nlm_SetHeapLimit() //NOTE: has no effect on OSF1(?)
*
* Revision 6.1  1997/09/09 23:42:22  vakatov
* struct Nlm_MemMap::  made "hMap" field be WIN32-only;  removed "hFile" field
* Nlm_MemMapInit():: use '/' instead of '\' in the file-mapping object
* name to allow specify(and then -- reuse) "pathed" filenames
*
* Revision 6.0  1997/08/25 18:16:43  madden
* Revision changed to 6.0
*
* Revision 5.5  1997/07/22 19:11:36  vakatov
* Separated Main() from GetArg[svc]() functions;  [WIN_MSWIN] converged
* console and GUI libraries; [for WIN32-DLL] encapsulated global variables
*
* Revision 5.4  1997/01/06 22:26:54  vakatov
* [WIN16,WIN32][_DEBUG]  MemFree --> Nlm_MemFreeTrace (included [_CONSOLE])
*
 * Revision 5.3  1996/12/13  19:12:42  epstein
 * add MMAP_AVAIL definition
 *
 * Revision 5.2  1996/12/03  21:48:33  vakatov
 * Adopted for 32-bit MS-Windows DLLs
 *
 * Revision 5.1  1996/08/29  20:50:44  madden
 * Added functions for memory-mapping.
 *
 * Revision 5.0  1996/05/28  13:18:57  ostell
 * Set to revision 5.0
 *
 * Revision 4.0  1995/07/26  13:46:50  ostell
 * force revision to 4.0
 *
 * Revision 2.15  1995/05/15  18:45:58  ostell
 * added Log line
 *
*
*
* ==========================================================================
*/
#ifndef _NCBIMEM_
#define _NCBIMEM_

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif


#ifdef __cplusplus
extern "C" {
#endif

/* ======== PROTOTYPES ======== */

NLM_EXTERN void * LIBCALL Nlm_MemNew PROTO((size_t size));
NLM_EXTERN void * LIBCALL Nlm_MemGet PROTO((size_t size, unsigned int flags));
NLM_EXTERN void * LIBCALL Nlm_MemMore PROTO((void *ptr, size_t size));
NLM_EXTERN void * LIBCALL Nlm_MemExtend PROTO((void *ptr, size_t size, size_t oldsize));
NLM_EXTERN void * LIBCALL Nlm_MemFree PROTO((void *ptr));
NLM_EXTERN void * LIBCALL Nlm_MemCopy PROTO((void *dst, const void *src, size_t bytes));
NLM_EXTERN void * LIBCALL Nlm_MemMove PROTO((void *dst, const void *src, size_t bytes));
NLM_EXTERN void * LIBCALL Nlm_MemFill PROTO((void *ptr, int value, size_t bytes));
NLM_EXTERN void * LIBCALL Nlm_MemDup PROTO((const void *orig, size_t size));

#if defined(OS_MAC) || defined(OS_MSWIN) || defined(MSC_VIRT)
NLM_EXTERN Nlm_Handle  LIBCALL Nlm_HandNew PROTO((size_t size));
NLM_EXTERN Nlm_Handle  LIBCALL Nlm_HandGet PROTO((size_t size, Nlm_Boolean clear_out));
NLM_EXTERN Nlm_Handle  LIBCALL Nlm_HandMore PROTO((Nlm_Handle hnd, size_t size));
NLM_EXTERN Nlm_Handle  LIBCALL Nlm_HandFree PROTO((Nlm_Handle hnd));
NLM_EXTERN void *LIBCALL Nlm_HandLock PROTO((Nlm_Handle hnd));
NLM_EXTERN void *LIBCALL Nlm_HandUnlock PROTO((Nlm_Handle hnd));
#endif

#ifdef WIN16
NLM_EXTERN void *  LIBCALL win16_Malloc (size_t bytes);
NLM_EXTERN void *  LIBCALL win16_Calloc (size_t items, size_t size);
NLM_EXTERN void *  LIBCALL win16_Realloc (void *ptr, size_t size);
NLM_EXTERN void    LIBCALL win16_Free (void *ptr);
#endif

#ifdef WIN_MAC
#ifdef USE_MAC_MEMORY
NLM_EXTERN void *mac_Malloc  PROTO((size_t size));
NLM_EXTERN void *mac_Calloc  PROTO((size_t nmemb, size_t size));
NLM_EXTERN void *mac_Realloc PROTO((void *ptr, size_t size));
NLM_EXTERN void  mac_Free    PROTO((void *ptr));
#endif
#endif


/* [UNIX only] set the limit for the heap size and its increase policy for
 *  NCBI memory allocation functions:
 *    MemGet, MemNew, MemMore, MemExtend 
 *  when the heap size reaches "curr", it ussues a warning, then it increases
 *  "curr" by "add" bytes -- unless "curr" has already reached "max".
 *  in the latter case, program ussues a FATAL_ERROR error message and
 *  the NCBI allocation function returns NULL
 *  NOTE: specifying "curr" == 0 switch off the heap restriction algorithm;
 *        and it is off by default(if Nlm_SetHeapLimit has not been called)
 */
NLM_EXTERN Nlm_Boolean Nlm_SetHeapLimit(size_t curr, size_t add, size_t max);


/* Do Nlm_Calloc by {Nlm_Malloc + Nlm_MemSet} rather than native {calloc}
 */
NLM_EXTERN void* Nlm_CallocViaMalloc(size_t n_elem, size_t item_size);


/* ========= MACROS ======== */

/* low-level ANSI-style functions */
#ifdef WIN16
#define Nlm_Malloc  win16_Malloc
#define Nlm_Calloc  win16_Calloc
#define Nlm_Realloc win16_Realloc
#define Nlm_Free    win16_Free
#define Nlm_MemSet  _fmemset
#define Nlm_MemCpy  _fmemcpy
#define Nlm_MemChr  _fmemchr
#define Nlm_MemCmp  _fmemcmp
#else
#ifdef USE_MAC_MEMORY
#define Nlm_Malloc  mac_Malloc
#define Nlm_Calloc  mac_Calloc
#define Nlm_Realloc mac_Realloc
#define Nlm_Free    mac_Free
#define Nlm_MemSet  memset
#define Nlm_MemCpy  memcpy
#define Nlm_MemChr  memchr
#define Nlm_MemCmp  memcmp
#else
#define Nlm_Malloc  malloc
#define Nlm_Calloc  calloc
#define Nlm_Realloc realloc
#define Nlm_Free    free
#define Nlm_MemSet  memset
#define Nlm_MemCpy  memcpy
#define Nlm_MemChr  memchr
#define Nlm_MemCmp  memcmp
#endif
#endif


#ifdef OS_UNIX_SOL
/* Kludge for Solaris("calloc()" sometimes fails in MT aplications) */
#undef Nlm_Calloc
#define Nlm_Calloc Nlm_CallocViaMalloc
#endif


#define Malloc  Nlm_Malloc
#define Calloc  Nlm_Calloc
#define Realloc Nlm_Realloc
#define Free    Nlm_Free
#define MemSet  Nlm_MemSet
#define MemCpy  Nlm_MemCpy
#define MemChr  Nlm_MemChr
#define MemCmp  Nlm_MemCmp

/*** High-level NCBI functions ***/

/* Fake handle functions with pointer functions */

#if !(defined(OS_MAC) || defined(OS_MSWIN) || defined(MSC_VIRT))
#define Nlm_HandNew(a)    Nlm_MemNew(a)
#define Nlm_HandGet(a,b)  Nlm_MemGet(a,b)
#define Nlm_HandMore(a,b) Nlm_MemMore(a,b)
#define Nlm_HandFree(a)   Nlm_MemFree(a)
#define Nlm_HandLock(a)   (a)
#define Nlm_HandUnlock(a) NULL
#endif

/* Pointer functions */
#define MemNew(x)     Nlm_MemGet(x,MGET_CLEAR|MGET_ERRPOST)
#define MemGet(x,y)   Nlm_MemGet(x,(unsigned int)(y))
#define MemFree       Nlm_MemFree
#define MemMore       Nlm_MemMore
#define MemExtend     Nlm_MemExtend
#define MemCopy       Nlm_MemCopy
#define MemMove       Nlm_MemMove
#define MemFill       Nlm_MemFill
#define MemDup        Nlm_MemDup

#define HandNew     Nlm_HandNew
#define HandGet     Nlm_HandGet
#define HandMore    Nlm_HandMore
#define HandFree    Nlm_HandFree
#define HandLock    Nlm_HandLock
#define HandUnlock  Nlm_HandUnlock

#if (defined(OS_UNIX_SYSV) || defined(OS_UNIX_SUN) || defined(OS_UNIX_OSF1) || defined(OS_UNIX_LINUX) || defined(OS_UNIX_AIX)) && !defined(OS_UNIX_ULTRIX)
#define MMAP_AVAIL
#endif



#if defined(_DEBUG)  &&  defined(OS_MSWIN)
NLM_EXTERN void * LIBCALL Nlm_MemFreeTrace (void *, const char*, const char*, int);
#undef MemFree
#define MemFree(_ptr_)  Nlm_MemFreeTrace(_ptr_,THIS_MODULE,THIS_FILE,__LINE__)
#endif


#ifdef _WINDLL
NLM_EXTERN void * dll_Malloc (size_t bytes);
NLM_EXTERN void   dll_Free (void *pMem);
#else
#define dll_Malloc(x)	(void*)Nlm_Malloc(x)
#define dll_Free(x)	Nlm_Free((void*)(x))
#endif



/* flags for MemGet */
#define MGET_CLEAR    0x0001
#define MGET_ERRPOST  0x0004


/* obsolete */
#define MG_CLEAR    MGET_CLEAR
#define MG_MAXALLOC 0x0002
#define MG_ERRPOST  MGET_ERRPOST

/* structures for memory-mapping. */

/* This structure is allocated and filled by Nlm_MemMapInit.
The Nlm_Handle's are used by WIN32, "file_size" is used by
UNIX memory mapping when the the files are unmapped. */

typedef struct _nlm_mem_map
{
#ifdef WIN32
  Nlm_Handle hMap;
#endif
  Nlm_Int4    file_size;
  Nlm_CharPtr mmp_begin;
} Nlm_MemMap, PNTR Nlm_MemMapPtr;

/* prototypes for memory-mapping. */

/* Used to determine if memory-mapping is supported by the NCBI toolkit. */
NLM_EXTERN Nlm_Boolean Nlm_MemMapAvailable PROTO((void));

/* Initializes the memory mapping on file "name"
 * Return NULL on error
 * NOTE:  return non-NULL zero-filled structure if the file has zero length
 */
NLM_EXTERN Nlm_MemMapPtr Nlm_MemMapInit PROTO((const Nlm_Char PNTR name));

/* Closes the memory mapping. */
NLM_EXTERN void Nlm_MemMapFini PROTO((Nlm_MemMapPtr mem_mapp));


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