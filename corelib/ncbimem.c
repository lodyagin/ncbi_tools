/*   ncbimem.c
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
* File Name:  ncbimem.c
*
* Author:  Gish, Kans, Ostell, Schuler
*
* Version Creation Date:   6/4/91
*
* $Revision: 6.14 $
*
* File Description:
*   	portable memory handlers for Mac, PC, Unix
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: ncbimem.c,v $
* Revision 6.14  1998/09/25 15:57:31  vakatov
* + #include <sys/fcntl.h> to provide O_RDONLY
*
* Revision 6.13  1998/08/21 16:37:42  vakatov
* [MMAP_AVAIL]  Nlm_MemMapInit():  use "open()/close()" instead of
* "fopen()/fclose()" to avoid problems when more than FOPEN_MAX files
* are opened at the same time
*
* Revision 6.12  1998/06/25 19:40:56  vakatov
* MemMapInit():  "FileLengthEx()" to catch non-existent files
*
* Revision 6.11  1998/06/23 16:30:02  vakatov
* MemMapInit(): return zero-filled structure(not a NULL) if file is empty
*
* Revision 6.10  1998/06/11 21:02:27  shavirin
* Fixed few warnings
*
* Revision 6.9  1998/05/28 19:10:45  vakatov
* Nlm_MemMapInit():  allow one to open a zero-length file
*
* Revision 6.8  1998/04/28 14:35:02  vakatov
* #define MAP_FAILED -- if not defined yet(e.g. for IRIX5)
*
* Revision 6.7  1998/04/27 21:40:19  vakatov
* Nlm_MemMapInit():  check for the mapping error(MAP_FAILED)
*
* Revision 6.6  1998/03/26 17:54:10  vakatov
* Added Nlm_CallocViaMalloc to replace native "calloc" on Solaris
*
* Revision 6.5  1998/03/20 17:09:37  vakatov
* [OS_UNIX] Added Nlm_SetHeapLimit() //NOTE: has no effect on OSF1(?)
*
* Revision 6.4  1997/12/12 20:32:37  kans
* DisposHandle now DisposeHandle on Mac
*
* Revision 6.3  1997/11/26 21:26:21  vakatov
* Fixed errors and warnings issued by C and C++ (GNU and Sun) compilers
*
* Revision 6.2  1997/10/29 02:43:48  vakatov
* Type castings to pass through the C++ compiler
*
* Revision 6.1  1997/09/09 23:42:19  vakatov
* struct Nlm_MemMap::  made "hMap" field be WIN32-only;  removed "hFile" field
* Nlm_MemMapInit():: use '/' instead of '\' in the file-mapping object
* name to allow specify(and then -- reuse) "pathed" filenames
*
* Revision 5.6  1997/07/22 19:11:34  vakatov
* Separated Main() from GetArg[svc]() functions;  [WIN_MSWIN] converged
* console and GUI libraries; [for WIN32-DLL] encapsulated global variables
*
* Revision 5.5  1997/01/06 22:28:04  vakatov
* [WIN16,WIN32][_DEBUG]  Nlm_MemFreeTrace(incl. [_CONSOLE], excl. [!_DEBUG])
*
* Revision 5.4  1996/12/13  20:35:32  epstein
* logic fix for WIN32
*
* Revision 5.3  1996/12/13  19:13:01  epstein
* use MMAP_AVAIL definition instead of cumbersome #defines
*
* Revision 5.2  1996/12/03  21:48:33  vakatov
* Adopted for 32-bit MS-Windows DLLs
*
* Revision 5.1  1996/08/29  20:50:44  madden
* Added functions for memory-mapping.
*
* 2/11/91  Kans        MemMore now frees old memory block
* 3/12/91  Kans        MemGet, HandGet will not clear if allocation failed
* 6/4/91   Kans        Macintosh version uses malloc, realloc, free
* 6/6/91   Schuler     New versions of MemGet, MemMore, MemFree for Windows
* 8/12/91  Ostell      Protection from NULL ptrs to MemFill, MemCopy
* 09-19-91 Schuler     Modified for closer resemblance to ANSI functions
* 09-19-91 Schuler     Changed all functions to _cdecl calling convention
* 04-15-93 Schuler     Changed _cdecl to LIBCALL
* 05-21-93 Schuler     win_Free checks now for invalid pointers
* 05-21-93 Schuler     Nlm_MemFreeTrace added for debugging MemFree
* 06-14-93 Schuler     Added dll_Malloc and dll_Free
* 12-20-93 Schuler     Converted ERRPOST to ErrPostEx
* ==========================================================================
*/

#define THIS_MODULE g_corelib
#define THIS_FILE  _this_file

#include <ncbi.h>

#ifdef OS_UNIX
#define USE_SETHEAPLIMIT
#endif

#ifdef USE_SETHEAPLIMIT
#include <ncbithr.h>
#endif

#include <ncbiwin.h>

/* Used for UNIX memory-mapping. */
#ifdef MMAP_AVAIL
#include <sys/mman.h>
#include <sys/fcntl.h>
#ifndef MAP_FAILED
#define MAP_FAILED ((void *) -1)
#endif
#endif

#ifdef USE_SETHEAPLIMIT
#include <sys/resource.h>
#endif



#ifdef WIN16

struct heapnode
{
    struct heapnode *next;
    HGLOBAL hSeg;
    WORD    wSeg;
    WORD    wFree;
    int     nItems;
    int     sn;
};

size_t	g_nMemHeap = ((size_t)4*KBYTE);
#define HEAPSIZE	g_nMemHeap
size_t	g_nMemThresh = ((size_t)128);
#define THRESHOLD	g_nMemThresh

struct heapnode  headFirst;
struct heapnode *headPtr = &headFirst;

static BOOL NEAR  HeapNew(struct heapnode *pNode);
static VOID NEAR  HeapFree(struct heapnode *pNode);

#endif  /* WIN16 */

short	g_bBadPtr;

extern char *g_corelib;
static char * _this_file = __FILE__;

static char * _msgMemory  = "Ran out of memory";


NLM_EXTERN void* Nlm_CallocViaMalloc(size_t n_elem, size_t item_size)
{
  size_t size = n_elem * item_size;
  void*  ptr = Nlm_Malloc(size);
  if ( ptr )
    Nlm_MemSet(ptr, 0, size);
  return ptr;
}



#ifdef USE_SETHEAPLIMIT
static size_t s_SetHeapLimit_Curr = 0;
static size_t s_SetHeapLimit_Add  = 0;
static size_t s_SetHeapLimit_Max  = 0;
static TNlmMutex s_SetHeapLimit_Mutex;

NLM_EXTERN Nlm_Boolean Nlm_SetHeapLimit(size_t curr, size_t add, size_t max)
{
  Nlm_Boolean ok = FALSE;
  struct rlimit rl;

  if (NlmMutexLockEx(&s_SetHeapLimit_Mutex) != 0)
    return FALSE;

  if ( curr ) {
    rl.rlim_cur = curr;
    rl.rlim_max = RLIM_INFINITY;
  }
  else {
    rl.rlim_cur = rl.rlim_max = RLIM_INFINITY;
  }

  if (setrlimit(RLIMIT_DATA, &rl) == 0) {
    s_SetHeapLimit_Curr = curr;
    s_SetHeapLimit_Add  = add;
    s_SetHeapLimit_Max  = max;
    ok = TRUE;
  }

  NlmMutexUnlock(s_SetHeapLimit_Mutex);
  return ok;
}

#else
NLM_EXTERN Nlm_Boolean Nlm_SetHeapLimit(size_t curr, size_t add, size_t max) {
  return FALSE;
}
#endif /* USE_SETHEAPLIMIT */


typedef enum {
  eA_Malloc,
  eA_Calloc,
  eA_Realloc
} EAllocator;

/****************************************************************************
 *
 * s_MemAllocator(ptr, size, flags, allocator)
 *   ptr       -- origin pointer(for eA_Realloc only)
 *   size      -- number of bytes to allocate
 *   flags     -- any of the following bits may be set
 *      MGET_CLEAR     clear to zeros
 *      MGET_ERRPOST   post error on allocaion failure
 *   allocator -- method
 *
 * It is a generic routine for:
 *   Nlm_MemGet(size, flags)
 *   Nlm_MemNew(size)
 *   Nlm_MemMore(ptr, size)
 *   Nlm_MemExtend(ptr, size, oldsize)
 *
 ****************************************************************************/

static void* s_MemAllocator(void *ptr, size_t size,
                            unsigned int flags, EAllocator allocator)
{
  void *x_ptr = 0;

  switch ( allocator ) {
  case eA_Malloc:
    if ( !size )
      return 0;
    x_ptr = Nlm_Malloc(size);
    break;
  case eA_Calloc:
    if ( !size )
      return 0;
    x_ptr = Nlm_Calloc(size, 1);
    break;
  case eA_Realloc:
    if ( !ptr ) {
      if (flags & MGET_ERRPOST)
        ErrPostEx(SEV_WARNING, E_Programmer, 0, "Attempt to realloc NULL");
      return 0;
    }
    if ( !size )
      return Nlm_MemFree(ptr);
      
    x_ptr = Nlm_Realloc(ptr, size);
    break;
  }

#ifdef USE_SETHEAPLIMIT
  if (!x_ptr  &&  s_SetHeapLimit_Curr) {
    NlmMutexLock(s_SetHeapLimit_Mutex);

    while (s_SetHeapLimit_Curr < s_SetHeapLimit_Max) {
      struct rlimit rl;
      size_t x_curr = s_SetHeapLimit_Curr + s_SetHeapLimit_Add;
      if (x_curr > s_SetHeapLimit_Max)
        x_curr = s_SetHeapLimit_Max;

      if (flags & MGET_ERRPOST) {
        ErrPostEx(SEV_WARNING, E_NoMemory, 0,
                  "Trying to allocate %ld bytes;  "
                  "adjusting max.avail. heap size from %ld to %ld",
                  (long)size, (long)s_SetHeapLimit_Curr, (long)x_curr);
      }

      rl.rlim_cur = x_curr;
      rl.rlim_max = RLIM_INFINITY;
      if (setrlimit(RLIMIT_DATA, &rl) != 0)
        break;

      s_SetHeapLimit_Curr = x_curr;

      switch ( allocator ) {
      case eA_Malloc:
        x_ptr = Nlm_Malloc(size);
        break;
      case eA_Calloc:
        x_ptr = Nlm_Calloc(size, 1);
        break;
      case eA_Realloc:
        x_ptr = Nlm_Realloc(ptr, size);
        break;
      }
      if ( x_ptr )
        break;
    }

    NlmMutexUnlock(s_SetHeapLimit_Mutex);
  }
#endif /* USE_SETHEAPLIMIT */

  if ( x_ptr ) {
    if (flags & MGET_CLEAR)
      memset(x_ptr, 0, size);
  }
  else if (flags & MGET_ERRPOST) {
    ErrPostEx(SEV_FATAL, E_NoMemory, 0,
              "Failed to allocate %ld bytes", (long)size);
  }

  return x_ptr;
}


NLM_EXTERN void* LIBCALL Nlm_MemGet(size_t size, unsigned int flags)
{
  return s_MemAllocator(0, size, flags, eA_Malloc);
}

NLM_EXTERN void* LIBCALL Nlm_MemNew(size_t size)
{
  return s_MemAllocator(0, size, MGET_ERRPOST, eA_Calloc);
}

NLM_EXTERN void* LIBCALL Nlm_MemMore(void *ptr, size_t size)
{
  return s_MemAllocator(ptr, size, MGET_ERRPOST, eA_Realloc);
}

NLM_EXTERN void* LIBCALL Nlm_MemExtend(void *ptr, size_t size, size_t oldsize)
{
	void *x_ptr = s_MemAllocator(ptr, size, MGET_ERRPOST, eA_Realloc);
	if (x_ptr  &&  size > oldsize)
		memset((char*)x_ptr + oldsize, 0, size - oldsize);

	return x_ptr;
}


/*****************************************************************************
*
*   Nlm_MemFree(ptr)
*   	frees allocated memory
*
*****************************************************************************/

NLM_EXTERN void * LIBCALL  Nlm_MemFree (void *ptr)
{
	if (ptr != NULL)
		Nlm_Free (ptr);

    return NULL;
}


#if defined(_DEBUG)  &&  defined(OS_MSWIN)
NLM_EXTERN void * LIBCALL  Nlm_MemFreeTrace (void *ptr, const char *module,
			const char *filename, int linenum)
{
	if (ptr != NULL)
	{
		Nlm_Free(ptr);
		if (g_bBadPtr)
		{
			Nlm_ErrSetContext(module,filename,linenum,TRUE,0,0,0);
			Nlm_ErrPostEx(SEV_WARNING,E_Programmer,0,
						"MemFree: attempt to free invalid pointer");
		}
	}
	return NULL;
}
#endif


/*****************************************************************************
*
*    void Nlm_MemCopy(Pointer to, Pointer from, Uint4 bytes)
*       WARNING: no check on overlapping regions
*
*****************************************************************************/

NLM_EXTERN void * LIBCALL  Nlm_MemCopy (void *dst, const void *src, size_t bytes)
{
    return (dst&&src) ? Nlm_MemCpy (dst, src, bytes) : NULL;
}

/*****************************************************************************
*
*    void Nlm_MemDup (Pointer orig, Uint4 bytes)
*       Duplicate the region of memory pointed to by 'orig' for 'size' length
*
*****************************************************************************/

NLM_EXTERN void * LIBCALL  Nlm_MemDup (const void *orig, size_t size)
{
	Nlm_VoidPtr	copy;

	if (orig == NULL || size == 0)
		return NULL;

	if ((copy = Nlm_Malloc (size)) == NULL)
		ErrPostEx(SEV_FATAL,E_NoMemory,0,_msgMemory);

	Nlm_MemCpy(copy, orig, size);
		return copy;
}

/*****************************************************************************
*
*    void Nlm_MemMove (Pointer to, Pointer from, Uint4 bytes)
*       This code will work on overlapping regions
*
*****************************************************************************/

NLM_EXTERN void * LIBCALL  Nlm_MemMove (void * dst, const void *src, size_t bytes)
{
	register char *dest = (char*) dst;
	register const char *sorc = (const char*) src;

	if (dest > sorc) {
		sorc += bytes;
		dest += bytes;
		while (bytes-- != 0) {
			*--dest = *--sorc;
		}
		} else {
		while (bytes-- != 0) {
			*dest++ = *sorc++;
		}
	}
	return dst;
}

/*****************************************************************************
*
*   void Nlm_MemFill(to, value, bytes)
*   	set a block of memory to a value
*
*****************************************************************************/

NLM_EXTERN void * LIBCALL  Nlm_MemFill (void *buf, int value, size_t bytes)
{
    return  buf ? Nlm_MemSet (buf, value, bytes) : NULL;
}


#if defined(OS_MAC) || defined(OS_MSWIN) || defined(MSC_VIRT)
/***** Handle functions are for Macintosh and Windows only *****/
/***** or Microsoft virtual memory manager ****/

static char * _msgNullHnd = "NULL handle passed as an argument";

#ifdef MSC_VIRT
Nlm_Boolean wrote_to_handle;   /* used by ncbibs write routines */
#endif

/*****************************************************************************
*
*   Nlm_HandGet(size, clear_out)
*   	returns handle to allocated memory
*       if (clear_out) clear memory to 0
*
*****************************************************************************/

NLM_EXTERN Nlm_Handle LIBCALL  Nlm_HandGet (size_t size, Nlm_Boolean clear_out)

{
    Nlm_VoidPtr ptr;
    Nlm_Handle  hnd;

    if (size == 0) return NULL;

#ifdef OS_MAC
    hnd = (Nlm_Handle) NewHandle (size);
#endif

#ifdef OS_MSWIN
    hnd = (Nlm_Handle) GlobalAlloc (GMEM_MOVEABLE, size);
#endif

#ifdef MSC_VIRT
	hnd = (Nlm_Handle) _vmalloc ((unsigned long)size);
#endif

    if (hnd == NULL)
    	ErrPostEx(SEV_FATAL,E_NoMemory,0,_msgMemory);

    else if (clear_out)	{
#ifdef MSC_VIRT
		wrote_to_handle = TRUE;
#endif
        if ((ptr = HandLock (hnd)) != NULL)
            Nlm_MemSet (ptr, 0, size);
        HandUnlock (hnd);
    }

    return  hnd;
}


/*****************************************************************************
*
*   Nlm_HandNew(size)
*
*****************************************************************************/

NLM_EXTERN Nlm_Handle LIBCALL  Nlm_HandNew (size_t size)

{
    Nlm_Handle  hnd;

    if (size == 0)  return NULL;

    if ((hnd = HandGet (size, TRUE)) == NULL)
    	ErrPostEx(SEV_FATAL,E_NoMemory,0,_msgMemory);

    return hnd;
}

/*****************************************************************************
*
*   Nlm_HandMore(hnd, size)
*
*****************************************************************************/

NLM_EXTERN Nlm_Handle LIBCALL  Nlm_HandMore (Nlm_Handle hnd, size_t size)

{
    Nlm_Handle  hnd2;

	if (size == 0) {
		Nlm_HandFree (hnd);
		return NULL;
	}

    if (hnd == NULL) {
    	ErrPostEx(SEV_WARNING,E_Programmer,0,"HandMore: %s", _msgNullHnd);
		return NULL;
    }

#ifdef OS_MAC
    SetHandleSize ((Handle)hnd, (Size)size);
    hnd2 = hnd;
    if (MemError() != noErr)
        hnd2 = NULL;
#endif

#ifdef OS_MSWIN
    hnd2 = (Nlm_Handle) GlobalReAlloc ((HANDLE)hnd, size, GHND);
#endif

#ifdef MSC_VIRT
	hnd2 = (Nlm_Handle) _vrealloc ((_vmhnd_t)hnd, (unsigned long)size);
#endif

    if (hnd2 == NULL)
    	ErrPostEx(SEV_FATAL,E_NoMemory,0,_msgMemory);

    return  hnd2;
}


/*****************************************************************************
*
*   Nlm_HandFree (hnd)
*
*****************************************************************************/

NLM_EXTERN Nlm_Handle LIBCALL  Nlm_HandFree (Nlm_Handle hnd)
{
#ifdef MSC_VIRT
	_vmhnd_t x;
#endif

    if (hnd) {

#ifdef OS_MAC
        DisposeHandle ((Handle) hnd);
#endif

#ifdef OS_MSWIN
        GlobalFree ((HANDLE) hnd);
#endif

#ifdef MSC_VIRT
		x = (_vmhnd_t)hnd;
		while (_vlockcnt(x))
			_vunlock(x, _VM_CLEAN);
		_vfree(x);
#endif
    }
    else
    	ErrPostEx(SEV_WARNING,E_Programmer,0,"HandFree: %s", _msgNullHnd);

    return NULL;
}


/*****************************************************************************
*
*   Nlm_HandLock (hnd)
*
*****************************************************************************/

NLM_EXTERN Nlm_VoidPtr LIBCALL  Nlm_HandLock (Nlm_Handle hnd)
{
    Nlm_VoidPtr ptr;

    if (hnd == NULL) {
    	ErrPostEx(SEV_WARNING,E_Programmer,0,"HandLock: %s", _msgNullHnd);
        return NULL;
    }

#ifdef OS_MAC
    HLock ((Handle) hnd);
    ptr = *((Handle) hnd);
#endif

#ifdef OS_MSWIN
    ptr = GlobalLock ((HANDLE) hnd);
#endif

#ifdef MSC_VIRT
	ptr = _vlock((_vmhnd_t) hnd);
#endif

    return  ptr;
}

/*****************************************************************************
*
*   Nlm_HandUnlock(hnd)
*
*****************************************************************************/

NLM_EXTERN Nlm_VoidPtr LIBCALL  Nlm_HandUnlock (Nlm_Handle hnd)
{
#ifdef MSC_VIRT
	int dirty = _VM_CLEAN;
#endif

    if (hnd == NULL)
    	ErrPostEx(SEV_WARNING,E_Programmer,0,"HandUnlock: %s", _msgNullHnd);
    else {
#ifdef OS_MAC
        HUnlock ((Handle) hnd);
#endif

#ifdef OS_MSWIN
        GlobalUnlock ((HANDLE) hnd);
#endif

#ifdef MSC_VIRT
		if (wrote_to_handle == TRUE)
			dirty = _VM_DIRTY;
		_vunlock ((_vmhnd_t) hnd, dirty);  /* always assume dirty */
		wrote_to_handle = FALSE;
#endif
    }

    return NULL;
}

#endif /* Mac or Win */




#ifdef WIN_MAC
#ifdef USE_MAC_MEMORY
/*****************************************************************************
*
*   Macintosh-specific functions (using Mac Toolbox calls):
*
*   mac_Malloc     Macintosh version of malloc
*   mac_Calloc     Macintosh version of calloc
*   mac_Realloc    Macintosh version of realloc
*   mac_Free       Macintosh version of free
*
*****************************************************************************/

void *mac_Malloc (size_t size)

{
    return NewPtr (size);
}

void *mac_Calloc (size_t nmemb, size_t size)

{
    return NewPtrClear (nmemb * size);
}

void *mac_Realloc (void *ptr, size_t size)

{
    void *ptr2;

    ptr2 = NewPtrClear (size);
    if (ptr2 != NULL && ptr != NULL) {
      Nlm_MemCpy (ptr2, ptr, size);
      DisposPtr (ptr);
    }
    return ptr2;
}

void mac_Free (void *ptr)

{
    DisposPtr (ptr);
}
#endif
#endif



#ifdef WIN16
/*****************************************************************************
*
*   Windows-specific functions:
*
*   win16_Malloc     Windows version of malloc
*   win16_Calloc     Windows version of calloc
*   win16_Realloc    Windows version of realloc
*   win16_Free       Windows version of free
*
*   The following fucntions are not exported.
*
*   HeapNew     allocates an auxilliary local heap
*   HeapFree    frees an auxilliary local heap
*
*****************************************************************************/

LPVOID LIBCALL  win16_Malloc (size_t size)
{
	struct heapnode *pNode;
	WORD      wSeg, wSize, wFree, wOfs;
	HGLOBAL   hSeg;
	LPVOID    pMem;

    if (size == 0)  return NULL;

   /*
    * Two allocation strategies are used depending on the reguested
    * size of the new block.  Large objects (size >= THRESHOLD) are
    * allocated as independent blocks from the global heap in the
    * usual way and then locked.  Small objects (lSize < THRESHOLD)
    * are allocated from auxilliary local heaps that are created as
    * needed.  References to these heaps are maintained in a linked
    * list structure (located in the default data segment).
    */

    if (size < THRESHOLD)
    {
        wSize = (WORD)size;

        for (pNode=headPtr; pNode; pNode=pNode->next)
        {

           /*
            * If the current node is empty, we allocate a new heap.
            * On failure, break out of the loop to do an error return.
            */

            if ((wSeg = pNode->wSeg) ==0)
            {
                if (!HeapNew (pNode))  break;
                wSeg = pNode->wSeg;
            }

            if (wSize <= pNode->wFree)
            {

               /*
                * Attempt to allocate memory from this node's heap.  To do
                * this, we need to switch the DS register to the heap's
                * segment, then do a LocalAlloc, and of course restore
                * the DS register afterwards (some assembly required).
                * If this is successful, we construct a far pointer to the
                * allocated block and return it as the function result.
                * Otherwise, we continue around the cycle...
                */

                _asm {
                	push  ds
                    mov   ax, wSeg
                    mov   ds, ax
                }

                wOfs = (WORD) LocalAlloc(LMEM_FIXED, wSize);
                wFree = LocalCompact(0);

                _asm {
               		pop   ds
                }

                if (wOfs)
                {
                    pNode->nItems +=1;
                    pNode->wFree = wFree;
                    pMem = (LPVOID) MAKELONG(wOfs, wSeg);
                    return  pMem;
                }
            }

           /*
            * If we are at the end of the chain, create a new node and
            * add it to the chain.  The new node becomes the current
            * node in the next loop iteration.  If the allocation fails
            * we will drop out of the loop and do an error return.
            */

            if ( ! pNode->next )
			{
                pNode->next = (struct heapnode*) calloc(sizeof(struct heapnode),1);
				pNode->next->sn = pNode->sn +1;
			}
        }

        return  NULL;
    }
	else 
	{
		/* 
	    	 *  Large memory objects are allocated from the global heap as
		 *  independent moveable blocks.  They are then locked to generate 
		 *  the far pointer and they remain locked until freed.
		 */

	    if ((hSeg = GlobalAlloc(GMEM_MOVEABLE,size)) != NULL) 
    	{
        	pMem = GlobalLock(hSeg);
			ASSERT(OFFSETOF(pMem)==0);
        	return  pMem;
	    }
	}
    return NULL;
}


LPVOID LIBCALL  win16_Calloc (size_t items, size_t size)
{
    LPVOID ptr;
    unsigned long bytes = (unsigned long) items * size;

     if (bytes > (unsigned long) SIZE_MAX)  return NULL;

    if ((ptr = win16_Malloc((size_t) bytes)) != NULL)
        _fmemset(ptr,0,(size_t) bytes);

    return ptr;
}


LPVOID LIBCALL  win16_Realloc (LPVOID pMem, size_t size)
{
    WORD wSeg1 = SELECTOROF(pMem);
    HGLOBAL   hSeg1, hSeg2;

    if (size ==0) 
    {
        win16_Free(pMem);
        return NULL;
    }

	if (OFFSETOF(pMem) ==0) 
	{
       /*
        *  The existing block was not allocated from one of the
        *  auxilliary local heaps, so it must be an independent
        *  global block.  Use the normal GlobalReAlloc strategy.
        */

		hSeg1 = (HGLOBAL) LOWORD(GlobalHandle(wSeg1));
		GlobalUnlock (hSeg1);
		if ((hSeg2 = GlobalReAlloc(hSeg1,size,GMEM_MOVEABLE)) != NULL)
			return GlobalLock(hSeg2);
			/*
		else
			pMem = NULL;
		return  pMem;
		*/
	}
	else
	{
	    struct heapnode *pNode;
	    WORD      wSize1, wOfs1;
	    WORD      wSize2, wOfs2, wFree;
	    LPVOID    pMem2;
        
       /*
        *  Scan the linked list to see if the segment of the pointer
        *  to be freed matches any of the auxilliary local heaps.
        *  If so, LocalFree from that heap.  Otherwise we assume that
        *  it is a normal global memory object and GlobalFree.
        */

        for (pNode=headPtr; pNode; pNode=pNode->next) 
        {
            if (pNode->wSeg == wSeg1)
			{
         	   wSize2 = (WORD) size;
            	wOfs1 = LOWORD((DWORD)pMem);

	            _asm {
    	        	push  ds
        	        mov   ax, wSeg1
            	    mov   ds, ax   
	            }

    	        wSize1 = LocalSize((HLOCAL)wSeg1);
        	    wOfs2 = (WORD) LocalReAlloc((HLOCAL)wOfs1,wSize2,LMEM_MOVEABLE);
            	wFree = LocalCompact(0);

	            _asm { 
    	        	pop   ds 
        	    }

            	if (wOfs2)     /* successfully realloc'ed within same heap */
				{
    	            pNode->wFree = wFree;
        	        pMem = (LPVOID) MAKELONG(wOfs2,wSeg1);
            	    return  pMem;
	            }

    	       /*
        	    *  We couldn't realloc the object within the same heap, so
            	*  we will now attempt to allocate a new pointer and copy
	            *  contents of the existing block to the new one.  If
    	        *  successful, the old block is freed.
        	    */

            	if ((pMem2 = win16_Malloc(size)) != NULL) 
	            {
    	            _fmemcpy(pMem2,pMem,MIN(wSize1,wSize2));
        	        win16_Free(pMem);
            	}
	            return  pMem2;
    	    }
		}
    }

    return NULL;
}


VOID LIBCALL  win16_Free (LPVOID pMem)
{
	WORD wSeg = SELECTOROF(pMem);
	WORD wOfs = OFFSETOF(pMem);

	g_bBadPtr = FALSE;

	if (wOfs ==0)
	{
		HGLOBAL hSeg = (HGLOBAL) LOWORD(GlobalHandle(wSeg));
		if (GlobalSize(hSeg))
		{
			GlobalUnlock(hSeg);
			GlobalFree(hSeg);
			return;
		}
	}
	else
	{
		struct heapnode *pNode;
        
        /*
         *  Scan the linked list to see if the segment of the pointer
         *  to be freed matches any of the auxilliary local heaps.
         *  If so, LocalFree from that heap.  Otherwise we assume that
         *  it is a normal global memory object and GlobalFree.
         */

        for (pNode=headPtr; pNode; pNode=pNode->next) 
        {
            if (pNode->wSeg == wSeg)
			{
		    	WORD wFree;

	           	_asm {
	           		push  ds
    	       		mov   ax, wSeg
        	        mov   ds, ax   
           		}

				if (LocalSize((HLOCAL)wOfs))
				{
		            wOfs = (WORD) LocalFree((HLOCAL)wOfs);
    		        wFree = LocalCompact(0);
				}

        	   	_asm { 
           			pop   ds 
	           	}

    	        if (wOfs == NULL)         /* NULL means success here */
				{
            	    pNode->nItems -= 1;
                	pNode->wFree = wFree;
	                if (pNode->nItems ==0)  
					{
        	            HeapFree(pNode);  /* free the heap, but keep the node */
					}
					return;
	            }
	            goto BAD_POINTER;
			}
        }
	    goto BAD_POINTER;
    }
    
    
BAD_POINTER :
	g_bBadPtr = TRUE;
	TRACE("win16_Free(%04X:%04X)  *** invalid pointer ***\n",wSeg,wOfs);
}


static BOOL NEAR  HeapNew (struct heapnode *pNode)
{
	HGLOBAL hMem;

	pNode->nItems = 0;
	pNode->hSeg = NULL;
	pNode->wSeg = NULL;
    if ((hMem = GlobalAlloc(GMEM_MOVEABLE,HEAPSIZE)) != NULL) 
    {
	    LPSTR pMem = GlobalLock(hMem);
        WORD  wSeg = SELECTOROF(pMem);
        WORD  wSize = (WORD) GlobalSize(hMem) -16;

        if (LocalInit(wSeg,0,wSize)) 
        {
            pNode->hSeg = hMem;
            pNode->wSeg = wSeg;
            pNode->wFree = wSize;
	        GlobalUnlock(hMem);    /* remove LocalInit's lock */

			return TRUE;
        }
		GlobalUnlock(hMem);
		GlobalFree(hMem);	
    }
	TRACE("HeapNew: ** unable to allocate %ld bytes\n",HEAPSIZE);
    return FALSE;
}

static VOID NEAR  HeapFree (struct heapnode *pNode)
{   
    HGLOBAL hMem;

    if (pNode->nItems==0 && ((hMem=pNode->hSeg) != NULL)) 
    {
        GlobalUnlock(hMem);
        GlobalFree(hMem);
        pNode->hSeg = NULL;
        pNode->wSeg = 0;
    }
}


#endif /* WIN16 */


#ifdef _WINDLL
/*****************************************************************************
*
*   Windows DLL-specific functions (shared memory)
*   
*   dll_Malloc
*   dll_Calloc	(not yet)
*   dll_Realloc	(not yet)
*   dll_Free   
*   
*****************************************************************************/


void * dll_Malloc (size_t bytes)
{
	HGLOBAL hMem;
	void *pMem;

	if (bytes >0 && (hMem = GlobalAlloc(GMEM_DDESHARE,bytes)))
	{
		if (pMem = GlobalLock(hMem))
			return pMem;
		else
			GlobalFree(hMem);
	}

	TRACE("dll_Malloc(%ld) failed\n",bytes);
	return NULL;
}

void   dll_Free (void *pMem)
{
	HGLOBAL hMem;
#ifdef WIN16
	DWORD dwHandle = GlobalHandle(SELECTOROF(pMem));
	hMem = (HGLOBAL)LOWORD(dwHandle);
#else
	hMem = GlobalHandle(pMem);
#endif

	GlobalUnlock(hMem);
	GlobalFree(hMem);
}

#endif


/*********************************************************************
*	Function to test whether memory-mapping is available. 
*
*	returns TRUE if it is supported by NCBI routines.
*********************************************************************/

NLM_EXTERN Nlm_Boolean Nlm_MemMapAvailable(void)
{
#if defined(MMAP_AVAIL) || defined(WIN32)
  return TRUE;
#else
  return FALSE;
#endif
}


NLM_EXTERN Nlm_MemMapPtr Nlm_MemMapInit(const Nlm_Char PNTR name)
{
  Nlm_MemMapPtr mem_mapp;
  if (!Nlm_MemMapAvailable()  ||  !name  ||  !*name  ||
      (mem_mapp = (Nlm_MemMapPtr)Nlm_MemNew(sizeof(Nlm_MemMap))) == NULL)
    return NULL;

  for (;;) {{ /* (quasi-TRY block) */
    if ((mem_mapp->file_size = Nlm_FileLengthEx(name)) < 0)
      break;

    if (mem_mapp->file_size == 0) /* Special case */
      return mem_mapp;
	
#ifdef WIN32
    {{
      char x_name[MAX_PATH], *str;
      Nlm_StringNCpy_0(x_name, name, sizeof(x_name));
      for (str = x_name;  *str;  str++)
        if (*str == '\\')
          *str = '/';  /* name of a file-mapping object cannot contain '\' */

      if ( !(mem_mapp->hMap =
             OpenFileMapping(FILE_MAP_READ, FALSE, x_name)) ) 
        { /* If failed to attach to an existing file-mapping object then
           * create a new one(based on the specified file) */
          HANDLE hFile= CreateFile(name, GENERIC_READ, FILE_SHARE_READ, NULL,
                                   OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
          if (hFile == INVALID_HANDLE_VALUE)
            break;

          mem_mapp->hMap = CreateFileMapping(hFile, NULL, PAGE_READONLY,
                                             0, 0, x_name);
          CloseHandle( hFile );
          if ( !mem_mapp->hMap )
            break;
        }

      if ( !(mem_mapp->mmp_begin =
             MapViewOfFile(mem_mapp->hMap, FILE_MAP_READ,
                           0, 0, mem_mapp->file_size)) ) {
        CloseHandle( mem_mapp->hMap );
        break;
      }
    }}

#elif defined(MMAP_AVAIL)
    {{  /* UNIX memory mapping. */
      int fd = open(name, O_RDONLY);
      if (fd < 0)
        break;

      mem_mapp->mmp_begin = mmap(NULL, mem_mapp->file_size, PROT_READ,
                                 MAP_SHARED, fd, 0);
      close(fd);
      if (mem_mapp->mmp_begin == MAP_FAILED)
        break;
    }}
#endif

    /* Success */
    return mem_mapp;
  }}

  /* Error;  cleanup */
  Nlm_MemFree(mem_mapp);
  return NULL;
}


NLM_EXTERN void Nlm_MemMapFini(Nlm_MemMapPtr mem_mapp)
{
  if ( !mem_mapp )
    return;

#ifdef WIN32
  UnmapViewOfFile( mem_mapp->mmp_begin );
  if ( mem_mapp->hMap )
    CloseHandle( mem_mapp->hMap );
#elif defined(MMAP_AVAIL)
  munmap(mem_mapp->mmp_begin, mem_mapp->file_size);
#endif

  Nlm_MemFree( mem_mapp );
}
