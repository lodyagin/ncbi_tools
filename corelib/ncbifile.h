#ifndef _NCBIFILE_
#define _NCBIFILE_

/*   ncbifile.h
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
* File Name:  ncbifile.h
*
* Author:  Gish, Kans, Ostell, Schuler
*
* Version Creation Date:   1/1/91
*
* $Revision: 6.8 $
*
* File Description:
*   	prototypes for portable file routines
*
* Modifications:
* --------------------------------------------------------------------------
*
* ==========================================================================
*/

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

#ifdef __cplusplus
extern "C" {
#endif


typedef FILE* (LIBCALLBACK *Nlm_FileOpenHook)
       (const char *filename, const char *mode);

NLM_EXTERN FILE* LIBCALL Nlm_FileOpen(const char *filename, const char *mode);
NLM_EXTERN void LIBCALL Nlm_FileClose(FILE *stream);
NLM_EXTERN size_t LIBCALL Nlm_FileRead(void *ptr, size_t size, size_t n, FILE *stream);
NLM_EXTERN size_t LIBCALL Nlm_FileWrite(const void *ptr, size_t size, size_t n, FILE *stream);
NLM_EXTERN int LIBCALL Nlm_FilePuts(const char *ptr, FILE * fp);
NLM_EXTERN char * LIBCALL Nlm_FileGets(char *ptr, size_t size, FILE * fp);
NLM_EXTERN Nlm_CharPtr LIBCALL Nlm_FileBuildPath(Nlm_CharPtr root, Nlm_CharPtr sub_path, Nlm_CharPtr filename);
NLM_EXTERN Nlm_CharPtr LIBCALL Nlm_FileNameFind(Nlm_CharPtr pathname);
NLM_EXTERN Nlm_CharPtr LIBCALL Nlm_FilePathFind(const Nlm_Char* fullname);
NLM_EXTERN Nlm_Int8 LIBCALL Nlm_FileLength(Nlm_CharPtr fileName);
/* Nlm_FileLengthEx() returns -1 if the file does not exist) */
NLM_EXTERN Nlm_Int8 LIBCALL Nlm_FileLengthEx(const Nlm_Char* fileName);
NLM_EXTERN Nlm_Boolean LIBCALL Nlm_FileRemove(Nlm_CharPtr fileName);
NLM_EXTERN Nlm_Boolean LIBCALL Nlm_FileRename(Nlm_CharPtr oldFileName, Nlm_CharPtr newFileName);
NLM_EXTERN void LIBCALL Nlm_FileCreate(Nlm_CharPtr fileName, Nlm_CharPtr type, Nlm_CharPtr creator);
NLM_EXTERN Nlm_Boolean LIBCALL Nlm_CreateDir(Nlm_CharPtr pathname);
NLM_EXTERN ValNodePtr LIBCALL Nlm_DirCatalog (Nlm_CharPtr pathname);
NLM_EXTERN Nlm_CharPtr LIBCALL Nlm_TmpNam(Nlm_CharPtr s);
NLM_EXTERN void LIBCALL Nlm_SetFileOpenHook(Nlm_FileOpenHook hook);

/* FileCache provides buffered text read for handling Unix, Mac, and DOS line endings gracefully */

typedef struct nlm_filecachedata {
  FILE         *fp;
  Nlm_Char     buf [516];
  Nlm_Int2     ctr;
  Nlm_Int2     total;
  Nlm_Int4     offset;
  Nlm_Boolean  failed;
} Nlm_FileCache, PNTR Nlm_FileCachePtr;

NLM_EXTERN Nlm_Boolean Nlm_FileCacheSetup (Nlm_FileCache PNTR fcp, FILE *fp);
NLM_EXTERN Nlm_CharPtr Nlm_FileCacheGetString (Nlm_FileCache PNTR fcp, Nlm_CharPtr str, size_t size);
NLM_EXTERN Nlm_CharPtr Nlm_FileCacheReadLine (Nlm_FileCache PNTR fcp, Nlm_CharPtr str, size_t size, Nlm_BoolPtr nonewline);
NLM_EXTERN void Nlm_FileCacheSeek (Nlm_FileCache PNTR fcp, Nlm_Int4 pos);
NLM_EXTERN Nlm_Int4 Nlm_FileCacheTell (Nlm_FileCache PNTR fcp);
NLM_EXTERN Nlm_Boolean Nlm_FileCacheFree (Nlm_FileCache PNTR fcp, Nlm_Boolean restoreFilePos);

/* general file recursion function - directory must not be empty, proc callback function must not be NULL */

typedef void (*Nlm_DirExpProc) (Nlm_CharPtr filename, Nlm_VoidPtr userdata);

NLM_EXTERN Nlm_Int4 Nlm_DirExplore (
  Nlm_CharPtr directory,
  Nlm_CharPtr filter,
  Nlm_CharPtr suffix,
  Nlm_Boolean recurse,
  Nlm_DirExpProc proc,
  Nlm_VoidPtr userdata
);

#define FileOpen Nlm_FileOpen
#define FileClose Nlm_FileClose
#define FileRead Nlm_FileRead
#define FileWrite Nlm_FileWrite
#define FilePuts Nlm_FilePuts
#define FileGets Nlm_FileGets
#define FileBuildPath Nlm_FileBuildPath
#define FileNameFind Nlm_FileNameFind
#define FilePathFind Nlm_FilePathFind
#define FileLength Nlm_FileLength
#define FileLengthEx Nlm_FileLengthEx
#define FileRemove Nlm_FileRemove
#define FileRename Nlm_FileRename
#define FileCreate Nlm_FileCreate
#define CreateDir Nlm_CreateDir
#define DirCatalog Nlm_DirCatalog
#define TmpNam Nlm_TmpNam

#define FileCache Nlm_FileCache
#define FileCacheSetup Nlm_FileCacheSetup
#define FileCachePtr Nlm_FileCachePtr

#define FileCacheSetup Nlm_FileCacheSetup
#define FileCacheGetString Nlm_FileCacheGetString
#define FileCacheReadLine Nlm_FileCacheReadLine
#define FileCacheSeek Nlm_FileCacheSeek
#define FileCacheTell Nlm_FileCacheTell
#define FileCacheFree Nlm_FileCacheFree

#define DirExpProc Nlm_DirExpProc
#define DirExplore Nlm_DirExplore

#define EjectCd(sVolume, deviceName, rawDeviceName, mountPoint, mountCmd)  FALSE
#define MountCd(sVolume, deviceName, mountPoint, mountCmd)                 FALSE


#ifdef __cplusplus
} /* extern "C" */
#endif

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif /* _NCBIFILE_ */
