/*   ncbifile.c
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
* File Name:  ncbifile.c
*
* Author:  Gish, Kans, Ostell, Schuler
*
* Version Creation Date:   3/4/91
*
* $Revision: 6.12 $
*
* File Description: 
*     portable file routines
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
* 04-15-93 Schuler     Changed _cdecl to LIBCALL
* 12-20-93 Schuler     Converted ErrPost to ErrPostEx
* 11-27-94 Ostell      moved includes to ncbiwin.h to avoid conflict MSC
*
* $Log: ncbifile.c,v $
* Revision 6.12  1998/06/26 20:39:41  vakatov
* Added FilePathFind() -- complimentary to FileNameFind()
*
* Revision 6.11  1998/06/26 14:11:56  madden
* tempnam called with NULLs rather than empty strings
*
* Revision 6.10  1998/06/25 19:39:17  vakatov
* Added "FileLengthEx()" -- it returns -1(not 0!) if the file does not exist
*
* Revision 6.9  1998/05/28 15:59:55  vakatov
* [WIN16-BORLAND]  Nlm_DirCatalog:  "_find_t" -> "find_t"(Borland-specific)
*
* Revision 6.8  1998/05/28 14:21:53  vakatov
* [WIN16,WIN32] Nlm_DirCatalog -- tested, fixed
*
* Revision 6.7  1998/05/27 18:50:04  kans
* DirCat UNIX version redirects stderr to /dev/null to avoid printing unwanted error message when no directory exists
*
* Revision 6.6  1998/05/27 11:43:24  kans
* implemented DirCatalog for WIN16
*
* Revision 6.5  1998/05/26 20:28:17  kans
* stripped newline at end of DirCat strings in UNIX
*
* Revision 6.4  1998/05/26 17:42:45  vakatov
* [WIN32] Nlm_DirCatalog -- implemented, but not tested yet
*
* Revision 6.3  1998/05/26 15:17:14  kans
* implemented DirCatalog for OS_UNIX
*
* Revision 6.2  1998/05/24 19:20:56  kans
* added Nlm_DirCatalog (Mac implementation only so far)
*
* Revision 6.1  1998/03/31 20:31:33  vakatov
* FileRead()/FileWrite(): get rid of the redundant type cast that caused
* Int4 overflow under Solaris 2.6
*
* Revision 6.0  1997/08/25 18:15:30  madden
* Revision changed to 6.0
*
* Revision 5.3  1997/07/22 19:11:26  vakatov
* Separated Main() from GetArg[svc]() functions;  [WIN_MSWIN] converged
* console and GUI libraries; [for WIN32-DLL] encapsulated global variables
*
* Revision 5.2  1997/01/14 21:57:14  vakatov
* Fixed inaccurate string copying -- <mostly potential> 1-byte exceeding of
* the string size by StringNCat;  missing terminating '\0' by StringNCpy.
*
 * Revision 5.1  1996/12/03  21:48:33  vakatov
 * Adopted for 32-bit MS-Windows DLLs
 *
 * Revision 5.0  1996/05/28  13:18:57  ostell
 * Set to revision 5.0
 *
 * Revision 4.2  1996/02/29  14:40:29  kans
 * added USE_MPW_FILE_OPEN symbol to simplify ifdefs in FileOpen
 *
 * Revision 4.1  1996/02/28  21:55:54  kans
 * "MPW" file open also used for CodeWarrior, possible path limit otherwise
 *
 * Revision 4.0  1995/07/26  13:46:50  ostell
 * force revision to 4.0
 *
 * Revision 2.41  1995/05/15  18:45:58  ostell
 * added Log line
 *
*
*
* ==========================================================================
*/

#define THIS_MODULE g_corelib
#define THIS_FILE  _this_file

#include "corepriv.h"

#ifdef OS_UNIX_SUN
#define DEFAULT_CDROM "/dev/sr0"
#define DEFAULT_RAW_CDROM "/dev/rsr0"
#endif

#ifdef PROC_MIPS
#define DEFAULT_CDROM "/dev/scsi/sc0d4l0"
#endif

#ifdef OS_UNIX
#ifndef DEFAULT_CDROM
#define DEFAULT_CDROM "/dev/cdrom"
#endif
#endif

#if defined(OS_MSWIN) || defined (OS_NT)
#ifdef COMP_MSC
#ifndef mkdir
#define mkdir _mkdir
#endif
#ifndef stat
#define stat _stat
#endif
#endif
#endif

#ifdef OS_VMS
#ifndef DEFAULT_CDROM
#define DEFAULT_CDROM "cdrom:"
#endif
#endif

extern char *g_corelib;
static char * _this_file = __FILE__;


/*****************************************************************************
*
*   Macintosh file utilities
*
*****************************************************************************/

#ifdef OS_MAC
static short Nlm_MacGetVRefNum (Nlm_CharPtr pathname, OSErr *errptr)

{
  OSErr           err;
  Nlm_Char        filename [FILENAME_MAX];
  Nlm_Char        path [256];
  HParamBlockRec  pbh;
  Nlm_CharPtr     ptr;

  memset (&pbh, 0, sizeof (HParamBlockRec));
  Nlm_StringNCpy_0(path, pathname, sizeof(path));
  ptr = Nlm_StringRChr (path, (int) DIRDELIMCHR);
  if (ptr != NULL) {
    ptr++;
    Nlm_StringNCpy_0(filename, ptr, sizeof(filename));
    *ptr = '\0';
    Nlm_CtoPstr ((Nlm_CharPtr) path);
    pbh.volumeParam.ioNamePtr = (StringPtr) path;
    pbh.volumeParam.ioVolIndex = -1;
    err = PBHGetVInfo (&pbh, FALSE);
    if (errptr != NULL) {
      *errptr = err;
    }
    return pbh.volumeParam.ioVRefNum;
  } else {
    if (errptr != NULL) {
      *errptr = noErr;
    }
    return 0;
  }
}

static long Nlm_MacGetDirID (Nlm_CharPtr pathname, short newVRefNum, OSErr *errptr)

{
  OSErr           err;
  Nlm_Char        path [256];
  CInfoPBRec      pbc;
  Nlm_CharPtr     ptr;

  memset (&pbc, 0, sizeof (CInfoPBRec));
  Nlm_StringNCpy_0(path, pathname, sizeof(path));
  ptr = Nlm_StringRChr (path, (int) DIRDELIMCHR);
  if (ptr != NULL) {
    ptr++;
    Nlm_StringNCpy_0(path, pathname, sizeof(path));
    *ptr = '\0';
    Nlm_CtoPstr ((Nlm_CharPtr) path);
    pbc.dirInfo.ioNamePtr = (StringPtr) path;
    pbc.dirInfo.ioVRefNum = newVRefNum;
    err = PBGetCatInfo (&pbc, FALSE);
    if (errptr != NULL) {
      *errptr = err;
    }
    return pbc.dirInfo.ioDrDirID;
  } else {
    if (errptr != NULL) {
      *errptr = noErr;
    }
    return 0;
  }
}

static OSErr Nlm_SetDefault (short newVRefNum, long newDirID, short *oldVRefNum, long *oldDirID)

{
  OSErr  error;

  error = HGetVol (NULL, oldVRefNum, oldDirID);
  if (error == noErr) {
    error = HSetVol (NULL, newVRefNum, newDirID);
  }
  return (error);
}

static OSErr Nlm_RestoreDefault (short oldVRefNum, long oldDirID)

{
  OSErr  error;
  short  defaultVRefNum;
  long   defaultDirID;
  long   defaultProcID;

  error = GetWDInfo (oldVRefNum, &defaultVRefNum, &defaultDirID, &defaultProcID);
  if (error == noErr) {
    if (defaultDirID != fsRtDirID) {
      error = SetVol (NULL, oldVRefNum);
    } else {
      error = HSetVol (NULL, oldVRefNum, oldDirID);
    }
  }
  return (error);
}
#endif

/*****************************************************************************
*
*   FileOpen(filename, mode)
*     if (filename == "stdin" or "stdout" or "stderr"
*           returns those predefined
*           streams on non-windowing systems)
*
*****************************************************************************/

#ifdef OS_MAC
#ifndef COMP_THINKC
#define USE_MPW_FILE_OPEN
#endif
#endif

static Nlm_FileOpenHook _hookFile = NULL;

#ifdef USE_MPW_FILE_OPEN
/*
*  MPWOptimizationErrorBypass was called in order to avoid an apparent MPC C
*  compiler optimization problem that resulted in the newDirID value sometimes
*  appearing to be 0.  Placing any debugging statement (or this dummy function)
*  after the statement that gets newDirID and before the statement that uses it
*  originally appeared to fix the problem.  Upon further investigation, it turns
*  out to be a problem with a UNIX to Mac file server product.  The problem did
*  not occur when compiling under THINK C.
*/

static Nlm_Int2 Nlm_MPWOptimizationErrorBypass (short newVRefNum, long newDirID)
{
  return 0;
}

/*
*  In MPW, if a temporary file is written first, without being created, it is
*  created on the hard disk, rather than in the appropriate path.  This code
*  creates the file in the desired location.
*/

static void Nlm_MPWCreateOutputFile (Nlm_CharPtr pathname, Nlm_CharPtr filename)

{
  OSErr     err;
  FILE      *f;
  Nlm_Char  temp [256];

  f = fopen (filename, "r");
  if (f == NULL) {
    Nlm_StringNCpy_0(temp, pathname, sizeof(temp));
    Nlm_CtoPstr ((Nlm_CharPtr) temp);
    err = Create ((StringPtr) temp, 0, '    ', 'TEXT');
  } else {
    fclose (f);
  }
}

static FILE * LIBCALL  Nlm_MPWFileOpen (Nlm_CharPtr pathname, Nlm_CharPtr mode)

{
  Nlm_Boolean  createfile;
  OSErr        err;
  FILE         *f;
  long         newDirID;
  short        newVRefNum;
  long         oldDirID;
  short        oldVRefNum;
  Nlm_CharPtr  ptr;

  newVRefNum = Nlm_MacGetVRefNum (pathname, &err);
  newDirID = Nlm_MacGetDirID (pathname, newVRefNum, &err);
  ptr = Nlm_StringRChr (pathname, (int) DIRDELIMCHR);
  createfile = (Nlm_Boolean) (strchr (mode, 'w') != NULL);
  if (ptr != NULL) {
    ptr++;
	Nlm_MPWOptimizationErrorBypass (newVRefNum, newDirID);
    err = Nlm_SetDefault (newVRefNum, newDirID, &oldVRefNum, &oldDirID);
	if (createfile) {
	  Nlm_MPWCreateOutputFile (pathname, ptr);
	}
    f = fopen (ptr, mode);
    err = Nlm_RestoreDefault (oldVRefNum, oldDirID);
  } else {
    if (createfile) {
	  Nlm_MPWCreateOutputFile (pathname, pathname);
	}
    f = fopen (pathname, mode);
  }
  return f;
}
#endif

NLM_EXTERN FILE * LIBCALL Nlm_FileOpen(const char *filename, const char *mode)
{
  FILE *f = NULL;

  if ( _hookFile )
    return _hookFile(filename, mode);

#if defined(WIN_DUMB)
    if ( Nlm_HasConsole )
      {
        if      ( !StringCmp("stdin",  filename))
          f = stdin;
        else if ( !StringCmp("stdout", filename) )
          f = stdout;
        else if ( !StringCmp("stderr", filename))
          f = stderr;
        else
          f = fopen(filename, mode);

#ifdef WIN32
        if (strchr(mode, 'b')  &&
            (f == stdin  ||  f == stdout  ||  f == stderr))
          setmode(fileno(f), O_BINARY);
#endif
      }
    else
      f = fopen(filename, mode);

#elif defined(OS_MAC)
  {{
    OSType    fCreator;
    Nlm_Int2  fError;
    FInfo     fInfo;
    OSType    fType;
    Nlm_Char  temp [256];
    Nlm_StringNCpy_0(temp, filename, sizeof(temp));
    Nlm_CtoPstr ((Nlm_CharPtr) temp);
    fError = GetFInfo ((StringPtr) temp, 0, &fInfo);
    if (fError == 0) {
      fCreator = fInfo.fdCreator;
      fType = fInfo.fdType;
    } else {
      if (strchr(mode, 'b') != NULL)
        fType = '    ';
      else
        fType = 'TEXT';
      fCreator = '    ';
    }
#ifdef USE_MPW_FILE_OPEN
    {{
      Nlm_Char localmode [16];
      Nlm_Char path [PATH_MAX];
      Nlm_StringNCpy_0(path, filename, sizeof(path));
      Nlm_StringNCpy_0(localmode, mode, sizeof(localmode));
      f = Nlm_MPWFileOpen(path, localmode);
    }}
#else
    f = fopen (filename,mode);
#endif /* USE_MPW_FILE_OPEN */

    fError = GetFInfo ((StringPtr) temp, 0, &fInfo);
    if (fError == 0) {
      fInfo.fdCreator = fCreator;
      fInfo.fdType = fType;
      fError = SetFInfo ((StringPtr) temp, 0, &fInfo);
    }
  }} /* def OS_MAC */

#elif defined(OS_VMS) && defined(DCC4DW12)
  /* never used */ 
  f = fopen (filename, mode);
  if (f  &&
      fstat(fileno(f), &statbuf) == 0  &&
      statbuf.st_fab_rfm == FAB$C_UDF)
    {
      fclose(f);
      f = fopen(filename,mode,"ctx=stm");
    }

#else
    f = fopen(filename, mode);  
#endif

  if (f == NULL)
    ErrPostEx(SEV_INFO, E_File, E_FOpen, "FileOpen(\"%s\",\"%s\") failed",
              filename, mode);
		
  return f;
}

/*****************************************************************************
*
*   SetFileOpenHook(hook)
*
*****************************************************************************/

NLM_EXTERN void LIBCALL Nlm_SetFileOpenHook (Nlm_FileOpenHook hook)
{
	_hookFile = hook;
}

/*****************************************************************************
*
*   FileClose(fp)
*
*****************************************************************************/

NLM_EXTERN void LIBCALL  Nlm_FileClose (FILE *stream)
{
  if (stream == NULL)
    return;

#ifdef WIN_DUMB    
  if (stream == stdin  ||  stream == stdout  ||  stream == stderr)
    {
#ifdef WIN32
      setmode(fileno(stream), O_TEXT);
#endif
      return;
    }
#endif

  fclose(stream);
}

/*****************************************************************************
*   FileRead(buf, size, fp)
*****************************************************************************/
#ifdef WIN16
#include <dos.h> /* dos.h defines the FP_SEG macro */
#endif

NLM_EXTERN size_t LIBCALL Nlm_FileRead
(void *ptr, size_t size, size_t n, FILE *stream)
{
  if (n  &&  (SIZE_MAX / n) < size) {
    ErrPostEx(SEV_WARNING,E_Programmer,0,"FileRead: size > SIZE_MAX");
    return 0;
  }
  if (!ptr  ||  !stream)
    return 0;

  return fread(ptr,size,n,stream);
}

/*****************************************************************************
*   FileWrite(buf, size, fp)
*****************************************************************************/
NLM_EXTERN size_t LIBCALL Nlm_FileWrite
(const void *ptr, size_t size, size_t n, FILE *stream)
{
  size_t cnt;
  if (n   &&  (SIZE_MAX / n)  <  size) {
    ErrPostEx(SEV_WARNING,E_Programmer,0,"FileWrite:  size > SIZE_MAX");
    return 0;
  }
  if (!ptr  ||  !stream)
    return 0;

  cnt = fwrite(ptr,size,n,stream);
  if (cnt != n)
    ErrPostEx(SEV_FATAL,E_File,E_FWrite,"File write error");

  return cnt;
}

/*****************************************************************************
*
*   FilePuts(ptr, fp)
*
*****************************************************************************/
NLM_EXTERN int LIBCALL  Nlm_FilePuts (const char *ptr, FILE *fp)
{
	int retval;

	if ((ptr == NULL) || (fp == NULL))
    	return EOF;
	if ((retval = fputs(ptr,fp)) ==EOF)
		ErrPostEx(SEV_FATAL,E_File,E_FWrite,"File write error");
	return retval;
}

/*****************************************************************************
*
*   FileGets()
*
*****************************************************************************/
NLM_EXTERN char * LIBCALL  Nlm_FileGets (Nlm_CharPtr ptr, size_t size, FILE *fp)
{
	if ((ptr == NULL) || (size <= 0) || (fp == NULL))
		return NULL;
	return fgets(ptr,size,fp);
}


/*****************************************************************************
*
*   FileBuildPath()
*
*****************************************************************************/
NLM_EXTERN Nlm_CharPtr LIBCALL  Nlm_FileBuildPath (Nlm_CharPtr root, Nlm_CharPtr sub_path, Nlm_CharPtr filename)

{
    Nlm_CharPtr tmp;
    Nlm_Boolean dir_start = FALSE;
#ifdef OS_VMS
  Nlm_Boolean had_root = FALSE;
#endif

    if (root == NULL)              /* no place to put it */
        return NULL;

    tmp = root;
    if (*tmp != '\0')                /* if not empty */
    {
#ifndef OS_VMS
        dir_start = TRUE;
#else
        had_root = TRUE;
#endif
        while (*tmp != '\0')
        {
#ifdef OS_VMS
            if (*tmp == '[')
                dir_start = TRUE;
#endif
            tmp++;
        }

        if ((*(tmp - 1) != DIRDELIMCHR) && (dir_start))
        {
            *tmp = DIRDELIMCHR;
            tmp++; *tmp = '\0';
        }
    }

    if (sub_path != NULL)
    {
#ifdef OS_VMS
        if (dir_start)
        {
            *(tmp-1) = '.';
            if (*sub_path == '[')
                sub_path++;
        }
        else if ((had_root) && (*sub_path != '['))
        {
            *tmp = '[';
            tmp++; *tmp = '\0';
        }
#else
        if ((dir_start) && (*sub_path == DIRDELIMCHR))
            sub_path++;
#endif
        tmp = StringMove(tmp, sub_path);
        if (*(tmp-1) != DIRDELIMCHR)
        {
            *tmp = DIRDELIMCHR;
            tmp++; *tmp = '\0';
        }
    }

    if (filename != NULL)
        StringMove(tmp, filename);

    return root;
}

/*****************************************************************************
*
*   FileNameFind()
*
*****************************************************************************/
NLM_EXTERN Nlm_CharPtr LIBCALL Nlm_FileNameFind (Nlm_CharPtr pathname)

{
  Nlm_CharPtr  filename;
  Nlm_Int2     len;

  if (pathname != NULL) {
    len = Nlm_StringLen (pathname);
    filename = &(pathname [len]);
    while (len > 0 && pathname [len - 1] != DIRDELIMCHR) {
      len--;
      filename--;
    }
    return filename;
  } else {
    return NULL;
  }
}


/*****************************************************************************
*
*   FilePathFind()
*
*****************************************************************************/
NLM_EXTERN Nlm_CharPtr LIBCALL Nlm_FilePathFind(const Nlm_Char* fullname)
{
  Nlm_CharPtr str;
  size_t	     len = Nlm_StringLen(fullname);
  if ( !len )
    return 0;

  while (len  &&  fullname[len] != DIRDELIMCHR)
    len--;

  str = (Nlm_Char*)Nlm_MemGet(len + 1, MGET_ERRPOST);
  Nlm_MemCpy(str, fullname, len);
  str[len] = '\0';
  return str;
}


/*****************************************************************************
*
*   FileLength()
*
*****************************************************************************/
NLM_EXTERN Nlm_Int4 LIBCALL Nlm_FileLength(Nlm_CharPtr fileName)
{
  Nlm_Int4 file_len = Nlm_FileLengthEx(fileName);
  return (file_len > 0) ? file_len : 0;
}


/*****************************************************************************
*
*   FileLengthEx()
*
*****************************************************************************/
NLM_EXTERN Nlm_Int4 LIBCALL Nlm_FileLengthEx(const Nlm_Char* fileName)
{
  if (!fileName  ||  !*fileName)
    return -1;

#ifdef OS_MAC
  {{
    ParamBlockRec params;
    Nlm_Char      path[256];

    Nlm_StringNCpy_0(path, fileName, sizeof(path));
    params.fileParam.ioNamePtr = (StringPtr)path;
    params.fileParam.ioVRefNum = 0;
    params.fileParam.ioFDirIndex = 0;
    Nlm_CtoPstr((Nlm_CharPtr) path);
    return (PBGetFInfo(&params, FALSE) == noErr) ?
      params.fileParam.ioFlLgLen : -1;
  }}
#else
  {{
    struct stat sbuf;
    return (stat(fileName, &sbuf) == 0) ? sbuf.st_size : -1;
  }}
#endif
}


/*****************************************************************************
*
*   FileDelete()
*
*****************************************************************************/
NLM_EXTERN Nlm_Boolean LIBCALL Nlm_FileRemove (Nlm_CharPtr fileName)

{
  Nlm_Char  local [256];

  if (fileName != NULL && fileName [0] != '\0') {
    Nlm_StringNCpy_0(local, fileName, sizeof(local));
    return (Nlm_Boolean) (remove (local) == 0);
  } else {
    return FALSE;
  }
}

/*****************************************************************************
*
*   FileRename()
*
*****************************************************************************/
NLM_EXTERN Nlm_Boolean LIBCALL Nlm_FileRename (Nlm_CharPtr oldFileName, Nlm_CharPtr newFileName)

{
  Nlm_Char  localnew [256];
  Nlm_Char  localold [256];

  if (oldFileName != NULL && oldFileName [0] != '\0'
    && newFileName != NULL && newFileName [0] != '\0') {
    Nlm_StringNCpy_0(localold, oldFileName, sizeof(localold));
    Nlm_StringNCpy_0(localnew, newFileName, sizeof(localnew));
    return (Nlm_Boolean) (rename (localold, localnew) == 0);
  } else {
    return FALSE;
  }
}

/*****************************************************************************
*
*   FileCreate()
*
*****************************************************************************/
#ifdef WIN_MAC
static OSType Nlm_GetOSType (Nlm_CharPtr str, OSType dfault)

{
  OSType  rsult;

  rsult = dfault;
  if (str != NULL && str [0] != '\0') {
    rsult = *(OSType*) str;
  }
  return rsult;
}
#endif

NLM_EXTERN void LIBCALL Nlm_FileCreate (Nlm_CharPtr fileName, Nlm_CharPtr type, Nlm_CharPtr creator)

{
  FILE      *fp;
#ifdef WIN_MAC
  OSType    fCreator;
  Nlm_Int2  fError;
  OSType    fType;
  Nlm_Char  temp [256];
#endif

  if (fileName != NULL && fileName [0] != '\0') {
#ifdef WIN_MAC
    if (type != NULL || creator != NULL) {
      fp = Nlm_FileOpen (fileName, "r");
      if (fp == NULL) {
        fType = Nlm_GetOSType (type, 'TEXT');
        fCreator = Nlm_GetOSType (creator, '    ');
        Nlm_StringNCpy_0(temp, fileName, sizeof(temp));
        Nlm_CtoPstr ((Nlm_CharPtr) temp);
        fError = Create ((StringPtr) temp, 0, fCreator, fType);
      } else {
        Nlm_FileClose (fp);
      }
    }
#else
    fp = Nlm_FileOpen (fileName, "w");
    if (fp != NULL) {
      Nlm_FileClose (fp);
    }
#endif
  }
}

/*****************************************************************************
*
*   CreateDir(pathname)
*
*****************************************************************************/

NLM_EXTERN Nlm_Boolean LIBCALL  Nlm_CreateDir (Nlm_CharPtr pathname)
{
#ifndef OS_VMS
  size_t          len;
  Nlm_Char        path[PATH_MAX];
#endif
#ifdef OS_MAC
  long            dirID;
  Nlm_Char        dirname [FILENAME_MAX];
  OSErr           err;
  HParamBlockRec  pbh;
  Nlm_CharPtr     ptr;
  short           vRefNum;
#endif
#ifdef OS_UNIX
  mode_t          oldmask;
#endif
  Nlm_Boolean     rsult = FALSE;

  if (pathname != NULL && pathname [0] != '\0') {
#ifdef OS_MAC
    Nlm_StringNCpy_0(path, pathname, sizeof(path));
    len = Nlm_StringLen (path);
    if (len > 0 && path [len - 1] == DIRDELIMCHR) {
        path [len - 1] = '\0';
    }
    memset (&pbh, 0, sizeof (HParamBlockRec));
    vRefNum = Nlm_MacGetVRefNum (path, &err);
    if (err == noErr) {
      dirID = Nlm_MacGetDirID (path, vRefNum, &err);
      if (err == noErr) {
        ptr = Nlm_StringRChr (path, (int) DIRDELIMCHR);
        if (ptr != NULL) {
          ptr++;
          Nlm_StringNCpy_0(dirname, ptr, sizeof(dirname));
          Nlm_CtoPstr ((Nlm_CharPtr) dirname);
          pbh.fileParam.ioNamePtr = (StringPtr) dirname;
          pbh.fileParam.ioVRefNum = vRefNum;
          pbh.fileParam.ioDirID = dirID;
          err = PBDirCreate (&pbh, FALSE);
          rsult = (Nlm_Boolean) (err == noErr || err == dupFNErr);
        }
      }
    }
#endif
#if defined(OS_MSWIN) || defined(OS_NT)
    Nlm_StringNCpy_0(path, pathname, sizeof(path));
    len = Nlm_StringLen (path);
    if (len > 0 && path [len - 1] == DIRDELIMCHR) {
        path [len - 1] = '\0';
    }
    rsult = (Nlm_Boolean) (mkdir ((char *) path) == 0);
    if (errno == EACCES) { /* it's O.K. if it was already there */
	rsult = TRUE;
    }
#endif
#ifdef OS_UNIX
    oldmask = umask (0000);
    Nlm_StringNCpy_0(path, pathname, sizeof(path));
    len = Nlm_StringLen (path);
    if (len > 0 && path [len - 1] == DIRDELIMCHR) {
        path [len - 1] = '\0';
    }
    rsult = (Nlm_Boolean) (mkdir ((char *) path, 0755) == 0);
    if (errno == EEXIST) { /* it's O.K. if it was already there */
	rsult = TRUE;
    }
    umask (oldmask);
#endif
#ifdef OS_VMS
    rsult = (Nlm_Boolean) (mkdir ((char *) pathname, 0755) == 0);
#endif
  }
  return rsult;
}

/*****************************************************************************
*
*   DirectoryContents()
*
*****************************************************************************/

NLM_EXTERN ValNodePtr LIBCALL Nlm_DirCatalog (Nlm_CharPtr pathname)

{
#ifdef OS_MAC
  long            dirID;
  OSErr           err;
  short           index;
  unsigned short  num;
  Nlm_Char        path[PATH_MAX];
  CInfoPBRec      pbc;
  HParamBlockRec  pbh;
  short           vRefNum;
#endif
#ifdef OS_UNIX
  Nlm_Char        buf [256];
  Nlm_Char        ch;
  Nlm_Uint1       choice;
  Nlm_Char        cmmd [PATH_MAX + 20];
  FILE            *fp;
  Nlm_CharPtr     ptr;
#endif
  ValNodePtr      vnp = NULL;

  if (pathname != NULL && pathname [0] != '\0') {
#ifdef OS_MAC
    Nlm_StringNCpy_0 (path, pathname, sizeof (path));
    Nlm_CtoPstr ((Nlm_CharPtr) path);
    Nlm_MemSet ((Nlm_VoidPtr) (&pbh), 0, sizeof (HParamBlockRec));
    pbh.volumeParam.ioNamePtr = (StringPtr) path;
    pbh.volumeParam.ioVolIndex = -1;
    err = PBHGetVInfo (&pbh, FALSE);
    if (err != noErr) return NULL;
    vRefNum = pbh.volumeParam.ioVRefNum;
    Nlm_StringNCpy_0 (path, pathname, sizeof (path));
    Nlm_CtoPstr ((Nlm_CharPtr) path);
    Nlm_MemSet ((Nlm_VoidPtr) (&pbc), 0, sizeof (CInfoPBRec));
    pbc.dirInfo.ioNamePtr = (StringPtr) path;
    pbc.dirInfo.ioVRefNum = vRefNum;
    err = PBGetCatInfo (&pbc, FALSE);
    if (err != noErr) return NULL;
    if (pbc.dirInfo.ioFlAttrib & 16) {
      num = pbc.dirInfo.ioDrNmFls;
      dirID = pbc.dirInfo.ioDrDirID;
      for (index = 1; index <= num; index++) {
        Nlm_MemSet ((Nlm_VoidPtr) (&pbc), 0, sizeof (CInfoPBRec));
        pbc.dirInfo.ioNamePtr = (StringPtr) path;
        pbc.dirInfo.ioVRefNum = vRefNum;
        pbc.dirInfo.ioFDirIndex = index;
        pbc.dirInfo.ioDrDirID = dirID;
        pbc.dirInfo.ioACUser = 0;
        err = PBGetCatInfo (&pbc, FALSE);
        if (err == noErr) {
          Nlm_PtoCstr ((Nlm_CharPtr) path);
          if (pbc.dirInfo.ioFlAttrib & 16) {
            ValNodeCopyStr (&vnp, 1, path);
          } else {
            ValNodeCopyStr (&vnp, 0, path);
          }
        }
      }
    }
#endif
#if defined(WIN32)
    {{
      Nlm_Char x_path[PATH_MAX];
      WIN32_FIND_DATA fData;
      HANDLE hFindFile;
      Nlm_StringNCpy_0(x_path, pathname, sizeof(x_path) - 5);
      Nlm_StringCat(x_path, "\\*.*");
      hFindFile = FindFirstFile(x_path, &fData);
      if (hFindFile == INVALID_HANDLE_VALUE)
        return 0;
      do {
        if (fData.cFileName[0] != '.'  ||
            (fData.cFileName[1] != '.'  &&  fData.cFileName[1] != '\0'))
          ValNodeCopyStr
            (&vnp, (fData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) ? 1 : 0,
             fData.cFileName);
      } while ( FindNextFile(hFindFile, &fData) );
      FindClose(hFindFile);
    }}
#endif
#if defined(WIN16)
    {{
      Nlm_Char x_path[PATH_MAX];
      struct find_t fData;
      Nlm_StringNCpy_0(x_path, pathname, sizeof(x_path) - 5);
      Nlm_StringCat(x_path, "\\*.*");
      if (_dos_findfirst(x_path, _A_SUBDIR|_A_RDONLY|_A_NORMAL|_A_ARCH, &fData)
          != 0)
        return 0;
      do {
        if (fData.name[0] != '.'  ||
            (fData.name[1] != '.'  &&  fData.name[1] != '\0'))
          ValNodeCopyStr(&vnp, (fData.attrib & _A_SUBDIR) ? 1 : 0, fData.name);
      } while (_dos_findnext(&fData) == 0);
    }}
#endif
#ifdef OS_UNIX
    sprintf (cmmd, "ls -1p %s 2>/dev/null", pathname);
    fp = popen (cmmd, "r");
    if (fp == NULL) return NULL;
    while (Nlm_FileGets (buf, sizeof (buf), fp) != NULL) {
      ptr = buf;
      ch = *ptr;
      while (ch != '\0' && ch != '\n' && ch != '\r') {
        ptr++;
        ch = *ptr;
      }
      *ptr = '\0';
      choice = 0;
      ptr = Nlm_StringChr (buf, '/');
      if (ptr != NULL) {
        *ptr = '\0';
        choice = 1;
      }
      ValNodeCopyStr (&vnp, choice, buf);
    }
    pclose (fp);
#endif
#ifdef OS_VMS
#endif
  }
  return vnp;
}

/*****************************************************************************
*
*   TmpNam()
*
*****************************************************************************/
NLM_EXTERN Nlm_CharPtr LIBCALL Nlm_TmpNam (Nlm_CharPtr s)

{
#ifdef TEMPNAM_AVAIL
    char *filename;
    static Nlm_Char save_filename[L_tmpnam+30];

    /* emulate tmpnam(), except get the benefits of tempnam()'s ability to */
    /* place the files in another directory specified by the environment   */
    /* variable TMPDIR                                                     */

    filename = tempnam(NULL, NULL);

    if (s == NULL)
    { /* return pointer to static string */
        if (filename != NULL) {
          strcpy ((char *) save_filename, (char *) filename);
          free ((void *) filename);
        } else {
          save_filename [0] = '\0';
        }
        return save_filename;
    } else {
        if (filename != NULL) {
          strcpy ((char *) save_filename, (char *) filename);
          Nlm_StrCpy (s, save_filename);
          free ((void *) filename);
        } else {
          *s = '\0';
        }
        return s;
    }
#else
#ifdef OS_MAC
    static Nlm_Char  directory [PATH_MAX];
    OSErr        err;
    long         gesResponse;
    long         newDirID;
    short        newVRefNum;
    long         oldDirID;
    short        oldVRefNum;
    CInfoPBRec   params;
    Nlm_Char     temp [PATH_MAX];
    Nlm_CharPtr  tmp;
    Nlm_Boolean  useTempFolder;
    char * filename;

    useTempFolder = FALSE;
    if (! Gestalt (gestaltFindFolderAttr, &gesResponse) &&
        (gesResponse & (1 << gestaltFindFolderPresent))) {
      err = FindFolder(kOnSystemDisk, kTemporaryFolderType,
                       kCreateFolder, &newVRefNum, &newDirID);
      if (err == noErr) {
        useTempFolder = TRUE;
        err = Nlm_SetDefault (newVRefNum, newDirID, &oldVRefNum, &oldDirID);
      }
    }
    filename = tmpnam (NULL);
    if (useTempFolder) {
      err = Nlm_RestoreDefault (oldVRefNum, oldDirID);
      temp [0] = '\0';
      params.dirInfo.ioNamePtr = (StringPtr) directory;
      params.dirInfo.ioDrParID = newDirID;
      do {
        params.dirInfo.ioVRefNum = newVRefNum;
        params.dirInfo.ioFDirIndex = -1;
        params.dirInfo.ioDrDirID = params.dirInfo.ioDrParID;
        err = PBGetCatInfo (&params, FALSE);
        Nlm_PtoCstr ((Nlm_CharPtr) directory);
        Nlm_StringCat (directory, DIRDELIMSTR);
        Nlm_StringCat (directory, temp);
        Nlm_StringCpy (temp, directory);
      } while (params.dirInfo.ioDrDirID != fsRtDirID);
      tmp = Nlm_StringMove (directory, temp);
      tmp = Nlm_StringMove (tmp, (Nlm_CharPtr) filename);
      if (s == NULL) {
          return (Nlm_CharPtr) directory;
      } else {
          s [0] = '\0';
          Nlm_StringCpy (s, directory);
          return s;
      }
    } else {
      if (s == NULL) {
          return (Nlm_CharPtr) filename;
      } else {
          s [0] = '\0';
          Nlm_StringCpy (s, filename);
          return s;
      }
    }
#else
    char * filename;

    filename = tmpnam (NULL);
    if (s == NULL) {
        return (Nlm_CharPtr) filename;
    } else {
        s [0] = '\0';
        Nlm_StringCpy (s, filename);
        return s;
    }
#endif
#endif
}

/*****************************************************************************
*
*   CD-ROM Ejection Routines
*
*****************************************************************************/

NLM_EXTERN Nlm_Boolean LIBCALL  Nlm_EjectCd(Nlm_CharPtr sVolume, Nlm_CharPtr deviceName,
			Nlm_CharPtr rawDeviceName, 
			Nlm_CharPtr mountPoint,
			Nlm_CharPtr mountCmd)
{
    Nlm_Boolean retval = FALSE;
#ifdef OS_MAC
    OSErr err;
    Nlm_CharPtr prob_area = "Ejection";
    Nlm_Char    temp [64];

    
    Nlm_StringNCpy_0(temp, sVolume, sizeof(temp) - 1);
    Nlm_StringCat (temp, ":");
    if ((err = Eject((StringPtr) NULL, Nlm_MacGetVRefNum(temp, NULL))) == noErr)
    {
        if ((err = UnmountVol((StringPtr) NULL, Nlm_MacGetVRefNum(temp, NULL))) == noErr)
        	return TRUE;
        
        /* We should still return TRUE if we Eject() successfully but failed to    */
        /* unmount the volume; however, we need to warn them, because a subsequent */
        /* GetFInfo() will result in a bus error, at least with System 7.0.        */
        retval = TRUE;
        prob_area = "Unmounting";
    }

    switch (err) {
    case bdNamErr:
	ErrPostEx(SEV_ERROR,E_File,E_CdEject,"%s error - bad volume name %s", prob_area, sVolume);
	break;
    case extFSErr:
	ErrPostEx(SEV_ERROR,E_File,E_CdEject,"%s error - external file system %s", prob_area, sVolume);
	break;
    case ioErr:
	ErrPostEx(SEV_ERROR,E_File,E_CdEject,"%s error - I/O error %s", prob_area, sVolume);
	break;
    case nsDrvErr:
	ErrPostEx(SEV_ERROR,E_File,E_CdEject,"%s error - No such drive %s", prob_area, sVolume);
	break;
    case nsvErr:
	ErrPostEx(SEV_ERROR,E_File,E_CdEject,"%s error - No such volume %s", prob_area, sVolume);
	break;
    case paramErr:
	ErrPostEx(SEV_ERROR,E_File,E_CdEject,"%s error - No default volume %s", prob_area, sVolume);
	break;
    }
    
    return retval;
#endif /* OS_MAC */

#ifdef OS_UNIX
	char cmd[100];
#endif
#ifdef OS_UNIX_SUN
	int fd;

	if (deviceName == NULL)
	{
		deviceName = DEFAULT_CDROM;
	}

	if (rawDeviceName == NULL)
	{
		rawDeviceName = DEFAULT_RAW_CDROM;
	}

	/* Open the CD-ROM character-based device */
	if ((fd = open(rawDeviceName, O_RDONLY, 0)) < 0)
	{
		ErrPostEx(SEV_ERROR,E_File,E_CdEject,"Ejection error - Unable to open device %s", rawDeviceName);
		return FALSE;
	}

	retval = ioctl(fd, CDROMEJECT, 0) >= 0;
	close (fd);

	if (! retval)
	{
		ErrPostEx(SEV_ERROR,E_File,E_CdEject,"Ejection error - Ioctl failure for %s", rawDeviceName);
    	return FALSE;
	}
#endif /* OS_UNIX_SUN */

#ifdef OS_UNIX
	/* Now try to unmount device using (un)mount-script */
	if (mountCmd != NULL)
	{
		sprintf(cmd, "%s -u %s >/dev/null 2>/dev/null", mountCmd,
				deviceName);
		retval = system(cmd) == 0;
	}
	else {
		if (deviceName != NULL)
		{
			retval = Message(MSG_OKC,
				            "Unmount device <%s> now; select OK when completed",
						    deviceName) != ANS_CANCEL;
		}
		else
		if (sVolume != NULL)
		{
			retval = Message(MSG_OKC,
						    "Unmount volume <%s> now; select OK when completed",
						    sVolume) != ANS_CANCEL;
		}
		else
		{
			retval = Message(MSG_OKC,
							"Unmount CD-ROM now; select OK when completed") !=
							ANS_CANCEL;
		}
	}
#endif /* OS_UNIX */

#ifdef OS_VMS
	char  cmd[100];
	char  tmp[100];
	char* cPtr;


	if ( mountPoint == NULL || *mountPoint == '\0' ) 
		strcpy(tmp,DEFAULT_CDROM); 
	else {
		strcpy(tmp,mountPoint);
		if ( cPtr = strchr(tmp,':') ) *(cPtr+1) = '\0';
	}
	/* 
	** Try to mount device using mount-script 
	*/

	sprintf(cmd, "CD_DISMOUNT/UNLOAD %s",tmp);
	retval = (system(cmd) == 0);

	 Message(MSG_OK,
		"Press the eject button on <%s>.",tmp);

#endif
	
    return retval;
}

NLM_EXTERN Nlm_Boolean LIBCALL  Nlm_MountCd(Nlm_CharPtr sVolume, Nlm_CharPtr deviceName,
			Nlm_CharPtr mountPoint, Nlm_CharPtr mountCmd)
{
	Nlm_Boolean retval = FALSE;

#ifdef OS_UNIX
	char cmd[100];

	if (deviceName == NULL)
	{
		deviceName = DEFAULT_CDROM;
	}

	/* Try to mount device using mount-script */
	if (mountCmd != NULL)
	{
		sprintf(cmd, "%s -m %s %s >/dev/null 2>/dev/null", mountCmd, deviceName,
				mountPoint != NULL ? mountPoint : "");
		retval = system(cmd) == 0;
	}
	else {
		if (deviceName != NULL)
		{
		}
		else
		{
			retval = Message(MSG_OKC,
							"Mount CD-ROM now; select OK when completed") !=
							ANS_CANCEL;
		}
	}
#endif

#ifdef OS_VMS
	char  cmd[100];
	char  tmp[100];
	char* cPtr;


	if ( mountPoint == NULL || *mountPoint == '\0' ) 
		strcpy(tmp,DEFAULT_CDROM); 
	else {
		strcpy(tmp,mountPoint);
		if ( cPtr = strchr(tmp,':') ) *(cPtr+1) = '\0';
	}


	/* Try to mount device using mount-script */

	sprintf(cmd, "CD_MOUNT/MEDIA=CDROM/OVERRIDE=IDENTIFICATION/NOASSIST %s",
          tmp);

	retval = (system(cmd) == 0);

#endif

	return retval;
}


