/*   ncbienv.c
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
* File Name:  ncbienv.c
*
* Author:  Ostell
*
* Version Creation Date:   7/7/91
*
* $Revision: 6.6 $
*
* File Description:
*       portable environment functions, companions for ncbimain.c
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: ncbienv.c,v $
* Revision 6.6  1998/12/10 17:04:06  vakatov
* Fixed to compile under LINUX(Red Hat 2.XX, gcc, with POSIX threads)
*
* Revision 6.5  1997/11/26 21:26:10  vakatov
* Fixed errors and warnings issued by C and C++ (GNU and Sun) compilers
*
* Revision 6.4  1997/10/29 02:43:31  vakatov
* Type castings to pass through the C++ compiler
*
* Revision 6.3  1997/10/27 21:58:11  vakatov
* Added Nlm_FreeArgs() to reset args earlier set by GetArgs[Silent]()
*
* Revision 6.2  1997/09/09 00:05:43  vakatov
* Nlm_x_HasConsole() -- "fileno" instead of "_fileno" to pass Borland compiler
*
* Revision 6.1  1997/08/27 16:16:38  vakatov
* [WIN32] Nlm_x_HasConsole() -- fixed for the case when the
* Vibrant-based(GUI) application is run from the console(command prompt)
*
* Revision 5.18  1997/07/22 19:11:18  vakatov
* Separated Main() from GetArg[svc]() functions;  [WIN_MSWIN] converged
* console and GUI libraries; [for WIN32-DLL] encapsulated global variables
*
* Revision 5.17  1997/06/02 15:28:36  vakatov
* [WIN32]  Read/write config(*.ini) files via WinSDK calls([WIN_MSWIN]-like)
*
* Revision 5.16  1997/05/08 19:37:08  vakatov
* [WIN_MSWIN]  Nlm_SetupArguments_ST() -- interpret all text(including
* space symbols) embraced by a pair of ""(or '') quote marks as a
* single command argument;  strip the embracing quote marks
*
* Revision 5.15  1997/03/05 20:32:55  vakatov
* Nlm_WorkGetAppParam: now free Nlm_lastParamFile before reset, avoid mem.leak
*
 * Revision 5.14  1997/03/03  16:44:00  kans
 * cleanup of memory leaks accidentally freed the cached file name,
 * so now only the Nlm_envList is freed with Nlm_FreeEnvData
 *
 * Revision 5.13  1997/02/27  19:33:48  vakatov
 * Nlm_ReadConfigFile():  free the obsolete config data --> fight memory leaks
 *
 * Revision 5.12  1997/01/28  21:19:12  kans
 * <Desk.h>, <OSEvents.h> and <GestaltEqu.h> are obsolete in CodeWarrior
 *
 * Revision 5.11  1997/01/24  17:03:49  epstein
 * make threaded version compatible with OSF/1
 *
 * Revision 5.10  1997/01/03  16:12:07  vakatov
 * Fixed inaccurate string copying -- <mostly potential> 1-byte exceeding of
 * the string size by StringNCat;  missing terminating '\0' by StringNCpy.
 *
 * Revision 5.9  1996/12/30  15:13:12  vakatov
 * [WIN_MSWIN]  Command-line parsing implemented inside Nlm_SetupArguments()
 *
 * Revision 5.8  1996/12/16  22:38:53  vakatov
 * Rolled back the changes made in "* Revision 5.4.  1996/11/27  20:38:14
 * epstein"
 *
 * Revision 5.7  1996/12/04  21:44:47  vakatov
 * [OS_UNIX][POSIX_THREADS_AVAIL]  Added _POSIX1C case (see Rev.5.5)
 *
 * Revision 5.6  1996/12/03  21:48:33  vakatov
 * Adopted for 32-bit MS-Windows DLLs
 *
 * Revision 5.5  1996/11/27  21:55:17  vakatov
 * [OS_UNIX][POSIX_THREADS_AVAIL]  Added (_POSIX_C_SOURCE - 0 >= 199506L)
 * preprocessor condition to match POSIX.1c function interface.
 *
 * Revision 5.4  1996/11/27  20:38:14  epstein
 * disable error logging when opening files in parameter fetching functions
 *
 * Revision 5.3  1996/11/25  19:04:26  vakatov
 * Wrapped all basic functions to MT-safe wrappers(named by adding '_ST'
 * suffix to the original function names).
 *
 * Revision 5.2  1996/08/19  18:46:06  vakatov
 * [WIN32]  Made modifications to let one create console applications
 *
* Revision 5.1  1996/07/16 19:56:00  vakatov
* Just castings
*
 * Revision 5.0  1996/05/28  13:18:57  ostell
 * Set to revision 5.0
 *
 * Revision 4.9  1996/02/15  22:00:49  kans
 * changed platform symbol back to OS_NT
 *
 * Revision 4.8  1996/01/31  20:29:57  epstein
 * fix uninitialized ptr in Nlm_GetHome()
 *
 * Revision 4.7  1996/01/29  22:33:37  epstein
 * Added GetHome() changes per Mr. G.M. Lack (gml4410@ggr.co.uk)
 *
 * Revision 4.6  1995/12/13  17:18:51  kans
 * fixed caching bug (JE & JK)
 *
 * Revision 4.5  1995/10/28  15:03:20  ostell
 * added casts to quiet DOS compile warnings
 *
 * Revision 4.4  1995/10/11  13:53:05  kans
 * made some variables static
 *
 * Revision 4.3  1995/10/06  19:59:18  epstein
 * more performance fixes
 *
 * Revision 4.2  1995/10/06  15:53:31  epstein
 * add CacheAppParam() and FlushAppParam() for improved performance
 *
 * Revision 4.1  1995/10/03  15:59:21  epstein
 * add NCBI_DONT_USE_LOCAL_CONFIG environment variable to avoid local config files
 *
*  7/7/91   Kans        Multiple configuration files, get and set functions
*  9-20-91  Schuler     GetAppParam takes default value as an argument
* 01-14-94  Epstein     Merged ncbienv.{unx,vms,dos,msw,mac} into a single file
* 06-14-94  Schuler     Put SetAppPropery("ProgramName",..) in SetupArguments
* 06-14-94  Schuler     Add LIBCALL to SetupArguments
* 08-23-94  Schuler     Add SetupArguments case for OS_NT/WIN_DUMB
* 01-29-96  Epstein     Added GetHome() changes per Mr. G.M. Lack
*
* ==========================================================================
*/

#ifdef OS_MAC
#ifdef PROC_MC680X0
#define OBSOLETE
#endif
#endif

#include "corepriv.h"

#ifdef OS_UNIX
#include <pwd.h>
#endif /* OS_UNIX */
#ifdef OS_MAC
#include <Gestalt.h>
#include <Folders.h>
#include <Processes.h>
#endif /* OS_MAC */


typedef struct nlm_env_item {
  struct nlm_env_item  PNTR next;
  Nlm_CharPtr          name;
  Nlm_CharPtr          comment;
  Nlm_CharPtr          value;
} Nlm_env_item, PNTR Nlm_env_itemPtr;

typedef struct nlm_env_sect {
  struct nlm_env_sect  PNTR next;
  Nlm_CharPtr          name;
  Nlm_CharPtr          comment;
  Nlm_Boolean          transientOnly; /* this field used only by Transient fns */
  struct nlm_env_item  PNTR children;
} Nlm_env_sect, PNTR Nlm_env_sectPtr;

typedef struct nlm_env_file {
  struct nlm_env_file  PNTR next;
  Nlm_CharPtr          name;
  Nlm_env_sectPtr      envList;
} Nlm_env_file, PNTR Nlm_env_filePtr;


static Nlm_env_filePtr Nlm_transientFileList = NULL;
static Nlm_Boolean caching = FALSE;
static Nlm_Boolean dirty = FALSE;
static Nlm_Boolean mustCreateFile = TRUE;

static Nlm_Boolean Nlm_Qualified PROTO((Nlm_CharPtr path));
static Nlm_Boolean Nlm_TransientLookup PROTO((Nlm_CharPtr file, Nlm_CharPtr section, Nlm_CharPtr type, Nlm_CharPtr dflt, Nlm_CharPtr buf, Nlm_Int2 buflen));
static void Nlm_TransientLogSetApp PROTO((Nlm_CharPtr file, Nlm_CharPtr section, Nlm_CharPtr type, Nlm_CharPtr value));
static void Nlm_FreeEnvData PROTO((Nlm_env_sectPtr esp));
static void Nlm_FreeTransientData PROTO((void));

#ifndef OS_MSWIN
static FILE *Nlm_OpenConfigFile PROTO((Nlm_CharPtr file, Nlm_Boolean writeMode, Nlm_Boolean create));
static Nlm_CharPtr Nlm_TrimString PROTO((Nlm_CharPtr str));
static Nlm_Boolean Nlm_ReadConfigFile PROTO((FILE *fp));
static Nlm_env_sectPtr Nlm_FindConfigSection PROTO((Nlm_CharPtr section));
static Nlm_env_itemPtr Nlm_FindConfigItem PROTO((Nlm_CharPtr section, Nlm_CharPtr type, Nlm_Boolean create));
static Nlm_Boolean Nlm_WriteConfigFile PROTO((FILE *fp));
static void Nlm_PutComment PROTO((Nlm_CharPtr s, FILE *fp));
static void Nlm_FreeConfigFileData PROTO((void));
static void Nlm_FlushConfigFile PROTO((Nlm_CharPtr file, Nlm_Boolean create));

static Nlm_env_sectPtr Nlm_envList = NULL;
static Nlm_CharPtr Nlm_lastParamFile = NULL;
static Nlm_CharPtr Nlm_bottomComment = NULL;

/* always FALSE, because this file is trying to emulating MS Windows's  */
/* handling of comments in Param files; however, just change this value */
/* to TRUE to turn this approach around                                 */
static Nlm_Boolean destroyDeadComments = FALSE;

/*****************************************************************************
*
* The "guts" of:
*   Nlm_GetAppParam (file, section, type, dflt, buf, buflen)
*      finds parameters from configuration files
*      this version, searching for configuration file(s) in a
*      platform-dependent basis as handled by Nlm_OpenConfigFile()
*
*      if configuration file is found, tries to read the parameter from it.
*
*****************************************************************************/

static Nlm_Int2 Nlm_WorkGetAppParam (Nlm_CharPtr file, Nlm_CharPtr section, Nlm_CharPtr type, Nlm_CharPtr dflt, Nlm_CharPtr buf, Nlm_Int2 buflen, Nlm_Boolean searchTransient)
{
  Nlm_env_itemPtr  eip;
  FILE             *fp;

  if (buf == NULL  ||  buflen <= 0)
    return 0;

  *buf = '\0';
  if (searchTransient  &&
      Nlm_TransientLookup(file, section, type, dflt, buf, buflen))
    {
      return (Nlm_Int2)Nlm_StringLen(buf);
    }

  if ( dflt )
    Nlm_StringNCpy_0(buf, dflt, buflen);

  if (file    == NULL   ||  *file    == '\0'  ||
      section == NULL   ||  *section == '\0')
    return (Nlm_Int2)Nlm_StringLen( buf );

  if (Nlm_lastParamFile == NULL  ||
      Nlm_StringICmp(Nlm_lastParamFile, file) != 0)
    {
      mustCreateFile = TRUE;
      if ( caching )
        Nlm_FlushAppParam();
      Nlm_FreeConfigFileData();
      fp = Nlm_OpenConfigFile(file, FALSE, FALSE);
      if (fp != NULL)
        {
          MemFree( Nlm_lastParamFile );
          Nlm_lastParamFile = Nlm_StringSave( file );
          Nlm_ReadConfigFile( fp );
          Nlm_FileClose( fp );
        }
    }

  if (type != NULL  &&  *type != '\0')
    {
      eip = Nlm_FindConfigItem(section, type, FALSE);
      if (eip != NULL)
        Nlm_StringNCpy_0(buf, eip->value, buflen);
    }
  else
    { /* return all the types in that section */
      Nlm_env_sectPtr  esp    = Nlm_FindConfigSection( section );
      Nlm_Int2         totlen = 0;
      *buf = '\0';
      if (esp == NULL)
        return 0;

      /* traverse the children, allowing the null chars to be inserted */
      /* in between each type-name                                     */
      for (eip = esp->children;
           eip != NULL  &&  totlen != buflen;  eip = eip->next)
        {
          Nlm_Int2 bytesToAppend = StrLen(eip->name) + 1;
          bytesToAppend = MIN(bytesToAppend, buflen - totlen);
          StrNCpy(&buf[totlen], eip->name, bytesToAppend - 1);
          totlen += bytesToAppend;
          buf[totlen - 1] = '\0';
        }
    }

  return (Nlm_Int2)Nlm_StringLen(buf);
}


/*****************************************************************************
*
*   Nlm_SetAppParam (file, section, type, value)
*      finds paths for types of data and fills in path in buf
*      this version
*      1)  looks in the current directory for ".filerc", but will not
*          create a new file in this directory.
*      2)  then looks in the home directory for ".filerc".
*      3)  then looks for an environment variable "NCBI" and takes its
*          value as a complete path to a directory containing the
*          configuration file ".filerc".
*      if configuration file is found, tries to write the parameter to it.
*
*****************************************************************************/

static Nlm_Boolean Nlm_SetAppParam_ST PROTO ((Nlm_CharPtr file, Nlm_CharPtr section, Nlm_CharPtr type, Nlm_CharPtr value));
static Nlm_Boolean Nlm_SetAppParam_ST (Nlm_CharPtr file, Nlm_CharPtr section, Nlm_CharPtr type, Nlm_CharPtr value)
{
  Nlm_env_itemPtr  eip;
  Nlm_env_sectPtr  esp;
  FILE             *fp = NULL;
  Nlm_Boolean      rsult;

  rsult = FALSE;
  if (file != NULL && *file != '\0' && section != NULL && *section != '\0') {
    Nlm_TransientLogSetApp (file, section, type, value);

    if (Nlm_lastParamFile == NULL || Nlm_StringICmp(Nlm_lastParamFile, file) != 0) {
      mustCreateFile = TRUE;
    }
    if (mustCreateFile)
      fp = Nlm_OpenConfigFile (file, FALSE, TRUE);

    if (TRUE) {
      if (fp != NULL) {
        Nlm_FlushAppParam();
        Nlm_FreeConfigFileData();
        Nlm_lastParamFile = Nlm_StringSave(file);
        Nlm_ReadConfigFile (fp);
        Nlm_FileClose (fp);
        mustCreateFile = FALSE;
      }
      if (type != NULL && *type != '\0')
      {
        eip = Nlm_FindConfigItem (section, type, TRUE);
        if (eip != NULL) {
          if (eip->value != NULL) {
            eip->value = (Nlm_CharPtr) Nlm_MemFree (eip->value);
          }
          eip->value = Nlm_StringSave (value);
          rsult = TRUE;
        }
      }
      else { /* wipe out that section */
        esp = Nlm_FindConfigSection (section);
        if (esp != NULL) { /* kill section by deleting name (leave comments)*/
          esp->name = (Nlm_CharPtr) Nlm_MemFree(esp->name);
          rsult = TRUE;
        }
      }

      if (rsult) {
        dirty = TRUE;
      }
      if (! caching) {
        Nlm_FlushConfigFile(file, TRUE);
      }
    }
  }

  return rsult;
}

/*****************************************************************************
*
*   Nlm_FlushAppParam()
*      flush the current parameter file's parameters to disk
*
*****************************************************************************/
static void Nlm_FlushAppParam_ST PROTO ((void));
static void Nlm_FlushAppParam_ST(void)
{
  if (Nlm_lastParamFile != NULL)
    Nlm_FlushConfigFile(Nlm_lastParamFile, TRUE);
}

/*****************************************************************************
*
*   Nlm_CacheAppParam()
*      Indicates whether data should be flushed to disk after each call
*      to SetAppParam()
*
*****************************************************************************/
static Nlm_Boolean Nlm_CacheAppParam_ST PROTO ((Nlm_Boolean value));
static Nlm_Boolean Nlm_CacheAppParam_ST(Nlm_Boolean value)
{
  Nlm_Boolean oldvalue = caching;

  caching = value;
  if (! value)
    Nlm_FlushAppParam();

  return oldvalue;
}


#ifdef OS_UNIX
#define NLM_POSIX1B (_POSIX1B || _POSIX1C || \
                     (_POSIX_C_SOURCE - 0 >= 199309L) || \
                     defined(_POSIX_PTHREAD_SEMANTICS))

/*****************************************************************************
*
*   Nlm_GetHome (buf, buflen)
*      returns the path of the home directory
*
*****************************************************************************/

static Nlm_Boolean Nlm_GetHome PROTO((Nlm_CharPtr buf, Nlm_Int2 buflen));
static Nlm_Boolean Nlm_GetHome (Nlm_CharPtr buf, Nlm_Int2 buflen)
{
  static Nlm_Char saveHome[PATH_MAX + 1];

  /* Have we passed this way before?  If not, then try get the info. */
  if (saveHome[0] == '\0'  &&  saveHome[1] == '\0')
    {
#if  defined(SOLARIS_THREADS_AVAIL)  ||  defined(POSIX_THREADS_AVAIL)
    static TNlmMutex saveHomeMutex;
    VERIFY_HARD ( NlmMutexLockEx( &saveHomeMutex ) == 0 );
    if (saveHome[0] == '\0'  &&  saveHome[1] == '\0')
    {
#endif

      struct passwd *pwd_ptr = NULL;

      /* Get the user login name */
#if  defined(SOLARIS_THREADS_AVAIL)  ||  defined(POSIX_THREADS_AVAIL)
      struct passwd pwd;
#ifndef LOGNAME_MAX
#if defined(MAXLOGNAME)
#define LOGNAME_MAX MAXLOGNAME
#elif defined(_POSIX_LOGIN_NAME_MAX)
#define LOGNAME_MAX _POSIX_LOGIN_NAME_MAX
#endif
#endif /* ndef LOGNAME_MAX */
      Nlm_Char      login_name[LOGNAME_MAX + 1];
      Nlm_Char      pwd_buffer[LOGNAME_MAX + PATH_MAX + 1024 + 1];
      Nlm_Boolean   ok = getlogin_r(login_name, sizeof(login_name)) ?
#if NLM_POSIX1B
                        FALSE : TRUE;
#else
                        TRUE : FALSE;
#endif
#else

      Nlm_CharPtr login_name = getlogin();
      Nlm_Boolean ok         = (login_name != NULL);
#endif

      /* Get the info using user-login-name */
      if ( ok )
        {
#if  defined(SOLARIS_THREADS_AVAIL)  ||  defined(POSIX_THREADS_AVAIL)
          pwd_ptr = &pwd;
          ok = (getpwnam_r(login_name, pwd_ptr, pwd_buffer, sizeof(pwd_buffer)
#if NLM_POSIX1B
                           , &pwd_ptr) == 0);
#else
                           ) != NULL);
#endif
#else
          pwd_ptr = getpwnam( login_name );
#endif
        }


      /* Get the info using user-ID */
      if ( !ok )
        {
#if  defined(SOLARIS_THREADS_AVAIL)  ||  defined(POSIX_THREADS_AVAIL)
          pwd_ptr = &pwd;
          ok = (getpwuid_r(getuid(), pwd_ptr, pwd_buffer, sizeof(pwd_buffer)
#if NLM_POSIX1B
                           , &pwd_ptr) == 0);
#else
                           ) != NULL);
#endif
#else
          pwd_ptr = getpwuid( getuid() );
#endif
        }

      /* Now we either *do* have it or we *never* will. */
      /* Remember failure by setting first two bytes to "\000\001" */
      saveHome[1] = (Nlm_Char)1;
      if ( ok )
        Nlm_StringNCpy_0(saveHome, pwd_ptr->pw_dir, sizeof(saveHome));

#if  defined(SOLARIS_THREADS_AVAIL)  ||  defined(POSIX_THREADS_AVAIL)
    VERIFY_HARD ( NlmMutexUnlock(saveHomeMutex) == 0 );
    }
#endif
    }

  /* Return the now-saved value and success code */
  Nlm_StringNCpy_0(buf, saveHome, buflen);
  return  (Nlm_Boolean)(*saveHome != '\0');
}


/*****************************************************************************
*
*   Nlm_OpenConfigFile (file, writeMode, create)
*      returns a file pointer to the specified configuration file.
*      1)  looks in the current directory for ".filerc", but will not
*          create a new file in this directory.
*      2)  then looks in the home directory for ".filerc".
*      3)  then looks for an environment variable "NCBI" and takes its
*          value as a complete path to a directory containing the
*          configuration file "filerc" or ".filerc".
*
*      Steps (1) and (2) above are omitted if the NCBI_DONT_USE_LOCAL_CONFIG
*      environment variable is set.  This can be used to allow specific
*      production applications to avoid stray .ncbirc files which may have
*      been erroneously generated.
*
*****************************************************************************/

static FILE *Nlm_OpenConfigFile (Nlm_CharPtr file, Nlm_Boolean writeMode, Nlm_Boolean create)

{
  FILE      *fp;
  Nlm_Int2  i;
  Nlm_Int2  len;
  FILE      *newfp;
  Nlm_Char  path [PATH_MAX+1];
  char      *pth;
  Nlm_Char  str [FILENAME_MAX+1];
  Nlm_Boolean dontUseLocalConfig;

  fp = NULL;
  if (file != NULL) {
    dontUseLocalConfig = getenv("NCBI_DONT_USE_LOCAL_CONFIG") != NULL;
    newfp = NULL;
    Nlm_StringMove(str, ".");
    Nlm_StringNCat(str, file, sizeof(str) - 4);
    if ( ! Nlm_Qualified (str))
    { /* use the user's extension instead of the "rc" extension */
      Nlm_StringCat(str, "rc");
    }
    len = (Nlm_Int2) Nlm_StringLen (str);
    for (i = 0; i < len; i++) {
      str [i] = TO_LOWER (str [i]);
    }
    Nlm_StringNCpy_0(path, str, sizeof(path));

    if (! dontUseLocalConfig)
      fp = Nlm_FileOpen (path, "r");
    if (fp == NULL) {
      if (Nlm_GetHome (path, sizeof (path))) {
        Nlm_FileBuildPath(path, NULL, str);
      } else {
        Nlm_StringNCpy_0(path, str, sizeof(path));
      }
      if (! dontUseLocalConfig)
        fp = Nlm_FileOpen (path, "r");
      if (fp == NULL && create) {
        newfp = Nlm_FileOpen (path, "w");
        Nlm_FileClose (newfp);
        newfp = Nlm_FileOpen (path, "r");
      }
    }
    if (fp == NULL) {
      path[0] = '\0';
      pth = getenv ("NCBI");
      if (pth != NULL) {
        Nlm_FileBuildPath(path, pth, str + 1);
        fp = Nlm_FileOpen (path, "r");
        if (fp == NULL) {
          path[0] = '\0';
          Nlm_FileBuildPath(path, pth, str);
          fp = Nlm_FileOpen (path, "r");
        }
      }
    }
    if (newfp != NULL) {
      if (fp != NULL) {
        Nlm_FileClose (newfp);
        newfp = NULL;
      } else {
        fp = newfp;
      }
    }
    if (writeMode && fp != NULL) {
      Nlm_FileClose (fp);
      fp = Nlm_FileOpen (path, "w");
    }
  }
  return fp;
}

#endif /* OS_UNIX */


#ifdef OS_MAC
/*****************************************************************************
*
*   Nlm_OpenConfigFile (file, writeMode, create)
*      returns a file pointer to the specified configuration file.
*      1)  looks in the System Folder for "file.cnf"
*      2)  then looks in System Folder:Preferences for "file.cnf"
*
*****************************************************************************/

static FILE *Nlm_OpenConfigFile (Nlm_CharPtr file, Nlm_Boolean writeMode, Nlm_Boolean create)

{
  WDPBRec      block;
  Nlm_Char     directory [PATH_MAX];
  long         dirID;
  OSErr        err;
  OSType       fCreator;
  Nlm_Int2     fError;
  FILE         *fp;
  FInfo        finfo;
  OSType       fType;
  long         gesResponse;
  Nlm_Int2     i;
  Nlm_Int2     len;
  CInfoPBRec   params;
  Nlm_Char     str [FILENAME_MAX+1];
  SysEnvRec    sysenv;
  Nlm_Char     temp [PATH_MAX+1];
  Nlm_CharPtr  tmp;
  short        vRefNum;

  fp = NULL;
  if (file != NULL) {
    Nlm_StringNCpy_0(str, file, sizeof(str) - 4);
    if ( ! Nlm_Qualified (str) ) {
      Nlm_StringCat(str, ".cnf");
    }
    len = (Nlm_Int2) Nlm_StringLen (str);
    for (i = 0; i < len; i++) {
      str [i] = TO_LOWER (str [i]);
    }
    if (SysEnvirons (curSysEnvVers, &sysenv) == noErr) {
      block.ioNamePtr = NULL;
      block.ioVRefNum = sysenv.sysVRefNum;
      block.ioWDIndex = 0;
      block.ioWDProcID = 0;
      PBGetWDInfo (&block, FALSE);
      dirID = block.ioWDDirID;
      vRefNum = block.ioWDVRefNum;
      temp [0] = '\0';
      params.dirInfo.ioNamePtr = (StringPtr) directory;
      params.dirInfo.ioDrParID = dirID;
      do {
        params.dirInfo.ioVRefNum = vRefNum;
        params.dirInfo.ioFDirIndex = -1;
        params.dirInfo.ioDrDirID = params.dirInfo.ioDrParID;
        err = PBGetCatInfo (&params, FALSE);
        Nlm_PtoCstr ((Nlm_CharPtr) directory);
        Nlm_StringCat (directory, DIRDELIMSTR);
        Nlm_StringCat (directory, temp);
        Nlm_StringCpy (temp, directory);
      } while (params.dirInfo.ioDrDirID != fsRtDirID);
      tmp = Nlm_StringMove (directory, temp);
      tmp = Nlm_StringMove (tmp, str);
      fp = Nlm_FileOpen (directory, "r");
      if (fp == NULL) {
        if (! Gestalt (gestaltFindFolderAttr, &gesResponse) &&
            (gesResponse & (1 << gestaltFindFolderPresent))) {
          err = FindFolder(kOnSystemDisk, kPreferencesFolderType,
                           kCreateFolder, &vRefNum, &dirID);
          if (err == noErr) {
            params.dirInfo.ioNamePtr = (StringPtr) directory;
            params.dirInfo.ioDrDirID = dirID;
            params.dirInfo.ioVRefNum = vRefNum;
            params.dirInfo.ioFDirIndex = -1;
            err = PBGetCatInfo (&params, FALSE);
            Nlm_PtoCstr ((Nlm_CharPtr) directory);
            Nlm_StringCat (temp, directory);
            Nlm_StringCat (temp, DIRDELIMSTR);
            tmp = Nlm_StringMove (directory, temp);
            tmp = Nlm_StringMove (tmp, str);
          } else {
            tmp = Nlm_StringMove (directory, temp);
            tmp = Nlm_StringMove (tmp, "Preferences");
            tmp = Nlm_StringMove (tmp, DIRDELIMSTR);
            tmp = Nlm_StringMove (tmp, str);
          }
        } else {
          tmp = Nlm_StringMove (directory, temp);
          tmp = Nlm_StringMove (tmp, "Preferences");
          tmp = Nlm_StringMove (tmp, DIRDELIMSTR);
          tmp = Nlm_StringMove (tmp, str);
        }
        fp = Nlm_FileOpen (directory, "r");
      }
      if (fp == NULL && create) {
        tmp = Nlm_StringMove (directory, temp);
        tmp = Nlm_StringMove (tmp, str);
        fp = Nlm_FileOpen (directory, "w");
        Nlm_StringCpy (temp, directory);
        Nlm_CtoPstr ((Nlm_CharPtr) temp);
        fError = GetFInfo ((StringPtr) temp, 0, &finfo);
        if (fError == 0) {
          finfo.fdCreator = 'ttxt';
          finfo.fdType = 'TEXT';
          fError = SetFInfo ((StringPtr) temp, 0, &finfo);
        }
        Nlm_FileClose (fp);
        fp = Nlm_FileOpen (directory, "r");
      }
      Nlm_StringCpy (temp, directory);
      if (writeMode && fp != NULL) {
        Nlm_FileClose (fp);
        Nlm_CtoPstr ((Nlm_CharPtr) temp);
        fType = 'TEXT';
        fCreator = '    ';
        fError = GetFInfo ((StringPtr) temp, 0, &finfo);
        if (fError == 0) {
          fCreator = finfo.fdCreator;
          fType = finfo.fdType;
        }
        fp = Nlm_FileOpen (directory, "w");
        fError = GetFInfo ((StringPtr) temp, 0, &finfo);
        if (fError == 0) {
          finfo.fdCreator = fCreator;
          finfo.fdType = fType;
          fError = SetFInfo ((StringPtr) temp, 0, &finfo);
        }
      }
    }
  }
  return fp;
}
#endif /* OS_MAC */

#ifdef OS_VMS
/*****************************************************************************
*
*   Nlm_GetHome (buf, buflen)
*      returns the path of the home directory
*
*****************************************************************************/

static Nlm_Boolean Nlm_GetHome PROTO((Nlm_CharPtr buf, Nlm_Int2 buflen));
static Nlm_Boolean Nlm_GetHome (Nlm_CharPtr buf, Nlm_Int2 buflen)
{
  Nlm_StringNCpy_0(buf, getenv("SYS$LOGIN"), buflen);
  return TRUE;
}

/*****************************************************************************
*
*   Nlm_OpenConfigFile (file, writeMode, create)
*      returns a file pointer to the specified configuration file.
*      1)  looks in the current directory for "file.cfg", but will not
*          create a new file in this directory.
*      2)  then looks in the home directory for "file.cfg".
*      3)  then looks for an environment variable "NCBI" and takes its
*          value as a complete path to a directory containing the
*          configuration file "file.cfg".
*
*****************************************************************************/

static FILE *Nlm_OpenConfigFile (Nlm_CharPtr file, Nlm_Boolean writeMode, Nlm_Boolean create)

{
  FILE      *fp;
  Nlm_Int2  i;
  Nlm_Int2  len;
  FILE      *newfp;
  Nlm_Char  path [PATH_MAX+1];
  char      *pth;
  Nlm_Char  str [FILENAME_MAX+1];

  fp = NULL;

  if (file != NULL) {
    newfp = NULL;
    Nlm_StringNCpy_0(str, file, sizeof(str) - 4);
    if ( ! Nlm_Qualified (str) ) {
      Nlm_StringCat (str, ".cfg");
    }
    len = (Nlm_Int2) Nlm_StringLen (str);
    for (i = 0; i < len; i++) {
      str [i] = TO_LOWER (str [i]);
    }
    Nlm_StringNCpy_0(path, str, sizeof(path));

    fp = Nlm_FileOpen (path, "r");  /* File exists? */
    if (fp == NULL) {
      if (Nlm_GetHome (path, sizeof (path))) {
        Nlm_FileBuildPath(path, NULL, str);
      } else {
        Nlm_StringNCpy_0(path, str, sizeof(path));
      }
      fp = Nlm_FileOpen (path, "r");   /* File exists? */
      if (fp == NULL && create) {
        newfp = Nlm_FileOpen (path, "w");
        Nlm_FileClose (newfp);
        newfp = Nlm_FileOpen (path, "r");
      }
    }

    if (fp == NULL) {
      path[0] = '\0';
      pth = getenv ("NCBI");
      if (pth != NULL) {
        Nlm_FileBuildPath(path, pth, str);
        fp = Nlm_FileOpen (path, "r");
      }
    }

    if (newfp != NULL) {
      if (fp != NULL) {
        Nlm_FileClose (newfp);
        newfp = NULL;
      } else {
        fp = newfp;
      }
    }

    /*
    ** On VMS if a file is opened for write a new version is created.
    ** This section of code check for "writeMode" and an existing file
    ** if both are true.  Get the currently open file's name and delete
    ** it.  Open a new one in write mode.
    **
    ** Side effects: This will replace the highest existing file version,
    ** but not older version.  There exists the possibility that a user's
    ** custom change may get lost.  A possible workaround for this would
    ** be to have the calling program make a new copy (one higher version)
    ** of the existing file before doing extensive write to the params
    ** file OR keep a static flag in this routine which does  delete the
    ** first time time.
    */

    if (writeMode && fp != NULL) {
      char temp[256];
      fgetname(fp,temp);
      Nlm_FileClose (fp);
      delete(temp);
      fp = Nlm_FileOpen (path, "w");
    }
  }
  return fp;
}

#endif /* OS_VMS */


/*****************************************************************************
*
*   Nlm_TrimString (str)
*      strips trailing spaces, \r, \n
*
*****************************************************************************/

static Nlm_CharPtr Nlm_TrimString (Nlm_CharPtr str)

{
  Nlm_Char     ch;
  Nlm_CharPtr  spc;
  Nlm_CharPtr  tmp;

  if (str != NULL) {
    ch = *str;
    while (ch == ' ' || ch == '\t') {
      str++;
      ch = *str;
    }
    tmp = str;
    spc = NULL;
    ch = *tmp;
    while (ch != '\0' && ch != '\r' && ch != '\n') {
      if (ch == ' ' || ch == '\t') {
        if (spc == NULL) {
          spc = tmp;
        }
      } else {
        spc = NULL;
      }
      tmp++;
      ch = *tmp;
    }
    *tmp = '\0';
    if (spc != NULL) {
      *spc = '\0';
    }
  }
  return str;
}

/*****************************************************************************
*
*   Nlm_ReadConfigFile (fp)
*      reads parameters from configuration file to memory structure
*
*****************************************************************************/

static Nlm_Boolean Nlm_ReadConfigFile (FILE *fp)

{
  Nlm_Char         ch;
  Nlm_env_itemPtr  eip;
  Nlm_env_sectPtr  esp;
  Nlm_env_itemPtr  lastEip;
  Nlm_env_sectPtr  lastEsp;
  Nlm_CharPtr      mid;
  Nlm_Char         str [256];
  Nlm_CharPtr      tmp;
  Nlm_CharPtr      comment;

  if (fp != NULL) {
    if (Nlm_envList != NULL) {
      Nlm_FreeEnvData (Nlm_envList);
      Nlm_envList = NULL;
    }
    esp = NULL;
    lastEsp = NULL;
    eip = NULL;
    lastEip = NULL;
    comment = NULL;
    while (fgets (str, sizeof (str), fp)) {
      ch = *str;
      if (ch != '\n' && ch != '\r') {
        if (ch == ';') { /* comment */
          if (comment == NULL) { /* first comment */
             comment = Nlm_StringSave(str);
          }
          else { /* append to existing comment */
             tmp = (Nlm_CharPtr) Nlm_MemNew(StrLen(comment) + StrLen(str) + 1);
             StrCpy(tmp, comment);
             StrCat(tmp, str);
             comment = (Nlm_CharPtr) Nlm_MemFree(comment);
             comment = tmp;
          }
        } else if (ch == '[') {
          if (esp == NULL) {
            esp = (Nlm_env_sectPtr) Nlm_MemNew (sizeof (Nlm_env_sect));
            lastEsp = esp;
            Nlm_envList = esp;
          } else {
            esp = (Nlm_env_sectPtr) Nlm_MemNew (sizeof (Nlm_env_sect));
            lastEsp->next = esp;
            lastEsp = esp;
          }
          esp->comment = comment;
          comment = NULL;
          tmp = str;
          ch = *tmp;
          while (ch != '\0' && ch != ']') {
            tmp++;
            ch = *tmp;
          }
          *tmp = '\0';
          esp->name = Nlm_StringSave (str + 1);
          eip = NULL;
          lastEip = NULL;
        } else if (esp != NULL) {
          if (eip == NULL) {
            eip = (Nlm_env_itemPtr) Nlm_MemNew (sizeof (Nlm_env_item));
            lastEip = eip;
            esp->children = eip;
          } else {
            eip = (Nlm_env_itemPtr) Nlm_MemNew (sizeof (Nlm_env_item));
            lastEip->next = eip;
            lastEip = eip;
          }
          eip->comment = comment;
          comment = NULL;
          tmp = str;
          mid = str;
          ch = *tmp;
          while (ch != '\0' && ch != '\n' && ch != '\r') {
            if (ch == '=' && mid == str) {
              mid = tmp;
              *mid++ = '\0';
            }
            tmp++;
            ch = *tmp;
          }
          *tmp = '\0';
          eip->name = Nlm_StringSave (Nlm_TrimString (str));
          eip->value = Nlm_StringSave (Nlm_TrimString (mid));
        }
      }
    }

    /* any comments which appeared after the final key of the final section */
    Nlm_bottomComment = comment;
  }
  return TRUE;
}

static Nlm_env_sectPtr Nlm_FindConfigSection (Nlm_CharPtr section)
{
  Nlm_env_sectPtr esp;

  if (section == NULL)
    return NULL;

  for (esp = Nlm_envList; esp != NULL; esp = esp->next)
  {
    if (esp->name != NULL && Nlm_StringICmp(section, esp->name) == 0)
       return esp;
  }

  return NULL;
}

/*****************************************************************************
*
*   Nlm_FindConfigItem (section, type, create)
*      finds parameter in memory structure
*
*****************************************************************************/

static Nlm_env_itemPtr Nlm_FindConfigItem (Nlm_CharPtr section, Nlm_CharPtr type, Nlm_Boolean create)

{
  Nlm_env_itemPtr  eip;
  Nlm_env_sectPtr  esp;
  Nlm_Boolean      goOn;
  Nlm_env_itemPtr  lastEip;
  Nlm_env_sectPtr  lastEsp;

  eip = NULL;
  if (section != NULL && type != NULL) {
    goOn = TRUE;
    esp = Nlm_envList;
    lastEsp = esp;
    while (esp != NULL && goOn) {
      if (esp->name != NULL && Nlm_StringICmp (section, esp->name) == 0) {
        goOn = FALSE;
      } else {
        lastEsp = esp;
        esp = esp->next;
      }
    }
    if (goOn && create) {
      if (Nlm_envList != NULL) {
        esp = (Nlm_env_sectPtr) Nlm_MemNew (sizeof (Nlm_env_sect));
        if (esp != NULL) {
          esp->name = Nlm_StringSave (section);
          lastEsp->next = esp;
        }
      } else {
        esp = (Nlm_env_sectPtr) Nlm_MemNew (sizeof (Nlm_env_sect));
        if (esp != NULL) {
          esp->name = Nlm_StringSave (section);
        }
        Nlm_envList = esp;
      }
    }
    if (esp != NULL) {
      eip = esp->children;
      if (eip != NULL) {
        goOn = TRUE;
        lastEip = eip;
        while (eip != NULL && goOn) {
          if (eip->name != NULL && Nlm_StringICmp (type, eip->name) == 0) {
            goOn = FALSE;
          } else {
            lastEip = eip;
            eip = eip->next;
          }
        }
        if (goOn && create) {
          eip = (Nlm_env_itemPtr) Nlm_MemNew (sizeof (Nlm_env_item));
          if (eip != NULL) {
            eip->name = Nlm_StringSave (type);
            lastEip->next = eip;
          }
        }
      } else if (create) {
        eip = (Nlm_env_itemPtr) Nlm_MemNew (sizeof (Nlm_env_item));
        if (eip != NULL) {
          eip->name = Nlm_StringSave (type);
          esp->children = eip;
        }
      }
    }
  }
  return eip;
}

/*****************************************************************************
*
*   Nlm_WriteConfigFile (fp)
*      writes parameters to configuration file from memory structure
*
*****************************************************************************/

static Nlm_Boolean Nlm_WriteConfigFile (FILE *fp)

{
  Nlm_env_itemPtr  eip;
  Nlm_env_sectPtr  esp;

  if (Nlm_envList != NULL && fp != NULL) {
    esp = Nlm_envList;
    while (esp != NULL) {
      if (! destroyDeadComments || esp->name != NULL)
      {
        Nlm_PutComment (esp->comment, fp);
      }
      if (esp->name != NULL)
      {
        fputc ('[', fp);
        fputs (esp->name, fp);
        fputs ("]\n", fp);
      }
      eip = esp->children;
      while (eip != NULL) {
        if (! destroyDeadComments)
        {
          Nlm_PutComment (eip->comment, fp);
        }
        if (esp->name != NULL && eip->name != NULL && eip->value != NULL) {
          if (destroyDeadComments)
          {
            Nlm_PutComment (eip->comment, fp);
          }
          fputs (eip->name, fp);
          fputc ('=', fp);
          fputs (eip->value, fp);
          fputc ('\n', fp);
        }
        eip = eip->next;
      }
      if (esp->name != NULL)
      {
        fputc ('\n', fp);
      }
      esp = esp->next;
    }
  }

  if (fp != NULL)
    Nlm_PutComment(Nlm_bottomComment, fp);

  return TRUE;
}

/*****************************************************************************
*
*   Nlm_FreeConfigFileData ()
*      frees parameter structure in memory
*
*****************************************************************************/

static void Nlm_FreeConfigFileData (void)

{
  mustCreateFile = TRUE;
  Nlm_bottomComment = (Nlm_CharPtr) Nlm_MemFree(Nlm_bottomComment);
  if (Nlm_lastParamFile != NULL)
    Nlm_lastParamFile = (Nlm_CharPtr) Nlm_MemFree(Nlm_lastParamFile);

  if (Nlm_envList != NULL) {
    Nlm_FreeEnvData (Nlm_envList);
    Nlm_envList = NULL;
  }
}


/*****************************************************************************
*
*   Nlm_FlushConfigFile()
*      flush the specified file's parameters to disk
*
*****************************************************************************/

static void Nlm_FlushConfigFile(Nlm_CharPtr file, Nlm_Boolean create)
{
  FILE *fp;

  if (dirty && file != NULL)
  {
    fp = Nlm_OpenConfigFile (file, TRUE, create);
    if (fp != NULL) {
      Nlm_WriteConfigFile (fp);
      Nlm_FileClose (fp);
    }
  }
  dirty = FALSE;
}

/*****************************************************************************
*
*   Nlm_FreeConfigStruct ()
*      frees parameter structure in memory, and perform other cleanup
*
*****************************************************************************/

static void Nlm_FreeConfigStruct_ST PROTO ((void));
static void Nlm_FreeConfigStruct_ST (void)
{
  Nlm_FlushAppParam();
  Nlm_FreeConfigFileData ();
  Nlm_FreeTransientData ();
}


/*****************************************************************************
*
*   Nlm_PutComment()
*      output a comment to the config file
*
*****************************************************************************/

static void Nlm_PutComment (Nlm_CharPtr s, FILE *fp)

{
  if (s != NULL)
    fputs(s, fp);
}

#else /* ndef OS_MSWIN */

static void        Nlm_FlushAppParam_ST(void) {}
static Nlm_Boolean Nlm_CacheAppParam_ST(Nlm_Boolean value)  { return TRUE; }

/*****************************************************************************
*
* The "guts" of:
*   Nlm_GetAppParam (file, section, type, buf, buflen)
*      finds parameters from configuration files
*      if configuration file is found, trys to read the parameter from it.
*
*****************************************************************************/
static Nlm_Int2 Nlm_WorkGetAppParam (Nlm_CharPtr file, Nlm_CharPtr section, Nlm_CharPtr type, Nlm_CharPtr dflt, Nlm_CharPtr buf, Nlm_Int2 buflen, Nlm_Boolean searchTransient)
{
  static char _empty_string[] = "";
  Nlm_Char path[PATH_MAX + 1];

  if (buf == NULL  ||  buflen <= 0)
    return 0;

  *buf = '\0';
  if (searchTransient  &&
      Nlm_TransientLookup(file, section, type, dflt, buf, buflen))
    {
      return (Nlm_Int2)Nlm_StringLen(buf);
    }

  if ( dflt )
    Nlm_StringNCpy_0(buf, dflt, buflen);

  if (file != NULL && *file != '\0' && section != NULL && *section != '\0')
    {
      Nlm_StringNCpy_0(path, file, sizeof(path) - 4);
      if ( !Nlm_Qualified( path ) )
        Nlm_StringCat(path, ".INI");
      if (dflt == NULL)
        dflt = _empty_string;  /* can't use NULL, must be empty string */
      return (Nlm_Int2)GetPrivateProfileString(section,
                                               type, dflt, buf, buflen, path);
    }
  else
    return (Nlm_Int2)Nlm_StringLen( buf );
}

/*****************************************************************************
*
*   Nlm_SetAppParam (file, section, type, value)
*      finds paths for types of data and fills in path in buf
*      if configuration file is found, trys to write the parameter to it.
*
*****************************************************************************/

static Nlm_Boolean Nlm_SetAppParam_ST PROTO ((Nlm_CharPtr file, Nlm_CharPtr section, Nlm_CharPtr type, Nlm_CharPtr value));
static Nlm_Boolean Nlm_SetAppParam_ST (Nlm_CharPtr file, Nlm_CharPtr section, Nlm_CharPtr type, Nlm_CharPtr value)    /* value */
{
  Nlm_Char     path [PATH_MAX+1];
  Nlm_Boolean  rsult;

  rsult = FALSE;
  if (file != NULL && *file != '\0' && section != NULL && *section != '\0') {
    Nlm_StringNCpy_0(path, file, sizeof(path) - 4);
    if ( ! Nlm_Qualified (path) ) {
      Nlm_StringCat (path, ".INI");
    }
    Nlm_TransientLogSetApp (file, section, type, value);
    if (WritePrivateProfileString (section, type, value, path)) {
      rsult = TRUE;
    }
  }
  return rsult;
}

/*****************************************************************************
*
*   Nlm_FreeConfigStruct ()
*      frees parameter structure in memory
*
*****************************************************************************/

static void Nlm_FreeConfigStruct_ST PROTO ((void));
static void Nlm_FreeConfigStruct_ST (void)
{
  Nlm_FreeTransientData ();
}

#endif /* else !OS_MSWIN */


static Nlm_Boolean Nlm_Qualified( Nlm_CharPtr path )
{
  Nlm_Int4 l,k;
  Nlm_CharPtr  p;

  l = Nlm_StrLen(path);
  p = path+l;
  k = 0;
  while (k < l && k <= 4) {
     if (*p-- == '.') return TRUE;
     k++;
  }
  return FALSE;
}


/*****************************************************************************
*
*   Nlm_FindPath (file, section, type, buf, buflen)
*      finds paths for types of data from configuration files.
*      if configuration file is found, tries to read the parameter from it,
*      then appends a directory delimiter character, if necessary.
*
*****************************************************************************/

NLM_EXTERN Nlm_Boolean LIBCALL Nlm_FindPath (Nlm_CharPtr file, Nlm_CharPtr section, Nlm_CharPtr type, Nlm_CharPtr buf, Nlm_Int2 buflen)   /* length of path buffer */
{
  Nlm_Boolean rsult = FALSE;

  if (file == NULL  ||  section == 0  ||  type == NULL  ||
      buf == NULL  ||  buflen <= 0)
    return FALSE;

  NlmMutexLockEx( &corelibMutex );

  *buf = '\0';
  if (*file != '\0'  &&  *section != '\0'  &&  *type != '\0'  &&
      Nlm_GetAppParam(file, section, type, "", buf, (Nlm_Int2)(buflen - 1))  &&
      *buf != '\0')
    {
      Nlm_FileBuildPath(buf, NULL, NULL);
      rsult = TRUE;
    }

  NlmMutexUnlock( corelibMutex );
  return rsult;
}


static Nlm_Boolean Nlm_TransientSetAppParam_ST PROTO ((Nlm_CharPtr file, Nlm_CharPtr section, Nlm_CharPtr type, Nlm_CharPtr value));
static Nlm_Boolean Nlm_TransientSetAppParam_ST (Nlm_CharPtr file, Nlm_CharPtr section, Nlm_CharPtr type, Nlm_CharPtr value)
{
  Nlm_env_filePtr  theFile;
  Nlm_env_itemPtr  eip;
  Nlm_env_sectPtr  esp;
  Nlm_env_itemPtr  nextEip;

  if (file == NULL || *file == '\0' || section == NULL || *section == '\0')
    return FALSE;

  for (theFile = Nlm_transientFileList; theFile != NULL; theFile = theFile->next)
  {
    if (StringICmp(theFile->name, file) == 0)
    {
      for (esp = theFile->envList; esp != NULL; esp = esp->next)
      {
        if (esp->name != NULL && StringICmp(esp->name, section) == 0)
        {
          if (type == NULL || type[0] == '\0')
          {
            /* free all children */
            for (eip = esp->children; eip != NULL; eip = nextEip)
            {
              nextEip = eip->next;
              eip->next    = (Nlm_env_itemPtr)(-1);
              eip->name    = (Nlm_CharPtr) Nlm_MemFree (eip->name);
              eip->comment = (Nlm_CharPtr) Nlm_MemFree (eip->comment);
              eip->value   = (Nlm_CharPtr) Nlm_MemFree (eip->value);
              Nlm_MemFree (eip);
            }
            esp->children = NULL;
            esp->transientOnly = TRUE;
          } else { /* append this type to the section */
            eip = (Nlm_env_itemPtr) MemNew(sizeof(*eip));
            eip->name = StringSave(type);
            eip->value = StringSave(value);
            eip->next = esp->children;
            esp->children = eip;
          }
          return TRUE;
        }
      }
      break;
    }
  }

  /* create the file data structure if needed */
  if (theFile == NULL)
  {
    theFile = (Nlm_env_filePtr) MemNew(sizeof(*theFile));
    theFile->name = StringSave(file);
    theFile->next = Nlm_transientFileList;
    Nlm_transientFileList = theFile;
  }

  /* create the section and type */
  esp = (Nlm_env_sectPtr) MemNew(sizeof(*esp));
  esp->name = StringSave(section);
  esp->next = theFile->envList;
  theFile->envList = esp;
  if (type == NULL || type[0] == '\0')
  {
    esp->transientOnly = TRUE;
  } else { /* create the section */
    esp->transientOnly = FALSE;
    eip = (Nlm_env_itemPtr) MemNew(sizeof(*eip));
    eip->name = StringSave(type);
    eip->value = StringSave(value);
    eip->next = NULL;
    esp->children = eip;
  }

  return TRUE;
}


/* SetAppParam is writing a value to the real config file, so log this value,
   if necessary, into the "transient" data structures */
static void Nlm_TransientLogSetApp (Nlm_CharPtr file, Nlm_CharPtr section, Nlm_CharPtr type, Nlm_CharPtr value)

{
  Nlm_env_filePtr  theFile;
  Nlm_env_itemPtr  eip;
  Nlm_env_sectPtr  esp;

  if (file == NULL || *file == '\0' || section == NULL || *section == '\0')
    return;

  if (type == NULL || type[0] == '\0')
  {
    for (theFile = Nlm_transientFileList; theFile != NULL; theFile = theFile->next)
    {
      if (StringICmp(theFile->name, file) == 0)
      {
        for (esp = theFile->envList; esp != NULL; esp = esp->next)
        {
          if (esp->name != NULL && StringICmp(esp->name, section) == 0)
          { /* delete the section by removing section name */
            esp->name = (Nlm_CharPtr) MemFree(esp->name);
          }
        }
      }
    }
  } else {
    for (theFile = Nlm_transientFileList; theFile != NULL; theFile = theFile->next)
    {
      if (StringICmp(theFile->name, file) == 0)
      {
        for (esp = theFile->envList; esp != NULL; esp = esp->next)
        {
          if (esp->name != NULL && StringICmp(esp->name, section) == 0 &&
              esp->transientOnly)
          { /* append this type to the section */
            eip = (Nlm_env_itemPtr) MemNew(sizeof(*eip));
            eip->name = StringSave(type);
            eip->value = StringSave(value);
            eip->next = esp->children;
            esp->children = eip;
          }
        }
      }
    }
  }
}

static Nlm_Boolean Nlm_TransientLookup (Nlm_CharPtr file, Nlm_CharPtr section, Nlm_CharPtr type, Nlm_CharPtr dflt, Nlm_CharPtr buf, Nlm_Int2 buflen)

{
  Nlm_env_filePtr  theFile;
  Nlm_env_itemPtr  eip;
  Nlm_env_sectPtr  esp;
  Nlm_Int2         totlen;
  Nlm_Int4         bytesToAppend;

  if (file == NULL || *file == '\0' || section == NULL || *section == '\0')
    return FALSE;

  for (theFile = Nlm_transientFileList; theFile != NULL; theFile = theFile->next)
  {
    if (StringICmp(theFile->name, file) == 0)
    {
      for (esp = theFile->envList; esp != NULL; esp = esp->next)
      {
        if (esp->name != NULL && StringICmp(esp->name, section) == 0)
        {
          if (type == NULL || type[0] == '\0')
          { /* concatenate all types (keys) within section */
            *buf = '\0';
            totlen = 0;
            for (eip = esp->children; eip != NULL; eip = eip->next)
            {
              bytesToAppend = StrLen(eip->name) + 1;
              bytesToAppend = MIN(bytesToAppend, buflen - totlen);
              StrNCpy(&buf[totlen], eip->name, bytesToAppend);
              totlen += bytesToAppend;
            }
            if (totlen > 0 && buf[totlen] == '\0')
            {
                totlen--; /* account for final null character */
            }
            /* now append the GetAppParam() data */
            if (! esp->transientOnly)
            { /* GetAppParam data can be trusted ... append it to buf */
              Nlm_WorkGetAppParam(file, section, NULL, "", &buf[totlen],
                                  (Nlm_Int2)(buflen - totlen), FALSE);
            }
            return TRUE;
          } else {
            for (eip = esp->children; eip != NULL; eip = eip->next)
            {
              if (StringICmp(eip->name, type) == 0)
              {
                Nlm_StringNCpy_0(buf, eip->value, buflen);
                return TRUE;
              }
            }
            if (esp->transientOnly)
            { /* GetAppParam data cannot be trusted ... use the default */
              Nlm_StringNCpy_0(buf, dflt, buflen);
              return TRUE;
            }
          }
        }
      }
    }
  }

  /* not found ... GetAppParam() should search the real config file */
  return FALSE;
}

static void Nlm_FreeEnvData (Nlm_env_sectPtr esp)

{
  Nlm_env_itemPtr  eip;
  Nlm_env_itemPtr  nextEip;
  Nlm_env_sectPtr  nextEsp;

  while (esp != NULL) {
    nextEsp = esp->next;
    eip = esp->children;
    while (eip != NULL) {
      nextEip = eip->next;
      eip->next    = (Nlm_env_itemPtr)(-1);
      eip->name    = (Nlm_CharPtr) Nlm_MemFree (eip->name);
      eip->comment = (Nlm_CharPtr) Nlm_MemFree (eip->comment);
      eip->value   = (Nlm_CharPtr) Nlm_MemFree (eip->value);
      Nlm_MemFree (eip);
      eip = nextEip;
    }
    esp->next    = (Nlm_env_sectPtr)(-1);
    esp->name    = (Nlm_CharPtr) Nlm_MemFree (esp->name);
    esp->comment = (Nlm_CharPtr) Nlm_MemFree (esp->comment);
    Nlm_MemFree (esp);
    esp = nextEsp;
  }
}


static void Nlm_FreeTransientData (void)
{
  Nlm_env_filePtr efp, nextEfp;

  efp = Nlm_transientFileList;
  while (efp != NULL) {
    nextEfp = efp->next;
    Nlm_FreeEnvData (efp->envList);
    efp->envList = NULL;
    efp->next = (Nlm_env_filePtr)(-1);
    efp->name = (Nlm_CharPtr) Nlm_MemFree (efp->name);
    Nlm_MemFree (efp);
    efp = nextEfp;
  }
  Nlm_transientFileList = NULL;
}

/*****************************************************************************
*
*   Command-line arguments
*     &
*   Nlm_ProgramPath (buf, maxsize)
*   	returns full path to executing program
*
*****************************************************************************/

static Nlm_Boolean wasSetup = FALSE;
static int    targc = 0;
static char **targv = NULL;


#ifdef WIN_MAC
#ifdef __CONDITIONALMACROS__
static FSSpec       apFileSpec;
#endif
static Str255       apName;
static Handle       apParam;
static short        apRefNum;

static Nlm_Boolean Nlm_SetupArguments_ST_Mac (void)
{
#ifdef __CONDITIONALMACROS__
  long                 gval;
  ProcessInfoRec       pirec;
  ProcessSerialNumber  psn;

  if (Gestalt (gestaltSystemVersion, &gval) == noErr && (short) gval >= 7 * 256) {
    GetCurrentProcess (&psn);
    pirec.processInfoLength = sizeof (ProcessInfoRec);
    pirec.processName = NULL;
    pirec.processAppSpec = &apFileSpec;
    GetProcessInformation (&psn, &pirec);
    Nlm_PtoCstr ((Nlm_CharPtr) apFileSpec.name);
  } else {
#ifdef PROC_MC680X0
    GetAppParms (apName, &apRefNum, &apParam);
    Nlm_PtoCstr ((Nlm_CharPtr) apName);
#else
    return FALSE;
#endif
  }
#else
  GetAppParms (apName, &apRefNum, &apParam);
  Nlm_PtoCstr ((Nlm_CharPtr) apName);
#endif

  SetAppProperty("ProgramName",(void*)apName);
  return TRUE;
}


static void Nlm_ProgramPath_ST PROTO ((Nlm_CharPtr buf, size_t maxsize));
static void Nlm_ProgramPath_ST (Nlm_CharPtr buf, size_t maxsize)
{
#ifdef __CONDITIONALMACROS__
  CInfoPBRec  block;
  Nlm_Char    path [256];
  Nlm_Char    temp [256];
  short nErr;

  if (buf != NULL && maxsize > 0) {
    *buf = '\0';
    if (wasSetup) {
      memset (&block, 0, sizeof (CInfoPBRec));
      Nlm_StringNCpy_0(path, (Nlm_CharPtr)apFileSpec.name, sizeof (path));

      block.dirInfo.ioNamePtr = (StringPtr) path;
      block.dirInfo.ioDrParID = apFileSpec.parID;

      do {
        Nlm_StringCpy (temp, path);
        block.dirInfo.ioVRefNum = apFileSpec.vRefNum;
        block.dirInfo.ioFDirIndex = -1;
        block.dirInfo.ioDrDirID = block.dirInfo.ioDrParID;
        nErr = PBGetCatInfo (&block, FALSE);
        if (nErr != noErr) break;
        Nlm_PtoCstr ((Nlm_CharPtr) path);
        Nlm_StringCat (path, DIRDELIMSTR);
        Nlm_StringCat (path, temp);
      } while (block.dirInfo.ioDrDirID != fsRtDirID);

      Nlm_StringNCpy_0(buf, path, maxsize);
    }
  }
#else
  WDPBRec      block;
  Nlm_Char     directory [PATH_MAX];
  Nlm_Int4     dirID;
  OSErr        err;
  CInfoPBRec   params;
  Nlm_Char     temp [PATH_MAX];
  Nlm_CharPtr  tmp;
  short        vRefNum;

  if (buf != NULL && maxsize > 0) {
    *buf = '\0';
    if (wasSetup) {
      memset (&block, 0, sizeof (WDPBRec));
      block.ioNamePtr = NULL;
      block.ioVRefNum = apRefNum;
      block.ioWDIndex = 0;
      block.ioWDProcID = 0;
      PBGetWDInfo (&block, FALSE);
      dirID = block.ioWDDirID;
      vRefNum = block.ioWDVRefNum;
      temp [0] = '\0';
      params.dirInfo.ioNamePtr = (StringPtr) directory;
      params.dirInfo.ioDrParID = dirID;
      do {
        params.dirInfo.ioVRefNum = vRefNum;
        params.dirInfo.ioFDirIndex = -1;
        params.dirInfo.ioDrDirID = params.dirInfo.ioDrParID;
        err = PBGetCatInfo (&params, FALSE);
        Nlm_PtoCstr ((Nlm_CharPtr) directory);
        Nlm_StringCat (directory, DIRDELIMSTR);
        Nlm_StringCat (directory, temp);
        Nlm_StringCpy (temp, directory);
      } while (params.dirInfo.ioDrDirID != fsRtDirID);
      tmp = Nlm_StringMove (directory, temp);
      tmp = Nlm_StringMove (tmp, (Nlm_CharPtr) apName);
      Nlm_StringNCpy_0(buf, directory, maxsize);
    }
  }
#endif
}
#endif /* WIN_MAC */


#if defined(OS_MSWIN) || defined(OS_VMS)
static void Nlm_ProgramPath_ST (Nlm_CharPtr buf, size_t maxsize)
{
  if (!buf  ||  maxsize <= 0)
    return;

  *buf = '\0';
  if (wasSetup  &&  targv  &&  targv[0])
    Nlm_StringNCpy_0(buf, targv[0], maxsize);
}
#endif  /* OS_MSWIN || OS_VMS */


#ifdef OS_UNIX
static void Nlm_ProgramPath_ST (Nlm_CharPtr buf, size_t maxsize)
{
  Nlm_Char     path [PATH_MAX];
  Nlm_CharPtr  pth;
  Nlm_CharPtr  ptr;

  if (buf != NULL && maxsize > 0) {
    *buf = '\0';
    if (wasSetup) {
      ptr = targv [0];
      if (ptr [0] == DIRDELIMCHR) {
        Nlm_StringNCpy_0(buf, targv[0], maxsize);
      } else if (getcwd (path, sizeof (path)) != NULL) {
        ptr = targv [0];
        while (ptr [0] == '.' || ptr [0] == DIRDELIMCHR) {
          if (ptr [0] == '.') {
            if (ptr [1] == '.' && ptr [2] == DIRDELIMCHR) {
              ptr += 3;
              pth = StringRChr (path, DIRDELIMCHR);
              if (pth != NULL) {
                *pth = '\0';
              }
            } else if (ptr [1] == DIRDELIMCHR) {
              ptr += 2;
            } else {
              ptr++;
            }
          } else if (ptr [0] == DIRDELIMCHR) {
            ptr++;
          } else {
            ptr++;
          }
        }
        FileBuildPath (path, NULL, ptr);
        Nlm_StringNCpy_0(buf, path, maxsize);
      } else {
        Nlm_StringNCpy_0(buf, targv[0], maxsize);
      }
    }
  }
}
#endif


/*****************************************************************************
* Multi-Thread protected external functions
*****************************************************************************/

NLM_EXTERN Nlm_Boolean LIBCALL Nlm_TransientSetAppParam (Nlm_CharPtr file, Nlm_CharPtr section, Nlm_CharPtr type, Nlm_CharPtr value)
{
  Nlm_Boolean rsult;
  NlmMutexLockEx( &corelibMutex );

  rsult = Nlm_TransientSetAppParam_ST(file, section, type, value);

  NlmMutexUnlock( corelibMutex );
  return rsult;
}

NLM_EXTERN void LIBCALL Nlm_FreeConfigStruct(void)
{
  NlmMutexLockEx( &corelibMutex );
  Nlm_FreeConfigStruct_ST();
  NlmMutexUnlock( corelibMutex );
}


NLM_EXTERN void LIBCALL Nlm_SetupArguments(int argc, char *argv[])
{
  NlmMutexLockEx( &corelibMutex );
  wasSetup = TRUE;
#if defined(OS_UNIX)
  {{
    char *p;
    if ((p = strrchr(argv[0],DIRDELIMCHR)) != NULL)  p++;
    else
      p = argv[0];
    SetAppProperty("ProgramName", (void*)p);  
  }}
#elif defined(WIN_MAC)
  wasSetup = Nlm_SetupArguments_ST_Mac();
#endif
  targc = argc;
  targv = argv;
  NlmMutexUnlock( corelibMutex );
}

NLM_EXTERN Nlm_CharPtr PNTR Nlm_GetArgv(void)
{
  return targv;
}

NLM_EXTERN Nlm_Int4 Nlm_GetArgc(void)
{
  return targc;
}

NLM_EXTERN void LIBCALL Nlm_ProgramPath(Nlm_CharPtr buf, size_t maxsize)
{
  NlmMutexLockEx( &corelibMutex );
  Nlm_ProgramPath_ST(buf, maxsize);
  NlmMutexUnlock( corelibMutex );
}

NLM_EXTERN void LIBCALL Nlm_FlushAppParam (void)
{
  NlmMutexLockEx( &corelibMutex );
  Nlm_FlushAppParam_ST();
  NlmMutexUnlock( corelibMutex );
}

NLM_EXTERN Nlm_Boolean LIBCALL Nlm_CacheAppParam(Nlm_Boolean value)
{
  Nlm_Boolean rsult;
  NlmMutexLockEx( &corelibMutex );

  rsult = Nlm_CacheAppParam_ST( value );

  NlmMutexUnlock( corelibMutex );
  return rsult;
}

NLM_EXTERN Nlm_Int2 LIBCALL Nlm_GetAppParam (Nlm_CharPtr file, Nlm_CharPtr section, Nlm_CharPtr type, Nlm_CharPtr dflt, Nlm_CharPtr buf, Nlm_Int2 buflen)
{
  Nlm_Int2 rsult;
  NlmMutexLockEx( &corelibMutex );

  rsult = Nlm_WorkGetAppParam(file, section, type, dflt, buf, buflen, TRUE);

  NlmMutexUnlock( corelibMutex );
  return rsult;
}

NLM_EXTERN Nlm_Boolean LIBCALL Nlm_SetAppParam (Nlm_CharPtr file, Nlm_CharPtr section, Nlm_CharPtr type, Nlm_CharPtr value)
{
  Nlm_Boolean rsult;
  NlmMutexLockEx( &corelibMutex );

  rsult = Nlm_SetAppParam_ST(file, section, type, value);

  NlmMutexUnlock( corelibMutex );
  return rsult;
}


/*****************************************************************************
*
*   GetAppParamBoolean()
*       SetAppParamBoolean()
*   GetAppParamLong()
*       SetAppParamLong()
*
*****************************************************************************/

NLM_EXTERN Nlm_Boolean LIBCALL GetAppParamBoolean (const char *filebase, const char *sect,
			const char *key, Nlm_Boolean dflt)
{
	char buffer[32];
	if (GetAppParam((char*)filebase,(char*)sect,(char*)key,"",buffer,sizeof buffer))
	{
		if (strchr("1yYtT",buffer[0]))
			return TRUE;
		if (strchr("0nNfF",buffer[0]))
			return FALSE;
	}
	return dflt;
}

NLM_EXTERN Nlm_Boolean LIBCALL SetAppParamBoolean (const char *filebase, const char *sect,
			const char *key, Nlm_Boolean value)
{
	return SetAppParam((char*)filebase,(char*)sect,(char*)key,(value)?"Yes":"No");
}

NLM_EXTERN long LIBCALL GetAppParamLong (const char *filebase, const char *sect,
			const char *key, long dflt)
{
	char buffer[32];
	return (GetAppParam((char*)filebase,(char*)sect,(char*)key,"",buffer,sizeof buffer)) ?
				atol(buffer) : dflt;
}

NLM_EXTERN Nlm_Boolean LIBCALL SetAppParamLong (const char *filebase, const char *sect,
			const char *key, long value)
{
	char buffer[32];
	sprintf(buffer,"%ld",value);
	return SetAppParam((char*)filebase,(char*)sect,(char*)key,buffer);
}


#ifdef WIN32
extern int Nlm_x_HasConsole(void)
{
  static int has_console = -1;
  if (has_console == -1)
    has_console = fileno(stdin) >= 0 ? 1 : 0;

  return has_console;
}
#endif


NLM_EXTERN Nlm_Boolean Nlm_FreeArgs(Nlm_Int2 numargs, Nlm_ArgPtr ap)
{
  Nlm_Int2 i;
  for (i = 0;  i < numargs;  i++, ap++)
    {
      switch ( ap->type )
        {
        case ARG_BOOLEAN:
          ap->intvalue = 0;
          break;
        case ARG_INT:
          ap->intvalue = 0;
          break;
        case ARG_FLOAT:
          ap->floatvalue = 0.0;
          break;
        case ARG_STRING:
        case ARG_FILE_IN:
        case ARG_FILE_OUT:
        case ARG_DATA_IN:
        case ARG_DATA_OUT:
          ap->strvalue = (Nlm_CharPtr) Nlm_MemFree( ap->strvalue );
          break;
        default:
          ASSERT ( FALSE );
          return FALSE;
        }
    }
  return TRUE;
}
