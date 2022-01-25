/*   cn3dmain.c
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*            National Center for Biotechnology Information (NCBI)
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government do not place any restriction on its use or reproduction.
*  We would, however, appreciate having the NCBI and the author cited in
*  any work or product based on this material
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
* ===========================================================================
*
* File Name:  cn3dmain.c
*
* Author:  Christopher Hogue
*
* Version Creation Date:   1/31/96
*
* File Description: Main entry point for Cn3d
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: cn3dmain.c,v $
* Revision 6.27  1998/08/03 18:33:12  lewisg
* added netentcf to cn3d
*
* Revision 6.26  1998/06/17 22:41:26  kans
* now uses UseLocalAsnloadDataAndErrMsg from sqnutils.h
*
* Revision 6.25  1998/06/15 14:26:06  ywang
* automatic launch salsa when mime data got in either through command line or via local reading in
*
 * Revision 6.24  1998/06/11  20:07:38  chappey
 * Added SetAppProperty ("SeqEditDisplayForm") to pass flags to sequence editor
 *
* Revision 6.23  1998/06/05 21:14:22  ywang
* A solution for ensure getting bsp
*
 * Revision 6.22  1998/06/04  16:48:34  ywang
 * fix bug triggered by automatic salsa launch
 *
 * Revision 6.20  1998/05/06  23:50:21  lewisg
 * fixed launching problem with sequin
 *
* Revision 6.19  1998/04/28 19:39:58  lewisg
* codewarrior fixes
*
* Revision 6.18  1998/04/28 19:38:35  lewisg
* codewarrior fixes
*
* Revision 6.17  1998/04/28 18:54:01  ywang
* slight modification
*
 * Revision 6.15  1998/04/28  15:14:25  lewisg
 * moved OpenMimeFileWithDeletion to cn3dopen
 *
* Revision 6.14  1998/04/27 23:23:02  lewisg
* added ability to open mime files
*
* Revision 6.13  1998/04/21 23:00:56  lewisg
* added show aligned/unaligned
*
* Revision 6.12  1998/04/15 00:51:36  lewisg
* bug fixes for multiple alignment mode and alignment pane
*
* Revision 6.11  1998/04/04 18:07:45  lewisg
* get rid of typo
*
* Revision 6.10  1998/04/04 05:57:52  lewisg
* got rid of dos line breaks
*
* Revision 6.9  1998/04/04 00:53:45  lewisg
* added support for multiple alignments
*
* Revision 6.8  1998/04/01 23:26:16  lewisg
* added new startup mode + fixed slave rendering
*
* Revision 6.7  1998/03/26 22:42:10  lewisg
* added seqentry and seq annot to msd
*
* Revision 6.6  1998/03/13 22:30:34  lewisg
* fix neighbor mode
*
* Revision 6.5  1998/03/07 20:43:51  kans
* moved Cn3D_fEntrezOn to cn3dwin.c
*
* Revision 6.4  1998/03/06 23:19:14  lewisg
* codewarrior fixes
*
* Revision 6.3  1998/03/06 01:19:32  lewisg
* merge
*
* Revision 6.2  1997/10/09 13:01:54  epstein
*  add return values for OpenMimeFile
*
* Revision 6.1  1997/09/30 20:09:21  epstein
* ADD ABILITY TO PERFORM cN3d demos driven from configuration files
*
* Revision 6.0  1997/08/25 18:13:31  madden
* Revision changed to 6.0
*
* Revision 5.20  1997/07/16 20:55:48  vakatov
* Use "Nlm_GetArg[vc]()" instead of "arg[vc]"
*
* Revision 5.19  1997/03/20 19:04:08  vakatov
* Now contains only standalone-specific code;  the generic Cn3D code has
* been moved to "cn3dwin.c", and the Entrez-specific code -- to "cn3dentr.c".
*
*
* ==========================================================================
*/

#include <ncbi.h>
#include <vibrant.h>  /* for netentcf */
#include <netcnfg.h>
#include <cn3dmain.h>
#include <objmime.h>
#include <accentr.h>
#include <objalign.h>
#include <objseq.h>
#include <objmgr.h>
#include <sequtil.h>
#include <saledit.h>
#include <lsqfetch.h>
#include <cn3dpane.h>
#include <algorend.h>
#include <cn3dopen.h>
#include <cn3dmsg.h>
#include <sqnutils.h>

extern Boolean viewalign_only;

static Boolean LIBCALLBACK OpenMimeFile(CharPtr filename)
{
#ifdef OS_MAC
  /* the Web browsers on other platforms get upset if you delete the file, */
  /* but apparently on the Mac you must delete it */
  return OpenMimeFileWithDeletion(filename, TRUE);
#else
  return OpenMimeFileWithDeletion(filename, FALSE);
#endif
}

/******* SEQUENCE EDITOR *********/
static SeqEditViewProcs    seqedprocs;
/******** END ************/

/* used by netentcf */

static IteM  cutItem = NULL;
static IteM  copyItem = NULL;
static IteM  pasteItem = NULL;
static IteM  deleteItem = NULL;

/* no idea what this mac stuff is for */
#ifdef WIN_MAC
static void ConfigFormActivated (WindoW w)

{
  currentFormDataPtr = (VoidPtr) GetObjectExtra (w);
  RepeatProcOnHandles (Enable,
                   (HANDLE) cutItem,
                   (HANDLE) copyItem,
                   (HANDLE) pasteItem,
                   (HANDLE) deleteItem,
                   NULL);
}
#endif

#ifndef WIN_MAC
#define ConfigFormActivated NULL
#endif

static void ConfigFormMessage (ForM f, Int2 mssg)

{
  BaseFormPtr  bfp;

  bfp = (BaseFormPtr) GetObjectExtra (f);
  if (bfp != NULL) {
    switch (mssg) {
      case VIB_MSG_CUT :
        StdCutTextProc (NULL);
        break;
      case VIB_MSG_COPY :
        StdCopyTextProc (NULL);
        break;
      case VIB_MSG_PASTE :
        StdPasteTextProc (NULL);
        break;
      case VIB_MSG_DELETE :
        StdDeleteTextProc (NULL);
        break;
      default :
        break;
    }
  }
}

static void ConfigAccepted (void)

{
  Message (MSG_OK, "Please rerun Cn3D now that it is configured");
  QuitProgram ();
}

static void ConfigCancelled (void)

{
  Message (MSG_OK, "Configuration cancelled.  Exiting program.");
  QuitProgram ();
}

static CharPtr configMessage =
"Entrez cannot find ncbi configuration file.\n\
Do you wish to run the configuration program?";

/* end of stuff used by netentcf */



Int2 Main(void)
{
  char buffer [PATH_MAX];
  WindoW www;
  Char           str [16];



  ErrSetFatalLevel( SEV_MAX );

if (! GetAppParam ("NCBI", "NCBI", NULL, NULL, str, sizeof (str) - 1)) {
    if (Message (MSG_YN, configMessage) == ANS_YES) {
      ShowNetConfigForm (ConfigFormActivated, ConfigFormMessage,
                         ConfigAccepted, ConfigCancelled);
      ProcessEvents ();
      return 0;
    }
}


  UseLocalAsnloadDataAndErrMsg();
  if ( !OpenMMDBAPI(0, NULL) )
    return 1;

/******* SEQUENCE EDITOR *********/
  MemSet ((Pointer) (&seqedprocs), 0, sizeof (SeqEditViewProcs));
  SetAppProperty ("SeqEditDisplayForm", &seqedprocs);
/******** END ************/
  

#ifdef WIN_MAC
    DeskAccGroup( AppleMenu( NULL ) );
#endif

  www = Cn3DWin(NULL, NULL);
  if ( www )
    {
      if (GetAppParam("Cn3D","demo","mandatory_file","", buffer, sizeof(buffer)) > 0)
	{
	  OpenMimeFileWithDeletion(buffer, FALSE);
	} else {
#if defined(OS_UNIX) || defined(WIN_MSWIN)
	  if (GetArgc() == 2)
	    OpenMimeFile( GetArgv()[1] );
	  else {
	    if (GetAppParam("Cn3D","demo","optional_file","", buffer, sizeof(buffer)) > 0)
	      {
		OpenMimeFileWithDeletion(buffer, FALSE);
	      }
	  }
#endif
#if defined(WIN_MAC)  ||  defined(WIN_MSWIN)
	  RegisterDropProc( OpenMimeFile );
#endif
	}

      if(!viewalign_only) Show( www );
/*    if(Mime_ReadIn) LaunchSequenceWindow(); */    /*  yanli */
      ProcessEvents();
  }

  CloseMMDBAPI();
  
  if(Cn3D_fEntrezOn)
  {
    EntrezBioseqFetchDisable();
/*  BioseqFetchDisable();
    EntrezFini();   */
  }
  
  return 0;
}

