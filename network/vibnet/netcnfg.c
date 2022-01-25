/*   netcnfg.c
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
 * RCS $Id: netcnfg.c,v 6.7 1998/08/28 18:53:25 vakatov Exp $
 *
 * Author:  Kans, Epstein
 *
 * Version Creation Date:   9/10/96
 *
 * File Description:
 *       Network Entrez configuration
 *
 * Modifications:
 * --------------------------------------------------------------------------
 * Date     Name        Description of modification
 * -------  ----------  -----------------------------------------------------
 */

#include <vibrant.h>
#include <document.h>
#include <accentr.h>
#include <ncbinet.h>
#include <netentr.h>
#include <netcnfg.h>

#ifdef NETP_INET_NEWT
#define SIN_ADDR        sin_addr.S_un.S_addr
#else
#define SIN_ADDR        sin_addr
#endif

#define FIRST_PAGE   0
#define SECOND_PAGE  1
#define THIRD_PAGE   2
#define FOURTH_PAGE  3
#define FIFTH_PAGE   4
#define NUM_PAGES    5

#define IDENTIFYING_STRING "NCBI Configuration"
#define SERVICE_TYPE "Entrez"


struct DispAddr {
  CharPtr         fqdn;
  CharPtr         addr;
};

struct DispAddr dispAddrStrs[] = {
  "dispatch1.nlm.nih.gov", "130.14.25.211",
  "dispatch2.nlm.nih.gov", "130.14.25.47",
  "dispatch3.nlm.nih.gov", "130.14.25.1"
};

#define NUMDISP ((sizeof(dispAddrStrs)/sizeof(dispAddrStrs[0])) * 2)

typedef struct netconfigdata {
  FORM_MESSAGE_BLOCK

  GrouP pages[NUM_PAGES];
  Int2 currentPage;
  Boolean movingForward;

  TexT            asnloadPath;
  TexT            dataPath;
  TexT            errmsgPath;
  Boolean         page1ok;

  TexT            username;
  GrouP           ifCantFindDisp;
  ButtoN          wantEncryption;
  ButtoN          outgoingOnly;

  GrouP           dispList;
  TexT            dispAddr;
  ButtoN          testDisp;
  ValNodePtr      dispTitles;

  LisT            servList;
  PrompT          servDesc;
  PrompT          datasetDesc;
  ButtoN          testServ;
  ValNodePtr      servTitles;

  PrompT          fyi[3];

  ButtoN          nextBtn;
  ButtoN          prevBtn;

  VoidProc        accepted;
  VoidProc        cancelled;

  NICatalog      *catalog;
  NI_DispatcherPtr dispatcher;
  NI_DispInfoPtr  dispInfoPtr;

  NIToolset      *toolset[16];
  NIService      *svc[16];
  NIResource     *rsc;

  Int2            svc_ct;
  Int2            svc_sel;
  Int2            verSvcMin;
  Int2            verSvcMax;
  Int2            verRscMin;
  Int2            verRscMax;
  CharPtr         svname;
  CharPtr         resname;

  CharPtr         adminInfo;

  ENIInterface    prev_interface;
}               NetConfigData, PNTR NetConfigPtr;

static void     ConfigMessageProc (ForM f, Int2 mssg)
{
  VoidProc        cancelled;
  NetConfigPtr    ncp;

  ncp = (NetConfigPtr) GetObjectExtra (f);
  if (ncp != NULL) {
    switch (mssg) {
      case VIB_MSG_CLOSE:
        cancelled = ncp->cancelled;
        Remove (f);
        if (cancelled != NULL) {
          cancelled ();
        }
        NI_SetInterface(ncp->prev_interface);
        break;
      case VIB_MSG_CUT:
        StdCutTextProc (NULL);
        break;
      case VIB_MSG_COPY:
        StdCopyTextProc (NULL);
        break;
      case VIB_MSG_PASTE:
        StdPasteTextProc (NULL);
        break;
      case VIB_MSG_DELETE:
        StdDeleteTextProc (NULL);
        break;
      default:
        if (ncp->appmessage != NULL) {
          ncp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static void     ConfigFormActivate (WindoW w)
{
  NetConfigPtr    ncp;

  ncp = (NetConfigPtr) GetObjectExtra (w);
  if (ncp != NULL) {
    if (ncp->activate != NULL) {
      ncp->activate (w);
    }
  }
}

#ifdef OS_MAC
#define ASNLOAD_NEEDED 1
#endif

#if defined(OS_DOS) || defined(WIN16)
#define ASNLOAD_NEEDED 1
#endif


static Boolean  FileExists (CharPtr dirname, CharPtr subname, CharPtr filename)
{
  Char            path[PATH_MAX];

  StringNCpy (path, dirname, sizeof (path));
  FileBuildPath (path, subname, NULL);
  FileBuildPath (path, NULL, filename);
  return (Boolean) (FileLength (path) > 0);
}

static Boolean  CheckAsnloadPath (CharPtr dirname, CharPtr subdir)
{

#ifdef ASNLOAD_NEEDED
  Char            fname[16];
  int             i;

  for (i = 60; i <= 69; ++i) {
    sprintf (fname, "asnmedli.l%02d", (int) i);
    if (FileExists (dirname, subdir, fname)) {
      return TRUE;
    }
  }
  return FALSE;
#else
  return TRUE;
#endif
}

static Boolean  CheckDataPath (CharPtr dirname, CharPtr subdir)
{
  return (Boolean) (FileExists (dirname, subdir, "seqcode.val"));
}

static Boolean  CheckErrMsgPath (CharPtr dirname, CharPtr subdir)
{
  return (Boolean) (FileExists (dirname, subdir, "valid.msg"));
}

static void     SetPathText (CharPtr dirname, CharPtr subname, TexT t)
{
  Char            path[PATH_MAX];

  StringNCpy (path, dirname, sizeof (path));
  FileBuildPath (path, subname, NULL);
  SafeSetTitle (t, path);
}

#ifdef ASNLOAD_NEEDED
static void     FindAsnloadPath (ButtoN b)
{
  NetConfigPtr    ncp;
  Char            path[PATH_MAX];
  CharPtr         ptr;

  ncp = (NetConfigPtr) GetObjectExtra (b);
  if (ncp == NULL)
    return;
  if (GetInputFileName (path, sizeof (path), "", "")) {
    ptr = StringRChr (path, DIRDELIMCHR);
    if (ptr != NULL) {
      ptr++;
      *ptr = '\0';
      SetTitle (ncp->asnloadPath, path);
    }
  }
}
#endif

static void     FindDataPath (ButtoN b)
{
  NetConfigPtr    ncp;
  Char            path[PATH_MAX];
  CharPtr         ptr;

  ncp = (NetConfigPtr) GetObjectExtra (b);
  if (ncp == NULL)
    return;
  if (GetInputFileName (path, sizeof (path), "", "")) {
    ptr = StringRChr (path, DIRDELIMCHR);
    if (ptr != NULL) {
      ptr++;
      *ptr = '\0';
      SetTitle (ncp->dataPath, path);
    }
  }
}

static void     FindErrMsgPath (ButtoN b)
{
  NetConfigPtr    ncp;
  Char            path[PATH_MAX];
  CharPtr         ptr;

  ncp = (NetConfigPtr) GetObjectExtra (b);
  if (ncp == NULL)
    return;
  if (GetInputFileName (path, sizeof (path), "", "")) {
    ptr = StringRChr (path, DIRDELIMCHR);
    if (ptr != NULL) {
      ptr++;
      *ptr = '\0';
      SetTitle (ncp->errmsgPath, path);
    }
  }
}

static void     GetDispName (NetConfigPtr ncp, CharPtr buf, size_t bufsiz)
{
  Int2            val;
  ValNodePtr      vnp;

  val = GetValue (ncp->dispList);
  if (val > 0 && val < (NUMDISP + 1)) {
    vnp = ncp->dispTitles;
    while (val > 1) {
      val--;
      vnp = vnp->next;
    }
    if (vnp != NULL) {
      StringNCpy (buf, vnp->data.ptrvalue, bufsiz);
    }
  } else {
    GetTitle (ncp->dispAddr, buf, bufsiz);
  }
}

static Boolean  EstabAnonConnection (NetConfigPtr ncp)
{
  char            buf[64];

  GetDispName (ncp, buf, sizeof (buf));
  ncp->dispatcher = NI_SetDispatcher (NULL, buf, NULL, 0, 0, NULL, FALSE);
  GetTitle (ncp->username, buf, sizeof (buf));
  return NI_InitServices (ncp->dispatcher, buf, "ANONYMOUS", NULL, &ncp->dispInfoPtr) >= 0;
}

static CharPtr  page1Msg = "\
Entrez needs to use some system files, which are \
normally located in folders called \"asnload\" and \
\"data\" within the application program folder.\n";

static void     CreatePage1 (GrouP h, NetConfigPtr ncp)
{
  Boolean         asnFound;
  ButtoN          b;
  Boolean         dataFound;
  Boolean         errmsgFound;
  GrouP           g;
  GrouP           j;
  GrouP           k;
  Char            path[PATH_MAX];
  CharPtr         ptr;

  g = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (g, 3, 10);
  StaticPrompt (g, "NCBI System Files", 0, 0, systemFont, 'l');
  MultiLinePrompt (g, page1Msg, stdCharWidth * 28, programFont);

#ifdef ASNLOAD_NEEDED
  j = HiddenGroup (g, -1, 0, NULL);
  StaticPrompt (j, "Folder for ASN.1 parser files:", 0, 0, programFont, 'l');
  k = HiddenGroup (j, 2, 0, NULL);
  ncp->asnloadPath = DialogText (k, "", 20, NULL);
  SetObjectExtra (ncp->asnloadPath, ncp, NULL);
  b = PushButton (k, "Browse", FindAsnloadPath);
  SetObjectExtra (b, ncp, NULL);
#endif

  j = HiddenGroup (g, -1, 0, NULL);
  StaticPrompt (j, "Folder for common data files:", 0, 0, programFont, 'l');
  k = HiddenGroup (j, 2, 0, NULL);
  ncp->dataPath = DialogText (k, "", 20, NULL);
  SetObjectExtra (ncp->dataPath, ncp, NULL);
  b = PushButton (k, "Browse", FindDataPath);
  SetObjectExtra (b, ncp, NULL);

  j = HiddenGroup (g, -1, 0, NULL);
  StaticPrompt (j, "Folder for error message files:", 0, 0, programFont, 'l');
  k = HiddenGroup (j, 2, 0, NULL);
  ncp->errmsgPath = DialogText (k, "", 20, NULL);
  SetObjectExtra (ncp->errmsgPath, ncp, NULL);
  b = PushButton (k, "Browse", FindErrMsgPath);
  SetObjectExtra (b, ncp, NULL);

  ncp->page1ok = FALSE;
  asnFound = FALSE;
  dataFound = FALSE;

  GetAppParam ("NCBI", "NCBI", "ASNLOAD", NULL, path, sizeof (path));
  asnFound = CheckAsnloadPath (path, NULL);
  if (asnFound) {
    SetTitle (ncp->asnloadPath, path);
  }
  GetAppParam ("NCBI", "NCBI", "DATA", NULL, path, sizeof (path));
  dataFound = CheckDataPath (path, NULL);
  if (dataFound) {
    SetTitle (ncp->dataPath, path);
  }
  GetAppParam ("NCBI", "ErrorProcessing", "MsgPath", NULL, path, sizeof (path));
  errmsgFound = CheckErrMsgPath (path, NULL);
  if (errmsgFound) {
    SetTitle (ncp->errmsgPath, path);
  }
  if (asnFound && dataFound) {
    ncp->page1ok = TRUE;
  }
  if (TextHasNoText (ncp->asnloadPath) && TextHasNoText (ncp->dataPath)) {
    ProgramPath (path, sizeof (path));
    ptr = StringRChr (path, DIRDELIMCHR);
    if (ptr != NULL) {
      ptr++;
      *ptr = '\0';
    }
    asnFound = CheckAsnloadPath (path, "asnload");
    dataFound = CheckDataPath (path, "data");
    if (asnFound && dataFound) {
      SetPathText (path, "asnload", ncp->asnloadPath);
      SetPathText (path, "data", ncp->dataPath);
      ncp->page1ok = TRUE;
      errmsgFound = CheckErrMsgPath (path, "errmsg");
      if (errmsgFound) {
        SetPathText (path, "errmsg", ncp->errmsgPath);
      }
    }
  }
}

static Boolean  LeavePage1 (NetConfigPtr ncp)
{
  Char            path[PATH_MAX];

#ifdef ASNLOAD_NEEDED
  GetTitle (ncp->asnloadPath, path, sizeof (path));
  if (StringHasNoText (path)) {
    MsgAlert (KEY_OK, SEV_ERROR, IDENTIFYING_STRING,
              "You must enter a valid path for "
              "the ASNLOAD setting in order to continue");
    return FALSE;
  }
  if (!CheckAsnloadPath (path, NULL)) {
    MsgAlert (KEY_OK, SEV_ERROR, IDENTIFYING_STRING,
              "The path you have entered for the ASNLOAD setting\n(%s)\n"
            "does not contain the expected files. Please try again.", path);
    return FALSE;
  }
#endif

  GetTitle (ncp->dataPath, path, sizeof (path));
  if (StringHasNoText (path)) {
    MsgAlert (KEY_OK, SEV_ERROR, IDENTIFYING_STRING,
              "You must enter a valid path for "
              "the DATA setting in order to continue");
    return FALSE;
  }
  if (!CheckDataPath (path, NULL)) {
    MsgAlert (KEY_OK, SEV_ERROR, IDENTIFYING_STRING,
              "The path you have entered for the DATA setting\n(%s)\n"
           "does not contain the expected files.  Please try again.", path);
    return FALSE;
  }

#ifdef ASNLOAD_NEEDED
  GetTitle (ncp->asnloadPath, path, sizeof (path));
  TransientSetAppParam ("NCBI", "NCBI", "ASNLOAD", path);
#endif

  GetTitle (ncp->dataPath, path, sizeof (path));
  TransientSetAppParam ("NCBI", "NCBI", "DATA", path);
  return TRUE;
}

static CharPtr  page2Msg = "\
Enter a username to identify this user or this site; this \
information is not used for validation, but helps to track \
usage.  Also indicate what action your client software should \
perform when it is unable to contact the primary Dispatcher \
computer.\n";

static void     CreatePage2 (GrouP h, NetConfigPtr ncp)
{
  Char            buf[64];
  GrouP           g;
  GrouP           j;
  GrouP           k;

  g = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (g, 3, 10);
  StaticPrompt (g, "Internet Setup Preferences", 0, 0, systemFont, 'l');
  MultiLinePrompt (g, page2Msg, stdCharWidth * 28, programFont);

  j = HiddenGroup (g, -1, 0, NULL);
  k = HiddenGroup (j, 2, 0, NULL);
  StaticPrompt (k, "Username:", 0, dialogTextHeight, programFont, 'l');
  ncp->username = DialogText (k, "", 20, NULL);
  SetObjectExtra (ncp->username, ncp, NULL);

  GetAppParam ("NCBI", "NET_SERV", "DISP_USERNAME", NULL, buf, sizeof (buf));

#ifdef OS_UNIX
  if (buf[0] == '\0') {
    StringCpy (buf, getlogin ());
  }
#endif

#ifdef OS_VMS
  if (buf[0] == '\0') {
    StringCpy (buf, getenv ("USER"));
  }
#endif

  SetTitle (ncp->username, buf);

  ncp->ifCantFindDisp = NormalGroup (g, 0, 5,
                                     "When cannot find primary dispatcher",
                                     programFont, NULL);
  SetObjectExtra (ncp->ifCantFindDisp, ncp, NULL);
  RadioButton (ncp->ifCantFindDisp, "Try alternate dispatchers");
  RadioButton (ncp->ifCantFindDisp, "Ask the user what to do");
  RadioButton (ncp->ifCantFindDisp, "Give up");

  GetAppParam ("NCBI", "NET_SERV", "DISP_RECONN_ACTION", "CONT", buf, sizeof (buf));
  if (StringICmp (buf, "CONT") == 0) {
    SetValue (ncp->ifCantFindDisp, 1);
  } else if (StringICmp (buf, "ASK") == 0) {
    SetValue (ncp->ifCantFindDisp, 2);
  } else if (StringICmp (buf, "QUIT") == 0) {
    SetValue (ncp->ifCantFindDisp, 3);
  } else {
    SetValue (ncp->ifCantFindDisp, 1);
  }

  ncp->outgoingOnly = CheckBox (g, "Outgoing connections only", NULL);
  GetAppParam ("NCBI", "NET_SERV", "DIRECT_SVC_CON", "FALSE", buf, sizeof (buf));
  if (StringICmp (buf, "TRUE") == 0) {
    SetStatus (ncp->outgoingOnly, TRUE);
  }
  if (NI_EncrAvailable ()) {
    ncp->wantEncryption = CheckBox (g, "Use encryption", NULL);
    GetAppParam ("NCBI", "NET_SERV", "ENCRYPTION_DESIRED", "FALSE", buf, sizeof (buf));
    if (StringICmp (buf, "TRUE") == 0) {
      SetStatus (ncp->wantEncryption, TRUE);
    }
  }
}

static Boolean  LeavePage2 (NetConfigPtr ncp)
{
  return TRUE;
}

static CharPtr  page3Msg = "\
The client software needs to know either the fully qualified \
domain name (FQDN) or the dotted Internet address of an NCBI \
dispatcher.  Please select one of the possibilities below \
by clicking on that address, or enter an address manually.\n";


static Boolean  ResourceCompatWService (NIResPtr resp, NISvcPtr service)
{
  NodePtr         np = service->typeL;

  if (np == NULL)
    return FALSE;

  do {
    np = np->next;
    if (StrICmp (np->elem, MATCHES_ANY_TYPE) == 0 ||
        StrICmp (np->elem, resp->type) == 0)
      return TRUE;
  } while (np != NULL && np != service->typeL);

  return FALSE;
}

static void     SetService (NetConfigPtr ncp, Int2 n)
{
  ncp->rsc = NULL;

  if (n > 0 && n <= ncp->svc_ct) {
    NIToolset      *tset = ncp->toolset[n];

    NIService      *svc = ncp->svc[n];

    NIResource     *rsc;

    Node           *rn;

    /* find the "best" resource for this service */
    for (rn = tset->resources; rn; rn = rn->next) {
      rsc = (NIResource *) rn->elem;
      if (!ResourceCompatWService (rsc, svc))
        continue;
      if (ncp->rsc == NULL)
        ncp->rsc = rsc;
      else if (ncp->rsc->maxVersion != 0) {
        if (rsc->maxVersion == 0)
          ncp->rsc = rsc;
      }
    }

    SetTitle (ncp->servDesc, svc->descrip);
    SetTitle (ncp->datasetDesc, ncp->rsc->descrip);
  } else {
    n = -1;
    SetTitle (ncp->servDesc, NULL);
    SetTitle (ncp->datasetDesc, NULL);
  }
}

static void     ServListProc (LisT l)
{
  NetConfigPtr    ncp;

  ncp = (NetConfigPtr) GetObjectExtra (l);
  if (ncp == NULL)
    return;
  SetService (ncp, GetValue (l));
}

static void     DispListProc (GrouP g)
{
  NetConfigPtr    ncp;

  ncp = (NetConfigPtr) GetObjectExtra (g);
  if (ncp == NULL)
    return;
  if (GetValue (g) > NUMDISP) {
    SafeShow (ncp->dispAddr);
  } else {
    SafeHide (ncp->dispAddr);
  }
}

static void     CreatePage3 (GrouP h, NetConfigPtr ncp)
{
  Char            buf[64];
  GrouP           g;
  Int2            i;
  ValNodePtr      vnp;

  g = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (g, 3, 10);
  StaticPrompt (g, "Dispatcher Internet Address", 0, 0, systemFont, 'l');
  MultiLinePrompt (g, page3Msg, stdCharWidth * 28, programFont);

  ncp->dispList = NormalGroup (g, -2, 0,
                               "Primary dispatcher address",
                               programFont, DispListProc);
  SetGroupSpacing (ncp->dispList, 10, 2);
  SetObjectExtra (ncp->dispList, ncp, NULL);

  for (i = 0; i < (Int2) DIM (dispAddrStrs); i++) {
    StringNCpy (buf, dispAddrStrs[i].fqdn, sizeof (buf) - 1);
    RadioButton (ncp->dispList, buf);
    ValNodeCopyStr (&ncp->dispTitles, 0, buf);
    if (NI_FqdnToIpaddr((CharPtr) dispAddrStrs[i].fqdn, buf, sizeof(buf))) {
      RadioButton (ncp->dispList, buf);
      ValNodeCopyStr (&ncp->dispTitles, 0, buf);
    } else {
      StringNCpy (buf, dispAddrStrs[i].addr, sizeof (buf) - 1);
      RadioButton (ncp->dispList, buf);
      ValNodeCopyStr (&ncp->dispTitles, 0, buf);
    }
  }
  RadioButton (ncp->dispList, "Other:");
  ncp->dispAddr = DialogText (ncp->dispList, "", 15, NULL);
  SetObjectExtra (ncp->dispAddr, ncp, NULL);

  GetAppParam ("NCBI", "NET_SERV", "DISPATCHER", "", buf, sizeof (buf));
  i = 1;
  vnp = ncp->dispTitles;
  while (vnp != NULL && StringICmp (buf, vnp->data.ptrvalue) != 0) {
    i++;
    vnp = vnp->next;
  }
  if (vnp != NULL) {
    SetValue (ncp->dispList, i);
    Hide (ncp->dispAddr);
  } else {
    if (StringCmp (buf, "") == 0) {
      SetValue (ncp->dispList, 2);
      Hide (ncp->dispAddr);
    } else {
      SetValue (ncp->dispList, i);
      SetTitle (ncp->dispAddr, buf);
      Show (ncp->dispAddr);
    }
  }

  ncp->testDisp = CheckBox (g, "Test connection during configuration", NULL);
  if (GetAppParamBoolean ("NCBI", "Net-Disp", "Connect", TRUE)) {
    SetStatus (ncp->testDisp, TRUE);
  }
}

static Boolean  LeavePage3 (NetConfigPtr ncp)
{
  Char            buf[64];

  if (! ncp->movingForward)
      return TRUE;

  buf[0] = '\0';
  GetDispName (ncp, buf, sizeof (buf));
  if (StringHasNoText (buf)) {
    MsgAlert (KEY_OK, SEV_ERROR, IDENTIFYING_STRING,
              "No Dispatcher has been specified");
    return FALSE;
  }
  if (GetStatus (ncp->testDisp)) {
    if (!EstabAnonConnection (ncp)) {
      MsgAlert (KEY_OK, SEV_ERROR, IDENTIFYING_STRING,
                "Unable to connect to dispatcher");
      NI_EndServices (ncp->dispatcher);
      ncp->dispatcher = NULL;
      return FALSE;
    }
  }
  return TRUE;
}

static CharPtr  page4Msg = "\
Please select the desired Entrez service from the list below.\n";

static void     CreatePage4 (GrouP h, NetConfigPtr ncp)
{
  GrouP           g;
  GrouP           j;
  GrouP           k;

  g = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (g, 3, 10);
  StaticPrompt (g, "Entrez Service Selection", 0, 0, systemFont, 'l');
  MultiLinePrompt (g, page4Msg, stdCharWidth * 28, programFont);

  j = HiddenGroup (g, -1, 0, NULL);
  ncp->servList = SingleList (j, 28, 3, ServListProc);
  SetObjectExtra (ncp->servList, ncp, NULL);

  k = NormalGroup (g, -1, 0, "Service Description", programFont, NULL);
  ncp->servDesc = StaticPrompt (k, "", stdCharWidth * 25, 0, programFont, 'l');

  k = NormalGroup (g, -1, 0, "Default Dataset Description", programFont, NULL);
  ncp->datasetDesc = StaticPrompt (k, "", stdCharWidth * 25, 0, programFont, 'l');

  j = HiddenGroup (g, -1, 0, NULL);
  ncp->testServ = CheckBox (j, "Test connection during configuration", NULL);
#ifndef WIN32
  if (GetAppParamBoolean ("NCBI", "Net-Serv", "Connect", TRUE)) {
    SetStatus (ncp->testServ, TRUE);
  }
#else
 SetStatus (ncp->testServ, FALSE);
#endif
/* the above ifndef works around the problem on windows where the initial select()
on the server fails.  lyg */

}

extern Int2     NI_DestroyMsgCatalog (NICatalogPtr cp);

static Int2     PopulateServiceList (NetConfigPtr ncp)
{
  Char            buf[64];
  Int2            localSelService;
  Node           *sn;
  NIService      *svc;
  Int2            svc_ct;
  Node           *tn;
  NIToolset      *tset;

  svc_ct = 0;
  GetAppParam ("NCBI", "ENTREZ_NET", "SERVICE_NAME", "", buf, sizeof (buf));
  if (ncp->catalog != NULL) {
    localSelService = 0;
    for (tn = ncp->catalog->toolsetL; tn != NULL; tn = tn->next) {
      tset = (NIToolset *) tn->elem;
      if (tset != NULL) {
        for (sn = tset->services; sn != NULL; sn = sn->next) {
          svc = (NIService *) sn->elem;
          if (svc != NULL && svc->typeL &&
              StringCmp ((CharPtr) svc->typeL->elem, SERVICE_TYPE) == 0) {
            svc_ct++;
            ncp->toolset[svc_ct] = tset;
            ncp->svc[svc_ct] = svc;
            if (svc->name != NULL) {
              ListItem (ncp->servList, svc->name);
              ValNodeCopyStr (&ncp->servTitles, 0, svc->name);
              if (buf[0] != '\0' && StringCmp (svc->name, buf) == 0) {
                localSelService = svc_ct;
              }
            }
          }
        }
      }
    }
    if (localSelService > 0) {
      SetValue (ncp->servList, localSelService);
      ncp->svc_ct = svc_ct;
      SetService (ncp, localSelService);
    }
  }
  return svc_ct;
}

static NICatalog *
                Catalog_New (const char *motd)
{
  NICatalog      *cat = (NICatalog *) MemNew (sizeof (NICatalog));

  if (cat != NULL) {
    cat->motd = motd ? StrSave (motd) : NULL;
  }
  return cat;
}

static void     Catalog_AddToolset (NICatalog * cat, NIToolset * tset)
{
  Node           *node = (Node *) MemNew (sizeof (Node));

  if (node != NULL) {
    ASSERT (cat != NULL);
    ASSERT (tset != NULL);

    node->next = cat->toolsetL;
    if (node->next)
      node->next->last = node;
    cat->toolsetL = node;
    node->elem = (void *) tset;
  }
}

static NIToolset *
                Toolset_New (const char *name, int minv, int maxv, const char *descr)
{
  NIToolset      *tset = (NIToolset *) MemNew (sizeof (NIToolset));

  if (tset != NULL) {
    NIService      *svc = MemNew (sizeof (NIService));

    if (svc == NULL) {
      MemFree ((void *) tset);
      tset = NULL;
    } else {
      svc->name = name ? StrSave (name) : NULL;
      svc->minVersion = (Uint2) minv;
      svc->maxVersion = (Uint2) maxv;
      svc->descrip = descr ? StrSave (descr) : NULL;
      svc->typeL = (Node *) MemNew (sizeof (Node));
      svc->typeL->elem = (void *) StrSave (SERVICE_TYPE);
      svc->typeL->next = svc->typeL;        /* this is a single-element circular
                                         * list */
      svc->typeL->last = svc->typeL;        /* this is a single-element circular
                                         * list */
    }
    tset->services = (Node *) MemNew (sizeof (Node));
    tset->services->elem = (void *) svc;
  }
  return tset;
}

static void     Toolset_AddResource (NIToolset * tset, NIResource * rsrc)
{
  Node           *node = (Node *) MemNew (sizeof (Node));

  if (node != NULL) {
    ASSERT (tset != NULL);
    ASSERT (rsrc != NULL);

    node->next = tset->resources;
    if (node->next)
      node->next->last = node;
    tset->resources = node;
    node->elem = (void *) rsrc;
  }
}

static NIResource *
                Resource_New (const char *name, int minv, int maxv, const char *descr)
{
  NIResource     *rsrc = (NIResource *) MemNew (sizeof (NIResource));

  if (rsrc != NULL) {
    rsrc->name = name ? StrSave (name) : NULL;
    rsrc->minVersion = (Uint2) minv;
    rsrc->maxVersion = (Uint2) maxv;
    rsrc->descrip = descr ? StrSave (descr) : NULL;
    rsrc->type = StringSave (SERVICE_TYPE);
  }
  return rsrc;
}

static NICatalog *
                FabricateCatalog (void)
{
  NICatalog      *cat = Catalog_New (NULL);

  ErrLogPrintf ("Using fabricated service catalog\n");
  if (cat != NULL) {
    NIToolset      *tset1 = Toolset_New ("Entrez", 1, 0, "Network Entrez");

    Toolset_AddResource (tset1, Resource_New ("Entrez", 130, 130, "Entrez 13.0"));
    Toolset_AddResource (tset1, Resource_New ("Entrez", 140, 140, "Entrez 14.0"));
    Toolset_AddResource (tset1, Resource_New ("Entrez", 1, 0, "Most recent Entrez release"));
    Catalog_AddToolset (cat, tset1);
  }
  return cat;
}


static void     EnterPage3 (NetConfigPtr ncp)
{
  if (!ncp->movingForward && GetStatus (ncp->testDisp)) {
      NI_EndServices (ncp->dispatcher);
      ncp->dispatcher = NULL;
  }
}


static void     EnterPage4 (NetConfigPtr ncp)
{
  Reset (ncp->servList);
  if (!ncp->movingForward && GetStatus (ncp->testDisp)) {
    if (!EstabAnonConnection (ncp)) {
      MsgAlert (KEY_OK, SEV_ERROR, IDENTIFYING_STRING, "Unable to connect to dispatcher");
      NI_EndServices (ncp->dispatcher);
      ncp->dispatcher = NULL;
      return;
    }
  }
  if (ncp->dispatcher != NULL) {
    ncp->catalog = NI_GetCatalog (ncp->dispatcher);
    if (ncp->catalog != NULL) {
      Enable (ncp->testServ);
    } else {
      MsgAlert (KEY_OK, SEV_ERROR, IDENTIFYING_STRING,
                "Unable to retrieve service catalog from the dispatcher");
      Disable (ncp->testServ);
    }
  }
  ncp->svc_ct = PopulateServiceList (ncp);
  if (ncp->svc_ct <= 0) {
    NI_DestroyMsgCatalog (ncp->catalog);
    ncp->catalog = NULL;
    /* fabricate and try again */
    ncp->catalog = FabricateCatalog ();
    ncp->svc_ct = PopulateServiceList (ncp);
    MsgAlert (KEY_OK, SEV_WARNING, IDENTIFYING_STRING, "Using a fabricated service catalog since no suitable services are currently available");

  }
}

static Boolean  LeavePage4 (NetConfigPtr ncp)
{
  char            buffer[64];
  NI_HandPtr      shSvc;
  NIService      *svc;
  Boolean         outgoingConns;


  ncp->svc_sel = GetValue (ncp->servList);

  if (ncp->movingForward) {
    if (ncp->svc_sel <= 0 || ncp->svc_sel > ncp->svc_ct) {
      MsgAlert (KEY_OK, SEV_ERROR, IDENTIFYING_STRING, "No service has been specified\n");
      return FALSE;
    }
    svc = ncp->svc[ncp->svc_sel];
    MemFree (ncp->svname);
    ncp->svname = StringSave (svc->name);
    TransientSetAppParam ("NCBI", "ENTREZ_NET", "SERVICE_NAME", ncp->svname);
    ncp->verSvcMin = svc->minVersion;
    ncp->verSvcMax = svc->maxVersion;
    MemFree (ncp->resname);
    ncp->resname = StringSave ("Entrez");
    ncp->verRscMin = 1;
    ncp->verRscMax = 0;
    if (ncp->rsc != NULL) {
      MemFree (ncp->resname);
      ncp->resname = StringSave (ncp->rsc->name);
      ncp->verRscMin = ncp->rsc->minVersion;
      ncp->verRscMax = ncp->rsc->maxVersion;
    }
    if (GetStatus (ncp->testServ)) {
      WatchCursor ();
      /* must disconnect as anonymous and reconnect as GUEST */
      NI_EndServices (ncp->dispatcher);
      GetDispName (ncp, buffer, sizeof (buffer));
      outgoingConns = GetStatus (ncp->outgoingOnly);
      ncp->dispatcher = NI_SetDispatcher (NULL, buffer, NULL, 0, 0, NULL,
                                          outgoingConns);
      GetTitle (ncp->username, buffer, sizeof (buffer));

      if (NI_InitServices (ncp->dispatcher, buffer, "GUEST", NULL, &ncp->dispInfoPtr) < 0) {
        ArrowCursor ();
        ErrShow ();
        NI_DestroyMsgCatalog (ncp->catalog);        /* JAE ??? */
        ncp->catalog = NULL;
        ncp->svc_ct = 0;
        NI_EndServices (ncp->dispatcher);
        ncp->dispatcher = NULL;
        EstabAnonConnection (ncp);
        return FALSE;
      }
      WatchCursor ();
      if (ncp->dispatcher->adminInfo != NULL && ncp->dispatcher->adminInfo[0] != NULLB) {
        ncp->adminInfo = StringSave (ncp->dispatcher->adminInfo);
      }
      if ((shSvc = NI_ServiceGet (ncp->dispatcher, ncp->svname, ncp->verSvcMin, ncp->verSvcMax, ncp->resname, SERVICE_TYPE, ncp->verRscMin, ncp->verRscMax)) == NULL) {
        ArrowCursor ();
        MsgAlert (KEY_OK, SEV_ERROR, IDENTIFYING_STRING, "Unable to get service [%s]: %s", ncp->svname, ni_errlist[ni_errno]);
        return FALSE;
      } else {
        NI_DestroyMsgCatalog (ncp->catalog);
        ncp->catalog = NULL;
        ncp->svc_ct = 0;
        NI_ServiceDisconnect (shSvc);
        NI_EndServices (ncp->dispatcher);
        ncp->dispatcher = NULL;
        ArrowCursor ();
        return TRUE;
      }
    } else {
      NI_DestroyMsgCatalog (ncp->catalog);
      ncp->catalog = NULL;
      ncp->svc_ct = 0;
      if (ncp->dispatcher) {
         NI_EndServices (ncp->dispatcher);
         ncp->dispatcher = NULL;
      }
      return TRUE;
    }
  } else {                        /* backward */
    NI_EndServices (ncp->dispatcher);
    ncp->dispatcher = NULL;
    NI_DestroyMsgCatalog (ncp->catalog);
    ncp->catalog = NULL;
    ncp->svc_ct = 0;
    return TRUE;
  }
}

static CharPtr  page5Msg = "Press Accept to save the settings.\n\n";

static void     CreatePage5 (GrouP h, NetConfigPtr ncp)
{
  GrouP           g;
  GrouP           j;

  g = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (g, 3, 10);
  StaticPrompt (g, "Configuration Complete!", 0, 0, systemFont, 'l');
  MultiLinePrompt (g, page5Msg, stdCharWidth * 28, programFont);

  j = NormalGroup (g, -1, 0, "For your information", programFont, NULL);
  ncp->fyi[0] = StaticPrompt (j, "", stdCharWidth * 25, 0, programFont, 'l');
  ncp->fyi[1] = StaticPrompt (j, "", stdCharWidth * 25, 0, programFont, 'l');
  ncp->fyi[2] = StaticPrompt (j, "", stdCharWidth * 25, 0, programFont, 'l');
}

static void     EnterPage5 (NetConfigPtr ncp)
{
  char            buffer[300];
  CharPtr         buf2;
  CharPtr         ptr;

  if (ncp->adminInfo != NULL) {
    SetTitle (ncp->fyi[0], "Your Network Entrez administrator is:");
    buf2 = StringSave (ncp->adminInfo);
    if ((ptr = StrStr (buf2, " Email:")) != NULL) {        /* break into two lines */
      *ptr = NULLB;
      sprintf (buffer, "  %s", buf2);
      SetTitle (ncp->fyi[1], buffer);
      sprintf (buffer, "        %s", ptr + 1);
      SetTitle (ncp->fyi[2], buffer);
    } else {
      sprintf (buffer, "  %s", buf2);
      SetTitle (ncp->fyi[1], buffer);
      SetTitle (ncp->fyi[2], "");
    }
    MemFree (buf2);
  } else {
    SetTitle (ncp->fyi[0], "Your Network Entrez administrator is unidentified");
    SetTitle (ncp->fyi[1], "");
    SetTitle (ncp->fyi[2], "");
  }
}


static void     EnterConfigPage (NetConfigPtr ncp, Int2 newpage, Int2 oldpage)
{
  if (ncp != NULL) {
    ncp->movingForward = newpage > oldpage;
    switch (newpage) {
      case FIRST_PAGE:
        SafeDisable (ncp->prevBtn);
        SafeSetTitle (ncp->nextBtn, "Next Page >>");

#ifdef ASNLOAD_NEEDED
        Select (ncp->asnloadPath);
#else
        Select (ncp->dataPath);
#endif

        break;
      case SECOND_PAGE:
        SafeEnable (ncp->prevBtn);
        SafeSetTitle (ncp->nextBtn, "Next Page >>");
        Select (ncp->username);
        break;
      case THIRD_PAGE:
        SafeEnable (ncp->prevBtn);
        SafeSetTitle (ncp->nextBtn, "Next Page >>");
        EnterPage3 (ncp);
        break;
      case FOURTH_PAGE:
        SafeEnable (ncp->prevBtn);
        SafeSetTitle (ncp->nextBtn, "Next Page >>");
        EnterPage4 (ncp);
        break;
      case FIFTH_PAGE:
        SafeEnable (ncp->prevBtn);
        SafeSetTitle (ncp->nextBtn, "Accept");
        EnterPage5 (ncp);
        break;
      default:
        break;
    }
  }
}

static void     ChangeConfigFormPage (VoidPtr data, Int2 newval, Int2 oldval)
{
  NetConfigPtr    ncp;

  ncp = (NetConfigPtr) data;
  if (ncp != NULL) {
    ncp->currentPage = newval;
    SafeHide (ncp->pages[oldval]);
    Update ();
    EnterConfigPage (ncp, newval, oldval);
    SafeShow (ncp->pages[newval]);
    Update ();
  }
}

static Boolean  TestConfigPage (NetConfigPtr ncp, Int2 newpage, Int2 oldpage)
{
  ncp->movingForward = newpage > oldpage;
  if (ncp != NULL) {
    switch (oldpage) {
      case FIRST_PAGE:
        return LeavePage1 (ncp);
      case SECOND_PAGE:
        return LeavePage2 (ncp);
      case THIRD_PAGE:
        return LeavePage3 (ncp);
      case FOURTH_PAGE:
        return LeavePage4 (ncp);
      default:
        break;
    }
  }
  return TRUE;
}

static void     SaveSettings (NetConfigPtr ncp)
{
  char            buffer[PATH_MAX];

  SetAppParamInt ("entrez", "Preferences", "MaxLoad", 1000);

#ifdef ASNLOAD_NEEDED
  GetTitle (ncp->asnloadPath, buffer, sizeof (buffer));
  SetAppParam ("NCBI", "NCBI", "ASNLOAD", buffer);
#endif

  GetTitle (ncp->dataPath, buffer, sizeof (buffer));
  SetAppParam ("NCBI", "NCBI", "DATA", buffer);

  GetTitle (ncp->errmsgPath, buffer, sizeof (buffer));
  if (CheckErrMsgPath (buffer, NULL)) {
    SetAppParam ("NCBI", "ErrorProcessing", "MsgPath", buffer);
    SetAppParam ("NCBI", "ErrorProcessing", "EO_BEEP", "No");
  }

  GetTitle (ncp->username, buffer, sizeof (buffer));
  SetAppParam ("NCBI", "NET_SERV", "DISP_USERNAME", buffer);

  if (NI_EncrAvailable ()) {
    SetAppParam ("NCBI", "NET_SERV", "ENCRYPTION_DESIRED", GetStatus (ncp->wantEncryption) ? "TRUE" : "FALSE");
  }
  SetAppParam ("NCBI", "NET_SERV", "DIRECT_SVC_CON", GetStatus (ncp->outgoingOnly) ? "TRUE" : "FALSE");

  switch (GetValue (ncp->ifCantFindDisp)) {
    case 2:
      StrCpy (buffer, "ASK");
      break;
    case 3:
      StrCpy (buffer, "QUIT");
      break;
    case 1:
    default:
      StrCpy (buffer, "CONT");
      break;
  }
  SetAppParam ("NCBI", "NET_SERV", "DISP_RECONN_ACTION", buffer);

  GetDispName (ncp, buffer, sizeof (buffer));
  if (ncp->dispInfoPtr == NULL)
    SetAppParam ("NCBI", "NET_SERV", "DISPATCHER", buffer);
  else
    NI_SetDispConfig(&ncp->dispInfoPtr, buffer, sizeof (buffer));

  SetAppParam ("NCBI", "NCBI", "MEDIA", "ENTREZ_NET");
  SetAppParam ("NCBI", "REFERENCE_FROM_NET", "MEDIA", "ENTREZ_NET");
  SetAppParam ("NCBI", "SEQUENCE_FROM_NET", "MEDIA", "ENTREZ_NET");
  SetAppParam ("NCBI", "LINKS_FROM_NET", "MEDIA", "ENTREZ_NET");
  SetAppParamInt ("NCBI", "LINKS_FROM_NET", "ENTR_SEQ__ENTR_SEQ", 1);
  SetAppParamInt ("NCBI", "LINKS_FROM_NET", "ENTR_REF__ENTR_SEQ", 1);
  SetAppParamInt ("NCBI", "LINKS_FROM_NET", "ENTR_SEQ__ENTR_REF", 1);
  SetAppParamInt ("NCBI", "LINKS_FROM_NET", "ENTR_REF__ENTR_REF", 1);
  SetAppParam ("NCBI", "ENTR_LINK", "CHANNELS", "LINKS_FROM_NET");
  SetAppParam ("NCBI", "ENTR_REF", "CHANNELS", "REFERENCE_FROM_NET");
  SetAppParam ("NCBI", "ENTR_SEQ", "CHANNELS", "SEQUENCE_FROM_NET");
  SetAppParam ("NCBI", "ENTREZ_NET", "TYPE", "NET");
  SetAppParam ("NCBI", "ENTREZ_NET", "SERVICE_NAME", ncp->svname);
  SetAppParam ("NCBI", "ENTREZ_NET", "RESOURCE_NAME", ncp->resname);
  SetAppParam ("NCBI", "ENTREZ_NET", "RESOURCE_TYPE", SERVICE_TYPE);
  SetAppParamInt ("NCBI", "ENTREZ_NET", "SERV_VERS_MIN", ncp->verSvcMin);
  SetAppParamInt ("NCBI", "ENTREZ_NET", "SERV_VERS_MAX", ncp->verSvcMax);
  SetAppParamInt ("NCBI", "ENTREZ_NET", "RES_VERS_MIN", ncp->verRscMin);
  SetAppParamInt ("NCBI", "ENTREZ_NET", "RES_VERS_MAX", ncp->verRscMax);
}

static void     NextConfigFormBtn (ButtoN b)
{
  VoidProc        accepted;
  NetConfigPtr    ncp;

  ncp = (NetConfigPtr) GetObjectExtra (b);
  if (ncp != NULL) {
    if (ncp->currentPage + 1 < NUM_PAGES) {
      if (TestConfigPage (ncp, ncp->currentPage + 1, ncp->currentPage)) {
        ChangeConfigFormPage (ncp, ncp->currentPage + 1, ncp->currentPage);
      }
    } else {
      accepted = ncp->accepted;
      Hide (ncp->form);
      SaveSettings (ncp);
      Remove (ncp->form);
      if (accepted != NULL) {
        accepted ();
      }
    }
  }
}

static void     PrevConfigFormBtn (ButtoN b)
{
  NetConfigPtr    ncp;

  ncp = (NetConfigPtr) GetObjectExtra (b);
  if (ncp != NULL) {
    if (ncp->currentPage > 0) {
      ChangeConfigFormPage (ncp, ncp->currentPage - 1, ncp->currentPage);
    }
  }
}

static void     CleanupConfigForm (GraphiC g, VoidPtr data)
{
  NetConfigPtr    ncp;
  Int2            sockCount;

  ncp = (NetConfigPtr) data;
  if (ncp != NULL) {
    if (ncp->dispatcher != NULL) {
      NI_EndServices (ncp->dispatcher);
      ncp->dispatcher = NULL;
    }
    if ((sockCount = NI_SocketsOpen ()) > 0) {
      TRACE ("At termination time, %d open sockets\n", (int) sockCount);
    }
    ncp->dispTitles = ValNodeFreeData (ncp->dispTitles);
    ncp->servTitles = ValNodeFreeData (ncp->servTitles);
    /* make sure we're disconnected from the net (user may have cancelled) */
  }
  StdCleanupFormProc (g, data);
}

extern void     ShowNetConfigForm (WndActnProc activate, FormMessageFunc messages,
                                      VoidProc accepted, VoidProc cancelled)
{
  ButtoN          b;
  GrouP           c;
  GrouP           g;
  GrouP           h;
  NetConfigPtr    ncp;
  WindoW          w;

  /* It is valid for the old-fashioned NCBI dispatcher only... */
  if ( !NI_IsInterfaceSupported(eNII_Dispatcher) ) {
    ErrPostEx(SEV_FATAL, 0, 0,
              "The eNII_Dispatcher interface is not supported.\n"
              "Please quit the configuration program and report.");
    return;
  }

  w = NULL;
  ncp = (NetConfigPtr) MemNew (sizeof (NetConfigData));
  if (ncp != NULL) {
    ncp->prev_interface = NI_SetInterface(eNII_Dispatcher);

    w = FixedWindow (-50, -33, -10, -10, "Network Configuration",
                     StdSendCloseWindowMessageProc);
    SetObjectExtra (w, ncp, CleanupConfigForm);
    ncp->form = (ForM) w;
    ncp->formmessage = ConfigMessageProc;

    ncp->appmessage = messages;
    ncp->activate = activate;
    SetActivate (w, ConfigFormActivate);

    ncp->accepted = accepted;
    ncp->cancelled = cancelled;

    ncp->currentPage = FIRST_PAGE;
    ncp->page1ok = FALSE;

    h = HiddenGroup (w, 0, 0, NULL);

    g = HiddenGroup (h, -1, 0, NULL);
    CreatePage1 (g, ncp);
    ncp->pages[FIRST_PAGE] = g;
    Hide (g);

    g = HiddenGroup (h, -1, 0, NULL);
    CreatePage2 (g, ncp);
    ncp->pages[SECOND_PAGE] = g;
    Hide (g);

    g = HiddenGroup (h, -1, 0, NULL);
    CreatePage3 (g, ncp);
    ncp->pages[THIRD_PAGE] = g;
    Hide (g);

    g = HiddenGroup (h, -1, 0, NULL);
    CreatePage4 (g, ncp);
    ncp->pages[FOURTH_PAGE] = g;
    Hide (g);

    g = HiddenGroup (h, -1, 0, NULL);
    CreatePage5 (g, ncp);
    ncp->pages[FIFTH_PAGE] = g;
    Hide (g);

    c = HiddenGroup (w, 4, 0, NULL);
    SetGroupSpacing (c, 10, 2);
    b = PushButton (c, "Cancel", StdSendCancelButtonMessageProc);
    SetObjectExtra (b, ncp, NULL);
    StaticPrompt (c, "", 30, 0, systemFont, 'l');
    ncp->prevBtn = PushButton (c, " << Prev Page ", PrevConfigFormBtn);
    SetObjectExtra (ncp->prevBtn, ncp, NULL);
    ncp->nextBtn = PushButton (c, " Next Page >> ", NextConfigFormBtn);
    SetObjectExtra (ncp->nextBtn, ncp, NULL);

    AlignObjects (ALIGN_CENTER,
                  (HANDLE) ncp->pages[FIRST_PAGE],
                  (HANDLE) ncp->pages[SECOND_PAGE],
                  (HANDLE) ncp->pages[THIRD_PAGE],
                  (HANDLE) ncp->pages[FOURTH_PAGE],
                  (HANDLE) ncp->pages[FIFTH_PAGE],
                  (HANDLE) c, NULL);

    RealizeWindow (w);

    Show (ncp->pages[ncp->currentPage]);
    EnterConfigPage (ncp, ncp->currentPage, ncp->currentPage);
    if (ncp->page1ok) {
      NextConfigFormBtn (ncp->nextBtn);
    }
    Show (w);
    Select (w);
  }
}
