/* $Id: cddfind.c,v 1.3 2000/07/21 21:35:15 bauer Exp $
*===========================================================================
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
* File Name:  cddfind.c
*
* Author:  Aron Marchler-Bauer
*
* Initial Version Creation Date: 10/19/1999
*
* $Revision: 1.3 $
*
* File Description: CDD locator/finder
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: cddfind.c,v $
* Revision 1.3  2000/07/21 21:35:15  bauer
* added verbose message to report that no hits are found
*
* Revision 1.2  2000/07/19 19:56:17  bauer
* added modification logging
*
*
* ==========================================================================
*/



#include <stdio.h>
#include <ncbi.h>
#include <accentr.h>
#include <lsqfetch.h>
#include <netentr.h>
#include <www.h>
#include <sys/resource.h>
#include <asn.h>
#include <accutils.h>
#include <mmdbapi.h>
#include <mmdbapi1.h>
#include <asnmime.h>
#include <objmime.h>
#include <strimprt.h>
#include "cddsrv.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* report Errors in processing and exit immediately                          */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
static void CddHtmlError(CharPtr cErrTxt) 
{
  printf("Content-type: text/html\n\n");
  printf("<h2>CDDfind Error:</h2>\n");
  printf("<h3>%s</h3>\n",cErrTxt);
  exit(1);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* read parameters from configuration file                                   */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
static Boolean CddGetParams()
{
  URLBase[0] = URLcgi[0] = ENTREZurl[0] = DOCSUMurl[0] = MAILto[0] = '\0';
  MMDBpath[0] = gunzip[0] = '\0';

  GetAppParam("cdd", "CDDSRV", "URLBase", "", URLBase, PATH_MAX);
  if (URLBase[0] == '\0') {
    ErrPostEx(SEV_FATAL,0,0,"CDD config file\nCDDSRV section has no URLBase...\n");
    return FALSE;
  }
  GetAppParam("cdd", "CDDSRV", "URLcgi", "", URLcgi, PATH_MAX);
  if (URLcgi[0] == '\0') {
                ErrPostEx(SEV_FATAL,0,0,"CDD config file\nCDDSRV section has no URLcgi...\n");
                return FALSE;
  }
  GetAppParam("cdd", "CDDSRV", "ENTREZurl", "", ENTREZurl, PATH_MAX);
  if (ENTREZurl[0] == '\0') {
    ErrPostEx(SEV_FATAL,0,0,"CDD config file\nCDDSRV section has no ENTREZurl...\n");
    return FALSE;
  }
  GetAppParam("cdd", "CDDSRV", "DOCSUMurl", "", DOCSUMurl, PATH_MAX);
  if (DOCSUMurl[0] == '\0') {
    ErrPostEx(SEV_FATAL,0,0,"CDD config file\nCDDSRV section has no DOCSUMurl...\n");
    return FALSE;
  }
  GetAppParam("cdd", "CDDSRV", "Gunzip", "", gunzip, (size_t) 256*(sizeof(Char)));
  if (gunzip[0] == '\0') {
    ErrPostEx(SEV_FATAL,0,0,"CDD config file\nCDDSRV section has no Gunzip...\n");
    return FALSE;
  }
  GetAppParam("cdd", "CDDSRV", "Database", "", MMDBpath, PATH_MAX);
  if (MMDBpath[0] == '\0') {
    ErrPostEx(SEV_FATAL,0,0,"MMDB config file\nMMDBSRV section has no Database...\n");
    return FALSE;
  }
  GetAppParam("cdd", "CDDSRV", "MAILto", "", MAILto, PATH_MAX);
  if (MAILto[0] == '\0') {
    ErrPostEx(SEV_FATAL,0,0,"CDD config file\nCDDSRV section has no MAILto...\n");
    return FALSE;
  }
  GetAppParam("cdd", "CDDSRV", "DATApath", "", DATApath, PATH_MAX);
  if (DATApath[0] == '\0') {
    ErrPostEx(SEV_FATAL,0,0,"CDD config file\nCDDSRV section has no VAST Data path...\n");
    return FALSE;
  }
  GetAppParam("cdd", "CDDSRV", "CDDatabase", "", CDDdpath, PATH_MAX);
  if (CDDdpath[0] == '\0') {
    ErrPostEx(SEV_FATAL,0,0,"CDD config file\nCDDSRV section has no CDD data path...\n");
    return FALSE;
  }
  GetAppParam("cdd", "CDDSRV", "CVDatabase", "", CDDvpath, PATH_MAX);
  if (CDDvpath[0] == '\0') {
    ErrPostEx(SEV_FATAL,0,0,"CDD config file\nCDDSRV section has no CDD/VAST data path...\n");
    return FALSE;
  }
  GetAppParam("cdd", "CDDSRV", "CDDextens", "", CDDextens, PATH_MAX);
  if (CDDextens[0] == '\0') {
    ErrPostEx(SEV_FATAL,0,0,"CDD config file\nCDDSRV section has no CDD file name extension...\n");
    return FALSE;
  }
  GetAppParam("cdd", "CDDSRV", "RAWextens", "", RAWextens, PATH_MAX);
  if (RAWextens[0] == '\0') {
    ErrPostEx(SEV_FATAL,0,0,"CDD config file\nCDDSRV section has no RAW file name extension...\n");
    return FALSE;
  }
  GetAppParam("cdd", "CDDSRV", "CDDdescr", "", CDDdescr, PATH_MAX);
  if (CDDdescr[0] == '\0') {
    ErrPostEx(SEV_FATAL,0,0,"CDD config file\nCDDSRV section has no description file name extension...\n");
    return FALSE;
  }
  GetAppParam("cdd", "CDDSRV", "Database", "", database, PATH_MAX);
  if (database[0] == '\0') {
    ErrPostEx(SEV_FATAL,0,0,"CDD config file\nCDDSRV section has no VAST data path...\n");
    return FALSE;
  }
  return TRUE;
}                                                       /* end GetVastParams */

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* allocate a new CddDesc linked list entry                                  */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
static CddDescPtr CddDescNew()
{
  CddDescPtr pcdd;

  pcdd = (CddDescPtr)MemNew(sizeof(CddDesc));
  if (pcdd == NULL) return pcdd;
  pcdd->cCddId[0] = '\0';
  pcdd->cDescr[0] = '\0';
  pcdd->cSourc[0] = '\0';
  pcdd->next = NULL;
  return(pcdd);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* Free a CddDesc linked list                                                */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
static CddDescPtr CddDescFree(CddDescPtr pcdd)
{
  CddDescPtr    next;
  
  while (pcdd) {
    next = pcdd->next;
    Nlm_MemFree(pcdd);
    pcdd = next;
  }
  return NULL;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* adds a to a linked list of CddDescPtr, always returns the beginning of the*/
/* list and always adds to the end of the list!!                             */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
static CddDescPtr CddDescLink(CddDescPtr PNTR head, CddDescPtr newnode)
{
  CddDescPtr     pcdd;
 
  if (head == NULL) return newnode;
  pcdd = *head;
  if (pcdd != NULL) {
    while(pcdd->next != NULL) pcdd = pcdd->next;
    pcdd->next = newnode;
  } else *head = newnode;
  return *head;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* Read the table of Cdd Names and descriptions                              */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
static CddDescPtr CddReadDescr() {
  FILE             *fp;
  Char              pcBuf[2048];
  CharPtr           pcTest;
  CddDescPtr        pcdd = NULL;
  CddDescPtr        pcddThis;

  if (!(fp = FileOpen(CDDdescr, READ))) 
    CddHtmlError("Can not read description file!");
  do {
    pcBuf[0]='\0';
    pcTest = fgets(pcBuf, (size_t)2048, fp);
    if (pcTest) if (pcTest[0] != ' ') {
      pcddThis = CddDescNew();
      strcpy(pcddThis->cCddId,strtok(pcTest,"\t"));
      strcpy(pcddThis->cDescr,strtok(NULL,"\t"));
      strcpy(pcddThis->cSourc,strtok(NULL,"\t"));
      CddDescLink(&(pcdd),pcddThis);
    }
  } while (pcTest);
  FileClose(fp);
  return(pcdd);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* convert strings into upper case                                           */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
static void CddUpCase(CharPtr p)
{
  Int4     i = 0;

  while (p[i] != '\0') {
    if (p[i]> 96 && p[i]< 123) p[i]-= 32;
    i++;
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* For faster (?) string comparisons                                         */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
static CharPtr BoyerMooreSearch(CharPtr txt, CharPtr pat)
{
  int    ch;
  Int2   d[256];
  Int2   i, j, k;
  size_t N, M;
  
  if (txt == NULL || pat == NULL) return NULL;
  N = StringLen(txt); M = StringLen(pat); if (N < M) return NULL;
  for (ch=0;ch<256;ch++) d[ch]=M;
  for (j=0;j<M-1;j++) {
    ch = (int) pat[j];
    d[ch] = M-j-1;
  }
  for (i=0;i<=N-M;ch=(int)txt[i+M-1],i+=d[ch]) {
    for (j=M-1,k=i+M-1;j>=0 && txt[k]==pat[j];j--,k--) continue;
    if (j < 0) return (txt+i);
  }
  return NULL;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* MAIN Function for cddfind                                                 */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Int2 Main()
{
  CddDescPtr               pcdd = NULL, pcddThis;
  CharPtr                  www_arg;
  WWWInfoPtr               www_info;
  Boolean                  bExecMode = FALSE;
  Int4                     NumLabels = 0, indx;
  Int4                     nhits     = 0;
  struct rlimit            rl;
  Char                     cDescr[CDD_MAX_DESCR];
  Char                     cSourc[PATH_MAX];
  Char                     cCddId[PATH_MAX];
  Char                     cSrch[PATH_MAX];
  Char                     cAll[CDD_MAX_DESCR+PATH_MAX];
  Boolean                  bFound;

/*---------------------------------------------------------------------------*/
/* this sets up the unix time limit                                          */
/*---------------------------------------------------------------------------*/
  getrlimit(RLIMIT_CPU, &rl);
  rl.rlim_max = rl.rlim_cur = CPUTIME_MAX;
  setrlimit(RLIMIT_CPU, &rl);

/*---------------------------------------------------------------------------*/
/* Begin processing www information block                                    */
/*---------------------------------------------------------------------------*/
  if (WWWGetArgs(&www_info) != WWWErrOk)
    CddHtmlError("Failed to process posting - check your get/post syntax.");    
  if ((NumLabels = WWWGetNumEntries(www_info)) == 0) { 
    printf("Content-type: text/html\n\n");
	  printf("<HTML><TITLE>CD finder</title>\n");
    printf("<BODY BGCOLOR=\"#FFFFFF\">\n");
    printf("<A HREF=\"cdd_form.map\">\n");
    printf("<IMG SRC=\"cdbrowse.gif\" BORDER=0 ISMAP></A>\n");
    printf("<FORM METHOD=\"POST\" ACTION=\"%scddfind.cgi\">\n",URLBase);
    printf("<INPUT TYPE=\"SUBMIT\" VALUE=\"Search CDD\"> for \n");
    printf("<INPUT TYPE=\"TEXT\" NAME=\"query\" SIZE=\"40\">\n");
    printf("<INPUT TYPE=\"RESET\">\n");
    printf("</FORM>\n");
    printf("</BODY></HTML>\n");
    exit(0);
  }

/*---------------------------------------------------------------------------*/
/* retrieve the query string                                                 */
/*---------------------------------------------------------------------------*/
  if ((indx = WWWFindName(www_info, "query")) >= 0) {
    www_arg =  WWWGetValueByIndex(www_info, indx);
    strcpy(cSrch,www_arg);
    CddUpCase(cSrch);
    /* printf("Query: %s\n",www_arg); */
  }
  if (www_arg == NULL) CddHtmlError("Can't interpret query string");

/*---------------------------------------------------------------------------*/
/* retrieve names for directories etc.                                       */
/*---------------------------------------------------------------------------*/
  if (!CddGetParams()) CddHtmlError("Couldn't read from config file...");
  pcdd = CddReadDescr();
  if (!pcdd) CddHtmlError("Could not read CDD descriptions!");
  
  printf("Content-type: text/html\n\n");
  printf("<HTML><TITLE>CD finder</title>\n");
  printf("<BODY BGCOLOR=\"#FFFFFF\">\n");
  printf("<A HREF=\"cdd_form.map\">\n");
  printf("<IMG SRC=\"cdbrowse.gif\" BORDER=0 ISMAP></A>\n");
  printf("<FORM METHOD=\"POST\" ACTION=\"%scddfind.cgi\">\n",URLBase);
  printf("<INPUT TYPE=\"SUBMIT\" VALUE=\"Search CDD\"> for \n");
  printf("<INPUT TYPE=\"TEXT\" NAME=\"query\" SIZE=\"40\">\n");
  printf("<INPUT TYPE=\"RESET\">\n");
  printf("</FORM>\n");
  printf("<BR>\n");
  printf("<p>Search results for: %s</p>\n",www_arg);
  printf("<TABLE BORDER=1 CELLPADDING=3 CELLSPACING=3>\n");
  if (strlen(cSrch) == 0 || strcmp(cSrch," ") == 0) {
    pcddThis = pcdd;
    while (pcddThis) {
      nhits++;
      printf("<TR>\n");
      printf("<TD NOWRAP VALIGN=TOP>\n");
      printf("<A HREF=\"%scddsrv.cgi?uid=%s\">\n",URLBase,pcddThis->cCddId);
      printf("%s</A>\n",pcddThis->cCddId);
      printf("</TD>\n");
      printf("<TD VALIGN=TOP>\n");
      printf("%s\n",pcddThis->cDescr);
      printf("</TD>\n");
      printf("<TD VALIGN=TOP>\n");
      printf("%s\n",pcddThis->cSourc);
      printf("</TD>\n");
      printf("</TR>\n");
      pcddThis = pcddThis->next;
    }
  } else { 
    pcddThis = pcdd;
    while (pcddThis) {
      strcpy(cCddId, pcddThis->cCddId);  CddUpCase(cCddId); 
      strcpy(cDescr, pcddThis->cDescr);  CddUpCase(cDescr); 
      strcpy(cSourc, pcddThis->cDescr);  CddUpCase(cSourc); 
      if (strstr(cCddId,cSrch) ||
          strstr(cDescr,cSrch) ||
	  strstr(cSourc,cSrch)) { 
	nhits++;
        printf("<TR>\n");
        printf("<TD NOWRAP VALIGN=TOP>\n");
        printf("<A HREF=\"%scddsrv.cgi?uid=%s\">\n",URLBase,pcddThis->cCddId);
        printf("%s</A>\n",pcddThis->cCddId);
        printf("</TD>\n");
        printf("<TD VALIGN=TOP>\n");
        printf("%s\n",pcddThis->cDescr);
        printf("</TD>\n");
        printf("<TD VALIGN=TOP>\n");
        printf("%s\n",pcddThis->cSourc);
        printf("</TD>\n");
        printf("</TR>\n");
      }
      pcddThis = pcddThis->next;
    }
  }
  printf("</TABLE>\n");
  if (nhits) {
    printf("<p>Number of CDs with matching terms: %d</p>\n",nhits);
  } else {
    printf("<H2>No matching terms found!</H2>\n");
  
  }
  printf("</BODY></HTML>\n");
  CddDescFree(pcdd);
  exit(0);
}

