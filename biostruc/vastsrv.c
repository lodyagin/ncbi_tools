/* ===========================================================================
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
 * File Name: vastsrv.c
 *
 * Author: Christopher Hogue, Tom Madej
 *
 * Version Creation Date: 10 March 1998
 *
 * $Log: vastsrv.c,v $
 * Revision 6.17  1999/05/07 14:02:13  zimmerma
 *  Removed local copy of MMDBBiostrucGet
 *
 * Revision 6.16  1999/02/09 15:14:00  addess
 * Modified by DZ to extract html path dependencies: created two new static char vars: DATApath and VASTpath, and text file getCn3D.txt in data and modified .vastrc to include new vars.
 *
 * Revision 6.15  1998/12/22 18:01:51  addess
 * changes relevant to reading new type of annot-set data
 *
 * Revision 6.14  1998/12/01  15:12:56  addess
 * removed unused variable pbsaShort
 *
 * Revision 6.13  1998/11/20  20:03:08  addess
 * related to platform independence of VAST Search
 *
 * Revision 6.12  1998/10/23  19:48:50  addess
 * fix pagination bug
 *
 * Revision 6.11  1998/10/21  14:07:13  addess
 * add error statement related to NR subset filtering
 *
 * Revision 6.10  1998/10/14  17:12:54  addess
 * pagination and aligned chain
 *
 * Revision 6.9  1998/08/06  17:48:52  madej
 * Fix Display/Show hits button for Vast Search.
 *
 * Revision 6.8  1998/07/27  16:08:15  madej
 * Minor change having to do with gifs.
 *
 * Revision 6.7  1998/07/17  18:48:14  madej
 * Miscellaneous changes.
 *
 * Revision 6.6  1998/06/11  19:13:57  madej
 * Modifications to html.
 *
 * Revision 6.5  1998/05/19  20:26:13  madej
 * Comment out VastViewAlign feature.
 *
 * Revision 6.4  1998/05/19  20:17:12  madej
 * Add general WWW routines for running on Sun servers.
 *
 * Revision 6.3  1998/03/30  19:11:50  madej
 * Added subset filtering, changes to neighbor page layout.
 *
 * Revision 6.2  1998/03/16  17:59:59  lewisg
 * added temporary hooks to communicate to cn3d
 *
 * Revision 6.1  1998/03/10 16:29:12  madej
 * First official version of vastsrv.c
 *
 */
 


/* Main program for VAST structure neighbor server. */

/* define DEBUG_1 1 */

#include <stdio.h>
#include <ncbi.h>
#include <accentr.h>
#include <netentr.h>
#include <www.h>
#include <sys/resource.h>
#include <asn.h>
#include <accutils.h>
#include <mmdbapi.h>
#include "vastlocl.h"
#include "mmdblocl.h"
#include "mmdbdata.h"
#include "vastsrv.h"



#define CPUTIME_MAX		120
#define DOMID_SIZE		6
#define MAX_MMDBIDS		4096
#define DEFAULT_SUBSET_NUM    	2		/* default to the NR set BLAST 10e-7 */
#define ABRIDGED_DISPLAY	0		/* display only Rmsd, Nres, %Id */
#define FULL_DISPLAY		1		/* display all columns in table */
#define SORT_BY_SCORE		1		/* sort table by Score */
#define SORT_BY_PVAL		2		/* sort table by P-value */
#define SORT_BY_RMSD		3		/* sort table by Rmsd */
#define SORT_BY_NRES		4		/* sort table by Nres */
#define SORT_BY_ID		5		/* sort table by %Id */
#define VIEW_ALIGNMENT		4		/* selects action "View Alignment" */
#define LOG10_500		2.69897		/* -log10(500); database size correction */
#define LOG_10			2.302585	/* log(10.0) */
#define NUM_HITS_PER_PAGE       20
#define DEFAULT_PAGE            1

static Char URLBase[PATH_MAX];
static Char URLcgi[PATH_MAX];
static Char DATApath[PATH_MAX];
static Char VASTpath[PATH_MAX];
static Char ENTREZurl[PATH_MAX];
static Char DOCSUMurl[PATH_MAX];
static Char MAILto[PATH_MAX];
static FILE *OutputFile = NULL;
static char OutputName[200];
static Char gunzip[PATH_MAX];
static Char MMDBpath[PATH_MAX];
static Int2 SortOn = 0;
static Int4 cnt_MMDBid;
static long save_MMDBid[MAX_MMDBIDS];
Boolean Chain;
Char VSPATH[PATH_MAX];
Char SlaveChain[2];
Char MasterChain[2];
 
static void
VnpHeapSort (ValNodePtr PNTR vnp, int (LIBCALLBACK *compar )PROTO ((Nlm_VoidPtr, Nlm_VoidPtr )))	
{
	Int2 index, total;
	ValNodePtr vnp1;
	ValNodePtr PNTR temp;

	total=0;
	for (vnp1 = *vnp; vnp1; vnp1=vnp1->next)
		total++;

	temp = (ValNodePtr PNTR) MemNew(total*sizeof(ValNodePtr));

	index=0;
	for (vnp1 = *vnp; vnp1; vnp1=vnp1->next)
	{
		temp[index] = vnp1;
		index++;
	}

	HeapSort ((VoidPtr) temp, (size_t) index, sizeof(ValNodePtr), compar);

	*vnp = temp[0];
	for (vnp1 = *vnp, index=0; index<(total-1); vnp1=vnp1->next, index++)
	{
		vnp1->next = temp[index+1];
	}
	vnp1 = temp[total-1];
	vnp1->next = NULL;

	temp = MemFree(temp);
}

 

/* Common code for the vastsrv page header is collected here.  This generates the html for
 * the title bar, the domain name, and the defline.
 */

static void
VastPageHeader(FILE *table, CharPtr pcPDB, Char cChain, int iDomain, Int4 iMMDBid, CharPtr JobID)
{
	DocSumPtr dsp;

	fprintf(table, "Content-type: text/html\n\n");
	fprintf(table, "<html><title>Vast Results</title>\n");
	fprintf(table, "<body bgcolor = \"ffffff\">\n");
	fprintf(table, "<base href=\"%s\">\n", URLBase);
	fprintf(table, "<img src=\"%s/vast2.gif\" alt=\"VAST Structure Neighbors\" usemap=#map ISMAP>\n", VASTpath);
	fprintf(table, "<map name=map>\n");
	fprintf(table, "<area shape=rect coords=4,5,42,18 href=\"http://www.ncbi.nlm.nih.gov\">\n");
	fprintf(table, "<area shape=rect coords=43,3,438,18 href=\"%s/vast.html\">\n", VASTpath);
	fprintf(table, "<area shape=rect coords=439,3,485,18 href=\"%s\">\n", ENTREZurl);
	fprintf(table, "<area shape=rect coords=490,3,507,18 href=\"%s/vasthelp.html\">\n", VASTpath);
	fprintf(table, "</map><p>\n");
	fprintf(table, "\n\n<HR SIZE=5 NOSHADE>\n");
	if (JobID == NULL)
        {
          fprintf(table, "<h1>Structures similar to MMDB\n");
	  fprintf(table, "<a href=\"%smmdbsrv?uid=%s", URLcgi, pcPDB);
	  fprintf(table, "&form=6&db=t&Dopt=s\">\n%d</a>", (long)iMMDBid);
	  fprintf(table, ",&nbsp\n");
        }
        else
        {
          fprintf(table, "<h1>Structures similar to VAST Search\n");
          fprintf(table, " %s\n", JobID);
          fprintf(table, ",&nbsp\n");
        }
	fprintf(table, "%s", pcPDB); 
  
	if (cChain != ' ')
		fprintf(table, " chain %c", cChain);

	if (iDomain > 0)
		fprintf(table, " domain %d", (int) iDomain);

	fprintf(table, "</h1>\n");
	dsp = NetDocSum(TYP_ST, (DocUid) iMMDBid);

	if (dsp) {
		if (dsp->title)
			fprintf(table, "%s\n", dsp->title);
 
		DocSumFree(dsp);
	}

} /* end VastPageHeader */


/************************** DZ: extracted from mmdbsrv.c ************************
 * WWWPrintFileData looks in the current CGI-BIN directory 
 *  or the "data" subdirectory for the data file.
 *  and prints it out to pFile
 */
 
static void WWWPrintFileData(CharPtr FName,  FILE *pFile)
{
 
   FILE *f = NULL;
   Char fullpath [PATH_MAX];
   CharPtr ptr;  
   Char pcBuf[1024];
   
   fullpath[0] = '\0';
   StringCpy(fullpath,  DATApath); /* look in DATApath */
   StringCat(fullpath,  FName);
   f = FileOpen (fullpath, "r");
   if (f == NULL) {
       f = FileOpen (FName, "r");  /* look in curent */
       if (f == NULL)  {  /* look in ./data/ */
         ProgramPath (fullpath, sizeof (fullpath) - 1);
         ptr = StringRChr (fullpath, DIRDELIMCHR);
         if (ptr != NULL) {
	  *ptr = '\0';
         }
         FileBuildPath (fullpath, "data", FName);  
         f = FileOpen (fullpath, "r");
         if (f == NULL)  {
           return;  
         } 
       }
   }
      
   do {
     pcBuf[0] = '\0';
     ptr = fgets(pcBuf, (size_t)1024, f);
     if (ptr) fprintf(pFile, ptr);
   } while (ptr);
  
   FileClose(f);
   return;
}
 

/* Note: A lot of this common html text could be put into a header file, read in, and then
 * spewed out.  Cf. the way this is done in mmdbsrv with WWWPrintFile.  Sometime I'll take
 * the time to write a more general version of the latter function that can be used for this
 * purpose.  T.M.
 */

/* 
   DZ: most of the HTML file URLs are inside the new VAST subdirectory, and references to these
   are interwoven with many query-specific references. Becasue of this interleaving, using 
   WWWPrintFileData would require several invocations and corresponding text files. Instead.
   WWWPrintFileData is only used to define the new path to CN3D, and all vast html references 
   are now prefixed with the VASTpath string. This variable has been added to .vastrc and to 
   the list of parameters to be extracted using GetVastParams. 
*/

static void
VastTableBegin (FILE *table, CharPtr pcPDB, CharPtr JobID, CharPtr pcPass, 
 Char cChain, int iDomain, Int4 iMMDBid, Int4 iFSID, 
 Int2 iFull, Int4 numhits, Int4 upper, Int4 lower, Int4 numpages, Int4 HitsPerPage, Int4 pagenum)
{
    Int4 DisplayOpt;

    VastPageHeader(table, pcPDB, cChain, iDomain, iMMDBid, JobID);
    fprintf(table,"<FORM METHOD=\"POST\" ACTION=\"%svastsrv\">\n", URLcgi);
    fprintf(table,"<INPUT TYPE=\"hidden\" NAME=\"uid\" VALUE=\"%ld\">\n", (long)iMMDBid);
    fprintf(table,"<INPUT TYPE=\"hidden\" NAME=\"chaindom\" VALUE=\"%ld\">\n", (long) iFSID);
    if (JobID)
    {
      fprintf(table,"<INPUT TYPE=\"hidden\" NAME=\"vsid\" VALUE=\"%s\">\n", JobID);
      fprintf(table,"<INPUT TYPE=\"hidden\" NAME=\"pass\" VALUE=\"%s\">\n", pcPass);
    } 

    fprintf(table, "<TABLE CELLPADDING=0 CELLSPACING=4>\n");
    /*****
    fprintf(table, "<TR><TH COLSPAN=2 ALIGN=LEFT><strong><INPUT TYPE=SUBMIT VALUE=\"View / Save Alignments\"></strong></TH></TR>\n");
    *****/
    fprintf(table, "<TR><TH COLSPAN=2 ALIGN=LEFT>\n");
    fprintf(table, "<strong><INPUT TYPE=SUBMIT VALUE=\"View / Save Alignments\"></strong>\n");
    fprintf(table, "&nbsp<img src=\"%s/new.gif\" alt=\"New\">\n", VASTpath);
    WWWPrintFileData("getCn3D.txt", table);
    fprintf(table, "</TH></TR>");
    fprintf(table, "<tr><td colspan = 1><strong>Options:</strong></td>\n");
    fprintf(table, "<td colspan = 1><strong>Viewer:</strong></td>\n"); 
    fprintf(table, "<td colspan = 2><strong>Complexity:</strong></td></tr>\n");
    /*fprintf(table, "<td colspan = 1><strong>Cn3D Atom Complexity:</strong></td></tr>\n");*/
    fprintf(table,"<TD VALIGN=TOP NOWRAP>\n");
    fprintf(table,"<INPUT TYPE=\"radio\" NAME=\"action\" value=\"0\"CHECKED> Launch Viewer<BR>\n");
    fprintf(table,"<INPUT TYPE=\"radio\" NAME=\"action\" value=\"1\"> See File<BR>\n");
    fprintf(table,"<INPUT TYPE=\"radio\" NAME=\"action\" value=\"2\"> Save File<BR></TD>\n");
    fprintf(table,"<TD VALIGN=TOP NOWRAP>\n");
    fprintf(table,"<INPUT TYPE=\"radio\" NAME=\"calltype\" value=\"a\"CHECKED> Cn3D v2.0 (asn.1)<BR>\n");
    fprintf(table,"<INPUT TYPE=\"radio\" NAME=\"calltype\" value=\"m\"> Mage (Kinemage)<BR>\n");
    fprintf(table,"<INPUT TYPE=\"radio\" NAME=\"calltype\" value=\"p\"> (PDB)<BR></TD>\n");
    fprintf(table, "<TD VALIGN=TOP NOWRAP>\n");
    fprintf(table, "<INPUT TYPE=\"radio\" NAME=\"chn_complexity\" value=\"1\"CHECKED> Aligned Chains only<BR>\n");
    fprintf(table, "<INPUT TYPE=\"radio\" NAME=\"chn_complexity\" value=\"0\"> All Chains<BR></TD>\n");
    fprintf(table, "<TD VALIGN=TOP NOWRAP>\n");
    fprintf(table, "<INPUT TYPE=\"radio\" NAME=\"atm_complexity\" value=\"1\"CHECKED> Alpha Carbons only<BR>\n");
    fprintf(table, "<INPUT TYPE=\"radio\" NAME=\"atm_complexity\" value=\"0\"> All Atoms<BR></TD>\n");
    fprintf(table, "</TR>\n</TABLE>\n");
    
    DisplayOpt = upper - lower;
    if (!DisplayOpt)
      fprintf(table, "<H4>Structure neighbor %d out of %d displayed. Page %d of %d.</H4>\n", lower, numhits, pagenum, numpages);
    else
      fprintf(table, "<H4>Structure neighbors %d-%d out of %d displayed. Page %d of %d.</H4>\n", lower, upper, numhits, pagenum, numpages);
    /* VAST data table begins here */
    fprintf(table,"<table cellspacing=3 cellpadding=2 width=100%% border=1>\n");
    fprintf(table,"<tr valign=middle>\n");
    fprintf(table,"<th>&nbsp</th>\n");
    fprintf(table,"<th align=left><pre> <a href=\"%s/vasthelp.html#Structure\">PDB</a>", VASTpath);
    fprintf(table," <a href=\"%s/vasthelp.html#C\">C</a>", VASTpath);
    fprintf(table," <a href=\"%s/vasthelp.html#D\">D</a></pre></th>\n", VASTpath);

    if (iFull) {
      fprintf(table,"<th align=right><pre><a href=\"%s/vasthelp.html#SCORE\">SCO</a></pre></th>\n", VASTpath); 
      fprintf(table,"<th align=right><pre><a href=\"%s/vasthelp.html#P-VAL\">P-VAL</a></pre></th>\n", VASTpath);
    }
    
    fprintf(table,"<th align=right><pre><a href=\"%s/vasthelp.html#RMSD\">RMSD</a></pre></th>\n", VASTpath);
    fprintf(table,"<th align=right><pre><a href=\"%s/vasthelp.html#NRES\">NRES</a></pre></th>\n",VASTpath);
    fprintf(table,"<th align=right><pre><a href=\"%s/vasthelp.html#Id\">%s</a></pre></th>\n", VASTpath, "%Id");
    fprintf(table,"<th align=left><pre><a href=\"%s/vasthelp.html#Contents\">Description</pre></th>\n",VASTpath);
    fprintf(table,"</tr><br>\n");
    fflush(table);

} /* end of VastTableBegin */



static void
VastTableRows(FILE *table, BiostrucFeatureSetPtr pbsfs, Int4 iMMDBid1, Int4 iFSID, Int2 iFull, ValNodePtr pvnBools)
{
   BiostrucFeaturePtr pbsf = NULL;
   ChemGraphAlignmentPtr pcga = NULL;
   BiostrucIdPtr pbsidThis = NULL;
   AlignStatsPtr pasp = NULL;
   ValNodePtr pvn =  NULL;
   DocSumPtr dsp = NULL;  
   CharPtr pcPDB;
   CharPtr pcSlaveName;
   Int4 iMMDBid;
   Char cChain;
   int iDomain;
   Int4 iAlign;
   float f;
   Boolean page;

   /* use cnt_MMDBid and save_MMDBid[] for docsum display of the correct subset hits */
   cnt_MMDBid = 0;

   for (pbsf = pbsfs->features; pbsf != NULL, pvnBools != NULL; pbsf = pbsf->next, pvnBools = pvnBools->next) {
       /* Filter Hits By Page*/
       page = pvnBools->data.boolvalue;
       if (page == FALSE) continue;
  
       /* get the embedded PDB code of the hit */
       pcPDB = StringSave(PDBNAME_DEFAULT);
       pcSlaveName = NULL;
       if (pbsf->name[14] != '\0') pcSlaveName = StringSave(&pbsf->name[14]);
       iDomain = 0;
       cChain = '-';
       if (StringLen(pbsf->name) >= 13)
        {
           pcPDB[0] = pbsf->name[7];
	   pcPDB[1] = pbsf->name[8];
	   pcPDB[2] = pbsf->name[9];
	   pcPDB[3] = pbsf->name[10];
	   cChain = pbsf->name[11];
	   iDomain = atoi((char *) &pbsf->name[12]);  
        }

       fprintf(table, "<tr>\n");
       fprintf(table, "<td VALIGN=TOP><INPUT TYPE=\"checkbox\" NAME=\"hit\"");
       fprintf(table, "VALUE=\"%ld\"></td>\n", (long) pbsf->id);
       pvn = ValNodeFindNext(pbsf->Location_location,NULL,Location_location_alignment);
       if (pvn) pcga = (ChemGraphAlignmentPtr) pvn->data.ptrvalue;
       iMMDBid = 0;
       if (pcga)
         {
	   
 	   pbsidThis = ValNodeFindNext(pcga->biostruc_ids,NULL,BiostrucId_mmdb_id);
	   if (pbsidThis)
	     {
		if (pbsidThis->next)  /* want only the second one - the slave */
		 iMMDBid = (long) pbsidThis->next->data.intvalue;
 	     }
         }

      /* save the MMDB id for later output for the docsum format */
      if (cnt_MMDBid < MAX_MMDBIDS)
         save_MMDBid[cnt_MMDBid++] = (long) iMMDBid;

      /* PDB and Kinemage file generators */
      /*****
      fprintf(table, "(<a href=\"%svastsrv?calltype=p&uid=%ld&fid=%ld&fsid=%ld&pdb=1\">P</a>) ",
		                       URLcgi, (long) iMMDBid1, pbsf->id, iFSID);
      fprintf(table, "(<a href=\"%svastsrv?calltype=m&uid=%ld&fid=%ld&fsid=%ld&pdb=1\">K</a>) ",
		                       URLcgi, (long) iMMDBid1, pbsf->id, iFSID);
      *****/
      /* Please do not delete.  for testing.  lyg */
      /*****
      fprintf(table, "(<a href=\"%svastsrv?calltype=c&uid=%ld&fid=%ld&fsid=%ld&pdb=1\">C</a>) ",
		                       URLcgi, (long) iMMDBid1, pbsf->id, iFSID);
      *****/

      /*****
      fprintf(table,"<a href=\"%smmdbsrv?uid=%ld&form=6&db=t&Dopt=s\">%ld</a></pre></td>\n",
          URLcgi, (long) iMMDBid, iMMDBid);
      *****/

      fprintf(table,"<td VALIGN=TOP><pre>");
      fprintf(table, "<a href=\"%smmdbsrv?uid=%ld&form=6&db=t&Dopt=s\">%s</a>", URLcgi,
         (long) iMMDBid, pcPDB);

      if (cChain == ' ')
         fprintf(table, "&nbsp ");
      else
         fprintf(table, " <a href=\"%smmdbsrv?uid=%ld&form=6&db=t&Dopt=s\">%c</a>",
            URLcgi, (long) iMMDBid, cChain);

      /* get the alignment number from the id */
      iDomain = (int) (pbsf->id/10) % 100;

      if (iDomain > 0)
         fprintf(table, " <a href=\"%smmdbsrv?uid=%ld&form=6&db=t&Dopt=s\">%d</a></pre></td>\n",
            URLcgi, (long) iMMDBid, iDomain);
      else
         fprintf(table, " </pre></td>\n");

       /* dig into aligndata */
       if (pcga)
         {
	     pasp = pcga->aligndata;
             if (iFull)
             {
               if (pasp->vast_score)
	         {
		   f = (FloatLo) pasp->vast_score;
		   f = f/(FloatLo) pasp->scale_factor;
		   fprintf(table,"<td VALIGN=top ALIGN=right>%.1f</td>\n",(float)f);
	         }
	       else
	         fprintf(table,"<td> </td>\n");  
	       if (pasp->vast_mlogp)
	         {
		   f = (float) pasp->vast_mlogp;
		   f = f/(float) pasp->scale_factor;

		   /* adjust for database size */
		   f -= LOG10_500;

		   if (f <= 4.0) {
		     f = (float) exp(-LOG_10*f);
		     fprintf(table, "<td VALIGN=top ALIGN=right>%.4f</td>\n", f);
		   }
		   else
		     fprintf(table,"<td VALIGN=top ALIGN=right>10e-%.1f</td>\n", f);
	         }
	       else
	          fprintf(table,"<td> </td>\n");
             }  
	     if (pasp->rmsd)
	       {
		   f = (FloatLo) pasp->rmsd;
		   f = f/(FloatLo) pasp->scale_factor;
		   fprintf(table,"<td VALIGN=top ALIGN=right>%.1f</td>\n",(float)f);
	       }
	     else
	        fprintf(table,"<td> </td>\n");  
	     if (pasp->align_res)
	      {
		  fprintf(table,"<td VALIGN=top ALIGN=center>%d</td>\n",(int) pasp->align_res);
	      }
             else
	       fprintf(table,"<td> </td>\n");  
	     if (pasp->other_score)
	      {
		   f = (FloatLo) pasp->other_score;
		   f = f/(FloatLo) pasp->scale_factor;
		   fprintf(table,"<td VALIGN=top ALIGN=right>%.1f</td>\n",(float)f * 100.0); 
	      }
             else
	       fprintf(table,"<td VALIGN=top ALIGN=right>0.0</td>\n");  
	 }
       else
        {
	    fprintf(table, "<td> </td><td> </td><td> </td><td> </td><td> </td>\n");
	}

      /* get the Entrez docsum */
       fprintf(table,"<td VALIGN=top>\n");
      /* Names are read directly from Annot-set. DocSums are used no longer */
       if(StringLen(pbsf->name) <= 13)
       {
         dsp = NULL;
         dsp = NetDocSum(TYP_ST, (DocUid) iMMDBid);
         if (dsp)
         {
           if (dsp->title)
 	   {
	    fprintf(table,"%s" , dsp->title);
	   }
	   else
	   {
	     fprintf(table,"  ");
	   }
	   DocSumFree(dsp);
         } 
         else
          fprintf(table, "   ");
       }  
       else
       {
         if (pcSlaveName) 
          fprintf(table,"%s" , pcSlaveName);
         else
          fprintf(table, "   ");
       }
       fprintf(table,"<BR>\n");
       fprintf(table,"</td>\n</tr>\n");
       fflush(table);
       MemFree(pcPDB);
       MemFree(pcSlaveName);

    } /* end of for */
 
} /* end of VastTableRows */



static void
VastTableEnd(FILE *table, Int4 iMMDBid, Int4 FSID, BiostrucAnnotSetPtr pbsas, Int4 subsetnum,
	Int4 iSort, Int2 iFull, CharPtr JobID, CharPtr pcPass, Int4 numpages, Int4 pagenum, Int4 HitsPerPage)
{
    BiostrucFeatureSetPtr pbsfs = NULL;
    BiostrucFeaturePtr pbsf = NULL;
    ChemGraphAlignmentPtr pcga = NULL;
    BiostrucIdPtr pbsidThis = NULL;
    ValNodePtr pvn = NULL;
    Int4 i, n;

    fprintf(table,"</table><hr>\n");
    fprintf(table,"</form>\n\n");
    fprintf(table,"<FORM METHOD=\"POST\" ACTION=\"%svastsrv\">\n", URLcgi);
    fprintf(table,"<INPUT TYPE=\"hidden\" NAME=\"uid\" VALUE=\"%ld\">\n", (long) iMMDBid);
    fprintf(table,"<INPUT TYPE=\"hidden\" NAME=\"chaindom\" VALUE=\"%ld\">\n",(long) FSID);
    if (JobID)
    {
      fprintf(table,"<INPUT TYPE=\"hidden\" NAME=\"vsid\" VALUE=\"%s\">\n", JobID);
      fprintf(table,"<INPUT TYPE=\"hidden\" NAME=\"pass\" VALUE=\"%s\">\n", pcPass);
    } 
    fprintf(table, "<TABLE CELLPADDING=0 CELLSPACING=4>\n");
    fprintf(table, "<TR><TD COLSPAN=3 ALIGN=LEFT> <strong><INPUT TYPE=SUBMIT VALUE=\"Display / Sort Hits\"></strong>\n");
    fprintf(table, "    &nbsp &nbsp page number:\n");
    fprintf(table, "<SELECT name=\"doclistpage\">\n");
    pagenum++;
    for (i = 1; i<= numpages; i++) {
      if (pagenum < numpages) {
        if (i == pagenum)
          fprintf(table,"<OPTION VALUE = \"%d\" SELECTED>%d\n", i, i);
        else 
          fprintf(table,"<OPTION VALUE = \"%d\">%d\n", i, i);
      }
      else {
        if (i == numpages)  
          fprintf(table,"<OPTION VALUE = \"%d\" SELECTED>%d\n", i, i);
        else 
          fprintf(table,"<OPTION VALUE = \"%d\">%d\n", i, i);
      }
    }

    fprintf(table, "</SELECT>\n");
    fprintf(table, "&nbsp <small>Hits to display per page:\n");
    fprintf(table, "<INPUT name = \"dispmax\" size=6 Value=%d> choose between 20-100 neighbors per page.</small></TD></TR>\n", HitsPerPage);
    fprintf(table, "<TR>\n");
    fprintf(table, "<TD><a href=\"%s/vasthelp.html#NRSet\"><strong>Display Subset:</strong></a></TD>\n", VASTpath);
    fprintf(table, "<TD COLSPAN=1><strong>Sorted by:</strong></TD>\n");
    fprintf(table, "<TD COLSPAN=1><strong>Column Format:</strong></TD>\n");
    fprintf(table, "</TR>\n");
    fprintf(table, "<TR><TD VALIGN=TOP NOWRAP>\n");
    n = GetNumberOfSubsets();

    /* subset 1 should be "All of PDB"; put it last */
    for (i = 2; i <= n; i++) {
      fprintf(table, "<INPUT TYPE=radio NAME=subset VALUE=\"%s\"", GetSubsetName(i));

      if (i == subsetnum)
	fprintf(table, " CHECKED");

      fprintf(table, "> %s<BR>\n", GetSubsetName(i));
    }

    fprintf(table, "<INPUT TYPE=radio NAME=subset VALUE=\"%s\"", GetSubsetName(1));

    if (subsetnum == 1)
      fprintf(table, " CHECKED");

    fprintf(table, "> %s<BR>\n", GetSubsetName(1));
    fprintf(table, "<TD VALIGN=TOP NOWRAP>\n");
    fprintf(table, "<INPUT TYPE=radio NAME=sort VALUE=\"1\"");
    if ((iSort == SORT_BY_SCORE) || (iSort == 0)) fprintf(table, " CHECKED");
    fprintf(table, "> VAST Score<BR>\n");
    fprintf(table, "<INPUT TYPE=radio NAME=sort VALUE=\"2\"");
    if (iSort == SORT_BY_PVAL) fprintf(table, " CHECKED");
    fprintf(table, "> VAST P-value<BR>\n");
    fprintf(table, "<INPUT TYPE=radio NAME=sort VALUE=\"3\"");
    if (iSort == SORT_BY_RMSD) fprintf(table, " CHECKED");
    fprintf(table, "> Rmsd<BR>\n");
    /*****
    fprintf(table, "<TD VALIGN=TOP NOWRAP>\n");
    *****/
    fprintf(table, "<INPUT TYPE=radio NAME=sort VALUE=\"4\"");
    if (iSort == SORT_BY_NRES) fprintf(table, " CHECKED");
    fprintf(table, "> Aligned residues<BR>\n");
    fprintf(table, "<INPUT TYPE=radio NAME=sort VALUE=\"5\"");
    if (iSort == SORT_BY_ID) fprintf(table, " CHECKED");
    fprintf(table, "> Identities<BR>\n");
    fprintf(table, "<TD VALIGN=TOP NOWRAP>\n");
    fprintf(table, "<INPUT TYPE=radio NAME=version VALUE=\"0\"");
    if (iFull == ABRIDGED_DISPLAY) fprintf(table, " CHECKED");
    fprintf(table, "> RMSD, NRES, %s<BR>\n", "%Id");
    fprintf(table, "<INPUT TYPE=radio NAME=version VALUE=\"1\"");
    if (iFull == FULL_DISPLAY) fprintf(table, " CHECKED");
    fprintf(table, "> All values<BR>\n");
    fprintf(table, "</TR>\n");
    fprintf(table, "</TABLE>\n");
    fprintf(table, "</FORM>\n");

    /***** Removing docsum button:
    fprintf(table,"<FORM METHOD=\"POST\" ACTION=\"%s/query\">\n", DOCSUMurl);
    fprintf(table,"<INPUT TYPE=\"hidden\" NAME=\"form\" VALUE=\"6\">\n");
    fprintf(table,"<INPUT TYPE=\"hidden\" NAME=\"db\" VALUE=\"t\">\n");
    fprintf(table,"<INPUT TYPE=\"hidden\" NAME=\"uid\" VALUE=\"%ld,", (long) iMMDBid);
    *****/
   
    /* output the MMDB id's for the slaves */
    /***** Removing docsum button:
    if (cnt_MMDBid > 0) {
      for (i = 0; i < cnt_MMDBid - 1; i++)
	fprintf(table, "%ld,", save_MMDBid[i]);

      fprintf(table, "%ld ", save_MMDBid[cnt_MMDBid - 1]);
    }

    fprintf(table,"\">\n");
    fprintf(table,"<INPUT TYPE=\"submit\" VALUE=\"Display\"> table as an Entrez Document Summary form.</form>\n");
    *****/
    fprintf(table, "\n<HR SIZE=5 NOSHADE>\n");
    fprintf(table,"</body></html>\n");

} /* end of VastTableEnd */

/* I decided to filter by domain subset first because it make writing the pagination code much easier Ken */

static BiostrucFeaturePtr FilterHitsByDomainSubset(BiostrucFeaturePtr pbsf, Int4 subsetnum)
{
  BiostrucFeaturePtr current, pbsfHead = NULL, pbsfTail;
  Int4 h, i, n;
  Int4 gn, gr, hcnt, *min_ranks, *group_num, *group_rank;
  Char domid[DOMID_SIZE + 1];
 
/* The next bit of code is used for filtering the hit lists.  When we go through a hit
 * list we skip domains that do not belong to the subset of interest, or which belong to
 * the subset but for which a group representative has already been encountered.
*/
  for (i = 0; i <= DOMID_SIZE; i++)
    domid[i] = '\0';

  n = GetNumberOfDomains();
  min_ranks = (Int4 *) MemNew(n*sizeof(Int4));
  group_num = (Int4 *) MemNew(n*sizeof(Int4));
  group_rank = (Int4 *) MemNew(n*sizeof(Int4));

/* use a first pass to select those hits belonging to the specified subset */
  if (group_rank != NULL) {
      /* initializations */
    for (i = 0; i < n; i++) {
         min_ranks[i] = n + 1;
         group_num[i] = 0;
         group_rank[i] = n + 1;
    }

    for (current = pbsf, hcnt = 0; current != NULL; current = current->next, hcnt++) {
       /* copy domain identifier into domid[] */
       domid[0] = current->name[7];
       domid[1] = current->name[8];
       domid[2] = current->name[9];
       domid[3] = current->name[10];
       domid[4] = current->name[11];
       domid[5] = current->name[12];

       if (domid[5] == '0') domid[5] = ' ';

       /* skip over domains that do not belong to the subset */
       if (BelongsToSubset(domid, subsetnum, &gn, &gr) <= 0)
          continue;

       /* record group data for this hit */
       group_num[hcnt] = gn;
       group_rank[hcnt] = gr;

       /* is this the minimum group rank for this group? */
       if (gr < min_ranks[gn - 1])
	  min_ranks[gn - 1] = gr;
    }
  }
  current = pbsf;
  hcnt = 0;
  while (current) {
    if (group_rank != NULL) 
    {
      /* subset filtering done here */
      if (group_num[hcnt] == 0)
      {
        hcnt++;
        current = current->next;
        continue;
      }
   
      gn = group_num[hcnt];

      if (group_rank[hcnt] != min_ranks[gn - 1])
      {
        hcnt++;
        current = current->next;
        continue;
      }
    }
    
    if (pbsfHead == NULL)
    {
      pbsfHead = BiostrucFeatureNew();
      pbsfHead = current;
      pbsfTail = pbsfHead;
      hcnt++;
      current = current->next;
      pbsfTail->next = NULL;
    }
    else
    {
      pbsfTail->next = current;
      pbsfTail = current;
      hcnt++;
      current = current->next;
      pbsfTail->next = NULL;
    }
  }
  
  MemFree(min_ranks);
  MemFree(group_num);
  MemFree(group_rank);
  return pbsfHead;
}
 
static ValNodePtr FilterHitsByPage(BiostrucFeatureSetPtr pbsfs, Int4 PageNum, Int4 HitsPerPage, Int4 *numhits, Int4 *numpages, Int4 *upper, Int4 *lower)
{
  
  BiostrucFeaturePtr pbsf;
  Int4 index, FidCount, RemainFids, CompleteFidSet; 
  Int4 UpperLimit, LowerLimit;
  ValNodePtr pvnBools = NULL;
  float n;
  
  pbsf = pbsfs->features;
  FidCount = 0;
  
  while (pbsf)
  {
    FidCount++;
    pbsf = pbsf->next;
  }
  *numhits = FidCount;

  if (FidCount <= HitsPerPage) 
    *numpages = 1;
  else
  { 
    RemainFids = FidCount % HitsPerPage;
  
    if (RemainFids)
    {
     CompleteFidSet = FidCount - RemainFids;
     *numpages = (CompleteFidSet/HitsPerPage) + 1;
    }
    else *numpages = FidCount/HitsPerPage;
  }
  
  UpperLimit = HitsPerPage * PageNum;
  LowerLimit = UpperLimit - HitsPerPage + 1;
    
  if ((FidCount < HitsPerPage ) || (LowerLimit > FidCount)) {
    UpperLimit = FidCount;
    LowerLimit = 1;
  }
  
  if (UpperLimit > FidCount) UpperLimit = FidCount;
  
  *lower = LowerLimit;
  *upper = UpperLimit;
  
  for (index = 1; index <= FidCount; index++)
  {
    if ((index >= LowerLimit) && (index <= UpperLimit)) ValNodeAddBoolean(&pvnBools, 0, TRUE);
    else ValNodeAddBoolean(&pvnBools, 0, FALSE);
  }
  
  return pvnBools;
}

static void
MakeVastTable(Int4 FSID, BiostrucAnnotSetPtr pbsas, Int2 iSort, Int4 subsetnum, Int4 pagenum, Int4 HitsPerPage, Int2 iFull, CharPtr JobID, CharPtr pcPass)
{
  BiostrucFeatureSetPtr pbsfs = NULL, pbsfs2 = NULL;
  FILE *table = NULL;
  char tableName[200];
  CharPtr pcPDB = NULL;
  char cChain = '-';
  int iDomain = 0;
  Int4 iMMDBid = 0;
  Int4 numhits, numpages, upper, lower;
  BiostrucIdPtr pbsidThis = NULL;
  BiostrucDescrPtr pbsdrThis = NULL;
  CharPtr pcVast = NULL;
  ValNodePtr pvnBools;

  if ((!pbsas) || (!FSID) ) return;
	
  /* pull values out of BiostrucAnnotSet */
  pbsidThis = ValNodeFindNext(pbsas->id,NULL,BiostrucId_mmdb_id);
  if (pbsidThis)
     {
 	 iMMDBid = (Int4) pbsidThis->data.intvalue;  /* Get MMDB id no (only the first one) */ 
     }
  else 
     {
 	printf("Content-type: text/html\n\n");
        printf("<h2>Error</h2>\n");
        printf("Internal VASTSERV Error.  No MMDB-ID in Data on Server.<p>\n");
        return;
     }

  pbsfs = pbsas->features;
  while (pbsfs)
   {
     if (pbsfs->id == FSID)
      { 
        /* got the right one - make the table */
        
	/* pull out the PDB chain code and domain number from the name string */
        pbsidThis = ValNodeFindNext(pbsfs->descr,NULL,BiostrucFeatureSetDescr_name);
        if (pbsidThis)
         {
 	   pcVast = (CharPtr) pbsidThis->data.ptrvalue;  
	   pcPDB = StringSave(PDBNAME_DEFAULT);
	   iDomain = 0;
	   cChain = '-';
           if (StringLen(pcVast) >= 6)
	     { 
	       pcPDB[0] = pcVast[0];
	       pcPDB[1] = pcVast[1];
	       pcPDB[2] = pcVast[2];
	       pcPDB[3] = pcVast[3];
	       cChain = pcVast[4];
	       iDomain = (Int4) FSID % 100;
             }
         }

	if (EntrezInit("vastsrv", NULL, FALSE) == FALSE)
      	  {
	    printf("Content-type: text/html\n\n");
	    printf("<h2>Error</h2>\n");
	    printf("VASTSERV: EntrezInit Failed.<p>\n");
	    if (pcPDB) MemFree(pcPDB);
	    return;      
	  }
 
 	strcpy(tableName,GetTempName("vastsrv")); 
	if  (!(table = FileOpen(tableName,WRITE)))
	   {
	     printf("Content-type: text/html\n\n");
	     printf("<h2>Error</h2>\n");
	     printf("VASTSERV: Temp File Open Failed on Server<p>\n");
             goto oot;
	   }  

        if (pbsfs->features != NULL) {
            pbsfs2 = BiostrucFeatureSetNew();
            pbsfs2->id = pbsfs->id;
            pbsfs->id = NULL;
            pbsfs2->descr = pbsfs->descr;
            pbsfs->descr = NULL;
            pbsfs2->features = FilterHitsByDomainSubset(pbsfs->features, subsetnum);
            if (pbsfs2->features == NULL)
            {
              printf("Content-type: text/html\n\n");
	      printf("<h2>Error</h2>\n");
	      printf("VASTSERV: Hits are not present in non-redundant subset<p>\n");
              goto oot;
            }  
            pbsfs->features = NULL;
            pvnBools = ValNodeNew(NULL);
            pvnBools = FilterHitsByPage(pbsfs2, pagenum, HitsPerPage, &numhits, &numpages, &upper, &lower);
            if ((numhits < HitsPerPage) || (lower > numhits)) pagenum = DEFAULT_PAGE;
            VastTableSort(pbsfs2, iSort);
 	    VastTableBegin(table, pcPDB, JobID, pcPass, cChain, iDomain, iMMDBid, FSID, iFull, numhits, upper, lower, numpages, HitsPerPage, pagenum);
 	    VastTableRows(table, pbsfs2, iMMDBid, FSID, iFull, pvnBools);
            VastTableEnd(table, iMMDBid, FSID, pbsas, subsetnum, iSort, iFull, JobID, pcPass, numpages, pagenum, HitsPerPage);
	}
	else {
	    VastPageHeader(table, pcPDB, cChain, iDomain, iMMDBid, JobID);
	    fprintf(table, "<h1><a href=\"%s/vasthelp.html#NoNeighbor\">%s</a></h1>\n", VASTpath,
		"VAST did not find any structure neighbors.");
	    fprintf(table, "<HR SIZE=5 NOSHADE>\n");
	    fprintf(table, "</body></html>\n");
	}

	fflush(table);

	if (table != stdout)
	  {
	    fclose(table);
	    PrintFile(tableName);
          }
       oot:
        EntrezFini();
	RemoveTempFiles();  
	if (pcPDB) MemFree(pcPDB);
        return;
      }     
     pbsfs = pbsfs->next; 
   }
  printf("Content-type: text/html\n\n");
  printf("<h2>Error</h2>\n");
  printf("VASTSERV: Could Not Create VAST Table because Data Not Found.<p>\n");
  if (pcPDB) MemFree(pcPDB);
  return;       
}



BiostrucAnnotSetPtr LIBCALL LocalGetBiostrucAnnotSet(Int4 mmdbid, CharPtr JobID)
{
    FILE *pFile;
    AsnIoPtr aip = NULL;
    AsnTypePtr atp = NULL;
    BiostrucAnnotSetPtr pbsa = NULL;
    Char path[PATH_MAX];
    Char pcId[20];    
    Int2 iFileExists = 0;

   sprintf(pcId, "/%ld", (long) mmdbid);
   pFile = NULL;
   StringCpy(path, VSPATH);
   StringCat(path, JobID);
   StringCat(path, pcId);
   StringCat(path, ".bas");
   iFileExists = FileLength(path);
    if (iFileExists == 0)
      {
        return NULL;
      }
      
    aip = AsnIoOpen(path, "r");
    pbsa = BiostrucAnnotSetAsnRead(aip, NULL);  
    AsnIoClose (aip);
    if (!pbsa) return NULL;  
    return pbsa;
} 



BiostrucAnnotSetPtr LIBCALL
LocalGetFeatureSet(Int4 mmdbid, Int4 feature_set_id, CharPtr JobID)
{
    BiostrucAnnotSetPtr basp2 = NULL;
    BiostrucFeatureSetPtr pbsfs = NULL;
    BiostrucAnnotSetPtr basp = NULL;
    BiostrucFeatureSetPtr pbsfsLast = NULL;
    
    if (IsVASTData(mmdbid))
       basp = VASTBsAnnotSetGet(mmdbid);
    else if (IsVASTData(feature_set_id)) {
       basp = VASTBsAnnotSetGet(feature_set_id);
       if (basp != NULL) return basp;
    } 
    else if (JobID)
       basp = LocalGetBiostrucAnnotSet(mmdbid, JobID);

    if (basp == NULL)
	return NULL;
 
    pbsfs = basp->features;
    pbsfsLast  = NULL;
    basp2 = NULL;
    while (pbsfs)
     {
       if (pbsfs->id == feature_set_id)
        {
          basp2 = BiostrucAnnotSetNew();
          basp2->id = basp->id;
          basp->id = NULL; /* unlink the id valnode from basp object */
          basp2->descr = basp->descr; 
          basp->descr = NULL;  /* unlink the descr from basp object */
          basp2->features = pbsfs;
	  if (pbsfsLast) /* relink next to prev */
	    pbsfsLast->next = pbsfs->next;
          else 
	    basp->features = pbsfs->next;	  
	  basp2->features->next = NULL;
	  BiostrucAnnotSetFree(basp);
          return basp2;
        }
       pbsfsLast = pbsfs;
       pbsfs = pbsfs->next;
     }   
    BiostrucAnnotSetFree(basp);
    return basp2;
}



BiostrucAnnotSetPtr  PruneBiostrucAnnotHits( 
    BiostrucAnnotSetPtr basp, Int4 FSID, ValNodePtr pvnFids)
{
    BiostrucAnnotSetPtr basp2 = NULL;
    BiostrucFeatureSetPtr pbsfs = NULL;
    BiostrucFeaturePtr pbsf = NULL;
    BiostrucFeaturePtr pbsfHold = NULL;
    BiostrucFeaturePtr pbsf2 = NULL;
    BiostrucFeaturePtr pbsfDie = NULL;
 
    ValNodePtr pvnFid = NULL;
    Boolean found = FALSE;

    if ((basp == NULL) || (pvnFids == NULL) || (FSID == 0))
	return NULL;
 
    pbsfs = basp->features;
    while (pbsfs)
     {
       if (pbsfs->id == FSID)
        {
	  basp2 = BiostrucAnnotSetNew();
     	  basp2->id = basp->id;
	  basp->id = NULL; /* unlink the id valnode from basp object */
    	  basp2->descr = basp->descr; 
    	  basp->descr = NULL;  /* unlink the descr from basp object */
	  basp2->features = BiostrucFeatureSetNew();
          basp2->features->id = pbsfs->id;
	  pbsfs->id = NULL;
          basp2->features->descr = pbsfs->descr;
          pbsfs->descr = NULL; /* unlink the feature-set descr from basp  object */
	  pbsfHold = pbsfs->features;
	  pbsfs->features = NULL;
	  BiostrucAnnotSetFree(basp);
	  pbsfs = NULL;
          pbsf = pbsfHold;
	  while (pbsf)
	    {
	        found = FALSE;
		pvnFid = pvnFids;
	        while (pvnFid)
		 { 
		     if (pbsf->id == (Int4) pvnFid->data.intvalue)
		       {
			 found = TRUE;
			 break;
		       }
		     pvnFid = pvnFid->next;
		 }
	       if (!found) 
	         {
		   pbsfDie = pbsf;
	           pbsf = pbsf->next;
		   pbsfDie->next = NULL;
		 }
	       else
	         {
		   if (!basp2->features->features)
		     {
			 basp2->features->features = pbsf;
			 pbsf2 = basp2->features->features;
			 pbsf = pbsf->next; /* keep next for loop */
			 pbsf2->next = NULL; /* chop it */
		     }
		    else
		     {
		         pbsf2->next = pbsf;  
			 pbsf2 = pbsf; 
			 pbsf = pbsf->next; /* keep next for loop */
			 pbsf2->next = NULL; /* chop it */
		     } 
		 }	 
	    }  /* while pbsf */
	 }
	if(pbsfs) pbsfs = pbsfs->next;
      }
   return basp2;
}



BiostrucAnnotSetPtr LIBCALL  BiostrucAnnotSetGetByFid (
    BiostrucAnnotSetPtr basp, Int4 feature_id, Int4 feature_set_id)
{
     BiostrucAnnotSetPtr basp2 = NULL;
    BiostrucFeatureSetPtr pbsfs = NULL;
    BiostrucFeaturePtr pbsf = NULL;

    if (basp == NULL)
	return NULL;
 
    pbsfs = basp->features;
    while (pbsfs)
     {
       if (pbsfs->id == feature_set_id)
        {
          pbsf =  pbsfs->features;
          while(pbsf)
            {
              if (pbsf->id == feature_id)
                {  /* found it */
                     basp2 = BiostrucAnnotSetNew();
     		     basp2->id = basp->id;
		     basp->id = NULL; /* unlink the id valnode from basp object */
    		     basp2->descr = basp->descr; 
    		     basp->descr = NULL;  /* unlink the descr from basp object */
    		     basp2->features = BiostrucFeatureSetNew();
                     basp2->features->id = pbsfs->id;
                     basp2->features->descr = pbsfs->descr;
                     pbsfs->descr = NULL; /* unlink the feature-set descr from basp  object */
                     basp2->features->features = BiostrucFeatureNew();
                     basp2->features->features->id = pbsf->id;
                     basp2->features->features->name = StringSave(pbsf->name);
		     basp2->features->features->type = pbsf->type;
		     basp2->features->features->Property_property = pbsf->Property_property;
		     pbsf->Property_property = NULL; /* unlink the property from basp  object */
		     basp2->features->features->Location_location = pbsf->Location_location;
		     pbsf->Location_location = NULL; /* unlink the location from basp  object */ 
		     BiostrucAnnotSetFree(basp);
                     return basp2;
                }
               pbsf = pbsf->next;
            }
        }
       pbsfs = pbsfs->next;
     }
   
    BiostrucAnnotSetFree(basp);
    return basp2;
}



Int2 LIBCALL Check_VastSearch_Password(CharPtr pcPassNew, CharPtr JobID)
{
 
  Char pcPassFile[24];
  Char PassPath[PATH_MAX];
  FILE *passwdfile;
  CharPtr pcPassOld;
  Int4 iPassLen;

  iPassLen = StringLen(pcPassNew);
  pcPassOld = StringSave(pcPassNew);
   
  sprintf(pcPassFile, "/%s.passwd", JobID);
  PassPath[0]='\0';
  StringCpy(PassPath, VSPATH);
  StringCat(PassPath, JobID);
  StringCat(PassPath, pcPassFile); 
  
  if ((passwdfile = FileOpen(PassPath, "r")) == NULL)
  {
   
    printf("Content-type: text/html\n\n");
    printf("<h2>Error</h2>\n");
    printf("<body><h2>Cannot examine password</h2>\n");
    printf("Please alert info@ncbi.nlm.nih.gov\n");
    printf("of this problem\n");
    printf("</BODY>\n</HTML>\n");
    return 2;
  }
  
  fscanf(passwdfile, "%s", pcPassOld);
  FileClose(passwdfile);

  if (!StringNCmp(pcPassOld, pcPassNew, iPassLen))
    return 1;
  else 
    return 0;
}

/* Extract vastsrv parameters from the config file. */

static Boolean
GetVastParams()
{
	URLBase[0] = URLcgi[0] = ENTREZurl[0] = DOCSUMurl[0] = MAILto[0] = '\0';
	MMDBpath[0] = gunzip[0] = '\0';

	GetAppParam("vast", "VASTSRV", "URLBase", "", URLBase, PATH_MAX);

	if (URLBase[0] == '\0') {
		ErrPostEx(SEV_FATAL, 0, 0, "VAST config file\nVASTSRV section has no URLBase...\n");
		return FALSE;
	}
 
	GetAppParam("vast", "VASTSRV", "URLcgi", "", URLcgi, PATH_MAX);

	if (URLcgi[0] == '\0') {
		ErrPostEx(SEV_FATAL, 0, 0, "VAST config file\nVASTSRV section has no URLcgi...\n");
		return FALSE;
	}
 
	GetAppParam("vast", "VASTSRV", "ENTREZurl", "", ENTREZurl, PATH_MAX);

	if (ENTREZurl[0] == '\0') {
		ErrPostEx(SEV_FATAL, 0, 0, "VAST config file\nVASTSRV section has no ENTREZurl...\n");
		return FALSE;
	}
 
	GetAppParam("vast", "VASTSRV", "DOCSUMurl", "", DOCSUMurl, PATH_MAX);

	if (DOCSUMurl[0] == '\0') {
		ErrPostEx(SEV_FATAL, 0, 0, "VAST config file\nVASTSRV section has no DOCSUMurl...\n");
		return FALSE;
	}

	GetAppParam("vast", "VASTSRV", "Gunzip", "", gunzip, (size_t) 256*(sizeof(char)));

	if (gunzip[0] == '\0') {
		ErrPostEx(SEV_FATAL, 0, 0, "VAST config file\nVASTSRV section has no Gunzip...\n");
		return FALSE;
	}

	GetAppParam("mmdb", "MMDB", "Database", "", MMDBpath, PATH_MAX);

	if (MMDBpath[0] == '\0') {
		ErrPostEx(SEV_FATAL, 0, 0, "MMDB config file\nMMDBSRV section has no Database...\n");
		return FALSE;
	}

	GetAppParam("vast", "VASTSRV", "MAILto", "", MAILto, PATH_MAX);

	if (MAILto[0] == '\0') {
		ErrPostEx(SEV_FATAL, 0, 0, "VAST config file\nVASTSRV section has no MAILto...\n");
		return FALSE;
	}
 
	 GetAppParam("vast", "VASTSRV", "VSPATH", "", VSPATH, PATH_MAX);

        if (VSPATH[0] == '\0') {
		ErrPostEx(SEV_FATAL, 0, 0, "VAST config file\nVASTSRV section has no VAST Search path...\n");
		return FALSE;
	}
        
	GetAppParam("vast", "VASTSRV", "DATApath", "", DATApath, PATH_MAX);
        if (DATApath[0] == '\0') {
		ErrPostEx(SEV_FATAL, 0, 0, "VAST config file\nVASTSRV section has no VAST Data path...\n");
		return FALSE;
	}

	GetAppParam("vast", "VASTSRV", "VASTpath", "", VASTpath, PATH_MAX);
        if (DATApath[0] == '\0') {
		ErrPostEx(SEV_FATAL, 0, 0, "VAST config file\nVASTSRV section has no VAST html path...\n");
		return FALSE;
	}

        return TRUE;

} /* end GetVastParams */



Int2
Main()
{
	FILE *pFile = NULL, *pIn = NULL;
	Char pcBuf[100];
	CharPtr pcTest, pcL1 = NULL;
	Int4 GetGi, Fid, Fsid, iFileExists = 0, NumLabels = 0;
	BiostrucAnnotSetPtr pbsa = NULL;
	PDNMS pdnmsMaster = NULL, pdnmsSlave = NULL;
	AsnIoPtr aip = NULL;
	Int2 iTest = 0, iPDB = 0, iSort = 0, action = 0, viewer = 0, level = 0;
	CharPtr Name, Value, IPAddress = getenv("REMOTE_HOST");
	struct rlimit rl;
	ValNodePtr pvnFid = NULL, pvnFids = NULL;
	Int4 iFidCount = 0, count = 0, subsetnum, pagenum, HitsPerPage, indx;
	Char subsetname[256];
	CharPtr JobID = NULL, pcPass, www_arg;
	Int2 ret, iFull = 0;
	WWWInfoPtr www_info;

 

	/* this sets up the unix time limit */
	getrlimit(RLIMIT_CPU, &rl);
	rl.rlim_max = rl.rlim_cur = CPUTIME_MAX;
	setrlimit(RLIMIT_CPU, &rl);

	if (!GetVastParams()) {
		printf("Content-type: text/html\n\n");
		printf("<h2>VASTSERV Error</h2>\n");
		printf("<h3>Couldn't read from config file...</h3>\n");
		exit(1);
	}
    
	if (WWWGetArgs(&www_info) != WWWErrOk) {
		printf("Content-type: text/html\n\n");
		printf("<h2>VASTSERV</h2>\n");
		printf("<h3>Failed to process posting - check your get/post syntax.</h3>\n");
		exit(1);
	}

	if ((NumLabels = WWWGetNumEntries(www_info)) == 0) {
		printf("Content-type: text/html\n\n");
		printf("<h2>VASTSERV</h2>\n");
		printf("<h3>No input - nothing to report.</h3>\n");
		exit(1);
	}

	if ((indx = WWWFindName(www_info, "action")) >= 0) {
		www_arg = WWWGetValueByIndex(www_info, indx);

		if (isdigit(www_arg[0]))
			action = (Int2) atoi(www_arg);
		else
			/* default to asn.1 text */
			action = 3;
	}
	else
		action = 3;

	/***** We may not add this feature, comment it out for now.
	if (action == VIEW_ALIGNMENT) {
		(void) VastViewAlign(www_info);
		exit(0);
	}
	*****/

	/* check whether or not to launch a viewer */
	if ((indx = WWWFindName(www_info, "calltype")) >= 0) {
		www_arg = WWWGetValueByIndex(www_info, indx);
	
		switch (www_arg[0]) {
		case 'm':
			(void) VastToMage(www_info);
			exit(0);
		case 'p':
			(void) VastToPDB(www_info);
			exit(0);
		case 'a':
			(void) VastToCn3D(www_info);
			exit(0);
		default:
			printf("Content-type: text/html\n\n");
			printf("<h2>VASTSERV Error</h2>\n");
			printf("<h3>Internal failure (bad calltype).\nContact %s</h3>\n", MAILto);
			exit(1);
		}
	}
	
	if (!VASTInit()) {
		printf("Content-type: text/html\n\n");
		printf("<h2>VASTSERV Error</h2>\n");
		printf("<h3>Cannot find VAST data on server.\nContact %s</h3>\n", MAILto);
		exit(1);
	}
 
	if ((indx = WWWFindName(www_info, "chaindom")) < 0) {
		printf("Content-type: text/html\n\n");
		printf("<h2>VASTSERV Error</h2>\n");
		printf("<h3>Internal failure (no chaindom).\nContact %s</h3>\n", MAILto);
		exit(1);
	}

	www_arg = WWWGetValueByIndex(www_info, indx);

	if (isdigit(www_arg[0]))
		Fsid = (Int4) atol(www_arg);
	else {
		printf("Content-type: text/html\n\n");
		printf("<h2>VASTSERV Error</h2>\n");
		printf("<h3>Invalid feature set id input; no results.</h3>\n");
		exit(1);
	}

	GetGi = Fsid/10000;

	if ((indx = WWWFindName(www_info, "vsid")) >= 0) {
		/* we have a VAST Search job */
		www_arg = WWWGetValueByIndex(www_info, indx);
		JobID = StringSave(www_arg);

		if ((indx = WWWFindName(www_info, "pass")) < 0) {
			printf("Content-type: text/html\n\n");
			printf("<body bgcolor = \"#f0f0f0\"\n");
			printf("<h2>VAST SEARCH</h2>\n");
			printf("<h3>Password required.</h3>\n");
			exit(0);
		}

		www_arg = WWWGetValueByIndex(www_info, indx);
		pcPass = StringSave(www_arg);

		if ((ret = Check_VastSearch_Password(pcPass, JobID)) != 1) {
			if (ret == 2) exit(0);
			printf("Content-type: text/html\n\n");
			printf("<body bgcolor = \"#f0f0f0\"\n");
			printf("<h2>VAST SEARCH</h2>\n");
			printf("<h3>Incorrect password.</h3>\n");
			exit(0);
		}
	}

	/* load in the chaindom into memory */
	objmmdb1AsnLoad();
	objmmdb2AsnLoad();
	objmmdb3AsnLoad();
	pbsa = LocalGetFeatureSet(GetGi, Fsid, JobID);	 

	if (pbsa == NULL) {
		printf("Content-type: text/html\n\n");
		printf("<body bgcolor = \"#f0f0f0\">\n");
		printf("<br>\n<h2>VAST structure neighbor calculations for this entry are in progress.</h2>\n");
		exit(0);
	}

	/* at this point, there is a valid feature set id and pbsa in memory */
   
	/* subset filtering; identify which subset we're working with */
	if ((indx = WWWFindName(www_info, "subset")) < 0)
		subsetnum = DEFAULT_SUBSET_NUM;
	else {
		www_arg = WWWGetValueByIndex(www_info, indx);
		StringCpy(subsetname, www_arg);
		subsetnum = GetSubsetNum(subsetname);
	}

	if ((indx = WWWFindName(www_info, "doclistpage")) < 0)
		pagenum = DEFAULT_PAGE;
	else {
		www_arg = WWWGetValueByIndex(www_info, indx);
		if (isdigit(www_arg[0]))
                  pagenum = (Int4) atoi(www_arg);
                else
                  pagenum = DEFAULT_PAGE;
        }
           
        if ((indx = WWWFindName(www_info, "dispmax")) < 0)
		HitsPerPage = NUM_HITS_PER_PAGE;
	else {
		www_arg = WWWGetValueByIndex(www_info, indx);
		if (isdigit(www_arg[0])) {
                  HitsPerPage = (Int4) atoi(www_arg);
                  if ((HitsPerPage < NUM_HITS_PER_PAGE) || (HitsPerPage > 100)) {
                    printf("Content-type: text/html\n\n");
		    printf("<h2>VASTSERV Error</h2>\n");
		    printf("<h3>Can only display between 20 and 100 neighbors per page!</h3>");
		    exit(1);
                  }
                } 
                else
                  HitsPerPage = NUM_HITS_PER_PAGE;
	}

        if ((indx = WWWFindName(www_info, "sort")) < 0)
		iSort = 0;
	else {
		www_arg = WWWGetValueByIndex(www_info, indx);

		if (isdigit(www_arg[0]))
			iSort = (Int2) atoi(www_arg);
		else
			iSort = 0;
	}

	if ((indx = WWWFindName(www_info, "version")) < 0)
		iFull = ABRIDGED_DISPLAY;
	else {
		www_arg = WWWGetValueByIndex(www_info, indx);

		if (isdigit(www_arg[0]))
			iFull = (Int2) atoi(www_arg);
		else
			iFull = ABRIDGED_DISPLAY;
	}

	/* only two possibilities (so far)! */
	if ((iFull != FULL_DISPLAY) && (iFull != ABRIDGED_DISPLAY))
		iFull = ABRIDGED_DISPLAY;

	if ((indx = WWWFindName(www_info, "hit")) < 0) {
		MakeVastTable(Fsid, pbsa, iSort, subsetnum, pagenum, HitsPerPage, iFull, JobID, pcPass);
		BiostrucAnnotSetFree(pbsa);
		exit(0);
	}

	/* if we get to here then something's wrong! */
	exit(1);

} /* end Main */

