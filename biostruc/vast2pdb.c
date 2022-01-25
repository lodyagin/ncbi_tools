/* vast2pdb.c
 *
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
 * File Name: vast2pdb.c
 *
 * Author: Tom Madej
 *
 * Version Creation Date: 03/12/98
 *
 * $Log: vast2pdb.c,v $
 * Revision 6.5  1998/07/17 18:47:44  madej
 * Allow files to be seen with Vast Search.
 *
 * Revision 6.4  1998/06/11  19:10:39  madej
 * Minor cosmetic change.
 *
 * Revision 6.3  1998/05/19  20:22:26  madej
 * Add general WWW routines for running on Sun servers.
 *
 * Revision 6.2  1998/03/30  19:13:41  madej
 * Changes by Ken Addess.
 *
 * Revision 6.1  1998/03/12  17:10:26  madej
 * First official version of vast2pdb.c
 *
 */

#include <stdio.h>
#include <ncbi.h>
#include <accentr.h>
#include <netentr.h>
#include <www.h>
#include <sys/resource.h>
#include <mmdbapi.h>
#include "vastlocl.h"
#include "mmdblocl.h"
#include "mmdbdata.h"
#include "vast2mage.h"
#include "vast2pdb.h"
#include "vastsrv.h"

static char *BaseURL;
static FILE *OutputFile = NULL;
static char OutputName[200];



/* Generate a PDB-format file with the reference protein rotated in the frame
 * of the master.  This file can then be loaded, e.g. into Kinemage, and the
 * structural alignment with the master protein displayed.
 */

Boolean LIBCALL VastToPDB(WWWInfoPtr www_info)
{
  FILE *pFile = NULL;
  FILE *pIn = NULL;
  Char pcBuf[100];
  CharPtr pcTest;
  Int4 GetGi, Fid, Fsid;
  Int4 iFileExists = 0, indx;
  Char pcLine[256];
  CharPtr pcL1 = NULL, www_arg;
  CharPtr JobID = NULL, pcPass;
  BiostrucAnnotSetPtr pbsa = NULL; 
  BiostrucAnnotSetPtr pbsaShort = NULL;
  PDNMS pdnmsMaster = NULL;
  PDNMS pdnmsSlave = NULL;
  Int2 iTest = 0, iPDB = 0, ret;
  AsnIoPtr aip = NULL; 
  Char giBuf[20], URL[200];
  char *IPAddress = getenv("REMOTE_HOST");


 
  if ((indx = WWWFindName(www_info, "uid")) < 0) {
    printf("Content-type: text/html\n\n");
    printf("<h2>VASTSERV (VastToPDB)</h2>\n");
    printf("<h3>No accession (PDB ID) was input - nothing to report.</h3>\n");
    return 0;
  }

  www_arg = WWWGetValueByIndex(www_info, indx);

  if (isdigit(www_arg[0]))
    GetGi = (Int4) atoi(www_arg);
  else {
    printf("Content-type: text/html\n\n");
    printf("<h2>VASTSERV Error (VastToPDB)</h2>\n");
    printf("<h3>Non-numeric MMDB-id input - no results.</h3>\n"); 
    return 0;
  }
  
  /* vsid and pass are to look at alignments from VAST Search */
	if ((indx = WWWFindName(www_info, "vsid")) >= 0) {
		www_arg = WWWGetValueByIndex(www_info, indx);
		JobID = StringSave(www_arg);

		if ((indx = WWWFindName(www_info, "pass")) < 0) {
			printf("Content-type: text/html\n\n");
			printf("<body bgcolor = \"#f0f0f0\">\n");
			printf("<h2>VAST SEARCH</h2>\n");
			printf("<h3>Password required.</h3>\n");
			return 0;
		}
		else {
			www_arg = WWWGetValueByIndex(www_info, indx);
			pcPass = StringSave(www_arg);

			if ((ret = Check_VastSearch_Password(pcPass, JobID)) != 1) {
				if (ret == 2) return 0;
				printf("Content-type: text/html\n\n");
				printf("<body bgcolor = \"#f0f0f0\">\n");
				printf("<h2>VAST SEARCH</h2>\n");
				printf("<h3>Incorrect password.</h3>\n"); 
				return 0;
			}
		}
	}
  if ((indx = WWWFindName(www_info, "hit")) < 0) {
    printf("Content-type: text/html\n\n");
    printf("<body bgcolor = \"#f0f0f0\">\n");
    printf("<br>\n<h2>No alignment was selected!</h2>\n");
    printf("<h3>Please click on a box in the leftmost column of the table.</h3>\n");
    return 0;
  }

  www_arg = WWWGetValueByIndex(www_info, indx);

  if (isdigit(www_arg[0]))
    Fid = (Int4) atoi(www_arg);
  else {
    printf("Content-type: text/html\n\n");
    printf("<h2>VASTSERV Error (VastToPDB)</h2>\n");
    printf("<h3>Non-numeric slave alignment code - no results.</h3>\n"); 
    return 0;
  }

  if ((indx = WWWFindName(www_info, "chaindom")) < 0) {
    printf("Content-type: text/html\n\n");
    printf("<h2>VASTSERV Error (VastToPDB)</h2>\n");
    printf("<h3>No feature set ID (master alignment code) - nothing to report.</h3>\n"); 
    return 0;
  }

  www_arg = WWWGetValueByIndex(www_info, indx);

  if (isdigit(www_arg[0]))
    Fsid = (Int4) atoi(www_arg);
  else {
       printf("Content-type: text/html\n\n");
       printf("<h2>Error</h2>\n");
       printf("<h3>Non-numeric master alignment code - no results</h3>\n"); 
       return 0;
   }

  /* action == 0 indicates MIME; action == 1 is text; action == 2 is save */
  if ((indx = WWWFindName(www_info, "action")) < 0)
    iPDB = 0;
  else {
    www_arg = WWWGetValueByIndex(www_info, indx);

    if (isdigit(www_arg[0]))
      iPDB = (Int4) atoi(www_arg);
    else
      iPDB = 0;
  }

  if (VASTInit() != TRUE) {
    printf("Content-type: text/html\n\n");
    printf("<h2>VASTSERV Error (VastToPDB)</h2>\n");
    printf("<h3>Can't find VAST data on server.\n");
    printf("Contact info@ncbi.nlm.nih.gov</h3>\n");
    return 0;
  }

  OpenMMDBAPI((POWER_VIEW /* ^ FETCH_ENTREZ */), NULL);
  if (JobID == NULL)
    pbsa = VASTBsAnnotSetGet(GetGi);
  else
    pbsa = LocalGetFeatureSet(GetGi, Fsid, JobID);

  if (pbsa == NULL) {
    printf("Content-type: text/html\n\n");
    printf("<h2>VASTSERV Error (VastToPDB)</h2>\n");
    printf("<h3>No master alignment record exists for %ld.</h3>\n", (long) GetGi);
    return 0; 
  }

  pbsaShort =  BiostrucAnnotSetGetByFid(pbsa, Fid, Fsid);

  if (pbsaShort == NULL) {
    printf("Content-type: text/html\n\n");
    printf("<h2>VASTSERV Error (VastToPDB)</h2>\n");
    printf("<h3>Can't find alignment record.</h3>\n");
    return 0; 
  }
   
  if (MMDBInit() == FALSE) {
    printf("Content-type: text/html\n\n");
    printf("<h2>VASTSERV Error (VastToPDB)</h2>\n");
    printf("<h3>MMDBInit Failed.</h3>\n");
    return 0;
  }
 
 /*
   aip = AsnIoOpen("pbsashort", "w");  
   BiostrucAnnotSetAsnWrite(pbsaShort, aip, NULL);
   AsnIoClose(aip);
 */
 
  if (JobID == NULL)
     InstBSAnnotSet(pbsaShort);
  else
     InstBSAnnotSetVS(pbsaShort, JobID);
  pdnmsMaster = GetSelectedModelstruc();
  pdnmsSlave = GetSlaveModelstruc();
  if (!pdnmsMaster)
    {
      printf("Content-type: text/html\n\n");
      printf("<h2>Error</h2>\n");
      printf("<h3>Unable to load structures on server, - no VAST results</h3>\n");
      return 0; 
     }		
 
    								       
  strcpy(OutputName,GetTempName("vast3")); 
  if(!(OutputFile = FileOpen(OutputName,WRITE)))
    {
       printf("Content-type: text/html\n\n");
       printf("<h2>Error</h2>\n");
       printf("Temp File Open Failed At Vast3 WWW-Server<p>\n");
       CloseMMDBAPI();
       MMDBFini();
       VASTFini();
      exit(1);
    }
    
 /* PDB FILE GENERATOR */   

   if (iPDB == 2)
     {
       fprintf(OutputFile, "Content-type: application/octet-stream\n\n");
     }
   else if (iPDB == 1)
     {
       fprintf(OutputFile, "Content-type: text/html\n\n");  
       fprintf(OutputFile, "<HTML><body><pre>\n");
     }
   else
     { /* MIME */
       fprintf(OutputFile,"Content-type: chemical/x-pdb\n\n");
     }
  iTest = WritePDBOneModel(pdnmsSlave, OutputFile,  0);
  
  if (!iTest)
    {
      printf("Content-type: text/html\n\n");
      printf("<h2>Error</h2>\n");
      printf("PDB File write failed on Server.<p>\n");
      RemoveTempFiles();   
      CloseMMDBAPI();
      MMDBFini();
      VASTFini();
      exit(1);
    }
  if (iPDB == 1)
     {
       fprintf(OutputFile, "</pre></body></HTML>\n"); 

     }

  
  fflush(OutputFile);
  if (OutputFile != stdout)
    {
      fclose(OutputFile);
      PrintFile(OutputName);
    }
    
  CloseMMDBAPI();
  MMDBFini();
  VASTFini();
  RemoveTempFiles();   
  return 0;
}

