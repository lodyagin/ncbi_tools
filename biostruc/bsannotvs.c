/* bsannotvs.c
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
 * File Name: bsannotvs.c
 *
 * Author: Ken Addess
 *
 * $Log: bsannotvs.c,v $
 * Revision 6.1  1998/07/17 18:51:13  madej
 * Allows Vast Search results to be viewed with RasMol and Kinemage.
 *
 */

/* This is a copy from mmdbapi so that vastsearch alignment can be viewed with */
/* Rasmol and Kinemage                                                         */

#include <stdio.h>
#include <ncbi.h>
#include <accentr.h>
#include <objalign.h>
#include <objseq.h>
#include <objmgr.h>
#include <lsqfetch.h>
#include <netentr.h>
#include <www.h>
#include <sys/resource.h>
#include <mmdbapi.h>
#include <mmdbapi1.h>
#include "vastlocl.h"
#include "mmdblocl.h"
#include "mmdbdata.h"
#include "vast2mage.h"
#include "vastsrv.h"
#include <asnmime.h>
#include <objmime.h>
#include <strimprt.h>

#define VSPATH "/usr/people7/addess/vastsearch/"

Boolean InstBSAnnotSetVS(BiostrucAnnotSetPtr pbsasThis, CharPtr JobID)
{
  Int2 iTest;
  Int4 iMMDBId;
  PDNMS pdnmsThis = NULL;
  PMSD  pmsdThis = NULL;
  BiostrucPtr pbsThis = NULL;
  BiostrucIdPtr pbsidThis = NULL;
  BiostrucFeatureSetPtr pbsfsThis = NULL;
  BiostrucFeaturePtr pbsfThis = NULL;
  PSFS psfsThis = NULL;
  CharPtr szTemp;
  Char szName[5];
  Char AsnPath[PATH_MAX];
  Char AsnName[10];
 
 /* the Feature set either attaches to an in-memory Modelstruc */
 /* or, if it specifies a new Modelstruc - attempts a retrieval */
 
  if (pbsasThis == NULL) return FALSE;

/* grab the id out of the Biostruc-descr */

/*  pbsidThis = ValNodeFindNext(pbsasThis->id,NULL,BiostrucId_mmdb_id); */
/*  if (pbsidThis)   only deal with the first one */
/*    {                                           */
/*      iMMDBId = (Int4) pbsidThis->data.intvalue; */
      /*printf("iMMDBID=%ld\n", iMMDBId); master MMDB_ID*/
/*    } */
  
  szTemp = pbsasThis->features->features->name;  
  szName[0] = szTemp[0];
  szName[1] = szTemp[1];
  szName[2] = szTemp[2];
  szName[3] = szTemp[3];
  szName[4] = '\0';
  
  AsnName[0]='\0';
  StringCpy(AsnName, "/b");
  StringCat(AsnName, szName);
	
  AsnPath[0]='\0';
  StringCpy(AsnPath, VSPATH);
  StringCat(AsnPath, JobID);
  StringCat(AsnPath, AsnName);
   
  pbsThis = FetchBS(AsnPath, 0, ALLSIMPLEMDL, 3, POWER_VIEW);
  /*pbsThis = MMDBBiostrucGet(iMMDBId, ALLSIMPLEMDL, 3);*/
  if (!pbsThis) goto nogomem;
  pdnmsThis =  MakeAModelstruc(pbsThis);
  /* side effect is that this is now the selected modelstruc too */
  if (!pdnmsThis) goto nogomem;  

  pmsdThis = (PMSD) pdnmsThis->data.ptrvalue;

  /****** WARNING: fix this if this code is ever used to recognize Neighbor mode *********/
  SetSelectedModelstruc(pdnmsThis);  /* set it to selected one */
  /* at this point we have annot-set, and its associated Modelstruc */

  pbsfsThis = pbsasThis->features;

  while (pbsfsThis)  /* walk through each feature-set */
    {
    iTest = BiostrucAddFeature(pbsfsThis,pdnmsThis);
    if (iTest < 0) goto nogomem; /* a malloc error */
      if (!iTest) goto nogo;  /* bad feature table error - fatal */
      pbsfsThis = pbsfsThis->next;
    }
  /* given a sucessfully registered Biostruc-annot-set
     we can now free it */
  BiostrucAnnotSetFree(pbsasThis);
  return TRUE;
  
nogomem:
nogo:
  BiostrucAnnotSetFree(pbsasThis);
  return FALSE;
}   
