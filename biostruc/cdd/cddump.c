/* $Id: cddump.c,v 1.8 2000/08/11 19:54:00 hurwitz Exp $
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
* File Name:  cddump.c
*
* Author:  Aron Marchler-Bauer
*
* Initial Version Creation Date: 12/2/1999
*
* $Revision: 1.8 $
*
* File Description: CD-dumper, made from scrap parts of the prototype CDD
*                   server            
*         
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: cddump.c,v $
* Revision 1.8  2000/08/11 19:54:00  hurwitz
* restored CddDenDiagCposComputation and CddCposComputation to original, added CddCposComp which combines the 2
*
* Revision 1.7  2000/08/10 22:41:16  bauer
* replaced ID1init with ID1BiostrucFetchEnable
*
* Revision 1.6  2000/08/09 21:29:08  hurwitz
* adding cddutil.c to VC++ build
*
* Revision 1.5  2000/08/03 22:03:41  bauer
* added support for 3D-structure link highlighting
*
* Revision 1.4  2000/08/01 21:25:11  bauer
* initial changes for consensus
*
* Revision 1.3  2000/07/28 18:00:58  bauer
* fixed typecasts
*
* Revision 1.2  2000/07/19 19:39:17  bauer
* added modification logging
*
*
* ==========================================================================
*/


#include <stdio.h>
#include <ncbi.h>
#include <lsqfetch.h>
#include <netentr.h>
#include <www.h>
#include <sys/resource.h>
#include <asn.h>
#include <accid1.h>
#include <accentr.h>
#include <accutils.h>
#include <mmdbapi.h>
#include <mmdbapi1.h>
#include <objmmdb1.h>
#include <objmmdb2.h>
#include <objmmdb3.h>
#include <objmime.h>
#include <strimprt.h>
#include "objcdd.h"
#include "cddsrv.h"
#include "cdd.h"
#include "cddutil.h"
#include <posit.h>
#include <medutil.h>

#define NUMARGS 10
static Args myargs[NUMARGS] = {
  {"Cd-Name",                                                /*0*/
   "RHO", NULL, NULL, FALSE, 'c', ARG_STRING,  0.0, 0, NULL},
  {"Extension for ASN.1 output file name",                   /*1*/
	 "acd", NULL, NULL, FALSE, 'e', ARG_STRING,  0.0, 0, NULL},
  { "Binary output",                                         /*2*/
   "F",   NULL, NULL, FALSE, 'b', ARG_BOOLEAN, 0.0, 0, NULL},
  { "Use Kludge Checkpoint computations",                    /*3*/
   "F",   NULL, NULL, FALSE, 'k', ARG_BOOLEAN, 0.0, 0, NULL},
  { "Source Identifier",                                     /*4*/
   NULL, NULL, NULL, FALSE,  's', ARG_STRING,  0.0, 0, NULL},
  { "Convert to multiple Alignment",                         /*5*/
   "F",   NULL, NULL, FALSE, 'm', ARG_BOOLEAN, 0.0, 0, NULL},
  { "File extension for tree file",                          /*6*/
   "act", NULL, NULL, FALSE, 't', ARG_STRING,  0.0, 0, NULL},
  { "Reference file extension",                              /*7*/
   "REF", NULL, NULL, FALSE, 'r', ARG_STRING,  0.0, 0, NULL},
  { "Status flag for CDD",                                   /*8*/
    "2",  NULL, NULL, FALSE, 'f', ARG_INT,     0.0, 0, NULL},
  { "Calculate a consensus sequence",                        /*9*/
   "F",   NULL, NULL, FALSE, 'C', ARG_BOOLEAN, 0.0, 0, NULL}
};


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
  GetAppParam("cdd", "CDDSRV", "CDDpath", "", CDDpath, PATH_MAX);
  if (CDDpath[0] == '\0') {
    ErrPostEx(SEV_FATAL,0,0,"CDD config file\nCDDSRV section has no VAST html path...\n");
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
  GetAppParam("cdd", "CDDSRV", "REFpath", "", REFpath, PATH_MAX);
  if (database[0] == '\0') {
    ErrPostEx(SEV_FATAL,0,0,"CDD config file\nCDDSRV section has no REFpath...\n");
    return FALSE;
  }
  return TRUE;
}                                                       /* end GetVastParams */

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* Check whether a particular gi has occured previously in the CddSum lnklst */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
static Boolean UniqueUid(Int4 uid, Boolean bIsPdb, CddSumPtr pcds) {
  
  CddSumPtr    pcdsThis;

  if (!pcds) return TRUE;
  pcdsThis = pcds;
  while (pcdsThis) {
    if (pcdsThis->uid == uid && pcdsThis->bIsPdb == bIsPdb) return FALSE;
    pcdsThis = pcdsThis->next;
  }
  return TRUE;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* Stolen from mmdbsrv.c - converts PDB-Id to numerical mmdb-id              */
/* modified - does no longer check whether the string corresponds to an      */
/* integer in the first place, so that PDB-Id's like "1914" can be read as   */
/* well ..                                                                   */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
static Int4 ConvertMMDBUID(CharPtr pcString)  
{
  Int4    iUID;
  CharPtr pcTemp = NULL;
        
  if (pcString == NULL) return 0;
  iUID = 0;
  pcTemp = StringSave(pcString);
  CleanSpaces(pcTemp);
	iUID = MMDBEvalPDB(pcTemp);
  MemFree(pcTemp);
  return iUID; 
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* allocate a new CddSum linked list entry                                   */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
static CddSumPtr CddSumNew()
{
  CddSumPtr   pcds;
  
  pcds=(CddSumPtr)MemNew(sizeof(CddSum));
  if (pcds==NULL) return(pcds);
  pcds->bIsPdb      = FALSE;
  pcds->bIsMaster   = FALSE;
  pcds->cPdbId[0]   = '\0';
  pcds->cChainId[0] = '\0';
  pcds->cPKBMDom[0] = '\0';
  pcds->cPKBDom[0]  = '\0';
  pcds->cDefLine[0] = '\0';
  pcds->iFsid       = -1;
  pcds->iFid        = -1;
  pcds->iMMDBId     = -1;
  pcds->iCddIdx     = -1;
  pcds->uid         = 0;
  pcds->sip         = NULL;
  pcds->next        = NULL;
  return pcds;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* Free a CddSum linked list                                                 */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
static CddSumPtr CddSumFree(CddSumPtr pcds)
{
  CddSumPtr    next;
  
  while (pcds) {
    next = pcds->next;
    Nlm_MemFree(pcds);
    pcds = next;
  }
  return NULL;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* adds a to a linked list of CddSumPtr, always returns the beginning of the */
/* list and always adds to the end of the list!!                             */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
static CddSumPtr CddSumLink(CddSumPtr PNTR head, CddSumPtr newnode)
{
  CddSumPtr     pcds;
 
  if (head == NULL) return newnode;
  pcds = *head;
  if (pcds != NULL) {
    while(pcds->next != NULL) pcds = pcds->next;
    pcds->next = newnode;
  } else *head = newnode;
  return *head;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* move alignments to PDB-derived sequences up in a list                     */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
static SeqAlignPtr CddAlignSort(SeqAlignPtr salp, CddSumPtr pcds)
{
  CddSumPtr     pcdsThis;
  DenseDiagPtr  ddp;
  DenseSegPtr   dsp;
  SeqIdPtr      sip;
  SeqAlignPtr   salpThis;
  SeqAlignPtr   salpStruct = NULL;
  SeqAlignPtr   salpStructHead = NULL;
  SeqAlignPtr   salpSeq = NULL;
  SeqAlignPtr   salpSeqHead = NULL;
  Boolean       bCont;

  salpThis = salp;
  while (salpThis) {
    ddp = (DenseDiagPtr) salpThis->segs;
    sip = ddp->id;
    sip = sip->next;
    
    bCont = FALSE;
    pcdsThis = pcds;
    while (pcdsThis) {
      if (pcdsThis->bIsPdb && CddSameSip(pcdsThis->sip,sip)) {
        bCont = TRUE; break;
      }
      pcdsThis = pcdsThis->next;
    }
    if (bCont && sip->choice == 15) {
      if (!salpStruct) {
        salpStruct = salpThis; salpStructHead = salpThis;
      } else {
        salpStruct->next = salpThis;
        salpStruct = salpStruct->next;
      }
    } else {
      if (!salpSeq) {
        salpSeq = salpThis; salpSeqHead = salpThis;
      } else {
        salpSeq->next = salpThis;
        salpSeq = salpSeq->next;
      }
    }
    salpThis = salpThis->next;
  }
  if (salpSeqHead) {
    if (!salpStructHead) return(salpSeqHead);
    salpStruct->next = salpSeqHead;
    salpSeq->next = NULL;
  } else {
    salpStruct->next = NULL;
  }
  return(salpStructHead);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* move all pdb-derived sequences up in the list                             */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
static CddSumPtr CddSumSort(CddSumPtr pcds)
{
  CddSumPtr pcdsThis;
  CddSumPtr pcdsStruct = NULL;
  CddSumPtr pcdsStructHead = NULL;
  CddSumPtr pcdsSeq = NULL;
  CddSumPtr pcdsSeqHead = NULL;

  pcdsThis = pcds;
  while (pcdsThis) {
    if (pcdsThis->bIsPdb) {
      if (!pcdsStruct) {
        pcdsStruct = pcdsThis; pcdsStructHead = pcdsThis;
      } else {
        pcdsStruct->next = pcdsThis;
        pcdsStruct = pcdsStruct->next;
      }
    } else {
      if (!pcdsSeq) {
        pcdsSeq = pcdsThis; pcdsSeqHead = pcdsThis;
      } else {
        pcdsSeq->next = pcdsThis;
        pcdsSeq = pcdsSeq->next;
      }
    }
    pcdsThis = pcdsThis->next;
  }
  if (pcdsSeq) {
    pcdsStruct->next = pcdsSeqHead;
    pcdsSeq->next = NULL;
  } else {
    pcdsStruct->next = NULL;
  }
  return(pcdsStructHead);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* Read information about which VAST data have to be retrieved               */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
static void CddReadVASTInfo(CddSumPtr pcds)
{
  FILE             *fp;
  Char             pcBuf[100];
  Char             path[PATH_MAX];
  Char             cSlave[7];
  Char             cMaster[7];
  Char             cPKBMDom[7];
  Char             cPKBDom[9];
  CharPtr          pcTest;
  CddSumPtr        pcdsThis;

  strcpy(path,CDDvpath);
  strcat(path,cCDDid);
  strcat(path,".VASTinfo");
/*  if (!(fp=FileOpen(path, READ))) CddSevError("Can't read VAST information"); */
  if (!(fp=FileOpen(path, READ))) return;
  do {
    pcBuf[0]='\0';
    pcTest = fgets(pcBuf, (size_t)100,fp);
    if (pcTest) {
      strcpy(cMaster,strtok(pcTest,"\t"));  cMaster[6]  = '\0';
      strcpy(cSlave,strtok(NULL,"\t"));     cSlave[6]   = '\0';
      strcpy(cPKBMDom,strtok(NULL,"\t"));   cPKBMDom[6] = '\0';
      strncpy(cPKBDom,strtok(NULL,"\t"),8); cPKBDom[8]  = '\0';
      pcdsThis = pcds;
      while (pcdsThis) {
        if (pcdsThis->bIsPdb && !(pcdsThis->bIsMaster)) {
          if (strncmp(cSlave,pcdsThis->cPdbId,4)==0) {
            if (pcdsThis->cChainId[0]==cSlave[5]) {
              strcpy(pcdsThis->cPKBMDom,cPKBMDom);
              strcpy(pcdsThis->cPKBDom,cPKBDom);
            }
          } 
        }
        pcdsThis = pcdsThis->next;
      }
    }      
  } while (pcTest);
  FileClose(fp);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* stolen from vastlocl.c                                                    */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
static BiostrucAnnotSetPtr CddVASTBsAnnotSetGet (Int4 uid)
{
  AsnIoPtr            aip = NULL;
  AsnTypePtr          atp = NULL;
  Char                path[PATH_MAX];
  Char                compath[PATH_MAX];
  Char                tempfile[PATH_MAX];
  Char                pcId[20];    
  Int2                iFileExists = 0;
  BiostrucAnnotSetPtr pbsa = NULL;
  int                 iAvail = 1;
  FILE                *pipe;

  sprintf(pcId, "%ld", (long) uid);
  path[0] = '\0';
  StringCpy(path, database);
  StringCat(path, pcId);
  StringCat(path, ".bas");

#ifdef MMDB_UNIXCOMPRESSED
  compath[0] = '\0';
  sprintf(compath, "%s -c %s.gz ", gunzip, path);
  pipe = popen(compath, "rb");
  if (pipe == NULL) {
    ErrPostEx(SEV_FATAL,0,0, "VASTBsAnnotSetGet failed: Can't find gunzip in path.\n");
    return NULL;
  }
  aip = AsnIoNew(ASNIO_BIN_IN, pipe, NULL, NULL, NULL);
#else
  iFileExists = FileLength(path);
  if (iFileExists == 0) {
    return NULL;
  }
  aip = AsnIoOpen(path, "rb");
#endif
  if (aip) {
    pbsa = BiostrucAnnotSetAsnRead(aip, NULL);
    AsnIoClose (aip);
  }
#ifdef MMDB_UNIXCOMPRESSED 
   pclose(pipe);
#endif
   if (!pbsa) return NULL;  
   return pbsa;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* stolen from vastlocl.c                                                    */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
static Boolean CddIsVASTData(Int4 uid)
{
  AsnIoPtr   aip = NULL;
  AsnTypePtr atp = NULL;
  Char       path[PATH_MAX];
  Char       pcId[30];

  sprintf(pcId, "%ld", (long) uid);
  path[0] = '\0';
  StringCpy(path, database);
  StringCat(path, pcId);
  StringCat(path, ".bas");

#ifdef MMDB_UNIXCOMPRESSED 
  StringCat(path, ".gz");
  if (FileLength(path) != 0) return TRUE;
#else
  if (FileLength(path) != 0) return TRUE;
#endif
   return FALSE;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* stolen from vastsrv.c                                                     */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
static BiostrucAnnotSetPtr CddLocalGetFeatureSet(Int4 mmdbid, Int4 feature_set_id)
{
  BiostrucAnnotSetPtr   basp = NULL;
  BiostrucAnnotSetPtr   basp2 = NULL;
  BiostrucFeatureSetPtr pbsfs = NULL;
  BiostrucFeatureSetPtr pbsfsLast = NULL;
    
  if (CddIsVASTData(mmdbid))
    basp = (BiostrucAnnotSetPtr) CddVASTBsAnnotSetGet(mmdbid);
  else if (CddIsVASTData(feature_set_id)) {
    basp = (BiostrucAnnotSetPtr) CddVASTBsAnnotSetGet(feature_set_id);
    if (basp != NULL) return basp;
  } 

  if (basp == NULL) return NULL;
 
  pbsfs = basp->features;
  pbsfsLast  = NULL;
  basp2 = NULL;
  while (pbsfs) {
    if (pbsfs->id == feature_set_id) {
      basp2 = BiostrucAnnotSetNew();
      basp2->id = basp->id;
      basp->id = NULL; /* unlink the id valnode from basp object */
      basp2->descr = basp->descr; 
      basp->descr = NULL;  /* unlink the descr from basp object */
      basp2->features = pbsfs;
      if (pbsfsLast) /* relink next to prev */
        pbsfsLast->next = pbsfs->next;
      else basp->features = pbsfs->next;         
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

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* stolen from vastsrv.c                                                     */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
static BiostrucAnnotSetPtr CddBiostrucAnnotSetGetByFid (BiostrucAnnotSetPtr basp, Int4 feature_id, Int4 feature_set_id)
{
  BiostrucAnnotSetPtr   basp2 = NULL;
  BiostrucFeatureSetPtr pbsfs = NULL;
  BiostrucFeaturePtr    pbsf = NULL;

  if (basp == NULL) return NULL;
 
  pbsfs = basp->features;
  while (pbsfs) {
    if (pbsfs->id == feature_set_id) {
      pbsf =  pbsfs->features;
      while(pbsf) {
        if (pbsf->id == feature_id) {  /* found it */
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

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* From the VAST info get the feature set id's required for slaves           */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
static void CddDetermineFsids(CddSumPtr pcds, PDNMS pModelStruc)
{
  CddSumPtr       pcdsThis;
  PDNMG           pdnmg;
  PDNMM           pdnmm;
  PMGD            pmgd;
  PMMD            pmmd;
  PMSD            pmsd;
  Int4            uid;
  Int4            ichn, idom, domcnt;
  Int4            chndom0, chndom;
  Int2            cnt;
  Char            cDid[3];

  pmsd = (PMSD) pModelStruc->data.ptrvalue;
  uid = (Int4) pmsd->iMMDBid;

  for (pdnmm=pmsd->pdnmmHead,cnt=ichn=0;pdnmm!=NULL;pdnmm=pdnmm->next) {
    pmmd = (PMMD) pdnmm->data.ptrvalue;
/*---------------------------------------------------------------------------*/
/* getting the chain id from pmmd seems to be a safer way than to increment  */
/*---------------------------------------------------------------------------*/
    ichn = pmmd->iChainId;
    if ((pmmd->bWhat) & AM_PROT) {
      if ((pmmd->iResCount <= 1) || (pmmd->iGi <= 0)) continue;
      chndom0 = 10000 * uid + 100 * (Int4) ichn;
      pcdsThis = pcds;
      while (pcdsThis) {
        if (pcdsThis->bIsPdb && !pcdsThis->bIsMaster) {
          if (pcdsThis->cPKBMDom[4]==pmmd->pcMolName[0]) { 
            if (pcdsThis->cPKBMDom[5]=='0') pcdsThis->iFsid = chndom0;
          }
        }
        pcdsThis = pcdsThis->next;
      }
      idom = 0; domcnt = 0;
      for (pdnmg=pmmd->pdnmgHead; pdnmg != NULL; pdnmg=pdnmg->next) {
        pmgd = pdnmg->data.ptrvalue;
        if (pmgd->iDomain > idom) {
          idom = (Int4) pmgd->iDomain;
          chndom = chndom0+idom;
          domcnt++;
          pcdsThis = pcds;
          while (pcdsThis) {
            if (pcdsThis->bIsPdb && !pcdsThis->bIsMaster) {
              if (pcdsThis->cPKBMDom[4]==pmmd->pcMolName[0]) { 
                strcpy(cDid,&(pcdsThis->cPKBMDom[5]));
                if (atoi(cDid)==domcnt) pcdsThis->iFsid = chndom;
              }
            }
            pcdsThis = pcdsThis->next;
          }
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* get the Feature id's identifying the correct slaves in the VAST data      */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
static void CddDetermineFids(CddSumPtr pcds, PDNMS pModelStruc, CharPtr szName)
{
  CddSumPtr       pcdsThis;
  PDNMG           pdnmg;
  PDNMM           pdnmm;
  PMGD            pmgd;
  PMMD            pmmd;
  PMSD            pmsd;
  Int4            uid;
  Int4            ichn, idom, domcnt;
  Int4            chndom0, chndom;
  Int2            cnt;
  Char            cDid[3];

  pmsd = (PMSD) pModelStruc->data.ptrvalue;
  uid = (Int4) pmsd->iMMDBid;

  for (pdnmm=pmsd->pdnmmHead,cnt=ichn=0;pdnmm!=NULL;pdnmm=pdnmm->next) {
    pmmd = (PMMD) pdnmm->data.ptrvalue;
    ichn = pmmd->iChainId;
    if ((pmmd->bWhat) & AM_PROT) {
      if ((pmmd->iResCount <= 1) || (pmmd->iGi <= 0)) continue;
      chndom0 = 100000 * uid + 1000 * (Int4) ichn;
      pcdsThis = pcds;
      while (pcdsThis) {
        if (pcdsThis->bIsPdb && !pcdsThis->bIsMaster) {
          if (strncmp(pcdsThis->cPdbId,szName,4)==0) {
            if (pcdsThis->cPKBDom[5]==pmmd->pcMolName[0]) { 
              if (pcdsThis->cPKBDom[7]==' ') {
                chndom = chndom0 + 1;
                if (pcdsThis->iFid == -1) pcdsThis->iFid = chndom;
              }
            }
          }
        }
        pcdsThis = pcdsThis->next;
      }
      idom = 0; domcnt = 0;
      for (pdnmg=pmmd->pdnmgHead; pdnmg != NULL; pdnmg=pdnmg->next) {
        pmgd = pdnmg->data.ptrvalue;
        if (pmgd->iDomain > idom) {
          idom = (Int4) pmgd->iDomain;
          chndom = chndom0 + idom * 10;
          domcnt++;
          pcdsThis = pcds;
          while (pcdsThis) {
            if (pcdsThis->bIsPdb && !pcdsThis->bIsMaster) {
              if (strncmp(pcdsThis->cPdbId,szName,4)==0) {
                if (pcdsThis->cPKBDom[5]==pmmd->pcMolName[0]) { 
                  strcpy(cDid,&(pcdsThis->cPKBDom[7]));
                  chndom = chndom + 1;
                  if (atoi(cDid)==domcnt) {
                     if (pcdsThis->iFid == -1) pcdsThis->iFid = chndom;
                  }
                }
              }
            }
            pcdsThis = pcdsThis->next;
          }
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* Read the table of Cdd Names and descriptions                              */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
static CddDescPtr CddReadDescr() {
  FILE             *fp;
  Char              pcBuf[CDD_MAX_DESCR];
  CharPtr           pcTest;
  CddDescPtr        pcdd = NULL;
  CddDescPtr        pcddThis;

  if (!(fp = FileOpen(CDDdescr, READ))) 
    CddSevError("Can not read description file!");
  do {
    pcBuf[0]='\0';
    pcTest = fgets(pcBuf, (size_t)CDD_MAX_DESCR, fp);
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
/* Main Function the name of the CD to be dumped is a command-line parameter */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Int2 Main()
{
  AsnIoPtr                 paiFile, aipr, aip = NULL;
  BioseqPtr                bsp, bspMaster, bspTrunc, bspFake, bspNew;
  BioseqSetPtr             bssp;
  BiostrucAlignPtr         pbsaStruct;
  BiostrucAlignSeqPtr      pbsaSeq;
  BiostrucAnnotSetPtr      pbsa = NULL;
  BiostrucAnnotSetPtr      pbsaShort = NULL;
  BiostrucFeaturePtr       pbsf;
  BiostrucIdPtr            pbsi;
  BiostrucPtr              pbsMaster, pbsSlave, pbsSlaveHead = NULL, pbsSlaveTail;
  BiostrucPtr              pbsTemp, pbsXtra;
  BiostrucSeqsPtr          pbsaStrSeq;
  CddDescPtr               pcddesc;
  CddPtr                   pcdd;
  CddSumPtr                pcds = NULL;
  CddSumPtr                pcdsThis = NULL;
  CddTreePtr               pcddt;
  CharPtr                  Name, pcPDB, pcTest;
  CharPtr                  cCategory;
  ChemGraphAlignmentPtr    cgap;
  ChemGraphPntrsPtr        cgpp;
  DbtagPtr                 dbtp;
  DenseDiagPtr             ddp;
  GlobalIdPtr              pGid;
  ObjectIdPtr              oidp;
  PDBSeqIdPtr              pdb_seq_id;
  ResidueIntervalPntrPtr   masterrip, slaverip, mrip, srip, mtailrip, stailrip;
  SeqAlignPtr              salpHead, salpTail, salpCopy, salpNew;
  SeqAnnotPtr              psaCAlignHead = NULL;
  SeqAnnotPtr              psaOAlignHead = NULL;
  SeqAnnotPtr              psaVAlignHead = NULL;
  SeqAnnotPtr              psaVAlignTail = NULL;
  SeqEntryPtr              sep, sepNew;
  SeqEntryPtr              sepThis;
  SeqIdPtr                 sip, mastersip, slavesip, pdbsip, sip_master;
  SeqIdPtr                 sipNew, sipFake;
  SeqIntPtr                sintp;
  TextSeqIdPtr             tsip;
  TrianglePtr              pTri = NULL;
  ValNodePtr               pvnAlignment, vnpThis, pub;
  ValNodePtr               pvnId, vnp, desc = NULL;
  PDNMS                    pModelStruc;
  Boolean                  is_network;
  Boolean                  bWriteOK = FALSE;
  Boolean                  bAll  = TRUE;
  Boolean                  bChain = TRUE, bRemove;
  Boolean                  bShrunk = FALSE;
  Int2                     iSeqStrMode = CDDSEVSTRUC, iPDB = 1, iOal = 0;
  Int2                     retcode = 2;
  Int4                     uidmaster = 0, uid;
  Int4                     iCddSize = 0;
  Int4                     Gi, nGi = 0, nPdb = 0;
  Int4                     iModelComplexity = ONECOORDRES;
  Int4                     iMMDBid, iDomain;
  Int4                     i, iMmid, iSmid;
  Int4                     iPcount, muid;
  Char                     chain[2], cChain;
  Char                     CDDalign[PATH_MAX];
  Char                     cOutFile[PATH_MAX];
  Char                     CDDref[PATH_MAX];
  Char                     szName[5], pcBuf[100];
  Char                     cLink[23] = "linked to 3D-structure";
  FILE                     *fp;
  
  cLink[22] = '\0';

/*---------------------------------------------------------------------------*/
/* retrieve command-line parameters                                          */
/*---------------------------------------------------------------------------*/
  if (! GetArgs ("cddump", NUMARGS, myargs)) return (1);

/*---------------------------------------------------------------------------*/
/* assign CD-Id                                                              */
/*---------------------------------------------------------------------------*/
  strcpy(cCDDid,myargs[0].strvalue);

/*---------------------------------------------------------------------------*/
/* retrieve names for directories etc.                                       */
/*---------------------------------------------------------------------------*/
  if (!CddGetParams()) CddSevError("Couldn't read from config file...");

/*---------------------------------------------------------------------------*/
/* initialize Entrez interface - this is needed to retrieve sequences        */
/*---------------------------------------------------------------------------*/
/*  if (!ID1Init()) */
  if (!ID1BioseqFetchEnable("cddump",TRUE))
    CddSevError("Unable to initialize ID1");
  if ( !EntrezInit("Cn3D", FALSE, &is_network))
    CddSevError("Unable to start Entrez");
    
/*---------------------------------------------------------------------------*/
/* Initialize the data structures needed to collect alignments and sequences */
/*---------------------------------------------------------------------------*/
  pbsaSeq    = BiostrucAlignSeqNew();
  pbsaStrSeq = BiostrucSeqsNew();
  pbsaStruct = BiostrucAlignNew();

/*---------------------------------------------------------------------------*/
/* access the information contained in the master CDD alignment, this is     */
/* the number of sequences, whether these are PDB-derived or not, and their  */
/* PDB-Id's/GI's. No need to continue if this file can not be found/accessed */
/*---------------------------------------------------------------------------*/
  strcpy(CDDalign,CDDdpath);
  strcat(CDDalign,cCDDid);
  strcat(CDDalign,CDDextens);
  aip = AsnIoOpen(CDDalign,"r");
  psaCAlignHead = SeqAnnotAsnRead(aip, NULL);
  AsnIoClose(aip);
  if (!psaCAlignHead) CddSevError("Could not access CDD alignment!");

  salpHead = (SeqAlignPtr) psaCAlignHead->data;
  salpCopy = NULL;
  salpTail = NULL;

  sip_master = NULL;

  while (salpHead) {
    iPcount = 1;
    ddp = (DenseDiagPtr) salpHead->segs;
    sip = ddp->id;
    while (sip != NULL) {
      uid = ID1FindSeqId(sip);
#ifdef DEBUG
      printf(" DEBUG: Found gi:%d\n",uid);
#endif
/*---------------------------------------------------------------------------*/
/* skip processing the alignment if sequence id is not found in Entrez!      */
/* this effectively removes "outdated" sequence entries from what is         */
/* listed in and displayed as the CDD                                        */
/* now the processing should also allow alignments between masters & slaves  */
/* having the same id                                                        */
/*---------------------------------------------------------------------------*/
      if (uid) {
        pcdsThis = CddSumNew();
        pcdsThis->uid = uid;
        pcdsThis->iCddIdx = iCddSize;
        pcdsThis->sip = sip;
        if (sip->choice == 15) {
          pcdsThis->bIsPdb = TRUE;
          pdb_seq_id = (PDBSeqIdPtr) CddGetPdbSeqId(sip);
          pcdsThis->cChainId[0] = pdb_seq_id->chain;
          pcdsThis->cChainId[1] = '\0';
          pcdsThis->cPdbId[0] = pdb_seq_id->mol[0];
          pcdsThis->cPdbId[1] = pdb_seq_id->mol[1];
          pcdsThis->cPdbId[2] = pdb_seq_id->mol[2];
          pcdsThis->cPdbId[3] = pdb_seq_id->mol[3];
          pcdsThis->cPdbId[4] = '\0';
        } else pcdsThis->bIsPdb = FALSE;
        if (!uidmaster && iPcount == 1) {
          if (!sip_master) sip_master = sip;
          uidmaster = uid;
          pcdsThis->bIsMaster = TRUE;
          if (pcdsThis->bIsPdb) nPdb++;
          iCddSize++;
          if (UniqueUid(uid,pcdsThis->bIsPdb,pcds)) {
            sep = ID1SeqEntryGet(uid,retcode);
            if (sep == NULL) CddSevError("Unable to get MasterSeqEntry from Entrez");
            bspNew = CddExtractBioseq(sep,sip);
	    /* bspNew = BioseqLockById(sip); */
            sepNew = SeqEntryNew();
            sepNew->data.ptrvalue = bspNew;
            sepNew->choice = 1;
            ValNodeLink(&(pbsaSeq->sequences), sepNew);
            SeqEntryFree(sep);
          }
          CddSumLink(&(pcds),pcdsThis);
        } else if (uid == uidmaster && iPcount == 1) {
          sip = sip->next;
          iPcount++;
          CddSumFree(pcdsThis);
        } else {
          sep = ID1SeqEntryGet(uid,retcode);
          if (sep == NULL) {
            printf("Unable to get SeqEntry %d from ID1, skipping ..\n",uid);
          } else {
            if (pcdsThis->bIsPdb) nPdb++;
            iCddSize++;
            if (!salpCopy) {
              salpCopy = salpHead;
              salpTail = salpCopy;
            } else {
              salpTail->next = salpHead;
              salpTail = salpTail->next;
            }
            if (UniqueUid(uid,pcdsThis->bIsPdb,pcds)) {
              bspNew = CddExtractBioseq(sep,sip);
	      /* bspNew = BioseqLockById(sip); */
              sepNew = SeqEntryNew();
              sepNew->data.ptrvalue = bspNew;
              sepNew->choice = 1;
              ValNodeLink(&(pbsaSeq->sequences), sepNew);
            }
            SeqEntryFree(sep);
          }
          sip = sip->next;
          CddSumLink(&(pcds),pcdsThis);
        }
      }
      else {
        printf("Warning: %s - could not find uid #%d[%d] in ID1, removed from Cdd\n",cCDDid,iCddSize,sip->data.intvalue);
        sip=sip->next;
        iPcount++;
      }
    }
    salpHead=salpHead->next;
  }
  pbsaStruct->sequences = pbsaSeq->sequences;
  pbsaStrSeq->sequences = pbsaSeq->sequences;

/*---------------------------------------------------------------------------*/
/* the purged list of alignments must point to nothing after the last member */
/* and alignments must be sorted so that all pdb-derived slaves appear first */
/*---------------------------------------------------------------------------*/
  salpTail->next = NULL;
  if (nPdb > 1) salpCopy = CddAlignSort(salpCopy,pcds);
  psaCAlignHead->data = salpCopy;

/*---------------------------------------------------------------------------*/
/* now sort the pcds linked list so that structures appear in first place    */
/*---------------------------------------------------------------------------*/
  if (nPdb > 1) {
    pcds = CddSumSort(pcds);
  } else if (nPdb == 1) {
    iSeqStrMode = CDDONESTRUC;
  } else {
    iSeqStrMode = CDDSEQUONLY;
  }

/*---------------------------------------------------------------------------*/
/* if more than one structure present, VAST results have to be retrieved     */
/*---------------------------------------------------------------------------*/
  if (nPdb > 1) {
    objmmdb1AsnLoad();
    objmmdb2AsnLoad();
    objmmdb3AsnLoad();
    OpenMMDBAPI((POWER_VIEW /* ^ FETCH_ENTREZ */), NULL);
/*---------------------------------------------------------------------------*/
/* read the information contained in the VAST info file. This is used to     */
/* figure out which biostruc annot sets have to be retrieved from the VAST   */
/* data base                                                                 */
/*---------------------------------------------------------------------------*/
    CddReadVASTInfo(pcds);
/*---------------------------------------------------------------------------*/
/* Now identify the feature set id's that go with each structurally aligned  */
/* slave - these need not be the same                                        */
/*---------------------------------------------------------------------------*/
    iMMDBid = ConvertMMDBUID(pcds->cPdbId); pcds->iMMDBId = iMMDBid;
#ifdef DEBUG
    printf(" DEBUG: MMDB-id of master is: %d\n",iMMDBid);
#endif
/*---------------------------------------------------------------------------*/
/* MakeAModelstruc sets the chemical_graph ptr to NULL, kludge gets the bio- */
/* struc again                                                               */
/*---------------------------------------------------------------------------*/
    strcpy(szName,pcds->cPdbId);
    pbsXtra = FetchBiostrucPDB(szName,iModelComplexity,1);
    if (bChain) {
      strcpy(chain,pcds->cChainId);
      if (chain[0] != ' ') {
        pbsTemp = (BiostrucPtr)PruneBiostruc(pbsXtra,chain);
        pbsXtra = NULL;
        pbsXtra = pbsTemp;
      }
    }
    pModelStruc = MakeAModelstruc(pbsXtra);
    CddDetermineFsids(pcds, pModelStruc);
    ClearStructures();
/*    FreeAModelstruc(pModelStruc);
    BiostrucFree(pbsXtra); */
/*---------------------------------------------------------------------------*/
/* identify feature id's that go with each aligned slave                     */
/*---------------------------------------------------------------------------*/
    pcdsThis = pcds->next;
    while (pcdsThis) {
      if (pcdsThis->bIsPdb) {
        strcpy(szName,pcdsThis->cPdbId);
        pbsXtra = FetchBiostrucPDB(szName,iModelComplexity,1);
        if (bChain) {
          strcpy(chain,pcdsThis->cChainId);
          if (chain[0] != ' ') {
            pbsTemp = (BiostrucPtr)PruneBiostruc(pbsXtra,chain);
             pbsXtra = NULL;
            pbsXtra = pbsTemp;
          }
        }
        iMMDBid = ConvertMMDBUID(pcdsThis->cPdbId);
        pcdsThis->iMMDBId = iMMDBid;
        pModelStruc = MakeAModelstruc(pbsXtra);
        if (pModelStruc) CddDetermineFids(pcds,pModelStruc,szName);
        ClearStructures();
/*        FreeAModelstruc(pModelStruc);
        BiostrucFree(pbsXtra); */
      }
      pcdsThis = pcdsThis->next;
    }
/*---------------------------------------------------------------------------*/
/* load the biostruc annot sets (and trim them) for each aligned slave struc.*/
/*---------------------------------------------------------------------------*/
    pcdsThis = pcds;
    while (pcdsThis) {
      if (pcdsThis->bIsPdb && !pcdsThis->bIsMaster) {
        pbsa = (BiostrucAnnotSetPtr) CddLocalGetFeatureSet(pcds->iMMDBId,pcdsThis->iFsid);
        pcdsThis->pbsaShort = CddBiostrucAnnotSetGetByFid(pbsa, pcdsThis->iFid, pcdsThis->iFsid);
        if (pcdsThis->pbsaShort) {
          if (pbsaShort == NULL) {
            pbsaShort = pcdsThis->pbsaShort;
            pbsf = pbsaShort->features->features;
            pbsf->next = NULL;
          } else {
            pbsf->next = pcdsThis->pbsaShort->features->features;
            pbsf = pbsf->next;
            pbsf->next = NULL;
          }
        } else {
          bShrunk = TRUE;
          pcdsThis->bIsPdb = FALSE;
          nPdb--;
/*---------------------------------------------------------------------------*/
/* if a particular slave is not a VAST neighbor, the slave structure must be */
/* removed - otherwise Cn3D will barf                                        */
/*---------------------------------------------------------------------------*/
          iMMDBid = ConvertMMDBUID(pcdsThis->cPdbId);
          pbsSlaveHead = pbsaStruct->slaves;
          pbsSlave = pbsSlaveHead;
          pbsSlaveTail = NULL;
          while (pbsSlaveHead) {
            bRemove = FALSE;
            pbsi=pbsSlaveHead->id;
            while (pbsi) {
              if (pbsi->choice == BiostrucId_mmdb_id) {
                if (pbsi->data.intvalue == iMMDBid) {
                  bRemove = TRUE;
                  if (!pbsSlaveTail) {
                    pbsSlave = pbsSlaveHead->next;
                    pbsSlaveTail = pbsSlave;
                  } else {
                    pbsSlaveTail->next = pbsSlaveHead->next;
                    pbsSlaveTail = pbsSlaveTail->next;
                  }

                  break;
                }
              }
              pbsi = pbsi->next;
            }
            if (!bRemove) pbsSlaveTail = pbsSlaveHead;
            pbsSlaveHead = pbsSlaveHead->next;
          }
          pbsaStruct->slaves = pbsSlave;
        }
      }
      pcdsThis = pcdsThis->next;
      
    }
/*---------------------------------------------------------------------------*/
/* need to reorder the SeqAligns if some PDB-derived seqs. are not VAST ngb. */
/*---------------------------------------------------------------------------*/
    if (bShrunk) {
      if (nPdb > 1) {
        salpCopy = CddAlignSort(psaCAlignHead->data,pcds);
        psaCAlignHead->data = salpCopy;
        if (iSeqStrMode == SEVSTRUC) {
          salpCopy = CddAlignSort(psaOAlignHead->data,pcds);
          psaOAlignHead->data = salpCopy;
        }
      }
    }
    if (pbsaShort) pbsaStruct->alignments = pbsaShort;  
    if (nPdb == 1) {
      if (iSeqStrMode == SEVSTRUC) iSeqStrMode = ONESTRUC;
      if (iSeqStrMode == CDDSEVSTRUC) iSeqStrMode = CDDONESTRUC;
    }

/*---------------------------------------------------------------------------*/
/* if required, change structure alignment to reflect SMART sequence align-  */
/*  ment, i.e. "fix coloring" in the Cn3D display                            */
/*---------------------------------------------------------------------------*/
    if (pbsaStruct->alignments) {
      pbsf = pbsaStruct->alignments->features->features;
      while (pbsf) {
        pcPDB=StringSave(PDBNAME_DEFAULT);
        iDomain = 0; cChain = '-';
        pcPDB[0] = pbsf->name[0]; pcPDB[1] = pbsf->name[1];
        pcPDB[2] = pbsf->name[2]; pcPDB[3] = pbsf->name[3];
        cChain = pbsf->name[4];
        iDomain = atoi ((char *) &pbsf->name[5]);
        mastersip = MakePDBSeqId2(pcPDB,cChain,iDomain,FALSE);
        pdbsip = MakePDBSeqId2(pcPDB,cChain,iDomain,FALSE);
        pcPDB[0] = pbsf->name[7]; pcPDB[1] = pbsf->name[8];
        pcPDB[2] = pbsf->name[9]; pcPDB[3] = pbsf->name[10];
        cChain = pbsf->name[11];
        iDomain = atoi ((char *) &pbsf->name[12]);
        slavesip = MakePDBSeqId2(pcPDB,cChain,iDomain,FALSE);
        MemFree(pcPDB);
        pvnAlignment = ValNodeFindNext(pbsf->Location_location,NULL,Location_location_alignment);
        cgap = (ChemGraphAlignmentPtr) pvnAlignment->data.ptrvalue;
        cgpp = (ChemGraphPntrsPtr) cgap->alignment->data.ptrvalue;
        masterrip = (ResidueIntervalPntrPtr) cgpp->data.ptrvalue;
        iMmid = masterrip->molecule_id;
        cgpp = (ChemGraphPntrsPtr) cgap->alignment->next->data.ptrvalue;
        slaverip = (ResidueIntervalPntrPtr) cgpp->data.ptrvalue;
        iSmid = slaverip->molecule_id;
/*---------------------------------------------------------------------------*/
/* now identify the CDD alignment corresponding to the current structure ali.*/
/*---------------------------------------------------------------------------*/
        salpHead = (SeqAlignPtr) psaCAlignHead->data;
        while (salpHead) {
          ddp = (DenseDiagPtr) salpHead->segs;
          sip = ddp->id;
          if (CddSameSip(sip,mastersip)) {
            sip = sip->next; if (CddSameSip(sip,slavesip)) {
/*---------------------------------------------------------------------------*/
/* if the alignment is found, start building up a new set of residue interv. */
/*---------------------------------------------------------------------------*/
              masterrip = NULL; slaverip = NULL;
              while (ddp) {
                mrip = ResidueIntervalPntrNew();
                srip = ResidueIntervalPntrNew();
                mrip->from = ddp->starts[0]+1;
                srip->from = ddp->starts[1]+1;
                mrip->to = ddp->starts[0]+ddp->len;
                srip->to = ddp->starts[1]+ddp->len;
                mrip->molecule_id = iMmid;
                srip->molecule_id = iSmid;
                if (!masterrip) {
                  masterrip = mrip;
                  slaverip = srip;
                  mtailrip = masterrip;
                  stailrip = slaverip;
                } else {
                  mtailrip->next = mrip;
                  stailrip->next = srip;
                  mtailrip = mtailrip->next;
                  stailrip = stailrip->next;
                }
                ddp = ddp->next;
              }

            }
          }
          salpHead = salpHead->next;
        }
        cgpp = (ChemGraphPntrsPtr) cgap->alignment->data.ptrvalue;
        cgpp->data.ptrvalue = masterrip;
        cgpp = (ChemGraphPntrsPtr) cgap->alignment->next->data.ptrvalue;
        cgpp->data.ptrvalue = slaverip;
        pbsf = pbsf->next;
      }
    }
  }


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* allocate the CDD data structure                                           */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
  pcdd = (CddPtr) CddNew();
  
  pcdd->name = cCDDid;
    
  pvnId = ValNodeNew(NULL);
  pGid = (GlobalIdPtr) GlobalIdNew();
  pGid->accession = cCDDid;
  pvnId->choice = CddId_gid;
  pvnId->data.ptrvalue = (GlobalIdPtr) pGid;
  pcdd->id = pvnId;
  pcdd->seqannot = (SeqAnnotPtr) psaCAlignHead;

/*---------------------------------------------------------------------------*/
/* Assign pointer to bioseqs                                                 */
/*---------------------------------------------------------------------------*/
  bssp = BioseqSetNew();
  bssp->seq_set = pbsaSeq->sequences;

  pcdd->sequences = ValNodeNew(NULL);
  pcdd->sequences->choice = 2;
  pcdd->sequences->data.ptrvalue = bssp;

/*---------------------------------------------------------------------------*/
/* Assign pointer to BiostrucFeatureSet (holding the VAST alignments)        */
/*---------------------------------------------------------------------------*/
  if (nPdb > 1) { 
    pcdd->features = (BiostrucAnnotSetPtr) pbsaStruct->alignments;
  }

/*---------------------------------------------------------------------------*/
/* fill in as much of the descriptive information as available               */
/* use a set of generic calls to CddAssignDescr                              */
/*---------------------------------------------------------------------------*/
  pcddesc = CddReadDescr();
  while (pcddesc) {
    if (strcmp(cCDDid,pcddesc->cCddId)==0) {
      CddAssignDescr(pcdd,(CharPtr) pcddesc->cDescr,CddDescr_comment,0);
      if (strncmp(pcddesc->cSourc," ",1)!=0) {
        cCategory = strdup(pcddesc->cSourc);
        cCategory[strlen(cCategory)-1]='\0';
        CddAssignDescr(pcdd, (CharPtr) cCategory,CddDescr_category,0);
      }
      break;
    }  
    pcddesc = pcddesc->next;
  }
  if (nPdb) CddAssignDescr(pcdd, cLink, CddDescr_comment,0);
  CddAssignDescr(pcdd, NULL, CddDescr_status , myargs[8].intvalue);

/*---------------------------------------------------------------------------*/
/* assign the source identifier if present                                   */
/*---------------------------------------------------------------------------*/
  CddAssignDescr(pcdd,myargs[4].strvalue,CddDescr_source,0);

/*---------------------------------------------------------------------------*/
/* assign create date as the current date                                    */
/*---------------------------------------------------------------------------*/
  CddAssignDescr(pcdd,(DatePtr) DateCurr(),CddDescr_create_date,0);

/*---------------------------------------------------------------------------*/
/* Assign references if they can be found                                    */
/*---------------------------------------------------------------------------*/
  strcpy(CDDref,REFpath);
  strcat(CDDref,cCDDid);
  strcat(CDDref,".");
  strcat(CDDref,myargs[7].strvalue);
  fp = FileOpen(CDDref,"r");
  if (fp) {
    MedArchInit();
    do {
      pcBuf[0]='\0';
      pcTest = fgets(pcBuf, (size_t)100, fp);
      if (pcTest) {
        pcTest[strlen(pcTest)-1]='\0';
        muid = (Int4) atoi(pcTest);
        pub  = FetchPub(muid);
        if (pub) CddAssignDescr(pcdd, pub, CddDescr_reference, 0);
      }
    } while (pcTest);
    FileClose(fp);
    MedArchFini();
  }


/*---------------------------------------------------------------------------*/
/* assuming this is a set of pairwise master-slave dendiag alignments,       */
/* calculate the interval on the master that is aligned                      */
/*---------------------------------------------------------------------------*/
  CddAssignProfileRange(pcdd, sip_master);
/*---------------------------------------------------------------------------*/
/* create a truncated version of the master bioseq which corresponds to the  */
/* interval defined above                                                    */
/*---------------------------------------------------------------------------*/
  tsip = (TextSeqIdPtr) TextSeqIdNew();
  tsip->name = cCDDid;
  oidp = (ObjectIdPtr) ObjectIdNew();
  oidp->str = cCDDid;
  dbtp = DbtagNew();
  dbtp->tag = oidp;
  dbtp->db = myargs[4].strvalue;

  sintp = (SeqIntPtr) pcdd->profile_range;
  sipFake = (SeqIdPtr) ValNodeNew(NULL);
  sipFake->choice = 11;
  sipFake->data.ptrvalue = dbtp;
  bspTrunc = (BioseqPtr) BioseqCopy(sipFake,sip_master,sintp->from,sintp->to,0,FALSE);
  bspFake  = (BioseqPtr) BioseqCopy(NULL,sip_master,sintp->from,sintp->to,0,FALSE);
  vnp = pcdd->description;
  while (vnp) {
    if (vnp->choice == CddDescr_comment) {
      desc = ValNodeNew(NULL);
      desc->choice = Seq_descr_title;
      desc->data.ptrvalue = vnp->data.ptrvalue;
      break;
    }
    vnp = vnp->next;
  }
  if (desc) bspTrunc->descr = desc;
  pcdd->trunc_master = (struct bioseq PNTR) bspTrunc;
  
/*---------------------------------------------------------------------------*/
/* Calculate a consensus sequence and make it the master for the SeqAlign    */
/*---------------------------------------------------------------------------*/
  if (myargs[9].intvalue > 0) {
    salpNew = (SeqAlignPtr) CddConsensus(salpCopy,bssp->seq_set,
                                         bspTrunc,pcdd->profile_range);
  
  
  }

  if (myargs[3].intvalue > 0) {
/*---------------------------------------------------------------------------*/
/* converting to multiple DenseSeg using a kludge converter                  */
/*---------------------------------------------------------------------------*/
    salpCopy = (SeqAlignPtr)CddMSLDenDiagToMSLDenSeg(pcdd->seqannot->data);
/*---------------------------------------------------------------------------*/
/* now the alignment should reflect the new n-terminal offset for the master */
/*---------------------------------------------------------------------------*/
    CddReindexMSLDenSegMaster(salpCopy, sintp->from); 
/*---------------------------------------------------------------------------*/
/* using kludge alignment to calculate position-specific scoring matrix      */
/*---------------------------------------------------------------------------*/
    CddCposComputation(salpCopy, bspTrunc, pcdd);
  } else {
/*---------------------------------------------------------------------------*/
/* or do the pssm calculation directly on the DenseDiag Alignment            */
/*---------------------------------------------------------------------------*/
    salpCopy = (SeqAlignPtr) CddCopyMSLDenDiag(pcdd->seqannot->data);
    CddReindexMSLDenDiagMaster(salpCopy, sintp->from);
    CddDenDiagCposComputation((SeqAlignPtr)salpCopy,bspTrunc,bspFake,pcdd);
  }

/*---------------------------------------------------------------------------*/
/* Calculate pairwise percent identities between the sequences, for alignment*/
/* formatting purposes                                                       */
/*---------------------------------------------------------------------------*/
  pTri = CddCalculateTriangle(pcdd);
  if (pTri) pcdd->distance = pTri;


/*---------------------------------------------------------------------------*/
/* for output, convert the pairwise master-slave densediag to a multiple     */
/* dense-diag alignment if the block-structure is consistent - and if        */
/* specified in the command line                                             */
/*---------------------------------------------------------------------------*/
  if (myargs[5].intvalue > 0) {
    pcdd->seqannot->data = (SeqAlignPtr)CddMSLDenDiagToMULDenDiag(pcdd->seqannot->data);
  }

  strcpy(cOutFile,cCDDid);
  strcat(cOutFile,".");
  strcat(cOutFile,myargs[1].strvalue);

  if (!CddWriteToFile(pcdd,cOutFile,(Boolean) myargs[2].intvalue))
    CddSevError("Could not write ASN.1 output");

/*---------------------------------------------------------------------------*/
/* output the corresponding Cdd tree-file                                    */
/*---------------------------------------------------------------------------*/
  strcpy(cOutFile,cCDDid);
  strcat(cOutFile,".");
  strcat(cOutFile,myargs[6].strvalue);
  
  pcddt              = CddTreeNew();
  pcddt->name        = pcdd->name;
  pcddt->id          = pcdd->id;
  pcddt->description = pcdd->description;

  if (!CddTreeWriteToFile(pcddt,cOutFile,(Boolean) myargs[2].intvalue))
    CddSevError("Could not write CDD-tree");



/*  ID1BioseqFetchDisable(); */
  CloseMMDBAPI();
  MMDBFini();
  VASTFini();
  return 0;
}
