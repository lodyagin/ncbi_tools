/* $Id: cddaddtax.c,v 1.2 2000/07/19 19:51:32 bauer Exp $
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
* File Name:  cddaddtax.c
*
* Author:  Aron Marchler-Bauer
*
* Initial Version Creation Date: 3/21/2000
*
* $Revision: 1.2 $
*
* File Description: add taxonomy information to an existing CD
*
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: cddaddtax.c,v $
* Revision 1.2  2000/07/19 19:51:32  bauer
* added modification logging
*
*
* ==========================================================================
*/

#include <ncbi.h>
#include <asn.h>
#include <objseq.h>
#include <objalign.h>
#include "objcdd.h"
#include <cddutil.h>
#include "cdd.h"
#include <taxinc.h>
#include <objsset.h>
#include <tax_cmmn.h>

#define NUMARGS 5
static Args myargs[NUMARGS] = {
  {"Cd-Name",                                                      /*0*/
   NULL, NULL, NULL, FALSE, 'c', ARG_STRING,  0.0, 0, NULL},
  {"Binary output",                                                /*1*/
   "F",  NULL, NULL, FALSE, 'b', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Binary input",                                                 /*2*/
   "F",  NULL, NULL, FALSE, 'e', ARG_BOOLEAN, 0.0, 0, NULL},
  {"New Cd-Name",                                                  /*3*/
   NULL, NULL, NULL, TRUE,  'n', ARG_STRING,  0.0, 0, NULL},
  {"Overwrite existing taxonomy annotation",                       /*4*/
   "F",  NULL, NULL, FALSE, 'o', ARG_BOOLEAN, 0.0, 0, NULL}    
};

static _taxDBCtlPtr Cdd_DBctl = NULL;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* local version of tax1_join                                                */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
static Int4 Cdd_tax1_join(Int4 taxid1, Int4 taxid2)
{
  Int4 lin1[64];
  Int2 i= 1, j;

  lin1[0]= taxid1;
  if(taxid1 > 1) {
    for(i= 1; i < 64; i++) {
      lin1[i]= tax_getParent(Cdd_DBctl, lin1[i-1]);
      if(lin1[i] <= 1) break;
      if(lin1[i] == taxid2) return taxid2;
    }
  }
  while(taxid2 > 1) {
    for(j= 0; j <= i; j++) {
      if(taxid2 == lin1[j]) return taxid2;
    }
    taxid2= tax_getParent(Cdd_DBctl, taxid2);
  }
  return 1;
}


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* Main Function                                                             */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Int2 Main()
{
  Char          cCDDid[PATH_MAX];
  Char          cNEWid[PATH_MAX];
  CddPtr        pcdd;
  ValNodePtr    vnp, pub;
  Char          pcBuf[100];
  CharPtr       pcTest;
  Int4          muid;
  OrgRefPtr     pOrgRef;
  Int4          iTxid1 = -1, iTxid2;
  BioseqSetPtr  bssp;
  SeqEntryPtr   sep;
  BioseqPtr     bsp;
  ValNodePtr    descr, descrLast;
  BioSourcePtr  bsop;
  ObjectIdPtr   oidp;
  DbtagPtr      dbtp;
  Taxon1DataPtr t1dp;
  CharPtr       DB_PATH;
  CharPtr       orgname;
  Int4Ptr       Ids;
  Boolean       bDone = FALSE;

/*---------------------------------------------------------------------------*/
/* Yanli's fix for making binary reading work                                */
/*---------------------------------------------------------------------------*/
  objmmdb1AsnLoad();
  objmmdb2AsnLoad();
  objmmdb3AsnLoad();


/*---------------------------------------------------------------------------*/
/* retrieve command-line parameters                                          */
/*---------------------------------------------------------------------------*/
  if (! GetArgs ("cddaddtax", NUMARGS, myargs)) return (1);

/*---------------------------------------------------------------------------*/
/* assign CD-Id                                                              */
/*---------------------------------------------------------------------------*/
  strcpy(cCDDid,myargs[0].strvalue);
  if (myargs[3].strvalue) {
    strcpy(cNEWid, myargs[3].strvalue);
  } else {
    strcpy(cNEWid, myargs[0].strvalue);
  }
  
/*---------------------------------------------------------------------------*/
/* read Cdd into memory                                                      */
/*---------------------------------------------------------------------------*/
   if (myargs[2].intvalue > 0) {
     pcdd = CddReadFromFile(cCDDid,TRUE);
   } else pcdd = CddReadFromFile(cCDDid,FALSE);
   if (!pcdd) CddSevError("Could not read CDD from disk!");

   if (!pcdd->sequences) CddSevError("Cdd doesn't have sequences!");

   descr = pcdd->description;
   while (descr) {
     if (descr->choice == CddDescr_tax_source) {
       if (myargs[4].intvalue == 0) {
         CddSevError("Cdd already has taxonomy assignments! Aborting ..");
       } else {

       }
     }
     descr = descr->next;
   }


   if (!tax1_init()) CddSevError("Can't connect taxonomy server");
   DB_PATH = getenv("TAXDBPATH");
   Cdd_DBctl = tax_loadDBCtl(DB_PATH,"dbctl.tax");
   tax_loadTree(Cdd_DBctl);

   bssp = (BioseqSetPtr) pcdd->sequences->data.ptrvalue;
   sep = bssp->seq_set;
   while (sep) {
     if (sep->choice != 1) {
       printf("Warning: Error in Seq-entry!\n");
     } else {
       bsp = sep->data.ptrvalue;
       descr = bsp->descr;
       while (descr) {
         if (descr->choice == Seq_descr_source) {
           bsop = descr->data.ptrvalue;
           pOrgRef = bsop->org;
           if (pOrgRef->db) {
             dbtp = pOrgRef->db->data.ptrvalue;
             oidp = dbtp->tag;
             iTxid2 = oidp->id;
           } else {
             iTxid2 = tax1_getTaxId4Str(pOrgRef->taxname,&orgname,&Ids);
           }
           if (iTxid2 >= 1) {
             if (iTxid1 == -1) {
               iTxid1 = iTxid2;
             } else {
               iTxid1 = Cdd_tax1_join(iTxid2, iTxid1);
             }
           }
         }
         descr = descr->next;
       }
     }
     sep = sep->next;
   }
   if (iTxid1 >= 1) {
     t1dp = tax1_getbyid(iTxid1);
     pOrgRef = t1dp->org;
     if (pOrgRef) CddAssignDescr(pcdd,pOrgRef,CddDescr_tax_source,0);   
   }
   tax1_fini();

/*---------------------------------------------------------------------------*/
/* write Cdd out again                                                       */
/*---------------------------------------------------------------------------*/
   if (myargs[1].intvalue > 0) {
     if (!CddWriteToFile(pcdd,cNEWid,TRUE)) {
       CddSevError("Could not write CDD to disk!");
       return(1);
     }
   } else {
     if (!CddWriteToFile(pcdd,cNEWid,FALSE)) {
       CddSevError("Could not write CDD to disk!");
       return(1);
     }
   }
   return (0);

/*---------------------------------------------------------------------------*/
/* End of Main                                                               */
/*---------------------------------------------------------------------------*/
}

