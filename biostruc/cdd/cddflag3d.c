/* $Id: cddflag3d.c,v 1.1 2000/08/03 22:03:40 bauer Exp $
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
* File Name:  cddflag3d.c
*
* Author:  Aron Marchler-Bauer
*
* Initial Version Creation Date: 8/02/2000
*
* $Revision: 1.1 $
*
* File Description: detects and flags CD's with 3D-structure links  
*
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: cddflag3d.c,v $
* Revision 1.1  2000/08/03 22:03:40  bauer
* added support for 3D-structure link highlighting
*
*
* ==========================================================================
*/

#include <ncbi.h>
#include <asn.h>
#include <objseq.h>
#include <objalign.h>
#include <objcdd.h>
#include <cddutil.h>

#define NUMARGS 7
static Args myargs[NUMARGS] = {
  {"Cd-Name",                                                      /*0*/
   NULL, NULL, NULL, FALSE, 'c', ARG_STRING,  0.0, 0, NULL},
  {"Flag to be added",                                      /*1*/
   "linked to 3D-structure", NULL, NULL, TRUE,  'f', ARG_STRING,  0.0, 0, NULL},
  {"Binary output",                                                /*2*/
   "F",  NULL, NULL, FALSE, 'b', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Binary input",                                                 /*3*/
   "F",  NULL, NULL, FALSE, 'e', ARG_BOOLEAN, 0.0, 0, NULL},
  {"New Cd-Name",                                                  /*4*/
   NULL, NULL, NULL, TRUE,  'n', ARG_STRING,  0.0, 0, NULL},
  {"Check whether flag already exists",                            /*5*/
   "T",  NULL, NULL, FALSE, 'x', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Tree-file extension",                                          /*6*/
   "act",NULL, NULL, FALSE, 't', ARG_STRING,  0.0, 0, NULL}
   
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* Main Function                                                             */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Int2 Main()
{
  Char         cCDDid[PATH_MAX];
  Char         cNEWid[PATH_MAX];
  CharPtr      cCurrDesc;
  CddPtr       pcdd;
  ValNodePtr   vnp;
  SeqAlignPtr  salpHead;
  DenseDiagPtr ddp;
  SeqIdPtr     sip1, sip2;
  CddTreePtr   pcddt;

/*---------------------------------------------------------------------------*/
/* Yanli's fix for making binary reading work                                */
/*---------------------------------------------------------------------------*/
  objmmdb1AsnLoad();
  objmmdb2AsnLoad();
  objmmdb3AsnLoad();


/*---------------------------------------------------------------------------*/
/* retrieve command-line parameters                                          */
/*---------------------------------------------------------------------------*/
  if (! GetArgs ("cddadddesc", NUMARGS, myargs)) return (1);

/*---------------------------------------------------------------------------*/
/* assign CD-Id and Id for output                                            */
/*---------------------------------------------------------------------------*/
  strcpy(cCDDid,myargs[0].strvalue);
  if (myargs[4].strvalue) {
    strcpy(cNEWid, myargs[4].strvalue);
  } else {
    strcpy(cNEWid, myargs[0].strvalue);
  }

/*---------------------------------------------------------------------------*/
/* read Cdd into memory                                                      */
/*---------------------------------------------------------------------------*/
   if (myargs[3].intvalue > 0) {
     pcdd = CddReadFromFile(cCDDid,TRUE);
   } else pcdd = CddReadFromFile(cCDDid,FALSE);
   if (!pcdd) CddSevError("Could not read CDD from disk!");

/*---------------------------------------------------------------------------*/
/* check whether CD has 3D structure!                                        */
/*---------------------------------------------------------------------------*/
  if (!pcdd->seqannot) CddSevError("No SeqAlign!");
  salpHead = (SeqAlignPtr) pcdd->seqannot->data;
  ddp = (DenseDiagPtr) salpHead->segs;
  sip1 = ddp->id;
  sip2 = ddp->id->next;
  if (sip1->choice!=SEQID_PDB && sip2->choice!=SEQID_PDB) CddSevError("No 3D-structure!");

/*---------------------------------------------------------------------------*/
/* Add flag to the Cdd - but only if flat has been supplied and isn't there! */
/*---------------------------------------------------------------------------*/
   if (myargs[1].strvalue) {
     if (myargs[5].intvalue > 0) {
       vnp = pcdd->description;
       while (vnp) {
         if (vnp->choice == CddDescr_comment) {
	   cCurrDesc = (CharPtr) vnp->data.ptrvalue;
	   if (strcmp(cCurrDesc,myargs[1].strvalue) == 0) {
             CddSevError("Cdd already flagged!");
	   }
	 }
         vnp = vnp->next;
       }
     }
     CddAssignDescr(pcdd,myargs[1].strvalue,CddDescr_comment,0);
   }

/*---------------------------------------------------------------------------*/
/* write Cdd out again                                                       */
/*---------------------------------------------------------------------------*/
   if (myargs[2].intvalue > 0) {
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
   
/*---------------------------------------------------------------------------*/
/* assign name for Cdd-tree object                                           */
/*---------------------------------------------------------------------------*/
   strcpy(cNEWid, pcdd->name);
   strcat(cNEWid,".");
   strcat(cNEWid,myargs[6].strvalue);

/*---------------------------------------------------------------------------*/
/* Create Cdd-tree object                                                    */
/*---------------------------------------------------------------------------*/
   pcddt              = CddTreeNew();
   pcddt->name        = pcdd->name;
   pcddt->id          = pcdd->id;
   pcddt->description = pcdd->description;

/*---------------------------------------------------------------------------*/
/* write Cdd out again                                                       */
/*---------------------------------------------------------------------------*/
   if (myargs[2].intvalue > 0) {
     if (!CddTreeWriteToFile(pcddt,cNEWid,TRUE)) {
       CddSevError("Could not write CDD-tree to disk!");
       return(1);
     }
   } else {
     if (!CddTreeWriteToFile(pcddt,cNEWid,FALSE)) {
       CddSevError("Could not write CDD-tree to disk!");
       return(1);
     }
   }
   return (0);

/*---------------------------------------------------------------------------*/
/* End of Main                                                               */
/*---------------------------------------------------------------------------*/
}




