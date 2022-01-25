/* $Id: cddadddesc.c,v 1.3 2000/08/14 18:14:21 bauer Exp $
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
* File Name:  cddadddesc.c
*
* Author:  Aron Marchler-Bauer
*
* Initial Version Creation Date: 1/20/2000
*
* $Revision: 1.3 $
*
* File Description: example file for CDD api use, adds a description line
*                   of type "comment" to an existing CDD   
*
* Usage: cddaddcdesc -c CDD.name -d "Description Text" -b T/F (Binary Out)
*                    -e T/F (Binary in) -n "New Name" 
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: cddadddesc.c,v $
* Revision 1.3  2000/08/14 18:14:21  bauer
* fixed assignment of new name
*
* Revision 1.2  2000/07/19 19:48:46  bauer
* added modification logging
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

#define NUMARGS 6
static Args myargs[NUMARGS] = {
  {"Cd-Name",                                                      /*0*/
   NULL, NULL, NULL, FALSE, 'c', ARG_STRING,  0.0, 0, NULL},
  {"Description to be added",                                      /*1*/
	 NULL, NULL, NULL, TRUE,  'd', ARG_STRING,  0.0, 0, NULL},
  {"Binary output",                                                /*2*/
   "F",  NULL, NULL, FALSE, 'b', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Binary input",                                                 /*3*/
   "F",  NULL, NULL, FALSE, 'e', ARG_BOOLEAN, 0.0, 0, NULL},
  {"New Cd-Name",                                                  /*4*/
   NULL, NULL, NULL, TRUE,  'n', ARG_STRING,  0.0, 0, NULL},
  {"Check whether comment already exists",                         /*5*/
   "T",  NULL, NULL, FALSE, 'x', ARG_BOOLEAN, 0.0, 0, NULL}
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* Main Function                                                             */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Int2 Main()
{
  Char       cCDDid[PATH_MAX];
  Char       cNEWid[PATH_MAX];
  CddPtr     pcdd;
  ValNodePtr vnp;

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
  if (myargs[4].strvalue != NULL) {
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
/* Add description to the Cdd - but only if description has been supplied!   */
/*---------------------------------------------------------------------------*/
   if (myargs[1].strvalue) {
     if (myargs[5].intvalue > 0) {
       vnp = pcdd->description;
       while (vnp) {
         if (vnp->choice == CddDescr_comment) 
           CddSevError("Cdd already has comment added!");
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
   return (0);

/*---------------------------------------------------------------------------*/
/* End of Main                                                               */
/*---------------------------------------------------------------------------*/
}




