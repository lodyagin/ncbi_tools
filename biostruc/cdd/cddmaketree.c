/* $Id: cddmaketree.c,v 1.2 2000/07/19 20:02:06 bauer Exp $
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
* File Name:  cddmaketree.c
*
* Author:  Aron Marchler-Bauer
*
* Initial Version Creation Date: 1/20/2000
*
* $Revision: 1.2 $
*
* File Description:  example file for CDD api use, creates a compact
*                    Cdd-Tree structure from a Cdd and saves it to disk
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: cddmaketree.c,v $
* Revision 1.2  2000/07/19 20:02:06  bauer
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

#define NUMARGS 4
static Args myargs[NUMARGS] = {
  {"Cd-Name",                                                      /*0*/
   NULL, NULL, NULL, FALSE, 'c', ARG_STRING,  0.0, 0, NULL},
  {"Tree-file extension",                                          /*1*/
   "act",NULL, NULL, FALSE, 't', ARG_STRING,  0.0, 0, NULL} ,
  {"Binary output",                                                /*2*/
   "F",   NULL, NULL, FALSE, 'b', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Binary input",                                                 /*3*/
   "F",   NULL, NULL, FALSE, 'e', ARG_BOOLEAN, 0.0, 0, NULL},
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
  CddTreePtr pcddt;
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
  if (! GetArgs ("cddmaketree", NUMARGS, myargs)) return (1);

/*---------------------------------------------------------------------------*/
/* assign CD-Id                                                              */
/*---------------------------------------------------------------------------*/
  strcpy(cCDDid,myargs[0].strvalue);
  
/*---------------------------------------------------------------------------*/
/* read Cdd into memory                                                      */
/*---------------------------------------------------------------------------*/
   if (myargs[3].intvalue > 0) {
     pcdd = CddReadFromFile(cCDDid,TRUE);
   } else pcdd = CddReadFromFile(cCDDid,FALSE);
   if (!pcdd) CddSevError("Could not read CDD from disk!");

/*---------------------------------------------------------------------------*/
/* assign name for Cdd-tree object                                           */
/*---------------------------------------------------------------------------*/
   strcpy(cNEWid, pcdd->name);
   strcat(cNEWid,".");
   strcat(cNEWid,myargs[1].strvalue);

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

