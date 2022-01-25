/* $Id: cddistances.c,v 1.2 2000/07/19 19:59:11 bauer Exp $
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
* File Name:  cddistances.c
*
* Author:  Aron Marchler-Bauer
*
* Initial Version Creation Date: 6/6/2000
*
* $Revision: 1.2 $
*
* File Description: calculate pairwise distances/similarities in a CD
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: cddistances.c,v $
* Revision 1.2  2000/07/19 19:59:11  bauer
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

#define NUMARGS 3
static Args myargs[NUMARGS] = {
  {"Cd-Name",                                                      /*0*/
   NULL, NULL, NULL, FALSE, 'c', ARG_STRING,  0.0, 0, NULL},
  {"Binary output",                                                /*1*/
   "F",   NULL, NULL, FALSE, 'b', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Binary input",                                                 /*2*/
   "F",   NULL, NULL, FALSE, 'e', ARG_BOOLEAN, 0.0, 0, NULL},
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* Main Function                                                             */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Int2 Main()
{
  Char        cCDDid[PATH_MAX];
  CddPtr      pcdd;
  ValNodePtr  vnp;
  TrianglePtr pTri = NULL;
  Int4Ptr     piList;
  ScorePtr    psc;
  Int4        i;

/*---------------------------------------------------------------------------*/
/* Yanli's fix for making binary reading work                                */
/*---------------------------------------------------------------------------*/
  if (!ID1Init()) CddSevError("Unable to initialize ID1");
  objmmdb1AsnLoad();
  objmmdb2AsnLoad();
  objmmdb3AsnLoad();


/*---------------------------------------------------------------------------*/
/* retrieve command-line parameters                                          */
/*---------------------------------------------------------------------------*/
  if (! GetArgs ("cddistances", NUMARGS, myargs)) return (1);

/*---------------------------------------------------------------------------*/
/* assign CD-Id                                                              */
/*---------------------------------------------------------------------------*/
  strcpy(cCDDid,myargs[0].strvalue);
  
/*---------------------------------------------------------------------------*/
/* read Cdd into memory                                                      */
/*---------------------------------------------------------------------------*/
   if (myargs[2].intvalue > 0) {
     pcdd = CddReadFromFile(cCDDid,TRUE);
   } else pcdd = CddReadFromFile(cCDDid,FALSE);
   if (!pcdd) CddSevError("Could not read CDD from disk!");

/*---------------------------------------------------------------------------*/
/* Call the distance calculator                                              */
/*---------------------------------------------------------------------------*/
  pTri = CddCalculateTriangle(pcdd);
/*  piList = CddMostDiverse(pTri,10); */
/*  psc = pTri->scores; i=0;
  while (psc) {
    printf("%8d: ",++i);
    printf("%5.2f  %s\n",psc->value.realvalue, psc->id->str);
    psc = psc->next;
  }
*/

  if (pTri) pcdd->distance = pTri;

/*---------------------------------------------------------------------------*/
/* write Cdd out again                                                       */
/*---------------------------------------------------------------------------*/
   if (myargs[1].intvalue > 0) {
     if (!CddWriteToFile(pcdd,cCDDid,TRUE)) {
       CddSevError("Could not write CD to disk!");
       return(1);
     }
   } else {
     if (!CddWriteToFile(pcdd,cCDDid,FALSE)) {
       CddSevError("Could not write CD to disk!");
       return(1);
     }
   }
   return (0);

/*---------------------------------------------------------------------------*/
/* End of Main                                                               */
/*---------------------------------------------------------------------------*/
}

