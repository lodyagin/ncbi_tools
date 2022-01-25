/* $Id: thrdrsmp.c,v 1.2 2000/08/16 21:18:57 hurwitz Exp $
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
* File Name:  thrdrsmp.c
*
* Author:  Stephen Bryant
*
* Initial Version Creation Date: 08/16/2000
*
* $Revision: 1.2 $
*
* File Description: threader
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: thrdrsmp.c,v $
* Revision 1.2  2000/08/16 21:18:57  hurwitz
* fix dangerous warnings found by MS Visual C++, replace rand48 functions with toolkit functions
*
* Revision 1.1  2000/08/16 20:45:21  hurwitz
* initial check in of threading routines
*
* ==========================================================================
*/



/* Choose randomly from a multinomial distribution */

#include <thrdatd.h>
#include <thrddecl.h>
#include <ncbimath.h>
#include <stdio.h>
#include <math.h>

int rsmp(Rnd_Smp* pvl) {
/*----------------------------------------------------*/
/* pvl:  Number and probabilities of parameter values */
/*----------------------------------------------------*/

int	i;		/* The i-th offset parameter value will be the choice */
float	c;		/* Cumulative probabilites across parameter values */
float	r;		/* A uniform random number on the interval 0 - 1 */

/* for(i=0;i<pvl->n;i++) printf("%.4f ",pvl->p[i]); printf("pvl->p\n");  */

/* r=drand48(); c=0.; */

  /* 0 <= RandomNum() <= 2**31 - 1 */
  /* 0x7fffffff = 2**31 - 1 */
  r = (float) ((double)RandomNum()/(double)0x7fffffff);
  c = 0.;

/* printf("r: %.4f\n",r); */

for(i=0; i<pvl->n; i++) {
	c=c+pvl->p[i];
	/* printf("c: %.4f\n",c); */
	if(c>=r) return(i); }

return(i-1);

}
