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
* ===========================================================================*/
/*****************************************************************************

File name: blastasn.h

Author: Tom Madden

Contents: prototypes for functions that produce blast18 ASN.1 from
	BLAST structures.

******************************************************************************/

/* $Revision: 6.1 $ 
* $Log: blastasn.h,v $
* Revision 6.1  1999/03/17 16:59:19  madden
* remove comment in comment
*
* Revision 6.0  1997/08/25 18:34:00  madden
* Revision changed to 6.0
*
 * Revision 1.3  1996/11/27 18:04:28  madden
 * Added prototype for MakeBLAST0DbDesc
 *
 * Revision 1.2  1996/09/25  20:03:06  madden
 * Added prototype for GetParameterStack.
 *
 * Revision 1.1  1996/08/07  14:07:13  madden
 * Initial revision
 *
 * */
#ifndef __BLASTASN__
#define __BLASTASN__

#ifdef __cplusplus
extern "C" {
#endif

#include <blast.h>
#include <blastkar.h>
#define NLM_GENERATED_CODE_PROTO     /* Needed to get Scores_scaled_intsPtr*/
#include <blast18p.h>                   /* Patch for score-set */
#include <objblst2.h>


BLAST0ResultPtr LIBCALL MakeBLAST0Result PROTO((BlastSearchBlkPtr search, Boolean get_query_seq, Boolean get_db_seq));

BLAST0MatrixPtr LIBCALL GetBLAST0Matrix PROTO((BLAST_ScoreBlkPtr sbp, Boolean fullreport));

BLAST0KABlkPtr LIBCALL GetBLAST0KABlk PROTO((BLAST_ScoreBlkPtr sbp));

ValNodePtr LIBCALL GetParameterStack PROTO((BlastSearchBlkPtr search, ValNodePtr stp, Boolean old, Boolean stats));

/* Produce a BLAST0DbDesc from the ReadDBFILEPtr. */
BLAST0DbDescPtr LIBCALL MakeBLAST0DbDesc PROTO((ReadDBFILEPtr rdfp));

#ifdef __cplusplus
}
#endif
#endif /* !__BLASTASN__ */
