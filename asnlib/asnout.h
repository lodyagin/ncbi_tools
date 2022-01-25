/* asnout.h
* ===========================================================================
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
* File Name:  asnout.h
*
* Author:  James Ostell
*
* Version Creation Date: 1/1/91
*
* $Revision: 6.0 $
*
* File Description:
*   typedefs and prototypes used internally by asnout.c
*
* Modifications:
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
*
* $Log: asnout.h,v $
* Revision 6.0  1997/08/25 18:10:16  madden
* Revision changed to 6.0
*
* Revision 5.1  1996/12/03 21:43:48  vakatov
* Adopted for 32-bit MS-Windows DLLs
*
 * Revision 5.0  1996/05/28  14:00:29  ostell
 * Set to revision 5.0
 *
 * Revision 4.0  1995/07/26  13:47:38  ostell
 * force revision to 4.0
 *
 * Revision 2.2  1995/05/15  18:38:28  ostell
 * added Log line
 *
*
* ==========================================================================
*/

#ifndef _ASNOUT_
#define _ASNOUT_

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif


#define SYMBOLLEN_MAX	31	/* maximum length of defined symbol in chars */

typedef struct tmptype {
	AsnTypePtr ptr;
	Int2 index;
	struct tmptype PNTR next;
} AsnTmpType, PNTR AsnTmpTypePtr;

typedef struct tmpvalue {
	AsnValxNodePtr val;
	Int2 index;
	struct tmpvalue PNTR next;
} AsnTmpValue, PNTR AsnTmpValuePtr;

typedef struct asnproc {
	AsnTmpValue root_tmpvalue;
	AsnTmpValuePtr curr_tmpvalue;
	AsnTmpType root_tmptype;
	AsnTmpTypePtr curr_tmptype;
} AsnProc, PNTR AsnProcPtr;


/*****************************************************************************
*
*   prototypes
*
*****************************************************************************/

NLM_EXTERN void AsnOutput PROTO((CharPtr filename, AsnModulePtr amp, Boolean loader, Int2 maxDefineLength));
NLM_EXTERN Int2 AsnOutFindValue PROTO((AsnProcPtr app, AsnValxNodePtr avnp));
NLM_EXTERN AsnTypePtr AsnOutAddType PROTO((AsnProcPtr app, AsnTypePtr atp));
NLM_EXTERN Boolean AsnOutNewType PROTO((AsnProcPtr app, AsnTypePtr atp));
NLM_EXTERN AsnTmpTypePtr AsnOutFindType PROTO((AsnProcPtr app, AsnTypePtr atp));
NLM_EXTERN void AsnOutNewValueChain PROTO((AsnProcPtr app, AsnValxNodePtr avnp));
NLM_EXTERN void AsnOutDefineType PROTO((FILE *fp, AsnProcPtr app, AsnTypePtr atp, Int2 maxDefineLength));
NLM_EXTERN void AsnOutDefineElement PROTO((FILE *fp, AsnProcPtr app, AsnTypePtr atp, CharPtr buf, CharPtr pnt, Int2 maxDefineLength));


#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif
