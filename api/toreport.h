/*  toreport.h
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
* File Name:  toreport.h
*
* Author:  Jonathan Kans
*   
* Version Creation Date: 9/24/91
*
* $Revision: 6.0 $
*
* File Description:  Generates a sequence report from a SeqView array
*
* Modifications:  
* --------------------------------------------------------------------------
* Date	   Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
*
* $Log: toreport.h,v $
* Revision 6.0  1997/08/25 18:08:05  madden
* Revision changed to 6.0
*
* Revision 5.1  1997/06/19 18:39:36  vakatov
* [WIN32,MSVC++]  Adopted for the "NCBIOBJ.LIB" DLL'ization
*
* Revision 5.0  1996/05/28 13:23:23  ostell
* Set to revision 5.0
*
 * Revision 4.0  1995/07/26  13:49:01  ostell
 * force revision to 4.0
 *
 * Revision 2.8  1995/05/15  21:46:05  ostell
 * added Log line
 *
*
*
* ==========================================================================
*/

#ifndef _TOREPORT_
#define _TOREPORT_

#include <sequtil.h>

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* ----- Entry Tag Definitions ----- */

#define BEGENTRY      1
#define PUBSET        2
#define NPSET         3
#define SEGSET        4
#define OTHERSET      5
#define ENDSET        6
#define RAWNUCL       7
#define RAWPROT       8
#define VIRTSEQ       9
#define CONSTSEQ      10
#define MEMBERS       11
#define SEGMENTS      12 
#define PDBSET        13
#define NCBI		  15
#define EMBL		  16
#define DDBJ		  17
#define PRF           18
#define PDB			  19
#define GENBANK       20
#define PIR           21
#define SWISSPROT     22
#define PATENT        23
#define XREF          24
#define ORGANISM      25
#define PUBINFO       26
#define REFERENCE     27
#define COMMENT       28
#define SOURCE        29
#define ORIGIN        30
#define RELDATE		  31
#define METHOD        32

#define FEAT_GENE    40 
#define FEAT_ORG     41
#define FEAT_CDRGN   42
#define FEAT_PROT    43 
#define FEAT_RNA     44
#define FEAT_PUB     45  
#define FEAT_SEQ     46
#define FEAT_IMP     47
#define FEAT_REGION  48
#define FEAT_COMMENT 49
#define FEAT_BOND    50
#define FEAT_SITE    51

#define FEAT_NUM    52
#define FEAT_PSEC   53
#define FEAT_NSR    54
#define FEAT_HET    55

#define FEAT_OTHER  128

#define SEQUENCE    200

#define OTHER       256

#define BADTYPE     512


/* ----- Flattened Array Type Definition ----- */

typedef struct seqView {
    Int2    nType;
    Int2    nHeight;
    Pointer pDat;     
    Pointer pAux;
} SeqView, PNTR SeqViewPtr;

/* ----- Formatting Callback Function Prototype ----- */

typedef Boolean (*ReportProc) PROTO((VoidPtr inst, Int2 _class, CharPtr string));

/* ----- Standard Save-To-File Callback Function Prototype ----- */

NLM_EXTERN Boolean StdReportProc PROTO((VoidPtr inst, Int2 _class, CharPtr string));

/* ----- Primitive Function Prototypes ----- */

NLM_EXTERN SeqViewPtr FlatSeqEntryNew PROTO((SeqEntryPtr pEntry, Int2Ptr nItems));
NLM_EXTERN Boolean FlatSeqToReport PROTO((SeqViewPtr pView, Int2 nItems, ReportProc proc, VoidPtr inst, Boolean showSeq, Int2 charsPerLine));
NLM_EXTERN SeqViewPtr FlatSeqEntryFree PROTO((SeqViewPtr pView));

/* ----- Compound Function Prototypes ----- */

NLM_EXTERN Boolean SeqEntryToReport PROTO((SeqEntryPtr pEntry, ReportProc proc, VoidPtr inst, Boolean showSeq, Int2 charsPerLine));
NLM_EXTERN Boolean SeqEntryToFile PROTO((SeqEntryPtr pEntry, FILE *fp, Boolean showSeq, Int2 charsPerLine, Boolean checkform));

#ifdef __cplusplus
}
#endif

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif

