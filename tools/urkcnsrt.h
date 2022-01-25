/*
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
* File Name: urkcnsrt.h
*
* Author(s): John Kuzio
*
* Version Creation Date: 98-01-01
*
* $Revision: 6.6 $
*
* File Description: consort header
*
* Modifications:
* --------------------------------------------------------------------------
* Date       Name        Description of modification
* --------------------------------------------------------------------------
* $Log: urkcnsrt.h,v $
* Revision 6.6  1998/09/28 16:36:11  kuzio
* no met orf check
*
* Revision 6.5  1998/09/16 17:46:44  kuzio
* cvs logging
*
*
* ==========================================================================
*/

#ifndef _CONSORT__
#define _CONSORT__

#include <ncbi.h>
#include <accentr.h>
#include <urktree.h>

#ifdef __cplusplus
extern "C" {
#endif

#define CONMATSIZE    512

/* codon usage tree from genome */

extern TreeNodePtr ConsortSeqEntry (SeqEntryPtr sep);
extern Int4Ptr     ConformSeqEntry (SeqEntryPtr sep);

/* exploration:  ORFs related by codon usage */

extern Int4Ptr NewCodonTable (void);
extern Int4Ptr FreeCodonTable (Int4Ptr cutp);
extern Int4Ptr CodonTableFromSeqLoc (BioseqPtr bsp, SeqLocPtr slp);
extern void AddSeqLocToCodonTable (Int4Ptr cutp, BioseqPtr bsp,
                                   SeqLocPtr slp, Boolean flagAdd);
extern Int4Ptr MergeCodonTables (Int4Ptr cutp1, Int4Ptr cutp2);

extern FloatHi Confide (Int4Ptr cutgene, Int4Ptr cutgbl);
extern void Conform (Int4Ptr freq, FILE *fn);

extern ValNodePtr ClearNonMetOrfs (ValNodePtr orflist);

#ifdef __cplusplus
}
#endif

#endif
