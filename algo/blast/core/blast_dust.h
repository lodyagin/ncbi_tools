/* $Id: blast_dust.h,v 1.8 2003/09/09 14:21:14 coulouri Exp $
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's offical duties as a United States Government employee and
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

File name: blast_filter.h

Author: Ilya Dondoshansky

Contents: DUST filtering functions.

Detailed Contents: 

******************************************************************************
 * $Revision: 1.8 $
 * */
#ifndef __BLAST_DUST__
#define __BLAST_DUST__

#ifdef __cplusplus
extern "C" {
#endif

#include <algo/blast/core/blast_def.h>
#include <algo/blast/core/blast_stat.h>  /* Needed for NCBI4NA_TO_BLASTNA definition only */

Int2 SeqBufferDust (Uint1* sequence, Int4 length, Int4 offset,
                    Int2 level, Int2 window, Int2 minwin, Int2 linker,
                    BlastSeqLoc** dust_loc);

#ifdef __cplusplus
}
#endif
#endif /* !__BLAST_FILTER__ */