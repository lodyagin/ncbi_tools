/* $Id: phi_extend.h,v 1.1 2003/09/09 22:04:02 dondosha Exp $

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
* ===========================================================================
*****************************************************************************

File name: phi_extend.h

Author: Ilya Dondoshansky

Contents: Word finder for PHI-BLAST

Detailed Contents: 

******************************************************************************
 * $Revision: 1.1 $
 * */

#include <algo/blast/core/blast_extend.h>
#include <algo/blast/core/blast_util.h>

#ifndef PHI_EXTEND__H
#define PHI_EXTEND__H

#ifdef __cplusplus
extern "C" {
#endif

Int4 PHIBlastWordFinder(BLAST_SequenceBlk* subject, 
        BLAST_SequenceBlk* query, LookupTableWrap* lookup_wrap,
        Int4** matrix, BlastInitialWordParameters* word_params,
        BLAST_ExtendWord* ewp, Uint4* q_offsets, Uint4* s_offsets,
        Int4 max_hits, BlastInitHitList* init_hitlist);

#ifdef __cplusplus
}
#endif

#endif /* PHI_LOOKUP__H */