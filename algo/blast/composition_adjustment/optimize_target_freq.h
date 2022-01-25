/* $Id: optimize_target_freq.h,v 1.6 2005/12/01 13:54:04 gertz Exp $
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
 * ===========================================================================*/

/**
 * @file optimize_target_freq.h
 * @author E. Michael Gertz
 *
 * Exports for optimized_target_freq.c
 */

#ifndef __OPTIMIZE_TARGET_FREQ__
#define __OPTIMIZE_TARGET_FREQ__

#include <algo/blast/core/blast_export.h>

#ifdef __cplusplus
extern "C" {
#endif

NCBI_XBLAST_EXPORT
int
Blast_OptimizeTargetFrequencies(double x[],
                                int alphsize,
                                int * iterations,
                                const double q[],
                                const double row_sums[],
                                const double col_sums[],
                                int constrain_rel_entropy,
                                double relative_entropy,
                                double tol,
                                int maxits);

#ifdef __cplusplus
}
#endif

#endif
