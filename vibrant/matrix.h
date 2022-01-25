/*  $Id: matrix.h,v 6.0 1997/08/25 18:55:56 madden Exp $
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*            National Center for Biotechnology Information (NCBI)
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government do not place any restriction on its use or reproduction.
*  We would, however, appreciate having the NCBI and the author cited in
*  any work or product based on this material
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
* ===========================================================================
*
* Author:  Denis Vakatov
*
* File Description: 
*   Basic matrix & vector handling library.
*   NOTE:  most of advanced functions are valid for square matrices only!
*
* ===========================================================================
* $Log: matrix.h,v $
* Revision 6.0  1997/08/25 18:55:56  madden
* Revision changed to 6.0
*
* Revision 5.1  1997/05/20 14:40:02  vakatov
* Initial revision
*
*
* ==========================================================================
*/

#ifndef MATRIX_H
#define MATRIX_H

#include <ncbilcl.h>
#include <ncbistd.h>

/* The following preprocessor variables can be set on the compiling stage:
 * #TEST_MODULE_MATRIX -- to build test application (needs CORELIB library)
 * #CHECK_RANGE_MATRIX -- to check for the matrix indexes be in-range (slow) 
 */

#ifdef __cplusplus
extern "C" {
#endif

struct Nlm_MatrixTag;
typedef struct Nlm_MatrixTag *Nlm_Matrix;

extern Nlm_Matrix Nlm_MatrixNew       (Nlm_Uint4 n_row, Nlm_Uint4 n_column);
extern void       Nlm_MatrixDelete    (Nlm_Matrix mat);

extern void       Nlm_MatrixSetNode   (Nlm_Matrix mat,
                                       Nlm_Uint4 row, Nlm_Uint4 column,
                                       Nlm_FloatHi value);
extern void       Nlm_MatrixSetRow    (Nlm_Matrix       mat,
                                       Nlm_Uint4        row,
                                       const Nlm_Matrix mat_src,
                                       Nlm_Uint4        row_src);
extern void       Nlm_MatrixSetColumn (Nlm_Matrix       mat,
                                       Nlm_Uint4        column,
                                       const Nlm_Matrix mat_src,
                                       Nlm_Uint4        column_src);


extern Nlm_FloatHi Nlm_MatrixNode     (const Nlm_Matrix mat,
                                       Nlm_Uint4 row, Nlm_Uint4 column);
extern Nlm_Matrix  Nlm_MatrixRow      (const Nlm_Matrix mat, Nlm_Uint4 row);
extern Nlm_Matrix  Nlm_MatrixColumn   (const Nlm_Matrix mat, Nlm_Uint4 column);

extern Nlm_Boolean Nlm_MatrixCompare  (const Nlm_Matrix mat1,
                                       const Nlm_Matrix mat2);
extern Nlm_Matrix  Nlm_MatrixCopy     (const Nlm_Matrix mat);
extern Nlm_Matrix  Nlm_MatrixTranspose(const Nlm_Matrix mat);
extern Nlm_Matrix  Nlm_MatrixMultiply (const Nlm_Matrix mat_left,
                                       const Nlm_Matrix mat_right);
extern Nlm_Matrix  Nlm_MatrixSolve    (const Nlm_Matrix mat,
                                       const Nlm_Matrix vector);
extern Nlm_Matrix  Nlm_MatrixInvert   (const Nlm_Matrix mat);
extern void        Nlm_MatrixPrint    (const Nlm_Matrix mat,
                                       FILE *fd, const Char *descr);

#ifdef __cplusplus
}
#endif

#endif  /* MATRIX_H */
