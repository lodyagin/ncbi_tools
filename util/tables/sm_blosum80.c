/*  $Id: sm_blosum80.c,v 1.2 2005/05/16 12:17:43 madden Exp $
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
* Author:  Aaron Ucko (via ./convert_scoremat.pl)
*
* File Description:
*   Protein alignment score matrices; shared between the two toolkits.
*
* ===========================================================================
*/

#include <util/tables/raw_scoremat.h>

/* Matrix made by matblas from blosum80.iij */
/* * column uses minimum score */
/* BLOSUM Clustered Scoring Matrix in 1/2 Bit Units */
/* Blocks Database = /data/blocks_5.0/blocks.dat */
/* Cluster Percentage: >= 80 */
/* Entropy =   0.9868, Expected =  -0.7442 */

static const TNCBIScore s_Blosum80PSM[24][24] = {
    /*       A,  R,  N,  D,  C,  Q,  E,  G,  H,  I,  L,  K,
             M,  F,  P,  S,  T,  W,  Y,  V,  B,  Z,  X,  * */
    /*A*/ {  5, -2, -2, -2, -1, -1, -1,  0, -2, -2, -2, -1,
            -1, -3, -1,  1,  0, -3, -2,  0, -2, -1, -1, -6 },
    /*R*/ { -2,  6, -1, -2, -4,  1, -1, -3,  0, -3, -3,  2,
            -2, -4, -2, -1, -1, -4, -3, -3, -2,  0, -1, -6 },
    /*N*/ { -2, -1,  6,  1, -3,  0, -1, -1,  0, -4, -4,  0,
            -3, -4, -3,  0,  0, -4, -3, -4,  4,  0, -1, -6 },
    /*D*/ { -2, -2,  1,  6, -4, -1,  1, -2, -2, -4, -5, -1,
            -4, -4, -2, -1, -1, -6, -4, -4,  4,  1, -1, -6 },
    /*C*/ { -1, -4, -3, -4,  9, -4, -5, -4, -4, -2, -2, -4,
            -2, -3, -4, -2, -1, -3, -3, -1, -4, -4, -1, -6 },
    /*Q*/ { -1,  1,  0, -1, -4,  6,  2, -2,  1, -3, -3,  1,
             0, -4, -2,  0, -1, -3, -2, -3,  0,  3, -1, -6 },
    /*E*/ { -1, -1, -1,  1, -5,  2,  6, -3,  0, -4, -4,  1,
            -2, -4, -2,  0, -1, -4, -3, -3,  1,  4, -1, -6 },
    /*G*/ {  0, -3, -1, -2, -4, -2, -3,  6, -3, -5, -4, -2,
            -4, -4, -3, -1, -2, -4, -4, -4, -1, -3, -1, -6 },
    /*H*/ { -2,  0,  0, -2, -4,  1,  0, -3,  8, -4, -3, -1,
            -2, -2, -3, -1, -2, -3,  2, -4, -1,  0, -1, -6 },
    /*I*/ { -2, -3, -4, -4, -2, -3, -4, -5, -4,  5,  1, -3,
             1, -1, -4, -3, -1, -3, -2,  3, -4, -4, -1, -6 },
    /*L*/ { -2, -3, -4, -5, -2, -3, -4, -4, -3,  1,  4, -3,
             2,  0, -3, -3, -2, -2, -2,  1, -4, -3, -1, -6 },
    /*K*/ { -1,  2,  0, -1, -4,  1,  1, -2, -1, -3, -3,  5,
            -2, -4, -1, -1, -1, -4, -3, -3, -1,  1, -1, -6 },
    /*M*/ { -1, -2, -3, -4, -2,  0, -2, -4, -2,  1,  2, -2,
             6,  0, -3, -2, -1, -2, -2,  1, -3, -2, -1, -6 },
    /*F*/ { -3, -4, -4, -4, -3, -4, -4, -4, -2, -1,  0, -4,
             0,  6, -4, -3, -2,  0,  3, -1, -4, -4, -1, -6 },
    /*P*/ { -1, -2, -3, -2, -4, -2, -2, -3, -3, -4, -3, -1,
            -3, -4,  8, -1, -2, -5, -4, -3, -2, -2, -1, -6 },
    /*S*/ {  1, -1,  0, -1, -2,  0,  0, -1, -1, -3, -3, -1,
            -2, -3, -1,  5,  1, -4, -2, -2,  0,  0, -1, -6 },
    /*T*/ {  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -2, -1,
            -1, -2, -2,  1,  5, -4, -2,  0, -1, -1, -1, -6 },
    /*W*/ { -3, -4, -4, -6, -3, -3, -4, -4, -3, -3, -2, -4,
            -2,  0, -5, -4, -4, 11,  2, -3, -5, -4, -1, -6 },
    /*Y*/ { -2, -3, -3, -4, -3, -2, -3, -4,  2, -2, -2, -3,
            -2,  3, -4, -2, -2,  2,  7, -2, -3, -3, -1, -6 },
    /*V*/ {  0, -3, -4, -4, -1, -3, -3, -4, -4,  3,  1, -3,
             1, -1, -3, -2,  0, -3, -2,  4, -4, -3, -1, -6 },
    /*B*/ { -2, -2,  4,  4, -4,  0,  1, -1, -1, -4, -4, -1,
            -3, -4, -2,  0, -1, -5, -3, -4,  4,  0, -1, -6 },
    /*Z*/ { -1,  0,  0,  1, -4,  3,  4, -3,  0, -4, -3,  1,
            -2, -4, -2,  0, -1, -4, -3, -3,  0,  4, -1, -6 },
    /*X*/ { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -6 },
    /***/ { -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6,
            -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6,  1 }
};
const SNCBIPackedScoreMatrix NCBISM_Blosum80 = {
    "ARNDCQEGHILKMFPSTWYVBZX*",
    s_Blosum80PSM[0],
    -6
};
