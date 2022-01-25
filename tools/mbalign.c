/* $Id: mbalign.c,v 6.14 2000/05/24 19:51:28 dondosha Exp $
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
* File Name:  $RCSfile: mbalign.c,v $
*
* Author:  Webb Miller and Co.
* Adopted for NCBI standard libraries by Sergey Shavirin
*
* Initial Creation Date: 10/27/1999
*
* $Revision: 6.14 $
*
* File Description:
*        Alignment functions for Mega Blast program
*
* $Log: mbalign.c,v $
* Revision 6.14  2000/05/24 19:51:28  dondosha
* Do not use readdb when megablast is used for two sequences
*
* Revision 6.13  2000/05/17 17:12:39  dondosha
* Removed unused variable
*
* Revision 6.12  2000/05/12 19:45:32  dondosha
* Subject sequence is now passed in ncbi4na encoding
*
* Revision 6.11  2000/05/01 21:23:40  dondosha
* Changed greedy_align_basic to MegaBlastGreedyAlign, added MegaBlastAffineGreedyAlign
*
* Revision 6.10  2000/04/19 18:53:31  dondosha
* Bug fix: allow to extend through gap and masked characters, but not end-points
*
* Revision 6.9  2000/04/07 18:52:32  dondosha
* Do not allocate new space for flast_d[0,1] every time
*
* Revision 6.8  2000/04/05 16:30:07  dondosha
* Fixed a bug that allowed crossing the end of a strand or query
*
* Revision 6.7  2000/03/30 20:38:33  dondosha
* Returned one modification from 6.5 that was lost in 6.6
*
* Revision 6.6  2000/03/29 22:03:14  dondosha
* Made edit_script_new and edit_script_append public for use in MegaBlastNtWordExtend
*
* Revision 6.5  2000/02/01 22:42:11  dondosha
* Moved some definitions to mbalign.h; added 3 new routines
*
* Revision 6.4  2000/01/26 16:39:18  beloslyu
* change space_t to space_type because HP-UX compiler got confused with ANSI C space_t
*
* Revision 6.3  1999/11/24 20:38:09  shavirin
* Added possibility to produce Traditional Blast Output.
*
* Revision 6.2  1999/11/03 19:53:57  shavirin
* Fixed problem with Realloc() function.
*
* Revision 6.1  1999/10/29 15:37:12  shavirin
* Initial revision.
*
*
* Initial revision.
*
* ==========================================================================
*/

#include <mbutils.h>
#include <mbalign.h>
#include <blastdef.h>

/* -------- From original file edit.c ------------- */

static Uint4 edit_val_get(edit_op_t op)
{
    return op >> 2;
}
static Uint4 edit_opc_get(edit_op_t op)
{
    return op & EDIT_OP_MASK;
}

static edit_op_t edit_op_cons(Uint4 op, Uint4 val)
{
    return (val << 2) | (op & EDIT_OP_MASK);
}

static edit_op_t edit_op_inc(edit_op_t op, Uint4 n)
{
    return edit_op_cons(edit_opc_get(op), edit_val_get(op) + n);
}

static edit_op_t edit_op_inc_last(edit_script_t *es, Uint4 n)
{
    edit_op_t *last;
    ASSERT (es->num > 0);
    last = &(es->op[es->num-1]);
    *last = edit_op_inc(*last, n);
    return *last;
}

static Int4 edit_script_ready(edit_script_t *es, Uint4 n)
{
    edit_op_t *p;
    Uint4 m = n + n/2;
    
    if (es->size <= n) {
        p = Realloc(es->op, m*sizeof(edit_op_t));
        if (p == 0) {
            return 0;
        } else {
            es->op = p;
            es->size = m;
        }
    }
    return 1;
}

static Int4 edit_script_readyplus(edit_script_t *es, Uint4 n)
{
    if (es->size - es->num <= n)
        return edit_script_ready(es, n + es->num);
    return 1;
}

static Int4 edit_script_put(edit_script_t *es, Uint4 op, Uint4 n)
{
    if (!edit_script_readyplus(es, 2))
        return 0;
    es->last = op;
    ASSERT(op != 0);
    es->op[es->num] = edit_op_cons(op, n);
    es->num += 1;
    es->op[es->num] = 0; /* sentinal */

    return 1;
}

static edit_script_t *edit_script_init(edit_script_t *es)
{
	es->op = 0;
	es->size = es->num = 0;
	es->last = 0;
	edit_script_ready(es, 8);
	return es;
}
static edit_op_t *edit_script_first(edit_script_t *es)
{
    return es->num > 0 ? &es->op[0] : 0;
}

static edit_op_t *edit_script_next(edit_script_t *es, edit_op_t *op)
{
    /* XXX - assumes flat address space */
    if (&es->op[0] <= op && op < &es->op[es->num-1])
        return op+1;
    else
        return 0;
}
static Int4 edit_script_more(edit_script_t *data, Uint4 op, Uint4 k)
{
    if (op == EDIT_OP_ERR) {
        ErrPostEx(SEV_FATAL, 0, 0, 
                  "edit_script_more: bad opcode %d:%d", op, k);
        return -1;
    }
    
    if (edit_opc_get(data->last) == op)
        edit_op_inc_last(data, k);
    else
        edit_script_put(data, op, k);

    return 0;
}
/* External */
edit_script_t *edit_script_append(edit_script_t *es, edit_script_t *et)
{
    edit_op_t *op;
    
    for (op = edit_script_first(et); op; op = edit_script_next(et, op))
        edit_script_more(es, edit_opc_get(*op), edit_val_get(*op));

    return es;
}

/* External */
edit_script_t *edit_script_new(void)
{
    edit_script_t *es = MemNew(sizeof(*es));
    if (!es)
        return 0;

    return edit_script_init(es);
}

/* External */
edit_script_t *edit_script_free(edit_script_t *es)
{
    if (es) {
        if (es->op)
            MemFree(es->op);
        MemSet(es, 0, sizeof(*es));
        MemFree(es);
    }
    return 0;
}

/* External */
Int4 es_rep_len(edit_script_t *S, Int4Ptr n, UcharPtr p, 
                UcharPtr q, Int4Ptr match)
{
    edit_op_t op;
    Int4 len;

    len = *match = 0;
    while ((*n < S->num) && (op=edit_opc_get(S->op[*n])) == EDIT_OP_REP) {
        Int4 num = edit_val_get(S->op[*n]);
        len += num;
        while (num-- > 0)
            /* masked regions are lower case */
            if (toupper(*p++) == toupper(*q++)) ++(*match);
        *n += 1;
    }
    return len;
}

/* External. Used in mbutils.c */
Int4 es_indel_len(edit_script_t *S, Int4Ptr n, Int4Ptr i, Int4Ptr j)
{
    if (S->num <= *n)
        return 0;
    else {
        edit_op_t op = S->op[*n];
        Int4 len = edit_val_get(op);
        
        switch (edit_opc_get(op)) {
        case EDIT_OP_INS:
            *j += len;
            break;
        case EDIT_OP_DEL:
            *i += len;
            break;
        default:
            ErrPostEx(SEV_ERROR, 0, 0, "es_indel_len: cannot happen!");
        }
        *n += 1;
        return len;
    }
}

Int4 edit_script_del(edit_script_t *data, Uint4 k)
{
    return edit_script_more(data, EDIT_OP_DEL, k);
}

Int4 edit_script_ins(edit_script_t *data, Uint4 k)
{
    return edit_script_more(data, EDIT_OP_INS, k);
}
Int4 edit_script_rep(edit_script_t *data, Uint4 k)
{
    return edit_script_more(data, EDIT_OP_REP, k);
}

edit_script_t *edit_script_reverse_inplace(edit_script_t *es)
{
    Uint4 i;
    const Uint4 num = es->num;
    const Uint4 mid = num/2;
    const Uint4 end = num-1;
    
    for (i = 0; i < mid; ++i) {
        const edit_op_t t = es->op[i];
        es->op[i] = es->op[end-i];
        es->op[end-i] = t;
    }
    return es;
}

/* -------- From original file galign.c --------- */

space_ptr new_space(Int4 MAX_D)
{
    space_ptr p;
    Int4 amount;
    
    p = Nlm_Malloc(sizeof(space_type));
    amount = MAX_SPACE;
    p->space_array = Nlm_Malloc(sizeof(ThreeVal)*amount);
    p->used = 0; 
    p->size = amount;

    return p;
}

void free_space(space_ptr sp)
{
    sp->space_array = MemFree(sp->space_array);
    sp = MemFree(sp);
}

ThreeValPtr get_space(space_ptr S, Int4 amount)
{
    ThreeValPtr s = S->space_array+S->used;
    if (amount < 0) 
        return NULL;  

    if (S->used+amount > S->size) {
        ErrPostEx(SEV_FATAL, 0, 0, "out of space");
        return NULL;
    }

    S->used += amount;
    
    return s;
}

/* ----- */

static Int4 gcd(Int4 a, Int4 b)
{
    Int4 c;
    if (a < b) {
        c = a;
        a = b; b = c;
    }
    while ((c = a % b) != 0) {
        a = b; b = c;
    }

    return b;
}


static Int4 gdb3(Int4Ptr a, Int4Ptr b, Int4Ptr c)
{
    Int4 g;
    if (*b == 0) g = gcd(*a, *c);
    else g = gcd(*a, gcd(*b, *c));
    if (g > 1) {
        *a /= g;
        *b /= g;
        *c /= g;
    }
    return g;
}

static Int4 get_lastC(ThreeValPtr PNTR flast_d, Int4Ptr lower, Int4Ptr upper, 
                      Int4Ptr d, Int4 diag, Int4 Mis_cost, Int4Ptr row1)
{
    Int4 row;
    
    if (diag >= lower[(*d)-Mis_cost] && diag <= upper[(*d)-Mis_cost]) {
        row = flast_d[(*d)-Mis_cost][diag].C;
        if (row >= MAX(flast_d[*d][diag].I, flast_d[*d][diag].D)) {
            *d = *d-Mis_cost;
            *row1 = row;
            return sC;
        }
    }
    if (flast_d[*d][diag].I > flast_d[*d][diag].D) {
        *row1 = flast_d[*d][diag].I;
        return sI;
    } else {
        *row1 = flast_d[*d][diag].D;
        return sD;
    }
}

static Int4 get_last_ID(ThreeValPtr PNTR flast_d, Int4Ptr lower, Int4Ptr upper, 
                        Int4Ptr d, Int4 diag, Int4 GO_cost, 
                        Int4 GE_cost, Int4 IorD)
{
    Int4 ndiag; 
    Int4 row;

    if (IorD == sI) { ndiag = diag -1;} else ndiag = diag+1;
    if (ndiag >= lower[(*d)-GE_cost] && ndiag <=upper[(*d)-GE_cost]) 
        row = (IorD == sI)? flast_d[(*d)-GE_cost][ndiag].I: flast_d[(*d)-GE_cost][ndiag].D;
    else row = -100;
    if (ndiag >= lower[(*d)-GO_cost-GE_cost] && ndiag <= upper[(*d)-GO_cost-GE_cost] && row < flast_d[(*d)-GO_cost-GE_cost][ndiag].C) {
        *d = (*d)-GO_cost-GE_cost;
        return sC;
    }
    *d = (*d)-GE_cost;
    return IorD;
}

static Int4 get_lastI(ThreeValPtr PNTR flast_d, Int4Ptr lower, Int4Ptr upper, 
                      Int4Ptr d, Int4 diag, Int4 GO_cost, Int4 GE_cost)
{
    return get_last_ID(flast_d, lower, upper, d, diag, GO_cost, GE_cost, sI);
}


static int get_lastD(ThreeValPtr PNTR flast_d, Int4Ptr lower, Int4Ptr upper, 
                     Int4Ptr d, Int4 diag, Int4 GO_cost, Int4 GE_cost)
{
    return get_last_ID(flast_d, lower, upper, d, diag, GO_cost, GE_cost, sD);
}

/* --- From file align.c --- */
/* ----- */

Int4 get_last(Int4 **flast_d, Int4 d, Int4 diag, Int4 *row1)
{
    if (flast_d[d-1][diag-1] > MAX(flast_d[d-1][diag], flast_d[d-1][diag+1])) {
        *row1 = flast_d[d-1][diag-1];
        return diag-1;
    } 
    if (flast_d[d-1][diag] > flast_d[d-1][diag+1]) {
        *row1 = flast_d[d-1][diag];
        return diag;
    }
    *row1 = flast_d[d-1][diag+1];
    return diag+1;
}

/* Basic O(ND) time, O(N) space, alignment function. 
 * Parameters:
 *   s1, len1        - first sequence and its length
 *   s2, len2        - second sequence and its length
 *   reverse         - direction of alignment
 *   xdrop_threshold -
 *   mismatch_cost   -
 *   e1, e2          - endpoint of the computed alignment
 *   edit_script     -
 */
Int4 align_basic(
                 const UcharPtr s1, Int4 len1,
                 const UcharPtr s2, Int4 len2,
                 Int4 reverse, Int4 xdrop_threshold, Int4 match_cost, 
                 Int4 mismatch_cost,
                 Int4 *e1, Int4 *e2, edit_script_t *S)
{
    Int4 col,			/* column number */
        d,				/* current distance */
        k,				/* current diagonal */
        flower, fupper,            /* boundaries for searching diagonals */
        row,		        /* row number */
        MAX_D, 			/* maximum cost */
        ORIGIN,
        return_val = 0;
    Int4 **flast_d;		/* rows containing the last d */
    Int4 *max_row;		/* reached for cost d=0, ... len1.  */
    
    Int4 X_pen = xdrop_threshold;
    Int4 M_half = match_cost/2;
    Int4 Op_cost = mismatch_cost + M_half*2;
    Int4 D_diff = ICEIL(X_pen+M_half, Op_cost);
    
    Int4 *max_row_free, x, cur_max, b_diag = 0, best_diag = INT_MAX/2;
    Char nlower = 0, nupper = 0;
    space_ptr space;
    
    MAX_D = (Int4) (len1/ERROR_FRACTION + 1);
    ORIGIN = MAX_D + 2;
    *e1 = *e2 = 0;
    
    if (reverse)
	for (row = 0; row < len2 && row < len1 && (s1[len1-1-row] == s2[len2-1-row]); row++)
	    /*empty*/ ;
    else
	for (row = 0; row < len2 && row < len1 && (s1[row] == s2[row]); row++)
	    /*empty*/ ;
    
    *e1 = row;
    *e2 = row;
    if (row == len1) {
        if (S != NULL)
            edit_script_rep(S, row);
	/* hit last row; stop search */
	return 0;
    }
    flast_d = MemNew((MAX_D + 2) * sizeof(Int4 *));
    if (S==NULL) {
        space = 0;
        flast_d[0] = MemNew((MAX_D + MAX_D + 6) * sizeof(Int4) * 2);
        flast_d[1] = flast_d[0] + MAX_D + MAX_D + 6;
    } else {
        space = new_space(MAX_D);
        flast_d[0] = (Int4 *) get_space(space, 7)-ORIGIN+3;
        flast_d[1] = (Int4 *) get_space(space, 9)-ORIGIN+4;
    }
    max_row_free = MemNew(sizeof(Int4) * (MAX_D + 1 + D_diff));
    max_row = max_row_free + X_pen/Op_cost;
    for (k = 0; k < D_diff; k++)
	max_row_free[k] = 0;
    
    flast_d[0][ORIGIN] = row;
    max_row[0] = (row + row)*M_half;
    
    flower = ORIGIN - 1;
    fupper = ORIGIN + 1;

    d = 1;
    while (d <= MAX_D) {
	Int4 fl0, fu0;
	flast_d[d - 1][flower - 1] = flast_d[d - 1][flower] = -1;
	flast_d[d - 1][fupper] = flast_d[d - 1][fupper + 1] = -1;
	x = max_row[d - D_diff] + Op_cost * d - X_pen;
	x = ICEIL(x, M_half);	
	cur_max = 0;
	fl0 = flower;
	fu0 = fupper;
	for (k = fl0; k <= fu0; k++) {
	    /* get space for the next edit instruction */
	    row = MAX(flast_d[d - 1][k + 1], flast_d[d - 1][k]) + 1;
	    row = MAX(row, flast_d[d - 1][k - 1]);
	    col = row + k - ORIGIN;
	    if (row + col >= x)
		fupper = k;
	    else {
		if (k == flower)
		    flower++;
		else
		    flast_d[d][k] = -1;
		continue;
	    }
	    /* slide down the diagonal */
	    if (reverse) {
		while (row < len1 && col < len2 && s1[len1-1-row] == s2[len2-1-col]) {
		    ++row;
		    ++col;
		}
	    } else {
		while (row < len1 && col < len2 && s1[row] == s2[col]) {
		    ++row;
		    ++col;
		}
	    }
	    flast_d[d][k] = row;
	    if (row + col > cur_max) {
		cur_max = row + col;
		b_diag = k;
	    }
	    if (row == len1) {
		flower = k+1; nlower = 1;
	    }
	    if (col == len2) {
		fupper = k-1; nupper = 1;
	    }
	}
	k = cur_max*M_half - d * Op_cost;
	if (max_row[d - 1] < k) {
	    max_row[d] = k;
	    return_val = d;
	    best_diag = b_diag;
	    *e1 = flast_d[d][b_diag];
	    *e2 = (*e1)+b_diag-ORIGIN;
	} else {
	    max_row[d] = max_row[d - 1];
	}
	if (flower > fupper)
	    break;
	d++;
	if (!nlower) flower--; 
	if (!nupper) fupper++;
	if (S==NULL)
            flast_d[d] = flast_d[d - 2];
	else 
            flast_d[d] = (Int4 *) get_space(space, fupper-flower+5)-flower+2;
    }
    
    /*printf("return\n"); */
    if (S!=NULL) { /*trace back*/
        Int4 row1, col1, diag1, diag;
        d = return_val; diag = best_diag;
        row = *e1; col = *e2;
        while (d > 0) {
            diag1 = get_last(flast_d, d, diag, &row1);
            col1 = row1+diag1-ORIGIN;
            if (diag1 == diag) {
                if (row-row1 > 0) edit_script_rep(S, row-row1);
            } else if (diag1 < diag) {
                if (row-row1 > 0) edit_script_rep(S, row-row1);
                edit_script_ins(S,1);
            } else {
                if (row-row1-1> 0) edit_script_rep(S, row-row1-1);
                edit_script_del(S, 1);
            }
            d--; diag = diag1; col = col1; row = row1;
        }
        edit_script_rep(S, flast_d[0][ORIGIN]);
        if (!reverse) 
            edit_script_reverse_inplace(S);
        free_space(space);
    } else {
        free(flast_d[0]);
    }
    free(max_row_free);
    free(flast_d);
    return return_val;
}

/* Basic O(ND) time, O(N) space, alignment function. 
 * Parameters:
 *   s1, len1        - first sequence and its length
 *   s2, len2        - second sequence and its length
 *   reverse         - direction of alignment
 *   xdrop_threshold -
 *   mismatch_cost   -
 *   e1, e2          - endpoint of the computed alignment
 *   edit_script     -
 */

int gap_align_basic(
	const UcharPtr s1, Int4 len1,
	const UcharPtr s2, Int4 len2,
	Int4 reverse, Int4 xdrop_threshold, 
	Int4 match_score,
	Int4 mismatch_score,
	Int4 gap_open,
	Int4 gap_extend,
	Int4 *e1, Int4 *e2, edit_script_t *S)
{
    Int4 col,			/* column number */
        d,			/* current distance */
        k,			/* current diagonal */
        flower, fupper,         /* boundaries for searching diagonals */
        row,		        /* row number */
        MAX_D, 			/* maximum cost */
        ORIGIN,
        return_val = 0;
    ThreeValPtr PNTR flast_d;	/* rows containing the last d */
    Int4Ptr max_row;		/* reached for cost d=0, ... len1.  */
    Int4 Mis_cost, GO_cost, GE_cost;
    Int4 D_diff, gd;
    Int4 M_half;
    Int4 max_cost;
    Int4 *lower, *upper;
    
    Int4 *max_row_free, x, cur_max, b_diag = 0, best_diag = INT_MAX/2;
    Char nlower = 0, nupper = 0;
    space_ptr space;
    Int4 stop_condition;
    Int4 max_d;
    Int4Ptr uplow_free;
    
    if (match_score % 2 == 1) {
        match_score *= 2;
        mismatch_score *= 2;
        xdrop_threshold *= 2;
        gap_open *= 2;
        gap_extend *= 2;
    }
    
    M_half = match_score/2;
    if (gap_open == 0 && gap_extend == mismatch_score+M_half) {
        return align_basic(s1, len1, s2, len2, reverse, xdrop_threshold,
                           match_score, mismatch_score, e1, e2, S);
    }

    Mis_cost = mismatch_score + match_score;
    GO_cost = gap_open;
    GE_cost = gap_extend+M_half;
    gd = gdb3(&Mis_cost, &GO_cost, &GE_cost);
    D_diff = ICEIL(xdrop_threshold+M_half, gd);
    
    
    MAX_D = (Int4) (len1/ERROR_FRACTION + 1);
    max_d = MAX_D*GE_cost;
    ORIGIN = MAX_D + 2;
    max_cost = MAX(Mis_cost, GO_cost+GE_cost);
    *e1 = *e2 = 0;
    
    if (reverse)
	for (row = 0; row < len2 && row < len1 && (s1[len1-1-row] == s2[len2-1-row]); row++)
	    /*empty*/ ;
    else
	for (row = 0; row < len2 && row < len1 && (s1[row] == s2[row]); row++)
	    /*empty*/ ;
    
    *e1 = row;
    *e2 = row;
    if (row == len1 || row == len2) {
        if (S != NULL)
            edit_script_rep(S, row);
	/* hit last row; stop search */
	return row*match_score;
    }
    flast_d = MemNew((MAX(max_d, max_cost) + 2) * sizeof(ThreeValPtr));
    if (S==NULL) {
        space = 0;
        flast_d[0] = MemNew((MAX_D + MAX_D + 6) * sizeof(ThreeVal) * (max_cost+1));
        for (k = 1; k <= max_cost; k++)
            flast_d[k] = flast_d[k-1] + MAX_D + MAX_D + 6;
    } else {
        space = new_space(MAX_D);
        flast_d[0] = get_space(space, 1)-ORIGIN;
        flast_d[1] = get_space(space, 3)-ORIGIN+1;
    }
    max_row_free = MemNew(sizeof(Int4) * (max_d + 1 + D_diff));
    max_row = max_row_free + D_diff;
    for (k = 0; k < D_diff; k++)
	max_row_free[k] = 0;
    uplow_free = MemNew(sizeof(Int4)*2*(max_d+1+max_cost));
    lower = uplow_free;
    upper = uplow_free+max_d+1+max_cost;
    /* next 3 lines set boundary for -1,-2,...,-max_cost+1*/
    for (k = 0; k < max_cost; k++) {lower[k] =LARGE;  upper[k] = -LARGE;}
    lower += max_cost;
    upper += max_cost; 
    
    flast_d[0][ORIGIN].C = row;
    flast_d[0][ORIGIN].I = flast_d[0][ORIGIN].D = -2;
    max_row[0] = (row + row)*M_half;
    lower[0] = upper[0] = ORIGIN;
    
    flower = ORIGIN - 1;
    fupper = ORIGIN + 1;
    
    d = 1;
    stop_condition = 1;
    while (d <= max_d) {
	Int4 fl0, fu0;
	x = max_row[d - D_diff] + gd * d - xdrop_threshold;
	x = ICEIL(x, M_half);
	if (x < 0) x=0;
	cur_max = 0;
	fl0 = flower;
	fu0 = fupper;
	for (k = fl0; k <= fu0; k++) {
	    row = -2;
	    if (k+1 <= upper[d-GO_cost-GE_cost] && k+1 >= lower[d-GO_cost-GE_cost]) 
                row = flast_d[d-GO_cost-GE_cost][k+1].C;
	    if (k+1  <= upper[d-GE_cost] && k+1 >= lower[d-GE_cost] &&
		row < flast_d[d-GE_cost][k+1].D) 
                row = flast_d[d-GE_cost][k+1].D;
	    row++;
	    if (row+row+k-ORIGIN >= x) 
	      flast_d[d][k].D = row;
	    else flast_d[d][k].D = -2;
	    row = -1; 
	    if (k-1 <= upper[d-GO_cost-GE_cost] && k-1 >= lower[d-GO_cost-GE_cost]) 
                row = flast_d[d-GO_cost-GE_cost][k-1].C;
	    if (k-1  <= upper[d-GE_cost] && k-1 >= lower[d-GE_cost] &&
		row < flast_d[d-GE_cost][k-1].I) 
                row = flast_d[d-GE_cost][k-1].I;
	    if (row+row+k-ORIGIN >= x) 
                flast_d[d][k].I = row;
	    else flast_d[d][k].I = -2;
            
	    row = MAX(flast_d[d][k].I, flast_d[d][k].D);
	    if (k <= upper[d-Mis_cost] && k >= lower[d-Mis_cost]) 
                row = MAX(flast_d[d-Mis_cost][k].C+1,row);
            
	    col = row + k - ORIGIN;
	    if (row + col >= x)
		fupper = k;
	    else {
		if (k == flower)
		    flower++;
		else
		    flast_d[d][k].C = -2;
		continue;
	    }
	    /* slide down the diagonal */
	    if (reverse) {
		while (row < len1 && col < len2 && s1[len1-1-row] == s2[len2-1-col]) {
		    ++row;
		    ++col;
		}
	    } else {
		while (row < len1 && col < len2 && s1[row] == s2[col]) {
		    ++row;
		    ++col;
		}
	    }
	    flast_d[d][k].C = row;
	    if (row + col > cur_max) {
		cur_max = row + col;
		b_diag = k;
	    }
	    if (row == len1) {
		flower = k; nlower = k+1;
	    }
	    if (col == len2) {
		fupper = k; nupper = k-1;
	    }
	}
	k = cur_max*M_half - d * gd;
	if (max_row[d - 1] < k) {
	    max_row[d] = k;
	    return_val = d;
	    best_diag = b_diag;
	    *e1 = flast_d[d][b_diag].C;
	    *e2 = (*e1)+b_diag-ORIGIN;
	} else {
	    max_row[d] = max_row[d - 1];
	}
	if (flower <= fupper) {
            stop_condition++;
            lower[d] = flower; upper[d] = fupper;
	} else {
            lower[d] = LARGE; upper[d] = -LARGE;
	}
	if (lower[d-max_cost] <= upper[d-max_cost]) stop_condition--;
	if (stop_condition == 0) break;
	d++;
	flower = MIN(lower[d-Mis_cost], MIN(lower[d-GO_cost-GE_cost], lower[d-GE_cost])-1);
	if (nlower) flower = MAX(flower, nlower);
	fupper = MAX(upper[d-Mis_cost], MAX(upper[d-GO_cost-GE_cost], upper[d-GE_cost])+1);
	if (nupper) fupper = MIN(fupper, nupper);
	if (S==NULL) {
            if (d > max_cost)
                flast_d[d] = flast_d[d - max_cost-1];
	} else 
            flast_d[d] = get_space(space, fupper-flower+1)-flower;
    }
    
    /*printf("return\n"); */
    if (S!=NULL) { /*trace back*/
        Int4 row1, diag, state;
        d = return_val; diag = best_diag;
        row = *e1; state = sC;
        while (d > 0) {
            if (state == sC) {
                /* diag will not be changed*/
                state = get_lastC(flast_d, lower, upper, &d, diag, Mis_cost, &row1);
                if (row-row1 > 0) edit_script_rep(S, row-row1);
                row = row1;
            } else {
                if (state == sI) {
                    /*row unchanged */
                    state = get_lastI(flast_d, lower, upper, &d, diag, GO_cost, GE_cost);
                    diag--;
                    edit_script_ins(S,1);
                } else {
                    edit_script_del(S,1);
                    state = get_lastD(flast_d, lower, upper, &d, diag, GO_cost, GE_cost);
                    diag++;
                    row--;
                }
            }
        }
        edit_script_rep(S, flast_d[0][ORIGIN].C);
        if (!reverse) 
            edit_script_reverse_inplace(S);
        free_space(space);
    } else {
        free(flast_d[0]);
    }
    return_val = max_row[return_val];
    free(max_row_free);
    free(flast_d);
    free(uplow_free);
    return return_val;
}

/* -------- From original file mbalign.c ------------- */

static Int4 extend_gap_align_internal(gal_t *m, Int4 X, Int4 M, Int4 R, 
                                      Int4 G, Int4 E, Int4 l1, 
                                      const UcharPtr s1, Int4 l2, 
                                      const UcharPtr s2, Int4 scr)
{
    Int4 b1, b2, e1, e2, i;
    edit_script_t *scr1 = 0, *scr2 = 0;
    
    if (scr) {
        scr1 = edit_script_new();
        scr2 = edit_script_new();
    }
    
    /* nb. msp is 1-indexed; align is 1-indexed */
    
    /* Align forward from end of the msp, and backwards from the beginning. */
    i = gap_align_basic(s1+m->e1, l1-m->e1, s2+m->e2, l2-m->e2, 0,X,M,R,G,E, &e1,&e2,scr1)
        + gap_align_basic(s1,       m->b1-1,  s2,       m->b2-1,  1,X,M,R,G,E, &b1,&b2,scr2);
    
    if (scr) {
	/* The improved alignment is the computed scripts, joined by the msp. */
        edit_script_rep(scr2, m->e1 - m->b1 + 1);
        edit_script_append(scr2, scr1);
        edit_script_free(scr1);
    }
    
    if (G!=0 || E*2!=R*2+M) i += (m->e1-m->b1+1)*M;
    m->e1 += e1;
    m->e2 += e2;
    m->b1 -= b1;
    m->b2 -= b2;
    m->S = scr2;
    return i;
}

Int4 mb_extend_gap_align_script(gal_t *mq, Int4 X, Int4 M, Int4 R, Int4 G, 
                                Int4 E, const UcharPtr s1, Int4 l1, 
                                const UcharPtr s2, Int4 l2)
{
    Int4 i = extend_gap_align_internal(mq, X, M, R, G, E, l1, s1, l2, s2, 1);
    /*print_align_lav(mq->diag, s1, s2, 
      mq->b1, mq->e1, mq->b2, mq->e2, mq->S);*/
    return i;
}

Int4 mb_extend_gap_align(gal_t *mq, Int4 X, Int4 M, Int4 R, Int4 G, 
                         Int4 E, const UcharPtr s1, Int4 l1, 
                         const UcharPtr s2, Int4 l2)
{
    return extend_gap_align_internal(mq, X, M, R, G, E, l1, s1, l2, s2, 0);
}

/* ------ Functions, that create SeqAlignPtr from gap_align_ptr */
/*
	Get the current position.
*/

static Int4 get_current_pos(Int4Ptr pos, Int4 length)
{
    Int4 val;
    if(*pos < 0)
        val = -(*pos + length -1);
    else
        val = *pos;
    *pos += length;
    return val;
}

static SeqAlignPtr MBMakeSeqAlign(SeqIdPtr query_id, SeqIdPtr subject_id, 
                                  Int4 numseg, Int4Ptr length, Int4Ptr start, 
                                  Uint1Ptr strands)
{
    SeqAlignPtr sap;
    DenseSegPtr dsp;
    
    sap = SeqAlignNew();
    
    sap->dim = 2; /**only two dimention alignment**/
    
    /** make the Denseg Object for SeqAlign **/
    sap->segtype = SAS_DENSEG; /** use denseg to store the alignment **/
    sap->type = SAT_PARTIAL;   /* partial for gapped translating search.*/
    dsp = DenseSegNew();
    dsp->dim = 2;
    dsp->numseg = numseg;
    dsp->ids = SeqIdDup(query_id);
    dsp->ids->next = SeqIdDup(subject_id);
    dsp->starts = start;
    dsp->strands = strands;
    dsp->lens = length;
    sap->segs = dsp;
    sap->next = NULL;
    
    return sap;
}

SeqAlignPtr MBCreateSeqAlign(gal_t *galp, SeqIdPtr subject_id, 
                             SeqIdPtr query_id)
{
    Int4 begin1, begin2, index, length1, length2;
    Int4Ptr length, start;
    Uint1Ptr strands;
    Uint1 strand1, strand2;
    Uint4 esp_number, esp_operation;
    edit_script_t *esp;
    Int4 i, numseg, start1, start2;
    SeqAlignPtr sap;
    
    if(galp->S == NULL)
        return NULL;
    
    numseg = galp->S->num;
#if 0
    start1 = galp->b1;
    start2 = galp->b2;

    length1 = galp->e1 - galp->b1;
    length2 = galp->e2 - galp->b2;
    esp = galp->S;

    strand1 = galp->query_strand;
    strand2 = Seq_strand_plus; 
#else
    start1 = galp->b2 - 1;  /* Pos from 1 to 2 eq 0 to 1 in seq-align */ 
    start2 = galp->b1 - 1;
    
    length1 = galp->e2 - galp->b2;
    length2 = galp->e1 - galp->b1;
    esp = galp->S;
    
    strand1 = Seq_strand_plus; 
    /* strand2 = galp->query_strand; */
    strand2 = Seq_strand_plus; 
#endif
    
    start = MemNew((2*numseg+1)*sizeof(Int4));
    length = MemNew((numseg+1)*sizeof(Int4));
    strands = MemNew((2*numseg+1)*sizeof(Uint1));
    
    index = 0;
    for (i = 0; i < numseg; i++) {
        
        esp_number = edit_val_get(esp->op[i]); 
        esp_operation = edit_opc_get(esp->op[i]);
        
        switch(esp_operation) {
        case EDIT_OP_REP:
            if (strand1 != Seq_strand_minus) {
                begin1 = get_current_pos(&start1, esp_number);
            } else {
                begin1 = length1 - 
                    get_current_pos(&start1, esp_number) - esp_number;
            }
            
            if (strand2 != Seq_strand_minus) {
                begin2 = get_current_pos(&start2, esp_number);
            } else {
                begin2 = length2 - 
                    get_current_pos(&start2, esp_number) - esp_number;
            }
            
            strands[2*index] = strand1;
            strands[2*index+1] = strand2;
            start[2*index] = begin1;
            start[2*index+1] = begin2;
            
            break;
            
        case EDIT_OP_DEL:
            begin1 = -1;
            if (strand2 != Seq_strand_minus) {
                begin2 = get_current_pos(&start2, esp_number);
            } else {
                begin2 = length2 - 
                    get_current_pos(&start2, esp_number) - esp_number;
            }
            
            if (index > 0)
                strands[2*index] = strands[2*(index-1)];
            else
                strands[2*index] = Seq_strand_unknown;
            
            strands[2*index+1] = strand2;
            start[2*index] = begin1;
            start[2*index+1] = begin2;
            
            break;
            
        case EDIT_OP_INS:
            if (strand1 != Seq_strand_minus) {
                begin1 = get_current_pos(&start1, esp_number);
            } else {
                begin1 = length1 - 
                    get_current_pos(&start1, esp_number) - esp_number;
            }
            begin2 = -1;
            
            strands[2*index] = strand1;
            if (index > 0)
                strands[2*index+1] = strands[2*(index-1)+1];
            else
                strands[2*index+1] = Seq_strand_unknown;
            
            start[2*index] = begin1;
            start[2*index+1] = begin2;
            
            break;
        default:
            break;
        }

        length[index] = esp_number;
        index++;
    }    

    /* Now building SeqAlign itself */
    
    sap = MBMakeSeqAlign(query_id, subject_id, numseg, length, start, strands);
    
    return sap;
}

BlastSearchBlkPtr GreedyAlignMemAlloc(BlastSearchBlkPtr search)
{
   Int4 max_d, max_d_1, Xdrop, d_diff, max_len, max_cost, gd, i;
   Int4 reward, penalty, gap_open, gap_extend;
   Int4 Mis_cost, GE_cost;
   
   if (search == NULL) 
      return search;
   
   if (search->sbp->reward % 2 == 1) {
      reward = 2*search->sbp->reward;
      penalty = -2*search->sbp->penalty;
      Xdrop = 2*search->pbp->gap_x_dropoff;
      gap_open = 2*search->pbp->gap_open;
      gap_extend = 2*search->pbp->gap_extend;
   } else {
      reward = search->sbp->reward;
      penalty = -search->sbp->penalty;
      Xdrop = search->pbp->gap_x_dropoff;
      gap_open = search->pbp->gap_open;
      gap_extend = search->pbp->gap_extend;
   }
   if (search->rdfp)
      max_len = readdb_get_maxlen(search->rdfp);
   else
      max_len = search->dblen;
   max_d = (Int4) (max_len / ERROR_FRACTION + 1);

   search->abmp = (GreedyAlignMemPtr) Malloc(sizeof(GreedyAlignMem));

   if (search->pbp->gap_open==0 && search->pbp->gap_extend==0) {
      d_diff = ICEIL(Xdrop+reward/2, penalty+reward);
   
      search->abmp->flast_d = (Int4Ptr PNTR) Malloc((max_d + 2) * sizeof(Int4Ptr));
      search->abmp->flast_d[0] = Malloc((max_d + max_d + 6) * sizeof(Int4) * 2);
      search->abmp->flast_d[1] = search->abmp->flast_d[0] + max_d + max_d + 6;
      search->abmp->flast_d_affine = NULL;
      search->abmp->uplow_free = NULL;
   } else {
      search->abmp->flast_d = NULL;
      Mis_cost = reward + penalty;
      GE_cost = gap_extend+reward/2;
      max_d_1 = max_d;
      max_d *= GE_cost;
      max_cost = MAX(Mis_cost, gap_open+GE_cost);
      gd = gdb3(&Mis_cost, &gap_open, &GE_cost);
      d_diff = ICEIL(Xdrop+reward/2, gd);
      search->abmp->uplow_free = MemNew(sizeof(Int4)*2*(max_d+1+max_cost));
      search->abmp->flast_d_affine = (ThreeValPtr PNTR) 
	 Malloc((MAX(max_d, max_cost) + 2) * sizeof(ThreeValPtr));
      search->abmp->flast_d_affine[0] = 
	 MemNew((2*max_d_1 + 6) * sizeof(ThreeVal) * (max_cost+1));
      for (i = 1; i <= max_cost; i++)
	 search->abmp->flast_d_affine[i] = 
	    search->abmp->flast_d_affine[i-1] + 2*max_d_1 + 6;
   }
   search->abmp->max_row_free = Malloc(sizeof(Int4) * (max_d + 1 + d_diff));
   
   return search;
}

GreedyAlignMemPtr GreedyAlignMemFree(GreedyAlignMemPtr abmp)
{
   if (abmp->flast_d) {
      MemFree(abmp->flast_d[0]);
      MemFree(abmp->flast_d);
   } else {
      MemFree(abmp->flast_d_affine[0]);
      MemFree(abmp->flast_d_affine);
      MemFree(abmp->uplow_free);
   }
   MemFree(abmp->max_row_free);
   MemFree(abmp);
   return abmp;
}

/*
	Version to search a packed nucl. sequence against
an unpacked sequence.  s2 is the packed nucl. sequence.
len2 corresponds to the unpacked (true) length.

 * Basic O(ND) time, O(N) space, alignment function. 
 * Parameters:
 *   s1, len1        - first sequence and its length
 *   s2, len2        - second sequence and its length
 *   reverse         - direction of alignment
 *   xdrop_threshold -
 *   mismatch_cost   -
 *   e1, e2          - endpoint of the computed alignment
 *   edit_script     -
 */
Int4 MegaBlastGreedyAlign(const UcharPtr s1, Int4 len1,
			  const UcharPtr s2, Int4 len2,
			  Boolean reverse, Int4 xdrop_threshold, 
			  Int4 match_cost, Int4 mismatch_cost,
			  Int4Ptr e1, Int4Ptr e2, 
			  GreedyAlignMemPtr abmp, edit_script_t *S)
{
    Int4 col,			/* column number */
        d,				/* current distance */
        k,				/* current diagonal */
        flower, fupper,            /* boundaries for searching diagonals */
        row,		        /* row number */
        MAX_D, 			/* maximum cost */
        ORIGIN,
        return_val = 0;
    Int4Ptr PNTR flast_d = abmp->flast_d; /* rows containing the last d */
    Int4Ptr max_row;		/* reached for cost d=0, ... len1.  */
    
    Int4 X_pen = xdrop_threshold;
    Int4 M_half = match_cost/2;
    Int4 Op_cost = mismatch_cost + M_half*2;
    Int4 D_diff = ICEIL(X_pen+M_half, Op_cost);
    
    Int4 x, cur_max, b_diag = 0, best_diag = INT_MAX/2;
    Int4Ptr max_row_free = abmp->max_row_free;
    Char nlower = 0, nupper = 0;
    space_ptr space;
    
    MAX_D = (Int4) (len1/ERROR_FRACTION + 1);
    ORIGIN = MAX_D + 2;
    *e1 = *e2 = 0;
    
    if (reverse)
       /*for (row = 0; row < len2 && row < len1 && (s2[len2-1-row] == READDB_UNPACK_BASE_N(s1[(len1-1-row)>>2], row&3)); row++)*/
	for (row = 0; row < len2 && row < len1 && (s2[len2-1-row] == ncbi4na_to_blastna[s1[len1-1-row]]); row++)
	    /*empty*/ ;
    else
       /*for (row = 0; row < len2 && row < len1 && (s2[row] == READDB_UNPACK_BASE_N(s1[row>>2], (3-(row&3)))); row++)*/
	for (row = 0; row < len2 && row < len1 && (s2[row] == ncbi4na_to_blastna[s1[row]]); row++)
	    /*empty*/ ;
    
    *e1 = row;
    *e2 = row;
    if (row == len1) {
        if (S != NULL)
            edit_script_rep(S, row);
	/* hit last row; stop search */
	return 0;
    }
    if (S==NULL) 
       space = 0;
    else 
       space = new_space(MAX_D);
    
    max_row = max_row_free + X_pen/Op_cost;
    for (k = 0; k < D_diff; k++)
	max_row_free[k] = 0;
    
    flast_d[0][ORIGIN] = row;
    max_row[0] = (row + row)*M_half;
    
    flower = ORIGIN - 1;
    fupper = ORIGIN + 1;

    d = 1;
    while (d <= MAX_D) {
	Int4 fl0, fu0;
	flast_d[d - 1][flower - 1] = flast_d[d - 1][flower] = -1;
	flast_d[d - 1][fupper] = flast_d[d - 1][fupper + 1] = -1;
	x = max_row[d - D_diff] + Op_cost * d - X_pen;
	x = ICEIL(x, M_half);	
	cur_max = 0;
	fl0 = flower;
	fu0 = fupper;
	for (k = fl0; k <= fu0; k++) {
	   row = MAX(flast_d[d - 1][k + 1], flast_d[d - 1][k]) + 1;
	   row = MAX(row, flast_d[d - 1][k - 1]);
	    col = row + k - ORIGIN;
	    if (row + col >= x)
		fupper = k;
	    else {
		if (k == flower)
		    flower++;
		else
		    flast_d[d][k] = -1;
		continue;
	    }
	    /* Slide down the diagonal. Don't do this if reached 
	       the end point, which has value 0x08 */
	    if (reverse) {
	       if (s2[len2 - row] != 0x08) {
		  /*while (row < len2 && col < len1 && s2[len2-1-row] == READDB_UNPACK_BASE_N(s1[(len1-1-col)>>2], col&3)) {*/
		  while (row < len2 && col < len1 && s2[len2-1-row] == ncbi4na_to_blastna[s1[len1-1-col]]) {
		     ++row;
		     ++col;
		  }
	       } else {
		  row--;
		  col--;
	       }
	    } else if (s2[row-1] != 0x08) { 
	       /* while (row < len2 && col < len1 && s2[row] == READDB_UNPACK_BASE_N(s1[col>>2], (3-col&3))) {*/
	       while (row < len2 && col < len1 && s2[row] == ncbi4na_to_blastna[s1[col]]) {
		    ++row;
		    ++col;
		}
	    } else {
	       row--;
	       col--;
	    }
	    flast_d[d][k] = row;
	    if (row + col > cur_max) {
		cur_max = row + col;
		b_diag = k;
	    }
	    if (row == len2) {
		flower = k+1; nlower = 1;
	    }
	    if (col == len1) {
		fupper = k-1; nupper = 1;
	    }
	}
	k = cur_max*M_half - d * Op_cost;
	if (max_row[d - 1] < k) {
	    max_row[d] = k;
	    return_val = d;
	    best_diag = b_diag;
	    *e2 = flast_d[d][b_diag];
	    *e1 = (*e2)+b_diag-ORIGIN;
	} else {
	    max_row[d] = max_row[d - 1];
	}
	if (flower > fupper)
	    break;
	d++;
	if (!nlower) flower--; 
	if (!nupper) fupper++;
	if (S==NULL)
	   flast_d[d] = flast_d[d - 2];
	else 
	   flast_d[d] = (Int4Ptr) get_space(space, fupper-flower+5)-flower+2;
    }
    
    if (S!=NULL) { /*trace back*/
        Int4 row1, col1, diag1, diag;
        d = return_val; diag = best_diag;
        row = *e2; col = *e1;
        while (d > 0) {
            diag1 = get_last(flast_d, d, diag, &row1);
            col1 = row1+diag1-ORIGIN;
            if (diag1 == diag) {
                if (row-row1 > 0) edit_script_rep(S, row-row1);
            } else if (diag1 < diag) {
                if (row-row1 > 0) edit_script_rep(S, row-row1);
                edit_script_ins(S,1);
            } else {
                if (row-row1-1> 0) edit_script_rep(S, row-row1-1);
                edit_script_del(S, 1);
            }
            d--; diag = diag1; col = col1; row = row1;
        }
        edit_script_rep(S, flast_d[0][ORIGIN]);
        if (!reverse) 
            edit_script_reverse_inplace(S);
        free_space(space);
    }
    return return_val;
}


Int4 MegaBlastAffineGreedyAlign (const UcharPtr s1, Int4 len1, 
				 const UcharPtr s2, Int4 len2,
				 Boolean reverse, Int4 xdrop_threshold, 
				 Int4 match_score, Int4 mismatch_score, 
				 Int4 gap_open, Int4 gap_extend,
				 Int4Ptr e1, Int4Ptr e2, 
				 GreedyAlignMemPtr abmp, edit_script_t *S)
{
    Int4 col,			/* column number */
        d,			/* current distance */
        k,			/* current diagonal */
        flower, fupper,         /* boundaries for searching diagonals */
        row,		        /* row number */
        MAX_D, 			/* maximum cost */
        ORIGIN,
        return_val = 0;
    ThreeValPtr PNTR flast_d;	/* rows containing the last d */
    Int4Ptr max_row_free = abmp->max_row_free;
    Int4Ptr max_row;		/* reached for cost d=0, ... len1.  */
    Int4 Mis_cost, GO_cost, GE_cost;
    Int4 D_diff, gd;
    Int4 M_half;
    Int4 max_cost;
    Int4 *lower, *upper;
    
    Int4 x, cur_max, b_diag = 0, best_diag = INT_MAX/2;
    Char nlower = 0, nupper = 0;
    space_ptr space;
    Int4 stop_condition;
    Int4 max_d;
    Int4Ptr uplow_free;
    
    if (match_score % 2 == 1) {
        match_score *= 2;
        mismatch_score *= 2;
        xdrop_threshold *= 2;
        gap_open *= 2;
	gap_extend *= 2;
    }
    M_half = match_score/2;
    if (gap_extend == 0)
       gap_extend = mismatch_score + M_half;
    
    if (gap_open == 0 && gap_extend == mismatch_score+M_half) {
        return MegaBlastGreedyAlign(s1, len1, s2, len2, reverse, 
				    xdrop_threshold, match_score, 
				    mismatch_score, e1, e2, abmp, S);
    }

    Mis_cost = mismatch_score + match_score;
    GO_cost = gap_open;
    GE_cost = gap_extend+M_half;
    gd = gdb3(&Mis_cost, &GO_cost, &GE_cost);
    D_diff = ICEIL(xdrop_threshold+M_half, gd);
    
    
    MAX_D = (Int4) (len1/ERROR_FRACTION + 1);
    max_d = MAX_D*GE_cost;
    ORIGIN = MAX_D + 2;
    max_cost = MAX(Mis_cost, GO_cost+GE_cost);
    *e1 = *e2 = 0;
    
    if (reverse)
       /*for (row = 0; row < len2 && row < len1 && (s2[len2-1-row] == READDB_UNPACK_BASE_N(s1[(len1-1-row)>>2], row&3)); row++)*/
	for (row = 0; row < len2 && row < len1 && (s2[len2-1-row] == ncbi4na_to_blastna[s1[len1-1-row]]); row++)
	    /*empty*/ ;
    else
       /*for (row = 0; row < len2 && row < len1 && (s2[row] == READDB_UNPACK_BASE_N(s1[row>>2], (3-(row&3)))); row++)*/
	for (row = 0; row < len2 && row < len1 && (s2[row] == ncbi4na_to_blastna[s1[row]]); row++)
	    /*empty*/ ;
    
    *e1 = row;
    *e2 = row;
    if (row == len1 || row == len2) {
        if (S != NULL)
            edit_script_rep(S, row);
	/* hit last row; stop search */
	return row*match_score;
    }
    /*flast_d = MemNew((MAX(max_d, max_cost) + 2) * sizeof(ThreeValPtr));*/
    flast_d = abmp->flast_d_affine;
    if (S==NULL) {
        space = 0;
        /*flast_d[0] = MemNew((MAX_D + MAX_D + 6) * sizeof(ThreeVal) * (max_cost+1));
        for (k = 1; k <= max_cost; k++) 
	flast_d[k] = flast_d[k-1] + MAX_D + MAX_D + 6;*/
    } else { 
        space = new_space(MAX_D);
        /*flast_d[0] = get_space(space, 1)-ORIGIN; 
        flast_d[1] = get_space(space, 3)-ORIGIN+1;*/
    }
    max_row = max_row_free + D_diff;
    for (k = 0; k < D_diff; k++)
	max_row_free[k] = 0;
    /*uplow_free = MemNew(sizeof(Int4)*2*(max_d+1+max_cost));*/
    uplow_free = abmp->uplow_free;
    lower = uplow_free;
    upper = uplow_free+max_d+1+max_cost;
    /* next 3 lines set boundary for -1,-2,...,-max_cost+1*/
    for (k = 0; k < max_cost; k++) {lower[k] =LARGE;  upper[k] = -LARGE;}
    lower += max_cost;
    upper += max_cost; 
    
    flast_d[0][ORIGIN].C = row;
    flast_d[0][ORIGIN].I = flast_d[0][ORIGIN].D = -2;
    max_row[0] = (row + row)*M_half;
    lower[0] = upper[0] = ORIGIN;
    
    flower = ORIGIN - 1;
    fupper = ORIGIN + 1;
    
    d = 1;
    stop_condition = 1;
    while (d <= max_d) {
	Int4 fl0, fu0;
	x = max_row[d - D_diff] + gd * d - xdrop_threshold;
	x = ICEIL(x, M_half);
	if (x < 0) x=0;
	cur_max = 0;
	fl0 = flower;
	fu0 = fupper;
	for (k = fl0; k <= fu0; k++) {
	    row = -2;
	    if (k+1 <= upper[d-GO_cost-GE_cost] && k+1 >= lower[d-GO_cost-GE_cost]) 
                row = flast_d[d-GO_cost-GE_cost][k+1].C;
	    if (k+1  <= upper[d-GE_cost] && k+1 >= lower[d-GE_cost] &&
		row < flast_d[d-GE_cost][k+1].D) 
                row = flast_d[d-GE_cost][k+1].D;
	    row++;
	    if (row+row+k-ORIGIN >= x) 
	      flast_d[d][k].D = row;
	    else flast_d[d][k].D = -2;
	    row = -1; 
	    if (k-1 <= upper[d-GO_cost-GE_cost] && k-1 >= lower[d-GO_cost-GE_cost]) 
                row = flast_d[d-GO_cost-GE_cost][k-1].C;
	    if (k-1  <= upper[d-GE_cost] && k-1 >= lower[d-GE_cost] &&
		row < flast_d[d-GE_cost][k-1].I) 
                row = flast_d[d-GE_cost][k-1].I;
	    if (row+row+k-ORIGIN >= x) 
                flast_d[d][k].I = row;
	    else flast_d[d][k].I = -2;
            
	    row = MAX(flast_d[d][k].I, flast_d[d][k].D);
	    if (k <= upper[d-Mis_cost] && k >= lower[d-Mis_cost]) 
                row = MAX(flast_d[d-Mis_cost][k].C+1,row);
            
	    col = row + k - ORIGIN;
	    if (row + col >= x)
		fupper = k;
	    else {
		if (k == flower)
		    flower++;
		else
		    flast_d[d][k].C = -2;
		continue;
	    }
	    /* slide down the diagonal */
	    if (reverse) {
	       if (s2[len2 - row] != 0x08) {
		  /* while (row < len2 && col < len1 && s2[len2-1-row] ==
		     READDB_UNPACK_BASE_N(s1[(len1-1-col)>>2], col&3)) {*/
		  while (row < len2 && col < len1 && s2[len2-1-row] ==
			 ncbi4na_to_blastna[s1[len1-1-col]]) {
		    ++row;
		    ++col;
		  }
	       } else {
		  row--;
		  col--;
	       }
	    } else if (s2[row-1] != 0x08) {
	       /*while (row < len2 && col < len1 && s2[row] == READDB_UNPACK_BASE_N(s1[col>>2], (3-col&3))) {*/
		while (row < len2 && col < len1 && s2[row] == ncbi4na_to_blastna[s1[col]]) {
		    ++row;
		    ++col;
		}
	    } else {
	       row--;
	       col--;
	    }
	    flast_d[d][k].C = row;
	    if (row + col > cur_max) {
		cur_max = row + col;
		b_diag = k;
	    }
	    if (row == len2) {
		flower = k; nlower = k+1;
	    }
	    if (col == len1) {
		fupper = k; nupper = k-1;
	    }
	}
	k = cur_max*M_half - d * gd;
	if (max_row[d - 1] < k) {
	    max_row[d] = k;
	    return_val = d;
	    best_diag = b_diag;
	    *e2 = flast_d[d][b_diag].C;
	    *e1 = (*e2)+b_diag-ORIGIN;
	} else {
	    max_row[d] = max_row[d - 1];
	}
	if (flower <= fupper) {
            stop_condition++;
            lower[d] = flower; upper[d] = fupper;
	} else {
            lower[d] = LARGE; upper[d] = -LARGE;
	}
	if (lower[d-max_cost] <= upper[d-max_cost]) stop_condition--;
	if (stop_condition == 0) break;
	d++;
	flower = MIN(lower[d-Mis_cost], MIN(lower[d-GO_cost-GE_cost], lower[d-GE_cost])-1);
	if (nlower) flower = MAX(flower, nlower);
	fupper = MAX(upper[d-Mis_cost], MAX(upper[d-GO_cost-GE_cost], upper[d-GE_cost])+1);
	if (nupper) fupper = MIN(fupper, nupper);
	if (d > max_cost) {
	   if (S==NULL) {
	      /*if (d > max_cost)*/
	      flast_d[d] = flast_d[d - max_cost-1];
	   } else 
	      flast_d[d] = get_space(space, fupper-flower+1)-flower;
	}
    }
    
    if (S!=NULL) { /*trace back*/
        Int4 row1, diag, state;
        d = return_val; diag = best_diag;
        row = *e2; state = sC;
        while (d > 0) {
            if (state == sC) {
                /* diag will not be changed*/
                state = get_lastC(flast_d, lower, upper, &d, diag, Mis_cost, &row1);
                if (row-row1 > 0) edit_script_rep(S, row-row1);
                row = row1;
            } else {
                if (state == sI) {
                    /*row unchanged */
                    state = get_lastI(flast_d, lower, upper, &d, diag, GO_cost, GE_cost);
                    diag--;
                    edit_script_ins(S,1);
                } else {
                    edit_script_del(S,1);
                    state = get_lastD(flast_d, lower, upper, &d, diag, GO_cost, GE_cost);
                    diag++;
                    row--;
                }
            }
        }
        edit_script_rep(S, flast_d[0][ORIGIN].C);
        if (!reverse) 
            edit_script_reverse_inplace(S);
        free_space(space);
    }/* else {
        MemFree(flast_d[0]);
	}*/
    return_val = max_row[return_val];
    /*MemFree(flast_d);
      MemFree(uplow_free);*/
    return return_val;
}



/* ------------------------------------------------------------ */
