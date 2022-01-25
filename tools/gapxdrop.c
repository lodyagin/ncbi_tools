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

File name: gapxdrop.c

Author: Gennadiy Savchuk, Jinqhui Zhang, Tom Madden
[From code originally supplied by Webb Miller's lab].

Contents: Functions to perform a gapped alignment on two sequences.

****************************************************************************/
/* $Revision: 6.18 $ 
* $Log: gapxdrop.c,v $
* Revision 6.18  1999/08/27 18:07:34  shavirin
* Passed parameter decline_align from top to the engine.
*
* Revision 6.17  1999/05/12 18:48:52  madden
* lift assumption that no deletion follows insertion
*
* Revision 6.16  1999/05/03 18:59:33  madden
* Removed set, but unused, variable count
*
* Revision 6.15  1999/03/17 18:39:11  madden
* Removed unneccessary memset
*
* Revision 6.14  1999/02/18 21:18:49  madden
* MINIINT is INT4_MIN/2
*
* Revision 6.13  1999/02/17 19:42:42  madden
* change INT4_MIN to -9999999 again
*
* Revision 6.11  1998/12/18 16:20:42  madden
* Removed unnecessary memsets
*
 * Revision 6.10  1998/11/19 14:03:38  madden
 * minor efficiency
 *
 * Revision 6.9  1998/11/17 13:39:03  madden
 * Made ALIGN non-static
 *
 * Revision 6.8  1998/09/24 15:26:40  egorov
 * Fix lint complaints
 *
 * Revision 6.7  1998/08/05 21:28:16  madden
 * Protection against core-dump when gap_extend is zero
 *
 * Revision 6.6  1998/05/22 19:16:52  egorov
 * Another try to fix time problem (increase chunk size to 2Mb)
 *
 * Revision 6.5  1998/05/21 12:49:26  egorov
 * Remove memory allocation optimization to reduce time
 *
 * Revision 6.4  1998/04/21 18:39:25  madden
 * Zhengs suggestion to save memory on long sequences
 *
 * Revision 6.3  1998/04/17 19:41:17  madden
 * Zhengs changes for decline to align
 *
 * Revision 6.2  1998/01/06 18:26:48  madden
 * Gaps have same strand as surrounding sequence
 *
 * Revision 6.1  1997/09/22 17:36:29  madden
 * MACROS for position-specific matrices from Andy Neuwald
 *
 * Revision 6.0  1997/08/25 18:53:14  madden
 * Revision changed to 6.0
 *
 * Revision 1.32  1997/06/24 21:03:38  madden
 * TracebackToGapXEditBlock() check for last GapXEditScript fixed
 *
 * Revision 1.31  1997/06/24 19:12:12  tatiana
 * TracebackToGapXEditBlock() check for last GapXEditScript added
 *
 * Revision 1.30  1997/06/24 18:14:33  tatiana
 * TracebackToGapXEditBlock() fixed
 *
 * Revision 1.29  1997/05/13 20:45:50  madden
 * fixed offset problem in PerformGappedAlignment
 *
 * Revision 1.28  1997/04/21 15:35:26  madden
 * Fixes for 'gapped' StdSegs.
 *
 * Revision 1.27  1997/04/16  16:21:38  madden
 * Set seqalign_type.
 *
 * Revision 1.26  1997/04/15  22:01:53  madden
 * Added original_length[12] for translating searches.
 *
 * Revision 1.25  1997/04/14  14:10:37  madden
 * Fix for the case of dropoff less than the opening penalty.
 *
 * Revision 1.24  1997/04/11  19:02:49  madden
 * Changes for in-frame blastx, tblastn gapping.
 *
 * Revision 1.23  1997/03/14  21:01:59  madden
 * Changed to use less memory in ALIGN, with GapXDropStateArrayStructPtr.
 *
 * Revision 1.22  1997/03/11  19:24:08  madden
 * Check return value of MemNew for state_array in ALIGN.
 *
 * Revision 1.21  1997/03/01  18:25:33  madden
 * Sequences reversed on SeqAlign if reverse flag set on EditBlockPtr.
 *
 * Revision 1.20  1997/02/24  16:40:38  madden
 * Change to GapXEditBlockToSeqAlign to use first SeqIdPtr, duplicate.
 *
 * Revision 1.19  1997/02/24  15:09:38  madden
 * Replaced MemSet by setting some elements of an array to zero,.
 *
 * Revision 1.18  1997/02/23  16:44:47  madden
 * Memory, saved on GapAlignBlkPtr, reuses.
 *
 * Revision 1.17  1997/02/20  22:58:49  madden
 * Support for gapped translated results.
 *
 * Revision 1.16  1997/02/20  21:50:24  madden
 * Added frame and translation information to GapAlignBlk, assigned it.
 *
 * Revision 1.15  1997/02/12  22:19:08  madden
 * Changes for position-based comparisons.
 *
 * Revision 1.14  1997/02/10  15:25:33  madden
 * Changed args. of ALIGN and SEMI_G_ALIGN to pass in gap_align block,
 * also whether the query was reversed and the query_offset.
 * Also General cleanup.
 *
 * Revision 1.13  1997/02/04  16:22:32  madden
 * Changes to enable gapped alignments on the reverse strand.
 *
 * Revision 1.12  1997/01/22  17:45:08  madden
 * Added positionBased Boolean to ALIGN and SEMI_G_ALIGN.
 *
 * Revision 1.11  1997/01/21  18:55:17  madden
 * Fixed bug with translation of EditScript to SeqAlign.
 *
 * Revision 1.10  1997/01/17  14:29:34  madden
 * Protection against NULL args added in GapXEditBlockDelete.
 *
 * Revision 1.9  1997/01/16  20:20:49  madden
 * TracebackToGapXEditBlock made non-static.
 *
 * Revision 1.8  1997/01/06  22:40:55  madden
 * Added function SimpleIntervalToGapXEditBlock.
 *
 * Revision 1.7  1997/01/06  17:22:59  madden
 * Used GapXEditScriptToSeqAlign to find SeqAlign.
 *
 * Revision 1.6  1996/12/30  21:45:28  madden
 * UMR fix.
 *
 * Revision 1.5  1996/12/30  17:14:06  madden
 * Fixes for changes for "require a portion of the query sequence".
 *
 * Revision 1.4  1996/12/30  15:44:25  madden
 * Added capability to require a portion of the query sequence.
 *
 * Revision 1.3  1996/12/16  15:29:12  madden
 * Changed gapalign.h to gapxdrop.h
 *
 * Revision 1.2  1996/12/12  16:45:03  madden
 * GapAlignBlkPtr used instead of arguments in functions.
 *
 * Revision 1.1  1996/12/12  14:02:51  madden
 * Initial revision
 *
*/


#include "gapxdrop.h"


/* A PACKAGE FOR LOCALLY ALIGNING TWO SEQUENCES WITHIN A BAND:

   To invoke, call SEMI_G_ALIGN(A,B,M,N,W,G,H,S,X,&EI,&EJ, score_only, positionBased).
   The parameters are explained as follows:
	A, B : two sequences to be aligned
	M : the length of sequence A
	N : the length of sequence B
	W : scoring table for matches and mismatches
	G : gap-opening penalty
	H : gap-extension penalty
	X : maximum drop off
	S : script for DISPLAY routine
	*EI : ending position of sequence A in the optimal local alignment
	*EJ : ending position of sequence B in the optimal local alignment
	score_only: indicate if ony score is needed. 1--score only 0--alignment
	            is also needed.

	Use best_score to cut unfeasible points. quadratic space and row-wise
*/


/* Append "Delete k" op */
#define DEL_(k) \
data.last = (data.last < 0) ? (data.sapp[-1] -= (k)) : (*data.sapp++ = -(k));

/* Append "Insert k" op */
#define INS_(k) \
data.last = (data.last > 0) ? (data.sapp[-1] += (k)) : (*data.sapp++ = (k));

/* Append "Replace" op */
#define REP_ \
{ data.last = *data.sapp++ = 0; }

/* Divide by two to prevent underflows. */
#define MININT INT4_MIN/2
#define REPP_ \
{ *data.sapp++ = MININT; data.last = 0; }


/* Dynamic Programming structure. */
typedef struct DP {
  Int4 CC, DD, FF;	/* Values for gap opening and extensions (?). */
} PNTR dp_ptr, dp_node;

typedef struct {
  Int4Ptr sapp;			/* Current script append ptr */
  Int4  last;	
} data_t;

/* #define	CHUNKSIZE	1048576 */
#define	CHUNKSIZE	2097152

static GapXDropStateArrayStructPtr
GapXDropGetState(GapXDropStateArrayStructPtr PNTR head, Int4 length)

{
	GapXDropStateArrayStructPtr	retval, var, last;
	Int4				chunksize = MAX(CHUNKSIZE, length + length/3);

	length += length/3;	/* Add on about 30% so the end will get reused. */
	retval = NULL;
	if (*head == NULL)
	{
		retval = (GapXDropStateArrayStructPtr) MemNew(sizeof(GapXDropStateArrayStruct));
		retval->state_array = MemNew(chunksize*sizeof(Uint1));
		retval->length = chunksize;
		retval->used = 0;
		*head = retval;
	}
	else
	{
		var = *head;
		last = *head;
		while (var)
		{
			if (length < (var->length - var->used))
			{
				retval = var;
				break;
			}
			else if (var->used == 0)
			{ /* If it's empty and too small, replace. */
				var->state_array = MemFree(var->state_array);
				var->state_array = MemNew(chunksize*sizeof(Uint1));
				var->length = chunksize;
				retval = var;
				break;
			}
			last = var;
			var = var->next;
		}

		if (var == NULL)
		{
			retval = (GapXDropStateArrayStructPtr) MemNew(sizeof(GapXDropStateArrayStruct));
			retval->state_array = MemNew(chunksize*sizeof(Uint1));
			retval->length = chunksize;
			retval->used = 0;
			last->next = retval;
		}
	}

	if (retval->state_array == NULL)
		ErrPostEx(SEV_ERROR, 0, 0, "state array is NULL");
		
	return retval;

}

static Boolean
GapXDropPurgeState(GapXDropStateArrayStructPtr state_struct)

{
	while (state_struct)
	{
		MemSet(state_struct->state_array, 0, state_struct->used);
		state_struct->used = 0;
		state_struct = state_struct->next;
	}

	return TRUE;
}

GapXDropStateArrayStructPtr 
GapXDropStateDestroy(GapXDropStateArrayStructPtr state_struct)

{
	GapXDropStateArrayStructPtr next;

	while (state_struct)
	{
		next = state_struct->next;
		MemFree(state_struct->state_array);
		MemFree(state_struct);
		state_struct = next;
		
	}

	return NULL;
}

Int4 LIBCALL 
ALIGN(Uint1Ptr A, Uint1Ptr B, Int4 M, Int4 N,
		Int4Ptr S, Int4Ptr pei, Int4Ptr pej, Int4Ptr PNTR sapp, 
		GapAlignBlkPtr gap_align, Int4 query_offset, Boolean reversed)
	
{ 
  data_t data;
  Int4 i, j, cb,  j_r, s, k;
  Uint1 st, std, ste;
  Int4 gap_open, gap_extend, decline_penalty;
  register Int4 c, d, e, m,t, tt, f, tt_start;
  Int4 best_score = 0;
  register Int4Ptr wa;
  register dp_ptr dp, dyn_prog;
  Uint1Ptr PNTR state, stp, tmp;
  Uint1Ptr state_array;
  Int4Ptr *matrix;
  register Int4 X;
  GapXDropStateArrayStructPtr state_struct;
  
  matrix = gap_align->matrix;
  *pei = *pej = 0;
  data.sapp = *sapp = S;
  data.last= 0;
  m = gap_align->gap_open + gap_align->gap_extend;
  decline_penalty = gap_align->decline_align;
  gap_open = gap_align->gap_open;
  gap_extend = gap_align->gap_extend;
  X = gap_align->x_parameter;

  if (X < m)
	X = m;

  if(N <= 0 || M <= 0) { 
    *pei = *pej;
    return 0;
  }

  GapXDropPurgeState(gap_align->state_struct);

  j = (N + 2) * sizeof(dp_node);
  dyn_prog = (dp_ptr)Nlm_Malloc(j);

  state = Nlm_Malloc(sizeof(Uint1Ptr)*(M+1));
  dyn_prog[0].CC = 0;
  dyn_prog[0].FF = -m;
  c = dyn_prog[0].DD = -m;

  /* Protection against divide by zero. */
  if (gap_extend > 0)
  	state_struct = GapXDropGetState(&gap_align->state_struct, X/gap_extend+5);
  else
	state_struct = GapXDropGetState(&gap_align->state_struct, N+3);

  state_array = state_struct->state_array;
  state[0] = state_array;
  stp  = state[0];
  for(i = 1; i <= N; i++) {
    if(c < -X) break;
    dyn_prog[i].CC = c;
    dyn_prog[i].DD = c - m; 
    dyn_prog[i].FF = c-m;
    c -= gap_extend;
    stp[i] = 1;
  }
  state_struct->used = i+1;
  
  tt = 0;  j = i;
  for(j_r = 1; j_r <= M; j_r++) {
     /* Protection against divide by zero. */
    if (gap_extend > 0)
    	state_struct = GapXDropGetState(&gap_align->state_struct, j-tt+5+X/gap_extend);
    else
	state_struct = GapXDropGetState(&gap_align->state_struct, N-tt+3);
    state_array = state_struct->state_array + state_struct->used + 1;
    state[j_r] = state_array - tt + 1;
    stp = state[j_r];
    tt_start = tt; 
    if (!(gap_align->positionBased)) /*AAS*/
      wa = matrix[A[j_r]]; 
    else {
      if(reversed)
        wa = gap_align->posMatrix[M - j_r];
      else
        wa = gap_align->posMatrix[j_r + query_offset];
    }
    e = c = f= MININT;
    
    for (cb = i = tt, dp = &dyn_prog[i-1]; i < j; i++) {
	d = (++dp)->DD;
	st = 0;
	if (c < f) {c = f; st = 3;}
	if (f > d) {d = f; std = 60;} 
	else {
	  std = 30;
	  if (c < d) { c= d;st = 2;}
	}
	if (f > e) {e = f; ste = 20;} 
	else {
	  ste = 10;
	  if (c < e) {c=e; st=1;}
	}
	if (best_score - c > X){
	  c = dp->CC+wa[B[i+1]]; f= dp->FF;
	  if (tt == i) tt++;
	  else { dp->CC = dp->FF =  MININT; }
	} else {
	  cb = i;
	  if (c > best_score) {
	    best_score = c;
	    *pei = j_r; *pej = i;
	  }
	  if ((c-=m) > (d-=gap_extend)) {
	    dp->DD = c; 
	  } else {
	    dp->DD = d;
	    st+=std;
	  } 
	  if (c > (e-=gap_extend)) {
	    e = c; 
	  }  else {
	    st+=ste;
	  }
	  c+=m; 
	  d = dp->FF;
	  if (f < c-gap_open) { 
	    dp->FF = c-gap_open-decline_penalty; 
	  } else {
	    dp->FF = f-decline_penalty; st+= 5; 
	  }
	  f = d;
	  d = dp->CC+wa[B[i+1]]; dp->CC = c; c = d;      
	}
	stp[i] = st;
    }
    if (tt == j) { j_r++; break;}
    if (cb < j-1) { j = cb+1;}
    else {
	while (e >= best_score-X && j <= N) {
	    dyn_prog[j].CC = e; dyn_prog[j].DD = e-m; dyn_prog[j].FF = e-gap_open-decline_penalty;
	    e -= gap_extend; stp[j] = 1;
	    j++;
	}
    }
    if (j <= N) {
	dyn_prog[j].DD = dyn_prog[j].CC= dyn_prog[j].FF = MININT;
	j++; 
    }
    state_struct->used += (MAX(i, j) - tt_start + 1);
  }
  i = *pei; j = *pej;
  tmp = Nlm_Malloc(i+j);
  for (s=0, c = 0; i> 0 || j > 0; c++) {
      t = state[i][j];
      k  = t %5;
      if (s == 1) 
	  if ((t/10)%3 == 1) k = 1;
	  else if ((t/10)%3 == 2) k = 3;
      if (s == 2)
	  if ((t/30) == 1) k = 2;
	  else if ((t/30) == 2) k = 3;
      if (s == 3 && ((t/5)%2) == 1) k = 3;
      if (k == 1) { j--;}
      else if (k == 2) {i--;}
      else {j--; i--;}
      tmp[c] = s = k;
  }
  c--;
  while (c >= 0) {
      if (tmp[c] == 0) REP_
      else if (tmp[c] == 1) INS_(1)
      else if (tmp[c] == 3) REPP_			  
      else DEL_(1)     
      c--;
  }

  MemFree(tmp);

  MemFree(state);

  MemFree(dyn_prog);
  *sapp = data.sapp;

  return best_score;
}

Int4 SEMI_G_ALIGN(Uint1Ptr A, Uint1Ptr B, Int4 M, Int4 N,
		Int4Ptr S, Int4Ptr pei, Int4Ptr pej, 
		Boolean score_only, Int4Ptr PNTR sapp, GapAlignBlkPtr gap_align,
		Int4 query_offset, Boolean reversed)
		
{ 
  dp_ptr dyn_prog;
  Int4 i, j, cb, j_r, g, decline_penalty;
  register Int4 c, d, e, m, tt, h, X, f;
  Int4 best_score = 0;
  Int4Ptr *matrix;
  register Int4Ptr wa;
  register dp_ptr dp;
  
  if(!score_only)
    return ALIGN(A, B, M, N, S, pei, pej, sapp, gap_align, query_offset, reversed);
  
  matrix = gap_align->matrix;
  *pei = *pej = 0;
  m = (g=gap_align->gap_open) + gap_align->gap_extend;
  h = gap_align->gap_extend;
  decline_penalty = gap_align->decline_align;
  X = gap_align->x_parameter;

  if (X < m)
	X = m;

  if(N <= 0 || M <= 0) return 0;

  j = (N + 2) * sizeof(dp_node);
  dyn_prog = (dp_ptr)Nlm_Malloc(j);

  dyn_prog[0].CC = 0; c = dyn_prog[0].DD = -m;
  dyn_prog[0].FF = -m;
  for(i = 1; i <= N; i++) {
    if(c < -X) break;
    dyn_prog[i].CC = c;
    dyn_prog[i].DD = c - m; 
    dyn_prog[i].FF = c-m;
    c -= h;
  }

  tt = 0;  j = i;
  for (j_r = 1; j_r <= M; j_r++) {
      if (!(gap_align->positionBased)) /*AAS*/
	  wa = matrix[A[j_r]]; 
      else {
	  if(reversed)
	      wa = gap_align->posMatrix[M - j_r];
	  else
	      wa = gap_align->posMatrix[j_r + query_offset];
      }
      e = c =f = MININT;
      for (cb = i = tt, dp = &dyn_prog[i]; i < j; i++) {
	  d = dp->DD;
	  if (e < f) e = f;
	  if (d < f) d = f;
	  if (c < d || c < e) {
	      if (d < e) {
		  c = e;
	      } else {
		  c = d; 
	      }
	      if (best_score - c > X) {
		  c = dp->CC+wa[B[i+1]]; f = dp->FF;
		  if (tt == i) tt++;
		  else { dp->CC =dp->FF= MININT;}
	      } else {
		  cb = i;
                  if ((c-=m) > (d-=h)) {
                      dp->DD = c;
                  } else {
                      dp->DD = d;
                  }
                  if (c > (e-=h)) {
                      e = c;
                  }
                  c+=m;
		  d = dp->CC+wa[B[i+1]]; dp->CC = c; c=d;
		  d = dp->FF; dp->FF = f-decline_penalty; f = d;
	      }
	  } else {
	      if (best_score - c > X){
		  c = dp->CC+wa[B[i+1]]; f= dp->FF;
		  if (tt == i) tt++;
		  else { dp->CC =dp->FF= MININT;}
	      } else {
		  cb = i; 
		  if (c > best_score) {
		      best_score = c;
		      *pei = j_r; *pej = i;
		  } 
		  if ((c-=m) > (d-=h)) {
		      dp->DD = c; 
		  } else {
		      dp->DD = d;
		  } 
		  if (c > (e-=h)) {
		      e = c;
		  } 
		  c+=m;
		  d = dp->FF;
		  if (c-g>f) dp->FF = c-g-decline_penalty; else dp->FF = f-decline_penalty;
		  f = d;
		  d = dp->CC+wa[B[i+1]]; dp->CC = c; c = d;
	      }
	  }
	  dp++;
      }
      if (tt == j) break;
      if (cb < j-1) { j = cb+1;}
      else while (e >= best_score-X && j <= N) {
	  dyn_prog[j].CC = e; dyn_prog[j].DD = e-m;dyn_prog[j].FF = e-g-decline_penalty;
	  e -= h; j++;
      }
      if (j <= N) {
	  dyn_prog[j].DD = dyn_prog[j].CC = dyn_prog[j].FF = MININT; j++;
      }
  }
  

  MemFree(dyn_prog);

  return best_score;
}


/*
	Allocates an EditScriptPtr and puts it on the end of
	the chain of EditScriptPtr's.  Returns a pointer to the
	new one.

*/
static GapXEditScriptPtr 
GapXEditScriptNew(GapXEditScriptPtr old)

{
	GapXEditScriptPtr new;

	new = (GapXEditScriptPtr) MemNew(sizeof(GapXEditScript));

	if (old == NULL)
		return new;

	while (old->next != NULL)
	{
		old = old->next;
	}

	old->next = new;

	return new;
}

static GapXEditScriptPtr
GapXEditScriptDelete(GapXEditScriptPtr old)
{
	GapXEditScriptPtr next;

	while (old)
	{
		next = old->next;
		old = MemFree(old);
		old = next;
	}
	return old;
}

GapXEditBlockPtr LIBCALL
GapXEditBlockNew(Int4 start1, Int4 start2)

{
	GapXEditBlockPtr edit_block;

	edit_block = (GapXEditBlockPtr) MemNew(sizeof(GapXEditBlock));
	edit_block->start1 = start1;
	edit_block->start2 = start2;

	return edit_block;
}

GapXEditBlockPtr LIBCALL
GapXEditBlockDelete(GapXEditBlockPtr edit_block)

{
	if (edit_block == NULL)
		return NULL;

	edit_block->esp = GapXEditScriptDelete(edit_block->esp);

	edit_block = MemFree(edit_block);

	return edit_block;
}

/* Alignment display routine */

GapXEditBlockPtr LIBCALL
TracebackToGapXEditBlock(Uint1Ptr A, Uint1Ptr B, Int4 M, Int4 N, Int4Ptr S, Int4 start1, Int4 start2)
{

  register Int4 i, j, op, number_of_subs, number_of_decline;
  GapXEditScriptPtr e_script;
  GapXEditBlockPtr edit_block;

  if (S == NULL)
	return NULL;

  i = j = op = number_of_subs = number_of_decline = 0;

  edit_block = GapXEditBlockNew(start1, start2);
  edit_block->esp = e_script = GapXEditScriptNew(NULL);

  while(i < M || j < N) 
  {
	op = *S;
	if (op != MININT && number_of_decline > 0) 
	{
               e_script->op_type = GAPALIGN_DECLINE;
               e_script->num = number_of_decline;
               e_script = GapXEditScriptNew(e_script);
		number_of_decline = 0;
	}
        if (op != 0 && number_of_subs > 0) 
	{
                        e_script->op_type = GAPALIGN_SUB;
                        e_script->num = number_of_subs;
                        e_script = GapXEditScriptNew(e_script);
                        number_of_subs = 0;
        }      
	if (op == 0) {
		i++; j++; number_of_subs++;
	} else if (op == MININT) {
		i++; j++;
		number_of_decline++; 
	}	
	else
	{
		if(op > 0) 
		{
			e_script->op_type = GAPALIGN_DEL;
			e_script->num = op;
			j += op;
			if (i < M || j < N)
                                e_script = GapXEditScriptNew(e_script);
		}
		else
		{
			e_script->op_type = GAPALIGN_INS;
			e_script->num = ABS(op);
			i += ABS(op);
			if (i < M || j < N)
                                e_script = GapXEditScriptNew(e_script);
		}
    	}
	S++;
  }

  if (number_of_subs > 0)
  {
	e_script->op_type = GAPALIGN_SUB;
	e_script->num = number_of_subs;
  } else if (number_of_decline > 0) {
        e_script->op_type = GAPALIGN_DECLINE;
        e_script->num = number_of_decline;
  }

  return edit_block;


}

static Int4
reverse_seq(Uint1Ptr seq, Uint1Ptr pos, Uint1Ptr target) 
{
    Uint1Ptr c;

    for (c = pos; c >= seq; c--) 
	*target++ = *c;
    *target = '\0';

    return (Int4) (pos-c);
}
/*
	Destruct Function for GapAlignBlk.  If "state" is not NULL, then
	it's deallocated.
*/
GapAlignBlkPtr LIBCALL
GapAlignBlkDelete(GapAlignBlkPtr gap_align)

{
	if (gap_align == NULL)
		return NULL;

	gap_align->state_struct = GapXDropStateDestroy(gap_align->state_struct);
	gap_align = MemFree(gap_align);

	return gap_align;
}

/*
	Allocates GapAlignBlkPtr and "state", if state_column_length and 
	state_row_length are not NULL.

	For "call and forget" applications, state_column_length and
	state_row_length should both be set to zero.  ALIGN will
	then allocate and deallocate this memory.
*/

GapAlignBlkPtr LIBCALL
GapAlignBlkNew(Int4 state_column_length, Int4 state_row_length)

{
	GapAlignBlkPtr gap_align;

	gap_align = MemNew(sizeof(GapAlignBlk));
	if (gap_align == NULL)
		return NULL;

	/* gap_align->decline_align = INT2_MAX; */
	return gap_align;
}


Boolean LIBCALL
PerformGappedAlignment(GapAlignBlkPtr gap_align)

{
	Boolean found_start, found_end;
	Int4 q_length=0, s_length=0, middle_score, score_right, score_left, private_q_start, private_s_start;
	Int4 include_query, index;
	Uint1Ptr q_left=NULL, s_left=NULL;
	Uint1Ptr query, subject, query_var, subject_var;

	if (gap_align == NULL)
		return FALSE;

	found_start = FALSE;
	found_end = FALSE;

	query = gap_align->query;
	subject = gap_align->subject;
	include_query = gap_align->include_query;
	score_left = 0;
	if (gap_align->q_start != 0 && gap_align->s_start != 0)
	{
		found_start = TRUE;
		q_left = (Uint1Ptr) Nlm_Malloc((gap_align->q_start+2)*sizeof(Uint1));
		s_left = (Uint1Ptr) Nlm_Malloc((gap_align->s_start+2)*sizeof(Uint1));

		q_length = reverse_seq(query, query+gap_align->q_start, q_left);
		s_length = reverse_seq(subject, subject+gap_align->s_start, s_left);

		score_left = SEMI_G_ALIGN(q_left-1, s_left-1, q_length, s_length, NULL, &private_q_start, &private_s_start, TRUE, NULL, gap_align, gap_align->q_start, TRUE);
		gap_align->query_start = q_length - private_q_start;
		gap_align->subject_start = s_length - private_s_start;
	}
	else
	{
		q_length = gap_align->q_start;
		s_length = gap_align->s_start;
	}

	middle_score = 0;
	query_var = query+gap_align->q_start;
	subject_var = subject+gap_align->s_start;
	for (index=0; index<include_query; index++)
	{
		query_var++;
		subject_var++;
		if (!(gap_align->positionBased))  /*AAS*/
		  middle_score += gap_align->matrix[*query_var][*subject_var];
		else 
		  middle_score += MtrxScoreGapAlign(gap_align,
				gap_align->q_start+1 + index,*subject_var);
	}

	score_right = 0;
	if (gap_align->q_start+include_query < gap_align->query_length && gap_align->s_start+include_query < gap_align->subject_length)
	{
		found_end = TRUE;
		score_right = SEMI_G_ALIGN(query+gap_align->q_start+include_query, subject+gap_align->s_start+include_query, gap_align->query_length-q_length-include_query, gap_align->subject_length-s_length-include_query, NULL, &(gap_align->query_stop), &(gap_align->subject_stop), TRUE, NULL, gap_align, gap_align->q_start+include_query, FALSE);
		gap_align->query_stop += gap_align->q_start+include_query;
		gap_align->subject_stop += gap_align->s_start+include_query;
	}

	if (found_start == FALSE)
	{	/* Start never found */
		gap_align->query_start = gap_align->q_start;
		gap_align->subject_start = gap_align->s_start;
	}

	if (found_end == FALSE)
	{
		gap_align->query_stop = gap_align->q_start + include_query - 1;
		gap_align->subject_stop = gap_align->s_start + include_query - 1;
	}

	q_left = MemFree(q_left);
	s_left = MemFree(s_left);

	gap_align->score = score_right+score_left+middle_score;

	return TRUE;
}

/*
	Perform a gapped alignment and return the score.  A SeqAlignPtr with
	the traceback is also returned.
*/
Boolean LIBCALL
PerformGappedAlignmentWithTraceback(GapAlignBlkPtr gap_align)

{
	Boolean found_start, found_end;
	Int4 q_length=0, s_length=0, score_right, middle_score, score_left, private_q_length, private_s_length, tmp;
	Int4 include_query, index;
	Int4Ptr tback, tback1, p, q;
	Uint1Ptr q_left=NULL, s_left=NULL;
	Uint1Ptr query, subject, query_var, subject_var;

	if (gap_align == NULL)
		return FALSE;

	found_start = FALSE;
	found_end = FALSE;

	query = gap_align->query;
	subject = gap_align->subject;

	tback = tback1 = Nlm_Malloc((gap_align->subject_length+gap_align->query_length)*sizeof(Int4));
	include_query = gap_align->include_query;

	score_left = 0;
	if (gap_align->q_start != 0 && gap_align->s_start != 0)
	{
		found_start = TRUE;
		q_left = (Uint1Ptr) Nlm_Malloc((gap_align->q_start+2)*sizeof(Uint1));
		s_left = (Uint1Ptr) Nlm_Malloc((gap_align->s_start+2)*sizeof(Uint1));

		q_length = reverse_seq(query, query+gap_align->q_start, q_left);
		s_length = reverse_seq(subject, subject+gap_align->s_start, s_left);

		score_left = SEMI_G_ALIGN(q_left-1, s_left-1, q_length, s_length, tback, &private_q_length, &private_s_length, FALSE, &tback1, gap_align, gap_align->q_start, TRUE);
	        for(p = tback, q = tback1 - 1; p < q; p++, q--) 
		{
		        tmp = *p;
		  	*p = *q;
			*q = tmp;
        	}
		gap_align->query_start = q_length - private_q_length;
		gap_align->subject_start = s_length - private_s_length;
	}
	else
	{
		q_length = gap_align->q_start;
		s_length = gap_align->s_start;
	}

	middle_score = 0;
	query_var = query+gap_align->q_start;
	subject_var = subject+gap_align->s_start;
	for (index=0; index<include_query; index++)
	{
		query_var++;
		subject_var++;
		if (!(gap_align->positionBased))  /*AAS*/
		  middle_score += gap_align->matrix[*query_var][*subject_var];
		else 
		  middle_score += MtrxScoreGapAlign(gap_align,
			gap_align->q_start+1 + index,*subject_var);
		*tback1 = 0;
		tback1++;
	}
	
	score_right = 0;
	if ((gap_align->q_start+include_query) < gap_align->query_length && (gap_align->s_start+include_query) < gap_align->subject_length)
	{
		found_end = TRUE;
		score_right = SEMI_G_ALIGN(query+gap_align->q_start+include_query, subject+gap_align->s_start+include_query, gap_align->query_length-q_length-include_query, gap_align->subject_length-s_length-include_query, tback1, &private_q_length, &private_s_length, FALSE, &tback1, gap_align, gap_align->q_start+include_query, FALSE);
		gap_align->query_stop = gap_align->q_start + private_q_length+include_query;
		gap_align->subject_stop = gap_align->s_start + private_s_length+include_query;
	}

	if (found_start == FALSE)
	{	/* Start never found */
		gap_align->query_start = gap_align->q_start;
		gap_align->subject_start = gap_align->s_start;
	}

	if (found_end == FALSE)
	{
		gap_align->query_stop = gap_align->q_start + include_query - 1;
		gap_align->subject_stop = gap_align->s_start + include_query - 1;
	}

	gap_align->edit_block = TracebackToGapXEditBlock(query, subject, gap_align->query_stop-gap_align->query_start+1, gap_align->subject_stop-gap_align->subject_start+1, tback, gap_align->query_start, gap_align->subject_start);

	gap_align->edit_block->frame1 = gap_align->query_frame;
	gap_align->edit_block->frame2 = gap_align->subject_frame;
	gap_align->edit_block->length1 = gap_align->query_length;
	gap_align->edit_block->length2 = gap_align->subject_length;
	gap_align->edit_block->translate1 = gap_align->translate1;
	gap_align->edit_block->translate2 = gap_align->translate2;

	tback = MemFree(tback);
	q_left = MemFree(q_left);
	s_left = MemFree(s_left);

	gap_align->score = score_right+score_left+middle_score;

	return TRUE;
}

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

/* 
	Convert an EditScript chain to a SeqAlign of type DenseSeg.
	Used for a non-simple interval (i.e., one without subs. or 
	deletions).  

	The first SeqIdPtr in the chain of subject_id and query_id is duplicated for
	the SeqAlign.
*/

SeqAlignPtr LIBCALL
GapXEditBlockToSeqAlign(GapXEditBlockPtr edit_block, SeqIdPtr subject_id, SeqIdPtr query_id)

{
	Boolean reverse, translate1, translate2;
	GapXEditScriptPtr curr, esp;
	Int2 frame1, frame2, numseg;
	Int4 begin1, begin2, index, start1, start2, length1, length2;
	Int4 original_length1, original_length2;
	Int4Ptr length, start;
	DenseSegPtr dsp;
	SeqAlignPtr sap;
	SeqIntPtr seq_int1;
	SeqLocPtr slp, slp1, slp2;
	StdSegPtr sseg, sseg_head, sseg_old;
	Uint1 strand1, strand2;
	Uint1Ptr strands;

	reverse = edit_block->reverse;	

	numseg=0;
	start1 = edit_block->start1;
	start2 = edit_block->start2;
	length1 = edit_block->length1;
	length2 = edit_block->length2;
	original_length1 = edit_block->original_length1;
	original_length2 = edit_block->original_length2;
	translate1 = edit_block->translate1;
	translate2 = edit_block->translate2;
	frame1 = edit_block->frame1;
	frame2 = edit_block->frame2;

	if (frame1 > 0)
		strand1 = Seq_strand_plus; 
	else if (frame1 < 0)
		strand1 = Seq_strand_minus; 
	else
		strand1 = Seq_strand_unknown; 
		
	if (frame2 > 0)
		strand2 = Seq_strand_plus; 
	else if (frame2 < 0)
		strand2 = Seq_strand_minus; 
	else
		strand2 = Seq_strand_unknown; 

	esp = edit_block->esp;
	for (curr=esp; curr; curr=curr->next)
	{
		numseg++;
	}

	start = MemNew((2*numseg+1)*sizeof(Int4));
	length = MemNew((numseg+1)*sizeof(Int4));
	strands = MemNew((2*numseg+1)*sizeof(Uint1));

	index=0;
	for (curr=esp; curr; curr=curr->next)
	{
		switch(curr->op_type)
		{
			case GAPALIGN_SUB:
			case GAPALIGN_DECLINE:
				if (strand1 != Seq_strand_minus)
				{
				    if(translate1 == FALSE)
					begin1 = get_current_pos(&start1, curr->num);
				    else
					begin1 = frame1 - 1 + CODON_LENGTH*get_current_pos(&start1, curr->num);
				}
				else
				{
				    if(translate1 == FALSE)
					begin1 = length1 - get_current_pos(&start1, curr->num) - curr->num;
				    else
					begin1 = original_length1 - CODON_LENGTH*(get_current_pos(&start1, curr->num)+curr->num) + frame1 + 1;
				}
					
				if (strand2 != Seq_strand_minus)
				{
				    if(translate2 == FALSE)
					begin2 = get_current_pos(&start2, curr->num);
				    else
					begin2 = frame2 - 1 + CODON_LENGTH*get_current_pos(&start2, curr->num);
				}
				else
				{
				    if(translate2 == FALSE)
					begin2 = length2 - get_current_pos(&start2, curr->num) - curr->num;
				    else
					begin2 = original_length2 - CODON_LENGTH*(get_current_pos(&start2, curr->num)+curr->num) + frame2 + 1;
				}

				if (reverse)
				{
					strands[2*index] = strand2;
					strands[2*index+1] = strand1;
					start[2*index] = begin2;
					start[2*index+1] = begin1;
				}
				else
				{
					strands[2*index] = strand1;
					strands[2*index+1] = strand2;
					start[2*index] = begin1;
					start[2*index+1] = begin2;
				}

				break;

			case GAPALIGN_DEL:
				begin1 = -1;
				if (strand2 != Seq_strand_minus)
				{
				    if(translate2 == FALSE)
					begin2 = get_current_pos(&start2, curr->num);
				    else
					begin2 = frame2 - 1 + CODON_LENGTH*get_current_pos(&start2, curr->num);
				}
				else
				{
				    if(translate2 == FALSE)
					begin2 = length2 - get_current_pos(&start2, curr->num) - curr->num;
				    else
					begin2 = original_length2 - CODON_LENGTH*(get_current_pos(&start2, curr->num)+curr->num) + frame2 + 1;
				}

				if (reverse)
				{
					strands[2*index] = strand2;
					if (index > 0)
						strands[2*index+1] = strands[2*(index-1)+1];
					else
						strands[2*index+1] = Seq_strand_unknown;
					start[2*index] = begin2;
					start[2*index+1] = begin1;
				}
				else
				{
					if (index > 0)
						strands[2*index] = strands[2*(index-1)];
					else
						strands[2*index] = Seq_strand_unknown;
					strands[2*index+1] = strand2;
					start[2*index] = begin1;
					start[2*index+1] = begin2;
				}

				break;

			case GAPALIGN_INS:
				if (strand1 != Seq_strand_minus)
				{
				    if(translate1 == FALSE)
					begin1 = get_current_pos(&start1, curr->num);
				    else
					begin1 = frame1 - 1 + CODON_LENGTH*get_current_pos(&start1, curr->num);
				}
				else
				{
				    if(translate1 == FALSE)
					begin1 = length1 - get_current_pos(&start1, curr->num) - curr->num;
				    else
					begin1 = original_length1 - CODON_LENGTH*(get_current_pos(&start1, curr->num)+curr->num) + frame1 + 1;
				}
				begin2 = -1;
				if (reverse)
				{
					if (index > 0)
						strands[2*index] = strands[2*(index-1)];
					else
						strands[2*index] = Seq_strand_unknown;
					strands[2*index+1] = strand1;
					start[2*index] = begin2;
					start[2*index+1] = begin1;
				}
				else
				{
					strands[2*index] = strand1;
					if (index > 0)
						strands[2*index+1] = strands[2*(index-1)+1];
					else
						strands[2*index+1] = Seq_strand_unknown;
					start[2*index] = begin1;
					start[2*index+1] = begin2;
				}

				break;
			default:
				break;
		}
		length[index] = curr->num;
		index++;
	}

	sap = SeqAlignNew();

	sap->dim =2; /**only two dimention alignment**/

	/**make the Denseg Object for SeqAlign**/
	if (translate1 == FALSE && translate2 == FALSE)
	{
		sap->segtype = 2; /**use denseg to store the alignment**/
		sap->type = 3; /**partial for gapped translating search. */
		dsp = DenseSegNew();
		dsp->dim = 2;
		dsp->numseg = numseg;
		if (reverse)
		{
			dsp->ids = SeqIdDup(subject_id);
			dsp->ids->next = SeqIdDup(query_id);
		}
		else
		{
			dsp->ids = SeqIdDup(query_id);
			dsp->ids->next = SeqIdDup(subject_id);
		}
		dsp->starts = start;
		dsp->strands = strands;
		dsp->lens = length;
		sap->segs = dsp;
		sap->next = NULL;
	}
	else
	{
		sap->type =3; /**partial for gapped translating search. */
		sap->segtype =3; /**use denseg to store the alignment**/
		sseg_head = NULL;
		sseg_old = NULL;
		for (index=0; index<numseg; index++)
		{
			sseg = StdSegNew();
			sseg->dim = 2;
			if (sseg_head == NULL)
			{
				sseg_head = sseg;
			}
				if (reverse)
				{
					sseg->ids = SeqIdDup(subject_id);
					sseg->ids->next = SeqIdDup(query_id);
				}
				else
				{
					sseg->ids = SeqIdDup(query_id);
					sseg->ids->next = SeqIdDup(subject_id);
				}
			slp1 = NULL;
			if (start[2*index] != -1)
			{
				seq_int1 = SeqIntNew();
				seq_int1->from = start[2*index];
				if (translate1)
					seq_int1->to = start[2*index] + CODON_LENGTH*length[index] - 1;
				else
					seq_int1->to = start[2*index] + length[index] - 1;
				seq_int1->strand = strands[2*index];
				seq_int1->id = SeqIdDup(query_id);
				ValNodeAddPointer(&slp1, SEQLOC_INT, seq_int1);
			}
			else
			{
				ValNodeAddPointer(&slp1, SEQLOC_EMPTY, SeqIdDup(query_id));
			}
			slp2 = NULL;
			if (start[2*index+1] != -1)
			{
				seq_int1 = SeqIntNew();
				seq_int1->from = start[2*index+1];
				if (translate2)
					seq_int1->to = start[2*index+1] + CODON_LENGTH*length[index] - 1;
				else
					seq_int1->to = start[2*index+1] + length[index] - 1;
				seq_int1->strand = strands[2*index+1];
				seq_int1->id = SeqIdDup(subject_id);
				ValNodeAddPointer(&slp2, SEQLOC_INT, seq_int1);
			}
			else
			{
				ValNodeAddPointer(&slp2, SEQLOC_EMPTY, SeqIdDup(subject_id));
			}

			if (reverse)
			{
				slp = slp2;
				slp2->next = slp1;
			}
			else
			{
				slp = slp1;
				slp1->next = slp2;
			}
			sseg->loc = slp;

			if (sseg_old)
				sseg_old->next = sseg;
			sseg_old = sseg;
		}
		sap->segs = sseg_head;
		sap->next = NULL;
	}

	return sap;
}

/*
	SimpleIntervalToGapXEditBlock(Int4 start1, Int4 start2, Int4 length)

	Int4 start1, start2: offset of this interval in sequence 1 and 2.
	Int4 length: length of the interval.

	May be used to produce a gapXEditBlock when there are no subs. or deletions
	in the interval (e.g., ungapped BLAST HSP).
*/

GapXEditBlockPtr LIBCALL
SimpleIntervalToGapXEditBlock (Int4 start1, Int4 start2, Int4 length)

{
	GapXEditBlockPtr edit_block;
  	GapXEditScriptPtr e_script;

	edit_block = GapXEditBlockNew(start1, start2);

	edit_block->esp = e_script = GapXEditScriptNew(NULL);

	e_script->op_type = GAPALIGN_SUB;
	e_script->num = length;

	return edit_block;
}

