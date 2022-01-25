#ifndef _MBALIGN_H_
#define _MBALIGN_H_
#ifdef __cplusplus
extern "C" {
#endif

#include <ncbi.h>
#include <mbutils.h>

enum {
    EDIT_OP_MASK = 0x3,
    EDIT_OP_ERR  = 0x0,
    EDIT_OP_INS  = 0x1,
    EDIT_OP_DEL  = 0x2,
    EDIT_OP_REP  = 0x3
};

enum {         /* half of the (fixed) match score */
    ERROR_FRACTION=2,  /* 1/this */
    MAX_SPACE=1000000,
    sC = 0, sI = 1, sD = 2, LARGE=100000000
};

#define ICEIL(x,y) ((((x)-1)/(y))+1)

/* ----- pool allocator ----- */
typedef struct three {
    Int4 I, C, D;
} three_val;

typedef struct space_struct {
    three_val *space_array;
    Int4 used, size;
} space_type, *space_ptr;

#define EDIT_VAL(op) (op >> 2)

#define EDIT_OPC(op) (op & EDIT_OP_MASK)

space_ptr new_space(Int4 MAX_D);
void free_space(space_ptr sp);
three_val *get_space(space_ptr S, Int4 amount);
Int4 get_last(Int4 **flast_d, Int4 d, Int4 diag, Int4 *row1);

typedef struct greedy_align_mem {
   Int4 **flast_d;
   Int4 *max_row_free;
} GreedyAlignMem, PNTR GreedyAlignMemPtr;

Int4 
greedy_gapped_align PROTO((const UcharPtr s1, Int4 len1,
			     const UcharPtr s2, Int4 len2,
			     Boolean reverse, Int4 xdrop_threshold, 
			     Int4 match_cost, Int4 mismatch_cost,
			     Int4 *e1, Int4 *e2, GreedyAlignMemPtr abmp, 
			     edit_script_t *S));

GreedyAlignMemPtr 
GreedyAlignMemFree PROTO((GreedyAlignMemPtr abmp));

#ifdef __cplusplus
}
#endif
#endif /* _MBALIGN_H_ */
