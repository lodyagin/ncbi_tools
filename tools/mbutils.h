/* $Id: mbutils.h,v 6.8 2000/03/29 21:58:01 dondosha Exp $
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
* File Name:  $RCSfile: mbutils.h,v $
*
* Author:  Webb Miller and Co.
* Adopted for NCBI standard libraries by Sergey Shavirin
*
* Initial Creation Date: 10/27/1999
*
* $Revision: 6.8 $
*
* File Description:
*         Main file for Mega Blast program
*
* $Log: mbutils.h,v $
* Revision 6.8  2000/03/29 21:58:01  dondosha
* Added prototypes for edit_script_new and edit_script_append
*
* Revision 6.7  1999/11/26 19:42:21  shavirin
* Fixed warnings of windows compiler.
*
* Revision 6.6  1999/11/24 20:38:10  shavirin
* Added possibility to produce Traditional Blast Output.
*
* Revision 6.5  1999/11/03 21:27:55  shavirin
* Added option to print output in real time.
*
* Revision 6.4  1999/11/03 19:55:11  shavirin
* All printing functions distinguished from search functions.
*
* Revision 6.3  1999/10/29 22:02:34  kans
* added prototypes for es_rep_len and es_indel_len required by Mac compiler
*
* Revision 6.2  1999/10/29 18:33:57  shavirin
* Removed usage of function snprintf().
*
* Revision 6.1  1999/10/29 15:39:30  shavirin
* Initial revision.
*
*
* Initial revision.
*
* ==========================================================================
*/

#ifndef _MBUTILS_H_
#define _MBUTILS_H_

#include <ncbi.h>
#include <readdb.h>

/* ----- From original file db.h ------- */

#ifdef SMALL
#define BLOCKLIM  (10)
#define SEQLIM    (3)
#define HDRSAVE   (50)
#else
#define BLOCKLIM  (5000000)
#define SEQLIM    (50000)
#define HDRSAVE   (50)
#endif
#define BLOCKSIZE (BLOCKLIM*2+1)
#define SEQSIZE   (SEQLIM*2+1)

/* This query database stores the concatenation of a number of query
 * sequences, seperated by 'X's, along with the file position, byte
 * offset, and saved header for each sequence.  The reverse-compliment of
 * each sequence is also stored.
 */

typedef struct qdb {
    Int4 fpos[SEQSIZE];
    Char hdr[SEQSIZE][HDRSAVE];
    Int4 beg[SEQSIZE];
    Uchar sequence[BLOCKSIZE];
    Uchar *p; /* next free position in sequence */
    Int4 nseq;
} qdb_t;

/* -----  End from original file db.h ----- */

/* ----- From original file edit.h ------- */

typedef Uint4 edit_op_t; /* 32 bits */
typedef struct {
    edit_op_t *op;                  /* array of edit operations */
    Uint4 size, num;         /* size of allocation, number in use */
    edit_op_t last;                 /* most recent operation added */
} edit_script_t;
edit_script_t *edit_script_free(edit_script_t *es);
edit_script_t *edit_script_new();
edit_script_t *edit_script_append(edit_script_t *es, edit_script_t *et);

/* ----  End from original file edit.h ---- */

/* ----- From original file mbalign.h ------- */

#define GAL_CLOSE_ENOUGH 10 /* XXX - arbitrary */

#define GAL_CLOSE(a,b) \
        (abs(((a)->b1 - (a)->b2) - ((b)->b1 - (b)->b2)) < GAL_CLOSE_ENOUGH)

#define GAL_CONTAINS(outer,inner) \
        ((outer)->b1 <= (inner)->b1 && (outer)->b2 <= (inner)->b2 && \
         (outer)->e1 >= (inner)->e1 && (outer)->e2 >= (inner)->e2)

typedef struct gapped_align {
    struct gapped_align *next;
    Int4 b1, b2, e1, e2;  /* bounding box */
    Int4 diag, length;
    Int4 query_id;       /* Sequence number of query seq in qdb storage */
    Int4 subject_id;     /* Sequence number of subject if rdfp database */
    Int4 query_strand; /* First part of qdb == Seq_strand_plus */
    Int4 max, prune;
    edit_script_t *S;
} gal_t;

Int4 mb_extend_gap_align_script(gal_t *mq, Int4 X, Int4 M, Int4 R, Int4 G, 
                                Int4 E, const UcharPtr s1, Int4 l1, 
                                const UcharPtr s2, Int4 l2);

Int4 mb_extend_gap_align(gal_t *mq, Int4 X, Int4 M, Int4 R, Int4 G, Int4 E, 
                         const UcharPtr s1, Int4 l1, 
                         const UcharPtr s2, Int4 l2);


/* ----  End from original file mbalign.h ---- */

/* ----- From original file mb.c ------- */
enum {
    ESTACK_SIZE = 5000
};

struct stack {
    Int4 diag, level, length;
};

#define NACHARS 128
typedef Int4 ss_t[NACHARS][NACHARS];

typedef struct mb_table {
    UcharPtr seq;
    Int4 slen;
    
    Int4 W;		/* width */
    Int4 K;		/* blastz threshold */
    Int4 X;		/* xdrop condition */
    Int4 k;		/* gap_align_cutoff */
    Int4 x;		/* xdrop_gapfree */
    Int4 lpm;		/* length of perfect match */
    Int4 surfeit;	/* of 12mer entries */
    ss_t ss;
    
    Int4    HASH_SIZE;
    Int4Ptr hashtab;
    Int4Ptr next_pos;
    Uint4   mask;
    
    Int4 stack_index;
    struct stack estack[ESTACK_SIZE];
    
    Int4 numMSPs;
    ReadDBFILEPtr rdfp;
    gal_t *gal_list;
} mbt_t;

typedef struct align_params {
    Int4 prog_flag;
    Int4 gap_align_cutoff;
    Int4 gap_extension_penalty;
    Int4 gap_open_penalty;
    Int4 idlen;
    Int4 len_perf_match;
    Int4 match_score;
    Int4 transition_score;
    Int4 transversion_score;
    Int4 max_12mer_ents;
    Int4 mismatch_penalty_mb;
    Int4 normalize_revcomp_flag;
    Int4 xdrop_ext;
    Int4 xdrop_gapfree;
    Int4 output_al;
    Int4 revcomp_flag;
    Int4 real_time_print;
    FILE *outfp;
} align_params_t;
/* ----  End from original file mb.c ---- */

typedef struct _MBHitlist {
    Int4 hitcount;
    Int4 allocated;
    gal_t **galpp; /* Array of alignment structures */
} MBHitlist, PNTR MBHitlistPtr;

MBHitlistPtr MBEngineCore(mbt_t *mbt, qdb_t *qdb, align_params_t *P);
void MBPrintHitList(mbt_t *mbt, MBHitlistPtr mbhlp, qdb_t *qdb, 
                    align_params_t *mbap);
void MBHitlistFree(MBHitlistPtr mbhlp);

ValNodePtr MBSeqAlignFromHitList(mbt_t *mbt, qdb_t *qdb,
                                 MBHitlistPtr mbhlp, ValNodePtr *bsp_vnp,
                                 Int4Ptr num_seqs, align_params_t *mbap);

Int4 est_search(mbt_t *);

Int4 qdb_fill_block(FILE *fp, qdb_t *qdb, Int4 do_revcomp);
Int4 qdb_total_length(const qdb_t *qdb);

mbt_t *mb_init(CharPtr blast_database);
void mb_start(mbt_t *mbt, UcharPtr dna, Int4 ld, Int4 w, Int4 p);
void mb_end(mbt_t *mbt);

void free_gal_list(gal_t *list);
SeqAlignPtr MBCreateSeqAlign(gal_t * galp, SeqIdPtr subject_id, 
                             SeqIdPtr query_id);
void mb_print_sequence(UcharPtr seq, Int4 len, CharPtr);

Int4 es_rep_len(edit_script_t *S, Int4Ptr n, UcharPtr p, 
                UcharPtr q, Int4Ptr match);
Int4 es_indel_len(edit_script_t *S, Int4Ptr n, Int4Ptr i, Int4Ptr j);

enum {
    MB_DO_EST    = 0,
    MB_DO_GREEDY = 1
};

enum {
    MB_ALIGN_NONE = 0,
    MB_ALIGN_SUM  = 1,
    MB_ALIGN_LAV  = 2
};

enum {
    MATCH_SCORE_SCALE = 10
};

#endif /* _MBUTILS_H_ */

