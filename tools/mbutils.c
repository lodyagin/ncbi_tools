/* $Id: mbutils.c,v 6.7 1999/11/26 19:42:21 shavirin Exp $
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
* File Name:  $RCSfile: mbutils.c,v $
*
* Author:  Webb Miller and Co.
* Adopted for NCBI standard libraries by Sergey Shavirin
*
* Initial Creation Date: 10/27/1999
*
* $Revision: 6.7 $
*
* File Description:
*        All core functions for Mega Blast program
*
* $Log: mbutils.c,v $
* Revision 6.7  1999/11/26 19:42:21  shavirin
* Fixed warnings of windows compiler.
*
* Revision 6.6  1999/11/24 20:38:10  shavirin
* Added possibility to produce Traditional Blast Output.
*
* Revision 6.5  1999/11/08 19:10:33  shavirin
* Fixed minor SGI warning.
*
* Revision 6.4  1999/11/03 21:27:54  shavirin
* Added option to print output in real time.
*
* Revision 6.3  1999/11/03 19:55:11  shavirin
* All printing functions distinguished from search functions.
*
* Revision 6.2  1999/10/29 18:33:56  shavirin
* Removed usage of function snprintf().
*
* Revision 6.1  1999/10/29 15:38:16  shavirin
* Initial revision.
*
*
* Initial revision.
*
* ==========================================================================
*/

#include <mbutils.h>
#include <blastkar.h>

/* ----- From original file encoding.h ------- */
#define _ (-1)
static const Int1 fasta_encoding[256] = {
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, 0, _, 1, _, _, _, 2, _, _, _, _, _, _, _, _,
 _, _, _, _, 3, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
 _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
};

static const Uchar fasta_decoding[] = "ACGT";

static const Uchar dna_complement[256] = 
"                                                                "
" TVGH  CD  M KN   YSA BWXR       tvgh  cd  m kn   ysa bwxr      "
"                                                                "
"                                                                ";
/* ................................................................ */
/* @ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~. */
/* ................................................................ */
/* ................................................................ */

/* ----  End from original file encoding.h ---- */

/* -------- Start of functions from various original files -------- */

/* -------- From original file prnt.c ------------- */

/* XXX */
static Int4 offset1;
static Int4 offset2;

static Int4 align_match_percent(Int4 run, Int4 match)
{
    return (run > 0) ? ((100*match + run/2)/run) : 0; /* round */
}

static void prefn_lav(Int4 b1, Int4 b2, Int4 e1, Int4 e2, Int4 s)
{
    printf("a {\n  s %d\n  b %d %d\n  e %d %d\n", s, b1, b2, e1, e2);
}
static void perfn_lav(Int4 b1, Int4 b2, Int4 e1, Int4 e2, Int4 p)
{
    printf("  l %d %d %d %d %d\n", b1, b2, e1, e2, p);
}

/*ARGSUSED*/
static void postfn_lav(Int4 t, Int4 m, Int4 n, Int4 l)
{
    printf("}\n");
}
static void prefn_sum(Int4 b1, Int4 b2, Int4 e1, Int4 e2, Int4 s) { /* */ }
static void perfn_sum(Int4 b1, Int4 b2, Int4 e1, Int4 e2, Int4 p) { /* */ }
static void postfn_sum(Int4 t, Int4 m, Int4 n, Int4 l)
{
    /* match, mismatch, #Gaps, total gap length */
    printf("%d %d %d %d\n", m, t-m, n, l);
}

static void print_align_lav_fn(
    const UcharPtr seq1, const UcharPtr seq2,
    Int4 beg1, Int4 end1,
    Int4 beg2, Int4 end2,
    edit_script_t *S,
    Int4 score,
    void (*prefn)(Int4 b1, Int4 b2, Int4 e1, Int4 e2, Int4 s),
    void (*perfn)(Int4 b1, Int4 b2, Int4 e1, Int4 e2, Int4 p),
    void (*postfn)(Int4 t, Int4 m, Int4 n, Int4 l)
)
{
    Int4 M, N, i, j, o, nG = 0, lG = 0;
    Int4 total_run = 0, total_match = 0;
    
    beg1 += offset1;
    end1 += offset1;
    beg2 += offset2;
    end2 += offset2;
    
    M = end1 - beg1 + 1;
    N = end2 - beg2 + 1;
    
    prefn(beg1, beg2, end1, end2, score);
    for (o = i = j = 0; i < M || j < N; ) {
        Int4 start_i = i;
        Int4 start_j = j; 
        Int4 match = 0;
        Int4 run;
        
        run = es_rep_len(S, &o, seq1+beg1+i-1, seq2+beg2+j-1, &match);
        i += run; j += run;
        total_match += match;
        total_run += run;
        
        perfn(beg1+start_i, beg2+start_j, beg1+i-1, beg2+j-1, 
              align_match_percent(run, match));
        
        if (i < M || j < N) {
            nG++;
            lG += es_indel_len(S, &o, &i, &j);
        }
        
    }
    postfn(total_run, total_match, nG, lG);
}
static void print_align_lav(Int4 score, const UcharPtr seq1, const UcharPtr seq2,
                     Int4 beg1, Int4 end1, Int4 beg2, Int4 end2, edit_script_t *S)
{

    print_align_lav_fn(seq1, seq2, beg1, end1, beg2, end2, S, score,
                       prefn_lav, perfn_lav, postfn_lav);
}
static void print_align_summary(Int4 score, const UcharPtr seq1, const UcharPtr seq2,
                         Int4 beg1, Int4 end1, Int4 beg2, Int4 end2, 
                         edit_script_t *S)
{
    
    print_align_lav_fn(seq1, seq2, beg1, end1, beg2, end2, S, score,
                       prefn_sum, perfn_sum, postfn_sum);
}

/* -------- From original file unpack.c ------------- */

static UcharPtr readdb_unpack(const UcharPtr seq, Int4 len)
{
    UcharPtr s, t;
    Int4 i, j;
    const Uchar *const code = fasta_decoding;
    
    t = s = MemNew(len + 2 + 1);
    for (i = 0; i < len / 4; i++) {
	*t++ = code[READDB_UNPACK_BASE_1(seq[i])];
	*t++ = code[READDB_UNPACK_BASE_2(seq[i])];
	*t++ = code[READDB_UNPACK_BASE_3(seq[i])];
	*t++ = code[READDB_UNPACK_BASE_4(seq[i])];
    }
    if (len % 4 != 0) {
	j = seq[i] % 4;
	i = seq[i];
	/*debugf("%d %d\n",j, i); */
	for (; j > 0; j--) {
	    *t++ = code[(i & 192) >> 6];
	    i = (i << 2);
	}
    }
    *t = '\0';
    return s;
}


/* -------- From original file print.c ------------- */

static CharPtr skipword(CharPtr s)
{
    s += strspn(s, " \t");
    s += strcspn(s, " \t");
    return s;
}

static Int4 nwordlen(Int4 n, CharPtr s)
{
    CharPtr t = s;
    while (n-- > 0) {
        t = skipword(t);
    }
    return t-s;
}

/* External */
void mb_print_sequence(UcharPtr seq, Int4 len, CharPtr defline)
{
    Int4 i, col, remain;
    Uchar ch;
    const Uchar *const code = fasta_decoding;
    
    printf(">%s\n", defline);
    for (i = col = 0; i < len / 4; i++) {
	printf("%c%c%c%c",
               code[READDB_UNPACK_BASE_1(seq[i])],
               code[READDB_UNPACK_BASE_2(seq[i])],
               code[READDB_UNPACK_BASE_3(seq[i])],
               code[READDB_UNPACK_BASE_4(seq[i])]);
	col += 4;
	if (col >= 50) {
	    printf("\n");
	    col = 0;
	}
    }
    remain = len % 4;

    if (col == 0 && remain == 0)
	return;
    
    ch = seq[i];
    for (i = 0; i < remain; ++i) {
	printf("%c", code[READDB_UNPACK_BASE_1(ch)]);
	ch <<= 2;
    }
    printf("\n");
}

static void print_mb_result(Int4 idlen, const CharPtr def, const CharPtr hdr,
                            Int4 b1, Int4 b2, Int4 e1, Int4 e2, Int4 diag, Int4 nl)
{
    Int4 nd = idlen+1, nh = idlen+1; /* +1 for the ">" */
    
    if (idlen < 0) {
	nd = nwordlen(-idlen, def);
	nh = nwordlen(-idlen, hdr);
        /* nh used to get -idlen+1, but Lukas said he likes this better. */
    }
    printf("'%.*s'=='%.*s' (%d %d %d %d) %d%c",
	nd, def, nh, hdr, b1, b2, e1, e2, diag, nl);
}

/* -------- From original file db.c ------------- */

enum { Ok, Full, Eof };

static void fmt_id(CharPtr h, Int4 c, CharPtr p)
{
    sprintf(h, "%c%.*s", c, HDRSAVE-2, p);
}

static void save_id(qdb_t *qdb, Int4 i, Int4 c, CharPtr p)
{
    fmt_id(qdb->hdr[i], c, p);
}

/* read_id -- read the id line from fp, and save some of it. */
static Int4 read_id(FILE *fp, qdb_t *qdb, Int4 i)
{
    Char buf[1024]; /* XXX - arbitrary limit */
    Int4 n;
    
    if (fgets(buf, sizeof(buf), fp) == 0)
        return Eof;

    buf[strlen(buf)-1] = 0;
    n = strspn(buf+1, " \t");
    save_id(qdb, i, '+', buf+1+n);

    return Ok;
}

static UcharPtr qdb_seq(qdb_t *qdb, Int4 i)
{
    ASSERT(0 <= i);
    ASSERT(i < qdb->nseq);
    
    return qdb->sequence + qdb->beg[i];
}
static Int4 qdb_length(const qdb_t *qdb, Int4 i)
{
    ASSERT(0 <= i);
    ASSERT(i < qdb->nseq);
    
    return qdb->beg[i+1] - qdb->beg[i];
}
static Int4 qdb_is_revcomp(qdb_t *qdb, Int4 i)
{
    ASSERT(0 <= i);
    ASSERT(i < qdb->nseq);

    return qdb->hdr[i][0] == '-';  /* XXX -- kludge */
}

/* copy sequence data from fp into qdb->p[].
 * return pointer to end, or 0 if overflowed.
 * rely on caller to update qdb->p for next time.
 */
static UcharPtr read_seq(FILE *fp, qdb_t *qdb)
{
    UcharPtr p = qdb->p;
    
#if 0
    *p++ = 'X';
#endif
    while (p < qdb->sequence+BLOCKLIM) {
        Int4 c = getc(fp);
	/* test most common case first.  
           the manpage says EOF is an acceptable argument for is*(). */
        if (isalpha(c)) {
            *p++ = c & 0xFF;
        } else if (isspace(c)) {
            ; /* skip it */
        } else if (ispunct(c)) {
            ungetc(c, fp);
            return p;
        } else if (c == EOF) {
            return p;
        } else {
            /* can't happen, hopefully */
            ErrPostEx(SEV_ERROR, 0,0 , "bad character '%c' in read_seq", c);
            exit(1);
        }
    }
    return 0;
}

/* read_one --- Read one sequence from fp.
 * Return Ok for success, Full if the qdb table is full,
 * and Eof if we ran out of input.  
 * The file must be seekable.
 */
static Int4 read_one(FILE *fp, qdb_t *qdb, Int4Ptr len)
{
    UcharPtr q;
    Int4 i = qdb->nseq;
    
    if (i >= SEQLIM)
        return Full;
    
    if (feof(fp))
        return Eof;
    
    qdb->fpos[i] = ftell(fp);
    qdb->beg[i] = qdb->p - qdb->sequence;
    qdb->hdr[i][0] = 0;
    
    if (i >= SEQLIM)
        return Full;
    
    read_id(fp, qdb, i);
    if ((q = read_seq(fp, qdb)) != 0) {
	/* got one */
        *q++ = 'X';
	*q = 0;
	*len = q - qdb->p;
	qdb->p = q;
        qdb->beg[i+1] = qdb->p - qdb->sequence;
        return Ok;
    } else {
	/* rewind for next time */
	fseek(fp, qdb->fpos[i], SEEK_SET);
	*qdb->p = 0;
        
        return Full;
    }
}

/* make_revcomp_seq -- given a query sequence in qdb, construct
 * its reverse complement and append it to the original.
 * assume that there is enough space.
 */
static void make_revcomp_seq(qdb_t *qdb)
{
    UcharPtr p = qdb->p;
    UcharPtr rp = p;

    --p; /* p's nul and preceeding 'X' go at end. */

    while (--p >= qdb->sequence)
        *rp++ = dna_complement[*p];
    *rp++ = 'X';
    *rp = 0;
}

/* make_revcomp_tab -- add new header information to a qdb
 * to reflect the recently appended reverse complement.
 * assume that there is enough space.
 */
static void make_revcomp_tab(qdb_t *qdb)
{
    Int4 i, j;
    i = qdb->nseq;

    for (j = 1; j <= i; ++j) {
        qdb->beg[i+j] = qdb->beg[i]+(qdb->beg[i]-qdb->beg[i-j]);
        save_id(qdb, i+j-1, '-', qdb->hdr[i-j]+1);
    }

    qdb->fpos[2*i] = -1;
    qdb->hdr[2*i][0] = 0;
}

/* make_revcomp -- modify qdb in place by appending the query
 * sequence's reverse complement, and updating the headers to match.
 * return the new size of the table.
 */
static Int4 make_revcomp(qdb_t *qdb)
{
    make_revcomp_seq(qdb);
    make_revcomp_tab(qdb);
    qdb->nseq *= 2;
    
    return qdb->nseq;
}

/* read_fasta -- read as many as possible sequences from
 * a FASTA format file, and store them in the qdb.
 */
static Int4 read_fasta(FILE *fp, qdb_t *qdb)
{
    Int4 len = 0;
    qdb->nseq = 0;
    qdb->p = qdb->sequence;
    while (read_one(fp, qdb, &len) == Ok)
        ++qdb->nseq;

    return qdb->nseq;
}

/* qdb_fill_block -- Read a number of sequences from a FASTA
 * format file and store them and their reverse-complements
 * in a query database structure.  The file must be seekable.
 */
Int4 qdb_fill_block(FILE *fp, qdb_t *qdb, Int4 do_revcomp)
{
    Int4 i = read_fasta(fp, qdb);
    if (do_revcomp)
        i = make_revcomp(qdb);
    return i;
}
Int4 qdb_total_length(const qdb_t *qdb)
{
    return qdb->beg[qdb->nseq];
}

/* -------- End of functions from various original files -------- */


/* normalize -- transform reverse compliment coord to normal coord */
static Int4 normalize(Int4 want, Int4 state, Int4 len, Int4 x)
{
    return (want && state) ? len-x : x; /* len+1? */
}

void free_gal_list(gal_t *list)
{
    while (list) {
        gal_t *p = list->next;	
        edit_script_free(list->S);
        MemFree(list);
        list = p;
    }
    return;
}

static void add_gal(mbt_t *mbt, Int4 level, Int4 diag, Int4 length)
{
    gal_t *mp;
    
    if (length <= mbt->lpm*4)
	return;
    
    mp = MemNew(sizeof(*mp)); 
    mp->b1 = level - length + 1;
    mp->e1 = level;
    mp->b2 = mp->b1 - diag;
    mp->e2 = mp->e1 - diag;
    mp->length = length;
    
    mp->next = mbt->gal_list;
    mbt->gal_list = mp;
    
    mbt->numMSPs++;
}

static Int4 extend_hit(mbt_t *mbt, Int4 i, Int4 p, Int4 prog_flag)
{
    Int4 j, len = mbt->W*4;
    const Int4 EST = prog_flag == MB_DO_EST;

    for (j = 0; j < mbt->stack_index;) {
	if (mbt->estack[j].diag == i - p) {
	    if (mbt->estack[j].level < i - 4) {
		if (!EST)
		    add_gal(mbt, mbt->estack[j].level, mbt->estack[j].diag, mbt->estack[j].length);
		mbt->estack[j].length = len;
		mbt->estack[j].level = i;
		return 0;
	    }
	    mbt->estack[j].length += 4;
	    if (EST)
		if (mbt->estack[j].length > mbt->lpm*4)
		    return 1;
	    mbt->estack[j].level = i;
	    return 0;
	} else if (mbt->estack[j].level >= i - 4)
	    j++;
	else {
	    mbt->stack_index--;
	    if (!EST)
		add_gal(mbt, mbt->estack[j].level, mbt->estack[j].diag, mbt->estack[j].length);
	    mbt->estack[j] = mbt->estack[mbt->stack_index];
	}
    }

    if (mbt->stack_index >= ESTACK_SIZE) {
	ErrPostEx(SEV_FATAL, 0, 0, "stack overflow in extend_hit");
        return -1;
    }
    
    mbt->estack[mbt->stack_index].diag = i - p;
    mbt->estack[mbt->stack_index].level = i;
    mbt->estack[mbt->stack_index].length = len;
    mbt->stack_index++;
    
    return 0;
}

static gal_t *search(mbt_t *mbt, UcharPtr s, Int4 len, 
                     Int4 prog_flag, Int4 N[])
{
    register Int4 h;
    register UcharPtr t;
    register Uint4 ecode;
    register Uint4 mask;
    Int4 i;
    const Int4 EST = prog_flag == MB_DO_EST;

    mbt->numMSPs = 0;
    mbt->gal_list = 0;
    t = s - 1;
    mbt->stack_index = 0;
    ecode = 0L;
    for (i = 0; i < (mbt->W - 1)*4; i += 4) {
	ecode = (ecode << 8) + *++t;
    }
    i += 4;
    mask = mbt->mask;
    while (i < len) {
	ecode = ((ecode & mask) << 8) + *++t;
	for (h = mbt->hashtab[ecode]; h; h = mbt->next_pos[h])
	    if (!EST)
		(void) extend_hit(mbt, i, h, prog_flag);
	    else if (extend_hit(mbt, i, h, prog_flag))
		return (gal_t *) 1; /* XXX */
	i += 4;
    }

    if (EST) {
	return 0;
    } else {
	for (i = 0; i < mbt->stack_index; i++)
	    add_gal(mbt, mbt->estack[i].level, mbt->estack[i].diag, mbt->estack[i].length);

	*N = mbt->numMSPs;
	return mbt->gal_list;
    }
}

static Int4 binary_search(Int4 n, Int4 A[], Int4 i)
{
    Int4 m, b, e;

    b = 0;
    e = i;
    while (b < e - 1) {
	m = (b + e) / 2;
	if (A[m] > n)
	    e = m;
	else
	    b = m;
    }
    return b;
}

static void readdb_get_stuff(ReadDBFILEPtr rdpt, Int4 dbidx,
                             UcharPtr PNTR s, Int4Ptr slen, 
                             CharPtr PNTR defline)
{
    SeqIdPtr sip = 0;
    *slen = readdb_get_sequence(rdpt, dbidx, s);
    readdb_get_descriptor(rdpt, dbidx, &sip, defline);
    sip = SeqIdSetFree(sip);
}

void MBPrintResult(mbt_t *mbt, gal_t *gal_head, qdb_t *qdb, align_params_t *P)
{
    gal_t *gp;
    Char def[HDRSAVE];
    Int4 tlen, q_len, cmpl, b1, e1, b2, e2;
    UcharPtr s_seq, t_seq;
    CharPtr defline;

    if(gal_head == NULL || qdb == NULL)
        return;
    
    for (gp = gal_head; gp; gp = gp->next) {
        
        if (gp->prune != 0)
            continue;

        readdb_get_stuff(mbt->rdfp, gp->subject_id, &t_seq, 
                         &tlen, &defline);
        
        fmt_id(def, '>', defline);
        s_seq = readdb_unpack(t_seq, tlen);

        q_len = qdb_length(qdb, gp->query_id);
        cmpl = (gp->query_strand == Seq_strand_minus);
        
        b1 = gp->b1;
        e1 = gp->e1;
        b2 = normalize(P->normalize_revcomp_flag, cmpl, q_len, gp->b2);
        e2 = normalize(P->normalize_revcomp_flag, cmpl, q_len, gp->e2);
        
        if (P->output_al == MB_ALIGN_LAV) putchar('#');
        
        print_mb_result(P->idlen, def, qdb->hdr[gp->query_id],
                        b1, b2, e1, e2, gp->diag, ' ');
        if (gp->S) {
            if (P->output_al == MB_ALIGN_LAV) {
                putchar('\n');
                print_align_lav(gp->diag, s_seq, 
                                qdb_seq(qdb,gp->query_id), 
                                gp->b1, gp->e1, gp->b2, gp->e2, gp->S);
            } else {
                putchar(' ');
                print_align_summary(gp->diag, s_seq, 
                                    qdb_seq(qdb,gp->query_id),
                                    gp->b1, gp->e1, gp->b2, gp->e2, gp->S);
            }
            
        } else {
            putchar('\n');
        }
        
        MemFree(s_seq);
        MemFree(defline);
    }
}
void MBPrintHitList(mbt_t *mbt, MBHitlistPtr mbhlp, qdb_t *qdb, 
                    align_params_t *mbap)
{
    Int4 i;
    
    for(i = 0; i <  mbhlp->hitcount; i++) {
        MBPrintResult(mbt, mbhlp->galpp[i], qdb, mbap);
    }    
}

/* Public functions belong to the object MBHitlistPtr */

#define MBHLP_ALLOC 1024
MBHitlistPtr MBHitlistNew(void)
{
    MBHitlistPtr mbhlp;
    
    mbhlp = MemNew(sizeof(MBHitlist));
    
    if(mbhlp == NULL)
        return NULL;
    
    mbhlp->allocated = MBHLP_ALLOC;
    mbhlp->galpp = MemNew(sizeof(gal_t *) * mbhlp->allocated);
    
    if(mbhlp->galpp == NULL)
        return NULL;
    
    mbhlp->hitcount = 0;
    
    return mbhlp;
}

static Boolean MBHitlistAddOne(MBHitlistPtr mbhlp, gal_t *gp)
{
    
    if(mbhlp->hitcount >= mbhlp->allocated) {
        mbhlp->allocated += MBHLP_ALLOC;
        mbhlp->galpp = Realloc(mbhlp->galpp, 
                               sizeof(gal_t *) * mbhlp->allocated);
        
        if(mbhlp->galpp == NULL)
            return FALSE;   
    }
    mbhlp->galpp[mbhlp->hitcount] = gp;
    mbhlp->hitcount++;
    
    return TRUE;
}

void MBHitlistFree(MBHitlistPtr mbhlp)
{
    Int4 i;
    
    for(i = 0; i <  mbhlp->hitcount; i++) {
        free_gal_list(mbhlp->galpp[i]);
    }
    
    MemFree(mbhlp->galpp);
    MemFree(mbhlp);
}

/* --------------------------------------------- */

static SeqEntryPtr MBReadQuerySeqEntry(qdb_t *qdb, Int4 id)
{
    Int4 length;
    CharPtr buffer, next_char;
    SeqEntryPtr sep;
    
    length = qdb_length(qdb, id) + sizeof(qdb->hdr[id]) + 2;
    buffer = MemNew(length);
    length = StringLen(qdb->hdr[id]);

    MemCpy(buffer + 1, qdb->hdr[id], length);
    *buffer = '>'; /* Replacing '+' or '-' */
    *(buffer + length + 1) = '\n'; /* Immediately after defline */
    MemCpy(buffer+length+2, qdb_seq(qdb, id), qdb_length(qdb, id));

    sep = FastaToSeqBuffEx(buffer, &next_char, TRUE, NULL, TRUE);
    
    MemFree(buffer);
    
    return sep;
}

static int LIBCALLBACK MBSortCallback(VoidPtr i, VoidPtr j)
{
    gal_t *gal1, *gal2;
    
    gal1 =  *((gal_t **) i);
    gal2 =  *((gal_t **) j);


    
    if (gal1->query_id > gal2->query_id)
        return (1);
    
    if (gal1->query_id < gal2->query_id)
        return (-1);

    if(gal1->query_id == gal2->query_id) {
        if (gal2->diag > gal1->diag)
            return 1;
        if (gal2->diag < gal1->diag)
            return -1;
    }
    
    return (0);
}
Boolean MBSortHitList(MBHitlistPtr mbhlp)
{
    HeapSort(mbhlp->galpp, mbhlp->hitcount, 
             sizeof(mbhlp->galpp[0]), MBSortCallback);
    
    return TRUE;
}

ScorePtr MBMakeScore(Int4 score, Nlm_FloatHi bit_score, Nlm_FloatHi e_value)
{
    ScorePtr scorep = NULL;
    
    MakeBlastScore(&scorep, "score", 0.0, score);
    MakeBlastScore(&scorep,"bit_score", bit_score, 0);
    MakeBlastScore(&scorep, "e_value", e_value, 0);
    
    return scorep;
}

static Int4 MBCalculateScore(gal_t *galp, align_params_t *mbap)
{
    Int4 score;

    if(mbap->gap_extension_penalty != 
       mbap->mismatch_penalty_mb + MATCH_SCORE_SCALE/2) {
        score = galp->diag;
    } else {
        score = ((galp->e1 - galp->b1) - galp->diag) * 10 - 
            galp->diag * mbap->mismatch_penalty_mb;
    }
    
    score = (1.37*score + 0.34)/0.69;
    
    if(score > 0)
        return score;
    
    return 0;
}

static Nlm_FloatHi MBCalculateEvalue(gal_t *galp, Int4 score)
{
    Nlm_FloatHi evalue;
    
    evalue = (galp->e1 - galp->b1)*(galp->e2 - galp->b2)*pow(2, (-1.0*score));
    
    return evalue;
}

ValNodePtr MBSeqAlignFromHitList(mbt_t *mbt, qdb_t *qdb, 
                                 MBHitlistPtr mbhlp, 
                                 ValNodePtr *sep_vnp, Int4Ptr num_seqs,
                                 align_params_t *mbap)
{
    SeqIdPtr query_sip, subject_sip;
    SeqAlignPtr sap, sap_head = NULL, sap_tmp;
    Int4 old_query_id, i;
    CharPtr defline;
    gal_t *galp, *galp_tmp;
    BioseqPtr query_bsp;
    ValNodePtr vnp_sap = NULL;
    SeqEntryPtr qsep;
    Int4 align_score;
    Nlm_FloatHi align_evalue;

    *num_seqs = 0;
    for(old_query_id = -1, i = 0; i <  mbhlp->hitcount; i++) {
        
        galp = mbhlp->galpp[i];
        
        if(old_query_id == -1) {
            old_query_id = galp->query_id;
        } else {
            if(old_query_id != galp->query_id) { /* New SeqAlign is coming */
                old_query_id = galp->query_id;
                sap_head = NULL;
                i--; /* We want to pass this entry anyway */
                continue;
            }
        }

        qsep  = MBReadQuerySeqEntry(qdb, galp->query_id);
        query_bsp = (BioseqPtr) qsep->data.ptrvalue;
        
        ValNodeAddPointer(sep_vnp, 0, qsep); 
        query_sip = query_bsp->id;
        
        readdb_get_descriptor(mbt->rdfp, galp->subject_id, 
                              &subject_sip, &defline);
        MemFree(defline);
        
        for(galp_tmp = galp; galp_tmp != NULL; galp_tmp = galp_tmp->next) {
            
            sap = MBCreateSeqAlign(galp_tmp, subject_sip, query_sip);

            align_score = MBCalculateScore(galp_tmp, mbap);
            align_evalue = MBCalculateEvalue(galp_tmp, align_score);
            
            sap->score = MBMakeScore(align_score, 
                                     (Nlm_FloatHi) align_score, align_evalue);
            
            if(sap_head == NULL) {
                sap_head = sap;    
                ValNodeAddPointer(&vnp_sap, 0, sap_head);
                (*num_seqs)++;
            } else {
                sap_tmp->next = sap;
            }
            sap_tmp = sap;
        }
        SeqIdSetFree(subject_sip);
    }
    return vnp_sap;
}

MBHitlistPtr MBEngineCore(mbt_t *mbt, qdb_t *qdb, align_params_t *mbap)
{
    Int4 N, dbidx, tlen;
    gal_t *gal_orig, *gp, *gal_new = NULL;
    gal_t *gal_tmp, *gal_prev;
    UcharPtr s_seq, t_seq;
    CharPtr defline;
    SeqAlignPtr sap;
    SeqIdPtr subject_sip, query_sip;
    MBHitlistPtr mbhlp = NULL;
    BioseqPtr bsp;

    if(mbt == NULL)
        return 0;

    if(!mbap->real_time_print)
        mbhlp = MBHitlistNew();
        
    /* Loop over all sequences in the blast database */
    for (dbidx = 0; dbidx < mbt->rdfp->num_seqs; ++dbidx) {

        tlen = readdb_get_sequence(mbt->rdfp, dbidx, &t_seq);
        readdb_get_descriptor(mbt->rdfp, dbidx, &subject_sip, &defline);
        
        gal_orig = search(mbt, t_seq, tlen, mbap->prog_flag, &N);
        
        if (N <= 0) {
            MemFree(defline);
            SeqIdSetFree(subject_sip);
            continue; /* No hits found for the sequence */
        }
        
        s_seq = readdb_unpack(t_seq, tlen);
        
        gal_new = NULL;  /* For next sequence should be NULL */
        
        for (gp = gal_orig; gp;) {
            gal_t *ngp;
            Int4 i = binary_search(gp->b2, qdb->beg, qdb->nseq);
            if (gp->b2 < qdb->beg[i])
                ErrPostEx(SEV_WARNING, 0, 0, "mb: warning %d %ld %d\n", 
                          gp->b2, qdb->beg[i], i);
            
            gp->query_id = i;
            gp->subject_id = dbidx;
            
            /* Adjustment coordinates of alignment relatively to start of 
               sequence */
            gp->b2 -= qdb->beg[i];
            gp->e2 -= qdb->beg[i];
            
            /* "Strand" of the query sequence */

            gp->query_strand = 
                qdb_is_revcomp(qdb, i) ? Seq_strand_minus : Seq_strand_plus; 
            
            /* printf("<<***>> %d %d %d %d\n", gp->b1, gp->b2, gp->e1, gp->e2); */
            for (ngp = gal_new; ngp; ngp = ngp->next) {
                if (ngp->query_id == i)
                    if (GAL_CONTAINS(ngp, gp) && GAL_CLOSE(ngp, gp))
                        break;
            }
            
            if (ngp != 0) {
                /* unlink and free */
                gal_t *t = gp->next;
                MemFree(gp);
                gp = t;
            } else {

                Int4 q_len, cmpl;
                
                Int4 j = (mbap->output_al != MB_ALIGN_NONE ?
                          mb_extend_gap_align_script :
                          mb_extend_gap_align) \
                    (gp, mbap->xdrop_ext, mbap->match_score,
                     mbap->mismatch_penalty_mb, mbap->gap_open_penalty, 
                     mbap->gap_extension_penalty, s_seq, tlen,
                     qdb_seq(qdb,i), qdb_length(qdb,i));
                
                gp->prune = 0;
                gp->diag = j;
                for (ngp = gal_new; ngp; ngp = ngp->next) {
                    if (ngp->prune == 0)
                        if (ngp->query_id == i)
                            if (GAL_CONTAINS(gp, ngp) && GAL_CLOSE(ngp, gp))
                                ngp->prune = 1;
                }
                

                ngp = gp->next;
                gp->next= gal_new;
                gal_new = gp;
                gp = ngp;
                
            }
        }
                
        if(mbap->real_time_print) {
            MBPrintResult(mbt, gal_new, qdb, mbap);
            free_gal_list(gal_new);
            gal_new = NULL;
        } else {
            gal_prev = NULL;
            for (gal_tmp = gal_new; gal_tmp != NULL; gal_tmp = gal_tmp->next) {
                
                if(gal_prev != NULL)
                    gal_prev->next = NULL;
                
                MBHitlistAddOne(mbhlp, gal_tmp);
                gal_prev = gal_tmp;
            }
        }

        s_seq = MemFree(s_seq);       /* Plain sequence */
        defline = MemFree(defline);
        SeqIdSetFree(subject_sip);
        
    } /* Loop over all sequences in the blast database */

    if(mbhlp != NULL) {
        MBSortHitList(mbhlp);
    }

    return mbhlp;
}

Int4 est_search(mbt_t *mbt)
{
    Int4 i, dbidx;
    gal_t *gal_orig;
    UcharPtr s;
    Int4 slen;
    CharPtr defline;
    
    if(mbt == NULL)
        return FALSE;
    
    for (dbidx = 0; dbidx < mbt->rdfp->num_seqs; ++dbidx) {
        
        readdb_get_stuff(mbt->rdfp, dbidx, &s, &slen, &defline);

        /* Here is gal_orig == 0 or == 1 !! */
        gal_orig = search(mbt, s, slen, MB_DO_EST, &i);
        
        if (gal_orig != NULL) {
            /* free_gal_list(gal_orig); */
            mb_print_sequence(s, slen, defline);
        }
        defline = (CharPtr)  MemFree(defline);
    }

    return TRUE;
}

static void sieve_mbt(mbt_t *mbt)
{
    const Uchar *const code = fasta_decoding;
    Int4 i;

    for (i = 0; i <= mbt->HASH_SIZE; i++) {
	Int4 h;
	Int4 j = 0;
	for (h = mbt->hashtab[i]; h; h = mbt->next_pos[h])
	    j++;
	if (j > mbt->surfeit) {
	    fprintf(stderr, "%d ", j);
	    for (j = 0; j < mbt->W*4; j++) {
		fprintf(stderr, "%c", code[(i >> (mbt->W*8 - 2 - j*2)) % 4]);
	    }
	    fprintf(stderr, "\n");
	    mbt->hashtab[i] = 0;
	}
    }   
}

static void bld_table(mbt_t *mbt)
{
    const Int1 *const encoding = fasta_encoding;
    register Uint4 ecode;
    register Uint4 mask;
    Int4 j, h;
    register Int4 i;
    register Uchar *t, *p;
    Uchar *s;
    
    /* perform initializations */
    mbt->HASH_SIZE = (1 << (mbt->W*8)) - 1;
    mask = mbt->mask = (1 << (mbt->W*8 - 2)) - 1;
    if (mbt->hashtab == 0) {
	mbt->hashtab = MemNew(sizeof(Int4Ptr) * (mbt->HASH_SIZE + 1));
	mbt->next_pos = MemNew(sizeof(Int4Ptr) * BLOCKSIZE);
    }

    s = mbt->seq;
    t = s - 1;
    ecode = 0L;
    for (i = 1; i < mbt->W*4; ++i)
	ecode = (ecode << 2) + (Uint4) encoding[*++t];
    p = s;
    while (*++t) {
	Uint4 e = encoding[*t];
	if (e == ~0U) {
	    ecode = 0;
	    p = t + mbt->W*4;
	    *t = toupper(*t);
	} else {
	    ecode = ((ecode & mask) << 2) + e;
	    if (t >= p) {
		mbt->next_pos[i] = mbt->hashtab[ecode];
		mbt->hashtab[ecode] = i;
	    }
	}
	i++;
    }
    mbt->mask = (1 << (mbt->W*8 - 8)) - 1;
    /* sieve_mbt(mbt); */
}

mbt_t *mb_init(CharPtr blast_database)
{
    mbt_t *mbt = MemNew(sizeof(*mbt));
    mbt->seq = 0;
    mbt->slen = 0;

    mbt->W = 3;
    mbt->lpm = 20;
    mbt->surfeit = 0;

    mbt->numMSPs = 0;
    mbt->gal_list = 0;
    
    mbt->HASH_SIZE = (1 << 24) - 1;
    mbt->hashtab = 0;
    mbt->next_pos = 0;
    mbt->mask = 0;
    mbt->stack_index = 0;
    
    if ((mbt->rdfp = readdb_new(blast_database, 0)) == 0) {
	ErrPostEx(SEV_ERROR, 0, 0, 
                  "cannot open database '%s'", blast_database);
        return NULL;
    }
    
    return mbt;
}

void mb_start(mbt_t *mbt, UcharPtr dna, Int4 len, Int4 w, Int4 p)
{
    mbt->seq = dna;
    mbt->slen = len;
    mbt->lpm = w;
    mbt->surfeit = p;
    
    bld_table(mbt);
}

void mb_end(mbt_t *mbt)
{

    readdb_destruct(mbt->rdfp);
    mbt->hashtab = MemFree(mbt->hashtab);
    mbt->next_pos = MemFree(mbt->next_pos);
    
    mbt = MemFree(mbt);
}

static Int4 gap_free_extend(mbt_t *mbt, gal_t *mp, const UcharPtr seq1, 
                            Int4 len1, const UcharPtr seq2, Int4 len2)
{
    const Int4 X = mbt->x;
    const Uchar *s, *q;
    Int4 left_sum, right_sum, sum, score;
    const Uchar *const seq1_end = seq1 + len1;
    const Uchar *const seq2_end = seq2 + len2;
    
    score = (mp->e1 - mp->b1 + 1)*mbt->ss['A']['A'];
    if (score > mbt->k)
	return 1;

    score = mbt->k - score;
    right_sum = sum = 0;
    q = seq1 + mp->e1;
    s = seq2 + mp->e2;
    while (s < seq2_end && q < seq1_end && sum >= right_sum - X)
	if ((sum += mbt->ss[*s++][*q++]) > right_sum) {
	    right_sum = sum;
	    if (right_sum > score)
		return 1;
	}
    score -= right_sum;
    /* extend to the left */
    left_sum = sum = 0;
    q = seq1 + mp->b1;
    s = seq2 + mp->b2;
    while (s > seq2 && q > seq1 && sum >= left_sum - X)
	if ((sum += mbt->ss[*--s][*--q]) > left_sum) {
	    left_sum = sum;
	    if (left_sum > score)
		return 1;
	}
    return 0;
}
