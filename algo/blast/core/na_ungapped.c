/* $Id: na_ungapped.c,v 1.21 2009/07/30 19:34:30 kazimird Exp $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
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
 */

/** @file na_ungapped.c
 * Nucleotide ungapped extension routines
 */

#ifndef SKIP_DOXYGEN_PROCESSING
static char const rcsid[] =
    "$Id: na_ungapped.c,v 1.21 2009/07/30 19:34:30 kazimird Exp $";
#endif                          /* SKIP_DOXYGEN_PROCESSING */

#include <algo/blast/core/na_ungapped.h>
#include <algo/blast/core/blast_nalookup.h>
#include <algo/blast/core/blast_nascan.h>
#include <algo/blast/core/mb_indexed_lookup.h>
#include <algo/blast/core/blast_util.h> /* for NCBI2NA_UNPACK_BASE macros */

#include "index_ungapped.h"
#include "masksubj.inl"

/** Perform ungapped extension of a word hit, using a score
 *  matrix and extending one base at a time
 * @param query The query sequence [in]
 * @param subject The subject sequence [in]
 * @param matrix The scoring matrix [in]
 * @param q_off The offset of a word in query [in]
 * @param s_off The offset of a word in subject [in]
 * @param X The drop-off parameter for the ungapped extension [in]
 * @param ungapped_data The ungapped extension information [out]
 */
static void
s_NuclUngappedExtendExact(BLAST_SequenceBlk * query,
                          BLAST_SequenceBlk * subject, Int4 ** matrix,
                          Int4 q_off, Int4 s_off, Int4 X,
                          BlastUngappedData * ungapped_data)
{
    Uint1 *q;
    Int4 sum, score;
    Uint1 ch;
    Uint1 *subject0, *sf, *q_beg, *q_end, *s, *start;
    Int2 remainder, base;
    Int4 q_avail, s_avail;

    base = 3 - (s_off % 4);

    subject0 = subject->sequence;
    q_avail = query->length - q_off;
    s_avail = subject->length - s_off;

    q = q_beg = q_end = query->sequence + q_off;
    s = subject0 + s_off / COMPRESSION_RATIO;
    if (q_off < s_off) {
        start = subject0 + (s_off - q_off) / COMPRESSION_RATIO;
        remainder = 3 - ((s_off - q_off) % COMPRESSION_RATIO);
    } else {
        start = subject0;
        remainder = 3;
    }

    score = 0;
    sum = 0;

    /* extend to the left */
    while ((s > start) || (s == start && base < remainder)) {
        if (base == 3) {
            s--;
            base = 0;
        } else {
            ++base;
        }
        ch = *s;
        if ((sum += matrix[*--q][NCBI2NA_UNPACK_BASE(ch, base)]) > 0) {
            q_beg = q;
            score += sum;
            sum = 0;
        } else if (sum < X) {
            break;
        }
    }

    ungapped_data->q_start = q_beg - query->sequence;
    ungapped_data->s_start = s_off - (q_off - ungapped_data->q_start);

    if (q_avail < s_avail) {
        sf = subject0 + (s_off + q_avail) / COMPRESSION_RATIO;
        remainder = 3 - ((s_off + q_avail) % COMPRESSION_RATIO);
    } else {
        sf = subject0 + (subject->length) / COMPRESSION_RATIO;
        remainder = 3 - ((subject->length) % COMPRESSION_RATIO);
    }

    /* extend to the right */
    q = query->sequence + q_off;
    s = subject0 + s_off / COMPRESSION_RATIO;
    sum = 0;
    base = 3 - (s_off % COMPRESSION_RATIO);

    while (s < sf || (s == sf && base > remainder)) {
        ch = *s;
        if ((sum += matrix[*q++][NCBI2NA_UNPACK_BASE(ch, base)]) > 0) {
            q_end = q;
            score += sum;
            sum = 0;
        } else if (sum < X)
            break;
        if (base == 0) {
            base = 3;
            s++;
        } else
            base--;
    }

    ungapped_data->length = q_end - q_beg;
    ungapped_data->score = score;
}

/** Perform ungapped extension of a word hit. Use an approximate method
 * and revert to rigorous ungapped alignment if the approximate score
 * is high enough
 * @param query The query sequence [in]
 * @param subject The subject sequence [in]
 * @param matrix The scoring matrix [in]
 * @param q_off The offset of a word in query [in]
 * @param s_match_end The first offset in the subject sequence that 
 *              is not known to exactly match the query sequence [in]
 * @param s_off The offset of a word in subject [in]
 * @param X The drop-off parameter for the ungapped extension [in]
 * @param ungapped_data The ungapped extension information [out]
 * @param score_table Array containing precompute sums of
 *                    match and mismatch scores [in]
 * @param reduced_cutoff Score beyond which a rigorous extension is used [in]
 */
static void
s_NuclUngappedExtend(BLAST_SequenceBlk * query,
                     BLAST_SequenceBlk * subject, Int4 ** matrix,
                     Int4 q_off, Int4 s_match_end, Int4 s_off,
                     Int4 X, BlastUngappedData * ungapped_data,
                     const Int4 * score_table, Int4 reduced_cutoff)
{
    Uint1 *q_start = query->sequence;
    Uint1 *s_start = subject->sequence;
    Uint1 *q;
    Uint1 *s;
    Int4 sum, score;
    Int4 i, len;
    Uint1 *new_q;
    Int4 q_ext, s_ext;

    /* The left extension begins behind (q_ext,s_ext); this is the first
       4-base boundary after s_off. */

    len = (COMPRESSION_RATIO - (s_off % COMPRESSION_RATIO)) %
        COMPRESSION_RATIO;
    q_ext = q_off + len;
    s_ext = s_off + len;
    q = q_start + q_ext;
    s = s_start + s_ext / COMPRESSION_RATIO;
    len = MIN(q_ext, s_ext) / COMPRESSION_RATIO;

    score = 0;
    sum = 0;
    new_q = q;

    for (i = 0; i < len; s--, q -= 4, i++) {
        Uint1 s_byte = s[-1];
        Uint1 q_byte = (q[-4] << 6) | (q[-3] << 4) | (q[-2] << 2) | q[-1];

        sum += score_table[q_byte ^ s_byte];
        if (sum > 0) {
            new_q = q - 4;
            score += sum;
            sum = 0;
        }
        if (sum < X) {
            break;
        }
    }

    /* record the start point of the extension */

    ungapped_data->q_start = new_q - q_start;
    ungapped_data->s_start = s_ext - (q_ext - ungapped_data->q_start);

    /* the right extension begins at the first bases not examined by the
       previous loop */

    q = q_start + q_ext;
    s = s_start + s_ext / COMPRESSION_RATIO;
    len = MIN(query->length - q_ext, subject->length - s_ext) /
        COMPRESSION_RATIO;
    sum = 0;
    new_q = q;

    for (i = 0; i < len; s++, q += 4, i++) {
        Uint1 s_byte = s[0];
        Uint1 q_byte = (q[0] << 6) | (q[1] << 4) | (q[2] << 2) | q[3];

        sum += score_table[q_byte ^ s_byte];
        if (sum > 0) {
            new_q = q + 3;
            score += sum;
            sum = 0;
        }
        if (sum < X) {
            break;
        }
    }

    if (score >= reduced_cutoff) {
        /* the current score is high enough; throw away the alignment and
           recompute it rigorously */
        s_NuclUngappedExtendExact(query, subject, matrix, q_off,
                                  s_off, X, ungapped_data);
    } else {
        /* record the length and score of the extension. Make sure the
           alignment extends at least to s_match_end */
        ungapped_data->score = score;
        ungapped_data->length = MAX(s_match_end - ungapped_data->s_start,
                                    (new_q - q_start) -
                                    ungapped_data->q_start + 1);
    }
}

/**
 * Attempt to retrieve information associated with diagonal diag.
 * @param table The hash table [in]
 * @param diag The diagonal to be retrieved [in]
 * @param level The offset of the last hit on the specified diagonal [out]
 * @param hit_saved Whether or not the last hit on the specified diagonal was saved [out]
 * @param hit_length length of the last hit on the specified diagonal [out]
 * @return 1 if successful, 0 if no hit was found on the specified diagonal.
 */
static NCBI_INLINE Int4 s_BlastDiagHashRetrieve(BLAST_DiagHash * table,
                                                Int4 diag, Int4 * level,
                                                Int4 * hit_len,
                                                Int4 * hit_saved)
{
    /* see http://lxr.linux.no/source/include/linux/hash.h */
    /* mod operator will be strength-reduced to an and by the compiler */
    Uint4 bucket = ((Uint4) diag * 0x9E370001) % DIAGHASH_NUM_BUCKETS;
    Uint4 index = table->backbone[bucket];

    while (index) {
        if (table->chain[index].diag == diag) {
            *level = table->chain[index].level;
            *hit_len = table->chain[index].hit_len;
            *hit_saved = table->chain[index].hit_saved;
            return 1;
        }

        index = table->chain[index].next;
    }
    return 0;
}

/**
 * Attempt to store information associated with diagonal diag.
 * Cleans up defunct entries along the way or allocates more space if necessary.
 * @param table The hash table [in]
 * @param diag The diagonal to be stored [in]
 * @param level The offset of the hit to be stored [in]
 * @param len The length of the hit to be stored [in]
 * @param hit_saved Whether or not this hit was stored [in]
 * @param s_end Needed to clean up defunct entries [in]
 * @param window_size Needed to clean up defunct entries [in]
 * @return 1 if successful, 0 if memory allocation failed.
 */
static NCBI_INLINE Int4 s_BlastDiagHashInsert(BLAST_DiagHash * table,
                                              Int4 diag, Int4 level,
                                              Int4 len,
                                              Int4 hit_saved,
                                              Int4 s_off,
                                              Int4 window_size)
{
    Uint4 bucket = ((Uint4) diag * 0x9E370001) % DIAGHASH_NUM_BUCKETS;
    Uint4 index = table->backbone[bucket];
    DiagHashCell *cell = NULL;

    while (index) {
        /* if we find what we're looking for, save into it */
        if (table->chain[index].diag == diag) {
            table->chain[index].level = level;
            table->chain[index].hit_len = len;
            table->chain[index].hit_saved = hit_saved;
            return 1;
        }
        /* otherwise, if this hit is stale, save into it. */
        else {
            /* if this hit is stale, save into it. */
            if ( s_off - table->chain[index].level > window_size) {
                table->chain[index].diag = diag;
                table->chain[index].level = level;
                table->chain[index].hit_len = len;
                table->chain[index].hit_saved = hit_saved;
                return 1;
            }
        }
        index = table->chain[index].next;
    }

    /* if we got this far, we were unable to replace any existing entries. */

    /* if there's no more room, allocate more */
    if (table->occupancy == table->capacity) {
        table->capacity *= 2;
        table->chain =
            realloc(table->chain, table->capacity * sizeof(DiagHashCell));
        if (table->chain == NULL)
            return 0;
    }

    cell = table->chain + table->occupancy;
    cell->diag = diag;
    cell->level = level;
    cell->hit_len = len;
    cell->hit_saved = hit_saved;
    cell->next = table->backbone[bucket];
    table->backbone[bucket] = table->occupancy;
    table->occupancy++;

    return 1;
}

/** Extends the 1st hit to the right.  Used only in 2-hit algorithm.
 * This is called when 1st hit is identified.  It will scan areas that
 * are relevant to 2-hit algorithm (overlap_size, window_size, etc.),
 * which are nevertheless skipped by the strided scan.
 * @param query Query sequence data [in]
 * @param subject Subject sequence data [in]
 * @param q_end Query ending offset (starting point for right extension) [in]
 * @param s_end Subject ending offset (starting point for right extension) [in]
 * @param query_info Structure containing query context ranges [in]
 * @param max_ext  The maximum number of bases to extend [in]
 * @return Actual number of bases extended 
 */
static Uint4
s_BlastRightExtend(BLAST_SequenceBlk * query,
                   BLAST_SequenceBlk * subject,
                   Uint4 q_end, Uint4 s_end,
                   BlastQueryInfo * query_info,
                   Uint4 max_ext)
{
    Uint1 *q = query->sequence + q_end;
    Uint1 *s = subject->sequence + s_end / COMPRESSION_RATIO;
    Uint4 extended;

    /* reduce the max_ext to prevent hitting the sequence boundary */
    /* note: we probably can return immediately if full extent can not be achieved*/
    /* TODO: subject masking is ignored, to be fixed */
    if (query_info) {
        Uint4 context = BSearchContextInfo(q_end, query_info);
        Uint4 context_i = query_info->contexts[context].query_offset;
        Uint4 context_f = context_i + query_info->contexts[context].query_length;
        Uint4 subject_length = subject->length;
        max_ext = MIN(max_ext, MIN(context_f - q_end, subject_length - s_end));
    }
 
    for (extended = 0; extended < max_ext; extended++, q++) {
        if (((Uint1) (*s << (2 * (s_end % COMPRESSION_RATIO))) >> 6) != *q)
            break;
        s_end++;
        if (s_end % COMPRESSION_RATIO == 0) s++;
    }
    return extended;
}

/** Perform ungapped extension given an offset pair, and save the initial 
 * hit information if the hit qualifies. This function assumes that the
 * exact match has already been extended to the word size parameter. It also
 * supports a two-hit version, where extension is performed only after two hits
 * are detected within a window.
 * @param query The query sequence [in]
 * @param subject The subject sequence [in]
 * @param contiguous Is this a contiguous blastn? [in]
 * @param word_params The parameters related to initial word extension [in]
 * @param matrix the substitution matrix for ungapped extension [in]
 * @param query_info Structure containing query context ranges [in]
 * @param diag_table Structure containing diagonal table with word extension 
 *                   information [in]
 * @param q_off The offset in the query sequence [in]
 * @param s_off The offset in the subject sequence [in]
 * @param word_length The length of the hit [in]
 * @param init_hitlist The structure containing information about all 
 *                     initial hits [in] [out]
 * @return 1 if hit was extended, 0 if not
 */
static Int4
s_BlastnDiagTableExtendInitialHit(BLAST_SequenceBlk * query,
                             BLAST_SequenceBlk * subject, 
                             Boolean contiguous,
                             const BlastInitialWordParameters * word_params,
                             Int4 ** matrix, BlastQueryInfo * query_info,
                             BLAST_DiagTable * diag_table,
                             Int4 q_off, Int4 s_off, Int4 word_length,
                             BlastInitHitList * init_hitlist)
{
    Int4 diag, real_diag;
    Int4 s_end, s_off_pos, s_end_pos;
    Int4 ext_right = 0;
    BlastUngappedData *ungapped_data;
    BlastUngappedData dummy_ungapped_data;
    Int4 window_size = word_params->options->window_size;
    Int4 hit_ready = 1;
    Int4 last_hit, hit_saved;
    DiagStruct *hit_level_array;
    BlastUngappedCutoffs *cutoffs = NULL;
    Boolean two_hits = (window_size > 0);
    Boolean found = FALSE;
    Int4 Delta = MIN(5, window_size - word_length);

    hit_level_array = diag_table->hit_level_array;
    ASSERT(hit_level_array);

    diag = s_off + diag_table->diag_array_length - q_off;
    real_diag = diag & diag_table->diag_mask;
    last_hit = hit_level_array[real_diag].last_hit;
    hit_saved = hit_level_array[real_diag].flag;
    s_end = s_off + word_length;
    s_off_pos = s_off + diag_table->offset;
    s_end_pos = s_end + diag_table->offset;

    /* hit within the explored area should be rejected*/
    if (s_off_pos < last_hit) return 0;

    if (two_hits && (hit_saved || s_end_pos > last_hit + window_size )) {
        /* check to see if it can be extended to the right by
           word_length and therefore qualifies for a double-hit */
        if (contiguous) {
            ext_right = s_BlastRightExtend(query, subject,
                  q_off + word_length, s_end, query_info, word_length);
            /* update the right end*/
            s_end += ext_right;
            s_end_pos += ext_right;
        }

        if (ext_right < word_length) {
            /* try off-diagonals */
            Int4 orig_diag = real_diag + diag_table->diag_array_length;
            Int4 s_a = s_off_pos + word_length - window_size;
            Int4 s_b = s_end_pos - 2 * word_length;
            Int4 delta;
            if (Delta < 0) Delta = 0;
            for (delta = 1; delta < Delta ; ++delta) {
                Int4 off_diag  = (orig_diag + delta) & diag_table->diag_mask;
                Int4 off_s_end = hit_level_array[off_diag].last_hit;
                Int4 off_s_l   = diag_table->hit_len_array[off_diag];
                if ( off_s_l
                 && off_s_end - delta >= s_a 
                 && off_s_end - off_s_l <= s_b) {
                    found = TRUE;
                    break;
                }
                off_diag  = (orig_diag - delta) & diag_table->diag_mask;
                off_s_end = hit_level_array[off_diag].last_hit;
                off_s_l   = diag_table->hit_len_array[off_diag];
                if ( off_s_l
                 && off_s_end >= s_a 
                 && off_s_end - off_s_l + delta <= s_b) {
                    found = TRUE;
                    break;
                }
            }
            if (!found) {
                /* This is a new hit */
                hit_ready = 0;
            }
        }
    }

    if (hit_ready) {
        if (word_params->ungapped_extension) {
            Int4 context = BSearchContextInfo(q_off, query_info);
            cutoffs = word_params->cutoffs + context;
            ungapped_data = &dummy_ungapped_data;
            s_NuclUngappedExtend(query, subject, matrix, q_off, s_end, s_off,
                                 -(cutoffs->x_dropoff), ungapped_data,
                                 word_params->nucl_score_table,
                                 cutoffs->reduced_nucl_cutoff_score);

            if (found || ungapped_data->score >= cutoffs->cutoff_score) {
                BlastUngappedData *final_data =
                    (BlastUngappedData *) malloc(sizeof(BlastUngappedData));
                *final_data = *ungapped_data;
                BLAST_SaveInitialHit(init_hitlist, q_off, s_off, final_data);
                s_end_pos = ungapped_data->length + ungapped_data->s_start
                          + diag_table->offset;
            } else {
                hit_ready = 0;
            }
        } else {
            ungapped_data = NULL;
            BLAST_SaveInitialHit(init_hitlist, q_off, s_off, ungapped_data);
        }
    } 

    hit_level_array[real_diag].last_hit = s_end_pos;
    hit_level_array[real_diag].flag = hit_ready;
    if (two_hits) {
        diag_table->hit_len_array[real_diag] = (hit_ready) ? 0 : s_end_pos - s_off_pos;
    }

    return hit_ready;
}

/** Perform ungapped extension given an offset pair, and save the initial 
 * hit information if the hit qualifies. This function assumes that the
 * exact match has already been extended to the word size parameter. It also
 * supports a two-hit version, where extension is performed only after two hits
 * are detected within a window.
 * @param query The query sequence [in]
 * @param subject The subject sequence [in]
 * @param contiguous Is this a contiguous blastn? [in]
 * @param word_params The parameters related to initial word extension [in]
 * @param matrix the substitution matrix for ungapped extension [in]
 * @param query_info Structure containing query context ranges [in]
 * @param hash_table Structure containing initial hits [in] [out]
 * @param q_off The offset in the query sequence [in]
 * @param s_off The offset in the subject sequence [in]
 * @param word_length The length of the hit [in]
 * @param init_hitlist The structure containing information about all 
 *                     initial hits [in] [out]
 * @return 1 if hit was extended, 0 if not
 */
static Int4
s_BlastnDiagHashExtendInitialHit(BLAST_SequenceBlk * query,
                               BLAST_SequenceBlk * subject, 
                               Boolean contiguous, 
                               const BlastInitialWordParameters * word_params,
                               Int4 ** matrix, BlastQueryInfo * query_info,
                               BLAST_DiagHash * hash_table,
                               Int4 q_off, Int4 s_off, Int4 word_length,
                               BlastInitHitList * init_hitlist)
{
    Int4 diag;
    Int4 s_end, s_off_pos, s_end_pos, s_l;
    Int4 ext_right = 0;
    BlastUngappedData *ungapped_data;
    BlastUngappedData dummy_ungapped_data;
    Int4 window_size = word_params->options->window_size;
    Int4 hit_ready = 1; 
    Int4 last_hit, hit_saved = 0;
    BlastUngappedCutoffs *cutoffs = NULL;
    Boolean two_hits = (window_size > 0);
    Boolean found = FALSE;
    Int4 Delta = MIN(5, window_size - word_length);
    Int4 rc;

    diag = s_off - q_off;
    s_end = s_off + word_length;
    s_off_pos = s_off + hash_table->offset;
    s_end_pos = s_end + hash_table->offset;

    rc = s_BlastDiagHashRetrieve(hash_table, diag, &last_hit, &s_l, &hit_saved);

    /* if there is no record in hashtable, we set last_hit to be a very negative number */
    if(!rc)  last_hit = 0;

    /* hit within the explored area should be rejected*/
    if (s_off_pos < last_hit) return 0;

    if (two_hits && (hit_saved || s_end_pos > last_hit + window_size )) {
        /* this must be the 1st hit */
        /* check to see if it can be extended to the right by
           word_length and therefore qualifies for a double-hit */
        if (contiguous) {
            ext_right = s_BlastRightExtend(query, subject,
                  q_off + word_length, s_end, query_info, word_length);
            /* update the right end*/
            s_end += ext_right;
            s_end_pos += ext_right;
        }

        if (ext_right < word_length) {
            /* try off-diagonal */
            Int4 s_a = s_off_pos + word_length - window_size;
            Int4 s_b = s_end_pos - 2 * word_length;
            Int4 delta;
            if (Delta < 0) Delta = 0;
            for (delta = 1; delta < Delta; ++delta) {
                Int4 off_s_end = 0;
                Int4 off_s_l = 0;
                Int4 off_hit_saved = 0;
                Int4 off_rc = s_BlastDiagHashRetrieve(hash_table, diag + delta, 
                              &off_s_end, &off_s_l, &off_hit_saved);
                if ( off_rc
                  && off_s_l
                  && off_s_end - delta >= s_a
                  && off_s_end - off_s_l <= s_b) {
                     found = TRUE;
                     break;
                }
                off_rc = s_BlastDiagHashRetrieve(hash_table, diag - delta, 
                              &off_s_end, &off_s_l, &off_hit_saved);
                if ( off_rc
                  && off_s_l
                  && off_s_end >= s_a
                  && off_s_end - off_s_l + delta <= s_b) {
                     found = TRUE;
                     break;
                }
            }
            if (!found) {
                /* This is a new hit */
                hit_ready = 0;
            }
        }
    }

    if (hit_ready) {
        if (word_params->ungapped_extension) {
            /* Perform ungapped extension */
            Int4 context = BSearchContextInfo(q_off, query_info);
            cutoffs = word_params->cutoffs + context;
            ungapped_data = &dummy_ungapped_data;
            s_NuclUngappedExtend(query, subject, matrix, q_off, s_end,
                                 s_off, -(cutoffs->x_dropoff),
                                 ungapped_data,
                                 word_params->nucl_score_table,
                                 cutoffs->reduced_nucl_cutoff_score);
            if (found || ungapped_data->score >= cutoffs->cutoff_score) {
                BlastUngappedData *final_data =
                    (BlastUngappedData *) malloc(sizeof(BlastUngappedData));
                *final_data = *ungapped_data;
                BLAST_SaveInitialHit(init_hitlist, q_off, s_off, final_data);
                s_end_pos = ungapped_data->length + ungapped_data->s_start 
                          + hash_table->offset;
            } else {
                hit_ready = 0;
            }
        } else {
            ungapped_data = NULL;
            BLAST_SaveInitialHit(init_hitlist, q_off, s_off, ungapped_data);
        }
    } 
    
    s_BlastDiagHashInsert(hash_table, diag, s_end_pos, 
                          (hit_ready) ? 0 : s_end_pos - s_off_pos,
                          hit_ready, s_off_pos, window_size + Delta);

    return hit_ready;
}

/** Check that word just found is not part of a soft-masked region.
 * This check only happens for soft masking for certain word sizes
 * @param locations not masked regions [in]
 * @param query query sequence [in]
 * @param subject subject sequence [in]
 * @param q_off query start [in]
 * @param word_length length of word [in]
 * @param s_end subject end (starting pos to extend) [in]
 * @param max_ext maxium possible right extension [in]
 */
static Boolean 
s_IsWordSeedMasked(BlastSeqLoc* locations, 
                   BLAST_SequenceBlk * query,
                   BLAST_SequenceBlk * subject,
                   Uint4 q_off, 
                   Uint4 word_length,
                   Uint4 s_end,
                   Uint4 max_ext)
{
    Int4 q_end = q_off + word_length;
    while (locations) {
       if (q_end-1 >= locations->ssr->left && q_end-1 <= locations->ssr->right) {
          /* since we have exhausted left extension, this hit is hopeless */
          return TRUE;
       } 
       if (q_off >= locations->ssr->left && q_off <= locations->ssr->right) {
          /* determine the right boundary for extension*/
          if (locations->next) {
              max_ext = MIN (max_ext, locations->next->ssr->left - q_end);
          } 
          if (locations->ssr->right+1 - q_off > max_ext) {
              /*  query range too short to contain a word  */
              return TRUE;
          }
          /* try to extend to the right */
          max_ext = locations->ssr->right+1 - q_off;
          if (s_BlastRightExtend(query, subject, q_end, s_end, NULL, max_ext) < max_ext) {
              return TRUE;
          }
          /* extension sucessful, (TODO should we update new q_off? )*/
          return FALSE;
       }
       locations = locations->next;
    }
    return FALSE;
}

/** Perform ungapped extensions on the hits retrieved from
 * blastn/megablast lookup tables, skipping the mini-extension process
 * @param offset_pairs Array of query and subject offsets. [in]
 * @param num_hits Size of the above arrays [in]
 * @param word_params Parameters for word extension [in]
 * @param lookup_wrap Lookup table wrapper structure [in]
 * @param query Query sequence data [in]
 * @param subject Subject sequence data [in]
 * @param matrix Scoring matrix for ungapped extension [in]
 * @param query_info Structure containing query context ranges [in]
 * @param ewp Word extension structure containing information about the 
 *            extent of already processed hits on each diagonal [in]
 * @param init_hitlist Structure to keep the extended hits. 
 *                     Must be allocated outside of this function [in] [out]
 * @param range the subject range [in]
 * @return Number of hits extended. 
 */
static Int4
s_BlastNaExtendDirect(const BlastOffsetPair * offset_pairs, Int4 num_hits,
                      const BlastInitialWordParameters * word_params,
                      LookupTableWrap * lookup_wrap,
                      BLAST_SequenceBlk * query,
                      BLAST_SequenceBlk * subject, Int4 ** matrix,
                      BlastQueryInfo * query_info,
                      Blast_ExtendWord * ewp,
                      BlastInitHitList * init_hitlist,
                      Int4 * range)
{
    Int4 index;
    Int4 hits_extended = 0;
    Int4 word_length;
    Boolean contiguous = TRUE;

    if (lookup_wrap->lut_type == eMBLookupTable) {
        BlastMBLookupTable *lut = (BlastMBLookupTable *) lookup_wrap->lut;
        word_length = lut->word_length;
        if (lut->discontiguous) {
            word_length = lut->template_length;
            contiguous = FALSE;
        }
        ASSERT(word_length == lut->lut_word_length || lut->discontiguous);
    } 
    else if (lookup_wrap->lut_type == eSmallNaLookupTable) {
        BlastSmallNaLookupTable *lut = 
                        (BlastSmallNaLookupTable *) lookup_wrap->lut;
        word_length = lut->word_length;
    } 
    else {
        BlastNaLookupTable *lut = (BlastNaLookupTable *) lookup_wrap->lut;
        word_length = lut->word_length;
    }

    if (word_params->container_type == eDiagHash) {
        for (index = 0; index < num_hits; ++index) {
            Uint4 s_offset = offset_pairs[index].qs_offsets.s_off;
            Uint4 q_offset = offset_pairs[index].qs_offsets.q_off;

            hits_extended += s_BlastnDiagHashExtendInitialHit(query, subject, 
                                                   contiguous, 
                                                   word_params, matrix,
                                                   query_info, ewp->hash_table,
                                                   q_offset, 
                                                   s_offset,  word_length,
                                                   init_hitlist);
        }
    } 
    else {
        for (index = 0; index < num_hits; ++index) {
            Uint4 s_offset = offset_pairs[index].qs_offsets.s_off;
            Uint4 q_offset = offset_pairs[index].qs_offsets.q_off;

            hits_extended += s_BlastnDiagTableExtendInitialHit(query, subject, 
                                                 contiguous,
                                                 word_params, matrix, 
                                                 query_info, ewp->diag_table, 
                                                 q_offset,
                                                 s_offset,  word_length,
                                                 init_hitlist);
        }
    }
    return hits_extended;
}

/** Perform exact match extensions on the hits retrieved from
 * blastn/megablast lookup tables, assuming an arbitrary number of bases 
 * in a lookup and arbitrary start offset of each hit. Also 
 * update the diagonal structure.
 * @param offset_pairs Array of query and subject offsets [in]
 * @param num_hits Size of the above arrays [in]
 * @param word_params Parameters for word extension [in]
 * @param lookup_wrap Lookup table wrapper structure [in]
 * @param query Query sequence data [in]
 * @param subject Subject sequence data [in]
 * @param matrix Scoring matrix for ungapped extension [in]
 * @param query_info Structure containing query context ranges [in]
 * @param ewp Word extension structure containing information about the 
 *            extent of already processed hits on each diagonal [in]
 * @param init_hitlist Structure to keep the extended hits. 
 *                     Must be allocated outside of this function [in] [out]
 * @param range The subject range
 * @return Number of hits extended. 
 */
static Int4
s_BlastNaExtend(const BlastOffsetPair * offset_pairs, Int4 num_hits,
                const BlastInitialWordParameters * word_params,
                LookupTableWrap * lookup_wrap,
                BLAST_SequenceBlk * query,
                BLAST_SequenceBlk * subject, Int4 ** matrix,
                BlastQueryInfo * query_info,
                Blast_ExtendWord * ewp,
                BlastInitHitList * init_hitlist,
                Int4 * range)
{
    Int4 index;
    Uint1 *q_start = query->sequence;
    Uint1 *s_start = subject->sequence;
    Uint4 word_length, lut_word_length;
    Uint4 q_off, s_off;
    Uint4 max_bases_left, max_bases_right;
    Uint4 extended_left, extended_right;
    Uint4 extra_bases;
    Uint1 *q, *s;
    Int4 hits_extended = 0;
    BlastSeqLoc* masked_locations = NULL;
    Int4 max_ext;

    if (lookup_wrap->lut_type == eMBLookupTable) {
        BlastMBLookupTable *lut = (BlastMBLookupTable *) lookup_wrap->lut;
        word_length = lut->word_length;
        lut_word_length = lut->lut_word_length;
        masked_locations = lut->masked_locations;
    } 
    else {
        BlastNaLookupTable *lut = (BlastNaLookupTable *) lookup_wrap->lut;
        word_length = lut->word_length;
        lut_word_length = lut->lut_word_length;
        masked_locations = lut->masked_locations;
    }
    extra_bases = word_length - lut_word_length;

    /* We trust that the bases of the hit itself are exact matches, 
       and look only for exact matches before and after the hit.

       Most of the time, the lookup table width is close to the word size 
       so only a few bases need examining. Also, most of the time (for
       random sequences) extensions will fail almost immediately (the
       first base examined will not match about 3/4 of the time). Thus it 
       is critical to reduce as much as possible all work that is not the 
       actual examination of sequence data */

    for (index = 0; index < num_hits; ++index) {
        Uint4 s_offset = offset_pairs[index].qs_offsets.s_off;
        Uint4 q_offset = offset_pairs[index].qs_offsets.q_off;

        /* begin with the left extension; the initialization is slightly
           faster. Point to the first base of the lookup table hit and
           work backwards */

        Int4 context = BSearchContextInfo(q_offset, query_info);
        Int4 context_i = query_info->contexts[context].query_offset;
        Int4 context_f = context_i + query_info->contexts[context].query_length;
        max_bases_left = MIN(extra_bases, MIN(q_offset - context_i, s_offset - range[0]));
        q = q_start + q_offset;
        s = s_start + s_offset / COMPRESSION_RATIO;
        s_off = s_offset;

        for (extended_left = 0; extended_left < max_bases_left;
             extended_left++) {
            s_off--;
            q--;
            if (s_off % COMPRESSION_RATIO == 3)
                s--;
            if (((Uint1) (*s << (2 * (s_off % COMPRESSION_RATIO))) >> 6)
                != *q)
                break;
        }

        /* do the right extension if the left extension did not find all
           the bases required. Begin at the first base beyond the lookup
           table hit and move forwards */

        q_off = q_offset + lut_word_length;
        s_off = s_offset + lut_word_length;
        max_ext = MIN(context_f - q_off, range[1] - s_off);

        if (extended_left < extra_bases) {
            q = q_start + q_off;
            s = s_start + s_off / COMPRESSION_RATIO;
            max_bases_right = MIN(extra_bases - extended_left, max_ext);

            for (extended_right = 0; extended_right < max_bases_right;
                 extended_right++) {
                if (((Uint1) (*s << (2 * (s_off % COMPRESSION_RATIO))) >>
                     6) != *q)
                    break;
                s_off++;
                q++;
                if (s_off % COMPRESSION_RATIO == 0)
                    s++;
            }

            /* check if enough extra matches were found */

            if (extended_left + extended_right < extra_bases)
                continue;
            max_ext -= extended_right;
        }
        
        /* Check that the exact match was not masked out. */
        if (word_length > lut_word_length && 
            s_IsWordSeedMasked(masked_locations, query, subject, q_offset-extended_left,
                               word_length, s_off, max_ext)) {
            continue; 
        }

        /* check the diagonal on which the hit lies. The boundaries
           extend from the first match of the hit to one beyond the last
           match */

        if (word_params->container_type == eDiagHash) {
            hits_extended += s_BlastnDiagHashExtendInitialHit(query, subject, 
                                               TRUE,
                                               word_params, matrix,
                                               query_info, ewp->hash_table,
                                               q_offset - extended_left,
                                               s_offset - extended_left,
                                               word_length,
                                               init_hitlist);
        } else {
            hits_extended += s_BlastnDiagTableExtendInitialHit(query, subject, 
                                               TRUE,
                                               word_params, matrix, query_info,
                                               ewp->diag_table,
                                               q_offset - extended_left,
                                               s_offset - extended_left,
                                               word_length,
                                               init_hitlist);
        }
    }
    return hits_extended;
}

/** Perform exact match extensions on the hits retrieved from
 * blastn/megablast lookup tables, assuming the number of bases in a lookup
 * table word, and the start offset of each hit, is a multiple
 * of 4. Also update the diagonal structure.
 * @param offset_pairs Array of query and subject offsets. [in]
 * @param num_hits Size of the above arrays [in]
 * @param word_params Parameters for word extension [in]
 * @param lookup_wrap Lookup table wrapper structure [in]
 * @param query Query sequence data [in]
 * @param subject Subject sequence data [in]
 * @param matrix Scoring matrix for ungapped extension [in]
 * @param query_info Structure containing query context ranges [in]
 * @param ewp Word extension structure containing information about the 
 *            extent of already processed hits on each diagonal [in]
 * @param init_hitlist Structure to keep the extended hits. 
 *                     Must be allocated outside of this function [in] [out]
 * @param range The subject range [in]
 * @return Number of hits extended. 
 */
static Int4
s_BlastNaExtendAligned(const BlastOffsetPair * offset_pairs, Int4 num_hits,
                       const BlastInitialWordParameters * word_params,
                       LookupTableWrap * lookup_wrap,
                       BLAST_SequenceBlk * query, BLAST_SequenceBlk * subject,
                       Int4 ** matrix, BlastQueryInfo * query_info,
                       Blast_ExtendWord * ewp,
                       BlastInitHitList * init_hitlist,
                       Int4 * range)
{
    Int4 index;
    Uint1 *q_start = query->sequence;
    Uint1 *s_start = subject->sequence;
    Uint4 word_length, lut_word_length;
    Uint4 q_off, s_off;
    Uint4 max_bases_left, max_bases_right;
    Uint4 extended_left, extended_right;
    Uint4 extra_bases;
    Uint1 *q, *s;
    Int4 hits_extended = 0;
    BlastSeqLoc* masked_locations = NULL;
    Int4 max_ext;

    if (lookup_wrap->lut_type == eMBLookupTable) {
        BlastMBLookupTable *lut = (BlastMBLookupTable *) lookup_wrap->lut;
        word_length = lut->word_length;
        lut_word_length = lut->lut_word_length;
        masked_locations = lut->masked_locations;
    } 
    else {
        BlastNaLookupTable *lut = 
                             (BlastNaLookupTable *) lookup_wrap->lut;
        word_length = lut->word_length;
        lut_word_length = lut->lut_word_length;
        masked_locations = lut->masked_locations;
    }
    extra_bases = word_length - lut_word_length;

    /* We trust that the bases of the hit itself are exact matches, 
       and look only for exact matches before and after the hit.

       Most of the time, the lookup table width is close to the word size 
       so only a few bases need examining. Also, most of the time (for
       random sequences) extensions will fail almost immediately (the
       first base examined will not match about 3/4 of the time). Thus it 
       is critical to reduce as much as possible all work that is not the 
       actual examination of sequence data */

    for (index = 0; index < num_hits; ++index) {
        Uint4 s_offset = offset_pairs[index].qs_offsets.s_off;
        Uint4 q_offset = offset_pairs[index].qs_offsets.q_off;

        /* begin with the left extension; the initialization is slightly
           faster. q below points to the first base of the lookup table hit
           and s points to the first four bases of the hit (which is
           guaranteed to be aligned on a byte boundary) */

        Int4 context = BSearchContextInfo(q_offset, query_info);
        Int4 context_i = query_info->contexts[context].query_offset;
        Int4 context_f = context_i + query_info->contexts[context].query_length;
        max_bases_left = MIN(extra_bases, MIN(q_offset - context_i, s_offset - range[0]));
        q = q_start + q_offset;
        s = s_start + s_offset / COMPRESSION_RATIO;

        for (extended_left = 0; extended_left < max_bases_left;
             s--, q -= 4, extended_left++) {
            Uint1 byte = s[-1];

            if ((byte & 3) != q[-1] || ++extended_left == max_bases_left)
                break;
            if (((byte >> 2) & 3) != q[-2]
                || ++extended_left == max_bases_left)
                break;
            if (((byte >> 4) & 3) != q[-3]
                || ++extended_left == max_bases_left)
                break;
            if ((byte >> 6) != q[-4])
                break;
        }

        /* do the right extension if the left extension did not find all the
           bases required. Begin at the first base past the lookup table hit 
           and move forwards */

        s_off = s_offset + lut_word_length;
        extended_right = 0;

        if (extended_left < extra_bases) {
            q_off = q_offset + lut_word_length;
            q = q_start + q_off;
            s = s_start + s_off / COMPRESSION_RATIO;
            max_bases_right = MIN(extra_bases - extended_left,
                                  MIN(context_f - q_off,
                                      range[1]  - s_off));

            for (; extended_right < max_bases_right;
                 s++, q += 4, extended_right++) {
                Uint1 byte = s[0];

                if ((byte >> 6) != q[0]
                    || ++extended_right == max_bases_right)
                    break;
                if (((byte >> 4) & 3) != q[1]
                    || ++extended_right == max_bases_right)
                    break;
                if (((byte >> 2) & 3) != q[2]
                    || ++extended_right == max_bases_right)
                    break;
                if ((byte & 3) != q[3])
                    break;
            }

            /* check if enough extra matches were found */

            if (extended_left + extended_right < extra_bases)
                continue;
        }

        q_off = q_offset - extended_left;
        s_off = s_offset - extended_left;
        max_ext = MIN(context_f - q_off, range[1] - s_off) - word_length;
        /* Check that the exact match was not masked out. */
        if (word_length > lut_word_length && 
            s_IsWordSeedMasked(masked_locations, query, subject, q_off,
                               word_length, s_off + word_length, max_ext)) {
            continue; 
        }

        /* check the diagonal on which the hit lies. The boundaries extend
           from the first match of the hit to one beyond the last match */

        if (word_params->container_type == eDiagHash) {
            hits_extended += s_BlastnDiagHashExtendInitialHit(query, subject, 
                                               TRUE,
                                               word_params, matrix, query_info,
                                               ewp->hash_table,
                                               q_off,
                                               s_off,
                                               word_length,
                                               init_hitlist);
        } else {
            hits_extended += s_BlastnDiagTableExtendInitialHit(query, subject, 
                                               TRUE,
                                               word_params, matrix,
                                               query_info,
                                               ewp->diag_table,
                                               q_off,
                                               s_off,
                                               word_length,
                                               init_hitlist);
        }
    }
    return hits_extended;
}

/** Entry i of this list gives the number of pairs of
 * bits that are zero in the bit pattern of i, looking 
 * from right to left
 */
static const Uint1 s_ExactMatchExtendLeft[256] = {
4, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
3, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
3, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
3, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
};

/** Entry i of this list gives the number of pairs of
 * bits that are zero in the bit pattern of i, looking 
 * from left to right
 */
static const Uint1 s_ExactMatchExtendRight[256] = {
4, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
};

/** Perform exact match extensions on the hits retrieved from
 * small-query lookup tables. Assumes the number of bases in a lookup
 * table word and the start offset of each hit is a multiple
 * of 4, and also that the word size is within 4 bases of the lookup
 * table width
 * @param offset_pairs Array of query and subject offsets. [in]
 * @param num_hits Size of the above arrays [in]
 * @param word_params Parameters for word extension [in]
 * @param lookup_wrap Lookup table wrapper structure [in]
 * @param query Query sequence data [in]
 * @param subject Subject sequence data [in]
 * @param matrix Scoring matrix for ungapped extension [in]
 * @param query_info Structure containing query context ranges [in]
 * @param ewp Word extension structure containing information about the 
 *            extent of already processed hits on each diagonal [in]
 * @param init_hitlist Structure to keep the extended hits. 
 *                     Must be allocated outside of this function [in] [out]
 * @param range The subject range
 * @return Number of hits extended. 
 */
static Int4
s_BlastSmallNaExtendAlignedOneByte(const BlastOffsetPair * offset_pairs, 
                       Int4 num_hits,
                       const BlastInitialWordParameters * word_params,
                       LookupTableWrap * lookup_wrap,
                       BLAST_SequenceBlk * query, BLAST_SequenceBlk * subject,
                       Int4 ** matrix, BlastQueryInfo * query_info,
                       Blast_ExtendWord * ewp,
                       BlastInitHitList * init_hitlist,
                       Int4 * range)
{
    BlastSmallNaLookupTable *lut = (BlastSmallNaLookupTable *) lookup_wrap->lut;
    Int4 i;
    Int4 subject_length = subject->length;
    Uint1 *q = query->compressed_nuc_seq;
    Uint1 *s = subject->sequence;
    Int4 lut_word_length = lut->lut_word_length;
    Int4 word_length = lut->word_length;
    Int4 extra_bases = word_length - lut_word_length;
    Int4 hits_extended = 0;
    Int4 max_ext;

    for (i = 0; i < num_hits; i++) {
        Int4 s_off = offset_pairs[i].qs_offsets.s_off;
        Int4 q_off = offset_pairs[i].qs_offsets.q_off;
        Int4 extended_left = 0;
        Int4 extended_right = 0;

        /* the seed is assumed to start on a multiple of 4 bases
           in the subject sequence. Look for up to 4 exact matches
           to the left. The index into q[] below can
           technically be negative, but the compressed version
           of the query has extra pad bytes before q[0] */

        Int4 context = BSearchContextInfo(q_off, query_info);
        Int4 context_i = query_info->contexts[context].query_offset;
        Int4 context_f = context_i + query_info->contexts[context].query_length;
        if (q_off > context_i  && s_off > 0) {
            Uint1 q_byte = q[q_off - 4];
            Uint1 s_byte = s[s_off / COMPRESSION_RATIO - 1];
            extended_left = s_ExactMatchExtendLeft[q_byte ^ s_byte];
            extended_left = MIN(extended_left, extra_bases);
            extended_left = MIN(extended_left, q_off - context_i);
        }

        /* look for up to 4 exact matches to the right of the seed */

        if (extended_left < extra_bases) {
            Int4 q_end = q_off + lut_word_length;
            Int4 s_end = s_off + lut_word_length;
            Uint1 q_byte, s_byte;

            if (q_end == context_f || s_end == subject_length)
                continue;

            q_byte = q[q_end];
            s_byte = s[s_end / COMPRESSION_RATIO];
            extended_right = s_ExactMatchExtendRight[q_byte ^ s_byte];
            extended_right = MIN(extended_right, 
                                 extra_bases - extended_left);
            extended_right = MIN(extended_right, context_f - q_end);
            extended_right = MIN(extended_right, subject_length - s_end);
            if (extended_left + extended_right < extra_bases)
                continue;
        }

        q_off -= extended_left;
        s_off -= extended_left;
        max_ext = MIN(context_f - q_off, subject_length - s_off) - word_length;
        /* Check that the exact match was not masked out. */
        if (word_length > lut_word_length && 
            s_IsWordSeedMasked(lut->masked_locations, query, subject, q_off,
                               word_length, s_off + word_length, max_ext)) {
                continue; 
        }
        
        if (word_params->container_type == eDiagHash) {
            hits_extended += s_BlastnDiagHashExtendInitialHit(query, 
                                           subject, TRUE,  word_params, 
                                           matrix, query_info,
                                           ewp->hash_table,
                                           q_off,
                                           s_off,
                                           word_length,
                                           init_hitlist);
        }
        else {
            hits_extended += s_BlastnDiagTableExtendInitialHit(query, 
                                           subject, TRUE,  word_params, 
                                           matrix, query_info,
                                           ewp->diag_table,
                                           q_off,
                                           s_off,
                                           word_length,
                                           init_hitlist);
        }
    }
    return hits_extended;
}


/** Perform exact match extensions on the hits retrieved from
 * small-query blastn lookup tables, assuming an arbitrary number of bases 
 * in a lookup and arbitrary start offset of each hit. Also 
 * update the diagonal structure.
 * @param offset_pairs Array of query and subject offsets [in]
 * @param num_hits Size of the above arrays [in]
 * @param word_params Parameters for word extension [in]
 * @param lookup_wrap Lookup table wrapper structure [in]
 * @param query Query sequence data [in]
 * @param subject Subject sequence data [in]
 * @param matrix Scoring matrix for ungapped extension [in]
 * @param query_info Structure containing query context ranges [in]
 * @param ewp Word extension structure containing information about the 
 *            extent of already processed hits on each diagonal [in]
 * @param init_hitlist Structure to keep the extended hits. 
 *                     Must be allocated outside of this function [in] [out]
 * @param range The subject range [in]
 * @return Number of hits extended. 
 */
static Int4
s_BlastSmallNaExtend(const BlastOffsetPair * offset_pairs, Int4 num_hits,
                     const BlastInitialWordParameters * word_params,
                     LookupTableWrap * lookup_wrap,
                     BLAST_SequenceBlk * query,
                     BLAST_SequenceBlk * subject, Int4 ** matrix,
                     BlastQueryInfo * query_info,
                     Blast_ExtendWord * ewp,
                     BlastInitHitList * init_hitlist,
                     Int4 * range)
{
    BlastSmallNaLookupTable *lut = (BlastSmallNaLookupTable *) lookup_wrap->lut;
    Int4 i;
    Uint1 *q = query->compressed_nuc_seq;
    Uint1 *s = subject->sequence;
    Int4 word_length = lut->word_length; 
    Int4 lut_word_length = lut->lut_word_length; 
    Int4 hits_extended = 0;
    Int4 max_ext;

    for (i = 0; i < num_hits; i++) {
        Int4 s_off = offset_pairs[i].qs_offsets.s_off;
        Int4 q_off = offset_pairs[i].qs_offsets.q_off;

        Int4 context = BSearchContextInfo(q_off, query_info);
        Int4 context_i = query_info->contexts[context].query_offset;
        Int4 context_f = context_i + query_info->contexts[context].query_length;

        Int4 s_start;
        Int4 q_start;
        Int4 extended_left = 0;
        Int4 extended_right = 0;
        Int4 extra_bases = COMPRESSION_RATIO - (s_off % COMPRESSION_RATIO);

        s_off += extra_bases;
        q_off += extra_bases;

        /* Start the extension at the first multiple of 4 bases in
           the subject sequence to the right of the seed.
           Collect exact matches in groups of four, until a
           mismatch is encountered or the expected number of
           matches is found. The index into q[] below can
           technically be negative, but the compressed version
           of the query has extra pad bytes before q[0] */

        extra_bases = MIN(word_length - lut_word_length + extra_bases, 
                          MIN(q_off - context_i, s_off - range[0]));
        s_start = s_off;
        q_start = q_off;
        while (extended_left < extra_bases) {
            Uint1 q_byte = q[q_start - 4];
            Uint1 s_byte = s[s_start / COMPRESSION_RATIO - 1];
            Uint1 bases = s_ExactMatchExtendLeft[q_byte ^ s_byte];
            extended_left += bases;
            if (bases < 4)
                break;
            q_start -= 4;
            s_start -= 4;
        }
        extended_left = MIN(extended_left, extra_bases);

        /* extend to the right. The extension begins at the first
           base not examined by the left extension */

        extra_bases = word_length - extended_left;
        extra_bases = MIN(extra_bases, context_f - q_off);
        extra_bases = MIN(extra_bases, range[1] - s_off);
        s_start = s_off;
        q_start = q_off;
        while (extended_right < extra_bases) {
            Uint1 q_byte = q[q_start];
            Uint1 s_byte = s[s_start / COMPRESSION_RATIO];
            Uint1 bases = s_ExactMatchExtendRight[q_byte ^ s_byte];
            extended_right += bases;
            if (bases < 4)
                break;
            q_start += 4;
            s_start += 4;
        }
        extended_right = MIN(extended_right, extra_bases);

        if (extended_left + extended_right < word_length)
            continue;

        q_off -= extended_left;
        s_off -= extended_left;
        max_ext = MIN(context_f - q_off, range[1] - s_off) - word_length;
        /* Check that the exact match was not masked out. */
        if (word_length > lut_word_length && 
            s_IsWordSeedMasked(lut->masked_locations, query, subject, q_off,
                               word_length, s_off + word_length, max_ext)) {
            continue; 
        }
        

        if (word_params->container_type == eDiagHash) {
            hits_extended += s_BlastnDiagHashExtendInitialHit(query, subject, 
                                               TRUE, word_params, matrix,
                                               query_info, ewp->hash_table,
                                               q_off,
                                               s_off,
                                               word_length,
                                               init_hitlist);
        } else {
            hits_extended += s_BlastnDiagTableExtendInitialHit(query, subject, 
                                               TRUE, word_params, matrix,
                                               query_info, ewp->diag_table,
                                               q_off,
                                               s_off,
                                               word_length,
                                               init_hitlist);
        }
    }
    return hits_extended;
}

/* Description in na_ungapped.h */

Int2 BlastNaWordFinder(BLAST_SequenceBlk * subject,
                       BLAST_SequenceBlk * query,
                       BlastQueryInfo * query_info,
                       LookupTableWrap * lookup_wrap,
                       Int4 ** matrix,
                       const BlastInitialWordParameters * word_params,
                       Blast_ExtendWord * ewp,
                       BlastOffsetPair * offset_pairs,
                       Int4 max_hits,
                       BlastInitHitList * init_hitlist,
                       BlastUngappedStats * ungapped_stats)
{
    Int4 hitsfound, total_hits = 0;
    Int4 hits_extended = 0;
    TNaScanSubjectFunction scansub = NULL;
    TNaExtendFunction extend = NULL;
    Int4 scan_range[5];
    Int4 word_length;
    Int4 lut_word_length;

    if (lookup_wrap->lut_type == eSmallNaLookupTable) {
        BlastSmallNaLookupTable *lookup = 
                                (BlastSmallNaLookupTable *) lookup_wrap->lut;
        word_length = lookup->word_length;
        lut_word_length = lookup->lut_word_length;
        scansub = (TNaScanSubjectFunction)lookup->scansub_callback;
        extend = (TNaExtendFunction)lookup->extend_callback;
    }
    else if (lookup_wrap->lut_type == eMBLookupTable) {
        BlastMBLookupTable *lookup = 
                                (BlastMBLookupTable *) lookup_wrap->lut;
        if (lookup->discontiguous) {
            word_length = lookup->template_length;
            lut_word_length = lookup->template_length;
        } else {
            word_length = lookup->word_length;
            lut_word_length = lookup->lut_word_length;
        }
        scansub = (TNaScanSubjectFunction)lookup->scansub_callback;
        extend = (TNaExtendFunction)lookup->extend_callback;
    }
    else {
        BlastNaLookupTable *lookup = 
                                (BlastNaLookupTable *) lookup_wrap->lut;
        word_length = lookup->word_length;
        lut_word_length = lookup->lut_word_length;
        scansub = (TNaScanSubjectFunction)lookup->scansub_callback;
        extend = (TNaExtendFunction)lookup->extend_callback;
    }

    scan_range[0] = 0;  /* subject seq mask index */
    scan_range[1] = 0;	/* start pos of scan */
    scan_range[2] = subject->length - lut_word_length; /*end pos of scan*/
    scan_range[3] = 0;  /* current subject seq range start */
    scan_range[4] = subject->length;  /* current subject seq range end */

    /* if sequence is masked, fall back to generic scanner and extender */
    if (subject->num_seq_ranges) {
        scansub = (TNaScanSubjectFunction) 
                  BlastChooseNucleotideScanSubjectAny(lookup_wrap);
        if (extend != s_BlastNaExtendDirect) {
             extend = (lookup_wrap->lut_type == eSmallNaLookupTable) 
                    ? s_BlastSmallNaExtend
                    : s_BlastNaExtend;
        }
        scan_range[1] = word_length - lut_word_length;
        scan_range[2] = subject->length - word_length; 
    }

    ASSERT(scansub);
    ASSERT(extend);

    while(s_DetermineScanningOffsets(subject, word_length, lut_word_length, scan_range)) {
        Int4 next_start;

        /* Pass the last word ending offset */
        next_start = scan_range[2];

        hitsfound = scansub(lookup_wrap, subject, scan_range[1], 
                            offset_pairs, max_hits, &next_start);
        scan_range[1] = next_start;
        if (hitsfound == 0)
            continue;

        total_hits += hitsfound;
        hits_extended += extend(offset_pairs, hitsfound, word_params,
                                lookup_wrap, query, subject, matrix, 
                                query_info, ewp, init_hitlist, &scan_range[3]);
    }

    Blast_ExtendWordExit(ewp, subject->length);

    Blast_UngappedStatsUpdate(ungapped_stats, total_hits, hits_extended,
                              init_hitlist->total);

    if (word_params->ungapped_extension)
        Blast_InitHitListSortByScore(init_hitlist);

    return 0;
}

Int2 MB_IndexedWordFinder( 
        BLAST_SequenceBlk * subject,
        BLAST_SequenceBlk * query,
        BlastQueryInfo * query_info,
        LookupTableWrap * lookup_wrap,
        Int4 ** matrix,
        const BlastInitialWordParameters * word_params,
        Blast_ExtendWord * ewp,
        BlastOffsetPair * offset_pairs,
        Int4 max_hits,
        BlastInitHitList * init_hitlist,
        BlastUngappedStats * ungapped_stats)
{ 
    BlastInitHSP * hsp, * new_hsp, * hsp_end;
    BlastUngappedData dummy_ungapped_data;
    BlastUngappedData * ungapped_data = 0;
    ir_diag_hash * hash = 0;
    ir_hash_entry * e = 0;
    Uint4 word_size;
    Uint4 q_off, s_off;
    Uint4 diag, key;
    Int4 oid = subject->oid;
    Int4 chunk = subject->chunk;
    Int4 context;
    BlastUngappedCutoffs *cutoffs;
    T_MB_IdbGetResults get_results = 
                        (T_MB_IdbGetResults)lookup_wrap->read_indexed_db;
    ASSERT(get_results);
    word_size = get_results(lookup_wrap->lut, oid, chunk, init_hitlist);

    if( word_size > 0 && word_params->ungapped_extension ) {
        hash = ir_hash_create();
        new_hsp = hsp = init_hitlist->init_hsp_array;
        hsp_end = hsp + init_hitlist->total;

        for( ; hsp < hsp_end; ++hsp ) {
            q_off = hsp->offsets.qs_offsets.q_off;
            s_off = hsp->offsets.qs_offsets.s_off;
            diag = IR_DIAG( q_off, s_off );
            key  = IR_KEY( diag );
            e = IR_LOCATE( hash, diag, key );
            if( e != 0 ) {
                if( q_off + word_size - 1 > e->diag_data.qend ) {
                    context = BSearchContextInfo(q_off, query_info);
                    cutoffs = word_params->cutoffs + context;
                    s_NuclUngappedExtend( 
                            query, subject, matrix, 
                            q_off, s_off + word_size, s_off,
                            -(cutoffs->x_dropoff), &dummy_ungapped_data,
                            word_params->nucl_score_table,
                            cutoffs->reduced_nucl_cutoff_score);

                    if( dummy_ungapped_data.score >= cutoffs->cutoff_score ) {
                        ungapped_data = 
                            (BlastUngappedData *)malloc(sizeof(BlastUngappedData));
                        *ungapped_data = dummy_ungapped_data;
                        if( new_hsp != hsp ) *new_hsp = *hsp;
                        new_hsp->ungapped_data = ungapped_data;
                        ++new_hsp;
                    }

                    if( e->diag_data.diag != diag ) e->diag_data.diag = diag;
                    e->diag_data.qend = dummy_ungapped_data.q_start + dummy_ungapped_data.length - 1;
                }
            }
            else {
                if( new_hsp != hsp ) *new_hsp = *hsp;
                ++new_hsp;
            }
        }

        init_hitlist->total = new_hsp - init_hitlist->init_hsp_array;
        hash = ir_hash_destroy( hash );
    }

    if (word_params->ungapped_extension)
        Blast_InitHitListSortByScore(init_hitlist);

    return 0;
}

void BlastChooseNaExtend(LookupTableWrap * lookup_wrap)
{
    if (lookup_wrap->lut_type == eMBLookupTable) {
        BlastMBLookupTable *lut = (BlastMBLookupTable *) lookup_wrap->lut;

        if (lut->lut_word_length == lut->word_length || lut->discontiguous)
            lut->extend_callback = (void *)s_BlastNaExtendDirect;
        else if (lut->lut_word_length % COMPRESSION_RATIO == 0 &&
                 lut->scan_step % COMPRESSION_RATIO == 0)
            lut->extend_callback = (void *)s_BlastNaExtendAligned;
        else
            lut->extend_callback = (void *)s_BlastNaExtend;
    } 
    else if (lookup_wrap->lut_type == eSmallNaLookupTable) {
        BlastSmallNaLookupTable *lut = 
                             (BlastSmallNaLookupTable *) lookup_wrap->lut;

        if (lut->lut_word_length == lut->word_length)
            lut->extend_callback = (void *)s_BlastNaExtendDirect;
        else if (lut->lut_word_length % COMPRESSION_RATIO == 0 &&
                 lut->scan_step % COMPRESSION_RATIO == 0 &&
                 lut->word_length - lut->lut_word_length <= 4)
            lut->extend_callback = (void *)s_BlastSmallNaExtendAlignedOneByte;
        else
            lut->extend_callback = (void *)s_BlastSmallNaExtend;
    }
    else {
        BlastNaLookupTable *lut = (BlastNaLookupTable *) lookup_wrap->lut;

        if (lut->lut_word_length == lut->word_length)
            lut->extend_callback = (void *)s_BlastNaExtendDirect;
        else if (lut->lut_word_length % COMPRESSION_RATIO == 0 &&
                 lut->scan_step % COMPRESSION_RATIO == 0)
            lut->extend_callback = (void *)s_BlastNaExtendAligned;
        else
            lut->extend_callback = (void *)s_BlastNaExtend;
    }
}
