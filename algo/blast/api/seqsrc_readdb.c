/*  $Id: seqsrc_readdb.c,v 1.35 2004/10/06 14:59:28 dondosha Exp $
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
* Author:  Christiam Camacho
*
* File Description:
*   Implementation of the BlastSeqSrc interface using readdb
*
*/

static char const rcsid[] = "$Id: seqsrc_readdb.c,v 1.35 2004/10/06 14:59:28 dondosha Exp $";

#include <algo/blast/api/seqsrc_readdb.h>
#include <algo/blast/core/blast_def.h>
#include <algo/blast/core/blast_util.h>

/** Retrieves the length of the longest sequence in the BlastSeqSrc.
 * @param readdb_handle Pointer to initialized ReadDBFILEPtr structure [in]
 * @param ignoreme Unused by this implementation [in]
 */
static Int4 ReaddbGetMaxLength(void* readdb_handle, void* ignoreme)
{
    ReadDBFILEPtr rdfp = (ReadDBFILEPtr) readdb_handle;
    Int4 retval = 0;

    for (; rdfp; rdfp = rdfp->next)
        retval = MAX(retval, readdb_get_maxlen(rdfp));

    return retval;
}

/** Retrieves the number of sequences in the BlastSeqSrc.
 * @param readdb_handle Pointer to initialized ReadDBFILEPtr structure [in]
 * @param ignoreme Unused by this implementation [in]
 */
static Int4 ReaddbGetNumSeqs(void* readdb_handle, void* ignoreme)
{
   ReadDBFILEPtr rdfp = (ReadDBFILEPtr) readdb_handle;
   Int4 dbnseqs = 0;
   Int8 dblength = 0;
   
   readdb_get_totals_ex(rdfp, &dblength, &dbnseqs, TRUE);
   return dbnseqs;
}

/** Retrieves the total length of all sequences in the BlastSeqSrc.
 * @param readdb_handle Pointer to initialized ReadDBFILEPtr structure [in]
 * @param ignoreme Unused by this implementation [in]
 */
static Int8 ReaddbGetTotLen(void* readdb_handle, void* ignoreme)
{
    ReadDBFILEPtr rdfp = (ReadDBFILEPtr) readdb_handle;
    Int4 dbnseqs = 0;
    Int8 dblength = 0;

    readdb_get_totals_ex(rdfp, &dblength, &dbnseqs, TRUE);
    return dblength;
}

/** Retrieves the average length of sequences in the BlastSeqSrc.
 * @param readdb_handle Pointer to initialized ReadDBFILEPtr structure [in]
 * @param ignoreme Unused by this implementation [in]
 */
static Int4 ReaddbGetAvgLength(void* readdb_handle, void* ignoreme)
{
   Int8 total_length = ReaddbGetTotLen(readdb_handle, ignoreme);
   Int4 num_seqs = MAX(1, ReaddbGetNumSeqs(readdb_handle, ignoreme));

   return (Int4) (total_length/num_seqs);
}

/** Retrieves the name of the BLAST database.
 * @param readdb_handle Pointer to initialized ReadDBFILEPtr structure [in]
 * @param ignoreme Unused by this implementation [in]
 */
static const char* ReaddbGetName(void* readdb_handle, void* ignoreme)
{
    ReadDBFILEPtr rdfp = (ReadDBFILEPtr) readdb_handle;

    return readdb_get_filename(rdfp);
}

/** Retrieves the date of the BLAST database.
 * @param readdb_handle Pointer to initialized ReadDBFILEPtr structure [in]
 * @param ignoreme Unused by this implementation [in]
 */
static Boolean ReaddbGetIsProt(void* readdb_handle, void* ignoreme)
{
    ReadDBFILEPtr rdfp = (ReadDBFILEPtr) readdb_handle;

    return readdb_is_prot(rdfp);
}

/** Retrieves the sequence meeting the criteria defined by its second argument.
 * @param readdb_handle Pointer to initialized ReadDBFILEPtr structure [in]
 * @param args Pointer to GetSeqArg structure [in]
 * @return return codes defined in blast_seqsrc.h
 */
static Int2 ReaddbGetSequence(void* readdb_handle, void* args)
{
    ReadDBFILEPtr rdfp = (ReadDBFILEPtr) readdb_handle;
    GetSeqArg* readdb_args = (GetSeqArg*) args;
    Int4 oid = -1, len = 0, buflen = 0;
    Uint1 *buf = NULL, encoding;
    Boolean has_sentinel_byte;
    Boolean buffer_allocated;

    if (!rdfp || !readdb_args)
        return BLAST_SEQSRC_ERROR;

    oid = readdb_args->oid;
    encoding = readdb_args->encoding;
    has_sentinel_byte = (encoding == BLASTNA_ENCODING);
    buffer_allocated = 
       (encoding == BLASTNA_ENCODING || encoding == NCBI4NA_ENCODING);

    /* free buffers if necessary */
    if (readdb_args->seq)
        BlastSequenceBlkClean(readdb_args->seq);

    /* TODO: this should be cached somewhere */
    if (oid >= readdb_get_num_entries_total(rdfp))
        return BLAST_SEQSRC_EOF;

    if (!buffer_allocated) 
        len = readdb_get_sequence(rdfp, oid, &buf);
    else
        len = readdb_get_sequence_ex(rdfp, oid, &buf, &buflen, has_sentinel_byte);
       
    if (len <= 0) {
        sfree(buf);
        return BLAST_SEQSRC_ERROR;
    }

    BlastSetUp_SeqBlkNew(buf, len, 0, &readdb_args->seq, buffer_allocated);
    /* If there is no sentinel byte, and buffer is allocated, i.e. this is
       the traceback stage of a translated search, set "sequence" to the same 
       position as "sequence_start". */
    if (buffer_allocated && !has_sentinel_byte)
       readdb_args->seq->sequence = readdb_args->seq->sequence_start;

    readdb_args->seq->oid = oid;

    return BLAST_SEQSRC_SUCCESS;
}

/** Deallocates uncompressed sequence buffer, obtained by ReaddbGetSequence.
 * @param readdb_handle Pointer to initialized ReadDBFILEPtr structure [in]
 * @param args Pointer to GetSeqArg structure [in]
 * @return return codes defined in blast_seqsrc.h
 */
static Int2 ReaddbRetSequence(void* readdb_handle, void* args)
{
    GetSeqArg* readdb_args = (GetSeqArg*) args;

    ASSERT(readdb_args);
    return BlastSequenceBlkClean(readdb_args->seq);
}

/** Retrieve length of a given database sequence.
 * @param readdb_handle Pointer to initialized ReadDBFILEPtr structure [in]
 * @param args Pointer to integer indicating ordinal id [in]
 * @return Length of the database sequence or BLAST_SEQSRC_ERROR.
 */
static Int4 ReaddbGetSeqLen(void* readdb_handle, void* args)
{
    ReadDBFILEPtr rdfp = (ReadDBFILEPtr) readdb_handle;
    Int4* oid = (Int4*) args;

    if (!rdfp || !oid)
       return BLAST_SEQSRC_ERROR;

    return readdb_get_sequence_length(rdfp, *oid);
}

/* There are no error messages saved in the ReadDBFILE structure, so the 
 * following getter function is implemented as always returning NULL.
 * @todo Should more meaningful error reporting be implemented?
 */
static Blast_Message* ReaddbGetError(void* readdb_handle, void* args)
{
   return NULL;
}

/** Mutex for retrieving ordinal id chunks from ReadDB in a multi-threaded
 * search.
 */
static TNlmMutex ReaddbMutex;


/** Retrieve next chunk of ordinal ids from a ReadDBFILE structure, in case
 * it contains an oidlist.
 * NB: this function is not MT-safe: ReaddbMutex must be locked/unlocked around
 * any call to this function.
 * @param rdfp List of ReadDBFILE structures [in]
 * @param itr BLAST sequence source iterator [in]
 * @return Status
 */ 
static Int2
ReadDbGetNextOidListChunk(ReadDBFILEPtr rdfp, BlastSeqSrcIterator* itr,
                          Uint4* last_oid_assigned)
     
{
   Int2 status = BLAST_SEQSRC_SUCCESS;
   OIDListPtr oidlist;
   Uint4  gi_start, gi_end;
   Uint4* id_list;
   Uint4 oidindex  = 0;

   if (!itr || !last_oid_assigned)
      return BLAST_SEQSRC_ERROR;

   for ( ; rdfp; rdfp = rdfp->next) {
      oidlist = rdfp->oidlist;
   
      /* If there is no OID list, go to the next readdb structure. */
      if (!oidlist)
         continue;

      gi_start = MAX(*last_oid_assigned, (Uint4)rdfp->start) - rdfp->start;
      gi_end = (Uint4)oidlist->total + 1;
      id_list = itr->oid_list;

      if (gi_start < gi_end) {
         Uint4 bit_start = gi_start % MASK_WORD_SIZE;
         Uint4 gi;

         for(gi = gi_start; (gi < gi_end) && (oidindex < itr->chunk_sz);) {
            Int4 bit_end = ((gi_end - gi + bit_start) < MASK_WORD_SIZE) ?
               (gi_end - gi + bit_start) : MASK_WORD_SIZE;
            Int4 bit;
            
            Uint4 mask_index = gi / MASK_WORD_SIZE;
            Uint4 mask_word  = Nlm_SwapUint4(oidlist->list[mask_index]);
            
            if ( mask_word ) {
               for(bit = bit_start; bit<bit_end && oidindex<itr->chunk_sz; bit++) {
                  Uint4 bitshift = (MASK_WORD_SIZE-1)-bit;
                  
                  if ((mask_word >> bitshift) & 1) {
                     id_list[ oidindex++ ] = rdfp->start + (gi - bit_start) + bit;
                  }
               }
               gi += bit - bit_start;
            } else {
               gi += bit_end - bit_start;
            }
            
            bit_start = 0;
         }

         if (oidindex == itr->chunk_sz || !rdfp->next) {
            itr->itr_type = eOidList;
            itr->current_pos = 0;
            *last_oid_assigned = rdfp->start + gi;
            itr->chunk_sz = oidindex;
            break;
         }
      } /* End if (gi_start < gi_end) */
   } /* End loop over ReadDBFILE's */
   
   if (!rdfp) 
      status = BLAST_SEQSRC_EOF;

   return status;
}



static Int2 ReaddbGetNextChunk(void* readdb_handle, BlastSeqSrcIterator* itr)
{
    ReadDBFILEPtr rdfp = (ReadDBFILEPtr) readdb_handle;
    ReadDBFILEPtr rdfp_head = rdfp;
    unsigned int nseqs;
    Uint4 current_oid;
    Int2 status = BLAST_SEQSRC_SUCCESS;
    Uint4 real_readdb_entries;

    if (!rdfp || !itr)
        return BLAST_SEQSRC_ERROR;

    real_readdb_entries = readdb_get_num_entries_total_real(rdfp);
    
    /* Lock the mutex before retrieving the next chunk */
    NlmMutexLockEx(&ReaddbMutex);
    ASSERT(rdfp->shared_info);

    current_oid = rdfp->shared_info->last_oid_assigned;

    if (current_oid < (unsigned int) rdfp->start)
       current_oid = (unsigned int) rdfp->start;
    
    for ( ; rdfp && !rdfp->oidlist; rdfp = rdfp->next) {
       if (rdfp->stop > 0) {
          nseqs = rdfp->stop + 1;
       } else if (rdfp->aliasnseq) {
          nseqs = rdfp->aliasnseq;
       } else {
          nseqs = rdfp->num_seqs;
       }
       
       if (current_oid < nseqs)
          break;
    }

    if (!rdfp) {
       status = BLAST_SEQSRC_EOF;
    } else if (!rdfp->oidlist) {
       itr->itr_type = eOidRange;
       itr->current_pos = itr->oid_range[0] = current_oid;
       itr->oid_range[1] = MIN(current_oid + itr->chunk_sz, nseqs);
       rdfp_head->shared_info->last_oid_assigned = itr->oid_range[1];
    } else {
       status = ReadDbGetNextOidListChunk(rdfp, itr, 
                   &rdfp_head->shared_info->last_oid_assigned);
    }

    NlmMutexUnlock(ReaddbMutex);

    return status;
}

static Int4 ReaddbIteratorNext(void* seqsrc, BlastSeqSrcIterator* itr)
{
    BlastSeqSrc* bssp = (BlastSeqSrc*) seqsrc;
    Int4 retval = BLAST_SEQSRC_EOF;
    Int4 status = BLAST_SEQSRC_SUCCESS;
    Uint4 last_pos = 0;

    ASSERT(bssp);
    ASSERT(itr);

    /* If iterator is uninitialized/invalid, retrieve the next chunk from the
     * BlastSeqSrc */
    if (itr->current_pos == UINT4_MAX) {
        status = BLASTSeqSrcGetNextChunk(bssp, itr);
        if (status != BLAST_SEQSRC_SUCCESS) {
            return status;
        }
    }

    if (itr->itr_type == eOidRange) {
        retval = itr->current_pos;
        last_pos = itr->oid_range[1];
    } else if (itr->itr_type == eOidList) {
        retval = itr->oid_list[itr->current_pos];
        last_pos = itr->chunk_sz;
    } else {
        /* Unsupported/invalid iterator type! */
        fprintf(stderr, "Invalid iterator type: %d\n", itr->itr_type);
        retval = BLAST_SEQSRC_ERROR;
    }

    ++itr->current_pos;
    if (itr->current_pos >= last_pos) {
        itr->current_pos = UINT4_MAX;  /* invalidate internal iteration */
    }

    return retval;
}

BlastSeqSrc* ReaddbSeqSrcNew(BlastSeqSrc* retval, void* args)
{
    ReaddbNewArgs* rargs = (ReaddbNewArgs*) args;
    ReadDBFILEPtr rdfp = NULL;

    if (!retval)
        return NULL;

    /* Initialize the rdfp */
    if ( !(rdfp = readdb_new(rargs->dbname, rargs->is_protein)))
        return NULL;

    /* Initialize the BlastSeqSrc structure fields with user-defined function
     * pointers and rdfp */
    SetDeleteFnPtr(retval, &ReaddbSeqSrcFree);
    SetCopyFnPtr(retval, &ReaddbSeqSrcCopy);
    SetDataStructure(retval, (void*) rdfp);
    SetGetNumSeqs(retval, &ReaddbGetNumSeqs);
    SetGetMaxSeqLen(retval, &ReaddbGetMaxLength);
    SetGetAvgSeqLen(retval, &ReaddbGetAvgLength);
    SetGetTotLen(retval, &ReaddbGetTotLen);
    SetGetAvgSeqLen(retval, &ReaddbGetAvgLength);
    SetGetName(retval, &ReaddbGetName);
    SetGetIsProt(retval, &ReaddbGetIsProt);
    SetGetSequence(retval, &ReaddbGetSequence);
    SetGetSeqLen(retval, &ReaddbGetSeqLen);
    SetGetNextChunk(retval, &ReaddbGetNextChunk);
    SetIterNext(retval, &ReaddbIteratorNext);
    SetGetError(retval, &ReaddbGetError);
    SetRetSequence(retval, &ReaddbRetSequence);

    /* Set the range, if it is specified */
    if (rargs->first_db_seq > 0) {
       while (rdfp && rdfp->stop < rargs->first_db_seq) {
          /* Make this rdfp's range empty */
          rdfp->start = rdfp->stop + 1;
          rdfp = rdfp->next;
       }
       rdfp->start = rargs->first_db_seq;
    }
    if (rargs->final_db_seq > 0) {
       while (rdfp && rdfp->stop < rargs->final_db_seq)
          rdfp = rdfp->next;
       /* Set last sequence for this and all subsequent rdfp's to the one
          in the arguments, making the subsequent rdfp's ranges empty. 
          Note that final_db_seq in arguments is 1 beyond the last sequence
          number to search. */
       for ( ; rdfp; rdfp = rdfp->next)
          rdfp->stop = rargs->final_db_seq - 1;
    }

    return retval;
}

BlastSeqSrc* ReaddbSeqSrcFree(BlastSeqSrc* bssp)
{
    if (!bssp) 
        return NULL;
    readdb_destruct((ReadDBFILEPtr)GetDataStructure(bssp));
    sfree(bssp);
    return NULL;
}

BlastSeqSrc* ReaddbSeqSrcCopy(BlastSeqSrc* bssp)
{
   ReadDBFILE* rdfp = NULL;

   if (!bssp) 
      return NULL;

   rdfp = readdb_attach((ReadDBFILEPtr)GetDataStructure(bssp));

   SetDataStructure(bssp, (void*) rdfp);
    
   return bssp;
}

BlastSeqSrc* 
ReaddbBlastSeqSrcInit(const char* dbname, Boolean is_prot, int first_seq, 
                      int last_seq, void* extra_arg)
{
    BlastSeqSrcNewInfo bssn_info;
    BlastSeqSrc* seq_src = NULL;
    ReaddbNewArgs* readdb_args = 
        (ReaddbNewArgs*) calloc(1, sizeof(ReaddbNewArgs));;
    readdb_args->dbname = strdup(dbname);
    readdb_args->is_protein = is_prot;
    readdb_args->first_db_seq = first_seq;
    readdb_args->final_db_seq = last_seq;
    bssn_info.constructor = &ReaddbSeqSrcNew;
    bssn_info.ctor_argument = (void*) readdb_args;

    seq_src = BlastSeqSrcNew(&bssn_info);
    sfree(readdb_args->dbname);
    sfree(readdb_args);
    return seq_src;
}
