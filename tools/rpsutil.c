/* $Id: rpsutil.c,v 6.22 2000/05/08 20:52:20 shavirin Exp $
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
* File Name:  $RCSfile: rpsutil.c,v $
*
* Author:  Sergei Shavirin
*
* Initial Version Creation Date: 12/14/1999
*
* $Revision: 6.22 $
*
* File Description:
*         Reversed PSI BLAST utilities file
*
* $Log: rpsutil.c,v $
* Revision 6.22  2000/05/08 20:52:20  shavirin
* Updated query_gapped_start in case when hsp hits 2 boundaries.
*
* Revision 6.21  2000/05/02 17:57:19  shavirin
* Corrected path to RPS Databases changed definition of RPSInit() function.
*
* Revision 6.20  2000/04/13 18:49:59  shavirin
* Fixed serious memory leaks.
*
* Revision 6.19  2000/04/12 14:27:04  shavirin
* Updated case when query sequence hit more, than one boundary in PSSM
* database.
*
* Revision 6.18  2000/03/28 20:32:27  shavirin
* Added functions RPSInfoDetach and RPSInfoAttach.
*
* Revision 6.17  2000/02/28 14:21:21  madden
* Fix bugs with reuse of index and NULLB for copyMatrix
*
* Revision 6.16  2000/02/25 17:00:58  madden
* Added RPSimpalaStatCorrections to correct statistics in the IMPALA manner
*
* Revision 6.15  2000/02/25 14:43:26  shavirin
* Changes, that fixed 64bit Solaris sparc and partually fixed DEC Alpha.
*
* Revision 6.14  2000/02/23 22:19:19  shavirin
* Fixed bug with optimized version in RPSAlignTraceBack().
*
* Revision 6.13  2000/02/23 21:51:08  shavirin
* Empty message.
*
* Revision 6.12  2000/02/23 21:02:50  shavirin
* Fixed -z and -Y options in rpsblast.
*
* Revision 6.11  2000/02/17 19:18:23  shavirin
* Removed reference to theCasheSize.
*
* Revision 6.10  2000/02/17 18:36:44  shavirin
* Added conversion of protein filter SeqLoc into DNA in case of
* translated RPS Blast.
*
* Revision 6.9  2000/02/15 16:15:51  shavirin
* Adjusted effective search space.
*
* Revision 6.8  2000/02/11 22:10:51  shavirin
* Checked seqalign for NULL, removed some memory leaks.
*
* Revision 6.7  2000/02/11 20:44:07  shavirin
* Added possibility to search PSSM database against DNA sequence.
*
* Revision 6.6  2000/02/08 17:36:42  shavirin
* Fixed problem with extentions, that cross boundaries.
*
* Revision 6.5  2000/02/03 17:45:08  shavirin
* Fixed problem with filtering.
*
* Revision 6.4  2000/01/14 18:32:44  shavirin
* Added adjustment of ExtendWord* structure.
*
* Revision 6.3  2000/01/07 22:35:02  shavirin
* Major changes - to use single lookup table for complete database set.
*
* Revision 6.2  1999/12/30 18:35:20  shavirin
* Adapted to the new format of big matrix file.
*
* Revision 6.1  1999/12/29 19:39:47  shavirin
* Initial revision.
*
*
* ==========================================================================
*/

#include <rpsutil.h>

#define RPS_MAGIC_NUMBER 7702

typedef struct _RPSSap {
    SeqAlignPtr sap;
    Nlm_FloatHi e_value;
} RPSap, PNTR RPSapPtr;

typedef struct _RPSapSort {
    RPSapPtr PNTR rpsap;
    Int4 count;
    Int4 allocated;
} RPSapSort, PNTR RPSapSortPtr;

void RPSFreeLookup(RPSLookupPtr lookup)
{
    Nlm_MemMapFini(lookup->mmLookup);
    MemFree(lookup);

    return;
}

RPSLookupPtr RPSInitLookup(CharPtr LookupFile)
{
    RPSLookupPtr lookup;
    
    lookup = MemNew(sizeof(RPSLookup));

    if((lookup->mmLookup = Nlm_MemMapInit(LookupFile)) == NULL) {
        ErrPostEx(SEV_FATAL, 0, 0, "RPSInit: mmap of lookup failed");
        return NULL;
    }
    
    lookup->header = (Int4Ptr) lookup->mmLookup->mmp_begin;
    
    if(lookup->header[0] != RPS_MAGIC_NUMBER) {
        ErrPostEx(SEV_FATAL, 0, 0, "RPSInit: invalid lookup file");
        return NULL;
    }

    lookup->offsets = lookup->header + 8;
    lookup->entries = lookup->header[1];

    lookup->looktbl = lookup->offsets + (lookup->entries + 1);
    
    return lookup;
}
RPSInfoPtr RPSInfoAttach(RPSInfoPtr rpsinfo_old)
{
    RPSInfoPtr rpsinfo_new;
    
    rpsinfo_new = MemNew(sizeof(RPSInfo));
    
    MemCpy(rpsinfo_new, rpsinfo_old, sizeof(RPSInfo));
    
    rpsinfo_new->rdfp = readdb_attach(rpsinfo_old->rdfp);
    
    return rpsinfo_new;
}

void RPSInfoDetach(RPSInfoPtr rpsinfo)
{
    readdb_destruct(rpsinfo->rdfp);
    MemFree(rpsinfo);
    return;
}

RPSInfoPtr RPSInit(CharPtr dbname, Int4 query_is_prot)
{
    RPSInfoPtr rpsinfo;
    Int4Ptr header;
    Char rps_matrix[128], rps_lookup[128];

    rpsinfo = MemNew(sizeof(RPSInfo));
    
    if((rpsinfo->rdfp = readdb_new(dbname, TRUE)) == NULL)
        return NULL;

    rpsinfo->query_is_prot = query_is_prot;

    ReadDBBioseqFetchEnable ("rpsblast", dbname, FALSE, TRUE); 
    

    sprintf(rps_matrix, "%s.rps", rpsinfo->rdfp->full_filename);
    sprintf(rps_lookup, "%s.loo", rpsinfo->rdfp->full_filename);
    
    
    if((rpsinfo->mmMatrix = Nlm_MemMapInit(rps_matrix)) == NULL) {
        ErrPostEx(SEV_FATAL, 0, 0, "RPSInit: mmap of matrix failed");
        return (NULL);
    }

    header = (Int4Ptr) rpsinfo->mmMatrix->mmp_begin;
    
    if(header[0] != RPS_MAGIC_NUMBER) {
        ErrPostEx(SEV_FATAL, 0, 0, "RPSInit: invalid matrix memmap file");
        return (NULL);
    }
    
    rpsinfo->matrixCount = header[1];
    
    rpsinfo->offsets = header + 2; /* Strarting from 3rd integer */

    /* Matrix  is started from this position */
    rpsinfo->bigMatrix = (RPScoreRow *)(header + 2 + rpsinfo->matrixCount + 1);


    /* Now initializing lookup tables */

    rpsinfo->lookup = RPSInitLookup(rps_lookup);

    
    return rpsinfo;
}

void RPSClose(RPSInfoPtr rpsinfo)
{

    readdb_destruct(rpsinfo->rdfp);
    Nlm_MemMapFini(rpsinfo->mmMatrix);

    RPSFreeLookup(rpsinfo->lookup);

    MemFree(rpsinfo);
    
    ReadDBBioseqFetchDisable();
    
    return;
}

static Int4 RPSBinary_Search(Int4 n, Uint4Ptr A, Int4 i)
{
    Int4 m, b, e;

    n++;
    b = 0;
    e = i;

    while (b < e - 1) {
	m = (b + e) / 2;
	if (Nlm_SwapUint4(A[m]) > n)
	    e = m;
	else
	    b = m;
    }

    return b;
}

Int4 RPSFindSequence(ReadDBFILEPtr rdfp, Int4 start, Int4Ptr sequence_start,
                     Int4Ptr seq_length)
{
    Int4 seq_num, num_seqs;
    
    num_seqs =  readdb_get_num_entries(rdfp);
    
    seq_num = RPSBinary_Search(start, rdfp->sequence_index, num_seqs);
    
    *sequence_start = Nlm_SwapUint4(rdfp->sequence_index[seq_num]) - 1;
    *seq_length = Nlm_SwapUint4(rdfp->sequence_index[seq_num+1]) - 
        Nlm_SwapUint4(rdfp->sequence_index[seq_num]) -1;
    
    return seq_num;
}

RPSequencePtr RPSGetSequenceEx(RPSInfoPtr rpsinfo, Int4 seqnum, Boolean all_seqs)
{
    RPSequencePtr rpsp;
    Int4 offset, seqlength, i, row, num_seqs;

    if(rpsinfo == NULL)
        return NULL;

    if(all_seqs) seqnum = 0;
    
    rpsp = MemNew(sizeof(RPSequence));

    readdb_get_descriptor(rpsinfo->rdfp, seqnum, 
                          &rpsp->seqid, &rpsp->description);
    
    rpsp->number = seqnum;
    rpsp->seqlen = readdb_get_sequence (rpsinfo->rdfp, seqnum, 
                                        &rpsp->sequence);
    
    if(all_seqs) {
        /* Length of sequence updated to all length */
        num_seqs = readdb_get_num_entries(rpsinfo->rdfp);   
        rpsp->seqlen = Nlm_SwapUint4(rpsinfo->rdfp->sequence_index[num_seqs]) - Nlm_SwapUint4(rpsinfo->rdfp->sequence_index[0]) - 1;
    }
    
    /* Now reading the posMatrix */
    
    if(all_seqs) {
        offset = 0;
        seqlength = rpsp->seqlen;
    } else {
        offset = rpsinfo->offsets[seqnum];
        seqlength = rpsinfo->offsets[seqnum + 1] - rpsinfo->offsets[seqnum] - 1;
    }
    
    rpsp->posMatrix = MemNew((seqlength + 1) * sizeof(BLAST_Score *));
    for (i = 0, row = offset; i < seqlength+1; i++, row++) {
        rpsp->posMatrix[i] = (BLAST_Score *) &(rpsinfo->bigMatrix[row][0]);
    }

    if(all_seqs) {        
        /* This is pointer to the start of the lookup table */
        rpsp->mod_lt = (ModLAEntry *) (rpsinfo->lookup->mmLookup->mmp_begin + rpsinfo->lookup->offsets[seqnum]); /* Only one lookup table */
        
        rpsp->mod_lookup_table_memory = 
            (ModLookupPositionPtr) (rpsp->mod_lt + 1 + RPS_ARRAY_SIZE);
        
        rpsp->num_pos_added = rpsinfo->lookup->header[2];
        rpsp->num_unique_pos_added = rpsinfo->lookup->header[3];
        rpsp->mod_lookup_table_size = rpsinfo->lookup->header[4];
    } 

    return rpsp;
}

RPSequencePtr RPSGetSequence(RPSInfoPtr rpsinfo, Int4 seqnum)
{
    return RPSGetSequenceEx(rpsinfo, seqnum, FALSE);
}

RPSequencePtr RPSGetBIGSequence(RPSInfoPtr rpsinfo, BioseqPtr PNTR bsp_out)
{
    RPSequencePtr rpsp;
    ValNodePtr vnp;
    BioseqPtr bsp;

    rpsp = RPSGetSequenceEx(rpsinfo, 0, TRUE);
    
    /* Creating BIG Bioseq - for testing only */
#if 0
    bsp = BioseqNew();
    
    bsp->mol = Seq_mol_aa;
    bsp->seq_data_type = Seq_code_ncbistdaa;
    bsp->repr = Seq_repr_raw;
    bsp->length = rpsp->seqlen;
    
    bsp->seq_data = BSNew(1024);
    BSWrite(bsp->seq_data, rpsp->sequence, bsp->length);
    
    /* BIG SeqId */
    /* bsp->id = SeqIdDup(rpsp->seqid); */

    bsp->id = MakeNewProteinSeqId (NULL, NULL);
    
    SeqIdSetFree(rpsp->seqid);
    rpsp->seqid = SeqIdDup(bsp->id);

    /* BIG Description */
    vnp = SeqDescrNew(NULL);
    vnp->choice = Seq_descr_title;
    vnp->data.ptrvalue = StringSave (rpsp->description);
    bsp->descr = vnp;    
    
    SeqMgrAddToBioseqIndex (bsp);
    
    *bsp_out = bsp;
#endif
    return rpsp;
}

void RPSequenceFree(RPSequencePtr rpseqp)
{
    
    MemFree(rpseqp->posMatrix);
    SeqIdSetFree(rpseqp->seqid);
    MemFree(rpseqp->description);
    
    MemFree(rpseqp);
    
    return;
}

void RPSUpdateDbSize(BLAST_OptionsBlkPtr options, RPSInfoPtr rpsinfo, Int4 query_length)
{
    Int4 db_length;

    options->dbseq_num = rpsinfo->matrixCount;

    if(options->db_length == 0) {
        options->db_length = rpsinfo->offsets[rpsinfo->matrixCount];
    }

    if(options->searchsp_eff == 0) {
        if(rpsinfo->query_is_prot)
            options->searchsp_eff = query_length * options->db_length;
        else
            options->searchsp_eff = query_length * 6 * options->db_length;
    }
    
    return;
}

/* ------ Aux. functions to sort alignments --------- */
#define SAP_SORT_CHUNK 512
static RPSapSortPtr RPSapSortInit(void)
{
    RPSapSortPtr ssp;
    
    ssp = MemNew(sizeof(RPSapSort));
    
    ssp->allocated = SAP_SORT_CHUNK;
    ssp->rpsap = MemNew(sizeof(RPSapPtr) * ssp->allocated);
    ssp->count = 0;
    return ssp;
}

static void RPSAddSap(RPSapSortPtr ssp, SeqAlignPtr sap)
{
    RPSapPtr rpsp;
    ScorePtr thisScorePtr;

    if(sap == NULL)
        return;
    
    if(ssp->count >= ssp->allocated) {
        ssp->allocated += SAP_SORT_CHUNK;
        ssp->rpsap = Realloc(ssp->rpsap, sizeof(Int4Ptr) * ssp->allocated);
    }
    
    rpsp = MemNew(sizeof(RPSap));
    
    rpsp->sap = sap;
        
    /* Extracting e_value from SeqAlignPtr */
    thisScorePtr = sap->score;
    while ((thisScorePtr != NULL) &&
           (StringICmp(thisScorePtr->id->str, "e_value") != 0) &&
           (StringICmp(thisScorePtr->id->str, "sum_e") != 0)) {
        thisScorePtr = thisScorePtr->next;
    }

    if(NULL == thisScorePtr)
        rpsp->e_value = 10.0;
    else
        rpsp->e_value = (Nlm_FloatHi) (thisScorePtr->value.realvalue);

    ssp->rpsap[ssp->count] = rpsp;
    ssp->count++;

    return;
}
static int LIBCALLBACK RPSortCallback(VoidPtr i, VoidPtr j)
{
    RPSapPtr rpsp1, rpsp2;
    
    rpsp1 =  *((RPSapPtr *) i);
    rpsp2 =  *((RPSapPtr *) j);
    
    if (rpsp1->e_value > rpsp2->e_value)
        return (1);
    
    if (rpsp1->e_value < rpsp2->e_value)
        return (-1);
    
    return (0);    
}

static SeqAlignPtr RPSReadSapSort(RPSapSortPtr ssp)
{
    SeqAlignPtr head = NULL, seqalign_var, sap_last = NULL;
    Int4 i;

    HeapSort(ssp->rpsap, ssp->count, sizeof(RPSapPtr), RPSortCallback);
    
    for(i = 0; i < ssp->count; i++) {
        
        if (head == NULL) {
            head = ssp->rpsap[i]->sap;
        } else {
            for (seqalign_var = sap_last; seqalign_var->next != NULL;) {
                seqalign_var = seqalign_var->next;
            }
            seqalign_var->next = ssp->rpsap[i]->sap;
        }
        sap_last = ssp->rpsap[i]->sap;
    }
    
    return head;
}

static void RPSapSortFree(RPSapSortPtr ssp)
{
    Int4 i;

    for(i = 0; i < ssp->count; i++) {
        MemFree(ssp->rpsap[i]);
    }
    MemFree(ssp->rpsap);
    
    MemFree(ssp);
    return;
}

Boolean RPSReturnQuery(BlastSearchBlkPtr search, BioseqPtr query_bsp,
                       Uint1Ptr query_seq_start)
{
    SeqPortPtr spp;
    Int4 index, query_length;
    /*    Uint1Ptr query_seq, query_seq_start;
          Uint1 residue; */
  
    if(search == NULL || query_bsp == NULL)
        return FALSE;

    SeqIdSetFree(search->query_id);
    search->query_id = SeqIdDup(query_bsp->id);

    search->query_slp = SeqLocFree(search->query_slp);
    query_length = query_bsp->length;
    
    ValNodeAddPointer(&search->query_slp, SEQLOC_WHOLE, 
                      SeqIdDup(SeqIdFindBest(query_bsp->id, SEQID_GI)));
    
    BlastSequenceAddSequence(search->context[0].query, NULL, query_seq_start, 
                             query_length, query_length, 0);
    
    search->context[0].query_allocated = TRUE;    
    
    search->positionBased = FALSE; /* Now query is simple sequence back */
    search->sbp->posMatrix = NULL;
    
    return TRUE;
}

Boolean RPSubstituteQueryLookup(BlastSearchBlkPtr search, 
                                RPSequencePtr rpseq, Boolean update_lookup)
{
    Int4 hitlist_max;
    SeqLocPtr slp;
    LookupTablePtr lookup;
    
    if(search == NULL || rpseq == NULL)
        return FALSE;
    
    SeqIdSetFree(search->query_id);
    search->query_id = SeqIdDup(rpseq->seqid);
    
    if(search->allocated & BLAST_SEARCH_ALLOC_QUERY_SLP)
        search->query_slp = SeqLocFree(search->query_slp);
    else 
        search->allocated += BLAST_SEARCH_ALLOC_QUERY_SLP;

    ValNodeAddPointer(&search->query_slp, SEQLOC_WHOLE, 
                      SeqIdDup(SeqIdFindBest(rpseq->seqid, SEQID_GI)));

    search->sbp->posMatrix = rpseq->posMatrix;
    search->positionBased = TRUE;
    search->sbp->kbp = search->sbp->kbp_psi;
    search->sbp->kbp_gap = search->sbp->kbp_gap_psi;
    hitlist_max = search->result_struct->hitlist_max;
    
    BlastSequenceAddSequence(search->context[0].query, NULL, 
                             rpseq->sequence-1, rpseq->seqlen, 
                             rpseq->seqlen, 0);
    search->context[0].query->effective_length = rpseq->seqlen;

    /* search->context[0].query_allocated = FALSE; */
    
    if(update_lookup) {

        search->result_struct = 
            BLASTResultsStructDelete(search->result_struct);
        
        if(search->result_struct == NULL) {
            search->result_struct = BLASTResultsStructNew(hitlist_max, search->pbp->max_pieces, search->pbp->hsp_range_max);
        }

        if (!(search->allocated & BLAST_SEARCH_ALLOC_WFP_FIRST)) {
            search->wfp = search->wfp_first = search->wfp_second = 
                BLAST_WordFinderNew(search->sbp->alphabet_size, 
                                    RPS_WORD_SIZE, 1, FALSE);
        }
        
        lookup = search->wfp->lookup;
        
        lookup->mod_lt = rpseq->mod_lt;
        lookup->mod_lookup_table_memory = rpseq->mod_lookup_table_memory;
        
        lookup->num_pos_added = rpseq->num_pos_added;
        lookup->num_unique_pos_added = rpseq->num_unique_pos_added;
        lookup->mod_lookup_table_size = rpseq->mod_lookup_table_size;

        MemFree(search->ewp_params);
        search->ewp_params = BLAST_ExtendWordParamsNew(rpseq->seqlen, (Boolean) search->pbp->window_size != 0, search->pbp->window_size);
        
        /* Only one context used in this program */
        BLAST_ExtendWordDestruct(search->context[0].ewp);
        search->context[0].ewp = BLAST_ExtendWordNew(search->ewp_params);

        /* Adjusting database length, that should be set to query length */
        /* search->dblen = rpseq->seqlen; */

    }
    
    return TRUE;
}

void RPSLookupCleanUp(LookupTablePtr lookup)
{
    MemFree(lookup->pv_array);
    lookup->pv_array = NULL;

    MemFree(lookup->mod_lt);
    lookup->mod_lt = NULL;

    MemFree(lookup->mod_lookup_table_memory);
    lookup->mod_lookup_table_memory = NULL;
    
    /* lookup_deallocate_memory(lookup); */
    /* lookup->mem_struct       = NULL;
       lookup->mem_struct_start = NULL; */

    lookup->num_pos_added = 0;
    lookup->num_unique_pos_added = 0;
    lookup->mod_lookup_table_size = 0;

    return;
}
void RPSExchangeInt(Int4 *a, Int4 *b)
{
    Int4 value;
    
    value = *a;
    *a = *b;
    *b = value;

    return;
}

Boolean RPSUpdateCoordinates(ReadDBFILEPtr rdfp, BLASTResultHitlistPtr result,
                             Boolean reverse)
{
    Int4 hspcnt, index, index2;
    BLASTResultHspPtr hsp_array;
    Int4 seq_no1, seq_no2, sequence_start1, seq_length1; 
    Int4 sequence_start2, seq_length2;
    Int4 num_regions, query_length_real, subject_length_real;
    
    hspcnt = result->hspcnt;
    hsp_array = result->hsp_array;


    for (index = 0; index < result->hspcnt; index++) {

        seq_no1 = RPSFindSequence(rdfp, hsp_array[index].query_offset, 
                                  &sequence_start1, &seq_length1);
        seq_no2 = RPSFindSequence(rdfp, hsp_array[index].query_offset +
                                  hsp_array[index].query_length, 
                                  &sequence_start2, &seq_length2);

        /* if seq_no1 != seq_no2 we crossed the boundary between sequences
           and therefore we need to expand hsp_array */

        if((num_regions = (seq_no2 - seq_no1)) > 0) {
            
            hspcnt = result->hspcnt; /* Count to be filled */

            /* Adjusting existing hit */
            query_length_real = hsp_array[index].query_length;
            subject_length_real = hsp_array[index].subject_length;

            hsp_array[index].query_length = sequence_start1 + seq_length1 - 
                hsp_array[index].query_offset;
            
            if(subject_length_real > hsp_array[index].query_length + 1) {
                hsp_array[index].subject_length = hsp_array[index].query_length; 
                hsp_array[index].query_gapped_start = hsp_array[index].query_offset;
                /* If more regions exist - reallocating hsp_array */
                
                num_regions = 1; /* Assuming only obe boundary for now */
                
                result->hsp_array = Realloc(hsp_array, (hspcnt+num_regions)*sizeof(BLASTResultHsp));
                hsp_array = result->hsp_array;
                result->hspcnt++;
                
                /* for(index2 = 0; index2 < hspcnt+num_regions; index2++) {
                   hsp_array[index2].point_back = result;
                   } */
                
                /* Last tail of the hit */
                
                MemCpy(&result->hsp_array[hspcnt], 
                       &result->hsp_array[index], sizeof(BLASTResultHsp));
                
                hsp_array[hspcnt].query_offset = sequence_start1 + seq_length1+1;
                
                hsp_array[hspcnt].query_length = query_length_real - hsp_array[index].query_length - 1;
                
                hsp_array[hspcnt].query_gapped_start = hsp_array[hspcnt].query_offset + hsp_array[hspcnt].query_length;
                
                hsp_array[hspcnt].subject_offset = hsp_array[index].subject_offset + hsp_array[index].subject_length;
                
                hsp_array[hspcnt].subject_length = subject_length_real - hsp_array[index].subject_length;
                
                hsp_array[hspcnt].subject_gapped_start = hsp_array[hspcnt].subject_offset + hsp_array[hspcnt].subject_length;
                
            
                hsp_array[hspcnt].number = seq_no1 + 1;

                /* hsp_array[hspcnt].query_offset -= sequence_start2;
                   hsp_array[hspcnt].query_gapped_start -= sequence_start2; */
                
                hspcnt = result->hspcnt;
            } else {
                hsp_array[index].query_length = subject_length_real;
                hsp_array[index].query_gapped_start = hsp_array[index].query_offset + hsp_array[index].query_length;
            }
        }
        
        hsp_array[index].number = seq_no1;
        
        
        /* Assigning coordinates in sequence own coordinates */
        hsp_array[index].query_offset -= sequence_start1;
        hsp_array[index].query_gapped_start -= sequence_start1;
        
        if(reverse) {
            RPSExchangeInt(&hsp_array[index].query_offset, 
                           &hsp_array[index].subject_offset);
            RPSExchangeInt(&hsp_array[index].query_length, 
                           &hsp_array[index].subject_length);
            RPSExchangeInt(&hsp_array[index].query_gapped_start, 
                           &hsp_array[index].subject_gapped_start);
        }
    }
    
    return TRUE;
}

BLASTResultHitlistPtr RPSExtractNewResult(BLASTResultHitlistPtr result)
{
    BLASTResultHitlistPtr new_result;
    Int4 i, index, hspcnt, seq_no;
    Nlm_FloatHi best_evalue = 1000.0; /* any large will work... */
    Int4	high_score = 0; 
    
    /* Pre-calculating number of HSPs for next new_result */
    hspcnt = 0;
    seq_no = -1;
    for(index = 0; index < result->hspcnt; index++) {
        if(result->hsp_array[index].number < 0)
            continue;
        if(seq_no < 0) {
            seq_no = result->hsp_array[index].number;
            hspcnt++;
            best_evalue = MIN(best_evalue, result->hsp_array[index].e_value);
            high_score = MAX(high_score, result->hsp_array[index].score);
        } else if (seq_no == result->hsp_array[index].number) {
            hspcnt++;
            best_evalue = MIN(best_evalue, result->hsp_array[index].e_value);
            high_score = MAX(high_score, result->hsp_array[index].score);
        }
        continue;
    }

    if(seq_no == -1) /* Nothing new found... */
        return NULL;
    
    new_result = BLASTResultHitlistNew(hspcnt);
    
    new_result->best_evalue = best_evalue;
    new_result->high_score = high_score;
    new_result->subject_id = seq_no; /*  Fake ... in fact */
    
    for(i = 0, index = 0; index < result->hspcnt; index++) {
        if(result->hsp_array[index].number == seq_no) {
            
            MemCpy(&new_result->hsp_array[i], &result->hsp_array[index], 
                   sizeof(BLASTResultHsp));

            result->hsp_array[index].number = -1;
            
            new_result->hsp_array[i].point_back = new_result;
            /* new_result->hsp_array[i].gap_info = NULL; */
            new_result->hsp_array[i].back_left = NULL;
            new_result->hsp_array[i].back_right = NULL;
            i++;
        }
    }

    /* Here we have to set correctly subject ID and query ID */
    
    return new_result;
}

static Nlm_FloatHi getEvalueFromSeqAlign(SeqAlignPtr thisSeqAlign)
{
    ScorePtr thisScorePtr;
    
    thisScorePtr = thisSeqAlign->score;
    while ((thisScorePtr != NULL) &&
           (StringICmp(thisScorePtr->id->str, "e_value") != 0) &&
           (StringICmp(thisScorePtr->id->str, "sum_e") != 0))
        thisScorePtr = thisScorePtr->next;
    if(NULL == thisScorePtr)
        return(10.0);
    else
        return((Nlm_FloatHi) (thisScorePtr->value.realvalue));
}

Boolean RPSUpdateResult(BlastSearchBlkPtr search, ReadDBFILEPtr rdfp)
{
    BLASTResultsStructPtr result_struct, result_struct_new;
    Int4 index;
    BLASTResultHitlistPtr new_result;
    
    result_struct = search->result_struct;

    if(result_struct->hitlist_count == 0) /* No brain no pain */
        return TRUE;
    
    result_struct_new = BLASTResultsStructNew(search->result_size, search->pbp->max_pieces, search->pbp->hsp_range_max);
    
    if(result_struct->hitlist_count != 1) {
        ErrPostEx(SEV_ERROR, 0,0, "RPSUpdateResult: This function works only for hitlist_count == 1");
        return FALSE;
    }

    /* First we have to update coordinates */
    RPSUpdateCoordinates(rdfp, result_struct->results[0], TRUE);
    
    /* Now we will extract correct "results" one by one */
    for(index = 0; index < result_struct_new->hitlist_max; index++) {
        
        new_result = RPSExtractNewResult(result_struct->results[0]);
        if(new_result == NULL)
            break;
        
        result_struct_new->results[index] = new_result;
        result_struct_new->hitlist_count++;
    }

    search->result_struct = result_struct_new;
    BLASTResultsStructDelete(result_struct);

    return TRUE;
}

/*
	this function calculates the corrected lambda ratio and rescales the
	matrix.  This is done in the IMPALA manner.

	The arguments subject_length and subject_seq are REALLY the length and sequence
	for the query, as the 'BLAST' terminology is used.
*/


static Boolean
RPSimpalaStatCorrections(RPSequencePtr rpseq, Nlm_FloatHiPtr LambdaRatio, Nlm_FloatHi scalingFactor, Int4 subject_length, Uint1Ptr subject_seq)

{
   Nlm_FloatHi *scoreArray; /*array of score probabilities*/
   Nlm_FloatHi *resProb; /*array of probabilities for each residue*/
   BLAST_ScoreFreqPtr this_sfp, return_sfp; /*score frequency pointers to compute lambda*/
   Int4Ptr PNTR posMatrix; /* position-specific matrix. */
   Int4Ptr PNTR copyMatrix; /* copy of position-specific matrix. */
   Nlm_FloatHi initialUngappedLambda, scaledInitialUngappedLambda, correctUngappedLambda;
   Nlm_FloatHi temp1; /*intermediate variable for adjusting matrix*/
   Int4 temp2; /*intermediate variable for adjusting matrix*/
   Int4 seqlength; /* length of posMatrix (or target sequence). */
   Int4 index, inner_index; 

   if (rpseq == NULL)
	return FALSE;

   posMatrix = rpseq->posMatrix;

   resProb = (Nlm_FloatHi *) MemNew (PRO_ALPHABET_SIZE * sizeof(Nlm_FloatHi));
   scoreArray = (Nlm_FloatHi *) MemNew(scoreRange * sizeof(Nlm_FloatHi));
   return_sfp = (BLAST_ScoreFreqPtr) MemNew(1 * sizeof(BLAST_ScoreFreq));

   seqlength = rpseq->seqlen;

  IMPALAfillResidueProbability(subject_seq, subject_length, resProb);
   this_sfp =  IMPALAfillSfp(posMatrix, seqlength, resProb, scoreArray, return_sfp, scoreRange);
   initialUngappedLambda = IMPALAfindUngappedLambda("BLOSUM62");
   scaledInitialUngappedLambda = initialUngappedLambda/scalingFactor;
   correctUngappedLambda = impalaKarlinLambdaNR(this_sfp, scaledInitialUngappedLambda);
   *LambdaRatio = correctUngappedLambda/scaledInitialUngappedLambda;

   resProb = MemFree(resProb);
   scoreArray = MemFree(scoreArray);
   return_sfp = MemFree(return_sfp);


   copyMatrix = (BLAST_Score **) MemNew((seqlength+1) * sizeof(BLAST_Score *));
   for (index = 0; index < (seqlength+1); index++)
     copyMatrix[index] = (BLAST_Score *) MemNew(PRO_ALPHABET_SIZE * sizeof(BLAST_Score ));

/*
   if (correctUngappedLambda > scaledInitialUngappedLambda) 
   {
*/
	   for (index = 0; index < seqlength+1; index++) {
	     for (inner_index = 0; inner_index < PRO_ALPHABET_SIZE; inner_index++) {
	       if ((posMatrix[index][inner_index] == BLAST_SCORE_MIN) || (Xchar == inner_index))
		 copyMatrix[index][inner_index] = posMatrix[index][inner_index];
	       else {
		 temp1 = ((Nlm_FloatHi) (posMatrix[index][inner_index]));
		 temp1 = temp1 * (*LambdaRatio);
		 temp2 = Nlm_Nint(temp1);
		 copyMatrix[index][inner_index] = temp2;
	       }
	     }
	   }
/*
    }
*/

    rpseq->copyMatrix = copyMatrix;

    return TRUE;
}
SeqAlignPtr RPSAlignTraceBack(BlastSearchBlkPtr search, RPSInfoPtr rpsinfo,
                              Uint1Ptr subject_seq, SeqLocPtr slp, BioseqPtr
                              subject_bsp)
{
    SeqAlignPtr seqalign;
    BLASTResultsStructPtr result_struct, result_struct_new;
    Int4 index, index1, subject_length, rev_subject_length;
    Int4Ptr PNTR posMatrix;
    BLASTResultHitlistPtr new_result;
    RPSapSortPtr rpssp;
    RPSequencePtr rpseq;
    Nlm_FloatHi e_value, lambda_ratio;
    Uint1Ptr subject = NULL, rev_subject = NULL;
    SeqPortPtr spp;

    result_struct = search->result_struct;

    if(result_struct->hitlist_count == 0) /* No brain no pain */
        return NULL;
    
    result_struct_new = BLASTResultsStructNew(search->result_size, search->pbp->max_pieces, search->pbp->hsp_range_max);
    
    if(result_struct->hitlist_count != 1) {
        ErrPostEx(SEV_ERROR, 0,0, "RPSUpdateResult: This function works only for hitlist_count == 1 (%d)", result_struct->hitlist_count);
        return NULL;
    }

    subject_length = SeqLocLen(slp);

    search->result_struct = result_struct_new;

    search->subject_info = BLASTSubjectInfoNew(SeqIdDup(SeqIdFindBest(subject_bsp->id, SEQID_GI)), StringSave(BioseqGetTitle(subject_bsp)), subject_length);    
    rpssp =  RPSapSortInit();
    
    /* First we have to update coordinates */
    RPSUpdateCoordinates(rpsinfo->rdfp, result_struct->results[0], FALSE);
    
    /* Now we will extract correct "results" one by one 
       and make alignment then */

    result_struct_new->hitlist_count = 1; /* Always will be 1 */
    
    /* We have to create 4na buffer for our query, that is subject */
    if(!rpsinfo->query_is_prot) {
        spp = SeqPortNew(subject_bsp, 0, -1, Seq_strand_plus, 
                         Seq_code_ncbi4na);
        /* make one longer to "protect" ALIGN. */
        subject = MemNew((1+subject_bsp->length)*sizeof(Uint1));
        for (index1 = 0; index1 < subject_bsp->length; index1++) {
            subject[index1] = SeqPortGetResidue(spp);
        }
        
        /* Gap character in last space. */
        subject[subject_bsp->length] = NULLB;
        subject_length = subject_bsp->length;
        spp = SeqPortFree(spp);
        
        spp = SeqPortNew(subject_bsp, 0, -1, Seq_strand_minus, Seq_code_ncbi4na);
        /* make one longer to "protect" ALIGN. */
        rev_subject = MemNew((1+subject_bsp->length)*sizeof(Uint1));
        for (index1=0; index1<subject_bsp->length; index1++) {
            rev_subject[index1] = SeqPortGetResidue(spp);
        }
        /* Gap character in last space. */
        rev_subject[subject_bsp->length] = NULLB;
        rev_subject_length = subject_bsp->length;
        spp = SeqPortFree(spp);
    }

    /* ---------------------- */
    for(index = 0; index < result_struct_new->hitlist_max; index++) {
        
        new_result = RPSExtractNewResult(result_struct->results[0]);
        if(new_result == NULL)
            break;

        new_result->subject_info = BLASTSubjectInfoNew(SeqIdDup(SeqIdFindBest(subject_bsp->id, SEQID_GI)), StringSave(BioseqGetTitle(subject_bsp)), subject_length);        
        result_struct_new->results[0] = new_result;
        
        /* subject_id == query_id in rdfp database */
        rpseq = RPSGetSequence(rpsinfo, new_result->subject_id);
        
        RPSubstituteQueryLookup(search, rpseq, FALSE);
	/* Correct statistics for proteins, should also be doen for nucleotides. */
	if (rpsinfo->query_is_prot == TRUE) {
            RPSimpalaStatCorrections(rpseq, &lambda_ratio, 1.0, 
                                     subject_length, subject_seq);
            if (rpseq->copyMatrix) {
                posMatrix = search->sbp->posMatrix;
                search->sbp->posMatrix = rpseq->copyMatrix;
            }
	}

        if(rpsinfo->query_is_prot) {

            if(search->pbp->gapped_calculation == TRUE) {
                seqalign = BlastGetGapAlgnTbck(search, 0, TRUE, FALSE, 
                                               subject_seq, 
                                               subject_length, NULL, 0);
                
            } else {
                seqalign = GetSeqAlignForResultHitList(search, TRUE, FALSE, 
                                                       search->pbp->discontinuous, 
                                                       TRUE, FALSE);
            }
        } else { /* tblastn */
            
            seqalign = BlastGetGapAlgnTbck (search, 0, TRUE, FALSE, 
                                            subject, subject_length, 
                                            rev_subject, rev_subject_length);
        }

	/* Deallocate copyMatrix to prevent a leak. */
	if (rpseq->copyMatrix)
	{
		for (index1=0; index1<(1+rpseq->seqlen); index1++)
			MemFree(rpseq->copyMatrix[index1]);
		rpseq->copyMatrix = MemFree(rpseq->copyMatrix);
    		search->sbp->posMatrix = posMatrix;
	}
        
        if(seqalign != NULL && (e_value = getEvalueFromSeqAlign(seqalign)) < 
           search->pbp->cutoff_e) {
            AdjustOffSetsInSeqAlign(seqalign, slp, search->query_slp);
            RPSAddSap(rpssp, seqalign);
        } else {
            SeqAlignSetFree(seqalign);
        }
        
        seqalign = NULL; 
        
        RPSequenceFree(rpseq);
        BLASTResultHitlistFree(new_result);
       
        result_struct_new->results[0] = NULL;
    }

    result_struct_new->hitlist_count = 0;
    
    BLASTResultsStructDelete(result_struct);

    seqalign = RPSReadSapSort(rpssp); 
    
    MemFree(rev_subject);
    MemFree(subject);

    RPSapSortFree(rpssp);

    return seqalign;
}
#if 0
void RRRPrintStatistics(RPSequencePtr rpseq)
{
    ModLAEntry *mod_lt;
    Int4 i, index, zero_count = 0;
    FILE *fd;
    Int4 sss[120], sss1[100];

    mod_lt =  rpseq->mod_lt;
    
    fd = FileOpen("statistics.out", "w");
    memset(sss, 0, sizeof(sss));
    memset(sss1, 0, sizeof(sss1));

    for(i = 0; i < RPS_ARRAY_SIZE; i++) {
        fprintf(fd, "%d %d\n", i, mod_lt[i].num_used);

        if(mod_lt[i].num_used != 0) {
            index = mod_lt[i].num_used / 100;
            
            if(index == 0)
                sss1[mod_lt[i].num_used]++;
            
            sss[index] ++;
        } else {
            zero_count++;
        }
    }

    FileClose(fd);

    fd = FileOpen("distribution.out", "w");

    for(i=0; i < 120; i++) {

        fprintf(fd, "%d %d\n", i, sss[i]);
        
    }

    FileClose(fd);
    printf("Zero count = %d\n", zero_count);

    fd = FileOpen("distribution1.out", "w");
    
    for(i=0; i < 100; i++) {
        
        fprintf(fd, "%d %d\n", i, sss1[i]);
        
    }

    FileClose(fd);
    
    exit(1);
    return;
}
#endif
SeqAlignPtr RPSBlastSearch (BlastSearchBlkPtr search,
                            BioseqPtr query_bsp, RPSInfoPtr rpsinfo)
{
    Int2 status;
    Int4 index, MaxProfiles;
    SeqAlignPtr seqalign=NULL, head = NULL, seqalign_var = NULL;
    SeqPortPtr spp;
    Uint1Ptr subject_seq, subject_seq_start;
    Uint1 residue;
    Int4 subject_length;
    SeqLocPtr slp;
    BioseqPtr subject_bsp;
    RPSequencePtr rpseq;
    BioseqPtr bsp = NULL;
    SPCompressPtr spc = NULL;
    ValNodePtr vnp;
    SeqLocPtr filter_slp;
    
    if (search == NULL || search->query_invalid)
        return NULL;
    

    subject_bsp = query_bsp;    /* Reversed ... */
    
    slp = NULL;
    ValNodeAddPointer(&slp, SEQLOC_WHOLE, 
                      SeqIdDup(SeqIdFindBest(subject_bsp->id, SEQID_GI)));
    
    subject_length = SeqLocLen(slp);
    
    search->subject_info = BLASTSubjectInfoDestruct(search->subject_info);

    /* Extracting sequence from the Bioseq */
    
    if(rpsinfo->query_is_prot) {         /* For blastp search */
        subject_seq_start = subject_seq = NULL;
        
        subject_seq_start = (Uint1Ptr) MemNew(((subject_bsp->length)+2)*sizeof(Uint1));
        /* The first residue is the sentinel. */
        subject_seq_start[0] = NULLB;
        subject_seq = subject_seq_start+1;
        index = 0;
        spp = SeqPortNewByLoc(slp, Seq_code_ncbistdaa);
        while ((residue=SeqPortGetResidue(spp)) != SEQPORT_EOF) {
            if (IS_residue(residue)) {
                subject_seq[index] = residue;
                index++;
            }
        }
        subject_seq[index] = NULLB;
        spp = SeqPortFree(spp);
    } else { /* tblastn */
        spp = SeqPortNewByLoc(slp, Seq_code_ncbi4na);
        subject_bsp = BioseqFindCore(SeqLocId(slp));
        if (subject_bsp != NULL && subject_bsp->repr == Seq_repr_delta)
            SeqPortSet_do_virtual(spp, TRUE);
        spc = SPCompressDNA(spp);
        if (spc == NULL)
            return NULL;
        subject_seq_start = subject_seq = spc->buffer;
        spp = SeqPortFree(spp);

        /* Adjusting translation buffer */
        MemFree(spc);
        MemFree(search->translation_buffer);
        search->translation_buffer = MemNew((3+(subject_length/3))*sizeof(Uint1));
        search->translation_buffer_size = 1 + (subject_length/3);
    }
    
    if(search->context[0].query->sequence_start != NULL) {
        MemFree(search->context[0].query->sequence_start);
        search->context[0].query->sequence_start = NULL;
    }

    /* Cleaning up word finder */

    RPSLookupCleanUp(search->wfp->lookup);
    
    search->subject_info = BLASTSubjectInfoNew(SeqIdDup(SeqIdFindBest(subject_bsp->id, SEQID_GI)), StringSave(BioseqGetTitle(subject_bsp)), subject_length);        
    rpseq = RPSGetBIGSequence(rpsinfo, &bsp);

    /*    RRRPrintStatistics(rpseq); */
  
    RPSubstituteQueryLookup(search, rpseq, TRUE);

    /* Filter the subject sequence, that is actually query */    
    if(rpsinfo->query_is_prot && search->mask != NULL) {
        BlastMaskTheResidues(subject_seq, subject_length, 21, 
                             search->mask->data.ptrvalue, FALSE, 0);
    }

    /* Performing real search here */    
    search = BLASTPerformSearch(search, subject_length, subject_seq);

    if(!rpsinfo->query_is_prot) {
        /* status = BlastLinkHsps(search); */
    } else if (search->pbp->gapped_calculation == FALSE) {
        if (search->pbp->do_sum_stats == TRUE)
            status = BlastLinkHsps(search);
        else
            status = BlastGetNonSumStatsEvalue(search);
    }
    
    status = BlastReapHitlistByEvalue(search);
    
    BlastSaveCurrentHitlist(search);
    
    seqalign = RPSAlignTraceBack(search, rpsinfo, subject_seq, 
                                 slp, subject_bsp); 
    
    /* Now converting protein filter SeqLoc-s into DNA seqlocs */
    if(!rpsinfo->query_is_prot) {  
        for(vnp = search->mask; vnp != NULL; vnp = vnp->next) {
            filter_slp = (SeqLocPtr) vnp->data.ptrvalue;
            BlastConvertProteinSeqLoc(filter_slp, DefineToFrame(vnp->choice), 
                                      subject_bsp->length);
        }
    }

    search->context[0].query->sequence_start = NULL;
    search->context[0].query->sequence = NULL;

    MemFree(subject_seq_start);
    RPSequenceFree(rpseq);
    BioseqFree(bsp);

    search->sbp->posMatrix = NULL;
    search->wfp->lookup->mod_lt = NULL;
    search->wfp->lookup->mod_lookup_table_memory = NULL;

    SeqLocFree(slp);
    
    return seqalign;
}

/* These functions may be never be used ... */
Int4Ptr PNTR RPSReadPSMatrix(CharPtr filename, Int4Ptr mat_len)
{
    Int4Ptr PNTR psmatrix;
    FILE *fd;
    Int4 i, length, num, bytes;
    
    if((length = FileLength(filename)) <= 0)
        return NULL;
    
    num = RPS_ALPHABET_SIZE*sizeof(Uint4);
    
    if(length%num) {
        ErrPostEx(SEV_ERROR, 0,0, "Invalid size of the matrix %s", filename);
        return NULL;
    }
    
    if((fd = FileOpen(filename, "r")) == NULL)
        return NULL;

    *mat_len = length/num;
    
    psmatrix = Nlm_Malloc((*mat_len + 1) * sizeof(Int4Ptr));
    
    for(i = 0 ; i < *mat_len; i++ ){
        psmatrix[i] = MemNew(RPS_ALPHABET_SIZE * sizeof(Int4));
        
        if((bytes = FileRead(psmatrix[i], 1, num, fd)) != num) {
            ErrPostEx(SEV_ERROR, 0,0, 
                      "Failure to read matrix from %s", filename);
            MemFree(psmatrix);
            return NULL;
        }
    }
    
    psmatrix[*mat_len] = MemNew(RPS_ALPHABET_SIZE * sizeof(Int4));
    for(i = 0; i < RPS_ALPHABET_SIZE; i++)
        psmatrix[*mat_len][i] = -INT2_MAX;
    
    *mat_len = length/num;

    return psmatrix;
}
BioseqPtr RPSGetBioseqFromMatrix(Int4Ptr PNTR psmatrix, Int4 length)
{
    BioseqPtr bsp;
    ByteStorePtr seqbs;
    CharPtr buffer;
    Int4 i, j, score, old_score, j_max;
    ValNodePtr     vnp = NULL;

    buffer = MemNew(length);

    old_score = -100; j_max = 0;
    for(i = 0; i < length; i++) {
        for(j = 0; j < RPS_ALPHABET_SIZE; j++) {
            score = psmatrix[i][j];
            if(score > old_score)
                j_max = j;
        }
        buffer[i] = j_max;
    }

    
    if((bsp = BioseqNew()) == NULL) {
        ErrPostEx(SEV_ERROR, 0,0, "Failure to allocate Bioseq");
        MemFree(buffer);
        return NULL;
    }
    
    bsp->mol = Seq_mol_aa;
    bsp->seq_data_type = Seq_code_ncbistdaa;
    bsp->repr = Seq_repr_raw;

    bsp->id = MakeNewProteinSeqId (NULL, NULL);

    if((vnp = SeqDescrNew(NULL)) != NULL) {
        vnp->choice = Seq_descr_title;
        vnp->data.ptrvalue = StringSave ("This is generated sequence");
    }
    bsp->descr = vnp;
    bsp->length = length;
 
    BSWrite(bsp->seq_data, buffer, length);
    
    SeqMgrAddToBioseqIndex (bsp);
    
    return bsp;  
}
/* End of functions may be never be used ... */

