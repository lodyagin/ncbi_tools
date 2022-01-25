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

File name: mblast.c

Author: Ilya Dondoshansky

Contents: Mega BLAST functions

Detailed Contents: 

        - Mega BLAST versions of some functions from blast.c returning array 
	of SeqAlignPtrs

	- Functions specific to Mega BLAST

******************************************************************************
 * $Revision: 6.64 $
 *
 * $Log: mblast.c,v $
 * Revision 6.64  2000/07/06 17:25:49  dondosha
 * Pass megablast_full_deflines option to search parameters
 *
 * Revision 6.63  2000/07/03 21:51:05  dondosha
 * For small wordsizes look for 2 hits in a window
 *
 * Revision 6.62  2000/06/30 15:07:11  dondosha
 * Fixed single strand search
 *
 * Revision 6.61  2000/06/29 21:10:37  dondosha
 * Sort seqaligns by evalue only among hits from same database sequence
 *
 * Revision 6.60  2000/06/29 14:52:10  dondosha
 * Typo in previous change
 *
 * Revision 6.59  2000/06/29 14:45:43  dondosha
 * Mask queries with N, not lower case for -Q option output
 *
 * Revision 6.58  2000/06/27 22:21:53  dondosha
 * Set no_check_score to TRUE; Added functions for masked query output
 *
 * Revision 6.57  2000/06/27 14:48:39  dondosha
 * Changed the procedure of HSP inclusion checking before gapped extension to avoid quadratic complexity
 *
 * Revision 6.56  2000/06/26 18:21:19  dondosha
 * Set awake_index to FALSE before any return from MegaBlastSetUpSearchInternalByLoc
 *
 * Revision 6.55  2000/06/26 17:35:30  dondosha
 * If some of the queries are completely masked, the search can still proceed
 *
 * Revision 6.54  2000/06/22 22:22:19  dondosha
 * Fixed a memory leak in MegaBlastNtWordExtend
 *
 * Revision 6.53  2000/06/21 18:03:20  dondosha
 * Use the same hsp_array for gapped extension
 *
 * Revision 6.52  2000/06/19 20:08:06  madden
 * Add last parameter to readdb_get_sequence_ex
 *
 * Revision 6.51  2000/06/19 19:42:35  dondosha
 * Fixed a bug causing megablast to crash when database is not found
 *
 * Revision 6.50  2000/06/19 19:26:56  dondosha
 * Do not go through a double loop on HSPs for small wordsize
 *
 * Revision 6.49  2000/06/07 17:39:51  dondosha
 * Fixed bug causing extra progress messages in nucleotide neighboring
 *
 * Revision 6.48  2000/05/31 14:07:14  dondosha
 * Reallocate memory for the hit stack if it is overflowing
 *
 * Revision 6.47  2000/05/26 19:22:28  dondosha
 * For neighboring test HSP length and percent identity before saving
 *
 * Revision 6.46  2000/05/25 21:00:34  dondosha
 * Set redundant HSPs to NULL in BlastSortUniqHspArray
 *
 * Revision 6.45  2000/05/24 19:51:03  dondosha
 * Initialize qid_array during search set-up; do not use readdb when megablast is used for two sequences
 *
 * Revision 6.44  2000/05/22 18:48:07  dondosha
 * Conform to change in ReadDBFILE structure
 *
 * Revision 6.43  2000/05/17 21:27:02  dondosha
 * Small improvement in handling edit script
 *
 * Revision 6.42  2000/05/17 17:11:13  dondosha
 * Removed some unused variables; fixed comparison function for sorting hsps
 *
 * Revision 6.41  2000/05/12 19:44:19  dondosha
 * Do not use SeqPort for retrieving queries; keep query ids in array to allow binary search; use ncbi4na encoding for subject sequences during extension; added routine BinarySearchInt4
 *
 * Revision 6.40  2000/05/03 20:23:04  dondosha
 * Added function MegaBlastGappedAlign to do the gapped extension after all perfect matches are saved - this significantly improves speed
 *
 * Revision 6.39  2000/05/01 21:25:54  dondosha
 * Changed greedy_gapped_align calls to MegaBlastAffineGreedyAlign; bug fix for neighboring; memory leak fix
 *
 * Revision 6.38  2000/04/25 21:51:45  dondosha
 * Little clean-up in MegaBlastNtWordExtend
 *
 * Revision 6.37  2000/04/20 15:14:36  dondosha
 * Fixed MegaBlastGetFirstAndLastContext and BlastFillQueryOffsets for one strand only search
 *
 * Revision 6.36  2000/04/19 18:54:33  dondosha
 * Bug fix: assign values for end-points different from gap or masked characters
 *
 * Revision 6.35  2000/04/18 19:09:03  dondosha
 * Set the X dropoff correctly for the megablast search
 *
 * Revision 6.34  2000/04/14 19:22:52  dondosha
 * Remove HSPs contained in other HSPs
 *
 * Revision 6.33  2000/04/12 19:11:28  dondosha
 * Create index_callback thread right after the creation of BlastSearchBlk
 *
 * Revision 6.32  2000/04/12 18:24:28  dondosha
 * Added MegaBlastGetPercentIdentity; fixed bugs with memory handling
 *
 * Revision 6.31  2000/04/11 21:00:19  dondosha
 * Fixed memory bug in BlastNeedHumanRepeatFiltering
 *
 * Revision 6.30  2000/04/11 19:51:50  dondosha
 * Removed setting of queue_callback to NULL
 *
 * Revision 6.29  2000/04/10 18:03:26  dondosha
 * Fixed a bug in previous change
 *
 * Revision 6.28  2000/04/10 17:44:52  dondosha
 * Find correct path for Homo_sapiens.n.gil, get this gilist only once
 *
 * Revision 6.27  2000/04/10 15:22:35  dondosha
 * Moved call to BlastFillQueryOffsets to MegaBlastSetUpSearchInternalByLoc
 *
 * Revision 6.26  2000/04/07 16:51:57  dondosha
 * Include callback argument for handling results in BioseqMegaBlastEngine, to be initialized in Main
 *
 * Revision 6.25  2000/04/05 18:14:48  dondosha
 * Moved SeqIdSetDup to objloc.c, slightly changed format of BlastSearchHandleResults
 *
 * Revision 6.24  2000/04/04 20:50:46  dondosha
 * Cleaned from non-megablast (non-blastn) stuff
 *
 * Revision 6.23  2000/04/04 16:16:35  dondosha
 * Fixed some memory leaks in MegaBlast traceback
 *
 * Revision 6.22  2000/04/03 23:37:08  dondosha
 * Added test for human repeat filtering for neighboring
 *
 * Revision 6.21  2000/03/31 19:10:32  dondosha
 * Changed some names related to MegaBlast
 *
 * Revision 6.20  2000/03/29 23:32:21  dondosha
 * Fixed minor bugs from previous change
 *
 * Revision 6.19  2000/03/29 22:13:43  dondosha
 * Enabled processing gap information for MegaBlast code
 *
 * Revision 6.18  2000/03/27 20:56:43  madden
 * Set query frame properly in BlastAdjustHitOffsets (for ungapped blast)
 *
 * Revision 6.17  2000/03/23 20:02:02  dondosha
 * BlastSearchHandleResults exits immediately if hspcnt==0
 *
 * Revision 6.16  2000/03/22 17:57:12  dondosha
 * Improved memory management in MegaBlast
 *
 * Revision 6.15  2000/03/16 18:14:00  dondosha
 * Added routine SeqIdSetDup, improved handling of megablast results
 *
 * Revision 6.14  2000/03/13 21:06:59  dondosha
 * Routine printing results rewritten, plus minor improvements
 *
 * Revision 6.13  2000/03/08 20:34:53  madden
 * Remove ungapped psi-blast stuff
 *
 * Revision 6.12  2000/03/03 18:12:22  dondosha
 * Added routine MegaBlastWordFinderDeallocate to fix multithreaded MegaBlast
 *
 * Revision 6.11  2000/03/02 18:29:00  dondosha
 * Minor bug fix in setting context offsets
 *
 * Revision 6.10  2000/03/02 17:23:21  dondosha
 * Added lower case masking plus bug fixes
 *
 * Revision 6.9  2000/02/24 23:20:15  dondosha
 * Fixed a bug in BlastAdjustHitOffsets
 *
 * Revision 6.8  2000/02/24 17:53:38  dondosha
 * Added routine BlastAdjustHitOffsets
 *
 * Revision 6.7  2000/02/23 20:28:38  dondosha
 * Bug fix for ungapped alignment search, plus style improvements
 *
 * Revision 6.6  2000/02/17 20:01:53  dondosha
 * Removed references to theCacheSize
 *
 * Revision 6.5  2000/02/12 21:18:57  kans
 * added prototype for MegaBlastBuildLookupTable - implemented in lookup.c, called from mblast.c
 *
 * Revision 6.4  2000/02/11 21:03:36  dondosha
 * Added new word finder and extension routines to work with new lokup table
 *
 * Revision 6.3  2000/02/03 21:21:10  dondosha
 * Added header; few bug fixes
 *
 * Revision 6.2  2000/02/02  15:03:44  dondosha
 * Removed unused routine ReapHitlistByContext
 *
 * */

#include <time.h>
#include <ncbi.h>
#include <blastpri.h>
#include <lookup.h>
#include <objcode.h>
#include <objseq.h>
#include <sequtil.h>
#include <tofasta.h>
#include <seqport.h>
#include <readdb.h>
#include <ncbithr.h>
#include <gapxdrop.h>
#include <dust.h>
#include <mbutils.h>
#include <mbalign.h>
#include <mblast.h>

SeqAlignPtr PNTR
BioseqMegaBlastEngine (BioseqPtr PNTR bspp, CharPtr progname, CharPtr database,
		       BLAST_OptionsBlkPtr options, ValNodePtr *other_returns,
		       ValNodePtr *error_returns, int (LIBCALLBACK
						       *callback)(Int4 done,
								  Int4
								  positives),
		       SeqIdPtr seqid_list, BlastDoubleInt4Ptr gi_list, 
		       Int4 gi_list_total, SeqLocPtr PNTR mask_slpp, 
		       int (LIBCALLBACK *results_callback)PROTO((VoidPtr Ptr)))
{
	SeqLocPtr slp;

	Boolean options_allocated=FALSE;
	BlastSearchBlkPtr search;
	Int2 status;
	SeqAlignPtr PNTR head;
	SeqLocPtr whole_slp=NULL, slp_var;
	Int4 index;


	slp = NULL;
	for (index=0; bspp[index] != NULL; index++) 
	   ValNodeAddPointer(&slp, SEQLOC_WHOLE,
			  SeqIdDup(SeqIdFindBest(bspp[index]->id,
			  SEQID_GI)));

	head = NULL;

	if (error_returns)
	{
		*error_returns = NULL;
	}

	if (other_returns)
	{
		*other_returns = NULL;
	}

	if (progname == NULL)
		return NULL;

	/* If no options, use default. */
        if (options == NULL)
	{
		options = BLASTOptionNew(progname, FALSE);
		options_allocated = TRUE;
	}

	status = BLASTOptionValidateEx(options, progname, error_returns);
	if (status != 0)
	{	/* error messages in other_returns? */
		return NULL;
	}

	if (slp == NULL || database == NULL)
		return NULL;

	search = MegaBlastSetUpSearchWithReadDbInternal(slp, NULL, progname,
							0,
							database, options, NULL,
							seqid_list, gi_list,
							gi_list_total, NULL, mask_slpp);

	if (search == NULL)
	{
		return NULL;
	}

	search->thr_info->tick_callback = callback;
	search->thr_info->star_callback = callback;
	search->handle_results = results_callback;

	head = BioseqMegaBlastEngineCore(search, options, NULL);
	
	if (search->error_return)
	{
		ValNodeLink(error_returns, search->error_return);
		search->error_return = NULL;
	}

	if (other_returns)
	{ /* format dbinfo etc.  */
		*other_returns = BlastOtherReturnsPrepare(search);
	}

	if (options_allocated)
	{
		options = BLASTOptionDelete(options);
	}
 
	search->rdfp = ReadDBCloseMHdrAndSeqFiles(search->rdfp);
	search->rdfp = ReadDBFreeSharedInfo(search->rdfp);
	search = BlastSearchBlkDestruct(search);

	/* Adjsut the offset if the query does not cover the entire sequence. */
	slp_var = slp;
	for (index=0; slp_var!=NULL; index++,slp_var=slp_var->next) {
	if (slp_var->choice != SEQLOC_WHOLE)
	{
		ValNodeAddPointer(&whole_slp, SEQLOC_WHOLE, SeqIdFindBest(SeqLocId(slp_var), SEQID_GI));
		if (SeqLocAinB(whole_slp, slp_var) != 0)
		{
			AdjustOffSetsInSeqAlign(head[index], slp_var, NULL);
		}
		ValNodeFree(whole_slp);
	}
	}

	SeqLocSetFree(slp);
	return head;
}

SeqAlignPtr PNTR
BioseqMegaBlastEngineCore(BlastSearchBlkPtr search, BLAST_OptionsBlkPtr options,
			  Int4Ptr *pos_matrix)
{
	BLASTResultsStructPtr result_struct;
	Int4 hitlist_count, index;
	SeqAlignPtr head, seqalign, seqalign_var, PNTR seqalignp;
	SeqIdPtr gi_list=NULL;

	head = seqalign = NULL;
	seqalignp = NULL;
start_timer;
	if (search == NULL || search->query_invalid)
		return NULL;

	/* Starting awake thread if multithreaded. */
	if (search->searchsp_eff > AWAKE_THR_MIN_SIZE)
		BlastStartAwakeThread(search->thr_info);

	stop_timer("BioseqBlastEngineCore: before do_the_blast_run()");
	do_the_blast_run(search);
	
start_timer;
	head = NULL;

	search->sbp->kbp_gap[search->first_context] = search->sbp->kbp[search->first_context];

	search->pbp->gap_open = options->gap_open;
	search->pbp->gap_extend = options->gap_extend;
	
	result_struct = search->result_struct;
	hitlist_count = result_struct->hitlist_count;
	search->sbp->kbp_gap[search->first_context] = NULL;
	
	HeapSort(result_struct->results, hitlist_count, sizeof(BLASTResultHitlistPtr), evalue_compare_hits);
	
	/* 
	   The next loop organizes the SeqAligns (and the alignments in 
	   the BLAST report) in the same order as the deflines.
	*/
	head = NULL;
	for (index=0; index<hitlist_count; index++) {
	   seqalign = result_struct->results[index]->seqalign;
	   if (seqalign) {
	      if (head == NULL)
		 head = seqalign;
	      else {
		 for (seqalign_var=head; seqalign_var->next; )
		    seqalign_var = seqalign_var->next;
		 seqalign_var->next = seqalign;
	      }
	   }
	}

	/* Now rearrange the seqaligns into an array so each element
	   contains alignments from a single query sequence */
	seqalignp = MegaBlastPackAlignmentsByQuery(search, head);

	gi_list = SeqIdSetFree(gi_list);

	if (seqalignp==NULL) {
	   seqalignp = (SeqAlignPtr PNTR) MemNew(sizeof(SeqAlignPtr));
	   *seqalignp = head;
	}
	/* Stop the awake thread. */
	BlastStopAwakeThread(search->thr_info);

stop_timer("BioseqBlastEngineCore: after do_the_blast_run()");
	return seqalignp;
}

#define INDEX_THR_MIN_SIZE 20000
BlastSearchBlkPtr
MegaBlastSetUpSearchWithReadDbInternal (SeqLocPtr query_slp, BioseqPtr
					query_bsp, CharPtr prog_name, Int4 qlen, 
					CharPtr dbname, BLAST_OptionsBlkPtr
					options, int (LIBCALLBACK *callback)
					PROTO((Int4 done, Int4 positives)), 
					SeqIdPtr seqid_list, 
					BlastDoubleInt4Ptr gi_list, Int4
					gi_list_total, ReadDBFILEPtr rdfp, 
					SeqLocPtr PNTR mask_slpp)

{

	BlastSearchBlkPtr search;
	Boolean multiple_hits, options_alloc=FALSE;
	Int2 status, first_context, last_context;
        Int8	dblen;
	Int4	query_length;
        ValNodePtr	vnp;
        Int4		i;
	Nlm_FloatHi	searchsp_eff=0;
	Boolean		use_private_gilist = FALSE;
	OIDListPtr	alias_oidlist;
	Int4		mask_index, virtual_mask_index;
	Uint4		oid_bit, virtual_oid_bit;
	ReadDBFILEPtr	tmprdfp;
        
	/* Allocate default options if none are allocated yet. */
	if (options == NULL)
	{
		options = BLASTOptionNew(prog_name, FALSE);
		options_alloc = TRUE;
	}

     
	/* non-NULL gi_list means that standalone program called the function */
	/* non-NULL options->gilist means that server got this gilist from client */
        if(options->gilist) {
                /* translate list of gis from ValNodePtr to BlastDoubleInt4Ptr */

            gi_list = Nlm_Malloc(ValNodeLen(options->gilist) * sizeof(BlastDoubleInt4));
            for (vnp=options->gilist, i=0; vnp; vnp = vnp->next, ++i) {
                gi_list[i].gi = vnp->data.intvalue;
            }
            gi_list_total = i;
        }
        
	/* Using "options->gifile" file and gi_list,
	   construct new gi_list with all needed gis */


        if (options->gifile && StringCmp(options->gifile, "")) {
            Int4	gi_list_total_2;
            BlastDoubleInt4Ptr	gi_list_2, tmptr;
            Int4	size = sizeof(BlastDoubleInt4);
            Char	buf[PATH_MAX], blast_dir[PATH_MAX];

                /**
                 * first looking in current directory, then checking .ncbirc,
                 * then $BLASTDB and then assuming BLASTDB_DIR
                 */
	    if (FileLength(options->gifile) > 0) {
	    	char *path = Nlm_FilePathFind(options->gifile);
	    	if (StringLen(path) > 0) {
	    		StringCpy(blast_dir, path);
	    	}
	    	else {
	    		StringCpy(blast_dir, ".");
	    	}
	    	MemFree(path);
	    }
	    else {
#ifdef OS_UNIX
		if (getenv("BLASTDB"))
		    Nlm_GetAppParam("NCBI", "BLAST", "BLASTDB", getenv("BLASTDB"), blast_dir, PATH_MAX);
		else
#endif
		    Nlm_GetAppParam ("NCBI", "BLAST", "BLASTDB", BLASTDB_DIR, blast_dir, PATH_MAX);
	    }
            sprintf(buf, "%s%s%s", blast_dir, DIRDELIMSTR, options->gifile);

            gi_list_2 = GetGisFromFile(buf, &gi_list_total_2);

                /* replace or append this list to main one */
            if (gi_list && gi_list_2) {
                    /* append */
                tmptr = Malloc((gi_list_total+gi_list_total_2)*size);
                MemCpy(tmptr, gi_list, gi_list_total * size);
                MemCpy(tmptr+gi_list_total, gi_list_2, gi_list_total_2 * size);

                MemFree(gi_list);
                MemFree(gi_list_2);
                
                gi_list = tmptr;
                gi_list_total += gi_list_total_2;
            }
            else if (gi_list_2) {
                    /* replace */
                gi_list = gi_list_2;
                gi_list_total = gi_list_total_2;
            }
        }
	if (options->window_size != 0)
		multiple_hits = TRUE;
	else
		multiple_hits = FALSE;

	/* last context is the total number of contexts in concatenated 
	   queries */
	MegaBlastGetFirstAndLastContext(prog_name, query_slp, &first_context, &last_context, options->strand_option);

	if (query_slp)
	   query_length = SeqLocTotalLen(prog_name, query_slp);
	else
	   query_length = query_bsp->length;
		
	/* Pass 0 length, since we don't need to allocate ewp here */
	search = BlastSearchBlkNewExtra(options->wordsize, 0,
					dbname, multiple_hits,
					options->threshold_first,
					options->threshold_second,
					options->hitlist_size, prog_name, NULL,
					first_context, last_context, rdfp,
					options->window_size);

	if (search) {
	   if (NlmThreadsAvailable() && query_length > INDEX_THR_MIN_SIZE) {
	      search->thr_info->awake_index = TRUE;
	      search->thr_info->last_tick = Nlm_GetSecs();
	      search->thr_info->index_thr = 
		 NlmThreadCreate(index_proc, search->thr_info);
	      search->thr_info->index_callback = callback;
	   }
	   readdb_get_totals(search->rdfp, &(dblen), &(search->dbseq_num));

	   if (seqid_list)
	      BlastAdjustDbNumbers(search->rdfp, &(dblen), 
				   &(search->dbseq_num), seqid_list, NULL, 
				   NULL, NULL, 0);

		if (gi_list) {
			/* transform the list into OID mask */

		    Int4		i;
		    Int4		maxoid, virtual_oid, oid;
		    OIDListPtr		oidlist;
		    Int4		total;
		    Int4		start;
		    Boolean		done;

		    BlastDoubleInt4Ptr PNTR gi_list_pointers;

		    use_private_gilist = TRUE;
		    gi_list_pointers = Nlm_Malloc(gi_list_total*sizeof(BlastDoubleInt4Ptr));
		    maxoid = 0;
		    for (i=0; i < gi_list_total; i++) {
			/* get virtual OID and start position for the 
			   database this gi is in */
			gi_list[i].ordinal_id = readdb_gi2seq(search->rdfp, gi_list[i].gi, &start);
			gi_list[i].start = start;
			maxoid = MAX(maxoid, gi_list[i].ordinal_id);
			gi_list_pointers[i] = &(gi_list[i]);
		    }

		    /* allocate space for mask for virtual database */
		    oidlist = (OIDListPtr) MemNew(sizeof(OIDList));
		    oidlist->total = maxoid + 1;
		    total = maxoid/MASK_WORD_SIZE + 2;
		    oidlist->list = (Uint4Ptr) MemNew (total*sizeof(Int4));
		    oidlist->memory = oidlist->list;
		    /* Merge this list with virtual database (OID list) */

		    for (i=0; i < gi_list_total; i++) {
			/* get start possition in that database */
			start = gi_list[i].start;

			/* find out if this is an mask database */

			done = FALSE;
			tmprdfp = search->rdfp;

			alias_oidlist = NULL;
			while (tmprdfp && !done) {
			    if (tmprdfp->start == start) {
				alias_oidlist = tmprdfp->oidlist;
				done = TRUE;
			    } else if (tmprdfp->start > start) {
				done = TRUE;
			    } else {
				tmprdfp = tmprdfp->next;
			    }
			}

			/* populate the mask */
			virtual_oid = gi_list[i].ordinal_id;

			if (virtual_oid >= 0) {
			    virtual_mask_index = virtual_oid/MASK_WORD_SIZE;
			    if (alias_oidlist) {
				mask_index = (virtual_oid - start) / MASK_WORD_SIZE;
				oid_bit = 0x1 << (MASK_WORD_SIZE - 1 - (virtual_oid-start) % MASK_WORD_SIZE);
			    }

			    virtual_oid_bit = 0x1 << (MASK_WORD_SIZE - 1 - virtual_oid % MASK_WORD_SIZE);

			    if ((!alias_oidlist) ||
				    (alias_oidlist && alias_oidlist->list && 
				     (Nlm_SwapUint4(alias_oidlist->list[mask_index])) & oid_bit)) { 
				oidlist->list[virtual_mask_index] |= virtual_oid_bit;
			    }
			}
		    }
		    for (i=0; i<total; i++) {
			oidlist->list[i] = Nlm_SwapUint4(oidlist->list[i]);
		    }

		    search->rdfp->oidlist = oidlist;

		    search->rdfp->parameters |= READDB_CONTENTS_ALLOCATED;
		    /*search->rdfp->contents_allocated = TRUE;*/

		    /* in this case, the case when we have .gil file, the only database mask
		       should be used in Blast Search, so set number of sequences for the first
		       database in rdfp list to 0 avoiding search this real database: */
		    search->rdfp->num_seqs = 0;

		    /* Adjust db size; the size should be equal to size of the
		       .gil */
		    if (!options->use_real_db_size)
		       BlastAdjustDbNumbers(search->rdfp, &(dblen), &(search->dbseq_num), 
		      NULL, NULL, oidlist, NULL, 0);

		    /* keep list of gi's (needed for formating) */
		    if (options->sort_gi_list)
		       HeapSort(gi_list_pointers, gi_list_total, sizeof(BlastDoubleInt4Ptr PNTR), compare);
		    search->thr_info->blast_gi_list = BlastGiListNew(gi_list, gi_list_pointers, gi_list_total);
		} else {
		    /* Ok, we do not have a gi-list specified, but maybe
		       we have an a mask database in the list of databases,
		       we need to create one mask for all such databases */
		    OIDListPtr		virtual_oidlist = NULL;
		    Int4		final_virtual_db_seq=0, final_db_seq=0;
		    Int4		mask, oid, virtual_oid, maskindex,
					virtual_mask_index, total_virtual_mask,
					base;
		    Uint4		virtual_oid_bit;

		    tmprdfp = search->rdfp;
		    while (tmprdfp) {

			final_virtual_db_seq = tmprdfp->stop;
			if (!tmprdfp->oidlist)
			    final_db_seq = tmprdfp->stop;
			tmprdfp = tmprdfp->next;
		    }

		    tmprdfp = search->rdfp;
		    while (tmprdfp) {
			if (tmprdfp->oidlist) {
			    if (!virtual_oidlist) {
				/* create new oidlist for virtual database */
				virtual_oidlist = (OIDListPtr) MemNew(sizeof(OIDList));
				virtual_oidlist->total = final_virtual_db_seq + 1;
				total_virtual_mask = final_virtual_db_seq/MASK_WORD_SIZE + 2;
				virtual_oidlist->list = (Uint4Ptr) MemNew (total_virtual_mask*sizeof(Int4));
			    }
			    /* Now populate the virtual_oidlist */
			    maskindex = 0;
			    base = 0;

			    while (maskindex < (tmprdfp->oidlist->total/MASK_WORD_SIZE +1)) {
				/* for each long-word mask */
				mask = Nlm_SwapUint4(tmprdfp->oidlist->list[maskindex]);

				i = 0;
				while (mask) {
				    if (mask & (((Uint4)0x1)<<(MASK_WORD_SIZE-1))) {
					oid = base + i;
					virtual_oid = oid + tmprdfp->start;

					virtual_mask_index = virtual_oid/MASK_WORD_SIZE;
					virtual_oid_bit = 0x1 << (MASK_WORD_SIZE - 1 - virtual_oid % MASK_WORD_SIZE);
					virtual_oidlist->list[virtual_mask_index] |= virtual_oid_bit;
				    }
				    mask <<= 1;
				    i++;
				}
				maskindex++;
				base += MASK_WORD_SIZE;
			    }

			    /* free old mask */
			    tmprdfp->oidlist = OIDListFree(tmprdfp->oidlist);
			}
			tmprdfp = tmprdfp->next;
		    }
		    if (virtual_oidlist) {
			for (i=0; i<total_virtual_mask; i++) {
			    virtual_oidlist->list[i] = Nlm_SwapUint4(virtual_oidlist->list[i]);
			}
		    }
		    search->rdfp->oidlist = virtual_oidlist;

		    readdb_get_totals_ex(search->rdfp, &(dblen), &(search->dbseq_num), TRUE);

		}
		/* Intended for use when a list of gi's is sent in, but the real size is needed. */
		/* It's probably still necessary to call BlastAdjustDbNumbers, but it would be nice
			if this were not required. */
		if (options->use_real_db_size && gi_list)
		      search->dbseq_num = MIN(search->dbseq_num, gi_list_total);

#if 0
		/* use length and num of seqs of the database from alias file */
		if (search->rdfp->aliaslen && !gi_list)
		    dblen = search->rdfp->aliaslen;
		if (search->rdfp->aliasnseq && !gi_list) 
		    search->dbseq_num = search->rdfp->aliasnseq;
#endif
		/* command-line/options trump alias file. */
		if (options->db_length > 0)
			dblen = options->db_length;
		if (options->dbseq_num > 0)
			search->dbseq_num = options->dbseq_num;
		if (options->searchsp_eff > 0)
			searchsp_eff = options->searchsp_eff;

		search->dblen = dblen;
		search->searchsp_eff = searchsp_eff;
		status = MegaBlastSetUpSearchInternalByLoc (search, query_slp, query_bsp,
							prog_name, qlen, options,
							callback, mask_slpp);
		if (status != 0)
		{
	  		ErrPostEx(SEV_WARNING, 0, 0, "MegaBlastSetUpBlastSearch failed.");
			search->query_invalid = TRUE;
		}

		if (search->pbp->is_megablast_search) 
		   search = GreedyAlignMemAlloc(search);
		else 
		   search->abmp = NULL;
		search->rdfp = ReadDBCloseMHdrAndSeqFiles(search->rdfp);
	}
	
	if (options_alloc)
		options = BLASTOptionDelete(options);

	return search;
}

#define DROPOFF_NUMBER_OF_BITS 10.0
#define TOP_HALFBYTE 0xf0
#define BOTTOM_HALFBYTE 0x0f
static Uint1 inverse_halfbyte[16] = 
{ 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };

Int2
MegaBlastSetUpSearchInternalByLoc (BlastSearchBlkPtr search, SeqLocPtr
				   query_slp, BioseqPtr query_bsp, CharPtr
				   prog_name, Int4 qlen, BLAST_OptionsBlkPtr
				   options, int (LIBCALLBACK
						 *callback)PROTO((Int4 done,
								  Int4
								  positives)),
				   SeqLocPtr PNTR mask_slpp)

{
   BioseqPtr bsp;
   Boolean mask_at_hash=FALSE, private_slp_delete;
   Boolean query_is_na, db_is_na;
   Char buffer[128];
   Int2 retval, status, rem;
   Int4 effective_query_length, query_length,
      index, index1, length, length_adjustment=0, last_length_adjustment,
      min_query_length, full_query_length=0;
   Int4 context, last_context, num_queries;
   Nlm_FloatHi avglen;
   SeqIdPtr query_id, qid;
   /*SeqPortPtr spp=NULL, spp_reverse=NULL;*/
   SeqLocPtr filter_slp=NULL, private_slp=NULL, private_slp_rev=NULL, slp,
      tmp_slp;
   Uint1 residue, strand, seq_data_type;
   Uint1Ptr query_seq, query_seq_start, query_seq_rev, query_seq_start_rev;
   Uint1Ptr query_seq_combined, sequence;
   CharPtr filter_string;
   Uint4 query_gi;
   Int4 homo_gilist_size;
   BlastDoubleInt4Ptr homo_gilist;
   CharPtr homo_gifile;
   Nlm_BSUnitPtr seq_chain, seq_chain_start;
   
   if (options == NULL) {
      ErrPostEx(SEV_FATAL, 0, 0, "BLAST_OptionsBlkPtr is NULL\n");
      search->thr_info->awake_index = FALSE;
      return 1;
   }
   
   if (query_slp == NULL && query_bsp == NULL) {
      ErrPostEx(SEV_FATAL, 0, 0, "Query is NULL\n");
      search->thr_info->awake_index = FALSE;
      return 1;
   }
   
   query_seq = NULL;	/* Gets rid of warning. */
   query_seq_rev = NULL;	/* Gets rid of warning. */
   query_seq_start = NULL;	/* Gets rid of warning. */
   query_seq_start_rev = NULL;	/* Gets rid of warning. */
   
   context = search->first_context; 
   
   search = BlastFillQueryOffsets(search, query_slp);
   query_length = SeqLocTotalLen(search->prog_name, query_slp);
   num_queries = search->last_context / 2 + 1;
   search->qid_array = (SeqIdPtr PNTR) Malloc(num_queries*sizeof(SeqIdPtr));

   search->context[search->first_context].query->sequence_start = 
      (Uint1Ptr) Malloc((query_length+2)*sizeof(Uint1));

   if (!query_slp) {
      private_slp = SeqLocIntNew(0, query_bsp->length-1 , Seq_strand_plus, SeqIdFindBest(query_bsp->id, SEQID_GI));
      private_slp_rev = SeqLocIntNew(0, query_bsp->length-1 , Seq_strand_minus, SeqIdFindBest(query_bsp->id, SEQID_GI));
      private_slp_delete = FALSE;
   }

   if (query_slp) {
      search->query_slp = query_slp;
   }	else {
      search->query_slp = private_slp;
      search->allocated += BLAST_SEARCH_ALLOC_QUERY_SLP;
   }
   
   search->translation_buffer = NULL;
   search->translation_buffer_size = 0;
   
   /*
     Set the context_factor, which specifies how many different 
     ways the query or db is examined (e.g., blastn looks at both
     stands of query, context_factor is 2).
   */
   /* All strands of all queries concatenated into one */
   search->context_factor = 1;
   
   /* Set the ambiguous residue before the ScoreBlk is filled. */
   if(options->matrix!=NULL) {
      search->sbp->read_in_matrix = TRUE;
   } else 
      search->sbp->read_in_matrix = FALSE;
   BlastScoreSetAmbigRes(search->sbp, 'N');
   
   search->sbp->penalty = options->penalty;
   search->sbp->reward = options->reward;
   
   /* option is to use alignments chosen by user in PSM computation API (used in WWW PSI-Blast); */
   search->pbp->use_best_align = options->use_best_align;
   
   /* Should culling be used at all? */
   search->pbp->perform_culling = options->perform_culling;
   search->pbp->hsp_range_max = options->hsp_range_max;
   /* This assures that search->pbp->max_pieces is at least one wide. */
   
   search->pbp->block_width = MIN(query_length, options->block_width);
   if (search->pbp->block_width > 0)
      search->pbp->max_pieces = query_length/search->pbp->block_width;
   
   search->sbp->query_length = query_length;
   
   search->result_struct = BLASTResultsStructNew(search->result_size, search->pbp->max_pieces, search->pbp->hsp_range_max);
   if (options->matrix != NULL)
      status = BlastScoreBlkMatFill(search->sbp, options->matrix);
   else
      status = BlastScoreBlkMatFill(search->sbp, "BLOSUM62");
   if (status != 0) {
      ErrPostEx(SEV_WARNING, 0, 0, "BlastScoreBlkMatFill returned non-zero status");
      search->thr_info->awake_index = FALSE;
      return 1;
   }
   
   /* This is used right below. */
   search->pbp->gapped_calculation = options->gapped_calculation;
   search->pbp->do_not_reevaluate = options->do_not_reevaluate;
   search->pbp->do_sum_stats = options->do_sum_stats;
   search->pbp->first_db_seq = options->first_db_seq;
   search->pbp->final_db_seq = options->final_db_seq;
   search->pbp->is_neighboring = options->is_neighboring;
   search->pbp->one_line_results = options->one_line_results;
   search->pbp->is_megablast_search = options->is_megablast_search;
   search->pbp->megablast_full_deflines = options->megablast_full_deflines;

   retval = 0;
   
   context = search->first_context;
   
   if (private_slp && private_slp_delete)
      private_slp = SeqLocFree(private_slp);
   if (private_slp_rev)
      private_slp_rev = SeqLocFree(private_slp_rev);
   
   search->query_id = NULL;
   slp = query_slp;
   
   if (search->pbp->is_neighboring) {
      homo_gifile = FindBlastDBFile("Homo_sapiens.n.gil");
      homo_gilist = GetGisFromFile (homo_gifile, &homo_gilist_size);
      MemFree(homo_gifile);
   }

   while(context<=search->last_context && slp != NULL) {
      if (search->first_context == 0) {
	 private_slp = SeqLocIntNew(SeqLocStart(slp), SeqLocStop(slp), 
				    Seq_strand_plus, SeqLocId(slp));
      }
      if (search->last_context % 2 == 1) {
	 private_slp_rev = SeqLocIntNew(SeqLocStart(slp), SeqLocStop(slp), 
					Seq_strand_minus, SeqLocId(slp));
      }
      
      if (private_slp)
	 tmp_slp = private_slp;
      else 
	 tmp_slp = private_slp_rev;

      query_length = 0;
      query_length = SeqLocLen(tmp_slp);
      if (query_length == 0) {
	 sprintf(buffer, "No valid query sequence");
	 BlastConstructErrorMessage("Blast", buffer, 2, &(search->error_return));
	 continue;
      }

      if (options->filter && !options->filter_string)
	 options->filter_string = BlastConstructFilterString(options->filter);
      if (search->pbp->is_neighboring) {
	 /* Test if we have to do human repeat filtering on this query */
	 if (private_slp)
	    qid = SeqLocId(private_slp);
	 else 
	    qid = SeqLocId(private_slp_rev);
	 query_gi = GetGIForSeqId(qid);
	 if (BlastNeedHumanRepeatFiltering(homo_gilist, homo_gilist_size, 
					   query_gi)) {
	    filter_string = 
	       (CharPtr) Malloc(StringLen(options->filter_string)+3);
	    sprintf(filter_string, "%s%s", options->filter_string, ";R");
	 } else
	    filter_string = StringSave(options->filter_string);
      } else
	 filter_string = options->filter_string;
      
      if (private_slp)
	 filter_slp = BlastSeqLocFilterEx(private_slp, filter_string, &mask_at_hash);
      else if (private_slp_rev)
	 filter_slp = BlastSeqLocFilterEx(private_slp_rev, filter_string, &mask_at_hash);
      if (search->pbp->is_neighboring)
	 MemFree(filter_string);
	/* 
           Dusting of query sequence. Only needed for blastn, optional
        */
	
      if (filter_slp && !mask_at_hash)
	 ValNodeAddPointer(&(search->mask), SEQLOC_MASKING_NOTSET, filter_slp);

      /*
      if (private_slp) {
	 spp = SeqPortNewByLoc(private_slp, Seq_code_ncbi4na);
	 SeqPortSet_do_virtual(spp, TRUE);
      }
      if (private_slp_rev) {
	 spp_reverse = SeqPortNewByLoc(private_slp_rev, Seq_code_ncbi4na);
	 SeqPortSet_do_virtual(spp_reverse, TRUE);
      }
      */
      if (search->first_context==0 && search->last_context % 2 == 1)
	 query_seq_combined = 
	    (Uint1Ptr) Malloc((2*query_length+3)*sizeof(Uint1));
      else
	 query_seq_combined = 
	    (Uint1Ptr) Malloc((query_length+2)*sizeof(Uint1));

      bsp = BioseqLockById(SeqLocId(tmp_slp));
      
      if (bsp == NULL) {
	 ErrPostEx(SEV_WARNING, 0, 0, "No valid query sequence, BioseqLockById returned NULL\n");
         search->thr_info->awake_index = FALSE;
	 return 1;
      }
      seq_data_type = bsp->seq_data_type;
      if (seq_data_type != Seq_code_ncbi4na && 
	  seq_data_type != Seq_code_ncbi2na) {
	 bsp->seq_data = BSConvertSeq(bsp->seq_data, Seq_code_ncbi4na, 
				      bsp->seq_data_type, 
				      bsp->length);
      }
      seq_chain = seq_chain_start = bsp->seq_data->chain;
      sequence = seq_chain->str;
      BioseqUnlock(bsp);

      query_seq_combined[0] = 0x08;

      if (private_slp) {
	 query_seq_start = (Uint1Ptr) Malloc(((query_length)+2)*sizeof(Char));
	 query_seq = query_seq_start+1;

	 index=0;
	 /*
	 while ((residue=SeqPortGetResidue(spp)) != SEQPORT_EOF) {
	    if (IS_residue(residue)) {
	       query_seq[index] = residue;
	       index++;
	    }
	    }*/
	 if (seq_data_type == Seq_code_ncbi4na) {
	    for (index=0, index1=0; index<query_length/2; index++, index1++) {
	       if (index1 >= seq_chain->len_avail) {
		  seq_chain = seq_chain->next;
		  sequence = seq_chain->str;
		  index1 = 0;
	       }
	       query_seq[2*index] =
		  ncbi4na_to_blastna[(sequence[index1] & TOP_HALFBYTE)>>4];
	       query_seq[2*index+1] =
		  ncbi4na_to_blastna[sequence[index1] & BOTTOM_HALFBYTE];
	    }
	    if (query_length & 1)
	       query_seq[query_length-1] = 
		  ncbi4na_to_blastna[(sequence[index1] & TOP_HALFBYTE)>>4];
	 } else {/* if (seq_data_type == Seq_code_ncbi2na) */
	    rem = 3;
	    for (index=0, index1=0; index<query_length; index++, index1++) { 
	       if (index1 >= 4*seq_chain->len_avail) {
		  seq_chain = seq_chain->next;
		  sequence = seq_chain->str;
		  index1 = 0;
	       }
	       query_seq[index] = READDB_UNPACK_BASE_N(sequence[index1/4], rem);
	       if (rem>0) rem--;
	       else rem = 3;
	    }
	 }


	 query_seq_start[0] = query_seq[query_length] = 0x08;
	 if (mask_slpp && mask_slpp[context/2] != NULL)
	    MegaBlastMaskTheResidues(query_seq, query_length,
				     15, mask_slpp[context/2], FALSE,
				     SeqLocStart(private_slp), mask_at_hash);
	 if (filter_slp) 
	    MegaBlastMaskTheResidues(query_seq, query_length,
				     15, filter_slp, FALSE,
				     SeqLocStart(private_slp), mask_at_hash);
	 MemCpy(&query_seq_combined[1], query_seq, query_length+1);
	 full_query_length = query_length + 1;
      }

      if (private_slp_rev) {
	 query_seq_start_rev = (Uint1Ptr)
	    Malloc(((query_length)+2)*sizeof(Char));
	 query_seq_rev = query_seq_start_rev+1;
	 index=0;
	 /*
	 while ((residue=SeqPortGetResidue(spp_reverse)) != SEQPORT_EOF) {
	    if (IS_residue(residue)) {
	       query_seq_rev[index] = residue;
	       index++;
	    }
	    }*/
	 seq_chain = seq_chain_start;
	 sequence = seq_chain->str;
	 if (seq_data_type == Seq_code_ncbi4na) {
	    for (index=0, index1=0; index<query_length/2; index++, index1++) {
	       if (index1 >= seq_chain->len_avail) {
		  seq_chain = seq_chain->next;
		  sequence = seq_chain->str;
		  index1 = 0;
	       }
	       query_seq_rev[query_length-2*(index+1)] =
		  ncbi4na_to_blastna[inverse_halfbyte[sequence[index1] & BOTTOM_HALFBYTE]];
	       query_seq_rev[query_length-2*index-1] =
		  ncbi4na_to_blastna[inverse_halfbyte[(sequence[index1] & TOP_HALFBYTE)>>4]];
	    }
	    if (query_length & 1)
	       query_seq_rev[0] = 
		  ncbi4na_to_blastna[inverse_halfbyte[(sequence[index1] & TOP_HALFBYTE)>>4]];
	 } else {
	    rem = 3;
	    for (index=0, index1=0; index<query_length; index++, index1++) {
	       if (index1 >= 4*seq_chain->len_avail) {
		  seq_chain = seq_chain->next;
		  sequence = seq_chain->str;
		  index1 = 0;
	       }
	       query_seq_rev[query_length-index-1] = 3 -
		  READDB_UNPACK_BASE_N(sequence[index1/4], rem);
	       if (rem>0) rem--;
	       else rem = 3;
	    }
	 }

	 query_seq_start_rev[0] = query_seq_rev[query_length] = 0x08;
	 if (mask_slpp && mask_slpp[context/2] != NULL)
	    MegaBlastMaskTheResidues(query_seq_rev, query_length,
				     15, mask_slpp[context/2], TRUE,
				     query_length -
				     SeqLocStop(private_slp_rev) - 1, mask_at_hash);
	 if (filter_slp)	
	    MegaBlastMaskTheResidues(query_seq_rev, query_length, 15,
				     filter_slp, TRUE, query_length -
				     SeqLocStop(private_slp_rev) - 1, mask_at_hash);
	 
	 MemCpy(query_seq_combined+full_query_length+1, 
		query_seq_rev,query_length+1);
	 full_query_length += query_length + 1;
      }
      
      if (!search->query_id) {
	 search->query_id = SeqIdDup(SeqIdFindBest(SeqLocId(tmp_slp), SEQID_GI));
	 query_id = search->query_id;
      } else {
	 query_id->next =  SeqIdDup(SeqIdFindBest(SeqLocId(tmp_slp), SEQID_GI));
	 query_id = query_id->next;
      }

      MegaBlastSequenceAddSequence(search->context[search->first_context].query,
				   NULL, query_seq_combined, full_query_length,
				   full_query_length, 0);
      
      query_seq_combined = MemFree(query_seq_combined);
      if (context==search->first_context)
	 query_seq_start = MemFree(query_seq_start);
      
      if (private_slp && context != search->first_context) 
	 MegaBlastSequenceAddSequence(search->context[context].query, query_seq,
				      query_seq_start, query_length, 
				      query_length, 0);
      if (context % 2 == 0)
	 context++;
      /*spp = SeqPortFree(spp);*/
      if (private_slp_rev && context != search->first_context) {
	 MegaBlastSequenceAddSequence(search->context[context].query, query_seq_rev,
				      query_seq_start_rev, query_length,
				      query_length, 0);
	 /*spp_reverse = SeqPortFree(spp_reverse);*/
      }
      context++;

      if (mask_at_hash)
	 /* No longer needed. */
	 filter_slp = SeqLocSetFree(filter_slp);
      
      if (private_slp)
	 private_slp = SeqLocFree(private_slp);
      if (private_slp_rev)
	 private_slp_rev = SeqLocFree(private_slp_rev);
      slp = slp->next;
   } /* End of loop over query contexts (strands) */
   
   for (qid=search->query_id, index=0; qid && index<num_queries; 
	qid=qid->next, index++)
      search->qid_array[index] = qid;

   if (search->pbp->is_neighboring) 
      MemFree(homo_gilist);

   if (!search->pbp->one_line_results) 
      last_context = search->last_context;
   else
      last_context = search->first_context;
	for (index=search->first_context; index<=last_context; index++)
	{
	   length = search->context[index].query->length;

	   if (index>search->first_context || 
	       search->first_context==last_context) 
		status = BlastScoreBlkFill(search->sbp, (CharPtr)
					   search->context[index].query->sequence,search->context[index].query->length, index);
	   else {
	      length = search->context[index+1].query->length;
		status = BlastScoreBlkFill(search->sbp, (CharPtr)
					   search->context[index].query->sequence, length, index);
	   }
		if (status != 0 && index == search->first_context)
		{
			sprintf(buffer, "Unable to calculate Karlin-Altschul params, check query sequence");
			BlastConstructErrorMessage("BLASTSetUpSearch", buffer, 2, &(search->error_return));
			retval = 1;
		}
	}

	search->sbp->kbp_gap = search->sbp->kbp_gap_std;
        search->sbp->kbp = search->sbp->kbp_std;

	/* If retval was set non-zero above (by the routines calculating Karlin-Altschul params),
	   return here before these values are used.
	*/
	if (retval) {
           search->thr_info->awake_index = FALSE;
           return retval;
        }

	min_query_length = (Int4) 1/(search->sbp->kbp[search->first_context]->K);

	last_length_adjustment = 0;
	for (index=0; index<5; index++)
	{
		length_adjustment = (Int4) ((search->sbp->kbp[search->first_context]->logK)+log((Nlm_FloatHi)(length-last_length_adjustment)*(Nlm_FloatHi)(MAX(1, (search->dblen)-(search->dbseq_num*last_length_adjustment)))))/(search->sbp->kbp[search->first_context]->H);
		if (length_adjustment >= length-min_query_length)
		{
			length_adjustment = length-min_query_length;
			break;
		}
	
		if (ABS(last_length_adjustment-length_adjustment) <= 1)
			break;
		last_length_adjustment = length_adjustment;
	}
	search->length_adjustment = MAX(length_adjustment, 0);

	search->dblen_eff = MAX(1, search->dblen - search->dbseq_num*search->length_adjustment);
	effective_query_length = MAX(length - search->length_adjustment, min_query_length);
	
	for (index=search->first_context; index<=last_context; index++)
	{
		search->context[index].query->effective_length = effective_query_length;
	}

	if (search->searchsp_eff == 0)
		search->searchsp_eff = ((Nlm_FloatHi) search->dblen_eff)*((Nlm_FloatHi) effective_query_length);

	/* The default is that cutoff_s was not set and is zero. */
	if (options->cutoff_s == 0)
	{
		search->pbp->cutoff_e = options->expect_value;
		search->pbp->cutoff_e_set = TRUE;
		search->pbp->cutoff_s = options->cutoff_s;
		search->pbp->cutoff_s_set = FALSE;
	}
	else
	{
		search->pbp->cutoff_e = options->expect_value;
		search->pbp->cutoff_e_set = FALSE;
		search->pbp->cutoff_s = options->cutoff_s;
		search->pbp->cutoff_s_set = TRUE;
	}
/* For now e2 is set to 0.5 and cutoff_e2_set is FALSE.  This is then
changed to the proper values in blast_set_parameters.  In the final version
of this program (where more blast programs and command-line options are
available) this needs to be set higher up. */
	if (options->cutoff_s2 == 0)
	{
		search->pbp->cutoff_e2 = options->e2;
		search->pbp->cutoff_e2_set = FALSE;
		search->pbp->cutoff_s2 = options->cutoff_s2;
		search->pbp->cutoff_s2_set = FALSE;
	}
	else
	{
		search->pbp->cutoff_e2 = options->e2;
		search->pbp->cutoff_e2_set = FALSE;
		search->pbp->cutoff_s2 = options->cutoff_s2;
		search->pbp->cutoff_s2_set = TRUE;
	}
	
	search->pbp->discontinuous = options->discontinuous;

	
	/* For postion based blast. */
	search->pbp->ethresh = options->ethresh;
	search->pbp->maxNumPasses = options->maxNumPasses;
	search->pbp->pseudoCountConst = options->pseudoCountConst;

	search->pbp->process_num = options->number_of_cpus;
	search->pbp->cpu_limit = options->cpu_limit;
	search->pbp->gap_decay_rate = options->gap_decay_rate;
	search->pbp->gap_size = options->gap_size;
	search->pbp->gap_prob = options->gap_prob;
	search->pbp->old_stats = options->old_stats;
	search->pbp->use_large_gaps = options->use_large_gaps;
	search->pbp->number_of_bits = options->number_of_bits;
	search->pbp->two_pass_method = options->two_pass_method;
	search->pbp->multiple_hits_only = options->multiple_hits_only;
	search->pbp->gap_open = options->gap_open;
	search->pbp->gap_extend = options->gap_extend;
        search->pbp->decline_align = options->decline_align;

	search->pbp->hsp_num_max = options->hsp_num_max;
	/*search->pbp->gap_x_dropoff = (BLAST_Score)
	  (options->gap_x_dropoff*NCBIMATH_LN2 /
	  search->sbp->kbp[search->first_context]->Lambda);*/
	search->pbp->gap_x_dropoff = options->gap_x_dropoff;
	search->pbp->gap_x_dropoff_final = (BLAST_Score) (options->gap_x_dropoff_final*NCBIMATH_LN2 / search->sbp->kbp[search->first_context]->Lambda);
	search->pbp->gap_trigger = (BLAST_Score) ((options->gap_trigger*NCBIMATH_LN2+search->sbp->kbp[search->first_context]->logK)/ search->sbp->kbp[search->first_context]->Lambda);
	/* Ensures that gap_x_dropoff_final is at least as large as gap_x_dropoff. */
	search->pbp->gap_x_dropoff_final = MAX(search->pbp->gap_x_dropoff_final, search->pbp->gap_x_dropoff);

/* "threshold" (first and second) must be set manually for two-pass right now.*/
	search->pbp->threshold_set = TRUE;
	search->pbp->threshold_first = options->threshold_first;
	search->pbp->threshold_second = options->threshold_second;

	search->pbp->window_size = options->window_size;
	search->pbp->window_size_set = TRUE;

	search->whole_query = TRUE;
	if (options->required_start != 0 || options->required_end != -1)
	{
		search->whole_query = FALSE;
		search->required_start = options->required_start;
		if (options->required_end != -1)
		   search->required_end = options->required_end;
		else 
		   search->required_end = qlen;
	}

	/* Use DROPOFF_NUMBER_OF_BITS as the default if it's set to zero. */
	if (options->dropoff_1st_pass == 0)
		options->dropoff_1st_pass = (Int4) DROPOFF_NUMBER_OF_BITS;

	if (options->dropoff_2nd_pass == 0)
		options->dropoff_2nd_pass = (Int4) DROPOFF_NUMBER_OF_BITS;

        search->pbp->no_check_score = TRUE;

	avglen = BLAST_NT_AVGLEN;
	/* Use only one type of gap for blastn */
	search->pbp->ignore_small_gaps = TRUE;

	/* Only do this if this is not a pattern search. */
	if (options->isPatternSearch == FALSE)
	{
		search->wfp = search->wfp_second;
		MegaBlastBuildLookupTable(search);
	}
	/* 
	Turn off the index thread by setting this flag.  Don't wait for a join, as the
	search will take much longer than the one second for this to die.
	*/
	search->thr_info->awake_index = FALSE;

	return 0;
}

Boolean 
MegaBlastGetFirstAndLastContext(CharPtr prog_name, SeqLocPtr query_slp, Int2Ptr first_context, Int2Ptr last_context, Uint1 strand_options)
{
   SeqLocPtr tmp_slp = query_slp;
   Int2 tmp_first, tmp_last;
   Int2 index;

   if (StringCmp(prog_name, "blastn"))
      return BlastGetFirstAndLastContext(prog_name, query_slp, first_context, 
					 last_context, strand_options);
   
   if (!BlastGetFirstAndLastContext(prog_name, query_slp, first_context, 
				    last_context, strand_options))
      return FALSE;
   
   for (index=0; tmp_slp->next != NULL; index++, tmp_slp=tmp_slp->next);

   if (!BlastGetFirstAndLastContext(prog_name, tmp_slp, &tmp_first, 
				    &tmp_last, strand_options))
      return FALSE;
   /* For all intermediate sequences count 2 contexts */
   *last_context = 2 * index + tmp_last;
   return TRUE;
}


Int4 SeqLocTotalLen(CharPtr prog_name, SeqLocPtr slp)
{
   Int4 total_length = 0;
   SeqLocPtr tmp_slp = slp;

   if (StringCmp(prog_name, "blastn"))
      return SeqLocLen(slp);

   while (tmp_slp != NULL) {
      total_length += 2*(SeqLocLen(tmp_slp) + 1);
      tmp_slp = tmp_slp->next;
   }
   return total_length - 1;
}

BlastSearchBlkPtr 
BlastFillQueryOffsets(BlastSearchBlkPtr search, SeqLocPtr query_slp)
{
   SeqLocPtr slp = query_slp;
   Int2 num_contexts;
   Int2 i = 0;
   Int4 length;
   Uint1 strand;

   if (slp == NULL) return search;
   /* Note: the last offset will be equal to the combined query length + 1, for
      convenience of computing each individual query length later */
   for (slp=query_slp, num_contexts=0; slp; slp=slp->next, num_contexts+=2);
      

   search->query_context_offsets = (Int4Ptr) Malloc((num_contexts+1)*sizeof(Int4));
   
   search->query_context_offsets[0] = 0;
   for (slp = query_slp; slp; slp = slp->next) {
      length = SeqLocLen(slp) + 1;
      if (search->first_context == 0) 
	 search->query_context_offsets[i+1] = 
	    search->query_context_offsets[i] + length;
      else /* Top strand is not searched */
	 search->query_context_offsets[i+1] = 
	    search->query_context_offsets[i];
      i++;
      if (search->last_context % 2 == 1) 
	    search->query_context_offsets[i+1] = 
	       search->query_context_offsets[i] + length;
      else 
	 /* Bottom strand is not searched */
	 search->query_context_offsets[i+1] = 
	    search->query_context_offsets[i];
      i++;
   }	    
   

   return search;
}

static int LIBCALLBACK
evalue_compare_seqaligns(VoidPtr v1, VoidPtr v2)
{
   SeqAlignPtr sap1, sap2, PNTR sapp1, PNTR sapp2;
   SeqIdPtr id1, id2;

   sapp1 = (SeqAlignPtr PNTR) v1;
   sapp2 = (SeqAlignPtr PNTR) v2;
   
   sap1 = *sapp1;
   sap2 = *sapp2;
   
   if (sap1->score->value.realvalue > sap2->score->value.realvalue)
      return -1;
   else if (sap1->score->value.realvalue < sap2->score->value.realvalue)
      return 1;
   else if (sap1->score->value.intvalue < sap2->score->value.intvalue)
      return -1;
   else if (sap1->score->value.intvalue > sap2->score->value.intvalue)
      return 1;
   else 
      return 0;

}


SeqAlignPtr PNTR MegaBlastPackAlignmentsByQuery(BlastSearchBlkPtr search,
						SeqAlignPtr seqalign)
{
   SeqAlignPtr PNTR seqalignp, PNTR sapp;
   SeqAlignPtr seqalign_var, sap, seqalign_prev;
   Int2 index, i, hit_count;
   SeqIdPtr query_id, seg_query_id;
   SeqLocPtr slp;
   Int4 num_queries = search->last_context / 2 + 1;

   slp = search->query_slp;
   
   seqalignp = (SeqAlignPtr PNTR)
      MemNew(num_queries*sizeof(SeqAlignPtr));

   if (!seqalign) { 
      *seqalignp = NULL;
      return seqalignp;
   }

   for (index=0; index<num_queries; index++) {
      query_id = search->qid_array[index];
      hit_count = 0;
      seqalign_var = seqalign;
      seqalign_prev = NULL;
      while (seqalign_var != NULL) {
	 seg_query_id = ((DenseSegPtr)seqalign_var->segs)->ids;
	 if (SeqIdComp(query_id, seg_query_id) == SIC_YES) {
	    /* Remove this link from the list - we won't need it any more */
	    if (seqalign_prev != NULL)
	       seqalign_prev->next = seqalign_var->next;
	    else /* first link in the list */
	       seqalign = seqalign_var->next;
	    
	    seqalign_var->next = NULL;
	    
	    /* Add this seqalign to the list for this query */
	    if (seqalignp[index]==NULL) {
	       seqalignp[index] = seqalign_var;
	       sap = seqalignp[index];
	    } else {
	       sap->next = seqalign_var;
	       sap = sap->next;
	    }
	    hit_count++;
	    /* Move to the next link */
	    if (seqalign_prev != NULL)
	       seqalign_var = seqalign_prev->next;
	    else
	       seqalign_var = seqalign;
	 } else { /* hit from different query */
	    seqalign_prev = seqalign_var;
	    seqalign_var = seqalign_var->next;
	 }
      }
      
      /* Now sort seqaligns for this query based on evalue */
      if (hit_count>1) {
	 Int4 start, end;
	 SeqIdPtr sid;

	 sapp = Nlm_Malloc(hit_count*sizeof(SeqAlignPtr));
	 sapp[0] = seqalignp[index];

	 for (i=1; i<hit_count; i++)
	    sapp[i] = sapp[i-1]->next;

	 for (start=0; start<hit_count; ) {
	    end = start; 
	    sid = ((DenseSegPtr)sapp[start]->segs)->ids->next;
	    for(end=start+1; end < hit_count &&
		SeqIdComp(sid, ((DenseSegPtr)sapp[end]->segs)->ids->next)
		   ==SIC_YES;
		end++);

	    if (end > start)
	       HeapSort(&sapp[start], end-start, sizeof(SeqAlignPtr),
			evalue_compare_seqaligns);
	    start = end;
	 }
	 for (i=0; i<hit_count-1; i++)
	    sapp[i]->next = sapp[i+1];
	 sapp[hit_count-1]->next = NULL;
	 
	 seqalignp[index] = sapp[0];
	 MemFree(sapp);
      }
   }
   /* Shouldn't be necessary, but clean just in case */
   SeqAlignSetFree(seqalign);
   return seqalignp;
}

/* Attach the "sequence" pointer to the BlastSequenceBlkPtr. sequence_start may be the
actual start of the sequence (this pointer is kept for deallocation purposes).  The
sequence may start before "sequence" starts as there may be a sentinel (i.e., NULLB)
before the start of the sequence.  When the extension function extends this way it
can tell that there is a NULLB there and stop the extension.

*/

Int2 LIBCALL
MegaBlastSequenceAddSequence (BlastSequenceBlkPtr sequence_blk, Uint1Ptr sequence, 
			      Uint1Ptr sequence_start, Int4 length, 
			      Int4 original_length, Int4 effective_length)
{
   Uint1Ptr seqptr;

   if (sequence_blk == NULL)
      return 1;
   /* Assume that memory for sequence_blk->sequence is
      allocated elsewhere */
   if (sequence == NULL && sequence_start != NULL) {
      if (sequence_blk->sequence!=NULL && sequence_blk->length>0) {
	 seqptr = sequence_blk->sequence_start + sequence_blk->length;
	 MemCpy(seqptr, sequence_start, length+1);
      } else {
	 MemCpy(sequence_blk->sequence_start, sequence_start, length+1);
	 sequence_blk->sequence = sequence_blk->sequence_start + 1;
	 sequence_blk->effective_length = effective_length;
      }
   }
   else if (sequence != NULL) {
      sequence_blk->sequence = sequence;
      sequence_blk->sequence_start = sequence_start;
   }

   sequence_blk->length += length;
   sequence_blk->original_length += original_length;
   
   return 0;
}

#ifndef BUFFER_LENGTH
#define BUFFER_LENGTH 255
#endif
#define MAX_LINE 48

Int4 
MegaBlastWordFinder(BlastSearchBlkPtr search, LookupTablePtr lookup)
{
   register Uint1Ptr subject;
   register Int4 s_off, ecode, mask, q_off;
   Int4 index;
   Int4 subj_length = search->subject->length;

   if (search->pbp->is_neighboring && 
       subj_length < MIN_NEIGHBOR_HSP_LENGTH)
      return 0;
   if (search->current_hitlist == NULL)
      search->current_hitlist = BlastHitListNew(search);
   
   mask = lookup->mb_lt->mask;
   subject = search->subject->sequence - 1;
   ecode = 0;
   lookup->mb_lt->stack_index = 0;

   for (s_off = 0; s_off < (lookup->mb_lt->width - 1)*4; s_off += 4) {
      ecode = (ecode << 8) + *++subject;
   }
   s_off += 4;
   while (s_off<subj_length) {
      ecode = ((ecode & mask) << 8) + *++subject;
      for (q_off = lookup->mb_lt->hashtable[ecode]; q_off>0; 
	   q_off = lookup->mb_lt->next_pos[q_off]) 
	 MegaBlastExtendHit(search, lookup, s_off, q_off);
      s_off += 4;
   }

   /* Save hits remaining on the stack */
   for (index=0; index<lookup->mb_lt->stack_index; index++) {
      q_off = lookup->mb_lt->estack[index].level - 
	 lookup->mb_lt->estack[index].diag;
      s_off = lookup->mb_lt->estack[index].level;
      if (lookup->mb_lt->estack[index].length >= lookup->mb_lt->lpm)
	 BlastSaveCurrentHsp(search, 0, q_off, s_off, 
			     lookup->mb_lt->estack[index].length, 0);
   }

   /* Do greedy gapped extension */
   if (search->current_hitlist->hspcnt > 0)
      MegaBlastGappedAlign(search);

   if (search->current_hitlist)
      return search->current_hitlist->hspcnt;
   else return 0;
}

static int LIBCALLBACK
diag_compare_hsps(VoidPtr v1, VoidPtr v2)
{
   BLAST_HSPPtr h1, h2;

   h1 = *((BLAST_HSPPtr PNTR) v1);
   h2 = *((BLAST_HSPPtr PNTR) v2);
   
   return SIGN((h1->query.offset - h1->subject.offset) - 
      (h2->query.offset - h2->subject.offset));
}

Int4 MegaBlastGappedAlign(BlastSearchBlkPtr search)
{
   BLAST_HSPPtr hsp, PNTR hsp_array;
   Int4 index, i, hspcnt, new_hspcnt=0, buf_len=0;
   Uint1Ptr subject0 = NULL, seq_4na, seq_2na;
   Uint1 rem;
   Boolean delete_hsp;

   hspcnt = search->current_hitlist->hspcnt;
   search->second_pass_extends += hspcnt;
   hsp_array = search->current_hitlist->hsp_array;
   /* Make current hitlist available for rewriting of extended hsps 
      without freeing the hsp_array since it's used in this function */
   search->current_hitlist->hspcnt = 0;

   /* Get subject sequence in Seq_code_ncbi4na format with all 
      ambiguities information */
   if (search->rdfp)
      readdb_get_sequence_ex(search->rdfp, search->subject_id, 
			     &subject0, &buf_len, FALSE);
   else { /* Convert ncbi2na to ncbi4na encoding */
      subject0 = (Uint1Ptr) MemNew(search->subject->length + 1);
      seq_4na = subject0;
      seq_2na = search->subject->sequence;
      rem = 3;
      for (i=0; i<search->subject->length; i++) {
	 *seq_4na = (Uint1) (1 << READDB_UNPACK_BASE_N(*seq_2na, rem));
	 seq_4na++;
	 if (rem>0) rem--;
	 else {
	    rem = 3;
	    seq_2na++;
	 }
      }
   }
   
   HeapSort(hsp_array, hspcnt, sizeof(BLAST_HSPPtr), diag_compare_hsps);

   for (index=0; index<hspcnt; index++) {
      hsp = hsp_array[index];
      delete_hsp = FALSE;
      for (i = search->current_hitlist->hspcnt-1; 
	   i >= 0 && MB_HSP_CLOSE(hsp->query.offset, search->current_hitlist->hsp_array[i]->query.offset, hsp->subject.offset, search->current_hitlist->hsp_array[i]->subject.offset, MB_DIAG_NEAR);
	   i--) {
	if (MB_HSP_CONTAINED(hsp->query.offset, search->current_hitlist->hsp_array[i]->query.offset, search->current_hitlist->hsp_array[i]->query.end, hsp->subject.offset, search->current_hitlist->hsp_array[i]->subject.offset, search->current_hitlist->hsp_array[i]->subject.end, MB_DIAG_CLOSE)) {
	  delete_hsp = TRUE;
	  break;
	}
      }
      if (!delete_hsp) 
	MegaBlastNtWordExtend(search, subject0, hsp->query.offset, hsp->subject.offset);
      MemFree(hsp);
   }
   MemFree(subject0);
   return 0;
}

Int4
MegaBlastExtendHit(BlastSearchBlkPtr search, LookupTablePtr lookup, 
		   Int4 s_off, Int4 q_off) {
   Int4 index, len = 4*lookup->mb_lt->width; 
   Int4 query_offset;
   Int2 skip, min_hit_length;


   if (lookup->mb_lt->lpm >= 12) {
      skip = 0;
      min_hit_length = lookup->mb_lt->lpm;
   } else {
      skip = 12;
      min_hit_length = lookup->mb_lt->lpm + 4;
   }
   
   for (index=0; index<lookup->mb_lt->stack_index; ) {
      if (lookup->mb_lt->estack[index].diag == s_off - q_off) {
	 if (lookup->mb_lt->estack[index].level < s_off - 4 - skip) {
	    /* That hit doesn't go further, save it and substitute by
	       the new one */
	    if (lookup->mb_lt->estack[index].length >= min_hit_length) {
	       query_offset = lookup->mb_lt->estack[index].level -
		  lookup->mb_lt->estack[index].diag;
	       BlastSaveCurrentHsp(search, 0, query_offset, 
				   lookup->mb_lt->estack[index].level, 
				   lookup->mb_lt->estack[index].length, 0);
	    }
	    lookup->mb_lt->estack[index].length = len;
	    lookup->mb_lt->estack[index].level = s_off;
	 } else if (lookup->mb_lt->estack[index].level < s_off - 4) {
	    /* We had a hit on this diagonal within a window */
	    lookup->mb_lt->estack[index].length += len;
	    lookup->mb_lt->estack[index].level = s_off;
	 } else {
	    /* We had a hit on this diagonal which is extended by current hit */
	    lookup->mb_lt->estack[index].length += 4;
	    lookup->mb_lt->estack[index].level = s_off;
	 }
	 return 0;   
      } else if (lookup->mb_lt->estack[index].level >= s_off - 4 - skip) 
	 /* Just skip - it's on a different diagonal */
	 index++;
      else { 
	 /* Hit from different diagonal, and it is finished, so save it 
	    and remove from stack */
	 lookup->mb_lt->stack_index--;
	 if (lookup->mb_lt->estack[index].length >= min_hit_length) {
	    query_offset = lookup->mb_lt->estack[index].level -
	       lookup->mb_lt->estack[index].diag;
	    BlastSaveCurrentHsp(search, 0, query_offset, 
				lookup->mb_lt->estack[index].level, 
				lookup->mb_lt->estack[index].length, 0);
	 }

	 lookup->mb_lt->estack[index] = 
	    lookup->mb_lt->estack[lookup->mb_lt->stack_index];
      }
   }
   /* Need an extra slot on the stack for this hit */
   if (lookup->mb_lt->stack_index>=lookup->mb_lt->stack_size) {
      /* Stack about to overflow - reallocate memory */
      MbStackPtr ptr;
      if (!(ptr = (MbStackPtr)Realloc(lookup->mb_lt->estack, 
				      2*lookup->mb_lt->stack_size*sizeof(MbStack)))) {
	 ErrPostEx(SEV_WARNING, 0, 0, "Unable to allocate %ld extra stack spaces",
		   2*lookup->mb_lt->stack_size);
	 return 1; 
      } else {
	 lookup->mb_lt->stack_size *= 2;
	 lookup->mb_lt->estack = ptr;
      }
   }
   lookup->mb_lt->estack[lookup->mb_lt->stack_index].diag = s_off - q_off;
   lookup->mb_lt->estack[lookup->mb_lt->stack_index].level = s_off;
   lookup->mb_lt->estack[lookup->mb_lt->stack_index].length = len;
   lookup->mb_lt->stack_index++;
   
   return 0;
}
 
Int2
MegaBlastNtWordExtend(BlastSearchBlkPtr search, Uint1Ptr subject0, 
		      Int4 q_off, Int4 s_off)
{
	register Uint1Ptr	q;
	register BLAST_Score    score;
	Uint1Ptr query0, s;
	BLAST_Score	X;
        BLAST_ParameterBlkPtr   pbp;
        BLAST_ScoreBlkPtr       sbp;
	Int4 q_avail, s_avail;
	/* Variables for the call to MegaBlastGreedyAlign */
	Int4 ndiff, q_ext_l, q_ext_r, s_ext_l, s_ext_r;
	edit_script_t *ed_script_fwd=NULL, *ed_script_rev=NULL;
	Boolean good_hit;
	
	good_hit = TRUE;
        sbp=search->sbp;
        pbp=search->pbp;

	query0 = (Uint1Ptr) search->context[search->first_context].query->sequence;
        q_avail = search->context[search->first_context].query->length - q_off;
        s_avail = search->subject->length - s_off;

	q = query0 + q_off;
	s = subject0 + s_off;

	/* Find where positive scoring starts & ends within the word hit */
	score = 0;

	X = pbp->gap_x_dropoff;

	if (!search->pbp->one_line_results) {
	   ed_script_fwd = edit_script_new();
	   ed_script_rev = edit_script_new();
	}

	
	/* extend to the right */
	score = MegaBlastAffineGreedyAlign(s, s_avail, q, q_avail, FALSE, X,
					   sbp->reward, -sbp->penalty,
					   pbp->gap_open, pbp->gap_extend,
					   &s_ext_r, &q_ext_r, search->abmp, 
					   ed_script_fwd);
	
	/* extend to the left */
	score += MegaBlastAffineGreedyAlign(subject0, s_off, query0, q_off, 
					    TRUE, X, sbp->reward, 
					    -sbp->penalty, pbp->gap_open, 
					    pbp->gap_extend, &s_ext_l, 
					    &q_ext_l, search->abmp, 
					    ed_script_rev);
	
	/* For neighboring we have a stricter criterion to keep an HSP */
	if (search->pbp->is_neighboring) {
	   FloatHi perc_identity;
	   Int4 hsp_length;
	   
	   hsp_length = MIN(q_ext_l+q_ext_r, s_ext_l+s_ext_r);
	   perc_identity = 100 * (1 - ((FloatHi) score) / hsp_length);
	   if (hsp_length < MIN_NEIGHBOR_HSP_LENGTH || 
	       perc_identity < MIN_NEIGHBOR_PERC_IDENTITY)
	      good_hit = FALSE;
	}

	/* In basic case the greedy algorithm returns number of 
	   differences, hence we need to convert it to score */
	if (pbp->gap_open==0 && pbp->gap_extend==0) 
	   score = 
	      (q_ext_r + s_ext_r + q_ext_l + s_ext_l)*sbp->reward/2 - 
	      score*(sbp->reward - sbp->penalty);
	else if (sbp->reward % 2 == 1)
	   score /= 2;

	if (good_hit && score >= pbp->cutoff_s2) { /* Score is reportable */
	   if (!search->pbp->one_line_results) 
	      edit_script_append(ed_script_rev, ed_script_fwd);
	   
	   search->second_pass_good_extends++;
           
	   BlastSaveCurrentHspGapped(search, score, q_off-q_ext_l,
				     s_off-s_ext_l, q_ext_l+q_ext_r,
				     s_ext_l+s_ext_r, search->first_context,
				     ed_script_rev); 
	}
        edit_script_free(ed_script_fwd);
        edit_script_free(ed_script_rev);

	return 0;
}

/* A routine to mask the residues from a mask SeqLoc by either converting them 
 * to lowercase or changing to a given value.
 * The buffer should contain an already decoded sequence. 
 */
void MegaBlastMaskTheResidues(Uint1Ptr buffer, Int4 max_length, Uint1
			      mask_residue, SeqLocPtr mask_slp, Boolean reverse,
			      Int4 offset, Boolean mask_at_hash)

{
   SeqLocPtr slp=NULL;
   Int4 index, start, stop;
   
   while (mask_slp) {
      slp=NULL;
      while((slp = SeqLocFindNext(mask_slp, slp))!=NULL) {
	 if (reverse) {
	    start = max_length - 1 - SeqLocStop(slp);
	    stop = max_length - 1 - SeqLocStart(slp);
	 } else {
	    start = SeqLocStart(slp);
	    stop = SeqLocStop(slp);
	 }
	 
	 start -= offset;
	 stop  -= offset;
	 
	 if (mask_at_hash) 
	    for (index=start; index<=stop; index++)
	       /* Set the 3rd bit to mark that this is a masked residue */
	       buffer[index] = buffer[index] | 0x04;
	 else
	    for (index=start; index<=stop; index++)
	       buffer[index] = mask_residue;
      }
      mask_slp = mask_slp->next;
   }
}

/* The following function can only be used for IUPAC encoding */
void
BlastLCaseMaskTheResidues(Uint1Ptr buffer, Int4 max_length, SeqLocPtr mask_slp, Boolean reverse, Int4 offset)

{
   SeqLocPtr slp=NULL;
   Int4 index, start, stop;
   
   while (mask_slp) {
      slp=NULL;
      while((slp = SeqLocFindNext(mask_slp, slp))!=NULL) {
	 if (reverse) {
	    start = max_length - 1 - SeqLocStop(slp);
	    stop = max_length - 1 - SeqLocStart(slp);
	 }	else {
	    start = SeqLocStart(slp);
	    stop = SeqLocStop(slp);
	 }
	 
	 start -= offset;
	 stop  -= offset;
	 
	 for (index=start; index<=stop; index++)
	    buffer[index] = tolower(buffer[index]);
      }
      mask_slp = mask_slp->next;
   }
}

BLAST_WordFinderPtr
MegaBlastWordFinderDeallocate(BLAST_WordFinderPtr wfp)
{
   LookupTablePtr lookup;
   
   if (wfp == NULL || (lookup = wfp->lookup) == NULL || 
       lookup->mb_lt == NULL) 
      return wfp;
   if (lookup->mb_lt->estack)
      MemFree(lookup->mb_lt->estack);
   MemFree(lookup->mb_lt);
   MemFree(lookup);
   wfp = MemFree(wfp);
   return wfp;
}

static int LIBCALLBACK
fwd_compare_hsps(VoidPtr v1, VoidPtr v2)

{
	BLAST_HSPPtr h1, h2;
	BLAST_HSPPtr PNTR hp1, PNTR hp2;

	hp1 = (BLAST_HSPPtr PNTR) v1;
	hp2 = (BLAST_HSPPtr PNTR) v2;
	h1 = *hp1;
	h2 = *hp2;
	
	if (h1==NULL && h2==NULL) return 0;
	else if (h1==NULL) return 1;
	else if (h2==NULL) return -1;

	/* If the two HSP's have same coordinates, they are equal */
	if (h1->query.offset == h2->query.offset && 
	    h1->query.end == h2->query.end && 
	    h1->subject.offset == h2->subject.offset &&
	    h1->subject.end == h2->subject.end)
	   return 0;

	/* Check if one HSP is contained in the other, if so, 
	   leave only the longer one */
	if (h1->query.offset >= h2->query.offset && 
	    h1->query.end <= h2->query.end &&  
	    h1->subject.offset >= h2->subject.offset && 
	    h1->subject.end <= h2->subject.end) { 
	   GapXEditBlockDelete(h1->gap_info); 
	   *hp1 = MemFree(*hp1);
	   return 1; 
	} else if (h1->query.offset <= h2->query.offset &&  
		   h1->query.end >= h2->query.end &&  
		   h1->subject.offset <= h2->subject.offset && 
		   h1->subject.end >= h2->subject.end) { 
	   GapXEditBlockDelete(h2->gap_info); 
	   *hp2 = MemFree(*hp2);
	   return -1; 
	}

	if (h1->query.offset < h2->query.offset) 
		return -1;
	if (h1->query.offset > h2->query.offset) 
		return 1;
	/* Necessary in case both HSP's have the same query offset. */
	if (h1->subject.offset < h2->subject.offset) 
		return -1;
	if (h1->subject.offset > h2->subject.offset) 
		return 1;

	return 0;
}

void
BlastSortUniqHspArray(BLAST_HitListPtr hitlist)
{
   Int4 index, new_hspcnt;
   BLAST_HSPPtr PNTR hsp_array = hitlist->hsp_array;

   HeapSort(hitlist->hsp_array, hitlist->hspcnt, sizeof(BLAST_HSPPtr), fwd_compare_hsps);
   for (index=1, new_hspcnt=0; index<hitlist->hspcnt; index++) {
      if (hsp_array[index]==NULL) 
	 continue;
      if (fwd_compare_hsps(&hsp_array[new_hspcnt], &hsp_array[index])) {
	 if (hsp_array[index] != NULL) {
	    new_hspcnt++;
	    if (index > new_hspcnt)
	       hsp_array[new_hspcnt] = hsp_array[index];
	 }
      } else {
	 GapXEditBlockDelete(hsp_array[index]->gap_info);
	 hsp_array[index] = MemFree(hsp_array[index]);
      }
   }
   /* Set all HSP pointers that will not be used any more to NULL */
   MemSet(&hsp_array[new_hspcnt+1], 0, 
	  (hitlist->hspcnt - new_hspcnt - 1)*sizeof(BLAST_HSPPtr));
   hitlist->hspcnt = new_hspcnt + 1;
   hitlist->hspcnt_max = hitlist->hspcnt;
}

void
MegaBlastFillHspGapInfo(BLAST_HSPPtr hsp, edit_script_t PNTR ed_script)
{
   hsp->gap_info = GapXEditBlockNew(hsp->query.offset, hsp->subject.offset);

   hsp->gap_info->start1 = hsp->query.offset;
   hsp->gap_info->start2 = hsp->subject.offset;
   hsp->gap_info->length1 = hsp->query.length;
   hsp->gap_info->length2 = hsp->subject.length;
   hsp->gap_info->frame1 = hsp->gap_info->frame2 = 1;
   hsp->gap_info->reverse = 0;
   
   hsp->gap_info->esp = MBToGapXEditScript(ed_script);
}

static Uint4 opc_trans[4] = {0,2,1,3};

GapXEditScriptPtr
MBToGapXEditScript (edit_script_t PNTR ed_script)
{
   GapXEditScriptPtr esp_start, esp, esp_prev;
   Int4 i;

   if (ed_script==NULL || ed_script->num == 0)
      return NULL;

   for (i=0; i<ed_script->num; i++) {
      /*esp_var = GapXEditScriptNew(esp_var);*/
      esp = (GapXEditScriptPtr) MemNew(sizeof(GapXEditScript));
      esp->num = EDIT_VAL(ed_script->op[i]); 
      esp->op_type = opc_trans[EDIT_OPC(ed_script->op[i])];
      if (i==0) 
	 esp_start = esp_prev = esp;
      else {
	 esp_prev->next = esp;
	 esp_prev = esp;
      }
   }

   return esp_start;

}

SeqAlignPtr
MegaBlastSeqAlignFromResultHitlist(BlastSearchBlkPtr search,
				   BLASTResultHitlistPtr result_hitlist,
				   SeqIdPtr subject_id)
{
   SeqAlignPtr seqalign, seqalign_var, sap;
   BLASTResultHspPtr hsp_array;
   Int4 index, i;
   SeqIdPtr query_id, new_subject_seqid;
   StdSegPtr ssp;
   ValNodePtr gi_list=NULL;

   hsp_array = result_hitlist->hsp_array;

   gi_list = BlastGetAllowedGis(search, result_hitlist->subject_id, &new_subject_seqid);
   
   for (index=0; index<result_hitlist->hspcnt; index++) { 
      query_id = search->qid_array[hsp_array[index].context/2];
      if (index==0) {
	 if (new_subject_seqid)
	    seqalign = seqalign_var =
	       GapXEditBlockToSeqAlign(hsp_array[index].gap_info, 
				       new_subject_seqid, query_id); 
	 else
	    seqalign = seqalign_var =
	       GapXEditBlockToSeqAlign(hsp_array[index].gap_info, 
				       subject_id, query_id);
      } else {
	 if (new_subject_seqid)
	    seqalign_var->next = 
	       GapXEditBlockToSeqAlign(hsp_array[index].gap_info, 
				       new_subject_seqid, query_id); 
	 else
	    seqalign_var->next = 
	       GapXEditBlockToSeqAlign(hsp_array[index].gap_info, 
				       subject_id, query_id); 
	    
	 seqalign_var = seqalign_var->next;
      }
      seqalign_var->score = GetScoreSetFromBlastResultHsp(&hsp_array[index], gi_list);
      if (seqalign_var->segtype == 3) {
	 ssp = seqalign_var->segs;
	 while (ssp) {
	    ssp->scores = GetScoreSetFromBlastResultHsp(&hsp_array[index], gi_list);
	    ssp = ssp->next;
	 }
      } else if (seqalign_var->segtype == 5) { /* Discontinuous */
	 sap = (SeqAlignPtr) seqalign_var->segs;
	 for(;sap != NULL; sap = sap->next) {
	    sap->score = GetScoreSetFromBlastResultHsp(&hsp_array[index], gi_list);
	 }
      }
   }
   return seqalign;
}
/* Function below assumes that the gi file is sorted in increasing order */
Boolean 
BlastNeedHumanRepeatFiltering(BlastDoubleInt4Ptr gi_list, 
			      Int4 gi_list_size, Uint4 query_gi)
{
   Int4 begin, end, middle;
   
   begin = 0; 
   end = gi_list_size - 1;
   while (begin <= end) {
      middle = (begin + end) / 2;
      if (gi_list[middle].gi > query_gi)
	 end = middle - 1;
      else if (gi_list[middle].gi < query_gi)
	 begin = middle + 1;
      else {
	 return TRUE;
      }
   }
   return FALSE;
}

/* Subject sequence is assumed to be packed */
FloatHi
MegaBlastGetPercentIdentity(Uint1Ptr query, Uint1Ptr subject, Int4 q_start, 
			    Int4 s_start, Int4 length, Boolean reverse)
{
   Int4 i, ident = 0;
   FloatHi perc_ident;
   Uint1Ptr q, s;
   Int2 rem;

   q = &query[q_start];
   if (!reverse) {
      s = &subject[s_start / READDB_COMPRESSION_RATIO];
      rem = 3 - s_start % READDB_COMPRESSION_RATIO;
   } else {
      s = &subject[(s_start + length - 1) / READDB_COMPRESSION_RATIO];
      rem = 3 - (s_start + length - 1) % READDB_COMPRESSION_RATIO;
   }

   for (i=0; i<length; i++) {
      if (*q == READDB_UNPACK_BASE_N(*s, rem))
	 ident++;
      q++;
      if (!reverse) {
	 if (rem > 0)
	    rem--;
	 else {
	    rem = READDB_COMPRESSION_RATIO - 1;
	    s++;
	 }
      } else {
	 if (rem < READDB_COMPRESSION_RATIO - 1)
	    rem++;
	 else {
	    rem = 0;
	    s--;
	 }
      }
   }
   perc_ident = ((FloatHi) ident) / length * 100.0;

   return perc_ident;
}

/* The following binary search routine assumes that array A is filled. */
Int4 BinarySearchInt4(Int4 n, Int4Ptr A, Int4 size)
{
    Int4 m, b, e;

    b = 0;
    e = size;
    while (b < e - 1) {
	m = (b + e) / 2;
	if (A[m] > n)
	    e = m;
	else
	    b = m;
    }
    return b;
}

SeqLocPtr MaskSeqLocFromSeqAlign(SeqAlignPtr seqalign)
{
   SeqLocPtr slp = NULL, new_slp, slp_var;
   SeqAlignPtr sap;
   DenseSegPtr dsp;
   Int4 from, to;

   for (sap = seqalign; sap; sap = sap->next) {
      dsp = (DenseSegPtr) sap->segs;
      from = dsp->starts[0];
      to = dsp->starts[2*(dsp->numseg - 1)] + dsp->lens[dsp->numseg - 1] - 1;
      new_slp = SeqLocIntNew(from, to, Seq_strand_plus, dsp->ids);
      if (slp==NULL)
	 slp = slp_var = new_slp;
      else {
	 slp_var->next = new_slp;
	 slp_var = slp_var->next;
      }
   }
   
   return slp;
}

void PrintMaskedSequence(BioseqPtr query_bsp, SeqLocPtr mask_slp, 
				CharPtr file_name)
{
   SeqLocPtr query_slp = NULL;
   SeqPortPtr spp;
   FILE *fp;
   Uint1Ptr query_seq, seq_line;
   BioseqPtr masked_bsp;
   ByteStorePtr byte_sp;
   Int4 index, id_length, defline_len;
   Uint1 residue;
   CharPtr buffer = NULL;

   ValNodeAddPointer(&query_slp, SEQLOC_WHOLE, SeqIdDup(SeqIdFindBest(query_bsp->id, SEQID_GI)));
   

   spp = SeqPortNewByLoc(query_slp, Seq_code_iupacna);
   SeqPortSet_do_virtual(spp, TRUE);
   query_seq = (Uint1Ptr) MemNew((query_bsp->length + 1)*sizeof(Uint1));
   index=0;
   while ((residue=SeqPortGetResidue(spp)) != SEQPORT_EOF) {
      if (IS_residue(residue)) {
	 query_seq[index] = residue;
	 index++;
      }
   }
   BlastMaskTheResidues(query_seq, query_bsp->length, 'N', mask_slp, FALSE,
			     SeqLocStart(query_slp));

   fp = fopen(file_name, "w");
   if (query_bsp->descr)
      defline_len = StringLen(BioseqGetTitle(query_bsp));
   else
      defline_len = 0;
   defline_len += 255;     /* Sufficient for an ID. */
   buffer = MemNew((defline_len+1)*sizeof(Char));
   SeqIdWrite(query_bsp->id, buffer, PRINTID_FASTA_LONG, BUFFER_LENGTH);
   id_length = StringLen(buffer);
   buffer[id_length] = ' ';
   id_length++;
   StringCpy(&buffer[id_length], BioseqGetTitle(query_bsp));
   fprintf(fp, ">%s\n", buffer);
   
   seq_line = query_seq;
   for (index=0; index<=(query_bsp->length-1)/60; index++) {
      fprintf(fp, "%.60s\n", seq_line);
      seq_line += 60;
   }
      
   FileClose(fp);
   spp = SeqPortFree(spp);
   query_seq = MemFree(query_seq);

}
