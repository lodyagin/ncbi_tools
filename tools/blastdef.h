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

File name: blastdef.h

Author: Tom Madden

Contents: #defines and definitions for structures used by BLAST.

******************************************************************************/

/* $Revision: 6.43 $ */
/* $Log: blastdef.h,v $
/* Revision 6.43  1999/01/05 13:57:19  madden
/* Changed version and release date
/*
 * Revision 6.42  1998/12/31 18:17:03  madden
 * Added strand option
 *
 * Revision 6.41  1998/12/29 17:45:06  madden
 * Add do_sum_stats flag
 *
 * Revision 6.40  1998/12/21 13:09:53  madden
 * Changed version and release date
 *
 * Revision 6.39  1998/11/04 01:36:05  egorov
 * Add support for entrez-query and org-name to blast3
 *
 * Revision 6.38  1998/09/16 18:58:57  madden
 * Changed release number and date
 *
 * Revision 6.37  1998/09/14 15:11:15  egorov
 * Add support for Int8 length databases; remove unused variables
 *
 * Revision 6.36  1998/07/30 19:00:32  madden
 * Change to allow search of subset of database
 *
 * Revision 6.35  1998/07/28 21:17:59  madden
 * Added do_not_reevaluate
 *
 * Revision 6.34  1998/07/25 14:26:38  madden
 * Added comments
 *
 * Revision 6.33  1998/07/22 12:16:25  madden
 * Added handle_results
 *
 * Revision 6.32  1998/07/21 20:58:04  madden
 * Changes to allow masking at hash only
 *
 * Revision 6.31  1998/07/17 15:39:56  madden
 * Changes for Effective search space.
 *
 * Revision 6.30  1998/07/14 20:17:05  egorov
 * Add two new parameters (gilist and gifile) to BLAST_OptionsBlk
 *
 * Revision 6.29  1998/06/17 18:10:07  madden
 * Added isPatternSearch to Options
 *
 * Revision 6.28  1998/06/12 16:08:49  madden
 * BlastHitRange stuff
 *
 * Revision 6.27  1998/05/28 19:59:16  madden
 * Added typedef for BLASTHeapStruct
 *
 * Revision 6.26  1998/05/17 16:28:43  madden
 * Allow changes to filter options and cc filtering.
 *
 * Revision 6.25  1998/05/05 13:56:38  madden
 * Raised version to 2.0.5 and changed date
 *
 * Revision 6.24  1998/04/24 19:27:05  madden
 * Added BlastMatrixRescalePtr
 *
 * Revision 6.23  1998/04/01 22:47:14  madden
 * Added query_invalid flag
 *
 * Revision 6.22  1998/03/24 15:38:22  madden
 * Use BlastDoubleInt4Ptr to keep track of gis and ordinal_ids
 *
 * Revision 6.21  1998/03/18 14:14:20  madden
 * Support random access by gi list
 *
 * Revision 6.20  1998/03/14 18:29:21  madden
 * Added BlastSeqIdListPtr
 *
 * Revision 6.19  1998/02/26 22:34:37  madden
 * Changes for 16 bit windows
 *
 * Revision 6.18  1998/02/26 19:10:37  madden
 * Removed elements with BLAST_COLLECT_SPECIAL_STATS defines
 *
 * Revision 6.17  1998/02/24 22:46:29  madden
 * Added perform_culling Boolean and changed release date
 *
 * Revision 6.16  1998/02/19 17:17:10  madden
 * Use of Int4 rather than Int2 when pruning SeqAlign
 *
 * Revision 6.15  1998/01/05 16:46:52  madden
 * One or both strands can be searched, as opposed to only both, changes to number of contexts
 *
 * Revision 6.14  1997/12/23 19:14:14  madden
 * release version to 2.0.4
 *
 * Revision 6.13  1997/12/23 18:12:32  madden
 * Changes for range-dependent blast
 *
 * Revision 6.12  1997/12/12 20:38:02  madden
 * Fix to comments
 *
 * Revision 6.11  1997/12/11 22:20:16  madden
 * Corrected blast_type defines
 *
 * Revision 6.10  1997/12/10 22:41:40  madden
 * program number defines
 *
 * Revision 6.9  1997/11/14 21:30:16  madden
 * Changed version and date
 *
 * Revision 6.8  1997/10/26 17:26:59  madden
 * Changes for range dependent limits
 *
 * Revision 6.7  1997/10/01 13:35:28  madden
 * Changed BLAST_VERSION to BLAST_ENGINE_VERSION
 *
 * Revision 6.6  1997/09/22 17:36:24  madden
 * MACROS for position-specific matrices from Andy Neuwald
 *
 * Revision 6.5  1997/09/18 22:22:12  madden
 * Added prune functions
 *
 * Revision 6.4  1997/09/11 18:49:26  madden
 * Changes to enable searches against multiple databases.
 *
 * Revision 6.3  1997/09/10 21:27:57  madden
 * Changes to set CPU limits
 *
 * Revision 6.2  1997/09/03 19:06:35  madden
 * changed BLAST_VERSION and BLAST_RELEASE_DATE
 *
 * Revision 6.1  1997/08/27 14:46:48  madden
 * Changes to enable multiple DB searches
 *
 * Revision 6.0  1997/08/25 18:52:32  madden
 * Revision changed to 6.0
 *
 * Revision 1.63  1997/08/20 21:43:10  madden
 * Updated release date
 *
 * Revision 1.62  1997/07/21 17:37:15  madden
 * Added define for BLAST_RELEASE_DATE
 *
 * Revision 1.61  1997/07/18 20:55:45  madden
 * Added BLAST_VERSION
 *
 * Revision 1.60  1997/07/15 20:36:43  madden
 * Added ValNodePtr mask
 *
 * Revision 1.59  1997/07/14 15:33:00  madden
 * typedef for BlastErrorMsg
 *
 * Revision 1.58  1997/05/22 21:24:52  madden
 * Added support for final gapX dropoff value
 *
 * Revision 1.57  1997/05/20 17:51:33  madden
 * Added element SeqLocPtr query_slp to BlastSearch
 *
 * Revision 1.56  1997/05/06 22:19:35  madden
 * Added use_large_gaps and subject_length
 *
 * Revision 1.55  1997/04/09  20:01:53  madden
 * Added seqid_list to SearchBlk
 *
 * Revision 1.54  1997/04/03  19:48:13  madden
 * Changes to use effective database length instead of the length of each
 * sequence in statistical calculations.
 *
 * Revision 1.53  1997/03/31  17:07:57  madden
 * Added BLAST_COLLECT_STATS define.
 *
 * Revision 1.52  1997/03/20  22:56:24  madden
 * Added gap_info to hsp.
 *
 * Revision 1.51  1997/03/14  22:06:11  madden
 * fixed MT bug in BlastReevaluateWithAmbiguities.
 *
 * Revision 1.50  1997/03/08  16:52:16  madden
 * y
 * Added discontinuous option to ParameterBlk.
 *
 * Revision 1.49  1997/02/25  19:17:05  madden
 * Added discontinuous flag to options.
 *
 * Revision 1.48  1997/02/23  16:44:47  madden
 * GapAlignBlkPtr added to search structure.
 *
 * Revision 1.47  1997/02/20  18:38:34  madden
 * Added Int4 db_length to Options block.
 *
 * Revision 1.46  1997/02/18  21:03:00  madden
 * Added #define FILTER_NONE 0.
 *
 * Revision 1.45  1997/02/17  17:40:18  madden
 * Added seqalign to ResultHitlistptr
 *
 * Revision 1.44  1997/02/11  19:30:54  madden
 * Added program_name to Options.
 *
 * Revision 1.43  1997/02/10  20:27:01  madden
 * Changed some CharPtr's into Uint1Ptr's.
 *
 * Revision 1.42  1997/02/10  20:14:23  madden
 * replaced doubles by Nlm_FloatHi's.
 *
 * Revision 1.41  1997/02/10  20:03:58  madden
 * Added specific to BlastAllWordsPtr.
 *
 * Revision 1.40  1997/02/10  15:36:40  madden
 * added posConverged to the BlastSearchBlk.
 *
 * Revision 1.39  1997/02/06  14:27:15  madden
 * Addition of BlastAllWord structure.
 *
 * Revision 1.38  1997/02/03  13:02:12  madden
 * Added length to BLASTSubjectInfo.
 *
 * Revision 1.37  1997/01/17  17:41:44  madden
 * Added flags for position based BLAST.
 *
 * Revision 1.36  1997/01/13  15:37:05  madden
 * Changed prototypes for star_callback and tick_callback.
 *
 * Revision 1.35  1997/01/11  18:22:10  madden
 * Changes to allow S2 to be set.
 *
 * Revision 1.34  1997/01/09  17:44:35  madden
 * Added "bit_score" to BLASTResultHsp.
 *
 * Revision 1.33  1996/12/27  20:44:10  madden
 * Chnages to require that part of the query be included.
 *
 * Revision 1.32  1996/12/23  14:04:44  madden
 * Added gap_trigger.
 *
 * Revision 1.31  1996/12/20  21:11:40  madden
 * Changes to allow multiple hits runs only.
 *
 * Revision 1.30  1996/12/18  14:33:13  madden
 * Added high_score element.
 *
 * Revision 1.29  1996/12/17  17:27:03  madden
 * Count number of attempted gappings.
 *
 * Revision 1.28  1996/12/17  13:47:57  madden
 * Added star_proc.
 *
 * Revision 1.27  1996/12/16  14:35:48  madden
 * Added gapped_calculation Boolean
 *
 * Revision 1.26  1996/12/13  22:00:23  madden
 * Corrected starting point for gapped extension with traceback.
 *
 * Revision 1.25  1996/12/13  18:13:56  madden
 * Added tick callback functions
 *
 * Revision 1.24  1996/12/13  15:09:31  madden
 * Changes to parameters used for gapped extensions.
 *
 * Revision 1.23  1996/12/09  23:24:05  madden
 * Added parameters to control which sequences get a gapped alignment.
 *
 * Revision 1.22  1996/12/08  15:19:59  madden
 * Added parameters for gapped alignments.
 *
 * Revision 1.21  1996/11/27  21:56:57  madden
 * Removed define for XNU.
 *
 * Revision 1.20  1996/11/18  18:07:57  madden
 * *** empty log message ***
 *
 * Revision 1.19  1996/11/18  17:28:13  madden
 * Added BLAST_SEARCH_ALLOC_TRANS_INFO define.
 *
 * Revision 1.18  1996/11/18  15:45:40  madden
 * Defines for filter type added (by S. Shavirin),.
 *
 * Revision 1.17  1996/11/15  17:54:54  madden
 * Added support for alternate genetic codes for blastx, tblast[nx].
 *
 * Revision 1.16  1996/11/13  22:35:18  madden
 * Added genetic_code and db_genetic_code elements to blastdef.h
 *
 * Revision 1.15  1996/11/12  16:21:53  madden
 * Added context_factor
 *
 * Revision 1.14  1996/11/06  22:10:01  madden
 * translation_buffer changed from CharPtr to Uint1Ptr.
 *
 * Revision 1.13  1996/11/04  16:59:43  madden
 * Added translation_table and translation_table_rc elements
 * to BlastSearchBlk.
 *
 * Revision 1.12  1996/10/03  20:49:29  madden
 * Added xsum member to HSP_Link structure.
 * ,.
 *
 * Revision 1.11  1996/10/01  21:24:02  madden
 * Added e2.
 *
 * Revision 1.10  1996/09/26  13:02:32  madden
 * Removed ifdef for BLAST_COLLECT_STATS with counters.
 *
 * Revision 1.9  1996/09/12  21:13:46  madden
 * *** empty log message ***
 *
 * Revision 1.8  1996/09/11  22:21:51  madden
 * *** empty log message ***
 *
 * Revision 1.7  1996/09/11  19:14:09  madden
 * Added BLAST_OptionsBlkPtr structure and use thereof.
 *
 * Revision 1.6  1996/08/14  18:16:13  madden
 * removed frame from Context.
 *
 * Revision 1.5  1996/08/14  17:19:02  madden
 * Added frame to BlastSeqBlkPtr.
 *
 * Revision 1.4  1996/08/13  15:26:29  madden
 * Changes for tblastn.
 *
 * Revision 1.3  1996/08/09  22:11:12  madden
 * Added original_sequence to BlastSequenceBlk.
 *
 * Revision 1.2  1996/08/07  14:24:42  madden
 * Removed include for blast18p.h and objblst2.h
 *
 * Revision 1.1  1996/08/05  20:32:18  madden
 * Initial revision
 *
 * Revision 1.51  1996/08/02  14:20:06  madden
 * Removed multiproc strucutre.
 *
 * Revision 1.50  1996/07/31  13:09:17  madden
 * Changes for threaded blast.
 *
 * Revision 1.49  1996/07/24  12:01:28  madden
 * Changes for blastx
 *
 * Revision 1.48  1996/07/18  22:00:49  madden
 * Addition of BLAST_ExtendWordParams structure.
 *
 * Revision 1.47  1996/07/18  13:36:34  madden
 * Addition of the BLASTContextStructPtr.
 *
 * Revision 1.46  1996/07/16  14:37:42  madden
 * Removed _blast_link_structure .
 *
 * Revision 1.45  1996/07/11  16:03:58  madden
 * SaveCurrentHitlist keeps track of which set an HSP belongs to.
 *
 * Revision 1.44  1996/07/02  14:33:16  madden
 * Added hspcnt_max.
 *
 * Revision 1.43  1996/07/02  12:04:15  madden
 * HSP's saved on array, rather than linked list.
 *
 * Revision 1.42  1996/06/26  19:38:12  madden
 * Removed ifdef.
 *
 * Revision 1.41  1996/06/24  20:26:46  madden
 * Added dropoff_1st_pass and dropoff_2nd_pass to ParameterBlkPtr.
 *
 * Revision 1.40  1996/06/24  17:58:21  madden
 * Removed X_set parameter, added right and left dropoff's.
 *
 * Revision 1.39  1996/06/20  16:15:57  madden
 * Replaced int's with Int4's.
 *
 * Revision 1.38  1996/06/19  14:19:53  madden
 * Added define for BLASTSubjectInfoPtr.
 *
 * Revision 1.37  1996/06/17  19:03:07  madden
 * Rmoved unused structure.
 *
 * Revision 1.36  1996/06/14  17:58:13  madden
 * Changes to avoid nulling out arrays for every sequence.
 *
 * Revision 1.35  1996/06/13  21:03:06  madden
 * Added actual_window element to ExtendWord structure.
 *
 * Revision 1.34  1996/06/11  17:58:31  madden
 * Changes to allow shorter arrays for multiple hits type blast.
 *
 * Revision 1.33  1996/06/10  16:52:16  madden
 * Use bit-shifting and masking instead of dividing and remainder.
 *
 * Revision 1.32  1996/06/10  13:44:07  madden
 * Changes to reduce the size of the "already visited" array.
 *
 * Revision 1.31  1996/06/06  17:55:16  madden
 * Added number_of_bits to ParameterBlkPtr.
 *
 * Revision 1.30  1996/06/06  13:23:17  madden
 * Added elements cutoff_big_gap and ignore_small_gaps to ParameterBlkPt.
 *
 * Revision 1.29  1996/05/29  12:44:04  madden
 * Added structure BlastTimeKeeper.
 *
 * Revision 1.28  1996/05/28  14:16:32  madden
 * Added Int4's to collect statistics info.
 *
 * Revision 1.27  1996/05/23  21:55:04  madden
 * Removed unused variable initlen
 *
 * Revision 1.26  1996/05/23  21:48:23  madden
 * Removed unused defines.
 *
 * Revision 1.25  1996/05/16  19:51:09  madden
 * Added documentation block.
 *
 * Revision 1.24  1996/05/16  13:29:38  madden
 * Added defines for contiguous or discontiguous calls.
 *
 * Revision 1.23  1996/05/01  15:00:00  madden
 * Added BlastResults sturcture defs.
 *
 * Revision 1.22  1996/04/24  16:17:26  madden
 * Added new structure, BLAST_Link.
 *
 * Revision 1.21  1996/04/24  12:52:48  madden
 * ID's for sequences simplified.
 *
 * Revision 1.20  1996/04/03  19:14:35  madden
 * Removed defunct HSP ptr's.
 *
 * Revision 1.19  1996/03/29  21:27:43  madden
 * "hitlist" now kept on SeqAlign rather than HitList.
 *
 * Revision 1.17  1996/03/27  19:51:53  madden
 * "current_hitlist" added to Search Structure.
 *
 * Revision 1.16  1996/03/26  19:36:59  madden
 * Added  ReadDBFILEPtr to Search structure.
 *
 * Revision 1.15  1996/03/25  16:35:18  madden
 * Added old_stats.
 *
 * Revision 1.14  1996/02/28  21:37:43  madden
 * Added "trim" variables to segments for HSP.
 *
 * Revision 1.13  1996/02/06  22:51:13  madden
 * Added "prelim" to BlastSearch
 *
 * Revision 1.12  1996/02/02  19:25:32  madden
 * Added wfp_first and wfp_second to BlastParameterBlk for first and second pass.
 *
 * Revision 1.11  1996/01/29  21:12:07  madden
 * *** empty log message ***
 *
 * Revision 1.10  1996/01/23  16:31:47  madden
 * e_cutoff changed from BLAST_Score to double in ParameterBlk.
 *
 * Revision 1.9  1996/01/17  17:00:40  madden
 * Added gap parameters to ParameterBlk, dblen to SearchBlk.
 *
 * Revision 1.8  1996/01/17  13:45:58  madden
 * Added gap_prob and gap_decay_rate to ParameterBlk.
 *
 * Revision 1.7  1996/01/11  15:17:36  madden
 * Added process_num to ParameterBlk.
 *
 * Revision 1.6  1996/01/08  23:23:55  madden
 * removed "len" from HSP.
 *
 * Revision 1.5  1996/01/06  18:57:47  madden
 * Added BLAST_HSP_LINK structure.
 *
 * Revision 1.4  1995/12/28  21:26:05  madden
 * *** empty log message ***
 *
 * Revision 1.3  1995/12/26  23:04:14  madden
 * Added parameters to BlastParameterBlk.
 *
 * Revision 1.2  1995/12/21  23:10:41  madden
 * BLAST_Score prototypes moved to blastkar.h.
 *
 * Revision 1.1  1995/12/19  22:33:06  madden
 * Initial revision
 *
 * Revision 1.1  1995/12/08  15:48:23  madden
 * Initial revision
 *
 * */
#ifndef __BLASTSTR__
#define __BLASTSTR__
#ifdef __cplusplus
extern "C" {
#endif

#include <ncbi.h>
#include <lookup.h>
#include <blastkar.h>
#include <objalign.h>
#include <sequtil.h>
#include <readdb.h>
#include <gapxdrop.h>

/* the version of BLAST. */
#define BLAST_ENGINE_VERSION "2.0.8"
#define BLAST_RELEASE_DATE "Jan-05-1999"

/* Defines for program numbers. (Translated in BlastGetProgramNumber). */
#define blast_type_undefined 0
#define blast_type_blastn 1
#define blast_type_blastp 2
#define blast_type_blastx 3
#define blast_type_tblastn 4
#define blast_type_tblastx 5

/* defines for strand_option, determines which strand of query to compare. */
#define BLAST_TOP_STRAND 1
#define BLAST_BOTTOM_STRAND 2
#define BLAST_BOTH_STRAND 3

/* Defines that specify whether or not BLAST should delete some memory, or
	leave it up to the caller.
*/
#define BLAST_OWN 0
#define BLAST_NOT_OWN 1

/********************************************************************
*
*	define for collecting BLAST stats.
*
***********************************************************************/

#define BLAST_COLLECT_STATS

/********************************************************************
*
*	Structure to save timing info. in.  Right now this only
*	works for UNIX.
*
********************************************************************/

typedef struct _blast_time_keeper {
                FloatLo	user, 	/* CPU time in user space of the process. */
			system, /* CPU time used by system. */
			total;	/* total CPU time (i.e., both of the above). */
        } BlastTimeKeeper, PNTR BlastTimeKeeperPtr;


/*****************************************************************
*
*	Used for pruing SeqALigns that are too big.
*
********************************************************************/

typedef struct _blast_prune_hits_from_sap {
		SeqAlignPtr sap;
		Int4 	original_number,	/* how may unique hits were there originally. */
			number;		/* How many hits on SeqALignPtr above. */
		Boolean allocated; /* If FALSE, SeqAlignPtr above does NOT belong to this struc.*/
        } BlastPruneSapStruct, PNTR BlastPruneSapStructPtr;

/***************************************************************************
  Macros added by Andy Neuwald in order to allow easy modification of matrices.
***************************************************************************/

#define  MtrxScorePosSearch(S,x,y)	((S)->posMatrix[(x)][(y)])
#define  PtrMtrxScorePosSearch(S,x)	((S)->posMatrix[(x)])

/*****
#define  MtrxScorePosSearchi2(S,x,y)	\
	((S)->posMatrix[( (x) %(S)->query_length)][(y)])
#define  PtrMtrxScorePosSearch2(S,x)	\
	((S)->posMatrix[( (x) %(S)->query_length)])
*****/

/********************************************************************

	Defines for discontiguous word hits on 1st and 2nd pass.

********************************************************************/

#define BLAST_NO_PASS_DISCONTIG 0
#define BLAST_1ST_PASS_DISCONTIG 1
#define BLAST_2ND_PASS_DISCONTIG 2
#define BLAST_BOTH_PASS_DISCONTIG 3

#define CODON_LENGTH 3  /* three is always the codon length. */

#define BLAST_SMALL_GAPS 0
#define BLAST_LARGE_GAPS 1

/*********************************************************************
    Filter types definitions
*********************************************************************/

#define FILTER_NONE 0
#define FILTER_DUST 1
#define FILTER_SEG  2

/**********************************************************************
	Structure for the blast options (available to user/programmer).
	This should be filled in by the "Main" program before blast
	is called.
***********************************************************************/

typedef struct _blast_optionsblk {
	Nlm_FloatHi gap_decay_rate,	/* decay rate. */
		    gap_prob;	/* Prob of decay. */
	Int4	gap_size,	/* Small gap size. */
		window_size, 	/* Multiple Hits window size (zero for single hit algorithm) */
		threshold_first, /* Threshold for extending hits (preliminary pass), zero if one-pass algorithm is used. */ 
		threshold_second;/* Threshold for extending hits (second pass) */
	Nlm_FloatHi	expect_value, 	/* Expectation value (E) */
			e2; 	  	/* Expect value for a single HSP */
	/* These two scores are zero, unless they've been set, then they set
	the expect_value and e2 above. */
	Int4		cutoff_s,	/* score corresponds to expect_value above.*/
			cutoff_s2;	/* score corresponds to e2 above. */
	Boolean two_pass_method; /* should two passes be used? */
	Boolean	multiple_hits_only; /* Only the multiple hits alg. used. */
	Int4	hitlist_size;	/* How many hits should be returned. */
	Nlm_FloatHi number_of_bits; /* Number of bits to initiate 2nd pass (default is used if zero) */
	Int4	dropoff_1st_pass, /* dropoff ("X") used for 1st pass. */
                dropoff_2nd_pass; /* dropoff ("X") used for 2nd pass. */
	Int2	number_of_cpus;	/* How many CPU's. */
	CharPtr matrix;		/* name of matrix to use. */
	Boolean old_stats; /* Use old stats (option may disappear later) */
	Boolean do_sum_stats;   /* Should sum statistics be used? */
	Boolean use_large_gaps; /* Use only large gaps for linking HSP's with sum stats. */
	Int2	wordsize;	/* size of word used to find hits. */
	Int2	penalty, reward; /* penalty and reward, only for blastn */
	/* The ID numbers from gc.prt are used for the genetic codes. */
	Int4	genetic_code,		/* genetic code for query (blastx, tblastx) */
		db_genetic_code;	/* genetic code for db (tblast[nx]). */
        Int4 filter;          /* filter type 0 mean no filter
                                 non-zero value indicate filer type */
	CharPtr filter_string;	/* String specifying the type of filtering and filter options. */
	Boolean		gapped_calculation; /* Is a gapped calc. being done? */
	/* The next three are used ONLY for gapped alignments. */
	Int4		gap_open,	/* Cost to open a gap (NO extension). */
			gap_extend,	/* Cost to extend a gap one letter. */
			gap_x_dropoff,	/* X-dropoff (in bits) used by Gapped align routine. */
			gap_x_dropoff_final;	/* X-dropoff (in bits) used by Gapped align routine for FINAL alignment. */
	Nlm_FloatHi	gap_trigger; /* Score (in bits) to gap, if an HSP gaps well. */

	Boolean		discontinuous;	/* Should the SeqAlign be discontinuous.*/
	/* What region of the query is required for the alignment.  If start is
	zero and end is -1 (the entire query), then these are not checked. */
	Int4		required_start,
			required_end;
	Int4		db_length;	/* database size used for stat. calcul. */
	Int4		dbseq_num;	/* number of database sequences used for stat. calcul. */
	Nlm_FloatHi	searchsp_eff;	/* Effective search space to be used. */
			
	/* Options for postion based blast. */
	Nlm_FloatHi	ethresh;
	Int4		maxNumPasses,
			pseudoCountConst;
	CharPtr program_name;		/* program name, for reference. */
	Int4 cpu_limit;	/* timeout total. */
	/* Used for region-dependent limits when storing hits. */
        Int4    hsp_range_max,          /* maximum hits for a range */
                block_width;            /* width of a block */
	Boolean perform_culling;	/* Should results be culled at all? */
        Boolean isPatternSearch;        /* Is this a use of PHI-BLAST?*/
	CharPtr		gifile;		/* name of file containing list of gis on server */
	ValNodePtr	gilist;		/* list of gis specified by client */
	Boolean		do_not_reevaluate;	/* Don't perform BlastReevaluateWithAmbiguities. */
	/* These options allow a subset of the database to be examined.  IF they
		are set to zero, then the entire database is examined. */
	Int4		first_db_seq,		/* 1st sequence in db to be compared. */
			final_db_seq;		/* Final sequence to be compared. */
	CharPtr		entrez_query;	/* user specified Entrez query to make selection from databases */
	CharPtr		org_name;	/* user specified name of organizm;  corresponding .gil file will be used */
	Uint1		strand_option;	/* BLAST_TOP_STRAND, BLAST_BOTTOM_STRAND, or BLAST_BOTH_STRAND.  used by blast[nx] and tblastx */
	} BLAST_OptionsBlk, PNTR BLAST_OptionsBlkPtr;

/****************************************************************************

	PARAMETER BLOCK: parameters for the BLAST search entered by on
	command line by user.

*****************************************************************************/

typedef struct _blast_parameterblk {
        BLAST_Score     threshold,      /* threshold for extending a word hit*/
        		threshold_first, /* threshold for 1st pass. */
        		threshold_second, /* threshold for 2nd pass. */
                        X,              /* drop-off score for extension. */
			dropoff_1st_pass, /* dropoff ("X") used for 1st pass. */
			dropoff_2nd_pass, /* dropoff ("X") used for 2nd pass. */
                        cutoff_s,	/* Score to report a hit. */
                        cutoff_s2,	/* Score to report a hsp. */
			cutoff_s_first, /* Score (S2) to use on 1st pass */
			cutoff_s_second, /* Score (S2) to use on 2nd pass and
			   for "small" gaps in link_hsps (in blast.c) */
	/* Max value of s2, used if s2 is set or s2 becomes larger than s. */
			cutoff_s2_max,	
			cutoff_big_gap; /* cutoff value for a "big" gap in
			   link_hsps (in blast.c). */
	Nlm_FloatHi	cutoff_e,	/* Expect value to report a hit. */
                        cutoff_e2,	/* Expect value to report a hsp. */
			number_of_bits; /* number of bits of significance, used
			   to calculate cutoff_s_first (above). */
	Boolean		threshold_set, /*TRUE if threshold set on command-line*/
			cutoff_s_set,	/* TRUE if cutoff score set on c-l */
			cutoff_s2_set,	/* TRUE if cutoff score2 set on c-l */
			cutoff_e_set,	/* TRUE if cutoff expect set on c-l */
			cutoff_e2_set,	/* TRUE if cutoff expect2 set on c-l */
			ignore_small_gaps, /* ignore small gaps if TRUE, set by
			   CalculateSecondCutoffScore in blast.c if the search 
			   space is smalled than 8*gap_size*gap_size. */
			window_size_set;/* TRUE if window size set for MHBLAST*/
        Boolean         sump_option;    /* TRUE if sump is used. */
	Int4		gap_size,	/* max. gap allowed for small gaps.*/
			window_size;	/* used for multiple hits BLAST. */
	Nlm_FloatHi	gap_prob; 	/* prob. of gap of size "gap" (above).*/
	Nlm_FloatHi	gap_decay_rate; /* prob. of only one HSP */
	Int2		process_num;	/* max # processrs permitted (for MP).*/
	Boolean		old_stats;	/* Use "old" stats if TRUE. */
	Boolean 	do_sum_stats;   /* Should sum statistics be used? */
	Boolean         use_large_gaps; /* Use only large gaps for linking HSP's with sum stats. */
	Boolean		two_pass_method; /* should two passes be used? */
	Boolean		multiple_hits_only; /* Only the multiple hits alg. used. */
	Boolean		discontinuous;	/* Should discontinuous SeqAlign's be produced? */
	Boolean		gapped_calculation; /* Is a gapped calc. being done? */
	Boolean		do_not_reevaluate;	/* Don't perform BlastReevaluateWithAmbiguities. */
	/* The next three are used ONLY for gapped alignments. */
	Int4		gap_open,	/* Cost to open a gap (NO extension). */
			gap_extend,	/* Cost to extend a gap one letter. */
			gap_x_dropoff,	/* X-dropoff used by Gapped align routine. */
			gap_x_dropoff_final;	/* X-dropoff (in bits) used by Gapped align routine for FINAL alignment. */
	Nlm_FloatHi	gap_trigger; /* Score (in bits) to gap, if an HSP gaps well.*/

	/* Options for postion based blast. */
	Nlm_FloatHi	ethresh;
	Int4		maxNumPasses,
			pseudoCountConst;
	Int4 cpu_limit;	/* timeout total. */
        Int4    hsp_range_max,          /* maximum hits for a range */
                block_width,            /* width of a block */
		max_pieces;		/* Max number of pieces allowed (query_length/block_width) */
	Boolean perform_culling;	/* determines whether culling should be used or not.
					If not, then hsp_range_max, block_width, and max_pieces are ignored. */
	/* These options allow a subset of the database to be examined.  IF they
		are set to zero, then the entire database is examined. */
	Int4		first_db_seq,		/* 1st sequence in db to be compared. */
			final_db_seq;		/* Final sequence to be compared. */
        } BLAST_ParameterBlk, PNTR BLAST_ParameterBlkPtr;

typedef Nlm_Int4	BLAST_Diag, PNTR BLAST_DiagPtr;

/*
	BLAST_ExtendWord contains information about which diagonals
	have been extended over (i.e., which diagonals have been 
	tested).  This structure will be duplicated once for each
	context as every context is different.
*/
typedef struct _blast_extend_word {
		Int4Ptr	_buffer, /* The "real" buffer for diag_level, version,
				and last_hit arrays. */
			diag_level, /* How far the "diagonal" has been traversed. */
			version, /* "version" of the diagonal. */
			last_hit;/* keeps track of last hit for Multiple hits. */
		Int4	actual_window; /* The actual window used if the multiple
				hits method was used and a hit was found. */	
	} BLAST_ExtendWord, PNTR BLAST_ExtendWordPtr;

/*
	BLAST_ExtendWordParams contains parameters about the extensions.
	Only one copy of this structure is needed, regardless of how many
	contexts there are.
*/
typedef struct _blast_extend_word_params {
		Int4	bits_to_shift; /* how many bits should the diagonal be
				shifted to get the "version" */
		Int4	min_diag_length, /* Min. length of diagonal, actuall
				2**bits_to_shift. */
			min_diag_mask; /* Used to mask off everything above
				min_diag_length (mask = min_diag_length-1). */
		Int4	offset; /* "offset" added to query and subject position
				so that "diag_level" and "last_hit" don't have
				to be zeroed out every time. */
		Int4	window;	/* The "window" size, within which two (or more)
				hits must be found in order to be extended. */
		/* Used by BLAST_ExtendWordNew to decide whether or not
		to prepare the structure for multiple-hit type searches.
		If TRUE, multiple hits are not neccessary, but possible. */
		Boolean multiple_hits;  
	} BLAST_ExtendWordParams, PNTR BLAST_ExtendWordParamsPtr;
/*
	Data block to describe a single sequence.
*/

typedef struct blast_sequence_block {
	Uint1Ptr	sequence,	/* Actual (perhaps transl.) sequence. */
		sequence_start; /* Start of sequence, used if the sequence is preceded by a NULLB.  Sequences
				starting with a NULLB are used by BlastWordExtend_L1. */
	Int4	length,		/* length of sequence. */
		original_length,/* length before translation. */
		effective_length;/* effective length, used only by query. */
	Int2 frame;		/* frame of the sequence. */
} BlastSequenceBlk, PNTR BlastSequenceBlkPtr;


typedef struct _blast_seg {
		Int2		frame;
		Int4		offset;	/* start of hsp */
		Int4		length;	/* length of hsp */
		Int4		end;	/* end of HSP */
		Int4		offset_trim;	/* start of trimmed hsp */
		Int4		end_trim;	/* end of trimmed HSP */
		/* Where the gapped extension (with X-dropoff) started. */
		Int4		gapped_start;
	} BLAST_Seg, PNTR BLAST_SegPtr;

#define BLAST_NUMBER_OF_ORDERING_METHODS 2

/*
	The following structure is used in "link_hsps" to decide between
	two different "gapping" models.  Here link is used to hook up
	a chain of HSP's (this is a VoidPtr as _blast_hsp is not yet
	defined), num is the number of links, and sum is the sum score.
	Once the best gapping model has been found, this information is
	transferred up to the BLAST_HSP.  This structure should not be
	used outside of the function link_hsps.
*/
typedef struct _blast_hsp_link {
		/* Used to order the HSP's (i.e., hook-up w/o overlapping). */ 
	VoidPtr	link[BLAST_NUMBER_OF_ORDERING_METHODS]; 
		/* number of HSP in the ordering. */
	Int2	num[BLAST_NUMBER_OF_ORDERING_METHODS];
		/* Sum-Score of HSP. */
	Int4	sum[BLAST_NUMBER_OF_ORDERING_METHODS]; 
		/* Sum-Score of HSP, multiplied by the appropriate Lambda. */
	Nlm_FloatHi	xsum[BLAST_NUMBER_OF_ORDERING_METHODS]; 
	} BLAST_HSP_LINK, PNTR BLAST_HSP_LINKPtr;
/*
	BLAST_NUMBER_OF_ORDERING_METHODS tells how many methods are used
	to "order" the HSP's.
*/

typedef struct _blast_hsp {
		struct _blast_hsp PNTR next, /* the next HSP */
				  PNTR prev; /* the previous one. */
		BLAST_HSP_LINK	hsp_link;
/* Is this HSp part of a linked set? */
		Boolean		linked_set;
/* which method (max or no max for gaps) was used? */
		Int2		ordering_method; 
/* how many HSP's make up this (sum) segment */
		Int4		num;
/* sumscore of a set of "linked" HSP's. */
		BLAST_Score	sumscore;
		/* If TRUE this HSP starts a chain along the "link" pointer. */
		Boolean 	start_of_chain;
		BLAST_Score	score;
		Nlm_FloatHi	pvalue;
		Nlm_FloatHi	evalue;
		BLAST_Seg query,	/* query sequence info. */
			subject;	/* subject sequence info. */
		Int2		context;	/* Context number of query */
                GapXEditBlockPtr gap_info; /* ALL gapped alignment is here */
		Int4 num_ref;
	} BLAST_HSP, PNTR BLAST_HSPPtr;

typedef struct _blast_hitlist {
	struct _blast_hitlist	PNTR next;
	BLAST_HSPPtr PNTR	hsp_array; /* head of linked list of HSPs */
	Int4		hspmax, /* max no. of HSPs allowed per hit list */
			hspcnt, /* no. of HSPs in hit list */
			hspcnt_max; /* no. of HSPs in hitlist, before reaping */
	Boolean		further_process; /* This sequence has been found interesting,
					    it should be further processed by a gapped
					    alignment etc. */
	} BLAST_HitList, PNTR BLAST_HitListPtr;

/*
	The next two structures are the final output produced by BLAST.  Formatters should then
	convert the data into SeqAligns or the BLAST ASN.1 spec.  
*/

typedef struct _blast_results_hsp {
		Int2		ordering_method;/* determines whether large or small gap was used. */
		Int4 		number;	/* number of HSP's used to calculate the p-value. */
		BLAST_Score	score;	/* score of this HSP. */
		Nlm_FloatHi	e_value,/* expect value of this set of HSP's. */
				p_value,/* p-value of this set of HSP's. */
				bit_score; /* above score * lambda/ln2 */
		Int2		context;	/* context number of query. */
		Int2		query_frame, /* frame of query, non-zero if transl. */
				subject_frame; /* frame of subject, non-zero if transl. */
		Int4 		query_offset,	/* Start of the query HSP. */
				query_length,	/* Length of the query HSP. */
				subject_offset,	/* Start of the subject HSP. */
				subject_length, /* Length of the subject HSP.*/
				hspset_cnt;	/* which set of HSP's? */
	/* Starting points (on original HSP) for a gapped extension with X dropoff. */
		Int4		query_gapped_start,
				subject_gapped_start;

                GapXEditBlockPtr gap_info; /* ALL gapped alignment is here */
		struct _blast_result_hitlist PNTR point_back;
		struct _blast_heap_struct PNTR back_left, PNTR back_right;
		} BLASTResultHsp, PNTR BLASTResultHspPtr;

/*
	The following structure contains the subject info, if the readdb
	facility is not being used.  Then the subject information is
	kept here.  Otherwise this structure is NULL.
*/
typedef struct _blast_subject_info {
		SeqIdPtr sip;	/* ID of the subject. */
		CharPtr defline; /* Defline of the subject. */
		Int4 length; 	/* untranslated length of the database sequence. */
		} BLASTSubjectInfo, PNTR BLASTSubjectInfoPtr;

typedef struct _blast_result_hitlist {
		BLASTResultHspPtr hsp_array;	/* An array holding the HSP's. */
		Nlm_FloatHi	best_evalue;	/* best evalue in all the HSP's. */
		Int4	high_score; 	/* HSP with highest score. */
		Int4	hspcnt,		/* Number of HSP's. */
			subject_id;	/* ID of the subject. */
		Int2    db_id;          /* ID (0,1,2...) of the db if multiple db's searched. */
		Int4    subject_length; /* length of the database sequence. */
		BLASTSubjectInfoPtr subject_info; /* Subject info if the readdb facility is not being used. */
		SeqAlignPtr seqalign; /* alignment, if this a gapped calculation. */
		Int4 num_ref;
		} BLASTResultHitlist, PNTR BLASTResultHitlistPtr;


typedef struct _blast_heap_struct {
  Int4 cutvalue;	/* start of a region? */
  BLASTResultHspPtr PNTR heap;
  Int4 num_in_heap;	/* Number in 'heap' */
  Int4 num_of_ref;
  struct _blast_heap_struct PNTR next, PNTR prev;
} BLASTHeapStruct, PNTR BLASTHeapPtr;

/*
	Holds the results already saved.
*/

typedef struct _blast_results_struct {

		BLASTResultHitlistPtr PNTR results;
		Int4	hitlist_count,	/* Number of hitlists saved on results array already. */
			hitlist_max, 	/* Length of results array. */
			max_pieces;	/* For range-dependent limits. */
		BLASTResultHspPtr **heap;
        	Int4 *num_in_heap;
		BLASTHeapPtr heap_ptr;
		} BLASTResultsStruct, PNTR BLASTResultsStructPtr;

/*
	Holds the data for all possible words that might be used by BLAST.
*/

typedef struct _blast_all_words {
		Uint1Ptr *array;	/* All the possible words */
		Int4 	num_of_cols, 
			wordsize;
		Boolean rows_allocated,	/* are the rows (of length the wordsize) alloc.*/
			specific;	/* specific (limited) words are to be indexed. */
	} BlastAllWord, *BlastAllWordPtr;
		
typedef struct _blast_seqid_list {
		SeqIdPtr 	seqid_list;	/* A list of SeqId's (may or may not be in the database) that should
						serve as subject sequences for BLAST'ing. */
		Uint1 		mode;	/* Either BLAST_OWN, or BLAST_NOT_OWN.  Specifies whether they should be deleted. */
	} BlastSeqIdList, *BlastSeqIdListPtr;
		

/*
	Contains gi and ordinal number for use by random access BLAST.
*/
typedef struct _double_int4 {
        Int4    gi,
                ordinal_id;
} BlastDoubleInt4, *BlastDoubleInt4Ptr;


typedef struct _blast_gi_list {
		BlastDoubleInt4Ptr	gi_list;	/* List of gi's. */
		BlastDoubleInt4Ptr	*gi_list_pointer;	/* Pointer to above list. */
		Int4		current,	/* Current position in gi list. */
				total;		/* total number of gi's. */
	} BlastGiList, *BlastGiListPtr;


/*
	used for keeping start and stop of hits to query, for ALU filtering.
*/
typedef struct _blast_hit_range {
	BlastDoubleInt4Ptr      range_list;        /* ranges. */
        BlastDoubleInt4Ptr      *range_list_pointer;       /* Pointer to above list. */
	Int4		current,	/* current position in list. */
			total;		/* total number in list. */
	} BlastHitRange, *BlastHitRangePtr;

/*
	Contains BLAST error messages.
*/

typedef struct _blast_error_msg {
	Uint1 level;/* corresponds to levels of ErrPostEx [none(0), info(1), warn(2), error(3) and fatal(4)] */
	CharPtr msg;
} BlastErrorMsg, *BlastErrorMsgPtr;

/*
	Holds data for each "context" (which is generally equal to
	one frame of the query).  blastx would have six contexts,
	blastp would have one.
*/

typedef struct _blast_context_structure {
	Boolean query_allocated;/* The BlastSequenceBlkPtr IS allocated. */
	BlastSequenceBlkPtr query;  /* query sequence. */
	BLAST_ExtendWordPtr ewp;/* keep track of diagonal etc. for each frame */
	ValNodePtr location;    /* Where to start/stop masking. */
		} BLASTContextStruct, PNTR BLASTContextStructPtr;


/*
	Structure used for matrix rescaling. 
*/

typedef struct _blast_matrix_rescale {
	Int4 		alphabet_size,
			query_length;	/* length of query. */
	Uint1Ptr	query;
	Nlm_FloatHi 	*standardProb;
	Int4Ptr  	*matrix;
	Int4Ptr  	*private_matrix;
	BLAST_KarlinBlkPtr 	*kbp_std, 
				*kbp_psi, 
				*kbp_gap_std, 
				*kbp_gap_psi;
	Nlm_FloatHi	lambda_ideal,
                	K_ideal;
} BlastMatrixRescale, *BlastMatrixRescalePtr;
	
		
		
/*
	The central structure for the BLAST search.  This structure
	should contain data (or pointers to data) for all the
	information in a BLAST search.
*/


#define BLAST_SEARCH_ALLOC_QUERY 1
#define BLAST_SEARCH_ALLOC_SUBJECT 2
#define BLAST_SEARCH_ALLOC_PBP 4
#define BLAST_SEARCH_ALLOC_SBP 8
#define BLAST_SEARCH_ALLOC_WFP_FIRST 16
#define BLAST_SEARCH_ALLOC_WFP_SECOND 32
#define BLAST_SEARCH_ALLOC_EWPPARAMS 64
#define BLAST_SEARCH_ALLOC_CONTEXT 128
#define BLAST_SEARCH_ALLOC_RESULTS 256
#define BLAST_SEARCH_ALLOC_READDB 512
#define BLAST_SEARCH_ALLOC_TRANS_INFO 1024
#define BLAST_SEARCH_ALLOC_ALL_WORDS 2048
#define BLAST_SEARCH_ALLOC_QUERY_SLP 4096

typedef struct blast_search_block {
	Int4		allocated; 
/* bit fields specify which structures from below are allocated.  If 
a field is allocated, then it's bit is non-zero. 

		structure     		bit-field (define)
		-----------------------------------------
		query			BLAST_SEARCH_ALLOC_QUERY
		subject			BLAST_SEARCH_ALLOC_SUBJECT
		pbp			BLAST_SEARCH_ALLOC_PBP
		sbp			BLAST_SEARCH_ALLOC_SBP
		wfp_first       	BLAST_SEARCH_ALLOC_WFP_FIRST
		wfp_second      	BLAST_SEARCH_ALLOC_WFP_SECOND
		ewp_params		BLAST_SEARCH_ALLOC_EWPPARAMS
		context			BLAST_SEARCH_ALLOC_CONTEXT
		result_struct		BLAST_SEARCH_ALLOC_RESULTS
		rdfp	        	BLAST_SEARCH_ALLOC_READDB
		translation_table       BLAST_SEARCH_ALLOC_TRANS_INFO
		translation_table_rc
		all_words		BLAST_SEARCH_ALLOC_ALL_WORDS
		query_slp		BLAST_SEARCH_ALLOC_QUERY_SLP
*/

/*
	Specifies whether the search is position based or not.
*/
	Boolean positionBased;
	Boolean posConverged;
/*
	Specifies that the query sequence was invalid (e.g., XXXXXXXXXXXXXXXXXXXXXX).
*/
	Boolean query_invalid;
/*
	The BLASTContextStructPtr is an array and each element contains
	information about the query sequence and the frame number.
	If there are six frames (e.g., blastx) then the BLASTContextStructPtr
	is six elements long; if there's one frame (e.g., blastp) then
	BLASTContextStructPtr is one element long.

	number_of_contexts states how long the context array is.
*/	
	BLASTContextStructPtr context;
	Int2	first_context,
		last_context;
/* 
	The GapAlignBlkPtr used by ALIGN (in gapxdrop.c) for gapped alignments.
*/
	
	GapAlignBlkPtr gap_align;

/*
	All the possible words.
*/
	BlastAllWordPtr all_words;
/*
        Set the context_factor, which specifies how many different 
        ways the query or db is examined (e.g., blastn looks at both
        stands of query, context_factor is 2).
*/
	Int2 context_factor;

/*
	What type of search (e.g., blastp, blastx, etc.)?
*/
	CharPtr prog_name;
	Uint1 prog_number;
/*
	translation_table and translation_table_rc holds the translation
	from ncbi2na to ncbistdaa for normal and reverse-complement
	translations.  Only used and initialized with tblast[nx].
	Initialized by GetPrivatTranslationTable
*/
	Uint1Ptr translation_table,
		 translation_table_rc;

/*
	ValNodePtr containing error messages. 
*/
	ValNodePtr error_return;

/*
	ValNodePtr containing masking SeqLocPtr's
*/
	ValNodePtr mask;
/*
	What genetic codes are we using to translate the query or database
	when needed.  Based upon NCBI genetic codes.
*/
	CharPtr genetic_code,		/* genetic code used for query. */
		db_genetic_code;	/* genetic code used for database. */

/* 	
	The BlastSequenceBlk's subject hold info about the subject.  
	Info about the original sequence is in original_seq.  This will
	be NULL if the sequence was not translated. 
*/
	Uint1Ptr translation_buffer;	/* Buffer for (tblast[nx]) db translations*/
	Int4 translation_buffer_size;	/* size of translation_buffer. */
	CharPtr original_seq;	/* Original (i.e.,  untransl.) sequence. */
	BlastSequenceBlkPtr	subject;/* subject sequence. */

/*
	SeqLocPtr for the query, owned by the called and not by BLAST.
*/
	SeqLocPtr query_slp;

/* Id's for the query and subject. */
	SeqIdPtr		query_id;	/* ID for the query, any form. */
	Int4			subject_id;	/* the number of the subject, in the DB. */
	BLAST_ParameterBlkPtr pbp;	/* options selected. */
	BLAST_ScoreBlkPtr sbp;		/* info on scoring. */
	BLAST_ExtendWordParamsPtr ewp_params; /* parameters for extensions.*/

/* 	For the two-pass method two BLAST_WordFinderPtr's are required.
	The actual wfp's are in wfp_first and wfp_second.  "wfp" is just
	a pointer to one of those two.  If they have been allocated (at all)
	is signified by setting the bit-fields above. 
*/
	BLAST_WordFinderPtr     wfp, 	/* find initial words. */
				wfp_first, /* words for first pass. */
				wfp_second;/* words for second pass. */
/*	For the two-pass this should be set to TRUE on the first (preliminary)
	pass and FALSE on the second pass.
*/
	Boolean			prelim;
/*
	The "current" hit, that is the one being worked on right now.  
	If a hitlist is deemed significant, then "current_hitlist" is 
	moved to "seqalign".  current_hitlist_purge specifies 
	whether the hitlist should be purged after each call to a
	WordFinder; it will generally be purged except for non-initial
	frames of tblast[nx].
*/
	Boolean			current_hitlist_purge;
	BLAST_HitListPtr	current_hitlist;
/*
	The worst evalue seen by this thread so far.
	Only filled in if the hitlist is already full, otherwise
	it should be DBL_MAX.
*/
	Nlm_FloatHi	worst_evalue;
/*
	Size of the HSP array on the "current_hitlist"
*/
	Int4 hsp_array_size;
/*
	Contains hits that are significant. 
*/
	Int4			result_size;
	BLASTResultsStructPtr	result_struct;

	Int8			dblen;	/* total length of the database. */
	Int8                    dblen_eff;      /* effective length of the database. */
	Int8                    dblen_eff_real;      /* effective length of the database. */
	Int4                    dbseq_num;      /* number of sequences in the database. */
	Int4                    length_adjustment; /* amount removed from end of query and db sequences. */
	Nlm_FloatHi		searchsp_eff;	/* Effective search space (used for statistics). */
	ReadDBFILEPtr		rdfp, /* I/O PTR for database files. */
				rdfp_list;	/* linked rdfp list of all databases. */
/* The subject info (id and defline) is kept here for the current sequence
	if the readdb facility is not used.  This structure should only
	be used if rdfp is NULL.
*/
	BLASTSubjectInfoPtr subject_info;
/*
	A list of SeqId's from the database to use for the BLAST run.
*/
	BlastSeqIdListPtr blast_seqid_list;
	BlastGiListPtr	blast_gi_list;
/*
	start and stop of query that must be included for an alignment
	to be counted.  The Boolean whole_query specifies whether these
	are valid (i.e., have been set) or not.
*/
	Boolean whole_query;
	Int4 required_start, required_end;

/*
	Callback functions to indicate progress, or lack thereof.
*/
	int (LIBCALLBACK *tick_callback)PROTO((Int4 done, Int4 positives));
	int (LIBCALLBACK *star_callback)PROTO((Int4 done, Int4 positives));
/*
        Callback function to handle results (e.g., print them out for neighboring)
        in place of BlastSaveCurrentHitlist.
*/
        int (LIBCALLBACK *handle_results)PROTO((VoidPtr search));
/*
	These "counters" keep track of how often certain operations
	were performed.

	This counting is performed only if BLAST_COLLECT_STATS is defined.
*/
	Int4	first_pass_hits,	/* no. of hits on 1st pass. */
		second_pass_hits,	/* no. of hits on 2nd pass. */
		second_pass_trys,	/* no. of seqs that made it to 2nd pass. */
		first_pass_extends,	/* no. extended on 1st pass. */
		second_pass_extends,	/* no. extended on 2nd pass. */
		first_pass_good_extends,/* no. successfully extended on 1st pass. */
		second_pass_good_extends,/* no. successfully extended on 2nd pass. */
		number_of_seqs_better_E,/* how many sequences were better than E. */
		prelim_gap_no_contest,	/* No. of HSP's under E=10 alone. */
		prelim_gap_passed,	/* No. of HSP's that passed prelim gapping. */
		prelim_gap_attempts,	/* No. of HSP's we attempted to gap. */
		real_gap_number_of_hsps; /* How many HSP's were gapped in BlastGetGappedScore. */
} BlastSearchBlk, PNTR BlastSearchBlkPtr;



#ifdef __cplusplus
}
#endif
#endif /* !__BLASTSTR__ */