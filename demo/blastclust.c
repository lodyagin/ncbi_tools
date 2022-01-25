/*  $RCSfile: blastclust.c,v $  $Revision: 6.6 $  $Date: 2000/05/23 18:53:55 $
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
* Author: Ilya Dondoshansky 
*
* File Description:
*   This is server which does neighboring searches
*
* ---------------------------------------------------------------------------
* $Log: blastclust.c,v $
* Revision 6.6  2000/05/23 18:53:55  dondosha
* Function GetAccessionFromSeqId moved to sequtil.c; fixed minor bugs
*
* Revision 6.5  2000/05/22 18:19:59  dondosha
* In case search set up fails, destruct all necessary stuff before going to next query
*
* Revision 6.4  2000/05/17 17:44:50  dondosha
* Cleaned from unused variables
*
* Revision 6.3  2000/05/05 18:16:56  dondosha
* Enhanced to process all types of SeqIds
*
* Revision 6.2  2000/05/03 16:50:40  dondosha
* Removed system calls, added sorting of clusters by size, sequences within clusters by length
*
* Revision 6.1  2000/04/27 14:47:18  dondosha
* Clustering of protein neighbours - initial revision
*
* Revision 1.1  2000/04/12 12:53:04   dondosha
* Clustering of protein neighbors
*
* ===========================================================================
*/
#include <ncbi.h>
#include <objseq.h>
#include <objsset.h>
#include <sequtil.h>
#include <seqport.h>
#include <tofasta.h>
#include <blast.h>
#include <blastpri.h>
#include <simutil.h>
#include <txalign.h>
#include <gapxdrop.h>
#include <sqnutils.h>
#include <mbalign.h>
#include <mblast.h>

#define DEBUG 0

/* Used by the callback function. */
FILE *global_fp=NULL;
static Int4Ptr root;
static TNlmMutex root_mutex;

static int LIBCALLBACK
tick_callback(Int4 sequence_number, Int4 number_of_positive_hits)
{
   return 0;
}

typedef struct cluster_parameters
{
   Boolean bidirectional;
   FloatHi length_threshold;
   FloatHi score_threshold;
   FILE *logfp;
} ClusterParameters, PNTR ClusterParametersPtr;

typedef struct cluster_log_header
{
   Boolean numeric_id_type;
   Int4 size;
} ClusterLogHeader, PNTR ClusterLogHeaderPtr;

static ClusterParametersPtr global_parameters;

typedef struct cluster_log_info
{
   Int4 id1, id2;
   Int4 length1, length2;
   Int4 hsp_length1, hsp_length2;
   FloatHi bit_score;
} ClusterLogInfo, PNTR ClusterLogInfoPtr;

typedef struct blast_cluster_element
{
   Int4 index;
   Int4 len;
} BlastClusterElement, PNTR BlastClusterElementPtr;

typedef struct blast_cluster
{
   Int4 size;
   BlastClusterElementPtr PNTR elements;
} BlastCluster, PNTR BlastClusterPtr;

static int LIBCALLBACK
compare_cluster_elements(VoidPtr v1, VoidPtr v2)
{
   BlastClusterElementPtr e1, e2;
   
   e1 = *(BlastClusterElementPtr PNTR) v1;
   e2 = *(BlastClusterElementPtr PNTR) v2;

   if (e1->len > e2->len) return -1;
   else if (e1->len < e2->len) return 1; 
   else return 0;
}

static int LIBCALLBACK
compare_clusters(VoidPtr v1, VoidPtr v2)
{
   BlastClusterPtr c1, c2;

   c1 = *(BlastClusterPtr PNTR) v1;
   c2 = *(BlastClusterPtr PNTR) v2;

   if (c1->size > c2->size) return -1;
   else if (c1->size < c2->size) return 1;
   else return 0;
}

#define ORIGINAL_CLUSTER_SIZE 10
static int BlastClusterNeighbours(Int4 num_queries, Int4Ptr seq_len, 
				  CharPtr PNTR id_list, Int4Ptr gi_list)
{
   BlastClusterPtr PNTR cluster;
   Int4 num_clusters, available_size, index, i;
   Boolean numeric_id_type;

   if (gi_list)
      numeric_id_type = TRUE;
   else if (id_list)
      numeric_id_type = FALSE;
   else 
      return 0;

    for (index=1; index<num_queries; index++) {
       i = index;
       while (root[i] != i)
	  i = root[i];
       root[index] = i;
    }

    cluster = (BlastClusterPtr PNTR) 
       MemNew(num_queries*sizeof(BlastClusterPtr));
    num_clusters = 0;

    for (index=0; index<num_queries; index++) {
       if (root[index]==index) {
	  cluster[num_clusters] = 
	     (BlastClusterPtr) MemNew(sizeof(BlastCluster)); 
	  cluster[num_clusters]->size = 1;
	  cluster[num_clusters]->elements = 
	     (BlastClusterElementPtr PNTR) 
	     Malloc(ORIGINAL_CLUSTER_SIZE*sizeof(BlastClusterElementPtr));
	  cluster[num_clusters]->elements[0] = 
	     (BlastClusterElementPtr) Malloc(sizeof(BlastClusterElement));
	  cluster[num_clusters]->elements[0]->index = index;
	  cluster[num_clusters]->elements[0]->len = seq_len[index];
	  available_size = ORIGINAL_CLUSTER_SIZE;

	  for (i=index+1; i<num_queries; i++) {
	     if (root[i]==index) {
		if (cluster[num_clusters]->size >= available_size) {
		      available_size *= 2;
		      cluster[num_clusters]->elements = 
			 (BlastClusterElementPtr PNTR) 
			 Realloc(cluster[num_clusters]->elements, 
				 available_size*sizeof(BlastClusterElementPtr));
		} 
		cluster[num_clusters]->elements[cluster[num_clusters]->size] = 
		   (BlastClusterElementPtr) Malloc(sizeof(BlastClusterElement));
		cluster[num_clusters]->elements[cluster[num_clusters]->size]->index = i;
		cluster[num_clusters]->elements[cluster[num_clusters]->size]->len = seq_len[i];
		cluster[num_clusters]->size++;
	     }
	  }
	  num_clusters++;
       }
    }

    MemFree(seq_len);

    /* Sort clusters in decreasing order of sizes */
    HeapSort(cluster, num_clusters, sizeof(BlastClusterPtr), compare_clusters);
    /* Sort each cluster in decreasing order of sequence lengths */
    for (index=0; index<num_clusters; index++) {
       HeapSort(cluster[index]->elements, cluster[index]->size, 
		sizeof(BlastClusterElementPtr), compare_cluster_elements);
       for (i=0; i<cluster[index]->size; i++) {
	  if (numeric_id_type)
	     fprintf(global_fp, "%ld ", 
		     gi_list[cluster[index]->elements[i]->index]); 
	  else 
	     fprintf(global_fp, "%s ", 
		     id_list[cluster[index]->elements[i]->index]);
	  MemFree(cluster[index]->elements[i]);
       }
       fprintf(global_fp, "\n");
       MemFree(cluster[index]->elements);
       MemFree(cluster[index]);
    }
    MemFree(cluster);
    return 1;
}

#define INFO_LIST_SIZE 1000
static int ReclusterFromFile(FILE *infofp, FILE *outfp)
{
   ClusterLogInfoPtr info;
   Int4 num_queries, num_hits, i, root1, root2, total_id_len;
   Int4Ptr gi_list = NULL, seq_len;
   CharPtr id_string, PNTR id_list = NULL, ptr;
   FloatHi length_coverage, score_coverage; 
   ClusterLogHeader header;

   FileRead(&header, sizeof(ClusterLogHeader), 1, infofp);

   if (header.numeric_id_type) {
      num_queries = header.size;
      gi_list = (Int4Ptr) MemNew(num_queries*sizeof(Int4));
      FileRead(gi_list, sizeof(Int4), num_queries, infofp);
   } else {
      total_id_len = header.size;
      num_queries = 0;
      id_string = (CharPtr) MemNew(total_id_len+1);
      FileRead(id_string, sizeof(Char), total_id_len, infofp);
      ptr = id_string;
      /* Count the ID's and change delimiter from space to null character */
      for (i=0; i<total_id_len; i++) {
	 if (*ptr==' ') {
	    num_queries++;
	    *ptr = '\0';
	 }
	 ptr++;
      }
      id_list = (CharPtr PNTR) MemNew(num_queries*sizeof(CharPtr));
      for (i=0; i<num_queries; i++) {
	 /* No need for allocation of new memory */
	 id_list[i] = id_string;
	 id_string += StringLen(id_string) + 1;
      }
      id_string = id_list[0];
   }

   seq_len = (Int4Ptr) Malloc(num_queries*sizeof(Int4));
   FileRead(seq_len, sizeof(Int4), num_queries, infofp);

   info = (ClusterLogInfoPtr) MemNew(INFO_LIST_SIZE*sizeof(ClusterLogInfo));

   root = (Int4Ptr) Malloc(num_queries*sizeof(Int4));
   for (i=0; i<num_queries; i++)
      root[i] = i;

   while ((num_hits = FileRead(info, sizeof(ClusterLogInfo), INFO_LIST_SIZE,
			       infofp)) > 0) {
      for (i=0; i<num_hits; i++) {
	 if (global_parameters->bidirectional)
	    length_coverage = MIN(((FloatHi)info[i].hsp_length1) / 
				  info[i].length1, 
				  ((FloatHi)info[i].hsp_length2) / 
				  info[i].length2);
	 else
	    length_coverage = MAX(((FloatHi)info[i].hsp_length1) / 
				  info[i].length1, 
				  ((FloatHi)info[i].hsp_length2) / 
				  info[i].length2);
	 score_coverage = info[i].bit_score / 
	    (MAX(info[i].hsp_length1, info[i].hsp_length2));
	 
	 if (length_coverage >= global_parameters->length_threshold && 
	     score_coverage >= global_parameters->score_threshold) {
	    root1 = info[i].id1;
	    while (root[root1] != root1)
	       root1 = root[root1];
	    root2 = info[i].id2;
	    while (root[root2] != root2)
	       root2 = root[root2];
	    
	    if (root1 < root2)
	       root[root2] = root1;
	    else if (root1 > root2)
	       root[root1] = root2;
	 }
      } /* End loop on hits from a chunk */	 
   } /* End loop on chunks of hits */
   
   global_fp = outfp;
   BlastClusterNeighbours(num_queries, seq_len, id_list, gi_list);

   if (header.numeric_id_type)
      MemFree(gi_list);
   else {
      MemFree(id_string);
      MemFree(id_list);
   }
   MemFree(info);
   MemFree(root);
   return 0;
}

/* The following function prints only those hits which correspond to 
   an almost identical match of two sequences */
#define BUFFER_SIZE 80
static int LIBCALLBACK
PrintProteinNeighbors(VoidPtr ptr)
{
    BLAST_HitListPtr current_hitlist;
    BLAST_HSPPtr hsp; 
    Int4 index;
    BlastSearchBlkPtr search;
    Int4 id1, id2, root1, root2, hspcnt;
    Int4 high_score=0;
    Nlm_FloatHi current_evalue=DBL_MAX;
    Int4 subject_length;
    Uint1Ptr subject;
    ClusterLogInfoPtr loginfo = NULL;
    
    if (ptr == NULL)
        return 0;	
    
    search = (BlastSearchBlkPtr) ptr;
    
    if (search->current_hitlist == NULL || search->current_hitlist->hspcnt <= 0) {
        /* No hits to save. */
        search->subject_info = BLASTSubjectInfoDestruct(search->subject_info);
        return 0;
    }
    
    current_hitlist = search->current_hitlist;
    if (search->prog_number == blast_type_blastn) 
        return 0;
    
    subject_length = readdb_get_sequence(search->rdfp, search->subject_id, &subject);
    
    hspcnt = current_hitlist->hspcnt;
    
    id1 = SeqId2OrdinalId(search->rdfp, search->query_id);
    
    id2 = search->subject_id;
    
    if (id1 < id2) { /* Must be always true */
#define BUF_CHUNK_SIZE 1024
        Int4 query_length, q_length, s_length;
        FloatHi length_coverage, bit_score, score_coverage; 
	BLAST_KarlinBlkPtr kbp;

        query_length = search->context[0].query->length;

	if (global_parameters->logfp)
	   loginfo = (ClusterLogInfoPtr) MemNew(hspcnt*sizeof(ClusterLogInfo));
        
        for (index=0; index<hspcnt; index++) {
            hsp = current_hitlist->hsp_array[index];
            if (hsp) {
                /* Test if this hsp corresponds to an almost identical hit */
                q_length = hsp->query.end - hsp->query.offset;
                s_length = hsp->subject.end - hsp->subject.offset;
		
		if (global_parameters->bidirectional)
		   length_coverage = MIN(((FloatHi)q_length) / query_length, 
                                      ((FloatHi)s_length) / subject_length);
		else
		   length_coverage = MAX(((FloatHi)q_length) / query_length, 
                                      ((FloatHi)s_length) / subject_length);
		   
		if (search->pbp->gapped_calculation)
		   kbp = search->sbp->kbp_gap[search->first_context];
		else
		   kbp = search->sbp->kbp[search->first_context];
		bit_score = ((hsp->score*kbp->Lambda) -
			     kbp->logK)/NCBIMATH_LN2;

		if (global_parameters->logfp) {
		   loginfo[index].id1 = id1;
		   loginfo[index].id2 = id2;
		   loginfo[index].length1 = query_length;
		   loginfo[index].length2 = subject_length;
		   loginfo[index].hsp_length1 = q_length;
		   loginfo[index].hsp_length2 = s_length;
		   loginfo[index].bit_score = bit_score;
		}

                score_coverage = bit_score / (MAX(q_length, s_length));
                
                if (length_coverage >= global_parameters->length_threshold && 
                    score_coverage >= global_parameters->score_threshold) {
		   root1 = id1;

		   NlmMutexLockEx(&root_mutex);
		   while (root[root1] != root1)
		      root1 = root[root1];
		   root2 = id2;
		   while (root[root2] != root2)
		      root2 = root[root2];

		   if (root1 < root2)
		      root[root2] = root1;
		   else if (root1 > root2)
		      root[root1] = root2;
		   NlmMutexUnlock(root_mutex);
		}
            }
        }
	if (global_parameters->logfp) {
	   FileWrite(loginfo, sizeof(ClusterLogInfo), hspcnt, 
		     global_parameters->logfp);
	   fflush(global_parameters->logfp);
	}
    } else
       fprintf(stderr, "Error: this can't happen!\n");
    return 1;
}


#define VERSION_NO	5

#ifdef NUMARG
#undef NUMARG
#endif
#define NUMARG (sizeof(myargs)/sizeof(myargs[0]))

static Args myargs [] = {
   { "FASTA input file (program will format the database and remove files in the end)",                                                         /* 0 */
     "stdin", NULL, NULL, FALSE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
   { "Number of CPU's to use (do not use)",                      /* 1 */
      "1", NULL, NULL, FALSE, 'a', ARG_INT, 0.0, 0, NULL},       
   { "Output file for list of clusters",                         /* 2 */
      "stdout", NULL, NULL, TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
   { "Length coverage threshold",                                /* 3 */
      "0.9", NULL, NULL, FALSE, 'L', ARG_FLOAT, 0.0, 0, NULL},   
   { "Score coverage threshold",                                 /* 4 */
      "1.75", NULL, NULL, FALSE, 'S', ARG_FLOAT, 0.0, 0, NULL},  
   { "Require coverage on both neighbours?",                     /* 5 */
      "TRUE", NULL, NULL, FALSE, 'b', ARG_BOOLEAN, 0.0, 0, NULL},
   { "File to save all neighbours",                              /* 6 */
      "", NULL, NULL, TRUE, 's', ARG_FILE_OUT, 0.0, 0, NULL},
   { "File to restore neighbors for reclustering",               /* 7 */
      "", NULL, NULL, TRUE, 'r', ARG_FILE_IN, 0.0, 0, NULL}, 
   { "Input as a database",                                      /* 8 */
      "", NULL, NULL, TRUE, 'd', ARG_FILE_IN, 0.0, 0, NULL}, 
   { "Print progress messages (verbose mode)",                   /* 9 */
     "FALSE", NULL, NULL, FALSE, 'v', ARG_BOOLEAN, 0.0, 0, NULL}, 
   
};

#define MAX_DB_SIZE 100000
#define MAX_GI_LENGTH 9
#define PROGRESS_INTERVAL 1000

#ifndef TMPDIR
#define TMPDIR "/tmp"
#endif

Int2 Main (void)
{
    BLAST_OptionsBlkPtr options;
    BlastSearchBlkPtr search;
    Boolean db_is_na, query_is_na;
    Int4 qsize, dbsize;
    ReadDBFILEPtr rdfp, rdfp_var;
    Uint1 align_type;
    SeqIdPtr sip;
    SeqLocPtr query_slp = NULL;
    CharPtr blast_program, blast_inputfile, blast_outputfile, blast_database;
    CharPtr logfile, info_file, defline, input_name;
    Int4 total_id_len = 0;
    FILE *outfp, *infofp;
    Int4 index, i, num_queries;
    Int8 total_length;
    Int4Ptr gi_list, seq_len;
    CharPtr PNTR id_list, id_string;
    Boolean db_formatted = FALSE, numeric_id_type = TRUE;
    Char db_file[BUFFER_SIZE];
    ClusterLogHeader header;
    Boolean print_progress;
    FDB_optionsPtr fdb_options;
    
    if (! GetArgs ("blastclust", NUMARG, myargs))
       return (1);
    
    global_parameters = (ClusterParametersPtr) MemNew(sizeof(ClusterParameters));

    global_parameters->length_threshold = myargs[3].floatvalue;
    global_parameters->score_threshold = myargs[4].floatvalue;
    global_parameters->bidirectional = (Boolean) myargs[5].intvalue;
    
    print_progress = (Boolean) myargs[9].intvalue;
    
    blast_outputfile = myargs[2].strvalue;
    outfp = NULL;
    if (blast_outputfile != NULL) {
        if ((outfp = FileOpen(blast_outputfile, "w")) == NULL) {
            ErrPostEx(SEV_FATAL, 0, 0, "blast: Unable to open output file %s\n", blast_outputfile);
            return (1);
        }
    }

    info_file = myargs[7].strvalue;
    if (*info_file) {
       /* Non-empty string means only retrieve neighbors for reclustering */
       if ((infofp = FileOpen(info_file, "rb")) == NULL) { 
	  ErrPostEx(SEV_FATAL, 0, 0, "blast: Unable to open neighbors file %s for reading\n", logfile);
	  return (1);
       }
       /* No need for another search, simply get all the neighbours
	  and reculster them using new thresholds */
       ReclusterFromFile(infofp, outfp);
       FileClose(infofp);
       FileClose(outfp);
       return 0;
    } else
       infofp = NULL;

    blast_program = StringSave("blastp");
    if (*myargs[8].strvalue)
       db_formatted = TRUE;

    if (db_formatted) {
       blast_database = myargs[8].strvalue;
    } else {
       /* Need to format the database */
       blast_inputfile = myargs[0].strvalue;
       input_name = FileNameFind(blast_inputfile);
#ifdef TMPDIR
       blast_database = 
	  Malloc(StringLen(input_name) + StringLen(TMPDIR) + 1);
       sprintf(blast_database, "%s/%s", TMPDIR, input_name);
#else
       blast_database = blast_inputfile;
#endif

       fdb_options = MemNew(sizeof(FDB_options));
       fdb_options->db_file = StringSave(blast_inputfile);
       fdb_options->is_protein = 1; /* TRUE */
       fdb_options->parse_mode = 1; /* TRUE */
       fdb_options->base_name = StringSave(blast_database);
       FastaToBlastDB(fdb_options, blast_database, 0);

       MemFree(fdb_options->db_file);
       MemFree(fdb_options->base_name);
       MemFree(fdb_options);
    }

    logfile = myargs[6].strvalue;
    if (*logfile) { /* Empty string means do not write log information */
       if ((global_parameters->logfp = FileOpen(logfile, "wb")) == NULL) { 
	  ErrPostEx(SEV_FATAL, 0, 0, "blast: Unable to open log file %s for writing\n", logfile);
	  return (1);
       }
    }

    global_fp = outfp;

    align_type = BlastGetTypes(blast_program, &query_is_na, &db_is_na);
    
    rdfp = readdb_new(blast_database, READDB_DB_IS_PROT);
    
    options = BLASTOptionNew(blast_program, TRUE);
    
    if (options == NULL)
	return 3;
    
    options->use_real_db_size = TRUE;
    options->sort_gi_list = FALSE;
    
    options->expect_value  = 1e-6;	
    qsize = 300;
    dbsize = (20*1000*1000);
    options->searchsp_eff = (FloatHi) qsize * (FloatHi) dbsize;
    options->perform_culling = FALSE;
    options->do_not_reevaluate = TRUE;
    options->do_sum_stats = FALSE;
    options->number_of_cpus = myargs[1].intvalue;
    
    readdb_get_totals_ex(rdfp, &total_length, &num_queries, TRUE);

    root = (Int4Ptr) Malloc(num_queries*sizeof(Int4));
    for (index=0; index<num_queries; index++)
       root[index] = index;
    
    ReadDBBioseqFetchEnable ("blastclust", blast_database, db_is_na, TRUE);
       
    gi_list = (Int4Ptr) MemNew(num_queries*sizeof(Int4));
    id_list = (CharPtr PNTR) MemNew(num_queries*sizeof(CharPtr));
    seq_len = (Int4Ptr) MemNew(num_queries*sizeof(Int4));

    for (index=0; index<num_queries; index++) {
       readdb_get_descriptor(rdfp, index, &sip, &defline);
       seq_len[index] = readdb_get_sequence_length(rdfp, index);

       if (!GetAccessionFromSeqId(sip, &gi_list[index], &id_list[index])) 
	  numeric_id_type = FALSE;
       sip = SeqIdSetFree(sip);
       defline = MemFree(defline);
    }
    header.numeric_id_type = numeric_id_type;

    if (numeric_id_type) {
       id_list = MemFree(id_list);
       header.size = num_queries;
    } else {
       total_id_len = 0;
       /* Check if some ids were gis and convert them to strings */
       for (i=0; i<num_queries; i++) {
	  if (gi_list[i] > 0) {
	     id_list[i] = (CharPtr) MemNew(10);
	     sprintf(id_list[i], "%ld", gi_list[i]);
	  }
	  total_id_len += StringLen(id_list[i]) + 1;
       }
       gi_list = MemFree(gi_list);
       id_string = (CharPtr) MemNew(total_id_len+1);
       for (i=0; i<num_queries; i++) {
	  StringCat(id_string, id_list[i]);
	  StringCat(id_string, " ");
       }
       header.size = total_id_len;
    }

    FileWrite(&header, sizeof(ClusterLogHeader), 1, 
	      global_parameters->logfp);
    
    if (numeric_id_type) 
       FileWrite(gi_list, sizeof(Int4), num_queries, 
		 global_parameters->logfp);
    else {
       FileWrite(id_string, sizeof(Char), total_id_len, 
		 global_parameters->logfp);
       MemFree(id_string);
    }

    FileWrite(seq_len, sizeof(Int4), num_queries, global_parameters->logfp);
    fflush(global_parameters->logfp);

    for (index=0; index<num_queries; index++) {
       rdfp = readdb_new(blast_database, READDB_DB_IS_PROT);

       readdb_get_descriptor(rdfp, index, &sip, &defline);
       ValNodeAddPointer(&query_slp, SEQLOC_WHOLE, SeqIdDup(sip));
       SeqIdSetFree(sip);
       MemFree(defline);
       /* Set up the search */
       options->first_db_seq = index + 1;

       search = BLASTSetUpSearchWithReadDbInternal(query_slp, NULL, blast_program, seq_len[index], blast_database, options, tick_callback, NULL, NULL, 0, rdfp);
       if (search != NULL && !search->query_invalid) {
	  search->handle_results = PrintProteinNeighbors;
            
	  /* Run BLAST. */
	  search->queue_callback = NULL;
	  search->thr_info->tick_callback = tick_callback;
	  
	  do_the_blast_run(search);
       }
       search = BlastSearchBlkDestruct(search);
       query_slp = SeqLocSetFree(query_slp);
       for (rdfp_var=rdfp; rdfp_var; rdfp_var=rdfp_var->next) 
	  rdfp_var->shared_info = NULL;
       rdfp = readdb_destruct(rdfp);
       
       if (print_progress && (index + 1)%PROGRESS_INTERVAL == 0)
	  fprintf(stdout, "Finished processing of %ld queries\n", index+1);
    } /* End of loop on queries */

    BlastClusterNeighbours(num_queries, seq_len, id_list, gi_list);


    if (numeric_id_type)
       MemFree(gi_list);
    else {
       for (i=0; i<num_queries; i++)
	  MemFree(id_list[i]);
       MemFree(id_list);
    }
    MemFree(root);
    
    fflush(global_fp);
    options = BLASTOptionDelete(options);
    blast_program = (CharPtr) MemFree(blast_program);
    FileClose(outfp);
    if (global_parameters->logfp)
       FileClose(global_parameters->logfp);

    if (!db_formatted && StringLen(blast_database) > 0) {
       sprintf(db_file, "%s.phr", blast_database);
       FileRemove(db_file);
       sprintf(db_file, "%s.pin", blast_database);
       FileRemove(db_file);
       sprintf(db_file, "%s.pnd", blast_database);
       FileRemove(db_file);
       sprintf(db_file, "%s.pni", blast_database);
       FileRemove(db_file);
       sprintf(db_file, "%s.psd", blast_database);
       FileRemove(db_file);
       sprintf(db_file, "%s.psi", blast_database);
       FileRemove(db_file);
       sprintf(db_file, "%s.psq", blast_database);
       FileRemove(db_file);
    }
    return 0;
}

