#ifndef _objblst3_ 
#define _objblst3_ 

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif


#ifdef __cplusplus
extern "C" { /* } */
#endif


/**************************************************
*
*    Generated objects for Module NCBI-Blast
*    Generated using ASNCODE Revision: 6.0 at Jul 17, 1998  3:34 PM
*
**************************************************/

NLM_EXTERN Boolean LIBCALL
objblst3AsnLoad PROTO((void));
typedef ValNodePtr BlastRequestPtr;
typedef ValNode BlastRequest;
#define BlastRequest_init 1
#define BlastRequest_motd 2
#define BlastRequest_db_info 3
#define BlastRequest_db_info_specific 4
#define BlastRequest_matrix_get 5
#define BlastRequest_search 6
#define BlastRequest_db_seq_get 7
#define BlastRequest_db_redundant_ids_get 8
#define BlastRequest_db_redundant_descr_get 9
#define BlastRequest_fini 10


NLM_EXTERN BlastRequestPtr LIBCALL BlastRequestFree PROTO ((BlastRequestPtr ));
NLM_EXTERN BlastRequestPtr LIBCALL BlastRequestAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL BlastRequestAsnWrite PROTO (( BlastRequestPtr , AsnIoPtr, AsnTypePtr));

typedef ValNodePtr BlastResponsePtr;
typedef ValNode BlastResponse;
#define BlastResponse_init 1
#define BlastResponse_motd 2
#define BlastResponse_error 3
#define BlastResponse_db_seq_get 4
#define BlastResponse_db_redundant_ids_get 5
#define BlastResponse_db_redundant_descr_get 6
#define BlastResponse_db_info 7
#define BlastResponse_db_info_specific 8
#define BlastResponse_matrix 9
#define BlastResponse_alignment 10
#define BlastResponse_mask 11
#define BlastResponse_kablk 12
#define BlastResponse_parameters 13
#define BlastResponse_queued 14
#define BlastResponse_start 15
#define BlastResponse_progress 16
#define BlastResponse_done 17
#define BlastResponse_fini 18


NLM_EXTERN BlastResponsePtr LIBCALL BlastResponseFree PROTO ((BlastResponsePtr ));
NLM_EXTERN BlastResponsePtr LIBCALL BlastResponseAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL BlastResponseAsnWrite PROTO (( BlastResponsePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    BlastParameters
*
**************************************************/
typedef struct struct_Blast_parameters {
   Int4   first_threshold;
   Int4   second_threshold;
   ValNodePtr   Cutoff_cutoff;
   ValNodePtr   Cutoff2_cutoff2;
   Int4   hitlist_size;
   Int4   nucl_penalty;
   Int4   nucl_reward;
   Int4   genetic_code;
   Int4   db_genetic_code;
   Int4   low_complexity_filtering;
   Uint1   gapped_alignment;
   Int4   gap_open;
   Int4   gap_extend;
   Int4   required_start;
   Int4   required_end;
   FloatHi   ethresh;
   Int4   max_num_passes;
   Int4   pseudo_count_const;
   CharPtr   other_options;
   ValNodePtr   gilist;
   CharPtr   gifile;
   CharPtr   matrix;
   CharPtr   filter_string;
} BlastParameters, PNTR BlastParametersPtr;


NLM_EXTERN BlastParametersPtr LIBCALL BlastParametersFree PROTO ((BlastParametersPtr ));
NLM_EXTERN BlastParametersPtr LIBCALL BlastParametersNew PROTO (( void ));
NLM_EXTERN BlastParametersPtr LIBCALL BlastParametersAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL BlastParametersAsnWrite PROTO (( BlastParametersPtr , AsnIoPtr, AsnTypePtr));


#ifdef NLM_GENERATED_CODE_PROTO

typedef ValNodePtr Cutoff2_cutoff2Ptr;
typedef ValNode Cutoff2_cutoff2;

#endif /* NLM_GENERATED_CODE_PROTO */

#define Cutoff2_cutoff2_evalue 1
#define Cutoff2_cutoff2_score 2

#ifdef NLM_GENERATED_CODE_PROTO

static Cutoff2_cutoff2Ptr LIBCALL Cutoff2_cutoff2Free PROTO ((Cutoff2_cutoff2Ptr ));
static Cutoff2_cutoff2Ptr LIBCALL Cutoff2_cutoff2AsnRead PROTO (( AsnIoPtr, AsnTypePtr));
static Boolean LIBCALL Cutoff2_cutoff2AsnWrite PROTO (( Cutoff2_cutoff2Ptr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */


#ifdef NLM_GENERATED_CODE_PROTO

typedef ValNodePtr Cutoff_cutoffPtr;
typedef ValNode Cutoff_cutoff;

#endif /* NLM_GENERATED_CODE_PROTO */

#define Cutoff_cutoff_evalue 1
#define Cutoff_cutoff_score 2

#ifdef NLM_GENERATED_CODE_PROTO

static Cutoff_cutoffPtr LIBCALL Cutoff_cutoffFree PROTO ((Cutoff_cutoffPtr ));
static Cutoff_cutoffPtr LIBCALL Cutoff_cutoffAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
static Boolean LIBCALL Cutoff_cutoffAsnWrite PROTO (( Cutoff_cutoffPtr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */



/**************************************************
*
*    BlastDbinfoGet
*
**************************************************/
typedef struct struct_Blast_dbinfo_get {
   CharPtr   name;
   Uint2   type;
   /* following #defines are for enumerated type, not used by object loaders */
#define Blast_dbinfo_get_type_unknown 0
#define Blast_dbinfo_get_type_protein 1
#define Blast_dbinfo_get_type_nucleotide 2

} BlastDbinfoGet, PNTR BlastDbinfoGetPtr;


NLM_EXTERN BlastDbinfoGetPtr LIBCALL BlastDbinfoGetFree PROTO ((BlastDbinfoGetPtr ));
NLM_EXTERN BlastDbinfoGetPtr LIBCALL BlastDbinfoGetNew PROTO (( void ));
NLM_EXTERN BlastDbinfoGetPtr LIBCALL BlastDbinfoGetAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL BlastDbinfoGetAsnWrite PROTO (( BlastDbinfoGetPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    BlastSearch
*
**************************************************/
typedef struct struct_Blast_search {
   Uint2   program;
   /* following #defines are for enumerated type, not used by object loaders */
#define Blast_search_program_blastn 0
#define Blast_search_program_blastp 1
#define Blast_search_program_blastx 2
#define Blast_search_program_tblastn 3
#define Blast_search_program_tblastx 4

   struct struct_Bioseq PNTR   query;
   CharPtr   database;
   struct struct_Blast_parameters PNTR   parameters;
   ValNodePtr   mask;
   struct struct_Blast_matrix PNTR   matrix;
} BlastSearch, PNTR BlastSearchPtr;


NLM_EXTERN BlastSearchPtr LIBCALL BlastSearchFree PROTO ((BlastSearchPtr ));
NLM_EXTERN BlastSearchPtr LIBCALL BlastSearchNew PROTO (( void ));
NLM_EXTERN BlastSearchPtr LIBCALL BlastSearchAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL BlastSearchAsnWrite PROTO (( BlastSearchPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    BlastSeqId
*
**************************************************/
typedef struct struct_Blast_seq_id {
   Uint1   is_protein;
   CharPtr   database;
   ValNodePtr   id;
} BlastSeqId, PNTR BlastSeqIdPtr;


NLM_EXTERN BlastSeqIdPtr LIBCALL BlastSeqIdFree PROTO ((BlastSeqIdPtr ));
NLM_EXTERN BlastSeqIdPtr LIBCALL BlastSeqIdNew PROTO (( void ));
NLM_EXTERN BlastSeqIdPtr LIBCALL BlastSeqIdAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL BlastSeqIdAsnWrite PROTO (( BlastSeqIdPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    BlastMatrix
*
**************************************************/
typedef struct struct_Blast_matrix {
   Uint1   is_protein;
   CharPtr   name;
   ValNodePtr   comments;
   Int4   row_length;
   Int4   column_length;
   ValNodePtr   scores;
   FloatHi   karlinK;
} BlastMatrix, PNTR BlastMatrixPtr;


NLM_EXTERN BlastMatrixPtr LIBCALL BlastMatrixFree PROTO ((BlastMatrixPtr ));
NLM_EXTERN BlastMatrixPtr LIBCALL BlastMatrixNew PROTO (( void ));
NLM_EXTERN BlastMatrixPtr LIBCALL BlastMatrixAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL BlastMatrixAsnWrite PROTO (( BlastMatrixPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    BlastDbinfo
*
**************************************************/
typedef struct struct_Blast_dbinfo {
   struct struct_Blast_dbinfo PNTR next;
   Uint1   is_protein;
   CharPtr   name;
   CharPtr   definition;
   CharPtr   date;
   Int4   total_length;
   Int4   number_seqs;
} BlastDbinfo, PNTR BlastDbinfoPtr;


NLM_EXTERN BlastDbinfoPtr LIBCALL BlastDbinfoFree PROTO ((BlastDbinfoPtr ));
NLM_EXTERN BlastDbinfoPtr LIBCALL BlastDbinfoNew PROTO (( void ));
NLM_EXTERN BlastDbinfoPtr LIBCALL BlastDbinfoAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL BlastDbinfoAsnWrite PROTO (( BlastDbinfoPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    BlastQueued
*
**************************************************/
typedef struct struct_Blast_Queued {
   Int4   length;
} BlastQueued, PNTR BlastQueuedPtr;


NLM_EXTERN BlastQueuedPtr LIBCALL BlastQueuedFree PROTO ((BlastQueuedPtr ));
NLM_EXTERN BlastQueuedPtr LIBCALL BlastQueuedNew PROTO (( void ));
NLM_EXTERN BlastQueuedPtr LIBCALL BlastQueuedAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL BlastQueuedAsnWrite PROTO (( BlastQueuedPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    BlastProgress
*
**************************************************/
typedef struct struct_Blast_Progress {
   Int4   completed;
} BlastProgress, PNTR BlastProgressPtr;


NLM_EXTERN BlastProgressPtr LIBCALL BlastProgressFree PROTO ((BlastProgressPtr ));
NLM_EXTERN BlastProgressPtr LIBCALL BlastProgressNew PROTO (( void ));
NLM_EXTERN BlastProgressPtr LIBCALL BlastProgressAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL BlastProgressAsnWrite PROTO (( BlastProgressPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    BlastKABlk
*
**************************************************/
typedef struct struct_Blast_KABlk {
   FloatHi   lambda;
   FloatHi   k;
   FloatHi   h;
   Uint1   gapped;
} BlastKABlk, PNTR BlastKABlkPtr;


NLM_EXTERN BlastKABlkPtr LIBCALL BlastKABlkFree PROTO ((BlastKABlkPtr ));
NLM_EXTERN BlastKABlkPtr LIBCALL BlastKABlkNew PROTO (( void ));
NLM_EXTERN BlastKABlkPtr LIBCALL BlastKABlkAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL BlastKABlkAsnWrite PROTO (( BlastKABlkPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    BlastDefline
*
**************************************************/
typedef struct struct_Blast_defline {
   struct struct_Blast_defline PNTR next;
   ValNodePtr   id;
   CharPtr   defline;
} BlastDefline, PNTR BlastDeflinePtr;


NLM_EXTERN BlastDeflinePtr LIBCALL BlastDeflineFree PROTO ((BlastDeflinePtr ));
NLM_EXTERN BlastDeflinePtr LIBCALL BlastDeflineNew PROTO (( void ));
NLM_EXTERN BlastDeflinePtr LIBCALL BlastDeflineAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL BlastDeflineAsnWrite PROTO (( BlastDeflinePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    BlastMask
*
**************************************************/
typedef struct struct_Blast_mask {
   ValNodePtr   location;
   Uint2   frame;
   /* following #defines are for enumerated type, not used by object loaders */
#define Blast_mask_frame_notset 0
#define Blast_mask_frame_plus1 1
#define Blast_mask_frame_plus2 2
#define Blast_mask_frame_plus3 3
#define Blast_mask_frame_minus1 4
#define Blast_mask_frame_minus2 5
#define Blast_mask_frame_minus3 6

} BlastMask, PNTR BlastMaskPtr;


NLM_EXTERN BlastMaskPtr LIBCALL BlastMaskFree PROTO ((BlastMaskPtr ));
NLM_EXTERN BlastMaskPtr LIBCALL BlastMaskNew PROTO (( void ));
NLM_EXTERN BlastMaskPtr LIBCALL BlastMaskAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL BlastMaskAsnWrite PROTO (( BlastMaskPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    BlastVersion
*
**************************************************/
typedef struct struct_Blast_version {
   CharPtr   version;
   CharPtr   date;
} BlastVersion, PNTR BlastVersionPtr;


NLM_EXTERN BlastVersionPtr LIBCALL BlastVersionFree PROTO ((BlastVersionPtr ));
NLM_EXTERN BlastVersionPtr LIBCALL BlastVersionNew PROTO (( void ));
NLM_EXTERN BlastVersionPtr LIBCALL BlastVersionAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL BlastVersionAsnWrite PROTO (( BlastVersionPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    BlastError
*
**************************************************/
typedef struct struct_Blast_error {
   Uint2   level;
   /* following #defines are for enumerated type, not used by object loaders */
#define Blast_error_level_none 0
#define Blast_error_level_info 1
#define Blast_error_level_warn 2
#define Blast_error_level_error 3
#define Blast_error_level_fatal 4

   CharPtr   msg;
} BlastError, PNTR BlastErrorPtr;


NLM_EXTERN BlastErrorPtr LIBCALL BlastErrorFree PROTO ((BlastErrorPtr ));
NLM_EXTERN BlastErrorPtr LIBCALL BlastErrorNew PROTO (( void ));
NLM_EXTERN BlastErrorPtr LIBCALL BlastErrorAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL BlastErrorAsnWrite PROTO (( BlastErrorPtr , AsnIoPtr, AsnTypePtr));

#ifdef __cplusplus
/* { */ }
#endif

#endif /* _objblst3_ */

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

