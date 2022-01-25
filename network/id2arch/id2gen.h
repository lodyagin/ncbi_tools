#ifndef _id2gen_ 
#define _id2gen_ 

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
*    Generated objects for Module NCBI-ID2Access
*    Generated using ASNCODE Revision: 6.0 at Dec 15, 2003  5:08 PM
*
**************************************************/

NLM_EXTERN Boolean LIBCALL
id2genAsnLoad PROTO((void));


/**************************************************
*
*    ID2RequestPacket
*
**************************************************/
typedef struct struct_ID2Request ID2RequestPacket;
typedef struct struct_ID2Request PNTR ID2RequestPacketPtr;
#define ID2RequestPacketNew() ID2RequestNew() 

#ifdef NLM_GENERATED_CODE_PROTO

NLM_EXTERN ID2RequestPacketPtr LIBCALL ID2RequestPacketFree PROTO ((ID2RequestPacketPtr ));
NLM_EXTERN ID2RequestPacketPtr LIBCALL ID2RequestPacketNew PROTO (( void ));
NLM_EXTERN ID2RequestPacketPtr LIBCALL ID2RequestPacketAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2RequestPacketAsnWrite PROTO (( ID2RequestPacketPtr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */



/**************************************************
*
*    ID2Request
*
**************************************************/
typedef struct struct_ID2_Request {
   struct struct_ID2_Request PNTR next;
   Int4   serial_number;
   struct struct_ID2_Param PNTR   params;
   ValNodePtr   Request_request;
} ID2Request, PNTR ID2RequestPtr;


NLM_EXTERN ID2RequestPtr LIBCALL ID2RequestFree PROTO ((ID2RequestPtr ));
NLM_EXTERN ID2RequestPtr LIBCALL ID2RequestNew PROTO (( void ));
NLM_EXTERN ID2RequestPtr LIBCALL ID2RequestAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2RequestAsnWrite PROTO (( ID2RequestPtr , AsnIoPtr, AsnTypePtr));


#ifdef NLM_GENERATED_CODE_PROTO

typedef ValNodePtr Request_requestPtr;
typedef ValNode Request_request;

#endif /* NLM_GENERATED_CODE_PROTO */

#define Request_request_init 1
#define Request_request_get_packages 2
#define Request_request_string_to_gi 3
#define Request_request_seq_id_to_gi 4
#define Request_request_gi_to_tse_id 5
#define Request_request_get_tse 6
#define Request_request_reget_tse 7
#define Request_request_get_chunks 8

#ifdef NLM_GENERATED_CODE_PROTO

static Request_requestPtr LIBCALL Request_requestFree PROTO ((Request_requestPtr ));
static Request_requestPtr LIBCALL Request_requestAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
static Boolean LIBCALL Request_requestAsnWrite PROTO (( Request_requestPtr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */



/**************************************************
*
*    ID2Params
*
**************************************************/
typedef struct struct_ID2Param ID2Params;
typedef struct struct_ID2Param PNTR ID2ParamsPtr;
#define ID2ParamsNew() ID2ParamNew() 

#ifdef NLM_GENERATED_CODE_PROTO

NLM_EXTERN ID2ParamsPtr LIBCALL ID2ParamsFree PROTO ((ID2ParamsPtr ));
NLM_EXTERN ID2ParamsPtr LIBCALL ID2ParamsNew PROTO (( void ));
NLM_EXTERN ID2ParamsPtr LIBCALL ID2ParamsAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2ParamsAsnWrite PROTO (( ID2ParamsPtr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */



/**************************************************
*
*    ID2RequestGetPackages
*
**************************************************/
typedef struct struct_ID2_Request_Get_Packages {
   ValNodePtr   names;
   Uint1   no_contents;
} ID2RequestGetPackages, PNTR ID2RequestGetPackagesPtr;


NLM_EXTERN ID2RequestGetPackagesPtr LIBCALL ID2RequestGetPackagesFree PROTO ((ID2RequestGetPackagesPtr ));
NLM_EXTERN ID2RequestGetPackagesPtr LIBCALL ID2RequestGetPackagesNew PROTO (( void ));
NLM_EXTERN ID2RequestGetPackagesPtr LIBCALL ID2RequestGetPackagesAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2RequestGetPackagesAsnWrite PROTO (( ID2RequestGetPackagesPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ID2RequestStringToGi
*
**************************************************/
typedef struct struct_ID2_Request_String_To_Gi {
   CharPtr   id;
} ID2RequestStringToGi, PNTR ID2RequestStringToGiPtr;


NLM_EXTERN ID2RequestStringToGiPtr LIBCALL ID2RequestStringToGiFree PROTO ((ID2RequestStringToGiPtr ));
NLM_EXTERN ID2RequestStringToGiPtr LIBCALL ID2RequestStringToGiNew PROTO (( void ));
NLM_EXTERN ID2RequestStringToGiPtr LIBCALL ID2RequestStringToGiAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2RequestStringToGiAsnWrite PROTO (( ID2RequestStringToGiPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ID2RequestSeqIdToGi
*
**************************************************/
typedef struct struct_ID2_Request_Seq_id_To_Gi {
   ValNodePtr   seq_id;
} ID2RequestSeqIdToGi, PNTR ID2RequestSeqIdToGiPtr;


NLM_EXTERN ID2RequestSeqIdToGiPtr LIBCALL ID2RequestSeqIdToGiFree PROTO ((ID2RequestSeqIdToGiPtr ));
NLM_EXTERN ID2RequestSeqIdToGiPtr LIBCALL ID2RequestSeqIdToGiNew PROTO (( void ));
NLM_EXTERN ID2RequestSeqIdToGiPtr LIBCALL ID2RequestSeqIdToGiAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2RequestSeqIdToGiAsnWrite PROTO (( ID2RequestSeqIdToGiPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ID2RequestGiToTSEId
*
**************************************************/
typedef struct struct_ID2_Request_Gi_To_TSE_Id {
   ValNodePtr   Gi_gi;
   ValNodePtr   sources;
   Uint1   external;
   Uint1   current_gis;
} ID2RequestGiToTSEId, PNTR ID2RequestGiToTSEIdPtr;


NLM_EXTERN ID2RequestGiToTSEIdPtr LIBCALL ID2RequestGiToTSEIdFree PROTO ((ID2RequestGiToTSEIdPtr ));
NLM_EXTERN ID2RequestGiToTSEIdPtr LIBCALL ID2RequestGiToTSEIdNew PROTO (( void ));
NLM_EXTERN ID2RequestGiToTSEIdPtr LIBCALL ID2RequestGiToTSEIdAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2RequestGiToTSEIdAsnWrite PROTO (( ID2RequestGiToTSEIdPtr , AsnIoPtr, AsnTypePtr));


#ifdef NLM_GENERATED_CODE_PROTO

typedef ValNodePtr Gi_giPtr;
typedef ValNode Gi_gi;

#endif /* NLM_GENERATED_CODE_PROTO */

#define Gi_gi_gi 1
#define Gi_gi_string 2
#define Gi_gi_seq_id 3

#ifdef NLM_GENERATED_CODE_PROTO

static Gi_giPtr LIBCALL Gi_giFree PROTO ((Gi_giPtr ));
static Gi_giPtr LIBCALL Gi_giAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
static Boolean LIBCALL Gi_giAsnWrite PROTO (( Gi_giPtr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */



/**************************************************
*
*    ID2RequestGetTSE
*
**************************************************/
typedef struct struct_ID2_Request_Get_TSE {
   ValNodePtr   TseId_tse_id;
   struct struct_ID2_Get_TSE_Details PNTR   details;
} ID2RequestGetTSE, PNTR ID2RequestGetTSEPtr;


NLM_EXTERN ID2RequestGetTSEPtr LIBCALL ID2RequestGetTSEFree PROTO ((ID2RequestGetTSEPtr ));
NLM_EXTERN ID2RequestGetTSEPtr LIBCALL ID2RequestGetTSENew PROTO (( void ));
NLM_EXTERN ID2RequestGetTSEPtr LIBCALL ID2RequestGetTSEAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2RequestGetTSEAsnWrite PROTO (( ID2RequestGetTSEPtr , AsnIoPtr, AsnTypePtr));


#ifdef NLM_GENERATED_CODE_PROTO

typedef ValNodePtr TseId_tse_idPtr;
typedef ValNode TseId_tse_id;

#endif /* NLM_GENERATED_CODE_PROTO */

#define TseId_tse_id_tse_id 1
#define TseId_tse_id_TseId_Gi 2

#ifdef NLM_GENERATED_CODE_PROTO

static TseId_tse_idPtr LIBCALL TseId_tse_idFree PROTO ((TseId_tse_idPtr ));
static TseId_tse_idPtr LIBCALL TseId_tse_idAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
static Boolean LIBCALL TseId_tse_idAsnWrite PROTO (( TseId_tse_idPtr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */



/**************************************************
*
*    TseId_gi
*
**************************************************/

#ifdef NLM_GENERATED_CODE_PROTO

typedef struct struct_TseId_Gi {
   struct struct_ID2_Request_Gi_To_TSE_Id PNTR   request;
   struct struct_ID2_TSE_Id PNTR   exclude_tses;
} TseId_gi, PNTR TseId_giPtr;
#endif /* NLM_GENERATED_CODE_PROTO */

#ifdef NLM_GENERATED_CODE_PROTO

static TseId_giPtr LIBCALL TseId_giFree PROTO ((TseId_giPtr ));
static TseId_giPtr LIBCALL TseId_giNew PROTO (( void ));
static TseId_giPtr LIBCALL TseId_giAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
static Boolean LIBCALL TseId_giAsnWrite PROTO (( TseId_giPtr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */



/**************************************************
*
*    ID2RequestReGetTSE
*
**************************************************/
typedef struct struct_ID2_Request_ReGet_TSE {
   struct struct_ID2_TSE_Id PNTR   tse_id;
   struct struct_ID2_Get_TSE_Details PNTR   details;
   Int4   offset;
} ID2RequestReGetTSE, PNTR ID2RequestReGetTSEPtr;


NLM_EXTERN ID2RequestReGetTSEPtr LIBCALL ID2RequestReGetTSEFree PROTO ((ID2RequestReGetTSEPtr ));
NLM_EXTERN ID2RequestReGetTSEPtr LIBCALL ID2RequestReGetTSENew PROTO (( void ));
NLM_EXTERN ID2RequestReGetTSEPtr LIBCALL ID2RequestReGetTSEAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2RequestReGetTSEAsnWrite PROTO (( ID2RequestReGetTSEPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ID2SRequestGetChunks
*
**************************************************/
typedef struct struct_ID2S_Request_Get_Chunks {
   struct struct_ID2_TSE_Id PNTR   tse_id;
   ValNodePtr   chunks;
} ID2SRequestGetChunks, PNTR ID2SRequestGetChunksPtr;


NLM_EXTERN ID2SRequestGetChunksPtr LIBCALL ID2SRequestGetChunksFree PROTO ((ID2SRequestGetChunksPtr ));
NLM_EXTERN ID2SRequestGetChunksPtr LIBCALL ID2SRequestGetChunksNew PROTO (( void ));
NLM_EXTERN ID2SRequestGetChunksPtr LIBCALL ID2SRequestGetChunksAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2SRequestGetChunksAsnWrite PROTO (( ID2SRequestGetChunksPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ID2TSEId
*
**************************************************/
typedef struct struct_ID2_TSE_Id {
   struct struct_ID2_TSE_Id PNTR next;
   Int4   sat;
   Int4   sat_key;
} ID2TSEId, PNTR ID2TSEIdPtr;


NLM_EXTERN ID2TSEIdPtr LIBCALL ID2TSEIdFree PROTO ((ID2TSEIdPtr ));
NLM_EXTERN ID2TSEIdPtr LIBCALL ID2TSEIdNew PROTO (( void ));
NLM_EXTERN ID2TSEIdPtr LIBCALL ID2TSEIdAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2TSEIdAsnWrite PROTO (( ID2TSEIdPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ID2GetTSEDetails
*
**************************************************/
typedef struct struct_ID2_Get_TSE_Details {
   ValNodePtr   location;
   Int4   seq_class_level;
   Int4   descr_level;
   Int4   descr_type_mask;
   Int4   annot_type_mask;
   Int4   feat_type_mask;
   Uint2   sequence_level;
   /* following #defines are for enumerated type, not used by object loaders */
#define ID2_Get_TSE_Details_sequence_level_none 0
#define ID2_Get_TSE_Details_sequence_level_seq_map 1
#define ID2_Get_TSE_Details_sequence_level_seq_data 2

} ID2GetTSEDetails, PNTR ID2GetTSEDetailsPtr;


NLM_EXTERN ID2GetTSEDetailsPtr LIBCALL ID2GetTSEDetailsFree PROTO ((ID2GetTSEDetailsPtr ));
NLM_EXTERN ID2GetTSEDetailsPtr LIBCALL ID2GetTSEDetailsNew PROTO (( void ));
NLM_EXTERN ID2GetTSEDetailsPtr LIBCALL ID2GetTSEDetailsAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2GetTSEDetailsAsnWrite PROTO (( ID2GetTSEDetailsPtr , AsnIoPtr, AsnTypePtr));

typedef ValNodePtr ID2SeqLocPtr;
typedef ValNode ID2SeqLoc;
#define ID2SeqLoc_whole 1
#define ID2SeqLoc_int__ 2
#define ID2SeqLoc_int_set 3
#define ID2SeqLoc_whole_range 4
#define ID2SeqLoc_loc_set 5


NLM_EXTERN ID2SeqLocPtr LIBCALL ID2SeqLocFree PROTO ((ID2SeqLocPtr ));
NLM_EXTERN ID2SeqLocPtr LIBCALL ID2SeqLocAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2SeqLocAsnWrite PROTO (( ID2SeqLocPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ID2Reply
*
**************************************************/
typedef struct struct_ID2_Reply {
   Int4   serial_number;
   struct struct_ID2_Param PNTR   params;
   ValNodePtr   Reply_reply;
   struct struct_ID2_Error PNTR   error;
   Uint1   end_of_reply;
} ID2Reply, PNTR ID2ReplyPtr;


NLM_EXTERN ID2ReplyPtr LIBCALL ID2ReplyFree PROTO ((ID2ReplyPtr ));
NLM_EXTERN ID2ReplyPtr LIBCALL ID2ReplyNew PROTO (( void ));
NLM_EXTERN ID2ReplyPtr LIBCALL ID2ReplyAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2ReplyAsnWrite PROTO (( ID2ReplyPtr , AsnIoPtr, AsnTypePtr));


#ifdef NLM_GENERATED_CODE_PROTO

typedef ValNodePtr Reply_replyPtr;
typedef ValNode Reply_reply;

#endif /* NLM_GENERATED_CODE_PROTO */

#define Reply_reply_init 1
#define Reply_reply_get_package 2
#define Reply_reply_seq_id_to_gi 3
#define Reply_reply_gi_to_tse_id 4
#define Reply_reply_get_tse 5
#define Reply_reply_get_tse_info 6
#define Reply_reply_get_chunk 7

#ifdef NLM_GENERATED_CODE_PROTO

static Reply_replyPtr LIBCALL Reply_replyFree PROTO ((Reply_replyPtr ));
static Reply_replyPtr LIBCALL Reply_replyAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
static Boolean LIBCALL Reply_replyAsnWrite PROTO (( Reply_replyPtr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */



/**************************************************
*
*    ID2ReplyGetPackage
*
**************************************************/
typedef struct struct_ID2_Reply_Get_Package {
   CharPtr   name;
   struct struct_ID2_Param PNTR   params;
} ID2ReplyGetPackage, PNTR ID2ReplyGetPackagePtr;


NLM_EXTERN ID2ReplyGetPackagePtr LIBCALL ID2ReplyGetPackageFree PROTO ((ID2ReplyGetPackagePtr ));
NLM_EXTERN ID2ReplyGetPackagePtr LIBCALL ID2ReplyGetPackageNew PROTO (( void ));
NLM_EXTERN ID2ReplyGetPackagePtr LIBCALL ID2ReplyGetPackageAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2ReplyGetPackageAsnWrite PROTO (( ID2ReplyGetPackagePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ID2ReplySeqIdToGi
*
**************************************************/
typedef struct struct_ID2_Reply_Seq_id_To_Gi {
   ValNodePtr   seq_id;
   Int4   gi;
} ID2ReplySeqIdToGi, PNTR ID2ReplySeqIdToGiPtr;


NLM_EXTERN ID2ReplySeqIdToGiPtr LIBCALL ID2ReplySeqIdToGiFree PROTO ((ID2ReplySeqIdToGiPtr ));
NLM_EXTERN ID2ReplySeqIdToGiPtr LIBCALL ID2ReplySeqIdToGiNew PROTO (( void ));
NLM_EXTERN ID2ReplySeqIdToGiPtr LIBCALL ID2ReplySeqIdToGiAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2ReplySeqIdToGiAsnWrite PROTO (( ID2ReplySeqIdToGiPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ID2ReplyGiToTSEId
*
**************************************************/
typedef struct struct_ID2_Reply_Gi_To_TSE_Id {
   Int4   gi;
   CharPtr   source;
   struct struct_ID2_TSE_Id_Info PNTR   tses;
} ID2ReplyGiToTSEId, PNTR ID2ReplyGiToTSEIdPtr;


NLM_EXTERN ID2ReplyGiToTSEIdPtr LIBCALL ID2ReplyGiToTSEIdFree PROTO ((ID2ReplyGiToTSEIdPtr ));
NLM_EXTERN ID2ReplyGiToTSEIdPtr LIBCALL ID2ReplyGiToTSEIdNew PROTO (( void ));
NLM_EXTERN ID2ReplyGiToTSEIdPtr LIBCALL ID2ReplyGiToTSEIdAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2ReplyGiToTSEIdAsnWrite PROTO (( ID2ReplyGiToTSEIdPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ID2ReplyGetTSE
*
**************************************************/
typedef struct struct_ID2_Reply_Get_TSE {
   struct struct_ID2_TSE_Id PNTR   tse_id;
   struct struct_ID2_Reply_Data PNTR   data;
} ID2ReplyGetTSE, PNTR ID2ReplyGetTSEPtr;


NLM_EXTERN ID2ReplyGetTSEPtr LIBCALL ID2ReplyGetTSEFree PROTO ((ID2ReplyGetTSEPtr ));
NLM_EXTERN ID2ReplyGetTSEPtr LIBCALL ID2ReplyGetTSENew PROTO (( void ));
NLM_EXTERN ID2ReplyGetTSEPtr LIBCALL ID2ReplyGetTSEAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2ReplyGetTSEAsnWrite PROTO (( ID2ReplyGetTSEPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ID2SReplyGetTSEInfo
*
**************************************************/
typedef struct struct_ID2S_Reply_Get_TSE_Info {
   struct struct_ID2_TSE_Id PNTR   tse_id;
   Int4   split_version;
   struct struct_ID2_Reply_Data PNTR   info;
} ID2SReplyGetTSEInfo, PNTR ID2SReplyGetTSEInfoPtr;


NLM_EXTERN ID2SReplyGetTSEInfoPtr LIBCALL ID2SReplyGetTSEInfoFree PROTO ((ID2SReplyGetTSEInfoPtr ));
NLM_EXTERN ID2SReplyGetTSEInfoPtr LIBCALL ID2SReplyGetTSEInfoNew PROTO (( void ));
NLM_EXTERN ID2SReplyGetTSEInfoPtr LIBCALL ID2SReplyGetTSEInfoAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2SReplyGetTSEInfoAsnWrite PROTO (( ID2SReplyGetTSEInfoPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ID2SReplyGetChunk
*
**************************************************/
typedef struct struct_ID2S_Reply_Get_Chunk {
   struct struct_ID2_TSE_Id PNTR   tse_id;
   Int4   chunk_id;
   struct struct_ID2_Reply_Data PNTR   data;
} ID2SReplyGetChunk, PNTR ID2SReplyGetChunkPtr;


NLM_EXTERN ID2SReplyGetChunkPtr LIBCALL ID2SReplyGetChunkFree PROTO ((ID2SReplyGetChunkPtr ));
NLM_EXTERN ID2SReplyGetChunkPtr LIBCALL ID2SReplyGetChunkNew PROTO (( void ));
NLM_EXTERN ID2SReplyGetChunkPtr LIBCALL ID2SReplyGetChunkAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2SReplyGetChunkAsnWrite PROTO (( ID2SReplyGetChunkPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ID2Error
*
**************************************************/
typedef struct struct_ID2_Error {
   struct struct_ID2_Error PNTR next;
   Uint2   severity;
   /* following #defines are for enumerated type, not used by object loaders */
#define ID2_Error_severity_warning 1
#define ID2_Error_severity_failed_command 2
#define ID2_Error_severity_failed_connection 3
#define ID2_Error_severity_failed_server 4
#define ID2_Error_severity_no_data 5
#define ID2_Error_severity_restricted_data 6
#define ID2_Error_severity_unsupported_command 7
#define ID2_Error_severity_invalid_arguments 8

   Int4   retry_delay;
   CharPtr   message;
} ID2Error, PNTR ID2ErrorPtr;


NLM_EXTERN ID2ErrorPtr LIBCALL ID2ErrorFree PROTO ((ID2ErrorPtr ));
NLM_EXTERN ID2ErrorPtr LIBCALL ID2ErrorNew PROTO (( void ));
NLM_EXTERN ID2ErrorPtr LIBCALL ID2ErrorAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2ErrorAsnWrite PROTO (( ID2ErrorPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ID2TSEIdInfo
*
**************************************************/
typedef struct struct_ID2_TSE_Id_Info {
   struct struct_ID2_TSE_Id_Info PNTR next;
   struct struct_ID2_TSE_Id PNTR   tse_id;
   Int4   split_version;
} ID2TSEIdInfo, PNTR ID2TSEIdInfoPtr;


NLM_EXTERN ID2TSEIdInfoPtr LIBCALL ID2TSEIdInfoFree PROTO ((ID2TSEIdInfoPtr ));
NLM_EXTERN ID2TSEIdInfoPtr LIBCALL ID2TSEIdInfoNew PROTO (( void ));
NLM_EXTERN ID2TSEIdInfoPtr LIBCALL ID2TSEIdInfoAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2TSEIdInfoAsnWrite PROTO (( ID2TSEIdInfoPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ID2ReplyData
*
**************************************************/
typedef struct struct_ID2_Reply_Data {
   Int4   data_type;
   Int4   data_format;
   Int4   data_compression;
   ValNodePtr   data;
} ID2ReplyData, PNTR ID2ReplyDataPtr;


NLM_EXTERN ID2ReplyDataPtr LIBCALL ID2ReplyDataFree PROTO ((ID2ReplyDataPtr ));
NLM_EXTERN ID2ReplyDataPtr LIBCALL ID2ReplyDataNew PROTO (( void ));
NLM_EXTERN ID2ReplyDataPtr LIBCALL ID2ReplyDataAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2ReplyDataAsnWrite PROTO (( ID2ReplyDataPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ID2SSplitInfo
*
**************************************************/
typedef struct struct_ID2S_Split_Info {
   struct struct_ID2S_Bioseqs_Info PNTR   bioseqs_info;
   struct struct_ID2S_Chunk_Info PNTR   chunks;
} ID2SSplitInfo, PNTR ID2SSplitInfoPtr;


NLM_EXTERN ID2SSplitInfoPtr LIBCALL ID2SSplitInfoFree PROTO ((ID2SSplitInfoPtr ));
NLM_EXTERN ID2SSplitInfoPtr LIBCALL ID2SSplitInfoNew PROTO (( void ));
NLM_EXTERN ID2SSplitInfoPtr LIBCALL ID2SSplitInfoAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2SSplitInfoAsnWrite PROTO (( ID2SSplitInfoPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ID2SBioseqsInfo
*
**************************************************/
typedef struct struct_ID2S_Bioseqs_Info {
   struct struct_ID2S_Bioseqs_Info PNTR next;
   struct struct_ID2S_Bioseq_Info PNTR   info;
   struct struct_ID2_Id_Range PNTR   bioseqs;
} ID2SBioseqsInfo, PNTR ID2SBioseqsInfoPtr;


NLM_EXTERN ID2SBioseqsInfoPtr LIBCALL ID2SBioseqsInfoFree PROTO ((ID2SBioseqsInfoPtr ));
NLM_EXTERN ID2SBioseqsInfoPtr LIBCALL ID2SBioseqsInfoNew PROTO (( void ));
NLM_EXTERN ID2SBioseqsInfoPtr LIBCALL ID2SBioseqsInfoAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2SBioseqsInfoAsnWrite PROTO (( ID2SBioseqsInfoPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ID2SChunkInfo
*
**************************************************/
typedef struct struct_ID2S_Chunk_Info {
   struct struct_ID2S_Chunk_Info PNTR next;
   Int4   id;
   ValNodePtr   content;
} ID2SChunkInfo, PNTR ID2SChunkInfoPtr;


NLM_EXTERN ID2SChunkInfoPtr LIBCALL ID2SChunkInfoFree PROTO ((ID2SChunkInfoPtr ));
NLM_EXTERN ID2SChunkInfoPtr LIBCALL ID2SChunkInfoNew PROTO (( void ));
NLM_EXTERN ID2SChunkInfoPtr LIBCALL ID2SChunkInfoAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2SChunkInfoAsnWrite PROTO (( ID2SChunkInfoPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ID2SBioseqInfo
*
**************************************************/
typedef struct struct_ID2S_Bioseq_Info {
   Int4   gap_count;
   Uint1   seq_map_has_ref;
   struct struct_ID2S_Sequence_Split_Info PNTR   sequence_split;
} ID2SBioseqInfo, PNTR ID2SBioseqInfoPtr;


NLM_EXTERN ID2SBioseqInfoPtr LIBCALL ID2SBioseqInfoFree PROTO ((ID2SBioseqInfoPtr ));
NLM_EXTERN ID2SBioseqInfoPtr LIBCALL ID2SBioseqInfoNew PROTO (( void ));
NLM_EXTERN ID2SBioseqInfoPtr LIBCALL ID2SBioseqInfoAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2SBioseqInfoAsnWrite PROTO (( ID2SBioseqInfoPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ID2IdRange
*
**************************************************/
typedef struct struct_ID2_Id_Range {
   struct struct_ID2_Id_Range PNTR next;
   Int4   start;
   Int4   count;
} ID2IdRange, PNTR ID2IdRangePtr;


NLM_EXTERN ID2IdRangePtr LIBCALL ID2IdRangeFree PROTO ((ID2IdRangePtr ));
NLM_EXTERN ID2IdRangePtr LIBCALL ID2IdRangeNew PROTO (( void ));
NLM_EXTERN ID2IdRangePtr LIBCALL ID2IdRangeAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2IdRangeAsnWrite PROTO (( ID2IdRangePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ID2SSequenceSplitInfo
*
**************************************************/
typedef struct struct_ID2S_Sequence_Split_Info {
   Int4   block_size;
   Int4   chunk_start;
   ValNodePtr   chunk_blocks;
} ID2SSequenceSplitInfo, PNTR ID2SSequenceSplitInfoPtr;


NLM_EXTERN ID2SSequenceSplitInfoPtr LIBCALL ID2SSequenceSplitInfoFree PROTO ((ID2SSequenceSplitInfoPtr ));
NLM_EXTERN ID2SSequenceSplitInfoPtr LIBCALL ID2SSequenceSplitInfoNew PROTO (( void ));
NLM_EXTERN ID2SSequenceSplitInfoPtr LIBCALL ID2SSequenceSplitInfoAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2SSequenceSplitInfoAsnWrite PROTO (( ID2SSequenceSplitInfoPtr , AsnIoPtr, AsnTypePtr));

typedef ValNodePtr ID2SChunkContentPtr;
typedef ValNode ID2SChunkContent;
#define ID2SChunkContent_seq_descr 1
#define ID2SChunkContent_seq_annot 2
#define ID2SChunkContent_seq_assembly 3
#define ID2SChunkContent_seq_map 4
#define ID2SChunkContent_seq_data 5


NLM_EXTERN ID2SChunkContentPtr LIBCALL ID2SChunkContentFree PROTO ((ID2SChunkContentPtr ));
NLM_EXTERN ID2SChunkContentPtr LIBCALL ID2SChunkContentAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2SChunkContentAsnWrite PROTO (( ID2SChunkContentPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ID2SSeqDescrInfo
*
**************************************************/
typedef struct struct_ID2S_Seq_descr_Info {
   Int4   type_mask;
   struct struct_ID2_Id_Range PNTR   bioseqs;
   struct struct_ID2_Id_Range PNTR   bioseq_sets;
} ID2SSeqDescrInfo, PNTR ID2SSeqDescrInfoPtr;


NLM_EXTERN ID2SSeqDescrInfoPtr LIBCALL ID2SSeqDescrInfoFree PROTO ((ID2SSeqDescrInfoPtr ));
NLM_EXTERN ID2SSeqDescrInfoPtr LIBCALL ID2SSeqDescrInfoNew PROTO (( void ));
NLM_EXTERN ID2SSeqDescrInfoPtr LIBCALL ID2SSeqDescrInfoAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2SSeqDescrInfoAsnWrite PROTO (( ID2SSeqDescrInfoPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ID2SSeqAnnotInfo
*
**************************************************/
typedef struct struct_ID2S_Seq_annot_Info {
   CharPtr   name;
   Uint1   align;
   Uint1   graph;
   struct struct_ID2S_Feat_type_Info PNTR   feat;
   ValNodePtr   seq_loc;
} ID2SSeqAnnotInfo, PNTR ID2SSeqAnnotInfoPtr;


NLM_EXTERN ID2SSeqAnnotInfoPtr LIBCALL ID2SSeqAnnotInfoFree PROTO ((ID2SSeqAnnotInfoPtr ));
NLM_EXTERN ID2SSeqAnnotInfoPtr LIBCALL ID2SSeqAnnotInfoNew PROTO (( void ));
NLM_EXTERN ID2SSeqAnnotInfoPtr LIBCALL ID2SSeqAnnotInfoAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2SSeqAnnotInfoAsnWrite PROTO (( ID2SSeqAnnotInfoPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ID2SSeqAssemblyInfo
*
**************************************************/
typedef struct struct_ID2S_Seq_assembly_Info {
   struct struct_ID2_Id_Range PNTR   bioseqs;
} ID2SSeqAssemblyInfo, PNTR ID2SSeqAssemblyInfoPtr;


NLM_EXTERN ID2SSeqAssemblyInfoPtr LIBCALL ID2SSeqAssemblyInfoFree PROTO ((ID2SSeqAssemblyInfoPtr ));
NLM_EXTERN ID2SSeqAssemblyInfoPtr LIBCALL ID2SSeqAssemblyInfoNew PROTO (( void ));
NLM_EXTERN ID2SSeqAssemblyInfoPtr LIBCALL ID2SSeqAssemblyInfoAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2SSeqAssemblyInfoAsnWrite PROTO (( ID2SSeqAssemblyInfoPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ID2SSeqMapInfo
*
**************************************************/
#define ID2SSeqMapInfo  ValNode
#define ID2SSeqMapInfoPtr  ValNodePtr
#define ID2SSeqMapInfoFree  ValNodeFree
#define ID2SSeqMapInfoNew  ValNodeNew
#define ID2SSeqMapInfoAsnRead  ValNodeAsnRead
#define ID2SSeqMapInfoAsnWrite  ValNodeAsnWrite


/**************************************************
*
*    ID2SSeqDataInfo
*
**************************************************/
#define ID2SSeqDataInfo  ValNode
#define ID2SSeqDataInfoPtr  ValNodePtr
#define ID2SSeqDataInfoFree  ValNodeFree
#define ID2SSeqDataInfoNew  ValNodeNew
#define ID2SSeqDataInfoAsnRead  ValNodeAsnRead
#define ID2SSeqDataInfoAsnWrite  ValNodeAsnWrite


/**************************************************
*
*    ID2SFeatTypeInfo
*
**************************************************/
typedef struct struct_ID2S_Feat_type_Info {
   struct struct_ID2S_Feat_type_Info PNTR next;
   Int4   type;
   ValNodePtr   subtypes;
} ID2SFeatTypeInfo, PNTR ID2SFeatTypeInfoPtr;


NLM_EXTERN ID2SFeatTypeInfoPtr LIBCALL ID2SFeatTypeInfoFree PROTO ((ID2SFeatTypeInfoPtr ));
NLM_EXTERN ID2SFeatTypeInfoPtr LIBCALL ID2SFeatTypeInfoNew PROTO (( void ));
NLM_EXTERN ID2SFeatTypeInfoPtr LIBCALL ID2SFeatTypeInfoAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2SFeatTypeInfoAsnWrite PROTO (( ID2SFeatTypeInfoPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ID2SChunk
*
**************************************************/
typedef struct struct_ID2S_Chunk {
   struct struct_ID2S_Chunk_Data PNTR   data;
} ID2SChunk, PNTR ID2SChunkPtr;


NLM_EXTERN ID2SChunkPtr LIBCALL ID2SChunkFree PROTO ((ID2SChunkPtr ));
NLM_EXTERN ID2SChunkPtr LIBCALL ID2SChunkNew PROTO (( void ));
NLM_EXTERN ID2SChunkPtr LIBCALL ID2SChunkAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2SChunkAsnWrite PROTO (( ID2SChunkPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ID2SChunkData
*
**************************************************/
typedef struct struct_ID2S_Chunk_Data {
   struct struct_ID2S_Chunk_Data PNTR next;
   ValNodePtr   Id_id;
   ValNodePtr   descrs;
   struct struct_Seq_annot PNTR   annots;
   struct struct_Seq_align PNTR   assembly;
   struct struct_Seq_literal PNTR   seq_map;
   struct struct_Seq_literal PNTR   seq_data;
} ID2SChunkData, PNTR ID2SChunkDataPtr;


NLM_EXTERN ID2SChunkDataPtr LIBCALL ID2SChunkDataFree PROTO ((ID2SChunkDataPtr ));
NLM_EXTERN ID2SChunkDataPtr LIBCALL ID2SChunkDataNew PROTO (( void ));
NLM_EXTERN ID2SChunkDataPtr LIBCALL ID2SChunkDataAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2SChunkDataAsnWrite PROTO (( ID2SChunkDataPtr , AsnIoPtr, AsnTypePtr));


#ifdef NLM_GENERATED_CODE_PROTO

typedef ValNodePtr Id_idPtr;
typedef ValNode Id_id;

#endif /* NLM_GENERATED_CODE_PROTO */

#define Id_id_bioseq_set 1
#define Id_id_gi 2

#ifdef NLM_GENERATED_CODE_PROTO

static Id_idPtr LIBCALL Id_idFree PROTO ((Id_idPtr ));
static Id_idPtr LIBCALL Id_idAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
static Boolean LIBCALL Id_idAsnWrite PROTO (( Id_idPtr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */



/**************************************************
*
*    ID2Interval
*
**************************************************/
typedef struct struct_ID2_Interval {
   Int4   gi;
   Int4   start;
   Int4   length;
} ID2Interval, PNTR ID2IntervalPtr;


NLM_EXTERN ID2IntervalPtr LIBCALL ID2IntervalFree PROTO ((ID2IntervalPtr ));
NLM_EXTERN ID2IntervalPtr LIBCALL ID2IntervalNew PROTO (( void ));
NLM_EXTERN ID2IntervalPtr LIBCALL ID2IntervalAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2IntervalAsnWrite PROTO (( ID2IntervalPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ID2PackedSeqInts
*
**************************************************/
typedef struct struct_ID2_Packed_Seq_ints {
   Int4   gi;
   struct struct_ID2_Seq_range PNTR   ints;
} ID2PackedSeqInts, PNTR ID2PackedSeqIntsPtr;


NLM_EXTERN ID2PackedSeqIntsPtr LIBCALL ID2PackedSeqIntsFree PROTO ((ID2PackedSeqIntsPtr ));
NLM_EXTERN ID2PackedSeqIntsPtr LIBCALL ID2PackedSeqIntsNew PROTO (( void ));
NLM_EXTERN ID2PackedSeqIntsPtr LIBCALL ID2PackedSeqIntsAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2PackedSeqIntsAsnWrite PROTO (( ID2PackedSeqIntsPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ID2SeqRange
*
**************************************************/
typedef struct struct_ID2_Seq_range {
   struct struct_ID2_Seq_range PNTR next;
   Int4   start;
   Int4   length;
} ID2SeqRange, PNTR ID2SeqRangePtr;


NLM_EXTERN ID2SeqRangePtr LIBCALL ID2SeqRangeFree PROTO ((ID2SeqRangePtr ));
NLM_EXTERN ID2SeqRangePtr LIBCALL ID2SeqRangeNew PROTO (( void ));
NLM_EXTERN ID2SeqRangePtr LIBCALL ID2SeqRangeAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2SeqRangeAsnWrite PROTO (( ID2SeqRangePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ID2Param
*
**************************************************/
typedef struct struct_ID2_Param {
   struct struct_ID2_Param PNTR next;
   CharPtr   name;
   ValNodePtr   value;
   Uint2   type;
   /* following #defines are for enumerated type, not used by object loaders */
#define ID2_Param_type_set_value 1
#define ID2_Param_type_get_value 2
#define ID2_Param_type_force_value 3
#define ID2_Param_type_use_package 4

} ID2Param, PNTR ID2ParamPtr;


NLM_EXTERN ID2ParamPtr LIBCALL ID2ParamFree PROTO ((ID2ParamPtr ));
NLM_EXTERN ID2ParamPtr LIBCALL ID2ParamNew PROTO (( void ));
NLM_EXTERN ID2ParamPtr LIBCALL ID2ParamAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ID2ParamAsnWrite PROTO (( ID2ParamPtr , AsnIoPtr, AsnTypePtr));

#ifdef __cplusplus
/* { */ }
#endif

#endif /* _id2gen_ */

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

