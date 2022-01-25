/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "seqsplit.h15";
static AsnValxNode avnx[3] = {
    {3,NULL,1,0.0,NULL } ,
    {3,NULL,1,0.0,NULL } ,
    {3,NULL,1,0.0,NULL } };

static AsnType atx[123] = {
  {401, "ID2S-Chunk-Id" ,1,0,0,0,0,1,0,0,NULL,&atx[1],NULL,0,&atx[2]} ,
  {302, "INTEGER" ,0,2,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {402, "ID2S-Seq-annot-Info" ,1,0,0,0,0,1,0,0,NULL,&atx[15],&atx[3],0,&atx[41]} ,
  {0, "name" ,128,0,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[5]} ,
  {323, "VisibleString" ,0,26,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "align" ,128,1,0,1,0,0,0,0,NULL,&atx[6],NULL,0,&atx[7]} ,
  {305, "NULL" ,0,5,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "graph" ,128,2,0,1,0,0,0,0,NULL,&atx[6],NULL,0,&atx[8]} ,
  {0, "feat" ,128,3,0,1,0,0,0,0,NULL,&atx[14],&atx[9],0,&atx[16]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[10],NULL,0,NULL} ,
  {425, "ID2S-Feat-type-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[11],0,&atx[17]} ,
  {0, "type" ,128,0,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[12]} ,
  {0, "subtypes" ,128,1,0,1,0,0,0,0,NULL,&atx[14],&atx[13],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[1],NULL,0,NULL} ,
  {314, "SET OF" ,0,17,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {311, "SEQUENCE" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "seq-loc" ,128,4,0,1,0,0,0,0,NULL,&atx[17],NULL,0,NULL} ,
  {426, "ID2-Seq-loc" ,1,0,0,0,0,0,0,0,NULL,&atx[40],&atx[18],0,&atx[99]} ,
  {0, "gi-whole" ,128,0,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[19]} ,
  {0, "interval" ,128,1,0,0,0,0,0,0,NULL,&atx[20],NULL,0,&atx[24]} ,
  {430, "ID2-Interval" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[21],0,&atx[25]} ,
  {0, "gi" ,128,0,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[22]} ,
  {0, "start" ,128,1,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[23]} ,
  {0, "length" ,128,2,0,0,1,0,0,0,&avnx[0],&atx[1],NULL,0,NULL} ,
  {0, "packed-ints" ,128,2,0,0,0,0,0,0,NULL,&atx[25],NULL,0,&atx[32]} ,
  {431, "ID2-Packed-Seq-ints" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[26],0,&atx[29]} ,
  {0, "gi" ,128,0,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[27]} ,
  {0, "intervals" ,128,1,0,0,0,0,0,0,NULL,&atx[14],&atx[28],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[29],NULL,0,NULL} ,
  {432, "ID2-Seq-range" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[30],0,NULL} ,
  {0, "start" ,128,0,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[31]} ,
  {0, "length" ,128,1,0,0,1,0,0,0,&avnx[1],&atx[1],NULL,0,NULL} ,
  {0, "gi-whole-range" ,128,3,0,0,0,0,0,0,NULL,&atx[33],NULL,0,&atx[36]} ,
  {414, "ID2-Id-Range" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[34],0,&atx[57]} ,
  {0, "start" ,128,0,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[35]} ,
  {0, "count" ,128,1,0,0,1,0,0,0,&avnx[2],&atx[1],NULL,0,NULL} ,
  {0, "loc-set" ,128,4,0,0,0,0,0,0,NULL,&atx[14],&atx[37],0,&atx[38]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[17],NULL,0,NULL} ,
  {0, "seq-loc" ,128,5,0,0,0,0,0,0,NULL,&atx[39],NULL,0,NULL} ,
  {404, "Seq-loc" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[42]} ,
  {315, "CHOICE" ,0,-1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {403, "Seq-id" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[39]} ,
  {405, "Bioseq" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[43]} ,
  {406, "Seq-annot" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[44]} ,
  {407, "Seq-descr" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[45]} ,
  {408, "Seq-literal" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[46]} ,
  {409, "Seq-align" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[47]} ,
  {410, "ID2S-Split-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[48],0,&atx[50]} ,
  {0, "bioseqs-info" ,128,0,0,1,0,0,0,0,NULL,&atx[14],&atx[49],0,&atx[64]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[50],NULL,0,NULL} ,
  {411, "ID2S-Bioseqs-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[51],0,&atx[66]} ,
  {0, "info" ,128,0,0,0,0,0,0,0,NULL,&atx[52],NULL,0,&atx[63]} ,
  {413, "ID2S-Bioseq-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[53],0,&atx[33]} ,
  {0, "gap-count" ,128,0,0,1,0,0,0,0,NULL,&atx[1],NULL,0,&atx[54]} ,
  {0, "seq-map-has-ref" ,128,1,0,1,0,0,0,0,NULL,&atx[55],NULL,0,&atx[56]} ,
  {301, "BOOLEAN" ,0,1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "sequence-split" ,128,2,0,1,0,0,0,0,NULL,&atx[57],NULL,0,NULL} ,
  {415, "ID2S-Sequence-Split-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[58],0,&atx[70]} ,
  {0, "block-size" ,128,0,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[59]} ,
  {0, "chunk-start" ,128,1,0,0,0,0,0,0,NULL,&atx[0],NULL,0,&atx[60]} ,
  {0, "chunk-blocks" ,128,2,0,0,0,0,0,0,NULL,&atx[62],&atx[61],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[1],NULL,0,NULL} ,
  {312, "SEQUENCE OF" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "bioseqs" ,128,1,0,0,0,0,0,0,NULL,&atx[33],NULL,0,NULL} ,
  {0, "chunks" ,128,1,0,0,0,0,0,0,NULL,&atx[14],&atx[65],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[66],NULL,0,NULL} ,
  {412, "ID2S-Chunk-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[67],0,&atx[52]} ,
  {0, "id" ,128,0,0,0,0,0,0,0,NULL,&atx[0],NULL,0,&atx[68]} ,
  {0, "content" ,128,1,0,0,0,0,0,0,NULL,&atx[14],&atx[69],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[70],NULL,0,NULL} ,
  {416, "ID2S-Chunk-Content" ,1,0,0,0,0,0,0,0,NULL,&atx[40],&atx[71],0,&atx[72]} ,
  {0, "seq-descr" ,128,0,0,0,0,0,0,0,NULL,&atx[72],NULL,0,&atx[80]} ,
  {417, "ID2S-Seq-descr-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[73],0,&atx[82]} ,
  {0, "type-mask" ,128,0,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[74]} ,
  {0, "bioseqs" ,128,1,0,1,0,0,0,0,NULL,&atx[75],NULL,0,&atx[77]} ,
  {423, "ID2-Bioseq-Ids" ,1,0,0,0,0,0,0,0,NULL,&atx[14],&atx[76],0,&atx[78]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[33],NULL,0,NULL} ,
  {0, "bioseq-sets" ,128,2,0,1,0,0,0,0,NULL,&atx[78],NULL,0,NULL} ,
  {424, "ID2-Bioseq-set-Ids" ,1,0,0,0,0,0,0,0,NULL,&atx[14],&atx[79],0,&atx[10]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[1],NULL,0,NULL} ,
  {0, "seq-annot" ,128,1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[81]} ,
  {0, "seq-assembly" ,128,2,0,0,0,0,0,0,NULL,&atx[82],NULL,0,&atx[84]} ,
  {418, "ID2S-Seq-assembly-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[83],0,&atx[85]} ,
  {0, "bioseqs" ,128,0,0,0,0,0,0,0,NULL,&atx[75],NULL,0,NULL} ,
  {0, "seq-map" ,128,3,0,0,0,0,0,0,NULL,&atx[85],NULL,0,&atx[86]} ,
  {419, "ID2S-Seq-map-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[17],NULL,0,&atx[87]} ,
  {0, "seq-data" ,128,4,0,0,0,0,0,0,NULL,&atx[87],NULL,0,&atx[88]} ,
  {420, "ID2S-Seq-data-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[17],NULL,0,&atx[89]} ,
  {0, "seq-annot-place" ,128,5,0,0,0,0,0,0,NULL,&atx[89],NULL,0,&atx[93]} ,
  {421, "ID2S-Seq-annot-place-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[90],0,&atx[95]} ,
  {0, "name" ,128,0,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[91]} ,
  {0, "bioseqs" ,128,1,0,1,0,0,0,0,NULL,&atx[75],NULL,0,&atx[92]} ,
  {0, "bioseq-sets" ,128,2,0,1,0,0,0,0,NULL,&atx[78],NULL,0,NULL} ,
  {0, "bioseq-place" ,128,6,0,0,0,0,0,0,NULL,&atx[14],&atx[94],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[95],NULL,0,NULL} ,
  {422, "ID2S-Bioseq-place-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[96],0,&atx[75]} ,
  {0, "bioseq-set" ,128,0,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[97]} ,
  {0, "seq-ids" ,128,1,0,0,0,0,0,0,NULL,&atx[14],&atx[98],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[41],NULL,0,NULL} ,
  {427, "ID2S-Chunk" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[100],0,&atx[102]} ,
  {0, "data" ,128,0,0,0,0,0,0,0,NULL,&atx[14],&atx[101],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[102],NULL,0,NULL} ,
  {428, "ID2S-Chunk-Data" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[103],0,&atx[115]} ,
  {0, "id" ,128,0,0,0,0,0,0,0,NULL,&atx[40],&atx[104],0,&atx[107]} ,
  {0, "bioseq-set" ,128,0,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[105]} ,
  {0, "gi" ,128,1,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[106]} ,
  {0, "seq-id" ,128,2,0,0,0,0,0,0,NULL,&atx[41],NULL,0,NULL} ,
  {0, "descrs" ,128,1,0,1,0,0,0,0,NULL,&atx[14],&atx[108],0,&atx[109]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[44],NULL,0,NULL} ,
  {0, "annots" ,128,2,0,1,0,0,0,0,NULL,&atx[14],&atx[110],0,&atx[111]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[43],NULL,0,NULL} ,
  {0, "assembly" ,128,3,0,1,0,0,0,0,NULL,&atx[14],&atx[112],0,&atx[113]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[46],NULL,0,NULL} ,
  {0, "seq-map" ,128,4,0,1,0,0,0,0,NULL,&atx[62],&atx[114],0,&atx[119]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[115],NULL,0,NULL} ,
  {429, "ID2S-Sequence-Piece" ,1,0,0,0,0,0,0,0,NULL,&atx[15],&atx[116],0,&atx[20]} ,
  {0, "start" ,128,0,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[117]} ,
  {0, "data" ,128,1,0,0,0,0,0,0,NULL,&atx[62],&atx[118],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[45],NULL,0,NULL} ,
  {0, "seq-data" ,128,5,0,1,0,0,0,0,NULL,&atx[62],&atx[120],0,&atx[121]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[115],NULL,0,NULL} ,
  {0, "bioseqs" ,128,6,0,1,0,0,0,0,NULL,&atx[14],&atx[122],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[42],NULL,0,NULL} };

static AsnModule ampx[1] = {
  { "NCBI-Seq-split" , "seqsplit.h15",&atx[0],NULL,NULL,0,0} };

static AsnValxNodePtr avn = avnx;
static AsnTypePtr at = atx;
static AsnModulePtr amp = ampx;



/**************************************************
*
*    Defines for Module NCBI-Seq-split
*
**************************************************/

#define ID2S_CHUNK_ID &at[0]

#define ID2S_SEQ_ANNOT_INFO &at[2]
#define ID2S_SEQ_ANNOT_INFO_name &at[3]
#define ID2S_SEQ_ANNOT_INFO_align &at[5]
#define ID2S_SEQ_ANNOT_INFO_graph &at[7]
#define ID2S_SEQ_ANNOT_INFO_feat &at[8]
#define ID2S_SEQ_ANNOT_INFO_feat_E &at[9]
#define ID2S_SEQ_ANNOT_INFO_seq_loc &at[16]

#define ID2S_SPLIT_INFO &at[47]
#define ID2S_SPLIT_INFO_bioseqs_info &at[48]
#define ID2S_SPLIT_INFO_bioseqs_info_E &at[49]
#define ID2S_SPLIT_INFO_chunks &at[64]
#define ID2S_SPLIT_INFO_chunks_E &at[65]

#define ID2S_BIOSEQS_INFO &at[50]
#define ID2S_BIOSEQS_INFO_info &at[51]
#define ID2S_BIOSEQS_INFO_bioseqs &at[63]

#define ID2S_CHUNK_INFO &at[66]
#define ID2S_CHUNK_INFO_id &at[67]
#define ID2S_CHUNK_INFO_content &at[68]
#define ID2S_CHUNK_INFO_content_E &at[69]

#define ID2S_BIOSEQ_INFO &at[52]
#define ID2S_BIOSEQ_INFO_gap_count &at[53]
#define BIOSEQ_INFO_seq_map_has_ref &at[54]
#define ID2S_BIOSEQ_INFO_sequence_split &at[56]

#define ID2_ID_RANGE &at[33]
#define ID2_ID_RANGE_start &at[34]
#define ID2_ID_RANGE_count &at[35]

#define ID2S_SEQUENCE_SPLIT_INFO &at[57]
#define SEQUENCE_SPLIT_INFO_block_size &at[58]
#define SPLIT_INFO_chunk_start &at[59]
#define SPLIT_INFO_chunk_blocks &at[60]
#define SPLIT_INFO_chunk_blocks_E &at[61]

#define ID2S_CHUNK_CONTENT &at[70]
#define ID2S_CHUNK_CONTENT_seq_descr &at[71]
#define ID2S_CHUNK_CONTENT_seq_annot &at[80]
#define ID2S_CHUNK_CONTENT_seq_assembly &at[81]
#define ID2S_CHUNK_CONTENT_seq_map &at[84]
#define ID2S_CHUNK_CONTENT_seq_data &at[86]
#define CHUNK_CONTENT_seq_annot_place &at[88]
#define ID2S_CHUNK_CONTENT_bioseq_place &at[93]
#define CHUNK_CONTENT_bioseq_place_E &at[94]

#define ID2S_SEQ_DESCR_INFO &at[72]
#define ID2S_SEQ_DESCR_INFO_type_mask &at[73]
#define ID2S_SEQ_DESCR_INFO_bioseqs &at[74]
#define ID2S_SEQ_DESCR_INFO_bioseq_sets &at[77]

#define ID2S_SEQ_ASSEMBLY_INFO &at[82]
#define ID2S_SEQ_ASSEMBLY_INFO_bioseqs &at[83]

#define ID2S_SEQ_MAP_INFO &at[85]

#define ID2S_SEQ_DATA_INFO &at[87]

#define ID2S_SEQ_ANNOT_PLACE_INFO &at[89]
#define ID2S_SEQ_ANNOT_PLACE_INFO_name &at[90]
#define SEQ_ANNOT_PLACE_INFO_bioseqs &at[91]
#define ANNOT_PLACE_INFO_bioseq_sets &at[92]

#define ID2S_BIOSEQ_PLACE_INFO &at[95]
#define BIOSEQ_PLACE_INFO_bioseq_set &at[96]
#define ID2S_BIOSEQ_PLACE_INFO_seq_ids &at[97]
#define BIOSEQ_PLACE_INFO_seq_ids_E &at[98]

#define ID2_BIOSEQ_IDS &at[75]
#define ID2_BIOSEQ_IDS_E &at[76]

#define ID2_BIOSEQ_SET_IDS &at[78]
#define ID2_BIOSEQ_SET_IDS_E &at[79]

#define ID2S_FEAT_TYPE_INFO &at[10]
#define ID2S_FEAT_TYPE_INFO_type &at[11]
#define ID2S_FEAT_TYPE_INFO_subtypes &at[12]
#define ID2S_FEAT_TYPE_INFO_subtypes_E &at[13]

#define ID2_SEQ_LOC &at[17]
#define ID2_SEQ_LOC_gi_whole &at[18]
#define ID2_SEQ_LOC_interval &at[19]
#define ID2_SEQ_LOC_packed_ints &at[24]
#define ID2_SEQ_LOC_gi_whole_range &at[32]
#define ID2_SEQ_LOC_loc_set &at[36]
#define ID2_SEQ_LOC_loc_set_E &at[37]
#define ID2_SEQ_LOC_seq_loc &at[38]

#define ID2S_CHUNK &at[99]
#define ID2S_CHUNK_data &at[100]
#define ID2S_CHUNK_data_E &at[101]

#define ID2S_CHUNK_DATA &at[102]
#define ID2S_CHUNK_DATA_id &at[103]
#define ID2S_CHUNK_DATA_id_bioseq_set &at[104]
#define ID2S_CHUNK_DATA_id_gi &at[105]
#define ID2S_CHUNK_DATA_id_seq_id &at[106]
#define ID2S_CHUNK_DATA_descrs &at[107]
#define ID2S_CHUNK_DATA_descrs_E &at[108]
#define ID2S_CHUNK_DATA_annots &at[109]
#define ID2S_CHUNK_DATA_annots_E &at[110]
#define ID2S_CHUNK_DATA_assembly &at[111]
#define ID2S_CHUNK_DATA_assembly_E &at[112]
#define ID2S_CHUNK_DATA_seq_map &at[113]
#define ID2S_CHUNK_DATA_seq_map_E &at[114]
#define ID2S_CHUNK_DATA_seq_data &at[119]
#define ID2S_CHUNK_DATA_seq_data_E &at[120]
#define ID2S_CHUNK_DATA_bioseqs &at[121]
#define ID2S_CHUNK_DATA_bioseqs_E &at[122]

#define ID2S_SEQUENCE_PIECE &at[115]
#define ID2S_SEQUENCE_PIECE_start &at[116]
#define ID2S_SEQUENCE_PIECE_data &at[117]
#define ID2S_SEQUENCE_PIECE_data_E &at[118]

#define ID2_INTERVAL &at[20]
#define ID2_INTERVAL_gi &at[21]
#define ID2_INTERVAL_start &at[22]
#define ID2_INTERVAL_length &at[23]

#define ID2_PACKED_SEQ_INTS &at[25]
#define ID2_PACKED_SEQ_INTS_gi &at[26]
#define ID2_PACKED_SEQ_INTS_intervals &at[27]
#define ID2_PACKED_SEQ_INTS_intervals_E &at[28]

#define ID2_SEQ_RANGE &at[29]
#define ID2_SEQ_RANGE_start &at[30]
#define ID2_SEQ_RANGE_length &at[31]
