/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "blstspc.h17";
static AsnValxNode avnx[20] = {
    {20,"blastn" ,0,0.0,&avnx[1] } ,
    {20,"blastp" ,1,0.0,&avnx[2] } ,
    {20,"blastx" ,2,0.0,&avnx[3] } ,
    {20,"tblastn" ,3,0.0,&avnx[4] } ,
    {20,"tblastx" ,4,0.0,NULL } ,
    {20,"unknown" ,0,0.0,&avnx[6] } ,
    {20,"protein" ,1,0.0,&avnx[7] } ,
    {20,"nucleotide" ,2,0.0,NULL } ,
    {20,"none" ,0,0.0,&avnx[9] } ,
    {20,"info" ,1,0.0,&avnx[10] } ,
    {20,"warn" ,2,0.0,&avnx[11] } ,
    {20,"error" ,3,0.0,&avnx[12] } ,
    {20,"fatal" ,4,0.0,NULL } ,
    {20,"notset" ,0,0.0,&avnx[14] } ,
    {20,"plus1" ,1,0.0,&avnx[15] } ,
    {20,"plus2" ,2,0.0,&avnx[16] } ,
    {20,"plus3" ,3,0.0,&avnx[17] } ,
    {20,"minus1" ,4,0.0,&avnx[18] } ,
    {20,"minus2" ,5,0.0,&avnx[19] } ,
    {20,"minus3" ,6,0.0,NULL } };

static AsnType atx[130] = {
  {401, "Blast-search" ,1,0,0,0,0,1,0,0,NULL,&atx[44],&atx[1],0,&atx[58]} ,
  {0, "program" ,128,0,0,0,0,0,0,0,NULL,&atx[2],&avnx[0],0,&atx[3]} ,
  {310, "ENUMERATED" ,0,10,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "query" ,128,1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[5]} ,
  {408, "Bioseq" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[74]} ,
  {0, "database" ,128,2,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[7]} ,
  {323, "VisibleString" ,0,26,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "parameters" ,128,3,0,1,0,0,0,0,NULL,&atx[8],NULL,0,&atx[45]} ,
  {404, "Blast-parameters" ,1,0,0,0,0,1,0,0,NULL,&atx[44],&atx[9],0,&atx[98]} ,
  {0, "first-threshold" ,128,0,0,1,0,0,0,0,NULL,&atx[10],NULL,0,&atx[11]} ,
  {302, "INTEGER" ,0,2,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "second-threshold" ,128,1,0,1,0,0,0,0,NULL,&atx[10],NULL,0,&atx[12]} ,
  {0, "cutoff" ,128,2,0,1,0,0,0,0,NULL,&atx[16],&atx[13],0,&atx[17]} ,
  {0, "evalue" ,128,0,0,0,0,0,0,0,NULL,&atx[14],NULL,0,&atx[15]} ,
  {309, "REAL" ,0,9,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "score" ,128,1,0,0,0,0,0,0,NULL,&atx[10],NULL,0,NULL} ,
  {315, "CHOICE" ,0,-1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "cutoff2" ,128,3,0,1,0,0,0,0,NULL,&atx[16],&atx[18],0,&atx[20]} ,
  {0, "evalue" ,128,0,0,0,0,0,0,0,NULL,&atx[14],NULL,0,&atx[19]} ,
  {0, "score" ,128,1,0,0,0,0,0,0,NULL,&atx[10],NULL,0,NULL} ,
  {0, "hitlist-size" ,128,4,0,1,0,0,0,0,NULL,&atx[10],NULL,0,&atx[21]} ,
  {0, "nucl-penalty" ,128,5,0,1,0,0,0,0,NULL,&atx[10],NULL,0,&atx[22]} ,
  {0, "nucl-reward" ,128,6,0,1,0,0,0,0,NULL,&atx[10],NULL,0,&atx[23]} ,
  {0, "genetic-code" ,128,7,0,1,0,0,0,0,NULL,&atx[10],NULL,0,&atx[24]} ,
  {0, "db-genetic-code" ,128,8,0,1,0,0,0,0,NULL,&atx[10],NULL,0,&atx[25]} ,
  {0, "low-complexity-filtering" ,128,9,0,1,0,0,0,0,NULL,&atx[10],NULL,0,&atx[26]} ,
  {0, "gapped-alignment" ,128,10,0,0,0,0,0,0,NULL,&atx[27],NULL,0,&atx[28]} ,
  {301, "BOOLEAN" ,0,1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "gap-open" ,128,11,0,1,0,0,0,0,NULL,&atx[10],NULL,0,&atx[29]} ,
  {0, "gap-extend" ,128,12,0,1,0,0,0,0,NULL,&atx[10],NULL,0,&atx[30]} ,
  {0, "required-start" ,128,13,0,1,0,0,0,0,NULL,&atx[10],NULL,0,&atx[31]} ,
  {0, "required-end" ,128,14,0,1,0,0,0,0,NULL,&atx[10],NULL,0,&atx[32]} ,
  {0, "ethresh" ,128,15,0,1,0,0,0,0,NULL,&atx[14],NULL,0,&atx[33]} ,
  {0, "max-num-passes" ,128,16,0,1,0,0,0,0,NULL,&atx[10],NULL,0,&atx[34]} ,
  {0, "pseudo-count-const" ,128,17,0,1,0,0,0,0,NULL,&atx[10],NULL,0,&atx[35]} ,
  {0, "other-options" ,128,18,0,1,0,0,0,0,NULL,&atx[6],NULL,0,&atx[36]} ,
  {0, "gilist" ,128,19,0,1,0,0,0,0,NULL,&atx[38],&atx[37],0,&atx[39]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[10],NULL,0,NULL} ,
  {312, "SEQUENCE OF" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "gifile" ,128,20,0,1,0,0,0,0,NULL,&atx[6],NULL,0,&atx[40]} ,
  {0, "matrix" ,128,21,0,1,0,0,0,0,NULL,&atx[6],NULL,0,&atx[41]} ,
  {0, "filter-string" ,128,22,0,1,0,0,0,0,NULL,&atx[6],NULL,0,&atx[42]} ,
  {0, "entrez-query" ,128,23,0,1,0,0,0,0,NULL,&atx[6],NULL,0,&atx[43]} ,
  {0, "word-size" ,128,24,0,1,0,0,0,0,NULL,&atx[10],NULL,0,NULL} ,
  {311, "SEQUENCE" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "mask" ,128,4,0,1,0,0,0,0,NULL,&atx[46],NULL,0,&atx[47]} ,
  {410, "Seq-loc" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[108]} ,
  {0, "matrix" ,128,5,0,1,0,0,0,0,NULL,&atx[48],NULL,0,NULL} ,
  {414, "Blast-matrix" ,1,0,0,0,0,0,0,0,NULL,&atx[44],&atx[49],0,&atx[122]} ,
  {0, "is-protein" ,128,0,0,0,0,0,0,0,NULL,&atx[27],NULL,0,&atx[50]} ,
  {0, "name" ,128,1,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[51]} ,
  {0, "comments" ,128,2,0,1,0,0,0,0,NULL,&atx[38],&atx[52],0,&atx[53]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[6],NULL,0,NULL} ,
  {0, "row-length" ,128,3,0,0,0,0,0,0,NULL,&atx[10],NULL,0,&atx[54]} ,
  {0, "column-length" ,128,4,0,0,0,0,0,0,NULL,&atx[10],NULL,0,&atx[55]} ,
  {0, "scores" ,128,5,0,0,0,0,0,0,NULL,&atx[38],&atx[56],0,&atx[57]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[10],NULL,0,NULL} ,
  {0, "karlinK" ,128,6,0,0,0,0,0,0,NULL,&atx[14],NULL,0,NULL} ,
  {402, "Blast-request" ,1,0,0,0,0,1,0,0,NULL,&atx[16],&atx[59],0,&atx[78]} ,
  {0, "init" ,128,0,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[60]} ,
  {0, "motd" ,128,1,0,0,0,0,0,0,NULL,&atx[61],NULL,0,&atx[62]} ,
  {305, "NULL" ,0,5,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "db-info" ,128,2,0,0,0,0,0,0,NULL,&atx[61],NULL,0,&atx[63]} ,
  {0, "db-info-specific" ,128,3,0,0,0,0,0,0,NULL,&atx[64],NULL,0,&atx[67]} ,
  {412, "Blast-dbinfo-get" ,1,0,0,0,0,0,0,0,NULL,&atx[44],&atx[65],0,&atx[70]} ,
  {0, "name" ,128,0,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[66]} ,
  {0, "type" ,128,1,0,0,0,0,0,0,NULL,&atx[2],&avnx[5],0,NULL} ,
  {0, "matrix-get" ,128,4,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[68]} ,
  {0, "search" ,128,5,0,0,0,0,0,0,NULL,&atx[0],NULL,0,&atx[69]} ,
  {0, "db-seq-get" ,128,6,0,0,0,0,0,0,NULL,&atx[70],NULL,0,&atx[75]} ,
  {413, "Blast-seq-id" ,1,0,0,0,0,0,0,0,NULL,&atx[44],&atx[71],0,&atx[48]} ,
  {0, "is-protein" ,128,0,0,0,0,0,0,0,NULL,&atx[27],NULL,0,&atx[72]} ,
  {0, "database" ,128,1,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[73]} ,
  {0, "id" ,128,2,0,0,0,0,0,0,NULL,&atx[74],NULL,0,NULL} ,
  {409, "Seq-id" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[46]} ,
  {0, "db-redundant-ids-get" ,128,7,0,0,0,0,0,0,NULL,&atx[70],NULL,0,&atx[76]} ,
  {0, "db-redundant-descr-get" ,128,8,0,0,0,0,0,0,NULL,&atx[70],NULL,0,&atx[77]} ,
  {0, "fini" ,128,9,0,0,0,0,0,0,NULL,&atx[61],NULL,0,NULL} ,
  {403, "Blast-response" ,1,0,0,0,0,1,0,0,NULL,&atx[16],&atx[79],0,&atx[8]} ,
  {0, "init" ,128,0,0,0,0,0,0,0,NULL,&atx[80],NULL,0,&atx[83]} ,
  {418, "Blast-version" ,1,0,0,0,0,0,0,0,NULL,&atx[44],&atx[81],0,&atx[85]} ,
  {0, "version" ,128,0,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[82]} ,
  {0, "date" ,128,1,0,0,0,0,0,0,NULL,&atx[6],NULL,0,NULL} ,
  {0, "motd" ,128,1,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[84]} ,
  {0, "error" ,128,2,0,0,0,0,0,0,NULL,&atx[85],NULL,0,&atx[88]} ,
  {419, "Blast-error" ,1,0,0,0,0,0,0,0,NULL,&atx[44],&atx[86],0,NULL} ,
  {0, "level" ,128,0,0,0,0,0,0,0,NULL,&atx[2],&avnx[8],0,&atx[87]} ,
  {0, "msg" ,128,1,0,1,0,0,0,0,NULL,&atx[6],NULL,0,NULL} ,
  {0, "db-seq-get" ,128,3,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[89]} ,
  {0, "db-redundant-ids-get" ,128,4,0,0,0,0,0,0,NULL,&atx[38],&atx[90],0,&atx[91]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[74],NULL,0,NULL} ,
  {0, "db-redundant-descr-get" ,128,5,0,0,0,0,0,0,NULL,&atx[38],&atx[92],0,&atx[96]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[93],NULL,0,NULL} ,
  {417, "Blast-defline" ,1,0,0,0,0,0,0,0,NULL,&atx[44],&atx[94],0,&atx[80]} ,
  {0, "id" ,128,0,0,0,0,0,0,0,NULL,&atx[74],NULL,0,&atx[95]} ,
  {0, "defline" ,128,1,0,0,0,0,0,0,NULL,&atx[6],NULL,0,NULL} ,
  {0, "db-info" ,128,6,0,0,0,0,0,0,NULL,&atx[38],&atx[97],0,&atx[105]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[98],NULL,0,NULL} ,
  {405, "Blast-dbinfo" ,1,0,0,0,0,1,0,0,NULL,&atx[44],&atx[99],0,&atx[110]} ,
  {0, "is-protein" ,128,0,0,0,0,0,0,0,NULL,&atx[27],NULL,0,&atx[100]} ,
  {0, "name" ,128,1,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[101]} ,
  {0, "definition" ,128,2,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[102]} ,
  {0, "date" ,128,3,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[103]} ,
  {0, "total-length" ,128,4,0,0,0,0,0,0,NULL,&atx[10],NULL,0,&atx[104]} ,
  {0, "number-seqs" ,128,5,0,0,0,0,0,0,NULL,&atx[10],NULL,0,NULL} ,
  {0, "db-info-specific" ,128,7,0,0,0,0,0,0,NULL,&atx[98],NULL,0,&atx[106]} ,
  {0, "matrix" ,128,8,0,0,0,0,0,0,NULL,&atx[48],NULL,0,&atx[107]} ,
  {0, "alignment" ,128,9,0,0,0,0,0,0,NULL,&atx[108],NULL,0,&atx[109]} ,
  {411, "Seq-align-set" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[64]} ,
  {0, "mask" ,128,10,0,0,0,0,0,0,NULL,&atx[110],NULL,0,&atx[114]} ,
  {406, "Blast-mask" ,1,0,0,0,0,1,0,0,NULL,&atx[44],&atx[111],0,&atx[115]} ,
  {0, "location" ,128,0,0,0,0,0,0,0,NULL,&atx[38],&atx[112],0,&atx[113]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[46],NULL,0,NULL} ,
  {0, "frame" ,128,1,0,0,0,0,0,0,NULL,&atx[2],&avnx[13],0,NULL} ,
  {0, "kablk" ,128,11,0,0,0,0,0,0,NULL,&atx[115],NULL,0,&atx[120]} ,
  {407, "Blast-KABlk" ,1,0,0,0,0,1,0,0,NULL,&atx[44],&atx[116],0,&atx[4]} ,
  {0, "lambda" ,128,0,0,0,0,0,0,0,NULL,&atx[14],NULL,0,&atx[117]} ,
  {0, "k" ,128,1,0,0,0,0,0,0,NULL,&atx[14],NULL,0,&atx[118]} ,
  {0, "h" ,128,2,0,0,0,0,0,0,NULL,&atx[14],NULL,0,&atx[119]} ,
  {0, "gapped" ,128,3,0,0,0,0,0,0,NULL,&atx[27],NULL,0,NULL} ,
  {0, "parameters" ,128,12,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[121]} ,
  {0, "queued" ,128,13,0,0,0,0,0,0,NULL,&atx[122],NULL,0,&atx[124]} ,
  {415, "Blast-Queued" ,1,0,0,0,0,0,0,0,NULL,&atx[44],&atx[123],0,&atx[125]} ,
  {0, "length" ,128,0,0,0,0,0,0,0,NULL,&atx[10],NULL,0,NULL} ,
  {0, "start" ,128,14,0,0,0,0,0,0,NULL,&atx[125],NULL,0,&atx[127]} ,
  {416, "Blast-Progress" ,1,0,0,0,0,0,0,0,NULL,&atx[44],&atx[126],0,&atx[93]} ,
  {0, "completed" ,128,0,0,0,0,0,0,0,NULL,&atx[10],NULL,0,NULL} ,
  {0, "progress" ,128,15,0,0,0,0,0,0,NULL,&atx[125],NULL,0,&atx[128]} ,
  {0, "done" ,128,16,0,0,0,0,0,0,NULL,&atx[125],NULL,0,&atx[129]} ,
  {0, "fini" ,128,17,0,0,0,0,0,0,NULL,&atx[61],NULL,0,NULL} };

static AsnModule ampx[1] = {
  { "NCBI-Blast" , "blstspc.h17",&atx[0],NULL,NULL,0,0} };

static AsnValxNodePtr avn = avnx;
static AsnTypePtr at = atx;
static AsnModulePtr amp = ampx;



/**************************************************
*
*    Defines for Module NCBI-Blast
*
**************************************************/

#define BLAST_SEARCH &at[0]
#define BLAST_SEARCH_program &at[1]
#define BLAST_SEARCH_query &at[3]
#define BLAST_SEARCH_database &at[5]
#define BLAST_SEARCH_parameters &at[7]
#define BLAST_SEARCH_mask &at[45]
#define BLAST_SEARCH_matrix &at[47]

#define BLAST_REQUEST &at[58]
#define BLAST_REQUEST_init &at[59]
#define BLAST_REQUEST_motd &at[60]
#define BLAST_REQUEST_db_info &at[62]
#define BLAST_REQUEST_db_info_specific &at[63]
#define BLAST_REQUEST_matrix_get &at[67]
#define BLAST_REQUEST_search &at[68]
#define BLAST_REQUEST_db_seq_get &at[69]
#define BLAST_REQUEST_db_redundant_ids_get &at[75]
#define BLAST_REQUEST_db_redundant_descr_get &at[76]
#define BLAST_REQUEST_fini &at[77]

#define BLAST_RESPONSE &at[78]
#define BLAST_RESPONSE_init &at[79]
#define BLAST_RESPONSE_motd &at[83]
#define BLAST_RESPONSE_error &at[84]
#define BLAST_RESPONSE_db_seq_get &at[88]
#define BLAST_RESPONSE_db_redundant_ids_get &at[89]
#define BLAST_RESPONSE_db_redundant_ids_get_E &at[90]
#define BLAST_RESPONSE_db_redundant_descr_get &at[91]
#define BLAST_RESPONSE_db_redundant_descr_get_E &at[92]
#define BLAST_RESPONSE_db_info &at[96]
#define BLAST_RESPONSE_db_info_E &at[97]
#define BLAST_RESPONSE_db_info_specific &at[105]
#define BLAST_RESPONSE_matrix &at[106]
#define BLAST_RESPONSE_alignment &at[107]
#define BLAST_RESPONSE_mask &at[109]
#define BLAST_RESPONSE_kablk &at[114]
#define BLAST_RESPONSE_parameters &at[120]
#define BLAST_RESPONSE_queued &at[121]
#define BLAST_RESPONSE_start &at[124]
#define BLAST_RESPONSE_progress &at[127]
#define BLAST_RESPONSE_done &at[128]
#define BLAST_RESPONSE_fini &at[129]

#define BLAST_PARAMETERS &at[8]
#define BLAST_PARAMETERS_first_threshold &at[9]
#define BLAST_PARAMETERS_second_threshold &at[11]
#define BLAST_PARAMETERS_cutoff &at[12]
#define BLAST_PARAMETERS_cutoff_evalue &at[13]
#define BLAST_PARAMETERS_cutoff_score &at[15]
#define BLAST_PARAMETERS_cutoff2 &at[17]
#define BLAST_PARAMETERS_cutoff2_evalue &at[18]
#define BLAST_PARAMETERS_cutoff2_score &at[19]
#define BLAST_PARAMETERS_hitlist_size &at[20]
#define BLAST_PARAMETERS_nucl_penalty &at[21]
#define BLAST_PARAMETERS_nucl_reward &at[22]
#define BLAST_PARAMETERS_genetic_code &at[23]
#define BLAST_PARAMETERS_db_genetic_code &at[24]
#define BLAST_PARAMETERS_low_complexity_filtering &at[25]
#define BLAST_PARAMETERS_gapped_alignment &at[26]
#define BLAST_PARAMETERS_gap_open &at[28]
#define BLAST_PARAMETERS_gap_extend &at[29]
#define BLAST_PARAMETERS_required_start &at[30]
#define BLAST_PARAMETERS_required_end &at[31]
#define BLAST_PARAMETERS_ethresh &at[32]
#define BLAST_PARAMETERS_max_num_passes &at[33]
#define BLAST_PARAMETERS_pseudo_count_const &at[34]
#define BLAST_PARAMETERS_other_options &at[35]
#define BLAST_PARAMETERS_gilist &at[36]
#define BLAST_PARAMETERS_gilist_E &at[37]
#define BLAST_PARAMETERS_gifile &at[39]
#define BLAST_PARAMETERS_matrix &at[40]
#define BLAST_PARAMETERS_filter_string &at[41]
#define BLAST_PARAMETERS_entrez_query &at[42]
#define BLAST_PARAMETERS_word_size &at[43]

#define BLAST_DBINFO &at[98]
#define BLAST_DBINFO_is_protein &at[99]
#define BLAST_DBINFO_name &at[100]
#define BLAST_DBINFO_definition &at[101]
#define BLAST_DBINFO_date &at[102]
#define BLAST_DBINFO_total_length &at[103]
#define BLAST_DBINFO_number_seqs &at[104]

#define BLAST_MASK &at[110]
#define BLAST_MASK_location &at[111]
#define BLAST_MASK_location_E &at[112]
#define BLAST_MASK_frame &at[113]

#define BLAST_KABLK &at[115]
#define BLAST_KABLK_lambda &at[116]
#define BLAST_KABLK_k &at[117]
#define BLAST_KABLK_h &at[118]
#define BLAST_KABLK_gapped &at[119]

#define BLAST_DBINFO_GET &at[64]
#define BLAST_DBINFO_GET_name &at[65]
#define BLAST_DBINFO_GET_type &at[66]

#define BLAST_SEQ_ID &at[70]
#define BLAST_SEQ_ID_is_protein &at[71]
#define BLAST_SEQ_ID_database &at[72]
#define BLAST_SEQ_ID_id &at[73]

#define BLAST_MATRIX &at[48]
#define BLAST_MATRIX_is_protein &at[49]
#define BLAST_MATRIX_name &at[50]
#define BLAST_MATRIX_comments &at[51]
#define BLAST_MATRIX_comments_E &at[52]
#define BLAST_MATRIX_row_length &at[53]
#define BLAST_MATRIX_column_length &at[54]
#define BLAST_MATRIX_scores &at[55]
#define BLAST_MATRIX_scores_E &at[56]
#define BLAST_MATRIX_karlinK &at[57]

#define BLAST_QUEUED &at[122]
#define BLAST_QUEUED_length &at[123]

#define BLAST_PROGRESS &at[125]
#define BLAST_PROGRESS_completed &at[126]

#define BLAST_DEFLINE &at[93]
#define BLAST_DEFLINE_id &at[94]
#define BLAST_DEFLINE_defline &at[95]

#define BLAST_VERSION &at[80]
#define BLAST_VERSION_version &at[81]
#define BLAST_VERSION_date &at[82]

#define BLAST_ERROR &at[85]
#define BLAST_ERROR_level &at[86]
#define BLAST_ERROR_msg &at[87]