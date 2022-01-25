/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "blstspc.l17";
static AsnValxNodePtr avn;
static AsnTypePtr at;
static AsnModulePtr amp;


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
