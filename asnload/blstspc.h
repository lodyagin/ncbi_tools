/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "blstspc.l18";
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
#define BLAST_SEARCH_mask &at[51]
#define BLAST_SEARCH_matrix &at[53]

#define BLAST_REQUEST &at[64]
#define BLAST_REQUEST_init &at[65]
#define BLAST_REQUEST_motd &at[66]
#define BLAST_REQUEST_db_info &at[68]
#define BLAST_REQUEST_db_info_specific &at[69]
#define BLAST_REQUEST_matrix_get &at[73]
#define BLAST_REQUEST_search &at[74]
#define BLAST_REQUEST_db_seq_get &at[75]
#define BLAST_REQUEST_db_redundant_ids_get &at[81]
#define BLAST_REQUEST_db_redundant_descr_get &at[82]
#define BLAST_REQUEST_fini &at[83]

#define BLAST_RESPONSE &at[84]
#define BLAST_RESPONSE_init &at[85]
#define BLAST_RESPONSE_motd &at[89]
#define BLAST_RESPONSE_error &at[90]
#define BLAST_RESPONSE_db_seq_get &at[94]
#define BLAST_RESPONSE_db_redundant_ids_get &at[95]
#define BLAST_RESPONSE_db_redundant_ids_get_E &at[96]
#define BLAST_RESPONSE_db_redundant_descr_get &at[97]
#define BLAST_RESPONSE_db_redundant_descr_get_E &at[98]
#define BLAST_RESPONSE_db_info &at[102]
#define BLAST_RESPONSE_db_info_E &at[103]
#define BLAST_RESPONSE_db_info_specific &at[111]
#define BLAST_RESPONSE_matrix &at[112]
#define BLAST_RESPONSE_alignment &at[113]
#define BLAST_RESPONSE_mask &at[115]
#define BLAST_RESPONSE_kablk &at[120]
#define BLAST_RESPONSE_parameters &at[126]
#define BLAST_RESPONSE_queued &at[127]
#define BLAST_RESPONSE_start &at[130]
#define BLAST_RESPONSE_progress &at[133]
#define BLAST_RESPONSE_done &at[134]
#define BLAST_RESPONSE_fini &at[135]

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
#define BLAST_PARAMETERS_db_length &at[44]
#define BLAST_PARAMETERS_searchsp_eff &at[45]
#define BLAST_PARAMETERS_hsp_range_max &at[46]
#define BLAST_PARAMETERS_block_width &at[47]
#define BLAST_PARAMETERS_perform_culling &at[48]
#define BLAST_PARAMETERS_strand_option &at[49]

#define BLAST_DBINFO &at[104]
#define BLAST_DBINFO_is_protein &at[105]
#define BLAST_DBINFO_name &at[106]
#define BLAST_DBINFO_definition &at[107]
#define BLAST_DBINFO_date &at[108]
#define BLAST_DBINFO_total_length &at[109]
#define BLAST_DBINFO_number_seqs &at[110]

#define BLAST_MASK &at[116]
#define BLAST_MASK_location &at[117]
#define BLAST_MASK_location_E &at[118]
#define BLAST_MASK_frame &at[119]

#define BLAST_KABLK &at[121]
#define BLAST_KABLK_lambda &at[122]
#define BLAST_KABLK_k &at[123]
#define BLAST_KABLK_h &at[124]
#define BLAST_KABLK_gapped &at[125]

#define BLAST_DBINFO_GET &at[70]
#define BLAST_DBINFO_GET_name &at[71]
#define BLAST_DBINFO_GET_type &at[72]

#define BLAST_SEQ_ID &at[76]
#define BLAST_SEQ_ID_is_protein &at[77]
#define BLAST_SEQ_ID_database &at[78]
#define BLAST_SEQ_ID_id &at[79]

#define BLAST_MATRIX &at[54]
#define BLAST_MATRIX_is_protein &at[55]
#define BLAST_MATRIX_name &at[56]
#define BLAST_MATRIX_comments &at[57]
#define BLAST_MATRIX_comments_E &at[58]
#define BLAST_MATRIX_row_length &at[59]
#define BLAST_MATRIX_column_length &at[60]
#define BLAST_MATRIX_scores &at[61]
#define BLAST_MATRIX_scores_E &at[62]
#define BLAST_MATRIX_karlinK &at[63]

#define BLAST_QUEUED &at[128]
#define BLAST_QUEUED_length &at[129]

#define BLAST_PROGRESS &at[131]
#define BLAST_PROGRESS_completed &at[132]

#define BLAST_DEFLINE &at[99]
#define BLAST_DEFLINE_id &at[100]
#define BLAST_DEFLINE_defline &at[101]

#define BLAST_VERSION &at[86]
#define BLAST_VERSION_version &at[87]
#define BLAST_VERSION_date &at[88]

#define BLAST_ERROR &at[91]
#define BLAST_ERROR_level &at[92]
#define BLAST_ERROR_msg &at[93]
