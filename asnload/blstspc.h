/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "blstspc.l14";
static AsnValxNodePtr avn;
static AsnTypePtr at;
static AsnModulePtr amp;


/**************************************************
*
*    Defines for Module NCBI-Blast
*
**************************************************/

#define BLAST_REQUEST &at[0]
#define BLAST_REQUEST_init &at[1]
#define BLAST_REQUEST_motd &at[3]
#define BLAST_REQUEST_db_info &at[5]
#define BLAST_REQUEST_db_info_specific &at[6]
#define BLAST_REQUEST_matrix_get &at[12]
#define BLAST_REQUEST_search &at[13]
#define BLAST_REQUEST_db_seq_get &at[67]
#define BLAST_REQUEST_db_redundant_ids_get &at[73]
#define BLAST_REQUEST_db_redundant_descr_get &at[74]
#define BLAST_REQUEST_fini &at[75]

#define BLAST_RESPONSE &at[76]
#define BLAST_RESPONSE_init &at[77]
#define BLAST_RESPONSE_motd &at[81]
#define BLAST_RESPONSE_error &at[82]
#define BLAST_RESPONSE_db_seq_get &at[86]
#define BLAST_RESPONSE_db_redundant_ids_get &at[87]
#define BLAST_RESPONSE_db_redundant_ids_get_E &at[88]
#define BLAST_RESPONSE_db_redundant_descr_get &at[89]
#define BLAST_RESPONSE_db_redundant_descr_get_E &at[90]
#define BLAST_RESPONSE_db_info &at[94]
#define BLAST_RESPONSE_db_info_E &at[95]
#define BLAST_RESPONSE_db_info_specific &at[103]
#define BLAST_RESPONSE_matrix &at[104]
#define BLAST_RESPONSE_alignment &at[105]
#define BLAST_RESPONSE_mask &at[107]
#define BLAST_RESPONSE_kablk &at[112]
#define BLAST_RESPONSE_parameters &at[118]
#define BLAST_RESPONSE_queued &at[119]
#define BLAST_RESPONSE_start &at[122]
#define BLAST_RESPONSE_progress &at[125]
#define BLAST_RESPONSE_done &at[126]
#define BLAST_RESPONSE_fini &at[127]

#define BLAST_PARAMETERS &at[20]
#define BLAST_PARAMETERS_first_threshold &at[21]
#define BLAST_PARAMETERS_second_threshold &at[23]
#define BLAST_PARAMETERS_cutoff &at[24]
#define BLAST_PARAMETERS_cutoff_evalue &at[25]
#define BLAST_PARAMETERS_cutoff_score &at[27]
#define BLAST_PARAMETERS_cutoff2 &at[29]
#define BLAST_PARAMETERS_cutoff2_evalue &at[30]
#define BLAST_PARAMETERS_cutoff2_score &at[31]
#define BLAST_PARAMETERS_hitlist_size &at[32]
#define BLAST_PARAMETERS_nucl_penalty &at[33]
#define BLAST_PARAMETERS_nucl_reward &at[34]
#define BLAST_PARAMETERS_genetic_code &at[35]
#define BLAST_PARAMETERS_db_genetic_code &at[36]
#define BLAST_PARAMETERS_low_complexity_filtering &at[37]
#define BLAST_PARAMETERS_gapped_alignment &at[38]
#define BLAST_PARAMETERS_gap_open &at[40]
#define BLAST_PARAMETERS_gap_extend &at[41]
#define BLAST_PARAMETERS_required_start &at[42]
#define BLAST_PARAMETERS_required_end &at[43]
#define BLAST_PARAMETERS_ethresh &at[44]
#define BLAST_PARAMETERS_max_num_passes &at[45]
#define BLAST_PARAMETERS_pseudo_count_const &at[46]
#define BLAST_PARAMETERS_other_options &at[47]
#define BLAST_PARAMETERS_gilist &at[48]
#define BLAST_PARAMETERS_gilist_E &at[49]
#define BLAST_PARAMETERS_gifile &at[51]
#define BLAST_PARAMETERS_matrix &at[52]
#define BLAST_PARAMETERS_filter_string &at[53]

#define BLAST_DBINFO_GET &at[7]
#define BLAST_DBINFO_GET_name &at[8]
#define BLAST_DBINFO_GET_type &at[9]

#define BLAST_SEARCH &at[14]
#define BLAST_SEARCH_program &at[15]
#define BLAST_SEARCH_query &at[16]
#define BLAST_SEARCH_database &at[18]
#define BLAST_SEARCH_parameters &at[19]
#define BLAST_SEARCH_mask &at[54]
#define BLAST_SEARCH_matrix &at[56]

#define BLAST_SEQ_ID &at[68]
#define BLAST_SEQ_ID_is_protein &at[69]
#define BLAST_SEQ_ID_database &at[70]
#define BLAST_SEQ_ID_id &at[71]

#define BLAST_MATRIX &at[57]
#define BLAST_MATRIX_is_protein &at[58]
#define BLAST_MATRIX_name &at[59]
#define BLAST_MATRIX_comments &at[60]
#define BLAST_MATRIX_comments_E &at[61]
#define BLAST_MATRIX_row_length &at[62]
#define BLAST_MATRIX_column_length &at[63]
#define BLAST_MATRIX_scores &at[64]
#define BLAST_MATRIX_scores_E &at[65]
#define BLAST_MATRIX_karlinK &at[66]

#define BLAST_DBINFO &at[96]
#define BLAST_DBINFO_is_protein &at[97]
#define BLAST_DBINFO_name &at[98]
#define BLAST_DBINFO_definition &at[99]
#define BLAST_DBINFO_date &at[100]
#define BLAST_DBINFO_total_length &at[101]
#define BLAST_DBINFO_number_seqs &at[102]

#define BLAST_QUEUED &at[120]
#define BLAST_QUEUED_length &at[121]

#define BLAST_PROGRESS &at[123]
#define BLAST_PROGRESS_completed &at[124]

#define BLAST_KABLK &at[113]
#define BLAST_KABLK_lambda &at[114]
#define BLAST_KABLK_k &at[115]
#define BLAST_KABLK_h &at[116]
#define BLAST_KABLK_gapped &at[117]

#define BLAST_DEFLINE &at[91]
#define BLAST_DEFLINE_id &at[92]
#define BLAST_DEFLINE_defline &at[93]

#define BLAST_MASK &at[108]
#define BLAST_MASK_location &at[109]
#define BLAST_MASK_location_E &at[110]
#define BLAST_MASK_frame &at[111]

#define BLAST_VERSION &at[78]
#define BLAST_VERSION_version &at[79]
#define BLAST_VERSION_date &at[80]

#define BLAST_ERROR &at[83]
#define BLAST_ERROR_level &at[84]
#define BLAST_ERROR_msg &at[85]
