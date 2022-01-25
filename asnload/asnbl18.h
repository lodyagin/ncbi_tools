/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "asnbl18.l60";
static AsnValxNodePtr avn;
static AsnTypePtr at;
static AsnModulePtr amp;


/**************************************************
*
*    Defines for Module NCBI-BLAST-1
*
**************************************************/

#define BLAST0_PREFACE &at[0]
#define BLAST0_PREFACE_program &at[1]
#define BLAST0_PREFACE_desc &at[3]
#define BLAST0_PREFACE_version &at[4]
#define BLAST0_PREFACE_dev_date &at[5]
#define BLAST0_PREFACE_bld_date &at[6]
#define BLAST0_PREFACE_cit &at[7]
#define BLAST0_PREFACE_cit_E &at[8]
#define BLAST0_PREFACE_notice &at[10]
#define BLAST0_PREFACE_notice_E &at[11]
#define BLAST0_PREFACE_prog_usage &at[12]
#define BLAST0_PREFACE_prog_usage_E &at[13]
#define BLAST0_PREFACE_susage &at[14]
#define BLAST0_PREFACE_qusage &at[21]

#define BLAST0_JOB_DESC &at[22]
#define BLAST0_JOB_DESC_jid &at[23]
#define BLAST0_JOB_DESC_desc &at[24]
#define BLAST0_JOB_DESC_size &at[25]

#define BLAST0_JOB_PROGRESS &at[27]
#define BLAST0_JOB_PROGRESS_done &at[28]
#define BLAST0_JOB_PROGRESS_positives &at[29]

#define BLAST0_SEQUENCE &at[30]
#define BLAST0_SEQUENCE_desc &at[31]
#define BLAST0_SEQUENCE_desc_E &at[32]
#define BLAST0_SEQUENCE_length &at[42]
#define BLAST0_SEQUENCE_gcode &at[43]
#define BLAST0_SEQUENCE_seq &at[44]

#define BLAST0_KA_BLK &at[50]
#define BLAST0_KA_BLK_matid &at[51]
#define BLAST0_KA_BLK_frames &at[52]
#define BLAST0_KA_BLK_frames_E &at[53]
#define BLAST0_KA_BLK_lambda &at[54]
#define BLAST0_KA_BLK_k &at[56]
#define BLAST0_KA_BLK_h &at[57]

#define BLAST0_DB_DESC &at[58]
#define BLAST0_DB_DESC_name &at[59]
#define BLAST0_DB_DESC_type &at[60]
#define BLAST0_DB_DESC_def &at[61]
#define BLAST0_DB_DESC_rel_date &at[62]
#define BLAST0_DB_DESC_bld_date &at[63]
#define BLAST0_DB_DESC_count &at[64]
#define BLAST0_DB_DESC_totlen &at[65]
#define BLAST0_DB_DESC_maxlen &at[66]

#define BLAST0_RESULT &at[67]
#define BLAST0_RESULT_hist &at[68]
#define BLAST0_RESULT_count &at[79]
#define BLAST0_RESULT_hitlists &at[80]
#define BLAST0_RESULT_hitlists_E &at[81]

#define BLAST0_MATRIX &at[105]
#define BLAST0_MATRIX_matid &at[106]
#define BLAST0_MATRIX_name &at[107]
#define BLAST0_MATRIX_comments &at[108]
#define BLAST0_MATRIX_comments_E &at[109]
#define BLAST0_MATRIX_qalpha &at[110]
#define BLAST0_MATRIX_salpha &at[112]
#define BLAST0_MATRIX_scores &at[113]
#define MATRIX_scores_scaled_ints &at[114]
#define scores_scaled_ints_scale &at[115]
#define MATRIX_scores_scaled_ints_ints &at[116]
#define scores_scaled_ints_ints_E &at[117]
#define BLAST0_MATRIX_scores_reals &at[118]
#define BLAST0_MATRIX_scores_reals_E &at[119]

#define BLAST0_WARNING &at[120]

#define BLAST0_STATUS &at[121]
#define BLAST0_STATUS_code &at[122]
#define BLAST0_STATUS_reason &at[123]

#define BLAST0_OUTBLK &at[124]
#define BLAST0_OUTBLK_E &at[125]
#define BLAST0_OUTBLK_E_preface &at[126]
#define BLAST0_OUTBLK_E_query &at[127]
#define BLAST0_OUTBLK_E_dbdesc &at[128]
#define BLAST0_OUTBLK_E_matrix &at[129]
#define BLAST0_OUTBLK_E_matrix_E &at[130]
#define BLAST0_OUTBLK_E_kablk &at[131]
#define BLAST0_OUTBLK_E_kablk_E &at[132]
#define BLAST0_OUTBLK_E_job_start &at[133]
#define BLAST0_OUTBLK_E_job_progress &at[134]
#define BLAST0_OUTBLK_E_job_done &at[135]
#define BLAST0_OUTBLK_E_result &at[136]
#define BLAST0_OUTBLK_E_parms &at[137]
#define BLAST0_OUTBLK_E_parms_E &at[138]
#define BLAST0_OUTBLK_E_stats &at[139]
#define BLAST0_OUTBLK_E_stats_E &at[140]
#define BLAST0_OUTBLK_E_warning &at[141]
#define BLAST0_OUTBLK_E_status &at[142]

#define BLAST0_REQUEST &at[144]
#define BLAST0_REQUEST_hello &at[145]
#define BLAST0_REQUEST_motd &at[146]
#define BLAST0_REQUEST_prog_info &at[148]
#define BLAST0_REQUEST_usage_info &at[149]
#define BLAST0_REQUEST_db_info &at[150]
#define BLAST0_REQUEST_matrix_info &at[151]
#define BLAST0_REQUEST_matrix_get &at[152]
#define BLAST0_REQUEST_search &at[153]
#define BLAST0_REQUEST_goodbye &at[166]

#define BLAST0_SEARCH &at[154]
#define BLAST0_SEARCH_program &at[155]
#define BLAST0_SEARCH_database &at[156]
#define BLAST0_SEARCH_query &at[157]
#define BLAST0_SEARCH_options &at[158]
#define BLAST0_SEARCH_options_E &at[159]
#define BLAST0_SEARCH_return_matrix &at[160]
#define BLAST0_SEARCH_return_query &at[162]
#define SEARCH_return_BLAST0result &at[163]
#define SEARCH_return_query_seq_in_seg &at[164]
#define SEARCH_return_db_seq_in_seg &at[165]

#define BLAST0_RESPONSE &at[167]
#define BLAST0_RESPONSE_hello &at[168]
#define BLAST0_RESPONSE_motd &at[169]
#define BLAST0_RESPONSE_prog_info &at[170]
#define BLAST0_RESPONSE_prog_info_E &at[171]
#define BLAST0_RESPONSE_usage_info &at[172]
#define BLAST0_RESPONSE_db_info &at[173]
#define BLAST0_RESPONSE_db_info_E &at[174]
#define BLAST0_RESPONSE_matrix_info &at[175]
#define BLAST0_RESPONSE_matrix_info_E &at[176]
#define BLAST0_RESPONSE_ack &at[177]
#define BLAST0_RESPONSE_goodbye &at[181]
#define BLAST0_RESPONSE_queued &at[182]
#define BLAST0_RESPONSE_preface &at[186]
#define BLAST0_RESPONSE_query &at[187]
#define BLAST0_RESPONSE_dbdesc &at[188]
#define BLAST0_RESPONSE_matrix &at[189]
#define BLAST0_RESPONSE_matrix_E &at[190]
#define BLAST0_RESPONSE_kablk &at[191]
#define BLAST0_RESPONSE_kablk_E &at[192]
#define BLAST0_RESPONSE_job_start &at[193]
#define BLAST0_RESPONSE_job_progress &at[194]
#define BLAST0_RESPONSE_job_done &at[195]
#define BLAST0_RESPONSE_score_defs &at[196]
#define BLAST0_RESPONSE_score_defs_E &at[197]
#define BLAST0_RESPONSE_result &at[202]
#define BLAST0_RESPONSE_seqalign &at[203]
#define BLAST0_RESPONSE_parms &at[204]
#define BLAST0_RESPONSE_parms_E &at[205]
#define BLAST0_RESPONSE_stats &at[206]
#define BLAST0_RESPONSE_stats_E &at[207]
#define BLAST0_RESPONSE_warning &at[208]
#define BLAST0_RESPONSE_status &at[209]

#define BLAST0_ACK &at[178]
#define BLAST0_ACK_code &at[179]
#define BLAST0_ACK_reason &at[180]

#define BLAST0_QUEUED &at[183]
#define BLAST0_QUEUED_name &at[184]
#define BLAST0_QUEUED_length &at[185]

#define BLAST0_SCORE_INFO &at[198]
#define BLAST0_SCORE_INFO_sid &at[199]
#define BLAST0_SCORE_INFO_tag &at[200]
#define BLAST0_SCORE_INFO_desc &at[201]

#define BLAST0_SEQ_USAGE &at[15]
#define BLAST0_SEQ_USAGE_raw &at[16]
#define BLAST0_SEQ_USAGE_cooked &at[19]

#define BLAST0_ALPHATYPE &at[17]

#define BLAST0_HISTOGRAM &at[69]
#define BLAST0_HISTOGRAM_expect &at[70]
#define BLAST0_HISTOGRAM_observed &at[71]
#define BLAST0_HISTOGRAM_base &at[72]
#define BLAST0_HISTOGRAM_nbars &at[73]
#define BLAST0_HISTOGRAM_bar &at[74]
#define BLAST0_HISTOGRAM_bar_E &at[75]

#define BLAST0_HITLIST &at[82]
#define BLAST0_HITLIST_count &at[83]
#define BLAST0_HITLIST_kablk &at[84]
#define BLAST0_HITLIST_kablk_E &at[85]
#define BLAST0_HITLIST_hsps &at[86]
#define BLAST0_HITLIST_hsps_E &at[87]
#define BLAST0_HITLIST_seqs &at[103]
#define BLAST0_HITLIST_seqs_E &at[104]

#define BLAST0_HISTOGRAM_BAR &at[76]
#define BLAST0_HISTOGRAM_BAR_x &at[77]
#define BLAST0_HISTOGRAM_BAR_n &at[78]

#define BLAST0_HSP &at[88]
#define BLAST0_HSP_matid &at[89]
#define BLAST0_HSP_scores &at[90]
#define BLAST0_HSP_len &at[92]
#define BLAST0_HSP_segs &at[93]
#define BLAST0_HSP_segs_E &at[94]

#define BLAST0_SEGMENT &at[95]
#define BLAST0_SEGMENT_loc &at[96]
#define BLAST0_SEGMENT_str &at[101]
#define BLAST0_SEGMENT_str_raw &at[102]

#define BLAST0_SEQ_INTERVAL &at[97]
#define BLAST0_SEQ_INTERVAL_strand &at[98]
#define BLAST0_SEQ_INTERVAL_from &at[99]
#define BLAST0_SEQ_INTERVAL_to &at[100]

#define BLAST0_SEQ_DATA &at[45]
#define BLAST0_SEQ_DATA_ncbistdaa &at[46]
#define BLAST0_SEQ_DATA_ncbi2na &at[48]
#define BLAST0_SEQ_DATA_ncbi4na &at[49]

#define BLAST0_SEQ_DESC &at[33]
#define BLAST0_SEQ_DESC_id &at[34]
#define BLAST0_SEQ_DESC_defline &at[41]

#define BLAST0_SEQ_ID &at[35]
#define BLAST0_SEQ_ID_E &at[36]
#define BLAST0_SEQ_ID_E_giid &at[37]
#define BLAST0_SEQ_ID_E_textid &at[38]

#define BLAST0_ALPHA_ID &at[111]
