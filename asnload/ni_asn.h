/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "ni_asn.l60";
static AsnValxNodePtr avn;
static AsnTypePtr at;
static AsnModulePtr amp;


/**************************************************
*
*    Defines for Module NCBI-MESSAGE
*
**************************************************/

#define SVC_ENTRY &at[0]
#define SVC_ENTRY_name &at[1]
#define SVC_ENTRY_minvers &at[3]
#define SVC_ENTRY_maxvers &at[5]
#define SVC_ENTRY_id &at[6]
#define SVC_ENTRY_priority &at[7]
#define SVC_ENTRY_group &at[8]
#define SVC_ENTRY_description &at[9]
#define SVC_ENTRY_types &at[10]
#define SVC_ENTRY_types_E &at[11]
#define SVC_ENTRY_priority_timeout &at[13]
#define SVC_ENTRY_priority_penalty &at[14]
#define SVC_ENTRY_encryption_supported &at[15]
#define SVC_ENTRY_tracking_period &at[17]
#define SVC_ENTRY_tracking_count &at[18]

#define RES_ENTRY &at[20]
#define RES_ENTRY_name &at[21]
#define RES_ENTRY_type &at[22]
#define RES_ENTRY_minvers &at[23]
#define RES_ENTRY_maxvers &at[24]
#define RES_ENTRY_id &at[25]
#define RES_ENTRY_group &at[26]
#define RES_ENTRY_description &at[27]

#define TOOLSET &at[28]
#define TOOLSET_host &at[29]
#define TOOLSET_motd &at[30]
#define TOOLSET_services &at[31]
#define TOOLSET_services_E &at[32]
#define TOOLSET_resources &at[33]
#define TOOLSET_resources_E &at[34]
#define TOOLSET_regions &at[35]
#define TOOLSET_regions_E &at[36]

#define IDENTITY &at[40]
#define IDENTITY_username &at[41]
#define IDENTITY_group &at[42]
#define IDENTITY_domain &at[43]

#define REQUEST &at[44]
#define REQUEST_address &at[45]
#define REQUEST_port &at[46]
#define REQUEST_svcentry &at[47]
#define REQUEST_resentry &at[48]
#define REQUEST_resentry_E &at[49]

#define MSG_ACK &at[50]
#define MSG_ACK_seqno &at[51]
#define MSG_ACK_disp_info &at[52]
#define MSG_ACK_admin_info &at[66]
#define MSG_ACK_motd &at[67]

#define MSG_NACK &at[68]
#define MSG_NACK_seqno &at[69]
#define MSG_NACK_code &at[70]
#define MSG_NACK_reason &at[71]
#define MSG_NACK_disp_info &at[72]

#define MSG_LOGIN &at[73]
#define MSG_LOGIN_seqno &at[74]
#define MSG_LOGIN_uid &at[75]
#define MSG_LOGIN_password &at[76]
#define MSG_LOGIN_disp_serial_no &at[77]
#define MSG_LOGIN_encryption_desired &at[78]
#define MSG_LOGIN_pub_key &at[79]
#define MSG_LOGIN_des_key &at[80]
#define MSG_LOGIN_connect_delay &at[81]
#define MSG_LOGIN_server_port &at[82]

#define MSG_SVC_LIST &at[83]
#define MSG_SVC_LIST_seqno &at[84]
#define MSG_SVC_LIST_toollist &at[85]
#define MSG_SVC_LIST_knows_tracking &at[86]

#define MSG_SVC_REQUEST &at[87]
#define MSG_SVC_REQUEST_seqno &at[88]
#define MSG_SVC_REQUEST_conid &at[89]
#define MSG_SVC_REQUEST_uid &at[90]
#define MSG_SVC_REQUEST_request &at[91]
#define MSG_SVC_REQUEST_platform &at[92]
#define MSG_SVC_REQUEST_appl_id &at[93]
#define MSG_SVC_REQUEST_des_key &at[94]
#define SVC_REQUEST_want_pre_response &at[95]
#define MSG_SVC_REQUEST_server_ip &at[96]
#define MSG_SVC_REQUEST_server_port &at[97]
#define MSG_SVC_REQUEST_want_ticket &at[98]
#define MSG_SVC_REQUEST_ticket &at[99]

#define MSG_SVC_RESPONSE &at[109]
#define MSG_SVC_RESPONSE_seqno &at[110]
#define MSG_SVC_RESPONSE_request &at[111]

#define MSG_CMD &at[112]
#define MSG_CMD_seqno &at[113]
#define MSG_CMD_command &at[114]

#define MSG_ACCT &at[115]
#define MSG_ACCT_seqno &at[116]
#define MSG_ACCT_conid &at[117]
#define MSG_ACCT_jobname &at[118]
#define MSG_ACCT_usertime &at[119]
#define MSG_ACCT_systemtime &at[120]

#define MSG_CATALOG &at[121]
#define MSG_CATALOG_seqno &at[122]
#define MSG_CATALOG_motd &at[123]
#define MSG_CATALOG_toollists &at[124]
#define MSG_CATALOG_toollists_E &at[125]

#define MESSAGE &at[126]
#define MESSAGE_ack &at[127]
#define MESSAGE_nack &at[128]
#define MESSAGE_login &at[129]
#define MESSAGE_svc_list &at[130]
#define MESSAGE_svc_request &at[131]
#define MESSAGE_svc_response &at[132]
#define MESSAGE_command &at[133]
#define MESSAGE_acct &at[134]
#define MESSAGE_catalog &at[135]
#define MESSAGE_svc_pre_response &at[136]
#define MESSAGE_load_status &at[140]

#define REGION_DESCR &at[37]
#define REGION_DESCR_region_name &at[38]
#define REGION_DESCR_priority_delta &at[39]

#define RSA_PUBKEY &at[61]
#define RSA_PUBKEY_bits &at[62]
#define RSA_PUBKEY_modulus &at[63]
#define RSA_PUBKEY_exponent &at[65]

#define DISPATCHER_INFO &at[53]
#define DISPATCHER_INFO_serial_no &at[54]
#define INFO_is_alternate_list &at[55]
#define DISPATCHER_INFO_num_dispatchers &at[56]
#define DISPATCHER_INFO_disp_list &at[57]
#define DISPATCHER_INFO_disp_list_E &at[58]
#define DISPATCHER_INFO_pub_key &at[60]

#define TICKET &at[100]
#define TICKET_seqno &at[101]
#define TICKET_confounding_rand_num &at[102]
#define TICKET_client_ip_1 &at[103]
#define TICKET_client_ip_2 &at[104]
#define TICKET_server_ip &at[105]
#define TICKET_client_des_key &at[106]
#define TICKET_ticket_expiration &at[107]
#define TICKET_checksum &at[108]

#define MSG_SVC_PRE_RESPONSE &at[137]
#define MSG_SVC_PRE_RESPONSE_seqno &at[138]
#define MSG_SVC_PRE_RESPONSE_server_ip &at[139]

#define MSG_LOAD_STATUS &at[141]
#define MSG_LOAD_STATUS_load &at[142]
#define MSG_LOAD_STATUS_power &at[144]
#define MSG_LOAD_STATUS_light_thresh &at[145]
#define MSG_LOAD_STATUS_heavy_thresh &at[146]
#define MSG_LOAD_STATUS_job_penalty &at[147]
