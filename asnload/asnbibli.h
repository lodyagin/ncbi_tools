/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "asnbibli.l60";
static AsnValxNodePtr avn;
static AsnTypePtr at;
static AsnModulePtr amp;


/**************************************************
*
*    Defines for Module NCBI-Biblio
*
**************************************************/

#define CIT_ART &at[0]
#define CIT_ART_title &at[1]
#define CIT_ART_authors &at[17]
#define CIT_ART_from &at[51]
#define CIT_ART_from_journal &at[52]
#define CIT_ART_from_book &at[73]
#define CIT_ART_from_proc &at[79]

#define CIT_JOUR &at[53]
#define CIT_JOUR_title &at[54]
#define CIT_JOUR_imp &at[55]

#define CIT_BOOK &at[74]
#define CIT_BOOK_title &at[75]
#define CIT_BOOK_coll &at[76]
#define CIT_BOOK_authors &at[77]
#define CIT_BOOK_imp &at[78]

#define CIT_PAT &at[87]
#define CIT_PAT_title &at[88]
#define CIT_PAT_authors &at[89]
#define CIT_PAT_country &at[90]
#define CIT_PAT_doc_type &at[91]
#define CIT_PAT_number &at[92]
#define CIT_PAT_date_issue &at[93]
#define CIT_PAT_class &at[94]
#define CIT_PAT_class_E &at[95]
#define CIT_PAT_app_number &at[96]
#define CIT_PAT_app_date &at[97]
#define CIT_PAT_applicants &at[98]
#define CIT_PAT_assignees &at[99]
#define CIT_PAT_priority &at[100]
#define CIT_PAT_priority_E &at[101]
#define CIT_PAT_abstract &at[106]

#define CIT_LET &at[107]
#define CIT_LET_cit &at[108]
#define CIT_LET_man_id &at[109]
#define CIT_LET_type &at[110]

#define ID_PAT &at[111]
#define ID_PAT_country &at[112]
#define ID_PAT_id &at[113]
#define ID_PAT_id_number &at[114]
#define ID_PAT_id_app_number &at[115]
#define ID_PAT_doc_type &at[116]

#define CIT_GEN &at[117]
#define CIT_GEN_cit &at[118]
#define CIT_GEN_authors &at[119]
#define CIT_GEN_muid &at[120]
#define CIT_GEN_journal &at[122]
#define CIT_GEN_volume &at[123]
#define CIT_GEN_issue &at[124]
#define CIT_GEN_pages &at[125]
#define CIT_GEN_date &at[126]
#define CIT_GEN_serial_number &at[127]
#define CIT_GEN_title &at[128]
#define CIT_GEN_pmid &at[129]

#define CIT_PROC &at[80]
#define CIT_PROC_book &at[81]
#define CIT_PROC_meet &at[82]

#define CIT_SUB &at[131]
#define CIT_SUB_authors &at[132]
#define CIT_SUB_imp &at[133]
#define CIT_SUB_medium &at[134]
#define CIT_SUB_date &at[135]
#define CIT_SUB_descr &at[136]

#define TITLE &at[2]
#define TITLE_E &at[3]
#define TITLE_E_name &at[4]
#define TITLE_E_tsub &at[6]
#define TITLE_E_trans &at[7]
#define TITLE_E_jta &at[8]
#define TITLE_E_iso_jta &at[9]
#define TITLE_E_ml_jta &at[10]
#define TITLE_E_coden &at[11]
#define TITLE_E_issn &at[12]
#define TITLE_E_abr &at[13]
#define TITLE_E_isbn &at[14]

#define AUTHOR &at[22]
#define AUTHOR_name &at[23]
#define AUTHOR_level &at[25]
#define AUTHOR_role &at[27]
#define AUTHOR_affil &at[28]
#define AUTHOR_is_corr &at[43]

#define PUBMEDID &at[130]

#define AUTH_LIST &at[18]
#define AUTH_LIST_names &at[19]
#define AUTH_LIST_names_std &at[20]
#define AUTH_LIST_names_std_E &at[21]
#define AUTH_LIST_names_ml &at[46]
#define AUTH_LIST_names_ml_E &at[47]
#define AUTH_LIST_names_str &at[48]
#define AUTH_LIST_names_str_E &at[49]
#define AUTH_LIST_affil &at[50]

#define IMPRINT &at[56]
#define IMPRINT_date &at[57]
#define IMPRINT_volume &at[59]
#define IMPRINT_issue &at[60]
#define IMPRINT_pages &at[61]
#define IMPRINT_section &at[62]
#define IMPRINT_pub &at[63]
#define IMPRINT_cprt &at[64]
#define IMPRINT_part_sup &at[65]
#define IMPRINT_language &at[66]
#define IMPRINT_prepub &at[67]
#define IMPRINT_part_supi &at[68]
#define IMPRINT_retract &at[69]

#define MEETING &at[83]
#define MEETING_number &at[84]
#define MEETING_date &at[85]
#define MEETING_place &at[86]

#define PATENT_PRIORITY &at[102]
#define PATENT_PRIORITY_country &at[103]
#define PATENT_PRIORITY_number &at[104]
#define PATENT_PRIORITY_date &at[105]

#define AFFIL &at[29]
#define AFFIL_str &at[30]
#define AFFIL_std &at[31]
#define AFFIL_std_affil &at[32]
#define AFFIL_std_div &at[33]
#define AFFIL_std_city &at[34]
#define AFFIL_std_sub &at[35]
#define AFFIL_std_country &at[36]
#define AFFIL_std_street &at[37]
#define AFFIL_std_email &at[38]
#define AFFIL_std_fax &at[39]
#define AFFIL_std_phone &at[40]
#define AFFIL_std_postal_code &at[41]

#define CITRETRACT &at[70]
#define CITRETRACT_type &at[71]
#define CITRETRACT_exp &at[72]
