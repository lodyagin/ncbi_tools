/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "allpub.l60";
static AsnValxNodePtr avn;
static AsnTypePtr at;
static AsnModulePtr amp;


/**************************************************
*
*    Defines for Module NCBI-General
*
**************************************************/

#define DATE &at[0]
#define DATE_str &at[1]
#define DATE_std &at[3]

#define PERSON_ID &at[12]
#define PERSON_ID_dbtag &at[13]
#define PERSON_ID_name &at[20]
#define PERSON_ID_ml &at[29]
#define PERSON_ID_str &at[30]

#define OBJECT_ID &at[17]
#define OBJECT_ID_id &at[18]
#define OBJECT_ID_str &at[19]

#define DBTAG &at[14]
#define DBTAG_db &at[15]
#define DBTAG_tag &at[16]

#define INT_FUZZ &at[31]
#define INT_FUZZ_p_m &at[32]
#define INT_FUZZ_range &at[33]
#define INT_FUZZ_range_max &at[34]
#define INT_FUZZ_range_min &at[35]
#define INT_FUZZ_pct &at[36]
#define INT_FUZZ_lim &at[37]
#define INT_FUZZ_alt &at[39]
#define INT_FUZZ_alt_E &at[40]

#define USER_OBJECT &at[42]
#define USER_OBJECT_class &at[43]
#define USER_OBJECT_type &at[44]
#define USER_OBJECT_data &at[45]
#define USER_OBJECT_data_E &at[46]

#define DATE_STD &at[4]
#define DATE_STD_year &at[5]
#define DATE_STD_month &at[7]
#define DATE_STD_day &at[8]
#define DATE_STD_season &at[9]

#define NAME_STD &at[21]
#define NAME_STD_last &at[22]
#define NAME_STD_first &at[23]
#define NAME_STD_middle &at[24]
#define NAME_STD_full &at[25]
#define NAME_STD_initials &at[26]
#define NAME_STD_suffix &at[27]
#define NAME_STD_title &at[28]

#define USER_FIELD &at[47]
#define USER_FIELD_label &at[48]
#define USER_FIELD_num &at[49]
#define USER_FIELD_data &at[50]
#define USER_FIELD_data_str &at[51]
#define USER_FIELD_data_int &at[52]
#define USER_FIELD_data_real &at[53]
#define USER_FIELD_data_bool &at[55]
#define USER_FIELD_data_os &at[57]
#define USER_FIELD_data_object &at[59]
#define USER_FIELD_data_strs &at[60]
#define USER_FIELD_data_strs_E &at[61]
#define USER_FIELD_data_ints &at[63]
#define USER_FIELD_data_ints_E &at[64]
#define USER_FIELD_data_reals &at[65]
#define USER_FIELD_data_reals_E &at[66]
#define USER_FIELD_data_oss &at[67]
#define USER_FIELD_data_oss_E &at[68]
#define USER_FIELD_data_fields &at[69]
#define USER_FIELD_data_fields_E &at[70]
#define USER_FIELD_data_objects &at[71]
#define USER_FIELD_data_objects_E &at[72]


/**************************************************
*
*    Defines for Module NCBI-Pub
*
**************************************************/

#define PUB &at[73]
#define PUB_gen &at[74]
#define PUB_sub &at[131]
#define PUB_medline &at[155]
#define PUB_muid &at[226]
#define PUB_article &at[227]
#define PUB_journal &at[229]
#define PUB_book &at[231]
#define PUB_proc &at[233]
#define PUB_patent &at[235]
#define PUB_pat_id &at[257]
#define PUB_man &at[265]
#define PUB_equiv &at[271]
#define PUB_pmid &at[274]

#define PUB_SET &at[276]
#define PUB_SET_pub &at[277]
#define PUB_SET_pub_E &at[278]
#define PUB_SET_medline &at[279]
#define PUB_SET_medline_E &at[280]
#define PUB_SET_article &at[281]
#define PUB_SET_article_E &at[282]
#define PUB_SET_journal &at[283]
#define PUB_SET_journal_E &at[284]
#define PUB_SET_book &at[285]
#define PUB_SET_book_E &at[286]
#define PUB_SET_proc &at[287]
#define PUB_SET_proc_E &at[288]
#define PUB_SET_patent &at[289]
#define PUB_SET_patent_E &at[290]

#define PUB_EQUIV &at[272]
#define PUB_EQUIV_E &at[273]


/**************************************************
*
*    Defines for Module NCBI-Biblio
*
**************************************************/

#define CIT_ART &at[163]
#define CIT_ART_title &at[164]
#define CIT_ART_authors &at[165]
#define CIT_ART_from &at[166]
#define CIT_ART_from_journal &at[167]
#define CIT_ART_from_book &at[171]
#define CIT_ART_from_proc &at[177]

#define CIT_JOUR &at[168]
#define CIT_JOUR_title &at[169]
#define CIT_JOUR_imp &at[170]

#define CIT_BOOK &at[172]
#define CIT_BOOK_title &at[173]
#define CIT_BOOK_coll &at[174]
#define CIT_BOOK_authors &at[175]
#define CIT_BOOK_imp &at[176]

#define CIT_PAT &at[237]
#define CIT_PAT_title &at[238]
#define CIT_PAT_authors &at[239]
#define CIT_PAT_country &at[240]
#define CIT_PAT_doc_type &at[241]
#define CIT_PAT_number &at[242]
#define CIT_PAT_date_issue &at[243]
#define CIT_PAT_class &at[244]
#define CIT_PAT_class_E &at[245]
#define CIT_PAT_app_number &at[246]
#define CIT_PAT_app_date &at[247]
#define CIT_PAT_applicants &at[248]
#define CIT_PAT_assignees &at[249]
#define CIT_PAT_priority &at[250]
#define CIT_PAT_priority_E &at[251]
#define CIT_PAT_abstract &at[256]

#define CIT_LET &at[267]
#define CIT_LET_cit &at[268]
#define CIT_LET_man_id &at[269]
#define CIT_LET_type &at[270]

#define ID_PAT &at[259]
#define ID_PAT_country &at[260]
#define ID_PAT_id &at[261]
#define ID_PAT_id_number &at[262]
#define ID_PAT_id_app_number &at[263]
#define ID_PAT_doc_type &at[264]

#define CIT_GEN &at[76]
#define CIT_GEN_cit &at[77]
#define CIT_GEN_authors &at[78]
#define CIT_GEN_muid &at[108]
#define CIT_GEN_journal &at[109]
#define CIT_GEN_volume &at[122]
#define CIT_GEN_issue &at[123]
#define CIT_GEN_pages &at[124]
#define CIT_GEN_date &at[125]
#define CIT_GEN_serial_number &at[127]
#define CIT_GEN_title &at[128]
#define CIT_GEN_pmid &at[129]

#define CIT_PROC &at[178]
#define CIT_PROC_book &at[179]
#define CIT_PROC_meet &at[180]

#define CIT_SUB &at[133]
#define CIT_SUB_authors &at[134]
#define CIT_SUB_imp &at[135]
#define CIT_SUB_medium &at[152]
#define CIT_SUB_date &at[153]
#define CIT_SUB_descr &at[154]

#define TITLE &at[110]
#define TITLE_E &at[111]
#define TITLE_E_name &at[112]
#define TITLE_E_tsub &at[113]
#define TITLE_E_trans &at[114]
#define TITLE_E_jta &at[115]
#define TITLE_E_iso_jta &at[116]
#define TITLE_E_ml_jta &at[117]
#define TITLE_E_coden &at[118]
#define TITLE_E_issn &at[119]
#define TITLE_E_abr &at[120]
#define TITLE_E_isbn &at[121]

#define AUTHOR &at[83]
#define AUTHOR_name &at[84]
#define AUTHOR_level &at[86]
#define AUTHOR_role &at[87]
#define AUTHOR_affil &at[88]
#define AUTHOR_is_corr &at[102]

#define PUBMEDID &at[130]

#define AUTH_LIST &at[79]
#define AUTH_LIST_names &at[80]
#define AUTH_LIST_names_std &at[81]
#define AUTH_LIST_names_std_E &at[82]
#define AUTH_LIST_names_ml &at[103]
#define AUTH_LIST_names_ml_E &at[104]
#define AUTH_LIST_names_str &at[105]
#define AUTH_LIST_names_str_E &at[106]
#define AUTH_LIST_affil &at[107]

#define IMPRINT &at[136]
#define IMPRINT_date &at[137]
#define IMPRINT_volume &at[138]
#define IMPRINT_issue &at[139]
#define IMPRINT_pages &at[140]
#define IMPRINT_section &at[141]
#define IMPRINT_pub &at[142]
#define IMPRINT_cprt &at[143]
#define IMPRINT_part_sup &at[144]
#define IMPRINT_language &at[145]
#define IMPRINT_prepub &at[146]
#define IMPRINT_part_supi &at[147]
#define IMPRINT_retract &at[148]

#define MEETING &at[181]
#define MEETING_number &at[182]
#define MEETING_date &at[183]
#define MEETING_place &at[184]

#define PATENT_PRIORITY &at[252]
#define PATENT_PRIORITY_country &at[253]
#define PATENT_PRIORITY_number &at[254]
#define PATENT_PRIORITY_date &at[255]

#define AFFIL &at[89]
#define AFFIL_str &at[90]
#define AFFIL_std &at[91]
#define AFFIL_std_affil &at[92]
#define AFFIL_std_div &at[93]
#define AFFIL_std_city &at[94]
#define AFFIL_std_sub &at[95]
#define AFFIL_std_country &at[96]
#define AFFIL_std_street &at[97]
#define AFFIL_std_email &at[98]
#define AFFIL_std_fax &at[99]
#define AFFIL_std_phone &at[100]
#define AFFIL_std_postal_code &at[101]

#define CITRETRACT &at[149]
#define CITRETRACT_type &at[150]
#define CITRETRACT_exp &at[151]


/**************************************************
*
*    Defines for Module NCBI-Medline
*
**************************************************/

#define MEDLINE_ENTRY &at[157]
#define MEDLINE_ENTRY_uid &at[158]
#define MEDLINE_ENTRY_em &at[159]
#define MEDLINE_ENTRY_cit &at[161]
#define MEDLINE_ENTRY_abstract &at[185]
#define MEDLINE_ENTRY_mesh &at[186]
#define MEDLINE_ENTRY_mesh_E &at[187]
#define MEDLINE_ENTRY_substance &at[196]
#define MEDLINE_ENTRY_substance_E &at[197]
#define MEDLINE_ENTRY_xref &at[202]
#define MEDLINE_ENTRY_xref_E &at[203]
#define MEDLINE_ENTRY_idnum &at[207]
#define MEDLINE_ENTRY_idnum_E &at[208]
#define MEDLINE_ENTRY_gene &at[209]
#define MEDLINE_ENTRY_gene_E &at[210]
#define MEDLINE_ENTRY_pmid &at[211]
#define MEDLINE_ENTRY_pub_type &at[213]
#define MEDLINE_ENTRY_pub_type_E &at[214]
#define MEDLINE_ENTRY_mlfield &at[215]
#define MEDLINE_ENTRY_mlfield_E &at[216]
#define MEDLINE_ENTRY_status &at[225]

#define MEDLINE_SI &at[204]
#define MEDLINE_SI_type &at[205]
#define MEDLINE_SI_cit &at[206]

#define MEDLINE_MESH &at[188]
#define MEDLINE_MESH_mp &at[189]
#define MEDLINE_MESH_term &at[190]
#define MEDLINE_MESH_qual &at[191]
#define MEDLINE_MESH_qual_E &at[192]

#define MEDLINE_RN &at[198]
#define MEDLINE_RN_type &at[199]
#define MEDLINE_RN_cit &at[200]
#define MEDLINE_RN_name &at[201]

#define MEDLINE_FIELD &at[217]
#define MEDLINE_FIELD_type &at[218]
#define MEDLINE_FIELD_str &at[219]
#define MEDLINE_FIELD_ids &at[220]
#define MEDLINE_FIELD_ids_E &at[221]

#define MEDLINE_QUAL &at[193]
#define MEDLINE_QUAL_mp &at[194]
#define MEDLINE_QUAL_subh &at[195]

#define DOCREF &at[222]
#define DOCREF_type &at[223]
#define DOCREF_uid &at[224]
