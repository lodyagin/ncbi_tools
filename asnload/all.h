/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "all.l60";
static AsnValxNodePtr avn;
static AsnTypePtr at;
static AsnModulePtr amp;


/**************************************************
*
*    Defines for Module NCBI-Access
*
**************************************************/

#define LINK_SET &at[0]
#define LINK_SET_num &at[1]
#define LINK_SET_uids &at[3]
#define LINK_SET_uids_E &at[4]
#define LINK_SET_weights &at[6]
#define LINK_SET_weights_E &at[7]


/**************************************************
*
*    Defines for Module NCBI-Biblio
*
**************************************************/

#define CIT_ART &at[9]
#define CIT_ART_title &at[10]
#define CIT_ART_authors &at[26]
#define CIT_ART_from &at[77]
#define CIT_ART_from_journal &at[78]
#define CIT_ART_from_book &at[107]
#define CIT_ART_from_proc &at[113]

#define CIT_JOUR &at[79]
#define CIT_JOUR_title &at[80]
#define CIT_JOUR_imp &at[81]

#define CIT_BOOK &at[108]
#define CIT_BOOK_title &at[109]
#define CIT_BOOK_coll &at[110]
#define CIT_BOOK_authors &at[111]
#define CIT_BOOK_imp &at[112]

#define CIT_PAT &at[121]
#define CIT_PAT_title &at[122]
#define CIT_PAT_authors &at[123]
#define CIT_PAT_country &at[124]
#define CIT_PAT_doc_type &at[125]
#define CIT_PAT_number &at[126]
#define CIT_PAT_date_issue &at[127]
#define CIT_PAT_class &at[128]
#define CIT_PAT_class_E &at[129]
#define CIT_PAT_app_number &at[130]
#define CIT_PAT_app_date &at[131]
#define CIT_PAT_applicants &at[132]
#define CIT_PAT_assignees &at[133]
#define CIT_PAT_priority &at[134]
#define CIT_PAT_priority_E &at[135]
#define CIT_PAT_abstract &at[140]

#define CIT_LET &at[141]
#define CIT_LET_cit &at[142]
#define CIT_LET_man_id &at[143]
#define CIT_LET_type &at[144]

#define ID_PAT &at[145]
#define ID_PAT_country &at[146]
#define ID_PAT_id &at[147]
#define ID_PAT_id_number &at[148]
#define ID_PAT_id_app_number &at[149]
#define ID_PAT_doc_type &at[150]

#define CIT_GEN &at[151]
#define CIT_GEN_cit &at[152]
#define CIT_GEN_authors &at[153]
#define CIT_GEN_muid &at[154]
#define CIT_GEN_journal &at[155]
#define CIT_GEN_volume &at[156]
#define CIT_GEN_issue &at[157]
#define CIT_GEN_pages &at[158]
#define CIT_GEN_date &at[159]
#define CIT_GEN_serial_number &at[160]
#define CIT_GEN_title &at[161]
#define CIT_GEN_pmid &at[162]

#define CIT_PROC &at[114]
#define CIT_PROC_book &at[115]
#define CIT_PROC_meet &at[116]

#define CIT_SUB &at[164]
#define CIT_SUB_authors &at[165]
#define CIT_SUB_imp &at[166]
#define CIT_SUB_medium &at[167]
#define CIT_SUB_date &at[168]
#define CIT_SUB_descr &at[169]

#define TITLE &at[11]
#define TITLE_E &at[12]
#define TITLE_E_name &at[13]
#define TITLE_E_tsub &at[15]
#define TITLE_E_trans &at[16]
#define TITLE_E_jta &at[17]
#define TITLE_E_iso_jta &at[18]
#define TITLE_E_ml_jta &at[19]
#define TITLE_E_coden &at[20]
#define TITLE_E_issn &at[21]
#define TITLE_E_abr &at[22]
#define TITLE_E_isbn &at[23]

#define AUTHOR &at[31]
#define AUTHOR_name &at[32]
#define AUTHOR_level &at[53]
#define AUTHOR_role &at[55]
#define AUTHOR_affil &at[56]
#define AUTHOR_is_corr &at[70]

#define PUBMEDID &at[163]

#define AUTH_LIST &at[27]
#define AUTH_LIST_names &at[28]
#define AUTH_LIST_names_std &at[29]
#define AUTH_LIST_names_std_E &at[30]
#define AUTH_LIST_names_ml &at[72]
#define AUTH_LIST_names_ml_E &at[73]
#define AUTH_LIST_names_str &at[74]
#define AUTH_LIST_names_str_E &at[75]
#define AUTH_LIST_affil &at[76]

#define IMPRINT &at[82]
#define IMPRINT_date &at[83]
#define IMPRINT_volume &at[93]
#define IMPRINT_issue &at[94]
#define IMPRINT_pages &at[95]
#define IMPRINT_section &at[96]
#define IMPRINT_pub &at[97]
#define IMPRINT_cprt &at[98]
#define IMPRINT_part_sup &at[99]
#define IMPRINT_language &at[100]
#define IMPRINT_prepub &at[101]
#define IMPRINT_part_supi &at[102]
#define IMPRINT_retract &at[103]

#define MEETING &at[117]
#define MEETING_number &at[118]
#define MEETING_date &at[119]
#define MEETING_place &at[120]

#define PATENT_PRIORITY &at[136]
#define PATENT_PRIORITY_country &at[137]
#define PATENT_PRIORITY_number &at[138]
#define PATENT_PRIORITY_date &at[139]

#define AFFIL &at[57]
#define AFFIL_str &at[58]
#define AFFIL_std &at[59]
#define AFFIL_std_affil &at[60]
#define AFFIL_std_div &at[61]
#define AFFIL_std_city &at[62]
#define AFFIL_std_sub &at[63]
#define AFFIL_std_country &at[64]
#define AFFIL_std_street &at[65]
#define AFFIL_std_email &at[66]
#define AFFIL_std_fax &at[67]
#define AFFIL_std_phone &at[68]
#define AFFIL_std_postal_code &at[69]

#define CITRETRACT &at[104]
#define CITRETRACT_type &at[105]
#define CITRETRACT_exp &at[106]


/**************************************************
*
*    Defines for Module NCBI-FeatDef
*
**************************************************/

#define FEATDEF &at[170]
#define FEATDEF_typelabel &at[171]
#define FEATDEF_menulabel &at[172]
#define FEATDEF_featdef_key &at[173]
#define FEATDEF_seqfeat_key &at[174]
#define FEATDEF_entrygroup &at[175]
#define FEATDEF_displaygroup &at[176]
#define FEATDEF_molgroup &at[177]

#define FEATDEFSET &at[179]
#define FEATDEFSET_E &at[180]

#define FEATDISPGROUP &at[181]
#define FEATDISPGROUP_groupkey &at[182]
#define FEATDISPGROUP_groupname &at[183]

#define FEATDISPGROUPSET &at[184]
#define FEATDISPGROUPSET_E &at[185]

#define FEATMOLTYPE &at[178]

#define FEATDEFGROUPSET &at[186]
#define FEATDEFGROUPSET_groups &at[187]
#define FEATDEFGROUPSET_defs &at[188]


/**************************************************
*
*    Defines for Module NCBI-General
*
**************************************************/

#define DATE &at[85]
#define DATE_str &at[86]
#define DATE_std &at[87]

#define PERSON_ID &at[34]
#define PERSON_ID_dbtag &at[35]
#define PERSON_ID_name &at[42]
#define PERSON_ID_ml &at[51]
#define PERSON_ID_str &at[52]

#define OBJECT_ID &at[39]
#define OBJECT_ID_id &at[40]
#define OBJECT_ID_str &at[41]

#define DBTAG &at[36]
#define DBTAG_db &at[37]
#define DBTAG_tag &at[38]

#define INT_FUZZ &at[189]
#define INT_FUZZ_p_m &at[190]
#define INT_FUZZ_range &at[191]
#define INT_FUZZ_range_max &at[192]
#define INT_FUZZ_range_min &at[193]
#define INT_FUZZ_pct &at[194]
#define INT_FUZZ_lim &at[195]
#define INT_FUZZ_alt &at[196]
#define INT_FUZZ_alt_E &at[197]

#define USER_OBJECT &at[198]
#define USER_OBJECT_class &at[199]
#define USER_OBJECT_type &at[200]
#define USER_OBJECT_data &at[201]
#define USER_OBJECT_data_E &at[202]

#define DATE_STD &at[88]
#define DATE_STD_year &at[89]
#define DATE_STD_month &at[90]
#define DATE_STD_day &at[91]
#define DATE_STD_season &at[92]

#define NAME_STD &at[43]
#define NAME_STD_last &at[44]
#define NAME_STD_first &at[45]
#define NAME_STD_middle &at[46]
#define NAME_STD_full &at[47]
#define NAME_STD_initials &at[48]
#define NAME_STD_suffix &at[49]
#define NAME_STD_title &at[50]

#define USER_FIELD &at[203]
#define USER_FIELD_label &at[204]
#define USER_FIELD_num &at[205]
#define USER_FIELD_data &at[206]
#define USER_FIELD_data_str &at[207]
#define USER_FIELD_data_int &at[208]
#define USER_FIELD_data_real &at[209]
#define USER_FIELD_data_bool &at[211]
#define USER_FIELD_data_os &at[212]
#define USER_FIELD_data_object &at[214]
#define USER_FIELD_data_strs &at[215]
#define USER_FIELD_data_strs_E &at[216]
#define USER_FIELD_data_ints &at[217]
#define USER_FIELD_data_ints_E &at[218]
#define USER_FIELD_data_reals &at[219]
#define USER_FIELD_data_reals_E &at[220]
#define USER_FIELD_data_oss &at[221]
#define USER_FIELD_data_oss_E &at[222]
#define USER_FIELD_data_fields &at[223]
#define USER_FIELD_data_fields_E &at[224]
#define USER_FIELD_data_objects &at[225]
#define USER_FIELD_data_objects_E &at[226]


/**************************************************
*
*    Defines for Module NCBI-Medlars
*
**************************************************/

#define MEDLARS_ENTRY &at[227]
#define MEDLARS_ENTRY_pmid &at[228]
#define MEDLARS_ENTRY_muid &at[230]
#define MEDLARS_ENTRY_recs &at[231]
#define MEDLARS_ENTRY_recs_E &at[232]

#define MEDLARS_RECORD &at[233]
#define MEDLARS_RECORD_code &at[234]
#define MEDLARS_RECORD_abbr &at[235]
#define MEDLARS_RECORD_data &at[236]


/**************************************************
*
*    Defines for Module NCBI-Medline
*
**************************************************/

#define MEDLINE_ENTRY &at[237]
#define MEDLINE_ENTRY_uid &at[238]
#define MEDLINE_ENTRY_em &at[239]
#define MEDLINE_ENTRY_cit &at[241]
#define MEDLINE_ENTRY_abstract &at[243]
#define MEDLINE_ENTRY_mesh &at[244]
#define MEDLINE_ENTRY_mesh_E &at[245]
#define MEDLINE_ENTRY_substance &at[254]
#define MEDLINE_ENTRY_substance_E &at[255]
#define MEDLINE_ENTRY_xref &at[260]
#define MEDLINE_ENTRY_xref_E &at[261]
#define MEDLINE_ENTRY_idnum &at[265]
#define MEDLINE_ENTRY_idnum_E &at[266]
#define MEDLINE_ENTRY_gene &at[267]
#define MEDLINE_ENTRY_gene_E &at[268]
#define MEDLINE_ENTRY_pmid &at[269]
#define MEDLINE_ENTRY_pub_type &at[271]
#define MEDLINE_ENTRY_pub_type_E &at[272]
#define MEDLINE_ENTRY_mlfield &at[273]
#define MEDLINE_ENTRY_mlfield_E &at[274]
#define MEDLINE_ENTRY_status &at[283]

#define MEDLINE_SI &at[262]
#define MEDLINE_SI_type &at[263]
#define MEDLINE_SI_cit &at[264]

#define MEDLINE_MESH &at[246]
#define MEDLINE_MESH_mp &at[247]
#define MEDLINE_MESH_term &at[248]
#define MEDLINE_MESH_qual &at[249]
#define MEDLINE_MESH_qual_E &at[250]

#define MEDLINE_RN &at[256]
#define MEDLINE_RN_type &at[257]
#define MEDLINE_RN_cit &at[258]
#define MEDLINE_RN_name &at[259]

#define MEDLINE_FIELD &at[275]
#define MEDLINE_FIELD_type &at[276]
#define MEDLINE_FIELD_str &at[277]
#define MEDLINE_FIELD_ids &at[278]
#define MEDLINE_FIELD_ids_E &at[279]

#define MEDLINE_QUAL &at[251]
#define MEDLINE_QUAL_mp &at[252]
#define MEDLINE_QUAL_subh &at[253]

#define DOCREF &at[280]
#define DOCREF_type &at[281]
#define DOCREF_uid &at[282]


/**************************************************
*
*    Defines for Module NCBI-Mime
*
**************************************************/

#define NCBI_MIME_ASN1 &at[284]
#define NCBI_MIME_ASN1_entrez &at[285]
#define NCBI_MIME_ASN1_alignstruc &at[1058]
#define NCBI_MIME_ASN1_alignseq &at[1069]
#define NCBI_MIME_ASN1_strucseq &at[1075]

#define ENTREZ_GENERAL &at[286]
#define ENTREZ_GENERAL_title &at[287]
#define ENTREZ_GENERAL_data &at[288]
#define ENTREZ_GENERAL_data_ml &at[289]
#define ENTREZ_GENERAL_data_prot &at[291]
#define ENTREZ_GENERAL_data_nuc &at[1049]
#define ENTREZ_GENERAL_data_genome &at[1050]
#define ENTREZ_GENERAL_data_structure &at[1051]
#define ENTREZ_GENERAL_data_strucAnnot &at[1053]
#define ENTREZ_GENERAL_style &at[1055]
#define ENTREZ_GENERAL_location &at[1057]

#define BIOSTRUC_ALIGN &at[1059]
#define BIOSTRUC_ALIGN_master &at[1060]
#define BIOSTRUC_ALIGN_slaves &at[1061]
#define BIOSTRUC_ALIGN_slaves_E &at[1062]
#define BIOSTRUC_ALIGN_alignments &at[1063]
#define BIOSTRUC_ALIGN_sequences &at[1064]
#define BIOSTRUC_ALIGN_sequences_E &at[1065]
#define BIOSTRUC_ALIGN_seqalign &at[1066]
#define BIOSTRUC_ALIGN_seqalign_E &at[1067]

#define BIOSTRUC_ALIGN_SEQ &at[1070]
#define BIOSTRUC_ALIGN_SEQ_sequences &at[1071]
#define BIOSTRUC_ALIGN_SEQ_sequences_E &at[1072]
#define BIOSTRUC_ALIGN_SEQ_seqalign &at[1073]
#define BIOSTRUC_ALIGN_SEQ_seqalign_E &at[1074]

#define BIOSTRUC_SEQ &at[1076]
#define BIOSTRUC_SEQ_structure &at[1077]
#define BIOSTRUC_SEQ_sequences &at[1078]
#define BIOSTRUC_SEQ_sequences_E &at[1079]

#define ENTREZ_STYLE &at[1056]


/**************************************************
*
*    Defines for Module NCBI-ObjPrt
*
**************************************************/

#define PRINTTEMPLATE &at[1080]
#define PRINTTEMPLATE_name &at[1081]
#define PRINTTEMPLATE_labelfrom &at[1083]
#define PRINTTEMPLATE_format &at[1084]

#define PRINTTEMPLATESET &at[1114]
#define PRINTTEMPLATESET_E &at[1115]

#define TEMPLATENAME &at[1082]

#define PRINTFORMAT &at[1085]
#define PRINTFORMAT_asn1 &at[1086]
#define PRINTFORMAT_label &at[1087]
#define PRINTFORMAT_prefix &at[1088]
#define PRINTFORMAT_suffix &at[1089]
#define PRINTFORMAT_form &at[1090]

#define PRINTFORM &at[1091]
#define PRINTFORM_block &at[1092]
#define PRINTFORM_boolean &at[1097]
#define PRINTFORM_enum &at[1101]
#define PRINTFORM_text &at[1105]
#define PRINTFORM_use_template &at[1108]
#define PRINTFORM_user &at[1109]
#define PRINTFORM_null &at[1113]

#define PRINTFORMBLOCK &at[1093]
#define PRINTFORMBLOCK_separator &at[1094]
#define PRINTFORMBLOCK_components &at[1095]
#define PRINTFORMBLOCK_components_E &at[1096]

#define PRINTFORMBOOLEAN &at[1098]
#define PRINTFORMBOOLEAN_true &at[1099]
#define PRINTFORMBOOLEAN_false &at[1100]

#define PRINTFORMENUM &at[1102]
#define PRINTFORMENUM_values &at[1103]
#define PRINTFORMENUM_values_E &at[1104]

#define PRINTFORMTEXT &at[1106]
#define PRINTFORMTEXT_textfunc &at[1107]

#define USERFORMAT &at[1110]
#define USERFORMAT_printfunc &at[1111]
#define USERFORMAT_defaultfunc &at[1112]


/**************************************************
*
*    Defines for Module NCBI-Project
*
**************************************************/

#define PROJECT &at[1116]
#define PROJECT_descr &at[1117]
#define PROJECT_data &at[1132]

#define PROJECT_ITEM &at[1133]
#define PROJECT_ITEM_pmuid &at[1134]
#define PROJECT_ITEM_pmuid_E &at[1135]
#define PROJECT_ITEM_protuid &at[1136]
#define PROJECT_ITEM_protuid_E &at[1137]
#define PROJECT_ITEM_nucuid &at[1138]
#define PROJECT_ITEM_nucuid_E &at[1139]
#define PROJECT_ITEM_sequid &at[1140]
#define PROJECT_ITEM_sequid_E &at[1141]
#define PROJECT_ITEM_genomeuid &at[1142]
#define PROJECT_ITEM_genomeuid_E &at[1143]
#define PROJECT_ITEM_structuid &at[1144]
#define PROJECT_ITEM_structuid_E &at[1145]
#define PROJECT_ITEM_pmid &at[1146]
#define PROJECT_ITEM_pmid_E &at[1147]
#define PROJECT_ITEM_protid &at[1149]
#define PROJECT_ITEM_protid_E &at[1150]
#define PROJECT_ITEM_nucid &at[1152]
#define PROJECT_ITEM_nucid_E &at[1153]
#define PROJECT_ITEM_seqid &at[1154]
#define PROJECT_ITEM_seqid_E &at[1155]
#define PROJECT_ITEM_genomeid &at[1156]
#define PROJECT_ITEM_genomeid_E &at[1157]
#define PROJECT_ITEM_structid &at[1158]
#define PROJECT_ITEM_pment &at[1159]
#define PROJECT_ITEM_pment_E &at[1160]
#define PROJECT_ITEM_protent &at[1174]
#define PROJECT_ITEM_protent_E &at[1175]
#define PROJECT_ITEM_nucent &at[1177]
#define PROJECT_ITEM_nucent_E &at[1178]
#define PROJECT_ITEM_seqent &at[1179]
#define PROJECT_ITEM_seqent_E &at[1180]
#define PROJECT_ITEM_genomeent &at[1181]
#define PROJECT_ITEM_genomeent_E &at[1182]
#define PROJECT_ITEM_structent &at[1183]
#define PROJECT_ITEM_seqannot &at[1184]
#define PROJECT_ITEM_seqannot_E &at[1185]
#define PROJECT_ITEM_loc &at[1187]
#define PROJECT_ITEM_loc_E &at[1188]
#define PROJECT_ITEM_proj &at[1190]
#define PROJECT_ITEM_proj_E &at[1191]

#define PROJECT_DESCR &at[1118]
#define PROJECT_DESCR_id &at[1119]
#define PROJECT_DESCR_id_E &at[1120]
#define PROJECT_DESCR_name &at[1122]
#define PROJECT_DESCR_descr &at[1123]
#define PROJECT_DESCR_descr_E &at[1124]

#define PROJECT_ID &at[1121]

#define PROJDESC &at[1125]
#define PROJDESC_pub &at[1126]
#define PROJDESC_date &at[1128]
#define PROJDESC_comment &at[1130]
#define PROJDESC_title &at[1131]


/**************************************************
*
*    Defines for Module NCBI-Pub
*
**************************************************/

#define PUB &at[580]
#define PUB_gen &at[581]
#define PUB_sub &at[583]
#define PUB_medline &at[585]
#define PUB_muid &at[587]
#define PUB_article &at[588]
#define PUB_journal &at[590]
#define PUB_book &at[592]
#define PUB_proc &at[594]
#define PUB_patent &at[596]
#define PUB_pat_id &at[598]
#define PUB_man &at[600]
#define PUB_equiv &at[602]
#define PUB_pmid &at[603]

#define PUB_SET &at[909]
#define PUB_SET_pub &at[910]
#define PUB_SET_pub_E &at[911]
#define PUB_SET_medline &at[912]
#define PUB_SET_medline_E &at[913]
#define PUB_SET_article &at[914]
#define PUB_SET_article_E &at[915]
#define PUB_SET_journal &at[916]
#define PUB_SET_journal_E &at[917]
#define PUB_SET_book &at[918]
#define PUB_SET_book_E &at[919]
#define PUB_SET_proc &at[920]
#define PUB_SET_proc_E &at[921]
#define PUB_SET_patent &at[922]
#define PUB_SET_patent_E &at[923]

#define PUB_EQUIV &at[578]
#define PUB_EQUIV_E &at[579]


/**************************************************
*
*    Defines for Module NCBI-PubMed
*
**************************************************/

#define PUBMED_ENTRY &at[1162]
#define PUBMED_ENTRY_pmid &at[1163]
#define PUBMED_ENTRY_medent &at[1165]
#define PUBMED_ENTRY_publisher &at[1167]
#define PUBMED_ENTRY_urls &at[1168]
#define PUBMED_ENTRY_urls_E &at[1169]
#define PUBMED_ENTRY_pubid &at[1173]

#define PUBMED_URL &at[1170]
#define PUBMED_URL_location &at[1171]
#define PUBMED_URL_url &at[1172]


/**************************************************
*
*    Defines for Module NCBI-Sequence
*
**************************************************/

#define BIOSEQ &at[296]
#define BIOSEQ_id &at[297]
#define BIOSEQ_id_E &at[298]
#define BIOSEQ_descr &at[337]
#define BIOSEQ_inst &at[716]
#define BIOSEQ_annot &at[957]
#define BIOSEQ_annot_E &at[958]

#define SEQ_ANNOT &at[959]
#define SEQ_ANNOT_id &at[960]
#define SEQ_ANNOT_id_E &at[961]
#define SEQ_ANNOT_db &at[967]
#define SEQ_ANNOT_name &at[968]
#define SEQ_ANNOT_desc &at[969]
#define SEQ_ANNOT_data &at[987]
#define SEQ_ANNOT_data_ftable &at[988]
#define SEQ_ANNOT_data_ftable_E &at[989]
#define SEQ_ANNOT_data_align &at[990]
#define SEQ_ANNOT_data_align_E &at[991]
#define SEQ_ANNOT_data_graph &at[992]
#define SEQ_ANNOT_data_graph_E &at[993]
#define SEQ_ANNOT_data_ids &at[1027]
#define SEQ_ANNOT_data_ids_E &at[1028]
#define SEQ_ANNOT_data_locs &at[1029]
#define SEQ_ANNOT_data_locs_E &at[1030]

#define PUBDESC &at[575]
#define PUBDESC_pub &at[576]
#define PUBDESC_name &at[605]
#define PUBDESC_fig &at[606]
#define PUBDESC_num &at[607]
#define PUBDESC_numexc &at[608]
#define PUBDESC_poly_a &at[609]
#define PUBDESC_maploc &at[610]
#define PUBDESC_seq_raw &at[611]
#define PUBDESC_align_group &at[613]
#define PUBDESC_comment &at[614]
#define PUBDESC_reftype &at[615]

#define SEQ_DESCR &at[338]
#define SEQ_DESCR_E &at[339]

#define SEQDESC &at[340]
#define SEQDESC_mol_type &at[341]
#define SEQDESC_modif &at[343]
#define SEQDESC_modif_E &at[344]
#define SEQDESC_method &at[346]
#define SEQDESC_name &at[348]
#define SEQDESC_title &at[349]
#define SEQDESC_org &at[350]
#define SEQDESC_comment &at[393]
#define SEQDESC_num &at[394]
#define SEQDESC_maploc &at[539]
#define SEQDESC_pir &at[541]
#define SEQDESC_genbank &at[560]
#define SEQDESC_pub &at[574]
#define SEQDESC_region &at[616]
#define SEQDESC_user &at[617]
#define SEQDESC_sp &at[619]
#define SEQDESC_dbxref &at[640]
#define SEQDESC_embl &at[641]
#define SEQDESC_create_date &at[663]
#define SEQDESC_update_date &at[665]
#define SEQDESC_prf &at[666]
#define SEQDESC_pdb &at[678]
#define SEQDESC_het &at[694]
#define SEQDESC_source &at[696]
#define SEQDESC_molinfo &at[710]

#define NUMBERING &at[395]
#define NUMBERING_cont &at[396]
#define NUMBERING_enum &at[401]
#define NUMBERING_ref &at[406]
#define NUMBERING_real &at[534]

#define HETEROGEN &at[695]

#define SEQ_HIST &at[945]
#define SEQ_HIST_assembly &at[946]
#define SEQ_HIST_assembly_E &at[947]
#define SEQ_HIST_replaces &at[948]
#define SEQ_HIST_replaced_by &at[953]
#define SEQ_HIST_deleted &at[954]
#define SEQ_HIST_deleted_bool &at[955]
#define SEQ_HIST_deleted_date &at[956]

#define SEQ_INST &at[717]
#define SEQ_INST_repr &at[718]
#define SEQ_INST_mol &at[719]
#define SEQ_INST_length &at[720]
#define SEQ_INST_fuzz &at[721]
#define SEQ_INST_topology &at[723]
#define SEQ_INST_strand &at[724]
#define SEQ_INST_seq_data &at[725]
#define SEQ_INST_ext &at[747]
#define SEQ_INST_hist &at[944]

#define GIBB_MOL &at[342]

#define GIBB_MOD &at[345]

#define GIBB_METHOD &at[347]

#define MOLINFO &at[711]
#define MOLINFO_biomol &at[712]
#define MOLINFO_tech &at[713]
#define MOLINFO_techexp &at[714]
#define MOLINFO_completeness &at[715]

#define NUM_CONT &at[397]
#define NUM_CONT_refnum &at[398]
#define NUM_CONT_has_zero &at[399]
#define NUM_CONT_ascending &at[400]

#define NUM_ENUM &at[402]
#define NUM_ENUM_num &at[403]
#define NUM_ENUM_names &at[404]
#define NUM_ENUM_names_E &at[405]

#define NUM_REF &at[407]
#define NUM_REF_type &at[408]
#define NUM_REF_aligns &at[409]

#define NUM_REAL &at[535]
#define NUM_REAL_a &at[536]
#define NUM_REAL_b &at[537]
#define NUM_REAL_units &at[538]

#define SEQ_DATA &at[726]
#define SEQ_DATA_iupacna &at[727]
#define SEQ_DATA_iupacaa &at[729]
#define SEQ_DATA_ncbi2na &at[731]
#define SEQ_DATA_ncbi4na &at[733]
#define SEQ_DATA_ncbi8na &at[735]
#define SEQ_DATA_ncbipna &at[737]
#define SEQ_DATA_ncbi8aa &at[739]
#define SEQ_DATA_ncbieaa &at[741]
#define SEQ_DATA_ncbipaa &at[743]
#define SEQ_DATA_ncbistdaa &at[745]

#define SEQ_EXT &at[748]
#define SEQ_EXT_seg &at[749]
#define SEQ_EXT_ref &at[753]
#define SEQ_EXT_map &at[755]
#define SEQ_EXT_delta &at[934]

#define SEG_EXT &at[750]
#define SEG_EXT_E &at[751]

#define REF_EXT &at[754]

#define MAP_EXT &at[756]
#define MAP_EXT_E &at[757]

#define DELTA_EXT &at[935]
#define DELTA_EXT_E &at[936]

#define DELTA_SEQ &at[937]
#define DELTA_SEQ_loc &at[938]
#define DELTA_SEQ_literal &at[939]

#define SEQ_LITERAL &at[940]
#define SEQ_LITERAL_length &at[941]
#define SEQ_LITERAL_fuzz &at[942]
#define SEQ_LITERAL_seq_data &at[943]

#define SEQ_HIST_REC &at[949]
#define SEQ_HIST_REC_date &at[950]
#define SEQ_HIST_REC_ids &at[951]
#define SEQ_HIST_REC_ids_E &at[952]

#define IUPACNA &at[728]

#define IUPACAA &at[730]

#define NCBI2NA &at[732]

#define NCBI4NA &at[734]

#define NCBI8NA &at[736]

#define NCBIPNA &at[738]

#define NCBI8AA &at[740]

#define NCBIEAA &at[742]

#define NCBIPAA &at[744]

#define NCBISTDAA &at[746]

#define ANNOT_ID &at[962]
#define ANNOT_ID_local &at[963]
#define ANNOT_ID_ncbi &at[965]
#define ANNOT_ID_general &at[966]

#define ANNOT_DESCR &at[970]
#define ANNOT_DESCR_E &at[971]

#define ANNOTDESC &at[972]
#define ANNOTDESC_name &at[973]
#define ANNOTDESC_title &at[974]
#define ANNOTDESC_comment &at[975]
#define ANNOTDESC_pub &at[976]
#define ANNOTDESC_user &at[977]
#define ANNOTDESC_create_date &at[978]
#define ANNOTDESC_update_date &at[979]
#define ANNOTDESC_src &at[980]
#define ANNOTDESC_align &at[981]
#define ANNOTDESC_region &at[986]

#define ALIGN_DEF &at[982]
#define ALIGN_DEF_align_type &at[983]
#define ALIGN_DEF_ids &at[984]
#define ALIGN_DEF_ids_E &at[985]


/**************************************************
*
*    Defines for Module NCBI-Seqalign
*
**************************************************/

#define SEQ_ALIGN &at[411]
#define SEQ_ALIGN_type &at[412]
#define SEQ_ALIGN_dim &at[413]
#define SEQ_ALIGN_score &at[414]
#define SEQ_ALIGN_score_E &at[415]
#define SEQ_ALIGN_segs &at[422]
#define SEQ_ALIGN_segs_dendiag &at[423]
#define SEQ_ALIGN_segs_dendiag_E &at[424]
#define SEQ_ALIGN_segs_denseg &at[439]
#define SEQ_ALIGN_segs_std &at[453]
#define SEQ_ALIGN_segs_std_E &at[454]
#define SEQ_ALIGN_segs_packed &at[514]
#define SEQ_ALIGN_segs_disc &at[529]
#define SEQ_ALIGN_bounds &at[532]
#define SEQ_ALIGN_bounds_E &at[533]

#define SCORE &at[416]
#define SCORE_id &at[417]
#define SCORE_value &at[419]
#define SCORE_value_real &at[420]
#define SCORE_value_int &at[421]

#define SCORE_SET &at[1192]
#define SCORE_SET_E &at[1193]

#define SEQ_ALIGN_SET &at[530]
#define SEQ_ALIGN_SET_E &at[531]

#define DENSE_DIAG &at[425]
#define DENSE_DIAG_dim &at[426]
#define DENSE_DIAG_ids &at[427]
#define DENSE_DIAG_ids_E &at[428]
#define DENSE_DIAG_starts &at[430]
#define DENSE_DIAG_starts_E &at[431]
#define DENSE_DIAG_len &at[432]
#define DENSE_DIAG_strands &at[433]
#define DENSE_DIAG_strands_E &at[434]
#define DENSE_DIAG_scores &at[437]
#define DENSE_DIAG_scores_E &at[438]

#define DENSE_SEG &at[440]
#define DENSE_SEG_dim &at[441]
#define DENSE_SEG_numseg &at[442]
#define DENSE_SEG_ids &at[443]
#define DENSE_SEG_ids_E &at[444]
#define DENSE_SEG_starts &at[445]
#define DENSE_SEG_starts_E &at[446]
#define DENSE_SEG_lens &at[447]
#define DENSE_SEG_lens_E &at[448]
#define DENSE_SEG_strands &at[449]
#define DENSE_SEG_strands_E &at[450]
#define DENSE_SEG_scores &at[451]
#define DENSE_SEG_scores_E &at[452]

#define STD_SEG &at[455]
#define STD_SEG_dim &at[456]
#define STD_SEG_ids &at[457]
#define STD_SEG_ids_E &at[458]
#define STD_SEG_loc &at[459]
#define STD_SEG_loc_E &at[460]
#define STD_SEG_scores &at[512]
#define STD_SEG_scores_E &at[513]

#define PACKED_SEG &at[515]
#define PACKED_SEG_dim &at[516]
#define PACKED_SEG_numseg &at[517]
#define PACKED_SEG_ids &at[518]
#define PACKED_SEG_ids_E &at[519]
#define PACKED_SEG_starts &at[520]
#define PACKED_SEG_starts_E &at[521]
#define PACKED_SEG_present &at[522]
#define PACKED_SEG_lens &at[523]
#define PACKED_SEG_lens_E &at[524]
#define PACKED_SEG_strands &at[525]
#define PACKED_SEG_strands_E &at[526]
#define PACKED_SEG_scores &at[527]
#define PACKED_SEG_scores_E &at[528]


/**************************************************
*
*    Defines for Module EMBL-General
*
**************************************************/

#define EMBL_DBNAME &at[657]
#define EMBL_DBNAME_code &at[658]
#define EMBL_DBNAME_name &at[659]

#define EMBL_XREF &at[655]
#define EMBL_XREF_dbname &at[656]
#define EMBL_XREF_id &at[660]
#define EMBL_XREF_id_E &at[661]

#define EMBL_BLOCK &at[643]
#define EMBL_BLOCK_class &at[644]
#define EMBL_BLOCK_div &at[645]
#define EMBL_BLOCK_creation_date &at[646]
#define EMBL_BLOCK_update_date &at[648]
#define EMBL_BLOCK_extra_acc &at[649]
#define EMBL_BLOCK_extra_acc_E &at[650]
#define EMBL_BLOCK_keywords &at[651]
#define EMBL_BLOCK_keywords_E &at[652]
#define EMBL_BLOCK_xref &at[653]
#define EMBL_BLOCK_xref_E &at[654]


/**************************************************
*
*    Defines for Module SP-General
*
**************************************************/

#define SP_BLOCK &at[621]
#define SP_BLOCK_class &at[622]
#define SP_BLOCK_extra_acc &at[623]
#define SP_BLOCK_extra_acc_E &at[624]
#define SP_BLOCK_imeth &at[625]
#define SP_BLOCK_plasnm &at[626]
#define SP_BLOCK_plasnm_E &at[627]
#define SP_BLOCK_seqref &at[628]
#define SP_BLOCK_seqref_E &at[629]
#define SP_BLOCK_dbref &at[631]
#define SP_BLOCK_dbref_E &at[632]
#define SP_BLOCK_keywords &at[634]
#define SP_BLOCK_keywords_E &at[635]
#define SP_BLOCK_created &at[636]
#define SP_BLOCK_sequpd &at[638]
#define SP_BLOCK_annotupd &at[639]


/**************************************************
*
*    Defines for Module PIR-General
*
**************************************************/

#define PIR_BLOCK &at[543]
#define PIR_BLOCK_had_punct &at[544]
#define PIR_BLOCK_host &at[545]
#define PIR_BLOCK_source &at[546]
#define PIR_BLOCK_summary &at[547]
#define PIR_BLOCK_genetic &at[548]
#define PIR_BLOCK_includes &at[549]
#define PIR_BLOCK_placement &at[550]
#define PIR_BLOCK_superfamily &at[551]
#define PIR_BLOCK_keywords &at[552]
#define PIR_BLOCK_keywords_E &at[553]
#define PIR_BLOCK_cross_reference &at[554]
#define PIR_BLOCK_date &at[555]
#define PIR_BLOCK_seq_raw &at[556]
#define PIR_BLOCK_seqref &at[557]
#define PIR_BLOCK_seqref_E &at[558]


/**************************************************
*
*    Defines for Module GenBank-General
*
**************************************************/

#define GB_BLOCK &at[562]
#define GB_BLOCK_extra_accessions &at[563]
#define GB_BLOCK_extra_accessions_E &at[564]
#define GB_BLOCK_source &at[565]
#define GB_BLOCK_keywords &at[566]
#define GB_BLOCK_keywords_E &at[567]
#define GB_BLOCK_origin &at[568]
#define GB_BLOCK_date &at[569]
#define GB_BLOCK_entry_date &at[570]
#define GB_BLOCK_div &at[572]
#define GB_BLOCK_taxonomy &at[573]


/**************************************************
*
*    Defines for Module PRF-General
*
**************************************************/

#define PRF_BLOCK &at[668]
#define PRF_BLOCK_extra_src &at[669]
#define PRF_BLOCK_keywords &at[676]
#define PRF_BLOCK_keywords_E &at[677]

#define PRF_EXTRASRC &at[670]
#define PRF_EXTRASRC_host &at[671]
#define PRF_EXTRASRC_part &at[672]
#define PRF_EXTRASRC_state &at[673]
#define PRF_EXTRASRC_strain &at[674]
#define PRF_EXTRASRC_taxon &at[675]


/**************************************************
*
*    Defines for Module PDB-General
*
**************************************************/

#define PDB_BLOCK &at[680]
#define PDB_BLOCK_deposition &at[681]
#define PDB_BLOCK_class &at[683]
#define PDB_BLOCK_compound &at[684]
#define PDB_BLOCK_compound_E &at[685]
#define PDB_BLOCK_source &at[686]
#define PDB_BLOCK_source_E &at[687]
#define PDB_BLOCK_exp_method &at[688]
#define PDB_BLOCK_replace &at[689]

#define PDB_REPLACE &at[690]
#define PDB_REPLACE_date &at[691]
#define PDB_REPLACE_ids &at[692]
#define PDB_REPLACE_ids_E &at[693]


/**************************************************
*
*    Defines for Module NCBI-SeqCode
*
**************************************************/

#define SEQ_CODE_TABLE &at[1194]
#define SEQ_CODE_TABLE_code &at[1195]
#define SEQ_CODE_TABLE_num &at[1197]
#define SEQ_CODE_TABLE_one_letter &at[1198]
#define SEQ_CODE_TABLE_start_at &at[1199]
#define SEQ_CODE_TABLE_table &at[1200]
#define SEQ_CODE_TABLE_table_E &at[1201]
#define SEQ_CODE_TABLE_table_E_symbol &at[1202]
#define SEQ_CODE_TABLE_table_E_name &at[1203]
#define SEQ_CODE_TABLE_comps &at[1204]
#define SEQ_CODE_TABLE_comps_E &at[1205]

#define SEQ_MAP_TABLE &at[1206]
#define SEQ_MAP_TABLE_from &at[1207]
#define SEQ_MAP_TABLE_to &at[1208]
#define SEQ_MAP_TABLE_num &at[1209]
#define SEQ_MAP_TABLE_start_at &at[1210]
#define SEQ_MAP_TABLE_table &at[1211]
#define SEQ_MAP_TABLE_table_E &at[1212]

#define SEQ_CODE_SET &at[1213]
#define SEQ_CODE_SET_codes &at[1214]
#define SEQ_CODE_SET_codes_E &at[1215]
#define SEQ_CODE_SET_maps &at[1216]
#define SEQ_CODE_SET_maps_E &at[1217]

#define SEQ_CODE_TYPE &at[1196]


/**************************************************
*
*    Defines for Module NCBI-Seqfeat
*
**************************************************/

#define SEQ_FEAT &at[759]
#define SEQ_FEAT_id &at[760]
#define SEQ_FEAT_data &at[761]
#define SEQ_FEAT_partial &at[895]
#define SEQ_FEAT_except &at[896]
#define SEQ_FEAT_comment &at[897]
#define SEQ_FEAT_product &at[898]
#define SEQ_FEAT_location &at[899]
#define SEQ_FEAT_qual &at[900]
#define SEQ_FEAT_qual_E &at[901]
#define SEQ_FEAT_title &at[905]
#define SEQ_FEAT_ext &at[906]
#define SEQ_FEAT_cit &at[907]
#define SEQ_FEAT_exp_ev &at[924]
#define SEQ_FEAT_xref &at[925]
#define SEQ_FEAT_xref_E &at[926]
#define SEQ_FEAT_dbxref &at[930]
#define SEQ_FEAT_dbxref_E &at[931]
#define SEQ_FEAT_pseudo &at[932]
#define SEQ_FEAT_except_text &at[933]

#define FEAT_ID &at[504]
#define FEAT_ID_gibb &at[505]
#define FEAT_ID_giim &at[506]
#define FEAT_ID_local &at[508]
#define FEAT_ID_general &at[510]

#define GENETIC_CODE &at[786]
#define GENETIC_CODE_E &at[787]
#define GENETIC_CODE_E_name &at[788]
#define GENETIC_CODE_E_id &at[789]
#define GENETIC_CODE_E_ncbieaa &at[790]
#define GENETIC_CODE_E_ncbi8aa &at[791]
#define GENETIC_CODE_E_ncbistdaa &at[792]
#define GENETIC_CODE_E_sncbieaa &at[793]
#define GENETIC_CODE_E_sncbi8aa &at[794]
#define GENETIC_CODE_E_sncbistdaa &at[795]

#define SEQFEATDATA &at[762]
#define SEQFEATDATA_gene &at[763]
#define SEQFEATDATA_org &at[776]
#define SEQFEATDATA_cdregion &at[778]
#define SEQFEATDATA_prot &at[806]
#define SEQFEATDATA_rna &at[820]
#define SEQFEATDATA_pub &at[838]
#define SEQFEATDATA_seq &at[840]
#define SEQFEATDATA_imp &at[841]
#define SEQFEATDATA_region &at[846]
#define SEQFEATDATA_comment &at[847]
#define SEQFEATDATA_bond &at[848]
#define SEQFEATDATA_site &at[849]
#define SEQFEATDATA_rsite &at[850]
#define SEQFEATDATA_user &at[856]
#define SEQFEATDATA_txinit &at[858]
#define SEQFEATDATA_num &at[887]
#define SEQFEATDATA_psec_str &at[889]
#define SEQFEATDATA_non_std_residue &at[890]
#define SEQFEATDATA_het &at[891]
#define SEQFEATDATA_biosrc &at[893]

#define GB_QUAL &at[902]
#define GB_QUAL_qual &at[903]
#define GB_QUAL_val &at[904]

#define SEQFEATXREF &at[927]
#define SEQFEATXREF_id &at[928]
#define SEQFEATXREF_data &at[929]

#define CDREGION &at[779]
#define CDREGION_orf &at[780]
#define CDREGION_frame &at[781]
#define CDREGION_conflict &at[782]
#define CDREGION_gaps &at[783]
#define CDREGION_mismatch &at[784]
#define CDREGION_code &at[785]
#define CDREGION_code_break &at[796]
#define CDREGION_code_break_E &at[797]
#define CDREGION_stops &at[805]

#define IMP_FEAT &at[842]
#define IMP_FEAT_key &at[843]
#define IMP_FEAT_loc &at[844]
#define IMP_FEAT_descr &at[845]

#define CODE_BREAK &at[798]
#define CODE_BREAK_loc &at[799]
#define CODE_BREAK_aa &at[801]
#define CODE_BREAK_aa_ncbieaa &at[802]
#define CODE_BREAK_aa_ncbi8aa &at[803]
#define CODE_BREAK_aa_ncbistdaa &at[804]

#define GENETIC_CODE_TABLE &at[1218]
#define GENETIC_CODE_TABLE_E &at[1219]


/**************************************************
*
*    Defines for Module NCBI-Rsite
*
**************************************************/

#define RSITE_REF &at[852]
#define RSITE_REF_str &at[853]
#define RSITE_REF_db &at[854]


/**************************************************
*
*    Defines for Module NCBI-RNA
*
**************************************************/

#define RNA_REF &at[822]
#define RNA_REF_type &at[823]
#define RNA_REF_pseudo &at[824]
#define RNA_REF_ext &at[825]
#define RNA_REF_ext_name &at[826]
#define RNA_REF_ext_tRNA &at[827]

#define TRNA_EXT &at[828]
#define TRNA_EXT_aa &at[829]
#define TRNA_EXT_aa_iupacaa &at[830]
#define TRNA_EXT_aa_ncbieaa &at[831]
#define TRNA_EXT_aa_ncbi8aa &at[832]
#define TRNA_EXT_aa_ncbistdaa &at[833]
#define TRNA_EXT_codon &at[834]
#define TRNA_EXT_codon_E &at[835]
#define TRNA_EXT_anticodon &at[836]


/**************************************************
*
*    Defines for Module NCBI-Gene
*
**************************************************/

#define GENE_REF &at[765]
#define GENE_REF_locus &at[766]
#define GENE_REF_allele &at[767]
#define GENE_REF_desc &at[768]
#define GENE_REF_maploc &at[769]
#define GENE_REF_pseudo &at[770]
#define GENE_REF_db &at[771]
#define GENE_REF_db_E &at[772]
#define GENE_REF_syn &at[774]
#define GENE_REF_syn_E &at[775]


/**************************************************
*
*    Defines for Module NCBI-Organism
*
**************************************************/

#define ORG_REF &at[352]
#define ORG_REF_taxname &at[353]
#define ORG_REF_common &at[354]
#define ORG_REF_mod &at[355]
#define ORG_REF_mod_E &at[356]
#define ORG_REF_db &at[357]
#define ORG_REF_db_E &at[358]
#define ORG_REF_syn &at[360]
#define ORG_REF_syn_E &at[361]
#define ORG_REF_orgname &at[362]

#define ORGNAME &at[363]
#define ORGNAME_name &at[364]
#define ORGNAME_name_binomial &at[365]
#define ORGNAME_name_virus &at[370]
#define ORGNAME_name_hybrid &at[371]
#define ORGNAME_name_namedhybrid &at[374]
#define ORGNAME_name_partial &at[375]
#define ORGNAME_attrib &at[382]
#define ORGNAME_mod &at[383]
#define ORGNAME_mod_E &at[384]
#define ORGNAME_lineage &at[389]
#define ORGNAME_gcode &at[390]
#define ORGNAME_mgcode &at[391]
#define ORGNAME_div &at[392]

#define BINOMIALORGNAME &at[366]
#define BINOMIALORGNAME_genus &at[367]
#define BINOMIALORGNAME_species &at[368]
#define BINOMIALORGNAME_subspecies &at[369]

#define MULTIORGNAME &at[372]
#define MULTIORGNAME_E &at[373]

#define PARTIALORGNAME &at[376]
#define PARTIALORGNAME_E &at[377]

#define ORGMOD &at[385]
#define ORGMOD_subtype &at[386]
#define ORGMOD_subname &at[387]
#define ORGMOD_attrib &at[388]

#define TAXELEMENT &at[378]
#define TAXELEMENT_fixed_level &at[379]
#define TAXELEMENT_level &at[380]
#define TAXELEMENT_name &at[381]


/**************************************************
*
*    Defines for Module NCBI-BioSource
*
**************************************************/

#define BIOSOURCE &at[698]
#define BIOSOURCE_genome &at[699]
#define BIOSOURCE_origin &at[700]
#define BIOSOURCE_org &at[701]
#define BIOSOURCE_subtype &at[703]
#define BIOSOURCE_subtype_E &at[704]
#define BIOSOURCE_is_focus &at[709]

#define SUBSOURCE &at[705]
#define SUBSOURCE_subtype &at[706]
#define SUBSOURCE_name &at[707]
#define SUBSOURCE_attrib &at[708]


/**************************************************
*
*    Defines for Module NCBI-Protein
*
**************************************************/

#define PROT_REF &at[808]
#define PROT_REF_name &at[809]
#define PROT_REF_name_E &at[810]
#define PROT_REF_desc &at[811]
#define PROT_REF_ec &at[812]
#define PROT_REF_ec_E &at[813]
#define PROT_REF_activity &at[814]
#define PROT_REF_activity_E &at[815]
#define PROT_REF_db &at[816]
#define PROT_REF_db_E &at[817]
#define PROT_REF_processed &at[819]


/**************************************************
*
*    Defines for Module NCBI-TxInit
*
**************************************************/

#define TXINIT &at[860]
#define TXINIT_name &at[861]
#define TXINIT_syn &at[862]
#define TXINIT_syn_E &at[863]
#define TXINIT_gene &at[864]
#define TXINIT_gene_E &at[865]
#define TXINIT_protein &at[867]
#define TXINIT_protein_E &at[868]
#define TXINIT_rna &at[870]
#define TXINIT_rna_E &at[871]
#define TXINIT_expression &at[872]
#define TXINIT_txsystem &at[873]
#define TXINIT_txdescr &at[874]
#define TXINIT_txorg &at[875]
#define TXINIT_mapping_precise &at[877]
#define TXINIT_location_accurate &at[878]
#define TXINIT_inittype &at[879]
#define TXINIT_evidence &at[880]
#define TXINIT_evidence_E &at[881]

#define TX_EVIDENCE &at[882]
#define TX_EVIDENCE_exp_code &at[883]
#define TX_EVIDENCE_expression_system &at[884]
#define TX_EVIDENCE_low_prec_data &at[885]
#define TX_EVIDENCE_from_homolog &at[886]


/**************************************************
*
*    Defines for Module NCBI-Seqloc
*
**************************************************/

#define SEQ_ID &at[300]
#define SEQ_ID_local &at[301]
#define SEQ_ID_gibbsq &at[303]
#define SEQ_ID_gibbmt &at[304]
#define SEQ_ID_giim &at[305]
#define SEQ_ID_genbank &at[310]
#define SEQ_ID_embl &at[316]
#define SEQ_ID_pir &at[317]
#define SEQ_ID_swissprot &at[318]
#define SEQ_ID_patent &at[319]
#define SEQ_ID_other &at[324]
#define SEQ_ID_general &at[325]
#define SEQ_ID_gi &at[327]
#define SEQ_ID_ddbj &at[328]
#define SEQ_ID_prf &at[329]
#define SEQ_ID_pdb &at[330]

#define SEQ_LOC &at[462]
#define SEQ_LOC_null &at[463]
#define SEQ_LOC_empty &at[465]
#define SEQ_LOC_whole &at[466]
#define SEQ_LOC_int &at[467]
#define SEQ_LOC_packed_int &at[476]
#define SEQ_LOC_pnt &at[479]
#define SEQ_LOC_packed_pnt &at[485]
#define SEQ_LOC_mix &at[492]
#define SEQ_LOC_equiv &at[495]
#define SEQ_LOC_bond &at[498]
#define SEQ_LOC_feat &at[502]

#define SEQ_INTERVAL &at[468]
#define SEQ_INTERVAL_from &at[469]
#define SEQ_INTERVAL_to &at[470]
#define SEQ_INTERVAL_strand &at[471]
#define SEQ_INTERVAL_id &at[472]
#define SEQ_INTERVAL_fuzz_from &at[473]
#define SEQ_INTERVAL_fuzz_to &at[475]

#define PACKED_SEQINT &at[477]
#define PACKED_SEQINT_E &at[478]

#define SEQ_POINT &at[480]
#define SEQ_POINT_point &at[481]
#define SEQ_POINT_strand &at[482]
#define SEQ_POINT_id &at[483]
#define SEQ_POINT_fuzz &at[484]

#define PACKED_SEQPNT &at[486]
#define PACKED_SEQPNT_strand &at[487]
#define PACKED_SEQPNT_id &at[488]
#define PACKED_SEQPNT_fuzz &at[489]
#define PACKED_SEQPNT_points &at[490]
#define PACKED_SEQPNT_points_E &at[491]

#define NA_STRAND &at[436]

#define GIIMPORT_ID &at[306]
#define GIIMPORT_ID_id &at[307]
#define GIIMPORT_ID_db &at[308]
#define GIIMPORT_ID_release &at[309]

#define TEXTSEQ_ID &at[311]
#define TEXTSEQ_ID_name &at[312]
#define TEXTSEQ_ID_accession &at[313]
#define TEXTSEQ_ID_release &at[314]
#define TEXTSEQ_ID_version &at[315]

#define PATENT_SEQ_ID &at[320]
#define PATENT_SEQ_ID_seqid &at[321]
#define PATENT_SEQ_ID_cit &at[322]

#define PDB_SEQ_ID &at[331]
#define PDB_SEQ_ID_mol &at[332]
#define PDB_SEQ_ID_chain &at[334]
#define PDB_SEQ_ID_rel &at[335]

#define PDB_MOL_ID &at[333]

#define SEQ_LOC_MIX &at[493]
#define SEQ_LOC_MIX_E &at[494]

#define SEQ_LOC_EQUIV &at[496]
#define SEQ_LOC_EQUIV_E &at[497]

#define SEQ_BOND &at[499]
#define SEQ_BOND_a &at[500]
#define SEQ_BOND_b &at[501]


/**************************************************
*
*    Defines for Module NCBI-Seqres
*
**************************************************/

#define SEQ_GRAPH &at[995]
#define SEQ_GRAPH_title &at[996]
#define SEQ_GRAPH_comment &at[997]
#define SEQ_GRAPH_loc &at[998]
#define SEQ_GRAPH_title_x &at[1000]
#define SEQ_GRAPH_title_y &at[1001]
#define SEQ_GRAPH_comp &at[1002]
#define SEQ_GRAPH_a &at[1003]
#define SEQ_GRAPH_b &at[1004]
#define SEQ_GRAPH_numval &at[1005]
#define SEQ_GRAPH_graph &at[1006]
#define SEQ_GRAPH_graph_real &at[1007]
#define SEQ_GRAPH_graph_int &at[1014]
#define SEQ_GRAPH_graph_byte &at[1021]

#define REAL_GRAPH &at[1008]
#define REAL_GRAPH_max &at[1009]
#define REAL_GRAPH_min &at[1010]
#define REAL_GRAPH_axis &at[1011]
#define REAL_GRAPH_values &at[1012]
#define REAL_GRAPH_values_E &at[1013]

#define INT_GRAPH &at[1015]
#define INT_GRAPH_max &at[1016]
#define INT_GRAPH_min &at[1017]
#define INT_GRAPH_axis &at[1018]
#define INT_GRAPH_values &at[1019]
#define INT_GRAPH_values_E &at[1020]

#define BYTE_GRAPH &at[1022]
#define BYTE_GRAPH_max &at[1023]
#define BYTE_GRAPH_min &at[1024]
#define BYTE_GRAPH_axis &at[1025]
#define BYTE_GRAPH_values &at[1026]


/**************************************************
*
*    Defines for Module NCBI-Seqset
*
**************************************************/

#define BIOSEQ_SET &at[1032]
#define BIOSEQ_SET_id &at[1033]
#define BIOSEQ_SET_coll &at[1035]
#define BIOSEQ_SET_level &at[1037]
#define BIOSEQ_SET_class &at[1038]
#define BIOSEQ_SET_release &at[1039]
#define BIOSEQ_SET_date &at[1040]
#define BIOSEQ_SET_descr &at[1042]
#define BIOSEQ_SET_seq_set &at[1044]
#define BIOSEQ_SET_seq_set_E &at[1045]
#define BIOSEQ_SET_annot &at[1046]
#define BIOSEQ_SET_annot_E &at[1047]

#define SEQ_ENTRY &at[293]
#define SEQ_ENTRY_seq &at[294]
#define SEQ_ENTRY_set &at[1031]


/**************************************************
*
*    Defines for Module NCBI-Submit
*
**************************************************/

#define SEQ_SUBMIT &at[1220]
#define SEQ_SUBMIT_sub &at[1221]
#define SEQ_SUBMIT_data &at[1249]
#define SEQ_SUBMIT_data_entrys &at[1250]
#define SEQ_SUBMIT_data_entrys_E &at[1251]
#define SEQ_SUBMIT_data_annots &at[1253]
#define SEQ_SUBMIT_data_annots_E &at[1254]
#define SEQ_SUBMIT_data_delete &at[1256]
#define SEQ_SUBMIT_data_delete_E &at[1257]

#define SUBMIT_BLOCK &at[1222]
#define SUBMIT_BLOCK_contact &at[1223]
#define SUBMIT_BLOCK_cit &at[1240]
#define SUBMIT_BLOCK_hup &at[1242]
#define SUBMIT_BLOCK_reldate &at[1243]
#define SUBMIT_BLOCK_subtype &at[1245]
#define SUBMIT_BLOCK_tool &at[1246]
#define SUBMIT_BLOCK_user_tag &at[1247]
#define SUBMIT_BLOCK_comment &at[1248]

#define CONTACT_INFO &at[1224]
#define CONTACT_INFO_name &at[1225]
#define CONTACT_INFO_address &at[1226]
#define CONTACT_INFO_address_E &at[1227]
#define CONTACT_INFO_phone &at[1228]
#define CONTACT_INFO_fax &at[1229]
#define CONTACT_INFO_email &at[1230]
#define CONTACT_INFO_telex &at[1231]
#define CONTACT_INFO_owner_id &at[1232]
#define CONTACT_INFO_password &at[1234]
#define CONTACT_INFO_last_name &at[1235]
#define CONTACT_INFO_first_name &at[1236]
#define CONTACT_INFO_middle_initial &at[1237]
#define CONTACT_INFO_contact &at[1238]
