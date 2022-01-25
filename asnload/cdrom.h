/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "cdrom.l60";
static AsnValxNodePtr avn;
static AsnTypePtr at;
static AsnModulePtr amp;


/**************************************************
*
*    Defines for Module NCBI-CdRom
*
**************************************************/

#define CDROM_INF &at[0]
#define CDROM_INF_volume_label &at[1]
#define CDROM_INF_version &at[3]
#define CDROM_INF_issue &at[5]
#define CDROM_INF_format &at[6]
#define CDROM_INF_descr &at[7]
#define CDROM_INF_no_compression &at[8]
#define CDROM_INF_huff_count &at[10]
#define CDROM_INF_huff_left &at[11]
#define CDROM_INF_huff_left_E &at[12]
#define CDROM_INF_huff_right &at[14]
#define CDROM_INF_huff_right_E &at[15]
#define CDROM_INF_type_count &at[16]
#define CDROM_INF_type_names &at[17]
#define CDROM_INF_type_names_E &at[18]
#define CDROM_INF_type_bucket_size &at[19]
#define CDROM_INF_field_count &at[20]
#define CDROM_INF_field_names &at[21]
#define CDROM_INF_field_names_E &at[22]
#define CDROM_INF_field_bucket_size &at[23]
#define CDROM_INF_types &at[24]
#define CDROM_INF_types_E &at[25]
#define CDROM_INF_release_date &at[38]
#define CDROM_INF_close_date &at[43]
#define CDROM_INF_type_info &at[44]
#define CDROM_INF_type_info_E &at[45]
#define CDROM_INF_field_info &at[52]
#define CDROM_INF_field_info_E &at[53]
#define CDROM_INF_div_count &at[66]
#define CDROM_INF_div_info &at[67]
#define CDROM_INF_div_info_E &at[68]
#define CDROM_INF_link_count &at[76]
#define CDROM_INF_link_info &at[77]
#define CDROM_INF_link_info_E &at[78]

#define DOCSUM &at[88]
#define DOCSUM_no_abstract &at[89]
#define DOCSUM_translated_title &at[90]
#define DOCSUM_no_authors &at[91]
#define DOCSUM_caption &at[92]
#define DOCSUM_title &at[93]
#define DOCSUM_extra &at[94]
#define DOCSUM_non_document &at[95]
#define DOCSUM_is_segmented &at[96]
#define DOCSUM_is_partial &at[97]
#define DOCSUM_create &at[98]
#define DOCSUM_modify &at[99]
#define DOCSUM_link_count &at[100]
#define DOCSUM_link_count_E &at[101]
#define DOCSUM_uid &at[102]
#define DOCSUM_secondaryUid &at[103]
#define DOCSUM_not_yet_neighbored &at[104]

#define TYPEDATA &at[26]
#define TYPEDATA_num &at[27]
#define TYPEDATA_num_uids &at[28]
#define TYPEDATA_minuid &at[29]
#define TYPEDATA_maxuid &at[30]
#define TYPEDATA_num_bucket &at[31]
#define TYPEDATA_fields &at[32]
#define TYPEDATA_fields_E &at[33]

#define CD_DATE &at[39]
#define CD_DATE_year &at[40]
#define CD_DATE_month &at[41]
#define CD_DATE_day &at[42]

#define TYPE_INFO &at[46]
#define TYPE_INFO_id &at[47]
#define TYPE_INFO_tag &at[48]
#define TYPE_INFO_name &at[49]
#define TYPE_INFO_descr &at[50]
#define TYPE_INFO_asntype &at[51]

#define FIELD_INFO &at[54]
#define FIELD_INFO_id &at[55]
#define FIELD_INFO_tag &at[56]
#define FIELD_INFO_name &at[57]
#define FIELD_INFO_descr &at[58]
#define FIELD_INFO_single_token &at[59]
#define FIELD_INFO_has_special &at[60]
#define FIELD_INFO_hier_avail &at[61]
#define FIELD_INFO_hier_id &at[62]
#define FIELD_INFO_post_type &at[63]

#define DIV_INFO &at[69]
#define DIV_INFO_tag &at[70]
#define DIV_INFO_descr &at[71]
#define DIV_INFO_reldate &at[72]
#define DIV_INFO_date &at[73]
#define DIV_INFO_docs &at[74]
#define DIV_INFO_docs_E &at[75]

#define LINK_INFO &at[79]
#define LINK_INFO_id &at[80]
#define LINK_INFO_tag &at[81]
#define LINK_INFO_name &at[82]
#define LINK_INFO_descr &at[83]
#define LINK_INFO_dbfrom &at[84]
#define LINK_INFO_dbto &at[85]
#define LINK_INFO_datasize &at[86]
#define LINK_INFO_reciprocal &at[87]

#define FIELDDATA &at[34]
#define FIELDDATA_num_terms &at[35]
#define FIELDDATA_num_bucket &at[36]

#define POST_TYPE &at[64]

#define DOCSUM_SET &at[105]
#define DOCSUM_SET_E &at[106]
