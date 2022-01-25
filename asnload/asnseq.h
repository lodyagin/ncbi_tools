/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "asnseq.l60";
static AsnValxNodePtr avn;
static AsnTypePtr at;
static AsnModulePtr amp;


/**************************************************
*
*    Defines for Module NCBI-Sequence
*
**************************************************/

#define BIOSEQ &at[0]
#define BIOSEQ_id &at[1]
#define BIOSEQ_id_E &at[2]
#define BIOSEQ_descr &at[5]
#define BIOSEQ_inst &at[97]
#define BIOSEQ_annot &at[164]
#define BIOSEQ_annot_E &at[165]

#define SEQ_ANNOT &at[166]
#define SEQ_ANNOT_id &at[167]
#define SEQ_ANNOT_id_E &at[168]
#define SEQ_ANNOT_db &at[174]
#define SEQ_ANNOT_name &at[175]
#define SEQ_ANNOT_desc &at[176]
#define SEQ_ANNOT_data &at[194]
#define SEQ_ANNOT_data_ftable &at[195]
#define SEQ_ANNOT_data_ftable_E &at[196]
#define SEQ_ANNOT_data_align &at[197]
#define SEQ_ANNOT_data_align_E &at[198]
#define SEQ_ANNOT_data_graph &at[199]
#define SEQ_ANNOT_data_graph_E &at[200]
#define SEQ_ANNOT_data_ids &at[202]
#define SEQ_ANNOT_data_ids_E &at[203]
#define SEQ_ANNOT_data_locs &at[204]
#define SEQ_ANNOT_data_locs_E &at[205]

#define PUBDESC &at[58]
#define PUBDESC_pub &at[59]
#define PUBDESC_name &at[61]
#define PUBDESC_fig &at[62]
#define PUBDESC_num &at[63]
#define PUBDESC_numexc &at[64]
#define PUBDESC_poly_a &at[65]
#define PUBDESC_maploc &at[66]
#define PUBDESC_seq_raw &at[67]
#define PUBDESC_align_group &at[69]
#define PUBDESC_comment &at[70]
#define PUBDESC_reftype &at[71]

#define SEQ_DESCR &at[6]
#define SEQ_DESCR_E &at[7]

#define SEQDESC &at[8]
#define SEQDESC_mol_type &at[9]
#define SEQDESC_modif &at[12]
#define SEQDESC_modif_E &at[13]
#define SEQDESC_method &at[15]
#define SEQDESC_name &at[17]
#define SEQDESC_title &at[19]
#define SEQDESC_org &at[20]
#define SEQDESC_comment &at[22]
#define SEQDESC_num &at[23]
#define SEQDESC_maploc &at[51]
#define SEQDESC_pir &at[53]
#define SEQDESC_genbank &at[55]
#define SEQDESC_pub &at[57]
#define SEQDESC_region &at[72]
#define SEQDESC_user &at[73]
#define SEQDESC_sp &at[75]
#define SEQDESC_dbxref &at[77]
#define SEQDESC_embl &at[78]
#define SEQDESC_create_date &at[80]
#define SEQDESC_update_date &at[82]
#define SEQDESC_prf &at[83]
#define SEQDESC_pdb &at[85]
#define SEQDESC_het &at[87]
#define SEQDESC_source &at[89]
#define SEQDESC_molinfo &at[91]

#define NUMBERING &at[24]
#define NUMBERING_cont &at[25]
#define NUMBERING_enum &at[33]
#define NUMBERING_ref &at[39]
#define NUMBERING_real &at[44]

#define HETEROGEN &at[88]

#define SEQ_HIST &at[152]
#define SEQ_HIST_assembly &at[153]
#define SEQ_HIST_assembly_E &at[154]
#define SEQ_HIST_replaces &at[155]
#define SEQ_HIST_replaced_by &at[160]
#define SEQ_HIST_deleted &at[161]
#define SEQ_HIST_deleted_bool &at[162]
#define SEQ_HIST_deleted_date &at[163]

#define SEQ_INST &at[98]
#define SEQ_INST_repr &at[99]
#define SEQ_INST_mol &at[100]
#define SEQ_INST_length &at[101]
#define SEQ_INST_fuzz &at[102]
#define SEQ_INST_topology &at[104]
#define SEQ_INST_strand &at[105]
#define SEQ_INST_seq_data &at[106]
#define SEQ_INST_ext &at[129]
#define SEQ_INST_hist &at[151]

#define GIBB_MOL &at[10]

#define GIBB_MOD &at[14]

#define GIBB_METHOD &at[16]

#define MOLINFO &at[92]
#define MOLINFO_biomol &at[93]
#define MOLINFO_tech &at[94]
#define MOLINFO_techexp &at[95]
#define MOLINFO_completeness &at[96]

#define NUM_CONT &at[26]
#define NUM_CONT_refnum &at[27]
#define NUM_CONT_has_zero &at[29]
#define NUM_CONT_ascending &at[31]

#define NUM_ENUM &at[34]
#define NUM_ENUM_num &at[35]
#define NUM_ENUM_names &at[36]
#define NUM_ENUM_names_E &at[37]

#define NUM_REF &at[40]
#define NUM_REF_type &at[41]
#define NUM_REF_aligns &at[42]

#define NUM_REAL &at[45]
#define NUM_REAL_a &at[46]
#define NUM_REAL_b &at[48]
#define NUM_REAL_units &at[49]

#define SEQ_DATA &at[107]
#define SEQ_DATA_iupacna &at[108]
#define SEQ_DATA_iupacaa &at[110]
#define SEQ_DATA_ncbi2na &at[112]
#define SEQ_DATA_ncbi4na &at[115]
#define SEQ_DATA_ncbi8na &at[117]
#define SEQ_DATA_ncbipna &at[119]
#define SEQ_DATA_ncbi8aa &at[121]
#define SEQ_DATA_ncbieaa &at[123]
#define SEQ_DATA_ncbipaa &at[125]
#define SEQ_DATA_ncbistdaa &at[127]

#define SEQ_EXT &at[130]
#define SEQ_EXT_seg &at[131]
#define SEQ_EXT_ref &at[135]
#define SEQ_EXT_map &at[137]
#define SEQ_EXT_delta &at[141]

#define SEG_EXT &at[132]
#define SEG_EXT_E &at[133]

#define REF_EXT &at[136]

#define MAP_EXT &at[138]
#define MAP_EXT_E &at[139]

#define DELTA_EXT &at[142]
#define DELTA_EXT_E &at[143]

#define DELTA_SEQ &at[144]
#define DELTA_SEQ_loc &at[145]
#define DELTA_SEQ_literal &at[146]

#define SEQ_LITERAL &at[147]
#define SEQ_LITERAL_length &at[148]
#define SEQ_LITERAL_fuzz &at[149]
#define SEQ_LITERAL_seq_data &at[150]

#define SEQ_HIST_REC &at[156]
#define SEQ_HIST_REC_date &at[157]
#define SEQ_HIST_REC_ids &at[158]
#define SEQ_HIST_REC_ids_E &at[159]

#define IUPACNA &at[109]

#define IUPACAA &at[111]

#define NCBI2NA &at[113]

#define NCBI4NA &at[116]

#define NCBI8NA &at[118]

#define NCBIPNA &at[120]

#define NCBI8AA &at[122]

#define NCBIEAA &at[124]

#define NCBIPAA &at[126]

#define NCBISTDAA &at[128]

#define ANNOT_ID &at[169]
#define ANNOT_ID_local &at[170]
#define ANNOT_ID_ncbi &at[172]
#define ANNOT_ID_general &at[173]

#define ANNOT_DESCR &at[177]
#define ANNOT_DESCR_E &at[178]

#define ANNOTDESC &at[179]
#define ANNOTDESC_name &at[180]
#define ANNOTDESC_title &at[181]
#define ANNOTDESC_comment &at[182]
#define ANNOTDESC_pub &at[183]
#define ANNOTDESC_user &at[184]
#define ANNOTDESC_create_date &at[185]
#define ANNOTDESC_update_date &at[186]
#define ANNOTDESC_src &at[187]
#define ANNOTDESC_align &at[188]
#define ANNOTDESC_region &at[193]

#define ALIGN_DEF &at[189]
#define ALIGN_DEF_align_type &at[190]
#define ALIGN_DEF_ids &at[191]
#define ALIGN_DEF_ids_E &at[192]
