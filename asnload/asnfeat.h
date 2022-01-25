/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "asnfeat.l60";
static AsnValxNodePtr avn;
static AsnTypePtr at;
static AsnModulePtr amp;


/**************************************************
*
*    Defines for Module NCBI-Seqfeat
*
**************************************************/

#define SEQ_FEAT &at[0]
#define SEQ_FEAT_id &at[1]
#define SEQ_FEAT_data &at[12]
#define SEQ_FEAT_partial &at[207]
#define SEQ_FEAT_except &at[208]
#define SEQ_FEAT_comment &at[209]
#define SEQ_FEAT_product &at[210]
#define SEQ_FEAT_location &at[211]
#define SEQ_FEAT_qual &at[212]
#define SEQ_FEAT_qual_E &at[213]
#define SEQ_FEAT_title &at[217]
#define SEQ_FEAT_ext &at[218]
#define SEQ_FEAT_cit &at[219]
#define SEQ_FEAT_exp_ev &at[221]
#define SEQ_FEAT_xref &at[222]
#define SEQ_FEAT_xref_E &at[223]
#define SEQ_FEAT_dbxref &at[227]
#define SEQ_FEAT_dbxref_E &at[228]
#define SEQ_FEAT_pseudo &at[229]
#define SEQ_FEAT_except_text &at[230]

#define FEAT_ID &at[2]
#define FEAT_ID_gibb &at[3]
#define FEAT_ID_giim &at[5]
#define FEAT_ID_local &at[7]
#define FEAT_ID_general &at[9]

#define GENETIC_CODE &at[84]
#define GENETIC_CODE_E &at[85]
#define GENETIC_CODE_E_name &at[86]
#define GENETIC_CODE_E_id &at[87]
#define GENETIC_CODE_E_ncbieaa &at[88]
#define GENETIC_CODE_E_ncbi8aa &at[89]
#define GENETIC_CODE_E_ncbistdaa &at[91]
#define GENETIC_CODE_E_sncbieaa &at[92]
#define GENETIC_CODE_E_sncbi8aa &at[93]
#define GENETIC_CODE_E_sncbistdaa &at[94]

#define SEQFEATDATA &at[13]
#define SEQFEATDATA_gene &at[14]
#define SEQFEATDATA_org &at[31]
#define SEQFEATDATA_cdregion &at[75]
#define SEQFEATDATA_prot &at[105]
#define SEQFEATDATA_rna &at[119]
#define SEQFEATDATA_pub &at[137]
#define SEQFEATDATA_seq &at[139]
#define SEQFEATDATA_imp &at[140]
#define SEQFEATDATA_region &at[145]
#define SEQFEATDATA_comment &at[146]
#define SEQFEATDATA_bond &at[148]
#define SEQFEATDATA_site &at[149]
#define SEQFEATDATA_rsite &at[150]
#define SEQFEATDATA_user &at[156]
#define SEQFEATDATA_txinit &at[158]
#define SEQFEATDATA_num &at[187]
#define SEQFEATDATA_psec_str &at[189]
#define SEQFEATDATA_non_std_residue &at[190]
#define SEQFEATDATA_het &at[191]
#define SEQFEATDATA_biosrc &at[193]

#define GB_QUAL &at[214]
#define GB_QUAL_qual &at[215]
#define GB_QUAL_val &at[216]

#define SEQFEATXREF &at[224]
#define SEQFEATXREF_id &at[225]
#define SEQFEATXREF_data &at[226]

#define CDREGION &at[76]
#define CDREGION_orf &at[77]
#define CDREGION_frame &at[78]
#define CDREGION_conflict &at[80]
#define CDREGION_gaps &at[81]
#define CDREGION_mismatch &at[82]
#define CDREGION_code &at[83]
#define CDREGION_code_break &at[95]
#define CDREGION_code_break_E &at[96]
#define CDREGION_stops &at[104]

#define IMP_FEAT &at[141]
#define IMP_FEAT_key &at[142]
#define IMP_FEAT_loc &at[143]
#define IMP_FEAT_descr &at[144]

#define CODE_BREAK &at[97]
#define CODE_BREAK_loc &at[98]
#define CODE_BREAK_aa &at[100]
#define CODE_BREAK_aa_ncbieaa &at[101]
#define CODE_BREAK_aa_ncbi8aa &at[102]
#define CODE_BREAK_aa_ncbistdaa &at[103]

#define GENETIC_CODE_TABLE &at[231]
#define GENETIC_CODE_TABLE_E &at[232]


/**************************************************
*
*    Defines for Module NCBI-Rsite
*
**************************************************/

#define RSITE_REF &at[152]
#define RSITE_REF_str &at[153]
#define RSITE_REF_db &at[154]


/**************************************************
*
*    Defines for Module NCBI-RNA
*
**************************************************/

#define RNA_REF &at[121]
#define RNA_REF_type &at[122]
#define RNA_REF_pseudo &at[123]
#define RNA_REF_ext &at[124]
#define RNA_REF_ext_name &at[125]
#define RNA_REF_ext_tRNA &at[126]

#define TRNA_EXT &at[127]
#define TRNA_EXT_aa &at[128]
#define TRNA_EXT_aa_iupacaa &at[129]
#define TRNA_EXT_aa_ncbieaa &at[130]
#define TRNA_EXT_aa_ncbi8aa &at[131]
#define TRNA_EXT_aa_ncbistdaa &at[132]
#define TRNA_EXT_codon &at[133]
#define TRNA_EXT_codon_E &at[134]
#define TRNA_EXT_anticodon &at[135]


/**************************************************
*
*    Defines for Module NCBI-Gene
*
**************************************************/

#define GENE_REF &at[16]
#define GENE_REF_locus &at[17]
#define GENE_REF_allele &at[19]
#define GENE_REF_desc &at[20]
#define GENE_REF_maploc &at[21]
#define GENE_REF_pseudo &at[22]
#define GENE_REF_db &at[24]
#define GENE_REF_db_E &at[25]
#define GENE_REF_syn &at[28]
#define GENE_REF_syn_E &at[29]


/**************************************************
*
*    Defines for Module NCBI-Organism
*
**************************************************/

#define ORG_REF &at[33]
#define ORG_REF_taxname &at[34]
#define ORG_REF_common &at[35]
#define ORG_REF_mod &at[36]
#define ORG_REF_mod_E &at[37]
#define ORG_REF_db &at[38]
#define ORG_REF_db_E &at[39]
#define ORG_REF_syn &at[41]
#define ORG_REF_syn_E &at[42]
#define ORG_REF_orgname &at[43]

#define ORGNAME &at[44]
#define ORGNAME_name &at[45]
#define ORGNAME_name_binomial &at[46]
#define ORGNAME_name_virus &at[51]
#define ORGNAME_name_hybrid &at[52]
#define ORGNAME_name_namedhybrid &at[56]
#define ORGNAME_name_partial &at[57]
#define ORGNAME_attrib &at[64]
#define ORGNAME_mod &at[65]
#define ORGNAME_mod_E &at[66]
#define ORGNAME_lineage &at[71]
#define ORGNAME_gcode &at[72]
#define ORGNAME_mgcode &at[73]
#define ORGNAME_div &at[74]

#define BINOMIALORGNAME &at[47]
#define BINOMIALORGNAME_genus &at[48]
#define BINOMIALORGNAME_species &at[49]
#define BINOMIALORGNAME_subspecies &at[50]

#define MULTIORGNAME &at[53]
#define MULTIORGNAME_E &at[54]

#define PARTIALORGNAME &at[58]
#define PARTIALORGNAME_E &at[59]

#define ORGMOD &at[67]
#define ORGMOD_subtype &at[68]
#define ORGMOD_subname &at[69]
#define ORGMOD_attrib &at[70]

#define TAXELEMENT &at[60]
#define TAXELEMENT_fixed_level &at[61]
#define TAXELEMENT_level &at[62]
#define TAXELEMENT_name &at[63]


/**************************************************
*
*    Defines for Module NCBI-BioSource
*
**************************************************/

#define BIOSOURCE &at[195]
#define BIOSOURCE_genome &at[196]
#define BIOSOURCE_origin &at[197]
#define BIOSOURCE_org &at[198]
#define BIOSOURCE_subtype &at[200]
#define BIOSOURCE_subtype_E &at[201]
#define BIOSOURCE_is_focus &at[206]

#define SUBSOURCE &at[202]
#define SUBSOURCE_subtype &at[203]
#define SUBSOURCE_name &at[204]
#define SUBSOURCE_attrib &at[205]


/**************************************************
*
*    Defines for Module NCBI-Protein
*
**************************************************/

#define PROT_REF &at[107]
#define PROT_REF_name &at[108]
#define PROT_REF_name_E &at[109]
#define PROT_REF_desc &at[110]
#define PROT_REF_ec &at[111]
#define PROT_REF_ec_E &at[112]
#define PROT_REF_activity &at[113]
#define PROT_REF_activity_E &at[114]
#define PROT_REF_db &at[115]
#define PROT_REF_db_E &at[116]
#define PROT_REF_processed &at[118]


/**************************************************
*
*    Defines for Module NCBI-TxInit
*
**************************************************/

#define TXINIT &at[160]
#define TXINIT_name &at[161]
#define TXINIT_syn &at[162]
#define TXINIT_syn_E &at[163]
#define TXINIT_gene &at[164]
#define TXINIT_gene_E &at[165]
#define TXINIT_protein &at[167]
#define TXINIT_protein_E &at[168]
#define TXINIT_rna &at[170]
#define TXINIT_rna_E &at[171]
#define TXINIT_expression &at[172]
#define TXINIT_txsystem &at[173]
#define TXINIT_txdescr &at[174]
#define TXINIT_txorg &at[175]
#define TXINIT_mapping_precise &at[177]
#define TXINIT_location_accurate &at[178]
#define TXINIT_inittype &at[179]
#define TXINIT_evidence &at[180]
#define TXINIT_evidence_E &at[181]

#define TX_EVIDENCE &at[182]
#define TX_EVIDENCE_exp_code &at[183]
#define TX_EVIDENCE_expression_system &at[184]
#define TX_EVIDENCE_low_prec_data &at[185]
#define TX_EVIDENCE_from_homolog &at[186]
