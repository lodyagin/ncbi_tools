/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "mmdb1.l60";
static AsnValxNodePtr avn;
static AsnTypePtr at;
static AsnModulePtr amp;


/**************************************************
*
*    Defines for Module MMDB
*
**************************************************/

#define BIOSTRUC &at[0]
#define BIOSTRUC_id &at[1]
#define BIOSTRUC_id_E &at[2]
#define BIOSTRUC_descr &at[13]
#define BIOSTRUC_descr_E &at[14]
#define BIOSTRUC_chemical_graph &at[41]
#define BIOSTRUC_features &at[132]
#define BIOSTRUC_features_E &at[133]
#define BIOSTRUC_model &at[135]
#define BIOSTRUC_model_E &at[136]

#define BIOSTRUC_ID &at[3]
#define BIOSTRUC_ID_mmdb_id &at[4]
#define BIOSTRUC_ID_other_database &at[7]
#define BIOSTRUC_ID_local_id &at[9]

#define BIOSTRUC_SET &at[138]
#define BIOSTRUC_SET_id &at[139]
#define BIOSTRUC_SET_id_E &at[140]
#define BIOSTRUC_SET_descr &at[141]
#define BIOSTRUC_SET_descr_E &at[142]
#define BIOSTRUC_SET_biostrucs &at[143]
#define BIOSTRUC_SET_biostrucs_E &at[144]

#define BIOSTRUC_ANNOT_SET &at[145]
#define BIOSTRUC_ANNOT_SET_id &at[146]
#define BIOSTRUC_ANNOT_SET_id_E &at[147]
#define BIOSTRUC_ANNOT_SET_descr &at[148]
#define BIOSTRUC_ANNOT_SET_descr_E &at[149]
#define BIOSTRUC_ANNOT_SET_features &at[150]
#define BIOSTRUC_ANNOT_SET_features_E &at[151]

#define BIOSTRUC_RESIDUE_GRAPH_SET &at[152]
#define BIOSTRUC_RESIDUE_GRAPH_SET_id &at[153]
#define BIOSTRUC_RESIDUE_GRAPH_SET_id_E &at[154]
#define BIOSTRUC_RESIDUE_GRAPH_SET_descr &at[155]
#define BIOSTRUC_RESIDUE_GRAPH_SET_descr_E &at[156]
#define BIOSTRUC_RESIDUE_GRAPH_SET_residue_graphs &at[158]
#define BIOSTRUC_RESIDUE_GRAPH_SET_residue_graphs_E &at[159]

#define BIOSTRUC_DESCR &at[15]
#define BIOSTRUC_DESCR_name &at[16]
#define BIOSTRUC_DESCR_pdb_comment &at[18]
#define BIOSTRUC_DESCR_other_comment &at[19]
#define BIOSTRUC_DESCR_history &at[20]
#define BIOSTRUC_DESCR_attribution &at[39]

#define MMDB_ID &at[5]

#define BIOSTRUC_HISTORY &at[21]
#define BIOSTRUC_HISTORY_replaces &at[22]
#define BIOSTRUC_HISTORY_replaced_by &at[28]
#define BIOSTRUC_HISTORY_data_source &at[29]

#define BIOSTRUC_REPLACE &at[23]
#define BIOSTRUC_REPLACE_id &at[24]
#define BIOSTRUC_REPLACE_date &at[25]

#define BIOSTRUC_SOURCE &at[30]
#define BIOSTRUC_SOURCE_name_of_database &at[31]
#define BIOSTRUC_SOURCE_version_of_database &at[32]
#define BIOSTRUC_SOURCE_version_of_database_release_date &at[33]
#define BIOSTRUC_SOURCE_version_of_database_release_code &at[34]
#define BIOSTRUC_SOURCE_database_entry_id &at[35]
#define BIOSTRUC_SOURCE_database_entry_date &at[36]
#define BIOSTRUC_SOURCE_database_entry_history &at[37]
#define BIOSTRUC_SOURCE_database_entry_history_E &at[38]


/**************************************************
*
*    Defines for Module MMDB-Chemical-graph
*
**************************************************/

#define BIOSTRUC_GRAPH &at[43]
#define BIOSTRUC_GRAPH_descr &at[44]
#define BIOSTRUC_GRAPH_descr_E &at[45]
#define BIOSTRUC_GRAPH_molecule_graphs &at[58]
#define BIOSTRUC_GRAPH_molecule_graphs_E &at[59]
#define BIOSTRUC_GRAPH_inter_molecule_bonds &at[97]
#define BIOSTRUC_GRAPH_inter_molecule_bonds_E &at[98]
#define BIOSTRUC_GRAPH_residue_graphs &at[99]
#define BIOSTRUC_GRAPH_residue_graphs_E &at[100]

#define BIOMOL_DESCR &at[46]
#define BIOMOL_DESCR_name &at[47]
#define BIOMOL_DESCR_pdb_class &at[48]
#define BIOMOL_DESCR_pdb_source &at[49]
#define BIOMOL_DESCR_pdb_comment &at[50]
#define BIOMOL_DESCR_other_comment &at[51]
#define BIOMOL_DESCR_organism &at[52]
#define BIOMOL_DESCR_attribution &at[54]
#define BIOMOL_DESCR_assembly_type &at[56]
#define BIOMOL_DESCR_molecule_type &at[57]

#define RESIDUE_GRAPH &at[101]
#define RESIDUE_GRAPH_id &at[102]
#define RESIDUE_GRAPH_descr &at[103]
#define RESIDUE_GRAPH_descr_E &at[104]
#define RESIDUE_GRAPH_residue_type &at[105]
#define RESIDUE_GRAPH_iupac_code &at[106]
#define RESIDUE_GRAPH_iupac_code_E &at[107]
#define RESIDUE_GRAPH_atoms &at[108]
#define RESIDUE_GRAPH_atoms_E &at[109]
#define RESIDUE_GRAPH_bonds &at[118]
#define RESIDUE_GRAPH_bonds_E &at[119]
#define RESIDUE_GRAPH_chiral_centers &at[124]
#define RESIDUE_GRAPH_chiral_centers_E &at[125]

#define MOLECULE_ID &at[62]

#define RESIDUE_ID &at[71]

#define ATOM_ID &at[94]

#define MOLECULE_GRAPH &at[60]
#define MOLECULE_GRAPH_id &at[61]
#define MOLECULE_GRAPH_descr &at[63]
#define MOLECULE_GRAPH_descr_E &at[64]
#define MOLECULE_GRAPH_seq_id &at[65]
#define MOLECULE_GRAPH_residue_sequence &at[67]
#define MOLECULE_GRAPH_residue_sequence_E &at[68]
#define MOLECULE_GRAPH_inter_residue_bonds &at[86]
#define MOLECULE_GRAPH_inter_residue_bonds_E &at[87]

#define INTER_RESIDUE_BOND &at[88]
#define INTER_RESIDUE_BOND_atom_id_1 &at[89]
#define INTER_RESIDUE_BOND_atom_id_2 &at[95]
#define INTER_RESIDUE_BOND_bond_order &at[96]

#define RESIDUE &at[69]
#define RESIDUE_id &at[70]
#define RESIDUE_name &at[72]
#define RESIDUE_residue_graph &at[73]

#define RESIDUE_GRAPH_PNTR &at[74]
#define RESIDUE_GRAPH_PNTR_local &at[75]
#define RESIDUE_GRAPH_PNTR_biostruc &at[77]
#define RESIDUE_GRAPH_PNTR_standard &at[82]

#define RESIDUE_GRAPH_ID &at[76]

#define BIOSTRUC_GRAPH_PNTR &at[78]
#define BIOSTRUC_GRAPH_PNTR_biostruc_id &at[79]
#define BIOSTRUC_GRAPH_PNTR_residue_graph_id &at[81]

#define BIOSTRUC_RESIDUE_GRAPH_SET_PNTR &at[83]
#define BIOSTRUC_RESIDUE_GRAPH_SET_PNTR_biostruc_residue_graph_set_id &at[84]
#define BIOSTRUC_RESIDUE_GRAPH_SET_PNTR_residue_graph_id &at[85]

#define ATOM &at[110]
#define ATOM_id &at[111]
#define ATOM_name &at[112]
#define ATOM_iupac_code &at[113]
#define ATOM_iupac_code_E &at[114]
#define ATOM_element &at[115]
#define ATOM_ionizable_proton &at[117]

#define INTRA_RESIDUE_BOND &at[120]
#define INTRA_RESIDUE_BOND_atom_id_1 &at[121]
#define INTRA_RESIDUE_BOND_atom_id_2 &at[122]
#define INTRA_RESIDUE_BOND_bond_order &at[123]

#define CHIRAL_CENTER &at[126]
#define CHIRAL_CENTER_c &at[127]
#define CHIRAL_CENTER_n1 &at[128]
#define CHIRAL_CENTER_n2 &at[129]
#define CHIRAL_CENTER_n3 &at[130]
#define CHIRAL_CENTER_sign &at[131]

#define ATOM_PNTR &at[90]
#define ATOM_PNTR_molecule_id &at[91]
#define ATOM_PNTR_residue_id &at[92]
#define ATOM_PNTR_atom_id &at[93]

#define ATOM_PNTR_SET &at[161]
#define ATOM_PNTR_SET_E &at[162]
