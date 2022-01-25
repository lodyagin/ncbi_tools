/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "mmdb3.l61";
static AsnValxNodePtr avn;
static AsnTypePtr at;
static AsnModulePtr amp;


/**************************************************
*
*    Defines for Module MMDB-Features
*
**************************************************/

#define BIOSTRUC_FEATURE_SET &at[0]
#define BIOSTRUC_FEATURE_SET_id &at[1]
#define BIOSTRUC_FEATURE_SET_descr &at[4]
#define BIOSTRUC_FEATURE_SET_descr_E &at[5]
#define BIOSTRUC_FEATURE_SET_features &at[15]
#define BIOSTRUC_FEATURE_SET_features_E &at[16]

#define CHEM_GRAPH_PNTRS &at[100]
#define CHEM_GRAPH_PNTRS_atoms &at[101]
#define CHEM_GRAPH_PNTRS_residues &at[113]
#define CHEM_GRAPH_PNTRS_molecules &at[128]

#define ATOM_PNTRS &at[102]
#define ATOM_PNTRS_number_of_ptrs &at[103]
#define ATOM_PNTRS_molecule_ids &at[104]
#define ATOM_PNTRS_molecule_ids_E &at[105]
#define ATOM_PNTRS_residue_ids &at[107]
#define ATOM_PNTRS_residue_ids_E &at[108]
#define ATOM_PNTRS_atom_ids &at[110]
#define ATOM_PNTRS_atom_ids_E &at[111]

#define CHEM_GRAPH_ALIGNMENT &at[174]
#define CHEM_GRAPH_ALIGNMENT_dimension &at[175]
#define CHEM_GRAPH_ALIGNMENT_biostruc_ids &at[176]
#define CHEM_GRAPH_ALIGNMENT_biostruc_ids_E &at[177]
#define CHEM_GRAPH_ALIGNMENT_alignment &at[178]
#define CHEM_GRAPH_ALIGNMENT_alignment_E &at[179]
#define CHEM_GRAPH_ALIGNMENT_domain &at[180]
#define CHEM_GRAPH_ALIGNMENT_domain_E &at[181]
#define CHEM_GRAPH_ALIGNMENT_transform &at[182]
#define CHEM_GRAPH_ALIGNMENT_transform_E &at[183]
#define CHEM_GRAPH_ALIGNMENT_aligndata &at[184]
#define CHEM_GRAPH_ALIGNMENT_aligndata_E &at[185]

#define SPHERE &at[150]
#define SPHERE_center &at[151]
#define SPHERE_radius &at[152]

#define CONE &at[154]
#define CONE_axis_top &at[155]
#define CONE_axis_bottom &at[156]
#define CONE_radius_bottom &at[157]

#define CYLINDER &at[159]
#define CYLINDER_axis_top &at[160]
#define CYLINDER_axis_bottom &at[161]
#define CYLINDER_radius &at[162]

#define BRICK &at[164]
#define BRICK_corner_000 &at[165]
#define BRICK_corner_001 &at[166]
#define BRICK_corner_010 &at[167]
#define BRICK_corner_011 &at[168]
#define BRICK_corner_100 &at[169]
#define BRICK_corner_101 &at[170]
#define BRICK_corner_110 &at[171]
#define BRICK_corner_111 &at[172]

#define TRANSFORM &at[33]
#define TRANSFORM_id &at[34]
#define TRANSFORM_moves &at[35]
#define TRANSFORM_moves_E &at[36]

#define BIOSTRUC_FEATURE_SET_ID &at[2]

#define BIOSTRUC_FEATURE_ID &at[19]

#define BIOSTRUC_FEATURE_SET_DESCR &at[6]
#define BIOSTRUC_FEATURE_SET_DESCR_name &at[7]
#define BIOSTRUC_FEATURE_SET_DESCR_pdb_comment &at[9]
#define BIOSTRUC_FEATURE_SET_DESCR_other_comment &at[10]
#define BIOSTRUC_FEATURE_SET_DESCR_attribution &at[11]

#define BIOSTRUC_FEATURE &at[17]
#define BIOSTRUC_FEATURE_id &at[18]
#define BIOSTRUC_FEATURE_name &at[20]
#define BIOSTRUC_FEATURE_type &at[21]
#define BIOSTRUC_FEATURE_property &at[22]
#define BIOSTRUC_FEATURE_property_color &at[23]
#define BIOSTRUC_FEATURE_property_render &at[30]
#define BIOSTRUC_FEATURE_property_transform &at[32]
#define BIOSTRUC_FEATURE_property_camera &at[56]
#define BIOSTRUC_FEATURE_property_script &at[76]
#define BIOSTRUC_FEATURE_property_user &at[96]
#define BIOSTRUC_FEATURE_location &at[98]
#define BIOSTRUC_FEATURE_location_subgraph &at[99]
#define BIOSTRUC_FEATURE_location_region &at[133]
#define BIOSTRUC_FEATURE_location_alignment &at[173]
#define BIOSTRUC_FEATURE_location_similarity &at[196]
#define BIOSTRUC_FEATURE_location_indirect &at[205]

#define COLOR_PROP &at[24]
#define COLOR_PROP_r &at[25]
#define COLOR_PROP_g &at[26]
#define COLOR_PROP_b &at[27]
#define COLOR_PROP_name &at[28]

#define RENDER_PROP &at[31]

#define CAMERA &at[57]
#define CAMERA_mode &at[58]
#define CAMERA_x &at[59]
#define CAMERA_y &at[65]
#define CAMERA_z &at[66]
#define CAMERA_up &at[67]
#define CAMERA_fore &at[68]
#define CAMERA_norm &at[69]
#define CAMERA_center &at[70]
#define CAMERA_tooclose &at[71]
#define CAMERA_toofar &at[75]

#define BIOSTRUC_SCRIPT &at[77]
#define BIOSTRUC_SCRIPT_E &at[78]

#define REGION_PNTRS &at[134]
#define REGION_PNTRS_model_id &at[135]
#define REGION_PNTRS_region &at[137]
#define REGION_PNTRS_region_site &at[138]
#define REGION_PNTRS_region_site_E &at[139]
#define REGION_PNTRS_region_boundary &at[146]
#define REGION_PNTRS_region_boundary_E &at[147]

#define REGION_SIMILARITY &at[197]
#define REGION_SIMILARITY_dimension &at[198]
#define REGION_SIMILARITY_biostruc_ids &at[199]
#define REGION_SIMILARITY_biostruc_ids_E &at[200]
#define REGION_SIMILARITY_similarity &at[201]
#define REGION_SIMILARITY_similarity_E &at[202]
#define REGION_SIMILARITY_transform &at[203]
#define REGION_SIMILARITY_transform_E &at[204]

#define OTHER_FEATURE &at[85]
#define OTHER_FEATURE_biostruc_id &at[86]
#define OTHER_FEATURE_set &at[88]
#define OTHER_FEATURE_feature &at[89]

#define RESIDUE_PNTRS &at[114]
#define RESIDUE_PNTRS_explicit &at[115]
#define RESIDUE_PNTRS_interval &at[122]
#define RESIDUE_PNTRS_interval_E &at[123]

#define MOLECULE_PNTRS &at[129]
#define MOLECULE_PNTRS_number_of_ptrs &at[130]
#define MOLECULE_PNTRS_molecule_ids &at[131]
#define MOLECULE_PNTRS_molecule_ids_E &at[132]

#define RESIDUE_EXPLICIT_PNTRS &at[116]
#define RESIDUE_EXPLICIT_PNTRS_number_of_ptrs &at[117]
#define RESIDUE_EXPLICIT_PNTRS_molecule_ids &at[118]
#define RESIDUE_EXPLICIT_PNTRS_molecule_ids_E &at[119]
#define RESIDUE_EXPLICIT_PNTRS_residue_ids &at[120]
#define RESIDUE_EXPLICIT_PNTRS_residue_ids_E &at[121]

#define RESIDUE_INTERVAL_PNTR &at[124]
#define RESIDUE_INTERVAL_PNTR_molecule_id &at[125]
#define RESIDUE_INTERVAL_PNTR_from &at[126]
#define RESIDUE_INTERVAL_PNTR_to &at[127]

#define REGION_COORDINATES &at[140]
#define REGION_COORDINATES_model_coord_set_id &at[141]
#define REGION_COORDINATES_number_of_coords &at[143]
#define REGION_COORDINATES_coordinate_indices &at[144]
#define REGION_COORDINATES_coordinate_indices_E &at[145]

#define REGION_BOUNDARY &at[148]
#define REGION_BOUNDARY_sphere &at[149]
#define REGION_BOUNDARY_cone &at[153]
#define REGION_BOUNDARY_cylinder &at[158]
#define REGION_BOUNDARY_brick &at[163]

#define ALIGN_STATS &at[186]
#define ALIGN_STATS_descr &at[187]
#define ALIGN_STATS_scale_factor &at[188]
#define ALIGN_STATS_vast_score &at[189]
#define ALIGN_STATS_vast_mlogp &at[190]
#define ALIGN_STATS_align_res &at[191]
#define ALIGN_STATS_rmsd &at[192]
#define ALIGN_STATS_blast_score &at[193]
#define ALIGN_STATS_blast_mlogp &at[194]
#define ALIGN_STATS_other_score &at[195]

#define MODEL_SPACE_POINT &at[60]
#define MODEL_SPACE_POINT_scale_factor &at[61]
#define MODEL_SPACE_POINT_x &at[62]
#define MODEL_SPACE_POINT_y &at[63]
#define MODEL_SPACE_POINT_z &at[64]

#define REALVALUE &at[72]
#define REALVALUE_scale_factor &at[73]
#define REALVALUE_scaled_integer_value &at[74]

#define MOVE &at[37]
#define MOVE_rotate &at[38]
#define MOVE_translate &at[50]

#define ROT_MATRIX &at[39]
#define ROT_MATRIX_scale_factor &at[40]
#define ROT_MATRIX_rot_11 &at[41]
#define ROT_MATRIX_rot_12 &at[42]
#define ROT_MATRIX_rot_13 &at[43]
#define ROT_MATRIX_rot_21 &at[44]
#define ROT_MATRIX_rot_22 &at[45]
#define ROT_MATRIX_rot_23 &at[46]
#define ROT_MATRIX_rot_31 &at[47]
#define ROT_MATRIX_rot_32 &at[48]
#define ROT_MATRIX_rot_33 &at[49]

#define TRANS_MATRIX &at[51]
#define TRANS_MATRIX_scale_factor &at[52]
#define TRANS_MATRIX_tran_1 &at[53]
#define TRANS_MATRIX_tran_2 &at[54]
#define TRANS_MATRIX_tran_3 &at[55]

#define BIOSTRUC_SCRIPT_STEP &at[79]
#define BIOSTRUC_SCRIPT_STEP_step_id &at[80]
#define BIOSTRUC_SCRIPT_STEP_step_name &at[82]
#define BIOSTRUC_SCRIPT_STEP_feature_do &at[83]
#define BIOSTRUC_SCRIPT_STEP_feature_do_E &at[84]
#define BIOSTRUC_SCRIPT_STEP_camera_move &at[90]
#define BIOSTRUC_SCRIPT_STEP_pause &at[91]
#define BIOSTRUC_SCRIPT_STEP_waitevent &at[92]
#define BIOSTRUC_SCRIPT_STEP_extra &at[94]
#define BIOSTRUC_SCRIPT_STEP_jump &at[95]

#define STEP_ID &at[81]
