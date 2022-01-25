/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "mmdb2.l60";
static AsnValxNodePtr avn;
static AsnTypePtr at;
static AsnModulePtr amp;


/**************************************************
*
*    Defines for Module MMDB-Structural-model
*
**************************************************/

#define BIOSTRUC_MODEL &at[0]
#define BIOSTRUC_MODEL_id &at[1]
#define BIOSTRUC_MODEL_type &at[4]
#define BIOSTRUC_MODEL_descr &at[6]
#define BIOSTRUC_MODEL_descr_E &at[7]
#define BIOSTRUC_MODEL_model_space &at[19]
#define BIOSTRUC_MODEL_model_coordinates &at[33]
#define BIOSTRUC_MODEL_model_coordinates_E &at[34]

#define MODEL_ID &at[2]

#define MODEL_COORDINATE_SET_ID &at[37]

#define MODEL_TYPE &at[5]

#define MODEL_DESCR &at[8]
#define MODEL_DESCR_name &at[9]
#define MODEL_DESCR_pdb_reso &at[11]
#define MODEL_DESCR_pdb_method &at[12]
#define MODEL_DESCR_pdb_comment &at[13]
#define MODEL_DESCR_other_comment &at[14]
#define MODEL_DESCR_attribution &at[15]

#define MODEL_SPACE &at[20]
#define MODEL_SPACE_coordinate_units &at[21]
#define MODEL_SPACE_thermal_factor_units &at[23]
#define MODEL_SPACE_occupancy_factor_units &at[24]
#define MODEL_SPACE_density_units &at[25]
#define MODEL_SPACE_reference_frame &at[26]

#define MODEL_COORDINATE_SET &at[35]
#define MODEL_COORDINATE_SET_id &at[36]
#define MODEL_COORDINATE_SET_descr &at[38]
#define MODEL_COORDINATE_SET_descr_E &at[39]
#define MODEL_COORDINATE_SET_coordinates &at[40]
#define MODEL_COORDINATE_SET_coordinates_literal &at[41]
#define MODEL_COORDINATE_SET_coordinates_reference &at[149]

#define REFERENCE_FRAME &at[27]
#define REFERENCE_FRAME_biostruc_id &at[28]
#define REFERENCE_FRAME_rotation_translation &at[30]

#define COORDINATES &at[42]
#define COORDINATES_atomic &at[43]
#define COORDINATES_surface &at[94]
#define COORDINATES_density &at[137]

#define ATOMIC_COORDINATES &at[44]
#define ATOMIC_COORDINATES_number_of_points &at[45]
#define ATOMIC_COORDINATES_atoms &at[46]
#define ATOMIC_COORDINATES_sites &at[48]
#define ATOMIC_COORDINATES_temperature_factors &at[57]
#define ATOMIC_COORDINATES_occupancies &at[79]
#define ATOMIC_COORDINATES_alternate_conf_ids &at[84]
#define ATOMIC_COORDINATES_conf_ensembles &at[88]
#define ATOMIC_COORDINATES_conf_ensembles_E &at[89]

#define SURFACE_COORDINATES &at[95]
#define SURFACE_COORDINATES_contents &at[96]
#define SURFACE_COORDINATES_surface &at[98]
#define SURFACE_COORDINATES_surface_sphere &at[99]
#define SURFACE_COORDINATES_surface_cone &at[101]
#define SURFACE_COORDINATES_surface_cylinder &at[103]
#define SURFACE_COORDINATES_surface_brick &at[105]
#define SURFACE_COORDINATES_surface_tmesh &at[107]
#define SURFACE_COORDINATES_surface_triangles &at[120]

#define DENSITY_COORDINATES &at[138]
#define DENSITY_COORDINATES_contents &at[139]
#define DENSITY_COORDINATES_grid_corners &at[140]
#define DENSITY_COORDINATES_grid_steps_x &at[141]
#define DENSITY_COORDINATES_grid_steps_y &at[142]
#define DENSITY_COORDINATES_grid_steps_z &at[143]
#define DENSITY_COORDINATES_fastest_varying &at[144]
#define DENSITY_COORDINATES_slowest_varying &at[145]
#define DENSITY_COORDINATES_scale_factor &at[146]
#define DENSITY_COORDINATES_density &at[147]
#define DENSITY_COORDINATES_density_E &at[148]

#define MODEL_SPACE_POINTS &at[49]
#define MODEL_SPACE_POINTS_scale_factor &at[50]
#define MODEL_SPACE_POINTS_x &at[51]
#define MODEL_SPACE_POINTS_x_E &at[52]
#define MODEL_SPACE_POINTS_y &at[53]
#define MODEL_SPACE_POINTS_y_E &at[54]
#define MODEL_SPACE_POINTS_z &at[55]
#define MODEL_SPACE_POINTS_z_E &at[56]

#define ATOMIC_TEMPERATURE_FACTORS &at[58]
#define ATOMIC_TEMPERATURE_FACTORS_isotropic &at[59]
#define ATOMIC_TEMPERATURE_FACTORS_anisotropic &at[64]

#define ATOMIC_OCCUPANCIES &at[80]
#define ATOMIC_OCCUPANCIES_scale_factor &at[81]
#define ATOMIC_OCCUPANCIES_o &at[82]
#define ATOMIC_OCCUPANCIES_o_E &at[83]

#define ALTERNATE_CONFORMATION_IDS &at[85]
#define ALTERNATE_CONFORMATION_IDS_E &at[86]

#define CONFORMATION_ENSEMBLE &at[90]
#define CONFORMATION_ENSEMBLE_name &at[91]
#define CONFORMATION_ENSEMBLE_alt_conf_ids &at[92]
#define CONFORMATION_ENSEMBLE_alt_conf_ids_E &at[93]

#define ISOTROPIC_TEMPERATURE_FACTORS &at[60]
#define ISOTROPIC_TEMPERATURE_FACTORS_scale_factor &at[61]
#define ISOTROPIC_TEMPERATURE_FACTORS_b &at[62]
#define ISOTROPIC_TEMPERATURE_FACTORS_b_E &at[63]

#define ANISOTROPIC_TEMPERATURE_FACTORS &at[65]
#define ANISOTROPIC_TEMPERATURE_FACTORS_scale_factor &at[66]
#define ANISOTROPIC_TEMPERATURE_FACTORS_b_11 &at[67]
#define ANISOTROPIC_TEMPERATURE_FACTORS_b_11_E &at[68]
#define ANISOTROPIC_TEMPERATURE_FACTORS_b_12 &at[69]
#define ANISOTROPIC_TEMPERATURE_FACTORS_b_12_E &at[70]
#define ANISOTROPIC_TEMPERATURE_FACTORS_b_13 &at[71]
#define ANISOTROPIC_TEMPERATURE_FACTORS_b_13_E &at[72]
#define ANISOTROPIC_TEMPERATURE_FACTORS_b_22 &at[73]
#define ANISOTROPIC_TEMPERATURE_FACTORS_b_22_E &at[74]
#define ANISOTROPIC_TEMPERATURE_FACTORS_b_23 &at[75]
#define ANISOTROPIC_TEMPERATURE_FACTORS_b_23_E &at[76]
#define ANISOTROPIC_TEMPERATURE_FACTORS_b_33 &at[77]
#define ANISOTROPIC_TEMPERATURE_FACTORS_b_33_E &at[78]

#define ALTERNATE_CONFORMATION_ID &at[87]

#define T_MESH &at[108]
#define T_MESH_number_of_points &at[109]
#define T_MESH_scale_factor &at[110]
#define T_MESH_swap &at[111]
#define T_MESH_swap_E &at[112]
#define T_MESH_x &at[114]
#define T_MESH_x_E &at[115]
#define T_MESH_y &at[116]
#define T_MESH_y_E &at[117]
#define T_MESH_z &at[118]
#define T_MESH_z_E &at[119]

#define TRIANGLES &at[121]
#define TRIANGLES_number_of_points &at[122]
#define TRIANGLES_scale_factor &at[123]
#define TRIANGLES_x &at[124]
#define TRIANGLES_x_E &at[125]
#define TRIANGLES_y &at[126]
#define TRIANGLES_y_E &at[127]
#define TRIANGLES_z &at[128]
#define TRIANGLES_z_E &at[129]
#define TRIANGLES_number_of_triangles &at[130]
#define TRIANGLES_v1 &at[131]
#define TRIANGLES_v1_E &at[132]
#define TRIANGLES_v2 &at[133]
#define TRIANGLES_v2_E &at[134]
#define TRIANGLES_v3 &at[135]
#define TRIANGLES_v3_E &at[136]
