/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "cn3d.h14";
static AsnValxNode avnx[26] = {
    {20,"off" ,1,0.0,&avnx[1] } ,
    {20,"trace" ,2,0.0,&avnx[2] } ,
    {20,"partial" ,3,0.0,&avnx[3] } ,
    {20,"complete" ,4,0.0,NULL } ,
    {20,"wire" ,1,0.0,&avnx[5] } ,
    {20,"tubes" ,2,0.0,&avnx[6] } ,
    {20,"ball-and-stick" ,3,0.0,&avnx[7] } ,
    {20,"space-fill" ,4,0.0,&avnx[8] } ,
    {20,"wire-worm" ,5,0.0,&avnx[9] } ,
    {20,"tube-worm" ,6,0.0,&avnx[10] } ,
    {20,"with-arrows" ,7,0.0,&avnx[11] } ,
    {20,"without-arrows" ,8,0.0,NULL } ,
    {20,"element" ,1,0.0,&avnx[13] } ,
    {20,"object" ,2,0.0,&avnx[14] } ,
    {20,"molecule" ,3,0.0,&avnx[15] } ,
    {20,"domain" ,4,0.0,&avnx[16] } ,
    {20,"secondary-structure" ,5,0.0,&avnx[17] } ,
    {20,"user-select" ,6,0.0,&avnx[18] } ,
    {20,"aligned" ,7,0.0,&avnx[19] } ,
    {20,"identity" ,8,0.0,&avnx[20] } ,
    {20,"variety" ,9,0.0,&avnx[21] } ,
    {20,"weighted-variety" ,10,0.0,&avnx[22] } ,
    {20,"information-content" ,11,0.0,&avnx[23] } ,
    {20,"fit" ,12,0.0,NULL } ,
    {3,NULL,255,0.0,NULL } ,
    {3,NULL,255,0.0,NULL } };

static AsnType atx[80] = {
  {401, "Cn3d-style-dictionary" ,1,0,0,0,0,1,0,0,NULL,&atx[20],&atx[1],0,&atx[55]} ,
  {0, "global-style" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[48]} ,
  {412, "Cn3d-style-settings" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[3],0,&atx[52]} ,
  {0, "protein-backbone" ,128,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[21]} ,
  {410, "Cn3d-backbone-style" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[5],0,&atx[23]} ,
  {0, "type" ,128,0,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[8]} ,
  {406, "Cn3d-backbone-type" ,1,0,0,0,0,0,0,0,NULL,&atx[7],&avnx[0],0,&atx[9]} ,
  {310, "ENUMERATED" ,0,10,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "style" ,128,1,0,0,0,0,0,0,NULL,&atx[9],NULL,0,&atx[10]} ,
  {407, "Cn3d-drawing-style" ,1,0,0,0,0,0,0,0,NULL,&atx[7],&avnx[4],0,&atx[11]} ,
  {0, "color-scheme" ,128,2,0,0,0,0,0,0,NULL,&atx[11],NULL,0,&atx[12]} ,
  {408, "Cn3d-color-scheme" ,1,0,0,0,0,0,0,0,NULL,&atx[7],&avnx[12],0,&atx[13]} ,
  {0, "user-color" ,128,3,0,0,0,0,0,0,NULL,&atx[13],NULL,0,NULL} ,
  {409, "Cn3d-color" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[14],0,&atx[4]} ,
  {0, "scale-factor" ,128,0,0,0,1,0,0,0,&avnx[24],&atx[15],NULL,0,&atx[16]} ,
  {302, "INTEGER" ,0,2,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "red" ,128,1,0,0,0,0,0,0,NULL,&atx[15],NULL,0,&atx[17]} ,
  {0, "green" ,128,2,0,0,0,0,0,0,NULL,&atx[15],NULL,0,&atx[18]} ,
  {0, "blue" ,128,3,0,0,0,0,0,0,NULL,&atx[15],NULL,0,&atx[19]} ,
  {0, "alpha" ,128,4,0,0,1,0,0,0,&avnx[25],&atx[15],NULL,0,NULL} ,
  {311, "SEQUENCE" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "nucleotide-backbone" ,128,1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[22]} ,
  {0, "protein-sidechains" ,128,2,0,0,0,0,0,0,NULL,&atx[23],NULL,0,&atx[29]} ,
  {411, "Cn3d-general-style" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[24],0,&atx[2]} ,
  {0, "is-on" ,128,0,0,0,0,0,0,0,NULL,&atx[25],NULL,0,&atx[26]} ,
  {301, "BOOLEAN" ,0,1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "style" ,128,1,0,0,0,0,0,0,NULL,&atx[9],NULL,0,&atx[27]} ,
  {0, "color-scheme" ,128,2,0,0,0,0,0,0,NULL,&atx[11],NULL,0,&atx[28]} ,
  {0, "user-color" ,128,3,0,0,0,0,0,0,NULL,&atx[13],NULL,0,NULL} ,
  {0, "nucleotide-sidechains" ,128,3,0,0,0,0,0,0,NULL,&atx[23],NULL,0,&atx[30]} ,
  {0, "heterogens" ,128,4,0,0,0,0,0,0,NULL,&atx[23],NULL,0,&atx[31]} ,
  {0, "solvents" ,128,5,0,0,0,0,0,0,NULL,&atx[23],NULL,0,&atx[32]} ,
  {0, "connections" ,128,6,0,0,0,0,0,0,NULL,&atx[23],NULL,0,&atx[33]} ,
  {0, "helix-objects" ,128,7,0,0,0,0,0,0,NULL,&atx[23],NULL,0,&atx[34]} ,
  {0, "strand-objects" ,128,8,0,0,0,0,0,0,NULL,&atx[23],NULL,0,&atx[35]} ,
  {0, "virtual-disulfides-on" ,128,9,0,0,0,0,0,0,NULL,&atx[25],NULL,0,&atx[36]} ,
  {0, "virtual-disulfide-color" ,128,10,0,0,0,0,0,0,NULL,&atx[13],NULL,0,&atx[37]} ,
  {0, "hydrogens-on" ,128,11,0,0,0,0,0,0,NULL,&atx[25],NULL,0,&atx[38]} ,
  {0, "background-color" ,128,12,0,0,0,0,0,0,NULL,&atx[13],NULL,0,&atx[39]} ,
  {0, "scale-factor" ,128,13,0,0,0,0,0,0,NULL,&atx[15],NULL,0,&atx[40]} ,
  {0, "space-fill-proportion" ,128,14,0,0,0,0,0,0,NULL,&atx[15],NULL,0,&atx[41]} ,
  {0, "ball-radius" ,128,15,0,0,0,0,0,0,NULL,&atx[15],NULL,0,&atx[42]} ,
  {0, "stick-radius" ,128,16,0,0,0,0,0,0,NULL,&atx[15],NULL,0,&atx[43]} ,
  {0, "tube-radius" ,128,17,0,0,0,0,0,0,NULL,&atx[15],NULL,0,&atx[44]} ,
  {0, "tube-worm-radius" ,128,18,0,0,0,0,0,0,NULL,&atx[15],NULL,0,&atx[45]} ,
  {0, "helix-radius" ,128,19,0,0,0,0,0,0,NULL,&atx[15],NULL,0,&atx[46]} ,
  {0, "strand-width" ,128,20,0,0,0,0,0,0,NULL,&atx[15],NULL,0,&atx[47]} ,
  {0, "strand-thickness" ,128,21,0,0,0,0,0,0,NULL,&atx[15],NULL,0,NULL} ,
  {0, "style-table" ,128,1,0,1,0,0,0,0,NULL,&atx[54],&atx[49],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[50],NULL,0,NULL} ,
  {414, "Cn3d-style-table-item" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[51],0,&atx[75]} ,
  {0, "id" ,128,0,0,0,0,0,0,0,NULL,&atx[52],NULL,0,&atx[53]} ,
  {413, "Cn3d-style-table-id" ,1,0,0,0,0,0,0,0,NULL,&atx[15],NULL,0,&atx[50]} ,
  {0, "style" ,128,1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {312, "SEQUENCE OF" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {402, "Cn3d-user-annotations" ,1,0,0,0,0,1,0,0,NULL,&atx[20],&atx[56],0,&atx[67]} ,
  {0, "annotations" ,128,0,0,0,0,0,0,0,NULL,&atx[54],&atx[57],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[58],NULL,0,NULL} ,
  {418, "Cn3d-user-annotation" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[59],0,NULL} ,
  {0, "name" ,128,0,0,0,0,0,0,0,NULL,&atx[60],NULL,0,&atx[61]} ,
  {323, "VisibleString" ,0,26,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "description" ,128,1,0,1,0,0,0,0,NULL,&atx[60],NULL,0,&atx[62]} ,
  {0, "style-id" ,128,2,0,0,0,0,0,0,NULL,&atx[52],NULL,0,&atx[63]} ,
  {0, "residues" ,128,3,0,0,0,0,0,0,NULL,&atx[54],&atx[64],0,&atx[79]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[65],NULL,0,NULL} ,
  {417, "Cn3d-object-location" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[66],0,&atx[58]} ,
  {0, "structure-id" ,128,0,0,0,0,0,0,0,NULL,&atx[67],NULL,0,&atx[68]} ,
  {403, "Biostruc-id" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[72]} ,
  {0, "residues" ,128,1,0,0,0,0,0,0,NULL,&atx[54],&atx[69],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[70],NULL,0,NULL} ,
  {416, "Cn3d-molecule-location" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[71],0,&atx[65]} ,
  {0, "molecule-id" ,128,0,0,0,0,0,0,0,NULL,&atx[72],NULL,0,&atx[73]} ,
  {404, "Molecule-id" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[77]} ,
  {0, "residues" ,128,1,0,1,0,0,0,0,NULL,&atx[54],&atx[74],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[75],NULL,0,NULL} ,
  {415, "Cn3d-residue-range" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[76],0,&atx[70]} ,
  {0, "from" ,128,0,0,0,0,0,0,0,NULL,&atx[77],NULL,0,&atx[78]} ,
  {405, "Residue-id" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[6]} ,
  {0, "to" ,128,1,0,0,0,0,0,0,NULL,&atx[77],NULL,0,NULL} ,
  {0, "is-on" ,128,4,0,0,0,0,0,0,NULL,&atx[25],NULL,0,NULL} };

static AsnModule ampx[1] = {
  { "NCBI-Cn3d" , "cn3d.h14",&atx[0],NULL,NULL,0,0} };

static AsnValxNodePtr avn = avnx;
static AsnTypePtr at = atx;
static AsnModulePtr amp = ampx;



/**************************************************
*
*    Defines for Module NCBI-Cn3d
*
**************************************************/

#define CN3D_STYLE_DICTIONARY &at[0]
#define CN3D_STYLE_DICTIONARY_global_style &at[1]
#define CN3D_STYLE_DICTIONARY_style_table &at[48]
#define CN3D_STYLE_DICTIONARY_style_table_E &at[49]

#define CN3D_USER_ANNOTATIONS &at[55]
#define CN3D_USER_ANNOTATIONS_annotations &at[56]
#define CN3D_USER_ANNOTATIONS_annotations_E &at[57]

#define CN3D_BACKBONE_TYPE &at[6]

#define CN3D_DRAWING_STYLE &at[9]

#define CN3D_COLOR_SCHEME &at[11]

#define CN3D_COLOR &at[13]
#define CN3D_COLOR_scale_factor &at[14]
#define CN3D_COLOR_red &at[16]
#define CN3D_COLOR_green &at[17]
#define CN3D_COLOR_blue &at[18]
#define CN3D_COLOR_alpha &at[19]

#define CN3D_BACKBONE_STYLE &at[4]
#define CN3D_BACKBONE_STYLE_type &at[5]
#define CN3D_BACKBONE_STYLE_style &at[8]
#define CN3D_BACKBONE_STYLE_color_scheme &at[10]
#define CN3D_BACKBONE_STYLE_user_color &at[12]

#define CN3D_GENERAL_STYLE &at[23]
#define CN3D_GENERAL_STYLE_is_on &at[24]
#define CN3D_GENERAL_STYLE_style &at[26]
#define CN3D_GENERAL_STYLE_color_scheme &at[27]
#define CN3D_GENERAL_STYLE_user_color &at[28]

#define CN3D_STYLE_SETTINGS &at[2]
#define CN3D_STYLE_SETTINGS_protein_backbone &at[3]
#define CN3D_STYLE_SETTINGS_nucleotide_backbone &at[21]
#define CN3D_STYLE_SETTINGS_protein_sidechains &at[22]
#define CN3D_STYLE_SETTINGS_nucleotide_sidechains &at[29]
#define CN3D_STYLE_SETTINGS_heterogens &at[30]
#define CN3D_STYLE_SETTINGS_solvents &at[31]
#define CN3D_STYLE_SETTINGS_connections &at[32]
#define CN3D_STYLE_SETTINGS_helix_objects &at[33]
#define CN3D_STYLE_SETTINGS_strand_objects &at[34]
#define CN3D_STYLE_SETTINGS_virtual_disulfides_on &at[35]
#define CN3D_STYLE_SETTINGS_virtual_disulfide_color &at[36]
#define CN3D_STYLE_SETTINGS_hydrogens_on &at[37]
#define CN3D_STYLE_SETTINGS_background_color &at[38]
#define CN3D_STYLE_SETTINGS_scale_factor &at[39]
#define CN3D_STYLE_SETTINGS_space_fill_proportion &at[40]
#define CN3D_STYLE_SETTINGS_ball_radius &at[41]
#define CN3D_STYLE_SETTINGS_stick_radius &at[42]
#define CN3D_STYLE_SETTINGS_tube_radius &at[43]
#define CN3D_STYLE_SETTINGS_tube_worm_radius &at[44]
#define CN3D_STYLE_SETTINGS_helix_radius &at[45]
#define CN3D_STYLE_SETTINGS_strand_width &at[46]
#define CN3D_STYLE_SETTINGS_strand_thickness &at[47]

#define CN3D_STYLE_TABLE_ID &at[52]

#define CN3D_STYLE_TABLE_ITEM &at[50]
#define CN3D_STYLE_TABLE_ITEM_id &at[51]
#define CN3D_STYLE_TABLE_ITEM_style &at[53]

#define CN3D_RESIDUE_RANGE &at[75]
#define CN3D_RESIDUE_RANGE_from &at[76]
#define CN3D_RESIDUE_RANGE_to &at[78]

#define CN3D_MOLECULE_LOCATION &at[70]
#define CN3D_MOLECULE_LOCATION_molecule_id &at[71]
#define CN3D_MOLECULE_LOCATION_residues &at[73]
#define CN3D_MOLECULE_LOCATION_residues_E &at[74]

#define CN3D_OBJECT_LOCATION &at[65]
#define CN3D_OBJECT_LOCATION_structure_id &at[66]
#define CN3D_OBJECT_LOCATION_residues &at[68]
#define CN3D_OBJECT_LOCATION_residues_E &at[69]

#define CN3D_USER_ANNOTATION &at[58]
#define CN3D_USER_ANNOTATION_name &at[59]
#define CN3D_USER_ANNOTATION_description &at[61]
#define CN3D_USER_ANNOTATION_style_id &at[62]
#define CN3D_USER_ANNOTATION_residues &at[63]
#define CN3D_USER_ANNOTATION_residues_E &at[64]
#define CN3D_USER_ANNOTATION_is_on &at[79]
