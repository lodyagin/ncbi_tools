#ifndef _objmmdb3_ 
#define _objmmdb3_ 


#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

#ifdef __cplusplus
extern "C" { /* } */
#endif


/**************************************************
*
*    Generated objects for Module MMDB-features
*    Generated using ASNCODE Revision: 4.2 at Aug 14, 1996  5:26 PM
*
**************************************************/

NLM_EXTERN Boolean LIBCALL
objmmdb3AsnLoad PROTO((void));


/**************************************************
*
*    BiostrucFeatureSet
*
**************************************************/
typedef struct struct_Biostruc_feature_set {
   struct struct_Biostruc_feature_set PNTR next;
   Uint4 OBbits__;
   Int4   id;
   ValNodePtr   descr;
   struct struct_Biostruc_feature PNTR   features;
} BiostrucFeatureSet, PNTR BiostrucFeatureSetPtr;


BiostrucFeatureSetPtr LIBCALL BiostrucFeatureSetFree PROTO ((BiostrucFeatureSetPtr ));
NLM_EXTERN BiostrucFeatureSetPtr LIBCALL BiostrucFeatureSetNew PROTO ((void));
BiostrucFeatureSetPtr LIBCALL BiostrucFeatureSetAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL BiostrucFeatureSetAsnWrite PROTO (( BiostrucFeatureSetPtr , AsnIoPtr, AsnTypePtr));

typedef ValNodePtr ChemGraphPntrsPtr;
typedef ValNode ChemGraphPntrs;
#define ChemGraphPntrs_atoms 1
#define ChemGraphPntrs_residues 2
#define ChemGraphPntrs_molecules 3


ChemGraphPntrsPtr LIBCALL ChemGraphPntrsFree PROTO ((ChemGraphPntrsPtr ));
ChemGraphPntrsPtr LIBCALL ChemGraphPntrsAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL ChemGraphPntrsAsnWrite PROTO (( ChemGraphPntrsPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    AtomPntrs
*
**************************************************/
typedef struct struct_Atom_pntrs {
   Uint4 OBbits__;
   Int4   number_of_ptrs;
   ValNodePtr   molecule_ids;
   ValNodePtr   residue_ids;
   ValNodePtr   atom_ids;
} AtomPntrs, PNTR AtomPntrsPtr;


AtomPntrsPtr LIBCALL AtomPntrsFree PROTO ((AtomPntrsPtr ));
AtomPntrsPtr LIBCALL AtomPntrsNew PROTO (( void ));
AtomPntrsPtr LIBCALL AtomPntrsAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL AtomPntrsAsnWrite PROTO (( AtomPntrsPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ChemGraphAlignment
*
**************************************************/
typedef struct struct_Chem_graph_alignment {
   Uint4 OBbits__;
   Int4   dimension;
   ValNodePtr   biostruc_ids;
   ValNodePtr   alignment;
   ValNodePtr   domain;
   struct struct_Transform PNTR   transform;
   struct struct_Align_stats PNTR   aligndata;
} ChemGraphAlignment, PNTR ChemGraphAlignmentPtr;


ChemGraphAlignmentPtr LIBCALL ChemGraphAlignmentFree PROTO ((ChemGraphAlignmentPtr ));
ChemGraphAlignmentPtr LIBCALL ChemGraphAlignmentNew PROTO (( void ));
ChemGraphAlignmentPtr LIBCALL ChemGraphAlignmentAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL ChemGraphAlignmentAsnWrite PROTO (( ChemGraphAlignmentPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    Sphere
*
**************************************************/
typedef struct struct_Sphere {
   Uint4 OBbits__;
   struct struct_Model_space_point PNTR   center;
   struct struct_RealValue PNTR   radius;
} Sphere, PNTR SpherePtr;


SpherePtr LIBCALL SphereFree PROTO ((SpherePtr ));
SpherePtr LIBCALL SphereNew PROTO (( void ));
SpherePtr LIBCALL SphereAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL SphereAsnWrite PROTO (( SpherePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    Cone
*
**************************************************/
typedef struct struct_Cone {
   Uint4 OBbits__;
   struct struct_Model_space_point PNTR   axis_top;
   struct struct_Model_space_point PNTR   axis_bottom;
   struct struct_RealValue PNTR   radius_bottom;
} Cone, PNTR ConePtr;


ConePtr LIBCALL ConeFree PROTO ((ConePtr ));
ConePtr LIBCALL ConeNew PROTO (( void ));
ConePtr LIBCALL ConeAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL ConeAsnWrite PROTO (( ConePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    Cylinder
*
**************************************************/
typedef struct struct_Cylinder {
   Uint4 OBbits__;
   struct struct_Model_space_point PNTR   axis_top;
   struct struct_Model_space_point PNTR   axis_bottom;
   struct struct_RealValue PNTR   radius;
} Cylinder, PNTR CylinderPtr;


CylinderPtr LIBCALL CylinderFree PROTO ((CylinderPtr ));
CylinderPtr LIBCALL CylinderNew PROTO (( void ));
CylinderPtr LIBCALL CylinderAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL CylinderAsnWrite PROTO (( CylinderPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    Brick
*
**************************************************/
typedef struct struct_Brick {
   Uint4 OBbits__;
   struct struct_Model_space_point PNTR   corner_000;
   struct struct_Model_space_point PNTR   corner_001;
   struct struct_Model_space_point PNTR   corner_010;
   struct struct_Model_space_point PNTR   corner_011;
   struct struct_Model_space_point PNTR   corner_100;
   struct struct_Model_space_point PNTR   corner_101;
   struct struct_Model_space_point PNTR   corner_110;
   struct struct_Model_space_point PNTR   corner_111;
} Brick, PNTR BrickPtr;


BrickPtr LIBCALL BrickFree PROTO ((BrickPtr ));
BrickPtr LIBCALL BrickNew PROTO (( void ));
BrickPtr LIBCALL BrickAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL BrickAsnWrite PROTO (( BrickPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    Transform
*
**************************************************/
typedef struct struct_Transform {
   struct struct_Transform PNTR next;
   Uint4 OBbits__;
   Int4   id;
   ValNodePtr   moves;
} Transform, PNTR TransformPtr;


TransformPtr LIBCALL TransformFree PROTO ((TransformPtr ));
TransformPtr LIBCALL TransformNew PROTO (( void ));
TransformPtr LIBCALL TransformAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL TransformAsnWrite PROTO (( TransformPtr , AsnIoPtr, AsnTypePtr));

typedef ValNodePtr BiostrucFeatureSetDescrPtr;
typedef ValNode BiostrucFeatureSetDescr;
#define BiostrucFeatureSetDescr_name 1
#define BiostrucFeatureSetDescr_pdb_comment 2
#define BiostrucFeatureSetDescr_other_comment 3
#define BiostrucFeatureSetDescr_attribution 4


BiostrucFeatureSetDescrPtr LIBCALL BiostrucFeatureSetDescrFree PROTO ((BiostrucFeatureSetDescrPtr ));
BiostrucFeatureSetDescrPtr LIBCALL BiostrucFeatureSetDescrAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL BiostrucFeatureSetDescrAsnWrite PROTO (( BiostrucFeatureSetDescrPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    BiostrucFeature
*
**************************************************/
typedef struct struct_Biostruc_feature {
   struct struct_Biostruc_feature PNTR next;
   Uint4 OBbits__;
   Int4   id;
   CharPtr   name;
#define OB__Biostruc_feature_type 0

   Int4   type;
   ValNodePtr   Property_property;
   ValNodePtr   Location_location;
} BiostrucFeature, PNTR BiostrucFeaturePtr;


BiostrucFeaturePtr LIBCALL BiostrucFeatureFree PROTO ((BiostrucFeaturePtr ));
NLM_EXTERN BiostrucFeaturePtr LIBCALL BiostrucFeatureNew PROTO (( void ));
BiostrucFeaturePtr LIBCALL BiostrucFeatureAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL BiostrucFeatureAsnWrite PROTO (( BiostrucFeaturePtr , AsnIoPtr, AsnTypePtr));


#ifdef NLM_GENERATED_CODE_PROTO

typedef ValNodePtr Location_locationPtr;
typedef ValNode Location_location;

#endif /* NLM_GENERATED_CODE_PROTO */

#define Location_location_subgraph 1
#define Location_location_region 2
#define Location_location_alignment 3
#define Location_location_similarity 4
#define Location_location_indirect 5

#ifdef NLM_GENERATED_CODE_PROTO

static Location_locationPtr LIBCALL Location_locationFree PROTO ((Location_locationPtr ));
static Location_locationPtr LIBCALL Location_locationAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
static Boolean LIBCALL Location_locationAsnWrite PROTO (( Location_locationPtr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */


#ifdef NLM_GENERATED_CODE_PROTO

typedef ValNodePtr Property_propertyPtr;
typedef ValNode Property_property;

#endif /* NLM_GENERATED_CODE_PROTO */

#define Property_property_color 1
#define Property_property_render 2
#define Property_property_transform 3
#define Property_property_camera 4
#define Property_property_script 5
#define Property_property_user 6

#ifdef NLM_GENERATED_CODE_PROTO

static Property_propertyPtr LIBCALL Property_propertyFree PROTO ((Property_propertyPtr ));
static Property_propertyPtr LIBCALL Property_propertyAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
static Boolean LIBCALL Property_propertyAsnWrite PROTO (( Property_propertyPtr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */



/**************************************************
*
*    ColorProp
*
**************************************************/
typedef struct struct_Color_prop {
   Uint4 OBbits__;
#define OB__Color_prop_r 0

   Int4   r;
#define OB__Color_prop_g 1

   Int4   g;
#define OB__Color_prop_b 2

   Int4   b;
   CharPtr   name;
} ColorProp, PNTR ColorPropPtr;


ColorPropPtr LIBCALL ColorPropFree PROTO ((ColorPropPtr ));
ColorPropPtr LIBCALL ColorPropNew PROTO (( void ));
ColorPropPtr LIBCALL ColorPropAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL ColorPropAsnWrite PROTO (( ColorPropPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    Camera
*
**************************************************/
typedef struct struct_Camera {
   Uint4 OBbits__;
   Int4   mode;
   struct struct_Model_space_point PNTR   x;
   struct struct_Model_space_point PNTR   y;
   struct struct_Model_space_point PNTR   z;
   struct struct_Model_space_point PNTR   up;
   struct struct_Model_space_point PNTR   fore;
   struct struct_Model_space_point PNTR   norm;
   struct struct_Model_space_point PNTR   center;
   struct struct_RealValue PNTR   tooclose;
   struct struct_RealValue PNTR   toofar;
} Camera, PNTR CameraPtr;


CameraPtr LIBCALL CameraFree PROTO ((CameraPtr ));
CameraPtr LIBCALL CameraNew PROTO (( void ));
CameraPtr LIBCALL CameraAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL CameraAsnWrite PROTO (( CameraPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    BiostrucScript
*
**************************************************/
typedef struct struct_BiostrucScriptStep BiostrucScript;
typedef struct struct_BiostrucScriptStep PNTR BiostrucScriptPtr;
#define BiostrucScriptNew() BiostrucScriptStepNew() 

#ifdef NLM_GENERATED_CODE_PROTO

BiostrucScriptPtr LIBCALL BiostrucScriptFree PROTO ((BiostrucScriptPtr ));
BiostrucScriptPtr LIBCALL BiostrucScriptNew PROTO (( void ));
BiostrucScriptPtr LIBCALL BiostrucScriptAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL BiostrucScriptAsnWrite PROTO (( BiostrucScriptPtr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */



/**************************************************
*
*    RegionPntrs
*
**************************************************/
typedef struct struct_Region_pntrs {
   struct struct_Region_pntrs PNTR next;
   Uint4 OBbits__;
   Int4   model_id;
   ValNodePtr   Region_region;
} RegionPntrs, PNTR RegionPntrsPtr;


RegionPntrsPtr LIBCALL RegionPntrsFree PROTO ((RegionPntrsPtr ));
RegionPntrsPtr LIBCALL RegionPntrsNew PROTO (( void ));
RegionPntrsPtr LIBCALL RegionPntrsAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL RegionPntrsAsnWrite PROTO (( RegionPntrsPtr , AsnIoPtr, AsnTypePtr));


#ifdef NLM_GENERATED_CODE_PROTO

typedef ValNodePtr Region_regionPtr;
typedef ValNode Region_region;

#endif /* NLM_GENERATED_CODE_PROTO */

#define Region_region_site 1
#define Region_region_boundary 2

#ifdef NLM_GENERATED_CODE_PROTO

static Region_regionPtr LIBCALL Region_regionFree PROTO ((Region_regionPtr ));
static Region_regionPtr LIBCALL Region_regionAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
static Boolean LIBCALL Region_regionAsnWrite PROTO (( Region_regionPtr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */



/**************************************************
*
*    RegionSimilarity
*
**************************************************/
typedef struct struct_Region_similarity {
   Uint4 OBbits__;
   Int4   dimension;
   ValNodePtr   biostruc_ids;
   struct struct_Region_pntrs PNTR   similarity;
   struct struct_Transform PNTR   transform;
} RegionSimilarity, PNTR RegionSimilarityPtr;


RegionSimilarityPtr LIBCALL RegionSimilarityFree PROTO ((RegionSimilarityPtr ));
RegionSimilarityPtr LIBCALL RegionSimilarityNew PROTO (( void ));
RegionSimilarityPtr LIBCALL RegionSimilarityAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL RegionSimilarityAsnWrite PROTO (( RegionSimilarityPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    OtherFeature
*
**************************************************/
typedef struct struct_Other_feature {
   struct struct_Other_feature PNTR next;
   Uint4 OBbits__;
   ValNodePtr   biostruc_id;
   Int4   set;
   Int4   feature;
} OtherFeature, PNTR OtherFeaturePtr;


OtherFeaturePtr LIBCALL OtherFeatureFree PROTO ((OtherFeaturePtr ));
OtherFeaturePtr LIBCALL OtherFeatureNew PROTO (( void ));
OtherFeaturePtr LIBCALL OtherFeatureAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL OtherFeatureAsnWrite PROTO (( OtherFeaturePtr , AsnIoPtr, AsnTypePtr));

typedef ValNodePtr ResiduePntrsPtr;
typedef ValNode ResiduePntrs;
#define ResiduePntrs_explicit 1
#define ResiduePntrs_interval 2


ResiduePntrsPtr LIBCALL ResiduePntrsFree PROTO ((ResiduePntrsPtr ));
ResiduePntrsPtr LIBCALL ResiduePntrsAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL ResiduePntrsAsnWrite PROTO (( ResiduePntrsPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    MoleculePntrs
*
**************************************************/
typedef struct struct_Molecule_pntrs {
   Uint4 OBbits__;
   Int4   number_of_ptrs;
   ValNodePtr   molecule_ids;
} MoleculePntrs, PNTR MoleculePntrsPtr;


MoleculePntrsPtr LIBCALL MoleculePntrsFree PROTO ((MoleculePntrsPtr ));
MoleculePntrsPtr LIBCALL MoleculePntrsNew PROTO (( void ));
MoleculePntrsPtr LIBCALL MoleculePntrsAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL MoleculePntrsAsnWrite PROTO (( MoleculePntrsPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ResidueExplicitPntrs
*
**************************************************/
typedef struct struct_Residue_explicit_pntrs {
   Uint4 OBbits__;
   Int4   number_of_ptrs;
   ValNodePtr   molecule_ids;
   ValNodePtr   residue_ids;
} ResidueExplicitPntrs, PNTR ResidueExplicitPntrsPtr;


ResidueExplicitPntrsPtr LIBCALL ResidueExplicitPntrsFree PROTO ((ResidueExplicitPntrsPtr ));
ResidueExplicitPntrsPtr LIBCALL ResidueExplicitPntrsNew PROTO (( void ));
ResidueExplicitPntrsPtr LIBCALL ResidueExplicitPntrsAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL ResidueExplicitPntrsAsnWrite PROTO (( ResidueExplicitPntrsPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ResidueIntervalPntr
*
**************************************************/
typedef struct struct_Residue_interval_pntr {
   struct struct_Residue_interval_pntr PNTR next;
   Uint4 OBbits__;
   Int4   molecule_id;
   Int4   from;
   Int4   to;
} ResidueIntervalPntr, PNTR ResidueIntervalPntrPtr;


ResidueIntervalPntrPtr LIBCALL ResidueIntervalPntrFree PROTO ((ResidueIntervalPntrPtr ));
ResidueIntervalPntrPtr LIBCALL ResidueIntervalPntrNew PROTO (( void ));
ResidueIntervalPntrPtr LIBCALL ResidueIntervalPntrAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL ResidueIntervalPntrAsnWrite PROTO (( ResidueIntervalPntrPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    RegionCoordinates
*
**************************************************/
typedef struct struct_Region_coordinates {
   struct struct_Region_coordinates PNTR next;
   Uint4 OBbits__;
   Int4   model_coord_set_id;
#define OB__Region_coordinates_number_of_coords 0

   Int4   number_of_coords;
   ValNodePtr   coordinate_indices;
} RegionCoordinates, PNTR RegionCoordinatesPtr;


RegionCoordinatesPtr LIBCALL RegionCoordinatesFree PROTO ((RegionCoordinatesPtr ));
RegionCoordinatesPtr LIBCALL RegionCoordinatesNew PROTO (( void ));
RegionCoordinatesPtr LIBCALL RegionCoordinatesAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL RegionCoordinatesAsnWrite PROTO (( RegionCoordinatesPtr , AsnIoPtr, AsnTypePtr));

typedef ValNodePtr RegionBoundaryPtr;
typedef ValNode RegionBoundary;
#define RegionBoundary_sphere 1
#define RegionBoundary_cone 2
#define RegionBoundary_cylinder 3
#define RegionBoundary_brick 4


RegionBoundaryPtr LIBCALL RegionBoundaryFree PROTO ((RegionBoundaryPtr ));
RegionBoundaryPtr LIBCALL RegionBoundaryAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL RegionBoundaryAsnWrite PROTO (( RegionBoundaryPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    AlignStats
*
**************************************************/
typedef struct struct_Align_stats {
   struct struct_Align_stats PNTR next;
   Uint4 OBbits__;
   CharPtr   descr;
#define OB__Align_stats_scale_factor 0

   Int4   scale_factor;
#define OB__Align_stats_vast_score 1

   Int4   vast_score;
#define OB__Align_stats_vast_mlogp 2

   Int4   vast_mlogp;
#define OB__Align_stats_align_res 3

   Int4   align_res;
#define OB__Align_stats_rmsd 4

   Int4   rmsd;
#define OB__Align_stats_blast_score 5

   Int4   blast_score;
#define OB__Align_stats_blast_mlogp 6

   Int4   blast_mlogp;
#define OB__Align_stats_other_score 7

   Int4   other_score;
} AlignStats, PNTR AlignStatsPtr;


AlignStatsPtr LIBCALL AlignStatsFree PROTO ((AlignStatsPtr ));
AlignStatsPtr LIBCALL AlignStatsNew PROTO (( void ));
AlignStatsPtr LIBCALL AlignStatsAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL AlignStatsAsnWrite PROTO (( AlignStatsPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ModelSpacePoint
*
**************************************************/
typedef struct struct_Model_space_point {
   Uint4 OBbits__;
   Int4   scale_factor;
   Int4   x;
   Int4   y;
   Int4   z;
} ModelSpacePoint, PNTR ModelSpacePointPtr;


ModelSpacePointPtr LIBCALL ModelSpacePointFree PROTO ((ModelSpacePointPtr ));
ModelSpacePointPtr LIBCALL ModelSpacePointNew PROTO (( void ));
ModelSpacePointPtr LIBCALL ModelSpacePointAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL ModelSpacePointAsnWrite PROTO (( ModelSpacePointPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    RealValue
*
**************************************************/
typedef struct struct_RealValue {
   Uint4 OBbits__;
   Int4   scale_factor;
   Int4   scaled_integer_value;
} RealValue, PNTR RealValuePtr;


RealValuePtr LIBCALL RealValueFree PROTO ((RealValuePtr ));
RealValuePtr LIBCALL RealValueNew PROTO (( void ));
RealValuePtr LIBCALL RealValueAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL RealValueAsnWrite PROTO (( RealValuePtr , AsnIoPtr, AsnTypePtr));

typedef ValNodePtr MovePtr;
typedef ValNode Move;
#define Move_rotate 1
#define Move_translate 2


MovePtr LIBCALL MoveFree PROTO ((MovePtr ));
MovePtr LIBCALL MoveAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL MoveAsnWrite PROTO (( MovePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    RotMatrix
*
**************************************************/
typedef struct struct_Rot_matrix {
   Uint4 OBbits__;
   Int4   scale_factor;
   Int4   rot_11;
   Int4   rot_12;
   Int4   rot_13;
   Int4   rot_21;
   Int4   rot_22;
   Int4   rot_23;
   Int4   rot_31;
   Int4   rot_32;
   Int4   rot_33;
} RotMatrix, PNTR RotMatrixPtr;


RotMatrixPtr LIBCALL RotMatrixFree PROTO ((RotMatrixPtr ));
RotMatrixPtr LIBCALL RotMatrixNew PROTO (( void ));
RotMatrixPtr LIBCALL RotMatrixAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL RotMatrixAsnWrite PROTO (( RotMatrixPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    TransMatrix
*
**************************************************/
typedef struct struct_Trans_matrix {
   Uint4 OBbits__;
   Int4   scale_factor;
   Int4   tran_1;
   Int4   tran_2;
   Int4   tran_3;
} TransMatrix, PNTR TransMatrixPtr;


TransMatrixPtr LIBCALL TransMatrixFree PROTO ((TransMatrixPtr ));
TransMatrixPtr LIBCALL TransMatrixNew PROTO (( void ));
TransMatrixPtr LIBCALL TransMatrixAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL TransMatrixAsnWrite PROTO (( TransMatrixPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    BiostrucScriptStep
*
**************************************************/
typedef struct struct_Biostruc_script_step {
   struct struct_Biostruc_script_step PNTR next;
   Uint4 OBbits__;
   Int4   step_id;
   CharPtr   step_name;
   struct struct_Other_feature PNTR   feature_do;
   struct struct_Transform PNTR   camera_move;
   Int4   pause;
   Uint1   waitevent;
   Int4   extra;
#define OB__Biostruc_script_step_jump 0

   Int4   jump;
} BiostrucScriptStep, PNTR BiostrucScriptStepPtr;


BiostrucScriptStepPtr LIBCALL BiostrucScriptStepFree PROTO ((BiostrucScriptStepPtr ));
BiostrucScriptStepPtr LIBCALL BiostrucScriptStepNew PROTO (( void ));
BiostrucScriptStepPtr LIBCALL BiostrucScriptStepAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL BiostrucScriptStepAsnWrite PROTO (( BiostrucScriptStepPtr , AsnIoPtr, AsnTypePtr));

#ifdef __cplusplus
/* { */ }
#endif

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif /* _objmmdb3_ */
