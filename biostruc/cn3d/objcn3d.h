#ifndef _objcn3d_ 
#define _objcn3d_ 

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
*    Generated objects for Module NCBI-Cn3d
*    Generated using ASNCODE Revision: 6.10 at Jun 21, 2001 10:32 AM
*
**************************************************/

NLM_EXTERN Boolean LIBCALL
objcn3dAsnLoad PROTO((void));


/**************************************************
*
*    Cn3dStyleDictionary
*
**************************************************/
typedef struct struct_Cn3d_style_dictionary {
   struct struct_Cn3d_style_settings PNTR   global_style;
   struct struct_Cn3d_style_table_item PNTR   style_table;
} Cn3dStyleDictionary, PNTR Cn3dStyleDictionaryPtr;


NLM_EXTERN Cn3dStyleDictionaryPtr LIBCALL Cn3dStyleDictionaryFree PROTO ((Cn3dStyleDictionaryPtr ));
NLM_EXTERN Cn3dStyleDictionaryPtr LIBCALL Cn3dStyleDictionaryNew PROTO (( void ));
NLM_EXTERN Cn3dStyleDictionaryPtr LIBCALL Cn3dStyleDictionaryAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL Cn3dStyleDictionaryAsnWrite PROTO (( Cn3dStyleDictionaryPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    Cn3dUserAnnotations
*
**************************************************/
typedef struct struct_Cn3d_user_annotations {
   struct struct_Cn3d_user_annotation PNTR   annotations;
} Cn3dUserAnnotations, PNTR Cn3dUserAnnotationsPtr;


NLM_EXTERN Cn3dUserAnnotationsPtr LIBCALL Cn3dUserAnnotationsFree PROTO ((Cn3dUserAnnotationsPtr ));
NLM_EXTERN Cn3dUserAnnotationsPtr LIBCALL Cn3dUserAnnotationsNew PROTO (( void ));
NLM_EXTERN Cn3dUserAnnotationsPtr LIBCALL Cn3dUserAnnotationsAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL Cn3dUserAnnotationsAsnWrite PROTO (( Cn3dUserAnnotationsPtr , AsnIoPtr, AsnTypePtr));

/* following #defines are for enumerated type, not used by object loaders */
#define Cn3d_backbone_type_off 1
#define Cn3d_backbone_type_trace 2
#define Cn3d_backbone_type_partial 3
#define Cn3d_backbone_type_complete 4

/* following #defines are for enumerated type, not used by object loaders */
#define Cn3d_drawing_style_wire 1
#define Cn3d_drawing_style_tubes 2
#define Cn3d_drawing_style_ball_and_stick 3
#define Cn3d_drawing_style_space_fill 4
#define Cn3d_drawing_style_wire_worm 5
#define Cn3d_drawing_style_tube_worm 6
#define Cn3d_drawing_style_with_arrows 7
#define Cn3d_drawing_style_without_arrows 8

/* following #defines are for enumerated type, not used by object loaders */
#define Cn3d_color_scheme_element 1
#define Cn3d_color_scheme_object 2
#define Cn3d_color_scheme_molecule 3
#define Cn3d_color_scheme_domain 4
#define Cn3d_color_scheme_secondary_structure 5
#define Cn3d_color_scheme_user_select 6
#define Cn3d_color_scheme_aligned 7
#define Cn3d_color_scheme_identity 8
#define Cn3d_color_scheme_variety 9
#define Cn3d_color_scheme_weighted_variety 10
#define Cn3d_color_scheme_information_content 11
#define Cn3d_color_scheme_fit 12



/**************************************************
*
*    Cn3dColor
*
**************************************************/
typedef struct struct_Cn3d_color {
   Int4   scale_factor;
   Int4   red;
   Int4   green;
   Int4   blue;
   Int4   alpha;
} Cn3dColor, PNTR Cn3dColorPtr;


NLM_EXTERN Cn3dColorPtr LIBCALL Cn3dColorFree PROTO ((Cn3dColorPtr ));
NLM_EXTERN Cn3dColorPtr LIBCALL Cn3dColorNew PROTO (( void ));
NLM_EXTERN Cn3dColorPtr LIBCALL Cn3dColorAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL Cn3dColorAsnWrite PROTO (( Cn3dColorPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    Cn3dBackboneStyle
*
**************************************************/
typedef struct struct_Cn3d_backbone_style {
   Uint2   type;
   Uint2   style;
   Uint2   color_scheme;
   struct struct_Cn3d_color PNTR   user_color;
} Cn3dBackboneStyle, PNTR Cn3dBackboneStylePtr;


NLM_EXTERN Cn3dBackboneStylePtr LIBCALL Cn3dBackboneStyleFree PROTO ((Cn3dBackboneStylePtr ));
NLM_EXTERN Cn3dBackboneStylePtr LIBCALL Cn3dBackboneStyleNew PROTO (( void ));
NLM_EXTERN Cn3dBackboneStylePtr LIBCALL Cn3dBackboneStyleAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL Cn3dBackboneStyleAsnWrite PROTO (( Cn3dBackboneStylePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    Cn3dGeneralStyle
*
**************************************************/
typedef struct struct_Cn3d_general_style {
   Uint1   is_on;
   Uint2   style;
   Uint2   color_scheme;
   struct struct_Cn3d_color PNTR   user_color;
} Cn3dGeneralStyle, PNTR Cn3dGeneralStylePtr;


NLM_EXTERN Cn3dGeneralStylePtr LIBCALL Cn3dGeneralStyleFree PROTO ((Cn3dGeneralStylePtr ));
NLM_EXTERN Cn3dGeneralStylePtr LIBCALL Cn3dGeneralStyleNew PROTO (( void ));
NLM_EXTERN Cn3dGeneralStylePtr LIBCALL Cn3dGeneralStyleAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL Cn3dGeneralStyleAsnWrite PROTO (( Cn3dGeneralStylePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    Cn3dStyleSettings
*
**************************************************/
typedef struct struct_Cn3d_style_settings {
   struct struct_Cn3d_backbone_style PNTR   protein_backbone;
   struct struct_Cn3d_backbone_style PNTR   nucleotide_backbone;
   struct struct_Cn3d_general_style PNTR   protein_sidechains;
   struct struct_Cn3d_general_style PNTR   nucleotide_sidechains;
   struct struct_Cn3d_general_style PNTR   heterogens;
   struct struct_Cn3d_general_style PNTR   solvents;
   struct struct_Cn3d_general_style PNTR   connections;
   struct struct_Cn3d_general_style PNTR   helix_objects;
   struct struct_Cn3d_general_style PNTR   strand_objects;
   Uint1   virtual_disulfides_on;
   struct struct_Cn3d_color PNTR   virtual_disulfide_color;
   Uint1   hydrogens_on;
   struct struct_Cn3d_color PNTR   background_color;
   Int4   scale_factor;
   Int4   space_fill_proportion;
   Int4   ball_radius;
   Int4   stick_radius;
   Int4   tube_radius;
   Int4   tube_worm_radius;
   Int4   helix_radius;
   Int4   strand_width;
   Int4   strand_thickness;
} Cn3dStyleSettings, PNTR Cn3dStyleSettingsPtr;


NLM_EXTERN Cn3dStyleSettingsPtr LIBCALL Cn3dStyleSettingsFree PROTO ((Cn3dStyleSettingsPtr ));
NLM_EXTERN Cn3dStyleSettingsPtr LIBCALL Cn3dStyleSettingsNew PROTO (( void ));
NLM_EXTERN Cn3dStyleSettingsPtr LIBCALL Cn3dStyleSettingsAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL Cn3dStyleSettingsAsnWrite PROTO (( Cn3dStyleSettingsPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    Cn3dStyleTableItem
*
**************************************************/
typedef struct struct_Cn3d_style_table_item {
   struct struct_Cn3d_style_table_item PNTR next;
   Int4   id;
   struct struct_Cn3d_style_settings PNTR   style;
} Cn3dStyleTableItem, PNTR Cn3dStyleTableItemPtr;


NLM_EXTERN Cn3dStyleTableItemPtr LIBCALL Cn3dStyleTableItemFree PROTO ((Cn3dStyleTableItemPtr ));
NLM_EXTERN Cn3dStyleTableItemPtr LIBCALL Cn3dStyleTableItemNew PROTO (( void ));
NLM_EXTERN Cn3dStyleTableItemPtr LIBCALL Cn3dStyleTableItemAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL Cn3dStyleTableItemAsnWrite PROTO (( Cn3dStyleTableItemPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    Cn3dResidueRange
*
**************************************************/
typedef struct struct_Cn3d_residue_range {
   struct struct_Cn3d_residue_range PNTR next;
   Int4   from;
   Int4   to;
} Cn3dResidueRange, PNTR Cn3dResidueRangePtr;


NLM_EXTERN Cn3dResidueRangePtr LIBCALL Cn3dResidueRangeFree PROTO ((Cn3dResidueRangePtr ));
NLM_EXTERN Cn3dResidueRangePtr LIBCALL Cn3dResidueRangeNew PROTO (( void ));
NLM_EXTERN Cn3dResidueRangePtr LIBCALL Cn3dResidueRangeAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL Cn3dResidueRangeAsnWrite PROTO (( Cn3dResidueRangePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    Cn3dMoleculeLocation
*
**************************************************/
typedef struct struct_Cn3d_molecule_location {
   struct struct_Cn3d_molecule_location PNTR next;
   Int4   molecule_id;
   struct struct_Cn3d_residue_range PNTR   residues;
} Cn3dMoleculeLocation, PNTR Cn3dMoleculeLocationPtr;


NLM_EXTERN Cn3dMoleculeLocationPtr LIBCALL Cn3dMoleculeLocationFree PROTO ((Cn3dMoleculeLocationPtr ));
NLM_EXTERN Cn3dMoleculeLocationPtr LIBCALL Cn3dMoleculeLocationNew PROTO (( void ));
NLM_EXTERN Cn3dMoleculeLocationPtr LIBCALL Cn3dMoleculeLocationAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL Cn3dMoleculeLocationAsnWrite PROTO (( Cn3dMoleculeLocationPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    Cn3dObjectLocation
*
**************************************************/
typedef struct struct_Cn3d_object_location {
   struct struct_Cn3d_object_location PNTR next;
   ValNodePtr   structure_id;
   struct struct_Cn3d_molecule_location PNTR   residues;
} Cn3dObjectLocation, PNTR Cn3dObjectLocationPtr;


NLM_EXTERN Cn3dObjectLocationPtr LIBCALL Cn3dObjectLocationFree PROTO ((Cn3dObjectLocationPtr ));
NLM_EXTERN Cn3dObjectLocationPtr LIBCALL Cn3dObjectLocationNew PROTO (( void ));
NLM_EXTERN Cn3dObjectLocationPtr LIBCALL Cn3dObjectLocationAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL Cn3dObjectLocationAsnWrite PROTO (( Cn3dObjectLocationPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    Cn3dUserAnnotation
*
**************************************************/
typedef struct struct_Cn3d_user_annotation {
   struct struct_Cn3d_user_annotation PNTR next;
   CharPtr   name;
   CharPtr   description;
   Int4   style_id;
   struct struct_Cn3d_object_location PNTR   residues;
   Uint1   is_on;
} Cn3dUserAnnotation, PNTR Cn3dUserAnnotationPtr;


NLM_EXTERN Cn3dUserAnnotationPtr LIBCALL Cn3dUserAnnotationFree PROTO ((Cn3dUserAnnotationPtr ));
NLM_EXTERN Cn3dUserAnnotationPtr LIBCALL Cn3dUserAnnotationNew PROTO (( void ));
NLM_EXTERN Cn3dUserAnnotationPtr LIBCALL Cn3dUserAnnotationAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL Cn3dUserAnnotationAsnWrite PROTO (( Cn3dUserAnnotationPtr , AsnIoPtr, AsnTypePtr));

#ifdef __cplusplus
/* { */ }
#endif

#endif /* _objcn3d_ */

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

