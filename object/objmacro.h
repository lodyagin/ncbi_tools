#ifndef _objmacro_ 
#define _objmacro_ 

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
*    Generated objects for Module NCBI-Macro
*    Generated using ASNCODE Revision: 6.16 at Feb 22, 2008 12:57 PM
*
**************************************************/

NLM_EXTERN Boolean LIBCALL
objmacroAsnLoad PROTO((void));


/**************************************************
*
*    AECRAction
*
**************************************************/
typedef struct struct_AECR_action {
   ValNodePtr   action;
   ValNodePtr   constraint;
} AECRAction, PNTR AECRActionPtr;


NLM_EXTERN AECRActionPtr LIBCALL AECRActionFree PROTO ((AECRActionPtr ));
NLM_EXTERN AECRActionPtr LIBCALL AECRActionNew PROTO (( void ));
NLM_EXTERN AECRActionPtr LIBCALL AECRActionAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL AECRActionAsnWrite PROTO (( AECRActionPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ParseAction
*
**************************************************/
typedef struct struct_Parse_action {
   struct struct_Text_portion PNTR   portion;
   ValNodePtr   src;
   ValNodePtr   dest;
   Uint2   capitalization;
   Uint1   remove_from_parsed;
   Uint2   existing_text;
} ParseAction, PNTR ParseActionPtr;


NLM_EXTERN ParseActionPtr LIBCALL ParseActionFree PROTO ((ParseActionPtr ));
NLM_EXTERN ParseActionPtr LIBCALL ParseActionNew PROTO (( void ));
NLM_EXTERN ParseActionPtr LIBCALL ParseActionAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ParseActionAsnWrite PROTO (( ParseActionPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    MacroActionList
*
**************************************************/
typedef ValNode MacroActionList;
typedef ValNodePtr MacroActionListPtr;
#define MacroActionListNew() ValNodeNew(NULL) 

#ifdef NLM_GENERATED_CODE_PROTO

NLM_EXTERN MacroActionListPtr LIBCALL MacroActionListFree PROTO ((MacroActionListPtr ));
NLM_EXTERN MacroActionListPtr LIBCALL MacroActionListNew PROTO (( void ));
NLM_EXTERN MacroActionListPtr LIBCALL MacroActionListAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL MacroActionListAsnWrite PROTO (( MacroActionListPtr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */

/* following #defines are for enumerated type, not used by object loaders */
#define String_location_contains 1
#define String_location_equals 2
#define String_location_starts 3
#define String_location_ends 4
#define String_location_inlist 5



/**************************************************
*
*    StringConstraint
*
**************************************************/
typedef struct struct_String_constraint {
   CharPtr   match_text;
   Uint2   match_location;
   Uint1   case_sensitive;
   Uint1   whole_word;
   Uint1   not_present;
} StringConstraint, PNTR StringConstraintPtr;


NLM_EXTERN StringConstraintPtr LIBCALL StringConstraintFree PROTO ((StringConstraintPtr ));
NLM_EXTERN StringConstraintPtr LIBCALL StringConstraintNew PROTO (( void ));
NLM_EXTERN StringConstraintPtr LIBCALL StringConstraintAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL StringConstraintAsnWrite PROTO (( StringConstraintPtr , AsnIoPtr, AsnTypePtr));

/* following #defines are for enumerated type, not used by object loaders */
#define Strand_constraint_any 0
#define Strand_constraint_plus 1
#define Strand_constraint_minus 2

/* following #defines are for enumerated type, not used by object loaders */
#define Seqtype_constraint_any 0
#define Seqtype_constraint_nuc 1
#define Seqtype_constraint_prot 2



/**************************************************
*
*    LocationConstraint
*
**************************************************/
typedef struct struct_Location_constraint {
   Uint2   strand;
   Uint2   seq_type;
} LocationConstraint, PNTR LocationConstraintPtr;


NLM_EXTERN LocationConstraintPtr LIBCALL LocationConstraintFree PROTO ((LocationConstraintPtr ));
NLM_EXTERN LocationConstraintPtr LIBCALL LocationConstraintNew PROTO (( void ));
NLM_EXTERN LocationConstraintPtr LIBCALL LocationConstraintAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL LocationConstraintAsnWrite PROTO (( LocationConstraintPtr , AsnIoPtr, AsnTypePtr));

/* following #defines are for enumerated type, not used by object loaders */
#define Object_type_constraint_any 0
#define Object_type_constraint_feature 1
#define Object_type_constraint_descriptor 2

/* following #defines are for enumerated type, not used by object loaders */
#define Feature_type_any 0
#define Feature_type_gene 1
#define Feature_type_org 2
#define Feature_type_cds 3
#define Feature_type_prot 4
#define Feature_type_preRNA 5
#define Feature_type_mRNA 6
#define Feature_type_tRNA 7
#define Feature_type_rRNA 8
#define Feature_type_snRNA 9
#define Feature_type_scRNA 10
#define Feature_type_otherRNA 11
#define Feature_type_pub 12
#define Feature_type_seq 13
#define Feature_type_imp 14
#define Feature_type_allele 15
#define Feature_type_attenuator 16
#define Feature_type_c_region 17
#define Feature_type_caat_signal 18
#define Feature_type_imp_CDS 19
#define Feature_type_conflict 20
#define Feature_type_d_loop 21
#define Feature_type_d_segment 22
#define Feature_type_enhancer 23
#define Feature_type_exon 24
#define Feature_type_gC_signal 25
#define Feature_type_iDNA 26
#define Feature_type_intron 27
#define Feature_type_j_segment 28
#define Feature_type_ltr 29
#define Feature_type_mat_peptide 30
#define Feature_type_misc_binding 31
#define Feature_type_misc_difference 32
#define Feature_type_misc_feature 33
#define Feature_type_misc_recomb 34
#define Feature_type_misc_RNA 35
#define Feature_type_misc_signal 36
#define Feature_type_misc_structure 37
#define Feature_type_modified_base 38
#define Feature_type_mutation 39
#define Feature_type_n_region 40
#define Feature_type_old_sequence 41
#define Feature_type_polyA_signal 42
#define Feature_type_polyA_site 43
#define Feature_type_precursor_RNA 44
#define Feature_type_prim_transcript 45
#define Feature_type_primer_bind 46
#define Feature_type_promoter 47
#define Feature_type_protein_bind 48
#define Feature_type_rbs 49
#define Feature_type_repeat_region 50
#define Feature_type_repeat_unit 51
#define Feature_type_rep_origin 52
#define Feature_type_s_region 53
#define Feature_type_satellite 54
#define Feature_type_sig_peptide 55
#define Feature_type_source 56
#define Feature_type_stem_loop 57
#define Feature_type_sts 58
#define Feature_type_tata_signal 59
#define Feature_type_terminator 60
#define Feature_type_transit_peptide 61
#define Feature_type_unsure 62
#define Feature_type_v_region 63
#define Feature_type_v_segment 64
#define Feature_type_variation 65
#define Feature_type_virion 66
#define Feature_type_n3clip 67
#define Feature_type_n3UTR 68
#define Feature_type_n5clip 69
#define Feature_type_n5UTR 70
#define Feature_type_n10_signal 71
#define Feature_type_n35_signal 72
#define Feature_type_site_ref 73
#define Feature_type_region 74
#define Feature_type_comment 75
#define Feature_type_bond 76
#define Feature_type_site 77
#define Feature_type_rsite 78
#define Feature_type_user 79
#define Feature_type_txinit 80
#define Feature_type_num 81
#define Feature_type_psec_str 82
#define Feature_type_non_std_residue 83
#define Feature_type_het 84
#define Feature_type_biosrc 85
#define Feature_type_preprotein 86
#define Feature_type_mat_peptide_aa 87
#define Feature_type_sig_peptide_aa 88
#define Feature_type_transit_peptide_aa 89
#define Feature_type_snoRNA 90
#define Feature_type_gap 91
#define Feature_type_operon 92
#define Feature_type_oriT 93
#define Feature_type_ncRNA 94
#define Feature_type_tmRNA 95

/* following #defines are for enumerated type, not used by object loaders */
#define Feat_qual_legal_allele 1
#define Feat_qual_legal_activity 2
#define Feat_qual_legal_anticodon 3
#define Feat_qual_legal_bound_moiety 4
#define Feat_qual_legal_chromosome 5
#define Feat_qual_legal_citation 6
#define Feat_qual_legal_codon 7
#define Feat_qual_legal_codon_start 8
#define Feat_qual_legal_codons_recognized 9
#define Feat_qual_legal_compare 10
#define Feat_qual_legal_cons_splice 11
#define Feat_qual_legal_db_xref 12
#define Feat_qual_legal_description 13
#define Feat_qual_legal_direction 14
#define Feat_qual_legal_ec_number 15
#define Feat_qual_legal_environmental_sample 16
#define Feat_qual_legal_evidence 17
#define Feat_qual_legal_exception 18
#define Feat_qual_legal_experiment 19
#define Feat_qual_legal_focus 20
#define Feat_qual_legal_frequency 21
#define Feat_qual_legal_function 22
#define Feat_qual_legal_gene 23
#define Feat_qual_legal_gene_description 24
#define Feat_qual_legal_inference 25
#define Feat_qual_legal_label 26
#define Feat_qual_legal_locus_tag 27
#define Feat_qual_legal_map 28
#define Feat_qual_legal_mobile_element 29
#define Feat_qual_legal_mod_base 30
#define Feat_qual_legal_mol_type 31
#define Feat_qual_legal_ncRNA_class 32
#define Feat_qual_legal_note 33
#define Feat_qual_legal_number 34
#define Feat_qual_legal_old_locus_tag 35
#define Feat_qual_legal_operon 36
#define Feat_qual_legal_organism 37
#define Feat_qual_legal_organelle 38
#define Feat_qual_legal_partial 39
#define Feat_qual_legal_phenotype 40
#define Feat_qual_legal_plasmid 41
#define Feat_qual_legal_product 42
#define Feat_qual_legal_protein_id 43
#define Feat_qual_legal_pseudo 44
#define Feat_qual_legal_rearranged 45
#define Feat_qual_legal_replace 46
#define Feat_qual_legal_rpt_family 47
#define Feat_qual_legal_rpt_type 48
#define Feat_qual_legal_rpt_unit 49
#define Feat_qual_legal_rpt_unit_seq 50
#define Feat_qual_legal_rpt_unit_range 51
#define Feat_qual_legal_segment 52
#define Feat_qual_legal_sequenced_mol 53
#define Feat_qual_legal_standard_name 54
#define Feat_qual_legal_synonym 55
#define Feat_qual_legal_transcript_id 56
#define Feat_qual_legal_transgenic 57
#define Feat_qual_legal_translation 58
#define Feat_qual_legal_transl_except 59
#define Feat_qual_legal_transl_table 60
#define Feat_qual_legal_usedin 61



/**************************************************
*
*    FeatQualLegalVal
*
**************************************************/
typedef struct struct_Feat_qual_legal_val {
   Uint2   qual;
   CharPtr   val;
} FeatQualLegalVal, PNTR FeatQualLegalValPtr;


NLM_EXTERN FeatQualLegalValPtr LIBCALL FeatQualLegalValFree PROTO ((FeatQualLegalValPtr ));
NLM_EXTERN FeatQualLegalValPtr LIBCALL FeatQualLegalValNew PROTO (( void ));
NLM_EXTERN FeatQualLegalValPtr LIBCALL FeatQualLegalValAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL FeatQualLegalValAsnWrite PROTO (( FeatQualLegalValPtr , AsnIoPtr, AsnTypePtr));

typedef ValNodePtr FeatQualLegalValChoicePtr;
typedef ValNode FeatQualLegalValChoice;
#define FeatQualLegalValChoice_qual 1


NLM_EXTERN FeatQualLegalValChoicePtr LIBCALL FeatQualLegalValChoiceFree PROTO ((FeatQualLegalValChoicePtr ));
NLM_EXTERN FeatQualLegalValChoicePtr LIBCALL FeatQualLegalValChoiceAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL FeatQualLegalValChoiceAsnWrite PROTO (( FeatQualLegalValChoicePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    FeatQualLegalSet
*
**************************************************/
typedef ValNode FeatQualLegalSet;
typedef ValNodePtr FeatQualLegalSetPtr;
#define FeatQualLegalSetNew() ValNodeNew(NULL) 

#ifdef NLM_GENERATED_CODE_PROTO

NLM_EXTERN FeatQualLegalSetPtr LIBCALL FeatQualLegalSetFree PROTO ((FeatQualLegalSetPtr ));
NLM_EXTERN FeatQualLegalSetPtr LIBCALL FeatQualLegalSetNew PROTO (( void ));
NLM_EXTERN FeatQualLegalSetPtr LIBCALL FeatQualLegalSetAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL FeatQualLegalSetAsnWrite PROTO (( FeatQualLegalSetPtr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */

typedef ValNodePtr FeatQualChoicePtr;
typedef ValNode FeatQualChoice;
#define FeatQualChoice_legal_qual 1
#define FeatQualChoice_illegal_qual 2


NLM_EXTERN FeatQualChoicePtr LIBCALL FeatQualChoiceFree PROTO ((FeatQualChoicePtr ));
NLM_EXTERN FeatQualChoicePtr LIBCALL FeatQualChoiceAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL FeatQualChoiceAsnWrite PROTO (( FeatQualChoicePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    FeatureField
*
**************************************************/
typedef struct struct_Feature_field {
   Uint2   type;
   ValNodePtr   field;
} FeatureField, PNTR FeatureFieldPtr;


NLM_EXTERN FeatureFieldPtr LIBCALL FeatureFieldFree PROTO ((FeatureFieldPtr ));
NLM_EXTERN FeatureFieldPtr LIBCALL FeatureFieldNew PROTO (( void ));
NLM_EXTERN FeatureFieldPtr LIBCALL FeatureFieldAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL FeatureFieldAsnWrite PROTO (( FeatureFieldPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    FeatureFieldLegal
*
**************************************************/
typedef struct struct_Feature_field_legal {
   Uint2   type;
   Uint2   field;
} FeatureFieldLegal, PNTR FeatureFieldLegalPtr;


NLM_EXTERN FeatureFieldLegalPtr LIBCALL FeatureFieldLegalFree PROTO ((FeatureFieldLegalPtr ));
NLM_EXTERN FeatureFieldLegalPtr LIBCALL FeatureFieldLegalNew PROTO (( void ));
NLM_EXTERN FeatureFieldLegalPtr LIBCALL FeatureFieldLegalAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL FeatureFieldLegalAsnWrite PROTO (( FeatureFieldLegalPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    FeatureFieldPair
*
**************************************************/
typedef struct struct_Feature_field_pair {
   Uint2   type;
   ValNodePtr   field_from;
   ValNodePtr   field_to;
} FeatureFieldPair, PNTR FeatureFieldPairPtr;


NLM_EXTERN FeatureFieldPairPtr LIBCALL FeatureFieldPairFree PROTO ((FeatureFieldPairPtr ));
NLM_EXTERN FeatureFieldPairPtr LIBCALL FeatureFieldPairNew PROTO (( void ));
NLM_EXTERN FeatureFieldPairPtr LIBCALL FeatureFieldPairAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL FeatureFieldPairAsnWrite PROTO (( FeatureFieldPairPtr , AsnIoPtr, AsnTypePtr));

/* following #defines are for enumerated type, not used by object loaders */
#define Source_qual_acronym 1
#define Source_qual_anamorph 2
#define Source_qual_authority 3
#define Source_qual_bio_material 4
#define Source_qual_biotype 5
#define Source_qual_biovar 6
#define Source_qual_breed 7
#define Source_qual_cell_line 8
#define Source_qual_cell_type 9
#define Source_qual_chemovar 10
#define Source_qual_chromosome 11
#define Source_qual_clone 12
#define Source_qual_clone_lib 13
#define Source_qual_collected_by 14
#define Source_qual_collection_date 15
#define Source_qual_common 16
#define Source_qual_common_name 17
#define Source_qual_country 18
#define Source_qual_cultivar 19
#define Source_qual_culture_collection 20
#define Source_qual_dev_stage 21
#define Source_qual_division 22
#define Source_qual_dosage 23
#define Source_qual_ecotype 24
#define Source_qual_endogenous_virus_name 25
#define Source_qual_environmental_sample 26
#define Source_qual_forma 27
#define Source_qual_forma_specialis 28
#define Source_qual_frequency 29
#define Source_qual_fwd_primer_name 30
#define Source_qual_fwd_primer_seq 31
#define Source_qual_gb_acronym 32
#define Source_qual_gb_anamorph 33
#define Source_qual_gb_synonym 34
#define Source_qual_genotype 35
#define Source_qual_germline 36
#define Source_qual_group 37
#define Source_qual_haplotype 38
#define Source_qual_identified_by 39
#define Source_qual_insertion_seq_name 40
#define Source_qual_isolate 41
#define Source_qual_isolation_source 42
#define Source_qual_lab_host 43
#define Source_qual_lat_lon 44
#define Source_qual_lineage 45
#define Source_qual_map 46
#define Source_qual_metagenome_source 47
#define Source_qual_metagenomic 48
#define Source_qual_old_lineage 49
#define Source_qual_old_name 50
#define Source_qual_orgmod_note 51
#define Source_qual_nat_host 52
#define Source_qual_pathovar 53
#define Source_qual_plasmid_name 54
#define Source_qual_plastid_name 55
#define Source_qual_pop_variant 56
#define Source_qual_rearranged 57
#define Source_qual_rev_primer_name 58
#define Source_qual_rev_primer_seq 59
#define Source_qual_segment 60
#define Source_qual_serogroup 61
#define Source_qual_serotype 62
#define Source_qual_serovar 63
#define Source_qual_sex 64
#define Source_qual_specimen_voucher 65
#define Source_qual_strain 66
#define Source_qual_subclone 67
#define Source_qual_subgroup 68
#define Source_qual_subsource_note 69
#define Source_qual_sub_species 70
#define Source_qual_substrain 71
#define Source_qual_subtype 72
#define Source_qual_synonym 73
#define Source_qual_taxname 74
#define Source_qual_teleomorph 75
#define Source_qual_tissue_lib 76
#define Source_qual_tissue_type 77
#define Source_qual_transgenic 78
#define Source_qual_transposon_name 79
#define Source_qual_type 80
#define Source_qual_variety 81



/**************************************************
*
*    SourceQualPair
*
**************************************************/
typedef struct struct_Source_qual_pair {
   Uint2   field_from;
   Uint2   field_to;
} SourceQualPair, PNTR SourceQualPairPtr;


NLM_EXTERN SourceQualPairPtr LIBCALL SourceQualPairFree PROTO ((SourceQualPairPtr ));
NLM_EXTERN SourceQualPairPtr LIBCALL SourceQualPairNew PROTO (( void ));
NLM_EXTERN SourceQualPairPtr LIBCALL SourceQualPairAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL SourceQualPairAsnWrite PROTO (( SourceQualPairPtr , AsnIoPtr, AsnTypePtr));

/* following #defines are for enumerated type, not used by object loaders */
#define Source_location_unknown 0
#define Source_location_genomic 1
#define Source_location_chloroplast 2
#define Source_location_chromoplast 3
#define Source_location_kinetoplast 4
#define Source_location_mitochondrion 5
#define Source_location_plastid 6
#define Source_location_macronuclear 7
#define Source_location_extrachrom 8
#define Source_location_plasmid 9
#define Source_location_transposon 10
#define Source_location_insertion_seq 11
#define Source_location_cyanelle 12
#define Source_location_proviral 13
#define Source_location_virion 14
#define Source_location_nucleomorph 15
#define Source_location_apicoplast 16
#define Source_location_leucoplast 17
#define Source_location_proplastid 18
#define Source_location_endogenous_virus 19
#define Source_location_hydrogenosome 20
#define Source_location_chromosome 21
#define Source_location_chromatophore 22

/* following #defines are for enumerated type, not used by object loaders */
#define Source_origin_unknown 0
#define Source_origin_natural 1
#define Source_origin_natmut 2
#define Source_origin_mut 3
#define Source_origin_artificial 4
#define Source_origin_synthetic 5
#define Source_origin_other 255

typedef ValNodePtr SourceQualChoicePtr;
typedef ValNode SourceQualChoice;
#define SourceQualChoice_textqual 1
#define SourceQualChoice_location 2
#define SourceQualChoice_origin 3


NLM_EXTERN SourceQualChoicePtr LIBCALL SourceQualChoiceFree PROTO ((SourceQualChoicePtr ));
NLM_EXTERN SourceQualChoicePtr LIBCALL SourceQualChoiceAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL SourceQualChoiceAsnWrite PROTO (( SourceQualChoicePtr , AsnIoPtr, AsnTypePtr));

/* following #defines are for enumerated type, not used by object loaders */
#define CDSGeneProt_field_cds_comment 1
#define CDSGeneProt_field_gene_locus 2
#define CDSGeneProt_field_gene_description 3
#define CDSGeneProt_field_gene_comment 4
#define CDSGeneProt_field_gene_allele 5
#define CDSGeneProt_field_gene_maploc 6
#define CDSGeneProt_field_gene_locus_tag 7
#define CDSGeneProt_field_gene_synonym 8
#define CDSGeneProt_field_gene_old_locus_tag 9
#define CDSGeneProt_field_mrna_product 10
#define CDSGeneProt_field_mrna_comment 11
#define CDSGeneProt_field_prot_name 12
#define CDSGeneProt_field_prot_description 13
#define CDSGeneProt_field_prot_ec_number 14
#define CDSGeneProt_field_prot_activity 15
#define CDSGeneProt_field_prot_comment 16
#define CDSGeneProt_field_mat_peptide_name 17
#define CDSGeneProt_field_mat_peptide_description 18
#define CDSGeneProt_field_mat_peptide_ec_number 19
#define CDSGeneProt_field_mat_peptide_activity 20
#define CDSGeneProt_field_mat_peptide_comment 21



/**************************************************
*
*    CDSGeneProtFieldPair
*
**************************************************/
typedef struct struct_CDSGeneProt_field_pair {
   Uint2   field_from;
   Uint2   field_to;
} CDSGeneProtFieldPair, PNTR CDSGeneProtFieldPairPtr;


NLM_EXTERN CDSGeneProtFieldPairPtr LIBCALL CDSGeneProtFieldPairFree PROTO ((CDSGeneProtFieldPairPtr ));
NLM_EXTERN CDSGeneProtFieldPairPtr LIBCALL CDSGeneProtFieldPairNew PROTO (( void ));
NLM_EXTERN CDSGeneProtFieldPairPtr LIBCALL CDSGeneProtFieldPairAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL CDSGeneProtFieldPairAsnWrite PROTO (( CDSGeneProtFieldPairPtr , AsnIoPtr, AsnTypePtr));

/* following #defines are for enumerated type, not used by object loaders */
#define Molecule_type_unknown 0
#define Molecule_type_genomic 1
#define Molecule_type_precursor_RNA 2
#define Molecule_type_mRNA 3
#define Molecule_type_rRNA 4
#define Molecule_type_tRNA 5
#define Molecule_type_snRNA 6
#define Molecule_type_scRNA 7
#define Molecule_type_genomic_mRNA 8
#define Molecule_type_cRNA 9
#define Molecule_type_snoRNA 10
#define Molecule_type_transcribed_RNA 11
#define Molecule_type_ncRNA 12
#define Molecule_type_transfer_messenger_RNA 13
#define Molecule_type_other 14

/* following #defines are for enumerated type, not used by object loaders */
#define Technique_type_unknown 0
#define Technique_type_standard 1
#define Technique_type_est 2
#define Technique_type_sts 3
#define Technique_type_survey 4
#define Technique_type_genetic_map 5
#define Technique_type_physical_map 6
#define Technique_type_derived 7
#define Technique_type_concept_trans 8
#define Technique_type_seq_pept 9
#define Technique_type_both 10
#define Technique_type_seq_pept_overlap 11
#define Technique_type_seq_pept_homol 12
#define Technique_type_concept_trans_a 13
#define Technique_type_htgs_1 14
#define Technique_type_htgs_2 15
#define Technique_type_htgs_3 16
#define Technique_type_fli_cDNA 17
#define Technique_type_htgs_0 18
#define Technique_type_htc 19
#define Technique_type_wgs 20
#define Technique_type_barcode 21
#define Technique_type_composite_wgs_htgs 22
#define Technique_type_tsa 23
#define Technique_type_other 24

/* following #defines are for enumerated type, not used by object loaders */
#define Completedness_type_unknown 0
#define Completedness_type_complete 1
#define Completedness_type_partial 2
#define Completedness_type_no_left 3
#define Completedness_type_no_right 4
#define Completedness_type_no_ends 5
#define Completedness_type_has_left 6
#define Completedness_type_has_right 7
#define Completedness_type_other 6

/* following #defines are for enumerated type, not used by object loaders */
#define Molecule_class_type_unknown 0
#define Molecule_class_type_dna 1
#define Molecule_class_type_rna 2
#define Molecule_class_type_protein 3
#define Molecule_class_type_nucleotide 4
#define Molecule_class_type_other 5

/* following #defines are for enumerated type, not used by object loaders */
#define Topology_type_unknown 0
#define Topology_type_linear 1
#define Topology_type_circular 2
#define Topology_type_tandem 3
#define Topology_type_other 4

/* following #defines are for enumerated type, not used by object loaders */
#define Strand_type_unknown 0
#define Strand_type_single 1
#define Strand_type_double__ 2
#define Strand_type_mixed 3
#define Strand_type_mixed_rev 4
#define Strand_type_other 5

typedef ValNodePtr MolinfoFieldPtr;
typedef ValNode MolinfoField;
#define MolinfoField_molecule 1
#define MolinfoField_technique 2
#define MolinfoField_completedness 3
#define MolinfoField_mol_class 4
#define MolinfoField_topology 5
#define MolinfoField_strand 6


NLM_EXTERN MolinfoFieldPtr LIBCALL MolinfoFieldFree PROTO ((MolinfoFieldPtr ));
NLM_EXTERN MolinfoFieldPtr LIBCALL MolinfoFieldAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL MolinfoFieldAsnWrite PROTO (( MolinfoFieldPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    MolinfoMoleculePair
*
**************************************************/
typedef struct struct_Molinfo_molecule_pair {
   Uint2   from;
   Uint2   to;
} MolinfoMoleculePair, PNTR MolinfoMoleculePairPtr;


NLM_EXTERN MolinfoMoleculePairPtr LIBCALL MolinfoMoleculePairFree PROTO ((MolinfoMoleculePairPtr ));
NLM_EXTERN MolinfoMoleculePairPtr LIBCALL MolinfoMoleculePairNew PROTO (( void ));
NLM_EXTERN MolinfoMoleculePairPtr LIBCALL MolinfoMoleculePairAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL MolinfoMoleculePairAsnWrite PROTO (( MolinfoMoleculePairPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    MolinfoTechniquePair
*
**************************************************/
typedef struct struct_Molinfo_technique_pair {
   Uint2   from;
   Uint2   to;
} MolinfoTechniquePair, PNTR MolinfoTechniquePairPtr;


NLM_EXTERN MolinfoTechniquePairPtr LIBCALL MolinfoTechniquePairFree PROTO ((MolinfoTechniquePairPtr ));
NLM_EXTERN MolinfoTechniquePairPtr LIBCALL MolinfoTechniquePairNew PROTO (( void ));
NLM_EXTERN MolinfoTechniquePairPtr LIBCALL MolinfoTechniquePairAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL MolinfoTechniquePairAsnWrite PROTO (( MolinfoTechniquePairPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    MolinfoCompletednessPair
*
**************************************************/
typedef struct struct_Molinfo_completedness_pair {
   Uint2   from;
   Uint2   to;
} MolinfoCompletednessPair, PNTR MolinfoCompletednessPairPtr;


NLM_EXTERN MolinfoCompletednessPairPtr LIBCALL MolinfoCompletednessPairFree PROTO ((MolinfoCompletednessPairPtr ));
NLM_EXTERN MolinfoCompletednessPairPtr LIBCALL MolinfoCompletednessPairNew PROTO (( void ));
NLM_EXTERN MolinfoCompletednessPairPtr LIBCALL MolinfoCompletednessPairAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL MolinfoCompletednessPairAsnWrite PROTO (( MolinfoCompletednessPairPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    MolinfoMolClassPair
*
**************************************************/
typedef struct struct_Molinfo_mol_class_pair {
   Uint2   from;
   Uint2   to;
} MolinfoMolClassPair, PNTR MolinfoMolClassPairPtr;


NLM_EXTERN MolinfoMolClassPairPtr LIBCALL MolinfoMolClassPairFree PROTO ((MolinfoMolClassPairPtr ));
NLM_EXTERN MolinfoMolClassPairPtr LIBCALL MolinfoMolClassPairNew PROTO (( void ));
NLM_EXTERN MolinfoMolClassPairPtr LIBCALL MolinfoMolClassPairAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL MolinfoMolClassPairAsnWrite PROTO (( MolinfoMolClassPairPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    MolinfoTopologyPair
*
**************************************************/
typedef struct struct_Molinfo_topology_pair {
   Uint2   from;
   Uint2   to;
} MolinfoTopologyPair, PNTR MolinfoTopologyPairPtr;


NLM_EXTERN MolinfoTopologyPairPtr LIBCALL MolinfoTopologyPairFree PROTO ((MolinfoTopologyPairPtr ));
NLM_EXTERN MolinfoTopologyPairPtr LIBCALL MolinfoTopologyPairNew PROTO (( void ));
NLM_EXTERN MolinfoTopologyPairPtr LIBCALL MolinfoTopologyPairAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL MolinfoTopologyPairAsnWrite PROTO (( MolinfoTopologyPairPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    MolinfoStrandPair
*
**************************************************/
typedef struct struct_Molinfo_strand_pair {
   Uint2   from;
   Uint2   to;
} MolinfoStrandPair, PNTR MolinfoStrandPairPtr;


NLM_EXTERN MolinfoStrandPairPtr LIBCALL MolinfoStrandPairFree PROTO ((MolinfoStrandPairPtr ));
NLM_EXTERN MolinfoStrandPairPtr LIBCALL MolinfoStrandPairNew PROTO (( void ));
NLM_EXTERN MolinfoStrandPairPtr LIBCALL MolinfoStrandPairAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL MolinfoStrandPairAsnWrite PROTO (( MolinfoStrandPairPtr , AsnIoPtr, AsnTypePtr));

typedef ValNodePtr MolinfoFieldPairPtr;
typedef ValNode MolinfoFieldPair;
#define MolinfoFieldPair_molecule 1
#define MolinfoFieldPair_technique 2
#define MolinfoFieldPair_completedness 3
#define MolinfoFieldPair_mol_class 4
#define MolinfoFieldPair_topology 5
#define MolinfoFieldPair_strand 6


NLM_EXTERN MolinfoFieldPairPtr LIBCALL MolinfoFieldPairFree PROTO ((MolinfoFieldPairPtr ));
NLM_EXTERN MolinfoFieldPairPtr LIBCALL MolinfoFieldPairAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL MolinfoFieldPairAsnWrite PROTO (( MolinfoFieldPairPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    SourceConstraint
*
**************************************************/
typedef struct struct_Source_constraint {
   ValNodePtr   field1;
   ValNodePtr   field2;
   struct struct_String_constraint PNTR   constraint;
   Uint2   type_constraint;
} SourceConstraint, PNTR SourceConstraintPtr;


NLM_EXTERN SourceConstraintPtr LIBCALL SourceConstraintFree PROTO ((SourceConstraintPtr ));
NLM_EXTERN SourceConstraintPtr LIBCALL SourceConstraintNew PROTO (( void ));
NLM_EXTERN SourceConstraintPtr LIBCALL SourceConstraintAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL SourceConstraintAsnWrite PROTO (( SourceConstraintPtr , AsnIoPtr, AsnTypePtr));

/* following #defines are for enumerated type, not used by object loaders */
#define CDSGeneProt_feature_type_constraint_gene 1
#define CDSGeneProt_feature_type_constraint_mRNA 2
#define CDSGeneProt_feature_type_constraint_cds 3
#define CDSGeneProt_feature_type_constraint_prot 4
#define CDSGeneProt_feature_type_constraint_exon 5
#define CDSGeneProt_feature_type_constraint_mat_peptide 6



/**************************************************
*
*    CDSGeneProtPseudoConstraint
*
**************************************************/
typedef struct struct_CDSGeneProt_pseudo_constraint {
   Uint2   feature;
   Uint1   is_pseudo;
} CDSGeneProtPseudoConstraint, PNTR CDSGeneProtPseudoConstraintPtr;


NLM_EXTERN CDSGeneProtPseudoConstraintPtr LIBCALL CDSGeneProtPseudoConstraintFree PROTO ((CDSGeneProtPseudoConstraintPtr ));
NLM_EXTERN CDSGeneProtPseudoConstraintPtr LIBCALL CDSGeneProtPseudoConstraintNew PROTO (( void ));
NLM_EXTERN CDSGeneProtPseudoConstraintPtr LIBCALL CDSGeneProtPseudoConstraintAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL CDSGeneProtPseudoConstraintAsnWrite PROTO (( CDSGeneProtPseudoConstraintPtr , AsnIoPtr, AsnTypePtr));

typedef ValNodePtr CDSGeneProtConstraintFieldPtr;
typedef ValNode CDSGeneProtConstraintField;
#define CDSGeneProtConstraintField_field 1


NLM_EXTERN CDSGeneProtConstraintFieldPtr LIBCALL CDSGeneProtConstraintFieldFree PROTO ((CDSGeneProtConstraintFieldPtr ));
NLM_EXTERN CDSGeneProtConstraintFieldPtr LIBCALL CDSGeneProtConstraintFieldAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL CDSGeneProtConstraintFieldAsnWrite PROTO (( CDSGeneProtConstraintFieldPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    CDSGeneProtQualConstraint
*
**************************************************/
typedef struct struct_CDSGeneProt_qual_constraint {
   ValNodePtr   field1;
   ValNodePtr   field2;
   struct struct_String_constraint PNTR   constraint;
} CDSGeneProtQualConstraint, PNTR CDSGeneProtQualConstraintPtr;


NLM_EXTERN CDSGeneProtQualConstraintPtr LIBCALL CDSGeneProtQualConstraintFree PROTO ((CDSGeneProtQualConstraintPtr ));
NLM_EXTERN CDSGeneProtQualConstraintPtr LIBCALL CDSGeneProtQualConstraintNew PROTO (( void ));
NLM_EXTERN CDSGeneProtQualConstraintPtr LIBCALL CDSGeneProtQualConstraintAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL CDSGeneProtQualConstraintAsnWrite PROTO (( CDSGeneProtQualConstraintPtr , AsnIoPtr, AsnTypePtr));

/* following #defines are for enumerated type, not used by object loaders */
#define Sequence_constraint_rnamol_any 0
#define Sequence_constraint_rnamol_genomic 1
#define Sequence_constraint_rnamol_precursor_RNA 2
#define Sequence_constraint_rnamol_mRNA 3
#define Sequence_constraint_rnamol_rRNA 4
#define Sequence_constraint_rnamol_tRNA 5
#define Sequence_constraint_rnamol_snRNA 6
#define Sequence_constraint_rnamol_scRNA 7
#define Sequence_constraint_rnamol_genomic_mRNA 8
#define Sequence_constraint_rnamol_cRNA 9
#define Sequence_constraint_rnamol_snoRNA 10
#define Sequence_constraint_rnamol_transcribed_RNA 11
#define Sequence_constraint_rnamol_ncRNA 12
#define Sequence_constraint_rnamol_transfer_messenger_RNA 13

typedef ValNodePtr SequenceConstraintMolTypeConstraintPtr;
typedef ValNode SequenceConstraintMolTypeConstraint;
#define SequenceConstraintMolTypeConstraint_any 1
#define SequenceConstraintMolTypeConstraint_nucleotide 2
#define SequenceConstraintMolTypeConstraint_dna 3
#define SequenceConstraintMolTypeConstraint_rna 4
#define SequenceConstraintMolTypeConstraint_protein 5


NLM_EXTERN SequenceConstraintMolTypeConstraintPtr LIBCALL SequenceConstraintMolTypeConstraintFree PROTO ((SequenceConstraintMolTypeConstraintPtr ));
NLM_EXTERN SequenceConstraintMolTypeConstraintPtr LIBCALL SequenceConstraintMolTypeConstraintAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL SequenceConstraintMolTypeConstraintAsnWrite PROTO (( SequenceConstraintMolTypeConstraintPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    SequenceConstraint
*
**************************************************/
typedef struct struct_Sequence_constraint {
   ValNodePtr   seqtype;
   struct struct_String_constraint PNTR   id;
   Uint2   feature;
} SequenceConstraint, PNTR SequenceConstraintPtr;


NLM_EXTERN SequenceConstraintPtr LIBCALL SequenceConstraintFree PROTO ((SequenceConstraintPtr ));
NLM_EXTERN SequenceConstraintPtr LIBCALL SequenceConstraintNew PROTO (( void ));
NLM_EXTERN SequenceConstraintPtr LIBCALL SequenceConstraintAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL SequenceConstraintAsnWrite PROTO (( SequenceConstraintPtr , AsnIoPtr, AsnTypePtr));

typedef ValNodePtr ConstraintChoicePtr;
typedef ValNode ConstraintChoice;
#define ConstraintChoice_string 1
#define ConstraintChoice_location 2
#define ConstraintChoice_source 3
#define ConstraintChoice_cdsgeneprot_qual 4
#define ConstraintChoice_cdsgeneprot_pseudo 5
#define ConstraintChoice_sequence 6


NLM_EXTERN ConstraintChoicePtr LIBCALL ConstraintChoiceFree PROTO ((ConstraintChoicePtr ));
NLM_EXTERN ConstraintChoicePtr LIBCALL ConstraintChoiceAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ConstraintChoiceAsnWrite PROTO (( ConstraintChoicePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ConstraintChoiceSet
*
**************************************************/
typedef ValNode ConstraintChoiceSet;
typedef ValNodePtr ConstraintChoiceSetPtr;
#define ConstraintChoiceSetNew() ValNodeNew(NULL) 

#ifdef NLM_GENERATED_CODE_PROTO

NLM_EXTERN ConstraintChoiceSetPtr LIBCALL ConstraintChoiceSetFree PROTO ((ConstraintChoiceSetPtr ));
NLM_EXTERN ConstraintChoiceSetPtr LIBCALL ConstraintChoiceSetNew PROTO (( void ));
NLM_EXTERN ConstraintChoiceSetPtr LIBCALL ConstraintChoiceSetAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ConstraintChoiceSetAsnWrite PROTO (( ConstraintChoiceSetPtr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */



/**************************************************
*
*    TextPortion
*
**************************************************/
typedef struct struct_Text_portion {
   CharPtr   left_text;
   Uint1   include_left;
   CharPtr   right_text;
   Uint1   include_right;
   Uint1   inside;
   Uint1   case_sensitive;
   Uint1   whole_word;
} TextPortion, PNTR TextPortionPtr;


NLM_EXTERN TextPortionPtr LIBCALL TextPortionFree PROTO ((TextPortionPtr ));
NLM_EXTERN TextPortionPtr LIBCALL TextPortionNew PROTO (( void ));
NLM_EXTERN TextPortionPtr LIBCALL TextPortionAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL TextPortionAsnWrite PROTO (( TextPortionPtr , AsnIoPtr, AsnTypePtr));

/* following #defines are for enumerated type, not used by object loaders */
#define Field_edit_location_anywhere 0
#define Field_edit_location_beginning 1
#define Field_edit_location_end 2



/**************************************************
*
*    FieldEdit
*
**************************************************/
typedef struct struct_Field_edit {
   CharPtr   find_txt;
   CharPtr   repl_txt;
   Uint2   location;
} FieldEdit, PNTR FieldEditPtr;


NLM_EXTERN FieldEditPtr LIBCALL FieldEditFree PROTO ((FieldEditPtr ));
NLM_EXTERN FieldEditPtr LIBCALL FieldEditNew PROTO (( void ));
NLM_EXTERN FieldEditPtr LIBCALL FieldEditAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL FieldEditAsnWrite PROTO (( FieldEditPtr , AsnIoPtr, AsnTypePtr));

typedef ValNodePtr FieldTypePtr;
typedef ValNode FieldType;
#define FieldType_source_qual 1
#define FieldType_feature_field 2
#define FieldType_cds_gene_prot 3
#define FieldType_molinfo_field 4


NLM_EXTERN FieldTypePtr LIBCALL FieldTypeFree PROTO ((FieldTypePtr ));
NLM_EXTERN FieldTypePtr LIBCALL FieldTypeAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL FieldTypeAsnWrite PROTO (( FieldTypePtr , AsnIoPtr, AsnTypePtr));

typedef ValNodePtr FieldPairTypePtr;
typedef ValNode FieldPairType;
#define FieldPairType_source_qual 1
#define FieldPairType_feature_field 2
#define FieldPairType_cds_gene_prot 3
#define FieldPairType_molinfo_field 4


NLM_EXTERN FieldPairTypePtr LIBCALL FieldPairTypeFree PROTO ((FieldPairTypePtr ));
NLM_EXTERN FieldPairTypePtr LIBCALL FieldPairTypeAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL FieldPairTypeAsnWrite PROTO (( FieldPairTypePtr , AsnIoPtr, AsnTypePtr));

/* following #defines are for enumerated type, not used by object loaders */
#define ExistingTextOption_replace_old 1
#define ExistingTextOption_append_semi 2
#define ExistingTextOption_append_space 3
#define ExistingTextOption_append_colon 4
#define ExistingTextOption_append_none 5
#define ExistingTextOption_prefix_semi 6
#define ExistingTextOption_prefix_space 7
#define ExistingTextOption_prefix_colon 8
#define ExistingTextOption_prefix_none 9
#define ExistingTextOption_leave_old 10



/**************************************************
*
*    ApplyAction
*
**************************************************/
typedef struct struct_Apply_action {
   ValNodePtr   field;
   CharPtr   value;
   Uint2   existing_text;
} ApplyAction, PNTR ApplyActionPtr;


NLM_EXTERN ApplyActionPtr LIBCALL ApplyActionFree PROTO ((ApplyActionPtr ));
NLM_EXTERN ApplyActionPtr LIBCALL ApplyActionNew PROTO (( void ));
NLM_EXTERN ApplyActionPtr LIBCALL ApplyActionAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ApplyActionAsnWrite PROTO (( ApplyActionPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    EditAction
*
**************************************************/
typedef struct struct_Edit_action {
   struct struct_Field_edit PNTR   edit;
   ValNodePtr   field;
} EditAction, PNTR EditActionPtr;


NLM_EXTERN EditActionPtr LIBCALL EditActionFree PROTO ((EditActionPtr ));
NLM_EXTERN EditActionPtr LIBCALL EditActionNew PROTO (( void ));
NLM_EXTERN EditActionPtr LIBCALL EditActionAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL EditActionAsnWrite PROTO (( EditActionPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ConvertAction
*
**************************************************/
typedef struct struct_Convert_action {
   ValNodePtr   fields;
   Uint2   existing_text;
} ConvertAction, PNTR ConvertActionPtr;


NLM_EXTERN ConvertActionPtr LIBCALL ConvertActionFree PROTO ((ConvertActionPtr ));
NLM_EXTERN ConvertActionPtr LIBCALL ConvertActionNew PROTO (( void ));
NLM_EXTERN ConvertActionPtr LIBCALL ConvertActionAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ConvertActionAsnWrite PROTO (( ConvertActionPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    CopyAction
*
**************************************************/
typedef struct struct_Copy_action {
   ValNodePtr   fields;
   Uint2   existing_text;
} CopyAction, PNTR CopyActionPtr;


NLM_EXTERN CopyActionPtr LIBCALL CopyActionFree PROTO ((CopyActionPtr ));
NLM_EXTERN CopyActionPtr LIBCALL CopyActionNew PROTO (( void ));
NLM_EXTERN CopyActionPtr LIBCALL CopyActionAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL CopyActionAsnWrite PROTO (( CopyActionPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    SwapAction
*
**************************************************/
typedef struct struct_Swap_action {
   ValNodePtr   fields;
   ValNodePtr   field_to;
} SwapAction, PNTR SwapActionPtr;


NLM_EXTERN SwapActionPtr LIBCALL SwapActionFree PROTO ((SwapActionPtr ));
NLM_EXTERN SwapActionPtr LIBCALL SwapActionNew PROTO (( void ));
NLM_EXTERN SwapActionPtr LIBCALL SwapActionAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL SwapActionAsnWrite PROTO (( SwapActionPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    AECRParseAction
*
**************************************************/
typedef struct struct_AECRParse_action {
   struct struct_Text_portion PNTR   portion;
   ValNodePtr   fields;
   Uint1   remove_from_parsed;
   Uint2   existing_text;
} AECRParseAction, PNTR AECRParseActionPtr;


NLM_EXTERN AECRParseActionPtr LIBCALL AECRParseActionFree PROTO ((AECRParseActionPtr ));
NLM_EXTERN AECRParseActionPtr LIBCALL AECRParseActionNew PROTO (( void ));
NLM_EXTERN AECRParseActionPtr LIBCALL AECRParseActionAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL AECRParseActionAsnWrite PROTO (( AECRParseActionPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    RemoveAction
*
**************************************************/
typedef struct struct_Remove_action {
   ValNodePtr   field;
} RemoveAction, PNTR RemoveActionPtr;


NLM_EXTERN RemoveActionPtr LIBCALL RemoveActionFree PROTO ((RemoveActionPtr ));
NLM_EXTERN RemoveActionPtr LIBCALL RemoveActionNew PROTO (( void ));
NLM_EXTERN RemoveActionPtr LIBCALL RemoveActionAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL RemoveActionAsnWrite PROTO (( RemoveActionPtr , AsnIoPtr, AsnTypePtr));

typedef ValNodePtr ActionChoicePtr;
typedef ValNode ActionChoice;
#define ActionChoice_apply 1
#define ActionChoice_edit 2
#define ActionChoice_convert 3
#define ActionChoice_copy 4
#define ActionChoice_swap 5
#define ActionChoice_remove 6
#define ActionChoice_parse 7


NLM_EXTERN ActionChoicePtr LIBCALL ActionChoiceFree PROTO ((ActionChoicePtr ));
NLM_EXTERN ActionChoicePtr LIBCALL ActionChoiceAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ActionChoiceAsnWrite PROTO (( ActionChoicePtr , AsnIoPtr, AsnTypePtr));

/* following #defines are for enumerated type, not used by object loaders */
#define Cap_change_none 0
#define Cap_change_tolower 1
#define Cap_change_toupper 2
#define Cap_change_firstcap 3

typedef ValNodePtr ParseSrcOrgChoicePtr;
typedef ValNode ParseSrcOrgChoice;
#define ParseSrcOrgChoice_source_qual 1
#define ParseSrcOrgChoice_taxname_after_binomial 2


NLM_EXTERN ParseSrcOrgChoicePtr LIBCALL ParseSrcOrgChoiceFree PROTO ((ParseSrcOrgChoicePtr ));
NLM_EXTERN ParseSrcOrgChoicePtr LIBCALL ParseSrcOrgChoiceAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ParseSrcOrgChoiceAsnWrite PROTO (( ParseSrcOrgChoicePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ParseSrcOrg
*
**************************************************/
typedef struct struct_Parse_src_org {
   ValNodePtr   field;
   Uint2   type;
} ParseSrcOrg, PNTR ParseSrcOrgPtr;


NLM_EXTERN ParseSrcOrgPtr LIBCALL ParseSrcOrgFree PROTO ((ParseSrcOrgPtr ));
NLM_EXTERN ParseSrcOrgPtr LIBCALL ParseSrcOrgNew PROTO (( void ));
NLM_EXTERN ParseSrcOrgPtr LIBCALL ParseSrcOrgAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ParseSrcOrgAsnWrite PROTO (( ParseSrcOrgPtr , AsnIoPtr, AsnTypePtr));

typedef ValNodePtr ParseSrcPtr;
typedef ValNode ParseSrc;
#define ParseSrc_defline 1
#define ParseSrc_flatfile 2
#define ParseSrc_local_id 3
#define ParseSrc_org 4
#define ParseSrc_comment 5
#define ParseSrc_bankit_comment 6
#define ParseSrc_structured_comment 7
#define ParseSrc_file_id 8


NLM_EXTERN ParseSrcPtr LIBCALL ParseSrcFree PROTO ((ParseSrcPtr ));
NLM_EXTERN ParseSrcPtr LIBCALL ParseSrcAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ParseSrcAsnWrite PROTO (( ParseSrcPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ParseDstOrg
*
**************************************************/
typedef struct struct_Parse_dst_org {
   ValNodePtr   field;
   Uint2   type;
} ParseDstOrg, PNTR ParseDstOrgPtr;


NLM_EXTERN ParseDstOrgPtr LIBCALL ParseDstOrgFree PROTO ((ParseDstOrgPtr ));
NLM_EXTERN ParseDstOrgPtr LIBCALL ParseDstOrgNew PROTO (( void ));
NLM_EXTERN ParseDstOrgPtr LIBCALL ParseDstOrgAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ParseDstOrgAsnWrite PROTO (( ParseDstOrgPtr , AsnIoPtr, AsnTypePtr));

typedef ValNodePtr ParseDestPtr;
typedef ValNode ParseDest;
#define ParseDest_defline 1
#define ParseDest_org 2
#define ParseDest_featqual 3
#define ParseDest_dbxref 4


NLM_EXTERN ParseDestPtr LIBCALL ParseDestFree PROTO ((ParseDestPtr ));
NLM_EXTERN ParseDestPtr LIBCALL ParseDestAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ParseDestAsnWrite PROTO (( ParseDestPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    LocationInterval
*
**************************************************/
typedef struct struct_Location_interval {
   Int4   from;
   Int4   to;
} LocationInterval, PNTR LocationIntervalPtr;


NLM_EXTERN LocationIntervalPtr LIBCALL LocationIntervalFree PROTO ((LocationIntervalPtr ));
NLM_EXTERN LocationIntervalPtr LIBCALL LocationIntervalNew PROTO (( void ));
NLM_EXTERN LocationIntervalPtr LIBCALL LocationIntervalAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL LocationIntervalAsnWrite PROTO (( LocationIntervalPtr , AsnIoPtr, AsnTypePtr));

typedef ValNodePtr LocationChoicePtr;
typedef ValNode LocationChoice;
#define LocationChoice_interval 1
#define LocationChoice_whole_sequence 2


NLM_EXTERN LocationChoicePtr LIBCALL LocationChoiceFree PROTO ((LocationChoicePtr ));
NLM_EXTERN LocationChoicePtr LIBCALL LocationChoiceAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL LocationChoiceAsnWrite PROTO (( LocationChoicePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    SequenceList
*
**************************************************/
typedef ValNode SequenceList;
typedef ValNodePtr SequenceListPtr;
#define SequenceListNew() ValNodeNew(NULL) 

#ifdef NLM_GENERATED_CODE_PROTO

NLM_EXTERN SequenceListPtr LIBCALL SequenceListFree PROTO ((SequenceListPtr ));
NLM_EXTERN SequenceListPtr LIBCALL SequenceListNew PROTO (( void ));
NLM_EXTERN SequenceListPtr LIBCALL SequenceListAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL SequenceListAsnWrite PROTO (( SequenceListPtr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */

typedef ValNodePtr SequenceListChoicePtr;
typedef ValNode SequenceListChoice;
#define SequenceListChoice_list 1
#define SequenceListChoice_all 2


NLM_EXTERN SequenceListChoicePtr LIBCALL SequenceListChoiceFree PROTO ((SequenceListChoicePtr ));
NLM_EXTERN SequenceListChoicePtr LIBCALL SequenceListChoiceAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL SequenceListChoiceAsnWrite PROTO (( SequenceListChoicePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    ApplyFeatureAction
*
**************************************************/
typedef struct struct_Apply_feature_action {
   Uint2   type;
   Uint1   partial5;
   Uint1   partial3;
   Uint1   plus_strand;
   ValNodePtr   location;
   ValNodePtr   seq_list;
   Uint1   add_redundant;
   Uint1   add_mrna;
   Uint1   apply_to_parts;
   Int4   only_seg_num;
   ValNodePtr   fields;
} ApplyFeatureAction, PNTR ApplyFeatureActionPtr;


NLM_EXTERN ApplyFeatureActionPtr LIBCALL ApplyFeatureActionFree PROTO ((ApplyFeatureActionPtr ));
NLM_EXTERN ApplyFeatureActionPtr LIBCALL ApplyFeatureActionNew PROTO (( void ));
NLM_EXTERN ApplyFeatureActionPtr LIBCALL ApplyFeatureActionAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL ApplyFeatureActionAsnWrite PROTO (( ApplyFeatureActionPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    RemoveFeatureAction
*
**************************************************/
typedef struct struct_Remove_feature_action {
   Uint2   type;
   ValNodePtr   constraint;
} RemoveFeatureAction, PNTR RemoveFeatureActionPtr;


NLM_EXTERN RemoveFeatureActionPtr LIBCALL RemoveFeatureActionFree PROTO ((RemoveFeatureActionPtr ));
NLM_EXTERN RemoveFeatureActionPtr LIBCALL RemoveFeatureActionNew PROTO (( void ));
NLM_EXTERN RemoveFeatureActionPtr LIBCALL RemoveFeatureActionAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL RemoveFeatureActionAsnWrite PROTO (( RemoveFeatureActionPtr , AsnIoPtr, AsnTypePtr));

/* following #defines are for enumerated type, not used by object loaders */
#define Feature_location_strand_from_any 0
#define Feature_location_strand_from_plus 1
#define Feature_location_strand_from_minus 2
#define Feature_location_strand_from_unknown 3
#define Feature_location_strand_from_both 4

/* following #defines are for enumerated type, not used by object loaders */
#define Feature_location_strand_to_plus 1
#define Feature_location_strand_to_minus 2
#define Feature_location_strand_to_unknown 3
#define Feature_location_strand_to_both 4
#define Feature_location_strand_to_reverse 5



/**************************************************
*
*    EditLocationStrand
*
**************************************************/
typedef struct struct_Edit_location_strand {
   Uint2   strand_from;
   Uint2   strand_to;
} EditLocationStrand, PNTR EditLocationStrandPtr;


NLM_EXTERN EditLocationStrandPtr LIBCALL EditLocationStrandFree PROTO ((EditLocationStrandPtr ));
NLM_EXTERN EditLocationStrandPtr LIBCALL EditLocationStrandNew PROTO (( void ));
NLM_EXTERN EditLocationStrandPtr LIBCALL EditLocationStrandAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL EditLocationStrandAsnWrite PROTO (( EditLocationStrandPtr , AsnIoPtr, AsnTypePtr));

/* following #defines are for enumerated type, not used by object loaders */
#define Partial_5_set_constraint_all 0
#define Partial_5_set_constraint_at_end 1
#define Partial_5_set_constraint_bad_start 2
#define Partial_5_set_constraint_frame_not_one 3



/**************************************************
*
*    Partial5SetAction
*
**************************************************/
typedef struct struct_Partial_5_set_action {
   Uint2   constraint;
   Uint1   extend;
} Partial5SetAction, PNTR Partial5SetActionPtr;


NLM_EXTERN Partial5SetActionPtr LIBCALL Partial5SetActionFree PROTO ((Partial5SetActionPtr ));
NLM_EXTERN Partial5SetActionPtr LIBCALL Partial5SetActionNew PROTO (( void ));
NLM_EXTERN Partial5SetActionPtr LIBCALL Partial5SetActionAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL Partial5SetActionAsnWrite PROTO (( Partial5SetActionPtr , AsnIoPtr, AsnTypePtr));

/* following #defines are for enumerated type, not used by object loaders */
#define Partial_5_clear_constraint_all 0
#define Partial_5_clear_constraint_not_at_end 1
#define Partial_5_clear_constraint_good_start 2

/* following #defines are for enumerated type, not used by object loaders */
#define Partial_3_set_constraint_all 0
#define Partial_3_set_constraint_at_end 1
#define Partial_3_set_constraint_bad_end 2



/**************************************************
*
*    Partial3SetAction
*
**************************************************/
typedef struct struct_Partial_3_set_action {
   Uint2   constraint;
   Uint1   extend;
} Partial3SetAction, PNTR Partial3SetActionPtr;


NLM_EXTERN Partial3SetActionPtr LIBCALL Partial3SetActionFree PROTO ((Partial3SetActionPtr ));
NLM_EXTERN Partial3SetActionPtr LIBCALL Partial3SetActionNew PROTO (( void ));
NLM_EXTERN Partial3SetActionPtr LIBCALL Partial3SetActionAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL Partial3SetActionAsnWrite PROTO (( Partial3SetActionPtr , AsnIoPtr, AsnTypePtr));

/* following #defines are for enumerated type, not used by object loaders */
#define Partial_3_clear_constraint_all 0
#define Partial_3_clear_constraint_not_at_end 1
#define Partial_3_clear_constraint_good_end 2

/* following #defines are for enumerated type, not used by object loaders */
#define Convert_location_type_join 1
#define Convert_location_type_order 2
#define Convert_location_type_merge 3

typedef ValNodePtr LocationEditTypePtr;
typedef ValNode LocationEditType;
#define LocationEditType_strand 1
#define LocationEditType_set_5_partial 2
#define LocationEditType_clear_5_partial 3
#define LocationEditType_set_3_partial 4
#define LocationEditType_clear_3_partial 5
#define LocationEditType_convert 6


NLM_EXTERN LocationEditTypePtr LIBCALL LocationEditTypeFree PROTO ((LocationEditTypePtr ));
NLM_EXTERN LocationEditTypePtr LIBCALL LocationEditTypeAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL LocationEditTypeAsnWrite PROTO (( LocationEditTypePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    EditFeatureLocationAction
*
**************************************************/
typedef struct struct_Edit_feature_location_action {
   Uint2   type;
   ValNodePtr   action;
   ValNodePtr   constraint;
} EditFeatureLocationAction, PNTR EditFeatureLocationActionPtr;


NLM_EXTERN EditFeatureLocationActionPtr LIBCALL EditFeatureLocationActionFree PROTO ((EditFeatureLocationActionPtr ));
NLM_EXTERN EditFeatureLocationActionPtr LIBCALL EditFeatureLocationActionNew PROTO (( void ));
NLM_EXTERN EditFeatureLocationActionPtr LIBCALL EditFeatureLocationActionAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL EditFeatureLocationActionAsnWrite PROTO (( EditFeatureLocationActionPtr , AsnIoPtr, AsnTypePtr));

typedef ValNodePtr MacroActionChoicePtr;
typedef ValNode MacroActionChoice;
#define MacroActionChoice_aecr 1
#define MacroActionChoice_parse 2
#define MacroActionChoice_add_feature 3
#define MacroActionChoice_remove_feature 4
#define MacroActionChoice_edit_location 5


NLM_EXTERN MacroActionChoicePtr LIBCALL MacroActionChoiceFree PROTO ((MacroActionChoicePtr ));
NLM_EXTERN MacroActionChoicePtr LIBCALL MacroActionChoiceAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL MacroActionChoiceAsnWrite PROTO (( MacroActionChoicePtr , AsnIoPtr, AsnTypePtr));

#ifdef __cplusplus
/* { */ }
#endif

#endif /* _objmacro_ */

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

