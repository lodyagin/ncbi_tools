/*   macro_i.h
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*            National Center for Biotechnology Information (NCBI)
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government do not place any restriction on its use or reproduction.
*  We would, however, appreciate having the NCBI and the author cited in
*  any work or product based on this material
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
* ===========================================================================
*
* File Name:  macro_i.h
*
* Author:  Colleen Bollin
*
* Version Creation Date:   11/15/2007
*
* $Revision: 1.31 $
*
* File Description: 
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
*
* ==========================================================================
*/

#ifndef _macroapi_h_
#define _macroapi_h_

#ifdef __cplusplus
extern "C" {
#endif

NLM_EXTERN Int4 GetFeatdefFromFeatureType (Int4 feature_type);
NLM_EXTERN CharPtr GetFeatureNameFromFeatureType (Int4 feature_type);
NLM_EXTERN Int4 GetFeatureTypeByName (CharPtr feat_name);
NLM_EXTERN void AddImportFeaturesToChoiceList (ValNodePtr PNTR feature_type_list);
NLM_EXTERN void AddAllFeaturesToChoiceList (ValNodePtr PNTR feature_type_list);
NLM_EXTERN CharPtr GetFeatQualName (Int4 featqual); 
NLM_EXTERN Int4 GetFeatQualByName (CharPtr qualname); 
NLM_EXTERN Int4 GetNumFeatQual (void);
NLM_EXTERN void AddAllFeatureFieldsToChoiceList (ValNodePtr PNTR field_list);
NLM_EXTERN CharPtr GetSourceQualName (Int4 srcqual);
NLM_EXTERN Int4 GetSourceQualTypeByName (CharPtr qualname);
NLM_EXTERN ValNodePtr GetSourceQualList (void);
NLM_EXTERN Boolean IsNonTextSourceQual (Int4 srcqual);
NLM_EXTERN Int4 GenomeFromSrcLoc (Int4 srcloc);
NLM_EXTERN CharPtr LocNameFromGenome (Int4 genome);
NLM_EXTERN ValNodePtr GetLocationList (Boolean for_remove); 
NLM_EXTERN Int4 OriginFromSrcOrig (Int4 srcorig);
NLM_EXTERN CharPtr OriginNameFromOrigin (Int4 origin);
NLM_EXTERN ValNodePtr GetOriginList (Boolean for_remove);
NLM_EXTERN BioSourcePtr GetBioSourceFromObject (Uint1 choice, Pointer data);
NLM_EXTERN CharPtr CDSGeneProtNameFromField (Int4 field); 
NLM_EXTERN CharPtr CDSGeneProtFeatureNameFromFeatureType (Int4 feature_type);
NLM_EXTERN void AddAllCDSGeneProtFieldsToChoiceList (ValNodePtr PNTR field_list);
NLM_EXTERN void AddAllCDSGeneProtFeaturesToChoiceList (ValNodePtr PNTR field_list);
NLM_EXTERN FeatureFieldPtr FeatureFieldFromCDSGeneProtField (Uint2 cds_gene_prot_field);

NLM_EXTERN CharPtr BiomolNameFromBiomol (Int4 biomol);
NLM_EXTERN Int4 BiomolFromMoleculeType (Int4 molecule_type);
NLM_EXTERN ValNodePtr GetMoleculeTypeList (void);
NLM_EXTERN CharPtr TechNameFromTech (Int4 tech);
NLM_EXTERN Int4 TechFromTechniqueType (Int4 technique_type);
NLM_EXTERN ValNodePtr GetTechniqueTypeList (void);
NLM_EXTERN Int4 CompletenessFromCompletednessType (Int4 completedness_type);
NLM_EXTERN CharPtr CompletenessNameFromCompleteness (Int4 completeness); 
NLM_EXTERN ValNodePtr GetCompletednessTypeList (void);
NLM_EXTERN Int4 MolFromMoleculeClassType (Int4 moleculeclass_type);
NLM_EXTERN CharPtr MolNameFromMol (Int4 mol); 
NLM_EXTERN ValNodePtr GetMoleculeClassTypeList (void);
NLM_EXTERN Int4 TopologyFromTopologyType (Int4 topology_type);
NLM_EXTERN CharPtr TopologyNameFromTopology (Int4 topology);
NLM_EXTERN ValNodePtr GetTopologyTypeList (void);
NLM_EXTERN Int4 StrandFromStrandType (Int4 strand_type);
NLM_EXTERN CharPtr StrandNameFromStrand (Int4 strand);
NLM_EXTERN ValNodePtr GetStrandTypeList (void);



NLM_EXTERN FieldTypePtr GetFromFieldFromFieldPair (FieldPairTypePtr fieldpair);
NLM_EXTERN FieldTypePtr GetToFieldFromFieldPair (FieldPairTypePtr fieldpair);
NLM_EXTERN Uint1 FieldTypeFromAECRAction (AECRActionPtr action);
NLM_EXTERN Uint1 GetBiomolForRnaType (Int4 rnatype);
NLM_EXTERN CharPtr GetBiomolNameForRnaType (Int4 rnatype);
NLM_EXTERN void AddAllRNASubtypesToChoiceList (ValNodePtr PNTR field_list);

NLM_EXTERN CharPtr GetSourceQualFromBioSource (BioSourcePtr biop, SourceQualChoicePtr scp, StringConstraintPtr constraint);
NLM_EXTERN CharPtr GetQualFromFeature (SeqFeatPtr sfp, FeatureFieldPtr field, StringConstraintPtr scp);
NLM_EXTERN Boolean SetSourceQualInBioSource (BioSourcePtr biop, SourceQualChoicePtr scp, StringConstraintPtr constraint, CharPtr value, Uint2 existing_text);


NLM_EXTERN Boolean IsStringConstraintEmpty (StringConstraintPtr scp);
NLM_EXTERN Boolean DoesSingleStringMatchConstraint (CharPtr str, StringConstraintPtr scp);
NLM_EXTERN Boolean DoesStringMatchConstraint (CharPtr str, StringConstraintPtr scp);
NLM_EXTERN Boolean IsSourceConstraintEmpty (SourceConstraintPtr scp);
NLM_EXTERN Boolean DoesBiosourceMatchConstraint (BioSourcePtr biop, SourceConstraintPtr scp);
NLM_EXTERN Boolean IsSequenceConstraintEmpty (SequenceConstraintPtr constraint);
NLM_EXTERN Boolean IsCDSGeneProtQualConstraintEmpty (CDSGeneProtQualConstraintPtr constraint);
NLM_EXTERN ValNodePtr GetObjectListForAECRAction (SeqEntryPtr sep, AECRActionPtr action);
NLM_EXTERN Int4 DoApplyActionToObjectList (ApplyActionPtr action, ValNodePtr object_list, StringConstraintPtr scp);
NLM_EXTERN Int4 DoEditActionToObjectList (EditActionPtr action, ValNodePtr object_list);
NLM_EXTERN Int4 DoConvertActionToObjectList (ConvertActionPtr action, ValNodePtr object_list, StringConstraintPtr scp);
NLM_EXTERN Int4 DoCopyActionToObjectList (CopyActionPtr action, ValNodePtr object_list, StringConstraintPtr scp);
NLM_EXTERN Int4 DoSwapActionToObjectList (SwapActionPtr action, ValNodePtr object_list, StringConstraintPtr scp);
NLM_EXTERN Int4 DoRemoveActionToObjectList (RemoveActionPtr action, ValNodePtr object_list, StringConstraintPtr scp);
NLM_EXTERN Int4 DoParseActionToObjectList (AECRParseActionPtr action, ValNodePtr object_list, StringConstraintPtr scp);
NLM_EXTERN StringConstraintPtr FindStringConstraintInConstraintSetForField (FieldTypePtr field, ConstraintChoiceSetPtr csp);
NLM_EXTERN StringConstraintPtr StringConstraintFromFieldEdit (FieldEditPtr edit);




NLM_EXTERN void ApplyMacroToSeqEntry (SeqEntryPtr sep, ValNodePtr macro, Int4Ptr pNumFields, Int4Ptr pNumFeat);

/* for generating text representations of macro objects */
NLM_EXTERN CharPtr SummarizeSourceQual (ValNodePtr field);
NLM_EXTERN CharPtr FeatureFieldLabel (CharPtr feature_name, ValNodePtr field);
NLM_EXTERN Boolean IsFieldTypeEmpty (FieldTypePtr field);
NLM_EXTERN CharPtr SummarizeFieldType (ValNodePtr vnp);

typedef enum {
  eTableMatchFeatureID = 1,
  eTableMatchGeneLocusTag,
  eTableMatchProteinID,
  eTableMatchDbxref,
  eTableMatchNucID,
  eTableMatchBioSource
} ETableMatchType;



typedef struct tabcolumnconfig {
  Uint2 match_type;
  FieldTypePtr field;
  Uint2 existing_text;
  Boolean skip_blank;
  Boolean match_mrna;
} TabColumnConfigData, PNTR TabColumnConfigPtr;

NLM_EXTERN TabColumnConfigPtr TabColumnConfigNew (void);
NLM_EXTERN TabColumnConfigPtr TabColumnConfigFree (TabColumnConfigPtr t);
NLM_EXTERN TabColumnConfigPtr TabColumnConfigCopy (TabColumnConfigPtr orig);
NLM_EXTERN ValNodePtr TabColumnConfigListFree (ValNodePtr columns);
NLM_EXTERN ValNodePtr TabColumnConfigListCopy (ValNodePtr orig);
NLM_EXTERN ValNodePtr ValidateTabTableValues (ValNodePtr table, ValNodePtr columns);
NLM_EXTERN ValNodePtr ValidateFeatureFieldColumnNames (ValNodePtr header_line, ValNodePtr PNTR perr_list);
NLM_EXTERN ValNodePtr FreeObjectTableForTabTable (ValNodePtr table);
NLM_EXTERN ValNodePtr GetObjectTableForTabTable (SeqEntryPtr sep, ValNodePtr table, ValNodePtr columns, ValNodePtr PNTR p_err_list);
NLM_EXTERN ValNodePtr ApplyTableValuesToObjectTable (SeqEntryPtr sep, ValNodePtr table, ValNodePtr columns, ValNodePtr obj_table);
NLM_EXTERN ValNodePtr CheckObjTableForRowsThatApplyToTheSameDestination (ValNodePtr obj_table);
NLM_EXTERN ValNodePtr CheckObjTableForExistingText (SeqEntryPtr sep, ValNodePtr table, ValNodePtr columns, ValNodePtr obj_table);

NLM_EXTERN ValNodePtr ApplyTableToFeatures (SeqEntryPtr sep, ValNodePtr table, ValNodePtr columns);
NLM_EXTERN ValNodePtr CheckTableForExistingText (SeqEntryPtr sep, ValNodePtr table, ValNodePtr columns);

NLM_EXTERN SeqFeatPtr GetmRNAForFeature (SeqFeatPtr sfp);
NLM_EXTERN Boolean AdjustmRNAProductToMatchProteinProduct (SeqFeatPtr sfp);
NLM_EXTERN Boolean IsFieldTypeCDSProduct (FieldTypePtr ft);


#ifdef __cplusplus 
} 
#endif

#endif
