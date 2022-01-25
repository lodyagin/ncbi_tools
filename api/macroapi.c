/*   macro.c
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
* File Name:  macro.c
*
* Author:  Colleen Bollin
*
* Version Creation Date:   11/8/2007
*
* $Revision: 1.58 $
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

#include <asn.h>
#include <objmacro.h>
#include <objfeat.h>
#include <subutil.h>
#include <objmgr.h>
#include <objfdef.h>
#include <gbftdef.h>
#include <sqnutils.h>
#include <edutil.h>
#include <gather.h>
#include <asn2gnbi.h>
#define NLM_GENERATED_CODE_PROTO
#include <macroapi.h>
#include <seqport.h>

/* structure and create/free functions for CGPSet, used for handling CDS-Gene-Prot sets */
typedef struct cgpset 
{
  ValNodePtr cds_list;
  ValNodePtr gene_list;
  ValNodePtr prot_list;
  ValNodePtr mrna_list;
} CGPSetData, PNTR CGPSetPtr;



static CGPSetPtr CGPSetNew (void)
{
  CGPSetPtr c;

  c = (CGPSetPtr) MemNew (sizeof (CGPSetData));
  c->cds_list = NULL;
  c->gene_list = NULL;
  c->prot_list = NULL;
  c->mrna_list = NULL;
  return c;
}


static CGPSetPtr CGPSetFree (CGPSetPtr c)
{
  if (c != NULL) {
    c->cds_list = ValNodeFree (c->cds_list);
    c->gene_list = ValNodeFree (c->gene_list);
    c->prot_list = ValNodeFree (c->prot_list);
    c->mrna_list = ValNodeFree (c->mrna_list);
    c = MemFree (c);
  }
  return c;
}


static ValNodePtr FreeCGPSetList (ValNodePtr vnp)
{
  ValNodePtr vnp_next;
  
  while (vnp != NULL) {
    vnp_next = vnp->next;
    vnp->next = NULL;
    vnp->data.ptrvalue = CGPSetFree (vnp->data.ptrvalue);
    vnp = ValNodeFree (vnp);
    vnp = vnp_next;
  }
  return NULL;
}


/* generic functions for mapping constraints */

typedef struct feattypefeatdef {
  Int4 feattype;
  Int4 featdef;
  CharPtr featname;
} FeatTypeFeatDefData, PNTR FeatTypeFeatDefPtr;

static FeatTypeFeatDefData feattype_featdef[] = {
 { Feature_type_any , FEATDEF_ANY , "any" } , 
 { Feature_type_gene , FEATDEF_GENE , "gene" } , 
 { Feature_type_org , FEATDEF_ORG , "org" } , 
 { Feature_type_cds , FEATDEF_CDS , "CDS" } , 
 { Feature_type_prot , FEATDEF_PROT , "Protein" } , 
 { Feature_type_preRNA , FEATDEF_preRNA , "preRNA" } , 
 { Feature_type_mRNA , FEATDEF_mRNA , "mRNA" } , 
 { Feature_type_tRNA , FEATDEF_tRNA , "tRNA" } , 
 { Feature_type_rRNA , FEATDEF_rRNA , "rRNA" } , 
 { Feature_type_snRNA , FEATDEF_snRNA , "snRNA" } , 
 { Feature_type_scRNA , FEATDEF_scRNA , "scRNA" } , 
 { Feature_type_otherRNA , FEATDEF_otherRNA , "misc_RNA" } , 
 { Feature_type_pub , FEATDEF_PUB , "pub" } , 
 { Feature_type_seq , FEATDEF_SEQ , "seq" } , 
 { Feature_type_imp , FEATDEF_IMP , "imp" } , 
 { Feature_type_allele , FEATDEF_allele , "allele" } , 
 { Feature_type_attenuator , FEATDEF_attenuator , "attenuator" } , 
 { Feature_type_c_region , FEATDEF_C_region , "c_region" } , 
 { Feature_type_caat_signal , FEATDEF_CAAT_signal , "caat_signal" } , 
 { Feature_type_imp_CDS , FEATDEF_Imp_CDS , "imp_CDS" } , 
 { Feature_type_conflict , FEATDEF_conflict , "conflict" } , 
 { Feature_type_d_loop , FEATDEF_D_loop , "d_loop" } , 
 { Feature_type_d_segment , FEATDEF_D_segment , "d_segment" } , 
 { Feature_type_enhancer , FEATDEF_enhancer , "enhancer" } , 
 { Feature_type_exon , FEATDEF_exon , "exon" } , 
 { Feature_type_gC_signal , FEATDEF_GC_signal , "gC_signal" } , 
 { Feature_type_iDNA , FEATDEF_iDNA , "iDNA" } , 
 { Feature_type_intron , FEATDEF_intron , "intron" } , 
 { Feature_type_j_segment , FEATDEF_J_segment , "j_segment" } , 
 { Feature_type_ltr , FEATDEF_LTR , "ltr" } , 
 { Feature_type_mat_peptide , FEATDEF_mat_peptide , "mat_peptide" } , 
 { Feature_type_misc_binding , FEATDEF_misc_binding , "misc_binding" } , 
 { Feature_type_misc_difference , FEATDEF_misc_difference , "misc_difference" } , 
 { Feature_type_misc_feature , FEATDEF_misc_feature , "misc_feature" } , 
 { Feature_type_misc_recomb , FEATDEF_misc_recomb , "misc_recomb" } , 
 { Feature_type_misc_RNA , FEATDEF_misc_RNA , "misc_RNA" } , 
 { Feature_type_misc_signal , FEATDEF_misc_signal , "misc_signal" } , 
 { Feature_type_misc_structure , FEATDEF_misc_structure , "misc_structure" } , 
 { Feature_type_modified_base , FEATDEF_modified_base , "modified_base" } , 
 { Feature_type_mutation , FEATDEF_mutation , "mutation" } , 
 { Feature_type_n_region , FEATDEF_N_region , "n_region" } , 
 { Feature_type_old_sequence , FEATDEF_old_sequence , "old_sequence" } , 
 { Feature_type_polyA_signal , FEATDEF_polyA_signal , "polyA_signal" } , 
 { Feature_type_polyA_site , FEATDEF_polyA_site , "polyA_site" } , 
 { Feature_type_precursor_RNA , FEATDEF_precursor_RNA , "precursor_RNA" } , 
 { Feature_type_prim_transcript , FEATDEF_prim_transcript , "prim_transcript" } , 
 { Feature_type_primer_bind , FEATDEF_primer_bind , "primer_bind" } , 
 { Feature_type_promoter , FEATDEF_promoter , "promoter" } , 
 { Feature_type_protein_bind , FEATDEF_protein_bind , "protein_bind" } , 
 { Feature_type_rbs , FEATDEF_RBS , "rbs" } , 
 { Feature_type_repeat_region , FEATDEF_repeat_region , "repeat_region" } , 
 { Feature_type_repeat_unit , FEATDEF_repeat_unit , "repeat_unit" } , 
 { Feature_type_rep_origin , FEATDEF_rep_origin , "rep_origin" } , 
 { Feature_type_s_region , FEATDEF_S_region , "s_region" } , 
 { Feature_type_satellite , FEATDEF_satellite , "satellite" } , 
 { Feature_type_sig_peptide , FEATDEF_sig_peptide , "sig_peptide" } , 
 { Feature_type_source , FEATDEF_source , "source" } , 
 { Feature_type_stem_loop , FEATDEF_stem_loop , "stem_loop" } , 
 { Feature_type_sts , FEATDEF_STS , "sts" } , 
 { Feature_type_tata_signal , FEATDEF_TATA_signal , "tata_signal" } , 
 { Feature_type_terminator , FEATDEF_terminator , "terminator" } , 
 { Feature_type_transit_peptide , FEATDEF_transit_peptide , "transit_peptide" } , 
 { Feature_type_unsure , FEATDEF_unsure , "unsure" } , 
 { Feature_type_v_region , FEATDEF_V_region , "v_region" } , 
 { Feature_type_v_segment , FEATDEF_V_segment , "v_segment" } , 
 { Feature_type_variation , FEATDEF_variation , "variation" } , 
 { Feature_type_virion , FEATDEF_virion , "virion" } , 
 { Feature_type_n3clip , FEATDEF_3clip , "3clip" } , 
 { Feature_type_n3UTR , FEATDEF_3UTR , "3UTR" } , 
 { Feature_type_n5clip , FEATDEF_5clip , "5clip" } , 
 { Feature_type_n5UTR , FEATDEF_5UTR , "5UTR" } , 
 { Feature_type_n10_signal , FEATDEF_10_signal , "10_signal" } , 
 { Feature_type_n35_signal , FEATDEF_35_signal , "35_signal" } , 
 { Feature_type_site_ref , FEATDEF_site_ref , "site_ref" } , 
 { Feature_type_region , FEATDEF_REGION , "region" } , 
 { Feature_type_comment , FEATDEF_COMMENT , "comment" } , 
 { Feature_type_bond , FEATDEF_BOND , "bond" } , 
 { Feature_type_site , FEATDEF_SITE , "site" } , 
 { Feature_type_rsite , FEATDEF_RSITE , "rsite" } , 
 { Feature_type_user , FEATDEF_USER , "user" } , 
 { Feature_type_txinit , FEATDEF_TXINIT , "txinit" } , 
 { Feature_type_num , FEATDEF_NUM , "num" } , 
 { Feature_type_psec_str , FEATDEF_PSEC_STR , "psec_str" } , 
 { Feature_type_non_std_residue , FEATDEF_NON_STD_RESIDUE , "non_std_residue" } , 
 { Feature_type_het , FEATDEF_HET , "het" } , 
 { Feature_type_biosrc , FEATDEF_BIOSRC , "biosrc" } , 
 { Feature_type_preprotein , FEATDEF_preprotein , "preprotein" } , 
 { Feature_type_mat_peptide_aa , FEATDEF_mat_peptide_aa , "mat_peptide_aa" } , 
 { Feature_type_sig_peptide_aa , FEATDEF_sig_peptide_aa , "sig_peptide_aa" } , 
 { Feature_type_transit_peptide_aa , FEATDEF_transit_peptide_aa , "transit_peptide_aa" } , 
 { Feature_type_snoRNA , FEATDEF_snoRNA , "snoRNA" } , 
 { Feature_type_gap , FEATDEF_gap , "gap" } , 
 { Feature_type_operon , FEATDEF_operon , "operon" } , 
 { Feature_type_oriT , FEATDEF_oriT , "oriT" } , 
 { Feature_type_ncRNA , FEATDEF_ncRNA , "ncRNA" } , 
 { Feature_type_tmRNA , FEATDEF_tmRNA , "tmRNA" }};

#define NUM_feattype_featdef sizeof (feattype_featdef) / sizeof (FeatTypeFeatDefData)

NLM_EXTERN Int4 GetFeatdefFromFeatureType (Int4 feature_type) 
{
  Int4 i;

  for (i = 0; i < NUM_feattype_featdef; i++) {
    if (feature_type == feattype_featdef[i].feattype) {
      return feattype_featdef[i].featdef;
    }
  }
  return FEATDEF_BAD;
}


NLM_EXTERN CharPtr GetFeatureNameFromFeatureType (Int4 feature_type)
{
  CharPtr str = NULL;
  Int4 i;

  for (i = 0; i < NUM_feattype_featdef && str == NULL; i++) {
    if (feature_type == feattype_featdef[i].feattype) {
      str = feattype_featdef[feature_type].featname;
    }
  } 
  if (str == NULL) {
    str = "Unknown feature type";
  }
  return str;
}


static Boolean Matchnamestring (CharPtr name1, CharPtr name2)
{
  if (name1 == NULL && name2 == NULL) {
    return TRUE;
  } else if (name1 == NULL || name2 == NULL) {
    return FALSE;
  } else {
    while (*name1 != 0 && *name2 != 0) {
      while (*name1 == ' ' || *name1 == '-' || *name1 == '_') {
        name1++;
      }
      while (*name2 == ' ' || *name2 == '-' || *name2 == '_') {
        name2++;
      }
      if (*name1 != *name2) {
        return FALSE;
      }
      name1++;
      name2++;
    }
    if (*name1 == 0 && *name2 == 0) {
      return TRUE;
    } else {
      return FALSE;
    }
  }
}


NLM_EXTERN Int4 GetFeatureTypeByName (CharPtr feat_name)
{
  Int4 i;

  for (i = 0; i < NUM_feattype_featdef; i++) {
    if (Matchnamestring (feattype_featdef[i].featname, feat_name)) {
      return feattype_featdef[i].feattype;
    }
  }
  return -1;  
}


NLM_EXTERN void AddImportFeaturesToChoiceList (ValNodePtr PNTR feature_type_list)
{
  Int4 i, seqfeattype;
  CharPtr featname;
  ValNodePtr tmp_list = NULL;

  for (i = 1; i < NUM_feattype_featdef; i++) {
    if (feattype_featdef[i].feattype == Feature_type_gap) continue;
    seqfeattype = FindFeatFromFeatDefType (feattype_featdef[i].featdef);
    if (seqfeattype == SEQFEAT_IMP) {
      featname = GetFeatureNameFromFeatureType (feattype_featdef[i].feattype);
      if (featname != NULL) {
        ValNodeAddPointer (&tmp_list, feattype_featdef[i].feattype, StringSave (featname));
      }
    }
  }
  tmp_list = ValNodeSort (tmp_list, SortVnpByString);
  ValNodeLink (feature_type_list, tmp_list);
}



static Boolean IsMostUsedFeature (Uint1 val)
{
  if (val == Feature_type_gene
      || val == Feature_type_cds
      || val == Feature_type_prot
      || val == Feature_type_exon
      || val == Feature_type_intron
      || val == Feature_type_mRNA
      || val == Feature_type_rRNA
      || val == Feature_type_otherRNA) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static int LIBCALLBACK SortVnpByFeatureName (VoidPtr ptr1, VoidPtr ptr2)

{
  CharPtr     str1;
  CharPtr     str2;
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;
  Boolean     most_used1, most_used2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      most_used1 = IsMostUsedFeature (vnp1->choice);
      most_used2 = IsMostUsedFeature (vnp2->choice);
      if (most_used1 && !most_used2) {
        return -1;
      } else if (!most_used1 && most_used2) {
        return 1;
      } else {
        str1 = (CharPtr) vnp1->data.ptrvalue;
        str2 = (CharPtr) vnp2->data.ptrvalue;
        if (str1 != NULL && str2 != NULL) {
          return StringICmp (str1, str2);
        }
      }
    }
  }
  return 0;
}


NLM_EXTERN void AddAllFeaturesToChoiceList (ValNodePtr PNTR feature_type_list)
{
  Int4 i;
  CharPtr featname;
  ValNodePtr tmp_list = NULL;

  for (i = 1; i < NUM_feattype_featdef; i++) {
    if (feattype_featdef[i].feattype == Feature_type_gap) continue;
    featname = GetFeatureNameFromFeatureType (feattype_featdef[i].feattype);
    if (featname != NULL) {
      ValNodeAddPointer (&tmp_list, feattype_featdef[i].feattype, StringSave (featname));
    }
  }
  tmp_list = ValNodeSort (tmp_list, SortVnpByFeatureName);
  ValNodeLink (feature_type_list, tmp_list);
}


typedef struct featqualgbqual {
  Int4 featqual;
  Int4 gbqual;
  CharPtr qualname;
} FeatQualGBQualData, PNTR FeatQualGBQualPtr;

static FeatQualGBQualData featqual_gbqual[] = {
 { Feat_qual_legal_allele , GBQUAL_allele , "allele" } , 
 { Feat_qual_legal_anticodon , GBQUAL_anticodon , "anticodon" } , 
 { Feat_qual_legal_bound_moiety , GBQUAL_bound_moiety , "bound-moiety" } , 
 { Feat_qual_legal_chromosome , GBQUAL_chromosome , "chromosome" } , 
 { Feat_qual_legal_citation , GBQUAL_citation , "citation" } , 
 { Feat_qual_legal_codon , GBQUAL_codon , "codon" } , 
 { Feat_qual_legal_codon_start , GBQUAL_codon_start , "codon-start" } , 
 { Feat_qual_legal_compare , GBQUAL_compare , "compare" } , 
 { Feat_qual_legal_cons_splice , GBQUAL_cons_splice , "cons-splice" } , 
 { Feat_qual_legal_db_xref , GBQUAL_db_xref , "db-xref" } , 
 { Feat_qual_legal_direction , GBQUAL_direction , "direction" } , 
 { Feat_qual_legal_environmental_sample , GBQUAL_environmental_sample , "environmental-sample" } , 
 { Feat_qual_legal_evidence , GBQUAL_evidence , "evidence" } , 
 { Feat_qual_legal_exception , GBQUAL_exception , "exception" } , 
 { Feat_qual_legal_experiment , GBQUAL_experiment , "experiment" } , 
 { Feat_qual_legal_focus , GBQUAL_focus , "focus" } , 
 { Feat_qual_legal_frequency , GBQUAL_frequency , "frequency" } , 
 { Feat_qual_legal_function , GBQUAL_function , "function" } , 
 { Feat_qual_legal_gene , GBQUAL_gene , "locus" } , 
 { Feat_qual_legal_inference , GBQUAL_inference , "inference" } , 
 { Feat_qual_legal_label , GBQUAL_label , "label" } , 
 { Feat_qual_legal_locus_tag , GBQUAL_locus_tag , "locus-tag" } , 
 { Feat_qual_legal_map , GBQUAL_map , "map" } , 
 { Feat_qual_legal_mobile_element , GBQUAL_mobile_element , "mobile-element" } , 
 { Feat_qual_legal_mod_base , GBQUAL_mod_base , "mod-base" } , 
 { Feat_qual_legal_mol_type , GBQUAL_mol_type , "mol-type" } , 
 { Feat_qual_legal_ncRNA_class , GBQUAL_ncRNA_class , "ncRNA-class" } , 
 { Feat_qual_legal_note , GBQUAL_note , "note" } , 
 { Feat_qual_legal_number , GBQUAL_number , "number" } , 
 { Feat_qual_legal_old_locus_tag , GBQUAL_old_locus_tag , "old-locus-tag" } , 
 { Feat_qual_legal_operon , GBQUAL_operon , "operon" } , 
 { Feat_qual_legal_organism , GBQUAL_organism , "organism" } , 
 { Feat_qual_legal_organelle , GBQUAL_organelle , "organelle" } , 
 { Feat_qual_legal_partial , GBQUAL_partial , "partial" } , 
 { Feat_qual_legal_phenotype , GBQUAL_phenotype , "phenotype" } , 
 { Feat_qual_legal_plasmid , GBQUAL_plasmid , "plasmid" } , 
 { Feat_qual_legal_product , GBQUAL_product , "product" } , 
 { Feat_qual_legal_protein_id , GBQUAL_protein_id , "protein-id" } , 
 { Feat_qual_legal_pseudo , GBQUAL_pseudo , "pseudo" } , 
 { Feat_qual_legal_rearranged , GBQUAL_rearranged , "rearranged" } , 
 { Feat_qual_legal_replace , GBQUAL_replace , "replace" } , 
 { Feat_qual_legal_rpt_family , GBQUAL_rpt_family , "rpt-family" } , 
 { Feat_qual_legal_rpt_type , GBQUAL_rpt_type , "rpt-type" } , 
 { Feat_qual_legal_rpt_unit , GBQUAL_rpt_unit , "rpt-unit" } , 
 { Feat_qual_legal_rpt_unit_seq , GBQUAL_rpt_unit_seq , "rpt-unit-seq" } , 
 { Feat_qual_legal_rpt_unit_range , GBQUAL_rpt_unit_range , "rpt-unit-range" } , 
 { Feat_qual_legal_segment , GBQUAL_segment , "segment" } , 
 { Feat_qual_legal_sequenced_mol , GBQUAL_sequenced_mol , "sequenced-mol" } , 
 { Feat_qual_legal_standard_name , GBQUAL_standard_name , "standard-name" } , 
 { Feat_qual_legal_transcript_id , GBQUAL_transcript_id , "transcript-id" } , 
 { Feat_qual_legal_transgenic , GBQUAL_transgenic , "transgenic" } , 
 { Feat_qual_legal_translation , GBQUAL_translation , "translation" } , 
 { Feat_qual_legal_transl_except , GBQUAL_transl_except , "transl-except" } , 
 { Feat_qual_legal_transl_table , GBQUAL_transl_table , "transl-table" } , 
 { Feat_qual_legal_usedin , GBQUAL_usedin , "usedin" } };

#define NUM_featqual_gbqual sizeof (featqual_gbqual) / sizeof (FeatQualGBQualData)


NLM_EXTERN Int4 GetNumFeatQual (void)
{
  return NUM_featqual_gbqual;
}


static Int4 GetGBQualFromFeatQual (Int4 featqual) 
{
  Int4 i;

  for (i = 0; i < NUM_featqual_gbqual; i++) {
    if (featqual == featqual_gbqual[i].featqual) {
      return featqual_gbqual[i].gbqual;
    }
  }
  return -1;
}


NLM_EXTERN CharPtr GetFeatQualName (Int4 featqual) 
{
  Int4 i;

  for (i = 0; i < NUM_featqual_gbqual; i++) {
    if (featqual == featqual_gbqual[i].featqual) {
      return featqual_gbqual[i].qualname;
    }
  }
  return NULL;
}


NLM_EXTERN Int4 GetFeatQualByName (CharPtr qualname) 
{
  Int4 i;

  for (i = 0; i < NUM_featqual_gbqual; i++) {
    if (Matchnamestring (featqual_gbqual[i].qualname, qualname)) {
      return featqual_gbqual[i].featqual;
    }
  }
  return -1;  
}


NLM_EXTERN void AddAllFeatureFieldsToChoiceList (ValNodePtr PNTR field_list)
{
  Int4 i;

  for (i = 1; i < NUM_featqual_gbqual; i++) {
    ValNodeAddPointer (field_list, featqual_gbqual[i].featqual, StringSave (featqual_gbqual[i].qualname));
  }
}


#define IS_ORGMOD 1
#define IS_SUBSRC 2
#define IS_OTHER  3

typedef struct srcqualscqual {
  Int4 srcqual;
  Int4 subtype;
  Int4 typeflag;
  CharPtr qualname;
} SrcQualSCQualData, PNTR SrcQualSCQualPtr;

static SrcQualSCQualData srcqual_scqual[] = {
 { Source_qual_acronym , ORGMOD_acronym , IS_ORGMOD , "acronym" } , 
 { Source_qual_anamorph , ORGMOD_anamorph , IS_ORGMOD , "anamorph" } , 
 { Source_qual_authority , ORGMOD_authority , IS_ORGMOD , "authority" } , 
 { Source_qual_bio_material , ORGMOD_bio_material , IS_ORGMOD , "bio-material" } , 
 { Source_qual_biotype , ORGMOD_biotype , IS_ORGMOD , "biotype" } , 
 { Source_qual_biovar , ORGMOD_biovar , IS_ORGMOD , "biovar" } , 
 { Source_qual_breed , ORGMOD_breed , IS_ORGMOD , "breed" } , 
 { Source_qual_cell_line , SUBSRC_cell_line , IS_SUBSRC , "cell-line" } , 
 { Source_qual_cell_type , SUBSRC_cell_type , IS_SUBSRC , "cell-type" } , 
 { Source_qual_chemovar , ORGMOD_chemovar , IS_ORGMOD , "chemovar" } , 
 { Source_qual_chromosome , SUBSRC_chromosome , IS_SUBSRC , "chromosome" } , 
 { Source_qual_clone , SUBSRC_clone , IS_SUBSRC , "clone" } , 
 { Source_qual_clone_lib , SUBSRC_clone_lib , IS_SUBSRC , "clone-lib" } , 
 { Source_qual_collected_by , SUBSRC_collected_by , IS_SUBSRC , "collected-by" } , 
 { Source_qual_collection_date , SUBSRC_collection_date , IS_SUBSRC , "collection-date" } , 
 { Source_qual_common , ORGMOD_common , IS_ORGMOD , "common" } , 
 { Source_qual_common_name , 0 , IS_OTHER , "common name" } , 
 { Source_qual_country , SUBSRC_country , IS_SUBSRC , "country" } , 
 { Source_qual_cultivar , ORGMOD_cultivar , IS_ORGMOD , "cultivar" } , 
 { Source_qual_culture_collection , ORGMOD_culture_collection , IS_ORGMOD , "culture-collection" } , 
 { Source_qual_dev_stage , SUBSRC_dev_stage , IS_SUBSRC , "dev-stage" } , 
 { Source_qual_division , 0 , IS_OTHER, "divistion" } ,
 { Source_qual_dosage , ORGMOD_dosage , IS_ORGMOD , "dosage" } , 
 { Source_qual_ecotype , ORGMOD_ecotype , IS_ORGMOD , "ecotype" } , 
 { Source_qual_endogenous_virus_name , SUBSRC_endogenous_virus_name , IS_SUBSRC , "endogenous-virus-name" } , 
 { Source_qual_environmental_sample , SUBSRC_environmental_sample , IS_SUBSRC , "environmental-sample" } , 
 { Source_qual_forma , ORGMOD_forma , IS_ORGMOD , "forma" } , 
 { Source_qual_forma_specialis , ORGMOD_forma_specialis , IS_ORGMOD , "forma-specialis" } , 
 { Source_qual_frequency , SUBSRC_frequency , IS_SUBSRC , "frequency" } , 
 { Source_qual_fwd_primer_name , SUBSRC_fwd_primer_name , IS_SUBSRC , "fwd-primer-name" } , 
 { Source_qual_fwd_primer_seq , SUBSRC_fwd_primer_seq , IS_SUBSRC , "fwd-primer-seq" } , 
 { Source_qual_gb_acronym , ORGMOD_gb_acronym , IS_ORGMOD , "gb-acronym" } , 
 { Source_qual_gb_anamorph , ORGMOD_gb_anamorph , IS_ORGMOD , "gb-anamorph" } , 
 { Source_qual_gb_synonym , ORGMOD_gb_synonym , IS_ORGMOD , "gb-synonym" } , 
 { Source_qual_genotype , SUBSRC_genotype , IS_SUBSRC , "genotype" } , 
 { Source_qual_germline , SUBSRC_germline , IS_SUBSRC , "germline" } , 
 { Source_qual_group , ORGMOD_group , IS_ORGMOD , "group" } , 
 { Source_qual_haplotype , SUBSRC_haplotype , IS_SUBSRC , "haplotype" } , 
 { Source_qual_identified_by , SUBSRC_identified_by , IS_SUBSRC , "identified-by" } , 
 { Source_qual_insertion_seq_name , SUBSRC_insertion_seq_name , IS_SUBSRC , "insertion-seq-name" } , 
 { Source_qual_isolate , ORGMOD_isolate , IS_ORGMOD , "isolate" } , 
 { Source_qual_isolation_source , SUBSRC_isolation_source , IS_SUBSRC , "isolation-source" } , 
 { Source_qual_lab_host , SUBSRC_lab_host , IS_SUBSRC , "lab-host" } , 
 { Source_qual_lat_lon , SUBSRC_lat_lon , IS_SUBSRC , "lat-lon" } , 
 { Source_qual_lineage , 0, IS_OTHER, "lineage" } ,
 { Source_qual_map , SUBSRC_map , IS_SUBSRC , "map" } , 
 { Source_qual_metagenome_source , ORGMOD_metagenome_source , IS_ORGMOD , "metagenome-source" } , 
 { Source_qual_metagenomic , SUBSRC_metagenomic , IS_SUBSRC , "metagenomic" } , 
 { Source_qual_old_lineage , ORGMOD_old_lineage , IS_ORGMOD , "old-lineage" } , 
 { Source_qual_old_name , ORGMOD_old_name , IS_ORGMOD , "old-name" } , 
 { Source_qual_orgmod_note , ORGMOD_other, IS_ORGMOD, "orgmod note" } ,
 { Source_qual_nat_host , ORGMOD_nat_host , IS_ORGMOD , "nat-host" } , 
 { Source_qual_pathovar , ORGMOD_pathovar , IS_ORGMOD , "pathovar" } , 
 { Source_qual_plasmid_name , SUBSRC_plasmid_name , IS_SUBSRC , "plasmid-name" } , 
 { Source_qual_plastid_name , SUBSRC_plastid_name , IS_SUBSRC , "plastid-name" } , 
 { Source_qual_pop_variant , SUBSRC_pop_variant , IS_SUBSRC , "pop-variant" } , 
 { Source_qual_rearranged , SUBSRC_rearranged , IS_SUBSRC , "rearranged" } , 
 { Source_qual_rev_primer_name , SUBSRC_rev_primer_name , IS_SUBSRC , "rev-primer-name" } , 
 { Source_qual_rev_primer_seq , SUBSRC_rev_primer_seq , IS_SUBSRC , "rev-primer-seq" } , 
 { Source_qual_segment , SUBSRC_segment , IS_SUBSRC , "segment" } , 
 { Source_qual_serogroup , ORGMOD_serogroup , IS_ORGMOD , "serogroup" } , 
 { Source_qual_serotype , ORGMOD_serotype , IS_ORGMOD , "serotype" } , 
 { Source_qual_serovar , ORGMOD_serovar , IS_ORGMOD , "serovar" } , 
 { Source_qual_sex , SUBSRC_sex , IS_SUBSRC , "sex" } , 
 { Source_qual_specimen_voucher , ORGMOD_specimen_voucher , IS_ORGMOD , "specimen-voucher" } , 
 { Source_qual_strain , ORGMOD_strain , IS_ORGMOD , "strain" } , 
 { Source_qual_subclone , SUBSRC_subclone , IS_SUBSRC , "subclone" } , 
 { Source_qual_subgroup , ORGMOD_subgroup , IS_ORGMOD , "subgroup" } , 
 { Source_qual_subsource_note , SUBSRC_other , IS_SUBSRC , "subsource note" } ,
 { Source_qual_sub_species , ORGMOD_sub_species , IS_ORGMOD , "sub-species" } , 
 { Source_qual_substrain , ORGMOD_substrain , IS_ORGMOD , "substrain" } , 
 { Source_qual_subtype , ORGMOD_subtype , IS_ORGMOD , "subtype" } , 
 { Source_qual_synonym , ORGMOD_synonym , IS_ORGMOD , "synonym" } , 
 { Source_qual_taxname , 0 , IS_OTHER , "taxname" } , 
 { Source_qual_teleomorph , ORGMOD_teleomorph , IS_ORGMOD , "teleomorph" } , 
 { Source_qual_tissue_lib , SUBSRC_tissue_lib , IS_SUBSRC , "tissue-lib" } , 
 { Source_qual_tissue_type , SUBSRC_tissue_type , IS_SUBSRC , "tissue-type" } , 
 { Source_qual_transgenic , SUBSRC_transgenic , IS_SUBSRC , "transgenic" } , 
 { Source_qual_transposon_name , SUBSRC_transposon_name , IS_SUBSRC , "transposon-name" } , 
 { Source_qual_type , ORGMOD_type , IS_ORGMOD , "type" } , 
 { Source_qual_variety , ORGMOD_variety , IS_ORGMOD , "variety" } };

#define NUM_srcqual_scqual sizeof (srcqual_scqual) / sizeof (SrcQualSCQualData)

static Int4 GetSubSrcQualFromSrcQual (Int4 srcqual) 
{
  Int4 i;

  for (i = 0; i < NUM_srcqual_scqual; i++) {
    if (srcqual == srcqual_scqual[i].srcqual) {
      if (srcqual_scqual[i].typeflag == IS_SUBSRC) {
        return srcqual_scqual[i].subtype;
      } else {
        return -1;
      }
    }
  }
  return -1;
}


static Int4 GetOrgModQualFromSrcQual (Int4 srcqual) 
{
  Int4 i;

  for (i = 0; i < NUM_srcqual_scqual; i++) {
    if (srcqual == srcqual_scqual[i].srcqual) {
      if (srcqual_scqual[i].typeflag == IS_ORGMOD) {
        return srcqual_scqual[i].subtype;
      } else {
        return -1;
      }
    }
  }
  return -1;
}


NLM_EXTERN Boolean IsNonTextSourceQual (Int4 srcqual)
{
  if (srcqual == Source_qual_transgenic
      || srcqual == Source_qual_germline
      || srcqual == Source_qual_metagenomic
      || srcqual == Source_qual_environmental_sample
      || srcqual == Source_qual_rearranged)
  {
    return TRUE;  
  }
  else
  {
    return FALSE;
  }
}


NLM_EXTERN CharPtr GetSourceQualName (Int4 srcqual)
{
  CharPtr str = NULL;
  Int4    i;

  for (i = 0; i < NUM_srcqual_scqual && str == NULL; i++) {
    if (srcqual_scqual[i].srcqual == srcqual) {
      str = srcqual_scqual[i].qualname;
    }
  }
  if (str == NULL) {
    str = "Unknown source qualifier";
  }
  return str;
}


NLM_EXTERN Int4 GetSourceQualTypeByName (CharPtr qualname)
{
  Int4    i;

  for (i = 0; i < NUM_srcqual_scqual; i++) {
    if (Matchnamestring(srcqual_scqual[i].qualname, qualname)) {
      return srcqual_scqual[i].srcqual;
    }
  }
  return -1;
}


NLM_EXTERN ValNodePtr GetSourceQualList (void)
{
  ValNodePtr list = NULL;
  Int4 i;

  for (i = 0; i < NUM_srcqual_scqual; i++) {
    ValNodeAddPointer (&list, 0, StringSave (srcqual_scqual[i].qualname));
  }
  return list;
}

typedef struct srclocgenome {
  Int4 srcloc;
  Int4 genome;
  CharPtr name;
} SrcLocGenomeData, PNTR SrcLocGenomePtr;

static SrcLocGenomeData srcloc_genome[] = {
 { Source_location_unknown , GENOME_unknown , "unknown" } ,
 { Source_location_genomic , GENOME_genomic , "genomic" } ,
 { Source_location_chloroplast , GENOME_chloroplast , "chloroplast" } ,
 { Source_location_chromoplast , GENOME_chromoplast , "chromoplast" } ,
 { Source_location_kinetoplast , GENOME_kinetoplast , "kinetoplast" } ,
 { Source_location_mitochondrion , GENOME_mitochondrion , "mitochondrion" } ,
 { Source_location_plastid , GENOME_plastid , "plastid" } ,
 { Source_location_macronuclear , GENOME_macronuclear , "macronuclear" } ,
 { Source_location_extrachrom , GENOME_extrachrom , "extrachrom" } ,
 { Source_location_plasmid , GENOME_plasmid , "plasmid" } ,
 { Source_location_transposon , GENOME_transposon , "transposon" } ,
 { Source_location_insertion_seq , GENOME_insertion_seq , "insertion-seq" } ,
 { Source_location_cyanelle , GENOME_cyanelle , "cyanelle" } ,
 { Source_location_proviral , GENOME_proviral , "proviral" } ,
 { Source_location_virion , GENOME_virion , "virion" } ,
 { Source_location_nucleomorph , GENOME_nucleomorph , "nucleomorph" } ,
 { Source_location_apicoplast , GENOME_apicoplast , "apicoplast" } ,
 { Source_location_leucoplast , GENOME_leucoplast , "leucoplast" } ,
 { Source_location_proplastid , GENOME_proplastid , "proplastid" } ,
 { Source_location_endogenous_virus , GENOME_endogenous_virus , "endogenous-virus" } ,
 { Source_location_hydrogenosome , GENOME_hydrogenosome , "hydrogenosome" } ,
 { Source_location_chromosome , 21 , "chromosome" } ,
 { Source_location_chromatophore , 22 , "chromatophore" } };

#define NUM_srcloc_genome sizeof (srcloc_genome) / sizeof (SrcLocGenomeData)

NLM_EXTERN Int4 GenomeFromSrcLoc (Int4 srcloc) \
{
  Int4 i;

  for (i = 0; i < NUM_srcloc_genome; i++) {
    if (srcloc_genome[i].srcloc == srcloc) {
      return srcloc_genome[i].genome;
    }
  }
  return -1;
}


NLM_EXTERN CharPtr LocNameFromGenome (Int4 genome) 
{
  Int4 i;

  for (i = 0; i < NUM_srcloc_genome; i++) {
    if (srcloc_genome[i].genome == genome) {
      return srcloc_genome[i].name;
    }
  }
  return NULL;
}


static Int4 GenomeFromLocName (CharPtr loc_name)
{
  Int4 i;

  for (i = 0; i < NUM_srcloc_genome; i++) {
    if (StringICmp (srcloc_genome[i].name, loc_name) == 0) {
      return srcloc_genome[i].genome;
    }
  }
  return -1;
}


NLM_EXTERN ValNodePtr GetLocationList (Boolean for_remove)
{
  ValNodePtr list = NULL;
  Int4 i;

  for (i = 0; i < NUM_srcloc_genome; i++) {
    if (for_remove && srcloc_genome[i].srcloc == Source_location_unknown) {
      ValNodeAddPointer (&list, srcloc_genome[i].srcloc, StringSave ("any"));
    } else {
      ValNodeAddPointer (&list, srcloc_genome[i].srcloc, StringSave (srcloc_genome[i].name));
    }
  }
  return list;
}


typedef struct srcorigorigin {
  Int4 srcorig;
  Int4 origin;
  CharPtr name;
} SrcOrigOriginData, PNTR SrcrigOriginPtr;

static SrcOrigOriginData srcorig_origin[] = {
 { Source_origin_unknown , 0 , "unknown" } ,
 { Source_origin_natural , 1 , "natural" } ,
 { Source_origin_natmut , 2 , "natmut" } ,
 { Source_origin_mut , 3 , "mut" } ,
 { Source_origin_artificial , 4 , "artificial" } ,
 { Source_origin_synthetic , 5 , "synthetic" } ,
 { Source_origin_other , 255 , "other" } };

#define NUM_srcorig_origin sizeof (srcorig_origin) / sizeof (SrcOrigOriginData)

NLM_EXTERN Int4 OriginFromSrcOrig (Int4 srcorig) 
{
  Int4 i;

  for (i = 0; i < NUM_srcorig_origin; i++) {
    if (srcorig_origin[i].srcorig == srcorig) {
      return srcorig_origin[i].origin;
    }
  }
  return -1;
}


NLM_EXTERN CharPtr OriginNameFromOrigin (Int4 origin) 
{
  Int4 i;

  for (i = 0; i < NUM_srcorig_origin; i++) {
    if (srcorig_origin[i].origin == origin) {
      return srcorig_origin[i].name;
    }
  }
  return NULL;
}


static Int4 OriginFromOriginName (CharPtr origin_name)
{
  Int4 i;

  for (i = 0; i < NUM_srcorig_origin; i++) {
    if (StringCmp (srcorig_origin[i].name, origin_name) == 0) {
      return srcorig_origin[i].origin;
    }
  }
  return -1;
}


NLM_EXTERN ValNodePtr GetOriginList (Boolean for_remove)
{
  ValNodePtr list = NULL;
  Int4 i;

  for (i = 0; i < NUM_srcorig_origin; i++) {
    if (for_remove && srcorig_origin[i].srcorig == Source_origin_unknown) {
      ValNodeAddPointer (&list, srcorig_origin[i].srcorig, StringSave ("any"));
    } else {
      ValNodeAddPointer (&list, srcorig_origin[i].srcorig, StringSave (srcorig_origin[i].name));
    }
  }
  return list;
}


typedef struct cdsgeneprotfieldname {
  Int4 field;
  CharPtr name;
} CDSGeneProtFieldNameData, PNTR CDSGeneProtFieldNamePtr;

static CDSGeneProtFieldNameData cdsgeneprotfield_name[] = {
{ CDSGeneProt_field_cds_comment , "CDS comment" } ,
{ CDSGeneProt_field_gene_locus , "gene locus" } ,
{ CDSGeneProt_field_gene_description , "gene description" } ,
{ CDSGeneProt_field_gene_comment , "gene comment" } ,
{ CDSGeneProt_field_gene_allele , "allele" } ,
{ CDSGeneProt_field_gene_maploc , "maploc" } ,
{ CDSGeneProt_field_gene_locus_tag , "locus tag" } ,
{ CDSGeneProt_field_gene_synonym , "synonym" } ,
{ CDSGeneProt_field_gene_old_locus_tag , "old locus tag" } ,
{ CDSGeneProt_field_mrna_product , "mRNA product" } ,
{ CDSGeneProt_field_mrna_comment , "mRNA comment" } ,
{ CDSGeneProt_field_prot_name , "protein name" } ,
{ CDSGeneProt_field_prot_description , "protein description" } ,
{ CDSGeneProt_field_prot_ec_number , "protein EC number" } ,
{ CDSGeneProt_field_prot_activity , "protein activity" } ,
{ CDSGeneProt_field_prot_comment , "protein comment" } ,
{ CDSGeneProt_field_mat_peptide_name , "mat-peptide name" } ,
{ CDSGeneProt_field_mat_peptide_description ,  "mat-peptide description" } ,
{ CDSGeneProt_field_mat_peptide_ec_number , "mat-peptide EC number" } ,
{ CDSGeneProt_field_mat_peptide_activity , "mat-peptide activity" } ,
{ CDSGeneProt_field_mat_peptide_comment , "mat-peptide comment" } };

#define NUM_cdsgeneprotfield_name sizeof (cdsgeneprotfield_name) / sizeof (CDSGeneProtFieldNameData)

NLM_EXTERN CharPtr CDSGeneProtNameFromField (Int4 field) 
{
  Int4 i;

  for (i = 0; i < NUM_cdsgeneprotfield_name; i++) {
    if (cdsgeneprotfield_name[i].field == field) {
      return cdsgeneprotfield_name[i].name;
    }
  }
  return NULL;
}


NLM_EXTERN void AddAllCDSGeneProtFieldsToChoiceList (ValNodePtr PNTR field_list)
{
  Int4 i;

  for (i = 0; i < NUM_cdsgeneprotfield_name; i++) {
    ValNodeAddPointer (field_list, cdsgeneprotfield_name[i].field, StringSave (cdsgeneprotfield_name[i].name));
  }
}


typedef struct cdsgeneprotfeatname {
  Int4 feature_type;
  CharPtr name;
} CDSGeneProtFeatNameData, PNTR CDSGeneProtFeatNamePtr;

static CDSGeneProtFeatNameData cdsgeneprotfeat_name[] = {
{ CDSGeneProt_feature_type_constraint_gene , "gene" } ,
{ CDSGeneProt_feature_type_constraint_mRNA , "mRNA" } ,
{ CDSGeneProt_feature_type_constraint_cds , "CDS" } ,
{ CDSGeneProt_feature_type_constraint_prot , "protein" } ,
{ CDSGeneProt_feature_type_constraint_mat_peptide , "mat-peptide" }};

#define NUM_cdsgeneprotfeat_name sizeof (cdsgeneprotfeat_name) / sizeof (CDSGeneProtFeatNameData)

NLM_EXTERN CharPtr CDSGeneProtFeatureNameFromFeatureType (Int4 feature_type)
{
  Int4 i;

  for (i = 0; i < NUM_cdsgeneprotfeat_name; i++) {
    if (cdsgeneprotfeat_name[i].feature_type == feature_type) {
      return cdsgeneprotfeat_name[i].name;
    }
  }
  return NULL;
}


NLM_EXTERN void AddAllCDSGeneProtFeaturesToChoiceList (ValNodePtr PNTR field_list)
{
  Int4 i;

  for (i = 0; i < NUM_cdsgeneprotfeat_name; i++) {
    ValNodeAddPointer (field_list, cdsgeneprotfeat_name[i].feature_type, StringSave (cdsgeneprotfeat_name[i].name));
  }
}


NLM_EXTERN FeatureFieldPtr FeatureFieldFromCDSGeneProtField (Uint2 cds_gene_prot_field)
{
  FeatureFieldPtr f = NULL;

  switch (cds_gene_prot_field) {
    case CDSGeneProt_field_cds_comment:
      f = FeatureFieldNew ();
      f->type = Feature_type_cds;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_note;
      break;
    case CDSGeneProt_field_gene_locus:
      f = FeatureFieldNew ();
      f->type = Feature_type_gene;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_gene;
      break;
    case CDSGeneProt_field_gene_description:
      f = FeatureFieldNew ();
      f->type = Feature_type_gene;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_gene_description;
      break;
    case CDSGeneProt_field_gene_comment:
      f = FeatureFieldNew ();
      f->type = Feature_type_gene;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_note;
      break;
    case CDSGeneProt_field_gene_allele:
      f = FeatureFieldNew ();
      f->type = Feature_type_gene;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_allele;
      break;
    case CDSGeneProt_field_gene_maploc:
      f = FeatureFieldNew ();
      f->type = Feature_type_gene;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_map;
      break;
    case CDSGeneProt_field_gene_locus_tag:
      f = FeatureFieldNew ();
      f->type = Feature_type_gene;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_locus_tag;
      break;
    case CDSGeneProt_field_gene_synonym:
      f = FeatureFieldNew ();
      f->type = Feature_type_gene;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_synonym;
      break;
    case CDSGeneProt_field_gene_old_locus_tag:
      f = FeatureFieldNew ();
      f->type = Feature_type_gene;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_old_locus_tag;
      break;
    case CDSGeneProt_field_mrna_product:
      f = FeatureFieldNew ();
      f->type = Feature_type_mRNA;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_product;
      break;
    case CDSGeneProt_field_mrna_comment:
      f = FeatureFieldNew ();
      f->type = Feature_type_mRNA;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_note;
      break;
    case CDSGeneProt_field_prot_name:
      f = FeatureFieldNew ();
      f->type = Feature_type_prot;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_product;
      break;
    case CDSGeneProt_field_prot_description:
      f = FeatureFieldNew ();
      f->type = Feature_type_prot;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_description;
      break;
    case CDSGeneProt_field_prot_ec_number:
      f = FeatureFieldNew ();
      f->type = Feature_type_prot;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_ec_number;
      break;
    case CDSGeneProt_field_prot_activity:
      f = FeatureFieldNew ();
      f->type = Feature_type_prot;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_activity;
      break;
    case CDSGeneProt_field_prot_comment:
      f = FeatureFieldNew ();
      f->type = Feature_type_prot;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_note;
      break;
    case CDSGeneProt_field_mat_peptide_name:
      f = FeatureFieldNew ();
      f->type = Feature_type_mat_peptide;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_product;
      break;
    case CDSGeneProt_field_mat_peptide_description:
      f = FeatureFieldNew ();
      f->type = Feature_type_mat_peptide;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_description;
      break;
    case CDSGeneProt_field_mat_peptide_ec_number:
      f = FeatureFieldNew ();
      f->type = Feature_type_mat_peptide;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_ec_number;
      break;
    case CDSGeneProt_field_mat_peptide_activity:
      f = FeatureFieldNew ();
      f->type = Feature_type_mat_peptide;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_activity;
      break;
    case CDSGeneProt_field_mat_peptide_comment:
      f = FeatureFieldNew ();
      f->type = Feature_type_mat_peptide;
      f->field = ValNodeNew (NULL);
      f->field->choice = FeatQualChoice_legal_qual;
      f->field->data.intvalue = Feat_qual_legal_note;
      break;
  }
  return f;
}


/* Molinfo fields */
typedef struct moleculetypebiomol {
  Int4 molecule_type;
  Int4 biomol;
  CharPtr name;
} MoleculeTypeBiomolData, PNTR MoleculeTypeBiomolPtr;

static MoleculeTypeBiomolData moleculetype_biomol[] = {
 { Molecule_type_unknown , 0, "unknown" } ,
 { Molecule_type_genomic , MOLECULE_TYPE_GENOMIC , "genomic" } ,
 { Molecule_type_precursor_RNA , MOLECULE_TYPE_PRE_MRNA , "precursor RNA" } ,
 { Molecule_type_mRNA , MOLECULE_TYPE_MRNA , "mRNA" } ,
 { Molecule_type_rRNA , MOLECULE_TYPE_RRNA , "rRNA" } ,
 { Molecule_type_tRNA , MOLECULE_TYPE_TRNA , "tRNA" } ,
 { Molecule_type_snRNA , MOLECULE_TYPE_SNRNA , "snRNA" } ,
 { Molecule_type_scRNA , MOLECULE_TYPE_SCRNA , "scRNA" } ,
 { Molecule_type_genomic_mRNA , MOLECULE_TYPE_GENOMIC_MRNA_MIX , "genomic mRNA" } ,
 { Molecule_type_cRNA , MOLECULE_TYPE_CRNA , "cRNA" } ,
 { Molecule_type_scRNA , MOLECULE_TYPE_SCRNA , "scRNA" } ,
 { Molecule_type_scRNA , MOLECULE_TYPE_SCRNA , "scRNA" } ,
 { Molecule_type_scRNA , MOLECULE_TYPE_SCRNA , "scRNA" } ,
 { Molecule_type_snoRNA , MOLECULE_TYPE_SNORNA, "snoRNA" } ,
 { Molecule_type_transcribed_RNA, MOLECULE_TYPE_TRANSCRIBED_RNA, "transcribed RNA" } ,
 { Molecule_type_ncRNA, MOLECULE_TYPE_NCRNA, "ncRNA" } ,
 { Molecule_type_transfer_messenger_RNA, MOLECULE_TYPE_TMRNA, "tmRNA" } ,
 { Molecule_type_other, MOLECULE_TYPE_OTHER_GENETIC_MATERIAL, "other" }
};


#define NUM_moleculetype_biomol sizeof (moleculetype_biomol) / sizeof (MoleculeTypeBiomolData)

NLM_EXTERN Int4 BiomolFromMoleculeType (Int4 molecule_type) 
{
  Int4 i;

  for (i = 0; i < NUM_moleculetype_biomol; i++) {
    if (moleculetype_biomol[i].molecule_type == molecule_type) {
      return moleculetype_biomol[i].biomol;
    }
  }
  return -1;
}


NLM_EXTERN CharPtr BiomolNameFromBiomol (Int4 biomol) 
{
  Int4 i;

  for (i = 0; i < NUM_moleculetype_biomol; i++) {
    if (moleculetype_biomol[i].biomol == biomol) {
      return moleculetype_biomol[i].name;
    }
  }
  return NULL;
}


static Int4 BiomolFromBiomolName (CharPtr biomol_name)
{
  Int4 i;

  for (i = 0; i < NUM_moleculetype_biomol; i++) {
    if (StringCmp (moleculetype_biomol[i].name, biomol_name) == 0) {
      return moleculetype_biomol[i].biomol;
    }
  }
  return -1;
}


NLM_EXTERN ValNodePtr GetMoleculeTypeList (void)
{
  ValNodePtr list = NULL;
  Int4 i;

  for (i = 0; i < NUM_moleculetype_biomol; i++) {
    ValNodeAddPointer (&list, moleculetype_biomol[i].molecule_type, StringSave (moleculetype_biomol[i].name));
  }
  return list;
}


/* Technique fields */
typedef struct techniquetypetech {
  Int4 technique_type;
  Int4 tech;
  CharPtr name;
} TechniqueTypeTechData, PNTR TechniqueTypeTechPtr;

static TechniqueTypeTechData techniquetype_tech[] = {
 { Technique_type_unknown , MI_TECH_unknown , "unknown" } ,
 { Technique_type_standard , MI_TECH_standard , "standard" } ,
 { Technique_type_est , MI_TECH_est , "EST" } ,
 { Technique_type_sts , MI_TECH_sts , "STS" } ,
 { Technique_type_survey , MI_TECH_survey , "survey" } ,
 { Technique_type_genetic_map , MI_TECH_genemap , "genetic map" } ,
 { Technique_type_physical_map , MI_TECH_physmap , "physical map" } ,
 { Technique_type_derived , MI_TECH_derived , "derived" } ,
 { Technique_type_concept_trans , MI_TECH_concept_trans , "concept-trans" } ,
 { Technique_type_seq_pept , MI_TECH_seq_pept , "seq-pept" } ,
 { Technique_type_both , MI_TECH_both , "both" } ,
 { Technique_type_seq_pept_overlap , MI_TECH_seq_pept_overlap , "seq-pept-overlap" } ,
 { Technique_type_seq_pept_homol , MI_TECH_seq_pept_homol, "seq-pept-homol" } ,
 { Technique_type_concept_trans_a, MI_TECH_concept_trans_a, "concept-trans-a" } ,
 { Technique_type_htgs_1, MI_TECH_htgs_1, "HTGS-1" } ,
 { Technique_type_htgs_2, MI_TECH_htgs_2, "HTGS-2" } ,
 { Technique_type_htgs_3, MI_TECH_htgs_3, "HTGS-3" } ,
 { Technique_type_fli_cDNA, MI_TECH_fli_cdna, "fli-cDNA" } ,
 { Technique_type_htgs_0, MI_TECH_htgs_0, "HTGS-0" } ,
 { Technique_type_htc, MI_TECH_htc, "HTC" } ,
 { Technique_type_wgs, MI_TECH_wgs, "WGS" } ,
 { Technique_type_barcode, MI_TECH_barcode, "BARCODE" } ,
 { Technique_type_composite_wgs_htgs, MI_TECH_composite_wgs_htgs, "composite WGS-HTGS" } ,
 { Technique_type_tsa, MI_TECH_tsa, "TSA" } ,
 { Technique_type_other, MI_TECH_other, "other" } 
};


#define NUM_techniquetype_tech sizeof (techniquetype_tech) / sizeof (TechniqueTypeTechData)

NLM_EXTERN Int4 TechFromTechniqueType (Int4 technique_type) 
{
  Int4 i;

  for (i = 0; i < NUM_techniquetype_tech; i++) {
    if (techniquetype_tech[i].technique_type == technique_type) {
      return techniquetype_tech[i].tech;
    }
  }
  return -1;
}


NLM_EXTERN CharPtr TechNameFromTech (Int4 tech) 
{
  Int4 i;

  for (i = 0; i < NUM_techniquetype_tech; i++) {
    if (techniquetype_tech[i].tech == tech) {
      return techniquetype_tech[i].name;
    }
  }
  return NULL;
}


static Int4 TechFromTechName (CharPtr tech_name)
{
  Int4 i;

  for (i = 0; i < NUM_techniquetype_tech; i++) {
    if (StringCmp (techniquetype_tech[i].name, tech_name) == 0) {
      return techniquetype_tech[i].tech;
    }
  }
  return -1;
}


NLM_EXTERN ValNodePtr GetTechniqueTypeList (void)
{
  ValNodePtr list = NULL;
  Int4 i;

  for (i = 0; i < NUM_techniquetype_tech; i++) {
    ValNodeAddPointer (&list, techniquetype_tech[i].technique_type, StringSave (techniquetype_tech[i].name));
  }
  return list;
}


/* Completedness fields */
typedef struct completednesstypecompleteness {
  Int4 completedness_type;
  Int4 completeness;
  CharPtr name;
} CompletednessTypeCompletenessData, PNTR CompletednessTypeCompletenessPtr;

static CompletednessTypeCompletenessData completednesstype_completeness[] = {
 { Completedness_type_unknown, 0, "unknown" } ,
 { Completedness_type_complete, 1, "complete" } ,
 { Completedness_type_partial, 2, "partial" } ,
 { Completedness_type_no_left, 3, "no left" } ,
 { Completedness_type_no_right, 4, "no right" } ,
 { Completedness_type_no_ends, 5, "no ends" } ,
 { Completedness_type_has_left, 6, "has left" } ,
 { Completedness_type_has_right, 7, "has right" } ,
 { Completedness_type_other, 255, "other" }
};

#define NUM_completednesstype_completeness sizeof (completednesstype_completeness) / sizeof (CompletednessTypeCompletenessData)

NLM_EXTERN Int4 CompletenessFromCompletednessType (Int4 completedness_type) 
{
  Int4 i;

  for (i = 0; i < NUM_completednesstype_completeness; i++) {
    if (completednesstype_completeness[i].completedness_type == completedness_type) {
      return completednesstype_completeness[i].completeness;
    }
  }
  return -1;
}


NLM_EXTERN CharPtr CompletenessNameFromCompleteness (Int4 completeness) 
{
  Int4 i;

  for (i = 0; i < NUM_completednesstype_completeness; i++) {
    if (completednesstype_completeness[i].completeness == completeness) {
      return completednesstype_completeness[i].name;
    }
  }
  return NULL;
}


static Int4 CompletenessFromCompletenessName (CharPtr completeness_name)
{
  Int4 i;

  for (i = 0; i < NUM_completednesstype_completeness; i++) {
    if (StringCmp (completednesstype_completeness[i].name, completeness_name) == 0) {
      return completednesstype_completeness[i].completeness;
    }
  }
  return -1;
}


NLM_EXTERN ValNodePtr GetCompletednessTypeList (void)
{
  ValNodePtr list = NULL;
  Int4 i;

  for (i = 0; i < NUM_completednesstype_completeness; i++) {
    ValNodeAddPointer (&list, completednesstype_completeness[i].completedness_type, StringSave (completednesstype_completeness[i].name));
  }
  return list;
}


/* Molecule class fields */
typedef struct moleculeclasstypemol {
  Int4 moleculeclass_type;
  Int4 mol;
  CharPtr name;
} MoleculeClassTypeMolData, PNTR MoleculeClassTypeMolPtr;

static MoleculeClassTypeMolData moleculeclasstype_mol[] = {
 { Molecule_class_type_unknown, 0, "unknown" } ,
 { Molecule_class_type_dna, MOLECULE_CLASS_DNA, "DNA" } ,
 { Molecule_class_type_rna, MOLECULE_CLASS_RNA, "RNA" } ,
 { Molecule_class_type_protein, MOLECULE_CLASS_PROTEIN, "protein" } ,
 { Molecule_class_type_nucleotide, MOLECULE_CLASS_NUC, "nucleotide" } ,
 { Molecule_class_type_other, 255, "other" } 
};


#define NUM_moleculeclasstype_mol sizeof (moleculeclasstype_mol) / sizeof (MoleculeClassTypeMolData)

NLM_EXTERN Int4 MolFromMoleculeClassType (Int4 moleculeclass_type) 
{
  Int4 i;

  for (i = 0; i < NUM_moleculeclasstype_mol; i++) {
    if (moleculeclasstype_mol[i].moleculeclass_type == moleculeclass_type) {
      return moleculeclasstype_mol[i].mol;
    }
  }
  return -1;
}


NLM_EXTERN CharPtr MolNameFromMol (Int4 mol) 
{
  Int4 i;

  for (i = 0; i < NUM_moleculeclasstype_mol; i++) {
    if (moleculeclasstype_mol[i].mol == mol) {
      return moleculeclasstype_mol[i].name;
    }
  }
  return NULL;
}


static Int4 MolFromMolName (CharPtr mol_name)
{
  Int4 i;

  for (i = 0; i < NUM_moleculeclasstype_mol; i++) {
    if (StringCmp (moleculeclasstype_mol[i].name, mol_name) == 0) {
      return moleculeclasstype_mol[i].mol;
    }
  }
  return -1;
}


NLM_EXTERN ValNodePtr GetMoleculeClassTypeList (void)
{
  ValNodePtr list = NULL;
  Int4 i;

  for (i = 0; i < NUM_moleculeclasstype_mol; i++) {
    ValNodeAddPointer (&list, moleculeclasstype_mol[i].moleculeclass_type, StringSave (moleculeclasstype_mol[i].name));
  }
  return list;
}


/* Topology fields */
typedef struct topologytypetopology {
  Int4 topology_type;
  Int4 topology;
  CharPtr name;
} TopologyTypeTopologyData, PNTR TopologyTypeTopologyPtr;

static TopologyTypeTopologyData topologytype_topology[] = {
 { Topology_type_unknown, 0, "unknown" } ,
 { Topology_type_linear, TOPOLOGY_LINEAR, "linear" } ,
 { Topology_type_circular, TOPOLOGY_CIRCULAR, "circular" } ,
 { Topology_type_tandem, TOPOLOGY_TANDEM, "tandem" } ,
 { Topology_type_other, 255, "other" } 
};

#define NUM_topologytype_topology sizeof (topologytype_topology) / sizeof (TopologyTypeTopologyData)

NLM_EXTERN Int4 TopologyFromTopologyType (Int4 topology_type) 
{
  Int4 i;

  for (i = 0; i < NUM_topologytype_topology; i++) {
    if (topologytype_topology[i].topology_type == topology_type) {
      return topologytype_topology[i].topology;
    }
  }
  return -1;
}


NLM_EXTERN CharPtr TopologyNameFromTopology (Int4 topology) 
{
  Int4 i;

  for (i = 0; i < NUM_topologytype_topology; i++) {
    if (topologytype_topology[i].topology == topology) {
      return topologytype_topology[i].name;
    }
  }
  return NULL;
}


static Int4 TopologyFromTopologyName (CharPtr topology_name)
{
  Int4 i;

  for (i = 0; i < NUM_topologytype_topology; i++) {
    if (StringCmp (topologytype_topology[i].name, topology_name) == 0) {
      return topologytype_topology[i].topology;
    }
  }
  return -1;
}


NLM_EXTERN ValNodePtr GetTopologyTypeList (void)
{
  ValNodePtr list = NULL;
  Int4 i;

  for (i = 0; i < NUM_topologytype_topology; i++) {
    ValNodeAddPointer (&list, topologytype_topology[i].topology_type, StringSave (topologytype_topology[i].name));
  }
  return list;
}


/* strand fields */
typedef struct strandtypestrand {
  Int4 strand_type;
  Int4 strand;
  CharPtr name;
} StrandTypeStrandData, PNTR StrandTypeStrandPtr;

static StrandTypeStrandData strandtype_strand[] = {
 { Strand_type_unknown, 0, "unknown" } ,
 { Strand_type_single, STRANDEDNESS_SINGLE, "single" } ,
 { Strand_type_double__, STRANDEDNESS_DOUBLE, "double" } ,
 { Strand_type_mixed, 3, "mixed" } ,
 { Strand_type_mixed_rev, 4, "mixed-rev" } ,
 { Strand_type_other, 255, "other" } 
};

#define NUM_strandtype_strand sizeof (strandtype_strand) / sizeof (StrandTypeStrandData)

NLM_EXTERN Int4 StrandFromStrandType (Int4 strand_type) 
{
  Int4 i;

  for (i = 0; i < NUM_strandtype_strand; i++) {
    if (strandtype_strand[i].strand_type == strand_type) {
      return strandtype_strand[i].strand;
    }
  }
  return -1;
}


NLM_EXTERN CharPtr StrandNameFromStrand (Int4 strand) 
{
  Int4 i;

  for (i = 0; i < NUM_strandtype_strand; i++) {
    if (strandtype_strand[i].strand == strand) {
      return strandtype_strand[i].name;
    }
  }
  return NULL;
}


static Int4 StrandFromStrandName (CharPtr strand_name)
{
  Int4 i;

  for (i = 0; i < NUM_strandtype_strand; i++) {
    if (StringCmp (strandtype_strand[i].name, strand_name) == 0) {
      return strandtype_strand[i].strand;
    }
  }
  return -1;
}


NLM_EXTERN ValNodePtr GetStrandTypeList (void)
{
  ValNodePtr list = NULL;
  Int4 i;

  for (i = 0; i < NUM_strandtype_strand; i++) {
    ValNodeAddPointer (&list, strandtype_strand[i].strand_type, StringSave (strandtype_strand[i].name));
  }
  return list;
}


static CharPtr GetSequenceQualValName (ValNodePtr field)
{
  CharPtr val = NULL;

  if (field == NULL) return NULL;
  switch (field->choice) {
    case MolinfoField_molecule:
      val = BiomolNameFromBiomol (BiomolFromMoleculeType (field->data.intvalue));
      break;
    case MolinfoField_technique:
      val = TechNameFromTech (TechFromTechniqueType (field->data.intvalue));
      break;
    case MolinfoField_completedness:
      val = CompletenessNameFromCompleteness (CompletenessFromCompletednessType (field->data.intvalue));
      break;
    case MolinfoField_mol_class:
      val = MolNameFromMol (MolFromMoleculeClassType (field->data.intvalue));
      break;
    case MolinfoField_topology:
      val = TopologyNameFromTopology (TopologyFromTopologyType (field->data.intvalue));
      break;
    case MolinfoField_strand:
      val = StrandNameFromStrand (StrandFromStrandType (field->data.intvalue));
      break;
  }
  return val;
}


static CharPtr GetSequenceQualName (ValNodePtr field)
{
  CharPtr str = NULL, fieldname = "invalid field", val = "invalid value";
  CharPtr fmt = "%s %s";

  if (field == NULL) return NULL;
  switch (field->choice) {
    case MolinfoField_molecule:
      fieldname = "molecule";
      val = BiomolNameFromBiomol (BiomolFromMoleculeType (field->data.intvalue));
      break;
    case MolinfoField_technique:
      fieldname = "technique";
      val = TechNameFromTech (TechFromTechniqueType (field->data.intvalue));
      break;
    case MolinfoField_completedness:
      fieldname = "completeness";
      val = CompletenessNameFromCompleteness (CompletenessFromCompletednessType (field->data.intvalue));
      break;
    case MolinfoField_mol_class:
      fieldname = "class";
      val = MolNameFromMol (MolFromMoleculeClassType (field->data.intvalue));
      break;
    case MolinfoField_topology:
      fieldname = "topology";
      val = TopologyNameFromTopology (TopologyFromTopologyType (field->data.intvalue));
      break;
    case MolinfoField_strand:
      fieldname = "strand";
      val = StrandNameFromStrand (StrandFromStrandType (field->data.intvalue));
      break;
  }
  if (val == NULL) {
    val = "Invalid value";
  }
  str = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) + StringLen (fieldname) + StringLen (val)));
  sprintf (str, fmt, fieldname, val);
  return str;
}


/* Simple constraints */
static Boolean IsWholeWordMatch (CharPtr start, CharPtr found, Int4 match_len)
{
  Boolean rval = TRUE;
  Char    char_after;
  Char    char_before;
  
  if (match_len == 0)
  {
    rval = TRUE;
  }
  else if (start == NULL || found == NULL)
  {
    rval = FALSE;
  }
  else
  {
	  char_after = *(found + match_len);
    if (found != start)
	  {
	    char_before = *(found - 1);
	    if (isalpha ((Int4) char_before) || isdigit ((Int4) char_before))
	    {
	      rval = FALSE;
	    }
	  }
	  if (char_after != 0 && (isalpha ((Int4) char_after) || isdigit ((Int4)char_after)))
	  {
	    rval = FALSE;
	  }   
  }
  return rval;
}


NLM_EXTERN Boolean IsStringConstraintEmpty (StringConstraintPtr scp)
{
  if (scp == NULL || StringHasNoText (scp->match_text)) return TRUE;
  else return FALSE;
}


NLM_EXTERN Boolean DoesSingleStringMatchConstraint (CharPtr str, StringConstraintPtr scp)
{
  CharPtr pFound;
  Boolean rval = FALSE;
  Char    char_after = 0;
  
  if (IsStringConstraintEmpty (scp)) return TRUE;
  if (StringHasNoText (str)) return FALSE;

  switch (scp->match_location) 
  {
    case String_location_contains:
	    if (scp->case_sensitive)
	    {
	      pFound = StringSearch (str, scp->match_text);
	    }
	    else
	    {
	      pFound = StringISearch (str, scp->match_text);
	    }
      if (pFound == NULL) 
      {
        rval = FALSE;
      }
      else if (scp->whole_word) 
      {
        rval = IsWholeWordMatch (str, pFound, StringLen (scp->match_text));
        while (!rval && pFound != NULL) 
        {
	        if (scp->case_sensitive)
	        {
	          pFound = StringSearch (pFound + 1, scp->match_text);
	        }
	        else
	        {
	          pFound = StringISearch (pFound + 1, scp->match_text);
	        }
          if (pFound != NULL)
          {
            rval = IsWholeWordMatch (str, pFound, StringLen (scp->match_text));
          }
        }
      }
      else
      {
        rval = TRUE;
      }
      break;
    case String_location_starts:
	    if (scp->case_sensitive)
	    {
	      pFound = StringSearch (str, scp->match_text);
	    }
	    else
	    {
	      pFound = StringISearch (str, scp->match_text);
	    }
      if (pFound == str)
      {
        if (scp->whole_word) 
        {
          rval = IsWholeWordMatch (str, pFound, StringLen (scp->match_text));
        }
        else
        {
          rval = TRUE;
        }
      }
      break;
    case String_location_ends:
	    if (scp->case_sensitive)
	    {
	      pFound = StringSearch (str, scp->match_text);
	    }
	    else
	    {
	      pFound = StringISearch (str, scp->match_text);
	    }
      while (pFound != NULL && !rval) {
  	    char_after = *(pFound + StringLen (scp->match_text));
        if (char_after == 0)
        {
          if (scp->whole_word) 
          {
            rval = IsWholeWordMatch (str, pFound, StringLen (scp->match_text));
          }
          else
          {
            rval = TRUE;
          }
          /* stop the search, we're at the end of the string */
          pFound = NULL;
        }
        else
        {
	        if (scp->case_sensitive)
	        {
	          pFound = StringSearch (pFound + 1, scp->match_text);
	        }
	        else
	        {
	          pFound = StringISearch (pFound + 1, scp->match_text);
	        }
        }
      }
      break;
    case String_location_equals:
      if (scp->case_sensitive) 
      {
        if (StringCmp (str, scp->match_text) == 0) 
        {
          rval = TRUE;
        }
      }
      else
      {
        if (StringICmp (str, scp->match_text) == 0) 
        {
          rval = TRUE;
        }
      }
      break;
    case String_location_inlist:
	    if (scp->case_sensitive)
	    {
	      pFound = StringSearch (scp->match_text, str);
	    }
	    else
	    {
	      pFound = StringISearch (scp->match_text, str);
	    }
      if (pFound == NULL) 
      {
        rval = FALSE;
      }
      else
      {
        rval = IsWholeWordMatch (scp->match_text, pFound, StringLen (str));
        while (!rval && pFound != NULL) 
        {
	        if (scp->case_sensitive)
	        {
	          pFound = StringSearch (pFound + 1, str);
	        }
	        else
	        {
	          pFound = StringISearch (pFound + 1, str);
	        }
          if (pFound != NULL)
          {
            rval = IsWholeWordMatch (scp->match_text, pFound, StringLen (str));
          }
        }
      }
      if (!rval) {
        /* look for spans */
        rval = IsStringInSpanInList (str, scp->match_text);
      }
      break;
	}
	return rval;
}


NLM_EXTERN Boolean DoesStringMatchConstraint (CharPtr str, StringConstraintPtr scp)
{
  Boolean rval;

  rval = DoesSingleStringMatchConstraint (str, scp);
  if (scp != NULL && scp->not_present) {
    rval = !rval;
  }
  return rval;
}


static Boolean DoesStringListMatchConstraint (ValNodePtr list, StringConstraintPtr scp)
{
  Int4 len = 1;
  CharPtr tmp;
  Boolean rval = FALSE;
  ValNodePtr vnp;

  if (IsStringConstraintEmpty (scp)) {
    return TRUE;
  }
  if (list == NULL) return FALSE;

  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    len += StringLen (vnp->data.ptrvalue) + 2;
  }

  tmp = (CharPtr) MemNew (sizeof (Char) * len);
  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    StringCat (tmp, vnp->data.ptrvalue);
    if (vnp->next != NULL) {
      StringCat (tmp, "; ");
    }
  }

  rval = DoesStringMatchConstraint (tmp, scp);
  tmp = MemFree (tmp);
  return rval;  
}


static Boolean DoesStrandMatchConstraint (SeqLocPtr slp, LocationConstraintPtr lcp)
{
  Uint2 strand;
  Boolean rval = FALSE;
  
  if (slp == NULL)
  {
    rval = FALSE;
  }
  else if (lcp == NULL || lcp->strand == Strand_constraint_any)
  {
    rval = TRUE;
  }
  else
  {
    strand = SeqLocStrand (slp);
    if (strand == Seq_strand_minus)
    {
      if (lcp->strand == Strand_constraint_minus)
      {
        rval = TRUE;
      }
      else
      {
        rval = FALSE;
      }
    }
    else
    {
      if (lcp->strand == Strand_constraint_plus)
      {
        rval = TRUE;
      }
      else
      {
        rval = FALSE;
      }
    }
  }
  return rval;
}


static Boolean DoesBioseqMatchSequenceType (BioseqPtr bsp, Uint2 seq_type)
{
  Boolean rval = FALSE;

  if (bsp == NULL) return FALSE;
  if (seq_type == Seqtype_constraint_any) return TRUE;

  if (ISA_na (bsp->mol) && seq_type == Seqtype_constraint_nuc)
  {
    rval = TRUE;
  }
  else if (ISA_aa (bsp->mol) && seq_type == Seqtype_constraint_prot)
  {
    rval = TRUE;
  }
  return rval;
}


static Boolean DoesSequenceTypeMatchContraint (SeqLocPtr slp, LocationConstraintPtr lcp)
{
  Boolean   rval = FALSE;
  BioseqPtr bsp;
  
  if (slp == NULL)
  {
    rval = FALSE;
  }
  else if (lcp == NULL || lcp->seq_type == Seqtype_constraint_any)
  {
    rval = TRUE;
  }
  else
  {
    bsp = BioseqFindFromSeqLoc (slp);
    rval = DoesBioseqMatchSequenceType (bsp, lcp->seq_type);
  }
  return rval;
}

static Boolean DoesLocationMatchConstraint (SeqLocPtr slp, LocationConstraintPtr lcp)

{
  Boolean rval = FALSE;
  
  if (slp == NULL)
  {
    rval = FALSE;
  }  
  else if (lcp == NULL || (DoesStrandMatchConstraint (slp, lcp) && DoesSequenceTypeMatchContraint (slp, lcp)))
  {
    rval = TRUE;
  }
  return rval; 
}


static Boolean DoesObjectMatchLocationConstraint (Uint1 choice, Pointer data, LocationConstraintPtr constraint)
{
  SeqFeatPtr  sfp;
  SeqDescrPtr sdp;
  CGPSetPtr   cgp;
  BioseqPtr  bsp = NULL;
  BioseqSetPtr bssp;
  ValNodePtr    vnp;
  ObjValNodePtr ovp;
  SeqMgrFeatContext context;

  if (data == NULL) return FALSE;
 
  if (constraint == NULL 
      || (constraint->strand == Strand_constraint_any
          && constraint->seq_type == Seqtype_constraint_any)) {
    return TRUE;
  }

  if (choice == OBJ_SEQFEAT) {
    sfp = (SeqFeatPtr) data;
    bsp = BioseqFindFromSeqLoc (sfp->location);
  } else if (choice == OBJ_SEQDESC) {
    sdp = (SeqDescrPtr) data;
    if (sdp->extended != 0) {
      ovp = (ObjValNodePtr) sdp;
      if (ovp->idx.parenttype == OBJ_BIOSEQSET) {
        bssp = (BioseqSetPtr) ovp->idx.parentptr;
        if (bssp != NULL && bssp->seq_set != NULL && IS_Bioseq_set (bssp->seq_set)) {
          bsp = (BioseqPtr) bssp->seq_set->data.ptrvalue;
        }
      } else if (ovp->idx.parenttype == OBJ_BIOSEQ) {
        bsp = (BioseqPtr) ovp->idx.parentptr;
      }
    }
  } else if (choice == 0) {
    if (constraint->seq_type != Seqtype_constraint_any) {
      return FALSE;
    }
    cgp = (CGPSetPtr) data;
    if (cgp->cds_list != NULL && cgp->cds_list->data.ptrvalue != NULL) {
      sfp = (SeqFeatPtr) cgp->cds_list->data.ptrvalue;
      bsp = BioseqFindFromSeqLoc (sfp->location);
    } else if (cgp->gene_list != NULL && cgp->gene_list->data.ptrvalue != NULL) {
      sfp = (SeqFeatPtr) cgp->gene_list->data.ptrvalue;
      bsp = BioseqFindFromSeqLoc (sfp->location);
    } else if (cgp->mrna_list != NULL && cgp->mrna_list->data.ptrvalue != NULL) {
      sfp = (SeqFeatPtr) cgp->mrna_list->data.ptrvalue;
      bsp = BioseqFindFromSeqLoc (sfp->location);
    } else if (cgp->prot_list != NULL && cgp->prot_list->data.ptrvalue != NULL) {
      sfp = (SeqFeatPtr) cgp->prot_list->data.ptrvalue;
      bsp = BioseqFindFromSeqLoc (sfp->location);
    }
  }
  if (!DoesBioseqMatchSequenceType(bsp, constraint->seq_type)) {
    return FALSE;
  }
  if (constraint->strand != Strand_constraint_any && ISA_aa (bsp->mol)) {      
    sfp = SeqMgrGetCDSgivenProduct (bsp, &context);
    if (constraint->strand == Strand_constraint_minus && context.strand != Seq_strand_minus) {
      return FALSE;
    }
    if (constraint->strand == Strand_constraint_plus && context.strand == Seq_strand_minus) {
      return FALSE;
    }
  } else if (constraint->strand != Strand_constraint_any) {
    if (choice == 0) {
      /* strand for CDS-Gene-Prot group */
      cgp = (CGPSetPtr) data;
      for (vnp = cgp->cds_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && !DoesStrandMatchConstraint (sfp->location, constraint)) {
          return FALSE;
        }
      }
      for (vnp = cgp->gene_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && !DoesStrandMatchConstraint (sfp->location, constraint)) {
          return FALSE;
        }
      }
      for (vnp = cgp->mrna_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && !DoesStrandMatchConstraint (sfp->location, constraint)) {
          return FALSE;
        }
      }
    } else if (choice == OBJ_SEQFEAT) {
      sfp = (SeqFeatPtr) data;
      if (!DoesStrandMatchConstraint (sfp->location, constraint)) {
        return FALSE;
      }
    } else {
      /* descriptors can't meet strand constraints */
      return FALSE;
    }
  }
  return TRUE;
}


/* for parsing and editing */
static CharPtr GetTextPortionFromString (CharPtr str, TextPortionPtr text_portion)
{
  CharPtr portion = NULL;
  CharPtr found_start, found_end;
  Int4    found_len;

  if (StringHasNoText (str)) {
    return NULL;
  }
  if (text_portion == NULL) {
    return StringSave (str);
  }  
  
  if (text_portion->left_text == NULL || text_portion->left_text [0] == 0)
  {
    found_start = str;
  }
  else
  {
    if (text_portion->case_sensitive)
    {
      found_start = StringSearch (str, text_portion->left_text);
    }
    else
    {
      found_start = StringISearch (str, text_portion->left_text);
    }
    
    if (text_portion->whole_word && ! IsWholeWordMatch (str, found_start, StringLen (text_portion->left_text)))
    {
      found_start = NULL;
    }
  }
  
  if (found_start == NULL)
  {
    return NULL;
  }
  
  if (!text_portion->include_left)
  {
    found_start += StringLen (text_portion->left_text);
  }
  
  if (text_portion->right_text == NULL || text_portion->right_text [0] == 0)
  {
    found_len = StringLen (found_start);
  }
  else
  {
    if (text_portion->case_sensitive)
    {
      found_end = StringSearch (found_start, text_portion->right_text);
    }
    else
    {
      found_end = StringISearch (found_start, text_portion->right_text);
    }
    if (text_portion->whole_word && ! IsWholeWordMatch (str, found_end, StringLen (text_portion->right_text)))
    {
      found_end = NULL;
    }    
    
    if (found_end == NULL)
    {
      found_len = 0;
    }
    else if (text_portion->include_right)
    {
      found_len = (Int4)(found_end - found_start) + StringLen (text_portion->right_text);
    }
    else
    {
      found_len = found_end - found_start;
    }
  }

  if (found_len > 0)
  {
    portion = (CharPtr) MemNew (sizeof (Char) * (found_len + 1));
    StringNCpy (portion, found_start, found_len);
    portion[found_len] = 0;
  }
  return portion;
}



static CharPtr FindTextPortionLocationInString (CharPtr str, TextPortionPtr text_portion)
{
  CharPtr start, stop;

  if (str == NULL || text_portion == NULL) return FALSE;

  if (text_portion->left_text != NULL) {
    start = StringSearch (str, text_portion->left_text);
    if (start != NULL) {
      if (!text_portion->include_left) {
        start += StringLen (text_portion->left_text);
      }
    }
  } else {
    start = str;
  }
  if (start != NULL) {
    if (text_portion->right_text != NULL) { 
      stop = StringSearch (start, text_portion->right_text);
      if (stop == NULL) {
        start = NULL;
      }
    }
  }
  return start;
}


static void ReplaceStringForParse(CharPtr src_text, TextPortionPtr text_portion)
{
  CharPtr         src, dst;
  
  if (src_text == NULL || text_portion == NULL) {
    return;
  }

  dst = FindTextPortionLocationInString (src_text, text_portion);
  if (dst == NULL) return;
  if (text_portion->right_text == NULL) {
    *dst = 0;
  } else {
    src = StringSearch (src_text, text_portion->right_text);
    if (src != NULL) {
      if (text_portion->include_right) {
        src += StringLen (text_portion->right_text);
      }
      while (*src != 0) {
        *dst = *src;
        dst++;
        src++;
      }
      *dst = 0;
    }
  }
}


/* generic functions for getting string values */
static Int4 GetDbtagStringLen (DbtagPtr db_tag)
{
  Int4 len;
  
  if (db_tag == NULL)
  {
    return 0;
  }
  
  len = StringLen (db_tag->db) + 2;
  if (db_tag->tag != NULL)
  {
    if (db_tag->tag->str != NULL)
    {
      len += StringLen (db_tag->tag->str);
    }
    else
    {
      len += 10;
    }
  }
  return len;
}


static CharPtr GetDbtagString (DbtagPtr db_tag)
{
  Int4    len;
  CharPtr str;
  
  if (db_tag == NULL) {
    return NULL;
  }
  
  len = GetDbtagStringLen (db_tag);
  if (len == 0) {
    return NULL;
  }
  
  str = (CharPtr) MemNew (len * sizeof (Char));
  if (str != NULL) {
    StringCpy (str, db_tag->db);
    StringCat (str, ":");
    if (db_tag->tag != NULL) {
      if (db_tag->tag->str != NULL) {
        StringCat (str, db_tag->tag->str);
      } else {
        sprintf (str + StringLen (str), "%d", db_tag->tag->id);
      }
    }
  }
  return str;
}


/* generic functions for setting field values */
static Boolean SetStringValue (CharPtr PNTR existing_val, CharPtr new_val, Uint2 existing_text)
{
  Boolean rval = FALSE;
  Int4 len;
  CharPtr tmp;

  if (existing_val == NULL) {
    return FALSE;
  }

  if (StringHasNoText (*existing_val)) {
    *existing_val = MemFree (*existing_val);
    *existing_val = StringSave (new_val);
    rval = TRUE;
  } else {
    switch (existing_text) {
      case ExistingTextOption_replace_old :
        *existing_val = MemFree (*existing_val);
        *existing_val = StringSave (new_val);
        rval = TRUE;
        break;
      case ExistingTextOption_append_semi :
        len = StringLen (new_val) + StringLen (*existing_val) + 3;
        tmp = (CharPtr) MemNew (sizeof (Char) * len);
        if (tmp != NULL) {
          sprintf (tmp, "%s; %s", *existing_val, new_val);
          MemFree (*existing_val);
          *existing_val = tmp;
          rval = TRUE;
        }
        break;
      case ExistingTextOption_append_space :
        len = StringLen (new_val) + StringLen (*existing_val) + 2;
        tmp = (CharPtr) MemNew (sizeof (Char) * len);
        if (tmp != NULL) {
          sprintf (tmp, "%s %s", *existing_val, new_val);
          MemFree (*existing_val);
          *existing_val = tmp;
          rval = TRUE;
        }
        break;
      case ExistingTextOption_append_colon :
        len = StringLen (new_val) + StringLen (*existing_val) + 3;
        tmp = (CharPtr) MemNew (sizeof (Char) * len);
        if (tmp != NULL) {
          sprintf (tmp, "%s: %s", *existing_val, new_val);
          MemFree (*existing_val);
          *existing_val = tmp;
          rval = TRUE;
        }
        break;
      case ExistingTextOption_append_none :
        len = StringLen (new_val) + StringLen (*existing_val) + 1;
        tmp = (CharPtr) MemNew (sizeof (Char) * len);
        if (tmp != NULL) {
          sprintf (tmp, "%s%s", *existing_val, new_val);
          MemFree (*existing_val);
          *existing_val = tmp;
          rval = TRUE;
        }
        break;
      case ExistingTextOption_prefix_semi :
        len = StringLen (new_val) + StringLen (*existing_val) + 3;
        tmp = (CharPtr) MemNew (sizeof (Char) * len);
        if (tmp != NULL) {
          sprintf (tmp, "%s; %s", new_val, *existing_val);
          MemFree (*existing_val);
          *existing_val = tmp;
          rval = TRUE;
        }
        break;
      case ExistingTextOption_prefix_space :
        len = StringLen (new_val) + StringLen (*existing_val) + 2;
        tmp = (CharPtr) MemNew (sizeof (Char) * len);
        if (tmp != NULL) {
          sprintf (tmp, "%s %s", new_val, *existing_val);
          MemFree (*existing_val);
          *existing_val = tmp;
          rval = TRUE;
        }
        break;
      case ExistingTextOption_prefix_colon :
        len = StringLen (new_val) + StringLen (*existing_val) + 3;
        tmp = (CharPtr) MemNew (sizeof (Char) * len);
        if (tmp != NULL) {
          sprintf (tmp, "%s: %s", new_val, *existing_val);
          MemFree (*existing_val);
          *existing_val = tmp;
          rval = TRUE;
        }
        break;
      case ExistingTextOption_prefix_none :
        len = StringLen (new_val) + StringLen (*existing_val) + 1;
        tmp = (CharPtr) MemNew (sizeof (Char) * len);
        if (tmp != NULL) {
          sprintf (tmp, "%s%s", new_val, *existing_val);
          MemFree (*existing_val);
          *existing_val = tmp;
          rval = TRUE;
        }
        break;
      case ExistingTextOption_leave_old :
        rval = FALSE;
    }
  }
  return rval;
}


static Boolean SetStringsInValNodeStringList (ValNodePtr PNTR list, StringConstraintPtr scp, CharPtr new_val, Uint2 existing_text)
{
  ValNodePtr vnp;
  CharPtr    cp;
  Boolean rval = FALSE;
  
  if (list == NULL)
  {
    return FALSE;
  }

  if (*list == NULL && (scp == NULL || StringHasNoText (scp->match_text))) {
    ValNodeAddPointer (list, 0, StringSave (new_val));
    rval = TRUE;
  } else if (existing_text == ExistingTextOption_append_semi) {
    if (DoesStringListMatchConstraint (*list, scp)) {
      ValNodeAddPointer (list, 0, StringSave (new_val));
      rval = TRUE;
    }
  } else if (existing_text == ExistingTextOption_prefix_semi) {
    if (DoesStringListMatchConstraint (*list, scp)) {
      vnp = ValNodeNew (NULL);
      vnp->data.ptrvalue = StringSave (new_val);
      vnp->next = *list;
      *list = vnp;
      rval = TRUE;
    }
  } else if (existing_text == ExistingTextOption_replace_old) {
    if (DoesStringListMatchConstraint (*list, scp)) {
      *list = ValNodeFreeData (*list);
      vnp = ValNodeNew (NULL);
      vnp->data.ptrvalue = StringSave (new_val);
      *list = vnp;
      rval = TRUE;
    }
  } else if (existing_text == ExistingTextOption_leave_old) {
    rval = FALSE;
  } else {
    for (vnp = *list; vnp != NULL; vnp = vnp->next)
    {
      cp = (CharPtr) vnp->data.ptrvalue;
      if (DoesStringMatchConstraint (cp, scp)) {
        rval |= SetStringValue (&cp, new_val, existing_text);
        vnp->data.ptrvalue = cp;
      }
    }
  }
  return rval;
}


static Boolean SetStringInGBQualList (GBQualPtr PNTR list, ValNodePtr field, StringConstraintPtr scp, CharPtr new_val, Uint2 existing_text)
{
  Boolean rval = FALSE;
  Int4 gbqual;
  CharPtr qual_name = NULL;
  GBQualPtr gbq, last_gbq = NULL;

  if (field == NULL) return FALSE;

  if (field->choice == FeatQualChoice_legal_qual) 
  {
    gbqual = GetGBQualFromFeatQual (field->data.intvalue);
    if (gbqual > -1) {
      qual_name = ParFlat_GBQual_names [gbqual].name;
      for (gbq = *list; gbq != NULL; gbq = gbq->next) {
        if (StringCmp (gbq->qual, qual_name) == 0
            && DoesStringMatchConstraint (gbq->val, scp)) {
          rval |= SetStringValue (&(gbq->val), new_val, existing_text);
        }
        last_gbq = gbq;
      }
      if (!rval && (scp == NULL || scp->match_text == NULL)) {
        gbq = GBQualNew ();
        gbq->qual = StringSave (qual_name);
        gbq->val = StringSave (new_val);
        if (last_gbq == NULL) {
          *list = gbq;
        } else {
          last_gbq->next = gbq;
        }
        rval = TRUE;
      }
    }
  } else if (field->choice == FeatQualChoice_illegal_qual) {
    for (gbq = *list; gbq != NULL; gbq = gbq->next) {
      if (DoesStringMatchConstraint (gbq->qual, field->data.ptrvalue)
          && DoesStringMatchConstraint (gbq->val, scp)) {
        rval |= SetStringValue (&(gbq->val), new_val, existing_text);
      }
    }
  }

  return rval;
}


static Boolean IsAllDigits (CharPtr str)
{
  CharPtr cp;

  if (StringHasNoText (str)) return FALSE;

  cp = str;
  while (*cp != 0 && isdigit (*cp)) {
    cp++;
  }
  if (*cp == 0) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static Boolean SetObjectIdString (ObjectIdPtr oip, CharPtr value, Uint2 existing_text)
{
  Boolean rval = FALSE;
  Char    num[15];
  CharPtr tmp = NULL;

  if (oip == NULL) {
    return FALSE;
  }

  if (oip->id > 0) {
    sprintf (num, "%d", oip->id);
    tmp = StringSave (num);
  } else {
    tmp = StringSaveNoNull (oip->str);
  }
  if (SetStringValue (&tmp, value, existing_text)) {
    oip->str = MemFree (oip->str);        
    oip->id = 0;
    if (IsAllDigits (tmp)) {
      oip->id = atoi (tmp);
    } else {
      oip->str = tmp;
      tmp = NULL;
    }
    rval = TRUE;
  }
  tmp = MemFree (tmp);
  return rval;
}


static Boolean SetDbtagString (DbtagPtr db_tag, CharPtr value, Uint2 existing_text)
{
  Boolean rval = FALSE;
  CharPtr cp;
  Int4    dbxvalid;
  CharPtr tmp;
  CharPtr twoval;
  
  if (db_tag == NULL || StringHasNoText (value)) {
    return FALSE;
  }

  cp = StringChr (value, ':');
  if (cp == NULL) {
    tmp = StringSave (db_tag->db);
    if (SetStringValue (&tmp, value, existing_text)) {
      dbxvalid = IsDbxrefValid (tmp, NULL, NULL, TRUE, NULL);
      if (dbxvalid != 0) {
        db_tag->db = MemFree (db_tag->db);
        db_tag->db = tmp;
        tmp = NULL;
        rval = TRUE;
      }
    }
    if (!rval) {
      if (db_tag->tag == NULL) {
        db_tag->tag = ObjectIdNew();
      }
      rval = SetObjectIdString (db_tag->tag, value, existing_text);
    }
    tmp = MemFree (tmp);
  } else {
    twoval = StringSave (value);
    cp = StringChr (twoval, ':');
    *cp = 0;
    cp++;
    rval = SetStringValue (&(db_tag->db), twoval, existing_text);
    if (db_tag->tag == NULL) {
      db_tag->tag = ObjectIdNew ();
    }
    rval |= SetObjectIdString (db_tag->tag, cp, existing_text);
    twoval = MemFree (twoval);
  }
  return rval;
}


static Boolean SetDbxrefString (SeqFeatPtr sfp, StringConstraintPtr scp, CharPtr value, Uint2 existing_text)
{
  ValNodePtr vnp;
  Boolean    rval = FALSE, skip;
  DbtagPtr   dbtag;
  CharPtr    cp;
  
  if (sfp == NULL) {
    return FALSE;
  }

  if ((sfp->dbxref == NULL || existing_text == ExistingTextOption_append_semi) && (scp == NULL || StringHasNoText (scp->match_text))) {
    dbtag = DbtagNew ();
    rval = SetDbtagString (dbtag, value, existing_text);
    if (rval) {
      ValNodeAddPointer (&(sfp->dbxref), 0, dbtag);
    } else {
      dbtag = DbtagFree (dbtag);
    }
  } else {
    for (vnp = sfp->dbxref; vnp != NULL; vnp = vnp->next) {
      skip = FALSE;
      if (scp != NULL) {
        cp = GetDbtagString (vnp->data.ptrvalue);
        if (!DoesStringMatchConstraint (cp, scp)) {
          skip = TRUE;
        }
        cp = MemFree (cp);
      }
      if (!skip) {
        rval |= SetDbtagString (vnp->data.ptrvalue, value, existing_text);
      }
    }
  }
  return rval;
}



static CharPtr GetFirstValNodeStringMatch (ValNodePtr vnp, StringConstraintPtr scp)
{
  CharPtr str = NULL;
  while (vnp != NULL && str == NULL) {
    if (!StringHasNoText (vnp->data.ptrvalue)
        && DoesStringMatchConstraint (vnp->data.ptrvalue, scp)) {
      str = StringSave (vnp->data.ptrvalue);
    } 
    vnp = vnp->next;
  }
  return str;
}


static Boolean RemoveValNodeStringMatch (ValNodePtr PNTR list, StringConstraintPtr scp)
{
  ValNodePtr vnp_prev = NULL, vnp_next, vnp;
  Boolean    rval = FALSE;

  if (list == NULL) return FALSE;
  vnp = *list;
  while (vnp != NULL) {
    vnp_next = vnp->next;
    if (!StringHasNoText (vnp->data.ptrvalue) 
        && DoesStringMatchConstraint (vnp->data.ptrvalue, scp)) {
      if (vnp_prev == NULL) {
        *list = vnp->next;
      } else {
        vnp_prev->next = vnp->next;
      }
      vnp->next = NULL;
      vnp = ValNodeFreeData (vnp);
      rval = TRUE;
    } else {
      vnp_prev = vnp;
    }
    vnp = vnp_next;
  }
  return rval;
}


static CharPtr GetFirstGBQualMatch (GBQualPtr qual, CharPtr qual_name, StringConstraintPtr scp)
{
  CharPtr str = NULL;
  while (qual != NULL && str == NULL) {
    if (StringICmp (qual->qual, qual_name) == 0
        &&!StringHasNoText (qual->val)
        && DoesStringMatchConstraint (qual->val, scp)) {
      str = StringSave (qual->val);
    } 
    qual = qual->next;
  }
  return str;
}


static CharPtr GetFirstGBQualMatchConstraintName (GBQualPtr qual, StringConstraintPtr qual_name, StringConstraintPtr scp)
{
  CharPtr str = NULL;
  while (qual != NULL && str == NULL) {
    if (DoesStringMatchConstraint (qual->qual, qual_name)
        &&!StringHasNoText (qual->val)
        && DoesStringMatchConstraint (qual->val, scp)) {
      str = StringSave (qual->val);
    } 
    qual = qual->next;
  }
  return str;
}


static Boolean RemoveGBQualMatch (GBQualPtr PNTR list, CharPtr qual_name, StringConstraintPtr scp)
{
  GBQualPtr qual_prev = NULL, qual_next, qual;
  Boolean   rval = FALSE;

  if (list == NULL) return FALSE;
  qual = *list;
  while (qual != NULL) {
    qual_next = qual->next;
    if (StringICmp (qual->qual, qual_name) == 0
        && !StringHasNoText (qual->val) 
        && DoesStringMatchConstraint (qual->val, scp)) {
      if (qual_prev == NULL) {
        *list = qual->next;
      } else {
        qual_prev->next = qual->next;
      }
      qual->next = NULL;
      qual = GBQualFree (qual);
      rval = TRUE;
    } else {
      qual_prev = qual;
    }
    qual = qual_next;
  }
  return rval;
}


static Boolean RemoveGBQualMatchConstraintName (GBQualPtr PNTR list, StringConstraintPtr qual_name, StringConstraintPtr scp)
{
  GBQualPtr qual_prev = NULL, qual_next, qual;
  Boolean   rval = FALSE;

  if (list == NULL) return FALSE;
  qual = *list;
  while (qual != NULL) {
    qual_next = qual->next;
    if (DoesStringMatchConstraint (qual->qual, qual_name)
        && !StringHasNoText (qual->val) 
        && DoesStringMatchConstraint (qual->val, scp)) {
      if (qual_prev == NULL) {
        *list = qual->next;
      } else {
        qual_prev->next = qual->next;
      }
      qual->next = NULL;
      qual = GBQualFree (qual);
      rval = TRUE;
    } else {
      qual_prev = qual;
    }
    qual = qual_next;
  }
  return rval;
}


static CharPtr GetDbxrefString (SeqFeatPtr sfp, StringConstraintPtr scp)
{
  ValNodePtr vnp;
  Int4       len = 0;
  CharPtr    str = NULL, cp;
  
  if (sfp == NULL || sfp->dbxref == NULL) {
    return NULL;
  }
  
  for (vnp = sfp->dbxref; vnp != NULL; vnp = vnp->next) {
    cp = GetDbtagString (vnp->data.ptrvalue);
    if (cp != NULL && DoesStringMatchConstraint(cp, scp)) {
      len += StringLen (cp) + 1;
    }
    cp = MemFree (cp);
  }
  
  if (len == 0) {
    return NULL;
  }
  
  str = (CharPtr) MemNew ((len + 1) * sizeof (Char));
  if (str != NULL) {
    for (vnp = sfp->dbxref; vnp != NULL; vnp = vnp->next) {
      cp = GetDbtagString (vnp->data.ptrvalue);
      if (cp != NULL && DoesStringMatchConstraint(cp, scp)) {
        StringCat (str, cp);
        StringCat (str, ";");
      }
      cp = MemFree (cp);
    }
  }
  if (StringLen (str) >1) {
    /* remove final semicolon */
    str [StringLen (str) - 2] = 0;
  }
  return str;
}


static Boolean RemoveDbxrefString (SeqFeatPtr sfp, StringConstraintPtr scp)
{
  ValNodePtr vnp, vnp_prev = NULL, vnp_next;
  CharPtr    cp;
  Boolean    rval = FALSE;
  
  if (sfp == NULL || sfp->dbxref == NULL) {
    return FALSE;
  }
  
  vnp = sfp->dbxref;
  while (vnp != NULL) {
    vnp_next = vnp->next;
    cp = GetDbtagString (vnp->data.ptrvalue);
    if (DoesStringMatchConstraint(cp, scp)) {
      if (vnp_prev == NULL) {
        sfp->dbxref = vnp->next;
      } else {
        vnp_prev->next = vnp->next;
      }
      vnp->next = NULL;
      vnp->data.ptrvalue = DbtagFree (vnp->data.ptrvalue);
      vnp = ValNodeFree (vnp);
      rval = TRUE;
    } else {
      vnp_prev = vnp;
    }
  }
  return rval;  
}


static CharPtr GetRNAProductString (SeqFeatPtr sfp, StringConstraintPtr scp)
{
  RnaRefPtr  rrp;
  SeqMgrFeatContext context;
  CharPtr    str = NULL;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA || sfp->data.value.ptrvalue == NULL) {
    return NULL;
  }

  rrp = sfp->data.value.ptrvalue;
  if (rrp->ext.choice == 0 
      || (rrp->ext.choice == 1 && StringHasNoText (rrp->ext.value.ptrvalue))
      || (rrp->ext.choice == 1 
          && (StringCmp (rrp->ext.value.ptrvalue, "ncRNA") == 0 
              || StringCmp (rrp->ext.value.ptrvalue, "tmRNA") == 0
              || StringCmp (rrp->ext.value.ptrvalue, "misc_RNA") == 0))) {
    str = GetFirstGBQualMatch (sfp->qual, "product", scp);
  }

  if (str == NULL) {
    if (rrp->ext.choice == 1 && !StringHasNoText (rrp->ext.value.ptrvalue)
        && StringCmp (rrp->ext.value.ptrvalue, "ncRNA") != 0
        && StringCmp (rrp->ext.value.ptrvalue, "tmRNA") != 0
        && StringCmp (rrp->ext.value.ptrvalue, "misc_RNA") != 0) {
      str = StringSave (rrp->ext.value.ptrvalue);        
    } else if (rrp->ext.choice == 2 && rrp->ext.value.ptrvalue != NULL) {
      if (SeqMgrGetDesiredFeature (sfp->idx.entityID, NULL, 0, 0, sfp, &context) != NULL
          && !StringHasNoText (context.label)
          && StringCmp (context.label, "tRNA") != 0) {
        str = (CharPtr) MemNew (sizeof (Char) + (StringLen (context.label) + 6));
        sprintf (str, "tRNA-%s", context.label);
      }
    }
    if (!DoesStringMatchConstraint(str, scp)) {
      str = MemFree (str);
    }
  }
  return str;
}


static Boolean IsParseabletRNAName (CharPtr name_string)
{
  if (StringHasNoText(name_string)) 
  {
    return TRUE;
  }
  else if (StringNICmp (name_string, "trna-", 5) != 0)
  {
    return FALSE;
  }
  else if (StringLen (name_string) != 8)
  {
    return FALSE;
  }
  else if (ParseTRnaString (name_string, NULL, NULL, TRUE) == 0)
  {
    return FALSE;
  }
  else
  {
    return TRUE;
  }
}


static Boolean SetRNAProductString (SeqFeatPtr sfp, StringConstraintPtr scp, CharPtr new_val, Uint2 existing_text)
{
  RnaRefPtr  rrp;
  Boolean rval = FALSE;
  ValNode vn;
  CharPtr cp, tmp;
  tRNAPtr trp;
  Boolean justTrnaText = FALSE;
  Uint1   codon [6];

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA || sfp->data.value.ptrvalue == NULL) {
    return FALSE;
  }

  rrp = sfp->data.value.ptrvalue;
  if (rrp->ext.choice == 0 
      || (rrp->ext.choice == 1 && StringHasNoText (rrp->ext.value.ptrvalue))
      || (rrp->ext.choice == 1 
          && (StringCmp (rrp->ext.value.ptrvalue, "ncRNA") == 0 
              || StringCmp (rrp->ext.value.ptrvalue, "tmRNA") == 0
              || StringCmp (rrp->ext.value.ptrvalue, "misc_RNA") == 0))) {
    vn.choice = FeatQualChoice_legal_qual;
    vn.data.intvalue = Feat_qual_legal_product;

    rval = SetStringInGBQualList (&(sfp->qual), &vn, scp, new_val, existing_text);
  }

  if (!rval) {
    if ((rrp->ext.choice == 0 || (rrp->ext.choice == 1 && StringHasNoText (rrp->ext.value.ptrvalue)))
        && (scp == NULL || scp->match_text == NULL)) {
      rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
      rrp->ext.value.ptrvalue = StringSave (new_val);
      rrp->ext.choice = 1;
      rval = TRUE;
    } else if (rrp->ext.choice == 1 
                && StringCmp (rrp->ext.value.ptrvalue, "ncRNA") != 0
                && StringCmp (rrp->ext.value.ptrvalue, "tmRNA") != 0
                && StringCmp (rrp->ext.value.ptrvalue, "misc_RNA") != 0
                && DoesStringMatchConstraint (rrp->ext.value.ptrvalue, scp)) {
      cp = rrp->ext.value.ptrvalue;
      rval = SetStringValue (&cp, new_val, existing_text);
      rrp->ext.value.ptrvalue = cp;
      rval = TRUE;
    } else if (rrp->ext.choice == 2) {
      tmp = GetRNAProductString (sfp, NULL);
      if (DoesStringMatchConstraint (tmp, scp)
          && SetStringValue (&tmp, new_val, existing_text)) {
        trp = (tRNAPtr) rrp->ext.value.ptrvalue;
        if (trp == NULL) {
          trp = MemNew (sizeof (tRNA));
          trp->aatype = 0;
          MemSet (trp->codon, 255, sizeof (trp->codon));
          trp->anticodon = NULL;
          rrp->ext.value.ptrvalue = trp;
        }

        if (!IsParseabletRNAName(tmp))
        {
          if (trp->anticodon == NULL
              && trp->codon[0] == 255
              && trp->codon[1] == 255
              && trp->codon[2] == 255
              && trp->codon[3] == 255
              && trp->codon[4] == 255
              && trp->codon[5] == 255)
          {
            trp = MemFree (trp);
            rrp->ext.choice = 1;
            rrp->ext.value.ptrvalue = tmp;
            tmp = NULL;
            rval = TRUE;
          }
          else
          {
            vn.choice = FeatQualChoice_legal_qual;
            vn.data.intvalue = Feat_qual_legal_product;
            if (SetStringInGBQualList (&(sfp->qual), &vn, scp, new_val, existing_text)) {
              trp->aa = 0;
              rval = TRUE;
            }
          }
        }
        else
        {
          trp->aa = ParseTRnaString (tmp, &justTrnaText, codon, TRUE);
          trp->aatype = 2;
          rval = TRUE;
        }
        tmp = MemFree (tmp);
      }
    }
  }
  return rval;
}


static Boolean RemoveRNAProductString (SeqFeatPtr sfp, StringConstraintPtr scp)
{
  RnaRefPtr  rrp;
  Boolean    rval = FALSE;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA || sfp->data.value.ptrvalue == NULL) {
    return FALSE;
  }

  rrp = sfp->data.value.ptrvalue;
  if (rrp->ext.choice == 0 
      || (rrp->ext.choice == 1 && StringHasNoText (rrp->ext.value.ptrvalue))
      || (rrp->ext.choice == 1 
          && (StringCmp (rrp->ext.value.ptrvalue, "ncRNA") == 0 
              || StringCmp (rrp->ext.value.ptrvalue, "tmRNA") == 0
              || StringCmp (rrp->ext.value.ptrvalue, "misc_RNA") == 0))) {
    rval = RemoveGBQualMatch (&(sfp->qual), "product", scp);
  }

  if (!rval 
      && rrp->ext.choice == 1 && !StringHasNoText (rrp->ext.value.ptrvalue)
      && StringCmp (rrp->ext.value.ptrvalue, "ncRNA") != 0
      && StringCmp (rrp->ext.value.ptrvalue, "tmRNA") != 0
      && StringCmp (rrp->ext.value.ptrvalue, "misc_RNA") != 0
      && DoesStringMatchConstraint(rrp->ext.value.ptrvalue, scp)) {
    rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
    rrp->ext.choice = 0;
    rval = TRUE;
  }
  return rval;
}


static SeqFeatPtr GetProtFeature (BioseqPtr protbsp)
{
  SeqMgrFeatContext fcontext;
  SeqAnnotPtr sap;
  SeqFeatPtr prot_sfp;
  ProtRefPtr prp;

  if (protbsp == NULL) return NULL;

  prot_sfp = SeqMgrGetNextFeature (protbsp, NULL, 0, FEATDEF_PROT, &fcontext);
  if (prot_sfp == NULL) {
    sap = protbsp->annot;
    while (sap != NULL && prot_sfp == NULL) {
      if (sap->type == 1) {
        prot_sfp = sap->data;
        while (prot_sfp != NULL
               && (prot_sfp->data.choice != SEQFEAT_PROT
                   || (prp = prot_sfp->data.value.ptrvalue) == NULL
                   || prp->processed != 0)) {
          prot_sfp = prot_sfp->next;
        }
      }
      sap = sap->next;
    }
  }
  return prot_sfp;
}


static ProtRefPtr GetProtRefForFeature (SeqFeatPtr sfp)
{
  BioseqPtr  protbsp;
  SeqFeatPtr protsfp;
  ProtRefPtr prp = NULL;
  SeqFeatXrefPtr xref;

  if (sfp == NULL) return NULL;

  if (sfp->data.choice == SEQFEAT_PROT) {
    prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  } else if (sfp->data.choice == SEQFEAT_CDREGION) {
    xref = sfp->xref;
    while (xref != NULL && xref->data.choice != SEQFEAT_PROT) {
      xref = xref->next;
    }
    if (xref != NULL) {
      prp = xref->data.value.ptrvalue;
    }
    if (prp == NULL && sfp->product != NULL) {
      protbsp = BioseqFindFromSeqLoc (sfp->product);
      protsfp = GetProtFeature (protbsp);    
      if (protsfp != NULL) {
        prp = protsfp->data.value.ptrvalue;
      }
    }
  }
  return prp;
}


NLM_EXTERN CharPtr GetQualFromFeature (SeqFeatPtr sfp, FeatureFieldPtr field, StringConstraintPtr scp)
{
  CharPtr   str = NULL;
  GeneRefPtr grp = NULL;
  ProtRefPtr prp = NULL;
  Int4      gbqual;

  if (sfp == NULL || field == NULL || field->field == NULL)
  {
    return NULL;
  }
  if (field->type != Feature_type_any && sfp->idx.subtype != GetFeatdefFromFeatureType (field->type))
  {
    return NULL;
  }

  // for gene fields
  if (sfp->idx.subtype == FEATDEF_GENE) {
    grp = sfp->data.value.ptrvalue;
  } else {
    grp = SeqMgrGetGeneXref (sfp);
  }

  // for protein fields
  prp = GetProtRefForFeature (sfp);

  /* fields common to all features */
  /* note, also known as comment */
  if ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_note)
      || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("note", field->field->data.ptrvalue)))
  {
    if (!StringHasNoText (sfp->comment) && DoesStringMatchConstraint(sfp->comment, scp))
    {
      str = StringSave (sfp->comment);
    }
  }
  /* db-xref */
  if (str == NULL
      && ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_db_xref)
          || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("db_xref", field->field->data.ptrvalue))))
  {
    str = GetDbxrefString (sfp, scp);
  }
  /* exception */
  if (str == NULL 
      && ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_exception)
          || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("exception", field->field->data.ptrvalue))))
  {
    if (!StringHasNoText (sfp->except_text) && DoesStringMatchConstraint(sfp->except_text, scp))
    {
      str = StringSave (sfp->except_text);
    }
  }
  /* evidence */
  if (str == NULL
      && ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_evidence)
          || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("evidence", field->field->data.ptrvalue))))
  {
    if (sfp->exp_ev == 1)
    {
      str = StringSave ("experimental");
    }
    else if (sfp->exp_ev == 2)
    {
      str = StringSave ("non-experimental");
    }
    if (!DoesStringMatchConstraint(str, scp)) {
      str = MemFree (str);
    }
  }

  /* fields common to some features */
  /* product */
  if (str == NULL
      && ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_product)
          || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("product", field->field->data.ptrvalue))))
  {
    if (prp != NULL) {
      str = GetFirstValNodeStringMatch (prp->name, scp);
    } else if (sfp->data.choice == SEQFEAT_RNA) {
      str = GetRNAProductString (sfp, scp);
    }
  }

  /* Gene fields */
  /* locus */
  if (str == NULL 
       && ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_gene)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("locus", field->field->data.ptrvalue)))
       && grp != NULL)
  {
    if (!StringHasNoText (grp->locus) && DoesStringMatchConstraint(grp->locus, scp))
    {
      str = StringSave (grp->locus);
    }
  }
  /* description */
  if (str == NULL 
       && ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_gene_description)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("description", field->field->data.ptrvalue)))
       && grp != NULL)
  {
    if (!StringHasNoText (grp->desc) && DoesStringMatchConstraint(grp->desc, scp))
    {
      str = StringSave (grp->desc);
    }
  }
  /* maploc */
  if (str == NULL 
       && ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_map)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("map", field->field->data.ptrvalue)))
       && grp != NULL)
  {
    if (!StringHasNoText (grp->maploc) && DoesStringMatchConstraint(grp->maploc, scp))
    {
      str = StringSave (grp->maploc);
    }
  }
  /* allele */
  if (str == NULL 
       && ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_allele)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("allele", field->field->data.ptrvalue)))
       && grp != NULL)
  {
    if (!StringHasNoText (grp->allele) && DoesStringMatchConstraint(grp->allele, scp))
    {
      str = StringSave (grp->allele);
    }
  }
  /* locus_tag */
  if (str == NULL 
       && ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_locus_tag)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("locus_tag", field->field->data.ptrvalue)))
       && grp != NULL)
  {
    if (!StringHasNoText (grp->locus_tag) && DoesStringMatchConstraint(grp->locus_tag, scp))
    {
      str = StringSave (grp->locus_tag);
    }
  }
  /* synonym */
  if (str == NULL 
       && ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_synonym)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("synonym", field->field->data.ptrvalue)))
       && grp != NULL)
  {
    str = GetFirstValNodeStringMatch (grp->syn, scp);
  }


  /* protein fields */
  /* note - product handled above */
  /* description */
  if (str == NULL 
       && ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_description)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("description", field->field->data.ptrvalue)))
       && prp != NULL)
  {
    if (!StringHasNoText (prp->desc) && DoesStringMatchConstraint(prp->desc, scp)) {
      str = StringSave (prp->desc);
    }
  }
  /* ec_number */
  if (str == NULL 
       && ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_ec_number)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("ec_number", field->field->data.ptrvalue)))
       && prp != NULL)
  {
    str = GetFirstValNodeStringMatch (prp->ec, scp);
  }
  /* activity */
  if (str == NULL 
       && ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_activity)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("activity", field->field->data.ptrvalue)))
       && prp != NULL)
  {
    str = GetFirstValNodeStringMatch (prp->activity, scp);
  }
  

  /* actual GenBank qualifiers */
  if (str == NULL)
  {
    if (field->field->choice == FeatQualChoice_legal_qual) 
    {
      gbqual = GetGBQualFromFeatQual (field->field->data.intvalue);
      if (gbqual > -1) {
        str = GetFirstGBQualMatch (sfp->qual, ParFlat_GBQual_names [gbqual].name, scp);
      } else {
        /* need to do something with non-qualifier qualifiers */
      }
    } else {
      str = GetFirstGBQualMatchConstraintName (sfp->qual, field->field->data.ptrvalue, scp);
    }
  }
  return str;
}


static Boolean RemoveQualFromFeature (SeqFeatPtr sfp, FeatureFieldPtr field, StringConstraintPtr scp)
{
  Boolean rval = FALSE;
  GeneRefPtr grp = NULL;
  ProtRefPtr prp = NULL;
  RnaRefPtr  rrp;
  tRNAPtr trp;
  Int4      gbqual;

  if (sfp == NULL || field == NULL || field->field == NULL)
  {
    return FALSE;
  }
  if (field->type != Feature_type_any && sfp->idx.subtype != GetFeatdefFromFeatureType (field->type))
  {
    return FALSE;
  }

  // for gene fields
  if (sfp->idx.subtype == FEATDEF_GENE) {
    grp = sfp->data.value.ptrvalue;
  } else {
    grp = SeqMgrGetGeneXref (sfp);
  }

  // for protein fields
  prp = GetProtRefForFeature (sfp);

  // for RNA fields
  if (sfp->data.choice == SEQFEAT_RNA) {
    rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  }

  /* fields common to all features */
  /* note, also known as comment */
  if ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_note)
      || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("note", field->field->data.ptrvalue)))
  {
    if (!StringHasNoText (sfp->comment) && DoesStringMatchConstraint (sfp->comment, scp))
    {
      sfp->comment = MemFree (sfp->comment);
      rval = TRUE;
    }
  }
  /* db-xref */
  if ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_db_xref)
      || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("db_xref", field->field->data.ptrvalue)))
  {
    rval = RemoveDbxrefString (sfp, scp);
  }
  /* exception */
  if ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_exception)
          || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("exception", field->field->data.ptrvalue)))
  {
    if (!StringHasNoText (sfp->except_text) && DoesStringMatchConstraint (sfp->except_text, scp))
    {
      sfp->except_text = MemFree (sfp->except_text);
      rval = TRUE;
    }
  }
  /* evidence */
  if ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_evidence)
          || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("evidence", field->field->data.ptrvalue)))
  {
    if ((sfp->exp_ev == 1 && DoesStringMatchConstraint("experimental", scp))
        || (sfp->exp_ev == 2 && DoesStringMatchConstraint("non-experimental", scp))) {
      sfp->exp_ev = 0;
      rval = TRUE;
    }
  }

  /* fields common to some features */
  /* product */
  if ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_product)
          || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("product", field->field->data.ptrvalue)))
  {
    if (prp != NULL) {
      rval = RemoveValNodeStringMatch (&(prp->name), scp);
    } else if (sfp->data.choice == SEQFEAT_RNA) {
      rval = RemoveRNAProductString (sfp, scp);
    }
  }

  /* Gene fields */
  /* locus */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_gene)
       || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("locus", field->field->data.ptrvalue)))
      && grp != NULL)
  {
    if (!StringHasNoText (grp->locus) && DoesStringMatchConstraint (grp->locus, scp)) {
      grp->locus = MemFree (grp->locus);
      rval = TRUE;
    }
  }
  /* description */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_gene_description)
       || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("description", field->field->data.ptrvalue)))
      && grp != NULL)
  {
    if (!StringHasNoText (grp->desc) && DoesStringMatchConstraint(grp->desc, scp))
    {
      grp->desc = MemFree (grp->desc);
      rval = TRUE;
    }
  }
  /* maploc */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_map)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("map", field->field->data.ptrvalue)))
       && grp != NULL)
  {
    if (!StringHasNoText (grp->maploc) && DoesStringMatchConstraint(grp->maploc, scp))
    {
      grp->maploc = MemFree (grp->maploc);
      rval = TRUE;
    }
  }
  /* allele */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_allele)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("allele", field->field->data.ptrvalue)))
      && grp != NULL)
  {
    if (!StringHasNoText (grp->allele) && DoesStringMatchConstraint(grp->allele, scp))
    {
      grp->allele = MemFree (grp->allele);
      rval = TRUE;
    }
  }
  /* locus_tag */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_locus_tag)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("locus_tag", field->field->data.ptrvalue)))
       && grp != NULL)
  {
    if (!StringHasNoText (grp->locus_tag) && DoesStringMatchConstraint(grp->locus_tag, scp))
    {
      grp->locus_tag = MemFree (grp->locus_tag);
      rval = TRUE;
    }
  }
  /* synonym */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_synonym)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("synonym", field->field->data.ptrvalue)))
       && grp != NULL)
  {
    rval = RemoveValNodeStringMatch (&(grp->syn), scp);
  }

  /* protein fields */
  /* note - product handled above */
  /* description */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_description)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("description", field->field->data.ptrvalue)))
       && prp != NULL)
  {
    if (!StringHasNoText (prp->desc) && DoesStringMatchConstraint(prp->desc, scp)) {
      prp->desc = MemFree (prp->desc);
      rval = TRUE;
    }
  }
  /* ec_number */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_ec_number)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("ec_number", field->field->data.ptrvalue)))
       && prp != NULL)
  {
    rval = RemoveValNodeStringMatch (&(prp->ec), scp);
  }
  /* activity */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_activity)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("activity", field->field->data.ptrvalue)))
       && prp != NULL)
  {
    rval = RemoveValNodeStringMatch (&(prp->activity), scp);
  }
  
  /* RNA fields */
  /* anticodon */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_anticodon)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("anticodon", field->field->data.ptrvalue)))
       && rrp != NULL && rrp->ext.choice == 2)
  {
    trp = (tRNAPtr) rrp->ext.value.ptrvalue;
    if (trp != NULL && trp->anticodon != NULL) {
      trp->anticodon = SeqLocFree (trp->anticodon);
      rval = TRUE;
    }
  }
  /* codons recognized */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_anticodon)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("anticodon", field->field->data.ptrvalue)))
       && rrp != NULL && rrp->ext.choice == 2)
  {
    trp = (tRNAPtr) rrp->ext.value.ptrvalue;
    if (trp != NULL && (trp->codon[0] != 255 || trp->codon[1] != 255 || trp->codon[2] != 255
                        || trp->codon[3] != 255 || trp->codon[4] != 255 || trp->codon[5] != 255)) {
      trp->codon [0] = 255;
      trp->codon [1] = 255;
      trp->codon [2] = 255;
      trp->codon [3] = 255;
      trp->codon [4] = 255;
      trp->codon [5] = 255;
      rval = TRUE;
    }
  }

  if (!rval) {
    /* actual GenBank qualifiers */
    if (field->field->choice == FeatQualChoice_legal_qual) 
    {
      gbqual = GetGBQualFromFeatQual (field->field->data.intvalue);
      if (gbqual > -1) {
        rval = RemoveGBQualMatch (&(sfp->qual), ParFlat_GBQual_names [gbqual].name, scp);
      } else {
        /* need to do something with non-qualifier qualifiers */
      }
    } else {
      rval = RemoveGBQualMatchConstraintName (&(sfp->qual), field->field->data.ptrvalue, scp);
    }
  }

  return rval;
}


static Boolean ChooseBestFrame (SeqFeatPtr sfp)
{
  CdRegionPtr  crp;
  Uint1        new_frame = 0, i, orig_frame;
  ByteStorePtr bs;
  Int4         lens [3];
  Int4         max;
  Boolean      retval = TRUE;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_CDREGION) return FALSE;
  
  crp = sfp->data.value.ptrvalue;
  if (crp == NULL) return FALSE;
  orig_frame = crp->frame;

  max = 0;
  for (i = 1; i <= 3; i++) {
    crp->frame = i;
    bs = ProteinFromCdRegionEx (sfp, FALSE, FALSE);
    lens[i - 1] = BSLen (bs);
    BSFree (bs);
    if (lens[i - 1] > max) {
      max = lens[i - 1];
      new_frame = i;
    }
  }
  for (i = 1; i <= 3; i++) {
    if (lens [i - 1] == max && i != new_frame) {
      retval = FALSE;
    }
  }
  if (retval) {
    crp->frame = new_frame;
  } else {
    crp->frame = orig_frame;
  }
  return retval;
}


static Boolean SetQualOnFeature (SeqFeatPtr sfp, FeatureFieldPtr field, StringConstraintPtr scp, CharPtr value, Uint2 existing_text)
{
  Boolean rval = FALSE;
  GeneRefPtr grp = NULL;
  ProtRefPtr prp = NULL;
  CharPtr    tmp;
  CdRegionPtr crp;

  if (sfp == NULL || field == NULL || field->field == NULL)
  {
    return FALSE;
  }
  if (field->type != Feature_type_any && sfp->idx.subtype != GetFeatdefFromFeatureType (field->type))
  {
    return FALSE;
  }

  // for gene fields
  if (sfp->idx.subtype == FEATDEF_GENE) {
    grp = sfp->data.value.ptrvalue;
  } else {
    grp = SeqMgrGetGeneXref (sfp);
  }

  // for protein fields
  prp = GetProtRefForFeature (sfp);

  /* fields common to all features */
  /* note, also known as comment */
  if ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_note)
      || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("note", field->field->data.ptrvalue)))
  {
    if (DoesStringMatchConstraint(sfp->comment, scp))
    {
      rval = SetStringValue ( &(sfp->comment), value, existing_text);
    }
  }
  /* db-xref */
  if ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_db_xref)
          || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("db_xref", field->field->data.ptrvalue)))
  {
    rval = SetDbxrefString (sfp, scp, value, existing_text);
  }
  /* exception */
  if ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_exception)
          || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("exception", field->field->data.ptrvalue)))
  {
    if (DoesStringMatchConstraint(sfp->except_text, scp))
    {
      rval = SetStringValue ( &(sfp->except_text), value, existing_text);
    }
  }
  /* evidence */
  if ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_evidence)
          || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("evidence", field->field->data.ptrvalue)))
  {
    tmp = NULL;
    if (sfp->exp_ev == 1)
    {
      tmp = StringSave ("experimental");
    }
    else if (sfp->exp_ev == 2)
    {
      tmp = StringSave ("non-experimental");
    }
    if (DoesStringMatchConstraint(tmp, scp)) {
      rval = SetStringValue (&tmp, value, existing_text);
      if (rval) {
        rval = FALSE;
        if (StringICmp (tmp, "experimental") == 0) {
          sfp->exp_ev = 1;
          rval = TRUE;
        } else if (StringICmp (tmp, "non-experimental") == 0) {
          sfp->exp_ev = 2;
          rval = TRUE;
        } else if (StringHasNoText (tmp)) {
          sfp->exp_ev = 0;
          rval = TRUE;
        }
      }
    }
    tmp = MemFree (tmp);
  }
  

  /* fields common to some features */
  /* product */
  if ((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_product)
          || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("product", field->field->data.ptrvalue)))
  {
    if (prp != NULL) {
      rval = SetStringsInValNodeStringList (&(prp->name), scp, value, existing_text);
    } else if (sfp->data.choice == SEQFEAT_RNA) {
      rval = SetRNAProductString (sfp, scp, value, existing_text);
    }
  }

  /* Gene fields */
  /* locus */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_gene)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("locus", field->field->data.ptrvalue)))
       && grp != NULL)
  {
    if (DoesStringMatchConstraint(grp->locus, scp))
    {
      rval = SetStringValue (&(grp->locus), value, existing_text);
    }
  }
  /* description */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_gene_description)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("description", field->field->data.ptrvalue)))
       && grp != NULL)
  {
    if (DoesStringMatchConstraint(grp->desc, scp))
    {
      rval = SetStringValue (&(grp->desc), value, existing_text);
    }
  }
  /* maploc */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_map)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("map", field->field->data.ptrvalue)))
       && grp != NULL)
  {
    if (DoesStringMatchConstraint(grp->maploc, scp))
    {
      rval = SetStringValue (&(grp->maploc), value, existing_text);
    }
  }
  /* allele */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_allele)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("allele", field->field->data.ptrvalue)))
       && grp != NULL)
  {
    if (DoesStringMatchConstraint(grp->allele, scp))
    {
      rval = SetStringValue (&(grp->allele), value, existing_text);
    }
  }
  /* locus_tag */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_locus_tag)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("locus_tag", field->field->data.ptrvalue)))
       && grp != NULL)
  {
    if (DoesStringMatchConstraint(grp->locus_tag, scp))
    {
      rval = SetStringValue (&(grp->locus_tag), value, existing_text);
    }
  }
  /* synonym */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_synonym)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("synonym", field->field->data.ptrvalue)))
       && grp != NULL)
  {
    rval = SetStringsInValNodeStringList (&(grp->syn), scp, value, existing_text);
  }


  /* protein fields */
  /* note - product handled above */
  /* description */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_description)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("description", field->field->data.ptrvalue)))
       && prp != NULL)
  {
    if (DoesStringMatchConstraint(prp->desc, scp)) {
      rval = SetStringValue (&(prp->desc), value, existing_text);
    }
  }
  /* ec_number */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_ec_number)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("ec_number", field->field->data.ptrvalue)))
       && prp != NULL)
  {
    rval = SetStringsInValNodeStringList (&(prp->ec), scp, value, existing_text);
  }
  /* activity */
  if (((field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_activity)
           || (field->field->choice == FeatQualChoice_illegal_qual && DoesStringMatchConstraint ("activity", field->field->data.ptrvalue)))
       && prp != NULL)
  {
    rval = SetStringsInValNodeStringList (&(prp->activity), scp, value, existing_text);
  }
 
  if (field->field->choice == FeatQualChoice_legal_qual && field->field->data.intvalue == Feat_qual_legal_codon_start
      && sfp->data.choice == SEQFEAT_CDREGION) 
  {
    crp = (CdRegionPtr) sfp->data.value.ptrvalue;
    if (StringICmp (value, "best") == 0)
    {
      rval = ChooseBestFrame (sfp);
    }
    else if (StringCmp (value, "1") == 0) 
    {
      crp->frame = 1;
      rval = TRUE;
    }
    else if (StringCmp (value, "2") == 0) 
    {
      crp->frame = 2;
      rval = TRUE;
    }
    else if (StringCmp (value, "3") == 0)
    {
      crp->frame = 3;
      rval = TRUE;
    } 
  } 

  /* actual GenBank qualifiers */
  if (!rval)
  {
    rval = SetStringInGBQualList (&(sfp->qual), field->field, scp, value, existing_text);
  }
  return rval;
}


NLM_EXTERN CharPtr GetSourceQualFromBioSource (BioSourcePtr biop, SourceQualChoicePtr scp, StringConstraintPtr constraint)
{
  CharPtr str = NULL;
  SubSourcePtr ssp;
  OrgModPtr mod;
  Int4 orgmod_subtype = -1, subsrc_subtype = -1;

  if (biop == NULL || scp == NULL) return NULL;

  switch (scp->choice) 
  {
    case SourceQualChoice_textqual:
      if (scp->data.intvalue == Source_qual_taxname) {
        if (biop->org != NULL && !StringHasNoText (biop->org->taxname)
            && DoesStringMatchConstraint (biop->org->taxname, constraint)) {
          str = StringSave (biop->org->taxname);
        }
      } else if (scp->data.intvalue == Source_qual_common_name) {
        if (biop->org != NULL && !StringHasNoText (biop->org->common)
            && DoesStringMatchConstraint (biop->org->common, constraint)) {
          str = StringSave (biop->org->common);
        }
      } else if (scp->data.intvalue == Source_qual_lineage) {
        if (biop->org != NULL && biop->org->orgname != NULL && !StringHasNoText (biop->org->orgname->lineage)
            && DoesStringMatchConstraint (biop->org->orgname->lineage, constraint)) {
          str = StringSave (biop->org->orgname->lineage);
        }
      } else if (scp->data.intvalue == Source_qual_division) {
        if (biop->org != NULL && biop->org->orgname != NULL  && !StringHasNoText (biop->org->orgname->div)
            && DoesStringMatchConstraint (biop->org->orgname->div, constraint)) {
          str = StringSave (biop->org->orgname->div);
        }
      } else {
        orgmod_subtype = GetOrgModQualFromSrcQual (scp->data.intvalue);
        if (orgmod_subtype == -1) {
          subsrc_subtype = GetSubSrcQualFromSrcQual (scp->data.intvalue);
          for (ssp = biop->subtype; ssp != NULL && str == NULL; ssp = ssp->next) {
            if (ssp->subtype == subsrc_subtype) {
              if (StringHasNoText (ssp->name)) {
                if (IsNonTextSourceQual (scp->data.intvalue)
                    && DoesStringMatchConstraint ("TRUE", constraint)) {
                  str = StringSave ("TRUE");
                }
              } else {
                if (DoesStringMatchConstraint (ssp->name, constraint)) {
                  str = StringSave (ssp->name);
                }
              }
            }
          }
        } else {
          if (biop->org != NULL && biop->org->orgname != NULL) {
            for (mod = biop->org->orgname->mod; mod != NULL && str == NULL; mod = mod->next) {
              if (mod->subtype == orgmod_subtype) {
                if (StringHasNoText (mod->subname)) {
                  if (IsNonTextSourceQual (scp->data.intvalue)
                      && DoesStringMatchConstraint ("TRUE", constraint)) {
                    str = StringSave ("TRUE");
                  }
                } else {
                  if (DoesStringMatchConstraint (mod->subname, constraint)) {
                    str = StringSave (mod->subname);
                  }
                }
              }
            }
          }
        }
      }
      break;
    case SourceQualChoice_location:
      str = LocNameFromGenome (biop->genome);
      if (DoesStringMatchConstraint (str, constraint)) {
        str = StringSave (str);
      } else {
        str = NULL;
      }
      break;
    case SourceQualChoice_origin:
      str = OriginNameFromOrigin (biop->origin);
      if (DoesStringMatchConstraint (str, constraint)) {
        str = StringSave (str);
      } else {
        str = NULL;
      }
      break;
  }
  return str;
}


static Boolean RemoveSourceQualFromBioSource (BioSourcePtr biop, SourceQualChoicePtr scp, StringConstraintPtr constraint)
{
  SubSourcePtr ssp, ssp_prev = NULL, ssp_next;
  OrgModPtr mod, mod_prev = NULL, mod_next;
  Int4 orgmod_subtype = -1, subsrc_subtype = -1;
  CharPtr str;
  Boolean rval = FALSE;

  if (biop == NULL || scp == NULL) return FALSE;

  switch (scp->choice) 
  {
    case SourceQualChoice_textqual:
      if (scp->data.intvalue == Source_qual_taxname) {
        if (biop->org != NULL && !StringHasNoText (biop->org->taxname)
            && DoesStringMatchConstraint (biop->org->taxname, constraint)) {
          biop->org->taxname = MemFree (biop->org->taxname);
          rval = TRUE;
        }
      } else if (scp->data.intvalue == Source_qual_common_name) {
        if (biop->org != NULL && !StringHasNoText (biop->org->common)
            && DoesStringMatchConstraint (biop->org->common, constraint)) {
          biop->org->common = MemFree (biop->org->common);
          rval = TRUE;
        }
      } else if (scp->data.intvalue == Source_qual_lineage) {
        if (biop->org != NULL && biop->org->orgname != NULL && !StringHasNoText (biop->org->orgname->lineage)
            && DoesStringMatchConstraint (biop->org->orgname->lineage, constraint)) {
          biop->org->orgname->lineage = MemFree (biop->org->orgname->lineage);
          rval = TRUE;
        }
      } else if (scp->data.intvalue == Source_qual_division) {
        if (biop->org != NULL && biop->org->orgname != NULL  && !StringHasNoText (biop->org->orgname->div)
            && DoesStringMatchConstraint (biop->org->orgname->div, constraint)) {
          biop->org->orgname->div = MemFree (biop->org->orgname->div);
          rval = TRUE;
        }
      } else {
        orgmod_subtype = GetOrgModQualFromSrcQual (scp->data.intvalue);
        if (orgmod_subtype == -1) {
          subsrc_subtype = GetSubSrcQualFromSrcQual (scp->data.intvalue);
          ssp = biop->subtype;
          while (ssp != NULL) {
            ssp_next = ssp->next;
            if (ssp->subtype == subsrc_subtype 
                && DoesStringMatchConstraint (ssp->name, constraint)) {
              if (ssp_prev == NULL) {
                biop->subtype = ssp->next;
              } else {
                ssp_prev->next = ssp->next;
              }
              ssp->next = NULL;
              ssp = SubSourceFree (ssp);
              rval = TRUE;
            } else {
              ssp_prev = ssp;
            }
            ssp = ssp_next;
          }
        } else {
          if (biop->org != NULL && biop->org->orgname != NULL) {
            mod = biop->org->orgname->mod;
            while (mod != NULL) {
              mod_next = mod->next;
              if (mod->subtype == orgmod_subtype
                  && DoesStringMatchConstraint (mod->subname, constraint)) {
                if (mod_prev == NULL) {
                  biop->org->orgname->mod = mod->next;
                } else {
                  mod_prev->next = mod->next;
                }
                mod->next = NULL;
                mod = OrgModFree (mod);
                rval = TRUE;
              } else {
                mod_prev = mod;
              }
              mod = mod_next;
            }
          }
        }
      }
      break;
    case SourceQualChoice_location:
      str = LocNameFromGenome (biop->genome);
      if (DoesStringMatchConstraint (str, constraint)) {
        if (scp->data.intvalue == 0 || biop->genome == GenomeFromSrcLoc (scp->data.intvalue)) {
          biop->genome = 0;
          rval = TRUE;
        }
      }
      break;
    case SourceQualChoice_origin:
      str = OriginNameFromOrigin (biop->origin);
      if (DoesStringMatchConstraint (str, constraint)) {
        if (scp->data.intvalue == 0 || biop->origin == OriginFromSrcOrig (scp->data.intvalue)) {
          biop->origin = 0;
          rval = TRUE;
        }
      }
      break; 
  }
  return rval;
}


NLM_EXTERN Boolean SetSourceQualInBioSource (BioSourcePtr biop, SourceQualChoicePtr scp, StringConstraintPtr constraint, CharPtr value, Uint2 existing_text)
{
  SubSourcePtr ssp, ssp_prev = NULL, ssp_next;
  OrgModPtr mod, mod_prev = NULL, mod_next;
  Int4 orgmod_subtype = -1, subsrc_subtype = -1;
  CharPtr str;
  Boolean rval = FALSE, found = FALSE;

  if (biop == NULL || scp == NULL) return FALSE;

  switch (scp->choice) 
  {
    case SourceQualChoice_textqual:
      if (scp->data.intvalue == Source_qual_taxname) {
        if (biop->org != NULL
            && DoesStringMatchConstraint (biop->org->taxname, constraint)) {
          rval = SetStringValue (&(biop->org->taxname), value, existing_text);
        }
      } else if (scp->data.intvalue == Source_qual_common_name) {
        if (biop->org != NULL
            && DoesStringMatchConstraint (biop->org->common, constraint)) {
          rval = SetStringValue (&(biop->org->common), value, existing_text);
        }
      } else if (scp->data.intvalue == Source_qual_lineage) {
        if (biop->org != NULL && biop->org->orgname != NULL 
            && DoesStringMatchConstraint (biop->org->orgname->lineage, constraint)) {
          rval = SetStringValue (&(biop->org->orgname->lineage), value, existing_text);
        }
      } else if (scp->data.intvalue == Source_qual_division) {
        if (biop->org != NULL && biop->org->orgname != NULL
            && DoesStringMatchConstraint (biop->org->orgname->div, constraint)) {
          rval = SetStringValue (&(biop->org->orgname->div), value, existing_text);
        }
      } else {
        orgmod_subtype = GetOrgModQualFromSrcQual (scp->data.intvalue);
        if (orgmod_subtype == -1) {
          subsrc_subtype = GetSubSrcQualFromSrcQual (scp->data.intvalue);
          if (subsrc_subtype > -1) {
            ssp = biop->subtype;
            while (ssp != NULL) {
              ssp_next = ssp->next;
              if (ssp->subtype == subsrc_subtype
                  && DoesStringMatchConstraint (ssp->name, constraint)) {
                rval = SetStringValue (&(ssp->name), value, existing_text);
                found = TRUE;
                if (rval && StringHasNoText (ssp->name) && !IsNonTextSourceQual(scp->data.intvalue)) {
                  if (ssp_prev == NULL) {
                    biop->subtype = ssp->next;
                  } else {
                    ssp_prev->next = ssp->next;
                  }
                  ssp->next = NULL;
                  ssp = SubSourceFree (ssp);
                } else {
                  ssp_prev = ssp;
                }
              } else {
                ssp_prev = ssp;
              }
              ssp = ssp_next;
            }
            if (!found && IsStringConstraintEmpty (constraint)) {
              ssp = SubSourceNew ();
              ssp->subtype = subsrc_subtype;
              rval = SetStringValue (&(ssp->name), value, existing_text);
              if (ssp_prev == NULL) {
                biop->subtype = ssp;
              } else {
                ssp_prev->next = ssp;
              }
            }
          }
        } else {
          if (biop->org != NULL && biop->org->orgname != NULL) {
            mod = biop->org->orgname->mod;
            while (mod != NULL) {
              mod_next = mod->next;
              if (mod->subtype == orgmod_subtype
                  && DoesStringMatchConstraint (mod->subname, constraint)) {
                rval = SetStringValue (&(mod->subname), value, existing_text);
                found = TRUE;
                if (rval && StringHasNoText (mod->subname) && !IsNonTextSourceQual(scp->data.intvalue)) {
                  if (mod_prev == NULL) {
                    biop->org->orgname->mod = mod->next;
                  } else {
                    mod_prev->next = mod->next;
                  }
                  mod->next = NULL;
                  mod = OrgModFree (mod);
                } else {
                  mod_prev = mod;
                }
              } else {
                mod_prev = mod;
              }
              mod = mod_next;
            }
          }
          if (!found && IsStringConstraintEmpty (constraint)) {
            if (biop->org == NULL) {
              biop->org = OrgRefNew();
            }
            if (biop->org->orgname == NULL) {
              biop->org->orgname = OrgNameNew();
            }
            mod = OrgModNew ();
            mod->subtype = orgmod_subtype;
            rval = SetStringValue (&(mod->subname), value, existing_text);
            if (mod_prev == NULL) {
              biop->org->orgname->mod = mod;
            } else {
              mod_prev->next = mod;
            }
          }
        }
      }
      break;
    case SourceQualChoice_location:
      str = LocNameFromGenome (biop->genome);
      if (DoesStringMatchConstraint (str, constraint)) {
        biop->genome = GenomeFromSrcLoc (scp->data.intvalue);
        rval = TRUE;
      }
      break;
    case SourceQualChoice_origin:
      str = OriginNameFromOrigin (biop->origin);
      if (DoesStringMatchConstraint (str, constraint)) {
        biop->origin = OriginFromSrcOrig(scp->data.intvalue);
        rval = TRUE;
      }
      break; 
  }
  return rval;
}


static BioseqPtr GetSequenceForObject (Uint1 choice, Pointer data)
{
  BioseqPtr bsp = NULL;
  SeqFeatPtr sfp;
  SeqDescrPtr sdp;
  ObjValNodePtr ovp;
  CGPSetPtr cgp;
  ValNodePtr vnp;

  if (data == NULL) return NULL;

  switch (choice) {
    case OBJ_BIOSEQ:
      bsp = (BioseqPtr) data;
      break;
    case OBJ_SEQFEAT:
      sfp = (SeqFeatPtr) data;
      bsp = BioseqFindFromSeqLoc (sfp->location);
      break;
    case OBJ_SEQDESC:
      sdp = (SeqDescrPtr) data;
      if (sdp->extended) {
        ovp = (ObjValNodePtr) sdp;
        if (ovp->idx.parenttype == OBJ_BIOSEQ && ovp->idx.parentptr != NULL) {
          bsp = ovp->idx.parentptr;
        }
      }
      break;
    case 0:
      cgp = (CGPSetPtr) data;
      for (vnp = cgp->cds_list; vnp != NULL && bsp == NULL; vnp = vnp->next) {
        sfp = vnp->data.ptrvalue;
        if (sfp != NULL) {
          bsp = BioseqFindFromSeqLoc (sfp->location);
        }
      }
      for (vnp = cgp->mrna_list; vnp != NULL && bsp == NULL; vnp = vnp->next) {
        sfp = vnp->data.ptrvalue;
        if (sfp != NULL) {
          bsp = BioseqFindFromSeqLoc (sfp->location);
        }
      }
      break;
      for (vnp = cgp->gene_list; vnp != NULL && bsp == NULL; vnp = vnp->next) {
        sfp = vnp->data.ptrvalue;
        if (sfp != NULL) {
          bsp = BioseqFindFromSeqLoc (sfp->location);
        }
      }
      break;
  }
  return bsp;
}


NLM_EXTERN BioSourcePtr GetBioSourceFromObject (Uint1 choice, Pointer data)
{
  BioSourcePtr biop = NULL;
  SeqDescrPtr  sdp;
  SeqFeatPtr   sfp;
  BioseqPtr    bsp = NULL;
  SeqMgrDescContext context;

  if (data == NULL) return NULL;

  switch (choice)
  {
    case OBJ_SEQDESC:
      sdp = (SeqDescrPtr) data;
      if (sdp->choice == Seq_descr_source) {
        biop = sdp->data.ptrvalue;
      }
      break;
    case OBJ_SEQFEAT:
      sfp = (SeqFeatPtr) data;
      if (sfp->data.choice == SEQFEAT_BIOSRC) {
        biop = sfp->data.value.ptrvalue;
      }
      break;
  }
  if (biop == NULL) {
    bsp = GetSequenceForObject (choice, data);
    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
    if (sdp != NULL && sdp->choice == Seq_descr_source) {
      biop = sdp->data.ptrvalue;
    }
  }
  return biop;
}


/* functions for dealing with CDS-Gene-Prot sets */
static CharPtr GetFieldValueFromCGPSet (CGPSetPtr c, Uint2 field, StringConstraintPtr scp)
{
  CharPtr str = NULL;
  ValNodePtr vnp;
  SeqFeatPtr sfp;
  GeneRefPtr grp;
  RnaRefPtr  rrp;
  ProtRefPtr prp;
  
  if (c == NULL) return NULL;
  switch (field) {
    case CDSGeneProt_field_cds_comment:
      for (vnp = c->cds_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && !StringHasNoText (sfp->comment) && DoesStringMatchConstraint(sfp->comment, scp))
        {
          str = StringSave (sfp->comment);
        }
      }
      break;
    case CDSGeneProt_field_gene_locus:
      for (vnp = c->gene_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE
            && (grp = sfp->data.value.ptrvalue) != NULL
            && !StringHasNoText (grp->locus) 
            && DoesStringMatchConstraint(grp->locus, scp))
        {
          str = StringSave (grp->locus);
        }
      }
      break;
    case CDSGeneProt_field_gene_description:
      for (vnp = c->gene_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE
            && (grp = sfp->data.value.ptrvalue) != NULL
            && !StringHasNoText (grp->desc) 
            && DoesStringMatchConstraint(grp->desc, scp))
        {
          str = StringSave (grp->desc);
        }
      }
      break;
    case CDSGeneProt_field_gene_comment:
      for (vnp = c->gene_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && !StringHasNoText (sfp->comment) && DoesStringMatchConstraint(sfp->comment, scp))
        {
          str = StringSave (sfp->comment);
        }
      }
      break;
    case CDSGeneProt_field_gene_allele:
      for (vnp = c->gene_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE
            && (grp = sfp->data.value.ptrvalue) != NULL
            && !StringHasNoText (grp->allele) 
            && DoesStringMatchConstraint(grp->allele, scp))
        {
          str = StringSave (grp->allele);
        }
      }
      break;
    case CDSGeneProt_field_gene_maploc:
      for (vnp = c->gene_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE
            && (grp = sfp->data.value.ptrvalue) != NULL
            && !StringHasNoText (grp->maploc) 
            && DoesStringMatchConstraint(grp->maploc, scp))
        {
          str = StringSave (grp->maploc);
        }
      }
      break;
    case CDSGeneProt_field_gene_locus_tag:
      for (vnp = c->gene_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE
            && (grp = sfp->data.value.ptrvalue) != NULL
            && !StringHasNoText (grp->locus_tag) 
            && DoesStringMatchConstraint(grp->locus_tag, scp))
        {
          str = StringSave (grp->locus_tag);
        }
      }
      break;
    case CDSGeneProt_field_gene_synonym:
      for (vnp = c->gene_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE
            && (grp = sfp->data.value.ptrvalue) != NULL)
        {
          str = GetFirstValNodeStringMatch (grp->syn, scp);
        }
      }
      break;
    case CDSGeneProt_field_gene_old_locus_tag:
      for (vnp = c->gene_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL) {
          str = GetFirstGBQualMatch (sfp->qual, "old-locus-tag", scp);
        }
      }
      break;
    case CDSGeneProt_field_mrna_product:
      for (vnp = c->mrna_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_RNA
            && (rrp = sfp->data.value.ptrvalue) != NULL
            && rrp->ext.choice == 1
            && !StringHasNoText (rrp->ext.value.ptrvalue) 
            && DoesStringMatchConstraint(rrp->ext.value.ptrvalue, scp))
        {
          str = StringSave (rrp->ext.value.ptrvalue);
        }
      }
      break;
    case CDSGeneProt_field_mrna_comment:
      for (vnp = c->mrna_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && !StringHasNoText (sfp->comment) && DoesStringMatchConstraint(sfp->comment, scp))
        {
          str = StringSave (sfp->comment);
        }
      }
      break;
    case CDSGeneProt_field_prot_name:
      for (vnp = c->prot_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_PROT
            && (prp = sfp->data.value.ptrvalue) != NULL)
        {
          str = GetFirstValNodeStringMatch (prp->name, scp);
        }
      }
      break;
    case CDSGeneProt_field_prot_description:
      for (vnp = c->prot_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_PROT
            && (prp = sfp->data.value.ptrvalue) != NULL
            && !StringHasNoText (prp->desc) && DoesStringMatchConstraint(prp->desc, scp)) {
          str = StringSave (prp->desc);
        }
      }
      break;
    case CDSGeneProt_field_prot_ec_number:
      for (vnp = c->prot_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_PROT
            && (prp = sfp->data.value.ptrvalue) != NULL)
        {
          str = GetFirstValNodeStringMatch (prp->ec, scp);
        }
      }
      break;
    case CDSGeneProt_field_prot_activity:
      for (vnp = c->prot_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_PROT
            && (prp = sfp->data.value.ptrvalue) != NULL)
        {
          str = GetFirstValNodeStringMatch (prp->activity, scp);
        }
      }
      break;
    case CDSGeneProt_field_prot_comment:
      for (vnp = c->prot_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->idx.subtype == FEATDEF_PROT
            && !StringHasNoText (sfp->comment) && DoesStringMatchConstraint(sfp->comment, scp))
        {
          str = StringSave (sfp->comment);
        }
      }
      break;
    case CDSGeneProt_field_mat_peptide_name:
      for (vnp = c->prot_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_mat_peptide_aa
            && (prp = sfp->data.value.ptrvalue) != NULL)
        {
          str = GetFirstValNodeStringMatch (prp->name, scp);
        }
      }
      break;
    case CDSGeneProt_field_mat_peptide_description:
      for (vnp = c->prot_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_mat_peptide_aa
            && (prp = sfp->data.value.ptrvalue) != NULL
            && !StringHasNoText (prp->desc) && DoesStringMatchConstraint(prp->desc, scp)) {
          str = StringSave (prp->desc);
        }
      }
      break;
    case CDSGeneProt_field_mat_peptide_ec_number:
      for (vnp = c->prot_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_mat_peptide_aa
            && (prp = sfp->data.value.ptrvalue) != NULL)
        {
          str = GetFirstValNodeStringMatch (prp->ec, scp);
        }
      }
      break;
    case CDSGeneProt_field_mat_peptide_activity:
      for (vnp = c->prot_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_mat_peptide_aa
            && (prp = sfp->data.value.ptrvalue) != NULL)
        {
          str = GetFirstValNodeStringMatch (prp->activity, scp);
        }
      }
      break;
    case CDSGeneProt_field_mat_peptide_comment:
      for (vnp = c->prot_list; vnp != NULL && str == NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->idx.subtype == FEATDEF_mat_peptide_aa
            && !StringHasNoText (sfp->comment) && DoesStringMatchConstraint(sfp->comment, scp))
        {
          str = StringSave (sfp->comment);
        }
      }
      break;
  }
  return str;
}


static Boolean RemoveFieldValueFromCGPSet (CGPSetPtr c, Uint2 field, StringConstraintPtr scp)
{
  Boolean    rval = FALSE;
  ValNodePtr vnp;
  SeqFeatPtr sfp;
  GeneRefPtr grp;
  RnaRefPtr  rrp;
  ProtRefPtr prp;
  
  if (c == NULL) return FALSE;
  switch (field) {
    case CDSGeneProt_field_cds_comment:
      for (vnp = c->cds_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && !StringHasNoText (sfp->comment) && DoesStringMatchConstraint(sfp->comment, scp))
        {
          sfp->comment = MemFree (sfp->comment);
          rval = TRUE;
        }
      }
      break;
    case CDSGeneProt_field_gene_locus:
      for (vnp = c->gene_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE
            && (grp = sfp->data.value.ptrvalue) != NULL
            && !StringHasNoText (grp->locus) 
            && DoesStringMatchConstraint(grp->locus, scp))
        {
          grp->locus = MemFree (grp->locus);
          rval = TRUE;
        }
      }
      break;
    case CDSGeneProt_field_gene_description:
      for (vnp = c->gene_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE
            && (grp = sfp->data.value.ptrvalue) != NULL
            && !StringHasNoText (grp->desc) 
            && DoesStringMatchConstraint(grp->desc, scp))
        {
          grp->desc = MemFree(grp->desc);
          rval = TRUE;
        }
      }
      break;
    case CDSGeneProt_field_gene_comment:
      for (vnp = c->gene_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && !StringHasNoText (sfp->comment) && DoesStringMatchConstraint(sfp->comment, scp))
        {
          sfp->comment = MemFree (sfp->comment);
          rval = TRUE;
        }
      }
      break;
    case CDSGeneProt_field_gene_allele:
      for (vnp = c->gene_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE
            && (grp = sfp->data.value.ptrvalue) != NULL
            && !StringHasNoText (grp->allele) 
            && DoesStringMatchConstraint(grp->allele, scp))
        {
          grp->allele = MemFree (grp->allele);
          rval = TRUE;
        }
      }
      break;
    case CDSGeneProt_field_gene_maploc:
      for (vnp = c->gene_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE
            && (grp = sfp->data.value.ptrvalue) != NULL
            && !StringHasNoText (grp->maploc) 
            && DoesStringMatchConstraint(grp->maploc, scp))
        {
          grp->maploc = MemFree (grp->maploc);
          rval = TRUE;
        }
      }
      break;
    case CDSGeneProt_field_gene_locus_tag:
      for (vnp = c->gene_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE
            && (grp = sfp->data.value.ptrvalue) != NULL
            && !StringHasNoText (grp->locus_tag) 
            && DoesStringMatchConstraint(grp->locus_tag, scp))
        {
          grp->locus_tag = MemFree (grp->locus_tag);
          rval = TRUE;
        }
      }
      break;
    case CDSGeneProt_field_gene_synonym:
      for (vnp = c->gene_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE
            && (grp = sfp->data.value.ptrvalue) != NULL)
        {
          rval |= RemoveValNodeStringMatch (&(grp->syn), scp);
        }
      }
      break;
    case CDSGeneProt_field_gene_old_locus_tag:
      for (vnp = c->gene_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL) {
          rval |= RemoveGBQualMatch (&(sfp->qual), "old-locus-tag", scp);
        }
      }
      break;
    case CDSGeneProt_field_mrna_product:
      for (vnp = c->mrna_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_RNA
            && (rrp = sfp->data.value.ptrvalue) != NULL
            && rrp->ext.choice == 1
            && !StringHasNoText (rrp->ext.value.ptrvalue) 
            && DoesStringMatchConstraint(rrp->ext.value.ptrvalue, scp))
        {
          rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
          rrp->ext.choice = 0;
          rval = TRUE;
        }
      }
      break;
    case CDSGeneProt_field_mrna_comment:
      for (vnp = c->mrna_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && !StringHasNoText (sfp->comment) && DoesStringMatchConstraint(sfp->comment, scp))
        {
          sfp->comment = MemFree (sfp->comment);
          rval = TRUE;
        }
      }
      break;
    case CDSGeneProt_field_prot_name:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_PROT
            && (prp = sfp->data.value.ptrvalue) != NULL)
        {
          rval |= RemoveValNodeStringMatch (&(prp->name), scp);
        }
      }
      break;
    case CDSGeneProt_field_prot_description:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_PROT
            && (prp = sfp->data.value.ptrvalue) != NULL
            && !StringHasNoText (prp->desc) && DoesStringMatchConstraint(prp->desc, scp)) {
          prp->desc = MemFree (prp->desc);
          rval = TRUE;
        }
      }
      break;
    case CDSGeneProt_field_prot_ec_number:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_PROT
            && (prp = sfp->data.value.ptrvalue) != NULL)
        {
          rval |= RemoveValNodeStringMatch (&(prp->ec), scp);
        }
      }
      break;
    case CDSGeneProt_field_prot_activity:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_PROT
            && (prp = sfp->data.value.ptrvalue) != NULL)
        {
          rval |= RemoveValNodeStringMatch (&(prp->activity), scp);
        }
      }
      break;
    case CDSGeneProt_field_prot_comment:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->idx.subtype == FEATDEF_PROT
            && !StringHasNoText (sfp->comment) && DoesStringMatchConstraint(sfp->comment, scp))
        {
          sfp->comment = MemFree (sfp->comment);
          rval = TRUE;
        }
      }
      break;
    case CDSGeneProt_field_mat_peptide_name:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_mat_peptide_aa
            && (prp = sfp->data.value.ptrvalue) != NULL)
        {
          rval |= RemoveValNodeStringMatch (&(prp->name), scp);
        }
      }
      break;
    case CDSGeneProt_field_mat_peptide_description:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_mat_peptide_aa
            && (prp = sfp->data.value.ptrvalue) != NULL
            && !StringHasNoText (prp->desc) && DoesStringMatchConstraint(prp->desc, scp)) {
          prp->desc = MemFree (prp->desc);
          rval = TRUE;
        }
      }
      break;
    case CDSGeneProt_field_mat_peptide_ec_number:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_mat_peptide_aa
            && (prp = sfp->data.value.ptrvalue) != NULL)
        {
          rval |= RemoveValNodeStringMatch (&(prp->ec), scp);
        }
      }
      break;
    case CDSGeneProt_field_mat_peptide_activity:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_mat_peptide_aa
            && (prp = sfp->data.value.ptrvalue) != NULL)
        {
          rval |= RemoveValNodeStringMatch (&(prp->activity), scp);
        }
      }
      break;
    case CDSGeneProt_field_mat_peptide_comment:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->idx.subtype == FEATDEF_mat_peptide_aa
            && !StringHasNoText (sfp->comment) && DoesStringMatchConstraint(sfp->comment, scp))
        {
          sfp->comment = MemFree (sfp->comment);
          rval = TRUE;
        }
      }
      break;
  }
  return rval;
}


static SeqFeatPtr CreateGeneForCGPSet (CGPSetPtr c)
{
  SeqFeatPtr gene = NULL, sfp = NULL;
  BioseqPtr  bsp;
  ValNodePtr vnp;

  if (c == NULL) return NULL;

  for (vnp = c->cds_list; vnp != NULL && sfp == NULL; vnp = vnp->next) {
    sfp = vnp->data.ptrvalue;
  }
  for (vnp = c->mrna_list; vnp != NULL && sfp == NULL; vnp = vnp->next) {
    sfp = vnp->data.ptrvalue;
  }
  if (sfp != NULL) {
    bsp = BioseqFindFromSeqLoc (sfp->location);
    if (bsp != NULL) {
      gene = CreateNewFeatureOnBioseq (bsp, SEQFEAT_GENE, sfp->location);
      if (gene != NULL) {
        gene->data.value.ptrvalue = GeneRefNew();
      }
    }
  }
  return gene;
}


static Boolean SetFieldValueInCGPSet (CGPSetPtr c, Uint2 field, StringConstraintPtr scp, CharPtr value, Uint2 existing_text)
{
  Boolean    rval = FALSE;
  ValNodePtr vnp;
  SeqFeatPtr sfp;
  GeneRefPtr grp;
  ProtRefPtr prp;
  
  if (c == NULL) return FALSE;
  switch (field) {
    case CDSGeneProt_field_cds_comment:
      for (vnp = c->cds_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && DoesStringMatchConstraint(sfp->comment, scp))
        {
          rval |= SetStringValue ( &(sfp->comment), value, existing_text);
        }
      }
      break;
    case CDSGeneProt_field_gene_locus:
      if (c->gene_list == NULL && scp == NULL) {
        sfp = CreateGeneForCGPSet (c);
        if (sfp != NULL) {
          ValNodeAddPointer (&(c->gene_list), 0, sfp);
        }
      }
      for (vnp = c->gene_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE
            && (grp = sfp->data.value.ptrvalue) != NULL
            && DoesStringMatchConstraint(grp->locus, scp))
        {
          rval |= SetStringValue ( &(grp->locus), value, existing_text);
        }
      }
      break;
    case CDSGeneProt_field_gene_description:
      for (vnp = c->gene_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE
            && (grp = sfp->data.value.ptrvalue) != NULL
            && DoesStringMatchConstraint(grp->desc, scp))
        {
          rval |= SetStringValue ( &(grp->desc), value, existing_text);
        }
      }
      break;
    case CDSGeneProt_field_gene_comment:
      for (vnp = c->gene_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && DoesStringMatchConstraint(sfp->comment, scp))
        {
          rval |= SetStringValue ( &(sfp->comment), value, existing_text);
        }
      }
      break;
    case CDSGeneProt_field_gene_allele:
      for (vnp = c->gene_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE
            && (grp = sfp->data.value.ptrvalue) != NULL
            && DoesStringMatchConstraint(grp->allele, scp))
        {
          rval |= SetStringValue (&(grp->allele), value, existing_text);
        }
      }
      break;
    case CDSGeneProt_field_gene_maploc:
      for (vnp = c->gene_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE
            && (grp = sfp->data.value.ptrvalue) != NULL
            && DoesStringMatchConstraint(grp->maploc, scp))
        {
          rval |= SetStringValue ( &(grp->maploc), value, existing_text);
        }
      }
      break;
    case CDSGeneProt_field_gene_locus_tag:
      for (vnp = c->gene_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE
            && (grp = sfp->data.value.ptrvalue) != NULL
            && DoesStringMatchConstraint(grp->locus_tag, scp))
        {
          rval |= SetStringValue ( &(grp->locus_tag), value, existing_text);
        }
      }
      break;
    case CDSGeneProt_field_gene_synonym:
      for (vnp = c->gene_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_GENE
            && (grp = sfp->data.value.ptrvalue) != NULL)
        {
          rval |= SetStringsInValNodeStringList (&(grp->syn), scp, value, existing_text);
        }
      }
      break;
    case CDSGeneProt_field_gene_old_locus_tag:
      for (vnp = c->gene_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL) {
          rval |= RemoveGBQualMatch (&(sfp->qual), "old-locus-tag", scp);
        }
      }
      break;
    case CDSGeneProt_field_mrna_product:
      for (vnp = c->mrna_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        rval |= SetRNAProductString (sfp, scp, value, existing_text);
      }
      break;
    case CDSGeneProt_field_mrna_comment:
      for (vnp = c->mrna_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL&& DoesStringMatchConstraint(sfp->comment, scp))
        {
          rval |= SetStringValue ( &(sfp->comment), value, existing_text);
        }
      }
      break;
    case CDSGeneProt_field_prot_name:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_PROT
            && (prp = sfp->data.value.ptrvalue) != NULL)
        {
          rval |= SetStringsInValNodeStringList (&(prp->name), scp, value, existing_text);
        }
      }
      break;
    case CDSGeneProt_field_prot_description:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_PROT
            && (prp = sfp->data.value.ptrvalue) != NULL
            && DoesStringMatchConstraint(prp->desc, scp)) {
          rval |= SetStringValue ( &(prp->desc), value, existing_text);
        }
      }
      break;
    case CDSGeneProt_field_prot_ec_number:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_PROT
            && (prp = sfp->data.value.ptrvalue) != NULL)
        {
          rval |= SetStringsInValNodeStringList (&(prp->ec), scp, value, existing_text);
        }
      }
      break;
    case CDSGeneProt_field_prot_activity:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_PROT
            && (prp = sfp->data.value.ptrvalue) != NULL)
        {
          rval |= SetStringsInValNodeStringList (&(prp->activity), scp, value, existing_text);
        }
      }
      break;
    case CDSGeneProt_field_prot_comment:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->idx.subtype == FEATDEF_PROT
            && DoesStringMatchConstraint(sfp->comment, scp))
        {
          rval |= SetStringValue ( &(sfp->comment), value, existing_text);
        }
      }
      break;
    case CDSGeneProt_field_mat_peptide_name:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_mat_peptide_aa
            && (prp = sfp->data.value.ptrvalue) != NULL)
        {
          rval |= SetStringsInValNodeStringList (&(prp->name), scp, value, existing_text);
        }
      }
      break;
    case CDSGeneProt_field_mat_peptide_description:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_mat_peptide_aa
            && (prp = sfp->data.value.ptrvalue) != NULL
            && DoesStringMatchConstraint(prp->desc, scp)) {
          rval |= SetStringValue ( &(prp->desc), value, existing_text);
        }
      }
      break;
    case CDSGeneProt_field_mat_peptide_ec_number:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_mat_peptide_aa
            && (prp = sfp->data.value.ptrvalue) != NULL)
        {
          rval |= SetStringsInValNodeStringList (&(prp->ec), scp, value, existing_text);
        }
      }
      break;
    case CDSGeneProt_field_mat_peptide_activity:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->data.choice == SEQFEAT_PROT
            && sfp->idx.subtype == FEATDEF_mat_peptide_aa
            && (prp = sfp->data.value.ptrvalue) != NULL)
        {
          rval |= SetStringsInValNodeStringList (&(prp->activity), scp, value, existing_text);
        }
      }
      break;
    case CDSGeneProt_field_mat_peptide_comment:
      for (vnp = c->prot_list; vnp != NULL; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->idx.subtype == FEATDEF_mat_peptide_aa
            && DoesStringMatchConstraint(sfp->comment, scp))
        {
          rval |= SetStringValue ( &(sfp->comment), value, existing_text);
        }
      }
      break;
  }
  return rval;
}


static MolInfoPtr GetMolInfoForBioseq (BioseqPtr bsp)
{
  MolInfoPtr m = NULL;
  SeqDescrPtr sdp;

  if (bsp == NULL) return NULL;
  sdp = bsp->descr;
  while (sdp != NULL && sdp->choice != Seq_descr_molinfo) {
    sdp = sdp->next;
  }
  if (sdp != NULL) {
    m = (MolInfoPtr) sdp->data.ptrvalue;
  }
  return m;
}
  

static CharPtr GetSequenceQualFromBioseq (BioseqPtr bsp, ValNodePtr field)
{
  CharPtr rval = NULL;
  MolInfoPtr m;

  if (bsp == NULL || field == NULL) return NULL;

  switch (field->choice) {
    case MolinfoField_molecule:
      m = GetMolInfoForBioseq (bsp);
      if (m != NULL) {
        rval = BiomolNameFromBiomol (m->biomol);
      }
      break;
    case MolinfoField_technique:
      m = GetMolInfoForBioseq (bsp);
      if (m != NULL) {
        rval = TechNameFromTech (m->tech);
      }
      break;
    case MolinfoField_completedness:
      m = GetMolInfoForBioseq (bsp);
      if (m != NULL) {
        rval = CompletenessNameFromCompleteness (m->completeness);
      }
      break;
    case MolinfoField_mol_class:
      rval = MolNameFromMol (bsp->mol);
      break;
    case MolinfoField_topology:
      rval = TopologyNameFromTopology (bsp->topology);
      break;
    case MolinfoField_strand:
      rval = StrandNameFromStrand (bsp->strand);
      break;
  }
  if (rval != NULL) rval = StringSave (rval);
  return rval;
}


static Boolean RemoveSequenceQualFromBioseq (BioseqPtr bsp, ValNodePtr field)
{
  MolInfoPtr m;
  Boolean    rval = FALSE;

  if (bsp == NULL || field == NULL) return FALSE;

  switch (field->choice) {
    case MolinfoField_molecule:
      m = GetMolInfoForBioseq (bsp);
      if (m != NULL) {
        m->biomol = 0;
        rval = TRUE;
      }
      break;
    case MolinfoField_technique:
      m = GetMolInfoForBioseq (bsp);
      if (m != NULL) {
        m->tech = 0;
        rval = TRUE;
      }
      break;
    case MolinfoField_completedness:
      m = GetMolInfoForBioseq (bsp);
      if (m != NULL) {
        m->completeness = 0;
        rval = TRUE;
      }
      break;
    case MolinfoField_mol_class:
      bsp->mol = 0;
      rval = TRUE;
      break;
    case MolinfoField_topology:
      bsp->topology = 0;
      rval = TRUE;
      break;
    case MolinfoField_strand:
      bsp->strand = 0;
      rval = TRUE;
      break;
  }
  return rval;
}


static MolInfoPtr AddMolInfoToBioseq (BioseqPtr bsp)
{
  SeqDescrPtr sdp;
  MolInfoPtr  m;

  sdp = CreateNewDescriptorOnBioseq (bsp, Seq_descr_molinfo);
  m = MolInfoNew ();
  sdp->data.ptrvalue = m;
  return m;
}


static Boolean SetSequenceQualOnBioseq (BioseqPtr bsp, ValNodePtr field)
{
  MolInfoPtr m;
  Boolean    rval = FALSE;

  if (bsp == NULL || field == NULL) return FALSE;

  switch (field->choice) {
    case MolinfoField_molecule:
      m = GetMolInfoForBioseq (bsp);
      if (m == NULL) {
        m = AddMolInfoToBioseq (bsp);
      }
      m->biomol = BiomolFromMoleculeType (field->data.intvalue);
      rval = TRUE;
      break;
    case MolinfoField_technique:
      m = GetMolInfoForBioseq (bsp);
      if (m == NULL) {
        m = AddMolInfoToBioseq (bsp);
      }
      m->tech = TechFromTechniqueType (field->data.intvalue);
      rval = TRUE;
      break;
    case MolinfoField_completedness:
      m = GetMolInfoForBioseq (bsp);
      if (m == NULL) {
        m = AddMolInfoToBioseq (bsp);
      }
      m->completeness = CompletenessFromCompletednessType (field->data.intvalue);
      rval = TRUE;
      break;
    case MolinfoField_mol_class:
      bsp->mol = MolFromMoleculeClassType (field->data.intvalue);
      rval = TRUE;
      break;
    case MolinfoField_topology:
      bsp->topology = TopologyFromTopologyType (field->data.intvalue);
      rval = TRUE;
      break;
    case MolinfoField_strand:
      bsp->strand = StrandFromStrandType (field->data.intvalue);
      rval = TRUE;
      break;
  }
  return rval;
}


NLM_EXTERN FieldTypePtr GetFromFieldFromFieldPair (FieldPairTypePtr fieldpair)
{
  SourceQualChoicePtr ss = NULL;
  SourceQualPairPtr sqpp;
  FeatureFieldPairPtr fp;
  FeatureFieldPtr fs;
  FieldTypePtr f = NULL;
  CDSGeneProtFieldPairPtr cp;
  MolinfoFieldPairPtr mp;
  ValNodePtr vnp;

  if (fieldpair == NULL) return NULL;
  switch (fieldpair->choice) {
    case FieldPairType_source_qual:
      sqpp = (SourceQualPairPtr) fieldpair->data.ptrvalue;
      if (sqpp != NULL) {
        ss = ValNodeNew (NULL);
        ss->choice = SourceQualChoice_textqual;
        ss->data.intvalue = sqpp->field_from;
        f = ValNodeNew (NULL);
        f->choice = FieldType_source_qual;
        f->data.ptrvalue = ss;
      }
      break;
    case FieldPairType_feature_field:
      fp = (FeatureFieldPairPtr) fieldpair->data.ptrvalue;
      if (fp != NULL) {
        fs = FeatureFieldNew ();
        fs->type = fp->type;
        fs->field = (FeatQualChoicePtr) AsnIoMemCopy (fp->field_from, (AsnReadFunc) FeatQualChoiceAsnRead, (AsnWriteFunc) FeatQualChoiceAsnWrite);
        f = ValNodeNew (NULL);
        f->choice = FieldType_feature_field;
        f->data.ptrvalue = fs;
      }
      break;
    case FieldPairType_cds_gene_prot:
      cp = (CDSGeneProtFieldPairPtr) fieldpair->data.ptrvalue;
      if (cp != NULL) {
        f = ValNodeNew (NULL);
        f->choice = FieldType_cds_gene_prot;
        f->data.intvalue = cp->field_from;
      }
      break;
    case FieldPairType_molinfo_field:
      mp = (MolinfoFieldPairPtr) fieldpair->data.ptrvalue;
      if (mp != NULL && mp->data.ptrvalue != NULL) {
        vnp = NULL;
        switch (mp->choice) {
          case MolinfoFieldPair_molecule:
            vnp = ValNodeNew (NULL);
            vnp->choice = MolinfoField_molecule;
            vnp->data.intvalue = ((MolinfoMoleculePairPtr)mp->data.ptrvalue)->from;
            break;
          case MolinfoFieldPair_technique:
            vnp = ValNodeNew (NULL);
            vnp->choice = MolinfoField_technique;
            vnp->data.intvalue = ((MolinfoTechniquePairPtr)mp->data.ptrvalue)->from;
            break;
          case MolinfoFieldPair_completedness:
            vnp = ValNodeNew (NULL);
            vnp->choice = MolinfoField_completedness;
            vnp->data.intvalue = ((MolinfoCompletednessPairPtr)mp->data.ptrvalue)->from;
            break;
          case MolinfoFieldPair_mol_class:
            vnp = ValNodeNew (NULL);
            vnp->choice = MolinfoField_mol_class;
            vnp->data.intvalue = ((MolinfoMolClassPairPtr)mp->data.ptrvalue)->from;
            break;
          case MolinfoFieldPair_topology:
            vnp = ValNodeNew (NULL);
            vnp->choice = MolinfoField_topology;
            vnp->data.intvalue = ((MolinfoTopologyPairPtr)mp->data.ptrvalue)->from;
            break;
          case MolinfoFieldPair_strand:
            vnp = ValNodeNew (NULL);
            vnp->choice = MolinfoField_strand;
            vnp->data.intvalue = ((MolinfoStrandPairPtr)mp->data.ptrvalue)->from;
            break;
        }
        if (vnp != NULL) {
          f = ValNodeNew (NULL);
          f->choice = FieldType_molinfo_field;
          f->data.ptrvalue = vnp;
        }
      }
      break;     
  }
  return f;
}


NLM_EXTERN FieldTypePtr GetToFieldFromFieldPair (FieldPairTypePtr fieldpair)
{
  SourceQualChoicePtr ss = NULL;
  SourceQualPairPtr sqpp;
  FeatureFieldPairPtr fp;
  FeatureFieldPtr fs;
  FieldTypePtr f = NULL;
  CDSGeneProtFieldPairPtr cp;
  MolinfoFieldPairPtr     mp;
  ValNodePtr              vnp;

  if (fieldpair == NULL) return NULL;
  switch (fieldpair->choice) {
    case FieldPairType_source_qual:
      sqpp = (SourceQualPairPtr) fieldpair->data.ptrvalue;
      if (sqpp != NULL) {
        ss = ValNodeNew (NULL);
        ss->choice = SourceQualChoice_textqual;
        ss->data.intvalue = sqpp->field_to;
        f = ValNodeNew (NULL);
        f->choice = FieldType_source_qual;
        f->data.ptrvalue = ss;
      }
      break;
    case FieldPairType_feature_field:
      fp = (FeatureFieldPairPtr) fieldpair->data.ptrvalue;
      if (fp != NULL) {
        fs = FeatureFieldNew ();
        fs->type = fp->type;
        fs->field = (FeatQualChoicePtr) AsnIoMemCopy (fp->field_to, (AsnReadFunc) FeatQualChoiceAsnRead, (AsnWriteFunc) FeatQualChoiceAsnWrite);
        f = ValNodeNew (NULL);
        f->choice = FieldType_feature_field;
        f->data.ptrvalue = fs;
      }
      break;
    case FieldPairType_cds_gene_prot:
      cp = (CDSGeneProtFieldPairPtr) fieldpair->data.ptrvalue;
      if (cp != NULL) {
        f = ValNodeNew (NULL);
        f->choice = FieldType_cds_gene_prot;
        f->data.intvalue = cp->field_to;
      }
      break;
    case FieldPairType_molinfo_field:
      mp = (MolinfoFieldPairPtr) fieldpair->data.ptrvalue;
      if (mp != NULL && mp->data.ptrvalue != NULL) {
        vnp = NULL;
        switch (mp->choice) {
          case MolinfoFieldPair_molecule:
            vnp = ValNodeNew (NULL);
            vnp->choice = MolinfoField_molecule;
            vnp->data.intvalue = ((MolinfoMoleculePairPtr)mp->data.ptrvalue)->to;
            break;
          case MolinfoFieldPair_technique:
            vnp = ValNodeNew (NULL);
            vnp->choice = MolinfoField_technique;
            vnp->data.intvalue = ((MolinfoTechniquePairPtr)mp->data.ptrvalue)->to;
            break;
          case MolinfoFieldPair_completedness:
            vnp = ValNodeNew (NULL);
            vnp->choice = MolinfoField_completedness;
            vnp->data.intvalue = ((MolinfoCompletednessPairPtr)mp->data.ptrvalue)->to;
            break;
          case MolinfoFieldPair_mol_class:
            vnp = ValNodeNew (NULL);
            vnp->choice = MolinfoField_mol_class;
            vnp->data.intvalue = ((MolinfoMolClassPairPtr)mp->data.ptrvalue)->to;
            break;
          case MolinfoFieldPair_topology:
            vnp = ValNodeNew (NULL);
            vnp->choice = MolinfoField_topology;
            vnp->data.intvalue = ((MolinfoTopologyPairPtr)mp->data.ptrvalue)->to;
            break;
          case MolinfoFieldPair_strand:
            vnp = ValNodeNew (NULL);
            vnp->choice = MolinfoField_strand;
            vnp->data.intvalue = ((MolinfoStrandPairPtr)mp->data.ptrvalue)->to;
            break;
        }
        if (vnp != NULL) {
          f = ValNodeNew (NULL);
          f->choice = FieldType_molinfo_field;
          f->data.ptrvalue = vnp;
        }
      }
      break;     
  }
  return f;
}


static Uint1 FieldTypeChoiceFromFieldPairTypeChoice (Uint1 field_pair_choice)
{
  Uint1 field_type_choice = 0;

  switch (field_pair_choice) {
    case FieldPairType_source_qual:
      field_type_choice = FieldType_source_qual;
      break;
    case FieldPairType_feature_field:
      field_type_choice = FieldType_feature_field;
      break;
    case FieldPairType_cds_gene_prot:
      field_type_choice = FieldType_cds_gene_prot;
      break;
    case FieldPairType_molinfo_field:
      field_type_choice = FieldType_molinfo_field;
      break;
  }

  return field_type_choice;
}


NLM_EXTERN Uint1 FieldTypeFromAECRAction (AECRActionPtr action)
{
  Uint1 field_type = 0;
  ApplyActionPtr a;
  EditActionPtr  e;
  ConvertActionPtr v;
  CopyActionPtr c;
  SwapActionPtr s;
  RemoveActionPtr r;
  AECRParseActionPtr p;

  if (action == NULL || action->action == NULL || action->action->data.ptrvalue == NULL) {
    return 0;
  }
  switch (action->action->choice) {
    case ActionChoice_apply:
      a = (ApplyActionPtr) action->action->data.ptrvalue;
      if (a->field != NULL) {
        field_type = a->field->choice;
      }
      break;
    case ActionChoice_edit:
      e = (EditActionPtr) action->action->data.ptrvalue;
      if (e->field != NULL) {
        field_type = e->field->choice;
      }
      break;
    case ActionChoice_convert:
      v = (ConvertActionPtr) action->action->data.ptrvalue;
      field_type = FieldTypeChoiceFromFieldPairTypeChoice (v->fields->choice);
      break;
    case ActionChoice_copy:
      c = (CopyActionPtr) action->action->data.ptrvalue;
      field_type = FieldTypeChoiceFromFieldPairTypeChoice (c->fields->choice);
      break;
    case ActionChoice_swap:
      s = (SwapActionPtr) action->action->data.ptrvalue;
      field_type = FieldTypeChoiceFromFieldPairTypeChoice (s->fields->choice);
      break;
    case ActionChoice_remove:
      r = (RemoveActionPtr) action->action->data.ptrvalue;
      if (r->field != NULL) {
        field_type = r->field->choice;
      }
      break;
    case ActionChoice_parse:
      p = (AECRParseActionPtr) action->action->data.ptrvalue;
      field_type = FieldTypeChoiceFromFieldPairTypeChoice (p->fields->choice);
      break;
  }
  return field_type;
}


static CharPtr GetFieldValueForObject (Uint1 choice, Pointer data, FieldTypePtr field, StringConstraintPtr scp)
{
  CharPtr str = NULL;
  FeatureFieldPtr feature_field;

  if (data == NULL || field == NULL || field->data.ptrvalue == NULL) return FALSE;

  switch (field->choice) {
    case FieldType_source_qual :
      str = GetSourceQualFromBioSource (GetBioSourceFromObject (choice, data), (SourceQualChoicePtr) field->data.ptrvalue, scp);
      break;
    case FieldType_feature_field :
      if (choice == OBJ_SEQFEAT) {
        str = GetQualFromFeature ((SeqFeatPtr) data, (FeatureFieldPtr) field->data.ptrvalue, scp);
      }
      break;
    case FieldType_cds_gene_prot :
      if (choice == 0) {
        str = GetFieldValueFromCGPSet ((CGPSetPtr) data, field->data.intvalue, scp);
      } else if (choice == OBJ_SEQFEAT) {
        feature_field = FeatureFieldFromCDSGeneProtField (field->data.intvalue);
        str = GetQualFromFeature ((SeqFeatPtr) data, feature_field, scp);
        feature_field = FeatureFieldFree (feature_field);
      }
      break;
    case FieldType_molinfo_field :
      if (choice == OBJ_BIOSEQ) {
        str = GetSequenceQualFromBioseq ((BioseqPtr) data, field->data.ptrvalue);
      }
      break;
  }
  return str;
}


static Boolean RemoveFieldValueForObject (Uint1 choice, Pointer data, FieldTypePtr field, StringConstraintPtr scp)
{
  Boolean rval = FALSE;
  FeatureFieldPtr feature_field;

  if (data == NULL || field == NULL || field->data.ptrvalue == NULL) return FALSE;

  switch (field->choice) {
    case FieldType_source_qual :
      rval = RemoveSourceQualFromBioSource (GetBioSourceFromObject (choice, data), (SourceQualChoicePtr) field->data.ptrvalue, scp);
      break;
    case FieldType_feature_field :
      if (choice == OBJ_SEQFEAT) {
        rval = RemoveQualFromFeature ((SeqFeatPtr) data, (FeatureFieldPtr) field->data.ptrvalue, scp);
      }
      break;
    case FieldType_cds_gene_prot:
      if (choice == 0) {
        rval = RemoveFieldValueFromCGPSet ((CGPSetPtr) data, field->data.intvalue, scp);
      } else if (choice == OBJ_SEQFEAT) {
        feature_field = FeatureFieldFromCDSGeneProtField (field->data.intvalue);
        rval = RemoveQualFromFeature ((SeqFeatPtr) data, feature_field, scp);
        feature_field = FeatureFieldFree (feature_field);
      }
      break;
    case FieldType_molinfo_field :
      if (choice == OBJ_BIOSEQ) {
        rval = RemoveSequenceQualFromBioseq ((BioseqPtr) data, field->data.ptrvalue);
      }
      break;
  }
  return rval;
}


static Boolean SetFieldValueForObject (Uint1 choice, Pointer data, FieldTypePtr field, StringConstraintPtr scp, CharPtr value, Uint2 existing_text)
{
  Boolean rval = FALSE;
  FeatureFieldPtr feature_field;

  if (data == NULL || field == NULL || field->data.ptrvalue == NULL) return FALSE;

  switch (field->choice) {
    case FieldType_source_qual :
      rval = SetSourceQualInBioSource (GetBioSourceFromObject (choice, data), (SourceQualChoicePtr) field->data.ptrvalue, scp, value, existing_text);
      break;
    case FieldType_feature_field :
      if (choice == OBJ_SEQFEAT) {
        rval = SetQualOnFeature ((SeqFeatPtr) data, (FeatureFieldPtr) field->data.ptrvalue, scp, value, existing_text);
      }
      break;
    case FieldType_cds_gene_prot:
      if (choice == 0) {
        rval = SetFieldValueInCGPSet ((CGPSetPtr) data, field->data.intvalue, scp, value, existing_text);
      } else if (choice == OBJ_SEQFEAT) {
        feature_field = FeatureFieldFromCDSGeneProtField (field->data.intvalue);
        rval = SetQualOnFeature ((SeqFeatPtr) data, feature_field, scp, value, existing_text);
        feature_field = FeatureFieldFree (feature_field);
      }
      break;
    case FieldType_molinfo_field:
      if (choice == OBJ_BIOSEQ) {
        rval = SetSequenceQualOnBioseq ((BioseqPtr) data, field->data.ptrvalue);
      }
      break;
  }
  return rval;
}


static Boolean IsObjectAppropriateForFieldValue (Uint1 choice, Pointer data, FieldTypePtr field)
{
  SeqFeatPtr        sfp;
  SeqDescrPtr       sdp;
  FeatureFieldPtr   fp;
  Boolean rval = FALSE;

  if (data == NULL || field == NULL) return FALSE;

  switch (field->choice) {
    case FieldType_source_qual :
      if (choice == OBJ_SEQFEAT) {
        sfp = (SeqFeatPtr) data;
        if (sfp->data.choice == SEQFEAT_BIOSRC) {
          rval = TRUE;
        }
      } else if (choice == OBJ_SEQDESC) {
        sdp = (SeqDescrPtr) data;
        if (sdp->choice == Seq_descr_source) {
          rval = TRUE;
        }
      }
      break;
    case FieldType_feature_field :
      if (choice == OBJ_SEQFEAT) {
        sfp = (SeqFeatPtr) data;
        fp = (FeatureFieldPtr) field->data.ptrvalue;
        if (fp != NULL && (fp->type == Feature_type_any || GetFeatdefFromFeatureType (fp->type) == sfp->idx.subtype)) {
          rval = TRUE;
        }
      }
      break;
    case FieldType_cds_gene_prot :
      if (choice == 0) {
        rval = TRUE;
      }
      break;
    case FieldType_molinfo_field :
      if (choice == OBJ_BIOSEQ) {
        rval = TRUE;
      }
      break;
  }
  return rval;
}


static Boolean IsObjectAppropriateForFieldPair (Uint1 choice, Pointer data, FieldPairTypePtr fieldpair)
{
  FieldTypePtr f;
  Boolean rval;

  f = GetFromFieldFromFieldPair(fieldpair);
  rval = IsObjectAppropriateForFieldValue(choice, data, f);
  f = FieldTypeFree (f);
  return rval;
}


static Boolean DoFieldTypesMatch (FieldTypePtr field1, FieldTypePtr field2)
{
  Boolean rval = FALSE;
  SourceQualChoicePtr scp1, scp2;
  FeatureFieldPtr fp1, fp2;

  if (field1 == NULL || field2 == NULL) return FALSE;
  if (field1->choice != field2->choice) return FALSE;

  switch (field1->choice) {
    case FieldType_source_qual :
      scp1 = (SourceQualChoicePtr) field1->data.ptrvalue;
      scp2 = (SourceQualChoicePtr) field2->data.ptrvalue;
      if (scp1 != NULL && scp2 != NULL && scp1->choice == scp2->choice) {
        switch (scp1->choice) {
          case SourceQualChoice_textqual:
            if (scp1->data.intvalue == scp2->data.intvalue) {
              rval = TRUE;
            }
            break;
          case SourceQualChoice_location:
          case SourceQualChoice_origin:
            rval = TRUE;
            break;
        }
      }
      break;
    case FieldType_feature_field :
      fp1 = (FeatureFieldPtr) field1->data.ptrvalue;
      fp2 = (FeatureFieldPtr) field2->data.ptrvalue;
      if (fp1 != NULL && fp2 != NULL
          && (fp1->type == fp2->type || fp1->type == Feature_type_any || fp2->type == Feature_type_any)
          && fp1->field != NULL && fp2->field != NULL
          && fp1->field->choice == FeatQualChoice_legal_qual && fp2->field->choice == FeatQualChoice_legal_qual
          && fp1->field->data.intvalue == fp2->field->data.intvalue) {
        rval = TRUE;
      }
      break;
    case FieldType_cds_gene_prot :
      if (field1->data.intvalue == field2->data.intvalue) {
        rval = TRUE;
      }
      break;
  }
  return rval;
}


static Boolean IsNonTextSourceQualPresent (BioSourcePtr biop, Int4 srcqual)
{
  Int4 orgmod_subtype, subsrc_subtype;
  OrgModPtr mod;
  SubSourcePtr ssp;
  Boolean      rval = FALSE;

  if (biop == NULL) return FALSE;

  orgmod_subtype = GetOrgModQualFromSrcQual (srcqual);
  if (orgmod_subtype == -1) {
    subsrc_subtype = GetSubSrcQualFromSrcQual (srcqual);
    for (ssp = biop->subtype; ssp != NULL && !rval; ssp = ssp->next) {
      if (ssp->subtype == subsrc_subtype) {
        rval = TRUE;
      }
    }
  } else {
    if (biop->org != NULL && biop->org->orgname != NULL) {
      for (mod = biop->org->orgname->mod; mod != NULL && !rval; mod = mod->next) {
        if (mod->subtype == orgmod_subtype) {
          rval = TRUE;
        }
      }
    }
  }
  return rval;
}


static Boolean IsSourceQualPresent (BioSourcePtr biop, SourceQualChoicePtr scp)
{
  Boolean rval = FALSE;
  CharPtr   str;

  if (biop == NULL) return FALSE;
  if (scp == NULL) return TRUE;

  switch (scp->choice) {
    case SourceQualChoice_textqual:
      if (IsNonTextSourceQual (scp->data.intvalue)) {
        rval = IsNonTextSourceQualPresent (biop, scp->data.intvalue);
      } else {
        str = GetSourceQualFromBioSource (biop, scp, NULL);
        if (!StringHasNoText (str)) {
          rval = TRUE;
        }
        str = MemFree (str);
      }
      break;
    case SourceQualChoice_location:
      if (biop->genome != 0) {
        rval = TRUE;
      }
      break;
    case SourceQualChoice_origin:
      if (biop->origin != 0) {
        rval = TRUE;
      }
      break;
  }
  return rval;
}


typedef struct objecthasstring
{
  StringConstraintPtr scp;
  Boolean             found;
} ObjectHasStringData, PNTR ObjectHasStringPtr;


static void LIBCALLBACK AsnWriteConstraintCallBack (AsnExpOptStructPtr pAEOS)

{
  CharPtr            pchSource;
  ObjectHasStringPtr ohsp;

  ohsp = (ObjectHasStringPtr) pAEOS->data;
  if (ISA_STRINGTYPE (AsnFindBaseIsa (pAEOS->atp))) 
  {
	  pchSource = (CharPtr) pAEOS->dvp->ptrvalue;
	  ohsp->found |= DoesSingleStringMatchConstraint (pchSource, ohsp->scp);
  }
}


static Boolean DoesObjectMatchStringConstraint (Uint1 choice, Pointer data, StringConstraintPtr scp)

{
  ObjMgrPtr         omp;
  ObjMgrTypePtr     omtp;
  AsnIoPtr          aip;
  AsnExpOptPtr      aeop;
  ObjectHasStringData ohsd;
  SeqFeatPtr          sfp, prot;
  SeqMgrFeatContext   fcontext;
  CharPtr             search_txt;
  CGPSetPtr           c;
  ValNodePtr          vnp;
  Boolean             all_match = TRUE, any_match = FALSE, rval;
  BioseqPtr           protbsp;

  if (data == NULL) return FALSE;
  if (scp == NULL) return TRUE;

  if (choice == 0) {
    /* CDS-Gene-Prot set */
    c = (CGPSetPtr) data;
    for (vnp = c->gene_list; vnp != NULL && (!any_match || all_match); vnp = vnp->next) {
      if (DoesObjectMatchStringConstraint (OBJ_SEQFEAT, vnp->data.ptrvalue, scp)) {
        any_match = TRUE;
      } else {
        all_match = FALSE;
      }
    }
    for (vnp = c->cds_list; vnp != NULL && (!any_match || all_match); vnp = vnp->next) {
      if (DoesObjectMatchStringConstraint (OBJ_SEQFEAT, vnp->data.ptrvalue, scp)) {
        any_match = TRUE;
      } else {
        all_match = FALSE;
      }
    }
    for (vnp = c->mrna_list; vnp != NULL && (!any_match || all_match); vnp = vnp->next) {
      if (DoesObjectMatchStringConstraint (OBJ_SEQFEAT, vnp->data.ptrvalue, scp)) {
        any_match = TRUE;
      } else {
        all_match = FALSE;
      }
    }
    for (vnp = c->prot_list; vnp != NULL && (!any_match || all_match); vnp = vnp->next) {
      if (DoesObjectMatchStringConstraint (OBJ_SEQFEAT, vnp->data.ptrvalue, scp)) {
        any_match = TRUE;
      } else {
        all_match = FALSE;
      }
    }
    if (scp->not_present) {
      rval = all_match;
    } else {
      rval = any_match;
    }        
  } else {
    omp = ObjMgrGet ();
    omtp = ObjMgrTypeFind (omp, choice, NULL, NULL);
    if (omtp == NULL) return FALSE;
    aip = AsnIoNullOpen ();
    aeop = AsnExpOptNew (aip, NULL, NULL, AsnWriteConstraintCallBack);
    ohsd.found = FALSE;
    ohsd.scp = scp;
    if (aeop != NULL) {
      aeop->user_data = (Pointer) &ohsd;
    }
    
    (omtp->asnwrite) (data, aip, NULL);
    
    if (!ohsd.found && omtp->datatype == OBJ_SEQFEAT)
    {
      sfp = (SeqFeatPtr) data;
      if (sfp->data.choice == SEQFEAT_CDREGION) {
        protbsp = BioseqFindFromSeqLoc (sfp->product);
        prot = SeqMgrGetNextFeature (protbsp, NULL, 0, FEATDEF_PROT, &fcontext);
        if (prot != NULL) {
          (omtp->asnwrite) (prot, aip, NULL);
        }
      } else {
        if (SeqMgrFeaturesAreIndexed(sfp->idx.entityID) == 0) {
          SeqMgrIndexFeatures (sfp->idx.entityID, NULL);
        }
        if (sfp->idx.subtype == FEATDEF_tRNA) {
          sfp = SeqMgrGetDesiredFeature (sfp->idx.entityID, NULL, sfp->idx.itemID, 0, sfp, &fcontext);
          ohsd.found = DoesSingleStringMatchConstraint (fcontext.label, ohsd.scp);
          if (!ohsd.found && sfp != NULL && sfp->idx.subtype == FEATDEF_tRNA)
          {
            search_txt = (CharPtr) MemNew ((StringLen (fcontext.label) + 6) * sizeof (Char));
            if (search_txt != NULL)
            {
              sprintf (search_txt, "tRNA-%s", fcontext.label);
              ohsd.found = DoesSingleStringMatchConstraint (search_txt, ohsd.scp);
              search_txt = MemFree (search_txt);
            }
          }
        }
      }
    }
    AsnIoClose (aip);
    if (scp->not_present) {
      rval = !ohsd.found;
    } else {
      rval = ohsd.found;
    }
  }
  return rval;
}


NLM_EXTERN Boolean IsSourceConstraintEmpty (SourceConstraintPtr scp)
{
  if (scp == NULL) return TRUE;

  if (scp->field1 == NULL
      && scp->field2 == NULL
      && IsStringConstraintEmpty(scp->constraint)) {
    return TRUE;
  } else {
    return FALSE;
  }
}

NLM_EXTERN Boolean DoesBiosourceMatchConstraint (BioSourcePtr biop, SourceConstraintPtr scp)
{
  Boolean rval = FALSE;
  CharPtr str1, str2;
  ValNode vn;

  if (biop == NULL) return FALSE;
  if (scp == NULL) return TRUE;

  if (IsStringConstraintEmpty(scp->constraint)) {
    /* looking for qual present */
    if (scp->field1 != NULL && scp->field2 == NULL) {
      rval = IsSourceQualPresent (biop, scp->field1);
    } else if (scp->field2 != NULL && scp->field1 == NULL) {
      rval = IsSourceQualPresent (biop, scp->field2);
    /* looking for quals to match */
    } else if (scp->field1 != NULL && scp->field2 != NULL) {
      str1 = GetSourceQualFromBioSource (biop, scp->field1, NULL);
      str2 = GetSourceQualFromBioSource (biop, scp->field2, NULL);
      if (StringCmp (str1, str2) == 0) {
        rval = TRUE;
      }
      str1 = MemFree (str1);
      str2 = MemFree (str2);
    } else {
      /* nothing specified, automatic match */
      rval = TRUE;
    }
  } else {
    if (scp->field1 != NULL && scp->field2 == NULL) {
      str1 = GetSourceQualFromBioSource (biop, scp->field1, scp->constraint);
      if (str1 == NULL) {
        if (scp->constraint->not_present) {
          str1 = GetSourceQualFromBioSource (biop, scp->field1, NULL);
          if (str1 == NULL) {
            rval = TRUE;
          }
        }
      } else if (!StringHasNoText (str1)) {
        rval = TRUE;
      }
      str1 = MemFree (str1);
    } else if (scp->field2 != NULL && scp->field1 == NULL) {
      str2 = GetSourceQualFromBioSource (biop, scp->field2, scp->constraint);
      if (str2 == NULL) {
        if (scp->constraint->not_present) {
          str2 = GetSourceQualFromBioSource (biop, scp->field2, NULL);
          if (str2 == NULL) {
            rval = TRUE;
          }
        }
      } else if (!StringHasNoText (str2)) {
        rval = TRUE;
      }
      str2 = MemFree (str2);
    } else if (scp->field1 != NULL && scp->field2 != NULL) {
      str1 = GetSourceQualFromBioSource (biop, scp->field1, scp->constraint);
      str2 = GetSourceQualFromBioSource (biop, scp->field2, scp->constraint);
      if (StringCmp (str1, str2) == 0) {
        rval = TRUE;
      }
      str1 = MemFree (str1);
      str2 = MemFree (str2);
    } else {
      /* generic string constraint */
      vn.choice = Seq_descr_source;
      vn.next = NULL;
      vn.extended = 0;
      vn.data.ptrvalue = biop;
      rval = DoesObjectMatchStringConstraint (OBJ_SEQDESC, &vn, scp->constraint);
    }
  }
  return rval;
}


static Boolean DoesCGPSetMatchPseudoConstraint (CGPSetPtr c, CDSGeneProtPseudoConstraintPtr constraint)
{
  Boolean    any_pseudo = FALSE;
  ValNodePtr vnp;
  SeqFeatPtr sfp;
  Boolean    rval = FALSE;

  if (c == NULL) return FALSE;
  if (constraint == NULL) return TRUE;

  switch (constraint->feature) {
    case CDSGeneProt_feature_type_constraint_gene :
      for (vnp = c->gene_list; vnp != NULL && !any_pseudo; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->pseudo) {
          any_pseudo = TRUE;
        }
      }
      break;
    case CDSGeneProt_feature_type_constraint_mRNA :
      for (vnp = c->mrna_list; vnp != NULL && !any_pseudo; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->pseudo) {
          any_pseudo = TRUE;
        }
      }
      break;
    case CDSGeneProt_feature_type_constraint_cds :
      for (vnp = c->mrna_list; vnp != NULL && !any_pseudo; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->pseudo) {
          any_pseudo = TRUE;
        }
      }
      break;
    case CDSGeneProt_feature_type_constraint_prot :
      for (vnp = c->mrna_list; vnp != NULL && !any_pseudo; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->pseudo && sfp->idx.subtype == FEATDEF_PROT) {
          any_pseudo = TRUE;
        }
      }
      break;
    case CDSGeneProt_feature_type_constraint_mat_peptide :
      for (vnp = c->mrna_list; vnp != NULL && !any_pseudo; vnp = vnp->next) {
        sfp = (SeqFeatPtr) vnp->data.ptrvalue;
        if (sfp != NULL && sfp->pseudo && sfp->idx.subtype == FEATDEF_mat_peptide_aa) {
          any_pseudo = TRUE;
        }
      }
      break;
  }

  if ((any_pseudo && constraint->is_pseudo)
      || (!any_pseudo && !constraint->is_pseudo)) {
    rval = TRUE;
  }
  return rval;
}


NLM_EXTERN Boolean IsCDSGeneProtQualConstraintEmpty (CDSGeneProtQualConstraintPtr constraint)
{
  if (constraint == NULL) return TRUE;
  if (constraint->field1 == NULL && constraint->field2 == NULL && IsStringConstraintEmpty (constraint->constraint)) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static Boolean DoesCGPSetMatchQualConstraint (CGPSetPtr c, CDSGeneProtQualConstraintPtr constraint)
{
  Boolean rval = FALSE, any_match = FALSE, all_match = TRUE;
  CharPtr str, str1, str2;

  if (c == NULL) return FALSE;
  if (constraint == NULL) return TRUE;

  if (IsStringConstraintEmpty (constraint->constraint)) {
    /* looking for qual present */
    if (constraint->field1 != NULL && constraint->field2 == NULL) {
      str = GetFieldValueFromCGPSet (c, constraint->field1->data.intvalue, NULL);
      if (str != NULL) {
        rval = TRUE;
        str = MemFree (str);
      }
    } else if (constraint->field2 != NULL && constraint->field1 == NULL) {
      str = GetFieldValueFromCGPSet (c, constraint->field2->data.intvalue, NULL);
      if (str == NULL) {
        rval = FALSE;
      } else {
        str = MemFree (str);
      }
    /* looking for quals to match */
    } else if (constraint->field1 != NULL && constraint->field2 != NULL) {
      str1 = GetFieldValueFromCGPSet (c, constraint->field1->data.intvalue, NULL);
      str2 = GetFieldValueFromCGPSet (c, constraint->field2->data.intvalue, NULL);
      if (StringCmp (str1, str2) == 0) {
        rval = TRUE;
      }
      str1 = MemFree (str1);
      str2 = MemFree (str2);
    } else {
      /* nothing specified, automatic match */
      rval = TRUE;
    }
  } else {
    if (constraint->field1 != NULL && constraint->field2 == NULL) {
      str1 = GetFieldValueFromCGPSet (c, constraint->field1->data.intvalue, constraint->constraint);
      if (str1 == NULL) {
        if (constraint->constraint->not_present) {
          str1 = GetFieldValueFromCGPSet (c, constraint->field1->data.intvalue, NULL);
          if (str1 == NULL) {
            rval = TRUE;
          }
        }
      } else if (!StringHasNoText (str1)) {
        rval = TRUE;
      }
      str1 = MemFree (str1);
    } else if (constraint->field2 != NULL && constraint->field1 == NULL) {
      str2 = GetFieldValueFromCGPSet (c, constraint->field2->data.intvalue, constraint->constraint);
      if (str2 == NULL) {
        if (constraint->constraint->not_present) {
          str2 = GetFieldValueFromCGPSet (c, constraint->field2->data.intvalue, NULL);
          if (str2 == NULL) {
            rval = TRUE;
          }
        }
      } else if (!StringHasNoText (str2)) {
        rval = TRUE;
      }
      str2 = MemFree (str2);
    } else if (constraint->field1 != NULL && constraint->field2 != NULL) {
      str1 = GetFieldValueFromCGPSet (c, constraint->field1->data.intvalue, constraint->constraint);
      str2 = GetFieldValueFromCGPSet (c, constraint->field2->data.intvalue, constraint->constraint);
      if (StringCmp (str1, str2) == 0) {
        rval = TRUE;
      }
      str1 = MemFree (str1);
      str2 = MemFree (str2);
    } else {
      /* generic string constraint */
      rval = DoesObjectMatchStringConstraint (0, c, constraint->constraint);
    }
  }
  return rval;
}


NLM_EXTERN Boolean IsSequenceConstraintEmpty (SequenceConstraintPtr constraint)
{
  if (constraint == NULL) return TRUE;
  if (constraint->seqtype != NULL && constraint->seqtype->choice != SequenceConstraintMolTypeConstraint_any) return FALSE;
  if (constraint->feature != Feature_type_any) return FALSE;
  if (!IsStringConstraintEmpty (constraint->id)) return FALSE;
  return TRUE;
}


extern Boolean DoesSeqIDListMeetStringConstraint (SeqIdPtr sip, StringConstraintPtr string_constraint)
{
  Char       id [41];
  CharPtr    cp, cp_dst;
  SeqIdPtr   tmp;
  Boolean    match, changed;

  if (sip == NULL) 
  {
    return FALSE;
  }
  if (string_constraint == NULL)
  {
    return TRUE;
  }
  
  while (sip != NULL)
  {
    /* temporary disconnect ID from list */
    tmp = sip->next;
    sip->next = NULL;
    id [0] = '\0';
    SeqIdWrite (sip, id, PRINTID_FASTA_LONG, sizeof (id) - 1);
    match = DoesSingleStringMatchConstraint (id, string_constraint);
    if (!match) 
    {
      changed = FALSE;
      /* remove terminating pipe character */
      if (id[StringLen(id) - 1] == '|') 
      {
        id[StringLen(id) - 1] = 0;
        changed = TRUE;
      }
      /* remove leading pipe identifier */
      cp = StringChr (id, '|');
      if (cp != NULL)
      {
        changed = TRUE;
        cp++;
        cp_dst = id;
        while (*cp != 0) 
        {
          *cp_dst = *cp;
          cp_dst++;
          cp++;
        }
        *cp_dst = 0;
      }  
      if (changed) 
      {
        match = DoesSingleStringMatchConstraint (id, string_constraint);
      }

      /* if search text doesn't have ., try ID without version */
      if (!match && StringChr (string_constraint->match_text, '.') == NULL) 
      {
        cp = StringChr (id, '.');
        if (cp != NULL) 
        {
          *cp = 0;
          match = DoesSingleStringMatchConstraint (id, string_constraint);
        }
      }       
    }
    sip->next = tmp;

    if (match)
    {
      if (string_constraint->not_present)
      {
        return FALSE;
      }
      else
      {
        return TRUE;
      }
    }
    sip = sip->next;
  }
  if (string_constraint->not_present)
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}


typedef struct rnatypebiomol {
  Int4 rnatype;
  Uint1 biomol;
  CharPtr rnamolname;
} RnaTypeBiomolData, PNTR RnaTypeBiomolPtr;

static RnaTypeBiomolData rna_type_biomol[] = {
{ Sequence_constraint_rnamol_genomic , MOLECULE_TYPE_GENOMIC, "Genomic RNA" } ,
{ Sequence_constraint_rnamol_precursor_RNA , MOLECULE_TYPE_PRE_MRNA , "Precursor RNA" } ,
{ Sequence_constraint_rnamol_mRNA , MOLECULE_TYPE_MRNA , "mRNA [cDNA]" } ,
{ Sequence_constraint_rnamol_rRNA , MOLECULE_TYPE_RRNA , "Ribosomal RNA" } ,
{ Sequence_constraint_rnamol_tRNA , MOLECULE_TYPE_TRNA , "Transfer RNA" } ,
{ Sequence_constraint_rnamol_snRNA , MOLECULE_TYPE_SNRNA , "Small nuclear RNA" } ,
{ Sequence_constraint_rnamol_scRNA , MOLECULE_TYPE_SCRNA , "Small cytoplasmic RNA" } ,
{ Sequence_constraint_rnamol_genomic_mRNA , MOLECULE_TYPE_GENOMIC_MRNA_MIX , "Genomic-mRNA" } ,
{ Sequence_constraint_rnamol_cRNA , MOLECULE_TYPE_CRNA , "cRNA" } ,
{ Sequence_constraint_rnamol_snoRNA , MOLECULE_TYPE_SNORNA , "Small nucleolar RNA" } ,
{ Sequence_constraint_rnamol_transcribed_RNA , MOLECULE_TYPE_TRANSCRIBED_RNA , "Transcribed RNA" } ,
{ Sequence_constraint_rnamol_ncRNA , MOLECULE_TYPE_NCRNA , "Non-coding  RNA" } ,
{ Sequence_constraint_rnamol_transfer_messenger_RNA , MOLECULE_TYPE_TMRNA , "Transfer-messenger RNA" } } ;

#define NUM_rna_type_biomol sizeof (rna_type_biomol) / sizeof (RnaTypeBiomolData)


NLM_EXTERN Uint1 GetBiomolForRnaType (Int4 rnatype) 
{
  Int4 i;

  for (i = 0; i <  NUM_rna_type_biomol; i++) {
    if (rna_type_biomol[i].rnatype == rnatype) {
      return rna_type_biomol[i].biomol;
    }
  }
  return 0;
}


NLM_EXTERN CharPtr GetBiomolNameForRnaType (Int4 rnatype)
{
  Int4 i;

  for (i = 0; i <  NUM_rna_type_biomol; i++) {
    if (rna_type_biomol[i].rnatype == rnatype) {
      return rna_type_biomol[i].rnamolname;
    }
  }
  return "invalid RNA type";
}

NLM_EXTERN void AddAllRNASubtypesToChoiceList (ValNodePtr PNTR field_list)
{
  Int4 i;

  if (field_list == NULL) return;

  ValNodeAddPointer (field_list, Sequence_constraint_rnamol_any, StringSave ("Any RNA"));
  for (i = 0; i < NUM_rna_type_biomol; i++) {
    ValNodeAddPointer (field_list, rna_type_biomol[i].rnatype, StringSave (rna_type_biomol[i].rnamolname));
  }
}


static Boolean DoesSequenceMatchSequenceConstraint (BioseqPtr bsp, SequenceConstraintPtr constraint)
{
  SeqFeatPtr sfp;
  SeqMgrFeatContext fcontext;
  SeqDescrPtr sdp;
  SeqMgrDescContext dcontext;
  MolInfoPtr mip;
  
  if (bsp == NULL) return FALSE;
  if (IsSequenceConstraintEmpty (constraint)) return TRUE;

  if (constraint->seqtype != NULL && constraint->seqtype->choice != SequenceConstraintMolTypeConstraint_any) {
    switch (constraint->seqtype->choice) {
      case SequenceConstraintMolTypeConstraint_nucleotide :
        if (ISA_aa (bsp->mol)) {
          return FALSE;
        }
        break;
      case SequenceConstraintMolTypeConstraint_dna :
        if (bsp->mol != Seq_mol_dna) {
          return FALSE;
        }
        break;
      case SequenceConstraintMolTypeConstraint_rna :
        if (bsp->mol != Seq_mol_rna) {
          return FALSE;
        }
        if (constraint->seqtype->data.intvalue != Sequence_constraint_rnamol_any) {
          sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
          if (sdp == NULL || sdp->data.ptrvalue == NULL || sdp->choice != Seq_descr_molinfo) {
            return FALSE;
          }
          mip = (MolInfoPtr) sdp->data.ptrvalue;
          if (GetBiomolForRnaType (constraint->seqtype->data.intvalue) != mip->biomol) {
            return FALSE;
          }
        }
        break;
      case SequenceConstraintMolTypeConstraint_protein :
        if (!ISA_aa (bsp->mol)) {
          return FALSE;
        }
        break;
    }
  }

  if (constraint->feature != Feature_type_any) {
    sfp = SeqMgrGetNextFeature (bsp, NULL, 0, GetFeatdefFromFeatureType (constraint->feature), &fcontext);
    if (sfp == NULL) {
      return FALSE;
    }
  }

  if (!IsStringConstraintEmpty (constraint->id) && !DoesSeqIDListMeetStringConstraint (bsp->id, constraint->id)) {
    return FALSE;
  }
  return TRUE;
}

static Boolean DoesSequenceInSetMatchSequenceConstraint (BioseqSetPtr bssp, SequenceConstraintPtr constraint)
{
  Boolean       rval = FALSE;
  SeqEntryPtr   sep;

  if (bssp == NULL) return FALSE;
  if (IsSequenceConstraintEmpty (constraint)) return TRUE;
  
  for (sep = bssp->seq_set; sep != NULL && !rval; sep = sep->next) {
    if (IS_Bioseq (sep)) {
      rval = DoesSequenceMatchSequenceConstraint ((BioseqPtr) sep->data.ptrvalue, constraint);
    } else if (IS_Bioseq_set (sep)) {
      rval = DoesSequenceInSetMatchSequenceConstraint ((BioseqSetPtr) sep->data.ptrvalue, constraint);
    }
  }
  return rval;
}


static Boolean DoesObjectMatchSequenceConstraint (Uint1 choice, Pointer data, SequenceConstraintPtr constraint)
{
  BioseqPtr bsp;
  SeqDescrPtr sdp;
  ObjValNodePtr ovp;
  Boolean       rval = FALSE;

  if (data == NULL) return FALSE;
  if (IsSequenceConstraintEmpty (constraint)) return TRUE;

  bsp = GetSequenceForObject (choice, data);
  if (bsp == NULL) {
    if (choice == OBJ_SEQDESC) {
      sdp = (SeqDescrPtr) data;
      if (sdp->extended) {
        ovp = (ObjValNodePtr) sdp;
        if (ovp->idx.parenttype == OBJ_BIOSEQSET && ovp->idx.parentptr != NULL) {
          rval = DoesSequenceInSetMatchSequenceConstraint ((BioseqSetPtr) ovp->idx.parentptr, constraint);
        }
      }
    }
  } else {
    rval = DoesSequenceMatchSequenceConstraint (bsp, constraint);
  }
  return rval; 
}


static Boolean DoesObjectMatchConstraint (Uint1 choice, Pointer data, ConstraintChoicePtr constraint)
{
  Boolean rval = TRUE;

  if (data == NULL) return FALSE;
  if (constraint == NULL) return TRUE;

  switch (constraint->choice) {
    case ConstraintChoice_string :
      rval = DoesObjectMatchStringConstraint (choice, data, constraint->data.ptrvalue);
      break;
    case ConstraintChoice_location :
      rval = DoesObjectMatchLocationConstraint (choice, data, constraint->data.ptrvalue);
      break;
    case ConstraintChoice_source :
      rval = DoesBiosourceMatchConstraint (GetBioSourceFromObject (choice, data), constraint->data.ptrvalue);
      break;
    case ConstraintChoice_cdsgeneprot_qual :
      if (choice == 0) {
        rval = DoesCGPSetMatchQualConstraint (data, constraint->data.ptrvalue);
      } else {
        rval = FALSE;
      }
      break;
    case ConstraintChoice_cdsgeneprot_pseudo :
      if (choice == 0) {
        rval = DoesCGPSetMatchPseudoConstraint (data, constraint->data.ptrvalue);
      } else {
        rval = FALSE;
      }
      break;
    case ConstraintChoice_sequence :
      rval = DoesObjectMatchSequenceConstraint (choice, data, constraint->data.ptrvalue);
      break;
  }
  return rval;
}


static Boolean DoesObjectMatchConstraintChoiceSet (Uint1 choice, Pointer data, ConstraintChoiceSetPtr csp)
{
  Boolean rval = TRUE;

  if (data == NULL) return FALSE;

  while (csp != NULL && rval) {
    rval = DoesObjectMatchConstraint (choice, data, csp);
    csp = csp->next;
  }
  return rval;
}


NLM_EXTERN StringConstraintPtr FindStringConstraintInConstraintSetForField (FieldTypePtr field, ConstraintChoiceSetPtr csp)
{
  StringConstraintPtr scp = NULL;
  ConstraintChoicePtr constraint;
  SourceConstraintPtr source_constraint;
  CDSGeneProtQualConstraintPtr cgp_constraint;

  while (csp != NULL) {
    constraint = (ConstraintChoicePtr) csp->data.ptrvalue;
    switch (constraint->choice) {
      case ConstraintChoice_string :
        scp = constraint->data.ptrvalue;
        break;
      case ConstraintChoice_source :
        source_constraint = (SourceConstraintPtr) constraint->data.ptrvalue;
        if (source_constraint != NULL && source_constraint->constraint != NULL
            && ((source_constraint->field1 != NULL
                 && DoFieldTypesMatch (field, source_constraint->field1))
                || (source_constraint->field2 != NULL
                 && DoFieldTypesMatch (field, source_constraint->field2)))) {
            scp = source_constraint->constraint;
        } 
      break;
      case ConstraintChoice_cdsgeneprot_qual :
        cgp_constraint = (CDSGeneProtQualConstraintPtr) field->data.ptrvalue;
        if (field->choice == FieldType_cds_gene_prot
            && cgp_constraint != NULL && cgp_constraint->constraint != NULL
            && ((cgp_constraint->field1 != NULL && cgp_constraint->field1->data.intvalue == field->data.intvalue)
                || (cgp_constraint->field2 != NULL && cgp_constraint->field2->data.intvalue == field->data.intvalue))) {
          scp = cgp_constraint->constraint;
        }
        break;
    }
    csp = csp->next;
  }
  return scp;
}


NLM_EXTERN StringConstraintPtr FindStringConstraintInConstraintSetForFieldPair (FieldPairTypePtr fieldpair, ConstraintChoiceSetPtr csp)
{
  StringConstraintPtr scp;
  FieldTypePtr f;

  f = GetFromFieldFromFieldPair (fieldpair);
  scp = FindStringConstraintInConstraintSetForField (f, csp);
  f = FieldTypeFree (f);
  return scp;
}
 

NLM_EXTERN StringConstraintPtr StringConstraintFromFieldEdit (FieldEditPtr edit)
{
  StringConstraintPtr scp;

  if (edit == NULL || edit->find_txt == NULL) return NULL;
  scp = StringConstraintNew ();
  scp->match_text = StringSave (edit->find_txt);

  switch (edit->location) {
    case Field_edit_location_anywhere :
      scp->match_location = String_location_contains;
      break;
    case Field_edit_location_beginning :
      scp->match_location = String_location_starts;
      break;
    case Field_edit_location_end :
      scp->match_location = String_location_ends;
      break;
  }

  scp->case_sensitive = TRUE;
  scp->whole_word = FALSE;
  scp->not_present = FALSE;

  return scp;
}


static CharPtr ApplyEditToString (CharPtr str, FieldEditPtr edit)
{
  CharPtr cp_found, new_str;
  Int4 found_len, replace_len, new_len;

  if (edit == NULL) return StringSave (str);

  str = StringSave (str);
  cp_found = StringISearch (str, edit->find_txt);

  found_len = StringLen (edit->find_txt);
  replace_len = StringLen (edit->repl_txt);
  if (edit->location == Field_edit_location_beginning
      && cp_found != str) {
    cp_found = NULL;
  } 
  while (cp_found != NULL)
  {
    if (edit->location == Field_edit_location_end
        && cp_found != str + StringLen (str) - found_len) {
      cp_found = StringISearch (cp_found + found_len, edit->find_txt);
    } else {
      new_len = StringLen (str) + 1 - found_len + replace_len;
      new_str = (CharPtr) MemNew (new_len * sizeof (Char));
      if (new_str != NULL)
      {
        if (cp_found != str)
        {
          StringNCpy (new_str, str, cp_found - str);
        }
        StringCat (new_str, edit->repl_txt);
        StringCat (new_str, cp_found + found_len);
        cp_found = new_str + (cp_found - str) + replace_len;
        str = MemFree (str);
        str = new_str;
      }
      cp_found = StringISearch (cp_found, edit->find_txt);
    }
  }
  return str;
}


typedef struct objectcollection {
  AECRActionPtr action;
  ValNodePtr object_list;
} ObjectCollectionData, PNTR ObjectCollectionPtr;


static void AECRActionObjectCollectionItemCallback (Uint1 objecttype, Pointer objectdata, ObjectCollectionPtr o)
{
  ApplyActionPtr a;
  EditActionPtr e;
  ConvertActionPtr v;
  CopyActionPtr c;
  SwapActionPtr s;
  RemoveActionPtr r;
  AECRParseActionPtr p;
  CharPtr str, portion;
  StringConstraintPtr scp;
  FieldTypePtr field_from = NULL, field_to = NULL;

  if (objectdata == NULL || o == NULL) return;

  /* check to make sure object is appropriate for field and meets filter */
  switch (o->action->action->choice) {
    case ActionChoice_apply :
      a = (ApplyActionPtr) o->action->action->data.ptrvalue;
      if (a != NULL
          && IsObjectAppropriateForFieldValue (objecttype, objectdata, a->field)
          && DoesObjectMatchConstraintChoiceSet (objecttype, objectdata, o->action->constraint)) {
        ValNodeAddPointer (&(o->object_list), objecttype, objectdata);
      }
      break;
    case ActionChoice_edit :
      e = (EditActionPtr) o->action->action->data.ptrvalue;
      if (e != NULL
          && IsObjectAppropriateForFieldValue (objecttype, objectdata, e->field)
          && DoesObjectMatchConstraintChoiceSet (objecttype, objectdata, o->action->constraint)) {
        scp = StringConstraintFromFieldEdit (e->edit);
        str = GetFieldValueForObject (objecttype, objectdata, e->field, scp);
        if (!StringHasNoText (str)) {
          ValNodeAddPointer (&(o->object_list), objecttype, objectdata);
        }
        str = MemFree (str);
      }
      break;
    case ActionChoice_convert :
      v = (ConvertActionPtr) o->action->action->data.ptrvalue;
      if (v != NULL
          && (field_from = GetFromFieldFromFieldPair(v->fields)) != NULL
          && IsObjectAppropriateForFieldValue (objecttype, objectdata, field_from)
          && DoesObjectMatchConstraintChoiceSet (objecttype, objectdata, o->action->constraint)) {
        scp = FindStringConstraintInConstraintSetForField (field_from, o->action->constraint);
        str = GetFieldValueForObject (objecttype, objectdata, field_from, scp);
        if (!StringHasNoText (str)) {
          ValNodeAddPointer (&(o->object_list), objecttype, objectdata);
        }
        str = MemFree (str);
      }
      field_from = FieldTypeFree (field_from);
      break;
    case ActionChoice_copy :
      c = (CopyActionPtr) o->action->action->data.ptrvalue;
      if (c != NULL
          && (field_from = GetFromFieldFromFieldPair(c->fields)) != NULL
          && (field_to = GetFromFieldFromFieldPair(c->fields)) != NULL
          && IsObjectAppropriateForFieldValue (objecttype, objectdata, field_from)
          && IsObjectAppropriateForFieldValue (objecttype, objectdata, field_to)
          && DoesObjectMatchConstraintChoiceSet (objecttype, objectdata, o->action->constraint)) {
        ValNodeAddPointer (&(o->object_list), objecttype, objectdata);
      }
      field_from = FieldTypeFree (field_from);
      field_to = FieldTypeFree (field_to);
      break;
    case ActionChoice_swap :
      s = (SwapActionPtr) o->action->action->data.ptrvalue;
      if (s != NULL
          && (field_from = GetFromFieldFromFieldPair(s->fields)) != NULL
          && (field_to = GetFromFieldFromFieldPair(s->fields)) != NULL
          && IsObjectAppropriateForFieldValue (objecttype, objectdata, field_from)
          && IsObjectAppropriateForFieldValue (objecttype, objectdata, field_to)
          && DoesObjectMatchConstraintChoiceSet (objecttype, objectdata, o->action->constraint)) {
        ValNodeAddPointer (&(o->object_list), objecttype, objectdata);
      }
      field_from = FieldTypeFree (field_from);
      field_to = FieldTypeFree (field_to);
      break;
    case ActionChoice_remove :
      r = (RemoveActionPtr) o->action->action->data.ptrvalue;
      if (r != NULL
          && IsObjectAppropriateForFieldValue (objecttype, objectdata, r->field)
          && DoesObjectMatchConstraintChoiceSet (objecttype, objectdata, o->action->constraint)) {
        ValNodeAddPointer (&(o->object_list), objecttype, objectdata);
      }
      break;
    case ActionChoice_parse :
      p = (AECRParseActionPtr) o->action->action->data.ptrvalue;
      if (p != NULL
          && (field_from = GetFromFieldFromFieldPair(p->fields)) != NULL
          && (field_to = GetFromFieldFromFieldPair(p->fields)) != NULL
          && IsObjectAppropriateForFieldValue (objecttype, objectdata, field_from)
          && IsObjectAppropriateForFieldValue (objecttype, objectdata, field_to)
          && DoesObjectMatchConstraintChoiceSet (objecttype, objectdata, o->action->constraint)) {
        scp = FindStringConstraintInConstraintSetForField (field_from, o->action->constraint);
        str = GetFieldValueForObject (objecttype, objectdata, field_from, scp);
        portion = GetTextPortionFromString (str, p->portion);
        if (!StringHasNoText (portion)) {
          ValNodeAddPointer (&(o->object_list), objecttype, objectdata);
        }
        portion = MemFree (portion);
        str = MemFree (str);
      }
      field_from = FieldTypeFree (field_from);
      field_to = FieldTypeFree (field_to);
      break;
  }

}


static void AECRActionObjectCollectionFeatureCallback (SeqFeatPtr sfp, Pointer data)
{
  ObjectCollectionPtr o;
  if (sfp == NULL || data == NULL) return;

  o = (ObjectCollectionPtr) data;
  AECRActionObjectCollectionItemCallback (OBJ_SEQFEAT, sfp, o);

}


static void AECRActionObjectCollectionDescriptorCallback (SeqDescrPtr sdp, Pointer data)
{
  ObjectCollectionPtr o;

  if (sdp == NULL || data == NULL) return;

  o = (ObjectCollectionPtr) data;
  AECRActionObjectCollectionItemCallback (OBJ_SEQDESC, sdp, o);
}


static void AECRObjectCollectionBioseqCallback (BioseqPtr bsp, Pointer data)
{
  ObjectCollectionPtr o;

  if (bsp == NULL || data == NULL) return;

  o = (ObjectCollectionPtr) data;
  AECRActionObjectCollectionItemCallback (OBJ_BIOSEQ, bsp, o);
}


NLM_EXTERN ValNodePtr GetObjectListForAECRAction (SeqEntryPtr sep, AECRActionPtr action)
{
  ObjectCollectionData ocd;

  ocd.action = action;
  ocd.object_list = NULL;

  if (action == NULL) return NULL;
  if (FieldTypeFromAECRAction (action) == FieldType_molinfo_field) {
    VisitBioseqsInSep (sep, &ocd, AECRObjectCollectionBioseqCallback);
  } else {
    VisitFeaturesInSep (sep, &ocd, AECRActionObjectCollectionFeatureCallback);
    VisitDescriptorsInSep (sep, &ocd, AECRActionObjectCollectionDescriptorCallback);
  }
  return ocd.object_list;
}


typedef struct buildcgpset
{
  ValNodePtr cds_list;
  ValNodePtr mrna_list;
  ValNodePtr gene_list;
} BuildCGPSetData, PNTR BuildCGPSetPtr;

static void BuildCGPSetCallback (SeqFeatPtr sfp, Pointer userdata)
{
  BuildCGPSetPtr b;

  if (sfp == NULL || sfp->idx.deleteme || userdata == NULL) return;
  b = (BuildCGPSetPtr) userdata;
  if (sfp->data.choice == SEQFEAT_CDREGION)
  {
    ValNodeAddPointer (&(b->cds_list), OBJ_SEQFEAT, sfp);
  }
  else if (sfp->data.choice == SEQFEAT_GENE)
  {
    ValNodeAddPointer (&(b->gene_list), OBJ_SEQFEAT, sfp);
  }
  else if (sfp->idx.subtype == FEATDEF_mRNA)
  {
    ValNodeAddPointer (&(b->mrna_list), OBJ_SEQFEAT, sfp);
  }
  else if (SeqMgrGetGeneXref (sfp) != NULL)
  {
    ValNodeAddPointer (&(b->gene_list), OBJ_SEQFEAT, sfp);
  }
}


static CGPSetPtr BuildCGPSetFromCodingRegion (SeqFeatPtr cds, BoolPtr indexing_needed)
{
  SeqMgrFeatContext fcontext;
  SeqFeatPtr        gene = NULL, mrna, prot;
  BioseqPtr         protbsp;
  CGPSetPtr         cdsp;
  ProtRefPtr        prp;

  if (cds == NULL || cds->data.choice != SEQFEAT_CDREGION) return NULL;

  cdsp = (CGPSetPtr) MemNew (sizeof (CGPSetData));
  ValNodeAddPointer (&(cdsp->cds_list), 0, cds);

  gene = GetGeneForFeature (cds);
  if (gene != NULL)
  {
    ValNodeAddPointer (&(cdsp->gene_list), 0, gene);
    /* mark gene, so that we'll know it isn't lonely */
    gene->idx.deleteme = TRUE;
  }

  mrna = SeqMgrGetOverlappingmRNA (cds->location, &fcontext);
  if (mrna != NULL)
  {
    ValNodeAddPointer (&(cdsp->mrna_list), 0, mrna);
    /* mark mrna, so that we'll know it's already in a set */
    mrna->idx.deleteme = TRUE;
  }

  if (cds->product != NULL)
  {
    protbsp = BioseqFindFromSeqLoc (cds->product);
    if (protbsp != NULL)
    {
      prot = SeqMgrGetNextFeature (protbsp, NULL, SEQFEAT_PROT, FEATDEF_PROT, &fcontext);
      /* if there is no full-length protein feature, make one */
      if (prot == NULL)
      {
        prp = ProtRefNew ();
        prot = CreateNewFeatureOnBioseq (protbsp, SEQFEAT_PROT, NULL);
        if (prot != NULL)
        {
          prot->data.value.ptrvalue = prp;
          if (indexing_needed != NULL)
          {
            *indexing_needed = TRUE;
          }
        }
      }
      if (prot != NULL)
      {
        ValNodeAddPointer (&(cdsp->prot_list), 0, prot);
      }
      
      /* also add in mat_peptides from protein feature */
      prot = SeqMgrGetNextFeature (protbsp, NULL, SEQFEAT_PROT, FEATDEF_mat_peptide_aa, &fcontext);
      while (prot != NULL)
      {
        ValNodeAddPointer (&(cdsp->prot_list), 0, prot);
        prot = SeqMgrGetNextFeature (protbsp, prot, SEQFEAT_PROT, FEATDEF_mat_peptide_aa, &fcontext);
      }
    }
  }  
  return cdsp;
}


static CGPSetPtr BuildCGPSetFrommRNA (SeqFeatPtr mrna)
{
  SeqFeatPtr        gene;
  CGPSetPtr          cdsp;

  if (mrna == NULL || mrna->idx.deleteme || mrna->idx.subtype != FEATDEF_mRNA) return NULL;

  cdsp = (CGPSetPtr) MemNew (sizeof (CGPSetData));
  ValNodeAddPointer (&(cdsp->mrna_list), 0, mrna);

  gene = GetGeneForFeature (mrna);
  if (gene != NULL)
  {
    ValNodeAddPointer (&(cdsp->gene_list), 0, gene);
    /* mark gene, so that we'll know it isn't lonely */
    gene->idx.deleteme = TRUE;
  }

  return cdsp;
}


static void UnmarkFeatureList (ValNodePtr list)
{
  SeqFeatPtr sfp;

  while (list != NULL)
  {
    sfp = list->data.ptrvalue;
    if (sfp != NULL)
    {
      sfp->idx.deleteme = FALSE;
    }
    list = list->next;
  }
}


static ValNodePtr BuildCGPSetList (Uint2 entityID, ValNodePtr constraint)
{
  SeqEntryPtr    sep;
  BuildCGPSetData b;
  CGPSetPtr       cdsp;
  ValNodePtr     vnp, vnp_next, vnp_prev;
  ValNodePtr     cdset_list = NULL;
  SeqFeatPtr     cds, gene, mrna;
  Boolean        need_indexing = FALSE;
  
  sep = GetTopSeqEntryForEntityID (entityID);

  b.cds_list = NULL;
  b.gene_list = NULL;
  b.mrna_list = NULL;
  
  VisitFeaturesInSep (sep, &b, BuildCGPSetCallback);

  /* build cdsets that have coding regions */
  for (vnp = b.cds_list; vnp != NULL; vnp = vnp->next)
  {
    cds = (SeqFeatPtr) vnp->data.ptrvalue;
    if (cds == NULL) continue;
    cdsp = BuildCGPSetFromCodingRegion (cds, &need_indexing);
    if (cdsp != NULL)
    {
      ValNodeAddPointer (&cdset_list, 0, cdsp);
    }
  }
  if (need_indexing)
  {
    /* indexing because we have created full-length protein features */
    SeqMgrIndexFeatures (entityID, NULL);
  }

  /* build cdsets for mrna features that don't have coding regions */
  for (vnp = b.mrna_list; vnp != NULL; vnp = vnp->next)
  {
    mrna = (SeqFeatPtr) vnp->data.ptrvalue;
    if (mrna == NULL || mrna->idx.deleteme) continue;
    cdsp = BuildCGPSetFrommRNA (mrna);
    if (cdsp != NULL)
    {
      ValNodeAddPointer (&cdset_list, 0, cdsp);
    }
  }

  /* build cdsets for lonely genes / features with gene xrefs that are not coding regions or mrnas */
  for (vnp = b.gene_list; vnp != NULL; vnp = vnp->next)
  {
    gene = (SeqFeatPtr) vnp->data.ptrvalue;
    if (gene == NULL || gene->idx.deleteme) continue;
    cdsp = CGPSetNew ();
    ValNodeAddPointer (&(cdsp->gene_list), 0, gene);
    ValNodeAddPointer (&cdset_list, 0, cdsp);
  }

  /* now unmark features */
  UnmarkFeatureList (b.cds_list);
  UnmarkFeatureList (b.mrna_list);
  UnmarkFeatureList (b.gene_list);

  b.cds_list = ValNodeFree (b.cds_list);
  b.mrna_list = ValNodeFree (b.mrna_list);
  b.gene_list = ValNodeFree (b.gene_list);

  /* now remove sets that don't match our choice constraint */
  vnp_prev = NULL;
  for (vnp = cdset_list; vnp != NULL; vnp = vnp_next)
  {
    vnp_next = vnp->next;
    if (!DoesObjectMatchConstraintChoiceSet (0, vnp->data.ptrvalue, constraint))
    {
      if (vnp_prev == NULL)
      {
        cdset_list = vnp->next;
      }
      else
      {
        vnp_prev->next = vnp->next;
      }
      vnp->next = NULL;
      FreeCGPSetList (vnp);     
    }
    else
    {
      vnp_prev = vnp;
    }
  }
  
  return cdset_list;
}


NLM_EXTERN Int4 DoApplyActionToObjectList (ApplyActionPtr action, ValNodePtr object_list, StringConstraintPtr scp)
{
  ValNodePtr vnp;
  Int4       num_succeed = 0, num_fail = 0;

  if (action == NULL || object_list == NULL) return 0;

  for (vnp = object_list; vnp != NULL; vnp = vnp->next) {
    if (SetFieldValueForObject (vnp->choice, vnp->data.ptrvalue, action->field, scp, action->value, action->existing_text)) {
      num_succeed ++;
    } else {
      num_fail++;
    }
  }
  return num_succeed;
}


NLM_EXTERN Int4 DoEditActionToObjectList (EditActionPtr action, ValNodePtr object_list)
{
  ValNodePtr vnp;
  Int4       num_succeed = 0, num_fail = 0;
  StringConstraintPtr scp;
  CharPtr    str, new_str;

  if (action == NULL || object_list == NULL) return 0;
  scp = StringConstraintFromFieldEdit (action->edit);

  for (vnp = object_list; vnp != NULL; vnp = vnp->next) {
    str = GetFieldValueForObject (vnp->choice, vnp->data.ptrvalue, action->field, scp);
    new_str = ApplyEditToString (str, action->edit);
    if (StringCmp (str, new_str) != 0
        && SetFieldValueForObject (vnp->choice, vnp->data.ptrvalue, action->field, scp, new_str, ExistingTextOption_replace_old)) {
      num_succeed ++;
    } else {
      num_fail++;
    }
    new_str = MemFree (new_str);
    str = MemFree (str);
  }
  return num_succeed;
}


NLM_EXTERN Int4 DoConvertActionToObjectList (ConvertActionPtr action, ValNodePtr object_list, StringConstraintPtr scp)
{
  ValNodePtr vnp;
  Int4       num_succeed = 0, num_fail = 0;
  CharPtr    str, from_val;
  FieldTypePtr field_from, field_to;

  if (action == NULL || object_list == NULL || action->fields == NULL) return 0;

  field_from = GetFromFieldFromFieldPair (action->fields);
  field_to = GetToFieldFromFieldPair (action->fields);

  if (action->fields->choice == FieldPairType_molinfo_field) {
    for (vnp = object_list; vnp != NULL; vnp = vnp->next) {
      str = GetFieldValueForObject (vnp->choice, vnp->data.ptrvalue, field_from, NULL);
      from_val = GetSequenceQualValName (field_from->data.ptrvalue);
      if (StringCmp (str, from_val) == 0
          && SetFieldValueForObject (vnp->choice, vnp->data.ptrvalue, field_to, NULL, str, ExistingTextOption_replace_old)) {
        num_succeed ++;
      }
      str = MemFree (str);
    }
  } else {
    for (vnp = object_list; vnp != NULL; vnp = vnp->next) {
      str = GetFieldValueForObject (vnp->choice, vnp->data.ptrvalue, field_from, scp);
      if (SetFieldValueForObject (vnp->choice, vnp->data.ptrvalue, field_to, NULL, str, action->existing_text)
          && RemoveFieldValueForObject (vnp->choice, vnp->data.ptrvalue, field_from, scp)) {
        num_succeed ++;
      } else {
        num_fail++;
      }
      str = MemFree (str);
    }
  }

  field_from = FieldTypeFree (field_from);
  field_to = FieldTypeFree (field_to);

  return num_succeed;
}


NLM_EXTERN Int4 DoCopyActionToObjectList (CopyActionPtr action, ValNodePtr object_list, StringConstraintPtr scp)
{
  ValNodePtr vnp;
  Int4       num_succeed = 0, num_fail = 0;
  CharPtr    str;
  FieldTypePtr field_from, field_to;

  if (action == NULL || object_list == NULL) return 0;
  field_from = GetFromFieldFromFieldPair (action->fields);
  field_to = GetToFieldFromFieldPair (action->fields);

  for (vnp = object_list; vnp != NULL; vnp = vnp->next) {
    str = GetFieldValueForObject (vnp->choice, vnp->data.ptrvalue, field_from, scp);
    if (SetFieldValueForObject (vnp->choice, vnp->data.ptrvalue, field_to, NULL, str, action->existing_text)) {
      num_succeed ++;
    } else {
      num_fail++;
    }
    str = MemFree (str);
  }

  field_from = FieldTypeFree (field_from);
  field_to = FieldTypeFree (field_to);
  return num_succeed;
}


NLM_EXTERN Int4 DoSwapActionToObjectList (SwapActionPtr action, ValNodePtr object_list, StringConstraintPtr scp)
{
  ValNodePtr vnp;
  Int4       num_succeed = 0, num_fail = 0;
  CharPtr    str1, str2;
  FieldTypePtr field_from, field_to;

  if (action == NULL || object_list == NULL) return 0;
  field_from = GetFromFieldFromFieldPair (action->fields);
  field_to = GetToFieldFromFieldPair (action->fields);

  for (vnp = object_list; vnp != NULL; vnp = vnp->next) {
    str1 = GetFieldValueForObject (vnp->choice, vnp->data.ptrvalue, field_from, scp);
    str2 = GetFieldValueForObject (vnp->choice, vnp->data.ptrvalue, field_to, NULL);
    if (SetFieldValueForObject (vnp->choice, vnp->data.ptrvalue, field_to, NULL, str1, ExistingTextOption_replace_old)
        && SetFieldValueForObject (vnp->choice, vnp->data.ptrvalue, field_from, scp, str2, ExistingTextOption_replace_old)) {
      num_succeed ++;
    } else {
      num_fail++;
    }
    str1 = MemFree (str1);
    str2 = MemFree (str2);
  }
  field_from = FieldTypeFree (field_from);
  field_to = FieldTypeFree (field_to);
  return num_succeed;
}


NLM_EXTERN Int4 DoRemoveActionToObjectList (RemoveActionPtr action, ValNodePtr object_list, StringConstraintPtr scp)
{
  ValNodePtr vnp;
  Int4       num_succeed = 0, num_fail = 0;

  if (action == NULL || object_list == NULL) return 0;

  for (vnp = object_list; vnp != NULL; vnp = vnp->next) {
    if (RemoveFieldValueForObject (vnp->choice, vnp->data.ptrvalue, action->field, scp)) {
      num_succeed ++;
    } else {
      num_fail++;
    }
  }
  return num_succeed;
}


NLM_EXTERN Int4 DoParseActionToObjectList (AECRParseActionPtr action, ValNodePtr object_list, StringConstraintPtr scp)
{
  ValNodePtr vnp;
  CharPtr    str1, str2, cp;
  Int4       len, num_succeed = 0;
  FieldTypePtr field_from, field_to;

  if (action == NULL || object_list == NULL) return 0;
  field_from = GetFromFieldFromFieldPair (action->fields);
  field_to = GetToFieldFromFieldPair (action->fields);

  for (vnp = object_list; vnp != NULL; vnp = vnp->next) {
    str1 = GetFieldValueForObject (vnp->choice, vnp->data.ptrvalue, field_from, scp);
    str2 = GetTextPortionFromString (str1, action->portion);    
    if (str2 != NULL) {
      if (action->remove_from_parsed) {
        cp = StringSearch (str1, str2);
        len = StringLen (str2);
        StringCpy (cp, cp + len);
        SetFieldValueForObject (vnp->choice, vnp->data.ptrvalue, field_from, scp, str1, ExistingTextOption_replace_old);
      }
      if (SetFieldValueForObject (vnp->choice, vnp->data.ptrvalue, field_to, NULL, str1, action->existing_text)) {
        num_succeed++;
      }
    }
    str1 = MemFree (str1);
    str2 = MemFree (str2);
  }
  field_from = FieldTypeFree (field_from);
  field_to = FieldTypeFree (field_to);
  return num_succeed;
}


static Int4 ApplyAECRActionToSeqEntry (AECRActionPtr act, SeqEntryPtr sep)
{
  StringConstraintPtr scp;
  ApplyActionPtr      a;
  ValNodePtr          object_list = NULL;
  Uint1               field_type;
  Uint2               entityID;
  Int4                num_succeed = 0;

  if (act == NULL || act->action == NULL) return 0;
  field_type = FieldTypeFromAECRAction (act);
  if (field_type == FieldType_cds_gene_prot) {
    entityID = ObjMgrGetEntityIDForChoice(sep);
    object_list = BuildCGPSetList (entityID, act->constraint);
  } else {
    object_list = GetObjectListForAECRAction (sep, act);
  }

  switch (act->action->choice) {
    case ActionChoice_apply:
      a = (ApplyActionPtr) act->action->data.ptrvalue;
      scp = FindStringConstraintInConstraintSetForField (a->field, act->constraint);
      num_succeed = DoApplyActionToObjectList (act->action->data.ptrvalue, object_list, scp);
      scp = StringConstraintFree (scp);
      break;
    case ActionChoice_edit:
      num_succeed = DoEditActionToObjectList (act->action->data.ptrvalue, object_list);
      break;
    case ActionChoice_convert:
      num_succeed = DoConvertActionToObjectList (act->action->data.ptrvalue, object_list, NULL);
      break;
    case ActionChoice_swap:
      num_succeed = DoSwapActionToObjectList (act->action->data.ptrvalue, object_list, NULL);
      break;
    case ActionChoice_copy:
      num_succeed = DoCopyActionToObjectList (act->action->data.ptrvalue, object_list, NULL);
      break;
    case ActionChoice_remove:
      num_succeed = DoRemoveActionToObjectList (act->action->data.ptrvalue, object_list, NULL);
      break;
    case ActionChoice_parse:
      num_succeed = DoRemoveActionToObjectList (act->action->data.ptrvalue, object_list, NULL);
      break;
  }
  object_list = ValNodeFree (object_list);  
  return num_succeed;
}


/* This section handles parsing where the source field and destination field may not be on the same
 * group of objects. */
typedef struct parsesourceinfo 
{
  BioseqPtr   bsp;
  SeqFeatPtr  sfp;
  SeqDescrPtr sdp;
  SeqIdPtr    sip;
  ValNodePtr  dest_list;
  CharPtr     parse_src_txt;
} ParseSourceInfoData, PNTR ParseSourceInfoPtr;

static ParseSourceInfoPtr ParseSourceInfoNew (BioseqPtr bsp, SeqFeatPtr sfp, SeqDescrPtr sdp, SeqIdPtr sip, CharPtr parse_src_txt)
{
  ParseSourceInfoPtr psip;

  psip = (ParseSourceInfoPtr) MemNew (sizeof (ParseSourceInfoData));
  if (psip != NULL) {
    psip->bsp = bsp;
    psip->sdp = sdp;
    psip->sfp = sfp;
    psip->sip = sip;
    psip->dest_list = NULL;
    psip->parse_src_txt = parse_src_txt;
  } 
  return psip;
}


static ParseSourceInfoPtr ParseSourceInfoFree (ParseSourceInfoPtr psip)
{
  if (psip != NULL)
  {
    psip->dest_list = ValNodeFree (psip->dest_list);
    psip->parse_src_txt = MemFree (psip->parse_src_txt);
    psip = MemFree (psip);
  }
  return psip;
}

static ParseSourceInfoPtr ParseSourceInfoCopy (ParseSourceInfoPtr psip)
{
  ParseSourceInfoPtr pcopy = NULL;
  
  if (psip != NULL) 
  {
    pcopy = (ParseSourceInfoPtr) MemNew (sizeof (ParseSourceInfoData));
    if (pcopy != NULL) {
      pcopy->bsp = psip->bsp;
      pcopy->sfp = psip->sfp;
      pcopy->sdp = psip->sdp;
      pcopy->sip = psip->sip;
      pcopy->dest_list = NULL;
      pcopy->parse_src_txt = NULL;
    }
  }
  return pcopy;
}

static ValNodePtr ParseSourceListFree (ValNodePtr vnp)
{
  ValNodePtr vnp_next;
  while (vnp != NULL) {
    vnp_next = vnp->next;
    vnp->next = NULL;
    vnp->data.ptrvalue = ParseSourceInfoFree (vnp->data.ptrvalue);
    vnp = ValNodeFree (vnp);
    vnp = vnp_next;
  }
  return vnp;
}


static void 
GetDeflineSourcesForBioseq 
(BioseqPtr              bsp,
 TextPortionPtr         portion,
 ValNodePtr PNTR source_list)
{
  SeqDescrPtr        sdp;
  SeqMgrDescContext  dcontext;
  CharPtr            str;
  ParseSourceInfoPtr psip;
  
  if (bsp == NULL || source_list == NULL)
  {
    return;
  }
  
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_title, &dcontext);
  while (sdp != NULL)
  {
    str = GetTextPortionFromString (sdp->data.ptrvalue, portion);    
    if (str != NULL) {
      psip = ParseSourceInfoNew (bsp, NULL, sdp, NULL, str);
      if (psip != NULL) {
        ValNodeAddPointer (source_list, 0, psip);
      } else {
        str = MemFree (str);
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_title, &dcontext);
  }
}


static CharPtr GetIDSrc (SeqIdPtr sip, Uint1 id_type, CharPtr tag)
{
  DbtagPtr    dbt = NULL;
  ObjectIdPtr oip = NULL;
  Char        id_str[128];
  CharPtr     str_src = NULL;

  if (sip == NULL || sip->choice != id_type) return NULL;

  if (id_type == SEQID_GENERAL)
  {
    dbt = (DbtagPtr) sip->data.ptrvalue;
    if (dbt == NULL || (tag != NULL && StringCmp (dbt->db, tag) != 0)) return NULL;
    oip = dbt->tag;
  }
  else if (id_type == SEQID_LOCAL)
  {
    oip = sip->data.ptrvalue;
  }

  if (oip == NULL)
  {
    SeqIdWrite (sip, id_str, PRINTID_REPORT, sizeof (id_str));
    str_src = StringSave (id_str);
  }
  else
  {
    if (oip->str == NULL)
    {
      sprintf (id_str, "%d", oip->id);
      str_src = StringSave (id_str);
    }
    else
    {
      str_src = StringSave (oip->str);
    }
  }
  return str_src;
}


static void
GetIDSourcesForBioseq
(BioseqPtr       bsp,
 TextPortionPtr  portion,
 Uint1           id_type,
 CharPtr         tag,
 ValNodePtr PNTR source_list)
{
  SeqIdPtr           sip;
  ParseSourceInfoPtr psip;
  CharPtr            src_str = NULL, str;
  
  if (bsp == NULL || source_list == NULL)
  {
    return;
  }
  
  sip = bsp->id;
  while (sip != NULL)
  {
    if ((src_str = GetIDSrc (sip, id_type, tag)) != NULL) { 
      str = GetTextPortionFromString (src_str, portion); 
      if (str != NULL) {
        psip = ParseSourceInfoNew (bsp, NULL, NULL, sip, str);
        if (psip != NULL) {
          ValNodeAddPointer (source_list, 0, psip);
        } else {
          str = MemFree (str);
        }
      }
      src_str = MemFree (src_str);
    }
    sip = sip->next;
  }
}


static void
GetLocalIDSourcesForBioseq
(BioseqPtr       bsp,
 TextPortionPtr  tp,
 ValNodePtr PNTR source_list)
{
  GetIDSourcesForBioseq (bsp, tp, SEQID_LOCAL, NULL, source_list);
}


static void GetNcbiFileSourceForBioseq
(BioseqPtr       bsp,
 TextPortionPtr  tp,
 ValNodePtr PNTR source_list)
{
  GetIDSourcesForBioseq (bsp, tp, SEQID_GENERAL, "NCBIFILE", source_list);
}


static void StripBankitCommentForParse (SeqDescrPtr sdp, TextPortionPtr tp)
{
  UserObjectPtr      uop;
  ObjectIdPtr        oip;
  UserFieldPtr       ufp;
  
  if (sdp == NULL || sdp->choice != Seq_descr_user || tp == NULL) {
    return;
  }
  
  /* Bankit Comments */
  uop = (UserObjectPtr) sdp->data.ptrvalue;
  if (uop != NULL && StringCmp (uop->_class, "SMART_V1.0") != 0) {
    oip = uop->type;
    if (oip != NULL && StringCmp (oip->str, "Submission") == 0) {
      for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
        oip = ufp->label;
        if (oip != NULL && StringCmp (oip->str, "AdditionalComment") == 0) {
          ReplaceStringForParse (ufp->data.ptrvalue, tp);
        }
      }
    }
  }
}


static void StripStructuredCommentForParse (SeqDescrPtr sdp, CharPtr comment_field, TextPortionPtr tp)
{
  UserObjectPtr      uop;
  ObjectIdPtr        oip;
  UserFieldPtr       ufp;

  if (sdp == NULL || sdp->choice != Seq_descr_user || tp == NULL || StringHasNoText (comment_field)) {
    return;
  }
    
  uop = (UserObjectPtr) sdp->data.ptrvalue;
  oip = uop->type;
  if (oip != NULL && StringCmp (oip->str, "StructuredComment") == 0) {
    for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
      oip = ufp->label;
      if (oip != NULL && StringCmp (oip->str, comment_field) == 0) {
        ReplaceStringForParse (ufp->data.ptrvalue, tp);
      }
    }
  }
}


static void
GetBankitCommentSourcesForBioseq 
(BioseqPtr       bsp,
 TextPortionPtr  tp,
 ValNodePtr PNTR source_list)
{
  SeqDescrPtr        sdp;
  SeqMgrDescContext  dcontext;
  ParseSourceInfoPtr psip;
  UserObjectPtr      uop;
  ObjectIdPtr        oip;
  UserFieldPtr       ufp;
  CharPtr            str = NULL;
  
  if (bsp == NULL || source_list == NULL) {
    return;
  }
  
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &dcontext);
  while (sdp != NULL) {
    if (sdp->extended != 0) {
      /* Bankit Comments */
      uop = (UserObjectPtr) sdp->data.ptrvalue;
      if (uop != NULL && StringCmp (uop->_class, "SMART_V1.0") != 0) {
        oip = uop->type;
        if (oip != NULL && StringCmp (oip->str, "Submission") == 0) {
          for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
            oip = ufp->label;
            if (oip != NULL && StringCmp (oip->str, "AdditionalComment") == 0) {
              str = GetTextPortionFromString (ufp->data.ptrvalue, tp);
              if (str != NULL) {
                psip = ParseSourceInfoNew (bsp, NULL, sdp, NULL, str);
                if (psip == NULL) {
                  str = MemFree (str);
                } else {
                  ValNodeAddPointer (source_list, 0, psip);
                }
              }
            }
          }
        }
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &dcontext);
  }
}


static void 
GetCommentSourcesForBioseq 
(BioseqPtr       bsp,
 TextPortionPtr  tp,
 ValNodePtr PNTR source_list)
{
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp;
  SeqMgrFeatContext  fcontext;
  SeqMgrDescContext  dcontext;
  ParseSourceInfoPtr psip;
  CharPtr            str;
  
  if (bsp == NULL || source_list == NULL) {
    return;
  }
  
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_comment, &dcontext);
  while (sdp != NULL) {
    str = GetTextPortionFromString (sdp->data.ptrvalue, tp);
    if (str != NULL) {
      psip = ParseSourceInfoNew (bsp, NULL, sdp, NULL, str);
      if (psip == NULL) {
        str = MemFree (str);
      } else {
        ValNodeAddPointer (source_list, 0, psip);
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_comment, &dcontext);
  }
  
  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_COMMENT, 0, &fcontext);
  while (sfp != NULL) {
    str = GetTextPortionFromString (sfp->data.value.ptrvalue, tp);
    if (str != NULL) {
      psip = ParseSourceInfoNew (bsp, sfp, NULL, NULL, str);
      if (psip == NULL) {
        str = MemFree (str);
      } else {
        ValNodeAddPointer (source_list, 0, psip);
      }
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_COMMENT, 0, &fcontext);
  }
  GetBankitCommentSourcesForBioseq (bsp, tp, source_list);
}


static void 
GetStructuredCommentSourcesForBioseq 
(BioseqPtr       bsp,
 TextPortionPtr  tp,
 CharPtr         comment_field,
 ValNodePtr PNTR source_list)
{
  SeqDescrPtr        sdp;
  UserObjectPtr      uop;
  ObjectIdPtr        oip;
  UserFieldPtr       ufp;
  SeqMgrDescContext  dcontext;
  CharPtr            str;
  ParseSourceInfoPtr psip;
  
  if (bsp == NULL || source_list == NULL)
  {
    return;
  }
  
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &dcontext);
  while (sdp != NULL) {  
    if (sdp->extended != 0
        && sdp->data.ptrvalue != NULL) {
      uop = (UserObjectPtr) sdp->data.ptrvalue;
      oip = uop->type;
      if (oip != NULL && StringCmp (oip->str, "StructuredComment") == 0) {
        for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
          oip = ufp->label;
          if (oip != NULL && StringCmp (oip->str, comment_field) == 0) {
            str = GetTextPortionFromString (ufp->data.ptrvalue, tp);
            if (str != NULL) {
              psip = ParseSourceInfoNew (bsp, NULL, sdp, NULL, str);
              if (psip == NULL) {
                str = MemFree (str);
              } else {
                ValNodeAddPointer (source_list, 0, psip);
              }
            }
          }
        }
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &dcontext);
  }
}


const CharPtr nomial_keywords[] = {
"f. sp. ",
"var.",
"pv.",
"bv.",
"serovar",
"subsp." };

const Int4 num_nomial_keywords = sizeof(nomial_keywords) / sizeof (CharPtr);

static CharPtr GetTextAfterNomial (CharPtr taxname)

{
  CharPtr ptr, nomial_end;
  Int4    i;
  Boolean found_keyword = TRUE;
  
  ptr = StringChr (taxname, ' ');
  if (ptr == NULL) return NULL;
  /* skip over the first word and the spaces after it. */
  while (*ptr == ' ') {
    ptr++;
  }
  ptr = StringChr (ptr, ' ');
  /* if there are only two words, give up. */
  if (ptr == NULL) {
    return NULL;
  }
  nomial_end = ptr;
  while (*ptr == ' ') {
    ptr++;
  }
  
  while (found_keyword) {
    found_keyword = FALSE;
    /* if the next word is a nomial keyword, skip that plus the first word that follows it. */
    for (i = 0; i < num_nomial_keywords && *nomial_end != 0; i++) {
      if (StringNCmp (ptr, nomial_keywords[i], StringLen(nomial_keywords[i])) == 0) {
        ptr += StringLen(nomial_keywords[i]);
        while (*ptr == ' ' ) {
          ptr++;
        }
        nomial_end = StringChr (ptr, ' ');
        if (nomial_end == NULL) {
          nomial_end = ptr + StringLen (ptr);
        } else {          
          ptr = nomial_end;
          while (*ptr == ' ') {
            ptr++;
          }
          found_keyword = TRUE;
        }
      }
    }
  }
  return nomial_end;
}


static void 
GetOrgParseSourcesForBioSource 
(BioSourcePtr    biop,
 BioseqPtr       bsp,
 SeqDescrPtr     sdp,
 SeqFeatPtr      sfp,
 ParseSrcOrgPtr  o,
 TextPortionPtr  tp,
 ValNodePtr PNTR source_list)
{
  CharPtr str = NULL, portion, tmp;
  ValNode vn;
  ParseSourceInfoPtr psip;

  if (biop == NULL || o == NULL || o->field == NULL || source_list == NULL) return;

  switch (o->field->choice) {
    case ParseSrcOrgChoice_source_qual :
      vn.choice = SourceQualChoice_textqual;
      vn.data.intvalue = o->field->data.intvalue;
      vn.next = NULL;
      str = GetSourceQualFromBioSource (biop, &vn, NULL);
      break;
    case ParseSrcOrgChoice_taxname_after_binomial :
      vn.choice = SourceQualChoice_textqual;
      vn.data.intvalue = Source_qual_taxname;
      vn.next = NULL;
      str = GetSourceQualFromBioSource (biop, &vn, NULL);
      tmp = GetTextAfterNomial (str);
      tmp = StringSave (tmp);
      str = MemFree (str);
      str = tmp;
      break;
  }
  portion = GetTextPortionFromString (str, tp);
  if (portion != NULL) {
    psip = ParseSourceInfoNew (bsp, sfp, sdp, NULL, portion);
    if (psip == NULL) {
      portion = MemFree (portion);
    } else {
      ValNodeAddPointer (source_list, 0, psip);
    }
  }
  str = MemFree (str);
}


static void GetOrgParseSourcesForBioseq (BioseqPtr bsp, ParseSrcOrgPtr o, TextPortionPtr tp, ValNodePtr PNTR source_list)
{
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp;
  SeqMgrFeatContext  fcontext;
  SeqMgrDescContext  dcontext;

  if (bsp == NULL || o == NULL || source_list == NULL) return;

  if (o->type == Object_type_constraint_any || o->type == Object_type_constraint_descriptor) {
    for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
         sdp != NULL;
         sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_source, &dcontext)) {
      GetOrgParseSourcesForBioSource (sdp->data.ptrvalue, bsp, sdp, NULL, o, tp, source_list);
    }
  }

  if (o->type == Object_type_constraint_any || o->type == Object_type_constraint_feature) {
    for (sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_BIOSRC, 0, &fcontext);
         sfp != NULL;
         sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_BIOSRC, 0, &fcontext)) {
      GetOrgParseSourcesForBioSource (sfp->data.value.ptrvalue, bsp, NULL, sfp, o, tp, source_list);
    }
  }
}


typedef struct parsesrccollection {
  ParseSrcPtr src;
  TextPortionPtr portion;
  ValNodePtr src_list;
} ParseSrcCollectionData, PNTR ParseSrcCollectionPtr;


static void FindParseSourceBioseqCallback (BioseqPtr bsp, Pointer userdata)
{
  ParseSrcCollectionPtr psp;
  
  if (bsp == NULL || userdata == NULL)
  {
    return;
  }
  
  psp = (ParseSrcCollectionPtr) userdata;
  if (psp->src == NULL) return;

  switch (psp->src->choice)
  {
    case ParseSrc_defline:
      if (!ISA_aa (bsp->mol)) {
        GetDeflineSourcesForBioseq (bsp, psp->portion, &(psp->src_list));
      }
      break;
    case ParseSrc_local_id:
      if (! ISA_aa (bsp->mol) && bsp->repr != Seq_repr_seg) {
        GetLocalIDSourcesForBioseq (bsp, psp->portion, &(psp->src_list));
      }
      break;
    case ParseSrc_file_id:
      GetNcbiFileSourceForBioseq (bsp, psp->portion, &(psp->src_list));
      break;
    case ParseSrc_org:
      GetOrgParseSourcesForBioseq (bsp, psp->src->data.ptrvalue, psp->portion, &(psp->src_list));
      break;
    case ParseSrc_comment:
      GetCommentSourcesForBioseq (bsp, psp->portion, &(psp->src_list));
      break;
    case ParseSrc_structured_comment:
      GetStructuredCommentSourcesForBioseq(bsp, psp->portion, psp->src->data.ptrvalue, &(psp->src_list));
      break;
    case ParseSrc_bankit_comment:
      if (!ISA_aa (bsp->mol)) {
        GetBankitCommentSourcesForBioseq (bsp, psp->portion, &(psp->src_list));
      }
      break;
  }
}


static void GetOrgNamesInRecordCallback (BioSourcePtr biop, Pointer userdata)
{
  ValNodePtr PNTR org_names;
  
  if (biop == NULL || biop->org == NULL || StringHasNoText (biop->org->taxname)
      || userdata == NULL)
  {
    return;
  }
  
  org_names = (ValNodePtr PNTR) userdata;
  
  ValNodeAddPointer (org_names, 0, biop->org->taxname);
}


static void SetToUpper (CharPtr cp)
{
  if (cp == NULL) return;
  while (*cp != 0) {
    if (isalpha (*cp)) {
      *cp = toupper (*cp);
    }
    cp++;
  }
}


static void 
FixCapitalizationInString 
(CharPtr PNTR pTitle,
 Uint2 capitalization,
 ValNodePtr   org_names)
{
  if (pTitle == NULL || capitalization == Cap_change_none) return;

  switch (capitalization) {
    case Cap_change_tolower:
      ResetCapitalization (FALSE, *pTitle);
      FixAbbreviationsInElement (pTitle);
      FixOrgNamesInString (*pTitle, org_names);
      break;
    case Cap_change_toupper:
      SetToUpper (*pTitle);
      FixAbbreviationsInElement (pTitle);
      FixOrgNamesInString (*pTitle, org_names);
      break;
    case Cap_change_firstcap:
      ResetCapitalization (TRUE, *pTitle);
      FixAbbreviationsInElement (pTitle);
      FixOrgNamesInString (*pTitle, org_names);
      break;
  }
}


static void AddDeflineDestinationsForBioseq (BioseqPtr bsp, ValNodePtr PNTR dest_list)
{
  SeqDescrPtr        sdp;
  SeqMgrDescContext  dcontext;

  if (bsp == NULL || dest_list == NULL) {
    return;
  }
  
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_title, &dcontext);
  while (sdp != NULL) {
    ValNodeAddPointer (dest_list, OBJ_SEQDESC, sdp);
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_title, &dcontext);
  }
}


static void AddFeatureDestinationsForBioseq (BioseqPtr bsp, FeatureFieldLegalPtr featfield, ValNodePtr PNTR dest_list)
{
  SeqFeatPtr        sfp;
  SeqMgrFeatContext fcontext;
  Int4             featdef;

  if (bsp == NULL || featfield == NULL || dest_list == NULL) return;

  featdef = GetFeatdefFromFeatureType (featfield->type);
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, featdef, &fcontext);
       sfp != NULL;
       sfp = SeqMgrGetNextFeature (bsp, sfp, 0, featdef, &fcontext)) {
    ValNodeAddPointer (dest_list, OBJ_SEQFEAT, sfp);
  }
}


static void GetBioSourceDestinationsForBioseq (BioseqPtr bsp, Uint2 object_type, ValNodePtr PNTR dest_list)
{
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp;
  SeqMgrFeatContext  fcontext;
  SeqMgrDescContext  dcontext;

  if (bsp == NULL || dest_list == NULL)
  {
    return;
  }
  
  if (object_type == Object_type_constraint_any || object_type == Object_type_constraint_descriptor) 
  {
    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
    while (sdp != NULL)
    {
      ValNodeAddPointer (dest_list, OBJ_SEQDESC, sdp);
      sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_source, &dcontext);
    }
  }
  
  if (object_type == Object_type_constraint_any || object_type == Object_type_constraint_feature)
  {
    sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_BIOSRC, 0, &fcontext);
    while (sfp != NULL)
    {
      ValNodeAddPointer (dest_list, OBJ_SEQFEAT, sfp);
      sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_BIOSRC, 0, &fcontext);
    }  
  }
}


static void AddParseDestinations (ParseSourceInfoPtr psip, ParseDestPtr dst)
{
  ParseDstOrgPtr o;

  if (psip == NULL || dst == NULL) return;

  switch (dst->choice) {
    case ParseDest_defline :
      AddDeflineDestinationsForBioseq (psip->bsp, &(psip->dest_list));
      break;
    case ParseDest_org :
      o = (ParseDstOrgPtr) dst->data.ptrvalue;
      if ((o->type == Object_type_constraint_any || o->type == Object_type_constraint_descriptor)
          && psip->sdp != NULL && psip->sdp->choice == Seq_descr_source) {
        ValNodeAddPointer (&(psip->dest_list), OBJ_SEQDESC, psip->sdp);
      } else if ((o->type == Object_type_constraint_any || o->type == Object_type_constraint_feature)
                 && psip->sfp != NULL && psip->sfp->data.choice == SEQFEAT_BIOSRC) {
        ValNodeAddPointer (&(psip->dest_list), OBJ_SEQFEAT, psip->sfp);
      } else {
        GetBioSourceDestinationsForBioseq (psip->bsp, o->type, &(psip->dest_list));
      }
      break;
    case ParseDest_featqual :
      AddFeatureDestinationsForBioseq (psip->bsp, dst->data.ptrvalue, &(psip->dest_list));
      break;
    case ParseDest_dbxref :
      GetBioSourceDestinationsForBioseq (psip->bsp, Object_type_constraint_any, &(psip->dest_list));
      break;
  }
}


static Boolean SourceHasOneUndeletedDestination (ParseSourceInfoPtr source)
{
  Int4       num_seen = 0;
  ValNodePtr vnp;
  
  if (source == NULL
      || source->dest_list == NULL)
  {
    return FALSE;
  }
  
  vnp = source->dest_list;
  while (vnp != NULL && num_seen < 2)
  {
    if (vnp->choice > 1)
    {
      num_seen ++;
    }
    vnp = vnp->next;
  }
  if (num_seen == 1)
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}


static void CombineSourcesForDestinations (ValNodePtr PNTR source_list)
{
  ValNodePtr         source1_vnp, source2_vnp, dest1_vnp, dest2_vnp;
  ValNodePtr         source_new, del_vnp;
  ParseSourceInfoPtr psip1, psip2, new_psip;
  CharPtr            comb_txt;
  
  for (source1_vnp = *source_list;
       source1_vnp != NULL; 
       source1_vnp = source1_vnp->next)
  {
    psip1 = (ParseSourceInfoPtr) source1_vnp->data.ptrvalue;
    if (psip1 == NULL || psip1->dest_list == NULL)
    {
      continue;
    }
    for (source2_vnp = source1_vnp->next;
         source2_vnp != NULL; 
         source2_vnp = source2_vnp->next)
    {
      if (source2_vnp->choice > 0) 
      {
        /* already marked for deletion */
        continue;
      }
      psip2 = (ParseSourceInfoPtr) source2_vnp->data.ptrvalue;
      if (psip2 == NULL || psip2->dest_list == NULL)
      {
        continue;
      }
      for (dest1_vnp = psip1->dest_list;
           dest1_vnp != NULL; 
           dest1_vnp = dest1_vnp->next)
      {
        if (dest1_vnp->choice == 0)
        {
          /* already marked for deletion */
          continue;
        }
        for (dest2_vnp = psip2->dest_list;
             dest2_vnp != NULL;
             dest2_vnp = dest2_vnp->next)
        {
          if (dest2_vnp->choice == 0)
          {
            /* already marked for deletion */
            continue;
          }
          if (dest1_vnp->choice == dest2_vnp->choice
              && dest1_vnp->data.ptrvalue == dest2_vnp->data.ptrvalue)
          {
            comb_txt = (CharPtr) (MemNew (sizeof (Char) 
                                  * (StringLen (psip1->parse_src_txt)
                                     + StringLen (psip2->parse_src_txt)
                                     + 2)));
            StringCpy (comb_txt, psip1->parse_src_txt);
            StringCat (comb_txt, ";");
            StringCat (comb_txt, psip2->parse_src_txt);
            
            /* If the first source has a single destination, then we can 
             * add the text from the second source to the first and remove
             * the destination from the second source.
             */
            if (SourceHasOneUndeletedDestination (psip1))
            {
              
              psip1->parse_src_txt = MemFree (psip1->parse_src_txt);
              psip1->parse_src_txt = comb_txt;
              dest2_vnp->choice = 0;
            }             
            /* If the first source has more than one destination and
             * the second source has a single destination, then we can 
             * remove the repeated desination from the first source
             * and add the text from the first source to the second source.
             */
            else if (SourceHasOneUndeletedDestination (psip2))
            {
              psip2->parse_src_txt = MemFree (psip2->parse_src_txt);
              psip2->parse_src_txt = comb_txt;
              dest1_vnp->choice = 0;
            }
            /* If the first and second sources have multiple destinations,
             * we need to remove the repeated destination from both the first
             * and second source and create a new source with the combined 
             * text for just the repeated destination.
             */
            else
            {
              new_psip = ParseSourceInfoNew (NULL, NULL, NULL, NULL, comb_txt);
              ValNodeAddPointer (&(new_psip->dest_list), 
                                 dest1_vnp->choice, 
                                 dest1_vnp->data.ptrvalue);
              dest1_vnp->choice = 0;
              dest2_vnp->choice = 0;
              source_new = ValNodeNew (NULL);
              source_new->choice = 0;
              source_new->data.ptrvalue = new_psip;
              source_new->next = source1_vnp->next;
              source1_vnp->next = source_new;
            }
          }
        }
      }
      
      del_vnp = ValNodeExtractList (&(psip1->dest_list), 0);
      del_vnp = ValNodeFree (del_vnp);
      if (psip1->dest_list == NULL)
      {
        source1_vnp->choice = 1;
      }
      del_vnp = ValNodeExtractList (&(psip2->dest_list), 0);
      del_vnp = ValNodeFree (del_vnp);
      if (psip2->dest_list == NULL)
      {
        source2_vnp->choice = 1;
      }
    }
  }

  /* now remove sources deleted */
  del_vnp = ValNodeExtractList (source_list, 1);
  del_vnp = ParseSourceListFree (del_vnp); 
}


static BioseqSetPtr GetPartsForSourceDescriptorOnSegSet (SeqDescrPtr sdp)
{
  ObjValNodePtr ovp;
  BioseqSetPtr  bssp;
  SeqEntryPtr   sep;
  
  if (sdp == NULL || sdp->extended != 1) {
    return NULL;
  }
  ovp = (ObjValNodePtr) sdp;
  if (ovp->idx.parenttype != OBJ_BIOSEQSET || ovp->idx.parentptr == NULL) {
    return NULL;
  }
  bssp = (BioseqSetPtr) ovp->idx.parentptr;
  
  if (bssp->_class == BioseqseqSet_class_nuc_prot
      && IS_Bioseq_set (bssp->seq_set)
      && bssp->seq_set->data.ptrvalue != NULL) {
    bssp = (BioseqSetPtr) bssp->seq_set->data.ptrvalue;
  }
  
  if (bssp->_class == BioseqseqSet_class_segset) {
    sep = bssp->seq_set;
    while (sep != NULL) {
      if (IS_Bioseq_set (sep) && sep->data.ptrvalue != NULL) {
        bssp = (BioseqSetPtr) sep->data.ptrvalue;
        if (bssp->_class == BioseqseqSet_class_parts) {
          return bssp;
        }
      }
      sep = sep->next;
    }
  }

  return NULL;
}


static SeqDescrPtr FindSourceDescriptorInSeqEntry (SeqEntryPtr sep)
{
  BioseqPtr    bsp;
  BioseqSetPtr bssp;
  SeqDescrPtr  sdp = NULL;
  
  if (sep != NULL && sep->data.ptrvalue != NULL) {
    if (IS_Bioseq (sep)) {
      bsp = (BioseqPtr) sep->data.ptrvalue;
      sdp = bsp->descr;
    } else if (IS_Bioseq_set (sep)) {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      sdp = bssp->descr;
    }
    while (sdp != NULL && sdp->choice != Seq_descr_source)
    {
      sdp = sdp->next;
    }
  }
  return sdp;
}


static SeqDescrPtr PropagateToSeqEntry (SeqEntryPtr sep, SeqDescrPtr sdp)
{
  BioseqPtr    bsp;
  BioseqSetPtr bssp;
  SeqDescrPtr  new_sdp = NULL;
  
  if (sep != NULL && sep->data.ptrvalue != NULL) {
    if (IS_Bioseq (sep)) {
      bsp = (BioseqPtr) sep->data.ptrvalue;
      new_sdp = AsnIoMemCopy ((Pointer) sdp,
                              (AsnReadFunc) SeqDescrAsnRead,
                              (AsnWriteFunc) SeqDescrAsnWrite);
      ValNodeLink (&(bsp->descr), new_sdp);
    } else if (IS_Bioseq_set (sep)) {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      new_sdp = AsnIoMemCopy ((Pointer) sdp,
                              (AsnReadFunc) SeqDescrAsnRead,
                              (AsnWriteFunc) SeqDescrAsnWrite);
      ValNodeLink (&(bssp->descr), new_sdp);
    }
  }
  return new_sdp;
}


static void PropagateSourceOnSegSetForParse (ValNodePtr parse_source_list)
{
  ParseSourceInfoPtr psip;
  ValNodePtr         vnp_src, vnp_dst;
  SeqDescrPtr        sdp, other_sdp;
  SeqEntryPtr        sep;
  ValNodePtr         extra_dests = NULL;
  BioseqSetPtr       parts_bssp;
  
  for (vnp_src = parse_source_list; vnp_src != NULL; vnp_src = vnp_src->next) {
    psip = (ParseSourceInfoPtr) vnp_src->data.ptrvalue;
    if (psip != NULL) {
      for (vnp_dst = psip->dest_list; vnp_dst != NULL; vnp_dst = vnp_dst->next) {
        if (vnp_dst->choice == OBJ_SEQDESC) {
          sdp = (SeqDescrPtr) vnp_dst->data.ptrvalue;
          if (sdp != NULL && sdp->choice == Seq_descr_source) {
            parts_bssp = GetPartsForSourceDescriptorOnSegSet (sdp);
            if (parts_bssp != NULL) {
              for (sep = parts_bssp->seq_set; sep != NULL; sep = sep->next) {
                if (IS_Bioseq(sep) && sep->data.ptrvalue == psip->bsp) {
                  other_sdp = FindSourceDescriptorInSeqEntry (sep);
                  if (other_sdp == NULL) {
                    other_sdp = PropagateToSeqEntry (sep, sdp);
                    ValNodeAddPointer (&extra_dests, OBJ_SEQDESC, other_sdp);
                  }
                }
              }
            
              /* set choice to 0 so master won't be a destination */
              vnp_dst->choice = 0;
            
            }
          }
        }
      }
      /* add extra destinations to list */
      ValNodeLink (&psip->dest_list, extra_dests);
      extra_dests = NULL;
    }
  }
  
}


static Boolean SetDBxrefForBioSource (BioSourcePtr biop, CharPtr db_name, CharPtr str, Uint2 existing_text)
{
  ValNodePtr    dbx;
  DbtagPtr      dbtag;
  Boolean       found = FALSE;
  Char          buf[20];
  Boolean       rval = FALSE;

  if (biop == NULL || StringHasNoText (db_name) || StringHasNoText (str)) {
    return FALSE;
  }

  if (biop->org == NULL)
  {
    biop->org = OrgRefNew();
  }
  dbx = biop->org->db;
  while (dbx != NULL && !found)
  {
    dbtag = (DbtagPtr) dbx->data.ptrvalue;
    if (dbtag != NULL && dbtag->tag != NULL
        && StringCmp (dbtag->db, db_name) == 0)
    {
      found = TRUE;
    }
    if (!found)
    {
      dbx = dbx->next;
    }
  }
  if (!found)
  {
    dbtag = DbtagNew();
    dbtag->db = StringSave (db_name);      
    ValNodeAddPointer (&(biop->org->db), 0, dbtag);
  }
  if (dbtag->tag == NULL)
  {
    dbtag->tag = ObjectIdNew();
  }
  /* if it was a number before, make it a string now */
  if (dbtag->tag->id > 0 && dbtag->tag->str == NULL)
  {
    sprintf (buf, "%s", dbtag->tag->id);
    dbtag->tag->id = 0;
    dbtag->tag->str = StringSave (buf);
  }
  rval = SetStringValue (&(dbtag->tag->str), str, existing_text);
  return rval;
}


static Int4 SetFieldForDestList (ValNodePtr dest_list, ParseDestPtr field, CharPtr str, Uint2 existing_text)
{
  ValNodePtr vnp;
  SeqDescrPtr sdp;
  CharPtr     cp;
  BioSourcePtr biop;
  ParseDstOrgPtr o;
  FeatureFieldLegalPtr fl;
  FeatureField f;
  Int4         num_succeeded = 0;

  if (dest_list == NULL || field == NULL) return 0;

  switch (field->choice) {
    case ParseDest_defline :
      for (vnp = dest_list; vnp != NULL; vnp = vnp->next) {
        if (vnp->choice == OBJ_SEQDESC && vnp->data.ptrvalue != NULL) {
          sdp = (SeqDescrPtr) vnp->data.ptrvalue;
          if (sdp->choice == Seq_descr_title) {
            cp = sdp->data.ptrvalue;
            if (SetStringValue (&cp, str, existing_text)) {
              num_succeeded++;
            }
            sdp->data.ptrvalue = cp;
          }
        }
      }
      break;
    case ParseDest_org :
      o = (ParseDstOrgPtr) field->data.ptrvalue;
      if (o != NULL) {
        for (vnp = dest_list; vnp != NULL; vnp = vnp->next) {
          biop = GetBioSourceFromObject (vnp->choice, vnp->data.ptrvalue);
          if (SetSourceQualInBioSource (biop, o->field, NULL, str, existing_text)) {
            num_succeeded++;
          }
        }
      }
      break;
    case ParseDest_featqual:
      fl = (FeatureFieldLegalPtr) field->data.ptrvalue;
      if (fl != NULL) {
        f.type = fl->type;
        f.field = ValNodeNew(NULL);
        f.field->next = NULL;
        f.field->choice = FeatQualChoice_legal_qual;
        f.field->data.intvalue = fl->field;        
        for (vnp = dest_list; vnp != NULL; vnp = vnp->next) {
          if (SetQualOnFeature (vnp->data.ptrvalue, &f, NULL, str, existing_text)) {
            num_succeeded++;
          }
        }
        f.field = ValNodeFree (f.field);
      }
      break;
    case ParseDest_dbxref:
      if (!StringHasNoText (field->data.ptrvalue)) {
        for (vnp = dest_list; vnp != NULL; vnp = vnp->next) {
          biop = GetBioSourceFromObject (vnp->choice, vnp->data.ptrvalue);
          if (SetDBxrefForBioSource (biop, field->data.ptrvalue, str, existing_text)) {
            num_succeeded++;
          }
        }
      }
      break;
  }
  return num_succeeded;
}


static void StripFieldForSrcList (ParseSourceInfoPtr psip, ParseSrcPtr field, TextPortionPtr text_portion)
{
  CharPtr     str;
  ParseSrcOrgPtr o;
  BioSourcePtr biop;

  if (psip == NULL || field == NULL || text_portion == NULL) return;

  switch (field->choice) {
    case ParseSrc_defline :
      if (psip->sdp != NULL && psip->sdp->choice == Seq_descr_title) {
        ReplaceStringForParse (psip->sdp->data.ptrvalue, text_portion);
      }
      break;
    case ParseSrc_org :
      o = (ParseSrcOrgPtr) field->data.ptrvalue;
      if (o != NULL) {
        if (psip->sdp != NULL && psip->sdp->choice == Seq_descr_source) {
          biop = (BioSourcePtr) psip->sdp->data.ptrvalue;
          str = GetSourceQualFromBioSource (biop, o->field, NULL);
          ReplaceStringForParse (str, text_portion);
          SetSourceQualInBioSource (biop, o->field, NULL, str, ExistingTextOption_replace_old);
          str = MemFree (str);
        } else if (psip->sfp != NULL && psip->sfp->data.choice == SEQFEAT_BIOSRC) {
          biop = (BioSourcePtr) psip->sfp->data.value.ptrvalue;
          str = GetSourceQualFromBioSource (biop, o->field, NULL);
          ReplaceStringForParse (str, text_portion);
          SetSourceQualInBioSource (biop, o->field, NULL, str, ExistingTextOption_replace_old);
          str = MemFree (str);
        }
      }
      break;
    case ParseSrc_comment:
      if (psip->sdp != NULL) {
        if (psip->sdp->choice == Seq_descr_user) {
          StripBankitCommentForParse (psip->sdp, text_portion);
        } else if (psip->sdp->choice == Seq_descr_comment) {
          ReplaceStringForParse (psip->sdp->data.ptrvalue, text_portion);
        }
      }
      if (psip->sfp != NULL && psip->sfp->data.choice == SEQFEAT_COMMENT) {
        ReplaceStringForParse (psip->sfp->data.value.ptrvalue, text_portion);
      }
      break;
    case ParseSrc_bankit_comment:
      if (psip->sdp != NULL && psip->sdp->choice == Seq_descr_user) {
        StripBankitCommentForParse (psip->sdp, text_portion);
      }
      break;
    case ParseSrc_structured_comment:
      if (psip->sdp != NULL && psip->sdp->choice == Seq_descr_user) {
        StripStructuredCommentForParse (psip->sdp, field->data.ptrvalue, text_portion);
      }
      break;
  }
}


static Int4 ApplyParseActionToSeqEntry (ParseActionPtr action, SeqEntryPtr sep)
{
  ParseSrcCollectionData psd;
  ParseSourceInfoPtr     psip;
  ValNodePtr             orgnames = NULL, source_list_for_removal = NULL, vnp;
  Int4                   num_succeeded = 0;

  if (action == NULL || sep == NULL) return 0;

  psd.src = action->src;
  psd.portion = action->portion;
  psd.src_list = NULL;

  /* first, we need to get a list of the parse sources */  
  VisitBioseqsInSep (sep, &psd, FindParseSourceBioseqCallback);

  if (action->capitalization != Cap_change_none) {
    /* if we will be fixing capitalization, get org names to use in fixes */
    VisitBioSourcesInSep (sep, &orgnames, GetOrgNamesInRecordCallback);
  }

  /* for each parse source, we need to get a list of the destinations */
  for (vnp = psd.src_list; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->data.ptrvalue == NULL) continue;
    psip = (ParseSourceInfoPtr) vnp->data.ptrvalue;
    if (action->remove_from_parsed) {
        ValNodeAddPointer (&source_list_for_removal, 0, ParseSourceInfoCopy (psip));
    }
    /* fix source text */
    FixCapitalizationInString (&(psip->parse_src_txt), action->capitalization, orgnames);

    /* find destinations */
    AddParseDestinations (psip, action->dest);

  }

  /* free orgname list if we created it */
  orgnames = ValNodeFree (orgnames);

  CombineSourcesForDestinations (&(psd.src_list));

  if (action->dest->choice == ParseDest_org) {
    PropagateSourceOnSegSetForParse (psd.src_list);
  }
  
  /* now do the parsing */
  for (vnp = psd.src_list; vnp != NULL; vnp = vnp->next) {
    psip = (ParseSourceInfoPtr) vnp->data.ptrvalue;
    num_succeeded += SetFieldForDestList (psip->dest_list, action->dest, psip->parse_src_txt, action->existing_text);
  }

  /* now remove strings from sources */
  for (vnp = source_list_for_removal; vnp != NULL; vnp = vnp->next)
  {
    if (vnp->data.ptrvalue == NULL) continue;
    psip = (ParseSourceInfoPtr) vnp->data.ptrvalue;
    StripFieldForSrcList (psip, action->src, action->portion);
  }
  return num_succeeded;
}


static void SetCdRegionGeneticCode (SeqFeatPtr cds)
{
  CdRegionPtr crp;
  SeqEntryPtr parent_sep;
  BioseqPtr   bsp;
  Int4        genCode;
  ValNodePtr  code, vnp;

  if (cds == NULL || cds->data.choice != SEQFEAT_CDREGION) return;
  if (cds->data.value.ptrvalue == NULL) {
    cds->data.value.ptrvalue = CdRegionNew();
  }
  crp = (CdRegionPtr) cds->data.value.ptrvalue;
  bsp = BioseqFindFromSeqLoc (cds->location);
  if (bsp == NULL) return;
  parent_sep = GetBestTopParentForData (bsp->idx.entityID, bsp);
  genCode = SeqEntryToGeneticCode (parent_sep, NULL, NULL, 0);

  code = ValNodeNew (NULL);
  if (code != NULL) {
    code->choice = 254;
    vnp = ValNodeNew (NULL);
    code->data.ptrvalue = vnp;
    if (vnp != NULL) {
      vnp->choice = 2;
      vnp->data.intvalue = genCode;
    }
  }
  crp->genetic_code = code;
}

  
static void CreateDataForFeature (SeqFeatPtr sfp, Int4 feature_type)
{
  Int4 featdef, seqfeattype;
  CharPtr    label = NULL;
  RnaRefPtr  rrp;
  GBQualPtr  gbq;
  ImpFeatPtr ifp;

  featdef = GetFeatdefFromFeatureType (feature_type);
  sfp->idx.subtype = featdef;
  seqfeattype = FindFeatFromFeatDefType (featdef);
  switch (seqfeattype) {
    case SEQFEAT_GENE:
      sfp->data.value.ptrvalue = GeneRefNew();
      break;
    case SEQFEAT_CDREGION:
      sfp->data.value.ptrvalue = CdRegionNew();
      SetCdRegionGeneticCode (sfp);
      break;
    case SEQFEAT_RNA:
      rrp = RnaRefNew();
      rrp->ext.choice = 0;
      sfp->data.value.ptrvalue = rrp;
      switch (featdef) {
        case FEATDEF_preRNA:
          rrp->type = RNA_TYPE_premsg;
          break;
        case FEATDEF_mRNA:
          rrp->type = RNA_TYPE_mRNA;
          break;
        case FEATDEF_tRNA:
          rrp->type = RNA_TYPE_tRNA;
          break;
        case FEATDEF_rRNA:
          rrp->type = RNA_TYPE_rRNA;
          break;
        case FEATDEF_snRNA:
          rrp->type = RNA_TYPE_other;
          rrp->ext.choice = 1;
          rrp->ext.value.ptrvalue = StringSave ("ncRNA");
          gbq = GBQualNew ();
          gbq->qual = StringSave ("ncRNA_class");
          gbq->val = StringSave ("snRNA");
          break;
        case FEATDEF_scRNA:
          rrp->type = RNA_TYPE_other;
          rrp->ext.choice = 1;
          rrp->ext.value.ptrvalue = StringSave ("ncRNA");
          gbq = GBQualNew ();
          gbq->qual = StringSave ("ncRNA_class");
          gbq->val = StringSave ("scRNA");
          break;
        case FEATDEF_tmRNA:
          rrp->type = RNA_TYPE_other;
          rrp->ext.choice = 1;
          rrp->ext.value.ptrvalue = StringSave ("tmRNA");
          break;
        case FEATDEF_ncRNA:
          rrp->type = RNA_TYPE_other;
          rrp->ext.choice = 1;
          rrp->ext.value.ptrvalue = StringSave ("ncRNA");
          break;
      }
      break;
    case SEQFEAT_IMP:
      ifp = ImpFeatNew();
      sfp->data.value.ptrvalue = ifp;
      label = GetFeatureNameFromFeatureType (feature_type);
      ifp->key = StringSave (label);
      break;
  }
}


static void ExtraCDSCreationActions (SeqFeatPtr cds, SeqEntryPtr parent_sep)
{
  ByteStorePtr       bs;
  CharPtr            prot, ptr;
  BioseqPtr          bsp;
  Char               ch;
  Int4               i;
  SeqEntryPtr        psep, nsep;
  MolInfoPtr         mip;
  ValNodePtr         vnp, descr;
  SeqFeatPtr         prot_sfp;
  ProtRefPtr         prp;
  Boolean            partial5, partial3;

  if (cds == NULL) return;

  CheckSeqLocForPartial (cds->location, &partial5, &partial3);

  /* Create corresponding protein sequence data for the CDS */

  bs = ProteinFromCdRegionEx (cds, TRUE, FALSE);
  if (NULL == bs)
    return;

  prot = BSMerge (bs, NULL);
  bs = BSFree (bs);
  if (NULL == prot)
    return;

  ptr = prot;
  ch = *ptr;
  while (ch != '\0') {
    *ptr = TO_UPPER (ch);
    ptr++;
    ch = *ptr;
  }
  i = StringLen (prot);
  if (i > 0 && prot [i - 1] == '*') {
    prot [i - 1] = '\0';
  }
  bs = BSNew (1000);
  if (bs != NULL) {
    ptr = prot;
    BSWrite (bs, (VoidPtr) ptr, (Int4) StringLen (ptr));
  }

  /* Create the product protein Bioseq */
  
  bsp = BioseqNew ();
  if (NULL == bsp)
    return;
  
  bsp->repr = Seq_repr_raw;
  bsp->mol = Seq_mol_aa;
  bsp->seq_data_type = Seq_code_ncbieaa;
  bsp->seq_data = (SeqDataPtr) bs;
  bsp->length = BSLen (bs);
  bs = NULL;
  bsp->id = MakeNewProteinSeqId (cds->location, NULL);
  SeqMgrAddToBioseqIndex (bsp);
  
  /* Create a new SeqEntry for the Prot Bioseq */
  
  psep = SeqEntryNew ();
  if (NULL == psep)
    return;
  
  psep->choice = 1;
  psep->data.ptrvalue = (Pointer) bsp;
  SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, psep);
  
  /* Add a descriptor to the protein Bioseq */
  
  mip = MolInfoNew ();
  if (NULL == mip)
    return;
  
  mip->biomol = 8;
  mip->tech = 8;
  if (partial5 && partial3) {
    mip->completeness = 5;
  } else if (partial5) {
    mip->completeness = 3;
  } else if (partial3) {
    mip->completeness = 4;
  }
  vnp = CreateNewDescriptor (psep, Seq_descr_molinfo);
  if (NULL == vnp)
    return;
  
  vnp->data.ptrvalue = (Pointer) mip;
  
  /**/
  
  descr = ExtractBioSourceAndPubs (parent_sep);

  AddSeqEntryToSeqEntry (parent_sep, psep, TRUE);
  nsep = FindNucSeqEntry (parent_sep);
  ReplaceBioSourceAndPubs (parent_sep, descr);
  SetSeqFeatProduct (cds, bsp);
  
  prp = ProtRefNew ();
  
  if (prp != NULL) {
    prot_sfp = CreateNewFeature (psep, NULL, SEQFEAT_PROT, NULL);
    if (prot_sfp != NULL) {
      prot_sfp->data.value.ptrvalue = (Pointer) prp;
      SetSeqLocPartial (prot_sfp->location, partial5, partial3);
      prot_sfp->partial = (partial5 || partial3);
    }
  }
}


static SeqLocPtr LocationFromApplyFeatureAction (BioseqPtr bsp, ApplyFeatureActionPtr action)
{
  LocationIntervalPtr l;
  SeqLocPtr slp = NULL;
  Uint1 strand = Seq_strand_plus;
  Int4  from, to;

  if (bsp == NULL || action == NULL || action->location == NULL) return NULL;

  if (!action->plus_strand) {
    strand = Seq_strand_minus;
  }
  if (action->location->choice == LocationChoice_interval) {
    l = (LocationIntervalPtr) action->location->data.ptrvalue;
    if (l != NULL) {
      from = MIN (l->from, l->to);
      to = MAX (l->from, l->to);
      slp = SeqLocIntNew (from, to, strand, SeqIdFindWorst (bsp->id));
    }
  } else if (action->location->choice == LocationChoice_whole_sequence) {
    slp = SeqLocIntNew (0, bsp->length - 1, strand, SeqIdFindWorst (bsp->id));
  }
  SetSeqLocPartial (slp, action->partial5, action->partial3);
  return slp;
}


static Boolean OkToApplyToBioseq (ApplyFeatureActionPtr action, BioseqPtr bsp)
{
  SeqFeatPtr sfp;
  SeqMgrFeatContext context;
  Int4 featdef;
  Boolean rval = TRUE;

  if (action == NULL || bsp == NULL) return FALSE;

  if (!action->add_redundant) {
    featdef = GetFeatdefFromFeatureType (action->type);
    sfp = SeqMgrGetNextFeature (bsp, NULL, 0, featdef, &context);
    if (sfp != NULL) {
      rval = FALSE;
    }
  }
  return rval;
} 

static void AddParts (ApplyFeatureActionPtr action, BioseqSetPtr parts, ValNodePtr PNTR bsp_list)
{
  SeqEntryPtr sep;
  Int4         seg_num;

  if (action == NULL || !action->apply_to_parts
      || parts == NULL || parts->_class != BioseqseqSet_class_parts
      || bsp_list == NULL) {
    return;
  }

  if (action->only_seg_num > -1) {
    seg_num = 0;
    sep = parts->seq_set;
    while (seg_num < action->only_seg_num && sep != NULL) {
      sep = sep->next;
      seg_num++;
    }
    if (sep != NULL && IS_Bioseq (sep) && OkToApplyToBioseq (action, sep->data.ptrvalue)) {
      ValNodeAddPointer (bsp_list, OBJ_BIOSEQ, sep->data.ptrvalue);
    }
  } else {
    for (sep = parts->seq_set; sep != NULL; sep = sep->next) {
      if (IS_Bioseq (sep) && OkToApplyToBioseq (action, sep->data.ptrvalue)) {
        ValNodeAddPointer (bsp_list, OBJ_BIOSEQ, sep->data.ptrvalue);
      }
    }
  }  
}


static void AddSequenceOrParts (ApplyFeatureActionPtr action, BioseqPtr bsp, ValNodePtr PNTR bsp_list)
{
  BioseqSetPtr bssp, parts;
  SeqEntryPtr sep;

  if (action == NULL || bsp == NULL || bsp_list == NULL) return;

  if (bsp->idx.parenttype == OBJ_BIOSEQSET && bsp->idx.parentptr != NULL) {
    bssp = (BioseqSetPtr) bsp->idx.parentptr;
    if (bssp->_class == BioseqseqSet_class_segset) {
      if (action->apply_to_parts) {
        sep = bssp->seq_set;
        while (sep != NULL && !IS_Bioseq_set (sep)) {
          sep = sep->next;
        }
        if (sep != NULL) {
          AddParts (action, sep->data.ptrvalue, bsp_list);
        }
      } else {
        if (OkToApplyToBioseq (action, bsp)) {
          ValNodeAddPointer (bsp_list, OBJ_BIOSEQ, bsp);
        }
      }       
    } else if (bssp->_class == BioseqseqSet_class_parts) {
      if (action->apply_to_parts) {
        AddParts (action, bssp, bsp_list);
      } else {
        parts = bssp;
        if (parts->idx.parenttype == OBJ_BIOSEQSET && parts->idx.parentptr != NULL) {
          bssp = (BioseqSetPtr) parts->idx.parentptr;
          if (IS_Bioseq (bssp->seq_set) && OkToApplyToBioseq (action, bssp->seq_set->data.ptrvalue)) {
            ValNodeAddPointer (bsp_list, OBJ_BIOSEQ, bsp_list);
          }
        }
      }
    } else {
      if (OkToApplyToBioseq (action, bsp)) {
        ValNodeAddPointer (bsp_list, OBJ_BIOSEQ, bsp);
      }
    }
  } else {
    if (OkToApplyToBioseq (action, bsp)) {
      ValNodeAddPointer (bsp_list, OBJ_BIOSEQ, bsp);
    }
  }
}

static void AddSequenceOrPartsFromSeqEntry (ApplyFeatureActionPtr action, SeqEntryPtr sep, ValNodePtr PNTR bsp_list)
{
  BioseqSetPtr bssp;
  SeqEntryPtr  seq_set;

  if (action == NULL || sep == NULL) return;

  while (sep != NULL) {
    if (IS_Bioseq (sep)) {
      AddSequenceOrParts (action, sep->data.ptrvalue, bsp_list);
    } else if (IS_Bioseq_set (sep)) {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      if (bssp->_class == BioseqseqSet_class_segset) {
        /* find master segment */
        seq_set = bssp->seq_set;
        while (seq_set != NULL && !IS_Bioseq (seq_set)) {
          seq_set = seq_set->next;
        }
        if (seq_set != NULL) {
          AddSequenceOrParts (action, seq_set->data.ptrvalue, bsp_list);
        }
      } else if (bssp->_class == BioseqseqSet_class_nuc_prot) {
        /* find nucleotide sequence */
        seq_set = bssp->seq_set;
        if (seq_set != NULL) {
          if (IS_Bioseq_set (seq_set)) {
            /* nucleotide is segmented set */
            bssp = (BioseqSetPtr) seq_set->data.ptrvalue;
            if (bssp != NULL && bssp->_class == BioseqseqSet_class_segset
                && bssp->seq_set != NULL && IS_Bioseq (bssp->seq_set)) {
              AddSequenceOrParts (action, bssp->seq_set->data.ptrvalue, bsp_list);
            }
          } else if (IS_Bioseq (seq_set)) {
            AddSequenceOrParts (action, seq_set->data.ptrvalue, bsp_list);
          }
        }
      } else {
        /* add from set members */
        AddSequenceOrPartsFromSeqEntry (action, bssp->seq_set, bsp_list);
      }
    }
    sep = sep->next;
  }  
}  
  

static void AdjustProteinSequenceForReadingFrame (SeqFeatPtr cds)
{
  BioseqPtr protbsp, bsp;
  ByteStorePtr bs;
  SeqFeatPtr   prot_sfp;
  Boolean      partial5, partial3;

  if (cds == NULL || cds->data.choice != SEQFEAT_CDREGION) return;

  protbsp = BioseqFindFromSeqLoc (cds->product);

  if (protbsp == NULL) {
    bsp = BioseqFindFromSeqLoc (cds->location);
    if (bsp != NULL) {
      ExtraCDSCreationActions (cds, GetBestTopParentForData (bsp->idx.entityID, bsp));
    }
  } else {
    bs = ProteinFromCdRegionExWithTrailingCodonHandling (cds,
                                              TRUE,
                                              FALSE,
                                              FALSE);
    protbsp->seq_data = (SeqDataPtr) BSFree ((ByteStorePtr)(protbsp->seq_data));
    protbsp->seq_data = (SeqDataPtr) bs;
    protbsp->length = BSLen (bs);
    prot_sfp = GetProtFeature (protbsp);
    if (prot_sfp == NULL) {
      CheckSeqLocForPartial (cds->location, &partial5, &partial3);
      prot_sfp = CreateNewFeatureOnBioseq (protbsp, SEQFEAT_PROT, NULL);
      prot_sfp->data.value.ptrvalue = ProtRefNew ();
      SetSeqLocPartial (prot_sfp->location, partial5, partial3);
      prot_sfp->partial = (partial5 || partial3);
    } else {
      if (SeqLocLen (prot_sfp->location) != protbsp->length) {
        prot_sfp->location = SeqLocFree (prot_sfp->location);
        prot_sfp->location = SeqLocIntNew (0, protbsp->length - 1, Seq_strand_plus, SeqIdFindWorst (protbsp->id));   
      }
    }
  }
}


static Int4 ApplyApplyFeatureActionToSeqEntry (ApplyFeatureActionPtr action, SeqEntryPtr sep)
{
  ValNodePtr bsp_list = NULL, vnp, field_vnp;
  Int4       featdef, seqfeattype;
  BioseqPtr  bsp;
  SeqFeatPtr sfp;
  SeqLocPtr  slp;
  FeatQualLegalValPtr q;
  FeatureField f;
  SeqIdPtr   sip;
  SeqFeatPtr gene;
  Int4       num_created = 0;

  if (sep == NULL || action == NULL) return 0;

  /* first, get list of Bioseqs to apply features to */
  /* relevant values : seq_list, add_redundant, apply_to_parts, only_seg_num */
  if (action->seq_list != NULL && action->seq_list->choice == SequenceListChoice_list) {
    for (vnp = action->seq_list->data.ptrvalue; vnp != NULL; vnp = vnp->next) {
      sip = CreateSeqIdFromText (vnp->data.ptrvalue, sep);
      bsp = BioseqFind (sip);
      if (bsp != NULL) {
        AddSequenceOrParts (action, bsp, &bsp_list);
      }
    }  
  } else {
    AddSequenceOrPartsFromSeqEntry (action, sep, &bsp_list);
  }

  /* now add feature to each bioseq in list */
  for (vnp = bsp_list; vnp != NULL; vnp = vnp->next) {
    bsp = vnp->data.ptrvalue;
    if (bsp == NULL) continue;
    featdef = GetFeatdefFromFeatureType (action->type);
    seqfeattype = FindFeatFromFeatDefType (featdef);
    slp = LocationFromApplyFeatureAction (bsp, action);
    sfp = CreateNewFeatureOnBioseq (bsp, seqfeattype, slp);
    if (sfp == NULL) continue;
    CreateDataForFeature (sfp, action->type);
    /* any extra actions */
    switch (action->type) {
      case (Feature_type_cds) :
        ExtraCDSCreationActions (sfp, GetBestTopParentForData (bsp->idx.entityID, bsp));
        break;
    }
    gene = NULL;
    for (field_vnp = action->fields; field_vnp != NULL; field_vnp = field_vnp->next) {
      q = (FeatQualLegalValPtr) field_vnp->data.ptrvalue;
      if (q != NULL) {
        f.field = ValNodeNew(NULL);
        f.field->next = NULL;
        f.field->choice = FeatQualChoice_legal_qual;
        f.field->data.intvalue = q->qual;        
        if (sfp->data.choice != SEQFEAT_GENE
            && (q->qual == Feat_qual_legal_gene || q->qual == Feat_qual_legal_gene_description)) {
          if (gene == NULL) {
            gene = CreateNewFeatureOnBioseq (bsp, SEQFEAT_GENE, slp);
            CreateDataForFeature (gene, Feature_type_gene);
          }
          f.type = Feature_type_gene;
          SetQualOnFeature (gene, &f, NULL, q->val, ExistingTextOption_replace_old);
        } else {
          f.type = action->type;
          SetQualOnFeature (sfp, &f, NULL, q->val, ExistingTextOption_replace_old);
        }
      }
    }
    if (action->type == Feature_type_cds) {
      /* retranslate, to account for change in reading frame */
      AdjustProteinSequenceForReadingFrame (sfp);
      /* after the feature has been created, then adjust it for gaps */
      /* Note - this step may result in multiple coding regions being created. */
      AdjustCDSLocationsForUnknownGapsCallback (sfp, NULL);
    }
    num_created++;
  }  
  return num_created;
}


typedef struct convertandremovefeaturecollection {
  Uint1 featdef;
  ValNodePtr constraint_set;
  ValNodePtr feature_list;
} ConvertAndRemoveFeatureCollectionData, PNTR ConvertAndRemoveFeatureCollectionPtr;

static void ConvertAndRemoveFeatureCollectionCallback (SeqFeatPtr sfp, Pointer data)
{
  ConvertAndRemoveFeatureCollectionPtr p;  

  if (sfp == NULL || data == NULL) return;

  p = (ConvertAndRemoveFeatureCollectionPtr) data;
  if (sfp->idx.subtype == p->featdef && DoesObjectMatchConstraintChoiceSet (OBJ_SEQFEAT, sfp, p->constraint_set)) {
    ValNodeAddPointer (&(p->feature_list), OBJ_SEQFEAT, sfp);
  }
}


static Int4 ApplyRemoveFeatureActionToSeqEntry (RemoveFeatureActionPtr action, SeqEntryPtr sep)
{
  ConvertAndRemoveFeatureCollectionData d;
  ValNodePtr vnp;
  SeqFeatPtr sfp;
  Int4       num_deleted = 0;

  if (action == NULL) return 0;

  d.featdef = GetFeatdefFromFeatureType (action->type);
  d.constraint_set = action->constraint;
  d.feature_list = NULL;

  VisitFeaturesInSep (sep, &d, ConvertAndRemoveFeatureCollectionCallback);
  for (vnp = d.feature_list; vnp != NULL; vnp = vnp->next) {
    sfp = vnp->data.ptrvalue;
    if (sfp != NULL) {
      sfp->idx.deleteme = TRUE;
      num_deleted ++;
    }
  }
  DeleteMarkedObjects (ObjMgrGetEntityIDForChoice(sep), 0, NULL);
  return num_deleted;
}


static Boolean DoesStrandMatch (Int4 strand_choice, Uint1 strand_val)
{
  Boolean rval = FALSE;
  
  switch (strand_choice)
  {
    case Feature_location_strand_from_any:
      rval = TRUE;
      break;
    case Feature_location_strand_from_unknown:
      if (strand_val == Seq_strand_unknown)
      {
        rval = TRUE;
      }
      break;
    case Feature_location_strand_from_plus:
      if (strand_val != Seq_strand_minus)
      {
        rval = TRUE;
      }
      break;
    case Feature_location_strand_from_minus:
      if (strand_val == Seq_strand_minus)
      {
        rval = TRUE;
      }
      break;
    case Feature_location_strand_from_both:
      if (strand_val == Seq_strand_both)
      {
        rval = TRUE;
      }
      break;
  }
  return rval;
}


static Uint1 GetNewStrandValue (Int4 strand_choice, Uint1 strand_val)
{
  Uint1 rval = Seq_strand_unknown;
  
  switch (strand_choice)
  {
    case Feature_location_strand_to_reverse:
      switch (strand_val)
      {
        case Seq_strand_plus:
        case Seq_strand_unknown:
          rval = Seq_strand_minus;
          break;
        case Seq_strand_minus:
          rval = Seq_strand_plus;
          break;
        default:
          rval = strand_val;
          break;
      }
      break;
    case Feature_location_strand_to_unknown:
      rval = Seq_strand_unknown;
      break;
    case Feature_location_strand_to_plus:
      rval = Seq_strand_plus;
      break;
    case Feature_location_strand_to_minus:
      rval = Seq_strand_minus;
      break;
    case Feature_location_strand_to_both:
      rval = Seq_strand_both;
      break;
  }  
  return rval;
}


static Boolean ConvertLocationStrand (SeqLocPtr slp, Int4 fromStrand, Int4 toStrand)
{
  SeqLocPtr      loc;
  PackSeqPntPtr  psp;
  SeqBondPtr     sbp;
  SeqIntPtr      sinp;
  SeqPntPtr      spp;
  Boolean        rval = FALSE;
  Uint1          strand_orig;

  while (slp != NULL) {
    switch (slp->choice) {
      case SEQLOC_NULL :
        break;
      case SEQLOC_EMPTY :
      case SEQLOC_WHOLE :
        break;
      case SEQLOC_INT :
        sinp = (SeqIntPtr) slp->data.ptrvalue;
        if (sinp != NULL && DoesStrandMatch (fromStrand, sinp->strand)) 
        {
          strand_orig = sinp->strand;
          sinp->strand = GetNewStrandValue (toStrand, sinp->strand);
          if (strand_orig != sinp->strand) {
            rval = TRUE;
          }
        }
        break;
      case SEQLOC_PNT :
        spp = (SeqPntPtr) slp->data.ptrvalue;
        if (spp != NULL && DoesStrandMatch (fromStrand, spp->strand))
        {
          strand_orig = spp->strand;
          spp->strand = GetNewStrandValue (toStrand, spp->strand);
          if (strand_orig != spp->strand) {
            rval = TRUE;
          }
        }
        break;
      case SEQLOC_PACKED_PNT :
        psp = (PackSeqPntPtr) slp->data.ptrvalue;
        if (psp != NULL && DoesStrandMatch (fromStrand, psp->strand)) 
        {
          strand_orig = psp->strand;
          psp->strand = GetNewStrandValue (toStrand, psp->strand);
          if (strand_orig != psp->strand) {
            rval = TRUE;
          }
        }
        break;
      case SEQLOC_PACKED_INT :
      case SEQLOC_MIX :
      case SEQLOC_EQUIV :
        loc = (SeqLocPtr) slp->data.ptrvalue;
        while (loc != NULL) {
          rval |= ConvertLocationStrand (loc, fromStrand, toStrand);
          loc = loc->next;
        }
        break;
      case SEQLOC_BOND :
        sbp = (SeqBondPtr) slp->data.ptrvalue;
        if (sbp != NULL) {
          spp = (SeqPntPtr) sbp->a;
          if (spp != NULL && DoesStrandMatch (fromStrand, spp->strand)) 
          {
            strand_orig = spp->strand;
            spp->strand = GetNewStrandValue (toStrand, spp->strand);
            if (strand_orig != spp->strand) {
              rval = TRUE;
            }
          }
          spp = (SeqPntPtr) sbp->b;
          if (spp != NULL && DoesStrandMatch (fromStrand, spp->strand)) 
          {
            strand_orig = spp->strand;
            spp->strand = GetNewStrandValue (toStrand, spp->strand);
            if (strand_orig != spp->strand) {
              rval = TRUE;
            }
          }
        }
        break;
      case SEQLOC_FEAT :
        break;
      default :
        break;
    }
    slp = slp->next;
  }
  return rval;
}


static Boolean ApplyEditLocationStrandToSeqFeat (EditLocationStrandPtr edit, SeqFeatPtr sfp)
{
  Boolean rval = FALSE;

  if (edit == NULL || sfp == NULL) {
    return FALSE;
  }

  rval = ConvertLocationStrand (sfp->location, edit->strand_from, edit->strand_to);  
  return rval;
}


static Boolean At5EndOfSequence (SeqLocPtr slp, BioseqPtr bsp)
{
  Uint1 strand;
  Int4  start;
  Boolean at_end = FALSE;

  if (slp == NULL || bsp == NULL) return FALSE;

  strand = SeqLocStrand (slp);

  if (strand == Seq_strand_minus) {
    start = SeqLocStop (slp);
    if (start == bsp->length - 1) {
      at_end = TRUE;
    }
  } else {
    start = SeqLocStart (slp);
    if (start == 0) {
      at_end = TRUE;
    }
  }
  return at_end;
}


static Boolean HasGoodStartCodon (SeqFeatPtr sfp)
{
  ByteStorePtr bs;
  CharPtr      prot;
  Boolean     has_start = FALSE;

  if (sfp != NULL && sfp->data.choice == SEQFEAT_CDREGION) {
    bs = ProteinFromCdRegionEx (sfp, TRUE, FALSE);
    if (bs != NULL) {
      prot = BSMerge (bs, NULL);
      bs = BSFree (bs);
      if (prot != NULL && *prot == 'M') {
        has_start = TRUE;
      }
      prot = MemFree (prot);
    }
  }
  return has_start;
}


static Boolean ApplyPartial5SetActionToSeqFeat (Partial5SetActionPtr action, SeqFeatPtr sfp)
{
  Boolean      rval = FALSE;
  Boolean      make_partial = FALSE;
  Uint1        strand;
  BioseqPtr    bsp;
  CdRegionPtr  crp;
  Boolean      partial5, partial3;

  if (action == NULL || sfp == NULL) return FALSE;
  bsp = BioseqFindFromSeqLoc (sfp->location);
  strand = SeqLocStrand (sfp->location);

  switch (action->constraint) {
    case Partial_5_set_constraint_all:
      make_partial = TRUE;
      break;
    case Partial_5_set_constraint_at_end:
      make_partial = At5EndOfSequence (sfp->location, bsp);
      break;
    case Partial_5_set_constraint_bad_start:
      make_partial = HasGoodStartCodon (sfp);
      break;
    case Partial_5_set_constraint_frame_not_one:
      if (sfp->data.choice == SEQFEAT_CDREGION
          && (crp = sfp->data.value.ptrvalue) != NULL
          && crp->frame != 0 && crp->frame != 1) {
        make_partial = TRUE;
      }
      break;
  }

  if (make_partial) {
    CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
    if (!partial5) {
      SetSeqLocPartial (sfp->location, TRUE, partial3);
      if (action->extend && bsp != NULL) {
        ExtendSeqLocToEnd (sfp->location, bsp, TRUE);
      }
      rval = TRUE; 
    }
  }
  return rval;
}


static Boolean ApplyClear5PartialToSeqFeat (Int4 action, SeqFeatPtr sfp)
{
  Boolean rval = FALSE, clear_partial = FALSE;
  Boolean partial5, partial3;

  if (sfp == NULL) return FALSE;

  switch (action) {
    case Partial_5_clear_constraint_all:
      clear_partial = TRUE;
      break;
    case Partial_5_clear_constraint_not_at_end:
      clear_partial = !At5EndOfSequence(sfp->location, BioseqFindFromSeqLoc (sfp->location));
      break;
    case Partial_5_clear_constraint_good_start:
      clear_partial = !HasGoodStartCodon(sfp);
      break;
  }
  if (clear_partial) {
    CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
    if (partial5) {
      SetSeqLocPartial (sfp->location, FALSE, partial3);
      rval = TRUE;
    }
  }
  return rval;
}


static Boolean At3EndOfSequence (SeqLocPtr slp, BioseqPtr bsp)
{
  Uint1 strand;
  Int4  stop;
  Boolean at_end = FALSE;

  if (slp == NULL || bsp == NULL) return FALSE;

  strand = SeqLocStrand (slp);

  if (strand == Seq_strand_minus) {
    stop = SeqLocStart (slp);
    if (stop == 0) {
      at_end = TRUE;
    }
  } else {
    stop = SeqLocStop (slp);
    if (stop == bsp->length - 1) {
      at_end = TRUE;
    }
  }
  return at_end;
}


static Boolean HasGoodStopCodon (SeqFeatPtr sfp)
{
  ByteStorePtr bs;
  CharPtr      prot;
  Boolean      has_stop = FALSE;

  if (sfp != NULL && sfp->data.choice == SEQFEAT_CDREGION) {
    bs = ProteinFromCdRegionEx (sfp, TRUE, FALSE);
    if (bs != NULL) {
      prot = BSMerge (bs, NULL);
      bs = BSFree (bs);
      if (prot != NULL && prot[StringLen (prot) - 1] == '*') {
        has_stop = TRUE;
      }
      prot = MemFree (prot);
    }
  }
  return has_stop;
}


static Boolean ApplyPartial3SetActionToSeqFeat (Partial3SetActionPtr action, SeqFeatPtr sfp)
{
  Boolean      rval = FALSE;
  Boolean      make_partial = FALSE;
  Uint1        strand;
  BioseqPtr    bsp;
  Boolean      partial5, partial3;

  if (action == NULL || sfp == NULL) return FALSE;
  bsp = BioseqFindFromSeqLoc (sfp->location);
  strand = SeqLocStrand (sfp->location);

  switch (action->constraint) {
    case Partial_3_set_constraint_all:
      make_partial = TRUE;
      break;
    case Partial_3_set_constraint_at_end:
      make_partial = At3EndOfSequence (sfp->location, bsp);
      break;
    case Partial_3_set_constraint_bad_end:
      make_partial = HasGoodStopCodon (sfp);
      break;
  }

  if (make_partial) {
    CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
    if (!partial3) {
      SetSeqLocPartial (sfp->location, partial5, TRUE);
      if (action->extend && bsp != NULL) {
        ExtendSeqLocToEnd (sfp->location, bsp, FALSE);
      }
      rval = TRUE; 
    }
  }
  return rval;
}


static Boolean ApplyClear3PartialToSeqFeat (Int4 action, SeqFeatPtr sfp)
{
  Boolean rval = FALSE, clear_partial = FALSE;
  Boolean partial5, partial3;

  if (sfp == NULL) return FALSE;

  switch (action) {
    case Partial_3_clear_constraint_all:
      clear_partial = TRUE;
      break;
    case Partial_3_clear_constraint_not_at_end:
      clear_partial = !At3EndOfSequence(sfp->location, BioseqFindFromSeqLoc (sfp->location));
      break;
    case Partial_3_clear_constraint_good_end:
      clear_partial = !HasGoodStopCodon(sfp);
      break;
  }
  if (clear_partial) {
    CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
    if (partial3) {
      SetSeqLocPartial (sfp->location, partial5, FALSE);
      rval = TRUE;
    }
  }
  return rval;
}


static Boolean ApplyConvertLocationToSeqFeat (Int4 convert_location, SeqFeatPtr sfp)
{
  Boolean hasNulls, rval = FALSE;
  SeqLocPtr slp;
  BioseqPtr bsp;
  Boolean   partial5, partial3;

  if (sfp == NULL || (bsp = BioseqFindFromSeqLoc (sfp->location))== NULL) {
    return FALSE;
  }

  CheckSeqLocForPartial (sfp->location, &partial5, &partial3);
	hasNulls = LocationHasNullsBetween (sfp->location);
	switch (convert_location) 
	{
	  case Convert_location_type_join :
	    if (hasNulls) 
	    {
	      slp = SeqLocMerge (bsp, sfp->location, NULL, FALSE, FALSE, FALSE);
		    sfp->location = SeqLocFree (sfp->location);
		    sfp->location = slp;
		    if (bsp->repr == Seq_repr_seg) 
		    {
		      slp = SegLocToPartsEx (bsp, sfp->location, FALSE);
		      sfp->location = SeqLocFree (sfp->location);
		      sfp->location = slp;
		      hasNulls = LocationHasNullsBetween (sfp->location);
		      sfp->partial = (sfp->partial || hasNulls);
		    }
		    FreeAllFuzz (sfp->location);
		    SetSeqLocPartial (sfp->location, partial5, partial3);
        rval = TRUE;
	    }
	    break;
  	case Convert_location_type_order :
	    if (!hasNulls) 
	    {
		    slp = SeqLocMerge (bsp, sfp->location, NULL, FALSE, FALSE, TRUE);
        sfp->location = SeqLocFree (sfp->location);
		    sfp->location = slp;
		    if (bsp->repr == Seq_repr_seg) 
		    {
		      slp = SegLocToPartsEx (bsp, sfp->location, TRUE);
		      sfp->location = SeqLocFree (sfp->location);
		      sfp->location = slp;
		      hasNulls = LocationHasNullsBetween (sfp->location);
		      sfp->partial = (sfp->partial || hasNulls);
		    }
		    FreeAllFuzz (sfp->location);
		    SetSeqLocPartial (sfp->location, partial5, partial3);
        rval = TRUE;
	    }
	    break;
	  case Convert_location_type_merge :
      if (sfp->location->choice != SEQLOC_INT) {
	      slp = SeqLocMerge (bsp, sfp->location, NULL, TRUE, FALSE, FALSE);
	      sfp->location = SeqLocFree (sfp->location);
	      sfp->location = slp;
		    SetSeqLocPartial (sfp->location, partial5, partial3);
        rval = TRUE;
      }
	  default:
	    break;
	}
  return rval;
}


static Boolean ApplyLocationEditTypeToSeqFeat (ValNodePtr action, SeqFeatPtr sfp)
{
  Boolean rval = FALSE;

  if (action == NULL || sfp == NULL) {
    return FALSE;
  }

  switch (action->choice) {
    case LocationEditType_strand:
      rval = ApplyEditLocationStrandToSeqFeat (action->data.ptrvalue, sfp);
      break;
    case LocationEditType_set_5_partial:
      rval = ApplyPartial5SetActionToSeqFeat (action->data.ptrvalue, sfp);
      break;
    case LocationEditType_clear_5_partial:
      rval = ApplyClear5PartialToSeqFeat (action->data.intvalue, sfp);
      break;
    case LocationEditType_set_3_partial:
      rval = ApplyPartial3SetActionToSeqFeat (action->data.ptrvalue, sfp);
      break;
    case LocationEditType_clear_3_partial:
      rval = ApplyClear3PartialToSeqFeat (action->data.intvalue, sfp);
      break;
    case LocationEditType_convert:
      rval = ApplyConvertLocationToSeqFeat (action->data.intvalue, sfp);
      break;
  }
  return rval;
}


static Int4 ApplyEditFeatureLocationActionToSeqEntry (EditFeatureLocationActionPtr action, SeqEntryPtr sep)
{
  ConvertAndRemoveFeatureCollectionData d;
  ValNodePtr vnp;
  SeqFeatPtr sfp;
  Int4       num_affected = 0;

  if (action == NULL) return 0;

  d.featdef = GetFeatdefFromFeatureType (action->type);
  d.constraint_set = action->constraint;
  d.feature_list = NULL;

  VisitFeaturesInSep (sep, &d, ConvertAndRemoveFeatureCollectionCallback);
  for (vnp = d.feature_list; vnp != NULL; vnp = vnp->next) {
    sfp = vnp->data.ptrvalue;
    if (sfp != NULL && ApplyLocationEditTypeToSeqFeat (action->action, sfp)) {
      num_affected++;
    }
  }
  return num_affected;
}


NLM_EXTERN void ApplyMacroToSeqEntry (SeqEntryPtr sep, ValNodePtr macro, Int4Ptr pNumFields, Int4Ptr pNumFeat)
{
  Int4 num_AECR = 0, num_parse = 0, num_feature = 0, num_fields = 0;

  while (macro != NULL) {
    switch (macro->choice) {
      case MacroActionChoice_aecr:
        num_AECR += ApplyAECRActionToSeqEntry ((AECRActionPtr) macro->data.ptrvalue, sep);
        break;
      case MacroActionChoice_parse:
        num_parse += ApplyParseActionToSeqEntry ((ParseActionPtr) macro->data.ptrvalue, sep);
        break;
      case MacroActionChoice_add_feature:
        num_feature += ApplyApplyFeatureActionToSeqEntry ((ApplyFeatureActionPtr) macro->data.ptrvalue, sep);
        SeqMgrIndexFeatures (ObjMgrGetEntityIDForChoice(sep), NULL);
        break;
      case MacroActionChoice_remove_feature:
        num_feature += ApplyRemoveFeatureActionToSeqEntry ((RemoveFeatureActionPtr) macro->data.ptrvalue, sep);
        break;
      case MacroActionChoice_edit_location:
        num_fields += ApplyEditFeatureLocationActionToSeqEntry ((EditFeatureLocationActionPtr) macro->data.ptrvalue, sep);
        break;
    }
    macro = macro->next;
  }
  if (pNumFields != NULL) {
    *pNumFields = num_AECR + num_parse + num_fields;
  }
  if (pNumFeat != NULL) {
    *pNumFeat = num_feature;
  }
}


/* for generating text descriptions of macro objects */
NLM_EXTERN CharPtr SummarizeSourceQual (ValNodePtr field)
{
  CharPtr summ = NULL, locname, origname;
  Int4    genome, origin;
  CharPtr loc_fmt = "location %s";
  CharPtr orig_fmt = "origin %s";

  if (field == NULL) return NULL;
  switch (field->choice) {
    case SourceQualChoice_textqual:
      summ = StringSave (GetSourceQualName (field->data.intvalue));
      break;
    case SourceQualChoice_location:
      genome = GenomeFromSrcLoc (field->data.intvalue);
      locname = LocNameFromGenome (genome);
      summ = (CharPtr) MemNew (sizeof (Char) * (StringLen (loc_fmt) + StringLen (locname)));
      sprintf (summ, loc_fmt, locname);
      break;
    case SourceQualChoice_origin:
      origin = OriginFromSrcOrig (field->data.intvalue);
      origname = OriginNameFromOrigin (origin);
      summ = (CharPtr) MemNew (sizeof (Char) * (StringLen (orig_fmt) + StringLen (origname)));
      sprintf (summ, orig_fmt, origname);
      break;
  }
  return summ;
}


NLM_EXTERN CharPtr FeatureFieldLabel (CharPtr feature_name, ValNodePtr field)
{
  CharPtr cp;
  CharPtr label = NULL;
  CharPtr legal_fmt = "%s %s";
  CharPtr illegal_fmt = "constrained field on %s";
  
  if (feature_name == NULL) {
    feature_name = "Unknown feature";
  }

  if (field == NULL) {
    return StringSave ("missing field");
  } else if (field->choice == FeatQualChoice_legal_qual) {
    cp = GetFeatQualName (field->data.intvalue);
    if (cp == NULL) cp = "Unknown field type";
    label = (CharPtr) MemNew (sizeof (Char) * (StringLen (legal_fmt) + StringLen (feature_name) + StringLen (cp)));
    sprintf (label, legal_fmt, feature_name, cp);
  } else if (field->choice == FeatQualChoice_illegal_qual) {
    label = (CharPtr) MemNew (sizeof (Char) * (StringLen (illegal_fmt) + StringLen (feature_name)));
    sprintf (label, illegal_fmt, feature_name);
  } else {
    label = StringSave ("illegal field value");
  }
  return label;
}


NLM_EXTERN Boolean IsFeatureFieldEmpty (FeatureFieldPtr field)
{
  if (field == NULL) return TRUE;
  if (field->field == NULL) return TRUE;
  return FALSE;
}


NLM_EXTERN Boolean IsFieldTypeEmpty (FieldTypePtr field)
{
  Boolean rval = TRUE;

  if (field == NULL) return TRUE;
  switch (field->choice) {
    case FieldType_source_qual:
      if (field->data.ptrvalue != NULL) {
        rval = FALSE;
      }
      break;
    case FieldType_feature_field:
      if (!IsFeatureFieldEmpty (field->data.ptrvalue)) {
        rval = FALSE;
      }
      break;
    case FieldType_cds_gene_prot:
      rval = FALSE;
      break;
  }
  return rval;
}


NLM_EXTERN CharPtr SummarizeFieldType (ValNodePtr vnp)
{
  FeatureFieldPtr ffp;
  CharPtr str = NULL;
  CharPtr    label = NULL;

  if (vnp == NULL) {
    str = StringSave ("missing field");
  } else {
    switch (vnp->choice) {
      case FieldType_source_qual:
        str = SummarizeSourceQual (vnp->data.ptrvalue);
        break;
      case FieldType_feature_field:
        ffp = (FeatureFieldPtr) vnp->data.ptrvalue;
        if (ffp == NULL || ffp->field == NULL) {
          str = StringSave ("missing field");
        } else {
          label = GetFeatureNameFromFeatureType (ffp->type);
          str = FeatureFieldLabel (label, ffp->field);
        }
        break;
      case FieldType_cds_gene_prot:
        str = StringSaveNoNull (CDSGeneProtNameFromField (vnp->data.intvalue));
        if (str == NULL) {
          str = StringSave ("Invalid CDS-Gene-Prot Field");
        }
        break;
      case FieldType_molinfo_field:
        str = GetSequenceQualName (vnp->data.ptrvalue);
        if (str == NULL) {
          str = StringSave ("Invalid Sequence Qual Field");
        }
        break;
      default:
        str = StringSave ("Invalid field type");
        break;
    }
  }
  return str;
}


/* for table readers that use the macro language functions */
NLM_EXTERN TabColumnConfigPtr TabColumnConfigNew (void)
{
  TabColumnConfigPtr t;

  t = (TabColumnConfigPtr) MemNew (sizeof (TabColumnConfigData));
  t->field = NULL;
  t->existing_text = ExistingTextOption_replace_old;
  t->skip_blank = TRUE;
  return t;
}



NLM_EXTERN TabColumnConfigPtr TabColumnConfigFree (TabColumnConfigPtr t)
{
  if (t != NULL) {
    t->field = FieldTypeFree (t->field);
    t = MemFree (t);
  }
  return t;
}


NLM_EXTERN TabColumnConfigPtr TabColumnConfigCopy (TabColumnConfigPtr orig)
{
  TabColumnConfigPtr t = NULL;

  if (orig != NULL) {
    t = TabColumnConfigNew ();
    t->match_type = orig->match_type;
    t->existing_text = orig->existing_text;
    t->skip_blank = orig->skip_blank;
    t->match_mrna = orig->match_mrna;
    t->field = AsnIoMemCopy (orig->field, (AsnReadFunc) FieldTypeAsnRead, (AsnWriteFunc) FieldTypeAsnWrite);
  }
  return t;
}


NLM_EXTERN ValNodePtr TabColumnConfigListFree (ValNodePtr columns)
{
  ValNodePtr vnp_next;

  while (columns != NULL) {
    vnp_next = columns->next;
    columns->data.ptrvalue = TabColumnConfigFree (columns->data.ptrvalue);
    columns->next = NULL;
    columns = ValNodeFree (columns);
    columns = vnp_next;
  }
  return columns;
}


NLM_EXTERN ValNodePtr TabColumnConfigListCopy (ValNodePtr orig)
{
  ValNodePtr new_list = NULL;
  TabColumnConfigPtr t;

  while (orig != NULL) {
    t = TabColumnConfigCopy (orig->data.ptrvalue);
    ValNodeAddPointer (&new_list, 0, t);
    orig = orig->next;
  }
  return new_list;
}



/* This checks the column names and returns a list of the feature fields */
NLM_EXTERN ValNodePtr ValidateFeatureFieldColumnNames (ValNodePtr header_line, ValNodePtr PNTR perr_list)
{
  ValNodePtr         header_vnp;
  ValNodePtr         err_list = NULL, col_list = NULL;
  Boolean            rval = TRUE;
  TabColumnConfigPtr t;
  FeatureFieldPtr    field;
  Int4               featqual, feat_type;
  CharPtr            first_space;
  
  if (header_line == NULL)
  {
    return FALSE;
  }
  
  header_vnp = header_line->data.ptrvalue;
  if (header_vnp == NULL || header_vnp->next == NULL)
  {
    return FALSE;
  }
  
  /* skip ID column */
  header_vnp = header_vnp->next;
  while (header_vnp != NULL && rval)
  {
    first_space = StringChr (header_vnp->data.ptrvalue, ' ');
    if (first_space != NULL) {
      *first_space = 0;
      feat_type = GetFeatureTypeByName (header_vnp->data.ptrvalue);
      featqual = GetFeatQualByName (first_space + 1);
      *first_space = ' ';
      if (feat_type < 0 || featqual < 0) {
        /* unable to recognize column name */
        ValNodeAddPointer (&err_list, 0, StringSave (header_vnp->data.ptrvalue));
        /* if we're not able to send back a list of errors, just quit now */
        if (perr_list == NULL) {
          rval = FALSE;
        }
      } else if (err_list == NULL) {
        /* if we've already found errors, don't bother collecting more fields */
        field = FeatureFieldNew ();
        field->type = feat_type;
        field->field = ValNodeNew (NULL);
        field->field->choice = FeatQualChoice_legal_qual;
        field->field->data.intvalue = featqual;
        t = TabColumnConfigNew ();
        t->field = ValNodeNew (NULL);
        t->field->choice = FieldType_feature_field;
        t->field->data.ptrvalue = field;
        ValNodeAddPointer (&col_list, 0, t);
      }
    } else {
      featqual = GetFeatQualByName (header_vnp->data.ptrvalue);
      if (featqual < 0) {
        /* unable to recognize column name */
        ValNodeAddPointer (&err_list, 0, StringSave (header_vnp->data.ptrvalue));
        /* if we're not able to send back a list of errors, just quit now */
        if (perr_list == NULL) {
          rval = FALSE;
        }
      } else if (err_list == NULL) {
        /* if we've already found errors, don't bother collecting more fields */
        field = FeatureFieldNew ();
        field->type = Feature_type_any;
        field->field = ValNodeNew (NULL);
        field->field->choice = FeatQualChoice_legal_qual;
        field->field->data.intvalue = featqual;
        t = TabColumnConfigNew ();
        t->field = ValNodeNew (NULL);
        t->field->choice = FieldType_feature_field;
        t->field->data.ptrvalue = field;
        ValNodeAddPointer (&col_list, 0, t);
      }
    }
    header_vnp = header_vnp->next;
  }
  if (err_list != NULL) {
    col_list = TabColumnConfigListFree (col_list);
    if (perr_list != NULL) {
      *perr_list = err_list;
    } else {
      err_list = ValNodeFreeData (err_list);
    }
  }
  return col_list;
}

typedef struct findgenelocustag {
  CharPtr locus_tag;
  ValNodePtr gene_list;
} FindGeneLocusTagData, PNTR FindGeneLocusTagPtr;

static void FindGeneByLocusTagBioseqCallback (BioseqPtr bsp, Pointer userdata)
{
  FindGeneLocusTagPtr p;
  SeqFeatPtr          gene;
  SeqMgrFeatContext   fcontext;

  if (bsp == NULL || userdata == NULL || !ISA_na (bsp->mol)) {
    return;
  }

  p = (FindGeneLocusTagPtr) userdata;

  gene = SeqMgrGetGeneByLocusTag (bsp, p->locus_tag, &fcontext);
  if (gene != NULL) {
    ValNodeAddPointer (&p->gene_list, OBJ_SEQFEAT, gene);
  }
}


typedef struct objbystr {
  ValNodePtr obj_list;
  CharPtr    str;
} ObjByStrData, PNTR ObjByStrPtr;

static void GetFeaturesByDbxrefCallback (SeqFeatPtr sfp, Pointer userdata)
{
  ObjByStrPtr p;
  ValNodePtr    vnp;
  DbtagPtr      dbt;
  Char          buf[20];
  Boolean       found = FALSE;

  if (sfp == NULL || sfp->dbxref == NULL || userdata == NULL) return;
  p = (ObjByStrPtr) userdata;

  if (StringHasNoText (p->str)) return;

  for (vnp = sfp->dbxref; vnp != NULL && !found; vnp = vnp->next) {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt != NULL && dbt->tag != NULL) {
      if (dbt->tag->id > 0) {
        sprintf (buf, "%d", dbt->tag->id);
        if (StringCmp (buf, p->str) == 0) {
          found = TRUE;
        }
      } else if (StringCmp (dbt->tag->str, p->str) == 0) {
        found = TRUE;
      }
    }
  }
  if (found) {
    ValNodeAddPointer (&(p->obj_list), OBJ_SEQFEAT, sfp);
  }

}


static ValNodePtr GetFeaturesByDbxref (SeqEntryPtr sep, CharPtr dbxref)
{
  ObjByStrData d;

  d.str = dbxref;
  d.obj_list = NULL;
  VisitFeaturesInSep (sep, &d, GetFeaturesByDbxrefCallback);
  return d.obj_list;
}


static void GetBioSourcesByTaxNameDescriptorCallback (SeqDescrPtr sdp, Pointer userdata)
{
  ObjByStrPtr p;
  BioSourcePtr biop;

  if (sdp == NULL || sdp->choice != Seq_descr_source || userdata == NULL) return;
  p = (ObjByStrPtr) userdata;

  if (StringHasNoText (p->str)) return;

  biop = (BioSourcePtr) sdp->data.ptrvalue;
  if (biop != NULL && biop->org != NULL && StringCmp (biop->org->taxname, p->str) == 0) {
    ValNodeAddPointer (&(p->obj_list), OBJ_SEQDESC, sdp);
  }

}


static void GetBioSourcesByTaxNameFeatureCallback (SeqFeatPtr sfp, Pointer userdata)
{
  ObjByStrPtr p;
  BioSourcePtr biop;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_BIOSRC || userdata == NULL) return;
  p = (ObjByStrPtr) userdata;

  if (StringHasNoText (p->str)) return;

  biop = (BioSourcePtr) sfp->data.value.ptrvalue;
  if (biop != NULL && biop->org != NULL && StringCmp (biop->org->taxname, p->str) == 0) {
    ValNodeAddPointer (&(p->obj_list), OBJ_SEQFEAT, sfp);
  }

}


static ValNodePtr GetBioSourcesByTaxName (SeqEntryPtr sep, CharPtr taxname)
{
  ObjByStrData d;

  d.str = taxname;
  d.obj_list = NULL;
  VisitDescriptorsInSep (sep, &d, GetBioSourcesByTaxNameDescriptorCallback);

  VisitFeaturesInSep (sep, &d, GetBioSourcesByTaxNameFeatureCallback);
  return d.obj_list;
}



static ValNodePtr 
FindMatchForRow 
(ValNodePtr  match_type,
 Uint2       entityID,
 SeqEntryPtr sep)
{
  ValNodePtr match_list = NULL;
  SeqIdPtr   sip;
  BioseqPtr  bsp, nbsp = NULL;
  FindGeneLocusTagData fd;
  SeqFeatPtr           sfp;
  SeqMgrFeatContext    fcontext;

  if (match_type == NULL || sep == NULL) return NULL;

  switch (match_type->choice) {
    case eTableMatchFeatureID:
      sfp = SeqMgrGetFeatureByFeatID (entityID, NULL, match_type->data.ptrvalue, NULL, &fcontext);
      if (sfp != NULL) {
        ValNodeAddPointer (&match_list, OBJ_SEQFEAT, sfp);
      }
      break;
    case eTableMatchGeneLocusTag:
      fd.locus_tag = match_type->data.ptrvalue;
      fd.gene_list = NULL;
      VisitBioseqsInSep (sep, &fd, FindGeneByLocusTagBioseqCallback);
      ValNodeLink (&match_list, fd.gene_list);
      break;
    case eTableMatchProteinID:
    case eTableMatchNucID:
      sip = CreateSeqIdFromText (match_type->data.ptrvalue, sep);
      bsp = BioseqFind (sip);
      sip = SeqIdFree (sip);
      if (bsp != NULL) 
      {
        ValNodeAddPointer (&match_list, OBJ_BIOSEQ, bsp);
      }
      break;
    case eTableMatchDbxref:
      match_list = GetFeaturesByDbxref (sep, match_type->data.ptrvalue);
      break;
    case eTableMatchBioSource:
      match_list = GetBioSourcesByTaxName (sep, match_type->data.ptrvalue);
      break;
  }
  return match_list;
}


static ValNodePtr GetFeatureListForProteinBioseq (Uint1 featdef, BioseqPtr bsp)
{
  ValNodePtr feat_list = NULL;
  SeqFeatPtr sfp, cds;
  SeqMgrFeatContext fcontext;
  Int4              seqfeattype;

  if (bsp == NULL || !ISA_aa (bsp->mol)) 
  {
    return NULL;
  }

  seqfeattype = FindFeatFromFeatDefType (featdef);
  if (seqfeattype == SEQFEAT_PROT)
  {
    for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, featdef, &fcontext);
         sfp != NULL;
         sfp = SeqMgrGetNextFeature (bsp, sfp, 0, featdef, &fcontext))
    {
      ValNodeAddPointer (&feat_list, OBJ_SEQFEAT, sfp);
    }
  }
  else
  {
    cds = SeqMgrGetCDSgivenProduct (bsp, NULL);
    if (cds != NULL) 
    {
      if (featdef == FEATDEF_CDS)
      {
        sfp = cds;
      }
      else if (featdef == FEATDEF_GENE)
      {
        sfp = GetGeneForFeature (cds);
      }
      else if (featdef == FEATDEF_mRNA)
      {
        sfp = SeqMgrGetOverlappingmRNA (cds->location, &fcontext);
      }
      if (sfp != NULL)
      {
        ValNodeAddPointer (&feat_list, OBJ_SEQFEAT, sfp);
      }
    }
  }
  return feat_list;
}


static ValNodePtr GetFeatureListForNucleotideBioseq (Uint1 featdef, BioseqPtr bsp)
{
  ValNodePtr feat_list = NULL;
  SeqFeatPtr sfp;
  SeqMgrFeatContext fcontext;
  Int4              seqfeattype;
  BioseqPtr         prot_bsp;

  if (bsp == NULL || ISA_aa (bsp->mol)) 
  {
    return NULL;
  }

  seqfeattype = FindFeatFromFeatDefType (featdef);
  if (seqfeattype == SEQFEAT_PROT)
  {
    for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, FEATDEF_CDS, &fcontext);
         sfp != NULL;
         sfp = SeqMgrGetNextFeature (bsp, sfp, 0, FEATDEF_CDS, &fcontext))
    {
      prot_bsp = BioseqFindFromSeqLoc (sfp->product);
      ValNodeLink (&feat_list, GetFeatureListForProteinBioseq (featdef, prot_bsp));
    }
  }
  else
  {
    for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, featdef, &fcontext);
         sfp != NULL;
         sfp = SeqMgrGetNextFeature (bsp, sfp, 0, featdef, &fcontext))
    {
      ValNodeAddPointer (&feat_list, OBJ_SEQFEAT, sfp);
    }
  }
  return feat_list;
}


static ValNodePtr GetFeaturesForGene (SeqFeatPtr gene, Uint1 featdef)
{
  BioseqPtr bsp;
  SeqFeatPtr sfp;
  ValNodePtr feat_list = NULL;
  SeqMgrFeatContext fcontext;
  Int4              start, stop, swap;

  if (gene == NULL) return NULL;

  bsp = BioseqFindFromSeqLoc (gene->location);
  start = SeqLocStart (gene->location);
  stop = SeqLocStop (gene->location);
  if (stop < start) 
  {
    swap = start;
    start = stop;
    stop = swap;
  }
  for (sfp = SeqMgrGetNextFeature (bsp, NULL, 0, featdef, &fcontext);
       sfp != NULL && fcontext.left < stop;
       sfp = SeqMgrGetNextFeature (bsp, sfp, 0, featdef, &fcontext))
  {
    if (fcontext.right >= start && gene == GetGeneForFeature (sfp))
    {
      ValNodeAddPointer (&feat_list, OBJ_SEQFEAT, sfp);
    }
  }
  return feat_list;
}


static ValNodePtr GetFeatureListForGene (Uint1 featdef, SeqFeatPtr gene)
{
  ValNodePtr feat_list = NULL, cds_list, vnp;
  SeqFeatPtr sfp, cds;
  SeqMgrFeatContext fcontext;
  BioseqPtr         protbsp;

  if (gene == NULL) 
  {
    return NULL;
  }

  if (featdef == FEATDEF_GENE)
  {
    ValNodeAddPointer (&feat_list, OBJ_SEQFEAT, gene);
  }
  else if (FindFeatFromFeatDefType (featdef == SEQFEAT_PROT))
  {
    cds_list = GetFeaturesForGene (gene, FEATDEF_CDS);
    for (vnp = cds_list; vnp != NULL; vnp = vnp->next) 
    {
      cds = vnp->data.ptrvalue;
      if (cds != NULL)
      {
        protbsp = BioseqFindFromSeqLoc (cds->product);
        for (sfp = SeqMgrGetNextFeature (protbsp, NULL, 0, featdef, &fcontext);
             sfp != NULL;
             sfp = SeqMgrGetNextFeature (protbsp, sfp, 0, featdef, &fcontext))
        {
          ValNodeAddPointer (&feat_list, OBJ_SEQFEAT, sfp);
        }
      }
    }
    cds_list = ValNodeFree (cds_list);
  }
  else
  {
    feat_list = GetFeaturesForGene (gene, featdef);
  }

  return feat_list;
}


static ValNodePtr AddFeaturesFromBioseqSet (BioseqSetPtr bssp, Uint1 featdef)
{
  SeqEntryPtr sep;
  BioseqPtr   bsp;
  Int4        seqfeattype;
  ValNodePtr  item_list = NULL;

  if (bssp == NULL) return NULL;

  seqfeattype = FindFeatFromFeatDefType (featdef);
  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
    if (sep->data.ptrvalue == NULL) continue;
    if (IS_Bioseq (sep)) {
      bsp = sep->data.ptrvalue;
      if (seqfeattype == SEQFEAT_PROT) {
        if (ISA_aa (bsp->mol)) {
          ValNodeLink (&item_list, GetFeatureListForProteinBioseq (featdef, bsp));
        }
      } else if (!ISA_aa (bsp->mol)) {
        ValNodeLink (&item_list, GetFeatureListForNucleotideBioseq (featdef, bsp));
      }
    } else if (IS_Bioseq_set (sep)) {
      ValNodeLink (&item_list, AddFeaturesFromBioseqSet (sep->data.ptrvalue, featdef));
    }
  }
  return item_list;
}


static ValNodePtr GetFeatureListForBioSourceObjects (ValNodePtr item_list, FeatureFieldPtr field)
{
  ValNodePtr vnp;
  SeqFeatPtr sfp;
  SeqDescrPtr sdp;
  BioseqPtr   bsp;
  ObjValNodePtr ovp;
  ValNodePtr  feature_list = NULL;

  if (item_list == NULL || field == NULL) return NULL;

  for (vnp = item_list; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == OBJ_SEQFEAT) {
      sfp = vnp->data.ptrvalue;
      if (sfp != NULL) {
        bsp = BioseqFindFromSeqLoc (sfp->location);
        ValNodeLink (&feature_list, GetFeatureListForNucleotideBioseq (GetFeatdefFromFeatureType(field->type), bsp));
      }
    } else if (vnp->choice == OBJ_SEQDESC) {
      sdp = vnp->data.ptrvalue;
      if (sdp != NULL && sdp->extended != 0) {
        ovp = (ObjValNodePtr) sdp;
        if (ovp->idx.parenttype == OBJ_BIOSEQSET) {
          ValNodeLink (&feature_list, AddFeaturesFromBioseqSet (ovp->idx.parentptr, GetFeatdefFromFeatureType(field->type)));
        } else if (ovp->idx.parenttype == OBJ_BIOSEQ) {
          bsp = (BioseqPtr) ovp->idx.parentptr;
          ValNodeLink (&feature_list, GetFeatureListForNucleotideBioseq (GetFeatdefFromFeatureType(field->type), bsp));
        }
      }
    }
  }
  return feature_list;
}


static ValNodePtr ValNodeCopyPtr (ValNodePtr orig)
{
  ValNodePtr new_list = NULL, last_vnp = NULL, vnp;

  while (orig != NULL) {
    vnp = ValNodeNew (NULL);
    vnp->choice = orig->choice;
    vnp->data.ptrvalue = orig->data.ptrvalue;
    if (last_vnp == NULL) {
      new_list = vnp;
    } else {
      last_vnp->next = vnp;
    }
    last_vnp = vnp;
    orig = orig->next;
  }
  return new_list;
}


static ValNodePtr GetFeatureListForRowAndColumn (Uint1 match_type, ValNodePtr match_list, FeatureFieldPtr field)
{
  ValNodePtr feature_list = NULL, vnp;

  if (match_list == NULL || field == NULL) return NULL;

  switch (match_type) {
    case eTableMatchFeatureID:
      feature_list = ValNodeCopyPtr (match_list);
      break;
    case eTableMatchGeneLocusTag:
      for (vnp = match_list; vnp != NULL; vnp = vnp->next) {
        ValNodeLink (&feature_list, GetFeatureListForGene (GetFeatdefFromFeatureType(field->type), vnp->data.ptrvalue));
      }
      break;
    case eTableMatchProteinID:
      for (vnp = match_list; vnp != NULL; vnp = vnp->next) {
        ValNodeLink (&feature_list, GetFeatureListForProteinBioseq (GetFeatdefFromFeatureType(field->type), vnp->data.ptrvalue));
      }
      break;
    case eTableMatchDbxref:
      feature_list = ValNodeCopyPtr (match_list);
      break;
    case eTableMatchNucID:
      for (vnp = match_list; vnp != NULL; vnp = vnp->next) {
        ValNodeLink (&feature_list, GetFeatureListForNucleotideBioseq (GetFeatdefFromFeatureType(field->type), vnp->data.ptrvalue));
      }
      break;
    case eTableMatchBioSource:
      ValNodeLink (&feature_list, GetFeatureListForBioSourceObjects (match_list, field));
      break;
  }
  return feature_list;
}


static void AddBioSourcesForBioseq (BioseqPtr bsp, ValNodePtr PNTR feature_list)
{
  SeqDescrPtr sdp;
  SeqMgrDescContext context;

  if (feature_list == NULL) return;
  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
        sdp != NULL;
        sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_source, &context)) {
    ValNodeAddPointer (feature_list, OBJ_SEQDESC, sdp);
  }
}

static void AddBioSourcesForFeature (SeqFeatPtr sfp, ValNodePtr PNTR feature_list)
{
  BioseqPtr bsp;

  if (sfp == NULL || feature_list == NULL) return;

  if (sfp->data.choice == SEQFEAT_BIOSRC) {
    ValNodeAddPointer (feature_list, OBJ_SEQFEAT, sfp);
  } else {
    bsp = BioseqFindFromSeqLoc (sfp->location);
    AddBioSourcesForBioseq (bsp, feature_list);
  }
}


static ValNodePtr GetBioSourceListForRowAndColumn (Uint1 match_type, ValNodePtr match_list, FeatureFieldPtr field)
{
  ValNodePtr feature_list = NULL, vnp;

  if (match_list == NULL || field == NULL) return NULL;

  switch (match_type) {
    case eTableMatchFeatureID:
      for (vnp = match_list; vnp != NULL; vnp = vnp->next) {
        if (vnp->choice == OBJ_SEQFEAT && vnp->data.ptrvalue != NULL) {
          AddBioSourcesForFeature (vnp->data.ptrvalue, &feature_list);
        }
      }
      break;
    case eTableMatchGeneLocusTag:
      for (vnp = match_list; vnp != NULL; vnp = vnp->next) {
        if (vnp->choice == OBJ_SEQFEAT && vnp->data.ptrvalue != NULL) {
          AddBioSourcesForFeature (vnp->data.ptrvalue, &feature_list);
        }
      }
      break;
    case eTableMatchProteinID:
    case eTableMatchNucID:
      for (vnp = match_list; vnp != NULL; vnp = vnp->next) {
        if (vnp->choice == OBJ_BIOSEQ) {
          AddBioSourcesForBioseq (vnp->data.ptrvalue, &feature_list);
        }
      }
      break;
    case eTableMatchDbxref:
      for (vnp = match_list; vnp != NULL; vnp = vnp->next) {
        if (vnp->choice == OBJ_SEQFEAT && vnp->data.ptrvalue != NULL) {
          AddBioSourcesForFeature (vnp->data.ptrvalue, &feature_list);
        }
      }
      break;
    case eTableMatchBioSource:
      feature_list = ValNodeCopyPtr (match_list);
      break;
  }
  return feature_list;
}


static ValNodePtr GetTargetListForRowAndColumn (Uint1 match_type, ValNodePtr match_list, FieldTypePtr field)
{
  ValNodePtr target_list = NULL;
  FeatureFieldPtr feature_field;

  if (field == NULL) return NULL;
  switch (field->choice) {
    case FieldType_source_qual:
      target_list = GetBioSourceListForRowAndColumn (match_type, match_list, field->data.ptrvalue);
      break;
    case FieldType_feature_field:
      target_list = GetFeatureListForRowAndColumn (match_type, match_list, field->data.ptrvalue);
      break;
    case FieldType_cds_gene_prot:
      feature_field = FeatureFieldFromCDSGeneProtField (field->data.intvalue);
      target_list = GetFeatureListForRowAndColumn (match_type, match_list, feature_field);
      feature_field = FeatureFieldFree (feature_field);
      break;
  }
  return target_list;
}


static void ReportMissingTargets (ValNodePtr PNTR perr_list, FieldTypePtr ft, CharPtr match_val, Int4 col_num, Int4 line_num)
{
  CharPtr            feat_name;
  FeatureFieldPtr    field;
  CharPtr            no_feat_fmt = "No %s feature for %s (column %d, line %d)";
  CharPtr            no_src_fmt = "No biosource for %s (column %d, line %d)";
  CharPtr            err_msg;

  if (perr_list == NULL || ft == NULL || match_val == NULL) return;

  switch (ft->choice) {
    case FieldType_source_qual:
      err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (no_feat_fmt) 
                                                    + StringLen (match_val)
                                                    + 30));
      sprintf (err_msg, no_src_fmt, match_val, col_num, line_num);
      ValNodeAddPointer (perr_list, 0, err_msg);
      break;
    case FieldType_feature_field:
      field = (FeatureFieldPtr) ft->data.ptrvalue;
      if (field != NULL) {
        feat_name = GetFeatureNameFromFeatureType (field->type);
        err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (no_feat_fmt) 
                                                      + StringLen (feat_name)
                                                      + StringLen (match_val)
                                                      + 30));
        sprintf (err_msg, no_feat_fmt, feat_name, match_val, col_num, line_num);
        ValNodeAddPointer (perr_list, 0, err_msg);
      }
      break;
    case FieldType_cds_gene_prot:
      field = FeatureFieldFromCDSGeneProtField (ft->data.intvalue);
      if (field != NULL) {
        feat_name = GetFeatureNameFromFeatureType (field->type);
        err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (no_feat_fmt) 
                                                      + StringLen (feat_name)
                                                      + StringLen (match_val)
                                                      + 30));
        sprintf (err_msg, no_feat_fmt, feat_name, match_val, col_num, line_num);
        ValNodeAddPointer (perr_list, 0, err_msg);
      }
      field = FeatureFieldFree (field);
      break;
  }
}


static void ReportEmptyIDColumn (ValNodePtr PNTR perr_list, Int4 line_num)
{
  CharPtr            err_msg;
  CharPtr            missing_id_fmt = "No ID for line %d";

  err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (missing_id_fmt) + 15));
  sprintf (err_msg, missing_id_fmt, line_num);
  ValNodeAddPointer (perr_list, 0, err_msg);
}

static ValNodePtr FindMatchChoiceInLine (ValNodePtr val_vnp, ValNodePtr col_vnp)
{
  TabColumnConfigPtr t;

  while (val_vnp != NULL && col_vnp != NULL) {
    t = (TabColumnConfigPtr) col_vnp->data.ptrvalue;
    if (t != NULL && t->match_type > 0) {
      val_vnp->choice = (Uint1) t->match_type;
      return val_vnp;
    }
    val_vnp = val_vnp->next;
    col_vnp = col_vnp->next;
  }
  return NULL;
}


NLM_EXTERN SeqFeatPtr GetmRNAForFeature (SeqFeatPtr sfp)
{
  SeqMgrFeatContext fcontext;
  BioseqPtr         pbsp;

  if (sfp == NULL) return NULL;
  if (sfp->data.choice == SEQFEAT_PROT) 
  { 
    pbsp = BioseqFindFromSeqLoc (sfp->location);
    sfp = SeqMgrGetCDSgivenProduct (pbsp, NULL);
    if (sfp == NULL) return NULL;
  }
  return SeqMgrGetOverlappingmRNA (sfp->location, &fcontext);
}


NLM_EXTERN Boolean AdjustmRNAProductToMatchProteinProduct (SeqFeatPtr sfp)
{
  SeqFeatPtr mrna;
  ProtRefPtr prp;
  RnaRefPtr  rrp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_PROT) return FALSE;

  prp = (ProtRefPtr) sfp->data.value.ptrvalue;
  mrna = GetmRNAForFeature (sfp);

  if (mrna == NULL) return FALSE;

  rrp = (RnaRefPtr) mrna->data.value.ptrvalue;
  if (rrp == NULL) 
  {
    rrp = RnaRefNew();
    mrna->data.value.ptrvalue = rrp;
  }

  rrp->ext.value.ptrvalue = MemFree (rrp->ext.value.ptrvalue);
  if (prp == NULL || prp->name == NULL || StringHasNoText (prp->name->data.ptrvalue))
  {
    rrp->ext.choice = 0;
  }
  else
  {
    rrp->ext.choice = 1;
    rrp->ext.value.ptrvalue = StringSave (prp->name->data.ptrvalue);
  }
  return TRUE;
}


NLM_EXTERN Boolean IsFieldTypeCDSProduct (FieldTypePtr ft)
{
  FeatureFieldPtr field;
  Boolean         rval = FALSE;

  if (ft == NULL) return FALSE;
  if (ft->choice == FieldType_feature_field) {
    field = (FeatureFieldPtr) ft->data.ptrvalue;
    if (field != NULL && field->type == Feature_type_cds
        && field->field != NULL
        && field->field->choice == FeatQualChoice_legal_qual
        && field->field->data.intvalue == Feat_qual_legal_product) {
      rval = TRUE;
    }
  } else if (ft->choice == FieldType_cds_gene_prot) {
    if (ft->data.intvalue == CDSGeneProt_field_prot_name) {
      rval = TRUE;
    }
  }
  return rval;
}


static Boolean IsFieldTypeGeneLocusTag (FieldTypePtr ft)
{
  FeatureFieldPtr field;
  Boolean         rval = FALSE;

  if (ft == NULL) return FALSE;
  if (ft->choice == FieldType_feature_field) {
    field = (FeatureFieldPtr) ft->data.ptrvalue;
    if (field != NULL && field->type == Feature_type_gene
        && field->field != NULL
        && field->field->choice == FeatQualChoice_legal_qual
        && field->field->data.intvalue == Feat_qual_legal_locus_tag) {
      rval = TRUE;
    }
  } else if (ft->choice == FieldType_cds_gene_prot) {
    if (ft->data.intvalue == CDSGeneProt_field_gene_locus_tag) {
      rval = TRUE;
    }
  }
  return rval;
}



NLM_EXTERN ValNodePtr ValidateTabTableValues (ValNodePtr table, ValNodePtr columns)
{
  ValNodePtr err_list = NULL;
  ValNodePtr line_vnp, col_vnp, val_vnp;
  Int4       line_num, col_num;
  TabColumnConfigPtr t;
  ValNodePtr locus_tag_values = NULL, bad_locus_tags = NULL, vnp;
  CharPtr    bad_format_fmt = "Locus tag %s has incorrect format";
  CharPtr    dup_fmt = "Locus tag %s appears in the table more than once";
  CharPtr    inconsistent_fmt = "Locus tag prefix for %s is inconsistent";
  CharPtr    err_msg;

  if (table == NULL || columns == NULL) {
    return NULL;
  }

  for (line_vnp = table, line_num = 1;
       line_vnp != NULL;
       line_vnp = line_vnp->next, line_num++) {
    for (val_vnp = line_vnp->data.ptrvalue, col_vnp = columns, col_num = 1;
         val_vnp != NULL && col_vnp != NULL;
         val_vnp = val_vnp->next, col_vnp = col_vnp->next, col_num++) {
      t = (TabColumnConfigPtr) col_vnp->data.ptrvalue;
      if (t == NULL || t->match_type > 0 || val_vnp == NULL || StringHasNoText (val_vnp->data.ptrvalue)) {
        continue;
      }
      if (IsFieldTypeGeneLocusTag (t->field)) {
        ValNodeAddPointer (&locus_tag_values, 0, val_vnp->data.ptrvalue);
      }
    }
  }

  bad_locus_tags = FindBadLocusTagsInList (locus_tag_values);
  for (vnp = bad_locus_tags; vnp != NULL; vnp = vnp->next) {
    switch (vnp->choice) {
      case eLocusTagErrorBadFormat:
        err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_format_fmt) + StringLen (vnp->data.ptrvalue)));
        sprintf (err_msg, bad_format_fmt, vnp->data.ptrvalue);
        ValNodeAddPointer (&err_list, 0, err_msg);
        break;
      case eLocusTagErrorDuplicate:
        err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (dup_fmt) + StringLen (vnp->data.ptrvalue)));
        sprintf (err_msg, dup_fmt, vnp->data.ptrvalue);
        ValNodeAddPointer (&err_list, 0, err_msg);
        break;
      case eLocusTagErrorInconsistentPrefix:
        err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (inconsistent_fmt) + StringLen (vnp->data.ptrvalue)));
        sprintf (err_msg, inconsistent_fmt, vnp->data.ptrvalue);
        ValNodeAddPointer (&err_list, 0, err_msg);
        break;
    }
  }
  locus_tag_values = ValNodeFree (locus_tag_values);
  return err_list;
}


NLM_EXTERN ValNodePtr GetObjectTableForTabTable (SeqEntryPtr sep, ValNodePtr table, ValNodePtr columns, ValNodePtr PNTR p_err_list)
{
  ValNodePtr err_list = NULL;
  ValNodePtr line_vnp, val_vnp, col_vnp;
  ValNodePtr obj_table = NULL, obj_row;
  Int4       line_num = 1, col_num;
  Uint2      entityID;
  ValNodePtr match_list, match_choice, target_list;
  TabColumnConfigPtr t;
  CharPtr            err_msg;
  CharPtr            no_match_fmt = "No match for %s, line %d";
  CharPtr            bad_col_val_fmt = "Did not set value for column %d, line %d";
  CharPtr            num_affected_fmt = "%d fields affected";
  Int4               num_fields_affected = 0;

  if (sep == NULL) {
    ValNodeAddPointer (&err_list, 0, StringSave ("No SeqEntry"));
  }
  if (table == NULL) {
    ValNodeAddPointer (&err_list, 0, StringSave ("No table"));
  }
  if (columns == NULL) {
    ValNodeAddPointer (&err_list, 0, StringSave ("No column information"));
  }
  if (err_list != NULL) {
    if (p_err_list == NULL) {
      err_list = ValNodeFreeData (err_list);
    } else {
      *p_err_list = err_list;
    }
    return NULL;
  }

  entityID = SeqMgrGetEntityIDForSeqEntry (sep);

  for (line_vnp = table, line_num = 1; line_vnp != NULL; line_vnp = line_vnp->next, line_num++) {
    obj_row = NULL;
    match_choice = FindMatchChoiceInLine (line_vnp->data.ptrvalue, columns);
    if (match_choice == NULL || StringHasNoText (match_choice->data.ptrvalue)) {
      ReportEmptyIDColumn (&err_list, line_num);
    } else {
      match_list = FindMatchForRow (match_choice, entityID, sep);
      if (match_list == NULL) {
        err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (no_match_fmt) + StringLen (match_choice->data.ptrvalue) + 15));
        sprintf (err_msg, no_match_fmt, match_choice->data.ptrvalue, line_num);
        ValNodeAddPointer (&err_list, 0, err_msg);
      } else {
        for (val_vnp = line_vnp->data.ptrvalue, col_vnp = columns, col_num = 1;
             col_vnp != NULL;
             col_vnp = col_vnp->next, col_num++) {
          target_list = NULL;
          t = (TabColumnConfigPtr) col_vnp->data.ptrvalue;
          if (t == NULL || t->match_type > 0 
              || (t->skip_blank && (val_vnp == NULL || StringHasNoText (val_vnp->data.ptrvalue)))) {
            /* no targets */
          } else {         
            target_list = GetTargetListForRowAndColumn (match_choice->choice, match_list, t->field);
            if (target_list == NULL) {
              ReportMissingTargets (&err_list, t->field, match_choice->data.ptrvalue, col_num, line_num); 
            }
          }
          ValNodeAddPointer (&obj_row, 0, target_list);
          if (val_vnp != NULL) {
            val_vnp = val_vnp->next;
          }
        }
      }
    }
    ValNodeAddPointer (&obj_table, 0, obj_row);
  }

  if (err_list != NULL) {
    if (p_err_list == NULL) {
      err_list = ValNodeFreeData (err_list);
    } else {
      *p_err_list = err_list;
    }
  }  
  return obj_table;
}


NLM_EXTERN ValNodePtr FreeObjectTableForTabTable (ValNodePtr table)
{
  ValNodePtr vnp_next, vnp_row, vnp_row_next;

  while (table != NULL) {
    vnp_next = table->next;
    table->next = NULL;
    vnp_row = table->data.ptrvalue;
    while (vnp_row != NULL) {
      vnp_row_next = vnp_row->next;
      vnp_row->next = NULL;
      vnp_row->data.ptrvalue = ValNodeFree (vnp_row->data.ptrvalue);
      vnp_row = ValNodeFree (vnp_row);
      vnp_row = vnp_row_next;
    }
    table = ValNodeFree (table);
    table = vnp_next;
  }
  return table;
}


typedef struct countfeat {
  Uint1 featdef;
  Int4 num;
} CountFeatData, PNTR CountFeatPtr;


static void CountFeaturesCallback (SeqFeatPtr sfp, Pointer userdata)
{
  CountFeatPtr p;

  if (sfp == NULL || userdata == NULL) return;

  p = (CountFeatPtr) userdata;
  if (sfp->idx.subtype == p->featdef) {
    p->num++;
  }
}

static void CountBioSourceDescriptorsCallback (SeqDescrPtr sdp, Pointer userdata)
{
  Int4Ptr p;

  p = (Int4Ptr) userdata;
  if (sdp != NULL && p != NULL && sdp->choice == Seq_descr_source) {
    (*p)++;
  }
}


static ValNodePtr CountObjectsForColumnFields (SeqEntryPtr sep, ValNodePtr columns)
{
  ValNodePtr count_list = NULL, vnp;
  TabColumnConfigPtr t;
  CountFeatData d;
  FeatureFieldPtr f;
  Int4 num;
  Uint1 featdef = 0;

  d.featdef = 0;
  d.num = 0;
  for (vnp = columns; vnp != NULL; vnp = vnp->next) {
    num = 0;
    t = (TabColumnConfigPtr) vnp->data.ptrvalue;
    if (t != NULL && t->match_type == 0 && t->field != NULL) {
      switch (t->field->choice) {
        case FieldType_source_qual:
          if (featdef != FEATDEF_BIOSRC) {
            d.featdef = FEATDEF_BIOSRC;
            d.num = 0;
            VisitFeaturesInSep (sep, &d, CountFeaturesCallback);
            VisitDescriptorsInSep (sep, &(d.num), CountBioSourceDescriptorsCallback);
          }
          num = d.num;
          break;
        case FieldType_feature_field:
          f = (FeatureFieldPtr) t->field->data.ptrvalue;
          if (f != NULL) {
            featdef = GetFeatdefFromFeatureType(f->type);
            if (featdef != d.featdef) {
              d.featdef = featdef;
              d.num = 0;
              VisitFeaturesInSep (sep, &d, CountFeaturesCallback);
            }
            num = d.num;
          }
          break;
        case FieldType_cds_gene_prot:
          f = FeatureFieldFromCDSGeneProtField (t->field->data.intvalue);
          if (f != NULL) {
            featdef = GetFeatdefFromFeatureType(f->type);
            if (featdef != d.featdef) {
              d.featdef = featdef;
              d.num = 0;
              VisitFeaturesInSep (sep, &d, CountFeaturesCallback);
            }
            num = d.num;
          }
          f = FeatureFieldFree (f);
          break;
      }
    }
    ValNodeAddInt (&count_list, 0, num);
  }
  return count_list;
}


NLM_EXTERN ValNodePtr ApplyTableValuesToObjectTable (SeqEntryPtr sep, ValNodePtr table, ValNodePtr columns, ValNodePtr obj_table)
{
  ValNodePtr val_line_vnp, obj_line_vnp;
  ValNodePtr val_vnp, obj_vnp, col_vnp;
  ValNodePtr target_vnp;
  TabColumnConfigPtr t;
  CharPtr val, qual_name;
  ValNodePtr         err_list = NULL, count_list, count_affected_list = NULL, count_vnp, count_tot_vnp;
  CharPtr            err_msg;
  CharPtr            bad_col_val_fmt = "Did not set value for column %d, line %d";
  CharPtr            num_affected_fmt = "%d fields affected";
  CharPtr            col_num_affected_fmt = "For %s (column %d), %d items were affected out of %d total";
  Int4 num_fields_affected = 0, col_num, line_num, num_this_column;
  Boolean success;
  ValNodePtr count_msg = NULL;

  count_list = CountObjectsForColumnFields (sep, columns);

  for (val_line_vnp = table, obj_line_vnp = obj_table, line_num = 1;
       val_line_vnp != NULL && obj_line_vnp != NULL;
       val_line_vnp = val_line_vnp->next, obj_line_vnp = obj_line_vnp->next, line_num++) {
    val_vnp = val_line_vnp->data.ptrvalue;
    obj_vnp = obj_line_vnp->data.ptrvalue;
    col_vnp = columns;
    col_num = 1;
    count_vnp = count_affected_list;
    while (obj_vnp != NULL && col_vnp != NULL) {
      num_this_column = 0;
      if (obj_vnp->data.ptrvalue != NULL) {
        t = (TabColumnConfigPtr) col_vnp->data.ptrvalue;
        if (t == NULL || t->match_type > 0 
            || (t->skip_blank && (val_vnp == NULL || StringHasNoText (val_vnp->data.ptrvalue)))) {
          /* ignore column or skip blank value */
        } else {
          if (val_vnp == NULL || val_vnp->data.ptrvalue == NULL) {
            val = "";
          } else {
            val = val_vnp->data.ptrvalue;
          }
          for (target_vnp = obj_vnp->data.ptrvalue; target_vnp != NULL; target_vnp = target_vnp->next) {
            if (val[0] == 0) {
              success = RemoveFieldValueForObject (target_vnp->choice, target_vnp->data.ptrvalue, t->field, NULL);
            } else {
              success = SetFieldValueForObject (target_vnp->choice, target_vnp->data.ptrvalue, t->field, NULL,
                                                val_vnp->data.ptrvalue, t->existing_text);
            }
            if (success) {
              num_fields_affected++;
              num_this_column++;
              if (t->match_mrna && IsFieldTypeCDSProduct (t->field)
                  && target_vnp->choice == OBJ_SEQFEAT) {
                if (AdjustmRNAProductToMatchProteinProduct (target_vnp->data.ptrvalue)) {
                  num_fields_affected++;
                }
              }
            } else {
              err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_col_val_fmt) + 30));
              sprintf (err_msg, bad_col_val_fmt, col_num, line_num);
              ValNodeAddPointer (&err_list, 0, err_msg);
            }
          }
        }
      }
      if (val_vnp != NULL) {
        val_vnp = val_vnp->next;
      }
      if (count_vnp == NULL) {
        ValNodeAddInt (&count_affected_list, 0, num_this_column);
      } else {
        count_vnp->data.intvalue ++;
        count_vnp = count_vnp->next;
      }
      obj_vnp = obj_vnp->next;
      col_vnp = col_vnp->next;
      col_num++;
    }
  }

  /* put message at top of list for number of fields affected */
  err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (num_affected_fmt) + 15));
  sprintf (err_msg, num_affected_fmt, num_fields_affected);
  ValNodeAddPointer (&count_msg, 0, err_msg);

  /* if any affected, list number of fields per column, and the total in the record */
  if (num_fields_affected > 0) {
    for (count_vnp = count_affected_list, count_tot_vnp = count_list, col_vnp = columns, col_num = 1;
         count_vnp != NULL && count_tot_vnp != NULL && col_vnp != NULL;
         count_vnp = count_vnp->next, count_tot_vnp = count_tot_vnp->next, col_vnp = col_vnp->next, col_num++) {
      t = (TabColumnConfigPtr) col_vnp->data.ptrvalue;
      if (t != NULL && t->match_type == 0) {
        qual_name = SummarizeFieldType (t->field);
        err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (col_num_affected_fmt) + StringLen (qual_name) + 45));
        sprintf (err_msg, col_num_affected_fmt, qual_name, col_num, count_vnp->data.intvalue, count_tot_vnp->data.intvalue);      
        ValNodeAddPointer (&count_msg, 0, err_msg);
        qual_name = MemFree (qual_name);
      }
    }
  }

  ValNodeLink (&count_msg, err_list);

  count_list = ValNodeFree (count_list);
  count_affected_list = ValNodeFree (count_affected_list);

  return count_msg;
}


static int LIBCALLBACK SortVnpByChoiceAndPtrvalue (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      if (vnp1->choice > vnp2->choice) {
        return 1;
      } else if (vnp1->choice < vnp2->choice) {
        return -1;
      } else if (vnp1->data.ptrvalue > vnp2->data.ptrvalue) {
        return 1;
      } else if (vnp2->data.ptrvalue < vnp2->data.ptrvalue) {
        return -1;
      } else {
        return 0;
      }
    }
  }
  return 0;
}


NLM_EXTERN ValNodePtr CheckObjTableForRowsThatApplyToTheSameDestination (ValNodePtr obj_table)
{
  Int4 col_num;
  ValNodePtr line_vnp, col_vnp, obj_vnp, vnp;
  ValNodePtr col_list = NULL, col_obj_list;
  Boolean any_column_values_left;
  ValNodePtr err_list = NULL;
  Boolean found_multi;
  CharPtr multi_fmt = "Multiple rows apply to the same object for column %d";
  CharPtr err_msg;
  
  /* now, for each row, get pointer to first column */
  for (line_vnp = obj_table; line_vnp != NULL; line_vnp = line_vnp->next) {
    ValNodeAddPointer (&col_list, 0, line_vnp->data.ptrvalue);
  }

  /* now for each column, make a list of all features in the column, then sort to see if there are duplicates */
  any_column_values_left = TRUE;
  col_num = 1;
  while (any_column_values_left) {
    any_column_values_left = FALSE;
    col_obj_list = NULL;
    for (vnp = col_list; vnp != NULL; vnp = vnp->next) {
      col_vnp = vnp->data.ptrvalue;
      if (col_vnp != NULL) {
        obj_vnp = col_vnp->data.ptrvalue;
        ValNodeLink (&col_obj_list, ValNodeCopyPtr (obj_vnp));
        vnp->data.ptrvalue = col_vnp->next;
        any_column_values_left = TRUE;
      }
    }
    if (col_obj_list != NULL) {
      found_multi = FALSE;
      col_obj_list = ValNodeSort (col_obj_list, SortVnpByChoiceAndPtrvalue);
      for (vnp = col_obj_list; vnp != NULL && vnp->next != NULL && !found_multi; vnp = vnp->next) {
        if (vnp->choice == vnp->next->choice
            && vnp->data.ptrvalue == vnp->next->data.ptrvalue) {
          found_multi = TRUE;
        }
      }
      if (found_multi) {
        err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (multi_fmt)
                                                      + 30));
        sprintf (err_msg, multi_fmt, col_num);
        ValNodeAddPointer (&err_list, col_num, err_msg);
      }
      col_obj_list = ValNodeFree (col_obj_list);
    }
    col_num++;
  }
  col_list = ValNodeFree (col_list);
  return err_list;
}


NLM_EXTERN ValNodePtr CheckObjTableForExistingText (SeqEntryPtr sep, ValNodePtr table, ValNodePtr columns, ValNodePtr obj_table)
{
  ValNodePtr err_list = NULL, vnp;
  ValNodePtr val_line_vnp, obj_line_vnp;
  ValNodePtr val_vnp, obj_vnp, col_vnp;
  Int4       line_num = 1, col_num, num_existing_text = 0;
  Uint2      entityID;
  TabColumnConfigPtr t;
  CharPtr            err_msg, str, qual_name, val;
  CharPtr            already_has_val_fmt = "%s already has value '%s' (column %d), line %d.  Replacement is '%s'";
  CharPtr            num_existing_text_fmt = "%d fields already have text.";
  CharPtr            mrna_warn_fmt = "%d coding region features have mRNAs, but %d do not.";
  ValNodePtr         target_list, feat_vnp;
  Int4               num_with_mrna = 0, num_without_mrna = 0;

  if (sep == NULL) {
    ValNodeAddPointer (&err_list, 1, StringSave ("No SeqEntry"));
  }
  if (table == NULL) {
    ValNodeAddPointer (&err_list, 1, StringSave ("No table"));
  }
  if (columns == NULL) {
    ValNodeAddPointer (&err_list, 1, StringSave ("No column information"));
  }
  if (err_list != NULL) {
    return err_list;
  }

  entityID = SeqMgrGetEntityIDForSeqEntry (sep);

  for (val_line_vnp = table, obj_line_vnp = obj_table, line_num = 1;
       val_line_vnp != NULL && obj_line_vnp != NULL;
       val_line_vnp = val_line_vnp->next, obj_line_vnp = obj_line_vnp->next, line_num++) {
    val_vnp = val_line_vnp->data.ptrvalue;
    obj_vnp = obj_line_vnp->data.ptrvalue;
    col_vnp = columns;
    if (val_vnp == NULL || obj_vnp == NULL) continue;
    col_num = 1;
    while (obj_vnp != NULL && col_vnp != NULL) {
      if (obj_vnp->data.ptrvalue != NULL) {
        t = (TabColumnConfigPtr) col_vnp->data.ptrvalue;
        if (t == NULL || t->match_type > 0 
            || (t->skip_blank && (val_vnp == NULL || StringHasNoText (val_vnp->data.ptrvalue)))) {
          /* ignore column or skip blank value */
        } else {
          target_list = obj_vnp->data.ptrvalue;
          if (val_vnp == NULL || val_vnp->data.ptrvalue == NULL) {
            val = "";
          } else {
            val = val_vnp->data.ptrvalue;
          }
          for (feat_vnp = target_list; feat_vnp != NULL; feat_vnp = feat_vnp->next) {
            /* check for existing text */
            str = GetFieldValueForObject (feat_vnp->choice, feat_vnp->data.ptrvalue, t->field, NULL);
            if (!StringHasNoText (str)) {
              qual_name = SummarizeFieldType (t->field);
              err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (already_has_val_fmt)
                                                          + StringLen (qual_name) + StringLen (str)  
                                                          + StringLen (val)
                                                          + 30));
              sprintf (err_msg, already_has_val_fmt, qual_name, str, col_num, line_num, val);
              ValNodeAddPointer (&err_list, col_num, err_msg);
              num_existing_text ++;
            }
            str = MemFree (str);
            /* check for mrna if changing CDS product */
            if (IsFieldTypeCDSProduct (t->field) && feat_vnp->choice == OBJ_SEQFEAT) {
              if (GetmRNAForFeature (feat_vnp->data.ptrvalue) != NULL) {
                num_with_mrna++;
              } else {
                num_without_mrna++;
              }
            }
          }
        }
      }
      if (val_vnp != NULL) {
        val_vnp = val_vnp->next;
      }
      obj_vnp = obj_vnp->next;
      col_vnp = col_vnp->next;
      col_num++;
    }
  }          
  if (num_existing_text > 0) {
    err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (num_existing_text_fmt)
                                                + 15));
    sprintf (err_msg, num_existing_text_fmt, num_existing_text);
    vnp = ValNodeNew (NULL);
    vnp->choice = 0;
    vnp->data.ptrvalue = err_msg;
    vnp->next = err_list;
    err_list = vnp;
  }
  if (num_with_mrna > 0 && num_without_mrna > 0) {
    err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (mrna_warn_fmt)
                                                + 30));
    sprintf (err_msg, mrna_warn_fmt, num_with_mrna, num_without_mrna);
    vnp = ValNodeNew (NULL);
    vnp->choice = 0;
    vnp->data.ptrvalue = err_msg;
    vnp->next = err_list;
    err_list = vnp;
  }    

  return err_list;
}


NLM_EXTERN ValNodePtr ApplyTableToFeatures (SeqEntryPtr sep, ValNodePtr table, ValNodePtr columns)
{
  ValNodePtr err_list = NULL;
  ValNodePtr line_vnp, val_vnp, col_vnp;
  Int4       line_num = 1, col_num;
  Uint2      entityID;
  ValNodePtr match_list, match_choice, target_list, feat_vnp, vnp;
  TabColumnConfigPtr t;
  CharPtr            err_msg;
  CharPtr            no_match_fmt = "No match for %s, line %d";
  CharPtr            bad_col_val_fmt = "Did not set value for column %d, line %d";
  CharPtr            num_affected_fmt = "%d fields affected";
  Int4               num_fields_affected = 0;
  CharPtr            val;
  Boolean            success;

  if (sep == NULL) {
    ValNodeAddPointer (&err_list, 0, StringSave ("No SeqEntry"));
  }
  if (table == NULL) {
    ValNodeAddPointer (&err_list, 0, StringSave ("No table"));
  }
  if (columns == NULL) {
    ValNodeAddPointer (&err_list, 0, StringSave ("No column information"));
  }
  if (err_list != NULL) {
    return err_list;
  }

  entityID = SeqMgrGetEntityIDForSeqEntry (sep);

  for (line_vnp = table, line_num = 1; line_vnp != NULL; line_vnp = line_vnp->next, line_num++) {
    match_choice = FindMatchChoiceInLine (line_vnp->data.ptrvalue, columns);
    if (match_choice == NULL || StringHasNoText (match_choice->data.ptrvalue)) {
      ReportEmptyIDColumn (&err_list, line_num);
    } else {
      match_list = FindMatchForRow (match_choice, entityID, sep);
      if (match_list == NULL) {
        err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (no_match_fmt) + StringLen (match_choice->data.ptrvalue) + 15));
        sprintf (err_msg, no_match_fmt, match_choice->data.ptrvalue, line_num);
        ValNodeAddPointer (&err_list, 0, err_msg);
      } else {
        for (val_vnp = line_vnp->data.ptrvalue, col_vnp = columns, col_num = 1;
             col_vnp != NULL;
             col_vnp = col_vnp->next, col_num++) {
          t = (TabColumnConfigPtr) col_vnp->data.ptrvalue;
          if (t == NULL || t->match_type > 0 
              || (t->skip_blank && (val_vnp == NULL || StringHasNoText (val_vnp->data.ptrvalue)))) {
            if (val_vnp != NULL) {
              val_vnp = val_vnp->next;
            }            
            continue;
          }
          
          target_list = GetTargetListForRowAndColumn (match_choice->choice, match_list, t->field);
          if (target_list == NULL) {
            ReportMissingTargets (&err_list, t->field, match_choice->data.ptrvalue, col_num, line_num); 
          } else {
            if (val_vnp == NULL || val_vnp->data.ptrvalue == NULL) {
              val = "";
            } else {
              val = val_vnp->data.ptrvalue;
            }
            for (feat_vnp = target_list; feat_vnp != NULL; feat_vnp = feat_vnp->next) {
              if (val[0] == 0) {
                success = RemoveFieldValueForObject (feat_vnp->choice, feat_vnp->data.ptrvalue, t->field, NULL);
              } else {
                success = SetFieldValueForObject (feat_vnp->choice, feat_vnp->data.ptrvalue, t->field, NULL,
                                                  val_vnp->data.ptrvalue, t->existing_text);
              }
              if (success) {
                num_fields_affected++;
                if (t->match_mrna && IsFieldTypeCDSProduct (t->field)
                    && feat_vnp->choice == OBJ_SEQFEAT) {
                  if (AdjustmRNAProductToMatchProteinProduct (feat_vnp->data.ptrvalue)) {
                    num_fields_affected++;
                  }
                }
              } else {
                err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (bad_col_val_fmt) + 30));
                sprintf (err_msg, bad_col_val_fmt, col_num, line_num);
                ValNodeAddPointer (&err_list, 0, err_msg);
              }
            }
          }
          target_list = ValNodeFree (target_list);
          if (val_vnp != NULL) {
            val_vnp = val_vnp->next;
          }
        }
      }
      match_list = ValNodeFree (match_list);
    }
  }
  
  err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (num_affected_fmt) + 15));
  sprintf (err_msg, num_affected_fmt, num_fields_affected);
  vnp = ValNodeNew (NULL);
  vnp->data.ptrvalue = err_msg;
  vnp->next = err_list;
  err_list = vnp;

  return err_list;
}


NLM_EXTERN ValNodePtr CheckTableForExistingText (SeqEntryPtr sep, ValNodePtr table, ValNodePtr columns)
{
  ValNodePtr err_list = NULL, vnp;
  ValNodePtr line_vnp, val_vnp, col_vnp;
  Int4       line_num = 1, col_num, num_existing_text = 0;
  Uint2      entityID;
  TabColumnConfigPtr t;
  CharPtr            err_msg, str, qual_name, val;
  CharPtr            no_match_fmt = "No match for %s, line %d";
  CharPtr            no_feat_fmt = "No %s feature for %s (column %d, line %d)";
  CharPtr            already_has_val_fmt = "%s already has value '%s' (column %d), line %d.  Replacement is '%s'";
  CharPtr            num_existing_text_fmt = "%d fields already have text.";
  ValNodePtr         match_choice, match_list;
  ValNodePtr         target_list, feat_vnp;

  if (sep == NULL) {
    ValNodeAddPointer (&err_list, 1, StringSave ("No SeqEntry"));
  }
  if (table == NULL) {
    ValNodeAddPointer (&err_list, 1, StringSave ("No table"));
  }
  if (columns == NULL) {
    ValNodeAddPointer (&err_list, 1, StringSave ("No column information"));
  }
  if (err_list != NULL) {
    return err_list;
  }

  entityID = SeqMgrGetEntityIDForSeqEntry (sep);

  for (line_vnp = table, line_num = 1; line_vnp != NULL; line_vnp = line_vnp->next, line_num++) {
    match_choice = FindMatchChoiceInLine (line_vnp->data.ptrvalue, columns);
    if (match_choice == NULL || StringHasNoText (match_choice->data.ptrvalue)) {
      ReportEmptyIDColumn (&err_list, line_num);
    } else {
      match_list = FindMatchForRow (match_choice, entityID, sep);
      if (match_list == NULL) {
        err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (no_match_fmt) + StringLen (match_choice->data.ptrvalue) + 15));
        sprintf (err_msg, no_match_fmt, match_choice->data.ptrvalue, line_num);
        ValNodeAddPointer (&err_list, 0, err_msg);
      } else {
        for (val_vnp = line_vnp->data.ptrvalue, col_vnp = columns, col_num = 1;
             col_vnp != NULL;
             col_vnp = col_vnp->next, col_num++) {
          t = (TabColumnConfigPtr) col_vnp->data.ptrvalue;
          if (t == NULL || t->match_type > 0 
              || (t->skip_blank && (val_vnp == NULL || StringHasNoText (val_vnp->data.ptrvalue)))) {
            if (val_vnp != NULL) {
              val_vnp = val_vnp->next;
            }            
            continue;
          }
          target_list = GetTargetListForRowAndColumn (match_choice->choice, match_list, t->field);
          if (target_list == NULL) {
            ReportMissingTargets (&err_list, t->field, match_choice->data.ptrvalue, col_num, line_num); 
          } else {
            if (val_vnp == NULL || val_vnp->data.ptrvalue == NULL) {
              val = "";
            } else {
              val = val_vnp->data.ptrvalue;
            }
            for (feat_vnp = target_list; feat_vnp != NULL; feat_vnp = feat_vnp->next) {
              str = GetFieldValueForObject (feat_vnp->choice, feat_vnp->data.ptrvalue, t->field, NULL);
              if (!StringHasNoText (str)) {
                qual_name = SummarizeFieldType (t->field);
                err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (already_has_val_fmt)
                                                            + StringLen (qual_name) + StringLen (str)  
                                                            + StringLen (val)
                                                            + 30));
                sprintf (err_msg, already_has_val_fmt, qual_name, str, col_num, line_num, val);
                ValNodeAddPointer (&err_list, col_num, err_msg);
                num_existing_text ++;
              }
              str = MemFree (str);
            }
          }
          target_list = ValNodeFree (target_list);
          if (val_vnp != NULL) {
            val_vnp = val_vnp->next;
          }
        }
      }
      match_list = ValNodeFree (match_list);
    }
  }          
  if (num_existing_text > 0) {
    err_msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (num_existing_text_fmt)
                                                + 15));
    sprintf (err_msg, num_existing_text_fmt, num_existing_text);
    vnp = ValNodeNew (NULL);
    vnp->choice = 0;
    vnp->data.ptrvalue = err_msg;
    vnp->next = err_list;
    err_list = vnp;
  }

  return err_list;
}





