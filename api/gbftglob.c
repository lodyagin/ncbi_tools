/**************************************************************************
*   gbftglob.c:
*   -- all the globally variables for gbfeat.c
*   -- all the defined variables in the gbfeat.h
*
* $Log: gbftglob.c,v $
* Revision 6.2  1998/01/15 20:28:08  tatiana
* increased the size of optional quals array for source feature
*
* Revision 6.1  1997/12/23 22:01:07  tatiana
* focus and specimen_voucher
*
* Revision 6.0  1997/08/25 18:06:02  madden
* Revision changed to 6.0
*
* Revision 5.10  1997/07/29 20:59:54  vakatov
* Encapsulated 'ParFlat_GBQual_names' and 'ParFlat_GBFeat'(formerly
* global) variables into access functions. Made other global variables
* be "extern" instead of "NLM_EXTERN"(i.e. local within the PC DLL).
*
* Revision 5.9  1997/07/11 17:00:58  tatiana
* in misc_recomb GBQUAL_organism changed from mandatory ro optional
*
 * Revision 5.8  1997/01/08  21:08:56  tatiana
 * /clone added to misc_difference
 *
 * Revision 5.7  1996/09/17  14:52:17  tatiana
 * 'virion' added
 *
 * Revision 5.6  1996/08/02  16:50:45  tatiana
 * a typo fixed
 *
 * Revision 5.5  1996/08/01  17:15:47  kans
 * fixed a typo
 *
 * Revision 5.4  1996/07/30  17:28:07  kans
 * ParFlat_... arrays now external in header file
 *
 * Revision 5.3  1996/07/29  19:45:34  tatiana
 * GBQual_names changed to use a structure
 *
 * Revision 5.2  1996/07/25  19:35:34  kans
 * ParFlat_GBQual_class array was missing an item at the cultivar position
 *
 * Revision 5.1  1996/07/25  14:17:34  tatiana
 * added qualifiers: allele, exception replace
 *
 * Revision 4.5  1995/11/15  18:03:32  tatiana
 * a bug fixed.
 *
 * Revision 4.4  1995/11/13  15:53:50  tatiana
 * serotype added
 *
 * Revision 4.3  1995/11/08  22:55:24  tatiana
 * /serotype added
 *
 * Revision 4.2  1995/08/16  22:02:08  tatiana
 * changes for db_xref
 *
 * Revision 4.1  1995/08/15  22:07:16  tatiana
 * db_xref added
 *
 * ..
 *
 * Revision 1.7  1995/05/15  21:46:05  ostell
 * added Log line
 *
*
*                                                                 10-14-93
***************************************************************************/

#include <stdio.h>
#include <ncbi.h>
#include <gbftdef.h>

static GbFeatName STATIC__ParFlat_GBQual_names[ParFlat_TOTAL_GBQUAL] = { 
 {"allele", Class_text}, {"anticodon", Class_pos_aa}, 
 {"bound_moiety", Class_text}, {"cell_line", Class_text}, 
 {"cell_type", Class_text}, {"chromosome", Class_text}, 
 {"chloroplast", Class_none}, {"chromoplast", Class_none}, 
 {"citation", Class_bracket_int}, {"clone", Class_text}, 
 {"clone_lib", Class_text}, {"codon", Class_seq_aa}, 
 {"codon_start", Class_int_or}, {"cons_splice", Class_site}, 
 {"cultivar", Class_text}, {"cyanelle", Class_none},
 {"db_xref", Class_text}, {"dev_stage", Class_text},
 {"direction", Class_L_R_B}, {"EC_number", Class_ecnum}, 
 {"evidence", Class_exper}, {"exception", Class_text},
 {"frequency", Class_text}, {"function", Class_text}, 
 {"gene", Class_text}, {"gdb_xref", Class_text}, 
 {"germline", Class_none}, {"haplotype", Class_text},
 {"insertion_seq", Class_text}, {"isolate", Class_text}, 
 {"kinetoplast", Class_none}, {"label", Class_token}, 
 {"lab_host", Class_text}, {"map", Class_text}, 
 {"macronuclear", Class_none}, {"mitochondrion", Class_none}, 
 {"mod_base", Class_token}, {"note", Class_note},
 {"number", Class_int}, {"organism", Class_text}, 
 {"partial", Class_none}, {"PCR_conditions", Class_text}, 
 {"pop_variant", Class_text}, {"phenotype", Class_text},
  {"plasmid", Class_text}, {"product", Class_text}, 
 {"proviral", Class_none}, {"pseudo", Class_none},
 {"rearranged", Class_none}, { "replace", Class_text}, 
 {"rpt_family", Class_text}, {"rpt_type", Class_rpt},
 { "rpt_unit", Class_token}, { "sex", Class_text},
 {"sequenced_mol", Class_text}, { "serotype", Class_text},
 {"specific_host", Class_text}, {"standard_name", Class_text}, 
 {"strain", Class_text}, {"sub_clone", Class_text}, 
 {"sub_species", Class_text}, {"sub_strain", Class_text}, 
 {"tissue_lib", Class_text},  {"tissue_type", Class_text}, 
 {"translation", Class_text}, {"transl_except", Class_pos_aa}, 
 {"transl_table", Class_int}, {"transposon", Class_text}, 
 {"usedin", Class_token}, {"variety", Class_text}, {"virion", Class_none},
 { "specimen_voucher", Class_text}, 
 {"focus", Class_none}
 };

NLM_EXTERN GbFeatNamePtr x_ParFlat_GBQual_names(void) {
  return STATIC__ParFlat_GBQual_names;
}

CharPtr ParFlat_IntOrString[ParFlat_TOTAL_IntOr] = {"1", "2", "3"};

CharPtr ParFlat_LRBString[ParFlat_TOTAL_LRB] = {"LEFT", "RIGHT", "BOTH"};

CharPtr ParFlat_ExpString[ParFlat_TOTAL_Exp] = {
                    "EXPERIMENTAL", "NOT_EXPERIMENTAL"};

CharPtr ParFlat_RptString[ParFlat_TOTAL_Rpt] = {
       "tandem", "inverted", "flanking", "terminal", "direct",
       "dispersed", "other"};

static SematicFeat STATIC__ParFlat_GBFeat[ParFlat_TOTAL_GBFEAT] = {
   {"allele", 0, {-1, -1, -1, -1, -1}, 13,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_frequency, GBQUAL_gene,
      GBQUAL_label, GBQUAL_map, GBQUAL_note, GBQUAL_partial, 
      GBQUAL_phenotype, GBQUAL_product, GBQUAL_replace, GBQUAL_standard_name, 
      GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"attenuator", 0, {-1, -1, -1, -1, -1}, 10,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_label, 
     GBQUAL_gene, GBQUAL_map,
      GBQUAL_note, GBQUAL_partial, GBQUAL_phenotype, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"C_region", 0, {-1, -1, -1, -1, -1}, 12,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, 
     GBQUAL_gene, GBQUAL_label, GBQUAL_map, GBQUAL_note, GBQUAL_partial, 
     GBQUAL_product, GBQUAL_pseudo, GBQUAL_standard_name, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"CAAT_signal", 0, {-1, -1, -1, -1, -1}, 9,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_label, 
     GBQUAL_gene, GBQUAL_map, GBQUAL_note, GBQUAL_partial, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"CDS", 0, {-1, -1, -1, -1, -1}, 22,
     {GBQUAL_citation, GBQUAL_codon, GBQUAL_codon_start, GBQUAL_db_xref, 
     GBQUAL_EC_number,GBQUAL_evidence, GBQUAL_exception, GBQUAL_function, 
     GBQUAL_gdb_xref, GBQUAL_gene, GBQUAL_label, GBQUAL_map, GBQUAL_note, 
     GBQUAL_number, GBQUAL_partial, GBQUAL_product, GBQUAL_pseudo, 
     GBQUAL_standard_name, GBQUAL_translation, GBQUAL_transl_except, 
     GBQUAL_transl_table, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"conflict",  1, {GBQUAL_citation, -1, -1, -1, -1}, 6,
     {GBQUAL_db_xref, GBQUAL_map, GBQUAL_note, GBQUAL_gene, GBQUAL_usedin,
     GBQUAL_replace, 
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"D-loop",  0, {-1, -1, -1, -1, -1}, 8,
     {GBQUAL_citation, GBQUAL_label, GBQUAL_gene, GBQUAL_map, GBQUAL_note,
      GBQUAL_partial, GBQUAL_usedin, GBQUAL_db_xref, 
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"D_segment", 0, {-1, -1, -1, -1, -1}, 12,
     {GBQUAL_citation, GBQUAL_evidence, GBQUAL_gene, 
      GBQUAL_label, GBQUAL_map, GBQUAL_note, GBQUAL_partial, GBQUAL_product,
      GBQUAL_pseudo, GBQUAL_standard_name, GBQUAL_usedin, GBQUAL_db_xref, 
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"enhancer",  0, {-1, -1, -1, -1, -1}, 10,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_label, 
     GBQUAL_gene, GBQUAL_map, 
      GBQUAL_note, GBQUAL_partial, GBQUAL_standard_name,  GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"exon", 0, {-1, -1, -1, -1, -1}, 15,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_EC_number,
      GBQUAL_evidence, GBQUAL_function, GBQUAL_gene, GBQUAL_label,
      GBQUAL_map, GBQUAL_note, GBQUAL_number, GBQUAL_partial,
      GBQUAL_product, GBQUAL_pseudo, GBQUAL_standard_name, 
      GBQUAL_usedin, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"GC_signal", 0, {-1, -1, -1, -1, -1}, 9,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_label, 
     GBQUAL_gene, GBQUAL_map, GBQUAL_note, GBQUAL_partial, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"gene", 1, {GBQUAL_gene, -1, -1, -1, -1}, 12,
     {GBQUAL_allele, GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, 
     GBQUAL_function, GBQUAL_label, GBQUAL_map, GBQUAL_note, GBQUAL_partial, 
     GBQUAL_pseudo, GBQUAL_phenotype, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
	{"iDNA", 0, {-1, -1, -1, -1, -1}, 12,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_function, 
     GBQUAL_label,
      GBQUAL_gene, GBQUAL_map, GBQUAL_note, GBQUAL_number, GBQUAL_partial,
      GBQUAL_standard_name, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"intron", 0, {-1, -1, -1, -1, -1}, 13,
     {GBQUAL_citation, GBQUAL_cons_splice, GBQUAL_db_xref, GBQUAL_evidence, 
     GBQUAL_function, 
      GBQUAL_gene, GBQUAL_label, GBQUAL_map, GBQUAL_note, GBQUAL_number, 
      GBQUAL_partial, GBQUAL_standard_name, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"J_segment", 0, {-1, -1, -1, -1, -1}, 12,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_gene, 
      GBQUAL_label, GBQUAL_map, GBQUAL_note, GBQUAL_partial, GBQUAL_product, 
      GBQUAL_pseudo, GBQUAL_standard_name, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"LTR", 0, {-1, -1, -1, -1, -1}, 10,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_function, 
     GBQUAL_gene, 
      GBQUAL_label, GBQUAL_note, GBQUAL_partial, GBQUAL_standard_name,
      GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"mat_peptide",  0, {-1, -1, -1, -1, -1}, 14,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_EC_number,
      GBQUAL_evidence, GBQUAL_function, GBQUAL_gene, GBQUAL_label,
      GBQUAL_map, GBQUAL_note, GBQUAL_partial, GBQUAL_pseudo, 
      GBQUAL_product, GBQUAL_standard_name, GBQUAL_usedin,
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"misc_binding",  1, {GBQUAL_bound_moiety, -1, -1, -1, -1}, 10,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_function, 
     GBQUAL_gene,
      GBQUAL_label, GBQUAL_map, GBQUAL_note, GBQUAL_partial,
      GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"misc_difference", 0, {-1, -1, -1, -1, -1}, 11,
     {GBQUAL_citation, GBQUAL_clone, GBQUAL_db_xref, GBQUAL_gene, 
     GBQUAL_label, GBQUAL_map, GBQUAL_partial, GBQUAL_replace, 
      GBQUAL_note, GBQUAL_standard_name, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"misc_feature",  0, {-1, -1, -1, -1, -1}, 15,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_function, 
     GBQUAL_gene,
      GBQUAL_label, GBQUAL_map, GBQUAL_note, GBQUAL_partial, GBQUAL_number, 
      GBQUAL_phenotype, GBQUAL_product, GBQUAL_pseudo, GBQUAL_standard_name,
      GBQUAL_usedin, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"misc_recomb", 0, {-1, -1, -1, -1, -1}, 11,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_gene, 
     GBQUAL_label, GBQUAL_map, GBQUAL_note, GBQUAL_organism, 
     GBQUAL_partial, GBQUAL_standard_name, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"misc_RNA",  0, {-1, -1, -1, -1, -1}, 12,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_function, 
     GBQUAL_gene,
      GBQUAL_label, GBQUAL_map, GBQUAL_note, GBQUAL_partial, GBQUAL_product, 
      GBQUAL_standard_name, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"misc_signal",  0, {-1, -1, -1, -1, -1}, 12,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_function, 
     GBQUAL_gene,
      GBQUAL_label, GBQUAL_map, GBQUAL_note, GBQUAL_partial, GBQUAL_phenotype, 
      GBQUAL_standard_name, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"misc_structure",  0, {-1, -1, -1, -1, -1}, 11,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_function, 
     GBQUAL_gene,
      GBQUAL_label, GBQUAL_map, GBQUAL_note, GBQUAL_partial,
      GBQUAL_standard_name, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"modified_base",  1, {GBQUAL_mod_base, -1, -1, -1, -1}, 9,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_frequency, 
     GBQUAL_gene,
      GBQUAL_label, GBQUAL_map, GBQUAL_note, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"mRNA",  0, {-1, -1, -1, -1, -1}, 13,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_function, 
     GBQUAL_gene, GBQUAL_label, GBQUAL_map, GBQUAL_note, GBQUAL_partial, 	
     GBQUAL_product, GBQUAL_pseudo, GBQUAL_standard_name, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"mutation",  0, {-1, -1, -1, -1, -1}, 12,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_frequency, GBQUAL_gene, 
     GBQUAL_label, GBQUAL_map,
      GBQUAL_note, GBQUAL_phenotype, GBQUAL_product, GBQUAL_replace, 
      GBQUAL_standard_name, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"N_region",  0, {-1, -1, -1, -1, -1}, 11,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_gene,
      GBQUAL_label, GBQUAL_map, GBQUAL_note, GBQUAL_product, 
      GBQUAL_pseudo, GBQUAL_standard_name, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"old_sequence",   1, {GBQUAL_citation, -1, -1, -1, -1}, 7,
     {GBQUAL_db_xref, GBQUAL_gene, GBQUAL_map, GBQUAL_note, GBQUAL_partial,
      GBQUAL_replace, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"polyA_signal", 0, {-1, -1, -1, -1, -1}, 9,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_gene,
      GBQUAL_label, GBQUAL_map, 
      GBQUAL_note, GBQUAL_partial, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"polyA_site", 0, {-1, -1, -1, -1, -1}, 8,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_gene,
      GBQUAL_label, GBQUAL_map, GBQUAL_note, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"precursor_RNA",  0, {-1, -1, -1, -1, -1}, 11,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_function,
      GBQUAL_gene, GBQUAL_label, GBQUAL_map, GBQUAL_note, GBQUAL_partial,
      GBQUAL_standard_name, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"prim_transcript",  0, {-1, -1, -1, -1, -1}, 11,
     {GBQUAL_citation, GBQUAL_evidence, GBQUAL_function, GBQUAL_gene,
      GBQUAL_label, GBQUAL_map, GBQUAL_note, GBQUAL_partial,
      GBQUAL_standard_name, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"primer_bind",  0, {-1, -1, -1, -1, -1}, 11,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_gene,
      GBQUAL_label, GBQUAL_map,
      GBQUAL_note, GBQUAL_partial, GBQUAL_standard_name, GBQUAL_PCR_conditions,
      GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"promoter",  0, {-1, -1, -1, -1, -1}, 13,
     {GBQUAL_citation,GBQUAL_db_xref,  GBQUAL_evidence, GBQUAL_function,
      GBQUAL_gene,
      GBQUAL_label, GBQUAL_map, GBQUAL_note, GBQUAL_partial, GBQUAL_phenotype,
      GBQUAL_pseudo, GBQUAL_standard_name, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"protein_bind",  1, {GBQUAL_bound_moiety, -1, -1, -1, -1}, 11,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_function,
      GBQUAL_gene, 
      GBQUAL_label, GBQUAL_map, GBQUAL_note, GBQUAL_partial, 
      GBQUAL_standard_name, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"RBS",  0, {-1, -1, -1, -1, -1}, 10,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_gene,
      GBQUAL_label, GBQUAL_map,
      GBQUAL_note, GBQUAL_partial, GBQUAL_standard_name, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"repeat_region",  0, {-1, -1, -1, -1, -1}, 14,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_function,
      GBQUAL_gene, 
      GBQUAL_label, GBQUAL_map, GBQUAL_note, GBQUAL_partial, GBQUAL_rpt_type,
      GBQUAL_rpt_family, GBQUAL_rpt_unit, GBQUAL_standard_name, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"repeat_unit",  0, {-1, -1, -1, -1, -1}, 12,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_function,
      GBQUAL_gene, 
      GBQUAL_label, GBQUAL_map, GBQUAL_note, GBQUAL_partial, GBQUAL_rpt_family,
      GBQUAL_rpt_type, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"rep_origin",  0, {-1, -1, -1, -1, -1}, 11,
     {GBQUAL_citation, GBQUAL_direction, GBQUAL_db_xref, GBQUAL_evidence,
      GBQUAL_gene,
      GBQUAL_label, GBQUAL_map, GBQUAL_note, GBQUAL_partial,
      GBQUAL_standard_name, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"rRNA",  0, {-1, -1, -1, -1, -1}, 13,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_function, 
     GBQUAL_gene,
      GBQUAL_label, GBQUAL_map, GBQUAL_note, GBQUAL_partial, GBQUAL_product,
      GBQUAL_pseudo, GBQUAL_standard_name, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"S_region",  0, {-1, -1, -1, -1, -1}, 12,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_gene,
      GBQUAL_label, GBQUAL_map, GBQUAL_note, GBQUAL_partial, GBQUAL_product,
      GBQUAL_pseudo, GBQUAL_standard_name, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"satellite",  0, {-1, -1, -1, -1, -1}, 13,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_gene,
      GBQUAL_label, GBQUAL_map,
      GBQUAL_note, GBQUAL_partial, GBQUAL_rpt_type, GBQUAL_rpt_family,
      GBQUAL_rpt_unit, GBQUAL_standard_name, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"scRNA",  0, {-1, -1, -1, -1, -1}, 13,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_function,
      GBQUAL_gene,
      GBQUAL_label, GBQUAL_map, GBQUAL_note, GBQUAL_partial, GBQUAL_product,
      GBQUAL_pseudo, GBQUAL_standard_name, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"sig_peptide", 0, {-1, -1, -1, -1, -1}, 13,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, 
      GBQUAL_function, GBQUAL_gene, GBQUAL_label, GBQUAL_map, GBQUAL_note,
      GBQUAL_partial, GBQUAL_product, GBQUAL_pseudo, 
      GBQUAL_standard_name, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"snRNA",  0, {-1, -1, -1, -1, -1}, 13,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_function,
      GBQUAL_gene,
      GBQUAL_label, GBQUAL_map, GBQUAL_note, GBQUAL_partial, GBQUAL_product,
      GBQUAL_pseudo, GBQUAL_standard_name, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"source", 1, {GBQUAL_organism, -1, -1, -1, -1}, 43,
     {GBQUAL_cell_line, GBQUAL_cell_type, GBQUAL_chloroplast, 
      GBQUAL_chromoplast, GBQUAL_chromosome,
      GBQUAL_citation, GBQUAL_clone, GBQUAL_clone_lib, GBQUAL_cultivar,
      GBQUAL_cyanelle, GBQUAL_db_xref, GBQUAL_dev_stage, GBQUAL_focus,
      GBQUAL_germline, GBQUAL_haplotype,
      GBQUAL_lab_host, GBQUAL_insertion_seq, GBQUAL_isolate, GBQUAL_kinetoplast,
      GBQUAL_label, GBQUAL_macronuclear, GBQUAL_map, GBQUAL_mitochondrion,
      GBQUAL_note, GBQUAL_plasmid, GBQUAL_pop_variant, 
      GBQUAL_proviral, GBQUAL_rearranged, GBQUAL_sex, GBQUAL_sequenced_mol, 
      GBQUAL_serotype, GBQUAL_specific_host, GBQUAL_strain, GBQUAL_sub_clone, 
      GBQUAL_sub_species, GBQUAL_sub_strain, GBQUAL_tissue_lib, 
      GBQUAL_tissue_type, GBQUAL_transposon, GBQUAL_usedin, 
      GBQUAL_specimen_voucher, GBQUAL_variety, GBQUAL_virion,
      -1, -1, -1, -1, -1, -1, -1}},
   {"stem_loop",  0, {-1, -1, -1, -1, -1}, 11,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_function,
      GBQUAL_gene,
      GBQUAL_label, GBQUAL_map, GBQUAL_note, GBQUAL_partial, 
      GBQUAL_standard_name, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"STS",  0, {-1, -1, -1, -1, -1}, 9,
     {GBQUAL_citation, GBQUAL_standard_name, GBQUAL_db_xref, GBQUAL_gene,
      GBQUAL_label, GBQUAL_usedin, GBQUAL_note, GBQUAL_partial, GBQUAL_map, 
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"TATA_signal",  0, {-1, -1, -1, -1, -1}, 9,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_gene,
      GBQUAL_label, GBQUAL_map,
      GBQUAL_note, GBQUAL_partial, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"terminator",  0, {-1, -1, -1, -1, -1}, 10,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_gene,
      GBQUAL_label, GBQUAL_map,GBQUAL_note, GBQUAL_partial,
       GBQUAL_standard_name, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"transit_peptide",  0, {-1, -1, -1, -1, -1}, 13,
     {GBQUAL_citation, GBQUAL_db_xref,
      GBQUAL_evidence, GBQUAL_function, GBQUAL_gene, GBQUAL_label, GBQUAL_map,
       GBQUAL_note,GBQUAL_partial, GBQUAL_product, GBQUAL_pseudo, 	
       GBQUAL_standard_name, GBQUAL_usedin,
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"tRNA",  0, {-1, -1, -1, -1, -1}, 14,
     {GBQUAL_anticodon, GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence,
      GBQUAL_function, GBQUAL_gene, GBQUAL_label, GBQUAL_map, GBQUAL_note, 
      GBQUAL_partial, GBQUAL_product, GBQUAL_pseudo, GBQUAL_standard_name,
       GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"unsure",  0, {-1, -1, -1, -1, -1}, 8,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_gene, GBQUAL_usedin, 
     GBQUAL_label, GBQUAL_map, GBQUAL_note, GBQUAL_replace, 
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"V_region",  0, {-1, -1, -1, -1, -1}, 12,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence,
      GBQUAL_gene, GBQUAL_label, GBQUAL_map, GBQUAL_note, GBQUAL_partial, 
      GBQUAL_product, GBQUAL_pseudo, GBQUAL_standard_name, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"V_segment",  0, {-1, -1, -1, -1, -1}, 12,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence,
      GBQUAL_gene, GBQUAL_label, GBQUAL_map, GBQUAL_note, GBQUAL_partial, 
      GBQUAL_product, GBQUAL_pseudo, GBQUAL_standard_name, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"variation",  0, {-1, -1, -1, -1, -1}, 13,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_frequency, GBQUAL_gene,
      GBQUAL_label, GBQUAL_map, GBQUAL_note, GBQUAL_partial, GBQUAL_replace,
      GBQUAL_phenotype, GBQUAL_product, GBQUAL_standard_name, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"virion",  0, {-1, -1, -1, -1, -1}, 9,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_gene, GBQUAL_label, GBQUAL_map,
      GBQUAL_note, GBQUAL_organism, GBQUAL_partial, GBQUAL_usedin, 
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"3'clip",  0, {-1, -1, -1, -1, -1}, 11,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_function,
      GBQUAL_gene, GBQUAL_label, GBQUAL_map, GBQUAL_note, GBQUAL_partial, 
      GBQUAL_standard_name, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"3'UTR",  0, {-1, -1, -1, -1, -1}, 11,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_function,
      GBQUAL_gene, GBQUAL_label, GBQUAL_map, GBQUAL_note, GBQUAL_partial, 
      GBQUAL_standard_name, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"5'clip",  0, {-1, -1, -1, -1, -1}, 11,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_function, GBQUAL_gene,
      GBQUAL_evidence, GBQUAL_label, GBQUAL_map,  
      GBQUAL_note, GBQUAL_partial, GBQUAL_standard_name, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"5'UTR",  0, {-1, -1, -1, -1, -1}, 11,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_function,
      GBQUAL_gene, GBQUAL_label, GBQUAL_map, 
      GBQUAL_note, GBQUAL_partial, GBQUAL_standard_name, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"-10_signal",  0, {-1, -1, -1, -1, -1}, 10,
     {GBQUAL_citation,GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_gene,
     GBQUAL_label, GBQUAL_map,GBQUAL_note,GBQUAL_partial, 
     GBQUAL_standard_name, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
   {"-35_signal",  0, {-1, -1, -1, -1, -1}, 9,
     {GBQUAL_citation, GBQUAL_db_xref, GBQUAL_evidence, GBQUAL_gene,
      GBQUAL_label, GBQUAL_map, GBQUAL_note, GBQUAL_partial, GBQUAL_usedin,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}}                      
      };

NLM_EXTERN SematicFeatPtr x_ParFlat_GBFeat(void) {
  return STATIC__ParFlat_GBFeat;
}


