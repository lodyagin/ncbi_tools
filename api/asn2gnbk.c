/*   asn2gnbk.c
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
* File Name:  asn2gnbk.c
*
* Author:  Karl Sirotkin, Tom Madden, Tatiana Tatusov, Jonathan Kans,
*          Mati Shomrat
*
* Version Creation Date:   10/21/98
*
* $Revision: 6.592 $
*
* File Description:  New GenBank flatfile generator - work in progress
*
* Modifications:
* --------------------------------------------------------------------------
* ==========================================================================
*/

#include <ncbi.h>
#include <objall.h>
#include <objsset.h>
#include <objsub.h>
#include <objfdef.h>
#include <seqport.h>
#include <sequtil.h>
#include <sqnutils.h>
#include <subutil.h>
#include <tofasta.h>
#include <explore.h>
#include <gbfeat.h>
#include <gbftdef.h>
#include <edutil.h>
#include <alignmgr2.h>
#include <asn2gnbp.h>

#ifdef WIN_MAC
#if __profile__
#include <Profiler.h>
#endif
#endif

#define ASN2FF_EMBL_MAX 78
#define ASN2FF_GB_MAX   79
#define SEQID_MAX_LEN   41

#define TILDE_IGNORE     0
#define TILDE_TO_SPACES  1
#define TILDE_EXPAND     2
#define TILDE_OLD_EXPAND 3


/* flags set by mode to customize behavior */

typedef struct asn2gbflags {
  Boolean             suppressLocalID;
  Boolean             validateFeats;
  Boolean             ignorePatPubs;
  Boolean             dropShortAA;
  Boolean             avoidLocusColl;
  Boolean             iupacaaOnly;
  Boolean             dropBadCitGens;
  Boolean             noAffilOnUnpub;
  Boolean             dropIllegalQuals;
  Boolean             checkQualSyntax;
  Boolean             needRequiredQuals;
  Boolean             needOrganismQual;
  Boolean             needAtLeastOneRef;
  Boolean             citArtIsoJta;
  Boolean             dropBadDbxref;
  Boolean             useEmblMolType;
  Boolean             hideBankItComment;
  Boolean             checkCDSproductID;
  Boolean             suppressSegLoc;
  Boolean             srcQualsToNote;
  Boolean             hideEmptySource;
  Boolean             forGbRelease;
} Asn2gbFlags, PNTR Asn2gbFlagsPtr;

/* internal Asn2gbSect structure has fields on top of Asn2gbSect fields */

typedef struct int_Asn2gbSect {
  Asn2gbSect  asp;
  SeqPortPtr  spp;
} IntAsn2gbSect, PNTR IntAsn2gbSectPtr;

/* string structure */

#define STRING_BUF_LEN  1024

typedef struct stringitem {
  struct stringitem  *curr;
  struct stringitem  *next;
  Pointer            iajp;
  Char               buf [STRING_BUF_LEN];
  Int4               pos;
} StringItem, PNTR StringItemPtr;

/* internal asn2gbjob structure has fields on top of Asn2gbJob fields */

typedef struct int_asn2gb_job {
  Asn2gbJob      ajp;
  FmtType        format;
  Asn2gbFlags    flags;
  Boolean        showFarTransl;
  Boolean        transIfNoProd;
  Boolean        alwaysTranslCds;
  Boolean        transientSeqPort;
  Boolean        masterStyle;
  Boolean        newSourceOrg;
  ValNodePtr     lockedBspList;
  Boolean        relModeError;
  IndxPtr        index;
  GBSeqPtr       gbseq;
  StringItemPtr  pool;
  Boolean        www;
} IntAsn2gbJob, PNTR IntAsn2gbJobPtr;

/* structure for storing working parameters while building asn2gb_job structure */

typedef struct asn2gbwork {
  IntAsn2gbJobPtr  ajp;
  Uint2            entityID;

  FmtType          format;
  ModType          mode;
  StlType          style;

  ValNodePtr       pubhead;    /* for collecting publications */
  ValNodePtr       srchead;    /* for collecting biosources */

  /* linked lists of paragraphs, sections, blocks */

  ValNodePtr       sectionList;
  ValNodePtr       blockList;    /* reset for each new section */

  /* most recent node of linked lists, for quickly adding next node */

  ValNodePtr       lastsection;
  ValNodePtr       lastblock;    /* reset for each new section */

  Int4             currsection;

  /* section fields needed for populating blocks */

  Asn2gbSectPtr    asp;

  BioseqPtr        target;
  BioseqPtr        parent;
  BioseqPtr        bsp;
  BioseqPtr        refs;
  SeqLocPtr        slp;
  Uint2            seg;
  Int4             numsegs;
  Int4             partcount;
  Int4             from;
  Int4             to;
  Boolean          showAllFeats;

  Boolean          contig;
  Boolean          showconfeats;
  Boolean          showconsource;

  Boolean          onlyNearFeats;
  Boolean          farFeatsSuppress;

  Boolean          citSubsFirst;
  Boolean          hideGeneFeats;
  Boolean          newLocusLine;

  Boolean          hideImpFeats;
  Boolean          hideSnpFeats;

  Boolean          isGPS;
  Boolean          copyGpsCdsUp;
  Boolean          copyGpsGeneDown;

  Char             basename [SEQID_MAX_LEN];

  SeqFeatPtr       lastsfp;
  SeqAnnotPtr      lastsap;
  Int4             lastleft;
  Int4             lastright;

  Boolean          firstfeat;
  Boolean          featseen;

  SeqSubmitPtr     ssp;
  Boolean          hup;

  Boolean          failed;
} Asn2gbWork, PNTR Asn2gbWorkPtr;

/* array for assigning biosource and feature data fields to qualifiers */
/* should be allocated to MAX (ASN2GNBK_TOTAL_SOURCE, ASN2GNBK_TOTAL_FEATUR) */

typedef union qualval {
  CharPtr        str;
  Boolean        ble;
  Int4           num;
  ValNodePtr     vnp;
  GBQualPtr      gbq;
  OrgModPtr      omp;
  SubSourcePtr   ssp;
  CodeBreakPtr   cbp;
  SeqLocPtr      slp;
  SeqIdPtr       sip;
  tRNAPtr        trp;
  UserObjectPtr  uop;
} QualVal, PNTR QualValPtr;

/* structure passed to individual paragraph format functions */

typedef struct asn2gbformat {
  IntAsn2gbJobPtr ajp;
  Asn2gbSectPtr   asp;
  QualValPtr      qvp;

  FmtType         format;
} Asn2gbFormat, PNTR Asn2gbFormatPtr;


/* Seq-hist replacedBy is preformatted into string field, */
/* then comment descriptors, Map location:, and Region:, */
/* then comment features, finally HTGS */

typedef struct comment_block {
  ASN2GB_BASE_BLOCK
  Boolean           first;
} CommentBlock, PNTR CommentBlockPtr;

/* internal reference block has fields on top of RefBlock fields */

typedef struct int_ref_block {
  RefBlock   rb;
  DatePtr    date;     /* internal sorting use only */
  SeqLocPtr  loc;      /* final location on target bioseq */
  CharPtr    authstr;  /* author string */
  Uint2      index;    /* index if feature on target bioseq */
  Boolean    justuids; /* gibb pub with uids and Figure, etc. */
  CharPtr    fig;      /* figure string from equivalent gibb pub */
  CharPtr    maploc;   /* maploc string from equivalent gibb pub */
  Boolean    poly_a;   /* poly_a field from equivalent gibb pub */
} IntRefBlock, PNTR IntRefBlockPtr;

/* internal source block has fields on top of BaseBlock fields */

typedef struct int_src_block {
  BaseBlock         bb;
  Boolean           is_descriptor;
  Boolean           is_focus;
  Boolean           is_synthetic;
  BioSourcePtr      biop;
  Uint4             orghash;
  Uint4             modhash;
  Uint4             subhash;
  Uint4             xrfhash;
  SeqLocPtr         loc;     /* final location on target bioseq */
  CharPtr           taxname;
  CharPtr           comment;
  OrgModPtr         omp;
  SubSourcePtr      ssp;
  ValNodePtr        vnp;
  Int4              left;
  Int4              right;
} IntSrcBlock, PNTR IntSrcBlockPtr;

/* internal feature block has fields on top of FeatBlock fields */

typedef struct int_feat_block {
  FeatBlock  fb;
  Boolean     mapToNuc;
  Boolean     mapToProt;
  Boolean     mapToGen;
  Boolean     mapToMrna;
  Boolean     mapToPep;
  Boolean     isCDS;     /* set if using IntCdsBlock */
  Boolean     firstfeat;
} IntFeatBlock, PNTR IntFeatBlockPtr;

/* internal cds block has fields on top of IntFeatBlock fields */

typedef struct int_cds_block {
  IntFeatBlock  ifb;
  FeatBlock     fb;
  CharPtr       fig;    /* figure string from pub */
  CharPtr       maploc; /* maploc string from pub */
} IntCdsBlock, PNTR IntCdsBlockPtr;


/* enumerated qualifier category definitions */

typedef enum {
  Qual_class_ignore = 0,
  Qual_class_string,
  Qual_class_tilde,
  Qual_class_sgml,
  Qual_class_boolean,
  Qual_class_int,
  Qual_class_evidence,
  Qual_class_valnode,
  Qual_class_EC_valnode,
  Qual_class_xtraprds,
  Qual_class_quote,
  Qual_class_EC_quote,
  Qual_class_noquote,
  Qual_class_label,
  Qual_class_number,
  Qual_class_paren,
  Qual_class_region,
  Qual_class_replace,
  Qual_class_consplice,
  Qual_class_bond,
  Qual_class_site,
  Qual_class_L_R_B,
  Qual_class_rpt,
  Qual_class_organelle,
  Qual_class_orgmod,
  Qual_class_subsource,
  Qual_class_code_break,
  Qual_class_anti_codon,
  Qual_class_codon,
  Qual_class_method,
  Qual_class_pubset,
  Qual_class_db_xref,
  Qual_class_seq_id,
  Qual_class_seq_loc,
  Qual_class_its,
  Qual_class_sec_str,
  Qual_class_trna_codons,
  Qual_class_translation,
  Qual_class_protnames,
  Qual_class_illegal,
  Qual_class_note,
  Qual_class_rpt_unit,
  Qual_class_product,
  Qual_class_model_ev,
  Qual_class_gene_syn
}  QualType;

/* source 'feature' */

/* some qualifiers will require additional content verification not
   explicitly indicated by the class type */

typedef enum {
  SCQUAL_acronym = 1,
  SCQUAL_anamorph,
  SCQUAL_authority,
  SCQUAL_biotype,
  SCQUAL_biovar,
  SCQUAL_breed,
  SCQUAL_cell_line,
  SCQUAL_cell_type,
  SCQUAL_chemovar,
  SCQUAL_chromosome,
  SCQUAL_citation,
  SCQUAL_clone,
  SCQUAL_clone_lib,
  SCQUAL_common,
  SCQUAL_common_name,
  SCQUAL_country,
  SCQUAL_cultivar,
  SCQUAL_db_xref,
  SCQUAL_org_xref,
  SCQUAL_dev_stage,
  SCQUAL_dosage,
  SCQUAL_ecotype,
  SCQUAL_endogenous_virus_name,
  SCQUAL_environmental_sample,
  SCQUAL_extrachrom,
  SCQUAL_focus,
  SCQUAL_forma,
  SCQUAL_forma_specialis,
  SCQUAL_frequency,
  SCQUAL_gb_acronym,
  SCQUAL_gb_anamorph,
  SCQUAL_gb_synonym,
  SCQUAL_genotype,
  SCQUAL_germline,
  SCQUAL_group,
  SCQUAL_haplotype,
  SCQUAL_ins_seq_name,
  SCQUAL_isolate,
  SCQUAL_isolation_source,
  SCQUAL_lab_host,
  SCQUAL_label,
  SCQUAL_macronuclear,
  SCQUAL_map,
  SCQUAL_mol_type,
  SCQUAL_note,
  SCQUAL_old_lineage,
  SCQUAL_old_name,
  SCQUAL_organism,
  SCQUAL_organelle,
  SCQUAL_orgmod_note,
  SCQUAL_pathovar,
  SCQUAL_plasmid_name,
  SCQUAL_plastid_name,
  SCQUAL_pop_variant,
  SCQUAL_rearranged,
  SCQUAL_segment,
  SCQUAL_seqfeat_note,
  SCQUAL_sequenced_mol,
  SCQUAL_serogroup,
  SCQUAL_serotype,
  SCQUAL_serovar,
  SCQUAL_sex,
  SCQUAL_spec_or_nat_host,
  SCQUAL_specimen_voucher,
  SCQUAL_strain,
  SCQUAL_sub_clone,
  SCQUAL_sub_group,
  SCQUAL_sub_species,
  SCQUAL_sub_strain,
  SCQUAL_sub_type,
  SCQUAL_subsource_note,
  SCQUAL_synonym,
  SCQUAL_teleomorph,
  SCQUAL_tissue_lib,
  SCQUAL_tissue_type,
  SCQUAL_transgenic,
  SCQUAL_transposon_name,
  SCQUAL_type,
  SCQUAL_unstructured,
  SCQUAL_usedin,
  SCQUAL_variety,
  SCQUAL_zero_orgmod,
  SCQUAL_one_orgmod,
  SCQUAL_zero_subsrc,
  ASN2GNBK_TOTAL_SOURCE
}  SourceType;

/* ordering arrays for qualifiers and note components */

static Uint1 source_qual_order [] = {
  SCQUAL_organism,

  SCQUAL_organelle,

  SCQUAL_strain,
  SCQUAL_sub_strain,
  SCQUAL_variety,
  SCQUAL_serotype,
  SCQUAL_serovar,
  SCQUAL_cultivar,
  SCQUAL_isolate,
  SCQUAL_isolation_source,
  SCQUAL_spec_or_nat_host,
  SCQUAL_sub_species,
  SCQUAL_specimen_voucher,

  SCQUAL_db_xref,
  SCQUAL_org_xref,

  SCQUAL_chromosome,
  SCQUAL_map,
  SCQUAL_clone,
  SCQUAL_sub_clone,
  SCQUAL_haplotype,
  SCQUAL_sex,
  SCQUAL_cell_line,
  SCQUAL_cell_type,
  SCQUAL_tissue_type,
  SCQUAL_clone_lib,
  SCQUAL_dev_stage,
  SCQUAL_frequency,

  SCQUAL_germline,
  SCQUAL_rearranged,
  SCQUAL_transgenic,
  SCQUAL_environmental_sample,

  SCQUAL_lab_host,
  SCQUAL_pop_variant,
  SCQUAL_tissue_lib,

  SCQUAL_plasmid_name,
  SCQUAL_transposon_name,
  SCQUAL_ins_seq_name,

  SCQUAL_country,

  SCQUAL_focus,

  SCQUAL_note,

  SCQUAL_sequenced_mol,
  SCQUAL_label,
  SCQUAL_usedin,
  SCQUAL_citation,
  0
};

static Uint1 source_desc_note_order [] = {
  SCQUAL_seqfeat_note,
  SCQUAL_orgmod_note,
  SCQUAL_subsource_note,

  SCQUAL_type,
  SCQUAL_sub_type,
  SCQUAL_serogroup,
  SCQUAL_pathovar,
  SCQUAL_chemovar,
  SCQUAL_biovar,
  SCQUAL_biotype,
  SCQUAL_group,
  SCQUAL_sub_group,
  SCQUAL_common,
  SCQUAL_acronym,
  SCQUAL_dosage,

  SCQUAL_authority,
  SCQUAL_forma,
  SCQUAL_forma_specialis,
  SCQUAL_ecotype,
  SCQUAL_synonym,
  SCQUAL_anamorph,
  SCQUAL_teleomorph,
  SCQUAL_breed,

  SCQUAL_genotype,
  SCQUAL_plastid_name,

  SCQUAL_segment,
  SCQUAL_endogenous_virus_name,

  SCQUAL_common_name,

  SCQUAL_zero_orgmod,
  SCQUAL_one_orgmod,
  SCQUAL_zero_subsrc,

  /* SCQUAL_old_lineage, */

  /* SCQUAL_old_name, */
  0
};

static Uint1 source_feat_note_order [] = {
  SCQUAL_unstructured,

  SCQUAL_type,
  SCQUAL_sub_type,
  SCQUAL_serogroup,
  SCQUAL_pathovar,
  SCQUAL_chemovar,
  SCQUAL_biovar,
  SCQUAL_biotype,
  SCQUAL_group,
  SCQUAL_sub_group,
  SCQUAL_common,
  SCQUAL_acronym,
  SCQUAL_dosage,

  SCQUAL_authority,
  SCQUAL_forma,
  SCQUAL_forma_specialis,
  SCQUAL_ecotype,
  SCQUAL_synonym,
  SCQUAL_anamorph,
  SCQUAL_teleomorph,
  SCQUAL_breed,

  SCQUAL_genotype,
  SCQUAL_plastid_name,

  SCQUAL_segment,
  SCQUAL_endogenous_virus_name,

  SCQUAL_seqfeat_note,
  SCQUAL_orgmod_note,
  SCQUAL_subsource_note,

  SCQUAL_common_name,

  SCQUAL_zero_orgmod,
  SCQUAL_one_orgmod,
  SCQUAL_zero_subsrc,

  /* SCQUAL_old_lineage, */

  /* SCQUAL_old_name, */
  0
};

typedef struct sourcequal {
  CharPtr     name;
  Uint1       qualclass;
} SourceQual, PNTR SourceQualPtr;

static SourceQual asn2gnbk_source_quals [ASN2GNBK_TOTAL_SOURCE] = {
  { "",                     Qual_class_ignore    },
  { "acronym",              Qual_class_orgmod    },
  { "anamorph",             Qual_class_orgmod    },
  { "authority",            Qual_class_orgmod    },
  { "biotype",              Qual_class_orgmod    },
  { "biovar",               Qual_class_orgmod    },
  { "breed",                Qual_class_orgmod    },
  { "cell_line",            Qual_class_subsource },
  { "cell_type",            Qual_class_subsource },
  { "chemovar",             Qual_class_orgmod    },
  { "chromosome",           Qual_class_subsource },
  { "citation",             Qual_class_pubset    },
  { "clone",                Qual_class_subsource },
  { "clone_lib",            Qual_class_subsource },
  { "common",               Qual_class_orgmod    },
  { "common",               Qual_class_string    },
  { "country",              Qual_class_subsource },
  { "cultivar",             Qual_class_orgmod    },
  { "db_xref",              Qual_class_db_xref   },
  { "db_xref",              Qual_class_db_xref   },
  { "dev_stage",            Qual_class_subsource },
  { "dosage",               Qual_class_orgmod    },
  { "ecotype",              Qual_class_orgmod    },
  { "endogenous_virus",     Qual_class_subsource },
  { "environmental_sample", Qual_class_subsource },
  { "extrachromosomal",     Qual_class_boolean   },
  { "focus",                Qual_class_boolean   },
  { "forma",                Qual_class_orgmod    },
  { "forma_specialis",      Qual_class_orgmod    },
  { "frequency",            Qual_class_subsource },
  { "gb_acronym",           Qual_class_orgmod    },
  { "gb_anamorph",          Qual_class_orgmod    },
  { "gb_synonym",           Qual_class_orgmod    },
  { "genotype",             Qual_class_subsource },
  { "germline",             Qual_class_subsource },
  { "group",                Qual_class_orgmod    },
  { "haplotype",            Qual_class_subsource },
  { "insertion_seq",        Qual_class_subsource },
  { "isolate",              Qual_class_orgmod    },
  { "isolation_source",     Qual_class_subsource },
  { "lab_host",             Qual_class_subsource },
  { "label",                Qual_class_label     },
  { "macronuclear",         Qual_class_boolean   },
  { "map",                  Qual_class_subsource },
  { "mol_type",             Qual_class_string    },
  { "note",                 Qual_class_note      },
  { "old_lineage",          Qual_class_orgmod    },
  { "old_name",             Qual_class_orgmod    },
  { "organism",             Qual_class_string    },
  { "organelle",            Qual_class_organelle },
  { "orgmod_note",          Qual_class_orgmod    },
  { "pathovar",             Qual_class_orgmod    },
  { "plasmid",              Qual_class_subsource },
  { "plastid",              Qual_class_subsource },
  { "pop_variant",          Qual_class_subsource },
  { "rearranged",           Qual_class_subsource },
  { "segment",              Qual_class_subsource },
  { "seqfeat_note",         Qual_class_string    },
  { "sequenced_mol",        Qual_class_quote     },
  { "serogroup",            Qual_class_orgmod    },
  { "serotype",             Qual_class_orgmod    },
  { "serovar",              Qual_class_orgmod    },
  { "sex",                  Qual_class_subsource },
  { "specific_host",        Qual_class_orgmod    },
  { "specimen_voucher",     Qual_class_orgmod    },
  { "strain",               Qual_class_orgmod    },
  { "sub_clone",            Qual_class_subsource },
  { "subgroup",             Qual_class_orgmod    },
  { "sub_species",          Qual_class_orgmod    },
  { "sub_strain",           Qual_class_orgmod    },
  { "subtype",              Qual_class_orgmod    },
  { "subsource_note",       Qual_class_subsource },
  { "synonym",              Qual_class_orgmod    },
  { "teleomorph",           Qual_class_orgmod    },
  { "tissue_lib",           Qual_class_subsource },
  { "tissue_type",          Qual_class_subsource },
  { "transgenic",           Qual_class_subsource },
  { "transposon",           Qual_class_subsource },
  { "type",                 Qual_class_orgmod    },
  { "unstructured",         Qual_class_valnode   },
  { "usedin",               Qual_class_quote     },
  { "variety",              Qual_class_orgmod    },
  { "?",                    Qual_class_orgmod    },
  { "?",                    Qual_class_orgmod    },
  { "?",                    Qual_class_subsource },
};

static Int2 subSourceToSourceIdx [30] = {
  SCQUAL_zero_subsrc,
  SCQUAL_chromosome,
  SCQUAL_map,
  SCQUAL_clone,
  SCQUAL_sub_clone,
  SCQUAL_haplotype,
  SCQUAL_genotype,
  SCQUAL_sex,
  SCQUAL_cell_line,
  SCQUAL_cell_type,
  SCQUAL_tissue_type,
  SCQUAL_clone_lib,
  SCQUAL_dev_stage,
  SCQUAL_frequency,
  SCQUAL_germline,
  SCQUAL_rearranged,
  SCQUAL_lab_host,
  SCQUAL_pop_variant,
  SCQUAL_tissue_lib,
  SCQUAL_plasmid_name,
  SCQUAL_transposon_name,
  SCQUAL_ins_seq_name,
  SCQUAL_plastid_name,
  SCQUAL_country,
  SCQUAL_segment,
  SCQUAL_endogenous_virus_name,
  SCQUAL_transgenic,
  SCQUAL_environmental_sample,
  SCQUAL_isolation_source,
  SCQUAL_subsource_note
};

typedef enum {
  FTQUAL_allele = 1,
  FTQUAL_anticodon,
  FTQUAL_bond,
  FTQUAL_bond_type,
  FTQUAL_bound_moiety,
  FTQUAL_cds_product,
  FTQUAL_citation,
  FTQUAL_clone,
  FTQUAL_coded_by,
  FTQUAL_codon,
  FTQUAL_codon_start,
  FTQUAL_cons_splice,
  FTQUAL_db_xref,
  FTQUAL_direction,
  FTQUAL_EC_number,
  FTQUAL_evidence,
  FTQUAL_exception,
  FTQUAL_exception_note,
  FTQUAL_figure,
  FTQUAL_frequency,
  FTQUAL_function,
  FTQUAL_gene,
  FTQUAL_gene_desc,
  FTQUAL_gene_allele,
  FTQUAL_gene_map,
  FTQUAL_gene_syn,
  FTQUAL_gene_note,
  FTQUAL_gene_xref,
  FTQUAL_heterogen,
  FTQUAL_illegal_qual,
  FTQUAL_insertion_seq,
  FTQUAL_label,
  FTQUAL_locus_tag,
  FTQUAL_map,
  FTQUAL_maploc,
  FTQUAL_mod_base,
  FTQUAL_modelev,
  FTQUAL_note,
  FTQUAL_number,
  FTQUAL_organism,
  FTQUAL_partial,
  FTQUAL_PCR_conditions,
  FTQUAL_phenotype,
  FTQUAL_product,
  FTQUAL_product_quals,
  FTQUAL_prot_activity,
  FTQUAL_prot_comment,
  FTQUAL_prot_EC_number,
  FTQUAL_prot_note,
  FTQUAL_prot_method,
  FTQUAL_prot_conflict,
  FTQUAL_prot_desc,
  FTQUAL_prot_missing,
  FTQUAL_prot_name,
  FTQUAL_prot_names,
  FTQUAL_protein_id,
  FTQUAL_pseudo,
  FTQUAL_region,
  FTQUAL_region_name,
  FTQUAL_replace,
  FTQUAL_rpt_family,
  FTQUAL_rpt_type,
  FTQUAL_rpt_unit,
  FTQUAL_rrna_its,
  FTQUAL_sec_str_type,
  FTQUAL_selenocysteine,
  FTQUAL_seqfeat_note,
  FTQUAL_site,
  FTQUAL_site_type,
  FTQUAL_standard_name,
  FTQUAL_transcript_id,
  FTQUAL_transl_except,
  FTQUAL_transl_table,
  FTQUAL_translation,
  FTQUAL_transposon,
  FTQUAL_trna_aa,
  FTQUAL_trna_codons,
  FTQUAL_usedin,
  FTQUAL_xtra_prod_quals,
  ASN2GNBK_TOTAL_FEATUR
}  FtQualType;

/* ordering arrays for qualifiers and note components */

static Uint1 feat_qual_order [] = {
  FTQUAL_partial,
  FTQUAL_gene,

  FTQUAL_product,

  FTQUAL_prot_EC_number,
  FTQUAL_prot_activity,

  FTQUAL_standard_name,
  FTQUAL_coded_by,

  FTQUAL_prot_name,
  FTQUAL_region_name,
  FTQUAL_bond_type,
  FTQUAL_site_type,
  FTQUAL_sec_str_type,
  FTQUAL_heterogen,

  FTQUAL_note, 
  FTQUAL_citation,

  FTQUAL_number,

  FTQUAL_pseudo,

  FTQUAL_codon_start,

  FTQUAL_anticodon,
  FTQUAL_bound_moiety,
  FTQUAL_clone,
  FTQUAL_cons_splice,
  FTQUAL_direction,
  FTQUAL_function,
  FTQUAL_evidence,
  FTQUAL_exception,
  FTQUAL_frequency,
  FTQUAL_EC_number,
  FTQUAL_gene_map,
  FTQUAL_gene_allele,
  FTQUAL_allele,
  FTQUAL_map,
  FTQUAL_mod_base,
  FTQUAL_PCR_conditions,
  FTQUAL_phenotype,
  FTQUAL_rpt_family,
  FTQUAL_rpt_type,
  FTQUAL_rpt_unit,
  FTQUAL_insertion_seq,
  FTQUAL_transposon,
  FTQUAL_usedin,

  FTQUAL_illegal_qual,

  FTQUAL_replace,

  FTQUAL_transl_except,
  FTQUAL_transl_table,
  FTQUAL_codon,
  FTQUAL_organism,
  FTQUAL_label,
  FTQUAL_cds_product,
  FTQUAL_protein_id,
  FTQUAL_transcript_id,
  FTQUAL_db_xref, 
  FTQUAL_gene_xref,
  FTQUAL_translation,
  0
};

/*
prot_names after seqfeat_note - gi|4210642|emb|AJ005084.1|HBVAJ5084
prot_conflict after prot_desc - gi|61183|emb|V01135.1|PIVM02
figure after prot_desc - gi|400553|gb|S64006.1|
seqfeat_note after prot_desc - gi|431713|gb|L20354.1|STVPATPOLB
  but prot_desc after seqfeat_note - AF252556.1
prot_names after figure - gi|234022|gb|S56149.1|S56149
seqfeat_note after prot_conflict after figure - gi|234046|gb|S51392.1|S51392
prot_method after prot_comment (descriptor) after prot_note after prot_desc
region after seqfeat_note - gi|6554164|gb|AF043644.3|AF043644
prot_desc after prot_names - gi|6581069|gb|AF202541.1|AF202541 - cannot do !!!
gene_syn after gene_desc - gi|3386543|gb|AF079528.1|AF079528
pseudo after note - gi|6598562|gb|AC006419.3|AC006419
*/

static Uint1 feat_note_order [] = {
  FTQUAL_locus_tag,
  FTQUAL_gene_desc,
  FTQUAL_gene_syn,
  FTQUAL_trna_codons,
  FTQUAL_prot_desc,
  FTQUAL_prot_note,
  FTQUAL_prot_comment,
  FTQUAL_prot_method,
  FTQUAL_figure,
  FTQUAL_maploc,
  FTQUAL_prot_conflict,
  FTQUAL_prot_missing,
  FTQUAL_seqfeat_note,
  FTQUAL_exception_note,
  FTQUAL_region,
  /* FTQUAL_selenocysteine, */
  FTQUAL_prot_names,
  FTQUAL_bond,
  FTQUAL_site,
  FTQUAL_rrna_its,
  FTQUAL_xtra_prod_quals,
  FTQUAL_modelev,
  0
};

typedef struct featurqual {
  CharPtr     name;
  Uint1       qualclass;
} FeaturQual, PNTR FeaturQualPtr;

static FeaturQual asn2gnbk_featur_quals [ASN2GNBK_TOTAL_FEATUR] = {
  { "",               Qual_class_ignore       },
  { "allele",         Qual_class_quote        },
  { "anticodon",      Qual_class_anti_codon   },
  { "bond",           Qual_class_bond         },
  { "bond_type",      Qual_class_bond         },
  { "bound_moiety",   Qual_class_quote        },
  { "product",        Qual_class_string       },
  { "citation",       Qual_class_pubset       },
  { "clone",          Qual_class_quote        },
  { "coded_by",       Qual_class_seq_loc      },
  { "codon",          Qual_class_codon        },
  { "codon_start",    Qual_class_int          },
  { "cons_splice",    Qual_class_consplice    },
  { "db_xref",        Qual_class_db_xref      },
  { "direction",      Qual_class_L_R_B        },
  { "EC_number",      Qual_class_EC_quote     },
  { "evidence",       Qual_class_evidence     },
  { "exception",      Qual_class_string       },
  { "exception_note", Qual_class_string       },
  { "figure",         Qual_class_string       },
  { "frequency",      Qual_class_quote        },
  { "function",       Qual_class_quote        },
  { "gene",           Qual_class_sgml         },
  { "gene_desc",      Qual_class_string       },
  { "allele",         Qual_class_string       },
  { "map",            Qual_class_string       },
  { "gene_syn",       Qual_class_gene_syn     },
  { "gene_note",      Qual_class_string       },
  { "db_xref",        Qual_class_db_xref      },
  { "heterogen",      Qual_class_string       },
  { "illegal",        Qual_class_illegal      },
  { "insertion_seq",  Qual_class_quote        },
  { "label",          Qual_class_label        },
  { "locus_tag",      Qual_class_string       },
  { "map",            Qual_class_quote        },
  { "maploc",         Qual_class_string       },
  { "mod_base",       Qual_class_noquote      },
  { "model_evidence", Qual_class_model_ev     },
  { "note",           Qual_class_note         },
  { "number",         Qual_class_number       },
  { "organism",       Qual_class_quote        },
  { "partial",        Qual_class_boolean      },
  { "PCR_conditions", Qual_class_quote        },
  { "phenotype",      Qual_class_quote        },
  { "product",        Qual_class_product      },
  { "product",        Qual_class_quote        },
  { "function",       Qual_class_valnode      },
  { "prot_comment",   Qual_class_string       },
  { "EC_number",      Qual_class_EC_valnode   },
  { "prot_note",      Qual_class_string       },
  { "prot_method",    Qual_class_method       },
  { "prot_conflict",  Qual_class_string       },
  { "prot_desc",      Qual_class_string       },
  { "prot_missing",   Qual_class_string       },
  { "name",           Qual_class_tilde        },
  { "prot_names",     Qual_class_protnames    },
  { "protein_id",     Qual_class_seq_id       },
  { "pseudo",         Qual_class_boolean      },
  { "region",         Qual_class_region       },
  { "region_name",    Qual_class_string       },
  { "replace",        Qual_class_replace      },
  { "rpt_family",     Qual_class_quote        },
  { "rpt_type",       Qual_class_rpt          },
  { "rpt_unit",       Qual_class_rpt_unit     },
  { "rrna_its",       Qual_class_its          },
  { "sec_str_type",   Qual_class_sec_str      },
  { "selenocysteine", Qual_class_string       },
  { "seqfeat_note",   Qual_class_string       },
  { "site",           Qual_class_site         },
  { "site_type",      Qual_class_site         },
  { "standard_name",  Qual_class_quote        },
  { "transcript_id",  Qual_class_seq_id       },
  { "transl_except",  Qual_class_code_break   },
  { "transl_table",   Qual_class_int          },
  { "translation",    Qual_class_translation  },
  { "transposon",     Qual_class_quote        },
  { "trna_aa",        Qual_class_ignore       },
  { "trna_codons",    Qual_class_trna_codons  },
  { "usedin",         Qual_class_paren        },
  { "xtra_products",  Qual_class_xtraprds     }
};

static Boolean AllowedValQual (Uint2 featureKey, FtQualType qualKey);


/* ********************************************************************** */

/* URLs */

#define MAX_WWWBUF 328

static Char link_feat [MAX_WWWBUF];
#define DEF_LINK_FEAT  "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?"

static Char link_seq [MAX_WWWBUF];
#define DEF_LINK_SEQ  "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?"

static Char link_wgs [MAX_WWWBUF];
#define DEF_LINK_WGS  "http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?"

static Char link_omim [MAX_WWWBUF];
#define DEF_LINK_OMIM "http://www.ncbi.nlm.nih.gov/entrez/dispomim.cgi?id="

static Char ref_link [MAX_WWWBUF];
#define DEF_LINK_REF  "http://www.ncbi.nlm.nih.gov/LocusLink/refseq.html"

static Char nt_link [MAX_WWWBUF];
#define DEF_LINK_NT  "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?view=graph&val="

static Char doc_link [MAX_WWWBUF];
#define DEF_LINK_DOC  "http://www.ncbi.nlm.nih.gov/genome/guide/build.html"

static Char ev_link [MAX_WWWBUF];
#define DEF_LINK_EV  "http://www.ncbi.nlm.nih.gov/cgi-bin/Entrez/evv.cgi?"

static Char ec_link [MAX_WWWBUF];
#define DEF_LINK_EC "http://www.expasy.ch/cgi-bin/nicezyme.pl?"

static Char link_tax [MAX_WWWBUF];
#define DEF_LINK_TAX "/htbin-post/Taxonomy/wgetorg?"

static Char link_ff[MAX_WWWBUF];
#define DEF_LINK_FF  "/cgi-bin/Entrez/getfeat?"

static Char link_muid[MAX_WWWBUF];
#define DEF_LINK_MUID  "/entrez/utils/qmap.cgi?"

static Char link_ace[MAX_WWWBUF];
#define DEF_LINK_ACE  "http://www.ncbi.nlm.nih.gov/AceView/hs.cgi?"

static Char link_code[MAX_WWWBUF];
#define DEF_LINK_CODE "/htbin-post/Taxonomy/wprintgc?"

static Char link_fly[MAX_WWWBUF];
#define DEF_LINK_FLY "http://flybase.bio.indiana.edu/.bin/fbidq.html?"

static Char link_fly_fban[MAX_WWWBUF];
#define DEF_LINK_FBAN "http://www.fruitfly.org/cgi-bin/annot/fban?"

static Char link_fly_fbgn[MAX_WWWBUF];
#define DEF_LINK_FBGN "http://flybase.bio.indiana.edu/.bin/fbidq.html?"

static Char link_cog[MAX_WWWBUF];
#define DEF_LINK_COG "http://www.ncbi.nlm.nih.gov/cgi-bin/COG/palox?"

static Char link_sgd[MAX_WWWBUF];
#define DEF_LINK_SGD "/cgi-bin/Entrez/referer?http://genome-www4.stanford.edu/cgi-bin/SGD/locus.pl?locus="

static Char link_gdb[MAX_WWWBUF];
#define DEF_LINK_GDB "http://gdbwww.gdb.org/gdb-bin/genera/genera/hgd/DBObject/GDB:"

static Char link_ck[MAX_WWWBUF];
#define DEF_LINK_CK "http://flybane.berkeley.edu/cgi-bin/cDNA/CK_clone.pl?db=CK&dbid="

static Char link_rice[MAX_WWWBUF];
#define DEF_LINK_RICE "http://ars-genome.cornell.edu/cgi-bin/WebAce/webace?db=ricegenes&class=Marker&object="

static Char link_sp[MAX_WWWBUF];
#define DEF_LINK_SP "/cgi-bin/Entrez/referer?http://expasy.hcuge.ch/cgi-bin/sprot-search-ac%3f"

static Char link_pdb[MAX_WWWBUF];
#define DEF_LINK_PDB "/cgi-bin/Entrez/referer?http://expasy.hcuge.ch/cgi-bin/get-pdb-entry%3f"

static Char link_UniSTS[MAX_WWWBUF];
#define DEF_LINK_UniSTS "http://www.ncbi.nlm.nih.gov/genome/sts/sts.cgi?uid="

static Char link_dbSTS[MAX_WWWBUF];
#define DEF_LINK_dbSTS "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?"

static Char link_dbEST[MAX_WWWBUF];
#define DEF_LINK_dbEST "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?"

static Char link_omim[MAX_WWWBUF];
#define DEF_LINK_OMIM "http://www.ncbi.nlm.nih.gov/entrez/dispomim.cgi?id="

static Char link_locus[MAX_WWWBUF];
#define DEF_LINK_LOCUS "http://www.ncbi.nlm.nih.gov/LocusLink/LocRpt.cgi?l="

static Char link_snp[MAX_WWWBUF];
#define DEF_LINK_SNP "http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?type=rs&rs="

static Char link_ratmap[MAX_WWWBUF];
#define DEF_LINK_RATMAP "http://ratmap.gen.gu.se/action.lasso?-database=RATMAPfmPro&-layout=Detail&-response=/RM/Detail+Format.html&-search&-recid="

static Char link_rgd[MAX_WWWBUF];
#define DEF_LINK_RGD "http://rgd.mcw.edu/query/query.cgi?id="

static Char link_mgd[MAX_WWWBUF];
#define DEF_LINK_MGD "http://www.informatics.jax.org/searches/accession_report.cgi?id=MGI:"

static Char link_cdd[MAX_WWWBUF];
#define DEF_LINK_CDD "http://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid="

static Char link_niaest[MAX_WWWBUF];
#define DEF_LINK_NIAEST "http://lgsun.grc.nia.nih.gov/cgi-bin/pro3?sname1="

static Char link_worm_sequence[MAX_WWWBUF];
#define DEF_LINK_WORM_SEQUENCE "http://www.wormbase.org/db/seq/sequence?name="

static Char link_worm_locus[MAX_WWWBUF];
#define DEF_LINK_WORM_LOCUS "http://www.wormbase.org/db/gene/locus?name="

static Char link_imgt[MAX_WWWBUF];
#define DEF_LINK_IMGT "http://imgt.cines.fr:8104/cgi-bin/IMGTlect.jv?query=202+"

static Char link_ifo[MAX_WWWBUF];
#define DEF_LINK_IFO "http://www.ifo.or.jp/index_e.html"

static Char link_jcm[MAX_WWWBUF];
#define DEF_LINK_JCM "http://www.jcm.riken.go.jp/cgi-bin/jcm/jcm_number?JCM="

static Char link_isfinder[MAX_WWWBUF];
#define DEF_LINK_ISFINDER "http://www-is.biotoul.fr/scripts/is/is_spec.idc?name="

static Char link_gabi[MAX_WWWBUF];
#define DEF_LINK_GABI "https://gabi.rzpd.de/cgi-bin-protected/GreenCards.pl.cgi?Mode=ShowBioObject&BioObjectName="

static Char link_fantom[MAX_WWWBUF];
#define DEF_LINK_FANTOM "http://fantom.gsc.riken.go.jp/db/view/main.cgi?masterid="


/* utility functions */

static StringItemPtr FFGetString (IntAsn2gbJobPtr ajp)

{
  StringItemPtr  sip;

  if (ajp == NULL) return NULL;
  if (ajp->pool != NULL) {
    sip = ajp->pool;
    ajp->pool = sip->next;
    sip->next = NULL;
    MemSet ((Pointer) sip, 0, sizeof (StringItem));
  } else {
    sip = (StringItemPtr) MemNew (sizeof (StringItem));
    if (sip == NULL) return NULL;
  }
  sip->curr = sip;
  sip->iajp = ajp;
  sip->pos = 0;
  return sip;
}

static void FFRecycleString (IntAsn2gbJobPtr ajp, StringItemPtr ffstring)

{
  StringItemPtr  nxt;

  if (ajp == NULL || ffstring == NULL) return;
  if ( ffstring->pos == -1 ) return;
  
  nxt = ffstring;
  nxt->pos = -1;
  while (nxt->next != NULL) {
    nxt->pos = -1;
    nxt = nxt->next;
  }
  nxt->next = ajp->pool;
  ajp->pool = ffstring;

  ffstring->curr = NULL;
}

static void FFAddOneChar (
  StringItemPtr sip, 
  Char ch,
  Boolean convertQuotes
)
{
  StringItemPtr current = sip->curr;

  if ( current->pos == STRING_BUF_LEN ) {
    current->next = FFGetString(sip->iajp);
    current = current->next;
    current->pos = 0;
    sip->curr = current;
  }

  if ( convertQuotes && ch == '\"' ) {
    ch = '\'';
  }
  current->buf[current->pos] = ch;
  current->pos++;
}

static void FFAddNewLine(StringItemPtr ffstring) {
  FFAddOneChar(ffstring, '\n', FALSE);
}

static void FFAddNChar (
  StringItemPtr sip, 
  Char ch,
  Int4 n,
  Boolean convertQuotes
)
{
  Int4 i;

  for ( i = 0; i < n; ++i ) {
    FFAddOneChar(sip, ch, convertQuotes);
  }
}
  

static void FFExpandTildes (StringItemPtr sip, CharPtr PNTR cpp) {
  Char replace = **cpp;

  if ( **cpp == '~' ) {
    if ( *((*cpp) + 1) == '~' ) {     /* "~~" -> '~' */
      replace = '~';
      (*cpp)++;
    } else {
      replace = '\n';
    }
  } 

  FFAddOneChar(sip, replace, FALSE);
}


static void FFReplaceTildesWithSpaces (StringItemPtr ffstring, CharPtr PNTR cpp) {
  Char replace = **cpp, lookahead;
  CharPtr cptr = *cpp;
  
  if ( *cptr == '`' ) {
    FFAddOneChar(ffstring, replace, FALSE);
    return;
  }

  replace = ' ';
  lookahead = *(cptr + 1);

  if ( IS_DIGIT(lookahead) ) {
    replace = '~';
  }
  else {
    if ( (lookahead == ' ') || (lookahead == '(') ) {
      if ( IS_DIGIT(*(cptr + 2)) ) {
        replace = '~';
      }
    }
  }

  FFAddOneChar(ffstring, replace, FALSE);
}

static void FFOldExpand (StringItemPtr sip, CharPtr PNTR cpp) {
  /* "~" -> "\n", "~~" or "~~ ~~" -> "\n\n" */ 
  CharPtr cp = *cpp;
  Char current = *cp;
  Char next = *(cp + 1);
  
  /* handle "'~" */
  if ( current == '`' ) {
    if ( next != '~' ) {
        FFAddOneChar(sip, current, FALSE);
    } else {
        FFAddOneChar(sip, '~', FALSE);
        (*cpp)++;
    }
    return;
  }

  /* handle "~", "~~" or "~~ ~~" */
  FFAddOneChar(sip, '\n', FALSE);
  if ( next == '~' ) {
    FFAddOneChar(sip, '\n', FALSE);
    cp++;
    *cpp = cp;
    cp++;
    if ( *cp == ' ' ) {
      cp++;
      if ( *cp == '~' ) {
        cp++;
        if ( *cp == '~' ) { /* saw "~~ ~~" */
          *cpp = cp;
        }
      }
    }
  }
}

static void AddStringWithTildes (StringItemPtr ffstring, CharPtr string)
{
/* One "~" is a  new line, "~~" or "~~ ~~" means 2 returns */       

    while (*string != '\0') {
        if (*string == '`' && *(string+1) == '~') {
            FFAddOneChar(ffstring, '~', FALSE);
            string += 2;
        } else if (*string == '~') {
            FFAddOneChar(ffstring, '\n', FALSE);
            string++;
            if (*string == '~') {
                FFAddOneChar(ffstring, '\n', FALSE);
                string++;
        if (*string == ' ' && *(string+1) == '~' && *(string+2) == '~') {
                    string += 3;
        }
      }
        } else if (*string == '\"') {   
            *string = '\'';
            FFAddOneChar(ffstring, *string, FALSE);
            string++;
        } else {  
            FFAddOneChar(ffstring, *string, FALSE);
            string++;
        }
    }
}    /* AddStringWithTildes */


static void FFProcessTildes (StringItemPtr sip, CharPtr PNTR cpp, Int2 tildeAction) {
    
  switch (tildeAction) {

  case TILDE_EXPAND :
      FFExpandTildes(sip, cpp);
      break;

  case TILDE_OLD_EXPAND :
      FFOldExpand(sip, cpp);
      break;

  case TILDE_TO_SPACES :
      FFReplaceTildesWithSpaces (sip, cpp);
      break;

  case TILDE_IGNORE:
  default:
      FFAddOneChar(sip, **cpp, FALSE);
      break;
  }
}

static void FFAddPeriod (StringItemPtr sip) {
  Int4 i;
  Char ch  = '\0';
  StringItemPtr riter = sip->curr, prev;
  IntAsn2gbJobPtr ajp;

  if ( sip == NULL ) return;
  ajp = (IntAsn2gbJobPtr)sip->iajp;
  if ( ajp == NULL ) return;

  for ( i = riter->pos - 1; i >= 0; --i ) {
    ch = riter->buf[i];

    if ( (ch == ' ') || (ch == '\t')  || (ch == '~')  || (ch == '.') || (ch == '\n')) {
      riter->pos--;
      
      if ( i < 0 && riter != sip ) {
        for ( prev = sip; prev->next != NULL; prev = prev->next ) {
          if ( prev->next == riter ) {
            i = prev->pos - 1;
            FFRecycleString(ajp, riter);
            riter = prev;
            riter->next = NULL;
	          sip->curr = riter;
            break;
          }
        }
      }

    } else {
      break;
    }
  }

  if (ch != '.') {
    FFAddOneChar(sip, '.', FALSE);
  }
}

static void FFAddOneString (
  StringItemPtr sip, 
  CharPtr string,
  Boolean addPeriod, 
  Boolean convertQuotes,
  Int2 tildeAction
)
{
  CharPtr strp = string;

  if ( string == NULL ) return;
  
  while ( *strp != '\0' ) {
    if ( (*strp == '`') || (*strp == '~') ) {
      FFProcessTildes(sip, &strp, tildeAction);
    } else {
      FFAddOneChar(sip, *strp, convertQuotes);
    }
    strp++;
  }

  if ( addPeriod ) {
    FFAddPeriod(sip);
  }
}

static Boolean GetWWW (IntAsn2gbJobPtr ajp);
static Boolean FFIsStartOfLink (StringItemPtr iter, Int4 pos);

static void FFCatenateSubString (
  StringItemPtr dest,
  StringItemPtr start_sip, Int4 start_pos,
  StringItemPtr end_sip, Int4 end_pos
)
{
  Int4 max_i, min_i, i;
  StringItemPtr current;
  Boolean in_url = FALSE;
  IntAsn2gbJobPtr ajp = (IntAsn2gbJobPtr)dest->iajp;

  if ( GetWWW(ajp) ) {
    for ( current = start_sip, i = start_pos;
    current != NULL; 
    current = current->next ) {
      if ( current == start_sip ) {
        min_i = start_pos;
      } else {
        min_i = 0;
      }
      
      if ( current == end_sip ) {
        max_i = end_pos;
      } else {
        max_i = current->pos;
      }
      
      for ( i = min_i; i < max_i; ++i ) {
        if ( current->buf[i] == '<' ) {
          if ( !FFIsStartOfLink(current, i) ) {
            FFAddOneString(dest, "&lt;", FALSE, FALSE, TILDE_IGNORE);
            continue;
          } else {
            in_url = TRUE;
          }
        }
        if ( current->buf[i] == '>' ) {
          if ( !in_url ) {
            FFAddOneString(dest, "&gt;", FALSE, FALSE, TILDE_IGNORE);
            continue;
          } else {
            in_url = FALSE;
          }
        } 

        FFAddOneChar(dest, current->buf[i], FALSE);
      }

      if ( current == end_sip ) break;
    }
  } else {
    for ( current = start_sip, i = start_pos;
    current != NULL; 
    current = current->next ) {
      if ( current == start_sip ) {
        min_i = start_pos;
      } else {
        min_i = 0;
      }
      
      if ( current == end_sip ) {
        max_i = end_pos;
      } else {
        max_i = current->pos;
      }
      
      for ( i = min_i; i < max_i; ++i ) {
        FFAddOneChar(dest, current->buf[i], FALSE);
      }
      
      if ( current == end_sip ) break;
    }
  }
}


static CharPtr FFToCharPtr (StringItemPtr sip) {
  Int4 size = 0, i;
  StringItemPtr iter;
  CharPtr result, temp;

  for ( iter = sip; iter != NULL; iter = iter->next ) {
    size += iter->pos;
  }

  result = (CharPtr)MemNew(size + 2);
  temp = result;

  for ( iter = sip; iter != NULL; iter = iter->next ) {
    for ( i = 0; i < iter->pos; ++i ) {
      *temp = iter->buf[i];
      ++temp;
    }
  }

  *temp = '\0';

  return result;
}


/* word wrap functions */

static void FFSkipLink (StringItemPtr PNTR iterp, Int4Ptr ip) {
  StringItemPtr iter = *iterp;
  Int4 i = *ip;

  while ( (iter != NULL) && (iter->buf[i] != '>') ) {
    ++i;

    if ( i == iter->pos ) {
      iter = iter->next;
      i = 0;
    }
  }
  ++i;
  if ( i == iter->pos && iter->next != NULL ) {
    iter = iter->next;
    i = 0;
  }

  *iterp = iter;
  *ip = i;
}

static Boolean FFIsStartOfLink (StringItemPtr iter, Int4 pos)  {
  static CharPtr start_link = "<A HREF";
  static CharPtr end_link = "</A>";
  Int4 start_len = StringLen(start_link);
  Int4 end_len = StringLen(end_link);
  Char temp[10];
  Int4 i;
  StringItemPtr current = NULL;

  if ( iter == NULL || pos >= iter->pos ) return FALSE;
  if ( iter->buf[pos] != '<' ) return FALSE;

  MemSet(temp, 0, sizeof(temp));
  for ( i = 0; i < start_len && iter != NULL; ++i ) {
    if ( pos + i < iter->pos ) {
      temp[i] = iter->buf[pos+i];
      if ( i == end_len - 1 ) {
        if ( StringNICmp(temp, end_link, end_len) == 0 ) {
          return TRUE;
        }
      }
    } else {
      iter = iter->next;
      pos = -i;
      --i;
    }
  }

  if ( i == start_len ) {
    if ( StringNICmp(temp, start_link, start_len) == 0 ) {
        return TRUE;
    }
  }

  return FALSE;
}


static void FFSavePosition(StringItemPtr ffstring, StringItemPtr PNTR bufptr, Int4 PNTR posptr) {
  *bufptr = ffstring->curr;
  *posptr = ffstring->curr->pos;
}


static void FFTrim (
    StringItemPtr ffstring,
    StringItemPtr line_start,
    Int4 line_pos,
    Int4 line_prefix_len
)
{
  StringItemPtr riter, iter;
  Int4 i;
  IntAsn2gbJobPtr ajp = (IntAsn2gbJobPtr)ffstring->iajp;

  for ( i = 0; i < line_prefix_len; ++i ) {
    ++line_pos;
    if ( line_pos == STRING_BUF_LEN ) {
      line_pos = 0;
      line_start= line_start->next;
    }
  }

  riter = ffstring->curr;
  while ( riter != NULL ) {
    for ( i = riter->pos - 1;
          (i >= 0) && !(riter == line_start && i <= line_pos);
          --i ) {
      if ( !IS_WHITESP(riter->buf[i]) || (riter->buf[i] == '\n') ) {
        break;
      }
    }
    if ( i < 0 ) {
      i = STRING_BUF_LEN - 1;
      for ( iter = ffstring; iter != NULL; iter = iter->next ) {
        if ( iter->next == riter ) {
          break;
        }
      }
      if ( iter == NULL ){
        ffstring->pos = 0;
        break;
      } else {
        
        riter = iter;
        ffstring->curr = riter;
      }
    } else {
      riter->pos = i + 1;
      FFRecycleString(ajp, riter->next);
      riter->next = NULL;
      break;
    }
  }
}



/* A line is wrapped when the visble text in th eline exceeds the line size. */
/* Visible text is text that is not an HTML hyper-link.                      */
/* A line may be broken in one of the following characters:                  */
/* space, comma and dash                                                     */
/* the oredr of search is first spaces, then commas and then dashes.         */
/* We nee to take into account the possiblity that a 'new-line' character    */
/* already exists in the line, in such case we break at the 'new-line'       */
/* spaces, dashes and new-lines will be broken at that character wheras for  */
/* commas we break at the character following the comma.                     */
static void FFCalculateLineBreak (
  StringItemPtr PNTR break_sip, Int4 PNTR break_pos,
  Int4 init_indent, Int4 visible
)
{
  StringItemPtr iter, prev;
  Int4 i,
       done = FALSE,
       copied = 0, 
       start = *break_pos;
  Char ch;
  Boolean found_comma = FALSE, found_dash = FALSE;
  /* each candidate is a pair of buffer and position withingh this buffer */
  StringItemPtr candidate_sip_space = NULL,
                candidate_sip_comma = NULL,
                candidate_sip_dash  = NULL;
  Int4          candidate_int_space = -1,  
                candidate_int_comma = -1,
                candidate_int_dash  = -1;
  

  iter = *break_sip;
  prev = iter;

  /* skip the first 'init_indent' characters of the line */
  while ( iter != NULL && !done ) {
    for ( i = start; i < iter->pos && init_indent > 0; ++i ) {
      if ( iter->buf[i] == '\n' ) {
        candidate_sip_space = iter;
        candidate_int_space = i;
        done = TRUE;
        break;
      }
      if ( FFIsStartOfLink(iter, i) ) {
        FFSkipLink(&iter, &i);
        --i;
        continue;
      }

      --init_indent;
      ++copied;
    }
    if ( init_indent > 0 ) {
      start = 0;
      iter = iter->next;
    } else {
      break;
    }
  }
  start = i;

  while ( iter != NULL && !done ) {
    for ( i = start; i < iter->pos; ++i ) {
      if ( found_comma ) {
        candidate_sip_comma = iter;
        candidate_int_comma = i;
        found_comma = FALSE;
      }
      if ( found_dash ) {
        candidate_sip_dash = iter;
        candidate_int_dash = i;
        found_dash= FALSE;
      }

      ch = iter->buf[i];
      if ( ch == '\n' ) {
        candidate_sip_space = iter;
        candidate_int_space = i;
        done = TRUE;
        break;
      } else if ( ch == ' ' ) {
        candidate_sip_space = iter;
        candidate_int_space = i;
      } else if ( ch == ',' ) {
        found_comma = TRUE;
      } else if ( ch == '-' ) {
        found_dash = TRUE;
        /*candidate_sip_dash = iter;
        candidate_int_dash = i;*/
      }

      if ( FFIsStartOfLink(iter, i) ) {
        FFSkipLink(&iter, &i);
        --i;
        continue;
      }

      ++copied;
      if ( copied >= visible ) {
        if ( (candidate_sip_space == NULL) && (candidate_int_space == -1) &&
             (candidate_sip_comma == NULL) && (candidate_int_comma == -1) &&
             (candidate_sip_dash == NULL)  && (candidate_int_dash == -1)  ) {
          candidate_sip_space = iter;
          candidate_int_space = i;
        }
        done = TRUE;
        break;
      }      
    }
    start = 0;
    if ( !done ) {
      prev = iter;
      iter = iter->next;
    }
  }
  
  /* the order in which we examine the various candidate breaks is important */
  if ( iter == NULL && !done) { /* reached the end */
    *break_sip = prev;
    *break_pos = prev->pos;
  } else {
    if( candidate_sip_space != NULL ) {
        *break_sip = candidate_sip_space;
        *break_pos = candidate_int_space;
    } else if( candidate_sip_comma != NULL ) {
        *break_sip = candidate_sip_comma;
      *break_pos = candidate_int_comma;
    } else if( candidate_sip_dash != NULL ) {
      *break_sip = candidate_sip_dash;
      *break_pos = candidate_int_dash;
    }
  }
}

static void FFLineWrap (
  StringItemPtr dest, 
  StringItemPtr src, 
  Int4 init_indent,
  Int4 cont_indent, 
  Int4 line_max,
  CharPtr eb_line_prefix
)
{
  /* line break candidate is a pair <StringItepPtr, position> */
  StringItemPtr break_sip = src;
  Int4          break_pos = 0;
  StringItemPtr line_start = NULL;
  Int4          line_pos = 0;
  Int4          i, line_prefix_len = 0;
  StringItemPtr iter;

  FFSavePosition(dest, &line_start, &line_pos);

  for ( iter = src; iter != NULL; iter = iter->next ) {
    for ( i = 0; i < iter->pos; ) {
      break_pos = i;
      break_sip = iter;

      FFCalculateLineBreak(&break_sip, &break_pos, init_indent, line_max - line_prefix_len + 1);
      FFCatenateSubString(dest, iter, i, break_sip, break_pos);
      FFTrim(dest, line_start, line_pos, cont_indent);
      FFAddOneChar(dest, '\n', FALSE);
      
      FFSavePosition(dest, &line_start, &line_pos);

      i = break_pos;
      iter = break_sip;

      if ( iter->buf[i-1] == 'X' && iter->buf[i-2] == 'X') {
        if ( (i == 2) || ((i > 2) && (iter->buf[i-3] == '\n')) ) {
          ++i;
          continue;
        } 
      }

      if ( IS_WHITESP(iter->buf[i]) ) {
        i++;
      }
      if ( iter != src->curr || i < iter->pos ) {
        if ( eb_line_prefix != NULL ) {
          FFAddOneString(dest, eb_line_prefix, FALSE, FALSE, TILDE_IGNORE);
        }
        FFAddNChar(dest, ' ', cont_indent - StringLen(eb_line_prefix), FALSE);
        init_indent = 0;
        line_prefix_len = cont_indent;
        /*FFSkipGarbage(&iter, &i);*/
      }
    }
  }
}

/* === */

static void FFStartPrint (
  StringItemPtr sip,
  FmtType format,
  Int4 gb_init_indent,
  Int4 gb_cont_indent,
  CharPtr gb_label,
  Int4 gb_tab_to,
  Int4 eb_init_indent,
  Int4 eb_cont_indent,
  CharPtr eb_line_prefix,
  Boolean eb_print_xx 
)

{
  if (format == GENBANK_FMT || format == GENPEPT_FMT) {
    FFAddNChar(sip, ' ', gb_init_indent, FALSE);
    FFAddOneString(sip, gb_label, FALSE, FALSE, TILDE_IGNORE);
    FFAddNChar(sip, ' ', gb_tab_to - gb_init_indent - StringLen(gb_label), FALSE);
  } else if (format == EMBL_FMT || format == EMBLPEPT_FMT) {
    if ( eb_print_xx ) {
      FFAddOneString(sip, "XX\n", FALSE, FALSE, TILDE_IGNORE);
    }
    FFAddOneString(sip, eb_line_prefix, FALSE, FALSE, TILDE_IGNORE);
    FFAddNChar(sip, ' ', eb_init_indent - StringLen(eb_line_prefix), FALSE);
  }
}

static void FFAddTextToString (
  StringItemPtr ffstring, 
  CharPtr prefix,
  CharPtr string,
  CharPtr suffix,
  Boolean addPeriod,
  Boolean convertQuotes,
  Int2 tildeAction
)

{
  FFAddOneString (ffstring, prefix, FALSE, FALSE, TILDE_IGNORE);
  FFAddOneString (ffstring, string, FALSE, convertQuotes, tildeAction);
  FFAddOneString (ffstring, suffix, FALSE, FALSE, TILDE_IGNORE);

  if ( addPeriod ) {
    FFAddPeriod(ffstring);
  }
}
   

static CharPtr FFEndPrint (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  FmtType format,
  Int2 gb_init_indent,
  Int2 gb_cont_indent,
  Int2 eb_init_indent,
  Int2 eb_cont_indent,
  CharPtr eb_line_prefix
)
{
  StringItemPtr temp = FFGetString(ajp);
  CharPtr result;

  if ( (ffstring == NULL) || (ajp == NULL) ) return NULL;

  if (format == GENBANK_FMT || format == GENPEPT_FMT) {
    FFLineWrap(temp, ffstring, gb_init_indent, gb_cont_indent, ASN2FF_GB_MAX, NULL);
  } else {
    FFLineWrap(temp, ffstring, eb_init_indent, eb_cont_indent, ASN2FF_EMBL_MAX, eb_line_prefix);
  }

  result = FFToCharPtr(temp);
  FFRecycleString(ajp, temp);
  return result;
}

static Uint4 FFLength(StringItemPtr ffstring) {
  Uint4 len = 0;
  StringItemPtr current;

  for ( current = ffstring; current != NULL; current = current->next ) {
    len += current->pos;
  }

  return len;
}


static Char FFCharAt(StringItemPtr ffstring, Uint4 pos) {
  Uint4 count = 0, inbufpos;
  StringItemPtr current = NULL;

  inbufpos = pos % STRING_BUF_LEN;
  
  for ( current = ffstring; current != NULL; current = current->next ) {
    count += current->pos;
    if ( count > pos ) break;
  }

  if ( current != NULL && inbufpos <= pos )  {
    return current->buf[inbufpos];
  }

  return '\0';
}


static Char FFFindChar (
  StringItemPtr ffstring,   /* StringItem to search in */
  StringItemPtr start_buf,  /* the position of the last char searched for (buffer) */
  Uint4 start_pos,          /* the position of the last char searched for (pos) */
  Uint4 old_pos,         /* the global position searched for */
  Uint4 new_pos             /* new search position */
)
{
  Uint4 delta;
  Uint4 count;
  StringItemPtr current = NULL;

  Char result = '\0';

  if ( new_pos == old_pos ) {
    result = start_buf->buf[start_pos];
  } 

  if ( new_pos > old_pos ) {
    delta = new_pos - old_pos;
    current = start_buf;
    count = current->pos - start_pos - 1;
    current = current->next;
    
    while ( delta > count && current != NULL ) {
      current = current->next;
      count += current->pos;
    }
    
    if ( current != NULL  )  {
      result = current->buf[new_pos % STRING_BUF_LEN];
    }
    
  } else /* new_pos < old_pos */ {
    delta = old_pos - new_pos;
    if ( old_pos % STRING_BUF_LEN >= delta ) {
      result = start_buf->buf[new_pos % STRING_BUF_LEN];
    } else {
      result = FFCharAt(ffstring, new_pos);
    }
  }

  return result;
}

static Boolean FFEmpty(StringItemPtr ffstring) {
  if ( ffstring != NULL && ffstring->pos != 0 ) {
    return FALSE;
  }
  return TRUE;
}

/*
 * Compute the right-most position in the pattern at which character a occurs,
 * for each character a in the alphabet (assumed ASCII-ISO 8859-1)
 * 
 * The result is returned in the supplied vector.
 */
static void ComputeLastOccurance(const CharPtr pattern, Uint4 last_occurance[])
{
    Uint4 i;
    Uint4 pat_len;

    /* Initilalize vector */
    for ( i = 0; i < 256; ++i ) {
        last_occurance[i] = 0;
    }

    /* compute right-most occurance */
    pat_len = StringLen(pattern);
    for ( i = 0; i < pat_len; ++i ) {
        last_occurance[(Uint1)pattern[i]] = i;
    }
}

static void ComputePrefix(const CharPtr pattern, Uint4 longest_prefix[])
{
    Uint4 pat_len = StringLen(pattern);
    Uint4 k, q;

    longest_prefix[0] = 0;

    k = 0;
    for ( q = 1; q < pat_len; ++q ) {
        while ( k > 0 && pattern[k] != pattern[q] ) {
            k = longest_prefix[k - 1];
        }
        if ( pattern[k] == pattern[q] ) {
            ++k;
        }
        longest_prefix[q] = k;
    }
}


static void ComputeGoodSuffix(const CharPtr pattern, Uint4 good_suffix[])
{
    Uint4 pat_len = StringLen(pattern);
    Uint4Ptr longest_prefix, reverse_longest_prefix;
    CharPtr reverse_pattern;
    Uint4 i, j;

    /* allocate memory */
    longest_prefix = MemNew(pat_len * sizeof(Uint4));
    reverse_longest_prefix = MemNew(pat_len * sizeof(Uint4));
    reverse_pattern = MemNew((pat_len + 1) * sizeof(Char));

    if ( longest_prefix == NULL  ||
         reverse_longest_prefix == NULL  ||
         reverse_pattern == NULL ) {
      MemFree(longest_prefix);
      MemFree(reverse_longest_prefix);
      MemFree(reverse_pattern);
      return;
    }

    /* compute reverse pattern */
    for ( i = 0; i < pat_len; ++i ) {
      reverse_pattern[pat_len - i] = pattern[i];
    }

    ComputePrefix(pattern, longest_prefix);
    ComputePrefix(reverse_pattern, reverse_longest_prefix);

    for ( j = 0; j < pat_len; ++j) {
        good_suffix[j] = pat_len - longest_prefix[pat_len-1];
    }

    for ( i = 0; i < pat_len; ++i ) {
        j = pat_len - reverse_longest_prefix[i] - 1;
        if ( good_suffix[j] > i - reverse_longest_prefix[i] + 1) {
            good_suffix[j] = i - reverse_longest_prefix[i] + 1;
        }
    }
}


/*
 * searches for a pattern in a StringItem.
 * Using the Boyer-Moore algorithm for the search.
 */
static Int4 FFStringSearch (
  StringItemPtr text,
  const CharPtr pattern,
  Uint4 position )
{
  Uint4 text_len = FFLength(text);
  Uint4 pat_len = StringLen(pattern);
  Uint4 last_occurance[256];
  Uint4Ptr good_suffix;
  Uint4 shift;
  Int4 j;

  if ( pat_len == 0 ) return 0;
  if ( text_len == 0 || pat_len > text_len - position ) return -1;
  
  good_suffix = (Uint4Ptr)MemNew(pat_len * sizeof(Int4));
  if ( good_suffix == NULL ) return -1;

  ComputeLastOccurance(pattern, last_occurance);
  ComputeGoodSuffix(pattern, good_suffix);

  shift = position;
  while ( shift <= text_len - pat_len ) {
    j = pat_len - 1;
    while( j >= 0 && pattern[j] == FFCharAt(text,shift + j) ) {
      --j;
    }
    if ( j == -1 ) {
      return shift;
    } else {
        shift += MAX( (Int4)good_suffix[j],
		      (Int4)(j - last_occurance[FFCharAt(text,shift + j)]));
    }
  }

  return -1;
}


/*                                                                   */
/* IsWholeWordSubstr () -- Determines if a substring that is         */
/*                         contained in another string is a whole    */
/*                         word or phrase -- i.e. is it both         */
/*                         preceded and followed by white space.     */
/*                                                                   */
static Boolean IsWholeWordSubstr (
  StringItemPtr searchStr,
  Uint4 foundPos,
  CharPtr subStr
)
{
	Boolean left, right;
	Char ch;


	/* check on the left only if there is a character there */
	if (foundPos > 0) {
		ch = FFCharAt(searchStr, foundPos - 1);
		left = IS_WHITESP(ch) || ispunct(ch);
	} else {
		left = TRUE;
	}

	foundPos += StringLen(subStr);
  if ( foundPos == FFLength(searchStr) ) {
    right = TRUE;
  } else {
    ch = FFCharAt(searchStr, foundPos);
	  right = IS_WHITESP(ch) || ispunct(ch);
  }

	return left; /* see comment above */
  /* return left && right;  this is how it should be!*/
}


/* www utility functions */

static Boolean GetWWW (IntAsn2gbJobPtr ajp) {
    return ajp->www;
}

static void FiniWWW (IntAsn2gbJobPtr ajp) {
    ajp->www = FALSE;
}

static void InitWWW (IntAsn2gbJobPtr ajp)
{
  ajp->www = TRUE;

  GetAppParam ("NCBI", "WWWENTREZ", "LINK_FEAT", DEF_LINK_FEAT, link_feat, MAX_WWWBUF);
  GetAppParam ("NCBI", "WWWENTREZ", "LINK_SEQ", DEF_LINK_SEQ, link_seq, MAX_WWWBUF);
  GetAppParam ("NCBI", "WWWENTREZ", "LINK_WGS", DEF_LINK_WGS, link_wgs, MAX_WWWBUF);
  GetAppParam ("NCBI", "WWWENTREZ", "LINK_OMIM", DEF_LINK_OMIM, link_omim, MAX_WWWBUF);
  GetAppParam ("NCBI", "WWWENTREZ", "LINK_REF", DEF_LINK_REF, ref_link, MAX_WWWBUF);
  GetAppParam ("NCBI", "WWWENTREZ", "LINK_NT", DEF_LINK_NT, nt_link, MAX_WWWBUF);
  GetAppParam ("NCBI", "WWWENTREZ", "LINK_DOC", DEF_LINK_DOC, doc_link, MAX_WWWBUF);
  GetAppParam ("NCBI", "WWWENTREZ", "LINK_EV", DEF_LINK_EV, ev_link, MAX_WWWBUF);
  GetAppParam ("NCBI", "WWWENTREZ", "LINK_EC", DEF_LINK_EC, ec_link, MAX_WWWBUF);
  GetAppParam ("NCBI", "WWWENTREZ", "LINK_TAX", DEF_LINK_TAX, link_tax, MAX_WWWBUF);
  GetAppParam ("NCBI", "WWWENTREZ", "LINK_FF", DEF_LINK_FF, link_ff, MAX_WWWBUF);
  GetAppParam ("NCBI", "WWWENTREZ", "LINK_MUID", DEF_LINK_MUID, link_muid, MAX_WWWBUF);
  GetAppParam ("NCBI", "WWWENTREZ", "LINK_SEQ", DEF_LINK_SEQ, link_seq, MAX_WWWBUF);
  GetAppParam ("NCBI", "WWWENTREZ", "LINK_SEQ", DEF_LINK_SEQ, link_seq, MAX_WWWBUF);
  GetAppParam ("NCBI", "WWWENTREZ", "LINK_SEQ", DEF_LINK_SEQ, link_seq, MAX_WWWBUF);
  GetAppParam ("NCBI", "WWWENTREZ", "LINK_SEQ", DEF_LINK_SEQ, link_seq, MAX_WWWBUF);
  GetAppParam ("NCBI", "WWWENTREZ", "LINK_SEQ", DEF_LINK_SEQ, link_seq, MAX_WWWBUF);
  GetAppParam ("NCBI", "WWWENTREZ", "LINK_SEQ", DEF_LINK_SEQ, link_seq, MAX_WWWBUF);
  GetAppParam("NCBI", "WWWENTREZ", "LINK_FF", DEF_LINK_FF, link_ff, MAX_WWWBUF);
  GetAppParam("NCBI", "WWWENTREZ", "LINK_MUID", DEF_LINK_MUID, link_muid, MAX_WWWBUF);
  GetAppParam("NCBI", "WWWENTREZ", "LINK_ACE", DEF_LINK_ACE, link_ace, MAX_WWWBUF);
  GetAppParam("NCBI", "WWWENTREZ", "LINK_SEQ", DEF_LINK_SEQ, link_seq, MAX_WWWBUF);
  GetAppParam("NCBI", "WWWENTREZ", "LINK_TAX", DEF_LINK_TAX, link_tax, MAX_WWWBUF);
  GetAppParam("NCBI", "WWWENTREZ", "LINK_CODE", DEF_LINK_CODE, link_code, MAX_WWWBUF);
  GetAppParam("NCBI", "WWWENTREZ", "LINK_FLY", DEF_LINK_FLY, link_fly, MAX_WWWBUF);
  GetAppParam("NCBI", "WWWENTREZ", "LINK_COG", DEF_LINK_COG, link_cog, MAX_WWWBUF);
  GetAppParam("NCBI", "WWWENTREZ", "LINK_SGD", DEF_LINK_SGD, link_sgd, MAX_WWWBUF);
  GetAppParam("NCBI", "WWWENTREZ", "LINK_SGD", DEF_LINK_GDB, link_gdb, MAX_WWWBUF);
  GetAppParam("NCBI", "WWWENTREZ", "LINK_CK", DEF_LINK_CK, link_ck, MAX_WWWBUF);
  GetAppParam("NCBI", "WWWENTREZ", "LINK_RICE", DEF_LINK_RICE, link_rice, MAX_WWWBUF);
  GetAppParam("NCBI", "WWWENTREZ", "LINK_SP", DEF_LINK_SP, link_sp, MAX_WWWBUF);
  GetAppParam("NCBI", "WWWENTREZ", "LINK_PDB", DEF_LINK_PDB, link_pdb, MAX_WWWBUF);
  GetAppParam("NCBI", "WWWENTREZ", "LINK_OMIM", DEF_LINK_OMIM, link_omim, MAX_WWWBUF);
  GetAppParam("NCBI", "WWWENTREZ", "LINK_UniSTS", DEF_LINK_UniSTS, link_UniSTS, MAX_WWWBUF);
  GetAppParam("NCBI", "WWWENTREZ", "LINK_dbSTS", DEF_LINK_dbSTS, link_dbSTS, MAX_WWWBUF);
  GetAppParam("NCBI", "WWWENTREZ", "LINK_dbEST", DEF_LINK_dbEST, link_dbEST, MAX_WWWBUF);
  GetAppParam("NCBI", "WWWENTREZ", "LINK_LOCUS", DEF_LINK_LOCUS, link_locus, MAX_WWWBUF);
  GetAppParam("NCBI", "WWWENTREZ", "LINK_SNP", DEF_LINK_SNP, link_snp, MAX_WWWBUF);
  GetAppParam("NCBI", "WWWENTREZ", "LINK_RATMAP", DEF_LINK_RATMAP, link_ratmap, MAX_WWWBUF);
  GetAppParam("NCBI", "WWWENTREZ", "LINK_RGD", DEF_LINK_RGD, link_rgd, MAX_WWWBUF);
  GetAppParam("NCBI", "WWWENTREZ", "LINK_MGD", DEF_LINK_MGD, link_mgd, MAX_WWWBUF);
  GetAppParam("NCBI", "WWWENTREZ", "LINK_FBGN", DEF_LINK_FBGN, link_fly_fbgn, MAX_WWWBUF);
  GetAppParam("NCBI", "WWWENTREZ", "LINK_FBAN", DEF_LINK_FBAN, link_fly_fban, MAX_WWWBUF);
  GetAppParam("NCBI", "WWWENTREZ", "LINK_CDD", DEF_LINK_CDD, link_cdd, MAX_WWWBUF);
  GetAppParam("NCBI", "WWWENTREZ", "LINK_NIAEST", DEF_LINK_NIAEST, link_niaest, MAX_WWWBUF);
  GetAppParam("NCBI", "WWWENTREZ", "LINK_WORM_SEQUENCE", DEF_LINK_WORM_SEQUENCE, link_worm_sequence, MAX_WWWBUF);
  GetAppParam("NCBI", "WWWENTREZ", "LINK_WORM_LOCUS", DEF_LINK_WORM_LOCUS, link_worm_locus, MAX_WWWBUF);
  GetAppParam("NCBI", "WWWENTREZ", "LINK_IMGT", DEF_LINK_IMGT, link_imgt, MAX_WWWBUF);
  GetAppParam("NCBI", "WWWENTREZ", "LINK_IFO", DEF_LINK_IFO, link_ifo, MAX_WWWBUF);
  GetAppParam("NCBI", "WWWENTREZ", "LINK_JCM", DEF_LINK_JCM, link_jcm, MAX_WWWBUF);
  GetAppParam("NCBI", "WWWENTREZ", "LINK_ISFINDER", DEF_LINK_ISFINDER, link_isfinder, MAX_WWWBUF);
  GetAppParam("NCBI", "WWWENTREZ", "LINK_GABI", DEF_LINK_GABI, link_gabi, MAX_WWWBUF);
  GetAppParam("NCBI", "WWWENTREZ", "LINK_FANTOM", DEF_LINK_FANTOM, link_fantom, MAX_WWWBUF);
}


static void FF_www_gcode (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  CharPtr gcode
)
{

  if ( GetWWW(ajp) ) {
    FFAddOneString(ffstring, "<a href=", FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString(ffstring, link_code, FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString(ffstring, "mode=c#SG", FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString(ffstring, gcode, FALSE, FALSE, TILDE_IGNORE);
    FFAddOneChar(ffstring, '>', FALSE);  
    FFAddOneString(ffstring, gcode, FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString(ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
  } else {
    FFAddOneString(ffstring, gcode, FALSE, FALSE, TILDE_IGNORE);
  }
}

static void FF_AddECnumber (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  CharPtr str
)
{
  if ( GetWWW(ajp) ) {
    FFAddOneString(ffstring, "<a href=", FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString(ffstring, ec_link, FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString(ffstring, str, FALSE, FALSE, TILDE_IGNORE);
    FFAddOneChar(ffstring, '>', FALSE);  
    FFAddOneString(ffstring, str, FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString(ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
  } else {
    FFAddOneString(ffstring, str, FALSE, FALSE, TILDE_IGNORE);
  }
}



static void FF_www_featloc(StringItemPtr ffstring, CharPtr loc)
{
  CharPtr ptr;

  if (loc == NULL) return;

  for ( ptr = loc; *ptr != '\0'; ++ptr ) {
    switch (*ptr) {
    case '<' :
      /*FFAddOneString(ffstring, "<", FALSE, FALSE, TILDE_IGNORE);*/
      FFAddOneString(ffstring, "&lt;", FALSE, FALSE, TILDE_IGNORE);
      break;
    case '>' :
      /*FFAddOneString(ffstring, ">", FALSE, FALSE, TILDE_IGNORE);*/
      FFAddOneString(ffstring, "&gt;", FALSE, FALSE, TILDE_IGNORE);
      break;
    default:
      FFAddOneChar(ffstring, *ptr, FALSE);
      break;
    }
  }
}


static void FF_www_db_xref_std (
  StringItemPtr ffstring,
  CharPtr db,
  CharPtr identifier,
  CharPtr link
)
{
  while (*identifier == ' ')
    identifier++;

  FFAddTextToString(ffstring, NULL, db, ":", FALSE, FALSE, TILDE_IGNORE);
  FFAddTextToString(ffstring, "<a href=", link, identifier, FALSE, FALSE, TILDE_IGNORE);
  FFAddTextToString(ffstring, ">", identifier, "</a>", FALSE, FALSE, TILDE_IGNORE);
}

static void FF_www_db_xref_fly (
  StringItemPtr ffstring,
  CharPtr db, 
  CharPtr identifier
)
{
  CharPtr link = link_fly;
  
  if ( StringStr(identifier, "FBa") != NULL ) {
    link = link_fly_fban;
  }
  if ( StringStr(identifier, "FBg") != NULL ) {
    link = link_fly_fbgn;
  }

  FF_www_db_xref_std(ffstring, db, identifier, link);
}


static void FF_www_db_xref_pid(
  StringItemPtr ffstring,
  CharPtr db, 
  CharPtr identifier
)
{
  if ( *identifier != 'g' ) {
    FF_www_db_xref_std(ffstring, db, identifier, link_seq);
    return;
  }
  ++identifier;

  FFAddTextToString(ffstring, NULL, db, ":g", FALSE, FALSE, TILDE_IGNORE);
  FFAddTextToString(ffstring, "<a href=", link_seq, "val=", FALSE, FALSE, TILDE_IGNORE);
  FFAddTextToString(ffstring, identifier, ">", identifier, FALSE, FALSE, TILDE_IGNORE);
  FFAddOneString(ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
}

/*prefix = "<a href=%sval=gnl|dbest|%s>"; */
static void FF_www_db_xref_dbEST(
  StringItemPtr ffstring,
  CharPtr db, 
  CharPtr identifier
)
{
  while (*identifier == ' ')
    identifier++;

  FFAddTextToString(ffstring, NULL, db, ":", FALSE, FALSE, TILDE_IGNORE);
  FFAddTextToString(ffstring, "<a href=", link_dbEST, "val=gnl|dbest|", FALSE, FALSE, TILDE_IGNORE);
  FFAddTextToString(ffstring, identifier, ">", identifier, FALSE, FALSE, TILDE_IGNORE);
  FFAddOneString(ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
}

static void FF_www_db_xref_dbSTS(
  StringItemPtr ffstring,
  CharPtr db, 
  CharPtr identifier
)
{
  while (*identifier == ' ')
    identifier++;

  FFAddTextToString(ffstring, NULL, db, ":", FALSE, FALSE, TILDE_IGNORE);
  FFAddTextToString(ffstring, "<a href=", link_dbSTS, "val=gnl|dbsst|", FALSE, FALSE, TILDE_IGNORE);
  FFAddTextToString(ffstring, identifier, ">", identifier, FALSE, FALSE, TILDE_IGNORE);
  FFAddOneString(ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
}

static void FF_www_db_xref_niaEST(
  StringItemPtr ffstring,
  CharPtr db, 
  CharPtr identifier
)
{
  while (*identifier == ' ')
    identifier++;

  FFAddTextToString(ffstring, NULL, db, ":", FALSE, FALSE, TILDE_IGNORE);
  FFAddTextToString(ffstring, "<a href=", link_niaest, identifier, FALSE, FALSE, TILDE_IGNORE);
  FFAddOneString(ffstring, "&val=1>", FALSE, FALSE, TILDE_IGNORE);
  FFAddOneString(ffstring, identifier, FALSE, FALSE, TILDE_IGNORE);
  FFAddOneString(ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
}


static void FF_www_db_xref_worm(
  StringItemPtr ffstring,
  CharPtr db, 
  CharPtr identifier
)
{
  CharPtr link = NULL;
  CharPtr suffix;

  if ( StringStr (identifier, "unc") != NULL) {
    link = link_worm_locus;
    suffix = ";class=Locus>";
  } else {
    link = link_worm_sequence;
    suffix = ";class=Sequence>";
  }

  while (*identifier == ' ')
    identifier++;

  FFAddTextToString(ffstring, NULL, db, ":", FALSE, FALSE, TILDE_IGNORE);
  FFAddTextToString(ffstring, "<a href=", link, identifier, FALSE, FALSE, TILDE_IGNORE);
  FFAddOneString(ffstring, suffix, FALSE, FALSE, TILDE_IGNORE);
  FFAddOneString(ffstring, identifier, FALSE, FALSE, TILDE_IGNORE);
  FFAddOneString(ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
}


static void FF_www_db_xref_ifo(
  StringItemPtr ffstring,
  CharPtr db, 
  CharPtr identifier
)
{
  while (*identifier == ' ')
    identifier++;

  FFAddTextToString(ffstring, NULL, db, ":", FALSE, FALSE, TILDE_IGNORE);
  FFAddTextToString(ffstring, "<a href=", link_ifo, ">", FALSE, FALSE, TILDE_IGNORE);
  FFAddOneString(ffstring, identifier, FALSE, FALSE, TILDE_IGNORE);
  FFAddOneString(ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
}


static void FF_www_db_xref_gdb(
  StringItemPtr ffstring,
  CharPtr db, 
  CharPtr identifier
)
{
  CharPtr start;
  Char id[20], PNTR idp = id;

  FFAddTextToString(ffstring, NULL, db, ":", FALSE, FALSE, TILDE_IGNORE);

  if ( (start = StringStr(identifier, "G00-")) != NULL ) {
    /* G00-id-id */
    start += StringLen("G00-");
    while ( *start != '\0' ) {
      if ( *start != '-' ) {
        *idp++ = *start++;
      } else {
        *start++;
      }
    }
    *idp = '\0';
    FFAddTextToString(ffstring, "<a href=", link_gdb, id, FALSE, FALSE, TILDE_IGNORE);  
    FFAddTextToString(ffstring, ">", identifier, "</a>", FALSE, FALSE, TILDE_IGNORE);
  } else if ( IS_DIGIT(*identifier) ) {
    /* id */
    FFAddTextToString(ffstring, "<a href=", link_gdb, identifier, FALSE, FALSE, TILDE_IGNORE);  
    FFAddTextToString(ffstring, ">", identifier, "</a>", FALSE, FALSE, TILDE_IGNORE);
  } else {
    FFAddOneString(ffstring, identifier, FALSE, FALSE, TILDE_IGNORE);
  }
}



static void FF_www_db_xref(
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  CharPtr db, CharPtr identifier
)
{
  if ( ffstring == NULL || db == NULL || identifier == NULL ) return;
  
  if ( GetWWW(ajp) ) {
    if ( StringCmp(db, "FLYBASE") == 0) {
      FF_www_db_xref_fly(ffstring, db, identifier);
    } else if ( StringCmp(db , "COG") == 0) {
      FF_www_db_xref_std(ffstring, db, identifier, link_cog);
    } else if ( StringCmp(db , "UniSTS") == 0) {
      FF_www_db_xref_std(ffstring, db, identifier, link_UniSTS);
    } else if ( StringCmp(db , "LocusID") == 0) {
      FF_www_db_xref_std(ffstring, db, identifier, link_locus);
    } else if ( StringCmp(db , "InterimID") == 0) {
      FF_www_db_xref_std(ffstring, db, identifier, link_locus);
    } else if ( StringCmp(db , "MIM") == 0) {
      FF_www_db_xref_std(ffstring, db, identifier, link_omim);
    } else if ( StringCmp(db , "SGD") == 0) {
      FF_www_db_xref_std(ffstring, db, identifier, link_sgd);
    } else if ( StringCmp(db , "IMGT/LIGM") == 0) {
      FF_www_db_xref_std(ffstring, db, identifier, link_imgt);
    } else if ( StringCmp(db , "CK") == 0) {
      FF_www_db_xref_std(ffstring, db, identifier, link_ck);
    } else if ( StringCmp(db , "RiceGenes") == 0) {
      FF_www_db_xref_std(ffstring, db, identifier, link_rice);
    } else if ( StringCmp(db , "dbSNP") == 0) {
      FF_www_db_xref_std(ffstring, db, identifier, link_snp);
    } else if ( StringCmp(db , "RATMAP") == 0) {
      FF_www_db_xref_std(ffstring, db, identifier, link_ratmap);
    } else if ( StringCmp(db , "RGD") == 0) {
      FF_www_db_xref_std(ffstring, db, identifier, link_rgd);
    } else if ( StringCmp(db , "MGD") == 0) {
      FF_www_db_xref_std(ffstring, db, identifier, link_mgd);
    } else if ( StringCmp(db , "CDD") == 0) {
      FF_www_db_xref_std(ffstring, db, identifier, link_cdd);
    } else if ( StringCmp(db , "JCM") == 0) {
      FF_www_db_xref_std(ffstring, db, identifier, link_jcm);
    } else if ( StringCmp(db , "ISFinder") == 0) {
      FF_www_db_xref_std(ffstring, db, identifier, link_isfinder);
    } else if ( StringCmp(db , "GABI") == 0) {
      FF_www_db_xref_std(ffstring, db, identifier, link_gabi);
    } else if ( StringCmp(db , "FANTOM_DB") == 0) {
      FF_www_db_xref_std(ffstring, db, identifier, link_fantom);
    } else if ( StringCmp(db , "PID") == 0) {
      FF_www_db_xref_pid(ffstring, db, identifier);
    } else if ( StringCmp(db , "dbEST") == 0) {
      FF_www_db_xref_dbEST(ffstring, db, identifier);
    } else if ( StringCmp(db , "dbSTS") == 0) {
      FF_www_db_xref_dbSTS(ffstring, db, identifier);
    } else if ( StringCmp(db , "niaEST") == 0) {
      FF_www_db_xref_niaEST(ffstring, db, identifier);
    } else if ( StringCmp(db , "WormBase") == 0) {
      FF_www_db_xref_worm(ffstring, db, identifier);
    } else if ( StringCmp(db , "IFO") == 0) {
      FF_www_db_xref_ifo(ffstring, db, identifier);
    } else if ( StringCmp(db , "GDB") == 0) {
      FF_www_db_xref_gdb(ffstring, db, identifier);

    } else {  
      /* default: no link just the text */
      FFAddTextToString(ffstring, db, ":", identifier, FALSE, FALSE, TILDE_IGNORE);
    }             
  } else { /* not in www mode */
    FFAddTextToString(ffstring, db, ":", identifier, FALSE, FALSE, TILDE_IGNORE);
  }
}

static void FF_www_protein_id(
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  CharPtr seqid
)
{

  if ( GetWWW(ajp) ) {
    FFAddOneString(ffstring, "<a href=", FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString(ffstring, link_seq, FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString(ffstring, "val=", FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString(ffstring, seqid, FALSE, FALSE, TILDE_IGNORE);
    FFAddOneChar(ffstring, '>', FALSE);  
    FFAddOneString(ffstring, seqid, FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString(ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
  } else {
    FFAddOneString(ffstring, seqid, FALSE, FALSE, TILDE_IGNORE);
  }
}

static void  FF_www_muid(
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  Int4 muid
)
{
  Char numbuf[40];
  
  if ( GetWWW(ajp) ) {
    FFAddTextToString(ffstring, "<a href=", link_muid, NULL, FALSE, FALSE, TILDE_IGNORE);
    sprintf(numbuf, "%ld", (long)muid);
    FFAddTextToString(ffstring, "uid=", numbuf, "&form=6&db=m&Dopt=r>", FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString(ffstring, numbuf, FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString(ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
  } else {
    sprintf(numbuf, "%ld", (long)muid);
    FFAddOneString(ffstring, numbuf, FALSE, FALSE, TILDE_IGNORE);
  }
}

static void FF_www_accession (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  CharPtr cstring
)
{
  if (cstring == NULL || ffstring == NULL) return;

  if ( GetWWW(ajp) ) {
    FFAddTextToString(ffstring, "<a href=", link_seq, NULL, FALSE, FALSE, TILDE_IGNORE);
    FFAddTextToString(ffstring, "val=", cstring, ">", FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString(ffstring, cstring, FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString(ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
  } else {
    FFAddOneString(ffstring, cstring, FALSE, FALSE, TILDE_IGNORE);
  }
  return;
}

static CharPtr TxtSave (CharPtr text, size_t len)

{
   CharPtr str = NULL;

   if ((text == NULL) || (len == 0))
      return str;

   str = MemNew((size_t)(len + 1));
   MemCopy(str, text, (size_t)len);

   return (str);
}

static Boolean FF_www_dbsource(
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  CharPtr str,
  Boolean first,
  Uint1 choice
)
{
  CharPtr  temp, end, text, loc, link=NULL;
  Int2 j;

  if( GetWWW(ajp) ) {
    if (choice == SEQID_PIR /*|| choice == SEQID_SWISSPROT */) {
      link = link_seq;
    } else if (choice == SEQID_PDB || choice == SEQID_PRF) {
      link = link_seq;
    } else if (choice == SEQID_EMBL || choice == SEQID_GENBANK || 
        choice == SEQID_DDBJ || choice == SEQID_GIBBSQ || 
        choice == SEQID_GIBBMT || choice == SEQID_GI || 
        choice == SEQID_GIIM || choice == SEQID_OTHER ||
        choice == SEQID_TPG || choice == SEQID_TPE || choice == SEQID_TPD)  {
      link = link_seq;
    } else {
      AddStringWithTildes(ffstring, str);
      return TRUE;
    }
  
    if ((text = StringStr(str, "accession")) != NULL) {
      end = text + 9;
      j = 9;
      while (*end == ' ') {
        ++end;
        j++;
      }
      if (first == FALSE) {
        FFAddOneString(ffstring, ", ", FALSE, FALSE, TILDE_IGNORE);
      }
      loc = TxtSave (str, end-str - j);
      FFAddOneString(ffstring, loc, FALSE, FALSE, TILDE_IGNORE);
      MemFree (loc);
      for (; text != end; ++text ) {
        FFAddOneChar(ffstring, *text, FALSE);
      }

      temp = text;
      end += StringLen(text) - 1;
      if ( *end != ';' ) {
        ++end;
      }

      FFAddTextToString(ffstring, "<a href=", link, "val=", FALSE, FALSE, TILDE_IGNORE);
      for (text = temp; text != end; ++text ) {
        FFAddOneChar(ffstring, *text, FALSE);
      }
      FFAddOneString(ffstring, ">", FALSE, FALSE, TILDE_IGNORE);

      text = temp;
      FFAddOneString(ffstring, text, FALSE, FALSE, TILDE_IGNORE);
      FFAddOneString(ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
    } else {
      if (first == FALSE) {
        FFAddOneString(ffstring, ", ", FALSE, FALSE, TILDE_IGNORE);
      }
      FFAddOneString(ffstring, str, FALSE, FALSE, TILDE_IGNORE);
    }
  } else {
    AddStringWithTildes(ffstring, str);
  }
  return TRUE;
}


/* old utility functions */

static ValNodePtr ValNodeCopyStrToHead (ValNodePtr PNTR head, Int2 choice, CharPtr str)

{
  ValNodePtr newnode;

  if (head == NULL || str == NULL) return NULL;

  newnode = ValNodeNew (NULL);
  if (newnode == NULL) return NULL;

  newnode->choice = (Uint1) choice;
  newnode->data.ptrvalue = StringSave (str);

  newnode->next = *head;
  *head = newnode;

  return newnode;
}

/* the val node strings mechanism will be replaced by a more efficient method later  */

static CharPtr MergeValNodeStrings (
  ValNodePtr list
)

{
  size_t      len;
  CharPtr     ptr;
  CharPtr     str;
  CharPtr     tmp;
  ValNodePtr  vnp;


  if (list == NULL) return NULL;

  for (vnp = list, len = 0; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    len += StringLen (str);
  }
  if (len == 0) return NULL;

  ptr = MemNew (sizeof (Char) * (len + 2));
  if (ptr == NULL) return NULL;

  for (vnp = list, tmp = ptr; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    tmp = StringMove (tmp, str);
  }

  return ptr;
}


static void AddValNodeString (
  ValNodePtr PNTR head,
  CharPtr prefix,
  CharPtr string,
  CharPtr suffix
)

{
  Char     buf [256];
  CharPtr  freeme = NULL;
  size_t   len;
  CharPtr  newstr;
  CharPtr  strptr;

  len = StringLen (prefix) + StringLen (string) + StringLen (suffix);
  if (len == 0) return;

  if (len < sizeof (buf)) {

    /* if new string fits in stack buffer, no need to allocate */

    MemSet ((Pointer) buf, 0, sizeof (buf));
    newstr = buf;

  } else {

    /* new string bigger than stack buffer, so allocate sufficient string */

    newstr = (CharPtr) MemNew (sizeof (Char) * (len + 2));
    if (newstr == NULL) return;

    /* allocated string will be freed at end of function */

    freeme = newstr;
  }

  strptr = newstr;

  if (prefix != NULL) {
    strptr = StringMove (strptr, prefix);
  }

  if (string != NULL) {
    strptr = StringMove (strptr, string);
  }

  if (suffix != NULL) {
    strptr = StringMove (strptr, suffix);
  }

  /* currently just makes a valnode list, to be enhanced later */

  ValNodeCopyStr (head, 0, newstr);

  /* if large string was allocated, free it now */

  if (freeme != NULL) {
    MemFree (freeme);
  }
}


static void FFAddString_NoRedund (
  StringItemPtr unique,
  CharPtr prefix,
  CharPtr string,
  CharPtr suffix
)
{
  CharPtr    str = string;
  Int4       foundPos = 0;
  Boolean    wholeWord = FALSE;

  if ( StringHasNoText(prefix)  &&
       StringHasNoText(string)  &&
       StringHasNoText(suffix)  ) return;

  if (StringNICmp (string, "tRNA-", 5) == 0) {
    str = string+5;
	}

  while ( foundPos >= 0 && !wholeWord ) {
    foundPos = FFStringSearch(unique, str, foundPos);
    if ( foundPos >= 0 ) {
      wholeWord = IsWholeWordSubstr(unique, foundPos, str);
      foundPos += StringLen(str);
    }
  }

  if ( foundPos < 0 || !wholeWord ) {
      FFAddTextToString(unique, prefix, string, suffix, FALSE, FALSE, TILDE_IGNORE);
  }
}



/* s_AddPeriodToEnd () -- Adds a '.' to the end of a given string if */
/*                        there is not already one there.            */
/*                                                                   */
/*                        Note that this adds one character to the   */
/*                        length of the string, leading to a         */
/*                        memory overrun if space was not previously */
/*                        allocated for this.                        */

static void s_AddPeriodToEnd (CharPtr someString)
{
  Int4  len;

  if (StringHasNoText (someString)) return;
  len = StringLen (someString);
  if (len < 1) return;
  if (someString[len-1] != '.')
    {
      someString[len] = '.';
      someString[len+1] = '\0';
    }
}

/* s_RemovePeriodFromEnd () -- If the last character in a given      */
/*                             string is a '.', removes it.          */

static Boolean s_RemovePeriodFromEnd (CharPtr someString)
{
  Int4  len;

  if (StringHasNoText (someString)) return FALSE;
  len = StringLen (someString);
  if (len < 1) return FALSE;
  if (someString[len-1] == '.') {
    someString[len-1] = '\0';
    return TRUE;
  }
  return FALSE;
}

/**/
/*   isEllipsis () - Determines if a string ends in an ellipses */
/**/

static Boolean IsEllipsis (
  CharPtr str
)

{
  size_t   len;
  CharPtr  ptr;

  if (StringHasNoText (str)) return FALSE;
  len = StringLen (str);
  if (len < 3) return FALSE;
  ptr = str + len - 3;
  return (Boolean) (ptr [0] == '.' && ptr [1] == '.' && ptr [2] == '.');
}

static void A2GBSeqLocReplaceID (
  SeqLocPtr newloc,
  SeqLocPtr ajpslp
)

{
  BioseqPtr  bsp;
  SeqIdPtr   sip;

  bsp = BioseqFindFromSeqLoc (ajpslp);
  if (bsp == NULL) return;
  sip = SeqIdFindBest (bsp->id, 0);
  SeqLocReplaceID (newloc, sip);
}

static CharPtr asn2gb_PrintDate (
  DatePtr dp
)

{
  Char    buf [30];
  size_t  len;

  if (dp == NULL) return NULL;

  if (DatePrint (dp, buf)) {
    if (StringICmp (buf, "Not given") != 0) {
      len = StringLen (buf);
      if (len > 0) {
        if (buf [len - 1] == '\n') {
          if (buf [len - 2] == '.') {
            buf [len - 2] = '\0';
          } else {
            buf [len - 1] = '\0';
          }
        }
      }
      return StringSave (buf);
    }
  }

  return NULL;
}

static CharPtr month_names [] = {
  "JAN", "FEB", "MAR", "APR", "MAY", "JUN",
  "JUL", "AUG", "SEP", "OCT", "NOV", "DEC",
  "??"
};

static CharPtr DateToGB (
  CharPtr buf,
  DatePtr dp,
  Boolean citSub
)

{
  Int2  day;
  Int2  month;
  Int2  year;

  if (buf != NULL) {
    *buf = '\0';
  }
  if (dp == NULL) return NULL;

  if (dp->data [0] == 0) {

    StringCpy (buf, dp->str);

  } else if (dp->data [0] == 1) {

    year = 1900 + (Int2) dp->data [1];
    month = (Int2) dp->data [2];
    day = (Int2) dp->data [3];

    if (citSub) {
      if (month < 1 || month > 12) {
        month = 13;
      }
      if (day < 1 || day > 31) {
        day = 0;
      }
    } else {
      if (month < 1 || month > 12) {
        month = 1;
      }
      if (day < 1 || day > 31) {
        day = 1;
      }
    }

    if (day < 1) {
      sprintf (buf, "??-%s-%ld",
               month_names [month-1], (long) year);
    } else if (day < 10) {
      sprintf (buf, "0%ld-%s-%ld",
               (long) day, month_names [month-1], (long) year);
    } else {
      sprintf(buf, "%ld-%s-%ld",
               (long) day, month_names [month-1], (long) year);
    }
  }

  return buf;
}


/* ********************************************************************** */

/* format functions allocate printable string for given paragraph */

static CharPtr DefaultFormatBlock (
  Asn2gbFormatPtr afp,
  BaseBlockPtr bbp
)

{
  if (afp == NULL || bbp == NULL) return NULL;

  /* default format function assumes string pre-allocated by add block function */

  return StringSaveNoNull (bbp->string);
}

/* superset of http://www.ncbi.nlm.nih.gov/collab/db_xref.html and RefSeq db_xrefs */

static CharPtr legalDbXrefs [] = {
  "PIDe", "PIDd", "PIDg", "PID",
  "ATCC",
  "ATCC(in host)",
  "ATCC(dna)",
  "BDGP_EST",
  "BDGP_INS",
  "CDD",
  "CK",
  "COG",
  "dbEST",
  "dbSNP",
  "dbSTS",
  "ENSEMBL",
  "ESTLIB",
  "FANTOM_DB",
  "FLYBASE",
  "GABI",
  "GDB",
  "GeneID",
  "GI",
  "GO",
  "IFO",
  "IMGT/LIGM",
  "IMGT/HLA",
  "InterimID",
  "ISFinder",
  "JCM",
  "LocusID",
  "MaizeDB",
  "MGD",
  "MGI",
  "MIM",
  "niaEST",
  "PIR",
  "PSEUDO",
  "RATMAP",
  "RiceGenes",
  "REMTREMBL",
  "RGD",
  "RZPD",
  "SGD",
  "SoyBase",
  "SPTREMBL",
  "SWISS-PROT",
  "taxon",
  "UniGene",
  "UniSTS",
  "WormBase",
  NULL
};

static CharPtr legalRefSeqDbXrefs [] = {
  NULL
};

static CharPtr organellePrefix [] = {
  NULL,
  NULL,
  "Chloroplast ",
  "Chromoplast ",
  "Kinetoplast ",
  "Mitochondrion ",
  "Plastid ",
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  "Cyanelle ",
  NULL,
  NULL,
  "Nucleomorph ",
  "Apicoplast ",
  "Leucoplast ",
  "Proplastid ",
  NULL
};

static CharPtr newOrganellePrefix [] = {
  NULL,
  NULL,
  "chloroplast ",
  "chromoplast ",
  "kinetoplast ",
  "mitochondrion ",
  "plastid ",
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  "cyanelle ",
  NULL,
  NULL,
  "nucleomorph ",
  "apicoplast ",
  "leucoplast ",
  "proplastid ",
  NULL
};

static CharPtr FormatSourceBlock (
  Asn2gbFormatPtr afp,
  BaseBlockPtr bbp
)

{
  CharPtr            acr = NULL;
  Boolean            addPeriod = TRUE;
  IntAsn2gbJobPtr    ajp;
  CharPtr            ana = NULL;
  Asn2gbSectPtr      asp;
  BioSourcePtr       biop = NULL;
  CharPtr            com = NULL;
  CharPtr            common = NULL;
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext;
  CharPtr            gbacr = NULL;
  CharPtr            gbana = NULL;
  GBBlockPtr         gbp = NULL;
  GBSeqPtr           gbseq;
  CharPtr            gbsyn = NULL;
  Uint1              genome;
  ValNodePtr         mod = NULL;
  OrgModPtr          omp = NULL;
  OrgNamePtr         onp;
  CharPtr            organelle = NULL;
  OrgRefPtr          orp;
  CharPtr            prefix = " (";
  SeqDescrPtr        sdp;
  CharPtr            second = NULL;
  SeqFeatPtr         sfp;
  CharPtr            str;
  CharPtr            syn = NULL;
  CharPtr            taxname = NULL;
  Boolean            using_anamorph = FALSE;
  StringItemPtr      ffstring, temp;

  if (afp == NULL || bbp == NULL) return NULL;
  ajp = afp->ajp;
  if (ajp == NULL) return NULL;
  asp = afp->asp;
  if (asp == NULL) return NULL;

  if (! StringHasNoText (bbp->string)) return StringSave (bbp->string);

  ffstring = FFGetString(ajp);
  if ( ffstring == NULL ) return NULL;

  if (bbp->itemtype == OBJ_SEQDESC) {
    sdp = SeqMgrGetDesiredDescriptor (bbp->entityID, NULL, bbp->itemID, 0, NULL, &dcontext);
    if (sdp != NULL) {
      if (dcontext.seqdesctype == Seq_descr_source) {
        biop = (BioSourcePtr) sdp->data.ptrvalue;
      } else if (dcontext.seqdesctype == Seq_descr_genbank) {
        gbp = (GBBlockPtr) sdp->data.ptrvalue;
      }
    }
  } else if (bbp->itemtype == OBJ_SEQFEAT) {
    sfp = SeqMgrGetDesiredFeature (bbp->entityID, NULL, bbp->itemID, 0, NULL, &fcontext);
    if (sfp != NULL && fcontext.seqfeattype == SEQFEAT_BIOSRC) {
      biop = (BioSourcePtr) sfp->data.value.ptrvalue;
    }
  }
  if (gbp != NULL) {
    common = gbp->source;
  }

  if (biop != NULL) {
    genome = biop->genome;
    if (genome <= 19) {
      if (ajp->newSourceOrg && (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT)) {
        organelle = newOrganellePrefix [genome];
      } else {
        organelle = organellePrefix [genome];
      }
    }
    orp = biop->org;
    if (orp != NULL) {
      taxname = orp->taxname;
      common = orp->common;
      mod = orp->mod;
      onp = orp->orgname;
      if (onp != NULL) {

        if (ajp->newSourceOrg) {
          for (omp = onp->mod; omp != NULL; omp = omp->next) {
            switch (omp->subtype) {
              case ORGMOD_common :
                com = omp->subname;
                break;
              case ORGMOD_acronym :
                acr = omp->subname;
                break;
              case ORGMOD_synonym :
                syn = omp->subname;
                break;
              case ORGMOD_anamorph :
                ana = omp->subname;
                break;
              case ORGMOD_gb_acronym :
                gbacr = omp->subname;
                break;
              case ORGMOD_gb_anamorph :
                gbana = omp->subname;
                break;
              case ORGMOD_gb_synonym :
                gbsyn = omp->subname;
                break;
              default :
                break;
            }
          }

          if (StringHasNoText (second)) {
            second = syn;
            using_anamorph = FALSE;
          }
           if (StringHasNoText (second)) {
             second = acr;
             using_anamorph = FALSE;
          }
          if (StringHasNoText (second)) {
            second = ana;
            using_anamorph = TRUE;
          }
          if (StringHasNoText (second)) {
            second = com;
            using_anamorph = FALSE;
          }

          if (StringHasNoText (second)) {
            second = gbsyn;
            using_anamorph = FALSE;
          }
          if (StringHasNoText (second)) {
            second = gbacr;
            using_anamorph = FALSE;
          }
          if (StringHasNoText (second)) {
            second = gbana;
            using_anamorph = TRUE;
          }

          if (StringHasNoText (second)) {
            second = common;
            using_anamorph = FALSE;
          }
          if (using_anamorph) {
            prefix = " (anamorph: ";
          }
        }
      }
    }
  }

  /* If the organelle prefix is already on the */
  /* name, don't add it.                       */

  if (StringNCmp (organelle, taxname, StringLen (organelle)) == 0)
    organelle = "";

  if (StringHasNoText (common)) {
    common = taxname;
  }
  if (StringHasNoText (common)) {
    common = "Unknown.";
  }
  if (StringHasNoText (taxname)) {
    taxname = "Unknown.";
  }

  if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
    
    temp = FFGetString(ajp);

    if (ajp->newSourceOrg) {

      if (! StringHasNoText (organelle)) {
        FFAddTextToString(temp, NULL, organelle, NULL, FALSE, FALSE, TILDE_IGNORE);
      }
      FFAddTextToString(temp, NULL, taxname, NULL, FALSE, FALSE, TILDE_IGNORE);
      if (! StringHasNoText (second)) {
        FFAddTextToString(temp, prefix, second, ")", FALSE, FALSE, TILDE_IGNORE);
      }
      addPeriod = FALSE;

    } else {
      FFAddTextToString(temp, NULL, common, NULL, FALSE, FALSE, TILDE_IGNORE);
      while (mod != NULL) {
        str = (CharPtr) mod->data.ptrvalue;
        if (! StringHasNoText (str)) {
          FFAddTextToString(temp, " ", str, NULL, FALSE, FALSE, TILDE_IGNORE);
        }
        mod = mod->next;
      }
    }

    str = FFToCharPtr(temp);
    if (StringCmp (str, ".") == 0) {
      str = MemFree (str);
    }
    FFRecycleString(ajp, temp);
    /* optionally populate gbseq for XML-ized GenBank format */

    if (ajp->gbseq) {
      gbseq = &asp->gbseq;
    } else {
      gbseq = NULL;
    }

    if (gbseq != NULL) {
      gbseq->source = StringSave (str);
    }

    
    FFStartPrint(ffstring, afp->format, 0, 12, "SOURCE", 12, 5, 5, "OS", TRUE);
    if (str != NULL) {
      FFAddTextToString(ffstring, NULL, str, NULL, addPeriod, FALSE, TILDE_TO_SPACES);
    } else {
      FFAddOneChar(ffstring, '.', FALSE);
    }
    
    MemFree (str);

  } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {

    FFStartPrint(ffstring, afp->format, 0, 12, "SOURCE", 12, 5, 5, "OS", TRUE);
    FFAddTextToString(ffstring, organelle, taxname, NULL, FALSE, FALSE, TILDE_TO_SPACES);
    FFAddTextToString(ffstring, " (", common, ")", FALSE, FALSE, TILDE_TO_SPACES);
    
  }
  
  str = FFEndPrint(ajp, ffstring, afp->format, 12, 12, 0, 5, "OS");
  FFRecycleString(ajp, ffstring);
  return str;
}

static CharPtr FormatOrganismBlock (
  Asn2gbFormatPtr afp,
  BaseBlockPtr bbp
)

{
  IntAsn2gbJobPtr    ajp;
  Asn2gbSectPtr      asp;
  BioSourcePtr       biop = NULL;
  CharPtr            common = NULL;
  DbtagPtr           dbt;
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext;
  GBSeqPtr           gbseq;
  Uint1              genome;
  CharPtr            lineage = NULL;
  ObjectIdPtr        oip;
  OrgNamePtr         onp;
  CharPtr            organelle = NULL;
  OrgRefPtr          orp;
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp;
  CharPtr            str;
  Int4               taxid = -1;
  CharPtr            taxname = NULL;
  ValNodePtr         vnp;
  StringItemPtr      ffstring, temp;
  Char               buf [16];

  if (afp == NULL || bbp == NULL) return NULL;
  ajp = afp->ajp;
  if (ajp == NULL) return NULL;
  asp = afp->asp;
  if (asp == NULL) return NULL;


  if (! StringHasNoText (bbp->string)) return StringSave (bbp->string);

  if (bbp->itemtype == OBJ_SEQDESC) {
    sdp = SeqMgrGetDesiredDescriptor (bbp->entityID, NULL, bbp->itemID, 0, NULL, &dcontext);
    if (sdp != NULL && dcontext.seqdesctype == Seq_descr_source) {
      biop = (BioSourcePtr) sdp->data.ptrvalue;
    }
  } else if (bbp->itemtype == OBJ_SEQFEAT) {
    sfp = SeqMgrGetDesiredFeature (bbp->entityID, NULL, bbp->itemID, 0, NULL, &fcontext);
    if (sfp != NULL && fcontext.seqfeattype == SEQFEAT_BIOSRC) {
      biop = (BioSourcePtr) sfp->data.value.ptrvalue;
    }
  }
  if (biop != NULL) {
    genome = biop->genome;
    if (genome <= 19) {
      organelle = organellePrefix [genome];
    }
    orp = biop->org;
    if (orp != NULL) {
      taxname = orp->taxname;
      common = orp->common;
      onp = orp->orgname;
      if (onp != NULL) {
        lineage = onp->lineage;
      }
      for (vnp = orp->db; vnp != NULL; vnp = vnp->next) {
        dbt = (DbtagPtr) vnp->data.ptrvalue;
        if (dbt == NULL) continue;
        if (StringCmp (dbt->db, "taxon") == 0) {
          oip = dbt->tag;
          if (oip != NULL) {
            taxid = oip->id;
          }
        }
      }
    }
  }

  /* If the organelle prefix is already on the */
  /* name, don't add it.                       */

  if (StringNCmp (organelle, taxname, StringLen (organelle)) == 0)
    organelle = "";

  if (StringHasNoText (common)) {
    common = taxname;
  }
  if (StringHasNoText (common)) {
    common = "Unknown.";
  }
  if (StringHasNoText (taxname)) {
    taxname = "Unknown.";
  }
  if (StringHasNoText (lineage)) {
    lineage = "Unclassified.";
  }

  ffstring = FFGetString(ajp);
  temp = FFGetString(ajp);
  if ( ffstring == NULL || temp == NULL ) return NULL;

  if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
    
    FFStartPrint(temp, afp->format, 2, 12, "ORGANISM", 12, 5, 5, "OC", FALSE);
    if (! ajp->newSourceOrg) {
      FFAddOneString(temp, organelle, FALSE, FALSE, TILDE_IGNORE);
    }
    if (StringNICmp (taxname, "Unknown", 7) != 0) {
      if ( GetWWW(ajp) ) { 
        if (taxid != -1) {
          FFAddOneString(temp, "<a href=", FALSE, FALSE, TILDE_IGNORE);
          FFAddOneString(temp, link_tax, FALSE, FALSE, TILDE_IGNORE);
          FFAddOneString(temp, "id=", FALSE, FALSE, TILDE_IGNORE);
          sprintf (buf, "%ld", (long) taxid);
          FFAddOneString(temp, buf, FALSE, FALSE, TILDE_IGNORE);
          FFAddOneString(temp, ">", FALSE, FALSE, TILDE_IGNORE);
        } else {
          FFAddOneString(temp, "<a href=", FALSE, FALSE, TILDE_IGNORE);
          FFAddOneString(temp, link_tax, FALSE, FALSE, TILDE_IGNORE);
          FFAddOneString(temp, "name=", FALSE, FALSE, TILDE_IGNORE);
          sprintf (buf, "%ld", (long) taxid);
          FFAddOneString(temp, taxname, FALSE, FALSE, TILDE_IGNORE);
          FFAddOneString(temp, ">", FALSE, FALSE, TILDE_IGNORE);
        }
        FFAddOneString(temp, taxname, FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString(temp, "</a>", FALSE, FALSE, TILDE_IGNORE);
      } else {
        FFAddOneString(temp, taxname, FALSE, FALSE, TILDE_IGNORE);
      }
    } else {
      FFAddOneString(temp, taxname, FALSE, FALSE, TILDE_IGNORE);
    }
    FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
    FFRecycleString(ajp, temp);

    temp = FFGetString(ajp);
    FFStartPrint(temp, afp->format, 12, 12, NULL, 0, 5, 5, "OC", FALSE);
    FFAddTextToString(temp, NULL, lineage, NULL, TRUE, FALSE, TILDE_TO_SPACES);
    FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
    FFRecycleString(ajp, temp);
    /* optionally populate gbseq for XML-ized GenBank format */

    if (ajp->gbseq) {
      gbseq = &asp->gbseq;
    } else {
      gbseq = NULL;
    }

    if (gbseq != NULL) {
      temp = FFGetString(ajp);
      if (! ajp->newSourceOrg) {
        FFAddOneString(temp, organelle, FALSE, FALSE, TILDE_IGNORE);
      }
      FFAddOneString(temp, taxname, FALSE, FALSE, TILDE_IGNORE);
      gbseq->organism = FFToCharPtr(temp);
      gbseq->taxonomy = StringSave (lineage);
      FFRecycleString(ajp, temp);
    }

  } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {
    FFStartPrint(temp, afp->format, 12, 12, NULL, 0, 5, 5, "OC", FALSE);
    FFAddTextToString(temp, NULL, lineage, NULL, TRUE, FALSE, TILDE_TO_SPACES);
    FFLineWrap(ffstring, temp, 5, 5, ASN2FF_EMBL_MAX, "OC");
  }

  
  str = FFToCharPtr(ffstring);
  FFRecycleString(ajp, ffstring);
  return str;
}

/* format references section */

static AuthListPtr GetAuthListPtr (
  PubdescPtr pdp,
  CitSubPtr csp
)

{
  AuthListPtr  alp = NULL;
  CitArtPtr    cap;
  CitBookPtr   cbp;
  CitGenPtr    cgp;
  CitPatPtr    cpp;
  ValNodePtr   vnp;

  if (csp != NULL) {
    alp = csp->authors;
    if (alp != NULL) return alp;
  }
  if (pdp == NULL) return NULL;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    switch (vnp->choice) {
      case PUB_Gen :
        cgp = (CitGenPtr) vnp->data.ptrvalue;
        if (cgp != NULL) {
          alp = cgp->authors;
        }
        break;
      case PUB_Sub :
        csp = (CitSubPtr) vnp->data.ptrvalue;
        if (csp != NULL) {
          alp = csp->authors;
        }
        break;
      case PUB_Article :
        cap = (CitArtPtr) vnp->data.ptrvalue;
        if (cap != NULL) {
          alp = cap->authors;
        }
        break;
      case PUB_Book :
      case PUB_Proc :
      case PUB_Man :
        cbp = (CitBookPtr) vnp->data.ptrvalue;
        if (cbp != NULL) {
          alp = cbp->authors;
        }
        break;
      case PUB_Patent :
        cpp = (CitPatPtr) vnp->data.ptrvalue;
        if (cpp != NULL) {
          alp = cpp->authors;
        }
        break;
      default :
        break;
    }

    if (alp != NULL) return alp;
  }

  return NULL;
}

static CharPtr MakeSingleAuthorString (
  FmtType format,
  CharPtr prefix,
  CharPtr name,
  CharPtr initials,
  CharPtr suffix,
  IndxPtr index,
  GBReferencePtr gbref
)

{
  Char     ch;
  Char     dummy [10];
  size_t   len;
  CharPtr  nametoindex;
  CharPtr  ptr;
  CharPtr  str;
  CharPtr  tmp;

  if (name == NULL) return NULL;

  /* !!! clean up 'et al' as (presumably) last author !!! */

  /* !!! temporary to suppress diff !!! */
  {
  if (StringLen (name) <= 6 &&
      (StringNICmp (name, "et al", 5) == 0 || StringNICmp (name, "et,al", 5) == 0)) {
    if (StringCmp (prefix, " and ") == 0) {
      prefix = NULL;
      dummy [0] = ' ';
      StringNCpy_0 (dummy + 1, name, sizeof (dummy) - 1);
      name = dummy;
    }
  }
  }
  /*
  if (StringLen (name) <= 6 &&
      (StringNICmp (name, "et al", 5) == 0 || StringNICmp (name, "et,al", 5) == 0)) {
    name = "et al.";
    if (StringCmp (prefix, " and ") == 0) {
      prefix = ", ";
    }
  }
  */

  len = StringLen (name) + StringLen (initials) + StringLen (suffix) + StringLen (prefix);
  str = MemNew (sizeof (Char) * (len + 4));
  if (str == NULL) return NULL;

  ptr = str;
  if (! StringHasNoText (prefix)) {
    ptr = StringMove (ptr, prefix);
  }
  nametoindex = ptr;

  /* initials and suffix to support structured name fields */

  tmp = StringMove (ptr, name);
  if (! StringHasNoText (initials)) {
    tmp = StringMove (tmp, ",");
    tmp = StringMove (tmp, initials);
  }
  if (! StringHasNoText (suffix)) {
    tmp = StringMove (tmp, " ");
    tmp = StringMove (tmp, suffix);
  }

  /* optionally populate indexes for NCBI internal database */

  if (index != NULL) {
    ValNodeCopyStrToHead (&(index->authors), 0, nametoindex);
  }

  /* optionally populate gbseq for XML-ized GenBank format */

  if (gbref != NULL) {
    ValNodeCopyStr (&(gbref->authors), 0, nametoindex);
  }

  /* if embl, remove commas in individual names, starting after prefix */

  if (format == EMBL_FMT || format == EMBLPEPT_FMT) {
    tmp = ptr;
    ch = *tmp;
    while (ch != '\0') {
      if (ch == ',') {
        *tmp = ' ';
      }
      tmp++;
      ch = *tmp;
    }
  }

  return str;
}

static CharPtr GetAuthorsString (
  FmtType format,
  AuthListPtr alp,
  CharPtr PNTR consortP,
  IndxPtr index,
  GBReferencePtr gbref
)

{
  AuthorPtr    ap;
  ValNodePtr   conslist;
  Int2         count;
  ValNodePtr   head = NULL;
  ValNodePtr   names;
  ValNodePtr   next;
  NameStdPtr   nsp;
  PersonIdPtr  pid;
  ValNodePtr   pidlist;
  CharPtr      prefix = NULL;
  CharPtr      str;
  ValNodePtr   vnp;

  if (alp == NULL) return NULL;

  alp = AsnIoMemCopy ((Pointer) alp,
                      (AsnReadFunc) AuthListAsnRead,
                      (AsnWriteFunc) AuthListAsnWrite);
  if (alp == NULL) return NULL;

  count = 0;
  if (alp->choice == 1) {

    pidlist = NULL;
    conslist = NULL;

    for (names = alp->names; names != NULL; names = names->next) {
      ap = (AuthorPtr) names->data.ptrvalue;
      if (ap == NULL) continue;
      pid = ap->name;
      if (pid == NULL) continue;
      if (pid->choice == 2 || pid->choice == 3 || pid->choice == 4) {
        ValNodeAddPointer (&pidlist, 0, (Pointer) pid);
      } else if (pid->choice == 5) {
        ValNodeAddPointer (&conslist, 0, (Pointer) pid);
      }
    }

    for (vnp = pidlist; vnp != NULL; vnp = vnp->next) {
      next = vnp->next;
      if (next == NULL) {
        if (format == GENBANK_FMT || format == GENPEPT_FMT) {
          if (count == 0) {
            prefix = NULL;
          } else {
            prefix = " and ";
          }
        }
      }
      str = NULL;
      pid = (PersonIdPtr) vnp->data.ptrvalue;
      if (pid->choice == 2) {
        nsp = (NameStdPtr) pid->data;
        if (nsp != NULL) {
          if (! StringHasNoText (nsp->names [0])) {
            str = MakeSingleAuthorString (format, prefix, nsp->names [0], nsp->names [4], nsp->names [5], index, gbref);
          } else if (! StringHasNoText (nsp->names [3])) {
            str = MakeSingleAuthorString (format, prefix, nsp->names [3], NULL, NULL, index, gbref);
          }
        }
      } else if (pid->choice == 3 || pid->choice == 4) {
        str = MakeSingleAuthorString (format, prefix, (CharPtr) pid->data, NULL, NULL, index, gbref);
      }
      if (str != NULL) {
        ValNodeAddStr (&head, 0, str);
        count++;
      }
      prefix = ", ";
    }

    for (vnp = conslist; vnp != NULL; vnp = vnp->next) {
      str = NULL;
      pid = (PersonIdPtr) vnp->data.ptrvalue;
      if (pid->choice == 5) {
        str = MakeSingleAuthorString (format, NULL, (CharPtr) pid->data, NULL, NULL, index, NULL);
        if ((! StringHasNoText (str)) && consortP != NULL && *consortP == NULL) {
          *consortP = StringSave (str);
        }

        /* optionally populate gbseq for XML-ized GenBank format */

        if (gbref != NULL) {
          gbref->consortium = StringSave (str);
        }

        str = MemFree (str);
      }
    }

    ValNodeFree (pidlist);
    ValNodeFree (conslist);

  } else if (alp->choice == 2 || alp->choice == 3) {
    for (vnp = alp->names; vnp != NULL; vnp = vnp->next) {
      next = vnp->next;
      if (next == NULL) {
        if (format == GENBANK_FMT || format == GENPEPT_FMT) {
          if (count == 0) {
            prefix = NULL;
          } else {
            prefix = " and ";
          }
        }
      }
      str = MakeSingleAuthorString (format, prefix, (CharPtr) vnp->data.ptrvalue, NULL, NULL, index, gbref);
      if (str != NULL) {
        ValNodeAddStr (&head, 0, str);
        count++;
      }
      prefix = ", ";
    }
  }

  str = MergeValNodeStrings (head);

  ValNodeFreeData (head);

  AuthListFree (alp);

  return str;
}

/*
Strips all spaces in string in following manner. If the function
meet several spaces (spaces and tabs) in succession it replaces them
with one space. Strips all spaces after '(' and before ')'
*/

static void StrStripSpaces (
  CharPtr str
)

{
  CharPtr  new_str;

  if (str == NULL) return;

  new_str = str;
  while (*str != '\0') {
    *new_str++ = *str;
    if (*str == ' ' || *str == '\t' || *str == '(') {
      for (str++; *str == ' ' || *str == '\t'; str++) continue;
      if (*str == ')' || *str == ',') {
        new_str--;
      }
    } else {
      str++;
    }
  }
  *new_str = '\0';
}

static Boolean AllCaps (
  CharPtr p
)

{
  if (p == NULL) return FALSE;

  for (p++; p != NULL && *p != '\0'; p++) {
    if (IS_LOWER (*p)) return FALSE;
  }
  return TRUE;
}

static void CleanEquals (
  CharPtr p
)

{
  if (p == NULL) return;

  for (; *p != '\0'; p++) {
    if (*p == '\"') {
      *p = '\'';
    }
  }
}

static CharPtr GetPubTitle (
  FmtType format,
  PubdescPtr pdp,
  CitSubPtr csp
)

{
  CitArtPtr        cap;
  CitBookPtr       cbp;
  CitGenPtr        cgp;
  Char             ch;
  CitPatPtr        cpp;
  MedlineEntryPtr  mep;
  CharPtr          ptr;
  CharPtr          title = NULL;
  ValNodePtr       ttl = NULL;
  ValNodePtr       vnp;

  if (csp != NULL) {
    if (format == GENBANK_FMT || format == GENPEPT_FMT) {
      title = "Direct Submission";
      return StringSave (title);
    } else if (format == EMBL_FMT || format == EMBLPEPT_FMT) {
      return NULL;
    }
  }
  if (pdp == NULL) return NULL;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    switch (vnp->choice) {
      case PUB_Gen :
        cgp = (CitGenPtr) vnp->data.ptrvalue;
        if (cgp != NULL) {
          if (! StringHasNoText (cgp->title)) return StringSave (cgp->title);
          if (! StringHasNoText (cgp->cit)) {
            ptr = StringStr (cgp->cit, "Title=\"");
            if (ptr != NULL) {
              title = StringSave (ptr + 7);
              for (ptr = title; *ptr != '\0'; ptr++) {
                if (*ptr == '"') {
                  *ptr = '\0';
                  break;
                }
              }
              return title;
            }
          }
        }
        break;
      case PUB_Sub :
        csp = (CitSubPtr) vnp->data.ptrvalue;
        if (csp != NULL) {
          if (format == GENBANK_FMT || format == GENPEPT_FMT) {
            title = "Direct Submission";
            return StringSave (title);
          } else if (format == EMBL_FMT || format == EMBLPEPT_FMT) {
            return NULL;
          }
        }
        break;
      case PUB_Medline :
        mep = (MedlineEntryPtr) vnp->data.ptrvalue;
        if (mep != NULL) {
          cap = mep->cit;
          if (cap != NULL) {
            ttl = cap->title;
          }
        }
        break;
      case PUB_Article :
        cap = (CitArtPtr) vnp->data.ptrvalue;
        if (cap != NULL) {
          ttl = cap->title;
        }
        break;
      /* case PUB_Book : */
      case PUB_Proc :
      case PUB_Man :
        cbp = (CitBookPtr) vnp->data.ptrvalue;
        if (cbp != NULL) {
          ttl = cbp->title;
          if (ttl != NULL) {
            title = (CharPtr) ttl->data.ptrvalue;
            if (! StringHasNoText (title)) {
              title = StringSave (title);
              if (StringLen (title) > 3) {
                ch = *title;
                if (IS_LOWER (ch)) {
                  *title = TO_UPPER (ch);
                }
                ptr = title;
                if (AllCaps (ptr)) {
                  for (ptr++; ptr != NULL && *ptr != '\0'; ptr++) {
                    ch = *ptr;
                    *ptr = TO_LOWER (ch);
                  }
                }
              }
              return title;
            }
          }
        }
        break;
      case PUB_Patent :
        cpp = (CitPatPtr) vnp->data.ptrvalue;
        if (cpp != NULL) {
          title = cpp->title;
          if (! StringHasNoText (title)) {
            return StringSave (title);
          }
        }
        break;
      default :
        break;
    }

    if (ttl != NULL) {
      title = (CharPtr) ttl->data.ptrvalue;
      if (! StringHasNoText (title)) {
        return StringSave (title);
      }
    }
  }

  return NULL;
}

static void CleanPubTitle (
  CharPtr title
)

{
  CharPtr  p;
  Boolean  remove_it;

  if (title == NULL) return;

  CleanEquals (title);

  for (p = title + StringLen (title) - 1; p > title + 2; p--) {
    if (*p == ' ') {
      *p = '\0';
    } else if (*p == '.') {
      remove_it = FALSE;
      if (p > title + 5) {
        if (*(p - 1) != '.' || *(p - 2) != '.') {
          remove_it = TRUE;
        }
      }
      if (remove_it) {
        *p = '\0';
      }
      break;
    } else {
      break;
    }
  }
}

/*
medline type page numbering is expanded (e.g., 125-35 -> 125-135,
F124-34 -> F124-F134, 12a-c -> 12a-12c).
If only one page is given, this is output without a dash.
Expanded numbering is validated to ensure that the
first number is smaller than or equal to the second and
that the first letter is less than or identical to the second
(i.e., a < c).  If the input is all letters (i.e., roman numerals)
this is not validated.

Return values:
 0 : valid page numbering.
-1 : invalid page numbering.
*/

#define MAX_PAGE_DIGITS 12

static Int2 FixPages (
  CharPtr out_pages,
  CharPtr in_pages
)

{
  Boolean dash=TRUE, first_alpha;
  Char firstbegin[MAX_PAGE_DIGITS];
  Char secondbegin[MAX_PAGE_DIGITS];
  Char firstend[MAX_PAGE_DIGITS];
  Char secondend[MAX_PAGE_DIGITS];
  Char temp[MAX_PAGE_DIGITS];
  CharPtr alphabegin, numbegin, alphaend, numend, ptr, in=in_pages;
  Int2 diff, index, retval=0;
  Int2 length_nb, length_ab, length_ne, length_ae;
  Int4 num1=0, num2=0;

  if (in_pages == NULL) return retval;

  while (*in != '\0')
  {      /* Check for digits in input*/
    if (IS_DIGIT(*in))
      break;
    in++;
  }

  if (*in == '\0' || (in != in_pages && *(in-1) == ' '))
  {    /* if all letters (i.e. roman numerals), put out. */
    out_pages = StringCpy(out_pages, in_pages);
    return retval;
  }

  in = in_pages;
  if (IS_DIGIT(*in))
  {      /* Do digits come first? */
    first_alpha = FALSE;
    index=0;
    while (IS_DIGIT(*in) || *in == ' ')
    {
      firstbegin[index] = *in;
      if (*in != ' ')
        index++;
      in++;
      if (*in == '-')
        break;

    }
    firstbegin[index] = '\0';
    index=0;
    if (*in != '-')
    {    /* After digits look for letters. */
      while (IS_ALPHA(*in)  || *in == ' ')
      {
        secondbegin[index] = *in;
        index++;
        in++;
        if (*in == '-')
          break;
      }
    }
    secondbegin[index] = '\0';
    if (*in == '-')    /* if dash is not present, note */
      in++;
    else
      dash=FALSE;
    index=0;
    while (IS_DIGIT(*in) || *in == ' ')
    {      /* Look for digits.  */
      firstend[index] = *in;
      if (*in != ' ')
        index++;
      in++;
    }
    firstend[index] = '\0';
    index=0;
    if (*in != '\0')
    {      /* Look for letters again. */
      while (IS_ALPHA(*in)  || *in == ' ')
      {
        secondend[index] = *in;
        index++;
        in++;
      }
    }
    secondend[index] = '\0';
  }
  else
  {      /* Do letters come first? */
    first_alpha = TRUE;
    index=0;
    while (IS_ALPHA(*in) || *in == ' ')
    {
      firstbegin[index] = *in;
      index++;
      in++;
      if (*in == '-')
        break;
    }
    firstbegin[index] = '\0';
    index=0;
    if (*in != '-')
    {    /* After letters look for digits.   */
      while (IS_DIGIT(*in)  || *in == ' ')
      {
        secondbegin[index] = *in;
        if (*in != ' ')
          index++;
        in++;
        if (*in == '-')
          break;
      }
    }
    secondbegin[index] = '\0';
    if (*in == '-')    /* Note if dash is missing. */
      in++;
    else
      dash=FALSE;
    index=0;
    while (IS_ALPHA(*in) || *in == ' ')
    {    /* Look for letters again. */
      firstend[index] = *in;
      index++;
      in++;
    }
    firstend[index] = '\0';
    index=0;
    if (*in != '\0')
    {    /* Any digits here? */
      while (IS_DIGIT(*in)  || *in == ' ')
      {
        secondend[index] = *in;
        if (*in != ' ')
          index++;
        in++;
      }
    }
    secondend[index] = '\0';
  }

  if (first_alpha)
  {
    alphabegin = firstbegin;
    numbegin = secondbegin;
    alphaend = firstend;
    numend = secondend;
  }
  else
  {
    numbegin = firstbegin;
    alphabegin = secondbegin;
    numend = firstend;
    alphaend = secondend;
  }

  length_nb = StringLen(numbegin);
  length_ab = StringLen(alphabegin);
  length_ne = StringLen(numend);
  length_ae = StringLen(alphaend);

  /* If no dash, but second letters or numbers present, reject. */
  if (dash == FALSE)
  {
    if (length_ne != 0 || length_ae != 0)
      retval = -1;
  }
  /* Check for situations like "AAA-123" or "222-ABC". */
  if (dash == TRUE)
  {
    if (length_ne == 0 && length_ab == 0)
      retval = -1;
    else if (length_ae == 0 && length_nb == 0)
      retval = -1;
  }

  /* The following expands "F502-512" into "F502-F512" and
  checks, for entries like "12a-12c" that a > c.  "12aa-12ab",
  "125G-137A", "125-G137" would be rejected. */
  if (retval == 0)
  {
    if (length_ab > 0)
    {
      if (length_ae > 0)
      {
        if (StringCmp(alphabegin, alphaend) != 0)
        {
          if (length_ab != 1 || length_ae != 1)
            retval = -1;
          else if (*alphabegin > *alphaend)
            retval = -1;
        }
      }
      else
      {
        alphaend = alphabegin;
        length_ae = length_ab;
      }
    }
    else if (length_ae > 0)
      retval = -1;
  }

/* The following expands "125-37" into "125-137".  */
  if (retval == 0)
  {
    if (length_nb > 0)
    {
      if (length_ne > 0)
      {
        diff = length_nb - length_ne;
        if (diff > 0)
        {
          index=0;
          while (numend[index] != '\0')
          {
            temp[index+diff] = numend[index];
            index++;
          }
          temp[index+diff] = numend[index];
          for (index=0; index<diff; index++)
            temp[index] = numbegin[index];
          index=0;
          while (temp[index] != '\0')
          {
            numend[index] = temp[index];
            index++;
          }
          numend[index] = temp[index];
        }
      }
      else
      {
        numend = numbegin;
        length_ne = length_nb;
      }

    }
    else if (length_ne > 0)
      retval = -1;
  /* Check that the first number is <= the second (expanded) number. */
    if (retval == 0)
    {
  /*    sscanf(numbegin, "%ld", &num_type);
      num1 = (Int4) num_type;
      sscanf(  numend, "%ld", &num_type);
      num2 = (Int4) num_type;
  */
      num1 = (Int4) atol(numbegin);
      num2 = (Int4) atol(numend);
      if (num2 < num1)
        retval = -1;
    }
  }

  if (retval == -1)
  {
    out_pages = StringCpy(out_pages, in_pages);
  }
  else
  {
    ptr = out_pages;
  /* Place expanded and validated page numbers into "out_pages". */
    if (first_alpha)
    {
      while (*alphabegin != '\0')
      {
        *ptr = *alphabegin;
        alphabegin++;
        ptr++;
      }
      while (*numbegin != '\0')
      {
        *ptr = *numbegin;
        numbegin++;
        ptr++;
      }
      if (dash == TRUE)
      {
        *ptr = '-';
        ptr++;
        while (*alphaend != '\0')
        {
          *ptr = *alphaend;
          alphaend++;
          ptr++;
        }
        while (*numend != '\0')
        {
          *ptr = *numend;
          numend++;
          ptr++;
        }
      }
      *ptr = '\0';
    }
    else
    {
      while (*numbegin != '\0')
      {
        *ptr = *numbegin;
        numbegin++;
        ptr++;
      }
      while (*alphabegin != '\0')
      {
        *ptr = *alphabegin;
        alphabegin++;
        ptr++;
      }
      if (dash == TRUE)
      {
        *ptr = '-';
        ptr++;
        while (*numend != '\0')
        {
          *ptr = *numend;
          numend++;
          ptr++;
        }
        while (*alphaend != '\0')
        {
          *ptr = *alphaend;
          alphaend++;
          ptr++;
        }
      }
      *ptr = '\0';
    }
  }
  return retval;
}

/* !!! still need to add StripParanthesis equivalent !!! */

static void DoSup (
  ValNodePtr PNTR head,
  CharPtr issue,
  CharPtr part_sup,
  CharPtr part_supi
)

{
  size_t   len;
  CharPtr  str;
  CharPtr  temp;

  len = StringLen (issue) + StringLen (part_sup) + StringLen (part_supi) + 25;
  str = MemNew (sizeof (Char) * len);
  if (str == NULL) return;
  temp = str;

  if (! StringHasNoText (part_sup)) {
    *temp = ' ';
    temp++;
    temp = StringMove (temp, part_sup);
  }
  if (StringHasNoText (issue) && StringHasNoText (part_supi)) {
    ValNodeCopyStr (head, 0, str);
    MemFree (str);
    return;
  }
  *temp = ' ';
  temp++;
  *temp = '(';
  temp++;
  if (! StringHasNoText (issue)) {
    temp = StringMove (temp, issue);
  }
  if (! StringHasNoText (part_supi)) {
    *temp = ' ';
    temp++;
    temp = StringMove (temp, part_supi);
  }
  *temp = ')';
  temp++;
  ValNodeCopyStr (head, 0, str);
  MemFree (str);
}

static CharPtr FormatCitJour (
  FmtType format,
  Boolean citArtIsoJta,
  CitJourPtr cjp
)

{
  Char        buf [256];
  DatePtr     dp;
  ValNodePtr  head = NULL;
  ImprintPtr  imp;
  CharPtr     issue = NULL;
  Char        pages [128];
  CharPtr     part_sup = NULL;
  CharPtr     part_supi = NULL;
  CharPtr     rsult = NULL;
  CharPtr     title;
  ValNodePtr  ttl;
  CharPtr     volume;
  Char        year [8];

  if (cjp == NULL) return NULL;

  ttl = cjp->title;
  if (ttl == NULL) return NULL;

  /* always use iso_jta title if present */

  while (ttl != NULL && ttl->choice != Cit_title_iso_jta) {
    ttl = ttl->next;
  }

  /* release mode requires iso_jta title */

  if (ttl == NULL) {
    if (citArtIsoJta) return NULL;
    ttl = cjp->title;
  }

  imp = cjp->imp;
  if (imp == NULL) return NULL;

  dp = imp->date;
  year [0] = '\0';
  if (dp != NULL) {
    if (dp->data [0] == 1) {
      if (dp->data [1] != 0) {
        sprintf (year, " (%ld)", (long) (1900 + dp->data [1]));
      }
    } else {
      StringCpy (year, " (");
      StringNCat (year, dp->str, 4);
      StringCat (year, ")");
    }
  }

  if (imp->prepub == 1 || imp->prepub == 255) {
    sprintf (buf, "Unpublished %s", year);
    return StringSave (buf);
  }

  title = (CharPtr) ttl->data.ptrvalue;
  if (StringLen (title) < 3) return StringSave (".");

  ValNodeCopyStr (&head, 0, title);

  volume = imp->volume;
  if (format == GENBANK_FMT || format == GENPEPT_FMT) {
    issue = imp->issue;
    part_sup = imp->part_sup;
    part_supi = imp->part_supi;
  }
  pages [0] = '\0';
  FixPages (pages, imp->pages);

  if (! StringHasNoText (volume)) {
    AddValNodeString (&head, " ", volume, NULL);
  }

  if ((! StringHasNoText (volume)) || (! StringHasNoText (pages))) {
    DoSup (&head, issue, part_sup, part_supi);
  }

  if (format == GENBANK_FMT || format == GENPEPT_FMT) {
    if (! StringHasNoText (pages)) {
      AddValNodeString (&head, ", ", pages, NULL);
    }
  } else if (format == EMBL_FMT || format == EMBLPEPT_FMT) {
    if (! StringHasNoText (pages)) {
      AddValNodeString (&head, ":", pages, NULL);
    } else if (imp->prepub == 2 || (StringHasNoText (volume))) {
      ValNodeCopyStr (&head, 0, " 0:0-0");
    }
  }

  ValNodeCopyStr (&head, 0, year);

  if (format == GENBANK_FMT || format == GENPEPT_FMT) {
    if (imp->prepub == 2) {
      ValNodeCopyStr (&head, 0, " In press");
    }
  }

  rsult = MergeValNodeStrings (head);
  ValNodeFreeData (head);

  return rsult;
}

static CharPtr MakeAffilStr (
  AffilPtr afp
)

{
  ValNodePtr  head = NULL;
  CharPtr     prefix = "";
  CharPtr     rsult = NULL;

  if (afp == NULL) return NULL;

  if (! StringHasNoText (afp->affil)) {
    ValNodeCopyStr (&head, 0, afp->affil);
    prefix = ", ";
  }

  if (afp->choice == 2) {
    if (! StringHasNoText (afp->div)) {
      AddValNodeString (&head, prefix, afp->div, NULL);
      prefix = ", ";
    }
    if (! StringHasNoText (afp->street)) {
      AddValNodeString (&head, prefix, afp->street, NULL);
      prefix = ", ";
    }
    if (! StringHasNoText (afp->city)) {
      AddValNodeString (&head, prefix, afp->city, NULL);
      prefix = ", ";
    }
    if (! StringHasNoText (afp->sub)) {
      AddValNodeString (&head, prefix, afp->sub, NULL);
      prefix = ", ";
    }
    if (! StringHasNoText (afp->country)) {
      AddValNodeString (&head, prefix, afp->country, NULL);
      prefix = ", ";
    }
  }

  rsult = MergeValNodeStrings (head);
  ValNodeFreeData (head);

  return rsult;
}

static CharPtr GetAffil (
  AffilPtr afp
)

{
  Boolean need_comma=FALSE;
  CharPtr string=NULL, temp, ptr;
  Char ch;
  Int2 aflen=15;

  if (afp == NULL) return NULL;
  if (afp) {
    if (afp -> choice == 1){
      if (afp -> affil){
        aflen += StringLen(afp -> affil);
      }
    }else if (afp -> choice == 2){
      aflen += StringLen (afp -> affil) +
      StringLen (afp -> div) +
      StringLen (afp -> city) +
      StringLen (afp -> sub) +
      StringLen (afp -> street) +
      StringLen (afp -> country) + StringLen(afp->postal_code);
    }

    temp = string = MemNew(aflen);

    if ( afp -> choice == 1){
       if (afp -> affil){
        ptr = afp->affil;
        while ((*temp = *ptr) != '\0')
        {
          temp++; ptr++;
        }
       }
    }else if (afp -> choice == 2){

      if( afp -> div) {
        if (need_comma)
        {
          *temp = ','; temp++;
          *temp = ' '; temp++;
        }
        ptr = afp->div;
        while ((*temp = *ptr) != '\0')
        {
          temp++; ptr++;
        }
        need_comma = TRUE;
      }

      if(afp -> affil) {
        if (need_comma)
        {
          *temp = ','; temp++;
          *temp = ' '; temp++;
        }
        ptr = afp->affil;
        while ((*temp = *ptr) != '\0')
        {
          temp++; ptr++;
        }
        need_comma = TRUE;
      }

      if(afp -> street) {
        if (need_comma)
        {
          *temp = ','; temp++;
          *temp = ' '; temp++;
        }
        ptr = afp->street;
        while ((*temp = *ptr) != '\0')
        {
          temp++; ptr++;
        }
        need_comma = TRUE;
      }

      if( afp -> city) {
        if (need_comma)
        {
          *temp = ','; temp++;
          *temp = ' '; temp++;
        }
        ptr = afp->city;
        while ((*temp = *ptr) != '\0')
        {
          temp++; ptr++;
        }
        need_comma = TRUE;
      }

      if( afp -> sub) {
        if (need_comma)
        {
          *temp = ','; temp++;
          *temp = ' '; temp++;
        }
        ptr = afp->sub;
        while ((*temp = *ptr) != '\0')
        {
          temp++; ptr++;
        }
        need_comma = TRUE;
      }

      if( afp -> postal_code){
        *temp = ' ';
        temp++;
        ptr = afp->postal_code;
        while ((*temp = *ptr) != '\0')
        {
          temp++; ptr++;
        }
      }

      if( afp -> country){
        if (need_comma)
        {
          *temp = ','; temp++;
          *temp = ' '; temp++;
        }
        ptr = afp->country;
        while ((*temp = *ptr) != '\0')
        {
          temp++; ptr++;
        }
        need_comma = TRUE;
      }
    }
    temp++;
    *temp = '\0';
  }

    /* convert double quotes to single quotes */

    ptr = string;
    ch = *ptr;
    while (ch != '\0') {
      if (ch == '\"') {
        *ptr = '\'';
      }
      ptr++;
      ch = *ptr;
    }

  return string;
}

static CharPtr FormatCitBookArt (
  FmtType format,
  CitBookPtr cbp
)

{
  AffilPtr     afp;
  AuthListPtr  alp;
  CharPtr      book_title = NULL;
  Char         buf [256];
  Char         ch;
  DatePtr      dp;
  ValNodePtr   head = NULL;
  ImprintPtr   imp;
  CharPtr      issue = NULL;
  ValNodePtr   names = NULL;
  Char         pages [128];
  CharPtr      part_sup = NULL;
  CharPtr      part_supi = NULL;
  CharPtr      rsult = NULL;
  CharPtr      str;
  CharPtr      title;
  ValNodePtr   ttl;
  ValNodePtr   vnp;
  CharPtr      volume;
  Char         year [8];

  if (cbp == NULL) return NULL;

  ttl = cbp->title;
  if (ttl == NULL) return NULL;

  imp = cbp->imp;
  if (imp == NULL) return NULL;

  dp = imp->date;
  year [0] = '\0';
  if (dp != NULL) {
    if (dp->data [0] == 1) {
      if (dp->data [1] != 0) {
        sprintf (year, "(%ld)", (long) (1900 + dp->data [1]));
      }
    } else {
      StringCpy (year, "(");
      StringNCat (year, dp->str, 4);
      StringCpy (year, ")");
    }
  }

  if (imp->prepub == 1 || imp->prepub == 255) {
    sprintf (buf, "Unpublished %s", year);
    return StringSave (buf);
  }

  title = (CharPtr) ttl->data.ptrvalue;
  if (StringLen (title) < 3) return StringSave (".");

  ValNodeCopyStr (&head, 0, "(in) ");

  alp = cbp->authors;
  if (alp != NULL) {
    str = GetAuthorsString (format, alp, NULL, NULL, NULL);
    if (str != NULL) {
      ValNodeCopyStr (&head, 0, str);
      names = alp->names;
      if (names != NULL) {
        if (names->next != NULL) {
          ValNodeCopyStr (&head, 0, " (Eds.);");
        } else {
          ValNodeCopyStr (&head, 0, " (Ed.);");
        }
      }
      ValNodeCopyStr (&head, 0, "\n");
    }
    MemFree (str);
  }

  book_title = StringSaveNoNull (title);
  vnp = ValNodeAddStr (&head, 0, book_title);
  if (book_title != NULL) {

    /* make book title all caps */

    title = book_title;
    ch = *title;
    while (ch != '\0') {
      *title = TO_UPPER (ch);
      title++;
      ch = *title;
    }
  }

  volume = imp->volume;
  if (format == GENBANK_FMT || format == GENPEPT_FMT) {
    issue = imp->issue;
    part_sup = imp->part_sup;
    part_supi = imp->part_supi;
  }
  pages [0] = '\0';
  FixPages (pages, imp->pages);

  if ((! StringHasNoText (volume)) && (StringCmp (volume, "0") != 0)) {
    AddValNodeString (&head, ", Vol. ", volume, NULL);
    DoSup (&head, issue, part_sup, part_supi);
  }

  if (! StringHasNoText (pages)) {
    AddValNodeString (&head, ": ", pages, NULL);
  }

  if (book_title != NULL) {
    ValNodeCopyStr (&head, 0, ";\n");
  }

  afp = imp->pub;
  if (afp != NULL) {
    str = MakeAffilStr (afp);
    if (str != NULL) {
      ValNodeCopyStr (&head, 0, str);
      ValNodeCopyStr (&head, 0, " ");
      MemFree (str);
    }
  }

  AddValNodeString (&head, NULL, year, NULL);

  if (format == GENBANK_FMT || format == GENPEPT_FMT) {
    if (imp->prepub == 2) {
      ValNodeCopyStr (&head, 0, " In press");
    }
  }

  rsult = MergeValNodeStrings (head);
  ValNodeFreeData (head);

  return rsult;
}

static CharPtr FormatCitBook (
  FmtType format,
  CitBookPtr cbp
)

{
  AffilPtr   afp;
  char       year[5];
  CharPtr    bookTitle=NULL;
  CharPtr    retval = NULL;
  CharPtr    temp;
  DatePtr    dp;
  ImprintPtr ip;
  int        aflen = 0;
  CharPtr    p;
  CharPtr    affilStr = NULL;

  /* Check parameters */

  if (cbp == NULL)
    return NULL;

  if ( cbp -> othertype != 0)
    return NULL;

  ip = cbp -> imp;

  /* Format the year */

  dp = ip -> date;
  year[0] = '\0';

  if ( dp -> data[0] == 1)
    sprintf(year,"%ld",(long) ( 1900+dp -> data[1]));
  else
    {
      StringNCpy( (CharPtr) year, (CharPtr) dp -> str, (size_t) 4);
      year[4] = '\0';
    }

  /* Get the book title */

  if (cbp->title)
    bookTitle = StringSave(cbp -> title -> data.ptrvalue);

  /* Get the affiliation length */

  if ( ip -> pub){
    afp = ip -> pub;
    aflen = StringLen(afp -> affil)+ 5;
    if ( afp -> choice == 2){
      aflen += 3 + StringLen(afp -> div);
      aflen += 3 + StringLen(afp -> street);
      aflen += 3 + StringLen(afp -> city);
      aflen += 3 + StringLen(afp -> sub);
      aflen += 3 + StringLen(afp -> country);
    }
  } else{
    aflen = 22;
  }
  if (ip->prepub == 2)
    aflen += 10;

  /* Create a Char String big enough to hold */
  /* the title, year, and affiliation.       */

  temp = retval = MemNew( (size_t) (30+StringLen( bookTitle)+StringLen( year) + aflen) );

  /* Convert the title to upper case and */
  /* add it to the string.               */

  for ( p = bookTitle; *p; p++)
    *p = TO_UPPER(*p);

  /* temp = StringMove(temp, "Book: "); */
  temp = StringMove(temp, "(in) ");
  temp = StringMove(temp, bookTitle);
  temp = StringMove(temp, ".");

  /* Add the affiliation to the string */

  if ( ip -> pub)
    {
      afp = ip -> pub;
      *temp = ' ';
      temp++;
      affilStr = MakeAffilStr(afp);
      temp = StringMove(temp,affilStr);
    }

  /* Add the year to the string */

  if (year[0] != '\0')
    {
      if (affilStr != NULL)
        temp = StringMove(temp," (");
      else
        temp = StringMove(temp, "(");
      temp = StringMove(temp, year);
      temp = StringMove(temp, ")");
    }

  /* If in press, add note */

  if (ip->prepub == 2)
    temp = StringMove(temp, ", In press");

  /* Clean up and return */

  if (bookTitle)
    MemFree(bookTitle);

  return retval;

}

static CharPtr FormatThesis (
  FmtType format,
  CitBookPtr cbp
)

{
  AffilPtr     afp;
  Char         ch;
  DatePtr      dp;
  ValNodePtr   head = NULL;
  ImprintPtr   imp;
  CharPtr      ptr;
  CharPtr      rsult = NULL;
  CharPtr      str;
  CharPtr      suffix = NULL;
  Char         year [8];

  if (cbp == NULL) return NULL;
  if (cbp->othertype != 2 || cbp->let_type != 3) return NULL;

  imp = cbp->imp;
  if (imp == NULL) return NULL;

  dp = imp->date;
  year [0] = '\0';
  if (dp != NULL) {
    if (dp->data [0] == 1) {
      if (dp->data [1] != 0) {
        sprintf (year, "%ld", (long) (1900 + dp->data [1]));
      }
    } else {
      StringNCpy (year, dp->str, (size_t) 4);
      year [4] = '\0';
    }
  }

  AddValNodeString (&head, "Thesis (", year, ")");

  if (imp->prepub == 2) {
    suffix = ", In press";
  }

  str = NULL;
  afp = imp->pub;
  if (afp != NULL) {
    if (afp->choice == 1) {
      if (StringLen (afp->affil) > 7) {
        str = StringSave (afp->affil);
      }
    } else if (afp->choice == 2) {
      str = MakeAffilStr (afp);
    }
  }

  if (str != NULL) {

    /* convert double quotes to single quotes */

    ptr = str;
    ch = *ptr;
    while (ch != '\0') {
      if (ch == '\"') {
        *ptr = '\'';
      }
      ptr++;
      ch = *ptr;
    }
    AddValNodeString (&head, " ", str, suffix);
    MemFree (str);
  }

  rsult = MergeValNodeStrings (head);
  ValNodeFreeData (head);

  return rsult;
}

static CharPtr FormatCitArt (
  FmtType format,
  Boolean citArtIsoJta,
  CitArtPtr cap
)

{
  CitBookPtr  cbp;
  CitJourPtr  cjp;
  CharPtr     rsult = NULL;

  if (cap == NULL) return NULL;

  switch (cap->from) {
    case 1 :
      cjp = (CitJourPtr) cap->fromptr;
      if (cjp != NULL) {
        rsult = FormatCitJour (format, citArtIsoJta, cjp);
      }
      break;
    case 2 :
      cbp = (CitBookPtr) cap->fromptr;
      if (cbp != NULL) {
        rsult = FormatCitBookArt (format, cbp);
      }
      break;
    case 3 :
      cbp = (CitBookPtr) cap->fromptr;
      if (cbp != NULL) {
        rsult = FormatCitBookArt (format, cbp);
      }
      break;
    default :
      break;
  }

  return rsult;
}

static CharPtr FormatCitPat (
  FmtType   format,
  CitPatPtr cpp,
  SeqIdPtr  seqidp
)

{
  AffilPtr       afp;
  AuthListPtr    alp;
  Char           date [40];
  ValNodePtr     head = NULL;
  CharPtr        prefix = NULL;
  CharPtr        rsult = NULL;
  SeqIdPtr       sip;
  CharPtr        suffix = NULL;
  PatentSeqIdPtr psip;
  Int4           pat_seqid = 0;
  Char           buf[10];

  if (cpp == NULL) return NULL;

  if (format == GENBANK_FMT || format == GENPEPT_FMT) {
    ValNodeCopyStr (&head, 0, "Patent: ");
    suffix = " ";
  } else if (format == EMBL_FMT || format == EMBLPEPT_FMT) {
    ValNodeCopyStr (&head, 0, "Patent number ");
  }

  if (! StringHasNoText (cpp->country)) {
    AddValNodeString (&head, NULL, cpp->country, suffix);
  }

  if (! StringHasNoText (cpp->number)) {
    ValNodeCopyStr (&head, 0, cpp->number);
  } else if (! StringHasNoText (cpp->app_number)) {
    AddValNodeString (&head, "(", cpp->app_number, ")");
  }

  if (! StringHasNoText (cpp->doc_type)) {
    AddValNodeString (&head, "-", cpp->doc_type, NULL);
  }

  /* pat_seqid test */

  for (sip = seqidp; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_PATENT) {
      psip = (PatentSeqIdPtr) sip -> data.ptrvalue;
      if (psip != NULL) {
        pat_seqid = psip->seqid;
      }
    }
  }
  if (pat_seqid > 0) {
    if (format == EMBL_FMT) {
      sprintf(buf,"%s%ld%s", "/", (long) pat_seqid, ", ");
      ValNodeCopyStr (&head, 0, buf);
    } else {
      sprintf(buf,"%s%ld ", " ", (long) pat_seqid);
      ValNodeCopyStr (&head, 0, buf);
    }
  } else {
    ValNodeCopyStr (&head, 0, " ");
  }

  /* Date */

  date [0] = '\0';
  if (cpp->date_issue != NULL) {
    DateToGB (date, cpp->date_issue, FALSE);
  } else if (cpp->app_date != NULL) {
    DateToGB (date, cpp->app_date, FALSE);
  }
  if (! StringHasNoText (date)) {
    ValNodeCopyStr (&head, 0, date);
  }

  if (format == GENBANK_FMT || format == GENPEPT_FMT) {
    ValNodeCopyStr (&head, 0, ";");
  } else if (format == EMBL_FMT || format == EMBLPEPT_FMT) {
    ValNodeCopyStr (&head, 0, ".");
  }

  alp = cpp->authors;
  if (alp != NULL) {
    afp = alp->affil;
    if (afp != NULL) {
      suffix = NULL;
      if (afp->choice == 2) {
        suffix = ";";
      }

      /* If and of the affiliation fields are */
      /* non-blank, put them on a new line.   */

      if ((! StringHasNoText (afp->affil)) ||
          (! StringHasNoText (afp->street)) ||
          (! StringHasNoText (afp->div)) ||
          (! StringHasNoText (afp->city)) ||
          (! StringHasNoText (afp->sub)) ||
          (! StringHasNoText (afp->country)))
        ValNodeCopyStr (&head, 0, "\n");

      /* Write out the affiliation fields */

      if (! StringHasNoText (afp->affil)) {
        AddValNodeString (&head, NULL, afp->affil, suffix);
        prefix = " ";
      }
      if (! StringHasNoText (afp->street)) {
        AddValNodeString (&head, prefix, afp->street, ";");
        prefix = " ";
      }
      if (! StringHasNoText (afp->div)) {
        AddValNodeString (&head, prefix, afp->div, ";");
        prefix = " ";
      }
      if (! StringHasNoText (afp->city)) {
        AddValNodeString (&head, prefix, afp->city, NULL);
        prefix = ", ";
      }
      if (! StringHasNoText (afp->sub)) {
        AddValNodeString (&head, prefix, afp->sub, NULL);
      }
      if (! StringHasNoText (afp->country)) {
        AddValNodeString (&head, ";\n", afp->country, ";");
      }
    }
  }

  rsult = MergeValNodeStrings (head);
  ValNodeFreeData (head);

  /*
  s_StringCleanup(rsult);
  */

  return rsult;
}

static CharPtr FormatCitGen (
  FmtType format,
  Boolean dropBadCitGens,
  Boolean noAffilOnUnpub,
  CitGenPtr cgp
)

{
  CharPtr      affil = NULL;
  AuthListPtr  alp = NULL;
  Char         ch;
  DatePtr      dp;
  ValNodePtr   head = NULL;
  CharPtr      inpress = NULL;
  CharPtr      journal = NULL;
  Char         pages [128];
  CharPtr      prefix = NULL;
  CharPtr      ptr;
  CharPtr      rsult = NULL;
  Char         year [8];

  if (cgp == NULL) return NULL;

  if (cgp->journal == NULL && StringNICmp (cgp->cit, "unpublished", 11) == 0) {
    if (noAffilOnUnpub) {

      /* !!! temporarily put date in unpublished citation for QA !!! */

      if (dropBadCitGens) {
        year [0] = '\0';
        dp = cgp->date;
        if (dp != NULL) {
          if (dp->data [0] == 1) {
            if (dp->data [1] != 0) {
              sprintf (year, " (%ld)", (long) (1900 + dp->data [1]));
            }
          } else {
            StringCpy (year, " (");
            StringNCat (year, dp->str, 4);
            StringCat (year, ")");
          }
        }
        AddValNodeString (&head, NULL, "Unpublished", NULL);
        AddValNodeString (&head, NULL, year, NULL);
        rsult = MergeValNodeStrings (head);
        ValNodeFreeData (head);
        return rsult;
      }

      /* !!! remove above section once QA against asn2ff is done !!! */

      return StringSave ("Unpublished");
    }

    alp = cgp->authors;
    if (alp != NULL) {
      affil = GetAffil (alp->affil);
      if (! StringHasNoText (affil)) {
        rsult = MemNew ((size_t) StringLen (affil) + (size_t) StringLen (cgp->cit) + 15);
        StringCpy (rsult, "Unpublished ");
        StringCat (rsult, affil);
        TrimSpacesAroundString (rsult);
        return rsult;
      }
    }

    rsult = StringSave (cgp->cit);
    TrimSpacesAroundString (rsult);
    return rsult;
  }

  year [0] = '\0';
  dp = cgp->date;
  if (dp != NULL) {
    if (dp->data [0] == 1) {
      if (dp->data [1] != 0) {
        sprintf (year, " (%ld)", (long) (1900 + dp->data [1]));
      }
    } else {
      StringCpy (year, " (");
      StringNCat (year, dp->str, 4);
      StringCat (year, ")");
    }
  }

  pages [0] = '\0';
  if (cgp->pages != NULL) {
    FixPages (pages, cgp->pages);
  }

  if (cgp->journal != NULL) {
    journal = (CharPtr) cgp->journal->data.ptrvalue;
  }
  if (cgp->cit != NULL) {
    ptr = StringStr (cgp->cit, "Journal=\"");
    if (ptr != NULL) {
      journal = ptr + 9;
    } else if (StringNICmp (cgp->cit, "submitted", 8) == 0 ||
               StringNICmp (cgp->cit, "unpublished", 11) == 0) {

      if ((! dropBadCitGens) || journal != NULL) {
        inpress = cgp->cit;
      } else {
        inpress = "Unpublished";
      }
    } else if (StringNICmp (cgp->cit, "Online Publication", 18) == 0 ||
               StringNICmp (cgp->cit, "Published Only in DataBase", 26) == 0) {

      inpress = cgp->cit;
    } else if ((! dropBadCitGens) && journal == NULL) {
      journal = cgp->cit;
    }
  }
  if (journal != NULL) {
    journal = StringSave (journal);
    for (ptr = journal, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
      if (ch == '=' || ch == '\"') {
        *ptr = '\0';
      }
    }
    ValNodeAddStr (&head, 0, journal);
    prefix = " ";
  }

  if (! StringHasNoText (inpress)) {
    AddValNodeString (&head, prefix, inpress, NULL);
    prefix = " ";
  }

  if (! StringHasNoText (cgp->volume)) {
    AddValNodeString (&head, prefix, cgp->volume, NULL);
  }

  if (! StringHasNoText (pages)) {
    if (format == GENBANK_FMT) {
      AddValNodeString (&head, ", ", pages, NULL);
    } else if (format == EMBL_FMT) {
      AddValNodeString (&head, ":", pages, NULL);
    }
  }

  if (! StringHasNoText (year)) {
    AddValNodeString (&head, NULL, year, NULL);
  }

  rsult = MergeValNodeStrings (head);
  ValNodeFreeData (head);

  return rsult;
}

static CharPtr FormatCitSub (
  FmtType format,
  CitSubPtr csp
)

{
  CharPtr      affil;
  AffilPtr     afp;
  AuthListPtr  alp;
  Char         buf [256];
  Char         date [40];
  ValNodePtr   head = NULL;
  CharPtr      rsult = NULL;

  if (csp == NULL) return NULL;

  date [0] = '\0';
  if (csp->date != NULL) {
    DateToGB (date, csp->date, TRUE);
  }
  if (StringHasNoText (date)) {
    StringCpy (date, "??-???-????");
  }

  sprintf (buf, "Submitted (%s)", date);
  ValNodeCopyStr (&head, 0, buf);

  alp = csp->authors;
  if (alp != NULL) {
    afp = alp->affil;
    if (afp != NULL) {
      affil = GetAffil (afp);
      if (format == EMBL_FMT || format == EMBLPEPT_FMT) {
        if (StringNCmp(affil, " to the EMBL/GenBank/DDBJ databases.", 36) != 0) {
          ValNodeCopyStr (&head, 0, " to the EMBL/GenBank/DDBJ databases\n");
        } else {
          ValNodeCopyStr (&head, 0, " ");
        }
      } else {
        ValNodeCopyStr (&head, 0, " ");
      }
      ValNodeCopyStr (&head, 0, affil);
      MemFree (affil);
    } else if (format == EMBL_FMT || format == EMBLPEPT_FMT) {
      ValNodeCopyStr (&head, 0, " to the EMBL/GenBank/DDBJ databases\n");
    }
  }

  rsult = MergeValNodeStrings (head);
  ValNodeFreeData (head);

  return rsult;
}

static CharPtr GetPubJournal (
  FmtType format,
  Boolean dropBadCitGens,
  Boolean noAffilOnUnpub,
  Boolean citArtIsoJta,
  PubdescPtr pdp,
  CitSubPtr csp,
  SeqIdPtr  seqidp,
  IndxPtr index
)

{
  CitArtPtr        cap;
  CitBookPtr       cbp;
  CitGenPtr        cgp;
  CitPatPtr        cpp;
  CharPtr          journal = NULL;
  MedlineEntryPtr  mep;
  ValNodePtr       vnp;

  if (csp != NULL) {
    return FormatCitSub (format, csp);
  }
  if (pdp == NULL) return NULL;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    switch (vnp->choice) {
      case PUB_Gen :
        cgp = (CitGenPtr) vnp->data.ptrvalue;
        if (cgp != NULL) {
          if (StringNICmp ("BackBone id_pub", cgp->cit, 15) != 0) {
            if (cgp->cit == NULL && cgp->journal == NULL && cgp->date == NULL && cgp->serial_number) {
              break; /* skip just serial number */
            }
          }
          journal = FormatCitGen (format, dropBadCitGens, noAffilOnUnpub, cgp);
        }
        break;
      case PUB_Sub :
        csp = (CitSubPtr) vnp->data.ptrvalue;
        if (csp != NULL) {
          journal = FormatCitSub (format, csp);
        }
        break;
      case PUB_Medline :
        mep = (MedlineEntryPtr) vnp->data.ptrvalue;
        if (mep != NULL) {
          cap = mep->cit;
          if (cap != NULL) {
            journal = FormatCitArt (format, citArtIsoJta, cap);
          }
        }
        break;
      case PUB_Article :
        cap = (CitArtPtr) vnp->data.ptrvalue;
        if (cap != NULL) {
          journal = FormatCitArt (format, citArtIsoJta, cap);
        }
        break;
      case PUB_Book :
      case PUB_Proc :
        cbp = (CitBookPtr) vnp->data.ptrvalue;
        if (cbp != NULL) {
          journal = FormatCitBook (format, cbp);
        }
        break;
      case PUB_Man :
        cbp = (CitBookPtr) vnp->data.ptrvalue;
        if (cbp != NULL) {
          journal = FormatThesis (format, cbp);
        }
        break;
      case PUB_Patent :
        cpp = (CitPatPtr) vnp->data.ptrvalue;
        if (cpp != NULL) {
          journal = FormatCitPat (format, cpp, seqidp);
        }
        break;
      default :
        break;
    }

    /* optionally populate indexes for NCBI internal database */

    if (index != NULL && journal != NULL) {

      /* skip non-informative cit-gens */

      if (StringNICmp (journal, "submitted", 8) == 0 ||
          StringNICmp (journal, "unpublished", 11) == 0 ||
          StringNICmp (journal, "Online Publication", 18) == 0 ||
          StringNICmp (journal, "Published Only in DataBase", 26) == 0) {
      } else {
        ValNodeCopyStrToHead (&(index->journals), 0, journal);
      }
    }

    if (journal != NULL) return journal;
  }

  return NULL;
}

static Int4 GetMuid (
  PubdescPtr pdp
)

{
  MedlineEntryPtr  mep;
  ValNodePtr       vnp;

  if (pdp == NULL) return 0;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    switch (vnp->choice) {
      case PUB_Medline :
        mep = (MedlineEntryPtr) vnp->data.ptrvalue;
        if (mep != NULL) {
          return mep->uid;
        }
        break;
      case PUB_Muid :
        return vnp->data.intvalue;
      default :
        break;
    }
  }

  return 0;
}

static Int4 GetPmid (
  PubdescPtr pdp
)

{
  MedlineEntryPtr  mep;
  ValNodePtr       vnp;

  if (pdp == NULL) return 0;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    switch (vnp->choice) {
      case PUB_Medline :
        mep = (MedlineEntryPtr) vnp->data.ptrvalue;
        if (mep != NULL) {
          return mep->pmid;
        }
        break;
      case PUB_PMid :
        return vnp->data.intvalue;
      default :
        break;
    }
  }

  return 0;
}

static CharPtr remarksText [] = {
  "full automatic", "full staff_review", "full staff_entry",
  "simple staff_review", "simple staff_entry", "simple automatic",
  "unannotated automatic", "unannotated staff_review", "unannotated staff_entry",
  NULL
};

static void AddReferenceToGbseq (
  GBSeqPtr gbseq,
  GBReferencePtr gbref,
  CharPtr str
)

{
  CharPtr  copy;
  CharPtr  ptr;
  CharPtr  ref;

  if (gbseq == NULL || gbref == NULL || StringHasNoText (str)) return;

  copy = StringSave (str);

  /* link in reverse order, to be reversed in slash block */

  gbref->next = gbseq->references;
  gbseq->references = gbref;

  /* now parse or make ASN required default values for remaining fields */

  if (StringNCmp (copy, "REFERENCE   ", 12) == 0) {
    ref = copy + 12;
    ptr = StringStr (ref, "\n  AUTHORS");
    if (ptr == NULL) {
      ptr = StringStr (ref, ")\n");
    }
    if (ptr != NULL) {
      *ptr = '\0';
      gbref->reference = StringSave (ref);
    }
  }

  if (gbref->reference == NULL) {
    gbref->reference = StringSave ("?");
  }

  if (gbref->journal == NULL) {
    gbref->journal = StringSave ("?");
  }

  MemFree (copy);
}

static CharPtr FormatReferenceBlock (
  Asn2gbFormatPtr afp,
  BaseBlockPtr bbp
)

{
  IntAsn2gbJobPtr    ajp;
  AuthListPtr        alp;
  Asn2gbSectPtr      asp;
  BioseqPtr          bsp;
  Char               buf [150];
  CitArtPtr          cap;
  Char               ch;
  CitJourPtr         cjp;
  Boolean            citArtIsoJta;
  CharPtr            consortium;
  CitRetractPtr      crp;
  CitSubPtr          csp = NULL;
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext;
  Int4               gibbsq;
  GBReferencePtr     gbref = NULL;
  GBSeqPtr           gbseq;
  Int2               i;
  ImprintPtr         imp;
  IndxPtr            index;
  IntRefBlockPtr     irp;
  size_t             len;
  SeqLocPtr          loc = NULL;
  Int4               muid = 0;
  Boolean            needsPeriod = FALSE;
  SeqLocPtr          nextslp;
  Boolean            notFound;
  ObjMgrDataPtr      omdp;
  PubdescPtr         pdp = NULL;
  Int4               pmid = 0;
  CharPtr            prefix = NULL;
  RefBlockPtr        rbp;
  SubmitBlockPtr     sbp;
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp = NULL;
  SeqIdPtr           sip;
  SeqLocPtr          slp;
  SeqSubmitPtr       ssp;
  Int4               start;
  Int4               stop;
  CharPtr            str;
  Boolean            strict_isojta;
  CharPtr            suffix = NULL;
  CharPtr            tmp;
  Boolean            trailingPeriod = TRUE;
  ValNodePtr         vnp;
  StringItemPtr      ffstring, temp;

  if (afp == NULL || bbp == NULL) return NULL;
  rbp = (RefBlockPtr) bbp;
  ajp = afp->ajp;
  if (ajp == NULL) return NULL;
  asp = afp->asp;
  if (asp == NULL) return NULL;
  bsp = asp->bsp;
  if (bsp == NULL) return NULL;

  ffstring = FFGetString(ajp);
  if ( ffstring == NULL ) return NULL;

  if (ajp->index) {
    index = &asp->index;
  } else {
    index = NULL;
  }

  if (ajp->gbseq) {
    gbseq = &asp->gbseq;
  } else {
    gbseq = NULL;
  }

  if (! StringHasNoText (rbp->string)) return StringSave (rbp->string);

  /* could be descriptor, feature, or submit block citation */

  if (rbp->itemtype == OBJ_SEQDESC) {

    sdp = SeqMgrGetDesiredDescriptor (rbp->entityID, NULL, rbp->itemID, 0, NULL, &dcontext);
    if (sdp != NULL && dcontext.seqdesctype == Seq_descr_pub) {
      pdp = (PubdescPtr) sdp->data.ptrvalue;
    }

  } else if (rbp->itemtype == OBJ_SEQFEAT) {

    sfp = SeqMgrGetDesiredFeature (rbp->entityID, NULL, rbp->itemID, 0, NULL, &fcontext);
    if (sfp != NULL && fcontext.seqfeattype == SEQFEAT_PUB) {
      pdp = (PubdescPtr) sfp->data.value.ptrvalue;
    }

  } else if (rbp->itemtype == OBJ_SEQSUB_CIT) {

    omdp = ObjMgrGetData (rbp->entityID);
    if (omdp != NULL && omdp->datatype == OBJ_SEQSUB) {
      ssp = (SeqSubmitPtr) omdp->dataptr;
      if (ssp != NULL && ssp->datatype == 1) {
        sbp = ssp->sub;
        if (sbp != NULL) {
          csp = sbp->cit;
        }
      }
    }
  }

  if (pdp == NULL && csp == NULL) return NULL;

  temp = FFGetString(ajp);
  if ( temp == NULL ) {
    FFRecycleString(ajp, ffstring);
    return NULL;
  }

  /* print serial number */
  FFStartPrint(temp, afp->format, 0, 12, "REFERENCE", 12, 5, 5, "RN", TRUE);

  if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
    if (rbp->serial > 99) {
      sprintf (buf, "%d ", (int) rbp->serial);
    } else {
      sprintf (buf, "%d", (int) rbp->serial);
    }
  } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {
    sprintf (buf, "[%d]", (int) rbp->serial);
  }

  FFAddOneString(temp, buf, FALSE, FALSE, TILDE_TO_SPACES);

  /* print base range */

  if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {

    if (rbp->sites != 3) {
      FFAddNChar(temp, ' ', 15 - temp->pos, FALSE);
    }
  } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {

    if (rbp->sites == 0) {
      FFLineWrap(ffstring, temp, 0, 5, ASN2FF_EMBL_MAX, "RN");   
      FFRecycleString(ajp, temp);
      temp = FFGetString(ajp);
      FFStartPrint(temp, afp->format, 0, 0, NULL, 0, 5, 5, "RP", FALSE);
    }
  }

  if (rbp->sites == 1 || rbp->sites == 2) {

    if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
      FFAddOneString(temp, "(sites)", FALSE, FALSE, TILDE_TO_SPACES);
      FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
    } else {
      FFLineWrap(ffstring, temp, 5, 5, ASN2FF_EMBL_MAX, "RP");
    }
  } else if (rbp->sites == 3) {
    if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
      FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
    } else {
      FFLineWrap(ffstring, temp, 5, 5, ASN2FF_EMBL_MAX, "RP");
    }
  } else {
    if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
      FFAddNChar(temp, ' ', 15 - temp->pos, FALSE);
      if (afp->format == GENBANK_FMT) {
        FFAddOneString(temp, "(bases ", FALSE, FALSE, TILDE_TO_SPACES);
      } else {
        FFAddOneString(temp, "(residues ", FALSE, FALSE, TILDE_TO_SPACES);
      }
    }

    irp = (IntRefBlockPtr) rbp;
    loc = irp->loc;

    if (loc != NULL) {
      if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
        suffix = "; ";
      } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {
        suffix = ", ";
      }

      slp = SeqLocFindNext (loc, NULL);
      while (slp != NULL) {
        nextslp = SeqLocFindNext (loc, slp);
        start = SeqLocStart (slp) + 1;
        stop = SeqLocStop (slp) + 1;
        if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
          sprintf (buf, "%ld to %ld", (long) start, (long) stop);
        } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {
          sprintf (buf, "%ld-%ld", (long) start, (long) stop);
        }
        if (nextslp == NULL) {
          suffix = NULL;
        }
        FFAddTextToString(temp, NULL, buf, suffix, FALSE, FALSE, TILDE_TO_SPACES);
        slp = nextslp;
      }

    } else {

      /* code still used for ssp->cit */

      start = 1;
      stop = bsp->length;
      if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
        sprintf (buf, "%ld to %ld", (long) start, (long) stop);
      } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {
        sprintf (buf, "%ld-%ld", (long) start, (long) stop);
      }
      FFAddOneString(temp, buf, FALSE, FALSE, TILDE_TO_SPACES);
    }

    if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
      FFAddOneString(temp, ")", FALSE, FALSE, TILDE_TO_SPACES);
    }
    if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
      FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
    } else {
      FFLineWrap(ffstring, temp, 5, 5, ASN2FF_EMBL_MAX, "RP");
    }
  }

  if (gbseq != NULL) {
    gbref = GBReferenceNew ();
  }

  /* print author list */

  FFRecycleString(ajp, temp);
  temp = FFGetString(ajp);
  FFStartPrint(temp, afp->format, 2, 12, "AUTHORS", 12, 5, 5, "RA", FALSE);

  str = NULL;
  consortium = NULL;

  alp = GetAuthListPtr (pdp, csp);
  if (alp != NULL) {
    str = GetAuthorsString (afp->format, alp, &consortium, index, gbref);
    TrimSpacesAroundString (str);
  }

  if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
    suffix = NULL;
    trailingPeriod = TRUE;
  } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {
    trailingPeriod = FALSE;
    len = StringLen (str);
    if (len > 0 && str [len - 1] != '.') {
      suffix = ".;";
    } else {
      suffix = ";";
    }
  }

  /* if no authors were found, period will still be added by this call */
  if (str != NULL) {
    FFAddTextToString(temp, NULL, str, suffix, trailingPeriod, FALSE, TILDE_TO_SPACES);
  } else {
    if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
      FFAddOneChar(temp, '.', FALSE);
    } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {
      FFAddOneChar(temp, ';', FALSE);
    }    
  }

  MemFree (str);
  if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
      FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
  } else {
      FFLineWrap(ffstring, temp, 5, 5, ASN2FF_EMBL_MAX, "RA");
  }

  /* print consortium */

  FFRecycleString(ajp, temp);
  temp = FFGetString(ajp);
  if (! StringHasNoText (consortium)) {
    FFStartPrint (temp, afp->format, 2, 12, "CONSRTM", 12, 5, 5, "CM", FALSE);
    FFAddTextToString (temp, NULL, consortium, suffix, FALSE, FALSE, TILDE_TO_SPACES);
    if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
      FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
    } else {
      FFLineWrap(ffstring, temp, 5, 5, ASN2FF_EMBL_MAX, "CM");
    }
  }
  MemFree (consortium);

  /* print title */
  FFRecycleString(ajp, temp);
  temp = FFGetString(ajp);

  str = GetPubTitle (afp->format, pdp, csp);
  CleanPubTitle (str);
  StrStripSpaces (str);

  if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
    prefix = NULL;
    suffix = NULL;
  } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {
    if (str != NULL) {
      prefix = "\"";
      suffix = "\";";
    } else {
      prefix = NULL;
      suffix = ";";
    }
  }

  if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
    if (! StringHasNoText (str)) {
      FFStartPrint (temp, afp->format, 2, 12, "TITLE", 12, 5, 5, "RT", FALSE);

      FFAddTextToString (temp, prefix, str, suffix, FALSE, FALSE, TILDE_TO_SPACES);
      FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
    }
  } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {
    FFStartPrint (temp, afp->format, 2, 12, "TITLE", 12, 5, 5, "RT", FALSE);
    if (! StringHasNoText (str)) {

      FFAddTextToString (temp, prefix, str, suffix, FALSE, FALSE, TILDE_TO_SPACES);

    } else {
      FFAddOneChar (temp, ';', FALSE);
    }
    FFLineWrap(ffstring, temp, 5, 5, ASN2FF_EMBL_MAX, "RT");
  }

  if (gbseq != NULL) {
    if (gbref != NULL) {
      gbref->title = StringSaveNoNull (str);
    }
  }

  MemFree (str);

  /* print journal */
  FFRecycleString(ajp, temp);
  temp = FFGetString(ajp);

  FFStartPrint (temp, afp->format, 2, 12, "JOURNAL", 12, 5, 5, "RL", FALSE);

  /* Only GenBank/EMBL/DDBJ require ISO JTA in ENTREZ/RELEASE modes (RefSeq should later) */

  citArtIsoJta = ajp->flags.citArtIsoJta;
  strict_isojta = FALSE;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_GENBANK ||
        sip->choice == SEQID_EMBL ||
        sip->choice == SEQID_DDBJ ||
        /* sip->choice == SEQID_OTHER || */
        sip->choice == SEQID_TPG ||
        sip->choice == SEQID_TPE ||
        sip->choice == SEQID_TPD) {
      strict_isojta = TRUE;
    }
  }
  if (! strict_isojta) {
    citArtIsoJta = FALSE;
  }

  str = GetPubJournal (afp->format, ajp->flags.dropBadCitGens,
                       ajp->flags.noAffilOnUnpub, citArtIsoJta,
                       pdp, csp, bsp->id, index);
  if (str == NULL) {
    str = StringSave ("Unpublished");
  }
  StrStripSpaces (str);
  TrimSpacesAroundString (str);

  if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
    needsPeriod = FALSE;
  } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {
    needsPeriod = TRUE;
  }

  FFAddOneString (temp, str, FALSE, FALSE, TILDE_IGNORE);
  if (needsPeriod) {
    FFAddOneChar(temp, '.', FALSE);
  }

  if (gbseq != NULL) {
    if (gbref != NULL) {
      gbref->journal = StringSaveNoNull (str);
      tmp = gbref->journal;
      if (tmp != NULL) {
        ch = *tmp;
        while (ch != '\0') {
          if (ch == '\n' || ch == '\r' || ch == '\t') {
            *tmp = ' ';
          }
          tmp++;
          ch = *tmp;
        }
        TrimSpacesAroundString (gbref->journal);
      }
    }
  }

  MemFree (str);
  if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
    FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
  } else {
    FFLineWrap(ffstring, temp, 5, 5, ASN2FF_EMBL_MAX, "RL");
  }

  /* print muid */
  FFRecycleString(ajp, temp);
  temp = FFGetString(ajp);
  
  muid = GetMuid (pdp);
  if (muid > 0) {
    FFStartPrint (temp, afp->format, 2, 12, "MEDLINE", 12, 5, 5, "RX", FALSE);

    if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
      FF_www_muid (ajp, temp, muid);
      FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
    } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {
      sprintf (buf, "MEDLINE; %ld.", (long) muid);
      FFAddOneString (temp, buf, FALSE, FALSE, TILDE_TO_SPACES);
      FFLineWrap(ffstring, temp, 5, 5, ASN2FF_EMBL_MAX, "RX");
    }
  }

  FFRecycleString(ajp, temp);
  temp = FFGetString(ajp);

  if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
    pmid = GetPmid (pdp);
    if (pmid > 0) {
      FFStartPrint (temp, afp->format, 3, 12, "PUBMED", 12, 5, 5, "RX", FALSE);

      FF_www_muid (ajp, temp, pmid);
      FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
    }
  }

  if (gbseq != NULL) {
    if (gbref != NULL) {
      gbref->medline = muid;
      gbref->pubmed = pmid;
    }
  }

  if (pdp == NULL) {
    str = FFToCharPtr(ffstring);

    if (gbseq != NULL) {
      if (gbref != NULL) {
        AddReferenceToGbseq (gbseq, gbref, str);
      }
    }

    FFRecycleString(ajp, ffstring);
    FFRecycleString(ajp, temp);
    return str;
  }


  /* !!! remainder of fields are only for GenBank !!! */

  if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {

    prefix = "REMARK";

    if (pdp->comment != NULL) {
      for (i = 0, notFound = TRUE; notFound && remarksText [i] != NULL; i++) {
        if (StringCmp (pdp->comment, remarksText [i]) == 0) {
          notFound = FALSE;
        }
      }
      if (notFound) {
        FFRecycleString(ajp, temp);
        temp = FFGetString(ajp);

        FFStartPrint (temp, afp->format, 2, 12, prefix, 12, 5, 5, NULL, FALSE);
        FFAddOneString (temp, pdp->comment, FALSE, TRUE, TILDE_EXPAND);
        FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
        prefix = NULL;
      }
    }

    gibbsq = 0;
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_GIBBSQ) {
        gibbsq = sip->data.intvalue;
      }
    }
    csp = NULL;
    for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
      if (vnp->choice == PUB_Sub) {
        csp = (CitSubPtr) vnp->data.ptrvalue;
      }
    }
    if (gibbsq > 0 /* && csp == NULL */) {
      FFRecycleString(ajp, temp);
      temp = FFGetString(ajp);

      sprintf (buf, "GenBank staff at the National Library of Medicine created this entry [NCBI gibbsq %ld] from the original journal article.", (long) gibbsq);
      FFStartPrint (temp, afp->format, 2, 12, prefix, 12, 5, 5, NULL, FALSE);
      FFAddOneString (temp, buf, FALSE, FALSE, TILDE_EXPAND);
      FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
      prefix = NULL;

      /* gibbsq comment section (fields may be copied from degenerate pubdesc) */

      str = pdp->fig;
      if (StringHasNoText (str)) {
        str = irp->fig;
      }
      if (! StringHasNoText (str)) {
        FFRecycleString(ajp, temp);
        temp = FFGetString(ajp);

        sprintf (buf, "This sequence comes from %s", str);
        FFStartPrint (temp, afp->format, 2, 12, prefix, 12, 5, 5, NULL, FALSE);
        FFAddOneString (temp, buf, TRUE, TRUE, TILDE_EXPAND);
        FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
        prefix = NULL;
      }

      if (pdp->poly_a || irp->poly_a) {
        FFRecycleString(ajp, temp);
        temp = FFGetString(ajp);

        FFStartPrint (temp ,afp->format, 2, 12, prefix, 12, 5, 5, NULL, FALSE);
        FFAddOneString (temp, "Polyadenylate residues occurring in the figure were omitted from the sequence.", TRUE, TRUE, TILDE_EXPAND);
        FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
        prefix = NULL;
      }

      str = pdp->maploc;
      if (StringHasNoText (str)) {
        str = irp->maploc;
      }
      if (! StringHasNoText (str)) {
        FFRecycleString(ajp, temp);
        temp = FFGetString(ajp);

        sprintf (buf, "Map location: %s", str);
        FFStartPrint (temp, afp->format, 2, 12, prefix, 12, 5, 5, NULL, FALSE);
        FFAddOneString (temp, buf, TRUE, TRUE, TILDE_EXPAND);
        FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
        prefix = NULL;
      }

    }

    for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
      if (vnp->choice == PUB_Article) {
        cap = (CitArtPtr) vnp->data.ptrvalue;
        if (cap != NULL && cap->from == 1) {
          cjp = (CitJourPtr) cap->fromptr;
          if (cjp != NULL) {
            imp = cjp->imp;
            if (imp != NULL) {
              crp = imp->retract;
              if (crp != NULL && crp->type == 3) {
                FFRecycleString(ajp, temp);
                temp = FFGetString(ajp);

                FFStartPrint (temp, afp->format, 2, 12, prefix, 12, 5, 5, NULL, FALSE);
                FFAddOneString (temp, "Erratum:", FALSE, FALSE, TILDE_TO_SPACES);
                FFAddTextToString (temp, "[", crp->exp, "]", FALSE, TRUE, TILDE_EXPAND);
                FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
                prefix = NULL;
              }
            }
          }
        }
      } else if (vnp->choice == PUB_Sub) {
        csp = (CitSubPtr) vnp->data.ptrvalue;
        if (csp != NULL) {
          if (! StringHasNoText (csp->descr)) {
            FFRecycleString(ajp, temp);
            temp = FFGetString(ajp);

            FFStartPrint (temp, afp->format, 2, 12, prefix, 12, 5, 5, NULL, FALSE);
            FFAddOneString (temp, csp->descr, FALSE, TRUE, TILDE_EXPAND);
            FFLineWrap(ffstring, temp, 12, 12, ASN2FF_GB_MAX, NULL);
            prefix = NULL;
          }
        }
      }
    }

  }

  str = FFToCharPtr(ffstring);

  if (gbseq != NULL) {
    if (gbref != NULL) {
      AddReferenceToGbseq (gbseq, gbref, str);
    }
  }

  FFRecycleString(ajp, ffstring);
  FFRecycleString(ajp, temp);
  return str;
}

/* A tilde is not an EOL if it is found in a sring of the form:    */
/* /~alpahnumdot/ where alphanumdot is either alpha numeric or '.' */
/*                                                                 */
/* str points to the tilde in question.                            */
static Boolean IsTildeEOL(CharPtr str) {
  CharPtr ptr;

  if ( *(str - 1) != '/' ) return TRUE;

  ++str;

  
  for ( ptr = str; 
	IS_ALPHANUM(*ptr) || *ptr == '_' || *ptr == '-' || *ptr == '.';
	++ptr) continue;

  return *ptr == '/' ? FALSE : TRUE;
}


/* returns a pointer to the first character past the url */
static CharPtr FindUrlEnding(CharPtr str) {
  CharPtr ptr;

  for ( ptr = str;
        !IS_WHITESP(*ptr) && *ptr != '\0' && *ptr != '(' && *ptr != '\"';
        ++ptr  ) {
    if ( *ptr == '~' ) {
      if ( IsTildeEOL(ptr) ) break;
    }
  }

  --ptr;

  /* back up over any trailing periods, commas, or parentheses */
  while ( (*ptr == '.') || (*ptr == ',') || (*ptr == ')') ) {
    --ptr;
  }

  ++ptr;

  return ptr;
}

static void AddCommentWithURLlinks (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  CharPtr prefix,
  CharPtr str,
  CharPtr suffix
)

{
  Char     ch;
  CharPtr  ptr;

  while (! StringHasNoText (str)) {
    ptr = StringStr (str, "http://");
    if (ptr == NULL) {
      if (prefix != NULL) {
        FFAddOneString(ffstring, prefix, FALSE, FALSE, TILDE_IGNORE);
      }
      AddStringWithTildes (ffstring, str);
      if (suffix != NULL) {
        FFAddOneString(ffstring, suffix, FALSE, FALSE, TILDE_IGNORE);
      }
      return;
    }

    *ptr = '\0';
    AddStringWithTildes (ffstring, str); 
    *ptr = 'h';

    str = ptr;
    ptr = FindUrlEnding(str);


    ch = *ptr;
    *ptr = '\0';
    if ( GetWWW(ajp) ) {
      FFAddTextToString(ffstring, "<a href=", str, ">", FALSE, FALSE, TILDE_IGNORE);
      FFAddTextToString(ffstring, NULL, str, "</a>", FALSE, FALSE, TILDE_IGNORE);
    } else {
      FFAddOneString(ffstring, str, FALSE, FALSE, TILDE_IGNORE);
    }

    *ptr = ch;
    str = ptr;
  }
}

static CharPtr Asn2gnbkCompressSpaces (CharPtr str)

{
  Char     ch;
  CharPtr  dst;
  Char     last;
  CharPtr  ptr;

  if (str != NULL && str [0] != '\0') {
    dst = str;
    ptr = str;
    ch = *ptr;
    while (ch != '\0' && ch <= ' ') {
      ptr++;
      ch = *ptr;
    }
    while (ch != '\0') {
      *dst = ch;
      dst++;
      ptr++;
      last = ch;
      ch = *ptr;
      if (ch != '\0' && ch < ' ') {
        *ptr = ' ';
        ch = *ptr;
      }
      while (ch != '\0' && last <= ' ' && ch <= ' ') {
        ptr++;
        ch = *ptr;
      }
    }
    *dst = '\0';
    dst = NULL;
    ptr = str;
    ch = *ptr;
    while (ch != '\0') {
      if (ch != ' ') {
        dst = NULL;
      } else if (dst == NULL) {
        dst = ptr;
      }
      ptr++;
      ch = *ptr;
    }
    if (dst != NULL) {
      *dst = '\0';
    }
  }
  return str;
}

static void CatenateCommentInGbseq (
  GBSeqPtr gbseq,
  CharPtr str
)

{
  Char     ch;
  CharPtr  tmp;

  if (gbseq == NULL || StringHasNoText (str)) return;

  if (StringNCmp (str, "COMMENT     ", 12) == 0) {
    str += 12;
  }
  if (gbseq->comment == NULL) {
    gbseq->comment = StringSave (str);
  } else {
    tmp = (CharPtr) MemNew (StringLen (gbseq->comment) + StringLen (str) + 4);
    StringCpy (tmp, gbseq->comment);
    StringCat (tmp, "; ");
    StringCat (tmp, str);
    gbseq->comment = MemFree (gbseq->comment);
    gbseq->comment = tmp;
  }

  tmp = gbseq->comment;
  if (tmp == NULL) return;
  ch = *tmp;
  while (ch != '\0') {
    if (ch == '\n' || ch == '\r' || ch == '\t') {
      *tmp = ' ';
    }
    tmp++;
    ch = *tmp;
  }
  TrimSpacesAroundString (gbseq->comment);
  Asn2gnbkCompressSpaces (gbseq->comment);
}


static CharPtr FormatCommentBlock (
  Asn2gbFormatPtr afp,
  BaseBlockPtr bbp
)

{
  Boolean            add_period;
  IntAsn2gbJobPtr    ajp;
  Asn2gbSectPtr      asp;
  CommentBlockPtr    cbp;
  CharPtr            db;
  DbtagPtr           dbt;
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext;
  GBSeqPtr           gbseq;
  size_t             len;
  ObjectIdPtr        oip;
  CharPtr            prefix;
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp;
  Char               sfx [32];
  CharPtr            str;
  CharPtr            suffix;
  CharPtr            title;
  StringItemPtr      ffstring;

  if (afp == NULL || bbp == NULL) return NULL;
  ajp = afp->ajp;
  if (ajp == NULL) return NULL;
  asp = afp->asp;
  if (asp == NULL) return NULL;

  cbp = (CommentBlockPtr) bbp;

  /* optionally populate gbseq for XML-ized GenBank format */

  if (ajp->gbseq) {
    gbseq = &asp->gbseq;
  } else {
    gbseq = NULL;
  }

  /* some comments are allocated (along with possible first COMMENT label) */

  if (! StringHasNoText (bbp->string)) {
    str = StringSave (bbp->string);
    CatenateCommentInGbseq (gbseq, str);
    return str;
  }

  title = NULL;
  prefix = NULL;
  suffix = NULL;
  add_period = FALSE;
  sfx [0] = '\0';

  if (bbp->itemtype == OBJ_SEQDESC) {

    /* usually should reference comment, maploc, or region descriptor IDs */

    sdp = SeqMgrGetDesiredDescriptor (bbp->entityID, NULL, bbp->itemID, 0, NULL, &dcontext);
    if (sdp != NULL) {

      if (dcontext.seqdesctype == Seq_descr_comment) {

        title = (CharPtr) sdp->data.ptrvalue;

      } else if (dcontext.seqdesctype == Seq_descr_maploc) {

        dbt = (DbtagPtr) sdp->data.ptrvalue;
        if (dbt != NULL) {
          db = dbt->db;
          oip = dbt->tag;
          if (oip != NULL) {
            if (oip->str != NULL) {

              title = oip->str;
              prefix = ("Map location: ");

            } else if (db != NULL && oip->id != 0) {

              title = db;
              prefix = ("Map location: (Database ");
              sprintf (sfx, "; id # %ld).", (long) oip->id);
              suffix = sfx;

            }
          }
        }

      } else if (dcontext.seqdesctype == Seq_descr_region) {

        title = (CharPtr) sdp->data.ptrvalue;
        prefix = "Region: ";

      }
    }

  } else if (bbp->itemtype == OBJ_SEQFEAT) {

    /* also have to deal with comment feature across entire sequence */

    sfp = SeqMgrGetDesiredFeature (bbp->entityID, NULL, bbp->itemID, 0, NULL, &fcontext);
    if (sfp != NULL && fcontext.seqfeattype == SEQFEAT_COMMENT) {

      title = sfp->comment;
    }
  }

  if (title == NULL) return NULL;

  ffstring = FFGetString(ajp);
  if ( ffstring == NULL ) return NULL;

  if (cbp->first) {
    FFStartPrint (ffstring, afp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
  } else {
    FFStartPrint (ffstring, afp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
  }

  str = StringSave (title);
  TrimSpacesAndJunkFromEnds (str, TRUE);
  if (! IsEllipsis (str)) {
    s_RemovePeriodFromEnd (str);
    len = StringLen (str);
    if (len > 0 && str [len - 1] != '.') {
      add_period = TRUE;
    }
  }
  AddCommentWithURLlinks(ajp, ffstring, prefix, str, suffix);
  /*
  if ( GetWWW(ajp) && prefix == NULL && suffix == NULL) {
    
    AddCommentWithURLlinks (ffstring, str);
  } else {
    FFAddTextToString (ffstring, prefix, str, suffix, FALSE, TRUE, TILDE_OLD_EXPAND);
  }
  */
  if (add_period) {
    FFAddOneChar (ffstring, '.',FALSE);
  }
  MemFree (str);

  str = FFEndPrint(ajp, ffstring, afp->format, 12, 12, 5, 5, "CC");

  CatenateCommentInGbseq (gbseq, str);

  FFRecycleString(ajp, ffstring);
  return str;
}

/* format features section */

static Boolean is_real_id (
  SeqIdPtr sip,
  SeqIdPtr this_sip
)

{
  BioseqPtr  bsp;

  if (sip == NULL || this_sip == NULL) return FALSE;

  if (! SeqIdIn (sip, this_sip)) {
    bsp = BioseqFind (sip);
    if (bsp == NULL) return TRUE;  /* ??? */
    if (bsp->repr == Seq_repr_virtual) return FALSE;
  }

  return TRUE;
}

static Boolean FlatVirtLoc (
  BioseqPtr bsp,
  SeqLocPtr location
)

{
  SeqIntPtr  sintp;
  SeqIdPtr   sip;
  SeqPntPtr  spp;

  if (bsp == NULL || location == NULL) return FALSE;

  switch (location->choice) {
    case SEQLOC_WHOLE :
      sip = (SeqIdPtr) location->data.ptrvalue;
      if (sip == NULL) return TRUE;
      if (! is_real_id (sip, bsp->id)) return TRUE;
      break;
    case SEQLOC_INT :
      sintp = (SeqIntPtr) location->data.ptrvalue;
      if (sintp == NULL) return TRUE;
      sip = sintp->id;
      if (sip == NULL) return TRUE;
      if (! is_real_id (sip, bsp->id)) return TRUE;
      break;
    case SEQLOC_PNT :
      spp = (SeqPntPtr) location->data.ptrvalue;
      if (spp == NULL) return TRUE;
      sip = spp->id;
      if (sip == NULL) return TRUE;
      if (! is_real_id (sip, bsp->id)) return TRUE;
      break;
    default :
      break;
  }

  return FALSE;
}

static Uint1    order [NUM_SEQID];
static Boolean  order_initialized = FALSE;

static CharPtr lim_str [5] = { "", ">","<", ">", "<" };

static Boolean GetAccnVerFromServer (Int4 gi, CharPtr buf)

{
  AccnVerLookupFunc  func;
  SeqMgrPtr          smp;
  CharPtr            str;

  if (buf == NULL) return FALSE;
  *buf = '\0';
  smp = SeqMgrWriteLock ();
  if (smp == NULL) return FALSE;
  func = smp->accn_ver_lookup_func;
  SeqMgrUnlock ();
  if (func == NULL) return FALSE;
  str = (*func) (gi);
  if (str == NULL) return FALSE;
  if (StringLen (str) < 40) {
    StringCpy (buf, str);
  }
  MemFree (str);
  return TRUE;
}


/******************************************************************************/
/*                              FlatLoc functions  .                          */
/******************************************************************************/

static Boolean FF_FlatNullAhead (
  BioseqPtr bsp,
  ValNodePtr location
)

{
  SeqLocPtr  next;

  if (bsp == NULL || location == NULL) return FALSE;

  next = location->next;
  if (next == NULL) return TRUE;
  if (next->choice == SEQLOC_NULL) return TRUE;
  if (FlatVirtLoc (bsp, next)) return TRUE;

  return FALSE;
}



static void FlatLocSeqId (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  SeqIdPtr sip
)

{
  BioseqPtr    bsp;
  Char         buf [40];
  ObjectIdPtr  oip;
  SeqIdPtr     use_id = NULL;
  Boolean      was_lock = FALSE;

  if (ffstring == NULL || sip == NULL) return;

  buf [0] = '\0';
  bsp = BioseqFind (sip);
  if (bsp != NULL) {
    use_id = SeqIdSelect (bsp->id, order, NUM_SEQID);
  } else if (sip->choice == SEQID_GI) {
    if (GetAccnVerFromServer (sip->data.intvalue, buf)) {
      FFAddTextToString(ffstring, NULL, buf, ":", FALSE, FALSE, TILDE_IGNORE);
      /*AddValNodeString (head, NULL, buf, ":");*/
      return;
    }
    use_id = GetSeqIdForGI (sip->data.intvalue);
  }
  if (use_id == NULL && bsp == NULL) {
    bsp = BioseqLockById (sip);
    was_lock = TRUE;
    if (bsp != NULL) {
      use_id = SeqIdSelect (bsp->id, order, NUM_SEQID);
    }
  }
  if (use_id != NULL) {
    SeqIdWrite (use_id, buf, PRINTID_TEXTID_ACC_VER, sizeof (buf) - 1);
    if (use_id->choice == SEQID_GI) {
      ajp->relModeError = TRUE;
    }
  } else if (sip->choice == SEQID_GI) {
    SeqIdWrite (sip, buf, PRINTID_FASTA_LONG, sizeof (buf) - 1);
    ajp->relModeError = TRUE;
  } else {
    SeqIdWrite (sip, buf, PRINTID_TEXTID_ACC_VER, sizeof (buf) - 1);
    if (sip->choice == SEQID_GI) {
      ajp->relModeError = TRUE;
    }
  }
  if (was_lock) {
    BioseqUnlock (bsp);
  }
  if (StringHasNoText (buf)) {
    StringCpy (buf, "?00000");
    ajp->relModeError = TRUE;
    if (use_id != NULL && use_id->choice == SEQID_LOCAL) {
      oip = (ObjectIdPtr) use_id->data.ptrvalue;
      if (oip != NULL && (! StringHasNoText (oip->str))) {
        StringNCpy_0 (buf, oip->str, 13);
      }
    }
  }
  FFAddTextToString(ffstring, NULL, buf, ":", FALSE, FALSE, TILDE_IGNORE);
}



static void FlatLocCaret (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  SeqIdPtr sip,
  SeqIdPtr this_sip,
  Int4 point,
  IntFuzzPtr fuzz
)

{
  Char   buf [128];
  Uint1  index;

  if (ffstring == NULL) return;

  if (sip != NULL && (! SeqIdIn (sip, this_sip))) {
    FlatLocSeqId (ajp, ffstring, sip);
  }

  buf [0] = '\0';
  point++; /* orginal FlatLocHalfCaret was called with point + 1 */

  if (fuzz != NULL) {
    switch (fuzz->choice) {
      case 1 :
        sprintf (buf, "(%ld.%ld)..(%ld.%ld)",
                 (long) (point - fuzz->a),
                 (long) point,
                 (long) point,
                 (long) (point + fuzz->a));
        break;
      case 2 :
        sprintf (buf, "%ld^%ld",
                 (long) (1 + fuzz->b),
                 (long) (1 + fuzz->a));
        break;
      case 3 :
        sprintf (buf, "%ld^%ld",
                 (long) (point - point * ((double) fuzz->a / 1000.0)),
                 (long) (point + point * ((double) fuzz->a / 1000.0)));
        break;
      case 4 :
        if (fuzz->a == 3) { /* space to right */
          sprintf (buf, "%ld^%ld", (long) (point), (long) (point + 1));
        } else if (fuzz->a == 4 && point > 1) { /* space to left */
          sprintf (buf, "%ld^%ld", (long) (point - 1), (long) point);
        } else {
          index = (Uint1) fuzz->a;
          if (index > 4) {
            index = 0;
          }
          sprintf (buf, "%s%ld", lim_str [index], (long) point);
        }
        break;
      default :
        sprintf (buf, "%ld", (long) point);
        break;
    }
  } else {
    sprintf (buf, "%ld", (long) point);
  }

  FFAddOneString(ffstring, buf, FALSE, FALSE, TILDE_IGNORE);
}


static void FlatLocPoint (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  SeqIdPtr sip,
  SeqIdPtr this_sip,
  Int4 point,
  IntFuzzPtr fuzz
)

{
  Char   buf [128];
  Uint1  index;

  if (ffstring == NULL) return;

  if (sip != NULL && (! SeqIdIn (sip, this_sip))) {
    FlatLocSeqId (ajp, ffstring, sip);
  }

  buf [0] = '\0';
  point++;

  if (fuzz != NULL) {
    switch (fuzz->choice) {
      case 1 :
        sprintf (buf, "(%ld.%ld)",
                 (long) (point - fuzz->a),
                 (long) (point + fuzz->a));
        break;
      case 2 :
        sprintf (buf, "(%ld.%ld)",
                 (long) (1 + fuzz->b),
                 (long) (1 + fuzz->a));
        break;
      case 3 :
        sprintf (buf, "(%ld.%ld)",
                 (long) (point - point * ((double) fuzz->a / 1000.0)),
                 (long) (point + point * ((double) fuzz->a / 1000.0)));
        break;
      case 4 :
        index = (Uint1) fuzz->a;
        if (index > 4) {
          index = 0;
        }
        sprintf (buf, "%s%ld", lim_str [index], (long) point);
        break;
      default :
        sprintf (buf, "%ld", (long) point);
        break;
    }
  } else {
    sprintf (buf, "%ld", (long) point);
  }

  FFAddOneString(ffstring, buf, FALSE, FALSE, TILDE_IGNORE);
}


static void FlatLocElement (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  BioseqPtr bsp,
  SeqLocPtr location
)

{
  Boolean     minus_strand = FALSE;
  SeqBondPtr  sbp;
  SeqIntPtr   sintp;
  SeqIdPtr    sip;
  SeqPntPtr   spp;
  BioseqPtr   wholebsp;

  if (ffstring == NULL || bsp == NULL || location == NULL) return;

  switch (location->choice) {
    case SEQLOC_WHOLE :
      sip = (SeqIdPtr) location->data.ptrvalue;
      if (sip == NULL) return;
      wholebsp = BioseqFind (sip);
      if (wholebsp == NULL) return;
      if (is_real_id (sip, bsp->id)) {
        FlatLocPoint (ajp, ffstring, sip, bsp->id, 0, NULL);
        if (bsp->length > 0) {
          FFAddOneString(ffstring, "..", FALSE, FALSE, TILDE_IGNORE);
          FlatLocPoint (ajp, ffstring, NULL, bsp->id, bsp->length - 1, NULL);
        }
      }
      break;
    case SEQLOC_INT :
      sintp = (SeqIntPtr) location->data.ptrvalue;
      if (sintp == NULL) return;
      sip = sintp->id;
      if (sip == NULL) return;
      if (is_real_id (sip, bsp->id)) {
        minus_strand = (Boolean) (sintp->strand == Seq_strand_minus);
        if (minus_strand) {
          FFAddOneString(ffstring, "complement(", FALSE, FALSE, TILDE_IGNORE);
        }
        FlatLocPoint (ajp, ffstring, sip, bsp->id, sintp->from, sintp->if_from);
        if (sintp->to > 0 &&
            (sintp->to != sintp->from ||
             sintp->if_from != NULL ||
             sintp->if_to != NULL)) {
          FFAddOneString(ffstring, "..", FALSE, FALSE, TILDE_IGNORE);
          FlatLocPoint (ajp, ffstring, NULL, bsp->id, sintp->to, sintp->if_to);
        }
        if (minus_strand) {
          FFAddOneString(ffstring, ")", FALSE, FALSE, TILDE_IGNORE);
        }
      }
      break;
    case SEQLOC_PNT :
      spp = (SeqPntPtr) location->data.ptrvalue;
      if (spp == NULL) return;
      sip = spp->id;
      if (sip == NULL) return;
      if (is_real_id (sip, bsp->id)) {
        minus_strand = (Boolean) (spp->strand == Seq_strand_minus);
        if (minus_strand) {
          FFAddOneString(ffstring, "complement(", FALSE, FALSE, TILDE_IGNORE);
        }
        if (spp->fuzz != NULL) {
          FlatLocCaret (ajp, ffstring, sip, bsp->id, spp->point, spp->fuzz);
        } else {
          FlatLocPoint (ajp, ffstring, sip, bsp->id, spp->point, NULL);
        }
        if (minus_strand) {
          FFAddOneString(ffstring, ")", FALSE, FALSE, TILDE_IGNORE);
        }
      }
      break;
    case SEQLOC_BOND :
      sbp = (SeqBondPtr) location->data.ptrvalue;
      if (sbp == NULL) return;
      spp = sbp->a;
      if (spp == NULL) return;
      sip = spp->id;
      if (sip == NULL) return;
      FFAddOneString(ffstring, "bond(", FALSE, FALSE, TILDE_IGNORE);
      FlatLocPoint (ajp, ffstring, sip, bsp->id, spp->point, spp->fuzz);
      spp = sbp->b;
      if (spp != NULL) {
        FFAddOneString(ffstring, ",", FALSE, FALSE, TILDE_IGNORE);
        FlatLocPoint (ajp, ffstring, NULL, bsp->id, spp->point, spp->fuzz);
      }
      FFAddOneString(ffstring, ")", FALSE, FALSE, TILDE_IGNORE);
      break;
    default :
      /* unexpected internal complex type or unimplemented SEQLOC_FEAT */
      return;
  }
}



static void FF_FlatPackedPoint (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  PackSeqPntPtr pspp,
  BioseqPtr bsp
)

{
  Uint1  dex;

  if (ffstring == NULL || pspp == NULL || bsp == NULL) return;

  for (dex = 0; dex < pspp->used; dex++) {
    FlatLocPoint (ajp, ffstring, pspp->id, bsp->id, pspp->pnts [dex], pspp->fuzz);
  }
}


static void FF_DoFlatLoc (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  BioseqPtr bsp,
  SeqLocPtr location,
  Boolean ok_to_complement
);

static void FF_GroupFlatLoc (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  BioseqPtr bsp,
  SeqLocPtr location,
  CharPtr prefix,
  Boolean is_flat_order
)

{
  Boolean        found_non_virt = FALSE;
  SeqIdPtr       hold_next;
  Int2           parens = 1;
  PackSeqPntPtr  pspp;
  SeqLocPtr      slp;
  Boolean        special_mode = FALSE; /* join in order */

  if (ffstring == NULL || bsp == NULL || location == NULL) return;

  /* prefix will have the first parenthesis */

  FFAddOneString(ffstring, prefix, FALSE, FALSE, TILDE_IGNORE);

  for (slp = (SeqLocPtr) location->data.ptrvalue; slp != NULL; slp = slp->next) {

    if (slp->choice == SEQLOC_NULL || FlatVirtLoc (bsp, slp)) {
      if (slp != location && slp->next != NULL) {
        if (special_mode) {
          special_mode = FALSE;
          FFAddOneString(ffstring, ")", FALSE, FALSE, TILDE_IGNORE);
          parens--;
        }
      }
      continue;
    }

    if (found_non_virt && slp->choice != SEQLOC_EMPTY && slp->choice != SEQLOC_NULL) {
      FFAddOneString(ffstring, ",", FALSE, FALSE, TILDE_IGNORE);
    }

    switch (slp->choice) {
      case SEQLOC_WHOLE :
      case SEQLOC_PNT :
      case SEQLOC_BOND :
      case SEQLOC_FEAT :
        found_non_virt = TRUE;
        if (FlatVirtLoc (bsp, slp)) {
          if (slp != location && slp->next != NULL) {
            if (special_mode) {
              special_mode = FALSE;
              FFAddOneString(ffstring, ")", FALSE, FALSE, TILDE_IGNORE);
              parens--;
            }
          }
        } else {
          FlatLocElement (ajp, ffstring, bsp, slp);
        }
        break;
      case SEQLOC_INT :
        found_non_virt = TRUE;
        if (is_flat_order && (! FF_FlatNullAhead (bsp, slp))) {
          special_mode = TRUE;
          FFAddOneString(ffstring, "join(", FALSE, FALSE, TILDE_IGNORE);
          parens++;
        }
        FlatLocElement (ajp, ffstring, bsp, slp);
        break;
      case SEQLOC_PACKED_PNT :
        found_non_virt = TRUE;
        pspp = (PackSeqPntPtr) slp->data.ptrvalue;
        if (pspp != NULL) {
          FF_FlatPackedPoint (ajp, ffstring, pspp, bsp);
        }
        break;
      case SEQLOC_PACKED_INT :
      case SEQLOC_MIX :
      case SEQLOC_EQUIV :
        found_non_virt = TRUE;
        hold_next = slp->next;
        slp->next = NULL;
        FF_DoFlatLoc (ajp, ffstring, bsp, slp, FALSE);
        slp->next = hold_next;
        break;
      default :
        break;
    }

  }

  while (parens > 0) {
    FFAddOneString(ffstring, ")", FALSE, FALSE, TILDE_IGNORE);
    parens--;
  }
}




static void FF_DoFlatLoc (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  BioseqPtr bsp,
  SeqLocPtr location,
  Boolean ok_to_complement
)

{
  Boolean        found_null;
  SeqLocPtr      next_loc;
  PackSeqPntPtr  pspp;
  SeqLocPtr      slp;

  if (ffstring == NULL || bsp == NULL || location == NULL) return;

  /* deal with complement of entire location */

  if (ok_to_complement && SeqLocStrand (location) == Seq_strand_minus) {
    slp = AsnIoMemCopy ((Pointer) location,
                        (AsnReadFunc) SeqLocAsnRead,
                        (AsnWriteFunc) SeqLocAsnWrite);
    if (slp != NULL) {
      SeqLocRevCmp (slp);
      FFAddOneString(ffstring, "complement(", FALSE, FALSE, TILDE_IGNORE);
      FF_DoFlatLoc (ajp, ffstring, bsp, slp, FALSE);
      FFAddOneString(ffstring, ")", FALSE, FALSE, TILDE_IGNORE);
    }
    SeqLocFree (slp);
    return;
  }

  /* handle each location component */

  for (slp = location; slp != NULL; slp = slp->next) {

    if (slp->choice == SEQLOC_NULL || FlatVirtLoc (bsp, slp)) continue;

    /* print comma between components */

    if (slp != location) {
      FFAddOneString(ffstring, ",", FALSE, FALSE, TILDE_IGNORE);
    }

    switch (slp->choice) {
      case SEQLOC_MIX :
      case SEQLOC_PACKED_INT :
        found_null = FALSE;
        for (next_loc = (SeqLocPtr) slp->data.ptrvalue;
         next_loc != NULL;
         next_loc = next_loc->next) {
          if (next_loc->choice == SEQLOC_NULL ||
              FlatVirtLoc (bsp, next_loc) /* ||
              LocationHasNullsBetween (slp) */ )
            found_null = TRUE;
        }
        if (found_null) {
          FF_GroupFlatLoc (ajp, ffstring, bsp, slp, "order(", TRUE);
        } else {
          FF_GroupFlatLoc (ajp, ffstring, bsp, slp, "join(", FALSE);
        }
        break;
      case SEQLOC_EQUIV :
        FF_GroupFlatLoc (ajp, ffstring, bsp, slp, "one-of(", FALSE);
        break;
      case SEQLOC_PACKED_PNT :
        pspp = (PackSeqPntPtr) slp->data.ptrvalue;
        if (pspp != NULL) {
          FF_FlatPackedPoint (ajp, ffstring, pspp, bsp);
        }
        break;
      default :
        FlatLocElement (ajp, ffstring, bsp, slp);
        break;
    }

  }
}




static CharPtr FlatLoc (
  IntAsn2gbJobPtr ajp,
  BioseqPtr bsp,
  SeqLocPtr location,
  Boolean masterStyle
)

{
  Boolean     hasNulls;
  IntFuzzPtr  fuzz = NULL;
  SeqLocPtr   loc;
  Boolean     noLeft;
  Boolean     noRight;
  Uint1       num = 1;
  SeqPntPtr   spp;
  CharPtr     str;
  SeqLocPtr   tmp;
  StringItemPtr ffstring = NULL;

  if (bsp == NULL || location == NULL) return NULL;

  ffstring = FFGetString(ajp);

  if (! order_initialized) {
    order [SEQID_GENBANK] = num++;
    order [SEQID_EMBL] = num++;
    order [SEQID_DDBJ] = num++;
    order [SEQID_LOCAL] = num++;
    order [SEQID_OTHER] = num++;
    order [SEQID_TPG] = num++;
    order [SEQID_TPE] = num++;
    order [SEQID_TPD] = num++;
    order [SEQID_GIBBSQ] = num++;
    order [SEQID_GIBBMT] = num++;
    order [SEQID_PRF] = num++;
    order [SEQID_PDB] = num++;
    order [SEQID_PIR] = num++;
    order [SEQID_SWISSPROT] = num++;
    order [SEQID_PATENT] = num++;
    order [SEQID_GI] = num++;;
    order [SEQID_GENERAL] = num++;
    order [SEQID_GIIM] = num++;
    order_initialized = TRUE;
  }

  if (masterStyle) {

    /* map location from parts to segmented bioseq */

    if (location->choice == SEQLOC_PNT) {
      spp = (SeqPntPtr) location->data.ptrvalue;
      if (spp != NULL) {
        fuzz = spp->fuzz;
      }
    }

    CheckSeqLocForPartial (location, &noLeft, &noRight);
    hasNulls = LocationHasNullsBetween (location);
    loc = SeqLocMerge (bsp, location, NULL, FALSE, TRUE, hasNulls);
    if (loc == NULL) {
      tmp = TrimLocInSegment (bsp, location, &noLeft, &noRight);
      loc = SeqLocMerge (bsp, tmp, NULL, FALSE, TRUE, hasNulls);
      SeqLocFree (tmp);
    }
    if (loc == NULL) {
      return StringSave ("?");
    }
    FreeAllFuzz (loc);
    SetSeqLocPartial (loc, noLeft, noRight);

    if (loc->choice == SEQLOC_PNT && fuzz != NULL) {
      spp = (SeqPntPtr) loc->data.ptrvalue;
      if (spp != NULL && spp->fuzz == NULL) {
        spp->fuzz = AsnIoMemCopy ((Pointer) fuzz,
                                  (AsnReadFunc) IntFuzzAsnRead,
                                  (AsnWriteFunc) IntFuzzAsnWrite);
      }
    }

    FF_DoFlatLoc (ajp, ffstring, bsp, loc, TRUE);

    SeqLocFree (loc);

  } else {
    FF_DoFlatLoc (ajp, ffstring, bsp, location, TRUE);
  }

  str = FFToCharPtr(ffstring);
  FFRecycleString(ajp, ffstring);
  return str;
}


/******************************************************************************/
/*                            End FlatLoc functions.                          */
/******************************************************************************/



static void SubSourceToQualArray (
  SubSourcePtr ssp,
  QualValPtr qvp
)

{
  Int2   idx;
  Uint1  subtype;

  if (ssp == NULL || qvp == NULL) return;

  while (ssp != NULL) {
    subtype = ssp->subtype;
    if (subtype == 255) {
      subtype = 29;
    }
    if (subtype < 30) {
      idx = subSourceToSourceIdx [subtype];
      if (idx > 0 && idx < ASN2GNBK_TOTAL_SOURCE) {
        if (qvp [idx].ssp == NULL) {
          qvp [idx].ssp = ssp;
        }
      }
    }
    ssp = ssp->next;
  }
}

static Int2 orgModToSourceIdx [38] = {
  SCQUAL_zero_orgmod,
  SCQUAL_one_orgmod,
  SCQUAL_strain,
  SCQUAL_sub_strain,
  SCQUAL_type,
  SCQUAL_sub_type,
  SCQUAL_variety,
  SCQUAL_serotype,
  SCQUAL_serogroup,
  SCQUAL_serovar,
  SCQUAL_cultivar,
  SCQUAL_pathovar,
  SCQUAL_chemovar,
  SCQUAL_biovar,
  SCQUAL_biotype,
  SCQUAL_group,
  SCQUAL_sub_group,
  SCQUAL_isolate,
  SCQUAL_common,
  SCQUAL_acronym,
  SCQUAL_dosage,
  SCQUAL_spec_or_nat_host,
  SCQUAL_sub_species,
  SCQUAL_specimen_voucher,
  SCQUAL_authority,
  SCQUAL_forma,
  SCQUAL_forma_specialis,
  SCQUAL_ecotype,
  SCQUAL_synonym,
  SCQUAL_anamorph,
  SCQUAL_teleomorph,
  SCQUAL_breed,
  SCQUAL_gb_acronym,
  SCQUAL_gb_anamorph,
  SCQUAL_gb_synonym,
  SCQUAL_old_lineage,
  SCQUAL_old_name,
  SCQUAL_orgmod_note
};

static void OrgModToQualArray (
  OrgModPtr omp,
  QualValPtr qvp
)

{
  Int2   idx;
  Uint1  subtype;

  if (omp == NULL || qvp == NULL) return;

  while (omp != NULL) {
    subtype = omp->subtype;
    if (subtype == 253) {
      subtype = 35;
    } else if (subtype == 254) {
      subtype = 36;
    } else if (subtype == 255) {
      subtype = 37;
    }
    if (subtype < 38) {
      idx = orgModToSourceIdx [subtype];
      if (idx > 0 && idx < ASN2GNBK_TOTAL_SOURCE) {
        if (qvp [idx].omp == NULL) {
          qvp [idx].omp = omp;
        }
      }
    }
    omp = omp->next;
  }
}

static CharPtr organelleQual [] = {
  NULL,
  NULL,
  "/organelle=\"plastid:chloroplast\"",
  "/organelle=\"plastid:chromoplast\"",
  "/organelle=\"mitochondrion:kinetoplast\"",
  "/organelle=\"mitochondrion\"",
  "/organelle=\"plastid\"",
  "/macronuclear",
  NULL,
  "/plasmid=\"\"",
  "/transposon=\"\"",
  "/insertion_seq=\"\"",
  "/organelle=\"plastid:cyanelle\"",
  "/proviral",
  "/virion",
  "/organelle=\"nucleomorph\"",
  "/organelle=\"plastid:apicoplast\"",
  "/organelle=\"plastid:leucoplast\"",
  "/organelle=\"plastid:proplastid\"",
  NULL
};

static Boolean StringIsJustQuotes (
  CharPtr str
)

{
  Nlm_Uchar  ch;    /* to use 8bit characters in multibyte languages */

  if (str != NULL) {
    ch = *str;
    while (ch != '\0') {
      if (ch > ' ' && ch != '"' && ch != '\'') {
        return FALSE;
      }
      str++;
      ch = *str;
    }
  }
  return TRUE;
}

static CharPtr CleanQualValue (
  CharPtr str
)

{
  Char     ch;
  CharPtr  dst;
  CharPtr  ptr;

  if (str == NULL || str [0] == '\0') return NULL;

  dst = str;
  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    if (ch == '\n' || ch == '\r' || ch == '\t' || ch == '"') {
      *dst = ' ';
      dst++;
    } else {
      *dst = ch;
      dst++;
    }
    ptr++;
    ch = *ptr;
  }
  *dst = '\0';

  return str;
}

static CharPtr RemoveAllSpaces (
  CharPtr str
)

{
  Char     ch;
  CharPtr  dst;
  CharPtr  ptr;

  if (str == NULL || str [0] == '\0') return NULL;

  dst = str;
  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    if (ch != ' ') {
      *dst = ch;
      dst++;
    }
    ptr++;
    ch = *ptr;
  }
  *dst = '\0';

  return str;
}

static void AddFeatureToGbseq (
  GBSeqPtr gbseq,
  GBFeaturePtr gbfeat,
  CharPtr str
)

{
  Char            ch;
  CharPtr         copy;
  GBQualifierPtr  gbqual;
  GBQualifierPtr  last = NULL;
  CharPtr         ptr;
  CharPtr         qual;
  CharPtr         val;

  if (gbseq == NULL || gbfeat == NULL || StringHasNoText (str)) return;

  copy = StringSave (str);

  /* link in reverse order, to be reversed in slash block */

  gbfeat->next = gbseq->feature_table;
  gbseq->feature_table = gbfeat;

  /* now parse qualifiers */

  ptr = StringStr (copy, "                     /");
  while (ptr != NULL) {
    qual = ptr + 22;
    val = qual;
    ch = *val;
    while (ch != '=' && ch != '\n' && ch != '\0') {
      val++;
      ch = *val;
    }
    /*
    val = StringChr (qual, '=');
    if (val == NULL) {
      val = StringChr (qual, '\n');
    }
    */
    if (ch != '\0' /* val != NULL */) {
      *val = '\0';
      val++;
      if (ch == '=') {
        if (*val == '"') {
          val++;
        }
        ptr = StringStr (val, "\n                     /");
        if (ptr != NULL) {
          *ptr = '\0';
          ptr++;
        }
      } else {
        ptr = StringStr (val, "                     /");
        val = NULL;
      }
      gbqual = GBQualifierNew ();
      if (gbqual != NULL) {
        gbqual->name = StringSave (qual);
        if (! StringHasNoText (val)) {
          gbqual->value = StringSave (val);
          CleanQualValue (gbqual->value);
          Asn2gnbkCompressSpaces (gbqual->value);
          if (StringICmp (gbfeat->key, "CDS") == 0 &&
              StringICmp (qual, "translation") == 0) {
            RemoveAllSpaces (gbqual->value);
          }
        }
      }
    } else {
      gbqual = GBQualifierNew ();
      if (gbqual != NULL) {
        gbqual->name = StringSave (qual);
      }
    }
    if (gbfeat->quals == NULL) {
      gbfeat->quals = gbqual;
    } else if (last != NULL) {
      last->next = gbqual;
    }
    last = gbqual;
  }

  MemFree (copy);
}

static CharPtr FormatSourceFeatBlock (
  Asn2gbFormatPtr afp,
  BaseBlockPtr bbp
)

{
  Boolean            add_period;
  IntAsn2gbJobPtr    ajp;
  Asn2gbSectPtr      asp;
  BioSourcePtr       biop = NULL;
  BioseqPtr          bsp;
  Char               buf [80];
  CharPtr            common = NULL;
  DbtagPtr           dbt;
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext;
  GBFeaturePtr       gbfeat = NULL;
  GBSeqPtr           gbseq;
  Int2               i;
  Uint1              idx;
  IntSrcBlockPtr     isp;
  Boolean            is_desc = TRUE;
  Int2               j;
  Uint1              jdx;
  Uint1              lastomptype;
  Uint1              lastssptype;
  SeqLocPtr          location = NULL;
  CharPtr            notestr;
  Uint1Ptr           notetbl = NULL;
  Boolean            okay;
  ObjectIdPtr        oip;
  OrgModPtr          omp;
  OrgNamePtr         onp = NULL;
  OrgRefPtr          orp = NULL;
  CharPtr            prefix;
  Uint1Ptr           qualtbl = NULL;
  QualValPtr         qvp;
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp = NULL;
  SubSourcePtr       ssp;
  CharPtr            str;
  CharPtr            taxname = NULL;
  ValNodePtr         vnp;
  StringItemPtr      ffstring, unique;

  if (afp == NULL || bbp == NULL) return NULL;
  ajp = afp->ajp;
  if (ajp == NULL) return NULL;
  asp = afp->asp;
  if (asp == NULL) return NULL;
  bsp = asp->bsp;
  if (bsp == NULL) return NULL;
  qvp = afp->qvp;
  if (qvp == NULL) return NULL;

  if (ajp->gbseq) {
    gbseq = &asp->gbseq;
  } else {
    gbseq = NULL;
  }

  if (! StringHasNoText (bbp->string)) return StringSave (bbp->string);

  isp = (IntSrcBlockPtr) bbp;

  /* could be descriptor or feature */

  if (bbp->itemtype == OBJ_SEQDESC) {
    sdp = SeqMgrGetDesiredDescriptor (bbp->entityID, NULL, bbp->itemID, 0, NULL, &dcontext);
    if (sdp != NULL && dcontext.seqdesctype == Seq_descr_source) {
      biop = (BioSourcePtr) sdp->data.ptrvalue;
    }
  } else if (bbp->itemtype == OBJ_SEQFEAT) {
    sfp = SeqMgrGetDesiredFeature (bbp->entityID, NULL, bbp->itemID, 0, NULL, &fcontext);
    if (sfp != NULL && fcontext.seqfeattype == SEQFEAT_BIOSRC) {
      biop = (BioSourcePtr) sfp->data.value.ptrvalue;
    }
    is_desc = FALSE;
  }

  if (biop == NULL) return NULL;

  unique = FFGetString(ajp);
  if ( unique == NULL ) return NULL;

  ffstring = FFGetString(ajp);
  if ( ffstring == NULL ) return NULL;

  FFStartPrint (ffstring, afp->format, 5, 21, NULL, 0, 5, 21, "FT", FALSE);
  FFAddOneString (ffstring, "source", FALSE, FALSE, TILDE_IGNORE);
  FFAddNChar(ffstring, ' ', 21 - 5 - StringLen("source"), FALSE);

  if (gbseq != NULL) {
    gbfeat = GBFeatureNew ();
    if (gbfeat != NULL) {
      gbfeat->key = StringSave ("source");
    }
  }

  location = isp->loc;

  str = FlatLoc (ajp, bsp, location, ajp->masterStyle);
  if ( GetWWW(ajp) ) {
    FF_www_featloc (ffstring, str);
  } else {
    FFAddOneString (ffstring, str, FALSE, FALSE, TILDE_IGNORE);
  }
  FFAddOneChar(ffstring, '\n', FALSE);

  if (gbseq != NULL) {
    if (gbfeat != NULL) {
      gbfeat->location = StringSave (str);
    }
  }

  MemFree (str);

  orp = biop->org;
  if (orp != NULL) {
    taxname = orp->taxname;
    /* common = orp->common; */
  }
  if (StringHasNoText (taxname)) {
    if (ajp->flags.needOrganismQual) {
      taxname = "unknown";
      common = orp->common;
#ifdef ASN2GNBK_PRINT_UNKNOWN_ORG
    } else {
      taxname = "unknown";
      common = orp->common;
#endif
    }
  }

  /* populate qualifier table from biosource fields */

  qvp [SCQUAL_organism].str = taxname;
  qvp [SCQUAL_common_name].str = common;

  if (biop->is_focus) {
    qvp [SCQUAL_focus].ble = TRUE;
  }

  SubSourceToQualArray (biop->subtype, qvp);

  if (orp != NULL) {
    onp = orp->orgname;
    if (onp != NULL) {
      OrgModToQualArray (onp->mod, qvp);
    }

    if (! is_desc) {
      qvp [SCQUAL_unstructured].vnp = orp->mod;
    }
    qvp [SCQUAL_db_xref].vnp = orp->db;
  }

  if (sfp != NULL) {
    qvp [SCQUAL_org_xref].vnp = sfp->dbxref;
  }

  /* organelle currently prints /mitochondrion, /virion, etc. */

  qvp [SCQUAL_organelle].num = biop->genome;

  /* some qualifiers are flags in genome and names in subsource, print once with name */

  if (qvp [SCQUAL_ins_seq_name].ssp != NULL &&
      qvp [SCQUAL_organelle].num == GENOME_insertion_seq) {
    qvp [SCQUAL_organelle].num = 0;
  }
  if (qvp [SCQUAL_plasmid_name].ssp != NULL &&
      qvp [SCQUAL_organelle].num == GENOME_plasmid) {
    qvp [SCQUAL_organelle].num = 0;
  }
  /* AF095904.1
  if (qvp [SCQUAL_plastid_name].ssp != NULL &&
      qvp [SCQUAL_organelle].num == GENOME_plastid) {
    qvp [SCQUAL_organelle].num = 0;
  }
  */
  if (qvp [SCQUAL_transposon_name].ssp != NULL &&
      qvp [SCQUAL_organelle].num == GENOME_transposon) {
    qvp [SCQUAL_organelle].num = 0;
  }

  if (sfp != NULL) {
    qvp [SCQUAL_seqfeat_note].str = sfp->comment;
  }

  /* now print qualifiers from table */

  qualtbl = source_qual_order;
  if (is_desc) {
    notetbl = source_desc_note_order;
  } else {
    notetbl = source_feat_note_order;
  }

  for (i = 0, idx = qualtbl [i]; idx != 0; i++, idx = qualtbl [i]) {

    lastomptype = 0;
    lastssptype = 0;
    switch (asn2gnbk_source_quals [idx].qualclass) {

      case Qual_class_ignore :
        break;

      case Qual_class_string :
        if (! StringHasNoText (qvp [idx].str)) {
          FFAddTextToString(ffstring, "/", asn2gnbk_source_quals [idx].name, "=", 
                            FALSE, FALSE, TILDE_IGNORE);
          FFAddTextToString(ffstring, "\"", qvp [idx].str, "\"", 
                            FALSE, FALSE, TILDE_TO_SPACES);
          FFAddOneChar(ffstring, '\n', FALSE);
        }
        break;

      case Qual_class_boolean :
        if (qvp [idx].ble) {
          FFAddTextToString(ffstring, "/", asn2gnbk_source_quals [idx].name, "\n",
                            FALSE, TRUE, TILDE_IGNORE);
        }
        break;

      case Qual_class_organelle :
        j = (Int2) qvp [idx].num;
        if (organelleQual [j] != NULL) {
          FFAddTextToString(ffstring, NULL, organelleQual[j], "\n",
                            FALSE, FALSE, TILDE_IGNORE);
        }
        break;

      case Qual_class_orgmod :
        omp = qvp [idx].omp;
        if (lastomptype == 0 && omp != NULL) {
          lastomptype = omp->subtype;
        }
        while (omp != NULL && omp->subtype == lastomptype) {
          if (StringIsJustQuotes (omp->subname)) {
            FFAddTextToString(ffstring, "/", asn2gnbk_source_quals [idx].name, "=\"\"\n",
                              FALSE, TRUE, TILDE_IGNORE);
          } else if (! StringHasNoText (omp->subname)) {
            FFAddTextToString(ffstring, "/", asn2gnbk_source_quals [idx].name, "=",
                              FALSE, TRUE, TILDE_IGNORE);
            FFAddTextToString(ffstring, "\"", omp->subname, "\"\n",
                              FALSE, TRUE, TILDE_TO_SPACES);
          }
          omp = omp->next;
        }
        break;

      case Qual_class_subsource :
        ssp = qvp [idx].ssp;
        if (lastssptype == 0 && ssp != NULL) {
          lastssptype = ssp->subtype;
        }
        while (ssp != NULL && ssp->subtype == lastssptype) {
          if (ssp->subtype == SUBSRC_germline ||
              ssp->subtype == SUBSRC_rearranged ||
              ssp->subtype == SUBSRC_transgenic ||
              ssp->subtype == SUBSRC_environmental_sample) {
            FFAddTextToString(ffstring, "/", asn2gnbk_source_quals [idx].name, "\n",
                              FALSE, TRUE, TILDE_TO_SPACES);
          } else if (StringIsJustQuotes (ssp->name)) {
            FFAddTextToString(ffstring, "/", asn2gnbk_source_quals [idx].name, "=\"\"\n",
                              FALSE, TRUE, TILDE_IGNORE);
          } else if (! StringHasNoText (ssp->name)) {
            FFAddTextToString(ffstring, "/", asn2gnbk_source_quals [idx].name, "=",
                              FALSE, TRUE, TILDE_IGNORE);
            FFAddTextToString(ffstring, "\"", ssp->name, "\"\n",
                              FALSE, TRUE, TILDE_TO_SPACES);
          }
          ssp = ssp->next;
        }
        break;

      case Qual_class_pubset :
        break;

      case Qual_class_quote :
        break;

      case Qual_class_noquote :
        break;

      case Qual_class_label :
        break;

      case Qual_class_db_xref :
        for (vnp = qvp [idx].vnp; vnp != NULL; vnp = vnp->next) {
          buf [0] = '\0';
          dbt = (DbtagPtr) vnp->data.ptrvalue;
          if (dbt != NULL && (! StringHasNoText (dbt->db))) {
            oip = dbt->tag;
            if (oip != NULL) {

              okay = TRUE;
              if (ajp->flags.dropBadDbxref) {
                /* if RELEASE_MODE, drop unknown dbtag */
                okay = FALSE;
                for (j = 0; legalDbXrefs [j] != NULL; j++) {
                  if (StringCmp (dbt->db, legalDbXrefs [j]) == 0) {
                    okay = TRUE;
                  }
                }
              }

              if (okay) {
                if (! StringHasNoText (oip->str)) {
                  if (StringLen (dbt->db) + StringLen (oip->str) < 80) {
                    sprintf (buf, "%s", oip->str);
                  }
                } else {
                  sprintf (buf, "%ld", (long) oip->id);
                }
              }
            }
          }
          if (! StringHasNoText (buf)) {
            FFAddOneString(ffstring, "/db_xref=\"", FALSE, FALSE, TILDE_IGNORE);
            FF_www_db_xref(ajp, ffstring, dbt->db, buf);
            FFAddOneString(ffstring, "\"\n", FALSE, FALSE, TILDE_IGNORE);
          }
        }
        break;

      case Qual_class_illegal :
        break;

      case Qual_class_note :
        if (! ajp->flags.srcQualsToNote) {

          /* in sequin_mode and dump_mode, all orgmods and subsources show up as separate /qualifiers */

          for (j = 0, jdx = notetbl [j]; jdx != 0; j++, jdx = notetbl [j]) {

            lastomptype = 0;
            lastssptype = 0;
            switch (asn2gnbk_source_quals [jdx].qualclass) {

              case Qual_class_orgmod :
                if (jdx == SCQUAL_orgmod_note) break;
                omp = qvp [jdx].omp;
                if (lastomptype == 0 && omp != NULL) {
                  lastomptype = omp->subtype;
                }
                while (omp != NULL && omp->subtype == lastomptype) {
                  if (StringIsJustQuotes (omp->subname)) {
                    FFAddTextToString(ffstring, "/", asn2gnbk_source_quals [jdx].name, "=\"\"\n",
                              FALSE, TRUE, TILDE_IGNORE);
                  } else if (! StringHasNoText (omp->subname)) {
                    FFAddTextToString(ffstring, "/", asn2gnbk_source_quals [jdx].name, "=",
                                      FALSE, TRUE, TILDE_IGNORE);
                    FFAddTextToString(ffstring, "\"", omp->subname, "\"\n",
                                      FALSE, TRUE, TILDE_TO_SPACES);
                  }
                  omp = omp->next;
                }
                break;

              case Qual_class_subsource :
                if (jdx == SCQUAL_subsource_note) break;
                ssp = qvp [jdx].ssp;
                if (lastssptype == 0 && ssp != NULL) {
                  lastssptype = ssp->subtype;
                }
                while (ssp != NULL && ssp->subtype == lastssptype) {
                  if (StringIsJustQuotes (ssp->name)) {
                    FFAddTextToString(ffstring, "/", asn2gnbk_source_quals [jdx].name, "=\"\"\n",
                                      FALSE, TRUE, TILDE_IGNORE);

                  } else if (! StringHasNoText (ssp->name)) {
                    FFAddTextToString(ffstring, "/", asn2gnbk_source_quals [jdx].name, "=",
                                      FALSE, TRUE, TILDE_IGNORE);
                    FFAddTextToString(ffstring, "\"", ssp->name, "\"\n",
                                      FALSE, TRUE, TILDE_TO_SPACES);
                  }
                  ssp = ssp->next;
                }
                break;

              default :
                break;
            }
          }
        }

        notestr = NULL;
        prefix = "";
        add_period = FALSE;

        if (biop->genome == 8) {
          FFAddTextToString(unique, "", "extrachromosomal", NULL, FALSE, FALSE, TILDE_IGNORE); 
          prefix = "\n";
        }

        for (j = 0, jdx = notetbl [j]; jdx != 0; j++, jdx = notetbl [j]) {

          lastomptype = 0;
          lastssptype = 0;
          switch (asn2gnbk_source_quals [jdx].qualclass) {

            case Qual_class_string :
              if (! StringHasNoText (qvp [jdx].str)) {
                FFAddString_NoRedund (unique, prefix, qvp [jdx].str, NULL);
                add_period = FALSE;
                prefix = "\n";
              }
              break;

            case Qual_class_orgmod :
              if ((! ajp->flags.srcQualsToNote) && jdx != SCQUAL_orgmod_note) break;
              omp = qvp [jdx].omp;
              if (lastomptype == 0 && omp != NULL) {
                lastomptype = omp->subtype;
              }
              while (omp != NULL && omp->subtype == lastomptype) {
                if (! StringHasNoText (omp->subname)) {
                  if (jdx == SCQUAL_orgmod_note) {
                    sprintf (buf, "%s", prefix);
                  } else {
                    sprintf (buf, "%s%s: ", prefix, asn2gnbk_source_quals [jdx].name);
                  }

                  str = StringSave (omp->subname);
                  add_period = s_RemovePeriodFromEnd (str);
                  if (jdx == SCQUAL_orgmod_note) {
                    FFAddString_NoRedund (unique, buf, str, NULL);
                  } else {
                    FFAddTextToString(unique, buf, str, NULL, FALSE, FALSE, TILDE_IGNORE);
                  }
                  MemFree (str);

                  if (jdx == SCQUAL_orgmod_note) {
                    if (add_period) {
                      prefix = ".\n";
                    } else {
                      prefix = "\n";
                    }
                  } else {
                    prefix = "; ";
                  }
                }
                omp = omp->next;
              }
              break;

            case Qual_class_subsource :
              if ((! ajp->flags.srcQualsToNote) && jdx != SCQUAL_subsource_note) break;
              ssp = qvp [jdx].ssp;
              if (lastssptype == 0 && ssp != NULL) {
                lastssptype = ssp->subtype;
              }
              while (ssp != NULL && ssp->subtype == lastssptype) {
                if (! StringHasNoText (ssp->name)) {
                  if (jdx == SCQUAL_subsource_note) {
                    sprintf (buf, "%s", prefix);
                  } else {
                    sprintf (buf, "%s%s: ", prefix, asn2gnbk_source_quals [jdx].name);
                  }

                  str = StringSave (ssp->name);
                  add_period = s_RemovePeriodFromEnd (str);
                  if (jdx == SCQUAL_subsource_note) {
                    FFAddString_NoRedund (unique, buf, str, NULL);
                  } else {
                    FFAddTextToString(unique, buf, str, NULL, FALSE, FALSE, TILDE_IGNORE);
                  }
                  MemFree (str);

                  if (jdx == SCQUAL_subsource_note) {
                    if (add_period) {
                      prefix = ".\n";
                    } else {
                      prefix = "\n";
                    }
                  } else {
                    prefix = "; ";
                 }
                }
                ssp = ssp->next;
              }
              break;

            case Qual_class_valnode :
              for (vnp = qvp [jdx].vnp; vnp != NULL; vnp = vnp->next) {
                str = (CharPtr) vnp->data.ptrvalue;
                if (! StringHasNoText (str)) {
                  FFAddString_NoRedund (unique, prefix, str, NULL);
                  add_period = FALSE;
                  prefix = "; ";
                }
              }
              break;

            default :
              break;
          }
        }
        if ( !FFEmpty(unique) ) {
          notestr = FFToCharPtr(unique);
        
          if (add_period) {
            s_AddPeriodToEnd (notestr);
          }

#ifdef ASN2GNBK_STRIP_NOTE_PERIODS
          if (! IsEllipsis (notestr))
            s_RemovePeriodFromEnd (notestr);
#endif

          FFAddOneString (ffstring, "/note=\"", FALSE, FALSE, TILDE_IGNORE);
          if (is_desc) {
            /* AB055064.1 says TILDE_IGNORE on descriptors */
            FFAddOneString (ffstring, notestr, FALSE, TRUE, TILDE_IGNORE);
          } else {
            /* ASZ93724.1 says TILDE_EXPAND on features */
            FFAddOneString (ffstring, notestr, FALSE, TRUE, TILDE_EXPAND);
          }
          FFAddOneString (ffstring, "\"", FALSE, FALSE, TILDE_IGNORE);

          MemFree (notestr);
        }
        break;
      default :
        break;
    }
  }

  /* and then deal with the various note types separately (not in order table) */

  str = FFEndPrint(ajp, ffstring, afp->format, 21, 21, 5, 21, "FT"); 

  /* optionally populate gbseq for XML-ized GenBank format */

  if (gbseq != NULL) {
    if (gbfeat != NULL) {
      AddFeatureToGbseq (gbseq, gbfeat, str);
    }
  }

  FFRecycleString(ajp, unique);
  FFRecycleString(ajp, ffstring);
  return str;
}


typedef struct qualfeatur {
  CharPtr     name;
  Uint1       featurclass;
} QualFeatur, PNTR QualFeaturPtr;

#define NUM_GB_QUALS 25

static QualFeatur qualToFeature [NUM_GB_QUALS] = {
  { "allele",         FTQUAL_allele         },
  { "bound_moiety",   FTQUAL_bound_moiety   },
  { "clone",          FTQUAL_clone          },
  { "codon",          FTQUAL_codon          },
  { "cons_splice",    FTQUAL_cons_splice    },
  { "direction",      FTQUAL_direction      },
  { "EC_number",      FTQUAL_EC_number      },
  { "frequency",      FTQUAL_frequency      },
  { "function",       FTQUAL_function       },
  { "insertion_seq",  FTQUAL_insertion_seq  },
  { "label",          FTQUAL_label          },
  { "map",            FTQUAL_map            },
  { "mod_base",       FTQUAL_mod_base       },
  { "number",         FTQUAL_number         },
  { "organism",       FTQUAL_organism       },
  { "PCR_conditions", FTQUAL_PCR_conditions },
  { "phenotype",      FTQUAL_phenotype      },
  { "product",        FTQUAL_product_quals  },
  { "replace",        FTQUAL_replace        },
  { "rpt_family",     FTQUAL_rpt_family     },
  { "rpt_type",       FTQUAL_rpt_type       },
  { "rpt_unit",       FTQUAL_rpt_unit       },
  { "standard_name",  FTQUAL_standard_name  },
  { "transposon",     FTQUAL_transposon     },
  { "usedin",         FTQUAL_usedin         }
};

static Int2 GbqualToFeaturIndex (
  CharPtr qualname
)

{
  Int2  L, R, mid;

  if (qualname == NULL || *qualname == '\0') return 0;

  L = 0;
  R = NUM_GB_QUALS - 1;

  while (L < R) {
    mid = (L + R) / 2;
    if (StringICmp (qualToFeature [mid].name, qualname) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  if (StringICmp (qualToFeature [R].name, qualname) == 0) {
    return qualToFeature [R].featurclass;
  }

  return 0;
}

#define NUM_ILLEGAL_QUALS 14

static FeaturQual illegalGbqualList [NUM_ILLEGAL_QUALS] = {
  { "anticodon",      Qual_class_noquote },
  { "citation",       Qual_class_noquote },
  { "codon_start",    Qual_class_noquote },
  { "db_xref",        Qual_class_quote   },
  { "evidence",       Qual_class_noquote },
  { "exception",      Qual_class_quote   },
  { "gene",           Qual_class_quote   },
  { "note",           Qual_class_quote   },
  { "protein_id",     Qual_class_quote   },
  { "pseudo",         Qual_class_noquote },
  { "transcript_id",  Qual_class_quote   },
  { "transl_except",  Qual_class_noquote },
  { "transl_table",   Qual_class_noquote },
  { "translation",    Qual_class_quote   },
};

static Int2 IllegalGbqualToClass (
  CharPtr qualname
)

{
  Int2  L, R, mid;

  if (qualname == NULL || *qualname == '\0') return 0;

  L = 0;
  R = NUM_ILLEGAL_QUALS - 1;

  while (L < R) {
    mid = (L + R) / 2;
    if (StringICmp (illegalGbqualList [mid].name, qualname) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  if (StringICmp (illegalGbqualList [R].name, qualname) == 0) {
    return illegalGbqualList [R].qualclass;
  }

  return 0;
}

static CharPtr trnaList [] = {
  "tRNA-Gap",
  "tRNA-Ala",
  "tRNA-Asx",
  "tRNA-Cys",
  "tRNA-Asp",
  "tRNA-Glu",
  "tRNA-Phe",
  "tRNA-Gly",
  "tRNA-His",
  "tRNA-Ile",
  "tRNA-Lys",
  "tRNA-Leu",
  "tRNA-Met",
  "tRNA-Asn",
  "tRNA-Pro",
  "tRNA-Gln",
  "tRNA-Arg",
  "tRNA-Ser",
  "tRNA-Thr",
  "tRNA-Sec",
  "tRNA-Val",
  "tRNA-Trp",
  "tRNA-OTHER",
  "tRNA-Tyr",
  "tRNA-Glx",
  NULL
};

static CharPtr evidenceText [] = {
  NULL, "experimental", "not_experimental"
};

static CharPtr secStrText [] = {
  NULL, "helix", "sheet", "turn"
};

static CharPtr  oops = "?";

static CharPtr SeqCodeNameGet (
  SeqCodeTablePtr table,
  Uint1 residue
)

{
  Uint1  index;

  if (table != NULL) {
    index = residue - table->start_at;
    if ( /*index >= 0 && */ index < table->num) {
      return (table->names) [index];
    }
  }

  return oops;
}

static CharPtr Get3LetterSymbol (
  IntAsn2gbJobPtr  ajp,
  Uint1 seq_code,
  SeqCodeTablePtr table,
  Uint1 residue
)

{
  Uint1            code = Seq_code_ncbieaa;
  Int2             index;
  Uint1            new_residue;
  CharPtr          ptr;
  CharPtr          retval = NULL;
  SeqMapTablePtr   smtp;
  SeqCodeTablePtr  table_3aa;

  if (residue == 42) { /* stop codon in NCBIeaa */
    retval = "TERM";
    return retval;
  }

  if (ajp->flags.iupacaaOnly) {
    code = Seq_code_iupacaa;
  } else {
    code = Seq_code_ncbieaa;
  }

  if (seq_code != code) {
    /* if code and seq_code are identical, then smtp is NULL?? */
    smtp = SeqMapTableFind (code, seq_code);
    new_residue = SeqMapTableConvert (smtp, residue);
  } else {
    new_residue = residue;
  }

  /* The following looks for non-symbols (255) and "Undetermined" (88) */
  if ((int) new_residue == 255 || (int) new_residue == 88) {
    retval = "OTHER";
    return retval;
  } else {
    ptr = SeqCodeNameGet (table, residue);
    table_3aa=SeqCodeTableFind  (Seq_code_iupacaa3);
    if (ptr != NULL && table_3aa != NULL) {
      for (index=0; index < (int) table_3aa->num; index++) {
        if (StringCmp(ptr, (table_3aa->names) [index]) == 0) {
          retval = (table_3aa->symbols) [index];
          return retval;
        }
      }
    }
  }

  retval = "OTHER";
  return retval;
}

static Boolean MatchCit (
  ValNodePtr ppr,
  RefBlockPtr rbp
)

{
  Char        buf [121];
  size_t      len;
  Int4        uid;
  ValNodePtr  vnp;

  if (ppr == NULL || rbp == NULL) return FALSE;
  switch (ppr->choice) {
    case PUB_Muid :
      uid = ppr->data.intvalue;
      if (rbp->muid == uid) return TRUE;
      break;
    case PUB_PMid :
      uid = ppr->data.intvalue;
      if (rbp->pmid == uid) return TRUE;
      break;
    case PUB_Equiv :
      for (vnp = (ValNodePtr) ppr->data.ptrvalue; vnp != NULL; vnp = vnp->next) {
        if (MatchCit (vnp, rbp)) return TRUE;
      }
      break;
    default :
      PubLabelUnique (ppr, buf, sizeof (buf) - 1, OM_LABEL_CONTENT, TRUE);
      len = StringLen (buf);
      if (len > 0 && buf [len - 1] == '>') {
        buf [len - 1] = '\0';
        len--;
      }
      len = MIN (len, StringLen (rbp->uniquestr));
      if (StringNICmp (rbp->uniquestr, buf, len) == 0) return TRUE;
      break;
  }
  return FALSE;
}

static Int2 MatchRef (
  ValNodePtr ppr,
  RefBlockPtr PNTR rbpp,
  Int2 numReferences
)

{
  Int2         j;
  RefBlockPtr  rbp;

  if (ppr == NULL || rbpp == NULL) return 0;

  for (j = 0; j < numReferences; j++) {
    rbp = rbpp [j];
    if (rbp == NULL) continue;
    if (MatchCit (ppr, rbp)) return j + 1;
  }
  return 0;
}

/* taken from asn2ff4.c */

static Boolean LookForFuzz (SeqLocPtr head)
{
  Boolean retval=FALSE;
  IntFuzzPtr ifp;
  PackSeqPntPtr pspp;
  SeqIntPtr sip;
  SeqLocPtr slp;
  SeqPntPtr spp;

  if (head == NULL)
    return retval;

  slp=NULL;
  while ((slp = SeqLocFindNext(head, slp)) != NULL)
  {
    switch (slp->choice)
    {
      case SEQLOC_INT:
        sip = (SeqIntPtr)(slp->data.ptrvalue);
        ifp = sip->if_from;
        if (ifp != NULL)
        {
          if (ifp->choice == 4)
          {
            if (ifp->a != 0)
              retval=TRUE;
          }
          else
            retval = TRUE;
        }
        ifp = sip->if_to;
        if (ifp != NULL)
        {
          if (ifp->choice == 4)
          {
            if (ifp->a != 0)
              retval=TRUE;
          }
          else
            retval = TRUE;
        }
        break;
      case SEQLOC_PNT:
        spp = (SeqPntPtr)(slp->data.ptrvalue);
        ifp = spp->fuzz;
        if (ifp != NULL)
        {
          if (ifp->choice == 4)
          {
            if (ifp->a != 0)
              retval=TRUE;
          }
          else
            retval = TRUE;
        }
        break;
      case SEQLOC_PACKED_PNT:
        pspp = (PackSeqPntPtr)(slp->data.ptrvalue);
        ifp = pspp->fuzz;
        if (ifp != NULL)
        {
          if (ifp->choice == 4)
          {
            if (ifp->a != 0)
              retval=TRUE;
          }
          else
            retval = TRUE;
        }
        break;
      default:
        break;
    }
    if (retval == TRUE)
      break;
  }
  return retval;
}

static CharPtr bondList [] = {
  NULL,
  "disulfide",
  "thiolester",
  "xlink",
  "thioether",
  "unclassified"
};

static CharPtr siteList [] = {
  NULL,
  "active",
  "binding",
  "cleavage",
  "inhibit",
  "modified",
  "glycosylation",
  "myristoylation",
  "mutagenized",
  "metal-binding",
  "phosphorylation",
  "acetylation",
  "amidation",
  "methylation",
  "hydroxylation",
  "sulfatation",
  "oxidative-deamination",
  "pyrrolidone-carboxylic-acid",
  "gamma-carboxyglutamic-acid",
  "blocked",
  "lipid-binding",
  "np-binding",
  "DNA binding",
  "signal-peptide",
  "transit-peptide",
  "transmembrane-region",
  "unclassified"
};

static CharPtr conflict_msg =
"Protein sequence is in conflict with the conceptual translation";

static CharPtr no_protein_msg =
"Coding region translates with internal stops";

/**/
/*  s_DisplayQVP () -- Displays the strings in a QVP structure.   */
/*                     This is a debugging function only.         */
/**/

static void s_DisplayQVP(QualValPtr qvp, Uint1Ptr notetbl)
{
  Int2 j;
  Int2 jdx;

  fprintf(stderr,"\n");
  for (j = 0, jdx = notetbl [j]; jdx != 0; j++, jdx = notetbl [j])
    {
      if (((int) qvp[jdx].str != 0x1000000) &&
      ((int) qvp[jdx].str != 0x1) &&
      ((int) qvp[jdx].str != 0xb) &&
      (qvp[jdx].str != NULL))
    fprintf(stderr, "%d\t%-25s %s\n", j, asn2gnbk_featur_quals[jdx].name,
        qvp[jdx].str);
      else
    fprintf(stderr, "%d\t%-25s %s\n", j, asn2gnbk_featur_quals[jdx].name,
        "NULL");
    }
}

static Boolean NotInGeneSyn (
  CharPtr str,
  ValNodePtr gene_syn)

{
  CharPtr     syn;
  ValNodePtr  vnp;

  for (vnp = gene_syn; vnp != NULL; vnp = vnp->next) {
    syn = (CharPtr) vnp->data.ptrvalue;
    if (! StringHasNoText (syn)) {
      if (StringICmp (str, syn) == 0) return FALSE;
    }
  }
  return TRUE;
}

typedef struct valqualstruc {
  Uint2       featdef;
  FtQualType  ftqual;
} ValQual, PNTR ValQualPtr;

/*
   WARNING - This list MUST be kept sorted in FEATDEF order as the primary
   key, and within a FEATDEF group sorted by FTQUAL as the secondary key
*/

static ValQual legalGbqualList [] = {

  {FEATDEF_GENE ,  FTQUAL_allele},
  {FEATDEF_GENE ,  FTQUAL_function},
  {FEATDEF_GENE ,  FTQUAL_label},
  {FEATDEF_GENE ,  FTQUAL_map},
  {FEATDEF_GENE ,  FTQUAL_phenotype},
  {FEATDEF_GENE ,  FTQUAL_product},
  {FEATDEF_GENE ,  FTQUAL_standard_name},
  {FEATDEF_GENE ,  FTQUAL_usedin},

  {FEATDEF_CDS ,  FTQUAL_allele},
  {FEATDEF_CDS ,  FTQUAL_codon},
  {FEATDEF_CDS ,  FTQUAL_label},
  {FEATDEF_CDS ,  FTQUAL_map},
  {FEATDEF_CDS ,  FTQUAL_number},
  {FEATDEF_CDS ,  FTQUAL_standard_name},
  {FEATDEF_CDS ,  FTQUAL_usedin},

  {FEATDEF_PROT ,  FTQUAL_product},

  {FEATDEF_preRNA ,  FTQUAL_allele},
  {FEATDEF_preRNA ,  FTQUAL_function},
  {FEATDEF_preRNA ,  FTQUAL_label},
  {FEATDEF_preRNA ,  FTQUAL_map},
  {FEATDEF_preRNA ,  FTQUAL_product},
  {FEATDEF_preRNA ,  FTQUAL_standard_name},
  {FEATDEF_preRNA ,  FTQUAL_usedin},

  {FEATDEF_mRNA ,  FTQUAL_allele},
  {FEATDEF_mRNA ,  FTQUAL_function},
  {FEATDEF_mRNA ,  FTQUAL_label},
  {FEATDEF_mRNA ,  FTQUAL_map},
  {FEATDEF_mRNA ,  FTQUAL_product},
  {FEATDEF_mRNA ,  FTQUAL_standard_name},
  {FEATDEF_mRNA ,  FTQUAL_usedin},

  {FEATDEF_tRNA ,  FTQUAL_function},
  {FEATDEF_tRNA ,  FTQUAL_label},
  {FEATDEF_tRNA ,  FTQUAL_map},
  {FEATDEF_tRNA ,  FTQUAL_product},
  {FEATDEF_tRNA ,  FTQUAL_standard_name},
  {FEATDEF_tRNA ,  FTQUAL_usedin},

  {FEATDEF_rRNA ,  FTQUAL_function},
  {FEATDEF_rRNA ,  FTQUAL_label},
  {FEATDEF_rRNA ,  FTQUAL_map},
  {FEATDEF_rRNA ,  FTQUAL_product},
  {FEATDEF_rRNA ,  FTQUAL_standard_name},
  {FEATDEF_rRNA ,  FTQUAL_usedin},

  {FEATDEF_snRNA ,  FTQUAL_function},
  {FEATDEF_snRNA ,  FTQUAL_label},
  {FEATDEF_snRNA ,  FTQUAL_map},
  {FEATDEF_snRNA ,  FTQUAL_product},
  {FEATDEF_snRNA ,  FTQUAL_standard_name},
  {FEATDEF_snRNA ,  FTQUAL_usedin},

  {FEATDEF_scRNA ,  FTQUAL_function},
  {FEATDEF_scRNA ,  FTQUAL_label},
  {FEATDEF_scRNA ,  FTQUAL_map},
  {FEATDEF_scRNA ,  FTQUAL_product},
  {FEATDEF_scRNA ,  FTQUAL_standard_name},
  {FEATDEF_scRNA ,  FTQUAL_usedin},

  {FEATDEF_otherRNA ,  FTQUAL_function},
  {FEATDEF_otherRNA ,  FTQUAL_label},
  {FEATDEF_otherRNA ,  FTQUAL_map},
  {FEATDEF_otherRNA ,  FTQUAL_product},
  {FEATDEF_otherRNA ,  FTQUAL_standard_name},
  {FEATDEF_otherRNA ,  FTQUAL_usedin},

  {FEATDEF_attenuator ,  FTQUAL_label},
  {FEATDEF_attenuator ,  FTQUAL_map},
  {FEATDEF_attenuator ,  FTQUAL_phenotype},
  {FEATDEF_attenuator ,  FTQUAL_usedin},

  {FEATDEF_C_region ,  FTQUAL_label},
  {FEATDEF_C_region ,  FTQUAL_map},
  {FEATDEF_C_region ,  FTQUAL_product},
  {FEATDEF_C_region ,  FTQUAL_standard_name},
  {FEATDEF_C_region ,  FTQUAL_usedin},

  {FEATDEF_CAAT_signal ,  FTQUAL_label},
  {FEATDEF_CAAT_signal ,  FTQUAL_map},
  {FEATDEF_CAAT_signal ,  FTQUAL_usedin},

  {FEATDEF_Imp_CDS ,  FTQUAL_codon},
  {FEATDEF_Imp_CDS ,  FTQUAL_EC_number},
  {FEATDEF_Imp_CDS ,  FTQUAL_function},
  {FEATDEF_Imp_CDS ,  FTQUAL_label},
  {FEATDEF_Imp_CDS ,  FTQUAL_map},
  {FEATDEF_Imp_CDS ,  FTQUAL_number},
  {FEATDEF_Imp_CDS ,  FTQUAL_product},
  {FEATDEF_Imp_CDS ,  FTQUAL_standard_name},
  {FEATDEF_Imp_CDS ,  FTQUAL_usedin},

  {FEATDEF_conflict ,  FTQUAL_label},
  {FEATDEF_conflict ,  FTQUAL_map},
  {FEATDEF_conflict ,  FTQUAL_replace},
  {FEATDEF_conflict ,  FTQUAL_usedin},

  {FEATDEF_D_loop ,  FTQUAL_label},
  {FEATDEF_D_loop ,  FTQUAL_map},
  {FEATDEF_D_loop ,  FTQUAL_usedin},

  {FEATDEF_D_segment ,  FTQUAL_label},
  {FEATDEF_D_segment ,  FTQUAL_map},
  {FEATDEF_D_segment ,  FTQUAL_product},
  {FEATDEF_D_segment ,  FTQUAL_standard_name},
  {FEATDEF_D_segment ,  FTQUAL_usedin},

  {FEATDEF_enhancer ,  FTQUAL_label},
  {FEATDEF_enhancer ,  FTQUAL_map},
  {FEATDEF_enhancer ,  FTQUAL_standard_name},
  {FEATDEF_enhancer ,  FTQUAL_usedin},

  {FEATDEF_exon ,  FTQUAL_allele},
  {FEATDEF_exon ,  FTQUAL_EC_number},
  {FEATDEF_exon ,  FTQUAL_function},
  {FEATDEF_exon ,  FTQUAL_label},
  {FEATDEF_exon ,  FTQUAL_map},
  {FEATDEF_exon ,  FTQUAL_number},
  {FEATDEF_exon ,  FTQUAL_product},
  {FEATDEF_exon ,  FTQUAL_standard_name},
  {FEATDEF_exon ,  FTQUAL_usedin},

  {FEATDEF_GC_signal ,  FTQUAL_label},
  {FEATDEF_GC_signal ,  FTQUAL_map},
  {FEATDEF_GC_signal ,  FTQUAL_usedin},

  {FEATDEF_iDNA ,  FTQUAL_function},
  {FEATDEF_iDNA ,  FTQUAL_label},
  {FEATDEF_iDNA ,  FTQUAL_map},
  {FEATDEF_iDNA ,  FTQUAL_number},
  {FEATDEF_iDNA ,  FTQUAL_standard_name},
  {FEATDEF_iDNA ,  FTQUAL_usedin},

  {FEATDEF_intron ,  FTQUAL_allele},
  {FEATDEF_intron ,  FTQUAL_cons_splice},
  {FEATDEF_intron ,  FTQUAL_function},
  {FEATDEF_intron ,  FTQUAL_label},
  {FEATDEF_intron ,  FTQUAL_map},
  {FEATDEF_intron ,  FTQUAL_number},
  {FEATDEF_intron ,  FTQUAL_standard_name},
  {FEATDEF_intron ,  FTQUAL_usedin},

  {FEATDEF_J_segment ,  FTQUAL_label},
  {FEATDEF_J_segment ,  FTQUAL_map},
  {FEATDEF_J_segment ,  FTQUAL_product},
  {FEATDEF_J_segment ,  FTQUAL_standard_name},
  {FEATDEF_J_segment ,  FTQUAL_usedin},

  {FEATDEF_LTR ,  FTQUAL_function},
  {FEATDEF_LTR ,  FTQUAL_label},
  {FEATDEF_LTR ,  FTQUAL_map},
  {FEATDEF_LTR ,  FTQUAL_standard_name},
  {FEATDEF_LTR ,  FTQUAL_usedin},

  {FEATDEF_mat_peptide ,  FTQUAL_EC_number},
  {FEATDEF_mat_peptide ,  FTQUAL_function},
  {FEATDEF_mat_peptide ,  FTQUAL_label},
  {FEATDEF_mat_peptide ,  FTQUAL_map},
  {FEATDEF_mat_peptide ,  FTQUAL_product},
  {FEATDEF_mat_peptide ,  FTQUAL_standard_name},
  {FEATDEF_mat_peptide ,  FTQUAL_usedin},

  {FEATDEF_misc_binding ,  FTQUAL_bound_moiety},
  {FEATDEF_misc_binding ,  FTQUAL_function},
  {FEATDEF_misc_binding ,  FTQUAL_label},
  {FEATDEF_misc_binding ,  FTQUAL_map},
  {FEATDEF_misc_binding ,  FTQUAL_usedin},

  {FEATDEF_misc_difference ,  FTQUAL_clone},
  {FEATDEF_misc_difference ,  FTQUAL_label},
  {FEATDEF_misc_difference ,  FTQUAL_map},
  {FEATDEF_misc_difference ,  FTQUAL_phenotype},
  {FEATDEF_misc_difference ,  FTQUAL_replace},
  {FEATDEF_misc_difference ,  FTQUAL_standard_name},
  {FEATDEF_misc_difference ,  FTQUAL_usedin},

  {FEATDEF_misc_feature ,  FTQUAL_function},
  {FEATDEF_misc_feature ,  FTQUAL_label},
  {FEATDEF_misc_feature ,  FTQUAL_map},
  {FEATDEF_misc_feature ,  FTQUAL_number},
  {FEATDEF_misc_feature ,  FTQUAL_phenotype},
  {FEATDEF_misc_feature ,  FTQUAL_product},
  {FEATDEF_misc_feature ,  FTQUAL_standard_name},
  {FEATDEF_misc_feature ,  FTQUAL_usedin},

  {FEATDEF_misc_recomb ,  FTQUAL_label},
  {FEATDEF_misc_recomb ,  FTQUAL_map},
  {FEATDEF_misc_recomb ,  FTQUAL_organism},
  {FEATDEF_misc_recomb ,  FTQUAL_standard_name},
  {FEATDEF_misc_recomb ,  FTQUAL_usedin},

  {FEATDEF_misc_signal ,  FTQUAL_function},
  {FEATDEF_misc_signal ,  FTQUAL_label},
  {FEATDEF_misc_signal ,  FTQUAL_map},
  {FEATDEF_misc_signal ,  FTQUAL_phenotype},
  {FEATDEF_misc_signal ,  FTQUAL_standard_name},
  {FEATDEF_misc_signal ,  FTQUAL_usedin},

  {FEATDEF_misc_structure ,  FTQUAL_function},
  {FEATDEF_misc_structure ,  FTQUAL_label},
  {FEATDEF_misc_structure ,  FTQUAL_map},
  {FEATDEF_misc_structure ,  FTQUAL_standard_name},
  {FEATDEF_misc_structure ,  FTQUAL_usedin},

  {FEATDEF_modified_base ,  FTQUAL_frequency},
  {FEATDEF_modified_base ,  FTQUAL_label},
  {FEATDEF_modified_base ,  FTQUAL_map},
  {FEATDEF_modified_base ,  FTQUAL_mod_base},
  {FEATDEF_modified_base ,  FTQUAL_usedin},

  {FEATDEF_N_region ,  FTQUAL_label},
  {FEATDEF_N_region ,  FTQUAL_map},
  {FEATDEF_N_region ,  FTQUAL_product},
  {FEATDEF_N_region ,  FTQUAL_standard_name},
  {FEATDEF_N_region ,  FTQUAL_usedin},

  {FEATDEF_old_sequence ,  FTQUAL_label},
  {FEATDEF_old_sequence ,  FTQUAL_map},
  {FEATDEF_old_sequence ,  FTQUAL_replace},
  {FEATDEF_old_sequence ,  FTQUAL_usedin},

  {FEATDEF_polyA_signal ,  FTQUAL_label},
  {FEATDEF_polyA_signal ,  FTQUAL_map},
  {FEATDEF_polyA_signal ,  FTQUAL_usedin},

  {FEATDEF_polyA_site ,  FTQUAL_label},
  {FEATDEF_polyA_site ,  FTQUAL_map},
  {FEATDEF_polyA_site ,  FTQUAL_usedin},

  {FEATDEF_prim_transcript ,  FTQUAL_allele},
  {FEATDEF_prim_transcript ,  FTQUAL_function},
  {FEATDEF_prim_transcript ,  FTQUAL_label},
  {FEATDEF_prim_transcript ,  FTQUAL_map},
  {FEATDEF_prim_transcript ,  FTQUAL_standard_name},
  {FEATDEF_prim_transcript ,  FTQUAL_usedin},

  {FEATDEF_primer_bind ,  FTQUAL_label},
  {FEATDEF_primer_bind ,  FTQUAL_map},
  {FEATDEF_primer_bind ,  FTQUAL_PCR_conditions},
  {FEATDEF_primer_bind ,  FTQUAL_standard_name},
  {FEATDEF_primer_bind ,  FTQUAL_usedin},

  {FEATDEF_promoter ,  FTQUAL_function},
  {FEATDEF_promoter ,  FTQUAL_label},
  {FEATDEF_promoter ,  FTQUAL_map},
  {FEATDEF_promoter ,  FTQUAL_phenotype},
  {FEATDEF_promoter ,  FTQUAL_standard_name},
  {FEATDEF_promoter ,  FTQUAL_usedin},

  {FEATDEF_protein_bind ,  FTQUAL_bound_moiety},
  {FEATDEF_protein_bind ,  FTQUAL_function},
  {FEATDEF_protein_bind ,  FTQUAL_label},
  {FEATDEF_protein_bind ,  FTQUAL_map},
  {FEATDEF_protein_bind ,  FTQUAL_standard_name},
  {FEATDEF_protein_bind ,  FTQUAL_usedin},

  {FEATDEF_RBS ,  FTQUAL_label},
  {FEATDEF_RBS ,  FTQUAL_map},
  {FEATDEF_RBS ,  FTQUAL_standard_name},
  {FEATDEF_RBS ,  FTQUAL_usedin},

  {FEATDEF_repeat_region ,  FTQUAL_function},
  {FEATDEF_repeat_region ,  FTQUAL_insertion_seq},
  {FEATDEF_repeat_region ,  FTQUAL_label},
  {FEATDEF_repeat_region ,  FTQUAL_map},
  {FEATDEF_repeat_region ,  FTQUAL_rpt_family},
  {FEATDEF_repeat_region ,  FTQUAL_rpt_type},
  {FEATDEF_repeat_region ,  FTQUAL_rpt_unit},
  {FEATDEF_repeat_region ,  FTQUAL_standard_name},
  {FEATDEF_repeat_region ,  FTQUAL_transposon},
  {FEATDEF_repeat_region ,  FTQUAL_usedin},

  {FEATDEF_repeat_unit ,  FTQUAL_function},
  {FEATDEF_repeat_unit ,  FTQUAL_label},
  {FEATDEF_repeat_unit ,  FTQUAL_map},
  {FEATDEF_repeat_unit ,  FTQUAL_rpt_family},
  {FEATDEF_repeat_unit ,  FTQUAL_rpt_type},
  {FEATDEF_repeat_unit ,  FTQUAL_usedin},

  {FEATDEF_rep_origin ,  FTQUAL_direction},
  {FEATDEF_rep_origin ,  FTQUAL_label},
  {FEATDEF_rep_origin ,  FTQUAL_map},
  {FEATDEF_rep_origin ,  FTQUAL_standard_name},
  {FEATDEF_rep_origin ,  FTQUAL_usedin},

  {FEATDEF_S_region ,  FTQUAL_label},
  {FEATDEF_S_region ,  FTQUAL_map},
  {FEATDEF_S_region ,  FTQUAL_product},
  {FEATDEF_S_region ,  FTQUAL_standard_name},
  {FEATDEF_S_region ,  FTQUAL_usedin},

  {FEATDEF_satellite ,  FTQUAL_label},
  {FEATDEF_satellite ,  FTQUAL_map},
  {FEATDEF_satellite ,  FTQUAL_rpt_family},
  {FEATDEF_satellite ,  FTQUAL_rpt_type},
  {FEATDEF_satellite ,  FTQUAL_rpt_unit},
  {FEATDEF_satellite ,  FTQUAL_standard_name},
  {FEATDEF_satellite ,  FTQUAL_usedin},

  {FEATDEF_sig_peptide ,  FTQUAL_function},
  {FEATDEF_sig_peptide ,  FTQUAL_label},
  {FEATDEF_sig_peptide ,  FTQUAL_map},
  {FEATDEF_sig_peptide ,  FTQUAL_product},
  {FEATDEF_sig_peptide ,  FTQUAL_standard_name},
  {FEATDEF_sig_peptide ,  FTQUAL_usedin},

  {FEATDEF_stem_loop ,  FTQUAL_function},
  {FEATDEF_stem_loop ,  FTQUAL_label},
  {FEATDEF_stem_loop ,  FTQUAL_map},
  {FEATDEF_stem_loop ,  FTQUAL_standard_name},
  {FEATDEF_stem_loop ,  FTQUAL_usedin},

  {FEATDEF_STS ,  FTQUAL_label},
  {FEATDEF_STS ,  FTQUAL_map},
  {FEATDEF_STS ,  FTQUAL_standard_name},
  {FEATDEF_STS ,  FTQUAL_usedin},

  {FEATDEF_TATA_signal ,  FTQUAL_label},
  {FEATDEF_TATA_signal ,  FTQUAL_map},
  {FEATDEF_TATA_signal ,  FTQUAL_usedin},

  {FEATDEF_terminator ,  FTQUAL_label},
  {FEATDEF_terminator ,  FTQUAL_map},
  {FEATDEF_terminator ,  FTQUAL_standard_name},
  {FEATDEF_terminator ,  FTQUAL_usedin},

  {FEATDEF_transit_peptide ,  FTQUAL_function},
  {FEATDEF_transit_peptide ,  FTQUAL_label},
  {FEATDEF_transit_peptide ,  FTQUAL_map},
  {FEATDEF_transit_peptide ,  FTQUAL_product},
  {FEATDEF_transit_peptide ,  FTQUAL_standard_name},
  {FEATDEF_transit_peptide ,  FTQUAL_usedin},

  {FEATDEF_unsure ,  FTQUAL_label},
  {FEATDEF_unsure ,  FTQUAL_map},
  {FEATDEF_unsure ,  FTQUAL_replace},
  {FEATDEF_unsure ,  FTQUAL_usedin},

  {FEATDEF_V_region ,  FTQUAL_label},
  {FEATDEF_V_region ,  FTQUAL_map},
  {FEATDEF_V_region ,  FTQUAL_product},
  {FEATDEF_V_region ,  FTQUAL_standard_name},
  {FEATDEF_V_region ,  FTQUAL_usedin},

  {FEATDEF_V_segment ,  FTQUAL_label},
  {FEATDEF_V_segment ,  FTQUAL_map},
  {FEATDEF_V_segment ,  FTQUAL_product},
  {FEATDEF_V_segment ,  FTQUAL_standard_name},
  {FEATDEF_V_segment ,  FTQUAL_usedin},

  {FEATDEF_variation ,  FTQUAL_allele},
  {FEATDEF_variation ,  FTQUAL_frequency},
  {FEATDEF_variation ,  FTQUAL_label},
  {FEATDEF_variation ,  FTQUAL_map},
  {FEATDEF_variation ,  FTQUAL_phenotype},
  {FEATDEF_variation ,  FTQUAL_product},
  {FEATDEF_variation ,  FTQUAL_replace},
  {FEATDEF_variation ,  FTQUAL_standard_name},
  {FEATDEF_variation ,  FTQUAL_usedin},

  {FEATDEF_3clip ,  FTQUAL_allele},
  {FEATDEF_3clip ,  FTQUAL_function},
  {FEATDEF_3clip ,  FTQUAL_label},
  {FEATDEF_3clip ,  FTQUAL_map},
  {FEATDEF_3clip ,  FTQUAL_standard_name},
  {FEATDEF_3clip ,  FTQUAL_usedin},

  {FEATDEF_3UTR ,  FTQUAL_allele},
  {FEATDEF_3UTR ,  FTQUAL_function},
  {FEATDEF_3UTR ,  FTQUAL_label},
  {FEATDEF_3UTR ,  FTQUAL_map},
  {FEATDEF_3UTR ,  FTQUAL_standard_name},
  {FEATDEF_3UTR ,  FTQUAL_usedin},

  {FEATDEF_5clip ,  FTQUAL_allele},
  {FEATDEF_5clip ,  FTQUAL_function},
  {FEATDEF_5clip ,  FTQUAL_label},
  {FEATDEF_5clip ,  FTQUAL_map},
  {FEATDEF_5clip ,  FTQUAL_standard_name},
  {FEATDEF_5clip ,  FTQUAL_usedin},

  {FEATDEF_5UTR ,  FTQUAL_allele},
  {FEATDEF_5UTR ,  FTQUAL_function},
  {FEATDEF_5UTR ,  FTQUAL_label},
  {FEATDEF_5UTR ,  FTQUAL_map},
  {FEATDEF_5UTR ,  FTQUAL_standard_name},
  {FEATDEF_5UTR ,  FTQUAL_usedin},

  {FEATDEF_10_signal ,  FTQUAL_label},
  {FEATDEF_10_signal ,  FTQUAL_map},
  {FEATDEF_10_signal ,  FTQUAL_standard_name},
  {FEATDEF_10_signal ,  FTQUAL_usedin},

  {FEATDEF_35_signal ,  FTQUAL_label},
  {FEATDEF_35_signal ,  FTQUAL_map},
  {FEATDEF_35_signal ,  FTQUAL_standard_name},
  {FEATDEF_35_signal ,  FTQUAL_usedin},

  {FEATDEF_REGION ,  FTQUAL_function},
  {FEATDEF_REGION ,  FTQUAL_label},
  {FEATDEF_REGION ,  FTQUAL_map},
  {FEATDEF_REGION ,  FTQUAL_number},
  {FEATDEF_REGION ,  FTQUAL_phenotype},
  {FEATDEF_REGION ,  FTQUAL_product},
  {FEATDEF_REGION ,  FTQUAL_standard_name},
  {FEATDEF_REGION ,  FTQUAL_usedin},

  {FEATDEF_mat_peptide_aa ,  FTQUAL_label},
  {FEATDEF_mat_peptide_aa ,  FTQUAL_map},
  {FEATDEF_mat_peptide_aa ,  FTQUAL_product},
  {FEATDEF_mat_peptide_aa ,  FTQUAL_standard_name},
  {FEATDEF_mat_peptide_aa ,  FTQUAL_usedin},

  {FEATDEF_sig_peptide_aa ,  FTQUAL_label},
  {FEATDEF_sig_peptide_aa ,  FTQUAL_map},
  {FEATDEF_sig_peptide_aa ,  FTQUAL_product},
  {FEATDEF_sig_peptide_aa ,  FTQUAL_standard_name},
  {FEATDEF_sig_peptide_aa ,  FTQUAL_usedin},

  {FEATDEF_transit_peptide_aa ,  FTQUAL_label},
  {FEATDEF_transit_peptide_aa ,  FTQUAL_map},
  {FEATDEF_transit_peptide_aa ,  FTQUAL_product},
  {FEATDEF_transit_peptide_aa ,  FTQUAL_standard_name},
  {FEATDEF_transit_peptide_aa ,  FTQUAL_usedin},

  {FEATDEF_snoRNA ,  FTQUAL_function},
  {FEATDEF_snoRNA ,  FTQUAL_label},
  {FEATDEF_snoRNA ,  FTQUAL_map},
  {FEATDEF_snoRNA ,  FTQUAL_product},
  {FEATDEF_snoRNA ,  FTQUAL_standard_name},
  {FEATDEF_snoRNA ,  FTQUAL_usedin}

};

/* comparison of ValQual's -- first compare featdef then ftqual */

/* macro did not work properly on linux machine, so using function instead */
/* #define COMPARE_VALQUAL(av,aq,bv,bq) ( ((av)-(bv)) ? ((av)-(bv)) : ((aq)-(bq)) ) */

static Int2 CompareValQual (Uint2 av, FtQualType aq, Uint2 bv, FtQualType bq)

{
  if (av == bv) return (aq - bq);
  return (av - bv);
}

/* Returns TRUE if {featureKey, qualKey} exists in legalGbqualList */

static Boolean AllowedValQual (Uint2 featureKey, FtQualType qualKey)

{
  Int2 L, R, mid;

  L = 0;
  R = sizeof (legalGbqualList) / sizeof (ValQual) - 1;
  while (L < R) {
    mid = (L + R) / 2;
    if (CompareValQual (legalGbqualList [mid].featdef,
       legalGbqualList [mid].ftqual,
       featureKey, qualKey) < 0)
      L = mid + 1;
    else
      R = mid;
  }
  if (CompareValQual (legalGbqualList [R].featdef,
      legalGbqualList [R].ftqual,
      featureKey, qualKey)) {
    return 0;
  }
  return 1;
}


static CharPtr validRptString [] = {
  "tandem", "inverted", "flanking", "terminal", "direct", "dispersed", "other", NULL
};

static CharPtr validLRBString [] = {
  "LEFT", "RIGHT", "BOTH", NULL
};

static CharPtr validConsSpliceString [] = {
  "(5'site:YES, 3'site:YES)",
  "(5'site:YES, 3'site:NO)",
  "(5'site:YES, 3'site:ABSENT)",
  "(5'site:NO, 3'site:YES)",
  "(5'site:NO, 3'site:NO)",
  "(5'site:NO, 3'site:ABSENT)",
  "(5'site:ABSENT, 3'site:YES)",
  "(5'site:ABSENT, 3'site:NO)",
  "(5'site:ABSENT, 3'site:ABSENT)",
  NULL
};

static CharPtr validExceptionString [] = {
  "RNA editing", "reasons given in citation", NULL
};

static Boolean StringInStringList (CharPtr testString, CharPtr PNTR stringList) {
  Int2 i;
  i = 0;
  while (stringList [i] != NULL) {
    if (StringICmp (testString, stringList [i]) == 0)
      return 1;
    i++;
  }
  return 0;
}

/*
Functions now public and prototyped in sequtil.h
Return values are:
 0: no problem - Accession is in proper format
-1: Accession did not start with a letter (or two letters)
-2: Accession did not contain five numbers (or six numbers after 2 letters)
-3: the original Accession number to be validated was NULL
-4: the original Accession number is too long (>16)
*/

NLM_EXTERN Int2 ValidateAccn (
  CharPtr accession
)

{
  Char     ch;
  Int2     numAlpha = 0;
  Int2     numDigits = 0;
  Int2     numUndersc = 0;
  CharPtr  str;

  if (accession == NULL || accession [0] == '\0') return -3;

  if (StringLen (accession) >= 16) return -4;

  if (accession [0] < 'A' || accession [0] > 'Z') return -1;

  str = accession;
  if (StringNCmp (str, "NZ_", 3) == 0) {
    str += 3;
  }
  ch = *str;
  while (IS_ALPHA (ch)) {
    numAlpha++;
    str++;
    ch = *str;
  }
  while (ch == '_') {
    numUndersc++;
    str++;
    ch = *str;
  }
  while (IS_DIGIT (ch)) {
    numDigits++;
    str++;
    ch = *str;
  }
  if (ch != '\0' && ch != ' ' && ch != '.') return -2;

  if (numUndersc > 1) return -2;

  if (numUndersc == 0) {
    if (numAlpha == 1 && numDigits == 5) return 0;
    if (numAlpha == 2 && numDigits == 6) return 0;
    if (numAlpha == 3 && numDigits == 5) return 0;
    if (numAlpha == 4 && numDigits == 8) return 0;
  } else if (numUndersc == 1) {
    if (numAlpha != 2 || numDigits != 6) return -2;
    if (accession [0] == 'N' || accession [0] == 'X' || accession [0] == 'Z') {
      if (accession [1] == 'M' ||
          accession [1] == 'C' ||
          accession [1] == 'T' ||
          accession [1] == 'P' ||
          accession [1] == 'G' ||
          accession [1] == 'R' ||
          accession [1] == 'S' ||
          accession [1] == 'Z') {
        return 0;
      }
    }
  }

  return -2;
}

NLM_EXTERN Int2 ValidateSeqID (
  SeqIdPtr sip
)

{
  Char  buf [41];

  if (sip == NULL) return -3;
  SeqIdWrite (sip, buf, PRINTID_TEXTID_ACC_VER, sizeof (buf) - 1);
  return ValidateAccn (buf);
}

static CharPtr mrnaevtext1 = "Derived by automated computational analysis";
static CharPtr mrnaevtext2 = "using gene prediction method:";
static CharPtr mrnaevtext3 = "Supporting evidence includes similarity to:";

static void GetStrFormRNAEvidence (
  UserObjectPtr uop,
  Pointer userdata
)

{
  size_t        len;
  CharPtr       method = NULL;
  Int2          nm = 0;
  ObjectIdPtr   oip;
  CharPtr       str = NULL;
  CharPtr PNTR  strp;
  Char          tmp [20];
  UserFieldPtr  u, ufp, uu;

  if (uop == NULL) return;
  oip = uop->type;
  if (oip == NULL) return;
  if (StringCmp (oip->str, "ModelEvidence") != 0) return;
  strp = (CharPtr PNTR) userdata;

  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    oip = ufp->label;
    if (oip == NULL || ufp->data.ptrvalue == NULL) continue;
    if (StringCmp (oip->str, "Method") == 0) {
      method = StringSaveNoNull ((CharPtr) ufp->data.ptrvalue);
    }
    if (StringCmp (oip->str, "mRNA") == 0) {
      for (u = (UserFieldPtr) ufp->data.ptrvalue; u != NULL; u = u->next) {
        if (u->data.ptrvalue == NULL) continue;
        for (uu = (UserFieldPtr) u->data.ptrvalue; uu != NULL; uu = uu->next) {
          oip = uu->label;
          if (oip == NULL) continue;
          if (StringCmp (oip->str, "accession") == 0) {
            nm++;
          }
        }
      }
    }
  }

  len = StringLen (mrnaevtext1) + StringLen (mrnaevtext2) + StringLen (mrnaevtext3) + StringLen (method) + 30;
  str = (CharPtr) MemNew (len);
  if (str == NULL) return;

  if (method != NULL) {
    sprintf (str, "%s %s %s.", mrnaevtext1, mrnaevtext2, method);
  } else {
    sprintf (str, "%s.", mrnaevtext1);
  }
  if (nm > 0) {
    StringCat (str, " ");
    StringCat (str, mrnaevtext3);
    if (nm > 1) {
      sprintf (tmp, " %d mRNAs", (int) nm);
    } else {
      sprintf (tmp, " %d mRNA", (int) nm);
    }
    StringCat (str, tmp);
  }

  *strp = str;
}

static Boolean ValidateRptUnit (
  CharPtr buf
)

{
  CharPtr  str;
  Char     tmp [255];

  StringNCpy_0 (tmp, buf, sizeof (tmp));
  TrimSpacesAroundString (tmp);

  str = tmp;
  /* first check for sequence letters with optional semicolons */
  while (IS_ALPHA (*str) || *str == ';') str++;
  if (*str == '\0') return TRUE;
  /* next check for letters, digits, commas, parentheses, dashes, and underscores */
  str = tmp;
  while (IS_ALPHANUM (*str) || *str == '(' || *str == ')' || *str == ',' || *str == ';' || *str == '-' || *str == '_') str++;
  if (*str == '\0') return TRUE;
  /* now check for officially legal styles */
  str = tmp;
  while (IS_ALPHANUM (*str)) str++;
  if (*str != '\0') { /* wasn't pure alphanumeric; now check for xxx..yyy */
    str = buf;
    while (IS_DIGIT (*str)) str++; /* xxx */
    if (*str == '\0' /* must be something after the xxx */
      || StringLen (str) < 3  /* need at least 2 '.'s and a digit*/
      || str[0] != '.' || str[1] != '.') return FALSE;
    str+=2;
    while (IS_DIGIT (*str)) str++;
    if (*str != '\0') return FALSE;  /* mustn't be anything after the yyy */
  }
  return TRUE;
}



static void RecordUserObjectsInQVP (
  UserObjectPtr uop,
  Pointer userdata
)

{
  ObjectIdPtr  oip;
  QualValPtr   qvp;

  if (uop == NULL || userdata == NULL) return;
  qvp = (QualValPtr) userdata;
  oip = uop->type;
  if (oip == NULL) return;
  if (StringCmp (oip->str, "ModelEvidence") == 0) {
    qvp [FTQUAL_modelev].uop = uop;
  }
}

static SeqIdPtr SeqLocIdForProduct (
  SeqLocPtr product
)

{
  SeqIdPtr   sip;
  SeqLocPtr  slp;

  /* in case product is a SEQLOC_EQUIV */

  if (product == NULL) return NULL;
  sip = SeqLocId (product);
  if (sip != NULL) return sip;
  slp = SeqLocFindNext (product, NULL);
  while (slp != NULL) {
    sip = SeqLocId (slp);
    if (sip != NULL) return sip;
    slp = SeqLocFindNext (product, slp);
  }
  return NULL;
}

static void AddIntervalsToGbfeat (
  GBFeaturePtr gbfeat,
  SeqLocPtr location,
  BioseqPtr target
)

{
  Char           accn [41];
  SeqLocPtr      copy = NULL;
  Int4           from;
  GBIntervalPtr  gbint;
  Int4           gi;
  GBIntervalPtr  last = NULL;
  Int4           point;
  SeqIntPtr      sint;
  SeqIdPtr       sip;
  SeqLocPtr      slp;
  SeqPntPtr      spp;
  Int4           to;
  Int4           swap;

  if (gbfeat == NULL || location == NULL) return;
  if (target != NULL) {
    copy = SeqLocMerge (target, location, NULL, FALSE, TRUE, FALSE);
    location = copy;
  }

  slp = SeqLocFindNext (location, NULL);
  while (slp != NULL) {
    from = 0;
    to = 0;
    point = 0;
    sip = NULL;
    switch (slp->choice) {
      case SEQLOC_INT :
        sint = (SeqIntPtr) slp->data.ptrvalue;
        if (sint != NULL) {
          from = sint->from + 1;
          to = sint->to + 1;
          sip = sint->id;
          if (sint->strand == Seq_strand_minus && from < to) {
            swap = from;
            from = to;
            to = swap;
          }
        }
        break;
      case SEQLOC_PNT :
        spp = (SeqPntPtr) slp->data.ptrvalue;
        if (spp != NULL) {
          point = spp->point + 1;
          sip = spp->id;
        }
        break;
      default :
        break;
    }
    if (sip != NULL) {
      accn [0] = '\0';
      if (sip->choice == SEQID_GI) {
        gi = sip->data.intvalue;
        if (! GetAccnVerFromServer (gi, accn)) {
          accn [0] = '\0';
        }
        if (StringHasNoText (accn)) {
          sip = GetSeqIdForGI (gi);
          SeqIdWrite (sip, accn, PRINTID_TEXTID_ACC_VER, sizeof (accn));
          SeqIdFree (sip);
        }
      } else {
        SeqIdWrite (sip, accn, PRINTID_TEXTID_ACC_VER, sizeof (accn));
      }
      if (! StringHasNoText (accn)) {
        gbint = GBIntervalNew ();
        if (gbint != NULL) {
          gbint->from = from;
          gbint->to = to;
          gbint->point = point;
          gbint->accession = StringSave (accn);
          if (gbfeat->intervals == NULL) {
            gbfeat->intervals = gbint;
          } else {
            last->next = gbint;
          }
          last = gbint;
        }
      }
    }
    slp = SeqLocFindNext (location, slp);
  }

  SeqLocFree (copy);
}

static CharPtr FormatBasecountBlock (
  Asn2gbFormatPtr afp,
  BaseBlockPtr bbp
)

{
  IntAsn2gbJobPtr   ajp;
  Asn2gbSectPtr     asp;
  Int4              base_count [5];
  BioseqPtr         bsp;
  Char              buf [80];
  Byte              bases [400];
  Uint1             code = Seq_code_iupacna;
  Int2              ctr;
  Int2              i;
  Int4              len;
  Uint1             residue;
  SeqPortPtr        spp = NULL;
  Int4              total = 0;
  StringItemPtr ffstring;
  CharPtr str;

  if (afp == NULL || bbp == NULL) return NULL;
  ajp = afp->ajp;
  if (ajp == NULL) return NULL;

  asp = afp->asp;
  if (asp == NULL) return NULL;
  bsp = (asp->bsp);
  if (bsp == NULL) return NULL;

  /* after first formatting, result is cached into bbp->string */

  if (! StringHasNoText (bbp->string)) return StringSave (bbp->string);

  for (i = 0; i < 5; i++) {
    base_count [i] = 0;
  }

  if (ISA_aa (bsp->mol)) {
    code = Seq_code_ncbieaa;
  }

  if (ajp->ajp.slp != NULL) {
    spp = SeqPortNewByLoc (ajp->ajp.slp, code);
    len = SeqLocLen (ajp->ajp.slp);
  } else {
    spp = SeqPortNew (bsp, 0, -1, 0, code);
    len = bsp->length;
  }
  if (spp == NULL) return NULL;
  if (bsp->repr == Seq_repr_delta || bsp->repr == Seq_repr_virtual) {
    SeqPortSet_do_virtual (spp, TRUE);
  }

  /* use SeqPortRead rather than SeqPortGetResidue for faster performance */

  ctr = SeqPortRead (spp, bases, sizeof (bases));
  i = 0;
  residue = (Uint1) bases [i];
  while (residue != SEQPORT_EOF) {
    if (IS_residue (residue)) {
      total++;
      switch (residue) {
        case 'A' :
          (base_count [0])++;
          break;
        case 'C' :
          (base_count [1])++;
          break;
        case 'G' :
          (base_count [2])++;
          break;
        case 'T' :
          (base_count [3])++;
          break;
        default :
          (base_count [4])++;
          break;
      }
    }
    i++;
    if (i >= ctr) {
      i = 0;
      ctr = SeqPortRead (spp, bases, sizeof (bases));
      if (ctr < 0) {
        bases [0] = -ctr;
      } else if (ctr < 1) {
        bases [0] = SEQPORT_EOF;
      }
    }
    residue = (Uint1) bases [i];
  }

  SeqPortFree (spp);

  if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {

    if (base_count [4] == 0) {
      sprintf (buf, "%7ld a%7ld c%7ld g%7ld t",
               (long) base_count [0], (long) base_count [1],
               (long) base_count [2], (long) base_count [3]);
    } else {
      sprintf (buf, "%7ld a%7ld c%7ld g%7ld t%7ld others",
               (long) base_count [0], (long) base_count [1],
               (long) base_count [2], (long) base_count [3],
               (long) base_count [4]);
    }

  } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {

    sprintf (buf, "Sequence %ld BP; %ld A; %ld C; %ld G; %ld T; %ld other;",
             (long) len,
             (long) base_count [0], (long) base_count [1],
             (long) base_count [2], (long) base_count [3],
             (long) base_count [4]);
  }

  ffstring = FFGetString(ajp);
  if ( ffstring == NULL ) return NULL;

  if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {
    FFAddOneString(ffstring, "XX\n", FALSE, FALSE, TILDE_IGNORE);
  }
  FFStartPrint (ffstring, afp->format, 0, 0, "BASE COUNT", 12, 5, 5, "SQ", FALSE);
  FFAddOneString (ffstring, buf, FALSE, FALSE, TILDE_TO_SPACES);
  str = FFEndPrint(ajp, ffstring, afp->format, 12, 0, 5, 5, "SQ");
  FFRecycleString(ajp, ffstring);

  return str;
}

static void PrintSeqLine (
  StringItemPtr ffstring,
  FmtType format,
  CharPtr buf,
  Int4 start,
  Int4 stop
)

{
  size_t  len;
  Char    pos [16];
  Int4    pad;

  len = StringLen (buf);
  if (len > 0 && buf [len - 1] == ' ') {
    buf [len - 1] = '\0';
  }

  if (format == GENBANK_FMT || format == GENPEPT_FMT) {

    sprintf (pos, "%9ld", (long) (start + 1));
    FFAddOneString(ffstring, pos, FALSE, FALSE, TILDE_TO_SPACES);
    FFAddOneChar(ffstring, ' ', FALSE);
    FFAddOneString(ffstring, buf, FALSE, FALSE, TILDE_TO_SPACES);
    FFAddOneChar(ffstring, '\n', FALSE);
  } else if (format == EMBL_FMT || format == EMBLPEPT_FMT) {

    sprintf (pos, "%8ld", (long) (stop));
    FFAddNChar(ffstring, ' ', 5, FALSE);
    FFAddOneString(ffstring, buf, FALSE, FALSE, TILDE_TO_SPACES);
    pad = 72 - 5 - StringLen(buf);
    FFAddNChar(ffstring, ' ', pad, FALSE);
    FFAddOneString(ffstring, pos, FALSE, FALSE, TILDE_TO_SPACES);
    FFAddOneChar(ffstring, '\n', FALSE);
  }
}

static CharPtr CompressNonBases (CharPtr str)

{
  Char     ch;
  CharPtr  dst;
  CharPtr  ptr;

  if (str == NULL || str [0] == '\0') return NULL;

  dst = str;
  ptr = str;
  ch = *ptr;
  while (ch != '\0') {
    if (IS_ALPHA (ch)) {
      *dst = ch;
      dst++;
    }
    ptr++;
    ch = *ptr;
  }
  *dst = '\0';

  return str;
}

static void CatenateSequenceInGbseq (
  GBSeqPtr gbseq,
  CharPtr str
)

{
  Char     ch;
  CharPtr  tmp;

  if (gbseq == NULL || StringHasNoText (str)) return;

  if (gbseq->sequence == NULL) {
    gbseq->sequence = StringSave (str);
  } else {
    tmp = (CharPtr) MemNew (StringLen (gbseq->sequence) + StringLen (str) + 2);
    StringCpy (tmp, gbseq->sequence);
    StringCat (tmp, str);
    gbseq->sequence = MemFree (gbseq->sequence);
    gbseq->sequence = tmp;
  }

  tmp = gbseq->sequence;
  if (tmp == NULL) return;
  ch = *tmp;
  while (ch != '\0') {
    if (ch == '\n' || ch == '\r' || ch == '\t') {
      *tmp = ' ';
    }
    tmp++;
    ch = *tmp;
  }
  TrimSpacesAroundString (gbseq->sequence);
  CompressNonBases (gbseq->sequence);
}

  static Uint1 fasta_order [NUM_SEQID] = {
    33, /* 0 = not set */
    20, /* 1 = local Object-id */
    15, /* 2 = gibbsq */
    16, /* 3 = gibbmt */
    30, /* 4 = giim Giimport-id */
    10, /* 5 = genbank */
    10, /* 6 = embl */
    10, /* 7 = pir */
    10, /* 8 = swissprot */
    15, /* 9 = patent */
    20, /* 10 = other TextSeqId */
    20, /* 11 = general Dbtag */
    255, /* 12 = gi */
    10, /* 13 = ddbj */
    10, /* 14 = prf */
    12, /* 15 = pdb */
    10, /* 16 = tpg */
    10, /* 17 = tpe */
    10  /* 18 = tpd */
  };

static void PrintGenome (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  SeqLocPtr slp_head, 
  CharPtr prefix, 
  Boolean segWithParts
)
{
  Char         buf[40], val[166];
  Boolean      first = TRUE;
  SeqLocPtr    slp;
  Int4         from, to, start, stop;
  SeqIdPtr     sid, newid;
  BioseqPtr    bsp = NULL;
  Int2         p1=0, p2=0;

  for (slp = slp_head; slp; slp = slp->next) {
    from = to = 0;
    sid = SeqLocId(slp);
    if (slp->choice == SEQLOC_INT || slp->choice == SEQLOC_WHOLE) {
      start = from = SeqLocStart(slp);
      stop = to = SeqLocStop(slp);
    } else if (slp->choice == SEQLOC_NULL){
      sprintf(val, ",%s", "gap()");
      FFAddOneString(ffstring, val, FALSE, FALSE, TILDE_IGNORE);
      continue;
    } else {
      continue;
    }
    if (sid == NULL) {
      continue;
    }
    newid = NULL;
    buf [0] = '\0';
    if (sid->choice == SEQID_GI) {
      if (GetAccnVerFromServer (sid->data.intvalue, buf)) {
        /* no need to call GetSeqIdForGI */
      } else {
        newid = GetSeqIdForGI (sid->data.intvalue);
        if (newid != NULL && segWithParts) {
          if (newid->choice == SEQID_GIBBSQ ||
              newid->choice == SEQID_GIBBMT ||
              newid->choice == SEQID_GIIM) {
            bsp = BioseqFind (newid);
            if (bsp != NULL && bsp->repr == Seq_repr_virtual) {
              if (bsp->length > 0) {
                sprintf (val, ",gap(%ld)", (long) bsp->length);
              } else {
                sprintf(val, ",%s", "gap()");
              }
              FFAddOneString(ffstring, val, FALSE, FALSE, TILDE_IGNORE);
              continue;
            }
          }
        }
      }
    } else if (sid->choice == SEQID_GENERAL) {
      newid = sid;
    } else {
      newid = sid;
    }
    if (prefix != NULL) {
      FFAddOneString(ffstring, prefix, FALSE, FALSE, TILDE_IGNORE);
    }
    if (first) {
      first = FALSE;
    } else {
      FFAddOneChar(ffstring, ',', FALSE);
      /*ff_AddChar(',');*/
    }
    if (! StringHasNoText (buf)) {
      /* filled in by GetAccnVerFromServer */
    } else if (newid != NULL) {
      SeqIdWrite (SeqIdSelect (newid, fasta_order, NUM_SEQID),
                 buf, PRINTID_TEXTID_ACC_VER, sizeof(buf) -1 );
    } else if (sid->choice == SEQID_GI) {
      SeqIdWrite (sid, buf, PRINTID_FASTA_LONG, sizeof (buf) - 1);
    }

    if (SeqLocStrand (slp) == Seq_strand_minus) {
      FFAddOneString(ffstring, "complement(", FALSE, FALSE, TILDE_IGNORE);
    }
    if ( GetWWW(ajp) ) {
      if (newid == NULL) {
        newid = sid;
      }
      if (newid->choice != SEQID_GENERAL) {
        FFAddTextToString(ffstring, "<a href=", link_seq, NULL, FALSE, FALSE, TILDE_IGNORE);
        FFAddTextToString(ffstring, "val=", buf, ">", FALSE, FALSE, TILDE_IGNORE);
        FFAddTextToString(ffstring, NULL, buf, "</a>", FALSE, FALSE, TILDE_IGNORE);
      }
    } else {
      FFAddOneString(ffstring, buf, FALSE, FALSE, TILDE_IGNORE);
    }

    if (SeqLocStrand(slp) == Seq_strand_minus) {
      sprintf (val,":%ld..%ld)", (long) start+1, (long) stop+1);
    } else {
      sprintf (val,":%ld..%ld", (long) start+1, (long) stop+1);
    }
    FFAddOneString(ffstring, val, FALSE, FALSE, TILDE_IGNORE);
    p1 += StringLen (val);
    p2 += StringLen (val);
  }
}

static Boolean SegHasParts (
  BioseqPtr bsp
)

{
  BioseqSetPtr  bssp;
  SeqEntryPtr   sep;

  if (bsp == NULL || bsp->repr != Seq_repr_seg) return FALSE;
  sep = bsp->seqentry;
  if (sep == NULL) return FALSE;
  sep = sep->next;
  if (sep == NULL || (! IS_Bioseq_set (sep))) return FALSE;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp != NULL && bssp->_class == BioseqseqSet_class_parts) return TRUE;
  return FALSE;
}

static Boolean DeltaLitOnly (
  BioseqPtr bsp
)

{
  ValNodePtr  vnp;

  if (bsp == NULL || bsp->repr != Seq_repr_delta) return FALSE;
  for (vnp = (ValNodePtr)(bsp->seq_ext); vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == 1) return FALSE;
  }
  return TRUE;
}

static CharPtr FormatContigBlock (
  Asn2gbFormatPtr afp,
  BaseBlockPtr bbp
)

{
  IntAsn2gbJobPtr  ajp;
  Asn2gbSectPtr    asp;
  BioseqPtr        bsp;
  DeltaSeqPtr      dsp;
  GBSeqPtr         gbseq;
  SeqLitPtr        litp;
  CharPtr          prefix = NULL;
  Boolean          segWithParts = FALSE;
  SeqLocPtr        slp_head = NULL;
  CharPtr          str;
  Char             val [20];
  StringItemPtr    ffstring;
  CharPtr          label;

  if (afp == NULL || bbp == NULL) return NULL;
  ajp = afp->ajp;
  if (ajp == NULL) return NULL;
  asp = afp->asp;
  if (asp == NULL) return NULL;
  bsp = (asp->bsp);
  if (bsp == NULL) return NULL;

  ffstring = FFGetString(ajp);
  if ( ffstring == NULL ) return NULL;

  if ( GetWWW(ajp) ) {
    label = "CONTIG   ";
  } else {
    label = "CONTIG";
  }
  FFAddOneString(ffstring, label,  FALSE, FALSE, TILDE_IGNORE);  
  FFAddNChar(ffstring, ' ', 12 - StringLen(label), FALSE);

  FFAddOneString(ffstring, "join(", FALSE, FALSE, TILDE_IGNORE);

  if (bsp->seq_ext_type == 1) {

    if (bsp->repr == Seq_repr_seg && SegHasParts (bsp)) {
      segWithParts = TRUE;
    }

    slp_head = (SeqLocPtr) bsp->seq_ext;
    PrintGenome (ajp, ffstring, slp_head, prefix, segWithParts);

  } else if (bsp->seq_ext_type == 4) {

    for (dsp = (DeltaSeqPtr) bsp->seq_ext; dsp != NULL; dsp=dsp->next) {
      if (dsp->choice == 1) {

        slp_head = (SeqLocPtr) dsp->data.ptrvalue;
        PrintGenome (ajp, ffstring, slp_head, prefix, FALSE);

      } else {

        litp = (SeqLitPtr) dsp->data.ptrvalue;
        if (litp != NULL) {
          if (litp->seq_data != NULL) {
            if (litp->length == 0) {
              sprintf (val, "gap(%ld)", (long) litp->length);
              FFAddOneString(ffstring, val, FALSE, FALSE, TILDE_IGNORE);
            } else {
              /* don't know what to do here */
            }
          } else {
            sprintf (val, ",gap(%ld)", (long) litp->length);
            FFAddOneString(ffstring, val, FALSE, FALSE, TILDE_IGNORE);
          }
        }
      }

      prefix = ",";
    }
  }

  FFAddOneChar(ffstring, ')', FALSE);

  str = FFEndPrint(ajp, ffstring, afp->format, 12, 12, 12, 12, NULL);
  FFRecycleString(ajp, ffstring);

  /* optionally populate gbseq for XML-ized GenBank format */

  if (ajp->gbseq) {
    gbseq = &asp->gbseq;
  } else {
    gbseq = NULL;
  }

  if (gbseq != NULL) {
    if (StringLen (str) > 12) {
      gbseq->contig = StringSave (str + 12);
    } else {
      gbseq->contig = StringSave (str);
    }
  }

  return str;
}

static CharPtr FormatSlashBlock (
  Asn2gbFormatPtr afp,
  BaseBlockPtr bbp
)

{
  IntAsn2gbJobPtr    ajp;
  Asn2gbSectPtr      asp;
  GBFeaturePtr       currf, headf, nextf;
  GBReferencePtr     currr, headr, nextr;
  GBSeqPtr           gbseq;
  IndxPtr            index;

  if (afp == NULL || bbp == NULL) return NULL;
  ajp = afp->ajp;
  if (ajp == NULL) return NULL;
  asp = afp->asp;
  if (asp == NULL) return NULL;

  /* sort and unique indexes */

  index = ajp->index;

  if (index != NULL) {

    MemCopy (index, &asp->index, sizeof (IndxBlock));
    MemSet (&asp->index, 0, sizeof (IndxBlock));

    index->authors = ValNodeSort (index->authors, SortVnpByString);
    index->authors = UniqueValNode (index->authors);

    index->genes = ValNodeSort (index->genes, SortVnpByString);
    index->genes = UniqueValNode (index->genes);

    index->journals = ValNodeSort (index->journals, SortVnpByString);
    index->journals = UniqueValNode (index->journals);

    index->keywords = ValNodeSort (index->keywords, SortVnpByString);
    index->keywords = UniqueValNode (index->keywords);

    index->secondaries = ValNodeSort (index->secondaries, SortVnpByString);
    index->secondaries = UniqueValNode (index->secondaries);
  }

  /* adjust XML-ized GenBank format */

  gbseq = ajp->gbseq;

  if (gbseq != NULL) {

    MemCopy (gbseq, &asp->gbseq, sizeof (GBSeq));
    MemSet (&asp->gbseq, 0, sizeof (GBSeq));

    /* reverse order of references */

    headr = NULL;
    for (currr = gbseq->references; currr != NULL; currr = nextr) {
      nextr = currr->next;
      currr->next = headr;
      headr = currr;
    }
    gbseq->references = headr;

    /* reverse order of features */

    headf = NULL;
    for (currf = gbseq->feature_table; currf != NULL; currf = nextf) {
      nextf = currf->next;
      currf->next = headf;
      headf = currf;
    }
    gbseq->feature_table = headf;
  }

  /* slash always has string pre-allocated by add slash block function */

  return StringSaveNoNull (bbp->string);
}

/* ********************************************************************** */

/* functions to record sections or blocks in linked lists */

static Asn2gbSectPtr Asn2gbAddSection (
  Asn2gbWorkPtr awp
)

{
  Asn2gbSectPtr  asp;
  ValNodePtr     vnp;

  if (awp == NULL) return NULL;

  asp = (Asn2gbSectPtr) MemNew (sizeof (IntAsn2gbSect));
  if (asp == NULL) return NULL;

  vnp = ValNodeAddPointer (&(awp->lastsection), 0, asp);
  if (vnp == NULL) return asp;

  awp->lastsection = vnp;
  if (awp->sectionList == NULL) {
    awp->sectionList = vnp;
  }

  return asp;
}

static BaseBlockPtr Asn2gbAddBlock (
  Asn2gbWorkPtr awp,
  BlockType blocktype,
  size_t size
)

{
  BaseBlockPtr  bbp;
  ValNodePtr    vnp;

  if (awp == NULL || size < 1) return NULL;

  bbp = (BaseBlockPtr) MemNew (size);
  if (bbp == NULL) return NULL;
  bbp->blocktype = blocktype;
  bbp->section = awp->currsection;

  vnp = ValNodeAddPointer (&(awp->lastblock), 0, bbp);
  if (vnp == NULL) return bbp;

  awp->lastblock = vnp;
  if (awp->blockList == NULL) {
    awp->blockList = vnp;
  }

  return bbp;
}


/* ********************************************************************** */

/* add functions allocate specific blocks, populate with paragraph print info */

static CharPtr strd [4] = {
  "   ", "ss-", "ds-", "ms-"
};

static CharPtr gnbk_mol [14] = {
  "    ", "DNA ", "RNA ", "mRNA", "rRNA", "tRNA", "uRNA", "scRNA", " AA ", "DNA ", "DNA ", "RNA ", "snoRNA", "RNA "
};

/* EMBL_FMT in RELEASE_MODE or ENTREZ_MODE, otherwise use gnbk_mol */

static CharPtr embl_mol [14] = {
  "xxx", "DNA", "RNA", "RNA", "RNA", "RNA", "RNA", "RNA", "AA ", "DNA", "DNA", "RNA", "RNA", "RNA"
};

static CharPtr embl_divs [18] = {
  "FUN", "INV", "MAM", "ORG", "PHG", "PLN", "PRI", "PRO", "ROD"
  "SYN", "UNA", "VRL", "VRT", "PAT", "EST", "STS", "HUM", "HTC"
};

static Uint1 imolToMoltype [14] = {
  0, 1, 2, 5, 4, 3, 6, 7, 9, 1, 1, 2, 8, 2
};

static DatePtr GetBestDate (
  DatePtr a,
  DatePtr b
)

{
  Int2  status;

  if (a == NULL) return b;
  if (b == NULL) return a;

  status = DateMatch (a, b, FALSE);
  if (status == 1) return a;

  return b;
}

/*--------------------------------------------------------*/
/*                                                        */
/*  s_IsSeperatorNeeded()                                 */
/*                                                        */
/*--------------------------------------------------------*/

static Boolean s_IsSeperatorNeeded(CharPtr baseString, Int4 baseLength, Int2 suffixLength)
{
  Char lastChar;
  Char nextToLastChar;

  lastChar = baseString[baseLength - 1];
  nextToLastChar = baseString[baseLength - 2];

  /* This first check put here to emulate what may be a  */
  /* bug in the original code (in CheckLocusLength() )   */
  /* which adds an 'S' segment seperator only if it      */
  /* DOES make the string longer than the max.           */

  if (baseLength + suffixLength < 16)
    return FALSE;

  /* If the last character is not a digit */
  /* then don't use a seperator.          */

  if (!IS_DIGIT(lastChar))
    return FALSE;

  /* If the last two characters are a non-digit   */
  /* followed by a '0', then don't use seperator. */

  if ((lastChar == '0') && (!IS_DIGIT(nextToLastChar)))
    return FALSE;

  /* If we made it to here, use a seperator */

  return TRUE;
}

/*--------------------------------------------------------*/
/*                                                        */
/*  s_LocusAddSuffix() -                                  */
/*                                                        */
/*--------------------------------------------------------*/

static Boolean s_LocusAddSuffix (CharPtr locus, Asn2gbWorkPtr awp)
{
  size_t  buflen;
  Char    ch;
  Char    segCountStr[6];
  Int2    segCountStrLen;
  Char    segSuffix[5];

  buflen = StringLen (locus);

  /* If there's one or less segments, */
  /* no suffix is needed.             */

  if (awp->numsegs <= 1)
    return FALSE;

  /* If the basestring has one or less */
  /* characters, no suffix is needed.  */

  if (buflen <=1)
    return FALSE;

  /* Add the suffix */

  ch = locus[buflen-1];
  sprintf(segCountStr,"%d",awp->numsegs);
  segCountStrLen = StringLen(segCountStr);
  segSuffix[0] = '\0';

  if (s_IsSeperatorNeeded(locus,buflen,segCountStrLen) == TRUE)
    sprintf(segSuffix,"S%0*d",segCountStrLen,awp->seg);
  else
    sprintf(segSuffix,"%0*d",segCountStrLen,awp->seg);
  StringCat(locus,segSuffix);

  /* Return successfully */

  return TRUE;
}

/*--------------------------------------------------------*/
/*                                                        */
/*  s_LocusAdjustLength() -                               */
/*                                                        */
/*--------------------------------------------------------*/

static Boolean s_LocusAdjustLength(CharPtr locus, Int2 maxLength)
{
  Int2     trimCount;
  Int2     buflen;
  CharPtr  buftmp;

  buflen = StringLen (locus);
  if (buflen <= maxLength) return FALSE;

  buftmp = MemNew(maxLength + 1);

  /* If the sequence id is an NCBI locus of the */
  /* form HSU00001, then make sure that if      */
  /* there is trimming the HS gets trimmed off  */
  /* as a unit, never just the 'H'.             */

  trimCount = buflen - maxLength;
  if (trimCount == 1)
    if (IS_ALPHA(locus[0]) != 0 &&
        IS_ALPHA(locus[1]) != 0 &&
        IS_ALPHA(locus[2]) != 0 &&
        IS_DIGIT(locus[3]) != 0 &&
        IS_DIGIT(locus[4]) != 0 &&
        IS_DIGIT(locus[5]) != 0 &&
        IS_DIGIT(locus[6]) != 0 &&
        IS_DIGIT(locus[7]) != 0 &&
        locus[8] == 'S' &&
        locus[9] == '\0')
      trimCount++;

  /* Left truncate the sequence id */

  StringCpy(buftmp, &locus[trimCount]);
  StringCpy(locus, buftmp);

  MemFree(buftmp);
  return TRUE;
}

/*--------------------------------------------------------*/
/*                                                        */
/*  AddLocusBlock() -                                     */
/*                                                        */
/*--------------------------------------------------------*/

static DatePtr GetBestDateForBsp (
  BioseqPtr bsp
)

{
  DatePtr            best_date = NULL;
  SeqMgrDescContext  dcontext;
  DatePtr            dp;
  EMBLBlockPtr       ebp;
  GBBlockPtr         gbp;
  PdbBlockPtr        pdp;
  PdbRepPtr          prp;
  SeqDescrPtr        sdp;
  SPBlockPtr         spp;

  if (bsp == NULL) return NULL;

  dp = NULL;
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_update_date, &dcontext);
  if (sdp != NULL) {
    dp = (DatePtr) sdp->data.ptrvalue;
    best_date = GetBestDate (dp, best_date);
  }

  /* !!! temporarily also look at genbank block entry date !!! */

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_genbank, &dcontext);
  if (sdp != NULL) {
    gbp = (GBBlockPtr) sdp->data.ptrvalue;
    if (gbp != NULL) {
      dp = gbp->entry_date;
      best_date = GetBestDate (dp, best_date);
    }
  }

  /* more complicated code for dates from various objects goes here */

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_embl, &dcontext);
  if (sdp != NULL) {
    ebp = (EMBLBlockPtr) sdp->data.ptrvalue;
    if (ebp != NULL) {
      dp = ebp->creation_date;
      best_date = GetBestDate (dp, best_date);
      dp = ebp->update_date;
      best_date = GetBestDate (dp, best_date);
    }
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_sp, &dcontext);
  if (sdp != NULL) {
    spp = (SPBlockPtr) sdp->data.ptrvalue;
    if (spp != NULL) {
      dp = spp->created;
      if (dp != NULL && dp->data [0] == 1) {
        best_date = GetBestDate (dp, best_date);
      }
      dp = spp->sequpd;
      if (dp != NULL && dp->data [0] == 1) {
        best_date = GetBestDate (dp, best_date);
      }
      dp = spp->annotupd;
      if (dp != NULL && dp->data [0] == 1) {
        best_date = GetBestDate (dp, best_date);
      }
    }
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_pdb, &dcontext);
  if (sdp != NULL) {
    pdp = (PdbBlockPtr) sdp->data.ptrvalue;
    if (pdp != NULL) {
      dp = pdp->deposition;
      if (dp != NULL && dp->data [0] == 1) {
        best_date = GetBestDate (dp, best_date);
      }
      prp = pdp->replace;
      if (prp != NULL) {
        dp = prp->date;
        if (dp != NULL && dp->data[0] == 1) {
          best_date = GetBestDate (dp, best_date);
        }
      }
    }
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_create_date, &dcontext);
  if (sdp != NULL) {
    dp = (DatePtr) sdp->data.ptrvalue;
    if (dp != NULL) {
      best_date = GetBestDate (dp, best_date);
    }
  }

  return best_date;
}

static void AddLocusBlock (
  Asn2gbWorkPtr awp
)

{
  IntAsn2gbJobPtr    ajp;
  Asn2gbSectPtr      asp;
  BaseBlockPtr       bbp;
  DatePtr            best_date = NULL;
  BioSourcePtr       biop;
  Int2               bmol = 0;
  BioseqPtr          bsp;
  SeqFeatPtr         cds;
  Char               date [40];
  SeqMgrDescContext  dcontext;
  Char               div [10];
  BioseqPtr          dna;
  DatePtr            dp;
  CharPtr            ebmol;
  EMBLBlockPtr       ebp;
  SeqMgrFeatContext  fcontext;
  GBBlockPtr         gbp;
  Char               gene [32];
  Boolean            genome_view;
  GBSeqPtr           gbseq;
  Char               id [41];
  Int2               imol = 0;
  IndxPtr            index;
  Int2               istrand;
  Boolean            is_nm = FALSE;
  Boolean            is_np = FALSE;
  Boolean            is_transgenic = FALSE;
  Char               len [32];
  Int4               length;
  Char               locus [41];
  MolInfoPtr         mip;
  Char               mol [30];
  BioseqPtr          nm = NULL;
  OrgNamePtr         onp;
  Uint1              origin;
  OrgRefPtr          orp;
  BioseqPtr          parent;
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp;
  SeqIdPtr           sip;
  SubSourcePtr       ssp;
  Uint1              tech;
  Uint1              topology;
  TextSeqIdPtr       tsip;
  Boolean            wgsmaster = FALSE;
  StringItemPtr      ffstring;

  if (awp == NULL) return;
  ajp = awp->ajp;
  if (ajp == NULL) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;
  asp = awp->asp;
  if (asp == NULL) return;

  bbp = Asn2gbAddBlock (awp, LOCUS_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  ffstring = FFGetString(ajp);
  if ( ffstring == NULL ) return;

  mol [0] = '\0';
  len [0] = '\0';
  div [0] = '\0';
  date [0] = '\0';
  gene [0] = '\0';

  genome_view = FALSE;
  if (bsp->repr == Seq_repr_seg && (! SegHasParts (bsp))) {
    genome_view = TRUE;
  }
  if (bsp->repr == Seq_repr_delta && (! DeltaLitOnly (bsp))) {
    genome_view = TRUE;
  }

  /* locus id */

  sip = NULL;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_OTHER) {
      tsip = (TextSeqIdPtr) sip->data.ptrvalue;
      if (tsip != NULL) {
        if (StringNCmp (tsip->accession, "NM_", 3) == 0 ||
            StringNCmp (tsip->accession, "NR_", 3) == 0 ||
            StringNCmp (tsip->accession, "XM_", 3) == 0 ||
            StringNCmp (tsip->accession, "XR_", 3) == 0) {
          is_nm = TRUE;
          nm = bsp;
        } else if (StringNCmp (tsip->accession, "NP_", 3) == 0  ||
                   StringNCmp (tsip->accession, "XP_", 3) == 0) {
          is_np = TRUE;
        }
      }
      break;
    }
    if (sip->choice == SEQID_GENBANK ||
        sip->choice == SEQID_EMBL ||
        sip->choice == SEQID_DDBJ) break;
    if (sip->choice == SEQID_PIR ||
        sip->choice == SEQID_SWISSPROT ||
        sip->choice == SEQID_PRF ||
        sip->choice == SEQID_PDB) break;
    if (sip->choice == SEQID_TPG ||
        sip->choice == SEQID_TPE ||
        sip->choice == SEQID_TPD) break;
  }
  if (sip == NULL) {
    sip = SeqIdFindBest (bsp->id, SEQID_GENBANK);
  }

  if (genome_view) {
    SeqIdWrite (sip, locus, PRINTID_TEXTID_ACCESSION, sizeof (locus) - 1);
  } else {
    SeqIdWrite (sip, locus, PRINTID_TEXTID_LOCUS, sizeof (locus) - 1);
  }

  if (is_np) {
    sfp = SeqMgrGetCDSgivenProduct (bsp, &fcontext);
    if (sfp != NULL) {
      nm = fcontext.bsp;
      for (sip = nm->id; sip != NULL; sip = sip->next) {
        if (sip->choice == SEQID_OTHER) {
          tsip = (TextSeqIdPtr) sip->data.ptrvalue;
          if (tsip != NULL) {
            if (StringNCmp (tsip->accession, "NM_", 3) == 0 ||
                StringNCmp (tsip->accession, "XM_", 3) == 0) {
              is_nm = TRUE;
            }
          }
        }
      }
      if (! is_nm) {
        nm = NULL;
      }
    }
  }
  if (nm != NULL) {
    sfp = SeqMgrGetNextFeature (nm, NULL, SEQFEAT_GENE, 0, &fcontext);
    if (sfp != NULL) {
      StringNCpy_0 (gene, fcontext.label, sizeof (gene));
      if (SeqMgrGetNextFeature (nm, sfp, SEQFEAT_GENE, 0, &fcontext) != NULL) {
        gene [0] = '\0';
      }
      if (StringLen (gene) > 15) {
        gene [0] = '\0';
      }
    }
  }

  /* more complicated code to get parent locus, if segmented, goes here */

  if (awp->slp != NULL) {
    length = SeqLocLen (awp->slp);
  } else {
    length = bsp->length;
  }

  mip = NULL;
  tech = MI_TECH_standard;
  origin = 0;
  bmol = bsp->mol;
  if (bmol > Seq_mol_aa) {
    bmol = 0;
  }
  istrand = bsp->strand;
  if (istrand > 3) {
    istrand = 0;
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
  if (sdp != NULL) {
    bbp->entityID = dcontext.entityID;
    bbp->itemID = dcontext.itemID;
    bbp->itemtype = OBJ_SEQDESC;

    mip = (MolInfoPtr) sdp->data.ptrvalue;
    if (mip != NULL) {
      if (mip->biomol <= MOLECULE_TYPE_TRANSCRIBED_RNA) {
        imol = (Int2) mip->biomol;
      }
      tech = mip->tech;

      if (tech == MI_TECH_wgs && bsp->repr == Seq_repr_virtual) {

        /* check for WGS master record */

        for (sip = bsp->id; sip != NULL; sip = sip->next) {
          switch (sip->choice) {
            case SEQID_GENBANK :
            case SEQID_EMBL :
            case SEQID_DDBJ :
              tsip = (TextSeqIdPtr) sip->data.ptrvalue;
              if (tsip != NULL && tsip->accession != NULL) {
                if (StringLen (tsip->accession) == 12) {
                  if (StringCmp (tsip->accession + 6, "000000") == 0) {
                    wgsmaster = TRUE;
                  }
                }
              }
              break;
            case SEQID_OTHER :
              tsip = (TextSeqIdPtr) sip->data.ptrvalue;
              if (tsip != NULL && tsip->accession != NULL) {
                if (StringLen (tsip->accession) == 15) {
                  if (StringCmp (tsip->accession + 9, "000000") == 0) {
                    wgsmaster = TRUE;
                  }
                }
              }
              break;
            default :
              break;
          }
        }
      }
    }
  }

  /* check inst.mol if mol-type is not-set or genomic */

  if (imol <= MOLECULE_TYPE_GENOMIC) {
    if (bmol == Seq_mol_aa) {
      imol = MOLECULE_TYPE_PEPTIDE;
    } else if (bmol == Seq_mol_na) {
      imol = 0;
    } else if (bmol == Seq_mol_rna) {
      imol = 2;
    } else {
      imol = 1;
    }
  } else if (imol == MOLECULE_TYPE_OTHER_GENETIC_MATERIAL) {
    if (bmol == Seq_mol_rna) {
      imol = 2;
    }
  }

  /* if ds-DNA don't show ds */

  if (bmol == Seq_mol_dna && istrand == 2) {
    istrand = 0;
  }

  /* ss=any RNA don't show ss */

  if ((bmol > Seq_mol_rna ||
      (imol >= MOLECULE_TYPE_MRNA && imol <= MOLECULE_TYPE_PEPTIDE)) &&
      istrand == 1) {
    istrand = 0;
  }

  topology = bsp->topology;
  if (awp->slp != NULL) {
    topology = TOPOLOGY_LINEAR;
  }

  /* length, topology, and molecule type */

  if (awp->format == GENBANK_FMT) {

    if (awp->newLocusLine) {

      if (wgsmaster) {
        sprintf (len, "%ld rc", (long) length);
      } else {
        sprintf (len, "%ld bp", (long) length);
      }
      sprintf (mol, "%s%-4s", strd [istrand], gnbk_mol [imol]);

    } else {

      if (topology == TOPOLOGY_CIRCULAR) {
        sprintf (len, "%7ld bp", (long) length);
        sprintf (mol, "%s%-4s  circular", strd [istrand], gnbk_mol [imol]);
      } else {
        sprintf (len, "%7ld bp", (long) length);
        sprintf (mol, "%s%-4s          ", strd [istrand], gnbk_mol [imol]);
      }
    }

  } else if (awp->format == GENPEPT_FMT) {

    if (awp->newLocusLine) {
      sprintf (len, "%ld aa", (long) length);
    } else {
      sprintf (len, "%7ld aa", (long) length);
    }

  } else if (awp->format == EMBL_FMT) {

    if (imol < MOLECULE_TYPE_PEPTIDE) {
      if (ajp->flags.useEmblMolType) {
        ebmol = embl_mol [imol];
      } else {
        ebmol = gnbk_mol [imol];
      }

      if (topology == TOPOLOGY_CIRCULAR) {
        sprintf (mol, "circular %s", ebmol);
        sprintf (len, "%ld BP.", (long) length);
      } else {
        sprintf (mol, "%s", ebmol);
        sprintf (len, "%ld BP.", (long) length);
      }
    }
  }

  /* division */

  biop = NULL;
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp != NULL) {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
  } else {
    sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_BIOSRC, 0, &fcontext);
    if (sfp != NULL) {
      biop = (BioSourcePtr) sfp->data.value.ptrvalue;
    } else if (ISA_aa (bsp->mol)) {

      /* if protein with no sources, get sources applicable to DNA location of CDS */

      cds = SeqMgrGetCDSgivenProduct (bsp, &fcontext);
      if (cds != NULL) {
        sfp = SeqMgrGetOverlappingSource (cds->location, &fcontext);
        if (sfp != NULL) {
          biop = (BioSourcePtr) sfp->data.value.ptrvalue;
        } else {
          dna = BioseqFindFromSeqLoc (cds->location);
          if (dna != NULL) {
            sdp = SeqMgrGetNextDescriptor (dna, NULL, Seq_descr_source, &dcontext);
            if (sdp != NULL) {
              biop = (BioSourcePtr) sdp->data.ptrvalue;
            }
          }
        }
      }
    }
  }
  if (biop != NULL) {
    origin = biop->origin;
    orp = biop->org;
    if (orp != NULL) {
      onp = orp->orgname;
      if (onp != NULL) {
        StringNCpy_0 (div, onp->div, sizeof (div));
      }
    }
    for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
      if (ssp->subtype == SUBSRC_transgenic) {
        is_transgenic = TRUE;
      }
    }
  }

  switch (tech) {
    case MI_TECH_est :
      StringCpy (div, "EST");
      break;
    case MI_TECH_sts :
      StringCpy (div, "STS");
      break;
    case MI_TECH_survey :
      StringCpy (div, "GSS");
      break;
    case MI_TECH_htgs_0 :
    case MI_TECH_htgs_1 :
    case MI_TECH_htgs_2 :
      StringCpy (div, "HTG");
      break;
    case MI_TECH_htc :
      StringCpy (div, "HTC");
      break;
    default :
      break;
  }

  if (origin == 5 || is_transgenic) {
    StringCpy (div, "SYN");
  }

  sip = SeqIdFindBest (bsp->id, SEQID_PATENT);
  if (sip != NULL && sip->choice == SEQID_PATENT) {
    StringCpy (div, "PAT");
  }

  /* more complicated code for division, if necessary, goes here */

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_genbank, &dcontext);
  while (sdp != NULL) {
    gbp = (GBBlockPtr) sdp->data.ptrvalue;
    if (gbp != NULL) {
      if (StringHasNoText (div) && gbp->div != NULL) {
        StringCpy (div, gbp->div);
      }
      else if (StringCmp(gbp->div, "PAT") == 0 ||
               StringCmp(gbp->div, "SYN") == 0 ) {
        StringCpy (div, gbp->div);
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_genbank, &dcontext);
  }

  if (awp->format == EMBL_FMT || awp->format == EMBLPEPT_FMT) {

    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_embl, &dcontext);
    if (sdp != NULL) {
      ebp = (EMBLBlockPtr) sdp->data.ptrvalue;
      if (ebp != NULL) {
        if (ebp->div == 255) {
          if (mip == NULL) {
            StringCpy (div, "HUM");
          }
        } else if (ebp->div < 18)  {
          StringCpy (div, embl_divs [ebp->div]);
        }
      }
    }

    if (StringHasNoText (div)) {
      StringCpy (div, "UNA");
    }
  }

  /* empty division field if unable to find anything */

  if (StringHasNoText (div)) {
    StringCpy (div, "   ");
  }

  /* contig style (old genome_view flag) forces CON division */

  if (awp->contig) {
    StringCpy (div, "CON");
  }

  if (genome_view) {
    StringCpy (div, "CON");
  }

  /* date */

  best_date = GetBestDateForBsp (bsp);

  if (best_date == NULL) {

    /* if bsp is product of CDS or mRNA feature, get date from sfp->location bsp */

    sfp = NULL;
    if (ISA_na (bsp->mol)) {
      sfp = SeqMgrGetRNAgivenProduct (bsp, NULL);
    } else if (ISA_aa (bsp->mol)) {
      sfp = SeqMgrGetCDSgivenProduct (bsp, NULL);
    }
    if (sfp != NULL) {
      parent = BioseqFindFromSeqLoc (sfp->location);
      if (parent != NULL) {
        best_date = GetBestDateForBsp (parent);
      }
    }
  }

  /* convert best date */

  if (best_date != NULL) {
    DateToGB (date, best_date, FALSE);
  }
  if (StringHasNoText (date)) {
    StringCpy (date, "01-JAN-1900");
  }

  if (awp->format == GENBANK_FMT || awp->format == GENPEPT_FMT) {

    /* Create the proper locus name */

    parent = awp->parent;
    if (parent->repr == Seq_repr_seg) {

      if (! StringHasNoText (awp->basename)) {
        StringCpy (locus, awp->basename);
        s_LocusAddSuffix (locus, awp);
      }
    }

    /* Print the "LOCUS_NEW" line, if requested */

    if (awp->newLocusLine) {
      FFStartPrint (ffstring, awp->format, 0, 0, "LOCUS", 12, 0, 0, NULL, FALSE);
      parent = awp->parent;

      if (parent->repr == Seq_repr_seg)
        s_LocusAdjustLength (locus,16);

      if (is_nm && (! StringHasNoText (gene))) {
        FFAddOneString (ffstring, gene, FALSE, FALSE, TILDE_IGNORE);
      } else {
        FFAddOneString (ffstring, locus, FALSE, FALSE, TILDE_IGNORE);
      }
      FFAddNChar(ffstring, ' ', 43 - StringLen(len)- ffstring->curr->pos, FALSE);
      FFAddOneString (ffstring, len, FALSE, FALSE, TILDE_IGNORE);
      FFAddNChar(ffstring, ' ', 44 - ffstring->curr->pos, FALSE);
      FFAddOneString (ffstring, mol, FALSE, FALSE, TILDE_IGNORE);
      FFAddNChar(ffstring, ' ', 55 - ffstring->curr->pos, FALSE);
      if (topology == TOPOLOGY_CIRCULAR) {
        FFAddOneString (ffstring, "circular", FALSE, FALSE, TILDE_IGNORE);
      } else {
        FFAddOneString (ffstring, "linear  ", FALSE, FALSE, TILDE_IGNORE);
      }
      FFAddNChar(ffstring, ' ', 64 - ffstring->curr->pos, FALSE);
      FFAddOneString (ffstring, div, FALSE, FALSE, TILDE_IGNORE);
      FFAddNChar(ffstring, ' ', 68 - ffstring->curr->pos, FALSE);
      FFAddOneString (ffstring, date, FALSE, FALSE, TILDE_IGNORE);
    }

    /* Else print the "LOCUS" line */

    else {
      FFStartPrint (ffstring, awp->format, 0, 0, "LOCUS", 12, 0, 0, NULL, FALSE);

      if (parent->repr == Seq_repr_seg)
        s_LocusAdjustLength (locus,16);

      FFAddOneString (ffstring, locus, FALSE, FALSE, TILDE_IGNORE);
      FFAddNChar(ffstring, ' ', 32 - StringLen(len) - ffstring->curr->pos, FALSE);
      FFAddOneString (ffstring, len, FALSE, FALSE, TILDE_IGNORE);
      FFAddNChar(ffstring, ' ', 33 - ffstring->curr->pos, FALSE);
      FFAddOneString (ffstring, mol, FALSE, FALSE, TILDE_IGNORE);
      FFAddNChar(ffstring, ' ', 52 - ffstring->curr->pos, FALSE);
      FFAddOneString (ffstring, div, FALSE, FALSE, TILDE_IGNORE);
      FFAddNChar(ffstring, ' ', 62 - ffstring->curr->pos, FALSE);
      FFAddOneString (ffstring, date, FALSE, FALSE, TILDE_IGNORE);
    }

  } else if (awp->format == EMBL_FMT || awp->format == EMBLPEPT_FMT) {
    FFStartPrint (ffstring, awp->format, 0, 0, NULL, 0, 5, 0, "ID", FALSE);
    FFAddOneString (ffstring, locus, FALSE, FALSE, TILDE_IGNORE);
    FFAddNChar(ffstring, ' ', 15 - 5 - StringLen(locus), FALSE);
    if (awp->hup) {
      FFAddOneString (ffstring, " confidential; ", FALSE, FALSE, TILDE_IGNORE);
    } else {
      FFAddOneString (ffstring, " standard; ", FALSE, FALSE, TILDE_IGNORE);
    }
    FFAddOneString (ffstring, mol, FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString (ffstring, "; ", FALSE, FALSE, TILDE_IGNORE);

    /* conditional code to make div "UNA" goes here */

    FFAddOneString (ffstring, div, FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString (ffstring, "; ", FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString (ffstring, len, FALSE, FALSE, TILDE_IGNORE);
  }

  /* optionally populate indexes for NCBI internal database */

  if (ajp->index) {
    index = &asp->index;
  } else {
    index = NULL;
  }

  if (index != NULL) {
    Char  tmp [20];
    index->locus = StringSave (locus);
    index->div = StringSave (div);
    sprintf (tmp, "%ld", (long) length);
    index->base_cnt = StringSave (tmp);
  }

  /* optionally populate gbseq for XML-ized GenBank format */

  if (ajp->gbseq) {
    gbseq = &asp->gbseq;
  } else {
    gbseq = NULL;
  }

  if (gbseq != NULL) {
    gbseq->locus = StringSave (locus);
    gbseq->length = length;
    gbseq->division = StringSave (div);
    gbseq->strandedness = istrand;
    gbseq->moltype = imolToMoltype [imol];
    gbseq->topology = topology;

    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      SeqIdWrite (sip, id, PRINTID_FASTA_SHORT, sizeof (id));
      ValNodeCopyStr (&gbseq->other_seqids, 0, id);
    }

    date [0] = '\0';
    dp = NULL;
    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_create_date, &dcontext);
    if (sdp != NULL) {
      dp = (DatePtr) sdp->data.ptrvalue;
    }
    if (dp != NULL) {
      DateToGB (date, dp, FALSE);
    }
    if (StringHasNoText (date)) {
      StringCpy (date, "01-JAN-1900");
    }
    gbseq->create_date = StringSave (date);

    date [0] = '\0';
    dp = NULL;
    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_update_date, &dcontext);
    if (sdp != NULL) {
      dp = (DatePtr) sdp->data.ptrvalue;
    }
    if (dp != NULL) {
      DateToGB (date, dp, FALSE);
    }
    if (StringHasNoText (date)) {
      StringCpy (date, "01-JAN-1900");
    }
    gbseq->update_date = StringSave (date);
  }

  bbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 0, 5, 0, "ID");
  FFRecycleString(ajp, ffstring);
}

static void AddDeflineBlock (
  Asn2gbWorkPtr awp
)

{
  IntAsn2gbJobPtr    ajp;
  Asn2gbSectPtr      asp;
  BaseBlockPtr       bbp;
  BioseqPtr          bsp;
  Char               buf[1024]; 
  /*CharPtr          buf;
  size_t             buflen = 1024;*/ 
  SeqMgrDescContext  dcontext;
  GBSeqPtr           gbseq;
  ItemInfo           ii;
  MolInfoPtr         mip;
  SeqDescrPtr        sdp;
  Uint1              tech;
  StringItemPtr      ffstring;

  if (awp == NULL) return;
  ajp = awp->ajp;
  if (ajp == NULL) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;
  asp = awp->asp;
  if (asp == NULL) return;

  bbp = Asn2gbAddBlock (awp, DEFLINE_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  tech = 0;
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
  if (sdp != NULL) {
    mip = (MolInfoPtr) sdp->data.ptrvalue;
    if (mip != NULL) {
      tech = mip->tech;
    }
  }

  ffstring = FFGetString(ajp);
  if ( ffstring == NULL ) return;

  /*buf = MemNew (sizeof (Char) * (buflen + 1));*/
  MemSet ((Pointer) (&ii), 0, sizeof (ItemInfo));
  MemSet ((Pointer) buf, 0, sizeof (buf));

  /* create default defline */

  if ( CreateDefLine (&ii, bsp, buf, sizeof(buf), tech, NULL, NULL)) {
    bbp->entityID = ii.entityID;
    bbp->itemID = ii.itemID;
    bbp->itemtype = ii.itemtype;

    FFStartPrint (ffstring, awp->format, 0, 12, "DEFINITION", 12, 5, 5, "DE", TRUE);

    if (StringHasNoText (buf)) {
      FFAddOneChar (ffstring, '.', FALSE);
    } else {
      FFAddOneString (ffstring, buf, TRUE, TRUE, TILDE_IGNORE);
    }

    bbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "DE");
  }

  /* optionally populate gbseq for XML-ized GenBank format */

  if (ajp->gbseq) {
    gbseq = &asp->gbseq;
  } else {
    gbseq = NULL;
  }

  if (gbseq != NULL) {
    gbseq->definition = StringSave (buf);
  }

  FFRecycleString(ajp, ffstring);
}

/* !!! this definitely needs more work to support all classes, use proper SeqId !!! */

static void AddAccessionBlock (
  Asn2gbWorkPtr awp
)

{
  SeqIdPtr           accn = NULL;
  IntAsn2gbJobPtr    ajp;
  Asn2gbSectPtr      asp;
  BaseBlockPtr       bbp;
  BioseqPtr          bsp;
  Char               buf [41];
  SeqMgrDescContext  dcontext;
  EMBLBlockPtr       ebp;
  ValNodePtr         extra_access;
  CharPtr            flatloc;
  GBBlockPtr         gbp;
  SeqIdPtr           gi = NULL;
  GBSeqPtr           gbseq;
  IndxPtr            index;
  SeqIdPtr           lcl = NULL;
  size_t             len = 0;
  MolInfoPtr         mip;
  SeqDescrPtr        sdp;
  CharPtr            separator = " ";
  SeqIdPtr           sip;
  TextSeqIdPtr       tsip;
  ValNodePtr         vnp;
  CharPtr            wgsaccn = NULL;
  CharPtr            xtra;
  StringItemPtr      ffstring;

  if (awp == NULL) return;
  ajp = awp->ajp;
  if (ajp == NULL) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;
  asp = awp->asp;
  if (asp == NULL) return;
  
  ffstring = FFGetString(ajp);
  if ( ffstring == NULL ) return;

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    switch (sip->choice) {
      case SEQID_GI :
        gi = sip;
        break;
      case SEQID_GENBANK :
      case SEQID_EMBL :
      case SEQID_DDBJ :
        accn = sip;
        tsip = (TextSeqIdPtr) sip->data.ptrvalue;
        if (tsip != NULL) {
          if (StringLen (tsip->accession) == 12) {
            wgsaccn = tsip->accession;
            len = 12;
          }
        }
        break;
      case SEQID_OTHER :
        accn = sip;
        tsip = (TextSeqIdPtr) sip->data.ptrvalue;
        if (tsip != NULL) {
          if (StringLen (tsip->accession) == 15) {
            wgsaccn = tsip->accession;
            len = 15;
          }
        }
        break;
      case SEQID_PIR :
      case SEQID_SWISSPROT :
      case SEQID_PRF :
      case SEQID_PDB :
        accn = sip;
        break;
      case SEQID_TPG :
      case SEQID_TPE :
      case SEQID_TPD :
        accn = sip;
        break;
      case SEQID_GENERAL :
        /* should not override better accession */
        if (accn == NULL) {
          accn = sip;
        }
        break;
      case SEQID_LOCAL :
        lcl = sip;
        break;
      default :
        break;
    }
  }

  sip = NULL;
  if (accn != NULL) {
    sip = accn;
  } else if (lcl != NULL) {
    sip = lcl;
  } else if (gi != NULL) {
    sip = gi;
  }

  if (sip == NULL) return;

  SeqIdWrite (sip, buf, PRINTID_TEXTID_ACC_ONLY, sizeof (buf));

  bbp = Asn2gbAddBlock (awp, ACCESSION_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  FFStartPrint (ffstring, awp->format, 0, 12, "ACCESSION", 12, 5, 5, "AC", TRUE);

  if (awp->hup && accn != NULL) {
    FFAddOneString (ffstring, ";", FALSE, FALSE, TILDE_TO_SPACES);

  } else if (ajp->ajp.slp != NULL) {
    FF_www_accession (ajp, ffstring, buf);
    flatloc =  FlatLoc (ajp, bsp, ajp->ajp.slp, ajp->masterStyle);
    FFAddTextToString (ffstring, " REGION: ", flatloc, NULL, FALSE, FALSE, TILDE_TO_SPACES);
    MemFree (flatloc);
  } else {
    FFAddOneString (ffstring, buf, FALSE, FALSE, TILDE_TO_SPACES);
  }

  /* optionally populate indexes for NCBI internal database */

  if (ajp->index) {
    index = &asp->index;
  } else {
    index = NULL;
  }

  if (index != NULL) {
    index->accession = StringSave (buf);
  }

  /* optionally populate gbseq for XML-ized GenBank format */

  if (ajp->gbseq) {
    gbseq = &asp->gbseq;
  } else {
    gbseq = NULL;
  }

  if (gbseq != NULL) {
    gbseq->primary_accession = StringSave (buf);
  }

  if (awp->format == GENBANK_FMT || awp->format == GENPEPT_FMT) {
    separator = " ";
  } else if (awp->format == EMBL_FMT || awp->format == EMBLPEPT_FMT) {
    separator = ";";
  }

  if (ajp->ajp.slp == NULL) {
    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
    if (sdp != NULL && wgsaccn != NULL) {
      mip = (MolInfoPtr) sdp->data.ptrvalue;
      if (mip != NULL && mip->tech == MI_TECH_wgs) {
        StringNCpy_0 (buf, wgsaccn, sizeof (buf));
        if (StringCmp (buf + len - 6, "000000") != 0) {
          StringCpy (buf + len - 6, "000000");
          FFAddTextToString(ffstring, separator, buf, NULL, FALSE, FALSE, TILDE_TO_SPACES);
        } else if (StringCmp (buf + len - 8, "00000000") != 0) {
          StringCpy (buf + len - 8, "00000000");
          FFAddTextToString(ffstring, separator, buf, NULL, FALSE, FALSE, TILDE_TO_SPACES);
        }
      }
    }

    sdp = SeqMgrGetNextDescriptor (bsp, NULL, 0, &dcontext);
    while (sdp != NULL) {

      extra_access = NULL;

      switch (dcontext.seqdesctype) {
        case Seq_descr_genbank :
          gbp = (GBBlockPtr) sdp->data.ptrvalue;
          if (gbp != NULL) {
            extra_access = gbp->extra_accessions;
          }
          break;
        case Seq_descr_embl :
          ebp = (EMBLBlockPtr) sdp->data.ptrvalue;
          if (ebp != NULL) {
            extra_access = ebp->extra_acc;
          }
          break;
        default :
          break;
      }

      if (extra_access != NULL) {
        bbp->entityID = dcontext.entityID;
        bbp->itemID = dcontext.itemID;
        bbp->itemtype = OBJ_SEQDESC;
      }

      for (vnp = extra_access; vnp != NULL; vnp = vnp->next) {
        xtra = (CharPtr) vnp->data.ptrvalue;
        if (ValidateAccn (xtra) == 0) {
          FFAddTextToString(ffstring, separator, xtra, NULL, FALSE, FALSE, TILDE_TO_SPACES);

          /* optionally populate indexes for NCBI internal database */

          if (index != NULL) {
            ValNodeCopyStrToHead (&(index->secondaries), 0, xtra);
          }

          /* optionally populate gbseq for XML-ized GenBank format */

          if (gbseq != NULL) {
            ValNodeCopyStr (&(gbseq->secondary_accessions), 0, xtra);
          }
        }
      }

      sdp = SeqMgrGetNextDescriptor (bsp, sdp, 0, &dcontext);
    }
  }

  bbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "AC");
  FFRecycleString(ajp, ffstring);
}

static void AddVersionBlock (
  Asn2gbWorkPtr awp
)

{
  SeqIdPtr         accn = NULL;
  IntAsn2gbJobPtr  ajp;
  Asn2gbSectPtr    asp;
  BaseBlockPtr     bbp;
  BioseqPtr        bsp;
  Char             buf [41];
  GBSeqPtr         gbseq;
  Int4             gi = -1;
  IndxPtr          index;
  CharPtr          ptr;
  SeqIdPtr         sip;
  Char             tmp [41];
  Char             version [64];
  StringItemPtr    ffstring;

  if (awp == NULL) return;
  ajp = awp->ajp;
  if (ajp == NULL) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;
  asp = awp->asp;
  if (asp == NULL) return;

  ffstring = FFGetString(ajp);
  if ( ffstring == NULL ) return;

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    switch (sip->choice) {
      case SEQID_GI :
        gi = sip->data.intvalue;
        break;
      case SEQID_GENBANK :
      case SEQID_EMBL :
      case SEQID_DDBJ :
      case SEQID_OTHER :
        accn = sip;
        break;
      case SEQID_PIR :
      case SEQID_SWISSPROT :
      case SEQID_PRF :
      case SEQID_PDB :
        accn = sip;
        break;
      case SEQID_TPG :
      case SEQID_TPE :
      case SEQID_TPD :
        accn = sip;
        break;
      default :
        break;
    }
  }

  /* if (gi < 1 && accn == NULL) return; */

  if (awp->format == EMBL_FMT || awp->format == EMBLPEPT_FMT) {
    if (gi < 1) return;
  }

  bbp = Asn2gbAddBlock (awp, VERSION_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  /* no longer displaying NID */

  /*
  if (gi > 0) {
    sprintf (version, "g%ld", (long) gi);

    gb_StartPrint (awp->format, needInitBuff, 0, 12, "NID", 13, 5, 5, "NI", TRUE);
    needInitBuff = FALSE;

    gb_AddString (NULL, version, NULL, FALSE, FALSE, TILDE_TO_SPACES);

    ff_EndPrint();
    needEndPrint = FALSE;
  }
  */

  version [0] = '\0';

  if (awp->format == EMBL_FMT || awp->format == EMBLPEPT_FMT) {
    sprintf (version, "g%ld", (long) gi);

    FFStartPrint (ffstring, awp->format, 0, 12, "VERSION", 12, 5, 5, "NI", TRUE);

    FFAddOneString (ffstring, version, FALSE, FALSE, TILDE_TO_SPACES);

    FFAddOneChar(ffstring, '\n', FALSE);

    bbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "NI");
    FFRecycleString(ajp, ffstring);
    return;
  }

  if (accn != NULL) {

    buf [0] = '\0';
    SeqIdWrite (accn, buf, PRINTID_TEXTID_ACC_VER, sizeof (buf) - 1);

    if (gi > 0) {
      sprintf (version, "%s  GI:%ld", buf, (long) gi);
    } else {
      sprintf (version, "%s", buf);
    }

    FFStartPrint (ffstring, awp->format, 0, 12, "VERSION", 12, 5, 5, "SV", TRUE);

    FFAddTextToString (ffstring, NULL, version, "\n", FALSE, FALSE, TILDE_TO_SPACES);
    /* optionally populate indexes for NCBI internal database */

    if (ajp->index) {
      index = &asp->index;
    } else {
      index = NULL;
    }

    if (index != NULL) {
      ptr = StringChr (buf, '.');
      if (ptr != NULL) {
        ptr++;
        index->version = StringSave (ptr);
      }
      if (gi > 0) {
        sprintf (tmp, "%ld", (long) gi);
        index->gi = StringSave (tmp);
      }
    }

    /* optionally populate gbseq for XML-ized GenBank format */

    if (ajp->gbseq) {
      gbseq = &asp->gbseq;
    } else {
      gbseq = NULL;
    }

    if (gbseq != NULL) {
      ptr = StringChr (buf, '.');
      if (ptr != NULL) {
        gbseq->accession_version = StringSave (buf);
      }
    }

  } else if (gi > 0) {

    FFStartPrint (ffstring, awp->format, 0, 0, "VERSION", 12, 5, 5, "SV", TRUE);

    sprintf (version, "  GI:%ld", (long) gi);

    FFAddTextToString (ffstring, NULL, version, "\n", FALSE, FALSE, TILDE_TO_SPACES);

  } else {

    FFStartPrint (ffstring, awp->format, 0, 0, "VERSION", 0, 5, 5, "SV", TRUE);
    FFAddOneChar(ffstring, '\n', FALSE);
    bbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "SV");
  }
  bbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "SV");
  FFRecycleString(ajp, ffstring);
}

/* only displaying PID in GenPept format */

static void AddPidBlock (Asn2gbWorkPtr awp)

{
  IntAsn2gbJobPtr  ajp;
  BaseBlockPtr  bbp;
  BioseqPtr     bsp;
  Int4          gi = -1;
  SeqIdPtr      sip;
  Char          version [64];
  StringItemPtr ffstring;

  if (awp == NULL) return;
  ajp = awp->ajp;
  if (ajp == NULL) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    switch (sip->choice) {
      case SEQID_GI :
        gi = sip->data.intvalue;
        break;
      default :
        break;
    }
  }

  if (gi < 1) return;

  bbp = Asn2gbAddBlock (awp, PID_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  ffstring = FFGetString(ajp);
  if ( ffstring == NULL ) return;

  FFStartPrint (ffstring, awp->format, 0, 12, "PID", 12, 5, 5, NULL, TRUE);

  sprintf (version, "g%ld", (long) gi);
  FFAddOneString (ffstring, version, FALSE, FALSE, TILDE_TO_SPACES);

  bbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, NULL);
  FFRecycleString(ajp, ffstring);
}

static Uint1 dbsource_fasta_order [NUM_SEQID] = {
  33, /* 0 = not set */
  20, /* 1 = local Object-id */
  15, /* 2 = gibbsq */
  16, /* 3 = gibbmt */
  30, /* 4 = giim Giimport-id */
  10, /* 5 = genbank */
  10, /* 6 = embl */
  10, /* 7 = pir */
  10, /* 8 = swissprot */
  15, /* 9 = patent */
  18, /* 10 = other TextSeqId */
  20, /* 11 = general Dbtag */
  31, /* 12 = gi */
  10, /* 13 = ddbj */
  10, /* 14 = prf */
  12, /* 15 = pdb */
  10, /* 16 = tpg */
  10, /* 17 = tpe */
  10  /* 18 = tpd */
};

static void AddToUniqueSipList (
  ValNodePtr PNTR list,
  SeqIdPtr sip
)

{
  ValNodePtr  vnp;

  if (list == NULL || sip == NULL) return;
  for (vnp = *list; vnp != NULL; vnp = vnp->next) {
    if (SeqIdMatch (sip, (SeqIdPtr) vnp->data.ptrvalue)) return;
  }
  ValNodeAddPointer (list, 0, (Pointer) sip);
}

static Boolean WriteDbsourceID (
  SeqIdPtr sip,
  CharPtr str
)

{
  DbtagPtr      db;
  CharPtr       dt;
  Int4          gi;
  ObjectIdPtr   oip;
  CharPtr       pfx;
  PDBSeqIdPtr   psip = NULL;
  CharPtr       prefix;
  Boolean       rsult = FALSE;
  CharPtr       sfx;
  CharPtr       suffix;
  Char          tmp [32];
  TextSeqIdPtr  tsip = NULL;

  if (sip == NULL || str == NULL) return FALSE;
  *str = '\0';
  switch (sip->choice) {
    case SEQID_LOCAL :
      oip = (ObjectIdPtr) sip->data.ptrvalue;
      if (oip == NULL) return FALSE;
      if (! StringHasNoText (oip->str)) {
        StringCat (str, oip->str);
        return TRUE;
      } else if (oip->id > 0) {
        sprintf (tmp, "%ld", (long) oip->id);
        StringCat (str, tmp);
        return TRUE;
      }
      return FALSE;
    case SEQID_GI :
      gi = (Int4) sip->data.intvalue;
      if (gi == 0) return FALSE;
      sprintf (tmp, "gi: %ld", (long) gi);
      StringCat (str, tmp);
      return TRUE;
    case SEQID_GENERAL :
      db = (DbtagPtr) sip->data.ptrvalue;
      if (db == NULL) return FALSE;
      /* !!! still need to implement this !!! */
      return FALSE;
    case SEQID_GENBANK :
    case SEQID_EMBL :
    case SEQID_DDBJ :
    case SEQID_OTHER :
    case SEQID_PIR :
    case SEQID_SWISSPROT :
    case SEQID_PRF :
    case SEQID_TPG :
    case SEQID_TPE :
    case SEQID_TPD :
      tsip = (TextSeqIdPtr) sip->data.ptrvalue;
      if (tsip == NULL) return FALSE;
      break;
    case SEQID_PDB :
      psip = (PDBSeqIdPtr) sip->data.ptrvalue;
      if (psip == NULL) return FALSE;
      break;
    default :
      break;
  }
  prefix = " ";
  suffix = NULL;
  switch (sip->choice) {
    case SEQID_EMBL :
      StringCat (str, "embl ");
      suffix = ",";
      break;
    case SEQID_OTHER :
      StringCat (str, "REFSEQ: ");
      break;
    case SEQID_SWISSPROT :
      StringCat (str, "swissprot: ");
      suffix = ",";
      break;
    case SEQID_PIR :
      StringCat (str, "pir: ");
      break;
    case SEQID_PRF :
      StringCat (str, "prf: ");
      break;
    case SEQID_PDB :
      StringCat (str, "pdb: ");
      suffix = ",";
      break;
    default :
      break;
  }
  pfx = NULL;
  sfx = NULL;
  if (tsip != NULL) {
    if (! StringHasNoText (tsip->name)) {
      StringCat (str, sfx);
      StringCat (str, pfx);
      StringCat (str, "locus ");
      StringCat (str, tsip->name);
      sfx = suffix;
      pfx = prefix;
      rsult = TRUE;
    }
    if (! StringHasNoText (tsip->accession)) {
      StringCat (str, sfx);
      StringCat (str, pfx);
      StringCat (str, "accession ");
      StringCat (str, tsip->accession);
      sfx = suffix;
      pfx = prefix;
      rsult = TRUE;
    }
    if (tsip->version > 0) {
      sprintf (tmp, ".%d", (int) tsip->version);
      StringCat (str, tmp);
      sfx = suffix;
      pfx = prefix;
    }
    if (! StringHasNoText (tsip->release)) {
      StringCat (str, pfx);
      StringCat (str, "release ");
      StringCat (str, tsip->release);
      sfx = suffix;
      pfx = prefix;
    }
    if (sip->choice == SEQID_SWISSPROT || sip->choice == SEQID_PIR || sip->choice == SEQID_PRF) {
      StringCat (str, ";");
    }
    return rsult;
  }
  if (psip != NULL) {
    if (! StringHasNoText (psip->mol)) {
      StringCat (str, "molecule ");
      StringCat (str, psip->mol);
      sfx = suffix;
      pfx = prefix;
      rsult = TRUE;
    }
    if (psip->chain > 0) {
      StringCat (str, sfx);
      StringCat (str, pfx);
      sprintf (tmp, "chain %d", (int) psip->chain);
      StringCat (str, tmp);
      sfx = suffix;
      pfx = prefix;
      rsult = TRUE;
    }
    if (psip->rel != NULL) {
      StringCat (str, sfx);
      StringCat (str, pfx);
      StringCat (str, "release ");
      dt = asn2gb_PrintDate (psip->rel);
      StringCat (str, dt);
      MemFree (dt);
      sfx = suffix;
      pfx = prefix;
      rsult = TRUE;
    }
    StringCat (str, ";");
    return rsult;
  }
  return rsult;
}


static void AddSPBlock (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  BioseqPtr bsp
)

{
  CharPtr            acc;
  DbtagPtr           db;
  SeqMgrDescContext  dcontext;
  Boolean            first;
  Boolean            has_link;
  Char               id [41];
  ObjectIdPtr        oip;
  SeqDescrPtr        sdp;
  SeqIdPtr           sid;
  SPBlockPtr         spb;
  CharPtr            string;
  ValNodePtr         vnp;
  CharPtr            str;
  Char               numbuf[40];

  if (bsp == NULL) return;
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_sp, &dcontext);
  if (sdp == NULL) return;
  spb = (SPBlockPtr) sdp->data.ptrvalue;
  if (spb == NULL) return;

  if (spb->_class == 1) {
    FFAddOneString (ffstring, "class: standard.", FALSE, FALSE, TILDE_IGNORE);
    FFAddNewLine(ffstring);
  } else if (spb->_class == 2) {
    FFAddOneString (ffstring, "class: preliminary.", FALSE, FALSE, TILDE_IGNORE);
    FFAddNewLine(ffstring);
  }

  if (spb->extra_acc) {
    FFAddOneString (ffstring, "extra accessions:", FALSE, FALSE, TILDE_IGNORE);
    for (vnp = spb->extra_acc; vnp != NULL; vnp = vnp->next) {
      FFAddOneString (ffstring, (CharPtr) vnp->data.ptrvalue, FALSE, FALSE, TILDE_IGNORE);
      FFAddOneChar (ffstring, ',', FALSE );
    }
  }

  if (spb->imeth) {
    FFAddOneString (ffstring, "seq starts with Met", FALSE, FALSE, TILDE_IGNORE);
  }

  if (spb->plasnm != NULL) {
    FFAddOneString (ffstring, "plasmid:", FALSE, FALSE, TILDE_IGNORE);
    for (vnp = spb->plasnm; vnp != NULL; vnp = vnp->next) {
      FFAddOneString (ffstring, (CharPtr) vnp->data.ptrvalue, FALSE, FALSE, TILDE_IGNORE);
      FFAddOneChar (ffstring, ',', FALSE );
    }
  }

  if (spb->created) {
    string = PrintDate (spb->created);
    FFAddOneString (ffstring, "created: ", FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString (ffstring, string, FALSE, FALSE, TILDE_IGNORE);

    MemFree (string);
  }

  if (spb->sequpd) {
    string = PrintDate (spb->sequpd);
    FFAddOneString (ffstring, "sequence updated: ", FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString (ffstring, string, FALSE, FALSE, TILDE_IGNORE);
    MemFree (string);
  }

  if (spb->annotupd) {
    string = PrintDate (spb->annotupd);
    FFAddOneString (ffstring, "annotation updated: ", FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString (ffstring, string, FALSE, FALSE, TILDE_IGNORE);
    MemFree (string);
  }

  if (spb->seqref) {
    FFAddOneString (ffstring, "xrefs: ", FALSE, FALSE, TILDE_IGNORE);
    first = TRUE;
    for (sid = spb->seqref; sid != NULL; sid = sid->next) {
      acc = NULL;
      has_link = FALSE;
      if (first == FALSE) {
        FFAddOneString (ffstring, ", ", FALSE, FALSE, TILDE_IGNORE);
      }
      first = FALSE;
      SeqIdWrite (sid, id, PRINTID_TEXTID_ACC_VER, sizeof (id) - 1);
      if (sid->choice == SEQID_GI) {
        has_link = TRUE;
      }
      acc = id;
      if (acc != NULL) {
        switch (sid->choice) {
          case SEQID_GENBANK:
            FFAddOneString (ffstring, "genbank accession ", FALSE, FALSE, TILDE_IGNORE);
            break; 
          case SEQID_EMBL:
            FFAddOneString (ffstring, "embl accession ", FALSE, FALSE, TILDE_IGNORE);
            break; 
          case SEQID_PIR:
            FFAddOneString (ffstring, "pir locus ", FALSE, FALSE, TILDE_IGNORE);
            break; 
          case SEQID_SWISSPROT:
            FFAddOneString (ffstring, "swissprot accession ", FALSE, FALSE, TILDE_IGNORE);
            break; 
          case SEQID_DDBJ:
            FFAddOneString (ffstring, "ddbj accession ", FALSE, FALSE, TILDE_IGNORE);
            break; 
          case SEQID_PRF:
            FFAddOneString (ffstring, "prf accession ", FALSE, FALSE, TILDE_IGNORE);
            break;
          case SEQID_PDB:
            FFAddOneString (ffstring, "pdb accession ", FALSE, FALSE, TILDE_IGNORE);
            break;
          case SEQID_GI:
            FFAddOneString (ffstring, "gi: ", FALSE, FALSE, TILDE_IGNORE);
            break; 
          case SEQID_TPG:
            FFAddOneString (ffstring, "genbank third party accession ", FALSE, FALSE, TILDE_IGNORE);
            break; 
          case SEQID_TPE:
            FFAddOneString (ffstring, "embl third party accession ", FALSE, FALSE, TILDE_IGNORE);
            break; 
          case SEQID_TPD:
            FFAddOneString (ffstring, "ddbj third party accession ", FALSE, FALSE, TILDE_IGNORE);
            break; 
          default:
            acc = NULL;
            break; 
        }
      }
      if (acc != NULL) {
        if ( GetWWW(ajp) && has_link ) {
          sprintf(numbuf, "%ld", (long) sid->data.intvalue);
          FFAddTextToString(ffstring, "<a href=", link_seq, NULL, FALSE, FALSE, TILDE_IGNORE);
          FFAddTextToString(ffstring, "val=", numbuf, ">", FALSE, FALSE, TILDE_IGNORE);
          FFAddOneString(ffstring, acc, FALSE, FALSE, TILDE_IGNORE);
          FFAddOneString(ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
        } else {
          FFAddOneString(ffstring, acc, FALSE, FALSE, TILDE_IGNORE);
        }
      }
    }
  }

  first = TRUE;
  for (vnp = spb->dbref; vnp != NULL; vnp = vnp->next) {
    db = (DbtagPtr) vnp->data.ptrvalue;
    if (db == NULL) continue;
    oip = db->tag;
    if (oip == NULL) continue;
    has_link = FALSE;
    if (first) {
      FFAddNewLine(ffstring);
      FFAddOneString (ffstring, "xrefs (non-sequence databases): ", FALSE, FALSE, TILDE_IGNORE);
      first = FALSE;
    } else {
      FFAddOneString (ffstring, ", ", FALSE, FALSE, TILDE_IGNORE);
    }
    FFAddOneString (ffstring, db->db, FALSE, FALSE, TILDE_IGNORE);
    if (StringCmp (db->db, "MIM") == 0) {
      has_link = TRUE;
    }

    if ( oip->str != NULL ) {
      str = oip->str;
    } else if ( oip->id > 0 ) {
      sprintf(numbuf, "%d", oip->id);
      str = numbuf;
    }

    if ( !StringHasNoText(str) ) {
      if ( GetWWW(ajp) && has_link) {
        FFAddOneChar (ffstring, ' ', FALSE);
        FFAddTextToString(ffstring, "<a href=", link_omim, str, FALSE, FALSE, TILDE_IGNORE);
        FFAddTextToString(ffstring, ">", str, "</a>", FALSE, FALSE, TILDE_IGNORE);
      } else {
        FFAddOneString(ffstring, str, FALSE, FALSE, TILDE_IGNORE);
      }
    }
  }
}

static void AddPIRBlock (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  BioseqPtr bsp
)

{
  CharPtr            acc;
  SeqMgrDescContext  dcontext;
  Boolean            first;
  Char               id [41];
  CharPtr            prefix = NULL;
  SeqDescrPtr        sdp;
  SeqIdPtr           sid;
  PirBlockPtr        pbp;

  if (bsp == NULL) return;
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_pir, &dcontext);
  if (sdp == NULL) return;
  pbp = (PirBlockPtr) sdp->data.ptrvalue;
  if (pbp == NULL) return;

  if (pbp->host != NULL) {
    FFAddTextToString (ffstring, "host:", pbp->host, "\n", FALSE, TRUE, TILDE_IGNORE);
    prefix = ";";
  }

  if (pbp->source != NULL) {
    FFAddOneString (ffstring, prefix, FALSE, FALSE, TILDE_IGNORE);
    FFAddNewLine(ffstring);
    FFAddTextToString(ffstring, "source: ", pbp->source, "\n", FALSE, TRUE, TILDE_IGNORE);
    prefix = ";";
  }

  if (pbp->summary != NULL) {
    FFAddOneString (ffstring, prefix, FALSE, FALSE, TILDE_IGNORE);
    FFAddNewLine(ffstring);
    FFAddTextToString(ffstring, "summary: ", pbp->summary, "\n", FALSE, TRUE, TILDE_IGNORE);
    prefix = ";";
  }

  if (pbp->genetic != NULL) {
    FFAddOneString (ffstring, prefix, FALSE, FALSE, TILDE_IGNORE);
    FFAddNewLine(ffstring);
    FFAddTextToString(ffstring, "genetic: ", pbp->genetic, "\n", FALSE, TRUE, TILDE_IGNORE);
    prefix = ";";
  }

  if (pbp->includes != NULL) {
    FFAddOneString (ffstring, prefix, FALSE, FALSE, TILDE_IGNORE);
    FFAddNewLine(ffstring);
    FFAddTextToString(ffstring, "includes: ", pbp->includes, "\n", FALSE, TRUE, TILDE_IGNORE);
    prefix = ";";
  }

  if (pbp->placement != NULL) {
    FFAddOneString (ffstring, prefix, FALSE, FALSE, TILDE_IGNORE);
    FFAddNewLine(ffstring);
    FFAddTextToString(ffstring, "placement: ", pbp->placement, "\n", FALSE, TRUE, TILDE_IGNORE);
    prefix = ";";
  }

  if (pbp->superfamily != NULL) {
    FFAddOneString (ffstring, prefix, FALSE, FALSE, TILDE_IGNORE);
    FFAddNewLine(ffstring);
    FFAddTextToString(ffstring, "superfamily: ", pbp->superfamily, "\n", FALSE, TRUE, TILDE_IGNORE);
    prefix = ";";
  }

  if (pbp->cross_reference != NULL) {
    FFAddOneString (ffstring, prefix, FALSE, FALSE, TILDE_IGNORE);
    FFAddNewLine(ffstring);
    FFAddTextToString(ffstring, "xref: ", pbp->cross_reference, "\n", FALSE, TRUE, TILDE_IGNORE);
    prefix = ";";
  }

  if (pbp->date != NULL) {
    FFAddOneString (ffstring, prefix, FALSE, FALSE, TILDE_IGNORE);
    FFAddNewLine(ffstring);
    FFAddTextToString (ffstring, "PIR dates: ", pbp->date, "\n", FALSE, TRUE, TILDE_IGNORE);
    prefix = ";";
  }

  if (pbp->had_punct) {
    FFAddOneString (ffstring, prefix, FALSE, FALSE, TILDE_IGNORE);
    FFAddNewLine(ffstring);
    FFAddOneString (ffstring, "punctuation in sequence", FALSE, FALSE, TILDE_IGNORE);
    prefix = ";";
  }

  if (pbp->seqref) {
    FFAddOneString (ffstring, prefix, FALSE, FALSE, TILDE_IGNORE);
    FFAddNewLine(ffstring);
    FFAddOneString (ffstring, "xrefs: ", FALSE, FALSE, TILDE_IGNORE);
    first = TRUE;
    for (sid = pbp->seqref; sid != NULL; sid = sid->next) {
      acc = NULL;
      if (first == FALSE) {
        FFAddOneString (ffstring, ", ", FALSE, FALSE, TILDE_IGNORE);
      }
      first = FALSE;
      SeqIdWrite (sid, id, PRINTID_TEXTID_ACC_VER, sizeof (id) - 1);
      acc = id;
      if (acc != NULL) {
        switch (sid->choice) {
          case SEQID_GENBANK:
            FFAddOneString (ffstring, "genbank ", FALSE, FALSE, TILDE_IGNORE);
            break; 
          case SEQID_EMBL:
            FFAddOneString (ffstring, "embl ", FALSE, FALSE, TILDE_IGNORE);
            break; 
          case SEQID_PIR:
            FFAddOneString (ffstring, "pir ", FALSE, FALSE, TILDE_IGNORE);
            break; 
          case SEQID_SWISSPROT:
            FFAddOneString (ffstring, "swissprot ", FALSE, FALSE, TILDE_IGNORE);
            break; 
          case SEQID_DDBJ:
            FFAddOneString (ffstring, "ddbj ", FALSE, FALSE, TILDE_IGNORE);
            break; 
          case SEQID_PRF:
            FFAddOneString (ffstring, "prf ", FALSE, FALSE, TILDE_IGNORE);
            break; 
          case SEQID_GI:
            FFAddOneString (ffstring, "gi: ", FALSE, FALSE, TILDE_IGNORE);
            break; 
          default:
            acc = NULL;
            break; 
        }
      }
      if (acc != NULL) {
        FFAddOneString (ffstring, acc, FALSE, FALSE, TILDE_IGNORE);
      }
    }
  }
  FFAddOneString (ffstring, ".", FALSE, FALSE, TILDE_IGNORE);
}

static void AddPRFBlock (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  BioseqPtr bsp
)

{
  SeqMgrDescContext  dcontext;
  PrfExtSrcPtr       extra;
  CharPtr            prefix = NULL;
  SeqDescrPtr        sdp;
  PrfBlockPtr        prf;

  if (bsp == NULL) return;
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_prf, &dcontext);
  if (sdp == NULL) return;
  prf = (PrfBlockPtr) sdp->data.ptrvalue;
  if (prf == NULL) return;
  if ( ffstring == NULL ) return;

  extra = prf->extra_src;
  if (extra != NULL) {

    if (extra->host != NULL) {
      FFAddTextToString(ffstring, "host:", extra->host, NULL, FALSE, TRUE, TILDE_IGNORE);
      prefix = ";";
    }

    if (extra->part != NULL) {
      FFAddOneString(ffstring, prefix, FALSE, FALSE, TILDE_IGNORE);
      FFAddNewLine(ffstring);
      FFAddTextToString(ffstring, "part: ", extra->part, NULL, FALSE, TRUE, TILDE_IGNORE);
      prefix = ";";
    }
    if (extra->state != NULL) {
      FFAddOneString(ffstring, prefix, FALSE, FALSE, TILDE_IGNORE);
      FFAddNewLine(ffstring);
      FFAddTextToString(ffstring, "state: ", extra->state, NULL, FALSE, TRUE, TILDE_IGNORE);
      prefix = ";";
    }
    if (extra->strain != NULL) {
      FFAddOneString(ffstring, prefix, FALSE, FALSE, TILDE_IGNORE);
      FFAddNewLine(ffstring);
      FFAddTextToString(ffstring, "strain: ", extra->strain, NULL, FALSE, TRUE, TILDE_IGNORE);
      prefix = ";";
    }
    if (extra->taxon != NULL) {
      FFAddOneString(ffstring, prefix, FALSE, FALSE, TILDE_IGNORE);
      FFAddNewLine(ffstring);
      FFAddTextToString(ffstring, "taxonomy: ", extra->taxon, NULL, FALSE, TRUE, TILDE_IGNORE);
      prefix = ";";
    }

    FFAddOneChar(ffstring, '.', FALSE);
  }
}

static void AddPDBBlock (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  BioseqPtr bsp
)

{
  SeqMgrDescContext  dcontext;
  CharPtr            dt;
  CharPtr            prefix = NULL;
  SeqDescrPtr        sdp;
  PdbBlockPtr        pdb;
  PdbRepPtr          replace;
  CharPtr            str;
  ValNodePtr         vnp;

  if (bsp == NULL) return;
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_pdb, &dcontext);
  if (sdp == NULL) return;
  pdb = (PdbBlockPtr) sdp->data.ptrvalue;
  if (pdb == NULL) return;

  if (pdb->deposition != NULL) {
    dt = asn2gb_PrintDate (pdb->deposition);
    FFAddTextToString (ffstring, "deposition: ", dt, NULL, FALSE, TRUE, TILDE_IGNORE);
    MemFree (dt);
    prefix = ";";
  }
  if (pdb->pdbclass != NULL) {
    FFAddOneString(ffstring, prefix, FALSE, FALSE, TILDE_IGNORE);
    FFAddNewLine(ffstring);
    FFAddTextToString(ffstring, "class: ", pdb->pdbclass, NULL, FALSE, TRUE, TILDE_IGNORE);
    prefix = ";";
  }
  if (pdb->source != NULL) {
    FFAddOneString(ffstring, prefix, FALSE, FALSE, TILDE_IGNORE);
    FFAddNewLine(ffstring);
    FFAddOneString(ffstring, "source: ", FALSE, TRUE, TILDE_IGNORE);
    prefix = NULL;
    for (vnp = pdb->source; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      FFAddTextToString (ffstring, prefix, str, NULL, FALSE, TRUE, TILDE_IGNORE);
      prefix = ", ";
    }
    prefix = ";";
  }
  if (pdb->exp_method != NULL) {
    FFAddOneString(ffstring, prefix, FALSE, FALSE, TILDE_IGNORE);
    FFAddNewLine(ffstring);
    FFAddTextToString(ffstring, "Exp. method: ", pdb->exp_method, NULL, FALSE, TRUE, TILDE_IGNORE);
    prefix = ";";
  }
  replace = pdb->replace;
  if (replace != NULL) {
    if (replace->ids != NULL) {
      FFAddOneString(ffstring, prefix, FALSE, FALSE, TILDE_IGNORE);
      FFAddNewLine(ffstring);
      FFAddOneString(ffstring, "ids replaced: ", FALSE, TRUE, TILDE_IGNORE);

      prefix = NULL;
      for (vnp = replace->ids; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        FFAddTextToString (ffstring, prefix, str, NULL, FALSE, TRUE, TILDE_IGNORE);
        prefix = ", ";
      }
      prefix = ";";
    }
    if (replace->date != NULL) {
      FFAddOneString(ffstring, prefix, FALSE, FALSE, TILDE_IGNORE);
      FFAddNewLine(ffstring);

      dt = asn2gb_PrintDate (replace->date);
      FFAddTextToString(ffstring, "replacement date: ", dt, NULL, FALSE, TRUE, TILDE_IGNORE);
      MemFree (dt);
      prefix = ";";
    }
  }

  FFAddOneChar(ffstring, '.', FALSE);
}

static void AddDbsourceBlock (
  Asn2gbWorkPtr awp
)

{
  IntAsn2gbJobPtr  ajp;
  Asn2gbSectPtr    asp;
  BaseBlockPtr     bbp;
  BioseqPtr        bsp;
  Char             buf [256];
  SeqFeatPtr       cds;
  DbtagPtr         db;
  GBSeqPtr         gbseq;
  SeqIdPtr         id;
  ValNodePtr       list = NULL;
  BioseqPtr        nuc;
  SeqIdPtr         sip;
  SeqLocPtr        slp;
  CharPtr          str;
  Boolean          unknown = TRUE;
  ValNodePtr       vnp;
  StringItemPtr    ffstring;

  if (awp == NULL) return;
  ajp = awp->ajp;
  if (ajp == NULL) return;
  asp = awp->asp;
  if (asp == NULL) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;

  bbp = Asn2gbAddBlock (awp, DBSOURCE_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  ffstring = FFGetString(ajp);
  if ( ffstring == NULL ) return;

  FFStartPrint (ffstring, awp->format, 0, 12, "DBSOURCE", 12, 5, 5, NULL, TRUE);

  sip = SeqIdSelect (bsp->id, dbsource_fasta_order, NUM_SEQID);

  if (sip != NULL) {

    switch (sip->choice) {
      case SEQID_PIR :
      case SEQID_SWISSPROT :
      case SEQID_PRF :
      case SEQID_PDB :
        if (WriteDbsourceID (sip, buf)) {
          FF_www_dbsource (ajp, ffstring, buf, TRUE, sip->choice);
          FFAddNewLine(ffstring);
          unknown = FALSE;
        }
        break;
      case SEQID_GENERAL :
        db = sip->data.ptrvalue;
        if (db == NULL) {
          break;
        }
        if (StringNCmp (db->db, "PIDe", 4) != 0 &&
            StringNCmp (db->db, "PIDd", 4) != 0 &&
            StringNCmp (db->db, "PID", 3) != 0) {
          break;
        }
        /* if (ChoicePID) found, continue on to next set of cases */
      case SEQID_EMBL :
      case SEQID_GENBANK :
      case SEQID_DDBJ :
      case SEQID_GIBBSQ :
      case SEQID_GIBBMT :
      case SEQID_OTHER :
      case SEQID_TPG :
      case SEQID_TPE :
      case SEQID_TPD :
      case SEQID_GI :
      case SEQID_GIIM :
        cds = SeqMgrGetCDSgivenProduct (bsp, NULL);
        if (cds != NULL) {
          nuc = BioseqFindFromSeqLoc (cds->location);
          if (nuc != NULL) {
            slp = SeqLocFindNext (cds->location, NULL);
            while (slp != NULL) {
              sip = SeqLocId (slp);
              AddToUniqueSipList (&list, sip);
              slp = SeqLocFindNext (cds->location, slp);
            }
            for (vnp = list; vnp != NULL; vnp = vnp->next) {
              id = (SeqIdPtr) vnp->data.ptrvalue;
              nuc = BioseqFindCore (id);
              sip = NULL;
              if (nuc != NULL) {
                sip = SeqIdSelect (nuc->id, dbsource_fasta_order, NUM_SEQID);
              } else if (id != NULL && id->choice == SEQID_GI) {
                sip = GetSeqIdForGI (id->data.intvalue);
              }
              if (sip == NULL) {
                sip = id;
              }
              if (sip != NULL) {
                if (WriteDbsourceID (sip, buf)) {
                  FF_www_dbsource (ajp, ffstring, buf, TRUE, sip->choice);
                  FFAddNewLine(ffstring);
                  unknown = FALSE;
                }
              }
            }
            ValNodeFree (list);
          }
        } else {
          if (WriteDbsourceID (sip, buf)) {
            FF_www_dbsource (ajp, ffstring, buf, TRUE, sip->choice);
            FFAddNewLine(ffstring);
            unknown = FALSE;
          }
        }
        break;
      default :
        break;
    }

    switch (sip->choice) {
      case SEQID_PIR :
        AddPIRBlock (ajp, ffstring, bsp);
        break;
      case SEQID_SWISSPROT :
        AddSPBlock (ajp, ffstring, bsp);
        break;
      case SEQID_PRF :
        AddPRFBlock (ajp, ffstring, bsp);
        break;
      case SEQID_PDB :
        AddPDBBlock (ajp, ffstring, bsp);
        break;
      default :
        break;
    }
  }

  if (unknown) {
    FFAddOneString (ffstring, "UNKNOWN", FALSE, FALSE, TILDE_TO_SPACES);
  }

  str = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, NULL);

  /* optionally populate gbseq for XML-ized GenBank format */

  if (ajp->gbseq) {
    gbseq = &asp->gbseq;
  } else {
    gbseq = NULL;
  }

  if (gbseq != NULL) {
    if (StringNCmp (str, "DBSOURCE    ", 12) == 0) {
      gbseq->source_db = StringSave (str + 12);
    } else {
      gbseq->source_db = StringSave (str);
    }
    CleanQualValue (gbseq->source_db);
    Asn2gnbkCompressSpaces (gbseq->source_db);
  }

  bbp->string = str;
  FFRecycleString(ajp, ffstring);
}

static void AddDateBlock (
  Asn2gbWorkPtr awp
)

{
  IntAsn2gbJobPtr    ajp;
  BaseBlockPtr       bbp;
  BioseqPtr          bsp;
  Char               date [40];
  SeqMgrDescContext  dcontext;
  DatePtr            dp;
  SeqDescrPtr        sdp;
  StringItemPtr      ffstring;

  if (awp == NULL) return;
  ajp = awp->ajp;
  if (ajp == NULL) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;

  ffstring = FFGetString(ajp);
  if ( ffstring == NULL ) return;

  bbp = Asn2gbAddBlock (awp, DATE_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  date [0] = '\0';

  dp = NULL;
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_create_date, &dcontext);
  if (sdp != NULL) {
    dp = (DatePtr) sdp->data.ptrvalue;
  }
  if (dp != NULL) {
    DateToGB (date, dp, FALSE);
  }
  if (StringHasNoText (date)) {
    StringCpy (date, "01-JAN-1900");
  }

  FFStartPrint (ffstring, awp->format, 0, 0, NULL, 0, 5, 5, "DT", TRUE);
  FFAddOneString (ffstring, date, FALSE, FALSE, TILDE_IGNORE);

  bbp->string = FFEndPrint(ajp, ffstring, awp->format, 0, 0, 5, 5, "DT");
  FFRecycleString(ajp, ffstring);

  bbp = Asn2gbAddBlock (awp, DATE_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  ffstring = FFGetString(ajp);

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_update_date, &dcontext);
  if (sdp != NULL) {
    dp = (DatePtr) sdp->data.ptrvalue;
  }
  if (dp != NULL) {
    DateToGB (date, dp, FALSE);
  }

  FFStartPrint (ffstring, awp->format, 0, 0, NULL, 0, 5, 5, "DT", FALSE);
  FFAddOneString (ffstring, date, FALSE, FALSE, TILDE_IGNORE);

  bbp->string = FFEndPrint(ajp, ffstring, awp->format, 0, 0, 5, 5, "DT");
  FFRecycleString(ajp, ffstring);
}


#define TOTAL_ESTKW 11
#define TOTAL_STSKW 5
#define TOTAL_GSSKW 2

static CharPtr EST_kw_array[ TOTAL_ESTKW] = {
  "EST", "EST PROTO((expressed sequence tag)", "expressed sequence tag",
  "EST (expressed sequence tag)", "EST(expressed sequence tag)",
  "partial cDNA sequence", "transcribed sequence fragment", "TSR",
  "putatively transcribed partial sequence", "UK putts"
};

static CharPtr GSS_kw_array [TOTAL_GSSKW] = {
  "GSS", "trapped exon"
};
static CharPtr STS_kw_array[TOTAL_STSKW] = {
  "STS", "STS(sequence tagged site)", "STS (sequence tagged site)",
  "STS sequence", "sequence tagged site"
};

static Int2 MatchArrayString (
  CharPtr array_string [],
  Int2 totalstr,
  CharPtr text
)

{
  Int2 i;

  for (i = 0; i < totalstr && text != NULL; i++) {
    if (StringCmp (array_string [i], text) == 0) {
      return (i);
    }
  }

  return (-1);
}

static Boolean CheckSpecialKeyword (
  Boolean is_est,
  Boolean is_sts,
  Boolean is_gss,
  CharPtr kwd
)

{
  if (kwd == NULL) return FALSE;

  if (is_est) {
    if (MatchArrayString (STS_kw_array, TOTAL_STSKW, kwd) != -1) return FALSE;
    if (MatchArrayString (GSS_kw_array, TOTAL_GSSKW, kwd) != -1) return FALSE;
  }

  if (is_sts) {
    if (MatchArrayString (EST_kw_array, TOTAL_ESTKW, kwd) != -1) return FALSE;
    if (MatchArrayString (GSS_kw_array, TOTAL_GSSKW, kwd) != -1) return FALSE;
  }

  if (is_gss) {
    if (MatchArrayString (STS_kw_array, TOTAL_STSKW, kwd) != -1) return FALSE;
    if (MatchArrayString (EST_kw_array, TOTAL_ESTKW, kwd) != -1) return FALSE;
  }

  return TRUE;
}

static Boolean KeywordAlreadyInList (
  ValNodePtr head,
  CharPtr kwd
)

{
  ValNodePtr  vnp;

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    if (StringICmp ((CharPtr) vnp->data.ptrvalue, kwd) == 0) return TRUE;
  }

  return FALSE;
}


static void AddKeywordsBlock (
  Asn2gbWorkPtr awp
)

{
  IntAsn2gbJobPtr    ajp;
  Asn2gbSectPtr      asp;
  BaseBlockPtr       bbp;
  BioseqPtr          bsp;
  SeqMgrDescContext  dcontext;
  EMBLBlockPtr       ebp;
  GBBlockPtr         gbp;
  GBSeqPtr           gbseq;
  ValNodePtr         head = NULL;
  IndxPtr            index;
  Boolean            is_est = FALSE;
  Boolean            is_gss = FALSE;
  Boolean            is_sts = FALSE;
  ValNodePtr         keywords;
  CharPtr            kwd;
  MolInfoPtr         mip;
  PirBlockPtr        pir;
  PrfBlockPtr        prf;
  SeqDescrPtr        sdp;
  SeqIdPtr           sip;
  SPBlockPtr         sp;
  CharPtr            str;
  ValNodePtr         vnp;
  StringItemPtr      ffstring;

  if (awp == NULL) return;
  ajp = awp->ajp;
  if (ajp == NULL) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;
  asp = awp->asp;
  if (asp == NULL) return;

  bbp = (BaseBlockPtr) Asn2gbAddBlock (awp, KEYWORDS_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  ffstring = FFGetString(ajp);
  if ( ffstring == NULL ) return;

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_TPG || sip->choice == SEQID_TPE || sip->choice == SEQID_TPD) {
      ValNodeCopyStr (&head, 0, "Third Party Annotation");
      ValNodeCopyStr (&head, 0, "; ");
      ValNodeCopyStr (&head, 0, "TPA");
    }
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
  if (sdp != NULL) {
    bbp->entityID = dcontext.entityID;
    bbp->itemID = dcontext.itemID;
    bbp->itemtype = OBJ_SEQDESC;

    mip = (MolInfoPtr) sdp->data.ptrvalue;
    if (mip != NULL) {
      switch (mip->tech) {
        case MI_TECH_htgs_1 :
          if (head != NULL) {
            ValNodeCopyStr (&head, 0, "; ");
          }
          ValNodeCopyStr (&head, 0, "HTG");
          ValNodeCopyStr (&head, 0, "; ");
          ValNodeCopyStr (&head, 0, "HTGS_PHASE1");
          break;
        case MI_TECH_htgs_2 :
          if (head != NULL) {
            ValNodeCopyStr (&head, 0, "; ");
          }
          ValNodeCopyStr (&head, 0, "HTG");
          ValNodeCopyStr (&head, 0, "; ");
          ValNodeCopyStr (&head, 0, "HTGS_PHASE2");
          break;
        case MI_TECH_htgs_3 :
          if (head != NULL) {
            ValNodeCopyStr (&head, 0, "; ");
          }
          ValNodeCopyStr (&head, 0, "HTG");
          break;
        case MI_TECH_est :
          if (head != NULL) {
            ValNodeCopyStr (&head, 0, "; ");
          }
          is_est = TRUE;
          ValNodeCopyStr (&head, 0, "EST");
          break;
        case MI_TECH_sts :
          if (head != NULL) {
            ValNodeCopyStr (&head, 0, "; ");
          }
          is_sts = TRUE;
          ValNodeCopyStr (&head, 0, "STS");
          break;
        case MI_TECH_survey :
          if (head != NULL) {
            ValNodeCopyStr (&head, 0, "; ");
          }
          is_gss = TRUE;
          ValNodeCopyStr (&head, 0, "GSS");
          break;
        case MI_TECH_fli_cdna :
          if (head != NULL) {
            ValNodeCopyStr (&head, 0, "; ");
          }
          ValNodeCopyStr (&head, 0, "FLI_CDNA");
          break;
        case MI_TECH_htgs_0 :
          if (head != NULL) {
            ValNodeCopyStr (&head, 0, "; ");
          }
          ValNodeCopyStr (&head, 0, "HTG");
          ValNodeCopyStr (&head, 0, "; ");
          ValNodeCopyStr (&head, 0, "HTGS_PHASE0");
          break;
        case MI_TECH_htc :
          if (head != NULL) {
            ValNodeCopyStr (&head, 0, "; ");
          }
          ValNodeCopyStr (&head, 0, "HTC");
          break;
        case MI_TECH_wgs :
          if (head != NULL) {
            ValNodeCopyStr (&head, 0, "; ");
          }
          ValNodeCopyStr (&head, 0, "WGS");
          break;
        default :
          break;
      }
    }
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, 0, &dcontext);
  while (sdp != NULL) {

    keywords = NULL;

    switch (dcontext.seqdesctype) {
      case Seq_descr_genbank :
        gbp = (GBBlockPtr) sdp->data.ptrvalue;
        if (gbp != NULL) {
          keywords = gbp->keywords;
        }
        break;
      case Seq_descr_embl :
        ebp = (EMBLBlockPtr) sdp->data.ptrvalue;
        if (ebp != NULL) {
          keywords = ebp->keywords;
        }
        break;
      case Seq_descr_pir :
        pir = (PirBlockPtr) sdp->data.ptrvalue;
        if (pir != NULL) {
          keywords = pir->keywords;
        }
        break;
      case Seq_descr_prf :
        prf = (PrfBlockPtr) sdp->data.ptrvalue;
        if (prf != NULL) {
          keywords = prf->keywords;
        }
        break;
      case Seq_descr_sp :
        sp = (SPBlockPtr) sdp->data.ptrvalue;
        if (sp != NULL) {
          keywords = sp->keywords;
        }
        break;
      default :
        break;
    }

    if (keywords != NULL) {
      bbp->entityID = dcontext.entityID;
      bbp->itemID = dcontext.itemID;
      bbp->itemtype = OBJ_SEQDESC;
    }

    for (vnp = keywords; vnp != NULL; vnp = vnp->next) {
      kwd = (CharPtr) vnp->data.ptrvalue;
      if (CheckSpecialKeyword (is_est, is_sts, is_gss, kwd)) {
        if (! KeywordAlreadyInList (head, kwd)) {
          if (head != NULL) {
            ValNodeCopyStr (&head, 0, "; ");
          }
          ValNodeCopyStr (&head, 0, kwd);
        }
      }
    }

    sdp = SeqMgrGetNextDescriptor (bsp, sdp, 0, &dcontext);
  }

  FFStartPrint( ffstring, awp->format, 0, 12, "KEYWORDS", 12, 5, 5, "KW", TRUE);
  str = MergeValNodeStrings (head);
  
  /* if no keywords were found, period will still be added by this call */
  if ( str != NULL ) {
    FFAddOneString (ffstring, str, TRUE, FALSE, TILDE_TO_SPACES);
  } else {
    FFAddOneChar(ffstring, '.', FALSE);
  }

  MemFree (str);

  /* optionally populate indexes for NCBI internal database */

  if (ajp->index) {
    index = &asp->index;
  } else {
    index = NULL;
  }

  if (index != NULL) {
    for (vnp = head; vnp != NULL; vnp = vnp->next) {
      kwd = (CharPtr) vnp->data.ptrvalue;
      if (StringCmp (kwd, "; ") == 0) continue;
      ValNodeCopyStrToHead (&(index->keywords), 0, kwd);
    }
  }

  /* optionally populate gbseq for XML-ized GenBank format */

  if (ajp->gbseq) {
    gbseq = &asp->gbseq;
  } else {
    gbseq = NULL;
  }

  if (gbseq != NULL) {
    for (vnp = head; vnp != NULL; vnp = vnp->next) {
      kwd = (CharPtr) vnp->data.ptrvalue;
      if (StringCmp (kwd, "; ") == 0) continue;
      ValNodeCopyStr (&(gbseq->keywords), 0, kwd);
    }
  }

  ValNodeFreeData (head);

  bbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "KW");

  FFRecycleString(ajp, ffstring);
}

static void AddSegmentBlock (
  Asn2gbWorkPtr awp,
  Boolean onePartOfSeg
)

{
  Char             acc [41];
  IntAsn2gbJobPtr  ajp;
  Asn2gbSectPtr    asp;
  BaseBlockPtr     bbp;
  Char             buf [32];
  GBSeqPtr         gbseq;
  StringItemPtr    ffstring;

  if (awp == NULL) return;
  ajp = awp->ajp;
  if (ajp == NULL) return;
  asp = awp->asp;
  if (asp == NULL) return;

  if (awp->seg < 1 || awp->numsegs < 1) return;

  bbp = Asn2gbAddBlock (awp, SEGMENT_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  ffstring = FFGetString(ajp);
  if ( ffstring == NULL ) return;


  FFStartPrint (ffstring, awp->format, 0, 12, "SEGMENT", 12, 5, 5, "XX", FALSE);

  if ( GetWWW(ajp) && awp->parent != NULL && onePartOfSeg) {
    sprintf (buf, "%d of ", (int) awp->seg);
    FFAddOneString (ffstring, buf, FALSE, FALSE, TILDE_IGNORE);
    SeqIdWrite (awp->parent->id, acc, PRINTID_TEXTID_ACC_VER, sizeof (acc) - 1);

    FFAddTextToString(ffstring, "<a href=", link_seq, NULL, FALSE, FALSE, TILDE_IGNORE);
    FFAddTextToString(ffstring, "val=", acc, ">", FALSE, FALSE, TILDE_IGNORE);

    sprintf (buf, "%ld", (long) awp->numsegs);
    FFAddOneString (ffstring, buf, FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
  } else {
    sprintf (buf, "%d of %ld", (int) awp->seg, (long) awp->numsegs);
    FFAddOneString (ffstring, buf, FALSE, TRUE, TILDE_TO_SPACES);
  }

  /* optionally populate gbseq for XML-ized GenBank format */

  if (ajp->gbseq) {
    gbseq = &asp->gbseq;
  } else {
    gbseq = NULL;
  }

  if (gbseq != NULL) {
    sprintf (buf, "%d of %ld", (int) awp->seg, (long) awp->numsegs);
    gbseq->segment = StringSave (buf);
  }

  bbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "XX");
  FFRecycleString(ajp, ffstring);
}

static void AddSourceBlock (
  Asn2gbWorkPtr awp
)

{
  IntAsn2gbJobPtr    ajp;
  BaseBlockPtr       bbp;
  BioseqPtr          bsp;
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext;
  GBBlockPtr         gbp;
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp;

  if (awp == NULL) return;
  ajp = awp->ajp;
  if (ajp == NULL) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;

  bbp = Asn2gbAddBlock (awp, SOURCE_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_genbank, &dcontext);
  if (sdp != NULL && (! ajp->newSourceOrg)) {
    gbp = (GBBlockPtr) sdp->data.ptrvalue;
    if (gbp != NULL && (! StringHasNoText (gbp->source))) {
      bbp->entityID = dcontext.entityID;
      bbp->itemID = dcontext.itemID;
      bbp->itemtype = OBJ_SEQDESC;
      return;
    }
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp != NULL) {
    bbp->entityID = dcontext.entityID;
    bbp->itemID = dcontext.itemID;
    bbp->itemtype = OBJ_SEQDESC;
  } else {
    sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_BIOSRC, 0, &fcontext);
    if (sfp != NULL) {
      bbp->entityID = fcontext.entityID;
      bbp->itemID = fcontext.itemID;
      bbp->itemtype = OBJ_SEQFEAT;
    }
  }
}

static void AddOrganismBlock (
  Asn2gbWorkPtr awp
)

{
  BaseBlockPtr       bbp;
  BioseqPtr          bsp;
  SeqFeatPtr         cds;
  SeqMgrDescContext  dcontext;
  BioseqPtr          dna;
  SeqMgrFeatContext  fcontext;
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp;

  if (awp == NULL) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;

  bbp = Asn2gbAddBlock (awp, ORGANISM_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp != NULL) {
    bbp->entityID = dcontext.entityID;
    bbp->itemID = dcontext.itemID;
    bbp->itemtype = OBJ_SEQDESC;
  } else {
    sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_BIOSRC, 0, &fcontext);
    if (sfp != NULL) {
      bbp->entityID = fcontext.entityID;
      bbp->itemID = fcontext.itemID;
      bbp->itemtype = OBJ_SEQFEAT;
    } else if (ISA_aa (bsp->mol)) {

      /* if protein with no sources, get sources applicable to DNA location of CDS */

      cds = SeqMgrGetCDSgivenProduct (bsp, &fcontext);
      if (cds != NULL) {
        sfp = SeqMgrGetOverlappingSource (cds->location, &fcontext);
        if (sfp != NULL) {
          bbp->entityID = fcontext.entityID;
          bbp->itemID = fcontext.itemID;
          bbp->itemtype = OBJ_SEQFEAT;
        } else {
          dna = BioseqFindFromSeqLoc (cds->location);
          if (dna != NULL) {
            sdp = SeqMgrGetNextDescriptor (dna, NULL, Seq_descr_source, &dcontext);
            if (sdp != NULL) {
              bbp->entityID = dcontext.entityID;
              bbp->itemID = dcontext.itemID;
              bbp->itemtype = OBJ_SEQDESC;
            }
          }
        }
      }
    }
  }
}

static RefBlockPtr AddPub (
  Asn2gbWorkPtr awp,
  ValNodePtr PNTR head,
  PubdescPtr pdp
)

{
  Char            buf [121];
  CitArtPtr       cap;
  CitBookPtr      cbp;
  CitGenPtr       cgp;
  CitJourPtr      cjp;
  CitPatPtr       cpp;
  CitSubPtr       csp;
  DatePtr         dp = NULL;
  Boolean         justuids = TRUE;
  ImprintPtr      imp = NULL;
  IntRefBlockPtr  irp;
  RefBlockPtr     rbp;
  ValNodePtr      vnp;

  if (awp == NULL || head == NULL || pdp == NULL) return NULL;

  rbp = (RefBlockPtr) MemNew (sizeof (IntRefBlock));
  if (rbp == NULL) return NULL;
  rbp->blocktype = REFERENCE_BLOCK;
  rbp->section = awp->currsection;

  rbp->serial = INT2_MAX;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    switch (vnp->choice) {
      case PUB_Gen :
        /* may be unpublished, or may be serial number of swiss-prot reference */
        cgp = (CitGenPtr) vnp->data.ptrvalue;
        if (cgp != NULL) {
          if (StringNICmp ("BackBone id_pub", cgp->cit, 15) != 0) {
            rbp->category = REF_CAT_UNP;
            dp = cgp->date;
            if (cgp->serial_number > 0) {
              rbp->serial = cgp->serial_number;
            }
            if (cgp->cit != NULL) {
              if (StringNICmp ("unpublished", cgp->cit, 11) != 0 &&
                  StringNICmp ("submitted", cgp->cit, 8) != 0 &&
                  StringNICmp ("to be published", cgp->cit, 15) != 0 &&
                  StringNICmp ("in press", cgp->cit, 8) != 0 &&
                  StringStr (cgp->cit, "Journal") == NULL) {
                if (cgp->serial_number == 0) {
                  MemFree (rbp);
                  return NULL;
                }
              }
            } else if (cgp->journal == NULL || cgp->date == NULL) {
              if (cgp->serial_number == 0) {
                MemFree (rbp);
                return NULL;
              }
            }
          }
        }
        break;
      case PUB_Sub :
        rbp->category = REF_CAT_SUB;
        csp = (CitSubPtr) vnp->data.ptrvalue;
        if (csp != NULL) {
          imp = csp->imp;
          if (imp != NULL) {
            dp = imp->date;
          }
          if (csp->date != NULL) {
            dp = csp->date;
          }
        }
        break;
      case PUB_Article:
        cap = (CitArtPtr) vnp->data.ptrvalue;
        if (cap != NULL) {
          switch (cap->from) {
            case 1:
              cjp = (CitJourPtr) cap->fromptr;
              if (cjp != NULL) {
                imp = (ImprintPtr) cjp->imp;
                if (imp != NULL) {
                  dp = imp->date;
                }
              }
              break;
            case 2:
              cbp = (CitBookPtr) cap->fromptr;
              if (cbp != NULL) {
                imp = (ImprintPtr) cbp->imp;
                if (imp != NULL) {
                  dp = imp->date;
                }
              }
              break;
            case 3:
              cbp = (CitBookPtr) cap->fromptr;
              if (cbp != NULL) {
                imp = (ImprintPtr) cbp->imp;
                if (imp != NULL) {
                  dp = imp->date;
                }
              }
              break;
            default:
              break;
          }
        }
        break;
      case PUB_Book:
        cbp = (CitBookPtr) vnp->data.ptrvalue;
        if (cbp != NULL) {
          imp = (ImprintPtr) cbp->imp;
          if (imp != NULL) {
            dp = imp->date;
          }
        }
        break;
      case PUB_Proc:
        cbp = (CitBookPtr) vnp->data.ptrvalue;
        if (cbp != NULL) {
          imp = (ImprintPtr) cbp->imp;
          if (imp != NULL) {
            dp = imp->date;
          }
        }
        break;
      case PUB_Patent :
        rbp->category = REF_CAT_PUB;
        cpp = (CitPatPtr) vnp->data.ptrvalue;
        if (cpp != NULL) {
          if (cpp->date_issue != NULL) {
            dp = (DatePtr) cpp->date_issue;
          } else if (cpp->app_date != NULL) {
            dp = (DatePtr) cpp->app_date;
          }
        }
        break;
      case PUB_Man:
        cbp = (CitBookPtr) vnp->data.ptrvalue;
        if (cbp != NULL) {
          imp = (ImprintPtr) cbp->imp;
          if (imp != NULL) {
            dp = imp->date;
          }
        }
        break;
      case PUB_Muid :
        rbp->muid = vnp->data.intvalue;
        rbp->category = REF_CAT_PUB;
        break;
      case PUB_PMid :
        rbp->pmid = vnp->data.intvalue;
        rbp->category = REF_CAT_PUB;
        break;
      default :
        break;
    }
    if (vnp->choice != PUB_Muid && vnp->choice != PUB_PMid) {
      justuids = FALSE;
    }
  }

  /* check for submitted vs. in-press */

  if (imp != NULL) {
    rbp->category = REF_CAT_PUB;
    switch (imp->prepub) {
      case 1 :
        rbp->category = REF_CAT_UNP;
        break;
      case 2 :
        rbp->category = REF_CAT_PUB;
        break;
      default :
        break;
    }
  }

  /* check for sites reftype */

  if (pdp->reftype != 0) {
    rbp->sites = pdp->reftype;
  }

  if (rbp->muid == 0 && rbp->pmid == 0) {
    vnp = pdp->pub;

    /* skip over just serial number */

    if (vnp != NULL && vnp->choice == PUB_Gen && vnp->next != NULL) {
      cgp = (CitGenPtr) vnp->data.ptrvalue;
      if (cgp != NULL) {
        if (StringNICmp ("BackBone id_pub", cgp->cit, 15) != 0) {
          if (cgp->cit == NULL && cgp->journal == NULL && cgp->date == NULL && cgp->serial_number) {
            vnp = vnp->next;
          }
        }
      }
    }

    if (PubLabelUnique (vnp, buf, sizeof (buf) - 1, OM_LABEL_CONTENT, TRUE) > 0) {
      rbp->uniquestr = StringSaveNoNull (buf);
    }
  }

  irp = (IntRefBlockPtr) rbp;
  irp->date = DateDup (dp);
  irp->justuids = justuids;
  /* if (justuids) { */
    irp->fig = StringSaveNoNull (pdp->fig);
    irp->maploc = StringSaveNoNull (pdp->maploc);
    irp->poly_a = pdp->poly_a;
  /* } */

  /* if not rejected by now, link in */

  ValNodeAddPointer (head, 0, rbp);

  return rbp;
}

static int LIBCALLBACK SortReferences (
  VoidPtr ptr1,
  VoidPtr ptr2,
  Boolean serialFirst
)

{
  int             compare;
  IntRefBlockPtr  irp1;
  IntRefBlockPtr  irp2;
  RefBlockPtr     rbp1;
  RefBlockPtr     rbp2;
  Int2            status;
  ValNodePtr      vnp1;
  ValNodePtr      vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  rbp1 = (RefBlockPtr) vnp1->data.ptrvalue;
  rbp2 = (RefBlockPtr) vnp2->data.ptrvalue;
  if (rbp1 == NULL || rbp2 == NULL) return 0;

  if (serialFirst) {
    if (rbp1->serial > rbp2->serial) {
      return 1;
    } else if (rbp1->serial < rbp2->serial) {
      return -1;
    }
  }

  /* usual first sort by published, unpublished, and cit-subs */

  if (rbp1->category > rbp2->category) {
    return 1;
  } else if (rbp1->category < rbp2->category) {
    return -1;
  }

  /* within class, sort by date, older publications first */

  irp1 = (IntRefBlockPtr) rbp1;
  irp2 = (IntRefBlockPtr) rbp2;
  status = DateMatch (irp1->date, irp2->date, FALSE);
  if (status == 1 || status == -1) return status;

  /* if dates (e.g., years) match, try to distinguish by uids */

  if (rbp1->pmid != 0 && rbp2->pmid != 0) {
    if (rbp1->pmid > rbp2->pmid) {
      return 1;
    } else if (rbp1->pmid < rbp2->pmid) {
      return -1;
    }
  }
  if (rbp1->muid != 0 && rbp2->muid != 0) {
    if (rbp1->muid > rbp2->muid) {
      return 1;
    } else if (rbp1->muid < rbp2->muid) {
      return -1;
    }
  }

  /* if same uid, one with just uids goes last to be excised but remembered */

  if ((rbp1->pmid != 0 && rbp2->pmid != 0) || (rbp1->muid != 0 && rbp2->muid != 0)) {
    if (irp1->justuids && (! irp2->justuids)) {
      return 1;
    } else if ((! irp1->justuids) && irp2->justuids) {
      return -1;
    }
  }

  /* put sites after pubs that refer to all or a range of bases */

  if (rbp1->sites > 0) {
    return 1;
  } else if (rbp2->sites > 0) {
    return -1;
  }

  /* for publication features, sort in explore index order */

  if (irp1->index > irp2->index) {
    return 1;
  } else if (irp1->index < irp2->index) {
    return -1;
  }

  /* next use author string */

  if (irp1->authstr != NULL && irp2->authstr != NULL) {
    compare = StringICmp (irp1->authstr, irp2->authstr);
    if (compare > 0) {
      return 1;
    } else if (compare < 0) {
      return -1;
    }
  }

  /* use unique label string to determine sort order */

  if (rbp1->uniquestr != NULL && rbp2->uniquestr != NULL) {
    compare = StringICmp (rbp1->uniquestr, rbp2->uniquestr);
    if (compare > 0) {
      return 1;
    } else if (compare < 0) {
      return -1;
    }
  }

  /* last resort for equivalent publication descriptors, sort in itemID order */

  if (rbp1->itemtype == OBJ_SEQDESC && rbp2->itemtype == OBJ_SEQDESC) {
    if (rbp1->itemID > rbp2->itemID) {
      return 1;
    } else if (rbp1->itemID < rbp2->itemID) {
      return -1;
    }
  }

  if (! serialFirst) {
    if (rbp1->serial > rbp2->serial) {
      return 1;
    } else if (rbp1->serial < rbp2->serial) {
      return -1;
    }
  }

  return 0;
}

static int LIBCALLBACK SortReferencesA (
  VoidPtr ptr1,
  VoidPtr ptr2
)

{
  return SortReferences (ptr1, ptr2, FALSE);
}

static int LIBCALLBACK SortReferencesB (
  VoidPtr ptr1,
  VoidPtr ptr2
)

{
  return SortReferences (ptr1, ptr2, TRUE);
}

static CharPtr GetAuthorsPlusConsortium (
  FmtType format,
  AuthListPtr alp
)

{
  CharPtr  consortium;
  CharPtr  str;
  CharPtr  tmp;

  consortium = NULL;
  str = GetAuthorsString (format, alp, &consortium, NULL, NULL);
  if (str == NULL) return consortium;
  if (consortium == NULL) return str;
  tmp = MemNew (StringLen (str) + StringLen (consortium) + 5);
  if (tmp == NULL) return NULL;
  StringCpy (tmp, str);
  StringCat (tmp, "; ");
  StringCat (tmp, consortium);
  MemFree (str);
  MemFree (consortium);
  return tmp;
}

static void GetRefsOnBioseq (
  Asn2gbWorkPtr awp,
  BioseqPtr target,
  BioseqPtr bsp,
  Int4 from,
  Int4 to
)

{
  IntAsn2gbJobPtr    ajp;
  AuthListPtr        alp;
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext;
  Int2               i;
  Int2               idx;
  IntRefBlockPtr     irp;
  Int4Ptr            ivals;
  Int4               left;
  SeqLocPtr          newloc;
  Int2               numivals;
  Boolean            okay;
  PubdescPtr         pdp;
  RefBlockPtr        rbp;
  Int4               right;
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp;
  SeqInt             sint;
  SeqIdPtr           sip;
  Boolean            split;
  Int4               start;
  Int4               stop;
  Uint1              strand;
  Boolean            takeIt;
  ValNode            vn;
  ValNodePtr         vnp;

  if (awp == NULL || target == NULL || bsp == NULL) return;
  ajp = awp->ajp;
  if (ajp == NULL) return;

  /* full length loc for descriptors */

  sint.from = 0;
  if (ajp->ajp.slp != NULL) {
    from = SeqLocStart (ajp->ajp.slp); /* other features use awp->slp for from and to */
  }
  if (ajp->ajp.slp != NULL) {
    sint.to = SeqLocLen (ajp->ajp.slp) - 1;
    to = SeqLocStop (ajp->ajp.slp); /* other features use awp->slp for from and to */
  } else {
    sint.to = bsp->length - 1;
  }
  sint.strand = Seq_strand_plus;
  sint.id = SeqIdStripLocus (SeqIdDup (SeqIdFindBest (bsp->id, 0)));
  sint.if_from = NULL;
  sint.if_to = NULL;

  vn.choice = SEQLOC_INT;
  vn.data.ptrvalue = (Pointer) &sint;
  vn.next = NULL;

  sdp = SeqMgrGetNextDescriptor (target, NULL, Seq_descr_pub, &dcontext);
  while (sdp != NULL) {

    /* check if descriptor on part already added on segmented bioseq */

    okay = TRUE;
    for (vnp = awp->pubhead; vnp != NULL && okay; vnp = vnp->next) {
      rbp = (RefBlockPtr) vnp->data.ptrvalue;
      if (rbp != NULL) {
        if (rbp->entityID == dcontext.entityID &&
            rbp->itemID == dcontext.itemID &&
            rbp->itemtype == OBJ_SEQDESC) {
          okay = FALSE;
        }
      }
    }

    if (okay) {
      pdp = (PubdescPtr) sdp->data.ptrvalue;
      rbp = AddPub (awp, &(awp->pubhead), pdp);
      if (rbp != NULL) {

        rbp->entityID = dcontext.entityID;
        rbp->itemID = dcontext.itemID;
        rbp->itemtype = OBJ_SEQDESC;

        irp = (IntRefBlockPtr) rbp;
        irp->loc = SeqLocMerge (target, &vn, NULL, FALSE, TRUE, FALSE);
        alp = GetAuthListPtr (pdp, NULL);
        if (alp != NULL) {
          irp->authstr = GetAuthorsPlusConsortium (awp->format, alp);
        }
        irp->index = 0;
      }
    }
    sdp = SeqMgrGetNextDescriptor (target, sdp, Seq_descr_pub, &dcontext);
  }

  SeqIdFree (sint.id);

  /* features are indexed on parent if segmented */

  bsp = awp->parent;

  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_PUB, 0, &fcontext);
  while (sfp != NULL) {
    ivals = fcontext.ivals;
    numivals = fcontext.numivals;
    if (ivals != NULL && numivals > 0) {

      /*
      idx = (numivals - 1) * 2;
      start = ivals [idx];
      stop = ivals [idx + 1];
      */

      takeIt = FALSE;
      for (i = 0, idx = 0; i < numivals; i++, idx += 2) {
        start = ivals [idx];
        stop = ivals [idx + 1];
        if ((start <= from && stop > from) ||
            (start < to && stop >= to) ||
            (start >= from && stop <= to)) {
          takeIt = TRUE;
        }
      }

      if (takeIt /* stop >= from && stop <= to */) {

        /*
        start = ivals [0] + 1;
        stop = ivals [idx + 1] + 1;
        */
        pdp = (PubdescPtr) sfp->data.value.ptrvalue;
        rbp = AddPub (awp, &(awp->pubhead), pdp);
        if (rbp != NULL) {

          rbp->entityID = fcontext.entityID;
          rbp->itemID = fcontext.itemID;
          rbp->itemtype = OBJ_SEQFEAT;

          irp = (IntRefBlockPtr) rbp;
          irp->loc = SeqLocMerge (target, sfp->location, NULL, FALSE, TRUE, FALSE);
          if (ajp->ajp.slp != NULL) {
            sip = SeqIdParse ("lcl|dummy");
            left = GetOffsetInBioseq (ajp->ajp.slp, bsp, SEQLOC_LEFT_END);
            right = GetOffsetInBioseq (ajp->ajp.slp, bsp, SEQLOC_RIGHT_END);
            strand = SeqLocStrand (ajp->ajp.slp);
            split = FALSE;
            newloc = SeqLocReMap (sip, ajp->ajp.slp, irp->loc, 0, FALSE);
            /*
            newloc = SeqLocCopyRegion (sip, irp->loc, bsp, left, right, strand, &split);
            */
            SeqIdFree (sip);
            if (newloc != NULL) {
              A2GBSeqLocReplaceID (newloc, ajp->ajp.slp);
              irp->loc = SeqLocFree (irp->loc);
              irp->loc = newloc;
            }
          }
          alp = GetAuthListPtr (pdp, NULL);
          if (alp != NULL) {
            irp->authstr = GetAuthorsPlusConsortium (awp->format, alp);
          }
          irp->index = fcontext.index;
        }
      }
    }

    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_PUB, 0, &fcontext);
  }
}

static Boolean LIBCALLBACK GetRefsOnSeg (
  SeqLocPtr slp,
  SeqMgrSegmentContextPtr context
)

{
  Asn2gbWorkPtr  awp;
  BioseqPtr      bsp;
  Int4           from;
  SeqLocPtr      loc;
  SeqEntryPtr    oldscope;
  SeqEntryPtr    sep;
  SeqIdPtr       sip;
  Int4           to;

  if (slp == NULL || context == NULL) return FALSE;
  awp = (Asn2gbWorkPtr) context->userdata;

  from = context->cumOffset;
  to = from + context->to - context->from;

  sip = SeqLocId (slp);
  if (sip == NULL) {
    loc = SeqLocFindNext (slp, NULL);
    if (loc != NULL) {
      sip = SeqLocId (loc);
    }
  }
  if (sip == NULL) return TRUE;

  /* reference descriptors only on parts within entity */

  sep = GetTopSeqEntryForEntityID (awp->entityID);
  oldscope = SeqEntrySetScope (sep);
  bsp = BioseqFind (sip);
  SeqEntrySetScope (oldscope);

  if (bsp != NULL) {
    GetRefsOnBioseq (awp, awp->refs, bsp, from, to);
    return TRUE;
  }

  /* if we ever want to fetch remote references, code goes here */

  return TRUE;
}

static Boolean AddReferenceBlock (
  Asn2gbWorkPtr awp
)

{
  IntAsn2gbJobPtr    ajp;
  AuthListPtr        alp;
  Asn2gbSectPtr      asp;
  BioseqPtr          bsp;
  SeqFeatPtr         cds;
  Boolean            combine;
  SeqMgrFeatContext  context;
  CitSubPtr          csp;
  BioseqPtr          dna;
  Boolean            excise;
  ValNodePtr         head = NULL;
  Int2               i;
  IntRefBlockPtr     irp;
  Boolean            is_embl = FALSE;
  Boolean            is_patent = FALSE;
  IntRefBlockPtr     lastirp;
  RefBlockPtr        lastrbp;
  ValNodePtr         next;
  Int2               numReferences;
  ValNodePtr         PNTR prev;
  RefBlockPtr        rbp;
  RefBlockPtr        PNTR referenceArray;
  BioseqPtr          refs;
  SubmitBlockPtr     sbp;
  SeqIdPtr           sip;
  SeqLocPtr          slp;
  BioseqPtr          target;
  ValNodePtr         vnp;

  if (awp == NULL) return FALSE;
  ajp = awp->ajp;
  if (ajp == NULL) return FALSE;
  asp = awp->asp;
  if (asp == NULL) return FALSE;
  bsp = awp->bsp;
  refs = awp->refs;
  if (bsp == NULL || refs == NULL) return FALSE;

  /* collect publications on bioseq */

  awp->pubhead = NULL;
  GetRefsOnBioseq (awp, bsp, refs, awp->from, awp->to);
  target = bsp;

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_EMBL) {
      is_embl = TRUE;
    } else if (sip->choice == SEQID_PATENT) {
      is_patent = TRUE;
    }
  }

  if (bsp->repr == Seq_repr_seg) {

    /* collect publication descriptors on local parts */

    SeqMgrExploreSegments (bsp, (Pointer) awp, GetRefsOnSeg);
    target = awp->refs;
  }

  if (awp->pubhead == NULL && ISA_aa (bsp->mol)) {

    /* if protein with no pubs, get pubs applicable to DNA location of CDS */

    cds = SeqMgrGetCDSgivenProduct (bsp, &context);
    if (cds != NULL) {
      dna = BioseqFindFromSeqLoc (cds->location);
      if (dna != NULL) {
        GetRefsOnBioseq (awp, dna, dna, context.left, context.right);
        target = dna;
      }
    }
  }

  head = awp->pubhead;
  awp->pubhead = NULL;

  if (head == NULL && awp->ssp == NULL) return FALSE;

  /* sort by pub/unpub/sites/sub, then date, finally existing serial */

  head = SortValNode (head, SortReferencesA);

  if (awp->ssp != NULL) {

    /* add seq-submit citation */

    rbp = (RefBlockPtr) MemNew (sizeof (IntRefBlock));
    if (rbp != NULL) {
      irp = (IntRefBlockPtr) rbp;

      rbp->blocktype = REFERENCE_BLOCK;
      rbp->section = awp->currsection;
      rbp->serial = INT2_MAX;
      rbp->category = REF_CAT_SUB;

      rbp->entityID = ajp->ajp.entityID;
      rbp->itemID = 1;
      rbp->itemtype = OBJ_SEQSUB_CIT;

      sbp = awp->ssp->sub;
      if (sbp != NULL) {
        csp = sbp->cit;
        if (csp != NULL) {
          alp = GetAuthListPtr (NULL, csp);
          if (alp != NULL) {
            irp->authstr = GetAuthorsPlusConsortium (awp->format, alp);
          }
        }
      }

      if (awp->citSubsFirst) {

        /* for DDBJ, add seq-submit citation to beginning of list */

        vnp = ValNodeNew (NULL);
        if (vnp != NULL) {
          vnp->choice = 0;
          vnp->data.ptrvalue = (VoidPtr) rbp;
          vnp->next = head;
          head = vnp;
        }

      } else {

        /* for GENBANK and EMBL add seq-submit citation to end of list */

        ValNodeAddPointer (&head, 0, rbp);
      }
    }
  }

  /* unique references, excise duplicates from list */

  prev = &(head);
  vnp = head;
  lastrbp = NULL;
  while (vnp != NULL) {
    excise = FALSE;
    combine = TRUE;
    next = vnp->next;
    rbp = (RefBlockPtr) vnp->data.ptrvalue;
    if (lastrbp != NULL) {
      lastirp = (IntRefBlockPtr) lastrbp;
      if (rbp != NULL) {
        irp = (IntRefBlockPtr) rbp;
        if (lastrbp->pmid != 0 && rbp->pmid != 0) {
          if (lastrbp->pmid == rbp->pmid) {
            excise = TRUE;
          }
        } else if (lastrbp->muid != 0 && rbp->muid != 0) {
          if (lastrbp->muid == rbp->muid) {
            excise = TRUE;
          }
        } else if (lastrbp->uniquestr != NULL && rbp->uniquestr != NULL) {
          if (StringICmp (lastrbp->uniquestr, rbp->uniquestr) == 0) {
            if (SeqLocCompare (irp->loc, lastirp->loc) == SLC_A_EQ_B) {
              if (StringICmp (irp->authstr, lastirp->authstr) == 0) {

                /* L76496.1 - removing duplicate submission pubs */
                excise = TRUE;
              }
            }
          }
        }
        if (excise && lastrbp->sites == 0 && rbp->sites > 0) {
          /* real range trumps sites */
          combine = FALSE;
        }
      }
    }
    if (rbp != NULL) {
      irp = (IntRefBlockPtr) rbp;
      if (irp->justuids) {
        /* do not allow justuids reference to appear by itself - S79174.1 */
        excise = TRUE;
        /* justuids should still combine, even if no authors - S67070.1 */
      } else if (is_embl && is_patent) {
        /* EMBL patent records do not need author or title - A29528.1 */
      } else if (StringHasNoText (irp->authstr)) {
        /* do not allow no author reference to appear by itself - U07000.1 */
        excise = TRUE;
        combine = FALSE;
      }
    }
    if (excise) {
      *prev = vnp->next;
      vnp->next = NULL;

      /* combine locations of duplicate references */

      irp = (IntRefBlockPtr) rbp;
      lastirp = (IntRefBlockPtr) lastrbp;
      if (combine) {
        if (lastirp != NULL) {
          slp = SeqLocMerge (target, lastirp->loc, irp->loc, FALSE, TRUE, FALSE);
          lastirp->loc = SeqLocFree (lastirp->loc);
          lastirp->loc = slp;
        }
        if (irp != NULL && lastirp != NULL) {
          if ((rbp->muid == lastrbp->muid && rbp->muid != 0) ||
              (rbp->pmid == lastrbp->pmid && rbp->pmid != 0)) {
            if (lastirp->fig == NULL) {
              lastirp->fig = StringSaveNoNull (irp->fig);
            }
            if (lastirp->maploc == NULL) {
              lastirp->maploc = StringSaveNoNull (irp->maploc);
            }
            lastirp->poly_a = irp->poly_a;
          }
        }
      }

      /* and remove duplicate reference */

      MemFree (rbp->uniquestr);
      DateFree (irp->date);
      SeqLocFree (irp->loc);
      MemFree (irp->authstr);
      MemFree (irp->fig);
      MemFree (irp->maploc);
      MemFree (rbp);
      ValNodeFree (vnp);

    } else {

      prev = &(vnp->next);
      lastrbp = rbp;
    }
    vnp = next;
  }

  /* resort by existing serial, then pub/unpub/sites/sub, then date */

  head = SortValNode (head, SortReferencesB);

  if (head == NULL) return FALSE;

  /* assign serial numbers */

  for (vnp = head, i = 1; vnp != NULL; vnp = vnp->next, i++) {
    rbp = (RefBlockPtr) vnp->data.ptrvalue;
    if (rbp != NULL) {
      rbp->serial = i;
    }
  }

  /* allocate reference array for this section */

  numReferences = i - 1;
  asp->numReferences = numReferences;

  if (numReferences > 0) {
    referenceArray = (RefBlockPtr PNTR) MemNew (sizeof (RefBlockPtr) * (numReferences + 1));
    asp->referenceArray = referenceArray;

    if (referenceArray != NULL) {

      /* fill in reference array */

      for (vnp = head, i = 0; vnp != NULL && i < numReferences; vnp = vnp->next, i++) {
        referenceArray [i] = (RefBlockPtr) vnp->data.ptrvalue;
      }
    }
  }

  /* finally link into blocks for current section */

  ValNodeLink (&(awp->lastblock), head);
  vnp = awp->lastblock;
  if (vnp == NULL) return FALSE;
  while (vnp->next != NULL) {
    vnp = vnp->next;
  }

  awp->lastblock = vnp;
  if (awp->blockList == NULL) {
    awp->blockList = vnp;
  }

  return TRUE;
}

static void AddHistCommentString (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  CharPtr prefix,
  CharPtr suffix,
  DatePtr dp,
  SeqIdPtr ids
)

{
  Int2      count = 0;
  Char      buf [256];
  Boolean   first;
  Int4      gi = 0;
  SeqIdPtr  sip;
  CharPtr   strd;
  
  if (dp == NULL || ids == NULL || prefix == NULL || suffix == NULL || ffstring == NULL) return;

  strd = asn2gb_PrintDate (dp);
  if (strd == NULL) {
    strd = StringSave ("?");
  }

  for (sip = ids; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_GI) {
      gi = (long) sip->data.intvalue;
      count++;
    }
  }

  if (count > 1) {
    sprintf (buf, "%s or before %s %s", prefix, strd, suffix);
  } else {
    sprintf (buf, "%s %s %s", prefix, strd, suffix);
  }
  FFAddOneString (ffstring, buf, FALSE, FALSE, TILDE_EXPAND);

  MemFree (strd);

  if (gi == 0) {
    FFAddOneString (ffstring, " gi:?", FALSE, FALSE, TILDE_EXPAND);
    return;
  }

  first = TRUE;
  for (sip = ids; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_GI) {
      gi = (long) sip->data.intvalue;
      if (! first) {
        FFAddOneString (ffstring, ",", FALSE, FALSE, TILDE_IGNORE);
      }
      first = FALSE;
      if ( GetWWW(ajp) ) {
        FFAddOneString (ffstring, " gi:", FALSE, FALSE, TILDE_IGNORE);
        FFAddTextToString (ffstring, "<a href=", link_seq, NULL, FALSE, FALSE, TILDE_IGNORE);
        sprintf (buf, "%ld", (long) gi);
        FFAddTextToString (ffstring, "val=", buf, ">", FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, buf, FALSE, FALSE, TILDE_EXPAND);
        FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
      } else {
        sprintf (buf, " gi:%ld", (long) gi);
        FFAddOneString (ffstring, buf, FALSE, FALSE, TILDE_EXPAND);
      }
    }
  }

  FFAddOneString (ffstring, ".", FALSE, FALSE, TILDE_EXPAND);
}

static void AddHTGSCommentString (
  StringItemPtr ffstring,
  BioseqPtr bsp,
  MolInfoPtr mip
)

{
  CharPtr      buf = NULL;
  Char         buffer [256];
  Int4         buflen = 0;
  DeltaSeqPtr  dsp;
  ValNodePtr   head = NULL;
  Int4         num_s = 0;
  Int4         num_g = 0;
  CharPtr      str = NULL;

  if (bsp == NULL || mip == NULL || mip->tech < 2) return;

  if (bsp->repr == Seq_repr_delta) {
    for (dsp = (DeltaSeqPtr) bsp->seq_ext, buflen = 0; dsp != NULL; dsp = dsp->next) {
      buflen += 80;
    }
    if (buflen > 0) {
      buf = MemNew ((size_t) (buflen + 1));
      if (buf == NULL) return;
      CountGapsInDeltaSeq (bsp, &num_s, &num_g, NULL, NULL, buf, buflen);
    }
  }

  if (mip->tech == MI_TECH_htgs_0) {

    if (num_s > 0) {
      sprintf (buffer, "* NOTE: This record contains %ld individual~", (long) (num_g + 1));
      ValNodeCopyStr (&head, 0, buffer);
      ValNodeCopyStr (&head, 0, "* sequencing reads that have not been assembled into~");
      ValNodeCopyStr (&head, 0, "* contigs. Runs of N are used to separate the reads~");
      ValNodeCopyStr (&head, 0, "* and the order in which they appear is completely~");
      ValNodeCopyStr (&head, 0, "* arbitrary. Low-pass sequence sampling is useful for~");
      ValNodeCopyStr (&head, 0, "* identifying clones that may be gene-rich and allows~");
      ValNodeCopyStr (&head, 0, "* overlap relationships among clones to be deduced.~");
      ValNodeCopyStr (&head, 0, "* However, it should not be assumed that this clone~");
      ValNodeCopyStr (&head, 0, "* will be sequenced to completion. In the event that~");
      ValNodeCopyStr (&head, 0, "* the record is updated, the accession number will~");
      ValNodeCopyStr (&head, 0, "* be preserved.");
    }
    ValNodeCopyStr (&head, 0, "~");
    ValNodeCopyStr (&head, 0, buf);

  } else if (mip->tech == MI_TECH_htgs_1) {

    ValNodeCopyStr (&head, 0, "* NOTE: This is a \"working draft\" sequence.");
    if (num_s > 0) {
      sprintf (buffer, " It currently~* consists of %ld contigs. The true order of the pieces~", (long) (num_g + 1));
      ValNodeCopyStr (&head, 0, buffer);
      ValNodeCopyStr (&head, 0, "* is not known and their order in this sequence record is~");
      ValNodeCopyStr (&head, 0, "* arbitrary. Gaps between the contigs are represented as~");
      ValNodeCopyStr (&head, 0, "* runs of N, but the exact sizes of the gaps are unknown.");
    }
    ValNodeCopyStr (&head, 0, "~* This record will be updated with the finished sequence~");
    ValNodeCopyStr (&head, 0, "* as soon as it is available and the accession number will~");
    ValNodeCopyStr (&head, 0, "* be preserved.");
    ValNodeCopyStr (&head, 0, "~");
    ValNodeCopyStr (&head, 0, buf);

  } else if (mip->tech == MI_TECH_htgs_2) {

    ValNodeCopyStr (&head, 0, "* NOTE: This is a \"working draft\" sequence.");
    if (num_s > 0) {
      sprintf (buffer, " It currently~* consists of %ld contigs. Gaps between the contigs~", (long) (num_g + 1));
      ValNodeCopyStr (&head, 0, buffer);
      ValNodeCopyStr (&head, 0, "* are represented as runs of N. The order of the pieces~");
      ValNodeCopyStr (&head, 0, "* is believed to be correct as given, however the sizes~");
      ValNodeCopyStr (&head, 0, "* of the gaps between them are based on estimates that have~");
      ValNodeCopyStr (&head, 0, "* provided by the submittor.");
    }
    ValNodeCopyStr (&head, 0, "~* This sequence will be replaced~");
    ValNodeCopyStr (&head, 0, "* by the finished sequence as soon as it is available and~");
    ValNodeCopyStr (&head, 0, "* the accession number will be preserved.");
    ValNodeCopyStr (&head, 0, "~");
    ValNodeCopyStr (&head, 0, buf);

  } else if ((str = StringForSeqTech (mip->tech)) != NULL) {

      sprintf (buffer, "Method: %s.", str);
      ValNodeCopyStr (&head, 0, buffer);
  }

  MemFree (buf);

  str = MergeValNodeStrings (head);

  FFAddOneString (ffstring, str, TRUE, TRUE, TILDE_EXPAND);

  MemFree (str);
  ValNodeFreeData (head);
}

static void AddWGSMasterCommentString (
  StringItemPtr ffstring,
  BioseqPtr bsp,
  CharPtr wgsaccn,
  CharPtr wgsname
)

{
  BioSourcePtr       biop;
  Char               buf [256];
  SeqMgrDescContext  dcontext;
  CharPtr            first = NULL;
  ValNodePtr         head = NULL;
  CharPtr            last = NULL;
  ObjectIdPtr        oip;
  OrgRefPtr          orp;
  SeqDescrPtr        sdp;
  CharPtr            taxname = NULL;
  UserFieldPtr       ufp;
  UserObjectPtr      uop;
  Char               ver [16];

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp != NULL) {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
    if (biop != NULL) {
      orp = biop->org;
      if (orp != NULL) {
        taxname = orp->taxname;
      }
    }
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &dcontext);
  while (sdp != NULL) {
    uop = (UserObjectPtr) sdp->data.ptrvalue;
    if (uop != NULL) {
      oip = uop->type;
      if (oip != NULL && StringICmp (oip->str, "WGSProjects") == 0) {
        for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
          oip = ufp->label;
          if (oip == NULL || oip->str == NULL || ufp->choice != 1) continue;
          if (StringICmp (oip->str, "WGS_accession_first") == 0) {
            first = (CharPtr) ufp->data.ptrvalue;
          } else if (StringICmp (oip->str, "WGS_accession_last") == 0) {
            last = (CharPtr) ufp->data.ptrvalue;
          }
        }
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &dcontext);
  }

  if (StringHasNoText (taxname)) {
    taxname = "?";
  }
  if (StringHasNoText (first)) {
    first = "?";
  }
  if (StringHasNoText (last)) {
    last = "?";
  }
  ver [0] = '\0';
  if (StringLen (wgsname) == 12) {
    StringCpy (ver, wgsname + 4);
    ver [2] = '\0';
  } else if (StringLen (wgsname) == 15) {
    StringCpy (ver, wgsname + 7);
    ver [2] = '\0';
  }

  sprintf (buf, "The %s whole genome shotgun (WGS) project has the project accession %s.", taxname, wgsaccn);
  FFAddOneString(ffstring, buf, TRUE, FALSE, TILDE_EXPAND);

  sprintf (buf, "  This version of the project (%s) has the accession number %s,", ver, wgsname);
  FFAddOneString(ffstring, buf, TRUE, FALSE, TILDE_EXPAND);

  if (StringCmp (first, last) != 0) {
    sprintf (buf, " and consists of sequences %s-%s.", first, last);
    FFAddOneString(ffstring, buf, TRUE, FALSE, TILDE_EXPAND);
  } else {
    sprintf (buf, " and consists of sequence %s.", first);
    FFAddOneString(ffstring, buf, TRUE, FALSE, TILDE_EXPAND);
  }
}


static CharPtr GetMolInfoCommentString (
  BioseqPtr bsp,
  MolInfoPtr mip
)

{
  Boolean  is_aa;
  CharPtr  str = NULL;

  if (bsp == NULL || mip == NULL) return NULL;

  is_aa = ISA_aa (bsp->mol);
  switch (mip->completeness) {
    case 1 :
      str = "COMPLETENESS: full length";
      break;
    case 2 :
      str = "COMPLETENESS: not full length";
      break;
    case 3 :
      if (is_aa) {
        str = "COMPLETENESS: incomplete on the amino end";
      } else {
        str = "COMPLETENESS: incomplete on the 5' end";
      }
      break;
    case 4 :
      if (is_aa) {
        str = "COMPLETENESS: incomplete on the carboxy end";
      } else {
        str = "COMPLETENESS: incomplete on the 3' end";
      }
      break;
    case 5 :
      str = "COMPLETENESS: incomplete on both ends";
      break;
    case 6 :
      if (is_aa) {
        str = "COMPLETENESS: complete on the amino end";
      } else {
        str = "COMPLETENESS: complete on the 5' end";
      }
      break;
    case 7 :
      if (is_aa) {
        str = "COMPLETENESS: complete on the carboxy end";
      } else {
        str = "COMPLETENESS: complete on the 3' end";
      }
      break;
    default :
      str = "COMPLETENESS: unknown";
      break;
  }

  return str;
}

static CharPtr GetStrForBankit (
  UserObjectPtr uop
)

{
  CharPtr       bic = NULL, uvc = NULL, ptr;
  ObjectIdPtr   oip;
  UserFieldPtr  ufp;

  if (uop == NULL) return NULL;
  if ((oip = uop->type) == NULL) return NULL;
  if (StringCmp (oip->str, "Submission") != 0) return NULL;

  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    oip = ufp->label;
    if (StringCmp(oip->str, "UniVecComment") == 0) {
      uvc = ufp->data.ptrvalue;
    } else if (StringCmp(oip->str, "AdditionalComment") == 0) {
      bic = ufp->data.ptrvalue;
    }
  }

  if (uvc == NULL && bic == NULL) return NULL;

  ptr = (CharPtr) MemNew (StringLen (uvc) + StringLen (bic) + 45);
  if (uvc != NULL && bic != NULL) {
    sprintf (ptr, "Vector Explanation: %s~Bankit Comment: %s", uvc, bic);
  } else if (uvc != NULL) {
    sprintf (ptr, "Vector Explanation: %s", uvc);
  } else if (bic != NULL) {
    sprintf (ptr, "Bankit Comment: %s", bic);
  }

  return ptr;
}

static CharPtr reftxt0 = "The reference sequence was derived from ";
static CharPtr reftxt1 = " This record is predicted by genome sequence analysis and is not yet supported by experimental evidence. ";
static CharPtr reftxt2 = " This record has not yet been subject to final NCBI review. ";
static CharPtr reftxt3 = " The mRNA record is supported by experimental evidence; however, the coding sequence is predicted. ";
static CharPtr reftxt4 = " This record has undergone preliminary review of the sequence, but has not yet been subject to final NCBI review. ";
static CharPtr reftxt5 = " This record has been curated by ";
static CharPtr reftxt6 = " This RefSeq record is provided to represent a collection of whole genome shotgun sequences. ";

static CharPtr GetStatusForRefTrack (
  UserObjectPtr uop
)

{
  CharPtr       st;
  ObjectIdPtr   oip;
  UserFieldPtr  ufp, urf = NULL;

  if (uop == NULL) return NULL;
  if ((oip = uop->type) == NULL) return NULL;
  if (StringCmp (oip->str, "RefGeneTracking") != 0) return NULL;
  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    oip = ufp->label;
    if (StringCmp(oip->str, "Assembly") == 0) {
      urf = ufp;
    }
  }
  if (urf == NULL || urf->choice != 11) return NULL;
  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    oip = ufp->label;
    if (StringCmp (oip->str, "Status") == 0) {
      st = (CharPtr) ufp->data.ptrvalue;
      if (StringCmp (st, "Inferred") == 0) {
        return "INFERRED ";
      } else if (StringCmp (st, "Provisional") == 0) {
        return "PROVISIONAL ";
      } else if (StringCmp (st, "Predicted") == 0) {
        return "PREDICTED ";
      } else if (StringCmp (st, "Validated") == 0) {
        return "VALIDATED ";
      } else if (StringCmp (st, "Reviewed") == 0) {
        return "REVIEWED ";
      } else if (StringCmp (st, "WGS") == 0) {
        return "WGS ";
      }
    }
  }
  return NULL;
}


static void AddStrForRefTrack (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  UserObjectPtr uop
)

{
  CharPtr       accn, curator = NULL, st;
  ObjectIdPtr   oip;
  UserFieldPtr  ufp, tmp, u, urf = NULL;
  Int2          i = 0;
  Int2          review = 0,len;
  Boolean       is_accn;

  if ( uop == NULL || ffstring == NULL ) return;
  if ((oip = uop->type) == NULL) return;
  if (StringCmp (oip->str, "RefGeneTracking") != 0) return;

  len = StringLen (reftxt0);
  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    oip = ufp->label;
    if (StringCmp(oip->str, "Assembly") == 0) {
      urf = ufp;
    }
    if (StringCmp (oip->str, "Status") == 0) {
      st = (CharPtr) ufp->data.ptrvalue;
      if (StringCmp (st, "Inferred") == 0) {
        review = 1;
      } else if (StringCmp (st, "Provisional") == 0) {
        review = 2;
      } else if (StringCmp (st, "Predicted") == 0) {
        review = 3;
      } else if (StringCmp (st, "Validated") == 0) {
        review = 4;
      } else if (StringCmp (st, "Reviewed") == 0) {
        review = 5;
      } else if (StringCmp (st, "WGS") == 0) {
        review = 6;
      }
    } else if (StringCmp (oip->str, "Collaborator") == 0) {
      st = (CharPtr) ufp->data.ptrvalue;
      if (! StringHasNoText (st)) {
        curator = st;
      }
    }
  }
  if (urf != NULL && urf->choice == 11) {
    for (tmp = urf->data.ptrvalue; tmp != NULL; tmp = tmp->next) {
      for (u = tmp->data.ptrvalue; u != NULL; u = u->next) {
        oip = u->label;
        if (StringCmp (oip->str, "accession") == 0 ||
            StringCmp (oip->str, "name") == 0) {
          i++;
        }
      }
    }
    if ( GetWWW(ajp) ) {
      FFAddTextToString(ffstring, "<a href=", ref_link, ">", FALSE, FALSE, TILDE_IGNORE);
      FFAddOneString (ffstring, "REFSEQ", FALSE, FALSE, TILDE_IGNORE);
      FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
    } else {
      FFAddOneString (ffstring, "REFSEQ", FALSE, FALSE, TILDE_IGNORE);
    }
    FFAddOneString (ffstring, ":", FALSE, FALSE, TILDE_IGNORE);
    if (review == 1) {
      FFAddOneString (ffstring, reftxt1, FALSE, FALSE, TILDE_IGNORE);
    } else if (review == 2) {
      FFAddOneString (ffstring, reftxt2, FALSE, FALSE, TILDE_IGNORE);
    } else if (review == 3) {
      FFAddOneString (ffstring, reftxt3, FALSE, FALSE, TILDE_IGNORE);
    } else if (review == 4) {
      FFAddOneString (ffstring, reftxt4, FALSE, FALSE, TILDE_IGNORE);
    } else if (review == 5) {
      if (curator == NULL) {
        curator = "NCBI staff";
      }
      FFAddOneString (ffstring, reftxt5, FALSE, FALSE, TILDE_IGNORE);
      FFAddOneString (ffstring, curator, FALSE, FALSE, TILDE_IGNORE);
      FFAddOneString (ffstring, ". ", FALSE, FALSE, TILDE_IGNORE);
    } else if (review == 6) {
      FFAddOneString (ffstring, reftxt6, FALSE, FALSE, TILDE_IGNORE);
    }
    if (i > 0) {
      FFAddOneString (ffstring, reftxt0, FALSE, FALSE, TILDE_IGNORE);

      for (tmp = urf->data.ptrvalue; tmp != NULL; tmp = tmp->next) {
        is_accn = TRUE;
        for (u = tmp->data.ptrvalue; u != NULL; u = u->next) {
          oip = u->label;
          if (StringCmp (oip->str, "accession") == 0) break;
          if (StringCmp (oip->str, "name") == 0) {
            is_accn = FALSE;
            break;
          }
        }
        if (u == NULL) continue;
        accn = (CharPtr) u->data.ptrvalue;
        if (StringHasNoText (accn)) continue;
        if (is_accn && GetWWW(ajp) ) {
          FFAddTextToString(ffstring, "<a href=", link_seq, NULL, FALSE, FALSE, TILDE_IGNORE);
          FFAddTextToString(ffstring, "val=", accn, ">", FALSE, FALSE, TILDE_IGNORE);
          FFAddOneString (ffstring, accn, FALSE, FALSE, TILDE_IGNORE);
          FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
        } else {
          FFAddOneString (ffstring, accn, FALSE, FALSE, TILDE_IGNORE);
        }
        if (tmp->next != NULL) {
          ufp = tmp->next;
          if (ufp->next != NULL) {
            FFAddOneString (ffstring, ", ", FALSE, FALSE, TILDE_IGNORE);
          } else {
            FFAddOneString (ffstring, " and ", FALSE, FALSE, TILDE_IGNORE);
          }
        }
      }
      FFAddOneString (ffstring, ".", FALSE, FALSE, TILDE_EXPAND);
    }
  }
}


static CharPtr reftxt11 = "This model reference sequence was predicted from NCBI contig";
static CharPtr reftxt12 = "by automated computational analysis";
static CharPtr reftxt13 = "using gene prediction method:";

static void FindModelEvidenceUop (
  UserObjectPtr uop,
  Pointer userdata
)

{
  ObjectIdPtr         oip;
  UserObjectPtr PNTR  uopp;

  if (uop == NULL || userdata == NULL) return;
  uopp = (UserObjectPtr PNTR) userdata;
  oip = uop->type;
  if (oip == NULL) return;
  if (StringCmp (oip->str, "ModelEvidence") == 0) {
    *uopp = uop;
  }
}

static Boolean DoGetAnnotationComment (
   BioseqPtr bsp,
   CharPtr PNTR namep,
   CharPtr PNTR methodp,
   BoolPtr mrnaEv,
   BoolPtr estEv
)

{
  SeqMgrDescContext  dcontext;
  CharPtr            method = NULL;
  UserObjectPtr      moduop;
  CharPtr            name = NULL;
  ObjectIdPtr        oip;
  SeqDescrPtr        sdp;
  UserFieldPtr       ufp;
  UserObjectPtr      uop;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &dcontext);
  while (sdp != NULL) {
    uop = (UserObjectPtr) sdp->data.ptrvalue;
    if (uop != NULL) {
      moduop = NULL;
      VisitUserObjectsInUop (uop, (Pointer) &moduop, FindModelEvidenceUop);
      if (moduop != NULL) {
        oip = moduop->type;
        if (oip != NULL && StringCmp(oip->str, "ModelEvidence") == 0) {
          for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
            oip = ufp->label;
            if (oip == NULL) continue;
            if (StringCmp(oip->str, "Contig Name") == 0) {
              name = (CharPtr) ufp->data.ptrvalue;
            } else if (StringCmp(oip->str, "Method") == 0) {
              method = (CharPtr) ufp->data.ptrvalue;
            } else if (StringCmp(oip->str, "mRNA") == 0) {
              *mrnaEv = TRUE;
            } else if (StringCmp(oip->str, "EST") == 0) {
              *estEv = TRUE;
            }
          }
        }
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &dcontext);
  }
  if (StringHasNoText (name)) return FALSE;
  *namep = name;
  if (! StringHasNoText (method)) {
    *methodp = method;
  }
  return TRUE;
}

static Boolean GetAnnotationComment (
   BioseqPtr bsp,
   CharPtr PNTR namep,
   CharPtr PNTR methodp,
   BoolPtr mrnaEv,
   BoolPtr estEv
)

{
  SeqFeatPtr  cds;

  if (DoGetAnnotationComment (bsp, namep, methodp, mrnaEv, estEv)) return TRUE;
  if (ISA_aa (bsp->mol)) {
    cds = SeqMgrGetCDSgivenProduct (bsp, NULL);
    if (cds != NULL) {
      bsp = BioseqFindFromSeqLoc (cds->location);
      if (bsp != NULL) {
        return DoGetAnnotationComment (bsp, namep, methodp, mrnaEv, estEv);
      }
    }
  }
  return FALSE;
}

static void FindGeneFeat (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  SeqFeatPtr PNTR  sfpp;

  if (sfp->data.choice != SEQFEAT_GENE) return;
  sfpp = (SeqFeatPtr PNTR) userdata;
  *sfpp = sfp;
}

static void FindLocusId (
  ValNodePtr dbxref,
  CharPtr locusIDp
)

{
  DbtagPtr     dbt;
  ObjectIdPtr  oip;
  ValNodePtr   vnp;

  for (vnp = dbxref; vnp != NULL; vnp = vnp->next) {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt == NULL) continue;
    if (StringICmp (dbt->db, "LocusID") != 0 && StringICmp (dbt->db, "InterimID") != 0) continue;
    oip = dbt->tag;
    if (oip == NULL) continue;
    if (oip->str != NULL) {
      StringCpy (locusIDp, oip->str);
    } else if (oip->id > 0) {
      sprintf (locusIDp, "%ld", (long) oip->id);
    }
  }
}

static Boolean GetGeneAndLocus (
  BioseqPtr bsp,
  CharPtr PNTR genep,
  CharPtr locusIDp,
  CharPtr taxIDp
)

{
  BioSourcePtr       biop;
  DbtagPtr           dbt;
  SeqMgrDescContext  dcontext;
  SeqFeatPtr         gene = NULL;
  GeneRefPtr         grp;
  ObjectIdPtr        oip;
  OrgRefPtr          orp;
  SeqDescrPtr        sdp;
  SeqEntryPtr        sep;
  CharPtr            str;
  ValNodePtr         syn;
  ValNodePtr         vnp;

  sep = GetTopSeqEntryForEntityID (bsp->idx.entityID);
  if (sep == NULL) return FALSE;
  VisitFeaturesInSep (sep, (Pointer) &gene, FindGeneFeat);
  if (gene == NULL) return FALSE;

  grp = (GeneRefPtr) gene->data.value.ptrvalue;
  if (grp == NULL) return FALSE;
  if (! StringHasNoText (grp->locus)) {
    *genep = grp->locus;
  } else {
    syn = grp->syn;
    if (syn != NULL) {
      str = (CharPtr) syn->data.ptrvalue;
      if (! StringHasNoText (str)) {
        *genep = str;
      }
    }
  }
  FindLocusId (gene->dbxref, locusIDp);
  FindLocusId (grp->db, locusIDp);

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp != NULL) {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
    if (biop != NULL) {
      orp = biop->org;
      if (orp != NULL) {
        for (vnp = orp->db; vnp != NULL; vnp = vnp->next) {
          dbt = (DbtagPtr) vnp->data.ptrvalue;
          if (dbt == NULL) continue;
          if (StringCmp (dbt->db, "taxon") == 0) {
            oip = dbt->tag;
            if (oip == NULL) continue;
            if (oip->str != NULL) {
              StringCpy (taxIDp, oip->str);
            } else if (oip->id > 0) {
              sprintf (taxIDp, "%ld", (long) oip->id);
            }
          }
        }
      }
    }
  }

  if (genep == NULL || StringHasNoText (locusIDp)) return FALSE;

  return TRUE;
}

static CharPtr reftxt21 = "GENOME ANNOTATION REFSEQ:  NCBI contigs are derived from assembled genomic sequence data. They may include both draft and finished sequence.";

static CharPtr tpaString = "THIRD PARTY ANNOTATION DATABASE: This TPA record uses data from DDBJ/EMBL/GenBank ";
static CharPtr defTpaString = "THIRD PARTY ANNOTATION DATABASE: This entry contains annotation of previously submitted data.";
static CharPtr nsAreGapsString = "The strings of n's in this record represent gaps between contigs, and the length of each string corresponds to the length of the gap.";

static Boolean IsTpa (
  BioseqPtr bsp,
  BoolPtr isRefSeqP
)

{
  DbtagPtr  dbt;
  Boolean   has_bankit = FALSE;
  Boolean   has_genbank = FALSE;
  Boolean   has_gi = FALSE;
  Boolean   has_local = FALSE;
  Boolean   has_refseq = FALSE;
  Boolean   has_smart = FALSE;
  Boolean   has_tpa = FALSE;
  SeqIdPtr  sip;

  if (bsp == NULL || bsp->id == NULL) return FALSE;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    switch (sip->choice) {
      case SEQID_LOCAL :
        has_local = TRUE;
        break;
      case SEQID_GENBANK :
      case SEQID_EMBL :
      case SEQID_DDBJ :
        has_genbank = TRUE;
        break;
      case SEQID_OTHER :
        has_refseq = TRUE;
        if (isRefSeqP != NULL) {
          *isRefSeqP = TRUE;
        }
        break;
      case SEQID_GI :
        has_gi = TRUE;
        break;
      case SEQID_TPG :
      case SEQID_TPE :
      case SEQID_TPD :
        has_tpa = TRUE;
        break;
      case SEQID_GENERAL :
        dbt = (DbtagPtr) sip->data.ptrvalue;
        if (dbt != NULL) {
          if (StringICmp (dbt->db, "BankIt") == 0) {
            has_bankit = TRUE;
          }
          if (StringICmp (dbt->db, "TMSMART") == 0) {
            has_smart = TRUE;
          }
        }
        break;
      default :
        break;
    }
  }

  if (has_genbank) return FALSE;
  if (has_tpa) return TRUE;
  if (has_refseq) return TRUE;
  if (has_bankit) return TRUE;
  if (has_smart) return TRUE;
  if (has_gi) return FALSE;
  if (has_local) return TRUE;

  return FALSE;
}

static CharPtr GetStrForTpaOrRefSeqHist (
  UserObjectPtr uop,
  BioseqPtr bsp
)

{
  Boolean      accn;
  Char         buf [64];
  DbtagPtr     dbt;
  Int4         gi;
  ValNodePtr   head = NULL;
  SeqHistPtr   hist;
  SeqIdPtr     id;
  Boolean      isRefSeq = FALSE;
  Boolean      minus1;
  Boolean      minus2;
  ObjectIdPtr  oip;
  SeqAlignPtr  salp;
  SeqAlignPtr  salptmp;
  SeqIdPtr     sip;
  Int4         start;
  Int4         stop;
  CharPtr      str;
  Char         tmp [80];

  if (uop == NULL) return NULL;
  if ((oip = uop->type) == NULL) return NULL;
  if (StringCmp (oip->str, "TpaAssembly") != 0 &&
      StringCmp (oip->str, "RefGeneTracking") != 0) return NULL;
  if (! IsTpa (bsp, &isRefSeq)) return NULL;

  hist = bsp->hist;
  if (hist != NULL && hist->assembly != NULL) {
    salp = SeqAlignListDup (hist->assembly);
    AlnMgr2IndexLite (salp);
    AlnMgr2SortAlnSetByNthRowPos (salp, 1);
    salptmp = (SeqAlignPtr) (salp->segs);
    while (salptmp != NULL) {
      AlnMgr2GetNthSeqRangeInSA (salptmp, 1, &start, &stop);
      sip = AlnMgr2GetNthSeqIdPtr (salptmp, 2);
      if (sip != NULL) {
        id = NULL;
        accn = FALSE;
        if (sip->choice == SEQID_GI) {
          gi = (Int4) sip->data.intvalue;
          if (GetAccnVerFromServer (gi, buf)) {
            accn = TRUE;
          } else {
            id = GetSeqIdForGI (gi);
          }
        } else {
          id = SeqIdDup (sip);
        }
        if (id != NULL || accn) {
          if (head == NULL) {
            if (isRefSeq) {
              ValNodeCopyStr (&head, 0, "REFSEQ_SPAN         PRIMARY_IDENTIFIER PRIMARY_SPAN        COMP");
            } else {
              ValNodeCopyStr (&head, 0, "TPA_SPAN            PRIMARY_IDENTIFIER PRIMARY_SPAN        COMP");
            }
          }
          if (id != NULL) {
            SeqIdWrite (id, buf, PRINTID_TEXTID_ACC_VER, sizeof (buf) - 1);
            if (id->choice == SEQID_GENERAL) {
              dbt = (DbtagPtr) id->data.ptrvalue;
              if (dbt != NULL && StringICmp (dbt->db, "ti") == 0) {
                StringCpy (buf, "TI");
                SeqIdWrite (id, buf + 2, PRINTID_TEXTID_ACC_VER, sizeof (buf) - 3);
              }
            }
          }
          sprintf (tmp, "~%ld-%ld                                        ",
                   (long) (start + 1), (long) (stop + 1));
          /*
          i = 39 - StringLen (buf);
          if (i > 0) {
            tmp [i] = '\0';
          } else {
            tmp [21] = '\0';
          }
          */
          tmp [21] = '\0';
          StringCat (buf, "                                        ");
          buf [18] = '\0';
          StringCat (tmp, buf);
          AlnMgr2GetNthSeqRangeInSA (salptmp, 2, &start, &stop);
          sprintf (buf, " %ld-%ld                                        ",
                   (long) (start + 1), (long) (stop + 1));
          buf [21] = '\0';
          StringCat (tmp, buf);
          minus1 = (Boolean) (AlnMgr2GetNthStrand (salptmp, 1) == Seq_strand_minus);
          minus2 = (Boolean) (AlnMgr2GetNthStrand (salptmp, 2) == Seq_strand_minus);
          if (minus1 || minus2) {
            if (! (minus1 && minus2)) {
              StringCat (tmp, "c");
            }
          }
          ValNodeCopyStr (&head, 0, tmp);
        }
        SeqIdFree (id);
      }
      SeqIdFree (sip);
      salptmp = salptmp->next;
    }
    SeqAlignFree (salp);
  }

  if (head == NULL) return NULL;

  str = MergeValNodeStrings (head);
  ValNodeFreeData (head);

  return str;
}

static CharPtr GetStrForTPA (
  UserObjectPtr uop,
  BioseqPtr bsp
)

{
  Char          ch;
  UserFieldPtr  curr;
  Int2          i;
  Char          id [41];
  Boolean       isRefSeq = FALSE;
  Int2          j;
  size_t        len;
  ObjectIdPtr   oip;
  CharPtr       ptr;
  CharPtr       str;
  CharPtr       tmp;
  UserFieldPtr  ufp;

  if (uop == NULL) return NULL;
  if ((oip = uop->type) == NULL) return NULL;
  if (StringCmp (oip->str, "TpaAssembly") != 0) return NULL;
  if (! IsTpa (bsp, &isRefSeq)) return NULL;
  if (isRefSeq) return NULL;

  len = StringLen (tpaString) + StringLen ("entries ") + StringLen ("and ") + 5;
  i = 0;
  for (curr = uop->data; curr != NULL; curr = curr->next) {
    if (curr->choice != 11) continue;
    for (ufp = curr->data.ptrvalue; ufp != NULL; ufp = ufp->next) {
      if (ufp->choice != 1) continue;
      oip = ufp->label;
      if (oip == NULL || StringICmp (oip->str, "accession") != 0) continue;
      str = (CharPtr) ufp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      len += StringLen (str) + 2;
      i++;
    }
  }
  if (i == 0) return NULL;

  ptr = (CharPtr) MemNew (len);
  if (ptr == NULL) return NULL;
  StringCpy (ptr, tpaString);
  if (i > 1) {
    StringCat (ptr, "entries ");
  } else {
    StringCat (ptr, "entry ");
  }

  j = 0;
  for (curr = uop->data; curr != NULL; curr = curr->next) {
    if (curr->choice != 11) continue;
    for (ufp = curr->data.ptrvalue; ufp != NULL; ufp = ufp->next) {
      if (ufp->choice != 1) continue;
      oip = ufp->label;
      if (oip == NULL || StringICmp (oip->str, "accession") != 0) continue;
      str = (CharPtr) ufp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      StringNCpy_0 (id, str, sizeof (id));
      tmp = id;
      ch = *tmp;
      while (ch != '\0') {
        if (IS_LOWER (ch)) {
          *tmp = TO_UPPER (ch);
        }
        tmp++;
        ch = *tmp;
      }
      if (j == i - 1 && i > 1) {
        StringCat (ptr, " and ");
      } else if (j > 0) {
        StringCat (ptr, ", ");
      }
      StringCat (ptr, id);
      j++;
    }
  }

  return ptr;
}

static CharPtr GetStrForGenome (
  UserObjectPtr uop,
  BioseqPtr bsp
)

{
  ObjectIdPtr  oip;

  if (uop == NULL) return NULL;
  if ((oip = uop->type) == NULL) return NULL;
  if (StringCmp (oip->str, "GenomeInfo") != 0) return NULL;

  /* !!! need to implement !!! */

  return NULL;
}

static void AddPrimaryBlock (
  Asn2gbWorkPtr awp
)

{
  IntAsn2gbJobPtr    ajp;
  Asn2gbSectPtr      asp;
  BaseBlockPtr       bbp = NULL;
  BioseqPtr          bsp;
  SeqMgrDescContext  dcontext;
  Boolean            didTpaHist = FALSE;
  GBSeqPtr           gbseq;
  SeqDescrPtr        sdp;
  CharPtr            str;
  UserObjectPtr      uop;
  StringItemPtr      ffstring;

  if (awp == NULL) return;
  ajp = awp->ajp;
  if (ajp == NULL) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;
  asp = awp->asp;
  if (asp == NULL) return;

  ffstring = FFGetString(ajp);
  if ( ffstring == NULL ) return;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &dcontext);
  while (sdp != NULL) {

    uop = (UserObjectPtr) sdp->data.ptrvalue;
    if (uop != NULL) {

      if (! didTpaHist) {
        str = GetStrForTpaOrRefSeqHist (uop, bsp);
        if (str != NULL) {

          bbp = (BaseBlockPtr) Asn2gbAddBlock (awp, PRIMARY_BLOCK, sizeof (BaseBlock));
          if (bbp != NULL) {

            bbp->entityID = dcontext.entityID;
            bbp->itemID = dcontext.itemID;
            bbp->itemtype = OBJ_SEQDESC;

            FFStartPrint (ffstring, awp->format, 0, 12, "PRIMARY", 12, 5, 5, "PR", TRUE);

            FFAddOneString (ffstring, str, FALSE, FALSE, TILDE_EXPAND);

            bbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "PR");

            /* optionally populate gbseq for XML-ized GenBank format */

            if (ajp->gbseq) {
              gbseq = &asp->gbseq;
            } else {
              gbseq = NULL;
            }

            if (gbseq != NULL) {
              gbseq->primary = StringSave (str);
            }
          }
          MemFree (str);
          didTpaHist = TRUE;
        }
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &dcontext);
  }

  FFRecycleString(ajp, ffstring);
}

static void AddCommentBlock (
  Asn2gbWorkPtr awp
)

{
  IntAsn2gbJobPtr    ajp;
  BioseqPtr          bsp;
  Char               buf [128];
  CommentBlockPtr    cbp = NULL;
  Char               ch;
  Boolean            didGenome = FALSE;
  Boolean            didRefTrack = FALSE;
  Boolean            didTPA = FALSE;
  DbtagPtr           dbt;
  SeqMgrDescContext  dcontext;
  DeltaSeqPtr        dsp;
  Boolean            estEv = FALSE;
  SeqMgrFeatContext  fcontext;
  Boolean            first = TRUE;
  GBBlockPtr         gbp;
  CharPtr            geneName = NULL;
  Int4               gi = 0;
  CommentBlockPtr    gsdbcbp = NULL;
  Int4               gsdbid = 0;
  Boolean            has_gaps = FALSE;
  SeqHistPtr         hist;
  Boolean            is_other = FALSE;
  Boolean            is_tpa = FALSE;
  Boolean            is_wgs = FALSE;
  SeqLitPtr          litp;
  Char               locusID [32];
  CharPtr            method = NULL;
  MolInfoPtr         mip;
  Boolean            mrnaEv = FALSE;
  CharPtr            name = NULL;
  Boolean            okay;
  BioseqPtr          parent;
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp;
  Boolean            showGBBSource = FALSE;
  SeqIdPtr           sip;
  CharPtr            str;
  Char               taxID [32];
  TextSeqIdPtr       tsip;
  UserObjectPtr      uop;
  CharPtr            wgsaccn = NULL;
  CharPtr            wgsname = NULL;
  StringItemPtr      ffstring = NULL;

  if (awp == NULL) return;
  ajp = awp->ajp;
  if (ajp == NULL) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;

  ffstring = FFGetString(ajp);
  if ( ffstring ==  NULL ) return;

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_OTHER) {
      tsip = (TextSeqIdPtr) sip->data.ptrvalue;

      if (tsip != NULL) {
        is_other = TRUE;
        if (StringNCmp(tsip->accession, "NT_", 3) == 0 || StringNCmp(tsip->accession, "NW_", 3) == 0) {

          cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
          if (cbp != NULL) {

            cbp->first = first;
            first = FALSE;

            if (cbp->first) {
              FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
            } else {
              FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
            }

            FFAddOneString (ffstring, reftxt21, TRUE, FALSE, TILDE_EXPAND);

            cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
            FFRecycleString(ajp, ffstring);
            ffstring = FFGetString(ajp);
          }

        } else if (StringNCmp(tsip->accession, "XP_", 3) == 0 ||
                   StringNCmp(tsip->accession, "XM_", 3) == 0 ||
                   StringNCmp(tsip->accession, "XR_", 3) == 0 ||
                   StringNCmp(tsip->accession, "ZP_", 3) == 0) {

          name = NULL;
          method = NULL;
          mrnaEv = FALSE;
          estEv = FALSE;
          if (GetAnnotationComment (bsp, &name, &method, &mrnaEv, &estEv)) {

            cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
            if (cbp != NULL) {

              cbp->first = first;
              first = FALSE;

              if (cbp->first) {
                FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
              } else {
                FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
              }

              FFAddOneString (ffstring, "GENOME ANNOTATION ", FALSE, FALSE, TILDE_IGNORE);

              if ( GetWWW(ajp) ) {
                FFAddTextToString (ffstring, "<a href=", ref_link, ">", FALSE, FALSE, TILDE_IGNORE);
              }
              FFAddOneString (ffstring, "REFSEQ", FALSE, FALSE, TILDE_IGNORE);
              if ( GetWWW(ajp) ) {
                FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
              }
              FFAddOneString (ffstring, ":  ", FALSE, FALSE, TILDE_IGNORE);

              FFAddTextToString (ffstring, NULL, reftxt11, " ", FALSE, FALSE, TILDE_IGNORE);

              if ( GetWWW(ajp) ) {
                FFAddTextToString (ffstring, "<a href=", nt_link, name, FALSE, FALSE, TILDE_IGNORE);
                FFAddOneString (ffstring, ">", FALSE, FALSE, TILDE_IGNORE);
              }
              FFAddOneString (ffstring, name, FALSE, FALSE, TILDE_IGNORE);
              if ( GetWWW(ajp) ) {
                FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
              }

              FFAddOneString (ffstring, " ", FALSE, FALSE, TILDE_IGNORE);
              FFAddOneString (ffstring, reftxt12, FALSE, FALSE, TILDE_IGNORE);

              if (method != NULL) {
                FFAddOneString (ffstring, " ", FALSE, FALSE, TILDE_IGNORE);
                FFAddOneString (ffstring, reftxt13, FALSE, FALSE, TILDE_IGNORE);
                FFAddOneString (ffstring, " ", FALSE, FALSE, TILDE_IGNORE);
                FFAddOneString (ffstring, method, FALSE, FALSE, TILDE_IGNORE);
              }

              if (mrnaEv || estEv) {
                FFAddOneString (ffstring, ", supported by ", FALSE, FALSE, TILDE_IGNORE);
                if (mrnaEv && estEv) {
                  FFAddOneString (ffstring, "mRNA and EST ", FALSE, FALSE, TILDE_IGNORE);
                } else if (mrnaEv) {
                  FFAddOneString (ffstring, "mRNA ", FALSE, FALSE, TILDE_IGNORE);
                } else {
                  FFAddOneString (ffstring, "EST ", FALSE, FALSE, TILDE_IGNORE);
                }
                geneName = NULL;
                locusID [0] = '\0';
                taxID [0] = '\0';
                if ( GetWWW(ajp) && GetGeneAndLocus (bsp, &geneName, locusID, taxID)) {
                  FFAddTextToString (ffstring, "<a href=", ev_link, NULL, FALSE, FALSE, TILDE_IGNORE);
                  FFAddTextToString (ffstring, "contig=", name, NULL, FALSE, FALSE, TILDE_IGNORE);
                  FFAddTextToString (ffstring, "&gene=", geneName, NULL, FALSE, FALSE, TILDE_IGNORE);
                  FFAddTextToString (ffstring, "&lid=", locusID, NULL, FALSE, FALSE, TILDE_IGNORE);
                  if (! StringHasNoText (taxID)) {
                    FFAddTextToString (ffstring, "&taxid=", taxID, NULL, FALSE, FALSE, TILDE_IGNORE);
                  }
                  FFAddOneString (ffstring, ">", FALSE, FALSE, TILDE_IGNORE);
                  FFAddOneString (ffstring, "evidence", FALSE, FALSE, TILDE_IGNORE);
                  FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
                } else {
                  FFAddOneString (ffstring, "evidence", FALSE, FALSE, TILDE_IGNORE);
                }
              }

              FFAddOneString (ffstring, ".", FALSE, FALSE, TILDE_IGNORE);

              FFAddOneString (ffstring, "~Also see:~    ", FALSE, FALSE, TILDE_EXPAND);

              if ( GetWWW(ajp) ) {
                FFAddTextToString (ffstring, "<a href=", doc_link, ">", FALSE, FALSE, TILDE_IGNORE);
              }
              FFAddOneString (ffstring, "Documentation", FALSE, FALSE, TILDE_IGNORE);
              if ( GetWWW(ajp) ) {
                FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
              }

              FFAddOneString (ffstring, " of NCBI's Annotation Process~    ", FALSE, FALSE, TILDE_EXPAND);

              cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
              FFRecycleString(ajp, ffstring);
              ffstring = FFGetString(ajp);
            }
          }
        } else {
          if (StringLen (tsip->accession) == 15) {
            is_wgs = TRUE;
            if (StringCmp (tsip->accession + 9, "000000") == 0) {
              wgsaccn = tsip->accession;
              wgsname = tsip->name; /* master accession has 8 zeroes, name has project version plus 6 zeroes */
            }
          }
        }
      }

    } else if (sip->choice == SEQID_TPG || sip->choice == SEQID_TPE || sip->choice == SEQID_TPD) {

      is_tpa = TRUE;

    } else if (sip->choice == SEQID_GENBANK || sip->choice == SEQID_EMBL || sip->choice == SEQID_DDBJ) {

      tsip = (TextSeqIdPtr) sip->data.ptrvalue;
      if (tsip != NULL && tsip->accession != NULL) {
        if (StringLen (tsip->accession) == 12) {
          is_wgs = TRUE;
          if (StringCmp (tsip->accession + 6, "000000") == 0) {
            wgsaccn = tsip->accession;
            wgsname = tsip->name; /* master accession has 8 zeroes, name has project version plus 6 zeroes */
          }
        } else if (ajp->newSourceOrg && StringLen (tsip->accession) == 6) {
          ch = tsip->accession [0];
          if (ch == 'J' || ch == 'K' || ch == 'L' || ch == 'M') {
            showGBBSource = TRUE;
          }
        }
      }

    } else if (sip->choice == SEQID_GENERAL) {
      dbt = (DbtagPtr) sip->data.ptrvalue;

      /* show GSDB sequence identifier */

      if (dbt != NULL && StringCmp (dbt->db, "GSDB") == 0 && dbt->tag != NULL) {
        gsdbcbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
        if (gsdbcbp != NULL) {
          gsdbcbp->first = first;

          /* string will be created after we know if there are additional comments */

          gsdbid = dbt->tag->id;
          first = FALSE;
        }
      }

    } else if (sip->choice == SEQID_GI) {
      gi = (Int4) sip->data.intvalue;
    }
  }

  /* RefSeq results in allocated comment string */

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &dcontext);
  while (sdp != NULL) {

    uop = (UserObjectPtr) sdp->data.ptrvalue;
    if (uop != NULL) {

      if (! didTPA) {
        str = GetStrForTPA (uop, bsp);
        if (str != NULL) {

          cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
          if (cbp != NULL) {

            cbp->entityID = dcontext.entityID;
            cbp->itemID = dcontext.itemID;
            cbp->itemtype = OBJ_SEQDESC;
            cbp->first = first;
            first = FALSE;

            if (cbp->first) {
              FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
            } else {
              FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
            }

            FFAddOneString (ffstring, str, FALSE, FALSE, TILDE_EXPAND);

            cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12,5, 5, "CC");
            FFRecycleString(ajp, ffstring);
            ffstring = FFGetString(ajp);
          }
          MemFree (str);
          didTPA = TRUE;
        }
      }

      if (! ajp->flags.hideBankItComment) {
        str = GetStrForBankit (uop);
        if (str != NULL) {

          cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
          if (cbp != NULL) {

            cbp->entityID = dcontext.entityID;
            cbp->itemID = dcontext.itemID;
            cbp->itemtype = OBJ_SEQDESC;
            cbp->first = first;
            first = FALSE;

            if (cbp->first) {
              FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
            } else {
              FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
            }

            FFAddOneString (ffstring, str, TRUE, FALSE, TILDE_EXPAND);

            cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12,5, 5, "CC");
            FFRecycleString(ajp, ffstring);
            ffstring = FFGetString(ajp);
          }
          MemFree (str);
        }
      }

      if (! didRefTrack) {
        str = GetStatusForRefTrack (uop);
        if (str != NULL) {

          cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
          if (cbp != NULL) {

            cbp->entityID = dcontext.entityID;
            cbp->itemID = dcontext.itemID;
            cbp->itemtype = OBJ_SEQDESC;
            cbp->first = first;
            first = FALSE;

            if (cbp->first) {
              FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
            } else {
              FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
            }

            FFAddOneString (ffstring, str, FALSE, FALSE, TILDE_EXPAND);

            AddStrForRefTrack (ajp, ffstring, uop);

            cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12,5, 5, "CC");
            FFRecycleString(ajp, ffstring);
            ffstring = FFGetString(ajp);
          }
          /* do not free static str from GetStatusForRefTrack */
          didRefTrack = TRUE;
        }
      }

      if (! didGenome) {
        str = GetStrForGenome (uop, bsp);
        if (str != NULL) {

          cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
          if (cbp != NULL) {

            cbp->entityID = dcontext.entityID;
            cbp->itemID = dcontext.itemID;
            cbp->itemtype = OBJ_SEQDESC;
            cbp->first = first;
            first = FALSE;

            if (cbp->first) {
              FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
            } else {
              FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
            }

            FFAddOneString (ffstring, str, FALSE, FALSE, TILDE_EXPAND);

            cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
            FFRecycleString(ajp, ffstring);
            ffstring = FFGetString(ajp);
          }
          MemFree (str);
          didGenome = TRUE;
        }
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &dcontext);
  }

  if (is_tpa && (! didTPA)) {
    cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
    if (cbp != NULL) {

      cbp->first = first;
      first = FALSE;

      if (cbp->first) {
        FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
      } else {
        FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
      }

      FFAddOneString (ffstring, defTpaString, TRUE, FALSE, TILDE_EXPAND);

      cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
      FFRecycleString(ajp, ffstring);
      ffstring = FFGetString(ajp);
    }
  }

  if (bsp->repr == Seq_repr_delta && bsp->seq_ext_type == 4 && is_wgs) {
    has_gaps = FALSE;
    for (dsp = (DeltaSeqPtr) bsp->seq_ext; dsp; dsp=dsp->next) {
      if (dsp->choice == 2) {
        litp = (SeqLitPtr) dsp->data.ptrvalue;
        if (litp != NULL) {
          if (litp->seq_data == NULL && litp->length > 0) {
            has_gaps = TRUE;
          }
        }
      }
    }
    if (has_gaps) {
      cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
      if (cbp != NULL) {

        cbp->first = first;
        first = FALSE;

        if (cbp->first) {
          FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
        } else {
          FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
        }

        FFAddOneString (ffstring, nsAreGapsString, TRUE, FALSE, TILDE_EXPAND);

        cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
        FFRecycleString(ajp, ffstring);
        ffstring = FFGetString(ajp);
      }
    }
  }

  /* Seq-hist results in allocated comment string */

  hist = bsp->hist;
  if (hist != NULL) {

    if (hist->replaced_by_ids != NULL && hist->replaced_by_date != NULL) {

      okay = TRUE;
      for (sip = hist->replaced_by_ids; sip != NULL; sip = sip->next) {
        if (sip->choice == SEQID_GI) {
          if (gi == (Int4) sip->data.intvalue) {
            okay = FALSE;
          }
        }
      }

      if (okay) {
        cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
        if (cbp != NULL) {
          cbp->first = first;
          first = FALSE;

          if (cbp->first) {
            FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
          } else {
            FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
          }

          AddHistCommentString (ajp, ffstring, "[WARNING] On", "this sequence was replaced by a newer version",
                                hist->replaced_by_date, hist->replaced_by_ids);

          cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
          FFRecycleString(ajp, ffstring);
          ffstring = FFGetString(ajp);
        }
      }
    }

    if (hist->replace_ids != NULL && hist->replace_date != NULL) {

      okay = TRUE;
      for (sip = hist->replace_ids; sip != NULL; sip = sip->next) {
        if (sip->choice == SEQID_GI) {
          if (gi == (Int4) sip->data.intvalue) {
            okay = FALSE;
          }
        }
      }

      if (okay) {
        cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
        if (cbp != NULL) {
          cbp->first = first;
          first = FALSE;

          if (cbp->first) {
            FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
          } else {
            FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
          }

          AddHistCommentString (ajp, ffstring, "On", "this sequence version replaced",
                                hist->replace_date, hist->replace_ids);

          cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
          FFRecycleString(ajp, ffstring);
          ffstring = FFGetString(ajp);
        }
      }
    }

  }

  /* just save IDs for comment, maploc, and region descriptors */

  /*
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_comment, &dcontext);
  while (sdp != NULL) {
    if (sdp->data.ptrvalue != NULL) {
      cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
      if (cbp != NULL) {
        cbp->entityID = dcontext.entityID;
        cbp->itemID = dcontext.itemID;
        cbp->itemtype = OBJ_SEQDESC;
        cbp->first = first;
        first = FALSE;
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_comment, &dcontext);
  }
  */

  /* WGS master comment goes before comment descriptors */

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
  if (sdp != NULL) {

    mip = (MolInfoPtr) sdp->data.ptrvalue;
    if (mip != NULL) {
      if (mip->tech == MI_TECH_wgs) {

        if (wgsname != NULL) {

          cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
          if (cbp != NULL) {

            /*
            cbp->entityID = dcontext.entityID;
            cbp->itemID = dcontext.itemID;
            cbp->itemtype = OBJ_SEQDESC;
            */
            cbp->first = first;
            first = FALSE;

            if (cbp->first) {
              FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
            } else {
              FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
            }

            AddWGSMasterCommentString (ffstring,bsp, wgsaccn, wgsname);

            cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
            FFRecycleString(ajp, ffstring);
            ffstring = FFGetString(ajp);
          }
        }
      }
    }
  }

  if (showGBBSource) {
    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_genbank, &dcontext);
    if (sdp != NULL) {
      gbp = (GBBlockPtr) sdp->data.ptrvalue;
      if (gbp != NULL && (! StringHasNoText (gbp->source))) {
        cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
        if (cbp != NULL) {
          cbp->entityID = dcontext.entityID;
          cbp->itemID = dcontext.itemID;
          cbp->itemtype = OBJ_SEQDESC;
          cbp->first = first;
          first = FALSE;

          if (cbp->first) {
            FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
          } else {
            FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
          }

          FFAddOneString (ffstring, "Original source text: ", FALSE, FALSE, TILDE_EXPAND);
          FFAddOneString (ffstring, gbp->source, TRUE, TRUE, TILDE_EXPAND);

          cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
          FFRecycleString(ajp, ffstring);
          ffstring = FFGetString(ajp);
        }
      }
    }
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_comment, &dcontext);
  while (sdp != NULL) {
    if (sdp->data.ptrvalue != NULL) {
      cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
      if (cbp != NULL) {
        cbp->entityID = dcontext.entityID;
        cbp->itemID = dcontext.itemID;
        cbp->itemtype = OBJ_SEQDESC;
        cbp->first = first;
        first = FALSE;
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_comment, &dcontext);
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_maploc, &dcontext);
  while (sdp != NULL) {
    cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
    if (cbp != NULL) {
      cbp->entityID = dcontext.entityID;
      cbp->itemID = dcontext.itemID;
      cbp->itemtype = OBJ_SEQDESC;
      cbp->first = first;
      first = FALSE;
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_maploc, &dcontext);
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_region, &dcontext);
  while (sdp != NULL) {
    if (sdp->data.ptrvalue != NULL) {
      cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
      if (cbp != NULL) {
        cbp->entityID = dcontext.entityID;
        cbp->itemID = dcontext.itemID;
        cbp->itemtype = OBJ_SEQDESC;
        cbp->first = first;
        first = FALSE;
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_region, &dcontext);
  }

  /* HTGS results in allocated comment string */

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
  if (sdp != NULL) {

    mip = (MolInfoPtr) sdp->data.ptrvalue;
    if (mip != NULL) {
      if (mip->completeness != 0 && is_other) {

        str = GetMolInfoCommentString (bsp, mip);

        if (str != NULL) {
          cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
          if (cbp != NULL) {

            cbp->entityID = dcontext.entityID;
            cbp->itemID = dcontext.itemID;
            cbp->itemtype = OBJ_SEQDESC;
            cbp->first = first;
            first = FALSE;

            if (cbp->first) {
              FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
            } else {
              FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
            }

            FFAddOneString (ffstring, str, TRUE, FALSE, TILDE_EXPAND);

            cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
            FFRecycleString(ajp, ffstring);
            ffstring = FFGetString(ajp);
          }
        }

      }
      if (mip->tech == MI_TECH_htgs_0 ||
          mip->tech == MI_TECH_htgs_1 ||
          mip->tech == MI_TECH_htgs_2) {

        cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
        if (cbp != NULL) {

          /*
          cbp->entityID = dcontext.entityID;
          cbp->itemID = dcontext.itemID;
          cbp->itemtype = OBJ_SEQDESC;
          */
          cbp->first = first;
          first = FALSE;

          if (cbp->first) {
            FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
          } else {
            FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
          }
          
          AddHTGSCommentString (ffstring, bsp, mip);

          cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
          FFRecycleString(ajp, ffstring);
          ffstring = FFGetString(ajp);
        }

      } else {
        str = StringForSeqTech (mip->tech);
        if (! StringHasNoText (str)) {

          cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
          if (cbp != NULL) {

            /*
            cbp->entityID = dcontext.entityID;
            cbp->itemID = dcontext.itemID;
            cbp->itemtype = OBJ_SEQDESC;
            */
            cbp->first = first;
            first = FALSE;

            if (cbp->first) {
              FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
            } else {
              FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
            }

            FFAddTextToString (ffstring, "Method: ", str, NULL, TRUE, FALSE, TILDE_EXPAND);

            cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
            FFRecycleString(ajp, ffstring);
            ffstring = FFGetString(ajp);
          }
        }
      }
    }
  }

  /* add comment features that are full length on appropriate segment */

  parent = awp->parent;
  if (parent == NULL) return;

  sfp = SeqMgrGetNextFeature (parent, NULL, SEQFEAT_COMMENT, 0, &fcontext);
  while (sfp != NULL) {
    if (fcontext.left == awp->from && fcontext.right == awp->to) {
      cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
      if (cbp != NULL) {
        cbp->entityID = fcontext.entityID;
        cbp->itemID = fcontext.itemID;
        cbp->itemtype = OBJ_SEQFEAT;
        cbp->first = first;
        first = FALSE;
      }
    }
    sfp = SeqMgrGetNextFeature (parent, sfp, SEQFEAT_COMMENT, 0, &fcontext);
  }

  if (gsdbcbp != NULL) {

    /* if there were no subsequent comments, do not add period after GSDB id */

    if (cbp == NULL) {
      sprintf (buf, "GSDB:S:%ld", (long) gsdbid);
    } else {
      sprintf (buf, "GSDB:S:%ld.", (long) gsdbid);
    }

    if (gsdbcbp->first) {
      FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
    } else {
      FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
    }

    /* CheckEndPunctuation, ConvertDoubleQuotes, and ExpandTildes already taken into account */

    FFAddOneString (ffstring, buf, FALSE, FALSE, TILDE_IGNORE);

    gsdbcbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
    FFRecycleString(ajp, ffstring);
    ffstring = FFGetString(ajp);
  }

  FFRecycleString(ajp, ffstring);
}

static void AddFeatHeaderBlock (
  Asn2gbWorkPtr awp
)

{
  IntAsn2gbJobPtr ajp;
  BaseBlockPtr    bbp;
  StringItemPtr   ffstring;

  if (awp == NULL) return;
  ajp = awp->ajp;
  if (ajp == NULL) return;

  bbp = Asn2gbAddBlock (awp, FEATHEADER_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  if (awp->format == FTABLE_FMT) return;

  ffstring = FFGetString(ajp);
  if ( ffstring == NULL ) return;

  FFStartPrint (ffstring, awp->format, 0, 12, "FEATURES", 21, 5, 0, "FH", TRUE);

  if (awp->format == EMBL_FMT || awp->format == EMBLPEPT_FMT) {
    FFAddOneString (ffstring, "Key", FALSE, FALSE, TILDE_IGNORE);
    FFAddNChar(ffstring, ' ', 13 , FALSE);
  }

  FFAddOneString (ffstring, "Location/Qualifiers", FALSE, FALSE, TILDE_TO_SPACES);

  if (awp->format == EMBL_FMT || awp->format == EMBLPEPT_FMT) {
    FFAddNewLine(ffstring);
    FFAddNewLine(ffstring);
  }

  bbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 21, 5, 0, "FH");
  FFRecycleString(ajp, ffstring);
}

static Uint2 ComputeSourceHash (
  CharPtr key,
  Uint2 start
)

{
  Uint4  h;
  Uint2  M;
  Uint2  S;

  if (key == NULL) return start;

  M = 101; /* prime key */
  S = 256; /* size of alphabet */

  for (h = start; *key != '\0'; key++) {
    h = (S * h + *key) % M;
  }

  return (Uint2) h;
}

static BaseBlockPtr AddSource (
  Asn2gbWorkPtr awp,
  ValNodePtr PNTR head,
  BioSourcePtr biop,
  CharPtr comment
)

{
  BaseBlockPtr    bbp;
  DbtagPtr        dbt;
  Uint2           hash;
  Int2            idx;
  IntSrcBlockPtr  isp;
  ObjectIdPtr     oip;
  OrgModPtr       omp;
  OrgNamePtr      onp;
  OrgRefPtr       orp;
  SubSourcePtr    ssp;
  CharPtr         str;
  Uint1           subtype;
  Char            tmp [16];
  ValNodePtr      vnp;

  if (awp == NULL || head == NULL || biop == NULL) return NULL;

  bbp = (BaseBlockPtr) MemNew (sizeof (IntSrcBlock));
  if (bbp == NULL) return NULL;
  bbp->blocktype = SOURCEFEAT_BLOCK;
  bbp->section = awp->currsection;

  ValNodeAddPointer (head, 0, bbp);

  isp = (IntSrcBlockPtr) bbp;
  isp->biop = biop;
  isp->is_focus = biop->is_focus;
  if (biop->origin == 5) {
    isp->is_synthetic = TRUE;
  }

  orp = biop->org;
  if (orp == NULL) return bbp;

  isp->orghash = ComputeSourceHash (orp->taxname, 0);
  isp->taxname = orp->taxname;

  hash = 0;
  onp = orp->orgname;
  if (onp != NULL) {
    if (StringICmp (onp->div, "SYN") == 0) {
      isp->is_synthetic = TRUE;
    }
    isp->omp = onp->mod;
    for (omp = onp->mod; omp != NULL; omp = omp->next) {
      subtype = omp->subtype;
      if (subtype == 253) {
        subtype = 35;
      } else if (subtype == 254) {
        subtype = 36;
      } else if (subtype == 255) {
        subtype = 37;
      }
      if (subtype < 38) {
        idx = orgModToSourceIdx [subtype];
        if (idx > 0 && idx < ASN2GNBK_TOTAL_SOURCE) {
          str = asn2gnbk_source_quals [idx].name;
          hash = ComputeSourceHash (str, hash);
          hash = ComputeSourceHash (omp->subname, hash);
        }
      }
    }
  }
  if (comment != NULL) {
    hash = ComputeSourceHash ("note", hash);
    hash = ComputeSourceHash (comment, hash);
  }
  isp->modhash = hash;

  hash = 0;
  for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
    subtype = ssp->subtype;
    if (subtype == 255) {
      subtype = 29;
    }
    if (subtype < 30) {
      idx = subSourceToSourceIdx [subtype];
      if (idx > 0 && idx < ASN2GNBK_TOTAL_SOURCE) {
        str = asn2gnbk_source_quals [idx].name;
        hash = ComputeSourceHash (str, hash);
        hash = ComputeSourceHash (ssp->name, hash);
      }
    }
  }
  isp->subhash = hash;
  isp->ssp = biop->subtype;

  hash = 0;
  for (vnp = orp->db; vnp != NULL; vnp = vnp->next) {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt != NULL) {
      hash = ComputeSourceHash (dbt->db, hash);
      oip = dbt->tag;
      if (oip != NULL) {
        if (oip->str != NULL) {
          hash = ComputeSourceHash (oip->str, hash);
        } else {
          sprintf (tmp, "%ld", (long) oip->id);
          hash = ComputeSourceHash (tmp, hash);
        }
      }
    }
  }
  isp->xrfhash = hash;
  isp->vnp = orp->db;

  return bbp;
}

static int LIBCALLBACK SortSourcesByHash (
  VoidPtr ptr1,
  VoidPtr ptr2
)

{
  Int4            diff;
  IntSrcBlockPtr  isp1;
  IntSrcBlockPtr  isp2;
  ValNodePtr      vnp1;
  ValNodePtr      vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  isp1 = (IntSrcBlockPtr) vnp1->data.ptrvalue;
  isp2 = (IntSrcBlockPtr) vnp2->data.ptrvalue;
  if (isp1 == NULL || isp2 == NULL) return 0;

  if (isp1->is_focus && (! isp2->is_focus)) return -1;
  if (isp2->is_focus && (! isp1->is_focus)) return 1;

  diff = isp1->orghash - isp2->orghash;
  if (diff > 0) return -1;
  if (diff < 0) return 1;

  diff = isp1->xrfhash - isp2->xrfhash;
  if (diff > 0) return -1;
  if (diff < 0) return 1;

  /* sort so that sources with modifiers come first */

  diff = isp1->modhash - isp2->modhash;
  if (diff > 0) return -1;
  if (diff < 0) return 1;

  diff = isp1->subhash - isp2->subhash;
  if (diff > 0) return -1;
  if (diff < 0) return 1;

  /* if all hashes are equal, descriptor comes first */

  if (isp1->is_descriptor && (! isp2->is_descriptor)) {
    return -1;
  } else if (isp2->is_descriptor && (! isp1->is_descriptor)) {
    return 1;
  }

  return 0;
}

static int LIBCALLBACK SortSourcesByPos (
  VoidPtr ptr1,
  VoidPtr ptr2
)

{
  IntSrcBlockPtr  isp1;
  IntSrcBlockPtr  isp2;
  ValNodePtr      vnp1;
  ValNodePtr      vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  isp1 = (IntSrcBlockPtr) vnp1->data.ptrvalue;
  isp2 = (IntSrcBlockPtr) vnp2->data.ptrvalue;
  if (isp1 == NULL || isp2 == NULL) return 0;

  /* descriptor always goes first */

  if (isp1->is_descriptor && (! isp2->is_descriptor)) {
    return -1;
  } else if (isp2->is_descriptor && (! isp1->is_descriptor)) {
    return 1;
  }

  /* feature with smallest left extreme is first */

  if (isp1->left > isp2->left) {
    return 1;
  } else if (isp1->left < isp2->left) {
    return -1;
  }

  /* if same left extreme, shortest source feature is first just for flatfile */

  if (isp1->right > isp2->right) {
    return 1;
  } else if (isp1->right < isp2->right) {
    return -1;
  }

  return 0;
}

/*                                                                   */
/* s_isFuzzyLoc () -- Determines is a location has fuzzy coordinates */
/*                                                                   */

static Boolean s_isFuzzyLoc ( SeqLocPtr pLocation )
{
  SeqIntPtr pIntLocation;

  if (pLocation == NULL)
    return FALSE;

  if (pLocation->choice != SEQLOC_INT)
    return FALSE;

  if (pLocation->data.ptrvalue == NULL)
    return FALSE;

  pIntLocation = (SeqIntPtr) pLocation->data.ptrvalue;

  if ((pIntLocation->if_from != NULL) && (pIntLocation->if_from->choice == 2))
    return TRUE;

  if ((pIntLocation->if_to != NULL) && (pIntLocation->if_to->choice == 2))
    return TRUE;

  return FALSE;
}

static void GetSourcesOnBioseq (
  Asn2gbWorkPtr awp,
  BioseqPtr target,
  BioseqPtr bsp,
  Int4 from,
  Int4 to
)

{
  IntAsn2gbJobPtr    ajp;
  BaseBlockPtr       bbp;
  BioSourcePtr       biop;
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext;
  Boolean            hasNulls;
  Int4               left;
  Boolean            loop = FALSE;
  Int2               idx;
  IntSrcBlockPtr     isp;
  Int4Ptr            ivals;
  SeqLocPtr          newloc;
  Boolean            noLeft;
  Boolean            noRight;
  Int2               numivals;
  Boolean            okay;
  Int4               right;
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp;
  SeqInt             sint;
  SeqIdPtr           sip;
  Boolean            split;
  Int4               start;
  Int4               stop;
  Uint1              strand;
  ValNode            vn;
  ValNodePtr         vnp;

  if (awp == NULL || target == NULL || bsp == NULL) return;
  ajp = awp->ajp;
  if (ajp == NULL) return;

  /* full length loc for descriptors */

  sint.from = 0;
  if (ajp->ajp.slp != NULL) {
    sint.to = SeqLocLen (ajp->ajp.slp) - 1;
  } else {
    sint.to = bsp->length - 1;
  }
  sint.strand = Seq_strand_plus;
  sint.id = SeqIdStripLocus (SeqIdDup (SeqIdFindBest (bsp->id, 0)));
  sint.if_from = NULL;
  sint.if_to = NULL;

  vn.choice = SEQLOC_INT;
  vn.data.ptrvalue = (Pointer) &sint;
  vn.next = NULL;

  /* if SWISS-PROT, may have multiple source descriptors */

  if (ISA_aa (bsp->mol)) {
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_SWISSPROT) {
        loop = TRUE;
      }
    }
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  while (sdp != NULL) {

    /* check if descriptor on part already added on segmented bioseq */

    okay = TRUE;
    for (vnp = awp->srchead; vnp != NULL && okay; vnp = vnp->next) {
      bbp = (BaseBlockPtr) vnp->data.ptrvalue;
      if (bbp != NULL) {
        if (bbp->entityID == dcontext.entityID &&
            bbp->itemID == dcontext.itemID &&
            bbp->itemtype == OBJ_SEQDESC) {
          okay = FALSE;
        }
      }
    }

    if (okay) {
      biop = (BioSourcePtr) sdp->data.ptrvalue;
      bbp = AddSource (awp, &(awp->srchead), biop, NULL);
      if (bbp != NULL) {

        bbp->entityID = dcontext.entityID;
        bbp->itemID = dcontext.itemID;
        bbp->itemtype = OBJ_SEQDESC;

        isp = (IntSrcBlockPtr) bbp;
        isp->loc = SeqLocMerge (target, &vn, NULL, FALSE, TRUE, FALSE);
        isp->left = 0;
        isp->right = bsp->length - 1;
        isp->is_descriptor = TRUE;
      }
    }

    /* if SWISS-PROT, loop through multiple source descriptors */

    if (loop) {
      sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_source, &dcontext);
    } else {
      sdp = NULL;
    }
  }

  SeqIdFree (sint.id);

  if ((! awp->contig) || awp->showconsource) {

    /* features are indexed on parent if segmented */

    bsp = awp->parent;

    sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_BIOSRC, 0, &fcontext);
    while (sfp != NULL) {
      ivals = fcontext.ivals;
      numivals = fcontext.numivals;
      if (ivals != NULL && numivals > 0) {

        idx = (numivals - 1) * 2;
        start = ivals [idx];
        stop = ivals [idx + 1];
        if (stop >= from && stop <= to && (ajp->ajp.slp == NULL || SeqLocCompare (sfp->location, ajp->ajp.slp) > 0)) {

          biop = (BioSourcePtr) sfp->data.value.ptrvalue;
          bbp = AddSource (awp, &(awp->srchead), biop, sfp->comment);
          if (bbp != NULL) {

            bbp->entityID = fcontext.entityID;
            bbp->itemID = fcontext.itemID;
            bbp->itemtype = OBJ_SEQFEAT;

            isp = (IntSrcBlockPtr) bbp;
            if (sfp->location != NULL && sfp->location->choice == SEQLOC_PNT) {
              isp->loc = AsnIoMemCopy ((Pointer) sfp->location,
                                       (AsnReadFunc) SeqLocAsnRead,
                                       (AsnWriteFunc) SeqLocAsnWrite);
            } else if (s_isFuzzyLoc (sfp->location)) {
              isp->loc = AsnIoMemCopy ((Pointer) sfp->location,
                                      (AsnReadFunc) SeqLocAsnRead,
                                      (AsnWriteFunc) SeqLocAsnWrite);
            } else if (SeqLocId(sfp->location) == NULL) {
              isp->loc = AsnIoMemCopy ((Pointer) sfp->location,
                                       (AsnReadFunc) SeqLocAsnRead,
                                       (AsnWriteFunc) SeqLocAsnWrite);
            } else {
              CheckSeqLocForPartial (sfp->location, &noLeft, &noRight);
              hasNulls = LocationHasNullsBetween (sfp->location);
              isp->loc = SeqLocMerge (target, sfp->location, NULL, FALSE, TRUE, hasNulls);
              SetSeqLocPartial (isp->loc, noLeft, noRight);
            }
            isp->left = fcontext.left;
            isp->right = fcontext.right;
            isp->comment = sfp->comment;
            if (ajp->ajp.slp != NULL) {
              sip = SeqIdParse ("lcl|dummy");
              left = GetOffsetInBioseq (ajp->ajp.slp, bsp, SEQLOC_LEFT_END);
              right = GetOffsetInBioseq (ajp->ajp.slp, bsp, SEQLOC_RIGHT_END);
              strand = SeqLocStrand (ajp->ajp.slp);
              split = FALSE;
              newloc = SeqLocReMap (sip, ajp->ajp.slp, isp->loc, 0, FALSE);
              /*
              newloc = SeqLocCopyRegion (sip, isp->loc, bsp, left, right, strand, &split);
              */
              SeqIdFree (sip);
              if (newloc != NULL) {
                A2GBSeqLocReplaceID (newloc, ajp->ajp.slp);
                isp->loc = SeqLocFree (isp->loc);
                isp->loc = newloc;
                isp->left = left;
                isp->right = right;
              }
            }
          }
        }
      }

      sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_BIOSRC, 0, &fcontext);
    }
  }
}

static Boolean LIBCALLBACK GetSourcesOnSeg (
  SeqLocPtr slp,
  SeqMgrSegmentContextPtr context
)

{
  Asn2gbWorkPtr  awp;
  BioseqPtr      bsp;
  Int4           from;
  SeqLocPtr      loc;
  SeqEntryPtr    oldscope;
  SeqEntryPtr    sep;
  SeqIdPtr       sip;
  Int4           to;

  if (slp == NULL || context == NULL) return FALSE;
  awp = (Asn2gbWorkPtr) context->userdata;

  from = context->cumOffset;
  to = from + context->to - context->from;

  sip = SeqLocId (slp);
  if (sip == NULL) {
    loc = SeqLocFindNext (slp, NULL);
    if (loc != NULL) {
      sip = SeqLocId (loc);
    }
  }
  if (sip == NULL) return TRUE;

  /* biosource descriptors only on parts within entity */

  sep = GetTopSeqEntryForEntityID (awp->entityID);
  oldscope = SeqEntrySetScope (sep);
  bsp = BioseqFind (sip);
  SeqEntrySetScope (oldscope);

  if (bsp != NULL) {
    GetSourcesOnBioseq (awp, awp->target, bsp, from, to);
    return TRUE;
  }

  /* if we ever want to fetch remote sources, code goes here */

#if 0
  Uint2          entityID;

  /* may remote fetch genome component if not already in memory */

  bsp = BioseqLockById (sip);

  if (bsp == NULL) return TRUE;

  entityID = ObjMgrGetEntityIDForPointer (bsp);

  if (entityID != awp->entityID) {

    /* if segment not packaged in record, may need to feature index it */

    if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
      SeqMgrIndexFeatures (entityID, NULL);
    }

    /* collect features indexed on the remote bioseq */

    from = 0;
    to = bsp->length - 1;
  }

  GetSourcesOnBioseq (awp, awp->target, bsp, from, to);

  BioseqUnlock (bsp);
#endif

  return TRUE;
}

/* isIdenticalSource() -- Checks to see if two sources are identical */
/*                        by comparing the actual values in the      */
/*                        fields.  This only gets called if the two  */
/*                        sources hashed the same -- it's a double-  */
/*                        check since two non-identical things will  */
/*                        occassionally hash to the same value.      */

static Boolean isIdenticalSource (IntSrcBlockPtr isp1, IntSrcBlockPtr isp2)
{
  OrgModPtr     omp1;
  OrgModPtr     omp2;
  SubSourcePtr  ssp1;
  SubSourcePtr  ssp2;
  ValNodePtr    vnp1;
  ValNodePtr    vnp2;
  ObjectIdPtr   oip1;
  ObjectIdPtr   oip2;
  DbtagPtr      dbt1;
  DbtagPtr      dbt2;

  if (isp1->is_focus != isp2->is_focus)
    return FALSE;

  /* Compare the taxonomy names */

  if (StringICmp(isp1->taxname,isp2->taxname) != 0)
    return FALSE;

  /* Compare the comment */

  if (StringICmp(isp1->comment,isp2->comment) != 0)
    return FALSE;

  /* Compare the org mods */

  omp1 = isp1->omp;
  omp2 = isp2->omp;
  while (omp1 != NULL && omp2 != NULL)
    {
      if (omp1->subtype != omp2->subtype)
        return FALSE;
      if (StringICmp (omp1->subname, omp2->subname) != 0)
        return FALSE;
      omp1 = omp1->next;
      omp2 = omp2->next;
    }

  if (omp1 != NULL || omp2 != NULL)
    return FALSE;

  /* Compare the subtypes */

  ssp1 = isp1->ssp;
  ssp2 = isp2->ssp;

  while (ssp1 != NULL && ssp2 != NULL)
    {
      if (ssp1->subtype != ssp2->subtype)
        return FALSE;
      if (StringICmp(ssp1->name,ssp2->name) != 0)
        return FALSE;
      ssp1 = ssp1->next;
      ssp2 = ssp2->next;
    }

  if (ssp1 != NULL || ssp2 != NULL)
    return FALSE;

  /* Compare the DB tags */

  vnp1 = isp1->vnp;
  vnp2 = isp2->vnp;

  while (vnp1 != NULL && vnp2 != NULL)
    {
      dbt1 = (DbtagPtr) vnp1->data.ptrvalue;
      dbt2 = (DbtagPtr) vnp2->data.ptrvalue;

      if ((dbt1 != NULL) && (dbt2 != NULL)) {
        if (dbt1->db != dbt2->db)
          return FALSE;

        oip1 = dbt1->tag;
        oip2 = dbt2->tag;
        if ((oip1 != NULL) && (oip2 != NULL)) {
          if (oip1->str != NULL) {
            if (StringICmp(oip1->str, oip2->str) != 0)
              return FALSE;
          } else  {
            if (oip1->id != oip2->id)
              return FALSE;
          }
        }
        else if (oip1 != NULL)
          return FALSE;
        else if (oip2 != NULL)
          return FALSE;
      }
      else if (dbt1 != NULL)
        return FALSE;
      else if (dbt2 != NULL)
        return FALSE;

      vnp1 = vnp1->next;
      vnp2 = vnp2->next;
    }

  if (vnp1 != NULL || vnp2 != NULL)
    return FALSE;

  /* If it passed all checks, then they */
  /* are the same, so return true.      */

  return TRUE;
}


static void AddSourceFeatBlock (
  Asn2gbWorkPtr awp
)

{
  IntAsn2gbJobPtr    ajp;
  BaseBlockPtr       bbp;
  BioseqPtr          bsp;
  SeqFeatPtr         cds;
  SeqMgrFeatContext  context;
  BioseqPtr          dna;
  SeqLocPtr          duploc;
  Boolean            excise;
  ValNodePtr         head = NULL;
  IntSrcBlockPtr     isp;
  IntSrcBlockPtr     lastisp;
  IntSrcBlockPtr     descrIsp;
  ValNodePtr         next;
  ValNodePtr         PNTR prev;
  SeqInt             sint;
  SeqLocPtr          slp;
  CharPtr            str;
  BioseqPtr          target;
  ValNode            vn;
  ValNodePtr         vnp;
  Boolean            descHasFocus = FALSE;
  StringItemPtr      ffstring;

  if (awp == NULL) return;
  ajp = awp->ajp;
  if (ajp == NULL) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;

  ffstring = FFGetString(ajp);
  if ( ffstring == NULL ) return;

  /* collect biosources on bioseq */

  awp->srchead = NULL;
  GetSourcesOnBioseq (awp, bsp, bsp, awp->from, awp->to);
  target = bsp;

  if (bsp->repr == Seq_repr_seg) {

    /* collect biosource descriptors on local parts */

    SeqMgrExploreSegments (bsp, (Pointer) awp, GetSourcesOnSeg);
    target = awp->target;
  }

  if (awp->srchead == NULL && ISA_aa (bsp->mol)) {

    /* if protein with no sources, get sources applicable to DNA location of CDS */

    cds = SeqMgrGetCDSgivenProduct (bsp, &context);
    if (cds != NULL) {
      dna = BioseqFindFromSeqLoc (cds->location);
      if (dna != NULL) {
        GetSourcesOnBioseq (awp, dna, dna, context.left, context.right);
        target = dna;
      }
    }
  }

  head = awp->srchead;
  awp->srchead = NULL;

  if (head == NULL) {
    sint.from = 0;
    sint.to = bsp->length - 1;
    sint.strand = Seq_strand_plus;
    sint.id = SeqIdStripLocus (SeqIdDup (SeqIdFindBest (bsp->id, 0)));
    sint.if_from = NULL;
    sint.if_to = NULL;

    vn.choice = SEQLOC_INT;
    vn.data.ptrvalue = (Pointer) &sint;
    vn.next = NULL;

    FFStartPrint (ffstring, awp->format, 5, 21, NULL, 0, 5, 21, "FT", FALSE);
    FFAddOneString(ffstring, "source", FALSE, FALSE, TILDE_IGNORE);
    FFAddNChar(ffstring, ' ', 21 - 5 - StringLen("source"), FALSE);

    str = FlatLoc (ajp, bsp, &vn, (Boolean) (awp->style == MASTER_STYLE));
    if ( GetWWW(ajp) ) {
      FF_www_featloc (ffstring, str);
    } else {
      FFAddOneString (ffstring, str, FALSE, FALSE, TILDE_IGNORE);
    }
    MemFree (str);

    if (ajp->flags.needOrganismQual) {
      FFAddNewLine(ffstring);
      FFAddTextToString (ffstring, "/organism=\"", "unknown", "\"", FALSE, TRUE, TILDE_TO_SPACES);
#ifdef ASN2GNBK_PRINT_UNKNOWN_ORG
    } else {
      FFAddNewLine(ffstring);
      FFAddTextToString (ffstring, "/organism=\"", "unknown", "\"", FALSE, TRUE, TILDE_TO_SPACES);
#endif
    }

    str = FFEndPrint(ajp, ffstring, awp->format, 5, 21, 5, 21, "FT");

    bbp = (BaseBlockPtr) Asn2gbAddBlock (awp, SOURCEFEAT_BLOCK, sizeof (IntSrcBlock));
    if (bbp != NULL) {
      bbp->section = awp->currsection;
      bbp->string = str;
    } else {
      MemFree(str);
    }
    FFRecycleString(ajp, ffstring);
    return;
  }

  /* sort by hash values */

  head = SortValNode (head, SortSourcesByHash);

  /* unique sources, excise duplicates from list */

  prev = &(head);
  vnp = head;
  lastisp = NULL;
  while (vnp != NULL) {
    excise = FALSE;
    next = vnp->next;
    isp = (IntSrcBlockPtr) vnp->data.ptrvalue;
    if (isp->is_descriptor && isp->is_focus)
      descHasFocus = TRUE;
    if (lastisp != NULL) {
      if (isp != NULL) {
        if (lastisp->is_focus == isp->is_focus &&
            lastisp->orghash == isp->orghash &&
            lastisp->xrfhash == isp->xrfhash) {

          /* check for identical modifiers */

          if (lastisp->modhash == isp->modhash &&
              lastisp->subhash == isp->subhash) {

            excise = isIdenticalSource (isp, lastisp);

          /* or modifiers only in lastisp (e.g., on part bioseq) */

          } else if (isp->modhash == 0 && isp->subhash == 0) {
            excise = isIdenticalSource (isp, lastisp);
          }
        }
      }
    }
    if (excise) {
      *prev = vnp->next;
      vnp->next = NULL;

      /* combine locations of duplicate sources */

      if (lastisp != NULL) {
        slp = SeqLocMerge (target, lastisp->loc, isp->loc, FALSE, TRUE, FALSE);
        lastisp->loc = SeqLocFree (lastisp->loc);
        lastisp->loc = slp;
        lastisp->left = MIN (lastisp->left,isp->left);
        lastisp->right = MAX (lastisp->right, isp->right);
      }

      /* and remove duplicate source */

      SeqLocFree (isp->loc);
      MemFree (isp);
      ValNodeFree (vnp);

    } else {

      prev = &(vnp->next);
      lastisp = isp;
    }
    vnp = next;
  }

  /* Sort again, by location this time */

  head = SortValNode (head, SortSourcesByPos);

  /* If the descriptor has a focus, then subtract */
  /* out all the other source locations.          */

  descrIsp = (IntSrcBlockPtr) head->data.ptrvalue; /* Sorted 1st by now */

  if ((descHasFocus) && (! descrIsp->is_synthetic)) {

    vnp = head;
    duploc = AsnIoMemCopy ((Pointer) descrIsp->loc,
                           (AsnReadFunc) SeqLocAsnRead,
                           (AsnWriteFunc) SeqLocAsnWrite);
    vnp = vnp->next;
    while (vnp != NULL) {
      isp = (IntSrcBlockPtr) vnp->data.ptrvalue;
      if (SeqLocAinB (descrIsp->loc, isp->loc) >= 0) {
        vnp = NULL; /* break the chain */
        descrIsp->loc = SeqLocFree (descrIsp->loc);
        descrIsp->loc = duploc;
        duploc = NULL;
      } else {
        descrIsp->loc = SeqLocSubtract (descrIsp->loc, isp->loc);
        vnp = vnp->next;
      }
    }
    descrIsp->left  = SeqLocStart (descrIsp->loc);
    descrIsp->right = SeqLocStop (descrIsp->loc);
    SeqLocFree (duploc);
  }

  /* if features completely subtracted descriptor
     intervals, suppress in release, entrez modes */

  if (descrIsp->loc == NULL && ajp->flags.hideEmptySource && head->next != NULL) {
    vnp = head->next;
    head->next = NULL;
    ValNodeFreeData (head);
    head = vnp;
  }

  /* finally link into blocks for current section */

  ValNodeLink (&(awp->lastblock), head);
  vnp = awp->lastblock;
  if (vnp == NULL) return;
  while (vnp->next != NULL) {
    vnp = vnp->next;
  }

  awp->lastblock = vnp;
  if (awp->blockList == NULL) {
    awp->blockList = vnp;
  }
  FFRecycleString(ajp, ffstring);
}

static void GetFeatsOnCdsProduct (
  SeqFeatPtr sfp,
  BioseqPtr bsp,
  Asn2gbWorkPtr awp
)

{
  FeatBlockPtr       fbp;
  IntFeatBlockPtr    ifp;
  Int4               lastleft;
  Int4               lastright;
  SeqAnnotPtr        lastsap;
  SeqFeatPtr         lastsfp;
  SeqMgrFeatContext  pcontext;
  SeqFeatPtr         prt;
  Boolean            suppress;

  if (sfp == NULL || awp == NULL) return;
  if (bsp == NULL || (! ISA_aa (bsp->mol))) return;

  /* explore mat_peptides, sites, etc. */

  lastsfp = NULL;
  lastsap = NULL;
  lastleft = 0;
  lastright = 0;

  prt = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &pcontext);
  while (prt != NULL) {

    if (pcontext.featdeftype == FEATDEF_REGION ||
        pcontext.featdeftype == FEATDEF_SITE ||
        pcontext.featdeftype == FEATDEF_BOND ||
        pcontext.featdeftype == FEATDEF_mat_peptide_aa ||
        pcontext.featdeftype == FEATDEF_sig_peptide_aa ||
        pcontext.featdeftype == FEATDEF_transit_peptide_aa) {

      if (pcontext.dnaStop >= awp->from && pcontext.dnaStop <= awp->to) {

        /* suppress duplicate features (on protein) */

        suppress = FALSE;
        if (lastsfp != NULL && lastsap != NULL) {
          if (lastsfp->idx.subtype == prt->idx.subtype &&
              lastleft == pcontext.left &&
              lastright == pcontext.right) {
              if (lastsap == pcontext.sap ||
                  (lastsap->desc == NULL && pcontext.sap->desc == NULL)) {
              if (AsnIoMemComp (lastsfp, prt, (AsnWriteFunc) SeqFeatAsnWrite)) {
                suppress = TRUE;
              }
            }
          }
        }

        if (! suppress) {

          fbp = (FeatBlockPtr) Asn2gbAddBlock (awp, FEATURE_BLOCK, sizeof (IntFeatBlock));
          if (fbp != NULL) {

            fbp->entityID = pcontext.entityID;
            fbp->itemID = pcontext.itemID;
            fbp->itemtype = OBJ_SEQFEAT;
            fbp->featdeftype = pcontext.featdeftype;
            ifp = (IntFeatBlockPtr) fbp;
            ifp->mapToNuc = TRUE;
            ifp->mapToProt = FALSE;
            ifp->mapToGen = FALSE;
            ifp->mapToMrna = FALSE;
            ifp->mapToPep = FALSE;
            ifp->firstfeat = awp->firstfeat;
            awp->firstfeat = FALSE;
          }
        }

        lastsfp = prt;
        lastsap = pcontext.sap;
        lastleft = pcontext.left;
        lastright = pcontext.right;

      }
    }
    prt = SeqMgrGetNextFeature (bsp, prt, 0, 0, &pcontext);
  }
}

static Boolean NotEMBLorDDBJ (
  BioseqPtr bsp
)

{
  SeqIdPtr  sip;

  if (bsp == NULL) return TRUE;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_EMBL) return FALSE;
    if (sip->choice == SEQID_DDBJ) return FALSE;
  }
  return TRUE;
}

static Boolean LIBCALLBACK GetFeatsOnBioseq (
  SeqFeatPtr sfp,
  SeqMgrFeatContextPtr fcontext
)

{
  IntAsn2gbJobPtr    ajp;
  Asn2gbSectPtr      asp;
  Asn2gbWorkPtr      awp;
  BioseqPtr          bsp;
  Char               buf [41];
  SeqFeatPtr         cds;
  SeqMgrFeatContext  cdscontext;
  SeqMgrDescContext  dcontext;
  FeatBlockPtr       fbp;
  SeqLocPtr          firstslp;
  GBQualPtr          gbq;
  SeqFeatPtr         gene;
  Int4               gi;
  GeneRefPtr         grp;
  Boolean            juststop = FALSE;
  IntCdsBlockPtr     icp;
  Int2               idx;
  IntFeatBlockPtr    ifp;
  Int4Ptr            ivals;
  Int2               j;
  SeqAnnotPtr        lastsap;
  SeqFeatPtr         lastsfp;
  SeqLocPtr          lastslp;
  SeqLocPtr          newloc;
  Int2               numivals;
  Boolean            okay;
  SeqEntryPtr        oldscope;
  BioseqPtr          parent;
  Boolean            partial5;
  Boolean            partial3;
  ValNodePtr         ppr;
  PubdescPtr         pdp;
  BioseqPtr          prod;
  Boolean            pseudo = FALSE;
  SeqDescrPtr        sdp;
  SeqEntryPtr        sep;
  SeqIdPtr           sip;
  SeqLocPtr          slp;
  Int4               start;
  Int4               stop;
  TextSeqIdPtr       tsip;
  ValNodePtr         vnp;

  if (sfp == NULL || fcontext == NULL) return FALSE;
  awp = (Asn2gbWorkPtr) fcontext->userdata;
  if (awp == NULL) return FALSE;
  ajp = awp->ajp;
  if (ajp == NULL) return FALSE;
  asp = awp->asp;
  if (asp == NULL) return FALSE;
  bsp = asp->bsp;
  if (bsp == NULL) return FALSE;

  if (awp->hideImpFeats && sfp->data.choice == SEQFEAT_IMP) return TRUE;
  if (awp->hideSnpFeats && fcontext->featdeftype == FEATDEF_variation) return TRUE;

  if (fcontext->featdeftype == FEATDEF_PUB ||
      fcontext->featdeftype == FEATDEF_NON_STD_RESIDUE ||
      fcontext->featdeftype == FEATDEF_BIOSRC ||
      fcontext->featdeftype == FEATDEF_RSITE ||
      fcontext->featdeftype == FEATDEF_SEQ) return TRUE;

  if (ajp->flags.validateFeats &&
      (fcontext->featdeftype == FEATDEF_BAD ||
       fcontext->featdeftype == FEATDEF_virion)) {
    return TRUE;
  }

  if (ISA_na (bsp->mol) && fcontext->featdeftype == FEATDEF_HET) return TRUE;

  /* DDBJ does not want to show gene features */

  if (fcontext->seqfeattype == SEQFEAT_GENE && awp->hideGeneFeats) return TRUE;

  /* suppress comment features that are full length */

  if (fcontext->seqfeattype == SEQFEAT_COMMENT &&
      fcontext->left == awp->from && fcontext->right == awp->to) return TRUE;

  ivals = fcontext->ivals;
  numivals = fcontext->numivals;

  /* check to see if last interval is on this awp->from - awp->to range */

  if (ivals != NULL && numivals > 0) {
    idx = (numivals - 1) * 2;
    start = ivals [idx];
    stop = ivals [idx + 1];
    if (stop < awp->from || stop > awp->to) {

      /* may need to map sig_peptide on a different segment */

      if (fcontext->seqfeattype == SEQFEAT_CDREGION) {
        sip = SeqLocIdForProduct (sfp->product);
        bsp = BioseqFind (sip);
        GetFeatsOnCdsProduct (sfp, bsp, awp);
      }

      if (! awp->showAllFeats) return TRUE;

      /* if showing one segment, only show features covering this segment */

      if (fcontext->right < awp->from || fcontext->left > awp->to) return TRUE;

    } else if (fcontext->farloc && NotEMBLorDDBJ (awp->bsp)) {

      /* last interval may not have been mapped to bioseq if far */

      firstslp = NULL;
      lastslp = NULL;

      slp = SeqLocFindNext (sfp->location, NULL);
      while (slp != NULL) {
        if (slp->choice != SEQLOC_NULL) {
          lastslp = slp;
          if (firstslp == NULL) {
            firstslp = slp;
          }
        }
        slp = SeqLocFindNext (sfp->location, slp);
      }

      /* !!! EMBL may have different desired behavior on where to map !!! */

      if (firstslp != NULL && SeqLocStrand (firstslp) == Seq_strand_minus) {
        slp = firstslp;
      } else {
        slp = lastslp;
      }

      if (slp != NULL) {
        sip = SeqLocId (slp);
        if (sip != NULL) {
          bsp = BioseqFindCore (sip);
          if (bsp == NULL || (bsp != awp->parent && bsp != awp->bsp)) {

            return TRUE;
          }
        }
      }
    }
  }

  /* make sure feature is within sublocation */

  if (ajp->ajp.slp != NULL) {
    if (SeqLocCompare (sfp->location, ajp->ajp.slp) == SLC_NO_MATCH) {
      slp = SeqLocMerge (bsp, sfp->location, NULL, FALSE, TRUE, FALSE);
      if (slp == NULL) return TRUE;
      sip = SeqIdParse ("lcl|dummy");
      newloc = SeqLocReMap (sip, ajp->ajp.slp, slp, 0, FALSE);
      SeqIdFree (sip);
      SeqLocFree (slp);
      if (newloc == NULL) return TRUE;
      SeqLocFree (newloc);
    }
  }

  /* suppress duplicate features (on nucleotide) */

  lastsfp = awp->lastsfp;
  lastsap = awp->lastsap;
  if (lastsfp != NULL && lastsap != NULL) {
    if (lastsfp->idx.subtype == sfp->idx.subtype &&
        awp->lastleft == fcontext->left &&
        awp->lastright == fcontext->right) {
        if (lastsap == fcontext->sap ||
            (lastsap->desc == NULL && fcontext->sap->desc == NULL)) {
        if (AsnIoMemComp (lastsfp, sfp, (AsnWriteFunc) SeqFeatAsnWrite)) {
          return TRUE;
        }
      }
    }
  }

  /* if RELEASE_MODE, verify that features have all mandatory qualifiers */

  if (ajp->flags.needRequiredQuals) {
    okay = FALSE;

    switch (fcontext->featdeftype) {

    case FEATDEF_CDS:
      if (ajp->flags.checkCDSproductID) {
        /* non-pseudo CDS must have /product */
        if (sfp->pseudo) {
          pseudo = TRUE;
        }
        grp = SeqMgrGetGeneXref (sfp);
        if (grp == NULL) {
          sep = GetTopSeqEntryForEntityID (ajp->ajp.entityID);
          oldscope = SeqEntrySetScope (sep);
          gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
          SeqEntrySetScope (oldscope);
          if (gene != NULL) {
            grp = (GeneRefPtr) gene->data.value.ptrvalue;
            if (gene->pseudo) {
              pseudo = TRUE;
            }
          }
        }
        if (grp != NULL && grp->pseudo) {
          pseudo = TRUE;
        }
        if (sfp->location != NULL) {
          if (CheckSeqLocForPartial (sfp->location, &partial5, &partial3)) {
            if (partial5 && (! partial3)) {
              if (SeqLocLen (sfp->location) <= 5) {
                juststop = TRUE;
              }
            }
          }
        }
        if (pseudo || juststop) {
          okay = TRUE;
        } else if (sfp->product != NULL) {
          sip = SeqLocIdForProduct (sfp->product);
          if (sip != NULL) {
            if (sip->choice == SEQID_GI && sip->data.intvalue > 0) {
              sep = GetTopSeqEntryForEntityID (ajp->ajp.entityID);
              oldscope = SeqEntrySetScope (sep);
              prod = BioseqFind (sip);
              SeqEntrySetScope (oldscope);
              if (prod != NULL) {
                for (sip = prod->id; sip != NULL; sip = sip->next) {
                  if (sip->choice == SEQID_GENBANK ||
                     sip->choice == SEQID_EMBL ||
                      sip->choice == SEQID_DDBJ ||
                      sip->choice == SEQID_OTHER ||
                      sip->choice == SEQID_PATENT ||
                      sip->choice == SEQID_TPG ||
                      sip->choice == SEQID_TPE ||
                      sip->choice == SEQID_TPD) {
                    tsip = (TextSeqIdPtr) sip->data.ptrvalue;
                    if (tsip != NULL && (! StringHasNoText (tsip->accession))) {
                      if (ValidateAccn (tsip->accession) == 0)
                      okay = TRUE;
                    }
                  }
                }
              } else {
                /* RELEASE_MODE requires that /protein_id is an accession */
                gi = sip->data.intvalue;
                if (GetAccnVerFromServer (gi, buf)) {
                  okay = TRUE;
                } else {
                  sip = GetSeqIdForGI (gi);
                  if (sip != NULL) {
                    okay = TRUE;
                  }
                }
              }
            } else if (sip->choice == SEQID_GENBANK ||
                       sip->choice == SEQID_EMBL ||
                       sip->choice == SEQID_DDBJ ||
                       sip->choice == SEQID_OTHER ||
                       sip->choice == SEQID_PATENT ||
                       sip->choice == SEQID_TPG ||
                       sip->choice == SEQID_TPE ||
                       sip->choice == SEQID_TPD) {
              tsip = (TextSeqIdPtr) sip->data.ptrvalue;
              if (tsip != NULL && (! StringHasNoText (tsip->accession))) {
                if (ValidateAccn (tsip->accession) == 0)
                okay = TRUE;
              }
            }
          }
        }
      } else {
        okay = TRUE;
      }
      if (! okay) {
        ajp->relModeError = TRUE;
      }
      break;

    case FEATDEF_conflict:
      if (sfp->cit == NULL) {
        /* RefSeq allows conflict with accession in comment instead of sfp->cit */
        for (sip = bsp->id; sip != NULL; sip = sip->next) {
          if (sip->choice == SEQID_OTHER) {
            if (! StringHasNoText (sfp->comment)) {
              okay = TRUE;
            }
          }
        }
      }
      /* continue on to old_sequence */
    case FEATDEF_old_sequence:
      /* conflict and old_sequence require a publication printable on the segment */
      vnp = sfp->cit;

      if (vnp != NULL && asp->referenceArray != NULL) {
        for (ppr = vnp->data.ptrvalue; ppr != NULL; ppr = ppr->next) {
          j = MatchRef (ppr, asp->referenceArray, asp->numReferences);
          if (j > 0) {
            okay = TRUE;
            break;
          }
        }
      }
      break;

    case FEATDEF_GENE:
      /* gene requires /gene (later locus_tag suffices instead) */
      grp = (GeneRefPtr) sfp->data.value.ptrvalue;
      if (grp != NULL) {
        if (! StringHasNoText (grp->locus)) {
          okay = TRUE;
        } else if (! StringHasNoText (grp->desc)) {
          okay = TRUE;
        }  else if (! StringHasNoText (grp->locus_tag)) {
          okay = TRUE;
        } else {
          vnp = grp->syn;
          if (vnp != NULL) {
            if (! StringHasNoText (vnp->data.ptrvalue)) {
              okay = TRUE;
            }
          }
        }
      }
      break;

    case FEATDEF_protein_bind:
    case FEATDEF_misc_binding:
      /* protein_bind or misc_binding require FTQUAL_bound_moiety */
      gbq = sfp->qual;
      while (gbq != NULL) {
        if (StringICmp (gbq->qual, "bound_moiety") == 0 && (! StringHasNoText (gbq->val))) {
          okay = TRUE;
          break;
        }
        gbq = gbq->next;
      }
      break;

    case FEATDEF_modified_base:
      /* modified_base requires FTQUAL_mod_base */
      gbq = sfp->qual;
      while (gbq != NULL) {
        if (StringICmp (gbq->qual, "mod_base") == 0 && (! StringHasNoText (gbq->val))) {
          okay = TRUE;
          break;
        }
        gbq = gbq->next;
      }
      break;

    default:
      if (fcontext->featdeftype >= FEATDEF_GENE && fcontext->featdeftype < FEATDEF_MAX) {
        okay = TRUE;
      }
      break;
    }

    if (okay == FALSE) return TRUE;
  }

  /* if RELEASE_MODE, suppress features with location on segmented Bioseq */

  if (ajp->flags.suppressSegLoc) {
    bsp = awp->parent;
    if (bsp != NULL && bsp->repr == Seq_repr_seg) {
      slp = SeqLocFindNext (sfp->location, NULL);
      while (slp != NULL) {
        sip = SeqLocId (slp);
        if (sip != NULL) {
          if (SeqIdIn (sip, bsp->id)) return TRUE;
        }
        slp = SeqLocFindNext (sfp->location, slp);
      }
    }
  }

  awp->lastsfp = sfp;
  awp->lastsap = fcontext->sap;
  awp->lastleft = fcontext->left;
  awp->lastright = fcontext->right;

  if (fcontext->seqfeattype == SEQFEAT_CDREGION) {
    fbp = (FeatBlockPtr) Asn2gbAddBlock (awp, FEATURE_BLOCK, sizeof (IntCdsBlock));
  } else {
    fbp = (FeatBlockPtr) Asn2gbAddBlock (awp, FEATURE_BLOCK, sizeof (IntFeatBlock));
  }
  if (fbp == NULL) return TRUE;

  fbp->entityID = fcontext->entityID;
  fbp->itemID = fcontext->itemID;
  fbp->itemtype = OBJ_SEQFEAT;
  fbp->featdeftype = fcontext->featdeftype;
  ifp = (IntFeatBlockPtr) fbp;
  ifp->mapToNuc = FALSE;
  ifp->mapToProt = FALSE;
  ifp->mapToGen = FALSE;
  ifp->mapToMrna = FALSE;
  ifp->mapToPep = FALSE;
  ifp->firstfeat = awp->firstfeat;
  awp->firstfeat = FALSE;
  awp->featseen = TRUE;

  /* optionally map CDS from cDNA onto genomic */

  if (awp->isGPS && ISA_na (bsp->mol) && awp->copyGpsCdsUp &&
      fcontext->featdeftype == FEATDEF_mRNA) {
    sip = SeqLocIdForProduct (sfp->product);
    bsp = BioseqFind (sip);
    if (bsp != NULL && ISA_na (bsp->mol)) {
      cds = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &cdscontext);
      if (cds != NULL) {
        fbp = (FeatBlockPtr) Asn2gbAddBlock (awp, FEATURE_BLOCK, sizeof (IntCdsBlock));
        if (fbp != NULL) {

          fbp->entityID = cdscontext.entityID;
          fbp->itemID = cdscontext.itemID;
          fbp->itemtype = OBJ_SEQFEAT;
          fbp->featdeftype = cdscontext.featdeftype;
          ifp = (IntFeatBlockPtr) fbp;
          ifp->mapToNuc = FALSE;
          ifp->mapToProt = FALSE;
          ifp->mapToGen = TRUE;
          ifp->mapToMrna = FALSE;
          ifp->mapToPep = FALSE;
          ifp->firstfeat = awp->firstfeat;
          awp->firstfeat = FALSE;
        }
      }
    }
  }

  if (fcontext->seqfeattype != SEQFEAT_CDREGION) return TRUE;

  /* if feature table format, do not get features from protein product */

  if (awp->format == FTABLE_FMT) return TRUE;

  /* if CDS, collect more information from product protein bioseq - may be part */

  sip = SeqLocIdForProduct (sfp->product);
  bsp = BioseqFind (sip);
  if (bsp == NULL || (! ISA_aa (bsp->mol))) return TRUE;

  ifp->isCDS = TRUE;
  icp = (IntCdsBlockPtr) ifp;

  /* first explore pubs to pick up figure and maploc */

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_pub, &dcontext);
  while (sdp != NULL) {
    pdp = (PubdescPtr) sdp->data.ptrvalue;
    if (pdp != NULL) {
      if (icp->fig == NULL) {
        icp->fig = StringSaveNoNull (pdp->fig);
      }
      if (icp->maploc == NULL) {
        icp->maploc = StringSaveNoNull (pdp->maploc);
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_pub, &dcontext);
  }

  /* product may be segmented part, and remaining features are indexed on parent */

  parent = SeqMgrGetParentOfPart (bsp, NULL);
  if (parent != NULL) {
    bsp = parent;
  }

  /* then explore mat_peptides, sites, etc. */

  GetFeatsOnCdsProduct (sfp, bsp, awp);

  return TRUE;
}

static Boolean LIBCALLBACK GetFeatsOnSeg (
  SeqLocPtr slp,
  SeqMgrSegmentContextPtr context
)

{
  Asn2gbWorkPtr  awp;
  BioseqPtr      bsp;
  Uint2          entityID;
  Int4           from;
  SeqLocPtr      loc;
  SeqIdPtr       sip;
  Int4           to;

  if (slp == NULL || context == NULL) return FALSE;
  awp = (Asn2gbWorkPtr) context->userdata;

  from = awp->from;
  to = awp->to;

  sip = SeqLocId (slp);
  if (sip == NULL) {
    loc = SeqLocFindNext (slp, NULL);
    if (loc != NULL) {
      sip = SeqLocId (loc);
    }
  }
  if (sip == NULL) return TRUE;

  /* may want to remote fetch genome component if not already in memory */

  bsp = BioseqLockById (sip);

  if (bsp == NULL) return TRUE;

  entityID = ObjMgrGetEntityIDForPointer (bsp);

  if (entityID != awp->entityID) {

    /* if segment not packaged in record, may need to feature index it */

    if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
      SeqMgrIndexFeatures (entityID, NULL);
    }

    /* collect features indexed on the remote bioseq */

    awp->from = 0;
    awp->to = bsp->length - 1;
  }

  awp->lastsfp = NULL;
  awp->lastsap = NULL;
  awp->lastleft = 0;
  awp->lastright = 0;

  if (context->strand == Seq_strand_minus) {
    SeqMgrExploreFeaturesRev (bsp, (Pointer) awp, GetFeatsOnBioseq, /* awp->slp */ slp, NULL, NULL);
  } else {
    SeqMgrExploreFeatures (bsp, (Pointer) awp, GetFeatsOnBioseq, /* awp->slp */ slp, NULL, NULL);
  }

  /* restore original from and to */

  awp->from = from;
  awp->to = to;

  BioseqUnlock (bsp);

  return TRUE;
}

static void AddFeatureBlock (
  Asn2gbWorkPtr awp
)

{
  IntAsn2gbJobPtr    ajp;
  BioseqPtr          bsp;
  SeqFeatPtr         cds;
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext;
  FeatBlockPtr       fbp;
  SeqFeatPtr         gene;
  IntFeatBlockPtr    ifp;
  Boolean            is_other;
  MolInfoPtr         mip;
  SeqFeatPtr         mrna;
  SeqFeatPtr         prot;
  SeqDescrPtr        sdp;
  SeqIdPtr           sip;
  SeqLocPtr          slp;

  if (awp == NULL) return;
  ajp = awp->ajp;
  if (ajp == NULL) return;
  bsp = awp->parent;
  if (bsp == NULL) return;

  awp->lastsfp = NULL;
  awp->lastsap = NULL;
  awp->lastleft = 0;
  awp->lastright = 0;

  /* optionally map gene from genomic onto cDNA */

  if (awp->isGPS && ISA_na (bsp->mol) && awp->copyGpsGeneDown) {
    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
    if (sdp != NULL && sdp->choice == Seq_descr_molinfo) {
      mip = (MolInfoPtr) sdp->data.ptrvalue;
      if (mip != NULL) {
        if (mip->biomol == MOLECULE_TYPE_MRNA) {
          mrna = SeqMgrGetRNAgivenProduct (bsp, NULL);
          if (mrna != NULL) {
            gene = SeqMgrGetOverlappingGene (mrna->location, &fcontext);
            if (gene != NULL && gene->data.choice == SEQFEAT_GENE) {

              fbp = (FeatBlockPtr) Asn2gbAddBlock (awp, FEATURE_BLOCK, sizeof (IntCdsBlock));
              if (fbp != NULL) {

                fbp->entityID = fcontext.entityID;
                fbp->itemID = fcontext.itemID;
                fbp->itemtype = OBJ_SEQFEAT;
                fbp->featdeftype = fcontext.featdeftype;
                ifp = (IntFeatBlockPtr) fbp;
                ifp->mapToNuc = FALSE;
                ifp->mapToProt = FALSE;
                ifp->mapToGen = FALSE;
                ifp->mapToMrna = TRUE;
                ifp->mapToPep = FALSE;
                ifp->isCDS = TRUE;
                ifp->firstfeat = awp->firstfeat;
                awp->firstfeat = FALSE;
              }
            }
          }
        }
      }
    }
  }

  if (awp->farFeatsSuppress) {

    if (bsp->repr == Seq_repr_seg || bsp->repr == Seq_repr_delta) {

      /* if farFeatsSuppress first collect features on remote segments in MASTER_STYLE */

      SeqMgrExploreSegments (bsp, (Pointer) awp, GetFeatsOnSeg);
    }
  }

  if ((! awp->farFeatsSuppress) || (! awp->featseen)) {

    /* reminder - features on near parts are indexed on segmented Bioseq */

    slp = ajp->ajp.slp;
    if (slp != NULL && SeqLocStrand (slp) == Seq_strand_minus) {
      SeqMgrExploreFeaturesRev (bsp, (Pointer) awp, GetFeatsOnBioseq, awp->slp, NULL, NULL);
    } else {
      SeqMgrExploreFeatures (bsp, (Pointer) awp, GetFeatsOnBioseq, awp->slp, NULL, NULL);
    }
  }

  if (awp->format == GENPEPT_FMT && ISA_aa (bsp->mol)) {
    cds = SeqMgrGetCDSgivenProduct (bsp, &fcontext);
    if (cds != NULL && cds->data.choice == SEQFEAT_CDREGION) {

      fbp = (FeatBlockPtr) Asn2gbAddBlock (awp, FEATURE_BLOCK, sizeof (IntCdsBlock));
      if (fbp != NULL) {

        fbp->entityID = fcontext.entityID;
        fbp->itemID = fcontext.itemID;
        fbp->itemtype = OBJ_SEQFEAT;
        fbp->featdeftype = fcontext.featdeftype;
        ifp = (IntFeatBlockPtr) fbp;
        ifp->mapToNuc = FALSE;
        ifp->mapToProt = TRUE;
        ifp->mapToGen = FALSE;
        ifp->mapToMrna = FALSE;
        ifp->mapToPep = FALSE;
        ifp->isCDS = TRUE;
        ifp->firstfeat = awp->firstfeat;
        awp->firstfeat = FALSE;
      }
    }
    prot = SeqMgrGetPROTgivenProduct (bsp, &fcontext);
    if (prot != NULL && prot->data.choice == SEQFEAT_PROT) {

      is_other = FALSE;
      for (sip = bsp->id; sip != NULL; sip = sip->next) {
        if (sip->choice == SEQID_OTHER) {
          is_other = TRUE;
        }
      }

      /* for RefSeq records or GenBank not release_mode */
      if (is_other || (! ajp->flags.forGbRelease)) {

        fbp = (FeatBlockPtr) Asn2gbAddBlock (awp, FEATURE_BLOCK, sizeof (IntCdsBlock));
        if (fbp != NULL) {

          fbp->entityID = fcontext.entityID;
          fbp->itemID = fcontext.itemID;
          fbp->itemtype = OBJ_SEQFEAT;
          fbp->featdeftype = fcontext.featdeftype;
          ifp = (IntFeatBlockPtr) fbp;
          ifp->mapToNuc = FALSE;
          ifp->mapToProt = FALSE;
          ifp->mapToGen = FALSE;
          ifp->mapToMrna = FALSE;
          ifp->mapToPep = TRUE;
          ifp->firstfeat = awp->firstfeat;
          awp->firstfeat = FALSE;
        }
      }
    }
  }

  if (awp->onlyNearFeats) return;

  if (! awp->farFeatsSuppress) {

    if (bsp->repr == Seq_repr_seg || bsp->repr == Seq_repr_delta) {

      /* if not farFeatsSuppress now collect features on remote segments in MASTER_STYLE */

      SeqMgrExploreSegments (bsp, (Pointer) awp, GetFeatsOnSeg);
    }
  }
}

static void AddWGSBlock (
  Asn2gbWorkPtr awp
)

{
  IntAsn2gbJobPtr    ajp;
  BaseBlockPtr       bbp;
  BioseqPtr          bsp;
  Char               buf [80];
  SeqMgrDescContext  dcontext;
  CharPtr            first = NULL;
  CharPtr            last = NULL;
  ObjectIdPtr        oip;
  SeqDescrPtr        sdp;
  UserFieldPtr       ufp;
  UserObjectPtr      uop;
  StringItemPtr      ffstring;

  if (awp == NULL) return;
  ajp = awp->ajp;
  if ( ajp == NULL ) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;

  if (awp->format == EMBL_FMT || awp->format == EMBLPEPT_FMT) return;

  bbp = Asn2gbAddBlock (awp, WGS_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  ffstring = FFGetString(ajp);
  if ( ffstring == NULL ) return;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &dcontext);
  while (sdp != NULL) {
    uop = (UserObjectPtr) sdp->data.ptrvalue;
    if (uop != NULL) {
      oip = uop->type;
      if (oip != NULL && StringICmp (oip->str, "WGSProjects") == 0) {
        for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
          oip = ufp->label;
          if (oip == NULL || oip->str == NULL || ufp->choice != 1) continue;
          if (StringICmp (oip->str, "WGS_accession_first") == 0) {
            first = (CharPtr) ufp->data.ptrvalue;
          } else if (StringICmp (oip->str, "WGS_accession_last") == 0) {
            last = (CharPtr) ufp->data.ptrvalue;
          }
        }
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &dcontext);
  }

  FFStartPrint (ffstring, awp->format, 0, 12, "WGS", 12, 0, 0, NULL, FALSE);

  if (first != NULL && last != NULL) {

    if ( GetWWW(ajp) ) {
      if (StringCmp (first, last) != 0) {
        FFAddTextToString(ffstring, "<a href=", link_wgs, NULL, FALSE, FALSE, TILDE_IGNORE);
        FFAddTextToString(ffstring, "db=Nucleotide&cmd=Search&term=", first, NULL, FALSE, FALSE, TILDE_IGNORE);
        FFAddTextToString(ffstring, ":", last, "[ACCN]>", FALSE, FALSE, TILDE_IGNORE);
        sprintf (buf, "%s-%s", first, last);
        FFAddOneString (ffstring, buf, FALSE, FALSE, TILDE_TO_SPACES);
        FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_TO_SPACES);
      } else {
        FFAddTextToString(ffstring, "<a href=", link_seq, NULL, FALSE, FALSE, TILDE_IGNORE);
        FFAddTextToString(ffstring, "val=", first, ">", FALSE, FALSE, TILDE_IGNORE);
        sprintf (buf, "%s", first);
        FFAddOneString (ffstring, buf, FALSE, FALSE, TILDE_TO_SPACES);
        FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_TO_SPACES);
      }
    } else {
      if (StringCmp (first, last) != 0) {
        sprintf (buf, "%s-%s", first, last);
        FFAddOneString (ffstring, buf, FALSE, FALSE, TILDE_TO_SPACES);
      } else {
        sprintf (buf, "%s", first);
        FFAddOneString (ffstring, buf, FALSE, FALSE, TILDE_TO_SPACES);
      }
    }
  }

  bbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 0, 0, NULL);
  FFRecycleString(ajp, ffstring);
}

static void AddGenomeBlock (
  Asn2gbWorkPtr awp
)

{
  CharPtr            accn;
  IntAsn2gbJobPtr    ajp;
  BaseBlockPtr       bbp;
  BioseqPtr          bsp;
  SeqMgrDescContext  dcontext;
  Boolean            first = TRUE;
  CharPtr            moltype;
  ObjectIdPtr        oip;
  SeqDescrPtr        sdp;
  UserFieldPtr       ufp;
  UserObjectPtr      uop;
  UserFieldPtr       urf;
  StringItemPtr      ffstring;

  if (awp == NULL) return;
  ajp = awp->ajp;
  if ( ajp == NULL ) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;

  if (awp->format == EMBL_FMT || awp->format == EMBLPEPT_FMT) return;

  bbp = Asn2gbAddBlock (awp, GENOME_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  ffstring = FFGetString(ajp);
  if ( ffstring == NULL ) return;

  FFStartPrint (ffstring, awp->format, 0, 12, "GENOME", 12, 0, 0, NULL, FALSE);

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &dcontext);
  while (sdp != NULL) {
    uop = (UserObjectPtr) sdp->data.ptrvalue;
    if (uop != NULL) {
      oip = uop->type;
      if (oip != NULL && StringICmp (oip->str, "GenomeProject") == 0) {
        for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
          oip = ufp->label;
          if (oip == NULL || oip->str == NULL || ufp->choice != 11) continue;
          if (StringICmp (oip->str, "Chromosome") != 0) continue;
          accn = NULL;
          moltype = NULL;
          for (urf = (UserFieldPtr) ufp->data.ptrvalue; urf != NULL; urf = urf->next) {
            oip = urf->label;
            if (oip == NULL || oip->str == NULL || urf->choice != 1) continue;
            if (StringICmp (oip->str, "accession") == 0) {
              accn = (CharPtr) urf->data.ptrvalue;
            } else if (StringICmp (oip->str, "Moltype") == 0) {
              moltype = (CharPtr) urf->data.ptrvalue;
            }
          }
          if (! StringHasNoText (accn)) {
            if (! first) {
              FFAddNewLine(ffstring);
            }
            first = FALSE;
            FFAddOneString (ffstring, accn, FALSE, FALSE, TILDE_IGNORE);
            if (! StringHasNoText (moltype)) {
              FFAddTextToString (ffstring, " (", moltype, ")", FALSE, FALSE, TILDE_TO_SPACES);
            }
          }
        }
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &dcontext);
  }

  bbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 0, 0, NULL);
  FFRecycleString(ajp, ffstring);
}

static void AddBasecountBlock (
  Asn2gbWorkPtr awp
)

{
  BaseBlockPtr  bbp;

  if (awp == NULL) return;

  bbp = Asn2gbAddBlock (awp, BASECOUNT_BLOCK, sizeof (BaseBlock));
}

static void AddOriginBlock (
  Asn2gbWorkPtr awp
)

{
  IntAsn2gbJobPtr    ajp;
  BaseBlockPtr       bbp;
  BioseqPtr          bsp;
  Char               buf [67];
  SeqMgrDescContext  dcontext;
  GBBlockPtr         gbp;
  SeqDescrPtr        sdp;
  StringItemPtr      ffstring;

  if (awp == NULL) return;
  ajp = awp->ajp;
  if (ajp == NULL) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;

  if (awp->format == EMBL_FMT || awp->format == EMBLPEPT_FMT) return;
  
  ffstring = FFGetString(ajp);
  if ( ffstring == NULL ) return;

  bbp = Asn2gbAddBlock (awp, ORIGIN_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  if (awp->format == GENBANK_FMT || awp->format == GENPEPT_FMT) {

    buf [0] = '\0';

    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_genbank, &dcontext);
    if (sdp != NULL) {
      gbp = (GBBlockPtr) sdp->data.ptrvalue;
      if (gbp != NULL && (! StringHasNoText (gbp->origin))) {
        StringNCpy_0 (buf, gbp->origin, sizeof (buf));
        bbp->entityID = dcontext.entityID;
        bbp->itemID = dcontext.itemID;
        bbp->itemtype = OBJ_SEQDESC;
      }
    }

    FFStartPrint (ffstring, awp->format, 0, 12, "ORIGIN", 12, 0, 0, NULL, FALSE);

    if (! StringHasNoText (buf)) {
      FFAddOneString (ffstring, buf, TRUE, FALSE, TILDE_TO_SPACES);
    }
  }

  bbp->string = FFEndPrint(ajp, ffstring, awp->format, 0, 12, 0, 0, NULL);
  FFRecycleString(ajp, ffstring);
}

#define BASES_PER_BLOCK 1200

static void AddSequenceBlock (
  Asn2gbWorkPtr awp
)

{
  BioseqPtr    bsp;
  Int4         len;
  SeqBlockPtr  sbp;
  Int4         start;
  Int4         stop;

  if (awp == NULL) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;

  if (awp->slp != NULL) {
    len = SeqLocLen (awp->slp);
  } else {
    len = bsp->length;
  }

  /* populate individual sequence blocks for given range */

  for (start = 0; start < len; start += BASES_PER_BLOCK) {
    sbp = (SeqBlockPtr) Asn2gbAddBlock (awp, SEQUENCE_BLOCK, sizeof (SeqBlock));
    if (sbp == NULL) continue;

    sbp->entityID = bsp->idx.entityID;
    sbp->itemID = bsp->idx.itemID;
    sbp->itemtype = OBJ_BIOSEQ;

    stop = start + BASES_PER_BLOCK;
    if (stop >= len) {
      stop = len;
    }

    sbp->start = start;
    sbp->stop = stop;
  }
}

static void AddContigBlock (
  Asn2gbWorkPtr awp
)

{
  BaseBlockPtr  bbp;

  if (awp == NULL) return;

  bbp = Asn2gbAddBlock (awp, CONTIG_BLOCK, sizeof (BaseBlock));
}

static void AddSlashBlock (
  Asn2gbWorkPtr awp
)

{
  BaseBlockPtr  bbp;
  CharPtr str;

  if (awp == NULL) return;

  bbp = Asn2gbAddBlock (awp, SLASH_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  str = MemNew(sizeof(Char) * 4);
  StringNCpy(str, "//\n", 4);

  bbp->string = str;
}

/*--------------------------------------------------------*/
/*                                                        */
/*  s_LocusGetBaseName() -                                */
/*                                                        */
/*--------------------------------------------------------*/

static Boolean s_LocusGetBaseName (BioseqPtr parent, BioseqPtr segment, CharPtr baseName)
{
  Char          parentName[SEQID_MAX_LEN];
  Char          segName[SEQID_MAX_LEN];
  SeqIdPtr      sip;
  TextSeqIdPtr  tsip;
  Char          prefix[5];
  Char          bufTmp[SEQID_MAX_LEN];
  Int2          deleteChars;
  Int2          newLength;
  Int2          i;
  Uint2         segNameLen;

  /* Get the parent Sequence ID */

  parentName [0] = '\0';
  sip = NULL;
  for (sip = parent->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_GENBANK ||
        sip->choice == SEQID_EMBL ||
        sip->choice == SEQID_DDBJ) break;
    if (sip->choice == SEQID_TPG ||
        sip->choice == SEQID_TPE ||
        sip->choice == SEQID_TPD) break;
  }

  if (sip != NULL) {
    tsip = (TextSeqIdPtr) sip->data.ptrvalue;
    if (tsip != NULL && (! StringHasNoText (tsip->name))) {
      StringNCpy_0 (parentName, tsip->name, sizeof (parentName));
    }
  }

  if (StringHasNoText (parentName)) {
    StringNCpy_0 (parentName, baseName, sizeof (parentName));
  }

  /* Get segment id */

  segName [0] = '\0';
  segNameLen = 0;
  sip = NULL;
  for (sip = segment->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_GENBANK ||
        sip->choice == SEQID_EMBL ||
        sip->choice == SEQID_DDBJ) break;
    if (sip->choice == SEQID_TPG ||
        sip->choice == SEQID_TPE ||
        sip->choice == SEQID_TPD) break;
    }

  if (sip != NULL) {
    tsip = (TextSeqIdPtr) sip->data.ptrvalue;
    if (tsip != NULL && (! StringHasNoText (tsip->name))) {
      StringNCpy_0 (segName, tsip->name, sizeof (segName));
      segNameLen = StringLen(segName);
    }
  }

  /* If there's no "SEG_" prefix, then */
  /* just use the parent ID.           */

  StringNCpy_0 (prefix,parentName,sizeof (prefix));
  prefix[4] = '\0';
  if (StringCmp(prefix,"SEG_") != 0)
    {
      StringCpy(baseName,parentName);
      return FALSE;
    }

  /* Otherwise, eliminate the "SEG_" ... */

  StringCpy(bufTmp, &parentName[4]);
  StringCpy(parentName,bufTmp);

  /* ... And calculate a base name */

  if (segNameLen > 0 &&
      (segName[segNameLen-1] == '1') &&
      (StringLen(parentName) == segNameLen) &&
      (parentName[segNameLen-1] == segName[segNameLen-1]))
    {
      deleteChars = 1;
      for (i = segNameLen-2; i >= 0; i--)
    if (parentName[i] == '0')
      deleteChars++;
    else
      break;
      newLength = segNameLen - deleteChars;
      StringNCpy (parentName,segName,newLength); /* not StringNCpy_0 */
      parentName[newLength] = '\0';
    }

  /* Return the base name in the basename parameter */

  StringCpy(baseName,parentName);
  return TRUE;
}

/* ********************************************************************** */

/* DoOneSection builds a single report for one bioseq or segment */

static void DoOneSection (
  BioseqPtr target,
  BioseqPtr parent,
  BioseqPtr bsp,
  BioseqPtr refs,
  SeqLocPtr slp,
  Uint2 seg,
  Int4 from,
  Int4 to,
  Boolean contig,
  Boolean onePartOfSeg,
  Asn2gbWorkPtr awp
)

{
  IntAsn2gbJobPtr      ajp;
  Asn2gbSectPtr        asp;
  SeqMgrBioseqContext  bcontext;
  BaseBlockPtr         PNTR blockArray;
  SeqMgrDescContext    dcontext;
  Boolean              hasRefs;
  Int4                 i;
  IntAsn2gbSectPtr     iasp;
  Boolean              isRefSeq = FALSE;
  MolInfoPtr           mip;
  Boolean              nsgenome = FALSE;
  Int4                 numBlocks;
  Int4                 numsegs = 0;
  SeqDescrPtr          sdp;
  SeqIdPtr             sip;
  TextSeqIdPtr         tsip;
  ValNodePtr           vnp;
  Boolean              wgsmaster = FALSE;
  Boolean              wgstech = FALSE;

  if (target == NULL || parent == NULL || bsp == NULL || awp == NULL) return;
  ajp = awp->ajp;
  if (ajp == NULL) return;

  if (awp->mode == RELEASE_MODE && awp->style == CONTIG_STYLE) {
    if (bsp->repr == Seq_repr_seg) {
    } else if (bsp->repr == Seq_repr_delta && (! DeltaLitOnly (bsp))) {
    } else return;
  }

  if (ajp->flags.suppressLocalID) {
    sip = SeqIdSelect (bsp->id, fasta_order, NUM_SEQID);
    if (sip == NULL || sip->choice == SEQID_LOCAL) return;
  }

  if (seg == 0) {
    awp->basename[0] = '\0';
  } else if (seg == 1) {
    s_LocusGetBaseName (parent, bsp, awp->basename);
  }

  asp = Asn2gbAddSection (awp);
  if (asp == NULL) return;

  numsegs = awp->partcount;
  if (numsegs == 0 && SeqMgrGetBioseqContext (parent, &bcontext)) {
    numsegs = bcontext.numsegs;
  }

  /* set working data fields */

  awp->asp = asp;

  awp->target = target;
  awp->parent = parent;
  awp->bsp = bsp;
  awp->refs = refs;
  awp->slp = slp;
  awp->seg = seg;
  awp->numsegs = numsegs;
  awp->from = from;
  awp->to = to;
  awp->contig = contig;

  awp->firstfeat = TRUE;
  awp->featseen = FALSE;

  /* initialize empty blockList for this section */

  awp->blockList = NULL;
  awp->lastblock = NULL;

  /* and store section data into section fields */

  asp->target = target;
  asp->bsp = bsp;
  asp->slp = slp;
  asp->seg = seg;
  asp->numsegs = numsegs;
  asp->from = from;
  asp->to = to;

  iasp = (IntAsn2gbSectPtr) asp;
  iasp->spp = NULL;

  asp->blockArray = NULL;
  asp->numBlocks = 0;

  /* WGS master and NS_ virtual records treated differently */

  if (bsp->repr == Seq_repr_virtual) {

    /* check for certain ID types */

    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_GENBANK ||
          sip->choice == SEQID_EMBL ||
          sip->choice == SEQID_DDBJ) {
        tsip = (TextSeqIdPtr) sip->data.ptrvalue;
        if (tsip != NULL && tsip->accession != NULL) {
          if (StringLen (tsip->accession) == 12) {
            if (StringCmp (tsip->accession + 6, "000000") == 0) {
              wgsmaster = TRUE;
            }
          }
        }
     } else if (sip->choice == SEQID_OTHER) {
        tsip = (TextSeqIdPtr) sip->data.ptrvalue;
        if (tsip != NULL && tsip->accession != NULL) {
          if (StringNICmp (tsip->accession, "NC_", 3) == 0) {
            wgsmaster = TRUE;
          } else if (StringNICmp (tsip->accession, "NS_", 3) == 0) {
            nsgenome = TRUE;
          } else if (StringNICmp (tsip->accession, "NZ_", 3) == 0) {
            if (StringLen (tsip->accession) == 15) {
              if (StringCmp (tsip->accession + 9, "000000") == 0) {
                wgsmaster = TRUE;
              }
            }
          }
        }
      }
    }

    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
    if (sdp != NULL) {
      mip = (MolInfoPtr) sdp->data.ptrvalue;
      if (mip != NULL && mip->tech == MI_TECH_wgs) {
        wgstech = TRUE;
      }
    }
  }

  /* start exploring and populating paragraphs */

  if (awp->format == FTABLE_FMT) {

    AddFeatHeaderBlock (awp);
    AddFeatureBlock (awp);

  } else {

    AddLocusBlock (awp);

    if (awp->format == GENBANK_FMT || awp->format == GENPEPT_FMT) {

      AddDeflineBlock (awp);
      AddAccessionBlock (awp);

      if (ISA_aa (bsp->mol)) {
        /*
        AddPidBlock (awp);
        */
      }

      AddVersionBlock (awp);

      if (ISA_aa (bsp->mol)) {
        AddDbsourceBlock (awp);
      }

    } else if (awp->format == EMBL_FMT || awp->format == EMBLPEPT_FMT) {

      AddAccessionBlock (awp);

      if (ISA_na (bsp->mol)) {
        AddVersionBlock (awp);
      }

      if (ISA_aa (bsp->mol)) {
        /* AddPidBlock (awp); */
        /* AddDbsourceBlock (awp); */
      }

      AddDateBlock (awp);

      AddDeflineBlock (awp);
    }

    AddKeywordsBlock (awp);

    if (awp->format == GENBANK_FMT || awp->format == GENPEPT_FMT) {
      AddSegmentBlock (awp, onePartOfSeg);
    }

    AddSourceBlock (awp);
    AddOrganismBlock (awp);

    /* !!! RELEASE_MODE should check return value of AddReferenceBlock !!! */

    hasRefs = AddReferenceBlock (awp);
    if (! hasRefs) {
      if (ajp->flags.needAtLeastOneRef) {
        /* RefSeq does not require a publication */
        for (sip = bsp->id; sip != NULL; sip = sip->next) {
          if (sip->choice == SEQID_OTHER) {
            isRefSeq = TRUE;
          }
        }
        if (! isRefSeq) {
          awp->failed = TRUE;
        }
      }
    }
    AddCommentBlock (awp);
    AddPrimaryBlock (awp);

    AddFeatHeaderBlock (awp);
    AddSourceFeatBlock (awp);

    if (wgsmaster && wgstech) {

      AddWGSBlock (awp);

    } else if (nsgenome) {

      AddGenomeBlock (awp);

    } else if (contig) {

      if (awp->showconfeats) {
        AddFeatureBlock (awp);
      }
      AddContigBlock (awp);

    } else {
 
      AddFeatureBlock (awp);

      if (ISA_na (bsp->mol)) {
        AddBasecountBlock (awp);
      }
      AddOriginBlock (awp);

      AddSequenceBlock (awp);
    }

    AddSlashBlock (awp);
  }

  /* allocate block array for this section */

  numBlocks = ValNodeLen (awp->blockList);
  asp->numBlocks = numBlocks;

  if (numBlocks > 0) {
    blockArray = (BaseBlockPtr PNTR) MemNew (sizeof (BaseBlockPtr) * (numBlocks + 1));
    asp->blockArray = blockArray;

    if (blockArray != NULL) {
      for (vnp = awp->blockList, i = 0; vnp != NULL; vnp = vnp->next, i++) {
        blockArray [i] = (BaseBlockPtr) vnp->data.ptrvalue;
      }
    }
  }

  /* free blockList, but leave data, now pointed to by blockArray elements */

  awp->blockList = ValNodeFree (awp->blockList);
  awp->lastblock = NULL;

  (awp->currsection)++;
}


/*
the following functions handle various kinds of input, all calling
DoOneSection once for each component that gets its own report
*/

static Boolean LIBCALLBACK Asn2Seg (
  SeqLocPtr slp,
  SeqMgrSegmentContextPtr context
)

{
  Asn2gbWorkPtr  awp;
  BioseqPtr      bsp = NULL;
  Uint2          entityID;
  Int4           from;
  SeqLocPtr      loc;
  BioseqPtr      parent;
  SeqIdPtr       sip;
  Int4           to;

  if (slp == NULL || context == NULL) return FALSE;
  awp = (Asn2gbWorkPtr) context->userdata;

  parent = context->parent;

  from = context->cumOffset;
  to = from + context->to - context->from;

  sip = SeqLocId (slp);
  if (sip == NULL) {
    loc = SeqLocFindNext (slp, NULL);
    if (loc != NULL) {
      sip = SeqLocId (loc);
    }
  }
  if (sip == NULL) return TRUE;

  /* may remote fetch genome component if not already in memory */

  bsp = BioseqLockById (sip);

  if (bsp == NULL) return TRUE;

  entityID = ObjMgrGetEntityIDForPointer (bsp);

  if (entityID != awp->entityID) {

    /* if segment not packaged in record, may need to feature index it */

    if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
      SeqMgrIndexFeatures (entityID, NULL);
    }

    /* collect features indexed on the remote bioseq */

    parent = bsp;
    from = 0;
    to = bsp->length - 1;
  }

  if (bsp->repr != Seq_repr_virtual) {
    (awp->seg)++;
    DoOneSection (bsp, parent, bsp, bsp, /* slp */ NULL, awp->seg, from, to, FALSE, FALSE, awp);
  }

  BioseqUnlock (bsp);

  return TRUE;
}

static Int4 CountRealParts (
  SeqLocPtr slp_head
)

{
  SeqIdPtr   id;
  Int4       numparts;
  BioseqPtr  part;
  SeqIdPtr   sip;
  SeqLocPtr  slp;

  numparts = 0;
  for (slp = (SeqLocPtr) slp_head; slp != NULL; slp = slp->next) {
    sip = SeqLocId (slp);
    if (sip == NULL) continue;
    if (sip->choice == SEQID_GI) {
      part = BioseqFind (sip);
      if (part == NULL) continue;
      for (id = part->id; id != NULL; id = id->next) {
        if (id->choice == SEQID_GIBBSQ ||
            id->choice == SEQID_GIBBMT ||
            id->choice == SEQID_GIIM) break;
      }
      if (id != NULL && part->repr == Seq_repr_virtual) continue;
    }
    numparts++;
  }
  return numparts;
}

typedef struct findseg {
  BioseqPtr  bsp;
  Uint2      seg;
} FindSeg, PNTR FindSegPtr;

static Boolean LIBCALLBACK FindSegForPart (
  SeqLocPtr slp,
  SeqMgrSegmentContextPtr context
)

{
  FindSegPtr  fsp;
  BioseqPtr   bsp = NULL;
  SeqLocPtr   loc;
  SeqIdPtr    sip;

  if (slp == NULL || context == NULL) return TRUE;
  fsp = (FindSegPtr) context->userdata;

  sip = SeqLocId (slp);
  if (sip == NULL) {
    loc = SeqLocFindNext (slp, NULL);
    if (loc != NULL) {
      sip = SeqLocId (loc);
    }
  }
  if (sip == NULL) return TRUE;

  bsp = BioseqFind (sip);
  if (bsp == NULL) return TRUE;

  if (bsp->repr != Seq_repr_virtual) {
    (fsp->seg)++;
  }

  if (bsp != fsp->bsp) return TRUE;

  return FALSE;
}

static void DoOneBioseq (
  BioseqPtr bsp,
  Pointer userdata
)

{
  IntAsn2gbJobPtr       ajp;
  Asn2gbWorkPtr         awp;
  SeqMgrSegmentContext  context;
  Boolean               contig = FALSE;
  Int4                  from;
  FindSeg               fs;
  SeqEntryPtr           oldscope;
  BioseqPtr             parent;
  Boolean               segmented = FALSE;
  SeqEntryPtr           sep;
  Int4                  to;

  if (bsp == NULL) return;
  awp = (Asn2gbWorkPtr) userdata;
  if (awp == NULL) return;
  ajp = awp->ajp;
  if (ajp == NULL) return;

  /* return if molecule not right for format */

  if (ISA_na (bsp->mol)) {
    if (ajp->format == GENPEPT_FMT || ajp->format == EMBLPEPT_FMT) return;
  } else if (ISA_aa (bsp->mol)) {
    if (ajp->format == GENBANK_FMT || ajp->format == EMBL_FMT) return;
  }

  if (awp->style == SEGMENT_STYLE) {
    segmented = TRUE;
  }
  if (awp->style == CONTIG_STYLE) {
    contig = TRUE;
  }

  awp->partcount = 0;

  if (bsp->repr == Seq_repr_seg && awp->style == NORMAL_STYLE) {

    /* if bsp followed by parts set, then do not default to contig style */

    if (SegHasParts (bsp)) {
      segmented = TRUE;
      contig = FALSE;

      if (bsp->seq_ext_type == 1) {

        /* count only non-virtual parts */

        sep = GetTopSeqEntryForEntityID (awp->entityID);
        oldscope = SeqEntrySetScope (sep);
        awp->partcount = CountRealParts ((SeqLocPtr) bsp->seq_ext);
        SeqEntrySetScope (oldscope);
      }
    } else {
      segmented = FALSE;
      contig = TRUE;
    }
  }
  if (bsp->repr == Seq_repr_delta && awp->style == NORMAL_STYLE) {
    if (! DeltaLitOnly (bsp)) {
      contig = TRUE;
    }
  }

  if (bsp->repr == Seq_repr_seg) {

    /* this is a segmented bioseq */

    if (segmented) {

      /* show all segments individually */

      awp->seg = 0;
      SeqMgrExploreSegments (bsp, (Pointer) awp, Asn2Seg);

    } else {

      /* show as single bioseq */

      parent = bsp;
      from = 0;
      to = bsp->length - 1;

      DoOneSection (parent, parent, bsp, parent, ajp->ajp.slp, 0, from, to, contig, FALSE, awp);
    }

  } else if (bsp->repr == Seq_repr_raw ||
             bsp->repr == Seq_repr_const ||
             bsp->repr == Seq_repr_delta ||
             bsp->repr == Seq_repr_virtual) {

    parent = SeqMgrGetParentOfPart (bsp, &context);
    if (parent != NULL) {

      /* this is a part of an indexed segmented bioseq */

      from = context.cumOffset;
      to = from + context.to - context.from;

      s_LocusGetBaseName (parent, bsp, awp->basename);

      fs.bsp = bsp;
      fs.seg = 0;
      SeqMgrExploreSegments (parent, (Pointer) &fs, FindSegForPart);
      awp->showAllFeats = TRUE;

      DoOneSection (bsp, parent, bsp, parent, ajp->ajp.slp, fs.seg, from, to, contig, TRUE, awp);

    } else {

      /* this is a regular non-segmented bioseq */

      parent = bsp;
      from = 0;
      to = bsp->length - 1;

      DoOneSection (bsp, parent, bsp, parent, ajp->ajp.slp, 0, from, to, contig, FALSE, awp);
    }
  }
}

static void DoBioseqSetList (
  SeqEntryPtr seq_set,
  Asn2gbWorkPtr awp
)

{
  BioseqSetPtr  bssp;
  SeqEntryPtr   sep;

  if (seq_set == NULL || awp == NULL) return;

  /* iterate rather than recurse unless multiple nested sets > nuc-prot */

  for (sep = seq_set; sep != NULL; sep = sep->next) {

    if (IS_Bioseq_set (sep)) {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      if (bssp == NULL) continue;

      if (bssp->_class == BioseqseqSet_class_genbank ||
          bssp->_class == BioseqseqSet_class_mut_set ||
          bssp->_class == BioseqseqSet_class_pop_set ||
          bssp->_class == BioseqseqSet_class_phy_set ||
          bssp->_class == BioseqseqSet_class_eco_set ||
          bssp->_class == BioseqseqSet_class_wgs_set ||
          bssp->_class == BioseqseqSet_class_gen_prod_set) {

        /* if popset within genbank set, for example, recurse */

        DoBioseqSetList (bssp->seq_set, awp);

        continue;
      }
    }

    /* at most nuc-prot set, so do main bioseqs that fit the format */

    VisitSequencesInSep (sep, (Pointer) awp, VISIT_MAINS, DoOneBioseq);
  }
}

static void DoOneBioseqSet (
  SeqEntryPtr sep,
  Asn2gbWorkPtr awp
)

{
  BioseqSetPtr  bssp;

  if (sep == NULL || awp == NULL) return;

  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;

    if (bssp->_class == BioseqseqSet_class_genbank ||
        bssp->_class == BioseqseqSet_class_mut_set ||
        bssp->_class == BioseqseqSet_class_pop_set ||
        bssp->_class == BioseqseqSet_class_phy_set ||
        bssp->_class == BioseqseqSet_class_eco_set ||
        bssp->_class == BioseqseqSet_class_wgs_set ||
        bssp->_class == BioseqseqSet_class_gen_prod_set) {

      /* this is a pop/phy/mut/eco set, catenate separate reports */

      DoBioseqSetList (bssp->seq_set, awp);

      return;
    }
  }

  /* at most nuc-prot set, so do main bioseqs that fit the format */

  VisitSequencesInSep (sep, (Pointer) awp, VISIT_MAINS, DoOneBioseq);
}



/* FormatFeatureblockQuals should not be called directly,
   except from FormatFeatureBlock.  It performs no input
   validation.  (perhaps it should?) */

static void FormatFeatureBlockQuals (
  StringItemPtr    ffstring,
  IntAsn2gbJobPtr  ajp,
  Asn2gbFormatPtr  afp,
  Asn2gbSectPtr    asp,
  BioseqPtr        bsp,
  Uint1            featdeftype,
  ValNodePtr       gene_syn,
  CharPtr          lasttype,
  SeqLocPtr        location,
  BioseqPtr        prod,
  CharPtr          protein_pid_g,
  QualValPtr       qvp,
  Int4             left,
  Int4             right,
  Uint1            strand,
  SeqFeatPtr       sfp,
  BioseqPtr        target,
  Boolean          is_other,
  Boolean          is_journalscan,
  Boolean          is_gps
)

{
  Boolean            add_period;
  CharPtr            ascii;
  Int2               ascii_len;
  Boolean            at_end = FALSE;
  ByteStorePtr       bs;
  Char               buf[80];
  Choice             cbaa;
  CodeBreakPtr       cbp;
  Char               ch;
  Uint1              choice;
  Uint1              code = Seq_code_ncbieaa;
  Int4               gi;
  Boolean            hadProtDesc = FALSE;
  DbtagPtr           dbt;
  Int4               exp_ev;
  GBQualPtr          gbq;
  Int2               i;
  Uint1              idx;
  Boolean            isTRNA;
  Int2               j;
  Uint1              jdx;
  Int4               len;
  SeqLocPtr          newloc;
  CharPtr            notestr;
  Char               numbuf [32];
  Int2               numcodons;
  Int2               numsyns;
  ObjectIdPtr        oip;
  Boolean            okay;
  Boolean            only_digits;
  ValNodePtr         ppr;
  CharPtr            prefix;
  CharPtr            protein_seq = NULL;
  size_t             prtlen;
  CharPtr            ptr;
  Uint1              residue;
  SeqCodeTablePtr    sctp;
  Int4               sec_str;
  Uint1              seqcode;
  Char               seqid [50];
  SeqIntPtr          sintp;
  SeqIdPtr           sip;
  SeqLocPtr          slp;
  Boolean            split;
  SeqPortPtr         spp;
  CharPtr            start;
  CharPtr            str;
  Boolean            suppress_period;
  CharPtr            tmp;
  tRNAPtr            trna;
  UserObjectPtr      uop;
  ValNodePtr         vnp;
  StringItemPtr      unique;

  unique = FFGetString(ajp);
  if ( unique == NULL ) return;

  for (i = 0, idx = feat_qual_order [i]; idx != 0; i++, idx = feat_qual_order [i]) {

    lasttype = NULL;
    switch (asn2gnbk_featur_quals [idx].qualclass) {

      case Qual_class_ignore :
        break;

      case Qual_class_string :
        if (! StringHasNoText (qvp [idx].str)) {
          FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=",
            FALSE, TRUE, TILDE_TO_SPACES);
          FFAddTextToString(ffstring, "\"", qvp [idx].str, "\"",
            FALSE, TRUE, TILDE_TO_SPACES);
          FFAddOneChar(ffstring, '\n', FALSE);
        }
        break;

      case Qual_class_tilde :
        if (! StringHasNoText (qvp [idx].str)) {
          FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=",
            FALSE, TRUE, TILDE_EXPAND);
          FFAddTextToString(ffstring, "\"", qvp [idx].str, "\"",
            FALSE, TRUE, TILDE_EXPAND);
          FFAddOneChar(ffstring, '\n', FALSE);
        }
        break;

      case Qual_class_product :
        if (StringHasNoText (qvp [idx].str) ||
            (ajp->flags.dropIllegalQuals &&
             (! AllowedValQual (featdeftype, FTQUAL_product)))) break;
          FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=",
            FALSE, TRUE, TILDE_TO_SPACES);
          FFAddTextToString(ffstring, "\"", qvp [idx].str, "\"",
            FALSE, TRUE, TILDE_TO_SPACES);
          FFAddOneChar(ffstring, '\n', FALSE);
        break;

      case Qual_class_sgml :
        if (! StringHasNoText (qvp [idx].str)) {
          if (is_journalscan) {
            ascii_len = Sgml2AsciiLen (qvp [idx].str);
            start = ascii = MemNew ((size_t) (10 + ascii_len));
            if (start != NULL) {
              ascii = Sgml2Ascii (qvp [idx].str, ascii, ascii_len + 1);

              FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=",
                FALSE, TRUE, TILDE_TO_SPACES);
              FFAddTextToString(ffstring, "\"", start, "\"",
                FALSE, TRUE, TILDE_TO_SPACES);
              FFAddOneChar(ffstring, '\n', FALSE);

              MemFree (start);
            }
          } else {
              FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=",
                FALSE, TRUE, TILDE_TO_SPACES);
              FFAddTextToString(ffstring, "\"", qvp[idx].str, "\"",
                FALSE, TRUE, TILDE_TO_SPACES);
              FFAddOneChar(ffstring, '\n', FALSE);
          }
        }
        break;

      case Qual_class_boolean :
        if (qvp [idx].ble) {
          FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "\n",
                            FALSE, TRUE, TILDE_IGNORE);
        }
        break;

      case Qual_class_int :
        if (qvp [idx].num > 0) {
          if (idx == FTQUAL_transl_table) {
            sprintf (numbuf, "%ld", (long) qvp [idx].num);
            FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=",
              FALSE, TRUE, TILDE_IGNORE);
            FF_www_gcode (ajp, ffstring, numbuf);
          } else {
            sprintf (numbuf, "%ld", (long) qvp [idx].num);
            FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=",
              FALSE, TRUE, TILDE_IGNORE);
            FFAddTextToString(ffstring, NULL, numbuf, NULL,
              FALSE, TRUE, TILDE_IGNORE);
          }
          FFAddOneChar(ffstring, '\n', FALSE);
        }
        break;

      case Qual_class_evidence :
        exp_ev = qvp [idx].num;
        if (exp_ev > 0 && exp_ev <= 2) {
          FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=",
              FALSE, TRUE, TILDE_IGNORE);
          FFAddOneString(ffstring, evidenceText [exp_ev], FALSE, TRUE, TILDE_IGNORE);
          FFAddOneChar(ffstring, '\n', FALSE);
        }
        break;

      case Qual_class_valnode :
        for (vnp = qvp[idx].vnp; vnp != NULL; vnp = vnp->next) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (str != NULL) {
              FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=",
                  FALSE, TRUE, TILDE_TO_SPACES);
              FFAddTextToString(ffstring, "\"", str, "\"",
                  FALSE, TRUE, TILDE_TO_SPACES);
              FFAddOneChar(ffstring, '\n', FALSE);
          }
        }
        break;

      case Qual_class_gene_syn :
        for (vnp = qvp[idx].vnp; vnp != NULL; vnp = vnp->next) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (str != NULL) {
              FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=",
                  FALSE, TRUE, TILDE_TO_SPACES);
              FFAddTextToString(ffstring, "\"", str, "\"",
                  FALSE, TRUE, TILDE_TO_SPACES);
              FFAddOneChar(ffstring, '\n', FALSE);
          }
        }
        break;

      case Qual_class_EC_valnode :
        for (vnp = qvp [idx].vnp; vnp != NULL; vnp = vnp->next) {
          str = (CharPtr) vnp->data.ptrvalue;
          okay = TRUE;

          if (str == NULL) continue;

          if (ajp->flags.dropIllegalQuals) {
            tmp = str;
            while (*tmp != '\0' && *tmp == '\"')
              tmp++;
            for (; *tmp != '\0' && *tmp != '\"'; tmp++) {
              if (!IS_DIGIT(*tmp) && *tmp != '.' && *tmp != '-') {
                okay = FALSE;
              }
            }
          }
          if (!okay) continue;

          FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=",
                FALSE, TRUE, TILDE_TO_SPACES);
          FFAddOneChar(ffstring, '\"', FALSE);
          FF_AddECnumber(ajp, ffstring, str);
          FFAddOneChar(ffstring, '\"', FALSE);
          FFAddOneChar(ffstring, '\n', FALSE);
        }
        break;

      case Qual_class_EC_quote :
        gbq = qvp [idx].gbq;
        if (gbq == NULL || (ajp->flags.dropIllegalQuals && (! AllowedValQual (featdeftype, idx)))) break;
        if (lasttype == NULL) {
          lasttype = gbq->qual;
        }
        while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
          okay = TRUE;
          if (gbq->val == NULL) {
            okay = FALSE;
          }

          if (ajp->flags.dropIllegalQuals && okay) {
            tmp = gbq->val;
            while (*tmp != '\0' && *tmp == '\"')
              tmp++;
            for (; *tmp != '\0' && *tmp != '\"'; tmp++) {
              if (!IS_DIGIT(*tmp) && *tmp != '.' && *tmp != '-') {
                okay = FALSE;
              }
            }
          }

          if (StringHasNoText (gbq->val)) {
            okay = FALSE;
          }

          if (okay) {
            FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=",
                  FALSE, TRUE, TILDE_TO_SPACES);
            FFAddOneChar(ffstring, '\"', FALSE);
            if (!StringIsJustQuotes (gbq->val)) {
              FF_AddECnumber (ajp, ffstring, gbq->val);
            }
            FFAddOneChar(ffstring, '\"', FALSE);
            FFAddOneChar(ffstring, '\n', FALSE);
          }
          gbq = gbq->next;
        }
        break;

      case Qual_class_quote :
        gbq = qvp [idx].gbq;
        if (gbq == NULL || (ajp->flags.dropIllegalQuals && (! AllowedValQual (featdeftype, idx)))) break;
        if (lasttype == NULL) {
          lasttype = gbq->qual;
        }
        while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
          if (! StringHasNoText (gbq->val)) {
            FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=\"",
                  FALSE, TRUE, TILDE_IGNORE);
            if (!StringIsJustQuotes (gbq->val)) {
              FFAddOneString(ffstring, gbq->val, FALSE, TRUE, TILDE_IGNORE);
            }
            FFAddOneChar(ffstring, '\"', FALSE);
            FFAddOneChar(ffstring, '\n', FALSE);
          }
          gbq = gbq->next;
        }
        break;

      case Qual_class_noquote :
        gbq = qvp [idx].gbq;
        if (gbq == NULL || (ajp->flags.dropIllegalQuals && (! AllowedValQual (featdeftype, idx)))) break;
        if (lasttype == NULL) {
          lasttype = gbq->qual;
        }
        while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
          if (! StringHasNoText (gbq->val)) {
            FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=",
                  FALSE, TRUE, TILDE_IGNORE);
            FFAddOneString(ffstring, gbq->val, FALSE, TRUE, TILDE_TO_SPACES);
            FFAddOneChar(ffstring, '\n', FALSE);
          }
          gbq = gbq->next;
        }
        break;

      case Qual_class_label :
        gbq = qvp [idx].gbq;
        if (gbq == NULL || (ajp->flags.dropIllegalQuals && (! AllowedValQual (featdeftype, idx)))) break;
        if (lasttype == NULL) {
          lasttype = gbq->qual;
        }
        while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
          if (! StringHasNoText (gbq->val)) {
            if (ajp->flags.checkQualSyntax) { /* single token, not just numeric */
              str = gbq->val;
              ch = *str;
              only_digits = TRUE;
              while (ch != '\0') {
                if (IS_WHITESP (ch)) break; /* only single token allowed */
                if (! IS_DIGIT (ch)) {
                  only_digits = FALSE;
                }
                str++;
                ch = *str;
              }
              if (only_digits) break; /* must not be just numeric */
            }
            FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=",
                  FALSE, TRUE, TILDE_IGNORE);
            FFAddOneString(ffstring, gbq->val, FALSE, TRUE, TILDE_TO_SPACES);
            FFAddOneChar(ffstring, '\n', FALSE);
          }
          gbq = gbq->next;
        }
        break;

      case Qual_class_number :
        gbq = qvp [idx].gbq;
        if (gbq == NULL || (ajp->flags.dropIllegalQuals && (! AllowedValQual (featdeftype, idx)))) break;

        if (ajp->flags.checkQualSyntax) {
          str = gbq->val;

          if ( StringHasNoText (str) )
            break;
          while (!IS_WHITESP (*str) && *str != '\0')
            str++;
          if (! StringHasNoText (str) )
            break;
        }

        if (lasttype == NULL) {
          lasttype = gbq->qual;
        }
        while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
          if (! StringHasNoText (gbq->val)) {
            FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=",
                  FALSE, TRUE, TILDE_IGNORE);
            FFAddOneString(ffstring, gbq->val, FALSE, TRUE, TILDE_TO_SPACES);
            FFAddOneChar(ffstring, '\n', FALSE);
          }
          gbq = gbq->next;
        }
        break;

      case Qual_class_paren :
        gbq = qvp [idx].gbq;
        if (gbq == NULL || (ajp->flags.dropIllegalQuals && (! AllowedValQual (featdeftype, idx)))) break;
        if (lasttype == NULL) {
          lasttype = gbq->qual;
        }
        while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
          if (! StringHasNoText (gbq->val)) {
            tmp = StringSave (gbq->val);
            str = tmp;
            len = StringLen (str);
            if (len > 1 && *str == '(' && str [len - 1] == ')' &&
                StringChr (str, ',') != NULL) {
              str++;
              while (! StringHasNoText (str)) {
                ptr = StringChr (str, ',');
                if (ptr == NULL) {
                  ptr = StringChr (str, ')');
                }
                if (ptr != NULL) {
                  *ptr = '\0';
                  ptr++;
                }
                FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=",
                  FALSE, TRUE, TILDE_IGNORE);
                FFAddOneString(ffstring, str, FALSE, TRUE, TILDE_TO_SPACES);
                FFAddOneChar(ffstring, '\n', FALSE);
                str = ptr;
              }
            } else {
              FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=",
                  FALSE, TRUE, TILDE_IGNORE);
              FFAddOneString(ffstring, str, FALSE, TRUE, TILDE_TO_SPACES);
              FFAddOneChar(ffstring, '\n', FALSE);
            }
            MemFree (tmp);
          }
          gbq = gbq->next;
        }
        break;

      case Qual_class_rpt :
        gbq = qvp [idx].gbq;
        if (gbq == NULL || (ajp->flags.dropIllegalQuals && (! AllowedValQual (featdeftype, idx)))) break;

        if (lasttype == NULL) {
          lasttype = gbq->qual;
        }
        while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
          if (! StringHasNoText (gbq->val)) {
            tmp = StringSave (gbq->val);
            str = tmp;
            len = StringLen (str);
            if (len > 1 && *str == '(' && str [len - 1] == ')' &&
                StringChr (str, ',') != NULL) {
              str++;
              while (! StringHasNoText (str)) {
                ptr = StringChr (str, ',');
                if (ptr == NULL) {
                  ptr = StringChr (str, ')');
                }
                if (ptr != NULL) {
                  *ptr = '\0';
                  ptr++;
                }
                if ((! ajp->flags.checkQualSyntax) || (StringInStringList (str, validRptString))) {
                  FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=",
                                    FALSE, TRUE, TILDE_IGNORE);
                  FFAddOneString(ffstring, str, FALSE, TRUE, TILDE_TO_SPACES);
                  FFAddOneChar(ffstring, '\n', FALSE);
                }
                str = ptr;
              }
            } else {
              if ((! ajp->flags.checkQualSyntax) || (StringInStringList (str, validRptString))) {
                FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=",
                                  FALSE, TRUE, TILDE_IGNORE);
                FFAddOneString(ffstring, str, FALSE, TRUE, TILDE_TO_SPACES);
                FFAddOneChar(ffstring, '\n', FALSE);
              }
            }
            MemFree (tmp);
          }
          gbq = gbq->next;
        }
        break;

      case Qual_class_rpt_unit :
        gbq = qvp [idx].gbq;
        if (gbq == NULL || (ajp->flags.dropIllegalQuals && (! AllowedValQual (featdeftype, idx)))) break;
        if (lasttype == NULL) {
          lasttype = gbq->qual;
        }

        /* in release_mode, must be of the form 123..4567 or a single-token label,
           or (technically illegal but common) letters and semicolons */

        while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
          if (! StringHasNoText (gbq->val)) {
            tmp = StringSave (gbq->val);
            str = tmp;
            len = StringLen (str);
            if (len > 1 && *str == '(' && str [len - 1] == ')' /* && StringChr (str, ',') != NULL */) {
              str++;
              while (! StringHasNoText (str)) {
                ptr = StringChr (str, ',');
                if (ptr == NULL) {
                  ptr = StringChr (str, ')');
                }
                if (ptr != NULL) {
                  *ptr = '\0';
                  ptr++;
                }
                if ((! ajp->flags.checkQualSyntax) || (ValidateRptUnit (str))) {
                  FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=",
                                    FALSE, TRUE, TILDE_IGNORE);
                  FFAddOneString(ffstring, str, FALSE, TRUE, TILDE_TO_SPACES);
                  FFAddOneChar(ffstring, '\n', FALSE);
                }
                str = ptr;
              }
            } else {
              if ((! ajp->flags.checkQualSyntax) || (ValidateRptUnit (str))) {
                FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=",
                                  FALSE, TRUE, TILDE_IGNORE);
                FFAddOneString(ffstring, str, FALSE, TRUE, TILDE_TO_SPACES);
                FFAddOneChar(ffstring, '\n', FALSE);
              }
            }
            MemFree (tmp);
          }
          gbq = gbq->next;
        }
        break;

      case Qual_class_replace :
        gbq = qvp [idx].gbq;
        if (gbq == NULL || (ajp->flags.dropIllegalQuals && (! AllowedValQual (featdeftype, idx)))) break;
        if (lasttype == NULL) {
          lasttype = gbq->qual;
        }
        while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
          FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=",
                            FALSE, TRUE, TILDE_IGNORE);
          FFAddOneChar(ffstring, '\"', FALSE);
          if (!StringHasNoText (gbq->val)) {
            FFAddOneString(ffstring, gbq->val, FALSE, TRUE, TILDE_TO_SPACES);
          }
          FFAddOneChar(ffstring, '\"', FALSE);
          FFAddOneChar(ffstring, '\n', FALSE);
          gbq = gbq->next;
        }
        break;

      case Qual_class_consplice :
        gbq = qvp [idx].gbq;
        if (gbq == NULL || (ajp->flags.dropIllegalQuals && (! AllowedValQual (featdeftype, idx)))) break;

        if (ajp->flags.checkQualSyntax && (! StringInStringList (gbq->val, validConsSpliceString)) ) {
          break;
        }

        if (lasttype == NULL) {
          lasttype = gbq->qual;
        }
        while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
          if (! StringHasNoText (gbq->val)) {
            FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=",
                              FALSE, TRUE, TILDE_IGNORE);
            FFAddOneString(ffstring, gbq->val, FALSE, TRUE, TILDE_TO_SPACES);
            FFAddOneChar(ffstring, '\n', FALSE);
          }
          gbq = gbq->next;
        }
        break;

      case Qual_class_site :
        if (! StringHasNoText (qvp [idx].str)) {
          FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=",
                            FALSE, TRUE, TILDE_IGNORE);
          FFAddTextToString(ffstring, "\"", qvp[idx].str, "\"", FALSE, TRUE, TILDE_TO_SPACES);
          FFAddOneChar(ffstring, '\n', FALSE);
        }
        break;

      case Qual_class_bond :
        if (! StringHasNoText (qvp [idx].str)) {
          FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=",
                            FALSE, TRUE, TILDE_IGNORE);
          FFAddTextToString(ffstring, "\"", qvp[idx].str, "\"", FALSE, TRUE, TILDE_TO_SPACES);
          FFAddOneChar(ffstring, '\n', FALSE);
        }
        break;

      case Qual_class_L_R_B :
        gbq = qvp [idx].gbq;
        if (gbq == NULL || (ajp->flags.dropIllegalQuals && (! AllowedValQual (featdeftype, idx)))) break;

        if (ajp->flags.checkQualSyntax && (! StringInStringList (gbq->val, validLRBString)) ) {
          break;
        }

        if (lasttype == NULL) {
          lasttype = gbq->qual;
        }
        while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
          if (! StringHasNoText (gbq->val)) {
            FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=",
                              FALSE, TRUE, TILDE_IGNORE);
            FFAddOneString(ffstring, gbq->val, FALSE, TRUE, TILDE_TO_SPACES);
            FFAddOneChar(ffstring, '\n', FALSE);
          }
          gbq = gbq->next;
        }
        break;

      case Qual_class_sec_str :
        sec_str = qvp [idx].num;
        if (sec_str > 0 && sec_str <= 3) {
          FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=",
                            FALSE, TRUE, TILDE_IGNORE);
          FFAddTextToString(ffstring, "\"", secStrText[sec_str], "\"", 
                            FALSE, FALSE, TILDE_IGNORE);
          FFAddOneChar(ffstring, '\n', FALSE);
        }
        break;

      case Qual_class_seq_loc :
        slp = qvp [idx].slp;
        if (slp != NULL) {
          FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=",
                            FALSE, TRUE, TILDE_IGNORE);
          str = FlatLoc (ajp, target, slp, /* ajp->masterStyle */ FALSE);
          FFAddTextToString(ffstring, "\"", str, "\"", 
                            FALSE, TRUE, TILDE_TO_SPACES);
          FFAddOneChar(ffstring, '\n', FALSE);
          MemFree (str);
        }
        break;

      case Qual_class_code_break :
        cbp = qvp [idx].cbp;
        seqcode = 0;
        sctp = NULL;
        while (cbp != NULL) {
          cbaa = cbp->aa;
          switch (cbaa.choice) {
            case 1 :
              seqcode = Seq_code_ncbieaa;
              break;
            case 2 :
              seqcode = Seq_code_ncbi8aa;
              break;
            case 3 :
              seqcode = Seq_code_ncbistdaa;
              break;
            default :
              break;
          }
          if (seqcode != 0) {
            sctp = SeqCodeTableFind (seqcode);
            if (sctp != NULL) {
              slp = NULL;
              while ((slp = SeqLocFindNext (cbp->loc, slp)) != NULL) {
                str = NULL;
                if (ajp->ajp.slp != NULL) {
                  sip = SeqIdParse ("lcl|dummy");
                  split = FALSE;
                  newloc = SeqLocReMap (sip, ajp->ajp.slp, slp, 0, FALSE);
 
                  SeqIdFree (sip);
                  if (newloc != NULL) {
                    A2GBSeqLocReplaceID (newloc, ajp->ajp.slp);
                    str = FlatLoc (ajp, target, newloc, ajp->masterStyle);
                    SeqLocFree (newloc);
                  }
                } else {
                  str = FlatLoc (ajp, target, slp, ajp->masterStyle);
                }
                if (str != NULL) {
                  residue = cbaa.value.intvalue;
                  ptr = Get3LetterSymbol (ajp, seqcode, sctp, residue);
                  if (ptr == NULL) {
                    ptr = "OTHER";
                  }
                  FFAddOneString(ffstring, "/transl_except=", FALSE, FALSE, TILDE_IGNORE);
                  FFAddTextToString(ffstring, "(pos:", str, ",", FALSE, FALSE, TILDE_IGNORE);
                  FFAddTextToString(ffstring, "aa:", ptr, ")", FALSE, FALSE, TILDE_IGNORE);
                  FFAddOneChar(ffstring, '\n', FALSE);
                }
                MemFree (str);
              }
            }
          }
          cbp = cbp->next;
        }
        break;

      case Qual_class_anti_codon :
        slp = qvp [FTQUAL_anticodon].slp;
        newloc = NULL;
        if (ajp->ajp.slp != NULL) {
          sip = SeqIdParse ("lcl|dummy");
          split = FALSE;
          newloc = SeqLocReMap (sip, ajp->ajp.slp, slp, 0, FALSE);
          /*
          newloc = SeqLocCopyRegion (sip, slp, bsp, left, right, strand, &split);
          */
          SeqIdFree (sip);
          slp = newloc;
          if (newloc != NULL) {
            A2GBSeqLocReplaceID (newloc, ajp->ajp.slp);
          }
        }
        if (slp != NULL && slp->choice == SEQLOC_INT) {
          sintp = (SeqIntPtr) slp->data.ptrvalue;
          if (sintp != NULL) {
            str = qvp [FTQUAL_trna_aa].str;
            if (! StringHasNoText (str)) {
              sprintf(numbuf, "%ld", (long) sintp->from + 1);
              FFAddTextToString(ffstring, "/anticodon=(pos:", numbuf, "..",
                                FALSE, FALSE, TILDE_IGNORE);
              sprintf(numbuf, "%ld", (long) sintp->to + 1);
              FFAddTextToString(ffstring, NULL, numbuf, ",",
                                FALSE, FALSE, TILDE_IGNORE);
              FFAddTextToString(ffstring, "aa:", str, ")",
                                FALSE, FALSE, TILDE_IGNORE);
              FFAddOneChar(ffstring, '\n', FALSE);
            }
          }
        }
        if (newloc != NULL) {
          SeqLocFree (newloc);
        }
        break;

      case Qual_class_codon :
        gbq = qvp [idx].gbq;
        if (gbq == NULL || (ajp->flags.dropIllegalQuals && (! AllowedValQual (featdeftype, idx)))) break;
        if (lasttype == NULL) {
          lasttype = gbq->qual;
        }
        while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
          if (! StringHasNoText (gbq->val)) {
            FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals[idx].name, "=",
                              FALSE, FALSE, TILDE_IGNORE);
            FFAddOneString(ffstring, gbq->val, FALSE, FALSE, TILDE_TO_SPACES);
            FFAddOneChar(ffstring, '\n', FALSE);
          }
          gbq = gbq->next;
        }
        break;

      case Qual_class_pubset :
        vnp = qvp [idx].vnp;
        if (vnp != NULL && asp->referenceArray != NULL) {
          for (ppr = vnp->data.ptrvalue; ppr != NULL; ppr = ppr->next) {
            j = MatchRef (ppr, asp->referenceArray, asp->numReferences);
            if (j > 0) {
              sprintf (numbuf, "%d", (int) j);
              FFAddTextToString(ffstring, "/citation=[", numbuf, "]",
                                FALSE, TRUE, TILDE_TO_SPACES);
              FFAddOneChar(ffstring, '\n', FALSE);
            }
          }
        }
        break;

      case Qual_class_db_xref :
        for (vnp = qvp [idx].vnp; vnp != NULL; vnp = vnp->next) {
          buf [0] = '\0';
          dbt = (DbtagPtr) vnp->data.ptrvalue;
          if (dbt != NULL && (! StringHasNoText (dbt->db))) {
            oip = dbt->tag;
            if (oip != NULL) {

              okay = TRUE;
              if (StringCmp (dbt->db, "PID") == 0 || StringCmp (dbt->db, "GI") == 0) {
                okay = FALSE;
              }
              if (ajp->flags.dropBadDbxref) {
                /* if RELEASE_MODE, drop unknown dbtag */
                okay = FALSE;
                for (j = 0; legalDbXrefs [j] != NULL; j++) {
                  if (StringCmp (dbt->db, legalDbXrefs [j]) == 0) {
                    okay = TRUE;
                  }
                }
                if (! okay) {
                  if (is_gps || is_other) {
                    for (j = 0; legalRefSeqDbXrefs [j] != NULL; j++) {
                      if (StringCmp (dbt->db, legalRefSeqDbXrefs [j]) == 0) {
                        okay = TRUE;
                      }
                    }
                  }
                }
              }

              if (okay) {
                if (! StringHasNoText (oip->str)) {
                  if (StringLen (oip->str) < 80) {
                    sprintf (buf, "%s", oip->str);
                  }
                } else {
                  sprintf (buf, "%ld", (long) oip->id);
                }
              }
            }
          }
          if (! StringHasNoText (buf)) {
            if (StringICmp (buf, protein_pid_g) != 0) {
              /* already sorted and uniqued by BasicSeqEntryCleanup, per feature */
              if (StringICmp (dbt->db, "LocusID") == 0 || StringICmp (dbt->db, "InterimID") == 0) {
                if (FFStringSearch (ffstring, dbt->db, 0) >= 0) {
                  okay = FALSE;
                }
              }
              if (okay) {
                FFAddOneString(ffstring, "/db_xref=\"", FALSE, FALSE, TILDE_IGNORE);
                FF_www_db_xref(ajp, ffstring, dbt->db, buf);
                FFAddOneString(ffstring, "\"\n", FALSE, FALSE, TILDE_IGNORE);
              }
            }
          }
        }
        break;

      
      case Qual_class_seq_id :
        sip = qvp [idx].sip;
        if (sip != NULL) {
          /* should always be found above for protein_id or transcript_id
          prod = BioseqFind (sip);
          */
          if (prod != NULL) {
            choice = 0;
            for (sip = prod->id; sip != NULL; sip = sip->next) {
              if (sip->choice == SEQID_GENBANK ||
                  sip->choice == SEQID_EMBL ||
                  sip->choice == SEQID_DDBJ ||
                  sip->choice == SEQID_OTHER ||
                  sip->choice == SEQID_TPG ||
                  sip->choice == SEQID_TPE ||
                  sip->choice == SEQID_TPD) {
                choice = sip->choice;
                if (SeqIdWrite (sip, seqid, PRINTID_TEXTID_ACC_VER, sizeof (seqid)) != NULL) {
                  FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=\"",
                                    FALSE, FALSE, TILDE_IGNORE);
                  FF_www_protein_id(ajp, ffstring, seqid);
                  FFAddOneString(ffstring, "\"\n", FALSE, FALSE, TILDE_IGNORE);
                }
              } else if (sip->choice == SEQID_GI) {
                if (choice == 0) {
                  sprintf (seqid, "%ld", (long) sip->data.intvalue);
                  FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=\"",
                                    FALSE, FALSE, TILDE_IGNORE);
                  FF_www_protein_id(ajp, ffstring, seqid);
                  FFAddOneString(ffstring, "\"\n", FALSE, FALSE, TILDE_IGNORE);
                }
                sprintf (seqid, "%ld", (long) sip->data.intvalue);
                FFAddOneString(ffstring, "/db_xref=\"", FALSE, FALSE, TILDE_IGNORE);
                FF_www_db_xref(ajp, ffstring, "GI", seqid);
                FFAddOneString(ffstring, "\"\n", FALSE, FALSE, TILDE_IGNORE);
              } else if (sip->choice == SEQID_GENERAL) {
                dbt = (DbtagPtr) sip->data.ptrvalue;
                if (dbt != NULL && StringCmp (dbt->db, "PID") == 0) {
                  /*
                  oip = dbt->tag;
                  if (oip != NULL) {
                    if (! StringHasNoText (oip->str)) {
                      sprintf (seqid, "PID:%s", oip->str);
                      NewContLine ();
                      gb_AddString ("/db_xref=\"", seqid, "\"", FALSE, TRUE, TILDE_TO_SPACES);
                    }
                  }
                  */
                }
              }
            }
          } else {
            if (sip->choice == SEQID_GI) {
              gi = sip->data.intvalue;
              if (GetAccnVerFromServer (gi, seqid)) {
                if ((! ajp->flags.dropIllegalQuals) || ValidateAccn (seqid) == 0) {
                  FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=\"",
                                    FALSE, FALSE, TILDE_IGNORE);
                  FF_www_protein_id(ajp, ffstring, seqid);
                  FFAddOneString(ffstring, "\"\n", FALSE, FALSE, TILDE_IGNORE);
                } else {
                  ajp->relModeError = TRUE;
                }
              } else {
                sip = GetSeqIdForGI(gi);
                if (sip != NULL && SeqIdWrite (sip, seqid, PRINTID_TEXTID_ACC_VER, sizeof (seqid)) != NULL) {
                  if ((! ajp->flags.dropIllegalQuals) || ValidateAccn (seqid) == 0) {
                    FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=\"",
                                    FALSE, FALSE, TILDE_IGNORE);
                    FF_www_protein_id(ajp, ffstring, seqid);
                    FFAddOneString(ffstring, "\"\n", FALSE, FALSE, TILDE_IGNORE);
                  } else {
                    ajp->relModeError = TRUE;
                  }
                } else if (! ajp->flags.dropIllegalQuals) {
                  sprintf (seqid, "%ld", (long) gi);
                  FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=\"",
                                    FALSE, FALSE, TILDE_IGNORE);
                  FF_www_protein_id(ajp, ffstring, seqid);
                  FFAddOneString(ffstring, "\"\n", FALSE, FALSE, TILDE_IGNORE);
                } else {
                  ajp->relModeError = TRUE;
                }
              }

              sprintf (seqid, "%ld", (long) gi);
              FFAddOneString(ffstring, "/db_xref=\"", FALSE, FALSE, TILDE_IGNORE);
              FF_www_db_xref(ajp, ffstring, "GI", seqid);
              FFAddOneString(ffstring, "\"\n", FALSE, FALSE, TILDE_IGNORE);
            } else if (SeqIdWrite (sip, seqid, PRINTID_TEXTID_ACC_VER, sizeof (seqid)) != NULL) {
              if ((! ajp->flags.dropIllegalQuals) || ValidateAccn (seqid) == 0) {
                FFAddTextToString(ffstring, "/", asn2gnbk_featur_quals [idx].name, "=\"",
                                  FALSE, FALSE, TILDE_IGNORE);
                FF_www_protein_id(ajp, ffstring, seqid);
                FFAddOneString(ffstring, "\"\n", FALSE, FALSE, TILDE_IGNORE);
              } else {
                ajp->relModeError = TRUE;
              }

              gi = GetGIForSeqId (sip);
              if (gi > 0) {
                sprintf (seqid, "%ld", (long) gi);
                FFAddOneString(ffstring, "/db_xref=\"", FALSE, FALSE, TILDE_IGNORE);
                FF_www_db_xref(ajp, ffstring, "GI", seqid);
                FFAddOneString(ffstring, "\"\n", FALSE, FALSE, TILDE_IGNORE);
              }
            }
          }
        }
        break;

      case Qual_class_translation :
        if (qvp [idx].ble) {
          if ((prod == NULL && ajp->transIfNoProd) || ajp->alwaysTranslCds) {
            bs = ProteinFromCdRegionEx (sfp, TRUE, FALSE);
            if (bs != NULL) {
              str = BSMerge (bs, NULL);
              bs = BSFree (bs);
              if (str != NULL) {
                ptr = str;
                ch = *ptr;
                while (ch != '\0') {
                  *ptr = TO_UPPER (ch);
                  ptr++;
                  ch = *ptr;
                }
                prtlen = StringLen (str);
                if (prtlen > 1) {
                   if (str [prtlen - 1] == '*') {
                     str [prtlen - 1] = '\0';
                   }
                }
                if (! StringHasNoText (str)) {
                  FFAddTextToString(ffstring, "/translation=\"", str, "\"", 
                                    FALSE, TRUE, TILDE_TO_SPACES);
                  FFAddOneChar(ffstring, '\n', FALSE);
                }
                MemFree (str);
              }
            } else {
              ajp->relModeError = TRUE;
            }
          } else if (prod != NULL) {
            len = SeqLocLen (sfp->product);
            if (len > 0) {
              if (SeqLocStart (location) == 0 || SeqLocStop (location) == bsp->length - 1) {
                at_end = TRUE;
              }
              str = (CharPtr) MemNew ((size_t) (len + 1) * sizeof (Char));
              protein_seq = str;
              if (ajp->flags.iupacaaOnly) {
                code = Seq_code_iupacaa;
              } else {
                code = Seq_code_ncbieaa;
              }
              spp = SeqPortNewByLoc (sfp->product, code);
              if (spp != NULL) {
                SeqPortSet_do_virtual (spp, TRUE);
                while ((residue = SeqPortGetResidue (spp)) != SEQPORT_EOF) {
                  if (! (IS_residue (residue))) continue;
                  if (residue == INVALID_RESIDUE) {
                    residue = (Uint1) 'X';
                  }
                  *protein_seq = residue;
                  protein_seq++;
                }
                /*
                if (at_end && StringLen (str) < GENPEPT_MIN) {
                  str = MemFree (str);
                }
                */
                if (! StringHasNoText (str)) {
                  FFAddTextToString(ffstring, "/translation=\"", str, "\"", 
                                    FALSE, TRUE, TILDE_TO_SPACES);
                  FFAddOneChar(ffstring, '\n', FALSE);
                }
                MemFree (str);
              } else {
                ajp->relModeError = TRUE;
              }
              SeqPortFree (spp);
            } else {
              ajp->relModeError = TRUE;
            }
          }
        }
        break;
      
      case Qual_class_illegal :
        for (vnp = qvp [idx].vnp; vnp != NULL; vnp = vnp->next) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (str != NULL) {
            FFAddOneString(ffstring, str, FALSE, FALSE, TILDE_TO_SPACES);
            FFAddNewLine(ffstring);
          }
        }
        break;

      case Qual_class_note :
        /*head = NULL;*/
        notestr = NULL;
        prefix = NULL;
        add_period = FALSE;
        suppress_period = FALSE;
        lasttype = NULL;
        isTRNA = TRUE;
        

#ifdef DISPLAY_STRINGS
        s_DisplayQVP (qvp, feat_note_order);
#endif
        for (j = 0, jdx = feat_note_order [j]; jdx != 0; j++, jdx = feat_note_order [j]) {
          switch (asn2gnbk_featur_quals [jdx].qualclass) {

            case Qual_class_string :
              if (! StringHasNoText (qvp [jdx].str)) {
                if (jdx == FTQUAL_figure) {
                  if (!IsEllipsis (qvp [jdx].str))
                    s_RemovePeriodFromEnd (qvp [jdx].str);
                  sprintf (buf, "This sequence comes from %s", qvp [jdx].str);
                  FFAddString_NoRedund (unique, prefix, buf, NULL);
                  add_period = FALSE;
                } else if (jdx == FTQUAL_maploc) {
                  if (!IsEllipsis (qvp [jdx].str))
                    s_RemovePeriodFromEnd (qvp [jdx].str);
                  sprintf (buf, "Map location %s", qvp [jdx].str);
                  FFAddString_NoRedund (unique, prefix, buf, NULL);
                  add_period = FALSE;
                } else if (jdx == FTQUAL_seqfeat_note) {
                  str = StringSave (qvp [jdx].str);
                  TrimSpacesAndJunkFromEnds (str, TRUE);
                  if (! IsEllipsis (str))
                    add_period = s_RemovePeriodFromEnd (str);
                  /*  NOTE -- The following function call cleans up some strings
                      (i.e., U34661 & U31565) but should be commented back
                      in only if the problem can't be fixed upstream of here

                  s_StringCleanup(str);

                  */
                  FFAddString_NoRedund (unique, prefix, str, NULL);
                  MemFree (str);
                  if (hadProtDesc) {
                    suppress_period = TRUE;
                  }
                } else if (jdx == FTQUAL_prot_note) {
                  str = StringSave (qvp [jdx].str);
                  TrimSpacesAndJunkFromEnds (str, TRUE);
                  if (! IsEllipsis (str))
                    s_RemovePeriodFromEnd (str);
                  FFAddString_NoRedund (unique, prefix, str, NULL);
                  MemFree (str);
                  add_period = FALSE;
                } else if (jdx == FTQUAL_prot_desc) {
                  str = StringSave (qvp [jdx].str);
                  TrimSpacesAndJunkFromEnds (str, TRUE);
                  if (! IsEllipsis (str))
                    add_period = s_RemovePeriodFromEnd (str);
                  FFAddString_NoRedund (unique, prefix, str, NULL);
                  MemFree (str);
                  hadProtDesc = TRUE; /* gi|347886|gb|M96268.1|ECOUBIA */
                } else {
                  if (! IsEllipsis (qvp [jdx].str)) {
                    s_RemovePeriodFromEnd (qvp [jdx].str);
                  }
                  FFAddString_NoRedund (unique, prefix, qvp [jdx].str, NULL);
                  add_period = FALSE;
                }
                prefix = "; ";
              }
              break;

            case Qual_class_method :
              if (! StringHasNoText (qvp [jdx].str)) {
                if ( FFEmpty(unique) ) {
                  prefix = "Method: ";
                } else {
                  prefix = "; Method: ";
                }
                FFAddString_NoRedund (unique, prefix, qvp [jdx].str, NULL);
                prefix = "; ";
                add_period = FALSE;
              }
              break;

            case Qual_class_valnode :
              for (vnp = qvp [jdx].vnp; vnp != NULL; vnp = vnp->next) {
                str = (CharPtr) vnp->data.ptrvalue;
                if (! StringHasNoText (str)) {
                  FFAddString_NoRedund (unique, prefix, str, NULL);
                  prefix = "; ";
                  add_period = FALSE;
                }
              }
              break;

            case Qual_class_gene_syn :
              numsyns = 0;
              for (vnp = qvp [jdx].vnp; vnp != NULL; vnp = vnp->next) {
                str = (CharPtr) vnp->data.ptrvalue;
                if (! StringHasNoText (str)) {
                  numsyns++;
                }
              }
              if (numsyns > 0) {
                if (numsyns > 1) {
                  FFAddTextToString (unique, prefix, "synonyms: ", NULL, FALSE, FALSE, TILDE_IGNORE);
                } else {
                  FFAddTextToString (unique, prefix, "synonym: ", NULL, FALSE, FALSE, TILDE_IGNORE);
                }
                prefix = NULL;
                for (vnp = qvp [jdx].vnp; vnp != NULL; vnp = vnp->next) {
                  str = (CharPtr) vnp->data.ptrvalue;
                  if (! StringHasNoText (str)) {
                    FFAddTextToString (unique, prefix, str, NULL, FALSE, FALSE, TILDE_IGNORE);
                    prefix = ", ";
                  }
                }
                prefix = "; ";
                add_period = FALSE;
              }
              break;

            case Qual_class_region :
              if (! StringHasNoText (qvp [jdx].str)) {
                if ( FFEmpty(unique) ) {
                  prefix = "Region: ";
                } else {
                  prefix = "; Region: ";
                }
#ifdef ASN2GNBK_STRIP_NOTE_PERIODS
                FFAddTextToString(unique, prefix, qvp [jdx].str, NULL, FALSE, FALSE, TILDE_IGNORE);
#else
                FFAddString_NoRedund (unique, prefix, qvp [jdx].str, NULL);
#endif
                prefix = "; ";
                add_period = FALSE;
              }
              break;

            case Qual_class_site :
              if (! StringHasNoText (qvp [jdx].str)) {
                FFAddString_NoRedund (unique, prefix, qvp [jdx].str, " site");
                add_period = FALSE;
                prefix = "\n";
              }
              break;

            case Qual_class_bond :
              if (! StringHasNoText (qvp [jdx].str)) {
                FFAddString_NoRedund (unique, prefix, qvp [jdx].str, " bond");
                add_period = FALSE;
                prefix = "\n";
              }
              break;

            case Qual_class_protnames :
              /* process gene sgml for check against subsequent protein names */
              start = NULL;
              if (! StringHasNoText (qvp [FTQUAL_gene].str)) {
                if (is_journalscan) {
                  ascii_len = Sgml2AsciiLen (qvp [FTQUAL_gene].str);
                  start = ascii = MemNew ((size_t) (10 + ascii_len));
                  if (start != NULL) {
                    ascii = Sgml2Ascii (qvp [FTQUAL_gene].str, ascii, ascii_len + 1);
                  }
                } else {
                  start = StringSaveNoNull (qvp [FTQUAL_gene].str);
                }
              }
              for (vnp = qvp [jdx].vnp; vnp != NULL; vnp = vnp->next) {
                str = (CharPtr) vnp->data.ptrvalue;
                if (! StringHasNoText (str)) {
                  /* case sensitive - gi|4973426|gb|AF148501.1|AF148501 */
                  /* check with and without sgml conversion */
                  if (StringCmp (start, str) != 0 &&
                      StringCmp (qvp [FTQUAL_gene].str, str) != 0) {
                    if (! StringStr (qvp [FTQUAL_prot_desc].str, str)) {
                      if (NotInGeneSyn (str, gene_syn)) {
                        FFAddString_NoRedund (unique, prefix, str, NULL);
                        prefix = "; ";
                        add_period = FALSE;
                      }
                    }
                  }
                }
              }
              MemFree (start);
              break;

            case Qual_class_xtraprds :
              gbq = qvp [jdx].gbq;
              if (lasttype == NULL && gbq != NULL) {
                lasttype = gbq->qual;
              }
              while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
                if (! StringHasNoText (gbq->val)) {
                  if (StringCmp(gbq->val,qvp[FTQUAL_gene].str) != 0 &&
                      StringCmp(gbq->val,qvp[FTQUAL_product].str) != 0) {
                    if (!isTRNA || !StringStr (gbq->val, "RNA")) {
                      FFAddString_NoRedund (unique, prefix, gbq->val, NULL);
                      prefix = "; ";
                      add_period = FALSE;
                    }
                  }
                }
                gbq = gbq->next;
              }
              break;

            case Qual_class_its :
              str = qvp [jdx].str;
              if (! StringHasNoText (str)) {
                if (sfp->comment == NULL || StringStr (sfp->comment, str) == NULL) {
                  FFAddString_NoRedund (unique, prefix, str, NULL);
                  prefix = "; ";
                  add_period = FALSE;
                }
              }
              break;

            case Qual_class_trna_codons :
              trna = qvp [jdx].trp;
              if (trna) {
                numcodons = ComposeCodonsRecognizedString (trna, numbuf, sizeof (numbuf));
                if (numcodons < 1 || StringHasNoText (numbuf)) {
                } else if (numcodons == 1) {
                  isTRNA = TRUE;
                  sprintf (buf, "codon recognized: %s", numbuf);
                  if (StringStr (qvp [FTQUAL_seqfeat_note].str, buf) == NULL) {
                    FFAddString_NoRedund (unique, prefix, "codon recognized: ", numbuf);
                    prefix = "; ";
                  }
                } else {
                  isTRNA = TRUE;
                  FFAddString_NoRedund (unique, prefix, "codons recognized: ", numbuf);
                  prefix = "; ";
                  add_period = FALSE;
                }
              }
              break;

            case Qual_class_model_ev :
              uop = qvp [jdx].uop;
              if (uop != NULL) {
                str = NULL;
                VisitUserObjectsInUop (sfp->ext, (Pointer) &str, GetStrFormRNAEvidence);
                if (! StringHasNoText (str)) {
                  FFAddString_NoRedund (unique, prefix, str, NULL);
                  prefix = "; ";
                  add_period = FALSE;
                }
              }
              break;

            default :
              break;
          }
        }
        if ( !FFEmpty(unique) ) {
          notestr = FFToCharPtr(unique);
          TrimSpacesAroundString (notestr);
          if (add_period) {
            if (! suppress_period) {
              s_AddPeriodToEnd (notestr);
            }
          }

#ifdef ASN2GNBK_STRIP_NOTE_PERIODS
          if (! IsEllipsis (notestr))
            s_RemovePeriodFromEnd (notestr);
#endif
    
          FFAddOneString(ffstring, "/note=\"", FALSE, FALSE, TILDE_IGNORE);
          FFAddOneString(ffstring, notestr, FALSE, TRUE, TILDE_EXPAND);
          FFAddOneString(ffstring, "\"\n", FALSE, FALSE, TILDE_IGNORE);

          MemFree (notestr);
          /*ValNodeFreeData (head);*/
        }
        break;

      default:
        break;

    }
  }
  FFRecycleString(ajp, unique);
}



static void FF_asn2gb_www_featkey (
  StringItemPtr ffstring,
  CharPtr key,
  SeqLocPtr slp,
  Int4 from,
  Int4 to,
  Uint1 strand,
  Uint4 itemID
)

{
  BioseqPtr  bsp;
  Int4       gi = 0;
  SeqIdPtr   sip;
  Boolean    is_aa = FALSE;
  Char gi_buf[16];
  Char itemID_buf[16];

  bsp = BioseqFindFromSeqLoc (slp);
  if (bsp != NULL) {
    is_aa = ISA_aa (bsp->mol);
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_GI) {
        gi = (Int4) sip->data.intvalue;
      }
    }
  }

  sprintf(gi_buf, "%ld", (long)gi);
  sprintf(itemID_buf, "%ld", (long)itemID);


  FFAddOneString(ffstring, "<a href=", FALSE, FALSE, TILDE_IGNORE);
  FFAddOneString(ffstring, link_feat, FALSE, FALSE, TILDE_IGNORE);
  FFAddOneString(ffstring, "val=", FALSE, FALSE, TILDE_IGNORE);
  FFAddOneString(ffstring, gi_buf, FALSE, FALSE, TILDE_IGNORE);
  FFAddOneString(ffstring, "&itemID=", FALSE, FALSE, TILDE_IGNORE);
  FFAddOneString(ffstring, itemID_buf, FALSE, FALSE, TILDE_IGNORE);
  

  if ( is_aa ) {
    FFAddOneString(ffstring, "&view=gpwithparts>", FALSE, FALSE, TILDE_IGNORE);
  } else {
    FFAddOneString(ffstring, "&view=gbwithparts>", FALSE, FALSE, TILDE_IGNORE);
  }

  FFAddOneString(ffstring, key, FALSE, FALSE, TILDE_IGNORE);
  FFAddOneString(ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
}



static CharPtr FormatFeatureBlock (
  Asn2gbFormatPtr afp,
  BaseBlockPtr bbp
)

{
  Uint1              aa;
  IntAsn2gbJobPtr    ajp;
  Asn2gbSectPtr      asp;
  Int2               bondidx;
  BioseqPtr          bsp;
  BioseqSetPtr       bssp;
  Char               buf [80];
  Choice             cbaa;
  CodeBreakPtr       cbp;
  BioseqPtr          cdna;
  SeqFeatPtr         cds;
  Char               ch;
  Uint1              code = Seq_code_ncbieaa;
  CdRegionPtr        crp;
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext;
  Uint1              featdeftype;
  Uint1              from;
  GBQualPtr          gbq;
  GBFeaturePtr       gbfeat = NULL;
  GBSeqPtr           gbseq;
  SeqMgrFeatContext  gcontext;
  ValNodePtr         gcp;
  SeqFeatPtr         gene = NULL;
  ValNodePtr         gene_syn = NULL;
  GeneRefPtr         grp;
  IntCdsBlockPtr     icp;
  Uint2              idx;
  IntFeatBlockPtr    ifp;
  ValNodePtr         illegal = NULL;
  ImpFeatPtr         imp = NULL;
  IndxPtr            index;
  Boolean            is_gps = FALSE;
  Boolean            is_journalscan = FALSE;
  Boolean            is_other = FALSE;
  CharPtr            key;
  CharPtr            lasttype = NULL;
  Int4               left = -1;
  SeqLocPtr          loc = NULL;
  SeqLocPtr          location = NULL;
  SeqLocPtr          locforgene = NULL;
  SeqMgrFeatContext  mcontext;
  MolInfoPtr         mip;
  SeqFeatPtr         mrna;
  SeqLocPtr          newloc;
  Boolean            noLeft;
  Boolean            noRight;
  SeqEntryPtr        oldscope;
  Uint2              partial;
  SeqMgrFeatContext  pcontext;
  BioseqPtr          prd;
  BioseqPtr          prod = NULL;
  SeqFeatPtr         prot;
  Boolean            protein = FALSE;
  Char               protein_pid_g [32];
  ProtRefPtr         prp;
  ProtRefPtr         prpxref;
  Boolean            pseudo = FALSE;
  CharPtr            ptr;
  Int2               qualclass;
  QualValPtr         qvp;
  Uint1              residue;
  Int4               right = -1;
  RnaRefPtr          rrp;
  SeqCodeTablePtr    sctp;
  SeqDescrPtr        sdp;
  SeqEntryPtr        sep;
  Uint1              seqcode;
  SeqFeatPtr         sfp;
  Uint1              shift;
  SeqIdPtr           sip;
  Int2               siteidx;
  SeqMapTablePtr     smtp;
  Boolean            split;
  CharPtr            str;
  Uint1              strand = Seq_strand_unknown;
  BioseqPtr          target;
  CharPtr            tmp;
  tRNAPtr            trna;
  BioseqPtr          unlockme = NULL;
  ValNodePtr         vnp;
  StringItemPtr      ffstring;


  if (afp == NULL || bbp == NULL) return NULL;
  asp = afp->asp;
  if (asp == NULL) return NULL;
  target = asp->target;
  bsp = asp->bsp;
  if (target == NULL || bsp == NULL) return NULL;
  ajp = afp->ajp;
  if (ajp == NULL) return NULL;
  qvp = afp->qvp;
  if (qvp == NULL) return NULL;

  ffstring = FFGetString(ajp);
  if ( ffstring == NULL ) return NULL;

  if (ajp->index) {
    index = &asp->index;
  } else {
    index = NULL;
  }

  if (ajp->gbseq) {
    gbseq = &asp->gbseq;
  } else {
    gbseq = NULL;
  }
  
  protein_pid_g [0] = '\0';

  /* all features in this list are known to be valid for the designated mode */

  sfp = SeqMgrGetDesiredFeature (bbp->entityID, NULL, bbp->itemID, 0, NULL, &fcontext);
  if (sfp == NULL) return NULL;

  /* may need to map location between aa and dna */

  ifp = (IntFeatBlockPtr) bbp;
  if (ifp->mapToNuc) {

    /* map mat_peptide, etc., to nucleotide coordinates */

    sip = SeqLocId (sfp->location);
    prd = BioseqFind (sip);
    cds = SeqMgrGetCDSgivenProduct (prd, NULL);
    CheckSeqLocForPartial (sfp->location, &noLeft, &noRight);
    location = aaFeatLoc_to_dnaFeatLoc (cds, sfp->location);
    SetSeqLocPartial (location, noLeft, noRight);
    locforgene = location;
    loc = location;

  } else if (ifp->mapToProt) {

    /* map CDS to protein product coordinates */

    sip = SeqLocIdForProduct (sfp->product);
    prd = BioseqFind (sip);
    cds = SeqMgrGetCDSgivenProduct (prd, NULL);
    location = dnaLoc_to_aaLoc (cds, sfp->location, TRUE, NULL, FALSE);
    SetSeqLocPartial (location, FALSE, FALSE);
    locforgene = sfp->location;
    loc = location;

  } else if (ifp->mapToGen) {

    /* map CDS from cDNA to genomic Bioseq */

    cdna = BioseqFindFromSeqLoc (sfp->location);
    mrna = SeqMgrGetRNAgivenProduct (cdna, &mcontext);
    CheckSeqLocForPartial (sfp->location, &noLeft, &noRight);
    location = productLoc_to_locationLoc (mrna, sfp->location);
    SetSeqLocPartial (location, noLeft, noRight);
    locforgene = location;
    loc = location;

  } else if (ifp->mapToMrna) {

    /* map gene from genomic to cDNA Bioseq */

    sep = SeqMgrGetSeqEntryForData (bsp);
    location = CreateWholeInterval (sep);
    SetSeqLocPartial (location, FALSE, FALSE);
    locforgene = location;
    loc = location;

  } else if (ifp->mapToPep) {

    /* map protein processing from precursor to subpeptide Bioseq */

    sep = SeqMgrGetSeqEntryForData (bsp);
    location = CreateWholeInterval (sep);
    SetSeqLocPartial (location, FALSE, FALSE);
    locforgene = location;
    loc = location;

  } else {

    /* no aa-dna or dna-aa mapping, just use location */

    location = sfp->location;
    locforgene = sfp->location;
  }
  if (location == NULL) return NULL;

  sep = GetTopSeqEntryForEntityID (ajp->ajp.entityID);
  if (sep == NULL && IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && bssp->_class == BioseqseqSet_class_gen_prod_set) {
      is_gps = TRUE;
    }
  }

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_OTHER) {
      is_other = TRUE;
    } else if (sip->choice == SEQID_GIBBSQ ||
               sip->choice == SEQID_GIBBMT ||
               sip->choice == SEQID_GIIM) {
      is_journalscan = TRUE;
    }
  }

  featdeftype = fcontext.featdeftype;
  if (featdeftype < FEATDEF_GENE || featdeftype >= FEATDEF_MAX) {
    featdeftype = FEATDEF_BAD;
  }
  key = FindKeyFromFeatDefType (featdeftype, TRUE);

  if (afp->format == GENPEPT_FMT && ISA_aa (bsp->mol)) {
    if (featdeftype == FEATDEF_REGION) {
      key = "Region";
    } else if (featdeftype == FEATDEF_BOND) {
      key = "Bond";
    } else if (featdeftype == FEATDEF_SITE) {
      key = "Site";
    } else if (featdeftype == FEATDEF_preprotein) {
      key = "Proprotein";
    }
    if (ifp->mapToPep) {
      if (featdeftype >= FEATDEF_preprotein && featdeftype <= FEATDEF_transit_peptide_aa) {
        key = "Precursor";
      }
    }
  }

  /* deal with unmappable impfeats */

  if (featdeftype == FEATDEF_BAD && fcontext.seqfeattype == SEQFEAT_IMP) {
    imp = (ImpFeatPtr) sfp->data.value.ptrvalue;
    if (imp != NULL) {
      key = imp->key;
    }
  }

  FFStartPrint(ffstring, afp->format, 5, 21, NULL, 0, 5, 21, "FT", ifp->firstfeat);
  if (ajp->ajp.slp != NULL) {
    FFAddOneString(ffstring, key, FALSE, FALSE, TILDE_IGNORE);
  } else if ( GetWWW(ajp) /* && SeqMgrGetParentOfPart (bsp, NULL) == NULL */ ) {
    FF_asn2gb_www_featkey (ffstring, key, sfp->location, fcontext.left + 1, fcontext.right + 1, fcontext.strand, fcontext.itemID);
  } else {
    FFAddOneString(ffstring, key, FALSE, FALSE, TILDE_IGNORE);
  }
  FFAddNChar(ffstring, ' ', 21 - 5 - StringLen(key), FALSE);

  if (gbseq != NULL) {
    gbfeat = GBFeatureNew ();
    if (gbfeat != NULL) {
      gbfeat->key = StringSave (key);
    }
  }

  if (imp == NULL || StringHasNoText (imp->loc)) {

    
    if (ajp->ajp.slp != NULL) {
      sip = SeqIdParse ("lcl|dummy");
      left = GetOffsetInBioseq (ajp->ajp.slp, bsp, SEQLOC_LEFT_END);
      right = GetOffsetInBioseq (ajp->ajp.slp, bsp, SEQLOC_RIGHT_END);
      strand = SeqLocStrand (ajp->ajp.slp);
      split = FALSE;
      newloc = SeqLocReMap (sip, ajp->ajp.slp, location, 0, FALSE);
      /*
      newloc = SeqLocCopyRegion (sip, location, bsp, left, right, strand, &split);
      */
      SeqIdFree (sip);
      if (newloc == NULL) return NULL;
      A2GBSeqLocReplaceID (newloc, ajp->ajp.slp);
      str = FlatLoc (ajp, target, newloc, ajp->masterStyle);
      SeqLocFree (newloc);
    } else {
      str = FlatLoc (ajp, target, location, ajp->masterStyle);
    }
  } else {
    str = StringSave (imp->loc);
  }
  if ( GetWWW(ajp) ) {
    FF_www_featloc (ffstring, str);
  } else {
    FFAddOneString(ffstring, str, FALSE, FALSE, TILDE_IGNORE);
  }

  if (gbseq != NULL) {
    if (gbfeat != NULL) {
      gbfeat->location = StringSave (str);
      if (ajp->masterStyle) {
        AddIntervalsToGbfeat (gbfeat, location, target);
      } else {
        AddIntervalsToGbfeat (gbfeat, location, NULL);
      }
    }
  }

  MemFree (str);

  /* populate qualifier table from feature fields */

  /*
  if (sfp->partial == TRUE)
    sfp->partial = FlatAnnotPartial(sfp, use_product);
  */

  if (sfp->partial) {
    partial = SeqLocPartialCheck (location);
    if (partial == SLP_COMPLETE /* || partial > SLP_OTHER */ ) {
      qvp [FTQUAL_partial].ble = TRUE;
    }
    if (LookForFuzz (location)) {
      qvp [FTQUAL_partial].ble = FALSE;
    }
    if (imp != NULL) {
      if (StringChr (imp->loc, '<') != NULL || StringChr (imp->loc, '>') != NULL) {
        qvp [FTQUAL_partial].ble = FALSE;
      }
    }

    /* a few features cannot show /partial in RELEASE_MODE - later no features will */

    if (ajp->flags.checkQualSyntax) {
      switch (featdeftype) {
      case FEATDEF_conflict:
      case FEATDEF_mutation:
      case FEATDEF_N_region:
      case FEATDEF_polyA_site:
        qvp [FTQUAL_partial].ble = FALSE;
        break;
      default:
        break;
      }
    }
  }
  if (ifp->mapToProt) {
    qvp [FTQUAL_partial].ble = FALSE;
  }

  if (sfp->pseudo) {
    pseudo = TRUE;
  }

  if (fcontext.seqfeattype == SEQFEAT_GENE) {
    grp = (GeneRefPtr) sfp->data.value.ptrvalue;
    if (grp != NULL) {
      if (! StringHasNoText (grp->locus)) {
        qvp [FTQUAL_gene].str = grp->locus;
        qvp [FTQUAL_gene_desc].str = grp->desc;
        qvp [FTQUAL_gene_syn].vnp = grp->syn;
        qvp [FTQUAL_locus_tag].str = grp->locus_tag;
      } else if (! StringHasNoText (grp->desc)) {
        qvp [FTQUAL_gene].str = grp->desc;
        qvp [FTQUAL_gene_syn].vnp = grp->syn;
        qvp [FTQUAL_locus_tag].str = grp->locus_tag;
      } else if (grp->syn != NULL) {
        vnp = grp->syn;
        qvp [FTQUAL_gene].str = (CharPtr) vnp->data.ptrvalue;
        vnp = vnp->next;
        qvp [FTQUAL_gene_syn].vnp = vnp;
        qvp [FTQUAL_locus_tag].str = grp->locus_tag;
      } else if (grp->locus_tag != NULL) {
        /* - for now it can be used as /gene */
        qvp [FTQUAL_gene].str = grp->locus_tag;
      }
      qvp [FTQUAL_gene_map].str = grp->maploc;
      qvp [FTQUAL_gene_allele].str = grp->allele;
      /* - for later - right now it can become the /gene
      qvp [FTQUAL_locus_tag].str = grp->locus_tag;
      */
      qvp [FTQUAL_gene_xref].vnp = grp->db;
      if (grp->pseudo) {
        pseudo = TRUE;
      }
    }

  } else {

    grp = SeqMgrGetGeneXref (sfp);
    if (grp != NULL) {
      qvp [FTQUAL_gene_xref].vnp = grp->db;
    }
    if (grp == NULL) {
      sep = GetTopSeqEntryForEntityID (ajp->ajp.entityID);
      oldscope = SeqEntrySetScope (sep);
      gene = SeqMgrGetOverlappingGene (locforgene, &gcontext);
      SeqEntrySetScope (oldscope);
      if (gene != NULL) {
        qvp [FTQUAL_gene_note].str = gene->comment;
        grp = (GeneRefPtr) gene->data.value.ptrvalue;
        if (gene->pseudo) {
          pseudo = TRUE;
        }
        if (grp != NULL && grp->db != NULL) {
          qvp [FTQUAL_gene_xref].vnp = grp->db;
        } else {
          qvp [FTQUAL_gene_xref].vnp = gene->dbxref;
        }
      }
    }
    if (grp != NULL && grp->pseudo) {
      pseudo = TRUE;
    }
    if (grp != NULL && (! SeqMgrGeneIsSuppressed (grp)) &&
        (fcontext.featdeftype != FEATDEF_repeat_region || gene == NULL)) {
      if (! StringHasNoText (grp->locus)) {
        qvp [FTQUAL_gene].str = grp->locus;
        gene_syn = grp->syn;
      } else if (! StringHasNoText (grp->desc)) {
        qvp [FTQUAL_gene].str = grp->desc;
        gene_syn = grp->syn;
      } else if (grp->syn != NULL) {
        vnp = grp->syn;
        qvp [FTQUAL_gene].str = (CharPtr) vnp->data.ptrvalue;
        vnp = vnp->next;
        gene_syn = vnp;
      } else if (! StringHasNoText (grp->locus_tag)) {
        qvp [FTQUAL_gene].str = grp->locus_tag;
      }
    }
    if (fcontext.seqfeattype != SEQFEAT_CDREGION &&
        fcontext.seqfeattype != SEQFEAT_RNA) {
      qvp [FTQUAL_gene_xref].vnp = NULL;
    }

    /* specific fields set here */

    switch (fcontext.seqfeattype) {
      case SEQFEAT_CDREGION :
        if (! ifp->mapToProt) {
          crp = (CdRegionPtr) sfp->data.value.ptrvalue;
          if (crp != NULL) {

            qvp [FTQUAL_codon_start].num = crp->frame;
            if (qvp [FTQUAL_codon_start].num == 0) {
              qvp [FTQUAL_codon_start].num = 1;
            }
            qvp [FTQUAL_transl_except].cbp = crp->code_break;
            for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next) {
              seqcode = 0;
              sctp = NULL;
              cbaa = cbp->aa;
              switch (cbaa.choice) {
                case 1 :
                  seqcode = Seq_code_ncbieaa;
                  break;
                case 2 :
                  seqcode = Seq_code_ncbi8aa;
                  break;
                case 3 :
                  seqcode = Seq_code_ncbistdaa;
                  break;
                default :
                  break;
              }
              if (seqcode != 0) {
                sctp = SeqCodeTableFind (seqcode);
                if (sctp != NULL) {
                  residue = cbaa.value.intvalue;
                  if (residue != 42) {
                    if (seqcode != Seq_code_ncbieaa) {
                      smtp = SeqMapTableFind (seqcode, Seq_code_ncbieaa);
                      residue = SeqMapTableConvert (smtp, residue);
                    }
                    if (residue == 'U') {
                      qvp [FTQUAL_selenocysteine].str = "selenocysteine";
                    }
                  }
                }
              }
            }

            gcp = crp->genetic_code;
            if (gcp != NULL) {
              for (vnp = gcp->data.ptrvalue; vnp != NULL; vnp = vnp->next) {
                if (vnp->choice == 2 && vnp->data.intvalue != 0) {
                  qvp [FTQUAL_transl_table].num = vnp->data.intvalue;
                }
              }

              /* suppress table 1 */

              if (qvp [FTQUAL_transl_table].num == 1) {
                qvp [FTQUAL_transl_table].num = 0;
              }
            }

            if (sfp->product != NULL && SeqLocLen (sfp->product) != 0) {
              protein = TRUE;
            }
            if (crp->conflict && (protein || (! sfp->excpt))) {
              if (protein) {
                qvp [FTQUAL_prot_conflict].str = conflict_msg;
              } else {
                /*
                qvp [FTQUAL_prot_missing].str = no_protein_msg;
                */
              }
            }
          }

          sip = SeqLocIdForProduct (sfp->product);
          qvp [FTQUAL_protein_id].sip = sip;

          sep = GetTopSeqEntryForEntityID (ajp->ajp.entityID);

          if (! ajp->alwaysTranslCds) {

            /* by default only show /translation if product bioseq is within entity */

            oldscope = SeqEntrySetScope (sep);
            prod = BioseqFind (sip);
            SeqEntrySetScope (oldscope);

            if (prod == NULL && ajp->showFarTransl) {

              /* but flag can override and force far /translation */

              prod = BioseqLockById (sip);
              unlockme = prod;
            }
          }

          prp = NULL;

          if (prod != NULL) {
            for (sip = prod->id; sip != NULL; sip = sip->next) {
              if (sip->choice == SEQID_GI) {
                sprintf (protein_pid_g, "PID:g%ld", (long) sip->data.intvalue);
              }
            }
            sdp = SeqMgrGetNextDescriptor (prod, NULL, Seq_descr_comment, &dcontext);
            if (sdp != NULL && dcontext.level == 0) {
              if (! StringHasNoText ((CharPtr) sdp->data.ptrvalue)) {
                qvp [FTQUAL_prot_comment].str = (CharPtr) sdp->data.ptrvalue;
              }
            }
            sdp = SeqMgrGetNextDescriptor (prod, NULL, Seq_descr_molinfo, &dcontext);
            if (sdp != NULL && dcontext.level == 0) {
              mip = (MolInfoPtr) sdp->data.ptrvalue;
              if (mip != NULL && mip->tech > 1 &&
                  mip->tech != MI_TECH_concept_trans &&
                  mip->tech != MI_TECH_concept_trans_a) {
                str = StringForSeqTech (mip->tech);
                if (! StringHasNoText (str)) {
                  qvp [FTQUAL_prot_method].str = str;
                }
              }
            }
            prot = SeqMgrGetBestProteinFeature (prod, &pcontext);
            if (prot != NULL) {
              prp = (ProtRefPtr) prot->data.value.ptrvalue;
              if (prp != NULL && prp->processed < 2) {
                qvp [FTQUAL_prot_note].str = prot->comment;
              }
            }
          }

          /* protein xref overrides names, but should not prevent /protein_id, etc. */

          prpxref = SeqMgrGetProtXref (sfp);
          if (prpxref != NULL) {
            prp = prpxref;
          }
          if (prp != NULL) {
            vnp = prp->name;
            if (vnp != NULL && (! StringHasNoText ((CharPtr) vnp->data.ptrvalue))) {
              qvp [FTQUAL_cds_product].str = (CharPtr) vnp->data.ptrvalue;
              vnp = vnp->next;
              qvp [FTQUAL_prot_names].vnp = vnp;
            }
            qvp [FTQUAL_prot_desc].str = prp->desc;
            qvp [FTQUAL_prot_activity].vnp = prp->activity;
            qvp [FTQUAL_prot_EC_number].vnp = prp->ec;
          }

          if (! pseudo) {
            if (prod != NULL || ajp->transIfNoProd || ajp->alwaysTranslCds) {
              qvp [FTQUAL_translation].ble = TRUE;
            }
          }

          if (ifp->isCDS) {
            icp = (IntCdsBlockPtr) ifp;
            qvp [FTQUAL_figure].str = icp->fig;
            qvp [FTQUAL_maploc].str = icp->maploc;
          }
        } else {
          qvp [FTQUAL_coded_by].slp = sfp->location;

          crp = (CdRegionPtr) sfp->data.value.ptrvalue;
          if (crp != NULL) {
            gcp = crp->genetic_code;
            if (gcp != NULL) {
              for (vnp = gcp->data.ptrvalue; vnp != NULL; vnp = vnp->next) {
                if (vnp->choice == 2 && vnp->data.intvalue != 0) {
                  qvp [FTQUAL_transl_table].num = vnp->data.intvalue;
                }
              }

              /* suppress table 1 */

              if (qvp [FTQUAL_transl_table].num == 1) {
                qvp [FTQUAL_transl_table].num = 0;
              }
            }
          }
        }
        break;
      case SEQFEAT_PROT :
        if (! ifp->mapToPep) {
          prp = (ProtRefPtr) sfp->data.value.ptrvalue;
          if (prp != NULL) {
            vnp = prp->name;
            if (vnp != NULL && (! StringHasNoText ((CharPtr) vnp->data.ptrvalue))) {
              qvp [FTQUAL_product].str = (CharPtr) vnp->data.ptrvalue;
              vnp = vnp->next;
              qvp [FTQUAL_prot_names].vnp = vnp;
            }
            if (afp->format != GENPEPT_FMT) {
              qvp [FTQUAL_prot_desc].str = prp->desc;
            } else {
              qvp [FTQUAL_prot_name].str = prp->desc;
            }
            if (afp->format != GENPEPT_FMT || prp->processed != 2) {
              qvp [FTQUAL_prot_activity].vnp = prp->activity;
            }
            qvp [FTQUAL_prot_EC_number].vnp = prp->ec;
          }
          sip = SeqLocIdForProduct (sfp->product);
          if (sip != NULL) {
            /* for RefSeq records or GenBank not release_mode */
            if (is_other || (! ajp->flags.forGbRelease)) {
              qvp [FTQUAL_protein_id].sip = sip;
            }
            prod = BioseqFind (sip);
          }
        } else {
          qvp [FTQUAL_coded_by].slp = sfp->location;
        }
        break;
      case SEQFEAT_RNA :
        rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
        if (rrp != NULL) {
          if (rrp->pseudo) {
            pseudo = TRUE;
          }
          if (rrp->type == 2) {
            sip = SeqLocIdForProduct (sfp->product);
            if (sip != NULL) {
              /* for RefSeq records or GenBank not release_mode */
              if (is_other || (! ajp->flags.forGbRelease)) {
                qvp [FTQUAL_transcript_id].sip = sip;
              }
              prod = BioseqFind (sip);
            }
          }
          if (rrp->type == 3) {
            if (rrp->ext.choice == 1) {
              /* amino acid could not be parsed into structured form */
              if (! ajp->flags.dropIllegalQuals) {
                str = (CharPtr) rrp->ext.value.ptrvalue;
                qvp [FTQUAL_product].str = str;
              } else {
                qvp [FTQUAL_product].str = "tRNA-OTHER";
              }
            } else if (rrp->ext.choice == 2) {
              trna = (tRNAPtr) rrp->ext.value.ptrvalue;
              if (trna != NULL) {
                aa = 0;
                if (trna->aatype == 2) {
                  aa = trna->aa;
                } else {
                  from = 0;
                  switch (trna->aatype) {
                    case 0 :
                      from = 0;
                      break;
                    case 1 :
                      from = Seq_code_iupacaa;
                      break;
                    case 2 :
                      from = Seq_code_ncbieaa;
                      break;
                    case 3 :
                      from = Seq_code_ncbi8aa;
                      break;
                    case 4 :
                      from = Seq_code_ncbistdaa;
                      break;
                    default:
                      break;
                  }
                  if (ajp->flags.iupacaaOnly) {
                    code = Seq_code_iupacaa;
                  } else {
                    code = Seq_code_ncbieaa;
                  }
                  smtp = SeqMapTableFind (code, from);
                  if (smtp != NULL) {
                    aa = SeqMapTableConvert (smtp, trna->aa);
                    if (aa == 255 && from == Seq_code_iupacaa && trna->aa == 'U') {
                      aa = 'U';
                    }
                  }
                }
                if (ajp->flags.iupacaaOnly) {
                  smtp = SeqMapTableFind (Seq_code_iupacaa, Seq_code_ncbieaa);
                  if (smtp != NULL) {
                    aa = SeqMapTableConvert (smtp, trna->aa);
                  } else {
                    aa = 'X';
                  }
                }
                if (aa > 0 && aa != 255) {
                  if (aa <= 74) {
                    shift = 0;
                  } else if (aa > 79) {
                    shift = 2;
                  } else {
                    shift = 1;
                  }
                  idx = aa - (64 + shift);
                  if (idx > 0 && idx < 25) {
                    str = trnaList [idx];
                    qvp [FTQUAL_product].str = str;
                    if (StringNICmp (str, "tRNA-", 5) == 0) {
                      qvp [FTQUAL_trna_aa].str = str + 5;
                    }
                  }
                }
                qvp [FTQUAL_anticodon].slp = trna->anticodon;
                qvp [FTQUAL_trna_codons].trp = trna;
              }
            }
          } else {
            if (rrp->ext.choice == 1) {
              str = (CharPtr) rrp->ext.value.ptrvalue;
              qvp [FTQUAL_product].str = str;

              /*
              if (rrp->type == 255 && (! StringHasNoText (str))) {
                if        (StringICmp (str, "internal transcribed spacer 1") == 0 ||
                           StringICmp (str, "internal transcribed spacer ITS1") == 0 ||
                           StringICmp (str, "ITS1") == 0) {
                  qvp [FTQUAL_rrna_its].str = "ITS1";
                } else if (StringICmp (str, "internal transcribed spacer 2") == 0 ||
                           StringICmp (str, "internal transcribed spacer ITS2") == 0 ||
                           StringICmp (str, "ITS2") == 0) {
                  qvp [FTQUAL_rrna_its].str = "ITS2";
                } else if (StringICmp (str, "internal transcribed spacer 3") == 0 ||
                           StringICmp (str, "internal transcribed spacer ITS3") == 0 ||
                           StringICmp (str, "ITS3") == 0) {
                  qvp [FTQUAL_rrna_its].str = "ITS3";
                }
              }
              */
            }
          }
        }
        break;
      case SEQFEAT_REGION :
        if (afp->format == GENPEPT_FMT && featdeftype == FEATDEF_REGION && ISA_aa (bsp->mol)) {
          qvp [FTQUAL_region_name].str = (CharPtr) sfp->data.value.ptrvalue;
        } else {
          qvp [FTQUAL_region].str = (CharPtr) sfp->data.value.ptrvalue;
        }
        break;
      case SEQFEAT_COMMENT :
        break;
      case SEQFEAT_BOND :
        bondidx = (Int2) sfp->data.value.intvalue;
        if (bondidx == 255) {
          bondidx = 5;
        }
        if (bondidx > 0 && bondidx < 6) {
          if (afp->format == GENPEPT_FMT && ISA_aa (bsp->mol)) {
            qvp [FTQUAL_bond_type].str = bondList [bondidx];
          } else {
            qvp [FTQUAL_bond].str = bondList [bondidx];
          }
        }
        break;
      case SEQFEAT_SITE :
        siteidx = (Int2) sfp->data.value.intvalue;
        if (siteidx == 255) {
          siteidx = 26;
        }
        if (siteidx > 0 && siteidx < 27) {
          if (afp->format == GENPEPT_FMT && ISA_aa (bsp->mol)) {
            qvp [FTQUAL_site_type].str = siteList [siteidx];
          } else {
            qvp [FTQUAL_site].str = siteList [siteidx];
          }
        }
        break;
      case SEQFEAT_PSEC_STR :
        qvp [FTQUAL_sec_str_type].num = sfp->data.value.intvalue;
        break;
      case SEQFEAT_HET :
        qvp [FTQUAL_heterogen].str = (CharPtr) sfp->data.value.ptrvalue;
        break;
      default :
        break;
    }
  }

  /* common fields set here */

  VisitUserObjectsInUop (sfp->ext, (Pointer) qvp, RecordUserObjectsInQVP);

  if (fcontext.featdeftype == FEATDEF_repeat_region) {
    pseudo = FALSE;
  }

  qvp [FTQUAL_seqfeat_note].str = sfp->comment;

  qvp [FTQUAL_pseudo].ble = pseudo;

  /* if RELEASE_MODE, check list of features that can have /pseudo */

  if (ajp->flags.dropIllegalQuals && pseudo  &&
      (fcontext.seqfeattype == SEQFEAT_RNA || fcontext.seqfeattype == SEQFEAT_IMP) ) {
    switch (featdeftype) {

    case  FEATDEF_allele:
    case  FEATDEF_attenuator:
    case  FEATDEF_CAAT_signal:
    case  FEATDEF_conflict:
    case  FEATDEF_D_loop:
    case  FEATDEF_enhancer:
    case  FEATDEF_GC_signal:
    case  FEATDEF_iDNA:
    case  FEATDEF_intron:
    case  FEATDEF_LTR:
    case  FEATDEF_misc_binding:
    case  FEATDEF_misc_difference:
    case  FEATDEF_misc_recomb:
    case  FEATDEF_misc_RNA:
    case  FEATDEF_misc_signal:
    case  FEATDEF_misc_structure:
    case  FEATDEF_modified_base:
    case  FEATDEF_mutation:
    case  FEATDEF_old_sequence:
    case  FEATDEF_polyA_signal:
    case  FEATDEF_polyA_site:
    case  FEATDEF_precursor_RNA:
    case  FEATDEF_prim_transcript:
    case  FEATDEF_primer_bind:
    case  FEATDEF_protein_bind:
    case  FEATDEF_RBS:
    case  FEATDEF_repeat_region:
    case  FEATDEF_repeat_unit:
    case  FEATDEF_rep_origin:
    case  FEATDEF_satellite:
    case  FEATDEF_stem_loop:
    case  FEATDEF_STS:
    case  FEATDEF_TATA_signal:
    case  FEATDEF_terminator:
    case  FEATDEF_unsure:
    case  FEATDEF_variation:
    case  FEATDEF_3clip:
    case  FEATDEF_3UTR:
    case  FEATDEF_5clip:
    case  FEATDEF_5UTR:
    case  FEATDEF_10_signal:
    case  FEATDEF_35_signal:
      qvp [FTQUAL_pseudo].ble = FALSE;
        break;
    default:
        break;
    }
  }

  if (afp->format != GENPEPT_FMT) {
    qvp [FTQUAL_evidence].num = sfp->exp_ev;
  }

  /* exception currently legal only on CDS */

  if ((! ajp->flags.dropIllegalQuals) || fcontext.seqfeattype == SEQFEAT_CDREGION) {
    qvp [FTQUAL_exception].str = sfp->except_text;

    if (sfp->excpt && qvp [FTQUAL_exception].str == NULL) {
      if (qvp [FTQUAL_seqfeat_note].str != NULL) {
        if (ajp->flags.dropIllegalQuals &&
            (! StringInStringList (qvp [FTQUAL_seqfeat_note].str, validExceptionString))) {
        /* !!! if ajp->flags.dropIllegalQuals, check CDS list to avoid losing note !!! */
          qvp [FTQUAL_exception].str = NULL;
        } else {
          /* if no /exception text, use text in comment, remove from /note */

          qvp [FTQUAL_exception].str = qvp [FTQUAL_seqfeat_note].str;
          qvp [FTQUAL_seqfeat_note].str = NULL;
        }
      } else {
        qvp [FTQUAL_exception].str = "No explanation supplied";
      }

      /* !!! if ajp->flags.dropIllegalQuals, check CDS list here as well !!! */
      if (ajp->flags.dropIllegalQuals &&
          (! StringInStringList (qvp [FTQUAL_seqfeat_note].str, validExceptionString)) ) {
        qvp [FTQUAL_exception].str = NULL;
      }
    }
    if (ajp->flags.dropIllegalQuals &&
        (! StringInStringList (qvp [FTQUAL_exception].str, validExceptionString))) {
      qvp [FTQUAL_exception_note].str = qvp [FTQUAL_exception].str;
      qvp [FTQUAL_exception].str = NULL;
    }
  } else {
    qvp [FTQUAL_exception_note].str = sfp->except_text;
  }
  if (StringHasNoText (qvp [FTQUAL_exception].str)) {
    qvp [FTQUAL_exception].str = NULL;
  }
  if (StringHasNoText (qvp [FTQUAL_exception_note].str)) {
    qvp [FTQUAL_exception_note].str = NULL;
  }

  qvp [FTQUAL_db_xref].vnp = sfp->dbxref;
  qvp [FTQUAL_citation].vnp = sfp->cit;

  /* /product same as sfp->comment will suppress /note */

  if (! StringHasNoText (qvp [FTQUAL_product].str) &&
      StringICmp (sfp->comment, qvp [FTQUAL_product].str) == 0) {
    qvp [FTQUAL_seqfeat_note].str = NULL;
  }
  /* case sensitive AJ011317.1 */
  if (! StringHasNoText (qvp [FTQUAL_cds_product].str) &&
      StringCmp (sfp->comment, qvp [FTQUAL_cds_product].str) == 0) {
    qvp [FTQUAL_seqfeat_note].str = NULL;
  }

  /* /gene same as sfp->comment will suppress /note */
  /* case sensitive -gi|6572973|gb|AF195052.1|AF195052 */

  if (! StringHasNoText (qvp [FTQUAL_gene].str) &&
      StringCmp (sfp->comment, qvp [FTQUAL_gene].str) == 0) {
    qvp [FTQUAL_seqfeat_note].str = NULL;
  }

  /* gene /note same as sfp->comment will suppress /note - U92435.1 says do not do this */

  /*
  if (! StringHasNoText (qvp [FTQUAL_gene_note].str) &&
      StringICmp (sfp->comment, qvp [FTQUAL_gene_note].str) == 0) {
    qvp [FTQUAL_seqfeat_note].str = NULL;
  }
  */

  /* if site sfp->comment contains site name, suppress site in /note */

  if (! StringHasNoText (qvp [FTQUAL_site].str) &&
      StringStr (sfp->comment, qvp [FTQUAL_site].str) != NULL) {
    qvp [FTQUAL_site].str = NULL;
  }

  /* /EC_number same as sfp->comment will suppress /note */

  for (vnp = qvp [FTQUAL_prot_EC_number].vnp; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if ((! StringHasNoText (str)) &&
        StringICmp (sfp->comment, str) == 0) {
      qvp [FTQUAL_seqfeat_note].str = NULL;
    }
  }


  /* now go through gbqual list */

  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    idx = GbqualToFeaturIndex (gbq->qual);
    if (idx > 0 && idx < ASN2GNBK_TOTAL_FEATUR) {
      if (qvp [idx].gbq == NULL) {
        if (idx == FTQUAL_product_quals) {
          if (qvp [FTQUAL_product].str == NULL) {
            qvp [FTQUAL_product].str = gbq->val;
          } else if (qvp [FTQUAL_xtra_prod_quals].gbq == NULL) {
            /* chain will include remaining product gbquals */
            qvp [FTQUAL_xtra_prod_quals].gbq = gbq;
          }
        } else {
          qvp [idx].gbq = gbq;
        }
      }

    } else if (idx == 0) {

      qualclass = IllegalGbqualToClass (gbq->qual);
      if (qualclass == 0) {
        qualclass = Qual_class_quote;
      }
      tmp = StringSave (gbq->val);
      if (tmp != NULL) {
        str = MemNew (sizeof (Char) * (StringLen (gbq->val) + StringLen (tmp) + 10));
        if (str != NULL) {
          if (qualclass == Qual_class_quote) {
            if (StringIsJustQuotes (tmp)) {
              sprintf (buf, "/%s", gbq->qual);
            } else {
              ptr = tmp;
              ch = *ptr;
              while (ch != '\0') {
                if (ch == '"') {
                  *ptr = '\'';
                }
                ptr++;
                ch = *ptr;
              }
              sprintf (str, "/%s=\"%s\"", gbq->qual, tmp);
            }
            ValNodeCopyStr (&illegal, 0, str);
          } else if (qualclass == Qual_class_noquote || qualclass == Qual_class_label) {
            if (StringIsJustQuotes (tmp)) {
              sprintf (str, "/%s", gbq->qual);
            } else {
              sprintf (str, "/%s=%s", gbq->qual, tmp);
            }
            ValNodeCopyStr (&illegal, 0, str);
          }
          MemFree (str);
        }
        MemFree (tmp);
      }
    }
  }

  /* illegal qualifiers are copied and formatted in valnode chain */

  if (! ajp->flags.dropIllegalQuals) {
    qvp [FTQUAL_illegal_qual].vnp = illegal;
  }

  /* remove protein description that equals the gene name, case sensitive */

  if (StringCmp (qvp [FTQUAL_gene].str, qvp [FTQUAL_prot_desc].str) == 0) {
    qvp [FTQUAL_prot_desc].str = NULL;
  }

  /* remove protein description that equals the cds product, case sensitive */

  if (StringCmp (qvp [FTQUAL_cds_product].str, qvp [FTQUAL_prot_desc].str) == 0) {
    qvp [FTQUAL_prot_desc].str = NULL;
  }

  /* remove comment contained in prot_desc - gi|4530123|gb|AF071539.1|AF071539 */

  if (StringStr (qvp [FTQUAL_prot_desc].str, qvp [FTQUAL_seqfeat_note].str) != NULL) {
    qvp [FTQUAL_seqfeat_note].str = NULL;
  }

  /* remove protein description that equals the standard name */

  if (qvp [FTQUAL_standard_name].gbq != NULL && qvp [FTQUAL_prot_desc].str != NULL) {
    gbq = qvp [FTQUAL_standard_name].gbq;
    lasttype = gbq->qual;
    while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
      if (StringICmp (gbq->val, qvp [FTQUAL_prot_desc].str) == 0) {
        qvp [FTQUAL_prot_desc].str = NULL;
      }
      gbq = gbq->next;
    }
  }

  /* remove protein description that equals a gene synonym - case insensitive AF109216.1 */

  for (vnp = gene_syn; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if ((! StringHasNoText (str)) &&
        StringCmp (str, qvp [FTQUAL_prot_desc].str) == 0) {
      /* NC_001823 leave in prot_desc if no cds_product */
      if (qvp [FTQUAL_cds_product].str != NULL) {
        qvp [FTQUAL_prot_desc].str = NULL;
      }
    }
  }

  /* remove comment that equals a gene synonym */

  if (afp->format != GENPEPT_FMT && (! ifp->mapToProt)) {
    for (vnp = gene_syn; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if ((! StringHasNoText (str)) &&
          StringICmp (str, qvp [FTQUAL_seqfeat_note].str) == 0) {
        qvp [FTQUAL_seqfeat_note].str = NULL;
      }
    }
  }

  /* remove protein comment descriptor that equals the protein note */

  if (StringCmp (qvp [FTQUAL_prot_note].str, qvp [FTQUAL_prot_comment].str) == 0) {
    qvp [FTQUAL_prot_comment].str = NULL;
  }

  /* suppress cds comment if a subset of protein note - AF002218.1 */

  if (StringStr (qvp [FTQUAL_prot_note].str, qvp [FTQUAL_seqfeat_note].str) != NULL) {
    qvp [FTQUAL_seqfeat_note].str = NULL;
  }

  /* suppress selenocysteine note if already in comment */

  if (StringStr (sfp->comment, "selenocysteine") != NULL) {
    qvp [FTQUAL_selenocysteine].str = NULL;
  }

  /* now print qualifiers from table */

#ifdef DISPLAY_STRINGS
  s_DisplayQVP(qvp, feat_note_order);
#endif

  /* Strip duplicate notes */

  if ((StringCmp(qvp[FTQUAL_product].str,
         qvp[FTQUAL_seqfeat_note].str) == 0)) {
    qvp[FTQUAL_seqfeat_note].str = NULL;
  }

  if ((qvp[FTQUAL_standard_name].gbq != NULL) &&
      (qvp[FTQUAL_standard_name].gbq->val != NULL)) {
    if ((StringCmp(qvp[FTQUAL_seqfeat_note].str,
           qvp[FTQUAL_standard_name].gbq->val) == 0)) {
      qvp[FTQUAL_seqfeat_note].str = NULL;
    }
  }

  /* Display strings for debugging purposes */

#ifdef DISPLAY_STRINGS
  s_DisplayQVP(qvp, feat_qual_order);
#endif

  /* optionally populate indexes for NCBI internal database */

  if (index != NULL) {
    if (! StringHasNoText (qvp [FTQUAL_gene].str)) {
      ValNodeCopyStrToHead (&(index->genes), 0, qvp [FTQUAL_gene].str);
    }
  }

  FFAddOneChar(ffstring, '\n', FALSE);

  /* Build the flat file */
  FormatFeatureBlockQuals (ffstring, ajp, afp, asp, bsp, featdeftype, gene_syn,
                           lasttype, location, prod,
                           protein_pid_g, qvp,
                           left, right, strand,
                           sfp, target, is_other,
                           is_journalscan, is_gps);

  /* ??? and then deal with the various note types separately (not in order table) ??? */

  /* free aa-dna or dna-aa mapped location */

  SeqLocFree (loc);

  ValNodeFreeData (illegal);

  BioseqUnlock (unlockme);

  str = FFEndPrint (ajp, ffstring, afp->format, 21, 21, 21, 21, "FT");

  /* optionally populate gbseq for XML-ized GenBank format */

  if (gbseq != NULL) {
    if (gbfeat != NULL) {
      AddFeatureToGbseq (gbseq, gbfeat, str);
    }
  }

  FFRecycleString(ajp, ffstring);
  return str;
}



static CharPtr FormatSequenceBlock (
  Asn2gbFormatPtr afp,
  BaseBlockPtr bbp
)

{
  IntAsn2gbJobPtr   ajp;
  Asn2gbSectPtr     asp;
  Byte              bases [400];
  Int2              blk;
  BioseqPtr         bsp;
  Char              buf [80];
  Int2              cnt;
  Uint1             code = Seq_code_iupacna;
  Int2              count;
  Int2              ctr;
  GBSeqPtr          gbseq;
  Int2              i;
  IntAsn2gbSectPtr  iasp;
  Boolean           is_na;
  Int2              lin;
  Int4              pos;
  Uint1             residue;
  SeqBlockPtr       sbp;
  SeqPortPtr        spp;
  Int4              start;
  Int4              stop;
  CharPtr           str;
  StringItemPtr     ffstring;

  if (afp == NULL || bbp == NULL) return NULL;
  sbp = (SeqBlockPtr) bbp;
  ajp = afp->ajp;
  if (ajp == NULL) return NULL;
  asp = afp->asp;
  if (asp == NULL) return NULL;
  iasp = (IntAsn2gbSectPtr) asp;
  bsp = (asp->bsp);
  if (bsp == NULL) return NULL;

  spp = iasp->spp;
  if (spp == NULL) {

    /* if first time, create SeqPort for this section */

    if (ISA_aa (bsp->mol)) {
      if (ajp->flags.iupacaaOnly) {
        code = Seq_code_iupacaa;
      } else {
        code = Seq_code_ncbieaa;
      }
    }

    if (ajp->ajp.slp != NULL) {
      spp = SeqPortNewByLoc (ajp->ajp.slp, code);
    } else {
      spp = SeqPortNew (bsp, 0, -1, 0, code);
    }
    if (spp == NULL) return NULL;
    if (bsp->repr == Seq_repr_delta || bsp->repr == Seq_repr_virtual) {
      SeqPortSet_do_virtual (spp, TRUE);
    }

    iasp->spp = spp;
  }

  start = sbp->start;
  stop = sbp->stop;

  if (start != spp->curpos) {
    SeqPortSeek (spp, start, SEEK_SET);
  }

  pos = start;

  count = 0;
  cnt = 0;
  blk = 0;
  lin = 0;

  is_na = ISA_na (bsp->mol);

  ctr = (Int2) MIN ((Int4) (stop - pos), (Int4) sizeof (bases));
  ctr = SeqPortRead (spp, bases, ctr);

  i = 0;

  if (ctr < 0) {
    residue = -ctr;
  } else if (ctr < 1) {
    residue = SEQPORT_EOF;
  } else {
    residue = (Uint1) bases [i];
  }


  ffstring = FFGetString(ajp);
  while (pos < stop && residue != SEQPORT_EOF) {

    if (residue == INVALID_RESIDUE) {
      if (is_na) {
        residue = 'N';
      } else {
        residue = 'X';
      }
    }

    if (IS_residue (residue)) {

      buf [count] = (Char) (TO_LOWER (residue));
      count++;
      cnt++;
      pos++;

      blk++;
      lin++;
      if (lin >= 60) {

        buf [count] = '\0';
        PrintSeqLine (ffstring, afp->format, buf, start, start + cnt);
        count = 0;
        cnt = 0;
        blk = 0;
        lin = 0;
        start += 60;

      } else if (blk >= 10) {

        buf [count] = ' ';
        count++;
        blk = 0;

      }
    }

    i++;
    if (i >= ctr) {
      i = 0;
      ctr = (Int2) MIN ((Int4) (stop - pos), (Int4) sizeof (bases));
      ctr = SeqPortRead (spp, bases, ctr);
      if (ctr < 0) {
        bases [0] = -ctr;
      } else if (ctr < 1) {
        bases [0] = SEQPORT_EOF;
      }
    }
    residue = (Uint1) bases [i];
  }

  buf [count] = '\0';
  if (count > 0) {
    PrintSeqLine (ffstring, afp->format, buf, start, start + cnt);
  }

  if (ajp->transientSeqPort) {
    iasp->spp = SeqPortFree (iasp->spp);
  }

  str = FFToCharPtr(ffstring);

  /* optionally populate gbseq for XML-ized GenBank format */

  if (ajp->gbseq) {
    gbseq = &asp->gbseq;
  } else {
    gbseq = NULL;
  }

  if (gbseq != NULL) {
    CatenateSequenceInGbseq (gbseq, str);
  }

  FFRecycleString(ajp, ffstring);
  return str;
}



/* ********************************************************************** */

/* public functions */

static int LIBCALLBACK SortParagraphByIDProc (
  VoidPtr vp1,
  VoidPtr vp2
)

{
  BaseBlockPtr  bbp1, bbp2;

  if (vp1 == NULL || vp2 == NULL) return 0;
  bbp1 = *((BaseBlockPtr PNTR) vp1);
  bbp2 = *((BaseBlockPtr PNTR) vp2);
  if (bbp1 == NULL || bbp2 == NULL) return 0;

  if (bbp1->entityID > bbp2->entityID) return 1;
  if (bbp1->entityID < bbp2->entityID) return -1;

  if (bbp1->itemtype > bbp2->itemtype) return 1;
  if (bbp1->itemtype < bbp2->itemtype) return -1;

  if (bbp1->itemID > bbp2->itemID) return 1;
  if (bbp1->itemID < bbp2->itemID) return -1;

  if (bbp1->paragraph > bbp2->paragraph) return 1;
  if (bbp1->paragraph < bbp2->paragraph) return -1;

  return 0;
}

typedef struct modeflags {
  Boolean  flags [22];
} ModeFlags, PNTR ModeFlagsPtr;

static ModeFlags flagTable [] = {

  /* RELEASE_MODE */
  {TRUE,  TRUE, TRUE, TRUE, TRUE,
   TRUE, TRUE, TRUE, TRUE, TRUE,
   TRUE,  TRUE, TRUE, TRUE, TRUE,
   TRUE,  TRUE, TRUE, TRUE, TRUE,
   TRUE, TRUE},

  /* ENTREZ_MODE */
  {FALSE, TRUE, TRUE, TRUE, TRUE,
   FALSE, TRUE, TRUE, TRUE, TRUE,
   TRUE,  TRUE, FALSE, TRUE, TRUE,
   TRUE,  TRUE, FALSE, FALSE, TRUE,
   TRUE, FALSE},

  /* SEQUIN_MODE */
  {FALSE, FALSE, FALSE, FALSE, FALSE,
   FALSE, FALSE, TRUE, FALSE, FALSE,
   FALSE, FALSE, FALSE, FALSE, FALSE,
   FALSE, FALSE, FALSE, FALSE, FALSE,
   FALSE, FALSE},

  /* DUMP_MODE */
  {FALSE, FALSE, FALSE, FALSE, FALSE,
   FALSE, FALSE, FALSE, FALSE, FALSE,
   FALSE, FALSE, FALSE, FALSE, FALSE,
   FALSE, FALSE, FALSE, FALSE, FALSE,
   FALSE, FALSE}
};

static void SetFlagsFromMode (
  IntAsn2gbJobPtr  ajp,
  ModType mode
)

{
  BoolPtr       bp;
  ModeFlagsPtr  mfp;

  if (ajp == NULL) return;
  if (! (mode >= RELEASE_MODE && mode <= DUMP_MODE)) {
    mode = DUMP_MODE;
  }
  mfp = &(flagTable [(int) (mode - 1)]);
  bp = &(mfp->flags [0]);

  ajp->flags.suppressLocalID = *(bp++);
  ajp->flags.validateFeats = *(bp++);
  ajp->flags.ignorePatPubs = *(bp++);
  ajp->flags.dropShortAA = *(bp++);
  ajp->flags.avoidLocusColl = *(bp++);

  ajp->flags.iupacaaOnly = *(bp++);
  ajp->flags.dropBadCitGens = *(bp++);
  ajp->flags.noAffilOnUnpub = *(bp++);
  ajp->flags.dropIllegalQuals = *(bp++);
  ajp->flags.checkQualSyntax = *(bp++);

  ajp->flags.needRequiredQuals = *(bp++);
  ajp->flags.needOrganismQual = *(bp++);
  ajp->flags.needAtLeastOneRef = *(bp++);
  ajp->flags.citArtIsoJta = *(bp++);
  ajp->flags.dropBadDbxref = *(bp++);

  ajp->flags.useEmblMolType = *(bp++);
  ajp->flags.hideBankItComment = *(bp++);
  ajp->flags.checkCDSproductID = *(bp++);
  ajp->flags.suppressSegLoc = *(bp++);
  ajp->flags.srcQualsToNote = *(bp)++;

  ajp->flags.hideEmptySource = *(bp++);
  ajp->flags.forGbRelease = *(bp++);

  /* collaboration unapproved source quals on their own line only in indexer Sequin */

  if (GetAppProperty ("InternalNcbiSequin") == NULL) {
    ajp->flags.srcQualsToNote = TRUE;
  }
}

static void CheckVersionWithGi (BioseqPtr bsp, Pointer userdata)

{
  Boolean       hasGi = FALSE;
  BoolPtr       missingVersion;
  SeqIdPtr      sip;
  TextSeqIdPtr  tsip;
  Boolean       zeroVersion = FALSE;

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    switch (sip->choice) {
      case SEQID_TPG:
      case SEQID_TPE:
      case SEQID_TPD:
      case SEQID_GENBANK:
      case SEQID_EMBL:
      case SEQID_DDBJ:
        tsip = (TextSeqIdPtr) sip->data.ptrvalue;
        if (tsip != NULL && tsip->version == 0) {
          zeroVersion = TRUE;
        }
        break;
      case SEQID_GI :
        hasGi = TRUE;
        break;
      default :
        break;
    }
  }
  if (hasGi && zeroVersion) {
    missingVersion = (BoolPtr) userdata;
    *missingVersion = TRUE;
  }
}


#define FEAT_FETCH_MASK (ONLY_NEAR_FEATURES | FAR_FEATURES_SUPPRESS)
#define HTML_XML_ASN_MASK (CREATE_HTML_FLATFILE | CREATE_XML_GBSEQ_FILE | CREATE_ASN_GBSEQ_FILE)

NLM_EXTERN Asn2gbJobPtr asn2gnbk_setup (
  BioseqPtr bsp,
  BioseqSetPtr bssp,
  SeqLocPtr slp,
  FmtType format,
  ModType mode,
  StlType style,
  FlgType flags,
  LckType locks,
  XtraPtr extra
)

{
  IntAsn2gbJobPtr  ajp = NULL;
  Asn2gbSectPtr    asp;
  Asn2gbWork       aw;
  BaseBlockPtr     bbp;
  BaseBlockPtr     PNTR blockArray;
  Uint2            entityID = 0;
  GBSeqPtr         gbseq = NULL;
  Int4             i;
  IndxPtr          index = NULL;
  Int4             j;
  Int4             k;
  Boolean          lockFarComp;
  Boolean          lockFarLocs;
  Boolean          lockFarProd;
  Boolean          lookupFarComp;
  Boolean          lookupFarHist;
  Boolean          lookupFarLocs;
  Boolean          lookupFarProd;
  Boolean          missingVersion;
  Int4             numBlocks;
  Int4             numSections;
  SeqEntryPtr      oldscope;
  ObjMgrDataPtr    omdp;
  BaseBlockPtr     PNTR paragraphArray;
  BaseBlockPtr     PNTR paragraphByIDs;
  Int4             numParagraphs;
  Asn2gbSectPtr    PNTR sectionArray;
  SubmitBlockPtr   sbp;
  SeqEntryPtr      sep;
  SeqIntPtr        sintp;
  SeqSubmitPtr     ssp;
  BioseqSetPtr     topbssp;
  ValNodePtr       vnp;
  Boolean          is_html = FALSE;

  if (extra != NULL) {
    index = extra->index;
    gbseq = extra->gbseq;
  }

  if (slp != NULL) {
    bsp = BioseqFind (SeqLocId (slp));
    if (bsp == NULL) {
      bsp = BioseqFindFromSeqLoc (slp);
    }
    if (bsp == NULL) return NULL;

    /* if location is whole, generate normal bioseq report */

    if (slp->choice == SEQLOC_WHOLE) {
      slp = NULL;
    } else if (slp->choice == SEQLOC_INT) {
      sintp = (SeqIntPtr) slp->data.ptrvalue;
      if (sintp != NULL &&
          sintp->from == 0 &&
          sintp->to == bsp->length - 1 &&
          sintp->strand == Seq_strand_plus) {
        slp = NULL;
      }
    }
  }

  if (bsp != NULL) {
    bssp = NULL;
    entityID = ObjMgrGetEntityIDForPointer (bsp);
  } else if (bssp != NULL) {
    entityID = ObjMgrGetEntityIDForPointer (bssp);
  }

  if (entityID == 0) return NULL;

  if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
    SeqMgrIndexFeatures (entityID, NULL);
  }

  if (mode == RELEASE_MODE) {
    sep = GetTopSeqEntryForEntityID (entityID);
    missingVersion = FALSE;
    VisitBioseqsInSep (sep, (Pointer) &missingVersion, CheckVersionWithGi);
    if (missingVersion) return NULL;
  }

  ajp = (IntAsn2gbJobPtr) MemNew (sizeof (IntAsn2gbJob));
  if (ajp == NULL) return NULL;

  is_html = (Boolean) ((flags & HTML_XML_ASN_MASK) == CREATE_HTML_FLATFILE);
  if (is_html) {
    InitWWW(ajp);
  }

  init_buff ();

  asn2ff_set_output (NULL, "\n");

  ajp->ajp.entityID = entityID;
  ajp->ajp.bsp = bsp;
  ajp->ajp.bssp = bssp;
  if (slp != NULL) {
    ajp->ajp.slp = AsnIoMemCopy ((Pointer) slp,
                                 (AsnReadFunc) SeqLocAsnRead,
                                 (AsnWriteFunc) SeqLocAsnWrite);
  } else {
    ajp->ajp.slp = NULL;
  }

  /* if location specified, normal defaults to master style */

  if (ajp->ajp.slp != NULL && style == NORMAL_STYLE) {
    style = MASTER_STYLE;
  }

  ajp->format = format;

  SetFlagsFromMode (ajp, mode);

  ajp->transientSeqPort = (Boolean) ((locks & FREE_SEQPORT_EACH_TIME) != 0);

  lockFarComp = (Boolean) ((locks & LOCK_FAR_COMPONENTS) != 0);
  lockFarLocs = (Boolean) ((locks & LOCK_FAR_LOCATIONS) != 0);
  lockFarProd = (Boolean) ((locks & LOCK_FAR_PRODUCTS) != 0);

  if (lockFarComp || lockFarLocs || lockFarProd) {

    /* lock all bioseqs in advance, including remote genome components */

    sep = GetTopSeqEntryForEntityID (entityID);
    ajp->lockedBspList = LockFarComponentsEx (sep, lockFarComp, lockFarLocs, lockFarProd);
  }

  lookupFarComp = (Boolean) ((locks & LOOKUP_FAR_COMPONENTS) != 0);
  lookupFarLocs = (Boolean) ((locks & LOOKUP_FAR_LOCATIONS) != 0);
  lookupFarProd = (Boolean) ((locks & LOOKUP_FAR_PRODUCTS) != 0);
  lookupFarHist = (Boolean) ((locks & LOOKUP_FAR_HISTORY) != 0);

  if (lookupFarComp || lookupFarLocs || lookupFarProd || lookupFarHist) {

    /* lookukp all far SeqIDs in advance */

    sep = GetTopSeqEntryForEntityID (entityID);
    LookupFarSeqIDs (sep, lookupFarComp, lookupFarLocs, lookupFarProd, FALSE, lookupFarHist);
  }

  ajp->showFarTransl = (Boolean) ((flags & SHOW_FAR_TRANSLATION) != 0);
  ajp->transIfNoProd = (Boolean) ((flags & TRANSLATE_IF_NO_PRODUCT) != 0);
  ajp->alwaysTranslCds = (Boolean) ((flags & ALWAYS_TRANSLATE_CDS) != 0);

  ajp->masterStyle = (Boolean) (style == MASTER_STYLE);
  if (format == GENBANK_FMT || format == GENPEPT_FMT) {
    ajp->newSourceOrg = (Boolean) ((flags & USE_OLD_SOURCE_ORG) == 0);
  }

  ajp->relModeError = FALSE;

  ajp->index = index;
  ajp->gbseq = gbseq;

  MemSet ((Pointer) (&aw), 0, sizeof (Asn2gbWork));
  aw.ajp = ajp;
  aw.entityID = entityID;

  aw.sectionList = NULL;
  aw.lastsection = NULL;

  aw.currsection = 0;
  aw.showAllFeats = FALSE;

  aw.showconfeats = (Boolean) ((flags & SHOW_CONTIG_FEATURES) != 0);
  aw.showconsource = (Boolean) ((flags & SHOW_CONTIG_SOURCES) != 0);

  aw.onlyNearFeats = (Boolean) ((flags & FEAT_FETCH_MASK) == ONLY_NEAR_FEATURES);
  aw.farFeatsSuppress = (Boolean) ((flags & FEAT_FETCH_MASK) == FAR_FEATURES_SUPPRESS);

  if (aw.onlyNearFeats && aw.farFeatsSuppress) {
    aw.farFeatsSuppress = FALSE;
  }

  aw.hideImpFeats = (Boolean) ((flags & HIDE_IMP_FEATS) != 0);
  aw.hideSnpFeats = (Boolean) ((flags & HIDE_SNP_FEATS) != 0);

  aw.isGPS = FALSE;
  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep != NULL && IS_Bioseq_set (sep)) {
    topbssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (topbssp != NULL && topbssp->_class == BioseqseqSet_class_gen_prod_set) {
      aw.isGPS = TRUE;
      aw.copyGpsCdsUp = (Boolean) ((flags & COPY_GPS_CDS_UP) != 0);
      aw.copyGpsGeneDown = (Boolean) ((flags & COPY_GPS_GENE_DOWN) != 0);
    }
  }

  aw.newLocusLine = TRUE;

  if ((Boolean) (flags & DDBJ_VARIANT_FORMAT) != 0) {
    aw.citSubsFirst = TRUE;
    aw.hideGeneFeats = TRUE;
    aw.newLocusLine = FALSE;
    ajp->newSourceOrg = FALSE;
  }

  aw.format = format;
  aw.mode = mode;
  aw.style = style;

  aw.hup = FALSE;
  aw.ssp = NULL;

  aw.failed = FALSE;

  omdp = ObjMgrGetData (entityID);
  if (omdp != NULL && omdp->datatype == OBJ_SEQSUB) {
    ssp = (SeqSubmitPtr) omdp->dataptr;
    if (ssp != NULL && ssp->datatype == 1) {
      aw.ssp = ssp;
      sbp = ssp->sub;
      if (sbp != NULL) {
        aw.hup = sbp->hup;
      }
    }
  }

  oldscope = SeqEntrySetScope (sep);

  if (bssp != NULL) {

    /* handle all components of a pop/phy/mut/eco set */

    sep = SeqMgrGetSeqEntryForData (bssp);
    DoOneBioseqSet (sep, &aw);

  } else {

    /* handle single bioseq, which may be segmented or a local part */

    DoOneBioseq (bsp, &aw);
  }

  SeqEntrySetScope (oldscope);

  /* check for failure to populate anything */

  if (ajp->flags.needAtLeastOneRef && aw.failed) return asn2gnbk_cleanup ((Asn2gbJobPtr) ajp);

  numSections = ValNodeLen (aw.sectionList);
  ajp->ajp.numSections = numSections;

  if (numSections == 0) return asn2gnbk_cleanup ((Asn2gbJobPtr) ajp);

  /* allocate section array for this job */

  sectionArray = (Asn2gbSectPtr PNTR) MemNew (sizeof (Asn2gbSectPtr) * (numSections + 1));
  ajp->ajp.sectionArray = sectionArray;

  if (sectionArray == NULL) return asn2gnbk_cleanup ((Asn2gbJobPtr) ajp);

  /* fill in section and paragraph arrays */

  numParagraphs = 0;
  for (vnp = aw.sectionList, i = 0; vnp != NULL && i < numSections; vnp = vnp->next, i++) {
    asp = (Asn2gbSectPtr) vnp->data.ptrvalue;
    sectionArray [i] = asp;
    if (asp != NULL) {
      numParagraphs += asp->numBlocks;
    }
  }

  /* allocate paragraph array pointing to all blocks in all sections */

  ajp->ajp.numParagraphs = numParagraphs;
  if (numParagraphs == 0) return asn2gnbk_cleanup ((Asn2gbJobPtr) ajp);

  paragraphArray = (BaseBlockPtr PNTR) MemNew (sizeof (BaseBlockPtr) * (numParagraphs + 1));
  ajp->ajp.paragraphArray = paragraphArray;

  paragraphByIDs = (BaseBlockPtr PNTR) MemNew (sizeof (BaseBlockPtr) * (numParagraphs + 1));
  ajp->ajp.paragraphByIDs = paragraphByIDs;

  if (paragraphArray == NULL || paragraphByIDs == NULL) return asn2gnbk_cleanup ((Asn2gbJobPtr) ajp);

  k = 0;
  for (i = 0; i < numSections; i++) {
    asp = sectionArray [i];
    if (asp != NULL) {

      numBlocks = asp->numBlocks;
      blockArray = asp->blockArray;
      if (blockArray != NULL) {

        for (j = 0; j < numBlocks; j++) {
          bbp = blockArray [j];

          paragraphArray [k] = bbp;
          paragraphByIDs [k] = bbp;
          bbp->paragraph = k;
          k++;
        }
      }
    }
  }

  /* sort paragraphByIDs array by entityID/itemtype/itemID/paragraph */

  HeapSort (paragraphByIDs, (size_t) numParagraphs, sizeof (BaseBlockPtr), SortParagraphByIDProc);

  /* free sectionList, but leave data, now pointed to by sectionArray elements */

  ValNodeFree (aw.sectionList);

  return (Asn2gbJobPtr) ajp;
}

typedef CharPtr (*FormatProc) (Asn2gbFormatPtr afp, BaseBlockPtr bbp);

static FormatProc asn2gnbk_fmt_functions [27] = {
  NULL, NULL,
  DefaultFormatBlock, DefaultFormatBlock, DefaultFormatBlock,
  DefaultFormatBlock, DefaultFormatBlock, DefaultFormatBlock,
  DefaultFormatBlock, DefaultFormatBlock, DefaultFormatBlock,
  FormatSourceBlock, FormatOrganismBlock, FormatReferenceBlock,
  DefaultFormatBlock, FormatCommentBlock, DefaultFormatBlock,
  FormatSourceFeatBlock, FormatFeatureBlock, FormatBasecountBlock,
  DefaultFormatBlock, FormatSequenceBlock, FormatContigBlock,
  DefaultFormatBlock, DefaultFormatBlock, FormatSlashBlock,
  NULL
};

static void PrintFtableIntervals (
  ValNodePtr PNTR head,
  BioseqPtr target,
  SeqLocPtr location,
  CharPtr label
)

{
  Boolean    partial5;
  Boolean    partial3;
  SeqLocPtr  slp;
  Int4       start;
  Int4       stop;
  Char       str [64];
  Char       str1 [32];
  Char       str2 [32];

  if (head == NULL || target == NULL || location == NULL || label == NULL) return;

  slp = SeqLocFindNext (location, NULL);
  if (slp == NULL) return;

  start = GetOffsetInBioseq (slp, target, SEQLOC_START) + 1;
  stop = GetOffsetInBioseq (slp, target, SEQLOC_STOP) + 1;
  CheckSeqLocForPartial (slp, &partial5, &partial3);
  if (partial5) {
    sprintf (str1, "<%ld", (long) start);
  } else {
    sprintf (str1, "%ld", (long) start);
  }
  if (partial3) {
    sprintf (str2, ">%ld", (long) stop);
  } else {
    sprintf (str2, "%ld", (long) stop);
  }
  sprintf (str, "%s\t%s\t%s\n", str1, str2, label);
  ValNodeCopyStr (head, 0, str);

  while ((slp = SeqLocFindNext (location, slp)) != NULL) {
    start = GetOffsetInBioseq (slp, target, SEQLOC_START) + 1;
    stop = GetOffsetInBioseq (slp, target, SEQLOC_STOP) + 1;
    CheckSeqLocForPartial (slp, &partial5, &partial3);
    if (partial5) {
      sprintf (str1, "<%ld", (long) start);
    } else {
      sprintf (str1, "%ld", (long) start);
    }
    if (partial3) {
      sprintf (str2, ">%ld", (long) stop);
    } else {
      sprintf (str2, "%ld", (long) stop);
    }
    if (start != 0 && stop != 0) {
      sprintf (str, "%s\t%s\n", str1, str2);
      ValNodeCopyStr (head, 0, str);
    }
  }
}

static CharPtr goQualList [] = {
  "", "go_process", "go_component", "go_function", NULL
};

static CharPtr goQualType [] = {
  "", "Process", "Component", "Function", NULL
};

static CharPtr goFieldType [] = {
  "", "text string", "go id", "pubmed id", "evidence", NULL
};

static void PrintFTUserFld (
  UserFieldPtr ufp,
  Pointer userdata
)

{
  UserFieldPtr     entry;
  CharPtr          evidence = NULL;
  CharPtr          goid = NULL;
  ValNodePtr PNTR  head;
  Int2             i;
  Int2             j;
  size_t           len;
  ObjectIdPtr      oip;
  Int4             pmid = 0;
  CharPtr          str;
  CharPtr          textstr = NULL;
  Char             tmp [16];

  if (ufp == NULL || ufp->choice != 11) return;
  oip = ufp->label;
  if (oip == NULL) return;
  for (i = 0; goQualType [i] != NULL; i++) {
    if (StringICmp (oip->str, goQualType [i]) == 0) break;
  }
  if (goQualType [i] == NULL) return;

  entry = ufp->data.ptrvalue;
  if (entry == NULL || entry->choice != 11) return;

  for (ufp = (UserFieldPtr)  entry->data.ptrvalue; ufp != NULL; ufp = ufp->next) {
    oip = ufp->label;
    if (oip == NULL) continue;
    for (j = 0; goFieldType [j] != NULL; j++) {
      if (StringICmp (oip->str, goFieldType [j]) == 0) break;
    }
    if (goFieldType [j] == NULL) continue;
    switch (j) {
      case 1 :
        textstr = (CharPtr) ufp->data.ptrvalue;
        break;
      case 2 :
        goid = (CharPtr) ufp->data.ptrvalue;
        break;
      case 3 :
        pmid = (Int4) ufp->data.intvalue;
        break;
      case 4 :
        evidence = (CharPtr) ufp->data.ptrvalue;
        break;
      default :
        break;
    }
  }
  if (StringHasNoText (textstr)) return;

  str = (CharPtr) MemNew (StringLen (textstr) + StringLen (goid) + StringLen (evidence) + 40);
  if (str == NULL) return;
  StringCpy (str, "\t\t\t");
  StringCat (str, goQualList [i]);
  StringCat (str, "\t");
  StringCat (str, textstr);
  if (! StringHasNoText (goid)) {
    StringCat (str, "|");
    StringCat (str, goid);
  } else {
    StringCat (str, "|");
  }
  if (pmid != 0) {
    sprintf (tmp, "|%ld", (long) pmid);
    StringCat (str, tmp);
  } else {
    StringCat (str, "|");
  }
  if (! StringHasNoText (evidence)) {
    StringCat (str, "|");
    StringCat (str, evidence);
  }
  len = StringLen (str);
  while (len > 0 && str [len - 1] == '|') {
    str [len - 1] = '\0';
    len--;
  }

  head = (ValNodePtr PNTR) userdata;
  StringCat (str, "\n");
  ValNodeCopyStr (head, 0, str);
}

static void PrintFTUserObj (
  UserObjectPtr uop,
  Pointer userdata
)

{
  ObjectIdPtr  oip;

  if (uop == NULL) return;
  oip = uop->type;
  if (oip == NULL || StringICmp (oip->str, "GeneOntology") != 0) return;
  VisitUserFieldsInUop (uop, userdata, PrintFTUserFld);
}

static void PrintFTCodeBreak (
  ValNodePtr PNTR head,
  CodeBreakPtr cbp,
  BioseqPtr target
)

{
  Char             buf [80];
  Choice           cbaa;
  IntAsn2gbJob     iaj;
  CharPtr          ptr;
  Uint1            residue;
  SeqCodeTablePtr  sctp;
  Uint1            seqcode;
  SeqLocPtr        slp;
  CharPtr          str;

  seqcode = 0;
  sctp = NULL;
  cbaa = cbp->aa;
  switch (cbaa.choice) {
    case 1 :
      seqcode = Seq_code_ncbieaa;
      break;
    case 2 :
      seqcode = Seq_code_ncbi8aa;
      break;
    case 3 :
      seqcode = Seq_code_ncbistdaa;
      break;
    default :
      break;
  }
  if (seqcode == 0) return;
  sctp = SeqCodeTableFind (seqcode);
  if (sctp == NULL) return;

  MemSet ((Pointer) &iaj, 0, sizeof (IntAsn2gbJob));
  iaj.flags.iupacaaOnly = FALSE;
  iaj.relModeError = FALSE;

  slp = SeqLocFindNext (cbp->loc, NULL);
  while (slp != NULL) {
    str = FlatLoc (&iaj, target, slp, FALSE);
    if (str != NULL) {
      residue = cbaa.value.intvalue;
      ptr = Get3LetterSymbol (&iaj, seqcode, sctp, residue);
      if (ptr == NULL) {
        ptr = "OTHER";
      }
      sprintf (buf, "\t\t\ttransl_except\t(pos:%s,aa:%s)\n", str, ptr);
      ValNodeCopyStr (head, 0, buf);
      MemFree (str);
    }
    slp = SeqLocFindNext (cbp->loc, slp);
  }
}

static void PrintFtableLocAndQuals (
  IntAsn2gbJobPtr ajp,
  ValNodePtr PNTR head,
  BioseqPtr target,
  SeqFeatPtr sfp,
  SeqMgrFeatContextPtr context
)

{
  CodeBreakPtr  cbp;
  CdRegionPtr   crp;
  DbtagPtr      dbt;
  GBQualPtr     gbq;
  ValNodePtr    geneorprotdb;
  GeneRefPtr    grp;
  CharPtr       label;
  ObjectIdPtr   oip;
  BioseqPtr     prod;
  SeqFeatPtr    prot;
  ProtRefPtr    prp;
  Boolean       pseudo;
  RnaRefPtr     rrp;
  SeqIdPtr      sip;
  Char          str [256];
  Char          tmp [300];
  tRNAPtr       trp;
  ValNodePtr    vnp;

  if (head == NULL || target == NULL || sfp == NULL || context == NULL) return;
  label = (CharPtr) FeatDefTypeLabel (sfp);
  if (StringCmp (label, "Gene") == 0) {
    label = "gene";
  }
  if (StringHasNoText (label)) {
    label = "???";
  }

  PrintFtableIntervals (head, target, sfp->location, label);

  geneorprotdb = NULL;
  pseudo = sfp->pseudo;

  switch (context->seqfeattype) {
    case SEQFEAT_GENE :
      grp = (GeneRefPtr) sfp->data.value.ptrvalue;
      if (grp != NULL) {
        geneorprotdb = grp->db;
        pseudo |= grp->pseudo;

        StringNCpy_0 (str, (CharPtr) grp->locus, sizeof (str));
        if (! StringHasNoText (str)) {
          sprintf (tmp, "\t\t\tgene\t%s\n", str);
          ValNodeCopyStr (head, 0, tmp);
        }
        for (vnp = grp->syn; vnp != NULL; vnp = vnp->next) {
          StringNCpy_0 (str, (CharPtr) vnp->data.ptrvalue, sizeof (str));
          if (! StringHasNoText (str)) {
            sprintf (tmp, "\t\t\tgene_syn\t%s\n", str);
            ValNodeCopyStr (head, 0, tmp);
          }
        }
        if (! StringHasNoText (grp->maploc)) {
          sprintf (tmp, "\t\t\tmap\t%s\n", grp->maploc);
          ValNodeCopyStr (head, 0, tmp);
        }
        if (! StringHasNoText (grp->locus_tag)) {
          sprintf (tmp, "\t\t\tlocus_tag\t%s\n", grp->locus_tag);
          ValNodeCopyStr (head, 0, tmp);
        }
      }
      break;
    case SEQFEAT_CDREGION :
      prod = BioseqFind (SeqLocId (sfp->product));
      prot = SeqMgrGetBestProteinFeature (prod, NULL);
      if (prot != NULL) {
        prp = (ProtRefPtr) prot->data.value.ptrvalue;
        if (prp != NULL) {
          geneorprotdb = prp->db;
          if (prp->name != NULL) {
            for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
              StringNCpy_0 (str, (CharPtr) vnp->data.ptrvalue, sizeof (str));
              if (! StringHasNoText (str)) {
                sprintf (tmp, "\t\t\tproduct\t%s\n", str);
                ValNodeCopyStr (head, 0, tmp);
              }
            }
          }
          if (prp->desc != NULL) {
            StringNCpy_0 (str, prp->desc, sizeof (str));
            if (! StringHasNoText (str)) {
              sprintf (tmp, "\t\t\tprot_desc\t%s\n", str);
              ValNodeCopyStr (head, 0, tmp);
            }
          }
          for (vnp = prp->activity; vnp != NULL; vnp = vnp->next) {
            StringNCpy_0 (str, (CharPtr) vnp->data.ptrvalue, sizeof (str));
            if (! StringHasNoText (str)) {
              sprintf (tmp, "\t\t\tfunction\t%s\n", str);
              ValNodeCopyStr (head, 0, tmp);
            }
          }
          for (vnp = prp->ec; vnp != NULL; vnp = vnp->next) {
            StringNCpy_0 (str, (CharPtr) vnp->data.ptrvalue, sizeof (str));
            if (! StringHasNoText (str)) {
              sprintf (tmp, "\t\t\tEC_number\t%s\n", str);
              ValNodeCopyStr (head, 0, tmp);
            }
          }
        }
        StringNCpy_0 (str, prot->comment, sizeof (str));
        if (! StringHasNoText (str)) {
          sprintf (tmp, "\t\t\tprot_note\t%s\n", str);
          ValNodeCopyStr (head, 0, tmp);
        }
      }
      crp = (CdRegionPtr) sfp->data.value.ptrvalue;
      if (crp != NULL) {
        if (crp->frame > 1 && crp->frame <= 3) {
          sprintf (tmp, "\t\t\tcodon_start\t%d\n", (int) crp->frame);
          ValNodeCopyStr (head, 0, tmp);
        }
        for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next) {
          PrintFTCodeBreak (head, cbp, target);
        }
      }
      if (prod != NULL) {
        for (sip = prod->id; sip != NULL; sip = sip->next) {
          if (sip->choice == SEQID_GENBANK ||
              sip->choice == SEQID_EMBL ||
              sip->choice == SEQID_DDBJ ||
              sip->choice == SEQID_OTHER ||
              sip->choice == SEQID_TPG ||
              sip->choice == SEQID_TPE ||
              sip->choice == SEQID_TPD) {
            if (SeqIdWrite (sip, str, PRINTID_TEXTID_ACC_VER, sizeof (str)) != NULL) {
              sprintf (tmp, "\t\t\tprotein_id\t%s\n", str);
              ValNodeCopyStr (head, 0, tmp);
            }
          } else if (sip->choice == SEQID_LOCAL && (! ajp->flags.suppressLocalID)) {
            if (SeqIdWrite (sip, str, PRINTID_TEXTID_ACC_VER, sizeof (str)) != NULL) {
              sprintf (tmp, "\t\t\tprotein_id\tlcl|%s\n", str);
              ValNodeCopyStr (head, 0, tmp);
            }
          } else if (sip->choice == SEQID_GENERAL) {
            if (SeqIdWrite (sip, str, PRINTID_FASTA_GENERAL, sizeof (str)) != NULL) {
              sprintf (tmp, "\t\t\tprotein_id\t%s\n", str);
              ValNodeCopyStr (head, 0, tmp);
            }
          }
        }
      }
      break;
    case SEQFEAT_RNA :
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp != NULL) {
        switch (rrp->ext.choice) {
          case 1 :
            StringNCpy_0 (str, (CharPtr) rrp->ext.value.ptrvalue, sizeof (str));
            if (! StringHasNoText (str)) {
              sprintf (tmp, "\t\t\tproduct\t%s\n", str);
              ValNodeCopyStr (head, 0, tmp);
            }
            break;
          case 2 :
            trp = rrp->ext.value.ptrvalue;
            if (trp != NULL) {
              FeatDefLabel (sfp, str, sizeof (str) - 1, OM_LABEL_CONTENT);
              if (! StringHasNoText (str)) {
                sprintf (tmp, "\t\t\tproduct\t%s\n", str);
                ValNodeCopyStr (head, 0, tmp);
              }
            }
            break;
          default :
            break;
        }
      }
      break;
    default :
      break;
  }
  if (pseudo) {
    ValNodeCopyStr (head, 0, "\t\t\tpseudo\n");
  }
  grp = SeqMgrGetGeneXref (sfp);
  if (grp != NULL && SeqMgrGeneIsSuppressed (grp)) {
    ValNodeCopyStr (head, 0, "\t\t\tgene\t-\n");
  }
  if (! StringHasNoText (sfp->comment)) {
    ValNodeCopyStr (head, 0, "\t\t\tnote\t");
    ValNodeCopyStr (head, 0, sfp->comment);
    ValNodeCopyStr (head, 0, "\n");
  }
  switch (sfp->exp_ev) {
    case 1 :
      ValNodeCopyStr (head, 0, "\t\t\tevidence\texperimental\n");
      break;
    case 2 :
      ValNodeCopyStr (head, 0, "\t\t\tevidence\tnot_experimental\n");
      break;
    default :
      break;
  }
  if (! StringHasNoText (sfp->except_text)) {
    ValNodeCopyStr (head, 0, "\t\t\texception\t");
    ValNodeCopyStr (head, 0, sfp->except_text);
    ValNodeCopyStr (head, 0, "\n");
  } else if (sfp->excpt) {
    ValNodeCopyStr (head, 0, "\t\t\texception\n");
  }
  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    if (! StringHasNoText (gbq->qual)) {
      if (! StringHasNoText (gbq->val)) {
        sprintf (tmp, "\t\t\t%s\t%s\n", gbq->qual, gbq->val);
        ValNodeCopyStr (head, 0, tmp);
      }
    }
  }
  VisitUserObjectsInUop (sfp->ext, (Pointer) head, PrintFTUserObj);
  for (vnp = geneorprotdb; vnp != NULL; vnp = vnp->next) {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt != NULL) {
      if (! StringHasNoText (dbt->db)) {
        oip = dbt->tag;
        if (oip->str != NULL && (! StringHasNoText (oip->str))) {
          sprintf (tmp, "\t\t\tdb_xref\t%s:%s\n", dbt->db, oip->str);
          ValNodeCopyStr (head, 0, tmp);
        } else {
          sprintf (tmp, "\t\t\tdb_xref\t%s:%ld\n", dbt->db, (long) oip->id);
          ValNodeCopyStr (head, 0, tmp);
        }
      }
    }
  }
  for (vnp = sfp->dbxref; vnp != NULL; vnp = vnp->next) {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt != NULL) {
      if (! StringHasNoText (dbt->db)) {
        oip = dbt->tag;
        if (oip->str != NULL && (! StringHasNoText (oip->str))) {
          sprintf (tmp, "\t\t\tdb_xref\t%s:%s\n", dbt->db, oip->str);
          ValNodeCopyStr (head, 0, tmp);
        } else {
          sprintf (tmp, "\t\t\tdb_xref\t%s:%ld\n", dbt->db, (long) oip->id);
          ValNodeCopyStr (head, 0, tmp);
        }
      }
    }
  }
}

static BioseqPtr FindFirstBioseq (SeqEntryPtr sep)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;

  if (sep == NULL || sep->data.ptrvalue == NULL ||
      /* sep->choice < 0 || */ sep->choice > 2) return NULL;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    return bsp;
  }
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
    bsp = FindFirstBioseq (sep);
    if (bsp != NULL) return bsp;
  }
  return NULL;
}

static BioseqPtr BioseqLockAndIndexByEntity (Uint2 entityID)

{
  BioseqPtr    bsp;
  SeqEntryPtr  sep;
  SeqIdPtr     sip;

  if (entityID < 1) return NULL;

  sep = SeqMgrGetSeqEntryForEntityID (entityID);
  if (sep == NULL) return NULL;

  bsp = FindFirstBioseq (sep);
  if (bsp == NULL) return NULL;

  sip = SeqIdFindBest (bsp->id, 0);
  if (sip == NULL) return NULL;

  bsp = BioseqLockById (sip);
  if (bsp == NULL) return NULL;

  if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
    SeqMgrIndexFeatures (entityID, NULL);
  }

  return bsp;
}

NLM_EXTERN CharPtr asn2gnbk_format (
  Asn2gbJobPtr ajp,
  Int4 paragraph
)

{
  Asn2gbFormat       af;
  Asn2gbSectPtr      asp;
  BaseBlockPtr       bbp;
  BlockType          blocktype;
  BioseqPtr          bsp;
  SeqMgrFeatContext  fcontext;
  FormatProc         fmt;
  ValNodePtr         head;
  IntAsn2gbJobPtr    iajp;
  Char               id [42];
  size_t             max;
  SeqEntryPtr        oldscope;
  QualValPtr         qv;
  Int4               section;
  SeqEntryPtr        sep;
  SeqFeatPtr         sfp;
  SeqIdPtr           sip;
  CharPtr            str = NULL;
  BioseqPtr          target;
  Char               tmp [53];

  /* qv must hold MAX (ASN2GNBK_TOTAL_SOURCE, ASN2GNBK_TOTAL_FEATUR) */

  iajp = (IntAsn2gbJobPtr) ajp;
  if (iajp == NULL || ajp->sectionArray == NULL || ajp->paragraphArray == NULL) return NULL;
  if (paragraph < 0 || paragraph >= ajp->numParagraphs) return NULL;

  bbp = ajp->paragraphArray [paragraph];
  if (bbp == NULL) return NULL;

  section = bbp->section;
  if (section < 0 || section >= ajp->numSections) return NULL;

  asp = ajp->sectionArray [section];
  if (asp == NULL) return NULL;

  blocktype = bbp->blocktype;
  if (blocktype < LOCUS_BLOCK || blocktype > SLASH_BLOCK) return NULL;

  max = (size_t) (MAX (ASN2GNBK_TOTAL_SOURCE, ASN2GNBK_TOTAL_FEATUR));
  qv = MemNew (sizeof (QualVal) * (max + 5));
  if (qv == NULL) return NULL;

  af.ajp = (IntAsn2gbJobPtr) ajp;
  af.asp = asp;
  af.qvp = qv;
  af.format = iajp->format;

  sep = GetTopSeqEntryForEntityID (bbp->entityID);

  if (iajp->format != FTABLE_FMT) {
    fmt = asn2gnbk_fmt_functions [(int) blocktype];
    if (fmt == NULL) return NULL;

    bsp = BioseqLockAndIndexByEntity (bbp->entityID);
    oldscope = SeqEntrySetScope (sep);

    str = fmt (&af, bbp);

    SeqEntrySetScope (oldscope);
    BioseqUnlock (bsp);

  } else {

    target = asp->target;
    if (target != NULL) {

      bsp = BioseqLockAndIndexByEntity (bbp->entityID);
      oldscope = SeqEntrySetScope (sep);

      if (blocktype == FEATHEADER_BLOCK) {
        sip = SeqIdFindBest (target->id, 0);
        SeqIdWrite (sip, id, PRINTID_FASTA_LONG, sizeof (id) - 1);
        if (! StringHasNoText (id)) {
          sprintf (tmp, ">Feature %s\n", id);
          str = StringSave (tmp);
        }

      } else if (blocktype == FEATURE_BLOCK) {

        sfp = SeqMgrGetDesiredFeature (bbp->entityID, NULL, bbp->itemID, 0, NULL, &fcontext);
        if (sfp != NULL) {
          head = NULL;
          PrintFtableLocAndQuals (af.ajp, &head, target, sfp, &fcontext);
          str = MergeValNodeStrings (head);
          ValNodeFreeData (head);
        }
      }

      SeqEntrySetScope (oldscope);
      BioseqUnlock (bsp);
    }
  }

  if (str == NULL) {
    str = StringSave ("???\n");
  }

  MemFree (qv);

  return str;
}

NLM_EXTERN Asn2gbJobPtr asn2gnbk_cleanup (
  Asn2gbJobPtr ajp
)

{
  Asn2gbSectPtr     asp;
  BaseBlockPtr      bbp;
  BaseBlockPtr      PNTR blockArray;
  Int4              i;
  IntAsn2gbJobPtr   iajp;
  IntAsn2gbSectPtr  iasp;
  IntCdsBlockPtr    icp;
  IntFeatBlockPtr   ifp;
  IntRefBlockPtr    irp;
  IntSrcBlockPtr    isp;
  Int4              j;
  Int4              numBlocks;
  Int4              numSections;
  RefBlockPtr       rrp;
  Asn2gbSectPtr     PNTR sectionArray;
  StringItemPtr     sip, nxt;

  iajp = (IntAsn2gbJobPtr) ajp;
  if (iajp == NULL) return NULL;

  SeqLocFree (iajp->ajp.slp);

  numSections = ajp->numSections;
  sectionArray = ajp->sectionArray;

  if (sectionArray != NULL) {

    for (i = 0; i < numSections; i++) {
      asp = sectionArray [i];
      if (asp != NULL) {
        iasp = (IntAsn2gbSectPtr) asp;

        numBlocks = asp->numBlocks;
        blockArray = asp->blockArray;
        if (blockArray != NULL) {

          for (j = 0; j < numBlocks; j++) {
            bbp = blockArray [j];
            if (bbp != NULL) {

              MemFree (bbp->string);

              if (bbp->blocktype == REFERENCE_BLOCK) {
                rrp = (RefBlockPtr) bbp;
                MemFree (rrp->uniquestr);
                irp = (IntRefBlockPtr) rrp;
                DateFree (irp->date);
                SeqLocFree (irp->loc);
                MemFree (irp->authstr);
                MemFree (irp->fig);
                MemFree (irp->maploc);

              } else if (bbp->blocktype == SOURCEFEAT_BLOCK) {

                isp = (IntSrcBlockPtr) bbp;
                SeqLocFree (isp->loc);

              } else if (bbp->blocktype == FEATURE_BLOCK) {

                ifp = (IntFeatBlockPtr) bbp;
                if (ifp->isCDS) {
                  icp = (IntCdsBlockPtr) ifp;
                  MemFree (icp->fig);
                  MemFree (icp->maploc);
                }
              }

              MemFree (bbp);
            }
          }
        }
        MemFree (asp->blockArray);
        MemFree (asp->referenceArray);
        SeqPortFree (iasp->spp);
        MemFree (asp);
      }
    }
  }

  MemFree (ajp->sectionArray);
  MemFree (ajp->paragraphArray);
  MemFree (ajp->paragraphByIDs);

  sip = iajp->pool;
  while (sip != NULL) {
    nxt = sip->next;
    MemFree (sip);
    sip = nxt;
  }

  if (iajp->lockedBspList != NULL) {
    UnlockFarComponents (iajp->lockedBspList);
  }

  free_buff ();
  FiniWWW (iajp);

  MemFree (iajp);

  return NULL;
}

static CharPtr defHead = "\
Content-type: text/html\n\n\
<HTML>\n\
<HEAD><TITLE>GenBank entry</TITLE></HEAD>\n\
<BODY>\n\
<hr>\n\
<pre>";

static CharPtr defTail = "\
</pre>\n\
<hr>\n\
</BODY>\n\
</HTML>\n";

NLM_EXTERN void AsnPrintNewLine PROTO((AsnIoPtr aip));

NLM_EXTERN Boolean SeqEntryToGnbk (
  SeqEntryPtr sep,
  SeqLocPtr slp,
  FmtType format,
  ModType mode,
  StlType style,
  FlgType flags,
  LckType locks,
  XtraPtr extra,
  FILE *fp
)

{
  AsnIoPtr           aip = NULL;
  Asn2gbJobPtr       ajp;
  AsnTypePtr         atp = NULL;
  BaseBlockPtr       bbp;
  BlockType          block;
  BioseqPtr          bsp = NULL;
  BioseqSetPtr       bssp = NULL;
  Boolean            do_gbseq_asn = FALSE;
  Boolean            do_gbseq_xml = FALSE;
  CharPtr            ffhead = NULL;
  CharPtr            fftail = NULL;
  Asn2gbWriteFunc    ffwrite = NULL;
  GBSeqPtr           gbseq = NULL;
  GBSeq              gbsq;
  GBSeqPtr           gbtmp;
  Int4               i;
  IntAsn2gbJobPtr    iajp;
  IndxPtr            index = NULL;
  Boolean            is_html;
  Int4               numParagraphs;
  BaseBlockPtr PNTR  paragraphArray;
  Boolean            rsult = FALSE;
  CharPtr            str;
  Int1               type = ASNIO_TEXT_OUT;
  Pointer            userdata = NULL;
  XtraBlock          xtra;
#ifdef WIN_MAC
#if __profile__
  ValNodePtr         bsplist = NULL;
  Uint2              entityID;
  Boolean            lockFarComp;
  Boolean            lockFarLocs;
  Boolean            lockFarProd;
  Boolean            lookupFarComp;
  Boolean            lookupFarHist;
  Boolean            lookupFarLocs;
  Boolean            lookupFarProd;
#endif
#endif

  if (extra != NULL) {
    ffwrite = extra->ffwrite;
    ffhead = extra->ffhead;
    fftail = extra->fftail;
    index = extra->index;
    gbseq = extra->gbseq;
    aip = extra->aip;
    atp = extra->atp;
    userdata = extra->userdata;
  }
  if (fp == NULL && ffwrite == NULL && aip == NULL) return FALSE;
  if (sep == NULL && slp == NULL) return FALSE;
  if (sep != NULL) {
    if (IS_Bioseq (sep)) {
      bsp = (BioseqPtr) sep->data.ptrvalue;
    } else if (IS_Bioseq_set (sep)) {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
    }
  }

#ifdef WIN_MAC
#if __profile__
  /* this allows profiling of just the formatter, without feature indexing, on the Mac */

  if (sep != NULL) {
    entityID = ObjMgrGetEntityIDForPointer (sep->data.ptrvalue);
    if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
      SeqMgrIndexFeatures (entityID, NULL);
    }
  }

  lockFarComp = (Boolean) ((locks & LOCK_FAR_COMPONENTS) != 0);
  lockFarLocs = (Boolean) ((locks & LOCK_FAR_LOCATIONS) != 0);
  lockFarProd = (Boolean) ((locks & LOCK_FAR_PRODUCTS) != 0);

  if (lockFarComp || lockFarLocs || lockFarProd) {
    locks = locks ^ (LOCK_FAR_COMPONENTS | LOCK_FAR_LOCATIONS | LOCK_FAR_PRODUCTS);
    bsplist = LockFarComponentsEx (sep, lockFarComp, lockFarLocs, lockFarProd);
  }

  lookupFarComp = (Boolean) ((locks & LOOKUP_FAR_COMPONENTS) != 0);
  lookupFarLocs = (Boolean) ((locks & LOOKUP_FAR_LOCATIONS) != 0);
  lookupFarProd = (Boolean) ((locks & LOOKUP_FAR_PRODUCTS) != 0);
  lookupFarHist = (Boolean) ((locks & LOOKUP_FAR_HISTORY) != 0);

  if (lookupFarComp || lookupFarLocs || lookupFarProd || lookupFarHist) {
    locks = locks ^ (LOOKUP_FAR_COMPONENTS | LOOKUP_FAR_LOCATIONS | LOOKUP_FAR_PRODUCTS | LOOKUP_FAR_HISTORY);
    LookupFarSeqIDs (sep, lookupFarComp, lookupFarLocs, lookupFarProd, FALSE, lookupFarHist);
  }

  ProfilerSetStatus (TRUE);
#endif
#endif

  do_gbseq_xml = (Boolean) ((flags & HTML_XML_ASN_MASK) == CREATE_XML_GBSEQ_FILE);
  do_gbseq_asn = (Boolean) ((flags & HTML_XML_ASN_MASK) == CREATE_ASN_GBSEQ_FILE);

  if (do_gbseq_xml || do_gbseq_asn) {
    if (fp != NULL && aip == NULL) {
      if (do_gbseq_xml) {
        type |= ASNIO_XML;
      }
      aip = AsnIoNew (type, fp, NULL, NULL, NULL);
      fp = NULL;
    }
    if (extra == NULL) {
      MemSet ((Pointer) &xtra, 0, sizeof (XtraBlock));
      extra = &xtra;
    }
    if (extra->gbseq == NULL) {
      MemSet ((Pointer) &gbsq, 0, sizeof (GBSeq));
      extra->gbseq = &gbsq;
      gbseq = extra->gbseq;
    }
  }

  ajp = asn2gnbk_setup (bsp, bssp, slp, format, mode, style, flags, locks, extra);


  if (ajp != NULL) {
    rsult = TRUE;
    iajp = (IntAsn2gbJobPtr) ajp;

    /* send optional head string */
    is_html = (Boolean) ((flags & HTML_XML_ASN_MASK) == CREATE_HTML_FLATFILE);
    if (ffhead == NULL && is_html) {
      ffhead = defHead;
    }
    if (ffhead != NULL) {
      if (fp != NULL) {
        fprintf (fp, ffhead);
      }
    }
    if (ffwrite != NULL) {
      ffwrite (ffhead, userdata, HEAD_BLOCK);
    }

    /* send each paragraph */

    numParagraphs = ajp->numParagraphs;
    paragraphArray = ajp->paragraphArray;

    for (i = 0; i < numParagraphs; i++) {
      str = asn2gnbk_format (ajp, i);
      block = 0;
      if (paragraphArray != NULL) {
        bbp = paragraphArray [i];
        if (bbp != NULL) {
          block = bbp->blocktype;
        }
      }
      if (str != NULL) {
        if (fp != NULL) {
          fprintf (fp, "%s", str);
        }
        if (ffwrite != NULL) {
          ffwrite (str, userdata, block);
        }
      } else {
        if (fp != NULL) {
          fprintf (fp, "?\n");
        }
        if (ffwrite != NULL) {
          ffwrite ("?\n", userdata, block);
        }
      }

      /* if generating GBSeq XML/ASN, write at each slash block */

      if (block == SLASH_BLOCK && gbseq != NULL && aip != NULL) {
        GBSeqAsnWrite (gbseq, aip, atp);
        if (atp == NULL) {
          AsnPrintNewLine (aip);
        }
        AsnIoFlush (aip);

        /* clean up gbseq fields */

        gbtmp = GBSeqNew ();
        MemCopy (gbtmp, gbseq, sizeof (GBSeq));
        MemSet (gbseq, 0, sizeof (GBSeq));
        GBSeqFree (gbtmp);
      }
      MemFree (str);
    }

    /* if RELEASE_MODE, warn if unresolved gi numbers, missing translation, etc. */
    
    if (iajp->relModeError && mode == RELEASE_MODE) {
      rsult = FALSE;
    }

    /* send optional tail string */

    if (fftail == NULL && is_html) {
      fftail = defTail;
    }
    if (fftail != NULL) {
      if (fp != NULL) {
        fprintf (fp, fftail);
      }
    }
    if (ffwrite != NULL) {
      ffwrite (fftail, userdata, TAIL_BLOCK);
    }

    asn2gnbk_cleanup (ajp);
  }


#ifdef WIN_MAC
#if __profile__
  ProfilerSetStatus (FALSE);

  UnlockFarComponents (bsplist);
#endif
#endif

  return rsult;
}

NLM_EXTERN Boolean BioseqToGnbk (
  BioseqPtr bsp,
  SeqLocPtr slp,
  FmtType format,
  ModType mode,
  StlType style,
  FlgType flags,
  LckType locks,
  XtraPtr extra,
  FILE *fp
)

{
  SeqEntryPtr  sep = NULL;

  if (bsp == NULL && slp == NULL) return FALSE;
  if (bsp != NULL) {
    sep = SeqMgrGetSeqEntryForData (bsp);
  }
  return SeqEntryToGnbk (sep, slp, format, mode, style, flags, locks, extra, fp);
}

