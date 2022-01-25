#ifndef _objmime_ 
#define _objmime_ 

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
*    Generated objects for Module NCBI-Mime
*    Generated using ASNCODE Revision: 6.1 at Mar 23, 1998  6:44 PM
*
**************************************************/

NLM_EXTERN Boolean LIBCALL
objmimeAsnLoad PROTO((void));
typedef ValNodePtr NcbiMimeAsn1Ptr;
typedef ValNode NcbiMimeAsn1;
#define NcbiMimeAsn1_entrez 1
#define NcbiMimeAsn1_alignstruc 2
#define NcbiMimeAsn1_alignseq 3
#define NcbiMimeAsn1_strucseq 4


NLM_EXTERN NcbiMimeAsn1Ptr LIBCALL NcbiMimeAsn1Free PROTO ((NcbiMimeAsn1Ptr ));
NLM_EXTERN NcbiMimeAsn1Ptr LIBCALL NcbiMimeAsn1AsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL NcbiMimeAsn1AsnWrite PROTO (( NcbiMimeAsn1Ptr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    EntrezGeneral
*
**************************************************/
typedef struct struct_Entrez_general {
   CharPtr   title;
   ValNodePtr   Data_data;
   Uint2   style;
   CharPtr   location;
} EntrezGeneral, PNTR EntrezGeneralPtr;


NLM_EXTERN EntrezGeneralPtr LIBCALL EntrezGeneralFree PROTO ((EntrezGeneralPtr ));
NLM_EXTERN EntrezGeneralPtr LIBCALL EntrezGeneralNew PROTO (( void ));
NLM_EXTERN EntrezGeneralPtr LIBCALL EntrezGeneralAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL EntrezGeneralAsnWrite PROTO (( EntrezGeneralPtr , AsnIoPtr, AsnTypePtr));


#ifdef NLM_GENERATED_CODE_PROTO

typedef ValNodePtr Data_dataPtr;
typedef ValNode Data_data;

#endif /* NLM_GENERATED_CODE_PROTO */

#define Data_data_ml 1
#define Data_data_prot 2
#define Data_data_nuc 3
#define Data_data_genome 4
#define Data_data_structure 5
#define Data_data_strucAnnot 6

#ifdef NLM_GENERATED_CODE_PROTO

static Data_dataPtr LIBCALL Data_dataFree PROTO ((Data_dataPtr ));
static Data_dataPtr LIBCALL Data_dataAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
static Boolean LIBCALL Data_dataAsnWrite PROTO (( Data_dataPtr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */



/**************************************************
*
*    BiostrucAlign
*
**************************************************/
typedef struct struct_Biostruc_align {
   struct struct_Biostruc PNTR   master;
   struct struct_Biostruc PNTR   slaves;
   struct struct_Biostruc_annot_set PNTR   alignments;
   ValNodePtr   sequences;
   struct seqannot PNTR   seqalign;   /* hand change -- lyg */
} BiostrucAlign, PNTR BiostrucAlignPtr;


NLM_EXTERN BiostrucAlignPtr LIBCALL BiostrucAlignFree PROTO ((BiostrucAlignPtr ));
NLM_EXTERN BiostrucAlignPtr LIBCALL BiostrucAlignNew PROTO (( void ));
NLM_EXTERN BiostrucAlignPtr LIBCALL BiostrucAlignAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL BiostrucAlignAsnWrite PROTO (( BiostrucAlignPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    BiostrucAlignSeq
*
**************************************************/
typedef struct struct_Biostruc_align_seq {
   ValNodePtr   sequences;
  struct seqannot PNTR   seqalign;   /* hand change struct -- lyg */
} BiostrucAlignSeq, PNTR BiostrucAlignSeqPtr;


NLM_EXTERN BiostrucAlignSeqPtr LIBCALL BiostrucAlignSeqFree PROTO ((BiostrucAlignSeqPtr ));
NLM_EXTERN BiostrucAlignSeqPtr LIBCALL BiostrucAlignSeqNew PROTO (( void ));
NLM_EXTERN BiostrucAlignSeqPtr LIBCALL BiostrucAlignSeqAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL BiostrucAlignSeqAsnWrite PROTO (( BiostrucAlignSeqPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    BiostrucSeq
*
**************************************************/
typedef struct struct_Biostruc_seq {
   struct struct_Biostruc PNTR  structure; 
   ValNodePtr   sequences;
} BiostrucSeq, PNTR BiostrucSeqPtr;


NLM_EXTERN BiostrucSeqPtr LIBCALL BiostrucSeqFree PROTO ((BiostrucSeqPtr ));
NLM_EXTERN BiostrucSeqPtr LIBCALL BiostrucSeqNew PROTO (( void ));
NLM_EXTERN BiostrucSeqPtr LIBCALL BiostrucSeqAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL BiostrucSeqAsnWrite PROTO (( BiostrucSeqPtr , AsnIoPtr, AsnTypePtr));

/* following #defines are for enumerated type, not used by object loaders */
#define Entrez_style_docsum 1
#define Entrez_style_genbank 2
#define Entrez_style_genpept 3
#define Entrez_style_fasta 4
#define Entrez_style_asn1 5
#define Entrez_style_graphic 6
#define Entrez_style_alignment 7
#define Entrez_style_globalview 8
#define Entrez_style_report 9
#define Entrez_style_medlars 10
#define Entrez_style_embl 11
#define Entrez_style_pdb 12
#define Entrez_style_kinemage 13

#ifdef __cplusplus
/* { */ }
#endif

#endif /* _objmime_ */

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

