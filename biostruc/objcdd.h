#ifndef _objcdd_ 
#define _objcdd_ 

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
*    Generated objects for Module NCBI-Cdd
*    Generated using ASNCODE Revision: 6.8 at Oct 19, 1999  1:13 PM
*
**************************************************/

NLM_EXTERN Boolean LIBCALL
objcddAsnLoad PROTO((void));
typedef ValNodePtr CddIdPtr;
typedef ValNode CddId;
#define CddId_uid 1
#define CddId_gid 2


NLM_EXTERN CddIdPtr LIBCALL CddIdFree PROTO ((CddIdPtr ));
NLM_EXTERN CddIdPtr LIBCALL CddIdAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL CddIdAsnWrite PROTO (( CddIdPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    CddIdSet
*
**************************************************/
typedef ValNode CddIdSet;
typedef ValNodePtr CddIdSetPtr;
#define CddIdSetNew() ValNodeNew(NULL) 

#ifdef NLM_GENERATED_CODE_PROTO

NLM_EXTERN CddIdSetPtr LIBCALL CddIdSetFree PROTO ((CddIdSetPtr ));
NLM_EXTERN CddIdSetPtr LIBCALL CddIdSetNew PROTO (( void ));
NLM_EXTERN CddIdSetPtr LIBCALL CddIdSetAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL CddIdSetAsnWrite PROTO (( CddIdSetPtr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */



/**************************************************
*
*    Cdd
*
**************************************************/
typedef struct struct_Cdd {
   struct struct_Cdd PNTR next;
   CharPtr   name;
   ValNodePtr   othernames;
   ValNodePtr   categories;
   ValNodePtr   id;
   ValNodePtr   description;
/*   ValNodePtr   create_date;*/
   DatePtr   create_date;
   CharPtr   source;
/*   struct struct_Seq_annot PNTR   seqannot;*/
   struct seqannot PNTR   seqannot;
   struct struct_Biostruc_feature_set PNTR   strucalign;
} Cdd, PNTR CddPtr;


NLM_EXTERN CddPtr LIBCALL CddFree PROTO ((CddPtr ));
NLM_EXTERN CddPtr LIBCALL CddNew PROTO (( void ));
NLM_EXTERN CddPtr LIBCALL CddAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL CddAsnWrite PROTO (( CddPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    CddSet
*
**************************************************/
typedef struct struct_Cdd CddSet;
typedef struct struct_Cdd PNTR CddSetPtr;
#define CddSetNew() CddNew() 

#ifdef NLM_GENERATED_CODE_PROTO

NLM_EXTERN CddSetPtr LIBCALL CddSetFree PROTO ((CddSetPtr ));
NLM_EXTERN CddSetPtr LIBCALL CddSetNew PROTO (( void ));
NLM_EXTERN CddSetPtr LIBCALL CddSetAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL CddSetAsnWrite PROTO (( CddSetPtr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */



/**************************************************
*
*    CddTree
*
**************************************************/
typedef struct struct_Cdd_tree {
   struct struct_Cdd_tree PNTR next;
   ValNodePtr   id;
   ValNodePtr   parents;
   ValNodePtr   children;
   CharPtr   name;
} CddTree, PNTR CddTreePtr;


NLM_EXTERN CddTreePtr LIBCALL CddTreeFree PROTO ((CddTreePtr ));
NLM_EXTERN CddTreePtr LIBCALL CddTreeNew PROTO (( void ));
NLM_EXTERN CddTreePtr LIBCALL CddTreeAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL CddTreeAsnWrite PROTO (( CddTreePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    CddTreeSet
*
**************************************************/
typedef struct struct_CddTree CddTreeSet;
typedef struct struct_CddTree PNTR CddTreeSetPtr;
#define CddTreeSetNew() CddTreeNew() 

#ifdef NLM_GENERATED_CODE_PROTO

NLM_EXTERN CddTreeSetPtr LIBCALL CddTreeSetFree PROTO ((CddTreeSetPtr ));
NLM_EXTERN CddTreeSetPtr LIBCALL CddTreeSetNew PROTO (( void ));
NLM_EXTERN CddTreeSetPtr LIBCALL CddTreeSetAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL CddTreeSetAsnWrite PROTO (( CddTreeSetPtr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */



/**************************************************
*
*    GlobalId
*
**************************************************/
typedef struct struct_Global_id {
   CharPtr   accession;
   CharPtr   release;
   Int4   version;
   CharPtr   database;
} GlobalId, PNTR GlobalIdPtr;


NLM_EXTERN GlobalIdPtr LIBCALL GlobalIdFree PROTO ((GlobalIdPtr ));
NLM_EXTERN GlobalIdPtr LIBCALL GlobalIdNew PROTO (( void ));
NLM_EXTERN GlobalIdPtr LIBCALL GlobalIdAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL GlobalIdAsnWrite PROTO (( GlobalIdPtr , AsnIoPtr, AsnTypePtr));

typedef ValNodePtr CddDescrPtr;
typedef ValNode CddDescr;
#define CddDescr_comment 1
#define CddDescr_reference 2


NLM_EXTERN CddDescrPtr LIBCALL CddDescrFree PROTO ((CddDescrPtr ));
NLM_EXTERN CddDescrPtr LIBCALL CddDescrAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL CddDescrAsnWrite PROTO (( CddDescrPtr , AsnIoPtr, AsnTypePtr));

#ifdef __cplusplus
/* { */ }
#endif

#endif /* _objcdd_ */

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

