#ifndef _objvalid_ 
#define _objvalid_ 

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
*    Generated objects for Module NCBI-Structured-comment-validation
*    Generated using ASNCODE Revision: 6.16 at Apr 7, 2009  1:38 PM
*
**************************************************/

NLM_EXTERN Boolean LIBCALL
objvalidAsnLoad PROTO((void));


/**************************************************
*
*    FieldRule
*
**************************************************/
typedef struct struct_Field_rule {
   struct struct_Field_rule PNTR next;
   CharPtr   field_name;
   CharPtr   match_expression;
   Uint1   required;
} FieldRule, PNTR FieldRulePtr;


NLM_EXTERN FieldRulePtr LIBCALL FieldRuleFree PROTO ((FieldRulePtr ));
NLM_EXTERN FieldRulePtr LIBCALL FieldRuleNew PROTO (( void ));
NLM_EXTERN FieldRulePtr LIBCALL FieldRuleAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL FieldRuleAsnWrite PROTO (( FieldRulePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    FieldSet
*
**************************************************/
typedef struct struct_Field_rule FieldSet;
typedef struct struct_Field_rule PNTR FieldSetPtr;
#define FieldSetNew() Field_ruleNew() 

#ifdef NLM_GENERATED_CODE_PROTO

NLM_EXTERN FieldSetPtr LIBCALL FieldSetFree PROTO ((FieldSetPtr ));
NLM_EXTERN FieldSetPtr LIBCALL FieldSetNew PROTO (( void ));
NLM_EXTERN FieldSetPtr LIBCALL FieldSetAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL FieldSetAsnWrite PROTO (( FieldSetPtr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */



/**************************************************
*
*    CommentRule
*
**************************************************/
typedef struct struct_Comment_rule {
   struct struct_Comment_rule PNTR next;
   CharPtr   prefix;
   Uint1   updated;
   struct struct_Field_rule PNTR   fields;
} CommentRule, PNTR CommentRulePtr;


NLM_EXTERN CommentRulePtr LIBCALL CommentRuleFree PROTO ((CommentRulePtr ));
NLM_EXTERN CommentRulePtr LIBCALL CommentRuleNew PROTO (( void ));
NLM_EXTERN CommentRulePtr LIBCALL CommentRuleAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL CommentRuleAsnWrite PROTO (( CommentRulePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    CommentSet
*
**************************************************/
typedef struct struct_Comment_rule CommentSet;
typedef struct struct_Comment_rule PNTR CommentSetPtr;
#define CommentSetNew() Comment_ruleNew() 

#ifdef NLM_GENERATED_CODE_PROTO

NLM_EXTERN CommentSetPtr LIBCALL CommentSetFree PROTO ((CommentSetPtr ));
NLM_EXTERN CommentSetPtr LIBCALL CommentSetNew PROTO (( void ));
NLM_EXTERN CommentSetPtr LIBCALL CommentSetAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
NLM_EXTERN Boolean LIBCALL CommentSetAsnWrite PROTO (( CommentSetPtr , AsnIoPtr, AsnTypePtr));

#endif /* NLM_GENERATED_CODE_PROTO */

#ifdef __cplusplus
/* { */ }
#endif

#endif /* _objvalid_ */

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

