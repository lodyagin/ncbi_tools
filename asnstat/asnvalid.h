/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "asnvalid.h11";
static AsnValxNode avnx[2] = {
    {2,NULL,0,0.0,NULL } ,
    {2,NULL,0,0.0,NULL } };

static AsnType atx[16] = {
  {401, "Field-rule" ,1,0,0,0,0,1,0,0,NULL,&atx[6],&atx[1],0,&atx[7]} ,
  {0, "field-name" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[3]} ,
  {323, "VisibleString" ,0,26,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "match-expression" ,128,1,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[4]} ,
  {0, "required" ,128,2,0,0,1,0,0,0,&avnx[0],&atx[5],NULL,0,NULL} ,
  {301, "BOOLEAN" ,0,1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {311, "SEQUENCE" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {402, "Field-set" ,1,0,0,0,0,1,0,0,NULL,&atx[9],&atx[8],0,&atx[10]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[0],NULL,0,NULL} ,
  {314, "SET OF" ,0,17,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {403, "Comment-rule" ,1,0,0,0,0,1,0,0,NULL,&atx[6],&atx[11],0,&atx[14]} ,
  {0, "prefix" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[12]} ,
  {0, "updated" ,128,1,0,0,1,0,0,0,&avnx[1],&atx[5],NULL,0,&atx[13]} ,
  {0, "fields" ,128,2,0,0,0,0,0,0,NULL,&atx[7],NULL,0,NULL} ,
  {404, "Comment-set" ,1,0,0,0,0,1,0,0,NULL,&atx[9],&atx[15],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[10],NULL,0,NULL} };

static AsnModule ampx[1] = {
  { "NCBI-Structured-comment-validation" , "asnvalid.h11",&atx[0],NULL,NULL,0,0} };

static AsnValxNodePtr avn = avnx;
static AsnTypePtr at = atx;
static AsnModulePtr amp = ampx;



/**************************************************
*
*    Defines for Module NCBI-Structured-comment-validation
*
**************************************************/

#define FIELD_RULE &at[0]
#define FIELD_RULE_field_name &at[1]
#define FIELD_RULE_match_expression &at[3]
#define FIELD_RULE_required &at[4]

#define FIELD_SET &at[7]
#define FIELD_SET_E &at[8]

#define COMMENT_RULE &at[10]
#define COMMENT_RULE_prefix &at[11]
#define COMMENT_RULE_updated &at[12]
#define COMMENT_RULE_fields &at[13]

#define COMMENT_SET &at[14]
#define COMMENT_SET_E &at[15]
