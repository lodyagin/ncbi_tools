/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "cdd.h";
static AsnValxNode avnx[6] = {
    {20,"unassigned" ,0,0.0,&avnx[1] } ,
    {20,"finished-ok" ,1,0.0,&avnx[2] } ,
    {20,"pending-release" ,2,0.0,&avnx[3] } ,
    {20,"other-asis" ,3,0.0,&avnx[4] } ,
    {20,"matrix-only" ,4,0.0,&avnx[5] } ,
    {20,"other" ,255,0.0,NULL } };

static AsnType atx[71] = {
  {401, "Cdd-id" ,1,0,0,0,0,1,0,0,NULL,&atx[11],&atx[1],0,&atx[12]} ,
  {0, "uid" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[3]} ,
  {302, "INTEGER" ,0,2,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "gid" ,128,1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {416, "Global-id" ,1,0,0,0,0,0,0,0,NULL,&atx[10],&atx[5],0,&atx[21]} ,
  {0, "accession" ,128,0,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[7]} ,
  {323, "VisibleString" ,0,26,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "release" ,128,1,0,1,0,0,0,0,NULL,&atx[6],NULL,0,&atx[8]} ,
  {0, "version" ,128,2,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[9]} ,
  {0, "database" ,128,3,0,1,0,0,0,0,NULL,&atx[6],NULL,0,NULL} ,
  {311, "SEQUENCE" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {315, "CHOICE" ,0,-1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {402, "Cdd-id-set" ,1,0,0,0,0,1,0,0,NULL,&atx[14],&atx[13],0,&atx[15]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[0],NULL,0,NULL} ,
  {312, "SEQUENCE OF" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {403, "Cdd" ,1,0,0,0,0,1,0,0,NULL,&atx[10],&atx[16],0,&atx[60]} ,
  {0, "name" ,128,0,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[17]} ,
  {0, "id" ,128,1,0,0,0,0,0,0,NULL,&atx[12],NULL,0,&atx[18]} ,
  {0, "description" ,128,2,0,1,0,0,0,0,NULL,&atx[19],NULL,0,&atx[34]} ,
  {418, "Cdd-descr-set" ,1,0,0,0,0,0,0,0,NULL,&atx[33],&atx[20],0,&atx[46]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[21],NULL,0,NULL} ,
  {417, "Cdd-descr" ,1,0,0,0,0,0,0,0,NULL,&atx[11],&atx[22],0,&atx[19]} ,
  {0, "othername" ,128,0,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[23]} ,
  {0, "category" ,128,1,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[24]} ,
  {0, "comment" ,128,2,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[25]} ,
  {0, "reference" ,128,3,0,0,0,0,0,0,NULL,&atx[26],NULL,0,&atx[27]} ,
  {408, "Pub" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[38]} ,
  {0, "create-date" ,128,4,0,0,0,0,0,0,NULL,&atx[28],NULL,0,&atx[29]} ,
  {407, "Date" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[26]} ,
  {0, "tax-source" ,128,5,0,0,0,0,0,0,NULL,&atx[30],NULL,0,&atx[31]} ,
  {413, "Org-ref" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[42]} ,
  {0, "source" ,128,6,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[32]} ,
  {0, "status" ,128,7,0,0,0,0,0,0,NULL,&atx[2],&avnx[0],0,NULL} ,
  {314, "SET OF" ,0,17,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "seqannot" ,128,3,0,1,0,0,0,0,NULL,&atx[14],&atx[35],0,&atx[37]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[36],NULL,0,NULL} ,
  {411, "Seq-annot" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[40]} ,
  {0, "features" ,128,4,0,1,0,0,0,0,NULL,&atx[38],NULL,0,&atx[39]} ,
  {409, "Biostruc-annot-set" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[44]} ,
  {0, "sequences" ,128,5,0,1,0,0,0,0,NULL,&atx[40],NULL,0,&atx[41]} ,
  {412, "Seq-entry" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[30]} ,
  {0, "profile-range" ,128,6,0,1,0,0,0,0,NULL,&atx[42],NULL,0,&atx[43]} ,
  {414, "Seq-interval" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[59]} ,
  {0, "trunc-master" ,128,7,0,1,0,0,0,0,NULL,&atx[44],NULL,0,&atx[45]} ,
  {410, "Bioseq" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[36]} ,
  {0, "posfreq" ,128,8,0,1,0,0,0,0,NULL,&atx[46],NULL,0,&atx[54]} ,
  {419, "Matrix" ,1,0,0,0,0,0,0,0,NULL,&atx[10],&atx[47],0,&atx[56]} ,
  {0, "ncolumns" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[48]} ,
  {0, "nrows" ,128,1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[49]} ,
  {0, "row-labels" ,128,2,0,1,0,0,0,0,NULL,&atx[14],&atx[50],0,&atx[51]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[6],NULL,0,NULL} ,
  {0, "scale-factor" ,128,3,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[52]} ,
  {0, "columns" ,128,4,0,0,0,0,0,0,NULL,&atx[14],&atx[53],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "scoremat" ,128,9,0,1,0,0,0,0,NULL,&atx[46],NULL,0,&atx[55]} ,
  {0, "distance" ,128,10,0,1,0,0,0,0,NULL,&atx[56],NULL,0,NULL} ,
  {420, "Triangle" ,1,0,0,0,0,0,0,0,NULL,&atx[10],&atx[57],0,NULL} ,
  {0, "nelements" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[58]} ,
  {0, "scores" ,128,1,0,0,0,0,0,0,NULL,&atx[59],NULL,0,NULL} ,
  {415, "Score-set" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[4]} ,
  {404, "Cdd-set" ,1,0,0,0,0,1,0,0,NULL,&atx[33],&atx[61],0,&atx[62]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[15],NULL,0,NULL} ,
  {405, "Cdd-tree" ,1,0,0,0,0,1,0,0,NULL,&atx[10],&atx[63],0,&atx[69]} ,
  {0, "name" ,128,0,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[64]} ,
  {0, "id" ,128,1,0,0,0,0,0,0,NULL,&atx[12],NULL,0,&atx[65]} ,
  {0, "description" ,128,2,0,1,0,0,0,0,NULL,&atx[19],NULL,0,&atx[66]} ,
  {0, "parents" ,128,3,0,1,0,0,0,0,NULL,&atx[12],NULL,0,&atx[67]} ,
  {0, "children" ,128,4,0,1,0,0,0,0,NULL,&atx[12],NULL,0,&atx[68]} ,
  {0, "siblings" ,128,5,0,1,0,0,0,0,NULL,&atx[12],NULL,0,NULL} ,
  {406, "Cdd-tree-set" ,1,0,0,0,0,1,0,0,NULL,&atx[14],&atx[70],0,&atx[28]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[62],NULL,0,NULL} };

static AsnModule ampx[1] = {
  { "NCBI-Cdd" , "cdd.h",&atx[0],NULL,NULL,0,0} };

static AsnValxNodePtr avn = avnx;
static AsnTypePtr at = atx;
static AsnModulePtr amp = ampx;



/**************************************************
*
*    Defines for Module NCBI-Cdd
*
**************************************************/

#define CDD_ID &at[0]
#define CDD_ID_uid &at[1]
#define CDD_ID_gid &at[3]

#define CDD_ID_SET &at[12]
#define CDD_ID_SET_E &at[13]

#define CDD &at[15]
#define CDD_name &at[16]
#define CDD_id &at[17]
#define CDD_description &at[18]
#define CDD_seqannot &at[34]
#define CDD_seqannot_E &at[35]
#define CDD_features &at[37]
#define CDD_sequences &at[39]
#define CDD_profile_range &at[41]
#define CDD_trunc_master &at[43]
#define CDD_posfreq &at[45]
#define CDD_scoremat &at[54]
#define CDD_distance &at[55]

#define CDD_SET &at[60]
#define CDD_SET_E &at[61]

#define CDD_TREE &at[62]
#define CDD_TREE_name &at[63]
#define CDD_TREE_id &at[64]
#define CDD_TREE_description &at[65]
#define CDD_TREE_parents &at[66]
#define CDD_TREE_children &at[67]
#define CDD_TREE_siblings &at[68]

#define CDD_TREE_SET &at[69]
#define CDD_TREE_SET_E &at[70]

#define GLOBAL_ID &at[4]
#define GLOBAL_ID_accession &at[5]
#define GLOBAL_ID_release &at[7]
#define GLOBAL_ID_version &at[8]
#define GLOBAL_ID_database &at[9]

#define CDD_DESCR &at[21]
#define CDD_DESCR_othername &at[22]
#define CDD_DESCR_category &at[23]
#define CDD_DESCR_comment &at[24]
#define CDD_DESCR_reference &at[25]
#define CDD_DESCR_create_date &at[27]
#define CDD_DESCR_tax_source &at[29]
#define CDD_DESCR_source &at[31]
#define CDD_DESCR_status &at[32]

#define CDD_DESCR_SET &at[19]
#define CDD_DESCR_SET_E &at[20]

#define MATRIX &at[46]
#define MATRIX_ncolumns &at[47]
#define MATRIX_nrows &at[48]
#define MATRIX_row_labels &at[49]
#define MATRIX_row_labels_E &at[50]
#define MATRIX_scale_factor &at[51]
#define MATRIX_columns &at[52]
#define MATRIX_columns_E &at[53]

#define TRIANGLE &at[56]
#define TRIANGLE_nelements &at[57]
#define TRIANGLE_scores &at[58]
