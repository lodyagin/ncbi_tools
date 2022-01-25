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
static AsnValxNode avnx[15] = {
    {20,"unassigned" ,0,0.0,&avnx[1] } ,
    {20,"finished-ok" ,1,0.0,&avnx[2] } ,
    {20,"pending-release" ,2,0.0,&avnx[3] } ,
    {20,"other-asis" ,3,0.0,&avnx[4] } ,
    {20,"matrix-only" ,4,0.0,&avnx[5] } ,
    {20,"update-running" ,5,0.0,&avnx[6] } ,
    {20,"auto-updated" ,6,0.0,&avnx[7] } ,
    {20,"claimed" ,7,0.0,&avnx[8] } ,
    {20,"other" ,255,0.0,NULL } ,
    {20,"unassigned" ,0,0.0,&avnx[10] } ,
    {20,"update" ,1,0.0,&avnx[11] } ,
    {20,"update-3d" ,2,0.0,&avnx[12] } ,
    {20,"demoted" ,51,0.0,&avnx[13] } ,
    {20,"demoted-3d" ,52,0.0,&avnx[14] } ,
    {20,"other" ,255,0.0,NULL } };

static AsnType atx[124] = {
  {401, "Cdd-id" ,1,0,0,0,0,1,0,0,NULL,&atx[11],&atx[1],0,&atx[12]} ,
  {0, "uid" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[3]} ,
  {302, "INTEGER" ,0,2,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "gid" ,128,1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {421, "Global-id" ,1,0,0,0,0,0,0,0,NULL,&atx[10],&atx[5],0,&atx[21]} ,
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
  {403, "Cdd" ,1,0,0,0,0,1,0,0,NULL,&atx[10],&atx[16],0,&atx[112]} ,
  {0, "name" ,128,0,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[17]} ,
  {0, "id" ,128,1,0,0,0,0,0,0,NULL,&atx[12],NULL,0,&atx[18]} ,
  {0, "description" ,128,2,0,1,0,0,0,0,NULL,&atx[19],NULL,0,&atx[38]} ,
  {423, "Cdd-descr-set" ,1,0,0,0,0,0,0,0,NULL,&atx[37],&atx[20],0,&atx[50]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[21],NULL,0,NULL} ,
  {422, "Cdd-descr" ,1,0,0,0,0,0,0,0,NULL,&atx[11],&atx[22],0,&atx[19]} ,
  {0, "othername" ,128,0,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[23]} ,
  {0, "category" ,128,1,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[24]} ,
  {0, "comment" ,128,2,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[25]} ,
  {0, "reference" ,128,3,0,0,0,0,0,0,NULL,&atx[26],NULL,0,&atx[27]} ,
  {408, "Pub" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[42]} ,
  {0, "create-date" ,128,4,0,0,0,0,0,0,NULL,&atx[28],NULL,0,&atx[29]} ,
  {407, "Date" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[26]} ,
  {0, "tax-source" ,128,5,0,0,0,0,0,0,NULL,&atx[30],NULL,0,&atx[31]} ,
  {413, "Org-ref" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[91]} ,
  {0, "source" ,128,6,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[32]} ,
  {0, "status" ,128,7,0,0,0,0,0,0,NULL,&atx[2],&avnx[0],0,&atx[33]} ,
  {0, "update-date" ,128,8,0,0,0,0,0,0,NULL,&atx[28],NULL,0,&atx[34]} ,
  {0, "scrapbook" ,128,9,0,0,0,0,0,0,NULL,&atx[14],&atx[35],0,&atx[36]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[6],NULL,0,NULL} ,
  {0, "source-id" ,128,10,0,0,0,0,0,0,NULL,&atx[12],NULL,0,NULL} ,
  {314, "SET OF" ,0,17,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "seqannot" ,128,3,0,1,0,0,0,0,NULL,&atx[14],&atx[39],0,&atx[41]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[40],NULL,0,NULL} ,
  {411, "Seq-annot" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[44]} ,
  {0, "features" ,128,4,0,1,0,0,0,0,NULL,&atx[42],NULL,0,&atx[43]} ,
  {409, "Biostruc-annot-set" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[48]} ,
  {0, "sequences" ,128,5,0,1,0,0,0,0,NULL,&atx[44],NULL,0,&atx[45]} ,
  {412, "Seq-entry" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[30]} ,
  {0, "profile-range" ,128,6,0,1,0,0,0,0,NULL,&atx[46],NULL,0,&atx[47]} ,
  {415, "Seq-interval" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[78]} ,
  {0, "trunc-master" ,128,7,0,1,0,0,0,0,NULL,&atx[48],NULL,0,&atx[49]} ,
  {410, "Bioseq" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[40]} ,
  {0, "posfreq" ,128,8,0,1,0,0,0,0,NULL,&atx[50],NULL,0,&atx[58]} ,
  {424, "Matrix" ,1,0,0,0,0,0,0,0,NULL,&atx[10],&atx[51],0,&atx[60]} ,
  {0, "ncolumns" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[52]} ,
  {0, "nrows" ,128,1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[53]} ,
  {0, "row-labels" ,128,2,0,1,0,0,0,0,NULL,&atx[14],&atx[54],0,&atx[55]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[6],NULL,0,NULL} ,
  {0, "scale-factor" ,128,3,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[56]} ,
  {0, "columns" ,128,4,0,0,0,0,0,0,NULL,&atx[14],&atx[57],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "scoremat" ,128,9,0,1,0,0,0,0,NULL,&atx[50],NULL,0,&atx[59]} ,
  {0, "distance" ,128,10,0,1,0,0,0,0,NULL,&atx[60],NULL,0,&atx[66]} ,
  {425, "Triangle" ,1,0,0,0,0,0,0,0,NULL,&atx[10],&atx[61],0,&atx[75]} ,
  {0, "nelements" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[62]} ,
  {0, "scores" ,128,1,0,1,0,0,0,0,NULL,&atx[63],NULL,0,&atx[64]} ,
  {418, "Score-set" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[109]} ,
  {0, "div-ranks" ,128,2,0,1,0,0,0,0,NULL,&atx[14],&atx[65],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "parents" ,128,11,0,1,0,0,0,0,NULL,&atx[12],NULL,0,&atx[67]} ,
  {0, "children" ,128,12,0,1,0,0,0,0,NULL,&atx[12],NULL,0,&atx[68]} ,
  {0, "siblings" ,128,13,0,1,0,0,0,0,NULL,&atx[12],NULL,0,&atx[69]} ,
  {0, "neighbors" ,128,14,0,1,0,0,0,0,NULL,&atx[12],NULL,0,&atx[70]} ,
  {0, "pending" ,128,15,0,1,0,0,0,0,NULL,&atx[14],&atx[71],0,&atx[84]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[72],NULL,0,NULL} ,
  {427, "Update-align" ,1,0,0,0,0,0,0,0,NULL,&atx[10],&atx[73],0,&atx[86]} ,
  {0, "description" ,128,0,0,1,0,0,0,0,NULL,&atx[14],&atx[74],0,&atx[82]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[75],NULL,0,NULL} ,
  {426, "Update-comment" ,1,0,0,0,0,0,0,0,NULL,&atx[11],&atx[76],0,&atx[72]} ,
  {0, "comment" ,128,0,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[77]} ,
  {0, "addthis" ,128,1,0,0,0,0,0,0,NULL,&atx[78],NULL,0,&atx[79]} ,
  {416, "Seq-loc" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[107]} ,
  {0, "replaces" ,128,2,0,0,0,0,0,0,NULL,&atx[78],NULL,0,&atx[80]} ,
  {0, "reject-loc" ,128,3,0,0,0,0,0,0,NULL,&atx[78],NULL,0,&atx[81]} ,
  {0, "reference" ,128,4,0,0,0,0,0,0,NULL,&atx[26],NULL,0,NULL} ,
  {0, "seqannot" ,128,1,0,1,0,0,0,0,NULL,&atx[40],NULL,0,&atx[83]} ,
  {0, "type" ,128,2,0,0,0,0,0,0,NULL,&atx[2],&avnx[9],0,NULL} ,
  {0, "rejects" ,128,16,0,1,0,0,0,0,NULL,&atx[14],&atx[85],0,&atx[92]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[86],NULL,0,NULL} ,
  {428, "Reject-id" ,1,0,0,0,0,0,0,0,NULL,&atx[10],&atx[87],0,&atx[102]} ,
  {0, "description" ,128,0,0,1,0,0,0,0,NULL,&atx[14],&atx[88],0,&atx[89]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[75],NULL,0,NULL} ,
  {0, "ids" ,128,1,0,0,0,0,0,0,NULL,&atx[37],&atx[90],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[91],NULL,0,NULL} ,
  {414, "Seq-id" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[46]} ,
  {0, "master3d" ,128,17,0,1,0,0,0,0,NULL,&atx[37],&atx[93],0,&atx[94]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[91],NULL,0,NULL} ,
  {0, "alignannot" ,128,18,0,1,0,0,0,0,NULL,&atx[95],NULL,0,&atx[108]} ,
  {431, "Align-annot-set" ,1,0,0,0,0,0,0,0,NULL,&atx[14],&atx[96],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[97],NULL,0,NULL} ,
  {430, "Align-annot" ,1,0,0,0,0,0,0,0,NULL,&atx[10],&atx[98],0,&atx[95]} ,
  {0, "location" ,128,0,0,0,0,0,0,0,NULL,&atx[78],NULL,0,&atx[99]} ,
  {0, "description" ,128,1,0,1,0,0,0,0,NULL,&atx[6],NULL,0,&atx[100]} ,
  {0, "evidence" ,128,2,0,1,0,0,0,0,NULL,&atx[14],&atx[101],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[102],NULL,0,NULL} ,
  {429, "Feature-evidence" ,1,0,0,0,0,0,0,0,NULL,&atx[11],&atx[103],0,&atx[97]} ,
  {0, "comment" ,128,0,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[104]} ,
  {0, "reference" ,128,1,0,0,0,0,0,0,NULL,&atx[26],NULL,0,&atx[105]} ,
  {0, "bsannot" ,128,2,0,0,0,0,0,0,NULL,&atx[42],NULL,0,&atx[106]} ,
  {0, "seqfeat" ,128,3,0,0,0,0,0,0,NULL,&atx[107],NULL,0,NULL} ,
  {417, "Seq-feat" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[63]} ,
  {0, "style-dictionary" ,128,19,0,1,0,0,0,0,NULL,&atx[109],NULL,0,&atx[110]} ,
  {419, "Cn3d-style-dictionary" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[111]} ,
  {0, "user-annotations" ,128,20,0,1,0,0,0,0,NULL,&atx[111],NULL,0,NULL} ,
  {420, "Cn3d-user-annotations" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[4]} ,
  {404, "Cdd-set" ,1,0,0,0,0,1,0,0,NULL,&atx[37],&atx[113],0,&atx[114]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[15],NULL,0,NULL} ,
  {405, "Cdd-tree" ,1,0,0,0,0,1,0,0,NULL,&atx[10],&atx[115],0,&atx[122]} ,
  {0, "name" ,128,0,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[116]} ,
  {0, "id" ,128,1,0,0,0,0,0,0,NULL,&atx[12],NULL,0,&atx[117]} ,
  {0, "description" ,128,2,0,1,0,0,0,0,NULL,&atx[19],NULL,0,&atx[118]} ,
  {0, "parents" ,128,3,0,1,0,0,0,0,NULL,&atx[12],NULL,0,&atx[119]} ,
  {0, "children" ,128,4,0,1,0,0,0,0,NULL,&atx[12],NULL,0,&atx[120]} ,
  {0, "siblings" ,128,5,0,1,0,0,0,0,NULL,&atx[12],NULL,0,&atx[121]} ,
  {0, "neighbors" ,128,6,0,1,0,0,0,0,NULL,&atx[12],NULL,0,NULL} ,
  {406, "Cdd-tree-set" ,1,0,0,0,0,1,0,0,NULL,&atx[14],&atx[123],0,&atx[28]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[114],NULL,0,NULL} };

static AsnModule ampx[1] = {
  { "NCBI-Cdd" , "cdd.h18",&atx[0],NULL,NULL,0,0} };

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
#define CDD_seqannot &at[38]
#define CDD_seqannot_E &at[39]
#define CDD_features &at[41]
#define CDD_sequences &at[43]
#define CDD_profile_range &at[45]
#define CDD_trunc_master &at[47]
#define CDD_posfreq &at[49]
#define CDD_scoremat &at[58]
#define CDD_distance &at[59]
#define CDD_parents &at[66]
#define CDD_children &at[67]
#define CDD_siblings &at[68]
#define CDD_neighbors &at[69]
#define CDD_pending &at[70]
#define CDD_pending_E &at[71]
#define CDD_rejects &at[84]
#define CDD_rejects_E &at[85]
#define CDD_master3d &at[92]
#define CDD_master3d_E &at[93]
#define CDD_alignannot &at[94]
#define CDD_style_dictionary &at[108]
#define CDD_user_annotations &at[110]

#define CDD_SET &at[112]
#define CDD_SET_E &at[113]

#define CDD_TREE &at[114]
#define CDD_TREE_name &at[115]
#define CDD_TREE_id &at[116]
#define CDD_TREE_description &at[117]
#define CDD_TREE_parents &at[118]
#define CDD_TREE_children &at[119]
#define CDD_TREE_siblings &at[120]
#define CDD_TREE_neighbors &at[121]

#define CDD_TREE_SET &at[122]
#define CDD_TREE_SET_E &at[123]

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
#define CDD_DESCR_update_date &at[33]
#define CDD_DESCR_scrapbook &at[34]
#define CDD_DESCR_scrapbook_E &at[35]
#define CDD_DESCR_source_id &at[36]

#define CDD_DESCR_SET &at[19]
#define CDD_DESCR_SET_E &at[20]

#define MATRIX &at[50]
#define MATRIX_ncolumns &at[51]
#define MATRIX_nrows &at[52]
#define MATRIX_row_labels &at[53]
#define MATRIX_row_labels_E &at[54]
#define MATRIX_scale_factor &at[55]
#define MATRIX_columns &at[56]
#define MATRIX_columns_E &at[57]

#define TRIANGLE &at[60]
#define TRIANGLE_nelements &at[61]
#define TRIANGLE_scores &at[62]
#define TRIANGLE_div_ranks &at[64]
#define TRIANGLE_div_ranks_E &at[65]

#define UPDATE_COMMENT &at[75]
#define UPDATE_COMMENT_comment &at[76]
#define UPDATE_COMMENT_addthis &at[77]
#define UPDATE_COMMENT_replaces &at[79]
#define UPDATE_COMMENT_reject_loc &at[80]
#define UPDATE_COMMENT_reference &at[81]

#define UPDATE_ALIGN &at[72]
#define UPDATE_ALIGN_description &at[73]
#define UPDATE_ALIGN_description_E &at[74]
#define UPDATE_ALIGN_seqannot &at[82]
#define UPDATE_ALIGN_type &at[83]

#define REJECT_ID &at[86]
#define REJECT_ID_description &at[87]
#define REJECT_ID_description_E &at[88]
#define REJECT_ID_ids &at[89]
#define REJECT_ID_ids_E &at[90]

#define FEATURE_EVIDENCE &at[102]
#define FEATURE_EVIDENCE_comment &at[103]
#define FEATURE_EVIDENCE_reference &at[104]
#define FEATURE_EVIDENCE_bsannot &at[105]
#define FEATURE_EVIDENCE_seqfeat &at[106]

#define ALIGN_ANNOT &at[97]
#define ALIGN_ANNOT_location &at[98]
#define ALIGN_ANNOT_description &at[99]
#define ALIGN_ANNOT_evidence &at[100]
#define ALIGN_ANNOT_evidence_E &at[101]

#define ALIGN_ANNOT_SET &at[95]
#define ALIGN_ANNOT_SET_E &at[96]
