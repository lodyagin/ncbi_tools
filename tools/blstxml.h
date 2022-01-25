/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "blstxml.h";
static AsnType atx[66] = {
  {401, "BlastOutput" ,1,0,0,0,0,0,0,0,NULL,&atx[42],&atx[1],0,&atx[14]} ,
  {0, "program" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[3]} ,
  {323, "VisibleString" ,0,26,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "version" ,128,1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[4]} ,
  {0, "reference" ,128,2,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[5]} ,
  {0, "db" ,128,3,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[6]} ,
  {0, "query-ID" ,128,4,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[7]} ,
  {0, "query-def" ,128,5,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[8]} ,
  {0, "query-len" ,128,6,0,0,0,0,0,0,NULL,&atx[9],NULL,0,&atx[10]} ,
  {302, "INTEGER" ,0,2,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "query-seq" ,128,7,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[11]} ,
  {0, "iter-num" ,128,8,0,1,0,0,0,0,NULL,&atx[9],NULL,0,&atx[12]} ,
  {0, "hits" ,128,9,0,1,0,0,0,0,NULL,&atx[43],&atx[13],0,&atx[44]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[14],NULL,0,NULL} ,
  {402, "Hit" ,1,0,0,0,0,0,0,0,NULL,&atx[42],&atx[15],0,&atx[45]} ,
  {0, "num" ,128,0,0,0,0,0,0,0,NULL,&atx[9],NULL,0,&atx[16]} ,
  {0, "id" ,128,1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[17]} ,
  {0, "def" ,128,2,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[18]} ,
  {0, "accession" ,128,3,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[19]} ,
  {0, "len" ,128,4,0,0,0,0,0,0,NULL,&atx[9],NULL,0,&atx[20]} ,
  {0, "hsps" ,128,5,0,1,0,0,0,0,NULL,&atx[43],&atx[21],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[22],NULL,0,NULL} ,
  {405, "Hsp" ,1,0,0,0,0,0,0,0,NULL,&atx[42],&atx[23],0,NULL} ,
  {0, "num" ,128,0,0,0,0,0,0,0,NULL,&atx[9],NULL,0,&atx[24]} ,
  {0, "score" ,128,1,0,0,0,0,0,0,NULL,&atx[25],NULL,0,&atx[26]} ,
  {309, "REAL" ,0,9,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "evalue" ,128,2,0,0,0,0,0,0,NULL,&atx[25],NULL,0,&atx[27]} ,
  {0, "query-from" ,128,3,0,0,0,0,0,0,NULL,&atx[9],NULL,0,&atx[28]} ,
  {0, "query-to" ,128,4,0,0,0,0,0,0,NULL,&atx[9],NULL,0,&atx[29]} ,
  {0, "hit-from" ,128,5,0,0,0,0,0,0,NULL,&atx[9],NULL,0,&atx[30]} ,
  {0, "hit-to" ,128,6,0,0,0,0,0,0,NULL,&atx[9],NULL,0,&atx[31]} ,
  {0, "pattern-from" ,128,7,0,1,0,0,0,0,NULL,&atx[9],NULL,0,&atx[32]} ,
  {0, "pattern-to" ,128,8,0,1,0,0,0,0,NULL,&atx[9],NULL,0,&atx[33]} ,
  {0, "query-frame" ,128,9,0,1,0,0,0,0,NULL,&atx[9],NULL,0,&atx[34]} ,
  {0, "hit-frame" ,128,10,0,1,0,0,0,0,NULL,&atx[9],NULL,0,&atx[35]} ,
  {0, "identity" ,128,11,0,1,0,0,0,0,NULL,&atx[9],NULL,0,&atx[36]} ,
  {0, "positive" ,128,12,0,1,0,0,0,0,NULL,&atx[9],NULL,0,&atx[37]} ,
  {0, "gaps" ,128,13,0,1,0,0,0,0,NULL,&atx[9],NULL,0,&atx[38]} ,
  {0, "density" ,128,14,0,1,0,0,0,0,NULL,&atx[9],NULL,0,&atx[39]} ,
  {0, "qseq" ,128,15,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[40]} ,
  {0, "hseq" ,128,16,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[41]} ,
  {0, "midline" ,128,17,0,1,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {311, "SEQUENCE" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {312, "SEQUENCE OF" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "param" ,128,10,0,0,0,0,0,0,NULL,&atx[45],NULL,0,&atx[56]} ,
  {403, "Parameters" ,1,0,0,0,0,0,0,0,NULL,&atx[42],&atx[46],0,&atx[57]} ,
  {0, "matrix" ,128,0,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[47]} ,
  {0, "expect" ,128,1,0,0,0,0,0,0,NULL,&atx[25],NULL,0,&atx[48]} ,
  {0, "include" ,128,2,0,1,0,0,0,0,NULL,&atx[25],NULL,0,&atx[49]} ,
  {0, "sc-match" ,128,3,0,1,0,0,0,0,NULL,&atx[9],NULL,0,&atx[50]} ,
  {0, "sc-mismatch" ,128,4,0,1,0,0,0,0,NULL,&atx[9],NULL,0,&atx[51]} ,
  {0, "gap-open" ,128,5,0,0,0,0,0,0,NULL,&atx[9],NULL,0,&atx[52]} ,
  {0, "gap-extend" ,128,6,0,0,0,0,0,0,NULL,&atx[9],NULL,0,&atx[53]} ,
  {0, "filter" ,128,7,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[54]} ,
  {0, "pattern" ,128,8,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[55]} ,
  {0, "entrez-query" ,128,9,0,1,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "stat" ,128,11,0,1,0,0,0,0,NULL,&atx[57],NULL,0,&atx[65]} ,
  {404, "Statistics" ,1,0,0,0,0,0,0,0,NULL,&atx[42],&atx[58],0,&atx[22]} ,
  {0, "db-num" ,128,0,0,0,0,0,0,0,NULL,&atx[9],NULL,0,&atx[59]} ,
  {0, "db-len" ,128,1,0,0,0,0,0,0,NULL,&atx[9],NULL,0,&atx[60]} ,
  {0, "hsp-len" ,128,2,0,0,0,0,0,0,NULL,&atx[9],NULL,0,&atx[61]} ,
  {0, "eff-space" ,128,3,0,0,0,0,0,0,NULL,&atx[25],NULL,0,&atx[62]} ,
  {0, "kappa" ,128,4,0,0,0,0,0,0,NULL,&atx[25],NULL,0,&atx[63]} ,
  {0, "lambda" ,128,5,0,0,0,0,0,0,NULL,&atx[25],NULL,0,&atx[64]} ,
  {0, "entropy" ,128,6,0,0,0,0,0,0,NULL,&atx[25],NULL,0,NULL} ,
  {0, "message" ,128,12,0,1,0,0,0,0,NULL,&atx[2],NULL,0,NULL} };

static AsnModule ampx[1] = {
  { "NCBI-BlastOutput" , "blstxml.h",&atx[0],NULL,NULL,0,0} };

static AsnValxNodePtr avn = NULL;
static AsnTypePtr at = atx;
static AsnModulePtr amp = ampx;



/**************************************************
*
*    Defines for Module NCBI-BlastOutput
*
**************************************************/

#define BLASTOUTPUT &at[0]
#define BLASTOUTPUT_program &at[1]
#define BLASTOUTPUT_version &at[3]
#define BLASTOUTPUT_reference &at[4]
#define BLASTOUTPUT_db &at[5]
#define BLASTOUTPUT_query_ID &at[6]
#define BLASTOUTPUT_query_def &at[7]
#define BLASTOUTPUT_query_len &at[8]
#define BLASTOUTPUT_query_seq &at[10]
#define BLASTOUTPUT_iter_num &at[11]
#define BLASTOUTPUT_hits &at[12]
#define BLASTOUTPUT_hits_E &at[13]
#define BLASTOUTPUT_param &at[44]
#define BLASTOUTPUT_stat &at[56]
#define BLASTOUTPUT_message &at[65]

#define HIT &at[14]
#define HIT_num &at[15]
#define HIT_id &at[16]
#define HIT_def &at[17]
#define HIT_accession &at[18]
#define HIT_len &at[19]
#define HIT_hsps &at[20]
#define HIT_hsps_E &at[21]

#define PARAMETERS &at[45]
#define PARAMETERS_matrix &at[46]
#define PARAMETERS_expect &at[47]
#define PARAMETERS_include &at[48]
#define PARAMETERS_sc_match &at[49]
#define PARAMETERS_sc_mismatch &at[50]
#define PARAMETERS_gap_open &at[51]
#define PARAMETERS_gap_extend &at[52]
#define PARAMETERS_filter &at[53]
#define PARAMETERS_pattern &at[54]
#define PARAMETERS_entrez_query &at[55]

#define STATISTICS &at[57]
#define STATISTICS_db_num &at[58]
#define STATISTICS_db_len &at[59]
#define STATISTICS_hsp_len &at[60]
#define STATISTICS_eff_space &at[61]
#define STATISTICS_kappa &at[62]
#define STATISTICS_lambda &at[63]
#define STATISTICS_entropy &at[64]

#define HSP &at[22]
#define HSP_num &at[23]
#define HSP_score &at[24]
#define HSP_evalue &at[26]
#define HSP_query_from &at[27]
#define HSP_query_to &at[28]
#define HSP_hit_from &at[29]
#define HSP_hit_to &at[30]
#define HSP_pattern_from &at[31]
#define HSP_pattern_to &at[32]
#define HSP_query_frame &at[33]
#define HSP_hit_frame &at[34]
#define HSP_identity &at[35]
#define HSP_positive &at[36]
#define HSP_gaps &at[37]
#define HSP_density &at[38]
#define HSP_qseq &at[39]
#define HSP_hseq &at[40]
#define HSP_midline &at[41]
