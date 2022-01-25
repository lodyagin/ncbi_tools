/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "asnmime.h64";
static AsnValxNode avnx[13] = {
    {20,"docsum" ,1,0.0,&avnx[1] } ,
    {20,"genbank" ,2,0.0,&avnx[2] } ,
    {20,"genpept" ,3,0.0,&avnx[3] } ,
    {20,"fasta" ,4,0.0,&avnx[4] } ,
    {20,"asn1" ,5,0.0,&avnx[5] } ,
    {20,"graphic" ,6,0.0,&avnx[6] } ,
    {20,"alignment" ,7,0.0,&avnx[7] } ,
    {20,"globalview" ,8,0.0,&avnx[8] } ,
    {20,"report" ,9,0.0,&avnx[9] } ,
    {20,"medlars" ,10,0.0,&avnx[10] } ,
    {20,"embl" ,11,0.0,&avnx[11] } ,
    {20,"pdb" ,12,0.0,&avnx[12] } ,
    {20,"kinemage" ,13,0.0,NULL } };

static AsnType atx[52] = {
  {401, "Ncbi-mime-asn1" ,1,0,0,0,0,1,0,0,NULL,&atx[16],&atx[1],0,&atx[13]} ,
  {0, "entrez" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[22]} ,
  {407, "Entrez-general" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[3],0,&atx[23]} ,
  {0, "title" ,128,0,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[5]} ,
  {323, "VisibleString" ,0,26,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "data" ,128,1,0,0,0,0,0,0,NULL,&atx[16],&atx[6],0,&atx[17]} ,
  {0, "ml" ,128,0,0,0,0,0,0,0,NULL,&atx[7],NULL,0,&atx[8]} ,
  {406, "Medline-entry" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[2]} ,
  {0, "prot" ,128,1,0,0,0,0,0,0,NULL,&atx[9],NULL,0,&atx[10]} ,
  {404, "Seq-entry" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[33]} ,
  {0, "nuc" ,128,2,0,0,0,0,0,0,NULL,&atx[9],NULL,0,&atx[11]} ,
  {0, "genome" ,128,3,0,0,0,0,0,0,NULL,&atx[9],NULL,0,&atx[12]} ,
  {0, "structure" ,128,4,0,0,0,0,0,0,NULL,&atx[13],NULL,0,&atx[14]} ,
  {402, "Biostruc" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[15]} ,
  {0, "strucAnnot" ,128,5,0,0,0,0,0,0,NULL,&atx[15],NULL,0,NULL} ,
  {403, "Biostruc-annot-set" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[9]} ,
  {315, "CHOICE" ,0,-1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "style" ,128,2,0,0,0,0,0,0,NULL,&atx[18],NULL,0,&atx[20]} ,
  {412, "Entrez-style" ,1,0,0,0,0,0,0,0,NULL,&atx[19],&avnx[0],0,NULL} ,
  {310, "ENUMERATED" ,0,10,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "location" ,128,3,0,1,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {311, "SEQUENCE" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "alignstruc" ,128,1,0,0,0,0,0,0,NULL,&atx[23],NULL,0,&atx[34]} ,
  {408, "Biostruc-align" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[24],0,&atx[35]} ,
  {0, "master" ,128,0,0,0,0,0,0,0,NULL,&atx[13],NULL,0,&atx[25]} ,
  {0, "slaves" ,128,1,0,0,0,0,0,0,NULL,&atx[27],&atx[26],0,&atx[28]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[13],NULL,0,NULL} ,
  {314, "SET OF" ,0,17,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "alignments" ,128,2,0,0,0,0,0,0,NULL,&atx[15],NULL,0,&atx[29]} ,
  {0, "sequences" ,128,3,0,0,0,0,0,0,NULL,&atx[27],&atx[30],0,&atx[31]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[9],NULL,0,NULL} ,
  {0, "seqalign" ,128,4,0,0,0,0,0,0,NULL,&atx[27],&atx[32],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[33],NULL,0,NULL} ,
  {405, "Seq-annot" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[7]} ,
  {0, "alignseq" ,128,2,0,0,0,0,0,0,NULL,&atx[35],NULL,0,&atx[40]} ,
  {409, "Biostruc-align-seq" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[36],0,&atx[41]} ,
  {0, "sequences" ,128,0,0,0,0,0,0,0,NULL,&atx[27],&atx[37],0,&atx[38]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[9],NULL,0,NULL} ,
  {0, "seqalign" ,128,1,0,0,0,0,0,0,NULL,&atx[27],&atx[39],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[33],NULL,0,NULL} ,
  {0, "strucseq" ,128,3,0,0,0,0,0,0,NULL,&atx[41],NULL,0,&atx[45]} ,
  {410, "Biostruc-seq" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[42],0,&atx[46]} ,
  {0, "structure" ,128,0,0,0,0,0,0,0,NULL,&atx[13],NULL,0,&atx[43]} ,
  {0, "sequences" ,128,1,0,0,0,0,0,0,NULL,&atx[27],&atx[44],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[9],NULL,0,NULL} ,
  {0, "strucseqs" ,128,4,0,0,0,0,0,0,NULL,&atx[46],NULL,0,NULL} ,
  {411, "Biostruc-seqs" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[47],0,&atx[18]} ,
  {0, "structure" ,128,0,0,0,0,0,0,0,NULL,&atx[13],NULL,0,&atx[48]} ,
  {0, "sequences" ,128,1,0,0,0,0,0,0,NULL,&atx[27],&atx[49],0,&atx[50]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[9],NULL,0,NULL} ,
  {0, "seqalign" ,128,2,0,0,0,0,0,0,NULL,&atx[27],&atx[51],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[33],NULL,0,NULL} };

static AsnModule ampx[1] = {
  { "NCBI-Mime" , "asnmime.h64",&atx[0],NULL,NULL,0,0} };

static AsnValxNodePtr avn = avnx;
static AsnTypePtr at = atx;
static AsnModulePtr amp = ampx;



/**************************************************
*
*    Defines for Module NCBI-Mime
*
**************************************************/

#define NCBI_MIME_ASN1 &at[0]
#define NCBI_MIME_ASN1_entrez &at[1]
#define NCBI_MIME_ASN1_alignstruc &at[22]
#define NCBI_MIME_ASN1_alignseq &at[34]
#define NCBI_MIME_ASN1_strucseq &at[40]
#define NCBI_MIME_ASN1_strucseqs &at[45]

#define ENTREZ_GENERAL &at[2]
#define ENTREZ_GENERAL_title &at[3]
#define ENTREZ_GENERAL_data &at[5]
#define ENTREZ_GENERAL_data_ml &at[6]
#define ENTREZ_GENERAL_data_prot &at[8]
#define ENTREZ_GENERAL_data_nuc &at[10]
#define ENTREZ_GENERAL_data_genome &at[11]
#define ENTREZ_GENERAL_data_structure &at[12]
#define ENTREZ_GENERAL_data_strucAnnot &at[14]
#define ENTREZ_GENERAL_style &at[17]
#define ENTREZ_GENERAL_location &at[20]

#define BIOSTRUC_ALIGN &at[23]
#define BIOSTRUC_ALIGN_master &at[24]
#define BIOSTRUC_ALIGN_slaves &at[25]
#define BIOSTRUC_ALIGN_slaves_E &at[26]
#define BIOSTRUC_ALIGN_alignments &at[28]
#define BIOSTRUC_ALIGN_sequences &at[29]
#define BIOSTRUC_ALIGN_sequences_E &at[30]
#define BIOSTRUC_ALIGN_seqalign &at[31]
#define BIOSTRUC_ALIGN_seqalign_E &at[32]

#define BIOSTRUC_ALIGN_SEQ &at[35]
#define BIOSTRUC_ALIGN_SEQ_sequences &at[36]
#define BIOSTRUC_ALIGN_SEQ_sequences_E &at[37]
#define BIOSTRUC_ALIGN_SEQ_seqalign &at[38]
#define BIOSTRUC_ALIGN_SEQ_seqalign_E &at[39]

#define BIOSTRUC_SEQ &at[41]
#define BIOSTRUC_SEQ_structure &at[42]
#define BIOSTRUC_SEQ_sequences &at[43]
#define BIOSTRUC_SEQ_sequences_E &at[44]

#define BIOSTRUC_SEQS &at[46]
#define BIOSTRUC_SEQS_structure &at[47]
#define BIOSTRUC_SEQS_sequences &at[48]
#define BIOSTRUC_SEQS_sequences_E &at[49]
#define BIOSTRUC_SEQS_seqalign &at[50]
#define BIOSTRUC_SEQS_seqalign_E &at[51]

#define ENTREZ_STYLE &at[18]
