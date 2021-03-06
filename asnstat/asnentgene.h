/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "asnentgene.h36";
static AsnValxNode avnx[52] = {
    {20,"live" ,0,0.0,&avnx[1] } ,
    {20,"secondary" ,1,0.0,&avnx[2] } ,
    {20,"discontinued" ,2,0.0,NULL } ,
    {3,NULL,0,0.0,NULL } ,
    {20,"unknown" ,0,0.0,&avnx[5] } ,
    {20,"tRNA" ,1,0.0,&avnx[6] } ,
    {20,"rRNA" ,2,0.0,&avnx[7] } ,
    {20,"snRNA" ,3,0.0,&avnx[8] } ,
    {20,"scRNA" ,4,0.0,&avnx[9] } ,
    {20,"snoRNA" ,5,0.0,&avnx[10] } ,
    {20,"protein-coding" ,6,0.0,&avnx[11] } ,
    {20,"pseudo" ,7,0.0,&avnx[12] } ,
    {20,"transposon" ,8,0.0,&avnx[13] } ,
    {20,"miscRNA" ,9,0.0,&avnx[14] } ,
    {20,"ncRNA" ,10,0.0,&avnx[15] } ,
    {20,"biological-region" ,11,0.0,&avnx[16] } ,
    {20,"other" ,255,0.0,NULL } ,
    {20,"cyto" ,0,0.0,&avnx[18] } ,
    {20,"bp" ,1,0.0,&avnx[19] } ,
    {20,"cM" ,2,0.0,&avnx[20] } ,
    {20,"cR" ,3,0.0,&avnx[21] } ,
    {20,"min" ,4,0.0,NULL } ,
    {2,NULL,0,0.0,NULL } ,
    {2,NULL,0,0.0,NULL } ,
    {2,NULL,0,0.0,NULL } ,
    {20,"genomic" ,1,0.0,&avnx[26] } ,
    {20,"pre-RNA" ,2,0.0,&avnx[27] } ,
    {20,"mRNA" ,3,0.0,&avnx[28] } ,
    {20,"rRNA" ,4,0.0,&avnx[29] } ,
    {20,"tRNA" ,5,0.0,&avnx[30] } ,
    {20,"snRNA" ,6,0.0,&avnx[31] } ,
    {20,"scRNA" ,7,0.0,&avnx[32] } ,
    {20,"peptide" ,8,0.0,&avnx[33] } ,
    {20,"other-genetic" ,9,0.0,&avnx[34] } ,
    {20,"genomic-mRNA" ,10,0.0,&avnx[35] } ,
    {20,"cRNA" ,11,0.0,&avnx[36] } ,
    {20,"mature-peptide" ,12,0.0,&avnx[37] } ,
    {20,"pre-protein" ,13,0.0,&avnx[38] } ,
    {20,"miscRNA" ,14,0.0,&avnx[39] } ,
    {20,"snoRNA" ,15,0.0,&avnx[40] } ,
    {20,"property" ,16,0.0,&avnx[41] } ,
    {20,"reference" ,17,0.0,&avnx[42] } ,
    {20,"generif" ,18,0.0,&avnx[43] } ,
    {20,"phenotype" ,19,0.0,&avnx[44] } ,
    {20,"complex" ,20,0.0,&avnx[45] } ,
    {20,"compound" ,21,0.0,&avnx[46] } ,
    {20,"ncRNA" ,22,0.0,&avnx[47] } ,
    {20,"gene-group" ,23,0.0,&avnx[48] } ,
    {20,"assembly" ,24,0.0,&avnx[49] } ,
    {20,"assembly-unit" ,25,0.0,&avnx[50] } ,
    {20,"comment" ,254,0.0,&avnx[51] } ,
    {20,"other" ,255,0.0,NULL } };

static AsnType atx[105] = {
  {401, "Entrezgene" ,1,0,0,0,0,1,0,0,NULL,&atx[14],&atx[1],0,&atx[102]} ,
  {0, "track-info" ,128,0,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[15]} ,
  {403, "Gene-track" ,1,0,0,0,0,1,0,0,NULL,&atx[14],&atx[3],0,&atx[47]} ,
  {0, "geneid" ,128,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[5]} ,
  {302, "INTEGER" ,0,2,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "status" ,128,1,0,0,1,0,0,0,&avnx[3],&atx[4],&avnx[0],0,&atx[6]} ,
  {0, "current-id" ,128,2,0,1,0,0,0,0,NULL,&atx[9],&atx[7],0,&atx[10]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[8],NULL,0,NULL} ,
  {409, "Dbtag" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[11]} ,
  {312, "SEQUENCE OF" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "create-date" ,128,3,0,0,0,0,0,0,NULL,&atx[11],NULL,0,&atx[12]} ,
  {410, "Date" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[72]} ,
  {0, "update-date" ,128,4,0,0,0,0,0,0,NULL,&atx[11],NULL,0,&atx[13]} ,
  {0, "discontinue-date" ,128,5,0,1,0,0,0,0,NULL,&atx[11],NULL,0,NULL} ,
  {311, "SEQUENCE" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "type" ,128,1,0,0,0,0,0,0,NULL,&atx[4],&avnx[4],0,&atx[16]} ,
  {0, "source" ,128,2,0,0,0,0,0,0,NULL,&atx[17],NULL,0,&atx[18]} ,
  {407, "BioSource" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[23]} ,
  {0, "gene" ,128,3,0,0,0,0,0,0,NULL,&atx[19],NULL,0,&atx[20]} ,
  {405, "Gene-ref" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[21]} ,
  {0, "prot" ,128,4,0,1,0,0,0,0,NULL,&atx[21],NULL,0,&atx[22]} ,
  {406, "Prot-ref" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[17]} ,
  {0, "rna" ,128,5,0,1,0,0,0,0,NULL,&atx[23],NULL,0,&atx[24]} ,
  {408, "RNA-ref" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[8]} ,
  {0, "summary" ,128,6,0,1,0,0,0,0,NULL,&atx[25],NULL,0,&atx[26]} ,
  {323, "VisibleString" ,0,26,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "location" ,128,7,0,1,0,0,0,0,NULL,&atx[9],&atx[27],0,&atx[35]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[28],NULL,0,NULL} ,
  {413, "Maps" ,1,0,0,0,0,0,0,0,NULL,&atx[14],&atx[29],0,&atx[36]} ,
  {0, "display-str" ,128,0,0,0,0,0,0,0,NULL,&atx[25],NULL,0,&atx[30]} ,
  {0, "method" ,128,1,0,0,0,0,0,0,NULL,&atx[34],&atx[31],0,NULL} ,
  {0, "proxy" ,128,0,0,0,0,0,0,0,NULL,&atx[25],NULL,0,&atx[32]} ,
  {0, "map-type" ,128,1,0,0,0,0,0,0,NULL,&atx[33],&avnx[17],0,NULL} ,
  {310, "ENUMERATED" ,0,10,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {315, "CHOICE" ,0,-1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "gene-source" ,128,8,0,1,0,0,0,0,NULL,&atx[36],NULL,0,&atx[45]} ,
  {414, "Gene-source" ,1,0,0,0,0,0,0,0,NULL,&atx[14],&atx[37],0,&atx[56]} ,
  {0, "src" ,128,0,0,0,0,0,0,0,NULL,&atx[25],NULL,0,&atx[38]} ,
  {0, "src-int" ,128,1,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[39]} ,
  {0, "src-str1" ,128,2,0,1,0,0,0,0,NULL,&atx[25],NULL,0,&atx[40]} ,
  {0, "src-str2" ,128,3,0,1,0,0,0,0,NULL,&atx[25],NULL,0,&atx[41]} ,
  {0, "gene-display" ,128,4,0,0,1,0,0,0,&avnx[22],&atx[42],NULL,0,&atx[43]} ,
  {301, "BOOLEAN" ,0,1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "locus-display" ,128,5,0,0,1,0,0,0,&avnx[23],&atx[42],NULL,0,&atx[44]} ,
  {0, "extra-terms" ,128,6,0,0,1,0,0,0,&avnx[24],&atx[42],NULL,0,NULL} ,
  {0, "locus" ,128,9,0,1,0,0,0,0,NULL,&atx[9],&atx[46],0,&atx[84]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[47],NULL,0,NULL} ,
  {404, "Gene-commentary" ,1,0,0,0,0,1,0,0,NULL,&atx[14],&atx[48],0,&atx[19]} ,
  {0, "type" ,128,0,0,0,0,0,0,0,NULL,&atx[4],&avnx[25],0,&atx[49]} ,
  {0, "heading" ,128,1,0,1,0,0,0,0,NULL,&atx[25],NULL,0,&atx[50]} ,
  {0, "label" ,128,2,0,1,0,0,0,0,NULL,&atx[25],NULL,0,&atx[51]} ,
  {0, "text" ,128,3,0,1,0,0,0,0,NULL,&atx[25],NULL,0,&atx[52]} ,
  {0, "accession" ,128,4,0,1,0,0,0,0,NULL,&atx[25],NULL,0,&atx[53]} ,
  {0, "version" ,128,5,0,1,0,0,0,0,NULL,&atx[4],NULL,0,&atx[54]} ,
  {0, "xtra-properties" ,128,6,0,1,0,0,0,0,NULL,&atx[9],&atx[55],0,&atx[59]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[56],NULL,0,NULL} ,
  {415, "Xtra-Terms" ,1,0,0,0,0,0,0,0,NULL,&atx[14],&atx[57],0,&atx[64]} ,
  {0, "tag" ,128,0,0,0,0,0,0,0,NULL,&atx[25],NULL,0,&atx[58]} ,
  {0, "value" ,128,1,0,0,0,0,0,0,NULL,&atx[25],NULL,0,NULL} ,
  {0, "refs" ,128,7,0,1,0,0,0,0,NULL,&atx[9],&atx[60],0,&atx[62]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[61],NULL,0,NULL} ,
  {412, "Pub" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[28]} ,
  {0, "source" ,128,8,0,1,0,0,0,0,NULL,&atx[9],&atx[63],0,&atx[70]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[64],NULL,0,NULL} ,
  {416, "Other-source" ,1,0,0,0,0,0,0,0,NULL,&atx[14],&atx[65],0,NULL} ,
  {0, "src" ,128,0,0,1,0,0,0,0,NULL,&atx[8],NULL,0,&atx[66]} ,
  {0, "pre-text" ,128,1,0,1,0,0,0,0,NULL,&atx[25],NULL,0,&atx[67]} ,
  {0, "anchor" ,128,2,0,1,0,0,0,0,NULL,&atx[25],NULL,0,&atx[68]} ,
  {0, "url" ,128,3,0,1,0,0,0,0,NULL,&atx[25],NULL,0,&atx[69]} ,
  {0, "post-text" ,128,4,0,1,0,0,0,0,NULL,&atx[25],NULL,0,NULL} ,
  {0, "genomic-coords" ,128,9,0,1,0,0,0,0,NULL,&atx[9],&atx[71],0,&atx[73]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[72],NULL,0,NULL} ,
  {411, "Seq-loc" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[61]} ,
  {0, "seqs" ,128,10,0,1,0,0,0,0,NULL,&atx[9],&atx[74],0,&atx[75]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[72],NULL,0,NULL} ,
  {0, "products" ,128,11,0,1,0,0,0,0,NULL,&atx[9],&atx[76],0,&atx[77]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[47],NULL,0,NULL} ,
  {0, "properties" ,128,12,0,1,0,0,0,0,NULL,&atx[9],&atx[78],0,&atx[79]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[47],NULL,0,NULL} ,
  {0, "comment" ,128,13,0,1,0,0,0,0,NULL,&atx[9],&atx[80],0,&atx[81]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[47],NULL,0,NULL} ,
  {0, "create-date" ,128,14,0,1,0,0,0,0,NULL,&atx[11],NULL,0,&atx[82]} ,
  {0, "update-date" ,128,15,0,1,0,0,0,0,NULL,&atx[11],NULL,0,&atx[83]} ,
  {0, "rna" ,128,16,0,1,0,0,0,0,NULL,&atx[23],NULL,0,NULL} ,
  {0, "properties" ,128,10,0,1,0,0,0,0,NULL,&atx[9],&atx[85],0,&atx[86]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[47],NULL,0,NULL} ,
  {0, "refgene" ,128,11,0,1,0,0,0,0,NULL,&atx[9],&atx[87],0,&atx[88]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[47],NULL,0,NULL} ,
  {0, "homology" ,128,12,0,1,0,0,0,0,NULL,&atx[9],&atx[89],0,&atx[90]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[47],NULL,0,NULL} ,
  {0, "comments" ,128,13,0,1,0,0,0,0,NULL,&atx[9],&atx[91],0,&atx[92]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[47],NULL,0,NULL} ,
  {0, "unique-keys" ,128,14,0,1,0,0,0,0,NULL,&atx[9],&atx[93],0,&atx[94]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[8],NULL,0,NULL} ,
  {0, "xtra-index-terms" ,128,15,0,1,0,0,0,0,NULL,&atx[9],&atx[95],0,&atx[96]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[25],NULL,0,NULL} ,
  {0, "xtra-properties" ,128,16,0,1,0,0,0,0,NULL,&atx[9],&atx[97],0,&atx[98]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[56],NULL,0,NULL} ,
  {0, "xtra-iq" ,128,17,0,1,0,0,0,0,NULL,&atx[9],&atx[99],0,&atx[100]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[56],NULL,0,NULL} ,
  {0, "non-unique-keys" ,128,18,0,1,0,0,0,0,NULL,&atx[9],&atx[101],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[8],NULL,0,NULL} ,
  {402, "Entrezgene-Set" ,1,0,0,0,0,1,0,0,NULL,&atx[104],&atx[103],0,&atx[2]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[0],NULL,0,NULL} ,
  {314, "SET OF" ,0,17,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} };

static AsnModule ampx[1] = {
  { "NCBI-Entrezgene" , "asnentgene.h36",&atx[0],NULL,NULL,0,0} };

static AsnValxNodePtr avn = avnx;
static AsnTypePtr at = atx;
static AsnModulePtr amp = ampx;



/**************************************************
*
*    Defines for Module NCBI-Entrezgene
*
**************************************************/

#define ENTREZGENE &at[0]
#define ENTREZGENE_track_info &at[1]
#define ENTREZGENE_type &at[15]
#define ENTREZGENE_source &at[16]
#define ENTREZGENE_gene &at[18]
#define ENTREZGENE_prot &at[20]
#define ENTREZGENE_rna &at[22]
#define ENTREZGENE_summary &at[24]
#define ENTREZGENE_location &at[26]
#define ENTREZGENE_location_E &at[27]
#define ENTREZGENE_gene_source &at[35]
#define ENTREZGENE_locus &at[45]
#define ENTREZGENE_locus_E &at[46]
#define ENTREZGENE_properties &at[84]
#define ENTREZGENE_properties_E &at[85]
#define ENTREZGENE_refgene &at[86]
#define ENTREZGENE_refgene_E &at[87]
#define ENTREZGENE_homology &at[88]
#define ENTREZGENE_homology_E &at[89]
#define ENTREZGENE_comments &at[90]
#define ENTREZGENE_comments_E &at[91]
#define ENTREZGENE_unique_keys &at[92]
#define ENTREZGENE_unique_keys_E &at[93]
#define ENTREZGENE_xtra_index_terms &at[94]
#define ENTREZGENE_xtra_index_terms_E &at[95]
#define ENTREZGENE_xtra_properties &at[96]
#define ENTREZGENE_xtra_properties_E &at[97]
#define ENTREZGENE_xtra_iq &at[98]
#define ENTREZGENE_xtra_iq_E &at[99]
#define ENTREZGENE_non_unique_keys &at[100]
#define ENTREZGENE_non_unique_keys_E &at[101]

#define ENTREZGENE_SET &at[102]
#define ENTREZGENE_SET_E &at[103]

#define GENE_TRACK &at[2]
#define GENE_TRACK_geneid &at[3]
#define GENE_TRACK_status &at[5]
#define GENE_TRACK_current_id &at[6]
#define GENE_TRACK_current_id_E &at[7]
#define GENE_TRACK_create_date &at[10]
#define GENE_TRACK_update_date &at[12]
#define GENE_TRACK_discontinue_date &at[13]

#define GENE_COMMENTARY &at[47]
#define GENE_COMMENTARY_type &at[48]
#define GENE_COMMENTARY_heading &at[49]
#define GENE_COMMENTARY_label &at[50]
#define GENE_COMMENTARY_text &at[51]
#define GENE_COMMENTARY_accession &at[52]
#define GENE_COMMENTARY_version &at[53]
#define GENE_COMMENTARY_xtra_properties &at[54]
#define COMMENTARY_xtra_properties_E &at[55]
#define GENE_COMMENTARY_refs &at[59]
#define GENE_COMMENTARY_refs_E &at[60]
#define GENE_COMMENTARY_source &at[62]
#define GENE_COMMENTARY_source_E &at[63]
#define GENE_COMMENTARY_genomic_coords &at[70]
#define COMMENTARY_genomic_coords_E &at[71]
#define GENE_COMMENTARY_seqs &at[73]
#define GENE_COMMENTARY_seqs_E &at[74]
#define GENE_COMMENTARY_products &at[75]
#define GENE_COMMENTARY_products_E &at[76]
#define GENE_COMMENTARY_properties &at[77]
#define GENE_COMMENTARY_properties_E &at[78]
#define GENE_COMMENTARY_comment &at[79]
#define GENE_COMMENTARY_comment_E &at[80]
#define GENE_COMMENTARY_create_date &at[81]
#define GENE_COMMENTARY_update_date &at[82]
#define GENE_COMMENTARY_rna &at[83]

#define MAPS &at[28]
#define MAPS_display_str &at[29]
#define MAPS_method &at[30]
#define MAPS_method_proxy &at[31]
#define MAPS_method_map_type &at[32]

#define GENE_SOURCE &at[36]
#define GENE_SOURCE_src &at[37]
#define GENE_SOURCE_src_int &at[38]
#define GENE_SOURCE_src_str1 &at[39]
#define GENE_SOURCE_src_str2 &at[40]
#define GENE_SOURCE_gene_display &at[41]
#define GENE_SOURCE_locus_display &at[43]
#define GENE_SOURCE_extra_terms &at[44]

#define XTRA_TERMS &at[56]
#define XTRA_TERMS_tag &at[57]
#define XTRA_TERMS_value &at[58]

#define OTHER_SOURCE &at[64]
#define OTHER_SOURCE_src &at[65]
#define OTHER_SOURCE_pre_text &at[66]
#define OTHER_SOURCE_anchor &at[67]
#define OTHER_SOURCE_url &at[68]
#define OTHER_SOURCE_post_text &at[69]
