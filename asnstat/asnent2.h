/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "asnent2.h15";
static AsnValxNode avnx[10] = {
    {20,"and" ,1,0.0,&avnx[1] } ,
    {20,"or" ,2,0.0,&avnx[2] } ,
    {20,"butnot" ,3,0.0,&avnx[3] } ,
    {20,"range" ,4,0.0,&avnx[4] } ,
    {20,"left-paren" ,5,0.0,&avnx[5] } ,
    {20,"right-paren" ,6,0.0,NULL } ,
    {2,NULL,0,0.0,NULL } ,
    {2,NULL,0,0.0,NULL } ,
    {2,NULL,0,0.0,NULL } ,
    {2,NULL,0,0.0,NULL } };

static AsnType atx[188] = {
  {401, "Entrez2-dt" ,1,0,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[2]} ,
  {302, "INTEGER" ,0,2,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {402, "Entrez2-db-id" ,1,0,0,0,0,0,0,0,NULL,&atx[3],NULL,0,&atx[4]} ,
  {323, "VisibleString" ,0,26,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {403, "Entrez2-field-id" ,1,0,0,0,0,0,0,0,NULL,&atx[3],NULL,0,&atx[5]} ,
  {404, "Entrez2-link-id" ,1,0,0,0,0,0,0,0,NULL,&atx[3],NULL,0,&atx[6]} ,
  {405, "Entrez2-id-list" ,1,0,0,0,0,0,0,0,NULL,&atx[11],&atx[7],0,&atx[12]} ,
  {0, "db" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[8]} ,
  {0, "num" ,128,1,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[9]} ,
  {0, "uids" ,128,2,0,1,0,0,0,0,NULL,&atx[10],NULL,0,NULL} ,
  {304, "OCTET STRING" ,0,4,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {311, "SEQUENCE" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {406, "Entrez2-boolean-exp" ,1,0,0,0,0,0,0,0,NULL,&atx[11],&atx[13],0,&atx[16]} ,
  {0, "db" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[14]} ,
  {0, "exp" ,128,1,0,0,0,0,0,0,NULL,&atx[30],&atx[15],0,&atx[31]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[16],NULL,0,NULL} ,
  {407, "Entrez2-boolean-element" ,1,0,0,0,0,0,0,0,NULL,&atx[29],&atx[17],0,&atx[32]} ,
  {0, "str" ,128,0,0,0,0,0,0,0,NULL,&atx[3],NULL,0,&atx[18]} ,
  {0, "op" ,128,1,0,0,0,0,0,0,NULL,&atx[19],NULL,0,&atx[20]} ,
  {409, "Entrez2-operator" ,1,0,0,0,0,0,0,0,NULL,&atx[1],&avnx[0],0,&atx[21]} ,
  {0, "term" ,128,2,0,0,0,0,0,0,NULL,&atx[21],NULL,0,&atx[28]} ,
  {410, "Entrez2-boolean-term" ,1,0,0,0,0,0,0,0,NULL,&atx[11],&atx[22],0,&atx[40]} ,
  {0, "field" ,128,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[23]} ,
  {0, "term" ,128,1,0,0,0,0,0,0,NULL,&atx[3],NULL,0,&atx[24]} ,
  {0, "term-count" ,128,2,0,1,0,0,0,0,NULL,&atx[1],NULL,0,&atx[25]} ,
  {0, "do-not-explode" ,128,3,0,0,1,0,0,0,&avnx[6],&atx[26],NULL,0,&atx[27]} ,
  {301, "BOOLEAN" ,0,1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "do-not-translate" ,128,4,0,0,1,0,0,0,&avnx[7],&atx[26],NULL,0,NULL} ,
  {0, "ids" ,128,3,0,0,0,0,0,0,NULL,&atx[6],NULL,0,NULL} ,
  {315, "CHOICE" ,0,-1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {312, "SEQUENCE OF" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "limits" ,128,2,0,1,0,0,0,0,NULL,&atx[32],NULL,0,NULL} ,
  {408, "Entrez2-limits" ,1,0,0,0,0,0,0,0,NULL,&atx[11],&atx[33],0,&atx[19]} ,
  {0, "filter-date" ,128,0,0,1,0,0,0,0,NULL,&atx[34],NULL,0,&atx[38]} ,
  {419, "Entrez2-dt-filter" ,1,0,0,0,0,0,0,0,NULL,&atx[11],&atx[35],0,&atx[82]} ,
  {0, "begin-date" ,128,0,0,0,0,0,0,0,NULL,&atx[0],NULL,0,&atx[36]} ,
  {0, "end-date" ,128,1,0,0,0,0,0,0,NULL,&atx[0],NULL,0,&atx[37]} ,
  {0, "type-date" ,128,2,0,0,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {0, "max-UIDs" ,128,1,0,1,0,0,0,0,NULL,&atx[1],NULL,0,&atx[39]} ,
  {0, "offset-UIDs" ,128,2,0,1,0,0,0,0,NULL,&atx[1],NULL,0,NULL} ,
  {411, "Entrez2-request" ,1,0,0,0,0,0,0,0,NULL,&atx[11],&atx[41],0,&atx[42]} ,
  {0, "request" ,128,0,0,0,0,0,0,0,NULL,&atx[42],NULL,0,&atx[80]} ,
  {412, "E2Request" ,1,0,0,0,0,0,0,0,NULL,&atx[29],&atx[43],0,&atx[46]} ,
  {0, "get-info" ,128,0,0,0,0,0,0,0,NULL,&atx[44],NULL,0,&atx[45]} ,
  {305, "NULL" ,0,5,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "eval-boolean" ,128,1,0,0,0,0,0,0,NULL,&atx[46],NULL,0,&atx[50]} ,
  {413, "Entrez2-eval-boolean" ,1,0,0,0,0,0,0,0,NULL,&atx[11],&atx[47],0,&atx[52]} ,
  {0, "return-UIDs" ,128,0,0,0,1,0,0,0,&avnx[8],&atx[26],NULL,0,&atx[48]} ,
  {0, "return-parse" ,128,1,0,0,1,0,0,0,&avnx[9],&atx[26],NULL,0,&atx[49]} ,
  {0, "query" ,128,2,0,0,0,0,0,0,NULL,&atx[12],NULL,0,NULL} ,
  {0, "get-docsum" ,128,2,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[51]} ,
  {0, "get-term-pos" ,128,3,0,0,0,0,0,0,NULL,&atx[52],NULL,0,&atx[56]} ,
  {414, "Entrez2-term-query" ,1,0,0,0,0,0,0,0,NULL,&atx[11],&atx[53],0,&atx[57]} ,
  {0, "db" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[54]} ,
  {0, "field" ,128,1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[55]} ,
  {0, "term" ,128,2,0,0,0,0,0,0,NULL,&atx[3],NULL,0,NULL} ,
  {0, "get-term-list" ,128,4,0,0,0,0,0,0,NULL,&atx[57],NULL,0,&atx[62]} ,
  {415, "Entrez2-term-pos" ,1,0,0,0,0,0,0,0,NULL,&atx[11],&atx[58],0,&atx[63]} ,
  {0, "db" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[59]} ,
  {0, "field" ,128,1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[60]} ,
  {0, "first-term-pos" ,128,2,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[61]} ,
  {0, "number-of-terms" ,128,3,0,1,0,0,0,0,NULL,&atx[1],NULL,0,NULL} ,
  {0, "get-term-hierarchy" ,128,5,0,0,0,0,0,0,NULL,&atx[63],NULL,0,&atx[68]} ,
  {416, "Entrez2-hier-query" ,1,0,0,0,0,0,0,0,NULL,&atx[11],&atx[64],0,&atx[69]} ,
  {0, "db" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[65]} ,
  {0, "field" ,128,1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[66]} ,
  {0, "term" ,128,2,0,1,0,0,0,0,NULL,&atx[3],NULL,0,&atx[67]} ,
  {0, "txid" ,128,3,0,1,0,0,0,0,NULL,&atx[1],NULL,0,NULL} ,
  {0, "get-links" ,128,6,0,0,0,0,0,0,NULL,&atx[69],NULL,0,&atx[75]} ,
  {417, "Entrez2-get-links" ,1,0,0,0,0,0,0,0,NULL,&atx[11],&atx[70],0,&atx[77]} ,
  {0, "uids" ,128,0,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[71]} ,
  {0, "linktype" ,128,1,0,0,0,0,0,0,NULL,&atx[5],NULL,0,&atx[72]} ,
  {0, "max-UIDS" ,128,2,0,1,0,0,0,0,NULL,&atx[1],NULL,0,&atx[73]} ,
  {0, "count-only" ,128,3,0,1,0,0,0,0,NULL,&atx[26],NULL,0,&atx[74]} ,
  {0, "parents-persist" ,128,4,0,1,0,0,0,0,NULL,&atx[26],NULL,0,NULL} ,
  {0, "get-linked" ,128,7,0,0,0,0,0,0,NULL,&atx[69],NULL,0,&atx[76]} ,
  {0, "get-link-counts" ,128,8,0,0,0,0,0,0,NULL,&atx[77],NULL,0,NULL} ,
  {418, "Entrez2-id" ,1,0,0,0,0,0,0,0,NULL,&atx[11],&atx[78],0,&atx[34]} ,
  {0, "db" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[79]} ,
  {0, "uid" ,128,1,0,0,0,0,0,0,NULL,&atx[1],NULL,0,NULL} ,
  {0, "version" ,128,1,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[81]} ,
  {0, "tool" ,128,2,0,1,0,0,0,0,NULL,&atx[3],NULL,0,NULL} ,
  {420, "Entrez2-reply" ,1,0,0,0,0,0,0,0,NULL,&atx[11],&atx[83],0,&atx[84]} ,
  {0, "reply" ,128,0,0,0,0,0,0,0,NULL,&atx[84],NULL,0,&atx[185]} ,
  {421, "E2Reply" ,1,0,0,0,0,0,0,0,NULL,&atx[29],&atx[85],0,&atx[87]} ,
  {0, "error" ,128,0,0,0,0,0,0,0,NULL,&atx[3],NULL,0,&atx[86]} ,
  {0, "get-info" ,128,1,0,0,0,0,0,0,NULL,&atx[87],NULL,0,&atx[118]} ,
  {422, "Entrez2-info" ,1,0,0,0,0,0,0,0,NULL,&atx[11],&atx[88],0,&atx[119]} ,
  {0, "db-count" ,128,0,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[89]} ,
  {0, "build-date" ,128,1,0,0,0,0,0,0,NULL,&atx[0],NULL,0,&atx[90]} ,
  {0, "db-info" ,128,2,0,0,0,0,0,0,NULL,&atx[30],&atx[91],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[92],NULL,0,NULL} ,
  {429, "Entrez2-db-info" ,1,0,0,0,0,0,0,0,NULL,&atx[11],&atx[93],0,&atx[100]} ,
  {0, "db-name" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[94]} ,
  {0, "db-menu" ,128,1,0,0,0,0,0,0,NULL,&atx[3],NULL,0,&atx[95]} ,
  {0, "db-descr" ,128,2,0,0,0,0,0,0,NULL,&atx[3],NULL,0,&atx[96]} ,
  {0, "doc-count" ,128,3,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[97]} ,
  {0, "field-count" ,128,4,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[98]} ,
  {0, "fields" ,128,5,0,0,0,0,0,0,NULL,&atx[30],&atx[99],0,&atx[109]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[100],NULL,0,NULL} ,
  {430, "Entrez2-field-info" ,1,0,0,0,0,0,0,0,NULL,&atx[11],&atx[101],0,&atx[112]} ,
  {0, "field-name" ,128,0,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[102]} ,
  {0, "field-menu" ,128,1,0,0,0,0,0,0,NULL,&atx[3],NULL,0,&atx[103]} ,
  {0, "field-descr" ,128,2,0,0,0,0,0,0,NULL,&atx[3],NULL,0,&atx[104]} ,
  {0, "term-count" ,128,3,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[105]} ,
  {0, "is-date" ,128,4,0,1,0,0,0,0,NULL,&atx[26],NULL,0,&atx[106]} ,
  {0, "is-numerical" ,128,5,0,1,0,0,0,0,NULL,&atx[26],NULL,0,&atx[107]} ,
  {0, "single-token" ,128,6,0,1,0,0,0,0,NULL,&atx[26],NULL,0,&atx[108]} ,
  {0, "hierarchy-avail" ,128,7,0,1,0,0,0,0,NULL,&atx[26],NULL,0,NULL} ,
  {0, "link-count" ,128,6,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[110]} ,
  {0, "links" ,128,7,0,0,0,0,0,0,NULL,&atx[30],&atx[111],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[112],NULL,0,NULL} ,
  {431, "Entrez2-link-info" ,1,0,0,0,0,0,0,0,NULL,&atx[11],&atx[113],0,&atx[128]} ,
  {0, "link-name" ,128,0,0,0,0,0,0,0,NULL,&atx[5],NULL,0,&atx[114]} ,
  {0, "link-menu" ,128,1,0,0,0,0,0,0,NULL,&atx[3],NULL,0,&atx[115]} ,
  {0, "link-descr" ,128,2,0,0,0,0,0,0,NULL,&atx[3],NULL,0,&atx[116]} ,
  {0, "db-to" ,128,3,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[117]} ,
  {0, "data-size" ,128,4,0,1,0,0,0,0,NULL,&atx[1],NULL,0,NULL} ,
  {0, "eval-boolean" ,128,2,0,0,0,0,0,0,NULL,&atx[119],NULL,0,&atx[123]} ,
  {423, "Entrez2-boolean-reply" ,1,0,0,0,0,0,0,0,NULL,&atx[11],&atx[120],0,&atx[124]} ,
  {0, "count" ,128,0,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[121]} ,
  {0, "uids" ,128,1,0,1,0,0,0,0,NULL,&atx[6],NULL,0,&atx[122]} ,
  {0, "query" ,128,2,0,1,0,0,0,0,NULL,&atx[12],NULL,0,NULL} ,
  {0, "get-docsum" ,128,3,0,0,0,0,0,0,NULL,&atx[124],NULL,0,&atx[149]} ,
  {424, "Entrez2-docsum-list" ,1,0,0,0,0,0,0,0,NULL,&atx[11],&atx[125],0,&atx[151]} ,
  {0, "count" ,128,0,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[126]} ,
  {0, "list" ,128,1,0,0,0,0,0,0,NULL,&atx[30],&atx[127],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[128],NULL,0,NULL} ,
  {432, "Entrez2-docsum" ,1,0,0,0,0,0,0,0,NULL,&atx[11],&atx[129],0,&atx[156]} ,
  {0, "uid" ,128,0,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[130]} ,
  {0, "secondary-uid" ,128,1,0,1,0,0,0,0,NULL,&atx[1],NULL,0,&atx[131]} ,
  {0, "caption" ,128,2,0,1,0,0,0,0,NULL,&atx[3],NULL,0,&atx[132]} ,
  {0, "title" ,128,3,0,1,0,0,0,0,NULL,&atx[3],NULL,0,&atx[133]} ,
  {0, "short-citation" ,128,4,0,1,0,0,0,0,NULL,&atx[3],NULL,0,&atx[134]} ,
  {0, "entrez-date" ,128,5,0,0,0,0,0,0,NULL,&atx[0],NULL,0,&atx[135]} ,
  {0, "language" ,128,6,0,1,0,0,0,0,NULL,&atx[3],NULL,0,&atx[136]} ,
  {0, "create-date" ,128,7,0,1,0,0,0,0,NULL,&atx[3],NULL,0,&atx[137]} ,
  {0, "update-date" ,128,8,0,1,0,0,0,0,NULL,&atx[3],NULL,0,&atx[138]} ,
  {0, "seqlen" ,128,9,0,1,0,0,0,0,NULL,&atx[1],NULL,0,&atx[139]} ,
  {0, "author" ,128,10,0,1,0,0,0,0,NULL,&atx[3],NULL,0,&atx[140]} ,
  {0, "source" ,128,11,0,1,0,0,0,0,NULL,&atx[3],NULL,0,&atx[141]} ,
  {0, "volume" ,128,12,0,1,0,0,0,0,NULL,&atx[3],NULL,0,&atx[142]} ,
  {0, "pages" ,128,13,0,1,0,0,0,0,NULL,&atx[3],NULL,0,&atx[143]} ,
  {0, "pub-type" ,128,14,0,1,0,0,0,0,NULL,&atx[3],NULL,0,&atx[144]} ,
  {0, "record-status" ,128,15,0,1,0,0,0,0,NULL,&atx[3],NULL,0,&atx[145]} ,
  {0, "no-abstract" ,128,16,0,1,0,0,0,0,NULL,&atx[26],NULL,0,&atx[146]} ,
  {0, "translated-title" ,128,17,0,1,0,0,0,0,NULL,&atx[26],NULL,0,&atx[147]} ,
  {0, "no-authors" ,128,18,0,1,0,0,0,0,NULL,&atx[26],NULL,0,&atx[148]} ,
  {0, "taxid" ,128,19,0,1,0,0,0,0,NULL,&atx[1],NULL,0,NULL} ,
  {0, "get-term-pos" ,128,4,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[150]} ,
  {0, "get-term-list" ,128,5,0,0,0,0,0,0,NULL,&atx[151],NULL,0,&atx[161]} ,
  {425, "Entrez2-term-list" ,1,0,0,0,0,0,0,0,NULL,&atx[11],&atx[152],0,&atx[162]} ,
  {0, "pos" ,128,0,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[153]} ,
  {0, "num" ,128,1,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[154]} ,
  {0, "list" ,128,2,0,0,0,0,0,0,NULL,&atx[30],&atx[155],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[156],NULL,0,NULL} ,
  {433, "Entrez2-term" ,1,0,0,0,0,0,0,0,NULL,&atx[11],&atx[157],0,&atx[182]} ,
  {0, "term" ,128,0,0,0,0,0,0,0,NULL,&atx[3],NULL,0,&atx[158]} ,
  {0, "txid" ,128,1,0,1,0,0,0,0,NULL,&atx[1],NULL,0,&atx[159]} ,
  {0, "count" ,128,2,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[160]} ,
  {0, "is-leaf-node" ,128,3,0,1,0,0,0,0,NULL,&atx[26],NULL,0,NULL} ,
  {0, "get-term-hierarchy" ,128,6,0,0,0,0,0,0,NULL,&atx[162],NULL,0,&atx[171]} ,
  {426, "Entrez2-hier-node" ,1,0,0,0,0,0,0,0,NULL,&atx[11],&atx[163],0,&atx[172]} ,
  {0, "cannonical-form" ,128,0,0,0,0,0,0,0,NULL,&atx[3],NULL,0,&atx[164]} ,
  {0, "lineage-count" ,128,1,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[165]} ,
  {0, "lineage" ,128,2,0,1,0,0,0,0,NULL,&atx[30],&atx[166],0,&atx[167]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[156],NULL,0,NULL} ,
  {0, "child-count" ,128,3,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[168]} ,
  {0, "children" ,128,4,0,0,0,0,0,0,NULL,&atx[30],&atx[169],0,&atx[170]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[156],NULL,0,NULL} ,
  {0, "is-ambiguous" ,128,5,0,1,0,0,0,0,NULL,&atx[26],NULL,0,NULL} ,
  {0, "get-links" ,128,7,0,0,0,0,0,0,NULL,&atx[172],NULL,0,&atx[176]} ,
  {427, "Entrez2-link-set" ,1,0,0,0,0,0,0,0,NULL,&atx[11],&atx[173],0,&atx[178]} ,
  {0, "ids" ,128,0,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[174]} ,
  {0, "data-size" ,128,1,0,1,0,0,0,0,NULL,&atx[1],NULL,0,&atx[175]} ,
  {0, "data" ,128,2,0,1,0,0,0,0,NULL,&atx[10],NULL,0,NULL} ,
  {0, "get-linked" ,128,8,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[177]} ,
  {0, "get-link-counts" ,128,9,0,0,0,0,0,0,NULL,&atx[178],NULL,0,NULL} ,
  {428, "Entrez2-link-count-list" ,1,0,0,0,0,0,0,0,NULL,&atx[11],&atx[179],0,&atx[92]} ,
  {0, "link-type-count" ,128,0,0,0,0,0,0,0,NULL,&atx[1],NULL,0,&atx[180]} ,
  {0, "links" ,128,1,0,0,0,0,0,0,NULL,&atx[30],&atx[181],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[182],NULL,0,NULL} ,
  {434, "Entrez2-link-count" ,1,0,0,0,0,0,0,0,NULL,&atx[11],&atx[183],0,NULL} ,
  {0, "link-type" ,128,0,0,0,0,0,0,0,NULL,&atx[5],NULL,0,&atx[184]} ,
  {0, "link-count" ,128,1,0,0,0,0,0,0,NULL,&atx[1],NULL,0,NULL} ,
  {0, "dt" ,128,1,0,0,0,0,0,0,NULL,&atx[0],NULL,0,&atx[186]} ,
  {0, "server" ,128,2,0,0,0,0,0,0,NULL,&atx[3],NULL,0,&atx[187]} ,
  {0, "msg" ,128,3,0,1,0,0,0,0,NULL,&atx[3],NULL,0,NULL} };

static AsnModule ampx[1] = {
  { "NCBI-Entrez2" , "asnent2.h15",&atx[0],NULL,NULL,0,0} };

static AsnValxNodePtr avn = avnx;
static AsnTypePtr at = atx;
static AsnModulePtr amp = ampx;



/**************************************************
*
*    Defines for Module NCBI-Entrez2
*
**************************************************/

#define ENTREZ2_DT &at[0]

#define ENTREZ2_DB_ID &at[2]

#define ENTREZ2_FIELD_ID &at[4]

#define ENTREZ2_LINK_ID &at[5]

#define ENTREZ2_ID_LIST &at[6]
#define ENTREZ2_ID_LIST_db &at[7]
#define ENTREZ2_ID_LIST_num &at[8]
#define ENTREZ2_ID_LIST_uids &at[9]

#define ENTREZ2_BOOLEAN_EXP &at[12]
#define ENTREZ2_BOOLEAN_EXP_db &at[13]
#define ENTREZ2_BOOLEAN_EXP_exp &at[14]
#define ENTREZ2_BOOLEAN_EXP_exp_E &at[15]
#define ENTREZ2_BOOLEAN_EXP_limits &at[31]

#define ENTREZ2_BOOLEAN_ELEMENT &at[16]
#define ENTREZ2_BOOLEAN_ELEMENT_str &at[17]
#define ENTREZ2_BOOLEAN_ELEMENT_op &at[18]
#define ENTREZ2_BOOLEAN_ELEMENT_term &at[20]
#define ENTREZ2_BOOLEAN_ELEMENT_ids &at[28]

#define ENTREZ2_LIMITS &at[32]
#define ENTREZ2_LIMITS_filter_date &at[33]
#define ENTREZ2_LIMITS_max_UIDs &at[38]
#define ENTREZ2_LIMITS_offset_UIDs &at[39]

#define ENTREZ2_OPERATOR &at[19]

#define ENTREZ2_BOOLEAN_TERM &at[21]
#define ENTREZ2_BOOLEAN_TERM_field &at[22]
#define ENTREZ2_BOOLEAN_TERM_term &at[23]
#define ENTREZ2_BOOLEAN_TERM_term_count &at[24]
#define BOOLEAN_TERM_do_not_explode &at[25]
#define BOOLEAN_TERM_do_not_translate &at[27]

#define ENTREZ2_REQUEST &at[40]
#define ENTREZ2_REQUEST_request &at[41]
#define ENTREZ2_REQUEST_version &at[80]
#define ENTREZ2_REQUEST_tool &at[81]

#define E2REQUEST &at[42]
#define E2REQUEST_get_info &at[43]
#define E2REQUEST_eval_boolean &at[45]
#define E2REQUEST_get_docsum &at[50]
#define E2REQUEST_get_term_pos &at[51]
#define E2REQUEST_get_term_list &at[56]
#define E2REQUEST_get_term_hierarchy &at[62]
#define E2REQUEST_get_links &at[68]
#define E2REQUEST_get_linked &at[75]
#define E2REQUEST_get_link_counts &at[76]

#define ENTREZ2_EVAL_BOOLEAN &at[46]
#define EVAL_BOOLEAN_return_UIDs &at[47]
#define EVAL_BOOLEAN_return_parse &at[48]
#define ENTREZ2_EVAL_BOOLEAN_query &at[49]

#define ENTREZ2_TERM_QUERY &at[52]
#define ENTREZ2_TERM_QUERY_db &at[53]
#define ENTREZ2_TERM_QUERY_field &at[54]
#define ENTREZ2_TERM_QUERY_term &at[55]

#define ENTREZ2_TERM_POS &at[57]
#define ENTREZ2_TERM_POS_db &at[58]
#define ENTREZ2_TERM_POS_field &at[59]
#define ENTREZ2_TERM_POS_first_term_pos &at[60]
#define TERM_POS_number_of_terms &at[61]

#define ENTREZ2_HIER_QUERY &at[63]
#define ENTREZ2_HIER_QUERY_db &at[64]
#define ENTREZ2_HIER_QUERY_field &at[65]
#define ENTREZ2_HIER_QUERY_term &at[66]
#define ENTREZ2_HIER_QUERY_txid &at[67]

#define ENTREZ2_GET_LINKS &at[69]
#define ENTREZ2_GET_LINKS_uids &at[70]
#define ENTREZ2_GET_LINKS_linktype &at[71]
#define ENTREZ2_GET_LINKS_max_UIDS &at[72]
#define ENTREZ2_GET_LINKS_count_only &at[73]
#define GET_LINKS_parents_persist &at[74]

#define ENTREZ2_ID &at[77]
#define ENTREZ2_ID_db &at[78]
#define ENTREZ2_ID_uid &at[79]

#define ENTREZ2_DT_FILTER &at[34]
#define ENTREZ2_DT_FILTER_begin_date &at[35]
#define ENTREZ2_DT_FILTER_end_date &at[36]
#define ENTREZ2_DT_FILTER_type_date &at[37]

#define ENTREZ2_REPLY &at[82]
#define ENTREZ2_REPLY_reply &at[83]
#define ENTREZ2_REPLY_dt &at[185]
#define ENTREZ2_REPLY_server &at[186]
#define ENTREZ2_REPLY_msg &at[187]

#define E2REPLY &at[84]
#define E2REPLY_error &at[85]
#define E2REPLY_get_info &at[86]
#define E2REPLY_eval_boolean &at[118]
#define E2REPLY_get_docsum &at[123]
#define E2REPLY_get_term_pos &at[149]
#define E2REPLY_get_term_list &at[150]
#define E2REPLY_get_term_hierarchy &at[161]
#define E2REPLY_get_links &at[171]
#define E2REPLY_get_linked &at[176]
#define E2REPLY_get_link_counts &at[177]

#define ENTREZ2_INFO &at[87]
#define ENTREZ2_INFO_db_count &at[88]
#define ENTREZ2_INFO_build_date &at[89]
#define ENTREZ2_INFO_db_info &at[90]
#define ENTREZ2_INFO_db_info_E &at[91]

#define ENTREZ2_BOOLEAN_REPLY &at[119]
#define ENTREZ2_BOOLEAN_REPLY_count &at[120]
#define ENTREZ2_BOOLEAN_REPLY_uids &at[121]
#define ENTREZ2_BOOLEAN_REPLY_query &at[122]

#define ENTREZ2_DOCSUM_LIST &at[124]
#define ENTREZ2_DOCSUM_LIST_count &at[125]
#define ENTREZ2_DOCSUM_LIST_list &at[126]
#define ENTREZ2_DOCSUM_LIST_list_E &at[127]

#define ENTREZ2_TERM_LIST &at[151]
#define ENTREZ2_TERM_LIST_pos &at[152]
#define ENTREZ2_TERM_LIST_num &at[153]
#define ENTREZ2_TERM_LIST_list &at[154]
#define ENTREZ2_TERM_LIST_list_E &at[155]

#define ENTREZ2_HIER_NODE &at[162]
#define HIER_NODE_cannonical_form &at[163]
#define ENTREZ2_HIER_NODE_lineage_count &at[164]
#define ENTREZ2_HIER_NODE_lineage &at[165]
#define ENTREZ2_HIER_NODE_lineage_E &at[166]
#define ENTREZ2_HIER_NODE_child_count &at[167]
#define ENTREZ2_HIER_NODE_children &at[168]
#define ENTREZ2_HIER_NODE_children_E &at[169]
#define ENTREZ2_HIER_NODE_is_ambiguous &at[170]

#define ENTREZ2_LINK_SET &at[172]
#define ENTREZ2_LINK_SET_ids &at[173]
#define ENTREZ2_LINK_SET_data_size &at[174]
#define ENTREZ2_LINK_SET_data &at[175]

#define ENTREZ2_LINK_COUNT_LIST &at[178]
#define COUNT_LIST_link_type_count &at[179]
#define ENTREZ2_LINK_COUNT_LIST_links &at[180]
#define ENTREZ2_LINK_COUNT_LIST_links_E &at[181]

#define ENTREZ2_DB_INFO &at[92]
#define ENTREZ2_DB_INFO_db_name &at[93]
#define ENTREZ2_DB_INFO_db_menu &at[94]
#define ENTREZ2_DB_INFO_db_descr &at[95]
#define ENTREZ2_DB_INFO_doc_count &at[96]
#define ENTREZ2_DB_INFO_field_count &at[97]
#define ENTREZ2_DB_INFO_fields &at[98]
#define ENTREZ2_DB_INFO_fields_E &at[99]
#define ENTREZ2_DB_INFO_link_count &at[109]
#define ENTREZ2_DB_INFO_links &at[110]
#define ENTREZ2_DB_INFO_links_E &at[111]

#define ENTREZ2_FIELD_INFO &at[100]
#define ENTREZ2_FIELD_INFO_field_name &at[101]
#define ENTREZ2_FIELD_INFO_field_menu &at[102]
#define ENTREZ2_FIELD_INFO_field_descr &at[103]
#define ENTREZ2_FIELD_INFO_term_count &at[104]
#define ENTREZ2_FIELD_INFO_is_date &at[105]
#define ENTREZ2_FIELD_INFO_is_numerical &at[106]
#define ENTREZ2_FIELD_INFO_single_token &at[107]
#define FIELD_INFO_hierarchy_avail &at[108]

#define ENTREZ2_LINK_INFO &at[112]
#define ENTREZ2_LINK_INFO_link_name &at[113]
#define ENTREZ2_LINK_INFO_link_menu &at[114]
#define ENTREZ2_LINK_INFO_link_descr &at[115]
#define ENTREZ2_LINK_INFO_db_to &at[116]
#define ENTREZ2_LINK_INFO_data_size &at[117]

#define ENTREZ2_DOCSUM &at[128]
#define ENTREZ2_DOCSUM_uid &at[129]
#define ENTREZ2_DOCSUM_secondary_uid &at[130]
#define ENTREZ2_DOCSUM_caption &at[131]
#define ENTREZ2_DOCSUM_title &at[132]
#define ENTREZ2_DOCSUM_short_citation &at[133]
#define ENTREZ2_DOCSUM_entrez_date &at[134]
#define ENTREZ2_DOCSUM_language &at[135]
#define ENTREZ2_DOCSUM_create_date &at[136]
#define ENTREZ2_DOCSUM_update_date &at[137]
#define ENTREZ2_DOCSUM_seqlen &at[138]
#define ENTREZ2_DOCSUM_author &at[139]
#define ENTREZ2_DOCSUM_source &at[140]
#define ENTREZ2_DOCSUM_volume &at[141]
#define ENTREZ2_DOCSUM_pages &at[142]
#define ENTREZ2_DOCSUM_pub_type &at[143]
#define ENTREZ2_DOCSUM_record_status &at[144]
#define ENTREZ2_DOCSUM_no_abstract &at[145]
#define ENTREZ2_DOCSUM_translated_title &at[146]
#define ENTREZ2_DOCSUM_no_authors &at[147]
#define ENTREZ2_DOCSUM_taxid &at[148]

#define ENTREZ2_TERM &at[156]
#define ENTREZ2_TERM_term &at[157]
#define ENTREZ2_TERM_txid &at[158]
#define ENTREZ2_TERM_count &at[159]
#define ENTREZ2_TERM_is_leaf_node &at[160]

#define ENTREZ2_LINK_COUNT &at[182]
#define ENTREZ2_LINK_COUNT_link_type &at[183]
#define ENTREZ2_LINK_COUNT_link_count &at[184]
