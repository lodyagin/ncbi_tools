/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "asnbl18.h60";
static AsnValxNode avnx[20] = {
    {20,"not-set" ,0,0.0,&avnx[1] } ,
    {20,"amino-acid" ,1,0.0,&avnx[2] } ,
    {20,"nucleic-acid" ,2,0.0,&avnx[3] } ,
    {20,"other" ,255,0.0,NULL } ,
    {20,"not-set" ,0,0.0,&avnx[5] } ,
    {20,"neighborhood" ,1,0.0,&avnx[6] } ,
    {20,"search" ,2,0.0,&avnx[7] } ,
    {20,"threecomps" ,3,0.0,NULL } ,
    {20,"plus" ,1,0.0,&avnx[9] } ,
    {20,"minus" ,2,0.0,&avnx[10] } ,
    {20,"both" ,3,0.0,&avnx[11] } ,
    {20,"plus-rf" ,5,0.0,&avnx[12] } ,
    {20,"minus-rf" ,6,0.0,NULL } ,
    {20,"ncbi4na" ,4,0.0,&avnx[14] } ,
    {20,"ncbistdaa" ,11,0.0,NULL } ,
    {2,NULL,1,0.0,NULL } ,
    {2,NULL,1,0.0,NULL } ,
    {2,NULL,1,0.0,NULL } ,
    {2,NULL,1,0.0,NULL } ,
    {2,NULL,1,0.0,NULL } };

static AsnType atx[210] = {
  {401, "BLAST0-Preface" ,1,0,0,0,0,1,0,0,NULL,&atx[20],&atx[1],0,&atx[22]} ,
  {0, "program" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[3]} ,
  {323, "VisibleString" ,0,26,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "desc" ,128,1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[4]} ,
  {0, "version" ,128,2,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[5]} ,
  {0, "dev-date" ,128,3,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[6]} ,
  {0, "bld-date" ,128,4,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[7]} ,
  {0, "cit" ,128,5,0,1,0,0,0,0,NULL,&atx[9],&atx[8],0,&atx[10]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {312, "SEQUENCE OF" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "notice" ,128,6,0,1,0,0,0,0,NULL,&atx[9],&atx[11],0,&atx[12]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "prog-usage" ,128,7,0,1,0,0,0,0,NULL,&atx[9],&atx[13],0,&atx[14]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "susage" ,128,8,0,0,0,0,0,0,NULL,&atx[15],NULL,0,&atx[21]} ,
  {420, "BLAST0-Seq-usage" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[16],0,&atx[17]} ,
  {0, "raw" ,128,0,0,0,0,0,0,0,NULL,&atx[17],NULL,0,&atx[19]} ,
  {421, "BLAST0-Alphatype" ,1,0,0,0,0,0,0,0,NULL,&atx[18],&avnx[0],0,&atx[69]} ,
  {310, "ENUMERATED" ,0,10,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "cooked" ,128,1,0,0,0,0,0,0,NULL,&atx[17],NULL,0,NULL} ,
  {311, "SEQUENCE" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "qusage" ,128,9,0,0,0,0,0,0,NULL,&atx[15],NULL,0,NULL} ,
  {402, "BLAST0-Job-desc" ,1,0,0,0,0,1,0,0,NULL,&atx[20],&atx[23],0,&atx[27]} ,
  {0, "jid" ,128,0,0,0,0,0,0,0,NULL,&atx[18],&avnx[4],0,&atx[24]} ,
  {0, "desc" ,128,1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[25]} ,
  {0, "size" ,128,2,0,0,0,0,0,0,NULL,&atx[26],NULL,0,NULL} ,
  {302, "INTEGER" ,0,2,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {403, "BLAST0-Job-progress" ,1,0,0,0,0,1,0,0,NULL,&atx[20],&atx[28],0,&atx[30]} ,
  {0, "done" ,128,0,0,0,0,0,0,0,NULL,&atx[26],NULL,0,&atx[29]} ,
  {0, "positives" ,128,1,0,0,0,0,0,0,NULL,&atx[26],NULL,0,NULL} ,
  {404, "BLAST0-Sequence" ,1,0,0,0,0,1,0,0,NULL,&atx[20],&atx[31],0,&atx[50]} ,
  {0, "desc" ,128,0,0,0,0,0,0,0,NULL,&atx[9],&atx[32],0,&atx[42]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[33],NULL,0,NULL} ,
  {429, "BLAST0-Seq-desc" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[34],0,&atx[35]} ,
  {0, "id" ,128,0,0,0,0,0,0,0,NULL,&atx[35],NULL,0,&atx[41]} ,
  {430, "BLAST0-Seq-id" ,1,0,0,0,0,0,0,0,NULL,&atx[40],&atx[36],0,&atx[111]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[39],&atx[37],0,NULL} ,
  {0, "giid" ,128,0,0,0,0,0,0,0,NULL,&atx[26],NULL,0,&atx[38]} ,
  {0, "textid" ,128,1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {315, "CHOICE" ,0,-1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {314, "SET OF" ,0,17,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "defline" ,128,1,0,1,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "length" ,128,1,0,0,0,0,0,0,NULL,&atx[26],NULL,0,&atx[43]} ,
  {0, "gcode" ,128,2,0,1,0,0,0,0,NULL,&atx[26],NULL,0,&atx[44]} ,
  {0, "seq" ,128,3,0,1,0,0,0,0,NULL,&atx[45],NULL,0,NULL} ,
  {428, "BLAST0-Seq-data" ,1,0,0,0,0,0,0,0,NULL,&atx[39],&atx[46],0,&atx[33]} ,
  {0, "ncbistdaa" ,128,0,0,0,0,0,0,0,NULL,&atx[47],NULL,0,&atx[48]} ,
  {304, "OCTET STRING" ,0,4,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "ncbi2na" ,128,1,0,0,0,0,0,0,NULL,&atx[47],NULL,0,&atx[49]} ,
  {0, "ncbi4na" ,128,2,0,0,0,0,0,0,NULL,&atx[47],NULL,0,NULL} ,
  {405, "BLAST0-KA-Blk" ,1,0,0,0,0,1,0,0,NULL,&atx[20],&atx[51],0,&atx[58]} ,
  {0, "matid" ,128,0,0,0,0,0,0,0,NULL,&atx[26],NULL,0,&atx[52]} ,
  {0, "frames" ,128,1,0,0,0,0,0,0,NULL,&atx[9],&atx[53],0,&atx[54]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[26],NULL,0,NULL} ,
  {0, "lambda" ,128,2,0,0,0,0,0,0,NULL,&atx[55],NULL,0,&atx[56]} ,
  {309, "REAL" ,0,9,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "k" ,128,3,0,0,0,0,0,0,NULL,&atx[55],NULL,0,&atx[57]} ,
  {0, "h" ,128,4,0,0,0,0,0,0,NULL,&atx[55],NULL,0,NULL} ,
  {406, "BLAST0-Db-Desc" ,1,0,0,0,0,1,0,0,NULL,&atx[20],&atx[59],0,&atx[67]} ,
  {0, "name" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[60]} ,
  {0, "type" ,128,1,0,0,0,0,0,0,NULL,&atx[17],NULL,0,&atx[61]} ,
  {0, "def" ,128,2,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[62]} ,
  {0, "rel-date" ,128,3,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[63]} ,
  {0, "bld-date" ,128,4,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[64]} ,
  {0, "count" ,128,5,0,1,0,0,0,0,NULL,&atx[26],NULL,0,&atx[65]} ,
  {0, "totlen" ,128,6,0,1,0,0,0,0,NULL,&atx[26],NULL,0,&atx[66]} ,
  {0, "maxlen" ,128,7,0,1,0,0,0,0,NULL,&atx[26],NULL,0,NULL} ,
  {407, "BLAST0-Result" ,1,0,0,0,0,1,0,0,NULL,&atx[20],&atx[68],0,&atx[105]} ,
  {0, "hist" ,128,0,0,1,0,0,0,0,NULL,&atx[69],NULL,0,&atx[79]} ,
  {422, "BLAST0-Histogram" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[70],0,&atx[82]} ,
  {0, "expect" ,128,0,0,0,0,0,0,0,NULL,&atx[55],NULL,0,&atx[71]} ,
  {0, "observed" ,128,1,0,0,0,0,0,0,NULL,&atx[26],NULL,0,&atx[72]} ,
  {0, "base" ,128,2,0,0,0,0,0,0,NULL,&atx[26],NULL,0,&atx[73]} ,
  {0, "nbars" ,128,3,0,0,0,0,0,0,NULL,&atx[26],NULL,0,&atx[74]} ,
  {0, "bar" ,128,4,0,0,0,0,0,0,NULL,&atx[9],&atx[75],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[76],NULL,0,NULL} ,
  {424, "BLAST0-Histogram-bar" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[77],0,&atx[88]} ,
  {0, "x" ,128,0,0,0,0,0,0,0,NULL,&atx[55],NULL,0,&atx[78]} ,
  {0, "n" ,128,1,0,0,0,0,0,0,NULL,&atx[26],NULL,0,NULL} ,
  {0, "count" ,128,1,0,0,0,0,0,0,NULL,&atx[26],NULL,0,&atx[80]} ,
  {0, "hitlists" ,128,2,0,0,0,0,0,0,NULL,&atx[9],&atx[81],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[82],NULL,0,NULL} ,
  {423, "BLAST0-HitList" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[83],0,&atx[76]} ,
  {0, "count" ,128,0,0,0,0,0,0,0,NULL,&atx[26],NULL,0,&atx[84]} ,
  {0, "kablk" ,128,1,0,1,0,0,0,0,NULL,&atx[9],&atx[85],0,&atx[86]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[50],NULL,0,NULL} ,
  {0, "hsps" ,128,2,0,0,0,0,0,0,NULL,&atx[9],&atx[87],0,&atx[103]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[88],NULL,0,NULL} ,
  {425, "BLAST0-HSP" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[89],0,&atx[95]} ,
  {0, "matid" ,128,0,0,0,0,0,0,0,NULL,&atx[26],NULL,0,&atx[90]} ,
  {0, "scores" ,128,1,0,0,0,0,0,0,NULL,&atx[91],NULL,0,&atx[92]} ,
  {413, "Score-set" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[144]} ,
  {0, "len" ,128,2,0,0,0,0,0,0,NULL,&atx[26],NULL,0,&atx[93]} ,
  {0, "segs" ,128,3,0,0,0,0,0,0,NULL,&atx[9],&atx[94],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[95],NULL,0,NULL} ,
  {426, "BLAST0-Segment" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[96],0,&atx[97]} ,
  {0, "loc" ,128,0,0,0,0,0,0,0,NULL,&atx[97],NULL,0,&atx[101]} ,
  {427, "BLAST0-Seq-interval" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[98],0,&atx[45]} ,
  {0, "strand" ,128,0,0,1,0,0,0,0,NULL,&atx[18],&avnx[8],0,&atx[99]} ,
  {0, "from" ,128,1,0,0,0,0,0,0,NULL,&atx[26],NULL,0,&atx[100]} ,
  {0, "to" ,128,2,0,0,0,0,0,0,NULL,&atx[26],NULL,0,NULL} ,
  {0, "str" ,128,1,0,1,0,0,0,0,NULL,&atx[45],NULL,0,&atx[102]} ,
  {0, "str-raw" ,128,2,0,1,0,0,0,0,NULL,&atx[45],NULL,0,NULL} ,
  {0, "seqs" ,128,3,0,0,0,0,0,0,NULL,&atx[9],&atx[104],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[30],NULL,0,NULL} ,
  {408, "BLAST0-Matrix" ,1,0,0,0,0,1,0,0,NULL,&atx[20],&atx[106],0,&atx[120]} ,
  {0, "matid" ,128,0,0,0,0,0,0,0,NULL,&atx[26],NULL,0,&atx[107]} ,
  {0, "name" ,128,1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[108]} ,
  {0, "comments" ,128,2,0,1,0,0,0,0,NULL,&atx[9],&atx[109],0,&atx[110]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "qalpha" ,128,3,0,0,0,0,0,0,NULL,&atx[111],NULL,0,&atx[112]} ,
  {431, "BLAST0-Alpha-ID" ,1,0,0,0,0,0,0,0,NULL,&atx[18],&avnx[13],0,NULL} ,
  {0, "salpha" ,128,4,0,0,0,0,0,0,NULL,&atx[111],NULL,0,&atx[113]} ,
  {0, "scores" ,128,5,0,1,0,0,0,0,NULL,&atx[39],&atx[114],0,NULL} ,
  {0, "scaled-ints" ,128,0,0,0,0,0,0,0,NULL,&atx[20],&atx[115],0,&atx[118]} ,
  {0, "scale" ,128,0,0,0,0,0,0,0,NULL,&atx[55],NULL,0,&atx[116]} ,
  {0, "ints" ,128,1,0,0,0,0,0,0,NULL,&atx[9],&atx[117],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[26],NULL,0,NULL} ,
  {0, "reals" ,128,1,0,0,0,0,0,0,NULL,&atx[9],&atx[119],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[55],NULL,0,NULL} ,
  {409, "BLAST0-Warning" ,1,0,0,0,0,1,0,0,NULL,&atx[121],NULL,0,&atx[121]} ,
  {410, "BLAST0-Status" ,1,0,0,0,0,1,0,0,NULL,&atx[20],&atx[122],0,&atx[124]} ,
  {0, "code" ,128,0,0,0,0,0,0,0,NULL,&atx[26],NULL,0,&atx[123]} ,
  {0, "reason" ,128,1,0,1,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {411, "BLAST0-Outblk" ,1,0,0,0,0,1,0,0,NULL,&atx[40],&atx[125],0,&atx[143]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[39],&atx[126],0,NULL} ,
  {0, "preface" ,128,0,0,0,0,0,0,0,NULL,&atx[0],NULL,0,&atx[127]} ,
  {0, "query" ,128,1,0,0,0,0,0,0,NULL,&atx[30],NULL,0,&atx[128]} ,
  {0, "dbdesc" ,128,2,0,0,0,0,0,0,NULL,&atx[58],NULL,0,&atx[129]} ,
  {0, "matrix" ,128,3,0,0,0,0,0,0,NULL,&atx[9],&atx[130],0,&atx[131]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[105],NULL,0,NULL} ,
  {0, "kablk" ,128,4,0,0,0,0,0,0,NULL,&atx[9],&atx[132],0,&atx[133]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[50],NULL,0,NULL} ,
  {0, "job-start" ,128,5,0,0,0,0,0,0,NULL,&atx[22],NULL,0,&atx[134]} ,
  {0, "job-progress" ,128,6,0,0,0,0,0,0,NULL,&atx[27],NULL,0,&atx[135]} ,
  {0, "job-done" ,128,7,0,0,0,0,0,0,NULL,&atx[27],NULL,0,&atx[136]} ,
  {0, "result" ,128,8,0,0,0,0,0,0,NULL,&atx[67],NULL,0,&atx[137]} ,
  {0, "parms" ,128,9,0,0,0,0,0,0,NULL,&atx[9],&atx[138],0,&atx[139]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "stats" ,128,10,0,0,0,0,0,0,NULL,&atx[9],&atx[140],0,&atx[141]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "warning" ,128,11,0,0,0,0,0,0,NULL,&atx[120],NULL,0,&atx[142]} ,
  {0, "status" ,128,12,0,0,0,0,0,0,NULL,&atx[121],NULL,0,NULL} ,
  {412, "Seq-align-set" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[91]} ,
  {414, "BLAST0-Request" ,1,0,0,0,0,0,0,0,NULL,&atx[39],&atx[145],0,&atx[154]} ,
  {0, "hello" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[146]} ,
  {0, "motd" ,128,1,0,0,0,0,0,0,NULL,&atx[147],NULL,0,&atx[148]} ,
  {305, "NULL" ,0,5,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "prog-info" ,128,2,0,0,0,0,0,0,NULL,&atx[147],NULL,0,&atx[149]} ,
  {0, "usage-info" ,128,3,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[150]} ,
  {0, "db-info" ,128,4,0,0,0,0,0,0,NULL,&atx[147],NULL,0,&atx[151]} ,
  {0, "matrix-info" ,128,5,0,0,0,0,0,0,NULL,&atx[147],NULL,0,&atx[152]} ,
  {0, "matrix-get" ,128,6,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[153]} ,
  {0, "search" ,128,7,0,0,0,0,0,0,NULL,&atx[154],NULL,0,&atx[166]} ,
  {415, "BLAST0-Search" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[155],0,&atx[167]} ,
  {0, "program" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[156]} ,
  {0, "database" ,128,1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[157]} ,
  {0, "query" ,128,2,0,0,0,0,0,0,NULL,&atx[30],NULL,0,&atx[158]} ,
  {0, "options" ,128,3,0,1,0,0,0,0,NULL,&atx[9],&atx[159],0,&atx[160]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "return-matrix" ,128,4,0,0,1,0,0,0,&avnx[15],&atx[161],NULL,0,&atx[162]} ,
  {301, "BOOLEAN" ,0,1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "return-query" ,128,5,0,0,1,0,0,0,&avnx[16],&atx[161],NULL,0,&atx[163]} ,
  {0, "return-BLAST0result" ,128,6,0,0,1,0,0,0,&avnx[17],&atx[161],NULL,0,&atx[164]} ,
  {0, "return-query-seq-in-seg" ,128,7,0,0,1,0,0,0,&avnx[18],&atx[161],NULL,0,&atx[165]} ,
  {0, "return-db-seq-in-seg" ,128,8,0,0,1,0,0,0,&avnx[19],&atx[161],NULL,0,NULL} ,
  {0, "goodbye" ,128,8,0,0,0,0,0,0,NULL,&atx[147],NULL,0,NULL} ,
  {416, "BLAST0-Response" ,1,0,0,0,0,0,0,0,NULL,&atx[39],&atx[168],0,&atx[178]} ,
  {0, "hello" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[169]} ,
  {0, "motd" ,128,1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[170]} ,
  {0, "prog-info" ,128,2,0,0,0,0,0,0,NULL,&atx[9],&atx[171],0,&atx[172]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[0],NULL,0,NULL} ,
  {0, "usage-info" ,128,3,0,0,0,0,0,0,NULL,&atx[0],NULL,0,&atx[173]} ,
  {0, "db-info" ,128,4,0,0,0,0,0,0,NULL,&atx[9],&atx[174],0,&atx[175]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[58],NULL,0,NULL} ,
  {0, "matrix-info" ,128,5,0,0,0,0,0,0,NULL,&atx[9],&atx[176],0,&atx[177]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[105],NULL,0,NULL} ,
  {0, "ack" ,128,6,0,0,0,0,0,0,NULL,&atx[178],NULL,0,&atx[181]} ,
  {417, "BLAST0-Ack" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[179],0,&atx[183]} ,
  {0, "code" ,128,0,0,0,0,0,0,0,NULL,&atx[26],NULL,0,&atx[180]} ,
  {0, "reason" ,128,1,0,1,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "goodbye" ,128,7,0,0,0,0,0,0,NULL,&atx[178],NULL,0,&atx[182]} ,
  {0, "queued" ,128,8,0,0,0,0,0,0,NULL,&atx[183],NULL,0,&atx[186]} ,
  {418, "BLAST0-Queued" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[184],0,&atx[198]} ,
  {0, "name" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[185]} ,
  {0, "length" ,128,1,0,0,0,0,0,0,NULL,&atx[26],NULL,0,NULL} ,
  {0, "preface" ,128,9,0,0,0,0,0,0,NULL,&atx[0],NULL,0,&atx[187]} ,
  {0, "query" ,128,10,0,0,0,0,0,0,NULL,&atx[30],NULL,0,&atx[188]} ,
  {0, "dbdesc" ,128,11,0,0,0,0,0,0,NULL,&atx[58],NULL,0,&atx[189]} ,
  {0, "matrix" ,128,12,0,0,0,0,0,0,NULL,&atx[9],&atx[190],0,&atx[191]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[105],NULL,0,NULL} ,
  {0, "kablk" ,128,13,0,0,0,0,0,0,NULL,&atx[9],&atx[192],0,&atx[193]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[50],NULL,0,NULL} ,
  {0, "job-start" ,128,14,0,0,0,0,0,0,NULL,&atx[22],NULL,0,&atx[194]} ,
  {0, "job-progress" ,128,15,0,0,0,0,0,0,NULL,&atx[27],NULL,0,&atx[195]} ,
  {0, "job-done" ,128,16,0,0,0,0,0,0,NULL,&atx[27],NULL,0,&atx[196]} ,
  {0, "score-defs" ,128,17,0,0,0,0,0,0,NULL,&atx[9],&atx[197],0,&atx[202]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[198],NULL,0,NULL} ,
  {419, "BLAST0-Score-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[20],&atx[199],0,&atx[15]} ,
  {0, "sid" ,128,0,0,0,0,0,0,0,NULL,&atx[26],NULL,0,&atx[200]} ,
  {0, "tag" ,128,1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[201]} ,
  {0, "desc" ,128,2,0,1,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "result" ,128,18,0,0,0,0,0,0,NULL,&atx[67],NULL,0,&atx[203]} ,
  {0, "seqalign" ,128,19,0,0,0,0,0,0,NULL,&atx[143],NULL,0,&atx[204]} ,
  {0, "parms" ,128,20,0,0,0,0,0,0,NULL,&atx[9],&atx[205],0,&atx[206]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "stats" ,128,21,0,0,0,0,0,0,NULL,&atx[9],&atx[207],0,&atx[208]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "warning" ,128,22,0,0,0,0,0,0,NULL,&atx[120],NULL,0,&atx[209]} ,
  {0, "status" ,128,23,0,0,0,0,0,0,NULL,&atx[121],NULL,0,NULL} };

static AsnModule ampx[1] = {
  { "NCBI-BLAST-1" , "asnbl18.h60",&atx[0],NULL,NULL,0,0} };

static AsnValxNodePtr avn = avnx;
static AsnTypePtr at = atx;
static AsnModulePtr amp = ampx;



/**************************************************
*
*    Defines for Module NCBI-BLAST-1
*
**************************************************/

#define BLAST0_PREFACE &at[0]
#define BLAST0_PREFACE_program &at[1]
#define BLAST0_PREFACE_desc &at[3]
#define BLAST0_PREFACE_version &at[4]
#define BLAST0_PREFACE_dev_date &at[5]
#define BLAST0_PREFACE_bld_date &at[6]
#define BLAST0_PREFACE_cit &at[7]
#define BLAST0_PREFACE_cit_E &at[8]
#define BLAST0_PREFACE_notice &at[10]
#define BLAST0_PREFACE_notice_E &at[11]
#define BLAST0_PREFACE_prog_usage &at[12]
#define BLAST0_PREFACE_prog_usage_E &at[13]
#define BLAST0_PREFACE_susage &at[14]
#define BLAST0_PREFACE_qusage &at[21]

#define BLAST0_JOB_DESC &at[22]
#define BLAST0_JOB_DESC_jid &at[23]
#define BLAST0_JOB_DESC_desc &at[24]
#define BLAST0_JOB_DESC_size &at[25]

#define BLAST0_JOB_PROGRESS &at[27]
#define BLAST0_JOB_PROGRESS_done &at[28]
#define BLAST0_JOB_PROGRESS_positives &at[29]

#define BLAST0_SEQUENCE &at[30]
#define BLAST0_SEQUENCE_desc &at[31]
#define BLAST0_SEQUENCE_desc_E &at[32]
#define BLAST0_SEQUENCE_length &at[42]
#define BLAST0_SEQUENCE_gcode &at[43]
#define BLAST0_SEQUENCE_seq &at[44]

#define BLAST0_KA_BLK &at[50]
#define BLAST0_KA_BLK_matid &at[51]
#define BLAST0_KA_BLK_frames &at[52]
#define BLAST0_KA_BLK_frames_E &at[53]
#define BLAST0_KA_BLK_lambda &at[54]
#define BLAST0_KA_BLK_k &at[56]
#define BLAST0_KA_BLK_h &at[57]

#define BLAST0_DB_DESC &at[58]
#define BLAST0_DB_DESC_name &at[59]
#define BLAST0_DB_DESC_type &at[60]
#define BLAST0_DB_DESC_def &at[61]
#define BLAST0_DB_DESC_rel_date &at[62]
#define BLAST0_DB_DESC_bld_date &at[63]
#define BLAST0_DB_DESC_count &at[64]
#define BLAST0_DB_DESC_totlen &at[65]
#define BLAST0_DB_DESC_maxlen &at[66]

#define BLAST0_RESULT &at[67]
#define BLAST0_RESULT_hist &at[68]
#define BLAST0_RESULT_count &at[79]
#define BLAST0_RESULT_hitlists &at[80]
#define BLAST0_RESULT_hitlists_E &at[81]

#define BLAST0_MATRIX &at[105]
#define BLAST0_MATRIX_matid &at[106]
#define BLAST0_MATRIX_name &at[107]
#define BLAST0_MATRIX_comments &at[108]
#define BLAST0_MATRIX_comments_E &at[109]
#define BLAST0_MATRIX_qalpha &at[110]
#define BLAST0_MATRIX_salpha &at[112]
#define BLAST0_MATRIX_scores &at[113]
#define MATRIX_scores_scaled_ints &at[114]
#define scores_scaled_ints_scale &at[115]
#define MATRIX_scores_scaled_ints_ints &at[116]
#define scores_scaled_ints_ints_E &at[117]
#define BLAST0_MATRIX_scores_reals &at[118]
#define BLAST0_MATRIX_scores_reals_E &at[119]

#define BLAST0_WARNING &at[120]

#define BLAST0_STATUS &at[121]
#define BLAST0_STATUS_code &at[122]
#define BLAST0_STATUS_reason &at[123]

#define BLAST0_OUTBLK &at[124]
#define BLAST0_OUTBLK_E &at[125]
#define BLAST0_OUTBLK_E_preface &at[126]
#define BLAST0_OUTBLK_E_query &at[127]
#define BLAST0_OUTBLK_E_dbdesc &at[128]
#define BLAST0_OUTBLK_E_matrix &at[129]
#define BLAST0_OUTBLK_E_matrix_E &at[130]
#define BLAST0_OUTBLK_E_kablk &at[131]
#define BLAST0_OUTBLK_E_kablk_E &at[132]
#define BLAST0_OUTBLK_E_job_start &at[133]
#define BLAST0_OUTBLK_E_job_progress &at[134]
#define BLAST0_OUTBLK_E_job_done &at[135]
#define BLAST0_OUTBLK_E_result &at[136]
#define BLAST0_OUTBLK_E_parms &at[137]
#define BLAST0_OUTBLK_E_parms_E &at[138]
#define BLAST0_OUTBLK_E_stats &at[139]
#define BLAST0_OUTBLK_E_stats_E &at[140]
#define BLAST0_OUTBLK_E_warning &at[141]
#define BLAST0_OUTBLK_E_status &at[142]

#define BLAST0_REQUEST &at[144]
#define BLAST0_REQUEST_hello &at[145]
#define BLAST0_REQUEST_motd &at[146]
#define BLAST0_REQUEST_prog_info &at[148]
#define BLAST0_REQUEST_usage_info &at[149]
#define BLAST0_REQUEST_db_info &at[150]
#define BLAST0_REQUEST_matrix_info &at[151]
#define BLAST0_REQUEST_matrix_get &at[152]
#define BLAST0_REQUEST_search &at[153]
#define BLAST0_REQUEST_goodbye &at[166]

#define BLAST0_SEARCH &at[154]
#define BLAST0_SEARCH_program &at[155]
#define BLAST0_SEARCH_database &at[156]
#define BLAST0_SEARCH_query &at[157]
#define BLAST0_SEARCH_options &at[158]
#define BLAST0_SEARCH_options_E &at[159]
#define BLAST0_SEARCH_return_matrix &at[160]
#define BLAST0_SEARCH_return_query &at[162]
#define SEARCH_return_BLAST0result &at[163]
#define SEARCH_return_query_seq_in_seg &at[164]
#define SEARCH_return_db_seq_in_seg &at[165]

#define BLAST0_RESPONSE &at[167]
#define BLAST0_RESPONSE_hello &at[168]
#define BLAST0_RESPONSE_motd &at[169]
#define BLAST0_RESPONSE_prog_info &at[170]
#define BLAST0_RESPONSE_prog_info_E &at[171]
#define BLAST0_RESPONSE_usage_info &at[172]
#define BLAST0_RESPONSE_db_info &at[173]
#define BLAST0_RESPONSE_db_info_E &at[174]
#define BLAST0_RESPONSE_matrix_info &at[175]
#define BLAST0_RESPONSE_matrix_info_E &at[176]
#define BLAST0_RESPONSE_ack &at[177]
#define BLAST0_RESPONSE_goodbye &at[181]
#define BLAST0_RESPONSE_queued &at[182]
#define BLAST0_RESPONSE_preface &at[186]
#define BLAST0_RESPONSE_query &at[187]
#define BLAST0_RESPONSE_dbdesc &at[188]
#define BLAST0_RESPONSE_matrix &at[189]
#define BLAST0_RESPONSE_matrix_E &at[190]
#define BLAST0_RESPONSE_kablk &at[191]
#define BLAST0_RESPONSE_kablk_E &at[192]
#define BLAST0_RESPONSE_job_start &at[193]
#define BLAST0_RESPONSE_job_progress &at[194]
#define BLAST0_RESPONSE_job_done &at[195]
#define BLAST0_RESPONSE_score_defs &at[196]
#define BLAST0_RESPONSE_score_defs_E &at[197]
#define BLAST0_RESPONSE_result &at[202]
#define BLAST0_RESPONSE_seqalign &at[203]
#define BLAST0_RESPONSE_parms &at[204]
#define BLAST0_RESPONSE_parms_E &at[205]
#define BLAST0_RESPONSE_stats &at[206]
#define BLAST0_RESPONSE_stats_E &at[207]
#define BLAST0_RESPONSE_warning &at[208]
#define BLAST0_RESPONSE_status &at[209]

#define BLAST0_ACK &at[178]
#define BLAST0_ACK_code &at[179]
#define BLAST0_ACK_reason &at[180]

#define BLAST0_QUEUED &at[183]
#define BLAST0_QUEUED_name &at[184]
#define BLAST0_QUEUED_length &at[185]

#define BLAST0_SCORE_INFO &at[198]
#define BLAST0_SCORE_INFO_sid &at[199]
#define BLAST0_SCORE_INFO_tag &at[200]
#define BLAST0_SCORE_INFO_desc &at[201]

#define BLAST0_SEQ_USAGE &at[15]
#define BLAST0_SEQ_USAGE_raw &at[16]
#define BLAST0_SEQ_USAGE_cooked &at[19]

#define BLAST0_ALPHATYPE &at[17]

#define BLAST0_HISTOGRAM &at[69]
#define BLAST0_HISTOGRAM_expect &at[70]
#define BLAST0_HISTOGRAM_observed &at[71]
#define BLAST0_HISTOGRAM_base &at[72]
#define BLAST0_HISTOGRAM_nbars &at[73]
#define BLAST0_HISTOGRAM_bar &at[74]
#define BLAST0_HISTOGRAM_bar_E &at[75]

#define BLAST0_HITLIST &at[82]
#define BLAST0_HITLIST_count &at[83]
#define BLAST0_HITLIST_kablk &at[84]
#define BLAST0_HITLIST_kablk_E &at[85]
#define BLAST0_HITLIST_hsps &at[86]
#define BLAST0_HITLIST_hsps_E &at[87]
#define BLAST0_HITLIST_seqs &at[103]
#define BLAST0_HITLIST_seqs_E &at[104]

#define BLAST0_HISTOGRAM_BAR &at[76]
#define BLAST0_HISTOGRAM_BAR_x &at[77]
#define BLAST0_HISTOGRAM_BAR_n &at[78]

#define BLAST0_HSP &at[88]
#define BLAST0_HSP_matid &at[89]
#define BLAST0_HSP_scores &at[90]
#define BLAST0_HSP_len &at[92]
#define BLAST0_HSP_segs &at[93]
#define BLAST0_HSP_segs_E &at[94]

#define BLAST0_SEGMENT &at[95]
#define BLAST0_SEGMENT_loc &at[96]
#define BLAST0_SEGMENT_str &at[101]
#define BLAST0_SEGMENT_str_raw &at[102]

#define BLAST0_SEQ_INTERVAL &at[97]
#define BLAST0_SEQ_INTERVAL_strand &at[98]
#define BLAST0_SEQ_INTERVAL_from &at[99]
#define BLAST0_SEQ_INTERVAL_to &at[100]

#define BLAST0_SEQ_DATA &at[45]
#define BLAST0_SEQ_DATA_ncbistdaa &at[46]
#define BLAST0_SEQ_DATA_ncbi2na &at[48]
#define BLAST0_SEQ_DATA_ncbi4na &at[49]

#define BLAST0_SEQ_DESC &at[33]
#define BLAST0_SEQ_DESC_id &at[34]
#define BLAST0_SEQ_DESC_defline &at[41]

#define BLAST0_SEQ_ID &at[35]
#define BLAST0_SEQ_ID_E &at[36]
#define BLAST0_SEQ_ID_E_giid &at[37]
#define BLAST0_SEQ_ID_E_textid &at[38]

#define BLAST0_ALPHA_ID &at[111]
