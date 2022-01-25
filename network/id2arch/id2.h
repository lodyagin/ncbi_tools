/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "id2.h14";
static AsnValxNode avnx[27] = {
    {20,"set-value" ,1,0.0,&avnx[1] } ,
    {20,"get-value" ,2,0.0,&avnx[2] } ,
    {20,"force-value" ,3,0.0,&avnx[3] } ,
    {20,"use-package" ,4,0.0,NULL } ,
    {3,NULL,1,0.0,NULL } ,
    {3,NULL,1,0.0,NULL } ,
    {3,NULL,1,0.0,NULL } ,
    {3,NULL,1,0.0,NULL } ,
    {3,NULL,1,0.0,NULL } ,
    {3,NULL,1,0.0,NULL } ,
    {3,NULL,0,0.0,NULL } ,
    {3,NULL,0,0.0,NULL } ,
    {3,NULL,0,0.0,NULL } ,
    {20,"none" ,0,0.0,&avnx[14] } ,
    {20,"seq-map" ,1,0.0,&avnx[15] } ,
    {20,"seq-data" ,2,0.0,NULL } ,
    {3,NULL,0,0.0,NULL } ,
    {1,"" ,0,0.0,NULL } ,
    {3,NULL,0,0.0,NULL } ,
    {20,"warning" ,1,0.0,&avnx[20] } ,
    {20,"failed-command" ,2,0.0,&avnx[21] } ,
    {20,"failed-connection" ,3,0.0,&avnx[22] } ,
    {20,"failed-server" ,4,0.0,&avnx[23] } ,
    {20,"no-data" ,5,0.0,&avnx[24] } ,
    {20,"restricted-data" ,6,0.0,&avnx[25] } ,
    {20,"unsupported-command" ,7,0.0,&avnx[26] } ,
    {20,"invalid-arguments" ,8,0.0,NULL } };

static AsnType atx[217] = {
  {401, "Seq-id" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[1]} ,
  {402, "Seq-annot" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[2]} ,
  {403, "Seq-descr" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[3]} ,
  {404, "Seq-literal" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[4]} ,
  {405, "Seq-align" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[5]} ,
  {406, "ID2-Request-Packet" ,1,0,0,0,0,0,0,0,NULL,&atx[18],&atx[6],0,&atx[7]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[7],NULL,0,NULL} ,
  {407, "ID2-Request" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[8],0,&atx[11]} ,
  {0, "serial-number" ,128,0,0,1,0,0,0,0,NULL,&atx[9],NULL,0,&atx[10]} ,
  {302, "INTEGER" ,0,2,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "params" ,128,1,0,1,0,0,0,0,NULL,&atx[11],NULL,0,&atx[22]} ,
  {408, "ID2-Params" ,1,0,0,0,0,0,0,0,NULL,&atx[18],&atx[12],0,&atx[26]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[13],NULL,0,NULL} ,
  {448, "ID2-Param" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[14],0,NULL} ,
  {0, "name" ,128,0,0,0,0,0,0,0,NULL,&atx[15],NULL,0,&atx[16]} ,
  {323, "VisibleString" ,0,26,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "value" ,128,1,0,1,0,0,0,0,NULL,&atx[18],&atx[17],0,&atx[19]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[15],NULL,0,NULL} ,
  {312, "SEQUENCE OF" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "type" ,128,2,0,0,1,0,0,0,&avnx[4],&atx[20],&avnx[0],0,NULL} ,
  {310, "ENUMERATED" ,0,10,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {311, "SEQUENCE" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "request" ,128,2,0,0,0,0,0,0,NULL,&atx[42],&atx[23],0,NULL} ,
  {0, "init" ,128,0,0,0,0,0,0,0,NULL,&atx[24],NULL,0,&atx[25]} ,
  {305, "NULL" ,0,5,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "get-packages" ,128,1,0,0,0,0,0,0,NULL,&atx[26],NULL,0,&atx[30]} ,
  {409, "ID2-Request-Get-Packages" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[27],0,&atx[31]} ,
  {0, "names" ,128,0,0,1,0,0,0,0,NULL,&atx[18],&atx[28],0,&atx[29]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[15],NULL,0,NULL} ,
  {0, "no-contents" ,128,1,0,1,0,0,0,0,NULL,&atx[24],NULL,0,NULL} ,
  {0, "string-to-gi" ,128,2,0,0,0,0,0,0,NULL,&atx[31],NULL,0,&atx[33]} ,
  {410, "ID2-Request-String-To-Gi" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[32],0,&atx[34]} ,
  {0, "id" ,128,0,0,0,0,0,0,0,NULL,&atx[15],NULL,0,NULL} ,
  {0, "seq-id-to-gi" ,128,3,0,0,0,0,0,0,NULL,&atx[34],NULL,0,&atx[36]} ,
  {411, "ID2-Request-Seq-id-To-Gi" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[35],0,&atx[37]} ,
  {0, "seq-id" ,128,0,0,0,0,0,0,0,NULL,&atx[0],NULL,0,NULL} ,
  {0, "gi-to-tse-id" ,128,4,0,0,0,0,0,0,NULL,&atx[37],NULL,0,&atx[48]} ,
  {412, "ID2-Request-Gi-To-TSE-Id" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[38],0,&atx[49]} ,
  {0, "gi" ,128,0,0,0,0,0,0,0,NULL,&atx[42],&atx[39],0,&atx[43]} ,
  {0, "gi" ,128,0,0,0,0,0,0,0,NULL,&atx[9],NULL,0,&atx[40]} ,
  {0, "string" ,128,1,0,0,0,0,0,0,NULL,&atx[31],NULL,0,&atx[41]} ,
  {0, "seq-id" ,128,2,0,0,0,0,0,0,NULL,&atx[34],NULL,0,NULL} ,
  {315, "CHOICE" ,0,-1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "sources" ,128,1,0,1,0,0,0,0,NULL,&atx[45],&atx[44],0,&atx[46]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[15],NULL,0,NULL} ,
  {314, "SET OF" ,0,17,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "external" ,128,2,0,1,0,0,0,0,NULL,&atx[24],NULL,0,&atx[47]} ,
  {0, "current-gis" ,128,3,0,1,0,0,0,0,NULL,&atx[24],NULL,0,NULL} ,
  {0, "get-tse" ,128,5,0,0,0,0,0,0,NULL,&atx[49],NULL,0,&atx[89]} ,
  {413, "ID2-Request-Get-TSE" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[50],0,&atx[90]} ,
  {0, "tse-id" ,128,0,0,0,0,0,0,0,NULL,&atx[42],&atx[51],0,&atx[59]} ,
  {0, "tse-id" ,128,0,0,0,0,0,0,0,NULL,&atx[52],NULL,0,&atx[55]} ,
  {416, "ID2-TSE-Id" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[53],0,&atx[60]} ,
  {0, "sat" ,128,0,0,0,0,0,0,0,NULL,&atx[9],NULL,0,&atx[54]} ,
  {0, "sat-key" ,128,1,0,0,0,0,0,0,NULL,&atx[9],NULL,0,NULL} ,
  {0, "gi" ,128,1,0,0,0,0,0,0,NULL,&atx[21],&atx[56],0,NULL} ,
  {0, "request" ,128,0,0,0,0,0,0,0,NULL,&atx[37],NULL,0,&atx[57]} ,
  {0, "exclude-tses" ,128,1,0,1,0,0,0,0,NULL,&atx[45],&atx[58],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[52],NULL,0,NULL} ,
  {0, "details" ,128,1,0,1,0,0,0,0,NULL,&atx[60],NULL,0,NULL} ,
  {417, "ID2-Get-TSE-Details" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[61],0,&atx[99]} ,
  {0, "location" ,128,0,0,0,0,0,0,0,NULL,&atx[62],NULL,0,&atx[83]} ,
  {419, "ID2-Seq-loc" ,1,0,0,0,0,0,0,0,NULL,&atx[42],&atx[63],0,&atx[100]} ,
  {0, "whole" ,128,0,0,0,0,0,0,0,NULL,&atx[9],NULL,0,&atx[64]} ,
  {0, "int" ,128,1,0,0,0,0,0,0,NULL,&atx[65],NULL,0,&atx[69]} ,
  {445, "ID2-Interval" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[66],0,&atx[70]} ,
  {0, "gi" ,128,0,0,0,0,0,0,0,NULL,&atx[9],NULL,0,&atx[67]} ,
  {0, "start" ,128,1,0,0,0,0,0,0,NULL,&atx[9],NULL,0,&atx[68]} ,
  {0, "length" ,128,2,0,0,1,0,0,0,&avnx[5],&atx[9],NULL,0,NULL} ,
  {0, "int-set" ,128,2,0,0,0,0,0,0,NULL,&atx[70],NULL,0,&atx[77]} ,
  {446, "ID2-Packed-Seq-ints" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[71],0,&atx[74]} ,
  {0, "gi" ,128,0,0,0,0,0,0,0,NULL,&atx[9],NULL,0,&atx[72]} ,
  {0, "ints" ,128,1,0,0,0,0,0,0,NULL,&atx[45],&atx[73],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[74],NULL,0,NULL} ,
  {447, "ID2-Seq-range" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[75],0,&atx[13]} ,
  {0, "start" ,128,0,0,0,0,0,0,0,NULL,&atx[9],NULL,0,&atx[76]} ,
  {0, "length" ,128,1,0,0,1,0,0,0,&avnx[6],&atx[9],NULL,0,NULL} ,
  {0, "whole-range" ,128,3,0,0,0,0,0,0,NULL,&atx[78],NULL,0,&atx[81]} ,
  {434, "ID2-Id-Range" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[79],0,&atx[160]} ,
  {0, "start" ,128,0,0,0,0,0,0,0,NULL,&atx[9],NULL,0,&atx[80]} ,
  {0, "count" ,128,1,0,0,1,0,0,0,&avnx[7],&atx[9],NULL,0,NULL} ,
  {0, "loc-set" ,128,4,0,0,0,0,0,0,NULL,&atx[45],&atx[82],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[62],NULL,0,NULL} ,
  {0, "seq-class-level" ,128,1,0,0,1,0,0,0,&avnx[8],&atx[9],NULL,0,&atx[84]} ,
  {0, "descr-level" ,128,2,0,0,1,0,0,0,&avnx[9],&atx[9],NULL,0,&atx[85]} ,
  {0, "descr-type-mask" ,128,3,0,0,1,0,0,0,&avnx[10],&atx[9],NULL,0,&atx[86]} ,
  {0, "annot-type-mask" ,128,4,0,0,1,0,0,0,&avnx[11],&atx[9],NULL,0,&atx[87]} ,
  {0, "feat-type-mask" ,128,5,0,0,1,0,0,0,&avnx[12],&atx[9],NULL,0,&atx[88]} ,
  {0, "sequence-level" ,128,6,0,0,1,0,0,0,&avnx[16],&atx[20],&avnx[13],0,NULL} ,
  {0, "reget-tse" ,128,6,0,0,0,0,0,0,NULL,&atx[90],NULL,0,&atx[94]} ,
  {414, "ID2-Request-ReGet-TSE" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[91],0,&atx[95]} ,
  {0, "tse-id" ,128,0,0,0,0,0,0,0,NULL,&atx[52],NULL,0,&atx[92]} ,
  {0, "details" ,128,1,0,1,0,0,0,0,NULL,&atx[60],NULL,0,&atx[93]} ,
  {0, "offset" ,128,2,0,0,0,0,0,0,NULL,&atx[9],NULL,0,NULL} ,
  {0, "get-chunks" ,128,7,0,0,0,0,0,0,NULL,&atx[95],NULL,0,NULL} ,
  {415, "ID2S-Request-Get-Chunks" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[96],0,&atx[52]} ,
  {0, "tse-id" ,128,0,0,0,0,0,0,0,NULL,&atx[52],NULL,0,&atx[97]} ,
  {0, "chunks" ,128,1,0,0,0,0,0,0,NULL,&atx[45],&atx[98],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[99],NULL,0,NULL} ,
  {418, "ID2S-Chunk-Id" ,1,0,0,0,0,0,0,0,NULL,&atx[9],NULL,0,&atx[62]} ,
  {420, "ID2-Reply" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[101],0,&atx[106]} ,
  {0, "serial-number" ,128,0,0,1,0,0,0,0,NULL,&atx[9],NULL,0,&atx[102]} ,
  {0, "params" ,128,1,0,1,0,0,0,0,NULL,&atx[11],NULL,0,&atx[103]} ,
  {0, "reply" ,128,2,0,0,0,0,0,0,NULL,&atx[42],&atx[104],0,&atx[143]} ,
  {0, "init" ,128,0,0,0,0,0,0,0,NULL,&atx[24],NULL,0,&atx[105]} ,
  {0, "get-package" ,128,1,0,0,0,0,0,0,NULL,&atx[106],NULL,0,&atx[109]} ,
  {421, "ID2-Reply-Get-Package" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[107],0,&atx[110]} ,
  {0, "name" ,128,0,0,0,0,0,0,0,NULL,&atx[15],NULL,0,&atx[108]} ,
  {0, "params" ,128,1,0,1,0,0,0,0,NULL,&atx[11],NULL,0,NULL} ,
  {0, "seq-id-to-gi" ,128,2,0,0,0,0,0,0,NULL,&atx[110],NULL,0,&atx[113]} ,
  {422, "ID2-Reply-Seq-id-To-Gi" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[111],0,&atx[114]} ,
  {0, "seq-id" ,128,0,0,0,0,0,0,0,NULL,&atx[0],NULL,0,&atx[112]} ,
  {0, "gi" ,128,1,0,1,0,0,0,0,NULL,&atx[9],NULL,0,NULL} ,
  {0, "gi-to-tse-id" ,128,3,0,0,0,0,0,0,NULL,&atx[114],NULL,0,&atx[122]} ,
  {423, "ID2-Reply-Gi-To-TSE-Id" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[115],0,&atx[123]} ,
  {0, "gi" ,128,0,0,0,0,0,0,0,NULL,&atx[9],NULL,0,&atx[116]} ,
  {0, "source" ,128,1,0,0,1,0,0,0,&avnx[17],&atx[15],NULL,0,&atx[117]} ,
  {0, "tses" ,128,2,0,1,0,0,0,0,NULL,&atx[45],&atx[118],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[119],NULL,0,NULL} ,
  {428, "ID2-TSE-Id-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[120],0,&atx[126]} ,
  {0, "tse-id" ,128,0,0,0,0,0,0,0,NULL,&atx[52],NULL,0,&atx[121]} ,
  {0, "split-version" ,128,1,0,0,1,0,0,0,&avnx[18],&atx[9],NULL,0,NULL} ,
  {0, "get-tse" ,128,4,0,0,0,0,0,0,NULL,&atx[123],NULL,0,&atx[133]} ,
  {424, "ID2-Reply-Get-TSE" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[124],0,&atx[134]} ,
  {0, "tse-id" ,128,0,0,0,0,0,0,0,NULL,&atx[52],NULL,0,&atx[125]} ,
  {0, "data" ,128,1,0,1,0,0,0,0,NULL,&atx[126],NULL,0,NULL} ,
  {429, "ID2-Reply-Data" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[127],0,&atx[151]} ,
  {0, "data-type" ,128,0,0,0,0,0,0,0,NULL,&atx[9],NULL,0,&atx[128]} ,
  {0, "data-format" ,128,1,0,0,0,0,0,0,NULL,&atx[9],NULL,0,&atx[129]} ,
  {0, "data-compression" ,128,2,0,0,0,0,0,0,NULL,&atx[9],NULL,0,&atx[130]} ,
  {0, "data" ,128,3,0,0,0,0,0,0,NULL,&atx[18],&atx[131],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[132],NULL,0,NULL} ,
  {304, "OCTET STRING" ,0,4,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "get-tse-info" ,128,5,0,0,0,0,0,0,NULL,&atx[134],NULL,0,&atx[138]} ,
  {425, "ID2S-Reply-Get-TSE-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[135],0,&atx[139]} ,
  {0, "tse-id" ,128,0,0,0,0,0,0,0,NULL,&atx[52],NULL,0,&atx[136]} ,
  {0, "split-version" ,128,1,0,0,0,0,0,0,NULL,&atx[9],NULL,0,&atx[137]} ,
  {0, "info" ,128,2,0,1,0,0,0,0,NULL,&atx[126],NULL,0,NULL} ,
  {0, "get-chunk" ,128,6,0,0,0,0,0,0,NULL,&atx[139],NULL,0,NULL} ,
  {426, "ID2S-Reply-Get-Chunk" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[140],0,&atx[145]} ,
  {0, "tse-id" ,128,0,0,0,0,0,0,0,NULL,&atx[52],NULL,0,&atx[141]} ,
  {0, "chunk-id" ,128,1,0,0,0,0,0,0,NULL,&atx[99],NULL,0,&atx[142]} ,
  {0, "data" ,128,2,0,1,0,0,0,0,NULL,&atx[126],NULL,0,NULL} ,
  {0, "error" ,128,3,0,1,0,0,0,0,NULL,&atx[18],&atx[144],0,&atx[149]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[145],NULL,0,NULL} ,
  {427, "ID2-Error" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[146],0,&atx[119]} ,
  {0, "severity" ,128,0,0,0,0,0,0,0,NULL,&atx[20],&avnx[19],0,&atx[147]} ,
  {0, "retry-delay" ,128,1,0,1,0,0,0,0,NULL,&atx[9],NULL,0,&atx[148]} ,
  {0, "message" ,128,2,0,1,0,0,0,0,NULL,&atx[15],NULL,0,NULL} ,
  {0, "end-of-reply" ,128,4,0,0,0,0,0,0,NULL,&atx[150],NULL,0,NULL} ,
  {301, "BOOLEAN" ,0,1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {430, "ID2S-Split-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[152],0,&atx[154]} ,
  {0, "bioseqs-info" ,128,0,0,1,0,0,0,0,NULL,&atx[45],&atx[153],0,&atx[166]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[154],NULL,0,NULL} ,
  {431, "ID2S-Bioseqs-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[155],0,&atx[168]} ,
  {0, "info" ,128,0,0,0,0,0,0,0,NULL,&atx[156],NULL,0,&atx[165]} ,
  {433, "ID2S-Bioseq-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[157],0,&atx[78]} ,
  {0, "gap-count" ,128,0,0,1,0,0,0,0,NULL,&atx[9],NULL,0,&atx[158]} ,
  {0, "seq-map-has-ref" ,128,1,0,1,0,0,0,0,NULL,&atx[150],NULL,0,&atx[159]} ,
  {0, "sequence-split" ,128,2,0,1,0,0,0,0,NULL,&atx[160],NULL,0,NULL} ,
  {435, "ID2S-Sequence-Split-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[161],0,&atx[172]} ,
  {0, "block-size" ,128,0,0,0,0,0,0,0,NULL,&atx[9],NULL,0,&atx[162]} ,
  {0, "chunk-start" ,128,1,0,0,0,0,0,0,NULL,&atx[99],NULL,0,&atx[163]} ,
  {0, "chunk-blocks" ,128,2,0,0,0,0,0,0,NULL,&atx[18],&atx[164],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[9],NULL,0,NULL} ,
  {0, "bioseqs" ,128,1,0,0,0,0,0,0,NULL,&atx[78],NULL,0,NULL} ,
  {0, "chunks" ,128,1,0,0,0,0,0,0,NULL,&atx[45],&atx[167],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[168],NULL,0,NULL} ,
  {432, "ID2S-Chunk-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[169],0,&atx[156]} ,
  {0, "id" ,128,0,0,0,0,0,0,0,NULL,&atx[99],NULL,0,&atx[170]} ,
  {0, "content" ,128,1,0,0,0,0,0,0,NULL,&atx[45],&atx[171],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[172],NULL,0,NULL} ,
  {436, "ID2S-Chunk-Content" ,1,0,0,0,0,0,0,0,NULL,&atx[42],&atx[173],0,&atx[174]} ,
  {0, "seq-descr" ,128,0,0,0,0,0,0,0,NULL,&atx[174],NULL,0,&atx[180]} ,
  {437, "ID2S-Seq-descr-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[175],0,&atx[181]} ,
  {0, "type-mask" ,128,0,0,0,0,0,0,0,NULL,&atx[9],NULL,0,&atx[176]} ,
  {0, "bioseqs" ,128,1,0,1,0,0,0,0,NULL,&atx[45],&atx[177],0,&atx[178]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[78],NULL,0,NULL} ,
  {0, "bioseq-sets" ,128,2,0,1,0,0,0,0,NULL,&atx[45],&atx[179],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[78],NULL,0,NULL} ,
  {0, "seq-annot" ,128,1,0,0,0,0,0,0,NULL,&atx[181],NULL,0,&atx[192]} ,
  {438, "ID2S-Seq-annot-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[182],0,&atx[193]} ,
  {0, "name" ,128,0,0,1,0,0,0,0,NULL,&atx[15],NULL,0,&atx[183]} ,
  {0, "align" ,128,1,0,1,0,0,0,0,NULL,&atx[24],NULL,0,&atx[184]} ,
  {0, "graph" ,128,2,0,1,0,0,0,0,NULL,&atx[24],NULL,0,&atx[185]} ,
  {0, "feat" ,128,3,0,1,0,0,0,0,NULL,&atx[45],&atx[186],0,&atx[191]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[187],NULL,0,NULL} ,
  {442, "ID2S-Feat-type-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[188],0,&atx[200]} ,
  {0, "type" ,128,0,0,0,0,0,0,0,NULL,&atx[9],NULL,0,&atx[189]} ,
  {0, "subtypes" ,128,1,0,1,0,0,0,0,NULL,&atx[45],&atx[190],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[9],NULL,0,NULL} ,
  {0, "seq-loc" ,128,4,0,0,0,0,0,0,NULL,&atx[62],NULL,0,NULL} ,
  {0, "seq-assembly" ,128,2,0,0,0,0,0,0,NULL,&atx[193],NULL,0,&atx[196]} ,
  {439, "ID2S-Seq-assembly-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[194],0,&atx[197]} ,
  {0, "bioseqs" ,128,0,0,0,0,0,0,0,NULL,&atx[45],&atx[195],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[78],NULL,0,NULL} ,
  {0, "seq-map" ,128,3,0,0,0,0,0,0,NULL,&atx[197],NULL,0,&atx[198]} ,
  {440, "ID2S-Seq-map-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[62],NULL,0,&atx[199]} ,
  {0, "seq-data" ,128,4,0,0,0,0,0,0,NULL,&atx[199],NULL,0,NULL} ,
  {441, "ID2S-Seq-data-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[62],NULL,0,&atx[187]} ,
  {443, "ID2S-Chunk" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[201],0,&atx[203]} ,
  {0, "data" ,128,0,0,0,0,0,0,0,NULL,&atx[45],&atx[202],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[203],NULL,0,NULL} ,
  {444, "ID2S-Chunk-Data" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[204],0,&atx[65]} ,
  {0, "id" ,128,0,0,0,0,0,0,0,NULL,&atx[42],&atx[205],0,&atx[207]} ,
  {0, "bioseq-set" ,128,0,0,0,0,0,0,0,NULL,&atx[9],NULL,0,&atx[206]} ,
  {0, "gi" ,128,1,0,0,0,0,0,0,NULL,&atx[9],NULL,0,NULL} ,
  {0, "descrs" ,128,1,0,1,0,0,0,0,NULL,&atx[45],&atx[208],0,&atx[209]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "annots" ,128,2,0,1,0,0,0,0,NULL,&atx[45],&atx[210],0,&atx[211]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[1],NULL,0,NULL} ,
  {0, "assembly" ,128,3,0,1,0,0,0,0,NULL,&atx[45],&atx[212],0,&atx[213]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {0, "seq-map" ,128,4,0,1,0,0,0,0,NULL,&atx[18],&atx[214],0,&atx[215]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[3],NULL,0,NULL} ,
  {0, "seq-data" ,128,5,0,1,0,0,0,0,NULL,&atx[18],&atx[216],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[3],NULL,0,NULL} };

static AsnModule ampx[1] = {
  { "NCBI-ID2Access" , "id2.h14",&atx[0],NULL,NULL,0,0} };

static AsnValxNodePtr avn = avnx;
static AsnTypePtr at = atx;
static AsnModulePtr amp = ampx;



/**************************************************
*
*    Defines for Module NCBI-ID2Access
*
**************************************************/

#define ID2_REQUEST_PACKET &at[5]
#define ID2_REQUEST_PACKET_E &at[6]

#define ID2_REQUEST &at[7]
#define ID2_REQUEST_serial_number &at[8]
#define ID2_REQUEST_params &at[10]
#define ID2_REQUEST_request &at[22]
#define ID2_REQUEST_request_init &at[23]
#define REQUEST_request_get_packages &at[25]
#define REQUEST_request_string_to_gi &at[30]
#define REQUEST_request_seq_id_to_gi &at[33]
#define REQUEST_request_gi_to_tse_id &at[36]
#define ID2_REQUEST_request_get_tse &at[48]
#define ID2_REQUEST_request_reget_tse &at[89]
#define ID2_REQUEST_request_get_chunks &at[94]

#define ID2_PARAMS &at[11]
#define ID2_PARAMS_E &at[12]

#define ID2_REQUEST_GET_PACKAGES &at[26]
#define ID2_REQUEST_GET_PACKAGES_names &at[27]
#define REQUEST_GET_PACKAGES_names_E &at[28]
#define GET_PACKAGES_no_contents &at[29]

#define ID2_REQUEST_STRING_TO_GI &at[31]
#define ID2_REQUEST_STRING_TO_GI_id &at[32]

#define ID2_REQUEST_SEQ_ID_TO_GI &at[34]
#define ID2_REQUEST_SEQ_ID_TO_GI_seq_id &at[35]

#define ID2_REQUEST_GI_TO_TSE_ID &at[37]
#define ID2_REQUEST_GI_TO_TSE_ID_gi &at[38]
#define ID2_REQUEST_GI_TO_TSE_ID_gi_gi &at[39]
#define REQUEST_GI_TO_TSE_ID_gi_string &at[40]
#define REQUEST_GI_TO_TSE_ID_gi_seq_id &at[41]
#define REQUEST_GI_TO_TSE_ID_sources &at[43]
#define REQUEST_GI_TO_TSE_ID_sources_E &at[44]
#define REQUEST_GI_TO_TSE_ID_external &at[46]
#define GI_TO_TSE_ID_current_gis &at[47]

#define ID2_REQUEST_GET_TSE &at[49]
#define ID2_REQUEST_GET_TSE_tse_id &at[50]
#define REQUEST_GET_TSE_tse_id_tse_id &at[51]
#define ID2_REQUEST_GET_TSE_tse_id_gi &at[55]
#define GET_TSE_tse_id_gi_request &at[56]
#define GET_TSE_tse_id_gi_exclude_tses &at[57]
#define TSE_tse_id_gi_exclude_tses_E &at[58]
#define ID2_REQUEST_GET_TSE_details &at[59]

#define ID2_REQUEST_REGET_TSE &at[90]
#define ID2_REQUEST_REGET_TSE_tse_id &at[91]
#define ID2_REQUEST_REGET_TSE_details &at[92]
#define ID2_REQUEST_REGET_TSE_offset &at[93]

#define ID2S_REQUEST_GET_CHUNKS &at[95]
#define ID2S_REQUEST_GET_CHUNKS_tse_id &at[96]
#define ID2S_REQUEST_GET_CHUNKS_chunks &at[97]
#define REQUEST_GET_CHUNKS_chunks_E &at[98]

#define ID2_TSE_ID &at[52]
#define ID2_TSE_ID_sat &at[53]
#define ID2_TSE_ID_sat_key &at[54]

#define ID2_GET_TSE_DETAILS &at[60]
#define ID2_GET_TSE_DETAILS_location &at[61]
#define TSE_DETAILS_seq_class_level &at[83]
#define ID2_GET_TSE_DETAILS_descr_level &at[84]
#define TSE_DETAILS_descr_type_mask &at[85]
#define TSE_DETAILS_annot_type_mask &at[86]
#define GET_TSE_DETAILS_feat_type_mask &at[87]
#define GET_TSE_DETAILS_sequence_level &at[88]

#define ID2S_CHUNK_ID &at[99]

#define ID2_SEQ_LOC &at[62]
#define ID2_SEQ_LOC_whole &at[63]
#define ID2_SEQ_LOC_int &at[64]
#define ID2_SEQ_LOC_int_set &at[69]
#define ID2_SEQ_LOC_whole_range &at[77]
#define ID2_SEQ_LOC_loc_set &at[81]
#define ID2_SEQ_LOC_loc_set_E &at[82]

#define ID2_REPLY &at[100]
#define ID2_REPLY_serial_number &at[101]
#define ID2_REPLY_params &at[102]
#define ID2_REPLY_reply &at[103]
#define ID2_REPLY_reply_init &at[104]
#define ID2_REPLY_reply_get_package &at[105]
#define ID2_REPLY_reply_seq_id_to_gi &at[109]
#define ID2_REPLY_reply_gi_to_tse_id &at[113]
#define ID2_REPLY_reply_get_tse &at[122]
#define ID2_REPLY_reply_get_tse_info &at[133]
#define ID2_REPLY_reply_get_chunk &at[138]
#define ID2_REPLY_error &at[143]
#define ID2_REPLY_error_E &at[144]
#define ID2_REPLY_end_of_reply &at[149]

#define ID2_REPLY_GET_PACKAGE &at[106]
#define ID2_REPLY_GET_PACKAGE_name &at[107]
#define ID2_REPLY_GET_PACKAGE_params &at[108]

#define ID2_REPLY_SEQ_ID_TO_GI &at[110]
#define ID2_REPLY_SEQ_ID_TO_GI_seq_id &at[111]
#define ID2_REPLY_SEQ_ID_TO_GI_gi &at[112]

#define ID2_REPLY_GI_TO_TSE_ID &at[114]
#define ID2_REPLY_GI_TO_TSE_ID_gi &at[115]
#define ID2_REPLY_GI_TO_TSE_ID_source &at[116]
#define ID2_REPLY_GI_TO_TSE_ID_tses &at[117]
#define ID2_REPLY_GI_TO_TSE_ID_tses_E &at[118]

#define ID2_REPLY_GET_TSE &at[123]
#define ID2_REPLY_GET_TSE_tse_id &at[124]
#define ID2_REPLY_GET_TSE_data &at[125]

#define ID2S_REPLY_GET_TSE_INFO &at[134]
#define ID2S_REPLY_GET_TSE_INFO_tse_id &at[135]
#define GET_TSE_INFO_split_version &at[136]
#define ID2S_REPLY_GET_TSE_INFO_info &at[137]

#define ID2S_REPLY_GET_CHUNK &at[139]
#define ID2S_REPLY_GET_CHUNK_tse_id &at[140]
#define ID2S_REPLY_GET_CHUNK_chunk_id &at[141]
#define ID2S_REPLY_GET_CHUNK_data &at[142]

#define ID2_ERROR &at[145]
#define ID2_ERROR_severity &at[146]
#define ID2_ERROR_retry_delay &at[147]
#define ID2_ERROR_message &at[148]

#define ID2_TSE_ID_INFO &at[119]
#define ID2_TSE_ID_INFO_tse_id &at[120]
#define ID2_TSE_ID_INFO_split_version &at[121]

#define ID2_REPLY_DATA &at[126]
#define ID2_REPLY_DATA_data_type &at[127]
#define ID2_REPLY_DATA_data_format &at[128]
#define ID2_REPLY_DATA_data_compression &at[129]
#define ID2_REPLY_DATA_data &at[130]
#define ID2_REPLY_DATA_data_E &at[131]

#define ID2S_SPLIT_INFO &at[151]
#define ID2S_SPLIT_INFO_bioseqs_info &at[152]
#define ID2S_SPLIT_INFO_bioseqs_info_E &at[153]
#define ID2S_SPLIT_INFO_chunks &at[166]
#define ID2S_SPLIT_INFO_chunks_E &at[167]

#define ID2S_BIOSEQS_INFO &at[154]
#define ID2S_BIOSEQS_INFO_info &at[155]
#define ID2S_BIOSEQS_INFO_bioseqs &at[165]

#define ID2S_CHUNK_INFO &at[168]
#define ID2S_CHUNK_INFO_id &at[169]
#define ID2S_CHUNK_INFO_content &at[170]
#define ID2S_CHUNK_INFO_content_E &at[171]

#define ID2S_BIOSEQ_INFO &at[156]
#define ID2S_BIOSEQ_INFO_gap_count &at[157]
#define BIOSEQ_INFO_seq_map_has_ref &at[158]
#define ID2S_BIOSEQ_INFO_sequence_split &at[159]

#define ID2_ID_RANGE &at[78]
#define ID2_ID_RANGE_start &at[79]
#define ID2_ID_RANGE_count &at[80]

#define ID2S_SEQUENCE_SPLIT_INFO &at[160]
#define SEQUENCE_SPLIT_INFO_block_size &at[161]
#define SPLIT_INFO_chunk_start &at[162]
#define SPLIT_INFO_chunk_blocks &at[163]
#define SPLIT_INFO_chunk_blocks_E &at[164]

#define ID2S_CHUNK_CONTENT &at[172]
#define ID2S_CHUNK_CONTENT_seq_descr &at[173]
#define ID2S_CHUNK_CONTENT_seq_annot &at[180]
#define ID2S_CHUNK_CONTENT_seq_assembly &at[192]
#define ID2S_CHUNK_CONTENT_seq_map &at[196]
#define ID2S_CHUNK_CONTENT_seq_data &at[198]

#define ID2S_SEQ_DESCR_INFO &at[174]
#define ID2S_SEQ_DESCR_INFO_type_mask &at[175]
#define ID2S_SEQ_DESCR_INFO_bioseqs &at[176]
#define ID2S_SEQ_DESCR_INFO_bioseqs_E &at[177]
#define ID2S_SEQ_DESCR_INFO_bioseq_sets &at[178]
#define SEQ_DESCR_INFO_bioseq_sets_E &at[179]

#define ID2S_SEQ_ANNOT_INFO &at[181]
#define ID2S_SEQ_ANNOT_INFO_name &at[182]
#define ID2S_SEQ_ANNOT_INFO_align &at[183]
#define ID2S_SEQ_ANNOT_INFO_graph &at[184]
#define ID2S_SEQ_ANNOT_INFO_feat &at[185]
#define ID2S_SEQ_ANNOT_INFO_feat_E &at[186]
#define ID2S_SEQ_ANNOT_INFO_seq_loc &at[191]

#define ID2S_SEQ_ASSEMBLY_INFO &at[193]
#define ID2S_SEQ_ASSEMBLY_INFO_bioseqs &at[194]
#define SEQ_ASSEMBLY_INFO_bioseqs_E &at[195]

#define ID2S_SEQ_MAP_INFO &at[197]

#define ID2S_SEQ_DATA_INFO &at[199]

#define ID2S_FEAT_TYPE_INFO &at[187]
#define ID2S_FEAT_TYPE_INFO_type &at[188]
#define ID2S_FEAT_TYPE_INFO_subtypes &at[189]
#define ID2S_FEAT_TYPE_INFO_subtypes_E &at[190]

#define ID2S_CHUNK &at[200]
#define ID2S_CHUNK_data &at[201]
#define ID2S_CHUNK_data_E &at[202]

#define ID2S_CHUNK_DATA &at[203]
#define ID2S_CHUNK_DATA_id &at[204]
#define ID2S_CHUNK_DATA_id_bioseq_set &at[205]
#define ID2S_CHUNK_DATA_id_gi &at[206]
#define ID2S_CHUNK_DATA_descrs &at[207]
#define ID2S_CHUNK_DATA_descrs_E &at[208]
#define ID2S_CHUNK_DATA_annots &at[209]
#define ID2S_CHUNK_DATA_annots_E &at[210]
#define ID2S_CHUNK_DATA_assembly &at[211]
#define ID2S_CHUNK_DATA_assembly_E &at[212]
#define ID2S_CHUNK_DATA_seq_map &at[213]
#define ID2S_CHUNK_DATA_seq_map_E &at[214]
#define ID2S_CHUNK_DATA_seq_data &at[215]
#define ID2S_CHUNK_DATA_seq_data_E &at[216]

#define ID2_INTERVAL &at[65]
#define ID2_INTERVAL_gi &at[66]
#define ID2_INTERVAL_start &at[67]
#define ID2_INTERVAL_length &at[68]

#define ID2_PACKED_SEQ_INTS &at[70]
#define ID2_PACKED_SEQ_INTS_gi &at[71]
#define ID2_PACKED_SEQ_INTS_ints &at[72]
#define ID2_PACKED_SEQ_INTS_ints_E &at[73]

#define ID2_SEQ_RANGE &at[74]
#define ID2_SEQ_RANGE_start &at[75]
#define ID2_SEQ_RANGE_length &at[76]

#define ID2_PARAM &at[13]
#define ID2_PARAM_name &at[14]
#define ID2_PARAM_value &at[16]
#define ID2_PARAM_value_E &at[17]
#define ID2_PARAM_type &at[19]
