/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "id2.h26";
static AsnValxNode avnx[66] = {
    {20,"live" ,0,0.0,&avnx[1] } ,
    {20,"suppressed-temp" ,1,0.0,&avnx[2] } ,
    {20,"suppressed" ,2,0.0,&avnx[3] } ,
    {20,"dead" ,3,0.0,&avnx[4] } ,
    {20,"protected" ,4,0.0,&avnx[5] } ,
    {20,"withdrawn" ,5,0.0,NULL } ,
    {20,"set-value" ,1,0.0,&avnx[7] } ,
    {20,"get-value" ,2,0.0,&avnx[8] } ,
    {20,"force-value" ,3,0.0,&avnx[9] } ,
    {20,"use-package" ,4,0.0,NULL } ,
    {3,NULL,1,0.0,NULL } ,
    {20,"any" ,0,0.0,&avnx[12] } ,
    {20,"gi" ,1,0.0,&avnx[13] } ,
    {20,"text" ,2,0.0,&avnx[14] } ,
    {20,"general" ,4,0.0,&avnx[15] } ,
    {20,"all" ,127,0.0,&avnx[16] } ,
    {20,"label" ,128,0.0,&avnx[17] } ,
    {20,"taxid" ,256,0.0,&avnx[18] } ,
    {20,"hash" ,512,0.0,NULL } ,
    {3,NULL,0,0.0,NULL } ,
    {20,"main" ,0,0.0,&avnx[21] } ,
    {20,"snp" ,1,0.0,&avnx[22] } ,
    {20,"snp-graph" ,4,0.0,&avnx[23] } ,
    {20,"cdd" ,8,0.0,&avnx[24] } ,
    {20,"mgc" ,16,0.0,&avnx[25] } ,
    {20,"hprd" ,32,0.0,&avnx[26] } ,
    {20,"sts" ,64,0.0,&avnx[27] } ,
    {20,"trna" ,128,0.0,&avnx[28] } ,
    {20,"exon" ,512,0.0,NULL } ,
    {3,NULL,0,0.0,NULL } ,
    {3,NULL,1,0.0,NULL } ,
    {3,NULL,1,0.0,NULL } ,
    {3,NULL,0,0.0,NULL } ,
    {3,NULL,0,0.0,NULL } ,
    {3,NULL,0,0.0,NULL } ,
    {20,"none" ,0,0.0,&avnx[36] } ,
    {20,"seq-map" ,1,0.0,&avnx[37] } ,
    {20,"seq-data" ,2,0.0,NULL } ,
    {3,NULL,0,0.0,NULL } ,
    {20,"warning" ,1,0.0,&avnx[40] } ,
    {20,"failed-command" ,2,0.0,&avnx[41] } ,
    {20,"failed-connection" ,3,0.0,&avnx[42] } ,
    {20,"failed-server" ,4,0.0,&avnx[43] } ,
    {20,"no-data" ,5,0.0,&avnx[44] } ,
    {20,"restricted-data" ,6,0.0,&avnx[45] } ,
    {20,"unsupported-command" ,7,0.0,&avnx[46] } ,
    {20,"invalid-arguments" ,8,0.0,NULL } ,
    {3,NULL,0,0.0,NULL } ,
    {20,"seq-entry" ,0,0.0,&avnx[49] } ,
    {20,"seq-annot" ,1,0.0,&avnx[50] } ,
    {20,"id2s-split-info" ,2,0.0,&avnx[51] } ,
    {20,"id2s-chunk" ,3,0.0,NULL } ,
    {3,NULL,0,0.0,NULL } ,
    {20,"asn-binary" ,0,0.0,&avnx[54] } ,
    {20,"asn-text" ,1,0.0,&avnx[55] } ,
    {20,"xml" ,2,0.0,NULL } ,
    {3,NULL,0,0.0,NULL } ,
    {20,"none" ,0,0.0,&avnx[58] } ,
    {20,"gzip" ,1,0.0,&avnx[59] } ,
    {20,"nlmzip" ,2,0.0,&avnx[60] } ,
    {20,"bzip2" ,3,0.0,NULL } ,
    {3,NULL,0,0.0,NULL } ,
    {3,NULL,0,0.0,NULL } ,
    {20,"reply" ,0,0.0,&avnx[64] } ,
    {20,"last-octet-string" ,1,0.0,&avnx[65] } ,
    {20,"nothing" ,2,0.0,NULL } };

static AsnType atx[150] = {
  {401, "ID2-Blob-State" ,1,0,0,0,0,1,0,0,NULL,&atx[1],&avnx[0],0,&atx[2]} ,
  {310, "ENUMERATED" ,0,10,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {402, "Seq-id" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[3]} ,
  {403, "Seq-loc" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[4]} ,
  {404, "ID2S-Chunk-Id" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[5]} ,
  {405, "ID2S-Seq-annot-Info" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[6]} ,
  {406, "ID2-Request-Packet" ,1,0,0,0,0,0,0,0,NULL,&atx[19],&atx[7],0,&atx[8]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[8],NULL,0,NULL} ,
  {407, "ID2-Request" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[9],0,&atx[12]} ,
  {0, "serial-number" ,128,0,0,1,0,0,0,0,NULL,&atx[10],NULL,0,&atx[11]} ,
  {302, "INTEGER" ,0,2,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "params" ,128,1,0,1,0,0,0,0,NULL,&atx[12],NULL,0,&atx[22]} ,
  {408, "ID2-Params" ,1,0,0,0,0,0,0,0,NULL,&atx[19],&atx[13],0,&atx[26]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[14],NULL,0,NULL} ,
  {431, "ID2-Param" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[15],0,NULL} ,
  {0, "name" ,128,0,0,0,0,0,0,0,NULL,&atx[16],NULL,0,&atx[17]} ,
  {323, "VisibleString" ,0,26,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "value" ,128,1,0,1,0,0,0,0,NULL,&atx[19],&atx[18],0,&atx[20]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[16],NULL,0,NULL} ,
  {312, "SEQUENCE OF" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "type" ,128,2,0,0,1,0,0,0,&avnx[10],&atx[1],&avnx[6],0,NULL} ,
  {311, "SEQUENCE" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "request" ,128,2,0,0,0,0,0,0,NULL,&atx[36],&atx[23],0,NULL} ,
  {0, "init" ,128,0,0,0,0,0,0,0,NULL,&atx[24],NULL,0,&atx[25]} ,
  {305, "NULL" ,0,5,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "get-packages" ,128,1,0,0,0,0,0,0,NULL,&atx[26],NULL,0,&atx[30]} ,
  {409, "ID2-Request-Get-Packages" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[27],0,&atx[31]} ,
  {0, "names" ,128,0,0,1,0,0,0,0,NULL,&atx[19],&atx[28],0,&atx[29]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[16],NULL,0,NULL} ,
  {0, "no-contents" ,128,1,0,1,0,0,0,0,NULL,&atx[24],NULL,0,NULL} ,
  {0, "get-seq-id" ,128,2,0,0,0,0,0,0,NULL,&atx[31],NULL,0,&atx[38]} ,
  {410, "ID2-Request-Get-Seq-id" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[32],0,&atx[39]} ,
  {0, "seq-id" ,128,0,0,0,0,0,0,0,NULL,&atx[33],NULL,0,&atx[37]} ,
  {415, "ID2-Seq-id" ,1,0,0,0,0,0,0,0,NULL,&atx[36],&atx[34],0,&atx[48]} ,
  {0, "string" ,128,0,0,0,0,0,0,0,NULL,&atx[16],NULL,0,&atx[35]} ,
  {0, "seq-id" ,128,1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {315, "CHOICE" ,0,-1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "seq-id-type" ,128,1,0,0,1,0,0,0,&avnx[19],&atx[10],&avnx[11],0,NULL} ,
  {0, "get-blob-id" ,128,3,0,0,0,0,0,0,NULL,&atx[39],NULL,0,&atx[44]} ,
  {411, "ID2-Request-Get-Blob-Id" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[40],0,&atx[45]} ,
  {0, "seq-id" ,128,0,0,0,0,0,0,0,NULL,&atx[31],NULL,0,&atx[41]} ,
  {0, "sources" ,128,1,0,1,0,0,0,0,NULL,&atx[19],&atx[42],0,&atx[43]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[16],NULL,0,NULL} ,
  {0, "external" ,128,2,0,1,0,0,0,0,NULL,&atx[24],NULL,0,NULL} ,
  {0, "get-blob-info" ,128,4,0,0,0,0,0,0,NULL,&atx[45],NULL,0,&atx[67]} ,
  {412, "ID2-Request-Get-Blob-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[46],0,&atx[68]} ,
  {0, "blob-id" ,128,0,0,0,0,0,0,0,NULL,&atx[36],&atx[47],0,&atx[57]} ,
  {0, "blob-id" ,128,0,0,0,0,0,0,0,NULL,&atx[48],NULL,0,&atx[53]} ,
  {416, "ID2-Blob-Id" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[49],0,&atx[59]} ,
  {0, "sat" ,128,0,0,0,0,0,0,0,NULL,&atx[10],NULL,0,&atx[50]} ,
  {0, "sub-sat" ,128,1,0,0,1,0,0,0,&avnx[29],&atx[10],&avnx[20],0,&atx[51]} ,
  {0, "sat-key" ,128,2,0,0,0,0,0,0,NULL,&atx[10],NULL,0,&atx[52]} ,
  {0, "version" ,128,3,0,1,0,0,0,0,NULL,&atx[10],NULL,0,NULL} ,
  {0, "resolve" ,128,1,0,0,0,0,0,0,NULL,&atx[21],&atx[54],0,NULL} ,
  {0, "request" ,128,0,0,0,0,0,0,0,NULL,&atx[39],NULL,0,&atx[55]} ,
  {0, "exclude-blobs" ,128,1,0,1,0,0,0,0,NULL,&atx[19],&atx[56],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[48],NULL,0,NULL} ,
  {0, "get-seq-ids" ,128,1,0,1,0,0,0,0,NULL,&atx[24],NULL,0,&atx[58]} ,
  {0, "get-data" ,128,2,0,1,0,0,0,0,NULL,&atx[59],NULL,0,NULL} ,
  {417, "ID2-Get-Blob-Details" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[60],0,&atx[78]} ,
  {0, "location" ,128,0,0,1,0,0,0,0,NULL,&atx[3],NULL,0,&atx[61]} ,
  {0, "seq-class-level" ,128,1,0,0,1,0,0,0,&avnx[30],&atx[10],NULL,0,&atx[62]} ,
  {0, "descr-level" ,128,2,0,0,1,0,0,0,&avnx[31],&atx[10],NULL,0,&atx[63]} ,
  {0, "descr-type-mask" ,128,3,0,0,1,0,0,0,&avnx[32],&atx[10],NULL,0,&atx[64]} ,
  {0, "annot-type-mask" ,128,4,0,0,1,0,0,0,&avnx[33],&atx[10],NULL,0,&atx[65]} ,
  {0, "feat-type-mask" ,128,5,0,0,1,0,0,0,&avnx[34],&atx[10],NULL,0,&atx[66]} ,
  {0, "sequence-level" ,128,6,0,0,1,0,0,0,&avnx[38],&atx[1],&avnx[35],0,NULL} ,
  {0, "reget-blob" ,128,5,0,0,0,0,0,0,NULL,&atx[68],NULL,0,&atx[72]} ,
  {413, "ID2-Request-ReGet-Blob" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[69],0,&atx[73]} ,
  {0, "blob-id" ,128,0,0,0,0,0,0,0,NULL,&atx[48],NULL,0,&atx[70]} ,
  {0, "split-version" ,128,1,0,0,0,0,0,0,NULL,&atx[10],NULL,0,&atx[71]} ,
  {0, "offset" ,128,2,0,0,0,0,0,0,NULL,&atx[10],NULL,0,NULL} ,
  {0, "get-chunks" ,128,6,0,0,0,0,0,0,NULL,&atx[73],NULL,0,NULL} ,
  {414, "ID2S-Request-Get-Chunks" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[74],0,&atx[33]} ,
  {0, "blob-id" ,128,0,0,0,0,0,0,0,NULL,&atx[48],NULL,0,&atx[75]} ,
  {0, "chunks" ,128,1,0,0,0,0,0,0,NULL,&atx[19],&atx[76],0,&atx[77]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {0, "split-version" ,128,2,0,1,0,0,0,0,NULL,&atx[10],NULL,0,NULL} ,
  {418, "ID2-Reply" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[79],0,&atx[83]} ,
  {0, "serial-number" ,128,0,0,1,0,0,0,0,NULL,&atx[10],NULL,0,&atx[80]} ,
  {0, "params" ,128,1,0,1,0,0,0,0,NULL,&atx[12],NULL,0,&atx[81]} ,
  {0, "error" ,128,2,0,1,0,0,0,0,NULL,&atx[19],&atx[82],0,&atx[87]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[83],NULL,0,NULL} ,
  {419, "ID2-Error" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[84],0,&atx[92]} ,
  {0, "severity" ,128,0,0,0,0,0,0,0,NULL,&atx[1],&avnx[39],0,&atx[85]} ,
  {0, "retry-delay" ,128,1,0,1,0,0,0,0,NULL,&atx[10],NULL,0,&atx[86]} ,
  {0, "message" ,128,2,0,1,0,0,0,0,NULL,&atx[16],NULL,0,NULL} ,
  {0, "end-of-reply" ,128,3,0,1,0,0,0,0,NULL,&atx[24],NULL,0,&atx[88]} ,
  {0, "reply" ,128,4,0,0,0,0,0,0,NULL,&atx[36],&atx[89],0,&atx[144]} ,
  {0, "init" ,128,0,0,0,0,0,0,0,NULL,&atx[24],NULL,0,&atx[90]} ,
  {0, "empty" ,128,1,0,0,0,0,0,0,NULL,&atx[24],NULL,0,&atx[91]} ,
  {0, "get-package" ,128,2,0,0,0,0,0,0,NULL,&atx[92],NULL,0,&atx[95]} ,
  {420, "ID2-Reply-Get-Package" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[93],0,&atx[96]} ,
  {0, "name" ,128,0,0,0,0,0,0,0,NULL,&atx[16],NULL,0,&atx[94]} ,
  {0, "params" ,128,1,0,1,0,0,0,0,NULL,&atx[12],NULL,0,NULL} ,
  {0, "get-seq-id" ,128,3,0,0,0,0,0,0,NULL,&atx[96],NULL,0,&atx[101]} ,
  {421, "ID2-Reply-Get-Seq-id" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[97],0,&atx[102]} ,
  {0, "request" ,128,0,0,0,0,0,0,0,NULL,&atx[31],NULL,0,&atx[98]} ,
  {0, "seq-id" ,128,1,0,1,0,0,0,0,NULL,&atx[19],&atx[99],0,&atx[100]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "end-of-reply" ,128,2,0,1,0,0,0,0,NULL,&atx[24],NULL,0,NULL} ,
  {0, "get-blob-id" ,128,4,0,0,0,0,0,0,NULL,&atx[102],NULL,0,&atx[110]} ,
  {422, "ID2-Reply-Get-Blob-Id" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[103],0,&atx[111]} ,
  {0, "seq-id" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[104]} ,
  {0, "blob-id" ,128,1,0,1,0,0,0,0,NULL,&atx[48],NULL,0,&atx[105]} ,
  {0, "split-version" ,128,2,0,0,1,0,0,0,&avnx[47],&atx[10],NULL,0,&atx[106]} ,
  {0, "annot-info" ,128,3,0,1,0,0,0,0,NULL,&atx[19],&atx[107],0,&atx[108]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[5],NULL,0,NULL} ,
  {0, "end-of-reply" ,128,4,0,1,0,0,0,0,NULL,&atx[24],NULL,0,&atx[109]} ,
  {0, "blob-state" ,128,5,0,1,0,0,0,0,NULL,&atx[10],NULL,0,NULL} ,
  {0, "get-blob-seq-ids" ,128,5,0,0,0,0,0,0,NULL,&atx[111],NULL,0,&atx[121]} ,
  {423, "ID2-Reply-Get-Blob-Seq-ids" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[112],0,&atx[122]} ,
  {0, "blob-id" ,128,0,0,0,0,0,0,0,NULL,&atx[48],NULL,0,&atx[113]} ,
  {0, "ids" ,128,1,0,1,0,0,0,0,NULL,&atx[114],NULL,0,NULL} ,
  {428, "ID2-Reply-Data" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[115],0,&atx[145]} ,
  {0, "data-type" ,128,0,0,0,1,0,0,0,&avnx[52],&atx[10],&avnx[48],0,&atx[116]} ,
  {0, "data-format" ,128,1,0,0,1,0,0,0,&avnx[56],&atx[10],&avnx[53],0,&atx[117]} ,
  {0, "data-compression" ,128,2,0,0,1,0,0,0,&avnx[61],&atx[10],&avnx[57],0,&atx[118]} ,
  {0, "data" ,128,3,0,0,0,0,0,0,NULL,&atx[19],&atx[119],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[120],NULL,0,NULL} ,
  {304, "OCTET STRING" ,0,4,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "get-blob" ,128,6,0,0,0,0,0,0,NULL,&atx[122],NULL,0,&atx[127]} ,
  {424, "ID2-Reply-Get-Blob" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[123],0,&atx[128]} ,
  {0, "blob-id" ,128,0,0,0,0,0,0,0,NULL,&atx[48],NULL,0,&atx[124]} ,
  {0, "split-version" ,128,1,0,0,1,0,0,0,&avnx[62],&atx[10],NULL,0,&atx[125]} ,
  {0, "data" ,128,2,0,1,0,0,0,0,NULL,&atx[114],NULL,0,&atx[126]} ,
  {0, "blob-state" ,128,3,0,1,0,0,0,0,NULL,&atx[10],NULL,0,NULL} ,
  {0, "reget-blob" ,128,7,0,0,0,0,0,0,NULL,&atx[128],NULL,0,&atx[133]} ,
  {425, "ID2-Reply-ReGet-Blob" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[129],0,&atx[134]} ,
  {0, "blob-id" ,128,0,0,0,0,0,0,0,NULL,&atx[48],NULL,0,&atx[130]} ,
  {0, "split-version" ,128,1,0,0,0,0,0,0,NULL,&atx[10],NULL,0,&atx[131]} ,
  {0, "offset" ,128,2,0,0,0,0,0,0,NULL,&atx[10],NULL,0,&atx[132]} ,
  {0, "data" ,128,3,0,1,0,0,0,0,NULL,&atx[114],NULL,0,NULL} ,
  {0, "get-split-info" ,128,8,0,0,0,0,0,0,NULL,&atx[134],NULL,0,&atx[139]} ,
  {426, "ID2S-Reply-Get-Split-Info" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[135],0,&atx[140]} ,
  {0, "blob-id" ,128,0,0,0,0,0,0,0,NULL,&atx[48],NULL,0,&atx[136]} ,
  {0, "split-version" ,128,1,0,0,0,0,0,0,NULL,&atx[10],NULL,0,&atx[137]} ,
  {0, "data" ,128,2,0,1,0,0,0,0,NULL,&atx[114],NULL,0,&atx[138]} ,
  {0, "blob-state" ,128,3,0,1,0,0,0,0,NULL,&atx[10],NULL,0,NULL} ,
  {0, "get-chunk" ,128,9,0,0,0,0,0,0,NULL,&atx[140],NULL,0,NULL} ,
  {427, "ID2S-Reply-Get-Chunk" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[141],0,&atx[114]} ,
  {0, "blob-id" ,128,0,0,0,0,0,0,0,NULL,&atx[48],NULL,0,&atx[142]} ,
  {0, "chunk-id" ,128,1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,&atx[143]} ,
  {0, "data" ,128,2,0,1,0,0,0,0,NULL,&atx[114],NULL,0,NULL} ,
  {0, "discard" ,128,5,0,1,0,0,0,0,NULL,&atx[1],&avnx[63],0,NULL} ,
  {429, "ID2-Blob-Seq-ids" ,1,0,0,0,0,0,0,0,NULL,&atx[19],&atx[146],0,&atx[147]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[147],NULL,0,NULL} ,
  {430, "ID2-Blob-Seq-id" ,1,0,0,0,0,0,0,0,NULL,&atx[21],&atx[148],0,&atx[14]} ,
  {0, "seq-id" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[149]} ,
  {0, "replaced" ,128,1,0,1,0,0,0,0,NULL,&atx[24],NULL,0,NULL} };

static AsnModule ampx[1] = {
  { "NCBI-ID2Access" , "id2.h26",&atx[0],NULL,NULL,0,0} };

static AsnValxNodePtr avn = avnx;
static AsnTypePtr at = atx;
static AsnModulePtr amp = ampx;



/**************************************************
*
*    Defines for Module NCBI-ID2Access
*
**************************************************/

#define ID2_BLOB_STATE &at[0]

#define ID2_REQUEST_PACKET &at[6]
#define ID2_REQUEST_PACKET_E &at[7]

#define ID2_REQUEST &at[8]
#define ID2_REQUEST_serial_number &at[9]
#define ID2_REQUEST_params &at[11]
#define ID2_REQUEST_request &at[22]
#define ID2_REQUEST_request_init &at[23]
#define REQUEST_request_get_packages &at[25]
#define ID2_REQUEST_request_get_seq_id &at[30]
#define ID2_REQUEST_request_get_blob_id &at[38]
#define REQUEST_request_get_blob_info &at[44]
#define ID2_REQUEST_request_reget_blob &at[67]
#define ID2_REQUEST_request_get_chunks &at[72]

#define ID2_PARAMS &at[12]
#define ID2_PARAMS_E &at[13]

#define ID2_REQUEST_GET_PACKAGES &at[26]
#define ID2_REQUEST_GET_PACKAGES_names &at[27]
#define REQUEST_GET_PACKAGES_names_E &at[28]
#define GET_PACKAGES_no_contents &at[29]

#define ID2_REQUEST_GET_SEQ_ID &at[31]
#define ID2_REQUEST_GET_SEQ_ID_seq_id &at[32]
#define REQUEST_GET_SEQ_ID_seq_id_type &at[37]

#define ID2_REQUEST_GET_BLOB_ID &at[39]
#define ID2_REQUEST_GET_BLOB_ID_seq_id &at[40]
#define ID2_REQUEST_GET_BLOB_ID_sources &at[41]
#define REQUEST_GET_BLOB_ID_sources_E &at[42]
#define REQUEST_GET_BLOB_ID_external &at[43]

#define ID2_REQUEST_GET_BLOB_INFO &at[45]
#define REQUEST_GET_BLOB_INFO_blob_id &at[46]
#define GET_BLOB_INFO_blob_id_blob_id &at[47]
#define GET_BLOB_INFO_blob_id_resolve &at[53]
#define INFO_blob_id_resolve_request &at[54]
#define blob_id_resolve_exclude_blobs &at[55]
#define id_resolve_exclude_blobs_E &at[56]
#define GET_BLOB_INFO_get_seq_ids &at[57]
#define REQUEST_GET_BLOB_INFO_get_data &at[58]

#define ID2_REQUEST_REGET_BLOB &at[68]
#define ID2_REQUEST_REGET_BLOB_blob_id &at[69]
#define REGET_BLOB_split_version &at[70]
#define ID2_REQUEST_REGET_BLOB_offset &at[71]

#define ID2S_REQUEST_GET_CHUNKS &at[73]
#define ID2S_REQUEST_GET_CHUNKS_blob_id &at[74]
#define ID2S_REQUEST_GET_CHUNKS_chunks &at[75]
#define REQUEST_GET_CHUNKS_chunks_E &at[76]
#define GET_CHUNKS_split_version &at[77]

#define ID2_SEQ_ID &at[33]
#define ID2_SEQ_ID_string &at[34]
#define ID2_SEQ_ID_seq_id &at[35]

#define ID2_BLOB_ID &at[48]
#define ID2_BLOB_ID_sat &at[49]
#define ID2_BLOB_ID_sub_sat &at[50]
#define ID2_BLOB_ID_sat_key &at[51]
#define ID2_BLOB_ID_version &at[52]

#define ID2_GET_BLOB_DETAILS &at[59]
#define ID2_GET_BLOB_DETAILS_location &at[60]
#define BLOB_DETAILS_seq_class_level &at[61]
#define GET_BLOB_DETAILS_descr_level &at[62]
#define BLOB_DETAILS_descr_type_mask &at[63]
#define BLOB_DETAILS_annot_type_mask &at[64]
#define BLOB_DETAILS_feat_type_mask &at[65]
#define BLOB_DETAILS_sequence_level &at[66]

#define ID2_REPLY &at[78]
#define ID2_REPLY_serial_number &at[79]
#define ID2_REPLY_params &at[80]
#define ID2_REPLY_error &at[81]
#define ID2_REPLY_error_E &at[82]
#define ID2_REPLY_end_of_reply &at[87]
#define ID2_REPLY_reply &at[88]
#define ID2_REPLY_reply_init &at[89]
#define ID2_REPLY_reply_empty &at[90]
#define ID2_REPLY_reply_get_package &at[91]
#define ID2_REPLY_reply_get_seq_id &at[95]
#define ID2_REPLY_reply_get_blob_id &at[101]
#define REPLY_reply_get_blob_seq_ids &at[110]
#define ID2_REPLY_reply_get_blob &at[121]
#define ID2_REPLY_reply_reget_blob &at[127]
#define ID2_REPLY_reply_get_split_info &at[133]
#define ID2_REPLY_reply_get_chunk &at[139]
#define ID2_REPLY_discard &at[144]

#define ID2_ERROR &at[83]
#define ID2_ERROR_severity &at[84]
#define ID2_ERROR_retry_delay &at[85]
#define ID2_ERROR_message &at[86]

#define ID2_REPLY_GET_PACKAGE &at[92]
#define ID2_REPLY_GET_PACKAGE_name &at[93]
#define ID2_REPLY_GET_PACKAGE_params &at[94]

#define ID2_REPLY_GET_SEQ_ID &at[96]
#define ID2_REPLY_GET_SEQ_ID_request &at[97]
#define ID2_REPLY_GET_SEQ_ID_seq_id &at[98]
#define ID2_REPLY_GET_SEQ_ID_seq_id_E &at[99]
#define REPLY_GET_SEQ_ID_end_of_reply &at[100]

#define ID2_REPLY_GET_BLOB_ID &at[102]
#define ID2_REPLY_GET_BLOB_ID_seq_id &at[103]
#define ID2_REPLY_GET_BLOB_ID_blob_id &at[104]
#define GET_BLOB_ID_split_version &at[105]
#define REPLY_GET_BLOB_ID_annot_info &at[106]
#define REPLY_GET_BLOB_ID_annot_info_E &at[107]
#define REPLY_GET_BLOB_ID_end_of_reply &at[108]
#define REPLY_GET_BLOB_ID_blob_state &at[109]

#define ID2_REPLY_GET_BLOB_SEQ_IDS &at[111]
#define REPLY_GET_BLOB_SEQ_IDS_blob_id &at[112]
#define ID2_REPLY_GET_BLOB_SEQ_IDS_ids &at[113]

#define ID2_REPLY_GET_BLOB &at[122]
#define ID2_REPLY_GET_BLOB_blob_id &at[123]
#define REPLY_GET_BLOB_split_version &at[124]
#define ID2_REPLY_GET_BLOB_data &at[125]
#define ID2_REPLY_GET_BLOB_blob_state &at[126]

#define ID2_REPLY_REGET_BLOB &at[128]
#define ID2_REPLY_REGET_BLOB_blob_id &at[129]
#define REPLY_REGET_BLOB_split_version &at[130]
#define ID2_REPLY_REGET_BLOB_offset &at[131]
#define ID2_REPLY_REGET_BLOB_data &at[132]

#define ID2S_REPLY_GET_SPLIT_INFO &at[134]
#define REPLY_GET_SPLIT_INFO_blob_id &at[135]
#define GET_SPLIT_INFO_split_version &at[136]
#define ID2S_REPLY_GET_SPLIT_INFO_data &at[137]
#define GET_SPLIT_INFO_blob_state &at[138]

#define ID2S_REPLY_GET_CHUNK &at[140]
#define ID2S_REPLY_GET_CHUNK_blob_id &at[141]
#define ID2S_REPLY_GET_CHUNK_chunk_id &at[142]
#define ID2S_REPLY_GET_CHUNK_data &at[143]

#define ID2_REPLY_DATA &at[114]
#define ID2_REPLY_DATA_data_type &at[115]
#define ID2_REPLY_DATA_data_format &at[116]
#define ID2_REPLY_DATA_data_compression &at[117]
#define ID2_REPLY_DATA_data &at[118]
#define ID2_REPLY_DATA_data_E &at[119]

#define ID2_BLOB_SEQ_IDS &at[145]
#define ID2_BLOB_SEQ_IDS_E &at[146]

#define ID2_BLOB_SEQ_ID &at[147]
#define ID2_BLOB_SEQ_ID_seq_id &at[148]
#define ID2_BLOB_SEQ_ID_replaced &at[149]

#define ID2_PARAM &at[14]
#define ID2_PARAM_name &at[15]
#define ID2_PARAM_value &at[17]
#define ID2_PARAM_value_E &at[18]
#define ID2_PARAM_type &at[20]
