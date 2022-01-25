/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "allpub.h60";
static AsnValxNode avnx[58] = {
    {20,"unk" ,0,0.0,&avnx[1] } ,
    {20,"gt" ,1,0.0,&avnx[2] } ,
    {20,"lt" ,2,0.0,&avnx[3] } ,
    {20,"tr" ,3,0.0,&avnx[4] } ,
    {20,"tl" ,4,0.0,&avnx[5] } ,
    {20,"circle" ,5,0.0,&avnx[6] } ,
    {20,"other" ,255,0.0,NULL } ,
    {20,"primary" ,1,0.0,&avnx[8] } ,
    {20,"secondary" ,2,0.0,NULL } ,
    {20,"compiler" ,1,0.0,&avnx[10] } ,
    {20,"editor" ,2,0.0,&avnx[11] } ,
    {20,"patent-assignee" ,3,0.0,&avnx[12] } ,
    {20,"translator" ,4,0.0,NULL } ,
    {1,"ENG" ,0,0.0,NULL } ,
    {20,"submitted" ,1,0.0,&avnx[15] } ,
    {20,"in-press" ,2,0.0,&avnx[16] } ,
    {20,"other" ,255,0.0,NULL } ,
    {20,"retracted" ,1,0.0,&avnx[18] } ,
    {20,"notice" ,2,0.0,&avnx[19] } ,
    {20,"in-error" ,3,0.0,&avnx[20] } ,
    {20,"erratum" ,4,0.0,NULL } ,
    {20,"paper" ,1,0.0,&avnx[22] } ,
    {20,"tape" ,2,0.0,&avnx[23] } ,
    {20,"floppy" ,3,0.0,&avnx[24] } ,
    {20,"email" ,4,0.0,&avnx[25] } ,
    {20,"other" ,255,0.0,NULL } ,
    {2,NULL,0,0.0,NULL } ,
    {2,NULL,0,0.0,NULL } ,
    {20,"nameonly" ,0,0.0,&avnx[29] } ,
    {20,"cas" ,1,0.0,&avnx[30] } ,
    {20,"ec" ,2,0.0,NULL } ,
    {20,"ddbj" ,1,0.0,&avnx[32] } ,
    {20,"carbbank" ,2,0.0,&avnx[33] } ,
    {20,"embl" ,3,0.0,&avnx[34] } ,
    {20,"hdb" ,4,0.0,&avnx[35] } ,
    {20,"genbank" ,5,0.0,&avnx[36] } ,
    {20,"hgml" ,6,0.0,&avnx[37] } ,
    {20,"mim" ,7,0.0,&avnx[38] } ,
    {20,"msd" ,8,0.0,&avnx[39] } ,
    {20,"pdb" ,9,0.0,&avnx[40] } ,
    {20,"pir" ,10,0.0,&avnx[41] } ,
    {20,"prfseqdb" ,11,0.0,&avnx[42] } ,
    {20,"psd" ,12,0.0,&avnx[43] } ,
    {20,"swissprot" ,13,0.0,&avnx[44] } ,
    {20,"gdb" ,14,0.0,NULL } ,
    {20,"other" ,0,0.0,&avnx[46] } ,
    {20,"comment" ,1,0.0,&avnx[47] } ,
    {20,"erratum" ,2,0.0,NULL } ,
    {20,"medline" ,1,0.0,&avnx[49] } ,
    {20,"pubmed" ,2,0.0,&avnx[50] } ,
    {20,"ncbigi" ,3,0.0,NULL } ,
    {20,"publisher" ,1,0.0,&avnx[52] } ,
    {20,"premedline" ,2,0.0,&avnx[53] } ,
    {20,"medline" ,3,0.0,NULL } ,
    {3,NULL,3,0.0,NULL } ,
    {20,"manuscript" ,1,0.0,&avnx[56] } ,
    {20,"letter" ,2,0.0,&avnx[57] } ,
    {20,"thesis" ,3,0.0,NULL } };

static AsnType atx[291] = {
  {401, "Date" ,1,0,0,0,0,1,0,0,NULL,&atx[11],&atx[1],0,&atx[12]} ,
  {0, "str" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[3]} ,
  {323, "VisibleString" ,0,26,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "std" ,128,1,0,0,0,0,0,0,NULL,&atx[4],NULL,0,NULL} ,
  {407, "Date-std" ,1,0,0,0,0,0,0,0,NULL,&atx[10],&atx[5],0,&atx[21]} ,
  {0, "year" ,128,0,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[7]} ,
  {302, "INTEGER" ,0,2,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "month" ,128,1,0,1,0,0,0,0,NULL,&atx[6],NULL,0,&atx[8]} ,
  {0, "day" ,128,2,0,1,0,0,0,0,NULL,&atx[6],NULL,0,&atx[9]} ,
  {0, "season" ,128,3,0,1,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {311, "SEQUENCE" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {315, "CHOICE" ,0,-1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {402, "Person-id" ,1,0,0,0,0,1,0,0,NULL,&atx[11],&atx[13],0,&atx[17]} ,
  {0, "dbtag" ,128,0,0,0,0,0,0,0,NULL,&atx[14],NULL,0,&atx[20]} ,
  {404, "Dbtag" ,1,0,0,0,0,1,0,0,NULL,&atx[10],&atx[15],0,&atx[31]} ,
  {0, "db" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[16]} ,
  {0, "tag" ,128,1,0,0,0,0,0,0,NULL,&atx[17],NULL,0,NULL} ,
  {403, "Object-id" ,1,0,0,0,0,1,0,0,NULL,&atx[11],&atx[18],0,&atx[14]} ,
  {0, "id" ,128,0,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[19]} ,
  {0, "str" ,128,1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "name" ,128,1,0,0,0,0,0,0,NULL,&atx[21],NULL,0,&atx[29]} ,
  {408, "Name-std" ,1,0,0,0,0,0,0,0,NULL,&atx[10],&atx[22],0,&atx[47]} ,
  {0, "last" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[23]} ,
  {0, "first" ,128,1,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[24]} ,
  {0, "middle" ,128,2,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[25]} ,
  {0, "full" ,128,3,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[26]} ,
  {0, "initials" ,128,4,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[27]} ,
  {0, "suffix" ,128,5,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[28]} ,
  {0, "title" ,128,6,0,1,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "ml" ,128,2,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[30]} ,
  {0, "str" ,128,3,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {405, "Int-fuzz" ,1,0,0,0,0,1,0,0,NULL,&atx[11],&atx[32],0,&atx[42]} ,
  {0, "p-m" ,128,0,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[33]} ,
  {0, "range" ,128,1,0,0,0,0,0,0,NULL,&atx[10],&atx[34],0,&atx[36]} ,
  {0, "max" ,128,0,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[35]} ,
  {0, "min" ,128,1,0,0,0,0,0,0,NULL,&atx[6],NULL,0,NULL} ,
  {0, "pct" ,128,2,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[37]} ,
  {0, "lim" ,128,3,0,0,0,0,0,0,NULL,&atx[38],&avnx[0],0,&atx[39]} ,
  {310, "ENUMERATED" ,0,10,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "alt" ,128,4,0,0,0,0,0,0,NULL,&atx[41],&atx[40],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[6],NULL,0,NULL} ,
  {314, "SET OF" ,0,17,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {406, "User-object" ,1,0,0,0,0,1,0,0,NULL,&atx[10],&atx[43],0,&atx[4]} ,
  {0, "class" ,128,0,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[44]} ,
  {0, "type" ,128,1,0,0,0,0,0,0,NULL,&atx[17],NULL,0,&atx[45]} ,
  {0, "data" ,128,2,0,0,0,0,0,0,NULL,&atx[62],&atx[46],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[47],NULL,0,NULL} ,
  {409, "User-field" ,1,0,0,0,0,0,0,0,NULL,&atx[10],&atx[48],0,NULL} ,
  {0, "label" ,128,0,0,0,0,0,0,0,NULL,&atx[17],NULL,0,&atx[49]} ,
  {0, "num" ,128,1,0,1,0,0,0,0,NULL,&atx[6],NULL,0,&atx[50]} ,
  {0, "data" ,128,2,0,0,0,0,0,0,NULL,&atx[11],&atx[51],0,NULL} ,
  {0, "str" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[52]} ,
  {0, "int" ,128,1,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[53]} ,
  {0, "real" ,128,2,0,0,0,0,0,0,NULL,&atx[54],NULL,0,&atx[55]} ,
  {309, "REAL" ,0,9,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "bool" ,128,3,0,0,0,0,0,0,NULL,&atx[56],NULL,0,&atx[57]} ,
  {301, "BOOLEAN" ,0,1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "os" ,128,4,0,0,0,0,0,0,NULL,&atx[58],NULL,0,&atx[59]} ,
  {304, "OCTET STRING" ,0,4,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "object" ,128,5,0,0,0,0,0,0,NULL,&atx[42],NULL,0,&atx[60]} ,
  {0, "strs" ,128,6,0,0,0,0,0,0,NULL,&atx[62],&atx[61],0,&atx[63]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {312, "SEQUENCE OF" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "ints" ,128,7,0,0,0,0,0,0,NULL,&atx[62],&atx[64],0,&atx[65]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[6],NULL,0,NULL} ,
  {0, "reals" ,128,8,0,0,0,0,0,0,NULL,&atx[62],&atx[66],0,&atx[67]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[54],NULL,0,NULL} ,
  {0, "oss" ,128,9,0,0,0,0,0,0,NULL,&atx[62],&atx[68],0,&atx[69]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[58],NULL,0,NULL} ,
  {0, "fields" ,128,10,0,0,0,0,0,0,NULL,&atx[62],&atx[70],0,&atx[71]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[47],NULL,0,NULL} ,
  {0, "objects" ,128,11,0,0,0,0,0,0,NULL,&atx[62],&atx[72],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[42],NULL,0,NULL} ,
  {401, "Pub" ,1,0,0,0,0,1,0,0,NULL,&atx[11],&atx[74],0,&atx[276]} ,
  {0, "gen" ,128,0,0,0,0,0,0,0,NULL,&atx[75],NULL,0,&atx[131]} ,
  {411, "Cit-gen" ,1,0,0,0,0,0,1,0,NULL,&atx[76],NULL,0,&atx[266]} ,
  {407, "Cit-gen" ,1,0,0,0,0,1,0,0,NULL,&atx[10],&atx[77],0,&atx[178]} ,
  {0, "cit" ,128,0,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[78]} ,
  {0, "authors" ,128,1,0,1,0,0,0,0,NULL,&atx[79],NULL,0,&atx[108]} ,
  {415, "Auth-list" ,1,0,0,0,0,0,0,0,NULL,&atx[10],&atx[80],0,&atx[136]} ,
  {0, "names" ,128,0,0,0,0,0,0,0,NULL,&atx[11],&atx[81],0,&atx[107]} ,
  {0, "std" ,128,0,0,0,0,0,0,0,NULL,&atx[62],&atx[82],0,&atx[103]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[83],NULL,0,NULL} ,
  {411, "Author" ,1,0,0,0,0,1,0,0,NULL,&atx[10],&atx[84],0,&atx[130]} ,
  {0, "name" ,128,0,0,0,0,0,0,0,NULL,&atx[85],NULL,0,&atx[86]} ,
  {413, "Person-id" ,1,0,0,0,0,0,1,0,NULL,&atx[12],NULL,0,&atx[126]} ,
  {0, "level" ,128,1,0,1,0,0,0,0,NULL,&atx[38],&avnx[7],0,&atx[87]} ,
  {0, "role" ,128,2,0,1,0,0,0,0,NULL,&atx[38],&avnx[9],0,&atx[88]} ,
  {0, "affil" ,128,3,0,1,0,0,0,0,NULL,&atx[89],NULL,0,&atx[102]} ,
  {419, "Affil" ,1,0,0,0,0,0,0,0,NULL,&atx[11],&atx[90],0,&atx[149]} ,
  {0, "str" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[91]} ,
  {0, "std" ,128,1,0,0,0,0,0,0,NULL,&atx[10],&atx[92],0,NULL} ,
  {0, "affil" ,128,0,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[93]} ,
  {0, "div" ,128,1,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[94]} ,
  {0, "city" ,128,2,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[95]} ,
  {0, "sub" ,128,3,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[96]} ,
  {0, "country" ,128,4,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[97]} ,
  {0, "street" ,128,5,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[98]} ,
  {0, "email" ,128,6,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[99]} ,
  {0, "fax" ,128,7,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[100]} ,
  {0, "phone" ,128,8,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[101]} ,
  {0, "postal-code" ,128,9,0,1,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "is-corr" ,128,4,0,1,0,0,0,0,NULL,&atx[56],NULL,0,NULL} ,
  {0, "ml" ,128,1,0,0,0,0,0,0,NULL,&atx[62],&atx[104],0,&atx[105]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "str" ,128,2,0,0,0,0,0,0,NULL,&atx[62],&atx[106],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "affil" ,128,1,0,1,0,0,0,0,NULL,&atx[89],NULL,0,NULL} ,
  {0, "muid" ,128,2,0,1,0,0,0,0,NULL,&atx[6],NULL,0,&atx[109]} ,
  {0, "journal" ,128,3,0,1,0,0,0,0,NULL,&atx[110],NULL,0,&atx[122]} ,
  {410, "Title" ,1,0,0,0,0,1,0,0,NULL,&atx[41],&atx[111],0,&atx[83]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[11],&atx[112],0,NULL} ,
  {0, "name" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[113]} ,
  {0, "tsub" ,128,1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[114]} ,
  {0, "trans" ,128,2,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[115]} ,
  {0, "jta" ,128,3,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[116]} ,
  {0, "iso-jta" ,128,4,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[117]} ,
  {0, "ml-jta" ,128,5,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[118]} ,
  {0, "coden" ,128,6,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[119]} ,
  {0, "issn" ,128,7,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[120]} ,
  {0, "abr" ,128,8,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[121]} ,
  {0, "isbn" ,128,9,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "volume" ,128,4,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[123]} ,
  {0, "issue" ,128,5,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[124]} ,
  {0, "pages" ,128,6,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[125]} ,
  {0, "date" ,128,7,0,1,0,0,0,0,NULL,&atx[126],NULL,0,&atx[127]} ,
  {414, "Date" ,1,0,0,0,0,0,1,0,NULL,&atx[0],NULL,0,&atx[79]} ,
  {0, "serial-number" ,128,8,0,1,0,0,0,0,NULL,&atx[6],NULL,0,&atx[128]} ,
  {0, "title" ,128,9,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[129]} ,
  {0, "pmid" ,128,10,0,1,0,0,0,0,NULL,&atx[130],NULL,0,NULL} ,
  {412, "PubMedId" ,1,0,0,0,0,1,0,0,NULL,&atx[6],NULL,0,&atx[85]} ,
  {0, "sub" ,128,1,0,0,0,0,0,0,NULL,&atx[132],NULL,0,&atx[155]} ,
  {413, "Cit-sub" ,1,0,0,0,0,0,1,0,NULL,&atx[133],NULL,0,&atx[275]} ,
  {409, "Cit-sub" ,1,0,0,0,0,1,0,0,NULL,&atx[10],&atx[134],0,&atx[110]} ,
  {0, "authors" ,128,0,0,0,0,0,0,0,NULL,&atx[79],NULL,0,&atx[135]} ,
  {0, "imp" ,128,1,0,1,0,0,0,0,NULL,&atx[136],NULL,0,&atx[152]} ,
  {416, "Imprint" ,1,0,0,0,0,0,0,0,NULL,&atx[10],&atx[137],0,&atx[181]} ,
  {0, "date" ,128,0,0,0,0,0,0,0,NULL,&atx[126],NULL,0,&atx[138]} ,
  {0, "volume" ,128,1,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[139]} ,
  {0, "issue" ,128,2,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[140]} ,
  {0, "pages" ,128,3,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[141]} ,
  {0, "section" ,128,4,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[142]} ,
  {0, "pub" ,128,5,0,1,0,0,0,0,NULL,&atx[89],NULL,0,&atx[143]} ,
  {0, "cprt" ,128,6,0,1,0,0,0,0,NULL,&atx[126],NULL,0,&atx[144]} ,
  {0, "part-sup" ,128,7,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[145]} ,
  {0, "language" ,128,8,0,0,1,0,0,0,&avnx[13],&atx[2],NULL,0,&atx[146]} ,
  {0, "prepub" ,128,9,0,1,0,0,0,0,NULL,&atx[38],&avnx[14],0,&atx[147]} ,
  {0, "part-supi" ,128,10,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[148]} ,
  {0, "retract" ,128,11,0,1,0,0,0,0,NULL,&atx[149],NULL,0,NULL} ,
  {420, "CitRetract" ,1,0,0,0,0,0,0,0,NULL,&atx[10],&atx[150],0,NULL} ,
  {0, "type" ,128,0,0,0,0,0,0,0,NULL,&atx[38],&avnx[17],0,&atx[151]} ,
  {0, "exp" ,128,1,0,1,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "medium" ,128,2,0,1,0,0,0,0,NULL,&atx[38],&avnx[21],0,&atx[153]} ,
  {0, "date" ,128,3,0,1,0,0,0,0,NULL,&atx[126],NULL,0,&atx[154]} ,
  {0, "descr" ,128,4,0,1,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "medline" ,128,2,0,0,0,0,0,0,NULL,&atx[156],NULL,0,&atx[226]} ,
  {404, "Medline-entry" ,1,0,0,0,0,0,1,0,NULL,&atx[157],NULL,0,&atx[228]} ,
  {401, "Medline-entry" ,1,0,0,0,0,1,0,0,NULL,&atx[10],&atx[158],0,&atx[204]} ,
  {0, "uid" ,128,0,0,1,0,0,0,0,NULL,&atx[6],NULL,0,&atx[159]} ,
  {0, "em" ,128,1,0,0,0,0,0,0,NULL,&atx[160],NULL,0,&atx[161]} ,
  {405, "Date" ,1,0,0,0,0,0,1,0,NULL,&atx[0],NULL,0,&atx[188]} ,
  {0, "cit" ,128,2,0,0,0,0,0,0,NULL,&atx[162],NULL,0,&atx[185]} ,
  {403, "Cit-art" ,1,0,0,0,0,0,1,0,NULL,&atx[163],NULL,0,&atx[212]} ,
  {401, "Cit-art" ,1,0,0,0,0,1,0,0,NULL,&atx[10],&atx[164],0,&atx[168]} ,
  {0, "title" ,128,0,0,1,0,0,0,0,NULL,&atx[110],NULL,0,&atx[165]} ,
  {0, "authors" ,128,1,0,1,0,0,0,0,NULL,&atx[79],NULL,0,&atx[166]} ,
  {0, "from" ,128,2,0,0,0,0,0,0,NULL,&atx[11],&atx[167],0,NULL} ,
  {0, "journal" ,128,0,0,0,0,0,0,0,NULL,&atx[168],NULL,0,&atx[171]} ,
  {402, "Cit-jour" ,1,0,0,0,0,1,0,0,NULL,&atx[10],&atx[169],0,&atx[172]} ,
  {0, "title" ,128,0,0,0,0,0,0,0,NULL,&atx[110],NULL,0,&atx[170]} ,
  {0, "imp" ,128,1,0,0,0,0,0,0,NULL,&atx[136],NULL,0,NULL} ,
  {0, "book" ,128,1,0,0,0,0,0,0,NULL,&atx[172],NULL,0,&atx[177]} ,
  {403, "Cit-book" ,1,0,0,0,0,1,0,0,NULL,&atx[10],&atx[173],0,&atx[237]} ,
  {0, "title" ,128,0,0,0,0,0,0,0,NULL,&atx[110],NULL,0,&atx[174]} ,
  {0, "coll" ,128,1,0,1,0,0,0,0,NULL,&atx[110],NULL,0,&atx[175]} ,
  {0, "authors" ,128,2,0,0,0,0,0,0,NULL,&atx[79],NULL,0,&atx[176]} ,
  {0, "imp" ,128,3,0,0,0,0,0,0,NULL,&atx[136],NULL,0,NULL} ,
  {0, "proc" ,128,2,0,0,0,0,0,0,NULL,&atx[178],NULL,0,NULL} ,
  {408, "Cit-proc" ,1,0,0,0,0,1,0,0,NULL,&atx[10],&atx[179],0,&atx[133]} ,
  {0, "book" ,128,0,0,0,0,0,0,0,NULL,&atx[172],NULL,0,&atx[180]} ,
  {0, "meet" ,128,1,0,0,0,0,0,0,NULL,&atx[181],NULL,0,NULL} ,
  {417, "Meeting" ,1,0,0,0,0,0,0,0,NULL,&atx[10],&atx[182],0,&atx[252]} ,
  {0, "number" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[183]} ,
  {0, "date" ,128,1,0,0,0,0,0,0,NULL,&atx[126],NULL,0,&atx[184]} ,
  {0, "place" ,128,2,0,1,0,0,0,0,NULL,&atx[89],NULL,0,NULL} ,
  {0, "abstract" ,128,3,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[186]} ,
  {0, "mesh" ,128,4,0,1,0,0,0,0,NULL,&atx[41],&atx[187],0,&atx[196]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[188],NULL,0,NULL} ,
  {406, "Medline-mesh" ,1,0,0,0,0,0,0,0,NULL,&atx[10],&atx[189],0,&atx[198]} ,
  {0, "mp" ,128,0,0,0,1,0,0,0,&avnx[26],&atx[56],NULL,0,&atx[190]} ,
  {0, "term" ,128,1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[191]} ,
  {0, "qual" ,128,2,0,1,0,0,0,0,NULL,&atx[41],&atx[192],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[193],NULL,0,NULL} ,
  {409, "Medline-qual" ,1,0,0,0,0,0,0,0,NULL,&atx[10],&atx[194],0,&atx[222]} ,
  {0, "mp" ,128,0,0,0,1,0,0,0,&avnx[27],&atx[56],NULL,0,&atx[195]} ,
  {0, "subh" ,128,1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "substance" ,128,5,0,1,0,0,0,0,NULL,&atx[41],&atx[197],0,&atx[202]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[198],NULL,0,NULL} ,
  {407, "Medline-rn" ,1,0,0,0,0,0,0,0,NULL,&atx[10],&atx[199],0,&atx[217]} ,
  {0, "type" ,128,0,0,0,0,0,0,0,NULL,&atx[38],&avnx[28],0,&atx[200]} ,
  {0, "cit" ,128,1,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[201]} ,
  {0, "name" ,128,2,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "xref" ,128,6,0,1,0,0,0,0,NULL,&atx[41],&atx[203],0,&atx[207]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[204],NULL,0,NULL} ,
  {402, "Medline-si" ,1,0,0,0,0,1,0,0,NULL,&atx[10],&atx[205],0,&atx[162]} ,
  {0, "type" ,128,0,0,0,0,0,0,0,NULL,&atx[38],&avnx[31],0,&atx[206]} ,
  {0, "cit" ,128,1,0,1,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "idnum" ,128,7,0,1,0,0,0,0,NULL,&atx[41],&atx[208],0,&atx[209]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "gene" ,128,8,0,1,0,0,0,0,NULL,&atx[41],&atx[210],0,&atx[211]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "pmid" ,128,9,0,1,0,0,0,0,NULL,&atx[212],NULL,0,&atx[213]} ,
  {404, "PubMedId" ,1,0,0,0,0,0,1,0,NULL,&atx[130],NULL,0,&atx[160]} ,
  {0, "pub-type" ,128,10,0,1,0,0,0,0,NULL,&atx[41],&atx[214],0,&atx[215]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "mlfield" ,128,11,0,1,0,0,0,0,NULL,&atx[41],&atx[216],0,&atx[225]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[217],NULL,0,NULL} ,
  {408, "Medline-field" ,1,0,0,0,0,0,0,0,NULL,&atx[10],&atx[218],0,&atx[193]} ,
  {0, "type" ,128,0,0,0,0,0,0,0,NULL,&atx[6],&avnx[45],0,&atx[219]} ,
  {0, "str" ,128,1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[220]} ,
  {0, "ids" ,128,2,0,1,0,0,0,0,NULL,&atx[62],&atx[221],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[222],NULL,0,NULL} ,
  {410, "DocRef" ,1,0,0,0,0,0,0,0,NULL,&atx[10],&atx[223],0,NULL} ,
  {0, "type" ,128,0,0,0,0,0,0,0,NULL,&atx[6],&avnx[48],0,&atx[224]} ,
  {0, "uid" ,128,1,0,0,0,0,0,0,NULL,&atx[6],NULL,0,NULL} ,
  {0, "status" ,128,12,0,0,1,0,0,0,&avnx[54],&atx[6],&avnx[51],0,NULL} ,
  {0, "muid" ,128,3,0,0,0,0,0,0,NULL,&atx[6],NULL,0,&atx[227]} ,
  {0, "article" ,128,4,0,0,0,0,0,0,NULL,&atx[228],NULL,0,&atx[229]} ,
  {405, "Cit-art" ,1,0,0,0,0,0,1,0,NULL,&atx[163],NULL,0,&atx[230]} ,
  {0, "journal" ,128,5,0,0,0,0,0,0,NULL,&atx[230],NULL,0,&atx[231]} ,
  {406, "Cit-jour" ,1,0,0,0,0,0,1,0,NULL,&atx[168],NULL,0,&atx[232]} ,
  {0, "book" ,128,6,0,0,0,0,0,0,NULL,&atx[232],NULL,0,&atx[233]} ,
  {407, "Cit-book" ,1,0,0,0,0,0,1,0,NULL,&atx[172],NULL,0,&atx[234]} ,
  {0, "proc" ,128,7,0,0,0,0,0,0,NULL,&atx[234],NULL,0,&atx[235]} ,
  {408, "Cit-proc" ,1,0,0,0,0,0,1,0,NULL,&atx[178],NULL,0,&atx[236]} ,
  {0, "patent" ,128,8,0,0,0,0,0,0,NULL,&atx[236],NULL,0,&atx[257]} ,
  {409, "Cit-pat" ,1,0,0,0,0,0,1,0,NULL,&atx[237],NULL,0,&atx[258]} ,
  {404, "Cit-pat" ,1,0,0,0,0,1,0,0,NULL,&atx[10],&atx[238],0,&atx[267]} ,
  {0, "title" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[239]} ,
  {0, "authors" ,128,1,0,0,0,0,0,0,NULL,&atx[79],NULL,0,&atx[240]} ,
  {0, "country" ,128,2,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[241]} ,
  {0, "doc-type" ,128,3,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[242]} ,
  {0, "number" ,128,4,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[243]} ,
  {0, "date-issue" ,128,5,0,1,0,0,0,0,NULL,&atx[126],NULL,0,&atx[244]} ,
  {0, "class" ,128,6,0,1,0,0,0,0,NULL,&atx[62],&atx[245],0,&atx[246]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "app-number" ,128,7,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[247]} ,
  {0, "app-date" ,128,8,0,1,0,0,0,0,NULL,&atx[126],NULL,0,&atx[248]} ,
  {0, "applicants" ,128,9,0,1,0,0,0,0,NULL,&atx[79],NULL,0,&atx[249]} ,
  {0, "assignees" ,128,10,0,1,0,0,0,0,NULL,&atx[79],NULL,0,&atx[250]} ,
  {0, "priority" ,128,11,0,1,0,0,0,0,NULL,&atx[62],&atx[251],0,&atx[256]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[252],NULL,0,NULL} ,
  {418, "Patent-priority" ,1,0,0,0,0,0,0,0,NULL,&atx[10],&atx[253],0,&atx[89]} ,
  {0, "country" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[254]} ,
  {0, "number" ,128,1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[255]} ,
  {0, "date" ,128,2,0,0,0,0,0,0,NULL,&atx[126],NULL,0,NULL} ,
  {0, "abstract" ,128,12,0,1,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "pat-id" ,128,9,0,0,0,0,0,0,NULL,&atx[258],NULL,0,&atx[265]} ,
  {410, "Id-pat" ,1,0,0,0,0,0,1,0,NULL,&atx[259],NULL,0,&atx[75]} ,
  {406, "Id-pat" ,1,0,0,0,0,1,0,0,NULL,&atx[10],&atx[260],0,&atx[76]} ,
  {0, "country" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[261]} ,
  {0, "id" ,128,1,0,0,0,0,0,0,NULL,&atx[11],&atx[262],0,&atx[264]} ,
  {0, "number" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[263]} ,
  {0, "app-number" ,128,1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "doc-type" ,128,2,0,1,0,0,0,0,NULL,&atx[2],NULL,0,NULL} ,
  {0, "man" ,128,10,0,0,0,0,0,0,NULL,&atx[266],NULL,0,&atx[271]} ,
  {412, "Cit-let" ,1,0,0,0,0,0,1,0,NULL,&atx[267],NULL,0,&atx[132]} ,
  {405, "Cit-let" ,1,0,0,0,0,1,0,0,NULL,&atx[10],&atx[268],0,&atx[259]} ,
  {0, "cit" ,128,0,0,0,0,0,0,0,NULL,&atx[172],NULL,0,&atx[269]} ,
  {0, "man-id" ,128,1,0,1,0,0,0,0,NULL,&atx[2],NULL,0,&atx[270]} ,
  {0, "type" ,128,2,0,1,0,0,0,0,NULL,&atx[38],&avnx[55],0,NULL} ,
  {0, "equiv" ,128,11,0,0,0,0,0,0,NULL,&atx[272],NULL,0,&atx[274]} ,
  {403, "Pub-equiv" ,1,0,0,0,0,1,0,0,NULL,&atx[41],&atx[273],0,&atx[156]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[73],NULL,0,NULL} ,
  {0, "pmid" ,128,12,0,0,0,0,0,0,NULL,&atx[275],NULL,0,NULL} ,
  {414, "PubMedId" ,1,0,0,0,0,0,1,0,NULL,&atx[130],NULL,0,NULL} ,
  {402, "Pub-set" ,1,0,0,0,0,1,0,0,NULL,&atx[11],&atx[277],0,&atx[272]} ,
  {0, "pub" ,128,0,0,0,0,0,0,0,NULL,&atx[41],&atx[278],0,&atx[279]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[73],NULL,0,NULL} ,
  {0, "medline" ,128,1,0,0,0,0,0,0,NULL,&atx[41],&atx[280],0,&atx[281]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[156],NULL,0,NULL} ,
  {0, "article" ,128,2,0,0,0,0,0,0,NULL,&atx[41],&atx[282],0,&atx[283]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[228],NULL,0,NULL} ,
  {0, "journal" ,128,3,0,0,0,0,0,0,NULL,&atx[41],&atx[284],0,&atx[285]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[230],NULL,0,NULL} ,
  {0, "book" ,128,4,0,0,0,0,0,0,NULL,&atx[41],&atx[286],0,&atx[287]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[232],NULL,0,NULL} ,
  {0, "proc" ,128,5,0,0,0,0,0,0,NULL,&atx[41],&atx[288],0,&atx[289]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[234],NULL,0,NULL} ,
  {0, "patent" ,128,6,0,0,0,0,0,0,NULL,&atx[41],&atx[290],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[236],NULL,0,NULL} };

static AsnModule ampx[4] = {
  { "NCBI-General" , "allpub.h60",&atx[0],NULL,&ampx[1],0,0} ,
  { "NCBI-Pub" , NULL,&atx[73],NULL,&ampx[2],0,0} ,
  { "NCBI-Biblio" , NULL,&atx[163],NULL,&ampx[3],0,0} ,
  { "NCBI-Medline" , NULL,&atx[157],NULL,NULL,0,0} };

static AsnValxNodePtr avn = avnx;
static AsnTypePtr at = atx;
static AsnModulePtr amp = ampx;



/**************************************************
*
*    Defines for Module NCBI-General
*
**************************************************/

#define DATE &at[0]
#define DATE_str &at[1]
#define DATE_std &at[3]

#define PERSON_ID &at[12]
#define PERSON_ID_dbtag &at[13]
#define PERSON_ID_name &at[20]
#define PERSON_ID_ml &at[29]
#define PERSON_ID_str &at[30]

#define OBJECT_ID &at[17]
#define OBJECT_ID_id &at[18]
#define OBJECT_ID_str &at[19]

#define DBTAG &at[14]
#define DBTAG_db &at[15]
#define DBTAG_tag &at[16]

#define INT_FUZZ &at[31]
#define INT_FUZZ_p_m &at[32]
#define INT_FUZZ_range &at[33]
#define INT_FUZZ_range_max &at[34]
#define INT_FUZZ_range_min &at[35]
#define INT_FUZZ_pct &at[36]
#define INT_FUZZ_lim &at[37]
#define INT_FUZZ_alt &at[39]
#define INT_FUZZ_alt_E &at[40]

#define USER_OBJECT &at[42]
#define USER_OBJECT_class &at[43]
#define USER_OBJECT_type &at[44]
#define USER_OBJECT_data &at[45]
#define USER_OBJECT_data_E &at[46]

#define DATE_STD &at[4]
#define DATE_STD_year &at[5]
#define DATE_STD_month &at[7]
#define DATE_STD_day &at[8]
#define DATE_STD_season &at[9]

#define NAME_STD &at[21]
#define NAME_STD_last &at[22]
#define NAME_STD_first &at[23]
#define NAME_STD_middle &at[24]
#define NAME_STD_full &at[25]
#define NAME_STD_initials &at[26]
#define NAME_STD_suffix &at[27]
#define NAME_STD_title &at[28]

#define USER_FIELD &at[47]
#define USER_FIELD_label &at[48]
#define USER_FIELD_num &at[49]
#define USER_FIELD_data &at[50]
#define USER_FIELD_data_str &at[51]
#define USER_FIELD_data_int &at[52]
#define USER_FIELD_data_real &at[53]
#define USER_FIELD_data_bool &at[55]
#define USER_FIELD_data_os &at[57]
#define USER_FIELD_data_object &at[59]
#define USER_FIELD_data_strs &at[60]
#define USER_FIELD_data_strs_E &at[61]
#define USER_FIELD_data_ints &at[63]
#define USER_FIELD_data_ints_E &at[64]
#define USER_FIELD_data_reals &at[65]
#define USER_FIELD_data_reals_E &at[66]
#define USER_FIELD_data_oss &at[67]
#define USER_FIELD_data_oss_E &at[68]
#define USER_FIELD_data_fields &at[69]
#define USER_FIELD_data_fields_E &at[70]
#define USER_FIELD_data_objects &at[71]
#define USER_FIELD_data_objects_E &at[72]


/**************************************************
*
*    Defines for Module NCBI-Pub
*
**************************************************/

#define PUB &at[73]
#define PUB_gen &at[74]
#define PUB_sub &at[131]
#define PUB_medline &at[155]
#define PUB_muid &at[226]
#define PUB_article &at[227]
#define PUB_journal &at[229]
#define PUB_book &at[231]
#define PUB_proc &at[233]
#define PUB_patent &at[235]
#define PUB_pat_id &at[257]
#define PUB_man &at[265]
#define PUB_equiv &at[271]
#define PUB_pmid &at[274]

#define PUB_SET &at[276]
#define PUB_SET_pub &at[277]
#define PUB_SET_pub_E &at[278]
#define PUB_SET_medline &at[279]
#define PUB_SET_medline_E &at[280]
#define PUB_SET_article &at[281]
#define PUB_SET_article_E &at[282]
#define PUB_SET_journal &at[283]
#define PUB_SET_journal_E &at[284]
#define PUB_SET_book &at[285]
#define PUB_SET_book_E &at[286]
#define PUB_SET_proc &at[287]
#define PUB_SET_proc_E &at[288]
#define PUB_SET_patent &at[289]
#define PUB_SET_patent_E &at[290]

#define PUB_EQUIV &at[272]
#define PUB_EQUIV_E &at[273]


/**************************************************
*
*    Defines for Module NCBI-Biblio
*
**************************************************/

#define CIT_ART &at[163]
#define CIT_ART_title &at[164]
#define CIT_ART_authors &at[165]
#define CIT_ART_from &at[166]
#define CIT_ART_from_journal &at[167]
#define CIT_ART_from_book &at[171]
#define CIT_ART_from_proc &at[177]

#define CIT_JOUR &at[168]
#define CIT_JOUR_title &at[169]
#define CIT_JOUR_imp &at[170]

#define CIT_BOOK &at[172]
#define CIT_BOOK_title &at[173]
#define CIT_BOOK_coll &at[174]
#define CIT_BOOK_authors &at[175]
#define CIT_BOOK_imp &at[176]

#define CIT_PAT &at[237]
#define CIT_PAT_title &at[238]
#define CIT_PAT_authors &at[239]
#define CIT_PAT_country &at[240]
#define CIT_PAT_doc_type &at[241]
#define CIT_PAT_number &at[242]
#define CIT_PAT_date_issue &at[243]
#define CIT_PAT_class &at[244]
#define CIT_PAT_class_E &at[245]
#define CIT_PAT_app_number &at[246]
#define CIT_PAT_app_date &at[247]
#define CIT_PAT_applicants &at[248]
#define CIT_PAT_assignees &at[249]
#define CIT_PAT_priority &at[250]
#define CIT_PAT_priority_E &at[251]
#define CIT_PAT_abstract &at[256]

#define CIT_LET &at[267]
#define CIT_LET_cit &at[268]
#define CIT_LET_man_id &at[269]
#define CIT_LET_type &at[270]

#define ID_PAT &at[259]
#define ID_PAT_country &at[260]
#define ID_PAT_id &at[261]
#define ID_PAT_id_number &at[262]
#define ID_PAT_id_app_number &at[263]
#define ID_PAT_doc_type &at[264]

#define CIT_GEN &at[76]
#define CIT_GEN_cit &at[77]
#define CIT_GEN_authors &at[78]
#define CIT_GEN_muid &at[108]
#define CIT_GEN_journal &at[109]
#define CIT_GEN_volume &at[122]
#define CIT_GEN_issue &at[123]
#define CIT_GEN_pages &at[124]
#define CIT_GEN_date &at[125]
#define CIT_GEN_serial_number &at[127]
#define CIT_GEN_title &at[128]
#define CIT_GEN_pmid &at[129]

#define CIT_PROC &at[178]
#define CIT_PROC_book &at[179]
#define CIT_PROC_meet &at[180]

#define CIT_SUB &at[133]
#define CIT_SUB_authors &at[134]
#define CIT_SUB_imp &at[135]
#define CIT_SUB_medium &at[152]
#define CIT_SUB_date &at[153]
#define CIT_SUB_descr &at[154]

#define TITLE &at[110]
#define TITLE_E &at[111]
#define TITLE_E_name &at[112]
#define TITLE_E_tsub &at[113]
#define TITLE_E_trans &at[114]
#define TITLE_E_jta &at[115]
#define TITLE_E_iso_jta &at[116]
#define TITLE_E_ml_jta &at[117]
#define TITLE_E_coden &at[118]
#define TITLE_E_issn &at[119]
#define TITLE_E_abr &at[120]
#define TITLE_E_isbn &at[121]

#define AUTHOR &at[83]
#define AUTHOR_name &at[84]
#define AUTHOR_level &at[86]
#define AUTHOR_role &at[87]
#define AUTHOR_affil &at[88]
#define AUTHOR_is_corr &at[102]

#define PUBMEDID &at[130]

#define AUTH_LIST &at[79]
#define AUTH_LIST_names &at[80]
#define AUTH_LIST_names_std &at[81]
#define AUTH_LIST_names_std_E &at[82]
#define AUTH_LIST_names_ml &at[103]
#define AUTH_LIST_names_ml_E &at[104]
#define AUTH_LIST_names_str &at[105]
#define AUTH_LIST_names_str_E &at[106]
#define AUTH_LIST_affil &at[107]

#define IMPRINT &at[136]
#define IMPRINT_date &at[137]
#define IMPRINT_volume &at[138]
#define IMPRINT_issue &at[139]
#define IMPRINT_pages &at[140]
#define IMPRINT_section &at[141]
#define IMPRINT_pub &at[142]
#define IMPRINT_cprt &at[143]
#define IMPRINT_part_sup &at[144]
#define IMPRINT_language &at[145]
#define IMPRINT_prepub &at[146]
#define IMPRINT_part_supi &at[147]
#define IMPRINT_retract &at[148]

#define MEETING &at[181]
#define MEETING_number &at[182]
#define MEETING_date &at[183]
#define MEETING_place &at[184]

#define PATENT_PRIORITY &at[252]
#define PATENT_PRIORITY_country &at[253]
#define PATENT_PRIORITY_number &at[254]
#define PATENT_PRIORITY_date &at[255]

#define AFFIL &at[89]
#define AFFIL_str &at[90]
#define AFFIL_std &at[91]
#define AFFIL_std_affil &at[92]
#define AFFIL_std_div &at[93]
#define AFFIL_std_city &at[94]
#define AFFIL_std_sub &at[95]
#define AFFIL_std_country &at[96]
#define AFFIL_std_street &at[97]
#define AFFIL_std_email &at[98]
#define AFFIL_std_fax &at[99]
#define AFFIL_std_phone &at[100]
#define AFFIL_std_postal_code &at[101]

#define CITRETRACT &at[149]
#define CITRETRACT_type &at[150]
#define CITRETRACT_exp &at[151]


/**************************************************
*
*    Defines for Module NCBI-Medline
*
**************************************************/

#define MEDLINE_ENTRY &at[157]
#define MEDLINE_ENTRY_uid &at[158]
#define MEDLINE_ENTRY_em &at[159]
#define MEDLINE_ENTRY_cit &at[161]
#define MEDLINE_ENTRY_abstract &at[185]
#define MEDLINE_ENTRY_mesh &at[186]
#define MEDLINE_ENTRY_mesh_E &at[187]
#define MEDLINE_ENTRY_substance &at[196]
#define MEDLINE_ENTRY_substance_E &at[197]
#define MEDLINE_ENTRY_xref &at[202]
#define MEDLINE_ENTRY_xref_E &at[203]
#define MEDLINE_ENTRY_idnum &at[207]
#define MEDLINE_ENTRY_idnum_E &at[208]
#define MEDLINE_ENTRY_gene &at[209]
#define MEDLINE_ENTRY_gene_E &at[210]
#define MEDLINE_ENTRY_pmid &at[211]
#define MEDLINE_ENTRY_pub_type &at[213]
#define MEDLINE_ENTRY_pub_type_E &at[214]
#define MEDLINE_ENTRY_mlfield &at[215]
#define MEDLINE_ENTRY_mlfield_E &at[216]
#define MEDLINE_ENTRY_status &at[225]

#define MEDLINE_SI &at[204]
#define MEDLINE_SI_type &at[205]
#define MEDLINE_SI_cit &at[206]

#define MEDLINE_MESH &at[188]
#define MEDLINE_MESH_mp &at[189]
#define MEDLINE_MESH_term &at[190]
#define MEDLINE_MESH_qual &at[191]
#define MEDLINE_MESH_qual_E &at[192]

#define MEDLINE_RN &at[198]
#define MEDLINE_RN_type &at[199]
#define MEDLINE_RN_cit &at[200]
#define MEDLINE_RN_name &at[201]

#define MEDLINE_FIELD &at[217]
#define MEDLINE_FIELD_type &at[218]
#define MEDLINE_FIELD_str &at[219]
#define MEDLINE_FIELD_ids &at[220]
#define MEDLINE_FIELD_ids_E &at[221]

#define MEDLINE_QUAL &at[193]
#define MEDLINE_QUAL_mp &at[194]
#define MEDLINE_QUAL_subh &at[195]

#define DOCREF &at[222]
#define DOCREF_type &at[223]
#define DOCREF_uid &at[224]
