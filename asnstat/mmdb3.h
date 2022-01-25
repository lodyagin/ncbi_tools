/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "mmdb3.h61";
static AsnValxNode avnx[53] = {
    {20,"helix" ,1,0.0,&avnx[1] } ,
    {20,"strand" ,2,0.0,&avnx[2] } ,
    {20,"sheet" ,3,0.0,&avnx[3] } ,
    {20,"turn" ,4,0.0,&avnx[4] } ,
    {20,"site" ,5,0.0,&avnx[5] } ,
    {20,"footnote" ,6,0.0,&avnx[6] } ,
    {20,"comment" ,7,0.0,&avnx[7] } ,
    {20,"subgraph" ,100,0.0,&avnx[8] } ,
    {20,"region" ,101,0.0,&avnx[9] } ,
    {20,"core" ,102,0.0,&avnx[10] } ,
    {20,"supercore" ,103,0.0,&avnx[11] } ,
    {20,"color" ,150,0.0,&avnx[12] } ,
    {20,"render" ,151,0.0,&avnx[13] } ,
    {20,"label" ,152,0.0,&avnx[14] } ,
    {20,"transform" ,153,0.0,&avnx[15] } ,
    {20,"camera" ,154,0.0,&avnx[16] } ,
    {20,"script" ,155,0.0,&avnx[17] } ,
    {20,"alignment" ,200,0.0,&avnx[18] } ,
    {20,"similarity" ,201,0.0,&avnx[19] } ,
    {20,"multalign" ,202,0.0,&avnx[20] } ,
    {20,"indirect" ,203,0.0,&avnx[21] } ,
    {20,"cn3dstate" ,254,0.0,&avnx[22] } ,
    {20,"other" ,255,0.0,NULL } ,
    {20,"default" ,0,0.0,&avnx[24] } ,
    {20,"wire" ,1,0.0,&avnx[25] } ,
    {20,"space" ,2,0.0,&avnx[26] } ,
    {20,"stick" ,3,0.0,&avnx[27] } ,
    {20,"ballNStick" ,4,0.0,&avnx[28] } ,
    {20,"thickWire" ,5,0.0,&avnx[29] } ,
    {20,"hide" ,9,0.0,&avnx[30] } ,
    {20,"name" ,10,0.0,&avnx[31] } ,
    {20,"number" ,11,0.0,&avnx[32] } ,
    {20,"pdbNumber" ,12,0.0,&avnx[33] } ,
    {20,"objWireFrame" ,150,0.0,&avnx[34] } ,
    {20,"objPolygons" ,151,0.0,&avnx[35] } ,
    {20,"colorsetCPK" ,225,0.0,&avnx[36] } ,
    {20,"colorsetbyChain" ,226,0.0,&avnx[37] } ,
    {20,"colorsetbyTemp" ,227,0.0,&avnx[38] } ,
    {20,"colorsetbyRes" ,228,0.0,&avnx[39] } ,
    {20,"colorsetbyLen" ,229,0.0,&avnx[40] } ,
    {20,"colorsetbySStru" ,230,0.0,&avnx[41] } ,
    {20,"colorsetbyHydro" ,231,0.0,&avnx[42] } ,
    {20,"colorsetbyObject" ,246,0.0,&avnx[43] } ,
    {20,"colorsetbyDomain" ,247,0.0,&avnx[44] } ,
    {20,"other" ,255,0.0,NULL } ,
    {20,"examine" ,1,0.0,&avnx[46] } ,
    {20,"fly" ,2,0.0,&avnx[47] } ,
    {20,"walk" ,3,0.0,&avnx[48] } ,
    {20,"free" ,4,0.0,&avnx[49] } ,
    {20,"other" ,255,0.0,NULL } ,
    {3,NULL,10,0.0,NULL } ,
    {3,NULL,2,0.0,NULL } ,
    {3,NULL,2,0.0,NULL } };

static AsnType atx[206] = {
  {401, "Biostruc-feature-set" ,1,0,0,0,0,1,0,0,NULL,&atx[29],&atx[1],0,&atx[100]} ,
  {0, "id" ,128,0,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[4]} ,
  {410, "Biostruc-feature-set-id" ,1,0,0,0,0,1,0,0,NULL,&atx[3],NULL,0,&atx[19]} ,
  {302, "INTEGER" ,0,2,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "descr" ,128,1,0,1,0,0,0,0,NULL,&atx[14],&atx[5],0,&atx[15]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[6],NULL,0,NULL} ,
  {420, "Biostruc-feature-set-descr" ,1,0,0,0,0,0,0,0,NULL,&atx[13],&atx[7],0,&atx[17]} ,
  {0, "name" ,128,0,0,0,0,0,0,0,NULL,&atx[8],NULL,0,&atx[9]} ,
  {323, "VisibleString" ,0,26,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "pdb-comment" ,128,1,0,0,0,0,0,0,NULL,&atx[8],NULL,0,&atx[10]} ,
  {0, "other-comment" ,128,2,0,0,0,0,0,0,NULL,&atx[8],NULL,0,&atx[11]} ,
  {0, "attribution" ,128,3,0,0,0,0,0,0,NULL,&atx[12],NULL,0,NULL} ,
  {419, "Pub" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[6]} ,
  {315, "CHOICE" ,0,-1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {312, "SEQUENCE OF" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "features" ,128,2,0,0,0,0,0,0,NULL,&atx[14],&atx[16],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[17],NULL,0,NULL} ,
  {421, "Biostruc-feature" ,1,0,0,0,0,0,0,0,NULL,&atx[29],&atx[18],0,&atx[24]} ,
  {0, "id" ,128,0,0,0,0,0,0,0,NULL,&atx[19],NULL,0,&atx[20]} ,
  {411, "Biostruc-feature-id" ,1,0,0,0,0,1,0,0,NULL,&atx[3],NULL,0,&atx[87]} ,
  {0, "name" ,128,1,0,1,0,0,0,0,NULL,&atx[8],NULL,0,&atx[21]} ,
  {0, "type" ,128,2,0,1,0,0,0,0,NULL,&atx[3],&avnx[0],0,&atx[22]} ,
  {0, "property" ,128,3,0,1,0,0,0,0,NULL,&atx[13],&atx[23],0,&atx[98]} ,
  {0, "color" ,128,0,0,0,0,0,0,0,NULL,&atx[24],NULL,0,&atx[30]} ,
  {422, "Color-prop" ,1,0,0,0,0,0,0,0,NULL,&atx[29],&atx[25],0,&atx[31]} ,
  {0, "r" ,128,0,0,1,0,0,0,0,NULL,&atx[3],NULL,0,&atx[26]} ,
  {0, "g" ,128,1,0,1,0,0,0,0,NULL,&atx[3],NULL,0,&atx[27]} ,
  {0, "b" ,128,2,0,1,0,0,0,0,NULL,&atx[3],NULL,0,&atx[28]} ,
  {0, "name" ,128,3,0,1,0,0,0,0,NULL,&atx[8],NULL,0,NULL} ,
  {311, "SEQUENCE" ,0,16,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "render" ,128,1,0,0,0,0,0,0,NULL,&atx[31],NULL,0,&atx[32]} ,
  {423, "Render-prop" ,1,0,0,0,0,0,0,0,NULL,&atx[3],&avnx[23],0,&atx[57]} ,
  {0, "transform" ,128,2,0,0,0,0,0,0,NULL,&atx[33],NULL,0,&atx[56]} ,
  {409, "Transform" ,1,0,0,0,0,1,0,0,NULL,&atx[29],&atx[34],0,&atx[2]} ,
  {0, "id" ,128,0,0,0,0,0,0,0,NULL,&atx[3],NULL,0,&atx[35]} ,
  {0, "moves" ,128,1,0,0,0,0,0,0,NULL,&atx[14],&atx[36],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[37],NULL,0,NULL} ,
  {438, "Move" ,1,0,0,0,0,0,0,0,NULL,&atx[13],&atx[38],0,&atx[39]} ,
  {0, "rotate" ,128,0,0,0,0,0,0,0,NULL,&atx[39],NULL,0,&atx[50]} ,
  {439, "Rot-matrix" ,1,0,0,0,0,0,0,0,NULL,&atx[29],&atx[40],0,&atx[51]} ,
  {0, "scale-factor" ,128,0,0,0,0,0,0,0,NULL,&atx[3],NULL,0,&atx[41]} ,
  {0, "rot-11" ,128,1,0,0,0,0,0,0,NULL,&atx[3],NULL,0,&atx[42]} ,
  {0, "rot-12" ,128,2,0,0,0,0,0,0,NULL,&atx[3],NULL,0,&atx[43]} ,
  {0, "rot-13" ,128,3,0,0,0,0,0,0,NULL,&atx[3],NULL,0,&atx[44]} ,
  {0, "rot-21" ,128,4,0,0,0,0,0,0,NULL,&atx[3],NULL,0,&atx[45]} ,
  {0, "rot-22" ,128,5,0,0,0,0,0,0,NULL,&atx[3],NULL,0,&atx[46]} ,
  {0, "rot-23" ,128,6,0,0,0,0,0,0,NULL,&atx[3],NULL,0,&atx[47]} ,
  {0, "rot-31" ,128,7,0,0,0,0,0,0,NULL,&atx[3],NULL,0,&atx[48]} ,
  {0, "rot-32" ,128,8,0,0,0,0,0,0,NULL,&atx[3],NULL,0,&atx[49]} ,
  {0, "rot-33" ,128,9,0,0,0,0,0,0,NULL,&atx[3],NULL,0,NULL} ,
  {0, "translate" ,128,1,0,0,0,0,0,0,NULL,&atx[51],NULL,0,NULL} ,
  {440, "Trans-matrix" ,1,0,0,0,0,0,0,0,NULL,&atx[29],&atx[52],0,&atx[79]} ,
  {0, "scale-factor" ,128,0,0,0,0,0,0,0,NULL,&atx[3],NULL,0,&atx[53]} ,
  {0, "tran-1" ,128,1,0,0,0,0,0,0,NULL,&atx[3],NULL,0,&atx[54]} ,
  {0, "tran-2" ,128,2,0,0,0,0,0,0,NULL,&atx[3],NULL,0,&atx[55]} ,
  {0, "tran-3" ,128,3,0,0,0,0,0,0,NULL,&atx[3],NULL,0,NULL} ,
  {0, "camera" ,128,3,0,0,0,0,0,0,NULL,&atx[57],NULL,0,&atx[76]} ,
  {424, "Camera" ,1,0,0,0,0,0,0,0,NULL,&atx[29],&atx[58],0,&atx[77]} ,
  {0, "mode" ,128,0,0,0,0,0,0,0,NULL,&atx[3],&avnx[45],0,&atx[59]} ,
  {0, "x" ,128,1,0,0,0,0,0,0,NULL,&atx[60],NULL,0,&atx[65]} ,
  {436, "Model-space-point" ,1,0,0,0,0,0,0,0,NULL,&atx[29],&atx[61],0,&atx[72]} ,
  {0, "scale-factor" ,128,0,0,0,0,0,0,0,NULL,&atx[3],NULL,0,&atx[62]} ,
  {0, "x" ,128,1,0,0,0,0,0,0,NULL,&atx[3],NULL,0,&atx[63]} ,
  {0, "y" ,128,2,0,0,0,0,0,0,NULL,&atx[3],NULL,0,&atx[64]} ,
  {0, "z" ,128,3,0,0,0,0,0,0,NULL,&atx[3],NULL,0,NULL} ,
  {0, "y" ,128,2,0,0,0,0,0,0,NULL,&atx[60],NULL,0,&atx[66]} ,
  {0, "z" ,128,3,0,0,0,0,0,0,NULL,&atx[60],NULL,0,&atx[67]} ,
  {0, "up" ,128,4,0,0,0,0,0,0,NULL,&atx[60],NULL,0,&atx[68]} ,
  {0, "fore" ,128,5,0,1,0,0,0,0,NULL,&atx[60],NULL,0,&atx[69]} ,
  {0, "norm" ,128,6,0,1,0,0,0,0,NULL,&atx[60],NULL,0,&atx[70]} ,
  {0, "center" ,128,7,0,1,0,0,0,0,NULL,&atx[60],NULL,0,&atx[71]} ,
  {0, "tooclose" ,128,8,0,1,0,0,0,0,NULL,&atx[72],NULL,0,&atx[75]} ,
  {437, "RealValue" ,1,0,0,0,0,0,0,0,NULL,&atx[29],&atx[73],0,&atx[37]} ,
  {0, "scale-factor" ,128,0,0,0,0,0,0,0,NULL,&atx[3],NULL,0,&atx[74]} ,
  {0, "scaled-integer-value" ,128,1,0,0,0,0,0,0,NULL,&atx[3],NULL,0,NULL} ,
  {0, "toofar" ,128,9,0,1,0,0,0,0,NULL,&atx[72],NULL,0,NULL} ,
  {0, "script" ,128,4,0,0,0,0,0,0,NULL,&atx[77],NULL,0,&atx[96]} ,
  {425, "Biostruc-script" ,1,0,0,0,0,0,0,0,NULL,&atx[14],&atx[78],0,&atx[134]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[79],NULL,0,NULL} ,
  {441, "Biostruc-script-step" ,1,0,0,0,0,0,0,0,NULL,&atx[29],&atx[80],0,&atx[81]} ,
  {0, "step-id" ,128,0,0,0,0,0,0,0,NULL,&atx[81],NULL,0,&atx[82]} ,
  {442, "Step-id" ,1,0,0,0,0,0,0,0,NULL,&atx[3],NULL,0,NULL} ,
  {0, "step-name" ,128,1,0,1,0,0,0,0,NULL,&atx[8],NULL,0,&atx[83]} ,
  {0, "feature-do" ,128,2,0,1,0,0,0,0,NULL,&atx[14],&atx[84],0,&atx[90]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[85],NULL,0,NULL} ,
  {428, "Other-feature" ,1,0,0,0,0,0,0,0,NULL,&atx[29],&atx[86],0,&atx[114]} ,
  {0, "biostruc-id" ,128,0,0,0,0,0,0,0,NULL,&atx[87],NULL,0,&atx[88]} ,
  {412, "Biostruc-id" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[106]} ,
  {0, "set" ,128,1,0,0,0,0,0,0,NULL,&atx[2],NULL,0,&atx[89]} ,
  {0, "feature" ,128,2,0,0,0,0,0,0,NULL,&atx[19],NULL,0,NULL} ,
  {0, "camera-move" ,128,3,0,1,0,0,0,0,NULL,&atx[33],NULL,0,&atx[91]} ,
  {0, "pause" ,128,4,0,0,1,0,0,0,&avnx[50],&atx[3],NULL,0,&atx[92]} ,
  {0, "waitevent" ,128,5,0,0,0,0,0,0,NULL,&atx[93],NULL,0,&atx[94]} ,
  {301, "BOOLEAN" ,0,1,0,0,0,0,0,0,NULL,NULL,NULL,0,NULL} ,
  {0, "extra" ,128,6,0,0,0,0,0,0,NULL,&atx[3],NULL,0,&atx[95]} ,
  {0, "jump" ,128,7,0,1,0,0,0,0,NULL,&atx[81],NULL,0,NULL} ,
  {0, "user" ,128,5,0,0,0,0,0,0,NULL,&atx[97],NULL,0,NULL} ,
  {418, "User-object" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[12]} ,
  {0, "location" ,128,4,0,0,0,0,0,0,NULL,&atx[13],&atx[99],0,NULL} ,
  {0, "subgraph" ,128,0,0,0,0,0,0,0,NULL,&atx[100],NULL,0,&atx[133]} ,
  {402, "Chem-graph-pntrs" ,1,0,0,0,0,1,0,0,NULL,&atx[13],&atx[101],0,&atx[102]} ,
  {0, "atoms" ,128,0,0,0,0,0,0,0,NULL,&atx[102],NULL,0,&atx[113]} ,
  {403, "Atom-pntrs" ,1,0,0,0,0,1,0,0,NULL,&atx[29],&atx[103],0,&atx[174]} ,
  {0, "number-of-ptrs" ,128,0,0,0,0,0,0,0,NULL,&atx[3],NULL,0,&atx[104]} ,
  {0, "molecule-ids" ,128,1,0,0,0,0,0,0,NULL,&atx[14],&atx[105],0,&atx[107]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[106],NULL,0,NULL} ,
  {413, "Molecule-id" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[109]} ,
  {0, "residue-ids" ,128,2,0,0,0,0,0,0,NULL,&atx[14],&atx[108],0,&atx[110]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[109],NULL,0,NULL} ,
  {414, "Residue-id" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[112]} ,
  {0, "atom-ids" ,128,3,0,0,0,0,0,0,NULL,&atx[14],&atx[111],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[112],NULL,0,NULL} ,
  {415, "Atom-id" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[136]} ,
  {0, "residues" ,128,1,0,0,0,0,0,0,NULL,&atx[114],NULL,0,&atx[128]} ,
  {429, "Residue-pntrs" ,1,0,0,0,0,0,0,0,NULL,&atx[13],&atx[115],0,&atx[129]} ,
  {0, "explicit" ,128,0,0,0,0,0,0,0,NULL,&atx[116],NULL,0,&atx[122]} ,
  {431, "Residue-explicit-pntrs" ,1,0,0,0,0,0,0,0,NULL,&atx[29],&atx[117],0,&atx[124]} ,
  {0, "number-of-ptrs" ,128,0,0,0,0,0,0,0,NULL,&atx[3],NULL,0,&atx[118]} ,
  {0, "molecule-ids" ,128,1,0,0,0,0,0,0,NULL,&atx[14],&atx[119],0,&atx[120]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[106],NULL,0,NULL} ,
  {0, "residue-ids" ,128,2,0,0,0,0,0,0,NULL,&atx[14],&atx[121],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[109],NULL,0,NULL} ,
  {0, "interval" ,128,1,0,0,0,0,0,0,NULL,&atx[14],&atx[123],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[124],NULL,0,NULL} ,
  {432, "Residue-interval-pntr" ,1,0,0,0,0,0,0,0,NULL,&atx[29],&atx[125],0,&atx[140]} ,
  {0, "molecule-id" ,128,0,0,0,0,0,0,0,NULL,&atx[106],NULL,0,&atx[126]} ,
  {0, "from" ,128,1,0,0,0,0,0,0,NULL,&atx[109],NULL,0,&atx[127]} ,
  {0, "to" ,128,2,0,0,0,0,0,0,NULL,&atx[109],NULL,0,NULL} ,
  {0, "molecules" ,128,2,0,0,0,0,0,0,NULL,&atx[129],NULL,0,NULL} ,
  {430, "Molecule-pntrs" ,1,0,0,0,0,0,0,0,NULL,&atx[29],&atx[130],0,&atx[116]} ,
  {0, "number-of-ptrs" ,128,0,0,0,0,0,0,0,NULL,&atx[3],NULL,0,&atx[131]} ,
  {0, "molecule-ids" ,128,1,0,0,0,0,0,0,NULL,&atx[14],&atx[132],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[106],NULL,0,NULL} ,
  {0, "region" ,128,1,0,0,0,0,0,0,NULL,&atx[134],NULL,0,&atx[173]} ,
  {426, "Region-pntrs" ,1,0,0,0,0,0,0,0,NULL,&atx[29],&atx[135],0,&atx[197]} ,
  {0, "model-id" ,128,0,0,0,0,0,0,0,NULL,&atx[136],NULL,0,&atx[137]} ,
  {416, "Model-id" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[142]} ,
  {0, "region" ,128,1,0,0,0,0,0,0,NULL,&atx[13],&atx[138],0,NULL} ,
  {0, "site" ,128,0,0,0,0,0,0,0,NULL,&atx[14],&atx[139],0,&atx[146]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[140],NULL,0,NULL} ,
  {433, "Region-coordinates" ,1,0,0,0,0,0,0,0,NULL,&atx[29],&atx[141],0,&atx[148]} ,
  {0, "model-coord-set-id" ,128,0,0,0,0,0,0,0,NULL,&atx[142],NULL,0,&atx[143]} ,
  {417, "Model-coordinate-set-id" ,1,0,0,0,0,0,1,0,NULL,NULL,NULL,0,&atx[97]} ,
  {0, "number-of-coords" ,128,1,0,1,0,0,0,0,NULL,&atx[3],NULL,0,&atx[144]} ,
  {0, "coordinate-indices" ,128,2,0,1,0,0,0,0,NULL,&atx[14],&atx[145],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[3],NULL,0,NULL} ,
  {0, "boundary" ,128,1,0,0,0,0,0,0,NULL,&atx[14],&atx[147],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[148],NULL,0,NULL} ,
  {434, "Region-boundary" ,1,0,0,0,0,0,0,0,NULL,&atx[13],&atx[149],0,&atx[186]} ,
  {0, "sphere" ,128,0,0,0,0,0,0,0,NULL,&atx[150],NULL,0,&atx[153]} ,
  {405, "Sphere" ,1,0,0,0,0,1,0,0,NULL,&atx[29],&atx[151],0,&atx[154]} ,
  {0, "center" ,128,0,0,0,0,0,0,0,NULL,&atx[60],NULL,0,&atx[152]} ,
  {0, "radius" ,128,1,0,0,0,0,0,0,NULL,&atx[72],NULL,0,NULL} ,
  {0, "cone" ,128,1,0,0,0,0,0,0,NULL,&atx[154],NULL,0,&atx[158]} ,
  {406, "Cone" ,1,0,0,0,0,1,0,0,NULL,&atx[29],&atx[155],0,&atx[159]} ,
  {0, "axis-top" ,128,0,0,0,0,0,0,0,NULL,&atx[60],NULL,0,&atx[156]} ,
  {0, "axis-bottom" ,128,1,0,0,0,0,0,0,NULL,&atx[60],NULL,0,&atx[157]} ,
  {0, "radius-bottom" ,128,2,0,0,0,0,0,0,NULL,&atx[72],NULL,0,NULL} ,
  {0, "cylinder" ,128,2,0,0,0,0,0,0,NULL,&atx[159],NULL,0,&atx[163]} ,
  {407, "Cylinder" ,1,0,0,0,0,1,0,0,NULL,&atx[29],&atx[160],0,&atx[164]} ,
  {0, "axis-top" ,128,0,0,0,0,0,0,0,NULL,&atx[60],NULL,0,&atx[161]} ,
  {0, "axis-bottom" ,128,1,0,0,0,0,0,0,NULL,&atx[60],NULL,0,&atx[162]} ,
  {0, "radius" ,128,2,0,0,0,0,0,0,NULL,&atx[72],NULL,0,NULL} ,
  {0, "brick" ,128,3,0,0,0,0,0,0,NULL,&atx[164],NULL,0,NULL} ,
  {408, "Brick" ,1,0,0,0,0,1,0,0,NULL,&atx[29],&atx[165],0,&atx[33]} ,
  {0, "corner-000" ,128,0,0,0,0,0,0,0,NULL,&atx[60],NULL,0,&atx[166]} ,
  {0, "corner-001" ,128,1,0,0,0,0,0,0,NULL,&atx[60],NULL,0,&atx[167]} ,
  {0, "corner-010" ,128,2,0,0,0,0,0,0,NULL,&atx[60],NULL,0,&atx[168]} ,
  {0, "corner-011" ,128,3,0,0,0,0,0,0,NULL,&atx[60],NULL,0,&atx[169]} ,
  {0, "corner-100" ,128,4,0,0,0,0,0,0,NULL,&atx[60],NULL,0,&atx[170]} ,
  {0, "corner-101" ,128,5,0,0,0,0,0,0,NULL,&atx[60],NULL,0,&atx[171]} ,
  {0, "corner-110" ,128,6,0,0,0,0,0,0,NULL,&atx[60],NULL,0,&atx[172]} ,
  {0, "corner-111" ,128,7,0,0,0,0,0,0,NULL,&atx[60],NULL,0,NULL} ,
  {0, "alignment" ,128,2,0,0,0,0,0,0,NULL,&atx[174],NULL,0,&atx[196]} ,
  {404, "Chem-graph-alignment" ,1,0,0,0,0,1,0,0,NULL,&atx[29],&atx[175],0,&atx[150]} ,
  {0, "dimension" ,128,0,0,0,1,0,0,0,&avnx[51],&atx[3],NULL,0,&atx[176]} ,
  {0, "biostruc-ids" ,128,1,0,0,0,0,0,0,NULL,&atx[14],&atx[177],0,&atx[178]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[87],NULL,0,NULL} ,
  {0, "alignment" ,128,2,0,0,0,0,0,0,NULL,&atx[14],&atx[179],0,&atx[180]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[100],NULL,0,NULL} ,
  {0, "domain" ,128,3,0,1,0,0,0,0,NULL,&atx[14],&atx[181],0,&atx[182]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[100],NULL,0,NULL} ,
  {0, "transform" ,128,4,0,1,0,0,0,0,NULL,&atx[14],&atx[183],0,&atx[184]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[33],NULL,0,NULL} ,
  {0, "aligndata" ,128,5,0,1,0,0,0,0,NULL,&atx[14],&atx[185],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[186],NULL,0,NULL} ,
  {435, "Align-stats" ,1,0,0,0,0,0,0,0,NULL,&atx[29],&atx[187],0,&atx[60]} ,
  {0, "descr" ,128,0,0,1,0,0,0,0,NULL,&atx[8],NULL,0,&atx[188]} ,
  {0, "scale-factor" ,128,1,0,1,0,0,0,0,NULL,&atx[3],NULL,0,&atx[189]} ,
  {0, "vast-score" ,128,2,0,1,0,0,0,0,NULL,&atx[3],NULL,0,&atx[190]} ,
  {0, "vast-mlogp" ,128,3,0,1,0,0,0,0,NULL,&atx[3],NULL,0,&atx[191]} ,
  {0, "align-res" ,128,4,0,1,0,0,0,0,NULL,&atx[3],NULL,0,&atx[192]} ,
  {0, "rmsd" ,128,5,0,1,0,0,0,0,NULL,&atx[3],NULL,0,&atx[193]} ,
  {0, "blast-score" ,128,6,0,1,0,0,0,0,NULL,&atx[3],NULL,0,&atx[194]} ,
  {0, "blast-mlogp" ,128,7,0,1,0,0,0,0,NULL,&atx[3],NULL,0,&atx[195]} ,
  {0, "other-score" ,128,8,0,1,0,0,0,0,NULL,&atx[3],NULL,0,NULL} ,
  {0, "similarity" ,128,3,0,0,0,0,0,0,NULL,&atx[197],NULL,0,&atx[205]} ,
  {427, "Region-similarity" ,1,0,0,0,0,0,0,0,NULL,&atx[29],&atx[198],0,&atx[85]} ,
  {0, "dimension" ,128,0,0,0,1,0,0,0,&avnx[52],&atx[3],NULL,0,&atx[199]} ,
  {0, "biostruc-ids" ,128,1,0,0,0,0,0,0,NULL,&atx[14],&atx[200],0,&atx[201]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[87],NULL,0,NULL} ,
  {0, "similarity" ,128,2,0,0,0,0,0,0,NULL,&atx[14],&atx[202],0,&atx[203]} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[134],NULL,0,NULL} ,
  {0, "transform" ,128,3,0,0,0,0,0,0,NULL,&atx[14],&atx[204],0,NULL} ,
  {0, NULL,1,-1,0,0,0,0,0,0,NULL,&atx[33],NULL,0,NULL} ,
  {0, "indirect" ,128,4,0,0,0,0,0,0,NULL,&atx[85],NULL,0,NULL} };

static AsnModule ampx[1] = {
  { "MMDB-Features" , "mmdb3.h61",&atx[0],NULL,NULL,0,0} };

static AsnValxNodePtr avn = avnx;
static AsnTypePtr at = atx;
static AsnModulePtr amp = ampx;



/**************************************************
*
*    Defines for Module MMDB-Features
*
**************************************************/

#define BIOSTRUC_FEATURE_SET &at[0]
#define BIOSTRUC_FEATURE_SET_id &at[1]
#define BIOSTRUC_FEATURE_SET_descr &at[4]
#define BIOSTRUC_FEATURE_SET_descr_E &at[5]
#define BIOSTRUC_FEATURE_SET_features &at[15]
#define BIOSTRUC_FEATURE_SET_features_E &at[16]

#define CHEM_GRAPH_PNTRS &at[100]
#define CHEM_GRAPH_PNTRS_atoms &at[101]
#define CHEM_GRAPH_PNTRS_residues &at[113]
#define CHEM_GRAPH_PNTRS_molecules &at[128]

#define ATOM_PNTRS &at[102]
#define ATOM_PNTRS_number_of_ptrs &at[103]
#define ATOM_PNTRS_molecule_ids &at[104]
#define ATOM_PNTRS_molecule_ids_E &at[105]
#define ATOM_PNTRS_residue_ids &at[107]
#define ATOM_PNTRS_residue_ids_E &at[108]
#define ATOM_PNTRS_atom_ids &at[110]
#define ATOM_PNTRS_atom_ids_E &at[111]

#define CHEM_GRAPH_ALIGNMENT &at[174]
#define CHEM_GRAPH_ALIGNMENT_dimension &at[175]
#define CHEM_GRAPH_ALIGNMENT_biostruc_ids &at[176]
#define CHEM_GRAPH_ALIGNMENT_biostruc_ids_E &at[177]
#define CHEM_GRAPH_ALIGNMENT_alignment &at[178]
#define CHEM_GRAPH_ALIGNMENT_alignment_E &at[179]
#define CHEM_GRAPH_ALIGNMENT_domain &at[180]
#define CHEM_GRAPH_ALIGNMENT_domain_E &at[181]
#define CHEM_GRAPH_ALIGNMENT_transform &at[182]
#define CHEM_GRAPH_ALIGNMENT_transform_E &at[183]
#define CHEM_GRAPH_ALIGNMENT_aligndata &at[184]
#define CHEM_GRAPH_ALIGNMENT_aligndata_E &at[185]

#define SPHERE &at[150]
#define SPHERE_center &at[151]
#define SPHERE_radius &at[152]

#define CONE &at[154]
#define CONE_axis_top &at[155]
#define CONE_axis_bottom &at[156]
#define CONE_radius_bottom &at[157]

#define CYLINDER &at[159]
#define CYLINDER_axis_top &at[160]
#define CYLINDER_axis_bottom &at[161]
#define CYLINDER_radius &at[162]

#define BRICK &at[164]
#define BRICK_corner_000 &at[165]
#define BRICK_corner_001 &at[166]
#define BRICK_corner_010 &at[167]
#define BRICK_corner_011 &at[168]
#define BRICK_corner_100 &at[169]
#define BRICK_corner_101 &at[170]
#define BRICK_corner_110 &at[171]
#define BRICK_corner_111 &at[172]

#define TRANSFORM &at[33]
#define TRANSFORM_id &at[34]
#define TRANSFORM_moves &at[35]
#define TRANSFORM_moves_E &at[36]

#define BIOSTRUC_FEATURE_SET_ID &at[2]

#define BIOSTRUC_FEATURE_ID &at[19]

#define BIOSTRUC_FEATURE_SET_DESCR &at[6]
#define BIOSTRUC_FEATURE_SET_DESCR_name &at[7]
#define BIOSTRUC_FEATURE_SET_DESCR_pdb_comment &at[9]
#define BIOSTRUC_FEATURE_SET_DESCR_other_comment &at[10]
#define BIOSTRUC_FEATURE_SET_DESCR_attribution &at[11]

#define BIOSTRUC_FEATURE &at[17]
#define BIOSTRUC_FEATURE_id &at[18]
#define BIOSTRUC_FEATURE_name &at[20]
#define BIOSTRUC_FEATURE_type &at[21]
#define BIOSTRUC_FEATURE_property &at[22]
#define BIOSTRUC_FEATURE_property_color &at[23]
#define BIOSTRUC_FEATURE_property_render &at[30]
#define BIOSTRUC_FEATURE_property_transform &at[32]
#define BIOSTRUC_FEATURE_property_camera &at[56]
#define BIOSTRUC_FEATURE_property_script &at[76]
#define BIOSTRUC_FEATURE_property_user &at[96]
#define BIOSTRUC_FEATURE_location &at[98]
#define BIOSTRUC_FEATURE_location_subgraph &at[99]
#define BIOSTRUC_FEATURE_location_region &at[133]
#define BIOSTRUC_FEATURE_location_alignment &at[173]
#define BIOSTRUC_FEATURE_location_similarity &at[196]
#define BIOSTRUC_FEATURE_location_indirect &at[205]

#define COLOR_PROP &at[24]
#define COLOR_PROP_r &at[25]
#define COLOR_PROP_g &at[26]
#define COLOR_PROP_b &at[27]
#define COLOR_PROP_name &at[28]

#define RENDER_PROP &at[31]

#define CAMERA &at[57]
#define CAMERA_mode &at[58]
#define CAMERA_x &at[59]
#define CAMERA_y &at[65]
#define CAMERA_z &at[66]
#define CAMERA_up &at[67]
#define CAMERA_fore &at[68]
#define CAMERA_norm &at[69]
#define CAMERA_center &at[70]
#define CAMERA_tooclose &at[71]
#define CAMERA_toofar &at[75]

#define BIOSTRUC_SCRIPT &at[77]
#define BIOSTRUC_SCRIPT_E &at[78]

#define REGION_PNTRS &at[134]
#define REGION_PNTRS_model_id &at[135]
#define REGION_PNTRS_region &at[137]
#define REGION_PNTRS_region_site &at[138]
#define REGION_PNTRS_region_site_E &at[139]
#define REGION_PNTRS_region_boundary &at[146]
#define REGION_PNTRS_region_boundary_E &at[147]

#define REGION_SIMILARITY &at[197]
#define REGION_SIMILARITY_dimension &at[198]
#define REGION_SIMILARITY_biostruc_ids &at[199]
#define REGION_SIMILARITY_biostruc_ids_E &at[200]
#define REGION_SIMILARITY_similarity &at[201]
#define REGION_SIMILARITY_similarity_E &at[202]
#define REGION_SIMILARITY_transform &at[203]
#define REGION_SIMILARITY_transform_E &at[204]

#define OTHER_FEATURE &at[85]
#define OTHER_FEATURE_biostruc_id &at[86]
#define OTHER_FEATURE_set &at[88]
#define OTHER_FEATURE_feature &at[89]

#define RESIDUE_PNTRS &at[114]
#define RESIDUE_PNTRS_explicit &at[115]
#define RESIDUE_PNTRS_interval &at[122]
#define RESIDUE_PNTRS_interval_E &at[123]

#define MOLECULE_PNTRS &at[129]
#define MOLECULE_PNTRS_number_of_ptrs &at[130]
#define MOLECULE_PNTRS_molecule_ids &at[131]
#define MOLECULE_PNTRS_molecule_ids_E &at[132]

#define RESIDUE_EXPLICIT_PNTRS &at[116]
#define RESIDUE_EXPLICIT_PNTRS_number_of_ptrs &at[117]
#define RESIDUE_EXPLICIT_PNTRS_molecule_ids &at[118]
#define RESIDUE_EXPLICIT_PNTRS_molecule_ids_E &at[119]
#define RESIDUE_EXPLICIT_PNTRS_residue_ids &at[120]
#define RESIDUE_EXPLICIT_PNTRS_residue_ids_E &at[121]

#define RESIDUE_INTERVAL_PNTR &at[124]
#define RESIDUE_INTERVAL_PNTR_molecule_id &at[125]
#define RESIDUE_INTERVAL_PNTR_from &at[126]
#define RESIDUE_INTERVAL_PNTR_to &at[127]

#define REGION_COORDINATES &at[140]
#define REGION_COORDINATES_model_coord_set_id &at[141]
#define REGION_COORDINATES_number_of_coords &at[143]
#define REGION_COORDINATES_coordinate_indices &at[144]
#define REGION_COORDINATES_coordinate_indices_E &at[145]

#define REGION_BOUNDARY &at[148]
#define REGION_BOUNDARY_sphere &at[149]
#define REGION_BOUNDARY_cone &at[153]
#define REGION_BOUNDARY_cylinder &at[158]
#define REGION_BOUNDARY_brick &at[163]

#define ALIGN_STATS &at[186]
#define ALIGN_STATS_descr &at[187]
#define ALIGN_STATS_scale_factor &at[188]
#define ALIGN_STATS_vast_score &at[189]
#define ALIGN_STATS_vast_mlogp &at[190]
#define ALIGN_STATS_align_res &at[191]
#define ALIGN_STATS_rmsd &at[192]
#define ALIGN_STATS_blast_score &at[193]
#define ALIGN_STATS_blast_mlogp &at[194]
#define ALIGN_STATS_other_score &at[195]

#define MODEL_SPACE_POINT &at[60]
#define MODEL_SPACE_POINT_scale_factor &at[61]
#define MODEL_SPACE_POINT_x &at[62]
#define MODEL_SPACE_POINT_y &at[63]
#define MODEL_SPACE_POINT_z &at[64]

#define REALVALUE &at[72]
#define REALVALUE_scale_factor &at[73]
#define REALVALUE_scaled_integer_value &at[74]

#define MOVE &at[37]
#define MOVE_rotate &at[38]
#define MOVE_translate &at[50]

#define ROT_MATRIX &at[39]
#define ROT_MATRIX_scale_factor &at[40]
#define ROT_MATRIX_rot_11 &at[41]
#define ROT_MATRIX_rot_12 &at[42]
#define ROT_MATRIX_rot_13 &at[43]
#define ROT_MATRIX_rot_21 &at[44]
#define ROT_MATRIX_rot_22 &at[45]
#define ROT_MATRIX_rot_23 &at[46]
#define ROT_MATRIX_rot_31 &at[47]
#define ROT_MATRIX_rot_32 &at[48]
#define ROT_MATRIX_rot_33 &at[49]

#define TRANS_MATRIX &at[51]
#define TRANS_MATRIX_scale_factor &at[52]
#define TRANS_MATRIX_tran_1 &at[53]
#define TRANS_MATRIX_tran_2 &at[54]
#define TRANS_MATRIX_tran_3 &at[55]

#define BIOSTRUC_SCRIPT_STEP &at[79]
#define BIOSTRUC_SCRIPT_STEP_step_id &at[80]
#define BIOSTRUC_SCRIPT_STEP_step_name &at[82]
#define BIOSTRUC_SCRIPT_STEP_feature_do &at[83]
#define BIOSTRUC_SCRIPT_STEP_feature_do_E &at[84]
#define BIOSTRUC_SCRIPT_STEP_camera_move &at[90]
#define BIOSTRUC_SCRIPT_STEP_pause &at[91]
#define BIOSTRUC_SCRIPT_STEP_waitevent &at[92]
#define BIOSTRUC_SCRIPT_STEP_extra &at[94]
#define BIOSTRUC_SCRIPT_STEP_jump &at[95]

#define STEP_ID &at[81]
