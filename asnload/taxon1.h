/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "taxon1.l61";
static AsnValxNodePtr avn;
static AsnTypePtr at;
static AsnModulePtr amp;


/**************************************************
*
*    Defines for Module NCBI-Taxon1
*
**************************************************/

#define TAXON1_REQ &at[1]
#define TAXON1_REQ_init &at[2]
#define TAXON1_REQ_findname &at[4]
#define TAXON1_REQ_getdesignator &at[6]
#define TAXON1_REQ_getunique &at[7]
#define TAXON1_REQ_getidbyorg &at[8]
#define TAXON1_REQ_getorgnames &at[9]
#define TAXON1_REQ_getcde &at[11]
#define TAXON1_REQ_getranks &at[12]
#define TAXON1_REQ_getdivs &at[13]
#define TAXON1_REQ_getgcs &at[14]
#define TAXON1_REQ_getlineage &at[15]
#define TAXON1_REQ_getchildren &at[16]
#define TAXON1_REQ_getbyid &at[17]
#define TAXON1_REQ_lookup &at[18]
#define TAXON1_REQ_getorgmod &at[19]
#define TAXON1_REQ_fini &at[25]

#define TAXON1_INFO &at[20]
#define TAXON1_INFO_ival1 &at[21]
#define TAXON1_INFO_ival2 &at[22]
#define TAXON1_INFO_sval &at[23]

#define TAXON1_RESP &at[27]
#define TAXON1_RESP_error &at[28]
#define TAXON1_RESP_init &at[33]
#define TAXON1_RESP_findname &at[34]
#define TAXON1_RESP_findname_E &at[35]
#define TAXON1_RESP_getdesignator &at[42]
#define TAXON1_RESP_getunique &at[43]
#define TAXON1_RESP_getidbyorg &at[44]
#define TAXON1_RESP_getorgnames &at[45]
#define TAXON1_RESP_getorgnames_E &at[46]
#define TAXON1_RESP_getcde &at[47]
#define TAXON1_RESP_getcde_E &at[48]
#define TAXON1_RESP_getranks &at[49]
#define TAXON1_RESP_getranks_E &at[50]
#define TAXON1_RESP_getdivs &at[51]
#define TAXON1_RESP_getdivs_E &at[52]
#define TAXON1_RESP_getgcs &at[53]
#define TAXON1_RESP_getgcs_E &at[54]
#define TAXON1_RESP_getlineage &at[55]
#define TAXON1_RESP_getlineage_E &at[56]
#define TAXON1_RESP_getchildren &at[57]
#define TAXON1_RESP_getchildren_E &at[58]
#define TAXON1_RESP_getbyid &at[59]
#define TAXON1_RESP_lookup &at[66]
#define TAXON1_RESP_getorgmod &at[67]
#define TAXON1_RESP_getorgmod_E &at[68]
#define TAXON1_RESP_fini &at[69]

#define TAXON1_ERROR &at[29]
#define TAXON1_ERROR_level &at[30]
#define TAXON1_ERROR_msg &at[32]

#define TAXON1_NAME &at[36]
#define TAXON1_NAME_taxid &at[37]
#define TAXON1_NAME_cde &at[38]
#define TAXON1_NAME_oname &at[39]
#define TAXON1_NAME_uname &at[40]

#define TAXON1_DATA &at[60]
#define TAXON1_DATA_org &at[61]
#define TAXON1_DATA_div &at[62]
#define TAXON1_DATA_embl_code &at[63]
#define TAXON1_DATA_is_species_level &at[64]
