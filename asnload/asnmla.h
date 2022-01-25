/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "asnmla.l60";
static AsnValxNodePtr avn;
static AsnTypePtr at;
static AsnModulePtr amp;


/**************************************************
*
*    Defines for Module NCBI-MedArchive
*
**************************************************/

#define MLA_REQUEST &at[7]
#define MLA_REQUEST_init &at[8]
#define MLA_REQUEST_getmle &at[10]
#define MLA_REQUEST_getpub &at[12]
#define MLA_REQUEST_gettitle &at[13]
#define MLA_REQUEST_citmatch &at[20]
#define MLA_REQUEST_fini &at[21]
#define MLA_REQUEST_getmriuids &at[22]
#define MLA_REQUEST_getaccuids &at[23]
#define MLA_REQUEST_uidtopmid &at[24]
#define MLA_REQUEST_pmidtouid &at[25]
#define MLA_REQUEST_getmlepmid &at[26]
#define MLA_REQUEST_getpubpmid &at[27]
#define MLA_REQUEST_citmatchpmid &at[28]
#define MLA_REQUEST_getmripmids &at[29]
#define MLA_REQUEST_getaccpmids &at[30]
#define MLA_REQUEST_citlstpmids &at[31]
#define MLA_REQUEST_getmleuid &at[32]
#define MLA_REQUEST_getmlrpmid &at[33]
#define MLA_REQUEST_getmlruid &at[34]

#define TITLE_MSG &at[14]
#define TITLE_MSG_type &at[15]
#define TITLE_MSG_title &at[18]

#define TITLE_TYPE &at[16]

#define TITLE_MSG_LIST &at[36]
#define TITLE_MSG_LIST_num &at[37]
#define TITLE_MSG_LIST_titles &at[38]
#define TITLE_MSG_LIST_titles_E &at[39]

#define ERROR_VAL &at[41]

#define MLA_BACK &at[42]
#define MLA_BACK_init &at[43]
#define MLA_BACK_error &at[44]
#define MLA_BACK_getmle &at[45]
#define MLA_BACK_getpub &at[46]
#define MLA_BACK_gettitle &at[47]
#define MLA_BACK_citmatch &at[48]
#define MLA_BACK_fini &at[49]
#define MLA_BACK_getuids &at[50]
#define MLA_BACK_getuids_E &at[51]
#define MLA_BACK_getpmids &at[52]
#define MLA_BACK_getpmids_E &at[53]
#define MLA_BACK_outuid &at[54]
#define MLA_BACK_outpmid &at[55]
#define MLA_BACK_getpme &at[56]
#define MLA_BACK_getmlr &at[57]
