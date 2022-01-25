/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "sugmap.l60";
static AsnValxNodePtr avn;
static AsnTypePtr at;
static AsnModulePtr amp;


/**************************************************
*
*    Defines for Module NCBI-Suggest
*
**************************************************/

#define SUGGEST_REQUEST &at[2]
#define SUGGEST_REQUEST_init &at[3]
#define SUGGEST_REQUEST_intervals &at[5]
#define SUGGEST_REQUEST_fini &at[19]

#define SUGGEST_INTERVALS &at[6]
#define SUGGEST_INTERVALS_params &at[7]
#define SUGGEST_INTERVALS_dna &at[16]
#define SUGGEST_INTERVALS_protein &at[17]
#define SUGGEST_INTERVALS_code &at[18]

#define SUGGEST_PARAMETERS &at[8]
#define SUGGEST_PARAMETERS_size &at[9]
#define SUGGEST_PARAMETERS_begin_search &at[11]
#define SUGGEST_PARAMETERS_end_search &at[12]
#define SUGGEST_PARAMETERS_term_stop &at[13]

#define SUGGEST_ERROR &at[21]
#define SUGGEST_ERROR_level &at[22]
#define SUGGEST_ERROR_msg &at[24]

#define SUGGEST_RESPONSE &at[26]
#define SUGGEST_RESPONSE_init &at[27]
#define SUGGEST_RESPONSE_error &at[28]
#define SUGGEST_RESPONSE_intervals &at[29]
#define SUGGEST_RESPONSE_fini &at[30]
