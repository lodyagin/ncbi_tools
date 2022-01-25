/***********************************************************************
*
**
*        Automatic header module from ASNTOOL
*
************************************************************************/

#ifndef _ASNTOOL_
#include <asn.h>
#endif

static char * asnfilename = "spell.l";
static AsnValxNodePtr avn;
static AsnTypePtr at;
static AsnModulePtr amp;


/**************************************************
*
*    Defines for Module NCBI-SPELL
*
**************************************************/

#define SPELL_REQUEST &at[0]
#define SPELL_REQUEST_init &at[1]
#define SPELL_REQUEST_spelltext &at[3]
#define SPELL_REQUEST_fini &at[5]

#define SPELL_RESPONSE &at[7]
#define SPELL_RESPONSE_error &at[8]
#define SPELL_RESPONSE_init &at[10]
#define SPELL_RESPONSE_spelltext &at[11]
#define SPELL_RESPONSE_spelltext_E &at[12]
#define SPELL_RESPONSE_fini &at[14]
