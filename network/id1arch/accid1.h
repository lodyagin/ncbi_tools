#ifndef __ACCID1_H_
#define __ACCID1_H_

#include <id1arch.h>

NLM_EXTERN Boolean LIBCALL ID1Init PROTO((void));
NLM_EXTERN void LIBCALL ID1Fini PROTO((void));

/**** Look Up a Uid from a SeqId using ID1 lookup ****/
NLM_EXTERN Int4 LIBCALL ID1FindSeqId PROTO((SeqIdPtr sip));

/**** Look Up the source SeqId given a GI ****************/
NLM_EXTERN SeqIdPtr LIBCALL ID1SeqIdForGI PROTO ((Int4 gi));

NLM_EXTERN SeqEntryPtr LIBCALL ID1SeqEntryGet PROTO((Int4 uid, Int2 retcode));


/*****************************************************************************
*
*   The Following two functions allow access by BioseqFetch using the
*   SeqMgr.  The application should call ID1BioseqFetchEnable() at the start
*   of the application and ID1BioseqFetchDisable() at the end; This
*   will make ID1BioseqFetch() the "remote" access procedure for the
*   SeqMgr. ID1Init() will only be called on the first fetch unless "now"
*   is true;
*
*   If you add your own fetch function after calling ID1BioseqFetchEnable,
*     it will be called BEFORE ID1BioseqFetchEnable. Add yours after this
*     call, and yours will be call AFTER ID1.
*
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL ID1BioseqFetchEnable PROTO((CharPtr progname, Boolean now));
NLM_EXTERN void LIBCALL ID1BioseqFetchDisable PROTO((void));

#endif
