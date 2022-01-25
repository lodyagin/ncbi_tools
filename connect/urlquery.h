/*   urlquery.h
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*            National Center for Biotechnology Information (NCBI)
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government do not place any restriction on its use or reproduction.
*  We would, however, appreciate having the NCBI and the author cited in
*  any work or product based on this material
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
* ===========================================================================
*
* File Name:  urlquery.h
*
* Author:  Jonathan Kans
*
* Version Creation Date:   4/16/98
*
* $Revision: 6.1 $
*
* File Description: 
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
*
* ==========================================================================
*/

#ifndef _URLQUERY_
#define _URLQUERY_

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include <ncbi.h>
#include <con_url.h>


/* URL QUERY CONVENIENCE FUNCTIONS */

/* Content-type flag. */
typedef enum {
  wwwencoded = 1,
  multipart
}  UrlContentType;

/*
  Initializes a POST connector based on URL query parameters.  For example:

  QUERY_OpenUrlQuery (
    "cruncher.nlm.nih.gov", 80,
    "/cgi-bin/Sequin/testcgi.cgi",
    "request=seg&window=12&lowcut=2.3&hicut=2.6",
    "Sequin", 30, wwwencoded, URLC_SURE_FLUSH
  );

  The returned CONN value is then passed data before being sent to the cgi.
*/
NLM_EXTERN CONN QUERY_OpenUrlQuery (
  Nlm_CharPtr host_machine, Nlm_Uint2 host_port,
  Nlm_CharPtr host_path, Nlm_CharPtr arguments,
  Nlm_CharPtr appName, Nlm_Uint4 timeoutsec,
  UrlContentType ctype, URLC_Flags flags
);

/*
  Copies file to CONN by repeated calls to CONN_Write.  Writing of data must be
  done before the query is sent.  For the above seg example, a FASTA file of
  protein sequence data is appropriate.
*/
NLM_EXTERN void QUERY_CopyFileToQuery (
  CONN conn, FILE *fp
);

/*
  Calculates length of data written to connection, writes HTML header, and calls
  CONN_Wait with a timeout of 0 to send the query without waiting for response.
*/
NLM_EXTERN void QUERY_SendQuery (
  CONN conn
);

/*
  Copies results from CONN by repeated calls to CONN_Read.  Reading of data is
  only done when the connection is ready (CONN_Wait returns eCONN_Success).  In
  the query system below, the user's completion routine is called when results
  are ready to be read.  For the seg example, the results are in FASTA format
  with lower case x characters replacing amino acids in low-complexity regions.
*/
NLM_EXTERN void QUERY_CopyResultsToFile (
  CONN conn, FILE *fp
);


/* FUNCTIONS FOR MAINTAINING A QUEUE OF PENDING URL QUERIES */

/* Callback type for queued queries */
typedef Nlm_Boolean (LIBCALLBACK *QueryResultProc) (CONN conn, Nlm_VoidPtr userdata, EConnStatus status);

/* Opaque handle type.  Variable must be kept by application and initialized to NULL. */
struct SQueueTag;
typedef struct SQueueTag* QUEUE;  /* queue handle */

/*
  Records connection, completion routine, and user data in queue.
*/
NLM_EXTERN void QUERY_AddToQueue (
  QUEUE* queue, CONN conn, QueryResultProc resultproc, Nlm_VoidPtr userdata
);

/*
  Checks queued connections (with CONN_Wait), calls completion routine, then removes
  query from queue and closes connection.  Application is responsible for calling
  QUERY_CheckQueue every once in a while, typically with a timer.
*/
NLM_EXTERN Nlm_Int4 QUERY_CheckQueue (
  QUEUE* queue
);


#ifdef __cplusplus
}
#endif

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif /* _URLQUERY_ */

