#ifndef CON_FILE_H
#define CON_FILE_H
/*  $Id: con_file.h,v 6.8 1999/06/28 20:16:05 aleksey Exp $
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================
*
* Author: Vladimir Alekseyev 
*
* File Description:
*  Implement CONNECTOR for a file stream (based on the NCBI "FILE").
*
* --------------------------------------------------------------------------
* $Log: con_file.h,v $
* Revision 6.8  1999/06/28 20:16:05  aleksey
* Extra comments were added since last realese
*
* ==========================================================================
* Revision 6.7  1999/06/28 18:12:36  aleksey
* Fixed final release
* ==========================================================================
* Revision 6.6  1999/06/28 16:04:10  aleksey
* Fixed compiler warnings
* Testing application has been provided
* ==========================================================================
* Revision 6.5  1999/06/28 14:29:50  kans
* more fixes to Mac compiler warnings
* ==========================================================================
* Revision 6.4  1999/06/28 14:14:31  aleksey
* Intermediate revision with test features
* ==========================================================================
* Revision 6.3  1999/06/25 18:24:36  kans
* fixes to Mac compiler warnings
* ==========================================================================
* Revision 6.2  1999/06/24 20:28:31  aleksey
* Corrected initial revision
* ==========================================================================
*/

#include <connectr.h>
#include <connutil.h>


#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* An output connection to a file can overrite, write to
 * any position or to the end of file
 */
typedef enum
{
  eFCM_Truncate, /* create new or replace existing one */
  eFCM_Seek,     /* take start position in count */
  eFCM_Append    /* add to the end of file */
} EFileConnMode;

/* Access to file(s) connector attributes
 */
typedef struct
{
  Nlm_Int4      r_pos; /* begin to read from  r_pos in seek mode */
  Nlm_Int4      w_pos; /* begin to write from w_pos in seek mode */
  EFileConnMode w_mode;/* where to set the pointer / how to open output file */
} SFileConnAttr;

/* Create new CONNECTOR structure to handle connection to a file
 * with access from the begining.
 * Return: CONNECTOR if OK or NULL if error
 */
NLM_EXTERN CONNECTOR FILE_CreateConnector
(const Nlm_Char* in_file_name,
 const Nlm_Char* out_file_name
 );

/* Create new CONNECTOR structure to handle connection to a file
 * with access defined in passed attributes.
 * Return: CONNECTOR if OK or NULL if error
 */
NLM_EXTERN CONNECTOR FILE_CreateConnectorEx
(const Nlm_Char*      in_file_name,
 const Nlm_Char*      out_file_name,
 const SFileConnAttr* attr
 );

#ifdef __cplusplus
}
#endif /* CON_FILE_H */

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif
