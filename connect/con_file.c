/*  $Id: con_file.c,v 6.9 1999/06/29 16:17:42 aleksey Exp $
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
*  Implement CONNECTOR for a FILE stream (based on the NCBI "FILE").
*
* --------------------------------------------------------------------------
* $Log: con_file.c,v $
* Revision 6.9  1999/06/29 16:17:42  aleksey
* Updated test procedure, some complication was added
*
* ==========================================================================
* Revision 6.8  1999/06/29 15:18:15  aleksey
* Updated s_VT_Read and test procedures
* ==========================================================================
* Revision 6.7  1999/06/28 20:16:04  aleksey
* Extra comments were added since last realese
* ==========================================================================
* Revision 6.6  1999/06/28 18:12:35  aleksey
* Fixed final release
* ==========================================================================
* Revision 6.5  1999/06/28 16:04:10  aleksey
* Fixed compiler warnings
* Testing application has been provided
* ==========================================================================
* Revision 6.4  1999/06/28 14:29:50  kans
* more fixes to Mac compiler warnings
* ==========================================================================
* Revision 6.3  1999/06/28 14:14:30  aleksey
* Intermediate revision with test features
* ==========================================================================
* Revision 6.2  1999/06/25 18:24:36  kans
* fixes to Mac compiler warnings
* ==========================================================================
* Revision 6.1  1999/06/24 20:26:13  aleksey
* Initial revision
* ==========================================================================
*/

#include <ncbi.h>
#include <con_file.h>

typedef struct
{
  Nlm_Char*     in_file_name,
          *     out_file_name;
  FILE*         fin,
      *         fout;
  SFileConnAttr attr;
} SFileConnector;

/*
 *
 */
static const Nlm_Char* s_VT_GetType (void* connector) { return "FILE"; }

/*
 *
 */
static EConnStatus s_VT_Connect (void* connector, const STimeout* timeout)
{
SFileConnector* fcon = (SFileConnector*) connector;
Char*           mode;

  /* check if it is connected already,
   * and a user is trying to reconnect
   */
  if (fcon->fin && fcon->fout) return eCONN_Success;

  fcon->fin = Nlm_FileOpen (fcon->in_file_name, "rb");
  if (!fcon->fin) return eCONN_InvalidArg;

  switch (fcon->attr.w_mode)
  {
    case eFCM_Truncate : mode = "wb";  break;
    case eFCM_Seek     : /* mode = "rb+"; break; */
    case eFCM_Append   : mode = "ab";  break;
  }

  fcon->fout = Nlm_FileOpen (fcon->out_file_name, mode);

  if (!fcon->fout)
  {
    Nlm_FileClose (fcon->fin);
    fcon->fin = NULL;
    return eCONN_InvalidArg;
  }

  /* due to the shortage of portable 'fseek' call
   * ignore read/write positions so far,
   * only 0/EOF are in use for writing,
   * and only 0 for reading
   */

  return eCONN_Success;
}

/*
 *
 */
static EConnStatus s_VT_Write (void*           connector,
                               const void*     buf,
                               Nlm_Uint4       size,
                               Nlm_Uint4*      n_written,
                               const STimeout* timeout)
{
SFileConnector* fcon = (SFileConnector*) connector;

  ASSERT (fcon->fout);

  /* return buffer unchanged if user requested to write 0 bytes
   */
  if (!size)
  {
    *n_written = size;
    return eCONN_Success;
  }

  /* hard to determine errors using this call, otherwise
   * it might be not portable
   */
  *n_written = Nlm_FileWrite (buf, 1, size, fcon->fout);

  if (!*n_written) return eCONN_Unknown;

  /* return OK so far
   */
  return eCONN_Success;
}

/*
 *
 */
static EConnStatus s_VT_Read (void*           connector,
                              void*           buf,
                              Nlm_Uint4       size,
                              Nlm_Uint4*      n_read,
                              const STimeout* timeout)
{
SFileConnector* fcon = (SFileConnector*) connector;

  ASSERT (fcon->fin);

  /* return buffer unchanged if user requested to read 0 bytes
   */
  if (!size)
  {
    *n_read = size;
    return eCONN_Success;
  }

  /* hard to determine errors using this call, otherwise
   * it might be not portable
   */
  *n_read = Nlm_FileRead (buf, 1, size, fcon->fin);

  /* if nothing was read
   * check for errors
   */
  if (!*n_read) /* check if an error occured */
  {
    if (feof (fcon->fin)) /* pointer behind eof */
      return eCONN_Closed;
    else
      return eCONN_Unknown;
  }

  /* return OK so far
   */
  return eCONN_Success;
}

/* always return OK
 *
 */
static EConnStatus s_VT_Wait (void*           connector,
                              EConnDirection  direction,
                              const STimeout* timeout)
{
  return eCONN_Success;
}

/* always return OK
 *
 */
static EConnStatus s_VT_Flush (void* connector, const STimeout* timeout)
{
SFileConnector* fcon = (SFileConnector*) connector;

  ASSERT (fcon->fout);

  fflush (fcon->fout);

  return eCONN_Success;
}

/* always return OK
 *
 */
static EConnStatus s_VT_Close (CONNECTOR connector, const STimeout* timeout)
{
SFileConnector* fcon = (SFileConnector*) connector->handle;

  Nlm_FileClose (fcon->fin);
  Nlm_FileClose (fcon->fout);

  Nlm_MemFree (fcon->in_file_name);
  Nlm_MemFree (fcon->out_file_name);
  Nlm_MemFree (fcon);
  Nlm_MemFree (connector);

  return eCONN_Success;
}

/*
 *
 */
CONNECTOR FILE_CreateConnectorEx (const Nlm_Char*      in_file_name,
                                  const Nlm_Char*      out_file_name,
                                  const SFileConnAttr* attr)
{
CONNECTOR       con  = (SConnector    *) Nlm_MemNew (sizeof(SConnector    ));
SFileConnector* fcon = (SFileConnector*) Nlm_MemNew (sizeof(SFileConnector));

  ASSERT (con);
  ASSERT (fcon);

  /* initialize
   */
  fcon->in_file_name  = Nlm_StringSave (in_file_name);
  fcon->out_file_name = Nlm_StringSave (out_file_name);

  Nlm_MemCopy (&fcon->attr, attr, sizeof (SFileConnAttr));

  /* initialize handle
   */
  con->handle = fcon;

  /* initialize virtual table
   */ 
  con->vtable.get_type   = s_VT_GetType;
  con->vtable.connect    = s_VT_Connect;
  con->vtable.wait       = s_VT_Wait;
#ifdef IMPLEMENTED__CONN_WaitAsync
  con->vtable.wait_async = s_VT_Wait;
#endif
  con->vtable.write      = s_VT_Write;
  con->vtable.flush      = s_VT_Flush;
  con->vtable.read       = s_VT_Read;
  con->vtable.close      = s_VT_Close;

  return con;
}

/*
 *
 */
CONNECTOR FILE_CreateConnector (const Nlm_Char* in_file_name,
                                const Nlm_Char* out_file_name)
{
SFileConnAttr attr;

  Nlm_MemSet (&attr, 0, sizeof (SFileConnAttr));

  return FILE_CreateConnectorEx (in_file_name, out_file_name, &attr);
}

/********************************************************
 * testing
 ********************************************************/

#if TEST_MODULE_CONN_FILE

#define CONN_FILE_SIZE 256

#include "conntest.h"

extern
#ifdef __cplusplus
"C"
#endif
Int2 Main(void)
{
  STimeout  timeout = {0, 0};
  CONNECTOR connector;
  FILE*     log_file;

  Int4        argc = Nlm_GetArgc();
  CharPtr*    argv = Nlm_GetArgv();

  CharPtr     in_file_name  = "con_file.tst",
              out_file_name = "con_file.tst";

  Char        buffer [CONN_FILE_SIZE] = "This is the test\n"
                                        "of connector which is bound\n"
                                        "to an ordinary file\n";
  Nlm_Uint4   n_read;
  EConnStatus status;
  CONN        conn;
  
  Nlm_ErrSetOpts(ERR_TEE, ERR_LOG_ON);
  ErrSetLogLevel(SEV_INFO);
  ErrSetMessageLevel(SEV_INFO);
  log_file = FileOpen("con_file.log", "wb");

  ErrPostEx(SEV_INFO, 0, 0,
            "Starting the CON_FILE test...\n"
            "Input file name is %s\n"
            "Output file name is %s\n",
            in_file_name, out_file_name);

  {
    FILE* test = Nlm_FileOpen (in_file_name, "rb");

    if (!test)
    {
      test = Nlm_FileOpen (in_file_name, "wb");
      Nlm_FileWrite (buffer, 1, CONN_FILE_SIZE, test);
    }
    Nlm_FileClose (test);
  }

  /* the awful thing is that the file for reading is opened
   * prior to the file for writing and it MUST exist, otherwise
   * the library will core dump
   */
  connector = FILE_CreateConnector(in_file_name, out_file_name);
  ASSERT (connector);

  VERIFY (CONN_Create (connector, &conn)==eCONN_Success);

  {
    int i;
    for (i=0; i<10; ++i)
      VERIFY (CONN_Write (conn, buffer, CONN_FILE_SIZE, &n_read)==
              eCONN_Success);
  }

  Nlm_MemSet (buffer, 0, sizeof (buffer));

  {
    EConnStatus status;
    Char*       fmt;
    
    do
    {
      status = CONN_Read  (conn, buffer, CONN_FILE_SIZE, &n_read, eCR_Read);
      if (status == eCONN_Success)
        puts (buffer);
    } while (status == eCONN_Success);

    switch (status)
    {
      case eCONN_Closed    : fmt = "Status is [%d] Closed\n"; break;
      case eCONN_InvalidArg: fmt = "Status is [%d] Invalid Arg\n"; break;
      default              :
      case eCONN_Unknown   : fmt = "Status is [%d] Unknown\n"; break;
    }
    printf (fmt, (int)status);
  }

  VERIFY (CONN_Close (conn)==eCONN_Success);
  
  FileClose (log_file);

  return 0;
}

#endif /* TEST_MODULE_CONN_FILE */
