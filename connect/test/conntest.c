/*  $Id: conntest.c,v 6.3 1999/07/26 17:51:21 vakatov Exp $
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
* Author:  Denis Vakatov
*
* File Description:
*   Test suite for the NCBI connector(see in "connectn.[ch]", "connectr.h")
*
* --------------------------------------------------------------------------
* $Log: conntest.c,v $
* Revision 6.3  1999/07/26 17:51:21  vakatov
* A tiny cosmetic tweak
*
* Revision 6.2  1999/07/19 19:10:57  vakatov
* s_SingleBouncePrint():  printout BEFORE assert()
*
* Revision 6.1  1999/04/01 21:34:13  vakatov
* s_MultiBouncePrint(): decreased the # of ...SinglePrint() calls from 50 to 5
*
* Revision 6.0  1999/03/25 23:05:00  vakatov
* Initial revision
*
* ==========================================================================
*/

#ifndef _DEBUG
#define _DEBUG
#endif

#include <ncbi.h>
#include <conntest.h>


/***********************************************************************
 *  INTERNAL
 ***********************************************************************/


/* Standard error report
 */

#define s_ErrPost(status, descr)	\
 Nlm_ErrSetContext("",__FILE__,__LINE__,DBFLAG,0,0,0) ? 0 : \
 Nlm_ErrPostEx(status == eCONN_Success ? SEV_INFO : SEV_ERROR, 333, 444, \
 "%s (status: \"%s\")", descr, CONN_StatusString(status))


/* TESTs
 */

static void s_SingleBouncePrint
(CONN  conn,
 FILE* log_file)
{
  static const char write_str[] = "This is a s_*BouncePrint test string.\n";
  EConnStatus  status;
  Nlm_Uint4    n_written, n_read;
  Nlm_Char     buf[8192];

  s_ErrPost(eCONN_Success, "[s_SingleBouncePrint]  Starting...");

  /* WRITE */
  status = CONN_Write(conn, write_str, Nlm_StrLen(write_str), &n_written);
  if (status != eCONN_Success  ||  n_written != Nlm_StrLen(write_str))
    s_ErrPost(status, "[s_SingleBouncePrint] Write failed!");
  ASSERT( n_written == Nlm_StrLen(write_str) );
  ASSERT( status == eCONN_Success );


  /* READ the "bounced" data from the connection */
  status = CONN_Read(conn, buf, sizeof(buf)-1, &n_read, eCR_Persist);
  s_ErrPost(status, "[s_SingleBouncePrint] after READ");

  /* Printout to LOG file, if any */
  if (log_file  &&  n_read) {
    fprintf(log_file, "\ns_SingleBouncePrint(BEGIN PRINT)\n");
    ASSERT( Nlm_FileWrite(buf, n_read, 1, log_file) == 1 );
    fprintf(log_file, "\ns_SingleBouncePrint(END PRINT)\n");
  }

  /* Check-up */
  ASSERT( n_read >= n_written );
  buf[n_read] = '\0';
  ASSERT( Nlm_StrStr(buf, write_str) );
}


static void s_MultiBouncePrint
(CONN  conn,
 FILE* log_file)
{
  int i;

  s_ErrPost(eCONN_Success, "[s_MultiBouncePrint]  Starting...");
  if ( log_file )
    fprintf(log_file, "\ns_MultiBouncePrint(BEGIN)\n");
  for (i = 0;  i < 5;  i++)
    s_SingleBouncePrint(conn, log_file);
  s_ErrPost(eCONN_Success, "[s_MultiBouncePrint]  ...finished");
  if ( log_file )
    fprintf(log_file, "\ns_MultiBouncePrint(END)\n");
}


static void s_SingleBounceCheck
(CONN            conn,
 const STimeout* timeout,
 FILE*           log_file)
{
  static const char sym[] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};

#define TEST_N_LINES 200
#define TEST_BUF_SIZE (TEST_N_LINES * (TEST_N_LINES + 3) / 2)
  char  buf[TEST_BUF_SIZE];

  EConnStatus status;

  s_ErrPost(eCONN_Success, "[s_SingleBounceCheck]  Starting...");


  /* WRITE to the connection:  "0\n12\n345\n6789\n01234\n........"
   */

  {{
    size_t k = 0;
    Nlm_Uint4 i, j=0;
    for (i = 0;  k != sizeof(buf);  i++) {
      /* prepare output data */
      Nlm_Uint4 n_write, n_written;
      for (n_write = 0;  n_write < i;  n_write++, k++) {
        ASSERT( k != sizeof(buf) );
        buf[n_write] = sym[j++ % sizeof(sym)];
      }
      ASSERT( k != sizeof(buf) );
      if ( n_write ) {
        buf[n_write++] = '\n';  k++;
      }
      buf[n_write] = '\0';

      do { /* persistently */
        /* WAIT... sometimes */
        if (n_write % 5 == 3) {
          status = CONN_Wait(conn, eCONN_Write, timeout);
          if (status != eCONN_Success) {
            s_ErrPost(status, "The pre-WRITE CONN_Wait failed");
            ASSERT( status == eCONN_Timeout );
          }
        }

        /* WRITE */
        status = CONN_Write(conn, buf, n_write, &n_written);
        if (status != eCONN_Success) {
          s_ErrPost(status, "Write failed. Retrying...");
          ASSERT( n_written < n_write );
          ASSERT( status == eCONN_Timeout );
        }
        else {
          ASSERT( n_written == n_write );
        }
      } while (status != eCONN_Success);
    }
  }}

  /* READ the "bounced" data from the connection, the first TEST_BUF_SIZE
   * bytes must be:  "0\n12\n345\n6789\n01234\n........"
   */

  {{
    char* x_buf;
    Nlm_Uint4 n_read, n_to_read;

    Nlm_MemSet(buf, '\0', TEST_BUF_SIZE);

    /* PEEK until the 1st 1/3 of the "bounced" data is available */
    x_buf = buf;
    n_to_read = TEST_BUF_SIZE/3;

    do {
      status = CONN_Read(conn, x_buf, n_to_read, &n_read, eCR_Peek);
      if (status != eCONN_Success) {
        s_ErrPost(status, "The 1/3 PEEK failed. Retrying...");
        ASSERT( n_read < n_to_read );
        ASSERT( status == eCONN_Timeout );
      }
      if (n_read < n_to_read) {
        s_ErrPost(status, "Not all expected data is peeked yet. Continue...");
      }
    } while (n_read != n_to_read);

    /* READ 1st 1/3 of "bounced" data, compare it with the PEEKed data */
    status = CONN_Read(conn, x_buf + n_to_read, n_to_read, &n_read,
                       eCR_Read);
    ASSERT( status == eCONN_Success );
    ASSERT( n_read == n_to_read );
    ASSERT( Nlm_MemCmp(x_buf, x_buf + n_to_read, n_to_read) == 0 );
    Nlm_MemSet(x_buf + n_to_read, '\0', n_to_read);

    /* WAIT on read */
    status = CONN_Wait(conn, eCONN_Read, timeout);
    if (status != eCONN_Success) {
      s_ErrPost(status, "The 2/3 pre-READ CONN_Wait failed");
      ASSERT( status == eCONN_Timeout );
    }

    /* READ the 2nd 1/3 of "bounced" data */
    x_buf = buf + TEST_BUF_SIZE/3;
    n_to_read = TEST_BUF_SIZE/3;

    while ( n_to_read ) {
      s_ErrPost(status, "2/3 READ...");
      status = CONN_Read(conn, x_buf, n_to_read, &n_read, eCR_Read);
      if (status != eCONN_Success) {
        s_ErrPost(status, "The 2/3 READ failed. Retrying...");
        ASSERT( n_read < n_to_read );
        ASSERT( status == eCONN_Timeout );
      } else {
        ASSERT( n_read <= n_to_read );
      }
      n_to_read -= n_read;
      x_buf     += n_read;
    }
    ASSERT( status == eCONN_Success );

    /* Persistently READ the 3rd 1/3 of "bounced" data */
    n_to_read = TEST_BUF_SIZE - (x_buf - buf);

    status = CONN_Read(conn, x_buf, n_to_read, &n_read, eCR_Persist);
    if (status != eCONN_Success) {
      s_ErrPost(status, "The 3/3 (persistent) READ failed!");
      ASSERT( n_read < n_to_read );
      ASSERT( 0 );
    } else {
      ASSERT( n_read == n_to_read );
    }
  }}


  /* Check for the received "bounced" data is identical to the sent data
   */
  {{
    const char* x_buf = buf;
    Nlm_Uint4 i;
    size_t k=0, j=0;
    for (i = 1;  k != sizeof(buf);  i++) {
      Nlm_Uint4 n;
      for (n = 0;  n < i;  n++, k++) {
        if (k == sizeof(buf))
          break;
        ASSERT( *x_buf++ == sym[j++ % sizeof(sym)] );
      }
      ASSERT( *x_buf++ == '\n' );  k++;
    }
  }}

  /* Now when the "bounced" data is read and tested, READ an arbitrary extra
   * data sent by the peer and print it out to LOG file
   */

  if ( !log_file )
    return;

  fprintf(log_file, "\ns_SingleBounceCheck(BEGIN EXTRA DATA)\n");
  for (;;) {
    Nlm_Uint4 n_read;

    status = CONN_Read(conn, buf, sizeof(buf), &n_read, eCR_Persist);
    s_ErrPost(status, "s_SingleBounceCheck(The extra data READ...)");
    if ( n_read )
      ASSERT( Nlm_FileWrite(buf, n_read, 1, log_file) == 1 );
    if (status == eCONN_Closed  ||  status == eCONN_Timeout)
      break; /* okay */

    ASSERT( status == eCONN_Success  ||  status == eCONN_Timeout );
  }
  fprintf(log_file, "\ns_SingleBounceCheck(END EXTRA DATA)\n\n");
}



/***********************************************************************
 *  EXTERNAL
 ***********************************************************************/


NLM_EXTERN void Ncbi_TestConnector
(CONNECTOR       connector,
 const STimeout* timeout,
 FILE*           log_file,
 Nlm_Uint4       flags)
{
  EConnStatus status;
  CONN        conn;

  ErrSetLogLevel(SEV_INFO);
  ErrSetMessageLevel(SEV_INFO);
  ErrSetOptFlags(EO_SHOW_FILELINE | EO_SHOW_MSGTEXT);

  s_ErrPost(eCONN_Success, "[Ncbi_TestConnector]  Starting...");


  /* CREATE new connection on the base of the connector, set
   * TIMEOUTs, try to RECONNECT, WAIT for the connection is writeable
   */

  ASSERT( CONN_Create(connector, &conn) == eCONN_Success );

  CONN_SetTimeout(conn, eCONN_Connect,   timeout);
  CONN_SetTimeout(conn, eCONN_Read,      timeout);
  CONN_SetTimeout(conn, eCONN_ReadWrite, timeout);
  CONN_SetTimeout(conn, eCONN_Close,     timeout);

  ASSERT( CONN_Reconnect(conn, connector) == eCONN_Success );

  CONN_SetTimeout(conn, eCONN_Write, timeout);

  status = CONN_Wait(conn, eCONN_Write, timeout);
  if (status != eCONN_Success) {
    s_ErrPost(status, "First CONN_Wait failed");
    ASSERT( status == eCONN_Timeout );
  }


  /* Run the specified TESTs
   */

  if ( !flags )
    flags = 0xffffffff;
  if (flags & TESTCONN_SINGLE_BOUNCE_PRINT)
    s_SingleBouncePrint(conn, log_file);
  if (flags & TESTCONN_MULTI_BOUNCE_PRINT)
    s_MultiBouncePrint(conn, log_file);
  if (flags & TESTCONN_SINGLE_BOUNCE_CHECK)
    s_SingleBounceCheck(conn, timeout, log_file);


  /* And CLOSE the connection...
   */

  ASSERT( CONN_Close(conn) == eCONN_Success );
}
