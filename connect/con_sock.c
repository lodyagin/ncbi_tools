/*  $Id: con_sock.c,v 6.0 1999/04/01 22:09:31 vakatov Exp $
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
*  Implement CONNECTOR for a network socket(based on the NCBI "SOCK").
*
*  See in "connectr.h" for the detailed specification of the underlying
*  connector("CONNECTOR", "SConnectorTag") methods and structures.
*
* --------------------------------------------------------------------------
* $Log: con_sock.c,v $
* Revision 6.0  1999/04/01 22:09:31  vakatov
* Initial revision
*
* ==========================================================================
*/

#include <ncbi.h>
#include <ncbisock.h>
#include <con_sock.h>


/***********************************************************************
 *  INTERNAL -- Auxiliary types and static functions
 ***********************************************************************/

/* All internal data necessary to perform the (re)connect and i/o
 */
typedef struct {
  SOCK      sock;      /* socket;  NULL if not connected yet */
  Nlm_Char* host;      /* server:  host */
  Nlm_Uint2 port;      /* server:  service port */
  Nlm_Uint4 conn_try;  /* max.number of attempts to establish conn */
  void*     init_data; /* data to send to the server on connect */
  Nlm_Uint4 init_size; /* size of the "inst_str" buffer */
} SSockConnector;


/***********************************************************************
 *  INTERNAL -- "s_VT_*" functions for the "virt.table" of connector methods
 ***********************************************************************/

static const Nlm_Char* s_VT_GetType
(void* connector)
{
  return "SOCK";
}


static EConnStatus s_VT_Connect
(void*           connector,
 const STimeout* timeout)
{
  SSockConnector* xxx = (SSockConnector*)connector;
  ESOCK_ErrCode   err_code = eSOCK_ESuccess;

  Nlm_Uint4 i;
  for (i = 0;  i < xxx->conn_try;  i++) {
    /* connect */
    err_code = xxx->sock ?
      SOCK_Reconnect(xxx->sock, 0, 0, timeout) :
      SOCK_Create(xxx->host, xxx->port, timeout, &xxx->sock);

    if (err_code == eSOCK_ESuccess) {
      /* write init data, if any */
      Nlm_Uint4 n_written = 0;
      if ( !xxx->init_data )
        return eCONN_Success;
      err_code = SOCK_Write(xxx->sock, xxx->init_data, xxx->init_size,
                            &n_written);
      if (err_code == eSOCK_ESuccess)
        return eCONN_Success;
    }

    /* error: close socket and continue trying */
    if ( xxx->sock ) {
      SOCK_Close(xxx->sock);
      xxx->sock = 0;
    }
  }

  /* error: return status */
  return ESOCK2ECONN(err_code);
}


static EConnStatus s_VT_Wait
(void*           connector,
 EConnDirection  direction,
 const STimeout* timeout)
{
  SSockConnector* xxx = (SSockConnector*)connector;

  if (direction != eCONN_Read  &&  direction != eCONN_Write)
    return eCONN_InvalidArg;

  return ESOCK2ECONN( SOCK_Select(xxx->sock,
                                  direction == eCONN_Read ?
                                  eSOCK_OnRead : eSOCK_OnWrite,
                                  timeout) );
}


#ifdef IMPLEMENTED__CONN_WaitAsync
static EConnStatus s_VT_WaitAsync
(void*                   connector,
 FConnectorAsyncHandler  func,
 SConnectorAsyncHandler* data)
{
  return eCONN_NotSupported;
}
#endif


static EConnStatus s_VT_Write
(void*           connector,
 const void*     buf,
 Nlm_Uint4       size,
 Nlm_Uint4*      n_written,
 const STimeout* timeout)
{
  SSockConnector* xxx = (SSockConnector*)connector;

  SOCK_SetTimeout(xxx->sock, eSOCK_OnWrite, timeout, 0, 0);

  return ESOCK2ECONN( SOCK_Write(xxx->sock, buf, size, n_written) );
}


static EConnStatus s_VT_Flush
(void*           connector,
 const STimeout* timeout)
{
  return eCONN_Success;
}


static EConnStatus s_VT_Read
(void*           connector,
 void*           buf,
 Nlm_Uint4       size,
 Nlm_Uint4*      n_read,
 const STimeout* timeout)
{
  SSockConnector* xxx = (SSockConnector*)connector;

  SOCK_SetTimeout(xxx->sock, eSOCK_OnRead, timeout, 0, 0);

  return ESOCK2ECONN( SOCK_Read(xxx->sock, buf, size, n_read) );
}


static EConnStatus s_VT_Close
(CONNECTOR       connector,
 const STimeout* timeout)
{
  SSockConnector* xxx = (SSockConnector*)connector->handle;
  EConnStatus status = eCONN_Success;

  if ( xxx->sock ) {
    SOCK_SetTimeout(xxx->sock, eSOCK_OnWrite, timeout, 0, 0);
    status = ESOCK2ECONN(SOCK_Close(xxx->sock));
  }

  Nlm_MemFree(xxx->host);
  Nlm_MemFree(xxx->init_data);
  Nlm_MemFree(xxx);
  Nlm_MemFree(connector);
  return status;
}



/***********************************************************************
 *  EXTERNAL -- the connector's "constructor"
 ***********************************************************************/


NLM_EXTERN CONNECTOR SOCK_CreateConnector
(const Nlm_Char* host,
 Nlm_Uint2       port,
 Nlm_Uint4       conn_try)
{
  return SOCK_CreateConnectorEx(host, port, conn_try, 0, 0);
}


NLM_EXTERN CONNECTOR SOCK_CreateConnectorEx
(const Nlm_Char* host,
 Nlm_Uint2       port,
 Nlm_Uint4       conn_try,
 const void*     init_data,
 Nlm_Uint4       init_size)
{
  CONNECTOR       ccc = (SConnector    *)Nlm_MemNew(sizeof(SConnector    ));
  SSockConnector* xxx = (SSockConnector*)Nlm_MemNew(sizeof(SSockConnector));

  /* initialize internal data structures */
  xxx->sock = 0;
  xxx->host = Nlm_StringSave(host);
  xxx->port = port;
  xxx->conn_try  = conn_try  ? conn_try  : 1;

  xxx->init_size = init_data ? init_size : 0;
  if ( xxx->init_size ) {
    xxx->init_data = Nlm_MemNew(init_size);
    Nlm_MemCpy(xxx->init_data, init_data, xxx->init_size);
  } else
    xxx->init_data = 0;

  /* initialize handle */
  ccc->handle = xxx;

  /* initialize virtual table */ 
  ccc->vtable.get_type   = s_VT_GetType;
  ccc->vtable.connect    = s_VT_Connect;
  ccc->vtable.wait       = s_VT_Wait;
#ifdef IMPLEMENTED__CONN_WaitAsync
  ccc->vtable.wait_async = s_VT_WaitAsync;
#endif
  ccc->vtable.write      = s_VT_Write;
  ccc->vtable.flush      = s_VT_Flush;
  ccc->vtable.read       = s_VT_Read;
  ccc->vtable.close      = s_VT_Close;

  /* done */
  return ccc;
}


#ifdef TEST_MODULE__CON_SOCK

#include <conntest.h>

extern
#ifdef __cplusplus
"C"
#endif
Int2 Main(void)
{
  STimeout  timeout;
  CONNECTOR connector;
  FILE*     log_file;

#define MIN_PORT 5001
  Int4        argc = Nlm_GetArgc();
  CharPtr*    argv = Nlm_GetArgv();
  const Char* host;
  Uint2       port;
  Uint4       conn_try;

#ifdef WIN16
  {{ /* a kludge to make sure the "vibwndws.c"(MainProc) get linked */
    extern void Nlm_Metronome(Nlm_VoidPtr actn);  Nlm_Metronome(0);
  }}
#endif

  Nlm_ErrSetOpts(ERR_TEE, ERR_LOG_ON);
  ErrSetLogLevel(SEV_INFO);
  ErrSetMessageLevel(SEV_INFO);
  log_file = FileOpen("con_sock.out", "wb");

  /* defaults */
  host         = 0;
  port         = 0;
  conn_try     = 2;
  timeout.sec  = 5;
  timeout.usec = 123456;

  /* parse cmd.-line args */
  switch ( argc ) {
  case 5: { /* timeout */
    float fff = 0;
    if (sscanf(argv[4], "%f", &fff) != 1  ||  fff < 0)
      break;
    timeout.sec  = (Uint4)fff;
    timeout.usec = (Uint4)((fff - timeout.sec) * 1000000);
  }
  case 4: { /* conn_try  */
    long lll;
    if (sscanf(argv[3], "%ld", &lll) != 1  ||  lll <= 0)
      break;
    conn_try = (Uint4)lll;
  }
  case 3: { /* host, port */
    int iii;
    if (sscanf(argv[2], "%d", &iii) != 1  ||
        iii < MIN_PORT  ||  INT2_MAX <= iii)
      break;
    port = (Uint2)iii;

    if ( !*argv[1] )
      break;
    host = argv[1];
  }
  default:
    break;
  }

  if ( !host ) {
    ErrPostEx(SEV_ERROR, 0, 0,
              "Usage: %s <host> <port> [conn_try] [timeout]\n"
              "  where <port> not less than %d; timeout is a float(in sec)",
              argv[0], (int)MIN_PORT);
    return 1;
  }

  ErrPostEx(SEV_INFO, 0, 0,
            "Starting the CON_SOCK test...\n"
            "%s:%d,  timeout = %lu.%06lu, max # of retry = %lu\n",
            host, (int)port,
            (unsigned long)timeout.sec, (unsigned long)timeout.usec,
            (unsigned long)conn_try);

  connector = SOCK_CreateConnector(host, port, conn_try);
  Ncbi_TestConnector(connector, &timeout, log_file,
                     TESTCONN_SINGLE_BOUNCE_PRINT);

  connector = SOCK_CreateConnector(host, port, conn_try);
  Ncbi_TestConnector(connector, &timeout, log_file,
                     TESTCONN_SINGLE_BOUNCE_CHECK);

  connector = SOCK_CreateConnector(host, port, conn_try);
  Ncbi_TestConnector(connector, &timeout, log_file,
                     TESTCONN_ALL);

  FileClose(log_file);
  return 0;
}

#endif  /* TEST_MODULE__CON_SOCK */
