/*  $Id: ncbibuf.c,v 6.5 1999/08/17 22:30:21 vakatov Exp $
 * ==========================================================================
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
 * ==========================================================================
 *
 * Author:  Denis Vakatov
 *
 * File Description:
 *   Memory-resident FIFO storage area (to be used e.g. in I/O buffering)
 *
 *   This is just a back-compatibility interface("Nlm_*") to the real
 *   BUF API located in "ncbi_buffer.[ch]".
 *   Unlike the "real" BUF API, this API uses:
 *    a) "Nlm_" name prefix for structures, types and functions;
 *    b) "Nlm_*" fixed size integer types like "Nlm_Uint4";
 *    c) "Nlm_Boolean" rather than a native "int" for the boolean type;
 *    d) [MSWIN] "NLM_EXTERN" rather than just "extern" to ease the compilation
 *       for MSWIN DLL.
 *
 * --------------------------------------------------------------------------
 * $Log: ncbibuf.c,v $
 * Revision 6.5  1999/08/17 22:30:21  vakatov
 * Use "Nlm_BUF" to avoid name clash with "BUF" in "ncbi_buffer.h" when
 * compiling with MipsPro IRIX compiler
 *
 * Revision 6.4  1999/08/17 19:47:24  vakatov
 * Moved all real code from NCBIBUF to NCBI_BUFFER;  the code has been cleaned
 * from the NCBI C toolkit specific types and API calls.
 * NCBIBUF module still exists for the backward compatibility -- it
 * provides old NCBI-wise interface.
 *
 * ==========================================================================
 */


#include <ncbibuf.h>

/* ...and undefine "Nlm_BUF_*" to allow to access <ncbi_buffer.h> API */
#undef BUF
#undef BUF_SetChunkSize
#undef BUF_Size
#undef BUF_Write
#undef BUF_Peek
#undef BUF_Read
#undef BUF_PushBack
#undef BUF_Destroy

/* this is the only place where both <ncbibuf.h> and <ncbi_buffer.h> can
 * be #include'd in one source module! */
#undef NCBIBUF__H

#include <ncbi_buffer.h>


NLM_EXTERN Nlm_Uint4 Nlm_BUF_SetChunkSize
(BUF* pBuf, Nlm_Uint4 chunk_size)
{
  return (Nlm_Uint4) BUF_SetChunkSize(pBuf, chunk_size);
}

NLM_EXTERN Nlm_Uint4 Nlm_BUF_Size
(BUF buf)
{
  return (Nlm_Uint4) BUF_Size(buf);
}

NLM_EXTERN Nlm_Boolean Nlm_BUF_Write
(BUF* pBuf, const void* data, Nlm_Uint4 size)
{
  return (Nlm_Boolean) BUF_Write(pBuf, data, size);
}

NLM_EXTERN Nlm_Boolean Nlm_BUF_PushBack
(BUF* pBuf, const void* data, Nlm_Uint4 size)
{
  return (Nlm_Boolean) BUF_PushBack(pBuf, data, size);
}


NLM_EXTERN Nlm_Uint4 Nlm_BUF_Peek
(BUF buf, void* data, Nlm_Uint4 size)
{
  return (Nlm_Uint4) BUF_Peek(buf, data, size);
}


NLM_EXTERN Nlm_Uint4 Nlm_BUF_Read
(BUF buf, void* data, Nlm_Uint4 size)
{
   return (Nlm_Uint4) BUF_Read(buf, data, size);
}


NLM_EXTERN BUF Nlm_BUF_Destroy
(BUF buf)
{
  return BUF_Destroy(buf);
}





#ifdef TEST_MODULE__NCBIBUF

#include <ncbi.h>

static Nlm_Uint4 s_Rand(void)
{ /* a uniform random number generator */
  static Nlm_Uint4 s_Random = 1;
  s_Random = s_Random * 1103515245 + 12345;
  return (Nlm_Uint4)(s_Random / 65536) % 32768;
}


Int2 Main(void)
{ /* test application */
#define X_MAX_N_IO  4
#define X_MAX_READ  BUF_DEF_CHUNK_SIZE * 3
#define X_TIMES     ((Nlm_Uint4)(s_Rand() % X_MAX_N_IO))
#define X_BYTES     ((Nlm_Uint4)(s_Rand() % X_MAX_READ))

  BUF buf = 0;
  Nlm_Boolean do_loop = TRUE;

  FILE *fin  = FileOpen("stdin",  "rb");
  FILE *fout = FileOpen("stdout", "wb");

  /* setup the error posting */
  Nlm_ErrSetOpts(ERR_TEE, ERR_LOG_ON);
  ErrSetLogLevel(SEV_INFO);
  ErrSetMessageLevel(SEV_INFO);
  VERIFY( Nlm_ErrSetLog("ncbibuf.log") );

  /* a simple test */
  {{
    char  charbuf[128];
    VERIFY( Nlm_BUF_PushBack(&buf, (const char*)"0", 1) );
    VERIFY( Nlm_BUF_Write(&buf, (const char*)"1", 1) );
    VERIFY( Nlm_BUF_Peek(buf, charbuf, sizeof(charbuf)) == 2);
    VERIFY( Nlm_BUF_PushBack(&buf, (const char*)"BB", 2) );
    VERIFY( Nlm_BUF_PushBack(&buf, (const char*)"aa", 2) );
    VERIFY( Nlm_BUF_Write(&buf, (const char*)"23", 3) );
    VERIFY( Nlm_BUF_Read(buf, charbuf, sizeof(charbuf)) == 9);
    ASSERT( StrCmp(charbuf, (const char*)"aaBB0123") == 0 );
    buf = Nlm_BUF_Destroy(buf);
  }}

  /* read up to the very end of input stream */
  while ( do_loop ) {
    Char charbuf[X_MAX_READ];
    Nlm_Uint4 i;
    Nlm_Uint4 n_times;

    /* read from the input stream, write to the NCBI IO-buf */
    n_times = X_TIMES;
    for (i = 0;  i < n_times;  i++) {
      Nlm_Uint4 n_bytes = X_BYTES;
      if ( !n_bytes )
        continue;
      n_bytes = (Nlm_Uint4)FileRead(charbuf, 1, (size_t)n_bytes, fin);
      ErrPostEx(SEV_INFO, 0, 0, "[FileRead] %lu", (unsigned long)n_bytes);
      if ( !n_bytes ) {
        do_loop = FALSE; /* end of the input stream */
        break;
      }
      VERIFY( Nlm_BUF_Write(&buf, charbuf, n_bytes) );
      ErrPostEx(SEV_INFO, 0, 0, "[BUF_Write] %lu", (unsigned long)n_bytes);
    }

    /* peek & read from the NCBI IO-buf, write to the output stream */
    n_times = X_TIMES;
    for (i = 0;  i < n_times  &&  Nlm_BUF_Size(buf);  i++) {
      Nlm_Boolean do_peek = (Nlm_Boolean)(s_Rand() % 2 == 0);
      Nlm_Uint4   n_peek;
      Nlm_Uint4   n_bytes = X_BYTES;
      if ( !n_bytes )
        continue;

      /* peek from the NCBI IO-buf */
      if ( do_peek ) {
        Nlm_Uint4 j, n_peek_times = s_Rand() % 3 + 1;
        for (j = 0;  j < n_peek_times;  j++) {
          n_peek = Nlm_BUF_Peek(buf, charbuf, n_bytes);
          ErrPostEx(SEV_INFO, 0, 0, "[BUF_Peek] %lu", (unsigned long)n_peek);
        }
      }

      /* read(or just discard) the data */
      if (do_peek  &&  s_Rand() % 2 == 0)
        n_bytes = Nlm_BUF_Read(buf, 0, n_bytes);
      else
        n_bytes = Nlm_BUF_Read(buf, charbuf, n_bytes);

      ErrPostEx(SEV_INFO, 0, 0, "[BUF_Read] %lu", (unsigned long)n_bytes);
      ASSERT( !do_peek  ||  n_bytes == n_peek );

      /* fake push back */
      if (s_Rand() % 3 == 0) {
        Nlm_Uint4 n_pushback = s_Rand() % n_bytes;
        VERIFY(Nlm_BUF_PushBack(&buf, charbuf + n_bytes - n_pushback, n_pushback));
        VERIFY(Nlm_BUF_Read(buf, charbuf + n_bytes - n_pushback, n_pushback));
      }
      
      /* write the read data to the output stream */
      VERIFY( n_bytes == (Nlm_Uint4)FileWrite(charbuf, 1, (size_t)n_bytes, fout) );
      ErrPostEx(SEV_INFO, 0, 0, "[FileWrite] %lu", (unsigned long)n_bytes);
    }
  }

  /* flush the IO-buf to the output stream */
  while ( Nlm_BUF_Size(buf) ) {
    Char charbuf[256];
    Nlm_Uint4 n_bytes = Nlm_BUF_Read(buf, charbuf, sizeof(charbuf));
    {{
      char      tmp[sizeof(charbuf)];
      Nlm_Uint4 n_pushback = s_Rand() % 64;
      if (n_pushback > n_bytes)
        n_pushback = n_bytes;
      VERIFY( Nlm_BUF_PushBack(&buf, charbuf + n_bytes - n_pushback, n_pushback) );
      VERIFY( Nlm_BUF_Read(buf, tmp, n_pushback)
              == n_pushback );
      MemCpy(charbuf + n_bytes - n_pushback, tmp, n_pushback);
    }}
    ErrPostEx(SEV_INFO, 0, 0, "[BUF_Read/flush] %lu", (unsigned long)n_bytes);
    ASSERT( n_bytes );
    VERIFY( n_bytes == (Nlm_Uint4)FileWrite(charbuf, 1, (size_t)n_bytes, fout) );
    ErrPostEx(SEV_INFO, 0, 0, "[FileWrite/flush] %lu", (unsigned long)n_bytes);
  }

  /* cleanup */
  Nlm_BUF_Destroy(buf);
  FileClose(fout);
  FileClose(fin);
  
  return 0;
}
#endif /* TEST_MODULE__NCBIBUF */

/* EOF */
