/*  $RCSfile: ncbibuf.c,v $  $Revision: 6.1 $  $Date: 1998/04/14 15:34:22 $
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
*   Memory-resided FIFO storage area(to be used e.g. in I/O buffering)
*
* --------------------------------------------------------------------------
* $Log: ncbibuf.c,v $
* Revision 6.1  1998/04/14 15:34:22  vakatov
* Initial revision
*
* ==========================================================================
*/

#include <ncbi.h>
#include <ncbibuf.h>

#define MIN_CHUNK_SIZE 512

typedef struct BufChunkTag {
  struct BufChunkTag* next;
  Uint4 size;       /* of data */
  Uint4 alloc_size; /* maximum avail.(allocated) size of "data" */
  Char  data[1];
} BufChunk;

typedef struct Nlm_BufferTag {
  BufChunk* list;
  Uint4     n_skip; /* # of bytes already "removed" from the 1st data chunk */
} Nlm_Buffer;


NLM_EXTERN Uint4 BUF_Size(BUF buf)
{
  Uint4  size;
  BufChunk  *pChunk;
  if ( !buf )
    return 0;

  for (size = 0, pChunk = buf->list;  pChunk;  pChunk = pChunk->next) {
    size += pChunk->size;
  }
  ASSERT((size > buf->n_skip) || (size == buf->n_skip && !buf->list && !size));
  size -= buf->n_skip;

  return size;
}


NLM_EXTERN Boolean BUF_Write(BUF* pBuf, const void* data, Uint4 size)
{
  BufChunk *pChunk, *pTail;
  if ( !size )
    return TRUE;

  /* create the buffer internals, if not created yet */
  if ( !*pBuf ) {
    *pBuf = (Nlm_Buffer*)MemNew(sizeof(Nlm_Buffer));
    if ( !*pBuf )
      return FALSE;
  }

  /* find the last allocated chunk */
  for (pTail = (*pBuf)->list;  pTail  &&  pTail->next;  pTail = pTail->next);

  /* write to an unfilled space of the last allocated chunk, if any */
  if (pTail  &&  pTail->size != pTail->alloc_size) {
    Uint4 n_avail = pTail->alloc_size - pTail->size;
    Uint4 n_write = (size <= n_avail) ? size : n_avail;
    ASSERT( pTail->size < pTail->alloc_size );
    MemCpy(pTail->data + pTail->size, data, n_write);
    pTail->size += n_write;
    size -= n_write;
    data = (char*)data + n_write;
  }

  /* allocate and write to the new chunk, if necessary */
  if ( size ) {
    Uint4 alloc_size =
      ((size + MIN_CHUNK_SIZE - 1) / MIN_CHUNK_SIZE) * MIN_CHUNK_SIZE;
    pChunk = (BufChunk*)MemNew(sizeof(BufChunk) - 1 + alloc_size);
    if ( !pChunk )
      return FALSE;
    pChunk->alloc_size = alloc_size;
    MemCpy(pChunk->data, data, size);
    pChunk->size = size;

    /* add the new chunk to the buffer list */
    if ( pTail )
      pTail->next = pChunk;
    else
      (*pBuf)->list = pChunk;
  }
  return TRUE;
}


NLM_EXTERN Uint4 BUF_Peek(BUF buf, void* data, Uint4 size)
{
  Uint4     n_todo = size;
  Uint4     n_skip;
  BufChunk *pChunk;

  if (!data  ||  !size  ||  !buf  ||  !buf->list)
    return 0;

  n_skip = buf->n_skip;
  for (pChunk = buf->list, n_skip = buf->n_skip;
       n_todo  &&  pChunk;
       pChunk = pChunk->next, n_skip = 0) {
    Uint4 n_copy = pChunk->size - n_skip;
    if (n_copy > n_todo)
      n_copy = n_todo;

    ASSERT( n_skip < pChunk->size );
    MemCpy(data, (char*)pChunk->data + n_skip, n_copy);
    data = (char*)data + n_copy;
    n_todo -= n_copy;
  }

  return (Uint4)(size - n_todo);
}


NLM_EXTERN Uint4 BUF_Read(BUF buf, void* data, Uint4 size)
{
  Uint4 n_todo;
  if (!buf  ||  !size)
    return 0;

  /* peek to the callers data buffer, if non-NULL */
  if ( data )
    size = BUF_Peek(buf, data, size);

  /* remove the read data from the buffer */ 
  n_todo = size;
  while (n_todo  &&  buf->list) {
    Uint4 n_avail = buf->list->size - buf->n_skip;
    if (n_todo >= n_avail) { /* discard the whole chunk */
      BufChunk *pChunk = buf->list;
      buf->list = pChunk->next;
      MemFree(pChunk);
      n_todo -= n_avail;
      buf->n_skip = 0;
    } else { /* discard some of the chunk data */
      buf->n_skip += n_todo;
      n_todo = 0;
    }
  }

  return (Uint4)(size - n_todo);
}


NLM_EXTERN BUF BUF_Destroy(BUF buf)
{
  if ( !buf )
    return 0;

  while ( buf->list ) {
    BufChunk *pChunk = buf->list;
    buf->list = pChunk->next;
    MemFree(pChunk);
  }

  MemFree(buf);
  return 0;
}


#ifdef TEST_MODULE__NCBIBUF

static Uint4 s_Rand(void)
{ /* a uniform random number generator */
  static Uint4 s_Random = 1;
  s_Random = s_Random * 1103515245 + 12345;
  return (Uint4)(s_Random / 65536) % 32768;
}


Int2 Main()
{ /* test application */
#define X_MAX_N_IO  3
#define X_MAX_READ  MIN_CHUNK_SIZE * 3
#define X_TIMES     ((Uint4)(s_Rand() % X_MAX_N_IO))
#define X_BYTES     ((Uint4)(s_Rand() % X_MAX_READ))

  BUF buf = 0;
  Boolean do_loop = TRUE;

  FILE *fin  = FileOpen("stdin",  "rb");
  FILE *fout = FileOpen("stdout", "wb");

  /* setup the error posting */
  Nlm_ErrSetOpts(ERR_TEE, ERR_LOG_ON);
  ErrSetLogLevel(SEV_INFO);
  ErrSetMessageLevel(SEV_INFO);
  VERIFY( Nlm_ErrSetLog("ncbibuf.log") );

  /* read up to the very end of input stream */
  while ( do_loop ) {
    Char charbuf[X_MAX_READ];
    Uint4 i;
    Uint4 n_times;

    /* read from the input stream, write to the NCBI IO-buf */
    n_times = X_TIMES;
    for (i = 0;  i < n_times;  i++) {
      Uint4 n_bytes = X_BYTES;
      if ( !n_bytes )
        continue;
      n_bytes = (Uint4)FileRead(charbuf, 1, (size_t)n_bytes, fin);
      ErrPostEx(SEV_INFO, 0, 0, "[FileRead] %lu", (unsigned long)n_bytes);
      if ( !n_bytes ) {
        do_loop = FALSE; /* end of the input stream */
        break;
      }
      VERIFY( BUF_Write(&buf, charbuf, n_bytes) );
      ErrPostEx(SEV_INFO, 0, 0, "[BUF_Write] %lu", (unsigned long)n_bytes);
    }

    /* peek & read from the NCBI IO-buf, write to the output stream */
    n_times = X_TIMES;
    for (i = 0;  i < n_times  &&  BUF_Size(buf);  i++) {
      Boolean do_peek = (Boolean)(s_Rand() % 2 == 0);
      Uint4   n_peek;
      Uint4   n_bytes = X_BYTES;
      if ( !n_bytes )
        continue;

      /* peek from the NCBI IO-buf */
      if ( do_peek ) {
        Uint4 j, n_peek_times = s_Rand() % 3 + 1;
        for (j = 0;  j < n_peek_times;  j++) {
          n_peek = BUF_Peek(buf, charbuf, n_bytes);
          ErrPostEx(SEV_INFO, 0, 0, "[BUF_Peek] %lu", (unsigned long)n_peek);
        }
      }

      /* read(or just discard) the data */
      if (do_peek  &&  s_Rand() % 2 == 0)
        n_bytes = BUF_Read(buf, 0, n_bytes);
      else
        n_bytes = BUF_Read(buf, charbuf, n_bytes);

      ErrPostEx(SEV_INFO, 0, 0, "[BUF_Read] %lu", (unsigned long)n_bytes);
      ASSERT( !do_peek  ||  n_bytes == n_peek );

      /* write the read data to the output stream */
      VERIFY( n_bytes == (Uint4)FileWrite(charbuf, 1, (size_t)n_bytes, fout) );
      ErrPostEx(SEV_INFO, 0, 0, "[FileWrite] %lu", (unsigned long)n_bytes);
    }
  }

  /* flush the IO-buf to the output stream */
  while ( BUF_Size(buf) ) {
    Char charbuf[1024];
    Uint4 n_bytes = BUF_Read(buf, charbuf, sizeof(charbuf));
    ErrPostEx(SEV_INFO, 0, 0, "[BUF_Read/flush] %lu", (unsigned long)n_bytes);
    ASSERT( n_bytes );
    VERIFY( n_bytes == (Uint4)FileWrite(charbuf, 1, (size_t)n_bytes, fout) );
    ErrPostEx(SEV_INFO, 0, 0, "[FileWrite/flush] %lu", (unsigned long)n_bytes);
  }

  /* cleanup */
  BUF_Destroy(buf);
  FileClose(fout);
  FileClose(fin);
  
  return 0;
}
#endif /* TEST_MODULE__NCBIBUF */

/* EOF */

