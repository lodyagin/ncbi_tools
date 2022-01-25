/*  $Id: hspstream_queue.h,v 1.2 2004/06/08 17:46:35 dondosha Exp $
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
 * Author:  Ilya Dondoshansky
 *
 */

/** @file hspstream_queue.h
 * Implementation of the BlastHSPStream interface for producing results on the
 * fly.
 */

#ifndef HSPSTREAM_QUEUE_H
#define HSPSTREAM_QUEUE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <ncbithr.h>
#include <algo/blast/core/blast_options.h>
#include <algo/blast/core/blast_hits.h>
#include <algo/blast/core/blast_seqsrc.h>
#include <algo/blast/core/blast_hspstream.h>

/** Data structure for the queue implementation of BlastHSPStream */
typedef struct BlastHSPListQueueData {
   ListNode* m_queueStart;
   ListNode* m_queueEnd;
   Boolean   m_writingDone;
   TNlmMutex m_resultsMutex;
   TNlmSemaphore m_resultsSema;
} BlastHSPListQueueData;

/** Function to initialize the queue implementation of BlastHSPStream */
BlastHSPStream* Blast_HSPListQueueInit(void);

#ifdef __cplusplus
}
#endif

#endif /* HSPSTREAM_QUEUE_H */
