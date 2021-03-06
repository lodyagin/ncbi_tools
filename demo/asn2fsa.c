/*   asn2fsa.c
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
* File Name:  asn2fsa.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   3/4/04
*
* $Revision: 1.72 $
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

#include <ncbi.h>
#include <objall.h>
#include <objsset.h>
#include <objsub.h>
#include <objfdef.h>
#include <seqport.h>
#include <sequtil.h>
#include <sqnutils.h>
#include <subutil.h>
#include <tofasta.h>
#include <gather.h>
#include <explore.h>
#include <lsqfetch.h>
#include <readdb.h>
#include <pmfapi.h>
#ifdef INTERNAL_NCBI_ASN2FSA
#include <accpubseq.h>
#endif
#include <connect/ncbi_gnutls.h>

#define ASN2FSA_APP_VER "6.1"

CharPtr ASN2FSA_APPLICATION = ASN2FSA_APP_VER;

static ValNodePtr  requested_uid_list = NULL;
static TNlmMutex   requested_uid_mutex = NULL;

static ValNodePtr  locked_bsp_list = NULL;
static TNlmMutex   locked_bsp_mutex = NULL;

static void AddUidToQueue (
  SeqIdPtr sip
)

{
  ValNodePtr  last = NULL, vnp;
  Int4        ret;
  Int4        uid;

  if (sip == NULL || sip->choice != SEQID_GI) return;
  uid = (Int4) sip->data.intvalue;
  if (uid < 1) return;

  ret = NlmMutexLockEx (&requested_uid_mutex);
  if (ret) {
    ErrPostEx (SEV_FATAL, 0, 0, "AddUidToQueue mutex failed [%ld]", (long) ret);
    return;
  }

  /* check against uids already in queue */

  last = NULL;
  for (vnp = requested_uid_list; vnp != NULL; vnp = vnp->next) {
    last = vnp;
    if ((Int4) vnp->data.intvalue == uid) break;
  }

  /* add uid to queue */

  if (vnp == NULL) {
    if (last != NULL) {
      vnp = ValNodeAddInt (&last, 0, uid);
      last = vnp;
    } else {
      requested_uid_list = ValNodeAddInt (NULL, 0, uid);
      last = requested_uid_list;
    }
  }

  NlmMutexUnlock (requested_uid_mutex);
}

static Int4 RemoveUidFromQueue (
  void
)

{
  Int4        ret, uid = 0;
  ValNodePtr  vnp;

  ret = NlmMutexLockEx (&requested_uid_mutex);
  if (ret) {
    ErrPostEx (SEV_FATAL, 0, 0, "RemoveUidFromQueue mutex failed [%ld]", (long) ret);
    return 0;
  }

  /* extract next requested uid from queue */

  if (requested_uid_list != NULL) {
    vnp = requested_uid_list;
    requested_uid_list = vnp->next;
    vnp->next = NULL;
    uid = (Int4) vnp->data.intvalue;
    ValNodeFree (vnp);
  }

  NlmMutexUnlock (requested_uid_mutex);

  return uid;
}

static void QueueFarSegments (SeqLocPtr slp)

{
  BioseqPtr   bsp;
  SeqLocPtr   loc;
  SeqIdPtr    sip;
  ValNodePtr  vnp;

  if (slp == NULL) return;

  sip = SeqLocId (slp);
  if (sip == NULL) {
    loc = SeqLocFindNext (slp, NULL);
    if (loc != NULL) {
      sip = SeqLocId (loc);
    }
  }
  if (sip == NULL) return;

  /* if packaged in record, no need to fetch it */

  if (BioseqFindCore (sip) != NULL) return;

  /* check against currently locked records */

  for (vnp = locked_bsp_list; vnp != NULL; vnp = vnp->next) {
    bsp = (BioseqPtr) vnp->data.ptrvalue;
    if (bsp == NULL) continue;
    if (SeqIdIn (sip, bsp->id)) return;
  }

  AddUidToQueue (sip);
}

static void QueueFarBioseqs (BioseqPtr bsp, Pointer userdata)

{
  DeltaSeqPtr  dsp;
  SeqLocPtr    slp = NULL;
  ValNode      vn;

  if (bsp == NULL) return;

  if (bsp->repr == Seq_repr_seg) {
    vn.choice = SEQLOC_MIX;
    vn.extended = 0;
    vn.data.ptrvalue = bsp->seq_ext;
    vn.next = NULL;
    while ((slp = SeqLocFindNext (&vn, slp)) != NULL) {
      if (slp != NULL && slp->choice != SEQLOC_NULL) {
        QueueFarSegments (slp);
      }
    }
  } else if (bsp->repr == Seq_repr_delta) {
    for (dsp = (DeltaSeqPtr) (bsp->seq_ext); dsp != NULL; dsp = dsp->next) {
      if (dsp->choice == 1) {
        slp = (SeqLocPtr) dsp->data.ptrvalue;
        if (slp != NULL && slp->choice != SEQLOC_NULL) {
          QueueFarSegments (slp);
        }
      }
    }
  }
}

static void AddBspToList (
  BioseqPtr bsp
)

{
  Int4        ret;
  ValNodePtr  vnp;

  if (bsp == NULL) return;

  ret = NlmMutexLockEx (&locked_bsp_mutex);
  if (ret) {
    ErrPostEx (SEV_FATAL, 0, 0, "AddBspToList mutex failed [%ld]", (long) ret);
    return;
  }

  vnp = ValNodeAddPointer (&locked_bsp_list, 0, (Pointer) bsp);

  NlmMutexUnlock (locked_bsp_mutex);
}

static ValNodePtr ExtractBspList (
  void
)

{
  Int4        ret;
  ValNodePtr  vnp;

  ret = NlmMutexLockEx (&locked_bsp_mutex);
  if (ret) {
    ErrPostEx (SEV_FATAL, 0, 0, "ExtractBspList mutex failed [%ld]", (long) ret);
    return NULL;
  }

  vnp = locked_bsp_list;
  locked_bsp_list = NULL;

  NlmMutexUnlock (locked_bsp_mutex);

  return vnp;
}

typedef struct fastaflags {
  Boolean  master_style;
  Boolean  html_spans;
  Boolean  expand_gaps;
  Boolean  use_dashes;
  Boolean  extended_ids;
  Boolean  far_genomic_qual;
  Boolean  qual_gap_is_zero;
  Boolean  automatic;
  Boolean  batch;
  Boolean  binary;
  Boolean  compressed;
  Boolean  lock;
  Boolean  useThreads;
  Boolean  usePUBSEQ;
  Boolean  useBLAST;
  CharPtr  blastdbname;
  Int2     type;
  Int2     linelen;
  Boolean  failed;
  FILE     *nt;
  FILE     *aa;
  FILE     *ql;
  FILE     *fr;
  FILE     *logfp;
  Boolean  debugging;
} FastaFlagData, PNTR FastaFlagPtr;

static VoidPtr DoAsyncLookup (
  VoidPtr arg
)

{
  BioseqPtr     bsp;
  FastaFlagPtr  ffp;
  Int4          uid;
  ValNode       vn;

  ffp = (FastaFlagPtr) arg;
  if (ffp == NULL) return NULL;

#ifdef INTERNAL_NCBI_ASN2FSA
  if (ffp->usePUBSEQ) {
    PUBSEQInit ();
  }
#endif
  if (ffp->useBLAST) {
    ReadDBBioseqFetchEnable ("asn2fsa", ffp->blastdbname, TRUE, FALSE);
  }

  MemSet ((Pointer) &vn, 0, sizeof (ValNode));

  uid = RemoveUidFromQueue ();
  while (uid > 0) {

    vn.choice = SEQID_GI;
    vn.data.intvalue = uid;
    vn.next = NULL;

    bsp = BioseqLockById (&vn);
    if (bsp != NULL) {
      AddBspToList (bsp);
    }

    uid = RemoveUidFromQueue ();
  }

  if (ffp->useBLAST) {
    ReadDBBioseqFetchDisable ();
  }
#ifdef INTERNAL_NCBI_ASN2FSA
  if (ffp->usePUBSEQ) {
    PUBSEQFini ();
  }
#endif

  return NULL;
}

#define NUM_ASYNC_LOOKUP_THREADS 5

static void ProcessAsyncLookups (
  FastaFlagPtr ffp
)

{
  Int2        i;
  VoidPtr     status;
  TNlmThread  thds [NUM_ASYNC_LOOKUP_THREADS];

  /* spawn several threads for individual BioseqLockById requests */

  for (i = 0; i < NUM_ASYNC_LOOKUP_THREADS; i++) {
    thds [i] = NlmThreadCreate (DoAsyncLookup, (Pointer) ffp);
  }

  /* wait for all fetching threads to complete */

  for (i = 0; i < NUM_ASYNC_LOOKUP_THREADS; i++) {
    NlmThreadJoin (thds [i], &status);
  }
}

static ValNodePtr AsyncLockFarComponents (
  SeqEntryPtr sep,
  FastaFlagPtr ffp
)

{
  BioseqPtr    bsp;
  ValNodePtr   bsplist = NULL, sublist, vnp;
  SeqEntryPtr  oldsep;

  if (sep == NULL || ffp == NULL) return NULL;
  oldsep = SeqEntrySetScope (sep);

  /* add far uids to queue */

  VisitBioseqsInSep (sep, NULL, QueueFarBioseqs);

  /* fetching from uid list using several threads */

  ProcessAsyncLookups (ffp);

  sublist = ExtractBspList ();

  /* take list, look for seg or delta, recurse */

  while (sublist != NULL) {
    for (vnp = sublist; vnp != NULL; vnp = vnp->next) {
      bsp = (BioseqPtr) vnp->data.ptrvalue;
      if (bsp == NULL) continue;
      QueueFarBioseqs (bsp, NULL);
    }

    ValNodeLink (&bsplist, sublist);
    sublist = NULL;

    ProcessAsyncLookups (ffp);

    sublist = ExtractBspList ();
  }

  SeqEntrySetScope (oldsep);
  return bsplist;
}

static ValNodePtr DoLockFarComponents (
  SeqEntryPtr sep,
  FastaFlagPtr ffp
)

{
  ValNodePtr  rsult;
  time_t      start_time, stop_time;

  /*
  if (NlmThreadsAvailable () && ffp->useThreads) {
    return AsyncLockFarComponents (sep);
  }

  return LockFarComponents (sep);
  */

  start_time = GetSecs ();

  if (NlmThreadsAvailable () && ffp->useThreads) {
    rsult = AsyncLockFarComponents (sep, ffp);
  } else if (ffp->useThreads) {
    Message (MSG_POST, "Threads not available in this executable");
    rsult = LockFarComponents (sep);
  } else {
    rsult = LockFarComponents (sep);
  }

  stop_time = GetSecs ();

  return rsult;
}

static Boolean SegHasParts (
  BioseqPtr bsp
)

{
  BioseqSetPtr  bssp;
  SeqEntryPtr   sep;

  if (bsp == NULL || bsp->repr != Seq_repr_seg) return FALSE;
  sep = bsp->seqentry;
  if (sep == NULL) return FALSE;
  sep = sep->next;
  if (sep == NULL || (! IS_Bioseq_set (sep))) return FALSE;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp != NULL && bssp->_class == BioseqseqSet_class_parts) return TRUE;
  return FALSE;
}

static void CacheFarComponents (
  FastaFlagPtr ffp,
  ValNodePtr bsplist
)

{
  BioseqPtr   bsp;
  Uint2       entityID;
  ValNodePtr  vnp;

  if (ffp == NULL || ffp->fr == NULL || bsplist == NULL) return;

  for (vnp = bsplist; vnp != NULL; vnp = vnp->next) {
    bsp = (BioseqPtr) vnp->data.ptrvalue;
    if (bsp == NULL) continue;

    /* cache raw and constructed, near segmented, and delta literal */

    switch (bsp->repr) {
        case Seq_repr_raw :
        case Seq_repr_const :
          if (BioseqFastaStream (bsp, ffp->fr, 0, ffp->linelen, 0, 0, TRUE) < 0) {
            ffp->failed = TRUE;
          }
          break;
        case Seq_repr_seg :
          entityID = ObjMgrGetEntityIDForPointer (bsp);
          AssignIDsInEntity (entityID, 0, NULL);
          if (SegHasParts (bsp)) {
            if (BioseqFastaStream (bsp, ffp->fr, 0, ffp->linelen, 0, 0, TRUE) < 0) {
              ffp->failed = TRUE;
            }
          }
          break;
        case Seq_repr_delta :
          if (DeltaLitOnly (bsp)) {
            if (BioseqFastaStream (bsp, ffp->fr, 0, ffp->linelen, 0, 0, TRUE) < 0) {
              ffp->failed = TRUE;
            }
          }
          break;
        default :
          break;
    }
  }
}

static void PrintQualProc (
  CharPtr buf,
  Uint4 buflen,
  Pointer userdata
)

{
  FILE  *fp;

  fp = (FILE*) userdata;
  fprintf (fp, "%s", buf);
}

static void PrintQualScores (
  BioseqPtr bsp,
  Pointer userdata
)

{
  FastaFlagPtr  ffp;

  ffp = (FastaFlagPtr) userdata;
  if (bsp == NULL || ffp == NULL) return;
  if (! ISA_na (bsp->mol)) return;

  if (ffp->far_genomic_qual) {
    PrintQualityScoresForContig (bsp, ffp->qual_gap_is_zero, ffp->ql);
  } else {
    PrintQualityScoresToBuffer (bsp, ffp->qual_gap_is_zero, ffp->ql, PrintQualProc);
  }
}

static void ProcessSingleRecord (
  CharPtr directory,
  CharPtr base,
  CharPtr suffix,
  FastaFlagPtr ffp
)

{
  AsnIoPtr       aip;
  BioseqPtr      bsp;
  ValNodePtr     bsplist;
  BioseqSetPtr   bssp;
  Pointer        dataptr = NULL;
  Uint2          datatype, entityID = 0;
  Char           file /* [FILENAME_MAX] */ [PATH_MAX], path [PATH_MAX];
  StreamFlgType  flags = STREAM_CORRECT_INVAL;
  FILE           *fp;
  ObjMgrPtr      omp;
  SeqEntryPtr    sep;

  if (ffp == NULL) return;

  if (base == NULL) {
    base = "";
  }
  if (suffix == NULL) {
    suffix = "";
  }
  StringNCpy_0 (path, directory, sizeof (path));
  sprintf (file, "%s%s", base, suffix);
  FileBuildPath (path, NULL, file);

  if (StringHasNoText (path)) return;

  if (ffp->type == 1) {
    fp = FileOpen (path, "r");
    if (fp == NULL) {
      Message (MSG_POSTERR, "Failed to open '%s'", path);
      return;
    }

    dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, FALSE, FALSE, TRUE, FALSE);

    FileClose (fp);

    entityID = ObjMgrRegister (datatype, dataptr);

  } else if (ffp->type >= 2 && ffp->type <= 5) {
    aip = AsnIoOpen (path, ffp->binary? "rb" : "r");
    if (aip == NULL) {
      Message (MSG_POSTERR, "AsnIoOpen failed for input file '%s'", path);
      return;
    }

    SeqMgrHoldIndexing (TRUE);
    switch (ffp->type) {
      case 2 :
        dataptr = (Pointer) SeqEntryAsnRead (aip, NULL);
        datatype = OBJ_SEQENTRY;
        break;
      case 3 :
        dataptr = (Pointer) BioseqAsnRead (aip, NULL);
        datatype = OBJ_BIOSEQ;
        break;
      case 4 :
        dataptr = (Pointer) BioseqSetAsnRead (aip, NULL);
        datatype = OBJ_BIOSEQSET;
        break;
      case 5 :
        dataptr = (Pointer) SeqSubmitAsnRead (aip, NULL);
        datatype = OBJ_SEQSUB;
        break;
      default :
        break;
    }
    SeqMgrHoldIndexing (FALSE);

    AsnIoClose (aip);

    entityID = ObjMgrRegister (datatype, dataptr);

  } else {
    Message (MSG_POSTERR, "Input format type '%d' unrecognized", (int) ffp->type);
    return;
  }

  if (entityID < 1 || dataptr == NULL) {
    Message (MSG_POSTERR, "Data read failed for input file '%s'", path);
    return;
  }

  if (datatype == OBJ_SEQSUB || datatype == OBJ_SEQENTRY ||
        datatype == OBJ_BIOSEQ || datatype == OBJ_BIOSEQSET) {

    sep = GetTopSeqEntryForEntityID (entityID);

    if (sep == NULL) {
      sep = SeqEntryNew ();
      if (sep != NULL) {
        if (datatype == OBJ_BIOSEQ) {
          bsp = (BioseqPtr) dataptr;
          sep->choice = 1;
          sep->data.ptrvalue = bsp;
          SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, sep);
        } else if (datatype == OBJ_BIOSEQSET) {
          bssp = (BioseqSetPtr) dataptr;
          sep->choice = 2;
          sep->data.ptrvalue = bssp;
          SeqMgrSeqEntry (SM_BIOSEQSET, (Pointer) bssp, sep);
        } else {
          sep = SeqEntryFree (sep);
        }
      }
      sep = GetTopSeqEntryForEntityID (entityID);
    }

    if (sep != NULL) {
      BasicSeqEntryCleanup (sep);

      if (ffp->expand_gaps && ffp->use_dashes) {
        flags |= EXPAND_GAPS_TO_DASHES;
      } else if (ffp->expand_gaps) {
        flags |= STREAM_EXPAND_GAPS;
      } else if (ffp->use_dashes) {
        flags |= GAP_TO_SINGLE_DASH;
      }
      if (ffp->html_spans) {
        flags |= STREAM_HTML_SPANS;
      }
      if (ffp->extended_ids) {
        flags |= STREAM_ALL_FASTA_IDS;
      }

      bsplist = NULL;
      if (ffp->lock) {
        bsplist = DoLockFarComponents (sep, ffp);
        if (bsplist != NULL && ffp->fr != NULL) {
          CacheFarComponents (ffp, bsplist);
        }
      }

      if (ffp->nt != NULL) {
        if (SeqEntryFastaStream (sep, ffp->nt, flags, ffp->linelen, 0, 0,
                                 TRUE, FALSE, ffp->master_style) < 0) {
          ffp->failed = TRUE;
        }
      }
      if (ffp->aa != NULL) {
        if (SeqEntryFastaStream (sep, ffp->aa, flags, ffp->linelen, 0, 0,
                                 FALSE, TRUE, ffp->master_style) < 0) {
          ffp->failed = TRUE;
        }
      }
      if (ffp->ql != NULL) {
        VisitBioseqsInSep (sep, (Pointer) ffp, PrintQualScores);
      }

      bsplist = UnlockFarComponents (bsplist);

    }
  } else {
    Message (MSG_POSTERR, "Datatype %d not recognized", (int) datatype);
  }

  ObjMgrFree (datatype, dataptr);

  omp = ObjMgrGet ();
  ObjMgrReapOne (omp);
  SeqMgrClearBioseqIndex ();
  ObjMgrFreeCache (0);
  FreeSeqIdGiCache ();

  SeqEntrySetScope (NULL);
}

static void ProcessMultipleRecord (
  CharPtr directory,
  CharPtr base,
  CharPtr suffix,
  FastaFlagPtr ffp
)

{
  AsnIoPtr       aip;
  AsnModulePtr   amp;
  AsnTypePtr     atp, atp_bss, atp_desc, atp_se;
  BioseqPtr      bsp;
  ValNodePtr     bsplist;
  Char           buf [64], file /* [FILENAME_MAX] */ [PATH_MAX], path [PATH_MAX], longest [64];
  StreamFlgType  flags = STREAM_CORRECT_INVAL;
  FILE           *fp;
  Int4           numrecords = 0;
  SeqEntryPtr    fsep, sep;
  ObjMgrPtr      omp;
  time_t         starttime, stoptime, worsttime;
#ifdef OS_UNIX
  Char           cmmd [256];
  CharPtr        gzcatprog;
  int            ret;
  Boolean        usedPopen = FALSE;
#endif

  if (ffp == NULL) return;

  if (base == NULL) {
    base = "";
  }
  if (suffix == NULL) {
    suffix = "";
  }
  StringNCpy_0 (path, directory, sizeof (path));
  sprintf (file, "%s%s", base, suffix);
  FileBuildPath (path, NULL, file);

  if (StringHasNoText (path)) return;

#ifndef OS_UNIX
  if (ffp->compressed) {
    Message (MSG_POSTERR, "Can only decompress on-the-fly on UNIX machines");
    return;
  }
#endif

  amp = AsnAllModPtr ();
  if (amp == NULL) {
    Message (MSG_POSTERR, "Unable to load AsnAllModPtr");
    return;
  }

  atp_bss = AsnFind ("Bioseq-set");
  if (atp_bss == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Bioseq-set");
    return;
  }

  atp_desc = AsnFind ("Bioseq-set.descr");
  if (atp_desc == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Bioseq-set.descr");
    return;
  }

  atp_se = AsnFind ("Bioseq-set.seq-set.E");
  if (atp_se == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Bioseq-set.seq-set.E");
    return;
  }

#ifdef OS_UNIX
  if (ffp->compressed) {
    gzcatprog = getenv ("NCBI_UNCOMPRESS_BINARY");
    if (gzcatprog != NULL) {
      sprintf (cmmd, "%s %s", gzcatprog, path);
    } else {
      ret = system ("gzcat -h >/dev/null 2>&1");
      if (ret == 0) {
        sprintf (cmmd, "gzcat %s", path);
      } else if (ret == -1) {
        Message (MSG_POSTERR, "Unable to fork or exec gzcat in ScanBioseqSetRelease");
        return;
      } else {
        ret = system ("zcat -h >/dev/null 2>&1");
        if (ret == 0) {
          sprintf (cmmd, "zcat %s", path);
        } else if (ret == -1) {
          Message (MSG_POSTERR, "Unable to fork or exec zcat in ScanBioseqSetRelease");
          return;
        } else {
          Message (MSG_POSTERR, "Unable to find zcat or gzcat in ScanBioseqSetRelease - please edit your PATH environment variable");
          return;
        }
      }
    }
    fp = popen (cmmd, /* ffp->binary? "rb" : */ "r");
    usedPopen = TRUE;
  } else {
    fp = FileOpen (path, ffp->binary? "rb" : "r");
  }
#else
  fp = FileOpen (path, ffp->binary? "rb" : "r");
#endif
  if (fp == NULL) {
    Message (MSG_POSTERR, "FileOpen failed for input file '%s'", path);
    return;
  }

  aip = AsnIoNew (ffp->binary? ASNIO_BIN_IN : ASNIO_TEXT_IN, fp, NULL, NULL, NULL);
  if (aip == NULL) {
    Message (MSG_ERROR, "AsnIoNew failed for input file '%s'", path);
    return;
  }

  atp = atp_bss;

  if (ffp->expand_gaps && ffp->use_dashes) {
    flags |= EXPAND_GAPS_TO_DASHES;
  } else if (ffp->expand_gaps) {
    flags |= STREAM_EXPAND_GAPS;
  } else if (ffp->use_dashes) {
    flags |= GAP_TO_SINGLE_DASH;
  }
  if (ffp->html_spans) {
    flags |= STREAM_HTML_SPANS;
  }
  if (ffp->extended_ids) {
    flags |= STREAM_ALL_FASTA_IDS;
  }

  longest [0] = '\0';
  worsttime = 0;

  while ((atp = AsnReadId (aip, amp, atp)) != NULL) {
    if (atp == atp_se) {

      SeqMgrHoldIndexing (TRUE);
      sep = SeqEntryAsnRead (aip, atp);
      SeqMgrHoldIndexing (FALSE);

      BasicSeqEntryCleanup (sep);

      starttime = GetSecs ();
      buf [0] = '\0';

      if (ffp->logfp != NULL) {
        fsep = FindNthBioseq (sep, 1);
        if (fsep != NULL && fsep->choice == 1) {
          bsp = (BioseqPtr) fsep->data.ptrvalue;
          if (bsp != NULL) {
            SeqIdWrite (bsp->id, buf, PRINTID_FASTA_LONG, sizeof (buf));
            fprintf (ffp->logfp, "%s\n", buf);
            fflush (ffp->logfp);
          }
        }
      }

      bsplist = NULL;
      if (ffp->lock) {
        bsplist = DoLockFarComponents (sep, ffp);
        if (bsplist != NULL && ffp->fr != NULL) {
          CacheFarComponents (ffp, bsplist);
        }
      }

      if (ffp->nt != NULL) {
        SeqEntryFastaStream (sep, ffp->nt, flags, ffp->linelen, 0, 0, TRUE, FALSE, ffp->master_style);
      }
      if (ffp->aa != NULL) {
        SeqEntryFastaStream (sep, ffp->aa, flags, ffp->linelen, 0, 0, FALSE, TRUE, ffp->master_style);
      }
      if (ffp->ql != NULL) {
        VisitBioseqsInSep (sep, (Pointer) ffp, PrintQualScores);
      }

      bsplist = UnlockFarComponents (bsplist);

      stoptime = GetSecs ();
      if (stoptime - starttime > worsttime && StringDoesHaveText (buf)) {
        worsttime = stoptime - starttime;
        StringCpy (longest, buf);
      }
      numrecords++;

      SeqEntryFree (sep);
      omp = ObjMgrGet ();
      ObjMgrReapOne (omp);
      SeqMgrClearBioseqIndex ();
      ObjMgrFreeCache (0);
      FreeSeqIdGiCache ();

      SeqEntrySetScope (NULL);
    } else {
      AsnReadVal (aip, atp, NULL);
    }
  }

  AsnIoFree (aip, FALSE);

#ifdef OS_UNIX
  if (usedPopen) {
    pclose (fp);
  } else {
    FileClose (fp);
  }
#else
  FileClose (fp);
#endif

  if (ffp->logfp != NULL && (! StringHasNoText (longest))) {
    fprintf (ffp->logfp, "Longest processing time %ld seconds on %s\n",
             (long) worsttime, longest);
    fprintf (ffp->logfp, "Total number of records %ld\n", (long) numrecords);
    fflush (ffp->logfp);
  }
}

static void FastaWrapper (
  SeqEntryPtr sep,
  Pointer userdata
)

{
  ValNodePtr     bsplist;
  FastaFlagPtr   ffp;
  StreamFlgType  flags = STREAM_CORRECT_INVAL;

  if (sep == NULL) return;
  ffp = (FastaFlagPtr) userdata;
  if (ffp == NULL) return;

  BasicSeqEntryCleanup (sep);

  if (ffp->expand_gaps && ffp->use_dashes) {
    flags |= EXPAND_GAPS_TO_DASHES;
  } else if (ffp->expand_gaps) {
    flags |= STREAM_EXPAND_GAPS;
  } else if (ffp->use_dashes) {
    flags |= GAP_TO_SINGLE_DASH;
  }
  if (ffp->html_spans) {
    flags |= STREAM_HTML_SPANS;
  }
  if (ffp->extended_ids) {
    flags |= STREAM_ALL_FASTA_IDS;
  }

  bsplist = NULL;
  if (ffp->lock) {
    bsplist = DoLockFarComponents (sep, ffp);
    if (bsplist != NULL && ffp->fr != NULL) {
      CacheFarComponents (ffp, bsplist);
    }
  }

  if (ffp->nt != NULL) {
    if (SeqEntryFastaStream (sep, ffp->nt, flags, ffp->linelen, 0, 0,
                             TRUE, FALSE, ffp->master_style) < 0) {
      ffp->failed = TRUE;
    }
  }
  if (ffp->aa != NULL) {
    if (SeqEntryFastaStream (sep, ffp->aa, flags, ffp->linelen, 0, 0,
                             FALSE, TRUE, ffp->master_style) < 0) {
      ffp->failed = TRUE;
    }
  }
  if (ffp->ql != NULL) {
    VisitBioseqsInSep (sep, (Pointer) ffp, PrintQualScores);
  }

  bsplist = UnlockFarComponents (bsplist);
}

static void ProcessAutomaticRecord (
  CharPtr directory,
  CharPtr base,
  CharPtr suffix,
  FastaFlagPtr ffp
)

{
  Char  file /* [FILENAME_MAX] */ [PATH_MAX], path [PATH_MAX];

  if (ffp == NULL) return;

  if (base == NULL) {
    base = "";
  }
  if (suffix == NULL) {
    suffix = "";
  }

  StringNCpy_0 (path, directory, sizeof (path));
  sprintf (file, "%s%s", base, suffix);
  if (ffp->debugging) {
    printf ("FILENAME_MAX %d\n PATH_MAX %d\n", (int) FILENAME_MAX, (int) PATH_MAX);
    printf ("base '%s'\nsufx '%s'\nfile '%s' \n", base, suffix, file);
  }
  FileBuildPath (path, NULL, file);
  if (ffp->debugging) {
    printf ("path '%s' \n", path);
  }

  if (StringHasNoText (path)) return;

  ReadSequenceAsnFile (path, ffp->binary, ffp->compressed, (Pointer) ffp, FastaWrapper);
}

static void ProcessOneRecord (
  CharPtr directory,
  CharPtr base,
  CharPtr suffix,
  FastaFlagPtr ffp
)

{
  if (ffp == NULL) return;

  if (ffp->automatic) {
    ProcessAutomaticRecord (directory, base, suffix, ffp);
  } else if (ffp->batch) {
    ProcessMultipleRecord (directory, base, suffix, ffp);
  } else {
    ProcessSingleRecord (directory, base, suffix, ffp);
  }
}

static void ProcessOneSeqEntry (
  SeqEntryPtr sep,
  FastaFlagPtr ffp
)


{
  ValNodePtr     bsplist;
  StreamFlgType  flags = STREAM_CORRECT_INVAL;

  if (sep == NULL || ffp == NULL) return;

  BasicSeqEntryCleanup (sep);

  if (ffp->expand_gaps && ffp->use_dashes) {
    flags |= EXPAND_GAPS_TO_DASHES;
  } else if (ffp->expand_gaps) {
    flags |= STREAM_EXPAND_GAPS;
  } else if (ffp->use_dashes) {
    flags |= GAP_TO_SINGLE_DASH;
  }
  if (ffp->html_spans) {
    flags |= STREAM_HTML_SPANS;
  }
  if (ffp->extended_ids) {
    flags |= STREAM_ALL_FASTA_IDS;
  }

  bsplist = NULL;
  if (ffp->lock) {
    bsplist = DoLockFarComponents (sep, ffp);
    if (bsplist != NULL && ffp->fr != NULL) {
      CacheFarComponents (ffp, bsplist);
    }
  }

  if (ffp->nt != NULL) {
    if (SeqEntryFastaStream (sep, ffp->nt, flags, ffp->linelen, 0, 0,
                             TRUE, FALSE, ffp->master_style) < 0) {
      ffp->failed = TRUE;
    }
  }
  if (ffp->aa != NULL) {
    if (SeqEntryFastaStream (sep, ffp->aa, flags, ffp->linelen, 0, 0,
                             FALSE, TRUE, ffp->master_style) < 0) {
      ffp->failed = TRUE;
    }
  }
  if (ffp->ql != NULL) {
    VisitBioseqsInSep (sep, (Pointer) ffp, PrintQualScores);
  }

  bsplist = UnlockFarComponents (bsplist);
}

static void FileRecurse (
  CharPtr directory,
  CharPtr subdir,
  CharPtr suffix,
  Boolean dorecurse,
  FastaFlagPtr ffp
)

{
  Char        path [PATH_MAX];
  CharPtr     ptr, str;
  ValNodePtr  head, vnp;

  /* get list of all files in source directory */

  head = DirCatalog (directory);

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == 0) {
      if (StringHasNoText (subdir) || StringStr (directory, subdir) != NULL) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (! StringHasNoText (str)) {

          /* does filename have desired substring? */

          ptr = StringStr (str, suffix);
          if (ptr != NULL) {
            *ptr = '\0';

            /* process file that has desired suffix (usually .ent) */

            ProcessOneRecord (directory, str, suffix, ffp);
          }
        }
      }
    } else if (vnp->choice == 1 && dorecurse) {

      /* recurse into subdirectory */

      StringNCpy_0 (path, directory, sizeof (path));
      str = (CharPtr) vnp->data.ptrvalue;
      FileBuildPath (path, str, NULL);

      FileRecurse (path, str, suffix, dorecurse, ffp);
    }
  }

  /* clean up file list */

  ValNodeFreeData (head);
}

static SeqEntryPtr SeqEntryFromAccnOrGi (
  CharPtr str
)

{
  CharPtr      accn;
  BioseqPtr    bsp;
  Char         buf [64];
  Int4         flags = 0;
  Int2         retcode = 0;
  SeqEntryPtr  sep = NULL;
  SeqIdPtr     sip;
  CharPtr      tmp1 = NULL;
  CharPtr      tmp2 = NULL;
  long int     val;

  if (StringHasNoText (str)) return NULL;
  StringNCpy_0 (buf, str, sizeof (buf));
  TrimSpacesAroundString (buf);

  accn = buf;
  tmp1 = StringChr (accn, ',');
  if (tmp1 != NULL) {
    *tmp1 = '\0';
    tmp1++;
    tmp2 = StringChr (tmp1, ',');
    if (tmp2 != NULL) {
      *tmp2 = '\0';
      tmp2++;
      if (StringDoesHaveText (tmp2) && sscanf (tmp2, "%ld", &val) == 1) {
        flags = (Int4) val;
      }
    }
    if (StringDoesHaveText (tmp1) && sscanf (tmp1, "%ld", &val) == 1) {
      retcode = (Int2) val;
    }
  }

  sip = SeqIdFromPubSeqString (accn);
  sep = PubSeqSynchronousQueryId (sip, retcode, flags);

  if (sep != NULL) {
    bsp = BioseqFind (sip);
    if (bsp != NULL) {
      sep = SeqMgrGetSeqEntryForData ((Pointer) bsp);
    }
  }
  sip = SeqIdFree (sip);

  return sep;
}

/* Args structure contains command-line arguments */

#define p_argInputPath     0
#define i_argInputFile     1
#define o_argNtOutFile     2
#define v_argAaOutFile     3
#define q_argQlOutFile     4
#define x_argSuffix        5
#define u_argRecurse       6
#define m_argMaster        7
#define g_argExpandGaps    8
#define D_argUseDashes     9
#define E_argExtendedIDs  10
#define s_argGenomicQual  11
#define z_argZeroQualGap  12
#define a_argType         13
#define b_argBinary       14
#define c_argCompressed   15
#define r_argRemote       16
#define f_argFastaIdx     17
#define d_argBlastDB      18
#define k_argLocalFetch   19
#define l_argLockFar      20
#define h_argFarOutFile   21
#define e_argLineLength   22
#define T_argThreads      23
#define L_argLogFile      24
#define A_argAccession    25
#define H_argHtmlSpans    26
#define y_argDebugging    27

Args myargs [] = {
  {"Path to ASN.1 Files", NULL, NULL, NULL,
    TRUE, 'p', ARG_STRING, 0.0, 0, NULL},
  {"Single Input File", "stdin", NULL, NULL,
    TRUE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
  {"Nucleotide Output File Name", NULL, NULL, NULL,
    TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
  {"Protein Output File Name", NULL, NULL, NULL,
    TRUE, 'v', ARG_FILE_OUT, 0.0, 0, NULL},
  {"Quality Score Output File Name", NULL, NULL, NULL,
    TRUE, 'q', ARG_FILE_OUT, 0.0, 0, NULL},
  {"File Selection Substring", ".ent", NULL, NULL,
    TRUE, 'x', ARG_STRING, 0.0, 0, NULL},
  {"Recurse", "F", NULL, NULL,
    TRUE, 'u', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Master Style for Near Segmented Sequences", "F", NULL, NULL,
    TRUE, 'm', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Expand Delta Gaps into Ns", "F", NULL, NULL,
    TRUE, 'g', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Use Dash for Gap", "F", NULL, NULL,
    TRUE, 'D', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Extended Seq-IDs", "F", NULL, NULL,
    TRUE, 'E', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Far Genomic Contig for Quality Scores", "F", NULL, NULL,
    TRUE, 's', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Print Quality Score Gap as -1", "F", NULL, NULL,
    TRUE, 'z', ARG_BOOLEAN, 0.0, 0, NULL},
  {"ASN.1 Type (a Automatic, z Any, e Seq-entry, b Bioseq, s Bioseq-set, m Seq-submit, t Batch Processing)", "a", NULL, NULL,
    TRUE, 'a', ARG_STRING, 0.0, 0, NULL},
  {"Bioseq-set is Binary", "F", NULL, NULL,
    TRUE, 'b', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Bioseq-set is Compressed", "F", NULL, NULL,
    TRUE, 'c', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Remote Fetching from ID", "F", NULL, NULL,
    TRUE, 'r', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Path to Indexed FASTA Data", NULL, NULL, NULL,
    TRUE, 'f', ARG_STRING, 0.0, 0, NULL},
  {"Path to ReadDB Database", NULL, NULL, NULL,
    TRUE, 'd', ARG_STRING, 0.0, 0, NULL},
  {"Local Fetching", "F", NULL, NULL,
    TRUE, 'k', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Lock Components in Advance", "F", NULL, NULL,
    TRUE, 'l', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Far Component Cache Output File Name", NULL, NULL, NULL,
    TRUE, 'h', ARG_FILE_OUT, 0.0, 0, NULL},
  {"Line Length", "70", "10", "120",
    TRUE, 'e', ARG_INT, 0.0, 0, NULL},
  {"Use Threads", "F", NULL, NULL,
    TRUE, 'T', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Log File", NULL, NULL, NULL,
    TRUE, 'L', ARG_FILE_OUT, 0.0, 0, NULL},
  {"Accession to Fetch", NULL, NULL, NULL,
    TRUE, 'A', ARG_STRING, 0.0, 0, NULL},
  {"HTML Spans", "F", NULL, NULL,
    TRUE, 'H', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Debugging", "F", NULL, NULL,
    TRUE, 'y', ARG_BOOLEAN, 0.0, 0, NULL},
};

Int2 Main (void)

{
  Char           app [64], sfx [32];
  CharPtr        accn, base, blastdb, directory, fastaidx, ntout,
                 aaout, qlout, frout, logfile, ptr, str, suffix;
  Boolean        automatic, batch, binary, blast, compressed, dorecurse,
                 expandgaps, extendedids, fargenomicqual, fasta, htmlspans,
                 local, lock, masterstyle, qualgapzero, remote, usedashes,
                 usethreads, debugging;
  FastaFlagData  ffd;
  Int2           linelen, type = 0;
  time_t         run_time, start_time, stop_time;
  SeqEntryPtr    sep;

  /* standard setup */

  ErrSetFatalLevel (SEV_MAX);
  ErrClearOptFlags (EO_SHOW_USERSTR);
  ErrSetLogfile ("stderr", ELOG_APPEND);
  UseLocalAsnloadDataAndErrMsg ();
  ErrPathReset ();

  SOCK_SetupSSL(NcbiSetupGnuTls);

  if (! AllObjLoad ()) {
    Message (MSG_FATAL, "AllObjLoad failed");
    return 1;
  }
  if (! SubmitAsnLoad ()) {
    Message (MSG_FATAL, "SubmitAsnLoad failed");
    return 1;
  }
  if (! FeatDefSetLoad ()) {
    Message (MSG_FATAL, "FeatDefSetLoad failed");
    return 1;
  }
  if (! SeqCodeSetLoad ()) {
    Message (MSG_FATAL, "SeqCodeSetLoad failed");
    return 1;
  }
  if (! GeneticCodeTableLoad ()) {
    Message (MSG_FATAL, "GeneticCodeTableLoad failed");
    return 1;
  }

  /* process command line arguments */

  sprintf (app, "asn2fsa %s", ASN2FSA_APPLICATION);
  if (! GetArgs (app, sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  /* additional setup modifications */

  MemSet ((Pointer) &ffd, 0, sizeof (FastaFlagData));

  directory = (CharPtr) myargs [p_argInputPath].strvalue;
  suffix = (CharPtr) myargs [x_argSuffix].strvalue;
  base = (CharPtr) myargs [i_argInputFile].strvalue;
  accn = (CharPtr) myargs [A_argAccession].strvalue;
  dorecurse = (Boolean) myargs [u_argRecurse].intvalue;
  remote = (Boolean ) myargs [r_argRemote].intvalue;
  fastaidx = (CharPtr) myargs [f_argFastaIdx].strvalue;
  fasta = (Boolean) StringDoesHaveText (fastaidx);
  blastdb = (CharPtr) myargs [d_argBlastDB].strvalue;
  blast = (Boolean) StringDoesHaveText (blastdb);
  local = (Boolean) myargs [k_argLocalFetch].intvalue;
  lock = (Boolean) myargs [l_argLockFar].intvalue;
  linelen = (Int2) myargs [e_argLineLength].intvalue;
  usethreads = (Boolean) myargs [T_argThreads].intvalue;

  expandgaps = (Boolean) myargs [g_argExpandGaps].intvalue;
  htmlspans = (Boolean) myargs [H_argHtmlSpans].intvalue;
  usedashes = (Boolean) myargs [D_argUseDashes].intvalue;
  extendedids = (Boolean) myargs [E_argExtendedIDs].intvalue;
  masterstyle = (Boolean) myargs [m_argMaster].intvalue;
  fargenomicqual = (Boolean) myargs [s_argGenomicQual].intvalue;
  qualgapzero = (Boolean) myargs [z_argZeroQualGap].intvalue;
  automatic = FALSE;
  batch = FALSE;
  binary = (Boolean) myargs [b_argBinary].intvalue;
  compressed = (Boolean) myargs [c_argCompressed].intvalue;

  debugging = (Boolean) myargs [y_argDebugging].intvalue;

  str = myargs [a_argType].strvalue;
  if (StringICmp (str, "a") == 0) {
    type = 1;
    automatic = TRUE;
  } else if (StringICmp (str, "z") == 0) {
    type = 1;
  } else if (StringICmp (str, "e") == 0) {
    type = 2;
  } else if (StringICmp (str, "b") == 0) {
    type = 3;
  } else if (StringICmp (str, "s") == 0) {
    type = 4;
  } else if (StringICmp (str, "m") == 0) {
    type = 5;
  } else if (StringICmp (str, "t") == 0) {
    type = 1;
    batch = TRUE;
  } else {
    type = 1;
  }

  if ((binary || compressed) && (! batch)) {
    if (type == 1) {
      Message (MSG_FATAL, "-b or -c cannot be used without -t or -a");
      return 1;
    }
  }

  if (StringHasNoText (directory) && StringHasNoText (base)) {
    Message (MSG_FATAL, "Input path or input file must be specified");
    return 1;
  }

  ntout = (CharPtr) myargs [o_argNtOutFile].strvalue;
  aaout = (CharPtr) myargs [v_argAaOutFile].strvalue;
  qlout = (CharPtr) myargs [q_argQlOutFile].strvalue;
  frout = (CharPtr) myargs [h_argFarOutFile].strvalue;

  logfile = (CharPtr) myargs [L_argLogFile].strvalue;

  /* default to stdout for nucleotide output if nothing specified */

  if (StringHasNoText (ntout) &&
      StringHasNoText (aaout) &&
      StringHasNoText (qlout)) {
    ntout = "stdout";
  }

  start_time = GetSecs ();

  /* populate parameter structure */

  ffd.expand_gaps = expandgaps;
  ffd.html_spans = htmlspans;
  ffd.use_dashes = usedashes;
  ffd.extended_ids = extendedids;
  ffd.master_style = masterstyle;
  ffd.far_genomic_qual = fargenomicqual;
  ffd.qual_gap_is_zero = (Boolean) (! qualgapzero);
  ffd.automatic = automatic;
  ffd.batch = batch;
  ffd.binary = binary;
  ffd.compressed = compressed;
  ffd.lock = lock;
  ffd.useThreads = usethreads;
  ffd.type = type;
  ffd.linelen = linelen;
  ffd.failed = FALSE;
  ffd.nt = NULL;
  ffd.aa = NULL;
  ffd.ql = NULL;
  ffd.fr = NULL;
  ffd.logfp = NULL;
  ffd.debugging = debugging;

  if (! StringHasNoText (ntout)) {
    ffd.nt = FileOpen (ntout, "w");
    if (ffd.nt == NULL) {
      Message (MSG_FATAL, "Unable to open nucleotide output file");
      return 1;
    }
  }

  if (! StringHasNoText (aaout)) {
    ffd.aa = FileOpen (aaout, "w");
    if (ffd.aa == NULL) {
      Message (MSG_FATAL, "Unable to open protein output file");
      return 1;
    }
  }

  if (! StringHasNoText (qlout)) {
    ffd.ql = FileOpen (qlout, "w");
    if (ffd.ql == NULL) {
      Message (MSG_FATAL, "Unable to open quality score output file");
      return 1;
    }
  }

  if (! StringHasNoText (frout)) {
    ffd.fr = FileOpen (frout, "w");
    if (ffd.fr == NULL) {
      Message (MSG_FATAL, "Unable to open far component cache output file");
      return 1;
    }
    ffd.lock = TRUE;
  }

  if (! StringHasNoText (logfile)) {
    ffd.logfp = FileOpen (logfile, "w");
    if (ffd.logfp == NULL) {
      Message (MSG_FATAL, "Unable to open log file");
      return 1;
    }
  }

  /* register fetch functions */

  if (remote) {
#ifdef INTERNAL_NCBI_ASN2FSA
    if (! PUBSEQBioseqFetchEnable ("asn2fsa", FALSE)) {
      Message (MSG_POSTERR, "PUBSEQBioseqFetchEnable failed");
      return 1;
    }
    ffd.usePUBSEQ = TRUE;
    ffd.useThreads = FALSE;
#else
    PubSeqFetchEnable ();
#endif
  }

  if (blast) {
    ptr = StringRChr (blastdb, DIRDELIMCHR);
    if (ptr != NULL) {
      *ptr = '\0';
      ptr++;
      TransientSetAppParam ("NCBI", "BLAST", "BLASTDB", blastdb);
      if (StringDoesHaveText (ptr)) {
        ReadDBBioseqFetchEnable ("asn2fsa", ptr, TRUE, FALSE);
        ffd.blastdbname = ptr;
        ffd.useBLAST = TRUE;
      } else {
        ReadDBBioseqFetchEnable ("asn2fsa", "nr", TRUE, FALSE);
        ffd.blastdbname = "nr";
        ffd.useBLAST = TRUE;
      }
    } else {
      ReadDBBioseqFetchEnable ("asn2fsa", blastdb, TRUE, FALSE);
      ffd.blastdbname = blastdb;
      ffd.useBLAST = TRUE;
    }
  }

  if (fasta) {
    AltIndexedFastaLibFetchEnable (fastaidx);
  }

  if (local) {
    LocalSeqFetchInit (FALSE);
  }

  /* recurse through all files within source directory or subdirectories */

  if (StringDoesHaveText (accn)) {

    if (remote) {
      sep = SeqEntryFromAccnOrGi (accn);
      if (sep != NULL) {
        ProcessOneSeqEntry (sep, &ffd);
        SeqEntryFree (sep);
      }
    }

  } else if (StringDoesHaveText (directory)) {

    FileRecurse (directory, NULL, suffix, dorecurse, &ffd);

  } else if (StringDoesHaveText (base)) {

    ptr = StringRChr (base, '.');
    sfx[0] = '\0';
    /* check for file without suffix but path includes dotted URL name */
    if (ptr != NULL && StringChr (ptr, '/') == NULL) {
      StringNCpy_0 (sfx, ptr, sizeof (sfx));
      *ptr = '\0';
    }
    ProcessOneRecord (directory, base, sfx, &ffd);
  }

  if (ffd.nt != NULL) {
    FileClose (ffd.nt);
  }
  if (ffd.aa != NULL) {
    FileClose (ffd.aa);
  }
  if (ffd.ql != NULL) {
    FileClose (ffd.ql);
  }
  if (ffd.fr != NULL) {
    FileClose (ffd.fr);
    CreateFastaIndex (frout);
  }

  stop_time = GetSecs ();
  run_time = stop_time - start_time;

  if (ffd.logfp != NULL) {
    fprintf (ffd.logfp, "Finished in %ld seconds\n", (long) run_time);
    FileClose (ffd.logfp);
  }

  /* close fetch functions */

  if (local) {
    LocalSeqFetchDisable ();
  }

  if (fasta) {
    AltIndexedFastaLibFetchDisable ();
  }

  if (blast) {
    ReadDBBioseqFetchDisable ();
  }

  if (remote) {
#ifdef INTERNAL_NCBI_ASN2FSA
    PUBSEQBioseqFetchDisable ();
#else
    PubSeqFetchDisable ();
#endif
  }

  if (ffd.failed) {
    return 1;
  }

  return 0;
}

