/*   asnval.c
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
* File Name:  asnval.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   11/3/04
*
* $Revision: 1.185 $
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
#include <gather.h>
#include <explore.h>
#include <lsqfetch.h>
#include <valid.h>
#include <validerr.h>
#include <pmfapi.h>
#include <ent2api.h>
#include <gbftdef.h>
#include <objtax3.h>
#include <tax3api.h>
#ifdef INTERNAL_NCBI_ASN2VAL
#include <accpubseq.h>
#endif
#include <connect/ncbi_gnutls.h>

#define ASNVAL_APP_VER "15.4"

CharPtr ASNVAL_APPLICATION = ASNVAL_APP_VER;

typedef struct valflags {
  Int2     severity;
  Int2     lowCutoff;
  Int2     highCutoff;
  CharPtr  errcode;
  Boolean  validateAlignments;
  Boolean  alignFindRemoteBsp;
  Boolean  doSeqHistAssembly;
  Boolean  farIDsInAlignments;
  Boolean  alwaysRequireIsoJTA;
  Boolean  farFetchCDSproducts;
  Boolean  farFetchMRNAproducts;
  Boolean  locusTagGeneralMatch;
  Boolean  validateIDSet;
  Boolean  seqSubmitParent;
  Boolean  ignoreExceptions;
  Boolean  validateExons;
  Boolean  inferenceAccnCheck;
  Boolean  testLatLonSubregion;
  Boolean  strictLatLonCountry;
  Boolean  rubiscoTest;
  Boolean  disableSuppression;
  Boolean  genomeSubmission;
  Boolean  debugTestDuJour;
  Boolean  indexerVersion;
  Boolean  automatic;
  Boolean  catenated;
  Boolean  batch;
  Boolean  binary;
  Boolean  compressed;
  Boolean  lock;
  Boolean  useThreads;
  Boolean  usePUBSEQ;
  Boolean  validateBarcode;
  Boolean  taxfetch;
  Int2     verbosity;
  Int2     type;
  Int4     skipcount;
  Int4     maxcount;
  CharPtr  outpath;
  FILE     *outfp;
  FILE     *logfp;
  Int4     num_errors;
  Int4     fatal_errors;
  Boolean  has_errors;
  Boolean  io_failure;
  Char     longest [64];
  time_t   worsttime;
  Int4     numrecords;
  Char     path [PATH_MAX];
} ValFlagData, PNTR ValFlagPtr;

#ifdef INTERNAL_NCBI_ASN2VAL
static CharPtr dirsubfetchproc = "DirSubBioseqFetch";

static CharPtr dirsubfetchcmd = NULL;

extern Pointer ReadFromDirSub (CharPtr accn, Uint2Ptr datatype, Uint2Ptr entityID);
extern Pointer ReadFromDirSub (CharPtr accn, Uint2Ptr datatype, Uint2Ptr entityID)

{
  Char     cmmd [256];
  Pointer  dataptr;
  FILE*    fp;
  Char     path [PATH_MAX];

  if (datatype != NULL) {
    *datatype = 0;
  }
  if (entityID != NULL) {
    *entityID = 0;
  }
  if (StringHasNoText (accn)) return NULL;

  if (dirsubfetchcmd == NULL) {
    if (GetAppParam ("SEQUIN", "DIRSUB", "FETCHSCRIPT", NULL, cmmd, sizeof (cmmd))) {
        dirsubfetchcmd = StringSaveNoNull (cmmd);
    }
  }
  if (dirsubfetchcmd == NULL) return NULL;

  TmpNam (path);

#ifdef OS_UNIX
  sprintf (cmmd, "csh %s %s > %s", dirsubfetchcmd, accn, path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, "%s %s -o %s", dirsubfetchcmd, accn, path);
  system (cmmd);
#endif

  fp = FileOpen (path, "r");
  if (fp == NULL) {
    FileRemove (path);
    return NULL;
  }
  dataptr = ReadAsnFastaOrFlatFile (fp, datatype, entityID, FALSE, FALSE, TRUE, FALSE);
  FileClose (fp);
  FileRemove (path);
  return dataptr;
}


static Int2 LIBCALLBACK DirSubBioseqFetchFunc (Pointer data)

{
  BioseqPtr         bsp;
  Char              cmmd [256];
  Pointer           dataptr;
  Uint2             datatype;
  Uint2             entityID;
  FILE*             fp;
  OMProcControlPtr  ompcp;
  ObjMgrProcPtr     ompp;
  Char              path [PATH_MAX];
  SeqEntryPtr       sep = NULL;
  SeqIdPtr          sip;
  TextSeqIdPtr      tsip;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL) return OM_MSG_RET_ERROR;
  ompp = ompcp->proc;
  if (ompp == NULL) return OM_MSG_RET_ERROR;
  sip = (SeqIdPtr) ompcp->input_data;
  if (sip == NULL) return OM_MSG_RET_ERROR;

  if (sip->choice != SEQID_GENBANK) return OM_MSG_RET_ERROR;
  tsip = (TextSeqIdPtr) sip->data.ptrvalue;
  if (tsip == NULL || StringHasNoText (tsip->accession)) return OM_MSG_RET_ERROR;

  if (dirsubfetchcmd == NULL) {
    if (GetAppParam ("SEQUIN", "DIRSUB", "FETCHSCRIPT", NULL, cmmd, sizeof (cmmd))) {
        dirsubfetchcmd = StringSaveNoNull (cmmd);
    }
  }
  if (dirsubfetchcmd == NULL) return OM_MSG_RET_ERROR;

  TmpNam (path);

#ifdef OS_UNIX
  sprintf (cmmd, "csh %s %s > %s", dirsubfetchcmd, tsip->accession, path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, "%s %s -o %s", dirsubfetchcmd, tsip->accession, path);
  system (cmmd);
#endif

  fp = FileOpen (path, "r");
  if (fp == NULL) {
    FileRemove (path);
    return OM_MSG_RET_ERROR;
  }
  dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, &entityID, FALSE, FALSE, TRUE, FALSE);
  FileClose (fp);
  FileRemove (path);

  if (dataptr == NULL) return OM_MSG_RET_OK;

  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) return OM_MSG_RET_ERROR;
  bsp = BioseqFindInSeqEntry (sip, sep);
  ompcp->output_data = (Pointer) bsp;
  ompcp->output_entityID = ObjMgrGetEntityIDForChoice (sep);
  return OM_MSG_RET_DONE;
}

static Boolean DirSubFetchEnable (void)

{
  ObjMgrProcLoad (OMPROC_FETCH, dirsubfetchproc, dirsubfetchproc,
                  OBJ_SEQID, 0, OBJ_BIOSEQ, 0, NULL,
                  DirSubBioseqFetchFunc, PROC_PRIORITY_DEFAULT);
  return TRUE;
}

static CharPtr smartfetchproc = "SmartBioseqFetch";

static CharPtr smartfetchcmd = NULL;

extern Pointer ReadFromSmart (CharPtr accn, Uint2Ptr datatype, Uint2Ptr entityID);
extern Pointer ReadFromSmart (CharPtr accn, Uint2Ptr datatype, Uint2Ptr entityID)

{
  Char     cmmd [256];
  Pointer  dataptr;
  FILE*    fp;
  Char     path [PATH_MAX];

  if (datatype != NULL) {
    *datatype = 0;
  }
  if (entityID != NULL) {
    *entityID = 0;
  }
  if (StringHasNoText (accn)) return NULL;

  if (smartfetchcmd == NULL) {
    if (GetAppParam ("SEQUIN", "SMART", "FETCHSCRIPT", NULL, cmmd, sizeof (cmmd))) {
        smartfetchcmd = StringSaveNoNull (cmmd);
    }
  }
  if (smartfetchcmd == NULL) return NULL;

  TmpNam (path);

#ifdef OS_UNIX
  sprintf (cmmd, "csh %s %s > %s", smartfetchcmd, accn, path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, "%s %s -o %s", smartfetchcmd, accn, path);
  system (cmmd);
#endif

  fp = FileOpen (path, "r");
  if (fp == NULL) {
    FileRemove (path);
    return NULL;
  }
  dataptr = ReadAsnFastaOrFlatFile (fp, datatype, entityID, FALSE, FALSE, TRUE, FALSE);
  FileClose (fp);
  FileRemove (path);
  return dataptr;
}


static Int2 LIBCALLBACK SmartBioseqFetchFunc (Pointer data)

{
  BioseqPtr         bsp;
  Char              cmmd [256];
  Pointer           dataptr;
  Uint2             datatype;
  Uint2             entityID;
  FILE*             fp;
  OMProcControlPtr  ompcp;
  ObjMgrProcPtr     ompp;
  Char              path [PATH_MAX];
  SeqEntryPtr       sep = NULL;
  SeqIdPtr          sip;
  TextSeqIdPtr      tsip;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL) return OM_MSG_RET_ERROR;
  ompp = ompcp->proc;
  if (ompp == NULL) return OM_MSG_RET_ERROR;
  sip = (SeqIdPtr) ompcp->input_data;
  if (sip == NULL) return OM_MSG_RET_ERROR;

  if (sip->choice != SEQID_GENBANK) return OM_MSG_RET_ERROR;
  tsip = (TextSeqIdPtr) sip->data.ptrvalue;
  if (tsip == NULL || StringHasNoText (tsip->accession)) return OM_MSG_RET_ERROR;

  if (smartfetchcmd == NULL) {
    if (GetAppParam ("SEQUIN", "SMART", "FETCHSCRIPT", NULL, cmmd, sizeof (cmmd))) {
        smartfetchcmd = StringSaveNoNull (cmmd);
    }
  }
  if (smartfetchcmd == NULL) return OM_MSG_RET_ERROR;

  TmpNam (path);

#ifdef OS_UNIX
  sprintf (cmmd, "csh %s %s > %s", smartfetchcmd, tsip->accession, path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, "%s %s -o %s", smartfetchcmd, tsip->accession, path);
  system (cmmd);
#endif

  fp = FileOpen (path, "r");
  if (fp == NULL) {
    FileRemove (path);
    return OM_MSG_RET_ERROR;
  }
  dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, &entityID, FALSE, FALSE, TRUE, FALSE);
  FileClose (fp);
  FileRemove (path);

  if (dataptr == NULL) return OM_MSG_RET_OK;

  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) return OM_MSG_RET_ERROR;
  bsp = BioseqFindInSeqEntry (sip, sep);
  ompcp->output_data = (Pointer) bsp;
  ompcp->output_entityID = ObjMgrGetEntityIDForChoice (sep);
  return OM_MSG_RET_DONE;
}

static Boolean SmartFetchEnable (void)

{
  ObjMgrProcLoad (OMPROC_FETCH, smartfetchproc, smartfetchproc,
                  OBJ_SEQID, 0, OBJ_BIOSEQ, 0, NULL,
                  SmartBioseqFetchFunc, PROC_PRIORITY_DEFAULT);
  return TRUE;
}

static CharPtr tpasmartfetchproc = "TPASmartBioseqFetch";

static CharPtr tpasmartfetchcmd = NULL;

extern Pointer ReadFromTPASmart (CharPtr accn, Uint2Ptr datatype, Uint2Ptr entityID);
extern Pointer ReadFromTPASmart (CharPtr accn, Uint2Ptr datatype, Uint2Ptr entityID)

{
  Char     cmmd [256];
  Pointer  dataptr;
  FILE*    fp;
  Char     path [PATH_MAX];

  if (datatype != NULL) {
    *datatype = 0;
  }
  if (entityID != NULL) {
    *entityID = 0;
  }
  if (StringHasNoText (accn)) return NULL;

  if (tpasmartfetchcmd == NULL) {
    if (GetAppParam ("SEQUIN", "TPASMART", "FETCHSCRIPT", NULL, cmmd, sizeof (cmmd))) {
        tpasmartfetchcmd = StringSaveNoNull (cmmd);
    }
  }
  if (tpasmartfetchcmd == NULL) return NULL;

  TmpNam (path);

#ifdef OS_UNIX
  sprintf (cmmd, "csh %s %s > %s", tpasmartfetchcmd, accn, path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, "%s %s -o %s", tpasmartfetchcmd, accn, path);
  system (cmmd);
#endif

  fp = FileOpen (path, "r");
  if (fp == NULL) {
    FileRemove (path);
    return NULL;
  }
  dataptr = ReadAsnFastaOrFlatFile (fp, datatype, entityID, FALSE, FALSE, TRUE, FALSE);
  FileClose (fp);
  FileRemove (path);
  return dataptr;
}


static Int2 LIBCALLBACK TPASmartBioseqFetchFunc (Pointer data)

{
  BioseqPtr         bsp;
  Char              cmmd [256];
  Pointer           dataptr;
  Uint2             datatype;
  Uint2             entityID;
  FILE*             fp;
  OMProcControlPtr  ompcp;
  ObjMgrProcPtr     ompp;
  Char              path [PATH_MAX];
  SeqEntryPtr       sep = NULL;
  SeqIdPtr          sip;
  TextSeqIdPtr      tsip;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL) return OM_MSG_RET_ERROR;
  ompp = ompcp->proc;
  if (ompp == NULL) return OM_MSG_RET_ERROR;
  sip = (SeqIdPtr) ompcp->input_data;
  if (sip == NULL) return OM_MSG_RET_ERROR;

  if (sip->choice != SEQID_TPG) return OM_MSG_RET_ERROR;
  tsip = (TextSeqIdPtr) sip->data.ptrvalue;
  if (tsip == NULL || StringHasNoText (tsip->accession)) return OM_MSG_RET_ERROR;

  if (tpasmartfetchcmd == NULL) {
    if (GetAppParam ("SEQUIN", "TPASMART", "FETCHSCRIPT", NULL, cmmd, sizeof (cmmd))) {
        tpasmartfetchcmd = StringSaveNoNull (cmmd);
    }
  }
  if (tpasmartfetchcmd == NULL) return OM_MSG_RET_ERROR;

  TmpNam (path);

#ifdef OS_UNIX
  sprintf (cmmd, "csh %s %s > %s", tpasmartfetchcmd, tsip->accession, path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, "%s %s -o %s", tpasmartfetchcmd, tsip->accession, path);
  system (cmmd);
#endif

  fp = FileOpen (path, "r");
  if (fp == NULL) {
    FileRemove (path);
    return OM_MSG_RET_ERROR;
  }
  dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, &entityID, FALSE, FALSE, TRUE, FALSE);
  FileClose (fp);
  FileRemove (path);

  if (dataptr == NULL) return OM_MSG_RET_OK;

  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) return OM_MSG_RET_ERROR;
  bsp = BioseqFindInSeqEntry (sip, sep);
  ompcp->output_data = (Pointer) bsp;
  ompcp->output_entityID = ObjMgrGetEntityIDForChoice (sep);
  return OM_MSG_RET_DONE;
}

static Boolean TPASmartFetchEnable (void)

{
  ObjMgrProcLoad (OMPROC_FETCH, tpasmartfetchproc, tpasmartfetchproc,
                  OBJ_SEQID, 0, OBJ_BIOSEQ, 0, NULL,
                  TPASmartBioseqFetchFunc, PROC_PRIORITY_DEFAULT);
  return TRUE;
}
#endif

static void LookForBigFarSeqs (
  BioseqPtr bsp,
  Pointer userdata
)

{
  Int4         count = 0;
  DeltaSeqPtr  dsp;
  Boolean      is_ddbj = FALSE;
  SeqIdPtr     sip;
  BoolPtr      toomanyfarP;

  if (bsp == NULL || userdata == NULL) return;

  if (bsp->repr != Seq_repr_delta) return;
  if (bsp->seq_ext_type != 4) return;

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_DDBJ) {
      is_ddbj = TRUE;
    }
  }

  if (! is_ddbj) return;

  for (dsp = (DeltaSeqPtr) bsp->seq_ext; dsp != NULL; dsp = dsp->next) {
    if (dsp->choice == 1) {
      count++;
    }
  }

  if (count > 10000) {
    toomanyfarP = (BoolPtr) userdata;
    *toomanyfarP = TRUE;
  }
}

static Boolean TooManyFarComponents (
  SeqEntryPtr sep
)

{
  Boolean  toomanyfar = FALSE;

  if (sep == NULL) return FALSE;

  VisitBioseqsInSep (sep, (Pointer) &toomanyfar, LookForBigFarSeqs);

  return toomanyfar;
}

static ValNodePtr DoLockFarComponents (
  SeqEntryPtr sep,
  ValFlagPtr vfp
)

{
  Boolean     farFetch;
  ValNodePtr  rsult;
  time_t      start_time, stop_time;

  start_time = GetSecs ();

#ifdef INTERNAL_NCBI_ASN2VAL
  if (vfp->useThreads) {
    Message (MSG_POST, "Threads will not be used in this executable");
    vfp->useThreads = FALSE;;
  }
#endif

  farFetch = (Boolean) (vfp->farFetchCDSproducts);

  if (NlmThreadsAvailable () && vfp->useThreads) {
    rsult = AdvcLockFarComponents (sep, TRUE, farFetch, farFetch, NULL, TRUE);
  } else if (vfp->useThreads) {
    Message (MSG_POST, "Threads not available in this executable");
    rsult = AdvcLockFarComponents (sep, TRUE, farFetch, farFetch, NULL, FALSE);
  } else {
    rsult = AdvcLockFarComponents (sep, TRUE, farFetch, farFetch, NULL, FALSE);
  }

  stop_time = GetSecs ();

  return rsult;
}

static CharPtr severityLabel [] = {
  "NONE", "INFO", "WARNING", "ERROR", "REJECT", "FATAL", "MAX", NULL
};

static CharPtr compatSeverityLabel [] = {
  "NONE", "NOTE: valid", "WARNING: valid", "ERROR: valid", "REJECT: valid", "FATAL: valid", "MAX", NULL
};

typedef struct vcdaa {
  FILE        *ofp;
  Int2        verbosity;
  Int2        lowCutoff;
  Int2        highCutoff;
  CharPtr     errcode;
  ValFlagPtr  vfp;
} VCData, PNTR VCPtr;

static void XmlEncode (CharPtr dst, CharPtr src)

{
  Char  ch;

  if (dst == NULL || src == NULL) return;

  ch = *src;
  while (ch != '\0') {
    if (ch == '<') {
      *dst = '&';
      dst++;
      *dst = 'l';
      dst++;
      *dst = 't';
      dst++;
      *dst = ';';
      dst++;
    } else if (ch == '>') {
      *dst = '&';
      dst++;
      *dst = 'g';
      dst++;
      *dst = 't';
      dst++;
      *dst = ';';
      dst++;
    } else {
      *dst = ch;
      dst++;
    }
    src++;
    ch = *src;
  }
  *dst = '\0';
}


static CharPtr GetXmlHeaderText (ErrSev cutoff)
{
  CharPtr         xml_header = NULL;
  CharPtr         xml_4_fmt = "asnval version=\"%s\" severity_cutoff=\"%s\"";

  xml_header = (CharPtr) MemNew (sizeof (Char) * (10 + StringLen (xml_4_fmt) +
                StringLen (ASNVAL_APPLICATION) + StringLen (severityLabel[cutoff])));
  sprintf (xml_header, xml_4_fmt, ASNVAL_APPLICATION, severityLabel[cutoff]);
  return xml_header;
}


static void LIBCALLBACK ValidCallback (
  ErrSev severity,
  int errcode,
  int subcode,
  Uint2 entityID,
  Uint2 itemtype,
  Uint4 itemID,
  CharPtr accession,
  CharPtr seqid,
  CharPtr featureID,
  CharPtr message,
  CharPtr objtype,
  CharPtr label,
  CharPtr context,
  CharPtr location,
  CharPtr product,
  Pointer userdata
)

{
  Char        buf [256];
  CharPtr     catname, errname, seqidloc = NULL, urlmssg = NULL, urlloc = NULL;
  ErrSev      cutoff;
  FILE        *fp;
  Int4        gi = 0;
  size_t      len;
  SeqIdPtr    sip;
  VCPtr       vcp;
  ValFlagPtr  vfp;
  CharPtr     xml_header;

  vcp = (VCPtr) userdata;
  if (vcp == NULL) return;
  fp = vcp->ofp;
  if (fp == NULL) return;
  vfp = vcp->vfp;
  if (vfp == NULL) return;

  if (severity < SEV_NONE || severity > SEV_MAX) {
    severity = SEV_MAX;
  }

  if (severity < vcp->lowCutoff || severity > vcp->highCutoff) return;

  if (errcode == 0 && subcode == 0) {
    if (vcp->verbosity == 4) {
      catname = "Progress";
      errname = "Progress";
    } else {
      return;
    }
  } else {
    catname = GetValidCategoryName (errcode);
    errname = GetValidErrorName (errcode, subcode);

    if (catname == NULL) {
      catname = "?";
    }
    if (errname == NULL) {
      errname = "?";
    }

    if (StringDoesHaveText (vcp->errcode)) {
      if (StringICmp (vcp->errcode, errname) != 0) return;
    }
  }

  if (accession == NULL) {
    accession = "";
  }
  if (seqid == NULL) {
    seqid = "";
  }
  if (featureID == NULL) {
    featureID = "";
  }
  if (message == NULL) {
    message = "";
  }
  if (objtype == NULL) {
    objtype = "";
  }
  if (label == NULL) {
    label = "";
  }

  if (vcp->verbosity == 1 || vcp->verbosity == 6) {

    fprintf (fp, "%s [%s.%s] %s %s: %s",
             compatSeverityLabel [severity],
             catname, errname, message, objtype, label);
    if (StringDoesHaveText (featureID)) {
      fprintf (fp, " <%s>", featureID);
    }
    if (location != NULL) {
      fprintf (fp, " %s", location);
    }
    if (context != NULL) {
      fprintf (fp, " %s", context);
    }
    if (product != NULL) {
      fprintf (fp, " -> %s", product);
    }
    fprintf (fp, "\n");

  } else if (vcp->verbosity == 2) {

    StringCpy (buf, accession);
    StringCat (buf, "                    ");
    buf [15] = '\0';

    StringCat (buf, severityLabel [severity]);
    StringCat (buf, "                    ");
    buf [30] = '\0';

    StringCat (buf, catname);
    StringCat (buf, "_");
    StringCat (buf, errname);

    fprintf (fp, "%s\n", buf);

  } else if (vcp->verbosity == 3) {

    fprintf (fp, "%s\t%s\t%s_%s\n",
             accession, severityLabel [severity],
             catname, errname);

  } else if (vcp->verbosity == 4) {

    if (! vfp->has_errors) {
      cutoff = (ErrSev) vcp->lowCutoff;
      if (cutoff < SEV_NONE || cutoff > SEV_MAX) {
        cutoff = SEV_MAX;
      }

      xml_header = GetXmlHeaderText (cutoff);
      fprintf (fp, "<%s>\n", xml_header);
      xml_header = MemFree (xml_header);
    }

    if (location == NULL) {
      location = "";
    }
    len = StringLen (message);
    if (len > 0) {
      urlmssg = MemNew (len * 3 + 2);
      if (urlmssg != NULL) {
        XmlEncode (urlmssg, message);
        seqidloc = MemNew (StringLen (seqid) * 3 + 2);
        if (seqidloc != NULL) {
          XmlEncode (seqidloc, seqid);
          if (StringDoesHaveText (location)) {
            urlloc = MemNew (StringLen (location) * 3 + 2);
            if (urlloc != NULL) {
              XmlEncode (urlloc, location);
              if (StringDoesHaveText (featureID)) {
                fprintf (fp, "  <message severity=\"%s\" seq-id=\"%s\" feat-id=\"%s\" interval=\"%s\" code=\"%s_%s\">%s</message>\n",
                         severityLabel [severity], seqidloc, featureID, urlloc, catname, errname, urlmssg);
              } else {
                fprintf (fp, "  <message severity=\"%s\" seq-id=\"%s\" interval=\"%s\" code=\"%s_%s\">%s</message>\n",
                         severityLabel [severity], seqidloc, urlloc, catname, errname, urlmssg);
              }
            }
            MemFree (urlloc);
          } else {
            fprintf (fp, "  <message severity=\"%s\" seq-id=\"%s\" code=\"%s_%s\">%s</message>\n",
                     severityLabel [severity], seqidloc, catname, errname, urlmssg);
          }
          MemFree (seqidloc);
        }
        MemFree (urlmssg);
      }
    }
    fflush (fp);

  } else if (vcp->verbosity == 5) {

    sip = SeqIdFromAccessionDotVersion (accession);
    gi = GetGIForSeqId (sip);
    SeqIdFree (sip);

    fprintf (fp, "%s\t%ld\t%s\t%s_%s\n",
             accession, (long) gi, severityLabel [severity],
             catname, errname);
  }

  vfp->has_errors = TRUE;
}

typedef struct tax3val {
  Uint2       entityID;
  Uint4       itemID;
  Uint2       itemtype;
  Uint1       organelle;
  OrgRefPtr   orp;
  BioseqPtr   bsp;
  SeqFeatPtr  sfp;
} TaxVal, PNTR TaxValPtr;

typedef struct tax3lst {
  ValNodePtr  head;
  ValNodePtr  tail;
} TaxLst, PNTR TaxLstPtr;

static void RecordSrc (Uint2 entityID, Uint4 itemID, Uint2 itemtype, OrgRefPtr orp,
                       Uint1 organelle, TaxLstPtr tlp, SeqDescrPtr sdp, SeqFeatPtr sfp)

{
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  ObjValNodePtr  ovp;
  SeqEntryPtr    sep;
  TaxValPtr      tvp;
  ValNodePtr     vnp;

  if (orp == NULL || tlp == NULL) return;

  tvp = (TaxValPtr) MemNew (sizeof (TaxVal));
  if (tvp == NULL) return;

  vnp = ValNodeNew (tlp->tail);
  if (vnp == NULL) return;

  if (tlp->head == NULL) {
    tlp->head = vnp;
  }
  tlp->tail = vnp;

  tvp->entityID = entityID;
  tvp->itemID = itemID;
  tvp->itemtype = itemtype;
  tvp->organelle = organelle;
  tvp->orp = orp;
  if (sdp != NULL && sdp->extended != 0) {
    ovp = (ObjValNodePtr) sdp;
    if (ovp->idx.parenttype == OBJ_BIOSEQ) {
      bsp = (BioseqPtr) ovp->idx.parentptr;
      if (bsp != NULL) {
        tvp->bsp = bsp;
      }
    } else if (ovp->idx.parenttype == OBJ_BIOSEQSET) {
      bssp = (BioseqSetPtr) ovp->idx.parentptr;
      if (bssp != NULL) {
        sep = bssp->seqentry;
        if (sep != NULL) {
          sep = FindNthBioseq (sep, 1);
          if (sep != NULL) {
            bsp = (BioseqPtr) sep->data.ptrvalue;
            if (bsp != NULL) {
              tvp->bsp = bsp;
            }
          }
        }
      }
    }
  } else if (sfp != NULL) {
    tvp->sfp = sfp;
  }

  vnp->data.ptrvalue = tvp;
}

static void GetSrcDesc (SeqDescrPtr sdp, Pointer userdata)

{
  BioSourcePtr   biop;
  ObjValNodePtr  ovp;
  TaxLstPtr      tlp;

  if (sdp == NULL || sdp->choice != Seq_descr_source) return;
  tlp = (TaxLstPtr) userdata;

  biop = (BioSourcePtr) sdp->data.ptrvalue;
  if (biop == NULL) return;

  if (sdp->extended != 0) {
    ovp = (ObjValNodePtr) sdp;
    RecordSrc (ovp->idx.entityID, ovp->idx.itemID, OBJ_SEQDESC, biop->org, biop->genome, tlp, sdp, NULL);
  }
}

static void GetSrcFeat (SeqFeatPtr sfp, Pointer userdata)

{
  BioSourcePtr  biop;
  TaxLstPtr     tlp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_BIOSRC) return;
  tlp = (TaxLstPtr) userdata;

  biop = (BioSourcePtr) sfp->data.value.ptrvalue;
  if (biop == NULL) return;

  RecordSrc (sfp->idx.entityID, sfp->idx.itemID, OBJ_SEQFEAT, biop->org, biop->genome, tlp, NULL, sfp);
}

NLM_EXTERN void CDECL  ValidErr VPROTO((ValidStructPtr vsp, int severity, int code1, int code2, const char *fmt, ...));


static void ReportOneBadSpecificHost (ValNodePtr vnp, ValidStructPtr vsp, ErrSev sev, int code1, int code2, CharPtr msg_fmt)
{
  ObjValNodePtr ovp;
  BioSourcePtr  biop;
  OrgModPtr     mod;

  if (vnp == NULL || vsp == NULL || StringHasNoText (msg_fmt)) return;

  vsp->sfp = NULL;
  vsp->descr = NULL;
  vsp->bsp = NULL;
  vsp->bssp = NULL;
  biop = NULL;
  mod = NULL;

  if (vnp->choice == OBJ_SEQFEAT)
  {
    vsp->sfp = (SeqFeatPtr) vnp->data.ptrvalue;
    vsp->gcp->entityID = vsp->sfp->idx.entityID;
    vsp->gcp->itemID = vsp->sfp->idx.itemID;
    vsp->gcp->thistype = OBJ_SEQFEAT;
    if (vsp->sfp->idx.parenttype == OBJ_BIOSEQ)
    {
      vsp->bsp = vsp->sfp->idx.parentptr;        
    }
    else if (vsp->sfp->idx.parenttype == OBJ_BIOSEQSET)
    {
      vsp->bssp = vsp->sfp->idx.parentptr;        
    }
    biop = (BioSourcePtr) vsp->sfp->data.value.ptrvalue;
  } 
  else if (vnp->choice == OBJ_SEQDESC)
  {
    vsp->descr = (SeqDescrPtr) vnp->data.ptrvalue;
    if (vsp->descr != NULL && vsp->descr->extended != 0) 
    {
      ovp = (ObjValNodePtr) vsp->descr;
      vsp->gcp->entityID = ovp->idx.entityID;
      vsp->gcp->itemID = ovp->idx.itemID;
      vsp->gcp->thistype = OBJ_SEQDESC;

      if (ovp->idx.parenttype == OBJ_BIOSEQ)
      {
        vsp->bsp = ovp->idx.parentptr;        
      }
      else if (ovp->idx.parenttype == OBJ_BIOSEQSET)
      {
        vsp->bssp = ovp->idx.parentptr;        
      }
    }
    if (vsp->descr != NULL) {
      biop = vsp->descr->data.ptrvalue;
    }
  }
  
  if (biop != NULL && biop->org != NULL && biop->org->orgname != NULL)
  {
    mod = biop->org->orgname->mod;
    while (mod != NULL && mod->subtype != ORGMOD_nat_host)
    {
      mod = mod->next;
    }
    if (mod != NULL)
    {      
      ValidErr (vsp, sev, code1, code2, msg_fmt, mod->subname);
    }
  }
}


static void ReportBadSpecificHostValues (SeqEntryPtr sep, ValidStructPtr vsp)
{
  ValNodePtr    misspelled = NULL, bad_caps = NULL, ambiguous = NULL, unrecognized = NULL, vnp;

  Taxon3ValidateSpecificHostsInSeqEntry (sep, &misspelled, &bad_caps, &ambiguous, &unrecognized);

  for (vnp = misspelled; vnp != NULL; vnp = vnp->next) {
    ReportOneBadSpecificHost (vnp, vsp, SEV_WARNING, ERR_SEQ_DESCR_BadSpecificHost, "Specific host value is misspelled: %s");
  }
  for (vnp = bad_caps; vnp != NULL; vnp = vnp->next) {
    ReportOneBadSpecificHost (vnp, vsp, SEV_WARNING, ERR_SEQ_DESCR_BadSpecificHost, "Specific host value is incorrectly capitalized: %s");
  } 
  for (vnp = ambiguous; vnp != NULL; vnp = vnp->next) {
    ReportOneBadSpecificHost (vnp, vsp, SEV_INFO, ERR_SEQ_DESCR_AmbiguousSpecificHost, "Specific host value is ambiguous: %s");
  } 
  for (vnp = unrecognized; vnp != NULL; vnp = vnp->next) {
    ReportOneBadSpecificHost (vnp, vsp, SEV_WARNING, ERR_SEQ_DESCR_BadSpecificHost, "Invalid value for specific host: %s");
  } 

  misspelled = ValNodeFree (misspelled);
  bad_caps = ValNodeFree (bad_caps);
  unrecognized = ValNodeFree (unrecognized);
}

static void ReportBadTaxID (ValidStructPtr vsp, OrgRefPtr orig, OrgRefPtr reply)
{
  ValNodePtr vnp_o, vnp_r;
  DbtagPtr db_o = NULL, db_r = NULL;
  CharPtr tag1, tag2;
  Char    buf1[15];
  Char    buf2[15];

  if (vsp == NULL || orig == NULL || reply == NULL
      || orig->db == NULL || reply->db == NULL) {
    return;
  }

  for (vnp_o = orig->db; vnp_o != NULL && db_o == NULL; vnp_o = vnp_o->next) {
    if ((db_o = (DbtagPtr) vnp_o->data.ptrvalue) != NULL) {
      if (StringCmp (db_o->db, "taxon") != 0) {
        db_o = NULL;
      }
    }
  }

  if (db_o == NULL) {
    return;
  }

  for (vnp_r = reply->db; vnp_r != NULL && db_r == NULL; vnp_r = vnp_r->next) {
    if ((db_r = (DbtagPtr) vnp_r->data.ptrvalue) != NULL) {
      if (StringCmp (db_r->db, "taxon") != 0) {
        db_r = NULL;
      }
    }
  }

  if (db_r == NULL) {
    return;
  }
  if (db_o->tag->id > 0) {
    sprintf (buf1, "%d", db_o->tag->id);
    tag1 = buf1;
  } else if (db_o->tag->str == NULL) {
    sprintf (buf1, "%s", "");
    tag1 = buf1;
  } else {
    tag1 = db_o->tag->str;
  }
  if (db_r->tag->id > 0) {
    sprintf (buf2, "%d", db_r->tag->id);
    tag2 = buf2;
  } else if (db_r->tag->str == NULL) {
    sprintf (buf2, "%s", "");
    tag2 = buf2;
  } else {
    tag2 = db_r->tag->str;
  }
  if (!ObjectIdMatch (db_o->tag, db_r->tag)) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_TaxonomyLookupProblem, 
      "Organism name is '%s', taxonomy ID should be '%s' but is '%s'", orig->taxname == NULL ? "" : orig->taxname,
                                                                        tag2, tag1);
  }
}

static void RecordTentativeName (Uint2 entityID, Uint4 itemID, Uint2 itemtype, UserObjectPtr uop,
                                 TaxLstPtr tlp, SeqDescrPtr sdp, SeqFeatPtr sfp)

{
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  UserFieldPtr   curr;
  CharPtr        field;
  ObjectIdPtr    oip;
  OrgRefPtr      orp;
  ObjValNodePtr  ovp;
  SeqEntryPtr    sep;
  CharPtr        str;
  CharPtr        taxname = NULL;
  TaxValPtr      tvp;
  ValNodePtr     vnp;

  if (uop == NULL || tlp == NULL) return;

  oip = uop->type;
  if (oip == NULL || StringICmp (oip->str, "StructuredComment") != 0) return;
  for (curr = uop->data; curr != NULL; curr = curr->next) {
    if (curr->choice != 1) continue;
    oip = curr->label;
    if (oip == NULL) continue;
    field = oip->str;
    if (StringHasNoText (field)) continue;
    if (StringCmp (field, "Tentative Name") != 0) continue;
    str = (CharPtr) curr->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    if (StringCmp (str, "not provided") == 0) continue;
    taxname = str;
  }
  if (StringHasNoText (taxname)) return;

  tvp = (TaxValPtr) MemNew (sizeof (TaxVal));
  if (tvp == NULL) return;

  vnp = ValNodeNew (tlp->tail);
  if (vnp == NULL) return;

  if (tlp->head == NULL) {
    tlp->head = vnp;
  }
  tlp->tail = vnp;

  tvp->entityID = entityID;
  tvp->itemID = itemID;
  tvp->itemtype = itemtype;
  tvp->organelle = 0;

  orp = OrgRefNew ();
  if (orp == NULL) return;
  orp->taxname = StringSave (taxname);

  tvp->orp = orp;
  if (sdp != NULL && sdp->extended != 0) {
    ovp = (ObjValNodePtr) sdp;
    if (ovp->idx.parenttype == OBJ_BIOSEQ) {
      bsp = (BioseqPtr) ovp->idx.parentptr;
      if (bsp != NULL) {
        tvp->bsp = bsp;
      }
    } else if (ovp->idx.parenttype == OBJ_BIOSEQSET) {
      bssp = (BioseqSetPtr) ovp->idx.parentptr;
      if (bssp != NULL) {
        sep = bssp->seqentry;
        if (sep != NULL) {
          sep = FindNthBioseq (sep, 1);
          if (sep != NULL) {
            bsp = (BioseqPtr) sep->data.ptrvalue;
            if (bsp != NULL) {
              tvp->bsp = bsp;
            }
          }
        }
      }
    }
  } else if (sfp != NULL) {
    tvp->sfp = sfp;
  }

  vnp->data.ptrvalue = tvp;
}

static void GetTentativeNameDesc (SeqDescrPtr sdp, Pointer userdata)

{
  ObjValNodePtr  ovp;
  TaxLstPtr      tlp;
  UserObjectPtr  uop;

  if (sdp == NULL || sdp->choice != Seq_descr_user) return;
  tlp = (TaxLstPtr) userdata;

  uop = (UserObjectPtr) sdp->data.ptrvalue;
  if (uop == NULL) return;

  if (sdp->extended != 0) {
    ovp = (ObjValNodePtr) sdp;
    RecordTentativeName (ovp->idx.entityID, ovp->idx.itemID, OBJ_SEQDESC, uop, tlp, sdp, NULL);
  }
}

static void GetTentativeNameFeat (SeqFeatPtr sfp, Pointer userdata)

{
  TaxLstPtr      tlp;
  UserObjectPtr  uop;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_USER) return;
  tlp = (TaxLstPtr) userdata;

  uop = (UserObjectPtr) sfp->data.value.ptrvalue;
  if (uop == NULL) return;

  RecordTentativeName (sfp->idx.entityID, sfp->idx.itemID, OBJ_SEQFEAT, uop, tlp, NULL, sfp);
}

static void StructCommentTentativeNameValidate (SeqEntryPtr sep, ValidStructPtr vsp)

{
  GatherContext     gc;
  ValNodePtr        last = NULL;
  OrgNamePtr        onp;
  OrgRefPtr         orp;
  ErrSev            sev;
  TaxLst            srclist;
  CharPtr           str;
  T3ErrorPtr        t3ep;
  Taxon3RequestPtr  t3rq;
  Taxon3ReplyPtr    t3ry;
  T3ReplyPtr        trp;
  TaxValPtr         tvp;
  ValNodePtr        vnp;
  ValNodePtr        vnp2;

  if (sep == NULL || vsp == NULL) return;
  MemSet ((Pointer) &gc, 0, sizeof (GatherContext));
  vsp->gcp = &gc;

  srclist.head = NULL;
  srclist.tail = NULL;
  VisitDescriptorsInSep (sep, (Pointer) &srclist, GetTentativeNameDesc);
  VisitFeaturesInSep (sep, (Pointer) &srclist, GetTentativeNameFeat);
  if (srclist.head == NULL) return;

  t3rq = Taxon3RequestNew ();
  if (t3rq == NULL) return;

  for (vnp = srclist.head; vnp != NULL; vnp = vnp->next) {
    tvp = (TaxValPtr) vnp->data.ptrvalue;
    if (tvp == NULL) continue;
    orp = AsnIoMemCopy (tvp->orp,
                        (AsnReadFunc) OrgRefAsnRead,
                        (AsnWriteFunc) OrgRefAsnWrite);
    vnp2 = ValNodeAddPointer (&last, 3, (Pointer) orp);
    if (orp != NULL) {
      onp = orp->orgname;
      if (onp != NULL) {
        onp->pgcode = 0;
      }
    }
    if (t3rq->request == NULL) {
      t3rq->request = vnp2;
    }
    last = vnp2;
  }

  sev = ErrSetMessageLevel (SEV_WARNING);
  t3ry = Tax3SynchronousQuery (t3rq);
  ErrSetMessageLevel (sev);
  Taxon3RequestFree (t3rq);
  if (t3ry == NULL) return;

  for (trp = t3ry->reply, vnp = srclist.head;
       trp != NULL && vnp != NULL;
       trp = trp->next, vnp = vnp->next) {
    tvp = (TaxValPtr) vnp->data.ptrvalue;
    if (tvp == NULL) continue;
    if (trp->choice == T3Reply_error) {
      t3ep = (T3ErrorPtr) trp->data.ptrvalue;
      if (t3ep != NULL) {
        str = NULL;
        orp = (OrgRefPtr) t3ep->org;
        if (orp != NULL) {
          str = orp->taxname;
        }
        if (str == NULL) {
          str = "?";
        }

        vsp->bssp = NULL;
        vsp->bsp = tvp->bsp;
        vsp->sfp = tvp->sfp;
        vsp->descr = NULL;

        gc.entityID = tvp->entityID;
        gc.itemID = tvp->itemID;
        gc.thistype = tvp->itemtype;

        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadTentativeName, "Taxonomy lookup failed for Tentative Name '%s'", str);
      }
    }
  }

  Taxon3ReplyFree (t3ry);
  ValNodeFreeData (srclist.head);
}

static void IsINSDpatent (BioseqPtr bsp, Pointer userdata)

{
  BoolPtr   bp;
  Boolean   is_insd = FALSE;
  Boolean   is_patent = FALSE;
  SeqIdPtr  sip;

  if (bsp == NULL || userdata == NULL) return;

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    switch (sip->choice) {
      case SEQID_GENBANK :
      case SEQID_EMBL :
      case SEQID_DDBJ :
        is_insd = TRUE;
        break;
      case SEQID_PATENT :
        is_patent = TRUE;
        break;
      default :
        break;
    }
  }

  if (is_insd && is_patent) {
    bp = (BoolPtr) userdata;
    *bp = TRUE;
  }
}

static void TaxonValidate (SeqEntryPtr sep, ValidStructPtr vsp)

{
  GatherContext     gc;
  Boolean           has_nucleomorphs;
  Boolean           has_plastids;
  Boolean           is_insd_patent = FALSE;
  Boolean           is_nucleomorph;
  Boolean           is_species_level;
  Boolean           force_tax_consult;
  ValNodePtr        last = NULL;
  OrgNamePtr        onp;
  OrgRefPtr         orp;
  ErrSev            sev;
  TaxLst            srclist;
  CharPtr           str;
  T3ErrorPtr        t3ep;
  Taxon3RequestPtr  t3rq;
  Taxon3ReplyPtr    t3ry;
  T3DataPtr         tdp;
  T3StatusFlagsPtr  tfp;
  T3ReplyPtr        trp;
  TaxValPtr         tvp;
  ValNodePtr        val;
  ValNodePtr        vnp;
  ValNodePtr        vnp2;

  if (sep == NULL || vsp == NULL) return;
  MemSet ((Pointer) &gc, 0, sizeof (GatherContext));
  vsp->gcp = &gc;

  VisitBioseqsInSep (sep, (Pointer) &is_insd_patent, IsINSDpatent);

  srclist.head = NULL;
  srclist.tail = NULL;
  VisitDescriptorsInSep (sep, (Pointer) &srclist, GetSrcDesc);
  VisitFeaturesInSep (sep, (Pointer) &srclist, GetSrcFeat);
  if (srclist.head == NULL) return;

  t3rq = Taxon3RequestNew ();
  if (t3rq == NULL) return;

  for (vnp = srclist.head; vnp != NULL; vnp = vnp->next) {
    tvp = (TaxValPtr) vnp->data.ptrvalue;
    if (tvp == NULL) continue;
    orp = AsnIoMemCopy (tvp->orp,
                        (AsnReadFunc) OrgRefAsnRead,
                        (AsnWriteFunc) OrgRefAsnWrite);
    vnp2 = ValNodeAddPointer (&last, 3, (Pointer) orp);
    if (orp != NULL) {
      onp = orp->orgname;
      if (onp != NULL) {
        onp->pgcode = 0;
      }
    }
    if (t3rq->request == NULL) {
      t3rq->request = vnp2;
    }
    last = vnp2;
  }

  sev = ErrSetMessageLevel (SEV_WARNING);
  t3ry = Tax3SynchronousQuery (t3rq);
  ErrSetMessageLevel (sev);
  Taxon3RequestFree (t3rq);
  if (t3ry == NULL) return;

  for (trp = t3ry->reply, vnp = srclist.head;
       trp != NULL && vnp != NULL;
       trp = trp->next, vnp = vnp->next) {
    tvp = (TaxValPtr) vnp->data.ptrvalue;
    if (tvp == NULL) continue;
    if (trp->choice == T3Reply_error) {
      t3ep = (T3ErrorPtr) trp->data.ptrvalue;
      if (t3ep != NULL) {
        str = t3ep->message;
        if (str == NULL) {
          str = "?";
        }

        vsp->bssp = NULL;
        vsp->bsp = tvp->bsp;
        vsp->sfp = tvp->sfp;
        vsp->descr = NULL;

        gc.entityID = tvp->entityID;
        gc.itemID = tvp->itemID;
        gc.thistype = tvp->itemtype;

        if (StringCmp (str, "Organism not found") == 0) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_OrganismNotFound, "Organism not found in taxonomy database");
        } else {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_TaxonomyLookupProblem, "Taxonomy lookup failed with message '%s'", str);
        }
      }
    }
    if (trp->choice != T3Reply_data) continue;
    tdp = (T3DataPtr) trp->data.ptrvalue;
    if (tdp == NULL) continue;

    vsp->bssp = NULL;
    vsp->bsp = tvp->bsp;
    vsp->sfp = tvp->sfp;
    vsp->descr = NULL;

    ReportBadTaxID (vsp, tvp->orp, (OrgRefPtr) tdp->org);

    is_species_level = FALSE;
    has_nucleomorphs = FALSE;
    is_nucleomorph = FALSE;
    has_plastids = FALSE;

    for (tfp = tdp->status; tfp != NULL; tfp = tfp->next) {

      /*
      val = tfp->Value_value;
      if (val != NULL && val->choice == Value_value_bool) {
        str = tfp->property;
        if (str == NULL) {
          str = "?";
        }
        if (val->data.intvalue != 0) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_TaxonomyLookupProblem, "'%s' TRUE", str);
        } else {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_TaxonomyLookupProblem, "'%s' FALSE", str);
        }
      }
      */

      if (StringICmp (tfp->property, "is_species_level") == 0) {
        val = tfp->Value_value;
        if (val != NULL && val->choice == Value_value_bool) {
          is_species_level = (Boolean) (val->data.intvalue != 0);
          if (! is_species_level) {
            gc.entityID = tvp->entityID;
            gc.itemID = tvp->itemID;
            gc.thistype = tvp->itemtype;

            if (! vsp->is_wp_in_sep) {
              ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_TaxonomyIsSpeciesProblem, "Taxonomy lookup reports is_species_level FALSE");
            }
          }
        }
      } else if (StringICmp (tfp->property, "force_consult") == 0) {
        val = tfp->Value_value;
        if (val != NULL && val->choice == Value_value_bool) {
          force_tax_consult = (Boolean) (val->data.intvalue != 0);
          if (force_tax_consult && is_insd_patent) {
            orp = (OrgRefPtr) tdp->org;
            if (orp != NULL) {
              if (StringICmp (orp->taxname, "unidentified") == 0) {
                force_tax_consult = FALSE;
              }
            }
          }
          if (force_tax_consult) {
            gc.entityID = tvp->entityID;
            gc.itemID = tvp->itemID;
            gc.thistype = tvp->itemtype;

            ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_TaxonomyConsultRequired, "Taxonomy lookup reports taxonomy consultation needed");
          }
        }
      } else if (StringICmp (tfp->property, "has_nucleomorphs") == 0) {
        val = tfp->Value_value;
        if (val != NULL && val->choice == Value_value_bool) {
          has_nucleomorphs = (Boolean) (val->data.intvalue != 0);
          if (has_nucleomorphs) {
            is_nucleomorph = TRUE;
          }
        }
      } else if (StringICmp (tfp->property, "has_plastids") == 0) {
        val = tfp->Value_value;
        if (val != NULL && val->choice == Value_value_bool) {
          has_plastids = (Boolean) (val->data.intvalue != 0);
        }
      }
    }
    if (tvp->organelle == GENOME_nucleomorph && (! is_nucleomorph)) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_TaxonomyNucleomorphProblem, "Taxonomy lookup does not have expected nucleomorph flag");
    } else if (tvp->organelle == GENOME_plastid && (! has_plastids)) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_TaxonomyPlastidsProblem, "Taxonomy lookup does not have expected plastid flag");
    }

  }

  Taxon3ReplyFree (t3ry);
  ValNodeFreeData (srclist.head);

  /* also validate specific-host values */

  ReportBadSpecificHostValues (sep, vsp);

  StructCommentTentativeNameValidate (sep, vsp);
}

static void DoValidation (
  SeqEntryPtr sep,
  ValFlagPtr vfp,
  FILE *ofp
)

{
  Int2            i;
  VCData          vcd;
  ValidStructPtr  vsp;
  ErrSev          cutoff;
  CharPtr         xml_header = NULL;

  if (vfp == NULL) return;

  vsp = ValidStructNew ();
  if (vsp == NULL) return;

  MemSet ((Pointer) &vcd, 0, sizeof (VCData));

  vsp->useSeqMgrIndexes = TRUE;

  vsp->cutoff = vfp->lowCutoff;
  vsp->validateAlignments = vfp->validateAlignments;
  vsp->alignFindRemoteBsp = vfp->alignFindRemoteBsp;
  vsp->doSeqHistAssembly = vfp->doSeqHistAssembly;
  vsp->farIDsInAlignments = vfp->farIDsInAlignments;
  vsp->alwaysRequireIsoJTA = vfp->alwaysRequireIsoJTA;
  vsp->farFetchCDSproducts = vfp->farFetchCDSproducts;
  vsp->farFetchMRNAproducts = vfp->farFetchMRNAproducts;
  vsp->locusTagGeneralMatch = vfp->locusTagGeneralMatch;
  vsp->validateIDSet = vfp->validateIDSet;
  vsp->seqSubmitParent = vfp->seqSubmitParent;
  vsp->ignoreExceptions = vfp->ignoreExceptions;
  vsp->validateExons = vfp->validateExons;
  vsp->inferenceAccnCheck = vfp->inferenceAccnCheck;
  vsp->testLatLonSubregion = vfp->testLatLonSubregion;
  vsp->strictLatLonCountry = vfp->strictLatLonCountry;
  vsp->rubiscoTest = vfp->rubiscoTest;
  vsp->disableSuppression = vfp->disableSuppression;
  vsp->genomeSubmission = vfp->genomeSubmission;
  vsp->debugTestDuJour = vfp->debugTestDuJour;
  vsp->indexerVersion = vfp->indexerVersion;

  if (vfp->verbosity == 4) {
    vsp->use_heartbeat = TRUE;
  }

  if (ofp == NULL && vfp->outfp != NULL) {
    ofp = vfp->outfp;
  }
  if (ofp != NULL) {
    vcd.ofp = ofp;
    vcd.verbosity = vfp->verbosity;
    vcd.lowCutoff = vfp->lowCutoff;
    vcd.highCutoff = vfp->highCutoff;
    vcd.errcode = vfp->errcode;
    vcd.vfp = vfp;
    vsp->errfunc = ValidCallback;
    vsp->userdata = (Pointer) &vcd;
    vsp->convertGiToAccn = FALSE;
  }

  Heartbeat(vsp, "Beginning validation");

  ValidateSeqEntry (sep, vsp);

  if (vfp->taxfetch) {
    TaxonValidate (sep, vsp);
  }

  for (i = 0; i <= 4; i++) {
    vfp->num_errors += vsp->errors [i];
    if (i >= vfp->severity) {
      vfp->fatal_errors += vsp->errors [i];
    }
  }

  Heartbeat(vsp, "Finished Validation");

  ValidStructFree (vsp);
  if (vfp->validateBarcode) {
    if (vfp->verbosity == 4 && !vfp->has_errors) {
      cutoff = (ErrSev) vfp->lowCutoff;
      if (cutoff < SEV_NONE || cutoff > SEV_MAX) {
        cutoff = SEV_MAX;
      }
      xml_header = GetXmlHeaderText(cutoff);
    }
    if (!BarcodeValidateOneSeqEntry (ofp, sep, TRUE,
                                     vfp->verbosity == 4,
                                     !vfp->has_errors,
                                     xml_header)) {
      vfp->has_errors = TRUE;
    }
    xml_header = MemFree (xml_header);
  }
}

static void ValidateOneSep (
  SeqEntryPtr sep,
  FILE *ofp,
  ValFlagPtr vfp
)
{
  ValNodePtr bsplist = NULL;

  if (! TooManyFarComponents (sep)) {
    if (vfp->inferenceAccnCheck) {
      if (! TooManyInferenceAccessions (sep, NULL, NULL)) {
        LookupFarSeqIDs (sep, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE);
      }
    }
    if (vfp->lock) {
      bsplist = DoLockFarComponents (sep, vfp);
    }
  }

  DoValidation (sep, vfp, ofp);

  bsplist = UnlockFarComponents (bsplist);
}

static void ReportReadFailure (
  ValFlagPtr vfp
)

{
  FILE     *ofp;
  CharPtr  xml_header;

  if (vfp == NULL) return;
  ofp = vfp->outfp;
  if (ofp == NULL) return;

  if (vfp->verbosity == 4) {

    xml_header = GetXmlHeaderText (0);
    fprintf (ofp, "<%s>\n", xml_header);
    xml_header = MemFree (xml_header);

    fprintf (ofp, "  <message severity=\"REJECT\" code=\"GENERIC_InvalidAsn\">Unable to read invalid ASN.1</message>\n");

    fprintf (ofp, "</asnval>\n");

    return;
  }

  fprintf (ofp, "REJECT: valid [GENERIC.InvalidAsn] Unable to read invalid ASN.1\n");
}

static void ProcessSingleRecord (
  CharPtr filename,
  ValFlagPtr vfp
)

{
  AsnIoPtr       aip;
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  Char           buf [64], path [PATH_MAX];
  Pointer        dataptr = NULL;
  Uint2          datatype = 0, entityID = 0;
  FILE           *fp, *ofp = NULL;
  SeqEntryPtr    fsep, sep, tmp_sep;
  ObjMgrPtr      omp;
  CharPtr        ptr;
  time_t         starttime, stoptime;
  SeqSubmitPtr   ssp;

  if (StringHasNoText (filename)) return;
  if (vfp == NULL) return;

  if (vfp->type == 1) {
    fp = FileOpen (filename, "r");
    if (fp == NULL) {
      Message (MSG_POSTERR, "Failed to open '%s'", filename);
      return;
    }

    dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, FALSE, FALSE, TRUE, FALSE);

    FileClose (fp);

    entityID = ObjMgrRegister (datatype, dataptr);

  } else if (vfp->type >= 2 && vfp->type <= 5) {
    aip = AsnIoOpen (filename, vfp->binary? "rb" : "r");
    if (aip == NULL) {
      Message (MSG_POSTERR, "AsnIoOpen failed for input file '%s'", filename);
      return;
    }

    SeqMgrHoldIndexing (TRUE);
    switch (vfp->type) {
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
    Message (MSG_POSTERR, "Input format type '%d' unrecognized", (int) vfp->type);
    return;
  }

  if (entityID < 1 || dataptr == NULL) {
    Message (MSG_POSTERR, "Data read failed for input file '%s'", filename);
    vfp->fatal_errors++;
    ReportReadFailure (vfp);
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

      starttime = GetSecs ();
      buf [0] = '\0';

      if (vfp->logfp != NULL) {
        fsep = FindNthBioseq (sep, 1);
        if (fsep != NULL && fsep->choice == 1) {
          bsp = (BioseqPtr) fsep->data.ptrvalue;
          if (bsp != NULL) {
            SeqIdWrite (bsp->id, buf, PRINTID_FASTA_LONG, sizeof (buf));
            fprintf (vfp->logfp, "%s\n", buf);
            fflush (vfp->logfp);
          }
        }
      }

      StringNCpy_0 (path, filename, sizeof (path));
      ptr = StringRChr (path, '.');
      if (ptr != NULL) {
        *ptr = '\0';
      }
      StringCat (path, ".val");

      if (vfp->outpath != NULL) {
        ErrSetLogfile (vfp->outpath, ELOG_APPEND);
        ErrSetLogLevel (SEV_INFO);
      } else if (vfp->verbosity == 0 || vfp->verbosity == 6) {
        ErrSetLogfile (path, ELOG_APPEND);
        ErrSetLogLevel (SEV_INFO);
      } else if (vfp->outfp == NULL) {
        ofp = FileOpen (path, "w");
      }

      if (datatype == OBJ_SEQSUB && (ssp = (SeqSubmitPtr)dataptr) != NULL
          && ssp->datatype == 1) {
        for (sep = ssp->data; sep != NULL; sep = sep->next) {
          tmp_sep = (SeqEntryPtr) AsnIoMemCopy (sep, (AsnReadFunc) SeqEntryAsnRead, (AsnWriteFunc) SeqEntryAsnWrite);
          entityID = ObjMgrGetEntityIDForChoice (tmp_sep);
          SeqMgrIndexFeatures (entityID, NULL);
          ValidateOneSep (tmp_sep, ofp, vfp);
          tmp_sep = SeqEntryFree (tmp_sep);
        }                   
      } else {
        ValidateOneSep (sep, ofp, vfp);
      } 

      if (ofp != NULL) {
        if (vfp->has_errors) {
          if (vfp->verbosity == 4) {
            fprintf (ofp, "</asnval>\n");
          }
          vfp->has_errors = FALSE;
        }
        FileClose (ofp);
      }

      stoptime = GetSecs ();
      if (stoptime - starttime > vfp->worsttime && StringDoesHaveText (buf)) {
        vfp->worsttime = stoptime - starttime;
        StringCpy (vfp->longest, buf);
      }
      (vfp->numrecords)++;
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
  CharPtr filename,
  ValFlagPtr vfp
)

{
  AsnIoPtr        aip;
  AsnModulePtr    amp;
  AsnTypePtr      atp, atp_bss, atp_desc, atp_sbp, atp_se = NULL, atp_ssp;
  BioseqPtr       bsp;
  ValNodePtr      bsplist;
  BioseqSetPtr    bssp;
  Char            buf [64], path [PATH_MAX], longest [64];
  Int2            skipcount = 0, maxcount = 0;
  CitSubPtr       csp = NULL;
  FILE            *fp, *ofp = NULL;
  Int4            numrecords = 0;
  SeqEntryPtr     fsep, sep;
  ObjMgrPtr       omp;
  ObjValNode      ovn;
  Pubdesc         pd;
  CharPtr         ptr;
  SubmitBlockPtr  sbp = NULL;
  time_t          starttime, stoptime, worsttime;
  SeqDescrPtr     subcit = NULL;
  ValNode         vn;
#ifdef OS_UNIX
  Char            cmmd [256];
  Boolean         detailed_report = FALSE;
  CharPtr         gzcatprog;
  Boolean         memory_usage = FALSE;
  int             ret;
  Boolean         usedPopen = FALSE;
#endif

  if (StringHasNoText (filename)) return;
  if (vfp == NULL) return;

#ifndef OS_UNIX
  if (vfp->compressed) {
    Message (MSG_POSTERR, "Can only decompress on-the-fly on UNIX machines");
    return;
  }
#endif

  amp = AsnAllModPtr ();
  if (amp == NULL) {
    Message (MSG_POSTERR, "Unable to load AsnAllModPtr");
    return;
  }

  atp_ssp = AsnFind ("Seq-submit");
  if (atp_ssp == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Seq-submit");
    return;
  }

  atp_sbp = AsnFind ("Seq-submit.sub");
  if (atp_sbp == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type Seq-submit.sub");
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

  if (vfp->type == 4) {
    atp_se = AsnFind ("Bioseq-set.seq-set.E");
    if (atp_se == NULL) {
      Message (MSG_POSTERR, "Unable to find ASN.1 type Bioseq-set.seq-set.E");
      return;
    }
  } else if (vfp->type == 5) {
    /*
    atp_se = AsnFind ("Seq-submit.data.entrys.E");
    if (atp_se == NULL) {
      Message (MSG_POSTERR, "Unable to find ASN.1 type Seq-submit.data.entrys.E");
      return;
    }
    */
    /* current use has genbank set containing batch, so iterate genbank set within Seq-submit */
    atp_se = AsnFind ("Bioseq-set.seq-set.E");
    if (atp_se == NULL) {
      Message (MSG_POSTERR, "Unable to find ASN.1 type Bioseq-set.seq-set.E");
      return;
    }
  } else {
    Message (MSG_POSTERR, "Batch processing type not set properly");
    return;
  }

  if (atp_se == NULL) {
    Message (MSG_POSTERR, "Unable to find ASN.1 type for atp_se");
    return;
  }

#ifdef OS_UNIX
  if (getenv ("ASNVAL_LOG_OBJMGR_REPORT") != NULL) {
    detailed_report = TRUE;
  }
  if (getenv ("ASNVAL_LOG_MEMORY_REPORT") != NULL) {
    memory_usage = TRUE;
  }

  if (vfp->compressed) {
    gzcatprog = getenv ("NCBI_UNCOMPRESS_BINARY");
    if (gzcatprog != NULL) {
      sprintf (cmmd, "%s %s", gzcatprog, filename);
    } else {
      ret = system ("gzcat -h >/dev/null 2>&1");
      if (ret == 0) {
        sprintf (cmmd, "gzcat %s", filename);
      } else if (ret == -1) {
        Message (MSG_POSTERR, "Unable to fork or exec gzcat in ScanBioseqSetRelease");
        return;
      } else {
        ret = system ("zcat -h >/dev/null 2>&1");
        if (ret == 0) {
          sprintf (cmmd, "zcat %s", filename);
        } else if (ret == -1) {
          Message (MSG_POSTERR, "Unable to fork or exec zcat in ScanBioseqSetRelease");
          return;
        } else {
          Message (MSG_POSTERR, "Unable to find zcat or gzcat in ScanBioseqSetRelease - please edit your PATH environment variable");
          return;
        }
      }
    }
    fp = popen (cmmd, /* vfp->binary? "rb" : */ "r");
    usedPopen = TRUE;
  } else {
    fp = FileOpen (filename, vfp->binary? "rb" : "r");
  }
#else
  fp = FileOpen (filename, vfp->binary? "rb" : "r");
#endif
  if (fp == NULL) {
    Message (MSG_POSTERR, "FileOpen failed for input file '%s'", filename);
    return;
  }

  aip = AsnIoNew (vfp->binary? ASNIO_BIN_IN : ASNIO_TEXT_IN, fp, NULL, NULL, NULL);
  if (aip == NULL) {
    Message (MSG_POSTERR, "AsnIoNew failed for input file '%s'", filename);
    return;
  }

  if (vfp->type == 4) {
    atp = atp_bss;
  } else if (vfp->type == 5) {
    atp = atp_ssp;
  } else {
    Message (MSG_POSTERR, "Batch processing type not set properly");
    return;
  }

  longest [0] = '\0';
  worsttime = 0;

  StringNCpy_0 (path, filename, sizeof (path));
  ptr = StringRChr (path, '.');
  if (ptr != NULL) {
    *ptr = '\0';
  }
  StringCat (path, ".val");

  if (vfp->outpath != NULL) {
    ErrSetLogfile (vfp->outpath, ELOG_APPEND);
    ErrSetLogLevel (SEV_INFO);
  } else if (vfp->verbosity == 0 || vfp->verbosity == 6) {
    ErrSetLogfile (path, ELOG_APPEND);
    ErrSetLogLevel (SEV_INFO);
  } else if (vfp->outfp == NULL) {
    ofp = FileOpen (path, "w");
  }

  while ((! vfp->io_failure) && maxcount < vfp->maxcount &&
         (atp = AsnReadId (aip, amp, atp)) != NULL) {
    if (aip->io_failure) {
      vfp->io_failure = TRUE;
      aip->io_failure = FALSE;
    }
    if (atp == atp_se) {

      SeqMgrHoldIndexing (TRUE);
      sep = SeqEntryAsnRead (aip, atp);
      SeqMgrHoldIndexing (FALSE);

      if (sep == NULL) {
        vfp->fatal_errors++;
        ReportReadFailure (vfp);
      }

      /* propagate submission citation as descriptor onto each Seq-entry */

      if (subcit != NULL && sep != NULL && sep->data.ptrvalue != NULL) {
        if (sep->choice == 1) {
          bsp = (BioseqPtr) sep->data.ptrvalue;
          ValNodeLink (&(bsp->descr),
                       AsnIoMemCopy ((Pointer) subcit,
                                     (AsnReadFunc) SeqDescrAsnRead,
                                     (AsnWriteFunc) SeqDescrAsnWrite));
        } else if (sep->choice == 2) {
          bssp = (BioseqSetPtr) sep->data.ptrvalue;
          ValNodeLink (&(bssp->descr),
                       AsnIoMemCopy ((Pointer) subcit,
                                     (AsnReadFunc) SeqDescrAsnRead,
                                     (AsnWriteFunc) SeqDescrAsnWrite));
        }
      }

      if (sep != NULL) {
        if (skipcount < vfp->skipcount) {
          skipcount++;
        } else {

          starttime = GetSecs ();
          buf [0] = '\0';

          if (vfp->logfp != NULL) {
            fsep = FindNthBioseq (sep, 1);
            if (fsep != NULL && fsep->choice == 1) {
              bsp = (BioseqPtr) fsep->data.ptrvalue;
              if (bsp != NULL) {
                SeqIdWrite (bsp->id, buf, PRINTID_FASTA_LONG, sizeof (buf));
                fprintf (vfp->logfp, "%s\n", buf);
                fflush (vfp->logfp);
              }
            }
          }

          bsplist = NULL;

          if (! TooManyFarComponents (sep)) {
            if (vfp->inferenceAccnCheck) {
              if (! TooManyInferenceAccessions (sep, NULL, NULL)) {
                LookupFarSeqIDs (sep, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE);
              }
            }
            if (vfp->lock) {
              bsplist = DoLockFarComponents (sep, vfp);
            }
          }

          DoValidation (sep, vfp, ofp);

          bsplist = UnlockFarComponents (bsplist);

          stoptime = GetSecs ();
          if (stoptime - starttime > worsttime && StringDoesHaveText (buf)) {
            worsttime = stoptime - starttime;
            StringCpy (longest, buf);
          }
          numrecords++;
          maxcount++;
        }
      }

      SeqEntryFree (sep);
      omp = ObjMgrGet ();
      ObjMgrReapOne (omp);
      SeqMgrClearBioseqIndex ();
      ObjMgrFreeCache (0);
      FreeSeqIdGiCache ();

      SeqEntrySetScope (NULL);

#ifdef OS_UNIX
      if (detailed_report && vfp->logfp != NULL) {
        ObjMgrReportProc (vfp->logfp);
      }

      if (memory_usage && vfp->logfp != NULL) {
        Char     mbuf [512];
        FILE     *mufp;
        Char     ch;
        Int4     len;
        CharPtr  ptr1, ptr2;
        Int2     spaces;
        uid_t    uid;
        uid = getpid ();
        sprintf (cmmd, "cat /proc/%d/stat", (int) uid);
        mufp = popen (cmmd, "r");
        if (mufp != NULL) {
          len = FileRead ((Pointer) mbuf, sizeof (Char), sizeof (mbuf), mufp);
          if (len > 0) {
            mbuf [(int) len] = '\0';
            ptr1 = mbuf;
            ch = *ptr1;
            spaces = 0;
            while (ch != '\0' && spaces < 22) {
              if (ch == ' ') {
                spaces++;
              }
              ptr1++;
              ch = *ptr1;
            }
            if (ch != '\0') {
              ptr2 = StringChr (ptr1, ' ');
              if (ptr2 != NULL) {
                *ptr2 = '\0';
                fprintf (vfp->logfp, "Memory usage %s\n", ptr1);
              }
            }
          }
          pclose (mufp);
        }
      }
#endif

    } else if (atp == atp_sbp) {
      sbp = SubmitBlockAsnRead (aip, atp);
      if (sbp != NULL) {
        csp = sbp->cit;
        if (csp != NULL) {
          MemSet ((Pointer) &ovn, 0, sizeof (ObjValNode));
          MemSet ((Pointer) &pd, 0, sizeof (Pubdesc));
          MemSet ((Pointer) &vn, 0, sizeof (ValNode));
          vn.choice = PUB_Sub;
          vn.data.ptrvalue = (Pointer) csp;
          vn.next = NULL;
          pd.pub = &vn;
          ovn.vn.choice = Seq_descr_pub;
          ovn.vn.data.ptrvalue = (Pointer) &pd;
          ovn.vn.next = NULL;
          ovn.vn.extended = 1;
          subcit = (SeqDescrPtr) &ovn;
        }
      }
    } else {
      AsnReadVal (aip, atp, NULL);
    }

    if (aip->io_failure) {
      vfp->io_failure = TRUE;
      aip->io_failure = FALSE;
    }
  }

  if (aip->io_failure) {
    vfp->io_failure = TRUE;
  }

  if (vfp->io_failure) {
    Message (MSG_POSTERR, "Asn io_failure for input file '%s'", filename);
  }

  if (ofp != NULL) {
    if (vfp->has_errors) {
      if (vfp->verbosity == 4) {
        fprintf (ofp, "</asnval>\n");
      }
      vfp->has_errors = FALSE;
    }
    FileClose (ofp);
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

  if (vfp->logfp != NULL && (! StringHasNoText (longest))) {
    fprintf (vfp->logfp, "Longest processing time %ld seconds on %s\n",
             (long) worsttime, longest);
    fprintf (vfp->logfp, "Total number of records %ld\n", (long) numrecords);
    fflush (vfp->logfp);
  }
}

static void ValidWrapper (
  SeqEntryPtr sep,
  Pointer userdata
)

{
  BioseqPtr      bsp;
  ValNodePtr     bsplist;
  Char           buf [64];
  SeqEntryPtr    fsep;
  FILE           *ofp = NULL;
  CharPtr        ptr;
  ErrSev         sev;
  time_t         starttime, stoptime;
  ValFlagPtr     vfp;

  if (sep == NULL) return;
  vfp = (ValFlagPtr) userdata;
  if (vfp == NULL) return;

  starttime = GetSecs ();
  buf [0] = '\0';

  if (vfp->logfp != NULL) {
    fsep = FindNthBioseq (sep, 1);
    if (fsep != NULL && fsep->choice == 1) {
      bsp = (BioseqPtr) fsep->data.ptrvalue;
      if (bsp != NULL) {
        SeqIdWrite (bsp->id, buf, PRINTID_FASTA_LONG, sizeof (buf));
        fprintf (vfp->logfp, "%s\n", buf);
        fflush (vfp->logfp);
      }
    }
  }

  ptr = StringRChr (vfp->path, '.');
  if (ptr != NULL) {
    *ptr = '\0';
  }
  StringCat (vfp->path, ".val");

  if (vfp->outpath != NULL) {
    ErrSetLogfile (vfp->outpath, ELOG_APPEND);
    ErrSetLogLevel (SEV_INFO);
  } else if (vfp->verbosity == 0 || vfp->verbosity == 6) {
    ErrSetLogfile (vfp->path, ELOG_APPEND);
    ErrSetLogLevel (SEV_INFO);
  } else if (vfp->outfp == NULL) {
    ofp = FileOpen (vfp->path, "w");
  }

  bsplist = NULL;

  if (! TooManyFarComponents (sep)) {
    sev = ErrSetMessageLevel (SEV_WARNING);
    if (vfp->inferenceAccnCheck) {
      if (! TooManyInferenceAccessions (sep, NULL, NULL)) {
        LookupFarSeqIDs (sep, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE);
      }
    }
    if (vfp->lock) {
      bsplist = DoLockFarComponents (sep, vfp);
    }
    ErrSetMessageLevel (sev);
  }

  DoValidation (sep, vfp, ofp);

  bsplist = UnlockFarComponents (bsplist);

  if (ofp != NULL) {
    if (vfp->has_errors) {
      if (vfp->verbosity == 4) {
        fprintf (ofp, "</asnval>\n");
      }
      vfp->has_errors = FALSE;
    }
    FileClose (ofp);
  }

  stoptime = GetSecs ();
  if (stoptime - starttime > vfp->worsttime && StringDoesHaveText (buf)) {
    vfp->worsttime = stoptime - starttime;
    StringCpy (vfp->longest, buf);
  }
  (vfp->numrecords)++;
}

static void ProcessOneRecord (
  CharPtr filename,
  Pointer userdata
)

{
  Pointer      dataptr;
  Uint2        datatype;
  Uint2        entityID;
  FILE         *fp;
  ObjMgrPtr    omp;
  SeqEntryPtr  sep, tmp_sep;
  ValFlagPtr   vfp;
  SeqSubmitPtr ssp;

  vfp = (ValFlagPtr) userdata;
  if (vfp == NULL) return;

  if (vfp->logfp != NULL) {
    fprintf (vfp->logfp, "%s\n", filename);
    fflush (vfp->logfp);
  }

  if (vfp->automatic) {
    StringNCpy_0 (vfp->path, filename, sizeof (vfp->path));
    if (!ReadSequenceAsnFile (filename, vfp->binary, vfp->compressed, (Pointer) vfp, ValidWrapper)) {
      vfp->fatal_errors++;
      ReportReadFailure (vfp);
    }
  } else if (vfp->catenated) {
    fp = FileOpen (filename, "r");
    if (fp != NULL) {

      SeqMgrHoldIndexing (TRUE);
      dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, &entityID, FALSE, FALSE, TRUE, FALSE);
      SeqMgrHoldIndexing (FALSE);

      if (dataptr == NULL) {
        Message (MSG_FATAL, "Unable to read from file %s", filename);
        vfp->fatal_errors++;
        ReportReadFailure (vfp);
      }
      while (dataptr != NULL) {
        if (datatype == OBJ_SEQSUB && (ssp = (SeqSubmitPtr)dataptr) != NULL
            && ssp->datatype == 1) {
          for (sep = ssp->data; sep != NULL; sep = sep->next) {
            tmp_sep = (SeqEntryPtr) AsnIoMemCopy (sep, (AsnReadFunc) SeqEntryAsnRead, (AsnWriteFunc) SeqEntryAsnWrite);
            entityID = ObjMgrGetEntityIDForChoice (tmp_sep);
            SeqMgrIndexFeatures (entityID, NULL);
            ValidWrapper (tmp_sep, vfp);
            tmp_sep = SeqEntryFree (tmp_sep);
          }          
        } else {
          sep = GetTopSeqEntryForEntityID (entityID);
          ValidWrapper (sep, vfp);
        }

        ObjMgrFree (datatype, dataptr);
  
        omp = ObjMgrGet ();
        ObjMgrReapOne (omp);
        SeqMgrClearBioseqIndex ();
        ObjMgrFreeCache (0);
        FreeSeqIdGiCache ();
  
        SeqEntrySetScope (NULL);

        SeqMgrHoldIndexing (TRUE);
        dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, &entityID, FALSE, FALSE, TRUE, FALSE);
        SeqMgrHoldIndexing (FALSE);
      }
      FileClose (fp);
    }
  } else if (vfp->batch) {
    ProcessMultipleRecord (filename, vfp);
  } else {
    ProcessSingleRecord (filename, vfp);
  }
}

static QUEUE bouncequeue = NULL;

static Boolean LIBCALLBACK BounceProc (
  CONN conn, VoidPtr userdata, EIO_Status status
)

{
  BoolPtr  bp;

  if (NetTestReadReply (conn, status)) {
    bp = (BoolPtr) userdata;
    *bp = TRUE;
  }
  return TRUE;
}

/* Args structure contains command-line arguments */

typedef enum {
  p_argInputPath = 0,
  i_argInputFile,
  o_argOutputFile,
  f_argFilter,
  x_argSuffix,
  u_argRecurse,
  R_argSeverity,
  Q_argLowCutoff,
  P_argHighCutoff,
  E_argOnlyThisErr,
  A_argAlignments,
  J_argIsoJta,
  Z_argRemoteCDS,
  X_argExonSplice,
  G_argInfAccns,
  N_argLatLonStrict,
  M_argMatchTag,
  Y_argCheckOld,
  e_argIgnoreExcept,
  v_argVerbosity,
  a_argType,
  b_argBinary,
  c_argCompressed,
  r_argRemote,
  k_argLocalFetch,
  d_argAsnIdx,
  q_argTaxLookup,
  l_argLockFar,
  T_argThreads,
  F_argTestNetwork,
  L_argLogFile,
  K_argSummary,
  D_argDisableSuppress,
  U_argGenomeSubmission,
  S_argSkipCount,
  B_argBarcodeVal,
  C_argMaxCount,
  j_argDebugDuJour,
  w_argSeqSubParent,
  y_argAIndexer,
#ifdef INTERNAL_NCBI_ASN2VAL
  H_argAccessHUP,
#endif
} Arguments;

#define LAT_LON_STATE    1
#define LAT_LON_STRICT   2

Args myargs [] = {
  {"Path to ASN.1 Files", NULL, NULL, NULL,
    TRUE, 'p', ARG_STRING, 0.0, 0, NULL},
  {"Single Input File", "stdin", NULL, NULL,
    TRUE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
  {"Single Output File", NULL, NULL, NULL,
    TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
  {"Substring Filter", NULL, NULL, NULL,
    TRUE, 'f', ARG_STRING, 0.0, 0, NULL},
  {"File Selection Substring", ".ent", NULL, NULL,
    TRUE, 'x', ARG_STRING, 0.0, 0, NULL},
  {"Recurse", "F", NULL, NULL,
    TRUE, 'u', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Severity for Error in Return Code", "4", "0", "6",
    FALSE, 'R', ARG_INT, 0.0, 0, NULL},
  {"Lowest Severity for Error to Show", "3", "0", "4",
    FALSE, 'Q', ARG_INT, 0.0, 0, NULL},
  {"Highest Severity for Error to Show", "4", "0", "4",
    FALSE, 'P', ARG_INT, 0.0, 0, NULL},
  {"Only Error Code to Show", NULL, NULL, NULL,
    TRUE, 'E', ARG_STRING, 0.0, 0, NULL},
  {"Validate Alignments", "F", NULL, NULL,
    TRUE, 'A', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Require ISO-JTA?", "F", NULL, NULL,
    TRUE, 'J', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Remote CDS Product Fetch", "F", NULL, NULL,
    TRUE, 'Z', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Exon Splice Check", "F", NULL, NULL,
    TRUE, 'X', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Verify Inference Accessions", "F", NULL, NULL,
    TRUE, 'G', ARG_BOOLEAN, 0.0, 0, NULL},
  {"LatLon/Country Flags (1 Test State/Province, 2 Ignore Water Exception)", "0", "0", "3",
    TRUE, 'N', ARG_INT, 0.0, 0, NULL},
  {"Match locus_tag against General ID", "F", NULL, NULL,
    TRUE, 'M', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Check Against Old IDs", "F", NULL, NULL,
    TRUE, 'Y', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Ignore Transcription/Translation Exceptions", "F", NULL, NULL,
    TRUE, 'e', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Verbosity\n"
   "      1 Standard Report\n"
   "      2 Accession / Severity / Code (space delimited)\n"
   "      3 Accession / Severity / Code (tab delimited\n"
   "      4 XML Report\n"
   "      5 Accession / GI / Severity / Code (tab delimited)\n"
   "      6 Logged Report\n", "1", "0", "6",
    FALSE, 'v', ARG_INT, 0.0, 0, NULL},
  {"ASN.1 Type\n"
   "      a Automatic\n"
   "      c Catenated\n"
   "      z Any\n"
   "      e Seq-entry\n"
   "      b Bioseq\n"
   "      s Bioseq-set\n"
   "      m Seq-submit\n"
   "      t Batch Bioseq-set\n"
   "      u Batch Seq-submit\n", "a", NULL, NULL,
    TRUE, 'a', ARG_STRING, 0.0, 0, NULL},
  {"Batch File is Binary", "F", NULL, NULL,
    TRUE, 'b', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Batch File is Compressed", "F", NULL, NULL,
    TRUE, 'c', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Remote Fetching from ID", "F", NULL, NULL,
    TRUE, 'r', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Local Fetching", "F", NULL, NULL,
    TRUE, 'k', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Path to Indexed Binary ASN.1 Data", NULL, NULL, NULL,
    TRUE, 'd', ARG_STRING, 0.0, 0, NULL},
  {"Taxonomy Lookup", "F", NULL, NULL,
    TRUE, 'q', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Lock Components in Advance", "F", NULL, NULL,
    TRUE, 'l', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Use Threads", "F", NULL, NULL,
    TRUE, 'T', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Test Network Access", "F", NULL, NULL,
    TRUE, 'F', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Log File", NULL, NULL, NULL,
    TRUE, 'L', ARG_FILE_OUT, 0.0, 0, NULL},
  {"Summary to Error File", "F", NULL, NULL,
    TRUE, 'K', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Disable Message Suppression", "F", NULL, NULL,
    TRUE, 'D', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Genome Center Submission", "F", NULL, NULL,
    TRUE, 'U', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Skip Count", "0", NULL, NULL,
    TRUE, 'S', ARG_INT, 0.0, 0, NULL},
  {"Barcode Validate", "F", NULL, NULL,
    TRUE, 'B', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Max Count", "0", NULL, NULL,
    TRUE, 'C', ARG_INT, 0.0, 0, NULL},
  {"Validator Code Performance Test", "F", NULL, NULL,
    TRUE, 'j', ARG_BOOLEAN, 0.0, 0, NULL},
  {"SeqSubmitParent Flag", "F", NULL, NULL,
    TRUE, 'w', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Special Indexer Tests", "F", NULL, NULL,
    TRUE, 'y', ARG_BOOLEAN, 0.0, 0, NULL},
#ifdef INTERNAL_NCBI_ASN2VAL
  {"Internal Access to HUP", "F", NULL, NULL,
    TRUE, 'H', ARG_BOOLEAN, 0.0, 0, NULL},
#endif
};

Int2 Main (void)

{
  Char         app [64];
  CharPtr      asnidx, directory, filter, infile, logfile, outfile, str, suffix;
  Boolean      automatic, batch, binary, catenated, compressed, dorecurse,
               indexed, local, lock, remote, summary, usethreads, bouncefound = FALSE;
#ifdef INTERNAL_NCBI_ASN2VAL
  Boolean      hup = FALSE;
#endif
  ErrSev       oldsev;
  ValNodePtr   parflat_list, vnp;
  time_t       run_time, start_time, stop_time;
  STimeout     timeout = { 0, 100000 };
  Int2         type = 0, val;
  ValFlagData  vfd;

  /* standard setup */

  ErrSetFatalLevel (SEV_MAX);
  ErrSetMessageLevel (SEV_MAX);
  ErrSetLogLevel (SEV_ERROR);
  ErrClearOptFlags (EO_SHOW_USERSTR);
  ErrSetLogfile ("stderr", ELOG_APPEND);
  ErrSetOpts (ERR_IGNORE, ERR_LOG_ON);

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

  parflat_list = Validate_ParFlat_GBFeat ();
  if (parflat_list != NULL) {
    Message (MSG_POSTERR, "Validate_ParFlat_GBFeat warnings");
    for (vnp = parflat_list; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      Message (MSG_POSTERR, "%s", str);
    }
    ValNodeFreeData (parflat_list);
  }

  /* process command line arguments */

  sprintf (app, "asnval %s", ASNVAL_APPLICATION);
  if (! GetArgs (app, sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  /* test network connection if requested */

  if ((Boolean) myargs [F_argTestNetwork].intvalue) {
    bouncefound = FALSE;
    start_time = GetSecs ();
    oldsev = ErrSetMessageLevel (SEV_FATAL);
    if (! NetTestAsynchronousQuery (&bouncequeue, BounceProc, (Pointer) &bouncefound)) {
      ErrSetMessageLevel (oldsev);
      Message (MSG_POSTERR, "NetTestAsynchronousQuery failed");
      return 1;
    }

    /* busy wait here, would normally call NetTestCheckQueue from event loop timer */
    while (! bouncefound) {
      stop_time = GetSecs ();
      if (stop_time - start_time >= 30) {
        Message (MSG_POSTERR, "Internet connection attempt timed out, exiting");
        return 1;
      }
      /* wait 0.1 seconds between attempts to avoid hogging machine */
      SOCK_Poll (0, 0, &timeout, 0);
      NetTestCheckQueue (&bouncequeue);
    }
    QUERY_CloseQueue (&bouncequeue);
    ErrSetMessageLevel (oldsev);

    Message (MSG_POSTERR, "Internet connection attempt succeeded, exiting");
    return 0;
  }

  /* additional setup modifications */

  MemSet ((Pointer) &vfd, 0, sizeof (ValFlagData));

  directory = (CharPtr) myargs [p_argInputPath].strvalue;
  suffix = (CharPtr) myargs [x_argSuffix].strvalue;
  filter = (CharPtr) myargs [f_argFilter].strvalue;
  infile = (CharPtr) myargs [i_argInputFile].strvalue;
  outfile = (CharPtr) myargs [o_argOutputFile].strvalue;
  dorecurse = (Boolean) myargs [u_argRecurse].intvalue;
  remote = (Boolean ) myargs [r_argRemote].intvalue;
  local = (Boolean) myargs [k_argLocalFetch].intvalue;
#ifdef INTERNAL_NCBI_ASN2VAL
  hup = (Boolean) myargs [H_argAccessHUP].intvalue;
#endif
  asnidx = (CharPtr) myargs [d_argAsnIdx].strvalue;
  indexed = (Boolean) StringDoesHaveText (asnidx);
  lock = (Boolean) myargs [l_argLockFar].intvalue;
  usethreads = (Boolean) myargs [T_argThreads].intvalue;

  vfd.severity = (Int2) myargs [R_argSeverity].intvalue;
  vfd.lowCutoff = (Int2) myargs [Q_argLowCutoff].intvalue;
  vfd.highCutoff = (Int2) myargs [P_argHighCutoff].intvalue;
  vfd.errcode = (CharPtr) myargs [E_argOnlyThisErr].strvalue;
  vfd.validateAlignments = (Boolean) myargs [A_argAlignments].intvalue;
  vfd.alignFindRemoteBsp = (Boolean) (vfd.validateAlignments && remote);
  vfd.doSeqHistAssembly = (Boolean) myargs [A_argAlignments].intvalue;
  vfd.farIDsInAlignments = (Boolean) myargs [A_argAlignments].intvalue;
  vfd.alwaysRequireIsoJTA = (Boolean) myargs [J_argIsoJta].intvalue;
  vfd.farFetchCDSproducts = (Boolean) myargs [Z_argRemoteCDS].intvalue;
  vfd.farFetchMRNAproducts = (Boolean) myargs [Z_argRemoteCDS].intvalue;
  vfd.locusTagGeneralMatch = (Boolean) myargs [M_argMatchTag].intvalue;
  vfd.validateIDSet = (Boolean) myargs [Y_argCheckOld].intvalue;
  vfd.ignoreExceptions = (Boolean) myargs [e_argIgnoreExcept].intvalue;
  vfd.validateExons = (Boolean) myargs [X_argExonSplice].intvalue;
  vfd.inferenceAccnCheck = (Boolean) myargs [G_argInfAccns].intvalue;
  vfd.disableSuppression = (Boolean) myargs [D_argDisableSuppress].intvalue;
  vfd.genomeSubmission = (Boolean) myargs [U_argGenomeSubmission].intvalue;
  vfd.debugTestDuJour = (Boolean) myargs [j_argDebugDuJour].intvalue;
  vfd.validateBarcode = (Boolean) myargs[B_argBarcodeVal].intvalue;
  vfd.taxfetch = (Boolean) myargs [q_argTaxLookup].intvalue;

  val = (Int2) myargs [N_argLatLonStrict].intvalue;
  vfd.testLatLonSubregion = (Boolean) ((val & LAT_LON_STATE) != 0);
  vfd.strictLatLonCountry = (Boolean) ((val & LAT_LON_STRICT) != 0);

  vfd.verbosity = (Int2) myargs [v_argVerbosity].intvalue;

  vfd.skipcount = (Int4) myargs [S_argSkipCount].intvalue;
  vfd.maxcount = (Int4) myargs [C_argMaxCount].intvalue;
  if (vfd.maxcount < 1) {
    vfd.maxcount = INT4_MAX;
  }

  vfd.seqSubmitParent = (Boolean) myargs [w_argSeqSubParent].intvalue;
  vfd.indexerVersion = (Boolean) myargs [y_argAIndexer].intvalue;

#ifdef INTERNAL_NCBI_ASN2VAL
  SetAppProperty ("InternalNcbiSequin", (void *) 1024);
#endif

  automatic = FALSE;
  catenated = FALSE;
  batch = FALSE;
  binary = (Boolean) myargs [b_argBinary].intvalue;
  compressed = (Boolean) myargs [c_argCompressed].intvalue;

  str = myargs [a_argType].strvalue;
  if (StringICmp (str, "a") == 0) {
    type = 1;
    automatic = TRUE;
  } else if (StringICmp (str, "c") == 0) {
    type = 1;
    catenated = TRUE;
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
    type = 4;
    batch = TRUE;
  } else if (StringICmp (str, "u") == 0) {
    type = 5;
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

  if (StringHasNoText (directory) && StringHasNoText (infile)) {
    Message (MSG_FATAL, "Input path or input file must be specified");
    return 1;
  }

  logfile = (CharPtr) myargs [L_argLogFile].strvalue;
  summary = (Boolean) myargs [K_argSummary].intvalue;

  start_time = GetSecs ();

  /* populate parameter structure */

  vfd.automatic = automatic;
  vfd.catenated = catenated;
  vfd.batch = batch;
  vfd.binary = binary;
  vfd.compressed = compressed;
  vfd.lock = lock;
  vfd.useThreads = usethreads;
  vfd.type = type;
  vfd.logfp = NULL;
  vfd.num_errors = 0;
  vfd.fatal_errors = 0;
  vfd.has_errors = FALSE;
  vfd.io_failure = FALSE;
  vfd.longest [0] = '\0';
  vfd.worsttime = 0;
  vfd.numrecords = 0;

  if (! StringHasNoText (outfile)) {
    if (vfd.verbosity == 0 || vfd.verbosity == 6) {
      vfd.outpath = outfile;
    } else {
      vfd.outfp = FileOpen (outfile, "w");
      if (vfd.outfp == NULL) {
        Message (MSG_FATAL, "Unable to open single output file");
        return 1;
      }
    }
  }

  if (! StringHasNoText (logfile)) {
    vfd.logfp = FileOpen (logfile, "w");
    if (vfd.logfp == NULL) {
      Message (MSG_FATAL, "Unable to open log file");
      return 1;
    }
  }

  /* register fetch functions */

  if (remote) {
#ifdef INTERNAL_NCBI_ASN2VAL
    if (hup) {
      DirSubFetchEnable ();
      SmartFetchEnable ();
      TPASmartFetchEnable ();
    }

    if (! PUBSEQBioseqFetchEnable ("asnval", FALSE)) {
      Message (MSG_POSTERR, "PUBSEQBioseqFetchEnable failed");
      return 1;
    }
    vfd.usePUBSEQ = TRUE;
    vfd.useThreads = FALSE;
#else
    PubSeqFetchEnable ();
#endif
    if (vfd.inferenceAccnCheck) {
      SeqMgrSetPreCache (GiRevHistLookupFarSeqIDs);
    }
    if (vfd.validateIDSet) {
      SeqMgrSetSeqIdSetFunc (GiRevHistLookupSeqIdSet);
    }
  }

  if (local) {
    LocalSeqFetchInit (FALSE);
  }

  if (indexed) {
    AsnIndexedLibFetchEnable (asnidx, TRUE);
  }

  /* recurse through all files within source directory or subdirectories */

  if (StringDoesHaveText (directory)) {

    DirExplore (directory, filter, suffix, dorecurse, ProcessOneRecord, (Pointer) &vfd);

  } else if (StringDoesHaveText (infile)) {

    ProcessOneRecord (infile, (Pointer) &vfd);
  }

  stop_time = GetSecs ();
  run_time = stop_time - start_time;

  if (vfd.outfp != NULL) {
    if (vfd.has_errors) {
      if (vfd.verbosity == 4) {
        fprintf (vfd.outfp, "</asnval>\n");
      }
      vfd.has_errors = FALSE;
    }
    if (summary) {
      fprintf (vfd.outfp, "Finished in %ld seconds\n", (long) run_time);
      if (StringDoesHaveText (vfd.longest)) {
        fprintf (vfd.outfp, "Longest processing time %ld seconds on %s\n",
                 (long) vfd.worsttime, vfd.longest);
        fprintf (vfd.outfp, "Total number of records %ld\n", (long) vfd.numrecords);
      }
    }
    FileClose (vfd.outfp);
  }

  if (vfd.logfp != NULL) {
    fprintf (vfd.logfp, "Finished in %ld seconds\n", (long) run_time);
    if (StringDoesHaveText (vfd.longest)) {
      fprintf (vfd.logfp, "Longest processing time %ld seconds on %s\n",
               (long) vfd.worsttime, vfd.longest);
      fprintf (vfd.logfp, "Total number of records %ld\n", (long) vfd.numrecords);
    }
    FileClose (vfd.logfp);
  }

  /* close fetch functions */

  if (indexed) {
    AsnIndexedLibFetchDisable ();
  }

  if (local) {
    LocalSeqFetchDisable ();
  }

  if (remote) {
#ifdef INTERNAL_NCBI_ASN2VAL
    PUBSEQBioseqFetchDisable ();
#else
    PubSeqFetchDisable ();
#endif
    SeqMgrSetPreCache (NULL);
    SeqMgrSetSeqIdSetFunc (NULL);
  }

  TransTableFreeAll ();

  ECNumberFSAFreeAll ();

  if (vfd.fatal_errors > 0) return 1;
  if (vfd.io_failure) return 1;

  return 0;
}

