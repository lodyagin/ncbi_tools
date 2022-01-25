/*   aceread_tst.c
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
* File Name:  aceread_tst.c
*
* Author:  Colleen Bollin
*
* Version Creation Date:   7/22/08
*
* $Revision: 1.11 $
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
#include <pmfapi.h>
#ifdef INTERNAL_NCBI_ASNDISC
#include <accpubseq.h>
#include <tax3api.h>
#endif

#include "aceread.h"
#include "acerdapi.h"

typedef enum {
  i_argInputFile,
  o_argOutputFile,
  f_argFASTA,
  S_argIDSubstitutionFile,
  R_argSRRids,
  L_argSuppressIdLookup,
  Q_argMakeQualScores,
  X_argXMLFile,
  t_argTemplateFile,
  T_argTSAFields,
  C_argCenter,
  F_argFormat,
  G_argGapString,
  V_argValidateAgainstAsn1File,
  q_argReadQualScoresFile,
  r_argReadFASTAFile,
  N_argRecalculateConsensus,
  l_argLimitNumContigs
} EArgNum;

Args myargs [] = {
  {"Single Input File", "stdin", NULL, NULL,
    TRUE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
  {"Single Output File", NULL, NULL, NULL,
    TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
  {"FASTA Output", "F", NULL, NULL,
    TRUE, 'f', ARG_BOOLEAN, 0.0, 0, NULL},
  {"ID Substitution File", "", NULL, NULL,
    TRUE, 'S', ARG_FILE_IN, 0.0, 0, NULL},
  {"Replacement IDs are SRR", "F", NULL, NULL,
    TRUE, 'R', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Suppress ID Lookup", "F", NULL, NULL,
    TRUE, 'L', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Make Qual Scores", "T", NULL, NULL,
    TRUE, 'Q', ARG_BOOLEAN, 0.0, 0, NULL},
  {"XML Output File", "", NULL, NULL,
    TRUE, 'X', ARG_FILE_OUT, 0.0, 0, NULL },
  {"Template File", "", NULL, NULL,
    TRUE, 't', ARG_FILE_IN, 0.0, 0, NULL },
  {"TSA fields", NULL, NULL, NULL,
    TRUE, 'T', ARG_STRING, 0.0, 0, NULL },
  {"Genome Center Tag", NULL, NULL, NULL,
    TRUE, 'C', ARG_STRING, 0.0, 0, NULL},
  {"Assembly Format\n\tM MAQ\n\tE Standalone Eland\n\tA ACE", "A", NULL, NULL,
    TRUE, 'F', ARG_STRING, 0.0, 0, NULL},
  {"Gap String", NULL, NULL, NULL,
    TRUE, 'G', ARG_STRING, 0.0, 0, NULL},
  {"ASN.1 File to validate against", NULL, NULL, NULL,
    TRUE, 'V', ARG_FILE_IN, 0.0, 0, NULL},
  {"Quality score file for read sequences", NULL, NULL, NULL,
    TRUE, 'q', ARG_FILE_IN, 0.0, 0, NULL},
  {"FASTA file for read sequences (to use when trimming read quality scores)", NULL, NULL, NULL,
    TRUE, 'r', ARG_FILE_IN, 0.0, 0, NULL},
  {"Recalculate consensus sequence using read data\n\tW Whole Consensus\n\tN Ns Only", "", NULL, NULL,
    TRUE, 'N', ARG_STRING, 0.0, 0, NULL},
  {"Limit number of contigs to read", NULL, NULL, NULL,
    TRUE, 'l', ARG_INT, 0.0, 0, NULL},
};


static FILE *OpenAceFile (CharPtr infile)
{
  FILE        *f;
  Int4        len;
#ifdef OS_UNIX
  Char            cmmd [256];
  CharPtr         gzcatprog;
  int             ret;
  Boolean         usedPopen = FALSE;
#endif

  len = StringLen (infile);
  if (StringCmp (infile + len - 3, ".gz") == 0) {
#ifdef OS_UNIX
    gzcatprog = getenv ("NCBI_UNCOMPRESS_BINARY");
    if (gzcatprog != NULL) {
      sprintf (cmmd, "%s %s", gzcatprog, infile);
    } else {
      ret = system ("gzcat -h >/dev/null 2>&1");
      if (ret == 0) {
        sprintf (cmmd, "gzcat %s", infile);
      } else if (ret == -1) {
        Message (MSG_POSTERR, "Unable to fork or exec gzcat in ScanBioseqSetRelease");
        return NULL;
      } else {
        ret = system ("zcat -h >/dev/null 2>&1");
        if (ret == 0) {
          sprintf (cmmd, "zcat %s", infile);
        } else if (ret == -1) {
          Message (MSG_POSTERR, "Unable to fork or exec zcat in ScanBioseqSetRelease");
          return NULL;
        } else {
          Message (MSG_POSTERR, "Unable to find zcat or gzcat in ScanBioseqSetRelease - please edit your PATH environment variable");
          return NULL;
        }
      }
    }
    f = popen (cmmd, "r");
    usedPopen = TRUE;
#else
    Message (MSG_POSTERR, "Unable to read gzipped files when not running in UNIX");
    return NULL;
#endif
  } else {
    f = FileOpen (infile, "r");
  }
  return f;
}


static Boolean ValidateAgainstASNFile (TACEFilePtr ace_file, CharPtr filename, char *has_errors)
{
  Pointer      dataptr;
  Uint2        datatype;
  SeqEntryPtr  sep = NULL;
  SeqSubmitPtr ssp = NULL;
  Boolean      chars_stripped = FALSE;
  FILE *fp;
  Boolean      rval = FALSE;
  

  fp = FileOpen (filename, "r");
  if (fp == NULL) {
    printf ("Unable to open %s\n", filename);
    return FALSE;
  }

  /* Read in one sequence from the file */
  dataptr = ReadAsnFastaOrFlatFileEx (fp, &datatype, NULL, FALSE, FALSE,
		                   	                  TRUE, FALSE, &chars_stripped);      
  FileClose (fp);
  if (NULL == dataptr) 
  {
    printf ("Unable to read SeqEntry from %s\n", filename);
    return FALSE;
  }

  /* Convert the file data to a SeqEntry */
  
  if (datatype == OBJ_SEQENTRY)
    sep = (SeqEntryPtr) dataptr;
  else if (datatype == OBJ_BIOSEQ || datatype == OBJ_BIOSEQSET)
    sep = SeqMgrGetSeqEntryForData (dataptr);
  else if (datatype == OBJ_SEQSUB) 
  {
    ssp = (SeqSubmitPtr) dataptr;
    if (ssp != NULL && ssp->datatype == 1)
    {
      sep = (SeqEntryPtr) ssp->data;
    }
  }
  
  rval = ValidateACEFileAgainstSeqEntry (ace_file, sep, has_errors);

  if (ssp != NULL) {
    ssp = SeqSubmitFree (ssp);
  } else {
    sep = SeqEntryFree (sep);
  }
  return rval;
 
}


static Boolean StringNHasNoText (CharPtr str, Int4 n)
{
  CharPtr cp;
  Int4    i;
  if (str == NULL) return TRUE;
  cp = str;
  i = 0;
  while (i < n) {
    if (*cp == 0) return TRUE;
    if (!isspace (*cp)) return FALSE;
    cp++;
    i++;
  }
  return TRUE;
}


static Boolean BracketMatchesLabel (CharPtr cp, CharPtr cp_equal, CharPtr label) 
{
  Int4 len;

  if (cp == NULL || cp_equal == NULL || label == NULL) return FALSE;

  len = StringLen (label);
  if (StringNCmp (cp, label, len) == 0 
        && StringNHasNoText (cp + len, cp_equal - cp - len)) {
    return TRUE;
  } else {
    return FALSE;
  }
}


static CharPtr GetBracketValue (CharPtr cp, CharPtr cp_end)
{
  Int4 len;
  CharPtr val = NULL;

  if (cp == NULL || cp_end == NULL || cp_end <= cp) return NULL;

  cp += StringSpn (cp, " \t");
  len = (cp_end - cp) + 1;
  val = (CharPtr) MemNew (sizeof (Char) * len);
  StringNCpy (val, cp, len - 1); 
  val [len] = 0;
  while (len > 1 && isspace (val [len-1])) {
    len--;
    val[len] = 0;
  }
  return val;
}


static Boolean
GetTSAFieldsFromString
(CharPtr str,
 CharPtr PNTR p_submitter_reference,
 CharPtr PNTR p_archive_id,
 CharPtr PNTR p_description)
{
  CharPtr cp, cp_next, cp_equal, cp_end;
  CharPtr subref = NULL, arch_id = NULL, desc = NULL;
  Boolean is_bad = FALSE;

  if (p_submitter_reference != NULL) {
    *p_submitter_reference = NULL;
  }
  if (p_archive_id != NULL) {
    *p_archive_id = NULL;
  }
  if (p_submitter_reference != NULL) {
    *p_description = NULL;
  }
  if (StringHasNoText (str)) {
    return TRUE;
  }

  cp = StringChr (str, '[');
  while (cp != NULL && !is_bad) {
    cp++;
    cp_next = StringChr (cp + 1, '[');
    cp_equal = StringChr (cp, '=');
    cp_end = StringChr (cp, ']');
    if (cp_equal == NULL || cp_end == NULL) {
      is_bad = TRUE;
    } else if (cp_equal > cp_end) {
      is_bad = TRUE;
    } else if (cp_next != NULL && (cp_equal > cp_next || cp_end > cp_next)) {
      is_bad = TRUE;
    } else {
      cp += StringSpn (cp, " \t");
      if (BracketMatchesLabel (cp, cp_equal, "subref")) {
        if (subref == NULL) {
          subref = GetBracketValue (cp_equal + 1, cp_end);
        } else {
          is_bad = TRUE;
        }
      } else if (BracketMatchesLabel (cp, cp_equal, "archive_id")) {
        if (arch_id == NULL) {
          arch_id = GetBracketValue (cp_equal + 1, cp_end);
        } else {
          is_bad = TRUE;
        }
      } else if (BracketMatchesLabel (cp, cp_equal, "desc")) {
        if (desc == NULL) {
          desc = GetBracketValue (cp_equal + 1, cp_end);
        } else {
          is_bad = TRUE;
        }
      } else {
        is_bad = TRUE;
      }
    }
    cp = cp_next;
  }
  if (p_submitter_reference == NULL) {
    subref = MemFree (subref);
  } else {
    *p_submitter_reference = subref;
  }
  if (p_archive_id == NULL) {
    arch_id = MemFree (arch_id);
  } else {
    *p_archive_id = arch_id;
  }
  if (p_description == NULL) {
    desc = MemFree (desc);
  } else {
    *p_description = desc;
  }
  return TRUE;
}


static void PrintTraceGapsXML (TGapInfoPtr gap_info)
{
  Int4 i;

  if (gap_info != NULL) {
    printf ("    <ntracegaps>%d</ntracegaps>\n", gap_info->num_gaps);
    if (gap_info->num_gaps > 0) {
      printf ("      <tracegaps source=\"INLINE\">");
      for (i = 0; i < gap_info->num_gaps - 1; i++) {
        printf ("%d,", gap_info->gap_offsets[i]);
      }
      printf ("%d</tracegaps>\n", gap_info->gap_offsets[gap_info->num_gaps - 1]);
    }
  }
}


static void TestPosConversions (TGapInfoPtr gap_info)
{
  Int4 i, t_pos, s_pos = 0, r_pos;
  Int4 test_len = 0;

  if (gap_info != NULL && gap_info->num_gaps > 0) {
    for (i = 0; i < gap_info->num_gaps; i++) {
      test_len += gap_info->gap_offsets[i] + 1;
    }
    for (i = 0; i < test_len; i++) {
      s_pos = SeqPosFromTilingPos (i, gap_info);
      t_pos = TilingPosFromSeqPos (s_pos, gap_info);
      if (t_pos != i) {
        printf ("Failed!  %d -> SeqPosFromTilingPos -> %d -> TilingPosFromSeqPos -> %d\n",
                i, s_pos, t_pos);
      }
      r_pos = SeqPosFromTilingPos (t_pos, gap_info);
      if (r_pos != s_pos) {
        printf ("Failed!  %d -> TilingPosFromSeqPos -> %d -> SeqPosFromTilingPos -> %d\n",
                s_pos, t_pos, r_pos);
      }
      /* printf ("%d:%d:%d:%d\n", i, s_pos, t_pos, r_pos); */
    }
  }
}


static void PrintTraceReadXML (TContigReadPtr read)
{
  if (read == NULL) {
    printf ("Bad read\n");
  } else {
    printf ("<trace>\n");
    printf ("  <trace_name>%s</trace_name>\n", read->read_id == NULL ? "" : read->read_id);
    PrintTraceGapsXML (read->gaps);
    printf ("  <nbasecalls>%d</nbasecalls>\n", StringLen (read->read_seq));
    printf ("  <valid>\n");
    printf ("    <start>%d</start>\n", read->read_assem_start + 1);
    printf ("    <stop>%d</stop>\n", read->read_assem_stop + 1);
    printf ("  </valid>\n");
    printf ("  <tiling direction = \"%s\">\n", read->is_complement ? "REVERSE" : "FORWARD");
    printf ("    <start>%d</start>\n", read->cons_start + 1);
    printf ("    <start>%d</start>\n", read->cons_start + StringLen (read->read_seq) + 1);
    printf ("  </tiling>\n");
    printf ("  <consensus>\n");
    printf ("    <start>%d</start>\n", read->cons_start + 1);
    printf ("    <start>%d</start>\n", read->cons_start + StringLen (read->read_seq) + 1);
    printf ("  </consensus>\n");
    printf ("<trace>\n");
  }
}



static void TestGapInfoReading (CharPtr gap_string)
{
  TGapInfoPtr  gap_info;
  ValNodePtr   list, vnp;
  
  if (!StringHasNoText (gap_string)) {
    gap_info = GapInfoFromSequenceString(gap_string, "*");
    if (gap_info == NULL) {
      printf ("error reading");
    } else {
      PrintTraceGapsXML (gap_info);
      TestPosConversions (gap_info);
      list = GetTransitionsFromGapInfo (gap_info, 0, 0, 40);
      for (vnp = list; vnp != NULL; vnp = vnp->next) {
        printf ("%d\n", vnp->data.intvalue);
      }
    }
    GapInfoFree (gap_info);
  }
}


static void AddAlignmentToSeqEntry (DenseSegPtr dsp, SeqEntryPtr sep)
{
  SeqAnnotPtr  sap;
  SeqAlignPtr  salp;
  BioseqPtr    bsp;
  BioseqSetPtr bssp;

  if (dsp == NULL || sep == NULL) return;

  sap = SeqAnnotNew ();
  sap->type = 2;

  salp = SeqAlignNew ();
  salp->type = 3;
  salp->segtype = 2;
  salp->segs = (Pointer) dsp;
  salp->dim = dsp->dim;
  sap->data = (Pointer) salp;

  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap->next = bsp->annot;
    bsp->annot = sap;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap->next = bssp->annot;
    bssp->annot = sap;
  }
}


static void AddDescrToNucBioseqCallback (BioseqPtr bsp, Pointer data)
{
  SeqDescrPtr sdp, sdp_copy;

  if (bsp == NULL || !ISA_na (bsp->mol) || data == NULL) { 
    return;
  }
  sdp = (SeqDescrPtr) data;
  sdp_copy = (SeqDescrPtr) AsnIoMemCopy (sdp, (AsnReadFunc) SeqDescrAsnRead, (AsnWriteFunc) SeqDescrAsnWrite);
  sdp_copy->next = bsp->descr;
  bsp->descr = sdp_copy;
}

  
static SeqSubmitPtr AddSeqSubmitFromTemplate (SeqEntryPtr sep, CharPtr filename)
{
  SeqSubmitPtr   ssp = NULL;
  SubmitBlockPtr sbp;
  CitSubPtr      csp;
  FILE *fp = NULL;
  Pointer         dataptr;
  Uint2           datatype;

  if (StringHasNoText (filename)) {
    return NULL;
  }
    
  fp = FileOpen (filename, "r");
  if (fp == NULL) {
    printf ("Unable to read template file %s\n", filename);
    return NULL;
  }

  while ((dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, FALSE, FALSE, TRUE, FALSE)) != NULL) {
    if (datatype == OBJ_SEQSUB) {
      ssp = (SeqSubmitPtr) dataptr;
      ssp->datatype = 1;
      ssp->data = sep;
    } else if (datatype == OBJ_SUBMIT_BLOCK) {
      sbp = (SubmitBlockPtr) dataptr;
      ssp = SeqSubmitNew ();
      ssp->datatype = 1;
      ssp->data = sep;
      ssp->sub = sbp;
    } else if (datatype == OBJ_SEQDESC) {
      VisitBioseqsInSep (sep, dataptr, AddDescrToNucBioseqCallback);
      ObjMgrFree (datatype, dataptr);
    } else {
      ObjMgrFree (datatype, dataptr);
    }
  }
  FileClose (fp);
  if (ssp == NULL) {
    ssp = SeqSubmitNew ();
    ssp->datatype = 1;
    ssp->data = sep;
  }

  if (ssp->sub == NULL) {
    ssp->sub = SubmitBlockNew ();
  } 

  ssp->sub->tool = MemFree (ssp->sub->tool);
  ssp->sub->tool = StringSave ("aceread");
  ssp->sub->hup = FALSE;
  ssp->sub->reldate = DateFree (ssp->sub->reldate);
  csp = ssp->sub->cit;
  if (csp != NULL) {
    csp->date = DateFree (csp->date);
    csp->date = DateCurr ();
  }
  return ssp;
}


static Boolean AddReadQualityScores (TACEFilePtr afp, CharPtr qs_filename, CharPtr rd_filename)
{
  ReadBufferData q, r;
  Boolean use_fasta = FALSE;
  Boolean rval = FALSE;

  if (afp == NULL || StringHasNoText (qs_filename)) {
    return TRUE;
  }

  q.current_data = NULL;
  r.current_data = NULL;

  q.fp = FileOpen (qs_filename, "r");
  if (q.fp == NULL) {
    printf ("Unable to read quality score file\n");
    return FALSE;
  }

  if (!StringHasNoText (rd_filename)) {
    r.fp = FileOpen (rd_filename, "r");
    if (r.fp == NULL) {
      printf ("Unable to open read FASTA file\n");
      FileClose (q.fp);
      return FALSE;
    }
    use_fasta = TRUE;
  }

  if (AddReadQualScores (afp, AbstractReadFunction, &q, use_fasta ? AbstractReadFunction : NULL, &r) > 0) {
    rval = TRUE;
  }

  FileClose (q.fp);
  if (use_fasta) {
    FileClose (r.fp);
  }
  return rval;
}


Int2 Main (void)

{
  CharPtr      infile, outfile, xmlfile;

  ReadBufferData    rbd;
  TACEFilePtr afp;
  Int4        i, len;
  SeqEntryPtr sep;
  AsnIoPtr    aip;
  FILE *f = NULL;
  FILE *f2;
  CharPtr app = "aceread_tst";
  BioseqSetPtr bssp;
  SeqEntryPtr  last_sep = NULL;
  Uint2        entityID;
  Boolean      make_qual_scores, suppress_lookup, srr_ids, fasta_out;
  CharPtr      submitter_ref = NULL, archive_id = NULL, description = NULL;
  CharPtr      center_name = NULL;
  CharPtr      format = NULL;
  CharPtr      gap_string;
  CharPtr      asn_file = NULL;
  Int4         limit = 0;
  char         has_errors = 0;
  Boolean      recalculate_consensus = FALSE, recalculate_only_Ns = FALSE;
  CharPtr      recalculate_options;
  SeqSubmitPtr ssp;

  /* standard setup */

  ErrSetFatalLevel (SEV_MAX);
  ErrSetMessageLevel (SEV_MAX);
  ErrClearOptFlags (EO_SHOW_USERSTR);
  ErrSetLogfile ("stderr", ELOG_APPEND);
  ErrSetOpts (ERR_IGNORE, ERR_LOG_ON);

  UseLocalAsnloadDataAndErrMsg ();
  ErrPathReset ();

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
  PubSeqFetchEnable ();

  if (! GetArgs (app, sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  recalculate_options = (CharPtr) myargs[N_argRecalculateConsensus].strvalue;
  if (!StringHasNoText (recalculate_options)) {
    if (StringCmp (recalculate_options, "W") == 0) {
      recalculate_consensus = TRUE;
      recalculate_only_Ns = FALSE;
    } else if (StringCmp (recalculate_options, "N") == 0) {
      recalculate_consensus = TRUE;
      recalculate_only_Ns = TRUE;
    } else {
      Message (MSG_FATAL, "Invalid consensus sequence recalculation option");
      return 1;
    }
  }


  /* test gap info reading if provided */
  gap_string = (CharPtr) myargs[G_argGapString].strvalue;
  TestGapInfoReading (gap_string);

  /* limit number of contigs?  for debugging purposes */
  limit = myargs[l_argLimitNumContigs].intvalue;

  /* select format of input file */
  format = (CharPtr) myargs[F_argFormat].strvalue;
  if (StringHasNoText (format)) {
    format = "A";
  }

  infile = (CharPtr) myargs [i_argInputFile].strvalue;
  if (StringHasNoText (infile)) {
    Message (MSG_FATAL, "Must supply input file!");
    return 1;
  }
  outfile = (CharPtr) myargs [o_argOutputFile].strvalue;
  xmlfile = (CharPtr) myargs[X_argXMLFile].strvalue;
  make_qual_scores = (Boolean) myargs [Q_argMakeQualScores].intvalue;
  center_name = (CharPtr) myargs[C_argCenter].strvalue;
  suppress_lookup = (Boolean) myargs [L_argSuppressIdLookup].intvalue;
  srr_ids = (Boolean) myargs[R_argSRRids].intvalue;
  fasta_out = (Boolean) myargs[f_argFASTA].intvalue;

  /* ASN.1 file to validate against */
  asn_file = (CharPtr) myargs [V_argValidateAgainstAsn1File].strvalue;

  if (!GetTSAFieldsFromString ((CharPtr) myargs [T_argTSAFields].strvalue,
                               &submitter_ref,
                               &archive_id,
                               &description)) {
    Message (MSG_FATAL, "Error reading TSA fields");
    return 1;
  }

  len = StringLen (infile);
  if (StringHasNoText (outfile)) {
    if (len > 3 && StringCmp (infile + len - 4, ".ace") == 0) {
      outfile = StringSave (infile);
      StringCpy (outfile + len - 3, "sqn");
    } else if (len > 6 && StringCmp (infile + len - 7, ".ace.gz") == 0) {
      outfile = StringSave (infile);
      StringCpy (outfile + len - 6, "sqn");
    } else {
      outfile = (CharPtr) MemNew (sizeof (Char) * (len + 5));
      sprintf (outfile, "%s.sqn", infile);
    }
  }

  if (!StringHasNoText ((CharPtr) myargs [S_argIDSubstitutionFile].strvalue)) {
    f = FileOpen (myargs [S_argIDSubstitutionFile].strvalue, "r");
    if (f == NULL) {
      Message (MSG_FATAL, "Unable to open %s", myargs [S_argIDSubstitutionFile].strvalue);
      return 1;
    }
  }

  if (StringChr (format, 'M') != NULL) {
    rbd.fp = FileOpen (infile, "r");
    if (rbd.fp == NULL) {
      Message (MSG_FATAL, "Unable to open %s", infile);
      return 1;
    }

    rbd.current_data = NULL;
    afp = ReadMAQFile (AbstractReadFunction, &rbd);
  } else if (StringChr (format, 'E') != NULL) {
    rbd.fp = FileOpen (infile, "r");
    if (rbd.fp == NULL) {
      Message (MSG_FATAL, "Unable to open %s", infile);
      return 1;
    }

    rbd.current_data = NULL;
    afp = ReadElandStandaloneFile (AbstractReadFunction, &rbd);
  } else if (StringChr (format, 'A') != NULL) { 
    rbd.fp = OpenAceFile (infile);
    if (rbd.fp == NULL) {
      Message (MSG_FATAL, "Unable to open %s", infile);
      return 1;
    }
    rbd.current_data = NULL;
    afp = ReadACEFile ( AbstractReadFunction, &rbd, make_qual_scores, &has_errors);
  } else {
    Message (MSG_FATAL, "Unrecognized format: %s\n", format);
    return 1;
  }
  FileClose (rbd.fp);
  if (afp == NULL) {
    printf ("<message severity=\"ERROR\" seq-id=\"No ID\" code=\"bad_format\">Unable to read file</message>\n");
  } else {
    if (recalculate_consensus) {
        if (!AddReadQualityScores (afp, (CharPtr) myargs [q_argReadQualScoresFile].strvalue, (CharPtr) myargs [r_argReadFASTAFile].strvalue)) {
            printf ("<message severity=\"ERROR\" seq-id=\"No ID\" code=\"bad_format\">Failed to add read quality scores</message>\n");
        } else {
            RecalculateConsensusSequences (afp, recalculate_only_Ns);
        }
    }

    if (limit > 0) {
      for (i = limit; i < afp->num_contigs; i++) {
        ContigFree (afp->contigs[i]);
        afp->contigs[i] = NULL;
      }
      afp->num_contigs = limit;
    }

    if (f != NULL) {
      UpdateAceFileIds (afp, f, suppress_lookup, srr_ids, &has_errors);
      FileClose (f);
      f = NULL;
    }
    ValidateAceFileIds (afp, &has_errors);

    if (asn_file != NULL) {
      if (ValidateAgainstASNFile (afp, asn_file, &has_errors)) {
        printf ("Validation against %s succeeded\n", asn_file);
      }
    }

    if (!StringHasNoText (xmlfile)) {
      f2 = FileOpen (xmlfile, "w");
      WriteTraceAssemblyFromAceFile (afp, submitter_ref, center_name, 0, description, f2);
      FileClose (f2);
    }

    if (fasta_out) {
      f2 = FileOpen (outfile, "w");
      WriteFASTAFromAceFile (afp, f2);
      FileClose (f2);
    } else {
      aip = AsnIoOpen (outfile, "w");
      if (aip == NULL) {
        printf ("Unable to open %s\n", outfile);
      } else {
        bssp = BioseqSetNew ();
        bssp->_class = BioseqseqSet_class_genbank;

        for (i = 0; i < afp->num_contigs; i++) {
          sep = MakeSeqEntryFromContig (afp->contigs[i]);
          if (last_sep == NULL) {
            bssp->seq_set = sep;
          } else {
            last_sep->next = sep;
          }
          last_sep = sep;          
        }
        sep = ValNodeNew (NULL);
        sep->choice = 2;
        sep->data.ptrvalue = bssp;
        bssp->seqentry = sep;
        SeqMgrLinkSeqEntry (sep, 0, NULL);
        entityID = ObjMgrGetEntityIDForChoice (sep);
        AssignIDsInEntityEx (entityID, 0, NULL, NULL);
        SeqMgrIndexFeatures (entityID, sep);
        ssp = AddSeqSubmitFromTemplate (sep, (CharPtr) myargs[t_argTemplateFile].strvalue);
        if (ssp == NULL) {
          SeqEntryAsnWrite (sep, aip, NULL);
          sep = SeqEntryFree (sep);
        } else {
          SeqSubmitAsnWrite (ssp, aip, NULL);
          ssp = SeqSubmitFree (ssp);
        }
        AsnIoClose (aip);
      }
    }
  }

  if (has_errors) {
    printf ("</aceread>\n");
  }

  return 0;

}

