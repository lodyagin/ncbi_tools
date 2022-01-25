/*   replicon.c
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
* File Name:  replicon.c
*
* Author:  Colleen Bollin
*
* Version Creation Date:   Feb. 1, 2012
*
* $Revision: 1.4 $
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
#include <sequtil.h>
#include <salutil.h>
#include <edutil.h>
#include <seqport.h>
#include <gather.h>
#include <sqnutils.h>
#include <subutil.h>
#include <toasn3.h>
#include <valid.h>
#include <asn2gnbk.h>
#include <explore.h>
#include <tofasta.h>
#include <simple.h>
#include <suggslp.h>
#include <toporg.h>
#include <aliparse.h>
#include <util/creaders/alnread.h>
#include <pmfapi.h>
#include <tax3api.h>
#ifdef INTERNAL_NCBI_TBL2ASN
#include <accpubseq.h>
#endif
#define NLM_GENERATED_CODE_PROTO
#include <objmacro.h>
#include <macroapi.h>

#define REPLICON_APP_VER "1.0"

CharPtr REPLICON_APPLICATION = REPLICON_APP_VER;



static void PopulateRepliconIdBuf (BioseqPtr bsp, CharPtr buf, Int4 buf_size)
{
  SeqIdPtr sip_local = NULL, sip_general = NULL, sip, sip_next;

  sip = bsp->id;
  while (sip != NULL && sip_local == NULL) {
    if (sip->choice == SEQID_LOCAL) {
      sip_local = sip;
    } else if (sip->choice == SEQID_GENERAL) {
      sip_general = sip;
    }
    sip = sip->next;
  }
  sip = NULL;
  if (sip_local != NULL) {
    sip = sip_local;
  } else if (sip_general != NULL) {
    sip = sip_general;
  } else {
    sip = SeqIdFindBest (bsp->id, SEQID_GENBANK);
  }

  sip_next = sip->next;
  sip->next = NULL;
  SeqIdWrite (sip, buf, PRINTID_FASTA_LONG, buf_size - 1);
  sip->next = sip_next;
}


typedef struct tablefiles {
  FILE *complete;
  FILE *incomplete;
  ValNodePtr chr_list;
} TableFilesData, PNTR TableFilesPtr;


static void MakeTable (BioseqPtr bsp, Pointer data)
{
  SeqMgrDescContext context;
  SeqDescPtr sdp;
  BioSourcePtr biop;
  Char         buf[PATH_MAX];
  CharPtr      chr_name = "ANONYMOUS";
  CharPtr      loc_str  = "UNKNOWN";
  CharPtr      type_str = "UNKNOWN";
  MolInfoPtr   mip;
  TableFilesPtr t;
  CharPtr      col3fmt = "%s\t%s\t%s\n";
  CharPtr      col3;
  Int4         len;

  if (bsp == NULL || ISA_aa(bsp->mol)) {
    return;
  }
  t = (TableFilesPtr) data;

  PopulateRepliconIdBuf (bsp, buf, sizeof (buf));

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
  if (sdp == NULL || (biop = (BioSourcePtr) sdp->data.ptrvalue) == NULL) {
    printf ("ERROR! No BioSource for %s\n", buf);
  } else {
    chr_name = GetRepliconChromosomeName (biop);
    loc_str = GetRepliconLocation (biop);
    type_str = GetRepliconType (biop);
    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &context);
    if (sdp == NULL || (mip = (MolInfoPtr) sdp->data.ptrvalue) == NULL || mip->completeness != 1) {
      /* not complete - looking for organelles */
      if (chr_name != NULL) {
        if (loc_str == NULL || type_str == NULL) {
          printf ("ERROR! Unrecognized BioSource.genome value!\n");
        } else {
          fprintf (t->incomplete, "%s\t%s\n", buf, chr_name);
          len = StringLen (col3fmt) + StringLen (chr_name) + StringLen (loc_str) + StringLen (type_str);
          col3 = (CharPtr) MemNew (sizeof (Char) * len);
          sprintf (col3, col3fmt, chr_name, loc_str, type_str);
          ValNodeAddPointer (&(t->chr_list), 0, col3);
        }
      }
    } else {
      /* complete */
      if (chr_name == NULL || loc_str == NULL || type_str == NULL) {
        printf ("ERROR! Unrecognized BioSource.genome value!\n");
      } else if (t != NULL && t->complete != NULL) {  
        fprintf (t->complete, "%s\t%s\t%s\t%s\n", buf, chr_name, loc_str, type_str);
      } else {
        printf ("%s\t%s\t%s\t%s\n", buf, chr_name, loc_str, type_str);
      }
    }
    chr_name = MemFree (chr_name);
    loc_str = MemFree (loc_str);
    type_str = MemFree (type_str);
  }
}
/* Args structure contains command-line arguments */ 
typedef enum {
  i_argInputFile = 0,
  c_argCompleteOuputFile,
  o_argIncompleteOrgFile,
  s_argIncompleteSeqFile
} Arguments;


Args myargs [] = {
  {"File List File", NULL, NULL, NULL,
    TRUE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
  {"Complete Output File", NULL, NULL, NULL,
    TRUE, 'c', ARG_FILE_OUT, 0.0, 0, NULL},
  {"Incomplete Org Output File", NULL, NULL, NULL,
    TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
  {"Incomplete Seq Output File", NULL, NULL, NULL,
    TRUE, 's', ARG_FILE_OUT, 0.0, 0, NULL},
};

Int2 Main (void)

{
  Char        app [64];
  CharPtr     input_file, complete_output_file, incomplete_org_file, incomplete_seq_file;
  FILE *      fp;
  FILE *      fi;
  Pointer     dataptr;
  Uint2       datatype;
  TableFilesData t;
  ValNodePtr     vnp;
  Int4           i;
  ReadBufferData rbd;
  CharPtr        line;

  /* standard setup */

  ErrSetFatalLevel (SEV_MAX);
  ErrSetMessageLevel (SEV_MAX);
  ErrClearOptFlags (EO_SHOW_USERSTR);
  UseLocalAsnloadDataAndErrMsg ();
  ErrPathReset ();

  /* finish resolving internal connections in ASN.1 parse tables */

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
  sprintf (app, "replicon %s", REPLICON_APPLICATION);
  if (! GetArgs (app, sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  input_file = (CharPtr) myargs [i_argInputFile].strvalue;
  complete_output_file = (CharPtr) myargs [c_argCompleteOuputFile].strvalue;
  incomplete_org_file = (CharPtr) myargs [o_argIncompleteOrgFile].strvalue;
  incomplete_seq_file = (CharPtr) myargs [s_argIncompleteSeqFile].strvalue;

  if (StringHasNoText (input_file)) {
    Message (MSG_FATAL, "Must supply input file.");
    return 1;
  }
  if (StringHasNoText (complete_output_file)) {
    Message (MSG_FATAL, "Must supply filename for complete replicons.");
    return 1;
  }
  if (StringHasNoText (incomplete_org_file)) {
    Message (MSG_FATAL, "Must supply filename for list of incomplete replicon sources.");
    return 1;
  }
  if (StringHasNoText (incomplete_seq_file)) {
    Message (MSG_FATAL, "Must supply filename for list of incomplete replicon sequences.");
    return 1;
  }
  

  t.complete = FileOpen (complete_output_file, "w");
  if (t.complete == NULL) {
    Message (MSG_FATAL, "Unable to open %s", complete_output_file);
    return 1;
  }
  t.incomplete = FileOpen (incomplete_seq_file, "w");
  if (t.incomplete == NULL) {
    Message (MSG_FATAL, "Unable to open %s", incomplete_seq_file);
    return 1;
  }
  t.chr_list = NULL;
  
  fi = FileOpen (input_file, "r");
  if (fi == NULL) {
    Message (MSG_FATAL, "Unable to open %s", input_file);
    return 1;
  }

  rbd.fp = fi;
  rbd.current_data = NULL;

  line = AbstractReadFunction (&rbd);
  while (line != NULL) {
    fp = FileOpen (line, "r");
    if (fp == NULL) {
      Message (MSG_FATAL, "Unable to open %s", line);
      return 1;
    }
    while ((dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, FALSE, FALSE, TRUE, FALSE)) != NULL) {
      switch (datatype) {
        case OBJ_SEQENTRY:
          VisitBioseqsInSep ((SeqEntryPtr) dataptr, &t, MakeTable);
          break;
        case OBJ_BIOSEQSET:
          VisitBioseqsInSet ((BioseqSetPtr) dataptr, &t, MakeTable);
          break;
        case OBJ_BIOSEQ:
          MakeTable ((BioseqPtr) dataptr, &t);
          break;
        default:
          Message (MSG_ERROR, "Unrecognized data type %d", datatype);
          break;
      }
      ObjMgrFree (datatype, dataptr);
    }
    FileClose (fp);
    line = AbstractReadFunction (&rbd);
  }
  FileClose (fi);

  FileClose (t.complete);
  FileClose (t.incomplete);
  fp = FileOpen (incomplete_org_file, "w");
  if (fp == NULL) {
    Message (MSG_FATAL, "Unable to open %s", incomplete_org_file);
    return 1;
  }
  t.chr_list = ValNodeSort (t.chr_list, SortVnpByString);
  ValNodeUnique (&(t.chr_list), SortVnpByString, ValNodeFreeData);
  for (vnp = t.chr_list; vnp != NULL; vnp = vnp->next) {
    fprintf (fp, "%s", (CharPtr)vnp->data.ptrvalue);
  }
  FileClose (fp);
  return 0;
}

