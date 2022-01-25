/*   export_cmt.c
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
* File Name:  export_cmt.c
*
* Author:  Colleen Bollin
*
* Version Creation Date:   Feb. 1, 2012
*
* $Revision: 1.1 $
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

#define EXPORT_APP_VER "1.0"

CharPtr EXPORT_APPLICATION = EXPORT_APP_VER;


static Boolean IsGenomeAssembly (UserObjectPtr uop)
{
  UserFieldPtr  ufp;

  if (uop == NULL || uop->type == NULL || StringICmp (uop->type->str, "StructuredComment") != 0) {
    return FALSE;
  }

  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    if (ufp->label != NULL) {
      if (StringICmp (ufp->label->str, "StructuredCommentPrefix") == 0) {
        if (ufp->choice == 1 && StringICmp (ufp->data.ptrvalue, "##Genome-Assembly-Data-START##") == 0) {
          return TRUE;
        }
      }
    }
  }
  return FALSE;
}


static void MakeTable (BioseqPtr bsp, Pointer data)
{
  SeqMgrDescContext context;
  SeqDescPtr sdp;
  UserObjectPtr uop;
  UserFieldPtr  ufp;
  SeqIdPtr      sip, sip_next;
  Char          buf[PATH_MAX];

  if (bsp == NULL || ISA_aa(bsp->mol)) {
    return;
  }

  sip = bsp->id;
  while (sip != NULL && sip->choice != SEQID_LOCAL) {
    sip = sip->next;
  }
  if (sip == NULL) {
    sip = SeqIdFindBest (bsp->id, SEQID_GENBANK);
  }

  sip_next = sip->next;
  sip->next = NULL;
  SeqIdWrite (sip, buf, PRINTID_FASTA_LONG, sizeof (buf) - 1);
  sip->next = sip_next;

  for (sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &context);
       sdp != NULL;
       sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &context)) {
    uop = (UserObjectPtr) sdp->data.ptrvalue;
    if (IsGenomeAssembly (uop)) {
      for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
        if (ufp->label != NULL) {
          if (StringICmp (ufp->label->str, "StructuredCommentPrefix") == 0) {
            /* ignore */
          } else if (StringICmp (ufp->label->str, "StructuredCommentSuffix") == 0) {
            /* ignore */
          } else {
            printf ("%s\t%s\n", ufp->label->str == NULL ? "" : (CharPtr) ufp->label->str, ufp->data.ptrvalue == NULL ? "" : (CharPtr) ufp->data.ptrvalue);
          }
        }
      }
    }
  }

}


/* Args structure contains command-line arguments */

typedef enum {
  i_argInputFile = 0,
} Arguments;


Args myargs [] = {
  {"Single Input File", NULL, NULL, NULL,
    TRUE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
};

Int2 Main (void)

{
  Char         app [64];
  CharPtr      file;
  FILE *       fp;
  Pointer      dataptr;
  Uint2        datatype;
  SeqSubmitPtr ssp;
  SeqEntryPtr  sep;

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

  sprintf (app, "export_cmt %s", EXPORT_APPLICATION);
  if (! GetArgs (app, sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  file = (CharPtr) myargs [i_argInputFile].strvalue;

  if (StringHasNoText (file)) {
    Message (MSG_FATAL, "Must supply input file.");
    return 1;
  }
  
  fp = FileOpen (file, "r");
  if (fp == NULL) {
    Message (MSG_FATAL, "Unable to open %s", file);
    return 1;
  }
  while ((dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, FALSE, FALSE, TRUE, FALSE)) != NULL) {
    switch (datatype) {
      case OBJ_SEQENTRY:
        VisitBioseqsInSep ((SeqEntryPtr) dataptr, NULL, MakeTable);
        break;
      case OBJ_BIOSEQSET:
        VisitBioseqsInSet ((BioseqSetPtr) dataptr, NULL, MakeTable);
        break;
      case OBJ_BIOSEQ:
        MakeTable ((BioseqPtr) dataptr, NULL);
        break;
      case OBJ_SEQSUB:
        ssp = (SeqSubmitPtr) dataptr;
        if (ssp != NULL && ssp->datatype == 1) {
          for (sep = ssp->data; sep != NULL; sep = sep->next) {
            VisitBioseqsInSep (sep, NULL, MakeTable);
          }
        }
        break;
      default:
        Message (MSG_ERROR, "Unrecognized data type %d", datatype);
        break;
    }
    ObjMgrFree (datatype, dataptr);
  }
  return 0;
}

