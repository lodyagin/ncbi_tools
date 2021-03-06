/*   asn2gnb3.c
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
* File Name:  asn2gnb3.c
*
* Author:  Karl Sirotkin, Tom Madden, Tatiana Tatusov, Jonathan Kans,
*          Mati Shomrat
*
* Version Creation Date:   10/21/98
*
* $Revision: 1.238 $
*
* File Description:  New GenBank flatfile generator - work in progress
*
* Modifications:
* --------------------------------------------------------------------------
* ==========================================================================
*/

#include <ncbi.h>
#include <objall.h>
#include <objsset.h>
#include <objsub.h>
#include <objfdef.h>
#include <objpubme.h>
#include <seqport.h>
#include <sequtil.h>
#include <sqnutils.h>
#include <subutil.h>
#include <tofasta.h>
#include <explore.h>
#include <gbfeat.h>
#include <gbftdef.h>
#include <edutil.h>
#include <alignmgr2.h>
#include <asn2gnbi.h>

#ifdef WIN_MAC
#if __profile__
#include <Profiler.h>
#endif
#endif

static CharPtr ref_link = "https://www.ncbi.nlm.nih.gov/RefSeq/";

static CharPtr doc_link = "https://www.ncbi.nlm.nih.gov/genome/annotation_euk/process/";

static CharPtr ev_link = "https://www.ncbi.nlm.nih.gov/sutils/evv.cgi?";

static CharPtr link_encode = "https://www.genome.gov/10005107";

static CharPtr link_seqn = "https://www.ncbi.nlm.nih.gov/nuccore/";
static CharPtr link_seqp = "https://www.ncbi.nlm.nih.gov/protein/";


/* ********************************************************************** */

static void AddHistCommentString (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  CharPtr prefix,
  CharPtr suffix,
  DatePtr dp,
  SeqIdPtr ids,
  Boolean is_na,
  Boolean use_accn
)

{
  Int2      count = 0;
  Char      buf [256], id [42];
  Boolean   first, skip;
  BIG_ID    gi = 0;
  SeqIdPtr  sip, sip2;
  CharPtr   strd;
  
  if (dp == NULL || ids == NULL || prefix == NULL || suffix == NULL || ffstring == NULL) return;

  strd = asn2gb_PrintDate (dp);
  if (strd == NULL) {
    strd = StringSave ("?");
  }

  for (sip = ids; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_GI) {
      gi = (BIG_ID) sip->data.intvalue;
      count++;
    }
  }

  if (count > 1) {
    sprintf (buf, "%s or before %s %s", prefix, strd, suffix);
  } else {
    sprintf (buf, "%s %s %s", prefix, strd, suffix);
  }
  FFAddOneString (ffstring, buf, FALSE, FALSE, TILDE_EXPAND);

  MemFree (strd);

  if (gi == 0) {
    FFAddOneString (ffstring, " gi:?", FALSE, FALSE, TILDE_EXPAND);
    return;
  }

  first = TRUE;
  for (sip = ids; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_GI) {
      gi = (BIG_ID) sip->data.intvalue;
      if (! first) {
        FFAddOneString (ffstring, ",", FALSE, FALSE, TILDE_IGNORE);
      }
      first = FALSE;
      skip = FALSE;
      if (use_accn) {
        sip2 = GetSeqIdForGI (gi);
        if (sip2 != NULL) {
          SeqIdWrite (sip2, id, PRINTID_TEXTID_ACC_VER, sizeof (id) - 1);
          if (StringDoesHaveText (id)) {
            if ( GetWWW(ajp) ) {
              FFAddOneString (ffstring, " ", FALSE, FALSE, TILDE_IGNORE);
              FFAddOneString (ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
              if (is_na) {
                FF_Add_NCBI_Base_URL (ffstring, link_seqn);
              } else {
                FF_Add_NCBI_Base_URL (ffstring, link_seqp);
              }
              sprintf (buf, "%ld", (long) gi);
              FFAddTextToString (ffstring, /* "val=" */ NULL, buf, "\">", FALSE, FALSE, TILDE_IGNORE);
              FFAddOneString (ffstring, id, FALSE, FALSE, TILDE_EXPAND);
              FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
            } else {
              sprintf (buf, " %s", id);
              FFAddOneString (ffstring, buf, FALSE, FALSE, TILDE_EXPAND);
            }
            skip = TRUE;
          }
          SeqIdFree (sip2);
        }
      }
      if (! skip) {
        if ( GetWWW(ajp) ) {
          FFAddOneString (ffstring, " gi:", FALSE, FALSE, TILDE_IGNORE);
          FFAddOneString (ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
          if (is_na) {
            FF_Add_NCBI_Base_URL (ffstring, link_seqn);
          } else {
            FF_Add_NCBI_Base_URL (ffstring, link_seqp);
          }
          sprintf (buf, "%ld", (long) gi);
          FFAddTextToString (ffstring, /* "val=" */ NULL, buf, "\">", FALSE, FALSE, TILDE_IGNORE);
          FFAddOneString (ffstring, buf, FALSE, FALSE, TILDE_EXPAND);
          FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
        } else {
          sprintf (buf, " gi:%ld", (long) gi);
          FFAddOneString (ffstring, buf, FALSE, FALSE, TILDE_EXPAND);
        }
      }
    }
  }

  FFAddOneString (ffstring, ".", FALSE, FALSE, TILDE_EXPAND);
}

static void AddUnorderedCommentString (
  StringItemPtr ffstring,
  BioseqPtr bsp
)

{
  Char         buffer [256];
  DeltaSeqPtr  dsp;
  ValNodePtr   head = NULL;
  Int4         num_gaps = 0;
  SeqLitPtr    slitp;
  SeqLocPtr    slocp;
  CharPtr      str;

  if (bsp == NULL) return;

  if (bsp->repr == Seq_repr_delta) {
    for (dsp = (DeltaSeqPtr) bsp->seq_ext; dsp != NULL; dsp = dsp->next) {
      switch (dsp->choice) {
        case 1:
          slocp = (SeqLocPtr)(dsp->data.ptrvalue);
          if (slocp == NULL) break;
          if (slocp->choice == SEQLOC_NULL) {
            num_gaps++;
          }
          break;
        case 2:
          slitp = (SeqLitPtr)(dsp->data.ptrvalue);
          if (slitp == NULL) break;
          if (slitp->seq_data == NULL || slitp->seq_data_type == Seq_code_gap) {
            num_gaps++;
          }
          break;
        default:
          break;
      }
    }
  }

  ValNodeCopyStr (&head, 0, "* NOTE: This is a partial genome representation.");
  if (num_gaps > 0) {
    sprintf (buffer, " It currently~* consists of %ld contigs. The true order of the pieces~", (long) (num_gaps + 1));
    ValNodeCopyStr (&head, 0, buffer);
    ValNodeCopyStr (&head, 0, "* is not known and their order in this sequence record is~");
    ValNodeCopyStr (&head, 0, "* arbitrary. Gaps between the contigs are represented as~");
    ValNodeCopyStr (&head, 0, "* runs of N, but the exact sizes of the gaps are unknown.");
  }
  ValNodeCopyStr (&head, 0, "~");

  str = MergeFFValNodeStrs (head);

  FFAddOneString (ffstring, str, TRUE, TRUE, TILDE_EXPAND);

  MemFree (str);
  ValNodeFreeData (head);
}

static void AddHTGSCommentString (
  StringItemPtr ffstring,
  BioseqPtr bsp,
  MolInfoPtr mip
)

{
  CharPtr      buf = NULL;
  Char         buffer [256];
  Int4         buflen = 0;
  DeltaSeqPtr  dsp;
  ValNodePtr   head = NULL;
  Int4         num_s = 0;
  Int4         num_g = 0;
  CharPtr      str = NULL;

  if (bsp == NULL || mip == NULL || mip->tech < 2) return;

  if (bsp->repr == Seq_repr_delta) {
    for (dsp = (DeltaSeqPtr) bsp->seq_ext, buflen = 0; dsp != NULL; dsp = dsp->next) {
      buflen += 80;
    }
    if (buflen > 0) {
      buf = MemNew ((size_t) (buflen + 1));
      if (buf == NULL) return;
      CountGapsInDeltaSeq (bsp, &num_s, &num_g, NULL, NULL, buf, buflen);
    }
  }

  if (mip->tech == MI_TECH_htgs_0) {

    if (num_s > 0) {
      sprintf (buffer, "* NOTE: This record contains %ld individual~", (long) (num_g + 1));
      ValNodeCopyStr (&head, 0, buffer);
      ValNodeCopyStr (&head, 0, "* sequencing reads that have not been assembled into~");
      ValNodeCopyStr (&head, 0, "* contigs. Runs of N are used to separate the reads~");
      ValNodeCopyStr (&head, 0, "* and the order in which they appear is completely~");
      ValNodeCopyStr (&head, 0, "* arbitrary. Low-pass sequence sampling is useful for~");
      ValNodeCopyStr (&head, 0, "* identifying clones that may be gene-rich and allows~");
      ValNodeCopyStr (&head, 0, "* overlap relationships among clones to be deduced.~");
      ValNodeCopyStr (&head, 0, "* However, it should not be assumed that this clone~");
      ValNodeCopyStr (&head, 0, "* will be sequenced to completion. In the event that~");
      ValNodeCopyStr (&head, 0, "* the record is updated, the accession number will~");
      ValNodeCopyStr (&head, 0, "* be preserved.");
    }
    ValNodeCopyStr (&head, 0, "~");
    ValNodeCopyStr (&head, 0, buf);

  } else if (mip->tech == MI_TECH_htgs_1) {

    ValNodeCopyStr (&head, 0, "* NOTE: This is a \"working draft\" sequence.");
    if (num_s > 0) {
      sprintf (buffer, " It currently~* consists of %ld contigs. The true order of the pieces~", (long) (num_g + 1));
      ValNodeCopyStr (&head, 0, buffer);
      ValNodeCopyStr (&head, 0, "* is not known and their order in this sequence record is~");
      ValNodeCopyStr (&head, 0, "* arbitrary. Gaps between the contigs are represented as~");
      ValNodeCopyStr (&head, 0, "* runs of N, but the exact sizes of the gaps are unknown.");
    }
    ValNodeCopyStr (&head, 0, "~* This record will be updated with the finished sequence~");
    ValNodeCopyStr (&head, 0, "* as soon as it is available and the accession number will~");
    ValNodeCopyStr (&head, 0, "* be preserved.");
    ValNodeCopyStr (&head, 0, "~");
    ValNodeCopyStr (&head, 0, buf);

  } else if (mip->tech == MI_TECH_htgs_2) {

    ValNodeCopyStr (&head, 0, "* NOTE: This is a \"working draft\" sequence.");
    if (num_s > 0) {
      sprintf (buffer, " It currently~* consists of %ld contigs. Gaps between the contigs~", (long) (num_g + 1));
      ValNodeCopyStr (&head, 0, buffer);
      ValNodeCopyStr (&head, 0, "* are represented as runs of N. The order of the pieces~");
      ValNodeCopyStr (&head, 0, "* is believed to be correct as given, however the sizes~");
      ValNodeCopyStr (&head, 0, "* of the gaps between them are based on estimates that have~");
      ValNodeCopyStr (&head, 0, "* provided by the submitter.");
    }
    ValNodeCopyStr (&head, 0, "~* This sequence will be replaced~");
    ValNodeCopyStr (&head, 0, "* by the finished sequence as soon as it is available and~");
    ValNodeCopyStr (&head, 0, "* the accession number will be preserved.");
    ValNodeCopyStr (&head, 0, "~");
    ValNodeCopyStr (&head, 0, buf);

  } else if ((str = StringForSeqTech (mip->tech)) != NULL) {

      sprintf (buffer, "Method: %s.", str);
      ValNodeCopyStr (&head, 0, buffer);
  }

  MemFree (buf);

  str = MergeFFValNodeStrs (head);

  FFAddOneString (ffstring, str, TRUE, TRUE, TILDE_EXPAND);

  MemFree (str);
  ValNodeFreeData (head);
}

static void AddWGSMasterCommentString (
  StringItemPtr ffstring,
  BioseqPtr bsp,
  CharPtr wgsaccn,
  CharPtr wgsname
)

{
  size_t             acclen;
  BioSourcePtr       biop;
  Char               buf [256];
  SeqMgrDescContext  dcontext;
  CharPtr            first = NULL;
  CharPtr            last = NULL;
  ObjectIdPtr        oip;
  OrgRefPtr          orp;
  SeqDescrPtr        sdp;
  CharPtr            taxname = NULL;
  UserFieldPtr       ufp;
  UserObjectPtr      uop;
  Char               ver [16];

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp != NULL) {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
    if (biop != NULL) {
      orp = biop->org;
      if (orp != NULL) {
        taxname = orp->taxname;
      }
    }
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &dcontext);
  while (sdp != NULL) {
    uop = (UserObjectPtr) sdp->data.ptrvalue;
    if (uop != NULL) {
      oip = uop->type;
      if (oip != NULL && StringICmp (oip->str, "WGSProjects") == 0) {
        for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
          oip = ufp->label;
          if (oip == NULL || oip->str == NULL || ufp->choice != 1) continue;
          if (StringICmp (oip->str, "WGS_accession_first") == 0) {
            first = (CharPtr) ufp->data.ptrvalue;
          } else if (StringICmp (oip->str, "WGS_accession_last") == 0) {
            last = (CharPtr) ufp->data.ptrvalue;
          }
        }
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &dcontext);
  }

  if (StringHasNoText (taxname)) {
    taxname = "?";
  }
  ver [0] = '\0';
  acclen = StringLen (wgsname);
  if (acclen == 12) {
    StringCpy (ver, wgsname + 4);
    ver [2] = '\0';
  } else if (acclen == 13) {
    StringCpy (ver, wgsname + 4);
    ver [2] = '\0';
  } else if (acclen == 14) {
    StringCpy (ver, wgsname + 4);
    ver [2] = '\0';
  } else if (acclen == 15) {
    StringCpy (ver, wgsname + 7);
    ver [2] = '\0';
  } else if (acclen == 16) {
    StringCpy (ver, wgsname + 7);
    ver [2] = '\0';
  }

  sprintf (buf, "The %s whole genome shotgun (WGS) project has the project accession %s.", taxname, wgsaccn);
  FFAddOneString(ffstring, buf, TRUE, FALSE, TILDE_EXPAND);

  sprintf (buf, "  This version of the project (%s) has the accession number %s", ver, wgsname);
  FFAddOneString(ffstring, buf, FALSE, FALSE, TILDE_EXPAND);

  if (first == NULL && last == NULL) {
    sprintf (buf, ".");
    FFAddOneString(ffstring, buf, TRUE, FALSE, TILDE_EXPAND);
  } else {
    if (first != NULL && last == NULL) {
      last = first;
    } else if (first == NULL && last != NULL) {
      first = last;
    }
    if (StringDoesHaveText (first) && StringDoesHaveText (last)) {
      if (StringCmp (first, last) != 0) {
        sprintf (buf, ", and consists of sequences %s-%s.", first, last);
        FFAddOneString(ffstring, buf, TRUE, FALSE, TILDE_EXPAND);
      } else {
        sprintf (buf, ", and consists of sequence %s.", first);
        FFAddOneString(ffstring, buf, TRUE, FALSE, TILDE_EXPAND);
      }
    } else {
      sprintf (buf, ".");
      FFAddOneString(ffstring, buf, TRUE, FALSE, TILDE_EXPAND);
    }
  }
}

static void AddTSAMasterCommentString (
  StringItemPtr ffstring,
  BioseqPtr bsp,
  CharPtr tsaaccn,
  CharPtr tsaname
)

{
  size_t             acclen;
  BioSourcePtr       biop;
  Char               buf [256];
  SeqMgrDescContext  dcontext;
  CharPtr            first = NULL;
  CharPtr            last = NULL;
  ObjectIdPtr        oip;
  OrgRefPtr          orp;
  SeqDescrPtr        sdp;
  CharPtr            taxname = NULL;
  UserFieldPtr       ufp;
  UserObjectPtr      uop;
  Char               ver [16];

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp != NULL) {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
    if (biop != NULL) {
      orp = biop->org;
      if (orp != NULL) {
        taxname = orp->taxname;
      }
    }
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &dcontext);
  while (sdp != NULL) {
    uop = (UserObjectPtr) sdp->data.ptrvalue;
    if (uop != NULL) {
      oip = uop->type;
      if (oip != NULL) {
        if (StringICmp (oip->str, "TSA-mRNA-List") == 0 || StringICmp (oip->str, "TSA-RNA-List") == 0) {
          for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
            oip = ufp->label;
            if (oip == NULL || oip->str == NULL || ufp->choice != 1) continue;
            if (StringICmp (oip->str, "TSA_accession_first") == 0) {
              first = (CharPtr) ufp->data.ptrvalue;
            } else if (StringICmp (oip->str, "TSA_accession_last") == 0) {
              last = (CharPtr) ufp->data.ptrvalue;
            }
          }
        }
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &dcontext);
  }

  if (StringHasNoText (taxname)) {
    taxname = "?";
  }
  ver [0] = '\0';
  acclen = StringLen (tsaname);
  if (acclen == 12) {
    StringCpy (ver, tsaname + 4);
    ver [2] = '\0';
  } else if (acclen == 13) {
    StringCpy (ver, tsaname + 4);
    ver [2] = '\0';
  } else if (acclen == 14) {
    StringCpy (ver, tsaname + 4);
    ver [2] = '\0';
  } else if (acclen == 15) {
    StringCpy (ver, tsaname + 7);
    ver [2] = '\0';
  }

  sprintf (buf, "The %s transcriptome shotgun assembly (TSA) project has the project accession %s.", taxname, tsaaccn);
  FFAddOneString(ffstring, buf, TRUE, FALSE, TILDE_EXPAND);

  sprintf (buf, "  This version of the project (%s) has the accession number %s", ver, tsaname);
  FFAddOneString(ffstring, buf, FALSE, FALSE, TILDE_EXPAND);

  if (first == NULL && last == NULL) {
    sprintf (buf, ".");
    FFAddOneString(ffstring, buf, TRUE, FALSE, TILDE_EXPAND);
  } else {
    if (first != NULL && last == NULL) {
      last = first;
    } else if (first == NULL && last != NULL) {
      first = last;
    }
    if (StringDoesHaveText (first) && StringDoesHaveText (last)) {
      if (StringCmp (first, last) != 0) {
        sprintf (buf, ", and consists of sequences %s-%s.", first, last);
        FFAddOneString(ffstring, buf, TRUE, FALSE, TILDE_EXPAND);
      } else {
        sprintf (buf, ", and consists of sequence %s.", first);
        FFAddOneString(ffstring, buf, TRUE, FALSE, TILDE_EXPAND);
      }
    } else {
      sprintf (buf, ".");
      FFAddOneString(ffstring, buf, TRUE, FALSE, TILDE_EXPAND);
    }
  }
}

static void AddTLSMasterCommentString (
  StringItemPtr ffstring,
  BioseqPtr bsp,
  CharPtr tlsaccn,
  CharPtr tlsname
)

{
  size_t             acclen;
  BioSourcePtr       biop;
  Char               buf [256];
  SeqMgrDescContext  dcontext;
  CharPtr            first = NULL;
  CharPtr            last = NULL;
  ObjectIdPtr        oip;
  OrgRefPtr          orp;
  SeqDescrPtr        sdp;
  CharPtr            taxname = NULL;
  UserFieldPtr       ufp;
  UserObjectPtr      uop;
  Char               ver [16];

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp != NULL) {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
    if (biop != NULL) {
      orp = biop->org;
      if (orp != NULL) {
        taxname = orp->taxname;
      }
    }
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &dcontext);
  while (sdp != NULL) {
    uop = (UserObjectPtr) sdp->data.ptrvalue;
    if (uop != NULL) {
      oip = uop->type;
      if (oip != NULL && StringICmp (oip->str, "TLSProjects") == 0) {
        for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
          oip = ufp->label;
          if (oip == NULL || oip->str == NULL || ufp->choice != 1) continue;
          if (StringICmp (oip->str, "TLS_accession_first") == 0) {
            first = (CharPtr) ufp->data.ptrvalue;
          } else if (StringICmp (oip->str, "TLS_accession_last") == 0) {
            last = (CharPtr) ufp->data.ptrvalue;
          }
        }
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &dcontext);
  }

  if (StringHasNoText (taxname)) {
    taxname = "?";
  }
  ver [0] = '\0';
  acclen = StringLen (tlsname);
  if (acclen == 12) {
    StringCpy (ver, tlsname + 4);
    ver [2] = '\0';
  } else if (acclen == 13) {
    StringCpy (ver, tlsname + 4);
    ver [2] = '\0';
  } else if (acclen == 14) {
    StringCpy (ver, tlsname + 4);
    ver [2] = '\0';
  } else if (acclen == 15) {
    StringCpy (ver, tlsname + 7);
    ver [2] = '\0';
  } else if (acclen == 16) {
    StringCpy (ver, tlsname + 7);
    ver [2] = '\0';
  }

  sprintf (buf, "The %s targeted locus study (TLS) project has the project accession %s.", taxname, tlsaccn);
  FFAddOneString(ffstring, buf, TRUE, FALSE, TILDE_EXPAND);

  sprintf (buf, "  This version of the project (%s) has the accession number %s", ver, tlsname);
  FFAddOneString(ffstring, buf, FALSE, FALSE, TILDE_EXPAND);

  if (first == NULL && last == NULL) {
    sprintf (buf, ".");
    FFAddOneString(ffstring, buf, TRUE, FALSE, TILDE_EXPAND);
  } else {
    if (first != NULL && last == NULL) {
      last = first;
    } else if (first == NULL && last != NULL) {
      first = last;
    }
    if (StringDoesHaveText (first) && StringDoesHaveText (last)) {
      if (StringCmp (first, last) != 0) {
        sprintf (buf, ", and consists of sequences %s-%s.", first, last);
        FFAddOneString(ffstring, buf, TRUE, FALSE, TILDE_EXPAND);
      } else {
        sprintf (buf, ", and consists of sequence %s.", first);
        FFAddOneString(ffstring, buf, TRUE, FALSE, TILDE_EXPAND);
      }
    } else {
      sprintf (buf, ".");
      FFAddOneString(ffstring, buf, TRUE, FALSE, TILDE_EXPAND);
    }
  }
}

static CharPtr GetMolInfoCommentString (
  BioseqPtr bsp,
  MolInfoPtr mip
)

{
  Boolean  is_aa;
  CharPtr  str = NULL;

  if (bsp == NULL || mip == NULL) return NULL;

  is_aa = ISA_aa (bsp->mol);
  switch (mip->completeness) {
    case 1 :
      str = "COMPLETENESS: full length";
      break;
    case 2 :
      str = "COMPLETENESS: not full length";
      break;
    case 3 :
      if (is_aa) {
        str = "COMPLETENESS: incomplete on the amino end";
      } else {
        str = "COMPLETENESS: incomplete on the 5' end";
      }
      break;
    case 4 :
      if (is_aa) {
        str = "COMPLETENESS: incomplete on the carboxy end";
      } else {
        str = "COMPLETENESS: incomplete on the 3' end";
      }
      break;
    case 5 :
      str = "COMPLETENESS: incomplete on both ends";
      break;
    case 6 :
      if (is_aa) {
        str = "COMPLETENESS: complete on the amino end";
      } else {
        str = "COMPLETENESS: complete on the 5' end";
      }
      break;
    case 7 :
      if (is_aa) {
        str = "COMPLETENESS: complete on the carboxy end";
      } else {
        str = "COMPLETENESS: complete on the 3' end";
      }
      break;
    default :
      str = "COMPLETENESS: unknown";
      break;
  }

  return str;
}

static CharPtr GetStrForBankit (
  UserObjectPtr uop,
  Boolean dumpMode,
  Boolean showedLocalId
)

{
  CharPtr       bic = NULL, smc = NULL, uvc = NULL, pfx = NULL, ptr;
  ObjectIdPtr   oip;
  UserFieldPtr  ufp;

  if (uop == NULL) return NULL;
  if ((oip = uop->type) == NULL) return NULL;
  if (StringCmp (oip->str, "Submission") != 0) return NULL;

  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    oip = ufp->label;
    if (StringCmp(oip->str, "UniVecComment") == 0) {
      uvc = ufp->data.ptrvalue;
    } else if (StringCmp(oip->str, "AdditionalComment") == 0) {
      bic = ufp->data.ptrvalue;
    } else if (StringCmp(oip->str, "SmartComment") == 0 && dumpMode) {
      smc = ufp->data.ptrvalue;
    }
  }

  if (showedLocalId) {
    if (StringNICmp (bic, "LocalID:", 8) == 0) {
      bic = NULL;
    }
    if (StringNICmp (smc, "LocalID:", 8) == 0) {
      smc = NULL;
    }
  }

  if (uvc == NULL && bic == NULL && smc == NULL) return NULL;

  ptr = (CharPtr) MemNew (StringLen (uvc) + StringLen (bic) + StringLen (smc) + 45);
  if (uvc != NULL) {
    StringCat (ptr, pfx);
    StringCat (ptr, "Vector Explanation: ");
    StringCat (ptr, uvc);
    pfx = "~";
  }
  if (bic != NULL) {
    StringCat (ptr, pfx);
    StringCat (ptr, "Bankit Comment: ");
    StringCat (ptr, bic);
    pfx = "~";
  }
  if (smc != NULL) {
    StringCat (ptr, pfx);
    StringCat (ptr, "Bankit Comment: ");
    StringCat (ptr, smc);
    pfx = "~";
  }

  return ptr;
}

static CharPtr reftxt0 = " The reference sequence was derived from ";
static CharPtr reftxtg = " The reference sequence was generated based on analysis of ";
static CharPtr reftxti = " The reference sequence is identical to ";
static CharPtr reftxt1 = " This record is predicted by genome sequence analysis and is not yet supported by experimental evidence.";
static CharPtr reftxt2 = " This record has not yet been subject to final NCBI review.";
static CharPtr reftxt3 = " This record has not been reviewed and the function is unknown.";
static CharPtr reftxt4 = " This record has undergone validation or preliminary review.";
static CharPtr reftxt5 = " This record has been curated by ";
static CharPtr reftxt6 = " This record is predicted by automated computational analysis.";
static CharPtr reftxt7 = " This record is provided to represent a collection of whole genome shotgun sequences.";
static CharPtr reftxt9 = " This record is derived from an annotated genomic sequence (";
static CharPtr reftxt21 = " NCBI contigs are derived from assembled genomic sequence data.";
static CharPtr reftxt22 = " Features on this sequence have been produced for build ";
static CharPtr reftxt23 = " of the NCBI's genome annotation";
static CharPtr reftxt41 = " This record is based on preliminary annotation provided by ";
static CharPtr reftxt51 = " This record represents a single, non-redundant, protein sequence which may be annotated on many different RefSeq genomes from the same, or different, species";

static CharPtr GetStatusForRefTrack (
  UserObjectPtr uop
)

{
  CharPtr       st;
  ObjectIdPtr   oip;
  UserFieldPtr  ufp, urf = NULL;

  if (uop == NULL) return NULL;
  if ((oip = uop->type) == NULL) return NULL;
  if (StringCmp (oip->str, "RefGeneTracking") != 0) return NULL;
  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    oip = ufp->label;
    if (StringCmp(oip->str, "Assembly") == 0) {
      urf = ufp;
    }
  }
  /* if (urf == NULL || urf->choice != 11) return NULL; */
  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    oip = ufp->label;
    if (StringCmp (oip->str, "Status") == 0) {
      st = (CharPtr) ufp->data.ptrvalue;
      if (StringICmp (st, "Inferred") == 0) {
        return "INFERRED ";
      } else if (StringICmp (st, "Provisional") == 0) {
        return "PROVISIONAL ";
      } else if (StringICmp (st, "Predicted") == 0) {
        return "PREDICTED ";
      } else if (StringICmp (st, "Validated") == 0) {
        return "VALIDATED ";
      } else if (StringICmp (st, "Reviewed") == 0) {
        return "REVIEWED ";
      } else if (StringICmp (st, "Model") == 0) {
        return "MODEL ";
      } else if (StringICmp (st, "WGS") == 0) {
        return "WGS ";
      } else if (StringICmp (st, "Pipeline") == 0) {
        return "Pipeline ";
      }
    }
  }
  return NULL;
}


static Boolean URLHasSuspiciousHtml (
  IntAsn2gbJobPtr ajp,
  CharPtr searchString
)

{
  Char        ch;
  CharPtr     ptr;
  Int4        state;
  ValNodePtr  matches;

  if (StringHasNoText (searchString)) return FALSE;

  state = 0;
  ptr = searchString;
  ch = *ptr;

  while (ch != '\0') {
    matches = NULL;
    ch = TO_LOWER (ch);
    state = TextFsaNext (ajp->bad_html_fsa, state, ch, &matches);
    if (matches != NULL) {
      return TRUE;
    }
    ptr++;
    ch = *ptr;
  }

  return FALSE;
}

static Boolean GetGiFromAccnDotVer (CharPtr source, BIG_ID_PNTR gip)

{
  BIG_ID    gi = 0;
  SeqIdPtr  sip;

  if (StringHasNoText (source) || gip == NULL) return FALSE;
  *gip = 0;

  sip = SeqIdFromAccessionDotVersion (source);
  if (sip == NULL) return FALSE;
  gi = GetGIForSeqId (sip);
  SeqIdFree (sip);
  if (gi == 0) return FALSE;

  *gip = gi;
  return TRUE;
}

static void AddStrForRefTrack (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  UserObjectPtr uop,
  Boolean is_na,
  CharPtr genomeBuildNumber,
  CharPtr genomeVersionNumber
)

{
  CharPtr       accn, curator = NULL, name, source = NULL, st, url = NULL;
  Char          buf [64];
  ObjectIdPtr   oip;
  UserFieldPtr  ufp, tmp, u, urf = NULL;
  Int4          from, to;
  BIG_ID        gi;
  Int2          i = 0;
  Int2          review = 0;
  Boolean       generated = FALSE, identical = FALSE;

  if ( uop == NULL || ffstring == NULL ) return;
  if ((oip = uop->type) == NULL) return;
  if (StringCmp (oip->str, "RefGeneTracking") != 0) return;

  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    oip = ufp->label;
    if (StringCmp(oip->str, "Assembly") == 0) {
      urf = ufp;
    } else if (StringCmp(oip->str, "IdenticalTo") == 0) {
      urf = ufp;
      identical = TRUE;
    }
    if (StringCmp (oip->str, "Status") == 0) {
      st = (CharPtr) ufp->data.ptrvalue;
      if (StringICmp (st, "Inferred") == 0) {
        review = 1;
      } else if (StringICmp (st, "Provisional") == 0) {
        review = 2;
      } else if (StringICmp (st, "Predicted") == 0) {
        review = 3;
      } else if (StringICmp (st, "Validated") == 0) {
        review = 4;
      } else if (StringICmp (st, "Reviewed") == 0) {
        review = 5;
      } else if (StringICmp (st, "Model") == 0) {
        review = 6;
      } else if (StringICmp (st, "WGS") == 0) {
        review = 7;
      } else if (StringICmp (st, "Pipeline") == 0) {
        review = 8;
      }
    } else if (StringCmp (oip->str, "Generated") == 0) {
      generated = ufp->data.boolvalue;
    } else if (StringCmp (oip->str, "Collaborator") == 0) {
      st = (CharPtr) ufp->data.ptrvalue;
      if (! StringHasNoText (st)) {
        curator = st;
      }
    } else if (StringCmp (oip->str, "CollaboratorURL") == 0) {
      st = (CharPtr) ufp->data.ptrvalue;
      if (! StringHasNoText (st)) {
        url = st;
      }
    } else if (StringCmp (oip->str, "GenomicSource") == 0) {
      st = (CharPtr) ufp->data.ptrvalue;
      if (! StringHasNoText (st)) {
        source = st;
      }
    }
  }
  if (urf != NULL && urf->choice == 11) {
    for (tmp = urf->data.ptrvalue; tmp != NULL; tmp = tmp->next) {
      if (tmp->choice != 11) continue;
      for (u = tmp->data.ptrvalue; u != NULL; u = u->next) {
        oip = u->label;
        if (oip == NULL) continue;
        if (StringCmp (oip->str, "accession") == 0 ||
            StringCmp (oip->str, "name") == 0) {
          i++;
        }
      }
    }
  }
  if ( GetWWW(ajp) ) {
    FFAddOneString (ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
    FF_Add_NCBI_Base_URL (ffstring, ref_link);
    FFAddOneString (ffstring, "\">", FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString (ffstring, "REFSEQ", FALSE, FALSE, TILDE_IGNORE);
    if (review == 8) {
      FFAddOneString (ffstring, " INFORMATION", FALSE, FALSE, TILDE_IGNORE);
    }
    FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
  } else {
    FFAddOneString (ffstring, "REFSEQ", FALSE, FALSE, TILDE_IGNORE);
    if (review == 8) {
      FFAddOneString (ffstring, " INFORMATION", FALSE, FALSE, TILDE_IGNORE);
    }
  }
  FFAddOneString (ffstring, ":", FALSE, FALSE, TILDE_IGNORE);
  if (review == 1) {
    FFAddOneString (ffstring, reftxt1, FALSE, FALSE, TILDE_IGNORE);
  } else if (review == 2) {
    if (curator == NULL) {
      FFAddOneString (ffstring, reftxt2, FALSE, FALSE, TILDE_IGNORE);
    }
  } else if (review == 3) {
    FFAddOneString (ffstring, reftxt3, FALSE, FALSE, TILDE_IGNORE);
  } else if (review == 4) {
    FFAddOneString (ffstring, reftxt4, FALSE, FALSE, TILDE_IGNORE);
  } else if (review == 5) {
    if (curator == NULL) {
      curator = "NCBI staff";
    }
  } else if (review == 6) {
    FFAddOneString (ffstring, reftxt6, FALSE, FALSE, TILDE_IGNORE);
  } else if (review == 7) {
    FFAddOneString (ffstring, reftxt7, FALSE, FALSE, TILDE_IGNORE);
  } else if (review == 8) {
  }
  if (curator != NULL) {
    if (review == 2) {
      FFAddOneString (ffstring, reftxt41, FALSE, FALSE, TILDE_IGNORE);
    } else {
      FFAddOneString (ffstring, reftxt5, FALSE, FALSE, TILDE_IGNORE);
    }
    if (GetWWW (ajp) && url != NULL && (! URLHasSuspiciousHtml (ajp, url))) {
      if (StringNCmp (url, "http://", 7) == 0 || StringNCmp (url, "https://", 8) == 0) {
        FFAddTextToString(ffstring, "<a href=\"", url, "\">", FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, curator, FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
      } else if (StringNCmp (url, "www.", 4) == 0) {
        FFAddTextToString(ffstring, "<a href=http://\"", url, "\">", FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, curator, FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
      } else {
        FFAddOneString (ffstring, curator, FALSE, FALSE, TILDE_IGNORE);
      }
    } else {
      FFAddOneString (ffstring, curator, FALSE, FALSE, TILDE_IGNORE);
    }
    FFAddOneString (ffstring, ".", FALSE, FALSE, TILDE_IGNORE);
  }
  if (source != NULL) {
    FFAddOneString (ffstring, reftxt9, FALSE, FALSE, TILDE_IGNORE);
    gi = 0;
    if (GetWWW (ajp) && ValidateAccnDotVer (source) == 0 && GetGiFromAccnDotVer (source, &gi)) {
      FFAddOneString (ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
      if (is_na) {
        FF_Add_NCBI_Base_URL (ffstring, link_seqn);
      } else {
        FF_Add_NCBI_Base_URL (ffstring, link_seqp);
      }
      if (gi > 0) {
        sprintf (buf, "%ld", (long) gi);
        FFAddTextToString(ffstring, /* "val=" */ NULL, buf, "\">", FALSE, FALSE, TILDE_IGNORE);
      } else {
        FFAddTextToString(ffstring, /* "val=" */ NULL, source, "\">", FALSE, FALSE, TILDE_IGNORE);
      }
      FFAddOneString (ffstring, source, FALSE, FALSE, TILDE_IGNORE);
      FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
    } else {
      FFAddOneString (ffstring, source, FALSE, FALSE, TILDE_IGNORE);
    }
    FFAddOneString (ffstring, ").", FALSE, FALSE, TILDE_IGNORE);
  }
  if (i > 0) {
    if (review == 8 && (genomeBuildNumber != NULL || genomeVersionNumber != NULL)) {
      FFAddOneString (ffstring, reftxt22, FALSE, FALSE, TILDE_EXPAND);
      FFAddOneString (ffstring, genomeBuildNumber, FALSE, FALSE, TILDE_EXPAND);
      if (StringHasNoText (genomeVersionNumber)) {
        genomeVersionNumber = "1";
      }
      FFAddOneString (ffstring, " version ", FALSE, FALSE, TILDE_EXPAND);
      FFAddOneString (ffstring, genomeVersionNumber, FALSE, FALSE, TILDE_EXPAND);
      FFAddOneString (ffstring, reftxt23, FALSE, FALSE, TILDE_EXPAND);

      FFAddOneString (ffstring, " [see ", FALSE, FALSE, TILDE_EXPAND);

      if ( GetWWW(ajp) ) {
        FFAddOneString (ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
        FF_Add_NCBI_Base_URL (ffstring, doc_link);
        FFAddOneString (ffstring, "\">", FALSE, FALSE, TILDE_IGNORE);
      }
      FFAddOneString (ffstring, "documentation", FALSE, FALSE, TILDE_IGNORE);
      if ( GetWWW(ajp) ) {
        FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
      }

      FFAddOneString (ffstring, "].", FALSE, FALSE, TILDE_EXPAND);
    }
    if (generated) {
      FFAddOneString (ffstring, reftxtg, FALSE, FALSE, TILDE_IGNORE);
    } else if (identical) {
      FFAddOneString (ffstring, reftxti, FALSE, FALSE, TILDE_IGNORE);
    } else {
      FFAddOneString (ffstring, reftxt0, FALSE, FALSE, TILDE_IGNORE);
    }

    for (tmp = urf->data.ptrvalue; tmp != NULL; tmp = tmp->next) {
      accn = NULL;
      from = 0;
      to = 0;
      name = NULL;
      gi = 0;
      for (u = tmp->data.ptrvalue; u != NULL; u = u->next) {
        oip = u->label;
        if (oip != NULL && oip->str != NULL) {
          if (StringICmp (oip->str, "accession") == 0 && u->choice == 1) {
            accn = (CharPtr) u->data.ptrvalue;
          } else if (StringICmp (oip->str, "from") == 0 && u->choice == 2) {
            from = u->data.intvalue;
          } else if (StringICmp (oip->str, "to") == 0 && u->choice == 2) {
            to = u->data.intvalue;
          } else if (StringICmp (oip->str, "name") == 0 && u->choice == 1) {
            name = (CharPtr) u->data.ptrvalue;
          } else if (StringICmp (oip->str, "gi") == 0 && u->choice == 2) {
            gi = (BIG_ID) u->data.intvalue;
          }
        }
      }
      if (StringDoesHaveText (accn)) {
        if (GetWWW (ajp) && ValidateAccnDotVer (accn) == 0 && GetGiFromAccnDotVer (accn, &gi)) {
          FFAddOneString (ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
          if (is_na) {
            FF_Add_NCBI_Base_URL (ffstring, link_seqn);
          } else {
            FF_Add_NCBI_Base_URL (ffstring, link_seqp);
          }
          if (gi > 0) {
            sprintf (buf, "%ld", (long) gi);
            FFAddTextToString(ffstring, /* "val=" */ NULL, buf, "\">", FALSE, FALSE, TILDE_IGNORE);
          } else {
            FFAddTextToString(ffstring, /* "val=" */ NULL, accn, "\">", FALSE, FALSE, TILDE_IGNORE);
          }
          FFAddOneString (ffstring, accn, FALSE, FALSE, TILDE_IGNORE);
          FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
        } else if (GetWWW (ajp) && ValidateAccn (accn) == 0) {
          FFAddOneString (ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
          if (is_na) {
            FF_Add_NCBI_Base_URL (ffstring, link_seqn);
          } else {
            FF_Add_NCBI_Base_URL (ffstring, link_seqp);
          }
          FFAddTextToString(ffstring, /* "val=" */ NULL, accn, "\">", FALSE, FALSE, TILDE_IGNORE);
          FFAddOneString (ffstring, accn, FALSE, FALSE, TILDE_IGNORE);
          FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
        } else {
          FFAddOneString (ffstring, accn, FALSE, FALSE, TILDE_IGNORE);
        }
        if (from > 0 && to > 0) {
          sprintf (buf, " (range: %ld-%ld)", (long) from, (long) to);
          FFAddOneString (ffstring, buf, FALSE, FALSE, TILDE_IGNORE);
        }
      } else if (StringDoesHaveText (name)) {
        FFAddOneString (ffstring, name, FALSE, FALSE, TILDE_IGNORE);
      } else continue;
      if (tmp->next != NULL) {
        ufp = tmp->next;
        if (ufp->next != NULL) {
          FFAddOneString (ffstring, ", ", FALSE, FALSE, TILDE_IGNORE);
        } else {
          FFAddOneString (ffstring, " and ", FALSE, FALSE, TILDE_IGNORE);
        }
      }
    }
    FFAddOneString (ffstring, ".", FALSE, FALSE, TILDE_EXPAND);
  }
}

static void AddStrForRefSeqGenome (
  IntAsn2gbJobPtr ajp,
  StringItemPtr ffstring,
  UserObjectPtr uop
)

{
  CharPtr       category = NULL, calc = NULL, cca = NULL, cli = NULL, com = NULL,
                fgs = NULL, mod = NULL, phy = NULL, prt = NULL, qfo = NULL,
                tys = NULL, upr = NULL;
  ObjectIdPtr   oip;
  UserFieldPtr  ufp, tmp, urf = NULL;

  if ( uop == NULL || ffstring == NULL ) return;
  if ((oip = uop->type) == NULL) return;
  if (StringCmp (oip->str, "RefSeqGenome") != 0) return;

  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    oip = ufp->label;
    if (oip == NULL) continue;
    if (StringCmp (oip->str, "RefSeq Category") == 0) {
      category = (CharPtr) ufp->data.ptrvalue;
    } else if (StringCmp (oip->str, "Details") == 0) {
      urf = ufp;
    }
  }
  if (urf != NULL && urf->choice == 11) {
    for (tmp = urf->data.ptrvalue; tmp != NULL; tmp = tmp->next) {
      oip = tmp->label;
      if (StringCmp (oip->str, "CALC") == 0) {
        calc = (CharPtr) tmp->data.ptrvalue;
      } else if (StringCmp (oip->str, "CCA") == 0) {
        cca = (CharPtr) tmp->data.ptrvalue;
      } else if (StringCmp (oip->str, "CLI") == 0) {
        cli = (CharPtr) tmp->data.ptrvalue;
      } else if (StringCmp (oip->str, "COM") == 0) {
        com = (CharPtr) tmp->data.ptrvalue;
      } else if (StringCmp (oip->str, "FGS") == 0) {
        fgs = (CharPtr) tmp->data.ptrvalue;
      } else if (StringCmp (oip->str, "MOD") == 0) {
        mod = (CharPtr) tmp->data.ptrvalue;
      } else if (StringCmp (oip->str, "PHY") == 0) {
        phy = (CharPtr) tmp->data.ptrvalue;
      } else if (StringCmp (oip->str, "PRT") == 0) {
        prt = (CharPtr) tmp->data.ptrvalue;
      } else if (StringCmp (oip->str, "QfO") == 0) {
        qfo = (CharPtr) tmp->data.ptrvalue;
      } else if (StringCmp (oip->str, "TYS") == 0) {
        tys = (CharPtr) tmp->data.ptrvalue;
      } else if (StringCmp (oip->str, "UPR") == 0) {
        upr = (CharPtr) tmp->data.ptrvalue;
      }
    }
  }
  FFAddOneString (ffstring, "RefSeq Category: ", FALSE, FALSE, TILDE_IGNORE);
  FFAddOneString (ffstring, category, FALSE, FALSE, TILDE_IGNORE);
  if (calc != NULL) {
    FFAddOneString (ffstring, "\n           CALC: ", FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString (ffstring, calc, FALSE, FALSE, TILDE_IGNORE);
  }
  if (cca != NULL) {
    FFAddOneString (ffstring, "\n            CCA: ", FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString (ffstring, cca, FALSE, FALSE, TILDE_IGNORE);
  }
  if (cli != NULL) {
    FFAddOneString (ffstring, "\n            CLI: ", FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString (ffstring, cli, FALSE, FALSE, TILDE_IGNORE);
  }
  if (com != NULL) {
    FFAddOneString (ffstring, "\n            COM: ", FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString (ffstring, com, FALSE, FALSE, TILDE_IGNORE);
  }
  if (fgs != NULL) {
    FFAddOneString (ffstring, "\n            FGS: ", FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString (ffstring, fgs, FALSE, FALSE, TILDE_IGNORE);
  }
  if (mod != NULL) {
    FFAddOneString (ffstring, "\n            MOD: ", FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString (ffstring, mod, FALSE, FALSE, TILDE_IGNORE);
  }
  if (phy != NULL) {
    FFAddOneString (ffstring, "\n            PHY: ", FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString (ffstring, phy, FALSE, FALSE, TILDE_IGNORE);
  }
  if (prt != NULL) {
    FFAddOneString (ffstring, "\n            PRT: ", FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString (ffstring, prt, FALSE, FALSE, TILDE_IGNORE);
  }
  if (qfo != NULL) {
    FFAddOneString (ffstring, "\n            QfO: ", FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString (ffstring, qfo, FALSE, FALSE, TILDE_IGNORE);
  }
  if (tys != NULL) {
    FFAddOneString (ffstring, "\n            TYS: ", FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString (ffstring, tys, FALSE, FALSE, TILDE_IGNORE);
  }
  if (upr != NULL) {
    FFAddOneString (ffstring, "\n            UPR: ", FALSE, FALSE, TILDE_IGNORE);
    FFAddOneString (ffstring, upr, FALSE, FALSE, TILDE_IGNORE);
  }
}

static CharPtr GetGenomeBuildNumber (
  UserObjectPtr uop
)

{
  ObjectIdPtr   oip;
  CharPtr       str;
  UserFieldPtr  ufp;

  if (uop == NULL) return NULL;
  if ((oip = uop->type) == NULL) return NULL;
  if (StringCmp (oip->str, "GenomeBuild") != 0) return NULL;
  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    oip = ufp->label;
    if (StringCmp(oip->str, "NcbiAnnotation") == 0) {
      if (ufp->choice == 1) { /* string */
        str = ufp->data.ptrvalue;
        if (! StringHasNoText (str)) return str;
      }
    } else if (StringCmp (oip->str, "Annotation") == 0) {
      if (ufp->choice == 1) { /* string */
        str = ufp->data.ptrvalue;
        if (! StringHasNoText (str)) {
          if (StringNICmp (str, "NCBI build ", 11) == 0) {
            if (! StringHasNoText (str + 11)) {
              return (str + 11);
            }
          }
        }
      }
    }
  }
  return NULL;
}

static CharPtr GetGenomeVersionNumber (
  UserObjectPtr uop
)

{
  ObjectIdPtr   oip;
  CharPtr       str;
  UserFieldPtr  ufp;

  if (uop == NULL) return NULL;
  if ((oip = uop->type) == NULL) return NULL;
  if (StringCmp (oip->str, "GenomeBuild") != 0) return NULL;
  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    oip = ufp->label;
    if (StringCmp(oip->str, "NcbiVersion") == 0) {
      if (ufp->choice == 1) { /* string */
        str = ufp->data.ptrvalue;
        if (! StringHasNoText (str)) return str;
      }
    }
  }
  return NULL;
}


static CharPtr reftxt11 = "This record is predicted by automated computational analysis. This record is derived from a genomic sequence";
static CharPtr reftxt12 = "annotated using gene prediction method:";
static CharPtr reftxt13 = "and transcript sequence";

static void FindModelEvidenceUop (
  UserObjectPtr uop,
  Pointer userdata
)

{
  ObjectIdPtr         oip;
  UserObjectPtr PNTR  uopp;

  if (uop == NULL || userdata == NULL) return;
  uopp = (UserObjectPtr PNTR) userdata;
  oip = uop->type;
  if (oip == NULL) return;
  if (StringCmp (oip->str, "ModelEvidence") == 0) {
    *uopp = uop;
  }
}

static Boolean DoGetAnnotationComment (
   BioseqPtr bsp,
   CharPtr PNTR namep,
   UserFieldPtr PNTR assmp,
   BIG_ID_PNTR gip,
   Int4Ptr leftp,
   Int4Ptr rightp,
   CharPtr PNTR methodp,
   BoolPtr mrnaEv,
   BoolPtr estEv
)

{
  UserFieldPtr       assm = NULL;
  Int2               ce = 0, cm = 0;
  SeqMgrDescContext  dcontext;
  BIG_ID             gi = 0;
  Int4               left = 0, right = 0;
  Int4Ptr            ints;
  CharPtr            method = NULL;
  UserObjectPtr      moduop;
  CharPtr            name = NULL;
  ObjectIdPtr        oip;
  SeqDescrPtr        sdp;
  SeqIdPtr           sip;
  TextSeqIdPtr       tsip;
  UserFieldPtr       u;
  UserFieldPtr       ufp;
  UserObjectPtr      uop;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &dcontext);
  while (sdp != NULL) {
    uop = (UserObjectPtr) sdp->data.ptrvalue;
    if (uop != NULL) {
      moduop = NULL;
      VisitUserObjectsInUop (uop, (Pointer) &moduop, FindModelEvidenceUop);
      if (moduop != NULL) {
        oip = moduop->type;
        if (oip != NULL && StringCmp(oip->str, "ModelEvidence") == 0) {
          for (ufp = moduop->data; ufp != NULL; ufp = ufp->next) {
            oip = ufp->label;
            if (oip == NULL) continue;
            if (StringCmp (oip->str, "Contig Name") == 0) {
              name = (CharPtr) ufp->data.ptrvalue;
            } else if (StringCmp (oip->str, "Assembly") == 0) {
              assm = ufp;
            } else if (StringCmp (oip->str, "Contig Gi") == 0) {
              gi = (BIG_ID) ufp->data.intvalue;
            } else if (StringCmp (oip->str, "Contig Span") == 0 && ufp->choice == 8 && ufp->num >= 2) {
              ints = (Int4Ptr) ufp->data.ptrvalue;
              if (ints != NULL) {
                left = ints [0] + 1;
                right = ints [1] + 1;
              }
            } else if (StringCmp (oip->str, "Method") == 0) {
              method = (CharPtr) ufp->data.ptrvalue;
            } else if (StringCmp (oip->str, "mRNA") == 0) {
              *mrnaEv = TRUE;
            } else if (StringCmp (oip->str, "EST") == 0) {
              *estEv = TRUE;
            } else if (StringCmp (oip->str, "Counts") == 0) {
              for (u = (UserFieldPtr) ufp->data.ptrvalue; u != NULL; u = u->next) {
                if (u->data.ptrvalue == NULL) continue;
                if (u->choice != 2) continue;
                oip = u->label;
                if (oip == NULL) continue;
                if (StringCmp (oip->str, "mRNA") == 0) {
                  cm = (Int2) u->data.intvalue;
                  if (cm > 0) {
                    *mrnaEv = TRUE;
                  }
                } else if (StringCmp (oip->str, "EST") == 0) {
                  ce = (Int2) u->data.intvalue;
                  if (ce > 0) {
                    *estEv = TRUE;
                  }
                }
              }
            }
          }
          if (StringHasNoText (name) && bsp != NULL) {
            for (sip = bsp->id; sip != NULL; sip = sip->next) {
              if (sip->choice == SEQID_OTHER) {
                tsip = (TextSeqIdPtr) sip->data.ptrvalue;
                if (tsip != NULL) {
                  name = tsip->accession;
                }
              }
            }
          }
        }
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &dcontext);
  }

  if (StringHasNoText (name)) return FALSE;
  *namep = name;
  *assmp = assm;
  *gip = gi;
  *leftp = left;
  *rightp = right;
  if (! StringHasNoText (method)) {
    *methodp = method;
  }
  return TRUE;
}

static Boolean GetAnnotationComment (
   BioseqPtr bsp,
   CharPtr PNTR namep,
   UserFieldPtr PNTR assmp,
   BIG_ID_PNTR gip,
   Int4Ptr leftp,
   Int4Ptr rightp,
   CharPtr PNTR methodp,
   BoolPtr mrnaEv,
   BoolPtr estEv
)

{
  SeqFeatPtr  cds;

  if (DoGetAnnotationComment (bsp, namep, assmp, gip, leftp, rightp, methodp, mrnaEv, estEv)) return TRUE;
  if (ISA_aa (bsp->mol)) {
    cds = SeqMgrGetCDSgivenProduct (bsp, NULL);
    if (cds != NULL) {
      bsp = BioseqFindFromSeqLoc (cds->location);
      if (bsp != NULL) {
        return DoGetAnnotationComment (bsp, namep, assmp, gip, leftp, rightp, methodp, mrnaEv, estEv);
      }
    }
  }
  return FALSE;
}

static void FindGeneFeat (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  SeqFeatPtr PNTR  sfpp;

  if (sfp->data.choice != SEQFEAT_GENE) return;
  sfpp = (SeqFeatPtr PNTR) userdata;
  *sfpp = sfp;
}

static void FindLocusId (
  ValNodePtr dbxref,
  CharPtr locusIDp
)

{
  DbtagPtr     dbt;
  ObjectIdPtr  oip;
  ValNodePtr   vnp;

  for (vnp = dbxref; vnp != NULL; vnp = vnp->next) {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt == NULL) continue;
    if (StringICmp (dbt->db, "LocusID") != 0 && StringICmp (dbt->db, "InterimID") != 0) continue;
    oip = dbt->tag;
    if (oip == NULL) continue;
    if (oip->str != NULL) {
      StringCpy (locusIDp, oip->str);
    } else if (oip->id > 0) {
      sprintf (locusIDp, "%ld", (long) oip->id);
    }
  }
}

static Boolean GetGeneAndLocus (
  BioseqPtr bsp,
  CharPtr PNTR genep,
  CharPtr locusIDp,
  CharPtr taxIDp
)

{
  BioSourcePtr       biop;
  DbtagPtr           dbt;
  SeqMgrDescContext  dcontext;
  SeqFeatPtr         gene = NULL;
  GeneRefPtr         grp;
  ObjectIdPtr        oip;
  OrgRefPtr          orp;
  SeqDescrPtr        sdp;
  SeqEntryPtr        sep;
  CharPtr            str;
  ValNodePtr         syn;
  ValNodePtr         vnp;

  sep = GetTopSeqEntryForEntityID (bsp->idx.entityID);
  if (sep == NULL) return FALSE;
  VisitFeaturesInSep (sep, (Pointer) &gene, FindGeneFeat);
  if (gene == NULL) return FALSE;

  grp = (GeneRefPtr) gene->data.value.ptrvalue;
  if (grp == NULL) return FALSE;
  if (! StringHasNoText (grp->locus)) {
    *genep = grp->locus;
  } else {
    syn = grp->syn;
    if (syn != NULL) {
      str = (CharPtr) syn->data.ptrvalue;
      if (! StringHasNoText (str)) {
        *genep = str;
      }
    }
  }
  FindLocusId (gene->dbxref, locusIDp);
  FindLocusId (grp->db, locusIDp);

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp != NULL) {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
    if (biop != NULL) {
      orp = biop->org;
      if (orp != NULL) {
        for (vnp = orp->db; vnp != NULL; vnp = vnp->next) {
          dbt = (DbtagPtr) vnp->data.ptrvalue;
          if (dbt == NULL) continue;
          if (StringCmp (dbt->db, "taxon") == 0) {
            oip = dbt->tag;
            if (oip == NULL) continue;
            if (oip->str != NULL) {
              StringCpy (taxIDp, oip->str);
            } else if (oip->id > 0) {
              sprintf (taxIDp, "%ld", (long) oip->id);
            }
          }
        }
      }
    }
  }

  if (genep == NULL || StringHasNoText (locusIDp)) return FALSE;

  return TRUE;
}

static CharPtr nsAreGapsString = "The strings of n's in this record represent gaps between contigs, and the length of each string corresponds to the length of the gap.";
static CharPtr nsWGSGapsString = "The strings of n's in this record represent gaps between contigs or uncallable bases.";

static Boolean IsTpa (
  BioseqPtr bsp,
  Boolean has_tpa_assembly,
  BoolPtr isRefSeqP,
  BoolPtr isTsaP
)

{
  SeqMgrDescContext  dcontext;
  DbtagPtr           dbt;
  Boolean            has_bankit = FALSE;
  Boolean            has_genbank = FALSE;
  Boolean            has_gi = FALSE;
  Boolean            has_local = FALSE;
  Boolean            has_refseq = FALSE;
  Boolean            has_smart = FALSE;
  Boolean            has_tpa = FALSE;
  Boolean            is_tsa = FALSE;
  MolInfoPtr         mip;
  SeqDescrPtr        sdp;
  SeqIdPtr           sip;

  if (bsp == NULL || bsp->id == NULL) return FALSE;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    switch (sip->choice) {
      case SEQID_LOCAL :
        has_local = TRUE;
        break;
      case SEQID_GENBANK :
      case SEQID_EMBL :
      case SEQID_DDBJ :
        has_genbank = TRUE;
        break;
      case SEQID_OTHER :
        has_refseq = TRUE;
        if (isRefSeqP != NULL) {
          *isRefSeqP = TRUE;
        }
        break;
      case SEQID_GI :
        has_gi = TRUE;
        break;
      case SEQID_TPG :
      case SEQID_TPE :
      case SEQID_TPD :
        has_tpa = TRUE;
        break;
      case SEQID_GENERAL :
        dbt = (DbtagPtr) sip->data.ptrvalue;
        if (dbt != NULL) {
          if (StringICmp (dbt->db, "BankIt") == 0) {
            has_bankit = TRUE;
          }
          if (StringICmp (dbt->db, "TMSMART") == 0) {
            has_smart = TRUE;
          }
        }
        break;
      default :
        break;
    }
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
  if (sdp != NULL && sdp->choice == Seq_descr_molinfo) {
    mip = (MolInfoPtr) sdp->data.ptrvalue;
    if (mip != NULL) {
      if (mip->tech == MI_TECH_tsa) {
        is_tsa = TRUE;
        if (isTsaP != NULL) {
          *isTsaP = TRUE;
        }
      }
    }
  }

  if (is_tsa) return TRUE;
  if (has_genbank) return FALSE;
  if (has_tpa) return TRUE;
  if (has_refseq) return TRUE;
  if (has_bankit && has_tpa_assembly) return TRUE;
  if (has_smart && has_tpa_assembly) return TRUE;
  if (has_gi) return FALSE;
  if (has_local && has_tpa_assembly) return TRUE;

  return FALSE;
}

static CharPtr GetPrimaryStrForDelta (
  BioseqPtr bsp
)

{
  Boolean      accn;
  Char         buf [128], tmp [128];
  Int4         curr_start = 0, len, start0, start1;
  DbtagPtr     dbt;
  DeltaSeqPtr  deltasp;
  BIG_ID       gi;
  ValNodePtr   head = NULL;
  SeqIdPtr     id, sip;
  SeqIntPtr    intp;
  SeqLitPtr    litp;
  SeqLocPtr    slp;
  CharPtr      str;
  Uint1        strand;

  if (bsp == NULL || bsp->repr != Seq_repr_delta || bsp->seq_ext_type != 4) return NULL;

  for (deltasp = (DeltaSeqPtr) bsp->seq_ext; deltasp != NULL; deltasp = deltasp->next) {
    if (deltasp->choice == 1) {
      slp = (SeqLocPtr) deltasp->data.ptrvalue;
      if (slp != NULL && slp->choice == SEQLOC_INT) {
        intp = (SeqIntPtr) slp->data.ptrvalue;
        start0 = curr_start;
        start1 = intp->from;
        len = intp->to - intp->from + 1;
        curr_start += len;
        strand = intp->strand;
        sip = intp->id;
        if (sip == NULL) continue;
        id = NULL;
        accn = FALSE;
        if (sip->choice == SEQID_GI) {
          gi = (BIG_ID) sip->data.intvalue;
          if (GetAccnVerFromServer (gi, buf)) {
            accn = TRUE;
          } else {
            id = GetSeqIdForGI (gi);
          }
          if (id == NULL) {
            sprintf (buf, "%ld", (long) gi);
            accn = TRUE;
          }
        } else {
          id = SeqIdDup (sip);
        }
        if (id != NULL || accn) {
          if (head == NULL) {
            ValNodeCopyStr (&head, 0, "CONTIG_SPAN         PRIMARY_IDENTIFIER PRIMARY_SPAN        COMP");
          }
          if (id != NULL) {
            SeqIdWrite (id, buf, PRINTID_TEXTID_ACC_VER, sizeof (buf) - 1);
            if (id->choice == SEQID_GENERAL) {
              dbt = (DbtagPtr) id->data.ptrvalue;
              if (dbt != NULL && StringICmp (dbt->db, "ti") == 0) {
                StringCpy (buf, "TI");
                SeqIdWrite (id, buf + 2, PRINTID_TEXTID_ACC_VER, sizeof (buf) - 3);
              }
            }
          }
          sprintf (tmp, "~%ld-%ld                                        ",
                   (long) (start0 + 1), (long) (start0 + len));
          tmp [21] = '\0';
          StringCat (buf, "                                        ");
          buf [18] = '\0';
          StringCat (tmp, buf);
          sprintf (buf, " %ld-%ld                                        ",
                   (long) (start1 + 1), (long) (start1 + len));
          buf [21] = '\0';
          StringCat (tmp, buf);
          if (strand == Seq_strand_minus) {
            StringCat (tmp, "c");
          }
          ValNodeCopyStr (&head, 0, tmp);
        }
        SeqIdFree (id);
      }
    } else if (deltasp->choice == 2) {
      litp = (SeqLitPtr) deltasp->data.ptrvalue;
      if (litp != NULL) {
        curr_start += litp->length;
      }
    }
  }

  if (head == NULL) return NULL;

  str = MergeFFValNodeStrs (head);
  ValNodeFreeData (head);

  return str;
}

static CharPtr GetStrForTpaOrRefSeqHist (
  BioseqPtr bsp,
  Boolean isRefSeq,
  Boolean isTsa,
  Boolean forcePrimaryBlock
)

{
  Boolean      accn;
  Char         bfr [100];
  Char         buf [100];
  DbtagPtr     dbt;
  BIG_ID       gi;
  ValNodePtr   head = NULL;
  SeqHistPtr   hist;
  SeqIdPtr     id;
  Int2         j;
  int          k;
  Int2         max;
  Boolean      minus1;
  Boolean      minus2;
  Int4         oldstop = -1;
  Uint1        residue;
  SeqAlignPtr  salp;
  SeqAlignPtr  salptmp;
  StreamCache  sc;
  SeqIdPtr     sip;
  Int4         start;
  Int4         stop;
  CharPtr      str;
  Char         tmp [120];

  if (bsp == NULL) return NULL;
  hist = bsp->hist;
  if (hist != NULL && hist->assembly != NULL) {
    salp = SeqAlignListDup (hist->assembly);
    AlnMgr2IndexLite (salp);
    AlnMgr2SortAlnSetByNthRowPos (salp, 1);
    salptmp = (SeqAlignPtr) (salp->segs);
    while (salptmp != NULL) {
      AlnMgr2GetNthSeqRangeInSA (salptmp, 1, &start, &stop);
      sip = AlnMgr2GetNthSeqIdPtr (salptmp, 2);
      if (sip != NULL) {
        id = NULL;
        accn = FALSE;
        buf [0] = '\0';
        if (sip->choice == SEQID_GI) {
          gi = (BIG_ID) sip->data.intvalue;
          if (GetAccnVerFromServer (gi, buf)) {
            accn = TRUE;
          } else {
            id = GetSeqIdForGI (gi);
          }
          if (id == NULL && forcePrimaryBlock) {
            id = SeqIdDup (sip);
          }
        } else {
          id = SeqIdDup (sip);
        }
        if (id != NULL || accn) {
          if (head == NULL) {
            if (isRefSeq) {
              ValNodeCopyStr (&head, 0, "REFSEQ_SPAN         PRIMARY_IDENTIFIER PRIMARY_SPAN        COMP");
            } else if (isTsa) {
              ValNodeCopyStr (&head, 0, "TSA_SPAN            PRIMARY_IDENTIFIER PRIMARY_SPAN        COMP");
            } else {
              ValNodeCopyStr (&head, 0, "TPA_SPAN            PRIMARY_IDENTIFIER PRIMARY_SPAN        COMP");
            }
          }
          if (isRefSeq && oldstop > -1 && oldstop < start) {
            sprintf (tmp, "~%ld-%ld                                        ",
                     (long) (oldstop + 1), (long) (start));
            tmp [21] = '\0';
            StringCpy (bfr, "                                        ");
            k = 0;
            if (StreamCacheSetup (bsp, NULL, 0, &sc)) {
              if (start - oldstop < 15) {
                StreamCacheSetPosition (&sc, oldstop);
                bfr [k] = '"';
                k++;
                max = start - oldstop;
                for (j = 0; j < max; j++) {
                  residue = StreamCacheGetResidue (&sc);
                  bfr [k] = (Char) residue;
                  k++;
                }
                bfr [k] = '"';
                k++;
              } else {
                StreamCacheSetPosition (&sc, oldstop);
                bfr [k] = '"';
                k++;
                for (j = 0; j < 4; j++) {
                  residue = StreamCacheGetResidue (&sc);
                  bfr [k] = (Char) residue;
                  k++;
                }
                bfr [k] = '.';
                k++;
                bfr [k] = '.';
                k++;
                bfr [k] = '.';
                k++;
                StreamCacheSetPosition (&sc, start - 4);
                for (j = 0; j < 4; j++) {
                  residue = StreamCacheGetResidue (&sc);
                  bfr [k] = (Char) residue;
                  k++;
                }
                bfr [k] = '"';
                k++;
              }
            } else {
              /*
              StringCpy (bfr, "inserted base(s)");
              */
            }
            bfr [k] = '\0';
            StringCat (bfr, "                                        ");
            bfr [18] = '\0';
            StringCat (tmp, bfr);
            sprintf (bfr, " %ld-%ld                                        ",
                     (long) 1, (long) (start - oldstop));
            bfr [21] = '\0';
            StringCat (tmp, bfr);
            ValNodeCopyStr (&head, 0, tmp);
          }
          oldstop = stop + 1;
          if (id != NULL) {
            SeqIdWrite (id, buf, PRINTID_TEXTID_ACC_VER, sizeof (buf) - 1);
            if (id->choice == SEQID_GENERAL) {
              dbt = (DbtagPtr) id->data.ptrvalue;
              if (dbt != NULL && StringICmp (dbt->db, "ti") == 0) {
                StringCpy (buf, "TI");
                SeqIdWrite (id, buf + 2, PRINTID_TEXTID_ACC_VER, sizeof (buf) - 3);
              }
            }
          }
          sprintf (tmp, "~%ld-%ld                                        ",
                   (long) (start + 1), (long) (stop + 1));
          /*
          i = 39 - StringLen (buf);
          if (i > 0) {
            tmp [i] = '\0';
          } else {
            tmp [21] = '\0';
          }
          */
          tmp [21] = '\0';
          StringCat (buf, "                                        ");
          buf [18] = '\0';
          StringCat (tmp, buf);
          AlnMgr2GetNthSeqRangeInSA (salptmp, 2, &start, &stop);
          sprintf (buf, " %ld-%ld                                        ",
                   (long) (start + 1), (long) (stop + 1));
          buf [21] = '\0';
          StringCat (tmp, buf);
          minus1 = (Boolean) (AlnMgr2GetNthStrand (salptmp, 1) == Seq_strand_minus);
          minus2 = (Boolean) (AlnMgr2GetNthStrand (salptmp, 2) == Seq_strand_minus);
          if (minus1 || minus2) {
            if (! (minus1 && minus2)) {
              StringCat (tmp, "c");
            }
          }
          ValNodeCopyStr (&head, 0, tmp);
        }
        SeqIdFree (id);
      }
      SeqIdFree (sip);
      salptmp = salptmp->next;
    }
    SeqAlignFree (salp);
  }

  if (head == NULL) return NULL;

  str = MergeFFValNodeStrs (head);
  ValNodeFreeData (head);

  return str;
}

static CharPtr tpaString = "THIRD PARTY ANNOTATION DATABASE: This TPA record uses data from DDBJ/EMBL/GenBank ";

static CharPtr GetStrForTPA (
  UserObjectPtr uop,
  BioseqPtr bsp
)

{
  Char          ch;
  UserFieldPtr  curr;
  SeqHistPtr    hist;
  Int2          i;
  Char          id [41];
  Boolean       isRefSeq = FALSE;
  Boolean       isTsa = FALSE;
  Int2          j;
  size_t        len;
  ObjectIdPtr   oip;
  CharPtr       ptr;
  CharPtr       str;
  CharPtr       tmp;
  UserFieldPtr  ufp;

  if (uop == NULL) return NULL;
  if ((oip = uop->type) == NULL) return NULL;
  if (StringCmp (oip->str, "TpaAssembly") != 0) return NULL;
  if (bsp == NULL) return NULL;
  hist = bsp->hist;
  if (hist != NULL && hist->assembly != NULL) return NULL;
  if (! IsTpa (bsp, TRUE, &isRefSeq, &isTsa)) return NULL;
  if (isRefSeq) return NULL;

  len = StringLen (tpaString) + StringLen ("entries ") + StringLen ("and ") + 5;
  i = 0;
  for (curr = uop->data; curr != NULL; curr = curr->next) {
    if (curr->choice != 11) continue;
    for (ufp = curr->data.ptrvalue; ufp != NULL; ufp = ufp->next) {
      if (ufp->choice != 1) continue;
      oip = ufp->label;
      if (oip == NULL || StringICmp (oip->str, "accession") != 0) continue;
      str = (CharPtr) ufp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      len += StringLen (str) + 2;
      i++;
    }
  }
  if (i == 0) return NULL;

  ptr = (CharPtr) MemNew (len);
  if (ptr == NULL) return NULL;
  StringCpy (ptr, tpaString);
  if (i > 1) {
    StringCat (ptr, "entries ");
  } else {
    StringCat (ptr, "entry ");
  }

  j = 0;
  for (curr = uop->data; curr != NULL; curr = curr->next) {
    if (curr->choice != 11) continue;
    for (ufp = curr->data.ptrvalue; ufp != NULL; ufp = ufp->next) {
      if (ufp->choice != 1) continue;
      oip = ufp->label;
      if (oip == NULL || StringICmp (oip->str, "accession") != 0) continue;
      str = (CharPtr) ufp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      StringNCpy_0 (id, str, sizeof (id));
      tmp = id;
      ch = *tmp;
      while (ch != '\0') {
        if (IS_LOWER (ch)) {
          *tmp = TO_UPPER (ch);
        }
        tmp++;
        ch = *tmp;
      }
      if (j == i - 1 && i > 1) {
        StringCat (ptr, " and ");
      } else if (j > 0) {
        StringCat (ptr, ", ");
      }
      StringCat (ptr, id);
      j++;
    }
  }

  return ptr;
}

static CharPtr GetStrForGenome (
  UserObjectPtr uop,
  BioseqPtr bsp
)

{
  ObjectIdPtr  oip;

  if (uop == NULL) return NULL;
  if ((oip = uop->type) == NULL) return NULL;
  if (StringCmp (oip->str, "GenomeInfo") != 0) return NULL;

  /* !!! need to implement !!! */

  return NULL;
}

static void AddAltPrimaryBlock (
  Asn2gbWorkPtr awp
)

{
  IntAsn2gbJobPtr  ajp;
  Asn2gbSectPtr    asp;
  BaseBlockPtr     bbp = NULL;
  BioseqPtr        bsp;
  GBSeqPtr         gbseq;
  CharPtr          str;
  StringItemPtr    ffstring;

  if (awp == NULL) return;
  ajp = awp->ajp;
  if (ajp == NULL) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;
  asp = awp->asp;
  if (asp == NULL) return;

  ffstring = FFGetString(ajp);
  if ( ffstring == NULL ) return;

  str = GetPrimaryStrForDelta (bsp);
  if (str != NULL) {

    bbp = (BaseBlockPtr) Asn2gbAddBlock (awp, PRIMARY_BLOCK, sizeof (BaseBlock));
    if (bbp != NULL) {

      FFStartPrint (ffstring, awp->format, 0, 12, "PRIMARY", 12, 5, 5, "PR", TRUE);

      FFAddOneString (ffstring, str, FALSE, FALSE, TILDE_EXPAND);

      bbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "PR");

      /* optionally populate gbseq for XML-ized GenBank format */

      if (ajp->gbseq) {
        gbseq = &asp->gbseq;
      } else {
        gbseq = NULL;
      }

      if (gbseq != NULL) {
        gbseq->primary = StringSave (str);
      }

      if (awp->afp != NULL) {
        DoImmediateFormat (awp->afp, (BaseBlockPtr) bbp);
      }
    }
    MemFree (str);
  }

  FFRecycleString(ajp, ffstring);
}

static CharPtr GeStrForTSA (
  UserObjectPtr uop
)

{
  Int4          asf, ast, prf, prt;
  Char          buf [128], tmp [128];
  UserFieldPtr  curr;
  Boolean       has_asf, has_ast, has_prf, has_prt;
  ValNodePtr    head = NULL;
  ObjectIdPtr   oip;
  CharPtr       pid;
  CharPtr       str;
  UserFieldPtr  ufp;

  if (uop == NULL) return NULL;
  if ((oip = uop->type) == NULL) return NULL;
  if (StringCmp (oip->str, "TSA") != 0) return NULL;

  for (curr = uop->data; curr != NULL; curr = curr->next) {
    if (curr->choice != 11) continue;
    asf = 0;
    ast = 0;
    prf = 0;
    prt = 0;
    pid = NULL;
    has_asf = FALSE;
    has_ast = FALSE;
    has_prf = FALSE;
    has_prt = FALSE;
    for (ufp = curr->data.ptrvalue; ufp != NULL; ufp = ufp->next) {
      oip = ufp->label;
      if (oip == NULL) continue;
      if (StringICmp (oip->str, "assembly from") == 0 && ufp->choice == 2) {
        asf = (Int4) ufp->data.intvalue;
        has_asf = TRUE;
      } else if (StringICmp (oip->str, "assembly to") == 0 && ufp->choice == 2) {
        ast = (Int4) ufp->data.intvalue;
        has_ast = TRUE;
      } else if (StringICmp (oip->str, "primary from") == 0 && ufp->choice == 2) {
        prf = (Int4) ufp->data.intvalue;
        has_prf = TRUE;
      } else if (StringICmp (oip->str, "primary to") == 0 && ufp->choice == 2) {
        prt = (Int4) ufp->data.intvalue;
        has_prt = TRUE;
      } else if (StringICmp (oip->str, "primary ID") == 0 && ufp->choice == 1) {
        pid = (CharPtr) ufp->data.ptrvalue;
      }
    }
    if (has_asf && has_ast && has_prf && has_prt && pid != NULL) {
      if (head == NULL) {
        ValNodeCopyStr (&head, 0, "TSA_SPAN            PRIMARY_IDENTIFIER PRIMARY_SPAN        COMP");
      }
      StringCpy (buf, pid);
      if (StringNCmp (pid, "gnl|ti|", 7) == 0) {
        StringCpy (buf, "TI");
        StringCat (buf, pid + 7);
      }
      sprintf (tmp, "~%ld-%ld                                        ",
               (long) (asf + 1), (long) (ast + 1));
      tmp [21] = '\0';
      StringCat (buf, "                                        ");
      buf [18] = '\0';
      StringCat (tmp, buf);
      sprintf (buf, " %ld-%ld                                        ",
               (long) (prf + 1), (long) (prt + 1));
      buf [21] = '\0';
      StringCat (tmp, buf);
      if (prf > prt) {
        StringCat (tmp, "c");
      }
      ValNodeCopyStr (&head, 0, tmp);
    }
  }

  if (head == NULL) return NULL;

  str = MergeFFValNodeStrs (head);
  ValNodeFreeData (head);

  return str;
}

static void AddTsaBlock (
  Asn2gbWorkPtr awp,
  UserObjectPtr uop
)

{
  IntAsn2gbJobPtr  ajp;
  Asn2gbSectPtr    asp;
  BaseBlockPtr     bbp = NULL;
  BioseqPtr        bsp;
  GBSeqPtr         gbseq;
  CharPtr          str;
  StringItemPtr    ffstring;

  if (awp == NULL) return;
  ajp = awp->ajp;
  if (ajp == NULL) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;
  asp = awp->asp;
  if (asp == NULL) return;

  ffstring = FFGetString(ajp);
  if ( ffstring == NULL ) return;

  str = GeStrForTSA (uop);
  if (str != NULL) {

    bbp = (BaseBlockPtr) Asn2gbAddBlock (awp, PRIMARY_BLOCK, sizeof (BaseBlock));
    if (bbp != NULL) {

      FFStartPrint (ffstring, awp->format, 0, 12, "PRIMARY", 12, 5, 5, "PR", TRUE);

      FFAddOneString (ffstring, str, FALSE, FALSE, TILDE_EXPAND);

      bbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "PR");

      /* optionally populate gbseq for XML-ized GenBank format */

      if (ajp->gbseq) {
        gbseq = &asp->gbseq;
      } else {
        gbseq = NULL;
      }

      if (gbseq != NULL) {
        gbseq->primary = StringSave (str);
      }

      if (awp->afp != NULL) {
        DoImmediateFormat (awp->afp, (BaseBlockPtr) bbp);
      }
    }
    MemFree (str);
  }

  FFRecycleString(ajp, ffstring);
}

NLM_EXTERN void AddPrimaryBlock (
  Asn2gbWorkPtr awp
)

{
  IntAsn2gbJobPtr    ajp;
  Asn2gbSectPtr      asp;
  BaseBlockPtr       bbp = NULL;
  BioseqPtr          bsp;
  SeqMgrDescContext  dcontext;
  GBSeqPtr           gbseq;
  Boolean            has_tpa_assembly = FALSE;
  Boolean            has_tsa = FALSE;
  SeqHistPtr         hist;
  Boolean            isRefSeq = FALSE;
  Boolean            isTsa = FALSE;
  ObjectIdPtr        oip;
  SeqDescrPtr        sdp;
  CharPtr            str;
  UserObjectPtr      uop;
  StringItemPtr      ffstring;

  if (awp == NULL) return;
  ajp = awp->ajp;
  if (ajp == NULL) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;
  asp = awp->asp;
  if (asp == NULL) return;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &dcontext);
  while (sdp != NULL) {
    uop = (UserObjectPtr) sdp->data.ptrvalue;
    if (uop != NULL) {
      oip = uop->type;
      if (oip != NULL) {
        if (StringCmp (oip->str, "TpaAssembly") == 0) {
          has_tpa_assembly = TRUE;
        } else if (StringCmp (oip->str, "TSA") == 0) {
          has_tsa = TRUE;
        }
      }
    }
    if (has_tpa_assembly || has_tsa) {
      sdp = NULL;
    } else {
      sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &dcontext);
    }
  }

  if (has_tsa) {
    AddTsaBlock (awp, uop);
    return;
  }

  hist = bsp->hist;
  if ((! IsTpa (bsp, has_tpa_assembly, &isRefSeq, &isTsa)) ||
      hist == NULL || hist->assembly == NULL) {
    if (awp->forcePrimaryBlock) {
      AddAltPrimaryBlock (awp);
    }
    return;
  }

  ffstring = FFGetString(ajp);
  if ( ffstring == NULL ) return;

  str = GetStrForTpaOrRefSeqHist (bsp, isRefSeq, isTsa, awp->forcePrimaryBlock);
  if (str != NULL) {

    bbp = (BaseBlockPtr) Asn2gbAddBlock (awp, PRIMARY_BLOCK, sizeof (BaseBlock));
    if (bbp != NULL) {

      if (has_tpa_assembly) {
        bbp->entityID = dcontext.entityID;
        bbp->itemID = dcontext.itemID;
        bbp->itemtype = OBJ_SEQDESC;
      }

      FFStartPrint (ffstring, awp->format, 0, 12, "PRIMARY", 12, 5, 5, "PR", TRUE);

      FFAddOneString (ffstring, str, FALSE, FALSE, TILDE_EXPAND);

      bbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "PR");

      /* optionally populate gbseq for XML-ized GenBank format */

      if (ajp->gbseq) {
        gbseq = &asp->gbseq;
      } else {
        gbseq = NULL;
      }

      if (gbseq != NULL) {
        gbseq->primary = StringSave (str);
      }

      if (awp->afp != NULL) {
        DoImmediateFormat (awp->afp, (BaseBlockPtr) bbp);
      }
    }
    MemFree (str);
  }

  FFRecycleString(ajp, ffstring);
}

static CharPtr reftxt32 = "It is defined by coordinates on the sequence of chromosome";
static CharPtr reftxt33 = "from the";
static CharPtr reftxt34 = "assembly of the human genome (NCBI build";
static CharPtr reftxt35 = ").";

static CharPtr GetEncodeString (
  UserObjectPtr uop,
  BioseqPtr bsp
)

{
  CharPtr            assembly_date = NULL;
  BioSourcePtr       biop;
  CharPtr            chromosome = NULL;
  SeqMgrDescContext  dcontext;
  size_t             len;
  CharPtr            ncbi_annotation = NULL;
  ObjectIdPtr        oip;
  SeqDescrPtr        sdp;
  SubSourcePtr       ssp;
  CharPtr            str;
  UserFieldPtr       ufp;

  if (uop == NULL || bsp == NULL) return NULL;

  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    oip = ufp->label;
    if (oip == NULL || oip->str == NULL || ufp->choice != 1) continue;
    if (StringICmp (oip->str, "AssemblyDate") == 0) {
      assembly_date = (CharPtr) ufp->data.ptrvalue;
    } else if (StringICmp (oip->str, "NcbiAnnotation") == 0) {
      ncbi_annotation = (CharPtr) ufp->data.ptrvalue;
    }
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp != NULL) {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
    if (biop != NULL) {
      for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
        if (ssp->subtype == SUBSRC_chromosome) {
          chromosome = ssp->name;
        }
      }
    }
  }

  if (chromosome == NULL || assembly_date == NULL || ncbi_annotation == NULL) return NULL;

  if (StringHasNoText (chromosome)) {
    chromosome = "?";
  }
  if (StringHasNoText (assembly_date)) {
    assembly_date = "?";
  }
  if (StringHasNoText (ncbi_annotation)) {
    ncbi_annotation = "?";
  }

  len = StringLen (reftxt32) + StringLen (reftxt33) +
        StringLen (reftxt34) + StringLen (reftxt35) +
        StringLen (chromosome) +
        StringLen (assembly_date) +
        StringLen (ncbi_annotation);

  str = (CharPtr) MemNew (sizeof (Char) * (len + 10));
  if (str == NULL) return NULL;

  sprintf (str, "%s %s %s %s %s %s%s", reftxt32, chromosome, reftxt33,
           assembly_date, reftxt34, ncbi_annotation, reftxt35);

  return str;
}


typedef struct unverifiedtypeinfodata {
  CharPtr match_name;
  CharPtr comment_text;
} UnverifiedTypeInfoData, PNTR UnverifiedTypeInfoPtr;


static UnverifiedTypeInfoData s_UnverifiedTypeInfo[] = {
  { "Organism", "source organism" },
  { "Features", "sequence and/or annotation" },
  { "Misassembled", "sequence assembly" }
};


NLM_EXTERN CharPtr GetUnverifiedMatchName (Int4 unverified_type)
{
  if (unverified_type < 0 || unverified_type > eUnverifiedType_Max) {
    return NULL;
  } else {
    return s_UnverifiedTypeInfo[unverified_type].match_name;
  }
}


static void GetUnverifiedFlags (UserObjectPtr uop, BoolPtr unverified_flags)
{
  Int4 i;
  UserFieldPtr ufp;
  ObjectIdPtr  oip;
  CharPtr      str;
  Boolean any = FALSE;

  if (unverified_flags == NULL) {
    return;
  }
  for (i = 0; i < eUnverifiedType_Max; i++) {
    unverified_flags[i] = FALSE;
  }
  if (uop == NULL) {
    return;
  }

  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    oip = ufp->label;
    if (oip != NULL && StringCmp (oip->str, "Type") == 0 && ufp->choice == 1) {
      str = (CharPtr) ufp->data.ptrvalue;
      for (i = 0; i < eUnverifiedType_Max; i++) {
        if (StringICmp (str, s_UnverifiedTypeInfo[i].match_name) == 0) {
          unverified_flags[i] = TRUE;
          any = TRUE;
          break;
        }
      }
    }
  }
  if (!any) {
    /* default in the past was to use feature if not source */
    unverified_flags[eUnverifiedType_Features] = TRUE;
  }
}


static CharPtr CommentTextFromUnverifiedFlags(BoolPtr unverified_flags)
{
  Int4 i, len, num_items = 0, item;
  CharPtr comment_start = "GenBank staff is unable to verify ";
  CharPtr comment_end = " provided by the submitter.";
  CharPtr and = "and ";
  CharPtr comma = ", ";
  CharPtr comment = NULL;

  if (unverified_flags == NULL) {
    return NULL;
  }

  len = StringLen (comment_start) + StringLen (comment_end) + 1;
  for (i = 0; i < eUnverifiedType_Max; i++) {
    if (unverified_flags[i]) {
      num_items++;
      len += StringLen (s_UnverifiedTypeInfo[i].comment_text);
    }
  }
  if (num_items > 1) {
    len += StringLen (and);
    if (num_items > 2) {
      len += StringLen (comma) * (num_items - 1);
    } else {
      len += 1;
    }
  } else if (num_items == 0) {
    return NULL;
  }

  comment = (CharPtr) MemNew (sizeof (Char) * len);
  StringCpy (comment, comment_start);
  item = 0;
  for (i = 0; i < eUnverifiedType_Max; i++) {
    if (unverified_flags[i]) {
      if (item > 0) {
        if (num_items > 2) {
          StringCat (comment, comma);
        }
        if (item == num_items - 1) {
          if (num_items == 2) {
            StringCat (comment, " ");
          }
          StringCat (comment, and);
        }
      }
      StringCat (comment, s_UnverifiedTypeInfo[i].comment_text);
      item++;
    }
  }
  StringCat (comment, comment_end);
  return comment;
}

static Int4 GetFileTrackPoint (SeqPntPtr spp, PackSeqPntPtr psp, Int4 index)

{
  if (spp != NULL) {
    return spp->point;
  } else if (psp != NULL) {
    return PackSeqPntGet (psp, index);
  }
  return 0;
}  

static Boolean CommentsAreDifferent (CharPtr str, CharPtr last_name)

{
  size_t  lens, lenl;

  if (str == NULL && last_name == NULL) return FALSE;

  if (StringCmp (str, last_name) == 0) return FALSE;

  lens = StringLen (str);
  lenl = StringLen (last_name);

  if (lens == lenl + 1) {
    if (StringNCmp (str, last_name, lenl) == 0) {
      if (str [lens - 1] == '.') {
        return FALSE;
      }
    }
  } else if (lenl == lens + 1) {
    if (StringNCmp (str, last_name, lens) == 0) {
      if (last_name [lenl - 1] == '.') {
        return FALSE;
      }
    }
  }

  return TRUE;
}

NLM_EXTERN void AddCommentBlock (
  Asn2gbWorkPtr awp
)

{
  size_t             acclen;
  CharPtr            accn;
  SeqMgrAndContext   acontext;
  AnnotDescPtr       adp;
  Boolean            annotDescCommentToComment = FALSE;
  IntAsn2gbJobPtr    ajp;
  UserFieldPtr       assm = NULL;
  CharPtr            authaccessvalue = NULL;
  Int4               authaccess_itemID = 0;
  BioseqPtr          bsp;
  Char               buf [2048];
  CommentBlockPtr    cbp;
  Char               ch;
  Int2               chunk;
  Int2               count;
  CharPtr PNTR       cpp;
  Boolean            didGenome = FALSE;
  Boolean            didRefTrack = FALSE;
  Boolean            didTPA = FALSE;
  DbtagPtr           dbt;
  SeqMgrDescContext  dcontext;
  DeltaSeqPtr        dsp;
  UserObjectPtr      encodeUop = NULL;
  Boolean            estEv = FALSE;
  BioseqPtr          farbsp;
  Uint2              fareid;
  /*
  SeqMgrFeatContext  fcontext;
  */
  CharPtr            field;
  PackSeqPntPtr      filetrackpsp = NULL;
  SeqPntPtr          filetrackspp = NULL;
  CharPtr            filetrackURL = NULL;
  Int4               basemodNum = 0;
  CharPtr PNTR       basemodURLhead = NULL;
  CharPtr            basemodURL = NULL;
  Int4               filetrack_itemID = 0;
  Boolean            first = TRUE;
  UserObjectPtr      firstGenAnnotSCAD = NULL;
  CharPtr            firstGenAnnotSCStr = NULL;
  Int4               frags;
  GBBlockPtr         gbp;
  CharPtr            geneName = NULL;
  CharPtr            genomeBuildNumber = NULL;
  CharPtr            genomeVersionNumber = NULL;
  BIG_ID             gi = 0;
  Int4               gsdbid = 0;
  /*
  Boolean            has_gaps = FALSE;
  */
  Boolean            hasRefTrackStatus = FALSE;
  SeqHistPtr         hist;
  Int4               idx;
  Boolean            is_collab = FALSE;
  Boolean            is_encode = FALSE;
  Boolean            is_other = FALSE;
  Boolean            is_tpa = FALSE;
  Boolean            is_wgs = FALSE;
  Boolean            isRefSeqStandard = FALSE;
  Boolean            is_unverified = FALSE;
  Int4               j;
  Int4               last;
  Boolean            last_had_tilde = FALSE;
  CharPtr            last_name;
  Int4               left;
  size_t             len;
  /*
  SeqLitPtr          litp;
  */
  ObjectIdPtr        localID = NULL;
  Char               locusID [32];
  CharPtr            method = NULL;
  MolInfoPtr         mip;
  Boolean            mrnaEv = FALSE;
  SeqIdPtr           msip;
  CharPtr            name = NULL;
  ObjectIdPtr        ncbifileID = NULL;
  CharPtr            nm;
  Int4               num;
  ObjectIdPtr        oip;
  Boolean            okay;
  CharPtr            origLocalID = NULL;
  /*
  BioseqPtr          parent;
  */
  CharPtr            pfx;
  CharPtr            plural;
  Int4               pos;
  Int4               right;
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp;
  CharPtr            sfx;
  Boolean            showedLocalID = FALSE;
  Boolean            showGBBSource = FALSE;
  SeqIdPtr           sip;
  SeqLocPtr          slp;
  CharPtr            str;
  Char               taxID [64];
  CharPtr            tlsaccn = NULL;
  CharPtr            tlsname = NULL;
  Char               tmp [128];
  CharPtr            tsaaccn = NULL;
  CharPtr            tsaname = NULL;
  TextSeqIdPtr       tsip;
  TextSeqIdPtr       tlstsip = NULL;
  UserFieldPtr       tufp;
  UserFieldPtr       ufp;
  Boolean            unordered = FALSE;
  Int4               unverified_itemID = 0;
  UserObjectPtr      uop;
  Int4               version;
  ValNodePtr         vnp;
  CharPtr            wgsaccn = NULL;
  CharPtr            wgsname = NULL;
  StringItemPtr      ffstring = NULL;
  Boolean            unverified_flags[eUnverifiedType_Max];
  CharPtr            unverified_comment;

  if (awp == NULL) return;
  ajp = awp->ajp;
  if (ajp == NULL) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;

  if (GetWWW (ajp) && awp->mode == ENTREZ_MODE && awp->afp != NULL &&
      (awp->format == GENBANK_FMT || awp->format == GENPEPT_FMT)) {
    sprintf (buf, "<a name=\"comment_%s\"></a>", awp->currAccVerLabel);
    DoQuickLinkFormat (awp->afp, buf);
  }

  ffstring = FFGetString(ajp);
  if ( ffstring ==  NULL ) return;

  GetUnverifiedFlags(NULL, unverified_flags);
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &dcontext);
  while (sdp != NULL) {
    uop = (UserObjectPtr) sdp->data.ptrvalue;
    if (uop != NULL) {
      str = GetStatusForRefTrack (uop);
      if (str != NULL) {
        hasRefTrackStatus = TRUE;
      }
      if (genomeBuildNumber == NULL) {
        genomeBuildNumber = GetGenomeBuildNumber (uop);
      }
      if (genomeVersionNumber == NULL) {
        genomeVersionNumber = GetGenomeVersionNumber (uop);
      }
      oip = uop->type;
      if (oip != NULL) {
        if (StringICmp (oip->str, "Unverified") == 0) {
          is_unverified = TRUE;
          unverified_itemID = dcontext.itemID;
          GetUnverifiedFlags(uop, unverified_flags);
        }
        if (StringICmp (oip->str, "ENCODE") == 0) {
          is_encode = TRUE;
          encodeUop = uop;
        }
        if (StringICmp (oip->str, "FileTrack") == 0) {
          filetrack_itemID = dcontext.itemID;
          for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
            oip = ufp->label;
            if (oip == NULL) continue;
            if (StringCmp (oip->str, "FileTrackURL") == 0 || StringCmp (oip->str, "Map-FileTrackURL") == 0) {
              if (ufp->choice == 1 && ufp->data.ptrvalue != NULL) {
                filetrackURL = (CharPtr) ufp->data.ptrvalue;
              } else if (ufp->choice == 7 && ufp->data.ptrvalue != NULL && ufp->num > 0) {
                cpp = (CharPtr PNTR) ufp->data.ptrvalue;
                if (cpp != NULL) {
                  filetrackURL = cpp [0];
                }
              }
            } else if (StringCmp (oip->str, "BaseModification-FileTrackURL") == 0) {
              if (ufp->choice == 1 && ufp->data.ptrvalue != NULL) {
                basemodURL = (CharPtr) ufp->data.ptrvalue;
                basemodNum = 1;
              } else if (ufp->choice == 7 && ufp->data.ptrvalue != NULL && ufp->num > 0) {
                cpp = (CharPtr PNTR) ufp->data.ptrvalue;
                if (cpp != NULL) {
                  basemodURLhead = cpp;
                  basemodNum = ufp->num;
                }
              }
            }
          }
        }
        if (StringICmp (oip->str, "AuthorizedAccess") == 0) {
          authaccess_itemID = dcontext.itemID;
          for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
            oip = ufp->label;
            if (oip == NULL) continue;
            if (StringCmp (oip->str, "Study") != 0) continue;
            if (ufp->choice != 1 || ufp->data.ptrvalue == NULL) continue;
            authaccessvalue = (CharPtr) ufp->data.ptrvalue;
          }
        }
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &dcontext);
  }

  if (is_unverified) {
    cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
    if (cbp != NULL) {

      cbp->entityID = awp->entityID;
      cbp->itemID = unverified_itemID;
      cbp->itemtype = OBJ_SEQDESC;
      cbp->first = first;
      cbp->no_blank_before = last_had_tilde;
      first = FALSE;

      if (cbp->first) {
        FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
      } else {
        FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
      }

      unverified_comment = CommentTextFromUnverifiedFlags(unverified_flags);
      if (unverified_comment != NULL) {
          FFAddOneString (ffstring,
                          unverified_comment,
                          FALSE, FALSE, TILDE_IGNORE);
          unverified_comment = MemFree (unverified_comment);
      }

      cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
      FFRecycleString(ajp, ffstring);
      ffstring = FFGetString(ajp);

      last_had_tilde = FALSE;
      if (awp->afp != NULL) {
        DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
      }
    }
  }

  if (bsp->repr == Seq_repr_map && bsp->seq_ext_type == 3) {
    for (sfp = (SeqFeatPtr) bsp->seq_ext; sfp != NULL; sfp = sfp->next) {
      if (sfp->data.choice != SEQFEAT_RSITE) continue;
      slp = sfp->location;
      if (slp == NULL) continue;
      if (slp->choice == SEQLOC_PNT) {
        filetrackspp = (SeqPntPtr) slp->data.ptrvalue;
      } else if (slp->choice == SEQLOC_PACKED_PNT) {
        filetrackpsp = (PackSeqPntPtr) slp->data.ptrvalue;
      }
    }
  }

  if (authaccessvalue != NULL) {
    cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
    if (cbp != NULL) {

      cbp->entityID = awp->entityID;
      cbp->itemID = authaccess_itemID;
      cbp->itemtype = OBJ_SEQDESC;
      cbp->first = first;
      cbp->no_blank_before = last_had_tilde;
      first = FALSE;

      if (cbp->first) {
        FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
      } else {
        FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
      }

      FFAddOneString (ffstring, "These data are available through the dbGaP authorized access system. ", FALSE, FALSE, TILDE_IGNORE);
      if (GetWWW (ajp)) {
        FFAddOneString (ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, "https://dbgap.ncbi.nlm.nih.gov/aa/wga.cgi?adddataset=", FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, authaccessvalue, FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, "&page=login", FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, "\">", FALSE, FALSE, TILDE_IGNORE); 
        FFAddOneString (ffstring, "Request access", FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, " to Study ", FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, "https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=", FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, authaccessvalue, FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, "\">", FALSE, FALSE, TILDE_IGNORE); 
        FFAddOneString (ffstring, authaccessvalue, FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
      } else {
        FFAddOneString (ffstring, "Request access to Study ", FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, authaccessvalue, FALSE, FALSE, TILDE_IGNORE);
      }
      FFAddOneString (ffstring, ".", FALSE, FALSE, TILDE_IGNORE);

      cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
      FFRecycleString(ajp, ffstring);
      ffstring = FFGetString(ajp);

      last_had_tilde = FALSE;
      if (awp->afp != NULL) {
        DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
      }
    }
  }

  /*
  look for Seq-annot.desc.comment on annots packaged on current bioseq,
  Genome-Annotation structured comment will suppress GenomeBuild user object
  */

  adp = SeqMgrGetNextAnnotDesc (bsp, NULL, Annot_descr_user, &acontext);
  while (adp != NULL) {
    uop = (UserObjectPtr) adp->data.ptrvalue;
    if (uop != NULL) {
      oip = uop->type;
      if (oip != NULL) {
        if (StringCmp (oip->str, "AnnotDescCommentPolicy") == 0) {
          for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
            oip = ufp->label;
            if (oip == NULL || ufp->data.ptrvalue == NULL) continue;
            if (StringCmp (oip->str, "Policy") == 0) {
              if (StringICmp ((CharPtr) ufp->data.ptrvalue, "ShowInComment") == 0) {
                annotDescCommentToComment = TRUE;
              }
            }
          }
        } else if (StringICmp (oip->str, "StructuredComment") == 0) {
          if (firstGenAnnotSCAD == NULL) {
            for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
              if (ufp->choice != 1) continue;
              oip = ufp->label;
              if (oip == NULL) continue;
              field = oip->str;
              if (StringHasNoText (field)) continue;
              if (StringCmp (field, "StructuredCommentPrefix") == 0) {
                if (StringCmp ((CharPtr) ufp->data.ptrvalue, "##Genome-Annotation-Data-START##") == 0) {
                  firstGenAnnotSCAD = uop;
                  genomeBuildNumber = NULL;
                  genomeVersionNumber = NULL;
                  firstGenAnnotSCStr = GetStrForStructuredComment (ajp, firstGenAnnotSCAD);
                }
              }
            }
          }
        }
      }
    }
    adp = SeqMgrGetNextAnnotDesc (bsp, adp, Annot_descr_user, &acontext);
  }

  /*
  also look on first far sequence component of NCBI_GENOMES records
  */

  if (awp->isNCBIGenomes && firstGenAnnotSCAD == NULL && bsp->repr == Seq_repr_delta && bsp->seq_ext_type == 4) {
    for (dsp = (DeltaSeqPtr) (bsp->seq_ext); dsp != NULL; dsp = dsp->next) {
      if (dsp->choice != 1) continue;
      slp = (SeqLocPtr) dsp->data.ptrvalue;
      if (slp == NULL) continue;
      sip = SeqLocId (slp);
      if (sip == NULL) continue;
      farbsp = BioseqLockById (sip);
      if (farbsp == NULL) break;
      fareid = ObjMgrGetEntityIDForPointer (farbsp);
      SeqMgrIndexFeatures (fareid, NULL);
      adp = SeqMgrGetNextAnnotDesc (farbsp, NULL, Annot_descr_user, &acontext);
      while (adp != NULL && firstGenAnnotSCAD == NULL) {
        uop = (UserObjectPtr) adp->data.ptrvalue;
        if (uop != NULL) {
          oip = uop->type;
          if (oip != NULL) {
            if (StringICmp (oip->str, "StructuredComment") == 0) {
              if (firstGenAnnotSCAD == NULL) {
                for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
                  if (ufp->choice != 1) continue;
                  oip = ufp->label;
                  if (oip == NULL) continue;
                  field = oip->str;
                  if (StringHasNoText (field)) continue;
                  if (StringCmp (field, "StructuredCommentPrefix") == 0) {
                    if (StringCmp ((CharPtr) ufp->data.ptrvalue, "##Genome-Annotation-Data-START##") == 0) {
                      firstGenAnnotSCAD = uop;
                      genomeBuildNumber = NULL;
                      genomeVersionNumber = NULL;
                      firstGenAnnotSCStr = GetStrForStructuredComment (ajp, firstGenAnnotSCAD);
                    }
                  }
                }
              }
            }
          }
        }
        adp = SeqMgrGetNextAnnotDesc (farbsp, adp, Annot_descr_user, &acontext);
      }
      if (firstGenAnnotSCAD == NULL) {
        sdp = SeqMgrGetNextDescriptor (farbsp, NULL, Seq_descr_user, &dcontext);
        while (sdp != NULL) {
          uop = (UserObjectPtr) sdp->data.ptrvalue;
          if (uop != NULL) {
            oip = uop->type;
            if (oip != NULL && StringCmp (oip->str, "StructuredComment") == 0) {
              for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
                if (ufp->choice != 1) continue;
                oip = ufp->label;
                if (oip == NULL) continue;
                field = oip->str;
                if (StringHasNoText (field)) continue;
                if (StringCmp (field, "StructuredCommentPrefix") == 0) {
                  if (StringCmp ((CharPtr) ufp->data.ptrvalue, "##Genome-Annotation-Data-START##") == 0) {
                    firstGenAnnotSCAD = uop;
                    genomeBuildNumber = NULL;
                    genomeVersionNumber = NULL;
                    firstGenAnnotSCStr = GetStrForStructuredComment (ajp, firstGenAnnotSCAD);
                  }
                }
              }
            }
          }
          sdp = SeqMgrGetNextDescriptor (farbsp, sdp, Seq_descr_user, &dcontext);
        }
      }
      BioseqUnlock (farbsp);
      break;
    }
  }

  /*
  also look for Genome-Annotation structured comment descriptor to suppress GenomeBuild user object
  */

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &dcontext);
  while (sdp != NULL) {
    uop = (UserObjectPtr) sdp->data.ptrvalue;
    if (uop != NULL) {
      oip = uop->type;
      if (oip != NULL && StringCmp (oip->str, "StructuredComment") == 0) {
        for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
          if (ufp->choice != 1) continue;
          oip = ufp->label;
          if (oip == NULL) continue;
          field = oip->str;
          if (StringHasNoText (field)) continue;
          if (StringCmp (field, "StructuredCommentPrefix") == 0) {
            if (StringCmp ((CharPtr) ufp->data.ptrvalue, "##Genome-Annotation-Data-START##") == 0) {
              genomeBuildNumber = NULL;
              genomeVersionNumber = NULL;
              if (firstGenAnnotSCAD == NULL) {
                firstGenAnnotSCAD = uop;
                firstGenAnnotSCStr = GetStrForStructuredComment (ajp, firstGenAnnotSCAD);
              }
            }
          }
        }
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &dcontext);
  }

  gi = 0;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    tsip = NULL;
    if (sip->choice == SEQID_OTHER) {
      tsip = (TextSeqIdPtr) sip->data.ptrvalue;

      if (tsip != NULL) {
        is_other = TRUE;
        if (StringNCmp (tsip->accession, "NC_", 3) == 0 || StringNCmp (tsip->accession, "AC_", 3) == 0) {
          if (hasRefTrackStatus) {
            /* will print elsewhere */
          } else if (! StringHasNoText (genomeBuildNumber)) {
            cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
            if (cbp != NULL) {

              cbp->entityID = awp->entityID;
              cbp->first = first;
              cbp->no_blank_before = last_had_tilde;
              first = FALSE;

              if (cbp->first) {
                FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
              } else {
                FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
              }

              FFAddOneString (ffstring, "GENOME ANNOTATION ", FALSE, FALSE, TILDE_IGNORE);

              if ( GetWWW(ajp) ) {
                FFAddOneString (ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
                FF_Add_NCBI_Base_URL (ffstring, ref_link);
                FFAddOneString (ffstring, "\">", FALSE, FALSE, TILDE_IGNORE);
              }
              FFAddOneString (ffstring, "REFSEQ", FALSE, FALSE, TILDE_IGNORE);
              if ( GetWWW(ajp) ) {
                FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
              }
              FFAddOneString (ffstring, ":", FALSE, FALSE, TILDE_IGNORE);

              FFAddOneString (ffstring, reftxt22, FALSE, FALSE, TILDE_EXPAND);
              FFAddOneString (ffstring, genomeBuildNumber, FALSE, FALSE, TILDE_EXPAND);
              if (StringHasNoText (genomeVersionNumber)) {
                genomeVersionNumber = "1";
              }
              FFAddOneString (ffstring, " version ", FALSE, FALSE, TILDE_EXPAND);
              FFAddOneString (ffstring, genomeVersionNumber, FALSE, FALSE, TILDE_EXPAND);
              FFAddOneString (ffstring, reftxt23, FALSE, FALSE, TILDE_EXPAND);

              FFAddOneString (ffstring, " [see ", FALSE, FALSE, TILDE_EXPAND);

              if ( GetWWW(ajp) ) {
                FFAddOneString (ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
                FF_Add_NCBI_Base_URL (ffstring, doc_link);
                FFAddOneString (ffstring, "\">", FALSE, FALSE, TILDE_IGNORE);
              }
              FFAddOneString (ffstring, "documentation", FALSE, FALSE, TILDE_IGNORE);
              if ( GetWWW(ajp) ) {
                FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
              }

              FFAddOneString (ffstring, "].", FALSE, FALSE, TILDE_EXPAND);

              cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
              FFRecycleString(ajp, ffstring);
              ffstring = FFGetString(ajp);

              last_had_tilde = FALSE;
              if (awp->afp != NULL) {
                DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
              }
            }
          }

        } else if (StringNCmp(tsip->accession, "NT_", 3) == 0 || StringNCmp(tsip->accession, "NW_", 3) == 0) {

          if (is_encode) {
            cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
            if (cbp != NULL) {

              cbp->entityID = awp->entityID;
              cbp->first = first;
              cbp->no_blank_before = last_had_tilde;
              first = FALSE;

              if (cbp->first) {
                FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
              } else {
                FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
              }

              FFAddOneString (ffstring, "REFSEQ", FALSE, FALSE, TILDE_IGNORE);
              FFAddOneString (ffstring, ":  ", FALSE, FALSE, TILDE_IGNORE);

              FFAddOneString (ffstring, "This record was provided by the ", FALSE, FALSE, TILDE_EXPAND);
              if ( GetWWW(ajp) ) {
                FFAddOneString (ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
                FF_Add_NCBI_Base_URL (ffstring, link_encode);
                FFAddOneString (ffstring, "\">", FALSE, FALSE, TILDE_IGNORE);
              }
              FFAddOneString (ffstring, "ENCODE", FALSE, FALSE, TILDE_EXPAND);
              if ( GetWWW(ajp) ) {
                FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
              }
              FFAddOneString (ffstring, " project.", FALSE, FALSE, TILDE_EXPAND);

              str = GetEncodeString (encodeUop, bsp);
              if (str != NULL) {
                FFAddOneString (ffstring, "  ", FALSE, FALSE, TILDE_EXPAND);
                FFAddOneString (ffstring, str, FALSE, FALSE, TILDE_EXPAND);
              }
              MemFree (str);

              cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
              FFRecycleString(ajp, ffstring);
              ffstring = FFGetString(ajp);

              last_had_tilde = FALSE;
              if (awp->afp != NULL) {
                DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
              }
            }

          } else if (! hasRefTrackStatus) {

            cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
            if (cbp != NULL) {

              cbp->entityID = awp->entityID;
              cbp->first = first;
              cbp->no_blank_before = last_had_tilde;
              first = FALSE;

              if (cbp->first) {
                FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
              } else {
                FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
              }

              FFAddOneString (ffstring, "GENOME ANNOTATION ", FALSE, FALSE, TILDE_IGNORE);

              if ( GetWWW(ajp) ) {
                FFAddOneString (ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
                FF_Add_NCBI_Base_URL (ffstring, ref_link);
                FFAddOneString (ffstring, "\">", FALSE, FALSE, TILDE_IGNORE);
              }
              FFAddOneString (ffstring, "REFSEQ", FALSE, FALSE, TILDE_IGNORE);
              if ( GetWWW(ajp) ) {
                FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
              }
              FFAddOneString (ffstring, ":", FALSE, FALSE, TILDE_IGNORE);

              if (! StringHasNoText (genomeBuildNumber)) {
                FFAddOneString (ffstring, reftxt22, FALSE, FALSE, TILDE_EXPAND);
                FFAddOneString (ffstring, genomeBuildNumber, FALSE, FALSE, TILDE_EXPAND);
                if (StringHasNoText (genomeVersionNumber)) {
                  genomeVersionNumber = "1";
                }
                FFAddOneString (ffstring, " version ", FALSE, FALSE, TILDE_EXPAND);
                FFAddOneString (ffstring, genomeVersionNumber, FALSE, FALSE, TILDE_EXPAND);
                FFAddOneString (ffstring, reftxt23, FALSE, FALSE, TILDE_EXPAND);

                FFAddOneString (ffstring, " [see ", FALSE, FALSE, TILDE_EXPAND);

                if ( GetWWW(ajp) ) {
                  FFAddOneString (ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
                  FF_Add_NCBI_Base_URL (ffstring, doc_link);
                  FFAddOneString (ffstring, "\">", FALSE, FALSE, TILDE_IGNORE);
                }
                FFAddOneString (ffstring, "documentation", FALSE, FALSE, TILDE_IGNORE);
                if ( GetWWW(ajp) ) {
                  FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
                }

                FFAddOneString (ffstring, "].", FALSE, FALSE, TILDE_EXPAND);
              } else {

                FFAddOneString (ffstring, reftxt21, TRUE, FALSE, TILDE_EXPAND);

                FFAddOneString (ffstring, "~Also see:~    ", FALSE, FALSE, TILDE_EXPAND);

                if ( GetWWW(ajp) ) {
                  FFAddOneString (ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
                  FF_Add_NCBI_Base_URL (ffstring, doc_link);
                  FFAddOneString (ffstring, "\">", FALSE, FALSE, TILDE_IGNORE);
                }
                FFAddOneString (ffstring, "Documentation", FALSE, FALSE, TILDE_IGNORE);
                if ( GetWWW(ajp) ) {
                  FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
                }

                FFAddOneString (ffstring, " of NCBI's Annotation Process~    ", FALSE, FALSE, TILDE_EXPAND);
              }

              cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
              FFRecycleString(ajp, ffstring);
              ffstring = FFGetString(ajp);

              last_had_tilde = TRUE;
              if (awp->afp != NULL) {
                DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
              }
            }
          }

        } else if (StringNCmp(tsip->accession, "XP_", 3) == 0 ||
                   StringNCmp(tsip->accession, "XM_", 3) == 0 ||
                   StringNCmp(tsip->accession, "XR_", 3) == 0 ||
                   StringNCmp(tsip->accession, "ZP_", 3) == 0) {

          name = NULL;
          gi = 0;
          version = 0;
          left = 0;
          right = 0;
          method = NULL;
          mrnaEv = FALSE;
          estEv = FALSE;
          if (GetAnnotationComment (bsp, &name, &assm, &gi, &left, &right, &method, &mrnaEv, &estEv)) {

            cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
            if (cbp != NULL) {

              cbp->entityID = awp->entityID;
              cbp->first = first;
              cbp->no_blank_before = last_had_tilde;
              first = FALSE;

              if (cbp->first) {
                FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
              } else {
                FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
              }

              FFAddOneString (ffstring, "MODEL ", FALSE, FALSE, TILDE_IGNORE);

              if ( GetWWW(ajp) ) {
                FFAddOneString (ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
                FF_Add_NCBI_Base_URL (ffstring, ref_link);
                FFAddOneString (ffstring, "\">", FALSE, FALSE, TILDE_IGNORE);
              }
              FFAddOneString (ffstring, "REFSEQ", FALSE, FALSE, TILDE_IGNORE);
              if ( GetWWW(ajp) ) {
                FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
              }
              FFAddOneString (ffstring, ":  ", FALSE, FALSE, TILDE_IGNORE);

              FFAddTextToString (ffstring, NULL, reftxt11, " (", FALSE, FALSE, TILDE_IGNORE);

              if ( GetWWW(ajp) ) {
                FFAddOneString (ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
                if (IS_protdb_accession (name)) {
                  FF_Add_NCBI_Base_URL (ffstring, link_seqp);
                } else {
                  FF_Add_NCBI_Base_URL (ffstring, link_seqn);
                }
                if (gi > 0) {
                  sprintf (tmp, "%ld", (long) gi);
                  FFAddOneString (ffstring, tmp, FALSE, FALSE, TILDE_IGNORE);
                  msip = GetSeqIdForGI (gi);
                  if (msip != NULL) {
                    switch (msip->choice) {
                      case SEQID_GENBANK:
                      case SEQID_EMBL:
                      case SEQID_DDBJ:
                      case SEQID_OTHER:
                      case SEQID_TPG:
                      case SEQID_TPE:
                      case SEQID_TPD:
                      case SEQID_PIR:
                      case SEQID_SWISSPROT:
                        tsip = (TextSeqIdPtr) msip->data.ptrvalue;
                        if (tsip != NULL) {
                          if (StringICmp (name, tsip->accession) == 0) {
                            version = tsip->version;
                          }
                        }
                        break;
                      default:
                        break;
                    }
                  }
                } else if (ValidateAccnDotVer (name) == 0 && GetGiFromAccnDotVer (name, &gi)) {
                  sprintf (tmp, "%ld", (long) gi);
                  FFAddOneString (ffstring, tmp, FALSE, FALSE, TILDE_IGNORE);
                } else {
                  FFAddOneString (ffstring, name, FALSE, FALSE, TILDE_IGNORE);
                }
                gi = 0;
                FFAddOneString (ffstring, "?report=graph", FALSE, FALSE, TILDE_IGNORE);
                if (left > 0 && right > 0) {
                  if (left > 500) {
                    left -= 500;
                  } else {
                    left = 1;
                  }
                  right += 500;
                  sprintf (tmp, "&v=%ld:%ld", (long) left, (long) right);
                  FFAddOneString (ffstring, tmp, FALSE, FALSE, TILDE_IGNORE);
                }
                FFAddOneString (ffstring, "\">", FALSE, FALSE, TILDE_IGNORE);
                FFAddOneString (ffstring, name, FALSE, FALSE, TILDE_IGNORE);
                if (version > 0 && StringChr (name, '.') == NULL) {
                  sprintf (tmp, ".%ld", (long) version);
                  FFAddOneString (ffstring, tmp, FALSE, FALSE, TILDE_IGNORE);
                }
                FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
              } else if (StringChr (name, '.') == NULL && gi > 0) {
                msip = GetSeqIdForGI (gi);
                if (msip != NULL) {
                  switch (msip->choice) {
                    case SEQID_GENBANK:
                    case SEQID_EMBL:
                    case SEQID_DDBJ:
                    case SEQID_OTHER:
                    case SEQID_TPG:
                    case SEQID_TPE:
                    case SEQID_TPD:
                    case SEQID_PIR:
                    case SEQID_SWISSPROT:
                      tsip = (TextSeqIdPtr) msip->data.ptrvalue;
                      if (tsip != NULL) {
                        if (StringICmp (name, tsip->accession) == 0) {
                          version = tsip->version;
                        }
                      }
                      break;
                    default:
                      break;
                  }
                }
                FFAddOneString (ffstring, name, FALSE, FALSE, TILDE_IGNORE);
                if (version > 0) {
                  sprintf (tmp, ".%ld", (long) version);
                  FFAddOneString (ffstring, tmp, FALSE, FALSE, TILDE_IGNORE);
                }
              } else {
                FFAddOneString (ffstring, name, FALSE, FALSE, TILDE_IGNORE);
              }

              FFAddOneString (ffstring, ")", FALSE, FALSE, TILDE_IGNORE);

              if (assm != NULL) {

                plural = " (";
                count = 0;
                for (tufp = assm->data.ptrvalue; tufp != NULL; tufp = tufp->next)  {
                  ufp = tufp->data.ptrvalue;
                  if (ufp != NULL) {
                    oip = ufp->label;
                    if (oip != NULL && oip->str != NULL && StringICmp (oip->str, "accession") == 0 && ufp->choice == 1) {
                      accn = (CharPtr) ufp->data.ptrvalue;
                      if (StringDoesHaveText (accn)) {
                        count++;
                      }
                    }
                  }
                }
                if (count > 1) {
                  plural = "s (";
                }

                if (count > 0) {
                  FFAddTextToString (ffstring, " ", reftxt13, plural, FALSE, FALSE, TILDE_IGNORE);

                  for (tufp = assm->data.ptrvalue; tufp != NULL; tufp = tufp->next) {
                    accn = NULL;
                    ufp = tufp->data.ptrvalue;
                    if (ufp != NULL) {
                      oip = ufp->label;
                      if (oip != NULL && oip->str != NULL && StringICmp (oip->str, "accession") == 0 && ufp->choice == 1) {
                        accn = (CharPtr) ufp->data.ptrvalue;
                      }
                    }
                    if (StringDoesHaveText (accn)) {
                      if (GetWWW (ajp) && ValidateAccnDotVer (accn) == 0 && GetGiFromAccnDotVer (accn, &gi)) {
                        FFAddOneString (ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
                        if (IS_protdb_accession (nm)) {
                            FF_Add_NCBI_Base_URL (ffstring, link_seqp);
                        } else {
                            FF_Add_NCBI_Base_URL (ffstring, link_seqn);
                        }
                        if (gi > 0) {
                          sprintf (buf, "%ld", (long) gi);
                          FFAddTextToString(ffstring, /* "val=" */ NULL, buf, "\">", FALSE, FALSE, TILDE_IGNORE);
                        } else {
                          FFAddTextToString(ffstring, /* "val=" */ NULL, accn, "\">", FALSE, FALSE, TILDE_IGNORE);
                        }
                        FFAddOneString (ffstring, accn, FALSE, FALSE, TILDE_IGNORE);
                        FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
                      } else if (GetWWW (ajp) && ValidateAccn (accn) == 0) {
                        FFAddOneString (ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
                        if (IS_protdb_accession (nm)) {
                            FF_Add_NCBI_Base_URL (ffstring, link_seqp);
                          } else {
                          FF_Add_NCBI_Base_URL (ffstring, link_seqn);
                        }
                        FFAddTextToString(ffstring, /* "val=" */ NULL, accn, "\">", FALSE, FALSE, TILDE_IGNORE);
                        FFAddOneString (ffstring, accn, FALSE, FALSE, TILDE_IGNORE);
                        FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
                      } else {
                        FFAddOneString (ffstring, accn, FALSE, FALSE, TILDE_IGNORE);
                      }
                    } else if (StringDoesHaveText (nm)) {
                      FFAddOneString (ffstring, nm, FALSE, FALSE, TILDE_IGNORE);
                    } else continue;
                    if (tufp->next != NULL) {
                      ufp = tufp->next;
                      if (ufp->next != NULL) {
                        FFAddOneString (ffstring, ", ", FALSE, FALSE, TILDE_IGNORE);
                      } else {
                        FFAddOneString (ffstring, " and ", FALSE, FALSE, TILDE_IGNORE);
                      }
                    }
                  }

                  FFAddOneString (ffstring, ")", FALSE, FALSE, TILDE_IGNORE);
                }
              }
 
              if (method != NULL) {
                FFAddOneString (ffstring, " ", FALSE, FALSE, TILDE_IGNORE);
                FFAddOneString (ffstring, reftxt12, FALSE, FALSE, TILDE_IGNORE);
                FFAddOneString (ffstring, " ", FALSE, FALSE, TILDE_IGNORE);
                FFAddOneString (ffstring, method, FALSE, FALSE, TILDE_IGNORE);
              }

              if (mrnaEv || estEv) {
                FFAddOneString (ffstring, ", supported by ", FALSE, FALSE, TILDE_IGNORE);
                if (mrnaEv && estEv) {
                  FFAddOneString (ffstring, "mRNA and EST ", FALSE, FALSE, TILDE_IGNORE);
                } else if (mrnaEv) {
                  FFAddOneString (ffstring, "mRNA ", FALSE, FALSE, TILDE_IGNORE);
                } else {
                  FFAddOneString (ffstring, "EST ", FALSE, FALSE, TILDE_IGNORE);
                }
                geneName = NULL;
                locusID [0] = '\0';
                taxID [0] = '\0';
                if ( GetWWW(ajp) && GetGeneAndLocus (bsp, &geneName, locusID, taxID)) {
                  FFAddOneString (ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
                  FF_Add_NCBI_Base_URL (ffstring, ev_link);
                  FFAddTextToString (ffstring, "contig=", name, NULL, FALSE, FALSE, TILDE_IGNORE);
                  FFAddTextToString (ffstring, "&gene=", geneName, NULL, FALSE, FALSE, TILDE_IGNORE);
                  FFAddTextToString (ffstring, "&lid=", locusID, NULL, FALSE, FALSE, TILDE_IGNORE);
                  if (! StringHasNoText (taxID)) {
                    FFAddTextToString (ffstring, "&taxid=", taxID, NULL, FALSE, FALSE, TILDE_IGNORE);
                  }
                  FFAddOneString (ffstring, "\">", FALSE, FALSE, TILDE_IGNORE);
                  FFAddOneString (ffstring, "evidence", FALSE, FALSE, TILDE_IGNORE);
                  FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
                } else {
                  FFAddOneString (ffstring, "evidence", FALSE, FALSE, TILDE_IGNORE);
                }
              }

              FFAddOneString (ffstring, ".", FALSE, FALSE, TILDE_IGNORE);

              FFAddOneString (ffstring, "~Also see:~    ", FALSE, FALSE, TILDE_EXPAND);

              if ( GetWWW(ajp) ) {
                FFAddOneString (ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
                FF_Add_NCBI_Base_URL (ffstring, doc_link);
                FFAddOneString (ffstring, "\">", FALSE, FALSE, TILDE_IGNORE);
              }
              FFAddOneString (ffstring, "Documentation", FALSE, FALSE, TILDE_IGNORE);
              if ( GetWWW(ajp) ) {
                FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
              }

              FFAddOneString (ffstring, " of NCBI's Annotation Process~    ", FALSE, FALSE, TILDE_EXPAND);

              cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
              FFRecycleString(ajp, ffstring);
              ffstring = FFGetString(ajp);

              last_had_tilde = TRUE;
              if (awp->afp != NULL) {
                DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
              }
            }
          }
        } else if (StringNCmp(tsip->accession, "WP_", 3) == 0) {
          cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
          if (cbp != NULL) {
      
            cbp->entityID = awp->entityID;
            cbp->itemID = unverified_itemID;
            cbp->itemtype = OBJ_SEQDESC;
            cbp->first = first;
            cbp->no_blank_before = last_had_tilde;
            first = FALSE;
      
            if (cbp->first) {
              FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
            } else {
              FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
            }
      
            FFAddOneString (ffstring, "REFSEQ:", FALSE, FALSE, TILDE_IGNORE);
            FFAddOneString (ffstring, reftxt51, FALSE, FALSE, TILDE_IGNORE);
            FFAddOneString (ffstring, ".", FALSE, FALSE, TILDE_IGNORE);
      
            cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
            FFRecycleString(ajp, ffstring);
            ffstring = FFGetString(ajp);
      
            last_had_tilde = FALSE;
            if (awp->afp != NULL) {
              DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
            }
          }
        } else if (StringNCmp(tsip->accession, "NZ_", 3) == 0) {
          if (StringLen (tsip->accession) == 15) {
            is_wgs = TRUE;
            if (StringCmp (tsip->accession + 9, "000000") == 0) {
              wgsaccn = tsip->accession;
              wgsname = tsip->name;
            }
          } else if (StringLen (tsip->accession) == 16) {
            is_wgs = TRUE;
            if (StringCmp (tsip->accession + 10, "000000") == 0) {
              wgsaccn = tsip->accession;
              wgsname = tsip->name;
            }
          }
       } else {
          if (StringLen (tsip->accession) == 15) {
            is_wgs = TRUE;
            if (StringCmp (tsip->accession + 9, "000000") == 0) {
              wgsaccn = tsip->accession;
              wgsname = tsip->name; /* master accession has 8 zeroes, name has project version plus 6 zeroes */
            }
          }
        }
      }

    } else if (sip->choice == SEQID_TPG || sip->choice == SEQID_TPE || sip->choice == SEQID_TPD) {

      is_tpa = TRUE;

      tsip = (TextSeqIdPtr) sip->data.ptrvalue;
      if (tsip != NULL && tsip->accession != NULL) {
        acclen = StringLen (tsip->accession);
        tsaaccn = tsip->accession;
        tsaname = tsip->name;
        if (acclen == 12) {
          is_wgs = TRUE;
          if (StringCmp (tsip->accession + 6, "000000") == 0) {
            wgsaccn = tsip->accession;
            wgsname = tsip->name; /* master accession has 8 zeroes, name has project version plus 6 zeroes */
          }
        } else if (acclen == 13) {
          is_wgs = TRUE;
          if (StringCmp (tsip->accession + 6, "0000000") == 0) {
            wgsaccn = tsip->accession;
            wgsname = tsip->name; /* master accession has 9 zeroes, name has project version plus 7 zeroes */
          }
        } else if (acclen == 14) {
          is_wgs = TRUE;
          if (StringCmp (tsip->accession + 6, "00000000") == 0) {
            wgsaccn = tsip->accession;
            wgsname = tsip->name; /* master accession has 10 zeroes, name has project version plus 8 zeroes */
          }
        } else if (ajp->newSourceOrg && StringLen (tsip->accession) == 6) {
          ch = tsip->accession [0];
          if (ch == 'J' || ch == 'K' || ch == 'L' || ch == 'M') {
            showGBBSource = TRUE;
          }
        }
      }

    } else if (sip->choice == SEQID_GENBANK || sip->choice == SEQID_EMBL || sip->choice == SEQID_DDBJ) {

      is_collab = TRUE;

      tsip = (TextSeqIdPtr) sip->data.ptrvalue;
      if (tsip != NULL && tsip->accession != NULL) {
        acclen = StringLen (tsip->accession);
        tsaaccn = tsip->accession;
        tsaname = tsip->name;
        if (acclen == 12) {
          is_wgs = TRUE;
          if (StringCmp (tsip->accession + 6, "000000") == 0) {
            wgsaccn = tsip->accession;
            wgsname = tsip->name; /* master accession has 8 zeroes, name has project version plus 6 zeroes */
          }
        } else if (acclen == 13) {
          is_wgs = TRUE;
          if (StringCmp (tsip->accession + 6, "0000000") == 0) {
            wgsaccn = tsip->accession;
            wgsname = tsip->name; /* master accession has 9 zeroes, name has project version plus 7 zeroes */
          }
        } else if (acclen == 14) {
          is_wgs = TRUE;
          if (StringCmp (tsip->accession + 6, "00000000") == 0) {
            wgsaccn = tsip->accession;
            wgsname = tsip->name; /* master accession has 10 zeroes, name has project version plus 8 zeroes */
          }
        } else if (ajp->newSourceOrg && StringLen (tsip->accession) == 6) {
          ch = tsip->accession [0];
          if (ch == 'J' || ch == 'K' || ch == 'L' || ch == 'M') {
            showGBBSource = TRUE;
          }
        }
      }

    } else if (sip->choice == SEQID_GENERAL) {
      dbt = (DbtagPtr) sip->data.ptrvalue;

      /* show GSDB sequence identifier */

      if (dbt != NULL && StringCmp (dbt->db, "GSDB") == 0 && dbt->tag != NULL) {
        cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
        if (cbp != NULL) {

          cbp->entityID = awp->entityID;
          cbp->first = first;
          cbp->no_blank_before = last_had_tilde;
          first = FALSE;

          /* string will be created after we know if there are additional comments */

          gsdbid = dbt->tag->id;
          sprintf (buf, "GSDB:S:%ld.", (long) gsdbid);

          if (cbp->first) {
            FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
          } else {
            FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
          }

          /* CheckEndPunctuation, ConvertDoubleQuotes, and ExpandTildes already taken into account */

          FFAddOneString (ffstring, buf, FALSE, FALSE, TILDE_IGNORE);

          cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
          FFRecycleString(ajp, ffstring);
          ffstring = FFGetString(ajp);

          last_had_tilde = FALSE;
          if (awp->afp != NULL) {
            DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
          }
        }
      } else if (dbt != NULL && StringCmp (dbt->db, "NCBIFILE") == 0 && dbt->tag != NULL) {
        ncbifileID = dbt->tag;
      }

    } else if (sip->choice == SEQID_GI) {
      gi = (BIG_ID) sip->data.intvalue;

    } else if (sip->choice == SEQID_LOCAL) {
      localID = (ObjectIdPtr) sip->data.ptrvalue;
    }

    if (tsip != NULL) {
      tlstsip = tsip;
    }
  }

  origLocalID = FastaGetOriginalId (bsp);

  if (localID != NULL) {
    if (is_tpa || is_collab) {
      if (awp->mode == SEQUIN_MODE || awp->mode == DUMP_MODE) {
        buf [0] = '\0';
        if (StringDoesHaveText (origLocalID)) {
          if (StringLen (origLocalID) < 1000) {
            sprintf (buf, "LocalID: %s", origLocalID);
            showedLocalID = TRUE;
          } else {
            sprintf (buf, "LocalID string too large");
          }
        } else if (! StringHasNoText (localID->str)) {
          if (StringLen (localID->str) < 1000) {
            sprintf (buf, "LocalID: %s", localID->str);
            showedLocalID = TRUE;
          } else {
            sprintf (buf, "LocalID string too large");
          }
        } else {
          sprintf (buf, "LocalID: %ld", (long) localID->id);
          showedLocalID = TRUE;
        }

        cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
        if (cbp != NULL) {

          cbp->entityID = awp->entityID;
          cbp->first = first;
          cbp->no_blank_before = last_had_tilde;
          first = FALSE;

          if (cbp->first) {
            FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
          } else {
            FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
          }

          FFAddOneString (ffstring, buf, FALSE, FALSE, TILDE_EXPAND);

          cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12,5, 5, "CC");
          FFRecycleString(ajp, ffstring);
          ffstring = FFGetString(ajp);

          last_had_tilde = FALSE;
          if (awp->afp != NULL) {
            DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
          }
        }
      }
    }
  }

  if (ncbifileID != NULL) {
    if (is_tpa || is_collab) {
      if (awp->mode == SEQUIN_MODE || awp->mode == DUMP_MODE) {
        buf [0] = '\0';
        if (! StringHasNoText (ncbifileID->str)) {
          if (StringLen (ncbifileID->str) < 1000) {
            sprintf (buf, "FileID: %s", ncbifileID->str);
          } else {
            sprintf (buf, "FileID string too large");
          }
        } else {
          sprintf (buf, "FileID: %ld", (long) ncbifileID->id);
        }

        cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
        if (cbp != NULL) {

          cbp->entityID = awp->entityID;
          cbp->first = first;
          cbp->no_blank_before = last_had_tilde;
          first = FALSE;

          if (cbp->first) {
            FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
          } else {
            FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
          }

          FFAddOneString (ffstring, buf, FALSE, FALSE, TILDE_EXPAND);

          cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12,5, 5, "CC");
          FFRecycleString(ajp, ffstring);
          ffstring = FFGetString(ajp);

          last_had_tilde = FALSE;
          if (awp->afp != NULL) {
            DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
          }
        }
      }
    }
  }

  /* RefSeq results in allocated comment string */

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &dcontext);
  while (sdp != NULL) {

    uop = (UserObjectPtr) sdp->data.ptrvalue;
    if (uop != NULL) {

      if (! didTPA) {
        str = GetStrForTPA (uop, bsp);
        if (str != NULL) {

          cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
          if (cbp != NULL) {

            cbp->entityID = dcontext.entityID;
            cbp->itemID = dcontext.itemID;
            cbp->itemtype = OBJ_SEQDESC;
            cbp->first = first;
            cbp->no_blank_before = last_had_tilde;
            first = FALSE;

            if (cbp->first) {
              FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
            } else {
              FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
            }

            FFAddOneString (ffstring, str, FALSE, FALSE, TILDE_EXPAND);

            cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12,5, 5, "CC");
            FFRecycleString(ajp, ffstring);
            ffstring = FFGetString(ajp);

            last_had_tilde = FALSE;
            if (awp->afp != NULL) {
              DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
            }
          }
          MemFree (str);
          didTPA = TRUE;
        }
      }

      if (! ajp->flags.hideBankItComment) {
        str = GetStrForBankit (uop, (Boolean) (awp->mode == DUMP_MODE),
                               (Boolean) (showedLocalID && awp->mode == SEQUIN_MODE));
        if (str != NULL) {

          cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
          if (cbp != NULL) {

            cbp->entityID = dcontext.entityID;
            cbp->itemID = dcontext.itemID;
            cbp->itemtype = OBJ_SEQDESC;
            cbp->first = first;
            cbp->no_blank_before = last_had_tilde;
            first = FALSE;

            if (cbp->first) {
              FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
            } else {
              FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
            }

            FFAddOneString (ffstring, str, TRUE, FALSE, TILDE_EXPAND);

            cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12,5, 5, "CC");
            FFRecycleString(ajp, ffstring);
            ffstring = FFGetString(ajp);

            last_had_tilde = FALSE;
            if (awp->afp != NULL) {
              DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
            }
          }
          MemFree (str);
        }
      }

      if (! didRefTrack) {
        str = GetStatusForRefTrack (uop);
        if (str != NULL) {

          cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
          if (cbp != NULL) {

            cbp->entityID = dcontext.entityID;
            cbp->itemID = dcontext.itemID;
            cbp->itemtype = OBJ_SEQDESC;
            cbp->first = first;
            cbp->no_blank_before = last_had_tilde;
            first = FALSE;

            if (cbp->first) {
              FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
            } else {
              FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
            }

            if (StringICmp (str, "Pipeline ") != 0) {
              FFAddOneString (ffstring, str, FALSE, FALSE, TILDE_EXPAND);
            }

            AddStrForRefTrack (ajp, ffstring, uop, ISA_na (bsp->mol), genomeBuildNumber, genomeVersionNumber);

            cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12,5, 5, "CC");
            FFRecycleString(ajp, ffstring);
            ffstring = FFGetString(ajp);

            last_had_tilde = FALSE;
            if (awp->afp != NULL) {
              DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
            }
          }
          /* do not free static str from GetStatusForRefTrack */
          didRefTrack = TRUE;
        }
      }

      if (! didGenome) {
        str = GetStrForGenome (uop, bsp);
        if (str != NULL) {

          cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
          if (cbp != NULL) {

            cbp->entityID = dcontext.entityID;
            cbp->itemID = dcontext.itemID;
            cbp->itemtype = OBJ_SEQDESC;
            cbp->first = first;
            cbp->no_blank_before = last_had_tilde;
            first = FALSE;

            if (cbp->first) {
              FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
            } else {
              FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
            }

            FFAddOneString (ffstring, str, FALSE, FALSE, TILDE_EXPAND);

            cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
            FFRecycleString(ajp, ffstring);
            ffstring = FFGetString(ajp);

            last_had_tilde = FALSE;
            if (awp->afp != NULL) {
              DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
            }
          }
          MemFree (str);
          didGenome = TRUE;
        }
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &dcontext);
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &dcontext);
  while (sdp != NULL) {
    uop = (UserObjectPtr) sdp->data.ptrvalue;
    if (uop != NULL) {
      oip = uop->type;
      if (oip != NULL && StringCmp (oip->str, "RefSeqGene") == 0) {
        for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
          oip = ufp->label;
          if (oip != NULL && StringCmp(oip->str, "Status") == 0 && ufp->choice == 1) {
            str = (CharPtr) ufp->data.ptrvalue;
            if (str != NULL && StringICmp (str, "Reference Standard") == 0) {
              isRefSeqStandard = TRUE;
            }
          }
        }
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &dcontext);
  }
  if (isRefSeqStandard) {
    cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
    if (cbp != NULL) {

      cbp->entityID = awp->entityID;
      cbp->first = first;
      cbp->no_blank_before = last_had_tilde;
      first = FALSE;

      if (cbp->first) {
        FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
      } else {
        FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
      }

        FFAddOneString (ffstring, "This sequence is a reference standard in the ",
                        FALSE, FALSE, TILDE_IGNORE);
      if ( GetWWW(ajp) ) {
        FFAddOneString (ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
        FF_Add_NCBI_Base_URL (ffstring, "https://www.ncbi.nlm.nih.gov/refseq/rsg/");
        FFAddOneString (ffstring, "\">", FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, "RefSeqGene", FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
      } else {
        FFAddOneString (ffstring, "RefSeqGene", FALSE, FALSE, TILDE_IGNORE);
      }
      FFAddOneString (ffstring, " project.", FALSE, FALSE, TILDE_IGNORE);

      cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
      FFRecycleString(ajp, ffstring);
      ffstring = FFGetString(ajp);

      last_had_tilde = FALSE;
      if (awp->afp != NULL) {
        DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
      }
    }
  }

  /*
  if (bsp->repr == Seq_repr_delta && bsp->seq_ext_type == 4 && is_wgs) {
    has_gaps = FALSE;
    for (dsp = (DeltaSeqPtr) bsp->seq_ext; dsp; dsp=dsp->next) {
      if (dsp->choice == 2) {
        litp = (SeqLitPtr) dsp->data.ptrvalue;
        if (litp != NULL) {
          if ((litp->seq_data == NULL || litp->seq_data_type == Seq_code_gap) &&
              litp->length > 0) {
            has_gaps = TRUE;
          }
        }
      }
    }
    if (has_gaps) {
      cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
      if (cbp != NULL) {

        cbp->entityID = awp->entityID;
        cbp->first = first;
        cbp->no_blank_before = last_had_tilde;
        first = FALSE;

        if (cbp->first) {
          FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
        } else {
          FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
        }

        if (is_wgs) {
          FFAddOneString (ffstring, nsWGSGapsString, TRUE, FALSE, TILDE_EXPAND);
        } else {
          FFAddOneString (ffstring, nsAreGapsString, TRUE, FALSE, TILDE_EXPAND);
        }

        cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
        FFRecycleString(ajp, ffstring);
        ffstring = FFGetString(ajp);

        last_had_tilde = FALSE;
        if (awp->afp != NULL) {
          DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
        }
      }
    }
  }
  */

  /* Seq-hist results in allocated comment string */

  hist = bsp->hist;
  if (hist != NULL) {

    if (hist->replaced_by_ids != NULL && hist->replaced_by_date != NULL) {

      okay = TRUE;
      for (sip = hist->replaced_by_ids; sip != NULL; sip = sip->next) {
        if (sip->choice == SEQID_GI) {
          if (gi == (BIG_ID) sip->data.intvalue) {
            okay = FALSE;
          }
        }
      }

      if (okay) {
        cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
        if (cbp != NULL) {

          cbp->entityID = awp->entityID;
          cbp->first = first;
          cbp->no_blank_before = last_had_tilde;
          first = FALSE;

          if (cbp->first) {
            FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
          } else {
            FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
          }

          if (wgsaccn != NULL) {
            AddHistCommentString (ajp, ffstring, "[WARNING] On", "this project was updated. The new version is",
                                  hist->replaced_by_date, hist->replaced_by_ids, ISA_na (bsp->mol), TRUE);
          } else {
            AddHistCommentString (ajp, ffstring, "[WARNING] On", "this sequence was replaced by",
                                  hist->replaced_by_date, hist->replaced_by_ids, ISA_na (bsp->mol), FALSE);
          }

          cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
          FFRecycleString(ajp, ffstring);
          ffstring = FFGetString(ajp);

          last_had_tilde = FALSE;
          if (awp->afp != NULL) {
            DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
          }
        }
      }
    }

    if (hist->replace_ids != NULL && hist->replace_date != NULL && awp->mode != SEQUIN_MODE) {

      okay = TRUE;
      for (sip = hist->replace_ids; sip != NULL; sip = sip->next) {
        if (sip->choice == SEQID_GI) {
          if (gi == (BIG_ID) sip->data.intvalue) {
            okay = FALSE;
          }
        }
      }

      if (okay) {
        cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
        if (cbp != NULL) {

          cbp->entityID = awp->entityID;
          cbp->first = first;
          cbp->no_blank_before = last_had_tilde;
          first = FALSE;

          if (cbp->first) {
            FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
          } else {
            FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
          }

          AddHistCommentString (ajp, ffstring, "On", "this sequence version replaced",
                                hist->replace_date, hist->replace_ids, ISA_na (bsp->mol), FALSE);

          cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
          FFRecycleString(ajp, ffstring);
          ffstring = FFGetString(ajp);

          last_had_tilde = FALSE;
          if (awp->afp != NULL) {
            DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
          }
        }
      }
    }

  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &dcontext);
  while (sdp != NULL) {
    uop = (UserObjectPtr) sdp->data.ptrvalue;
    if (uop != NULL) {
      oip = uop->type;
      if (oip != NULL) {
        if (StringCmp (oip->str, "RefSeqGenome") == 0) {
          cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
          if (cbp != NULL) {

            cbp->entityID = dcontext.entityID;
            cbp->itemID = dcontext.itemID;
            cbp->itemtype = OBJ_SEQDESC;
            cbp->first = first;
            cbp->no_blank_before = last_had_tilde;
            first = FALSE;

            if (cbp->first) {
              FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
            } else {
              FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
            }

            AddStrForRefSeqGenome (ajp, ffstring, uop);

            cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
            FFRecycleString(ajp, ffstring);
            ffstring = FFGetString(ajp);

            last_had_tilde = FALSE;
            if (awp->afp != NULL) {
              DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
            }
          }
        }
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &dcontext);
  }



  /* just save IDs for comment, maploc, and region descriptors */

  /*
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_comment, &dcontext);
  while (sdp != NULL) {
    if (sdp->data.ptrvalue != NULL) {
      cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
      if (cbp != NULL) {
        cbp->entityID = dcontext.entityID;
        cbp->itemID = dcontext.itemID;
        cbp->itemtype = OBJ_SEQDESC;
        cbp->first = first;
        cbp->no_blank_before = last_had_tilde;
        first = FALSE;

        last_had_tilde = FALSE;
        if (awp->afp != NULL) {
          DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
        }
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_comment, &dcontext);
  }
  */

  /* WGS master comment goes before comment descriptors */

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
  if (sdp != NULL) {

    mip = (MolInfoPtr) sdp->data.ptrvalue;
    if (mip != NULL) {
      if (mip->tech == MI_TECH_wgs) {

        if (wgsname != NULL) {

          cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
          if (cbp != NULL) {

            /*
            cbp->entityID = dcontext.entityID;
            cbp->itemID = dcontext.itemID;
            cbp->itemtype = OBJ_SEQDESC;
            */
            cbp->entityID = awp->entityID;
            cbp->first = first;
            cbp->no_blank_before = last_had_tilde;
            first = FALSE;

            if (cbp->first) {
              FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
            } else {
              FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
            }

            AddWGSMasterCommentString (ffstring, bsp, wgsaccn, wgsname);

            cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
            FFRecycleString(ajp, ffstring);
            ffstring = FFGetString(ajp);

            cbp->itemID = dcontext.itemID;
            cbp->itemtype = OBJ_SEQDESC;
            last_had_tilde = FALSE;
            if (awp->afp != NULL) {              
              DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
            }
            cbp->itemID = 0;
            cbp->itemtype = 0;
          }
        }
      } else if (mip->tech == MI_TECH_tsa) {

        if (tsaname != NULL && bsp->repr == Seq_repr_virtual) {

          cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
          if (cbp != NULL) {

            /*
            cbp->entityID = dcontext.entityID;
            cbp->itemID = dcontext.itemID;
            cbp->itemtype = OBJ_SEQDESC;
            */
            cbp->entityID = awp->entityID;
            cbp->first = first;
            cbp->no_blank_before = last_had_tilde;
            first = FALSE;

            if (cbp->first) {
              FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
            } else {
              FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
            }

            AddTSAMasterCommentString (ffstring, bsp, tsaaccn, tsaname);

            cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
            FFRecycleString(ajp, ffstring);
            ffstring = FFGetString(ajp);

            cbp->itemID = dcontext.itemID;
            cbp->itemtype = OBJ_SEQDESC;
            last_had_tilde = FALSE;
            if (awp->afp != NULL) {              
              DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
            }
            cbp->itemID = 0;
            cbp->itemtype = 0;
          }
        }
      } else if (mip->tech == MI_TECH_targeted) {

        if (tlstsip != NULL) {
          tlsaccn = tlstsip->accession;
          tlsname = tlstsip->name;

          if (tlsname != NULL && bsp->repr == Seq_repr_virtual) {
    
            cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
            if (cbp != NULL) {

              /*
              cbp->entityID = dcontext.entityID;
              cbp->itemID = dcontext.itemID;
              cbp->itemtype = OBJ_SEQDESC;
              */
              cbp->entityID = awp->entityID;
              cbp->first = first;
              cbp->no_blank_before = last_had_tilde;
              first = FALSE;

              if (cbp->first) {
                FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
              } else {
                FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
              }

              AddTLSMasterCommentString (ffstring, bsp, tlsaccn, tlsname);

              cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
              FFRecycleString(ajp, ffstring);
              ffstring = FFGetString(ajp);

              cbp->itemID = dcontext.itemID;
              cbp->itemtype = OBJ_SEQDESC;
              last_had_tilde = FALSE;
              if (awp->afp != NULL) {              
                DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
              }
              cbp->itemID = 0;
              cbp->itemtype = 0;
            }
          }
        }
      }
    }
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_genbank, &dcontext);
  if (sdp != NULL) {
    gbp = (GBBlockPtr) sdp->data.ptrvalue;
    if (gbp != NULL) {
      unordered = FALSE;
      for (vnp = gbp->keywords; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringCmp (str, "UNORDERED") == 0) {
          unordered = TRUE;
        }
      }
      if (unordered) {
        cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
        if (cbp != NULL) {

          cbp->entityID = dcontext.entityID;
          cbp->itemID = dcontext.itemID;
          cbp->itemtype = OBJ_SEQDESC;
          cbp->first = first;
          cbp->no_blank_before = last_had_tilde;
          first = FALSE;

          if (cbp->first) {
            FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
          } else {
            FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
          }
          
          AddUnorderedCommentString (ffstring, bsp);

          cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
          FFRecycleString(ajp, ffstring);
          ffstring = FFGetString(ajp);

          last_had_tilde = FALSE;
          if (awp->afp != NULL) {
            DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
          }
        }
      }
    }
  }

  if (showGBBSource) {
    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_genbank, &dcontext);
    if (sdp != NULL) {
      gbp = (GBBlockPtr) sdp->data.ptrvalue;
      if (gbp != NULL && (! StringHasNoText (gbp->source))) {
        cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
        if (cbp != NULL) {

          cbp->entityID = dcontext.entityID;
          cbp->itemID = dcontext.itemID;
          cbp->itemtype = OBJ_SEQDESC;
          cbp->first = first;
          cbp->no_blank_before = last_had_tilde;
          first = FALSE;

          if (cbp->first) {
            FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
          } else {
            FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
          }

          FFAddOneString (ffstring, "Original source text: ", FALSE, FALSE, TILDE_EXPAND);
          FFAddOneString (ffstring, gbp->source, TRUE, TRUE, TILDE_EXPAND);

          cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
          FFRecycleString(ajp, ffstring);
          ffstring = FFGetString(ajp);

          last_had_tilde = FALSE;
          if (awp->afp != NULL) {
            DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
          }
        }
      }
    }
  }

  last_name = NULL;
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_comment, &dcontext);
  while (sdp != NULL) {
    str = (CharPtr) sdp->data.ptrvalue;
    if (StringDoesHaveText (str) && (last_name == NULL || CommentsAreDifferent (str, last_name) || awp->mode == DUMP_MODE)) {
      cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
      if (cbp != NULL) {

        last_name = (CharPtr) str;

        cbp->entityID = dcontext.entityID;
        cbp->itemID = dcontext.itemID;
        cbp->itemtype = OBJ_SEQDESC;
        cbp->first = first;
        cbp->no_blank_before = last_had_tilde;
        first = FALSE;

        last_had_tilde = FALSE;
        len = StringLen (str);
        if (len > 4 && str [len - 1] == '~' && str [len - 2] == '~') {
          last_had_tilde = TRUE;
        }
        if (awp->afp != NULL) {
          DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
        }
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_comment, &dcontext);
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_maploc, &dcontext);
  while (sdp != NULL) {
    cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
    if (cbp != NULL) {

      cbp->entityID = dcontext.entityID;
      cbp->itemID = dcontext.itemID;
      cbp->itemtype = OBJ_SEQDESC;
      cbp->first = first;
      cbp->no_blank_before = last_had_tilde;
      first = FALSE;

      last_had_tilde = FALSE;
      if (awp->afp != NULL) {
        DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_maploc, &dcontext);
  }

  last_name = NULL;
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_region, &dcontext);
  while (sdp != NULL) {
    str = (CharPtr) sdp->data.ptrvalue;
    if (StringDoesHaveText (str) &&
        ((last_name == NULL || StringCmp (str, last_name) != 0) || awp->mode == DUMP_MODE) &&
        (StringCmp (str, ".") != 0 || awp->mode == DUMP_MODE)) {
      cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
      if (cbp != NULL) {

        last_name = (CharPtr) str;

        cbp->entityID = dcontext.entityID;
        cbp->itemID = dcontext.itemID;
        cbp->itemtype = OBJ_SEQDESC;
        cbp->first = first;
        cbp->no_blank_before = last_had_tilde;
        first = FALSE;

        last_had_tilde = FALSE;
        if (awp->afp != NULL) {
          DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
        }
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_region, &dcontext);
  }

  last_name = NULL;
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_name, &dcontext);
  while (sdp != NULL) {
    str = (CharPtr) sdp->data.ptrvalue;
    if (StringDoesHaveText (str) &&
        ((last_name == NULL || StringCmp (str, last_name) != 0) || awp->mode == DUMP_MODE) &&
        (StringCmp (str, ".") != 0 || awp->mode == DUMP_MODE)) {
      cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
      if (cbp != NULL) {

        last_name = (CharPtr) str;

        cbp->entityID = dcontext.entityID;
        cbp->itemID = dcontext.itemID;
        cbp->itemtype = OBJ_SEQDESC;
        cbp->first = first;
        cbp->no_blank_before = last_had_tilde;
        first = FALSE;

        last_had_tilde = FALSE;
        if (awp->afp != NULL) {
          DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
        }
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_name, &dcontext);
  }

  if (basemodNum > 0 && (basemodURLhead != NULL || basemodURL != NULL)) {
    cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
    if (cbp != NULL) {

      cbp->entityID = awp->entityID;
      cbp->itemID = filetrack_itemID;
      cbp->itemtype = OBJ_SEQDESC;
      cbp->first = first;
      cbp->no_blank_before = last_had_tilde;
      first = FALSE;

      if (cbp->first) {
        FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
      } else {
        FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
      }

      if (! last_had_tilde && ! cbp->first) {
        FFAddOneString (ffstring, "\n", FALSE, FALSE, TILDE_EXPAND);
      }

      if (basemodNum == 1) {
        FFAddOneString (ffstring, "This genome has a ", FALSE, FALSE, TILDE_IGNORE);
        if (GetWWW (ajp)) {
          str = NULL;
          if (basemodURL != NULL) {
            str = basemodURL;
          } else if (basemodURLhead != NULL) {
            str = basemodURLhead [0];
          }
          if (StringDoesHaveText (str)) {
            FFAddOneString(ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
            FFAddOneString (ffstring, str, FALSE, FALSE, TILDE_IGNORE);
            FFAddOneString(ffstring, "\">", FALSE, FALSE, TILDE_IGNORE);
            FFAddOneString (ffstring, "base modification file", FALSE, FALSE, TILDE_IGNORE);
            FFAddOneString(ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
          }
        } else {
          FFAddOneString (ffstring, "base modification file", FALSE, FALSE, TILDE_IGNORE);
        }
        FFAddOneString (ffstring, " available.", FALSE, FALSE, TILDE_IGNORE);
      } else {
        FFAddOneString (ffstring, "There are ", FALSE, FALSE, TILDE_IGNORE);
        sprintf (buf, "%ld", (long) basemodNum);
        FFAddOneString (ffstring, buf, FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, " base modification files", FALSE, FALSE, TILDE_IGNORE);
        if (GetWWW (ajp)) {
          pfx = " (";
          sfx = "";
          for (j = 0; j < basemodNum; j++) {
            str = NULL;
            if (basemodURL != NULL) {
              str = basemodURL;
            } else if (basemodURLhead != NULL) {
              str = basemodURLhead [j];
            }
            if (StringHasNoText (str)) continue;
            FFAddOneString (ffstring, pfx, FALSE, FALSE, TILDE_IGNORE);
            FFAddOneString(ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
            FFAddOneString (ffstring, str, FALSE, FALSE, TILDE_IGNORE);
            FFAddOneString(ffstring, "\">", FALSE, FALSE, TILDE_IGNORE);
            sprintf (buf, "%ld", (long) (j + 1));
            FFAddOneString (ffstring, buf, FALSE, FALSE, TILDE_IGNORE);
            FFAddOneString(ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
            if (basemodNum == 2) {
              pfx = " and ";
            } else if (j == basemodNum - 2) {
              pfx = ", and ";
            } else {
              pfx = ", ";
            }
            sfx = ")";
          }
          FFAddOneString (ffstring, sfx, FALSE, FALSE, TILDE_IGNORE);
        }
        FFAddOneString (ffstring, " available for this genome.", FALSE, FALSE, TILDE_IGNORE);
      }

      cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
      FFRecycleString(ajp, ffstring);
      ffstring = FFGetString(ajp);

      last_had_tilde = FALSE;
      if (awp->afp != NULL) {
        DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
      }
    }
  }

  /* StructuredComment user object */

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &dcontext);
  while (sdp != NULL) {
    uop = (UserObjectPtr) sdp->data.ptrvalue;
    if (uop != NULL) {
      oip = uop->type;
      if (oip != NULL && StringCmp (oip->str, "StructuredComment") == 0) {
        for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
          if (ufp->choice != 1) continue;
          oip = ufp->label;
          if (oip == NULL) continue;
          field = oip->str;
          if (StringHasNoText (field)) continue;
          if (StringCmp (field, "StructuredCommentPrefix") == 0) {
            if (StringCmp ((CharPtr) ufp->data.ptrvalue, "##Genome-Annotation-Data-START##") == 0) {
              if (firstGenAnnotSCAD == NULL) {
                firstGenAnnotSCAD = uop;
                genomeBuildNumber = NULL;
                genomeVersionNumber = NULL;
                firstGenAnnotSCStr = GetStrForStructuredComment (ajp, firstGenAnnotSCAD);
                uop = NULL;
              } else {
                firstGenAnnotSCAD = NULL;
              }
              break;
            }
          }
        }
        if (uop != NULL) {
          cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
          if (cbp != NULL) {

            cbp->entityID = dcontext.entityID;
            cbp->itemID = dcontext.itemID;
            cbp->itemtype = OBJ_SEQDESC;
            cbp->first = first;
            cbp->no_blank_before = last_had_tilde;
            first = FALSE;

            last_had_tilde = FALSE;
            if (awp->afp != NULL) {
              DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
            }
          }
        }
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_user, &dcontext);
  }

  /* HTGS results in allocated comment string */

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
  if (sdp != NULL) {

    mip = (MolInfoPtr) sdp->data.ptrvalue;
    if (mip != NULL) {
      if (mip->completeness != 0 && is_other) {

        str = GetMolInfoCommentString (bsp, mip);

        if (str != NULL) {
          cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
          if (cbp != NULL) {

            cbp->entityID = dcontext.entityID;
            cbp->itemID = dcontext.itemID;
            cbp->itemtype = OBJ_SEQDESC;
            cbp->first = first;
            cbp->no_blank_before = last_had_tilde;
            first = FALSE;

            if (cbp->first) {
              FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
            } else {
              FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
            }

            FFAddOneString (ffstring, str, TRUE, FALSE, TILDE_EXPAND);

            cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
            FFRecycleString(ajp, ffstring);
            ffstring = FFGetString(ajp);

            last_had_tilde = FALSE;
            if (awp->afp != NULL) {
              DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
            }
          }
        }

      }
      if (mip->tech == MI_TECH_htgs_0 ||
          mip->tech == MI_TECH_htgs_1 ||
          mip->tech == MI_TECH_htgs_2) {

        cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
        if (cbp != NULL) {

          /*
          cbp->entityID = dcontext.entityID;
          cbp->itemID = dcontext.itemID;
          cbp->itemtype = OBJ_SEQDESC;
          */
          cbp->entityID = awp->entityID;
          cbp->first = first;
          cbp->no_blank_before = last_had_tilde;
          first = FALSE;

          if (cbp->first) {
            FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
          } else {
            FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
          }
          
          AddHTGSCommentString (ffstring, bsp, mip);

          cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
          FFRecycleString(ajp, ffstring);
          ffstring = FFGetString(ajp);

          last_had_tilde = FALSE;
          if (awp->afp != NULL) {
            DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
          }
        }

      } else {
        str = StringForSeqTech (mip->tech);
        if (! StringHasNoText (str)) {

          cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
          if (cbp != NULL) {

            /*
            cbp->entityID = dcontext.entityID;
            cbp->itemID = dcontext.itemID;
            cbp->itemtype = OBJ_SEQDESC;
            */
            cbp->entityID = awp->entityID;
            cbp->first = first;
            cbp->no_blank_before = last_had_tilde;
            first = FALSE;

            if (cbp->first) {
              FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
            } else {
              FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
            }

            FFAddTextToString (ffstring, "Method: ", str, NULL, TRUE, FALSE, TILDE_EXPAND);

            cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
            FFRecycleString(ajp, ffstring);
            ffstring = FFGetString(ajp);

            last_had_tilde = FALSE;
            if (awp->afp != NULL) {
              DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
            }
          }
        }
      }
    }
  }

  /* no longer adding comment features that are full length on appropriate segment */

  /*
  parent = awp->parent;
  if (parent == NULL) return;

  sfp = SeqMgrGetNextFeature (parent, NULL, SEQFEAT_COMMENT, 0, &fcontext);
  while (sfp != NULL) {
    if (fcontext.left == awp->from && fcontext.right == awp->to) {
      cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
      if (cbp != NULL) {

        cbp->entityID = fcontext.entityID;
        cbp->itemID = fcontext.itemID;
        cbp->itemtype = OBJ_SEQFEAT;
        cbp->first = first;
        cbp->no_blank_before = last_had_tilde;
        first = FALSE;

        last_had_tilde = FALSE;
        if (awp->afp != NULL) {
          DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
        }
      }
    }
    sfp = SeqMgrGetNextFeature (parent, sfp, SEQFEAT_COMMENT, 0, &fcontext);
  }
  */

  /*
  search for Seq-annot.desc.comment on annots packaged on current bioseq
  is now done earlier in order to suppress GenomeBuild user object comment
  */

  /*
  annotDescCommentToComment = FALSE;
  adp = SeqMgrGetNextAnnotDesc (bsp, NULL, Annot_descr_user, &acontext);
  while (adp != NULL) {
    uop = (UserObjectPtr) adp->data.ptrvalue;
    if (uop != NULL) {
      oip = uop->type;
      if (oip != NULL) {
        if (StringCmp (oip->str, "AnnotDescCommentPolicy") == 0) {
          for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
            oip = ufp->label;
            if (oip == NULL || ufp->data.ptrvalue == NULL) continue;
            if (StringCmp (oip->str, "Policy") == 0) {
              if (StringICmp ((CharPtr) ufp->data.ptrvalue, "ShowInComment") == 0) {
                annotDescCommentToComment = TRUE;
              }
            }
          }
        } else if (StringICmp (oip->str, "StructuredComment") == 0) {
          if (firstGenAnnotSCAD == NULL) {
            for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
              if (ufp->choice != 1) continue;
              oip = ufp->label;
              if (oip == NULL) continue;
              field = oip->str;
              if (StringHasNoText (field)) continue;
              if (StringCmp (field, "StructuredCommentPrefix") == 0) {
                if (StringCmp ((CharPtr) ufp->data.ptrvalue, "##Genome-Annotation-Data-START##") == 0) {
                  firstGenAnnotSCAD = uop;
                }
              }
            }
          }
        }
      }
    }
    adp = SeqMgrGetNextAnnotDesc (bsp, adp, Annot_descr_user, &acontext);
  }
  */

  if (annotDescCommentToComment) {
    adp = SeqMgrGetNextAnnotDesc (bsp, NULL, Annot_descr_comment, &acontext);
    while (adp != NULL) {
      str = (CharPtr) adp->data.ptrvalue;
      if (StringDoesHaveText (str)) {
        cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
        if (cbp != NULL) {

          cbp->entityID = awp->entityID;
          cbp->first = first;
          cbp->no_blank_before = last_had_tilde;
          first = FALSE;

          if (cbp->first) {
            FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
          } else {
            FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
          }

          FFAddOneString (ffstring, str, TRUE, FALSE, TILDE_EXPAND);

          cbp->string = FFEndPrint (ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
          FFRecycleString (ajp, ffstring);
          ffstring = FFGetString (ajp);

          last_had_tilde = FALSE;
          if (awp->afp != NULL) {
            DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
          }
        }
      }
      adp = SeqMgrGetNextAnnotDesc (bsp, adp, Annot_descr_comment, &acontext);
    }
  }

  if (firstGenAnnotSCAD != NULL) {
    if (StringDoesHaveText (firstGenAnnotSCStr)) {
      cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
      if (cbp != NULL) {

        cbp->entityID = awp->entityID;
        cbp->first = first;
        cbp->no_blank_before = last_had_tilde;

        if (cbp->first) {
          FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
        } else {
          FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
          if (! last_had_tilde) {
            FFAddOneString (ffstring, "\n", FALSE, FALSE, TILDE_EXPAND);
          }
        }

        first = FALSE;

        FFAddOneString (ffstring, firstGenAnnotSCStr, FALSE, FALSE, TILDE_EXPAND);

        cbp->string = FFEndPrint (ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
        FFRecycleString (ajp, ffstring);
        ffstring = FFGetString (ajp);

        if (awp->afp != NULL) {
          DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
        }
      }
    }
  }
  if (firstGenAnnotSCStr != NULL) {
    MemFree (firstGenAnnotSCStr);
  }

  num = 0;
  if (filetrackspp != NULL) {
    num = 1;
  } else if (filetrackpsp != NULL) {
    num = PackSeqPntNum (filetrackpsp);
  }
  if (num > 0) {
    cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
    if (cbp != NULL) {

      cbp->entityID = awp->entityID;
      cbp->itemID = filetrack_itemID;
      cbp->itemtype = OBJ_SEQDESC;
      cbp->first = first;
      cbp->no_blank_before = last_had_tilde;
      first = FALSE;

      if (cbp->first) {
        FFStartPrint (ffstring, awp->format, 0, 12, "COMMENT", 12, 5, 5, "CC", TRUE);
      } else {
        FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
      }

      FFAddOneString (ffstring, "This ", FALSE, FALSE, TILDE_IGNORE);
      if (GetWWW (ajp) && filetrackURL != NULL) {
        FFAddOneString(ffstring, "<a href=\"", FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString (ffstring, filetrackURL, FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString(ffstring, "\">", FALSE, FALSE, TILDE_IGNORE); 
        FFAddOneString (ffstring, "map", FALSE, FALSE, TILDE_IGNORE);
        FFAddOneString(ffstring, "</a>", FALSE, FALSE, TILDE_IGNORE);
      } else {
        FFAddOneString (ffstring, "map", FALSE, FALSE, TILDE_IGNORE);
      }
      FFAddOneString (ffstring, " has ", FALSE, FALSE, TILDE_IGNORE);
      frags = num;

      if (bsp->topology != TOPOLOGY_CIRCULAR) {
        if (num > 1 && GetFileTrackPoint (filetrackspp, filetrackpsp, num - 1) < bsp->length - 1 ) {
          frags = num + 1;
        }
      }

      sprintf (tmp, "%ld", (long) frags);
      FFAddOneString (ffstring, tmp, FALSE, FALSE, TILDE_IGNORE);
      if (frags > 1) {
        FFAddOneString (ffstring, " pieces:", FALSE, FALSE, TILDE_IGNORE);
      } else if (frags == 1) {
        FFAddOneString (ffstring, " piece:", FALSE, FALSE, TILDE_IGNORE);
      }

      last = 1;
      pos = GetFileTrackPoint (filetrackspp, filetrackpsp, 0) + 1;
      if (bsp->topology != TOPOLOGY_CIRCULAR) {

        FFAddNewLine (ffstring);
        sprintf (tmp, "*  %7ld %7ld: fragment of %ld bp in length",
                 (long) last, (long) pos, (long) (pos - last + 1));
        FFAddOneString (ffstring, tmp, FALSE, FALSE, TILDE_IGNORE);

      }
      last = pos + 1;

      chunk = 0;
      for (idx = 1; idx < num; idx++) {

        chunk++;
        if (chunk >= 100) {
          chunk = 0;

          cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
          FFRecycleString(ajp, ffstring);
          ffstring = FFGetString(ajp);

          last_had_tilde = FALSE;
          if (awp->afp != NULL) {
            DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
          }

          cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
          if (cbp != NULL) {
            cbp->entityID = awp->entityID;
            cbp->itemID = filetrack_itemID;
            cbp->itemtype = OBJ_SEQDESC;
            cbp->first = FALSE;
            FFStartPrint (ffstring, awp->format, 0, 12, NULL, 12, 5, 5, "CC", FALSE);
          }
        } else {
          FFAddNewLine (ffstring);
        }

        pos = GetFileTrackPoint (filetrackspp, filetrackpsp, idx) + 1;

        sprintf (tmp, "*  %7ld %7ld: fragment of %ld bp in length",
                 (long) last, (long) pos, (long) (pos - last + 1));
        FFAddOneString (ffstring, tmp, FALSE, FALSE, TILDE_IGNORE);

        last = pos + 1;
      }

      if (bsp->topology != TOPOLOGY_CIRCULAR) {
        pos = bsp->length;

        if (last < pos) {
          FFAddNewLine (ffstring);
          sprintf (tmp, "*  %7ld %7ld: fragment of %ld bp in length",
                   (long) last, (long) pos, (long) (pos - last + 1));
          FFAddOneString (ffstring, tmp, FALSE, FALSE, TILDE_IGNORE);
        }

      } else {
        pos = GetFileTrackPoint (filetrackspp, filetrackpsp, 0) + 1;

        FFAddNewLine (ffstring);
        sprintf (tmp, "*  %7ld %7ld: fragment of %ld bp in length",
                 (long) last, (long) pos, (long) (bsp->length + pos - last + 1));
        FFAddOneString (ffstring, tmp, FALSE, FALSE, TILDE_IGNORE);
      }

      cbp->string = FFEndPrint(ajp, ffstring, awp->format, 12, 12, 5, 5, "CC");
      FFRecycleString(ajp, ffstring);
      ffstring = FFGetString(ajp);

      last_had_tilde = FALSE;
      if (awp->afp != NULL) {
        DoImmediateFormat (awp->afp, (BaseBlockPtr) cbp);
      }
    }
  }

  FFRecycleString(ajp, ffstring);
}

NLM_EXTERN void AddFeatHeaderBlock (
  Asn2gbWorkPtr awp
)

{
  IntAsn2gbJobPtr ajp;
  BaseBlockPtr    bbp;
  Char            buf [128];
  StringItemPtr   ffstring;
  CharPtr         suffix = NULL;

  if (awp == NULL) return;
  ajp = awp->ajp;
  if (ajp == NULL) return;

  bbp = Asn2gbAddBlock (awp, FEATHEADER_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  bbp->entityID = awp->entityID;

  if (GetWWW (ajp) && awp->mode == ENTREZ_MODE && awp->afp != NULL &&
      (awp->format == GENBANK_FMT || awp->format == GENPEPT_FMT)) {
    sprintf (buf, "<a name=\"feature_%s\"></a>", awp->currAccVerLabel);
    DoQuickLinkFormat (awp->afp, buf);
  }

  if (awp->format != FTABLE_FMT) {
    ffstring = FFGetString(ajp);
    if ( ffstring == NULL ) return;

    FFStartPrint (ffstring, awp->format, 0, 12, "FEATURES", 21, 5, 0, "FH", TRUE);

    if (awp->format == EMBL_FMT || awp->format == EMBLPEPT_FMT) {
      FFAddOneString (ffstring, "Key", FALSE, FALSE, TILDE_IGNORE);
      FFAddNChar(ffstring, ' ', 13 , FALSE);
    }

    FFAddOneString (ffstring, "Location/Qualifiers", FALSE, FALSE, TILDE_TO_SPACES);

    if (awp->format == EMBL_FMT || awp->format == EMBLPEPT_FMT) {
      FFAddNewLine(ffstring);
      FFAddNewLine(ffstring);
    }

    suffix = FFEndPrint(ajp, ffstring, awp->format, 12, 21, 5, 0, "FH");
    FFRecycleString(ajp, ffstring);
  }

  bbp->string = suffix;

  if (awp->afp != NULL) {
    DoImmediateFormat (awp->afp, bbp);
  }
}

static Uint2 ComputeSourceHash (
  CharPtr key,
  Uint2 start
)

{
  Uint4  h;
  Uint2  M;
  Uint2  S;

  if (key == NULL) return start;

  M = 101; /* prime key */
  S = 256; /* size of alphabet */

  for (h = start; *key != '\0'; key++) {
    h = (S * h + *key) % M;
  }

  return (Uint2) h;
}

static BaseBlockPtr AddSource (
  Asn2gbWorkPtr awp,
  ValNodePtr PNTR head,
  BioSourcePtr biop,
  CharPtr comment
)

{
  BaseBlockPtr    bbp;
  DbtagPtr        dbt;
  Uint2           hash;
  SourceType      idx;
  IntSrcBlockPtr  isp;
  ObjectIdPtr     oip;
  OrgModPtr       omp;
  OrgNamePtr      onp;
  OrgRefPtr       orp;
  SubSourcePtr    ssp;
  CharPtr         str;
  Uint1           subtype;
  Char            tmp [16];
  ValNodePtr      vnp;

  if (awp == NULL || head == NULL || biop == NULL) return NULL;

  bbp = (BaseBlockPtr) MemNew (sizeof (IntSrcBlock));
  if (bbp == NULL) return NULL;
  bbp->blocktype = SOURCEFEAT_BLOCK;
  bbp->section = awp->currsection;

  ValNodeAddPointer (head, 0, bbp);

  isp = (IntSrcBlockPtr) bbp;
  isp->biop = biop;
  isp->is_focus = biop->is_focus;
  if (biop->origin == 5) {
    isp->is_synthetic = TRUE;
  }

  orp = biop->org;
  if (orp == NULL) return bbp;

  if (StringICmp (orp->taxname, "synthetic construct") == 0) {
    isp->is_synthetic = TRUE;
  }

  isp->orghash = ComputeSourceHash (orp->taxname, 0);
  isp->taxname = orp->taxname;

  hash = 0;
  onp = orp->orgname;
  if (onp != NULL) {
    if (StringICmp (onp->div, "SYN") == 0) {
      isp->is_synthetic = TRUE;
    }
    isp->omp = onp->mod;
    for (omp = onp->mod; omp != NULL; omp = omp->next) {
      subtype = omp->subtype;
      if (subtype == 253) {
        subtype = 39;
      } else if (subtype == 254) {
        subtype = 40;
      } else if (subtype == 255) {
        subtype = 41;
      }
      if (subtype < 42) {
        idx = orgModToSourceIdx [subtype];
        if (idx > 0 && idx < ASN2GNBK_TOTAL_SOURCE) {
          str = asn2gnbk_source_quals [idx].name;
          hash = ComputeSourceHash (str, hash);
          hash = ComputeSourceHash (omp->subname, hash);
        }
      }
    }
  }
  if (comment != NULL) {
    hash = ComputeSourceHash ("note", hash);
    hash = ComputeSourceHash (comment, hash);
  }
  isp->modhash = hash;

  hash = 0;
  for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
    subtype = ssp->subtype;
    if (subtype == 255) {
      subtype = 44;
    }
    if (subtype < 45) {
      idx = subSourceToSourceIdx [subtype];
      if (idx > 0 && idx < ASN2GNBK_TOTAL_SOURCE) {
        str = asn2gnbk_source_quals [idx].name;
        hash = ComputeSourceHash (str, hash);
        hash = ComputeSourceHash (ssp->name, hash);
      }
    }
  }
  isp->subhash = hash;
  isp->ssp = biop->subtype;

  hash = 0;
  for (vnp = orp->db; vnp != NULL; vnp = vnp->next) {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt != NULL) {
      hash = ComputeSourceHash (dbt->db, hash);
      oip = dbt->tag;
      if (oip != NULL) {
        if (oip->str != NULL) {
          hash = ComputeSourceHash (oip->str, hash);
        } else {
          sprintf (tmp, "%ld", (long) oip->id);
          hash = ComputeSourceHash (tmp, hash);
        }
      }
    }
  }
  isp->xrfhash = hash;
  isp->vnp = orp->db;

  return bbp;
}

static int LIBCALLBACK SortSourcesByHash (
  VoidPtr ptr1,
  VoidPtr ptr2
)

{
  Int4            diff;
  IntSrcBlockPtr  isp1;
  IntSrcBlockPtr  isp2;
  ValNodePtr      vnp1;
  ValNodePtr      vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  isp1 = (IntSrcBlockPtr) vnp1->data.ptrvalue;
  isp2 = (IntSrcBlockPtr) vnp2->data.ptrvalue;
  if (isp1 == NULL || isp2 == NULL) return 0;

  if (isp1->is_focus && (! isp2->is_focus)) return -1;
  if (isp2->is_focus && (! isp1->is_focus)) return 1;

  diff = isp1->orghash - isp2->orghash;
  if (diff > 0) return -1;
  if (diff < 0) return 1;

  diff = isp1->xrfhash - isp2->xrfhash;
  if (diff > 0) return -1;
  if (diff < 0) return 1;

  /* sort so that sources with modifiers come first */

  diff = isp1->modhash - isp2->modhash;
  if (diff > 0) return -1;
  if (diff < 0) return 1;

  diff = isp1->subhash - isp2->subhash;
  if (diff > 0) return -1;
  if (diff < 0) return 1;

  /* if all hashes are equal, descriptor comes first */

  if (isp1->is_descriptor && (! isp2->is_descriptor)) {
    return -1;
  } else if (isp2->is_descriptor && (! isp1->is_descriptor)) {
    return 1;
  }

  /* now sort identical sources by position, to only fuse abutting ones */
  /* feature with smallest left extreme is first */

  if (isp1->left > isp2->left) {
    return 1;
  } else if (isp1->left < isp2->left) {
    return -1;
  }

  /* if same left extreme, shortest source feature is first just for flatfile */

  if (isp1->right > isp2->right) {
    return 1;
  } else if (isp1->right < isp2->right) {
    return -1;
  }

  return 0;
}

static int LIBCALLBACK SortSourcesByPos (
  VoidPtr ptr1,
  VoidPtr ptr2
)

{
  IntSrcBlockPtr  isp1;
  IntSrcBlockPtr  isp2;
  ValNodePtr      vnp1;
  ValNodePtr      vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  isp1 = (IntSrcBlockPtr) vnp1->data.ptrvalue;
  isp2 = (IntSrcBlockPtr) vnp2->data.ptrvalue;
  if (isp1 == NULL || isp2 == NULL) return 0;

  /* descriptor always goes first */

  if (isp1->is_descriptor && (! isp2->is_descriptor)) {
    return -1;
  } else if (isp2->is_descriptor && (! isp1->is_descriptor)) {
    return 1;
  }

  /* feature with smallest left extreme is first */

  if (isp1->left > isp2->left) {
    return 1;
  } else if (isp1->left < isp2->left) {
    return -1;
  }

  /* if same left extreme, shortest source feature is first just for flatfile */

  if (isp1->right > isp2->right) {
    return 1;
  } else if (isp1->right < isp2->right) {
    return -1;
  }

  return 0;
}

/*                                                                   */
/* s_isFuzzyLoc () -- Determines is a location has fuzzy coordinates */
/*                                                                   */

static Boolean s_isFuzzyLoc ( SeqLocPtr pLocation )
{
  SeqIntPtr pIntLocation;

  if (pLocation == NULL)
    return FALSE;

  if (pLocation->choice != SEQLOC_INT)
    return FALSE;

  if (pLocation->data.ptrvalue == NULL)
    return FALSE;

  pIntLocation = (SeqIntPtr) pLocation->data.ptrvalue;

  if ((pIntLocation->if_from != NULL) && (pIntLocation->if_from->choice == 2))
    return TRUE;

  if ((pIntLocation->if_to != NULL) && (pIntLocation->if_to->choice == 2))
    return TRUE;

  return FALSE;
}

static void GetSourcesOnBioseq (
  Asn2gbWorkPtr awp,
  BioseqPtr target,
  BioseqPtr bsp,
  Int4 from,
  Int4 to,
  SeqFeatPtr cds
)

{
  IntAsn2gbJobPtr    ajp;
  BaseBlockPtr       bbp;
  BioSourcePtr       biop;
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext;
  Boolean            hasNulls;
  Int4               left;
  Boolean            loop = FALSE;
  Int2               idx;
  IntSrcBlockPtr     isp;
  Boolean            is_wp = FALSE;
  Int4Ptr            ivals;
  SeqLocPtr          newloc;
  Boolean            noLeft;
  Boolean            noRight;
  Int2               numivals;
  Int2               num_super_kingdom = 0;
  Boolean            okay;
  OrgNamePtr         onp;
  OrgRefPtr          orp;
  ObjValNodePtr      ovp;
  Int4               right;
  SeqDescrPtr        sdp;
  ValNodePtr         sdplist = NULL;
  SeqFeatPtr         sfp;
  SeqInt             sint;
  SeqIntPtr          sintp;
  SeqIdPtr           sip;
  SeqLocPtr          slp, slpx;
  Boolean            split;
  SeqPntPtr          spp;
  Int4               start;
  Int4               stop;
  Uint1              strand;
  Boolean            super_kingdoms_different = FALSE;
  CharPtr            super_kingdom_name = NULL;
  TaxElementPtr      tep;
  TextSeqIdPtr       tsip;
  ValNode            vn;
  ValNodePtr         vnp;
  ValNodePtr         vnp2;

  if (awp == NULL || target == NULL || bsp == NULL) return;
  ajp = awp->ajp;
  if (ajp == NULL) return;

  if (cds != NULL) {
    slp = AsnIoMemCopy ((Pointer) cds->location, (AsnReadFunc) SeqLocAsnRead, (AsnWriteFunc) SeqLocAsnWrite);
    if (slp != NULL) {
      for (slpx = SeqLocFindNext (slp, NULL); slpx != NULL; slpx = SeqLocFindNext (slp, slpx)) {
        if (slpx->choice == SEQLOC_INT) {
          sintp = (SeqIntPtr) slpx->data.ptrvalue;
          if (sintp != NULL) {
            sintp->strand = Seq_strand_both;
          }
        } else if (slpx->choice == SEQLOC_PNT) {
          spp = (SeqPntPtr) slpx->data.ptrvalue;
          if (spp != NULL) {
            spp->strand = Seq_strand_both;
          }
        }
      }
    }
    sfp = SeqMgrGetOverlappingSource (slp, &fcontext);
    SeqLocFree (slp);
    if (sfp != NULL) {
      biop = (BioSourcePtr) sfp->data.value.ptrvalue;
      bbp = AddSource (awp, &(awp->srchead), biop, sfp->comment);
      if (bbp != NULL) {

        bbp->entityID = sfp->idx.entityID;
        bbp->itemID = sfp->idx.itemID;
        bbp->itemtype = OBJ_SEQFEAT;

        isp = (IntSrcBlockPtr) bbp;
        CheckSeqLocForPartial (sfp->location, &noLeft, &noRight);
        hasNulls = LocationHasNullsBetween (sfp->location);
        isp->loc = SeqLocMerge (target, sfp->location, NULL, FALSE, TRUE, hasNulls);
        SetSeqLocPartial (isp->loc, noLeft, noRight);
        isp->left = fcontext.left;
        isp->right = fcontext.right;
        isp->comment = sfp->comment;
      }
    }

    return;
  }

  if (awp->format != FTABLE_FMT || awp->mode == DUMP_MODE) {

    /* full length loc for descriptors */
  
    sint.from = 0;
    if (ajp->ajp.slp != NULL) {
      sint.to = SeqLocLen (ajp->ajp.slp) - 1;
    } else {
      sint.to = bsp->length - 1;
    }
    sint.strand = Seq_strand_plus;
    sint.id = SeqIdStripLocus (SeqIdDup (SeqIdFindBest (bsp->id, 0)));
    sint.if_from = NULL;
    sint.if_to = NULL;
  
    vn.choice = SEQLOC_INT;
    vn.data.ptrvalue = (Pointer) &sint;
    vn.next = NULL;
  
    /* if SWISS-PROT, may have multiple source descriptors */
  
    if (ISA_aa (bsp->mol)) {
      for (sip = bsp->id; sip != NULL; sip = sip->next) {
        if (sip->choice == SEQID_SWISSPROT) {
          loop = TRUE;
        } else if (sip->choice == SEQID_OTHER) {
          tsip = (TextSeqIdPtr) sip->data.ptrvalue;
          if (tsip != NULL && StringNICmp (tsip->accession, "WP_", 3) == 0) {
            is_wp = TRUE;
          }
        }
      }
    }

    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
    while (sdp != NULL) {
      ValNodeAddPointer (&sdplist, 0, (Pointer) sdp);
      biop = (BioSourcePtr) sdp->data.ptrvalue;
      if (biop != NULL) {
        orp = biop->org;
        if (orp != NULL) {
          onp = orp->orgname;
          if (onp != NULL) {
            if (onp->choice == 5) {
              for (tep = (TaxElementPtr) onp->data; tep != NULL; tep = tep->next) {
                if (tep->fixed_level == 0 && StringICmp (tep->level, "superkingdom") == 0) {
                  num_super_kingdom++;
                  if (super_kingdom_name == NULL) {
                    super_kingdom_name = tep->name;
                  } else if (StringICmp (super_kingdom_name, tep->name) != 0) {
                    super_kingdoms_different = TRUE;
                  }
                }
              }
            }
          }
        }
      }
      sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_source, &dcontext);
    }

    vnp = sdplist;
    while (vnp != NULL) {
      sdp = (SeqDescrPtr) vnp->data.ptrvalue;

      if (sdp != NULL && sdp->extended != 0) {
        ovp = (ObjValNodePtr) sdp;

        /* check if descriptor on part already added on segmented bioseq */

        okay = TRUE;
        for (vnp2 = awp->srchead; vnp2 != NULL && okay; vnp2 = vnp2->next) {
          bbp = (BaseBlockPtr) vnp2->data.ptrvalue;
          if (bbp != NULL) {
            if (bbp->entityID == ovp->idx.entityID &&
                bbp->itemID == ovp->idx.itemID &&
                bbp->itemtype == OBJ_SEQDESC) {
              okay = FALSE;
            }
          }
        }
  
        if (okay) {
          biop = (BioSourcePtr) sdp->data.ptrvalue;
          bbp = AddSource (awp, &(awp->srchead), biop, NULL);
          if (bbp != NULL) {
  
            bbp->entityID = ovp->idx.entityID;
            bbp->itemID = ovp->idx.itemID;
            bbp->itemtype = OBJ_SEQDESC;

            isp = (IntSrcBlockPtr) bbp;
            isp->loc = SeqLocMerge (target, &vn, NULL, FALSE, TRUE, FALSE);
            isp->left = 0;
            isp->right = bsp->length - 1;
            isp->is_descriptor = TRUE;
          }
        }
      }

      if ((num_super_kingdom > 1 && super_kingdoms_different && is_wp) || loop) {
        vnp = vnp->next;
      } else {
        vnp = NULL;
      }
    }

    SeqIdFree (sint.id);
  }

  ValNodeFree (sdplist);

  if ((! awp->contig) || awp->showconsource) {

    /* features are indexed on parent if segmented */

    bsp = awp->parent;

    sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_BIOSRC, 0, &fcontext);
    while (sfp != NULL) {
      ivals = fcontext.ivals;
      numivals = fcontext.numivals;
      if (ivals != NULL && numivals > 0) {

        idx = (numivals - 1) * 2;
        start = ivals [idx];
        stop = ivals [idx + 1];
        if (stop >= from && stop <= to && (ajp->ajp.slp == NULL || SeqLocCompare (sfp->location, ajp->ajp.slp) > 0)) {

          biop = (BioSourcePtr) sfp->data.value.ptrvalue;
          bbp = AddSource (awp, &(awp->srchead), biop, sfp->comment);
          if (bbp != NULL) {

            bbp->entityID = fcontext.entityID;
            bbp->itemID = fcontext.itemID;
            bbp->itemtype = OBJ_SEQFEAT;

            isp = (IntSrcBlockPtr) bbp;
            if (sfp->location != NULL && sfp->location->choice == SEQLOC_PNT) {
              isp->loc = AsnIoMemCopy ((Pointer) sfp->location,
                                       (AsnReadFunc) SeqLocAsnRead,
                                       (AsnWriteFunc) SeqLocAsnWrite);
            } else if (s_isFuzzyLoc (sfp->location)) {
              isp->loc = AsnIoMemCopy ((Pointer) sfp->location,
                                      (AsnReadFunc) SeqLocAsnRead,
                                      (AsnWriteFunc) SeqLocAsnWrite);
            } else if (SeqLocId(sfp->location) == NULL) {
              isp->loc = AsnIoMemCopy ((Pointer) sfp->location,
                                       (AsnReadFunc) SeqLocAsnRead,
                                       (AsnWriteFunc) SeqLocAsnWrite);
            } else {
              CheckSeqLocForPartial (sfp->location, &noLeft, &noRight);
              hasNulls = LocationHasNullsBetween (sfp->location);
              isp->loc = SeqLocMerge (target, sfp->location, NULL, FALSE, TRUE, hasNulls);
              SetSeqLocPartial (isp->loc, noLeft, noRight);
            }
            isp->left = fcontext.left;
            isp->right = fcontext.right;
            isp->comment = sfp->comment;
            if (ajp->ajp.slp != NULL) {
              sip = SeqIdParse ("lcl|dummy");
              left = GetOffsetInBioseq (ajp->ajp.slp, bsp, SEQLOC_LEFT_END);
              right = GetOffsetInBioseq (ajp->ajp.slp, bsp, SEQLOC_RIGHT_END);
              strand = SeqLocStrand (ajp->ajp.slp);
              split = FALSE;
              newloc = SeqLocReMapEx (sip, ajp->ajp.slp, isp->loc, 0, FALSE, ajp->masterStyle, ajp->relaxedMapping);
              /*
              newloc = SeqLocCopyRegion (sip, isp->loc, bsp, left, right, strand, &split);
              */
              SeqIdFree (sip);
              if (newloc != NULL) {
                A2GBSeqLocReplaceID (newloc, ajp->ajp.slp);
                isp->loc = SeqLocFree (isp->loc);
                isp->loc = newloc;
                isp->left = left;
                isp->right = right;
              }
            }
          }
        }
      }

      sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_BIOSRC, 0, &fcontext);
    }
  }
}

static Boolean LIBCALLBACK GetSourcesOnSeg (
  SeqLocPtr slp,
  SeqMgrSegmentContextPtr context
)

{
  Asn2gbWorkPtr  awp;
  BioseqPtr      bsp;
  Int4           from;
  SeqLocPtr      loc;
  SeqEntryPtr    oldscope;
  SeqEntryPtr    sep;
  SeqIdPtr       sip;
  Int4           to;

  if (slp == NULL || context == NULL) return FALSE;
  awp = (Asn2gbWorkPtr) context->userdata;

  from = context->cumOffset;
  to = from + context->to - context->from;

  sip = SeqLocId (slp);
  if (sip == NULL) {
    loc = SeqLocFindNext (slp, NULL);
    if (loc != NULL) {
      sip = SeqLocId (loc);
    }
  }
  if (sip == NULL) return TRUE;

  /* biosource descriptors only on parts within entity */

  sep = GetTopSeqEntryForEntityID (awp->entityID);
  oldscope = SeqEntrySetScope (sep);
  bsp = BioseqFind (sip);
  SeqEntrySetScope (oldscope);

  if (bsp != NULL) {
    GetSourcesOnBioseq (awp, awp->target, bsp, from, to, NULL);
    return TRUE;
  }

  /* if we ever want to fetch remote sources, code goes here */

#if 0
  Uint2          entityID;

  /* may remote fetch genome component if not already in memory */

  bsp = BioseqLockById (sip);

  if (bsp == NULL) return TRUE;

  entityID = ObjMgrGetEntityIDForPointer (bsp);

  if (entityID != awp->entityID) {

    /* if segment not packaged in record, may need to feature index it */

    if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
      SeqMgrIndexFeatures (entityID, NULL);
    }

    /* collect features indexed on the remote bioseq */

    from = 0;
    to = bsp->length - 1;
  }

  GetSourcesOnBioseq (awp, awp->target, bsp, from, to, NULL);

  BioseqUnlock (bsp);
#endif

  return TRUE;
}

/* isIdenticalSource() -- Checks to see if two sources are identical */
/*                        by comparing the actual values in the      */
/*                        fields.  This only gets called if the two  */
/*                        sources hashed the same -- it's a double-  */
/*                        check since two non-identical things will  */
/*                        occassionally hash to the same value.      */
/*                        Now checks for adjacent or overlapping.    */

static Boolean isIdenticalSource (IntSrcBlockPtr isp1, IntSrcBlockPtr isp2)
{
  OrgModPtr     omp1;
  OrgModPtr     omp2;
  SubSourcePtr  ssp1;
  SubSourcePtr  ssp2;
  ValNodePtr    vnp1;
  ValNodePtr    vnp2;
  ObjectIdPtr   oip1;
  ObjectIdPtr   oip2;
  DbtagPtr      dbt1;
  DbtagPtr      dbt2;

  if (isp1->is_focus != isp2->is_focus)
    return FALSE;

  /* Compare the taxonomy names */

  if (StringICmp(isp1->taxname,isp2->taxname) != 0)
    return FALSE;

  /* Compare the comment */

  if (StringICmp(isp1->comment,isp2->comment) != 0)
    return FALSE;

  /* Compare the org mods */

  omp1 = isp1->omp;
  omp2 = isp2->omp;
  while (omp1 != NULL && omp2 != NULL)
    {
      if (omp1->subtype != omp2->subtype)
        return FALSE;
      if (StringICmp (omp1->subname, omp2->subname) != 0)
        return FALSE;
      omp1 = omp1->next;
      omp2 = omp2->next;
    }

  if (omp1 != NULL || omp2 != NULL)
    return FALSE;

  /* Compare the subtypes */

  ssp1 = isp1->ssp;
  ssp2 = isp2->ssp;

  while (ssp1 != NULL && ssp2 != NULL)
    {
      if (ssp1->subtype != ssp2->subtype)
        return FALSE;
      if (StringICmp(ssp1->name, ssp2->name) != 0)
        return FALSE;
      ssp1 = ssp1->next;
      ssp2 = ssp2->next;
    }

  if (ssp1 != NULL || ssp2 != NULL)
    return FALSE;

  /* Compare the DB tags */

  vnp1 = isp1->vnp;
  vnp2 = isp2->vnp;

  while (vnp1 != NULL && vnp2 != NULL)
    {
      dbt1 = (DbtagPtr) vnp1->data.ptrvalue;
      dbt2 = (DbtagPtr) vnp2->data.ptrvalue;

      if ((dbt1 != NULL) && (dbt2 != NULL)) {
        if (StringCmp (dbt1->db, dbt2->db) != 0)
          return FALSE;

        oip1 = dbt1->tag;
        oip2 = dbt2->tag;
        if ((oip1 != NULL) && (oip2 != NULL)) {
          if (oip1->str != NULL) {
            if (StringICmp(oip1->str, oip2->str) != 0)
              return FALSE;
          } else  {
            if (oip1->id != oip2->id)
              return FALSE;
          }
        }
        else if (oip1 != NULL)
          return FALSE;
        else if (oip2 != NULL)
          return FALSE;
      }
      else if (dbt1 != NULL)
        return FALSE;
      else if (dbt2 != NULL)
        return FALSE;

      vnp1 = vnp1->next;
      vnp2 = vnp2->next;
    }

  if (vnp1 != NULL || vnp2 != NULL)
    return FALSE;

  /* now check for not adjacent or overlapping */

  if (isp2->right + 1 < isp1->left) return FALSE;

  /* If it passed all checks, then they */
  /* are the same, so return true.      */

  return TRUE;
}

static void CleanupPackedSeqInt (SeqLocPtr location)

{
  SeqLocPtr  head = NULL;
  SeqIntPtr  loc;
  SeqIntPtr  sintp;
  SeqLocPtr  slp;

  if (location == NULL || location->choice != SEQLOC_PACKED_INT || location->data.ptrvalue == NULL) return;

  slp = SeqLocFindNext (location, NULL);
  while (slp != NULL) {
    if (slp->choice == SEQLOC_INT) {
      sintp = (SeqIntPtr) slp->data.ptrvalue;
      if (sintp != NULL) {
        loc = AsnIoMemCopy (sintp, (AsnReadFunc) SeqIntAsnRead,
                            (AsnWriteFunc) SeqIntAsnWrite);
        ValNodeAddPointer (&head, SEQLOC_INT, loc);
      }
    }
    slp = SeqLocFindNext (location, slp);
  }
  if (head == NULL) return;

  location->data.ptrvalue = SeqLocFree (location->data.ptrvalue);
  location->data.ptrvalue = head;

  slp = location->data.ptrvalue;
  if (slp == NULL || slp->next != NULL) return;
    /* here seqloc_packed_int points to a single location element, so no need for seqloc_packed_int parent */
    location->choice = slp->choice;
    location->data.ptrvalue = (Pointer) slp->data.ptrvalue;
    MemFree (slp);
}

static Boolean x_NotSpecialTaxName (
  CharPtr taxname
)

{
  if (StringHasNoText (taxname)) return TRUE;

  if (StringICmp (taxname, "synthetic construct") == 0) return FALSE;
  if (StringICmp (taxname, "artificial sequence") == 0) return FALSE;
  if (StringStr (taxname, "vector") != NULL) return FALSE;
  if (StringStr (taxname, "Vector") != NULL) return FALSE;

  return TRUE;
}

NLM_EXTERN void AddSourceFeatBlock (
  Asn2gbWorkPtr awp
)

{
  IntAsn2gbJobPtr    ajp;
  Asn2gbSectPtr      asp;
  BaseBlockPtr       bbp;
  BioSourcePtr       biop;
  BioseqPtr          bsp;
  SeqFeatPtr         cds;
  SeqMgrFeatContext  context;
  BioseqPtr          dna;
  SeqLocPtr          duploc;
  Boolean            excise;
  GBFeaturePtr       gbfeat = NULL;
  GBSeqPtr           gbseq;
  ValNodePtr         head = NULL;
  IntSrcBlockPtr     isp;
  IntSrcBlockPtr     lastisp;
  IntSrcBlockPtr     descrIsp;
  ValNodePtr         next;
  OrgRefPtr          orp;
  Char               pfx [128], sfx [128];
  ValNodePtr         PNTR prev;
  SeqDescrPtr        sdp;
  SeqInt             sint;
  SeqIdPtr           sip;
  SeqLocPtr          slp;
  Int4               source_count = 0;
  CharPtr            str;
  BioseqPtr          target = NULL;
  CharPtr            taxname;
  ValNode            vn;
  ValNodePtr         vnp;
  Boolean            descHasFocus = FALSE;
  StringItemPtr      ffstring;

  if (awp == NULL) return;
  ajp = awp->ajp;
  if (ajp == NULL) return;
  asp = awp->asp;
  if (asp == NULL) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;

  ffstring = FFGetString(ajp);
  if ( ffstring == NULL ) return;

  pfx [0] = '\0';
  sfx [0] = '\0';

  /* collect biosources on bioseq */

  awp->srchead = NULL;

  if (ISA_aa (bsp->mol)) {

    /* if protein, get sources applicable to DNA location of CDS */

    sdp = GetNextDescriptorUnindexed (bsp, Seq_descr_source, NULL);
    if (sdp != NULL && sdp->choice == Seq_descr_source) {
      biop = (BioSourcePtr) sdp->data.ptrvalue;
      if (biop != NULL) {
        orp = biop->org;
        if (orp != NULL) {
          taxname = orp->taxname;
          if (StringHasNoText (taxname) || x_NotSpecialTaxName (taxname)) {
            cds = SeqMgrGetCDSgivenProduct (bsp, &context);
            if (cds != NULL) {
              dna = BioseqFindFromSeqLoc (cds->location);
              if (dna != NULL) {
                GetSourcesOnBioseq (awp, dna, dna, context.left, context.right, cds);
                target = dna;
              }
            }
          }
        }
      }
    }
  }

  if (awp->srchead == NULL) {
    GetSourcesOnBioseq (awp, bsp, bsp, awp->from, awp->to, NULL);
    target = bsp;
  }

  if (bsp->repr == Seq_repr_seg) {

    /* collect biosource descriptors on local parts */

    SeqMgrExploreSegments (bsp, (Pointer) awp, GetSourcesOnSeg);
    target = awp->target;
  }

  head = awp->srchead;
  awp->srchead = NULL;

  if (head == NULL && (awp->format != FTABLE_FMT || awp->mode == DUMP_MODE)) {

    if (ajp->gbseq) {
      gbseq = &asp->gbseq;
    } else {
      gbseq = NULL;
    }

    sint.from = 0;
    sint.to = bsp->length - 1;
    sint.strand = Seq_strand_plus;
    sint.id = SeqIdStripLocus (SeqIdDup (SeqIdFindBest (bsp->id, 0)));
    sint.if_from = NULL;
    sint.if_to = NULL;

    vn.choice = SEQLOC_INT;
    vn.data.ptrvalue = (Pointer) &sint;
    vn.next = NULL;

    FFStartPrint (ffstring, awp->format, 5, 21, NULL, 0, 5, 21, "FT", FALSE);

    /*
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_GI) {
        currGi = (BIG_ID) sip->data.intvalue;
      }
    }
    */
    if (GetWWW (ajp) && ajp->mode == ENTREZ_MODE && ajp->seqspans &&
        (ajp->format == GENBANK_FMT || ajp->format == GENPEPT_FMT)) {
      sprintf (pfx, "<span id=\"feature_%s_source_0\" class=\"feature\">", awp->currAccVerLabel);
    }

    FFAddOneString(ffstring, "source", FALSE, FALSE, TILDE_IGNORE);
    FFAddNChar(ffstring, ' ', 21 - 5 - StringLen("source"), FALSE);

    if (gbseq != NULL) {
      gbfeat = GBFeatureNew ();
      if (gbfeat != NULL) {
        gbfeat->key = StringSave ("source");
      }
    }

    str = FFFlatLoc (ajp, bsp, &vn, (Boolean) (awp->style == MASTER_STYLE), FALSE);
    if ( GetWWW(ajp) ) {
      FF_www_featloc (ffstring, str);
    } else {
      FFAddOneString (ffstring, str, FALSE, FALSE, TILDE_IGNORE);
    }

    if (gbseq != NULL) {
      if (gbfeat != NULL) {
        if (! StringHasNoText (str)) {
          gbfeat->location = StringSave (str);
        } else {
          gbfeat->location = StringSave ("");
        }
      }
    }

    MemFree (str);

    if (ajp->flags.needOrganismQual) {
      FFAddNewLine(ffstring);
      FFAddTextToString (ffstring, "/organism=\"", "unknown", "\"", FALSE, TRUE, TILDE_TO_SPACES);
#ifdef ASN2GNBK_PRINT_UNKNOWN_ORG
    } else {
      FFAddNewLine(ffstring);
      FFAddTextToString (ffstring, "/organism=\"", "unknown", "\"", FALSE, TRUE, TILDE_TO_SPACES);
#endif
    }

    str = GetMolTypeQual (bsp);
    if (StringICmp (str, "ncRNA") == 0) {
      str = "other RNA";
    }
    if (str == NULL) {
      switch (bsp->mol) {
        case Seq_mol_dna :
          str = "unassigned DNA";
          break;
        case Seq_mol_rna :
          str = "unassigned RNA";
          break;
        case Seq_mol_aa :
          break;
        default :
          str = "unassigned DNA";
          break;
      }
    }
    if (str != NULL) {
      FFAddNewLine(ffstring);
      FFAddTextToString (ffstring, "/mol_type=\"", str, "\"", FALSE, TRUE, TILDE_TO_SPACES);
    }

    if (GetWWW (ajp) && ajp->mode == ENTREZ_MODE && ajp->seqspans &&
        (ajp->format == GENBANK_FMT || ajp->format == GENPEPT_FMT)) {
      sprintf (sfx, "</span>");
    }

    str = FFEndPrintEx (ajp, ffstring, awp->format, 5, 21, 5, 21, "FT", pfx, sfx);

    bbp = (BaseBlockPtr) Asn2gbAddBlock (awp, SOURCEFEAT_BLOCK, sizeof (IntSrcBlock));
    if (bbp != NULL) {
      bbp->section = awp->currsection;
      bbp->string = str;
    } else {
      MemFree(str);
    }
    FFRecycleString(ajp, ffstring);

    if (awp->afp != NULL) {
      DoImmediateFormat (awp->afp, (BaseBlockPtr) bbp);
    }

    /* optionally populate gbseq for XML-ized GenBank format */

    if (gbseq != NULL) {
      if (gbfeat != NULL) {
        AddFeatureToGbseq (gbseq, gbfeat, str, NULL);
      }
    }

    return;
  }

  if (head == NULL) return;

  /* sort by hash values */

  head = ValNodeSort (head, SortSourcesByHash);

  /* unique sources, excise duplicates from list */

  prev = &(head);
  vnp = head;
  lastisp = NULL;
  while (vnp != NULL) {
    excise = FALSE;
    next = vnp->next;
    isp = (IntSrcBlockPtr) vnp->data.ptrvalue;
    if (isp->is_descriptor && isp->is_focus)
      descHasFocus = TRUE;
    if (lastisp != NULL) {
      if (isp != NULL) {
        if (lastisp->is_focus == isp->is_focus &&
            lastisp->orghash == isp->orghash &&
            lastisp->xrfhash == isp->xrfhash) {

          /* check for identical modifiers */

          if (lastisp->modhash == isp->modhash &&
              lastisp->subhash == isp->subhash) {

            excise = isIdenticalSource (isp, lastisp);

          /* or modifiers only in lastisp (e.g., on part bioseq) */

          } else if (isp->modhash == 0 && isp->subhash == 0) {
            excise = isIdenticalSource (isp, lastisp);
          }
        }
      }
    }
    if (awp->mode == DUMP_MODE) {
      excise = FALSE;
    }
    /* does not fuse equivalent source features for local, general, refseq, and 2+6 genbank ids */
    if (excise && awp->sourcePubFuse) {
      *prev = vnp->next;
      vnp->next = NULL;

      /* combine locations of duplicate sources */

      if (lastisp != NULL) {
        slp = SeqLocMerge (target, lastisp->loc, isp->loc, FALSE, TRUE, FALSE);
        lastisp->loc = SeqLocFree (lastisp->loc);
        lastisp->loc = slp;
        lastisp->left = MIN (lastisp->left,isp->left);
        lastisp->right = MAX (lastisp->right, isp->right);
      }

      /* and remove duplicate source */

      SeqLocFree (isp->loc);
      MemFree (isp);
      ValNodeFree (vnp);

    } else {

      prev = &(vnp->next);
      lastisp = isp;
    }
    vnp = next;
  }

  /* Sort again, by location this time */

  head = ValNodeSort (head, SortSourcesByPos);

  /* If the descriptor has a focus, then subtract */
  /* out all the other source locations.          */

  descrIsp = (IntSrcBlockPtr) head->data.ptrvalue; /* Sorted 1st by now */

  if ((descHasFocus) && (! descrIsp->is_synthetic)) {

    vnp = head;
    duploc = AsnIoMemCopy ((Pointer) descrIsp->loc,
                           (AsnReadFunc) SeqLocAsnRead,
                           (AsnWriteFunc) SeqLocAsnWrite);
    vnp = vnp->next;
    while (vnp != NULL) {
      isp = (IntSrcBlockPtr) vnp->data.ptrvalue;
      if (SeqLocAinB (descrIsp->loc, isp->loc) >= 0) {
        vnp = NULL; /* break the chain */
        descrIsp->loc = SeqLocFree (descrIsp->loc);
        descrIsp->loc = duploc;
        duploc = NULL;
      } else {
        descrIsp->loc = SeqLocSubtract (descrIsp->loc, isp->loc);
        vnp = vnp->next;
      }
    }
    CleanupPackedSeqInt (descrIsp->loc);
    descrIsp->left  = SeqLocStart (descrIsp->loc);
    descrIsp->right = SeqLocStop (descrIsp->loc);
    SeqLocFree (duploc);
  }

  /* if features completely subtracted descriptor
     intervals, suppress in release, entrez modes */

  if (descrIsp->loc == NULL && ajp->flags.hideEmptySource && head->next != NULL) {
    vnp = head->next;
    head->next = NULL;
    ValNodeFreeData (head);
    head = vnp;
  }

  /* finally link into blocks for current section */

  ValNodeLink (&(awp->lastblock), head);
  vnp = awp->lastblock;
  if (vnp == NULL) return;
  while (vnp->next != NULL) {
    vnp = vnp->next;
  }

  awp->lastblock = vnp;
  if (awp->blockList == NULL) {
    awp->blockList = vnp;
  }
  FFRecycleString(ajp, ffstring);

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    isp = (IntSrcBlockPtr) vnp->data.ptrvalue;
    if (isp == NULL) continue;
    isp->source_count = source_count;
    source_count++;
  }

  if (awp->afp != NULL) {
    for (vnp = head; vnp != NULL; vnp = vnp->next) {
      isp = (IntSrcBlockPtr) vnp->data.ptrvalue;
      if (isp == NULL) continue;
      DoImmediateFormat (awp->afp, (BaseBlockPtr) isp);
    }
  }

}

static Boolean IsCDD (
  SeqFeatPtr sfp
)

{
  DbtagPtr    dbt;
  ValNodePtr  vnp;

  for (vnp = sfp->dbxref; vnp != NULL; vnp = vnp->next) {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt != NULL && StringCmp (dbt->db, "CDD") == 0) return TRUE;
  }

  return FALSE;
}

NLM_EXTERN void SetIfpFeatCount (
  IntFeatBlockPtr ifp,
  IntAsn2gbJobPtr ajp,
  Asn2gbWorkPtr awp,
  Boolean isProt
)

{
  FeatBlockPtr      fbp;
  Uint1             featdeftype;
  IntAsn2gbSectPtr  iasp;
  Boolean           is_other = FALSE;

  if (ifp == NULL || ajp == NULL || awp == NULL) return;
  iasp = (IntAsn2gbSectPtr) awp->asp;
  if (iasp == NULL) return;

  fbp = (FeatBlockPtr) ifp;

  featdeftype = fbp->featdeftype;

  if (featdeftype == FEATDEF_COMMENT) {
    featdeftype = FEATDEF_misc_feature;
  }

  if (! isProt) {
    if (featdeftype == FEATDEF_REGION || featdeftype == FEATDEF_BOND || featdeftype == FEATDEF_SITE) {
      featdeftype = FEATDEF_misc_feature;
    }
  }

  if (ajp->format == GENPEPT_FMT && isProt) {
    if (ifp->mapToPep) {
      if (featdeftype >= FEATDEF_preprotein && featdeftype <= FEATDEF_transit_peptide_aa) {
        featdeftype = FEATDEF_preprotein;
      }
    }
  }

  if (featdeftype == FEATDEF_Imp_CDS) {
    featdeftype = FEATDEF_CDS;
  }
  if (featdeftype == FEATDEF_preRNA) {
    featdeftype = FEATDEF_precursor_RNA;
  }
  if (featdeftype == FEATDEF_otherRNA) {
    featdeftype = FEATDEF_misc_RNA;
  }
  if (featdeftype == FEATDEF_mat_peptide_aa) {
    featdeftype = FEATDEF_mat_peptide;
  }
  if (featdeftype == FEATDEF_sig_peptide_aa) {
    featdeftype = FEATDEF_sig_peptide;
  }
  if (featdeftype == FEATDEF_transit_peptide_aa) {
    featdeftype = FEATDEF_transit_peptide;
  }

  if (ajp->refseqConventions || awp->isRefSeq) {
    is_other = TRUE;
  }

  if (! isProt) {
    if (featdeftype == FEATDEF_preprotein) {
      if (! is_other) {
        featdeftype = FEATDEF_misc_feature;
      }
    }
  }

  if (featdeftype == FEATDEF_CLONEREF) {
    if (ajp->mode == RELEASE_MODE || ajp->mode == ENTREZ_MODE) {
      featdeftype = FEATDEF_misc_feature;
    }
  }

  if (featdeftype == FEATDEF_repeat_unit && (ajp->mode == RELEASE_MODE || ajp->mode == ENTREZ_MODE)) {
    featdeftype = FEATDEF_repeat_region;
  }

  if (featdeftype < FEATDEF_MAX) {
    ifp->feat_count = iasp->feat_counts [featdeftype];
    (iasp->feat_counts [featdeftype])++;
  }
}

static void GetFeatsOnCdsProduct (
  SeqFeatPtr cds,
  BioseqPtr nbsp,
  BioseqPtr pbsp,
  IntAsn2gbJobPtr ajp,
  Asn2gbWorkPtr awp
)

{
  FeatBlockPtr       fbp;
  IntFeatBlockPtr    ifp;
  Boolean            isRefSeq;
  Int4               lastleft;
  Int4               lastright;
  SeqAnnotPtr        lastsap;
  SeqFeatPtr         lastsfp;
  SeqLocPtr          location;
  SeqLocPtr          newloc;
  SeqMgrFeatContext  pcontext;
  SeqFeatPtr         prt;
  SeqIdPtr           sip;
  SeqLocPtr          slp;
  Boolean            suppress;

  if (cds == NULL || ajp == NULL || awp == NULL) return;
  if (nbsp == NULL || pbsp == NULL || (! ISA_aa (pbsp->mol))) return;

  if (awp->hideCdsProdFeats) return;

  isRefSeq = FALSE;
  for (sip = nbsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_OTHER) {
      isRefSeq = TRUE;
    }
  }

  /* explore mat_peptides, sites, etc. */

  lastsfp = NULL;
  lastsap = NULL;
  lastleft = 0;
  lastright = 0;

  prt = SeqMgrGetNextFeature (pbsp, NULL, 0, 0, &pcontext);
  while (prt != NULL) {

    if (pcontext.featdeftype == FEATDEF_REGION ||
        pcontext.featdeftype == FEATDEF_SITE ||
        pcontext.featdeftype == FEATDEF_BOND ||
        pcontext.featdeftype == FEATDEF_mat_peptide_aa ||
        pcontext.featdeftype == FEATDEF_sig_peptide_aa ||
        pcontext.featdeftype == FEATDEF_transit_peptide_aa ||
        pcontext.featdeftype == FEATDEF_preprotein ||
        (pcontext.featdeftype == FEATDEF_propeptide_aa /* && isRefSeq */)) {

      if (awp->hideSitesBondsRegions && (pcontext.featdeftype == FEATDEF_REGION ||
                                         pcontext.featdeftype == FEATDEF_SITE ||
                                         pcontext.featdeftype == FEATDEF_BOND)) {

        /* hide site, bond, and region features */

      } else if (awp->hideCddFeats && pcontext.featdeftype == FEATDEF_REGION && IsCDD (prt)) {

        /* passing this test prevents mapping of COG CDD region features */

      } else if (pcontext.dnaStop >= awp->from && pcontext.dnaStop <= awp->to) {

        /* suppress duplicate features (on protein) */

        suppress = FALSE;
        if (lastsfp != NULL && lastsap != NULL) {
          if (lastsfp->idx.subtype == prt->idx.subtype &&
              lastleft == pcontext.left &&
              lastright == pcontext.right) {
              if (lastsap == pcontext.sap ||
                  (lastsap->desc == NULL && pcontext.sap->desc == NULL)) {
              if (AsnIoMemComp (lastsfp, prt, (AsnWriteFunc) SeqFeatAsnWrite)) {
                suppress = TRUE;
              }
            }
          }
        }

        /* make sure feature maps within nucleotide sublocation */

        if (! suppress) {
          if (ajp->ajp.slp != NULL) {
            location = aaFeatLoc_to_dnaFeatLoc (cds, prt->location);
            slp = SeqLocMerge (nbsp, location, NULL, FALSE, TRUE, FALSE);
            if (slp != NULL) {
              sip = SeqIdParse ("lcl|dummy");
              newloc = SeqLocReMapEx (sip, ajp->ajp.slp, slp, 0, FALSE, ajp->masterStyle, ajp->relaxedMapping);
              SeqIdFree (sip);
              SeqLocFree (slp);
              if (newloc == NULL) {
                suppress = TRUE;
              }
              SeqLocFree (newloc);
            } else {
              suppress = TRUE;
            }
            SeqLocFree (location);
          }
        }

        if (! suppress) {

          fbp = (FeatBlockPtr) Asn2gbAddBlock (awp, FEATURE_BLOCK, sizeof (IntFeatBlock));
          if (fbp != NULL) {

            fbp->entityID = pcontext.entityID;
            fbp->itemID = pcontext.itemID;
            fbp->itemtype = OBJ_SEQFEAT;
            fbp->featdeftype = pcontext.featdeftype;
            ifp = (IntFeatBlockPtr) fbp;
            ifp->mapToNuc = TRUE;
            ifp->mapToProt = FALSE;
            ifp->mapToGen = FALSE;
            ifp->mapToMrna = FALSE;
            ifp->mapToPep = FALSE;
            ifp->left = 0;
            ifp->right = 0;
            SetIfpFeatCount (ifp, ajp, awp, FALSE);
            ifp->firstfeat = awp->firstfeat;
            awp->firstfeat = FALSE;

            if (awp->afp != NULL) {
              DoImmediateFormat (awp->afp, (BaseBlockPtr) fbp);
            }
          }
        }

        lastsfp = prt;
        lastsap = pcontext.sap;
        lastleft = pcontext.left;
        lastright = pcontext.right;

      }
    }
    prt = SeqMgrGetNextFeature (pbsp, prt, 0, 0, &pcontext);
  }
}

static void GetRemoteFeatsOnCdsProduct (
  SeqFeatPtr cds,
  BioseqPtr nbsp,
  BioseqPtr pbsp,
  IntAsn2gbJobPtr ajp,
  Asn2gbWorkPtr awp
)

{
  BioseqPtr        bsp;
  FeatBlockPtr     fbp;
  ValNodePtr       head = NULL;
  IntFeatBlockPtr  ifp;
  Boolean          isRefSeq;
  Int4             lastleft;
  Int4             lastright;
  SeqAnnotPtr      lastsap;
  SeqFeatPtr       lastsfp;
  SeqLocPtr        location;
  SeqLocPtr        newloc;
  SeqFeatPtr       prt;
  ValNodePtr       publist;
  Asn2gbFreeFunc   remotefree;
  Asn2gbLockFunc   remotelock;
  ValNodePtr       remotevnp;
  SeqAnnotPtr      sap;
  SeqFeatPtr       sfp;
  SeqIdPtr         sip;
  SeqLocPtr        slp;
  Boolean          suppress;
  ValNodePtr       vnp;

  if (cds == NULL || ajp == NULL || awp == NULL) return;
  if (nbsp == NULL || pbsp == NULL || (! ISA_aa (pbsp->mol))) return;

  if (awp->hideCdsProdFeats) return;

  if (ajp->remotelock == NULL) return;

  remotelock = ajp->remotelock;
  remotefree = ajp->remotefree;

  sip = SeqIdFindBest (pbsp->id, SEQID_GI);
  if (sip == NULL) return;

  remotevnp = remotelock (sip, ajp->remotedata);
  if (remotevnp == NULL) return;

  /* do cleanup of remotely fetched feature tables */

  for (vnp = remotevnp; vnp != NULL; vnp = vnp->next) {
    bsp = (BioseqPtr) vnp->data.ptrvalue;
    if (bsp == NULL) continue;
    for (sap = bsp->annot; sap != NULL; sap = sap->next) {
      if (sap->type != 1) continue;
      for (sfp = (SeqFeatPtr) sap->data; sfp != NULL; sfp = sfp->next) {
        publist = NULL;
        CleanUpSeqFeat (sfp, FALSE, FALSE, TRUE, TRUE, &publist);
        sfp->idx.subtype = FindFeatDefType (sfp);
        ValNodeFreeData (publist);
        ValNodeAddPointer (&head, 0, (Pointer) sfp);
      }
    }
  }

  if (head == NULL) return;

  isRefSeq = FALSE;
  for (sip = nbsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_OTHER) {
      isRefSeq = TRUE;
    }
  }

  /* explore mat_peptides, sites, etc. */

  lastsfp = NULL;
  lastsap = NULL;
  lastleft = 0;
  lastright = 0;

  for (vnp = head; vnp != NULL; vnp = vnp->next) {

    prt = (SeqFeatPtr) vnp->data.ptrvalue;
    if (prt == NULL) continue;

    if (prt->idx.subtype == FEATDEF_REGION ||
        prt->idx.subtype == FEATDEF_SITE ||
        prt->idx.subtype == FEATDEF_BOND ||
        prt->idx.subtype == FEATDEF_mat_peptide_aa ||
        prt->idx.subtype == FEATDEF_sig_peptide_aa ||
        prt->idx.subtype == FEATDEF_transit_peptide_aa ||
        prt->idx.subtype == FEATDEF_preprotein ||
        (prt->idx.subtype == FEATDEF_propeptide_aa /* && isRefSeq */)) {

      if (awp->hideSitesBondsRegions && (prt->idx.subtype == FEATDEF_REGION ||
                                         prt->idx.subtype == FEATDEF_SITE ||
                                         prt->idx.subtype == FEATDEF_BOND)) {

        /* hide site, bond, and region features */

      } else if (awp->hideCddFeats && prt->idx.subtype == FEATDEF_REGION && IsCDD (prt)) {

        /* passing this test prevents mapping of COG CDD region features */

      } else {

        suppress = FALSE;

        /* make sure feature maps within nucleotide sublocation */

        if (! suppress) {
          if (ajp->ajp.slp != NULL) {
            location = aaFeatLoc_to_dnaFeatLoc (cds, prt->location);
            slp = SeqLocMerge (nbsp, location, NULL, FALSE, TRUE, FALSE);
            if (slp != NULL) {
              sip = SeqIdParse ("lcl|dummy");
              newloc = SeqLocReMapEx (sip, ajp->ajp.slp, slp, 0, FALSE, ajp->masterStyle, ajp->relaxedMapping);
              SeqIdFree (sip);
              SeqLocFree (slp);
              if (newloc == NULL) {
                suppress = TRUE;
              }
              SeqLocFree (newloc);
            } else {
              suppress = TRUE;
            }
            SeqLocFree (location);
          }
        }

        if (! suppress) {

          fbp = (FeatBlockPtr) Asn2gbAddBlock (awp, FEATURE_BLOCK, sizeof (IntFeatBlock));
          if (fbp != NULL) {

            fbp->entityID = 0;
            fbp->itemID = 0;
            fbp->itemtype = OBJ_SEQFEAT;
            fbp->featdeftype = prt->idx.subtype;
            ifp = (IntFeatBlockPtr) fbp;
            ifp->mapToNuc = TRUE;
            ifp->mapToProt = FALSE;
            ifp->mapToGen = FALSE;
            ifp->mapToMrna = FALSE;
            ifp->mapToPep = FALSE;
            ifp->left = 0;
            ifp->right = 0;
            SetIfpFeatCount (ifp, ajp, awp, FALSE);
            ifp->firstfeat = awp->firstfeat;
            awp->firstfeat = FALSE;

            if (awp->afp != NULL) {
              DoImmediateRemoteFeatureFormat (awp->afp, (BaseBlockPtr) fbp, prt);
            }
          }
        }
      }
    }
  }

  ValNodeFree (head);

  if (remotefree != NULL) {
    remotefree (remotevnp, ajp->remotedata);
  } else {
    /* otherwise free Bioseqs and ValNode chain ourselves */
    for (vnp = remotevnp; vnp != NULL; vnp = vnp->next) {
      bsp = (BioseqPtr) vnp->data.ptrvalue;
      if (bsp != NULL) {
        BioseqFree (bsp);
      }
    }
    ValNodeFree (remotevnp);
  }
}

static Boolean NotEMBLorDDBJ (
  BioseqPtr bsp
)

{
  SeqIdPtr  sip;

  if (bsp == NULL) return TRUE;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_EMBL || sip->choice == SEQID_TPE) return FALSE;
    if (sip->choice == SEQID_DDBJ || sip->choice == SEQID_TPD) return FALSE;
  }
  return TRUE;
}

/*
static Boolean EquivProtFeats (
  SeqFeatPtr prot1,
  SeqFeatPtr prot2
)

{
  ProtRefPtr  prp1, prp2;

  if (prot1 == NULL || prot2 == NULL) return FALSE;
  prp1 = (ProtRefPtr) prot1->data.value.ptrvalue;
  prp2 = (ProtRefPtr) prot2->data.value.ptrvalue;
  if (prp1 == NULL || prp2 == NULL) return FALSE;

  if (! AsnIoMemComp (prp1, prp2, (AsnWriteFunc) ProtRefAsnWrite)) return FALSE;

  if (StringDoesHaveText (prot1->comment) && StringDoesHaveText (prot2->comment)) {
    if (StringCmp (prot1->comment, prot2->comment) != 0) return FALSE;
  }

  return TRUE;
}
*/

/*
static Boolean EquivProtFeats (
  SeqFeatPtr prot1,
  SeqFeatPtr prot2
)

{
  SeqFeatPtr  cpy1, cpy2;
  Boolean     rsult = FALSE;
  SeqLocPtr   tmp;

  if (prot1 == NULL || prot2 == NULL) return FALSE;

  cpy1 = AsnIoMemCopy ((Pointer) prot1,
                       (AsnReadFunc) SeqFeatAsnRead,
                       (AsnWriteFunc) SeqFeatAsnWrite);
  cpy2 = AsnIoMemCopy ((Pointer) prot2,
                       (AsnReadFunc) SeqFeatAsnRead,
                       (AsnWriteFunc) SeqFeatAsnWrite);
  if (cpy1 == NULL || cpy2 == NULL) return FALSE;

  tmp = cpy1->location;
  cpy1->location = cpy2->location;

  rsult = AsnIoMemComp (cpy1, cpy2, (AsnWriteFunc) SeqFeatAsnWrite);

  cpy1->location = tmp;
  SeqFeatFree (cpy1);
  SeqFeatFree (cpy2);

  return rsult;
}
*/

static Boolean LocInBioseq (
  SeqLocPtr slp,
  BioseqPtr bsp
)

{
  SeqIdPtr  sip;

  if (slp == NULL || bsp == NULL) return FALSE;
  sip = SeqLocId (slp);
  if (sip == NULL) return FALSE;
  return SeqIdIn (sip, bsp->id);
}

static Boolean LIBCALLBACK GetFeatsOnBioseq (
  SeqFeatPtr sfp,
  SeqMgrFeatContextPtr fcontext
)

{
  IntAsn2gbJobPtr    ajp;
  Asn2gbSectPtr      asp;
  Asn2gbWorkPtr      awp;
  BioseqPtr          bsp;
  Char               buf [41];
  SeqFeatPtr         cds;
  SeqMgrFeatContext  cdscontext;
  FeatBlockPtr       fbp;
  SeqLocPtr          firstslp;
  SeqFeatPtr         gap;
  GBQualPtr          gbq;
  /*
  SeqFeatPtr         gene;
  */
  BIG_ID             gi;
  GeneRefPtr         grp;
  Boolean            has_est_len;
  Boolean            has_gap_type;
  IntCdsBlockPtr     icp;
  Int2               idx;
  IntFeatBlockPtr    ifp;
  IntPrtBlockPtr     ipp;
  Boolean            is_whole;
  Int4Ptr            ivals;
  Int2               j;
  Boolean            juststop = FALSE;
  SeqAnnotPtr        lastsap;
  SeqFeatPtr         lastsfp;
  SeqLocPtr          lastslp;
  SeqLocPtr          newloc;
  Int2               numivals;
  Boolean            okay;
  SeqEntryPtr        oldscope;
  BioseqPtr          parent;
  Boolean            partial5;
  Boolean            partial3;
  ValNodePtr         ppr;
  BioseqPtr          prod;
  ProtRefPtr         prp;
  Boolean            psdo;
  Boolean            pseudo = FALSE;
  RNAGenPtr          rgp;
  RnaRefPtr          rrp;
  SeqEntryPtr        sep;
  SeqIntPtr          sintp;
  SeqIdPtr           sip;
  SeqLocPtr          slp;
  Int4               start;
  Int4               stop;
  Boolean            supr;
  TextSeqIdPtr       tsip;
  ValNodePtr         vnp;
  /*
  SeqMgrDescContext  dcontext;
  PubdescPtr         pdp;
  SeqDescrPtr        sdp;
  */

  if (sfp == NULL || fcontext == NULL) return FALSE;
  awp = (Asn2gbWorkPtr) fcontext->userdata;
  if (awp == NULL) return FALSE;
  ajp = awp->ajp;
  if (ajp == NULL) return FALSE;
  asp = awp->asp;
  if (asp == NULL) return FALSE;
  bsp = asp->bsp;
  if (bsp == NULL) return FALSE;

  if (fcontext->featdeftype == FEATDEF_PUB ||
      fcontext->featdeftype == FEATDEF_NON_STD_RESIDUE ||
      fcontext->featdeftype == FEATDEF_RSITE ||
      fcontext->featdeftype == FEATDEF_SEQ) return TRUE;

  if (fcontext->featdeftype == FEATDEF_BIOSRC) return TRUE;

  if (ajp->flags.validateFeats &&
      (fcontext->featdeftype == FEATDEF_BAD ||
       fcontext->featdeftype == FEATDEF_virion)) {
    return TRUE;
  }

  if (ISA_na (bsp->mol) && fcontext->featdeftype == FEATDEF_HET) return TRUE;

  /* check feature customization flags */

  if (awp->hideImpFeats && sfp->data.choice == SEQFEAT_IMP && fcontext->featdeftype != FEATDEF_operon) return TRUE;
  if (awp->hideVariations && fcontext->featdeftype == FEATDEF_variation) return TRUE;
  if (awp->hideRepeatRegions && fcontext->featdeftype == FEATDEF_repeat_region) return TRUE;
  if (awp->hideRepeatRegions && fcontext->featdeftype == FEATDEF_mobile_element) return TRUE;
  if (awp->hideGaps && fcontext->featdeftype == FEATDEF_gap) return TRUE;
  if (ISA_aa (bsp->mol) && fcontext->featdeftype == FEATDEF_REGION &&
      awp->hideCddFeats && IsCDD (sfp)) return TRUE;
  if (awp->hideSitesBondsRegions && (fcontext->featdeftype == FEATDEF_REGION ||
                                     fcontext->featdeftype == FEATDEF_SITE ||
                                     fcontext->featdeftype == FEATDEF_BOND)) return TRUE;

  /* DDBJ does not want to show gene features */

  if (fcontext->seqfeattype == SEQFEAT_GENE && awp->hideGeneFeats) return TRUE;

  /* no longer suppressing comment features that are full length */

  /*
  if (fcontext->seqfeattype == SEQFEAT_COMMENT &&
      fcontext->left == awp->from && fcontext->right == awp->to) return TRUE;
  */

  /*
  if (ISA_aa (bsp->mol) && awp->format == GENPEPT_FMT && fcontext->seqfeattype == SEQFEAT_PROT) {
    if (fcontext->left == awp->from && fcontext->right == awp->to) {
      if (awp->bestprot != sfp) {
        if (EquivProtFeats (awp->bestprot, sfp)) return TRUE;
      }
    }
  }
  */

  ivals = fcontext->ivals;
  numivals = fcontext->numivals;

  /* check to see if last interval is on this awp->from - awp->to range */

  if (ivals != NULL && numivals > 0) {
    idx = (numivals - 1) * 2;
    start = ivals [idx];
    stop = ivals [idx + 1];
    if (stop < awp->from || stop > awp->to) {

      /* may need to map sig_peptide on a different segment */

      if (fcontext->seqfeattype == SEQFEAT_CDREGION) {
        sip = SeqLocIdForProduct (sfp->product);
        bsp = BioseqFind (sip);
        GetFeatsOnCdsProduct (sfp, asp->bsp, bsp, ajp, awp);
      }

      if (! awp->showAllFeats) return TRUE;

      /* if showing one segment, only show features covering this segment */

      if (fcontext->right < awp->from || fcontext->left > awp->to) return TRUE;

    } else if (fcontext->farloc && NotEMBLorDDBJ (awp->bsp)) {

      /* last interval may not have been mapped to bioseq if far */

      firstslp = NULL;
      lastslp = NULL;

      slp = SeqLocFindNext (sfp->location, NULL);
      while (slp != NULL) {
        if (slp->choice != SEQLOC_NULL) {
          lastslp = slp;
          if (firstslp == NULL) {
            firstslp = slp;
          }
        }
        slp = SeqLocFindNext (sfp->location, slp);
      }

      /* !!! EMBL may have different desired behavior on where to map !!! */

      if (firstslp != NULL && SeqLocStrand (firstslp) == Seq_strand_minus) {
        slp = firstslp;
      } else {
        slp = lastslp;
      }

      if (slp != NULL) {
        sip = SeqLocId (slp);
        if (sip != NULL) {
          bsp = BioseqFindCore (sip);
          if (bsp == NULL || (bsp != awp->parent && bsp != awp->bsp)) {

            return TRUE;
          }
        }
      }
    }
  }

  /* make sure feature is within sublocation */

  if (ajp->ajp.slp != NULL) {
    if (SeqLocCompare (sfp->location, ajp->ajp.slp) == SLC_NO_MATCH) {
      slp = SeqLocMerge (bsp, sfp->location, NULL, FALSE, TRUE, FALSE);
      if (slp == NULL) return TRUE;
      sip = SeqIdParse ("lcl|dummy");
      newloc = SeqLocReMapEx (sip, ajp->ajp.slp, slp, 0, FALSE, ajp->masterStyle, ajp->relaxedMapping);
      SeqIdFree (sip);
      SeqLocFree (slp);
      if (newloc == NULL) return TRUE;
      SeqLocFree (newloc);
    }
  }

  /* suppress duplicate features (on nucleotide) */

  lastsfp = awp->lastsfp;
  lastsap = awp->lastsap;
  if (lastsfp != NULL && lastsap != NULL) {
    if (lastsfp->idx.subtype == sfp->idx.subtype &&
        awp->lastleft == fcontext->left &&
        awp->lastright == fcontext->right) {
        if (lastsap == fcontext->sap ||
            (lastsap->desc == NULL && fcontext->sap->desc == NULL)) {
        if (AsnIoMemComp (lastsfp, sfp, (AsnWriteFunc) SeqFeatAsnWrite)) {
          return TRUE;
        }
      }
    }
  }

  /* if RELEASE_MODE, verify that features have all mandatory qualifiers */

  if (ajp->flags.needRequiredQuals) {
    okay = FALSE;

    switch (fcontext->featdeftype) {

    case FEATDEF_CDS:
      if (ajp->flags.checkCDSproductID) {
        /* non-pseudo CDS must have /product */
        if (sfp->pseudo) {
          pseudo = TRUE;
        }
        /*
        grp = SeqMgrGetGeneXref (sfp);
        */
        grp = GetGeneByFeat (sfp, &psdo, &supr);
        if (psdo) {
          pseudo = TRUE;
        }
        /*
        if (grp == NULL) {
          sep = GetTopSeqEntryForEntityID (ajp->ajp.entityID);
          oldscope = SeqEntrySetScope (sep);
          gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
          SeqEntrySetScope (oldscope);
          if (gene != NULL) {
            grp = (GeneRefPtr) gene->data.value.ptrvalue;
            if (gene->pseudo) {
              pseudo = TRUE;
            }
          }
        }
        */
        if (grp != NULL && grp->pseudo) {
          pseudo = TRUE;
        }
        if (sfp->location != NULL) {
          if (CheckSeqLocForPartial (sfp->location, &partial5, &partial3)) {
            if (partial5 && (! partial3)) {
              if (SeqLocLen (sfp->location) <= 5) {
                juststop = TRUE;
              }
            }
          }
        }
        if (pseudo || juststop) {
          okay = TRUE;
        } else if (sfp->product != NULL) {
          sip = SeqLocIdForProduct (sfp->product);
          if (sip != NULL) {
            if ((sip->choice == SEQID_GI && sip->data.intvalue > 0) ||
                sip->choice == SEQID_LOCAL) {
              sep = GetTopSeqEntryForEntityID (ajp->ajp.entityID);
              oldscope = SeqEntrySetScope (sep);
              prod = BioseqFind (sip);
              SeqEntrySetScope (oldscope);
              if (prod != NULL) {
                for (sip = prod->id; sip != NULL; sip = sip->next) {
                  if (sip->choice == SEQID_GENBANK ||
                     sip->choice == SEQID_EMBL ||
                      sip->choice == SEQID_DDBJ ||
                      sip->choice == SEQID_OTHER ||
                      sip->choice == SEQID_PATENT ||
                      sip->choice == SEQID_TPG ||
                      sip->choice == SEQID_TPE ||
                      sip->choice == SEQID_TPD ||
                      sip->choice == SEQID_GPIPE) {
                    tsip = (TextSeqIdPtr) sip->data.ptrvalue;
                    if (tsip != NULL && (StringDoesHaveText (tsip->accession))) {
                      if (ValidateAccn (tsip->accession) == 0)
                      okay = TRUE;
                    }
                  }
                }
              } else if (sip->choice == SEQID_GI && sip->data.intvalue > 0) {
                /* RELEASE_MODE requires that /protein_id is an accession */
                gi = sip->data.intvalue;
                if (GetAccnVerFromServer (gi, buf)) {
                  okay = TRUE;
                } else {
                  sip = GetSeqIdForGI (gi);
                  if (sip != NULL) {
                    okay = TRUE;
                  }
                }
              }
            } else if (sip->choice == SEQID_GENBANK ||
                       sip->choice == SEQID_EMBL ||
                       sip->choice == SEQID_DDBJ ||
                       sip->choice == SEQID_OTHER ||
                       sip->choice == SEQID_PATENT ||
                       sip->choice == SEQID_TPG ||
                       sip->choice == SEQID_TPE ||
                       sip->choice == SEQID_TPD ||
                       sip->choice == SEQID_GPIPE) {
              tsip = (TextSeqIdPtr) sip->data.ptrvalue;
              if (tsip != NULL && (StringDoesHaveText (tsip->accession))) {
                if (ValidateAccn (tsip->accession) == 0)
                okay = TRUE;
              }
            }
          }
        } else {
          if (sfp->excpt && (StringDoesHaveText (sfp->except_text))) {
            if (StringStr (sfp->except_text, "rearrangement required for product") != NULL) {
              okay = TRUE;
            }
          }
        }
      } else {
        okay = TRUE;
      }
      if (! okay) {
        ajp->relModeError = TRUE;
      }
      break;

    case FEATDEF_conflict:
      if (sfp->cit == NULL) {
        /* RefSeq allows conflict with accession in comment instead of sfp->cit */
        for (sip = bsp->id; sip != NULL; sip = sip->next) {
          if (sip->choice == SEQID_OTHER) {
            if (StringDoesHaveText (sfp->comment)) {
              okay = TRUE;
            }
          }
        }
      }
      /* continue on to old_sequence */
    case FEATDEF_old_sequence:
      /* conflict and old_sequence require a publication printable on the segment */
      vnp = sfp->cit;

      if (vnp != NULL && asp->referenceArray != NULL) {
        for (ppr = vnp->data.ptrvalue; ppr != NULL; ppr = ppr->next) {
          j = MatchRef (ppr, asp->referenceArray, asp->numReferences);
          if (j > 0) {
            okay = TRUE;
            break;
          }
        }
      }
      if (! okay) {
        /* compare qualifier can now substitute for citation qualifier */
        gbq = sfp->qual;
        while (gbq != NULL) {
          if (StringICmp (gbq->qual, "compare") == 0 && (StringDoesHaveText (gbq->val))) {
            okay = TRUE;
            break;
          }
          gbq = gbq->next;
        }
      }
      break;

    case FEATDEF_GENE:
      /* gene requires /gene or /locus_tag, but desc or syn can be mapped to /gene */
      grp = (GeneRefPtr) sfp->data.value.ptrvalue;
      if (grp != NULL) {
        if (StringDoesHaveText (grp->locus)) {
          okay = TRUE;
        }  else if (StringDoesHaveText (grp->locus_tag)) {
          okay = TRUE;
        } else if (StringDoesHaveText (grp->desc)) {
          okay = TRUE;
        } else {
          vnp = grp->syn;
          if (vnp != NULL) {
            if (StringDoesHaveText (vnp->data.ptrvalue)) {
              okay = TRUE;
            }
          }
        }
      }
      break;

    case FEATDEF_protein_bind:
    case FEATDEF_misc_binding:
      /* protein_bind or misc_binding require FTQUAL_bound_moiety */
      gbq = sfp->qual;
      while (gbq != NULL) {
        if (StringICmp (gbq->qual, "bound_moiety") == 0 && (StringDoesHaveText (gbq->val))) {
          okay = TRUE;
          break;
        }
        gbq = gbq->next;
      }
      break;

    case FEATDEF_modified_base:
      /* modified_base requires FTQUAL_mod_base */
      gbq = sfp->qual;
      while (gbq != NULL) {
        if (StringICmp (gbq->qual, "mod_base") == 0 && (StringDoesHaveText (gbq->val))) {
          okay = TRUE;
          break;
        }
        gbq = gbq->next;
      }
      break;

    case FEATDEF_gap:
      /* gap requires FTQUAL_estimated_length */
      gbq = sfp->qual;
      while (gbq != NULL) {
        if (StringICmp (gbq->qual, "estimated_length") == 0 && (StringDoesHaveText (gbq->val))) {
          okay = TRUE;
          break;
        }
        gbq = gbq->next;
      }
      break;

    case FEATDEF_operon:
      /* operon requires FTQUAL_operon */
      gbq = sfp->qual;
      while (gbq != NULL) {
        if (StringICmp (gbq->qual, "operon") == 0 && (StringDoesHaveText (gbq->val))) {
          okay = TRUE;
          break;
        }
        gbq = gbq->next;
      }
      break;

    case FEATDEF_ncRNA:
      /* ncRNA requires FTQUAL_ncRNA_class */
      gbq = sfp->qual;
      while (gbq != NULL) {
        if (StringICmp (gbq->qual, "ncRNA_class") == 0 && (StringDoesHaveText (gbq->val))) {
          okay = TRUE;
          break;
        }
        gbq = gbq->next;
      }
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp != NULL && rrp->ext.choice == 3) {
        rgp = (RNAGenPtr) rrp->ext.value.ptrvalue;
        if (rgp != NULL) {
          if (StringDoesHaveText (rgp->_class)) {
            okay = TRUE;
            break;
          }
        }
      }
      break;

    case FEATDEF_mobile_element:
      /* mobile_element requires FTQUAL_mobile_element_type */
      gbq = sfp->qual;
      while (gbq != NULL) {
        if (StringICmp (gbq->qual, "mobile_element_type") == 0 && (StringDoesHaveText (gbq->val))) {
          okay = TRUE;
          break;
        }
        gbq = gbq->next;
      }
      break;

    case FEATDEF_assembly_gap:
      /* assembly_gap requires FTQUAL_estimated_length and FTQUAL_gap_type */
      has_est_len = FALSE;
      has_gap_type = FALSE;
      gbq = sfp->qual;
      while (gbq != NULL) {
        if (StringDoesHaveText (gbq->val)) {
          if (StringICmp (gbq->qual, "estimated_length") == 0) {
            has_est_len = TRUE;
          } else if (StringICmp (gbq->qual, "gap_type") == 0) {
            has_gap_type = TRUE;
          }
        }
        gbq = gbq->next;
      }
      if (has_est_len && has_gap_type) {
        okay = TRUE;
      }
      break;

    case FEATDEF_regulatory:
      /* regulatory requires FTQUAL_regulatory_class */
      gbq = sfp->qual;
      while (gbq != NULL) {
        if (StringICmp (gbq->qual, "regulatory_class") == 0 && (StringDoesHaveText (gbq->val))) {
          okay = TRUE;
          break;
        }
        gbq = gbq->next;
      }
      break;

    default:
      if (fcontext->featdeftype >= FEATDEF_GENE && fcontext->featdeftype < FEATDEF_MAX) {
        okay = TRUE;
      }
      break;
    }

    if (okay == FALSE) return TRUE;
  }

  /* if RELEASE_MODE, suppress features with location on near segmented Bioseq */

  if (ajp->flags.suppressSegLoc) {
    bsp = awp->parent;
    if (bsp != NULL && bsp->repr == Seq_repr_seg && SegHasParts (bsp)) {
      slp = SeqLocFindNext (sfp->location, NULL);
      while (slp != NULL) {
        sip = SeqLocId (slp);
        if (sip != NULL) {
          if (SeqIdIn (sip, bsp->id)) return TRUE;
        }
        slp = SeqLocFindNext (sfp->location, slp);
      }
    }
  }

  gap = awp->currfargap;
  if (gap != NULL && awp->afp != NULL) {
    while (gap != NULL && LocInBioseq (gap->location, asp->bsp) && GetOffsetInBioseq (gap->location, asp->bsp, SEQLOC_LEFT_END) < fcontext->left) {

      fbp = (FeatBlockPtr) Asn2gbAddBlock (awp, FEATURE_BLOCK, sizeof (IntFeatBlock));
      if (fbp != NULL) {
        fbp->entityID = 0;
        fbp->itemID = 0;
        fbp->itemtype = OBJ_SEQFEAT;
        fbp->featdeftype = FEATDEF_gap;
        ifp = (IntFeatBlockPtr) fbp;
        ifp->mapToNuc = FALSE;
        ifp->mapToProt = FALSE;
        ifp->mapToGen = FALSE;
        ifp->mapToMrna = FALSE;
        ifp->mapToPep = FALSE;
        ifp->left = 0;
        ifp->right = 0;
        if (bsp != NULL) {
          SetIfpFeatCount (ifp, ajp, awp, ISA_aa (bsp->mol));
        }
        ifp->firstfeat = awp->firstfeat;
        awp->firstfeat = FALSE;
        if (awp->afp != NULL) {
          DoImmediateRemoteFeatureFormat (awp->afp, (BaseBlockPtr) fbp, gap);
        }
      }

      awp->currfargap = gap->next;
      gap = awp->currfargap;
    }
  }

  /* check for Imp-feat gap that is same as next Seq-lit gap - but need to check against scaffold coordinate */
  if (! NotEMBLorDDBJ (awp->bsp)) {
    if (gap != NULL && LocInBioseq (gap->location, asp->bsp) && fcontext->featdeftype == FEATDEF_gap &&
        GetOffsetInBioseq (gap->location, asp->bsp, SEQLOC_LEFT_END) == fcontext->left &&
        GetOffsetInBioseq (gap->location, asp->bsp, SEQLOC_RIGHT_END) == fcontext->right) {
      awp->currfargap = gap->next;
    }
  }

  awp->lastsfp = sfp;
  awp->lastsap = fcontext->sap;
  awp->lastleft = fcontext->left;
  awp->lastright = fcontext->right;

  if (fcontext->seqfeattype == SEQFEAT_CDREGION) {
    fbp = (FeatBlockPtr) Asn2gbAddBlock (awp, FEATURE_BLOCK, sizeof (IntCdsBlock));
  } else if (fcontext->seqfeattype == SEQFEAT_PROT) {
    fbp = (FeatBlockPtr) Asn2gbAddBlock (awp, FEATURE_BLOCK, sizeof (IntPrtBlock));
  } else {
    fbp = (FeatBlockPtr) Asn2gbAddBlock (awp, FEATURE_BLOCK, sizeof (IntFeatBlock));
  }
  if (fbp == NULL) return TRUE;

  fbp->entityID = fcontext->entityID;
  fbp->itemID = fcontext->itemID;
  fbp->itemtype = OBJ_SEQFEAT;
  fbp->featdeftype = fcontext->featdeftype;
  ifp = (IntFeatBlockPtr) fbp;
  ifp->mapToNuc = FALSE;
  ifp->mapToProt = FALSE;
  ifp->mapToGen = FALSE;
  ifp->mapToMrna = FALSE;
  ifp->mapToPep = FALSE;
  ifp->left = 0;
  ifp->right = 0;
  if (bsp != NULL) {
    SetIfpFeatCount (ifp, ajp, awp, ISA_aa (bsp->mol));
  }
  ifp->firstfeat = awp->firstfeat;
  awp->firstfeat = FALSE;

  /* local centromere, telomere, rep_origin, and region features (e.g, on eukaryotic NC record) do not contribute to test for far fetch suppression */
  if (sfp->idx.subtype != FEATDEF_centromere &&
      sfp->idx.subtype != FEATDEF_telomere &&
      sfp->idx.subtype != FEATDEF_rep_origin &&
      sfp->idx.subtype != FEATDEF_REGION) {

    /* this allows remote SNP, CDD, MGC, etc., not to be treated as local annotation */
    if (awp->entityID != fbp->entityID || fbp->itemID <= awp->localFeatCount) {
      awp->featseen = TRUE;
    }
    awp->featjustseen = TRUE;
  }

  if (fcontext->seqfeattype == SEQFEAT_PROT) {

    /* set calculated molecular weight flags for proteins */

    ifp->isPrt = TRUE;
    ipp = (IntPrtBlockPtr) fbp;
    prp = (ProtRefPtr) sfp->data.value.ptrvalue;
    if (prp != NULL) {
      if (prp->processed < 2) {
        is_whole = FALSE;
        slp = sfp->location;
        if (slp != NULL) {
          if (slp->choice == SEQLOC_WHOLE) {
            is_whole = TRUE;
          } else if (slp->choice == SEQLOC_INT) {
            sintp = (SeqIntPtr) slp->data.ptrvalue;
            if (sintp != NULL && 
                bsp != NULL &&
                sintp->from == 0 &&
                sintp->to == bsp->length - 1) {
              is_whole = TRUE;
            }
          }
        }
        if (is_whole) {
          ipp->is_whole_loc = TRUE;
          if (awp->has_sig_peptide) {
            if (awp->has_mat_peptide) {
              ipp->suppress_mol_wt = TRUE;
            } else if (awp->sig_pept_trim_len > 0) {
              ipp->sig_pept_trim_len = awp->sig_pept_trim_len;
            }
          } else {
            ipp->trim_initial_met = TRUE;
          }
        }
      }
    }
  }

  if (awp->afp != NULL) {
    DoImmediateFormat (awp->afp, (BaseBlockPtr) fbp);
  }

  /* optionally map CDS from cDNA onto genomic */

  if (awp->isGPS && bsp != NULL && ISA_na (bsp->mol) && awp->copyGpsCdsUp &&
      fcontext->featdeftype == FEATDEF_mRNA) {
    sip = SeqLocIdForProduct (sfp->product);
    bsp = BioseqFind (sip);
    if (bsp != NULL && ISA_na (bsp->mol)) {
      cds = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_CDREGION, 0, &cdscontext);
      if (cds != NULL) {
        fbp = (FeatBlockPtr) Asn2gbAddBlock (awp, FEATURE_BLOCK, sizeof (IntCdsBlock));
        if (fbp != NULL) {

          fbp->entityID = cdscontext.entityID;
          fbp->itemID = cdscontext.itemID;
          fbp->itemtype = OBJ_SEQFEAT;
          fbp->featdeftype = cdscontext.featdeftype;
          ifp = (IntFeatBlockPtr) fbp;
          ifp->mapToNuc = FALSE;
          ifp->mapToProt = FALSE;
          ifp->mapToGen = TRUE;
          ifp->mapToMrna = FALSE;
          ifp->mapToPep = FALSE;
          ifp->left = 0;
          ifp->right = 0;
          SetIfpFeatCount (ifp, ajp, awp, FALSE);
          ifp->firstfeat = awp->firstfeat;
          awp->firstfeat = FALSE;

          if (awp->afp != NULL) {
            DoImmediateFormat (awp->afp, (BaseBlockPtr) fbp);
          }
        }
      }
    }
  }

  if (fcontext->seqfeattype != SEQFEAT_CDREGION) return TRUE;

  /* if feature table format, do not get features from protein product */

  if (awp->format == FTABLE_FMT) return TRUE;

  /* if CDS, collect more information from product protein bioseq - may be part */

  sip = SeqLocIdForProduct (sfp->product);
  bsp = BioseqFind (sip);
  if (bsp == NULL || (! ISA_aa (bsp->mol))) return TRUE;

  ifp->isCDS = TRUE;
  icp = (IntCdsBlockPtr) ifp;

  /* first explore pubs to pick up figure and maploc - no longer shown */

  /*
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_pub, &dcontext);
  while (sdp != NULL) {
    pdp = (PubdescPtr) sdp->data.ptrvalue;
    if (pdp != NULL) {
      if (icp->fig == NULL) {
        icp->fig = StringSaveNoNull (pdp->fig);
      }
      if (icp->maploc == NULL) {
        icp->maploc = StringSaveNoNull (pdp->maploc);
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_pub, &dcontext);
  }
  */

  /* product may be segmented part, and remaining features are indexed on parent */

  parent = SeqMgrGetParentOfPart (bsp, NULL);
  if (parent != NULL) {
    bsp = parent;
  }

  /* then explore mat_peptides, sites, etc. */

  GetFeatsOnCdsProduct (sfp, asp->bsp, bsp, ajp, awp);

  GetRemoteFeatsOnCdsProduct (sfp, asp->bsp, bsp, ajp, awp);

  return TRUE;
}

/*
static Boolean TestGetAccnVerFromServer (BIG_ID gi, CharPtr buf)

{
  Char      accn [64];
  SeqIdPtr  sip;

  if (buf == NULL) return FALSE;
  *buf = '\0';
  sip = GetSeqIdForGI (gi);
  if (sip == NULL) return FALSE;
  SeqIdWrite (sip, accn, PRINTID_TEXTID_ACC_VER, sizeof (accn) - 1);
  SeqIdFree (sip);
  if (StringLen (accn) < 40) {
    StringCpy (buf, accn);
  }
  return TRUE;
}
*/

static WgsAccnPtr GetWgsNode (
  Asn2gbWorkPtr awp,
  CharPtr accn
)

{
  ValNodePtr  vnp;
  WgsAccnPtr  wap = NULL;

  if (awp == NULL || StringHasNoText (accn)) return NULL;

  for (vnp = awp->wgsaccnlist; vnp != NULL; vnp = vnp->next) {
    wap = (WgsAccnPtr) vnp->data.ptrvalue;
    if (wap == NULL) continue;
    if (StringCmp (accn, wap->accn) == 0) return wap;
  }
  wap = (WgsAccnPtr) MemNew (sizeof (WgsAccn));
  if (wap == NULL) return NULL;
  StringCpy (wap->accn, accn);
  ValNodeAddPointer (&(awp->wgsaccnlist), 0, (Pointer) wap);
  return wap;
}

static Boolean LIBCALLBACK GetFeatsOnSeg (
  SeqLocPtr slp,
  SeqMgrSegmentContextPtr context
)

{
  Char             accn [41];
  Uint4            accntype;
  IntAsn2gbJobPtr  ajp;
  Asn2gbWorkPtr    awp;
  BioseqPtr        bsp;
  time_t           currTime;
  Uint2            entityID;
  Int4             from;
  BIG_ID           gi;
  Int4             left;
  SeqLocPtr        loc;
  CharPtr          ptr;
  Int4             right;
  SeqIdPtr         sip;
  Int4             to;
  WgsAccnPtr       wap = NULL;

  if (slp == NULL || context == NULL) return FALSE;
  awp = (Asn2gbWorkPtr) context->userdata;
  if (awp == NULL) return FALSE;
  ajp = awp->ajp;
  if (ajp == NULL) return FALSE;

  /* do not fetch outside of desired component */

  if (ajp->ajp.slp != NULL) {
    left = GetOffsetInBioseq (ajp->ajp.slp, awp->parent, SEQLOC_LEFT_END);
    right = GetOffsetInBioseq (ajp->ajp.slp, awp->parent, SEQLOC_RIGHT_END);

    from = context->cumOffset;
    to = from + context->to - context->from;

    if (left > to) return TRUE;
    if (right < from) return TRUE;
  }
 
  from = awp->from;
  to = awp->to;

  sip = SeqLocId (slp);
  if (sip == NULL) {
    loc = SeqLocFindNext (slp, NULL);
    if (loc != NULL) {
      sip = SeqLocId (loc);
    }
  }
  if (sip == NULL) return TRUE;

  /* if Web Entrez WGS */

  if (awp->farFeatTimeLimit) {
    if (sip->choice == SEQID_GI) {
      gi = (BIG_ID) sip->data.intvalue;
      if (GetAccnVerFromServer (gi, accn)) {
        ptr = StringChr (accn, '.');
        if (ptr != NULL) {
          *ptr = '\0';
        }
        accntype = WHICH_db_accession (accn);
        if (ACCN_IS_WGS (accntype)) {
          accn [4] = '\0';
          wap = GetWgsNode (awp, accn);
          if (wap != NULL) {
            (wap->count)++;
            if (wap->count > 50) {
              if (! wap->hasfeats) return TRUE;
            }
          }
        }
      }
    }
    if (! awp->featseen) {
      currTime = GetSecs ();
      if (currTime - awp->farFeatStartTime > 25) return FALSE;
    }
  }

  /* may want to remote fetch genome component if not already in memory */

  bsp = BioseqLockById (sip);

  if (bsp == NULL) return TRUE;

  entityID = ObjMgrGetEntityIDForPointer (bsp);

  if (entityID != awp->entityID) {

    /* if segment not packaged in record, may need to feature index it */

    if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
      SeqMgrIndexFeatures (entityID, NULL);
    }

    /* collect features indexed on the remote bioseq */

    awp->from = 0;
    awp->to = bsp->length - 1;
  }

  awp->lastsfp = NULL;
  awp->lastsap = NULL;
  awp->lastleft = 0;
  awp->lastright = 0;

  awp->featjustseen = FALSE;

  if (context->strand == Seq_strand_minus) {
    SeqMgrExploreFeaturesRev (bsp, (Pointer) awp, GetFeatsOnBioseq, /* awp->slp */ slp, NULL, NULL);
  } else {
    SeqMgrExploreFeatures (bsp, (Pointer) awp, GetFeatsOnBioseq, /* awp->slp */ slp, NULL, NULL);
  }

  if (awp->featjustseen && wap != NULL) {
    wap->hasfeats = TRUE;
  }

  /* restore original from and to */

  awp->from = from;
  awp->to = to;

  BioseqUnlock (bsp);

  return TRUE;
}

NLM_EXTERN void AddFeatureBlock (
  Asn2gbWorkPtr awp
)

{
  IntAsn2gbJobPtr    ajp;
  BioseqPtr          bsp;
  SeqFeatPtr         cds;
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext;
  FeatBlockPtr       fbp;
  SeqFeatPtr         gene;
  IntFeatBlockPtr    ifp;
  Boolean            is_other;
  MolInfoPtr         mip;
  SeqFeatPtr         mrna;
  SeqMgrFeatContext  pcontext;
  SeqFeatPtr         prot;
  SeqDescrPtr        sdp;
  SeqIdPtr           sip;
  SeqLocPtr          slp;

  if (awp == NULL) return;
  ajp = awp->ajp;
  if (ajp == NULL) return;
  bsp = awp->parent;
  if (bsp == NULL) return;

  awp->lastsfp = NULL;
  awp->lastsap = NULL;
  awp->lastleft = 0;
  awp->lastright = 0;

  /* for protein molecular weight calculation, need sig_peptide, etc. */

  awp->has_mat_peptide = FALSE;
  awp->has_sig_peptide = FALSE;
  awp->sig_pept_trim_len = 0;

  if (awp->format == GENPEPT_FMT && ISA_aa (bsp->mol)) {
    awp->bestprot = SeqMgrGetBestProteinFeature (bsp, &pcontext);

    prot = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
    while (prot != NULL) {
      if (fcontext.featdeftype == FEATDEF_sig_peptide_aa ||
          fcontext.featdeftype == FEATDEF_transit_peptide_aa) {
        awp->has_sig_peptide = TRUE;
        if (fcontext.left == 0 && fcontext.right < bsp->length - 1) {
          awp->sig_pept_trim_len = fcontext.right + 1;
        }
      } else if (fcontext.featdeftype == FEATDEF_mat_peptide_aa) {
        awp->has_mat_peptide = TRUE;
      }

      prot = SeqMgrGetNextFeature (bsp, prot, 0, 0, &fcontext);
    }
  }

  /* optionally map gene from genomic onto cDNA */

  if (awp->isGPS && ISA_na (bsp->mol) && awp->copyGpsGeneDown) {
    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
    if (sdp != NULL && sdp->choice == Seq_descr_molinfo) {
      mip = (MolInfoPtr) sdp->data.ptrvalue;
      if (mip != NULL) {
        if (mip->biomol == MOLECULE_TYPE_MRNA) {
          mrna = SeqMgrGetRNAgivenProduct (bsp, NULL);
          if (mrna != NULL) {
            gene = SeqMgrGetOverlappingGene (mrna->location, &fcontext);
            if (gene != NULL && gene->data.choice == SEQFEAT_GENE) {

              fbp = (FeatBlockPtr) Asn2gbAddBlock (awp, FEATURE_BLOCK, sizeof (IntCdsBlock));
              if (fbp != NULL) {

                fbp->entityID = fcontext.entityID;
                fbp->itemID = fcontext.itemID;
                fbp->itemtype = OBJ_SEQFEAT;
                fbp->featdeftype = fcontext.featdeftype;
                ifp = (IntFeatBlockPtr) fbp;
                ifp->mapToNuc = FALSE;
                ifp->mapToProt = FALSE;
                ifp->mapToGen = FALSE;
                ifp->mapToMrna = TRUE;
                ifp->mapToPep = FALSE;
                ifp->isCDS = TRUE;
                ifp->left = 0;
                ifp->right = 0;
                SetIfpFeatCount (ifp, ajp, awp, FALSE);
                ifp->firstfeat = awp->firstfeat;
                awp->firstfeat = FALSE;

                if (awp->afp != NULL) {
                  DoImmediateFormat (awp->afp, (BaseBlockPtr) fbp);
                }
              }
            }
          }
        }
      }
    }
  }

  awp->farFeatTimeLimit = FALSE;
  if (bsp->repr == Seq_repr_seg || bsp->repr == Seq_repr_delta || bsp->repr == Seq_repr_ref) {
    if (awp->mode == ENTREZ_MODE) {
      awp->farFeatTimeLimit = TRUE;
    }
    /*
    if (GetWWW (ajp) && awp->mode == ENTREZ_MODE) {
      sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
      if (sdp != NULL && sdp->choice == Seq_descr_molinfo && sdp->data.ptrvalue != NULL) {
        mip = (MolInfoPtr) sdp->data.ptrvalue;
        if (mip->tech == MI_TECH_wgs || mip->tech == MI_TECH_composite_wgs_htgs) {
          awp->farFeatTimeLimit = TRUE;
        }
      }
    }
    */
  }

  if (! awp->onlyNearFeats) {
    if (awp->farFeatsSuppress) {

      if (bsp->repr == Seq_repr_seg || bsp->repr == Seq_repr_delta || bsp->repr == Seq_repr_ref) {

        /* get start time for 25 second timeout in Web Entrez far WGS records */

        if (awp->farFeatTimeLimit) {
          awp->farFeatStartTime = GetSecs ();
        }

        /* if farFeatsSuppress first collect features on remote segments in MASTER_STYLE */

        SeqMgrExploreSegments (bsp, (Pointer) awp, GetFeatsOnSeg);

        awp->wgsaccnlist = ValNodeFreeData (awp->wgsaccnlist);
      }
    }
  }

  if ((! awp->farFeatsSuppress) || (! awp->featseen)) {

    /* reminder - features on near parts are indexed on segmented Bioseq */

    slp = ajp->ajp.slp;
    if (slp != NULL && SeqLocStrand (slp) == Seq_strand_minus) {
      SeqMgrExploreFeaturesRev (bsp, (Pointer) awp, GetFeatsOnBioseq, awp->slp, NULL, NULL);
    } else {
      SeqMgrExploreFeatures (bsp, (Pointer) awp, GetFeatsOnBioseq, awp->slp, NULL, NULL);
    }
  }


  if (awp->format == GENPEPT_FMT && ISA_aa (bsp->mol)) {
    cds = SeqMgrGetCDSgivenProduct (bsp, &fcontext);
    if (cds != NULL && cds->data.choice == SEQFEAT_CDREGION) {

      if (fcontext.entityID > 0 && fcontext.itemID > 0) {

        fbp = (FeatBlockPtr) Asn2gbAddBlock (awp, FEATURE_BLOCK, sizeof (IntCdsBlock));
        if (fbp != NULL) {

          fbp->entityID = fcontext.entityID;
          fbp->itemID = fcontext.itemID;
          fbp->itemtype = OBJ_SEQFEAT;
          fbp->featdeftype = fcontext.featdeftype;
          ifp = (IntFeatBlockPtr) fbp;
          ifp->mapToNuc = FALSE;
          ifp->mapToProt = TRUE;
          ifp->mapToGen = FALSE;
          ifp->mapToMrna = FALSE;
          ifp->mapToPep = FALSE;
          ifp->isCDS = TRUE;
          ifp->left = 0;
          ifp->right = 0;
          SetIfpFeatCount (ifp, ajp, awp, TRUE);
          ifp->firstfeat = awp->firstfeat;
          awp->firstfeat = FALSE;

          if (awp->afp != NULL) {
            DoImmediateFormat (awp->afp, (BaseBlockPtr) fbp);
          }
        }
      } else if (cds->idx.entityID > 0 && cds->idx.itemID > 0) {

        /* if protein bioseq and cds feature but no nucleotide, handle as special case */

        fbp = (FeatBlockPtr) Asn2gbAddBlock (awp, FEATURE_BLOCK, sizeof (IntCdsBlock));
        if (fbp != NULL) {

          fbp->entityID = cds->idx.entityID;
          fbp->itemID = cds->idx.itemID;
          fbp->itemtype = OBJ_SEQFEAT;
          fbp->featdeftype = FEATDEF_CDS;
          ifp = (IntFeatBlockPtr) fbp;
          ifp->mapToNuc = FALSE;
          ifp->mapToProt = TRUE;
          ifp->mapToGen = FALSE;
          ifp->mapToMrna = FALSE;
          ifp->mapToPep = FALSE;
          ifp->isCDS = TRUE;
          ifp->left = 0;
          ifp->right = 0;
          SetIfpFeatCount (ifp, ajp, awp, TRUE);
          ifp->firstfeat = awp->firstfeat;
          awp->firstfeat = FALSE;

          if (awp->afp != NULL) {
            DoImmediateFormat (awp->afp, (BaseBlockPtr) fbp);
          }
        }
      }
    }
    prot = SeqMgrGetPROTgivenProduct (bsp, &fcontext);
    if (prot != NULL && prot->data.choice == SEQFEAT_PROT) {

      is_other = FALSE;
      for (sip = bsp->id; sip != NULL; sip = sip->next) {
        if (sip->choice == SEQID_OTHER) {
          is_other = TRUE;
        }
      }

      /* for RefSeq records or GenBank not release_mode */
      if (is_other || (! ajp->flags.forGbRelease)) {

        fbp = (FeatBlockPtr) Asn2gbAddBlock (awp, FEATURE_BLOCK, sizeof (IntCdsBlock));
        if (fbp != NULL) {

          fbp->entityID = fcontext.entityID;
          fbp->itemID = fcontext.itemID;
          fbp->itemtype = OBJ_SEQFEAT;
          fbp->featdeftype = fcontext.featdeftype;
          ifp = (IntFeatBlockPtr) fbp;
          ifp->mapToNuc = FALSE;
          ifp->mapToProt = FALSE;
          ifp->mapToGen = FALSE;
          ifp->mapToMrna = FALSE;
          ifp->mapToPep = TRUE;
          ifp->left = 0;
          ifp->right = 0;
          SetIfpFeatCount (ifp, ajp, awp, TRUE);
          ifp->firstfeat = awp->firstfeat;
          awp->firstfeat = FALSE;

          if (awp->afp != NULL) {
            DoImmediateFormat (awp->afp, (BaseBlockPtr) fbp);
          }
        }
      }
    }
  }

  if (awp->onlyNearFeats) return;

  if (awp->nearFeatsSuppress && awp->featseen) return;

  if (! awp->farFeatsSuppress) {

    if (bsp->repr == Seq_repr_seg || bsp->repr == Seq_repr_delta || bsp->repr == Seq_repr_ref) {

      /* get start time for 25 second timeout in Web Entrez far WGS records */

      if (awp->farFeatTimeLimit) {
        awp->farFeatStartTime = GetSecs ();
      }

      /* if not farFeatsSuppress now collect features on remote segments in MASTER_STYLE */

      SeqMgrExploreSegments (bsp, (Pointer) awp, GetFeatsOnSeg);

      awp->wgsaccnlist = ValNodeFreeData (awp->wgsaccnlist);
    }
  }
}

NLM_EXTERN void AddFeatStatsBlock (
  Asn2gbWorkPtr awp
)

{
  IntAsn2gbJobPtr  ajp;
  BaseBlockPtr     bbp;
  BioseqPtr        bsp;
  StringItemPtr    ffstring;

  if (awp == NULL) return;
  ajp = awp->ajp;
  if ( ajp == NULL ) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;

  if (awp->format == EMBL_FMT || awp->format == EMBLPEPT_FMT) return;

  bbp = Asn2gbAddBlock (awp, FEAT_STATS_BLOCK, sizeof (BaseBlock));
  if (bbp != NULL) {
    ffstring = FFGetString (ajp);
    if (ffstring != NULL) {
      FFStartPrint (ffstring, awp->format, 0, 12, "FEATSTATS", 12, 0, 0, NULL, FALSE);
    
      FFAddOneString (ffstring, "placeholder", FALSE, FALSE, TILDE_TO_SPACES);
  
      bbp->string = FFEndPrint (ajp, ffstring, awp->format, 12, 12, 0, 0, NULL);
      FFRecycleString(ajp, ffstring);
    }

    if (awp->afp != NULL) {
      DoImmediateFormat (awp->afp, bbp);
    }
  }
}

