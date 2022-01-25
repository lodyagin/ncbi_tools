/*   scantest.c
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
* File Name:  scantest.c
*
* Author:  Kans
*
* Version Creation Date:   1/20/95
*
* $Revision: 6.13 $
*
* File Description: 
*       template for custom scans of ASN.1 release files
*       (was - scans through sequence records on the Entrez discs)
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
#include <sqnutils.h>
#include <explore.h>
#include <toasn3.h>

typedef struct appflags {
  Boolean  binary;
  Boolean  compressed;
  Boolean  verbose;
  FILE     *fp;
  Char     id [64];
} AppFlagData, PNTR AppFlagPtr;

static CharPtr Se2Str (
  SeqEntryPtr sep
)

{
  AsnIoBSPtr    aibp;
  ByteStorePtr  bs;
  CharPtr       str;

  if (sep == NULL) return NULL;

  bs = BSNew (1000);
  if (bs == NULL) return NULL;
  aibp = AsnIoBSOpen ("w", bs);
  if (aibp == NULL) return NULL;

  SeqEntryAsnWrite (sep, aibp->aip, NULL);

  AsnIoFlush (aibp->aip);
  AsnIoBSClose (aibp);

  str = BSMerge (bs, NULL);
  BSFree (bs);

  return str;
}

typedef struct chgdata {
  Boolean       rubisco;
  Boolean       rbc;
  Boolean       its;
  Boolean       sgml;
  Boolean       rnaother;
  Boolean       trnanote;
  Boolean       oldbiomol;
  Boolean       badname;
  Int4          protdesc;
  Int4          sfpnote;
  Int4          gbsource;
  Int4          cdsconf;
  AppFlagPtr    afp;
} ChangeData, PNTR ChangeDataPtr;

static Boolean IsRubisco (
  CharPtr name
)

{
  return (StringICmp (name, "rubisco large subunit") == 0 ||
          StringICmp (name, "rubisco small subunit") == 0);
}

static Boolean IsRbc (
  CharPtr name
)

{
  return (StringICmp (name, "RbcL") == 0 ||
          StringICmp (name, "RbcS") == 0);
}

static Boolean IsITS (
  CharPtr name
)

{
  return (StringICmp (name, "its1") == 0 ||
          StringICmp (name, "its 1") == 0 ||
          StringICmp (name, "its2") == 0 ||
          StringICmp (name, "its 2") == 0 ||
          StringICmp (name, "its3") == 0 ||
          StringICmp (name, "its 3") == 0 ||
          StringICmp (name, "Ribosomal DNA internal transcribed spacer 1") == 0 ||
          StringICmp (name, "Ribosomal DNA internal transcribed spacer 2") == 0 ||
          StringICmp (name, "Ribosomal DNA internal transcribed spacer 3") == 0 ||
          StringICmp (name, "internal transcribed spacer 1 (ITS1)") == 0 ||
          StringICmp (name, "internal transcribed spacer 2 (ITS2)") == 0 ||
          StringICmp (name, "internal transcribed spacer 3 (ITS3)") == 0);
}

static Boolean HasSgml (
  CharPtr  str
)

{
  Int2  ascii_len;
  Char  buf [1024];

  if (StringHasNoText (str)) return FALSE;

  ascii_len = Sgml2AsciiLen (str);
  if (ascii_len + 2 > sizeof (buf)) return FALSE;

  Sgml2Ascii (str, buf, ascii_len + 1);
  if (StringCmp (str, buf) != 0) {
    return TRUE;
  }

  return FALSE;
}

static void ScoreFeature (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  ChangeDataPtr  cdp;
  CharPtr        comment;
  CdRegionPtr    crp;
  CharPtr        desc;
  GBQualPtr      gbq;
  GeneRefPtr     grp;
  CharPtr        name;
  ProtRefPtr     prp;
  Uint1          residue;
  RnaRefPtr      rrp;
  CharPtr        str;
  ValNodePtr     vnp;

  if (sfp == NULL) return;
   cdp = (ChangeDataPtr) userdata;
   if (cdp == NULL) return;

  comment = sfp->comment;
  if (StringDoesHaveText (comment)) {
    (cdp->sfpnote)++;
  }

  /* skip feature types that do not use data.value.ptrvalue */
  switch (sfp->data.choice) {
    case SEQFEAT_COMMENT:
    case SEQFEAT_BOND:
    case SEQFEAT_SITE:
    case SEQFEAT_PSEC_STR:
      return;
    default:
      break;
  }

  if (sfp->data.value.ptrvalue == NULL) return;

  switch (sfp->data.choice) {
    case SEQFEAT_GENE:
      grp = (GeneRefPtr) sfp->data.value.ptrvalue;
      if (HasSgml (grp->locus)) {
        cdp->sgml = TRUE;
      }
      if (HasSgml (grp->desc)) {
        cdp->sgml = TRUE;
      }
      for (vnp = grp->syn; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        if (HasSgml (str)) {
          cdp->sgml = TRUE;
        }
      }
      break;
    case SEQFEAT_CDREGION:
      crp = (CdRegionPtr) sfp->data.value.ptrvalue;
      if (crp->conflict) {
        (cdp->cdsconf)++;
      }
      break;
    case SEQFEAT_PROT:
      prp = (ProtRefPtr) sfp->data.value.ptrvalue;
      desc = prp->desc;
      if (StringDoesHaveText (desc)) {
        (cdp->protdesc)++;
      }
      for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        if (IsRubisco (str)) {
          cdp->rubisco = TRUE;
        }
        if (IsRbc (str)) {
          cdp->rbc = TRUE;
        }
      }
      break;
    case SEQFEAT_RNA :
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp->type == 255 && rrp->ext.choice == 1) {
        name = (CharPtr) rrp->ext.value.ptrvalue;
        if (StringCmp (name, "misc_RNA") == 0) {
          for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
            if (StringCmp (gbq->qual, "product") != 0) continue;
            name = gbq->val;
            if (StringHasNoText (name)) continue;
            if (IsITS (name)) {
              cdp->its = TRUE;
            }
          }
        } else if (StringCmp (name, "ncRNA") == 0 || StringCmp (name, "tmRNA") == 0) {
        } else {
          cdp->rnaother = TRUE;
          if (IsITS (name)) {
            cdp->its = TRUE;
          }
        }
      } else if (rrp->type == 3 && rrp->ext.choice == 2) {
        if (StringDoesHaveText (comment)) {
          if (StringNCmp (comment, "aa: ", 4) == 0) {
            comment += 4;
          }
          residue = FindTrnaAA3 (comment);
          if (residue > 0 && residue != 255) {
            cdp->trnanote = TRUE;
          }
          residue = FindTrnaAA (comment);
          if (residue > 0 && residue != 255) {
            cdp->trnanote = TRUE;
          }
        }
      }
      break;
    default:
      break;
  }
}

static void ScoreDescriptor (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  ChangeDataPtr  cdp;
  GBBlockPtr     gbp;
  MolInfoPtr     mip;

  if (sdp == NULL) return;
  cdp = (ChangeDataPtr) userdata;
  if (cdp == NULL) return;

  switch (sdp->choice) {
    case Seq_descr_genbank :
      gbp = (GBBlockPtr) sdp->data.ptrvalue;
      if (gbp != NULL) {
        if (StringDoesHaveText (gbp->source)) {
          (cdp->gbsource)++;
        }
      }
      break;
    case Seq_descr_molinfo :
      mip = (MolInfoPtr) sdp->data.ptrvalue;
      if (mip != NULL) {
        switch (mip->biomol) {
          case MOLECULE_TYPE_SNRNA:
          case MOLECULE_TYPE_SCRNA:
          case MOLECULE_TYPE_SNORNA:
            cdp->oldbiomol = TRUE;
            break;
          default :
            break;
        }
      }
      break;
    default :
      break;
  }
}

static void CheckForChanges (
  SeqEntryPtr sep,
  ChangeDataPtr cdp
)

{
  if (sep == NULL || cdp == NULL) return;

  VisitFeaturesInSep (sep, (Pointer) cdp, ScoreFeature);
  VisitDescriptorsInSep (sep, (Pointer) cdp, ScoreDescriptor);
}

static void ModGenes (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  ModernizeGeneFields (sfp);
}

static void ModRNAs (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  ModernizeRNAFields (sfp);
}

static void ModPCRs (
  BioSourcePtr biop,
  Pointer userdata
)

{
  BoolPtr         namP;
  PCRPrimerPtr    ppp;
  PCRReactionPtr  prp;

  if (biop == NULL) return;

  ModernizePCRPrimers (biop);

  namP = (BoolPtr) userdata;
  if (namP == NULL) return;

  for (prp = biop->pcr_primers; prp != NULL; prp = prp->next) {
    if (prp->forward == NULL || prp->reverse == NULL) {
      *namP = TRUE;
      return;
    }
    for (ppp = prp->forward; ppp != NULL; ppp = ppp->next) {
      if (StringHasNoText (ppp->seq) && StringDoesHaveText (ppp->name)) {
        *namP = TRUE;
        return;
      }
    }
    for (ppp = prp->reverse; ppp != NULL; ppp = ppp->next) {
      if (StringHasNoText (ppp->seq) && StringDoesHaveText (ppp->name)) {
        *namP = TRUE;
        return;
      }
    }
  }
}

static void TestForRubisco (
  CharPtr str,
  AppFlagPtr afp,
  CharPtr prefix,
  CharPtr remainder
)

{
  if (StringHasNoText (str)) return;
  if (afp == NULL || afp->fp == NULL) return;

  if (StringICmp (str, "ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit") == 0) return;
  if (StringICmp (str, "ribulose-1,5-bisphosphate carboxylase/oxygenase small subunit") == 0) return;
  if (StringStr (str, "ribulose") == NULL || StringStr (str, "bisphosphate") == NULL) return;

  if (StringHasNoText (prefix)) {
    prefix = "?";
  }

  if (StringStr (str, "methyltransferase") == NULL) {
    if (StringICmp (str, "ribulose 1,5-bisphosphate carboxylase/oxygenase large subunit") == 0 ||
        StringICmp (str, "ribulose 1,5-bisphosphate carboxylase large subunit") == 0 ||
        StringICmp (str, "ribulose bisphosphate carboxylase large subunit") == 0 ||
        StringICmp (str, "ribulose-bisphosphate carboxylase large subunit") == 0 ||
        StringICmp (str, "ribulose-1,5-bisphosphate carboxylase large subunit") == 0 ||
        StringICmp (str, "ribulose-1,5-bisphosphate carboxylase, large subunit") == 0 ||
        StringICmp (str, "large subunit of ribulose-1,5-bisphosphate carboxylase/oxygenase") == 0 ||
        StringICmp (str, "ribulose-1,5-bisphosphate carboxylase oxygenase large subunit") == 0 ||
        StringICmp (str, "ribulose bisphosphate carboxylase large chain") == 0 ||
        StringICmp (str, "ribulose 1,5-bisphosphate carboxylase-oxygenase large subunit") == 0 ||
        StringICmp (str, "ribulose bisphosphate carboxylase oxygenase large subunit") == 0 ||
        StringICmp (str, "ribulose 1,5 bisphosphate carboxylase large subunit") == 0 ||
        StringICmp (str, "ribulose-1,5-bisphosphate carboxylase/oxygenase, large subunit") == 0 ||
        StringICmp (str, "large subunit of ribulose-1,5-bisphosphate carboxylase/oxgenase") == 0 ||
        StringICmp (str, "ribulose bisphosphate carboxylase/oxygenase large subunit") == 0 ||
        StringICmp (str, "ribulose-1,5-bisphosphate carboxylase oxygenase, large subunit") == 0 ||
        StringICmp (str, "ribulose 5-bisphosphate carboxylase, large subunit") == 0 ||
        StringICmp (str, "ribulosebisphosphate carboxylase large subunit") == 0 ||
        StringICmp (str, "ribulose bisphosphate large subunit") == 0 ||
        StringICmp (str, "ribulose 1,5 bisphosphate carboxylase/oxygenase large subunit") == 0 ||
        StringICmp (str, "ribulose 1,5-bisphosphate carboxylase/oxygenase large chain") == 0 ||
        StringICmp (str, "large subunit ribulose-1,5-bisphosphate carboxylase/oxygenase") == 0 ||
        StringICmp (str, "ribulose-bisphosphate carboxylase, large subunit") == 0 ||
        StringICmp (str, "ribulose-1, 5-bisphosphate carboxylase/oxygenase large-subunit") == 0) {
      if (afp->verbose) {
        fprintf (afp->fp, "%s\t%s\t%s\n", prefix, afp->id, str);
      } else {
        fprintf (afp->fp, "%s %s\n", prefix, afp->id);
      }
      fflush (afp->fp);
      return;
    }
  }

  if (StringHasNoText (remainder)) {
    remainder = "?";
  }
  if (afp->verbose) {
    fprintf (afp->fp, "%s\t%s\t%s\n", remainder, afp->id, str);
  } else {
    fprintf (afp->fp, "%s %s\n", remainder, afp->id);
  }
  fflush (afp->fp);
}

static void TrailingCommaFix (
  CharPtr str,
  AppFlagPtr afp,
  CharPtr prefix
)

{
  Char    ch;
  size_t  len;

  if (StringHasNoText (str)) return;
  len = StringLen (str);
  if (len < 1) return;
  ch = str [len - 1];
  while (ch == ' ' && len > 2) {
    len--;
    ch = str [len - 1];
  }
  if (ch == ',') {
    if (afp != NULL && afp->verbose && afp->fp != NULL) {
      str [len] = '\0';
      if (StringHasNoText (prefix)) {
        prefix = "?";
      }
      fprintf (afp->fp, "%s\t%s\t%s\n", prefix, afp->id, str);
      fflush (afp->fp);
    }
    str [len - 1] = '_';
    str [len] = '\0';
  }
}

static void RnaProtCmntTrailingCommaFix (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  AppFlagPtr  afp;
  ProtRefPtr  prp;
  RnaRefPtr   rrp;
  CharPtr     str;
  ValNodePtr  vnp;

  if (sfp == NULL) return;
  afp = (AppFlagPtr) userdata;
  if (afp == NULL) return;

  str = sfp->comment;
  if (StringDoesHaveText (str)) {
    TrailingCommaFix (str, afp, "SFPCOMM");
  }

  if (sfp->data.choice == SEQFEAT_PROT) {
    prp = (ProtRefPtr) sfp->data.value.ptrvalue;
    /* turn trailing space into trailing underscore for validator */
    for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      TrailingCommaFix (str, afp, "PRTCOMM");
      TestForRubisco (str, afp, "RIBBIS", "RIBREM");
    }
  } else if (sfp->data.choice == SEQFEAT_RNA) {
    rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
    /* turn trailing space into trailing underscore for validator */
    if (rrp->ext.choice == 1) {
      str = rrp->ext.value.ptrvalue;
      if (StringDoesHaveText (str)) {
        TrailingCommaFix (str, afp, "RNACOMM");
      }
    }
  }
}

static void LookForBadAuth (
  NameStdPtr nsp,
  Pointer userdata
)

{
  AppFlagPtr     afp;
  ChangeDataPtr  cdp;
  Char           ch;
  Int2           i;
  Boolean        is_bad = FALSE;
  CharPtr        prefix = "\t";
  CharPtr        str;

  if (nsp == NULL) return;
  cdp = (ChangeDataPtr) userdata;
  if (cdp == NULL) return;
  afp = cdp->afp;
  if (afp == NULL) return;

  for (i = 0; i < 6; i++) {
    str = nsp->names [i];
    if (StringHasNoText (str)) continue;
    ch = *str;
    while (ch != '\0') {
      if (IS_DIGIT (ch)) {
        cdp->badname = TRUE;
        is_bad = TRUE;
      }
      str++;
      ch = *str;
    }
  }

  if (is_bad && afp->fp != NULL && afp->verbose) {
    fprintf (afp->fp, "%s\t%s", "AUTHOR", afp->id);
    for (i = 0; i < 6; i++) {
      str = nsp->names [i];
      if (StringHasNoText (str)) continue;
      fprintf (afp->fp, "%s%s", prefix, str);
      prefix = " | ";
    }
    fprintf (afp->fp, "\n");
    fflush (afp->fp);
  }
}

static void LookForBadPub (
  PubdescPtr pdp,
  Pointer userdata
)

{
  VisitAuthorsInPub (pdp, userdata, LookForBadAuth);
}

static void CommentDescrTrailingCommaFix (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  AppFlagPtr  afp;
  CharPtr     str;

  if (sdp == NULL || sdp->choice != Seq_descr_comment) return;
  afp = (AppFlagPtr) userdata;
  if (afp == NULL) return;

  str = (CharPtr) sdp->data.ptrvalue;
  if (StringDoesHaveText (str)) {
    TrailingCommaFix (str, afp, "DSCCOMM");
  }
}

static void DoReport (
  SeqEntryPtr sep,
  AppFlagPtr afp
)

{
  Boolean     bsec = FALSE, cma = FALSE, norm = FALSE, ssec = FALSE;
  Boolean     gen = FALSE, ncr = FALSE, pcr = FALSE, nam = FALSE;
  ChangeData  cdbefore, cdafter;
  CharPtr     str = NULL, tmp = NULL;

  if (sep == NULL || afp == NULL) return;

  MemSet ((Pointer) &cdbefore, 0, sizeof (ChangeData));
  MemSet ((Pointer) &cdafter, 0, sizeof (ChangeData));

  cdbefore.afp = afp;
  cdafter.afp = afp;

  CheckForChanges (sep, &cdbefore);
 
  str = Se2Str (sep);
  NormalizeDescriptorOrder (sep);
  tmp = Se2Str (sep);
  if (StringCmp (str, tmp) != 0) {
    norm = TRUE;
  }
  MemFree (str);
  str = tmp;

  VisitFeaturesInSep (sep, (Pointer) afp, RnaProtCmntTrailingCommaFix);
  VisitDescriptorsInSep (sep, (Pointer) afp, CommentDescrTrailingCommaFix);
  tmp = Se2Str (sep);
  if (StringCmp (str, tmp) != 0) {
    cma = TRUE;
  }
  MemFree (str);
  str = tmp;

  BasicSeqEntryCleanup (sep);
  tmp = Se2Str (sep);
  if (StringCmp (str, tmp) != 0) {
    bsec = TRUE;
  }
  MemFree (str);
  str = tmp;

  VisitPubdescsInSep (sep, (Pointer) &cdbefore, LookForBadPub);

  VisitFeaturesInSep (sep, NULL, ModGenes);
  tmp = Se2Str (sep);
  if (StringCmp (str, tmp) != 0) {
    gen = TRUE;
  }
  MemFree (str);
  str = tmp;

  VisitFeaturesInSep (sep, NULL, ModRNAs);
  tmp = Se2Str (sep);
  if (StringCmp (str, tmp) != 0) {
    ncr = TRUE;
  }
  MemFree (str);
  str = tmp;

  VisitBioSourcesInSep (sep, (Pointer) &nam, ModPCRs);
  tmp = Se2Str (sep);
  if (StringCmp (str, tmp) != 0) {
    pcr = TRUE;
  }
  MemFree (str);
  str = tmp;

  SeriousSeqEntryCleanup (sep, NULL, NULL);
  tmp = Se2Str (sep);
  if (StringCmp (str, tmp) != 0) {
    ssec = TRUE;
  }
  MemFree (str);
  str = tmp;

  CheckForChanges (sep, &cdafter);

  MemFree (str);

  if (ssec) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "SSEC %s\n", afp->id);
      fflush (afp->fp);
    }
  } else if (bsec) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "BSEC %s\n", afp->id);
      fflush (afp->fp);
    }
  } else if (norm) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "NORM %s\n", afp->id);
      fflush (afp->fp);
    }
  } else {
    /*
    if (afp->fp != NULL) {
      fprintf (afp->fp, "OKAY %s\n", afp->id);
      fflush (afp->fp);
    }
    */
  }

  if (cma) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "CMA %s\n", afp->id);
      fflush (afp->fp);
    }
  }

  if (gen) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "GEN %s\n", afp->id);
      fflush (afp->fp);
    }
  }
  if (ncr) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "NCR %s\n", afp->id);
      fflush (afp->fp);
    }
  }
  if (pcr) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "PCR %s\n", afp->id);
      fflush (afp->fp);
    }
  }
  if (nam) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "NAM %s\n", afp->id);
      fflush (afp->fp);
    }
  }

  if (cdbefore.rubisco) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "RUB %s\n", afp->id);
      fflush (afp->fp);
    }
  }
  if (cdbefore.rbc) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "RBC %s\n", afp->id);
      fflush (afp->fp);
    }
  }
  if (cdbefore.its) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "ITS %s\n", afp->id);
      fflush (afp->fp);
    }
  }
  if (cdbefore.sgml) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "SGM %s\n", afp->id);
      fflush (afp->fp);
    }
  }
  if (cdbefore.rnaother) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "RNA %s\n", afp->id);
      fflush (afp->fp);
    }
  }
  if (cdbefore.trnanote) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "TRN %s\n", afp->id);
      fflush (afp->fp);
    }
  }
  if (cdbefore.oldbiomol) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "MOL %s\n", afp->id);
      fflush (afp->fp);
    }
  }
  if (cdbefore.badname) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "AUT %s\n", afp->id);
      fflush (afp->fp);
    }
  }

  if (cdbefore.protdesc != cdafter.protdesc) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "PRT %s\n", afp->id);
      fflush (afp->fp);
    }
  }
  if (cdbefore.sfpnote != cdafter.sfpnote) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "COM %s\n", afp->id);
      fflush (afp->fp);
    }
  }
  if (cdbefore.gbsource != cdafter.gbsource) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "SRC %s\n", afp->id);
      fflush (afp->fp);
    }
  }
  if (cdbefore.cdsconf != cdafter.cdsconf) {
    if (afp->fp != NULL) {
      fprintf (afp->fp, "CNF %s\n", afp->id);
      fflush (afp->fp);
    }
  }
}

static void DoRecord (SeqEntryPtr sep, Pointer userdata)

{
  AppFlagPtr   afp;
  BioseqPtr    fbsp;
  SeqEntryPtr  fsep;
  SeqIdPtr     sip, siphead;

  if (sep == NULL) return;
  afp = (AppFlagPtr) userdata;
  if (afp == NULL) return;

  fsep = FindNthBioseq (sep, 1);
  if (fsep == NULL) return;
  fbsp = (BioseqPtr) fsep->data.ptrvalue;
  if (fbsp == NULL) return;

  siphead = SeqIdSetDup (fbsp->id);
  for (sip = siphead; sip != NULL; sip = sip->next) {
    SeqIdStripLocus (sip);
  }
  SeqIdWrite (siphead, afp->id, PRINTID_FASTA_LONG, sizeof (afp->id));
  SeqIdSetFree (siphead);

  DoReport (sep, afp);
}

static void ProcessOneRecord (
  CharPtr filename,
  Pointer userdata
)

{
  AppFlagPtr  afp;

  if (StringHasNoText (filename)) return;
  afp = (AppFlagPtr) userdata;
  if (afp == NULL) return;

  if (StringStr (filename, "gbcon") != NULL ||
      StringStr (filename, "gbest") != NULL ||
      StringStr (filename, "gbgss") != NULL ||
      StringStr (filename, "gbhtg") != NULL ||
      StringStr (filename, "gbsts") != NULL) {
    printf ("Skipping %s\n", filename);
    return;
  }

  printf ("%s\n", filename);
  fflush (stdout);

  fprintf (afp->fp, "%s\n", filename);
  fflush (afp->fp);

  ScanBioseqSetRelease (filename, afp->binary, afp->compressed, (Pointer) afp, DoRecord);

  fprintf (afp->fp, "\n");
  fflush (afp->fp);
}

#define p_argInputPath    0
#define i_argInputFile    1
#define o_argOutputFile   2
#define f_argFilter       3
#define x_argSuffix       4
#define u_argRecurse      5
#define b_argBinary       6
#define c_argCompressed   7
#define v_argVerbose      8

Args myargs [] = {
  {"Path to Files", NULL, NULL, NULL,
    TRUE, 'p', ARG_STRING, 0.0, 0, NULL},
  {"Input File Name", NULL, NULL, NULL,
    TRUE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
  {"Output File Name", NULL, NULL, NULL,
    TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
  {"Substring Filter", NULL, NULL, NULL,
    TRUE, 'f', ARG_STRING, 0.0, 0, NULL},
  {"File Selection Suffix", ".aso", NULL, NULL,
    TRUE, 'x', ARG_STRING, 0.0, 0, NULL},
  {"Recurse", "F", NULL, NULL,
    TRUE, 'u', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Bioseq-set is Binary", "F", NULL, NULL,
    TRUE, 'b', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Bioseq-set is Compressed", "F", NULL, NULL,
    TRUE, 'c', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Verbose Output", "F", NULL, NULL,
    TRUE, 'v', ARG_BOOLEAN, 0.0, 0, NULL},
};

extern Int2 Main (void)

{
  AppFlagData  afd;
  Boolean      dorecurse;
  CharPtr      filter, infile, outfile, directory, suffix;

  /* standard setup */

  ErrSetFatalLevel (SEV_MAX);
  ErrClearOptFlags (EO_SHOW_USERSTR);
  ErrSetLogfile ("stderr", ELOG_APPEND);
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
  if (! SeqCodeSetLoad ()) {
    Message (MSG_FATAL, "SeqCodeSetLoad failed");
    return 1;
  }
  if (! GeneticCodeTableLoad ()) {
    Message (MSG_FATAL, "GeneticCodeTableLoad failed");
    return 1;
  }

  /* process command line arguments */

  if (! GetArgs ("scantest", sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  MemSet ((Pointer) &afd, 0, sizeof (AppFlagData));

  directory = (CharPtr) myargs [p_argInputPath].strvalue;
  infile = (CharPtr) myargs [i_argInputFile].strvalue;
  outfile = (CharPtr) myargs [o_argOutputFile].strvalue;
  filter = (CharPtr) myargs [f_argFilter].strvalue;
  suffix = (CharPtr) myargs [x_argSuffix].strvalue;
  dorecurse = (Boolean) myargs [u_argRecurse].intvalue;
  afd.binary = (Boolean) myargs [b_argBinary].intvalue;
  afd.compressed = (Boolean) myargs [c_argCompressed].intvalue;
  afd.verbose = (Boolean) myargs [v_argVerbose].intvalue;

  afd.fp = FileOpen (outfile, "w");
  if (afd.fp == NULL) {
    return 0;
  }

  if (StringDoesHaveText (directory)) {

    DirExplore (directory, NULL, suffix, dorecurse, ProcessOneRecord, (Pointer) &afd);

  } else if (StringDoesHaveText (infile)) {

    ProcessOneRecord (infile, &afd);
  }

  FileClose (afd.fp);

  return 0;
}
