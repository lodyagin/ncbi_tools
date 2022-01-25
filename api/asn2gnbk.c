/*   asn2gnbk.c
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
* File Name:  asn2gnbk.c
*
* Author:  Karl Sirotkin, Tom Madden, Tatiana Tatusov, Jonathan Kans
*
* Version Creation Date:   10/21/98
*
* $Revision: 6.120 $
*
* File Description:  New GenBank flatfile generator
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
#include <seqport.h>
#include <sequtil.h>
#include <sqnutils.h>
#include <subutil.h>
#include <tofasta.h>
#include <explore.h>
#include <ffprint.h>
#include <utilpub.h>
#include <jzmisc.h>
#include <gbfeat.h>
#include <gbftdef.h>
#include <asn2gnbk.h>

#ifdef WIN_MAC
#if __profile__
#include <Profiler.h>
#endif
#endif

#define ASN2FF_EMBL_MAX 78
#define ASN2FF_GB_MAX 79


/* structure for storing working parameters while building asn2gb_job structure */

typedef struct asn2gbwork {
  Asn2gbJobPtr  ajp;
  Uint2         entityID;

  FmtType       format;
  ModType       mode;
  StlType       style;

  ValNodePtr    pubhead;    /* for collecting publications */
  ValNodePtr    srchead;    /* for collecting biosources */

  /* linked lists of paragraphs, sections, blocks */

  ValNodePtr    sectionList;
  ValNodePtr    blockList;    /* reset for each new section */

  /* most recent node of linked lists, for quickly adding next node */

  ValNodePtr    lastsection;
  ValNodePtr    lastblock;    /* reset for each new section */

  Int4          currsection;

  /* section fields needed for populating blocks */

  BioseqPtr     target;
  BioseqPtr     parent;
  BioseqPtr     bsp;
  SeqLocPtr     slp;
  Uint2         seg;
  Int4          numsegs;
  Int4          from;
  Int4          to;

  SeqFeatPtr    lastsfp;
  SeqAnnotPtr   lastsap;
  Int4          lastleft;
  Int4          lastright;

  SeqSubmitPtr  ssp;
  Boolean       hup;
} Asn2gbWork, PNTR Asn2gbWorkPtr;

/* array for assigning biosource and feature data fields to qualifiers */
/* should be allocated to MAX (ASN2GNBK_TOTAL_SOURCE, ASN2GNBK_TOTAL_FEATUR) */

typedef union qualval {
  CharPtr       str;
  Boolean       bool;
  Int4          num;
  ValNodePtr    vnp;
  GBQualPtr     gbq;
  OrgModPtr     omp;
  SubSourcePtr  ssp;
  CodeBreakPtr  cbp;
  SeqLocPtr     slp;
  SeqIdPtr      sip;
  tRNAPtr       trp;
} QualVal, PNTR QualValPtr;

/* structure passed to individual paragraph format functions */

typedef struct asn2gbformat {
  Asn2gbJobPtr      ajp;
  Asn2gbSectionPtr  asp;
  QualValPtr        qvp;

  FmtType           format;
  ModType           mode;
  StlType           style;
} Asn2gbFormat, PNTR Asn2gbFormatPtr;


/* Seq-hist replacedBy is preformatted into string field, */
/* then comment descriptors, Map location:, and Region:, */
/* then comment features, finally HTGS */

typedef struct comment_block {
  ASN2GB_BASE_BLOCK
  Boolean           first;
} CommentBlock, PNTR CommentBlockPtr;

/* internal reference block has fields on top of ReferenceBlock fields */

typedef struct int_ref_block {
  ReferenceBlock    rb;
  DatePtr           date;     /* internal sorting use only */
  SeqLocPtr         loc;      /* final location on target bioseq */
  CharPtr           authstr;  /* author string */
  Uint2             index;    /* index if feature on target bioseq */
  Boolean           justuids; /* gibb pub with uids and Figure, etc. */
  CharPtr           fig;      /* figure string from equivalent gibb pub */
  CharPtr           maploc;   /* maploc string from equivalent gibb pub */
  Boolean           poly_a;   /* poly_a field from equivalent gibb pub */
} IntRefBlock, PNTR IntRefBlockPtr;

/* internal source block has fields on top of BaseBlock fields */

typedef struct int_src_block {
  BaseBlock         bb;
  Boolean           is_focus;
  Uint4             orghash;
  Uint4             modhash;
  Uint4             subhash;
  Uint4             xrfhash;
  SeqLocPtr         loc;     /* final location on target bioseq */
} IntSrcBlock, PNTR IntSrcBlockPtr;

/* internal feature block has fields on top of FeatureBlock fields */

typedef struct int_feat_block {
  FeatureBlock      fb;
  Boolean           mapToNuc;
  Boolean           mapToProt;
  Boolean           isCDS;     /* set if using IntCdsBlock */
} IntFeatBlock, PNTR IntFeatBlockPtr;

/* internal cds block has fields on top of IntFeatBlock fields */

typedef struct int_cds_block {
  IntFeatBlock      ifb;
  FeatureBlock      fb;
  CharPtr           fig;    /* figure string from pub */
  CharPtr           maploc; /* maploc string from pub */
} IntCdsBlock, PNTR IntCdsBlockPtr;


/* ********************************************************************** */

/* utility functions */

static CharPtr MergeValNodeStrings (
  ValNodePtr list
)

{
  size_t      len;
  CharPtr     ptr;
  CharPtr     str;
  CharPtr     tmp;
  ValNodePtr  vnp;


  if (list == NULL) return NULL;

  for (vnp = list, len = 0; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    len += StringLen (str);
  }
  if (len == 0) return NULL;

  ptr = MemNew (sizeof (Char) * (len + 2));
  if (ptr == NULL) return NULL;

  for (vnp = list, tmp = ptr; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    tmp = StringMove (tmp, str);
  }

  return ptr;
}

static void AddValNodeString (
  ValNodePtr PNTR head,
  CharPtr prefix,
  CharPtr string,
  CharPtr suffix
)

{
  Char     buf [256];
  CharPtr  freeme = NULL;
  size_t   len;
  CharPtr  newstr;
  CharPtr  strptr;

  len = StringLen (prefix) + StringLen (string) + StringLen (suffix);
  if (len == 0) return;

  if (len < sizeof (buf)) {

    /* if new string fits in stack buffer, no need to allocate */

    MemSet ((Pointer) buf, 0, sizeof (buf));
    newstr = buf;

  } else {

    /* new string bigger than stack buffer, so allocate sufficient string */

    newstr = (CharPtr) MemNew (sizeof (Char) * (len + 2));
    if (newstr == NULL) return;

    /* allocated string will be freed at end of function */

    freeme = newstr;
  }

  strptr = newstr;

  if (prefix != NULL) {
    strptr = StringMove (strptr, prefix);
  }

  if (string != NULL) {
    strptr = StringMove (strptr, string);
  }

  if (suffix != NULL) {
    strptr = StringMove (strptr, suffix);
  }

  /* currently just makes a valnode list, to be enhanced later */

  ValNodeCopyStr (head, 0, newstr);

  /* if large string was allocated, free it now */

  if (freeme != NULL) {
    MemFree (freeme);
  }
}

static Int2 gb_StartPrint (
  FmtType format,
  Boolean call_init_buff,
  Int2 gb_init_indent,
  Int2 gb_cont_indent,
  CharPtr gb_label,
  Int2 gb_tab_to,
  Int2 eb_init_indent,
  Int2 eb_cont_indent,
  CharPtr eb_line_prefix,
  Boolean eb_print_xx
)

{
  Int2  result = 0;

  if (call_init_buff) {
    init_buff ();
  }

  if (format == GENBANK_FMT || format == GENPEPT_FMT) {

    result = ff_StartPrint (gb_init_indent, gb_cont_indent, ASN2FF_GB_MAX, NULL);
    if (gb_label != NULL) {
      ff_AddString (gb_label);
    }
    if (gb_tab_to > 0) {
      TabToColumn (gb_tab_to);
    }

  } else if (format == EMBL_FMT || format == EMBLPEPT_FMT) {

    if (eb_print_xx) {
      PrintXX ();
    }
    result = ff_StartPrint (eb_init_indent, eb_cont_indent, ASN2FF_EMBL_MAX, eb_line_prefix);
  }

  return result;
}

/* convertQuotes retains double quotes in the prefix and suffix */

static void gb_AddString (
  CharPtr prefix,
  CharPtr string,
  CharPtr suffix,
  Boolean addPeriod,
  Boolean convertQuotes,
  Boolean expandTildes
)

{
  Char     buf [256];
  Char     ch;
  CharPtr  convertHere;
  CharPtr  freeme = NULL;
  size_t   len;
  CharPtr  newstr;
  CharPtr  strptr;

  len = StringLen (prefix) + StringLen (string) + StringLen (suffix);

  if (len == 0) {
    if (addPeriod) {
      ff_AddString (".");
    }
    return;
  }

  if (len + 1 < sizeof (buf)) {

    /* if new string fits in stack buffer, no need to allocate */

    MemSet ((Pointer) buf, 0, sizeof (buf));
    newstr = buf;

  } else {

    /* new string bigger than stack buffer, so allocate sufficient string */

    newstr = (CharPtr) MemNew (sizeof (Char) * (len + 3));
    if (newstr == NULL) return;

    /* allocated string will be freed at end of function */

    freeme = newstr;
  }

  strptr = newstr;

  if (prefix != NULL) {
    strptr = StringMove (strptr, prefix);
  }

  convertHere = strptr;
  if (string != NULL) {
    strptr = StringMove (strptr, string);
  }

  if (convertQuotes) {
    strptr = convertHere;
    ch = *strptr;
    while (ch != '\0') {
      if (ch == '\"') {
        *strptr = '\'';
      }
      strptr++;
      ch = *strptr;
    }
  }
  if (suffix != NULL) {
    strptr = StringMove (strptr, suffix);
  }

  if (addPeriod) {
    for (strptr = newstr + len - 1; strptr > newstr; strptr--) {
      ch = *strptr;
      if (ch == ' ' || ch == '\t' || ch == '~' || ch == '.') {
        *strptr = '\0';
      } else {
        break;
      }
    }
    if (*strptr != '.') {
      strptr++;
      *strptr = '.';
      strptr++;
      *strptr = '\0';
    }
  }

  if (expandTildes) {
    ff_AddStringWithTildes (newstr);
  } else {
    ff_AddString (newstr);
  }

  /* if large string was allocated, free it now */

  if (freeme != NULL) {
    MemFree (freeme);
  }
}

static CharPtr gb_MergeString (
  Boolean call_end_print
)

{
  if (call_end_print) {
    ff_EndPrint ();
  }

  return ff_MergeString ();
}

static CharPtr month_names [] = {
  "JAN", "FEB", "MAR", "APR", "MAY", "JUN",
  "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"
};

static CharPtr DateToGB (
  CharPtr buf,
  DatePtr dp
)

{
  Int2  day;
  Int2  month;
  Int2  year;

  if (buf != NULL) {
    *buf = '\0';
  }
  if (dp == NULL) return NULL;

  if (dp->data [0] == 0) {

    StringCpy (buf, dp->str);

  } else if (dp->data [0] == 1) {

    year = 1900 + (Int2) dp->data [1];
    month = (Int2) dp->data [2];
    day = (Int2) dp->data [3];

    if (month < 1) {
      month = 1;
    }
    if (day < 1) {
      day = 1;
    }

    if (day < 10) {
      sprintf (buf, "0%ld-%s-%ld",
               (long) day, month_names [month-1], (long) year);
    } else {
      sprintf(buf, "%ld-%s-%ld", 
               (long) day, month_names [month-1], (long) year);
    }
  }

  return buf;
}


/* ********************************************************************** */

/* format functions allocate printable string for given paragraph */

static CharPtr DefaultFormatBlock (
  Asn2gbFormatPtr afp,
  BaseBlockPtr bbp
)

{
  if (afp == NULL || bbp == NULL) return NULL;

  /* default format function assumes string pre-allocated by add block function */

  return StringSaveNoNull (bbp->string);
}

static CharPtr organellePrefix [] = {
  NULL,
  NULL,
  "Chloroplast ",
  "Chromoplast ",
  "Kinetoplast ",
  "Mitochondrion ",
  "Plastid ",
  "Macronuclear ",
  "Extrachrom ",
  "Plasmid ", 
  NULL,
  NULL,
  "Cyanelle ",
  "Proviral ",
  "Virion ",
  "Nucleomorph ",
  "Apicoplast ",
  "Leucoplast ",
  "Proplastid "
};

static CharPtr FormatSourceBlock (
  Asn2gbFormatPtr afp,
  BaseBlockPtr bbp
)

{
  BioSourcePtr       biop = NULL;
  CharPtr            common = NULL;
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext;
  GBBlockPtr         gbp = NULL;
  Uint1              genome;
  CharPtr            lineage = NULL;
  OrgNamePtr         onp;
  CharPtr            organelle = NULL;
  OrgRefPtr          orp;
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp;
  CharPtr            taxname = NULL;

  if (afp == NULL || bbp == NULL) return NULL;

  if (! StringHasNoText (bbp->string)) return StringSave (bbp->string);

  if (bbp->itemtype == OBJ_SEQDESC) {
    sdp = SeqMgrGetDesiredDescriptor (bbp->entityID, NULL, bbp->itemID, 0, NULL, &dcontext);
    if (sdp != NULL) {
      if (dcontext.seqdesctype == Seq_descr_source) {
        biop = (BioSourcePtr) sdp->data.ptrvalue;
      } else if (dcontext.seqdesctype == Seq_descr_genbank) {
        gbp = (GBBlockPtr) sdp->data.ptrvalue;
      }
    }
  } else if (bbp->itemtype == OBJ_SEQFEAT) {
    sfp = SeqMgrGetDesiredFeature (bbp->entityID, NULL, bbp->itemID, 0, NULL, &fcontext);
    if (sfp != NULL && fcontext.seqfeattype == SEQFEAT_BIOSRC) {
      biop = (BioSourcePtr) sfp->data.value.ptrvalue;
    }
  }
  if (gbp != NULL) {
    common = gbp->source;
  }
  if (biop != NULL) {
    genome = biop->genome;
    if (genome < 6 || genome == 12) {
      organelle = organellePrefix [genome];
    }
    orp = biop->org;
    if (orp != NULL) {
      taxname = orp->taxname;
      common = orp->common;
      onp = orp->orgname;
      if (onp != NULL) {
        lineage = onp->lineage;
      }
    }
  }

  if (StringHasNoText (common)) {
    common = taxname;
  }
  if (StringHasNoText (common)) {
    common = "Unknown.";
  }
  if (StringHasNoText (taxname)) {
    taxname = "Unknown.";
  }
  if (StringHasNoText (lineage)) {
    lineage = "Unclassified.";
  }

  if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {

    gb_StartPrint (afp->format, TRUE, 0, 12, "SOURCE", 13, 5, 5, "OS", TRUE);
    gb_AddString (NULL, common, NULL, TRUE, FALSE, FALSE);
    ff_EndPrint();

  } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {

    gb_StartPrint (afp->format, TRUE, 0, 12, "SOURCE", 13, 5, 5, "OS", TRUE);
    gb_AddString (organelle, taxname, NULL, FALSE, FALSE, FALSE);
    gb_AddString (" (", common, ")", FALSE, FALSE, FALSE);
    ff_EndPrint();
  }

  return gb_MergeString (FALSE);
}

static CharPtr FormatOrganismBlock (
  Asn2gbFormatPtr afp,
  BaseBlockPtr bbp
)

{
  BioSourcePtr       biop = NULL;
  CharPtr            common = NULL;
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext;
  Uint1              genome;
  CharPtr            lineage = NULL;
  OrgNamePtr         onp;
  CharPtr            organelle = NULL;
  OrgRefPtr          orp;
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp;
  CharPtr            taxname = NULL;

  if (afp == NULL || bbp == NULL) return NULL;

  if (! StringHasNoText (bbp->string)) return StringSave (bbp->string);

  if (bbp->itemtype == OBJ_SEQDESC) {
    sdp = SeqMgrGetDesiredDescriptor (bbp->entityID, NULL, bbp->itemID, 0, NULL, &dcontext);
    if (sdp != NULL && dcontext.seqdesctype == Seq_descr_source) {
      biop = (BioSourcePtr) sdp->data.ptrvalue;
    }
  } else if (bbp->itemtype == OBJ_SEQFEAT) {
    sfp = SeqMgrGetDesiredFeature (bbp->entityID, NULL, bbp->itemID, 0, NULL, &fcontext);
    if (sfp != NULL && fcontext.seqfeattype == SEQFEAT_BIOSRC) {
      biop = (BioSourcePtr) sfp->data.value.ptrvalue;
    }
  }
  if (biop != NULL) {
    genome = biop->genome;
    if (genome < 6 || genome == 12) {
      organelle = organellePrefix [genome];
    }
    orp = biop->org;
    if (orp != NULL) {
      taxname = orp->taxname;
      common = orp->common;
      onp = orp->orgname;
      if (onp != NULL) {
        lineage = onp->lineage;
      }
    }
  }

  if (StringHasNoText (common)) {
    common = taxname;
  }
  if (StringHasNoText (common)) {
    common = "Unknown.";
  }
  if (StringHasNoText (taxname)) {
    taxname = "Unknown.";
  }
  if (StringHasNoText (lineage)) {
    lineage = "Unclassified.";
  }

  if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {

    gb_StartPrint (afp->format, FALSE, 2, 12, "ORGANISM", 13, 5, 5, "OC", FALSE);
    gb_AddString (organelle, taxname, NULL, FALSE, FALSE, FALSE);
    ff_EndPrint();

    gb_StartPrint (afp->format, FALSE, 12, 12, NULL, 0, 5, 5, "OC", FALSE);
    gb_AddString (NULL, lineage, NULL, TRUE, FALSE, FALSE);
    ff_EndPrint();

  } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {

    gb_StartPrint (afp->format, FALSE, 12, 12, NULL, 0, 5, 5, "OC", FALSE);
    gb_AddString (NULL, lineage, NULL, TRUE, FALSE, FALSE);
    ff_EndPrint();
  }

  return gb_MergeString (FALSE);
}

/* format references section */

static AuthListPtr GetAuthListPtr (
  PubdescPtr pdp,
  CitSubPtr csp
)

{
  AuthListPtr  alp = NULL;
  CitArtPtr    cap;
  CitBookPtr   cbp;
  CitGenPtr    cgp;
  CitPatPtr    cpp;
  ValNodePtr   vnp;

  if (csp != NULL) {
    alp = csp->authors;
    if (alp != NULL) return alp;
  }
  if (pdp == NULL) return NULL;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    switch (vnp->choice) {
      case PUB_Gen :
        cgp = (CitGenPtr) vnp->data.ptrvalue;
        if (cgp != NULL) {
          alp = cgp->authors;
        }
        break;
      case PUB_Sub :
        csp = (CitSubPtr) vnp->data.ptrvalue;
        if (csp != NULL) {
          alp = csp->authors;
        }
        break;
      case PUB_Article :
        cap = (CitArtPtr) vnp->data.ptrvalue;
        if (cap != NULL) {
          alp = cap->authors;
        }
        break;
      case PUB_Book :
      case PUB_Proc :
      case PUB_Man :
        cbp = (CitBookPtr) vnp->data.ptrvalue;
        if (cbp != NULL) {
          alp = cbp->authors;
        }
        break;
      case PUB_Patent :
        cpp = (CitPatPtr) vnp->data.ptrvalue;
        if (cpp != NULL) {
          alp = cpp->authors;
        }
        break;
      default :
        break;
    }

    if (alp != NULL) return alp;
  }

  return NULL;
}

static CharPtr MakeSingleAuthorString (
  FmtType format,
  CharPtr prefix,
  CharPtr name,
  CharPtr initials,
  CharPtr suffix
)

{
  Char     ch;
  size_t   len;
  CharPtr  ptr;
  CharPtr  str;
  CharPtr  tmp;

  if (name == NULL) return NULL;

  /* !!! clean up 'et al' as (presumably) last author !!! */

  /* !!! temporary to suppress diff !!! */
  {
  Char dummy [10];
  if (StringLen (name) <= 6 &&
      (StringNICmp (name, "et al", 5) == 0 || StringNICmp (name, "et,al", 5) == 0)) {
    if (StringCmp (prefix, " and ") == 0) {
      prefix = NULL;
      dummy [0] = ' ';
      StringNCpy_0 (dummy + 1, name, sizeof (dummy) - 1);
      name = dummy;
    }
  }
  }
  /*
  if (StringLen (name) <= 6 &&
      (StringNICmp (name, "et al", 5) == 0 || StringNICmp (name, "et,al", 5) == 0)) {
    name = "et al.";
    if (StringCmp (prefix, " and ") == 0) {
      prefix = ", ";
    }
  }
  */

  len = StringLen (name) + StringLen (initials) + StringLen (suffix) + StringLen (prefix);
  str = MemNew (sizeof (Char) * (len + 4));
  if (str == NULL) return NULL;

  ptr = str;
  if (! StringHasNoText (prefix)) {
    ptr = StringMove (ptr, prefix);
  }

  /* initials and suffix to support structured name fields */

  tmp = StringMove (ptr, name);
  if (! StringHasNoText (initials)) {
    tmp = StringMove (tmp, ",");
    tmp = StringMove (tmp, initials);
  }
  if (! StringHasNoText (suffix)) {
    tmp = StringMove (tmp, " ");
    tmp = StringMove (tmp, suffix);
  }

  /* if embl, remove commas in individual names, starting after prefix */

  if (format == EMBL_FMT || format == EMBLPEPT_FMT) {
    tmp = ptr;
    ch = *tmp;
    while (ch != '\0') {
      if (ch == ',') {
        *tmp = ' ';
      }
      tmp++;
      ch = *tmp;
    }
  }

  return str;
}

static CharPtr GetAuthorsString (
  FmtType format,
  AuthListPtr alp
)

{
  AuthorPtr    ap;
  Int2         count;
  ValNodePtr   head = NULL;
  ValNodePtr   names;
  ValNodePtr   next;
  NameStdPtr   nsp;
  PersonIdPtr  pid;
  CharPtr      prefix = NULL;
  CharPtr      str;
  ValNodePtr   vnp;

  if (alp == NULL) return NULL;

  count = 0;
  if (alp->choice == 1) {
    for (names = alp->names; names != NULL; names = names->next) {
      next = names->next;
      if (next == NULL) {
        if (format == GENBANK_FMT || format == GENPEPT_FMT) {
          if (count == 0) {
            prefix = NULL;
          } else {
            prefix = " and ";
          }
        }
      }
      str = NULL;
      ap = (AuthorPtr) names->data.ptrvalue;
      if (ap != NULL) {
        pid = ap->name;
        if (pid != NULL) {
          if (pid->choice == 2) {
            nsp = (NameStdPtr) pid->data;
            if (nsp != NULL) {
              if (! StringHasNoText (nsp->names [0])) {
                str = MakeSingleAuthorString (format, prefix, nsp->names [0], nsp->names [4], nsp->names [5]);
              } else if (! StringHasNoText (nsp->names [3])) {
                str = MakeSingleAuthorString (format, prefix, nsp->names [3], NULL, NULL);
              }
            }
          } else if (pid->choice == 3 || pid->choice == 4) {
            str = MakeSingleAuthorString (format, prefix, (CharPtr) pid->data, NULL, NULL);
          }
        }
      }
      if (str != NULL) {
        ValNodeCopyStr (&head, 0, str);
        count++;
      }
      prefix = ", ";
    }
  } else if (alp->choice == 2 || alp->choice == 3) {
    for (vnp = alp->names; vnp != NULL; vnp = vnp->next) {
      next = vnp->next;
      if (next == NULL) {
        if (format == GENBANK_FMT || format == GENPEPT_FMT) {
          if (count == 0) {
            prefix = NULL;
          } else {
            prefix = " and ";
          }
        }
      }
      str = MakeSingleAuthorString (format, prefix, (CharPtr) vnp->data.ptrvalue, NULL, NULL);
      if (str != NULL) {
        ValNodeCopyStr (&head, 0, str);
        count++;
      }
      prefix = ", ";
    }
  }

  str = MergeValNodeStrings (head);

  ValNodeFreeData (head);

  return str;
}

/*
Strips all spaces in string in following manner. If the function
meet several spaces (spaces and tabs) in succession it replaces them
with one space. Strips all spaces after '(' and before ')'
*/

static void StrStripSpaces (
  CharPtr str
)

{
  CharPtr  new_str;

  if (str == NULL) return;

  new_str = str;
  while (*str != '\0') {
    *new_str++ = *str;
    if (*str == ' ' || *str == '\t' || *str == '(') {
      for (str++; *str == ' ' || *str == '\t'; str++) continue;
      if (*str == ')' || *str == ',') {
        new_str--;
      }
    } else {
      str++;
    }
  }
  *new_str = '\0';
}

static Boolean AllCaps (
  CharPtr p
)

{
  if (p == NULL) return FALSE;

  for (p++; p != NULL && *p != '\0'; p++) {
    if (IS_LOWER (*p)) return FALSE;
  }
  return TRUE;
}

static void CleanEquals (
  CharPtr p
)

{
  if (p == NULL) return;

  for (; *p != '\0'; p++) {
    if (*p == '\"') {
      *p = '\'';
    }
  }
}

static CharPtr GetPubTitle (
  FmtType format,
  PubdescPtr pdp,
  CitSubPtr csp
)

{
  CitArtPtr        cap;
  CitBookPtr       cbp;
  CitGenPtr        cgp;
  Char             ch;
  CitPatPtr        cpp;
  MedlineEntryPtr  mep;
  CharPtr          ptr;
  CharPtr          title = NULL;
  ValNodePtr       ttl = NULL;
  ValNodePtr       vnp;

  if (csp != NULL) {
    if (format == GENBANK_FMT || format == GENPEPT_FMT) {
      title = "Direct Submission";
      return StringSave (title);
    } else if (format == EMBL_FMT || format == EMBLPEPT_FMT) {
      return NULL;
    }
  }
  if (pdp == NULL) return NULL;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    switch (vnp->choice) {
      case PUB_Gen :
        cgp = (CitGenPtr) vnp->data.ptrvalue;
        if (cgp != NULL) {
          if (! StringHasNoText (cgp->title)) return StringSave (cgp->title);
          if (! StringHasNoText (cgp->cit)) {
            ptr = StringStr (cgp->cit, "Title=\"");
            if (ptr != NULL) {
              title = StringSave (ptr + 7);
              for (ptr = title; *ptr != '\0'; ptr++) {
                if (*ptr == '"') {
                  *ptr = '\0';
                  break;
                }
              }
              return title;
            }
          }
        }
        break;
      case PUB_Sub :
        csp = (CitSubPtr) vnp->data.ptrvalue;
        if (csp != NULL) {
          if (format == GENBANK_FMT || format == GENPEPT_FMT) {
            title = "Direct Submission";
            return StringSave (title);
          } else if (format == EMBL_FMT || format == EMBLPEPT_FMT) {
            return NULL;
          }
        }
        break;
      case PUB_Medline :
        mep = (MedlineEntryPtr) vnp->data.ptrvalue;
        if (mep != NULL) {
          cap = mep->cit;
          if (cap != NULL) {
            ttl = cap->title;
          }
        }
        break;
      case PUB_Article :
        cap = (CitArtPtr) vnp->data.ptrvalue;
        if (cap != NULL) {
          ttl = cap->title;
        }
        break;
      case PUB_Book :
      case PUB_Proc :
      case PUB_Man :
        cbp = (CitBookPtr) vnp->data.ptrvalue;
        if (cbp != NULL) {
          ttl = cbp->title;
          if (ttl != NULL) {
            title = (CharPtr) ttl->data.ptrvalue;
            if (! StringHasNoText (title)) {
              title = StringSave (title);
              if (StringLen (title) > 3) {
                ch = *title;
                if (IS_LOWER (ch)) {
                  *title = TO_UPPER (ch);
                }
                ptr = title;
                if (AllCaps (ptr)) {
                  for (ptr++; ptr != NULL && *ptr != '\0'; ptr++) {
                    ch = *ptr;
                    *ptr = TO_LOWER (ch);
                  }
                }
              }
              return title;
            }
          }
        }
        break;
      case PUB_Patent :
        cpp = (CitPatPtr) vnp->data.ptrvalue;
        if (cpp != NULL) {
          title = cpp->title;
          if (! StringHasNoText (title)) {
            return StringSave (title);
          }
        }
        break;
      default :
        break;
    }

    if (ttl != NULL) {
      title = (CharPtr) ttl->data.ptrvalue;
      if (! StringHasNoText (title)) {
        return StringSave (title);
      }
    }
  }

  return NULL;
}

static void CleanPubTitle (
  CharPtr title
)

{
  CharPtr  p;
  Boolean  remove_it;

  if (title == NULL) return;

  CleanEquals (title);

  for (p = title + StringLen (title) - 1; p > title + 2; p--) {
    if (*p == ' ') {
      *p = '\0';
    } else if (*p == '.') {
      remove_it = FALSE;
      if (p > title + 5) {
        if (*(p - 1) != '.' || *(p - 2) != '.') {
          remove_it = TRUE;
        }
      }
      if (remove_it) {
        *p = '\0';
      }
      break;
    } else {
      break;
    }
  }
}

/* !!! GetPubJournal needs to be implemented !!! */

/*
medline type page numbering is expanded (e.g., 125-35 -> 125-135,
F124-34 -> F124-F134, 12a-c -> 12a-12c).
If only one page is given, this is output without a dash.
Expanded numbering is validated to ensure that the
first number is smaller than or equal to the second and
that the first letter is less than or identical to the second
(i.e., a < c).  If the input is all letters (i.e., roman numerals)
this is not validated.

Return values:
 0 : valid page numbering.
-1 : invalid page numbering.
*/

#define MAX_PAGE_DIGITS 12

static Int2 FixPages (
  CharPtr out_pages,
  CharPtr in_pages
)

{
	Boolean dash=TRUE, first_alpha;
	Char firstbegin[MAX_PAGE_DIGITS];
	Char secondbegin[MAX_PAGE_DIGITS];
	Char firstend[MAX_PAGE_DIGITS];
	Char secondend[MAX_PAGE_DIGITS];
	Char temp[MAX_PAGE_DIGITS];
	CharPtr alphabegin, numbegin, alphaend, numend, ptr, in=in_pages;
	Int2 diff, index, retval=0;
	Int2 length_nb, length_ab, length_ne, length_ae;
	Int4 num1=0, num2=0;

	if (in_pages == NULL) return retval;

	while (*in != '\0')
	{			/* Check for digits in input*/
		if (IS_DIGIT(*in))
			break;
		in++;
	}

	if (*in == '\0' || (in != in_pages && *(in-1) == ' '))
	{		/* if all letters (i.e. roman numerals), put out. */
		out_pages = StringCpy(out_pages, in_pages);
		return retval;
	}

	in = in_pages;
	if (IS_DIGIT(*in))
	{			/* Do digits come first? */
		first_alpha = FALSE;
		index=0;
		while (IS_DIGIT(*in) || *in == ' ')
		{
			firstbegin[index] = *in;
			if (*in != ' ')
				index++;
			in++;
			if (*in == '-')
				break;

		}
		firstbegin[index] = '\0';
		index=0;
		if (*in != '-')
		{		/* After digits look for letters. */
			while (IS_ALPHA(*in)  || *in == ' ')
			{
				secondbegin[index] = *in;
				index++;
				in++;
				if (*in == '-')
					break;
			}
		}
		secondbegin[index] = '\0';
		if (*in == '-')		/* if dash is not present, note */
			in++;
		else
			dash=FALSE;
		index=0;
		while (IS_DIGIT(*in) || *in == ' ')
		{			/* Look for digits.	*/
			firstend[index] = *in;
			if (*in != ' ')
				index++;
			in++;
		}
		firstend[index] = '\0';
		index=0;
		if (*in != '\0')
		{			/* Look for letters again. */
			while (IS_ALPHA(*in)  || *in == ' ')
			{
				secondend[index] = *in;
				index++;
				in++;
			}
		}
		secondend[index] = '\0';
	}
	else
	{			/* Do letters come first? */
		first_alpha = TRUE;
		index=0;
		while (IS_ALPHA(*in) || *in == ' ')
		{
			firstbegin[index] = *in;
			index++;
			in++;
			if (*in == '-')
				break;
		}
		firstbegin[index] = '\0';
		index=0;
		if (*in != '-')
		{		/* After letters look for digits. 	*/
			while (IS_DIGIT(*in)  || *in == ' ')
			{
				secondbegin[index] = *in;
				if (*in != ' ')
					index++;
				in++;
				if (*in == '-')
					break;
			}
		}
		secondbegin[index] = '\0';
		if (*in == '-')		/* Note if dash is missing. */
			in++;
		else
			dash=FALSE;
		index=0;
		while (IS_ALPHA(*in) || *in == ' ')
		{		/* Look for letters again. */
			firstend[index] = *in;
			index++;
			in++;
		}
		firstend[index] = '\0';
		index=0;
		if (*in != '\0')
		{		/* Any digits here? */
			while (IS_DIGIT(*in)  || *in == ' ')
			{
				secondend[index] = *in;
				if (*in != ' ')
					index++;
				in++;
			}
		}
		secondend[index] = '\0';
	}

	if (first_alpha)
	{
		alphabegin = firstbegin;
		numbegin = secondbegin;
		alphaend = firstend;
		numend = secondend;
	}
	else
	{
		numbegin = firstbegin;
		alphabegin = secondbegin;
		numend = firstend;
		alphaend = secondend;
	}

	length_nb = StringLen(numbegin);
	length_ab = StringLen(alphabegin);
	length_ne = StringLen(numend);
	length_ae = StringLen(alphaend);

	/* If no dash, but second letters or numbers present, reject. */
	if (dash == FALSE)
	{
		if (length_ne != 0 || length_ae != 0)
			retval = -1;
	}
	/* Check for situations like "AAA-123" or "222-ABC". */
	if (dash == TRUE)
	{
		if (length_ne == 0 && length_ab == 0)
			retval = -1;
		else if (length_ae == 0 && length_nb == 0)
			retval = -1;
	}

	/* The following expands "F502-512" into "F502-F512" and
	checks, for entries like "12a-12c" that a > c.  "12aa-12ab",
	"125G-137A", "125-G137" would be rejected. */
	if (retval == 0)
	{
		if (length_ab > 0)
		{
			if (length_ae > 0) 	
			{ 
				if (StringCmp(alphabegin, alphaend) != 0)
				{
					if (length_ab != 1 || length_ae != 1)
						retval = -1;
					else if (*alphabegin > *alphaend)
						retval = -1;
				}
			}
			else
			{
				alphaend = alphabegin;
				length_ae = length_ab;
			}
		} 
		else if (length_ae > 0) 
			retval = -1;
	}

/* The following expands "125-37" into "125-137".	*/
	if (retval == 0)
	{
		if (length_nb > 0)
		{
			if (length_ne > 0)
			{
				diff = length_nb - length_ne;
				if (diff > 0)
				{
					index=0;
					while (numend[index] != '\0')
					{
						temp[index+diff] = numend[index];
						index++;
					}
					temp[index+diff] = numend[index];
					for (index=0; index<diff; index++)
						temp[index] = numbegin[index];
					index=0;
					while (temp[index] != '\0')
					{
						numend[index] = temp[index];
						index++;
					}
					numend[index] = temp[index];
				}
			}
			else
			{
				numend = numbegin;
				length_ne = length_nb;
			}
		
		}
		else if (length_ne > 0)
			retval = -1;
	/* Check that the first number is <= the second (expanded) number. */
		if (retval == 0)
		{
	/*		sscanf(numbegin, "%ld", &num_type);
			num1 = (Int4) num_type;
			sscanf(	numend, "%ld", &num_type);
			num2 = (Int4) num_type;
	*/
			num1 = (Int4) atol(numbegin);
			num2 = (Int4) atol(numend);
			if (num2 < num1)
				retval = -1;
		}
	}

	if (retval == -1)
	{
		out_pages = StringCpy(out_pages, in_pages);
	}
	else
	{
		ptr = out_pages;
	/* Place expanded and validated page numbers into "out_pages". */
		if (first_alpha)
		{
			while (*alphabegin != '\0')
			{
				*ptr = *alphabegin;
				alphabegin++;
				ptr++;
			}
			while (*numbegin != '\0')
			{
				*ptr = *numbegin;
				numbegin++;
				ptr++;
			}
			if (dash == TRUE)
			{
				*ptr = '-';
				ptr++;
				while (*alphaend != '\0')
				{
					*ptr = *alphaend;
					alphaend++;
					ptr++;
				}
				while (*numend != '\0')
				{
					*ptr = *numend;
					numend++;
					ptr++;
				}
			}
			*ptr = '\0';
		}
		else 
		{
			while (*numbegin != '\0')
			{
				*ptr = *numbegin;
				numbegin++;
				ptr++;
			}
			while (*alphabegin != '\0')
			{
				*ptr = *alphabegin;
				alphabegin++;
				ptr++;
			}
			if (dash == TRUE)
			{
				*ptr = '-';
				ptr++;
				while (*numend != '\0')
				{
					*ptr = *numend;
					numend++;
					ptr++;
				}
				while (*alphaend != '\0')
				{
					*ptr = *alphaend;
					alphaend++;
					ptr++;
				}
			}
			*ptr = '\0';
		}
	}
	return retval;
}

/* !!! still need to add StripParanthesis equivalent !!! */

static void DoSup (
  ValNodePtr PNTR head,
  CharPtr issue,
  CharPtr part_sup,
  CharPtr part_supi
)

{
	size_t   len;
	CharPtr  str;
	CharPtr  temp;

	len = StringLen (issue) + StringLen (part_sup) + StringLen (part_supi) + 25;
	str = MemNew (sizeof (Char) * len);
	if (str == NULL) return;
	temp = str;

	if (! StringHasNoText (part_sup)) {
		*temp = ' ';
		temp++;
		temp = StringMove (temp, part_sup);
	}
	if (StringHasNoText (issue) && StringHasNoText (part_supi)) {
		ValNodeCopyStr (head, 0, str);
		MemFree (str);
		return;
	}
	*temp = ' ';
	temp++;
	*temp = '(';
	temp++;
	if (! StringHasNoText (issue)) {
		temp = StringMove (temp, issue);
	}
	if (! StringHasNoText (part_supi)) {
		*temp = ' ';
		temp++;
		temp = StringMove (temp, part_supi);
	}
	*temp = ')';
	temp++;
	ValNodeCopyStr (head, 0, str);
	MemFree (str);
}

static CharPtr FormatCitJour (
  FmtType format,
  CitJourPtr cjp
)

{
  Char        buf [256];
  DatePtr     dp;
  ValNodePtr  head = NULL;
  ImprintPtr  imp;
  CharPtr     issue = NULL;
  Char        pages [128];
  CharPtr     part_sup = NULL;
  CharPtr     part_supi = NULL;
  CharPtr     rsult = NULL;
  CharPtr     title;
  ValNodePtr  ttl;
  CharPtr     volume;
  Char        year [8];

  if (cjp == NULL) return NULL;

  ttl = cjp->title;
  if (ttl == NULL) return NULL;

  imp = cjp->imp;
  if (imp == NULL) return NULL;

  dp = imp->date;
  year [0] = '\0';
  if (dp != NULL) {
    if (dp->data [0] == 1) {
      if (dp->data [1] != 0) {
        sprintf (year, " (%ld)", (long) (1900 + dp->data [1]));
      }
    } else {
      StringCpy (year, " (");
      StringNCat (year, dp->str, 4);
      StringCat (year, ")");
    }
  }

  if (imp->prepub == 1 || imp->prepub == 255) {
    sprintf (buf, "Unpublished (%s)", year);
    return StringSave (buf);
  }

  title = (CharPtr) ttl->data.ptrvalue;
  if (StringLen (title) < 3) return StringSave (".");

  ValNodeCopyStr (&head, 0, title);

  volume = imp->volume;
  if (format == GENBANK_FMT || format == GENPEPT_FMT) {
    issue = imp->issue;
    part_sup = imp->part_sup;
    part_supi = imp->part_supi;
  }
  pages [0] = '\0';
  FixPages (pages, imp->pages);

  if (! StringHasNoText (volume)) {
    AddValNodeString (&head, " ", volume, NULL);
  }

  if ((! StringHasNoText (volume)) || (! StringHasNoText (pages))) {
    DoSup (&head, issue, part_sup, part_supi);
  }

  if (format == GENBANK_FMT || format == GENPEPT_FMT) {
    if (! StringHasNoText (pages)) {
      AddValNodeString (&head, ", ", pages, NULL);
    }
  } else if (format == EMBL_FMT || format == EMBLPEPT_FMT) {
    if (! StringHasNoText (pages)) {
      AddValNodeString (&head, ":", pages, NULL);
    } else if (imp->prepub == 2 || (StringHasNoText (volume))) {
      ValNodeCopyStr (&head, 0, " 0:0-0");
    }
  }

  ValNodeCopyStr (&head, 0, year);

  if (format == GENBANK_FMT || format == GENPEPT_FMT) {
    if (imp->prepub == 2) {
      ValNodeCopyStr (&head, 0, " In press");
    }
  }

  rsult = MergeValNodeStrings (head);
  ValNodeFreeData (head);

  return rsult;
}

static CharPtr MakeAffilStr (
  AffilPtr afp
)

{
  ValNodePtr  head = NULL;
  CharPtr     prefix = "";
  CharPtr     rsult = NULL;

  if (afp == NULL) return NULL;

  if (! StringHasNoText (afp->affil)) {
    ValNodeCopyStr (&head, 0, afp->affil);
    prefix = ", ";
  }

  if (afp->choice == 2) {
    if (! StringHasNoText (afp->div)) {
      AddValNodeString (&head, prefix, afp->div, NULL);
      prefix = ", ";
    }
    if (! StringHasNoText (afp->street)) {
      AddValNodeString (&head, prefix, afp->street, NULL);
      prefix = ", ";
    }
    if (! StringHasNoText (afp->city)) {
      AddValNodeString (&head, prefix, afp->city, NULL);
      prefix = ", ";
    }
    if (! StringHasNoText (afp->sub)) {
      AddValNodeString (&head, prefix, afp->sub, NULL);
      prefix = ", ";
    }
    if (! StringHasNoText (afp->country)) {
      AddValNodeString (&head, prefix, afp->country, NULL);
      prefix = ", ";
    }
  }

  rsult = MergeValNodeStrings (head);
  ValNodeFreeData (head);

  return rsult;
}

static CharPtr FormatCitBook (
  FmtType format,
  CitBookPtr cbp
)

{
  AffilPtr     afp;
  AuthListPtr  alp;
  CharPtr      book_title = NULL;
  Char         buf [256];
  Char         ch;
  DatePtr      dp;
  ValNodePtr   head = NULL;
  ImprintPtr   imp;
  CharPtr      issue = NULL;
  ValNodePtr   names = NULL;
  Char         pages [128];
  CharPtr      part_sup = NULL;
  CharPtr      part_supi = NULL;
  CharPtr      rsult = NULL;
  CharPtr      str;
  CharPtr      title;
  ValNodePtr   ttl;
  ValNodePtr   vnp;
  CharPtr      volume;
  Char         year [8];

  if (cbp == NULL) return NULL;

  ttl = cbp->title;
  if (ttl == NULL) return NULL;

  imp = cbp->imp;
  if (imp == NULL) return NULL;

  dp = imp->date;
  year [0] = '\0';
  if (dp != NULL) {
    if (dp->data [0] == 1) {
      if (dp->data [1] != 0) {
        sprintf (year, "(%ld)", (long) (1900 + dp->data [1]));
      }
    } else {
      StringCpy (year, "(");
      StringNCat (year, dp->str, 4);
      StringCpy (year, ")");
    }
  }

  if (imp->prepub == 1 || imp->prepub == 255) {
    sprintf (buf, "Unpublished %s", year);
    return StringSave (buf);
  }

  title = (CharPtr) ttl->data.ptrvalue;
  if (StringLen (title) < 3) return StringSave (".");

  ValNodeCopyStr (&head, 0, "(in) ");

  alp = cbp->authors;
  if (alp != NULL) {
    str = GetAuthorsString (format, alp);
    if (str != NULL) {
      ValNodeCopyStr (&head, 0, str);
      names = alp->names;
      if (names != NULL) {
        if (names->next != NULL) {
          ValNodeCopyStr (&head, 0, " (Eds.);");
        } else {
          ValNodeCopyStr (&head, 0, " (Ed.);");
        }
      }
      ValNodeCopyStr (&head, 0, "\n");
    }
    MemFree (str);
  }

  book_title = StringSaveNoNull (title);
  vnp = ValNodeAddStr (&head, 0, book_title);
  if (book_title != NULL) {

    /* make book title all caps */

    title = book_title;
    ch = *title;
    while (ch != '\0') {
      *title = TO_UPPER (ch);
      title++;
      ch = *title;
    }
  }

  volume = imp->volume;
  if (format == GENBANK_FMT || format == GENPEPT_FMT) {
    issue = imp->issue;
    part_sup = imp->part_sup;
    part_supi = imp->part_supi;
  }
  pages [0] = '\0';
  FixPages (pages, imp->pages);

  if ((! StringHasNoText (volume)) && (StringCmp (volume, "0") != 0)) {
    AddValNodeString (&head, ", Vol. ", volume, NULL);
    DoSup (&head, issue, part_sup, part_supi);
  }

  if (! StringHasNoText (pages)) {
    AddValNodeString (&head, ": ", pages, NULL);
  }

  if (book_title != NULL) {
    ValNodeCopyStr (&head, 0, ";\n");
  }

  afp = imp->pub;
  if (afp != NULL) {
    str = MakeAffilStr (afp);
    if (str != NULL) {
      ValNodeCopyStr (&head, 0, str);
    ValNodeCopyStr (&head, 0, " ");
      MemFree (str);
    }
  }

  AddValNodeString (&head, NULL, year, NULL);

  if (format == GENBANK_FMT || format == GENPEPT_FMT) {
    if (imp->prepub == 2) {
      ValNodeCopyStr (&head, 0, " In press");
    }
  }

  rsult = MergeValNodeStrings (head);
  ValNodeFreeData (head);

  return rsult;
}

static CharPtr FormatThesis (
  FmtType format,
  CitBookPtr cbp
)

{
  AffilPtr     afp;
  Char         ch;
  DatePtr      dp;
  ValNodePtr   head = NULL;
  ImprintPtr   imp;
  CharPtr      ptr;
  CharPtr      rsult = NULL;
  CharPtr      str;
  CharPtr      suffix = NULL;
  Char         year [8];

  if (cbp == NULL) return NULL;
  if (cbp->othertype != 2 || cbp->let_type != 3) return NULL;

  imp = cbp->imp;
  if (imp == NULL) return NULL;

  dp = imp->date;
  year [0] = '\0';
  if (dp != NULL) {
    if (dp->data [0] == 1) {
      if (dp->data [1] != 0) {
        sprintf (year, "%ld", (long) (1900 + dp->data [1]));
      }
    } else {
      StringNCpy (year, dp->str, (size_t) 4);
      year [4] = '\0';
    }
  }

  AddValNodeString (&head, "Thesis (", year, ")");

  if (imp->prepub == 2) {
    suffix = ", In press";
  }

  str = NULL;
  afp = imp->pub;
  if (afp != NULL) {
    if (afp->choice == 1) {
      if (StringLen (afp->affil) > 7) {
        str = StringSave (afp->affil);
      }
    } else if (afp->choice == 2) {
      str = MakeAffilStr (afp);
    }
  }

  if (str != NULL) {

    /* convert double quotes to single quotes */

    ptr = str;
    ch = *ptr;
    while (ch != '\0') {
      if (ch == '\"') {
        *ptr = '\'';
      }
      ptr++;
      ch = *ptr;
    }
    AddValNodeString (&head, " ", str, suffix);
    MemFree (str);
  }

  rsult = MergeValNodeStrings (head);
  ValNodeFreeData (head);

  return rsult;
}

static CharPtr FormatCitArt (
  FmtType format,
  CitArtPtr cap
)

{
  CitBookPtr  cbp;
  CitJourPtr  cjp;
  CharPtr     rsult = NULL;

  if (cap == NULL) return NULL;

  switch (cap->from) {
    case 1 :
      cjp = (CitJourPtr) cap->fromptr;
      if (cjp != NULL) {
        rsult = FormatCitJour (format, cjp);
      }
      break;
    case 2 :
      cbp = (CitBookPtr) cap->fromptr;
      if (cbp != NULL) {
        rsult = FormatCitBook (format, cbp);
      }
      break;
    case 3 :
      cbp = (CitBookPtr) cap->fromptr;
      if (cbp != NULL) {
        rsult = FormatCitBook (format, cbp);
      }
      break;
    default :
      break;
  }

  return rsult;
}

static CharPtr FormatCitPat (
  FmtType format,
  CitPatPtr cpp
)

{
  AffilPtr     afp;
  AuthListPtr  alp;
  Char         date [40];
  ValNodePtr   head = NULL;
  CharPtr      prefix = NULL;
  CharPtr      rsult = NULL;
  CharPtr      suffix = NULL;

  if (cpp == NULL) return NULL;

  if (format == GENBANK_FMT || format == GENPEPT_FMT) {
    ValNodeCopyStr (&head, 0, "Patent: ");
    suffix = " ";
  } else if (format == EMBL_FMT || format == EMBLPEPT_FMT) {
    ValNodeCopyStr (&head, 0, "Patent number ");
  }

  if (! StringHasNoText (cpp->country)) {
    AddValNodeString (&head, NULL, cpp->country, suffix);
  }

  if (! StringHasNoText (cpp->number)) {
    ValNodeCopyStr (&head, 0, cpp->number);
  } else if (! StringHasNoText (cpp->app_number)) {
    AddValNodeString (&head, "(", cpp->app_number, ")");
  }

  if (! StringHasNoText (cpp->doc_type)) {
    AddValNodeString (&head, "-", cpp->doc_type, NULL);
  }

  /* !!! pat_seqid test not yet implemented !!! */

  ValNodeCopyStr (&head, 0, " ");

  date [0] = '\0';
  if (cpp->date_issue != NULL) {
    DateToGB (date, cpp->date_issue);
  } else if (cpp->app_date != NULL) {
    DateToGB (date, cpp->app_date);
  }
  if (! StringHasNoText (date)) {
    ValNodeCopyStr (&head, 0, date);
  }

  if (format == GENBANK_FMT || format == GENPEPT_FMT) {
    ValNodeCopyStr (&head, 0, ";");
  } else if (format == EMBL_FMT || format == EMBLPEPT_FMT) {
    ValNodeCopyStr (&head, 0, ".");
  }

  alp = cpp->authors;
  if (alp != NULL) {
    afp = alp->affil;
    if (afp != NULL) {
      suffix = NULL;
      if (afp->choice == 2) {
        suffix = ";";
      }
      if (! StringHasNoText (afp->affil)) {
        AddValNodeString (&head, "\n", afp->affil, suffix);
      }
      if (! StringHasNoText (afp->street)) {
        AddValNodeString (&head, "\n", afp->street, ";");
      }
      if (! StringHasNoText (afp->div)) {
        AddValNodeString (&head, "\n", afp->div, ";");
      }
      if (! StringHasNoText (afp->city)) {
        AddValNodeString (&head, "\n", afp->city, ", ");
      }
      if (! StringHasNoText (afp->sub)) {
        ValNodeCopyStr (&head, 0, afp->sub);
      }
      if (! StringHasNoText (afp->country)) {
        AddValNodeString (&head, ";\n", afp->country, ";");
      }
    }
  }

  rsult = MergeValNodeStrings (head);
  ValNodeFreeData (head);

  return rsult;
}

static CharPtr FormatCitGen (
  FmtType format,
  CitGenPtr cgp
)

{
  Char        ch;
  DatePtr     dp;
  ValNodePtr  head = NULL;
  CharPtr     inpress = NULL;
  CharPtr     journal = NULL;
  Char        pages [128];
  CharPtr     ptr;
  CharPtr     rsult = NULL;
  Char        year [8];

  if (cgp == NULL) return NULL;

  if (StringNICmp (cgp->cit, "unpublished", 11) == 0) {
    return StringSave ("Unpublished");
  }

  year [0] = '\0';
  dp = cgp->date;
  if (dp != NULL) {
    if (dp->data [0] == 1) {
      if (dp->data [1] != 0) {
        sprintf (year, " (%ld)", (long) (1900 + dp->data [1]));
      }
    } else {
      StringCpy (year, " (");
      StringNCat (year, dp->str, 4);
      StringCat (year, ")");
    }
  }

  pages [0] = '\0';
  if (cgp->pages != NULL) {
    FixPages (pages, cgp->pages);
  }

  if (cgp->journal != NULL) {
    journal = (CharPtr) cgp->journal->data.ptrvalue;
  }
  if (cgp->cit != NULL) {
    ptr = StringStr (cgp->cit, "Journal=\"");
    if (ptr != NULL) {
      journal = ptr + 9;
    } else {
      if (StringNICmp (cgp->cit, "submitted", 8) == 0 ||
          StringNICmp (cgp->cit, "in press", 8) == 0 ||
          StringNICmp (cgp->cit, "to be published", 15) == 0) {
        inpress = cgp->cit;
      }
    }
  }
  if (journal != NULL) {
    journal = StringSave (journal);
    for (ptr = journal, ch = *ptr; ch != '\0'; ptr++, ch = *ptr) {
      if (ch == '=' || ch == '\"') {
        *ptr = '\0';
      }
    }
    ValNodeAddStr (&head, 0, journal);
  }

  if (! StringHasNoText (inpress)) {
    AddValNodeString (&head, " ", inpress, NULL);
  }

  if (! StringHasNoText (cgp->volume)) {
    AddValNodeString (&head, " ", cgp->volume, NULL);
  }

  if (! StringHasNoText (pages)) {
    if (format == GENBANK_FMT) {
      AddValNodeString (&head, ", ", pages, NULL);
    } else if (format == EMBL_FMT) {
      AddValNodeString (&head, ":", pages, NULL);
    }
  }

  if (! StringHasNoText (year)) {
    AddValNodeString (&head, NULL, year, NULL);
  }

  rsult = MergeValNodeStrings (head);
  ValNodeFreeData (head);

  return rsult;
}

static CharPtr GetAffil (
  AffilPtr afp
)

{
	Boolean need_comma=FALSE;
	CharPtr string=NULL, temp, ptr;
	Char ch;
	Int2 aflen=15;

	if (afp) {
		if (afp -> choice == 1){
			if (afp -> affil){
				aflen += StringLen(afp -> affil);
			}
		}else if (afp -> choice == 2){
			aflen += StringLen (afp -> affil) + 
			StringLen (afp -> div) + 
			StringLen (afp -> city) + 
			StringLen (afp -> sub) + 
			StringLen (afp -> street) + 
			StringLen (afp -> country) + StringLen(afp->postal_code);
		}

		temp = string = MemNew(aflen);

		if ( afp -> choice == 1){
			 if (afp -> affil){
				ptr = afp->affil;
				while ((*temp = *ptr) != '\0')
				{
					temp++; ptr++;
				}
			 }
		}else if (afp -> choice == 2){

			if( afp -> div) { 
				if (need_comma)
				{
					*temp = ','; temp++;
					*temp = ' '; temp++;
				}
				ptr = afp->div;
				while ((*temp = *ptr) != '\0')
				{
					temp++; ptr++;
				}
				need_comma = TRUE;
			}

			if(afp -> affil) { 
				if (need_comma)
				{
					*temp = ','; temp++;
					*temp = ' '; temp++;
				}
				ptr = afp->affil;
				while ((*temp = *ptr) != '\0')
				{
					temp++; ptr++;
				}
				need_comma = TRUE;
			}

			if(afp -> street) { 
				if (need_comma)
				{
					*temp = ','; temp++;
					*temp = ' '; temp++;
				}
				ptr = afp->street;
				while ((*temp = *ptr) != '\0')
				{
					temp++; ptr++;
				}
				need_comma = TRUE;
			}

			if( afp -> city) { 
				if (need_comma)
				{
					*temp = ','; temp++;
					*temp = ' '; temp++;
				}
				ptr = afp->city;
				while ((*temp = *ptr) != '\0')
				{
					temp++; ptr++;
				}
				need_comma = TRUE;
			}

			if( afp -> sub) { 
				if (need_comma)
				{
					*temp = ','; temp++;
					*temp = ' '; temp++;
				}
				ptr = afp->sub;
				while ((*temp = *ptr) != '\0')
				{
					temp++; ptr++;
				}
				need_comma = TRUE;
			}

			if( afp -> postal_code){
				*temp = ' '; 
				temp++;
				ptr = afp->postal_code;
				while ((*temp = *ptr) != '\0')
				{
					temp++; ptr++;
				}
			}	

			if( afp -> country){
				if (need_comma)
				{
					*temp = ','; temp++;
					*temp = ' '; temp++;
				}
				ptr = afp->country;
				while ((*temp = *ptr) != '\0')
				{
					temp++; ptr++;
				}
				need_comma = TRUE;
			}	
		}
		temp++;
		*temp = '\0';
	}

    /* convert double quotes to single quotes */

    ptr = string;
    ch = *ptr;
    while (ch != '\0') {
      if (ch == '\"') {
        *ptr = '\'';
      }
      ptr++;
      ch = *ptr;
    }

	return string;
}

static CharPtr FormatCitSub (
  FmtType format,
  CitSubPtr csp
)

{
  CharPtr      affil;
  AffilPtr     afp;
  AuthListPtr  alp;
  Char         buf [256];
  Char         date [40];
  ValNodePtr   head = NULL;
  CharPtr      rsult = NULL;

  if (csp == NULL) return NULL;

  date [0] = '\0';
  if (csp->date != NULL) {
    DateToGB (date, csp->date);
  }
  if (StringHasNoText (date)) {
    StringCpy (date, "??-???-????");
  }

  sprintf (buf, "Submitted (%s)", date);
  ValNodeCopyStr (&head, 0, buf);

  alp = csp->authors;
  if (alp != NULL) {
    afp = alp->affil;
    if (afp != NULL) {
      affil = GetAffil (afp);
      if (format == EMBL_FMT || format == EMBLPEPT_FMT) {
        if (StringNCmp(affil, " to the EMBL/GenBank/DDBJ databases.", 36) != 0) {
          ValNodeCopyStr (&head, 0, " to the EMBL/GenBank/DDBJ databases.\n");
        } else {
          ValNodeCopyStr (&head, 0, " ");
        }
      } else {
        ValNodeCopyStr (&head, 0, " ");
      }
      ValNodeCopyStr (&head, 0, affil);
      MemFree (affil);
    }
  }

  rsult = MergeValNodeStrings (head);
  ValNodeFreeData (head);

  return rsult;
}

static CharPtr GetPubJournal (
  FmtType format,
  PubdescPtr pdp,
  CitSubPtr csp
)

{
  CitArtPtr        cap;
  CitBookPtr       cbp;
  CitGenPtr        cgp;
  CitPatPtr        cpp;
  CharPtr          journal = NULL;
  MedlineEntryPtr  mep;
  ValNodePtr       vnp;

  if (csp != NULL) {
    return FormatCitSub (format, csp);
  }
  if (pdp == NULL) return NULL;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    switch (vnp->choice) {
      case PUB_Gen :
        cgp = (CitGenPtr) vnp->data.ptrvalue;
        if (cgp != NULL) {
          journal = FormatCitGen (format, cgp);
        }
        break;
      case PUB_Sub :
        csp = (CitSubPtr) vnp->data.ptrvalue;
        if (csp != NULL) {
          journal = FormatCitSub (format, csp);
        }
        break;
      case PUB_Medline :
        mep = (MedlineEntryPtr) vnp->data.ptrvalue;
        if (mep != NULL) {
          cap = mep->cit;
          if (cap != NULL) {
            journal = FormatCitArt (format, cap);
          }
        }
        break;
      case PUB_Article :
        cap = (CitArtPtr) vnp->data.ptrvalue;
        if (cap != NULL) {
          journal = FormatCitArt (format, cap);
        }
        break;
      case PUB_Book :
      case PUB_Proc :
        cbp = (CitBookPtr) vnp->data.ptrvalue;
        if (cbp != NULL) {
          journal = FormatCitBook (format, cbp);
        }
        break;
      case PUB_Man :
        cbp = (CitBookPtr) vnp->data.ptrvalue;
        if (cbp != NULL) {
          journal = FormatThesis (format, cbp);
        }
        break;
      case PUB_Patent :
        cpp = (CitPatPtr) vnp->data.ptrvalue;
        if (cpp != NULL) {
          journal = FormatCitPat (format, cpp);
        }
        break;
      default :
        break;
    }

    if (journal != NULL) return journal;
  }

  return NULL;
}

static Int4 GetMuid (
  PubdescPtr pdp
)

{
  MedlineEntryPtr  mep;
  ValNodePtr       vnp;

  if (pdp == NULL) return 0;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    switch (vnp->choice) {
      case PUB_Medline :
        mep = (MedlineEntryPtr) vnp->data.ptrvalue;
        if (mep != NULL) {
          return mep->uid;
        }
        break;
      case PUB_Muid :
        return vnp->data.intvalue;
        break;
      default :
        break;
    }
  }

  return 0;
}

static Int4 GetPmid (
  PubdescPtr pdp
)

{
  MedlineEntryPtr  mep;
  ValNodePtr       vnp;

  if (pdp == NULL) return 0;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    switch (vnp->choice) {
      case PUB_Medline :
        mep = (MedlineEntryPtr) vnp->data.ptrvalue;
        if (mep != NULL) {
          return mep->pmid;
        }
        break;
      case PUB_PMid :
        return vnp->data.intvalue;
        break;
      default :
        break;
    }
  }

  return 0;
}

static CharPtr remarksText [] = {
  "full automatic", "full staff_review", "full staff_entry",
  "simple staff_review", "simple staff_entry", "simple automatic",
  "unannotated automatic", "unannotated staff_review", "unannotated staff_entry",
  NULL
};

static CharPtr FormatReferenceBlock (
  Asn2gbFormatPtr afp,
  BaseBlockPtr bbp
)

{
  Asn2gbJobPtr       ajp;
  AuthListPtr        alp;
  Asn2gbSectionPtr   asp;
  BioseqPtr          bsp;
  Char               buf [150];
  CitArtPtr          cap;
  CitJourPtr         cjp;
  CitRetractPtr      crp;
  CitSubPtr          csp = NULL;
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext;
  Int4               gibbsq;
  Int2               i;
  ImprintPtr         imp;
  IntRefBlockPtr     irp;
  SeqLocPtr          loc = NULL;
  Int4               muid;
  Boolean            needsPeriod;
  SeqLocPtr          nextslp;
  Boolean            notFound;
  ObjMgrDataPtr      omdp;
  PubdescPtr         pdp = NULL;
  Int4               pmid;
  CharPtr            prefix;
  ReferenceBlockPtr  rbp;
  SubmitBlockPtr     sbp;
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp = NULL;
  SeqIdPtr           sip;
  SeqLocPtr          slp;
  SeqSubmitPtr       ssp;
  Int4               start;
  Int4               stop;
  CharPtr            str;
  CharPtr            suffix = NULL;
  Boolean            trailingPeriod;
  ValNodePtr         vnp;

  if (afp == NULL || bbp == NULL) return NULL;
  rbp = (ReferenceBlockPtr) bbp;
  ajp = afp->ajp;
  if (ajp == NULL) return NULL;
  asp = afp->asp;
  if (asp == NULL) return NULL;
  bsp = asp->bsp;
  if (bsp == NULL) return NULL;

  if (! StringHasNoText (rbp->string)) return StringSave (rbp->string);

  /* could be descriptor, feature, or submit block citation */

  if (rbp->itemtype == OBJ_SEQDESC) {

    sdp = SeqMgrGetDesiredDescriptor (rbp->entityID, NULL, rbp->itemID, 0, NULL, &dcontext);
    if (sdp != NULL && dcontext.seqdesctype == Seq_descr_pub) {
      pdp = (PubdescPtr) sdp->data.ptrvalue;
    }

  } else if (rbp->itemtype == OBJ_SEQFEAT) {

    sfp = SeqMgrGetDesiredFeature (rbp->entityID, NULL, rbp->itemID, 0, NULL, &fcontext);
    if (sfp != NULL && fcontext.seqfeattype == SEQFEAT_PUB) {
      pdp = (PubdescPtr) sfp->data.value.ptrvalue;
    }

  } else if (rbp->itemtype == OBJ_SEQSUB_CIT) {

    omdp = ObjMgrGetData (rbp->entityID);
    if (omdp != NULL && omdp->datatype == OBJ_SEQSUB) {
      ssp = (SeqSubmitPtr) omdp->dataptr;
      if (ssp != NULL && ssp->datatype == 1) {
        sbp = ssp->sub;
        if (sbp != NULL) {
          csp = sbp->cit;
        }
      }
    }
  }

  if (pdp == NULL && csp == NULL) return NULL;

  /* print serial number */

  gb_StartPrint (afp->format, TRUE, 0, 12, "REFERENCE", 13, 5, 5, "RN", TRUE);

  if (afp->format ==
   GENBANK_FMT || afp->format == GENPEPT_FMT) {
    sprintf (buf, "%d", (int) rbp->serial);
  } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {
    sprintf (buf, "[%d]", (int) rbp->serial);
  }

  gb_AddString (NULL, buf, NULL, FALSE, FALSE, FALSE);

  /* print base range */

  if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
    TabToColumn (16);
  } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {
    ff_EndPrint ();

    ff_StartPrint (5, 5, ASN2FF_EMBL_MAX, "RP");
  }

  if (rbp->sites) {

    gb_AddString (NULL, "(sites)", NULL, FALSE, FALSE, FALSE);
    ff_EndPrint ();

  } else {

    if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
      TabToColumn (16);
      if (afp->format == GENBANK_FMT) {
        gb_AddString (NULL, "(bases ", NULL, FALSE, FALSE, FALSE);
      } else {
        gb_AddString (NULL, "(residues ", NULL, FALSE, FALSE, FALSE);
      }
    }

    irp = (IntRefBlockPtr) rbp;
    loc = irp->loc;

    if (loc != NULL) {
      if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
        suffix = "; ";
      } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {
        suffix = ", ";
      }

      slp = SeqLocFindNext (loc, NULL);
      while (slp != NULL) {
        nextslp = SeqLocFindNext (loc, slp);
        start = SeqLocStart (slp) + 1;
        stop = SeqLocStop (slp) + 1;
        if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
          sprintf (buf, "%ld to %ld", (long) start, (long) stop);
        } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {
          sprintf (buf, "%ld-%ld", (long) start, (long) stop);
        }
        if (nextslp == NULL) {
          suffix = NULL;
        }
        gb_AddString (NULL, buf, suffix, FALSE, FALSE, FALSE);
        slp = nextslp;
      }

    } else {

      /* code still used for ssp->cit */

      start = 1;
      stop = bsp->length;
      if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
        sprintf (buf, "%ld to %ld", (long) start, (long) stop);
      } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {
        sprintf (buf, "%ld-%ld", (long) start, (long) stop);
      }
      gb_AddString (NULL, buf, NULL, FALSE, FALSE, FALSE);
    }

    if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
      gb_AddString (NULL, ")", NULL, FALSE, FALSE, FALSE);
    }
    ff_EndPrint ();
  }

  /* print author list */

  gb_StartPrint (afp->format, FALSE, 2, 12, "AUTHORS", 13, 5, 5, "RA", FALSE);

  str = NULL;

  alp = GetAuthListPtr (pdp, csp);
  if (alp != NULL) {
    str = GetAuthorsString (afp->format, alp);
  }

  if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
    suffix = NULL;
    trailingPeriod = TRUE;
  } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {
    suffix = ";";
    trailingPeriod = FALSE;
  }

  gb_AddString (NULL, str, suffix, trailingPeriod, FALSE, FALSE);

  MemFree (str);
  ff_EndPrint ();

  /* print title */

  str = GetPubTitle (afp->format, pdp, csp);
  CleanPubTitle (str);
  StrStripSpaces (str);

  if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
    prefix = NULL;
    suffix = NULL;
  } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {
    if (str != NULL) {
      prefix = "\"";
      suffix = "\";";
    } else {
      prefix = NULL;
      suffix = ";";
    }
  }

  if (! StringHasNoText (str)) {
    gb_StartPrint (afp->format, FALSE, 2, 12, "TITLE", 13, 5, 5, "RT", FALSE);

    gb_AddString (prefix, str, suffix, FALSE, FALSE, FALSE);

    ff_EndPrint ();
  }

  MemFree (str);

  /* print journal */

  gb_StartPrint (afp->format, FALSE, 2, 12, "JOURNAL", 13, 5, 5, "RL", FALSE);

  str = GetPubJournal (afp->format, pdp, csp);
  if (str == NULL) {
    str = StringSave ("Unpublished");
  }
  StrStripSpaces (str);

  if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
    needsPeriod = FALSE;
  } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {
    needsPeriod = TRUE;
  }

  gb_AddString (NULL, str, NULL, needsPeriod, FALSE, FALSE);

  MemFree (str);
  ff_EndPrint ();

  /* print muid */

  muid = GetMuid (pdp);
  if (muid > 0) {
    gb_StartPrint (afp->format, FALSE, 2, 12, "MEDLINE", 13, 5, 5, "RX", FALSE);

    if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
      www_muid (muid);
    } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {
      sprintf (buf, "MEDLINE; %ld.", (long) muid);
      gb_AddString (NULL, buf, NULL, FALSE, FALSE, FALSE);
    }

    ff_EndPrint ();
  }

  if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
    pmid = GetPmid (pdp);
    if (pmid > 0) {
      gb_StartPrint (afp->format, FALSE, 3, 12, "PUBMED", 13, 5, 5, "RX", FALSE);

      www_muid (pmid);

      ff_EndPrint ();
    }
  }

  if (pdp == NULL) return gb_MergeString (FALSE);

  /* !!! remainder of fields are only for GenBank !!! */

  if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {

    prefix = "REMARK";

    if (pdp->comment != NULL) {
      for (i = 0, notFound = TRUE; notFound && remarksText [i] != NULL; i++) {
        if (StringCmp (pdp->comment, remarksText [i]) == 0) {
          notFound = FALSE;
        }
      }
      if (notFound) {
        gb_StartPrint (afp->format, FALSE, 2, 12, prefix, 13, 5, 5, NULL, FALSE);
        gb_AddString (NULL, pdp->comment, NULL, FALSE, FALSE, TRUE);
        ff_EndPrint ();
        prefix = NULL;
      }
    }

    gibbsq = 0;
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      if (sip->choice == SEQID_GIBBSQ) {
        gibbsq = sip->data.intvalue;
      }
    }
    if (gibbsq > 0) {
      sprintf (buf, "GenBank staff at the National Library of Medicine created this entry [NCBI gibbsq %ld] from the original journal article.", (long) gibbsq);
      gb_StartPrint (afp->format, FALSE, 2, 12, prefix, 13, 5, 5, NULL, FALSE);
      gb_AddString (NULL, buf, NULL, FALSE, FALSE, TRUE);
      ff_EndPrint ();
      prefix = NULL;

      /* gibbsq comment section (fields may be copied from degenerate pubdesc) */

      str = pdp->fig;
      if (StringHasNoText (str)) {
        str = irp->fig;
      }
      if (! StringHasNoText (str)) {
        sprintf (buf, "This sequence comes from %s", str);
        gb_StartPrint (afp->format, FALSE, 2, 12, prefix, 13, 5, 5, NULL, FALSE);
        gb_AddString (NULL, buf, NULL, TRUE, TRUE, TRUE);
        ff_EndPrint ();
        prefix = NULL;
      }

      if (pdp->poly_a || irp->poly_a) {
        gb_StartPrint (afp->format, FALSE, 2, 12, prefix, 13, 5, 5, NULL, FALSE);
        gb_AddString (NULL, "Polyadenylate residues occurring in the figure were omitted from the sequence.", NULL, TRUE, TRUE, TRUE);
        ff_EndPrint ();
        prefix = NULL;
      }

      str = pdp->maploc;
      if (StringHasNoText (str)) {
        str = irp->maploc;
      }
      if (! StringHasNoText (str)) {
        sprintf (buf, "Map location: %s", str);
        gb_StartPrint (afp->format, FALSE, 2, 12, prefix, 13, 5, 5, NULL, FALSE);
        gb_AddString (NULL, buf, NULL, TRUE, TRUE, TRUE);
        ff_EndPrint ();
        prefix = NULL;
      }

    }

    for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
      if (vnp->choice == PUB_Article) {
        cap = (CitArtPtr) vnp->data.ptrvalue;
        if (cap != NULL && cap->from == 1) {
          cjp = (CitJourPtr) cap->fromptr;
          if (cjp != NULL) {
            imp = cjp->imp;
            if (imp != NULL) {
              crp = imp->retract;
              if (crp != NULL && crp->type == 3) {
                gb_StartPrint (afp->format, FALSE, 2, 12, prefix, 13, 5, 5, NULL, FALSE);
                gb_AddString (NULL, "Erratum:", NULL, FALSE, FALSE, FALSE);
                gb_AddString ("[", crp->exp, "]", FALSE, TRUE, TRUE);
                ff_EndPrint ();
                prefix = NULL;
              }
            }
          }
        }
      } else if (vnp->choice == PUB_Sub) {
        csp = (CitSubPtr) vnp->data.ptrvalue;
        if (csp != NULL) {
          if (! StringHasNoText (csp->descr)) {
            gb_StartPrint (afp->format, FALSE, 2, 12, prefix, 13, 5, 5, NULL, FALSE);
            gb_AddString (NULL, csp->descr, NULL, FALSE, TRUE, TRUE);
            ff_EndPrint ();
            prefix = NULL;
          }
        }
      }
    }

  }

  return gb_MergeString (FALSE);
}

static CharPtr FormatCommentBlock (
  Asn2gbFormatPtr afp,
  BaseBlockPtr bbp
)

{
  CommentBlockPtr    cbp;
  CharPtr            db;
  DbtagPtr           dbt;
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext;
  ObjectIdPtr        oip;
  CharPtr            prefix;
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp;
  Char               sfx [32];
  CharPtr            suffix;
  CharPtr            title;

  if (afp == NULL || bbp == NULL) return NULL;
  cbp = (CommentBlockPtr) bbp;

  /* some comments are allocated (along with possible first COMMENT label) */

  if (! StringHasNoText (bbp->string)) return StringSave (bbp->string);

  title = NULL;
  prefix = NULL;
  suffix = NULL;
  sfx [0] = '\0';

  if (bbp->itemtype == OBJ_SEQDESC) {

    /* usually should reference comment, maploc, or region descriptor IDs */

    sdp = SeqMgrGetDesiredDescriptor (bbp->entityID, NULL, bbp->itemID, 0, NULL, &dcontext);
    if (sdp != NULL) {

      if (dcontext.seqdesctype == Seq_descr_comment) {

        title = (CharPtr) sdp->data.ptrvalue;

      } else if (dcontext.seqdesctype == Seq_descr_maploc) {

        dbt = (DbtagPtr) sdp->data.ptrvalue;
        if (dbt != NULL) {
          db = dbt->db;
          oip = dbt->tag;
          if (oip != NULL) {
            if (oip->str != NULL) {

              title = oip->str;
              prefix = ("Map location: ");

            } else if (db != NULL && oip->id != 0) {

              title = db;
              prefix = ("Map location: (Database ");
              sprintf (sfx, "; id # %ld).", (long) oip->id);
              suffix = sfx;

            }
          }
        }

      } else if (dcontext.seqdesctype == Seq_descr_region) {

        title = (CharPtr) sdp->data.ptrvalue;
        prefix = "Region: ";

      }
    }

  } else if (bbp->itemtype == OBJ_SEQFEAT) {

    /* also have to deal with comment feature across entire sequence */

    sfp = SeqMgrGetDesiredFeature (bbp->entityID, NULL, bbp->itemID, 0, NULL, &fcontext);
    if (sfp != NULL && fcontext.seqfeattype == SEQFEAT_COMMENT) {

      title = sfp->comment;
    }
  }

  if (title == NULL) return NULL;

  if (cbp->first) {
    gb_StartPrint (afp->format, TRUE, 0, 12, "COMMENT", 13, 5, 5, "CC", TRUE);
  } else {
    gb_StartPrint (afp->format, TRUE, 0, 12, NULL, 13, 5, 5, "CC", FALSE);
  }

  gb_AddString (prefix, title, suffix, TRUE, TRUE, TRUE);

  return gb_MergeString (TRUE);
}

/* format features section */

static Boolean is_real_id (
  SeqIdPtr sip,
  SeqIdPtr this_sip
)

{
  BioseqPtr  bsp;

  if (sip == NULL || this_sip == NULL) return FALSE;

  if (! SeqIdIn (sip, this_sip)) {
    bsp = BioseqFind (sip);
    if (bsp == NULL) return TRUE;  /* ??? */
    if (bsp->repr == Seq_repr_virtual) return FALSE;
  }

  return TRUE;
}

static Boolean FlatVirtLoc (
  BioseqPtr bsp,
  SeqLocPtr location
)

{
  SeqIntPtr  sintp;
  SeqIdPtr   sip;
  SeqPntPtr  spp;

  if (bsp == NULL || location == NULL) return FALSE;

  switch (location->choice) {
    case SEQLOC_WHOLE :
      sip = (SeqIdPtr) location->data.ptrvalue;
      if (sip == NULL) return TRUE;
      if (! is_real_id (sip, bsp->id)) return TRUE;
      break;
    case SEQLOC_INT :
      sintp = (SeqIntPtr) location->data.ptrvalue;
      if (sintp == NULL) return TRUE;
      sip = sintp->id;
      if (sip == NULL) return TRUE;
      if (! is_real_id (sip, bsp->id)) return TRUE;
      break;
    case SEQLOC_PNT :
      spp = (SeqPntPtr) location->data.ptrvalue;
      if (spp == NULL) return TRUE;
      sip = spp->id;
      if (sip == NULL) return TRUE;
      if (! is_real_id (sip, bsp->id)) return TRUE;
      break;
    default :
      break;
  }

  return FALSE;
}

static Uint1    order [18];
static Boolean  order_initialized = FALSE;

static CharPtr lim_str [5] = { "", ">","<", ">", "<" };

static void FlatLocSeqId (
  ValNodePtr PNTR head,
  SeqIdPtr sip
)

{
  BioseqPtr    bsp;
  Char         buf [40];
  ObjectIdPtr  oip;
  SeqIdPtr     use_id = NULL;
  Boolean      was_lock = FALSE;

  if (head == NULL || sip == NULL) return;

  bsp = BioseqFind (sip);
  if (bsp == NULL) {
    bsp = BioseqLockById (sip);
    was_lock = TRUE;
  }
  if (bsp != NULL) {
    use_id = SeqIdSelect (bsp->id, order, 18);
  } else if (sip->choice == SEQID_GI) {
    use_id = GetSeqIdForGI (sip->data.intvalue);
  }
  if (use_id != NULL) {
    SeqIdWrite (use_id, buf, PRINTID_TEXTID_ACC_VER, sizeof (buf) - 1);
  } else if (sip->choice == SEQID_GI) {
    SeqIdWrite (sip, buf, PRINTID_FASTA_LONG, sizeof (buf) - 1);
  }
  if (was_lock) {
    BioseqUnlock (bsp);
  }
  if (StringHasNoText (buf)) {
    StringCpy (buf, "?00000");
    if (use_id != NULL && use_id->choice == SEQID_LOCAL) {
      oip = (ObjectIdPtr) use_id->data.ptrvalue;
      if (oip != NULL && (! StringHasNoText (oip->str))) {
        StringNCpy_0 (buf, oip->str, 13);
      }
    }
  }
  AddValNodeString (head, NULL, buf, ":");
}

static void FlatLocCaret (
  ValNodePtr PNTR head,
  SeqIdPtr sip,
  SeqIdPtr this_sip,
  Int4 point,
  IntFuzzPtr fuzz
)

{
  Char   buf [128];
  Uint1  index;

  if (head == NULL) return;

  if (sip != NULL && (! SeqIdIn (sip, this_sip))) {
    FlatLocSeqId (head, sip);
  }

  buf [0] = '\0';
  point++;

  if (fuzz != NULL) {
    switch (fuzz->choice) {
      case 1 :
        sprintf (buf, "(%ld.%ld)..(%ld.%ld)",
                 (long) (point - fuzz->a),
                 (long) point,
                 (long) point,
                 (long) (point + fuzz->a));
        break;
      case 2 :
        sprintf (buf, "%ld^%ld",
                 (long) (1 + fuzz->b),
                 (long) (1 + fuzz->a));
        break;
      case 3 :
        sprintf (buf, "%ld^%ld",
                 (long) (point - point * ((double) fuzz->a / 1000.0)),
                 (long) (point + point * ((double) fuzz->a / 1000.0)));
        break;
      case 4 :
        if (point > 1) {
          sprintf (buf, "%ld^%ld", (long) (point - 1), (long) point);
        } else {
          index = (Uint1) fuzz->a;
          if (index > 4) {
            index = 0;
          }
          sprintf (buf, "%s%ld", lim_str [index], (long) point);
        }
        break;
      default :
        sprintf (buf, "%ld", (long) point);
        break;
    }
  } else {
    sprintf (buf, "%ld", (long) point);
  }

  ValNodeCopyStr (head, 0, buf);
}

static void FlatLocPoint (
  ValNodePtr PNTR head,
  SeqIdPtr sip,
  SeqIdPtr this_sip,
  Int4 point,
  IntFuzzPtr fuzz
)

{
  Char   buf [128];
  Uint1  index;

  if (head == NULL) return;

  if (sip != NULL && (! SeqIdIn (sip, this_sip))) {
    FlatLocSeqId (head, sip);
  }

  buf [0] = '\0';
  point++;

  if (fuzz != NULL) {
    switch (fuzz->choice) {
      case 1 :
        sprintf (buf, "(%ld.%ld)",
                 (long) (point - fuzz->a),
                 (long) (point + fuzz->a));
        break;
      case 2 :
        sprintf (buf, "(%ld.%ld)",
                 (long) (1 + fuzz->b),
                 (long) (1 + fuzz->a));
        break;
      case 3 :
        sprintf (buf, "(%ld.%ld)",
                 (long) (point - point * ((double) fuzz->a / 1000.0)),
                 (long) (point + point * ((double) fuzz->a / 1000.0)));
        break;
      case 4 :
        index = (Uint1) fuzz->a;
        if (index > 4) {
          index = 0;
        }
        sprintf (buf, "%s%ld", lim_str [index], (long) point);
        break;
      default :
        sprintf (buf, "%ld", (long) point);
        break;
    }
  } else {
    sprintf (buf, "%ld", (long) point);
  }

  ValNodeCopyStr (head, 0, buf);
}

static void FlatLocElement (
  ValNodePtr PNTR head,
  BioseqPtr bsp,
  SeqLocPtr location
)

{
  Boolean     minus_strand = FALSE;
  SeqBondPtr  sbp;
  SeqIntPtr   sintp;
  SeqIdPtr    sip;
  SeqPntPtr   spp;
  BioseqPtr   wholebsp;

  if (head == NULL || bsp == NULL || location == NULL) return;

  switch (location->choice) {
    case SEQLOC_WHOLE :
      sip = (SeqIdPtr) location->data.ptrvalue;
      if (sip == NULL) return;
      wholebsp = BioseqFind (sip);
      if (wholebsp == NULL) return;
      if (is_real_id (sip, bsp->id)) {
        FlatLocPoint (head, sip, bsp->id, 0, NULL);
        if (bsp->length > 0) {
          ValNodeCopyStr (head, 0, "..");
          FlatLocPoint (head, NULL, bsp->id, bsp->length - 1, NULL);
        }
      }
      break;
    case SEQLOC_INT :
      sintp = (SeqIntPtr) location->data.ptrvalue;
      if (sintp == NULL) return;
      sip = sintp->id;
      if (sip == NULL) return;
      if (is_real_id (sip, bsp->id)) {
        minus_strand = (Boolean) (sintp->strand == Seq_strand_minus);
        if (minus_strand) {
          ValNodeCopyStr (head, 0, "complement(");
        }
        FlatLocPoint (head, sip, bsp->id, sintp->from, sintp->if_from);
        if (sintp->to > 0 &&
            (sintp->to != sintp->from ||
             sintp->if_from != NULL ||
             sintp->if_to != NULL)) {
          ValNodeCopyStr (head, 0, "..");
          FlatLocPoint (head, NULL, bsp->id, sintp->to, sintp->if_to);
        }
        if (minus_strand) {
          ValNodeCopyStr (head, 0, ")");
        }
      }
      break;
    case SEQLOC_PNT :
      spp = (SeqPntPtr) location->data.ptrvalue;
      if (spp == NULL) return;
      sip = spp->id;
      if (sip == NULL) return;
      if (is_real_id (sip, bsp->id)) {
        minus_strand = (Boolean) (spp->strand == Seq_strand_minus);
        if (minus_strand) {
          ValNodeCopyStr (head, 0, "complement(");
        }
        if (spp->fuzz != NULL) {
          FlatLocCaret (head, sip, bsp->id, spp->point, spp->fuzz);
        } else {
          FlatLocPoint (head, sip, bsp->id, spp->point, NULL);
        }
        if (minus_strand) {
          ValNodeCopyStr (head, 0, ")");
        }
      }
      break;
    case SEQLOC_BOND :
      sbp = (SeqBondPtr) location->data.ptrvalue;
      if (sbp == NULL) return;
      spp = sbp->a;
      if (spp == NULL) return;
      sip = spp->id;
      if (sip == NULL) return;
      ValNodeCopyStr (head, 0, "bond(");
      FlatLocPoint (head, sip, bsp->id, spp->point, spp->fuzz);
      spp = sbp->b;
      if (spp != NULL) {
        ValNodeCopyStr (head, 0, ",");
        FlatLocPoint (head, NULL, bsp->id, spp->point, spp->fuzz);
      }
      ValNodeCopyStr (head, 0, ")");
      break;
    default :
      /* unexpected internal complex type or unimplemented SEQLOC_FEAT */
      return;
  }
}

static Boolean FlatNullAhead (
  BioseqPtr bsp,
  ValNodePtr location
)

{
  SeqLocPtr  next;

  if (bsp == NULL || location == NULL) return FALSE;

  next = location->next;
  if (next == NULL) return TRUE;
  if (next->choice == SEQLOC_NULL) return TRUE;
  if (FlatVirtLoc (bsp, next)) return TRUE;

  return FALSE;
}

static void FlatPackedPoint (
  ValNodePtr PNTR head,
  PackSeqPntPtr pspp,
  BioseqPtr bsp
)

{
  Uint1  dex;

  if (head == NULL || pspp == NULL || bsp == NULL) return;

  for (dex = 0; dex < pspp->used; dex++) {
    FlatLocPoint (head, pspp->id, bsp->id, pspp->pnts [dex], pspp->fuzz);
  }
}

static void DoFlatLoc (
  ValNodePtr PNTR head,
  BioseqPtr bsp,
  SeqLocPtr location,
  Boolean ok_to_complement
);

static void GroupFlatLoc (
  ValNodePtr PNTR head,
  BioseqPtr bsp,
  SeqLocPtr location,
  CharPtr prefix,
  Boolean is_flat_order
)

{
  Boolean        found_non_virt = FALSE;
  SeqIdPtr       hold_next;
  Int2           parens = 1;
  PackSeqPntPtr  pspp;
  SeqLocPtr      slp;
  Boolean        special_mode = FALSE; /* join in order */

  if (head == NULL || bsp == NULL || location == NULL) return;

  /* prefix will have the first parenthesis */

  ValNodeCopyStr (head, 0, prefix);

  for (slp = (SeqLocPtr) location->data.ptrvalue; slp != NULL; slp = slp->next) {

    if (slp->choice == SEQLOC_NULL || FlatVirtLoc (bsp, slp)) {
      if (slp != location && slp->next != NULL) {
        if (special_mode) {
          special_mode = FALSE;
          ValNodeCopyStr (head, 0, ")");
          parens--;
        }
      }
      continue;
    }

    if (found_non_virt && slp->choice != SEQLOC_EMPTY && slp->choice != SEQLOC_NULL) {
      ValNodeCopyStr (head, 0, ",");
    }

    switch (slp->choice) {
      case SEQLOC_WHOLE :
      case SEQLOC_PNT :
      case SEQLOC_BOND :
      case SEQLOC_FEAT :
        found_non_virt = TRUE;
        if (FlatVirtLoc (bsp, slp)) {
          if (slp != location && slp->next != NULL) {
            if (special_mode) {
              special_mode = FALSE;
              ValNodeCopyStr (head, 0, ")");
              parens--;
            }
          }
        } else {
          FlatLocElement (head, bsp, slp);
        }
        break;
      case SEQLOC_INT :
        found_non_virt = TRUE;
        if (is_flat_order && (! FlatNullAhead (bsp, slp))) {
          special_mode = TRUE;
          ValNodeCopyStr (head, 0, "join(");
          parens++;
        }
        FlatLocElement (head, bsp, slp);
        break;
      case SEQLOC_PACKED_PNT :
        found_non_virt = TRUE;
        pspp = (PackSeqPntPtr) slp->data.ptrvalue;
        if (pspp != NULL) {
          FlatPackedPoint (head, pspp, bsp);
        }
        break;
      case SEQLOC_PACKED_INT :
      case SEQLOC_MIX :
      case SEQLOC_EQUIV :
        found_non_virt = TRUE;
        hold_next = slp->next;
        slp->next = NULL;
        DoFlatLoc (head, bsp, slp, FALSE);
        slp->next = hold_next;
        break;
      default :
        break;
    }

  }

  while (parens > 0) {
    ValNodeCopyStr (head, 0, ")");
    parens--;
  }
}

static void DoFlatLoc (
  ValNodePtr PNTR head,
  BioseqPtr bsp,
  SeqLocPtr location,
  Boolean ok_to_complement
)

{
  Boolean        found_null;
  SeqLocPtr      next_loc;
  PackSeqPntPtr  pspp;
  SeqLocPtr      slp;

  if (head == NULL || bsp == NULL || location == NULL) return;

  /* deal with complement of entire location */

  if (ok_to_complement && SeqLocStrand (location) == Seq_strand_minus) {
    slp = AsnIoMemCopy ((Pointer) location,
                        (AsnReadFunc) SeqLocAsnRead,
                        (AsnWriteFunc) SeqLocAsnWrite);
    if (slp != NULL) {
      SeqLocRevCmp (slp);
      ValNodeCopyStr (head, 0, "complement(");
      DoFlatLoc (head, bsp, slp, FALSE);
      ValNodeCopyStr (head, 0, ")");
    }
    SeqLocFree (slp);
    return;
  }

  /* handle each location component */

  for (slp = location; slp != NULL; slp = slp->next) {

    if (slp->choice == SEQLOC_NULL || FlatVirtLoc (bsp, slp)) continue;

    /* print comma between components */

    if (slp != location) {
      ValNodeCopyStr (head, 0, ",");
    }

    switch (slp->choice) {
      case SEQLOC_MIX :
      case SEQLOC_PACKED_INT :
        found_null = FALSE;
        for (next_loc = (SeqLocPtr) slp->data.ptrvalue; next_loc != NULL; next_loc = next_loc->next) {
          if (next_loc->choice == SEQLOC_NULL || FlatVirtLoc (bsp, next_loc)) {
            found_null = TRUE;
          }
        }
        if (found_null) {
          GroupFlatLoc (head, bsp, slp, "order(", TRUE);
        } else {
          GroupFlatLoc (head, bsp, slp, "join(", FALSE);
        }
        break;
      case SEQLOC_EQUIV :
        GroupFlatLoc (head, bsp, slp, "one-of(", FALSE);
        break;
      case SEQLOC_PACKED_PNT :
        pspp = (PackSeqPntPtr) slp->data.ptrvalue;
        if (pspp != NULL) {
          FlatPackedPoint (head, pspp, bsp);
        }
        break;
      default :
        FlatLocElement (head, bsp, slp);
        break;
    }

  }
}

static CharPtr FlatLoc (
  BioseqPtr bsp,
  SeqLocPtr location,
  StlType style
)

{
  ValNodePtr  head = NULL;
  SeqLocPtr   loc;
  Boolean     noLeft;
  Boolean     noRight;
  CharPtr     str;

  if (bsp == NULL || location == NULL) return NULL;

  if (! order_initialized) {
    order[SEQID_GENBANK ] = 1;
    order[SEQID_EMBL ] = 2;
    order[SEQID_DDBJ ] = 3;
    order[SEQID_LOCAL ] = 4;
    order[SEQID_GI ] = 5;
    order[SEQID_GIBBSQ ] = 6;
    order[SEQID_GIBBMT ] = 7;
    order[SEQID_PRF ] = 8;
    order[SEQID_PDB ] = 9;
    order[SEQID_PIR ] = 10;
    order[SEQID_SWISSPROT ] = 11;
    order[SEQID_PATENT ] = 12;
    order[SEQID_OTHER ] = 13;
    order[SEQID_GENERAL ] = 14;
    order[SEQID_GIIM ] = 15;
    order_initialized = TRUE;
  }

  if (style == MASTER_STYLE) {

    /* map location from parts to segmented bioseq */

    CheckSeqLocForPartial (location, &noLeft, &noRight);
    loc = SeqLocMerge (bsp, location, NULL, FALSE, TRUE, FALSE);
    if (loc == NULL) {
      return StringSave ("?");
    }
    FreeAllFuzz (loc);
    SetSeqLocPartial (loc, noLeft, noRight);

    DoFlatLoc (&head, bsp, loc, TRUE);

    SeqLocFree (loc);

  } else {

    DoFlatLoc (&head, bsp, location, TRUE);
  }

  str = MergeValNodeStrings (head);
  ValNodeFreeData (head);

  return str;
}

/* enumerated qualifier category definitions */

typedef enum {
  Qual_class_ignore = 0,
  Qual_class_string,
  Qual_class_boolean,
  Qual_class_int,
  Qual_class_evidence,
  Qual_class_valnode,
  Qual_class_quote,
  Qual_class_noquote,
  Qual_class_paren,
  Qual_class_region,
  Qual_class_replace,
  Qual_class_L_R_B,
  Qual_class_organelle,
  Qual_class_orgmod,
  Qual_class_subsource,
  Qual_class_code_break,
  Qual_class_anti_codon,
  Qual_class_codon,
  Qual_class_method,
  Qual_class_pubset,
  Qual_class_db_xref,
  Qual_class_seq_id,
  Qual_class_its,
  Qual_class_trna_codons,
  Qual_class_translation,
  Qual_class_protnames,
  Qual_class_illegal,
  Qual_class_note
}  QualType;

/* source 'feature' */

/* some qualifiers will require additional content verification not
   explicitly indicated by the class type */

typedef enum {
  SOURCE_acronym = 1,
  SOURCE_biotype,
  SOURCE_biovar,
  SOURCE_cell_line,
  SOURCE_cell_type,
  SOURCE_chemovar,
  SOURCE_chromosome,
  SOURCE_citation,
  SOURCE_clone,
  SOURCE_clone_lib,
  SOURCE_common,
  SOURCE_common_name,
  SOURCE_country,
  SOURCE_cultivar,
  SOURCE_db_xref,
  SOURCE_dev_stage,
  SOURCE_dosage,
  SOURCE_extrachrom,
  SOURCE_focus,
  SOURCE_frequency,
  SOURCE_genotype,
  SOURCE_germline,
  SOURCE_group,
  SOURCE_haplotype,
  SOURCE_ins_seq_name,
  SOURCE_isolate,
  SOURCE_lab_host,
  SOURCE_label,
  SOURCE_macronuclear,
  SOURCE_map,
  SOURCE_note,
  SOURCE_old_name,
  SOURCE_organism,
  SOURCE_organelle,
  SOURCE_orgmod_note,
  SOURCE_pathovar,
  SOURCE_plasmid_name,
  SOURCE_plastid_name,
  SOURCE_pop_variant,
  SOURCE_rearranged,
  SOURCE_seqfeat_note,
  SOURCE_sequenced_mol,
  SOURCE_serogroup,
  SOURCE_serotype,
  SOURCE_serovar,
  SOURCE_sex,
  SOURCE_spec_or_nat_host,
  SOURCE_specimen_voucher,
  SOURCE_strain,
  SOURCE_sub_clone,
  SOURCE_sub_group,
  SOURCE_sub_species,
  SOURCE_sub_strain,
  SOURCE_sub_type,
  SOURCE_subsource_note,
  SOURCE_tissue_lib,
  SOURCE_tissue_type,
  SOURCE_transposon_name,
  SOURCE_type,
  SOURCE_usedin,
  SOURCE_variety,
  ASN2GNBK_TOTAL_SOURCE
}  SourceType;

/* ordering arrays for qualifiers and note components */

static Uint1 relmode_source_qual_order [] = {
  SOURCE_organism,
  SOURCE_focus,

  SOURCE_organelle,

  SOURCE_strain,
  SOURCE_sub_strain,
  SOURCE_variety,
  SOURCE_serotype,
  SOURCE_cultivar,
  SOURCE_isolate,
  SOURCE_spec_or_nat_host,
  SOURCE_sub_species,
  SOURCE_specimen_voucher,

  SOURCE_db_xref,

  SOURCE_chromosome,
  SOURCE_map,
  SOURCE_clone,
  SOURCE_sub_clone,
  SOURCE_haplotype,
  SOURCE_sex,
  SOURCE_cell_line,
  SOURCE_cell_type,
  SOURCE_tissue_type,
  SOURCE_clone_lib,
  SOURCE_dev_stage,
  SOURCE_frequency,

  SOURCE_germline,
  SOURCE_rearranged,

  SOURCE_lab_host,
  SOURCE_pop_variant,
  SOURCE_tissue_lib,

  SOURCE_plasmid_name,
  SOURCE_transposon_name,
  SOURCE_ins_seq_name,

  SOURCE_country,

  SOURCE_sequenced_mol, SOURCE_label, SOURCE_usedin,
  SOURCE_citation, SOURCE_note, 0
};

static Uint1 relmode_source_note_order [] = {
  SOURCE_seqfeat_note, SOURCE_orgmod_note, SOURCE_subsource_note,

  SOURCE_type,
  SOURCE_sub_type,
  SOURCE_serogroup,
  SOURCE_serovar,
  SOURCE_pathovar,
  SOURCE_chemovar,
  SOURCE_biovar,
  SOURCE_biotype,
  SOURCE_group,
  SOURCE_sub_group,
  SOURCE_common,
  SOURCE_common_name,
  SOURCE_acronym,
  SOURCE_dosage,
  SOURCE_genotype,
  SOURCE_plastid_name,
  /* SOURCE_old_name, */
  0
};

static Uint1 seqmode_source_qual_order [] = {
  SOURCE_organism,
  SOURCE_focus,

  SOURCE_organelle,

  SOURCE_strain,
  SOURCE_sub_strain,
  SOURCE_variety,
  SOURCE_serotype,
  SOURCE_cultivar,
  SOURCE_isolate,
  SOURCE_spec_or_nat_host,
  SOURCE_sub_species,
  SOURCE_specimen_voucher,

  SOURCE_db_xref,

  SOURCE_chromosome,
  SOURCE_map,
  SOURCE_clone,
  SOURCE_sub_clone,
  SOURCE_haplotype,
  SOURCE_sex,
  SOURCE_cell_line,
  SOURCE_cell_type,
  SOURCE_tissue_type,
  SOURCE_clone_lib,
  SOURCE_dev_stage,
  SOURCE_frequency,

  SOURCE_germline,
  SOURCE_rearranged,

  SOURCE_lab_host,
  SOURCE_pop_variant,
  SOURCE_tissue_lib,

  SOURCE_plasmid_name,
  SOURCE_transposon_name,
  SOURCE_ins_seq_name,

  SOURCE_country,

  SOURCE_sequenced_mol, SOURCE_label, SOURCE_usedin,
  SOURCE_citation, SOURCE_note, 0
};

static Uint1 seqmode_source_note_order [] = {
  SOURCE_seqfeat_note, SOURCE_orgmod_note, SOURCE_subsource_note,

  SOURCE_type,
  SOURCE_sub_type,
  SOURCE_serogroup,
  SOURCE_serovar,
  SOURCE_pathovar,
  SOURCE_chemovar,
  SOURCE_biovar,
  SOURCE_biotype,
  SOURCE_group,
  SOURCE_sub_group,
  SOURCE_common,
  SOURCE_common_name,
  SOURCE_acronym,
  SOURCE_dosage,
  SOURCE_genotype,
  SOURCE_plastid_name,
  /* SOURCE_old_name, */
  0
};

typedef struct sourcequal {
  CharPtr     name;
  Uint1       qualclass;
} SourceQual, PNTR SourceQualPtr;

static SourceQual asn2gnbk_source_quals [ASN2GNBK_TOTAL_SOURCE] = {
  { "",                 Qual_class_ignore    },
  { "acronym",          Qual_class_orgmod    },
  { "biotype",          Qual_class_orgmod    },
  { "biovar",           Qual_class_orgmod    },
  { "cell_line",        Qual_class_subsource },
  { "cell_type",        Qual_class_subsource },
  { "chemovar",         Qual_class_orgmod    },
  { "chromosome",       Qual_class_subsource },
  { "citation",         Qual_class_pubset    },
  { "clone",            Qual_class_subsource },
  { "clone_lib",        Qual_class_subsource },
  { "common",           Qual_class_orgmod    },
  { "common",           Qual_class_string    },
  { "country",          Qual_class_subsource },
  { "cultivar",         Qual_class_orgmod    },
  { "db_xref",          Qual_class_db_xref   },
  { "dev_stage",        Qual_class_subsource },
  { "dosage",           Qual_class_orgmod    },
  { "extrachromosomal", Qual_class_boolean   },
  { "focus",            Qual_class_boolean   },
  { "frequency",        Qual_class_subsource },
  { "genotype",         Qual_class_subsource },
  { "germline",         Qual_class_subsource },
  { "group",            Qual_class_orgmod    },
  { "haplotype",        Qual_class_subsource },
  { "insertion_seq",    Qual_class_subsource },
  { "isolate",          Qual_class_orgmod    },
  { "lab_host",         Qual_class_subsource },
  { "label",            Qual_class_noquote   },
  { "macronuclear",     Qual_class_boolean   },
  { "map",              Qual_class_subsource },
  { "note",             Qual_class_note      },
  { "old_name",         Qual_class_orgmod    },
  { "organism",         Qual_class_string    },
  { "organelle",        Qual_class_organelle },
  { "orgmod_note",      Qual_class_orgmod    },
  { "pathovar",         Qual_class_orgmod    },
  { "plasmid",          Qual_class_subsource },
  { "plastid",          Qual_class_subsource },
  { "pop_variant",      Qual_class_subsource },
  { "rearranged",       Qual_class_subsource },
  { "seqfeat_note",     Qual_class_string    },
  { "sequenced_mol",    Qual_class_quote     },
  { "serogroup",        Qual_class_orgmod    },
  { "serotype",         Qual_class_orgmod    },
  { "serovar",          Qual_class_orgmod    },
  { "sex",              Qual_class_subsource },
  { "specific_host",    Qual_class_orgmod    },
  { "specimen_voucher", Qual_class_orgmod    },
  { "strain",           Qual_class_orgmod    },
  { "sub_clone",        Qual_class_subsource },
  { "subgroup",         Qual_class_orgmod    },
  { "sub_species",      Qual_class_orgmod    },
  { "sub_strain",       Qual_class_orgmod    },
  { "subtype",          Qual_class_orgmod    },
  { "subsource_note",   Qual_class_subsource },
  { "tissue_lib",       Qual_class_subsource },
  { "tissue_type",      Qual_class_subsource },
  { "transposon",       Qual_class_subsource },
  { "type",             Qual_class_orgmod    },
  { "usedin",           Qual_class_quote     },
  { "variety",          Qual_class_orgmod    },
};

static Int2 subSourceToSourceIdx [25] = {
  0,
  SOURCE_chromosome,
  SOURCE_map,
  SOURCE_clone,
  SOURCE_sub_clone,
  SOURCE_haplotype,
  SOURCE_genotype,
  SOURCE_sex,
  SOURCE_cell_line,
  SOURCE_cell_type,
  SOURCE_tissue_type,
  SOURCE_clone_lib,
  SOURCE_dev_stage,
  SOURCE_frequency,
  SOURCE_germline,
  SOURCE_rearranged,
  SOURCE_lab_host,
  SOURCE_pop_variant,
  SOURCE_tissue_lib,
  SOURCE_plasmid_name,
  SOURCE_transposon_name,
  SOURCE_ins_seq_name,
  SOURCE_plastid_name,
  SOURCE_country,
  SOURCE_subsource_note
};

static void SubSourceToQualArray (
  SubSourcePtr ssp,
  QualValPtr qvp
)

{
  Int2   idx;
  Uint1  subtype;

  if (ssp == NULL || qvp == NULL) return;

  while (ssp != NULL) {
    subtype = ssp->subtype;
    if (subtype == 255) {
      subtype = 24;
    }
    if (subtype < 25) {
      idx = subSourceToSourceIdx [subtype];
      if (idx > 0 && idx < ASN2GNBK_TOTAL_SOURCE) {
        if (qvp [idx].ssp == NULL) {
          qvp [idx].ssp = ssp;
        }
      }
    }
    ssp = ssp->next;
  }
}

static Int2 orgModToSourceIdx [26] = {
  0,
  0,
  SOURCE_strain,
  SOURCE_sub_strain,
  SOURCE_type,
  SOURCE_sub_type,
  SOURCE_variety,
  SOURCE_serotype,
  SOURCE_serogroup,
  SOURCE_serovar,
  SOURCE_cultivar,
  SOURCE_pathovar,
  SOURCE_chemovar,
  SOURCE_biovar,
  SOURCE_biotype,
  SOURCE_group,
  SOURCE_sub_group,
  SOURCE_isolate,
  SOURCE_common,
  SOURCE_acronym,
  SOURCE_dosage,
  SOURCE_spec_or_nat_host,
  SOURCE_sub_species,
  SOURCE_specimen_voucher,
  SOURCE_old_name,
  SOURCE_orgmod_note
};

static void OrgModToQualArray (
  OrgModPtr omp,
  QualValPtr qvp
)

{
  Int2   idx;
  Uint1  subtype;

  if (omp == NULL || qvp == NULL) return;

  while (omp != NULL) {
    subtype = omp->subtype;
    if (subtype == 254) {
      subtype = 24;
    } else if (subtype == 255) {
      subtype = 25;
    }
    if (subtype < 26) {
      idx = orgModToSourceIdx [subtype];
      if (idx > 0 && idx < ASN2GNBK_TOTAL_SOURCE) {
        if (qvp [idx].omp == NULL) {
          qvp [idx].omp = omp;
        }
      }
    }
    omp = omp->next;
  }
}

/*
static CharPtr organelleQual [] = {
  NULL,
  NULL,
  "/chloroplast",
  "/chromoplast",
  "/kinetoplast",
  "/mitochondrion",
  "/plastid",
  "/macronuclear",
  "/extrachrom",
  "/plasmid", 
  "/transposon",
  "/insertion_seq",
  "/cyanelle",
  "/proviral",
  "/virion",
  NULL,
  NULL,
  NULL,
  NULL
};
*/

static CharPtr organelleQual [] = {
  NULL,
  NULL,
  "/organelle=\"plastid:chloroplast\"",
  "/organelle=\"plastid:chromoplast\"",
  "/organelle=\"mitochondrion:kinetoplast\"",
  "/organelle=\"mitochondrion\"",
  "/organelle=\"plastid\"",
  "/macronuclear",
  "/extrachrom",
  "/plasmid", 
  "/transposon",
  "/insertion_seq",
  "/organelle=\"plastid:cyanelle\"",
  "/proviral",
  "/virion",
  "/organelle=\"nucleomorph\"",
  "/organelle=\"plastid:apicoplast\"",
  "/organelle=\"plastid:leucoplast\"",
  "/organelle=\"plastid:proplastid\""
};

static Boolean StringIsJustQuotes (
  CharPtr str
)

{
  Nlm_Uchar  ch;	/* to use 8bit characters in multibyte languages */

  if (str != NULL) {
    ch = *str;
    while (ch != '\0') {
      if (ch > ' ' && ch != '"' && ch != '\'') {
        return FALSE;
      }
      str++;
      ch = *str;
    }
  }
  return TRUE;
}

static CharPtr FormatSourceFeatBlock (
  Asn2gbFormatPtr afp,
  BaseBlockPtr bbp
)

{
  Asn2gbJobPtr       ajp;
  Asn2gbSectionPtr   asp;
  BioSourcePtr       biop = NULL;
  BioseqPtr          bsp;
  Char               buf [80];
  CharPtr            common = NULL;
  DbtagPtr           dbt;
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext;
  ValNodePtr         head;
  Int2               i;
  Uint1              idx;
  Int2               j;
  Uint1              jdx;
  Uint1              lastomptype;
  Uint1              lastssptype;
  SeqLocPtr          location = NULL;
  CharPtr            notestr;
  Uint1Ptr           notetbl = NULL;
  ObjectIdPtr        oip;
  OrgModPtr          omp;
  OrgNamePtr         onp = NULL;
  OrgRefPtr          orp = NULL;
  CharPtr            prefix;
  Uint1Ptr           qualtbl = NULL;
  QualValPtr         qvp;
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp = NULL;
  SeqInt             sin;
  SubSourcePtr       ssp;
  CharPtr            str;
  CharPtr            taxname = NULL;
  ValNode            vn;
  ValNodePtr         vnp;
  CharPtr            wwwbuf;

  if (afp == NULL || bbp == NULL) return NULL;
  ajp = afp->ajp;
  if (ajp == NULL) return NULL;
  asp = afp->asp;
  if (asp == NULL) return NULL;
  bsp = asp->bsp;
  if (bsp == NULL) return NULL;
  qvp = afp->qvp;
  if (qvp == NULL) return NULL;

  if (! StringHasNoText (bbp->string)) return StringSave (bbp->string);

  /* could be descriptor or feature */

  if (bbp->itemtype == OBJ_SEQDESC) {
    sdp = SeqMgrGetDesiredDescriptor (bbp->entityID, NULL, bbp->itemID, 0, NULL, &dcontext);
    if (sdp != NULL && dcontext.seqdesctype == Seq_descr_source) {
      biop = (BioSourcePtr) sdp->data.ptrvalue;
    }
  } else if (bbp->itemtype == OBJ_SEQFEAT) {
    sfp = SeqMgrGetDesiredFeature (bbp->entityID, NULL, bbp->itemID, 0, NULL, &fcontext);
    if (sfp != NULL && fcontext.seqfeattype == SEQFEAT_BIOSRC) {
      biop = (BioSourcePtr) sfp->data.value.ptrvalue;
    }
  }

  if (biop == NULL) return NULL;

  orp = biop->org;
  if (orp != NULL) {
    taxname = orp->taxname;
    /* common = orp->common; */
  }
  if (StringHasNoText (taxname)) {
    taxname = "Unknown.";
  }

  gb_StartPrint (afp->format, TRUE, 5, 21, NULL, 0, 5, 21, "FT", FALSE);
  ff_AddString ("source");
  TabToColumn (22);

  if (sfp == NULL) {
    sin.from = 0;
    sin.to = bsp->length - 1;
    sin.strand = Seq_strand_plus;
    sin.id = SeqIdFindBest (bsp->id, 0);
    sin.if_from = NULL;
    sin.if_to = NULL;

    vn.choice = SEQLOC_INT;
    vn.data.ptrvalue = &sin;
    vn.next = NULL;

    location = &vn;
  } else {
    location = sfp->location;
  }
  str = FlatLoc (bsp, location, afp->style);
  if (get_www ()) {
    wwwbuf = www_featloc (str);
    ff_AddString (wwwbuf);
    MemFree (wwwbuf);
  } else {
    ff_AddString (str);
  }
  MemFree (str);

  /* populate qualifier table from biosource fields */

  qvp [SOURCE_organism].str = taxname;
  qvp [SOURCE_common_name].str = common;

  if (biop->is_focus) {
    qvp [SOURCE_focus].bool = TRUE;
  }

  SubSourceToQualArray (biop->subtype, qvp);

  if (orp != NULL) {
    onp = orp->orgname;
    if (onp != NULL) {
      OrgModToQualArray (onp->mod, qvp);
    }

    qvp [SOURCE_db_xref].vnp = orp->db;
  }

  /* organelle currently prints /mitochondrion, /virion, etc. */

  qvp [SOURCE_organelle].num = biop->genome;

  /* some qualifiers are flags in genome and names in subsource, print once with name */

  if (qvp [SOURCE_ins_seq_name].ssp != NULL &&
      qvp [SOURCE_organelle].num == GENOME_insertion_seq) {
    qvp [SOURCE_organelle].num = 0;
  }
  if (qvp [SOURCE_plasmid_name].ssp != NULL &&
      qvp [SOURCE_organelle].num == GENOME_plasmid) {
    qvp [SOURCE_organelle].num = 0;
  }
  if (qvp [SOURCE_plastid_name].ssp != NULL &&
      qvp [SOURCE_organelle].num == GENOME_plastid) {
    qvp [SOURCE_organelle].num = 0;
  }
  if (qvp [SOURCE_transposon_name].ssp != NULL &&
      qvp [SOURCE_organelle].num == GENOME_transposon) {
    qvp [SOURCE_organelle].num = 0;
  }

  if (sfp != NULL) {
    qvp [SOURCE_seqfeat_note].str = sfp->comment;
  }

  /* now print qualifiers from table */

  if (afp->mode == RELEASE_MODE) {
    qualtbl = relmode_source_qual_order;
    notetbl = relmode_source_note_order;
  } else {
    qualtbl = seqmode_source_qual_order;
    notetbl = seqmode_source_note_order;
  }

  for (i = 0, idx = qualtbl [i]; idx != 0; i++, idx = qualtbl [i]) {

    lastomptype = 0;
    lastssptype = 0;
    switch (asn2gnbk_source_quals [idx].qualclass) {

      case Qual_class_ignore :
        break;

      case Qual_class_string :
        if (! StringHasNoText (qvp [idx].str)) {
          NewContLine ();
          sprintf (buf, "/%s=\"", asn2gnbk_source_quals [idx].name);
          gb_AddString (buf, qvp [idx].str, "\"", FALSE, TRUE, FALSE);
        }
        break;

      case Qual_class_boolean :
        if (qvp [idx].bool) {
          NewContLine ();
          sprintf (buf, "/%s", asn2gnbk_source_quals [idx].name);
          ff_AddString (buf);
        }
        break;

      case Qual_class_organelle :
        j = (Int2) qvp [idx].num;
        if (organelleQual [j] != NULL) {
          NewContLine ();
          sprintf (buf, "%s", organelleQual [j]);
          ff_AddString (buf);
        }
        break;

      case Qual_class_orgmod :
        omp = qvp [idx].omp;
        if (lastomptype == 0 && omp != NULL) {
          lastomptype = omp->subtype;
        }
        while (omp != NULL && omp->subtype == lastomptype) {
          if (StringIsJustQuotes (omp->subname)) {
            NewContLine ();
            sprintf (buf, "/%s=\"", asn2gnbk_source_quals [idx].name);
            gb_AddString (buf, NULL, "\"", FALSE, TRUE, FALSE);
          } else if (! StringHasNoText (omp->subname)) {
            NewContLine ();
            sprintf (buf, "/%s=\"", asn2gnbk_source_quals [idx].name);
            gb_AddString (buf, omp->subname, "\"", FALSE, TRUE, FALSE);
          }
          omp = omp->next;
        }
        break;

      case Qual_class_subsource :
        ssp = qvp [idx].ssp;
        if (lastssptype == 0 && ssp != NULL) {
          lastssptype = ssp->subtype;
        }
        while (ssp != NULL && ssp->subtype == lastssptype) {
          if (ssp->subtype == SUBSRC_germline || ssp->subtype == SUBSRC_rearranged) {
            NewContLine ();
            sprintf (buf, "/%s", asn2gnbk_source_quals [idx].name);
            gb_AddString (buf, NULL, NULL, FALSE, TRUE, FALSE);
          } else if (StringIsJustQuotes (ssp->name)) {
            NewContLine ();
            sprintf (buf, "/%s=\"", asn2gnbk_source_quals [idx].name);
            gb_AddString (buf, NULL, "\"", FALSE, TRUE, FALSE);
          } else if (! StringHasNoText (ssp->name)) {
            NewContLine ();
            sprintf (buf, "/%s=\"", asn2gnbk_source_quals [idx].name);
            gb_AddString (buf, ssp->name, "\"", FALSE, TRUE, FALSE);
          } else if (ssp->subtype == SUBSRC_germline ||
                     ssp->subtype == SUBSRC_rearranged) {
            NewContLine ();
            gb_AddString ("/", asn2gnbk_source_quals [idx].name, NULL, FALSE, TRUE, FALSE);
          }
          ssp = ssp->next;
        }
        break;

      case Qual_class_pubset :
        break;

      case Qual_class_quote :
        break;

      case Qual_class_noquote :
        break;

      case Qual_class_db_xref :
        for (vnp = qvp [idx].vnp; vnp != NULL; vnp = vnp->next) {
          buf [0] = '\0';
          dbt = (DbtagPtr) vnp->data.ptrvalue;
          if (dbt != NULL && (! StringHasNoText (dbt->db))) {
            oip = dbt->tag;
            if (oip != NULL) {
              /* if release mode, drop unknown dbtag */
              if (! StringHasNoText (oip->str)) {
                if (StringLen (dbt->db) + StringLen (oip->str) < 80) {
                  sprintf (buf, "%s:%s", dbt->db, oip->str);
                }
              } else {
                sprintf (buf, "%s:%ld", dbt->db, (long) oip->id);
              }
            }
          }
          if (! StringHasNoText (buf)) {
            NewContLine ();
            gb_AddString ("/db_xref=\"", buf, "\"", FALSE, TRUE, FALSE);
          }
        }
        break;

      case Qual_class_illegal :
        break;

      case Qual_class_note :
        head = NULL;
        notestr = NULL;
        prefix = "";
        for (j = 0, jdx = notetbl [j]; jdx != 0; j++, jdx = notetbl [j]) {

          lastomptype = 0;
          lastssptype = 0;
          switch (asn2gnbk_source_quals [jdx].qualclass) {

            case Qual_class_string :
              if (! StringHasNoText (qvp [jdx].str)) {
                AddValNodeString (&head, prefix, qvp [jdx].str, NULL);
                /* prefix = "; "; */
                prefix = "\n";
              }
              break;
            case Qual_class_orgmod :
              omp = qvp [jdx].omp;
              if (lastomptype == 0 && omp != NULL) {
                lastomptype = omp->subtype;
              }
              while (omp != NULL && omp->subtype == lastomptype) {
                if (! StringHasNoText (omp->subname)) {
                  if (jdx == SOURCE_orgmod_note) {
                    sprintf (buf, "%s", prefix);
                  } else {
                    sprintf (buf, "%s%s: ", prefix, asn2gnbk_source_quals [jdx].name);
                  }
                  AddValNodeString (&head, buf, omp->subname, NULL);
                  if (jdx == SOURCE_orgmod_note) {
                    prefix = "\n";
                  } else {
                    prefix = "; ";
                  }
                }
                omp = omp->next;
              }
              break;
            case Qual_class_subsource :
              ssp = qvp [jdx].ssp;
              if (lastssptype == 0 && ssp != NULL) {
                lastssptype = ssp->subtype;
              }
              while (ssp != NULL && ssp->subtype == lastssptype) {
                if (! StringHasNoText (ssp->name)) {
                  if (jdx == SOURCE_subsource_note) {
                    sprintf (buf, "%s", prefix);
                  } else {
                    sprintf (buf, "%s%s: ", prefix, asn2gnbk_source_quals [jdx].name);
                  }
                  AddValNodeString (&head, buf, ssp->name, NULL);
                  if (jdx == SOURCE_subsource_note) {
                    prefix = "\n";
                  } else {
                    prefix = "; ";
                  }
                }
                ssp = ssp->next;
              }
              break;
            default :
              break;
          }
        }
        if (head != NULL) {
          NewContLine ();
          notestr = MergeValNodeStrings (head);
          gb_AddString ("/note=\"", notestr, "\"", FALSE, TRUE, FALSE);
          MemFree (notestr);
          ValNodeFreeData (head);
        }
        break;
      default :
        break;
    }
  }

  /* and then deal with the various note types separately (not in order table) */

  return gb_MergeString (TRUE);
}

static CharPtr featurekeys [] = {
  "???" ,
  "gene" ,
  "Org" ,
  "CDS" ,
  "Protein" ,
  "precursor_RNA" ,
  "mRNA" ,
  "tRNA" ,
  "rRNA" ,
  "snRNA" ,
  "scRNA" ,
  "misc_RNA" ,
  "Cit" ,
  "Xref" ,
  "Imp" ,
  "allele" ,
  "attenuator" ,
  "C_region" ,
  "CAAT_signal" ,
  "CDS" ,
  "conflict" ,
  "D-loop" ,
  "D_segment" ,
  "enhancer" ,
  "exon" ,
  "GC_signal" ,
  "iDNA" ,
  "intron" ,
  "J_segment" ,
  "LTR" ,
  "mat_peptide" ,
  "misc_binding" ,
  "misc_difference" ,
  "misc_feature" ,
  "misc_recomb" ,
  "misc_RNA" ,
  "misc_signal" ,
  "misc_structure" ,
  "modified_base" ,
  "mutation" ,
  "N_region" ,
  "old_sequence" ,
  "polyA_signal" ,
  "polyA_site" ,
  "precursor_RNA" ,
  "prim_transcript" ,
  "primer_bind" ,
  "promoter" ,
  "protein_bind" ,
  "RBS" ,
  "repeat_region" ,
  "repeat_unit" ,
  "rep_origin" ,
  "S_region" ,
  "satellite" ,
  "sig_peptide" ,
  "source" ,
  "stem_loop" ,
  "STS" ,
  "TATA_signal" ,
  "terminator" ,
  "transit_peptide" ,
  "unsure" ,
  "V_region" ,
  "V_segment" ,
  "variation" ,
  "virion" ,
  "3'clip" ,
  "3'UTR" ,
  "5'clip" ,
  "5'UTR" ,
  "-10_signal" ,
  "-35_signal" ,
  "Site-ref" ,
  "misc_feature" ,
  "misc_feature" ,
  "Bond" ,
  "misc_feature" ,
  "Rsite" ,
  "User" ,
  "TxInit" ,
  "Num" ,
  "SecStr" ,
  "NonStdRes" ,
  "Het" ,
  "Src" ,
  "preprotein" ,
  "mat_peptide" ,
  "sig_peptide" ,
  "transit_peptide"
};

typedef enum {
  FEATUR_allele = 1,
  FEATUR_anticodon,
  FEATUR_bound_moiety,
  FEATUR_cds_product,
  FEATUR_citation,
  FEATUR_clone,
  FEATUR_codon,
  FEATUR_codon_start,
  FEATUR_cons_splice,
  FEATUR_db_xref,
  FEATUR_direction,
  FEATUR_EC_number,
  FEATUR_evidence,
  FEATUR_exception,
  FEATUR_figure,
  FEATUR_frequency,
  FEATUR_function,
  FEATUR_gene,
  FEATUR_gene_desc,
  FEATUR_gene_map,
  FEATUR_gene_syn,
  FEATUR_gene_note,
  FEATUR_illegal_qual,
  FEATUR_label,
  FEATUR_map,
  FEATUR_maploc,
  FEATUR_mod_base,
  FEATUR_note,
  FEATUR_number,
  FEATUR_organism,
  FEATUR_partial,
  FEATUR_PCR_conditions,
  FEATUR_phenotype,
  FEATUR_product,
  FEATUR_product_quals,
  FEATUR_prot_activity,
  FEATUR_prot_comment,
  FEATUR_prot_note,
  FEATUR_prot_method,
  FEATUR_prot_conflict,
  FEATUR_prot_desc,
  FEATUR_prot_missing,
  FEATUR_prot_names,
  FEATUR_protein_id,
  FEATUR_pseudo,
  FEATUR_region,
  FEATUR_replace,
  FEATUR_rpt_family,
  FEATUR_rpt_type,
  FEATUR_rpt_unit,
  FEATUR_rrna_its,
  FEATUR_seqfeat_note,
  FEATUR_site,
  FEATUR_standard_name,
  FEATUR_transl_except,
  FEATUR_transl_table,
  FEATUR_translation,
  FEATUR_trna_aa,
  FEATUR_trna_codons,
  FEATUR_usedin,
  ASN2GNBK_TOTAL_FEATUR
}  FeaturType;

/* ordering arrays for qualifiers and note components */

static Uint1 relmode_feat_qual_order [] = {
  FEATUR_partial,
  FEATUR_gene,

  FEATUR_EC_number,
  FEATUR_prot_activity,

  FEATUR_note, FEATUR_citation,

  FEATUR_pseudo,

  FEATUR_number,

  FEATUR_product,

  FEATUR_allele,
  FEATUR_anticodon,
  FEATUR_bound_moiety,
  FEATUR_clone,
  FEATUR_codon,
  FEATUR_cons_splice,
  FEATUR_direction,
  FEATUR_evidence,
  FEATUR_exception,
  FEATUR_frequency,
  FEATUR_function,
  FEATUR_gene_map,
  FEATUR_map,
  FEATUR_mod_base,
  FEATUR_organism,
  FEATUR_PCR_conditions,
  FEATUR_phenotype,
  FEATUR_rpt_family,
  FEATUR_rpt_type,
  FEATUR_rpt_unit,
  FEATUR_standard_name,
  FEATUR_usedin,

  /* FEATUR_illegal_qual, */

  FEATUR_replace,

  FEATUR_codon_start,
  FEATUR_transl_except,
  FEATUR_transl_table,
  FEATUR_label,
  FEATUR_cds_product,
  FEATUR_protein_id, FEATUR_db_xref,
  FEATUR_translation,
  0
};

static Uint1 relmode_feat_note_order [] = {
  FEATUR_gene_syn,
  FEATUR_gene_desc,
  FEATUR_trna_codons,
  FEATUR_prot_desc,
  FEATUR_prot_note,
  FEATUR_prot_comment,
  FEATUR_prot_method,
  FEATUR_figure,
  FEATUR_maploc,
  FEATUR_prot_conflict,
  FEATUR_prot_missing,
  FEATUR_seqfeat_note,
  FEATUR_region,
  FEATUR_prot_names,
  FEATUR_site,
  FEATUR_rrna_its,
  0
};

/*
pseudo after note - gi|6598562|gb|AC006419.3|AC006419
*/

static Uint1 seqmode_feat_qual_order [] = {
  FEATUR_partial,
  FEATUR_gene,

  FEATUR_EC_number,
  FEATUR_prot_activity,

  FEATUR_note, FEATUR_citation,

  FEATUR_pseudo,

  FEATUR_number,

  FEATUR_product,

  FEATUR_allele,
  FEATUR_anticodon,
  FEATUR_bound_moiety,
  FEATUR_clone,
  FEATUR_codon,
  FEATUR_cons_splice,
  FEATUR_direction,
  FEATUR_evidence,
  FEATUR_exception,
  FEATUR_frequency,
  FEATUR_function,
  FEATUR_gene_map,
  FEATUR_map,
  FEATUR_mod_base,
  FEATUR_organism,
  FEATUR_PCR_conditions,
  FEATUR_phenotype,
  FEATUR_rpt_family,
  FEATUR_rpt_type,
  FEATUR_rpt_unit,
  FEATUR_standard_name,
  FEATUR_usedin,

  FEATUR_illegal_qual,

  FEATUR_replace,

  FEATUR_codon_start,
  FEATUR_transl_except,
  FEATUR_transl_table,
  FEATUR_label,
  FEATUR_cds_product,
  FEATUR_protein_id, FEATUR_db_xref,
  FEATUR_translation,
  0
};

/*
prot_names after seqfeat_note - gi|4210642|emb|AJ005084.1|HBVAJ5084
prot_conflict after prot_desc - gi|61183|emb|V01135.1|PIVM02
figure after prot_desc - gi|400553|gb|S64006.1|
seqfeat_note after prot_desc - gi|431713|gb|L20354.1|STVPATPOLB
prot_names after figure - gi|234022|gb|S56149.1|S56149
seqfeat_note after prot_conflict after figure - gi|234046|gb|S51392.1|S51392
prot_method after prot_comment (descriptor) after prot_note after prot_desc
region after seqfeat_note - gi|6554164|gb|AF043644.3|AF043644
prot_desc after prot_names - gi|6581069|gb|AF202541.1|AF202541 - cannot do !!!
*/

static Uint1 seqmode_feat_note_order [] = {
  FEATUR_gene_syn,
  FEATUR_gene_desc,
  FEATUR_trna_codons,
  FEATUR_prot_desc,
  FEATUR_prot_note,
  FEATUR_prot_comment,
  FEATUR_prot_method,
  FEATUR_figure,
  FEATUR_maploc,
  FEATUR_prot_conflict,
  FEATUR_prot_missing,
  FEATUR_seqfeat_note,
  FEATUR_region,
  FEATUR_prot_names,
  FEATUR_site,
  FEATUR_rrna_its,
  0
};

typedef struct featurqual {
  CharPtr     name;
  Uint1       qualclass;
} FeaturQual, PNTR FeaturQualPtr;

static FeaturQual asn2gnbk_featur_quals [ASN2GNBK_TOTAL_FEATUR] = {
  { "",               Qual_class_ignore       },
  { "allele",         Qual_class_string       },
  { "anticodon",      Qual_class_anti_codon   },
  { "bound_moiety",   Qual_class_quote        },
  { "product",        Qual_class_string       },
  { "citation",       Qual_class_pubset       },
  { "clone",          Qual_class_quote        },
  { "codon",          Qual_class_codon        },
  { "codon_start",    Qual_class_int          },
  { "cons_splice",    Qual_class_noquote      },
  { "db_xref",        Qual_class_db_xref      },
  { "direction",      Qual_class_L_R_B        },
  { "EC_number",      Qual_class_valnode      },
  { "evidence",       Qual_class_evidence     },
  { "exception",      Qual_class_string       },
  { "figure",         Qual_class_string       },
  { "frequency",      Qual_class_quote        },
  { "function",       Qual_class_quote        },
  { "gene",           Qual_class_string       },
  { "gene_desc",      Qual_class_string       },
  { "map",            Qual_class_string       },
  { "gene_syn",       Qual_class_valnode      },
  { "gene_note",      Qual_class_string       },
  { "illegal",        Qual_class_illegal      },
  { "label",          Qual_class_noquote      },
  { "map",            Qual_class_quote        },
  { "maploc",         Qual_class_string       },
  { "mod_base",       Qual_class_noquote      },
  { "note",           Qual_class_note         },
  { "number",         Qual_class_noquote      },
  { "organism",       Qual_class_quote        },
  { "partial",        Qual_class_boolean      },
  { "PCR_conditions", Qual_class_quote        },
  { "phenotype",      Qual_class_quote        },
  { "product",        Qual_class_string       },
  { "product",        Qual_class_quote        },
  { "function",       Qual_class_valnode      },
  { "prot_comment",   Qual_class_string       },
  { "prot_note",      Qual_class_string       },
  { "prot_method",    Qual_class_method       },
  { "prot_conflict",  Qual_class_string       },
  { "prot_desc",      Qual_class_string       },
  { "prot_missing",   Qual_class_string       },
  { "prot_names",     Qual_class_protnames    },
  { "protein_id",     Qual_class_seq_id       },
  { "pseudo",         Qual_class_boolean      },  
  { "region",         Qual_class_region       },
  { "replace",        Qual_class_replace      },
  { "rpt_family",     Qual_class_quote        },
  { "rpt_type",       Qual_class_paren        },
  { "rpt_unit",       Qual_class_paren        },
  { "rrna_its",       Qual_class_its          },
  { "seqfeat_note",   Qual_class_string       },
  { "site",           Qual_class_string       },
  { "standard_name",  Qual_class_quote        },
  { "transl_except",  Qual_class_code_break   },
  { "transl_table",   Qual_class_int          },
  { "translation",    Qual_class_translation  },
  { "trna_aa",        Qual_class_ignore       },
  { "trna_codons",    Qual_class_trna_codons  },
  { "usedin",         Qual_class_paren        }
};

typedef struct qualfeatur {
  CharPtr     name;
  Uint1       featurclass;
} QualFeatur, PNTR QualFeaturPtr;

#define NUM_GB_QUALS 20

static QualFeatur qualToFeature [NUM_GB_QUALS] = {
  { "bound_moiety",   FEATUR_bound_moiety   },
  { "clone",          FEATUR_clone          },
  { "cons_splice",    FEATUR_cons_splice    },
  { "direction",      FEATUR_direction      },
  { "frequency",      FEATUR_frequency      },
  { "function",       FEATUR_function       },
  { "label",          FEATUR_label          },
  { "map",            FEATUR_map            },
  { "mod_base",       FEATUR_mod_base       },
  { "number",         FEATUR_number         },
  { "organism",       FEATUR_organism       },
  { "PCR_conditions", FEATUR_PCR_conditions },
  { "phenotype",      FEATUR_phenotype      },
  { "product",        FEATUR_product_quals  },
  { "replace",        FEATUR_replace        },
  { "rpt_family",     FEATUR_rpt_family     },
  { "rpt_type",       FEATUR_rpt_type       },
  { "rpt_unit",       FEATUR_rpt_unit       },
  { "standard_name",  FEATUR_standard_name  },
  { "usedin",         FEATUR_usedin         }
};

static Int2 GbqualToFeaturIndex (
  CharPtr qualname
)

{
  Int2  L, R, mid;

  if (qualname == NULL || *qualname == '\0') return 0;

  L = 0;
  R = NUM_GB_QUALS - 1;

  while (L < R) {
    mid = (L + R) / 2;
    if (StringICmp (qualToFeature [mid].name, qualname) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  if (StringICmp (qualToFeature [R].name, qualname) == 0) {
    return qualToFeature [R].featurclass;
  }

  return 0;
}

#define NUM_ILLEGAL_QUALS 18

static FeaturQual illegalGbqualList [NUM_ILLEGAL_QUALS] = {
  { "allele",         Qual_class_quote   },
  { "anticodon",      Qual_class_noquote },
  { "citation",       Qual_class_noquote },
  { "codon",          Qual_class_noquote },
  { "codon_start",    Qual_class_noquote },
  { "cons_splice",    Qual_class_noquote },
  { "db_xref",        Qual_class_quote   },
  { "EC_number",      Qual_class_quote   },
  { "evidence",       Qual_class_noquote },
  { "exception",      Qual_class_quote   },
  { "frequency",      Qual_class_quote   },
  { "gene",           Qual_class_quote   },
  { "note",           Qual_class_quote   },
  { "protein_id",     Qual_class_quote   },
  { "pseudo",         Qual_class_noquote },  
  { "transl_except",  Qual_class_noquote },
  { "transl_table",   Qual_class_noquote },
  { "translation",    Qual_class_quote   },
};

static Int2 IllegalGbqualToClass (
  CharPtr qualname
)

{
  Int2  L, R, mid;

  if (qualname == NULL || *qualname == '\0') return 0;

  L = 0;
  R = NUM_ILLEGAL_QUALS - 1;

  while (L < R) {
    mid = (L + R) / 2;
    if (StringICmp (illegalGbqualList [mid].name, qualname) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  if (StringICmp (illegalGbqualList [R].name, qualname) == 0) {
    return illegalGbqualList [R].qualclass;
  }

  return 0;
}

static CharPtr trnaList [] = {
  "tRNA-Gap",
  "tRNA-Ala",
  "tRNA-Asx",
  "tRNA-Cys",
  "tRNA-Asp",
  "tRNA-Glu",
  "tRNA-Phe",
  "tRNA-Gly",
  "tRNA-His",
  "tRNA-Ile",
  "tRNA-Lys",
  "tRNA-Leu",
  "tRNA-Met",
  "tRNA-Asn",
  "tRNA-Pro",
  "tRNA-Gln",
  "tRNA-Arg",
  "tRNA-Ser",
  "tRNA-Thr",
  "tRNA-Sec",
  "tRNA-Val",
  "tRNA-Trp",
  "tRNA-Xxx",
  "tRNA-Tyr",
  "tRNA-Glx",
  NULL
};

static CharPtr evidenceText [] = {
  NULL, "experimental", "not_experimental"
};

static CharPtr SeqCodeNameGet (
  SeqCodeTablePtr table,
  Uint1 residue
)

{
  Uint1   index;
  CharPtr  retval = "?";

  if (table != NULL) {
    index = residue - table->start_at;
    if (index >= 0 && index < table->num) {
      retval = (table->names) [index];
    }
  }

  return retval;
}

static CharPtr Get3LetterSymbol (
  Uint1 seq_code,
  SeqCodeTablePtr table,
  Uint1 residue
)

{
  Int2             index;
  Uint1            new_residue;
  CharPtr          ptr;
  CharPtr          retval = NULL;
  SeqMapTablePtr   smtp;
  SeqCodeTablePtr  table_3aa;

  if (residue == 42) { /* stop codon in NCBIeaa */
    retval = "TERM";
    return retval;
  }

  if (seq_code != Seq_code_ncbieaa) {
    /* if code and seq_code are identical, then smtp is NULL?? */
    smtp = SeqMapTableFind (seq_code, Seq_code_ncbieaa);
    new_residue = SeqMapTableConvert (smtp, residue);
  } else {
    new_residue = residue;
  }

  /* The following looks for non-symbols (255) and "Undetermined" (88) */
  if ((int) new_residue == 255 || (int) new_residue == 88) {
    retval = "OTHER";
    return retval;
  } else {
    ptr = SeqCodeNameGet (table, residue);
    table_3aa=SeqCodeTableFind  (Seq_code_iupacaa3);
    if (ptr != NULL && table_3aa != NULL) {
      for (index=0; index < (int) table->num; index++) {
        if (StringCmp(ptr, (table_3aa->names) [index]) == 0) {
          retval = (table_3aa->symbols) [index];
          return retval;
        }
      }
    }
  }

  retval = "OTHER";
  return retval;
}

static Int2 MatchRef (
  ValNodePtr ppr,
  ReferenceBlockPtr PNTR rbpp,
  Int2 numReferences
)

{
  Char               buf [121];
  Int2               j;
  size_t             len;
  ReferenceBlockPtr  rbp;
  Int4               uid;

  if (ppr == NULL || rbpp == NULL) return 0;

  for (j = 0; j < numReferences; j++) {
    rbp = rbpp [j];
    if (rbp == NULL) continue;
    switch (ppr->choice) {
      case PUB_Muid :
        uid = ppr->data.intvalue;
        if (rbp->muid == uid) return j + 1;
        break;
      case PUB_PMid :
        uid = ppr->data.intvalue;
        if (rbp->pmid == uid) return j + 1;
        break;
      default :
        PubLabelUnique (ppr, buf, sizeof (buf) - 1, OM_LABEL_CONTENT, TRUE);
        len = StringLen (buf);
        if (len > 0 && buf [len - 1] == '>') {
          buf [len - 1] = '\0';
          len--;
        }
        len = MIN (len, StringLen (rbp->uniquestr));
        if (StringNICmp (rbp->uniquestr, buf, len) == 0) return j + 1;
        break;
    }
  }
  return 0;
}

static CharPtr TrimSpacesAndJunkFromEnds (
  CharPtr str
)

{
  Uchar    ch;	/* to use 8bit characters in multibyte languages */
  CharPtr  dst;
  CharPtr  ptr;

  if (str != NULL && str [0] != '\0') {
    dst = str;
    ptr = str;
    ch = *ptr;
    if (ch != '\0' && ch <= ' ') {
      while (ch != '\0' && ch <= ' ') {
        ptr++;
        ch = *ptr;
      }
      while (ch != '\0') {
        *dst = ch;
        dst++;
        ptr++;
        ch = *ptr;
      }
      *dst = '\0';
    }
    dst = NULL;
    ptr = str;
    ch = *ptr;
    while (ch != '\0') {
      if (ch == ' ' || ch == '.' || ch == '~' || ch == ';') {
        if (dst == NULL) {
          dst = ptr;
        }
      } else {
        dst = NULL;
      }
      ptr++;
      ch = *ptr;
    }
    if (dst != NULL) {
      *dst = '\0';
    }
  }
  return str;
}

static CharPtr siteList [] = {
  NULL,
  "active site",
  "binding site",
  "cleavage site",
  "inhibit site",
  "modified site",
  "glycosylation site",
  "myristoylation site",
  "mutagenized site",
  "metal-binding site",
  "phosphorylation site",
  "acetylation site",
  "amidation site",
  "methylation site",
  "hydroxylation site",
  "sulfatation site",
  "oxidative-deamination site",
  "pyrrolidone-carboxylic-acid site",
  "gamma-carboxyglutamic-acid site",
  "blocked site",
  "lipid-binding site",
  "np-binding site",
  "DNA binding site",
  "signal-peptide site",
  "transit-peptide site",
  "transmembrane-region site",
  "unclassified site"
};

static CharPtr conflict_msg =
"Protein sequence is in conflict with the conceptual translation";

static CharPtr no_protein_msg =
"Coding region translates with internal stops";

static CharPtr FormatFeatureBlock (
  Asn2gbFormatPtr afp,
  BaseBlockPtr bbp
)

{
  Uint1              aa;
  Asn2gbSectionPtr   asp;
  Boolean            at_end = FALSE;
  BioseqPtr          bsp;
  Char               buf [80];
  Choice             cbaa;
  CodeBreakPtr       cbp;
  SeqFeatPtr         cds;
  Char               ch;
  Uint1              choice;
  CdRegionPtr        crp;
  SeqMgrDescContext  dcontext;
  DbtagPtr           dbt;
  Int4               exp_ev;
  SeqMgrFeatContext  fcontext;
  Uint1              featdeftype;
  Uint1              from;
  GBQualPtr          gbq;
  SeqMgrFeatContext  gcontext;
  ValNodePtr         gcp;
  SeqFeatPtr         gene = NULL;
  GeneRefPtr         grp;
  ValNodePtr         head;
  Int2               i;
  IntCdsBlockPtr     icp;
  Uint1              idx;
  IntFeatBlockPtr    ifp;
  ValNodePtr         illegal = NULL;
  ImpFeatPtr         imp = NULL;
  Int2               j;
  Uint1              jdx;
  CharPtr            key;
  CharPtr            lasttype;
  Int4               len;
  SeqLocPtr          loc = NULL;
  SeqLocPtr          location = NULL;
  MolInfoPtr         mip;
  CharPtr            notestr;
  Uint1Ptr           notetbl = NULL;
  Char               numbuf [32];
  Int2               numcodons;
  ObjectIdPtr        oip;
  Uint2              partial;
  SeqMgrFeatContext  pcontext;
  ValNodePtr         ppr;
  CharPtr            prefix;
  BioseqPtr          prod = NULL;
  SeqFeatPtr         prot;
  Boolean            protein = FALSE;
  Char               protein_pid_g [32];
  CharPtr            protein_seq = NULL;
  ProtRefPtr         prp;
  Boolean            pseudo = FALSE;
  CharPtr            ptr;
  Int2               qualclass;
  Uint1Ptr           qualtbl = NULL;
  QualValPtr         qvp;
  Uint1              residue;
  RnaRefPtr          rrp;
  SeqCodeTablePtr    sctp;
  SeqDescrPtr        sdp;
  Uint1              seqcode;
  Char               seqid [50];
  SeqFeatPtr         sfp;
  Uint1              shift;
  SeqIntPtr          sintp;
  SeqIdPtr           sip;
  Int2               siteidx;
  SeqLocPtr          slp;
  SeqMapTablePtr     smtp;
  SeqPortPtr         spp;
  CharPtr            str;
  BioseqPtr          target;
  CharPtr            tmp;
  tRNAPtr            trna;
  ValNodePtr         vnp;
  CharPtr            wwwbuf;

  if (afp == NULL || bbp == NULL) return NULL;
  asp = afp->asp;
  if (asp == NULL) return NULL;
  target = asp->target;
  bsp = asp->bsp;
  if (target == NULL || bsp == NULL) return NULL;
  qvp = afp->qvp;
  if (qvp == NULL) return NULL;

  protein_pid_g [0] = '\0';

  /* all features in this list are known to be valid for the designated mode */

  sfp = SeqMgrGetDesiredFeature (bbp->entityID, NULL, bbp->itemID, 0, NULL, &fcontext);
  if (sfp == NULL) return NULL;

  /* may need to map location between aa and dna */

  ifp = (IntFeatBlockPtr) bbp;
  if (ifp->mapToNuc) {

    /* map mat_peptide, etc., to nucleotide coordinates */

    sip = SeqLocId (sfp->location);
    prod = BioseqFind (sip);
    cds = SeqMgrGetCDSgivenProduct (prod, NULL);
    location = aaFeatLoc_to_dnaFeatLoc (cds, sfp->location);
    loc = location;

  } else if (ifp->mapToProt) {
  } else {

    /* no aa-dna or dna-aa mapping, just use location */

    location = sfp->location;
  }
  if (location == NULL) return NULL;

  featdeftype = fcontext.featdeftype;
  if (featdeftype < FEATDEF_GENE || featdeftype >= FEATDEF_MAX) {
    featdeftype = FEATDEF_BAD;
  }
  key = featurekeys [featdeftype];

  /* deal with unmappable impfeats */

  if (featdeftype == FEATDEF_BAD && fcontext.seqfeattype == SEQFEAT_IMP) {
    imp = (ImpFeatPtr) sfp->data.value.ptrvalue;
    if (imp != NULL) {
      key = imp->key;
    }
  }

  gb_StartPrint (afp->format, TRUE, 5, 21, NULL, 0, 5, 21, "FT", FALSE);
  ff_AddString (key);
  TabToColumn (22);

  if (imp == NULL || StringHasNoText (imp->loc)) {
    str = FlatLoc (target, location, afp->style);
  } else {
    str = StringSave (imp->loc);
  }
  if (get_www ()) {
    wwwbuf = www_featloc (str);
    ff_AddString (wwwbuf);
    MemFree (wwwbuf);
  } else {
    ff_AddString (str);
  }
  MemFree (str);

  /* populate qualifier table from feature fields */

  if (sfp->partial) {
    partial = SeqLocPartialCheck (location);
    if (partial == SLP_COMPLETE /* || partial > SLP_OTHER */ ) {
      qvp [FEATUR_partial].bool = TRUE;
    }
    if (imp != NULL) {
      if (StringChr (imp->loc, '<') != NULL || StringChr (imp->loc, '>') != NULL) {
        qvp [FEATUR_partial].bool = FALSE;
      }
    }
  }

  if (sfp->pseudo) {
    pseudo = TRUE;
  }

  if (fcontext.seqfeattype == SEQFEAT_GENE) {
    grp = (GeneRefPtr) sfp->data.value.ptrvalue;
    if (grp != NULL) {
      if (! StringHasNoText (grp->locus)) {
        qvp [FEATUR_gene].str = grp->locus;
        qvp [FEATUR_gene_desc].str = grp->desc;
        qvp [FEATUR_gene_syn].vnp = grp->syn;
      } else if (! StringHasNoText (grp->desc)) {
        qvp [FEATUR_gene].str = grp->desc;
        qvp [FEATUR_gene_syn].vnp = grp->syn;
      } else if (grp->syn != NULL) {
        vnp = grp->syn;
        qvp [FEATUR_gene].str = (CharPtr) vnp->data.ptrvalue;
        vnp = vnp->next;
        qvp [FEATUR_gene_syn].vnp = vnp;
      }
      qvp [FEATUR_allele].str = grp->allele;
      qvp [FEATUR_gene_map].str = grp->maploc;
      if (grp->pseudo) {
        pseudo = TRUE;
      }
    }

  } else {

    grp = SeqMgrGetGeneXref (sfp);
    if (grp == NULL) {
      gene = SeqMgrGetOverlappingGene (location, &gcontext);
      if (gene != NULL) {
        qvp [FEATUR_gene_note].str = gene->comment;
        grp = (GeneRefPtr) gene->data.value.ptrvalue;
        if (gene->pseudo) {
          pseudo = TRUE;
        }
      }
    }
    if (grp != NULL && grp->pseudo) {
      pseudo = TRUE;
    }
    if (grp != NULL && (! SeqMgrGeneIsSuppressed (grp)) &&
        (fcontext.featdeftype != FEATDEF_repeat_region || gene == NULL)) {
      if (! StringHasNoText (grp->locus)) {
        qvp [FEATUR_gene].str = grp->locus;
      } else if (! StringHasNoText (grp->desc)) {
        qvp [FEATUR_gene].str = grp->desc;
      } else if (grp->syn != NULL) {
        vnp = grp->syn;
        qvp [FEATUR_gene].str = (CharPtr) vnp->data.ptrvalue;
      }
    }

    /* specific fields set here */

    switch (fcontext.seqfeattype) {
      case SEQFEAT_CDREGION :
        crp = (CdRegionPtr) sfp->data.value.ptrvalue;
        if (crp != NULL) {
          qvp [FEATUR_codon_start].num = crp->frame;
          if (qvp [FEATUR_codon_start].num == 0) {
            qvp [FEATUR_codon_start].num = 1;
          }
          qvp [FEATUR_transl_except].cbp = crp->code_break;
          gcp = crp->genetic_code;
          if (gcp != NULL) {
            for (vnp = gcp->data.ptrvalue; vnp != NULL; vnp = vnp->next) {
              if (vnp->choice == 2 && vnp->data.intvalue != 0) {
                qvp [FEATUR_transl_table].num = vnp->data.intvalue;
              }
            }

            /* suppress table 1 */

            if (qvp [FEATUR_transl_table].num == 1) {
              qvp [FEATUR_transl_table].num = 0;
            }
          }

          if (sfp->product != NULL && SeqLocLen (sfp->product) != 0) {
            protein = TRUE;
          }
          if (crp->conflict && (protein || (! sfp->excpt))) {
            if (protein) {
              qvp [FEATUR_prot_conflict].str = conflict_msg;
            } else {
              qvp [FEATUR_prot_missing].str = no_protein_msg;
            }
          }
        }

        prp = SeqMgrGetProtXref (sfp);
        if (prp == NULL) {
          sip = SeqLocId (sfp->product);
          qvp [FEATUR_protein_id].sip = sip;
          prod = BioseqFind (sip);
          if (prod != NULL) {
            for (sip = prod->id; sip != NULL; sip = sip->next) {
              if (sip->choice == SEQID_GI) {
                sprintf (protein_pid_g, "PID:g%ld", (long) sip->data.intvalue);
              }
            }
            sdp = SeqMgrGetNextDescriptor (prod, NULL, Seq_descr_comment, &dcontext);
            if (sdp != NULL && dcontext.level == 0) {
              if (! StringHasNoText ((CharPtr) sdp->data.ptrvalue)) {
                qvp [FEATUR_prot_comment].str = (CharPtr) sdp->data.ptrvalue;
              }
            }
            sdp = SeqMgrGetNextDescriptor (prod, NULL, Seq_descr_molinfo, &dcontext);
            if (sdp != NULL && dcontext.level == 0) {
              mip = (MolInfoPtr) sdp->data.ptrvalue;
              if (mip != NULL && mip->tech > 1 &&
                  mip->tech != MI_TECH_concept_trans &&
                  mip->tech != MI_TECH_concept_trans_a) {
                str = StringForSeqTech (mip->tech);
                if (! StringHasNoText (str)) {
                  qvp [FEATUR_prot_method].str = str;
                }
              }
            }
            prot = SeqMgrGetBestProteinFeature (prod, &pcontext);
            if (prot != NULL) {
              prp = (ProtRefPtr) prot->data.value.ptrvalue;
              qvp [FEATUR_prot_note].str = prot->comment;
            }
          }
        }

        if (prp != NULL) {
          vnp = prp->name;
          if (vnp != NULL && (! StringHasNoText ((CharPtr) vnp->data.ptrvalue))) {
            qvp [FEATUR_cds_product].str = (CharPtr) vnp->data.ptrvalue;
            vnp = vnp->next;
            qvp [FEATUR_prot_names].vnp = vnp;
          }
          qvp [FEATUR_prot_desc].str = prp->desc;
          qvp [FEATUR_prot_activity].vnp = prp->activity;
          qvp [FEATUR_EC_number].vnp = prp->ec;
        }

        if (! pseudo) {
          qvp [FEATUR_translation].bool = TRUE;
        }

        if (ifp->isCDS) {
          icp = (IntCdsBlockPtr) ifp;
          qvp [FEATUR_figure].str = icp->fig;
          qvp [FEATUR_maploc].str = icp->maploc;
        }
        break;
      case SEQFEAT_PROT :
        prp = (ProtRefPtr) sfp->data.value.ptrvalue;
        if (prp != NULL) {
          vnp = prp->name;
          if (vnp != NULL && (! StringHasNoText ((CharPtr) vnp->data.ptrvalue))) {
            qvp [FEATUR_product].str = (CharPtr) vnp->data.ptrvalue;
            vnp = vnp->next;
            qvp [FEATUR_prot_names].vnp = vnp;
          }
          qvp [FEATUR_prot_desc].str = prp->desc;
          qvp [FEATUR_prot_activity].vnp = prp->activity;
          qvp [FEATUR_EC_number].vnp = prp->ec;
        }
        break;
      case SEQFEAT_RNA :
        rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
        if (rrp != NULL) {
          if (rrp->type == 3) {
            if (rrp->ext.choice == 2) {
              trna = (tRNAPtr) rrp->ext.value.ptrvalue;
              if (trna != NULL) {
                aa = 0;
                if (trna->aatype == 2) {
                  aa = trna->aa;
                } else {
                  from = 0;
                  switch (trna->aatype) {
                    case 0 :
                      from = 0;
                      break;
                    case 1 :
                      from = Seq_code_iupacaa;
                      break;
                    case 2 :
                      from = Seq_code_ncbieaa;
                      break;
                    case 3 :
                      from = Seq_code_ncbi8aa;
                      break;
                    case 4 :
                      from = Seq_code_ncbistdaa;
                      break;
                    default:
                      break;
                  }
                  smtp = SeqMapTableFind (Seq_code_ncbieaa, from);
                  if (smtp != NULL) {
                    aa = SeqMapTableConvert (smtp, trna->aa);
                  }
                }
                if (aa > 0 && aa != 255) {
                  if (aa <= 74) {
                    shift = 0;
                  } else if (aa > 79) {
                    shift = 2;
                  } else {
                    shift = 1;
                  } 
                  idx = aa - (64 + shift);
                  if (idx > 0 && idx < 25) {
                    str = trnaList [idx];
                    qvp [FEATUR_product].str = str;
                    if (StringNICmp (str, "tRNA-", 5) == 0) {
                      qvp [FEATUR_trna_aa].str = str + 5;
                    }
                  }
                }
                qvp [FEATUR_anticodon].slp = trna->anticodon;
                qvp [FEATUR_trna_codons].trp = trna;
              }
            }
          } else {
            if (rrp->ext.choice == 1) {
              str = (CharPtr) rrp->ext.value.ptrvalue;
              qvp [FEATUR_product].str = str;

              /* !!! no preRNA product allowed is probably oversight !!! */

              if (fcontext.featdeftype == FEATDEF_preRNA) {
                qvp [FEATUR_product].str = NULL;
              }

              /*
              if (rrp->type == 255 && (! StringHasNoText (str))) {
                if        (StringICmp (str, "internal transcribed spacer 1") == 0 ||
                           StringICmp (str, "internal transcribed spacer ITS1") == 0 ||
                           StringICmp (str, "ITS1") == 0) {
                  qvp [FEATUR_rrna_its].str = "ITS1";
                } else if (StringICmp (str, "internal transcribed spacer 2") == 0 ||
                           StringICmp (str, "internal transcribed spacer ITS2") == 0 ||
                           StringICmp (str, "ITS2") == 0) {
                  qvp [FEATUR_rrna_its].str = "ITS2";
                } else if (StringICmp (str, "internal transcribed spacer 3") == 0 ||
                           StringICmp (str, "internal transcribed spacer ITS3") == 0 ||
                           StringICmp (str, "ITS3") == 0) {
                  qvp [FEATUR_rrna_its].str = "ITS3";
                }
              }
              */
            }
          }
        }
        break;
      case SEQFEAT_REGION :
        qvp [FEATUR_region].str = (CharPtr) sfp->data.value.ptrvalue;
        break;
      case SEQFEAT_COMMENT :
        break;
      case SEQFEAT_SITE :
        siteidx = (Int2) sfp->data.value.intvalue;
        if (siteidx == 255) {
          siteidx = 26;
        }
        if (siteidx > 0 && siteidx < 27) {
          qvp [FEATUR_site].str = siteList [siteidx];
        }
        break;
      default :
        break;
    }
  }

  /* common fields set here */

  if (fcontext.featdeftype == FEATDEF_repeat_region) {
    pseudo = FALSE;
  }

  qvp [FEATUR_pseudo].bool = pseudo;

  qvp [FEATUR_evidence].num = sfp->exp_ev;
  qvp [FEATUR_exception].str = sfp->except_text;
  if (sfp->excpt && qvp [FEATUR_exception].str == NULL) {
    qvp [FEATUR_exception].str = "No explanation supplied";
  }

  qvp [FEATUR_db_xref].vnp = sfp->dbxref;
  qvp [FEATUR_citation].vnp = sfp->cit;

  qvp [FEATUR_seqfeat_note].str = sfp->comment;

  /* /product same as sfp->comment will suppress /note */

  if (! StringHasNoText (qvp [FEATUR_product].str) &&
      StringICmp (sfp->comment, qvp [FEATUR_product].str) == 0) {
    qvp [FEATUR_seqfeat_note].str = NULL;
  }
  if (! StringHasNoText (qvp [FEATUR_cds_product].str) &&
      StringICmp (sfp->comment, qvp [FEATUR_cds_product].str) == 0) {
    qvp [FEATUR_seqfeat_note].str = NULL;
  }

  /* /gene same as sfp->comment will suppress /note */
  /* case sensitive -gi|6572973|gb|AF195052.1|AF195052 */

  if (! StringHasNoText (qvp [FEATUR_gene].str) &&
      StringCmp (sfp->comment, qvp [FEATUR_gene].str) == 0) {
    qvp [FEATUR_seqfeat_note].str = NULL;
  }

  /* gene /note same as sfp->comment will suppress /note */

  if (! StringHasNoText (qvp [FEATUR_gene_note].str) &&
      StringICmp (sfp->comment, qvp [FEATUR_gene_note].str) == 0) {
    qvp [FEATUR_seqfeat_note].str = NULL;
  }

  /* if site sfp->comment contains site name, suppress site in /note */

  if (! StringHasNoText (qvp [FEATUR_site].str) &&
      StringStr (sfp->comment, qvp [FEATUR_site].str) != NULL) {
    qvp [FEATUR_site].str = NULL;
  }

  /* now go through gbqual list */

  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    idx = GbqualToFeaturIndex (gbq->qual);
    if (idx > 0 && idx < ASN2GNBK_TOTAL_FEATUR) {
      if (qvp [idx].gbq == NULL) {
        if (idx == FEATUR_product_quals) {
          if (qvp [FEATUR_product].str == NULL) {
            qvp [FEATUR_product].str = gbq->val;
          } else {
            qvp [idx].gbq = gbq;
          }
        } else {
          qvp [idx].gbq = gbq;
        }
      }

    } else if (idx == 0) {

      qualclass = IllegalGbqualToClass (gbq->qual);
      tmp = StringSave (gbq->val);
      if (tmp != NULL) {
        str = MemNew (sizeof (Char) * (StringLen (gbq->val) + StringLen (tmp) + 10));
        if (str != NULL) {
          if (qualclass == Qual_class_quote) {
	        if (StringIsJustQuotes (tmp)) {
              sprintf (buf, "/%s", gbq->qual);
	        } else {
              ptr = tmp;
              ch = *ptr;
              while (ch != '\0') {
                if (ch == '"') {
                  *ptr = '\'';
                }
                ptr++;
                ch = *ptr;
              }
              sprintf (str, "/%s=\"%s\"", gbq->qual, tmp);
	        }
            ValNodeCopyStr (&illegal, 0, str);
          } else if (qualclass == Qual_class_noquote) {
	        if (StringIsJustQuotes (tmp)) {
              sprintf (str, "/%s", gbq->qual);
	        } else {
              sprintf (str, "/%s=%s", gbq->qual, tmp);
	        }
            ValNodeCopyStr (&illegal, 0, str);
          }
          MemFree (str);
        }
        MemFree (tmp);
      }
    }
  }

  /* illegal qualifiers are copied and formatted in valnode chain */

  qvp [FEATUR_illegal_qual].vnp = illegal;

  /* remove protein description that equals the gene name, case sensitive */

  if (StringCmp (qvp [FEATUR_gene].str, qvp [FEATUR_prot_desc].str) == 0) {
    qvp [FEATUR_prot_desc].str = NULL;
  }

  /* remove protein description that equals the cds product, case sensitive */

  if (StringCmp (qvp [FEATUR_cds_product].str, qvp [FEATUR_prot_desc].str) == 0) {
    qvp [FEATUR_prot_desc].str = NULL;
  }

  /* remove comment contained in prot_desc - gi|4530123|gb|AF071539.1|AF071539 */

  if (StringStr (qvp [FEATUR_prot_desc].str, qvp [FEATUR_seqfeat_note].str) != NULL) {
    qvp [FEATUR_seqfeat_note].str = NULL;
  }

  /* remove protein description that equals the standard name */

  if (qvp [FEATUR_standard_name].gbq != NULL && qvp [FEATUR_prot_desc].str != NULL) {
    gbq = qvp [FEATUR_standard_name].gbq;
    lasttype = gbq->qual;
    while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
      if (StringICmp (gbq->val, qvp [FEATUR_prot_desc].str) == 0) {
        qvp [FEATUR_prot_desc].str = NULL;
      }
      gbq = gbq->next;
    }
  }

  /* remove protein comment descriptor that equals the protein note */

  if (StringCmp (qvp [FEATUR_prot_note].str, qvp [FEATUR_prot_comment].str) == 0) {
    qvp [FEATUR_prot_comment].str = NULL;
  }

  /* now print qualifiers from table */

  if (afp->mode == RELEASE_MODE) {
    qualtbl = relmode_feat_qual_order;
    notetbl = relmode_feat_note_order;
  } else {
    qualtbl = seqmode_feat_qual_order;
    notetbl = seqmode_feat_note_order;
  }

  for (i = 0, idx = qualtbl [i]; idx != 0; i++, idx = qualtbl [i]) {

    lasttype = NULL;
    switch (asn2gnbk_featur_quals [idx].qualclass) {

      case Qual_class_ignore :
        break;

      case Qual_class_string :
        if (! StringHasNoText (qvp [idx].str)) {
          sprintf (buf, "/%s=\"", asn2gnbk_featur_quals [idx].name);
          NewContLine ();
          gb_AddString (buf, qvp [idx].str, "\"", FALSE, TRUE, FALSE);
        }
        break;

      case Qual_class_boolean :
        if (qvp [idx].bool) {
          sprintf (buf, "/%s", asn2gnbk_featur_quals [idx].name);
          NewContLine ();
          ff_AddString (buf);
        }
        break;

      case Qual_class_int :
        if (qvp [idx].num > 0) {
          sprintf (numbuf, "/%s=%ld", asn2gnbk_featur_quals [idx].name, (long) qvp [idx].num);
          NewContLine ();
          ff_AddString (numbuf);
        }
        break;

      case Qual_class_evidence :
        exp_ev = qvp [idx].num;
        if (exp_ev > 0 && exp_ev <= 2) {
          sprintf (numbuf, "/%s=%s", asn2gnbk_featur_quals [idx].name, evidenceText [exp_ev]);
          NewContLine ();
          ff_AddString (numbuf);
        }
        break;

      case Qual_class_valnode :
        for (vnp = qvp [idx].vnp; vnp != NULL; vnp = vnp->next) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (str != NULL) {
            sprintf (buf, "/%s=\"", asn2gnbk_featur_quals [idx].name);
            NewContLine ();
            gb_AddString (buf, str, "\"", FALSE, TRUE, FALSE);
          }
        }
        break;

      case Qual_class_quote :
        gbq = qvp [idx].gbq;
        if (lasttype == NULL && gbq != NULL) {
          lasttype = gbq->qual;
        }
        while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
          if (! StringHasNoText (gbq->val)) {
            sprintf (buf, "/%s=\"", asn2gnbk_featur_quals [idx].name);
            NewContLine ();
            if (StringIsJustQuotes (gbq->val)) {
              gb_AddString (buf, NULL, "\"", FALSE, TRUE, FALSE);
            } else {
              gb_AddString (buf, gbq->val, "\"", FALSE, TRUE, FALSE);
            }
          }
          gbq = gbq->next;
        }
        break;

      case Qual_class_noquote :
        gbq = qvp [idx].gbq;
        if (lasttype == NULL && gbq != NULL) {
          lasttype = gbq->qual;
        }
        while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
          if (! StringHasNoText (gbq->val)) {
            sprintf (buf, "/%s=", asn2gnbk_featur_quals [idx].name);
            NewContLine ();
            gb_AddString (buf, gbq->val, NULL, FALSE, TRUE, FALSE);
          }
          gbq = gbq->next;
        }
        break;

      case Qual_class_paren :
        gbq = qvp [idx].gbq;
        if (lasttype == NULL && gbq != NULL) {
          lasttype = gbq->qual;
        }
        while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
          if (! StringHasNoText (gbq->val)) {
            tmp = StringSave (gbq->val);
            str = tmp;
            if (*str == '(') {
              str++;
            }
            while (! StringHasNoText (str)) {
              ptr = StringChr (str, ',');
              if (ptr == NULL) {
                ptr = StringChr (str, ')');
              }
              if (ptr != NULL) {
                *ptr = '\0';
                ptr++;
              }
              sprintf (buf, "/%s=", asn2gnbk_featur_quals [idx].name);
              NewContLine ();
              gb_AddString (buf, str, NULL, FALSE, TRUE, FALSE);
              str = ptr;
            }
            MemFree (tmp);
          }
          gbq = gbq->next;
        }
        break;

      case Qual_class_region :
        break;

      case Qual_class_replace :
        gbq = qvp [idx].gbq;
        if (lasttype == NULL && gbq != NULL) {
          lasttype = gbq->qual;
        }
        while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
          if (! StringHasNoText (gbq->val)) {
            sprintf (buf, "/%s=\"", asn2gnbk_featur_quals [idx].name);
            NewContLine ();
            gb_AddString (buf, gbq->val, "\"", FALSE, TRUE, FALSE);
          } else {
            sprintf (buf, "/%s=\"", asn2gnbk_featur_quals [idx].name);
            NewContLine ();
            gb_AddString (buf, NULL, "\"", FALSE, TRUE, FALSE);
          }
          gbq = gbq->next;
        }
        break;

      case Qual_class_L_R_B :
        gbq = qvp [idx].gbq;
        if (lasttype == NULL && gbq != NULL) {
          lasttype = gbq->qual;
        }
        while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
          if (! StringHasNoText (gbq->val)) {
            sprintf (buf, "/%s=", asn2gnbk_featur_quals [idx].name);
            NewContLine ();
            gb_AddString (buf, gbq->val, NULL, FALSE, TRUE, FALSE);
          }
          gbq = gbq->next;
        }
        break;

      case Qual_class_code_break :
        cbp = qvp [idx].cbp;
        seqcode = 0;
        sctp = NULL;
        while (cbp != NULL) {
          cbaa = cbp->aa;
          switch (cbaa.choice) {
            case 1 :
              seqcode = Seq_code_ncbieaa;
              break;
            case 2 :
              seqcode = Seq_code_ncbi8aa;
              break;
            case 3 :
              seqcode = Seq_code_ncbistdaa;
              break;
            default :
              break;
          }
          if (seqcode != 0) {
            sctp = SeqCodeTableFind (seqcode);
            if (sctp != NULL) {
              slp = NULL;
              while ((slp = SeqLocFindNext (cbp->loc, slp)) != NULL) {
                str = FlatLoc (target, slp, afp->style);
                if (str != NULL) {
                  residue = cbaa.value.intvalue;
                  ptr = Get3LetterSymbol (seqcode, sctp, residue);
                  if (ptr == NULL) {
                    ptr = "OTHER";
                  }
                  sprintf (buf, "/transl_except=(pos:%s,aa:%s)", str, ptr);
                  NewContLine ();
                  ff_AddString (buf);
                }
                MemFree (str);
              }
            }
          }
          cbp = cbp->next;
        }
        break;

      case Qual_class_anti_codon :
        slp = qvp [FEATUR_anticodon].slp;
        if (slp != NULL && slp->choice == SEQLOC_INT) {
          sintp = (SeqIntPtr) slp->data.ptrvalue;
          if (sintp != NULL) {
            str = qvp [FEATUR_trna_aa].str;
            if (! StringHasNoText (str)) {
              sprintf (buf, "/anticodon=(pos:%ld..%ld,aa:%s)",
                       (long) sintp->from + 1, (long) sintp->to + 1, str);
              NewContLine ();
              ff_AddString (buf);
            }
          }
        }
        break;

      case Qual_class_codon :
        break;

      case Qual_class_pubset :
        vnp = qvp [idx].vnp;
        if (vnp != NULL && asp->referenceArray != NULL) {
          for (ppr = vnp->data.ptrvalue; ppr != NULL; ppr = ppr->next) {
            j = MatchRef (ppr, asp->referenceArray, asp->numReferences);
            if (j > 0) {
              sprintf (buf, "%d", (int) j);
              NewContLine ();
              gb_AddString ("/citation=[", buf, "]", FALSE, TRUE, FALSE);
            }
          }
        }
        break;

      case Qual_class_db_xref :
        for (vnp = qvp [idx].vnp; vnp != NULL; vnp = vnp->next) {
          buf [0] = '\0';
          dbt = (DbtagPtr) vnp->data.ptrvalue;
          if (dbt != NULL && (! StringHasNoText (dbt->db))) {
            oip = dbt->tag;
            if (oip != NULL) {
              /* if release mode, drop unknown dbtag */
              if (! StringHasNoText (oip->str)) {
                if (StringLen (dbt->db) + StringLen (oip->str) < 80) {
                  sprintf (buf, "%s:%s", dbt->db, oip->str);
                }
              } else {
                sprintf (buf, "%s:%ld", dbt->db, (long) oip->id);
              }
            }
          }
          if (! StringHasNoText (buf)) {
            if (StringICmp (buf, protein_pid_g) != 0) {
              NewContLine ();
              gb_AddString ("/db_xref=\"", buf, "\"", FALSE, TRUE, FALSE);
            }
          }
        }
        break;

      case Qual_class_seq_id :
        sip = qvp [idx].sip;
        if (sip != NULL) {
          prod = BioseqFind (sip);
          if (prod != NULL) {
            choice = 0;
            for (sip = prod->id; sip != NULL; sip = sip->next) {
              if (sip->choice == SEQID_GENBANK ||
                  sip->choice == SEQID_EMBL ||
                  sip->choice == SEQID_DDBJ) {
                choice = sip->choice;
                if (SeqIdWrite (sip, seqid, PRINTID_TEXTID_ACC_VER, sizeof (seqid)) != NULL) {
                  sprintf (buf, "/%s=\"", asn2gnbk_featur_quals [idx].name);
                  NewContLine ();
                  gb_AddString (buf, seqid, "\"", FALSE, TRUE, FALSE);
                }
              }
            }
            for (sip = prod->id; sip != NULL; sip = sip->next) {
              if (sip->choice == SEQID_GI) {
                /*
                if (choice != 0) {
                  sprintf (seqid, "PID:g%ld", (long) sip->data.intvalue);
                  NewContLine ();
                  gb_AddString ("/db_xref=\"", seqid, "\"", FALSE, TRUE, FALSE);
                }
                */
                sprintf (seqid, "GI:%ld", (long) sip->data.intvalue);
                NewContLine ();
                gb_AddString ("/db_xref=\"", seqid, "\"", FALSE, TRUE, FALSE);
              } else if (sip->choice == SEQID_GENERAL) {
                dbt = (DbtagPtr) sip->data.ptrvalue;
                if (dbt != NULL && StringCmp (dbt->db, "PID") == 0) {
                  /*
                  oip = dbt->tag;
                  if (oip != NULL) {
                    if (! StringHasNoText (oip->str)) {
                      sprintf (seqid, "PID:%s", oip->str);
                      NewContLine ();
                      gb_AddString ("/db_xref=\"", seqid, "\"", FALSE, TRUE, FALSE);
                    }
                  }
                  */
                }
              }
            }
          }
        }
        break;

      case Qual_class_trna_codons :
        break;

      case Qual_class_translation :
        if (qvp [idx].bool) {
          if (prod != NULL) {
            len = SeqLocLen (sfp->product);
            if (len > 0) {
              if (SeqLocStart (location) == 0 || SeqLocStop (location) == bsp->length - 1) {
                at_end = TRUE;
              }
              str = (CharPtr) MemNew ((size_t) (len + 1) * sizeof (Char));
              protein_seq = str;
              spp = SeqPortNewByLoc (sfp->product, Seq_code_ncbieaa);
              if (spp != NULL) {
                spp->do_virtual = TRUE;
                while ((residue = SeqPortGetResidue (spp)) != SEQPORT_EOF) {
                  if (! (IS_residue (residue))) continue;
                  if (residue == INVALID_RESIDUE) {
                    residue = (Uint1) 'X';
                  }
                  *protein_seq = residue;
                  protein_seq++;
                }
                if (at_end && StringLen (str) < 0) {
                  str = MemFree (str);
                }
                if (! StringHasNoText (str)) {
                  NewContLine ();
                  gb_AddString ("/translation=\"", str, "\"", FALSE, TRUE, FALSE);
                }
                MemFree (str);
              }
              SeqPortFree (spp);
            }
          /*
          } else {
            bs = ProteinFromCdRegionEx (sfp, TRUE, FALSE);
            if (bs != NULL) {
              str = BSMerge (bs, NULL);
              bs = BSFree (bs);
              if (str != NULL) {
                ptr = str;
                ch = *ptr;
                while (ch != '\0') {
                  *ptr = TO_UPPER (ch);
                  ptr++;
                  ch = *ptr;
                }
                if (! StringHasNoText (str)) {
                  NewContLine ();
                  gb_AddString ("/translation=\"", str, "\"", FALSE, TRUE, FALSE);
                }
                MemFree (str);
              }
            }
            */
          }
        }
        break;

      case Qual_class_illegal :
        for (vnp = qvp [idx].vnp; vnp != NULL; vnp = vnp->next) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (str != NULL) {
            NewContLine ();
            gb_AddString (NULL, str, NULL, FALSE, FALSE, FALSE);
          }
        }
        break;

      case Qual_class_note :
        head = NULL;
        notestr = NULL;
        prefix = NULL;
        for (j = 0, jdx = notetbl [j]; jdx != 0; j++, jdx = notetbl [j]) {
          switch (asn2gnbk_featur_quals [jdx].qualclass) {
            case Qual_class_string :
              if (! StringHasNoText (qvp [jdx].str)) {
                if (jdx == FEATUR_figure) {
                  sprintf (buf, "This sequence comes from %s", qvp [jdx].str);
                  AddValNodeString (&head, prefix, buf, NULL);
                } else if (jdx == FEATUR_maploc) {
                  sprintf (buf, "Map location %s", qvp [jdx].str);
                  AddValNodeString (&head, prefix, buf, NULL);
                } else if (jdx == FEATUR_seqfeat_note || jdx == FEATUR_prot_note) {
                  str = StringSave (qvp [jdx].str);
                  TrimSpacesAndJunkFromEnds (str);
                  AddValNodeString (&head, prefix, str, NULL);
                  MemFree (str);
                } else {
                  AddValNodeString (&head, prefix, qvp [jdx].str, NULL);
                }
                prefix = "; ";
              }
              break;
            case Qual_class_method :
              if (! StringHasNoText (qvp [jdx].str)) {
                sprintf (buf, "Method: %s", qvp [jdx].str);
                AddValNodeString (&head, prefix, buf, NULL);
                prefix = "; ";
              }
              break;
            case Qual_class_valnode :
              for (vnp = qvp [jdx].vnp; vnp != NULL; vnp = vnp->next) {
                str = (CharPtr) vnp->data.ptrvalue;
                if (! StringHasNoText (str)) {
                  AddValNodeString (&head, prefix, str, NULL);
                  prefix = "; ";
                }
              }
              break;
            case Qual_class_region :
              if (! StringHasNoText (qvp [jdx].str)) {
                if (head == NULL) {
                  prefix = "Region: ";
                } else {
                  prefix = "; Region: ";
                }
                AddValNodeString (&head, prefix, qvp [jdx].str, NULL);
                prefix = "; ";
              }
              break;
            case Qual_class_protnames :
              for (vnp = qvp [jdx].vnp; vnp != NULL; vnp = vnp->next) {
                str = (CharPtr) vnp->data.ptrvalue;
                if (! StringHasNoText (str)) {
                  /* case sensitive - gi|4973426|gb|AF148501.1|AF148501 */
                  if (StringCmp (qvp [FEATUR_gene].str, str) != 0) {
                    AddValNodeString (&head, prefix, str, NULL);
                    prefix = "; ";
                  }
                }
              }
              break;
            case Qual_class_its :
              str = qvp [jdx].str;
              if (! StringHasNoText (str)) {
                if (sfp->comment == NULL || StringStr (sfp->comment, str) == NULL) {
                  AddValNodeString (&head, prefix, str, NULL);
                  prefix = "; ";
                }
              }
              break;
            case Qual_class_trna_codons :
              trna = qvp [jdx].trp;
              if (trna) {
                numcodons = ComposeCodonsRecognizedString (trna, numbuf, sizeof (numbuf));
                if (numcodons < 1 || StringHasNoText (numbuf)) {
                } else if (numcodons == 1) {
                  AddValNodeString (&head, prefix, "codon recognized: ", numbuf);
                  prefix = "; ";
                } else {
                  AddValNodeString (&head, prefix, "codons recognized: ", numbuf);
                  prefix = "; ";
                }
              }
              break;
            default :
              break;
          }
        }
        if (head != NULL) {
          notestr = MergeValNodeStrings (head);
          NewContLine ();
          gb_AddString ("/note=\"", notestr, "\"", FALSE, TRUE, FALSE);
          MemFree (notestr);
          ValNodeFreeData (head);
        }
        break;
      default :
        break;
    }
  }

  /* and then deal with the various note types separately (not in order table) */

  /* free aa-dna or dna-aa mapped location */

  SeqLocFree (loc);

  ValNodeFreeData (illegal);

  return gb_MergeString (TRUE);
}

static CharPtr FormatBasecountBlock (
  Asn2gbFormatPtr afp,
  BaseBlockPtr bbp
)

{
  Asn2gbJobPtr      ajp;
  Asn2gbSectionPtr  asp;
  Int4              base_count [5];
  BioseqPtr         bsp;
  Char              buf [80];
  Byte              bases [400];
  Uint1             code = Seq_code_iupacna;
  Int2              ctr;
  Int2              i;
  Int4              len;
  Uint1             residue;
  SeqPortPtr        spp = NULL;
  Int4              total = 0;

  if (afp == NULL || bbp == NULL) return NULL;
  ajp = afp->ajp;
  if (ajp == NULL) return NULL;

  asp = afp->asp;
  if (asp == NULL) return NULL;
  bsp = (asp->bsp);
  if (bsp == NULL) return NULL;

  /* after first formatting, result is cached into bbp->string */

  if (! StringHasNoText (bbp->string)) return StringSave (bbp->string);

  for (i = 0; i < 5; i++) {
    base_count [i] = 0;
  }

  if (ISA_aa (bsp->mol)) {
    code = Seq_code_ncbieaa;
  }

  if (ajp->slp != NULL) {
    spp = SeqPortNewByLoc (ajp->slp, code);
    len = SeqLocLen (ajp->slp);
  } else {
    spp = SeqPortNew (bsp, 0, -1, 0, code);
    len = bsp->length;
  }
  if (bsp->repr == Seq_repr_delta || bsp->repr == Seq_repr_virtual) {
    SeqPortSet_do_virtual (spp, TRUE);
  }

  /* use SeqPortRead rather than SeqPortGetResidue for faster performance */

  ctr = SeqPortRead (spp, bases, sizeof (bases));
  i = 0;
  residue = (Uint1) bases [i];
  while (residue != SEQPORT_EOF) {
    if (IS_residue (residue)) {
      total++;
      switch (residue) {
        case 'A' :
          (base_count [0])++;
          break;
        case 'C' :
          (base_count [1])++;
          break;
        case 'G' :
          (base_count [2])++;
          break;
        case 'T' :
          (base_count [3])++;
          break;
        default :
          (base_count [4])++;
          break;
      }
    }
    i++;
    if (i >= ctr) {
      i = 0;
      ctr = SeqPortRead (spp, bases, sizeof (bases));
      if (ctr < 0) {
        bases [0] = -ctr;
      } else if (ctr < 1) {
        bases [0] = SEQPORT_EOF;
      }
    }
    residue = (Uint1) bases [i];
  }

  SeqPortFree (spp);

  init_buff ();

  if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {

    if (base_count [4] == 0) {
      sprintf (buf, "%7ld a%7ld c%7ld g%7ld t",
               (long) base_count [0], (long) base_count [1], 
               (long) base_count [2], (long) base_count [3]);
    } else {
      sprintf (buf, "%7ld a%7ld c%7ld g%7ld t%7ld others",
               (long) base_count [0], (long) base_count [1], 
               (long) base_count [2], (long) base_count [3],
               (long) base_count [4]);
    }

  } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {

    sprintf (buf, "Sequence %ld BP; %ld A; %ld C; %ld G; %ld T; %ld other;",
             (long) len,
             (long) base_count [0], (long) base_count [1], 
             (long) base_count [2], (long) base_count [3],
             (long) base_count [4]);

    PrintXX ();
  }

  gb_StartPrint (afp->format, FALSE, 0, 0, "BASE COUNT", 13, 5, 5, "SQ", FALSE);

  gb_AddString (NULL, buf, NULL, FALSE, FALSE, FALSE);

  return gb_MergeString (TRUE);
}

static void PrintSeqLine (
  FmtType format,
  CharPtr buf,
  Int4 start,
  Int4 stop
)

{
  Char  pos [16];

  if (format == GENBANK_FMT || format == GENPEPT_FMT) {

    sprintf (pos, "%9ld", (long) (start + 1));
    ff_StartPrint (0, 0, ASN2FF_GB_MAX, NULL);
    gb_AddString (NULL, pos, NULL, FALSE, FALSE, FALSE);
    TabToColumn (11);
    gb_AddString (NULL, buf, NULL, FALSE, FALSE, FALSE);
    ff_EndPrint();

  } else if (format == EMBL_FMT || format == EMBLPEPT_FMT) {

    sprintf (pos, "%8ld", (long) (stop));
    ff_StartPrint (5, 5, 0, NULL);
    gb_AddString (NULL, buf, NULL, FALSE, FALSE, FALSE);
    TabToColumn (73);
    gb_AddString (NULL, pos, NULL, FALSE, FALSE, FALSE);
    ff_EndPrint();
  }
}

static CharPtr FormatSequenceBlock (
  Asn2gbFormatPtr afp,
  BaseBlockPtr bbp
)

{
  Asn2gbJobPtr      ajp;
  Asn2gbSectionPtr  asp;
  Byte              bases [400];
  Int2              blk;
  BioseqPtr         bsp;
  Char              buf [80];
  Int2              cnt;
  Uint1             code = Seq_code_iupacna;
  Int2              count;
  Int2              ctr;
  Int2              i;
  Boolean           is_na;
  Int2              lin;
  Int4              pos;
  Uint1             residue;
  SequenceBlockPtr  sbp;
  SeqPortPtr        spp;
  Int4              start;
  Int4              stop;

  if (afp == NULL || bbp == NULL) return NULL;
  sbp = (SequenceBlockPtr) bbp;
  ajp = afp->ajp;
  if (ajp == NULL) return NULL;

  asp = afp->asp;
  if (asp == NULL) return NULL;
  bsp = (asp->bsp);
  if (bsp == NULL) return NULL;

  spp = asp->spp;
  if (spp == NULL) {

    /* if first time, create SeqPort for this section */

    if (ISA_aa (bsp->mol)) {
      code = Seq_code_ncbieaa;
    }

    if (ajp->slp != NULL) {
      spp = SeqPortNewByLoc (ajp->slp, code);
    } else {
      spp = SeqPortNew (bsp, 0, -1, 0, code);
    }
    if (bsp->repr == Seq_repr_delta || bsp->repr == Seq_repr_virtual) {
      SeqPortSet_do_virtual (spp, TRUE);
    }

    asp->spp = spp;
  }

  init_buff ();

  start = sbp->start;
  stop = sbp->stop;

  if (start != spp->curpos) {
    SeqPortSeek (spp, start, SEEK_SET);
  }

  pos = start;

  count = 0;
  cnt = 0;
  blk = 0;
  lin = 0;

  is_na = ISA_na (bsp->mol);

  ctr = (Int2) MIN ((Int4) (stop - pos), (Int4) sizeof (bases));
  ctr = SeqPortRead (spp, bases, ctr);
  i = 0;
  residue = (Uint1) bases [i];
  while (pos < stop && residue != SEQPORT_EOF) {

    if (residue == INVALID_RESIDUE) {
      if (is_na) {
        residue = 'N';
      } else {
        residue = 'X';
      }
    }

    if (IS_residue (residue)) {

      buf [count] = (Char) (TO_LOWER (residue));
      count++;
      cnt++;
      pos++;

      blk++;
      lin++;
      if (lin >= 60) {

        buf [count] = '\0';
        PrintSeqLine (afp->format, buf, start, start + cnt);
        count = 0;
        cnt = 0;
        blk = 0;
        lin = 0;
        start += 60;

      } else if (blk >= 10) {

        buf [count] = ' ';
        count++;
        blk = 0;

      }
    }

    i++;
    if (i >= ctr) {
      i = 0;
      ctr = (Int2) MIN ((Int4) (stop - pos), (Int4) sizeof (bases));
      ctr = SeqPortRead (spp, bases, ctr);
      if (ctr < 0) {
        bases [0] = -ctr;
      } else if (ctr < 1) {
        bases [0] = SEQPORT_EOF;
      }
    }
    residue = (Uint1) bases [i];
  }

  buf [count] = '\0';
  if (count > 0) {
    PrintSeqLine (afp->format, buf, start, start + cnt);
  }

  if (! ajp->keepSeqPortOpen) {
    asp->spp = SeqPortFree (asp->spp);
  }

  return gb_MergeString (FALSE);
}

/* !!! not yet implemented !!! */

static void PrintGenome (
  SeqLocPtr slp_head,
  CharPtr prefix
)

{
  gb_AddString (prefix, "Not yet implemented", NULL, FALSE, FALSE, FALSE);
}

static CharPtr FormatContigBlock (
  Asn2gbFormatPtr afp,
  BaseBlockPtr bbp
)

{
  Asn2gbSectionPtr  asp;
  BioseqPtr         bsp;
  DeltaSeqPtr       dsp;
  SeqLitPtr         litp;
  CharPtr           prefix = NULL;
  SeqLocPtr         slp_head = NULL;
  Char              val [20];

  if (afp == NULL || bbp == NULL) return NULL;

  asp = afp->asp;
  if (asp == NULL) return NULL;
  bsp = (asp->bsp);
  if (bsp == NULL) return NULL;


  ff_StartPrint (0, 12, ASN2FF_GB_MAX, NULL);
  ff_AddString ("CONTIG");
  TabToColumn (13);
  ff_AddString("join(");

  if (bsp->seq_ext_type == 1) {

    slp_head = (SeqLocPtr) bsp->seq_ext;
    PrintGenome (slp_head, prefix);

  } else if (bsp->seq_ext_type == 4) {

    for (dsp = (DeltaSeqPtr) bsp->seq_ext; dsp; dsp=dsp->next) {
      if (dsp->choice == 1) {

        slp_head = (SeqLocPtr) dsp->data.ptrvalue;
        PrintGenome (slp_head, prefix);

      } else {

        litp = (SeqLitPtr) dsp->data.ptrvalue;
        if (litp != NULL) {
          if (litp->seq_data != NULL) {
            if (litp->length == 0) {
              sprintf (val, "gap(%ld)", (long) litp->length);
              ff_AddString (val);
            } else {
              /* don't know what to do here */
            }
          } else {
            sprintf (val, ",gap(%ld)", (long) litp->length);
            ff_AddString (val);
          }
        }
      }

      prefix = ",";
    }
  }

  ff_AddChar (')');
  ff_EndPrint ();

  return gb_MergeString (FALSE);
}


/* ********************************************************************** */

/* functions to record sections or blocks in linked lists */

static Asn2gbSectionPtr Asn2gbAddSection (
  Asn2gbWorkPtr awp
)

{
  Asn2gbSectionPtr  asp;
  ValNodePtr        vnp;

  if (awp == NULL) return NULL;

  asp = (Asn2gbSectionPtr) MemNew (sizeof (Asn2gbSection));
  if (asp == NULL) return NULL;

  vnp = ValNodeAddPointer (&(awp->lastsection), 0, asp);
  if (vnp == NULL) return asp;

  awp->lastsection = vnp;
  if (awp->sectionList == NULL) {
    awp->sectionList = vnp;
  }

  return asp;
}

static BaseBlockPtr Asn2gbAddBlock (
  Asn2gbWorkPtr awp,
  BlockType blocktype,
  size_t size
)

{
  BaseBlockPtr  bbp;
  ValNodePtr    vnp;

  if (awp == NULL || size < 1) return NULL;

  bbp = (BaseBlockPtr) MemNew (size);
  if (bbp == NULL) return NULL;
  bbp->blocktype = blocktype;
  bbp->section = awp->currsection;

  vnp = ValNodeAddPointer (&(awp->lastblock), 0, bbp);
  if (vnp == NULL) return bbp;

  awp->lastblock = vnp;
  if (awp->blockList == NULL) {
    awp->blockList = vnp;
  }

  return bbp;
}


/* ********************************************************************** */

/* add functions allocate specific blocks, populate with paragraph print info */

static CharPtr strd [4] = {
  "   ", "ss-", "ds-", "ms-"
};

static CharPtr gnbk_mol [9] = {
  "    ", "DNA ", "RNA ", "mRNA", "rRNA", "tRNA", "uRNA", "scRNA", " AA "
};

/* this looks wrong */

static CharPtr embl_mol [8] = {
  "xxx", "DNA", "RNA", "RNA", "RNA", "RNA", "RNA", "AA "
};

static CharPtr embl_divs [17] = {
  "FUN", "INV", "MAM", "ORG", "PHG", "PLN", "PRI", "PRO", "ROD"
  "SYN", "UNA", "VRL", "VRT", "PAT", "EST", "STS", "HUM"
};

static void AddLocusBlock (
  Asn2gbWorkPtr awp
)

{
  BaseBlockPtr       bbp;
  DatePtr            best_create_date = NULL;
  DatePtr            best_update_date = NULL;
  BioSourcePtr       biop;
  BioseqPtr          bsp;
  Char               date [40];
  SeqMgrDescContext  dcontext;
  Char               div [10];
  DatePtr            dp;
  EMBLBlockPtr       ebp;
  SeqMgrFeatContext  fcontext;
  GBBlockPtr         gbp;
  Int2               imol;
  Int2               istrand;
  Char               len [15];
  Int4               length;
  Char               locus [41];
  MolInfoPtr         mip;
  Char               mol [30];
  OrgNamePtr         onp;
  Uint1              origin;
  OrgRefPtr          orp;
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp;
  SeqIdPtr           sip;
  Int2               status;
  Uint1              tech;
  Uint1              topology;

  if (awp == NULL) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;

  bbp = Asn2gbAddBlock (awp, LOCUS_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  mol [0] = '\0';
  len [0] = '\0';
  div [0] = '\0';
  date [0] = '\0';

  /* locus id */

  sip = NULL;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_GENBANK ||
        sip->choice == SEQID_EMBL ||
        sip->choice == SEQID_DDBJ) break;
  }
  if (sip == NULL) {
    sip = SeqIdFindBest (bsp->id, SEQID_GENBANK);
  }
  SeqIdWrite (sip, locus, PRINTID_TEXTID_LOCUS, sizeof (locus) - 1);

  /* more complicated code to get parent locus, if segmented, goes here */

  if (awp->slp != NULL) {
    length = SeqLocLen (awp->slp);
  } else {
    length = bsp->length;
  }

  mip = NULL;
  tech = MI_TECH_standard;
  origin = 0;
  imol = bsp->mol;
  if (imol > Seq_mol_aa) {
    imol = 0;
  }
  istrand = bsp->strand;
  if (istrand > Seq_strand_both) {
    istrand = Seq_strand_unknown;
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
  if (sdp != NULL) {
    bbp->entityID = dcontext.entityID;
    bbp->itemID = dcontext.itemID;
    bbp->itemtype = OBJ_SEQDESC;

    mip = (MolInfoPtr) sdp->data.ptrvalue;
    if (mip != NULL) {
      if (mip->biomol <= MOLECULE_TYPE_PEPTIDE) {
        imol = (Int2) mip->biomol;
      }
      tech = mip->tech;
    }
  }

  /* check inst.mol if mol-type is not-set or genomic */

  if (imol < Seq_mol_rna) {
    imol = bsp->mol;
    if (imol == Seq_mol_aa) {
      imol = MOLECULE_TYPE_PEPTIDE;
    } else if (imol == Seq_mol_na) {
      imol = 0;
    }
  }

  /* if ds-DNA don't show ds */

  if (imol == Seq_mol_dna && istrand == Seq_strand_minus) {
    istrand = Seq_strand_unknown;
  }

  /* ss=any RNA don't show ss */

  if (imol > Seq_mol_rna && istrand == Seq_strand_plus) {
    istrand = Seq_strand_unknown;
  }

  topology = bsp->topology;
  if (awp->slp != NULL) {
    topology = TOPOLOGY_LINEAR;
  }

  /* length, topology, and molecule type */

  if (awp->format == GENBANK_FMT) {

    if (topology == TOPOLOGY_CIRCULAR) {
      sprintf (len, "%7ld bp", (long) length);
      sprintf (mol, "%s%-4s  circular", strd [istrand], gnbk_mol [imol]);
    } else {
      sprintf (len, "%7ld bp", (long) length);
      sprintf (mol, "%s%-4s          ", strd [istrand], gnbk_mol [imol]);
    }

  } else if (awp->format == GENPEPT_FMT) {

    sprintf (len, "%7ld aa", (long) length);

  } else if (awp->format == EMBL_FMT) {

    if (imol < MOLECULE_TYPE_PEPTIDE) {
      if (topology == TOPOLOGY_CIRCULAR) {
        sprintf (mol, "circular %s", embl_mol [imol]);
        sprintf (len, "%7ld BP.", (long) length);
      } else {
        sprintf (mol, "%s", embl_mol [imol]);
        sprintf (len, "%7ld BP.", (long) length);
      }
    }
  }

  /* division */

  biop = NULL;
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp != NULL) {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
  } else {
    sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_BIOSRC, 0, &fcontext);
    if (sfp != NULL) {
      biop = (BioSourcePtr) sfp->data.value.ptrvalue;
    }
  }
  if (biop != NULL) {
    origin = biop->origin;
    orp = biop->org;
    if (orp != NULL) {
      onp = orp->orgname;
      if (onp != NULL) {
        StringNCpy_0 (div, onp->div, sizeof (div));
      }
    }
  }

  switch (tech) {
    case MI_TECH_est :
      StringCpy (div, "EST");
      break;
    case MI_TECH_sts :
      StringCpy (div, "STS");
      break;
    case MI_TECH_survey :
      StringCpy (div, "GSS");
      break;
    case MI_TECH_htgs_0 :
    case MI_TECH_htgs_1 :
    case MI_TECH_htgs_2 :
      StringCpy (div, "HTG");
      break;
    default :
      break;
  }

  if (origin == 5) {
    StringCpy (div, "SYN");
  }

  sip = SeqIdFindBest (bsp->id, SEQID_PATENT);
  if (sip != NULL && sip->choice == SEQID_PATENT) {
    StringCpy (div, "PAT");
  }

  /* more complicated code for division, if necessary, goes here */

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_genbank, &dcontext);
  while (sdp != NULL) {
    gbp = (GBBlockPtr) sdp->data.ptrvalue;
    if (gbp != NULL) {
      if (StringHasNoText (div) && gbp->div != NULL) {
        StringCpy (div, gbp->div);
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_genbank, &dcontext);
  }

  if (awp->format == EMBL_FMT || awp->format == EMBLPEPT_FMT) {

    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_embl, &dcontext);
    if (sdp != NULL) {
      ebp = (EMBLBlockPtr) sdp->data.ptrvalue;
      if (ebp != NULL) {
        if (ebp->div == 255) {
          if (mip == NULL) {
            StringCpy (div, "HUM");
          }
        } else if (ebp->div < 16)  {
          StringCpy (div, embl_divs [ebp->div]);
        }
      }
    }

    if (StringHasNoText (div)) {
      StringCpy (div, "UNA");
    }
  }

  /* empty division field if unable to find anything */

  if (StringHasNoText (div)) {
    StringCpy (div, "   ");
  }

  /* date */

  dp = NULL;
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_update_date, &dcontext);
  while (sdp != NULL) {
    dp = (DatePtr) sdp->data.ptrvalue;
    if (dp != NULL) {
      if (best_update_date == NULL) {
        best_update_date = dp;
      } else {
        status = DateMatch (dp, best_update_date, FALSE);
        if (status == 1) {
          best_update_date = dp;
        }
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_update_date, &dcontext);
  }

  /* !!! temporarily also look at genbank block entry date !!! */

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_genbank, &dcontext);
  while (sdp != NULL) {
    gbp = (GBBlockPtr) sdp->data.ptrvalue;
    if (gbp != NULL) {
      dp = gbp->entry_date;
      if (dp != NULL) {
        if (best_update_date == NULL) {
          best_update_date = dp;
        } else {
          status = DateMatch (dp, best_update_date, FALSE);
          if (status == 1) {
            best_update_date = dp;
          }
        }
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_genbank, &dcontext);
  }

  dp = NULL;
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_create_date, &dcontext);
  while (sdp != NULL) {
    dp = (DatePtr) sdp->data.ptrvalue;
    if (dp != NULL) {
      if (best_create_date == NULL) {
        best_create_date = dp;
      } else {
        status = DateMatch (dp, best_create_date, FALSE);
        if (status == 1) {
          best_create_date = dp;
        }
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_create_date, &dcontext);
  }

  dp = NULL;
  if (best_update_date != NULL && best_create_date != NULL) {
    status = DateMatch (best_update_date, best_create_date, FALSE);
    if (status == 0 || status == 1) {
      dp = best_update_date;
    } else {
      dp = best_create_date;
    }
  } else if (best_update_date != NULL) {
    dp = best_update_date;
  } else if (best_create_date != NULL) {
    dp = best_create_date;
  }
  if (dp != NULL) {
    DateToGB (date, dp);
  }
  if (StringHasNoText (date)) {
    StringCpy (date, "01-JAN-1900");
  }

  /* more complicated code for dates from various objects goes here */

  init_buff ();

  if (awp->format == GENBANK_FMT || awp->format == GENPEPT_FMT) {

    ff_StartPrint (0, 0, ASN2FF_GB_MAX, NULL);
    ff_AddString ("LOCUS");
    TabToColumn (13);
    ff_AddString (locus);
    TabToColumn (33 - StringLen (len));
    ff_AddString (len);
    TabToColumn (34);
    ff_AddString (mol);
    TabToColumn (53);
    ff_AddString (div);
    TabToColumn (63);
    ff_AddString (date);

  } else if (awp->format == EMBL_FMT || awp->format == EMBLPEPT_FMT) {

    ff_StartPrint (5, 0, ASN2FF_EMBL_MAX, "ID");
    ff_AddString (locus);
    if (awp->hup) {
      ff_AddString (" confidential; ");
    } else {
      ff_AddString (" standard; ");
    }
    ff_AddString (mol);
    ff_AddString ("; ");

    /* conditional code to make div "UNA" goes here */
  
    ff_AddString (div);
    ff_AddString ("; ");
    ff_AddString (len);
  }

  bbp->string = gb_MergeString (TRUE);
}

static void AddDeflineBlock (
  Asn2gbWorkPtr awp
)

{
  BaseBlockPtr       bbp;
  BioseqPtr          bsp;
  CharPtr            buf;
  size_t             buflen = 1001;
  SeqMgrDescContext  dcontext;
  ItemInfo           ii;
  MolInfoPtr         mip;
  SeqDescrPtr        sdp;
  Uint1              tech;

  if (awp == NULL) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;

  bbp = Asn2gbAddBlock (awp, DEFLINE_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  tech = 0;
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
  if (sdp != NULL) {
    mip = (MolInfoPtr) sdp->data.ptrvalue;
    if (mip != NULL) {
      tech = mip->tech;
    }
  }

  buf = MemNew (sizeof (Char) * (buflen + 1));
  MemSet ((Pointer) (&ii), 0, sizeof (ItemInfo));

  /* create default defline */

  if (buf != NULL && CreateDefLine (&ii, bsp, buf, buflen, tech, NULL, NULL)) {
    bbp->entityID = ii.entityID;
    bbp->itemID = ii.itemID;
    bbp->itemtype = ii.itemtype;

    gb_StartPrint (awp->format, TRUE, 0, 12, "DEFINITION", 13, 5, 5, "DE", TRUE);

    gb_AddString (NULL, buf, NULL, TRUE, TRUE, FALSE);

    bbp->string = gb_MergeString (TRUE);
  }

  MemFree (buf);
}

/*
Return values are:
 0: no problem - Accession is in proper format
-1: Accession did not start with a letter (or two letters)
-2: Accession did not contain five numbers (or six numbers after 2 letters)
-3: the original Accession number to be validated was NULL
*/

static Int2 ValidateAccession (
  CharPtr accession
)

{
  Int2  count;

  if (accession == NULL || accession [0] == '\0') return -3;

  if (accession [0] < 'A' || accession [0] > 'Z') return -1;

  for (count=1; count < 5; count++) {
    if (! IS_DIGIT (accession [count])) break;
  }

  if (count == 5 && (accession [count + 1] == '\0' || accession [count + 1] == ' ')) {
    return 0;
  }

  if (IS_ALPHA (accession [1])) {
    if (accession [1] < 'A' || accession [1] > 'Z') return -1;
    for (count = 2; count < 7; count++) {
      if (! IS_DIGIT (accession [count])) break;
    }
    if (count == 7 && (accession [count + 1] == '\0' || accession [count + 1] == ' ')) return 0;
  }
  return -2;
}

/* this definitely needs more work to support all classes, use proper SeqId */

static void AddAccessionBlock (
  Asn2gbWorkPtr awp
)

{
  SeqIdPtr           accn = NULL;
  BaseBlockPtr       bbp;
  BioseqPtr          bsp;
  Char               buf [41];
  SeqMgrDescContext  dcontext;
  EMBLBlockPtr       ebp;
  ValNodePtr         extra_access;
  GBBlockPtr         gbp;
  SeqIdPtr           gi = NULL;
  ValNodePtr         head = NULL;
  SeqIdPtr           lcl = NULL;
  SeqDescrPtr        sdp;
  CharPtr            separator;
  SeqIdPtr           sip;
  CharPtr            str;
  ValNodePtr         vnp;
  CharPtr            xtra;

  if (awp == NULL) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    switch (sip->choice) {
      case SEQID_GI :
        gi = sip;
        break;
      case SEQID_GENBANK :
      case SEQID_EMBL :
      case SEQID_DDBJ :
      case SEQID_OTHER :
        accn = sip;
        break;
      case SEQID_LOCAL :
        lcl = sip;
        break;
      default :
        break;
    }
  }

  sip = NULL;
  if (accn != NULL) {
    sip = accn;
  } else if (lcl != NULL) {
    sip = lcl;
  } else if (gi != NULL) {
    sip = gi;
  }

  if (sip == NULL) return;

  SeqIdWrite (sip, buf, PRINTID_TEXTID_ACC_ONLY, sizeof (buf));

  bbp = Asn2gbAddBlock (awp, ACCESSION_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  gb_StartPrint (awp->format, TRUE, 0, 12, "ACCESSION", 13, 5, 5, "AC", TRUE);

  if (awp->hup && accn != NULL) {
    gb_AddString (NULL, ";", NULL, FALSE, FALSE, FALSE);
  } else {
    gb_AddString (NULL, buf, NULL, FALSE, FALSE, FALSE);
  }

  if (awp->format == GENBANK_FMT || awp->format == GENPEPT_FMT) {
    separator = " ";
  } else if (awp->format == EMBL_FMT || awp->format == EMBLPEPT_FMT) {
    separator = ";";
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, 0, &dcontext);
  while (sdp != NULL) {

    extra_access = NULL;

    switch (dcontext.seqdesctype) {
      case Seq_descr_genbank :
        gbp = (GBBlockPtr) sdp->data.ptrvalue;
        if (gbp != NULL) {
          extra_access = gbp->extra_accessions;
        }
        break;
      case Seq_descr_embl :
        ebp = (EMBLBlockPtr) sdp->data.ptrvalue;
        if (ebp != NULL) {
          extra_access = ebp->extra_acc;
        }
        break;
      default :
        break;
    }

    if (extra_access != NULL) {
      bbp->entityID = dcontext.entityID;
      bbp->itemID = dcontext.itemID;
      bbp->itemtype = OBJ_SEQDESC;
    }

    for (vnp = extra_access; vnp != NULL; vnp = vnp->next) {
      xtra = (CharPtr) vnp->data.ptrvalue;
      if (ValidateAccession (xtra) == 0) {
        ValNodeCopyStr (&head, 0, separator);
        ValNodeCopyStr (&head, 0, xtra);
      }
    }

    sdp = SeqMgrGetNextDescriptor (bsp, sdp, 0, &dcontext);
  }

  str = MergeValNodeStrings (head);

  gb_AddString (NULL, str, NULL, FALSE, FALSE, FALSE);

  MemFree (str);
  ValNodeFreeData (head);

  bbp->string = gb_MergeString (TRUE);
}

static void AddVersionBlock (
  Asn2gbWorkPtr awp
)

{
  SeqIdPtr      accn = NULL;
  BaseBlockPtr  bbp;
  BioseqPtr     bsp;
  Char          buf [41];
  Int4          gi = -1;
  Boolean       needEndPrint = TRUE;
  Boolean       needInitBuff = TRUE;
  SeqIdPtr      sip;
  Char          version [64];

  if (awp == NULL) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    switch (sip->choice) {
      case SEQID_GI :
        gi = sip->data.intvalue;
        break;
      case SEQID_GENBANK :
      case SEQID_EMBL :
      case SEQID_DDBJ :
      case SEQID_OTHER :
        accn = sip;
        break;
      default :
        break;
    }
  }

  /* if (gi < 1 && accn == NULL) return; */

  bbp = Asn2gbAddBlock (awp, VERSION_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  /* no longer displaying NID */

  /*
  if (gi > 0) {
    sprintf (version, "g%ld", (long) gi);

    gb_StartPrint (awp->format, needInitBuff, 0, 12, "NID", 13, 5, 5, "NI", TRUE);
    needInitBuff = FALSE;

    gb_AddString (NULL, version, NULL, FALSE, FALSE, FALSE);

	ff_EndPrint();
	needEndPrint = FALSE;
  }
  */

  if (SHOWVERSION && accn != NULL) {

    buf [0] = '\0';
    SeqIdWrite (accn, buf, PRINTID_TEXTID_ACC_VER, sizeof (buf) - 1);

    if (awp->format == GENBANK_FMT || awp->format == GENPEPT_FMT) {
      if (gi > 0) {
        sprintf (version, "%s  GI:%ld", buf, (long) gi);
      } else {
        sprintf (version, "%s", buf);
      }
    } else if (awp->format == EMBL_FMT || awp->format == EMBLPEPT_FMT) {
      if (gi > 0) {
        sprintf (version, "%s.%ld", buf, (long) gi);
      } else {
        sprintf (version, "%s", buf);
      }
    }

    gb_StartPrint (awp->format, needInitBuff, 0, 12, "VERSION", 13, 5, 5, "SV", TRUE);
    needInitBuff = FALSE;

    gb_AddString (NULL, version, NULL, FALSE, FALSE, FALSE);

	ff_EndPrint();
	needEndPrint = FALSE;

  } else if (SHOWVERSION) {

    gb_StartPrint (awp->format, needInitBuff, 0, 0, "VERSION", 0, 5, 5, "SV", TRUE);
    needInitBuff = FALSE;
	ff_EndPrint();
	needEndPrint = FALSE;
  }

  bbp->string = gb_MergeString (needEndPrint);
}

/* no longer displaying PID */

/*
static void AddPidBlock (Asn2gbWorkPtr awp)

{
  BaseBlockPtr  bbp;

  if (awp == NULL) return;

  bbp = Asn2gbAddBlock (awp, PID_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  bbp->string = StringSave ("PID\n");
}
*/

static void AddDbsourceBlock (
  Asn2gbWorkPtr awp
)

{
  BaseBlockPtr  bbp;

  if (awp == NULL) return;

  bbp = Asn2gbAddBlock (awp, DBSOURCE_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  bbp->string = StringSave ("DBSOURCE\n");
}

static void AddDateBlock (
  Asn2gbWorkPtr awp
)

{
  BaseBlockPtr       bbp;
  BioseqPtr          bsp;
  Char               date [40];
  SeqMgrDescContext  dcontext;
  DatePtr            dp;
  SeqDescrPtr        sdp;

  if (awp == NULL) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;

  bbp = Asn2gbAddBlock (awp, DATE_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  date [0] = '\0';

  dp = NULL;
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_update_date, &dcontext);
  if (sdp == NULL) {
    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_create_date, &dcontext);
  }
  if (sdp != NULL) {
    dp = (DatePtr) sdp->data.ptrvalue;
  }
  if (dp != NULL) {
    DateToGB (date, dp);
  }
  if (StringHasNoText (date)) {
    StringCpy (date, "01-JAN-1900");
  }

  init_buff ();
  PrintXX ();
  ff_StartPrint (5, 5, ASN2FF_EMBL_MAX, "DT");
  ff_AddString (date);

  bbp->string = gb_MergeString (TRUE);
}

#define TOTAL_ESTKW 11
#define TOTAL_STSKW 5
#define TOTAL_GSSKW 2

static CharPtr EST_kw_array[ TOTAL_ESTKW] = {
  "EST", "EST PROTO((expressed sequence tag)", "expressed sequence tag",
  "EST (expressed sequence tag)", "EST(expressed sequence tag)",
  "partial cDNA sequence", "transcribed sequence fragment", "TSR",
  "putatively transcribed partial sequence", "UK putts"
};

static CharPtr GSS_kw_array [TOTAL_GSSKW] = {
  "GSS", "trapped exon"
};
static CharPtr STS_kw_array[TOTAL_STSKW] = {
  "STS", "STS(sequence tagged site)", "STS (sequence tagged site)", 
  "STS sequence", "sequence tagged site"
};

static Int2 MatchArrayString (
  CharPtr array_string [],
  Int2 totalstr,
  CharPtr text
)

{
  Int2 i;

  for (i = 0; i < totalstr && text != NULL; i++) {
    if (StringCmp (array_string [i], text) == 0) {
      return (i);
    }
  }

  return (-1);
}

static Boolean CheckSpecialKeyword (
  Boolean is_est,
  Boolean is_sts,
  Boolean is_gss,
  CharPtr kwd
)

{
  if (kwd == NULL) return FALSE;

  if (is_est) {
    if (MatchArrayString (STS_kw_array, TOTAL_STSKW, kwd) != -1) return FALSE;
    if (MatchArrayString (GSS_kw_array, TOTAL_GSSKW, kwd) != -1) return FALSE;
  }

  if (is_sts) {
    if (MatchArrayString (EST_kw_array, TOTAL_ESTKW, kwd) != -1) return FALSE;
    if (MatchArrayString (GSS_kw_array, TOTAL_GSSKW, kwd) != -1) return FALSE;
  }

  if (is_gss) {
    if (MatchArrayString (STS_kw_array, TOTAL_STSKW, kwd) != -1) return FALSE;
    if (MatchArrayString (EST_kw_array, TOTAL_ESTKW, kwd) != -1) return FALSE;
  }

  return TRUE;
}

static Boolean KeywordAlreadyInList (
  ValNodePtr head,
  CharPtr kwd
)

{
  ValNodePtr  vnp;

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    if (StringICmp ((CharPtr) vnp->data.ptrvalue, kwd) == 0) return TRUE;
  }

  return FALSE;
}

static void AddKeywordsBlock (
  Asn2gbWorkPtr awp
)

{
  BaseBlockPtr       bbp;
  BioseqPtr          bsp;
  SeqMgrDescContext  dcontext;
  EMBLBlockPtr       ebp;
  GBBlockPtr         gbp;
  ValNodePtr         head = NULL;
  Boolean            is_est = FALSE;
  Boolean            is_gss = FALSE;
  Boolean            is_sts = FALSE;
  ValNodePtr         keywords;
  CharPtr            kwd;
  MolInfoPtr         mip;
  PirBlockPtr        pir;
  PrfBlockPtr        prf;
  SeqDescrPtr        sdp;
  SPBlockPtr         sp;
  CharPtr            str;
  ValNodePtr         vnp;

  if (awp == NULL) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;

  bbp = (BaseBlockPtr) Asn2gbAddBlock (awp, KEYWORDS_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
  if (sdp != NULL) {
    bbp->entityID = dcontext.entityID;
    bbp->itemID = dcontext.itemID;
    bbp->itemtype = OBJ_SEQDESC;

    mip = (MolInfoPtr) sdp->data.ptrvalue;
    if (mip != NULL) {
      switch (mip->tech) {
        case MI_TECH_htgs_1 :
          ValNodeCopyStr (&head, 0, "HTG; HTGS_PHASE1");
          break;
        case MI_TECH_htgs_2 :
          ValNodeCopyStr (&head, 0, "HTG; HTGS_PHASE2");
          break;
        case MI_TECH_htgs_3 :
          ValNodeCopyStr (&head, 0, "HTG");
          break;
        case MI_TECH_est :
          is_est = TRUE;
          ValNodeCopyStr (&head, 0, "EST");
          break;
        case MI_TECH_sts :
          is_sts = TRUE;
          ValNodeCopyStr (&head, 0, "STS");
          break;
        case MI_TECH_survey :
          is_gss = TRUE;
          ValNodeCopyStr (&head, 0, "GSS");
          break;
        case MI_TECH_fli_cdna :
          ValNodeCopyStr (&head, 0, "FLI_CDNA");
          break;
        case MI_TECH_htgs_0 :
          ValNodeCopyStr (&head, 0, "HTG; HTGS_PHASE0");
          break;
        default :
          break;
      }
    }
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, 0, &dcontext);
  while (sdp != NULL) {

    keywords = NULL;

    switch (dcontext.seqdesctype) {
      case Seq_descr_genbank :
        gbp = (GBBlockPtr) sdp->data.ptrvalue;
        if (gbp != NULL) {
          keywords = gbp->keywords;
        }
        break;
      case Seq_descr_embl :
        ebp = (EMBLBlockPtr) sdp->data.ptrvalue;
        if (ebp != NULL) {
          keywords = ebp->keywords;
        }
        break;
      case Seq_descr_pir :
        pir = (PirBlockPtr) sdp->data.ptrvalue;
        if (pir != NULL) {
          keywords = pir->keywords;
        }
        break;
      case Seq_descr_prf :
        prf = (PrfBlockPtr) sdp->data.ptrvalue;
        if (prf != NULL) {
          keywords = prf->keywords;
        }
        break;
      case Seq_descr_sp :
        sp = (SPBlockPtr) sdp->data.ptrvalue;
        if (sp != NULL) {
          keywords = sp->keywords;
        }
        break;
      default :
        break;
    }

    if (keywords != NULL) {
      bbp->entityID = dcontext.entityID;
      bbp->itemID = dcontext.itemID;
      bbp->itemtype = OBJ_SEQDESC;
    }

    for (vnp = keywords; vnp != NULL; vnp = vnp->next) {
      kwd = (CharPtr) vnp->data.ptrvalue;
      if (CheckSpecialKeyword (is_est, is_sts, is_gss, kwd)) {
        if (! KeywordAlreadyInList (head, kwd)) {
          if (head != NULL) {
            ValNodeCopyStr (&head, 0, "; ");
          }
          ValNodeCopyStr (&head, 0, kwd);
        }
      }
    }

    sdp = SeqMgrGetNextDescriptor (bsp, sdp, 0, &dcontext);
  }

  gb_StartPrint (awp->format, TRUE, 0, 12, "KEYWORDS", 13, 5, 5, "KW", TRUE);

  str = MergeValNodeStrings (head);

  /* if no keywords were found, period will still be added by this call */

  gb_AddString (NULL, str, NULL, TRUE, FALSE, FALSE);

  MemFree (str);
  ValNodeFreeData (head);

  bbp->string = gb_MergeString (TRUE);
}

static void AddSegmentBlock (
  Asn2gbWorkPtr awp
)

{
  BaseBlockPtr  bbp;
  Char          buf [32];

  if (awp == NULL) return;

  if (awp->seg < 1 || awp->numsegs < 1) return;

  bbp = Asn2gbAddBlock (awp, SEGMENT_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  sprintf (buf, "%d of %ld", (int) awp->seg, (long) awp->numsegs);

  gb_StartPrint (awp->format, TRUE, 0, 12, "SEGMENT", 13, 5, 5, "XX", FALSE);

  gb_AddString (NULL, buf, NULL, FALSE, TRUE, FALSE);

  bbp->string = gb_MergeString (TRUE);
}

static void AddSourceBlock (
  Asn2gbWorkPtr awp
)

{
  BaseBlockPtr       bbp;
  BioseqPtr          bsp;
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext;
  GBBlockPtr         gbp;
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp;

  if (awp == NULL) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;

  bbp = Asn2gbAddBlock (awp, SOURCE_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_genbank, &dcontext);
  if (sdp != NULL) {
    gbp = (GBBlockPtr) sdp->data.ptrvalue;
    if (gbp != NULL && (! StringHasNoText (gbp->source))) {
      bbp->entityID = dcontext.entityID;
      bbp->itemID = dcontext.itemID;
      bbp->itemtype = OBJ_SEQDESC;
      return;
    }
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp != NULL) {
    bbp->entityID = dcontext.entityID;
    bbp->itemID = dcontext.itemID;
    bbp->itemtype = OBJ_SEQDESC;
  } else {
    sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_BIOSRC, 0, &fcontext);
    if (sfp != NULL) {
      bbp->entityID = fcontext.entityID;
      bbp->itemID = fcontext.itemID;
      bbp->itemtype = OBJ_SEQFEAT;
    }
  }
}

static void AddOrganismBlock (
  Asn2gbWorkPtr awp
)

{
  BaseBlockPtr       bbp;
  BioseqPtr          bsp;
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext;
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp;

  if (awp == NULL) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;

  bbp = Asn2gbAddBlock (awp, ORGANISM_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp != NULL) {
    bbp->entityID = dcontext.entityID;
    bbp->itemID = dcontext.itemID;
    bbp->itemtype = OBJ_SEQDESC;
  } else {
    sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_BIOSRC, 0, &fcontext);
    if (sfp != NULL) {
      bbp->entityID = fcontext.entityID;
      bbp->itemID = fcontext.itemID;
      bbp->itemtype = OBJ_SEQFEAT;
    }
  }
}

static ReferenceBlockPtr AddPub (
  Asn2gbWorkPtr awp,
  ValNodePtr PNTR head,
  PubdescPtr pdp
)

{
  Char               buf [121];
  CitArtPtr          cap;
  CitBookPtr         cbp;
  CitGenPtr          cgp;
  CitJourPtr         cjp;
  CitPatPtr          cpp;
  CitSubPtr          csp;
  DatePtr            dp = NULL;
  Boolean            justuids = TRUE;
  ImprintPtr         imp = NULL;
  IntRefBlockPtr     irp;
  ReferenceBlockPtr  rbp;
  ValNodePtr         vnp;

  if (awp == NULL || head == NULL || pdp == NULL) return NULL;

  rbp = (ReferenceBlockPtr) MemNew (sizeof (IntRefBlock));
  if (rbp == NULL) return NULL;
  rbp->blocktype = REFERENCE_BLOCK;
  rbp->section = awp->currsection;

  rbp->serial = INT2_MAX;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    switch (vnp->choice) {
      case PUB_Gen :
        /* may be unpublished, or may be serial number of swiss-prot reference */
        cgp = (CitGenPtr) vnp->data.ptrvalue;
        if (cgp != NULL) {
          if (StringNICmp ("BackBone id_pub", cgp->cit, 15) != 0) {
            rbp->category = REF_CAT_UNP;
            dp = cgp->date;
            if (cgp->serial_number > 0) {
              rbp->serial = cgp->serial_number;
            }
            if (cgp->cit != NULL) {
              if (StringNICmp ("unpublished", cgp->cit, 11) != 0 &&
                  StringNICmp ("submitted", cgp->cit, 8) != 0 &&
                  StringNICmp ("to be published", cgp->cit, 15) != 0 &&
                  StringNICmp ("in press", cgp->cit, 8) != 0 &&
                  StringStr (cgp->cit, "Journal") == NULL) {
                MemFree (rbp);
                return NULL;
              }
            } else if (cgp->journal == NULL || cgp->date == NULL) {
              MemFree (rbp);
              return NULL;
            }
          }
        }
        break;
      case PUB_Sub :
        rbp->category = REF_CAT_SUB;
        csp = (CitSubPtr) vnp->data.ptrvalue;
        if (csp != NULL) {
          imp = csp->imp;
          if (imp != NULL) {
            dp = imp->date;
          }
          if (csp->date != NULL) {
            dp = csp->date;
          }
        }
        break;
      case PUB_Article:
        cap = (CitArtPtr) vnp->data.ptrvalue;
        if (cap != NULL) {
          switch (cap->from) {
            case 1:
              cjp = (CitJourPtr) cap->fromptr;
              if (cjp != NULL) {
                imp = (ImprintPtr) cjp->imp;
                if (imp != NULL) {
                  dp = imp->date;
                }
              }
              break;
            case 2:
              cbp = (CitBookPtr) cap->fromptr;
              if (cbp != NULL) {
                imp = (ImprintPtr) cbp->imp;
                if (imp != NULL) {
                  dp = imp->date;
                }
              }
              break;
            case 3:
              cbp = (CitBookPtr) cap->fromptr;
              if (cbp != NULL) {
                imp = (ImprintPtr) cbp->imp;
                if (imp != NULL) {
                  dp = imp->date;
                }
              }
              break;
            default:
              break;
          }
        }
        break;
      case PUB_Book:
        cbp = (CitBookPtr) vnp->data.ptrvalue;
        if (cbp != NULL) {
          imp = (ImprintPtr) cbp->imp;
          if (imp != NULL) {
            dp = imp->date;
          }
        }
        break;
      case PUB_Proc:
        cbp = (CitBookPtr) vnp->data.ptrvalue;
        if (cbp != NULL) {
          imp = (ImprintPtr) cbp->imp;
          if (imp != NULL) {
            dp = imp->date;
          }
        }
        break;
      case PUB_Patent :
        rbp->category = REF_CAT_PUB;
        cpp = (CitPatPtr) vnp->data.ptrvalue;
        if (cpp != NULL) {
          if (cpp->date_issue != NULL) {
            dp = (DatePtr) cpp->date_issue;
          } else if (cpp->app_date != NULL) {
            dp = (DatePtr) cpp->app_date;
          }
        }
        break;
      case PUB_Man:
        cbp = (CitBookPtr) vnp->data.ptrvalue;
        if (cbp != NULL) {
          imp = (ImprintPtr) cbp->imp;
          if (imp != NULL) {
            dp = imp->date;
          }
        }
        break;
      case PUB_Muid :
        rbp->muid = vnp->data.intvalue;
        rbp->category = REF_CAT_PUB;
        break;
      case PUB_PMid :
        rbp->pmid = vnp->data.intvalue;
        rbp->category = REF_CAT_PUB;
        break;
      default :
        break;
    }
    if (vnp->choice != PUB_Muid && vnp->choice != PUB_PMid) {
      justuids = FALSE;
    }
  }

  /* check for submitted vs. in-press */

  if (imp != NULL) {
    rbp->category = REF_CAT_PUB;
    switch (imp->prepub) {
      case 1 :
        rbp->category = REF_CAT_UNP;
        break;
      case 2 :
        rbp->category = REF_CAT_PUB;
        break;
      default :
        break;
    }
  }

  /* check for sites reftype */

  if (pdp->reftype == 1) {
    rbp->sites = TRUE;
  }

  if (rbp->muid == 0 && rbp->pmid == 0) {
	if (PubLabelUnique (pdp->pub, buf, sizeof (buf) - 1, OM_LABEL_CONTENT, TRUE) > 0) {
	  rbp->uniquestr = StringSaveNoNull (buf);
	}
  }

  irp = (IntRefBlockPtr) rbp;
  irp->date = DateDup (dp);
  irp->justuids = justuids;
  if (justuids) {
    irp->fig = StringSaveNoNull (pdp->fig);
    irp->maploc = StringSaveNoNull (pdp->maploc);
    irp->poly_a = pdp->poly_a;
  }

  /* if not rejected by now, link in */

  ValNodeAddPointer (head, 0, rbp);

  return rbp;
}

static int LIBCALLBACK SortReferences (
  VoidPtr ptr1,
  VoidPtr ptr2
)

{
  int                compare;
  IntRefBlockPtr     irp1;
  IntRefBlockPtr     irp2;
  ReferenceBlockPtr  rbp1;
  ReferenceBlockPtr  rbp2;
  Int2               status;
  ValNodePtr         vnp1;
  ValNodePtr         vnp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  vnp1 = *((ValNodePtr PNTR) ptr1);
  vnp2 = *((ValNodePtr PNTR) ptr2);
  if (vnp1 == NULL || vnp2 == NULL) return 0;
  rbp1 = (ReferenceBlockPtr) vnp1->data.ptrvalue;
  rbp2 = (ReferenceBlockPtr) vnp2->data.ptrvalue;
  if (rbp1 == NULL || rbp2 == NULL) return 0;

  if (rbp1->serial > rbp2->serial) {
    return 1;
  } else if (rbp1->serial < rbp2->serial) {
    return -1;
  }

  /* usual first sort by published, unpublished, and cit-subs */

  if (rbp1->category > rbp2->category) {
    return 1;
  } else if (rbp1->category < rbp2->category) {
    return -1;
  }

  /* within class, sort by date, older publications first */

  irp1 = (IntRefBlockPtr) rbp1;
  irp2 = (IntRefBlockPtr) rbp2;
  status = DateMatch (irp1->date, irp2->date, FALSE);
  if (status == 1 || status == -1) return status;

  /* if dates (e.g., years) match, try to distinguish by uids */

  if (rbp1->pmid != 0 && rbp2->pmid != 0) {
    if (rbp1->pmid > rbp2->pmid) {
      return 1;
    } else if (rbp1->pmid < rbp2->pmid) {
      return -1;
    }
  }
  if (rbp1->muid != 0 && rbp2->muid != 0) {
    if (rbp1->muid > rbp2->muid) {
      return 1;
    } else if (rbp1->muid < rbp2->muid) {
      return -1;
    }
  }

  /* if same uid, one with just uids goes last to be excised but remembered */

  if ((rbp1->pmid != 0 && rbp2->pmid != 0) || (rbp1->muid != 0 && rbp2->muid != 0)) {
    if (irp1->justuids && (! irp2->justuids)) {
      return 1;
    } else if ((! irp1->justuids) && irp2->justuids) {
      return -1;
    }
  }

  /* put sites after pubs that refer to all or a range of bases */

  if (rbp2->sites) {
    return 1;
  } else if (rbp1->sites) {
    return -1;
  }

  /* for publication features, sort in explore index order */

  if (irp1->index > irp2->index) {
    return 1;
  } else if (irp1->index < irp2->index) {
    return -1;
  }

  /* next use author string */

  if (irp1->authstr != NULL && irp2->authstr != NULL) {
    compare = StringICmp (irp1->authstr, irp2->authstr);
    if (compare > 0) {
      return 1;
    } else if (compare < 0) {
      return -1;
    }
  }

  /* use unique label string to determine sort order */

  if (rbp1->uniquestr != NULL && rbp2->uniquestr != NULL) {
    compare = StringICmp (rbp1->uniquestr, rbp2->uniquestr);
    if (compare > 0) {
      return 1;
    } else if (compare < 0) {
      return -1;
    }
  }

  /* last resort for equivalent publication descriptors, sort in itemID order */

  if (rbp1->itemtype == OBJ_SEQDESC && rbp2->itemtype == OBJ_SEQDESC) {
    if (rbp1->itemID > rbp2->itemID) {
      return 1;
    } else if (rbp1->itemID < rbp2->itemID) {
      return -1;
    }
  }

  return 0;
}

static void GetRefsOnBioseq (
  Asn2gbWorkPtr awp,
  BioseqPtr target,
  BioseqPtr bsp,
  Int4 from,
  Int4 to
)

{
  AuthListPtr        alp;
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext;
  Int2               idx;
  IntRefBlockPtr     irp;
  Int4Ptr            ivals;
  Int2               numivals;
  Boolean            okay;
  PubdescPtr         pdp;
  ReferenceBlockPtr  rbp;
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp;
  SeqInt             sint;
  Int4               start;
  Int4               stop;
  ValNode            vn;
  ValNodePtr         vnp;

  if (awp == NULL || target == NULL || bsp == NULL) return;

  /* full length loc for descriptors */

  sint.from = 0;
  sint.to = bsp->length - 1;
  sint.strand = Seq_strand_plus;
  sint.id = SeqIdStripLocus (SeqIdDup (SeqIdFindBest (bsp->id, 0)));
  vn.choice = SEQLOC_INT;
  vn.data.ptrvalue = (Pointer) &sint;
  vn.next = NULL;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_pub, &dcontext);
  while (sdp != NULL) {

    /* check if descriptor on part already added on segmented bioseq */

    okay = TRUE;
    for (vnp = awp->pubhead; vnp != NULL && okay; vnp = vnp->next) {
      rbp = (ReferenceBlockPtr) vnp->data.ptrvalue;
      if (rbp != NULL) {
        if (rbp->entityID == dcontext.entityID &&
            rbp->itemID == dcontext.itemID &&
            rbp->itemtype == OBJ_SEQDESC) {
          okay = FALSE;
        }
      }
    }

    if (okay) {
      pdp = (PubdescPtr) sdp->data.ptrvalue;
      rbp = AddPub (awp, &(awp->pubhead), pdp);
      if (rbp != NULL) {

        rbp->entityID = dcontext.entityID;
        rbp->itemID = dcontext.itemID;
        rbp->itemtype = OBJ_SEQDESC;

        irp = (IntRefBlockPtr) rbp;
        irp->loc = SeqLocMerge (target, &vn, NULL, FALSE, TRUE, FALSE);
        alp = GetAuthListPtr (pdp, NULL);
        if (alp != NULL) {
          irp->authstr = GetAuthorsString (awp->format, alp);
        }
        irp->index = 0;
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_pub, &dcontext);
  }

  SeqIdFree (sint.id);

  /* features are indexed on parent if segmented */

  bsp = awp->parent;

  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_PUB, 0, &fcontext);
  while (sfp != NULL) {
    ivals = fcontext.ivals;
    numivals = fcontext.numivals;
    if (ivals != NULL && numivals > 0) {

      idx = (numivals - 1) * 2;
      start = ivals [idx];
      stop = ivals [idx + 1];
      if (stop >= from && stop <= to) {

        /*
        start = ivals [0] + 1;
        stop = ivals [idx + 1] + 1;
        */
        pdp = (PubdescPtr) sfp->data.value.ptrvalue;
        rbp = AddPub (awp, &(awp->pubhead), pdp);
        if (rbp != NULL) {

          rbp->entityID = fcontext.entityID;
          rbp->itemID = fcontext.itemID;
          rbp->itemtype = OBJ_SEQFEAT;

          irp = (IntRefBlockPtr) rbp;
          irp->loc = SeqLocMerge (target, sfp->location, NULL, FALSE, TRUE, FALSE);
          alp = GetAuthListPtr (pdp, NULL);
          if (alp != NULL) {
            irp->authstr = GetAuthorsString (awp->format, alp);
          }
          irp->index = fcontext.index;
        }
      }
    }

    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_PUB, 0, &fcontext);
  }
}

static Boolean LIBCALLBACK GetRefsOnSeg (
  SeqLocPtr slp,
  SeqMgrSegmentContextPtr context
)

{
  Asn2gbWorkPtr  awp;
  BioseqPtr      bsp = NULL;
  Uint2          entityID;
  Int4           from;
  SeqLocPtr      loc;
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

  GetRefsOnBioseq (awp, awp->target, bsp, from, to);

  BioseqUnlock (bsp);

  return TRUE;
}

static void AddReferenceBlock (
  Asn2gbWorkPtr awp,
  Asn2gbSectionPtr asp
)

{
  Asn2gbJobPtr       ajp;
  BioseqPtr          bsp;
  Boolean            excise;
  ValNodePtr         head = NULL;
  Int2               i;
  IntRefBlockPtr     irp;
  IntRefBlockPtr     lastirp;
  ReferenceBlockPtr  lastrbp;
  ValNodePtr         next;
  Int2               numReferences;
  ValNodePtr         PNTR prev;
  ReferenceBlockPtr  rbp;
  ReferenceBlockPtr  PNTR referenceArray;
  SeqLocPtr          slp;
  BioseqPtr          target;
  ValNodePtr         vnp;

  if (awp == NULL || asp == NULL) return;
  ajp = awp->ajp;
  if (ajp == NULL) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;

  /* collect publications on bioseq */

  awp->pubhead = NULL;
  GetRefsOnBioseq (awp, bsp, bsp, awp->from, awp->to);
  target = bsp;

  if (bsp->repr == Seq_repr_seg) {

    /* collect publications on remote segments in MASTER_STYLE */

    SeqMgrExploreSegments (bsp, (Pointer) awp, GetRefsOnSeg);
    target = awp->target;
  }

  head = awp->pubhead;
  awp->pubhead = NULL;

  if (head == NULL) return;

  /* sort by existing serial, then pub/unpub/sites/sub, then date */

  head = SortValNode (head, SortReferences);

  if (awp->ssp != NULL) {

    /* add seq-submit citation */

    rbp = (ReferenceBlockPtr) MemNew (sizeof (IntRefBlock));
    if (rbp != NULL) {

      rbp->blocktype = REFERENCE_BLOCK;
      rbp->section = awp->currsection;
      rbp->serial = INT2_MAX;
      rbp->category = REF_CAT_SUB;

      rbp->entityID = ajp->entityID;
      rbp->itemID = 1;
      rbp->itemtype = OBJ_SEQSUB_CIT;

      if (awp->format == DDBJ_FMT || awp->format == DDBJPEPT_FMT) {

        /* for DDBJ, add seq-submit citation to beginning of list */

        vnp = ValNodeNew (NULL);
        if (vnp != NULL) {
          vnp->choice = 0;
          vnp->data.ptrvalue = (VoidPtr) rbp;
          vnp->next = head;
          head = vnp;
        }

      } else {

        /* for GENBANK and EMBL add seq-submit citation to end of list */

        ValNodeAddPointer (&head, 0, rbp);
      }
    }
  }

  /* unique references, excise duplicates from list */

  prev = &(head);
  vnp = head;
  lastrbp = NULL;
  while (vnp != NULL) {
    excise = FALSE;
    next = vnp->next;
    rbp = (ReferenceBlockPtr) vnp->data.ptrvalue;
    if (lastrbp != NULL) {
      lastirp = (IntRefBlockPtr) lastrbp;
      if (rbp != NULL) {
        irp = (IntRefBlockPtr) rbp;
        if (lastrbp->pmid != 0 && rbp->pmid != 0) {
          if (lastrbp->pmid == rbp->pmid) {
            excise = TRUE;
          }
        } else if (lastrbp->muid != 0 && rbp->muid != 0) {
          if (lastrbp->muid == rbp->muid) {
            excise = TRUE;
          }
        } else if (lastrbp->uniquestr != NULL && rbp->uniquestr != NULL) {
          if (StringICmp (lastrbp->uniquestr, rbp->uniquestr) == 0) {
            if (SeqLocCompare (irp->loc, lastirp->loc) == SLC_A_EQ_B) {
              if (StringICmp (irp->authstr, lastirp->authstr) == 0) {
                /*
                excise = TRUE;
                */
              }
            }
          }
        }
      }
    }
    if (excise) {
      *prev = vnp->next;
      vnp->next = NULL;

      /* combine locations of duplicate references */

      irp = (IntRefBlockPtr) rbp;
      lastirp = (IntRefBlockPtr) lastrbp;
      if (lastirp != NULL) {
        slp = SeqLocMerge (target, lastirp->loc, irp->loc, FALSE, TRUE, FALSE);
        lastirp->loc = SeqLocFree (lastirp->loc);
        lastirp->loc = slp;
      }
      if (irp != NULL && lastirp != NULL) {
        if (irp->justuids && (! lastirp->justuids)) {
          if (lastirp->fig == NULL) {
            lastirp->fig = StringSaveNoNull (irp->fig);
          }
          if (lastirp->maploc == NULL) {
            lastirp->maploc = StringSaveNoNull (irp->maploc);
          }
          lastirp->poly_a = irp->poly_a;
        }
      }

      /* and remove duplicate reference */

      MemFree (rbp->uniquestr);
      DateFree (irp->date);
      SeqLocFree (irp->loc);
      MemFree (irp->authstr);
      MemFree (irp->fig);
      MemFree (irp->maploc);
      MemFree (rbp);
      ValNodeFree (vnp);

    } else {

      prev = &(vnp->next);
      lastrbp = rbp;
    }
    vnp = next;
  }

  /* assign serial numbers */

  for (vnp = head, i = 1; vnp != NULL; vnp = vnp->next, i++) {
    rbp = (ReferenceBlockPtr) vnp->data.ptrvalue;
    if (rbp != NULL) {
      rbp->serial = i;
    }
  }

  /* allocate reference array for this section */

  numReferences = i - 1;
  asp->numReferences = numReferences;

  if (numReferences > 0) {
    referenceArray = (ReferenceBlockPtr PNTR) MemNew (sizeof (ReferenceBlockPtr) * (numReferences + 1));
    asp->referenceArray = referenceArray;

    if (referenceArray != NULL) {

      /* fill in reference array */

      for (vnp = head, i = 0; vnp != NULL && i < numReferences; vnp = vnp->next, i++) {
        referenceArray [i] = (ReferenceBlockPtr) vnp->data.ptrvalue;
      }
    }
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
}

static CharPtr PrintDate (
  DatePtr dp
)

{
  Char    buf [30];
  size_t  len;

  if (dp == NULL) return NULL;

  if (DatePrint (dp, buf)) {
    if (StringICmp (buf, "Not given") != 0) {
      len = StringLen (buf);
      if (len > 0) {
        if (buf [len - 1] == '\n') {
          if (buf [len - 2] == '.') {
            buf [len - 2] = '\0';
          } else {
            buf [len - 1] = '\0';
          }
        }
      }
      return StringSave (buf);
    }
  }

  return NULL;
}

static void AddHistCommentString (
  CharPtr prefix,
  CharPtr suffix,
  DatePtr dp,
  SeqIdPtr ids
)

{
  Char        buf [256];
  ValNodePtr  head = NULL;
  SeqIdPtr    sip;
  CharPtr     str;
  CharPtr     strd;

  if (dp == NULL || ids == NULL || prefix == NULL || suffix == NULL) return;

  strd = PrintDate (dp);
  if (strd == NULL) {
    strd = StringSave ("?");
  }

  sprintf (buf, "%s %s %s", prefix, strd, suffix);
  ValNodeCopyStr (&head, 0, buf);

  MemFree (strd);

  for (sip = ids; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_GI) {
     sprintf (buf, " gi:%ld", (long) sip->data.intvalue);
      ValNodeCopyStr (&head, 0, buf);
    }
  }

  str = MergeValNodeStrings (head);

  gb_AddString (NULL, str, NULL, TRUE, FALSE, TRUE);

  MemFree (str);
  ValNodeFreeData (head);
}

static void AddHTGSCommentString (
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
  CharPtr      str;

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
      sprintf (buffer, "* NOTE: This record contains %ld individual~", (long) (num_s - num_g));
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
      sprintf (buffer, " It currently~* consists of %ld contigs. The true order of the pieces~", (long) (num_s - num_g));
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
      sprintf (buffer, " It currently~* consists of %ld contigs. Gaps between the contigs~", (long) (num_s - num_g));
      ValNodeCopyStr (&head, 0, buffer);
      ValNodeCopyStr (&head, 0, "* are represented as runs of N. The order of the pieces~");
      ValNodeCopyStr (&head, 0, "* is believed to be correct as given, however the sizes~");
      ValNodeCopyStr (&head, 0, "* of the gaps between them are based on estimates that have~");
      ValNodeCopyStr (&head, 0, "* provided by the submittor.");
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

  str = MergeValNodeStrings (head);

  gb_AddString (NULL, str, NULL, TRUE, FALSE, TRUE);

  MemFree (str);
  ValNodeFreeData (head);
}

static CharPtr GetStrForUserObject (
  UserObjectPtr uop
)

{
  ObjectIdPtr   oip;
  UserFieldPtr  ufp, tmp, u;
  CharPtr       ptr = NULL;
  Int2          i = 0;
  Char          p [13];
	
  if (uop == NULL) return NULL;
  if ((oip = uop->type) == NULL) return NULL;
  if (StringCmp (oip->str, "RefGeneTracking") != 0) return NULL;
  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    oip = ufp->label;
    if (StringCmp (oip->str, "Assembly") == 0) break;
  }
  if (ufp && ufp->choice == 11) {
    for (tmp = ufp->data.ptrvalue; tmp != NULL; tmp = tmp->next, i++) continue;
    ptr = MemNew (StringLen ("This reference sequence was derived from ") + 10*i + 1);
    if (ptr == NULL) return NULL;
    sprintf (ptr, "This reference sequence was derived from ");
    for (tmp = ufp->data.ptrvalue; tmp != NULL; tmp = tmp->next) {
      for (u = tmp->data.ptrvalue; u != NULL; u = u->next) {
        oip = u->label;
        if (StringCmp (oip->str, "accession") == 0) break;
      }
      if (u != NULL && tmp->next) {
        sprintf (p, "%s, ", u->data.ptrvalue);
      } else {
        sprintf (p, "%s.~~", u->data.ptrvalue);
      }
      StringCat (ptr, p);
    }
  }
  return ptr;
}

static void AddCommentBlock (
  Asn2gbWorkPtr awp
)

{
  BioseqPtr          bsp;
  Char               buf [128];
  CommentBlockPtr    cbp;
  DbtagPtr           dbt;
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext;
  Boolean            first = TRUE;
  SeqHistPtr         hist;
  MolInfoPtr         mip;
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp;
  SeqIdPtr           sip;
  CharPtr            str;
  UserObjectPtr      uop;

  if (awp == NULL) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;

  /* first show GSDB sequence identifier */

  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_GENERAL) {
      dbt = (DbtagPtr) sip->data.ptrvalue;
      if (dbt != NULL && StringCmp (dbt->db, "GSDB") == 0 && dbt->tag != NULL) {
        cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
        if (cbp != NULL) {
          cbp->first = first;
          first = FALSE;

          sprintf (buf, "GSDB:S:%ld", (long) dbt->tag->id);

          if (cbp->first) {
            gb_StartPrint (awp->format, TRUE, 0, 12, "COMMENT", 13, 5, 5, "CC", TRUE);
          } else {
            gb_StartPrint (awp->format, TRUE, 0, 12, NULL, 13, 5, 5, "CC", FALSE);
          }

          /* CheckEndPunctuation, ConvertDoubleQuotes, and ExpandTildes already taken into account */

          ff_AddString (buf);

          cbp->string = gb_MergeString (TRUE);
        }
      }
    }
  }

  /* Seq-hist results in allocated comment string */

  hist = bsp->hist;
  if (hist != NULL) {

    if (hist->replaced_by_ids != NULL && hist->replaced_by_date != NULL) {

      cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
      if (cbp != NULL) {
        cbp->first = first;
        first = FALSE;

        if (cbp->first) {
          gb_StartPrint (awp->format, TRUE, 0, 12, "COMMENT", 13, 5, 5, "CC", TRUE);
        } else {
          gb_StartPrint (awp->format, TRUE, 0, 12, NULL, 13, 5, 5, "CC", FALSE);
        }

        AddHistCommentString ("[WARNING] On", "this sequence was replaced by a newer version",
                              hist->replaced_by_date, hist->replaced_by_ids);

        cbp->string = gb_MergeString (TRUE);
      }
    }

    if (hist->replace_ids != NULL && hist->replace_date != NULL) {

      cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
      if (cbp != NULL) {
        cbp->first = first;
        first = FALSE;

        if (cbp->first) {
          gb_StartPrint (awp->format, TRUE, 0, 12, "COMMENT", 13, 5, 5, "CC", TRUE);
        } else {
          gb_StartPrint (awp->format, TRUE, 0, 12, NULL, 13, 5, 5, "CC", FALSE);
        }

        AddHistCommentString ("On", "this sequence version replaced",
                              hist->replace_date, hist->replace_ids);

        cbp->string = gb_MergeString (TRUE);
      }
    }

  }

  /* just save IDs for comment, maploc, and region descriptors */

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_comment, &dcontext);
  while (sdp != NULL) {
    cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
    if (cbp != NULL) {
      cbp->entityID = dcontext.entityID;
      cbp->itemID = dcontext.itemID;
      cbp->itemtype = OBJ_SEQDESC;
      cbp->first = first;
      first = FALSE;
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
      first = FALSE;
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_maploc, &dcontext);
  }

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_region, &dcontext);
  while (sdp != NULL) {
    cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
    if (cbp != NULL) {
      cbp->entityID = dcontext.entityID;
      cbp->itemID = dcontext.itemID;
      cbp->itemtype = OBJ_SEQDESC;
      cbp->first = first;
      first = FALSE;
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_region, &dcontext);
  }

  /* RefSeq results in allocated comment string */

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_user, &dcontext);
  if (sdp != NULL) {

    uop = (UserObjectPtr) sdp->data.ptrvalue;
    if (uop != NULL) {

      str = GetStrForUserObject (uop);
      if (str != NULL) {

        cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
        if (cbp != NULL) {

          /*
          cbp->entityID = dcontext.entityID;
          cbp->itemID = dcontext.itemID;
          cbp->itemtype = OBJ_SEQDESC;
          */
          cbp->first = first;
          first = FALSE;

          if (cbp->first) {
            gb_StartPrint (awp->format, TRUE, 0, 12, "COMMENT", 13, 5, 5, "CC", TRUE);
          } else {
            gb_StartPrint (awp->format, TRUE, 0, 12, NULL, 13, 5, 5, "CC", FALSE);
          }

          gb_AddString (NULL, str, NULL, TRUE, FALSE, TRUE);

          cbp->string = gb_MergeString (TRUE);
        }
      }
    }
  }

  /* HTGS results in allocated comment string */

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
  if (sdp != NULL) {

    mip = (MolInfoPtr) sdp->data.ptrvalue;
    if (mip != NULL && (mip->tech == MI_TECH_htgs_0 ||
                        mip->tech == MI_TECH_htgs_1 ||
                        mip->tech == MI_TECH_htgs_2)) {

      cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
      if (cbp != NULL) {

        /*
        cbp->entityID = dcontext.entityID;
        cbp->itemID = dcontext.itemID;
        cbp->itemtype = OBJ_SEQDESC;
        */
        cbp->first = first;
        first = FALSE;

        if (cbp->first) {
          gb_StartPrint (awp->format, TRUE, 0, 12, "COMMENT", 13, 5, 5, "CC", TRUE);
        } else {
          gb_StartPrint (awp->format, TRUE, 0, 12, NULL, 13, 5, 5, "CC", FALSE);
        }

        AddHTGSCommentString (bsp, mip);

        cbp->string = gb_MergeString (TRUE);
      }
    }
  }

  /* add comment features that are full length */

  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_COMMENT, 0, &fcontext);
  while (sfp != NULL) {
    if (fcontext.left == 0 && fcontext.right == bsp->length - 1) {
      cbp = (CommentBlockPtr) Asn2gbAddBlock (awp, COMMENT_BLOCK, sizeof (CommentBlock));
      if (cbp != NULL) {
        cbp->entityID = fcontext.entityID;
        cbp->itemID = fcontext.itemID;
        cbp->itemtype = OBJ_SEQFEAT;
        cbp->first = first;
        first = FALSE;
      }
    }
    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_COMMENT, 0, &fcontext);
  }
}

static void AddFeatHeaderBlock (
  Asn2gbWorkPtr awp
)

{
  BaseBlockPtr  bbp;

  if (awp == NULL) return;

  bbp = Asn2gbAddBlock (awp, FEATHEADER_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  if (awp->format == FTABLE_FMT) return;

  gb_StartPrint (awp->format, TRUE, 0, 12, "FEATURES", 22, 5, 0, "FH", TRUE);

  if (awp->format == EMBL_FMT || awp->format == EMBLPEPT_FMT) {
    ff_AddString ("Key");
    TabToColumn (22);
  }

  gb_AddString (NULL, "Location/Qualifiers", NULL, FALSE, FALSE, FALSE);

  if (awp->format == EMBL_FMT || awp->format == EMBLPEPT_FMT) {
    NewContLine();
  }

  bbp->string = gb_MergeString (TRUE);
}

static int SourceCompare (
  BioSourcePtr biop1,
  BioSourcePtr biop2
)

{
  int           compare;
  OrgModPtr     omp1, omp2;
  OrgNamePtr    onp1, onp2;
  OrgRefPtr     orp1, orp2;
  SubSourcePtr  ssp1, ssp2;

  if (biop1 == NULL || biop2 == NULL) return 0;
  if (biop1->is_focus && (! biop2->is_focus)) return 1;
  if (biop2->is_focus && (! biop1->is_focus)) return -1;

  orp1 = biop1->org;
  orp2 = biop2->org;
  if (orp1 == NULL || orp2 == NULL) return 0;

  compare = StringICmp (orp1->taxname, orp2->taxname);
  if (compare != 0) return compare;

  onp1 = orp1->orgname;
  onp2 = orp2->orgname;
  if (onp1 != NULL && onp2 == NULL) return 1;
  if (onp2 != NULL && onp1 == NULL) return -1;

  for (omp1 = onp1->mod, omp2 = onp2->mod;
       omp1 != NULL && omp2 != NULL;
       omp1 = omp1->next, omp2 = omp2->next) {
    if (omp1->subtype < omp2->subtype) return 1;
    if (omp2->subtype < omp1->subtype) return -1;
    compare = StringICmp (omp1->subname, omp2->subname);
    if (compare != 0) return compare;
  }
  if (omp1 != NULL && omp2 == NULL) return 1;
  if (omp2 != NULL && omp1 == NULL) return -1;

  for (ssp1 = biop1->subtype, ssp2 = biop2->subtype;
       ssp1 != NULL && ssp2 != NULL;
       ssp1 = ssp1->next, ssp2 = ssp2->next) {
    if (ssp1->subtype < ssp2->subtype) return 1;
    if (ssp2->subtype < ssp1->subtype) return -1;
    compare = StringICmp (ssp1->name, ssp2->name);
    if (compare != 0) return compare;
  }
  if (ssp1 != NULL && ssp2 == NULL) return 1;
  if (ssp2 != NULL && ssp1 == NULL) return -1;

  return 0;
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
  BioSourcePtr biop
)

{
  BaseBlockPtr    bbp;
  DbtagPtr        dbt;
  Uint2           hash;
  Int2            idx;
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
  isp->is_focus = biop->is_focus;

  orp = biop->org;
  if (orp == NULL) return bbp;

  isp->orghash = ComputeSourceHash (orp->taxname, 0);

  hash = 0;
  onp = orp->orgname;
  if (onp != NULL) {
    for (omp = onp->mod; omp != NULL; omp = omp->next) {
      subtype = omp->subtype;
      if (subtype == 254) {
        subtype = 24;
      } else if (subtype == 255) {
        subtype = 25;
      }
      if (subtype < 26) {
        idx = orgModToSourceIdx [subtype];
        if (idx > 0 && idx < ASN2GNBK_TOTAL_SOURCE) {
          str = asn2gnbk_source_quals [idx].name;
          hash = ComputeSourceHash (str, hash);
          hash = ComputeSourceHash (omp->subname, hash);
        }
      }
    }
  }
  isp->modhash = hash;

  hash = 0;
  for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
    subtype = ssp->subtype;
    if (subtype == 255) {
      subtype = 24;
    }
    if (subtype < 25) {
      idx = subSourceToSourceIdx [subtype];
      if (idx > 0 && idx < ASN2GNBK_TOTAL_SOURCE) {
        str = asn2gnbk_source_quals [idx].name;
        hash = ComputeSourceHash (str, hash);
        hash = ComputeSourceHash (ssp->name, hash);
      }
    }
  }
  isp->subhash = hash;

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

  return bbp;
}

static int LIBCALLBACK SortSources (
  VoidPtr ptr1,
  VoidPtr ptr2
)

{
  Uint4           diff;
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

  return 0;
}

static void GetSourcesOnBioseq (
  Asn2gbWorkPtr awp,
  BioseqPtr target,
  BioseqPtr bsp,
  Int4 from,
  Int4 to
)

{
  BaseBlockPtr       bbp;
  BioSourcePtr       biop;
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext;
  Int2               idx;
  IntSrcBlockPtr     isp;
  Int4Ptr            ivals;
  Int2               numivals;
  Boolean            okay;
  SeqDescrPtr        sdp;
  SeqFeatPtr         sfp;
  SeqInt             sint;
  Int4               start;
  Int4               stop;
  ValNode            vn;
  ValNodePtr         vnp;

  if (awp == NULL || target == NULL || bsp == NULL) return;

  /* full length loc for descriptors */

  sint.from = 0;
  sint.to = bsp->length - 1;
  sint.strand = Seq_strand_plus;
  sint.id = SeqIdStripLocus (SeqIdDup (SeqIdFindBest (bsp->id, 0)));
  vn.choice = SEQLOC_INT;
  vn.data.ptrvalue = (Pointer) &sint;
  vn.next = NULL;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  while (sdp != NULL) {

    /* check if descriptor on part already added on segmented bioseq */

    okay = TRUE;
    for (vnp = awp->srchead; vnp != NULL && okay; vnp = vnp->next) {
      bbp = (BaseBlockPtr) vnp->data.ptrvalue;
      if (bbp != NULL) {
        if (bbp->entityID == dcontext.entityID &&
            bbp->itemID == dcontext.itemID &&
            bbp->itemtype == OBJ_SEQDESC) {
          okay = FALSE;
        }
      }
    }

    if (okay) {
      biop = (BioSourcePtr) sdp->data.ptrvalue;
      bbp = AddSource (awp, &(awp->srchead), biop);
      if (bbp != NULL) {

        bbp->entityID = dcontext.entityID;
        bbp->itemID = dcontext.itemID;
        bbp->itemtype = OBJ_SEQDESC;

        isp = (IntSrcBlockPtr) bbp;
        isp->loc = SeqLocMerge (target, &vn, NULL, FALSE, TRUE, FALSE);
      }
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_source, &dcontext);
  }

  SeqIdFree (sint.id);

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
      if (stop >= from && stop <= to) {

        /*
        start = ivals [0] + 1;
        stop = ivals [idx + 1] + 1;
        */
        biop = (BioSourcePtr) sfp->data.value.ptrvalue;
        bbp = AddSource (awp, &(awp->srchead), biop);
        if (bbp != NULL) {

          bbp->entityID = fcontext.entityID;
          bbp->itemID = fcontext.itemID;
          bbp->itemtype = OBJ_SEQFEAT;

          isp = (IntSrcBlockPtr) bbp;
          isp->loc = SeqLocMerge (target, sfp->location, NULL, FALSE, TRUE, FALSE);
        }
      }
    }

    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_BIOSRC, 0, &fcontext);
  }
}

static Boolean LIBCALLBACK GetSourcesOnSeg (
  SeqLocPtr slp,
  SeqMgrSegmentContextPtr context
)

{
  Asn2gbWorkPtr  awp;
  BioseqPtr      bsp = NULL;
  Uint2          entityID;
  Int4           from;
  SeqLocPtr      loc;
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

  GetSourcesOnBioseq (awp, awp->target, bsp, from, to);

  BioseqUnlock (bsp);

  return TRUE;
}

static void AddSourceFeatBlock (
  Asn2gbWorkPtr awp
)

{
  BioseqPtr       bsp;
  Boolean         excise;
  ValNodePtr      head = NULL;
  IntSrcBlockPtr  isp;
  IntSrcBlockPtr  lastisp;
  ValNodePtr      next;
  ValNodePtr      PNTR prev;
  SeqLocPtr       slp;
  BioseqPtr       target;
  ValNodePtr      vnp;

  if (awp == NULL) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;

  /* collect biosources on bioseq */

  awp->srchead = NULL;
  GetSourcesOnBioseq (awp, bsp, bsp, awp->from, awp->to);
  target = bsp;

  if (bsp->repr == Seq_repr_seg) {

    /* collect biosources on remote segments in MASTER_STYLE */

    SeqMgrExploreSegments (bsp, (Pointer) awp, GetSourcesOnSeg);
    target = awp->target;
  }

  head = awp->srchead;
  awp->srchead = NULL;

  if (head == NULL) return;

  /* sort by ??? */

  head = SortValNode (head, SortSources);

  /* unique sources, excise duplicates from list */

  prev = &(head);
  vnp = head;
  lastisp = NULL;
  while (vnp != NULL) {
    excise = FALSE;
    next = vnp->next;
    isp = (IntSrcBlockPtr) vnp->data.ptrvalue;
    if (lastisp != NULL) {
      if (isp != NULL) {
        if (lastisp->is_focus == isp->is_focus &&
            lastisp->orghash == isp->orghash &&
            lastisp->xrfhash == isp->xrfhash) {

          /* check for identical modifiers */

          if (lastisp->modhash == isp->modhash &&
              lastisp->subhash == isp->subhash) {
            excise = TRUE;

          /* or modifiers only in lastisp (e.g., on part bioseq) */

          } else if (isp->modhash == 0 && isp->subhash == 0) {
            excise = TRUE;
          }
        }
      }
    }
    if (excise) {
      *prev = vnp->next;
      vnp->next = NULL;

      /* combine locations of duplicate sources */

      if (lastisp != NULL) {
        slp = SeqLocMerge (target, lastisp->loc, isp->loc, FALSE, TRUE, FALSE);
        lastisp->loc = SeqLocFree (lastisp->loc);
        lastisp->loc = slp;
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
}

static void GetFeatsOnCdsProduct (
  SeqFeatPtr sfp,
  BioseqPtr bsp,
  Asn2gbWorkPtr awp
)

{
  FeatureBlockPtr    fbp;
  IntFeatBlockPtr    ifp;
  Int4               lastleft;
  Int4               lastright;
  SeqAnnotPtr        lastsap;
  SeqFeatPtr         lastsfp;
  SeqMgrFeatContext  pcontext;
  SeqFeatPtr         prt;
  Boolean            suppress;

  if (sfp == NULL || awp == NULL) return;
  if (bsp == NULL || (! ISA_aa (bsp->mol))) return;

  /* explore mat_peptides, sites, etc. */

  lastsfp = NULL;
  lastsap = NULL;
  lastleft = 0;
  lastright = 0;

  prt = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &pcontext);
  while (prt != NULL) {
    
    if (pcontext.featdeftype == FEATDEF_REGION ||
        pcontext.featdeftype == FEATDEF_SITE ||
        pcontext.featdeftype == FEATDEF_mat_peptide_aa ||
        pcontext.featdeftype == FEATDEF_sig_peptide_aa ||
        pcontext.featdeftype == FEATDEF_transit_peptide_aa) {

      if (pcontext.dnaStop >= awp->from && pcontext.dnaStop <= awp->to) {

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

        if (! suppress) {

          fbp = (FeatureBlockPtr) Asn2gbAddBlock (awp, FEATURE_BLOCK, sizeof (IntFeatBlock));
          if (fbp != NULL) {

            fbp->entityID = pcontext.entityID;
            fbp->itemID = pcontext.itemID;
            fbp->itemtype = OBJ_SEQFEAT;
            fbp->featdeftype = pcontext.featdeftype;
            ifp = (IntFeatBlockPtr) fbp;
            ifp->mapToNuc = TRUE;
            ifp->mapToProt = FALSE;
          }
        }

        lastsfp = prt;
        lastsap = pcontext.sap;
        lastleft = pcontext.left;
        lastright = pcontext.right;

      }
    }
    prt = SeqMgrGetNextFeature (bsp, prt, 0, 0, &pcontext);
  }
}

static Boolean LIBCALLBACK GetFeatsOnBioseq (
  SeqFeatPtr sfp,
  SeqMgrFeatContextPtr fcontext
)

{
  Asn2gbWorkPtr      awp;
  BioseqPtr          bsp;
  SeqMgrDescContext  dcontext;
  FeatureBlockPtr    fbp;
  IntCdsBlockPtr     icp;
  Int2               idx;
  IntFeatBlockPtr    ifp;
  Int4Ptr            ivals;
  SeqAnnotPtr        lastsap;
  SeqFeatPtr         lastsfp;
  Int2               numivals;
  PubdescPtr         pdp;
  SeqDescrPtr        sdp;
  Int4               start;
  Int4               stop;

  if (sfp == NULL || fcontext == NULL) return FALSE;
  awp = (Asn2gbWorkPtr) fcontext->userdata;
  if (awp == NULL) return FALSE;

  if (fcontext->seqfeattype == SEQFEAT_BIOSRC ||
      fcontext->seqfeattype == SEQFEAT_PUB) return TRUE;

  /* DDBJ does not want to show gene features */

  if (fcontext->seqfeattype == SEQFEAT_GENE &&
      (awp->format == DDBJ_FMT || awp->format == DDBJPEPT_FMT)) return TRUE;

  /* suppress comment features that are full length */

  if (fcontext->seqfeattype == SEQFEAT_COMMENT &&
      fcontext->left == awp->from && fcontext->right == awp->to) return TRUE;

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
        bsp = BioseqFindFromSeqLoc (sfp->product);
        GetFeatsOnCdsProduct (sfp, bsp, awp);
      }

      return TRUE;
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

  /* !!! if RELEASE_MODE, verify that features have all mandatory qualifiers !!! */

  awp->lastsfp = sfp;
  awp->lastsap = fcontext->sap;
  awp->lastleft = fcontext->left;
  awp->lastright = fcontext->right;

  if (fcontext->seqfeattype == SEQFEAT_CDREGION) {
    fbp = (FeatureBlockPtr) Asn2gbAddBlock (awp, FEATURE_BLOCK, sizeof (IntCdsBlock));
  } else {
    fbp = (FeatureBlockPtr) Asn2gbAddBlock (awp, FEATURE_BLOCK, sizeof (IntFeatBlock));
  }
  if (fbp == NULL) return TRUE;

  fbp->entityID = fcontext->entityID;
  fbp->itemID = fcontext->itemID;
  fbp->itemtype = OBJ_SEQFEAT;
  fbp->featdeftype = fcontext->featdeftype;
  ifp = (IntFeatBlockPtr) fbp;
  ifp->mapToNuc = FALSE;
  ifp->mapToProt = FALSE;

  if (fcontext->seqfeattype != SEQFEAT_CDREGION) return TRUE;

  /* if CDS, collect more information from product bioseq */

  bsp = BioseqFindFromSeqLoc (sfp->product);
  if (bsp == NULL || (! ISA_aa (bsp->mol))) return TRUE;

  ifp->isCDS = TRUE;
  icp = (IntCdsBlockPtr) ifp;

  /* first explore pubs to pick up figure and maploc */

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

  /* then explore mat_peptides, sites, etc. */

  GetFeatsOnCdsProduct (sfp, bsp, awp);

  return TRUE;
}

static Boolean LIBCALLBACK GetFeatsOnSeg (
  SeqLocPtr slp,
  SeqMgrSegmentContextPtr context
)

{
  Asn2gbWorkPtr  awp;
  BioseqPtr      bsp = NULL;
  Uint2          entityID;
  Int4           from;
  SeqLocPtr      loc;
  SeqIdPtr       sip;
  Int4           to;

  if (slp == NULL || context == NULL) return FALSE;
  awp = (Asn2gbWorkPtr) context->userdata;

  from = awp->from;
  to = awp->from;

  sip = SeqLocId (slp);
  if (sip == NULL) {
    loc = SeqLocFindNext (slp, NULL);
    if (loc != NULL) {
      sip = SeqLocId (loc);
    }
  }
  if (sip == NULL) return TRUE;

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

    awp->from = 0;
    awp->to = bsp->length - 1;
  }

  awp->lastsfp = NULL;
  awp->lastsap = NULL;
  awp->lastleft = 0;
  awp->lastright = 0;

  SeqMgrExploreFeatures (bsp, (Pointer) awp, GetFeatsOnBioseq, awp->slp, NULL, NULL);

  /* restore original from and to */

  awp->from = from;
  awp->to = to;

  BioseqUnlock (bsp);

  return TRUE;
}

static void AddFeatureBlock (
  Asn2gbWorkPtr awp
)

{
  BioseqPtr  bsp;

  if (awp == NULL) return;
  bsp = awp->parent;
  if (bsp == NULL) return;

  awp->lastsfp = NULL;
  awp->lastsap = NULL;
  awp->lastleft = 0;
  awp->lastright = 0;

  SeqMgrExploreFeatures (bsp, (Pointer) awp, GetFeatsOnBioseq, awp->slp, NULL, NULL);

  if (bsp->repr == Seq_repr_seg) {

    /* collect features on remote segments in MASTER_STYLE */

    SeqMgrExploreSegments (bsp, (Pointer) awp, GetFeatsOnSeg);
  }

}

static void AddBasecountBlock (
  Asn2gbWorkPtr awp
)

{
  BaseBlockPtr  bbp;

  if (awp == NULL) return;

  bbp = Asn2gbAddBlock (awp, BASECOUNT_BLOCK, sizeof (BaseBlock));
}

static void AddOriginBlock (
  Asn2gbWorkPtr awp
)

{
  BaseBlockPtr       bbp;
  BioseqPtr          bsp;
  Char               buf [67];
  SeqMgrDescContext  dcontext;
  GBBlockPtr         gbp;
  SeqDescrPtr        sdp;

  if (awp == NULL) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;

  if (awp->format == EMBL_FMT || awp->format == EMBLPEPT_FMT) return;

  bbp = Asn2gbAddBlock (awp, ORIGIN_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  if (awp->format == GENBANK_FMT || awp->format == GENPEPT_FMT) {

    buf [0] = '\0';

    sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_genbank, &dcontext);
    if (sdp != NULL) {
      gbp = (GBBlockPtr) sdp->data.ptrvalue;
      if (gbp != NULL && (! StringHasNoText (gbp->origin))) {
        StringNCpy_0 (buf, gbp->origin, sizeof (buf));
        bbp->entityID = dcontext.entityID;
        bbp->itemID = dcontext.itemID;
        bbp->itemtype = OBJ_SEQDESC;
      }
    }

    ff_StartPrint (0, 12, ASN2FF_GB_MAX, NULL);

    ff_AddString ("ORIGIN");
    TabToColumn (13);

    if (! StringHasNoText (buf)) {
      gb_AddString (NULL, buf, NULL, TRUE, FALSE, FALSE);
    }
  }

  bbp->string = gb_MergeString (TRUE);
}

#define BASES_PER_BLOCK 1200

static void AddSequenceBlock (
  Asn2gbWorkPtr awp
)

{
  BioseqPtr         bsp;
  Int4              len;
  SequenceBlockPtr  sbp;
  Int4              start;
  Int4              stop;

  if (awp == NULL) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;

  if (awp->slp != NULL) {
    len = SeqLocLen (awp->slp);
  } else {
    len = bsp->length;
  }

  /* populate individual sequence blocks for given range */

  for (start = 0; start < len; start += BASES_PER_BLOCK) {
    sbp = (SequenceBlockPtr) Asn2gbAddBlock (awp, SEQUENCE_BLOCK, sizeof (SequenceBlock));
    if (sbp == NULL) continue;

    sbp->entityID = bsp->idx.entityID;
    sbp->itemID = bsp->idx.itemID;
    sbp->itemtype = OBJ_BIOSEQ;

    stop = start + BASES_PER_BLOCK;
    if (stop >= len) {
      stop = len;
    }

    sbp->start = start;
    sbp->stop = stop;
  }
}

static void AddContigBlock (
  Asn2gbWorkPtr awp
)

{
  BaseBlockPtr  bbp;

  if (awp == NULL) return;

  bbp = Asn2gbAddBlock (awp, CONTIG_BLOCK, sizeof (BaseBlock));
}

static void AddSlashBlock (
  Asn2gbWorkPtr awp
)

{
  BaseBlockPtr  bbp;

  if (awp == NULL) return;

  bbp = Asn2gbAddBlock (awp, SLASH_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  gb_StartPrint (awp->format, TRUE, 0, 0, NULL, 0, 0, 0, NULL, FALSE);

  gb_AddString (NULL, "//", NULL, FALSE, FALSE, FALSE);

  bbp->string = gb_MergeString (TRUE);
}


/* ********************************************************************** */

/* DoOneSection builds a single report for one bioseq or segment */

static void DoOneSection (
  BioseqPtr target,
  BioseqPtr parent,
  BioseqPtr bsp,
  SeqLocPtr slp,
  Uint2 seg,
  Int4 from,
  Int4 to,
  Asn2gbWorkPtr awp
)

{
  Asn2gbSectionPtr     asp;
  SeqMgrBioseqContext  bcontext;
  BaseBlockPtr         PNTR blockArray;
  Int4                 i;
  Int4                 numBlocks;
  Int4                 numsegs = 0;
  ValNodePtr           vnp;

  if (target == NULL || parent == NULL || bsp == NULL || awp == NULL) return;

  asp = Asn2gbAddSection (awp);
  if (asp == NULL) return;

  if (SeqMgrGetBioseqContext (parent, &bcontext)) {
    numsegs = bcontext.numsegs;
  }

  /* set working data fields */

  awp->target = target;
  awp->parent = parent;
  awp->bsp = bsp;
  awp->slp = slp;
  awp->seg = seg;
  awp->numsegs = numsegs;
  awp->from = from;
  awp->to = to;

  /* initialize empty blockList for this section */

  awp->blockList = NULL;
  awp->lastblock = NULL;

  /* and store section data into section fields */

  asp->target = target;
  asp->bsp = bsp;
  asp->slp = slp;
  asp->seg = seg;
  asp->numsegs = numsegs;
  asp->from = from;
  asp->to = to;

  asp->spp = NULL;

  asp->blockArray = NULL;
  asp->numBlocks = 0;

  /* start exploring and populating paragraphs */

  if (awp->format == FTABLE_FMT) {

    AddFeatHeaderBlock (awp);
    AddFeatureBlock (awp);

  } else {

    AddLocusBlock (awp);

    if (awp->format == GENBANK_FMT || awp->format == GENPEPT_FMT) {

      AddDeflineBlock (awp);
      AddAccessionBlock (awp);

      if (ISA_na (bsp->mol)) {
        AddVersionBlock (awp);
      }

      if (ISA_aa (bsp->mol)) {
        /* AddPidBlock (awp); */
        AddDbsourceBlock (awp);
      }

    } else if (awp->format == EMBL_FMT || awp->format == EMBLPEPT_FMT) {

      AddAccessionBlock (awp);

      if (ISA_na (bsp->mol)) {
        AddVersionBlock (awp);
      }

      if (ISA_aa (bsp->mol)) {
        /* AddPidBlock (awp); */
        AddDbsourceBlock (awp);
      }

      AddDateBlock (awp);

      AddDeflineBlock (awp);
    }

    AddKeywordsBlock (awp);

    if (awp->format == GENBANK_FMT || awp->format == GENPEPT_FMT) {
      AddSegmentBlock (awp);
    }

    AddSourceBlock (awp);
    AddOrganismBlock (awp);

    AddReferenceBlock (awp, asp);
    AddCommentBlock (awp);

    AddFeatHeaderBlock (awp);
    AddSourceFeatBlock (awp);

    if (awp->style == CONTIG_STYLE) {
      AddContigBlock (awp);
    } else {
      AddFeatureBlock (awp);

      AddBasecountBlock (awp);
      AddOriginBlock (awp);

      AddSequenceBlock (awp);
    }

    AddSlashBlock (awp);
  }

  /* allocate block array for this section */

  numBlocks = ValNodeLen (awp->blockList);
  asp->numBlocks = numBlocks;

  if (numBlocks > 0) {
    blockArray = (BaseBlockPtr PNTR) MemNew (sizeof (BaseBlockPtr) * (numBlocks + 1));
    asp->blockArray = blockArray;

    if (blockArray != NULL) {
      for (vnp = awp->blockList, i = 0; vnp != NULL; vnp = vnp->next, i++) {
        blockArray [i] = (BaseBlockPtr) vnp->data.ptrvalue;
      }
    }
  }

  /* free blockList, but leave data, now pointed to by blockArray elements */

  awp->blockList = ValNodeFree (awp->blockList);
  awp->lastblock = NULL;

  (awp->currsection)++;
}


/*
the following functions handle various kinds of input, all calling
DoOneSection once for each component that gets its own report
*/

static Boolean LIBCALLBACK Asn2Seg (
  SeqLocPtr slp,
  SeqMgrSegmentContextPtr context
)

{
  Asn2gbWorkPtr  awp;
  BioseqPtr      bsp = NULL;
  Uint2          entityID;
  Int4           from;
  SeqLocPtr      loc;
  BioseqPtr      parent;
  SeqIdPtr       sip;
  Int4           to;

  if (slp == NULL || context == NULL) return FALSE;
  awp = (Asn2gbWorkPtr) context->userdata;

  parent = context->parent;

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

    parent = bsp;
    from = 0;
    to = bsp->length - 1;
  }

  DoOneSection (bsp, parent, bsp, /* slp */ NULL, context->index, from, to, awp);

  BioseqUnlock (bsp);

  return TRUE;
}

static void DoOneBioseq (
  BioseqPtr bsp,
  Asn2gbWorkPtr awp
)

{
  Asn2gbJobPtr          ajp;
  SeqMgrSegmentContext  context;
  Int4                  from;
  BioseqPtr             parent;
  Int4                  to;

  if (bsp == NULL || awp == NULL) return;
  ajp = awp->ajp;
  if (ajp == NULL) return;

  if (bsp->repr == Seq_repr_seg) {

    /* this is a segmented bioseq */

    if (awp->style == NORMAL_STYLE) {

      /* show all segments individually */

      SeqMgrExploreSegments (bsp, (Pointer) awp, Asn2Seg);

    } else {

      /* show as single bioseq */

      parent = bsp;
      from = 0;
      to = bsp->length - 1;

      DoOneSection (parent, parent, bsp, ajp->slp, 0, from, to, awp);
    }

  } else if (bsp->repr == Seq_repr_raw ||
             bsp->repr == Seq_repr_const ||
             bsp->repr == Seq_repr_delta ||
             bsp->repr == Seq_repr_virtual) {

    parent = SeqMgrGetParentOfPart (bsp, &context);
    if (parent != NULL) {

      /* this is a part of an indexed segmented bioseq */

      from = context.cumOffset;
      to = from + context.to - context.from;

    } else {

      /* this is a regular non-segmented bioseq */

      parent = bsp;
      from = 0;
      to = bsp->length - 1;
    }

    DoOneSection (parent, parent, bsp, ajp->slp, 0, from, to, awp);
  }
}

static void DoPopPhyMutSet (
  SeqEntryPtr sep,
  Asn2gbWorkPtr awp
)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;

  if (sep == NULL || awp == NULL) return;

  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;
    if (bssp->_class == 7 || bssp->_class == 13 ||
        bssp->_class == 14 || bssp->_class == 15) {

      /* this is a pop/phy/mut set, catenate separate reports */

      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        DoPopPhyMutSet (sep, awp);
      }
      return;
    }
  }

  /* at most nuc-prot set, so get first bioseq */

  sep = FindNthBioseq (sep, 1);
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    DoOneBioseq (bsp, awp);
  }
}

/* ********************************************************************** */

/* public functions */

static int LIBCALLBACK SortParagraphByIDProc (
  VoidPtr vp1,
  VoidPtr vp2
)

{
  BaseBlockPtr  bbp1, bbp2;

  if (vp1 == NULL || vp2 == NULL) return 0;
  bbp1 = *((BaseBlockPtr PNTR) vp1);
  bbp2 = *((BaseBlockPtr PNTR) vp2);
  if (bbp1 == NULL || bbp2 == NULL) return 0;

  if (bbp1->entityID > bbp2->entityID) return 1;
  if (bbp1->entityID < bbp2->entityID) return -1;

  if (bbp1->itemtype > bbp2->itemtype) return 1;
  if (bbp1->itemtype < bbp2->itemtype) return -1;

  if (bbp1->itemID > bbp2->itemID) return 1;
  if (bbp1->itemID < bbp2->itemID) return -1;

  if (bbp1->paragraph > bbp2->paragraph) return 1;
  if (bbp1->paragraph < bbp2->paragraph) return -1;

  return 0;
}

NLM_EXTERN Asn2gbJobPtr asn2gnbk_setup (
  BioseqPtr bsp,
  BioseqSetPtr bssp,
  SeqLocPtr slp,
  FmtType format,
  ModType mode,
  StlType style
)

{
  Asn2gbJobPtr      ajp = NULL;
  Asn2gbSectionPtr  asp;
  Asn2gbWork        aw;
  BaseBlockPtr      bbp;
  BaseBlockPtr      PNTR blockArray;
  Uint2             entityID = 0;
  Int4              i;
  Int4              j;
  Int4              k;
  Int4              numBlocks;
  Int4              numSections;
  ObjMgrDataPtr     omdp;
  BaseBlockPtr      PNTR paragraphArray;
  BaseBlockPtr      PNTR paragraphByIDs;
  Int4              numParagraphs;
  Asn2gbSectionPtr  PNTR sectionArray;
  SubmitBlockPtr    sbp;
  SeqEntryPtr       sep;
  SeqSubmitPtr      ssp;
  ValNodePtr        vnp;

  if (bssp != NULL) {

    entityID = ObjMgrGetEntityIDForPointer (bssp);

  } else {

    if (bsp == NULL && slp != NULL) {
      bsp = BioseqFindFromSeqLoc (slp);
    }
    if (bsp == NULL) return NULL;

    entityID = ObjMgrGetEntityIDForPointer (bsp);
  }

  if (entityID == 0) return NULL;

  if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
    SeqMgrIndexFeatures (entityID, NULL);
  }

  ajp = (Asn2gbJobPtr) MemNew (sizeof (Asn2gbJob));
  if (ajp == NULL) return NULL;

  asn2ff_set_output (NULL, "\n");

  ajp->entityID = entityID;
  ajp->bsp = bsp;
  ajp->bssp = bssp;
  if (slp != NULL) {
    ajp->slp = AsnIoMemCopy ((Pointer) slp,
                             (AsnReadFunc) SeqLocAsnRead,
                             (AsnWriteFunc) SeqLocAsnWrite);
  /*
  } else if (bsp != NULL) {
    sep = SeqMgrGetSeqEntryForData (bsp);
    ajp->slp = CreateWholeInterval (sep);
  */
  } else {
    ajp->slp = NULL;
  }

  ajp->format = format;
  ajp->mode = mode;
  ajp->style = style;

  ajp->keepSeqPortOpen = FALSE;

  MemSet ((Pointer) (&aw), 0, sizeof (Asn2gbWork));
  aw.ajp = ajp;
  aw.entityID = entityID;

  aw.sectionList = NULL;
  aw.lastsection = NULL;

  aw.currsection = 0;

  aw.format = format;
  aw.mode = mode;
  aw.style = style;

  aw.hup = FALSE;
  aw.ssp = NULL;

  omdp = ObjMgrGetData (entityID);
  if (omdp != NULL && omdp->datatype == OBJ_SEQSUB) {
    ssp = (SeqSubmitPtr) omdp->dataptr;
    if (ssp != NULL && ssp->datatype == 1) {
      aw.ssp = ssp;
      sbp = ssp->sub;
      if (sbp != NULL) {
        aw.hup = sbp->hup;
      }
    }
  }

  if (bssp != NULL) {

    /* handle all components of a pop/phy/mut set */

    sep = SeqMgrGetSeqEntryForData (bssp);
    DoPopPhyMutSet (sep, &aw);

  } else {

    /* handle single bioseq, which may be segmented or a local part */

    DoOneBioseq (bsp, &aw);
  }

  /* check for failure to populate anything */

  numSections = ValNodeLen (aw.sectionList);
  ajp->numSections = numSections;

  if (numSections == 0) return ajp;

  /* allocate section array for this job */

  sectionArray = (Asn2gbSectionPtr PNTR) MemNew (sizeof (Asn2gbSectionPtr) * (numSections + 1));
  ajp->sectionArray = sectionArray;

  if (sectionArray == NULL) return ajp;

  /* fill in section and paragraph arrays */

  numParagraphs = 0;
  for (vnp = aw.sectionList, i = 0; vnp != NULL && i < numSections; vnp = vnp->next, i++) {
    asp = (Asn2gbSectionPtr) vnp->data.ptrvalue;
    sectionArray [i] = asp;
    if (asp != NULL) {
      numParagraphs += asp->numBlocks;
    }
  }

  /* allocate paragraph array pointing to all blocks in all sections */

  ajp->numParagraphs = numParagraphs;
  if (numParagraphs == 0) return ajp;

  paragraphArray = (BaseBlockPtr PNTR) MemNew (sizeof (BaseBlockPtr) * (numParagraphs + 1));
  ajp->paragraphArray = paragraphArray;

  paragraphByIDs = (BaseBlockPtr PNTR) MemNew (sizeof (BaseBlockPtr) * (numParagraphs + 1));
  ajp->paragraphByIDs = paragraphByIDs;

  if (paragraphArray == NULL || paragraphByIDs == NULL) return ajp;

  k = 0;
  for (i = 0; i < numSections; i++) {
    asp = sectionArray [i];
    if (asp != NULL) {

      numBlocks = asp->numBlocks;
      blockArray = asp->blockArray;
      if (blockArray != NULL) {

        for (j = 0; j < numBlocks; j++) {
          bbp = blockArray [j];

          paragraphArray [k] = bbp;
          paragraphByIDs [k] = bbp;
          bbp->paragraph = k;
          k++;
        }
      }
    }
  }

  /* sort paragraphByIDs array by entityID/itemtype/itemID/paragraph */

  HeapSort (paragraphByIDs, (size_t) numParagraphs, sizeof (BaseBlockPtr), SortParagraphByIDProc);

  /* free sectionList, but leave data, now pointed to by sectionArray elements */

  ValNodeFree (aw.sectionList);

  return ajp;
}

typedef CharPtr (*FormatProc) (Asn2gbFormatPtr afp, BaseBlockPtr bbp);

static FormatProc asn2gnbk_fmt_functions [22] = {
  NULL,
  DefaultFormatBlock, DefaultFormatBlock, DefaultFormatBlock,
  DefaultFormatBlock, DefaultFormatBlock, DefaultFormatBlock,
  DefaultFormatBlock, DefaultFormatBlock, DefaultFormatBlock,
  FormatSourceBlock, FormatOrganismBlock, FormatReferenceBlock,
  FormatCommentBlock, DefaultFormatBlock, FormatSourceFeatBlock,
  FormatFeatureBlock, FormatBasecountBlock, DefaultFormatBlock,
  FormatSequenceBlock, FormatContigBlock, DefaultFormatBlock
};

static void PrintFtableIntervals (
  ValNodePtr PNTR head,
  BioseqPtr target,
  SeqLocPtr location,
  CharPtr label
)

{
  SeqLocPtr  slp;
  Int4       start;
  Int4       stop;
  Char       str [64];

  if (head == NULL || target == NULL || location == NULL || label == NULL) return;

  slp = SeqLocFindNext (location, NULL);
  if (slp == NULL) return;

  start = GetOffsetInBioseq (slp, target, SEQLOC_START) + 1;
  stop = GetOffsetInBioseq (slp, target, SEQLOC_STOP) + 1;
  sprintf (str, "%ld\t%ld\t%s\n", (long) start, (long) stop, label);
  ValNodeCopyStr (head, 0, str);

  while ((slp = SeqLocFindNext (location, slp)) != NULL) {
    start = GetOffsetInBioseq (slp, target, SEQLOC_START) + 1;
    stop = GetOffsetInBioseq (slp, target, SEQLOC_STOP) + 1;
    if (start != 0 && stop != 0) {
      sprintf (str, "%ld\t%ld\n", (long) start, (long) stop);
      ValNodeCopyStr (head, 0, str);
    }
  }
}

static void PrintFtableLocAndQuals (
  ValNodePtr PNTR head,
  BioseqPtr target,
  SeqFeatPtr sfp,
  SeqMgrFeatContextPtr context
)

{
  DbtagPtr     dbt;
  GBQualPtr    gbq;
  GeneRefPtr   grp;
  CharPtr      label;
  ObjectIdPtr  oip;
  BioseqPtr    prod;
  SeqFeatPtr   prot;
  ProtRefPtr   prp;
  RnaRefPtr    rrp;
  SeqIdPtr     sip;
  Char         str [256];
  Char         tmp [300];
  tRNAPtr      trna;
  ValNodePtr   vnp;

  if (head == NULL || target == NULL || sfp == NULL || context == NULL) return;
  label = (CharPtr) FeatDefTypeLabel (sfp);
  if (StringCmp (label, "Gene") == 0) {
    label = "gene";
  }
  if (StringHasNoText (label)) {
    label = "???";
  }

  PrintFtableIntervals (head, target, sfp->location, label);
  switch (context->seqfeattype) {
    case SEQFEAT_GENE :
      grp = (GeneRefPtr) sfp->data.value.ptrvalue;
      if (grp != NULL) {
        StringNCpy_0 (str, (CharPtr) grp->locus, sizeof (str));
        if (! StringHasNoText (str)) {
          sprintf (tmp, "\t\t\tgene\t%s\n", str);
          ValNodeCopyStr (head, 0, tmp);
        }
        for (vnp = grp->syn; vnp != NULL; vnp = vnp->next) {
          StringNCpy_0 (str, (CharPtr) vnp->data.ptrvalue, sizeof (str));
          if (! StringHasNoText (str)) {
            sprintf (tmp, "\t\t\tgene_syn\t%s\n", str);
            ValNodeCopyStr (head, 0, tmp);
          }
        }
      }
      break;
    case SEQFEAT_CDREGION :
      prod = BioseqFind (SeqLocId (sfp->product));
      prot = SeqMgrGetBestProteinFeature (prod, NULL);
      if (prot != NULL) {
        prp = (ProtRefPtr) prot->data.value.ptrvalue;
        if (prp != NULL) {
          if (prp->name != NULL) {
            for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
              StringNCpy_0 (str, (CharPtr) vnp->data.ptrvalue, sizeof (str));
              if (! StringHasNoText (str)) {
                sprintf (tmp, "\t\t\tproduct\t%s\n", str);
                ValNodeCopyStr (head, 0, tmp);
              }
            }
          } else if (prp->desc != NULL) {
            StringNCpy_0 (str, prp->desc, sizeof (str));
            if (! StringHasNoText (str)) {
              sprintf (tmp, "\t\t\tproduct\t%s\n", str);
              ValNodeCopyStr (head, 0, tmp);
            }
          }
          for (vnp = prp->activity; vnp != NULL; vnp = vnp->next) {
            StringNCpy_0 (str, (CharPtr) vnp->data.ptrvalue, sizeof (str));
            if (! StringHasNoText (str)) {
              sprintf (tmp, "\t\t\tfunction\t%s\n", str);
              ValNodeCopyStr (head, 0, tmp);
            }
          }
          for (vnp = prp->ec; vnp != NULL; vnp = vnp->next) {
            StringNCpy_0 (str, (CharPtr) vnp->data.ptrvalue, sizeof (str));
            if (! StringHasNoText (str)) {
              sprintf (tmp, "\t\t\tEC_number\t%s\n", str);
              ValNodeCopyStr (head, 0, tmp);
            }
          }
        }
      }
      if (prod != NULL) {
        for (sip = prod->id; sip != NULL; sip = sip->next) {
          if (sip->choice == SEQID_GENBANK ||
              sip->choice == SEQID_EMBL ||
              sip->choice == SEQID_DDBJ) {
            if (SeqIdWrite (sip, str, PRINTID_TEXTID_ACC_VER, sizeof (str)) != NULL) {
              sprintf (tmp, "\t\t\tprotein_id\t%s\n", str);
              ValNodeCopyStr (head, 0, tmp);
            }
          }
        }
      }
      break;
    case SEQFEAT_RNA :
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp != NULL) {
        switch (rrp->ext.choice) {
          case 1 :
            StringNCpy_0 (str, (CharPtr) rrp->ext.value.ptrvalue, sizeof (str));
            if (! StringHasNoText (str)) {
              sprintf (tmp, "\t\t\tproduct\t%s\n", str);
              ValNodeCopyStr (head, 0, tmp);
            }
            break;
          case 2 :
            trna = rrp->ext.value.ptrvalue;
            if (trna != NULL) {
              FeatDefLabel (sfp, str, sizeof (str) - 1, OM_LABEL_CONTENT);
              if (! StringHasNoText (str)) {
                sprintf (tmp, "\t\t\tproduct\t%s\n", str);
                ValNodeCopyStr (head, 0, tmp);
              }
            }
            break;
          default :
            break;
        }
      }
      break;
    default :
      break;
  }
  if (! StringHasNoText (sfp->comment)) {
    ValNodeCopyStr (head, 0, "\t\t\tnote\t");
    ValNodeCopyStr (head, 0, sfp->comment);
    ValNodeCopyStr (head, 0, "\n");
  }
  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    if (! StringHasNoText (gbq->qual)) {
      if (! StringHasNoText (gbq->val)) {
        sprintf (tmp, "\t\t\t%s\t%s\n", gbq->qual, gbq->val);
        ValNodeCopyStr (head, 0, tmp);
      }
    }
  }
  for (vnp = sfp->dbxref; vnp != NULL; vnp = vnp->next) {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt != NULL) {
      if (! StringHasNoText (dbt->db)) {
        oip = dbt->tag;
        if (oip->str != NULL && (! StringHasNoText (oip->str))) {
          sprintf (tmp, "\t\t\tdb_xref\t%s:%s\n", dbt->db, oip->str);
          ValNodeCopyStr (head, 0, tmp);
        } else {
          sprintf (tmp, "\t\t\tdb_xref\t%s:%ld\n", dbt->db, (long) oip->id);
          ValNodeCopyStr (head, 0, tmp);
        }
      }
    }
  }
}

NLM_EXTERN CharPtr asn2gnbk_format (
  Asn2gbJobPtr ajp,
  Int4 paragraph
)

{
  Asn2gbFormat       af;
  Asn2gbSectionPtr   asp;
  BaseBlockPtr       bbp;
  BlockType          blocktype;
  SeqMgrFeatContext  fcontext;
  FormatProc         fmt;
  ValNodePtr         head;
  Char               id [42];
  QualVal            qv [ASN2GNBK_TOTAL_SOURCE + 5];
  Int4               section;
  SeqFeatPtr         sfp;
  SeqIdPtr           sip;
  CharPtr            str = NULL;
  BioseqPtr          target;
  Char               tmp [53];

  /* qv must hold MAX (ASN2GNBK_TOTAL_SOURCE, ASN2GNBK_TOTAL_FEATUR) */

  if (ajp == NULL || ajp->sectionArray == NULL || ajp->paragraphArray == NULL) return NULL;
  if (paragraph < 0 || paragraph >= ajp->numParagraphs) return NULL;

  bbp = ajp->paragraphArray [paragraph];
  if (bbp == NULL) return NULL;

  section = bbp->section;
  if (section < 0 || section >= ajp->numSections) return NULL;

  asp = ajp->sectionArray [section];
  if (asp == NULL) return NULL;

  blocktype = bbp->blocktype;
  if (blocktype < LOCUS_BLOCK || blocktype > SLASH_BLOCK) return NULL;

  MemSet ((Pointer) qv, 0, sizeof (qv));

  af.ajp = ajp;
  af.asp = asp;
  af.qvp = qv;
  af.format = ajp->format;
  af.mode = ajp->mode;
  af.style = ajp->style;

  if (ajp->format != FTABLE_FMT) {
    fmt = asn2gnbk_fmt_functions [(int) blocktype];
    if (fmt == NULL) return NULL;
    str = fmt (&af, bbp);

  } else {

    target = asp->target;
    if (target != NULL) {

      if (blocktype == FEATHEADER_BLOCK) {
        sip = SeqIdFindBest (target->id, 0);
        SeqIdWrite (sip, id, PRINTID_FASTA_LONG, sizeof (id) - 1);
        if (! StringHasNoText (id)) {
          sprintf (tmp, ">Feature %s\n", id);
          str = StringSave (tmp);
        }

      } else if (blocktype == FEATURE_BLOCK) {

        sfp = SeqMgrGetDesiredFeature (bbp->entityID, NULL, bbp->itemID, 0, NULL, &fcontext);
        if (sfp != NULL) {
          head = NULL;
          PrintFtableLocAndQuals (&head, target, sfp, &fcontext);
          str = MergeValNodeStrings (head);
          ValNodeFreeData (head);
        }
      }
    }
  }

  if (str == NULL) {
    str = StringSave ("???\n");
  }

  return str;
}

NLM_EXTERN Asn2gbJobPtr asn2gnbk_cleanup (
  Asn2gbJobPtr ajp
)

{
  Asn2gbSectionPtr   asp;
  BaseBlockPtr       bbp;
  BaseBlockPtr       PNTR blockArray;
  Int4               i;
  IntCdsBlockPtr     icp;
  IntFeatBlockPtr    ifp;
  IntRefBlockPtr     irp;
  Int4               j;
  Int4               numBlocks;
  Int4               numSections;
  ReferenceBlockPtr  rrp;
  Asn2gbSectionPtr   PNTR sectionArray;

  if (ajp == NULL) return NULL;

  SeqLocFree (ajp->slp);

  numSections = ajp->numSections;
  sectionArray = ajp->sectionArray;

  if (sectionArray != NULL) {

    for (i = 0; i < numSections; i++) {
      asp = sectionArray [i];
      if (asp != NULL) {

        numBlocks = asp->numBlocks;
        blockArray = asp->blockArray;
        if (blockArray != NULL) {

          for (j = 0; j < numBlocks; j++) {
            bbp = blockArray [j];
            if (bbp != NULL) {

              MemFree (bbp->string);

              if (bbp->blocktype == REFERENCE_BLOCK) {
                rrp = (ReferenceBlockPtr) bbp;
                MemFree (rrp->uniquestr);
                irp = (IntRefBlockPtr) rrp;
                DateFree (irp->date);
                SeqLocFree (irp->loc);
                MemFree (irp->authstr);
                MemFree (irp->fig);
                MemFree (irp->maploc);

              } else if (bbp->blocktype == FEATURE_BLOCK) {

                ifp = (IntFeatBlockPtr) bbp;
                if (ifp->isCDS) {
                  icp = (IntCdsBlockPtr) ifp;
                  MemFree (icp->fig);
                  MemFree (icp->maploc);
                }
              }

              MemFree (bbp);
            }
          }
        }
        MemFree (asp->blockArray);
        MemFree (asp->referenceArray);
        SeqPortFree (asp->spp);
      }
    }
  }

  MemFree (ajp->sectionArray);
  MemFree (ajp->paragraphArray);
  MemFree (ajp->paragraphByIDs);

  MemFree (ajp);
  return NULL;
}

static Boolean LIBCALLBACK LockAllSegments (
  SeqLocPtr slp,
  SeqMgrSegmentContextPtr context
)

{
  BioseqPtr        bsp;
  SeqLocPtr        loc;
  SeqIdPtr         sip;
  ValNodePtr PNTR  vnpp;

  if (slp == NULL || context == NULL) return FALSE;
  vnpp = (ValNodePtr PNTR) context->userdata;
  if (vnpp == NULL) return TRUE;

  sip = SeqLocId (slp);
  if (sip == NULL) {
    loc = SeqLocFindNext (slp, NULL);
    if (loc != NULL) {
      sip = SeqLocId (loc);
    }
  }
  if (sip == NULL) return TRUE;

  bsp = BioseqLockById (sip);
  ValNodeAddPointer (vnpp, 0, (Pointer) bsp);

  return TRUE;
}

static Boolean LIBCALLBACK LockAllBioseqs (
  BioseqPtr bsp,
  SeqMgrBioseqContextPtr context
)

{
  ValNodePtr PNTR  vnpp;

  if (bsp == NULL || context == NULL) return FALSE;
  vnpp = (ValNodePtr PNTR) context->userdata;
  if (vnpp == NULL) return TRUE;

  BioseqLock (bsp);
  ValNodeAddPointer (vnpp, 0, (Pointer) bsp);

  if (bsp->repr == Seq_repr_seg) {
    SeqMgrExploreSegments (bsp, (Pointer) vnpp, LockAllSegments);
  }

  return TRUE;
}

NLM_EXTERN Boolean SeqEntryToGnbk (
  SeqEntryPtr sep,
  FmtType format,
  ModType mode,
  StlType style,
  FILE *fp
)

{
  Asn2gbJobPtr  ajp;
  BioseqPtr     bsp = NULL;
  BioseqSetPtr  bssp = NULL;
  Uint2         entityID;
  Int4          i;
  Int4          numParagraphs;
  CharPtr       str;
  ValNodePtr    bsplist = NULL;
  ValNodePtr    vnp;

  if (sep == NULL || fp == NULL) return FALSE;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
  }

  entityID = ObjMgrGetEntityIDForPointer (sep->data.ptrvalue);
  if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
    SeqMgrIndexFeatures (entityID, NULL);
  }

  /* lock all bioseqs in advance, including remote genome components */

  SeqMgrExploreBioseqs (0, sep->data.ptrvalue, (Pointer) &bsplist, LockAllBioseqs, TRUE, FALSE, FALSE);

  ajp = asn2gnbk_setup (bsp, bssp, NULL, format, mode, style);

  if (ajp != NULL) {

    ajp->keepSeqPortOpen = TRUE;

    numParagraphs = ajp->numParagraphs;
    for (i = 0; i < numParagraphs; i++) {
      str = asn2gnbk_format (ajp, i);
      if (str != NULL) {
        fprintf (fp, "%s", str);
      } else {
        fprintf (fp, "?\n");
      }
      MemFree (str);
    }

    asn2gnbk_cleanup (ajp);
  }

  /* unlock all pre-locked bioseqs, including remote genome components */

  for (vnp = bsplist; vnp != NULL; vnp = vnp->next) {
    bsp = (BioseqPtr) vnp->data.ptrvalue;
    if (bsp != NULL) {
      BioseqUnlock (bsp);
    }
  }
  bsplist = ValNodeFree (bsplist);

  return TRUE;
}


/* ********************************************************************** */

#define TEST_MODULE_ASN2GNBK

#ifdef TEST_MODULE_ASN2GNBK

NLM_EXTERN Boolean SeqEntryToFlat (
  SeqEntryPtr sep,
  FILE *fp,
  Uint1 format,
  Uint1 mode
);

extern void SeriousSeqEntryCleanup (
  SeqEntryPtr sep,
  SeqEntryFunc taxfun,
  SeqEntryFunc taxmerge
);

/*
static void SummarizeQuals (
  void
)

{
  FILE  *fp;
  Int2  i;
  Int2  j;
  Int2  k;

  fp = FileOpen ("QualSummary", "w");
  if (fp == NULL) return;

  for (i = 0; i < ParFlat_TOTAL_GBQUAL; i++) {
    fprintf (fp, "%s\n", ParFlat_GBQual_names[i].name);
    for (j = 0; j < ParFlat_TOTAL_GBFEAT; j++) {
      for (k = 0; k < ParFlat_GBFeat [j].opt_num; k++) {
        if (ParFlat_GBFeat [j].opt_qual [k] == i) {
          fprintf (fp, "  opt %s\n", ParFlat_GBFeat [j].key);
        }
      }
      for (k = 0; k < ParFlat_GBFeat [j].mand_num; k++) {
        if (ParFlat_GBFeat [j].mand_qual [k] == i) {
          fprintf (fp, "  mnd %s\n", ParFlat_GBFeat [j].key);
        }
      }
    }
    fprintf (fp, "\n");
  }

  FileClose (fp);
}
*/

static Int2 HandleSingleRecord (
  CharPtr inputFile,
  CharPtr outputFile,
  FmtType format,
  ModType mode,
  StlType style
)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Pointer       dataptr;
  Uint2         datatype;
  Uint2         entityID;
  FILE          *fp;
  SeqEntryPtr   sep;

  fp = FileOpen (inputFile, "r");
  if (fp == NULL) {
    Message (MSG_FATAL, "FileOpen failed for input file '%s'", inputFile);
    return 1;
  }

  dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, FALSE, FALSE, FALSE, FALSE);

  FileClose (fp);

  entityID = ObjMgrRegister (datatype, dataptr);
  if (datatype == OBJ_SEQSUB || datatype == OBJ_SEQENTRY ||
        datatype == OBJ_BIOSEQ || datatype == OBJ_BIOSEQSET) {

#ifdef WIN_MAC
#if __profile__
    ProfilerSetStatus (TRUE);
#endif
#endif
    entityID = SeqMgrIndexFeatures (entityID, NULL);
#ifdef WIN_MAC
#if __profile__
    ProfilerSetStatus (FALSE);
#endif
#endif

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
      FileRemove (outputFile);
#ifdef WIN_MAC
      FileCreate (outputFile, "TEXT", "ttxt");
#endif
      fp = FileOpen (outputFile, "w");
      SeqEntryToGnbk (sep, format, mode, style, fp);
      FileClose (fp);
    }
  } else {
    Message (MSG_FATAL, "Datatype %d not recognized", (int) datatype);
  }
  ObjMgrFree (datatype, dataptr);

  return 0;
}

static void SaveSeqEntry (
  SeqEntryPtr sep,
  CharPtr filename
)

{
  AsnIoPtr  aop;

  if (sep == NULL) return;
  aop = AsnIoOpen (filename, "w");
  if (aop != NULL) {
    SeqEntryAsnWrite (sep, aop, NULL);
  }
  AsnIoClose (aop);
}

static void SaveAsn2gnbk (
  SeqEntryPtr sep,
  CharPtr filename,
  FmtType format,
  ModType mode,
  StlType style
)

{
  FILE  *fp;

  if (sep == NULL) return;
  fp = FileOpen (filename, "w");
  if (fp != NULL) {
    SeqEntryToGnbk (sep, format, mode, style, fp);
  }
  FileClose (fp);
}

static void SaveAsn2ff (
  SeqEntryPtr sep,
  CharPtr filename
)

{
  FILE  *fp;

  if (sep == NULL) return;
  fp = FileOpen (filename, "w");
  if (fp != NULL) {
    SeqEntryToFlat (sep, fp, 0, 8);
  }
  FileClose (fp);
}

typedef struct hasgidata {
  Int4     gi;
  Boolean  found;
} HasGiData, PNTR HasGiPtr;

static void LookForGi (
  SeqEntryPtr sep,
  Pointer mydata,
  Int4 index,
  Int2 indent
)

{
  BioseqPtr  bsp;
  HasGiPtr   hgp;
  SeqIdPtr   sip;

  if (sep == NULL) return;
  if (! IS_Bioseq (sep)) return;
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL) return;
  hgp = (HasGiPtr) mydata;
  if (hgp == NULL) return;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_GI) {
      if (sip->data.intvalue == hgp->gi) {
        hgp->found = TRUE;
        return;
      }
    }
  }
}

static Boolean SeqEntryHasGi (
  SeqEntryPtr sep,
  Int4 gi
)

{
  HasGiData  hgd;

  if (sep == NULL || gi < 1) return FALSE;
  hgd.gi = gi;
  hgd.found = FALSE;
  SeqEntryExplore (sep, (Pointer) (&hgd), LookForGi);
  return hgd.found;
}

static void CompareFlatFiles (
  CharPtr path1,
  CharPtr path2,
  CharPtr path3,
  SeqEntryPtr sep,
  FILE* fp,
  FmtType format,
  ModType mode,
  StlType style,
  Boolean batch,
  Boolean diff,
  Boolean gbdjoin
)

{
#ifdef OS_UNIX
  Char    buf [256];
  Char    cmmd [256];
  size_t  ct;
  FILE    *fpo;

  if (batch) {

    SeriousSeqEntryCleanup (sep, NULL, NULL);

    SaveAsn2gnbk (sep, path1, format, mode, style);

    SaveAsn2ff (sep, path2);

  } else if (diff) {

    SaveAsn2ff (sep, path1);

    SeriousSeqEntryCleanup (sep, NULL, NULL);

    SaveAsn2ff (sep, path2);

  }

  if (gbdjoin) {
    sprintf (cmmd, "/netopt/genbank/subtool/bin/gbdjoin -o %s -n %s -d reports", path1, path2);
    system (cmmd);

    sprintf (cmmd, "rm %s; rm %s", path1, path2);
    system (cmmd);
  } else {
    sprintf (cmmd, "sort %s | uniq -c > %s.suc; rm %s", path1, path1, path1);
    system (cmmd);

    sprintf (cmmd, "sort %s | uniq -c > %s.suc; rm %s", path2, path2, path2);
    system (cmmd);

    sprintf (cmmd, "diff %s.suc %s.suc > %s", path1, path2, path3);
    system (cmmd);

    sprintf (cmmd, "cat %s", path3);
    fpo = popen (cmmd, "r");
    if (fpo != NULL) {
      while ((ct = fread (buf, 1, sizeof (buf), fpo)) > 0) {
        fwrite (buf, 1, ct, fp);
        fflush (stdout);
      }
      pclose (fpo);
    }

    sprintf (cmmd, "rm %s.suc; rm %s.suc; rm %s", path1, path2, path3);
    system (cmmd);
  }

#else
  SeqEntryToGnbk (sep, format, mode, style, fp);
#endif
}

static Int2 HandleMultipleRecords (
  CharPtr inputFile,
  CharPtr outputFile,
  CharPtr iomode,
  FmtType format,
  ModType mode,
  StlType style,
  Int4 gi,
  Boolean batch,
  Boolean diff,
  Boolean gbdjoin
)

{
  AsnIoPtr      aip;
  AsnModulePtr  amp;
  AsnTypePtr    atp, atp_bss, atp_se;
  BioseqPtr     bsp;
  Char          buf [40];
  FILE          *fp;
  SeqEntryPtr   fsep;
  Boolean       hasgi;
  Char          path1 [PATH_MAX];
  Char          path2 [PATH_MAX];
  Char          path3 [PATH_MAX];
  SeqEntryPtr   sep;

  amp = AsnAllModPtr ();
  if (amp == NULL) {
    Message (MSG_FATAL, "Unable to load AsnAllModPtr");
    return 1;
  }

  atp_bss = AsnFind ("Bioseq-set");
  if (atp_bss == NULL) {
    Message (MSG_FATAL, "Unable to find ASN.1 type Bioseq-set");
    return 1;
  }

  atp_se = AsnFind ("Bioseq-set.seq-set.E");
  if (atp_se == NULL) {
    Message (MSG_FATAL, "Unable to find ASN.1 type Bioseq-set.seq-set.E");
    return 1;
  }

  aip = AsnIoOpen (inputFile, iomode);
  if (aip == NULL) {
    Message (MSG_FATAL, "AsnIoOpen failed for input file '%s'", inputFile);
    return 1;
  }

  fp = FileOpen (outputFile, "w");
  if (fp == NULL) {
    AsnIoClose (aip);
    Message (MSG_FATAL, "FileOpen failed for output file '%s'", outputFile);
    return 1;
  }

  TmpNam (path1);
  TmpNam (path2);
  TmpNam (path3);

  atp = atp_bss;

  while ((atp = AsnReadId (aip, amp, atp)) != NULL) {
    if (atp == atp_se) {
      sep = SeqEntryAsnRead (aip, atp);

      fsep = FindNthBioseq (sep, 1);
      if (fsep != NULL && fsep->choice == 1) {
        bsp = (BioseqPtr) fsep->data.ptrvalue;
        if (bsp != NULL) {
          SeqIdWrite (bsp->id, buf, PRINTID_FASTA_LONG, sizeof (buf));
          if (SeqEntryHasNucs (sep)) {
            fprintf (fp, "%s\n", buf);
          } else {
            fprintf (fp, "%s -- ignored\n", buf);
          }
#ifdef OS_UNIX
          printf ("%s\n", buf);
#endif
        }
      }

      hasgi = SeqEntryHasGi (sep, gi);
      if (hasgi) {
        sprintf (buf, "%ld.before", (long) gi);
        SaveSeqEntry (sep, buf);
        sprintf (buf, "%ld.gbff.before", (long) gi);
        SaveAsn2ff (sep, buf);
        SeriousSeqEntryCleanup (sep, NULL, NULL);
        sprintf (buf, "%ld.after", (long) gi);
        SaveSeqEntry (sep, buf);
        sprintf (buf, "%ld.gbff.after", (long) gi);
        SaveAsn2ff (sep, buf);
        FileClose (fp);
        AsnIoClose (aip);
        return 0;
      }
      if (gi == 0 && SeqEntryHasNucs (sep)) {
        CompareFlatFiles (path1, path2, path3, sep, fp,
                          format, mode, style,
                          batch, diff, gbdjoin);
      }
      SeqEntryFree (sep);
    } else {
      AsnReadVal (aip, atp, NULL);
    }
  }

  FileClose (fp);
  AsnIoClose (aip);

  return 0;
}

#define i_argInputFile 0
#define o_argOutputFile 1
#define f_argFormat 2
#define m_argMode 3
#define s_argStyle 4
#define t_argBatch 5
#define d_argDiffFF 6
#define j_argGbdjoin 7
#define b_argBinary 8
#define g_argGiToSave 9
#define r_argRemote 10

Args myargs [] = {
  {"Input File Name", "stdin", NULL, NULL,
    FALSE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
  {"Output File Name", "stdout", NULL, NULL,
    FALSE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
  {"Format (1 GenBank, 2 EMBL, 3 Feature Table)", "1", NULL, NULL,
    TRUE, 'f', ARG_INT, 0.0, 0, NULL},
  {"Mode (1 Release, 2 Sequin, 3 Dump)", "2", NULL, NULL,
    TRUE, 'm', ARG_INT, 0.0, 0, NULL},
  {"Style (1 Normal, 2 Master, 3 Contig)", "1", NULL, NULL,
    TRUE, 's', ARG_INT, 0.0, 0, NULL},
  {"Batch Compare Asn2gnbk to Asn2ff", "F", NULL, NULL,
    TRUE, 't', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Batch Diff SeriousSeqEntryCleanup with Asn2ff", "F", NULL, NULL,
    FALSE, 'd', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Batch Use Gbdjoin for Diffs", "F", NULL, NULL,
    FALSE, 'j', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Bioseq-set is Binary", "T", NULL, NULL,
    TRUE, 'b', ARG_BOOLEAN, 0.0, 0, NULL},
  {"GI to save", "0", NULL, NULL,
    TRUE, 'g', ARG_STRING, 0.0, 0, NULL},
  {"Remote fetching", "F", NULL, NULL,
    FALSE, 'r', ARG_BOOLEAN, 0.0, 0, NULL},
};


#ifdef OS_UNIX
#include <accid1.h>
#include <lsqfetch.h>
#endif

Int2 Main (
  void
)

{
  Boolean  batch = FALSE;
  Boolean  diff = FALSE;
  FmtType  format;
  Boolean  gbdjoin = FALSE;
  Int4     gi = 0;
  CharPtr  iomode = "r";
  ModType  mode;
  Char     path [PATH_MAX];
  CharPtr  progname;
  Int2     rsult = 0;
  StlType  style;
  long     val;

  ErrSetFatalLevel (SEV_MAX);
  ErrClearOptFlags (EO_SHOW_USERSTR);
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
  if (! SeqCodeSetLoad ()) {
    Message (MSG_FATAL, "SeqCodeSetLoad failed");
    return 1;
  }
  if (! GeneticCodeTableLoad ()) {
    Message (MSG_FATAL, "GeneticCodeTableLoad failed");
    return 1;
  }

  ProgramPath (path, sizeof (path));
  progname = StringRChr (path, DIRDELIMCHR);
  if (progname != NULL) {
    progname++;
  } else {
    progname = "asn2gnbk";
  }

  if (! GetArgs (progname, sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  if (myargs [t_argBatch].intvalue) {
    batch = TRUE;
  } else {
    batch = FALSE;
  }

  if (myargs [d_argDiffFF].intvalue) {
    diff = TRUE;
  } else {
    diff = FALSE;
  }

  if (myargs [j_argGbdjoin].intvalue) {
    gbdjoin = TRUE;
  } else {
    gbdjoin = FALSE;
  }

  if (gbdjoin) {
    if ((! batch) && (! diff)) {
      diff = TRUE;
    }
  }

  if (myargs [b_argBinary].intvalue) {
    iomode = "rb";
  } else {
    iomode = "r";
  }

  if (! StringHasNoText (myargs [g_argGiToSave].strvalue)) {
    if (sscanf (myargs [g_argGiToSave].strvalue, "%ld", &val) == 1) {
      gi = (Int4) val;
    }
  }

  switch (myargs [f_argFormat].intvalue) {
    case 1 :
      format = GENBANK_FMT;
      break;
    case 2 :
      format = EMBL_FMT;
      break;
    case 3 :
      format = FTABLE_FMT;
      break;
    default :
      format = GENBANK_FMT;
      break;
  }

  switch (myargs [m_argMode].intvalue) {
    case 1 :
      mode = RELEASE_MODE;
      break;
    case 2 :
      mode = SEQUIN_MODE;
      break;
    case 3 :
      mode = DUMP_MODE;
      break;
    default :
      mode = RELEASE_MODE;
      break;
  }

  switch (myargs [s_argStyle].intvalue) {
    case 1 :
      style = NORMAL_STYLE;
      break;
    case 2 :
      style = MASTER_STYLE;
      break;
    case 3 :
      style = CONTIG_STYLE;
      break;
    default :
      style = NORMAL_STYLE;
      break;
  }

#ifdef OS_UNIX
  if (myargs [r_argRemote].intvalue) {
    ID1BioseqFetchEnable ("asn2gnbk", FALSE);
    LocalSeqFetchInit (FALSE);
  }
#endif

  if (batch || diff || gbdjoin || gi != 0) {
    rsult = HandleMultipleRecords (myargs [i_argInputFile].strvalue,
                                   myargs [o_argOutputFile].strvalue,
                                   iomode, format, mode, style, gi,
                                   batch, diff, gbdjoin);
  } else {
    rsult = HandleSingleRecord (myargs [i_argInputFile].strvalue,
                                myargs [o_argOutputFile].strvalue,
                                format, mode, style);
  }

#ifdef OS_UNIX
  if (myargs [r_argRemote].intvalue) {
    LocalSeqFetchDisable ();
    ID1BioseqFetchDisable ();
  }
#endif
  return rsult;
}
#endif /* TEST_MODULE_ASN2GNBK */

