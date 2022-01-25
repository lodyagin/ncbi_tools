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
* $Revision: 6.18 $
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
  Int2          format;
  Int2          mode;
  Int2          segstyle;
  Int2          seqstyle;

  /* linked lists of paragraphs, sections, blocks */

  ValNodePtr    sectionList;
  ValNodePtr    blockList;    /* reset for each new section */

  /* most recent node of linked lists, for quickly adding next node */

  ValNodePtr    lastsection;
  ValNodePtr    lastblock;    /* reset for each new section */

  Int2          currsection;

  /* section fields needed for populating blocks */

  BioseqPtr     parent;
  BioseqPtr     bsp;
  SeqLocPtr     slp;
  Uint2         seg;
  Int4          numsegs;
  Int4          from;
  Int4          to;

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
  Uint1Ptr      uip;
} QualVal, PNTR QualValPtr;

/* structure passed to individual paragraph format functions */

typedef struct asn2gbformat {
  Asn2gbJobPtr      ajp;
  Asn2gbSectionPtr  asp;
  QualValPtr        qvp;

  Int2              format;
  Int2              mode;
} Asn2gbFormat, PNTR Asn2gbFormatPtr;


/* Seq-hist replacedBy is preformatted into string field, */
/* then comment descriptors, Map location:, and Region:, */
/* then comment features, finally HTGS */

typedef struct comment_block {
  ASN2GB_BASE_BLOCK
  Boolean           first;
} CommentBlock, PNTR CommentBlockPtr;

/* internal reference block has date field on top of ReferenceBlock fields */

typedef struct int_ref_block {
  ReferenceBlock    rb;
  DatePtr           date;  /* internal sorting use only */
} IntRefBlock, PNTR IntRefBlockPtr;


/* ********************************************************************** */

/* utility functions */

static CharPtr MergeValNodeStrings (ValNodePtr list)

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

static void AddValNodeString (ValNodePtr PNTR head, CharPtr prefix, CharPtr string, CharPtr suffix)

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

static Int2 gb_StartPrint (Int2 format, Boolean call_init_buff,
                           Int2 gb_init_indent, Int2 gb_cont_indent,
                           CharPtr gb_label, Int2 gb_tab_to,
                           Int2 eb_init_indent, Int2 eb_cont_indent,
                           CharPtr eb_line_prefix, Boolean eb_print_xx)

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

static void gb_AddString (CharPtr prefix, CharPtr string, CharPtr suffix,
                          Boolean addPeriod, Boolean convertQuotes, Boolean expandTildes)

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
      if (ch == ' ' || ch == '\t' || ch == '~') {
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

  /*
  if (convertQuotes) {
    strptr = newstr;
    ch = *strptr;
    while (ch != '\0') {
      if (ch == '\"') {
        *strptr = '\'';
      }
      strptr++;
      ch = *strptr;
    }
  }
  */

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

static CharPtr gb_MergeString (Boolean call_end_print)

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

static CharPtr DateToGB (CharPtr buf, DatePtr dp)

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

static CharPtr DefaultFormatBlock (Asn2gbFormatPtr afp, BaseBlockPtr bbp)

{
  if (afp == NULL || bbp == NULL) return NULL;

  /* default format function assumes string pre-allocated by add block function */

  return StringSaveNoNull (bbp->string);
}

static CharPtr FormatDeflineBlock (Asn2gbFormatPtr afp, BaseBlockPtr bbp)

{
  SeqMgrDescContext  dcontext;
  ValNodePtr         sdp;
  CharPtr            title = NULL;

  if (afp == NULL || bbp == NULL) return NULL;

  /* CreateDefLine results saved in bbp->string */

  if (! StringHasNoText (bbp->string)) return StringSave (bbp->string);

  /* otherwise should reference title descriptor IDs */

  if (bbp->itemtype == OBJ_SEQDESC) {
    sdp = SeqMgrGetDesiredDescriptor (bbp->entityID, NULL, bbp->itemID, 0, NULL, &dcontext);
    if (sdp != NULL && dcontext.seqdesctype == Seq_descr_title) {
      title = (CharPtr) sdp->data.ptrvalue;
    }
  }

  if (title == NULL) {
    title = "No definition line found";
  }

  gb_StartPrint (afp->format, TRUE, 0, 12, "DEFINITION", 13, 5, 5, "DE", TRUE);

  gb_AddString (NULL, title, NULL, TRUE, TRUE, FALSE);

  return gb_MergeString (TRUE);
}

static CharPtr FormatOrganismBlock (Asn2gbFormatPtr afp, BaseBlockPtr bbp)

{
  BioSourcePtr       biop = NULL;
  CharPtr            common = NULL;
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext;
  CharPtr            lineage = NULL;
  OrgNamePtr         onp;
  OrgRefPtr          orp;
  ValNodePtr         sdp;
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

    gb_StartPrint (afp->format, FALSE, 2, 12, "ORGANISM", 13, 5, 5, "OC", FALSE);
    gb_AddString (NULL, taxname, NULL, FALSE, FALSE, FALSE);
    ff_EndPrint();

    gb_StartPrint (afp->format, FALSE, 12, 12, NULL, 0, 5, 5, "OC", FALSE);
    gb_AddString (NULL, lineage, NULL, TRUE, FALSE, FALSE);
    ff_EndPrint();

  } else if (afp->format == EMBL_FMT || afp->format == EMBLPEPT_FMT) {

    gb_StartPrint (afp->format, TRUE, 0, 12, "SOURCE", 13, 5, 5, "OS", TRUE);
    gb_AddString (NULL, taxname, NULL, FALSE, FALSE, FALSE);
    gb_AddString (" (", common, ")", FALSE, FALSE, FALSE);
    ff_EndPrint();

    gb_StartPrint (afp->format, FALSE, 12, 12, NULL, 0, 5, 5, "OC", FALSE);
    gb_AddString (NULL, lineage, NULL, TRUE, FALSE, FALSE);
    ff_EndPrint();
  }

  return gb_MergeString (FALSE);
}

/* format references section */

static AuthListPtr GetAuthListPtr (PubdescPtr pdp, CitSubPtr csp)

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

static CharPtr MakeSingleAuthorString (Int2 format, CharPtr prefix, CharPtr name, CharPtr initials, CharPtr suffix)

{
  Char     ch;
  size_t   len;
  CharPtr  ptr;
  CharPtr  str;
  CharPtr  tmp;

  if (name == NULL) return NULL;
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

static CharPtr GetAuthorsString (Int2 format, AuthListPtr alp)

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

static void StrStripSpaces(CharPtr str)

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

static Boolean AllCaps (CharPtr p)

{
  if (p == NULL) return FALSE;

  for (p++; p != NULL && *p != '\0'; p++) {
    if (IS_LOWER (*p)) return FALSE;
  }
  return TRUE;
}

static void CleanEquals (CharPtr p)

{
  if (p == NULL) return;

  for (; *p != '\0'; p++) {
    if (*p == '\"') {
      *p = '\'';
    }
  }
}

static CharPtr GetPubTitle (Int2 format, PubdescPtr pdp, CitSubPtr csp)

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

static void CleanPubTitle (CharPtr title)

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

static Int2 FixPages (CharPtr out_pages, CharPtr in_pages)

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

static void DoSup (ValNodePtr PNTR head, CharPtr issue, CharPtr part_sup, CharPtr part_supi)

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

static CharPtr FormatCitJour (Int2 format, CitJourPtr cjp)

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

  DoSup (&head, issue, part_sup, part_supi);

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

static CharPtr FormatCitArt (Int2 format, CitArtPtr cap)

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
      break;
    case 3 :
      cbp = (CitBookPtr) cap->fromptr;
      break;
    default :
      break;
  }

  return rsult;
}

static CharPtr GetAffil (AffilPtr afp)

{
	Boolean need_comma=FALSE;
	CharPtr string=NULL, temp, ptr;
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
	return string;
}

static CharPtr FormatCitSub (Int2 format, CitSubPtr csp)

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

static CharPtr GetPubJournal (Int2 format, PubdescPtr pdp, CitSubPtr csp)

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
          journal = StringSave ("CitGen");
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
      case PUB_Man :
        cbp = (CitBookPtr) vnp->data.ptrvalue;
        if (cbp != NULL) {
          journal = StringSave ("CitBook");
        }
        break;
      case PUB_Patent :
        cpp = (CitPatPtr) vnp->data.ptrvalue;
        if (cpp != NULL) {
          journal = StringSave ("CitPat");
        }
        break;
      default :
        break;
    }

    if (journal != NULL) return journal;
  }

  return NULL;
}

static Int4 GetMuid (PubdescPtr pdp)

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

static CharPtr FormatReferenceBlock (Asn2gbFormatPtr afp, BaseBlockPtr bbp)

{
  Asn2gbJobPtr       ajp;
  AuthListPtr        alp;
  Asn2gbSectionPtr   asp;
  BioseqPtr          bsp;
  Char               buf [32];
  CitSubPtr          csp = NULL;
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext;
  SeqLocPtr          loc = NULL;
  Int4               muid;
  Boolean            needsPeriod;
  SeqLocPtr          nextslp;
  ObjMgrDataPtr      omdp;
  PubdescPtr         pdp = NULL;
  CharPtr            prefix;
  ReferenceBlockPtr  rbp;
  SubmitBlockPtr     sbp;
  ValNodePtr         sdp;
  SeqFeatPtr         sfp = NULL;
  SeqLocPtr          slp;
  SeqSubmitPtr       ssp;
  Int4               start;
  Int4               stop;
  CharPtr            str;
  CharPtr            suffix = NULL;
  Boolean            trailingPeriod;

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

  if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
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

  if (rbp->category == REF_CAT_SIT) {

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

    if (sfp != NULL) {
      loc = sfp->location;
    } else if (ajp->slp != NULL) {
      loc = ajp->slp;
    }

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

  gb_StartPrint (afp->format, FALSE, 2, 12, "TITLE", 13, 5, 5, "RT", FALSE);

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

  gb_AddString (prefix, str, suffix, FALSE, FALSE, FALSE);

  MemFree (str);
  ff_EndPrint ();

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

  /* !!! remainder of fields are only for GenBank !!! */

  if (afp->format == GENBANK_FMT || afp->format == GENPEPT_FMT) {
  }

  return gb_MergeString (FALSE);
}

static CharPtr FormatCommentBlock (Asn2gbFormatPtr afp, BaseBlockPtr bbp)

{
  CommentBlockPtr    cbp;
  CharPtr            db;
  DbtagPtr           dbt;
  SeqMgrDescContext  dcontext;
  ObjectIdPtr        oip;
  CharPtr            prefix;
  ValNodePtr         sdp;
  Char               sfx [32];
  CharPtr            suffix;
  CharPtr            title;

  if (afp == NULL || bbp == NULL) return NULL;
  cbp = (CommentBlockPtr) bbp;

  /* some comments are allocated (along with possible first COMMENT label) */

  if (! StringHasNoText (bbp->string)) return StringSave (bbp->string);

  /* !!! also have to deal with comment feature across entire sequence !!! */

  /* otherwise should reference comment, maploc, or region descriptor IDs */

  title = NULL;
  prefix = NULL;
  suffix = NULL;
  sfx [0] = '\0';

  if (bbp->itemtype == OBJ_SEQDESC) {
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

static Boolean is_real_id (SeqIdPtr sip, SeqIdPtr this_sip)

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

static Boolean FlatVirtLoc (BioseqPtr bsp, SeqLocPtr location)

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

static void FlatLocSeqId (ValNodePtr PNTR head, SeqIdPtr sip)

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
    SeqIdWrite (use_id, buf, PRINTID_TEXTID_ACCESSION, sizeof (buf) - 1);
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

static void FlatLocCaret (ValNodePtr PNTR head, SeqIdPtr sip, SeqIdPtr this_sip, Int4 point, IntFuzzPtr fuzz)

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

static void FlatLocPoint (ValNodePtr PNTR head, SeqIdPtr sip, SeqIdPtr this_sip,
                          Int4 point, IntFuzzPtr fuzz)

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

static void FlatLocElement (ValNodePtr PNTR head, BioseqPtr bsp, SeqLocPtr location)

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

static Boolean FlatNullAhead (BioseqPtr bsp, ValNodePtr location)

{
  SeqLocPtr  next;

  if (bsp == NULL || location == NULL) return FALSE;

  next = location->next;
  if (next == NULL) return TRUE;
  if (next->choice == SEQLOC_NULL) return TRUE;
  if (FlatVirtLoc (bsp, next)) return TRUE;

  return FALSE;
}

static void FlatPackedPoint (ValNodePtr PNTR head, PackSeqPntPtr pspp, BioseqPtr bsp)

{
  Uint1  dex;

  if (head == NULL || pspp == NULL || bsp == NULL) return;

  for (dex = 0; dex < pspp->used; dex++) {
    FlatLocPoint (head, pspp->id, bsp->id, pspp->pnts [dex], pspp->fuzz);
  }
}

static void DoFlatLoc (ValNodePtr PNTR head, BioseqPtr bsp,
                       SeqLocPtr location, Boolean ok_to_complement);

static void GroupFlatLoc (ValNodePtr PNTR head, BioseqPtr bsp, SeqLocPtr location,
                          CharPtr prefix, Boolean is_flat_order)

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

static void DoFlatLoc (ValNodePtr PNTR head, BioseqPtr bsp,
                       SeqLocPtr location, Boolean ok_to_complement)

{
  Boolean        found_null;
  SeqLocPtr      next_loc;
  PackSeqPntPtr  pspp;
  SeqLocPtr      slp;

  if (head == NULL || bsp == NULL || location == NULL) return;

  /* deal with complement of entire location */

  if (ok_to_complement && SeqLocStrand (location) == Seq_strand_minus) {
    slp = AsnIoMemCopy (location,
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

static CharPtr FlatLoc (BioseqPtr bsp, SeqLocPtr location)

{
  ValNodePtr  head = NULL;
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

  DoFlatLoc (&head, bsp, location, TRUE);

  str = MergeValNodeStrings (head);
  ValNodeFreeData (head);

  return str;
}

static CharPtr featurekeys [] = {
  "???" ,
  "gene" ,
  "Org" ,
  "CDS" ,
  "Protein" ,
  "preRNA" ,
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
  "Region" ,
  "misc_feature" ,
  "Bond" ,
  "Site" ,
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

/* enumerated qualifier category definitions */

typedef enum {
  Qual_class_ignore = 0,
  Qual_class_string,
  Qual_class_boolean,
  Qual_class_boolstr,
  Qual_class_int,
  Qual_class_valnode,
  Qual_class_gbqual,
  Qual_class_L_R_B,
  Qual_class_rpt,
  Qual_class_orgmod,
  Qual_class_subsource,
  Qual_class_code_break,
  Qual_class_anti_codon,
  Qual_class_codon,
  Qual_class_pubset,
  Qual_class_db_xref,
  Qual_class_seq_id,
  Qual_class_its,
  Qual_class_trna_codons,
  Qual_class_translation,
  Qual_class_illegal_qual,
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
  SOURCE_chloroplast,
  SOURCE_chromoplast,
  SOURCE_chromosome,
  SOURCE_citation,
  SOURCE_clone,
  SOURCE_clone_lib,
  SOURCE_common,
  SOURCE_common_name,
  SOURCE_country,
  SOURCE_cultivar,
  SOURCE_cyanelle,
  SOURCE_db_xref,
  SOURCE_dev_stage,
  SOURCE_dosage,
  SOURCE_extrachrom,
  SOURCE_focus,
  SOURCE_frequency,
  SOURCE_genomic,
  SOURCE_genotype,
  SOURCE_germline,
  SOURCE_group,
  SOURCE_haplotype,
  SOURCE_ins_seq_name,
  SOURCE_insertion_seq,
  SOURCE_isolate,
  SOURCE_kinetoplast,
  SOURCE_lab_host,
  SOURCE_label,
  SOURCE_macronuclear,
  SOURCE_map,
  SOURCE_mitochondrion,
  SOURCE_note,
  SOURCE_old_name,
  SOURCE_organism,
  SOURCE_orgmod_note,
  SOURCE_pathovar,
  SOURCE_plasmid,
  SOURCE_plasmid_name,
  SOURCE_plastid,
  SOURCE_plastid_name,
  SOURCE_pop_variant,
  SOURCE_proviral,
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
  SOURCE_transposon,
  SOURCE_transposon_name,
  SOURCE_type,
  SOURCE_usedin,
  SOURCE_variety,
  SOURCE_virion,
  ASN2GNBK_TOTAL_SOURCE
}  SourceType;

/* ordering arrays for qualifiers and note components */

static Uint1 relmode_source_qual_order [] = {
  SOURCE_organism,
  SOURCE_focus,

  SOURCE_genomic,
  SOURCE_chloroplast,
  SOURCE_chromoplast,
  SOURCE_kinetoplast,
  SOURCE_mitochondrion,
  SOURCE_extrachrom,
  SOURCE_macronuclear,
  SOURCE_cyanelle,
  SOURCE_proviral,
  SOURCE_virion,
  SOURCE_plasmid,
  SOURCE_plasmid_name,
  SOURCE_plastid,
  SOURCE_plastid_name,
  SOURCE_ins_seq_name,
  SOURCE_insertion_seq,
  SOURCE_transposon,
  SOURCE_transposon_name,

  SOURCE_germline,
  SOURCE_rearranged,

  SOURCE_acronym,
  SOURCE_biotype,
  SOURCE_biovar,
  SOURCE_cell_line,
  SOURCE_cell_type,
  SOURCE_chemovar,
  SOURCE_chromosome,
  SOURCE_clone,
  SOURCE_clone_lib,
  SOURCE_common,
  SOURCE_common_name,
  SOURCE_country,
  SOURCE_cultivar,
  SOURCE_dev_stage,
  SOURCE_dosage,
  SOURCE_frequency,
  SOURCE_genotype,
  SOURCE_group,
  SOURCE_haplotype,
  SOURCE_isolate,
  SOURCE_lab_host,
  SOURCE_map,
  SOURCE_old_name,
  SOURCE_pathovar,
  SOURCE_pop_variant,
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
  SOURCE_tissue_lib,
  SOURCE_tissue_type,
  SOURCE_type,
  SOURCE_variety,

  SOURCE_sequenced_mol, SOURCE_label, SOURCE_usedin,
  SOURCE_db_xref, SOURCE_citation, SOURCE_note, 0
};

static Uint1 relmode_source_note_order [] = {
  SOURCE_seqfeat_note, SOURCE_subsource_note, SOURCE_orgmod_note, 0
};

static Uint1 seqmode_source_qual_order [] = {
  SOURCE_organism,
  SOURCE_focus,

  SOURCE_genomic,
  SOURCE_chloroplast,
  SOURCE_chromoplast,
  SOURCE_kinetoplast,
  SOURCE_mitochondrion,
  SOURCE_extrachrom,
  SOURCE_macronuclear,
  SOURCE_cyanelle,
  SOURCE_proviral,
  SOURCE_virion,
  SOURCE_plasmid,
  SOURCE_plasmid_name,
  SOURCE_plastid,
  SOURCE_plastid_name,
  SOURCE_ins_seq_name,
  SOURCE_insertion_seq,
  SOURCE_transposon,
  SOURCE_transposon_name,

  SOURCE_germline,
  SOURCE_rearranged,

  SOURCE_acronym,
  SOURCE_biotype,
  SOURCE_biovar,
  SOURCE_cell_line,
  SOURCE_cell_type,
  SOURCE_chemovar,
  SOURCE_chromosome,
  SOURCE_clone,
  SOURCE_clone_lib,
  SOURCE_common,
  SOURCE_common_name,
  SOURCE_country,
  SOURCE_cultivar,
  SOURCE_dev_stage,
  SOURCE_dosage,
  SOURCE_frequency,
  SOURCE_genotype,
  SOURCE_group,
  SOURCE_haplotype,
  SOURCE_isolate,
  SOURCE_lab_host,
  SOURCE_map,
  SOURCE_old_name,
  SOURCE_pathovar,
  SOURCE_pop_variant,
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
  SOURCE_tissue_lib,
  SOURCE_tissue_type,
  SOURCE_type,
  SOURCE_variety,

  SOURCE_sequenced_mol, SOURCE_label, SOURCE_usedin,
  SOURCE_db_xref, SOURCE_citation, SOURCE_note, 0
};

static Uint1 seqmode_source_note_order [] = {
  SOURCE_seqfeat_note, SOURCE_subsource_note, SOURCE_orgmod_note, 0
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
  { "chloroplast",      Qual_class_boolean   },
  { "chromoplast",      Qual_class_boolean   },
  { "chromosome",       Qual_class_subsource },
  { "citation",         Qual_class_pubset    },
  { "clone",            Qual_class_subsource },
  { "clone_lib",        Qual_class_subsource },
  { "common",           Qual_class_orgmod    },
  { "common",           Qual_class_string    },
  { "country",          Qual_class_subsource },
  { "cultivar",         Qual_class_orgmod    },
  { "cyanelle",         Qual_class_boolean   },
  { "db_xref",          Qual_class_db_xref   },
  { "dev_stage",        Qual_class_subsource },
  { "dosage",           Qual_class_orgmod    },
  { "extrachromosomal", Qual_class_boolean   },
  { "focus",            Qual_class_boolean   },
  { "frequency",        Qual_class_subsource },
  { "genomic",          Qual_class_ignore    },
  { "genotype",         Qual_class_subsource },
  { "germline",         Qual_class_boolean   },
  { "group",            Qual_class_orgmod    },
  { "haplotype",        Qual_class_subsource },
  { "insertion_seq",    Qual_class_subsource },
  { "insertion_seq",    Qual_class_boolstr   },
  { "isolate",          Qual_class_orgmod    },
  { "kinetoplast",      Qual_class_boolean   },
  { "lab_host",         Qual_class_subsource },
  { "label",            Qual_class_gbqual    },
  { "macronuclear",     Qual_class_boolean   },
  { "map",              Qual_class_subsource },
  { "mitochondrion",    Qual_class_boolean   },
  { "note",             Qual_class_note      },
  { "old_name",         Qual_class_orgmod    },
  { "organism",         Qual_class_string    },
  { "orgmod_note",      Qual_class_orgmod    },
  { "pathovar",         Qual_class_orgmod    },
  { "plasmid",          Qual_class_boolstr   },
  { "plasmid",          Qual_class_subsource },
  { "plastid",          Qual_class_boolstr   },
  { "plastid",          Qual_class_subsource },
  { "pop_variant",      Qual_class_subsource },
  { "proviral",         Qual_class_boolean   },
  { "rearranged",       Qual_class_boolean   },
  { "seqfeat_note",     Qual_class_string    },
  { "sequenced_mol",    Qual_class_gbqual    },
  { "serogroup",        Qual_class_orgmod    },
  { "serotype",         Qual_class_orgmod    },
  { "serovar",          Qual_class_orgmod    },
  { "sex",              Qual_class_subsource },
  { "specific_host",    Qual_class_orgmod    },
  { "specimen_voucher", Qual_class_orgmod    },
  { "strain",           Qual_class_orgmod    },
  { "sub_clone",        Qual_class_subsource },
  { "sub_group",        Qual_class_orgmod    },
  { "sub_species",      Qual_class_orgmod    },
  { "sub_strain",       Qual_class_orgmod    },
  { "sub_type",         Qual_class_orgmod    },
  { "subsource_note",   Qual_class_subsource },
  { "tissue_lib",       Qual_class_subsource },
  { "tissue_type",      Qual_class_subsource },
  { "transposon",       Qual_class_boolstr   },
  { "transposon",       Qual_class_subsource },
  { "type",             Qual_class_orgmod    },
  { "usedin",           Qual_class_gbqual    },
  { "variety",          Qual_class_orgmod    },
  { "virion",           Qual_class_boolean   }
};

static Int2 genomeToSourceIdx [15] = {
  0,
  SOURCE_genomic,
  SOURCE_chloroplast,
  SOURCE_chromoplast,
  SOURCE_kinetoplast,
  SOURCE_mitochondrion,
  SOURCE_plastid,
  SOURCE_macronuclear,
  SOURCE_extrachrom,
  SOURCE_plasmid,
  SOURCE_transposon,
  SOURCE_insertion_seq,
  SOURCE_cyanelle,
  SOURCE_proviral,
  SOURCE_virion
};

static void GenomeToQualArray (Uint1 genome, QualValPtr qvp)

{
  Int2  idx;

  if (qvp == NULL) return;

  if (genome < 15) {
    idx = genomeToSourceIdx [genome];
    if (idx > 0 && idx < ASN2GNBK_TOTAL_SOURCE) {
      qvp [idx].bool = TRUE;
    }
  }
}

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

static void SubSourceToQualArray (SubSourcePtr ssp, QualValPtr qvp)

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

static void OrgModToQualArray (OrgModPtr omp, QualValPtr qvp)

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

static CharPtr FormatSourceBlock (Asn2gbFormatPtr afp, BaseBlockPtr bbp)

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
  Uint1              lasttype;
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
  ValNodePtr         sdp;
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
  str = FlatLoc (bsp, location);
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

  GenomeToQualArray (biop->genome, qvp);

  SubSourceToQualArray (biop->subtype, qvp);

  if (orp != NULL) {
    onp = orp->orgname;
    if (onp != NULL) {
      OrgModToQualArray (onp->mod, qvp);
    }

    qvp [SOURCE_db_xref].vnp = orp->db;
  }

  /* some qualifiers are flags in genome and names in subsource, print once with name */

  if (qvp [SOURCE_ins_seq_name].ssp != NULL) {
    qvp [SOURCE_insertion_seq].bool = FALSE;
  }
  if (qvp [SOURCE_plasmid_name].ssp != NULL) {
    qvp [SOURCE_plasmid].bool = FALSE;
  }
  if (qvp [SOURCE_plastid_name].ssp != NULL) {
    qvp [SOURCE_plastid].bool = FALSE;
  }
  if (qvp [SOURCE_transposon_name].ssp != NULL) {
    qvp [SOURCE_transposon].bool = FALSE;
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

    lasttype = 0;
    switch (asn2gnbk_source_quals [idx].qualclass) {
      case Qual_class_ignore :
        break;
      case Qual_class_string :
        if (! StringHasNoText (qvp [idx].str)) {
          sprintf (buf, "/%s=\"", asn2gnbk_source_quals [idx].name);
          NewContLine ();
          gb_AddString (buf, qvp [idx].str, "\"", FALSE, TRUE, FALSE);
        }
        break;
      case Qual_class_boolean :
        if (qvp [idx].bool) {
          sprintf (buf, "/%s", asn2gnbk_source_quals [idx].name);
          NewContLine ();
          ff_AddString (buf);
        }
        break;
      case Qual_class_boolstr :
        if (qvp [idx].bool) {
          sprintf (buf, "/%s=\"\"", asn2gnbk_source_quals [idx].name);
          NewContLine ();
          ff_AddString (buf);
        }
        break;
      case Qual_class_orgmod :
        omp = qvp [idx].omp;
        if (lasttype == 0 && omp != NULL) {
          lasttype = omp->subtype;
        }
        while (omp != NULL && omp->subtype == lasttype) {
          if (! StringHasNoText (omp->subname)) {
            sprintf (buf, "/%s=\"", asn2gnbk_source_quals [idx].name);
            NewContLine ();
            gb_AddString (buf, omp->subname, "\"", FALSE, TRUE, FALSE);
          }
          omp = omp->next;
        }
        break;
      case Qual_class_subsource :
        ssp = qvp [idx].ssp;
        if (lasttype == 0 && ssp != NULL) {
          lasttype = ssp->subtype;
        }
        while (ssp != NULL && ssp->subtype == lasttype) {
          if (! StringHasNoText (ssp->name)) {
            sprintf (buf, "/%s=\"", asn2gnbk_source_quals [idx].name);
            NewContLine ();
            gb_AddString (buf, ssp->name, "\"", FALSE, TRUE, FALSE);
          }
          ssp = ssp->next;
        }
        break;
      case Qual_class_pubset :
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
      case Qual_class_illegal_qual :
        break;
      case Qual_class_note :
        head = NULL;
        notestr = NULL;
        prefix = NULL;
        for (j = 0, jdx = notetbl [j]; jdx != 0; j++, jdx = notetbl [j]) {
          switch (asn2gnbk_source_quals [jdx].qualclass) {
            case Qual_class_string :
              if (! StringHasNoText (qvp [jdx].str)) {
                AddValNodeString (&head, prefix, qvp [jdx].str, NULL);
                prefix = "; ";
              }
              break;
            case Qual_class_orgmod :
              omp = qvp [jdx].omp;
              if (omp != NULL) {
                if (! StringHasNoText (omp->subname)) {
                  AddValNodeString (&head, prefix, omp->subname, NULL);
                  prefix = "; ";
                }
              }
              break;
            case Qual_class_subsource :
              ssp = qvp [jdx].ssp;
              if (ssp != NULL) {
                if (! StringHasNoText (ssp->name)) {
                  AddValNodeString (&head, prefix, ssp->name, NULL);
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

  return gb_MergeString (TRUE);
}

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
  FEATUR_frequency,
  FEATUR_function,
  FEATUR_gene,
  FEATUR_gene_desc,
  FEATUR_gene_map,
  FEATUR_gene_syn,
  FEATUR_illegal_qual,
  FEATUR_label,
  FEATUR_map,
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
  FEATUR_prot_desc,
  FEATUR_prot_names,
  FEATUR_protein_id,
  FEATUR_pseudo,
  FEATUR_replace,
  FEATUR_rpt_family,
  FEATUR_rpt_type,
  FEATUR_rpt_unit,
  FEATUR_rrna_its,
  FEATUR_seqfeat_note,
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

  FEATUR_product,

  FEATUR_allele,
  FEATUR_anticodon,
  FEATUR_bound_moiety,
  FEATUR_clone,
  FEATUR_codon,
  FEATUR_codon_start,
  FEATUR_cons_splice,
  FEATUR_direction,
  FEATUR_EC_number,
  FEATUR_evidence,
  FEATUR_exception,
  FEATUR_frequency,
  FEATUR_function,
  FEATUR_gene_map,
  FEATUR_label,
  FEATUR_map,
  FEATUR_mod_base,
  FEATUR_number,
  FEATUR_organism,
  FEATUR_PCR_conditions,
  FEATUR_phenotype,
  FEATUR_prot_activity,
  FEATUR_protein_id,
  FEATUR_pseudo,
  FEATUR_replace,
  FEATUR_rpt_family,
  FEATUR_rpt_type,
  FEATUR_rpt_unit,
  FEATUR_standard_name,
  FEATUR_transl_except,
  FEATUR_transl_table,
  FEATUR_usedin,

  FEATUR_db_xref, FEATUR_citation, FEATUR_note,
  FEATUR_cds_product, FEATUR_translation,
  0
};

static Uint1 relmode_feat_note_order [] = {
  FEATUR_gene_syn,
  FEATUR_gene_desc,
  FEATUR_prot_comment,
  FEATUR_prot_names,
  FEATUR_prot_desc,
  FEATUR_rrna_its,
  FEATUR_trna_codons,
  FEATUR_seqfeat_note,
  0
};

static Uint1 seqmode_feat_qual_order [] = {
  FEATUR_partial,
  FEATUR_gene,

  FEATUR_product,

  FEATUR_allele,
  FEATUR_anticodon,
  FEATUR_bound_moiety,
  FEATUR_clone,
  FEATUR_codon,
  FEATUR_codon_start,
  FEATUR_cons_splice,
  FEATUR_direction,
  FEATUR_EC_number,
  FEATUR_evidence,
  FEATUR_exception,
  FEATUR_frequency,
  FEATUR_function,
  FEATUR_gene_map,
  FEATUR_label,
  FEATUR_map,
  FEATUR_mod_base,
  FEATUR_number,
  FEATUR_organism,
  FEATUR_PCR_conditions,
  FEATUR_phenotype,
  FEATUR_prot_activity,
  FEATUR_protein_id,
  FEATUR_pseudo,
  FEATUR_replace,
  FEATUR_rpt_family,
  FEATUR_rpt_type,
  FEATUR_rpt_unit,
  FEATUR_standard_name,
  FEATUR_transl_except,
  FEATUR_transl_table,
  FEATUR_usedin,

  FEATUR_db_xref, FEATUR_citation, FEATUR_note,
  FEATUR_cds_product, FEATUR_translation,
  0
};

static Uint1 seqmode_feat_note_order [] = {
  FEATUR_gene_syn,
  FEATUR_gene_desc,
  FEATUR_prot_comment,
  FEATUR_prot_names,
  FEATUR_prot_desc,
  FEATUR_trna_codons,
  FEATUR_rrna_its,
  FEATUR_seqfeat_note,
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
  { "bound_moiety",   Qual_class_gbqual       },
  { "product",        Qual_class_string       },
  { "citation",       Qual_class_pubset       },
  { "clone",          Qual_class_gbqual       },
  { "codon",          Qual_class_codon        },
  { "codon_start",    Qual_class_int          },
  { "cons_splice",    Qual_class_gbqual       },
  { "db_xref",        Qual_class_db_xref      },
  { "direction",      Qual_class_L_R_B        },
  { "EC_number",      Qual_class_valnode      },
  { "evidence",       Qual_class_int          },
  { "exception",      Qual_class_string       },
  { "frequency",      Qual_class_gbqual       },
  { "function",       Qual_class_gbqual       },
  { "gene",           Qual_class_string       },
  { "gene_desc",      Qual_class_string       },
  { "gene_map",       Qual_class_string       },
  { "gene_syn",       Qual_class_valnode      },
  { "illegal",        Qual_class_illegal_qual },
  { "label",          Qual_class_gbqual       },
  { "map",            Qual_class_gbqual       },
  { "mod_base",       Qual_class_gbqual       },
  { "note",           Qual_class_note         },
  { "number",         Qual_class_gbqual       },
  { "organism",       Qual_class_gbqual       },
  { "partial",        Qual_class_boolean      },
  { "PCR_conditions", Qual_class_gbqual       },
  { "phenotype",      Qual_class_gbqual       },
  { "product",        Qual_class_string       },
  { "product",        Qual_class_gbqual       },
  { "function",       Qual_class_valnode      },
  { "prot_comment",   Qual_class_string       },
  { "prot_desc",      Qual_class_string       },
  { "prot_names",     Qual_class_valnode      },
  { "protein_id",     Qual_class_seq_id       },
  { "pseudo",         Qual_class_boolean      },  
  { "replace",        Qual_class_gbqual       },
  { "rpt_family",     Qual_class_gbqual       },
  { "rpt_type",       Qual_class_rpt          },
  { "rpt_unit",       Qual_class_gbqual       },
  { "rrna_its",       Qual_class_its          },
  { "seqfeat_note",   Qual_class_string       },
  { "standard_name",  Qual_class_gbqual       },
  { "transl_except",  Qual_class_code_break   },
  { "transl_table",   Qual_class_int          },
  { "translation",    Qual_class_translation  },
  { "trna_aa",        Qual_class_ignore       },
  { "trna_codons",    Qual_class_trna_codons  },
  { "usedin",         Qual_class_gbqual       }
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

static Int2 GbqualToFeaturIndex (CharPtr qualname)

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

static CharPtr FormatFeatureBlock (Asn2gbFormatPtr afp, BaseBlockPtr bbp)

{
  Uint1              aa;
  Asn2gbSectionPtr   asp;
  Boolean            at_end = FALSE;
  BioseqPtr          bsp;
  Char               buf [80];
  Uint1              codon[4];
  CdRegionPtr        crp;
  DbtagPtr           dbt;
  SeqMgrFeatContext  fcontext;
  Uint1              featdeftype;
  Uint1              from;
  GBQualPtr          gbq;
  SeqMgrFeatContext  gcontext;
  ValNodePtr         gcp;
  SeqFeatPtr         gene;
  GeneRefPtr         grp;
  ValNodePtr         head;
  Int2               i;
  Uint1              idx;
  Int2               j;
  Uint1              jdx;
  Int2               k;
  CharPtr            key;
  CharPtr            lasttype;
  Int4               len;
  CharPtr            notestr;
  Uint1Ptr           notetbl = NULL;
  Char               numbuf [32];
  ObjectIdPtr        oip;
  Uint2              partial;
  SeqMgrFeatContext  pcontext;
  CharPtr            prefix;
  BioseqPtr          prod = NULL;
  SeqFeatPtr         prot;
  CharPtr            protein_seq = NULL;
  ProtRefPtr         prp;
  Boolean            pseudo = FALSE;
  Uint1Ptr           qualtbl = NULL;
  QualValPtr         qvp;
  Uint1              residue;
  RnaRefPtr          rrp;
  SeqFeatPtr         sfp;
  Uint1              shift;
  SeqIntPtr          sintp;
  SeqIdPtr           sip;
  SeqLocPtr          slp;
  SeqMapTablePtr     smtp;
  SeqPortPtr         spp;
  CharPtr            str;
  tRNAPtr            trna;
  Uint1Ptr           uip;
  ValNodePtr         vnp;
  CharPtr            wwwbuf;

  if (afp == NULL || bbp == NULL) return NULL;
  asp = afp->asp;
  if (asp == NULL) return NULL;
  bsp = (asp->bsp);
  if (bsp == NULL) return NULL;
  qvp = afp->qvp;
  if (qvp == NULL) return NULL;

  /* all features in this list are known to be valid for the designated mode */

  sfp = SeqMgrGetDesiredFeature (bbp->entityID, NULL, bbp->itemID, 0, NULL, &fcontext);
  if (sfp == NULL) return NULL;

  featdeftype = fcontext.featdeftype;
  if (featdeftype < FEATDEF_GENE || featdeftype >= FEATDEF_MAX) {
    featdeftype = FEATDEF_BAD;
  }
  key = featurekeys [featdeftype];

  gb_StartPrint (afp->format, TRUE, 5, 21, NULL, 0, 5, 21, "FT", FALSE);
  ff_AddString (key);
  TabToColumn (22);

  str = FlatLoc (bsp, sfp->location);
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
    partial = SeqLocPartialCheck (sfp->location);
    if (partial == SLP_COMPLETE || partial > SLP_OTHER) {
      qvp [FEATUR_partial].bool = TRUE;
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
      gene = SeqMgrGetOverlappingGene (sfp->location, &gcontext);
      if (gene != NULL) {
        grp = (GeneRefPtr) gene->data.value.ptrvalue;
        if (gene->pseudo) {
          pseudo = TRUE;
        }
      }
    }
    if (grp != NULL && grp->pseudo) {
      pseudo = TRUE;
    }
    if (grp != NULL && (! SeqMgrGeneIsSuppressed (grp))) {
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
        }

        if (! pseudo) {
          prp = SeqMgrGetProtXref (sfp);
          if (prp == NULL) {
            sip = SeqLocId (sfp->product);
            qvp [FEATUR_protein_id].sip = sip;
            prod = BioseqFind (sip);
            prot = SeqMgrGetBestProteinFeature (prod, &pcontext);
            if (prot != NULL) {
              prp = (ProtRefPtr) prot->data.value.ptrvalue;
              qvp [FEATUR_prot_comment].str = prot->comment;
            }
          }

          if (prp != NULL) {
            vnp = prp->name;
            if (vnp != NULL && (! StringHasNoText ((CharPtr) vnp->data.ptrvalue))) {
              qvp [FEATUR_cds_product].str = (CharPtr) vnp->data.ptrvalue;
              vnp = vnp->next;
              qvp [FEATUR_prot_names].vnp = vnp;
              qvp [FEATUR_prot_desc].str = prp->desc;
            } else if (! StringHasNoText (prp->desc)) {
              qvp [FEATUR_cds_product].str = prp->desc;
            }
            qvp [FEATUR_prot_activity].vnp = prp->activity;
            qvp [FEATUR_EC_number].vnp = prp->ec;
          }

          qvp [FEATUR_translation].bool = TRUE;
        }
        break;
      case SEQFEAT_PROT :
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
                qvp [FEATUR_trna_codons].uip = trna->codon;
              }
            }
          } else {
            if (rrp->ext.choice == 1) {
              str = (CharPtr) rrp->ext.value.ptrvalue;
              qvp [FEATUR_product].str = str;
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
            }
          }
        }
        break;
      case SEQFEAT_REGION :
        break;
      case SEQFEAT_COMMENT :
        break;
      default :
        break;
    }
  }

  /* common fields set here */

  qvp [FEATUR_pseudo].bool = pseudo;

  qvp [FEATUR_exception].str = sfp->except_text;
  if (sfp->excpt && qvp [FEATUR_exception].str == NULL) {
    qvp [FEATUR_exception].str = "No explanation supplied";
  }

  qvp [FEATUR_seqfeat_note].str = sfp->comment;
  qvp [FEATUR_db_xref].vnp = sfp->dbxref;
  qvp [FEATUR_citation].vnp = sfp->cit;

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
    }
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
      case Qual_class_valnode :
        break;
      case Qual_class_gbqual :
        gbq = qvp [idx].gbq;
        if (lasttype == NULL && gbq != NULL) {
          lasttype = gbq->qual;
        }
        while (gbq != NULL && StringICmp (gbq->qual, lasttype) == 0) {
          if (! StringHasNoText (gbq->val)) {
            sprintf (buf, "/%s=\"", asn2gnbk_featur_quals [idx].name);
            NewContLine ();
            gb_AddString (buf, gbq->val, "\"", FALSE, TRUE, FALSE);
          }
          gbq = gbq->next;
        }
        break;

      /* Qual_class_L_R_B AND Qual_class_rpt are gbquals with different output */

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
      case Qual_class_rpt :
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
      case Qual_class_seq_id :
        break;
      case Qual_class_trna_codons :
        break;
      case Qual_class_translation :
        if (qvp [idx].bool && prod != NULL) {
          len = SeqLocLen (sfp->product);
          if (len > 0) {
            if (SeqLocStart (sfp->location) == 0 || SeqLocStop (sfp->location) == bsp->length - 1) {
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
        }
        break;
      case Qual_class_illegal_qual :
        break;
      case Qual_class_note :
        head = NULL;
        notestr = NULL;
        prefix = NULL;
        for (j = 0, jdx = notetbl [j]; jdx != 0; j++, jdx = notetbl [j]) {
          switch (asn2gnbk_featur_quals [jdx].qualclass) {
            case Qual_class_string :
              if (! StringHasNoText (qvp [jdx].str)) {
                AddValNodeString (&head, prefix, qvp [jdx].str, NULL);
                prefix = "; ";
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
              uip = qvp [jdx].uip;
              if (uip != NULL && uip [0] != 255) {
                numbuf [0] = '\0';
                str = &(numbuf [0]);
                codon [3] = '\0';
                for (k = 0; k < 6; k++) {
                  if (uip [k] == 255) break;
                  if (CodonForIndex (uip [k], Seq_code_iupacna, codon)) {
                    StringCpy (str, (CharPtr) codon);
                    str += 3;
                  } else {
                    *str = '?';
                    str++;
                  }
                  if (k < 5 && uip [k + 1] != 255) {
                    StringCpy (str, ", ");
                    str += 2;
                  }
                }
                if (uip [1] == 255) {
                  AddValNodeString (&head, prefix, "codon recognized: ", numbuf);
                } else {
                  AddValNodeString (&head, prefix, "codons recognized: ", numbuf);
                }
                prefix = "; ";
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

  return gb_MergeString (TRUE);

  /*
  if (fcontext.seqfeattype == SEQFEAT_IMP) {

    for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
      gbqual_index = GBQualNameValid (gbq->qual);
      if (gbqual_index == -1) continue;
      NewContLine ();
      ff_AddChar ('/');
      ff_AddString (gbq->qual);
      class_qual = ParFlat_GBQual_names [gbqual_index].gbclass;
      if (class_qual == Class_none && CheckForEqualSign (gbq->qual) == 1) continue;
      ff_AddChar('=');
      if (class_qual == Class_text && StringCmp (gbq->val, "\"\"") == 0) {
        ff_AddString  (gbq->val);
        continue;
      }
      buf = NULL;
      if (get_www () && (class_qual == Class_text || class_qual == Class_note)) {
        buf = www_featloc (gbq->val);
      } else {
        buf = StringSave (gbq->val);
      }
      if (class_qual == Class_text || class_qual == Class_none ||
          class_qual == Class_ecnum || class_qual == Class_note) {
        ff_AddString ("\"");
      }
      if (class_qual == Class_note) {
        www_note_gi (buf);
      } else if (class_qual != Class_none) {
        if (StringCmp (gbq->qual, "transl_table") == 0) {
          www_gcode (buf);
        } else if (StringCmp (gbq->qual, "db_xref") == 0) {
          www_db_xref (buf);
        } else {
          ff_AddString (buf);
        }
      }
      if (class_qual == Class_text || class_qual == Class_none ||
          class_qual == Class_ecnum || class_qual == Class_note) {
        ff_AddString ("\"");
      }
      if (buf != NULL) {
        buf = MemFree (buf);
      }
    }

    if (! StringHasNoText (sfp->comment)) {
      NewContLine ();
      gb_AddString ("/note=\"", sfp->comment, "\"", FALSE, TRUE, FALSE);
    }

  } else {

    NewContLine ();
    gb_AddString (NULL, fcontext.label, NULL, FALSE, FALSE, FALSE);
  }

  return gb_MergeString (TRUE);
  */
}

static CharPtr FormatBasecountBlock (Asn2gbFormatPtr afp, BaseBlockPtr bbp)

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
  if (bsp->repr == Seq_repr_delta) {
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
      if (ctr < 1) {
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

static void PrintSeqLine (Int2 format, CharPtr buf, Int4 start, Int4 stop)

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

static CharPtr FormatSequenceBlock (Asn2gbFormatPtr afp, BaseBlockPtr bbp)

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
    if (bsp->repr == Seq_repr_delta) {
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

    if (! IS_residue (residue)) {
      if (residue == INVALID_RESIDUE) {
        if (is_na) {
          residue = 'N';
        } else {
          residue = 'X';
        }
      } else {
        residue = '?';
      }
    }

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
    i++;
    if (i >= ctr) {
      i = 0;
      ctr = (Int2) MIN ((Int4) (stop - pos), (Int4) sizeof (bases));
      ctr = SeqPortRead (spp, bases, ctr);
      if (ctr < 1) {
        bases [0] = SEQPORT_EOF;
      }
    }
    residue = (Uint1) bases [i];
  }

  buf [count] = '\0';
  if (count > 0) {
    PrintSeqLine (afp->format, buf, start, start + cnt);
  }

  return gb_MergeString (FALSE);
}

/* not yet implemented */

static CharPtr FormatContigBlock (Asn2gbFormatPtr afp, BaseBlockPtr bbp)

{
  SequenceBlockPtr  sbp;

  if (afp == NULL || bbp == NULL) return NULL;
  sbp = (SequenceBlockPtr) bbp;

  return StringSaveNoNull (bbp->string);
}


/* ********************************************************************** */

/* functions to record sections or blocks in linked lists */

static Asn2gbSectionPtr Asn2gbAddSection (Asn2gbWorkPtr awp)

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

static BaseBlockPtr Asn2gbAddBlock (Asn2gbWorkPtr awp, Int2 blocktype, size_t size)

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

static void AddLocusBlock (Asn2gbWorkPtr awp)

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
  ValNodePtr         sdp;
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

  SeqIdWrite (bsp->id, locus, PRINTID_TEXTID_LOCUS, sizeof (locus) - 1);

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

static void AddDeflineBlock (Asn2gbWorkPtr awp)

{
  BaseBlockPtr       bbp;
  BioseqPtr          bsp;
  CharPtr            buf;
  size_t             buflen = 1001;
  SeqMgrDescContext  dcontext;
  ItemInfo           ii;
  MolInfoPtr         mip;
  ValNodePtr         sdp;
  Uint1              tech;
  CharPtr            title;

  if (awp == NULL) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;

  bbp = Asn2gbAddBlock (awp, DEFLINE_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  title = NULL;
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_title, &dcontext);
  if (sdp != NULL) {
    bbp->entityID = dcontext.entityID;
    bbp->itemID = dcontext.itemID;
    bbp->itemtype = OBJ_SEQDESC;
    title = (CharPtr) sdp->data.ptrvalue;
  }

  /* if found title descriptor, just save IDs to it */

  if (! StringHasNoText (title)) return;

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

static Int2 ValidateAccession (CharPtr accession)

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

static void AddAccessionBlock (Asn2gbWorkPtr awp)

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
  ValNodePtr         sdp;
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

static void AddVersionBlock (Asn2gbWorkPtr awp)

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

  if (gi < 1 && accn == NULL) return;

  bbp = Asn2gbAddBlock (awp, VERSION_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  if (gi > 0) {
    sprintf (version, "g%ld", (long) gi);

    gb_StartPrint (awp->format, needInitBuff, 0, 12, "NID", 13, 5, 5, "NI", TRUE);
    needInitBuff = FALSE;

    gb_AddString (NULL, version, NULL, FALSE, FALSE, FALSE);

	ff_EndPrint();
	needEndPrint = FALSE;
  }

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
  }

  bbp->string = gb_MergeString (needEndPrint);
}

static void AddPidBlock (Asn2gbWorkPtr awp)

{
  BaseBlockPtr  bbp;

  if (awp == NULL) return;

  bbp = Asn2gbAddBlock (awp, PID_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  bbp->string = StringSave ("PID\n");
}

static void AddDbsourceBlock (Asn2gbWorkPtr awp)

{
  BaseBlockPtr  bbp;

  if (awp == NULL) return;

  bbp = Asn2gbAddBlock (awp, DBSOURCE_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  bbp->string = StringSave ("DBSOURCE\n");
}

static void AddDateBlock (Asn2gbWorkPtr awp)

{
  BaseBlockPtr       bbp;
  BioseqPtr          bsp;
  Char               date [40];
  SeqMgrDescContext  dcontext;
  DatePtr            dp;
  ValNodePtr         sdp;

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

static Int2 MatchArrayString (CharPtr array_string [], Int2 totalstr, CharPtr text)

{
  Int2 i;

  for (i = 0; i < totalstr && text != NULL; i++) {
    if (StringCmp (array_string [i], text) == 0) {
      return (i);
    }
  }

  return (-1);
}

static Boolean CheckSpecialKeyword (Boolean is_est, Boolean is_sts, Boolean is_gss, CharPtr kwd)

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

static void AddKeywordsBlock (Asn2gbWorkPtr awp)

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
  ValNodePtr         sdp;
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
        if (head != NULL) {
          ValNodeCopyStr (&head, 0, "; ");
        }
        ValNodeCopyStr (&head, 0, kwd);
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

static void AddSegmentBlock (Asn2gbWorkPtr awp)

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

static void AddOrganismBlock (Asn2gbWorkPtr awp)

{
  BaseBlockPtr       bbp;
  BioseqPtr          bsp;
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext;
  ValNodePtr         sdp;
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

static ReferenceBlockPtr AddPub (Asn2gbWorkPtr awp, ValNodePtr PNTR head, PubdescPtr pdp)

{
  Char               buf [121];
  CitArtPtr          cap;
  CitBookPtr         cbp;
  CitGenPtr          cgp;
  CitJourPtr         cjp;
  CitPatPtr          cpp;
  CitSubPtr          csp;
  DatePtr            dp = NULL;
  ImprintPtr         imp = NULL;
  IntRefBlockPtr     irp;
  ReferenceBlockPtr  rbp;
  ValNodePtr         vnp;

  if (awp == NULL || head == NULL || pdp == NULL) return NULL;

  rbp = (ReferenceBlockPtr) MemNew (sizeof (IntRefBlock));
  if (rbp == NULL) return NULL;
  rbp->blocktype = REFERENCE_BLOCK;
  rbp->section = awp->currsection;

  ValNodeAddPointer (head, 0, rbp);

  rbp->serial = INT2_MAX;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    switch (vnp->choice) {
      case PUB_Gen :
        /* may be unpublished, or may be serial number of swiss-prot reference */
        rbp->category = REF_CAT_UNP;
        cgp = (CitGenPtr) vnp->data.ptrvalue;
        if (cgp != NULL) {
          dp = cgp->date;
          if (cgp->serial_number > 0) {
            rbp->serial = cgp->serial_number;
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
        cpp = (CitPatPtr) vnp->data.ptrvalue;
        if (cpp != NULL) {
          dp = (DatePtr) cpp->date_issue;
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
    rbp->category = REF_CAT_SIT;
  }

  if (rbp->muid == 0 && rbp->pmid == 0) {
	if (PubLabelUnique (pdp->pub, buf, sizeof (buf) - 1, OM_LABEL_CONTENT, TRUE) > 0) {
	  rbp->uniquestr = StringSaveNoNull (buf);
	}
  }

  irp = (IntRefBlockPtr) rbp;
  irp->date = dp;

  return rbp;
}

static int LIBCALLBACK SortReferences (VoidPtr ptr1, VoidPtr ptr2)

{
  int                compare;
  IntRefBlockPtr     irp1;
  IntRefBlockPtr     irp2;
  ReferenceBlockPtr  rbp1;
  ReferenceBlockPtr  rbp2;
  Int2               status;
  ValNodePtr         vnp1;
  ValNodePtr         vnp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      rbp1 = (ReferenceBlockPtr) vnp1->data.ptrvalue;
      rbp2 = (ReferenceBlockPtr) vnp2->data.ptrvalue;
      if (rbp1 != NULL && rbp2 != NULL) {
        if (rbp1->serial > rbp2->serial) {
          return 1;
        } else if (rbp1->serial < rbp2->serial) {
          return -1;
        }
        if (rbp1->category > rbp2->category) {
          return 1;
        } else if (rbp1->category < rbp2->category) {
          return -1;
        }
        irp1 = (IntRefBlockPtr) rbp1;
        irp2 = (IntRefBlockPtr) rbp2;
        status = DateMatch (irp1->date, irp2->date, FALSE);
        if (status == 1 || status == -1) return status;
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
        if (rbp1->uniquestr != NULL && rbp2->uniquestr != NULL) {
          compare = StringICmp (rbp1->uniquestr, rbp2->uniquestr);
          if (compare > 0) {
            return 1;
          } else if (compare < 0) {
            return -1;
          }
        }
      }
    }
  }
  return 0;
}

static void AddReferenceBlock (Asn2gbWorkPtr awp, Asn2gbSectionPtr asp)

{
  Asn2gbJobPtr       ajp;
  BioseqPtr          bsp;
  SeqMgrDescContext  dcontext;
  Boolean            excise;
  SeqMgrFeatContext  fcontext;
  ValNodePtr         head = NULL;
  Int2               i;
  Int2               idx;
  IntRefBlockPtr     irp;
  Int4Ptr            ivals;
  ReferenceBlockPtr  lastrbp;
  ValNodePtr         next;
  Int2               numivals;
  Int2               numReferences;
  PubdescPtr         pdp;
  ValNodePtr         PNTR prev;
  ReferenceBlockPtr  rbp;
  ReferenceBlockPtr  PNTR referenceArray;
  ValNodePtr         sdp;
  SeqFeatPtr         sfp;
  Int4               start;
  Int4               stop;
  ValNodePtr         vnp;

  if (awp == NULL || asp == NULL) return;
  ajp = awp->ajp;
  if (ajp == NULL) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_pub, &dcontext);
  while (sdp != NULL) {
    pdp = (PubdescPtr) sdp->data.ptrvalue;
    rbp = AddPub (awp, &head, pdp);
    if (rbp != NULL) {

      rbp->entityID = dcontext.entityID;
      rbp->itemID = dcontext.itemID;
      rbp->itemtype = OBJ_SEQDESC;
    }
    sdp = SeqMgrGetNextDescriptor (bsp, sdp, Seq_descr_pub, &dcontext);
  }

  sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_PUB, 0, &fcontext);
  while (sfp != NULL) {
    ivals = fcontext.ivals;
    numivals = fcontext.numivals;
    if (ivals != NULL && numivals > 0) {

      idx = (numivals - 1) * 2;
      start = ivals [idx];
      stop = ivals [idx + 1];
      if (stop >= awp->from && stop <= awp->to) {

        /*
        start = ivals [0] + 1;
        stop = ivals [idx + 1] + 1;
        */
        pdp = (PubdescPtr) sfp->data.value.ptrvalue;
        rbp = AddPub (awp, &head, pdp);
        if (rbp != NULL) {

          rbp->entityID = fcontext.entityID;
          rbp->itemID = fcontext.itemID;
          rbp->itemtype = OBJ_SEQFEAT;
        }
      }
    }

    sfp = SeqMgrGetNextFeature (bsp, sfp, SEQFEAT_PUB, 0, &fcontext);
  }

  if (head == NULL) return;

  /* sort by existing serial, then pub/unpub/sites/sub, then date */

  head = SortValNode (head, SortReferences);

  if (awp->ssp != NULL) {

    /* add seq-submit citation to end of list */

    rbp = (ReferenceBlockPtr) MemNew (sizeof (IntRefBlock));
    if (rbp != NULL) {

      rbp->blocktype = REFERENCE_BLOCK;
      rbp->section = awp->currsection;
      rbp->serial = INT2_MAX;
      rbp->category = REF_CAT_SUB;

      rbp->entityID = ajp->entityID;
      rbp->itemID = 1;
      rbp->itemtype = OBJ_SEQSUB_CIT;

      ValNodeAddPointer (&head, 0, rbp);
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
      if (rbp != NULL) {
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
            excise = TRUE;
          }
        }
      }
    } else {
      lastrbp = rbp;
    }
    if (excise) {
      *prev = vnp->next;
      vnp->next = NULL;
      MemFree (rbp->uniquestr);
      MemFree (rbp);
      ValNodeFree (vnp);
    } else {
      prev = &(vnp->next);
    }
    vnp = next;
  }

  /* assign serial numbers, null out (temporary, hidden) date field */

  for (vnp = head, i = 1; vnp != NULL; vnp = vnp->next, i++) {
    rbp = (ReferenceBlockPtr) vnp->data.ptrvalue;
    if (rbp != NULL) {
      rbp->serial = i;
      irp = (IntRefBlockPtr) rbp;
      irp->date = NULL;
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

static CharPtr PrintDate (DatePtr dp)

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

static void AddHistCommentString (CharPtr prefix, CharPtr suffix, DatePtr dp, SeqIdPtr ids)

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

static void AddHTGSCommentString (BioseqPtr bsp, MolInfoPtr mip)

{
  CharPtr      buf;
  Char         buffer [256];
  Int2         buflen = 0;
  DeltaSeqPtr  dsp;
  ValNodePtr   head = NULL;
  Int4         num_s = 0;
  Int4         num_g = 0;
  CharPtr      str;

  if (bsp == NULL || bsp->repr != Seq_repr_delta || mip == NULL || mip->tech < 2) return;

  for (dsp = (DeltaSeqPtr) bsp->seq_ext, buflen = 0; dsp != NULL; dsp = dsp->next) {
    buflen += 80;
  }
  buf = MemNew (buflen + 1);
  if (buf == NULL) return;

  CountGapsInDeltaSeq (bsp, &num_s, &num_g, NULL, NULL, buf, buflen);

  if (mip->tech == MI_TECH_htgs_0) {

    sprintf (buffer, "* WARNING: This record contains %ld individual~", (long) (num_s - num_g));
    ValNodeCopyStr (&head, 0, buffer);
    ValNodeCopyStr (&head, 0, "* sequencing reads that have not been assembled into~");
    ValNodeCopyStr (&head, 0, "* contigs. Runs of N are used to separate the reads~");
    ValNodeCopyStr (&head, 0, "* and the order in which they appear is completely~");
    ValNodeCopyStr (&head, 0, "* arbitrary. Low-pass sequence sampling is useful for  ~");
    ValNodeCopyStr (&head, 0, "* identifying clones that may be gene-rich and allows  ~");
    ValNodeCopyStr (&head, 0, "* overlap relationships among clones to be deduced. ~");
    ValNodeCopyStr (&head, 0, "* However, it should not be assumed that this clone ~");
    ValNodeCopyStr (&head, 0, "* will be sequenced to completion. In the event that~");
    ValNodeCopyStr (&head, 0, "* the record is updated, the accession number will ~");
    ValNodeCopyStr (&head, 0, "* be preserved.~");
    ValNodeCopyStr (&head, 0, buf);

  } else if (mip->tech == MI_TECH_htgs_1) {

    ValNodeCopyStr (&head, 0, "* NOTE: This is a \"working draft\" sequence that currently~");
    sprintf (buffer, "* consists of %ld contigs. The true order of the pieces~", (long) (num_s - num_g));
    ValNodeCopyStr (&head, 0, buffer);
    ValNodeCopyStr (&head, 0, "* is not known and their order in this sequence record is~");
    ValNodeCopyStr (&head, 0, "* contigs. Runs of N are used to separate the reads~");
    ValNodeCopyStr (&head, 0, "* is not known and their order in this sequence record is~");
    ValNodeCopyStr (&head, 0, "* arbitrary. Gaps between the contigs are represented as~");
    ValNodeCopyStr (&head, 0, "* runs of N, but the exact sizes of the gaps are unknown.~");
    ValNodeCopyStr (&head, 0, "* This record will be updated with the finished sequence~");
    ValNodeCopyStr (&head, 0, "* as soon as it is available and the accession number will~");
    ValNodeCopyStr (&head, 0, "* be preserved.~");
    ValNodeCopyStr (&head, 0, buf);

  } else if (mip->tech == MI_TECH_htgs_2) {

    ValNodeCopyStr (&head, 0, "* NOTE: This is a \"working draft\" sequence that currently~");
    sprintf (buffer, "* consists of %ld contigs. Gaps between the contigs~", (long) (num_s - num_g));
    ValNodeCopyStr (&head, 0, buffer);
    ValNodeCopyStr (&head, 0, "* are represented as runs of N. The order of the pieces~");
    ValNodeCopyStr (&head, 0, "* is believed to be correct as given, however the sizes~");
    ValNodeCopyStr (&head, 0, "* of the gaps between them are based on estimates that have~");
    ValNodeCopyStr (&head, 0, "* provided by the submittor. This sequence will be replaced~");
    ValNodeCopyStr (&head, 0, "* by the finished sequence as soon as it is available and~");
    ValNodeCopyStr (&head, 0, "* the accession number will be preserved.~");
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

static CharPtr GetStrForUserObject (UserObjectPtr uop)

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

static void AddCommentBlock (Asn2gbWorkPtr awp)

{
  BioseqPtr          bsp;
  Char               buf [128];
  CommentBlockPtr    cbp;
  DbtagPtr           dbt;
  SeqMgrDescContext  dcontext;
  Boolean            first = TRUE;
  SeqHistPtr         hist;
  MolInfoPtr         mip;
  ValNodePtr         sdp;
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

          sprintf (buf, "GSDB:S:%ld.", (long) dbt->tag->id);

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
    if (mip != NULL && mip->tech >= MI_TECH_htgs_1 && mip->tech <= MI_TECH_fli_cdna) {

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

  /* !!! also need to add feature comments that are full length !!! */
}

static void AddFeatHeaderBlock (Asn2gbWorkPtr awp)

{
  BaseBlockPtr  bbp;

  if (awp == NULL) return;

  bbp = Asn2gbAddBlock (awp, FEATHEADER_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

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

static void AddSourceBlock (Asn2gbWorkPtr awp)

{
  BaseBlockPtr       bbp;
  BioseqPtr          bsp;
  SeqMgrDescContext  dcontext;
  SeqMgrFeatContext  fcontext;
  ValNodePtr         sdp;
  SeqFeatPtr         sfp;

  if (awp == NULL) return;
  bsp = awp->bsp;
  if (bsp == NULL) return;

  bbp = Asn2gbAddBlock (awp, SOURCE_BLOCK, sizeof (BaseBlock));
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

static Boolean LIBCALLBACK DoOneFeature (SeqFeatPtr sfp, SeqMgrFeatContextPtr fcontext)

{
  Asn2gbWorkPtr      awp;
  BaseBlockPtr       bbp;
  BioseqPtr          bsp;
  Int2               idx;
  Int4Ptr            ivals;
  Int2               numivals;
  SeqMgrFeatContext  pcontext;
  SeqFeatPtr         prt;
  Int4               start;
  Int4               stop;

  if (sfp == NULL || fcontext == NULL) return FALSE;
  awp = (Asn2gbWorkPtr) fcontext->userdata;
  if (awp == NULL) return FALSE;

  if (fcontext->seqfeattype == SEQFEAT_BIOSRC ||
      fcontext->seqfeattype == SEQFEAT_PUB) return TRUE;

  ivals = fcontext->ivals;
  numivals = fcontext->numivals;

  /* check to see if last interval is on this awp->from - awp->to range */

  if (ivals != NULL && numivals > 0) {
    idx = (numivals - 1) * 2;
    start = ivals [idx];
    stop = ivals [idx + 1];
    if (stop < awp->from || stop > awp->to) return TRUE;
  }

  /* !!! if RELEASE_MODE, verify that features have all mandatory qualifiers !!! */

  bbp = Asn2gbAddBlock (awp, FEATURE_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return TRUE;

  bbp->entityID = fcontext->entityID;
  bbp->itemID = fcontext->itemID;
  bbp->itemtype = OBJ_SEQFEAT;

  /* if CDS, explore mat_peptides, etc. */

  if (fcontext->seqfeattype != SEQFEAT_CDREGION) return TRUE;

  bsp = BioseqFindFromSeqLoc (sfp->product);
  if (bsp == NULL || (! ISA_aa (bsp->mol))) return TRUE;

  prt = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_PROT, 0, &pcontext);
  while (prt != NULL) {
    if (pcontext.left != 0 && pcontext.right != bsp->length - 1) {

      bbp = Asn2gbAddBlock (awp, FEATURE_BLOCK, sizeof (BaseBlock));
      if (bbp != NULL) {

        bbp->entityID = pcontext.entityID;
        bbp->itemID = pcontext.itemID;
        bbp->itemtype = OBJ_SEQFEAT;
      }
    }
    prt = SeqMgrGetNextFeature (bsp, prt, SEQFEAT_PROT, 0, &pcontext);
  }

  return TRUE;
}

static void AddFeatureBlock (Asn2gbWorkPtr awp)

{
  if (awp == NULL) return;

  SeqMgrExploreFeatures (awp->parent, (Pointer) awp, DoOneFeature, awp->slp, NULL, NULL);
}

static void AddBasecountBlock (Asn2gbWorkPtr awp)

{
  BaseBlockPtr  bbp;

  if (awp == NULL) return;

  bbp = Asn2gbAddBlock (awp, BASECOUNT_BLOCK, sizeof (BaseBlock));
}

static void AddOriginBlock (Asn2gbWorkPtr awp)

{
  BaseBlockPtr  bbp;

  if (awp == NULL) return;

  if (awp->format == EMBL_FMT || awp->format == EMBLPEPT_FMT) return;

  bbp = Asn2gbAddBlock (awp, ORIGIN_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  if (awp->format == GENBANK_FMT || awp->format == GENPEPT_FMT) {

    ff_StartPrint (0, 12, ASN2FF_GB_MAX, NULL);

    ff_AddString ("ORIGIN");

  }

  bbp->string = gb_MergeString (TRUE);
}

#define BASES_PER_BLOCK 1200

static void AddSequenceBlock (Asn2gbWorkPtr awp)

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

    stop = start + BASES_PER_BLOCK;
    if (stop >= len) {
      stop = len;
    }

    sbp->start = start;
    sbp->stop = stop;
  }
}

static void AddContigBlock (Asn2gbWorkPtr awp)

{
  BaseBlockPtr  bbp;

  if (awp == NULL) return;

  bbp = Asn2gbAddBlock (awp, CONTIG_BLOCK, sizeof (BaseBlock));
  if (bbp == NULL) return;

  bbp->string = StringSave ("CONTIG\n");
}

static void AddSlashBlock (Asn2gbWorkPtr awp)

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

static void DoOneSection (BioseqPtr parent, BioseqPtr bsp, SeqLocPtr slp,
                          Uint2 seg, Int4 from, Int4 to, Asn2gbWorkPtr awp)

{
  Asn2gbSectionPtr     asp;
  SeqMgrBioseqContext  bcontext;
  BaseBlockPtr         PNTR blockArray;
  Int2                 i;
  Int2                 numBlocks;
  Int4                 numsegs = 0;
  ValNodePtr           vnp;

  if (parent == NULL || bsp == NULL || awp == NULL) return;

  asp = Asn2gbAddSection (awp);
  if (asp == NULL) return;

  if (SeqMgrGetBioseqContext (parent, &bcontext)) {
    numsegs = bcontext.numsegs;
  }

  /* set working data fields */

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

  asp->parent = parent;
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

  AddLocusBlock (awp);

  if (awp->format == GENBANK_FMT || awp->format == GENPEPT_FMT) {

    AddDeflineBlock (awp);
    AddAccessionBlock (awp);

    if (ISA_na (bsp->mol)) {
      AddVersionBlock (awp);
    }

    if (ISA_aa (bsp->mol)) {
      AddPidBlock (awp);
      AddDbsourceBlock (awp);
    }

  } else if (awp->format == EMBL_FMT || awp->format == EMBLPEPT_FMT) {

    AddAccessionBlock (awp);

    if (ISA_na (bsp->mol)) {
      AddVersionBlock (awp);
    }

    if (ISA_aa (bsp->mol)) {
      AddPidBlock (awp);
      AddDbsourceBlock (awp);
    }

    AddDateBlock (awp);

    AddDeflineBlock (awp);
  }

  AddKeywordsBlock (awp);

  if (awp->format == GENBANK_FMT || awp->format == GENPEPT_FMT) {
    AddSegmentBlock (awp);
  }

  AddOrganismBlock (awp);

  AddReferenceBlock (awp, asp);
  AddCommentBlock (awp);

  AddFeatHeaderBlock (awp);
  AddSourceBlock (awp);
  AddFeatureBlock (awp);

  AddBasecountBlock (awp);
  AddOriginBlock (awp);

  if (awp->seqstyle == SEQUENCE_STYLE) {
    AddSequenceBlock (awp);
  } else if (awp->seqstyle == CONTIG_STYLE) {
    AddContigBlock (awp);
  }

  AddSlashBlock (awp);

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

static Boolean LIBCALLBACK Asn2Seg (SeqLocPtr slp, SeqMgrSegmentContextPtr context)

{
  Asn2gbWorkPtr  awp;
  BioseqPtr      bsp = NULL;
  Int4           from;
  SeqLocPtr      loc;
  SeqIdPtr       sip;
  Int4           to;

  if (slp == NULL || context == NULL) return FALSE;
  awp = (Asn2gbWorkPtr) context->userdata;

  sip = SeqLocId (slp);
  if (sip == NULL) {
    loc = SeqLocFindNext (slp, NULL);
    if (loc != NULL) {
      sip = SeqLocId (loc);
    }
  }
  if (sip != NULL) {
    bsp = BioseqFind (sip);
  }

  from = context->cumOffset;
  to = from + context->to - context->from;
  DoOneSection (context->parent, bsp, slp, context->index, from, to, awp);
  return TRUE;
}

static void DoOneBioseq (BioseqPtr bsp, Asn2gbWorkPtr awp)

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

    if (awp->segstyle == SEGMENTED_STYLE) {

      /* show all segments */

      SeqMgrExploreSegments (bsp, (Pointer) awp, Asn2Seg);

    } else if (awp->segstyle == MASTER_STYLE) {

      /* show as single bioseq */

      parent = bsp;
      from = 0;
      to = bsp->length - 1;

      DoOneSection (parent, bsp, ajp->slp, 0, from, to, awp);

    }

  } else if (bsp->repr == Seq_repr_raw || bsp->repr == Seq_repr_const || bsp->repr == Seq_repr_delta) {

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
    DoOneSection (parent, bsp, ajp->slp, 0, from, to, awp);
  }
}

static void DoPopPhyMutSet (SeqEntryPtr sep, Asn2gbWorkPtr awp)

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

NLM_EXTERN Asn2gbJobPtr asn2gnbk_setup (BioseqPtr bsp, BioseqSetPtr bssp,
                                        SeqLocPtr slp, Uint1 format, Uint1 mode,
                                        Uint1 segstyle, Uint1 seqstyle)

{
  Asn2gbJobPtr      ajp = NULL;
  Asn2gbSectionPtr  asp;
  Asn2gbWork        aw;
  BaseBlockPtr      bbp;
  BaseBlockPtr      PNTR blockArray;
  Uint2             entityID = 0;
  Int2              i;
  Int2              j;
  Int2              k;
  Int2              numBlocks;
  Int2              numSections;
  ObjMgrDataPtr     omdp;
  BaseBlockPtr      PNTR paragraphArray;
  Int2              numParagraphs;
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

  MemSet ((Pointer) (&aw), 0, sizeof (Asn2gbWork));
  aw.ajp = ajp;

  aw.sectionList = NULL;
  aw.lastsection = NULL;

  aw.currsection = 0;

  aw.format = format;
  aw.mode = mode;
  aw.segstyle = segstyle;
  aw.seqstyle = seqstyle;

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

  if (paragraphArray == NULL) return ajp;

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
          k++;
        }
      }
    }
  }

  /* free sectionList, but leave data, now pointed to by sectionArray elements */

  ValNodeFree (aw.sectionList);

  return ajp;
}

typedef CharPtr (*FormatProc) (Asn2gbFormatPtr afp, BaseBlockPtr bbp);

static FormatProc asn2gnbk_fmt_functions [21] = {
  NULL,
  DefaultFormatBlock, FormatDeflineBlock, DefaultFormatBlock,
  DefaultFormatBlock, DefaultFormatBlock, DefaultFormatBlock,
  DefaultFormatBlock, DefaultFormatBlock, DefaultFormatBlock,
  FormatOrganismBlock, FormatReferenceBlock, FormatCommentBlock,
  DefaultFormatBlock, FormatSourceBlock, FormatFeatureBlock,
  FormatBasecountBlock, DefaultFormatBlock, FormatSequenceBlock,
  FormatContigBlock, DefaultFormatBlock
};

NLM_EXTERN CharPtr asn2gnbk_format (Asn2gbJobPtr ajp, Int2 paragraph)

{
  Asn2gbFormat      af;
  Asn2gbSectionPtr  asp;
  BaseBlockPtr      bbp;
  Int2              blocktype;
  FormatProc        fmt;
  QualVal           qv [ASN2GNBK_TOTAL_SOURCE];
  Int2              section;
  CharPtr           str = NULL;

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

  fmt = asn2gnbk_fmt_functions [blocktype];
  if (fmt == NULL) return NULL;
  str = fmt (&af, bbp);

  if (str == NULL) {
    str = StringSave ("???\n");
  }

  return str;
}

NLM_EXTERN Asn2gbJobPtr asn2gnbk_cleanup (Asn2gbJobPtr ajp)

{
  Asn2gbSectionPtr   asp;
  BaseBlockPtr       bbp;
  BaseBlockPtr       PNTR blockArray;
  Int2               i;
  Int2               j;
  Int2               numBlocks;
  Int2               numSections;
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

  MemFree (ajp);
  return NULL;
}

NLM_EXTERN Boolean SeqEntryToGnbk (SeqEntryPtr sep, Uint1 format, Uint1 mode,
                                   Uint1 segstyle, Uint1 seqstyle, FILE *fp)

{
  Asn2gbJobPtr  ajp;
  BioseqPtr     bsp = NULL;
  BioseqSetPtr  bssp = NULL;
  Int2          i;
  Int2          numParagraphs;
  CharPtr       str;

  if (sep == NULL || fp == NULL) return FALSE;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
  }
  ajp = asn2gnbk_setup (bsp, bssp, NULL, format, mode, segstyle, seqstyle);
  if (ajp == NULL) return FALSE;

  if (fp != NULL) {
    numParagraphs = ajp->numParagraphs;
    for (i = 0; i < numParagraphs; i++) {
      str = asn2gnbk_format (ajp, i);
      fprintf (fp, "%s", str);
      MemFree (str);
    }
  }

  asn2gnbk_cleanup (ajp);
  return TRUE;
}


/* ********************************************************************** */

#define TEST_MODULE_ASN2GNBK

#ifdef TEST_MODULE_ASN2GNBK

/*
static void SummarizeQuals (void)

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

Args myargs [] = {
  {"Input File Name", "stdin", NULL, NULL,
    FALSE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
  {"EMBL format", "F", NULL, NULL,
    TRUE, 'e', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Master instead of segmented", "F", NULL, NULL,
    TRUE, 'm', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Contig instead of sequence", "F", NULL, NULL,
    TRUE, 'c', ARG_BOOLEAN, 0.0, 0, NULL},
};

Int2 Main (void)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Pointer       dataptr;
  Uint2         datatype;
  Uint2         entityID;
  Int2          format;
  FILE          *fp;
  Char          path [PATH_MAX];
  CharPtr       progname;
  Int2          segstyle;
  SeqEntryPtr   sep;
  Int2          seqstyle;

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
    Message (MSG_FATAL, "GetArgs failed");
    return 1;
  }

  fp = FileOpen (myargs [0].strvalue, "r");
  if (fp == NULL) {
    Message (MSG_FATAL, "FileOpen failed");
    return 1;
  }

  if (myargs [1].intvalue) {
    format = EMBL_FMT;
  } else {
    format = GENBANK_FMT;
  }

  if (myargs [2].intvalue) {
    segstyle = MASTER_STYLE;
  } else {
    segstyle = SEGMENTED_STYLE;
  }

  if (myargs [3].intvalue) {
    seqstyle = CONTIG_STYLE;
  } else {
    seqstyle = SEQUENCE_STYLE;
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
      fp = FileOpen ("newgenbankresults", "a");
      SeqEntryToGnbk (sep, format, RELEASE_MODE, segstyle, seqstyle, fp);
      FileClose (fp);
    }
  } else {
    Message (MSG_FATAL, "Datatype %d not recognized", (int) datatype);
  }
  ObjMgrFree (datatype, dataptr);

  return 0;
}

#endif /* TEST_MODULE_ASN2GNBK */

