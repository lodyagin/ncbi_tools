/*   valapi.c
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
* File Name:  valapi.c
*
* Author:  Colleen Bollin
*
* Version Creation Date:   4/7/2009
*
* $Revision: 1.33 $
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


/* TODO - create a superstructure for comment rules, so that they can be stored in a searchable list
 * and the match expressions can be precompiled.
 */

#include <asn.h>
#include <objfeat.h>
#include <subutil.h>
#include <objmgr.h>
#include <objfdef.h>
#include <gbftdef.h>
#include <sqnutils.h>
#include <edutil.h>
#include <gather.h>
#include <ffprint.h>
#include <asn2gnbi.h>
#include <findrepl.h>
#include <utilpub.h>
#include <regex.h>
#define NLM_GENERATED_CODE_PROTO
#include <objvalid.h>
#include <valapi.h>
#include "validrules.inc"

static CommentRulePtr CommentRules = NULL;


static int CompareFieldRules (FieldRulePtr r1, FieldRulePtr r2)
{
  int rval = 0;

  if (r1 == NULL && r2 == NULL) {
    rval = 0;
  } else if (r1 == NULL) {
    rval = -1;
  } else if (r2 == NULL) {
    rval = 1;
  } else {
    rval = StringCmp (r1->field_name, r2->field_name);
  }
  return rval;
}


NLM_EXTERN int LIBCALLBACK SortVnpByFieldRule (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      return CompareFieldRules(vnp1->data.ptrvalue, vnp2->data.ptrvalue);
    }
  }
  return 0;
}


static void SortFieldsInCommentRule (CommentRulePtr cr)
{
  FieldRulePtr r, prev = NULL;
  ValNodePtr list = NULL, vnp;

  if (cr == NULL || cr->fields == NULL) {
      //LCOV_EXCL_START
      //not valid, not testable
    return;
    //LCOV_EXCL_STOP
  }

  for (r = cr->fields; r != NULL; r = r->next) {
    ValNodeAddPointer (&list, 0, r);
  }

  list = ValNodeSort (list, SortVnpByFieldRule);

  prev = list->data.ptrvalue;
  cr->fields = prev;
  for (vnp = list->next; vnp != NULL; vnp = vnp->next) {
    prev->next = vnp->data.ptrvalue;
    prev = prev->next;
  }
  prev->next = NULL;
  list = ValNodeFree (list);
}


static void SortCommentRuleFields (CommentSetPtr cr_set)
{
  while (cr_set != NULL) {
    if (!cr_set->require_order) {
      SortFieldsInCommentRule(cr_set);
    }
    cr_set = cr_set->next;
  }
}


//LCOV_EXCL_START
//not testable in regression
static Boolean LoadCommentRulesFromLocalString (void)

{
#ifndef WIN16
  AsnIoMemPtr aimp;
  CharPtr     ptr;

  ptr = MergeStringArray ((CharPtr PNTR) s_Defaultvalidrules, sizeof (s_Defaultvalidrules) / sizeof (char*));
  if (ptr == NULL) return FALSE;

  aimp = AsnIoMemOpen ("r", (BytePtr) ptr, (Int4) StringLen (ptr));
  if (aimp == NULL || aimp->aip == NULL) return FALSE;

  CommentRules = CommentSetAsnRead (aimp->aip, NULL);
  SortCommentRuleFields(CommentRules);
  AsnIoMemClose (aimp);
  MemFree (ptr);
#endif
  return (Boolean) (CommentRules != NULL);
}
//LCOV_EXCL_STOP


NLM_EXTERN CommentRulePtr LoadCommentRuleSet (void)
{
    Char buf[PATH_MAX];
    AsnIoPtr aip;

    if (CommentRules != NULL)
        return CommentRules;

    if (! FindPath("ncbi", "ncbi", "data", buf, sizeof (buf)))
    {
        //LCOV_EXCL_START
        //not testable in regression
        if (LoadCommentRulesFromLocalString ()) {
            return CommentRules;
        }

        ErrPostEx(SEV_WARNING, 0, 0, "FindPath failed in LoadCommentRuleSet - ncbi configuration file missing or incorrect");
        return CommentRules;
        //LCOV_EXCL_STOP
    }

    StringCat(buf, "validrules.prt");
    if ((aip = AsnIoOpen(buf, "r")) == NULL)
    {
        //LCOV_EXCL_START
        //not testable in regression
        if (LoadCommentRulesFromLocalString ()) {
            return CommentRules;
        }

        ErrPostEx(SEV_WARNING, 0, 0, "Couldn't open [%s]", buf);
        return CommentRules;
        //LCOV_EXCL_STOP
    }

    CommentRules = CommentSetAsnRead(aip, NULL);
    AsnIoClose(aip);
    SortCommentRuleFields (CommentRules);
    return CommentRules;
}


NLM_EXTERN CommentRulePtr GetCommentRuleFromRuleSet (CharPtr prefix)
{
  CommentRulePtr cr;
  Boolean        found = FALSE;

  if (CommentRules == NULL) {
    cr = LoadCommentRuleSet ();
  } else {
    cr = CommentRules;
  }

  while (cr != NULL && !found) {
    if (StringICmp (cr->prefix, prefix) == 0) {
      found = TRUE;
    } else {
      cr = cr->next;
    }
  }
  return cr;
}


static Boolean DoesFieldNameMatchRuleName (UserFieldPtr ufp, CharPtr rule_field_name)
{
  Char buf[15];
  Boolean rval = FALSE;

  if (ufp == NULL) {
    return FALSE;
  }

  /* Does this match the rule field name? */
  if (ufp->label == NULL) {
    rval = FALSE;
  } else if (ufp->label->id > 0) {
    sprintf (buf, "%d", ufp->label->id);
    if (StringCmp (rule_field_name, buf) == 0) {
      rval = TRUE;
    }
  } else if (StringCmp (rule_field_name, ufp->label->str) == 0) {
    rval = TRUE;
  }
  return rval;
}


static Boolean DoesStringMatchExpression (CharPtr str, CharPtr expression)
{
  Boolean rval;
  regex_t buffer;
  const char *errstr;
  Int4  len;

  MemSet (&buffer, 0, sizeof (regex_t));
  errstr = re_compile_pattern (expression, StringLen (expression), &buffer);

  len = StringLen (str);
  if (re_match (&buffer, str, len, 0, NULL) > -1) {
    rval = TRUE;
  } else {
    rval = FALSE;
  }

  return rval;
}


static Boolean DoesFieldValueMatchExpression (UserFieldPtr ufp, CharPtr match_expression)
{
  Char buf[15];
  Boolean rval = FALSE;

  if (ufp->choice == 1) {
    /* compare string value */
    if (DoesStringMatchExpression(ufp->data.ptrvalue, match_expression)) {
      rval = TRUE;
    }
  } else if (ufp->choice == 2) {
    /* compare int value */
    sprintf (buf, "%d", ufp->data.intvalue);
    if (DoesStringMatchExpression(buf, match_expression)) {
      rval = TRUE;
    }
  } else {
    /* for now, only allowing matches for int or string fields */
    rval = FALSE;
  }
  return rval;
}


static Boolean DoesFieldValueMatchRule (UserFieldPtr ufp, FieldRulePtr rule)
{
  Boolean rval = FALSE;

  if (ufp == NULL || rule == NULL) {
    return FALSE;
  }
  /* now validate field value */
  if (rule->match_expression == NULL) {
    rval = TRUE;
  } else {
    rval = DoesFieldValueMatchExpression (ufp, rule->match_expression);
  }
  return rval;
}


static Boolean DoesStructuredCommentHavePrefix (UserObjectPtr uop, CharPtr prefix)
{
  UserFieldPtr ufp;
  Int4         prefix_len;
  Boolean      rval = FALSE;

  if (uop == NULL || uop->type == NULL || StringICmp (uop->type->str, "StructuredComment") != 0) {
    return FALSE;
  } else if (StringHasNoText (prefix)) {
    return TRUE;
  }

  prefix_len = StringLen (prefix);

  for (ufp = uop->data; ufp != NULL && !rval; ufp = ufp->next) {
    if (ufp->label != NULL 
        && (StringICmp (ufp->label->str, "StructuredCommentPrefix") == 0
            || StringICmp (ufp->label->str, "StructuredCommentSuffix") == 0)
        && ufp->choice == 1
        && StringNCmp (ufp->data.ptrvalue, prefix, prefix_len) == 0) {
      rval = TRUE;
    }
  }
  return rval;
}


static FieldRulePtr FindRuleForFieldName (UserFieldPtr ufp, FieldRulePtr field_rule)
{
  if (ufp == NULL) {
    return NULL;
  }

  while (field_rule != NULL && ! DoesFieldNameMatchRuleName(ufp, field_rule->field_name)) {
    field_rule = field_rule->next;
  }
  return field_rule;
}


static UserFieldPtr FindFieldForRuleName (UserFieldPtr ufp, CharPtr field_rule)
{
  if (StringHasNoText (field_rule)) {
    return NULL;
  }

  while (ufp != NULL && !DoesFieldNameMatchRuleName(ufp, field_rule)) {
    ufp = ufp->next;
  }
  return ufp;
}


static int VACompareUserFields (UserFieldPtr ufp1, UserFieldPtr ufp2)
{
  int rval = 0;
  CharPtr cp1, cp2;
  Char buf1[15];
  Char buf2[15];

  if (ufp1 == NULL && ufp2 == NULL) {
    rval = 0;
  } else if (ufp1 == NULL) {
    rval = -1;
  } else if (ufp2 == NULL) {
    rval = 1;
  } else if (ufp1->label == NULL) {
    rval = -1;
  } else if (ufp2->label == NULL) {
    rval = 1;
  } else if (ufp1->label->id > 0 && ufp2->label->id > 0) {
    if (ufp1->label->id > ufp2->label->id) {
      rval = 1;
    } else if (ufp1->label->id < ufp2->label->id) {
      rval = -1;
    } else {
      rval = 0;
    }
  } else {
    if (ufp1->label->id > 0) {
      sprintf (buf1, "%d", ufp1->label->id);
      cp1 = buf1;
    } else {
      cp1 = ufp1->label->str;
    }
    if (ufp2->label->id > 0) {
      sprintf (buf2, "%d", ufp2->label->id);
      cp2 = buf2;
    } else {
      cp2 = ufp2->label->str;
    }
    rval = StringCmp (cp1, cp2);
  }
  return rval;
}


NLM_EXTERN int LIBCALLBACK SortVnpByUserField (VoidPtr ptr1, VoidPtr ptr2)

{
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      return VACompareUserFields(vnp1->data.ptrvalue, vnp2->data.ptrvalue);
    }
  }
  return 0;
}


static void SortFieldsInUserObject (UserObjectPtr uop)
{
  UserFieldPtr ufp, prev = NULL;
  ValNodePtr list = NULL, vnp;

  if (uop == NULL || uop->data == NULL) {
    //LCOV_EXCL_START
    //invalid
    return;
    //LCOV_EXCL_STOP
  }

  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    ValNodeAddPointer (&list, 0, ufp);
  }

  list = ValNodeSort (list, SortVnpByUserField);

  prev = list->data.ptrvalue;
  uop->data = prev;
  for (vnp = list->next; vnp != NULL; vnp = vnp->next) {
    prev->next = vnp->data.ptrvalue;
    prev = prev->next;
  }
  prev->next = NULL;
  list = ValNodeFree (list);
}

static Boolean PrefixOrSuffixInList (CharPtr val, CharPtr before, CharPtr after)

{
  Char        buf [1024];
  size_t      len, l_before, l_after;
  ValNodePtr  list, vnp;
  Boolean     rsult = FALSE;
  CharPtr     str;

  if (val == NULL) return FALSE;
  len = StringLen (val);
  if (len < 10) return FALSE;
  l_before = StringLen (before);
  l_after = StringLen (after);
  if (StringNCmp (val, before, l_before) != 0) return FALSE;
  if (StringNCmp (val + len - l_after, after, l_after) != 0) return FALSE;

  if (len > sizeof (buf)) return FALSE;
  StringNCpy_0 (buf, val + l_before, sizeof (buf));
  buf [len - l_before - l_after] = '\0';

  list = GetStructuredCommentPrefixList ();
  if (list == NULL) return FALSE;

  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    if (StringCmp (buf, str) == 0) rsult = TRUE;
  }

  ValNodeFreeData (list);

  return rsult;
}


static EFieldValid FindForbiddenPhrases
(UserObjectPtr uop,
 CommentRulePtr comment_rule,
 StructuredCommentCallback s_callback,
 Pointer s_callback_data)
{
  UserFieldPtr ufp;
  ValNodePtr   vnp;
  EFieldValid  rval = eFieldValid_Valid;

  if (uop == NULL || comment_rule == NULL || comment_rule->forbidden_phrases == NULL) {
    return eFieldValid_Valid;
  }

  /* examine fields for forbidden phrases */
  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    if (ufp->label != NULL 
        && (StringICmp (ufp->label->str, "StructuredCommentPrefix") == 0
            || StringICmp (ufp->label->str, "StructuredCommentSuffix") == 0)) {
        /* skip suffix and prefix */
        continue;
    } else {
      if (ufp->choice == 1) {
        /* compare string value */
        for (vnp = comment_rule->forbidden_phrases; vnp != NULL; vnp = vnp->next) {
          if (StringISearch(ufp->data.ptrvalue, vnp->data.ptrvalue) != NULL) {
            rval = eFieldValid_Inappropriate;
            if (s_callback == NULL) {
              break;
            } else {
              s_callback (eFieldValid_Inappropriate, NULL, ufp, NULL, s_callback_data, uop);
            }
          }
        }
      }
    }
  }
  return rval;
}


static EFieldValid AreStructuredCommentContentsValidForRule
(UserObjectPtr uop,
 CommentRulePtr comment_rule,
 StructuredCommentCallback s_callback,
 Pointer s_callback_data)
{
  UserFieldPtr ufp, ufp_tmp, depend_ufp;
  FieldRulePtr field_rule, rule_tmp;
  DependentFieldRulePtr depend_rule;
  EFieldValid  rval = eFieldValid_Valid, tmp_val;
  Boolean      free_uop = FALSE;

  if (uop == NULL || uop->type == NULL || StringICmp (uop->type->str, "StructuredComment") != 0) {
    return eFieldValid_Invalid;
  }

  if (!comment_rule->require_order) {
    free_uop = TRUE;
    uop = (UserObjectPtr) AsnIoMemCopy (uop, (AsnReadFunc) UserObjectAsnRead, (AsnWriteFunc) UserObjectAsnWrite);
    SortFieldsInUserObject(uop);
  }

  /* now check individual fields */
  ufp = uop->data;
  field_rule = comment_rule->fields;

  while (field_rule != NULL && ufp != NULL) {
    if (ufp->label != NULL 
        && (StringICmp (ufp->label->str, "StructuredCommentPrefix") == 0
            || StringICmp (ufp->label->str, "StructuredCommentSuffix") == 0)) {
      /* skip suffix and prefix */
      ufp = ufp->next;
    } else if (DoesFieldNameMatchRuleName (ufp, field_rule->field_name)) {
      if (DoesFieldValueMatchRule (ufp, field_rule)) {
        /* all good */
      } else {
        if (s_callback == NULL) {
          rval = eFieldValid_Invalid;
          goto IsStructuredCommentValidForRule_exit;
        } else {
          if (rval == eFieldValid_Valid) {
            rval = eFieldValid_Invalid;
          }
          s_callback (eFieldValid_Invalid, field_rule, ufp, NULL, s_callback_data, uop);
        }
      }
      ufp = ufp->next;
      field_rule = field_rule->next;
    } else if ((rule_tmp = FindRuleForFieldName(ufp, comment_rule->fields)) == NULL) {
      if (!comment_rule->allow_unlisted) {
        if (s_callback == NULL) {
          rval = eFieldValid_Invalid;
          goto IsStructuredCommentValidForRule_exit;
        } else {
          if (rval == eFieldValid_Valid) {
            rval = eFieldValid_Invalid;
          }
          s_callback (eFieldValid_Invalid, NULL, ufp, NULL, s_callback_data, uop);
        }
      }
      ufp = ufp->next;
    } else {
      ufp_tmp = FindFieldForRuleName (ufp, field_rule->field_name);
      if (ufp_tmp == NULL) {
        if (field_rule->required) {
          if (s_callback == NULL) {
            rval = eFieldValid_MissingRequiredField;
            goto IsStructuredCommentValidForRule_exit;
          } else {
            if (rval == eFieldValid_Valid) {
              rval = eFieldValid_MissingRequiredField;
            }
            s_callback (eFieldValid_MissingRequiredField, field_rule, NULL, NULL, s_callback_data, uop);
          }
        } else {
          /* field wasn't required, it's ok */
        }
      } else if (comment_rule->require_order) {
        /* field is out of order */
        if (s_callback == NULL) {
          rval = eFieldValid_FieldOutOfOrder;
          goto IsStructuredCommentValidForRule_exit;
        } else {
          if (rval == eFieldValid_Valid) {
            rval = eFieldValid_FieldOutOfOrder;
          }
          s_callback (eFieldValid_FieldOutOfOrder, field_rule, ufp_tmp, NULL, s_callback_data, uop);
          if (!DoesFieldValueMatchRule (ufp_tmp, field_rule)) {
            s_callback (eFieldValid_Invalid, field_rule, ufp_tmp, NULL, s_callback_data, uop);
          }
        }
      } else {
        if (!DoesFieldValueMatchRule(ufp_tmp, field_rule)) {
            if (s_callback == NULL) {
                rval = eFieldValid_Invalid;
                goto IsStructuredCommentValidForRule_exit;
            } else {
                s_callback(eFieldValid_Invalid, field_rule, ufp_tmp, NULL, s_callback_data, uop);
            }
        }
      }
      field_rule = field_rule->next;
    }
  }

  /* validate remaining rules */
  while (field_rule != NULL) {
    if (field_rule->required) {
      if (s_callback == NULL) {
        rval = eFieldValid_MissingRequiredField;
        goto IsStructuredCommentValidForRule_exit;
      } else {
        if (rval == eFieldValid_Valid) {
          rval = eFieldValid_MissingRequiredField;
        }
        s_callback (eFieldValid_MissingRequiredField, field_rule, NULL, NULL, s_callback_data, uop);
      }
    }
    field_rule = field_rule->next;
  }

  while (ufp != NULL) {
    if (ufp->label != NULL 
      && (StringICmp (ufp->label->str, "StructuredCommentPrefix") == 0
          || StringICmp (ufp->label->str, "StructuredCommentSuffix") == 0)) {
      /* skip suffix and prefix */
      ufp = ufp->next;
      continue;
    }
    field_rule = FindRuleForFieldName (ufp, comment_rule->fields);
    if (field_rule == NULL && !comment_rule->allow_unlisted) {
      if (s_callback == NULL) {
        rval = eFieldValid_Invalid;
        goto IsStructuredCommentValidForRule_exit;
      } else {
        if (rval == eFieldValid_Valid) {
          rval = eFieldValid_Invalid;
        }
        s_callback (eFieldValid_Invalid, NULL, ufp, NULL, s_callback_data, uop);
      }
    } else if (field_rule != NULL && FindFieldForRuleName(uop->data, field_rule->field_name) != ufp) {
      if (s_callback == NULL) {
        rval = eFieldValid_DuplicateField;
        goto IsStructuredCommentValidForRule_exit;
      } else {
        if (rval == eFieldValid_Valid) {
          rval = eFieldValid_DuplicateField;
        }
        s_callback (eFieldValid_DuplicateField, field_rule, ufp, NULL, s_callback_data, uop);
      }
    }
    ufp = ufp->next;
  }

  /* now look at fields in context of other fields */
  for (depend_rule = comment_rule->dependent_rules; depend_rule != NULL; depend_rule = depend_rule->next) {
    depend_ufp = FindFieldForRuleName(uop->data, depend_rule->match_name);
    if (depend_ufp != NULL && 
        ((!depend_rule->invert_match && DoesFieldValueMatchExpression(depend_ufp, depend_rule->value_constraint)) 
        || (depend_rule->invert_match && !DoesFieldValueMatchExpression(depend_ufp, depend_rule->value_constraint)))) {
      for (field_rule = depend_rule->other_fields; field_rule != NULL; field_rule = field_rule->next) {
        ufp = FindFieldForRuleName (uop->data, field_rule->field_name);
        if (ufp == NULL) {
          if (s_callback == NULL) {
            rval = eFieldValid_MissingRequiredField;
            goto IsStructuredCommentValidForRule_exit;
          } else {
            if (rval == eFieldValid_Valid) {
              rval = eFieldValid_MissingRequiredField;
            }
            s_callback (eFieldValid_MissingRequiredField, field_rule, ufp, depend_ufp, s_callback_data, uop);
          }
        } else if (!DoesFieldValueMatchRule (ufp, field_rule)) {
          if (s_callback == NULL) {
            rval = eFieldValid_Invalid;
            goto IsStructuredCommentValidForRule_exit;
          } else {
            if (rval == eFieldValid_Valid) {
              rval = eFieldValid_MissingRequiredField;
            }
            s_callback (eFieldValid_Invalid, field_rule, ufp, depend_ufp, s_callback_data, uop);
          }
        }
      }
      for (field_rule = depend_rule->disallowed_fields; field_rule != NULL; field_rule = field_rule->next) {
          //LCOV_EXCL_START
          //no rules currently have disallowed fields
        ufp = FindFieldForRuleName (uop->data, field_rule->field_name);
        if (ufp != NULL && DoesFieldValueMatchRule (ufp, field_rule)) {
          if (s_callback == NULL) {
            rval = eFieldValid_Disallowed;
            goto IsStructuredCommentValidForRule_exit;
          } else {
            if (rval == eFieldValid_Valid) {
              rval = eFieldValid_Disallowed;
            }
            s_callback (eFieldValid_Disallowed, field_rule, ufp, depend_ufp, s_callback_data, uop);
          }
        }
        //LCOV_EXCL_STOP
      }
    }
  }

  tmp_val = FindForbiddenPhrases(uop, comment_rule, s_callback, s_callback_data);
  if (rval == eFieldValid_Valid) {
    rval = tmp_val;
  }

IsStructuredCommentValidForRule_exit:
  if (free_uop) {
    uop = UserObjectFree (uop);
  }
  return rval;
}


NLM_EXTERN EFieldValid 
IsStructuredCommentValidForRule 
(UserObjectPtr uop,
 CommentRulePtr comment_rule,
 StructuredCommentCallback s_callback,
 Pointer s_callback_data)
{
  UserFieldPtr ufp;
  EFieldValid  rval = eFieldValid_Valid;

  if (uop == NULL || uop->type == NULL || StringICmp (uop->type->str, "StructuredComment") != 0) {
    return eFieldValid_Invalid;
  }

  for (ufp = uop->data; ufp != NULL; ufp = ufp->next) {
    if (ufp->label != NULL && StringICmp (ufp->label->str, "StructuredCommentPrefix") == 0) {
      /* check prefix */
      if (! PrefixOrSuffixInList ((CharPtr) ufp->data.ptrvalue, "##", "-START##")) {
        if (s_callback == NULL) {
          return eFieldValid_Invalid;
        } else {
          if (rval == eFieldValid_Valid) {
            rval = eFieldValid_Invalid;
          }
          s_callback (eFieldValid_Invalid, NULL, ufp, NULL, s_callback_data, uop);
        }
      }
    } else if (ufp->label != NULL && StringICmp (ufp->label->str, "StructuredCommentSuffix") == 0) {
      /* check suffix */
      if (! PrefixOrSuffixInList ((CharPtr) ufp->data.ptrvalue, "##", "-END##")) {
        if (s_callback == NULL) {
          return eFieldValid_Invalid;
        } else {
          if (rval == eFieldValid_Valid) {
            rval = eFieldValid_Invalid;
          }
          s_callback (eFieldValid_Invalid, NULL, ufp, NULL, s_callback_data, uop);
        }
      }
    }
  }

  if (comment_rule == NULL) {
    return rval;
  }

  /* first, make sure comment rule prefix matches comment */
  if (!DoesStructuredCommentHavePrefix (uop, comment_rule->prefix)) {
    return eFieldValid_Invalid;
  }

  return AreStructuredCommentContentsValidForRule (uop, comment_rule, s_callback, s_callback_data);
}


NLM_EXTERN EFieldValid IsStructuredCommentValid (UserObjectPtr uop, StructuredCommentCallback s_callback, Pointer s_callback_data)
{
  CommentRulePtr cr;
  UserFieldPtr ufp;
  CharPtr prefix = NULL;

  for (ufp = uop->data; ufp != NULL && prefix == NULL; ufp = ufp->next) {
    if (ufp->label != NULL 
        && StringICmp (ufp->label->str, "StructuredCommentPrefix") == 0
        && ufp->choice == 1) {
      prefix = ufp->data.ptrvalue;
    }
  }

  if (prefix == NULL) {
    return TRUE;
  }
  cr = GetCommentRuleFromRuleSet (prefix);
  return IsStructuredCommentValidForRule (uop, cr, s_callback, s_callback_data);
}


static Boolean MovePrefixAndSuffixFieldsToFlank (UserObjectPtr uop)
{
  UserFieldPtr ufp, ufp_prev = NULL, ufp_next, ufp_last = NULL;
  Boolean changed = FALSE;

  if (uop == NULL) {
    return FALSE;
  }

  for (ufp = uop->data; ufp != NULL; ufp = ufp_next) {
    ufp_next = ufp->next;
    if (IsStructuredCommentPrefix(ufp)) {
      /* move prefix to first position */
      if (ufp_prev != NULL) {
        ufp_prev->next = ufp_next;
        ufp->next = uop->data;
        uop->data = ufp;
        changed = TRUE;
      } else  {
        ufp_prev = ufp;
      }
    } else if (IsStructuredCommentSuffix(ufp) && ufp->next != NULL) {
      changed = TRUE;
      ufp_last = ufp;
      if (ufp_prev == NULL) {
        uop->data = ufp_next;
      } else {
        ufp_prev->next = ufp_next;
      }
      ufp->next = NULL;
    } else {
      ufp_prev = ufp;
    }
  }
  if (ufp_last != NULL) {
    if (ufp_prev == NULL) {
      uop->data = ufp_last;
    } else {
      ufp_prev->next = ufp_last;
    }
  }
  return changed;
}


//LCOV_EXCL_START
//not used for validation
NLM_EXTERN Boolean ReorderStructuredCommentFields (UserObjectPtr uop)
{
  CommentRulePtr cr;
  FieldRulePtr rule;
  UserFieldPtr ufp, ufp_prev = NULL, new_list = NULL, new_prev = NULL;
  CharPtr prefix = NULL;
  UserObjectPtr uop_orig;
  Boolean changed = FALSE;

  if (uop == NULL || uop->type == NULL || StringICmp (uop->type->str, "StructuredComment") != 0) return FALSE;

  uop_orig = AsnIoMemCopy (uop, (AsnReadFunc) UserObjectAsnRead, (AsnWriteFunc) UserObjectAsnWrite);
  prefix = GetStructuredCommentPrefix (uop);
  if (prefix != NULL
      && (cr = GetCommentRuleFromRuleSet (prefix)) != NULL) {
    ufp = uop->data;
    for (rule = cr->fields; rule != NULL; rule = rule->next) {
      ufp_prev = NULL;
      for (ufp = uop->data; ufp != NULL && StringCmp (ufp->label->str, rule->field_name) != 0; ufp = ufp->next) {
        ufp_prev = ufp;
      }
      if (ufp != NULL) {
        if (ufp_prev == NULL) {
          uop->data = ufp->next;
        } else {
          ufp_prev->next = ufp->next;
        }
        ufp->next = NULL;
        if (new_prev == NULL) {
          new_list = ufp;
        } else {
          new_prev->next = ufp;
        }
        new_prev = ufp;
      }
    }
    if (new_prev != NULL) {
      new_prev->next = uop->data;
      uop->data = new_list;
    }
    changed = TRUE;
  } 
  changed |= MovePrefixAndSuffixFieldsToFlank (uop);

  changed = !AsnIoMemComp (uop, uop_orig, (AsnWriteFunc) UserObjectAsnWrite);
  uop_orig = UserObjectFree (uop_orig);

  return changed;
}

//not used for validation
static void ReorderStructuredCommentFieldsCallback (SeqDescrPtr sdp, Pointer userdata)
{
  BoolPtr r = (BoolPtr) userdata;
  if (sdp->choice == Seq_descr_user) {
    if (ReorderStructuredCommentFields ((UserObjectPtr) sdp->data.ptrvalue)) {
      if (r != NULL) {
        *r = TRUE;
      }
    }
  }
}


//not used for validation
NLM_EXTERN Boolean ReorderStructuredCommentsInSeqEntry (SeqEntryPtr sep)
{
  Boolean rval = FALSE;

  VisitDescriptorsInSep (sep, &rval, ReorderStructuredCommentFieldsCallback);
  return rval;
}


//not used for validation
static Boolean DoesStructuredCommentHaveAnyPrefixOrSuffix (UserObjectPtr uop)
{
  UserFieldPtr ufp;
  Boolean      rval = FALSE;

  if (uop == NULL || uop->type == NULL || StringICmp (uop->type->str, "StructuredComment") != 0) {
    return FALSE;
  }

  for (ufp = uop->data; ufp != NULL && !rval; ufp = ufp->next) {
    if (ufp->label != NULL 
        && (StringICmp (ufp->label->str, "StructuredCommentPrefix") == 0
          || StringICmp (ufp->label->str, "StructuredCommentSuffix") == 0)) {
      rval = TRUE;
    }
  }
  return rval;
}


//LCOV_EXCL_START
//not used in validation
static Boolean IsRuleOkForStructuredComment (UserObjectPtr uop, CommentRulePtr cr)
{
  UserFieldPtr ufp;
  FieldRulePtr field_rule;
  Boolean  rval = FALSE;
  CommentRulePtr cr_tmp = NULL;

  if (uop == NULL || uop->data == NULL || cr == NULL) {
    return FALSE;
  }
  ufp = uop->data;

  /* all field names must be recognized */
  while (ufp != NULL) {
    if (!IsStructuredCommentPrefix(ufp) && !IsStructuredCommentSuffix(ufp)) {
      /* ignore prefix and suffix if present */
      field_rule = FindRuleForFieldName(ufp, cr->fields);
      if (field_rule == NULL) {
        return FALSE;
      }
    }
    ufp = ufp->next;
  }

  /* don't be picky about order when assigning prefix */
  if (cr->require_order) {
    cr_tmp = (CommentRulePtr) AsnIoMemCopy (cr, (AsnReadFunc) CommentRuleAsnRead, (AsnWriteFunc) CommentRuleAsnWrite);
    cr_tmp->require_order = FALSE;
    cr = cr_tmp;
  }
  if (AreStructuredCommentContentsValidForRule (uop, cr, NULL, NULL) == eFieldValid_Valid) {
    rval = TRUE;
  }
  cr_tmp = CommentRuleFree (cr_tmp);
  return rval;
}


//not used for validation
NLM_EXTERN CharPtr AutoapplyStructuredCommentPrefix (UserObjectPtr uop)
{
  CommentRulePtr cr;
  CharPtr        prefix = NULL;

  if (uop == NULL || DoesStructuredCommentHaveAnyPrefixOrSuffix(uop)) {
    return NULL;
  }

  if (CommentRules == NULL) {
    cr = LoadCommentRuleSet ();
  } else {
    cr = CommentRules;
  }

  while (cr != NULL) {
    if (IsRuleOkForStructuredComment(uop, cr)) {
      if (prefix == NULL) {
        prefix = cr->prefix;
      } else {
        return NULL;
      }
    }
    cr = cr->next;
  }

  if (prefix != NULL) {
    AddItemStructuredCommentUserObject (uop, "StructuredCommentPrefix", prefix);
  }
  return prefix;
}


//not used for validation
NLM_EXTERN CommentRulePtr NewRuleForStructuredComment (UserObjectPtr uop)
{
  CommentRulePtr cr;
  CommentRulePtr new_cr = NULL;

  if (uop == NULL || uop->type == NULL 
      || StringICmp (uop->type->str, "StructuredComment") != 0
      || !DoesStructuredCommentHaveAnyPrefixOrSuffix (uop)
      || IsStructuredCommentValid (uop, NULL, NULL) == eFieldValid_Valid) {
    return NULL;
  }

  if (CommentRules == NULL) {
    cr = LoadCommentRuleSet ();
  } else {
    cr = CommentRules;
  }

  while (cr != NULL) {
    if (IsRuleOkForStructuredComment(uop, cr)) {
      if (new_cr == NULL) {
        new_cr = cr;
      } else {
        return NULL;
      }
    }
    cr = cr->next;
  }

  return new_cr;
}
//LCOV_EXCL_STOP

