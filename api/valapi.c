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

static CommentRulePtr CommentRules = NULL;

#ifndef WIN16
static CharPtr commentRulesStr = "Comment-set ::= {\n" \
"{ prefix \"foo\" , fields { { field-name \"bar\" , match-expression \"^s\" } , { field-name \"baz\" , match-expression \"\\d+\\s+z\" } } } , \n " \
"{ prefix \"foo2\" , fields { { field-name \"bar\" , match-expression \"b\\W+a\" } , { field-name \"baz\" , match-expression \"z\\s+\\d+\" , required TRUE } } } , \n " \
"}\n";
#endif


static Boolean LoadCommentRulesFromLocalString (void)

{
#ifndef WIN16
  AsnIoMemPtr aimp;

  aimp = AsnIoMemOpen ("r", (BytePtr) commentRulesStr, (Int4) StringLen (commentRulesStr));
  if (aimp == NULL || aimp->aip == NULL) return FALSE;

  CommentRules = CommentSetAsnRead (aimp->aip, NULL);
  AsnIoMemClose (aimp);
#endif
  return (Boolean) (CommentRules != NULL);
}


NLM_EXTERN CommentRulePtr LoadCommentRuleSet (void)
{
    Char buf[PATH_MAX];
    AsnIoPtr aip;

    if (CommentRules != NULL)
        return CommentRules;

#ifdef OS_UNIX
    if (getenv ("USE_VALIDRULES_FILE") == NULL) {
          if (LoadCommentRulesFromLocalString ()) {
              return CommentRules;
          }
    }
#endif

    if (! FindPath("ncbi", "ncbi", "data", buf, sizeof (buf)))
    {

        if (LoadCommentRulesFromLocalString ()) {
            return CommentRules;
        }

        ErrPostEx(SEV_WARNING, 0, 0, "FindPath failed in LoadCommentRuleSet - ncbi configuration file missing or incorrect");
        return CommentRules;
    }

    StringCat(buf, "validrules.val");
    if ((aip = AsnIoOpen(buf, "rb")) == NULL)
    {

        if (LoadCommentRulesFromLocalString ()) {
            return CommentRules;
        }

        ErrPostEx(SEV_WARNING, 0, 0, "Couldn't open [%s]", buf);
        return CommentRules;
    }

    CommentSetAsnRead(aip, NULL);

    AsnIoClose(aip);
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


static EFieldValid IsStructuredCommentFieldValid (UserFieldPtr ufp, FieldRulePtr rule)
{
  Char buf[15];
  EFieldValid rval = eFieldValid_Invalid;

  if (rule == NULL) {
    return eFieldValid_Valid;
  }
  if (ufp == NULL) {
    if (rule->required) {
      return eFieldValid_MissingRequiredField;
    } else {
      return eFieldValid_OptionalFieldSkipped;
    }
  }

  /* Does this match the rule field name? */
  if (ufp->label == NULL) {
    return eFieldValid_Invalid;
  } else if (ufp->label->id > 0) {
    sprintf (buf, "%d", ufp->label->id);
    if (StringCmp (rule->field_name, buf) != 0) {
      if (rule->required) {
        return eFieldValid_MissingRequiredField;
      } else {
        return eFieldValid_OptionalFieldSkipped;
      }
    }
  } else if (StringCmp (rule->field_name, ufp->label->str) != 0) {
    if (rule->required) {
      return eFieldValid_MissingRequiredField;
    } else {
      return eFieldValid_OptionalFieldSkipped;
    }
  }

  /* now validate field value */
  if (rule->match_expression == NULL) {
    return eFieldValid_Valid;
  }

  if (ufp->choice == 1) {
    /* compare string value */
    if (DoesStringMatchExpression(ufp->data.ptrvalue, rule->match_expression)) {
      rval = eFieldValid_Valid;
    }
  } else if (ufp->choice == 2) {
    /* compare int value */
    sprintf (buf, "%d", ufp->data.intvalue);
    if (DoesStringMatchExpression(buf, rule->match_expression)) {
      rval = eFieldValid_Valid;
    }
  } else {
    /* for now, only allowing matches for int or string fields */
    rval = eFieldValid_Invalid;
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


NLM_EXTERN EFieldValid IsStructuredCommentValidForRule (UserObjectPtr uop, CommentRulePtr comment_rule)
{
  UserFieldPtr ufp;
  FieldRulePtr field_rule;
  EFieldValid  valid;
  EFieldValid  rval = eFieldValid_Valid;

  if (uop == NULL || uop->type == NULL || StringICmp (uop->type->str, "StructuredComment") != 0) {
    return eFieldValid_Invalid;
  } else if (comment_rule == NULL) {
    return eFieldValid_Valid;
  }

  /* first, make sure comment rule prefix matches comment */
  if (!DoesStructuredCommentHavePrefix (uop, comment_rule->prefix)) {
    return eFieldValid_Invalid;
  }
  /* now check fields in order */
  ufp = uop->data;
  field_rule = comment_rule->fields;

  while (field_rule != NULL && ufp != NULL && rval == eFieldValid_Valid) {
    if (ufp->label != NULL 
        && (StringICmp (ufp->label->str, "StructuredCommentPrefix") == 0
            || StringICmp (ufp->label->str, "StructuredCommentSuffix") == 0)) {
      /* skip suffix and prefix */
      ufp = ufp->next;
    } else {
      valid = IsStructuredCommentFieldValid (ufp, field_rule);
      if (valid == eFieldValid_OptionalFieldSkipped) {
        field_rule = field_rule->next;
      } else if (valid == eFieldValid_Valid) {
        field_rule = field_rule->next;
        ufp = ufp->next;
      } else {
        rval = valid;
      }
    }
  }
  if (valid == eFieldValid_Valid) {
    while (field_rule != NULL && !field_rule->required) {
      field_rule = field_rule->next;
    }
    if (field_rule != NULL) {
      rval = eFieldValid_MissingRequiredField;
    }
  }

  return rval;
}


NLM_EXTERN EFieldValid IsStructuredCommentValid (UserObjectPtr uop)
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
  return IsStructuredCommentValidForRule (uop, cr);
}


