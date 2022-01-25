/*   sequin9.c
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
* File Name:  sequin9.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   4/20/99
*
* $Revision: 6.2 $
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

#include "sequin.h"

/* general text parsing and conversion */

static ENUM_ALIST(from_fld_alist)
  {" ",                    0},
  {"Gene-CDS-Prot-RNA",    1},
  {"Source-Modifiers",     2},
  {"Feature Qualifiers",   3},
  {"COMMENT",              4},
  {"GenBank",              5},
  {"EMBL",                 6},
  {"Def Line",             7},
  {"Local ID",             8},
END_ENUM_ALIST

static ENUM_ALIST(to_fld_alist)
  {" ",                    0},
  {"Gene-CDS-Prot-RNA",    1},
  {"Source-Modifiers",     2},
  {"Feature Qualifiers",   3},
  {"COMMENT",              4},
END_ENUM_ALIST

static ENUM_ALIST(gene_cds_prot_rna_fld_alist)
  {" ",                    0},
  {"CDS comment",          1},
  {"CDS gene xref",        2},
  {"CDS protein xref",     3},
  {"Gene locus",           4},
  {"Gene description",     5},
  {"Gene allele",          6},
  {"Gene maploc",          7},
  {"Gene synonym",         8},
  {"Gene comment",         9},
  {"Protein name",        10},
  {"Protein description", 11},
  {"Protein E.C. number", 12},
  {"Protein activity",    13},
  {"Protein comment",     14},
  {"RNA name",            15},
  {"RNA comment",         16},
  {"RNA gene xref",       17},
END_ENUM_ALIST

static ENUM_ALIST(source_modifiers_fld_alist)
  {" ",                  0},
  {"Acronym",           19},
  {"Biotype",           14},
  {"Biovar",            13},
  {"Cell-line",        108},
  {"Cell-type",        109},
  {"Chemovar",          12},
  {"Chromosome",       101},
  {"Clone",            103},
  {"Clone-lib",        111},
  {"Common",            18},
  {"Common Name",      202},
  {"Country",          123},
  {"Cultivar",          10},
  {"Dev-stage",        112},
  {"Division",         204},
  {"Dosage",            20},
  {"Frequency",        113},
  {"Genotype",         106},
  {"Germline",         114},
  {"Group",             15},
  {"Haplotype",        105},
  {"Ins-seq-name",     121},
  {"Isolate",           17},
  {"Lab-host",         116},
  {"Lineage",          203},
  {"Map",              102},
  {"Natural-host",      21},
  {"Old Name",          54},
  {"OrgMod Note",       55},
  {"Pathovar",          11},
  {"Plasmid-name",     119},
  {"Plastid-name",     122},
  {"Pop-variant",      117},
  {"Rearranged",       115},
  {"Scientific Name",  201},
  {"Serogroup",          8},
  {"Serotype",           7},
  {"Serovar",            9},
  {"Sex",              107},
  {"Specimen-voucher",  23},
  {"Strain",             2},
  {"Sub-species",       22},
  {"Subclone",         104},
  {"Subgroup",          16},
  {"SubSource Note",   155},
  {"Substrain",          3},
  {"Subtype",            5},
  {"Tissue-lib",       118},
  {"Tissue-type",      110},
  {"Transposon-name",  120},
  {"Type",               4},
  {"Variety",            6},
END_ENUM_ALIST

static ENUM_ALIST(feature_qualifiers_fld_alist)
  {" ",               0},
  {"bound_moiety",    1},
  {"clone",           2},
  {"cons_splice",     3},
  {"direction",       4},
  {"frequency",       5},
  {"function",        6},
  {"label",           7},
  {"map",             8},
  {"mod_base",        9},
  {"note",           10},
  {"number",         11},
  {"organism",       12},
  {"PCR_conditions", 13},
  {"phenotype",      14},
  {"product",        15},
  {"replace",        16},
  {"rpt_family",     17},
  {"rpt_type",       18},
  {"rpt_unit",       19},
  {"standard_name",  20},
  {"usedin",         21},
END_ENUM_ALIST

static Uint1 SourceModListToOrgModType (UIEnum val) 

{
  if (val > 0 && val < 24) return (Uint1) val;
  if (val == 55) return 255;
  if (val == 54) return 254;
  return 0;
}

static Uint1 SourceModListToSubSourceType (UIEnum val) 

{
  if (val > 100 && val < 124) return (Uint1) val - 100;
  if (val == 155) return 255;
  return 0;
}

static Uint1 SourceModListToBioSourceField (UIEnum val) 

{
  if (val > 200 && val < 205) return (Uint1) val - 200;
  return 0;
}

#define NUM_SUBTARGET_POPUPS 10

typedef struct targetdata {
  PopuP              target;
  EnumFieldAssocPtr  alist;
  PopuP              subtarget [NUM_SUBTARGET_POPUPS];
  EnumFieldAssocPtr  alists [NUM_SUBTARGET_POPUPS];
  Int2               type;
  Int2               subtype;
} TargetData, PNTR TargetDataPtr;

static void ChangeTarget (PopuP p)

{
  Int2           i;
  TargetDataPtr  tdp;
  UIEnum         val;

  tdp = (TargetDataPtr) GetObjectExtra (p);
  if (tdp == NULL) return;
  if (GetEnumPopup (tdp->target, tdp->alist, &val) && val > 0) {
    for (i = 0; i < NUM_SUBTARGET_POPUPS; i++) {
      if (i != (Int2) val) {
        SafeHide (tdp->subtarget [i]);
      }
    }
    SafeShow (tdp->subtarget [(Int2) val]);
  }
}

static void CreateTransformTargetControl (GrouP h, TargetDataPtr tdp, EnumFieldAssocPtr alist)

{
  Int2   j;
  GrouP  p;
  GrouP  q;

  if (h == NULL || tdp == NULL || alist == NULL) return;

  p = HiddenGroup (h, 2, 0, NULL);

  tdp->target = PopupList (p, TRUE, ChangeTarget);
  SetObjectExtra (tdp->target, tdp, NULL);
  InitEnumPopup (tdp->target, alist, NULL);
  SetEnumPopup (tdp->target, alist, 0);

  tdp->alists [1] = gene_cds_prot_rna_fld_alist;
  tdp->alists [2] = source_modifiers_fld_alist;
  tdp->alists [3] = feature_qualifiers_fld_alist;

  q = HiddenGroup (p, 0, 0, NULL);
  for (j = 1; j < 4; j++) {
    tdp->subtarget [j] = PopupList (q, TRUE, NULL);
    SetObjectExtra (tdp->subtarget [j], tdp, NULL);
    InitEnumPopup (tdp->subtarget [j], tdp->alists [j], NULL);
    SetEnumPopup (tdp->subtarget [j], tdp->alists [j], 0);
    Hide (tdp->subtarget [j]);
  }
}

typedef struct stringdata {
  TexT     atleft;
  TexT     atright;
  GrouP    leftbehav;
  GrouP    rightbehav;
  CharPtr  leftstr;
  CharPtr  rightstr;
  Boolean  includeleft;
  Boolean  includeright;
  TexT     replaceby;
  CharPtr  replacestr;
} StringData, PNTR StringDataPtr;

static void CreateTransformStringControl (GrouP h, StringDataPtr sdp)

{
  GrouP  g;

  if (h == NULL || sdp == NULL) return;

  g = HiddenGroup (h, 3, 0, NULL);

  StaticPrompt (g, "Select text", 0, dialogTextHeight, programFont, 'l');
  sdp->leftbehav = HiddenGroup (g, 2, 0, NULL);
  RadioButton (sdp->leftbehav, "just after");
  RadioButton (sdp->leftbehav, "starting at");
  SetValue (sdp->leftbehav, 1);
  sdp->atleft = DialogText (g, "", 10, NULL);
  SetObjectExtra (sdp->atleft, sdp, NULL);

  StaticPrompt (g, "and", 0, dialogTextHeight, programFont, 'l');
  sdp->rightbehav = HiddenGroup (g, 2, 0, NULL);
  RadioButton (sdp->rightbehav, "just before");
  RadioButton (sdp->rightbehav, "ending with");
  SetValue (sdp->rightbehav, 1);
  sdp->atright = DialogText (g, "", 10, NULL);
  SetObjectExtra (sdp->atright, sdp, NULL);
}

typedef struct transformdata {
  FEATURE_FORM_BLOCK

  TargetData         fromtarget;
  TargetData         totarget;
  StringData         strings;
  ButtoN             accept;
  CharPtr            foundstr;
  Boolean            replaceOldAsked;
  Boolean            doReplaceAll;
  Int2               index;
} TransFormData, PNTR TransFormPtr;

static CharPtr SaveStringFromTextNoStripSpaces (TexT t)

{
  size_t   len;
  CharPtr  str;

  len = TextLength (t);
  if (len > 0) {
    str = (CharPtr) MemNew(len + 1);
    if (str != NULL) {
      GetTitle (t, str, len + 1);
      return str;
    } else {
      return NULL;
    }
  } else {
    return NULL;
  }
}

static void DoOneTransformText (Uint2 entityID, SeqEntryPtr sep, TransFormPtr tfp, MonitorPtr mon)

{
  /*
  BioSourcePtr  biop;
  BioseqSetPtr  bssp;
  Char          str [64];

  if (sep == NULL || tfp == NULL) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && (bssp->_class == 7 || bssp->_class == 13 ||
                         bssp->_class == 14 || bssp->_class == 15)) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        DoOneTransformText (entityID, sep, tfp, mon);
      }
      return;
    }
  }
  (tfp->index)++;
  if (mon != NULL) {
    sprintf (str, "Processing component %d", (int) tfp->index);
    MonitorStrValue (mon, str);
  }
  switch (tfp->fromtarget.type) {
    case 1 :
    case 2 :
    case 3 :
    case 4 :
      SeqEntryExplore (sep, (Pointer) tfp, RemoveAFeatureText);
      break;
    case 5 :
      SeqEntryToBioSource (sep, NULL, NULL, 0, &biop);
      RemoveASourceText (biop, tfp);
      break;
    case 6 :
    case 7 :
      SeqEntryToBioSource (sep, NULL, NULL, 0, &biop);
      RemoveASourceText (biop, tfp);
      break;
    default :
      break;
  }
  */
}

static void DoTransformProc (ButtoN b)

{
  MonitorPtr    mon;
  SeqEntryPtr   sep;
  TransFormPtr  tfp;
  UIEnum        val;

  tfp = (TransFormPtr) GetObjectExtra (b);
  if (tfp == NULL || tfp->input_entityID == 0) return;
  sep = GetTopSeqEntryForEntityID (tfp->input_entityID);
  if (sep == NULL) return;
  Hide (tfp->form);
  WatchCursor ();
  Update ();
  if (GetEnumPopup (tfp->fromtarget.target, tfp->fromtarget.alist, &val) && val > 0) {
    tfp->fromtarget.type = (Int2) val;
    if (GetEnumPopup (tfp->fromtarget.subtarget [tfp->fromtarget.type],
                      tfp->fromtarget.alists [tfp->fromtarget.type], &val) && val > 0) {
      tfp->fromtarget.subtype = (Int2) val;
      tfp->strings.leftstr = SaveStringFromTextNoStripSpaces (tfp->strings.atleft);
      tfp->strings.rightstr = SaveStringFromTextNoStripSpaces (tfp->strings.atright);
      tfp->strings.includeleft = (Boolean) (GetValue (tfp->strings.leftbehav) == 2);
      tfp->strings.includeright = (Boolean) (GetValue (tfp->strings.rightbehav) == 2);
      mon = MonitorStrNewEx ("Removing Text From String", 20, FALSE);
      tfp->index = 0;
      DoOneTransformText (tfp->input_entityID, sep, tfp, mon);
      MonitorFree (mon);
      tfp->strings.leftstr = MemFree (tfp->strings.leftstr);
      tfp->strings.rightstr = MemFree (tfp->strings.rightstr);
    }
  }
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (tfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, tfp->input_entityID, 0, 0);
  Remove (tfp->form);
}

