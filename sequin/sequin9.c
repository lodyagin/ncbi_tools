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
* $Revision: 6.12 $
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

#include <subutil.h>
#include <explore.h>
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

/* automatic defline generator */

typedef struct deffeats {
  SeqFeatPtr  sfp;
  SeqFeatPtr  gene;
  SeqFeatPtr  prot;
  CharPtr     genename;
  CharPtr     allelename;
  CharPtr     protname;
  Boolean     alreadyTrimmed;
  Uint2       entityID;
  Uint2       itemID;
  Uint2       subtype;
  Boolean     lastInString;
  Boolean     lastInGroup;
  Boolean     lastInType;
  Boolean     lastInPenultimate;
  Boolean     pseudo;
  Boolean     ignore;
  Boolean     suppressprefix;
  Int2        altSplices;
  Int2        numUnknown;
} DefFeatsData, PNTR DefFeatsPtr;

static Boolean GetMolBioFeatsGatherFunc (GatherContextPtr gcp, Boolean getGene, Boolean getMrna)

{
  DefFeatsPtr  dfp;
  RnaRefPtr    rrp;
  SeqFeatPtr   sfp;
  CharPtr      str;
  Uint1        type;
  ValNodePtr   PNTR vnpp;

  if (gcp == NULL || gcp->thisitem == NULL || gcp->userdata == NULL) return TRUE;
  if (gcp->thistype != OBJ_SEQFEAT) return TRUE;
  vnpp = (ValNodePtr PNTR) gcp->userdata;
  sfp = (SeqFeatPtr) gcp->thisitem;
  switch (sfp->data.choice) {
    case SEQFEAT_GENE :
      if (getGene) {
        dfp = MemNew (sizeof (DefFeatsData));
        if (dfp == NULL) return TRUE;
        dfp->entityID = gcp->entityID;
        dfp->itemID = gcp->itemID;
        dfp->sfp = sfp;
        dfp->subtype = FEATDEF_GENE;
        ValNodeAddPointer (vnpp, 0, (Pointer) dfp);
      }
      break;
    case SEQFEAT_CDREGION :
      dfp = MemNew (sizeof (DefFeatsData));
      if (dfp == NULL) return TRUE;
      dfp->entityID = gcp->entityID;
      dfp->itemID = gcp->itemID;
      dfp->sfp = sfp;
      dfp->subtype = FEATDEF_CDS;
      dfp->altSplices = 1;
      ValNodeAddPointer (vnpp, 0, (Pointer) dfp);
      break;
    case SEQFEAT_RNA :
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp == NULL) return TRUE;
      switch (rrp->type) {
        case 2 :
          if (getMrna) {
            dfp = MemNew (sizeof (DefFeatsData));
            if (dfp == NULL) return TRUE;
            dfp->entityID = gcp->entityID;
            dfp->itemID = gcp->itemID;
            dfp->sfp = sfp;
            dfp->subtype = FEATDEF_mRNA;
            ValNodeAddPointer (vnpp, 0, (Pointer) dfp);
          }
          break;
        case 3 :
          dfp = MemNew (sizeof (DefFeatsData));
          if (dfp == NULL) return TRUE;
          dfp->entityID = gcp->entityID;
          dfp->itemID = gcp->itemID;
          dfp->sfp = sfp;
          dfp->subtype = FEATDEF_tRNA;
          ValNodeAddPointer (vnpp, 0, (Pointer) dfp);
          break;
        case 4 :
          dfp = MemNew (sizeof (DefFeatsData));
          if (dfp == NULL) return TRUE;
          dfp->entityID = gcp->entityID;
          dfp->itemID = gcp->itemID;
          dfp->sfp = sfp;
          dfp->subtype = FEATDEF_rRNA;
          ValNodeAddPointer (vnpp, 0, (Pointer) dfp);
          break;
        case 255 :
          if (rrp->ext.choice == 1) {
            str = (CharPtr) rrp->ext.value.ptrvalue;
            if (StringICmp (str, "internal transcribed spacer 1") == 0 ||
                StringICmp (str, "internal transcribed spacer 2") == 0 ||
                StringICmp (str, "internal transcribed spacer 3") == 0 ||
                StringICmp (str, "internal transcribed spacer ITS1") == 0 ||
                StringICmp (str, "internal transcribed spacer ITS2") == 0 ||
                StringICmp (str, "internal transcribed spacer ITS3") == 0 ||
                StringICmp (str, "ITS1") == 0 ||
                StringICmp (str, "ITS2") == 0 ||
                StringICmp (str, "ITS3") == 0) {
              dfp = MemNew (sizeof (DefFeatsData));
              if (dfp == NULL) return TRUE;
              dfp->entityID = gcp->entityID;
              dfp->itemID = gcp->itemID;
              dfp->sfp = sfp;
              dfp->subtype = FEATDEF_otherRNA;
              ValNodeAddPointer (vnpp, 0, (Pointer) dfp);
            }
          }
          break;
        default :
          break;
      }
      break;
    case SEQFEAT_IMP :
      type = FindFeatDefType (sfp);
      if (type == FEATDEF_LTR || type == FEATDEF_exon) {
        dfp = MemNew (sizeof (DefFeatsData));
        if (dfp == NULL) return TRUE;
        dfp->entityID = gcp->entityID;
        dfp->itemID = gcp->itemID;
        dfp->sfp = sfp;
        dfp->subtype = type;
        ValNodeAddPointer (vnpp, 0, (Pointer) dfp);
      }
      break;
    default :
      break;
  }
  return TRUE;
}

static Boolean GetCDStRNArRNAGatherFunc (GatherContextPtr gcp)

{
  return GetMolBioFeatsGatherFunc (gcp, FALSE, FALSE);
}

static Boolean GetGeneCDStRNArRNAGatherFunc (GatherContextPtr gcp)

{
  return GetMolBioFeatsGatherFunc (gcp, TRUE, FALSE);
}

static Boolean GetGeneCDStRNArRNAmRNAGatherFunc (GatherContextPtr gcp)

{
  return GetMolBioFeatsGatherFunc (gcp, TRUE, TRUE);
}

static void LabelAModifier (CharPtr str, CharPtr text, Boolean labelMods)

{
  Char     ch;
  CharPtr  ptr;

  if (str == NULL || text == NULL) return;
  if (StringHasNoText (text)) {
    str [0] = '\0';
    text [0] = '\0';
  } else {
    if (labelMods) {
      ptr = str;
      while (*ptr != '\0') {
        ch = *ptr;
        *ptr = TO_LOWER (ch);
        ptr++;
      }
      StringCat (str, " ");
    } else {
      str [0] = '\0';
    }
    StringCat (str, text);
  }
}

static int LIBCALLBACK SortCDStRNArRNAByLocation (VoidPtr ptr1, VoidPtr ptr2)

{
  BioseqPtr    bsp1;
  BioseqPtr    bsp2;
  DefFeatsPtr  dfp1;
  DefFeatsPtr  dfp2;
  Int4         leftend1;
  Int4         leftend2;
  Int4         rightend1;
  Int4         rightend2;
  SeqFeatPtr   sfp1;
  SeqFeatPtr   sfp2;
  ValNodePtr   vnp1;
  ValNodePtr   vnp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      dfp1 = (DefFeatsPtr) vnp1->data.ptrvalue;
      dfp2 = (DefFeatsPtr) vnp2->data.ptrvalue;
      if (dfp1 != NULL && dfp2 != NULL) {
        sfp1 = dfp1->sfp;
        sfp2 = dfp2->sfp;
        if (sfp1 != NULL && sfp2 != NULL) {
          bsp1 = GetBioseqGivenSeqLoc (sfp1->location, dfp1->entityID);
          bsp2 = GetBioseqGivenSeqLoc (sfp2->location, dfp2->entityID);
          if (bsp1 != NULL && bsp2 != NULL) {
            leftend1 = GetOffsetInBioseq (sfp1->location, bsp1, SEQLOC_LEFT_END);
            leftend2 = GetOffsetInBioseq (sfp2->location, bsp2, SEQLOC_LEFT_END);
            rightend1 = GetOffsetInBioseq (sfp1->location, bsp1, SEQLOC_RIGHT_END);
            rightend2 = GetOffsetInBioseq (sfp2->location, bsp2, SEQLOC_RIGHT_END);
            if (leftend1 > leftend2) {
              return 1;
            } else if (leftend1 < leftend2) {
              return -1;
            } else if (sfp2->data.choice == SEQFEAT_GENE) {
              return 1;
            } else if (sfp1->data.choice == SEQFEAT_GENE) {
              return -1;
            } else if (rightend1 > rightend2) {
              return 1;
            } else if (rightend1 < rightend2) {
              return -1;
            } else {
              return 0;
            }
          }
        }
      }
    }
  }
  return 0;
}

static int LIBCALLBACK SortCDSAfterExons (VoidPtr ptr1, VoidPtr ptr2)

{
  DefFeatsPtr  dfp1;
  DefFeatsPtr  dfp2;
  ValNodePtr   vnp1;
  ValNodePtr   vnp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      dfp1 = (DefFeatsPtr) vnp1->data.ptrvalue;
      dfp2 = (DefFeatsPtr) vnp2->data.ptrvalue;
      if (dfp1 != NULL && dfp2 != NULL) {
        if (dfp1->subtype == FEATDEF_CDS && dfp2->subtype == FEATDEF_exon) {
          return 1;
        } else if (dfp1->subtype == FEATDEF_exon && dfp2->subtype == FEATDEF_CDS) {
          return -1;
        } else {
          /* return 0; */
          return SortCDStRNArRNAByLocation (ptr1, ptr2);
        }
      }
    }
  }
  return 0;
}

extern EnumFieldAssoc  orgmod_subtype_alist [];
extern EnumFieldAssoc  subsource_subtype_alist [];

static Int2  orgmod_rank [24];
static Int2  subsource_rank [24];

static Boolean StrainAlreadyInParentheses (CharPtr taxname, CharPtr strain)

{
  size_t   len;
  CharPtr  ptr;

  ptr = StringChr (taxname, '(');
  if (ptr == NULL) return FALSE;
  ptr++;
  len = StringLen (strain);
  if (StringNCmp (taxname, strain, len) != 0) return FALSE;
  ptr += len;
  if (*ptr != ')') return FALSE;
  return TRUE;
}

static void AddOrgModsToDef (ValNodePtr PNTR stringsPtr, BioSourcePtr biop, Boolean labelMods)

{
  EnumFieldAssocPtr  ap;
  OrgModPtr          mod;
  OrgNamePtr         onp;
  OrgRefPtr          orp;
  CharPtr            ptr;
  SubSourcePtr       ssp;
  Char               str [128];
  Char               text [64];

  if (stringsPtr != NULL && biop != NULL) {
    orp = biop->org;
    if (orp != NULL) {
      onp = orp->orgname;
      if (onp != NULL) {
        mod = onp->mod;
        while (mod != NULL) {
          if (mod->subtype < 24 && orgmod_rank [mod->subtype] > 0) {
            if (mod->subtype == 2 && StrainAlreadyInParentheses (orp->taxname, mod->subname)) {
              /* do not add strain if already parenthetical in organism name */
            } else {
              text [0] = '\0';
              str [0] = '\0';
              StringNCpy_0 (text, mod->subname, sizeof (text));
              for (ap = orgmod_subtype_alist; ap->name != NULL; ap++) {
                if (ap->value == mod->subtype) {
                  StringNCpy_0 (str, ap->name, sizeof (str));
                }
              }
              LabelAModifier (str, text, labelMods);
              if (! StringHasNoText (str)) {
                ValNodeCopyStr (stringsPtr, orgmod_rank [mod->subtype], str);
              }
            }
          }
          mod = mod->next;
        }
      }
    }

    ssp = biop->subtype;
    while (ssp != NULL) {
      if (ssp->subtype < 24 && subsource_rank [ssp->subtype] > 0) {
        text [0] = '\0';
        str [0] = '\0';
        StringNCpy_0 (text, ssp->name, sizeof (text));
        for (ap = subsource_subtype_alist; ap->name != NULL; ap++) {
          if (ap->value == ssp->subtype) {
            StringNCpy_0 (str, ap->name, sizeof (str));
            ptr = StringStr (str, "-name");
            if (ptr != NULL) {
              *ptr = '\0';
            }
            if (ssp->subtype == 23) { /* country */
              ptr = StringStr (text, ":");
              if (ptr != NULL) {
                *ptr = '\0';
              }
            }
          }
        }
        if (ssp->subtype == 23 && (! labelMods)) {
          StringCpy (str, "from");
        }
        LabelAModifier (str, text, labelMods || (ssp->subtype == 23));
        if (! StringHasNoText (str)) {
          ValNodeCopyStr (stringsPtr, subsource_rank [ssp->subtype], str);
        }
      }
      ssp = ssp->next;
    }
  }
}

static CharPtr GetExonNumber (GBQualPtr gbq)

{
  while (gbq != NULL) {
    if (StringICmp (gbq->qual, "number") == 0) {
      return gbq->val;
    }
  }

  return NULL;
}

static Boolean NextIsExon (ValNodePtr vnp)

{
  DefFeatsPtr  nextdfp;
  ValNodePtr   nextvnp;

  if (vnp == NULL) return FALSE;
  nextvnp = vnp->next;
  if (nextvnp == NULL) return FALSE;
  nextdfp = (DefFeatsPtr) nextvnp->data.ptrvalue;
  if (nextdfp != NULL && nextdfp->subtype == FEATDEF_exon) return TRUE;
  return FALSE;
}

static void FinishAutoDefProc (Uint2 entityID, SeqEntryPtr sep,
                               ValNodePtr head, BioseqPtr target,
                               SeqEntryPtr nsep, MolInfoPtr mip,
                               ValNodePtr strings, BioSourcePtr biop,
                               Int2 mitochloroflag)

{
  Int2          count;
  Boolean       ddbjstyle = FALSE;
  DefFeatsPtr   dfp;
  CharPtr       exonnumber;
  Int2          mitocount;
  DefFeatsPtr   nextdfp;
  CharPtr       ptr;
  RnaRefPtr     rrp;
  SeqFeatPtr    sfp;
  Char          str [380];
  tRNAPtr       trna;
  ValNodePtr    ttl;
  Char          text [64];
  ValNodePtr    vnp;

  if (GetAppParam ("SEQUIN", "PREFERENCES", "DATABASE", NULL, str, sizeof (str))) {
    if (StringICmp (str, "DDBJ") == 0) {
      ddbjstyle = TRUE;
    }
  }
  vnp = head;
  count = 0;
  mitocount = 0;
  while (vnp != NULL) {
    str [0] = '\0';
    dfp = (DefFeatsPtr) vnp->data.ptrvalue;
    if (dfp != NULL) {
      sfp = dfp->sfp;
      if (sfp != NULL || dfp->numUnknown > 0) {
        count++;
        mitocount++;
        /* FindGeneAndProtForCDS (entityID, sfp, &(dfp->gene), &(dfp->prot)); */
        /* StringCpy (str, "unknown"); */
        text [0] = '\0';
        if (dfp->subtype == FEATDEF_CDS && dfp->prot != NULL) {
          if (dfp->suppressprefix) {
            str [0] = '\0';
            if (dfp->sfp->partial) {
              StringCat (str, "partial cds");
            } else {
              StringCat (str, "complete cds");
            }
            if (dfp->altSplices > 1) {
              StringCat (str, ", alternatively spliced");
            }
          } else if (dfp->protname != NULL) {
            if (StringICmp (dfp->protname, "unknown") == 0 && (! StringHasNoText (dfp->genename))) {
              StringNCpy_0 (text, dfp->genename, sizeof (text));
              /* StringCat (str, "("); */
              StringCat (str, text);
              /* StringCat (str, ")"); */
            } else {
              StringNCpy_0 (str, dfp->protname, sizeof (str) - 100);
              if (dfp->genename != NULL) {
                StringNCpy_0 (text, dfp->genename, sizeof (text));
                if (! StringHasNoText (text)) {
                  StringCat (str, " (");
                  StringCat (str, text);
                  StringCat (str, ")");
                }
              }
            }
            if (dfp->lastInGroup || dfp->lastInType) {
              if (mip != NULL) {
                ptr = StringISearch (str, "precursor");
                if (ptr != NULL && StringICmp (ptr, "precursor") == 0) {
                  StringCat (str, ",");
                }
                if (mip->biomol == MOLECULE_TYPE_MRNA) {
                  StringCat (str, " mRNA");
                } else {
                  StringCat (str, " gene");
                }
              }
              if (count > 1) {
                StringCat (str, "s");
              }
              if (count < 2 && /* dfp->altSplices < 2 && */ (! StringHasNoText (dfp->allelename))) {
                StringNCpy_0 (text, dfp->allelename, sizeof (text));
                StringCat (str, ", ");
                StringCat (str, text);
                StringCat (str, " allele");
              }
              if (dfp->sfp->partial) {
                StringCat (str, ", partial cds");
              } else {
                StringCat (str, ", complete cds");
              }
              if (dfp->altSplices > 1) {
                StringCat (str, ", alternatively spliced");
              }
            }
          }
        } else if (dfp->subtype == FEATDEF_CDS && dfp->pseudo) {
          if (dfp->genename != NULL) {
            StringNCpy_0 (str, dfp->genename, sizeof (str) - 50);
          }
          if (dfp->lastInGroup || dfp->lastInType) {
            StringCat (str, " pseudogene");
            if (count > 1) {
              StringCat (str, "s");
            }
            if (count < 2 && (! StringHasNoText (dfp->allelename))) {
              StringNCpy_0 (text, dfp->allelename, sizeof (text));
              StringCat (str, ", ");
              StringCat (str, text);
              StringCat (str, " allele");
            }
            if (dfp->sfp->partial) {
              StringCat (str, ", partial sequence");
            } else {
              StringCat (str, ", complete sequence");
            }
          }
        } else if (dfp->subtype == FEATDEF_CDS && dfp->prot == NULL && dfp->genename != NULL) {
          StringNCpy_0 (str, dfp->genename, sizeof (str) - 50);
          if (dfp->lastInGroup || dfp->lastInType) {
            if (mip != NULL) {
              ptr = StringISearch (str, "precursor");
              if (ptr != NULL && StringICmp (ptr, "precursor") == 0) {
                StringCat (str, ",");
              }
              if (mip->biomol == MOLECULE_TYPE_MRNA) {
                StringCat (str, " mRNA");
              } else {
                StringCat (str, " gene");
              }
            }
            if (count > 1) {
              StringCat (str, "s");
            }
            if (count < 2 && (! StringHasNoText (dfp->allelename))) {
              StringNCpy_0 (text, dfp->allelename, sizeof (text));
              StringCat (str, ", ");
              StringCat (str, text);
              StringCat (str, " allele");
            }
            if (dfp->sfp->partial) {
              StringCat (str, ", partial cds");
            } else {
              StringCat (str, ", complete cds");
            }
          }
        } else if (dfp->subtype == FEATDEF_GENE) {
          if (dfp->genename != NULL) {
            StringNCpy_0 (str, dfp->genename, sizeof (str) - 50);
          }
          if (dfp->lastInGroup || dfp->lastInType) {
            StringCat (str, " gene");
            if (count > 1) {
              StringCat (str, "s");
            }
            if (count < 2 && (! StringHasNoText (dfp->allelename))) {
              StringNCpy_0 (text, dfp->allelename, sizeof (text));
              StringCat (str, ", ");
              StringCat (str, text);
              StringCat (str, " allele");
            }
            if (dfp->sfp->partial) {
              StringCat (str, ", partial sequence");
            } else {
              StringCat (str, ", complete sequence");
            }
          }
        } else if (dfp->subtype == FEATDEF_rRNA || dfp->subtype == FEATDEF_otherRNA) {
          rrp = (RnaRefPtr) dfp->sfp->data.value.ptrvalue;
          if (rrp != NULL) {
            if (rrp->ext.choice == 1) {
              StringNCpy_0 (str, (CharPtr) rrp->ext.value.ptrvalue, sizeof (str) - 50);
              if (dfp->subtype == FEATDEF_rRNA) {
                ptr = StringISearch (str, " rRNA");
                if (ptr != NULL) {
                  *ptr = '\0';
                }
                ptr = StringISearch (str, " ribosomal RNA");
                if (ptr != NULL) {
                  *ptr = '\0';
                }
                if (! StringHasNoText (str)) {
                  StringCat (str, " ribosomal RNA");
                }
              } else if (dfp->subtype == FEATDEF_otherRNA) {
              }
              if (dfp->genename != NULL) {
                StringNCpy_0 (text, dfp->genename, sizeof (text));
                if (! StringHasNoText (text)) {
                  StringCat (str, " (");
                  StringCat (str, text);
                  StringCat (str, ")");
                }
              }
              if (dfp->lastInString || dfp->lastInGroup || dfp->lastInType) {
                if (dfp->subtype == FEATDEF_rRNA && mip->biomol == MOLECULE_TYPE_GENOMIC) {
                  StringCat (str, " gene");
                  if (count > 1) {
                    StringCat (str, "s");
                  }
                }
              }
              if (dfp->lastInGroup || dfp->lastInType) {
                if (dfp->sfp->partial) {
                  StringCat (str, ", partial sequence");
                } else {
                  StringCat (str, ", complete sequence");
                }
              }
            }
          }
        } else if (dfp->subtype == FEATDEF_tRNA) {
          rrp = (RnaRefPtr) dfp->sfp->data.value.ptrvalue;
          if (rrp != NULL) {
            if (rrp->ext.choice == 2) {
              trna = rrp->ext.value.ptrvalue;
              if (trna != NULL) {
                if (FeatDefLabel (dfp->sfp, str, sizeof (str) - 2, OM_LABEL_CONTENT) > 0) {
                  if (dfp->genename != NULL) {
                    StringNCpy_0 (text, dfp->genename, sizeof (text));
                    if (! StringHasNoText (text)) {
                      StringCat (str, " (");
                      StringCat (str, text);
                      StringCat (str, ")");
                    }
                  }
                  if (dfp->lastInGroup || dfp->lastInType) {
                    StringCat (str, " gene");
                    if (count > 1) {
                      StringCat (str, "s");
                    }
                    if (dfp->sfp->partial) {
                      StringCat (str, ", partial sequence");
                    } else {
                      StringCat (str, ", complete sequence");
                    }
                  }
                }
              }
            }
          }
        } else if (dfp->subtype == FEATDEF_LTR) {
          if (! StringHasNoText (sfp->comment)) {
            StringNCpy_0 (str, sfp->comment, sizeof (str));
            ptr = StringISearch (str, " long terminal repeat");
            if (ptr != NULL) {
              *ptr = '\0';
            }
            ptr = StringISearch (str, " long terminal repeat");
            if (ptr != NULL) {
              *ptr = '\0';
            }
            if (! StringHasNoText (str)) {
              StringCat (str, " long terminal repeat");
            }
          } else {
            /* StringCpy (str, "uncharacterized"); */
            StringCpy (str, "long terminal repeat");
          }
          if (dfp->lastInGroup || dfp->lastInType) {
            if (dfp->sfp->partial) {
              StringCat (str, ", partial sequence");
            } else {
              StringCat (str, ", complete sequence");
            }
          }
        } else if (dfp->subtype == FEATDEF_exon && target != NULL) {
          if (dfp->protname != NULL) {
            str [0] = '\0';
            if (! dfp->suppressprefix) {
              StringNCpy_0 (str, dfp->protname, sizeof (str) - 100);
              if (dfp->genename != NULL) {
                StringNCpy_0 (text, dfp->genename, sizeof (text));
                if (! StringHasNoText (text)) {
                  StringCat (str, " (");
                  StringCat (str, text);
                  StringCat (str, ") gene,");
                  if (/* count < 2 && */ /* dfp->altSplices < 2 && */ (! StringHasNoText (dfp->allelename))) {
                    StringNCpy_0 (text, dfp->allelename, sizeof (text));
                    StringCat (str, " ");
                    StringCat (str, text);
                    StringCat (str, " allele,");
                  }
                }
              }
            }
            if ((! StringHasNoText (str)) || dfp->suppressprefix) {
              exonnumber = GetExonNumber (sfp->qual);
              if (! dfp->suppressprefix) {
                if (exonnumber == NULL) {
                  if (StringStr (sfp->comment, "exon") != NULL && (! NextIsExon (vnp))) {
                  } else {
                    StringCat (str, " exon");
                    if (NextIsExon (vnp)) {
                      StringCat (str, "s");
                    }
                  }
                } else {
                  StringCat (str, " exon");
                  if (NextIsExon (vnp)) {
                    StringCat (str, "s");
                  }
                }
              }
              if (exonnumber != NULL) {
                if (! dfp->suppressprefix) {
                  StringCat (str, " ");
                }
                StringCat (str, exonnumber);
              } else {
                if (! StringHasNoText (sfp->comment)) {
                  if (! dfp->suppressprefix) {
                    StringCat (str, " ");
                  }
                  StringCat (str, sfp->comment);
                }
              }
            }
          } else if (dfp->genename != NULL) {
            str [0] = '\0';
            if (! dfp->suppressprefix) {
              StringNCpy_0 (str, dfp->genename, sizeof (str));
              StringCat (str, " gene,");
            }
            if ((! StringHasNoText (str)) || dfp->suppressprefix) {
              exonnumber = GetExonNumber (sfp->qual);
              if (! dfp->suppressprefix) {
                if (exonnumber == NULL) {
                  if (StringStr (sfp->comment, "exon") != NULL && (! NextIsExon (vnp))) {
                  } else {
                    StringCat (str, " exon");
                    if (NextIsExon (vnp)) {
                      StringCat (str, "s");
                    }
                  }
                } else {
                  StringCat (str, " exon");
                  if (NextIsExon (vnp)) {
                    StringCat (str, "s");
                  }
                }
              }
              if (exonnumber != NULL) {
                if (! dfp->suppressprefix) {
                  StringCat (str, " ");
                }
                StringCat (str, exonnumber);
              } else {
                if (! StringHasNoText (sfp->comment)) {
                  if (! dfp->suppressprefix) {
                    StringCat (str, " ");
                  }
                  StringCat (str, sfp->comment);
                }
              }
            }
          } else {
            StringCpy (str, "uncharacterized exon");
          }
          /*
          if (dfp->lastInGroup || dfp->lastInType) {
            if (dfp->sfp->partial) {
              StringCat (str, ", partial sequence");
            } else {
              StringCat (str, ", complete sequence");
            }
          }
          */
        } else if (dfp->numUnknown > 0) {
          if (mip != NULL && mip->biomol == MOLECULE_TYPE_MRNA) {
            StringCat (str, "unknown mRNA");
          } else {
            StringCat (str, "unknown gene");
          }
          if (dfp->numUnknown > 1) {
            StringCat (str, "s");
          }
        }
      }
    }
    vnp = vnp->next;
    if (! StringHasNoText (str)) {
      if (vnp == NULL) {
        if (biop != NULL) {
          if (biop->genome == GENOME_mitochondrion) {
            if (mitocount > 1) {
              StringCat (str, "; mitochondrial genes for mitochondrial products");
            } else {
              StringCat (str, "; mitochondrial gene for mitochondrial product");
            }
          } else if (biop->genome == GENOME_chloroplast) {
            if (mitocount > 1) {
              StringCat (str, "; chloroplast genes for chloroplast products");
            } else {
              StringCat (str, "; chloroplast gene for chloroplast product");
            }
          } else if (mitochloroflag > 0) {
            switch (mitochloroflag) {
              case 1 :
                if (mitocount > 1) {
                  StringCat (str, "; nuclear genes for mitochondrial products");
                } else {
                  StringCat (str, "; nuclear gene for mitochondrial product");
                }
                break;
              case 2 :
                if (mitocount > 1) {
                   StringCat (str, "; nuclear genes for chloroplast products");
                 } else {
                   StringCat (str, "; nuclear gene for chloroplast product");
                 }
                break;
              case 3 :
                if (mitocount > 1) {
                  StringCat (str, "; mitochondrial genes for mitochondrial products");
                } else {
                  StringCat (str, "; mitochondrial gene for mitochondrial product");
                }
                break;
              case 4 :
                if (mitocount > 1) {
                   StringCat (str, "; chloroplast genes for chloroplast products");
                 } else {
                   StringCat (str, "; chloroplast gene for chloroplast product");
                 }
                break;
              default :
                break;
            }
          }
        }
        StringCat (str, ".");
      } else if (vnp->next == NULL) {
        nextdfp = (DefFeatsPtr) vnp->data.ptrvalue;
        if (dfp->lastInPenultimate) {
          if ((dfp->subtype == FEATDEF_rRNA && nextdfp->subtype == FEATDEF_otherRNA &&
               dfp->sfp->partial == nextdfp->sfp->partial) ||
              (dfp->subtype == FEATDEF_otherRNA && nextdfp->subtype == FEATDEF_rRNA &&
               dfp->sfp->partial == nextdfp->sfp->partial)) {
            StringCat (str, ", and");
          } else {
            StringCat (str, "; and");
          }
        } else if (dfp->lastInType || dfp->lastInGroup) {
          if (count > 1) {
            StringCat (str, ", and");
          } else {
            StringCat (str, " and");
          }
        } else if (nextdfp->lastInType || nextdfp->lastInGroup) {
          if (count > 1) {
            StringCat (str, ", and");
          } else {
            StringCat (str, " and");
          }
        } else {
          if (count > 1) {
            StringCat (str, ", and");
          } else {
            StringCat (str, " and");
          }
        }
      } else {
        nextdfp = (DefFeatsPtr) vnp->data.ptrvalue;
        if (dfp->lastInPenultimate) {
          if ((dfp->subtype == FEATDEF_rRNA && nextdfp->subtype == FEATDEF_otherRNA &&
               dfp->sfp->partial == nextdfp->sfp->partial) ||
              (dfp->subtype == FEATDEF_otherRNA && nextdfp->subtype == FEATDEF_rRNA &&
               dfp->sfp->partial == nextdfp->sfp->partial)) {
            StringCat (str, ", and");
          } else {
            StringCat (str, "; and");
          }
        } else if (dfp->lastInType || dfp->lastInGroup) {
          StringCat (str, ";");
        } else if (nextdfp->lastInType || nextdfp->lastInGroup) {
          if (count > 1) {
            StringCat (str, ", and");
          } else {
            StringCat (str, " and");
          }
        } else {
          StringCat (str, ",");
        }
      }
      ValNodeCopyStr (&strings, 0, str);
    }
    if (dfp->lastInString || dfp->lastInGroup || dfp->lastInType) {
      count = 0;
    }
  }

  ptr = MergeValNodeStrings (strings, FALSE);
  if (nsep != NULL) {
    ttl = SeqEntryGetSeqDescr (nsep, Seq_descr_title, NULL);
    if (ttl == NULL) {
      ttl = SeqEntryGetSeqDescr (sep, Seq_descr_title, NULL);
    }
    if (ttl == NULL) {
      ttl = CreateNewDescriptor (nsep, Seq_descr_title);
    }
    if (ttl != NULL) {
      MemFree (ttl->data.ptrvalue);
      ttl->data.ptrvalue = ptr;
      ptr = NULL;
    }
  }
  MemFree (ptr);
  ValNodeFreeData (strings);
  ValNodeFreeData (head);
}

static void CombineProteinNames (DefFeatsPtr dfp1, DefFeatsPtr dfp2)

{
  Int2     i, j, lastspace, lastdash, lastcomma;
  size_t   len1, len2;
  CharPtr  str1 = NULL, str2 = NULL;

  if (dfp1 == NULL || dfp2 == NULL) return;
  if (dfp1->protname == NULL && dfp2->protname == NULL) return;
  len1 = StringLen (dfp1->protname);
  len2 = StringLen (dfp2->protname);
  if (len1 < 1 || len2 < 1) return;
  i = 0;
  j = 0;
  lastspace = 0;
  lastdash = 0;
  lastcomma = 0;
  while (i < len1 && j < len2 && dfp1->protname [i] == dfp2->protname [j]) {
    if (dfp1->protname [i] == ' ') {
      lastspace = i;
    }
    if (dfp1->protname [i] == '-') {
      lastdash = i;
    }
    if (dfp1->protname [i] == ',') {
      lastcomma = i;
    }
    i++;
    j++;
  }
  str1 = StringSave (dfp1->protname);
  if (str1 != NULL) {
    str1 [i] = '\0';
    if (! dfp1->alreadyTrimmed) {
      if (lastcomma > 0) {
        str1 [lastcomma] = '\0';
        dfp1->alreadyTrimmed = TRUE;
      } else if (lastdash > 0) {
        str1 [lastdash] = '\0';
        dfp1->alreadyTrimmed = TRUE;
      } else if (lastspace > 0) {
        str1 [lastspace] = '\0';
        dfp1->alreadyTrimmed = TRUE;
      }
    }
  }

  i = len1;
  j = len2;
  while (i > 0 && j > 0 && dfp1->protname [i - 1] == dfp2->protname [j - 1]) {
    i--;
    j--;
  }
  str2 = StringSave (dfp1->protname + i);
  TrimSpacesAroundString (str1);
  TrimSpacesAroundString (str2);
  len1 = StringLen (str1);
  len2 = StringLen (str2);
  if (len1 > len2) {
    dfp1->protname = str1;
    MemFree (str2);
  } else {
    dfp1->protname = str2;
    MemFree (str1);
  }
}

static Boolean AreAltSpliceGenes (DefFeatsPtr dfp1, DefFeatsPtr dfp2)

{
  Int2        comp;
  SeqFeatPtr  sfp1, sfp2;

  if (dfp1 == NULL || dfp2 == NULL) return FALSE;
  if (dfp1->prot == NULL || dfp2->prot == NULL) return FALSE;
  if (dfp1->pseudo || dfp2->pseudo) return FALSE;
  sfp1 = dfp1->sfp;
  sfp2 = dfp2->sfp;
  if (sfp1 == NULL || sfp2 == NULL) return FALSE;
  if (sfp1->partial != sfp2->partial) return FALSE;
  if (SeqLocStrand (sfp1->location) != SeqLocStrand (sfp2->location)) return FALSE;
  comp = SeqLocCompare (sfp1->location, sfp2->location);
  if (comp == SLC_NO_MATCH) return FALSE;
  if (dfp1->genename == NULL || dfp2->genename == NULL) return FALSE;
  if (StringICmp (dfp1->genename, dfp2->genename) != 0) return FALSE;
  CombineProteinNames (dfp1, dfp2);
  return TRUE;
}

static void MergeAltSpliceCDSs (ValNodePtr head)

{
  DefFeatsPtr  dfp1, dfp2;
  ValNodePtr   nextvnp;
  ValNodePtr   PNTR prevvnp;
  ValNodePtr   vnp1, vnp2;

  vnp1 = head;
  while (vnp1 != NULL) {
    dfp1 = (DefFeatsPtr) vnp1->data.ptrvalue;
    if (dfp1 != NULL && dfp1->subtype == FEATDEF_CDS) {
      vnp2 = vnp1->next;
      prevvnp = &(vnp1->next);
      while (vnp2 != NULL) {
        nextvnp = vnp2->next;
        dfp2 = (DefFeatsPtr) vnp2->data.ptrvalue;
        if (dfp2 != NULL && dfp2->subtype == FEATDEF_CDS) {
          if (AreAltSpliceGenes (dfp1, dfp2)) {
            (dfp1->altSplices)++;
            *prevvnp = vnp2->next;
            vnp2->next = NULL;
            ValNodeFreeData (vnp2);
          } else {
            prevvnp = (ValNodePtr PNTR) &(vnp2->next);
          }
        } else {
          prevvnp = (ValNodePtr PNTR) &(vnp2->next);
        }
        vnp2 = nextvnp;
      }
    }
    vnp1 = vnp1->next;
  }
}

static void AutoDefProc (Uint2 entityID, SeqEntryPtr sep, Boolean addMods,
                         Boolean labelMods, Int2 maxMods, Boolean leaveInParen,
                         BioseqPtr target, BioseqPtr seg, ValNodePtr nonUniqueOrgs,
                         Int2 mitochloroflag, BioseqPtr parent)

{
  Char            allele [256];
  BioseqContextPtr  bcp;
  BioSourcePtr    biop;
  BioseqPtr       bsp;
  BioseqSetPtr    bssp;
  SeqFeatPtr      cds;
  Boolean         change;
  CdRegionPtr     crp;
  DefFeatsData    df;
  DefFeatsPtr     dfp;
  DefFeatsPtr     dfpx;
  SeqMgrFeatContext  fcontext;
  GBQualPtr       gbq;
  SeqFeatPtr      gene;
  Int2            group;
  GeneRefPtr      grp;
  GatherScope     gs;
  ValNodePtr      head;
  SeqLocPtr       lastslp;
  Int4            left;
  size_t          lenallele;
  size_t          lenlocus;
  MolInfoPtr      mip;
  DefFeatsPtr     nextdfp;
  SeqLocPtr       nextslp;
  ValNodePtr      nextvnp;
  SeqEntryPtr     nsep;
  Int2            numUnknown;
  BioseqPtr       part;
  DefFeatsPtr     penult;
  ValNodePtr      PNTR prevvnp;
  ProtRefPtr      prp;
  CharPtr         ptr;
  Int4            right;
  ValNodePtr      sdp;
  SeqFeatPtr      sfp;
  SeqIdPtr        sip;
  SeqLocPtr       slp;
  Char            str [128];
  ValNodePtr      ttl;
  ValNodePtr      strings;
  ValNode         vn;
  ValNodePtr      vnp;
  ValNodePtr      vnpx;
  SeqFeatXrefPtr  xref;

  if (sep == NULL) return;
  if (target == NULL && seg == NULL) {
    if (IS_Bioseq_set (sep)) {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      if (bssp == NULL) return;
      if (bssp->_class == 7 || bssp->_class == 13 ||
          bssp->_class == 14 || bssp->_class == 15) {
        for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
          AutoDefProc (entityID, sep, addMods, labelMods, maxMods, leaveInParen, NULL, NULL, nonUniqueOrgs, mitochloroflag, NULL);
        }
        return;
      }
    }

    nsep = FindNucSeqEntry (sep);
    if (nsep != NULL) {
      bsp = (BioseqPtr) nsep->data.ptrvalue;
      if (bsp != NULL && bsp->repr == Seq_repr_seg && bsp->seq_ext != NULL && bsp->seq_ext_type == 1) {
        vn.choice = SEQLOC_MIX;
        vn.next = NULL;
        vn.data.ptrvalue = bsp->seq_ext;
        slp = SeqLocFindNext (&vn, NULL);
        while (slp != NULL) {
          nextslp = SeqLocFindNext (&vn, slp);
          sip = SeqLocId (slp);
          if (sip != NULL) {
            part = BioseqFind (sip);
            if (part != NULL) {
              AutoDefProc (entityID, sep, addMods, labelMods, maxMods, leaveInParen, part, NULL, nonUniqueOrgs, mitochloroflag, bsp);
            }
          }
          slp = nextslp;
        }
        AutoDefProc (entityID, sep, addMods, labelMods, maxMods, leaveInParen, NULL, bsp, nonUniqueOrgs, mitochloroflag, bsp);
        return;
      }
    }
  } else if (target != NULL) {
    nsep = SeqMgrGetSeqEntryForData (target);
  } else {
    nsep = SeqMgrGetSeqEntryForData (seg);
  }

  biop = NULL;
  strings = NULL;
  str [0] = '\0';

  SeqEntryToBioSource (sep, NULL, str, sizeof (str) - 1, &biop);
  if (! leaveInParen) {
    ptr = StringStr (str, "(");
    if (ptr != NULL) {
      *ptr = '\0';
    }
  }
  TrimSpacesAroundString (str);
  if (StringICmp (str, "Human immunodeficiency virus type 1") == 0) {
    StringCpy (str, "HIV-1");
  } else if (StringICmp (str, "Human immunodeficiency virus type 2") == 0) {
    StringCpy (str, "HIV-2");
  }
  str [0] = TO_UPPER (str [0]);
  ValNodeCopyStr (&strings, 0, str);

  mip = NULL;
  if (nsep != NULL) {
    bsp = (BioseqPtr) nsep->data.ptrvalue;
    bcp = BioseqContextNew (bsp);
    sdp = BioseqContextGetSeqDescr (bcp, Seq_descr_molinfo, NULL, NULL);
    BioseqContextFree (bcp);
    if (sdp != NULL) {
      mip = (MolInfoPtr) sdp->data.ptrvalue;
      if (mip != NULL) {
        if (mip->tech == MI_TECH_htgs_0 ||
            mip->tech == MI_TECH_htgs_1 ||
            mip->tech == MI_TECH_htgs_2 ||
            mip->tech == MI_TECH_est ||
            mip->tech == MI_TECH_sts ||
            mip->tech == MI_TECH_survey) {
          ttl = ValNodeExtract (&(bsp->descr), Seq_descr_title);
          if (ttl != NULL) {
            ttl = ValNodeFreeData (ttl);
          }
          return;
        }
      }
    }
  }

  /*
  nsep = FindNucSeqEntry (sep);
  if (nsep != NULL) {
    sdp = SeqEntryGetSeqDescr (nsep, Seq_descr_molinfo, NULL);
    if (sdp == NULL) {
      sdp = SeqEntryGetSeqDescr (sep, Seq_descr_molinfo, NULL);
    }
    if (sdp != NULL) {
      mip = (MolInfoPtr) sdp->data.ptrvalue;
    }
  }
  */

  if (nonUniqueOrgs != NULL) {
    for (vnp = nonUniqueOrgs; vnp != NULL; vnp = vnp->next) {
      if (StringICmp ((CharPtr) vnp->data.ptrvalue, str) == 0) {
        if (vnp->choice == 1) {
          addMods = FALSE; /* if only one organism in record, already unique defline */
        }
      }
    }
  }
  if (addMods) {
    AddOrgModsToDef (&strings, biop, labelMods);
    strings = SortValNode (strings, SortByVnpChoice);
    vnp = strings;
    while (vnp != NULL && maxMods > 0) {
      maxMods--;
      vnp = vnp->next;
    }
    if (vnp != NULL) {
      vnp->next = ValNodeFreeData (vnp->next);
    }
  }

  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = TRUE;
  MemSet ((Pointer) (gs.ignore), (int)(TRUE), (size_t) (OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.scope = sep;
  gs.target = NULL;
  if (target != NULL) {
    slp = ValNodeNew (NULL);
    slp->choice = SEQLOC_WHOLE;
    sip = SeqIdStripLocus (SeqIdDup (SeqIdFindBest (target->id, 0)));
    slp->data.ptrvalue = sip;
    gs.target = slp;
  }
  head = NULL;
  GatherEntity (entityID, (Pointer) (&head), GetGeneCDStRNArRNAGatherFunc, &gs);
  gs.target = SeqLocFree (gs.target);
  head = SortValNode (head, SortCDStRNArRNAByLocation);

  if (head == NULL && mip != NULL && mip->tech == MI_TECH_survey) {
    ValNodeCopyStr (&strings, 0, ", genome survey sequence.");
  }

  numUnknown = 0;
  vnp = head;
  while (vnp != NULL) {
    dfp = (DefFeatsPtr) vnp->data.ptrvalue;
    if (dfp != NULL) {
      sfp = dfp->sfp;
      if (sfp != NULL) {
        FindGeneAndProtForCDS (entityID, sfp, &(dfp->gene), &(dfp->prot));
        grp = SeqMgrGetGeneXref (sfp);
        if (SeqMgrGeneIsSuppressed (grp)) {
          dfp->gene = NULL;
        }
        dfp->pseudo = FALSE;
        dfp->genename = NULL;
        dfp->allelename = NULL;
        grp = NULL;
        if (dfp->gene != NULL) {
          grp = (GeneRefPtr) dfp->gene->data.value.ptrvalue;
          if (grp != NULL && grp->pseudo) {
            dfp->pseudo = TRUE;
          }
        }
        xref = sfp->xref;
        while (xref != NULL && xref->data.choice != SEQFEAT_GENE) {
          xref = xref->next;
        }
        if (xref != NULL) {
          grp = (GeneRefPtr) xref->data.value.ptrvalue;
        }
        if (grp != NULL) {
          dfp->genename = (CharPtr) grp->locus;
          if ((! StringHasNoText (grp->locus)) && (! StringHasNoText (grp->allele))) {
            lenallele = StringLen (grp->allele);
            lenlocus = StringLen (grp->locus);
            if (lenallele > lenlocus && StringNICmp (grp->locus, grp->allele, lenlocus) == 0) {
              sprintf (allele, "%s", grp->allele);
            } else if (StringNCmp (grp->allele, "-", 1) == 0) {
              sprintf (allele, "%s%s", grp->locus, grp->allele);
            } else {
              sprintf (allele, "%s-%s", grp->locus, grp->allele);
            }
            dfp->allelename = StringSave (allele);
          }
          if (grp->pseudo) {
            dfp->pseudo = TRUE;
          }
        }
        if (dfp->subtype == FEATDEF_CDS) {
          if (target != NULL) {
            lastslp = NULL;
            slp = SeqLocFindNext (sfp->location, NULL);
            while (slp != NULL) {
              lastslp = slp;
              slp = SeqLocFindNext (sfp->location, slp);
            }
            if (lastslp != NULL && nsep != NULL) {
              bsp = (BioseqPtr) nsep->data.ptrvalue;
              if (GetBioseqGivenSeqLoc (lastslp, entityID) != bsp) {
                dfp->ignore = TRUE;
              }
            }
          }
          crp = (CdRegionPtr) sfp->data.value.ptrvalue;
          if (crp != NULL) {
            if (crp->orf) {
              dfp->ignore = TRUE;
            }
          }
          gbq = sfp->qual;
          while (gbq != NULL) {
            if (StringICmp (gbq->qual, "pseudo") == 0) {
              dfp->pseudo = TRUE;
            }
            gbq = gbq->next;
          }
          if (dfp->pseudo) {
          } else if (dfp->prot == NULL) {
            if (dfp->gene == NULL) {
              dfp->ignore = TRUE;
            }
          } else {
            prp = (ProtRefPtr) dfp->prot->data.value.ptrvalue;
            if (prp != NULL) {
              if (prp->name == NULL || StringHasNoText ((CharPtr) prp->name->data.ptrvalue)) {
                if (prp->desc == NULL || StringHasNoText (prp->desc)) {
                  dfp->ignore = TRUE;
                } else {
                  dfp->protname = prp->desc;
                }
              } else {
                dfp->protname = (CharPtr) prp->name->data.ptrvalue;
              }
              if (! dfp->ignore) {
                if (StringICmp (dfp->protname, "unknown") == 0) {
                  if (StringHasNoText (dfp->genename)) {
                    numUnknown++;
                    dfp->ignore = TRUE;
                  }
                }
              }
            }
          }
        } else if (dfp->subtype == FEATDEF_exon && target == NULL) {
          dfp->ignore = TRUE;
        }
      }
    }
    vnp = vnp->next;
  }

  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    dfp = (DefFeatsPtr) vnp->data.ptrvalue;
    if (dfp != NULL) {
      sfp = dfp->sfp;
      if (sfp != NULL) {
        if (sfp->data.choice == SEQFEAT_GENE) {
          for (vnpx = vnp->next; vnpx != NULL; vnpx = vnpx->next) {
            dfpx = (DefFeatsPtr) vnpx->data.ptrvalue;
            if (dfpx != NULL) {
              if (sfp == dfpx->gene) {
                dfp->ignore = TRUE;
              }
            }
          }
        }
      }
    }
  }

  vnp = head;
  prevvnp = &head;
  while (vnp != NULL) {
    nextvnp = vnp->next;
    dfp = (DefFeatsPtr) vnp->data.ptrvalue;
    if (dfp->ignore) {
      *prevvnp = vnp->next;
      vnp->next = NULL;
      ValNodeFreeData (vnp);
    } else {
      prevvnp = (ValNodePtr PNTR) &(vnp->next);
    }
    vnp = nextvnp;
  }

  MergeAltSpliceCDSs (head);

  if (target != NULL) {
    head = SortValNode (head, SortCDSAfterExons);
  }

  if (numUnknown > 0) {
    dfp = (DefFeatsPtr) MemNew (sizeof (DefFeatsData));
    if (dfp != NULL) {
      dfp->entityID = entityID;
      dfp->subtype = 0;
      dfp->numUnknown = numUnknown;
      ValNodeAddPointer (&head, 0, (Pointer) dfp);
    }
  }

  if (target != NULL) {
    for (vnp = head; vnp != NULL; vnp = vnp->next) {
      dfp = (DefFeatsPtr) vnp->data.ptrvalue;
      if (dfp != NULL && dfp->subtype == FEATDEF_exon && dfp->gene != NULL) {
        gene = SeqMgrGetDesiredFeature (entityID, NULL, 0, 0, dfp->gene, &fcontext);
        if (gene == dfp->gene && (! gene->pseudo)) {
          left = fcontext.left;
          right = fcontext.right;
          cds = SeqMgrGetNextFeature (parent, NULL, SEQFEAT_CDREGION, 0, &fcontext);
          while (cds != NULL) {
            if (fcontext.left >= left && fcontext.right <= right && (! cds->pseudo)) {
              if (dfp->protname == NULL) {
                dfp->protname = fcontext.label; /* points to stable string */
              } else {
                df.protname = fcontext.label;
                CombineProteinNames (dfp, &df);
              }
            }
            cds = SeqMgrGetNextFeature (parent, cds, SEQFEAT_CDREGION, 0, &fcontext);
          }
        }
      }
    }
    vnp = head;
    while (vnp != NULL) {
      dfp = (DefFeatsPtr) vnp->data.ptrvalue;
      vnp = vnp->next;
      if (dfp != NULL && dfp->subtype == FEATDEF_exon && vnp != NULL) {
        nextdfp = (DefFeatsPtr) vnp->data.ptrvalue;
        if (nextdfp != NULL && (nextdfp->subtype == FEATDEF_exon || nextdfp->subtype == FEATDEF_CDS)) {
          if (StringCmp (dfp->genename, nextdfp->genename) == 0 &&
              StringCmp (dfp->protname, nextdfp->protname) == 0) {
            nextdfp->suppressprefix = TRUE;
          }
        }
      }
    }
  }

  vnp = head;
  group = 0;
  penult = NULL;
  while (vnp != NULL) {
    dfp = (DefFeatsPtr) vnp->data.ptrvalue;
    vnp->choice = (Uint1) group;
    vnp = vnp->next;
    if (vnp != NULL) {
      change = FALSE;
      nextdfp = (DefFeatsPtr) vnp->data.ptrvalue;
      if ((dfp->subtype == FEATDEF_rRNA && nextdfp->subtype == FEATDEF_otherRNA) ||
          (dfp->subtype == FEATDEF_otherRNA && nextdfp->subtype == FEATDEF_rRNA)) {
        dfp->lastInString = TRUE;
        if (dfp->sfp->partial != nextdfp->sfp->partial) {
          dfp->lastInGroup = TRUE;
        }
        change = TRUE;
      } else if (dfp->subtype != nextdfp->subtype) {
        if (dfp->subtype == FEATDEF_exon && nextdfp->subtype == FEATDEF_CDS && nextdfp->suppressprefix) {
          /* no separator between exons and appropriate cds */
        } else {
          dfp->lastInType = TRUE;
          change = TRUE;
        }
      } else if (dfp->sfp->partial != nextdfp->sfp->partial) {
        dfp->lastInGroup = TRUE;
        change = TRUE;
      } else if (dfp->pseudo != nextdfp->pseudo) {
        dfp->lastInGroup = TRUE;
        change = TRUE;
      } else if (dfp->allelename != NULL || nextdfp->allelename != NULL) {
        dfp->lastInGroup = TRUE;
        change = TRUE;
      } else if (dfp->altSplices > 1 || nextdfp->altSplices > 1) {
        dfp->lastInGroup = TRUE;
        change = TRUE;
      }
      if (change) {
        group++;
        penult = dfp;
      }
    } else {
      dfp->lastInString = TRUE;
      dfp->lastInGroup = TRUE;
      dfp->lastInType = TRUE;
      group++;
    }
  }
  if (penult != NULL) {
    penult->lastInPenultimate = TRUE;
  }

  FinishAutoDefProc (entityID, sep, head, target, nsep, mip, strings, biop, mitochloroflag);
}

static CharPtr sourceModRankList [] = {
  "Strain", "Isolate", "Clone", "Type", "Cultivar", "Haplotype",
  "Substrain", "Subclone", "Subtype", "Serotype", "Serogroup", "Serovar",
  "Variety", "Pathovar", "Chemovar", "Biovar", "Biotype", "Group", "Subgroup",
  "Cell-line", "Cell-type", "Tissue-type", "Clone-lib", "Tissue-lib", "Dev-stage",
  "Lab-host", "Pop-variant", "Frequency", "Germline", "Rearranged",
  "Chromosome", "Map", "Genotype", "Sex", "Plasmid-name", "Transposon-name",
  "Ins-seq-name", "Plastid-name", "Country", "Old Name", "Common", "Acronym",
  "Dosage", "Natural-host", "Sub-species", "Specimen-voucher",
  NULL
};

typedef struct deflineform {
  FORM_MESSAGE_BLOCK
  SeqEntryPtr    sep;
  ButtoN         addLabels;
  GrouP          customGrp;
  GrouP          sourceListGrp;
  ButtoN         PNTR sourceBoxList;
  PopuP          modLimit;
  ButtoN         onlyModifyTargeted;
  BioseqPtr      target;
  ButtoN         leaveInParentheses;
  GrouP          nucformitoorchloro;
  Boolean        smartMods;
} DeflineForm, PNTR DeflineFormPtr;

static int LIBCALLBACK SortByName (VoidPtr ptr1, VoidPtr ptr2)

{
  CharPtr     str1;
  CharPtr     str2;
  ValNodePtr  vnp1;
  ValNodePtr  vnp2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      str1 = (CharPtr) vnp1->data.ptrvalue;
      str2 = (CharPtr) vnp2->data.ptrvalue;
      if (str1 != NULL && str2 != NULL) {
        return StringICmp (str1, str2);
      } else {
        return 0;
      }
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}

static Boolean GatherOrgnamesFunc (GatherContextPtr gcp)

{
  BioSourcePtr     biop;
  OrgRefPtr        orp;
  ValNodePtr       sdp;
  SeqFeatPtr       sfp;
  CharPtr          str;
  ValNodePtr PNTR  vnpp;

  if (gcp == NULL || gcp->thisitem == NULL) return TRUE;
  if (gcp->thistype != OBJ_SEQFEAT  && gcp->thistype != OBJ_SEQDESC) return TRUE;
  vnpp = (ValNodePtr PNTR) gcp->userdata;
  if (vnpp == NULL) return TRUE;

  orp = NULL;
  biop = NULL;
  switch (gcp->thistype) {
    case OBJ_SEQFEAT :
      sfp = (SeqFeatPtr) gcp->thisitem;
      switch (sfp->data.choice) {
        case SEQFEAT_ORG :
          orp = (OrgRefPtr) sfp->data.value.ptrvalue;
          break;
        case SEQFEAT_BIOSRC :
          biop = (BioSourcePtr) sfp->data.value.ptrvalue;
          break;
        default :
          break;
      }
      break;
    case OBJ_SEQDESC :
      sdp = (ValNodePtr) gcp->thisitem;
      switch (sdp->choice) {
        case Seq_descr_org :
          orp = (OrgRefPtr) sdp->data.ptrvalue;
          break;
        case Seq_descr_source :
          biop = (BioSourcePtr) sdp->data.ptrvalue;
          break;
        default :
          break;
      }
      break;
    default :
      break;
  }

  if (orp == NULL && biop != NULL) {
    orp = biop->org;
  }
  if (orp != NULL) {
    str = StringSave (orp->taxname);
    if (str == NULL) return TRUE;
    TrimSpacesAroundString (str);
    if (StringICmp (str, "Human immunodeficiency virus type 1") == 0) {
      StringCpy (str, "HIV-1");
    } else if (StringICmp (str, "Human immunodeficiency virus type 2") == 0) {
      StringCpy (str, "HIV-2");
    }
    str [0] = TO_UPPER (str [0]);
    ValNodeAddStr (vnpp, 1, str);
  }
  return TRUE;
}

static void DefLineModFormAcceptProc (ButtoN b)

{
  EnumFieldAssocPtr  ap;
  Int2               count;
  DeflineFormPtr     dfp;
  GatherScope        gs;
  Boolean            labelMods;
  Boolean            leaveInParen;
  Int2               maxMods;
  Int2               mitochloroflag;
  ValNodePtr         nextvnp;
  ValNodePtr         nonUniqueOrgs;
  Int2               val;
  ValNodePtr         vnp;

  dfp = (DeflineFormPtr) GetObjectExtra (b);
  if (dfp == NULL) return;
  WatchCursor ();
  Hide (dfp->form);
  Update ();
  labelMods = GetStatus (dfp->addLabels);

  if (dfp->sourceBoxList != NULL && GetValue (dfp->customGrp) == 2) {
    count = 0;
    while (sourceModRankList [count] != NULL && dfp->sourceBoxList [count] != NULL) {
      if (! GetStatus (dfp->sourceBoxList [count])) {
        for (ap = orgmod_subtype_alist; ap->name != NULL; ap++) {
          if (StringICmp (ap->name, sourceModRankList [count]) == 0 &&
              ap->value > 0 && ap->value < 24) {
            orgmod_rank [ap->value] = 0;
          }
        }
        for (ap = subsource_subtype_alist; ap->name != NULL; ap++) {
          if (StringICmp (ap->name, sourceModRankList [count]) == 0 &&
              ap->value > 0 && ap->value < 24) {
            subsource_rank [ap->value] = 0;
          }
        }
      }
      count++;
    }
  }

  maxMods = INT2_MAX;
  val = GetValue (dfp->modLimit);
  if (val > 1) {
    maxMods = val - 1;
  }

  if (dfp->onlyModifyTargeted != NULL && GetStatus (dfp->onlyModifyTargeted)) {
    dfp->sep = GetBestTopParentForData (dfp->input_entityID, dfp->target);
  }

  leaveInParen = FALSE;
  if (dfp->leaveInParentheses != NULL && GetStatus (dfp->leaveInParentheses)) {
    leaveInParen = TRUE;
  }

  nonUniqueOrgs = NULL;
  if (dfp->smartMods) {
    MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
    gs.seglevels = 1;
    MemSet ((Pointer) (gs.ignore), (int)(TRUE), (size_t) (OBJ_MAX * sizeof(Boolean)));
    gs.ignore[OBJ_BIOSEQ] = FALSE;
    gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
    gs.ignore[OBJ_SEQANNOT] = FALSE;
    gs.ignore[OBJ_SEQFEAT] = FALSE;
    gs.ignore[OBJ_SEQDESC] = FALSE;
    GatherEntity (dfp->input_entityID, (Pointer) (&nonUniqueOrgs), GatherOrgnamesFunc, &gs);
    nonUniqueOrgs = SortValNode (nonUniqueOrgs, SortByName);
    vnp = nonUniqueOrgs;
    while (vnp != NULL) {
      nextvnp = vnp->next;
      if (nextvnp != NULL && StringICmp ((CharPtr) vnp->data.ptrvalue, (CharPtr) nextvnp->data.ptrvalue) == 0) {
        vnp->next = nextvnp->next;
        nextvnp->next = NULL;
        ValNodeFreeData (nextvnp);
        nextvnp = vnp;
        (vnp->choice)++;
      }
      vnp = nextvnp;
    }
  }

  mitochloroflag = GetValue (dfp->nucformitoorchloro) - 1;
  AutoDefProc (dfp->input_entityID, dfp->sep, TRUE, labelMods, maxMods, leaveInParen, NULL, NULL, nonUniqueOrgs, mitochloroflag, NULL);
  ValNodeFreeData (nonUniqueOrgs);
  ArrowCursor ();
  Remove (dfp->form);
  Update ();
  ObjMgrSetDirtyFlag (dfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, dfp->input_entityID, 0, 0);
}

static void ChangeCustomGrp (GrouP g)

{
  DeflineFormPtr  dfp;

  dfp = (DeflineFormPtr) GetObjectExtra (g);
  if (dfp == NULL) return;
  if (GetValue (g) == 1) {
    SafeDisable (dfp->sourceListGrp);
  } else {
    SafeEnable (dfp->sourceListGrp);
  }
}

static void CleanupDeflineForm (GraphiC g, VoidPtr data)

{
  DeflineFormPtr  dfp;

  dfp = (DeflineFormPtr) data;
  if (dfp != NULL) {
    MemFree (dfp->sourceBoxList);
  }
  StdCleanupFormProc (g, data);
}

static void DeflineMessageProc (ForM f, Int2 mssg)

{
  StdEditorProcsPtr  sepp;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    if (sepp->handleMessages != NULL) {
      sepp->handleMessages (f, mssg);
    }
  }
}

static ForM CreateDefLineModForm (Uint2 entityID, SeqEntryPtr sep, BioseqPtr target,
                                  Boolean smartMods, Int2 maxCount)

{
  ButtoN          b;
  GrouP           c;
  Int2            count;
  DeflineFormPtr  dfp;
  GrouP           h;
  Int2            i;
  GrouP           p;
  GrouP           q;
  Char            str [16];
  WindoW          w;

  w = NULL;
  dfp = MemNew (sizeof (DeflineForm));
  if (dfp != NULL) {
    w = FixedWindow (-50, -33, -10, -10, "Defline Modifier Customization", StdCloseWindowProc);
    SetObjectExtra (w, dfp, CleanupDeflineForm);
    dfp->form = (ForM) w;
    dfp->formmessage = DeflineMessageProc;
    dfp->input_entityID = entityID;
    dfp->sep = sep;
    dfp->target = target;
    dfp->smartMods = smartMods;

    h = HiddenGroup (w, -1, 0, NULL);
    SetGroupSpacing (h, 10, 10);

    dfp->addLabels = CheckBox (h, "Use Labels (e.g., 'strain BALB/c')", NULL);
    if (smartMods) {
      SetStatus (dfp->addLabels, TRUE);
    }

    p = NormalGroup (h, -1, 0, "Modifier Classes", programFont, ChangeCustomGrp);
    SetGroupSpacing (p, 10, 10);

    dfp->customGrp = HiddenGroup (p, 2, 0, ChangeCustomGrp);
    SetObjectExtra (dfp->customGrp, dfp, NULL);
    RadioButton (dfp->customGrp, "All");
    RadioButton (dfp->customGrp, "Custom");
    SetValue (dfp->customGrp, 1);

    dfp->sourceListGrp = HiddenGroup (p, 3, 0, NULL);
    count = 0;
    while (sourceModRankList [count] != NULL) {
      count++;
    }
    dfp->sourceBoxList = MemNew (sizeof (ButtoN) * (count + 3));
    count = 0;
    if (dfp->sourceBoxList != NULL) {
      while (sourceModRankList [count] != NULL && count < maxCount) {
        dfp->sourceBoxList [count] = CheckBox (dfp->sourceListGrp, sourceModRankList [count], NULL);
        count++;
      }
    }
    Disable (dfp->sourceListGrp);

    q = HiddenGroup (h, 4, 0, NULL);
    StaticPrompt (q, "Maximum modifiers per line", 0, popupMenuHeight, programFont, 'l');
    dfp->modLimit = PopupList (q, TRUE, NULL);
    PopupItem (dfp->modLimit, "no limit");
    for (i = 1; i <= count; i++) {
      sprintf (str, "%d", (int) i);
      PopupItem (dfp->modLimit, str);
    }
    if (smartMods) {
      SetValue (dfp->modLimit, 2);
    } else {
      SetValue (dfp->modLimit, 1);
    }

    dfp->nucformitoorchloro = HiddenGroup (h, -1, 0, NULL);
    RadioButton (dfp->nucformitoorchloro, "No mitochondrial or chloroplast suffix");
    RadioButton (dfp->nucformitoorchloro, "Nuclear gene(s) for mitochondrial product(s)");
    RadioButton (dfp->nucformitoorchloro, "Nuclear gene(s) for chloroplast product(s)");
    RadioButton (dfp->nucformitoorchloro, "Mitochondrial gene(s) for mitochondrial product(s)");
    RadioButton (dfp->nucformitoorchloro, "Chloroplast gene(s) for chloroplast product(s)");
    SetValue (dfp->nucformitoorchloro, 1);

    dfp->leaveInParentheses = CheckBox (w, "Leave in parenthetical organism info", NULL);
    SetStatus (dfp->leaveInParentheses, TRUE);

    dfp->onlyModifyTargeted = NULL;
    if (target != NULL) {
      dfp->onlyModifyTargeted = CheckBox (w, "Only modify targeted record", NULL);
    }

    c = HiddenGroup (w, 4, 0, NULL);
    b = PushButton (c, "Accept", DefLineModFormAcceptProc);
    SetObjectExtra (b, dfp, NULL);
    PushButton (c, "Cancel", StdCancelButtonProc);

    AlignObjects (ALIGN_CENTER, (HANDLE) dfp->addLabels, (HANDLE) dfp->customGrp,
                  (HANDLE) dfp->sourceListGrp, (HANDLE) q,
                  (HANDLE) dfp->leaveInParentheses,
                  (HANDLE) dfp->nucformitoorchloro, (HANDLE) c,
                  (HANDLE) dfp->onlyModifyTargeted, NULL);

    RealizeWindow (w);
  }
  if (smartMods) {
    DefLineModFormAcceptProc (b);
    return NULL;
  }
  return (ForM) w;
}

extern void GenerateAutomaticDefLinesCommon (IteM i, Boolean addMods, Boolean smartMods, ButtoN b)

{
  EnumFieldAssocPtr  ap;
  BaseFormPtr        bfp;
  Int2               count;
  ForM               f;
  Int2               maxCount;
  SeqEntryPtr        sep;
  BioseqPtr          target;

  if (b != NULL) {
    bfp = GetObjectExtra (b);
  } else {
#ifdef WIN_MAC
    bfp = currentFormDataPtr;
#else
    bfp = GetObjectExtra (i);
#endif
  }
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  MemSet ((Pointer) (orgmod_rank), (int)(0), sizeof(orgmod_rank));
  MemSet ((Pointer) (subsource_rank), (int)(0), sizeof(subsource_rank));

  count = 0;
  while (sourceModRankList [count] != NULL) {
    for (ap = orgmod_subtype_alist; ap->name != NULL; ap++) {
      if (StringICmp (ap->name, sourceModRankList [count]) == 0 &&
          ap->value > 0 && ap->value < 24) {
        orgmod_rank [ap->value] = count + 1;
      }
    }
    for (ap = subsource_subtype_alist; ap->name != NULL; ap++) {
      if (StringICmp (ap->name, sourceModRankList [count]) == 0 &&
          ap->value > 0 && ap->value < 24) {
        subsource_rank [ap->value] = count + 1;
      }
    }
    count++;
  }

  if (smartMods) {
    if (CountSeqEntryComponents (sep) == 1) {
      smartMods = FALSE;
    }
  }

  if (addMods || smartMods) {
    target =  GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID, bfp->input_itemtype);
    if (addMods) {
      maxCount = INT2_MAX;
    } else {
      maxCount = 6;
    }
    f = CreateDefLineModForm (bfp->input_entityID, sep, target, smartMods, maxCount);
    if (addMods && f != NULL) {
      Show (f);
      Select (f);
    } else if (bfp->activate != NULL) {
      bfp->activate ((WindoW) bfp->form);
    }
    return;
  }

  WatchCursor ();
  Update ();
  AutoDefProc (bfp->input_entityID, sep, FALSE, FALSE, INT2_MAX, TRUE, NULL, NULL, NULL, 0, NULL);
  ArrowCursor ();
  Update ();
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

void GenerateAutoDefLinesNoMods (IteM i)

{
  GenerateAutomaticDefLinesCommon (i, FALSE, FALSE, NULL);
}

void GenerateAutoDefLinesWithMods (IteM i)

{
  GenerateAutomaticDefLinesCommon (i, TRUE, FALSE, NULL);
}

void GenerateAutoDefLinesSmartMods (IteM i)

{
  GenerateAutomaticDefLinesCommon (i, FALSE, TRUE, NULL);
}

static void StringToLower (CharPtr str)

{
  Char  ch;

  if (str == NULL) return;
  ch = *str;
  while (ch != '\0') {
    *str = TO_LOWER (ch);
    str++;
    ch = *str;
  }
}

static void MakeNucleotideTitlesInSequinStyle (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  EnumFieldAssocPtr  ap;
  BioseqContextPtr   bcp;
  BioSourcePtr       biop;
  BioseqPtr          bsp;
  OrgModPtr          mod;
  OrgNamePtr         onp;
  OrgRefPtr          orp;
  SeqDescrPtr        sdp;
  SubSourcePtr       ssp;
  CharPtr            str;
  Char               text [256];
  ValNodePtr         ttl;
  ValNodePtr         vnp;

  if (sep == NULL) return;
  if (! IS_Bioseq (sep)) return;
  bsp = sep->data.ptrvalue;
  if (bsp == NULL) return;
  if (! ISA_na (bsp->mol)) return;
  bcp = BioseqContextNew (bsp);
  sdp = BioseqContextGetSeqDescr (bcp, Seq_descr_source, NULL, NULL);
  BioseqContextFree (bcp);
  if (sdp == NULL) return;
  biop = (BioSourcePtr) sdp->data.ptrvalue;
  if (biop == NULL) return;
  if (bsp->descr != NULL) {
    vnp = ValNodeExtract (&(bsp->descr), Seq_descr_title);
    vnp = ValNodeFreeData (vnp);
  }
  str = MemNew (2000);

  orp = biop->org;
  if (orp != NULL) {
    StringCpy (text, "[org=");
    StringCat (text, orp->taxname);
    StringCat (text, "] ");
    StringCat (str, text);
  }

  ssp = biop->subtype;
  while (ssp != NULL) {
    for (ap = subsource_subtype_alist; ap->name != NULL; ap++) {
      if (ssp->subtype == ap->value) {
        StringCpy (text, "[");
        if (ap->value != 255) {
          StringCat (text, ap->name);
        } else {
          StringCat (text, "subsource");
        }
        StringToLower (text);
        StringCat (text, "=");
        StringCat (text, ssp->name);
        StringCat (text, "] ");
        StringCat (str, text);
      }
    }
    ssp = ssp->next;
  }
  orp = biop->org;
  if (orp != NULL) {
    onp = orp->orgname;
    if (onp != NULL) {
      mod = onp->mod;
      while (mod != NULL) {
        for (ap = orgmod_subtype_alist; ap->name != NULL; ap++) {
          if (mod->subtype == ap->value) {
            StringCpy (text, "[");
            if (ap->value != 255) {
              StringCat (text, ap->name);
            } else {
              StringCat (text, "note");
            }
            StringToLower (text);
            StringCat (text, "=");
            StringCat (text, mod->subname);
            StringCat (text, "] ");
            StringCat (str, text);
          }
        }
        mod = mod->next;
      }
    }
  }

  TrimSpacesAroundString (str);
  if (! StringHasNoText (str)) {
    ttl = CreateNewDescriptor (sep, Seq_descr_title);
    if (ttl != NULL) {
      ttl->data.ptrvalue = StringSave (str);
    }
  }
  MemFree (str);
}

extern Int2 LIBCALLBACK MakeSequinProteinTitles (Pointer data);
extern Int2 LIBCALLBACK MakeSequinNucleotideTitles (Pointer data);
extern Int2 LIBCALLBACK MakeSequinFeatureTable (Pointer data);
extern Int2 LIBCALLBACK MakeContigBuildTable (Pointer data);

extern Int2 LIBCALLBACK MakeSequinNucleotideTitles (Pointer data)

{
  OMProcControlPtr  ompcp;
  SeqEntryPtr       sep;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  switch (ompcp->input_itemtype) {
    case OBJ_BIOSEQ :
      break;
    case OBJ_BIOSEQSET :
      break;
    case 0 :
      return OM_MSG_RET_ERROR;
    default :
      return OM_MSG_RET_ERROR;
  }
  if (ompcp->input_data == NULL) return OM_MSG_RET_ERROR;
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  if (sep == NULL) return OM_MSG_RET_ERROR;
  SeqEntryExplore (sep, NULL, MakeNucleotideTitlesInSequinStyle);
  ObjMgrSetDirtyFlag (ompcp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, ompcp->input_entityID, 0, 0);
  return OM_MSG_RET_DONE;
}

static void MakeProteinTitlesInSequinStyle (Uint2 entityID, SeqEntryPtr sep)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  DefFeatsPtr   dfp;
  GeneRefPtr    grp;
  GatherScope   gs;
  ValNodePtr    head;
  SeqEntryPtr   nsep;
  ProtRefPtr    prp;
  SeqEntryPtr   psep;
  Char          str [256];
  Char          text [256];
  ValNodePtr    ttl;
  ValNodePtr    vnp;

  if (sep == NULL) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;
    if (bssp->_class == 7 || bssp->_class == 13 ||
        bssp->_class == 14 || bssp->_class == 15) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        MakeProteinTitlesInSequinStyle (entityID, sep);
      }
      return;
    }
  }
  nsep = FindNucSeqEntry (sep);
  if (nsep == NULL) return;

  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = TRUE;
  MemSet ((Pointer) (gs.ignore), (int)(TRUE), (size_t) (OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.scope = sep;
  gs.target = NULL;
  head = NULL;
  GatherEntity (entityID, (Pointer) (&head), GetCDStRNArRNAGatherFunc, &gs);
  /* head = SortValNode (head, SortCDStRNArRNAByLocation); */
  if (head == NULL) return;

  vnp = head;
  while (vnp != NULL) {
    dfp = (DefFeatsPtr) vnp->data.ptrvalue;
    if (dfp != NULL && dfp->sfp != NULL && dfp->subtype == FEATDEF_CDS) {
      FindGeneAndProtForCDS (entityID, dfp->sfp, &(dfp->gene), &(dfp->prot));
      bsp = GetBioseqGivenSeqLoc (dfp->sfp->product, entityID);
      if (bsp != NULL) {
        str [0] = '\0';
        if (dfp->gene != NULL) {
          grp = (GeneRefPtr) dfp->gene->data.value.ptrvalue;
          if (grp != NULL) {
            StringNCpy_0 (text, (CharPtr) grp->locus, sizeof (text));
            if (! StringHasNoText (text)) {
              StringCat (str, "[gene=");
              StringCat (str, text);
              StringCat (str, "]");
            }
          }
        }
        if (dfp->prot != NULL) {
          prp = (ProtRefPtr) dfp->prot->data.value.ptrvalue;
          if (prp != NULL && prp->name != NULL) {
            StringNCpy_0 (text, (CharPtr) prp->name->data.ptrvalue, sizeof (text));
            if (! StringHasNoText (text)) {
              if (str [0] != '\0') {
                StringCat (str, " ");
              }
              StringCat (str, "[prot=");
              StringCat (str, text);
              StringCat (str, "]");
            }
          }
        }
        if (! StringHasNoText (str)) {
          psep = SeqMgrGetSeqEntryForData (bsp);
          if (psep != NULL) {
            ttl = CreateNewDescriptor (psep, Seq_descr_title);
            if (ttl != NULL) {
              ttl->data.ptrvalue = StringSave (str);
            }
          }
        }
      }
    }
    vnp = vnp->next;
  }
  ValNodeFreeData (head);
}

extern Int2 LIBCALLBACK MakeSequinProteinTitles (Pointer data)

{
  OMProcControlPtr  ompcp;
  SeqEntryPtr       sep;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  switch (ompcp->input_itemtype) {
    case OBJ_BIOSEQ :
      break;
    case OBJ_BIOSEQSET :
      break;
    case 0 :
      return OM_MSG_RET_ERROR;
    default :
      return OM_MSG_RET_ERROR;
  }
  if (ompcp->input_data == NULL) return OM_MSG_RET_ERROR;
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  if (sep == NULL) return OM_MSG_RET_ERROR;
  MakeProteinTitlesInSequinStyle (ompcp->input_entityID, sep);
  ObjMgrSetDirtyFlag (ompcp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, ompcp->input_entityID, 0, 0);
  return OM_MSG_RET_DONE;
}

static void PrintFtableIntervals (FILE *fp, SeqLocPtr location, Uint2 entityID, CharPtr label)

{
  BioseqPtr  bsp;
  SeqLocPtr  slp;
  Int4       start;
  Int4       stop;

  if (fp == NULL || location == NULL || StringHasNoText (label)) return;
  bsp = BioseqFindFromSeqLoc (location);
  if (bsp == NULL) return;
  slp = SeqLocFindNext (location, NULL);
  if (slp == NULL) return;
  start = GetOffsetInBioseq (slp, bsp, SEQLOC_START) + 1;
  stop = GetOffsetInBioseq (slp, bsp, SEQLOC_STOP) + 1;
  fprintf (fp, "%ld\t%ld\t%s\n", (long) start, (long) stop, label);
  while ((slp = SeqLocFindNext (location, slp)) != NULL) {
    start = GetOffsetInBioseq (slp, bsp, SEQLOC_START) + 1;
    stop = GetOffsetInBioseq (slp, bsp, SEQLOC_STOP) + 1;
    if (start != 0 && stop != 0) {
      fprintf (fp, "%ld\t%ld\n", (long) start, (long) stop);
    }
  }
}

static void MakeFeatureInSequinStyle (Uint2 entityID, SeqEntryPtr sep, FILE *fp)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  DefFeatsPtr   dfp;
  GeneRefPtr    grp;
  GatherScope   gs;
  ValNodePtr    head;
  SeqEntryPtr   nsep;
  ProtRefPtr    prp;
  RnaRefPtr     rrp;
  SeqIdPtr      sip;
  Char          str [256];
  ValNodePtr    syn;
  ValNodePtr    vnp;

  if (sep == NULL || fp == NULL) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;
    if (bssp->_class == 7 || bssp->_class == 13 ||
        bssp->_class == 14 || bssp->_class == 15) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        MakeFeatureInSequinStyle (entityID, sep, fp);
      }
      return;
    }
  }
  nsep = FindNucSeqEntry (sep);
  if (nsep == NULL) return;
  bsp = nsep->data.ptrvalue;
  if (bsp == NULL) return;

  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  gs.get_feats_location = TRUE;
  MemSet ((Pointer) (gs.ignore), (int)(TRUE), (size_t) (OBJ_MAX * sizeof(Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  gs.ignore[OBJ_SEQFEAT] = FALSE;
  gs.scope = sep;
  gs.target = NULL;
  head = NULL;
  GatherEntity (entityID, (Pointer) (&head), GetGeneCDStRNArRNAmRNAGatherFunc, &gs);
  head = SortValNode (head, SortCDStRNArRNAByLocation);
  if (head == NULL) return;

  sip = SeqIdFindBest (bsp->id, 0);
  SeqIdWrite (sip, str, PRINTID_FASTA_LONG, sizeof (str) - 1);
  if (! StringHasNoText (str)) {
    fprintf (fp, ">Feature %s\n", str);
  }

  vnp = head;
  while (vnp != NULL) {
    dfp = (DefFeatsPtr) vnp->data.ptrvalue;
    if (dfp != NULL && dfp->sfp != NULL) {
      str [0] = '\0';
      switch (dfp->subtype) {
        case FEATDEF_GENE :
          PrintFtableIntervals (fp, dfp->sfp->location, entityID, "gene");
          grp = (GeneRefPtr) dfp->sfp->data.value.ptrvalue;
          if (grp != NULL) {
            StringNCpy_0 (str, (CharPtr) grp->locus, sizeof (str));
            if (! StringHasNoText (str)) {
              fprintf (fp, "\t\t\tgene\t%s\n", str);
            }
            for (syn = grp->syn; syn != NULL; syn = syn->next) {
              StringNCpy_0 (str, (CharPtr) syn->data.ptrvalue, sizeof (str));
              if (! StringHasNoText (str)) {
                fprintf (fp, "\t\t\tgene\t%s\n", str);
              }
            }
          }
          break;
        case FEATDEF_CDS :
          PrintFtableIntervals (fp, dfp->sfp->location, entityID, "CDS");
          FindGeneAndProtForCDS (entityID, dfp->sfp, &(dfp->gene), &(dfp->prot));
          if (dfp->gene != NULL) {
            grp = (GeneRefPtr) dfp->gene->data.value.ptrvalue;
            if (grp != NULL) {
              StringNCpy_0 (str, (CharPtr) grp->locus, sizeof (str));
              if (! StringHasNoText (str)) {
                fprintf (fp, "\t\t\tgene\t%s\n", str);
              }
            }
          }
          if (dfp->prot != NULL) {
            prp = (ProtRefPtr) dfp->prot->data.value.ptrvalue;
            if (prp != NULL) {
              if (prp->name != NULL) {
                StringNCpy_0 (str, (CharPtr) prp->name->data.ptrvalue, sizeof (str));
              } else if (prp->desc != NULL) {
                StringNCpy_0 (str, prp->desc, sizeof (str));
              }
              if (! StringHasNoText (str)) {
                fprintf (fp, "\t\t\tproduct\t%s\n", str);
              }
            }
          }
          break;
        case FEATDEF_mRNA :
          PrintFtableIntervals (fp, dfp->sfp->location, entityID, "mRNA");
          FindGeneAndProtForCDS (entityID, dfp->sfp, &(dfp->gene), NULL);
          rrp = (RnaRefPtr) dfp->sfp->data.value.ptrvalue;
          if (rrp != NULL) {
            if (rrp->ext.choice == 1) {
              StringNCpy_0 (str, (CharPtr) rrp->ext.value.ptrvalue, sizeof (str));
            }
            if (! StringHasNoText (str)) {
              fprintf (fp, "\t\t\tproduct\t%s\n", str);
            }
          }
          break;
      }
    }
    vnp = vnp->next;
  }
  ValNodeFreeData (head);
}

static Boolean LIBCALLBACK SequinFTableFeature (SeqFeatPtr sfp, SeqMgrFeatContextPtr context)

{
  DbtagPtr     dbt;
  FILE         *fp;
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
  tRNAPtr      trna;
  ValNodePtr   vnp;

  if (sfp == NULL) return TRUE;
  fp = (FILE *) context->userdata;
  label = (CharPtr) FeatDefTypeLabel (sfp);
  if (StringCmp (label, "Gene") == 0) {
    label = "gene";
  }
  PrintFtableIntervals (fp, sfp->location, context->entityID, label);
  switch (context->seqfeattype) {
    case SEQFEAT_GENE :
      grp = (GeneRefPtr) sfp->data.value.ptrvalue;
      if (grp != NULL) {
        StringNCpy_0 (str, (CharPtr) grp->locus, sizeof (str));
        if (! StringHasNoText (str)) {
          fprintf (fp, "\t\t\tgene\t%s\n", str);
        }
        for (vnp = grp->syn; vnp != NULL; vnp = vnp->next) {
          StringNCpy_0 (str, (CharPtr) vnp->data.ptrvalue, sizeof (str));
          if (! StringHasNoText (str)) {
            fprintf (fp, "\t\t\tgene_syn\t%s\n", str);
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
                fprintf (fp, "\t\t\tproduct\t%s\n", str);
              }
            }
          } else if (prp->desc != NULL) {
            StringNCpy_0 (str, prp->desc, sizeof (str));
            if (! StringHasNoText (str)) {
              fprintf (fp, "\t\t\tproduct\t%s\n", str);
            }
          }
          for (vnp = prp->activity; vnp != NULL; vnp = vnp->next) {
            StringNCpy_0 (str, (CharPtr) vnp->data.ptrvalue, sizeof (str));
            if (! StringHasNoText (str)) {
              fprintf (fp, "\t\t\tfunction\t%s\n", str);
            }
          }
          for (vnp = prp->ec; vnp != NULL; vnp = vnp->next) {
            StringNCpy_0 (str, (CharPtr) vnp->data.ptrvalue, sizeof (str));
            if (! StringHasNoText (str)) {
              fprintf (fp, "\t\t\tEC_number\t%s\n", str);
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
              fprintf (fp, "\t\t\tprotein_id\t%s\n", str);
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
              fprintf (fp, "\t\t\tproduct\t%s\n", str);
            }
            break;
          case 2 :
            trna = rrp->ext.value.ptrvalue;
            if (trna != NULL) {
              FeatDefLabel (sfp, str, sizeof (str) - 1, OM_LABEL_CONTENT);
              if (! StringHasNoText (str)) {
                fprintf (fp, "\t\t\tproduct\t%s\n", str);
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
    fprintf (fp, "\t\t\tnote\t%s\n", sfp->comment);
  }
  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    if (! StringHasNoText (gbq->qual)) {
      if (! StringHasNoText (gbq->val)) {
        fprintf (fp, "\t\t\t%s\t%s\n", gbq->qual, gbq->val);
      }
    }
  }
  for (vnp = sfp->dbxref; vnp != NULL; vnp = vnp->next) {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt != NULL) {
      if (! StringHasNoText (dbt->db)) {
        oip = dbt->tag;
        if (oip->str != NULL && (! StringHasNoText (oip->str))) {
          fprintf (fp, "\t\t\tdb_xref\t%s:%s\n", dbt->db, oip->str);
        } else {
          fprintf (fp, "\t\t\tdb_xref\t%s:%ld\n", dbt->db, (long) oip->id);
        }
      }
    }
  }
  return TRUE;
}

static Boolean LIBCALLBACK SequinFTableBioseq (BioseqPtr bsp, SeqMgrBioseqContextPtr context)

{
  FILE      *fp;
  SeqIdPtr  sip;
  Char      str [42];

  if (bsp == NULL) return TRUE;
  fp = (FILE *) context->userdata;
  sip = SeqIdFindBest (bsp->id, 0);
  SeqIdWrite (sip, str, PRINTID_FASTA_LONG, sizeof (str) - 1);
  if (! StringHasNoText (str)) {
    fprintf (fp, ">Feature %s\n", str);
  }
  SeqMgrExploreFeatures (bsp, (Pointer) fp, SequinFTableFeature, NULL, NULL, NULL);
  return TRUE;
}

extern Int2 LIBCALLBACK MakeSequinFeatureTable (Pointer data)

{
  FILE              *fp;
  OMProcControlPtr  ompcp;
  Char              path [PATH_MAX];
  SeqEntryPtr       sep;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  switch (ompcp->input_itemtype) {
    case OBJ_BIOSEQ :
      break;
    case OBJ_BIOSEQSET :
      break;
    case 0 :
      return OM_MSG_RET_ERROR;
    default :
      return OM_MSG_RET_ERROR;
  }
  if (ompcp->input_data == NULL) return OM_MSG_RET_ERROR;
  sep = GetTopSeqEntryForEntityID (ompcp->input_entityID);
  if (sep == NULL) return OM_MSG_RET_ERROR;
  TmpNam (path);
  fp = FileOpen (path, "w");
  if (fp == NULL) return OM_MSG_RET_ERROR;
  if (SeqMgrFeaturesAreIndexed (ompcp->input_entityID) != 0) {
    SeqMgrExploreBioseqs (ompcp->input_entityID, NULL, (Pointer) fp, SequinFTableBioseq, TRUE, FALSE, FALSE);
  } else {
    MakeFeatureInSequinStyle (ompcp->input_entityID, sep, fp);
  }
  FileClose (fp);
  LaunchGeneralTextViewer (path, "Gene - CDS Feature Table");
  FileRemove (path);
  return OM_MSG_RET_DONE;
}

/* the following two functions are modified from PrintGenome in wprint.c */

static void SqLitPrintGenome(SeqLitPtr slp, FILE *fp)
{
	static Char		val[166];

	if (slp->seq_data != NULL)         /* not a gap */
	{
		if (slp->length == 0)  /* unknown length */
		{
			fprintf(fp, "gap\t0\t0\t0\n");
		} else {
/* don't know what to do here */
		}
	} else {                  /* gap length was set */
			fprintf(fp, "gap\t0\t0\t%ld\n", (long) slp->length);
	}
}

static void SqLocPrintGenome(SeqLocPtr slp_head, FILE *fp)
{
	SeqLocPtr	slp;
	static Char		buf[11];
	SeqIdPtr	sid, newid;
	Int4 		start, stop;
		
	for (slp = slp_head; slp; slp = slp->next) {
		sid = SeqLocId(slp);
		if (slp->choice == SEQLOC_INT || slp->choice == SEQLOC_WHOLE) {
			start = SeqLocStart(slp);
			stop = SeqLocStop(slp);
		} else if (slp->choice == SEQLOC_NULL){
			fprintf(fp, "gap\t0\t0\t0\n");
			continue;
		} else {
			continue;
		}
		if (sid->choice == SEQID_GI) {
			newid = GetSeqIdForGI(sid->data.intvalue);
		} else if (sid->choice == SEQID_GENERAL) {
			newid = sid;
		} else {
			newid = sid;
		}
		SeqIdWrite(newid, buf, PRINTID_TEXTID_ACCESSION, 10);
		fprintf(fp, "%s\t%ld\t%ld\n", buf, (long) start+1, (long) stop+1);
	}
}

static Boolean DeltaLitOnly (BioseqPtr bsp)

{
  ValNodePtr  vnp;

  if (bsp == NULL || bsp->repr != Seq_repr_delta) return FALSE;
  for (vnp = (ValNodePtr)(bsp->seq_ext); vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == 1) return FALSE;
  }
  return TRUE;
}

extern Int2 LIBCALLBACK MakeContigBuildTable (Pointer data)

{
  BioseqPtr         bsp;
  DeltaSeqPtr       dsp;
  FILE              *fp;
  SeqLitPtr 	    litp;
  OMProcControlPtr  ompcp;
  Char              path [PATH_MAX];
  SeqEntryPtr       sep;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  bsp = NULL;
  switch (ompcp->input_itemtype) {
    case OBJ_BIOSEQ :
      bsp = (BioseqPtr) ompcp->input_data;
      break;
    default :
      return OM_MSG_RET_ERROR;
  }
  if (bsp == NULL) return OM_MSG_RET_ERROR;
  sep = GetBestTopParentForData (ompcp->input_entityID, bsp);
  if (sep == NULL) return OM_MSG_RET_ERROR;
  if (IsAGenomeRecord (sep) ||
      IsSegmentedBioseqWithoutParts (sep)) {
  } else if (IsADeltaBioseq (sep) && (! DeltaLitOnly (bsp))) {
  } else return OM_MSG_RET_ERROR;
  TmpNam (path);
  fp = FileOpen (path, "w");
  if (fp == NULL) return OM_MSG_RET_ERROR;
/* the following code is modified from PrintGenome in wprint.c */
	if (bsp->seq_ext_type == 1) {
		SqLocPrintGenome((SeqLocPtr) bsp->seq_ext, fp);
	} else if (bsp->seq_ext_type == 4) {
		for (dsp = (DeltaSeqPtr) bsp->seq_ext; dsp; dsp=dsp->next) {
			if (dsp->choice == 1) {  /* SeqLoc */
				SqLocPrintGenome((SeqLocPtr)(dsp->data.ptrvalue), fp);
			} else {
				litp = (SeqLitPtr)(dsp->data.ptrvalue);
				if (litp == NULL) continue;
				SqLitPrintGenome(litp, fp);
			}
		}
	}
  FileClose (fp);
  LaunchGeneralTextViewer (path, "Contig Build Table");
  FileRemove (path);
  return OM_MSG_RET_DONE;
}

/* resolve existing colliding ids section */

typedef struct lclidlist {
  BioseqPtr  firstbsp;
  SeqIdPtr   firstsip;
  CharPtr    key;
  Int2       count;
  struct lclidlist PNTR left;
  struct lclidlist PNTR right;
} LclIdList, PNTR LclIdListPtr;

static void ReplaceLocalID (BioseqPtr bsp, SeqIdPtr sip, CharPtr key, Int2 count)

{
  ObjectIdPtr  oip;
  Char         str [64];
  Char         tmp [70];

  if (bsp == NULL || sip == NULL || StringHasNoText (key)) return;
  oip = (ObjectIdPtr) sip->data.ptrvalue;
  if (oip == NULL) return;
  StringNCpy_0 (str, key, sizeof (str));
  sprintf (tmp, "%s__%d", str, (int) count);
  oip->str = MemFree (oip->str);
  oip->str = StringSave (tmp);
  SeqMgrReplaceInBioseqIndex (bsp);
}

static void BuildLclTree (LclIdListPtr PNTR head, BioseqPtr bsp, CharPtr x, SeqIdPtr sip)

{
  Int2          comp;
  LclIdListPtr  idlist;

  if (*head != NULL) {
    idlist = *head;
    comp = StringICmp (idlist->key, x);
    if (comp < 0) {
      BuildLclTree (&(idlist->right), bsp, x, sip);
    } else if (comp > 0) {
      BuildLclTree (&(idlist->left), bsp, x, sip);
    } else {
      if (idlist->firstbsp != NULL && idlist->firstsip != NULL) {
        ReplaceLocalID (idlist->firstbsp, idlist->firstsip, x, 1);
        idlist->count = 2;
        idlist->firstbsp = NULL;
        idlist->firstsip = NULL;
      }
      ReplaceLocalID (bsp, sip, x, idlist->count);
      (idlist->count)++;
    }
  } else {
    idlist = MemNew (sizeof (LclIdList));
    if (idlist != NULL) {
      *head = idlist;
      idlist->firstbsp = bsp;
      idlist->firstsip = sip;
      idlist->count = 1;
      idlist->key = StringSave (x);
      idlist->left = NULL;
      idlist->right = NULL;
    }
  }
}

static void FreeLclTree (LclIdListPtr PNTR head)

{
  LclIdListPtr  idlist;

  if (head != NULL && *head != NULL) {
    idlist = *head;
    FreeLclTree (&(idlist->left));
    FreeLclTree (&(idlist->right));
    MemFree (idlist->key);
    MemFree (idlist);
  }
}

static void ResolveExistingIDsCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr          bsp;
  LclIdListPtr PNTR  head;
  SeqIdPtr           sip;
  Char               str [64];

  head = (LclIdListPtr PNTR) mydata;
  if (sep == NULL || head == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp != NULL) {
      for (sip = bsp->id; sip != NULL; sip = sip->next) {
        if (sip->choice == SEQID_LOCAL) {
          SeqIdWrite (sip, str, PRINTID_REPORT, sizeof (str));
          BuildLclTree (head, bsp, str, sip);
        }
      }
    }
  }
}

extern Int2 DoOneSegFixup (SeqEntryPtr sep, Boolean ask);

extern void ResolveExistingLocalIDs (IteM i);
extern void ResolveExistingLocalIDs (IteM i)

{
  MsgAnswer     ans;
  BaseFormPtr   bfp;
  Boolean       doParts = FALSE;
  LclIdListPtr  head = NULL;
  SeqEntryPtr   sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  ans = Message (MSG_YN, "Do you wish to also reconstruct segmented bioseqs?");
  doParts = (Boolean) (ans == ANS_YES);
  SeqEntryExplore (sep, (Pointer) &head, ResolveExistingIDsCallback);
  FreeLclTree (&head);
  if (doParts) {
    DoOneSegFixup (sep, FALSE);
  }
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

