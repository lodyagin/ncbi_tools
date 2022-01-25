/*   asn2graphic.c
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
* File Name:  asn2graphic.c
* Author:  Eric Northup
*
* Version Creation Date:   11/8/01
*
* $Revision: 6.70 $
*
* File Description:
*
* Modifications:
* --------------------------------------------------------------------------
*
* ==========================================================================
*/

#include <asn2graphicp.h>
#include <sequtil.h>
#include <sqnutils.h>
#include <explore.h>
#include <drawingp.h>
#include <alignmgr2.h>

#define RELEVANT_FEATS_PER_CHUNK 128

/* these do not apply to individual features; they are used, ie in "AddFeatureToFilter (FEATDEF_ANY_RNA. . .)
   that would make the filter include _all_ types of RNA */
#define FEATDEF_ANY_RNA (FEATDEF_MAX + 2)
#define FEATDEF_ANY_PROT (FEATDEF_MAX + 1)

typedef struct renderInput {
  AppearanceItemPtr AIP;
  FilterItemPtr     FIP;
  RelevantFeatureItemPtr RFIP;
  SegmenT         labelSeg;
  SegmenT         drawSeg;
  Int4            yStart;
  Int4            decorationLeft;
  Int4            decorationRight;
  Int4            decorationHeight;
  Uint2           featureOffset;
  Uint2           Height, rowHeight;
} RenderInput, PNTR RenderInputPtr;

typedef struct viewerInstanceData {
  BioseqPtr       BSP;
  Uint4           viewScale;
  Uint4           seqLength;
  Int4Ptr         ceiling;
  SeqAnnotPtr PNTR sapList;
  Uint2           sapCount;
  SegmenT         drawOnMe;
  SegmenT         topLevelSeg;
  Uint2           featureCount;
  ValNodePtr      featVNP;      /* data.ptrvalue == RelevantFeatureItem [RELEVANT_FEATS_PER_CHUNK] */
  ValNodePtr      BSPsegmentVNP;
  AppearancePtr   AppPtr;
  FilterPtr       FltPtr;
  LayoutAlgorithm overrideLayout;
  Int4            gphheight;
  SegmenT         gphseg;
  Int4            gphyOffset;
  Int4            from;
  Int4            to;
  Boolean         allFeatures;
} ViewerContext, PNTR ViewerContextPtr;

typedef struct featureFilterState {
  ValNodePtr       currentRFIPblockVNP;
  Uint4            featureBlockOffset;
  Uint2            indexInBlock;
} FeatureFilterState, PNTR FeatureFilterStatePtr;

typedef struct AlignmentFilterState {
  SeqAlignPtr      SAPhead, SAPcurrent;
  Uint2            align;  
} AlignmentFilterState, PNTR AlignmentFilterStatePtr;

/*
Bioseq's are still special-case for now.  todo: make bioseq segments use this mechanism
*/

typedef union filterState {
  FeatureFilterState feat;
  AlignmentFilterState align;
} FilterState, PNTR FilterStatePtr;  

typedef struct filterProcessState {
  FilterState      state;
  FilterItemPtr    currentFIP;
  ValNodePtr       currentFilterVNP;
  ValNodePtr       needFreeList; /* things to be freed with ValNodeFreeData() after finishing this filter */
  SegmenT          labelSegs [FEATDEF_MAX];
  SegmenT          drawSegs [FEATDEF_MAX];
  Int4             ceiling;
  RenderInput      renderParm;
  Uint4            featuresProcessedCount;
  BoolPtr          featuresProcessed;    /*points to an array of boolean [vContext->featureCount] */
  ViewerContextPtr vContext;
} FilterProcessState, PNTR FilterProcessStatePtr;

typedef void (*RenderFuncPtr)    (RenderInputPtr rip, ViewerContextPtr vContext);
typedef void (*GetDimensionsPtr) (RenderInputPtr, Int4Ptr start, Int4Ptr stop, Uint2Ptr height, ViewerContextPtr vContext);

typedef struct renderClass {
  RenderFuncPtr RenderFunc;
  GetDimensionsPtr GetDimensions;
} RenderClass, PNTR RenderClassPtr;

typedef struct internalRow {
  Uint2           rowFeatureCount;
  DataVal         layoutData;
  ValNodePtr      feats;        /* data.ptrvalue == RelevantFeatureItemPtr */
  struct internalRow PNTR next;
} InternalRow, PNTR InternalRowPtr;
  ValNodePtr      cleanupList;  /* */
/* returns the total number of _rows_ */
typedef Uint2   (PNTR LayoutFunction) (InternalRowPtr firstRow, FilterProcessStatePtr FPSP);


static CharPtr LayoutStrings [] = {
  "Inherit", /* Do not override layout in filter */
  "Diagonal", /* One feature per row */
  /* "DiagonalSawtooth" */
  "ByType", /* like "compact", except that no row will contain _different_ types of features (this may use more rows) */
  "Smear", /* one row contains all features of a given type (overlapping features will render as an anonymous smear */
  /*  "Single Row", -- not working currently*/ /* _all_ features are rendered in a single row (equivalent to "smear" if given only one type of features)*/
  "Compact", /* Overlapping features are drawn in different rows, but this trys to minimize the number of rows, */
  "GeneProducts", /* not yet implemented. . . */
  "GeneProductsX", /* not yet implemented. . . */
  NULL
};

static LayoutAlgorithm LayoutValues [] = {
  Layout_Inherit,
  Layout_Diagonal,
  /* Layout_DiagonalSawtooth, */
  Layout_FeatTypePerLineGroup,
  Layout_FeatTypePerLine,
  /*  Layout_AllInOneLine, not working currently (and not all that useful, either)*/
  Layout_PackUpward,
  Layout_GroupCorrespondingFeats,
  Layout_GroupCorrespondingFeatsRepeat
};

static CharPtr  LlocStrings [] = {
  "none",
  "inside",
  "above",
  "below",
  "left",
  "right",
  NULL
};

static LabelLocEnum LlocValues [] = {
  LabelNone,
  LabelInside,
  LabelAbove,
  LabelBelow,
  LabelLeft,
  LabelRight
};

static CharPtr  BoolStrings [] = {
  "yes",
  "true",
  "on",
  "no",
  "false",
  "off",
  NULL
};

static Boolean  BoolValues [] = {
  TRUE,
  TRUE,
  TRUE,
  FALSE,
  FALSE,
  FALSE
};

static CharPtr  RenderStrings [] = {
  "none",
  "line",
  "cappedline",
  "box",
  "outlinebox",
  NULL
};

static RenderAlgorithm RenderValues [] = {
  Do_Not_Render,
  Render_Line,
  Render_CappedLine,
  Render_Box,
  Render_OutlineBox
};

static CharPtr  GapStrings [] = {
  "none",
  "line",
  "angle",
  NULL
};

static GapEnum  GapValues [] = {
  NoGap,
  LineGap,
  AngleGap
};

static CharPtr GroupLabelLocations [] = {
  "above",
  "top",
  "left",
  "below",
  "none",
  NULL
};

static GroupLabelLocation GroupLabelLocationValues [] = {
  LabelOnTop,
  LabelOnTop,
  LabelOnSide,
  LabelOnBottom,
  NoLabel
};

static CharPtr BioseqFormat [] = {
  "accession",
  "accn",
  "fasta",
  "long",
  NULL
};

static Uint1 BioseqFormatValues [] = {
  PRINTID_TEXTID_ACC_VER,
  PRINTID_TEXTID_ACC_VER,
  PRINTID_FASTA_SHORT,
  PRINTID_FASTA_LONG
};

static CharPtr StrandStrings [] = {
  "both",
  "minus",
  "plus",
/* these could be added if people would find them useful. . .
   "-",
   "+",
*/
  NULL
};

static StrandChoice StrandValues [] = {
  BothStrands,
  MinusStrand,
  PlusStrand
};

static CharPtr ColorStrings [] = {
"black",    /*0,0,0*/
"blue",     /*0,0,255*/
"brown",    /*133,62,38*/
"coral",    /*255,127,80*/
"cyan",     /*0,255,255*/
"gray",     /*127,127,127*/
"green",    /*0,255,0*/
"grey",     /*127,127,127*/
"lavender", /*230,230,250*/
"magenta",  /*255,0,255*/
"maroon",   /*176,48,96*/
"orange",   /*255,165,0*/
"olive",    /*107,142,35*/
"pink",     /*255,192,203*/
"purple",   /*175,0,255*/
"red",      /*255,0,0*/
"white",    /*255,255,255*/
"yellow",   /*255,255,0*/
NULL
};

static Uint1 ColorValues[][3] = {
/*black*/    {0,   0,   0},
/*blue*/     {0,   0,   255},
/*brown*/    {133, 62,  38},
/*coral*/    {255, 127, 80},
/*cyan*/     {0,   255, 255},
/*gray*/     {127, 127, 127},
/*green*/    {0,   255, 0},
/*grey*/     {127, 127, 127},
/*lavender*/ {230, 230, 250},
/*magenta*/  {255, 0,   255},
/*maroon*/   {176, 48,  96},
/*orange*/   {255, 165, 0},
/*olive*/    {107, 142, 35},
/*pink*/     {255, 192, 203},
/*purple*/   {175, 0,   255},
/*red*/      {255, 0,   0},
/*white*/    {255, 255, 255},
/*yellow*/   {255, 255, 0}
};

static Char  config_filename [] = "asn2gph";

static Int1 StringIndexInStringList (CharPtr testString, CharPtr PNTR stringList) {

  Int1  i;

  if (testString == NULL || stringList == NULL) return -1;
  i = 0;
  while (stringList [i] != NULL) {
    if (StringICmp (testString, stringList [i]) == 0) return i;
    i++;
  }
  return -1;
}

/* this parses fonts either in the format used by Vibrant's ParseFont, or allows "small" "meduim", etc. . . */
static FonT LocalParseFont (
  CharPtr FontSpec
)

{
  if (FontSpec == NULL) return NULL;
  if (StringICmp (FontSpec, "small") == 0) {
    return SetSmallFont ();
  } else if (StringICmp (FontSpec, "medium") == 0) {
    return SetMediumFont ();
  } else if (StringICmp (FontSpec, "large") == 0) {
    return SetLargeFont ();
  } else if (StringICmp (FontSpec, "system") == 0) {
    return systemFont;
  } else if (StringICmp (FontSpec, "program") == 0) {
    return programFont;
  } else {
    return ParseFont (FontSpec);
  }
}

/* COLOR must point to an array of Uint1 [3]; */
static Boolean ParseColor (
  CharPtr string,
  Uint1Ptr color
)

{
  unsigned  sscanfOffset, localColor [3] = { 0 }; /* "unsigned", to match %u and %n*/
  Uint1  offset;
  Uint1  i;
  Boolean isLight = FALSE, isDark = FALSE;
  Char modifierBuffer [8];


  if (string == NULL || color == NULL) return FALSE;
  /* first try to parse a human-readable color (ie, light blue)*/
  StringNCpy_0 (modifierBuffer, string, sizeof (modifierBuffer));
  modifierBuffer [5] = '\0'; /* truncate "light xyz" -> "light" */
  if (StringICmp (modifierBuffer, "light") == 0) {
    isLight = TRUE;
    offset = 5;
  } else {
    modifierBuffer [4] = '\0';
    if (StringICmp (modifierBuffer, "dark") == 0) {
      isDark = TRUE;
      offset = 4;
    } else {
      modifierBuffer [3] = '\0';
      if (StringICmp (modifierBuffer, "drk") == 0) {
        isDark = TRUE;
        offset = 3;
      } else {
        modifierBuffer [2] = '\0';
        if (StringICmp (modifierBuffer, "lt") == 0) {
          isLight = TRUE;
          offset = 2;
        } else if (StringICmp (modifierBuffer, "dk") == 0) {
          isDark = TRUE;
          offset = 2;
        } else {
          offset = 0;
        }
      }
    }
  }
  while (string[offset] != '\0' && !IS_ALPHA (string[offset])) {
    offset ++;
  }
  i = StringIndexInStringList (string + offset, ColorStrings);
  if (i > 0 && i < DIM (ColorValues)) {
    color [0] = ColorValues [i] [0];
    color [1] = ColorValues [i] [1];
    color [2] = ColorValues [i] [2];
    if (isLight) {
      color [0] = MIN (((3 * color [0]) / 2), 255);
      color [1] = MIN (((3 * color [1]) / 2), 255);
      color [2] = MIN (((3 * color [2]) / 2), 255);
    } else if (isDark) {
      color [0] = (2 * color [0]) / 3;
      color [1] = (2 * color [1]) / 3;
      color [2] = (2 * color [2]) / 3;
    }
  } else {
    offset = 0;
    for (i = 0; i < 3; i++) {
      for (; string[offset] != '\0' && !IS_DIGIT (string [offset]); offset++) continue;
      if (string [offset] == '\0') return FALSE;
      if (sscanf (string + offset, "%u%n", localColor + i, &sscanfOffset) == 0)  return FALSE;
      offset += sscanfOffset;
    }
    color [0] = localColor [0];
    color [1] = localColor [1];
    color [2] = localColor [2];
  }
  return TRUE;
}


NLM_EXTERN FilterPtr FindFilterByName (
  CharPtr name,
  ViewerConfigsPtr VCP
)

{
  Uint1  i;

  if (VCP == NULL || StringHasNoText (name) || !VCP->ArraysPopulated) return NULL;
  i = StringIndexInStringList (name, VCP->FilterNameArray);
  if (i < VCP->FilterCount && i >= 0) return (VCP->FilterArray[i]);
  return NULL;
}

NLM_EXTERN AppearancePtr FindAppearanceByName (
  CharPtr name,
  ViewerConfigsPtr VCP
)

{
  Uint1  i;

  if (VCP == NULL || StringHasNoText (name) || !VCP->ArraysPopulated) return NULL;
  i = StringIndexInStringList (name, VCP->AppearanceNameArray);
  if (i < VCP->AppearanceCount && i>= 0) return (VCP->AppearanceArray[i]);
  return NULL;
}

NLM_EXTERN LayoutAlgorithm FindLayoutByName (
  CharPtr name
)

{
  Uint1  i;

  if (StringHasNoText (name)) return 0;
  i = StringIndexInStringList (name, LayoutStrings);
  if (i >= 0 && i < DIM(LayoutValues)) {
    return LayoutValues [i];
  }
  return Layout_Inherit;
}

static FilterPtr FindFilterByName_T (
  CharPtr name,
  ViewerConfigsPtr VCP
)

{
  ValNodePtr  nameVNP, filtVNP;

  if (VCP == NULL || StringHasNoText (name)) return NULL;
  for (nameVNP = VCP->FilterNameList, filtVNP = VCP->FilterList;
       nameVNP != NULL && filtVNP != NULL;
       nameVNP = nameVNP->next, filtVNP = filtVNP->next) {
    if (nameVNP->data.ptrvalue != NULL &&
        StringICmp (name, nameVNP->data.ptrvalue) == 0) {
      return ((FilterPtr) filtVNP->data.ptrvalue);
    }
  }
  return NULL;
}

static AppearancePtr FindAppearanceByName_T (
  CharPtr name,
  ViewerConfigsPtr VCP
)

{
  ValNodePtr  nameVNP, appVNP;

  if (VCP == NULL || StringHasNoText (name)) return NULL;
  for (nameVNP = VCP->AppearanceNameList, appVNP = VCP->AppearanceList;
       nameVNP != NULL && appVNP != NULL;
       nameVNP = nameVNP->next, appVNP = appVNP->next) {
    if (nameVNP->data.ptrvalue != NULL &&
        StringICmp (name, nameVNP->data.ptrvalue) == 0) {
      return ((AppearancePtr) appVNP->data.ptrvalue);
    }
  }
  return NULL;
}

NLM_EXTERN AppearancePtr CreateAppearance (
  CharPtr newname,
  ViewerConfigsPtr VCP
)

{
  AppearancePtr   AP;

  if (VCP == NULL || StringHasNoText (newname)) return NULL;
  if (FindAppearanceByName_T (newname, VCP) != NULL) return NULL;  /* don't allow duplicate names */
  AP = (AppearancePtr) MemNew (sizeof (Appearance));
  if (AP == NULL) return NULL;
  AP->name = StringSave (newname);
  VCP->AppearanceCount++;
  ValNodeAddPointer (&VCP->AppearanceList, VCP->AppearanceCount, AP);
  ValNodeAddPointer (&VCP->AppearanceNameList, VCP->AppearanceCount, AP->name);
  return AP;
}

static BioseqAppearanceItemPtr ParseBioseqAppearanceItem (
  CharPtr sectionName,
  ViewerConfigsPtr VCP
)

{
  Char                     inputBuffer [128];
  BioseqAppearanceItemPtr  bioseqAIP;
  Int2                     i;
  unsigned                 val; /* "unsigned" to match sscanf("%ud")*/

  bioseqAIP = MemNew (sizeof (BioseqAppearanceItem));
  if (bioseqAIP == NULL) return NULL;

  bioseqAIP->labelLoc = GroupLabelLocationValues [0];
  if (GetAppParam (config_filename, sectionName, "label", NULL, inputBuffer, sizeof (inputBuffer))) {
    i = StringIndexInStringList (inputBuffer, GroupLabelLocations);
    if (i >= 0 && i < DIM (GroupLabelLocationValues)) {
      bioseqAIP->labelLoc = GroupLabelLocationValues  [i];
    }
  }

  bioseqAIP->drawScale = TRUE;
  if (GetAppParam (config_filename, sectionName, "scale", NULL, inputBuffer, sizeof (inputBuffer))) {
    i = StringIndexInStringList (inputBuffer, BoolStrings);
    if (i >= 0 && i < DIM (BoolValues)) {
      bioseqAIP->drawScale = BoolValues [i];
    }
  }

  if (GetAppParam (config_filename, sectionName, "labelfont", NULL, inputBuffer, sizeof (inputBuffer))) {
    bioseqAIP->labelFont = LocalParseFont (inputBuffer);
  }
  if (bioseqAIP->labelFont == NULL) {
    bioseqAIP->labelFont = systemFont;
  }

  if (GetAppParam (config_filename, sectionName, "scalefont", NULL, inputBuffer, sizeof (inputBuffer))) {
    bioseqAIP->scaleFont = LocalParseFont (inputBuffer);
  }
  if (bioseqAIP->scaleFont == NULL) {
    bioseqAIP->scaleFont = SetSmallFont ();
  }

  if (GetAppParam (config_filename, sectionName, "height", NULL, inputBuffer, sizeof (inputBuffer))) {
    if (inputBuffer != NULL) {
      sscanf (inputBuffer, "%ud", &val);
      bioseqAIP->height = MIN (val, 16);
    }
  }
  if (bioseqAIP->height == 0) {
    bioseqAIP->height = 10;
  }

  if (GetAppParam (config_filename, sectionName, "scaleheight", NULL, inputBuffer, sizeof (inputBuffer))) {
    if (inputBuffer != NULL) {
      sscanf (inputBuffer, "%ud", &val);
      bioseqAIP->scaleHeight = MIN (val, 16);
    }
  }
  if (bioseqAIP->scaleHeight == 0) {
    bioseqAIP->scaleHeight = 10;
  }

  if (GetAppParam (config_filename, sectionName, "color", NULL, inputBuffer, sizeof (inputBuffer))) {
    ParseColor (inputBuffer, bioseqAIP->bioseqColor);
  }
  if (GetAppParam (config_filename, sectionName, "labelcolor", NULL, inputBuffer, sizeof (inputBuffer))) {
    ParseColor (inputBuffer, bioseqAIP->labelColor);
  }
  if (GetAppParam (config_filename, sectionName, "scalecolor", NULL, inputBuffer, sizeof (inputBuffer))) {
    ParseColor (inputBuffer, bioseqAIP->scaleColor);
  }

  bioseqAIP->format = BioseqFormatValues [0];
  if (GetAppParam (config_filename, sectionName, "format", NULL, inputBuffer, sizeof (inputBuffer))) {
    i = StringIndexInStringList (inputBuffer, BioseqFormat);
    if (i >= 0 && i < DIM (BioseqFormatValues)) {
      bioseqAIP->format = BioseqFormatValues  [i];
    }
  }

  return bioseqAIP;
}

static AppearanceItemPtr ParseFeatureAppearanceItem (
  CharPtr sectionName,
  AppearanceItemPtr inheritFromMe,
  Boolean recursing,
  ViewerConfigsPtr VCP
)

{
  Char               inputBuffer [128];
  AppearanceItemPtr  newAIP;
  Int2               i;
  Boolean            changed = FALSE;
  unsigned           val; /* "unsigned" to match sscanf("%ud")*/
  AppearanceItemPtr  namedAIP;

  if (! recursing) {
    if (GetAppParam (config_filename, sectionName, "usenamedstyle", NULL, inputBuffer, sizeof (inputBuffer))) {
      namedAIP = ParseFeatureAppearanceItem (inputBuffer, inheritFromMe, TRUE, VCP);
      if (namedAIP != NULL) {
        inheritFromMe = namedAIP;
        changed = TRUE; /* !!! this will use more memory than necessary */
      }
    }
  }
  newAIP = MemNew (sizeof (AppearanceItem));
  if (newAIP == NULL) return NULL;
  MemCopy (newAIP, inheritFromMe, sizeof (AppearanceItem));
  if (GetAppParam (config_filename, sectionName, "color", NULL, inputBuffer, sizeof (inputBuffer))) {
    changed = ParseColor (inputBuffer, newAIP->Color);
  }
  if (GetAppParam (config_filename, sectionName, "labelcolor", NULL, inputBuffer, sizeof (inputBuffer))) {
    changed = ParseColor (inputBuffer, newAIP->LabelColor);
  }

  newAIP->LabelLoc = LabelAbove;
  if (GetAppParam (config_filename, sectionName, "label", NULL, inputBuffer, sizeof (inputBuffer))) {
    i = StringIndexInStringList (inputBuffer, LlocStrings);
    if (i >= 0 && i < DIM (LlocValues)) {
      newAIP->LabelLoc = LlocValues [i];
    }
  }
  if (GetAppParam (config_filename, sectionName, "displaywith", NULL, inputBuffer, sizeof (inputBuffer))) {
    i = StringIndexInStringList (inputBuffer, RenderStrings);
    if (i >= 0 && i < DIM (RenderValues)) {
      newAIP->RenderChoice = RenderValues [i];
      changed = TRUE;
    }
  }
  if (GetAppParam (config_filename, sectionName, "showarrow", NULL, inputBuffer, sizeof (inputBuffer))) {
    i = StringIndexInStringList (inputBuffer, BoolStrings);
    if (i >= 0 && i < DIM (BoolValues)) {
      newAIP->ShowArrow = BoolValues [i];
      changed = TRUE;
    }
  }
  if (GetAppParam (config_filename, sectionName, "gap", NULL, inputBuffer, sizeof (inputBuffer))) {
    i = StringIndexInStringList (inputBuffer, GapStrings);
    if (i >= 0 && i < DIM (RenderValues)) {
      newAIP->GapChoice = GapValues [i];
      changed = TRUE;
    }
  }
  if (GetAppParam (config_filename, sectionName, "showtype", NULL, inputBuffer, sizeof (inputBuffer))) {
    i = StringIndexInStringList (inputBuffer, BoolStrings);
    if (i >= 0 && i < DIM (BoolValues)) {
      newAIP->AddTypeToLabel = BoolValues [i];
    }
  }
  if (GetAppParam (config_filename, sectionName, "showcontent", NULL, inputBuffer, sizeof (inputBuffer))) {
    i = StringIndexInStringList (inputBuffer, BoolStrings);
    if (i >= 0 && i < DIM (BoolValues)) {
      newAIP->AddDescToLabel = BoolValues [i];
    }
  }
  if (GetAppParam (config_filename, sectionName, "height", NULL, inputBuffer, sizeof (inputBuffer))) {
    if (inputBuffer != NULL) {
      sscanf (inputBuffer, "%ud", &val);
      newAIP->Height = MIN (val, 16);
      changed = TRUE;
    }
  }
  if (GetAppParam (config_filename, sectionName, "labelfont", NULL, inputBuffer, sizeof (inputBuffer))) {
    newAIP->LabelFont = LocalParseFont (inputBuffer);
  }
  if (newAIP->LabelFont == NULL) {
    newAIP->LabelFont = programFont;
  }
  if (newAIP->LabelFont != inheritFromMe->LabelFont) {
    changed = TRUE;
  }
  if (! changed) {
    MemFree (newAIP);
    return NULL;
  }
  return newAIP;
}

static AppearancePtr ParseAppearance (
  CharPtr appearanceNameInFile,
  ViewerConfigsPtr VCP
)

{
  AppearancePtr      AP;
  AppearanceItemPtr  AIP, impAIP, newAIP;
  Char               inputBuffer [128];
  Char               sectionName [128];
  Char               outputBuffer [128];
  AppearanceItem     DefaultAppearanceItem = {
    {0, 0, 0}, {64, 64, 64}, Render_Box, 0, 5, 0, FALSE, LineGap, TRUE, TRUE, NULL, LabelAbove
  };
  Uint1              i;
  unsigned           val;

  if (appearanceNameInFile == NULL) return NULL;
  DefaultAppearanceItem.LabelFont = programFont;
  sprintf (sectionName, "%s.master", appearanceNameInFile);
  /* require all styles to have a name, since high-level interface uses the name to identify Filters */
  if (! GetAppParam (config_filename, sectionName, "name", NULL, inputBuffer, sizeof (inputBuffer))) return NULL;
  if (StringHasNoText (inputBuffer)) return NULL;
  AP = CreateAppearance (inputBuffer, VCP);
  if (AP == NULL) return NULL;
  AIP = ParseFeatureAppearanceItem (sectionName, &DefaultAppearanceItem, FALSE, VCP); /*parse xyz.master */
  if (AIP == NULL) {            /* require a "master" style */
    DestroyAppearance (AP, VCP);
    return NULL;
  }
  val = VCP->DefaultMaxScaleForArrow;
  if (GetAppParam (config_filename, sectionName, "maxarrowscale", NULL, inputBuffer, sizeof (inputBuffer))) {
    sscanf (inputBuffer, "%ud", &val);
  }
  AP->MaxScaleForArrow = val;
  val = VCP->DefaultMinPixelsForArrow;
  if (GetAppParam (config_filename, sectionName, "minarrowpixels", NULL, inputBuffer, sizeof (inputBuffer))) {
    sscanf (inputBuffer, "%ud", &val);
  }
  AP->MinPixelsForArrow = val;
  AP->ShadeSmears = VCP->DefaultShadeSmears;
  if (GetAppParam (config_filename, sectionName, "shadesmears", NULL, inputBuffer, sizeof (inputBuffer))) {
    i = StringIndexInStringList (inputBuffer, BoolStrings);
    if (i >= 0 && i < DIM(BoolValues) && i >= 0) {
      AP->ShadeSmears = BoolValues [i];
    }
  }

  if (GetAppParam (config_filename, sectionName, "groupboxcolor", NULL, inputBuffer, sizeof (inputBuffer))) {
    ParseColor (inputBuffer, AP->GroupBoxColor);
  }

  if (GetAppParam (config_filename, sectionName, "grouplabelfont", NULL, inputBuffer, sizeof (inputBuffer))) {
    AP->GroupLabelFont = LocalParseFont (inputBuffer);
    if (AP->GroupLabelFont == NULL) {
      AP->GroupLabelFont = programFont;
    }
  }
  if (GetAppParam (config_filename, sectionName, "grouplabelcolor", NULL, inputBuffer, sizeof (inputBuffer))) {
    ParseColor (inputBuffer, AP->GroupLabelColor);
  }

  sprintf (outputBuffer, "%s.bioseq", appearanceNameInFile);
  AP->bioseqAIP = ParseBioseqAppearanceItem (outputBuffer, VCP);
  sprintf (outputBuffer, "%s.imp", appearanceNameInFile);
  impAIP = ParseFeatureAppearanceItem (outputBuffer, AIP, FALSE, VCP);
  if (impAIP == NULL) {
    impAIP = AIP;
  } else {
    AddAppearanceItemToAppearance (impAIP, AP, FEATDEF_IMP, VCP);
  }
  for (i = 1; i < FEATDEF_MAX; i++) {
    if (i == FEATDEF_IMP) continue;
    sprintf (outputBuffer, "%s.%s", appearanceNameInFile, FindKeyFromFeatDefType (i, FALSE));
    if (i >= FEATDEF_allele && i <= FEATDEF_35_signal) {        /* is it an imp-feat ? */
      newAIP = ParseFeatureAppearanceItem (outputBuffer, impAIP, FALSE, VCP);
      newAIP = newAIP ? newAIP : impAIP;
    } else {
      newAIP = ParseFeatureAppearanceItem (outputBuffer, AIP, FALSE, VCP);
      newAIP = newAIP ? newAIP : AIP;
    }
    if (newAIP != NULL) {
      AddAppearanceItemToAppearance (newAIP, AP, i, VCP);
    }
  }
  return AP;
}

NLM_EXTERN FilterPtr CreateFilter (
  CharPtr name,
  ViewerConfigsPtr VCP
)

{
  FilterPtr  FP;

  if (VCP == NULL || StringHasNoText (name)) return NULL;
  if (FindFilterByName_T (name, VCP) != NULL) return NULL;  /* don't allow duplicate names */
  FP = MemNew (sizeof (Filter));
  if (FP == NULL) return FP;
  FP->name = StringSave (name);
  VCP->FilterCount++;
  ValNodeAddPointer (&VCP->FilterList, VCP->FilterCount, FP);
  ValNodeAddPointer (&VCP->FilterNameList, VCP->FilterCount, FP->name);
  return FP;
}

static void ChangeFeatureInFilterItem (
  FilterItemPtr FIP,
  Uint1 newFeatdef,
  Boolean includeMe,
  ViewerConfigsPtr VCP
)

{
  Uint1  i;
  Uint1  order = 0;
  Uint1  orderIncrement = 0;

  if (FIP == NULL) return;

  if (includeMe) {
    orderIncrement = 1;
    for (i = 0; i < FEATDEF_MAX; i++) {
      order = MAX (order, FIP->IncludeFeature [i]);
    }
    order++; /* the lowest available order index */
  }

  if (newFeatdef == FEATDEF_ANY) {
    for (i = 0; i < FEATDEF_MAX; i++) {
      FIP->IncludeFeature [i] = order;
      order += orderIncrement;
    }
  } else if (newFeatdef == FEATDEF_ANY_RNA) {
    for (i = FEATDEF_preRNA; i <= FEATDEF_otherRNA; i++) {
      FIP->IncludeFeature [i] = order;
      order += orderIncrement;
    }
    FIP->IncludeFeature [FEATDEF_misc_RNA] = order;
    order += orderIncrement;
    FIP->IncludeFeature [FEATDEF_precursor_RNA] = order;
    order += orderIncrement;
    FIP->IncludeFeature [FEATDEF_snoRNA] = order;
  } else if (newFeatdef == FEATDEF_ANY_PROT) {
    FIP->IncludeFeature [FEATDEF_PROT] = order;
    order += orderIncrement;
    for (i = FEATDEF_preprotein; i <= FEATDEF_transit_peptide_aa; i++) {
      FIP->IncludeFeature [i] = order;
      order += orderIncrement;
    }
  } else if (newFeatdef == FEATDEF_IMP) {
    for (i = FEATDEF_allele; i <= FEATDEF_35_signal; i++) {
      FIP->IncludeFeature [i] = order;
      order += orderIncrement;
    }
  } else if (newFeatdef < FEATDEF_MAX) {
    FIP->IncludeFeature [newFeatdef] = order;
  }
}

NLM_EXTERN void AddFeatureToFilterItem (
  FilterItemPtr FIP,
  Uint1 newFeatdef,
  ViewerConfigsPtr VCP
)

{
  ChangeFeatureInFilterItem (FIP, newFeatdef, TRUE, VCP);
}

NLM_EXTERN void RemoveFeatureFromFilterItem (
  FilterItemPtr FIP,
  Uint1 newFeatdef,
  ViewerConfigsPtr VCP
)

{
  ChangeFeatureInFilterItem (FIP, newFeatdef, FALSE, VCP);
}

static void AddFilterItemToFilter (
 FilterItemPtr newFIP,
  FilterPtr parent,
  ViewerConfigsPtr VCP
)

{
  ValNodePtr  newVNP, lastVNP;

  for (lastVNP = parent->FilterItemList;
       lastVNP != NULL && lastVNP->next != NULL;
       lastVNP = lastVNP->next) continue;
  newVNP = ValNodeAdd (&parent->FilterItemList);
  newVNP->data.ptrvalue = newFIP;
}

NLM_EXTERN FilterItemPtr CreateNewFilterItemInFilter (
  CharPtr name,
  FilterPtr parent,
  ViewerConfigsPtr VCP
)

{
  FilterItemPtr  FIP;

  FIP = MemNew (sizeof (FilterItem));
  if (FIP == NULL) return FIP;
  FIP->GroupLabel = StringSaveNoNull (name);
  AddFilterItemToFilter (FIP, parent, VCP);
  return FIP;
}

static FilterItemPtr ParseFilterItem (
  CharPtr filterItemName,
  Uint2 defaultRowPadding,
  Uint2 defaultGroupPadding,
  LayoutAlgorithm defaultLayout,
  ViewerConfigsPtr VCP
)

{
  Char           sectionName [128];
  Char           featureNum [128];
  Char           inputBuffer [128];
  Uint4          featdeftype;
  Uint4          featureCount = 0;
  FilterItemPtr  FIP;
  Int2           i;
  unsigned       val;

  FIP = MemNew (sizeof (FilterItem));
  if (FIP == NULL) return FIP;
  FIP->DrawScale = TristateUnset;
  FIP->Type = FeatureFilter; /* this will get changed if a graph or alignment is disconvered instead */
  sprintf (sectionName, "filters.%s", filterItemName);
  GetAppParam (config_filename, sectionName, "layout", "inherit", inputBuffer, sizeof (inputBuffer));
  i = StringIndexInStringList (inputBuffer, LayoutStrings);
  if (i >=0 && i < DIM (LayoutStrings) && i >= 0) {
    FIP->LayoutChoice = LayoutValues [i];
  } else {
    FIP->LayoutChoice = defaultLayout;
  }

  FIP->GroupPadding = defaultGroupPadding;
  if (GetAppParam (config_filename, sectionName, "grouppadding", NULL, inputBuffer, sizeof (inputBuffer))) {
    sscanf (inputBuffer, "%ud", &val);
    val = MIN (val, 100);
    FIP->GroupPadding = val;
  }

  FIP->IntraRowPaddingPixels = defaultRowPadding;
  if (GetAppParam (config_filename, sectionName, "rowpadding", NULL, inputBuffer, sizeof (inputBuffer))) {
    sscanf (inputBuffer, "%ud", &val);
    val = MIN (val, 100);
    FIP->IntraRowPaddingPixels = val;
  }
  FIP->DrawItemRect = FALSE;
  FIP->FillItemRect = FALSE;
  if (GetAppParam (config_filename, sectionName, "groupbox", NULL, inputBuffer, sizeof (inputBuffer))) {
    i = StringIndexInStringList (inputBuffer, BoolStrings);
    if (i >= 0 && i < DIM (BoolValues) && i >=0) {
      FIP->DrawItemRect = BoolValues [i];
    }
  }
  if (FIP->DrawItemRect) {
    if (GetAppParam (config_filename, sectionName, "groupboxcolor", NULL, inputBuffer, sizeof (inputBuffer))) {
      FIP->GroupBoxColorSet = TRUE;
      ParseColor (inputBuffer, FIP->GroupBoxColor);
    }
    if (GetAppParam (config_filename, sectionName, "fillbox", NULL, inputBuffer, sizeof (inputBuffer))) {
      i = StringIndexInStringList (inputBuffer, BoolStrings);
      if (i >= 0 && i < DIM (BoolValues) && i >= 0) {
        FIP->FillItemRect = BoolValues[i];
      }
    }
  }

  FIP->GroupLabel = NoLabel;
  FIP->GroupLabelFont = programFont;
  if (GetAppParam (config_filename, sectionName, "name", NULL, inputBuffer, sizeof (inputBuffer))) {
    FIP->GroupLabel = StringSaveNoNull (inputBuffer);
    FIP->GroupLabelLoc = LabelOnTop;
    if (GetAppParam (config_filename, sectionName, "grouplabel", NULL, inputBuffer, sizeof (inputBuffer))) {
      i = StringIndexInStringList (inputBuffer, GroupLabelLocations);
      if (i >= 0 && i < DIM (GroupLabelLocationValues) && i >= 0) {
        FIP->GroupLabelLoc = GroupLabelLocationValues[i];
      }
    }
    if (GetAppParam (config_filename, sectionName, "grouplabelfont", NULL, inputBuffer, sizeof (inputBuffer))) {
      FIP->GroupLabelFontSet = TRUE;
      FIP->GroupLabelFont = LocalParseFont (inputBuffer);
      if (FIP->GroupLabelFont == NULL) {
        FIP->GroupLabelFont = programFont;
      }
    }
    if (GetAppParam (config_filename, sectionName, "grouplabelcolor", NULL, inputBuffer, sizeof (inputBuffer))) {
      FIP->GroupLabelColorSet = TRUE;
      ParseColor (inputBuffer, FIP->GroupLabelColor);
    }
  }

  FIP->LabelLoc = LabelUnset;
  if (GetAppParam (config_filename, sectionName, "label", NULL, inputBuffer, sizeof (inputBuffer))) {
    i = StringIndexInStringList (inputBuffer, LlocStrings);
    if (i >= 0 && i < DIM (LlocValues)) {
      FIP->LabelLoc = LlocValues [i];
    }
  }

  FIP->AddTypeToLabel = TristateUnset;
  if (FIP->LabelLoc != LabelNone && GetAppParam (config_filename, sectionName, "showtype", NULL, inputBuffer, sizeof (inputBuffer))) {
    i = StringIndexInStringList (inputBuffer, BoolStrings);
    if (i >= 0 && i < DIM (BoolValues)) {
      FIP->AddTypeToLabel = BOOL_TO_TRISTATE (BoolValues [i]);
    }
  }
  FIP->AddDescToLabel = TristateUnset;
  if (FIP->LabelLoc != LabelNone && GetAppParam (config_filename, sectionName, "showcontent", NULL, inputBuffer, sizeof (inputBuffer))) {
    i = StringIndexInStringList (inputBuffer, BoolStrings);
    if (i >= 0 && i < DIM (BoolValues)) {
      FIP->AddDescToLabel = BOOL_TO_TRISTATE (BoolValues [i]);
    }
  }

  FIP->MatchStrand = StrandValues [0];
  if (GetAppParam (config_filename, sectionName, "strand", NULL, inputBuffer, sizeof (inputBuffer))) {
    i = StringIndexInStringList (inputBuffer, StrandStrings);
    if (i >= 0 && i < DIM (StrandValues)) {
      FIP->MatchStrand = StrandValues [i];
    }
  }
  for (i = 1; i < FEATDEF_MAX; i++) {
    sprintf (featureNum, "feature%d", (unsigned) i);
    if (GetAppParam (config_filename, sectionName, featureNum, NULL, inputBuffer, sizeof (inputBuffer))) {
      featdeftype = FindFeatDefTypeFromKey (inputBuffer);
      if (featdeftype == FEATDEF_BAD) {
        /* insert special-case checks for types of features not found by FindFeatDefTypeFromKey () here */
        if (StringICmp (inputBuffer, "everything") == 0 ||
            StringICmp (inputBuffer, "all") == 0 ||
            StringICmp (inputBuffer, "every") == 0 ||
            StringICmp (inputBuffer, "any") == 0) {
          featdeftype = FEATDEF_ANY;
        } else if (StringICmp (inputBuffer, "rna") == 0) {
          featdeftype = FEATDEF_ANY_RNA;
        } else if (StringICmp (inputBuffer, "prot") == 0) {
          featdeftype = FEATDEF_ANY_PROT;
        } else if (StringICmp (inputBuffer, "bioseq") == 0) {
          FIP->Type = BioseqFilter;
        } else if (StringICmp (inputBuffer, "graph") == 0) {
          FIP->Type = GraphFilter;
        } else if (StringICmp (inputBuffer, "align") == 0) {
          FIP->Type = AlignmentFilter;
        } else continue; /* failed to find a match */
      }
      AddFeatureToFilterItem (FIP, featdeftype, VCP);
      featureCount++;
    }
  }
  if (FIP->Type == BioseqFilter) {
    FIP->DrawScale = TristateUnset;
    if (GetAppParam (config_filename, sectionName, "scale", NULL, inputBuffer, sizeof (inputBuffer))) {
      i = StringIndexInStringList (inputBuffer, BoolStrings);
      if (i >= 0 && i < DIM (BoolValues)) {
        FIP->DrawScale = BOOL_TO_TRISTATE (BoolValues [i]);
      }
    }
  }
  if (featureCount == 0) {
    MemFree (FIP);
    return NULL;
  }
  return FIP;
}

static FilterPtr ParseFilter (
  CharPtr filterNameInFile,
  ViewerConfigsPtr VCP
)

{

  FilterPtr       FP;
  FilterItemPtr   FIP;
  Int2            i;
  Uint1           filterItemCount = 0;
  Char            inputBuffer [128];     /* for input *from* GetAppParam */
  Char            outputBuffer [128];    /* paramater *to* GetAppParam */
  Char            sectionName [128];
  Boolean         foundBioseqFilter = FALSE;
  Boolean         foundGraphFilter = FALSE;
  Boolean         foundAlignmentFilter = FALSE;
  ValNodePtr      VNP;
  Boolean         createImplicitBioseq = TRUE;
  Boolean         createImplicitGraphs = TRUE;
  unsigned        val; /* to match sscanf ("%d"...)*/
  Uint2           defaultRowPadding;
  Uint2           defaultGroupPadding;
  LayoutAlgorithm defaultLayout;

  if (filterNameInFile == NULL) return NULL;
  sprintf (sectionName, "%s", filterNameInFile);
  /* require all styles to have a name, since high-level interface uses the name to identify Filters */
  if (! GetAppParam (config_filename, sectionName, "name", NULL, inputBuffer, sizeof (inputBuffer))) return NULL;
  FP = CreateFilter (inputBuffer, VCP); /* Createfilter will check for duplucate names */
  if (FP == NULL) return FP;
  val = VCP->DefaultMaxScaleWithLabels;
  if (GetAppParam (config_filename, sectionName, "maxlabelscale", NULL, inputBuffer, sizeof (inputBuffer))) {
    sscanf (inputBuffer, "%ud", &val);
  }
  FP->MaxScaleWithLabels = val;

  GetAppParam (config_filename, sectionName, "layout", NULL, inputBuffer, sizeof (inputBuffer));
  i = StringIndexInStringList (inputBuffer, LayoutStrings);
  if (i >= 0 && i < DIM (LayoutValues) && i >= 0) {
    defaultLayout = LayoutValues [i];
  } else {
    defaultLayout = Layout_Inherit;
  }

  val = VCP->DefaultGroupPadding;
  if (GetAppParam (config_filename, sectionName, "grouppadding", NULL, inputBuffer, sizeof (inputBuffer))) {
    sscanf (inputBuffer, "%ud", &val);
    val = MIN (val, 100);
  }
  defaultGroupPadding = val;

  val = VCP->DefaultRowPadding;
  if (GetAppParam (config_filename, sectionName, "rowpadding", NULL, inputBuffer, sizeof (inputBuffer))) {
    sscanf (inputBuffer, "%ud", &val);
  }
  defaultRowPadding = val;

  if (GetAppParam (config_filename, sectionName, "suppressbioseq", NULL, inputBuffer, sizeof (inputBuffer))) {
    i = StringIndexInStringList (inputBuffer, BoolStrings);
    if (i >= 0 && i < DIM (BoolStrings) && i >= 0) {
      createImplicitBioseq = ! (BoolValues [i]);
    }
  }
  if (GetAppParam (config_filename, sectionName, "suppressgraphs", NULL, inputBuffer, sizeof (inputBuffer))) {
    i = StringIndexInStringList (inputBuffer, BoolStrings);
    if (i >= 0 && i < DIM (BoolStrings) && i >= 0) {
      createImplicitGraphs = ! (BoolValues [i]);
    }
  }

  for (i = 1; i < FEATDEF_MAX; i++) {
    sprintf (outputBuffer, "%s%d", "group", (unsigned) i);
    if (GetAppParam (config_filename, sectionName, outputBuffer, NULL, inputBuffer, sizeof (inputBuffer))) {
      FIP = ParseFilterItem (inputBuffer, defaultRowPadding, defaultGroupPadding, defaultLayout, VCP);
      if (FIP != NULL && FIP->Type == BioseqFilter) {
        foundBioseqFilter = TRUE;
      }
      if (FIP != NULL && FIP->Type == GraphFilter) {
        foundGraphFilter = TRUE;
      }
      if (FIP != NULL && FIP->Type == AlignmentFilter) {
        foundAlignmentFilter = TRUE;
      }
      AddFilterItemToFilter (FIP, FP, VCP);
      filterItemCount++;
    }
  }
  if (filterItemCount == 0) {
    DestroyFilter (FP, VCP);
    return NULL;
  }

  if (createImplicitBioseq && ! foundBioseqFilter) {
    VNP = MemNew (sizeof (ValNode));
    FIP = MemNew (sizeof (FilterItem));
    if (VNP == NULL || FIP == NULL) {
      DestroyFilter (FP, VCP);
      return NULL;
    }
    /* insert a Bioseq filter at the head of the list */
    VNP->next = FP->FilterItemList;
    VNP->data.ptrvalue = FIP;
    FP->FilterItemList = VNP;

    FIP->Type = BioseqFilter;
    FIP->IntraRowPaddingPixels = 5;
  }

  if (createImplicitGraphs && ! foundAlignmentFilter) {
    /* insert a Graph filter at the_end of the list */
    FIP = MemNew (sizeof (FilterItem));
    VNP = ValNodeAddPointer (&FP->FilterItemList, 0, FIP);
    if (FIP == NULL || VNP == NULL ) {
      DestroyFilter (FP, VCP);
      return NULL;
    }
    FIP->Type = AlignmentFilter;
    FIP->IntraRowPaddingPixels = 5;
  }

  if (createImplicitGraphs && ! foundGraphFilter) {
    /* insert a Graph filter at the_end of the list */
    FIP = MemNew (sizeof (FilterItem));
    VNP = ValNodeAddPointer (&FP->FilterItemList, 0, FIP);
    if (FIP == NULL || VNP == NULL ) {
      DestroyFilter (FP, VCP);
      return NULL;
    }
    FIP->Type = GraphFilter;
    FIP->IntraRowPaddingPixels = 5;
  }

  return FP;
}

/* if this will be used by multiple threads in a multi-threaded application, there should be a lock around writing this */
static ViewerConfigsPtr newGraphicViewer_ConfigFileParse_Global = NULL;
static void InitializeDefaultStyle (
  CharPtr configFileName
);

NLM_EXTERN ViewerConfigsPtr GetGraphicConfigParseResults (
  void
)

{
  Uint1             AppearanceCount;
  Uint1             FilterCount;
  Uint1             i;
  void PNTR PNTR    ptr2;
  ViewerConfigsPtr  VCP;
  ValNodePtr        nameVNP;
  ValNodePtr        VNP;

  if (newGraphicViewer_ConfigFileParse_Global != NULL)
    return newGraphicViewer_ConfigFileParse_Global;

  InitializeDefaultStyle (config_filename);

  VCP = MemNew (sizeof (ViewerConfigs));
  if (VCP == NULL) return NULL;

  if (ParseConfigFile (VCP) == 0) return NULL;
  /* this should never happen, because of the static default style*/
  if (VCP->AppearanceCount == 0 || VCP->FilterCount == 0) return NULL;

  AppearanceCount = VCP->AppearanceCount;
  FilterCount = VCP->FilterCount;
  i = (AppearanceCount + FilterCount) * 2 + 2; /* total number of pointers needed (the extra 2 are NULL's to terminate the name lists)*/
  ptr2 = MemNew (i * sizeof (void PNTR));
  if (ptr2 == NULL) {
    MemFree (VCP);
    return NULL;
  }
  VCP->AppearanceArray = (AppearancePtr PNTR) ptr2;
  ptr2 += AppearanceCount;
  VCP->AppearanceNameArray = (CharPtr PNTR) ptr2;
  ptr2 += AppearanceCount + 1; /* add a NULL pointer to terminate the list */
  VCP->FilterArray = (FilterPtr PNTR) ptr2;
  ptr2 += FilterCount;
  VCP->FilterNameArray = (CharPtr PNTR) ptr2;
  VNP = VCP->AppearanceList;
  nameVNP = VCP->AppearanceNameList;
  for (i = 0; i < AppearanceCount; i++) {
    VCP->AppearanceArray[i] = VNP->data.ptrvalue;
    VCP->AppearanceNameArray[i] = nameVNP->data.ptrvalue;
    VNP = VNP->next;
    nameVNP = nameVNP->next;
  }

  VNP = VCP->FilterList;
  nameVNP = VCP->FilterNameList;
  for (i = 0; i < FilterCount; i++) {
    VCP->FilterArray[i] = VNP->data.ptrvalue;
    VCP->FilterNameArray[i] = nameVNP->data.ptrvalue;
    VNP = VNP->next;
    nameVNP = nameVNP->next;
  }
  VCP->ArraysPopulated = TRUE;
  newGraphicViewer_ConfigFileParse_Global = VCP;
  return VCP;
}

/* returns count of objects successfully parsed -- so 0 on failure*/
NLM_EXTERN Uint2 ParseConfigFile (
  ViewerConfigsPtr VCP
)

{
  Char     tagBuffer [32];
  Char     nameBuffer [128];
  Int2     i;
  Uint2    fCount = 0, aCount = 0;
  VoidPtr  tempPtr;
  unsigned val; /* to match scanf("%ud"...) */

  GetAppParam (config_filename, "filters", "maxlabelscale", NULL, tagBuffer, sizeof (tagBuffer));
  if (sscanf (tagBuffer, "%ud", &val) != 1) {
    val = 200;
  }
  VCP->DefaultMaxScaleWithLabels = val;

  GetAppParam (config_filename, "filters", "grouppadding", NULL, tagBuffer, sizeof (tagBuffer));
  if (sscanf (tagBuffer, "%ud", &val) != 1) {
    val = 3;
  }
  VCP->DefaultGroupPadding = val;

  GetAppParam (config_filename, "filters", "rowpadding", NULL, tagBuffer, sizeof (tagBuffer));
  if (sscanf (tagBuffer, "%ud", &val) != 1) {
    val = 5;
  }
  VCP->DefaultRowPadding = val;

  GetAppParam (config_filename, "styles", "maxarrowscale", NULL, tagBuffer, sizeof (tagBuffer));
  if (sscanf (tagBuffer, "%ud", &val) != 1) {
    val = 5;
  }
  VCP->DefaultMaxScaleForArrow = val;

  GetAppParam (config_filename, "styles", "minarrowpixels", NULL, tagBuffer, sizeof (tagBuffer));
  if (sscanf (tagBuffer, "%ud", &val) != 1) {
    val = 5;
  }
  VCP->DefaultMinPixelsForArrow = val;

  VCP->DefaultShadeSmears = FALSE;
  if (GetAppParam (config_filename, "styles", "shadesmears", NULL, tagBuffer, sizeof (tagBuffer))) {
    i = StringIndexInStringList (tagBuffer, BoolStrings);
    if (i >= 0 && i < DIM (BoolValues)) {
      VCP->DefaultShadeSmears = BoolValues [i];
    }
  }

  for (i = 0; i < 110; i++) {   /* do filters first */
    if (i < 10) {
      sprintf (tagBuffer, "filter0%d", (unsigned) i);
    } else {
      sprintf (tagBuffer, "filter%d", (unsigned) i - 9);
    }
    if (GetAppParam (config_filename, "filters", tagBuffer, NULL, nameBuffer, sizeof (nameBuffer))) {
      tempPtr = ParseFilter (nameBuffer, VCP);
      if (tempPtr == NULL) continue;
      fCount++;
    }
  }
  for (i = 0; i < 110; i++) {
    if (i < 10) {
      sprintf (tagBuffer, "style0%d", (unsigned) i);
    } else {
      sprintf (tagBuffer, "style%d", (unsigned) i - 9);
    }
    if (GetAppParam (config_filename, "styles", tagBuffer, NULL, nameBuffer, sizeof (nameBuffer))) {
      tempPtr = ParseAppearance (nameBuffer, VCP);
      if (tempPtr == NULL) continue;
      aCount++;
    }
  }
  return (aCount + fCount);
}

NLM_EXTERN FilterPtr DestroyFilter (
  FilterPtr FP,
  ViewerConfigsPtr VCP
)

{
  FilterItemPtr  FIP;
  ValNodePtr     VNP;
  Uint1          i;

  if (FP == NULL || VCP == NULL) {
    return NULL;
  }
  for (VNP = FP->FilterItemList; VNP; VNP = VNP->next) {        /* free all filterItems, and their labels */
    FIP = (FilterItemPtr) VNP->data.ptrvalue;
    if (FIP == NULL) {
      continue;
    }
    MemFree (FIP->GroupLabel);
    MemFree (FIP);
  }
  for (VNP = VCP->FilterList; VNP != NULL; VNP = VNP->next) {
    if (VNP->data.ptrvalue == FP) {
      i = VNP->choice;
      VNP = ValNodeExtract (&VCP->FilterList, i);
      break;
    }
  }
  if (VNP != NULL) {
    MemFree (VNP);
    VNP = ValNodeExtract (&VCP->FilterNameList, i);
    MemFree (VNP->data.ptrvalue);
    MemFree (VNP);
  }
  MemFree (FP);
  return NULL;
}

NLM_EXTERN void AddAppearanceItemToAppearance (
  AppearanceItemPtr AIP,
  AppearancePtr AP,
  Uint1 newFeatdef,
  ViewerConfigsPtr VCP
)

{
  Uint1       i;
  ValNodePtr  VNP;

  if (AIP == NULL || AP == NULL || VCP == NULL) return;
  if (newFeatdef == FEATDEF_ANY) {
    for (i = 0; i < FEATDEF_MAX; i++) {
      AP->FeaturesAppearanceItem [i] = AIP;
    }
  } else if (newFeatdef == FEATDEF_ANY_RNA) {
    for (i = FEATDEF_preRNA; i <= FEATDEF_otherRNA; i++) {
      AP->FeaturesAppearanceItem [i] = AIP;
    }
    AP->FeaturesAppearanceItem [FEATDEF_misc_RNA] = AIP;
    AP->FeaturesAppearanceItem [FEATDEF_precursor_RNA] = AIP;
    AP->FeaturesAppearanceItem [FEATDEF_snoRNA] = AIP;
  } else if (newFeatdef == FEATDEF_ANY_PROT) {
    AP->FeaturesAppearanceItem [FEATDEF_PROT] = AIP;
    for (i = FEATDEF_preprotein; i <= FEATDEF_transit_peptide_aa; i++) {
      AP->FeaturesAppearanceItem [i] = AIP;
    }
  } else if (newFeatdef == FEATDEF_IMP) {
    for (i = FEATDEF_allele; i <= FEATDEF_35_signal; i++) {
      AP->FeaturesAppearanceItem [i] = AIP;
    }
  } else if (newFeatdef < FEATDEF_MAX) {
    AP->FeaturesAppearanceItem [newFeatdef] = AIP;
  } else return;
  for (VNP = AP->AppearanceItemList; VNP != NULL && VNP->data.ptrvalue != AIP; VNP = VNP->next) continue;
  if (! (VNP != NULL && VNP->data.ptrvalue == AIP)) {    /* AIP was not previously in the AppearanceItemList */
    ValNodeAddPointer (&AP->AppearanceItemList, 0, AIP);
  }
}

NLM_EXTERN AppearancePtr DestroyAppearance (
  AppearancePtr AP,
  ViewerConfigsPtr VCP
)

{
  Uint1       i;
  ValNodePtr  VNP;

  if (AP == NULL || VCP == NULL) return NULL;
  if (AP->AppearanceItemList != NULL) {
    ValNodeFreeData (AP->AppearanceItemList);
  }
  for (VNP = VCP->AppearanceList; VNP != NULL; VNP = VNP->next) {
    if (VNP->data.ptrvalue == AP) {
      i = VNP->choice;
      VNP = ValNodeExtract (&VCP->AppearanceList, i);
      break;
    }
  }
  if (VNP != NULL) {
    MemFree (VNP);
    VNP = ValNodeExtract (&VCP->AppearanceNameList, i);
    MemFree (VNP);
  }
  MemFree (AP);
  return NULL;
}

static void getDim_do_not_render (
  RenderInputPtr RIP,
  Int4Ptr Start,
  Int4Ptr Stop,
  Uint2Ptr height,
  ViewerContextPtr vContext
)

{
  RelevantFeatureItemPtr RFIP;
  RFIP = RIP->RFIP;

  *Start = *Stop = RFIP->Left;
  *height = 1;
}

static void do_not_render (
  RenderInputPtr RIP,
  ViewerContextPtr vContext
)

{
  return;
}

static void getDim_render_with_line (
  RenderInputPtr RIP,
  Int4Ptr Start,
  Int4Ptr Stop,
  Uint2Ptr height,
  ViewerContextPtr vContext
)

{
  RelevantFeatureItemPtr RFIP;
  RFIP = RIP->RFIP;

  *height = 1;
  if (vContext->allFeatures) {
    *Start = RFIP->Left;
    *Stop = RFIP->Right;
  } else {
    *Start = MAX (RFIP->Left, vContext->from);
    *Stop = MAX (RFIP->Right, vContext->to);
  }
}

static void render_with_line (
  RenderInputPtr RIP,
  ViewerContextPtr vContext
)

{
  Uint4      StartY;
  PrimitivE  thisPrim;
  RelevantFeatureItemPtr RFIP;
  AppearanceItemPtr       AIP;
  Int4       start, stop;

  RFIP = RIP->RFIP;
  AIP = RIP->AIP;

  StartY = RIP->yStart - (RIP->featureOffset) - AIP->Height / 2;
  if (vContext->allFeatures) {
    start = RFIP->Left;
    stop = RFIP->Right;
  } else {
    start = MAX (RFIP->Left, vContext->from);
    stop = MAX (RFIP->Right, vContext->to);
  }

  thisPrim = AddLine (RIP->drawSeg, start, StartY, stop, StartY, 0, 0);
  SetPrimitiveIDs (thisPrim, RFIP->entityID, RFIP->itemID, RFIP->itemType, 0);
}

static void getDim_render_with_capped_line (
  RenderInputPtr RIP,
  Int4Ptr Start,
  Int4Ptr Stop,
  Uint2Ptr height,
  ViewerContextPtr vContext
)

{
  RelevantFeatureItemPtr RFIP;
  AppearanceItemPtr AIP;

  RFIP = RIP->RFIP;
  AIP = RIP->AIP;


  *height = AIP->Height;;
  if (vContext->allFeatures) {
    *Start = RFIP->Left;
    *Stop = RFIP->Right;
  } else {
    *Start = MAX (RFIP->Left, vContext->from);
    *Stop = MAX (RFIP->Right, vContext->to);
  }
}

static void render_with_capped_line (
  RenderInputPtr RIP,
  ViewerContextPtr vContext
)

{
  PrimitivE               thisPrim;
  AppearanceItemPtr       AIP;
  RelevantFeatureItemPtr  RFIP;
  Int4                    StartY;

  StartY = RIP->yStart - (RIP->featureOffset);
  render_with_line (RIP, vContext);
  RFIP = RIP->RFIP;
  AIP = RIP->AIP;
  thisPrim = AddLine (RIP->drawSeg, RFIP->Left, StartY, RFIP->Left, StartY - AIP->Height, 0, 0);
  SetPrimitiveIDs (thisPrim, RFIP->entityID, RFIP->itemID, RFIP->itemType, 0);
  thisPrim = AddLine (RIP->drawSeg, RFIP->Right, StartY, RFIP->Right, StartY - AIP->Height, 0, 0);
  SetPrimitiveIDs (thisPrim, RFIP->entityID, RFIP->itemID, RFIP->itemType, 0);
}

static void getDim_render_with_box (
  RenderInputPtr RIP,
  Int4Ptr Start,
  Int4Ptr Stop,
  Uint2Ptr height,
  ViewerContextPtr vContext
)

{
  RelevantFeatureItemPtr RFIP;
  AppearanceItemPtr AIP;

  RFIP = RIP->RFIP;
  AIP = RIP->AIP;
  *Start = RFIP->Left;
  *Stop = RFIP->Right;
  *height = AIP->Height;
}

static Boolean TestForSmearOverlap (
  Int4 PrevEnd,
  Int4 NewStart,
  ViewerContextPtr vContext
)

{
  Int4 PrevPixelStart;
  /* base-pair coordinates of the beginning of the last occupied pixel*/
  PrevPixelStart = PrevEnd - PrevEnd % vContext->viewScale;
  if (PrevPixelStart + vContext->viewScale > NewStart) return TRUE;
  return FALSE;
}

static Boolean TestForVisibleGap (
  Int4 x1,
  Int4 x2,
  ViewerContextPtr vContext
)

{
  Int4 xmin, xmax;

  if (x1 == x2) return FALSE;
  if (x1 > x2) {
    xmin = x2;
    xmax = x1;
  } else {
    xmin = x1;
    xmax = x2;
  }

  xmin -= xmin % vContext->viewScale; /* xmin is now the _beginning_ of the pixel it occupies in _base pair_ coordinates (I hope) */
  if (abs (xmax - xmin) < 2 * vContext->viewScale) return TRUE;
  return FALSE;
}

static Boolean TestForSmear (
  RelevantFeatureItemPtr RFIP1,
  RelevantFeatureItemPtr RFIP2,
  ViewerContextPtr vContext
)

{
  Uint4                  minSeperation;

  minSeperation = 5 * vContext->viewScale;  /* do not smear a feature more than 5 pixels wide */

  if (abs (RFIP1->Right - RFIP1->Left) >= minSeperation) return FALSE;
  if (abs (RFIP2->Right - RFIP2->Left) >= minSeperation) return FALSE;

  return (TestForVisibleGap (RFIP1->Right, RFIP2->Left, vContext)
          || TestForVisibleGap (RFIP1->Right, RFIP2->Right, vContext)
          || TestForVisibleGap (RFIP1->Left, RFIP2->Right, vContext)
          || TestForVisibleGap (RFIP1->Left, RFIP2->Left, vContext) );

}


static void render_with_box_master (
  RenderInputPtr RIP,
  Boolean fillBox,
  ViewerContextPtr vContext
)

{
  Uint4                   StartY;
  Uint2                   pieceIValStart;
  PrimitivE               thisPrim;
  Uint2                   i;
  Int4                    mid;
  AppearanceItemPtr       AIP;
  AppearancePtr           AP;
  RelevantFeatureItemPtr  RFIP;
  Boolean                 shade_p;
  Uint1                   arrow;

  RFIP = RIP->RFIP;
  AIP = RIP->AIP;
  AP = vContext->AppPtr;
  StartY = RIP->yStart - (RIP->featureOffset);
  if (RFIP->LeftEnd == EndClipped) {
    thisPrim = AddLine (RIP->drawSeg, vContext->from, StartY - AIP->Height / 2, RFIP->Left, StartY - AIP->Height / 2, 0, 0);
    SetPrimitiveIDs (thisPrim, RFIP->entityID, RFIP->itemID, RFIP->itemType, 0);
  }
  if (RFIP->RightEnd == EndClipped) {
    thisPrim = AddLine (RIP->drawSeg, vContext->to, StartY - AIP->Height / 2, RFIP->Right, StartY - AIP->Height / 2, 0, 0);
    SetPrimitiveIDs (thisPrim, RFIP->entityID, RFIP->itemID, RFIP->itemType, 0);
  }
  if (RFIP->numivals == 1) {
    shade_p = (RFIP->entityID == 0 && RFIP->itemID == 0) ? AP->ShadeSmears : FALSE;
    if (RFIP->entityID == 0 && RFIP->itemID == 0) { /* is this a multi-feature smear? */
      if (shade_p) {
        /*        AddAttribute (RIP->drawSeg, SHADING_ATT, 0, 0, MEDIUM_SHADING, 0, 0);*/
      }
      thisPrim = AddRectangle (RIP->drawSeg, RFIP->Left, StartY, RFIP->Right, StartY - AIP->Height, NO_ARROW, fillBox, 0);
      SetPrimitiveIDs (thisPrim, RFIP->entityID, RFIP->itemID, RFIP->itemType, 0);
      if (shade_p) {
        /*        AddAttribute (RIP->drawSeg, SHADING_ATT, 0, 0, NO_SHADING, 0, 0);*/
      }

    } else { /* nope */
      arrow = NO_ARROW;
      if (RFIP->plusstrand) {
        arrow = RIGHT_ARROW;
      } else {
        arrow = LEFT_ARROW;
      }
      if (! AIP->ShowArrow) {
        arrow = NO_ARROW;
      }
      if (ABS (RFIP->Right - RFIP->Left) / vContext->viewScale < AP->MinPixelsForArrow) {
        arrow = NO_ARROW;
      }
      if (vContext->viewScale > AP->MaxScaleForArrow) {
        arrow = NO_ARROW;
      }
      thisPrim = AddRectangle (RIP->drawSeg, RFIP->Left, StartY, RFIP->Right, StartY - AIP->Height, arrow, fillBox, 0);
      SetPrimitiveIDs (thisPrim, RFIP->entityID, RFIP->itemID, RFIP->itemType, 0);
    }
    return;
  } else {
    i = 0;
    while (i < RFIP->numivals) {
      /* collect a group of interval(s) which do not contain any pixels between them */
      pieceIValStart = i;
      /* this tests is the i-plus-1th feature is part of the smear, which goes from pieceIValStart to i, inclusive */
      while (i + 1 < RFIP->numivals &&
             TestForVisibleGap (vContext->viewScale + RFIP->ivals [2 * i + 1], RFIP->ivals [2 * i + 2], vContext)) {
        i++;
      }

      /* draw the segment and the gap -- drawing the gap first, so that it is overdrawn by the segment */
      if (i + 1 < RFIP->numivals) { /* a gap is present if there are more ivals to consider after i*/
        if (AIP->GapChoice == LineGap) {
          thisPrim = AddLine (RIP->drawSeg, RFIP->ivals [2 * i + 1], StartY - AIP->Height / 2, RFIP->ivals [2 * i + 2], StartY - AIP->Height / 2, 0, 0);
          SetPrimitiveIDs (thisPrim, RFIP->entityID, RFIP->itemID, RFIP->itemType, 0);
        } else if (AIP->GapChoice == AngleGap) {
          mid = (RFIP->ivals [2 * i + 2] + RFIP->ivals [2 * i + 1]) / 2;
          thisPrim = AddLine (RIP->drawSeg, RFIP->ivals [2 * i + 1], StartY - AIP->Height / 2, mid, StartY, 0, 0);
          SetPrimitiveIDs (thisPrim, RFIP->entityID, RFIP->itemID, RFIP->itemType, 0);
          thisPrim = AddLine (RIP->drawSeg, mid, StartY, RFIP->ivals [2 * i + 2], StartY - AIP->Height / 2, 0, 0);
          SetPrimitiveIDs (thisPrim, RFIP->entityID, RFIP->itemID, RFIP->itemType, 0);
        }
      }
      arrow = NO_ARROW;
      if (i == RFIP->numivals - 1) {
        if (RFIP->plusstrand) {
          arrow = RIGHT_ARROW;
        } else {
          arrow = LEFT_ARROW;
        }
      }
      if (! AIP->ShowArrow) {
        arrow = NO_ARROW;
      }
      if (ABS (RFIP->Right - RFIP->Left) / vContext->viewScale < AP->MinPixelsForArrow) {
        arrow = NO_ARROW;
      }
      if (vContext->viewScale > AP->MaxScaleForArrow) {
        arrow = NO_ARROW;
      }
      thisPrim = AddRectangle (RIP->drawSeg, RFIP->ivals [2 * pieceIValStart], StartY, RFIP->ivals [2 * i + 1], StartY - AIP->Height, arrow, fillBox, 0);
      SetPrimitiveIDs (thisPrim, RFIP->entityID, RFIP->itemID, RFIP->itemType, 0);
      i++;
    }
  }
}

static void render_with_box (
  RenderInputPtr   RIP,
  ViewerContextPtr vContext
)

{
  render_with_box_master (RIP, TRUE, vContext);

}

static void render_with_outline_box (
  RenderInputPtr   RIP,
  ViewerContextPtr vContext
)

{
  render_with_box_master (RIP, FALSE, vContext);
}

static const RenderClass RenderAlgorithmTable [] = {
  {do_not_render, getDim_do_not_render}, /* do_not_render */
  {render_with_line, getDim_render_with_line}, /* Render_Line */
  {render_with_capped_line, getDim_render_with_capped_line}, /* Render_CappedLine */
  {render_with_box, getDim_render_with_box},   /* Render_Box */
  {render_with_outline_box, getDim_render_with_box},   /* Render_OutlineBox */
  /* these do not exist right now */
  {render_with_line, getDim_render_with_line},
  {render_with_box, getDim_render_with_box},
  {render_with_line, getDim_render_with_line},
  {render_with_line, getDim_render_with_line},
  {render_with_line, getDim_render_with_line},
  {render_with_line, getDim_render_with_line},
  {render_with_line, getDim_render_with_line}
};

static void DrawFeatureAndLabel (
  RenderInputPtr RIP,
  ViewerContextPtr vContext
)

{
  AppearanceItemPtr       AIP;
  RelevantFeatureItemPtr  RFIP;
  FilterItemPtr           FIP;
  FilterPtr               FP;
  CharPtr                 label;
  Char                    tempStringBuffer [256];
  Char                    shortLabel [41];
  Uint1                   stringFlags;
  Uint4                   textWidthBP;
  Int4                    textStartX;
  Int4                    textStartY;
  Uint1                   labelAlign;
  PrimitivE               thisPrim;
  Boolean                 addType;
  Boolean                 addDesc;

  FIP = RIP->FIP;
  RFIP = RIP->RFIP;
  AIP = RIP->AIP;
  FP = vContext->FltPtr;
  addType = AIP->AddTypeToLabel;
  addDesc = AIP->AddDescToLabel;
  if (FIP->AddTypeToLabel != TristateUnset) {
    addType = BOOL_FROM_SET_TRISTATE (FIP->AddTypeToLabel);
  }
  if (FIP->AddDescToLabel != TristateUnset) {
    addDesc = BOOL_FROM_SET_TRISTATE (FIP->AddDescToLabel);
  }

  /*  RIP->drawSeg = CreateSegment (RIP->drawSeg, 0, 0);*/
  /* Place each feature in its own segment.  This is not necessary for simple features
     but perhaps the inefficiency is acceptable because it simplifies the algorithm. */
  (*RenderAlgorithmTable [AIP->RenderChoice].RenderFunc) (RIP, vContext);
  if (FIP->LabelLoc == LabelNone) return;
  if (vContext->viewScale > FP->MaxScaleWithLabels) return;
  stringFlags = 0;
  if (addDesc && RFIP->ContentLabel != NULL) {
    stringFlags |= 1;
  }
  if (addType) {
    stringFlags |= 2;
  }
  label = NULL;
  switch (stringFlags) {
    case 0:                    /* no label */
      return;
    case 1:                    /*comment but not type */
      label = RFIP->ContentLabel;
      break;
    case 2:                    /*only type */
      label = FindKeyFromFeatDefType (RFIP->featdeftype, FALSE);
      break;
    case 3:                    /*add both */
      if (StringCmp (FindKeyFromFeatDefType (RFIP->featdeftype, FALSE), RFIP->ContentLabel) != 0) {
        sprintf (tempStringBuffer, "%s: %s", FindKeyFromFeatDefType (RFIP->featdeftype, FALSE), RFIP->ContentLabel);
        label = tempStringBuffer;
      } else {
        label = RFIP->ContentLabel;
      }
      break;
  }
  if (StringHasNoText (label)) return;
  LabelCopy (shortLabel, label, sizeof (shortLabel));
  SelectFont (AIP->LabelFont);
  textWidthBP = StringWidth (shortLabel) * vContext->viewScale;
  switch (FIP->LabelLoc) {
    case LabelInside:
      if (textWidthBP + 2 * vContext->viewScale >= (RFIP->Right - RFIP->Left)) {
          /*!!! add string-chopper here, for now just don't render it */
        return;
      }
      textStartX = (RFIP->Left + RFIP->Right) / 2;
      /*      textStartY = RIP->yStart - RIP->rowHeight / 2; -- change "labelInside" to mean "above, but not wider than"*/
      textStartY = RIP->yStart;
      labelAlign = UPPER_CENTER;
      break;
    case LabelAbove:
      textStartX = (RFIP->Left + RFIP->Right) / 2;
      textStartY = RIP->yStart;
      labelAlign = UPPER_CENTER;
      break;
    case LabelBelow:
      textStartY = RIP->yStart - RIP->Height;
      textStartX = (RFIP->Left + RFIP->Right) / 2;
      labelAlign = LOWER_CENTER;
      break;
    case LabelLeft:
      textStartX = (RFIP->Left);
      textStartY = (RIP->yStart);
      labelAlign = LOWER_LEFT;
      break;
    case LabelRight:
      textStartX = (RFIP->Right);
      textStartY = (RIP->yStart);
      labelAlign = LOWER_RIGHT;
      break;
    default:
      return;
  }
  thisPrim = AddTextLabel (RIP->labelSeg, textStartX, textStartY, shortLabel, AIP->LabelFont, 1, labelAlign, 0);
  SetPrimitiveIDs (thisPrim, RFIP->entityID, RFIP->itemID, RFIP->itemType, 0);
}

static RelevantFeatureItemPtr BuildClippedRFIP (
  RelevantFeatureItemPtr inputRFIP,
  FilterProcessStatePtr FPSP,
  ViewerContextPtr vContext
)

{
  RelevantFeatureItemPtr newRFIP;
  ValNodePtr VNP = NULL;
  Uint2 i, newnumivals;
  Int4  from, to;
  Boolean useThis = TRUE, useLast;
  Boolean clippedLeft, clippedRight;
  Boolean lastClippedLeft, lastClippedRight;

  if ((newRFIP = MemNew (sizeof (RelevantFeatureItem))) == NULL) return NULL;
  if ((VNP = MemNew (sizeof (ValNode))) == NULL) {
    MemFree (newRFIP);
    return NULL;
  }
  VNP->data.ptrvalue = newRFIP;
  VNP->next = FPSP->needFreeList;
  FPSP->needFreeList = VNP;
  MemCopy (newRFIP, inputRFIP, sizeof (RelevantFeatureItem));
  if ((inputRFIP->Left <= vContext->from && inputRFIP->Right <= vContext->from)
      || (inputRFIP->Left >= vContext->to && inputRFIP->Right >= vContext->to)) {
    return NULL; /* entire feature removed by clipping */
  }
  newRFIP->LeftEnd = (inputRFIP->Left >= vContext->from) ? inputRFIP->LeftEnd : EndPartial;
  newRFIP->RightEnd  = (inputRFIP->Right <= vContext->to) ?  inputRFIP->RightEnd  : EndPartial;
  if (inputRFIP->numivals == 1) {
    newRFIP->Left  = MAX (vContext->from, MIN (vContext->to, inputRFIP->Left ));
    newRFIP->Right = MAX (vContext->from, MIN (vContext->to, inputRFIP->Right));
    newRFIP->LeftEnd = (inputRFIP->Left >= vContext->from) ? inputRFIP->LeftEnd : EndPartial;
    newRFIP->RightEnd  = (inputRFIP->Right <= vContext->to) ?  inputRFIP->RightEnd  : EndPartial;
  } else {
    newRFIP->Left = vContext->to;
    newRFIP->Right = vContext->from;
    newnumivals = 0;
    for (i = 0; i < inputRFIP->numivals; i++) {
      from = inputRFIP->ivals [2 * i];
      to = inputRFIP->ivals [2 * i + 1];
      if (from <= vContext->to && from >= vContext->from) {
        newnumivals ++;
      } else if (to <=  vContext->to && to >= vContext->from) {
        newnumivals ++;
      }
    }
    newRFIP->numivals = newnumivals;
    if ((newRFIP->ivals = MemNew (2 * newnumivals * sizeof (Int4))) == NULL) return NULL;
    if ((VNP = MemNew (sizeof (ValNode))) == NULL) {
      MemFree (newRFIP->ivals);
      return NULL;
    }
    VNP->data.ptrvalue = newRFIP->ivals;
    VNP->next = FPSP->needFreeList;
    FPSP->needFreeList = VNP;
    newnumivals = 0;
    for (i = 0; i < inputRFIP->numivals; i++) {
      from = inputRFIP->ivals [2 * i];
      to = inputRFIP->ivals [2 * i + 1];
      if (from <= vContext->to && from >= vContext->from) {
        useThis = TRUE;
        clippedLeft = clippedRight = FALSE;
      } else if (to <=  vContext->to && to >= vContext->from) {
        useThis = TRUE;
        clippedLeft = clippedRight = FALSE;
      } else {
        useThis = FALSE;
        clippedLeft = ((from + to) < (vContext->from + vContext->to));
        clippedRight = !clippedLeft;
      }
      if (i == 0) {
        useLast = useThis;
        lastClippedLeft = clippedLeft;
        lastClippedRight = clippedRight;
      }
      if (lastClippedLeft && useThis) {
        newRFIP->LeftEnd  = EndClipped;
      } else if (lastClippedRight && useThis) {
        newRFIP->RightEnd = EndClipped;
      } else if (useLast && clippedLeft) {
        newRFIP->LeftEnd  = EndClipped;
      } else if (useLast && clippedRight) {
        newRFIP->RightEnd = EndClipped;
      }
      if (useThis) {
        from = MAX (vContext->from, MIN (vContext->to, from));
        to   = MAX (vContext->from, MIN (vContext->to, to  ));
        newRFIP->Left = MIN (newRFIP->Left, from);
        newRFIP->Left = MIN (newRFIP->Left, to);
        newRFIP->Right = MAX (newRFIP->Right, to);
        newRFIP->Right = MAX (newRFIP->Right, from);
        newRFIP->ivals [newnumivals++] = from;
        newRFIP->ivals [newnumivals++] = to;
      }
      useLast = useThis;
    }
  }
  if (newRFIP->Left > newRFIP->Right) return NULL;
  return newRFIP;
}


static void GetFeatureAndDecorationDimensions (
  RenderInputPtr RIP,
  ViewerContextPtr vContext
)

{
  AppearanceItemPtr       AIP;
  RelevantFeatureItemPtr  RFIP;
  FilterItemPtr           FIP;
  FilterPtr               FP;
  Int4                    Start, Stop;
  Uint2                   Height;
  Int2                    lineHeight;
  Int4                    textStartX;
  Uint4                   textWidthBP;
  CharPtr                 label = NULL;
  Char                    tempStringBuffer [256];
  Uint1                   stringFlags;
  Uint2                   featureOffset = 0;
  Boolean                 addType;
  Boolean                 addDesc;

  RFIP = RIP->RFIP;
  AIP = RIP->AIP;
  FIP = RIP->FIP;
  FP = vContext->FltPtr;
  addType = AIP->AddTypeToLabel;
  addDesc = AIP->AddDescToLabel;
  if (FIP->AddTypeToLabel != TristateUnset) {
    addType = BOOL_FROM_SET_TRISTATE (FIP->AddTypeToLabel);
  }
  if (FIP->AddDescToLabel != TristateUnset) {
    addDesc = BOOL_FROM_SET_TRISTATE (FIP->AddDescToLabel);
  }

  (*RenderAlgorithmTable [AIP->RenderChoice].GetDimensions) (RIP, &Start, &Stop, &Height, vContext);
  RIP->Height = Height;
  RIP->decorationHeight = Height;
  RIP->decorationLeft = RFIP->Left;
  RIP->decorationRight = RFIP->Right;
  RIP->featureOffset = 0;

  if (FIP->LabelLoc != LabelNone && vContext->viewScale <= FP->MaxScaleWithLabels) {
    stringFlags = 0;
    if (addDesc && RFIP->ContentLabel != NULL) {
      stringFlags |= 1;
    }
    if (addType) {
      stringFlags |= 2;
    }

    switch (stringFlags) {
      case 0:                  /* no label */
        break;
      case 1:                  /*comment but not type */
        label = RFIP->ContentLabel;
        break;
      case 2:                  /*only type */
        label = FindKeyFromFeatDefType (RFIP->featdeftype, FALSE);
        break;
      case 3:                  /*add both */
        if (StringCmp (FindKeyFromFeatDefType (RFIP->featdeftype, FALSE), RFIP->ContentLabel) != 0) {
          sprintf (tempStringBuffer, "%s: %s", FindKeyFromFeatDefType (RFIP->featdeftype, FALSE), RFIP->ContentLabel);
          label = tempStringBuffer;
        } else {
          label = RFIP->ContentLabel;
        }
        break;
      default:
        return;
    }
  }
  if (! StringHasNoText (label)) {
    SelectFont (AIP->LabelFont);
    textWidthBP = StringWidth (label) * vContext->viewScale;
    lineHeight = LineHeight ();
    switch (FIP->LabelLoc) {
      case LabelInside:
        Height += lineHeight + 3;
        featureOffset = lineHeight + 1;
        break;
      case LabelAbove:
        textStartX = (Start + Stop) / 2;
        Start = MIN (Start, (signed)(textStartX - textWidthBP / 2));
        Stop = MAX (Stop, textStartX + textWidthBP / 2);
        featureOffset = lineHeight + 1;
        Height += lineHeight + 3;
        break;
      case LabelBelow:
        textStartX = (Start + Stop) / 2;
        Start = MIN (Start, (signed)(textStartX - textWidthBP / 2));
        Stop = MAX (Stop, textStartX + textWidthBP / 2);
        Height += lineHeight + 3;
        break;
      case LabelLeft:
        Start -= textWidthBP;
        Height = MAX (Height + lineHeight, Height);
        break;
      case LabelRight:
        Stop = RFIP->Right + textWidthBP;
        Height = MAX (Height + lineHeight, Height);
        break;
      default:
        return;
    }
  }
  RIP->decorationLeft = Start;
  RIP->decorationRight = Stop;
  RIP->decorationHeight = Height;
}


static Boolean BuildRenderInputFromRFIP (
  RenderInputPtr target,
  RelevantFeatureItemPtr RFIP,
  FilterProcessStatePtr FPSP
)

{
  AppearancePtr    AppPtr;
  FilterItemPtr    currentFIP;
  ViewerContextPtr vContext;
  RelevantFeatureItemPtr newRFIP;

  vContext = FPSP->vContext;

  if (target == NULL || RFIP == NULL) return FALSE;
  if (!vContext->allFeatures &&
      (
       (RFIP->Left < vContext->from && RFIP->Right < vContext->from)
       || (RFIP->Left > vContext->to && RFIP->Right > vContext->to)
       )) {
    return FALSE; /* this feature is outside the clipping seqloc */
  }
  if (! vContext->allFeatures && (
      (RFIP->Right >= vContext->from || RFIP->Left <= vContext->to)
      || (RFIP->Right < vContext->from && RFIP->Left > vContext->to)
      )) {
    newRFIP = BuildClippedRFIP (RFIP, FPSP, vContext);
    if (newRFIP == NULL) return FALSE;
    RFIP = newRFIP;
  }
  target->RFIP = RFIP;
  target->labelSeg = FPSP->labelSegs [RFIP->featdeftype];
  target->drawSeg = FPSP->drawSegs [RFIP->featdeftype];
  AppPtr = vContext->AppPtr;
  currentFIP = FPSP->currentFIP;
  target->AIP = AppPtr->FeaturesAppearanceItem [RFIP->featdeftype];
  target->FIP = FPSP->currentFIP;
  GetFeatureAndDecorationDimensions (target, vContext);
  target->rowHeight = MAX (target->Height, target->decorationHeight) + currentFIP->IntraRowPaddingPixels;
  return TRUE;
}

/* todo: perhaps switch to using Explore functions */
static Boolean GetAndCountFeatures (
  ViewerContextPtr vContext
)

{
  SeqFeatPtr              sfp;
  SeqMgrFeatContext       fContext;
  RelevantFeatureItemPtr  rFeats;
  ValNodePtr              sapList = NULL, VNP, VNPtail;
  Uint2                   i = 0;
  Int4                    swap;

  if (vContext == NULL) return FALSE;
  vContext->sapCount = 0;
  vContext->featureCount = 0;
  vContext->featVNP = NULL;
  vContext->sapList = NULL;

  rFeats = MemNew (RELEVANT_FEATS_PER_CHUNK * sizeof (RelevantFeatureItem));
  if (rFeats == NULL) return FALSE;
  ValNodeAddPointer (&vContext->featVNP, 0, rFeats);
  VNPtail = vContext->featVNP;
  sfp = SeqMgrGetNextFeature (vContext->BSP, NULL, 0, 0, &fContext);
  while (sfp != NULL) {
    vContext->featureCount++;

    rFeats [i].Left = fContext.left;
    rFeats [i].Right = fContext.right;
    rFeats [i].LeftEnd = fContext.partialL ? EndPartial : EndAbsolute;
    rFeats [i].RightEnd  = fContext.partialR ? EndPartial : EndAbsolute;
    rFeats [i].ContentLabel = fContext.label;
    rFeats [i].featdeftype = fContext.featdeftype;
    rFeats [i].entityID = fContext.entityID;
    rFeats [i].itemID = fContext.itemID;
    rFeats [i].itemType = OBJ_SEQFEAT;
    rFeats [i].numivals = fContext.numivals;
    rFeats [i].ivals = fContext.ivals;
    rFeats [i].plusstrand = (fContext.strand != Seq_strand_minus);
    if (rFeats [i].Right < rFeats [i].Left) {
      /* protection against (feature indexing vs. trans-spliced features) */
      swap = rFeats [i].Right;
      rFeats [i].Right = rFeats [i].Left;
      rFeats [i].Left = swap;
    }
    /* look at fContext.sap -- if unique, add to sapList */
    for (VNP = sapList; VNP != NULL; VNP = VNP->next) {
      if (VNP->data.ptrvalue == fContext.sap) break;
    }
    if (VNP == NULL) {          /* fContext.sap was not found */
      vContext->sapCount++;
      ValNodeAddPointer (&sapList, 0, fContext.sap);
    }
    i++;
    if (i >= RELEVANT_FEATS_PER_CHUNK) {
      i = 0;
      rFeats = MemNew (RELEVANT_FEATS_PER_CHUNK * sizeof (RelevantFeatureItem));
      VNPtail = ValNodeNew (VNPtail);
      if (rFeats == NULL || VNPtail == NULL) {
        ValNodeFreeData (vContext->featVNP);
        MemFree (rFeats);
        return FALSE;
      }
      VNPtail->data.ptrvalue = rFeats;
    }
    sfp = SeqMgrGetNextFeature (vContext->BSP, sfp, 0, 0, &fContext);
  }
  if (vContext->sapCount > 0) {
    vContext->sapList = MemNew (vContext->sapCount * sizeof (SeqAnnotPtr));
    if (vContext->sapList == NULL) {
      MemFree (rFeats);
      ValNodeFree (sapList);
      return FALSE;
    }
    for (i = 0, VNP = sapList; VNP != NULL && i < vContext->sapCount; VNP = VNP->next, i++) {
      vContext->sapList[i] = VNP->data.ptrvalue;
    }
    ValNodeFree (sapList);
  }
  if (vContext->featureCount == 0) {
    MemFree (rFeats);
    MemFree (vContext->featVNP);
    vContext->featVNP = NULL;
  }
  return TRUE;
}

NLM_EXTERN RelevantFeaturesPtr CollectFeatures (
  BioseqPtr bsp
)

{
  RelevantFeaturesPtr  RFP;
  ViewerContext        VC;

  RFP = MemNew (sizeof (RelevantFeatures));
  if (RFP == NULL) return NULL;
  VC.BSP = bsp;
  if (! GetAndCountFeatures (&VC)) return NULL;
  RFP->featureCount = VC.featureCount;
  RFP->featVNP = VC.featVNP;
  RFP->sapCount = VC.sapCount;
  RFP->sapList = VC.sapList;
  return RFP;
}

static Boolean EnsureFeatureHasSegment (
  FilterProcessStatePtr FPSP,
  Uint1 featdeftype,
  SegmenT parentSegment
)

{
  AppearancePtr     AppPtr;
  AppearanceItemPtr AIP;
  ViewerContextPtr vContext;

  vContext = FPSP->vContext;


  if (parentSegment == NULL) {
    parentSegment = vContext->drawOnMe;
  }

  AppPtr = vContext->AppPtr;
  if (FPSP->drawSegs [featdeftype] == NULL) {
    FPSP->drawSegs [featdeftype] = CreateSegment (parentSegment, 0, 0);
    FPSP->labelSegs [featdeftype] = CreateSegment (parentSegment, 0, 0);
    /* cleaup needed if program is supposed to recover from this !!! */
    if (FPSP->drawSegs [featdeftype] == NULL || FPSP->labelSegs [featdeftype] == NULL) return FALSE;
    AIP = AppPtr->FeaturesAppearanceItem [featdeftype];
    AddAttribute (FPSP->drawSegs [featdeftype],
                  COLOR_ATT | SHADING_ATT | STYLE_ATT | WIDTH_ATT,
                  AIP->Color, AIP->VibLinestyle, AIP->VibShading, 1, 0);
    AddAttribute (FPSP->labelSegs [featdeftype], COLOR_ATT, AIP->LabelColor, 0, 0, 0, 0);
  }
  return TRUE;
}


static RelevantFeatureItemPtr GetNextRFIPinAlignmentFilter (
  FilterProcessStatePtr FPSP
)

{
  AlignmentFilterStatePtr alignSP;
  SeqAlignPtr             SAlnP;
  Int4                    alignRow, start, stop;
  SeqIdPtr                SID;
  BioseqPtr               BSP;
  RelevantFeatureItemPtr  RFIP;
  ViewerContextPtr        vContext;
  FilterItemPtr           currentFIP;

  if (FPSP == NULL) return NULL;
  vContext = FPSP->vContext;
  currentFIP = FPSP->currentFIP;
  if (vContext == NULL || currentFIP == NULL) return NULL;

  alignSP = &FPSP->state.align;

  if (alignSP->SAPcurrent == NULL) return NULL;      
  alignSP->SAPcurrent = alignSP->SAPcurrent->next;
  SAlnP = alignSP->SAPcurrent;
  RFIP = MemNew (sizeof (RelevantFeatureItem));
  if (RFIP == NULL || (ValNodeAddPointer (&FPSP->needFreeList, 0, RFIP)) == NULL) {
    MemFree (RFIP);
    return NULL;
  }
  if (!AlnMgr2IndexSingleChildSeqAlign (SAlnP)) return NULL;
  BSP = vContext->BSP;
  SID = BSP->id;
  alignRow = AlnMgr2GetFirstNForSip (SAlnP, SID);
  if (alignRow == -1) return NULL;
  AlnMgr2GetNthSeqRangeInSA (SAlnP, alignRow, &start, &stop);
  if (start < 0 || stop < 0) return NULL;
  RFIP->Left = MIN (start, stop);
  RFIP->Right = MAX (start, stop);
  RFIP->plusstrand = (start < stop);
  RFIP->numivals = 1;
  return RFIP;
}

static RelevantFeatureItemPtr GetNextRFIPinFeatureFilter (
  FilterProcessStatePtr FPSP
)

{
  ValNodePtr              currentRFIPblockVNP;
  FeatureFilterStatePtr   featSP;
  RelevantFeatureItemPtr  RFIP;
  ViewerContextPtr        vContext;
  FilterItemPtr           currentFIP;

  if (FPSP == NULL) return NULL;
  vContext = FPSP->vContext;
  currentFIP = FPSP->currentFIP;
  if (vContext == NULL || currentFIP == NULL) return NULL;
  featSP = &FPSP->state.feat;
  for (; featSP->featureBlockOffset + featSP->indexInBlock < vContext->featureCount; featSP->indexInBlock++) {
    if (featSP->indexInBlock >= RELEVANT_FEATS_PER_CHUNK) {
      featSP->indexInBlock = 0;
      featSP->featureBlockOffset += RELEVANT_FEATS_PER_CHUNK;
      featSP->currentRFIPblockVNP = featSP->currentRFIPblockVNP->next;
      if (featSP->currentRFIPblockVNP == NULL) return NULL;
    }
    currentRFIPblockVNP = featSP->currentRFIPblockVNP;
    RFIP = (RelevantFeatureItemPtr) (currentRFIPblockVNP->data.ptrvalue) + featSP->indexInBlock;
    if (! vContext->allFeatures
        && (RFIP->Right < vContext->from || RFIP->Left > vContext->to)) {
      continue;
    }
    if (FPSP->featuresProcessed [featSP->featureBlockOffset + featSP->indexInBlock]
        || (! currentFIP->IncludeFeature [RFIP->featdeftype])) continue;
    if (currentFIP->MatchStrand != BothStrands) {
      if (currentFIP->MatchStrand == MinusStrand && RFIP->plusstrand) continue;
      if (currentFIP->MatchStrand == PlusStrand && !RFIP->plusstrand) continue;
    }
    FPSP->featuresProcessed [featSP->featureBlockOffset + featSP->indexInBlock] = TRUE;
    break;
  }
  if (featSP->featureBlockOffset + featSP->indexInBlock >= vContext->featureCount) return NULL;
  return RFIP;
}



static RelevantFeatureItemPtr GetNextRFIPinFilterItem (
  FilterProcessStatePtr FPSP
)

{
  AlignmentFilterStatePtr alignSP;
  RelevantFeatureItemPtr  RFIP;
  ViewerContextPtr        vContext;
  FilterItemPtr           currentFIP;

  /* called as an iterator by the rendering functions -- builds & returns the next feature (in an RFIP), or returns NULL if no more left in this Filter */
  if (FPSP == NULL) return NULL;
  vContext = FPSP->vContext;
  currentFIP = FPSP->currentFIP;
  if (vContext == NULL || currentFIP == NULL) return NULL;
  switch (currentFIP->Type) {
  case InvalidFilter:
  case GraphFilter:
  case BioseqFilter:
    return NULL;
  case AlignmentFilter:
    alignSP = &FPSP->state.align;
    do {
      RFIP = GetNextRFIPinAlignmentFilter (FPSP);
      if (RFIP != NULL) return RFIP;
    } while (alignSP->SAPcurrent != NULL);
    /* note: if control reaches here, then RFIP == NULL */
    break;
  case FeatureFilter:
    RFIP = GetNextRFIPinFeatureFilter (FPSP);
    break;
  }
  if (RFIP != NULL) {
    FPSP->featuresProcessedCount++;
  }
  return RFIP;
}

static Boolean AddFeatureToRow (
  InternalRowPtr row,
  RelevantFeatureItemPtr RFIP,
  Boolean SkipSmearTest,
  FilterProcessStatePtr FPSP
)

{
  ValNodePtr             VNP;
  RelevantFeatureItemPtr newRFIP; /* for representing a multi-feature smear */
  RelevantFeatureItemPtr oldRFIP;
  ViewerContextPtr       vContext;

  vContext = FPSP->vContext;

  if (row == NULL || RFIP == NULL) return FALSE;

  if (!SkipSmearTest && row->rowFeatureCount &&
      TestForSmear (RFIP, (RelevantFeatureItemPtr) row->feats->data.ptrvalue, vContext)) {
    /* if the last feature was not a smear-in-progress, allocate newRFIP, else re-use the current one */
    oldRFIP = (RelevantFeatureItemPtr) row->feats->data.ptrvalue;
    if (oldRFIP->entityID == 0 && oldRFIP->itemID == 0) {
      /* oldRFIP is already a smear-in-prorgess, just extend it */
      newRFIP = oldRFIP;
    } else {
      /* need to create a new smear */
      newRFIP = MemNew (sizeof (RelevantFeatureItem));
      VNP = MemNew (sizeof (ValNode));
      if (newRFIP == NULL || VNP == NULL) return FALSE;
      VNP->data.ptrvalue = newRFIP;
      VNP->next = FPSP->needFreeList;
      FPSP->needFreeList = VNP;
      row->feats->data.ptrvalue = newRFIP;
      newRFIP->featdeftype = oldRFIP->featdeftype;
      newRFIP->numivals = 1;
    }

    newRFIP->Left = MIN (oldRFIP->Left, RFIP->Left);
    newRFIP->Right = MAX (oldRFIP->Right, RFIP->Right);
  } else {

    VNP = MemNew (sizeof (ValNode));
    if (VNP == NULL) return FALSE;
    VNP->data.ptrvalue = RFIP;
    VNP->next = row->feats;
    row->feats = VNP;
    row->rowFeatureCount++;
  }
  return TRUE;
}

static InternalRowPtr AddARow (
  InternalRowPtr sourceRow
)

{
  InternalRowPtr  IRP;

  if (sourceRow == NULL) return NULL;
  IRP = MemNew (sizeof (InternalRow));
  if (IRP == NULL) return NULL;
  sourceRow->next = IRP;
  return IRP;
}

static Uint2 SimpleDiagonalLayout (
  InternalRowPtr firstRow,
  FilterProcessStatePtr FPSP
)

{
  Uint2                   rows = 0;
  InternalRowPtr          thisRow;
  RelevantFeatureItemPtr  RFIP;

  thisRow = firstRow;
  if (thisRow == NULL) return 0;
  thisRow->rowFeatureCount = 0;
  while ((RFIP = GetNextRFIPinFilterItem (FPSP)) != NULL) {
    AddFeatureToRow (thisRow, RFIP, TRUE, FPSP);
    thisRow = AddARow (thisRow);
    rows++;
  }
  return rows;
}

/* this was copied from seqmgr.c (6.181, 27-feb-2002) */
static Boolean CheckInternalExonBoundaries (Int2 numivalsCDS, Int4Ptr ivalsCDS, Int2 numivalsMRNA, Int4Ptr ivalsMRNA)

{
  Int2  i;
  Int2  j;

  if (numivalsCDS > numivalsMRNA) return FALSE;
  if (ivalsCDS == NULL || ivalsMRNA == NULL) return TRUE;

  /* scan first exon-intron boundary against candidate start positions */

  for (i = 0; i <= numivalsMRNA - numivalsCDS; i++) {
    if (ivalsCDS [1] == ivalsMRNA [2 * i + 1]) break;
  }
  if (i > numivalsMRNA - numivalsCDS) return FALSE;

  /* Addition by Eric: the first interval in the CDS must not be larger than the corresponding interval in the mRNA */
  if (ABS (ivalsCDS [0] - ivalsCDS [1]) > ABS (ivalsMRNA [2 * i] - ivalsMRNA [2 * i + 1])) return FALSE;

  /* scan subsequent exon-intron and intron-exon boundaries */

  for (j = 2; j <= 2 * numivalsCDS - 2; j++) {
    if (ivalsCDS [j] != ivalsMRNA [2 * i + j]) return FALSE;
  }

  /* Addition by Eric: the last interval in the CDS must not be larger than the corresponding interval in the mRNA */
  if (ABS (ivalsCDS [j - 1] - ivalsCDS [j]) > ABS (ivalsMRNA [2 * i + j - 1] - ivalsMRNA [2 * i + j])) return FALSE;

  return TRUE;
}

static Boolean mRNAmatchesCDS (
  RelevantFeatureItemPtr mRNA,
  RelevantFeatureItemPtr CDS
)

{
  /* the mRNA must be the larger feature */
  if (CDS->Left < mRNA->Left || CDS->Right > mRNA->Right) return FALSE;
  /* check strands */
  if (CDS->plusstrand != mRNA->plusstrand) return FALSE;
  /* trivial case */
  if (CDS->numivals == 1 && mRNA->numivals == 1) return TRUE;
  /* . . . and the intervals must line up */
  return CheckInternalExonBoundaries (CDS->numivals, CDS->ivals, mRNA->numivals, mRNA->ivals);
}

typedef struct rFIPentry {
  RelevantFeatureItemPtr RFIP;
  Int4                   decorationLeft;
  Int4                   decorationRight;
  Int4                   Right;
  Int4                   Left;
  Boolean                used;
} RFIPentry, PNTR RFIPentryPtr;

typedef struct rFIPgroup {
  Uint4       memberCount;
  Int4        decorationLeft;
  Int4        decorationRight;
  Int4        Left;
  Int4        Right;
  ValNodePtr  members; /* data.ptrvalue = RFIPentryPtr */
} RFIPgroup, PNTR RFIPgroupPtr;


static Uint2 GeneProductsLayoutInternal (
  InternalRowPtr firstRow,
  FilterProcessStatePtr FPSP,
  Boolean MultiRender
)

{
  Uint2                   rows = 1, i;
  RelevantFeatureItemPtr  RFIP, RFIPcds, RFIPmrna, tRFIP;
  InternalRowPtr          thisRow, thisRow2, lastRow;
  ValNodePtr              fifoHead = NULL, fifoTail = NULL;
  ValNodePtr              bumpedListHead, bumpedListTail;
  ValNodePtr              potCDSvnp, potRNAvnp;
  ValNodePtr              GroupsHead = NULL, GroupsTailVNP, MemberTailVNP;
  ValNodePtr              groupVNP, featVNP;
  RFIPentryPtr            tRFIPentry;
  RFIPentryPtr            bumpedListEntry;
  RFIPentryPtr            potCDS, potRNA, GroupsTailMemberTail;
  Uint2                   bumpedCount;
  Boolean                 foundRows, doneCollecting = FALSE;
  Int4                    newLeft, featureStart, rowMaxRight;
  RenderInput             dummyRI;
  RFIPgroupPtr            currentGroup;

  if (firstRow == NULL) return 0;

  lastRow = firstRow;
  firstRow->layoutData.intvalue = -2000000000; /* the first feature will _always_ fit in the first row */

  while (1) {
    if (doneCollecting) {
      break;
    }
    RFIP = GetNextRFIPinFilterItem (FPSP);
    if (RFIP == NULL) {
      doneCollecting = TRUE;
      bumpedListHead = fifoHead; /* no more incoming features, bump all features in the queue*/
      bumpedListTail = fifoTail;
      fifoHead = NULL;
    } else {
      if (! BuildRenderInputFromRFIP (&dummyRI, RFIP, FPSP)) {
        continue; /* either we're out of memory or this feature doesn't overlap the clipping SeqLoc */
      }

      if ((tRFIPentry = MemNew (sizeof (RFIPentry))) == NULL) goto bail_out;
      if (fifoHead == NULL) {
        fifoTail = ValNodeAddPointer (&fifoHead, 0, tRFIPentry);
      } else {
        fifoTail = ValNodeAddPointer (&fifoTail, 0, tRFIPentry);
      }
      if (fifoTail == NULL) goto bail_out;
      tRFIPentry->RFIP = RFIP;
      tRFIPentry->Left = RFIP->Left;
      tRFIPentry->Right = RFIP->Right;
      tRFIPentry->decorationLeft = dummyRI.decorationLeft;
      tRFIPentry->decorationRight = dummyRI.decorationRight;

    /*
      now find all features which can not overlap the next feature (in the sequence; decoration overlap doesn't matter)
    */
      newLeft = RFIP->Left;

      bumpedCount = 0;
      bumpedListHead = fifoHead;
      for (featVNP = fifoHead;
           featVNP != NULL;
           featVNP = featVNP->next) {
        bumpedListEntry = featVNP->data.ptrvalue;
        RFIP = bumpedListEntry->RFIP;
        if (RFIP->Right > newLeft) break;
        bumpedListTail = featVNP;
        fifoHead = featVNP->next;
        bumpedCount++;
      }
      if (bumpedCount == 0) continue;
    }
    if (bumpedListTail != NULL) {
      bumpedListTail->next = NULL;
    }
    GroupsHead = GroupsTailVNP = NULL;

    for (potCDSvnp = bumpedListHead; potCDSvnp != NULL; potCDSvnp = potCDSvnp->next) {
      /* this pass is only searching for CDS features (and correspondins mRNAs) */
      potCDS = potCDSvnp->data.ptrvalue;
      RFIPcds = potCDS->RFIP;
      if (RFIPcds->featdeftype != FEATDEF_CDS) continue;
      potCDS->used = TRUE;

      if ((currentGroup = MemNew (sizeof (RFIPgroup))) == NULL) goto bail_out;
      if (GroupsHead == NULL) {
        GroupsTailVNP = ValNodeAddPointer (&GroupsHead, 0, currentGroup);
      } else {
        GroupsTailVNP = ValNodeAddPointer (&GroupsTailVNP, 0, currentGroup);
      }
      if (GroupsTailVNP == NULL) goto bail_out;

      /* do not add the CDS to currentGroup->members yet, b/c we want mRNA features to appear first
         (but remember to add it after!) */

      currentGroup->memberCount = 1;
      currentGroup->members = NULL;

      currentGroup->Left = potCDS->Left;
      currentGroup->Right = potCDS->Right;
      currentGroup->decorationLeft = potCDS->decorationLeft;
      currentGroup->decorationRight = potCDS->decorationRight;

      for (potRNAvnp = bumpedListHead; potRNAvnp != NULL; potRNAvnp = potRNAvnp->next) {
        potRNA = potRNAvnp->data.ptrvalue;
        RFIPmrna = potRNA->RFIP;
        if (RFIPmrna->featdeftype != FEATDEF_mRNA) continue;
        if (!MultiRender && potRNA->used) continue;

        if (mRNAmatchesCDS (RFIPmrna, RFIPcds)) {
          potRNA->used = TRUE;

          if ((tRFIPentry = MemNew (sizeof (RFIPentry))) == NULL) goto bail_out;
          if (ValNodeAddPointer (&currentGroup->members, 0, tRFIPentry) == NULL) goto bail_out;

          MemCopy (tRFIPentry, potRNA, sizeof (RFIPentry));

          currentGroup->memberCount++;
          currentGroup->Left = MIN (currentGroup->Left, RFIPmrna->Left);
          currentGroup->Right = MAX (currentGroup->Right, RFIPmrna->Right);
          currentGroup->decorationLeft = MIN (currentGroup->decorationLeft, potRNA->decorationLeft);
          currentGroup->decorationRight = MAX (currentGroup->decorationRight, potRNA->decorationRight);
        }
      }

      if ((GroupsTailMemberTail = MemNew (sizeof (RFIPentry))) == NULL) goto bail_out;
      if (ValNodeAddPointer (&currentGroup->members, 0, GroupsTailMemberTail) == NULL) goto bail_out;
      MemCopy (GroupsTailMemberTail, potCDS, sizeof (RFIPentry));

    }
    /*
      append all non-matched elements to the Groups list
    */
    for (featVNP = bumpedListHead; featVNP != NULL; featVNP = featVNP->next) {
      tRFIPentry = featVNP->data.ptrvalue;
      /* skip feature if it's been matched already */
      if (tRFIPentry->used) {
        continue;
      } else {
        /* add another singleton entry to Groups */
        tRFIP = tRFIPentry->RFIP;

        if ((currentGroup = MemNew (sizeof (RFIPgroup))) == NULL) goto bail_out;
        if (GroupsHead == NULL) {
          GroupsTailVNP = ValNodeAddPointer (&GroupsHead, 0, currentGroup);
        } else {
          GroupsTailVNP = ValNodeAddPointer (&GroupsTailVNP, 0 ,currentGroup);
        }
        if (GroupsTailVNP == NULL) goto bail_out;

        if ((GroupsTailMemberTail = MemNew (sizeof (RFIPentry))) == NULL) goto bail_out;
        if ((MemberTailVNP = ValNodeAddPointer (&currentGroup->members, 0, GroupsTailMemberTail)) == NULL) goto bail_out;

        currentGroup->memberCount = 1;
        MemCopy (GroupsTailMemberTail, tRFIPentry, sizeof (RFIPentry));

        currentGroup->Left = tRFIP->Left;
        currentGroup->Right = tRFIP->Right;
        currentGroup->decorationLeft = tRFIPentry->decorationLeft;
        currentGroup->decorationRight = tRFIPentry->decorationRight;

      }
    }

    /*
      now, assign each element in Groups to a row.  algorithm is the same as in BubbleUpLayout,
      except that instead of looking for the first matching row, need to find (RFIPgroup->members)
      consecutive free rows.
    */
    for (groupVNP = GroupsHead; groupVNP != NULL; groupVNP = groupVNP->next) {
      currentGroup = groupVNP->data.ptrvalue;

      featureStart = currentGroup->decorationLeft;
      for (thisRow = firstRow; thisRow != NULL; thisRow = thisRow->next) {
        foundRows = TRUE;
        for (i = 0, thisRow2 = thisRow, featVNP = currentGroup->members;
             featVNP != NULL && i < currentGroup->memberCount;
             thisRow2 = thisRow2->next, featVNP = featVNP->next, i++) {
          tRFIPentry = featVNP->data.ptrvalue;
/*
  note: if thisRow2 ends up being NULL then it meas thisRow was a good starting point and we just need to add some rows
*/
          if (thisRow2 == NULL) {
            foundRows = TRUE;
            break;
          }
          rowMaxRight = thisRow2->layoutData.intvalue;
          if (tRFIPentry->decorationLeft < rowMaxRight) {
            foundRows = FALSE;
            if (thisRow2->next == NULL) {
              rows++;
              if ((lastRow = AddARow (lastRow)) == NULL) goto bail_out;
              lastRow->layoutData.intvalue = -2000000000;
            }
            break;
          }
        }
        if (foundRows) {
          for (i = 0,
                 thisRow2 = thisRow,
                 featVNP = currentGroup->members;
               tRFIPentry != NULL
                 && i < currentGroup->memberCount;
               featVNP = featVNP->next,
                 thisRow2 = thisRow2->next,
                 i++) {
            tRFIPentry = featVNP->data.ptrvalue;
            RFIP = tRFIPentry->RFIP;
            if (thisRow2 == NULL) {
              rows++;
              if ((lastRow = AddARow (lastRow)) == NULL) goto bail_out;
              thisRow2 = lastRow;
            }
            AddFeatureToRow (thisRow2, RFIP, FALSE, FPSP);
            rowMaxRight = tRFIPentry->decorationRight;
            thisRow2->layoutData.intvalue = rowMaxRight;
          }
          break;
        }
      }
    }
    bumpedListHead = ValNodeFreeData (bumpedListHead);
    for (groupVNP = GroupsHead; groupVNP != NULL; groupVNP = groupVNP->next) {
      currentGroup = groupVNP->data.ptrvalue;
      ValNodeFreeData (currentGroup->members);
    }
    GroupsHead = ValNodeFreeData (GroupsHead);
  }
bail_out:
  ValNodeFreeData (bumpedListHead);
  ValNodeFreeData (fifoHead);
  for (groupVNP = GroupsHead; groupVNP != NULL; groupVNP = groupVNP->next) {
    currentGroup = groupVNP->data.ptrvalue;
    ValNodeFreeData (currentGroup->members);
  }
  GroupsHead = ValNodeFreeData (GroupsHead);
  return rows;
}

static Uint2 GeneProductsLayout (
  InternalRowPtr firstRow,
  FilterProcessStatePtr FPSP
)

{
  return GeneProductsLayoutInternal (firstRow, FPSP, FALSE);
}

static Uint2 GeneProductsLayoutX (
  InternalRowPtr firstRow,
  FilterProcessStatePtr FPSP
)

{
  return GeneProductsLayoutInternal (firstRow, FPSP, TRUE);
}

static Uint2 BubbleUpLayout (
  InternalRowPtr firstRow,
  FilterProcessStatePtr FPSP
)

{
  Uint2                   rows = 1;
  Int4                    featureStart, rowMaxRight;
  RelevantFeatureItemPtr  RFIP;
  InternalRowPtr          thisRow, lastRow;
  RenderInput             dummyRI;
  /*
     This uses InternalRow.layoutData as an Int4, which stores the (x) offset of the last used pixel.
     so a feature starting at (internalRow.layoutData + 1) or greater can be placed in the same row
     without a collision.
   */

  /* add the 1st feature to the 1st row */
  do {
    RFIP = GetNextRFIPinFilterItem (FPSP);
    if (RFIP == NULL) return 0;
  } while (! BuildRenderInputFromRFIP (&dummyRI, RFIP, FPSP));

  AddFeatureToRow (firstRow, RFIP, FALSE, FPSP);
  lastRow = firstRow;
  firstRow->layoutData.intvalue = dummyRI.decorationRight;

  while ((RFIP = GetNextRFIPinFilterItem (FPSP)) != NULL) {
    if (! BuildRenderInputFromRFIP (&dummyRI, RFIP, FPSP)) {
      continue;
    }
    featureStart = dummyRI.decorationLeft;
    for (thisRow = firstRow; thisRow != NULL; thisRow = thisRow->next) {
      rowMaxRight = thisRow->layoutData.intvalue;
      if (featureStart >= rowMaxRight) {
        AddFeatureToRow (thisRow, RFIP, FALSE, FPSP);
        rowMaxRight = MAX (rowMaxRight, dummyRI.decorationRight);
        thisRow->layoutData.intvalue = rowMaxRight;
        RFIP = NULL;
        break;
      }
    }
    if (RFIP != NULL) {
      rows++;
      lastRow = AddARow (lastRow);
      AddFeatureToRow (lastRow, RFIP, FALSE, FPSP);
      rowMaxRight = dummyRI.decorationRight;
      lastRow->layoutData.intvalue = rowMaxRight;
    }
  }
  return rows;
}

static Uint2 SingleRowLayout (
  InternalRowPtr firstRow,
  FilterProcessStatePtr FPSP
)

{
  RelevantFeatureItemPtr RFIP, lastRFIP;
  ValNodePtr       freeVNP;
  Uint4            smearStart, smearStop;
  Uint1            smearFeatdeftype;
  ViewerContextPtr vContext;
  Boolean          adjacentNoGap;
  Boolean          tooNarrow;
  Boolean          smearing = FALSE, justSmeared = FALSE;

  vContext = FPSP->vContext;

  lastRFIP = GetNextRFIPinFilterItem (FPSP);
  if (lastRFIP == NULL)
    return 0;
  smearFeatdeftype = lastRFIP->featdeftype;

  while ((RFIP = GetNextRFIPinFilterItem (FPSP)) != NULL) {
    adjacentNoGap = TestForSmearOverlap (lastRFIP->Right + 2 * vContext->viewScale, RFIP->Left, vContext);
    tooNarrow = TestForSmearOverlap (RFIP->Left + 4 * vContext->viewScale, RFIP->Right, vContext);
    if (adjacentNoGap && tooNarrow) {
      if (smearing == FALSE) {
        tooNarrow = TestForSmearOverlap (lastRFIP->Left + 4 * vContext->viewScale, lastRFIP->Right, vContext);
        if (tooNarrow) {
          smearStart = lastRFIP->Left;
        } else {
          smearStart = RFIP->Left;
        }
        smearing = TRUE;
      }
      justSmeared = TRUE;
      if (smearFeatdeftype != RFIP->featdeftype) {
        smearFeatdeftype = FEATDEF_ANY;
      } else if (smearFeatdeftype == 0) {
        smearFeatdeftype = RFIP->featdeftype;
      }
      smearStop = RFIP->Right;
      continue;
    }
    /* VNP is non-NULL during a smear; this tests for regular features*/
    if (! justSmeared) {
      AddFeatureToRow (firstRow, lastRFIP, TRUE, FPSP);
      lastRFIP = RFIP;
    } else { /* we just finished a smear group */
      lastRFIP = MemNew (sizeof (RelevantFeatureItem));
      if (lastRFIP == NULL) return 1;
      lastRFIP->Left = smearStart;
      lastRFIP->Right = smearStop;
      lastRFIP->plusstrand = TRUE;
      lastRFIP->LeftEnd = lastRFIP->RightEnd = EndAbsolute;
      lastRFIP->featdeftype = smearFeatdeftype;
      lastRFIP->numivals = 1;
      ValNodeAddPointer (&FPSP->needFreeList, 0, lastRFIP);
      smearing = FALSE;
      AddFeatureToRow (firstRow, lastRFIP, TRUE, FPSP);
      justSmeared = FALSE;
      lastRFIP = RFIP;
    }
  }
  if (! justSmeared) { /* the last feature was not in a smear group */
    AddFeatureToRow (firstRow, lastRFIP, TRUE, FPSP);
  } else {
    RFIP = MemNew (sizeof (RelevantFeatureItem));
    if (RFIP == NULL) return 1;
    RFIP->Left = smearStart;
    RFIP->Right = smearStop;
    RFIP->plusstrand = TRUE;
    RFIP->LeftEnd = RFIP->RightEnd = EndAbsolute;
    RFIP->featdeftype = smearFeatdeftype;
    RFIP->numivals = 1;

    if ((freeVNP = MemNew (sizeof (ValNode))) == NULL) {
      MemFree (RFIP);
      return 1;
    }
    freeVNP->data.ptrvalue = RFIP;
    freeVNP->next = FPSP->needFreeList;
    FPSP->needFreeList = freeVNP;

    AddFeatureToRow (firstRow, RFIP, TRUE, FPSP);
  }
  return 1;
}


static const LayoutFunction LayoutAlgorithmTable [] = {
  BubbleUpLayout,             /* placeholder for Layout_Inherit */
  SimpleDiagonalLayout,       /* Layout_Diagonal */
  /*  SimpleDiagonalLayout,*/ /* Layout_DiagonalSawtooth (to be implemented) */
  SingleRowLayout,            /* Layout_FeatTypePerLine (same as single-row, but features are grouped by type before processing) */
  BubbleUpLayout,             /* Layout_FeatTypePerLineGroup (same as bubble-up, but features are grouped by type before processing) */
  /*SingleRowLayout,*/        /* Layout_AllInOneLine (which isn't working currently, and may be less useful than expected*/
  BubbleUpLayout,             /* Layout_PackUpward */
  GeneProductsLayout,         /* Layout_GroupCorrespondingFeats (?working?) */
  GeneProductsLayoutX         /* Layout_GroupCorrespondingFeatsRepeat (?working?) */
};

/* non-leaf segments contain SEGMENT_TREE_BASE other segments */
#define SEGMENT_TREE_BASE 16
#define SEGMENT_TREE_BASE2 (SEGMENT_TREE_BASE * SEGMENT_TREE_BASE)
#define SEGMENT_TREE_BASE3 (SEGMENT_TREE_BASE * SEGMENT_TREE_BASE * SEGMENT_TREE_BASE)

static Uint2 ProcessRows (
  LayoutAlgorithm layoutC,
  FilterProcessStatePtr FPSP,
  ViewerContextPtr vContext
)

{
  Uint2                   I, J, K, lastI, lastJ, lastK; /* for keeping track of which segment in the tree we're in*/
  Uint1                   featdeftype;
  Int4                    featMidPoint;
  Uint2                   row, rowCount, feat, rowHeight, totalHeight;
  InternalRowPtr          firstRow, thisRow;
  ValNodePtr              VNP;
  RelevantFeatureItemPtr  RFIP;
  RenderInput             RI;           /* dummy used while finding height of each row */
  SegmenT                 SegmentTreeTop;
  SegmenT                 SegmentTreeMid;
  SegmenT                 SegmentTreeBot;
  Boolean                 SegmentChanged = TRUE;
  Boolean                 emptyRow = TRUE;
  Boolean                 allEmpty = TRUE;
  Boolean                 NeedToDrawSegRect;
  FilterItemPtr           FIP;

  firstRow = MemNew (sizeof (InternalRow));
  if (firstRow == NULL) return 0;
  firstRow->feats = NULL;
  VNP = FPSP->currentFilterVNP;
  if (VNP == NULL) return 0;
  FIP = (FilterItemPtr) VNP->data.ptrvalue;
  NeedToDrawSegRect = FIP->DrawItemRect;

  rowCount = (*LayoutAlgorithmTable [layoutC]) (firstRow, FPSP);

  thisRow = firstRow;
  totalHeight = 0;

  for (row = 0; row < rowCount; row++) {
    if (thisRow == NULL) continue;
    /*First iterate through features to find the row's height */
    VNP = thisRow->feats;
    rowHeight = 0;
    if (thisRow->rowFeatureCount > 0) {
      allEmpty = FALSE;
    }
    for (feat = 0; feat < thisRow->rowFeatureCount; feat++) {
      RFIP = VNP->data.ptrvalue;
      if (VNP == NULL) break;
      VNP = VNP->next;
      if (RFIP == NULL) return 0;
      if (! BuildRenderInputFromRFIP (&RI, RFIP, FPSP)) {
        continue;
      }
      rowHeight = MAX (rowHeight, RI.decorationHeight);
    }
    if (!allEmpty && NeedToDrawSegRect) {
      AddAttribute (vContext->drawOnMe, COLOR_ATT, FIP->GroupBoxColor, 0, 0, 0, 0);
      AddSegRect (vContext->drawOnMe, FALSE, 0);
      NeedToDrawSegRect = FALSE;
    }

    /*Repeat, but actually draw them this time */
    VNP = thisRow->feats;
    RI.rowHeight = rowHeight;
    RI.yStart = FPSP->ceiling - totalHeight;
    if (thisRow->rowFeatureCount > 0) {
      emptyRow = FALSE;
    } else {
      emptyRow = TRUE;
    }
    I = J = K = 0;
    SegmentChanged = TRUE;
    for (feat = 0; feat < thisRow->rowFeatureCount; feat++) {
      if (VNP == NULL) break;
      if ((RFIP = VNP->data.ptrvalue) == NULL) return 0;
      VNP = VNP->next;
      featdeftype = RFIP->featdeftype;
      featMidPoint = (RFIP->Left + RFIP->Right) / 2;
      lastI = I;
      lastJ = J;
      lastK = K;
      I = (SEGMENT_TREE_BASE * featMidPoint) / vContext->seqLength;
      J = ((SEGMENT_TREE_BASE2 * featMidPoint) / vContext->seqLength) % SEGMENT_TREE_BASE;
      K = ((SEGMENT_TREE_BASE3 * featMidPoint) / vContext->seqLength) % SEGMENT_TREE_BASE;

      if (I != lastI || SegmentChanged) {
        SegmentTreeTop  = CreateSegment (vContext->drawOnMe, 0, 0);
        SegmentChanged = TRUE;
      }
      if (J != lastJ || SegmentChanged) {
        SegmentTreeMid = CreateSegment (SegmentTreeTop , 0, 0);
        SegmentChanged = TRUE;
      }
      if (K != lastK || SegmentChanged) {
        SegmentTreeBot = CreateSegment (SegmentTreeMid, 0, 0);
        SegmentChanged = TRUE;
      }
      if (SegmentChanged) {
        MemSet (FPSP->drawSegs, 0, sizeof (FPSP->drawSegs));
        MemSet (FPSP->labelSegs, 0, sizeof (FPSP->labelSegs));
        SegmentChanged = FALSE;
      }
      if (! EnsureFeatureHasSegment (FPSP, featdeftype, SegmentTreeBot)) return 0;
      if (! BuildRenderInputFromRFIP (&RI, RFIP, FPSP)) {
        continue;
      }

      DrawFeatureAndLabel (&RI, vContext);
    }
    if (! emptyRow) {
      totalHeight += rowHeight + FPSP->currentFIP->IntraRowPaddingPixels;
    }
    if (thisRow->next == NULL) return totalHeight;
    thisRow = thisRow->next;
  }
  if (allEmpty) {
    return 0;
  }
  /*  totalHeight += rowHeight;*/
  totalHeight += FPSP->currentFIP->GroupPadding;
  return totalHeight;
}

typedef struct scalesntdata {
  SegmenT                   seg;
  PrimitivE                 snt;
  BoxInfo                   box;
  Int4                      length;
  Int4                      labelOffset;
  BioseqAppearanceItemPtr   bioseqAIP;
  Boolean                   disableLastLabel;
  Boolean                   disableFirstLabel;
  Int4                      offset;
} ScaleSntData, PNTR ScaleSntPtr;

static void ScaleSntDrawProc (
  BigScalar calldata,
  PrimDrawContext pdc
)

{
  Int4                     curScale;
  DrawInfoPtr              dInfoPtr;
  Int4                     i, j;
  RecT                     r;
  PntInfo                  pnt;
  PoinT                    pt;
  ScaleSntPtr              ssp;
  BoxInfo                  tmpBox;
  Uint1                    tickCount;
  Char                     buffer[16];
  BioseqAppearanceItemPtr  bioseqAIP;
  Uint1                    scaleHeight;
  Int4                     from, to;

  ssp = (ScaleSntPtr) calldata;
  if (ssp == NULL) return;

  dInfoPtr = (DrawInfoPtr) pdc;
  tmpBox = ssp->box;
  bioseqAIP = ssp->bioseqAIP;
  scaleHeight = bioseqAIP->scaleHeight;

  from = tmpBox.left;
  to   = tmpBox.right;

  curScale = dInfoPtr->scale.scaleX;
  r.left   = (Int2) ((dInfoPtr->scale.offsetX + tmpBox.left)  / curScale);
  r.right  = (Int2) ((dInfoPtr->scale.offsetX + tmpBox.right)  / curScale);
  SelectFont (bioseqAIP->scaleFont);
  SelectColor (bioseqAIP->scaleColor[0], bioseqAIP->scaleColor[1], bioseqAIP->scaleColor[2]);

  /* if the right-most edge is visible, draw the final tick mark and a right-justified total length label */
  if (dInfoPtr->scale.worldWindow.right >= to - 100 * curScale) {
    pnt.x = to;
    pnt.y = tmpBox.top;
    MapWorldPointToPixel (&pt, &pnt, &dInfoPtr->scale);
    MoveTo (pt.x, pt.y);
    pnt.y = tmpBox.top - scaleHeight;
    MapWorldPointToPixel (&pt, &pnt, &dInfoPtr->scale);
    LineTo (pt.x, pt.y);
    sprintf (buffer, "%ld", (long)to + ssp->labelOffset + 1);
    pt.x += 3 - StringWidth (buffer);
    pt.y += 1 + Ascent ();
    PaintStringEx (buffer, pt.x, pt.y);
  }

  /* if the left-most edge is visible, draw the first label, left-justified */
  if (dInfoPtr->scale.worldWindow.left - 400 * curScale <= from) {
    pnt.x = from;
    pnt.y = tmpBox.top - scaleHeight;
    sprintf (buffer, "%ld", (long)((from == 0) ? 1 : from)  + ssp->labelOffset);
    MapWorldPointToPixel (&pt, &pnt, &dInfoPtr->scale);
    pt.y += 1 + Ascent ();
    PaintStringEx (buffer, pt.x, pt.y);
  }

  if (tmpBox.left > dInfoPtr->scale.worldWindow.right ||
      tmpBox.right < dInfoPtr->scale.worldWindow.left ||
      tmpBox.top < dInfoPtr->scale.worldWindow.bottom ||
      tmpBox.bottom > dInfoPtr->scale.worldWindow.top)
    return;

  if (dInfoPtr->checked == FALSE ) {
    if (tmpBox.left < dInfoPtr->scale.worldWindow16.left)
      tmpBox.left = dInfoPtr->scale.worldWindow16.left;
    if (tmpBox.right > dInfoPtr->scale.worldWindow16.right)
      tmpBox.right = dInfoPtr->scale.worldWindow16.right;
    if (tmpBox.top > dInfoPtr->scale.worldWindow16.top)
      tmpBox.top = dInfoPtr->scale.worldWindow16.top;
    if (tmpBox.bottom < dInfoPtr->scale.worldWindow16.bottom)
      tmpBox.bottom = dInfoPtr->scale.worldWindow16.bottom;
  }

  i = MAX (dInfoPtr->scale.worldWindow.left, tmpBox.left);
  /*  i = i + 10 * curScale - (i % (10 * curScale)) - 1;*/
  while (i % (10 * curScale) != 0) {
    i++; /* !!! do this the right way */
  }
  for (tickCount = (i / (10 * curScale)) % 10;
       i < MIN (dInfoPtr->scale.worldWindow.right, to);
       i += curScale * 10) {
    pnt.x = i;
    pnt.y = tmpBox.top;
    MapWorldPointToPixel (&pt, &pnt, &dInfoPtr->scale);
    MoveTo (pt.x, pt.y);

    if (tickCount == 0 || tickCount == 5) {
      pnt.y = tmpBox.top - scaleHeight; /* draw full-height tick */
      MapWorldPointToPixel (&pt, &pnt, &dInfoPtr->scale);
      LineTo (pt.x, pt.y);
      if (tickCount == 0
          /* if i + 100curScale > the right boundary, this is the last label */
          && (!((i + 100 * curScale >= to) && ssp->disableLastLabel))
          && (!((i - 100 * curScale <= from) && ssp->disableFirstLabel))
          && (i != from)) {
        sprintf (buffer, "%ld", (long)((i == 0) ? 1 : i) + ssp->labelOffset); /* put the origin at 1 instead of 0 */
        pt.x -= StringWidth (buffer) / 2; /* center text around the tick mark*/
        pt.y += 1 + Ascent ();
        PaintStringEx (buffer, pt.x, pt.y);
      }
      /* disable the right-end one in the case of an overlap */
    } else {
      pnt.y = tmpBox.top - scaleHeight / 2;
      MapWorldPointToPixel (&pt, &pnt, &dInfoPtr->scale);
      LineTo (pt.x, pt.y);
    }
    tickCount = (tickCount + 1) % 10;
  }

  /* also, we might need to redraw the labels just to the left or to the right of the update region. */
  i = MAX (dInfoPtr->scale.worldWindow.left, tmpBox.left);
  j = i; /* j = leftmost visible pixel in world coordinates */
  while (i % (100 * curScale) != 0) {
    i--;
  }
  pnt.x = i;
  sprintf (buffer, "%ld", (long)((i == 0) ? 1 : i) + ssp->labelOffset);
  i += curScale * ((StringWidth (buffer) / 2) + 1);
  if (i >= j
      && i > from
      && (!((pnt.x + 100 * curScale >= to) && ssp->disableLastLabel))
      && (!((pnt.x - 100 * curScale <= from) && ssp->disableFirstLabel))
      && pnt.x > from) {
    pnt.y = tmpBox.top - scaleHeight;
    MapWorldPointToPixel (&pt, &pnt, &dInfoPtr->scale);
    pt.x -= StringWidth (buffer) / 2; /* center text around the tick mark*/
    pt.y += 1 + Ascent ();
    PaintStringEx (buffer, pt.x, pt.y);
  }
  /* now check the right side*/

  i = MIN (dInfoPtr->scale.worldWindow.right, to);
  j = i; /* j = rightmost visible pixel in world coordinates */
  while (i % (100 * curScale) != 0) {
    i++;
  }
  pnt.x = i;
  sprintf (buffer, "%ld", (long)i + ssp->labelOffset);
  i -= curScale * ((StringWidth (buffer) / 2) + 1);
  if (i <= j
      && (!((pnt.x + 100 * curScale >= to) && ssp->disableLastLabel))
      && (!((pnt.x - 100 * curScale <= from) && ssp->disableFirstLabel))
      && pnt.x < to) {
    pnt.y = tmpBox.top - scaleHeight;
    MapWorldPointToPixel (&pt, &pnt, &dInfoPtr->scale);
    pt.x -= StringWidth (buffer) / 2; /* center text around the tick mark*/
    pt.y += 1 + Ascent ();
    PaintStringEx (buffer, pt.x, pt.y);
  }
}

static void ScaleSntCleanupProc (
  BigScalar calldata
)

{
  ScaleSntPtr  ssp;

  ssp = (ScaleSntPtr) calldata;
  MemFree (ssp);
}

static Boolean LIBCALLBACK Asn2gphSegmentExploreProc (
  SeqLocPtr slp,
  SeqMgrSegmentContextPtr context
)

{
  /* context->userdata a pointer to vContext */
  RelevantFeatureItemPtr RFIP;
  SeqIdPtr           sip, newid = NULL;
  Int4               gi;
  Char               labelBuf [100];
  Int4               left, right;
  ViewerContextPtr   vContext;
  AppearancePtr      AP;
  BioseqAppearanceItemPtr BioAIP;

  if ((RFIP = MemNew (sizeof (RelevantFeatureItem))) == NULL) return FALSE;
  vContext = context->userdata;
  ValNodeAddPointer (&vContext->BSPsegmentVNP, 0, RFIP);
  sip = SeqLocId (slp);
  if (sip != NULL && sip->choice == SEQID_GI) {
    gi = sip->data.intvalue;
    newid = SeqIdStripLocus (GetSeqIdForGI (gi));
    sip = newid;
    if (newid != NULL) {
      ValNodeAddInt (&newid, SEQID_GI, gi);
    }
  }
  
  AP = vContext->AppPtr;
  BioAIP = AP->bioseqAIP;

  left = context->from;
  right = context->to - left + context->cumOffset;
  left = context->cumOffset;
  if (sip != NULL) {
    if (BioAIP->format != PRINTID_FASTA_LONG) { 
      sip = SeqIdFindWorst (sip);
    }
    SeqIdWrite (sip, labelBuf, BioAIP->format, sizeof (labelBuf) - 1);
    RFIP->ContentLabel = StringSave (labelBuf); /* record this in a list of things to be freed so that we dont leak memory!!!!! */
  }
  SeqIdSetFree (newid);
  RFIP->Left = MIN (left, right);
  RFIP->Right = MAX (left, right);
  RFIP->entityID = context->entityID;
  RFIP->itemID = context->itemID;
  RFIP->itemType = OBJ_BIOSEQ_SEG;
  RFIP->numivals = 1;
  RFIP->featdeftype = 0;
  RFIP->plusstrand = (context->strand == Seq_strand_plus);
  return TRUE;
}

static Int4 DrawBioseqSegments (
  Int4 initialYOffset,
  FilterItemPtr FIP,
  ViewerContextPtr vContext
)

{
  ValNodePtr   listHead = NULL, vnpSegment = NULL, vnpRow = NULL;
  ValNodePtr   firstRow = NULL;
  Int4         segmentCount, rowCount;
  Uint2         rowHeight, row;
  AppearanceItem segmentAI;
  AppearanceItemPtr oldAIPzero;
  AppearancePtr AP;
  BioseqAppearanceItemPtr BioAIP;
  RelevantFeatureItemPtr RFIP;
  FilterProcessState FPS;
  Uint2         oldMaxScale;

  vContext->BSPsegmentVNP = NULL;
  SeqMgrExploreSegments (vContext->BSP, vContext, Asn2gphSegmentExploreProc);
  listHead = vContext->BSPsegmentVNP;
  segmentCount = 0;
  for (vnpSegment = listHead; vnpSegment != NULL; vnpSegment = vnpSegment->next) {
    segmentCount++;
  }

  if (segmentCount <= 1) {
    ValNodeFreeData (listHead);
    return 0;
  }

  oldMaxScale = vContext->FltPtr->MaxScaleWithLabels;
  vContext->FltPtr->MaxScaleWithLabels = 10000;
  AP = vContext->AppPtr;
  BioAIP = AP->bioseqAIP;
  MemCpy (segmentAI.LabelColor, BioAIP->labelColor, sizeof (segmentAI.LabelColor));
  MemCpy (segmentAI.Color, BioAIP->bioseqColor, sizeof (segmentAI.Color));
  segmentAI.RenderChoice = Render_Box;
  segmentAI.Height = BioAIP->height;
  segmentAI.AddDescToLabel = TRUE;
  segmentAI.AddTypeToLabel = FALSE;
  segmentAI.LabelFont = BioAIP->labelFont;
  segmentAI.LabelLoc = FIP->LabelLoc;
  FIP->LabelLoc = LabelAbove;

  segmentAI.ShowArrow = FALSE;
  segmentAI.VibLinestyle = SOLID_LINE;
  segmentAI.VibShading = SOLID_SHADING;
  oldAIPzero = vContext->AppPtr->FeaturesAppearanceItem[0];
  vContext->AppPtr->FeaturesAppearanceItem[0] = &segmentAI;

  MemSet (&FPS, 0, sizeof (FPS));
  FPS.currentFIP = FIP;
  FPS.ceiling = initialYOffset;
  FPS.vContext = vContext;

  FPS.renderParm.AIP = &segmentAI;
  FPS.renderParm.FIP = FIP;

  SelectFont (segmentAI.LabelFont);
  rowHeight = LineHeight () + segmentAI.Height + 3 + FIP->IntraRowPaddingPixels;
  rowCount = 0;

  for (vnpSegment = listHead; vnpSegment != NULL; vnpSegment = vnpSegment->next) {

    RFIP = vnpSegment->data.ptrvalue;
    EnsureFeatureHasSegment (&FPS, 0, vContext->drawOnMe);
    if (!BuildRenderInputFromRFIP (&FPS.renderParm, RFIP, &FPS)) {
      continue;
    }
    row = 0;
    for (vnpRow = firstRow; vnpRow != NULL; vnpRow = vnpRow->next) {
      if (vnpRow->data.intvalue < FPS.renderParm.decorationLeft) break;
      row++;
    }

    if (row >= rowCount) {
      rowCount++; 
      vnpRow = ValNodeAddInt (&firstRow, 0, 0);
    }
    vnpRow->data.intvalue = FPS.renderParm.decorationRight;
    FPS.renderParm.rowHeight += 2;
    FPS.renderParm.yStart = initialYOffset - (row * rowHeight);
    DrawFeatureAndLabel (&FPS.renderParm, vContext);
  }
  if (FPS.needFreeList != NULL) {
    ValNodeFreeData (FPS.needFreeList);
  }
  ValNodeFreeData (listHead);
  vContext->BSPsegmentVNP = NULL;
  vContext->AppPtr->FeaturesAppearanceItem[0] = oldAIPzero;
  vContext->FltPtr->MaxScaleWithLabels = oldMaxScale;
  return (rowCount) * rowHeight + 7 + FIP->IntraRowPaddingPixels + FIP->GroupPadding;;
}

static Int4 DrawBioseq (
  Int4 initialYOffset,
  FilterItemPtr FIP,
  ViewerContextPtr vContext
)

{
  PrimitivE               thisPrim;
  Char                    labelbuf[128];
  BioseqPtr               bsp;
  ScaleSntPtr             ssp;
  Uint2                   labelLineHeight, lastLabelWidth, firstLabelWidth;
  SegmenT                 seg;
  SeqIdPtr                sip;
  BioseqAppearanceItemPtr bioseqAIP;
  Uint1                   scaleHeight;
  Uint1                   height;
  Int4                    yOffset;
  Uint4                   i, j;
  Boolean                 drawScale;
  AppearancePtr           AP;

  AP = vContext->AppPtr;
  bioseqAIP = AP->bioseqAIP;
  yOffset = initialYOffset;
  bsp = vContext->BSP;

  drawScale = bioseqAIP->drawScale;
  if (FIP->DrawScale != TristateUnset) {
    drawScale = BOOL_FROM_SET_TRISTATE(FIP->DrawScale);
  }

  scaleHeight = bioseqAIP->scaleHeight;
  if (bioseqAIP->format != PRINTID_FASTA_LONG) {
    sip = SeqIdFindWorst (bsp->id);
  } else {
    sip = bsp->id;
  }
  SeqIdWrite (sip, labelbuf, bioseqAIP->format, sizeof (labelbuf) - 1);
  seg = CreateSegment (vContext->topLevelSeg, 0, 0);

  if (bioseqAIP->labelLoc == LabelOnTop) {
    AddAttribute (seg, COLOR_ATT, bioseqAIP->labelColor, 0, 0, 1, 0);
    thisPrim = AddTextLabel (seg, vContext->from, yOffset, labelbuf, bioseqAIP->labelFont, 0, LOWER_RIGHT, 0);
    SetPrimitiveIDs (thisPrim, bsp->idx.entityID, bsp->idx.itemID, OBJ_BIOSEQ, 0);
    SelectFont (bioseqAIP->labelFont);
    height = LineHeight() + 1;
    yOffset = initialYOffset - height;
  } else {
    height = 0;
  }

  AddAttribute (seg, COLOR_ATT, bioseqAIP->bioseqColor, 0, 0, 1, 0);
  thisPrim = AddRectangle (seg, vContext->from, yOffset, vContext->to, yOffset - bioseqAIP->height, NO_ARROW, TRUE, 0);
  SetPrimitiveIDs (thisPrim, bsp->idx.entityID, bsp->idx.itemID, OBJ_BIOSEQ, 0);
  AddAttribute (seg, COLOR_ATT, bioseqAIP->labelColor, 0, 0, 1, 0);

  height += bioseqAIP->height + 2;
  yOffset = initialYOffset - height;

  if (bioseqAIP->labelLoc == LabelOnSide) {
    AddAttribute (seg, COLOR_ATT, bioseqAIP->labelColor, 0, 0, 1, 0);
    thisPrim = AddTextLabel (seg, vContext->from, yOffset, labelbuf, bioseqAIP->labelFont, 1, UPPER_LEFT, 0);
    SetPrimitiveIDs (thisPrim, bsp->idx.entityID, bsp->idx.itemID, OBJ_BIOSEQ, 0);
  } else if (bioseqAIP->labelLoc == LabelOnBottom) {
    AddAttribute (seg, COLOR_ATT, bioseqAIP->labelColor, 0, 0, 1, 0);
    thisPrim = AddTextLabel (seg, 0, yOffset, labelbuf, bioseqAIP->labelFont, 1, LOWER_RIGHT, 0);
    SetPrimitiveIDs (thisPrim, bsp->idx.entityID, bsp->idx.itemID, OBJ_BIOSEQ, 0);
    SelectFont (bioseqAIP->labelFont);
    height += LineHeight() + 2;
    yOffset = initialYOffset - height;
  }

  if (drawScale) {

    ssp = (ScaleSntPtr) MemNew (sizeof (ScaleSntData));
    if (ssp == NULL) return height;
    SelectFont (bioseqAIP->scaleFont);
    labelLineHeight = LineHeight();
    ssp->seg = seg;
    ssp->labelOffset = 0;
    ssp->length = vContext->to - vContext->from;
    ssp->box.left = vContext->from;
    ssp->box.right = vContext->to;
    ssp->box.top = yOffset;
    ssp->box.bottom = yOffset - 2 - scaleHeight - labelLineHeight;
    ssp->bioseqAIP = bioseqAIP;
    
    i = vContext->to - vContext->to % (100 * vContext->viewScale); /* the last labelled tick-mark */
    
    sprintf (labelbuf, "%ld", (long)i + ssp->labelOffset);
    j = StringWidth (labelbuf); /* j is the width (pixels) of the last tick mark's label */
    
    sprintf (labelbuf, "%ld", (long)vContext->to + ssp->labelOffset);
    lastLabelWidth = (StringWidth (labelbuf) + 2) * vContext->viewScale;
    
    if (vContext->from > 0) {
      sprintf (labelbuf, "%ld", (long)vContext->from + ssp->labelOffset);
    } else {
      sprintf (labelbuf, "1");
    }
    firstLabelWidth = StringWidth (labelbuf) * vContext->viewScale;
    
    i += vContext->viewScale * j / 2; /* right-most edge of the centered label */
    j = MAX (j, 10); /* always leave at least 10 pixels between the last 2 labels */
    if (i + j * vContext->viewScale >= vContext->to - lastLabelWidth) {
      ssp->disableLastLabel = TRUE;
    } else {
      ssp->disableLastLabel = FALSE;
    }
    
    i = vContext->from + (100 * vContext->viewScale) - vContext->from % (100 * vContext->viewScale); /* the second label (first is at the origin) */
    sprintf (labelbuf, "%ld", (long)vContext->from + ssp->labelOffset);
    if (i - ((StringWidth (labelbuf) / 2 + 5) * vContext->viewScale) <= vContext->from + firstLabelWidth) {
      ssp->disableFirstLabel = TRUE;
    } else {
      ssp->disableFirstLabel = FALSE;
    }
    
    ssp->snt = AddSntRectangle (seg, vContext->from, yOffset,
                                vContext->to + lastLabelWidth, yOffset - scaleHeight - labelLineHeight,
                                ScaleSntDrawProc, (BigScalar) ssp, ScaleSntCleanupProc, 0);
    /*  height += scaleHeight + labelLineHeight + 4;*/
    yOffset = ssp->box.bottom - 1 - 15;
  } else {
    yOffset++;
  }
  yOffset -= DrawBioseqSegments (yOffset, FIP, vContext);
  return initialYOffset - yOffset;
}

/* -=-=-=-==-=-==-=-= GRAPH STUFF =-=-=-=-=-=-=-=-=-*/

typedef struct gphsentdata {
  Boolean        flagIsGUI;
  SegmenT        seg;
  PrimitivE      snt;
  BoxInfo        box;
  Int4           min;
  Int4           max;
  Int4           axis;
  Int4           bottom;
  FloatHi        a;
  FloatHi        b;
  SeqGraphPtr    sgp;
  Uint1          red, green, blue;
} GphSentData, PNTR GphSentPtr;

static void PlotTheWorkingArray (Int4Ptr xp, Int4 len, Int4 scaleX,
                                 Int4 gmin, RecT r, Uint1 uR, Uint1 uG,
                                 Uint1 uB, AttPData PNTR curattrib)
{
  PoinT          pt;
  Int4           i;
  Int4           j;
  Int4           max;
  Int4           min;
  Int4           val;
  Uint1          curR, curG, curB;

  if (curattrib != NULL)
  {
    curR = curattrib->color[0];
    curG = curattrib->color[1];
    curB = curattrib->color[2];
  }
  else
  {
    curR = 0;
    curG = 0;
    curB = 0;
  }

  SelectColor (uR, uG, uB);
  pt.x = r.left;
  i = 0;
  val = (Int4) *xp;
  pt.y = (Int2) (r.bottom - (Int2) (val - gmin));
  if (pt.y < r.top) {
    pt.y = r.top;
  }
  MoveTo (pt.x, pt.y);
  for (; i < len - scaleX; i += scaleX) {
    val = (Int4) *xp;
    min = (Int4) val;
    max = (Int4) val;
    for (j = 1; j < scaleX; j++) {
      xp++;
      val = (Int4) *xp;
      min = MIN (min, (Int4) val);
      max = MAX (max, (Int4) val);
    }
    xp++;
    pt.y = (Int2) (r.bottom - (Int2) (min - gmin));
    if (pt.y < r.top) {
      pt.y = r.top;
    }
    LineTo (pt.x, pt.y);
    pt.y = (Int2) (r.bottom - (Int2) (max - gmin));
    if (pt.y < r.top) {
      pt.y = r.top;
    }
    LineTo (pt.x, pt.y);
    (pt.x)++;
  }
  SelectColor (curR, curG, curB);
  return;
}

static Int4Ptr MakeWorkingSeqGraphInt4Array (Pointer data, Uint1 type,
                                             Int4 pos, Int4 len,
                                             Boolean usescaleflags,
                                             FloatHi a, FloatHi b)
{
  Int4          i;
  Int4Ptr       xpoints, xp;
  FloatHi       fval;
  Int4          ival;
  Byte          bval;
  ByteStorePtr  bs;
  BytePtr       bp;

  xpoints = (Int4Ptr) MemNew ((size_t) (sizeof (Int4) * (len + 10)));

  if (xpoints == NULL)
    return xpoints;

  xp = xpoints;
  switch (type) {
   default:
   case 1:
     for (i = 0; i < len; i++, pos++, xp++) {
      if (usescaleflags) {
        fval = a * (FloatHi) (((FloatHiPtr) data)[pos]) + b;
      } else {
        fval = (FloatHi) (((FloatHiPtr) data)[pos]);
      }
      *xp = (Int4) fval;
    }
    break;
   case 2:
    for (i = 0; i < len; i++, pos++, xp++) {
      if (usescaleflags) {
        ival = (Int4) (a * (FloatHi) (((Int4Ptr) data)[pos]) + b);
      } else {
        ival = (Int4) (((Int4Ptr) data)[pos]);
      }
      *xp = (Int4) ival;
    }
    break;
   case 3:
    bp = MemNew (sizeof (Byte) * (len + 10));
    bs = (ByteStorePtr) data;
    BSSeek (bs, pos, SEEK_SET);
    len = BSRead (bs, (Pointer) bp, len * sizeof (Byte));
    for (i = 0; i < len; i++, pos++, xp++) {
      if (usescaleflags) {
        bval = (Byte) (a * (FloatHi) (bp [i]) + b);
      } else {
        bval = (Byte) (bp [i]);
      }
      *xp = (Int4) bval;
    }
    MemFree (bp);
    break;
  }
  return xpoints;
}

static void GphSentProc (BigScalar calldata, PrimDrawContext pdc)
{
  BioseqPtr    bsp;
  Int4         curScale;
  DrawInfoPtr  dInfoPtr;
  GphSentPtr   gsp;
  Int4         left;
  Int4         len, pos;
  RecT         r;
  Int4         scaleX;
  Int4         scaleY;
  SeqGraphPtr  sgp;
  BoxInfo      tmpBox;
  Int4Ptr      xpoints;
  Int4         tmpscaleX;

  gsp = (GphSentPtr) calldata;
  if (gsp == NULL) return;
  sgp = gsp->sgp;
  if (sgp == NULL || sgp->values == NULL) return;

  dInfoPtr = (DrawInfoPtr) pdc;
  tmpBox = gsp->box;
  if ( (tmpBox.left   > dInfoPtr->scale.worldWindow.right ) ||
       (tmpBox.right  < dInfoPtr->scale.worldWindow.left  ) ||
       (tmpBox.top    < dInfoPtr->scale.worldWindow.bottom) ||
       (tmpBox.bottom > dInfoPtr->scale.worldWindow.top   ) )
    return;

  if ( dInfoPtr->checked == FALSE ) {
    if ( tmpBox.left < dInfoPtr->scale.worldWindow16.left )
      tmpBox.left = dInfoPtr->scale.worldWindow16.left;
    if ( tmpBox.right > dInfoPtr->scale.worldWindow16.right )
      tmpBox.right = dInfoPtr->scale.worldWindow16.right;
    if ( tmpBox.top > dInfoPtr->scale.worldWindow16.top )
      tmpBox.top = dInfoPtr->scale.worldWindow16.top;
    if ( tmpBox.bottom < dInfoPtr->scale.worldWindow16.bottom )
      tmpBox.bottom = dInfoPtr->scale.worldWindow16.bottom;
  }
  scaleX = dInfoPtr->scale.scaleX;
  scaleY = dInfoPtr->scale.scaleY;

  curScale = dInfoPtr->scale.scaleX;
  r.left   = (Int2)((dInfoPtr->scale.offsetX + tmpBox.left )  / curScale);
  r.right  = (Int2)((dInfoPtr->scale.offsetX + tmpBox.right)  / curScale);
  curScale = dInfoPtr->scale.scaleY;
  r.top    = (Int2)((dInfoPtr->scale.offsetY - tmpBox.top   ) / curScale);
  r.bottom = (Int2)((dInfoPtr->scale.offsetY - tmpBox.bottom) / curScale);

  bsp = BioseqFind (SeqLocId (sgp->loc));
  left = GetOffsetInBioseq (sgp->loc, bsp, SEQLOC_LEFT_END);

  if (gsp->flagIsGUI)
  {
    len = tmpBox.right - tmpBox.left;
    pos = tmpBox.left - left;
  }
  else
  {
/* for non-GUI, plot all the sgp values */
    len = sgp->numval;
    pos = 0;
  }

  xpoints = MakeWorkingSeqGraphInt4Array (sgp->values, sgp->flags[2],
                                          pos, len,
                                          (Boolean) (sgp->flags[1] != 0),
                                          gsp->a, gsp->b);
  if (!gsp->flagIsGUI)
  {
    tmpscaleX = sgp->numval / (Int4) (tmpBox.right - tmpBox.left);
    if (tmpscaleX > scaleX)
      scaleX = tmpscaleX;
  }

  PlotTheWorkingArray (xpoints, len, scaleX, gsp->min, r,
                       gsp->red, gsp->green, gsp->blue,
                       &(dInfoPtr->curattrib));

  MemFree (xpoints);
}

static void CleanGSP (BigScalar calldata)
{
  GphSentPtr  gsp;

  gsp = (GphSentPtr) calldata;
  MemFree (gsp);
}

static GphSentPtr AddGphSentinelToPicture (SeqGraphPtr sgp, BioseqPtr bsp,
                                           SegmenT pict, Int4 scaleX,
                                           Int4 top, Int2 start,
                                           Uint1Ptr uRGB)
{
  Int4        axis;
  GphSentPtr  gsp;
  Int4        i;
  Boolean     is_phrap;
  Int4        max;
  Int4        min;
  SegmenT     seg;
  Char        str [32];
  Int4        leftoff, rightoff;

  if (sgp == NULL || bsp == NULL || pict == NULL)
    return NULL;
  gsp = MemNew (sizeof (GphSentData));
  if (gsp == NULL)
    return NULL;
  gsp->flagIsGUI = VibrantIsGUI ();
  if (uRGB != NULL)
  {
    gsp->red   = uRGB[0];
    gsp->green = uRGB[1];
    gsp->blue  = uRGB[2];
  }
  else
  {
    gsp->red   = 0;
    gsp->green = 0;
    gsp->blue  = 0;
  }

  leftoff = GetOffsetInBioseq (sgp->loc, bsp, SEQLOC_LEFT_END);
  rightoff = GetOffsetInBioseq (sgp->loc, bsp, SEQLOC_RIGHT_END);
  if (!gsp->flagIsGUI)
  {
    leftoff /= scaleX;
    rightoff /= scaleX;
  }
  gsp->box.left = leftoff + start;
  gsp->box.right = rightoff - 1 + start;

  gsp->sgp = sgp;
  gsp->a = sgp->a;
  gsp->b = sgp->b;
  is_phrap = (Boolean) (StringICmp (sgp->title, "Phrap Quality") == 0 ||
                        StringICmp (sgp->title, "Phred Quality") == 0);
  switch (sgp->flags [2]) {
    case 1 :
      min = (Int4) sgp->min.realvalue;
      max = (Int4) sgp->max.realvalue;
      axis = (Int4) sgp->axis.realvalue;
      if (sgp->flags [1] != 0) {
        min = (Int4) (sgp->a * ((FloatHi) min) + sgp->b);
        max = (Int4) (sgp->a * ((FloatHi) max) + sgp->b);
      }
      break;
    case 2 :
      min = (Int4) sgp->min.intvalue;
      max = (Int4) sgp->max.intvalue;
      axis = (Int4) sgp->axis.intvalue;
      if (sgp->flags [1] != 0) {
        min = (Int4) (sgp->a * ((FloatHi) min) + sgp->b);
        max = (Int4) (sgp->a * ((FloatHi) max) + sgp->b);
      }
      break;
    case 3 :
      min = (Int4) sgp->min.intvalue;
      max = (Int4) sgp->max.intvalue;
      if (is_phrap) {
        min = MIN (0, min);
        max = MAX (100, max);
      }
      axis = (Int4) sgp->axis.intvalue;
      if (sgp->flags [1] != 0) {
        min = (Int4) (sgp->a * ((FloatHi) min) + sgp->b);
        max = (Int4) (sgp->a * ((FloatHi) max) + sgp->b);
      }
      break;
    default :
      min = (Int4) 0;
      max = (Int4) 100;
      axis = (Int4) 0;
      break;
  }
  gsp->seg = seg = CreateSegment (pict, 0, 0);
  gsp->bottom = top - (max - min) - 20;
  AddSegRect (seg, FALSE, 0);

  if (sgp->title != NULL)  /* StringHasNoText -- vibforms */
  {
    if (StrLen (sgp->title) > 0)
    {
      AddLabel (seg, (gsp->box.left + gsp->box.right) / 2, top,
                sgp->title, SMALL_TEXT, 0, MIDDLE_CENTER, 0);
    }
  }

  top -= 10;
  gsp->box.top = top;
  gsp->box.bottom = gsp->bottom;
  gsp->min = min;
  gsp->max = max;
  gsp->axis = axis;
  gsp->bottom += 10;

  if (is_phrap)
  {
    for (i = 0; i <=100; i += 20) {
      sprintf (str, "%ld", (long) i);
      AddLabel (seg, gsp->box.left, gsp->bottom + i, str,
                SMALL_TEXT, 5, MIDDLE_LEFT, 0);
    }
  }
  else
  {
    sprintf (str, "%ld", (long) max);
    AddLabel (seg, gsp->box.left, top-10, str,
              SMALL_TEXT, 5, MIDDLE_LEFT, 0);
    sprintf (str, "%ld", (long) min);
    AddLabel (seg, gsp->box.left, gsp->bottom-10, str,
              SMALL_TEXT, 5, MIDDLE_LEFT, 0);
    if (min < 0 && max > 0)
    {
      sprintf (str, "%ld", 0L);
      AddLabel (seg, gsp->box.left, gsp->bottom-min-10, str,
                SMALL_TEXT, 5, MIDDLE_LEFT, 0);
    }
  }

  gsp->snt = AddSntRectangle (seg, gsp->box.left, gsp->box.top,
                              gsp->box.right, gsp->box.bottom,
                              GphSentProc, (BigScalar) gsp, CleanGSP, 0);
  return gsp;
}

static void VisitAndListGraphs (SeqGraphPtr sgp, Pointer userdata)

{
  GphSentPtr        gsp;
  ViewerContextPtr  vContext;

  vContext = (ViewerContextPtr) userdata;
  if (vContext == NULL) return;

  if (vContext->gphseg == NULL) {
    vContext->gphseg = CreateSegment (vContext->drawOnMe, 0, 0);
  }

  gsp = AddGphSentinelToPicture (sgp, vContext->BSP, vContext->gphseg,
                                 vContext->viewScale, vContext->gphyOffset,
                                 0, NULL);
  if (gsp == NULL) return;
  vContext->gphheight = MAX (vContext->gphheight, gsp->box.top - gsp->box.bottom);
}

static Int4 DrawGraphs (
  Int4 yOffset,
  ViewerContextPtr vContext
)

{
  vContext->gphheight = 0;
  vContext->gphseg = NULL;
  vContext->gphyOffset = yOffset - 16; /* (workaround) drawing the graphs touches some pixels above gphyOffset  */

  VisitGraphsOnBsp (vContext->BSP, (Pointer) vContext, VisitAndListGraphs);

  return vContext->gphheight + 16;
}

static void ResetFilterState (
  FilterProcessStatePtr FPSP
)

{
  FilterItemPtr FIP;
  ViewerContextPtr vContext;


  if (FPSP == NULL) return;
  vContext = FPSP->vContext;
  FIP = FPSP->currentFIP;
  MemSet (&FPSP->state, 0, sizeof (FilterState));
  if (FIP->Type == FeatureFilter) {
    FPSP->state.feat.currentRFIPblockVNP = vContext->featVNP;
  }
}

static Boolean FilterAndLayout (
  ViewerContextPtr vContext
)

{
  FilterProcessState  FPS;
  FilterItemPtr       FIP;
  FilterPtr           FP;
  AppearancePtr       AP;
  FilterItem          tempFI;
  Int1                featdeftype;
  LayoutAlgorithm     layoutC;
  Int2                height;
  SegmenT             filterSeg, invisibleSeg;
  Uint1               featdefOrder;
  Boolean             emptyFilterGroup;
  Int4                undoCeiling;
  BioseqPtr           BSP;
  SeqAnnotPtr         SAnnP;
  SeqAlignPtr         SAlnP;

  if (vContext == NULL) return FALSE;
  MemSet (&FPS, 0, sizeof (FilterProcessState));
  FPS.vContext = vContext;
  if (vContext->ceiling != NULL) {
    FPS.ceiling = *vContext->ceiling;
  } else {
    FPS.ceiling = -20;
  }
  FPS.featuresProcessed = MemNew (vContext->featureCount * sizeof (Boolean));
  if (FPS.featuresProcessed == NULL && vContext->featureCount != 0) return FALSE;
  AP = vContext->AppPtr;

  FP = vContext->FltPtr;
  for (FPS.currentFilterVNP = FP->FilterItemList;
       FPS.currentFilterVNP != NULL && FPS.currentFilterVNP->data.ptrvalue != NULL;
       FPS.currentFilterVNP = FPS.currentFilterVNP->next) {

    filterSeg = CreateSegment (vContext->topLevelSeg, 0, 0);
    MemSet (FPS.labelSegs, 0, sizeof (FPS.labelSegs));
    MemSet (FPS.drawSegs, 0, sizeof (FPS.drawSegs));
    vContext->drawOnMe = filterSeg;
    emptyFilterGroup = TRUE;
    undoCeiling = FPS.ceiling;
    FPS.featuresProcessedCount = 0;

    FIP = (FilterItemPtr) FPS.currentFilterVNP->data.ptrvalue;

    if (! FIP->GroupBoxColorSet) {
      MemCopy (FIP->GroupBoxColor, AP->GroupBoxColor, sizeof (Uint1 [3]));
    }
    if (! FIP->GroupLabelColorSet) {
      MemCopy (FIP->GroupLabelColor, AP->GroupLabelColor, sizeof (Uint1 [3]));
    }
    if (! FIP->GroupLabelFontSet) {
      FIP->GroupLabelFont = AP->GroupLabelFont;
    }

    if (FIP->DrawItemRect) {
      FPS.ceiling -= 4;
    }
    if (FIP->GroupLabelLoc == LabelOnTop) {
      SelectFont (FIP->GroupLabelFont);
      FPS.ceiling -= LineHeight () + 4;
    }

    switch (FIP->Type) {
      case InvalidFilter:
        break;
      case BioseqFilter:
        height = DrawBioseq (FPS.ceiling, FIP, vContext);
        emptyFilterGroup = FALSE;
        break;
      case GraphFilter:
        height = DrawGraphs (FPS.ceiling, vContext);
        if (height != 0) {
          emptyFilterGroup = FALSE;
        }
        break;
      case AlignmentFilter:
        BSP = vContext->BSP;
        height = 0;
        FPS.currentFIP = FIP;
        for (SAnnP = BSP->annot; SAnnP != NULL; SAnnP = SAnnP->next) {
          if (SAnnP->type != 2) continue;  /* type 2 is an alignment */
          emptyFilterGroup = FALSE;
          SAlnP = (SeqAlignPtr) SAnnP->data;
          ResetFilterState (&FPS);
          FPS.state.align.SAPhead = FPS.state.align.SAPcurrent = SAlnP;
          /*
          fixme: something in here is thrashing the heap
          height = ProcessRows (Layout_Diagonal, &FPS, vContext);
          */
          height += 10;
          FPS.ceiling -= height;
        }
        break;
      case FeatureFilter:
        layoutC = (vContext->overrideLayout != Layout_Inherit) ? vContext->overrideLayout : FIP->LayoutChoice;
        /*
          Some layouts act to the user as if they are single FilterGroups, but internally use multiple
          consecutive FilterGroups
          (currently, Layout_FeatTypePerLine, Layout_FeatTypePerLineGroup, Layout_GroupCorrespondingFeats, and Layout_GroupCorrespondingFeatsRepeat).
        */
        /*
          FeatTypePerLineGroup is equiv. to using PackUpwards several times, with single-feature filteritems
          FeatTypePerLine is similar (but using AllInOneLine)
        */
        switch (layoutC) {
        case Layout_FeatTypePerLine:
        case Layout_FeatTypePerLineGroup:
          FPS.currentFIP = &tempFI;
          MemCopy (&tempFI, FIP, sizeof (FilterItem)); /* copy the filter . . .*/
          MemSet (&tempFI.IncludeFeature, 0, sizeof (tempFI.IncludeFeature));  /* but don't include any features*/
          tempFI.AddTypeToLabel = FALSE;
          for (featdefOrder = 1; featdefOrder < FEATDEF_MAX; featdefOrder++) {
            for (featdeftype = 1; featdeftype < FEATDEF_MAX; featdeftype++) {
              if (FIP->IncludeFeature [featdeftype] == featdefOrder) {
                ResetFilterState (&FPS);
                tempFI.IncludeFeature[featdeftype] = TRUE;
                height = ProcessRows (layoutC, &FPS, vContext);
                if (FPS.featuresProcessedCount != 0) {
                  emptyFilterGroup = FALSE;
                }
                FPS.ceiling -= height;
                tempFI.IncludeFeature[featdeftype] = FALSE;
              }
            }
          }
          height = 0; /* prevent FPS.ceiling from being bumped again */
          break;
        case Layout_GroupCorrespondingFeats:
        case Layout_GroupCorrespondingFeatsRepeat:
          /*
            This uses 3 FilterItems:
              - All gene features (compact)
              - CDS & mRNA (grouped by products)
              - anything else included in the filter as specified by the user (compact)
          */
          FPS.currentFIP = &tempFI;
          MemCopy (&tempFI, FIP, sizeof (FilterItem)); /* copy the filter . . .*/
          MemSet (&tempFI.IncludeFeature, 0, sizeof (tempFI.IncludeFeature));  /* but don't include any features*/
          tempFI.GroupPadding = 0;
          tempFI.IncludeFeature [FEATDEF_GENE] = FIP->IncludeFeature [FEATDEF_GENE];
          ResetFilterState (&FPS);
          height = ProcessRows (Layout_PackUpward, &FPS, vContext);
          FPS.ceiling -= height;

          ResetFilterState (&FPS);
          tempFI.IncludeFeature [FEATDEF_CDS] = FIP->IncludeFeature [FEATDEF_CDS];
          tempFI.IncludeFeature [FEATDEF_mRNA] = FIP->IncludeFeature [FEATDEF_mRNA];
          height = ProcessRows (layoutC, &FPS, vContext);
          FPS.ceiling -= height;

          ResetFilterState (&FPS);
          MemCopy (&tempFI.IncludeFeature, &FIP->IncludeFeature, sizeof (tempFI.IncludeFeature));
          tempFI.GroupPadding = FIP->GroupPadding;
          height = ProcessRows (Layout_PackUpward, &FPS, vContext);
          FPS.ceiling -= height;
          height = 0;
          if (FPS.featuresProcessedCount != 0) {
            emptyFilterGroup = FALSE;
          }

          break;
        default:
          FPS.currentFIP = FIP;
          ResetFilterState (&FPS);
          height = ProcessRows (layoutC, &FPS, vContext);
          if (FPS.featuresProcessedCount != 0) {
            emptyFilterGroup = FALSE;
          }
          break;
        } /* switch (layoutC) */
        break;
    } /* switch (FIP->type) */
    if (emptyFilterGroup) {
      FPS.ceiling = undoCeiling;
      continue;
    } else {
      if (FIP->DrawItemRect) {
        /* AddAttribute (filterSeg, COLOR_ATT, FIP->GroupBoxColor, 0, 0, 0, 0);*/
        /* !!!! use a better invisibility control here !!!! */
        invisibleSeg = CreateSegment (filterSeg, 0, 0);
        SetSegmentVisibleFlag (invisibleSeg, FALSE);
        AddLine (invisibleSeg, vContext->from - 1, undoCeiling - 2, vContext->to + 1, undoCeiling - 2, FALSE, 0);
        undoCeiling -= 4;
      }
      switch (FIP->GroupLabelLoc) {
      default: break;
      case LabelOnTop:
        if (StringHasNoText (FIP->GroupLabel)) break;
        AddAttribute (filterSeg, COLOR_ATT, FIP->GroupLabelColor, 0, 0, 0, 0);
        AddTextLabel (filterSeg, (vContext->from + vContext->to) / 2, undoCeiling,
                      FIP->GroupLabel, FIP->GroupLabelFont, 1, LOWER_CENTER, 0);
        break;
      case LabelOnSide:
        AddTextLabel (filterSeg, 0, FPS.ceiling - height / 2, FIP->GroupLabel,
                      programFont, 1, MIDDLE_RIGHT, 0);
        height = MAX (height, LineHeight () + 3);
        break;
      case LabelOnBottom:
        AddTextLabel (filterSeg, (vContext->from + vContext->to) / 2 , FPS.ceiling - height, FIP->GroupLabel,
                      programFont, 1, LOWER_CENTER, 0);
        height += LineHeight () + 3;
        break;
      }
    }
    if (FIP->DrawItemRect && !emptyFilterGroup) {
      FPS.ceiling -= 20;
    }
    FPS.ceiling -= height + FIP->IntraRowPaddingPixels;
    if (FPS.needFreeList != NULL) {
      ValNodeFreeData (FPS.needFreeList);
      FPS.needFreeList = NULL;
    }
  }
  if (FPS.needFreeList != NULL) {
    ValNodeFreeData (FPS.needFreeList);
    FPS.needFreeList = NULL;
  }

  if (vContext->ceiling != NULL) {
    *vContext->ceiling = FPS.ceiling;
  }
  return TRUE;
}

NLM_EXTERN RelevantFeaturesPtr FreeCollectedFeatures (
  RelevantFeaturesPtr RFP
)

{
  if (RFP == NULL) return NULL;
  if (RFP->sapList) {
    MemFree (RFP->sapList);
  }
  if (RFP->featVNP != NULL) {
    ValNodeFreeData (RFP->featVNP);
  }
  MemFree (RFP);
  return NULL;
}

NLM_EXTERN SegmenT CreateGraphicViewInternal (
  BioseqPtr bsp,
  Int4 from,
  Int4 to,
  Boolean allFeatures,
  RelevantFeaturesPtr feats,
  Int4 scale,
  Int4Ptr ceiling,
  SegmenT topLevel,
  AppearancePtr AP,
  FilterPtr FP,
  LayoutAlgorithm overrideLayout
)

{
  ViewerContext  VC;

  /*
    Removed checks feats->featureCount == 0 || feats->featVNP == NULL
    to allow display of BSP w/0 features on it.
   */

  if (FP == NULL) {
    return NULL;
  }
  if (AP == NULL) {
    return NULL;
  }
  VC.from = MIN (from, to);
  VC.to =   MAX (from, to);
  VC.allFeatures = allFeatures;
  VC.BSP = bsp;
  VC.viewScale = scale;
  VC.sapList = feats->sapList;
  VC.sapCount = feats->sapCount;
  if (topLevel == NULL) {
    VC.drawOnMe = VC.topLevelSeg = CreatePicture ();
  } else {
    VC.drawOnMe = VC.topLevelSeg = topLevel;
  }
  VC.featureCount = feats->featureCount;
  VC.featVNP = feats->featVNP;
  VC.AppPtr = AP;
  VC.FltPtr = FP;
  VC.overrideLayout = overrideLayout;
  VC.seqLength = bsp->length;
  VC.ceiling = ceiling;
  FilterAndLayout (&VC);
  return VC.topLevelSeg;
}

/* returns the 1st segment in a linked list. caller must deallocate it */
NLM_EXTERN SegmenT CreateGraphicViewFromBsp (
  BioseqPtr bsp,
  SeqLocPtr location,
  Int4 scale,
  Int4Ptr ceiling,
  SegmenT topLevel,
  AppearancePtr AP,
  FilterPtr FP,
  LayoutAlgorithm overrideLayout
)

{
  RelevantFeaturesPtr   RFP;
  SegmenT               seg;
  SeqIntPtr             sintp;
  BioseqPtr             parent;
  SeqMgrSegmentContext  context;
  Int4                  from = 0;
  Int4                  to = 0;
  Boolean               allFeatures = TRUE;

  if (location != NULL) {
    bsp = BioseqFindFromSeqLoc (location);
    if (bsp == NULL) return NULL;
    to = bsp->length - 1;

    if (location->choice == SEQLOC_WHOLE) {
      location = NULL; /* no special behavior needed if it's whole */
    } else if (location->choice == SEQLOC_INT) {
      sintp = (SeqIntPtr) location->data.ptrvalue;
      if (sintp != NULL && sintp->from == 0 && sintp->to == bsp->length - 1) {
        location = NULL;
      } else {
        from = sintp->from;
        to = sintp->to;
        allFeatures = FALSE;
      }
    }
  } else if (bsp != NULL) {
    to = bsp->length - 1;
  }
  if (bsp == NULL) return NULL;
  parent = SeqMgrGetParentOfPart (bsp, &context);
  if (parent != NULL) {
    from = context.cumOffset;
    to = from + context.to - context.from;
    allFeatures = FALSE;
    bsp = parent;
  }
  RFP = CollectFeatures (bsp);
  if (RFP == NULL) return NULL;
  seg = CreateGraphicViewInternal (bsp, from, to, allFeatures, RFP, scale, ceiling, topLevel, AP, FP, overrideLayout);
  FreeCollectedFeatures (RFP);
  return seg;
}

NLM_EXTERN SegmenT CreateGraphicView (
  BioseqPtr bsp,
  SeqLocPtr location,
  Int4 scale,
  CharPtr styleName,
  CharPtr filterName,
  CharPtr overrideLayoutName
)

{
  ViewerConfigsPtr  myVCP;
  FilterPtr         FP;
  AppearancePtr     AP;
  LayoutAlgorithm   overrideLayout;
  Uint1             i;

  if (bsp == NULL && location == NULL) return NULL;
  myVCP = GetGraphicConfigParseResults ();

  FP = FindFilterByName (filterName, myVCP);
  AP = FindAppearanceByName (styleName, myVCP);
  i = StringIndexInStringList (overrideLayoutName, LayoutStrings);
  if (i >= 0 && i < DIM(LayoutValues)) {
    overrideLayout = LayoutValues[i];
  } else {
    overrideLayout = Layout_Inherit;
  }
  return CreateGraphicViewFromBsp (bsp, location, scale, NULL, NULL, AP, FP, overrideLayout);
}

NLM_EXTERN Uint2 GetAppearanceCount (void)

{
  ViewerConfigsPtr VCP;

  VCP = GetGraphicConfigParseResults ();
  if (VCP == NULL) return 0;
  return VCP->AppearanceCount;
}

NLM_EXTERN Uint2 GetFilterCount (void)

{
  ViewerConfigsPtr VCP;

  VCP = GetGraphicConfigParseResults ();
  if (VCP == NULL) return 0;
  return VCP->FilterCount;
}

NLM_EXTERN Uint2 GetLayoutCount (void)

{
  return DIM (LayoutStrings);
}

NLM_EXTERN CharPtr PNTR GetStyleNameList (void)

{
  ViewerConfigsPtr VCP;

  VCP = GetGraphicConfigParseResults ();
  if (VCP == NULL) return NULL;
  return VCP->AppearanceNameArray;
}

NLM_EXTERN CharPtr PNTR GetFilterNameList (void)

{
  ViewerConfigsPtr VCP;

  VCP = GetGraphicConfigParseResults ();
  if (VCP == NULL) return NULL;
  return VCP->FilterNameArray;
}

NLM_EXTERN CharPtr PNTR GetLayoutNameList (void)

{
  return LayoutStrings;
}


/* -=-=-=-=-=-=-=- append default-style.c contents after this point -=-=-=-=-=-=-=-=- */

/* default-style.c : creates a default style for the new graphic viewer
  This is an automatically generated file, which came from an asn2gph configuration file (asn2gph.default)
  It might be better to edit the input file and then re-run the script create-default-style.tcl than to
  edit this file directly.
*/


typedef struct configFileLine {
  CharPtr key, value;
} ConfigFileLine, PNTR ConfigFileLinePtr;

typedef struct staticConfigFile {
  ConfigFileLinePtr lines;
  CharPtr sectionName;
} StaticConfigFile, PNTR StaticConfigFilePtr;


/* [Styles] */
static ConfigFileLine defaultStyleLines1 [] = {
  {"style00", "defaultStyle"},
  {NULL, NULL}
};

/* [defaultStyle.master] */
static ConfigFileLine defaultStyleLines2 [] = {
  {"name", "Default"},
  {"maxarrowscale", "200"},
  {"minarrowpixels", "5"},
  {"shadesmears", "false"},
  {"color", "black"},
  {"labelfont", "program"},
  {"labelcolor", "black"},
  {"label", "above"},
  {"grouplabelfont", "program"},
  {"grouplabelcolor", "dark gray"},
  {"groupboxcolor", "gray"},
  {"displaywith", "box"},
  {"height", "8"},
  {"gap", "line"},
  {"showarrow", "no"},
  {"showtype", "yes"},
  {"showcontent", "yes"},
  {"shadesmears", "false"},
  {NULL, NULL}
};

/* [defaultStyle.bioseq] */
static ConfigFileLine defaultStyleLines3 [] = {
  {"label", "top"},
  {"format", "accn"},
  {"scale", "true"},
  {"labelfont", "program"},
  {"scalefont", "small"},
  {"height", "10"},
  {"scaleheight", "10"},
  {"color", "0, 0, 0"},
  {"labelcolor", "64, 64, 255"},
  {"scalecolor", "32, 32, 32"},
  {NULL, NULL}
};

/* [defaultStyle.gene] */
static ConfigFileLine defaultStyleLines4 [] = {
  {"label", "above"},
  {"color", "blue"},
  {"labelcolor", "blue"},
  {"showarrow", "true"},
  {NULL, NULL}
};

/* [defaultStyle.mRNA] */
static ConfigFileLine defaultStyleLines5 [] = {
  {"label", "above"},
  {"color", "cyan"},
  {"labelcolor", "cyan"},
  {"showarrow", "true"},
  {"gap", "line"},
  {NULL, NULL}
};

/* [defaultStyle.cds] */
static ConfigFileLine defaultStyleLines6 [] = {
  {"label", "above"},
  {"color", "magenta"},
  {"labelcolor", "magenta"},
  {"showarrow", "true"},
  {"gap", "angle"},
  {NULL, NULL}
};

/* [defaultStyle.tRNA] */
static ConfigFileLine defaultStyleLines7 [] = {
  {"label", "above"},
  {"color", "green"},
  {"labelcolor", "green"},
  {"showarrow", "true"},
  {"gap", "line"},
  {NULL, NULL}
};

/* [defaultStyle.imp] */
static ConfigFileLine defaultStyleLines8 [] = {
  {"showcontent", "no"},
  {"color", "gray"},
  {"labelcolor", "gray"},
  {NULL, NULL}
};

/* [defaultStyle.exon] */
static ConfigFileLine defaultStyleLines9 [] = {
  {"showcontent", "no"},
  {"color", "dark cyan"},
  {"labelcolor", "dark cyan"},
  {NULL, NULL}
};

/* [defaultStyle.intron] */
static ConfigFileLine defaultStyleLines10 [] = {
  {"showcontent", "no"},
  {"color", "light gray"},
  {"labelcolor", "light gray"},
  {NULL, NULL}
};

/* [defaultStyle.bond] */
static ConfigFileLine defaultStyleLines11 [] = {
  {"displaywith", "cappedline"},
  {NULL, NULL}
};

/* [Filters] */
static ConfigFileLine defaultStyleLines12 [] = {
  {"filter00", "defaultFilt"},
  {"maxlabelscale", "200"},
  {"grouppadding", "2"},
  {"rowpadding", "2"},
  {NULL, NULL}
};

/* [defaultFilt] */
static ConfigFileLine defaultStyleLines13 [] = {
  {"name", "Default"},
  {"layout", "compact"},
  {"group1", "defaultFilt-gene-cds-prot-mrna"},
  {"group2", "defaultFilt-other-rnas"},
  {"group3", "defaultFilt-exon-intron-label"},
  {"group4", "defaultFilt-variations"},
  {"group5", "defaultFilt-conflicts"},
  {"group6", "defaultFilt-impfeats"},
  {"group7", "defaultFilt-everything-else-label"},
  {NULL, NULL}
};

/* [filters.defaultFilt-gene-cds-prot-mrna] */
static ConfigFileLine defaultStyleLines14 [] = {
  {"feature1", "gene"},
  {"feature2", "cds"},
  {"feature3", "prot"},
  {"feature4", "mrna"},
  {"name", "Gene-mRNA-CDS-Prots"},
  {"grouplabel", "none"},
  {"layout", "geneproducts"},
  {"showtype", "yes"},
  {"showcontent", "yes"},
  {"label", "above"},
  {NULL, NULL}
};

/* [filters.defaultFilt-other-rnas] */
static ConfigFileLine defaultStyleLines15 [] = {
  {"feature1", "rna"},
  {"name", "Structural RNAs"},
  {"grouplabel", "none"},
  {"label", "above"},
  {NULL, NULL}
};

/* [filters.defaultFilt-exon-intron-label] */
static ConfigFileLine defaultStyleLines16 [] = {
  {"feature1", "exon"},
  {"feature2", "intron"},
  {"name", "Introns and Exons"},
  {"grouplabel", "none"},
  {"label", "above"},
  {NULL, NULL}
};

/* [filters.defaultFilt-variations] */
static ConfigFileLine defaultStyleLines17 [] = {
  {"feature1", "variation"},
  {"name", "Variations"},
  {"groupbox", "true"},
  {"boxcolor", "red"},
  {"grouplabel", "above"},
  {"layout", "smear"},
  {"showtype", "no"},
  {"showcontent", "no"},
  {NULL, NULL}
};

/* [filters.defaultFilt-conflicts] */
static ConfigFileLine defaultStyleLines18 [] = {
  {"feature1", "conflict"},
  {"name", "Conflicts"},
  {"groupbox", "true"},
  {"boxcolor", "dark red"},
  {"grouplabel", "above"},
  {"layout", "smear"},
  {"showtype", "no"},
  {"showcontent", "no"},
  {NULL, NULL}
};

/* [filters.defaultFilt-impfeats] */
static ConfigFileLine defaultStyleLines19 [] = {
  {"feature1", "imp"},
  {"name", "Import Features"},
  {"grouplabel", "none"},
  {"label", "above"},
  {NULL, NULL}
};

 /* [filters.defaultFilt-everything-else-label] */
static ConfigFileLine defaultStyleLines20 [] = {
  {"feature1", "everything"},
  {"grouplabel", "none"},
  {"label", "above"},
  {NULL, NULL}
};


static StaticConfigFile defaultStyle [] = {
  {defaultStyleLines1, "Styles"},
  {defaultStyleLines2, "defaultStyle.master"},
  {defaultStyleLines3, "defaultStyle.bioseq"},
  {defaultStyleLines4, "defaultStyle.gene"},
  {defaultStyleLines5, "defaultStyle.mRNA"},
  {defaultStyleLines6, "defaultStyle.cds"},
  {defaultStyleLines7, "defaultStyle.tRNA"},
  {defaultStyleLines8, "defaultStyle.imp"},
  {defaultStyleLines9, "defaultStyle.exon"},
  {defaultStyleLines10, "defaultStyle.intron"},
  {defaultStyleLines11, "defaultStyle.bond"},
  {defaultStyleLines12, "Filters"},
  {defaultStyleLines13, "defaultFilt"},
  {defaultStyleLines14, "filters.defaultFilt-gene-cds-prot-mrna"},
  {defaultStyleLines15, "filters.defaultFilt-other-rnas"},
  {defaultStyleLines16, "filters.defaultFilt-exon-intron-label"},
  {defaultStyleLines17, "filters.defaultFilt-variations"},
  {defaultStyleLines18, "filters.defaultFilt-conflicts"},
  {defaultStyleLines19, "filters.defaultFilt-impfeats"},
  {defaultStyleLines20, "filters.defaultFilt-everything-else-label"},
  {NULL, NULL}
};


static void InitializeDefaultStyle (
  CharPtr configFileName
)

{
  Uint2             sectionNum, lineNum;
  ConfigFileLinePtr lines;
  CharPtr           sectionName;

  for (sectionNum = 0; defaultStyle [sectionNum].lines != NULL; sectionNum++) {
    lines = defaultStyle [sectionNum].lines;
    sectionName = defaultStyle [sectionNum].sectionName;
    for (lineNum = 0; lines [lineNum].key != NULL && lines [lineNum].value != NULL; lineNum++) {
      TransientSetAppParam (configFileName, sectionName, lines [lineNum].key, lines [lineNum].value);
    }
  }
}

/* End of automatically generated file. */

