/*   ingengraph.h
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
* File Name:  ingengraph.h
*
* Author:  Fasika Aklilu
*
* Version Creation Date:   4/26/01
*
* $Revision: 6.1 $
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

#ifndef _INGENGRAPH_
#define _INGENGRAPH_

#ifdef __cplusplus
extern "C" {
#endif

#include <ncbi.h>
#include <vibrant.h> 
#include <picture.h> 
#include <viewer.h> 
#include <objmgr.h> 
#include <sequtil.h>
#include <explore.h>
#include <gather.h>
#include <seqport.h>
#include <sqnutils.h>
#include <alignmgr.h>


/******************************************************
 *
 *  defines for the Genome Viewer of UDV
 *
******************************************************/
#define FEATURE_SEGMENT_MAXSIZE 500 /*a segment contains a max of 500 features*/
#define Ing_FEAT_LINEDECAL 20/* was 8 */ /*on the Feature Viewer, the line height is 8 pixels*/
#define Ing_FEAT_HEIGHT 2 /*on the Feature Viewer the height of a Feature is 3 pixels*/
#define LEVEL1_MAXNUM  20 /* top level segment contains max of 20 level 2 segments */
#define LEVEL2_MAXNUM  20 /* second level segment contains max of 20 level 3 segments */
#define COMPRESS_SCALE 500

#define Ing_black       121
#define Ing_red         122
#define Ing_blue        123
#define Ing_green       124
#define Ing_cyan        125
#define Ing_magenta     126
#define Ing_purple      127
#define Ing_yellow      128
#define Ing_grey        129
#define Ing_MAX         130

#define Ing_GENSCAN   110
#define ALIGN_ANNOT    111 /* alignment from running blast vs. db */
#define Ing_SPIDEY    112
#define ALIGN_FILE     113/* alignment from asn file */
#define ALIGN_BLAST    114/* alignment annotated on sequence */
#define ALIGN_BLAST2SEQ    115 /* alignment from running blast2seq */
#define Ing_ORF       116
#define Ing_SUMMARY   117
#define Ing_FEAT_TABLE 118
#define Ing_TBDL    119
#define Ing_SEQIDLISTFORSPIDEY 120
#define Ing_SEQIDLISTFORBLAST 121

/******************************************************
 *
 *  data structures for the Genome Viewer 
 *
******************************************************/


typedef struct gvpopfeat{/*data structure used to populate features viewer*/
  SegmenT   pictRuler;
  SegmenT   pictMain;
  SegmenT   CurrentSeg;
  SegmenT   seg1;
  SegmenT   seg2;
  Uint1Ptr  PNTR pClr;
  Boolean    showLabels;
  Uint4      nSegs; /* number of subsegments below Current segment */
  Uint4      nPrims;  /* prim counter */
  Uint4      nLevel3; /* level 3 segments counter */
  Uint4      nLevel2; /* level 2 segments counter */
  Uint4      nLevel1; /* level 1 segments counter */
  ValNodePtr PopRowsList;/* track features */
  Int4      nPopRows;      /* feature row counter */
  Int4      yBase;
  Int4      scaleX;         
  Int4      max_label;      /* width of largest label */
  SeqLocPtr  slp;
  Int4       left;
  Int4       right;
} IngPopFeat, PNTR IngPopFeatPtr;


typedef struct gvgraphdata{
	/*font size*/
	Int2 cxChar;
	Int2 cyChar;
	/*values used to draw the rulers*/
	Int2 SegBoxHeight;
	Int2 SegRulerTick;
	Int2 SegRulerIn;
	Int4 DecalSegRuler;
	Int4 DecalZoomRuler;
   Uint1Ptr PNTR pClr;
} IngGraphData, PNTR IngGraphDataPtr;


typedef struct gvpopnames{/*used to populate the Segments in Viewer*/
  VieweR          viewer;
	SegmenT        seg;
  Uint1Ptr        seqbuf;
  Int4            scaleX;
  Uint1           comp;
  SeqLocPtr       slp;
  BioseqPtr       bsp;
  ValNodePtr      slp_list;
  Boolean         bPopSlp;
  Boolean         bShowGC;
  Int4            idx;
  Int4            left, right;
  IngGraphDataPtr GrData;
} IngExploreSegs, PNTR IngExploreSegsPtr;



  typedef struct ingseqannotdata{
    Int4         index;
    CharPtr PNTR names;
    SeqLocPtr    slp;
    ValNodePtr   vnp;
  }IngSeqAnnotData;
  


  typedef struct ing_entitylist{
    ValNode *Sips;               /* list of sips in SeqEntry */
    Uint2    entityID;
    Uint2    itemID;
    Int4     bspcount;
    struct   ing_entitylist PNTR next;
  } IngEntity, PNTR IngEntityPtr;
  

  typedef struct win_data{
    MenU      Edit;
    MenU      Run;
    MenU      Options;
    MenU      Analysis;
    MenU      Annotate;
    MenU      Feat;
    GrouP     Goto;
    IteM      labels;
    IteM      genesorexons;
    IteM      GC;
  } IngWinData, PNTR IngWinDataPtr;

typedef struct ing_genomeviewer{
  WindoW      hMain;
  WindoW      hReport;
  GrouP       gMain;
  /*graphic elements*/
  VieweR      vRuler1;
  VieweR      vTop;
  VieweR      vRuler2;
  VieweR      vBottom;
  SegmenT     pictRuler1;
  SegmenT     pictTop;
  SegmenT     pictRuler2;
  SegmenT     pictBottom;
  VieweR      vZoom;
  SegmenT     pictZoom;
  SegmenT     segbox;
  SegmenT     segArrow;
  SegmenT     seglabel;
  SegmenT     segGCRect;
  SegmenT     segORFRect;
  PrompT      featInfo[3];
  CharPtr     defline;
  GrouP       deflineg;
  PopuP       pageControl;
  Uint2       sel_entityID;
  Uint2       sel_itemID;
  GrouP       Buttons;
  IteM        item_usenetwork;
  IngGraphData GrData;
  Boolean     update;
  Boolean     bOverviewSelected;
  /*data to display*/
  Uint2       procID;
  Uint2       userKey;
  Uint2       entityID;
  Uint2       itemID;
  Uint2       itemType;
  BioseqPtr   bsp;
  Boolean     bspIsLocked;
  Int4        from;
  Int4        to;
  Int4        SeqLength;/* length of user specified region*/
  Int4        bspLength; /* total sequence length */
  Int4        scale_index;
  Int4Ptr     scale_array;
  Int4        scaleX;
  Int4        maxScaleX;
  Boolean     bMaxScale;
  IngEntityPtr entity_list;
  Int4        numseqs;
  Int4        ZoomedVal;
  Int4        last_feat_position;
  Int4        xposition;
  Int2        cxChar;
  Char        title[50];
  Boolean     bLabels;
  Boolean     bGenes;
  Boolean     user_defined;
  SeqLocPtr   slp;

  /* variables for imported features -- from sequence analysis */
  Boolean     bDustExists;
  Boolean     bOrfExists;
  Boolean     bShowGC;
  Pointer     data;
  Uint1Ptr    seqbuf;
  Uint1       comp; 
  IngWinDataPtr d_Win;
  LisT bsp_list;          /* listbox of bsps */
  SeqIdPtr    cur_sip;
  Int2        cur_target;
  IteM        imemory;
  /* data info */
  Int4        filetype; /* reading in output files */
  Int2        nAnnotFeats;
  /* zoom */
  PrimitivE   zoomFbar;
	/*internal use only*/
  Int4       idx;
  Boolean bSegmented;
  /* progress notice vars */
  CharPtr    progress_string;
  Boolean    progress_running;
  Int4       progress_timer;
  PrompT     progress_prompt;
  Char       progress_counter[20];
} IngGenomeViewer , * IngGenomeViewerPtr;



Boolean IngfeatDefFilter[FEATDEF_MAX];/* filter parameter for SeqMgrExploreFeatures() */

Boolean IngfeatDefTrack2[FEATDEF_MAX];/* used by Ing_AddToOverviewPage()  when scale is greater than COMPRESS_SCALE */

Uint1   IngfeatDefTrack[FEATDEF_MAX];/* track features 1=exists 2=exists and is visible */

/*******************************************************************************

	Static Function Declarations

*******************************************************************************/
extern void Ing_AddGCRect(SegmenT seg, BioseqPtr bsp, Uint2 entityID, Uint2 itemID, Uint2 itemtype, Uint1Ptr seq, Int4 left, Int4 top, Int4 right, Int4 bottom, Int4 scaleX,  Uint1 strand, Boolean needs_label, Int4 idx, Boolean bShowGC);

extern SegmenT Ing_PopulateSequinGraphic(SegmenT seg, BioseqPtr bsp, Uint2 entityID, Uint2 itemID, Int4 scaleX);

extern void Ing_AddRuler(SegmenT seg, Int4 height, Int4 xstart,Int4 xstop, Int4 scaleX, Int4 scaleY, Boolean add_whitespace);

extern void Ing_InitGrData(IngGraphDataPtr gdp);

extern Uint4 Ing_PopFeaturesPage(BioseqPtr bsp, SegmenT pictBottom, Int4 left, Int4 right, Uint1Ptr PNTR pClr, Int4 scaleX, Int4 start_row, Boolean Labels, Boolean bGenes);
extern Uint4 Ing_PopOverviewPage(BioseqPtr bsp, SegmenT pictTop, Int4 left, Int4 right, Uint1Ptr PNTR pClr, Int4 maxScaleX);

extern Uint1Ptr PNTR Ing_BuildColorTable(void);
  extern Uint1Ptr PNTR Ing_FreeColorTable(Uint1Ptr PNTR pClr);

extern Uint1Ptr  Ing_PutColor(Uint1 r, Uint1 g, Uint1 b);

extern void Ing_GetSeqLocations(IngGenomeViewerPtr igvp, BioseqPtr bsp);

  extern Boolean LIBCALLBACK Ing_ExploreSegments(SeqLocPtr slp, SeqMgrSegmentContextPtr context);

extern void Ing_InitfeatDefFilter(void);
extern void Ing_SearchAli(SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent);
extern Uint8 Ing_BigEncodeIdxFeat (Uint4 val1,Uint4 val2);
extern void  Ing_BigDecodeIdxFeat (Uint8 index_g, Uint4Ptr val1, Uint4Ptr val2);
extern void Ing_PopOverviewRuler(IngGenomeViewerPtr igvp, BioseqPtr bsp);
extern Int4 Ing_GetValue (TexT t);



#ifdef __cplusplus
}
#endif

#endif /* ndef _INGENGRAPH_ */
