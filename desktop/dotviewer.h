/* dotplot.h*/

#ifndef _DOTVIEWER_
#define _DOTVIEWER_


  /****************************************************************************

      INCLUDE SECTION                                                                       
   ***************************************************************************/
#include <vibrant.h>
#include <picture.h>
#include <picturep.h>
#include <viewer.h>
#include <tofasta.h>
#include <ncbidraw.h>
#include <seqport.h>
#include <sequtil.h>
#include <blastkar.h>
#include <blastpri.h>
#include <blastdef.h>
#include <objseq.h>
#include <salpstat.h>
#include <seqmgr.h>
#include <seqgraph.h>
#include <salpstat.h>
#include <actutils.h>
#include <alignmgr.h>
#include <actutils.h>
#include <dotseq.h>


#ifdef __cplusplus
extern "C" {
#endif


/****************************************************************************

      DEFINES SECTION                                                                 
 ***************************************************************************/

#define REGISTER_RunDotPlot  ObjMgrProcLoad (OMPROC_EDIT, "RunDotPlot", "RunDotPlot", OBJ_MAX, 0, OBJ_VIBRANT_PICTURE, 0, NULL, Run_DotPlot, PROC_PRIORITY_HIGHEST);
#define REGISTER_RunDiagPlot  ObjMgrProcLoad (OMPROC_EDIT, "RunDiagPlot", "RunDiagPlot", OBJ_MAX, 0, OBJ_VIBRANT_PICTURE, 0, NULL, Run_DiagPlot, PROC_PRIORITY_HIGHEST);


 /****************************************************************************

      DATA STRUCTURE SECTION                                                               
 ***************************************************************************/

/*******************************************************
* temporary form structures                            *
********************************************************/
typedef struct dot_blastinfo
{
   BioseqPtr  bsp1;
   BioseqPtr  bsp2;
   GrouP      localorglobal;
   GrouP      progname;
   GrouP      gapped;
   TexT       eval;
   TexT       wordsize;
   ButtoN     newblast;
   ButtoN     maskrep;
   ButtoN     masksimple;
} DOTblastinfo, PNTR DOTblastinfoPtr;

typedef struct dot_displayinfo
{
  GrouP      BLASTorDotPlot;
  GrouP      grid;
} DOTdisplayinfo, PNTR DOTdisplayinfoPtr;


typedef struct dot_paramsinfo
{
  TexT       xstart;
  TexT       xstop;
  TexT       ystart;
  TexT       ystop;
  TexT       word_size;
  TexT       tree_limit;
} DOTparamsinfo, PNTR DOTparamsinfoPtr;



/*******************************************************
* storage for  alignment display                       *
********************************************************/

typedef struct saln {
  Int4 q_start;
  Int4 q_stop;
  Int4 s_start;
  Int4 s_stop;
  Int4 primID;
} DOTAln, PNTR DOTAlnPtr;

typedef struct algn{
  Int4  xlen;
  Int4  ylen;
  Int4  xstart;
  Int4  ystart;
  Int4  Fh;
  WindoW  w;
  VieweR  v;
  SeqAlignPtr sap;
  SegmenT pict;
  SegmenT seg1;
  PopuP scale;
  Int4  scaleValue;
  Int4  scaleIndex;
  Boolean showLabels;
  Boolean do_scale;
  PrompT      Infopanel;
  CharPtr title;
  SeqIdPtr  sip;
  Int4  HORZ_MARGIN;
  Int4  VERT_MARGIN;
  Int4  index; /* num. of alignments */
  DOTAlnPtr    PNTR  Alnlist;
} DOTAlignInfo, PNTR DOTAlignInfoPtr;


/*******************************************************
* scroll and threshold data                            *
********************************************************/

  typedef struct scrolldata {
    /* threshold bar info */
    BaR      ScrollBar; 
    Int4     TrampPos;
    /* main window scroll information */
    Int4     YScrlMax;
    Int4     XScrlMax;
    Int4     YScrlPos;
    Int4     XScrlPos;
    Int4     UnitY;
    Int4     UnitX;
    Int4     TotUnitsY;
    Int4     TotUnitsX;
    Int4     PgWdth;
    Int4     PgLen;
    Int4     YScrlPage;
    Int4     XScrlPage;
    /* scroll positions in pixel coordinates */
    Int4     HFrom; 
    Int4     VFrom;
  } DOTScrollData, PNTR DOTScrollDataPtr;


/*******************************************************
* store sequence data on selected region for new display*
*                                                       *
********************************************************/
typedef struct seqviewr{
  Boolean   do_scale;
  PopuP     scale;
  Int4      scaleIndex;
  Int4      scaleValue;
  GrouP     Labels;
  Boolean   showLabels;
  PoinT     old_pt;
  Int2     old_primID;
  SegmenT   pict1;
  SegmenT   pict2;
  SegmenT   seg1;
  DOTAlnPtr   salp;
  VieweR    v1;
  VieweR    v2;
  WindoW    w;
  Int4      barp;
} DOTSeqViewr, PNTR DOTSeqViewrPtr;


/*******************************************************
* main vibrant data structure                          *
********************************************************/
typedef struct dotvibdata{
  WindoW MainWin;
  WindoW ChildWin;
  Int4   curr_slen;
  Int4   curr_qlen;
  PaneL  panel;
  Boolean showGrid;
  Boolean     showDotPlot;
  Boolean     showALIGN;
  PrompT      Infopanel;
  Char        iInfo[255];
  Int4        comp; 
  Int4        originalcomp;
  /* Alignment options */
  Boolean     Blast2Seq_show;
  MenU        displayOpts1;
  ChoicE        displayOpts2;
  DOTAlignInfoPtr alp;
  /* second window elements */
  DOTSeqViewrPtr sv;
  DOTScrollData  sdp;
  Uint2       selectMode;
  VoidPtr     data;
  FonT        Fnt;
  Int4        Fh;
  Int4        charw; 
  Int4        VERT_MARGIN;
  Int4        HORZ_MARGIN;
  DOTMainDataPtr     mip;
} DOTVibData, PNTR DOTVibDataPtr;


/*******************************************************
* structures to store feature data for the feature viewer *
********************************************************/

typedef struct dot_feat_info {
   CharPtr  label;
   Int4     left;
   Int4     right;
   Uint1    strand;
   Int2     type;
   struct dot_feat_info PNTR next;
} DOTFeat, PNTR DOTFeatPtr;


typedef struct dot_row {
   DOTFeatPtr dfp;
} DOTRow, PNTR DOTRowPtr;

typedef struct dot_featindex {
  CharPtr   label;
  Boolean   present;
  Boolean   show;
  Int2      deref; 
} DOTFeatIndex, PNTR DOTFeatIndexPtr;


typedef struct dot_popfeat{
  VieweR TopParentView;
  SegmenT TopParentSeg;
  DOTRowPtr drp;
  DOTFeatIndexPtr featindex;
  DOTFeatPtr dfp_cur;
  Int4       nfeats;
  Int2       fontHeight;
} DOTPopFeat, PNTR DOTPopFeatPtr;


typedef struct dot_feats_list{
  WindoW FeatWin;
  WindoW hFeatDlg;
  VieweR Query;
  SegmenT segQuery;
  SegmenT segQName;
  SegmenT segQCursor;
  VieweR  Subject;
  SegmenT segSubject;
  SegmenT segSName;
  SegmenT segSCursor;
  Int4    qFeatscount;
  Int4    sFeatscount;
  Int2    fontHt;
  Int4    vert_Qpos;
  Int4    vert_Spos;
  Uint2     procID;
  Uint2     userKey;
  Uint2     entityID;
  Uint2     itemID;
  LisT      featList;
  PrompT     QInfo;
  PrompT     SInfo;
  DOTRowPtr  query_drp;
  DOTRowPtr  subject_drp;
  DOTFeatIndexPtr featindex; 
  Int4       numrows; /* list of feature types */
  VoidPtr     data;
  DOTMainDataPtr   mip;
} DOTFeatList, PNTR DOTFeatListPtr;

typedef struct dot_sel_feat{
  Int4      feat_num;
  struct dot_sel_feat PNTR next;
} DOTSelFeat, PNTR DOTSelFeatPtr;


/*******************************************************
* store sequence data on selected region for new display*
*                                                       *
********************************************************/


typedef struct selectdata{
  Int4      q_start;
  Int4      s_start;
  Int4      q_stop;
  Int4      s_stop;
  Int4      qlen;
  Int4      slen;
  RecT      rcS;  
  RecT      old_rcS; 
  RecT      rcP; 
  Boolean   selected;
  Boolean   rm_lastselected;
  Int2      H_pos;
  Int2      V_pos;
  DOTVibDataPtr vdp;
} DOTSelData, PNTR DOTSelDataPtr;



 /****************************************************************************

      FUNCTION DECLARATIONS                                                               
 ***************************************************************************/

Int2 LIBCALLBACK Run_DotPlot (Pointer data);
Int2 LIBCALLBACK Run_DiagPlot (Pointer data);

extern void DOT_MakeMainViewer(DOTMainDataPtr vdp, SeqAlignPtr sap);
extern void DOT_AlignPlotGivenSeqAlign(SeqAlignPtr sap);

#ifdef __cplusplus
}
#endif

#endif /* ndef _DOTVIEWER_ */

