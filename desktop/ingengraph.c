/*   ingengraph.c
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
* File Name:  ingengraph.c
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

#include <ingengraph.h>


#define ING_MAX_VIEWER_HEIGHT  600
#define ARROWS_SCALE   200
#define LABELS_SCALE   100

/*******************************************************************************

	Static Function Declarations

*******************************************************************************/

static void Ing_GetNamesZoomedVal(VieweR vNames,Int4Ptr from, Int4Ptr to,
	Int4Ptr scaleX,Int4 MaxLength);
static void Ing_GetLeftandRight(BioseqPtr bsp, SeqLocPtr slp, Int4Ptr left, Int4Ptr right);

/******************************************************************************

  Function : Ing_InitfeatDefFilter(), Ing_AddStringToSad(), Ing_SearchAli(),Ing_GetAlignColor(), Ing_GetCurrentSegment()
 
  Purpose :  function to retrieve data for drawing

*******************************************************************************/
extern void Ing_InitfeatDefFilter(void)
{
  MemSet(IngfeatDefFilter,1,sizeof(Boolean)*(FEATDEF_MAX));
  IngfeatDefFilter[FEATDEF_COMMENT]=FALSE;
  IngfeatDefFilter[FEATDEF_BIOSRC]=FALSE;
  IngfeatDefFilter[FEATDEF_PUB]=FALSE;
  IngfeatDefFilter[FEATDEF_ORG]=FALSE;
 
}

static Uint1Ptr Ing_GetAlignColor(CharPtr name, Uint1Ptr PNTR pClr)
{
  if (StringCmp(name, "Spidey")==0)
    return pClr[Ing_SPIDEY];
  if (StringCmp(name, "Blast 2 seqs")==0)
    return pClr[ALIGN_BLAST2SEQ];
  if (StringCmp(name, "Blast")==0)
    return pClr[ALIGN_BLAST];
  if (StringCmp(name, "Blast file")==0)
    return pClr[ALIGN_FILE];
  return pClr[ALIGN_ANNOT];
}

static CharPtr PNTR Ing_AddStringToSad(CharPtr PNTR names, SeqAnnotPtr sanp, Int4 index)
{
  Int4    i;
  CharPtr new_name;
  CharPtr PNTR head;

  if (sanp->name)
    new_name=sanp->name;
  else if (sanp->desc || sanp->desc->data.ptrvalue)
    new_name=sanp->desc->data.ptrvalue;
  else 
    return NULL;
  if (names==NULL){
    head=(CharPtr PNTR)MemNew(sizeof(CharPtr));
    *head=new_name;
    return head;
  }
  else {
    head=(CharPtr PNTR)MemNew(sizeof(CharPtr)*index);
    for (i=0; i<index-1; i++){
      head[i]=names[i];
    }
    
    head[index-1]=new_name;
    MemFree(names);
    return head;
  }
}

extern void Ing_SearchAli(SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)
{
ValNodePtr PNTR vnpp=NULL;
BioseqPtr     bsp=NULL;
BioseqSetPtr  bssp=NULL;
SeqAnnotPtr   sanp=NULL;
SeqAlignPtr   salp=NULL;
IngSeqAnnotData * sadp=NULL; 


 sadp=(IngSeqAnnotData *)mydata;
 if (IS_Bioseq(sep)) {
   bsp = (BioseqPtr)sep->data.ptrvalue;
   sanp = bsp->annot;
 } else if (IS_Bioseq_set(sep)) {
   bssp = (BioseqSetPtr)sep->data.ptrvalue;
   sanp = bssp->annot;
 } else return;
 

 while (sanp != NULL) {
   if (sanp->type == 2) {/*seqalign*/
     salp=(SeqAlignPtr)sanp->data;
     sadp->index++;
     sadp->names=Ing_AddStringToSad(sadp->names, sanp, sadp->index);
     if (!(sadp->vnp)) 
       sadp->vnp=ValNodeAddPointer(NULL,OBJ_SEQALIGN,(Pointer)salp);
     else 
       ValNodeAddPointer(&(sadp->vnp),OBJ_SEQALIGN,(Pointer)salp);
   }
   sanp=sanp->next;
 }
}


static SegmenT Ing_GetCurrentSegment(IngPopFeatPtr  pfp)
{
  SegmenT CurrentSeg=NULL;

  if (!(pfp->nPrims%FEATURE_SEGMENT_MAXSIZE))
    {
      if (!(pfp->nLevel2%LEVEL2_MAXNUM)){
        if (!(pfp->nLevel1%LEVEL1_MAXNUM)){
          pfp->nLevel1++;
          pfp->seg1=CreateSegment(pfp->pictMain, pfp->nLevel1, 0);
        }
        pfp->nLevel2++;
        pfp->seg2=CreateSegment(pfp->seg1, pfp->nLevel2, 0);
      }
      pfp->nLevel3++;
      CurrentSeg=CreateSegment(pfp->seg2, pfp->nLevel3, 0);
      pfp->CurrentSeg=CurrentSeg;
/*       pfp->nSegs=1; */
    }
  else 
    CurrentSeg=pfp->CurrentSeg;

  return CurrentSeg;
}

/*******************************************************************************

  Function : Ing_BigEncodeIdxFeat()
 
  Purpose :  index data for collision detection

*******************************************************************************/
extern Uint8 Ing_BigEncodeIdxFeat (Uint4 val1,Uint4 val2)
{
Uint4 index_g[2];
	
	index_g[0]=val1;
	index_g[1]=val2;
	
	return *((Uint8Ptr) index_g);
	
}

/*******************************************************************************

  Function : Ing_BigDecodeIdxFeat()
  
  Purpose : retrieve data from index

*******************************************************************************/
extern void  Ing_BigDecodeIdxFeat (Uint8 index_g, Uint4Ptr val1, Uint4Ptr val2)
{
Uint4Ptr  index_g2;

	index_g2 = (Uint4Ptr) (&index_g);
	if (val1) *val1 = (Uint4) index_g2 [0];
	if (val2) *val2 = (Uint4) index_g2 [1];
}

/*******************************************************************************

  Function : Ing_PopOverviewPage(), Ing_PopOverviewRuler(), Ing_AddToOverViewPage(), Ing_AddToOverviewPict(), Ing_AddAlignsToOverviewPage(), Ing_AddOneAlignToOverviewPage()
 
  Purpose :  Top (Overview) viewer drawing functions

*******************************************************************************/
extern void Ing_PopOverviewRuler(IngGenomeViewerPtr igvp, BioseqPtr bsp)
{
  IngExploreSegs     gpn;
  Int2               nSegments=0;

  gpn.seg=igvp->pictRuler1;
  gpn.viewer=igvp->vRuler1;
  gpn.idx=1;
  gpn.GrData=&igvp->GrData;
  gpn.left=igvp->from;
  gpn.right=igvp->to;
  gpn.seqbuf=igvp->seqbuf;
  gpn.bShowGC = FALSE;

  nSegments=SeqMgrExploreSegments(bsp, (Pointer)&gpn, Ing_ExploreSegments);
  /*this bioseq doesn't have any segments*/
  if (nSegments==0){
    AddAttribute(igvp->pictRuler1, COLOR_ATT, BLACK_COLOR, 0, 0, 0, 0);
    AddRectangle(igvp->pictRuler1,0,igvp->GrData.SegBoxHeight/2 ,igvp->SeqLength,-igvp->GrData.SegBoxHeight/2,NO_ARROW,TRUE,  2);
    igvp->bSegmented=FALSE;
  }
  else{
    igvp->bSegmented=TRUE;
  }
  Ing_AddRuler(igvp->pictRuler1, 2*igvp->GrData.SegBoxHeight, igvp->from, igvp->to, igvp->maxScaleX, 1, TRUE);

}


static void Ing_AddOneAlignToOverviewPage(SegmenT seg, SeqAlignPtr sap, Int4 row, Int4 start, Int4 stop)
{
  BioseqPtr bsp=NULL;
  Uint1     strand;
  Int4      j;
  PrimitivE prim;

    strand=AlnMgrGetNthStrand(sap, 2);
  if (strand == Seq_strand_minus)
    j = -6;
  else
    j = 6;
    prim = AddRectangle(seg, start, (row)*(-20), stop, (row)*(-20)+j, NO_ARROW, FALSE,0);
    SetPrimitiveIDs(prim, sap->idx.entityID, sap->idx.itemID, OBJ_SEQALIGN, OBJ_SEQALIGN);

}


static Uint4 Ing_AddAlignsToOverviewPage(SegmenT pictTop, SegmenT CurrentSeg, Uint2 nPrims, Uint2 nSegs, SeqEntryPtr sep, SeqAlignPtr sap, Uint4 rowNum, Int4 left, Int4 right, CharPtr name, Uint1Ptr Clr)
{
  SeqAlignPtr    salp=NULL;
  enumPrimAddOrder oldOrder;
  Boolean        NewLine=FALSE;
  Int4           start;
  Int4           stop;
  Int4           len;
  Int4           offset;
  Char           buf[41]={""};
  Char           label[41]={""};
  


  len=right-left+1;
  offset=left;
  if (sap->segtype==SAS_DISC){ 
    salp=(SeqAlignPtr)sap->segs; 
  } 
  else {
    salp=sap; 
  } 

  AddAttribute(CurrentSeg, COLOR_ATT, BLACK_COLOR, 0, 0, 0, 0);
  if (StringHasNoText(name))
    sprintf(label, "SeqAlign");
  else {
    sprintf(label,  "%s", name);
  }
  AddLabel(CurrentSeg, -15, (rowNum)*(-20)+6, label, SMALL_TEXT, 0, UPPER_RIGHT, 0);
  AddAttribute(CurrentSeg, COLOR_ATT, (Clr?Clr:BLACK_COLOR), 0, 0, 0, 0);
  AddLine(CurrentSeg, 0, (rowNum)*(-20), len, (rowNum)*(-20), 0, 0);
  oldOrder=ChangeAddPrimOrder(ADD_TO_HEAD);
  while(salp){
    AlnMgrGetNthSeqRangeInSA(salp, 1, &start, &stop);
    if (start>right || stop<left){
      salp=salp->next;
      continue;
    }
    start=MAX(left, start);
    stop=MIN(right, stop);
    if (!(nPrims%FEATURE_SEGMENT_MAXSIZE))
      {
        CurrentSeg=CreateSegment(pictTop, nSegs, 0);
        nSegs++;
        AddAttribute(CurrentSeg, COLOR_ATT, 
                     (Clr?Clr:BLACK_COLOR), 0, 0, 0, 0);
      }
    Ing_AddOneAlignToOverviewPage(CurrentSeg, salp, (Int4)rowNum, start-offset, stop-offset);
    nPrims++;
    NewLine=TRUE;
    salp=salp->next; 
  }
  ChangeAddPrimOrder(oldOrder);
  return rowNum;
}


static void Ing_AddToOverviewPict(SegmenT seg, SeqMgrFeatContextPtr sfc, Int4 start, Int4 stop, Int4 base)
{
  Int4  j;
  PrimitivE  prim;
  Int4  left;
  Int4  right;
  Int4  offset;
    
  if (seg == NULL || sfc == NULL)
     return;
  if (sfc->strand == Seq_strand_minus)
    j = -6;
  else
    j = 6;
  
  offset=start;
  left=MAX(start, sfc->left);
  right=MIN(stop, sfc->right);
  prim = AddRectangle(seg, left-offset, base, right-offset, base+j, NO_ARROW, FALSE,0);
  SetPrimitiveIDs(prim, sfc->entityID, sfc->itemID, OBJ_SEQFEAT, OBJ_SEQFEAT);

}


static Boolean LIBCALLBACK Ing_AddToOverviewPage(SeqFeatPtr sfp, SeqMgrFeatContextPtr sfc)
{
  ValNodePtr vnp=NULL;
  Int4       leftmargin;
  Int4       rightmargin;
  Uint4      rowNum, featdeftype;
  Boolean    bPop=FALSE;
  Boolean    Visible=FALSE;
  Uint8      idx;
  IngPopFeatPtr  pfp;
  Char       str[60]={""};
  Int4       scale;
  Int4       scaledLen;
  Int4       feat_len;
  Int4       pict_len;
  Int4       offset;
  
  pfp=(IngPopFeatPtr)sfc->userdata;
  if (!pfp) return (FALSE);
  leftmargin=pfp->left;
  offset=pfp->left;
  rightmargin=pfp->right;
  pict_len=rightmargin-leftmargin+1;
  feat_len=sfc->right-sfc->left;
  scale=pfp->scaleX;
   if (scale>COMPRESS_SCALE){ 
     scaledLen=feat_len/scale; 
     if (scaledLen<2  &&  IngfeatDefTrack2[sfc->featdeftype]) 
       return (TRUE); 
     else 
       IngfeatDefTrack2[sfc->featdeftype]=TRUE; 
   } 
  
  if (!(pfp->nPrims%FEATURE_SEGMENT_MAXSIZE)) 
    { 
      pfp->CurrentSeg=CreateSegment(pfp->pictMain, pfp->nLevel1, 0);
      pfp->nLevel1++;
    } 
  
  AddAttribute(pfp->CurrentSeg, COLOR_ATT, pfp->pClr[sfc->featdeftype] , 0, 0, 0, 0);
  if (!pfp->PopRowsList){
    rowNum = 0 ;
    idx=Ing_BigEncodeIdxFeat(rowNum, sfc->featdeftype);
    pfp->PopRowsList = ValNodeAddBigInt (NULL, 5, idx);
    bPop=FALSE;
  }
  else{
    vnp = pfp->PopRowsList;
    while(vnp){
      Ing_BigDecodeIdxFeat((Uint8)vnp->data.bigintvalue, &rowNum, &featdeftype);
      if (sfc->featdeftype == featdeftype){
        bPop = TRUE;
        break;
        }
      vnp = vnp->next;
    }
  }
  
  if (bPop==FALSE){
    rowNum = pfp->nPopRows;
    idx = Ing_BigEncodeIdxFeat(rowNum,sfc->featdeftype);
    ValNodeAddBigInt (&pfp->PopRowsList, 5, idx);
    pfp->nPopRows++;
    FeatDefLabel(sfp, str, sizeof(str)-1, OM_LABEL_TYPE);
    AddLabel(pfp->CurrentSeg, -15, (rowNum)*(-20)+6, str, SMALL_TEXT, 0, UPPER_RIGHT, 0);
    AddLine(pfp->CurrentSeg, 0, (rowNum)*(-20), pict_len, (rowNum)*(-20), 0, 0);
    pfp->nPopRows++;
  }
  Ing_AddToOverviewPict(pfp->CurrentSeg, sfc, leftmargin, rightmargin, rowNum*(-20));
  
  pfp->nPrims++;
  return (TRUE);
}

extern Uint4 Ing_PopOverviewPage(BioseqPtr bsp, SegmenT pictTop, Int4 left, Int4 right, Uint1Ptr PNTR pClr, Int4 maxScaleX)
{
  IngPopFeat   pf;
  Uint2        entityID;
  Int4         i;
  SeqEntryPtr  sep;
  ValNodePtr   vnp=NULL;
  SeqAnnotPtr  sanp=NULL;
  SeqAlignPtr  sap=NULL;
  SeqLocPtr    slp=NULL, whole_slp=NULL;
  Uint4        nPopRows=0;
  IngSeqAnnotData  sad;
  CharPtr      curr_name=NULL;
  SegmenT      CurrentSeg;
  Char         buf[50]={""};
  Int4         margin, len;
  Int4         start=0;
  ValNodePtr   slp_list=NULL, slp_head=NULL;
  Uint1Ptr     Clr=NULL;
 

  whole_slp=SeqLocIntNew(left, right, Seq_strand_plus, (SeqIdPtr)bsp->id);

  if (maxScaleX>COMPRESS_SCALE){
    start=left;
    len=right;
    while (start < len){
      slp=SeqLocIntNew(start, start+maxScaleX, Seq_strand_plus, (SeqIdPtr)bsp->id);
      ValNodeAddPointer(&slp_list, OBJ_BIOSEQ, (Pointer)slp);
      start+=maxScaleX;
    }
    slp=SeqLocIntNew(start, len-1, Seq_strand_plus, (SeqIdPtr)bsp->id);
    ValNodeAddPointer(&slp_list, OBJ_BIOSEQ, (Pointer)slp);
    MemSet(&pf,0,sizeof(pf));
    pf.pictMain=pictTop;
    pf.pClr=pClr;
    pf.scaleX=maxScaleX;
    pf.nPrims=0;
    pf.nLevel1=1;
    pf.left=left;
    pf.right=right;
    slp_head=slp_list;
    while(slp_list){
      slp=slp_list->data.ptrvalue;
      MemSet((Pointer)IngfeatDefTrack2, 0, sizeof(IngfeatDefTrack2));
      SeqMgrExploreFeatures (bsp, (Pointer) &pf,Ing_AddToOverviewPage, slp, NULL,IngfeatDefFilter);
      slp_list=slp_list->next;
      SeqLocFree(slp);
    }
    ValNodeFree(slp_head); 
  }
  else {
    MemSet(&pf,0,sizeof(pf));
    pf.pictMain=pictTop;
    pf.pClr=pClr;
    pf.scaleX=maxScaleX;
    pf.nPrims=0;
    pf.nLevel1=1;
    pf.left=left;
    pf.right=right;
    SeqMgrExploreFeatures (bsp, (Pointer) &pf,Ing_AddToOverviewPage, whole_slp, NULL,IngfeatDefFilter);
  }

  nPopRows=pf.nPopRows;
  ValNodeFree(pf.PopRowsList);
  pf.nLevel1++;
  CurrentSeg=CreateSegment(pictTop, pf.nLevel1, 0);
  if (IngfeatDefFilter[0]){ 
     entityID=ObjMgrGetEntityIDForPointer(bsp); 
     sep = GetTopSeqEntryForEntityID(entityID); 
     sad.slp=whole_slp;
     sad.vnp=NULL;
     sad.names=NULL;
     sad.index=0;
     SeqEntryExplore(sep,(Pointer)&sad, Ing_SearchAli);
     vnp=sad.vnp;
     if (!vnp) goto end;
     IngfeatDefTrack[0]=2;
     i=0;
     len=right-left+1;
     margin=10*maxScaleX;
     while (vnp){
       sap=(SeqAlignPtr)vnp->data.ptrvalue;
       if (sap->segtype!=SAS_DISC && StringCmp(sad.names[i], "Spidey")==0)
         goto skipindex;
       
       if (sap->saip != NULL)
         AMAlignIndexFree(sap);
       sap->saip=NULL;
       AlnMgrIndexLite(sap);
       AlnMgrSortAlnSetByNthRowPos(sap, 1);
        
     skipindex:
       AssignIDsInEntity(entityID, OBJ_SEQALIGN, (Pointer)sap);
       Clr=Ing_GetAlignColor(sad.names[i], pClr);
       nPopRows=Ing_AddAlignsToOverviewPage(pictTop, CurrentSeg, 0, pf.nLevel1, sep, sap, nPopRows, left, right, sad.names[i], Clr);
       vnp=vnp->next;
       nPopRows+=2;
       i++;
     }
     ValNodeFree(sad.vnp);
   } 

 end:
  SeqLocFree(whole_slp);
  return (nPopRows);
}




/*******************************************************************************

  Function : Ing_PopFeaturesPage(), Ing_AddToFeaturesPage(), Ing_AddToFeaturePict(), Ing_AddAlignsToFeaturesPage(), Ing_AddOneAlignToFeaturesPage()
 
  Purpose :  Bottom (Features) viewer drawing functions

*******************************************************************************/

static void Ing_AddOneAlignToFeaturesPage(SegmenT seg, SeqAlignPtr sap, Int4 row, Int4 start, Int4 stop, Int4 offset, Int4 scale, Uint2 alignID)
{
   Uint1     strand;
  Int4      i=0;
  PrimitivE prim;
  Int4      numivals;
  AlnMsgPtr amp=NULL;
  SegmenT   subSeg;
  Boolean   more;
  Uint1     ARROW=NO_ARROW;

 
  strand=AlnMgrGetNthStrand(sap, 2);
  numivals= AlnMgrGetNumSegments(sap);
  if (numivals>1 && scale<ARROWS_SCALE){
    subSeg=CreateSegment(seg, alignID, 0);
    amp = AlnMsgNew();
    amp->to_m=-1; 
    amp->row_num = 1;
    i=1;
    while((Boolean) (more = AlnMgrGetNextAlnBit(sap, amp))){
        if (amp->gap==0){
          if (strand==Seq_strand_minus && i==1)
            ARROW=LEFT_ARROW;
          else if (strand==Seq_strand_plus && i==numivals)
            ARROW=RIGHT_ARROW;
          else
            ARROW=NO_ARROW;
          if (amp->from_b>stop || amp->to_b<start){
            i++;
            continue;
          } 
          prim = AddRectangle(subSeg, MAX(start, amp->from_b)-offset, (row)*(-20), MIN(stop, amp->to_b)-offset, (row)*(-20)-Ing_FEAT_HEIGHT, ARROW, TRUE,0);
          SetPrimitiveIDs(prim, sap->idx.entityID, sap->idx.itemID, OBJ_SEQALIGN, 0);
        }
        i++;
    }
    AlnMsgFree(amp);
    AddLine(subSeg, start-offset, (row)*(-20), stop-offset, (row)*(-20), 0, 0);
  }
  else {
    if (scale<ARROWS_SCALE){
      prim = AddRectangle(seg, start-offset, (row)*(-20), stop-offset, (row)*(-20)-Ing_FEAT_HEIGHT, (strand==Seq_strand_minus?LEFT_ARROW:RIGHT_ARROW), TRUE, 0);
    }
    else {
      prim= AddRectangle(seg, start-offset, (row)*(-20), stop-offset, (row)*(-20)-2, NO_ARROW, TRUE, 0);
    }
    SetPrimitiveIDs(prim, sap->idx.entityID, sap->idx.itemID, OBJ_SEQALIGN, alignID);
  }
  
}

static Uint4 Ing_NextExon(Uint4 exon, Uint1 strand)
{
  if (strand==Seq_strand_plus)
    exon++;
  else
    exon--;

  return exon;
}

static Uint4 Ing_AddAlignsToFeaturesPage(SegmenT pictBottom, IngPopFeatPtr pfp, SeqEntryPtr sep, SeqAlignPtr sap, Uint4 rowNum, Int4 left, Int4 right, Uint1Ptr Clr, Uint4 scaleX, CharPtr name, Uint2 alignID, Boolean bGenes)
{
  Char             buf[41]={""};
  Char             label[41]={""};
  Uint4            top, bottom;
  enumPrimAddOrder oldOrder;
  ValNodePtr       vnp=NULL, tmp_vnp;
  SeqIdPtr         sip=NULL;
  Uint8            idx;
  Boolean          NewLine=FALSE;
  Uint4            rowCnt=0;
  Uint4            ColStart;
  Int4             start, stop;
  Int4             i;
  AMAlignIndexPtr  amaip=NULL;
  Int4             len;
  Int4             len_diff;
  Int4             bExonlabels=FALSE;
  Int4             label_len;
  Int4             feat_len;
  Boolean          labels=TRUE;
  Boolean          bDraw = FALSE;
  Int4             offset;
  Uint4            exon=0;
  Uint1            strand=0;
  SegmenT          CurrentSeg=NULL;


  amaip = (AMAlignIndexPtr)(sap->saip);
  if (!amaip) return rowNum; 
  len=right-left+1;
  offset=left;

  if ((StringCmp(name, "Spidey")==0) && !bGenes){
    bExonlabels=TRUE;
    strand=AlnMgrGetNthStrand(amaip->saps[0], 2);
    if (strand==Seq_strand_minus)
      exon=amaip->numsaps;
    else  
      exon=1;
  }

  if (scaleX > LABELS_SCALE)
    labels=FALSE;
  oldOrder=ChangeAddPrimOrder(ADD_TO_HEAD);
  top=(rowNum)*(-20);
  rowNum++;
  rowCnt=rowNum;

  for(i=0; i<amaip->numsaps; i++){
    AlnMgrGetNthSeqRangeInSA(amaip->saps[i], 1, &start, &stop);
    if (start>right || stop<left)
      continue;

    CurrentSeg = Ing_GetCurrentSegment(pfp);
    if (!CurrentSeg) return (FALSE);
    start=MAX(start, left);
    stop=MIN(stop, right);
    if (labels && pfp->showLabels) {
      SetSmallFont();
      if (bExonlabels){
        sprintf(label, "Exon %ld", exon);
        exon = Ing_NextExon(exon, strand);
      }
      else {
        sip=AlnMgrGetNthSeqIdPtr(amaip->saps[i], 2);
        SeqIdWrite(sip, buf, PRINTID_REPORT, 41);
        sprintf(label,"gi|%s", buf);
      }
      label_len=(StringWidth(label)+2)*scaleX;
      feat_len=stop-start-1;
      if (feat_len<label_len)
        len_diff= ((label_len-feat_len)/2)+(1*scaleX);
      else
        len_diff=1*scaleX;
    }
    else{
      len_diff=1*scaleX;
    }

    if (!vnp){
      ColStart=stop+len_diff;
      idx=Ing_BigEncodeIdxFeat(rowCnt, ColStart);
      vnp=ValNodeAddBigInt(NULL, 5, idx);
      NewLine=FALSE;
    }
    else{
      tmp_vnp=vnp;
      while(tmp_vnp){
        Ing_BigDecodeIdxFeat((Uint8)tmp_vnp->data.bigintvalue, &rowCnt, &ColStart);
        if (ColStart < MAX(0, start-len_diff)){
          ColStart=stop+len_diff;
          tmp_vnp->data.bigintvalue=Ing_BigEncodeIdxFeat(rowCnt,ColStart);
          NewLine=FALSE;
          break;
        }
        tmp_vnp=tmp_vnp->next;
      }
    }
    if (NewLine){
      rowNum++;
      rowCnt=rowNum;
      ColStart=stop+len_diff;
      idx=Ing_BigEncodeIdxFeat(rowNum, ColStart);
      ValNodeAddBigInt(&vnp, 5, idx);
    }
    if (labels && pfp->showLabels){
      AddAttribute(CurrentSeg, COLOR_ATT, BLACK_COLOR, 0, 0, 0, 0);
      AddLabel (CurrentSeg, start - offset + (stop-start)/2, rowCnt*(-20), label, SMALL_TEXT, 0,  UPPER_CENTER, 0); 
    }
    AddAttribute(CurrentSeg, COLOR_ATT, Clr, 0, 0, 0, 0);
    Ing_AddOneAlignToFeaturesPage(CurrentSeg, amaip->saps[i], (Int4)rowCnt, start, stop, offset, scaleX, alignID);
    pfp->nPrims++;
    NewLine = TRUE;
    bDraw = TRUE;
  }

  if (bDraw){
    CurrentSeg = Ing_GetCurrentSegment(pfp);
    if (!CurrentSeg) return (FALSE);
    
    AddAttribute(CurrentSeg, COLOR_ATT, BLACK_COLOR, 0, 0, 0, 0);
    if (StringHasNoText(name))
      sprintf(label, "SeqAlign: ");
    else {
      sprintf(label,  "%s: ", name);
    }
    sip=AlnMgrGetNthSeqIdPtr(amaip->saps[0], 2);
    SeqIdWrite(sip, buf, PRINTID_REPORT, 41);
    StringCat(label, buf);
    AddLabel(CurrentSeg, -15, top+6, label, SMALL_TEXT, 0, UPPER_RIGHT, 0);
    
    rowNum++;
    bottom=(rowNum)*(-20);
    AddRectangle(CurrentSeg, (-10)*scaleX, top, len+(10*scaleX), bottom, NO_ARROW, FALSE,0);
    ChangeAddPrimOrder(oldOrder);
  }
  else {
    rowNum-=3;
  }

  return rowNum;
}
                                      
static void Ing_AddToFeaturePict(SegmenT Seg, SeqMgrFeatContextPtr sfc, Int4 yBase, Uint4Ptr subSegID, Int4 start, Int4 stop, Int4 scale)
{
  PrimitivE  prim;
  Int4       left_end, right_end;
  Uint1      strand;
  Int4       fhalf;
  Int4       numivals, max, i;
  Int4       left, right;
  Int4Ptr    ivals;
  SegmenT    subSeg;
  SeqFeatPtr sfp;
  Boolean    arrows=TRUE;
  Int4       offset;

  if (Seg == NULL || sfc == NULL)
     return;
  
  if (scale>ARROWS_SCALE)
    arrows=FALSE;
  fhalf=Ing_FEAT_HEIGHT/2;
  strand=sfc->strand;
  offset=start;
  left=MAX(start, sfc->left);
  right=MIN(stop, sfc->right);
  
  ivals=sfc->ivals;
  numivals=sfc->numivals;
  
  if (numivals>1 && arrows){
    subSeg=CreateSegment(Seg, *subSegID, 0);
    (*subSegID)++;
    sfp=(SeqFeatPtr)sfc->sfp;
    /* draw a line thought segmented feature */
    left_end=(strand == Seq_strand_minus ? ivals[2*numivals-1] : ivals[0]);
    right_end=(strand == Seq_strand_minus ? ivals[0] : ivals[2*numivals-1]);
    if (!(right_end < offset)){
      left_end=MAX(left_end, offset);
      prim=AddLine(subSeg, left_end-offset, yBase-fhalf, right_end-offset, yBase-fhalf, FALSE, 0);
      SetPrimitiveIDs(prim, sfc->entityID, sfc->itemID, sfp->idx.itemtype, 0);
    }
      /* far left ivals */
    left_end=(strand == Seq_strand_minus ? ivals[2*numivals-1] : ivals[0]); 
    right_end=(strand == Seq_strand_minus ? ivals[2*numivals-2] : ivals[1]);
    if (!(right_end < offset)){
      left_end=MAX(left_end, offset);
      prim = AddRectangle(subSeg, left_end-offset, yBase, right_end-offset, yBase-Ing_FEAT_HEIGHT, (strand == Seq_strand_minus ? LEFT_ARROW : NO_ARROW), TRUE, 0);  
      SetPrimitiveIDs(prim, sfc->entityID, sfc->itemID, sfp->idx.itemtype, 0);
    }
    /* middle ivals */
    if (numivals>2){
      if (strand==Seq_strand_minus){
        for (i=2*numivals-4;i>1;i-=2){
          left_end=MIN(ivals[i], ivals[i+1]);
          right_end=MAX(ivals[i], ivals[i+1]);
          if ((!(left_end<offset && right_end<offset))&&(!(left_end>stop && right_end>stop))){
            prim = AddRectangle(subSeg,MAX(start, left_end)-offset, yBase, MIN(stop, right_end)-offset,yBase-Ing_FEAT_HEIGHT,NO_ARROW,TRUE,0);
            SetPrimitiveIDs(prim, sfc->entityID, sfc->itemID, sfp->idx.itemtype, 0);
          }
          
        }
      }
      else{
        max=2*numivals-2;
        for (i=2;i<max;i+=2){
          left_end=MIN(ivals[i], ivals[i+1]);
          right_end=MAX(ivals[i], ivals[i+1]);
          if ((!(left_end<offset && right_end<offset)) && 
              (!(left_end>stop && right_end>stop))){
          prim = AddRectangle(subSeg, MAX(start, left_end)-offset,yBase, MIN(stop, right_end)-offset, yBase-Ing_FEAT_HEIGHT, NO_ARROW, TRUE, 0);
          SetPrimitiveIDs(prim, sfc->entityID, sfc->itemID, sfp->idx.itemtype, 0);
          }
        }
      }               
    }
    /* far right ivals*/
    left_end=(strand == Seq_strand_minus ? ivals[0] : ivals[2*numivals-2]);
    right_end=(strand == Seq_strand_minus ?ivals[1] : ivals[2*numivals-1]);
    if (!(left_end>right && right_end>right)){
      prim = AddRectangle(subSeg, MAX(start, left_end)-offset, yBase, MIN(stop, right_end)-offset, yBase-Ing_FEAT_HEIGHT,  (strand == Seq_strand_minus ? NO_ARROW:RIGHT_ARROW), TRUE, 0);
      SetPrimitiveIDs(prim, sfc->entityID, sfc->itemID, sfp->idx.itemtype, 0);
    }
  }
  else {
    if (arrows){
      prim = AddRectangle(Seg, left-offset, yBase, right-offset, yBase-Ing_FEAT_HEIGHT, (strand == Seq_strand_minus ?LEFT_ARROW:RIGHT_ARROW), TRUE, 0);
    }
    else {
      prim= AddRectangle(Seg, left-offset, yBase, right-offset, yBase-2, NO_ARROW, TRUE, 0);
    }

    SetPrimitiveIDs(prim, sfc->entityID, sfc->itemID, OBJ_SEQFEAT, OBJ_SEQFEAT);
  }
}


static Boolean LIBCALLBACK Ing_AddToFeaturesPage(SeqFeatPtr sfp, SeqMgrFeatContextPtr sfc)
{
  ValNodePtr vnp=NULL;
  Int4       leftmargin;
  Int4       rightmargin;
  Int4       left;
  Int4       right;
  Uint4      rowNum;
  Boolean    bPop=FALSE;
  Uint8      idx;
  Uint4       ColStart;
  IngPopFeatPtr  pfp;
  Char       str[60]={""};
  Int4       scaleX;
  Int4       label_len, len_diff;
  Int4       feat_len;
  Boolean    labels=TRUE;
  SegmenT    CurrentSeg=NULL;

 
  pfp=(IngPopFeatPtr)sfc->userdata;
  if (!pfp) return (FALSE);

  scaleX=pfp->scaleX;
  leftmargin=pfp->left;
  rightmargin=pfp->right;
  left=MAX(leftmargin, sfc->left);
  right=MIN(rightmargin, sfc->right);
  feat_len=right-left;
  if (scaleX>LABELS_SCALE){ 
    labels=FALSE; 
   }
  if (labels && pfp->showLabels){
    SetSmallFont();
    FeatDefLabel (sfp, str, sizeof (str) - 1, OM_LABEL_CONTENT);
    label_len=(StringWidth(str))*scaleX;
    if (label_len > feat_len)
      len_diff= ((label_len-feat_len)/2)+(1*scaleX);
    else
      len_diff=1*scaleX;
  }
  else{
    len_diff=1*scaleX;
  }
  IngfeatDefTrack[sfc->featdeftype]=2;

  CurrentSeg=Ing_GetCurrentSegment(pfp);
  if (!CurrentSeg) return (FALSE);

  AddAttribute(pfp->CurrentSeg, COLOR_ATT, pfp->pClr[sfc->featdeftype] , 0, 0, 0, 0);
  if (!pfp->PopRowsList){
    rowNum = pfp->nPopRows;
    ColStart=sfc->right+len_diff;
    idx=Ing_BigEncodeIdxFeat(rowNum, ColStart);
    pfp->PopRowsList = ValNodeAddBigInt (NULL, 5, idx);
    bPop=TRUE;
    pfp->nPopRows++;
  }
  else{
    vnp=pfp->PopRowsList;
    while(vnp){
      Ing_BigDecodeIdxFeat((Uint8)vnp->data.bigintvalue,&rowNum,&ColStart);
      if (ColStart < MAX(0, left-len_diff)){
        ColStart=sfc->right+len_diff;
        vnp->data.bigintvalue=Ing_BigEncodeIdxFeat(rowNum,ColStart);
        bPop=TRUE;
        break;
        }
      vnp=vnp->next;
    }
  }
  if (bPop==FALSE){
    rowNum = pfp->nPopRows;
    ColStart=sfc->right+len_diff;
    idx = Ing_BigEncodeIdxFeat(rowNum, ColStart);
    ValNodeAddBigInt (&pfp->PopRowsList, 5, idx);
    pfp->nPopRows++;
  }
  if (labels && pfp->showLabels)
    AddLabel (pfp->CurrentSeg, left -leftmargin + (right-left)/2, rowNum*(-20), str, SMALL_TEXT, 0,  UPPER_CENTER, 0); 

  Ing_AddToFeaturePict(pfp->CurrentSeg, sfc, rowNum*(-20), &pfp->nSegs, leftmargin, rightmargin, scaleX);
  
  pfp->nPrims++;
  return (TRUE);
}


extern Uint4 Ing_PopFeaturesPage(BioseqPtr bsp, SegmenT pictBottom, Int4 left, Int4 right, Uint1Ptr PNTR pClr, Int4 scaleX, Int4 start_row, Boolean Labels, Boolean bGenes)
{
  IngPopFeat   pf;
  Uint2        entityID;
  SeqEntryPtr  sep;
  ValNodePtr   vnp=NULL;
  SeqAlignPtr  sap=NULL;
  SeqLocPtr    slp=NULL;
  Uint4        nPopRows=0;
  Int4         i;
  Int4         alignID;
  IngSeqAnnotData  sad;
  SegmenT      CurrentSeg=NULL;
  Uint1Ptr     Clr=NULL;


	memset(&pf,0,sizeof(pf));
	pf.pictMain=pictBottom;
	pf.pClr=pClr;
	pf.scaleX=scaleX;
   pf.nPrims=0;
   pf.nLevel1=0;
   pf.nLevel2=0;
   pf.nLevel3=0;
   pf.left=left;
   pf.right=right;
   pf.showLabels=Labels;
   pf.nPopRows=start_row;
   slp=SeqLocIntNew(left, right, Seq_strand_plus, (SeqIdPtr)bsp->id);
	SeqMgrExploreFeatures (bsp, (Pointer) &pf,Ing_AddToFeaturesPage, slp, NULL, IngfeatDefFilter);
   /* get annotated alignments */
   nPopRows=pf.nPopRows;
   ValNodeFree(pf.PopRowsList); 
   CurrentSeg=CreateSegment(pictBottom,pf.nLevel1+1, 0);

   if (IngfeatDefFilter[0]){
     nPopRows+=2;
     entityID=ObjMgrGetEntityIDForPointer(bsp);
     sep = GetTopSeqEntryForEntityID(entityID); 
     sad.slp=slp;
     sad.vnp=NULL;
     sad.names=NULL;
     sad.index=0;
     SeqEntryExplore(sep,(Pointer)&sad, Ing_SearchAli);
     vnp=sad.vnp;
     if (!vnp) goto end;
     IngfeatDefTrack[0]=2;
     i=0;
     alignID=1;
     /* alignID is same as segID in Spidey Report */
     while(vnp){
       sap=(SeqAlignPtr)vnp->data.ptrvalue;
       if (sap->saip != NULL)
         AMAlignIndexFree(sap);
       sap->saip=NULL;
       AlnMgrIndexLite(sap);
       AlnMgrSortAlnSetByNthRowPos(sap, 1);
       AssignIDsInEntity(entityID, OBJ_SEQALIGN, (Pointer)sap);
       Clr=Ing_GetAlignColor(sad.names[i], pClr);
       nPopRows=Ing_AddAlignsToFeaturesPage(pictBottom, &pf, sep, sap, nPopRows, left, right, Clr, scaleX, sad.names[i], alignID, bGenes);
       nPopRows+=2;
       i++;
       alignID++;
       vnp=vnp->next;
     }
     ValNodeFree(sad.vnp);
   }

 end:
   SeqLocFree(slp);
   return (nPopRows);
}





/*******************************************************************************

  Function : Ing_GetValue()
 
  Purpose :  returns int value given TexT

*******************************************************************************/

extern Int4 Ing_GetValue (TexT t)
{
  Char str[20];
  Int4 val;

  GetTitle (t,  str,  sizeof(str));
  if (StringHasNoText(str))
    {
      ErrPostEx (SEV_WARNING, 0, 0, "%s", "missing parameter(s)");
      return 0;
    }

  val=atoi(str);

  return val;
}


/*****************************************************************************

Function: Ing_PutColor(), Ing_FreeColor(), Ing_BuildColorTable(), Ing_FreeColorTable()

Purpose: color functions

*****************************************************************************/
extern Uint1Ptr  Ing_PutColor(Uint1 r, Uint1 g, Uint1 b)
{
  Uint1Ptr Clr;

  Clr = (Uint1Ptr)MemNew(sizeof(Uint1)*3);
  Clr[0]=r;
  Clr[1]=g;
  Clr[2]=b;
  return (Clr);
}

static Uint1Ptr Ing_FreeColor(Uint1Ptr Clr)
{
  MemFree(Clr);
  return NULL;
}

extern Uint1Ptr PNTR Ing_BuildColorTable(void)
{
Uint1Ptr PNTR pClr;


	pClr=(Uint1Ptr PNTR)MemNew(Ing_MAX*sizeof(Uint1Ptr));
	if (!pClr) return(NULL);

	pClr[FEATDEF_GENE]=Ing_PutColor(204,45,61); /*deep red*/
	
	pClr[FEATDEF_precursor_RNA]=Ing_PutColor(128,128,128);/*grey*/
	pClr[FEATDEF_misc_RNA]=Ing_PutColor(128,128,128);/*grey*/
	pClr[FEATDEF_preRNA]=Ing_PutColor(128, 128,128);/*grey*/
	pClr[FEATDEF_mRNA]=Ing_PutColor(0,127,0);/*green*/
	pClr[FEATDEF_tRNA]=Ing_PutColor(224,95,39);/* orange*/
	pClr[FEATDEF_rRNA]=Ing_PutColor(157,34, 28);/*sienna*/
	pClr[FEATDEF_snRNA]=Ing_PutColor(127,0,0);/* burnt sienna */
	pClr[FEATDEF_scRNA]=Ing_PutColor(119,45,0);/*brown*/
	pClr[FEATDEF_otherRNA]=Ing_PutColor(128,0,0);/*grey*/
	pClr[FEATDEF_prim_transcript]=Ing_PutColor(0,74, 0);/*green*/
	pClr[FEATDEF_polyA_signal]=Ing_PutColor(0,174, 0);/*green*/
	pClr[FEATDEF_polyA_site]=Ing_PutColor(0,200, 12);/*green*/
   pClr[FEATDEF_repeat_region]=Ing_PutColor(114,204,0)/* dk lime *//* (0,155,220) */;/*cyan*/

	pClr[FEATDEF_CDS]=Ing_PutColor(235, 150, 235);	/*pink*/
	pClr[FEATDEF_exon]=Ing_PutColor(192,50,150); /*fuschia */
	pClr[FEATDEF_intron]=Ing_PutColor(255,170,170);/* pale pink*/

	pClr[FEATDEF_PROT]=Ing_PutColor(0,0,128);		/*v dk blue*/
	pClr[FEATDEF_mat_peptide]=Ing_PutColor(64, 128, 192);/*pale blue*/
	pClr[FEATDEF_sig_peptide]=Ing_PutColor(204,0,61);/* dk red */
	pClr[FEATDEF_transit_peptide]=Ing_PutColor(224,224,0);/* lime */
	pClr[FEATDEF_preprotein]=Ing_PutColor(0, 19, 127);/*dk blue */
	pClr[FEATDEF_mat_peptide_aa]=Ing_PutColor(64, 128, 192);/*pale blue*/
	pClr[FEATDEF_sig_peptide_aa]=Ing_PutColor(63, 52, 90);/* indigo */
	pClr[FEATDEF_transit_peptide_aa]=Ing_PutColor(224,224,0);/* lime */

	pClr[FEATDEF_misc_feature]= Ing_PutColor(210, 154, 14 );/* curry */


	pClr[FEATDEF_SITE]=Ing_PutColor(255,0,0);		/*red*/

	pClr[FEATDEF_REGION]=Ing_PutColor(210,154,14);		/*orange*/
	pClr[FEATDEF_mutation]=Ing_PutColor(210,154,14);
	pClr[FEATDEF_variation]=Ing_PutColor(210,154,14);

	pClr[FEATDEF_PSEC_STR]=Ing_PutColor(104,201,220); /*cyan*/
   pClr[FEATDEF_STS]=Ing_PutColor(203, 52, 220); /*pale purple*/
 
	pClr[FEATDEF_HET]=Ing_PutColor(128,128,0);	/*yellow*/

	pClr[FEATDEF_BOND]=Ing_PutColor(255,92,255);	/*pink*/
   pClr[FEATDEF_unsure]=Ing_PutColor(104, 0, 40);/* dkred */
   pClr[FEATDEF_IMP]=Ing_PutColor(104, 0, 40);/* dkred */
   /* user defined */
	pClr[Ing_SPIDEY]=Ing_PutColor(0,200,112);/*green*/
   pClr[ALIGN_BLAST]=Ing_PutColor(155,145,0);
   pClr[ALIGN_FILE]=Ing_PutColor(145,155,255);
   pClr[ALIGN_ANNOT]=Ing_PutColor(255,92,55);
   pClr[ALIGN_BLAST2SEQ]=Ing_PutColor(25,92,255); 
   pClr[Ing_GENSCAN]=Ing_PutColor(204,45,61)/* (224,224,0) */;
   pClr[Ing_FEAT_TABLE]=Ing_PutColor(0,224,0);
   pClr[Ing_TBDL]=Ing_PutColor(0,224,155);

   /* extras */
   pClr[Ing_black]=Ing_PutColor(0,0,0);
   pClr[Ing_red]=Ing_PutColor(224, 0, 60);
   pClr[Ing_blue]=Ing_PutColor(0, 0, 255);
   pClr[Ing_cyan]=Ing_PutColor(0, 235, 245);
   pClr[Ing_yellow]=Ing_PutColor(235, 245, 0);
   pClr[Ing_grey]=Ing_PutColor(127, 127, 127);
   pClr[Ing_purple]=Ing_PutColor(163, 52, 190);
     

	return(pClr);
}


extern Uint1Ptr PNTR Ing_FreeColorTable(Uint1Ptr PNTR pClr)
{


	if (!pClr) return(NULL);

	Ing_FreeColor(pClr[FEATDEF_GENE]);
	Ing_FreeColor(pClr[FEATDEF_precursor_RNA]);
   Ing_FreeColor(pClr[FEATDEF_misc_RNA]);
   Ing_FreeColor(pClr[FEATDEF_preRNA]);
   Ing_FreeColor(pClr[FEATDEF_mRNA]);
   Ing_FreeColor(pClr[FEATDEF_tRNA]);
   Ing_FreeColor(pClr[FEATDEF_rRNA]);
   Ing_FreeColor(pClr[FEATDEF_snRNA]);
   Ing_FreeColor(pClr[FEATDEF_scRNA]);
   Ing_FreeColor(pClr[FEATDEF_otherRNA]);
   Ing_FreeColor(pClr[FEATDEF_prim_transcript]);
   Ing_FreeColor(pClr[FEATDEF_polyA_signal]);
   Ing_FreeColor(pClr[FEATDEF_polyA_site]);
   Ing_FreeColor(pClr[FEATDEF_repeat_region]);
   Ing_FreeColor(pClr[FEATDEF_CDS]);
   Ing_FreeColor(pClr[FEATDEF_exon]);
   Ing_FreeColor(pClr[FEATDEF_intron]);
   Ing_FreeColor(pClr[FEATDEF_PROT]);
   Ing_FreeColor(pClr[FEATDEF_mat_peptide]);
   Ing_FreeColor(pClr[FEATDEF_sig_peptide]);
   Ing_FreeColor(pClr[FEATDEF_transit_peptide]);
   Ing_FreeColor(pClr[FEATDEF_preprotein]);
   Ing_FreeColor(pClr[FEATDEF_mat_peptide_aa]);
   Ing_FreeColor(pClr[FEATDEF_sig_peptide_aa]);
   Ing_FreeColor(pClr[FEATDEF_transit_peptide_aa]);
   Ing_FreeColor(pClr[FEATDEF_misc_feature]);
   Ing_FreeColor(pClr[FEATDEF_SITE]);
   Ing_FreeColor(pClr[FEATDEF_REGION]);
   Ing_FreeColor(pClr[FEATDEF_mutation]);
   Ing_FreeColor(pClr[FEATDEF_variation]);
   Ing_FreeColor(pClr[FEATDEF_PSEC_STR]);
   Ing_FreeColor(pClr[FEATDEF_STS]);
   Ing_FreeColor(pClr[FEATDEF_HET]);
   Ing_FreeColor(pClr[FEATDEF_BOND]);
   Ing_FreeColor(pClr[FEATDEF_unsure]);
   Ing_FreeColor(pClr[Ing_ORF]);
   Ing_FreeColor(pClr[Ing_SPIDEY]);
   Ing_FreeColor(pClr[ALIGN_BLAST]);
   Ing_FreeColor(pClr[ALIGN_FILE]);
   Ing_FreeColor(pClr[ALIGN_ANNOT]);
   Ing_FreeColor(pClr[ALIGN_BLAST2SEQ]);
   Ing_FreeColor(pClr[Ing_GENSCAN]);
   Ing_FreeColor(pClr[Ing_FEAT_TABLE]);
   Ing_FreeColor(pClr[Ing_TBDL]);
   Ing_FreeColor(pClr[Ing_black]);
   Ing_FreeColor(pClr[Ing_red]);
   Ing_FreeColor(pClr[Ing_blue]);
   Ing_FreeColor(pClr[Ing_cyan]);
   Ing_FreeColor(pClr[Ing_yellow]);
   Ing_FreeColor(pClr[Ing_grey]);
   Ing_FreeColor(pClr[Ing_purple]);
   MemFree(pClr);
   return NULL;
}



/*************************************************

  Function : Ing_ExploreSegments()
  
  Purpose : draw callback for SeqMgrExploreSegments()

**************************************************/
extern Boolean LIBCALLBACK Ing_ExploreSegments(SeqLocPtr slp, 
	SeqMgrSegmentContextPtr context)
{
  IngExploreSegsPtr  gpnp;
  Int4               left,top,right,bottom, left_ex, right_ex;
  Boolean            is_inrange=TRUE; 
  Uint1              orange[3]/* , indigo[3] */;
  Uint1              strand;

   orange[0]= 210;
   orange[1]= 154;
   orange[2]= 14;

	gpnp=(IngExploreSegsPtr)context->userdata;
   if (gpnp->idx)
     gpnp->idx++;
	left=context->from;
	right=context->to-left+context->cumOffset;
	left=context->cumOffset;

   if (!(gpnp->left==0 && gpnp->right==0)){
     left_ex=gpnp->left;
     right_ex=gpnp->right;
     if ((left<left_ex && right<left_ex)||(left>right_ex && right>right_ex)){
       return (TRUE);
     }
     else {
       left=MAX(left, left_ex)-left_ex;
       right=MIN(right, right_ex)-left_ex;
     }
   }
   
   strand=SeqLocStrand(slp);
   if (strand==Seq_strand_minus){
     AddAttribute(gpnp->seg, COLOR_ATT, orange, 0, 0, 0, 0);
   }
   else{
     AddAttribute(gpnp->seg, COLOR_ATT, BLACK_COLOR, 0, 0, 0, 0);
   }

   top=gpnp->GrData->SegBoxHeight;
   if (gpnp->idx % 2){
     top-=gpnp->GrData->SegBoxHeight;
   }
   bottom=top-gpnp->GrData->SegBoxHeight;
   
   Ing_AddGCRect(gpnp->seg, gpnp->bsp, context->entityID, context->itemID, OBJ_BIOSEQ, gpnp->seqbuf, left, top, right, bottom, gpnp->scaleX, strand, FALSE, gpnp->idx, gpnp->bShowGC);

	return(TRUE);
}



static void Ing_GetLeftandRight(BioseqPtr bsp, SeqLocPtr slp, Int4Ptr left, Int4Ptr right)
{
  if (slp){
    *left=GetOffsetInBioseq(slp, bsp, SEQLOC_LEFT_END);
    *right=GetOffsetInBioseq(slp, bsp, SEQLOC_RIGHT_END);
  }
  else{
    *left=1;
    *right=bsp->length;
  }


}

extern void Ing_GetSeqLocations(IngGenomeViewerPtr igvp, BioseqPtr bsp)
{
  Int4 start, stop;

  start=igvp->from;
  stop=igvp->to;
  if (stop !=0 && stop<= bsp->length-1 && (start<stop)){
    igvp->slp=SeqLocIntNew(start, stop, (bsp->strand==Seq_strand_plus?1:2), bsp->id);
    igvp->SeqLength=stop-start+1;
  }
  else if (start==0 && stop==0){
    igvp->to=bsp->length-1;
    igvp->SeqLength = bsp->length;
  }
  igvp->bspLength=bsp->length;
  Ing_GetLeftandRight(bsp, igvp->slp, &igvp->from, &igvp->to);
   
}

extern void Ing_AddGCRect(SegmenT seg, BioseqPtr bsp, Uint2 entityID, Uint2 itemID, Uint2 itemtype, Uint1Ptr seq, Int4 left, Int4 top, Int4 right, Int4 bottom, Int4 scaleX,  Uint1 strand, Boolean needs_label, Int4 idx, Boolean bShowGC)
{
   Int4         a1, a2;
   Int4         c1, c2;
   Int4         g1, g2, t1, t2;
   Int4         gc_av, at_av, loopexit;
   Int4         i;
   Int4         j1, j2;
   Int4         len, length, start;
   Int4         n1, n2;
   FloatHi      comp2, comp_prev;
   FloatHi      comp_av;
   Char         szBuf[15]={""};
   Int4         pos, temp;
   Int4         window;
   Uint1Ptr     Clr;
   PrimitivE    prim;
   SeqLocPtr    slp=NULL;
   SeqIdPtr     sip=NULL;
   Boolean      is_loopexit =FALSE;


   a1 = a2 = c1 = c2 =  g1 = g2 = t1 = t2 = n1= n2 = 0;
   comp2 = comp_av = comp_prev = 0;
   window= scaleX*2 ; 
   start= left;
   length=right;
   if (start>length)
     {
       temp=start;
       start=length;
       length=temp;
     }
   if (needs_label){
     sip = GetSeqIdForGI(GetGIForSeqId(bsp->id));
     sip = SeqIdFindBest(sip, SEQID_GENBANK);
     if (!sip)
       sip = SeqIdFindBest(bsp->id, 0);
     SeqIdWrite(sip, szBuf, PRINTID_REPORT, 64);
     prim=AddLabel (seg, start+(length)/2, top, szBuf, SMALL_TEXT, 0,  UPPER_CENTER, 0); 
     SetPrimitiveIDs(prim, entityID, itemID, itemtype, itemtype);
   }
   if (seq==NULL || !bShowGC){ 
     prim=AddRectangle (seg, left, top , right, bottom , NO_ARROW, TRUE, idx);
     SetPrimitiveIDs(prim, entityID, itemID, itemtype, itemtype);
   }
   else{
     pos = j1 = j2 = 0;
     at_av = 0;
     gc_av = 0;
     loopexit=length;
     
     for (i=start; i<loopexit; i++)
       {
         
         if (seq[i] == 1){
           a1++;
           a2++;
         }
         else if (seq[i] == 2){
           c1++;
           c2++;
         }
         else if (seq[i] == 4){
           g1++;
           g2++;
         }
         else if (seq[i] == 8){
           t1++;
           t2++;
         }
         else{
           n1++;
           n2++;
         }
         
         j1++;
         pos++;
     
         
         if (i==loopexit-1)
           is_loopexit=TRUE;
         if ((j1 == window) || is_loopexit)
           {
             len = j1 - n1;
             gc_av = ((255)*(g1+c1))/len; 
             at_av = ((255)*(a1+t1))/len;
             
             gc_av=gc_av>165?255:gc_av; /* gc > 65% */
             Clr=Ing_PutColor(0, gc_av, 0);
             if(Clr==NULL) Clr=Ing_PutColor(128, 128, 128);
             AddAttribute(seg, COLOR_ATT, (Clr?Clr:BLACK_COLOR), 0, 0, 0, 0);
             prim=AddRectangle (seg, left+(pos-j1), top , left+pos, bottom , 0, TRUE, 0);
             SetPrimitiveIDs(prim, entityID, itemID, itemtype, idx);
             MemFree(Clr); 
           a1 = c1 = t1 = g1 = n1 = 0;
           
           gc_av = at_av = 0;
           j1 = 0;
           }
       }
     
   }
}

extern void Ing_AddRuler(SegmenT seg, Int4 height, Int4 xstart,Int4 xstop, Int4 scaleX, Int4 scaleY, Boolean add_whitespace)
{
  Int4         pos;
  Int4         xlen, x, y;
  Int4         scale_pos, i, j;
  Int4         bigtick, midtick, smalltick;
  Int4         line_Y;
  Int4         pos_2;
  Int4         pos_10;
  Int4         yBase;
  Char         scale_buf[35] = {""};	/*scale value*/
  Boolean      Decrement=FALSE;
  Int4         label;
  Int4         SeqLength;

  if (xstart>xstop)
    Decrement=TRUE;

  if (scaleX==0)
    scaleX=1;
  if (scaleY==0)
    scaleY=1;

  pos=100*scaleX;
  pos_2=pos/2;
  pos_10=pos/10;
  SeqLength=ABS(xstop-xstart)+1;
  bigtick=10*scaleY;
  midtick=7*scaleY;
  smalltick=5*scaleY;
  label=12*scaleY;
  line_Y=height+4*bigtick;/* -igvp->GrData.SegRulerIn */

  xlen=ABS(xstop-xstart)+1;
  yBase = -1*Ing_FEAT_LINEDECAL;

  AddAttribute(seg, COLOR_ATT, BLACK_COLOR, 0,0,0,0);
  AddLine(seg, 0, line_Y, SeqLength, line_Y, FALSE,-1);
  if (add_whitespace){
    AddAttribute(seg, COLOR_ATT, WHITE_COLOR, 0,0,0,0);
    AddLine(seg, 0, yBase, SeqLength, yBase, FALSE,0);
  }

  sprintf(scale_buf, "%ld", xstop);
  AddAttribute(seg, COLOR_ATT, RED_COLOR, 0,0,0,0);
  AddLabel(seg, xlen+4*scaleY, line_Y, scale_buf, SMALL_TEXT, 0, UPPER_LEFT, 0);
  
  AddAttribute(seg, COLOR_ATT, BLACK_COLOR, 0,0,0,0);
  if (!Decrement)
    {
      for (scale_pos = 0, i=0; scale_pos <= xlen; scale_pos+=pos_10, i+=pos_10)
        {
          
          if (!(scale_pos % pos))
            {
              x = i;
              y = line_Y;
              AddLine(seg, x, y, x, y-bigtick, FALSE,-1);
              sprintf(scale_buf, "%ld", (scale_pos+xstart));
              AddLabel(seg, x, y-label, scale_buf, SMALL_TEXT, 0, LOWER_CENTER, 0);
            }
          else if (!(scale_pos % (pos_2)))
            {
              x = i;
              y = line_Y;
              AddLine(seg, x, y, x, y-midtick, FALSE,-1);
              
            }
          else if (!(scale_pos % (pos_10)))
            {
              x = i;
              y = line_Y;
              AddLine(seg, x, y, x, y-smalltick, FALSE,-1);
            }
        }
    }
  else
    {
      for (scale_pos = 0, i=0, j=xstart; scale_pos <= xlen; scale_pos+=pos_10, i+=pos_10, j-=pos_10)
        {
          
          if (!(scale_pos % pos))
                {
                  x = i;
                  y = line_Y;
                  AddLine(seg, x, y, x, y-bigtick, FALSE,-1);
                  sprintf(scale_buf, "%ld", j);
                  AddLabel(seg, x, y-label, scale_buf, SMALL_TEXT, 0, LOWER_CENTER, 0);

                }
              else if (!(scale_pos % (pos_2)))
                {
                  x = i;
                  y = line_Y;
                  AddLine(seg, x, y, x, y-midtick, FALSE,-1);
                  
                }
              else if (!(scale_pos % (pos_10)))
                {
                  x = i;
                  y = line_Y;
                  AddLine(seg, x, y, x, y-smalltick, FALSE,-1);
                }
        }
    }
  
}


/*******************************************************************************

  Function : Ing_InitGrData()
  
  Purpose : initialize graphic objects

*******************************************************************************/
extern void Ing_InitGrData(IngGraphDataPtr gdp)
{
FonT f;

	/*font size*/
	f=SetSmallFont();
	if (f==NULL){
		gdp->cyChar=16;
		gdp->cxChar=8;
	}
	else{
		SelectFont(f);
		gdp->cxChar=MaxCharWidth();
		gdp->cyChar=LineHeight();
	}
	
	/*values used to draw the rulers*/
	gdp->SegBoxHeight=gdp->cyChar/2;
	gdp->SegRulerTick=gdp->cyChar/2;
	gdp->SegRulerIn=3;
	
	/*colors*/
   gdp->pClr=Ing_BuildColorTable();	
}

/*******************************************************************************

  Function : Ing_PopulateSequinGraphic()
  
  Purpose : populate sequin viewer

*******************************************************************************/
extern SegmenT Ing_PopulateSequinGraphic(SegmenT seg, BioseqPtr bsp, Uint2 entityID, Uint2 itemID, Int4 scaleX)
{
  IngExploreSegs     gpn;
  Int2               nSegments;
  PrimitivE          prim;
  IngGraphData       GrData;
  Int4               left, right;

  WatchCursor();

  Ing_InitGrData(&GrData);
  Ing_InitfeatDefFilter();
  gpn.seg=seg;
  gpn.viewer=NULL;
  gpn.idx=1;
  gpn.GrData=&GrData;
  Ing_GetLeftandRight(bsp, NULL, &left, &right);
  gpn.left=left;
  gpn.right=right;
  gpn.scaleX=scaleX;
  gpn.seqbuf=NULL;

  nSegments=SeqMgrExploreSegments(bsp, (Pointer)&gpn, Ing_ExploreSegments);

  /*this bioseq doesn't have any segments*/
  if (nSegments==0){
    AddAttribute(seg, COLOR_ATT, BLACK_COLOR, 0, 0, 0, 0);
    prim=AddRectangle(seg,0, GrData.SegBoxHeight, bsp->length, 0,NO_ARROW,TRUE, 2);
    SetPrimitiveIDs(prim, entityID, itemID, OBJ_BIOSEQ, 2);
  }

  Ing_AddRuler(seg, 2*GrData.SegBoxHeight, left, right, scaleX, 1, FALSE);
  Ing_PopFeaturesPage(bsp, seg, left, right, GrData.pClr, scaleX, 2, TRUE, FALSE);
  ArrowCursor();
  return (seg);
}

