/* dotmain.c */
#include <dotviewer.h>

/****************************************************************************

       GLOBAL VARS                                                                  
 ***************************************************************************/
#define NUMARGS  (sizeof(myargs)/sizeof(myargs[0]))

static Args myargs[] = {
    /* 0 */
    {"Query File",NULL,NULL,NULL,TRUE,'q',ARG_FILE_IN,0.0,0,NULL},
    /* 1*/
    {"Subject File",NULL, NULL,NULL,TRUE,'s',ARG_FILE_IN,0.0,0,NULL},
    /* 2 */
    {"plus-minus =0 plus-plus =1","1", NULL,NULL,TRUE,'r',ARG_INT,0.0,0,NULL},
    /* 3 */
    {"query start","0", NULL,NULL,TRUE,'a',ARG_INT,0.0,0,NULL},
    /* 4 */
    {"query stop, default if 5","5", NULL,NULL,TRUE,'b',ARG_INT,0.0,0,NULL},
    /* 5 */
    {"subject start","0", NULL,NULL,TRUE,'c',ARG_INT,0.0,0,NULL},
    /* 6 */
    {"subject stop, default if 5", "5", NULL,NULL,TRUE,'d',ARG_INT,0.0,0,NULL},
    /* 7 */
    { "Word size: nucleotide[4 - 11], protein[1 or 2], default if 0",  "0", NULL, NULL, TRUE, 'w', ARG_INT, 0.0, 0, NULL},
    /*8 */
    { "Number of  takehits to keep",  "100000", NULL, NULL, TRUE, 'k', ARG_INT, 0.0, 0, NULL},
    /* 9*/
    {"Alignment asn file",NULL, NULL,NULL,TRUE,'l',ARG_FILE_IN,0.0,0,NULL}
};

static FILE        *qfile=NULL, *sfile=NULL, *afile=NULL;
static Int4        q_start,q_stop,s_start, s_stop, tree_limit, word_size;




static void DOT_ProcessDialog (Int4 word_size, Int4 tree_limit, Int4 q_start, Int4 q_stop, Int4 s_start, Int4 s_stop, Boolean plus_plus)
{
  
  SeqEntryPtr   ssep, qsep;
  SeqAnnotPtr   sanp;
  SeqAlignPtr   sap;
  BioseqPtr     qbsp, sbsp;
  SeqLocPtr     slp1, slp2;
  DOTMainDataPtr    mip;
  Uint2         entityID1, entityID2, entityID3, datatype;
  Pointer       dataptr;
  

  if (! AllObjLoad ())
    {
      ErrPostEx(SEV_FATAL, 0, 0, "dotplot: SeqEntryLoad failed\n");
      return ;
    }

  if (qfile){
   while ((dataptr = ReadAsnFastaOrFlatFile (qfile, &datatype, &entityID1, FALSE, FALSE, TRUE, FALSE)) != NULL)
     {
       qsep= GetTopSeqEntryForEntityID(entityID1);
       entityID1 = SeqMgrIndexFeatures(0, qsep);
       qbsp = qsep->data.ptrvalue;

     }
  }
  if (sfile){
   while ((dataptr = ReadAsnFastaOrFlatFile (sfile, &datatype, &entityID2, FALSE, FALSE, TRUE, FALSE)) != NULL)
     {
       ssep= GetTopSeqEntryForEntityID(entityID2);
       entityID2 = SeqMgrIndexFeatures(0, ssep);
       sbsp = ssep->data.ptrvalue;

     }
  }

  if (qfile || sfile){
  if (qbsp == NULL || sbsp == NULL)
    {
      ErrPostEx(SEV_FATAL, 0, 0, "Unable to obtain bioseq\n");
      goto end;
    } 

   if (!((ISA_aa (qbsp->mol) && ISA_aa (sbsp->mol))||(ISA_na(qbsp->mol) && ISA_na (sbsp->mol))))
    {
      ErrPostEx(SEV_ERROR, 0, 0, "DOT -missmatched sequence types");
      goto end;
    }

  if (q_stop==5) 
    q_stop=qbsp->length;   
  if (s_stop==5) 
    s_stop=sbsp->length; 
  
  if (qbsp->length<q_start || sbsp->length<s_start || q_start>q_stop|| s_start>s_stop || s_start<0 || q_start<0)
    {
      ErrPostEx(SEV_FATAL, 0, 0, "Bad sequence offset values");
      goto end;
    }

  if (plus_plus)
    slp1= SeqLocIntNew(q_start, q_stop-1, 1, qbsp->id);
  else
    slp1= SeqLocIntNew(q_start, q_stop-1, 2, qbsp->id);
  slp2 = SeqLocIntNew(s_start, s_stop-1,1, sbsp->id);

  mip=DOT_CreateAndStorebyLoc(slp1, slp2, 8, 100000);

  if (afile){
    dataptr = ReadAsnFastaOrFlatFile (afile, &datatype, NULL, FALSE, FALSE, TRUE, FALSE);
    if (!dataptr) return;
    sanp = (SeqAnnotPtr)(dataptr);
    sap = (SeqAlignPtr)(sanp->data);
    AlnMgrIndexSeqAlign(sap);
    fclose(afile);
  }

  DOT_MakeMainViewer(mip, sap);
  end:
  fclose (qfile);
  fclose (sfile);
  }
  else if (afile){
    dataptr = ReadAsnFastaOrFlatFile (afile, &datatype, &entityID3, FALSE, FALSE, TRUE, FALSE);
    if (!dataptr){
      ErrPostEx(SEV_FATAL, 0, 0, "no seqalign found");
      return ;
    } 
    sanp = (SeqAnnotPtr)(dataptr);
    sap = (SeqAlignPtr)(sanp->data);
    if (!sap){
      ErrPostEx(SEV_FATAL, 0, 0, "no seqalign found");
      return ;
    }
    AlnMgrIndexSeqAlign(sap);
    DOT_AlignPlotGivenSeqAlign(sap);
    fclose(afile);
  }


  return ;
}



Int2 Main ()
{

  Int2  orient;
  CharPtr txt=NULL;

    if (! GetArgs ("dotplot", NUMARGS, myargs))
        {
            return (1);
        }
    else
        {    
          
          if((txt = myargs[0].strvalue)!=NULL){
            qfile = FileOpen(txt, "r");
          }

          if((txt = myargs[1].strvalue)!=NULL){
            sfile = FileOpen(txt, "r");
          }
          else{
            if (qfile)
              sfile=qfile;
          }
          orient = (Int2) myargs[2].intvalue; /* plus-minus=0,  plus-plus=1 */
          q_start = myargs[3].intvalue;
          q_stop = myargs[4].intvalue;
          s_start = myargs[5].intvalue;
          s_stop = myargs[6].intvalue;
          word_size= myargs[7].intvalue;
          tree_limit = myargs[8].intvalue;
          if((txt = myargs[9].strvalue)!=NULL){
            afile = FileOpen(txt, "r");
          }          

          DOT_ProcessDialog(word_size, tree_limit, q_start, q_stop, s_start, s_stop, orient);
        }
    return 1;

}
