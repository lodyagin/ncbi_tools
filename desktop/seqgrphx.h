#ifndef _SEQGRPHX_
#define _SEQGRPHX_

#ifdef __cplusplus
extern "C" {
#endif

/* structures */

typedef struct seqgraphview
{
  FORM_MESSAGE_BLOCK
  BioseqPtr          bsp;
  GrouP              ftype, mtype, ptype;
  Uint2              entityID, procID, userKEY;
  Uint4              itemID;
} SeqGraphView, PNTR SeqGraphViewPtr;

/* prototypes */

Int2 LIBCALLBACK PCCPredictFunc (Pointer data); 

/* sequinx defines */


#define REGISTER_GROUP_PCC ObjMgrProcLoadEx (OMPROC_FILTER, \
        "Coiled-coil", "Coiled-coil", \
        OBJ_SEQFEAT, 0, OBJ_SEQFEAT, 0, \
        NULL, PCCPredictFunc, PROC_PRIORITY_DEFAULT, "Analysis")

#ifdef __cplusplus
}
#endif

#endif
