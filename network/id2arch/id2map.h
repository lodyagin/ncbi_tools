#include <objsset.h>

#define NLM_EXTERN_LOADS {if (! SeqSetAsnLoad())return FALSE; \
                          }
#define struct_Seq_hist seqhist
#define SeqLiteralAsnRead  SeqLitAsnRead
#define SeqLiteralAsnWrite SeqLitAsnWrite
#define SeqLiteralFree     SeqLitFree
#define struct_ID2Param    struct_ID2_Param
#define ID2_SEQ_LOC_int__  ID2_SEQ_LOC_int
