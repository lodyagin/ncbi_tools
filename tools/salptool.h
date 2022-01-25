/* salpedit.h */
#ifndef __SALPTOOL__
#define __SALPTOOL__
#include <ncbi.h>
#include <objalign.h>
#include <blast.h>

typedef struct p_seqaligninfo{
	SeqAlignPtr sap;
	SeqIdPtr sip;
	Boolean used;
	struct p_seqaligninfo PNTR next;
} PSeqAlignInfo, PNTR PSeqAlignInfoPtr;

NLM_EXTERN void SeqAlignReverseOrder(SeqAlignPtr align);
NLM_EXTERN void SeqAlignSwapSeqs(SeqAlignPtr align);
NLM_EXTERN PSeqAlignInfoPtr SeqAlignToPSeqAlignInfo (SeqAlignPtr sap);
NLM_EXTERN SeqAlignPtr ReassembleSeqAlignFromPSeqAlignInfo(PSeqAlignInfoPtr alip);
NLM_EXTERN SeqAlignPtr SeqAlignSetGlobalFromLocal(SeqAlignPtr align,SeqLocPtr loc_1, SeqLocPtr loc_2, BioseqPtr bsp_1, BioseqPtr bsp_2, FILE *err_fp,Int4 MaxGap);

NLM_EXTERN Boolean TruncateAlignment (SeqAlignPtr salp, Int4 num_aln_pos, Boolean from_left);
NLM_EXTERN SeqAlignPtr MakeDiscontiguousAlignments (SeqAlignPtr salp);

#define ALNMGR_GAP           -2
#define ALNMGR_ROW_UNDEFINED -1

extern void 
AlignmentIntervalToString 
(SeqAlignPtr salp,
 Int4        row,
 Int4        start,
 Int4        stop,
 Int4        target_row,
 Boolean     view_whole_entity,
 Uint1Ptr    seqbuf,
 Uint1Ptr    alnbuf,
 Int4 PNTR   alnbuffer_len,
 Boolean     show_substitutions);

#endif
