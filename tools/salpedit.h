/* salpedit.h */
#ifndef __SALPEDIT__
#define __SALPEDIT__
#include <blast.h>

typedef struct p_seqaligninfo{
	SeqAlignPtr sap;
	SeqIdPtr sip;
	Boolean used;
	struct p_seqaligninfo PNTR next;
} PSeqAlignInfo, PNTR PSeqAlignInfoPtr;
NLM_EXTERN SeqAlignPtr SeqAlignSplitGappedBlast(SeqLocPtr slp1, CharPtr progname, CharPtr database, ValNodePtr *other_returns, ValNodePtr *error_returns, Int4 split_len, Int4 overlap_len, BLAST_OptionsBlkPtr options);

NLM_EXTERN Boolean MergeTwoAlignList(SeqAlignPtr h_align_1, SeqAlignPtr PNTR p_align_2, Int4 from, Int4 to, Int2 order);

NLM_EXTERN PSeqAlignInfoPtr SeqAlignToPSeqAlignInfo (SeqAlignPtr sap);
NLM_EXTERN SeqAlignPtr ReassembleSeqAlignFromPSeqAlignInfo(PSeqAlignInfoPtr alip);


#endif
