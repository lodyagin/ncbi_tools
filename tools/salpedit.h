/* salpedit.h */
#ifndef __SALPEDIT__
#define __SALPEDIT__
#include <blast.h>
NLM_EXTERN SeqAlignPtr SeqAlignSplitBlastTwoSeq(SeqLocPtr slp1, SeqLocPtr slp2, 
		Int4 split_len, Int4 overlap_len, BLAST_OptionsBlkPtr options);

NLM_EXTERN Boolean MergeTwoAlignList (SeqAlignPtr h_align_1, SeqAlignPtr PNTR p_align_2, Int4 from, Int4 to, Int2 order);
#endif
