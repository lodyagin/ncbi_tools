#ifndef _TOPORG_
#define _TOPORG_

#include <stdio.h>
#include <ncbi.h>
#include <seqport.h>

#define Seq_descr_GIBB_mod_dna        0
#define Seq_descr_GIBB_mod_rna        1
#define Seq_descr_GIBB_mod_partial    10
#define Seq_descr_GIBB_mod_complete    11
#define Seq_descr_GIBB_mod_est        20
#define Seq_descr_GIBB_mod_sts        21
#define Seq_descr_GIBB_mod_gss        22

void toporg PROTO((SeqEntryPtr sep));
void ChkSegset PROTO((SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent));
void ChkNucProt PROTO((SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent));
ValNodePtr GetDescr PROTO((ValNodePtr PNTR descr));
Boolean check_GIBB PROTO((ValNodePtr descr));
ValNodePtr SrchSegChoice PROTO((SeqEntryPtr sep, Uint1 choice));
void SrchSegSeqMol PROTO((SeqEntryPtr sep));
Boolean CheckSegDescrChoice PROTO((SeqEntryPtr sep, Uint1 choice));
void CleanUpSeqDescrChoice PROTO((SeqEntryPtr sep, Uint1 choice));
void StripProtXref PROTO((SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent));
void CheckMaps PROTO((SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent));
void StripMaps PROTO((SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent));
void MapsToGenref PROTO((SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent));
void FindNewLineage PROTO((SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent));
void FindOldLineage PROTO((SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent));
void NewLineage PROTO((SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent));
void OldLineage PROTO((SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent));
void OldPubs PROTO((SeqEntryPtr sep));
void NewPubs PROTO((SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent));
void MoveSetPubs PROTO((SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent));
void DeletePubs PROTO((SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent));
void ChangeCitSub PROTO((SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent));
void CheckCitSubNew PROTO((ValNodePtr vnp));
void CmpPub PROTO((SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent));
void MoveSegmPubs PROTO((SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent));
void MoveNPPubs PROTO((SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent));
void MovePopPhyMutPubs (SeqEntryPtr sep);
ValNodePtr AddToList PROTO((ValNodePtr list, ValNodePtr check, PubdescPtr pdp));
BioSourcePtr BioSourceMerge PROTO((BioSourcePtr host, BioSourcePtr guest));
BioSourcePtr BioSourceCommon PROTO((BioSourcePtr host, BioSourcePtr guest));
void StripBSfromTop PROTO((SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent));
void StripBSfromParts PROTO((SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent));
Boolean CmpOrgById PROTO((BioSourcePtr b1, BioSourcePtr b2));
extern void NormalizeSegSeqMolInfo PROTO((SeqEntryPtr sep));

#endif
