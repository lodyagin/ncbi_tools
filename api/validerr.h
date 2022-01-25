#ifndef __MODULE_valid__
#define __MODULE_valid__

#define ERR_SEQ_INST  1,0
#define ERR_SEQ_INST_ExtNotAllowed  1,1
#define ERR_SEQ_INST_ExtBadOrMissing  1,2
#define ERR_SEQ_INST_SeqDataNotFound  1,3
#define ERR_SEQ_INST_SeqDataNotAllowed  1,4
#define ERR_SEQ_INST_ReprInvalid  1,5
#define ERR_SEQ_INST_CircularProtein  1,6
#define ERR_SEQ_INST_DSProtein  1,7
#define ERR_SEQ_INST_MolNotSet  1,8
#define ERR_SEQ_INST_MolOther  1,9
#define ERR_SEQ_INST_FuzzyLen  1,10
#define ERR_SEQ_INST_InvalidLen  1,11
#define ERR_SEQ_INST_InvalidAlphabet  1,12
#define ERR_SEQ_INST_SeqDataLenWrong  1,13
#define ERR_SEQ_INST_SeqPortFail  1,14
#define ERR_SEQ_INST_InvalidResidue  1,15
#define ERR_SEQ_INST_StopInProtein  1,16
#define ERR_SEQ_INST_PartialInconsistent  1,17
#define ERR_SEQ_INST_ShortSeq  1,18
#define ERR_SEQ_INST_NoIdOnBioseq  1,19
#define ERR_SEQ_INST_BadDeltaSeq  1,20
#define ERR_SEQ_INST_LongHtgsSequence  1,21
#define ERR_SEQ_INST_LongLiteralSequence  1,22
#define ERR_SEQ_INST_SequenceExceeds350kbp  1,23
#define ERR_SEQ_INST_ConflictingIdsOnBioseq  1,24
#define ERR_SEQ_DESCR  2,0
#define ERR_SEQ_DESCR_BioSourceMissing  2,1
#define ERR_SEQ_DESCR_InvalidForType  2,2
#define ERR_SEQ_DESCR_Inconsistent  2,3
#define ERR_SEQ_DESCR_Unknown  2,4
#define ERR_SEQ_DESCR_NoPubFound  2,5
#define ERR_SEQ_DESCR_NoOrgFound  2,6
#define ERR_SEQ_DESCR_MultipleBioSources  2,7
#define ERR_SEQ_DESCR_NoMolInfoFound  2,8
#define ERR_SEQ_DESCR_BadCountryCode  2,9
#define ERR_SEQ_DESCR_NoTaxonID  2,10
#define ERR_SEQ_DESCR_InconsistentBioSources  2,11
#define ERR_GENERIC  3,0
#define ERR_GENERIC_NonAsciiAsn  3,1
#define ERR_GENERIC_Spell  3,2
#define ERR_GENERIC_AuthorListHasEtAl  3,3
#define ERR_GENERIC_MissingPubInfo  3,4
#define ERR_SEQ_PKG  4,0
#define ERR_SEQ_PKG_NoCdRegionPtr  4,1
#define ERR_SEQ_PKG_NucProtProblem  4,2
#define ERR_SEQ_PKG_SegSetProblem  4,3
#define ERR_SEQ_PKG_EmptySet  4,4
#define ERR_SEQ_PKG_NucProtNotSegSet  4,5
#define ERR_SEQ_PKG_SegSetNotParts  4,6
#define ERR_SEQ_PKG_SegSetMixedBioseqs  4,7
#define ERR_SEQ_PKG_PartsSetMixedBioseqs  4,8
#define ERR_SEQ_PKG_PartsSetHasSets  4,9
#define ERR_SEQ_FEAT  5,0
#define ERR_SEQ_FEAT_InvalidForType  5,1
#define ERR_SEQ_FEAT_PartialProblem  5,2
#define ERR_SEQ_FEAT_InvalidType  5,3
#define ERR_SEQ_FEAT_Range  5,4
#define ERR_SEQ_FEAT_MixedStrand  5,5
#define ERR_SEQ_FEAT_SeqLocOrder  5,6
#define ERR_SEQ_FEAT_CdTransFail  5,7
#define ERR_SEQ_FEAT_StartCodon  5,8
#define ERR_SEQ_FEAT_InternalStop  5,9
#define ERR_SEQ_FEAT_NoProtein  5,10
#define ERR_SEQ_FEAT_MisMatchAA  5,11
#define ERR_SEQ_FEAT_TransLen  5,12
#define ERR_SEQ_FEAT_NoStop  5,13
#define ERR_SEQ_FEAT_TranslExcept  5,14
#define ERR_SEQ_FEAT_NoProtRefFound  5,15
#define ERR_SEQ_FEAT_NotSpliceConsensus  5,16
#define ERR_SEQ_FEAT_OrfCdsHasProduct  5,17
#define ERR_SEQ_FEAT_GeneRefHasNoData  5,18
#define ERR_SEQ_FEAT_ExceptInconsistent  5,19
#define ERR_SEQ_FEAT_ProtRefHasNoData  5,20
#define ERR_SEQ_FEAT_GenCodeMismatch  5,21
#define ERR_SEQ_FEAT_RNAtype0  5,22
#define ERR_SEQ_FEAT_UnknownImpFeatKey  5,23
#define ERR_SEQ_FEAT_UnknownImpFeatQual  5,24
#define ERR_SEQ_FEAT_WrongQualOnImpFeat  5,25
#define ERR_SEQ_FEAT_MissingQualOnImpFeat  5,26
#define ERR_SEQ_FEAT_PsuedoCdsHasProduct  5,27
#define ERR_SEQ_FEAT_IllegalDbXref  5,28
#define ERR_SEQ_FEAT_FarLocation  5,29
#define ERR_SEQ_FEAT_DuplicateFeat  5,30
#define ERR_SEQ_FEAT_UnnecessaryGeneXref  5,31
#define ERR_SEQ_FEAT_TranslExceptPhase  5,32
#define ERR_SEQ_FEAT_TrnaCodonWrong  5,33

#endif
