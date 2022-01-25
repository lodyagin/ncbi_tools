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
#define ERR_SEQ_INST_MolNuclAcid  1,25
#define ERR_SEQ_INST_ConflictingBiomolTech  1,26
#define ERR_SEQ_INST_SeqIdNameHasSpace  1,27
#define ERR_SEQ_INST_IdOnMultipleBioseqs  1,28
#define ERR_SEQ_INST_DuplicateSegmentReferences  1,29
#define ERR_SEQ_INST_TrailingX  1,30
#define ERR_SEQ_INST_BadSeqIdFormat  1,31
#define ERR_SEQ_INST_PartsOutOfOrder  1,32
#define ERR_SEQ_INST_BadSecondaryAccn  1,33
#define ERR_SEQ_INST_ZeroGiNumber  1,34
#define ERR_SEQ_INST_RnaDnaConflict  1,35
#define ERR_SEQ_INST_HistoryGiCollision  1,36
#define ERR_SEQ_INST_GiWithoutAccession  1,37
#define ERR_SEQ_INST_MultipleAccessions  1,38
#define ERR_SEQ_INST_HistAssemblyMissing  1,39
#define ERR_SEQ_INST_TerminalNs  1,40
#define ERR_SEQ_INST_UnexpectedIdentifierChange  1,41
#define ERR_SEQ_INST_InternalNsInSeqLit  1,42
#define ERR_SEQ_INST_SeqLitGapLength0  1,43
#define ERR_SEQ_INST_TpaAssmeblyProblem  1,44
#define ERR_SEQ_INST_SeqLocLength  1,45
#define ERR_SEQ_INST_MissingGaps  1,46
#define ERR_SEQ_INST_CompleteTitleProblem  1,47
#define ERR_SEQ_INST_CompleteCircleProblem  1,48
#define ERR_SEQ_INST_BadHTGSeq  1,49
#define ERR_SEQ_INST_GapInProtein  1,50
#define ERR_SEQ_INST_BadProteinStart  1,51
#define ERR_SEQ_INST_TerminalGap  1,52
#define ERR_SEQ_DESCR  2,0
#define ERR_SEQ_DESCR_BioSourceMissing  2,1
#define ERR_SEQ_DESCR_InvalidForType  2,2
#define ERR_SEQ_DESCR_FileOpenCollision  2,3
#define ERR_SEQ_DESCR_Unknown  2,4
#define ERR_SEQ_DESCR_NoPubFound  2,5
#define ERR_SEQ_DESCR_NoOrgFound  2,6
#define ERR_SEQ_DESCR_MultipleBioSources  2,7
#define ERR_SEQ_DESCR_NoMolInfoFound  2,8
#define ERR_SEQ_DESCR_BadCountryCode  2,9
#define ERR_SEQ_DESCR_NoTaxonID  2,10
#define ERR_SEQ_DESCR_InconsistentBioSources  2,11
#define ERR_SEQ_DESCR_MissingLineage  2,12
#define ERR_SEQ_DESCR_SerialInComment  2,13
#define ERR_SEQ_DESCR_BioSourceNeedsFocus  2,14
#define ERR_SEQ_DESCR_BadOrganelle  2,15
#define ERR_SEQ_DESCR_MultipleChromosomes  2,16
#define ERR_SEQ_DESCR_BadSubSource  2,17
#define ERR_SEQ_DESCR_BadOrgMod  2,18
#define ERR_SEQ_DESCR_InconsistentProteinTitle  2,19
#define ERR_SEQ_DESCR_Inconsistent  2,20
#define ERR_SEQ_DESCR_ObsoleteSourceLocation  2,21
#define ERR_SEQ_DESCR_ObsoleteSourceQual  2,22
#define ERR_SEQ_DESCR_StructuredSourceNote  2,23
#define ERR_SEQ_DESCR_UnnecessaryBioSourceFocus  2,24
#define ERR_SEQ_DESCR_RefGeneTrackingWithoutStatus  2,25
#define ERR_SEQ_DESCR_UnwantedCompleteFlag  2,26
#define ERR_SEQ_DESCR_CollidingPublications  2,27
#define ERR_SEQ_DESCR_TransgenicProblem  2,28
#define ERR_SEQ_DESCR_TaxonomyLookupProblem  2,29
#define ERR_SEQ_DESCR_MultipleTitles  2,30
#define ERR_SEQ_DESCR_RefGeneTrackingOnNonRefSeq  2,31
#define ERR_SEQ_DESCR_BioSourceInconsistency  2,32
#define ERR_SEQ_DESCR_FastaBracketTitle  2,33
#define ERR_SEQ_DESCR_MissingText  2,34
#define ERR_GENERIC  3,0
#define ERR_GENERIC_NonAsciiAsn  3,1
#define ERR_GENERIC_Spell  3,2
#define ERR_GENERIC_AuthorListHasEtAl  3,3
#define ERR_GENERIC_MissingPubInfo  3,4
#define ERR_GENERIC_UnnecessaryPubEquiv  3,5
#define ERR_GENERIC_BadPageNumbering  3,6
#define ERR_GENERIC_MedlineEntryPub  3,7
#define ERR_GENERIC_BadDate  3,8
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
#define ERR_SEQ_PKG_FeaturePackagingProblem  4,10
#define ERR_SEQ_PKG_GenomicProductPackagingProblem  4,11
#define ERR_SEQ_PKG_InconsistentMolInfoBiomols  4,12
#define ERR_SEQ_PKG_ArchaicFeatureLocation  4,13
#define ERR_SEQ_PKG_ArchaicFeatureProduct  4,14
#define ERR_SEQ_PKG_GraphPackagingProblem  4,15
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
#define ERR_SEQ_FEAT_PseudoCdsHasProduct  5,27
#define ERR_SEQ_FEAT_IllegalDbXref  5,28
#define ERR_SEQ_FEAT_FarLocation  5,29
#define ERR_SEQ_FEAT_DuplicateFeat  5,30
#define ERR_SEQ_FEAT_UnnecessaryGeneXref  5,31
#define ERR_SEQ_FEAT_TranslExceptPhase  5,32
#define ERR_SEQ_FEAT_TrnaCodonWrong  5,33
#define ERR_SEQ_FEAT_BothStrands  5,34
#define ERR_SEQ_FEAT_CDSgeneRange  5,35
#define ERR_SEQ_FEAT_CDSmRNArange  5,36
#define ERR_SEQ_FEAT_OverlappingPeptideFeat  5,37
#define ERR_SEQ_FEAT_SerialInComment  5,38
#define ERR_SEQ_FEAT_MultipleCDSproducts  5,39
#define ERR_SEQ_FEAT_FocusOnBioSourceFeature  5,40
#define ERR_SEQ_FEAT_PeptideFeatOutOfFrame  5,41
#define ERR_SEQ_FEAT_InvalidQualifierValue  5,42
#define ERR_SEQ_FEAT_MultipleMRNAproducts  5,43
#define ERR_SEQ_FEAT_mRNAgeneRange  5,44
#define ERR_SEQ_FEAT_TranscriptLen  5,45
#define ERR_SEQ_FEAT_TranscriptMismatches  5,46
#define ERR_SEQ_FEAT_CDSproductPackagingProblem  5,47
#define ERR_SEQ_FEAT_DuplicateInterval  5,48
#define ERR_SEQ_FEAT_PolyAsiteNotPoint  5,49
#define ERR_SEQ_FEAT_ImpFeatBadLoc  5,50
#define ERR_SEQ_FEAT_LocOnSegmentedBioseq  5,51
#define ERR_SEQ_FEAT_UnnecessaryCitPubEquiv  5,52
#define ERR_SEQ_FEAT_ImpCDShasTranslation  5,53
#define ERR_SEQ_FEAT_ImpCDSnotPseudo  5,54
#define ERR_SEQ_FEAT_MissingMRNAproduct  5,55
#define ERR_SEQ_FEAT_AbuttingIntervals  5,56
#define ERR_SEQ_FEAT_CollidingGeneNames  5,57
#define ERR_SEQ_FEAT_MultiIntervalGene  5,58
#define ERR_SEQ_FEAT_FeatContentDup  5,59
#define ERR_SEQ_FEAT_BadProductSeqId  5,60
#define ERR_SEQ_FEAT_RnaProductMismatch  5,61
#define ERR_SEQ_FEAT_MissingCDSproduct  5,62
#define ERR_SEQ_FEAT_BadTrnaCodon  5,63
#define ERR_SEQ_FEAT_BadTrnaAA  5,64
#define ERR_SEQ_FEAT_OnlyGeneXrefs  5,65
#define ERR_SEQ_FEAT_UTRdoesNotAbutCDS  5,66
#define ERR_SEQ_FEAT_BadConflictFlag  5,67
#define ERR_SEQ_FEAT_ConflictFlagSet  5,68
#define ERR_SEQ_FEAT_LocusTagProblem  5,69
#define ERR_SEQ_FEAT_CollidingLocusTags  5,70
#define ERR_SEQ_FEAT_AltStartCodon  5,71
#define ERR_SEQ_FEAT_PartialsInconsistent  5,72
#define ERR_SEQ_FEAT_GenesInconsistent  5,73
#define ERR_SEQ_FEAT_DuplicateTranslExcept  5,74
#define ERR_SEQ_FEAT_TranslExceptAndRnaEditing  5,75
#define ERR_SEQ_FEAT_NoNameForProtein  5,76
#define ERR_SEQ_FEAT_TaxonDbxrefOnFeature  5,77
#define ERR_SEQ_FEAT_UnindexedFeature  5,78
#define ERR_SEQ_FEAT_CDSmRNAmismatch  5,79
#define ERR_SEQ_FEAT_UnnecessaryException  5,80
#define ERR_SEQ_FEAT_LocusTagProductMismatch  5,81
#define ERR_SEQ_FEAT_MrnaTransFail  5,82
#define ERR_SEQ_FEAT_PseudoCdsViaGeneHasProduct  5,83
#define ERR_SEQ_FEAT_MissingGeneXref  5,84
#define ERR_SEQ_FEAT_FeatureCitationProblem  5,85
#define ERR_SEQ_FEAT_NestedSeqLocMix  5,86
#define ERR_SEQ_FEAT_WrongQualOnFeature  5,87
#define ERR_SEQ_FEAT_MissingQualOnFeature  5,88
#define ERR_SEQ_FEAT_CodonQualifierUsed  5,89
#define ERR_SEQ_FEAT_UnknownFeatureQual  5,90
#define ERR_SEQ_FEAT_BadCharInAuthorName  5,91
#define ERR_SEQ_FEAT_PolyATail  5,92
#define ERR_SEQ_FEAT_ProteinNameEndsInBracket  5,93
#define ERR_SEQ_FEAT_CDSwithMultipleMRNAs  5,94
#define ERR_SEQ_FEAT_MultipleEquivBioSources  5,95
#define ERR_SEQ_FEAT_MultipleEquivPublications  5,96
#define ERR_SEQ_FEAT_BadFullLengthFeature  5,97
#define ERR_SEQ_FEAT_RedundantFields  5,98
#define ERR_SEQ_FEAT_CDSwithNoMRNAOverlap  5,99
#define ERR_SEQ_FEAT_FeatureProductInconsistency  5,100
#define ERR_SEQ_FEAT_ImproperBondLocation  5,101
#define ERR_SEQ_FEAT_GeneXrefWithoutGene  5,102
#define ERR_SEQ_FEAT_SeqFeatXrefProblem  5,103
#define ERR_SEQ_FEAT_ProductFetchFailure  5,104
#define ERR_SEQ_ALIGN  6,0
#define ERR_SEQ_ALIGN_SeqIdProblem  6,1
#define ERR_SEQ_ALIGN_StrandRev  6,2
#define ERR_SEQ_ALIGN_DensegLenStart  6,3
#define ERR_SEQ_ALIGN_StartLessthanZero  6,4
#define ERR_SEQ_ALIGN_StartMorethanBiolen  6,5
#define ERR_SEQ_ALIGN_EndLessthanZero  6,6
#define ERR_SEQ_ALIGN_EndMorethanBiolen  6,7
#define ERR_SEQ_ALIGN_LenLessthanZero  6,8
#define ERR_SEQ_ALIGN_LenMorethanBiolen  6,9
#define ERR_SEQ_ALIGN_SumLenStart  6,10
#define ERR_SEQ_ALIGN_AlignDimSeqIdNotMatch  6,11
#define ERR_SEQ_ALIGN_SegsDimSeqIdNotMatch  6,12
#define ERR_SEQ_ALIGN_FastaLike  6,13
#define ERR_SEQ_ALIGN_NullSegs  6,14
#define ERR_SEQ_ALIGN_SegmentGap  6,15
#define ERR_SEQ_ALIGN_SegsDimOne  6,16
#define ERR_SEQ_ALIGN_AlignDimOne  6,17
#define ERR_SEQ_ALIGN_Segtype  6,18
#define ERR_SEQ_ALIGN_BlastAligns  6,19
#define ERR_SEQ_GRAPH  7,0
#define ERR_SEQ_GRAPH_GraphMin  7,1
#define ERR_SEQ_GRAPH_GraphMax  7,2
#define ERR_SEQ_GRAPH_GraphBelow  7,3
#define ERR_SEQ_GRAPH_GraphAbove  7,4
#define ERR_SEQ_GRAPH_GraphByteLen  7,5
#define ERR_SEQ_GRAPH_GraphOutOfOrder  7,6
#define ERR_SEQ_GRAPH_GraphBioseqLen  7,7
#define ERR_SEQ_GRAPH_GraphSeqLitLen  7,8
#define ERR_SEQ_GRAPH_GraphSeqLocLen  7,9
#define ERR_SEQ_GRAPH_GraphStartPhase  7,10
#define ERR_SEQ_GRAPH_GraphStopPhase  7,11
#define ERR_SEQ_GRAPH_GraphDiffNumber  7,12
#define ERR_SEQ_GRAPH_GraphACGTScore  7,13
#define ERR_SEQ_GRAPH_GraphNScore  7,14
#define ERR_SEQ_GRAPH_GraphGapScore  7,15
#define ERR_SEQ_GRAPH_GraphOverlap  7,16

#endif
