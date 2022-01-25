MODULE valid
$$ SEQ_INST, 1
$^   ExtNotAllowed, 1
# This is a comment
A Bioseq "extension" is used for special classes of Bioseq. This class of
Bioseq should not have one but it does. This is probably a software error.

$^   ExtBadOrMissing, 2
This class of Bioseq requires an "extension" but it is missing or of the
wrong type. This is probably a software error.

$^   SeqDataNotFound, 3
No actual sequence data was found on this Bioseq. This is probably a
software problem.

$^   SeqDataNotAllowed, 4
The wrong type of sequence data was found on this Bioseq. This is probably a
software problem.

$^   ReprInvalid, 5
This Bioseq has an invalid representation class. This is probably a software
error.

$^   CircularProtein, 6
This protein Bioseq is represented as circular. Circular topology is
normally used only for certain DNA molecules, for example, plasmids.

$^   DSProtein, 7
This protein Bioseq has strandedness indicated. Strandedness is normally a
property only of DNA sequences. Please unset the strandedness.

$^   MolNotSet, 8
It is not clear whether this sequence is nucleic acid or protein. Please set
the appropriate molecule type (Bioseq.mol).

$^   MolOther, 9
Most sequences are either nucleic acid or protein. However, the molecule
type (Bioseq.mol) is set to "other". It should probably be set to nucleic
acid or a protein.

$^   FuzzyLen, 10
This sequence is marked as having an uncertain length, but the length is
known exactly.

$^   InvalidLen, 11
The length indicated for this sequence is invalid. This is probably a
software error.

$^   InvalidAlphabet, 12
This Bioseq has an invalid alphabet (e.g. protein codes on a nucleic acid or
vice versa). This is probably a software error.

$^   SeqDataLenWrong, 13
The length of this Bioseq does not agree with the length of the actual data.
This is probably a software error.

$^   SeqPortFail, 14
Something is very wrong with this entry. The validator cannot open a SeqPort
on the Bioseq. Further testing cannot be done.

$^   InvalidResidue, 15
Invalid residue codes were found in this Bioseq.

$^   StopInProtein, 16
Stop codon symbols were found in this protein Bioseq.

$^   PartialInconsistent, 17
This segmented sequence is described as complete or incomplete in several
places, but these settings are inconsistent.

$^   ShortSeq, 18
This Bioseq is unusually short (less than 4 amino acids or less than 11
nucleic acids). GenBank does not usually accept such short sequences.

$^   NoIdOnBioseq, 19
No SeqIds were found on this Bioseq. This is probably a software error.

$^   BadDeltaSeq, 20
Delta sequences should only be HTGS-1 or HTGS-2.

$^   LongHtgsSequence, 21
HTGS-1 or HTGS-2 sequences must be < 350 KB in length.

$^   LongLiteralSequence, 22
Delta literals must be < 350 KB in length.

$^   SequenceTooLong, 23
Individual sequences must be < 350 KB in length.

$^   ConflictingIdsOnBioseq, 24
Two SeqIds of the same class was found on this Bioseq. This is probably a software error.

$$ SEQ_DESCR, 2

$^   BioSourceMissing, 1
The biological source of this sequence has not been described correctly.  A
Bioseq must have a BioSource descriptor that covers the entire molecule.
Additional BioSource features may also be added to recombinant molecules,
natural or otherwise, to designate the parts of the molecule. Please add the
source information.

$^   InvalidForType, 2
This descriptor cannot be used with this Bioseq. A descriptor placed at the
BioseqSet level applies to all of the Bioseqs in the set. Please make sure
the descriptor is consistent with every sequence to which it applies.

$^   Inconsistent, 3
There are two descriptors of the same type which are inconsistent with each
other. Please make them consistent.

$^   Unknown, 4
An unknown or "other" modifier was used.

$^   NoPubFound, 5
No publications were found in this entry which refer to this Bioseq. If a
publication descriptor is added to a BioseqSet, it will apply to all of the
Bioseqs in the set. A publication feature should be used if the publication
applies only to a subregion of a sequence.

$^   NoOrgFound, 6
This entry does not specify the organism that was the source of the sequence.
Please name the organism.

$^   MultipleBioSources, 7
There are multiple BioSource or OrgRef descriptors in the same chain with
the same taxonomic name. Their information should be combined into a single
BioSource descriptor.

$^   NoMolInfoFound, 8
This sequence does not have a Mol-info descriptor applying to it.  This indicates
genomic vs. message, sequencing technique, and whether the sequence is incomplete.

$$ GENERIC, 3

$^   NonAsciiAsn, 1
There is a non-ASCII type character in this entry.

$^   Spell, 2
There is a potentially misspelled word in this entry.

$^   AuthorListHasEtAl, 3
The author list contains et al, which should be replaced with the
remaining author names.

$^   MissingPubInfo, 4
The publication is missing essential information, such as title or
authors.

$$ SEQ_PKG, 4

$^   NoCdRegionPtr, 1
A protein is found in this entry, but the coding region has not been
described. Please add a CdRegion feature to the nucleotide Bioseq.

$^   NucProtProblem, 2
Both DNA and protein sequences were expected, but one of the two seems to be
missing. Perhaps this is the wrong package to use.

$^   SegSetProblem, 3
A segmented sequence was expected, but it was not found. Perhaps this is the
wrong package to use.

$^   EmptySet, 4
No Bioseqs were found in this BioseqSet. Is that what was intended?

$^   NucProtNotSegSet, 5
A nuc-prot set should not contain any other BioseqSet except segset.

$^   SegSetNotParts, 6
A segset should not contain any other BioseqSet except parts.

$^   SegSetMixedBioseqs, 7
A segset should not contain both nucleotide and protein Bioseqs.

$^   PartsSetMixedBioseqs, 8
A parts set should not contain both nucleotide and protein Bioseqs.

$^   PartsSetHasSets, 9
A parts set should not contain BioseqSets.

$$ SEQ_FEAT, 5

$^   InvalidForType, 1
This feature type is illegal on this type of Bioseq.

$^   PartialProblem, 2
There are several places in an entry where a sequence can be described as
either partial or complete. In this entry, these settings are inconsistent.
Make sure that the location and product Seq-locs, the Bioseqs, and the
SeqFeat partial flag all agree in describing this SeqFeat as partial or
complete.

$^   InvalidType, 3
A feature with an invalid type has been detected. This is most likely a
software problem.

$^   Range, 4
The coordinates describing the location of a feature do not fall within the
sequence itself. A feature location or a product Seq-loc is out of range of
the Bioseq it points to.

$^   MixedStrand, 5
Mixed strands (plus and minus) have been found in the same location. While
this is biologically possible, it is very unusual. Please check that this is
really what you mean.

$^   SeqLocOrder, 6
This location has intervals that are out of order. While whis is biologically
possible, it is very unusual. Please check that this is really what you mean.

$^   CdTransFail, 7
A fundamental error occurred in software while attempting to translate this
coding region. It is either a software problem or sever data corruption.

$^   StartCodon, 8
An illegal start codon was used. Some possible explanations are: (1) the
wrong genetic code may have been selected; (2) the wrong reading frame may
be in use; or (3) the coding region may be incomplete at the 5' end, in
which case a partial location should be indicated.

$^   InternalStop, 9
Internal stop codons are found in the protein sequence. Some possible
explanations are: (1) the wrong genetic code may have been selected; (2) the
wrong reading frame may be in use; (3) the coding region may be incomplete
at the 5' end, in which case a partial location should be indicated; or (4)
the CdRegion feature location is incorrect.

$^   NoProtein, 10
Normally a protein sequence is supplied. This sequence can then be compared
with the translation of the coding region. In this entry, no protein Bioseq
was found, and the comparison could not be made.

$^   MisMatchAA, 11
The protein sequence that was supplied is not identical to the translation
of the coding region. Mismatching amino acids are found between these two
sequences.

$^   TransLen, 12
The protein sequence that was supplied is not the same length as the
translation of the coding region. Please determine why they are different.

$^   NoStop, 13
A coding region that is complete should have a stop codon at the 3'end.  A
stop codon was not found on this sequence, although one was expected.

$^   TranslExcept, 14
An unparsed \transl_except qualifier was found. This indicates a parser
problem.

$^   NoProtRefFound, 15
The name and description of the protein is missing from this entry. Every
protein Bioseq must have one full-length Prot-ref feature to provide this
information.

$^   NotSpliceConsensus, 16
Splice junctions typically have GT as the first two bases of the intron
(splice donor) and AG as the last two bases of the intron (splice acceptor).
This intron does not conform to that pattern.

$^ OrfCdsHasProduct, 17
A coding region flagged as orf has a protein product.  There should be no
protein product bioseq on an orf.

$^ GeneRefHasNoData, 18
A gene feature exists with no locus name or other fields filled in.

$^ ExceptInconsistent, 19
A coding region has an exception gbqual but the excpt flag is not set.

$^ ProtRefHasNoData, 20
A protein feature exists with no name or other fields filled in.

$^ GenCodeMismatch, 21
The genetic code stored in the BioSource is different than that for this CDS.

$^ RNAtype0, 22
RNA type 0 (unknown RNA) should be type 255 (other).

$^ UnknownImpFeatKey, 23
An import feature has an unrecognized key.

$^ UnknownImpFeatQual, 24
An import feature has an unrecognized qualifier.

$^ WrongQualOnImpFeat, 25
This qualifier is not legal for this feature.

$^ MissingQualOnImpFeat, 26
An essential qualifier for this feature is missing.

$^ PsuedoCdsHasProduct, 27
A coding region flagged as pseudo has a protein product.  There should be no
protein product bioseq on a pseudo CDS.

$^ IllegalDbXref, 28
The database in a cross-reference is not on the list of officially
recognized database abbreviations.

$^ FarLocation, 29
The location has a reference to a bioseq that is not packaged in this record.

