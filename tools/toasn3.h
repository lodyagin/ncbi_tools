#ifndef _TOASN3_
#define _TOASN3_

#include <stdio.h>
#include <ncbi.h>
#include <seqport.h>

#define INFO_ASNOLD 0
#define INFO_ASNNEW 1
#define ERR_REJECT  2
#define ERR_INPUT   4

typedef struct {
	ValNodePtr	pept;
	ValNodePtr	cds;
} SeqFeatArr, PNTR SeqFeatArrPtr;

typedef struct {
	CharPtr name;
	Uint1   num;
} ORGMOD;

typedef struct {
	BioseqSetPtr bioset;
	ValNodePtr   list;
} PubList, PNTR PubListPtr;

typedef struct {
	Boolean first;
	ValNodePtr   list;
} PubSetList, PNTR PubSetListPtr;

typedef struct qualmap{
	CharPtr name;
	Boolean same;
} QualMap, PNTR QualMapPtr;

typedef struct bsmap{
	BioSourcePtr bsp;
	Boolean same;
} BSMap, PNTR BSMapPtr;

typedef struct orgfix {
	SeqEntryPtr 	contains;
	Boolean 		desc;               /* descriptor containing org-ref */
	SeqFeatPtr 		sfp;                /* or feature containing the org-ref */
	SeqFeatPtr 		imp;                /* ImpFeat source */
	OrgRefPtr 		orp;
	ValNodePtr		modif;
	BioSourcePtr 	bsp;
	Int4 			index;
	struct orgfix PNTR next;
} OrgFix, PNTR OrgFixPtr;

typedef struct molfix {
	SeqEntryPtr		contains;
	Uint1 			mol;               		/*  mol_type */
	ValNodePtr 		modif;               /* descriptor containing modif */
	Uint1 			method;               	/*  method */
	MolInfoPtr 		molinfo;
	Int4 			index;
	struct molfix PNTR next;
} MolFix, PNTR MolFixPtr;

typedef struct toasn3 {
	OrgFixPtr ofp;
	MolFixPtr mfp;
	Boolean had_biosource;
	Boolean had_molinfo;
} ToAsn3, PNTR ToAsn3Ptr;

Int4 ToAsn4 PROTO((SeqEntryPtr sep));
Int4 SeqEntryPubsAsn4 PROTO((SeqEntryPtr sep));
Int4 SeqEntryPubsAsn3 PROTO((SeqEntryPtr sep));
Int4 SeqEntryToAsn3 PROTO((SeqEntryPtr sep, Boolean strip, Boolean correct, Boolean taxserver, SeqEntryFunc taxfun));
Int4 SeqEntryToAsn3Ex PROTO((SeqEntryPtr sep, Boolean strip, Boolean correct, Boolean taxserver, SeqEntryFunc taxfun, SeqEntryFunc taxmerge));
Int2 seq_loc_compare PROTO(( SeqLocPtr a, SeqLocPtr b));
void compare_quals PROTO((GBQualPtr PNTR qual1, GBQualPtr PNTR qual2));
Boolean feat_join PROTO((SeqFeatPtr f1, SeqFeatPtr f2, SeqFeatPtr head));
void count_join PROTO((SeqFeatPtr f1, SeqFeatPtr f2));
void CountSourceFeat PROTO((SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent));
void CorrectSourceFeat PROTO((SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent));

Boolean CheckLocWhole(BioseqPtr bsp, SeqLocPtr slp);
void FindOrg PROTO((SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent));
void CheckQuals PROTO((BioSourcePtr bsp, GBQualPtr sfp));
void CheckQualsWithComm PROTO((BioSourcePtr bsp, SeqFeatPtr sfp));
MolInfoPtr new_info PROTO((MolInfoPtr mfi));
MolInfoPtr ModToMolInfo PROTO((MolInfoPtr mfi, Uint1 mod));
void ModToBiosource PROTO((BioSourcePtr bsp, Uint1 mod));
void StripOld PROTO((SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent));
void CkOrg PROTO((SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent));
void MergeBSinDescr PROTO((SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent));
Int4 FixNucProtSet PROTO((SeqEntryPtr sep));
void CheckBS PROTO((SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent));
Int4 BSComparison PROTO((BioSourcePtr one, BioSourcePtr two));
Int4 BSComparisonEx PROTO((BioSourcePtr one, BioSourcePtr two, Boolean clone));
Int2 BioSourceToGeneticCode PROTO((BioSourcePtr biop));
ValNodePtr GetMultBiosource PROTO((SeqEntryPtr sep));

/****************************************************
* Does a SeqEntryExplore, calling ImpFeatToCdregion() on
*  each feature table found
*
****************************************************/
void EntryChangeImpFeat PROTO((SeqEntryPtr sep));

/****************************************************
*  Changes a SeqFeat of type Imp-feat CDS to a real
*  CdRegion in place (sfp does not change)
*  returns TRUE if the change happened
*  returns FALSE if no changes were made
*    (so also returns FALSE if not an Imp-feat of type CDS)
*
*****************************************************/
Boolean ImpFeatToCdregion PROTO((SeqFeatPtr sfp));

void CdCheck PROTO((SeqEntryPtr sep, FILE *fp));

/****************************************************
*  Creates source string from BioSource structure
*  Compare with GBBlock.source
*  deletes GBBlock.source if it's the same as from BioSource
*
*****************************************************/
void EntryChangeGBSource PROTO((SeqEntryPtr sep));

/****************************************************
*  Finds multiple biosource descriptors on the same chain with
*  the same taxonomic name, moves subsource, orgmod, and some
*  other blocks, conservatively, deletes second biosource
*
*****************************************************/
void EntryMergeDupBioSources PROTO((SeqEntryPtr sep));

/****************************************************
*  Checks GBBlock.source, .taxonomy, and .div, removes any empty
*  GBBlock descriptors, and returns TRUE if information (other
*  than PAT or SYN division) exists in the these fields
*
*****************************************************/
Boolean EntryCheckGBBlock PROTO((SeqEntryPtr sep));

void EntryChangeImpFeatToProt PROTO((SeqEntryPtr sep));

void CombineBSFeat PROTO((SeqEntryPtr sep));
void ChangeReplaceToQual PROTO((SeqFeatPtr sfp));
void AddReplaceQual PROTO((SeqFeatPtr sfp, CharPtr p));
Boolean SeqEntryMoveDbxrefs  PROTO((SeqEntryPtr sep));

/* functions moved from Sequin */

void NormalizePeriodsOnInitials (SeqEntryPtr sep);
void MoveRnaGBQualProductToName (SeqEntryPtr sep);
void MoveProtGBQualProductToName (SeqEntryPtr sep);
void MoveCdsGBQualProductToName (SeqEntryPtr sep);
void MoveFeatGBQualsToFields (SeqEntryPtr sep);
void StripTitleFromProtsInNucProts (SeqEntryPtr sep);

void GetRidOfEmptyFeatsDescCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent);

/* from move_cds (S. Bazhin) */

Uint2 move_cds PROTO((SeqEntryPtr sep));

/* more functions moved from Sequin, placed in toporg.c */

extern void CleanUpPseudoProducts (Uint2 entityID, SeqEntryPtr sep);
extern void CleanupGenbankCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent);
extern void MergeAdjacentAnnotsCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent);
extern void CleanupEmptyFeatCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent);
extern void RemoveBioSourceOnPopSet (SeqEntryPtr sep, OrgRefPtr master);
extern Boolean NoBiosourceOrTaxonId (SeqEntryPtr sep);
extern void ExtendGeneFeatIfOnMRNA (Uint2 entityID, SeqEntryPtr sep);
extern void ConvertFullLenSourceFeatToDesc (SeqEntryPtr sep);

/* SeriousSeqEntryCleanup combines many of the above cleanups */

extern void SeriousSeqEntryCleanup (SeqEntryPtr sep, SeqEntryFunc taxfun, SeqEntryFunc taxmerge);

#endif

