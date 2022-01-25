/*  valid.c
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE                          
*               National Center for Biotechnology Information
*                                                                          
*  This software/database is a "United States Government Work" under the   
*  terms of the United States Copyright Act.  It was written as part of    
*  the author's official duties as a United States Government employee and 
*  thus cannot be copyrighted.  This software/database is freely available 
*  to the public for use. The National Library of Medicine and the U.S.    
*  Government have not placed any restriction on its use or reproduction.  
*                                                                          
*  Although all reasonable efforts have been taken to ensure the accuracy  
*  and reliability of the software and data, the NLM and the U.S.          
*  Government do not and cannot warrant the performance or results that    
*  may be obtained by using this software or data. The NLM and the U.S.    
*  Government disclaim all warranties, express or implied, including       
*  warranties of performance, merchantability or fitness for any particular
*  purpose.                                                                
*                                                                          
*  Please cite the author in any work or product based on this material.   
*
* ===========================================================================
*
* File Name:  valid.c
*
* Author:  James Ostell
*   
* Version Creation Date: 1/1/94
*
* $Revision: 6.76 $
*
* File Description:  Sequence editing utilities
*
* Modifications:  
* --------------------------------------------------------------------------
* Date	   Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
* $Log: valid.c,v $
* Revision 6.76  1999/09/06 21:36:03  kans
* ValidateSeqEntry sets scope
*
* Revision 6.75  1999/08/24 17:44:01  kans
* removed Wagad from country list
*
* Revision 6.74  1999/08/24 15:22:17  kans
* added Galapagos Islands and Wagad to the country list
*
* Revision 6.73  1999/08/18 20:24:49  kans
* self-recursive call of CheckForInconsistentBiosources was not using tmp, but original sep, resulting in stack overflow in complex records
*
* Revision 6.72  1999/08/17 19:46:12  kans
* ValidatePopSet posts ERR_SEQ_DESCR_InconsistentBioSources
*
* Revision 6.71  1999/08/03 00:13:02  kans
* vsp->suppressContext now causes simplified locations to be written, seqidworst fastashort no locus
*
* Revision 6.70  1999/07/29 15:41:48  kans
* changed Serbia and Montenegro to Yugoslavia
*
* Revision 6.69  1999/07/22 22:04:35  kans
* added suppressContext flag
*
* Revision 6.68  1999/07/15 22:37:32  kans
* ValidateBioSource called once per biosource, not once per bioseq
*
* Revision 6.67  1999/07/15 20:39:22  kans
* suppress no pub warning if seq-submit, which has a cit-sub
*
* Revision 6.66  1999/06/24 19:33:24  kans
* corrected country list
*
* Revision 6.65  1999/06/22 17:15:49  kans
* added ERR_SEQ_DESCR_NoTaxonID
*
* Revision 6.64  1999/06/18 20:57:46  kans
* using collab approved country list
*
* Revision 6.63  1999/06/18 20:21:04  kans
* implemented ERR_SEQ_DESCR_BadCountryCode, indexed descr callback sets proper itemtype, itemID for click responsiveness
*
* Revision 6.62  1999/06/15 20:04:03  kans
* no org or pub anywhere on record now reports context of first bioseq for batch processing
*
* Revision 6.61  1999/06/15 19:45:42  kans
* changed SequenceTooLong to SequenceExceeds350kbp
*
* Revision 6.60  1999/06/14 16:14:20  kans
* added ERR_SEQ_FEAT_TrnaCodonWrong check
*
* Revision 6.59  1999/06/11 18:31:16  kans
* added ERR_SEQ_FEAT_TranslExceptPhase
*
* Revision 6.58  1999/06/09 21:34:29  kans
* stop in protein message gives gene and protein name for reading report later
*
* Revision 6.57  1999/05/07 15:31:20  kans
* added ERR_SEQ_FEAT_UnnecessaryGeneXref
*
* Revision 6.56  1999/05/05 19:11:41  kans
* for no pubs or biosource anywhere, needed to set vsp->gcp for ValidErr/ErrPostItem
*
* Revision 6.55  1999/05/05 13:03:14  kans
* no org or pub anywhere after clearing error counts
*
* Revision 6.54  1999/05/03 20:06:35  kans
* if no pubs or no biosource, report only once, not once per bioseq
*
* Revision 6.53  1999/03/31 20:57:48  kans
* htgs phase 1 and 2 messages also check for phase 0
*
* Revision 6.52  1999/03/04 19:55:49  kans
* inconsistent create_date messages now sev_warning
*
* Revision 6.51  1999/02/25 21:53:58  kans
* relax duplicate feature severity to warning if label or comment are different, or if FEATDEF_PUB
*
* Revision 6.50  1999/02/16 22:19:02  kans
* fixed interval comparison in duplicate feature detection
*
* Revision 6.49  1999/02/02 16:39:10  kans
* added ERR_SEQ_FEAT_DuplicateFeat
*
* Revision 6.48  1999/01/05 23:20:50  kans
* SpliceCheckEx does not check exon junction if partial
*
* Revision 6.47  1998/12/14 22:27:28  kans
* CdTransCheck now deals with termination by polyA
*
* Revision 6.46  1998/12/07 20:00:56  kans
* meant to set bcp = NULL, not bsp = NULL, crashed with segmented protein
*
* Revision 6.45  1998/10/26 20:57:45  kans
* check gene and prot db fields for IllegalDbXref
*
* Revision 6.44  1998/10/23 15:25:57  kans
* added FarLocation warning
*
* Revision 6.43  1998/10/22 16:05:57  kans
* removed labeltype parameter from SeqMgrIndexFeatures, changed index parameter/field to Uint2
*
* Revision 6.42  1998/10/21 14:32:11  kans
* on invalid feature for bioseq, restore itemid itemid and itemtype to avoid weird(er) click association - need to rewrite valid with new index functions, which will give proper items
*
* Revision 6.41  1998/10/20 20:18:10  kans
* mRNA feature is invalid on an mRNA (cDNA) bioseq
*
* Revision 6.40  1998/10/20 18:12:54  kans
* invalid for type (e.g., intron on mRNA) now coerces gcp to have feature itemtype, itemID for selection
*
* Revision 6.39  1998/10/15 17:29:18  kans
* import feature of mat_, sig_, and transit_peptide now flagged as invalid for type
*
* Revision 6.38  1998/09/22 13:12:01  kans
* locationFilter parameter to explore features function
*
* Revision 6.37  1998/09/21 17:29:35  kans
* precursor rna can have intron feature
*
* Revision 6.36  1998/09/17 16:38:14  kans
* added ERR_SEQ_DESCR_NoMolInfoFound
*
* Revision 6.35  1998/09/01 19:25:27  kans
* context parameter in get best protein, get cds/rna given product
*
* Revision 6.34  1998/08/28 22:25:56  kans
* keep track of last biomol, tech, completeness in multiple molinfo descriptors
*
* Revision 6.33  1998/08/26 21:07:48  kans
* added check for ERR_SEQ_INST_ConflictingIdsOnBioseq
*
* Revision 6.32  1998/08/10 16:05:15  kans
* copy some old descriptor checks to Molinfo
*
* Revision 6.31  1998/07/23 14:25:38  kans
* intron and CAAT_signal are illegal on mRNA - first checks molinfo, then resorts to Seq_mol_rna as mRNA criterion
*
* Revision 6.30  1998/07/16 16:06:56  kans
* use ObjMgrGetEntityIDForChoice instead of ObjMgrGetEntityIDForPointer for SeqEntryPtr
*
* Revision 6.29  1998/07/14 18:10:33  kans
* invalid feature for nucleotide now says nucleotide, not protein
*
* Revision 6.28  1998/07/06 18:01:52  kans
* added LIBCALLBACK to SeqMgrExplore function callbacks
*
* Revision 6.27  1998/07/02 17:53:43  kans
* useSeqMgrIndexes field added to ValidStructPtr, validator can use either old (nested gathers) or new (SeqMgr indexing) method
*
* Revision 6.26  1998/06/24 18:49:15  kans
* added missing BioseqContextFree
*
* Revision 6.25  1998/06/22 20:13:21  kans
* gencode mismatch reports biosource and cds codes
*
* Revision 6.24  1998/06/12 20:05:53  kans
* fixed unix compiler warnings
*
* Revision 6.23  1998/04/16 15:12:15  kans
* slight fix to frame > 1 and not at splice site test
*
* Revision 6.22  1998/04/15 21:59:25  kans
* added ERR_SEQ_FEAT_IllegalDbXref
*
* Revision 6.21  1998/04/14 20:57:36  kans
* check for mixed bioseqs in segset, parts set, and for sets within parts set
*
* Revision 6.20  1998/04/14 19:11:25  kans
* improvements to PartialAtSpliceSite and frame > 1 check
*
* Revision 6.19  1998/04/14 18:55:56  kans
* cds frame > 1 but not 5prime partial now also checks for PartialAtSpliceSite
*
* Revision 6.18  1998/04/13 18:10:38  kans
* warn if CDS frame > 1 but not 5prime partial
*
* Revision 6.17  1998/04/02 15:45:51  kans
* MolInfoPtr had not been obtained for Seq_repr_raw for HTGS test on long sequences
*
* Revision 6.16  1998/03/30 17:35:22  kans
* check raw bioseq for htgs flags if greater than 350kb
*
* Revision 6.15  1998/03/18 20:41:50  kans
* SpliceCheck only on mRNA (not all RNAs) and CDS
*
* Revision 6.14  1998/03/09 17:48:46  kans
* OBJ_SEQSUB_CIT now satisfies need for publication
*
* Revision 6.13  1998/02/19 17:21:15  shavirin
* Added check for NULL in ValidErr() function
*
* Revision 6.12  1998/02/18 20:34:55  kans
* added ERR_GENERIC_MissingPubInfo
*
* Revision 6.11  1998/02/09 20:35:35  kans
* calls ERR_SEQ_FEAT_PseudoCdsHasProduct
*
* Revision 6.10  1998/01/30 21:05:54  kans
* check for ERR_SEQ_DESCR_MultipleBioSources
*
* Revision 6.9  1998/01/30 20:29:48  kans
* added PartialAtSpliceSite check
*
* Revision 6.8  1998/01/13 15:34:50  kans
* gbqual_citation satisfied by sfp->cit
*
* Revision 6.7  1998/01/10 00:05:36  kans
* added ValidateImpFeat
*
* Revision 6.6  1998/01/06 03:07:57  ostell
* in comparison of cdregion genetic code to biosource genetic code, set defaults
* to 0 instead of -1 to fix default behavior on building submission.
*
* Revision 6.5  1997/12/18 21:51:43  kans
* warn on cds/biosource genetic code conflict, rna type 0
*
* Revision 6.4  1997/11/14 17:10:13  kans
* added checks for bioseq length > 350K (based on Cavanaugh request)
*
* Revision 6.3  1997/08/27 20:11:02  kans
* order gene should in fact have partial flag set
*
* Revision 6.2  1997/08/27 19:48:32  kans
* print feature product seqloc
*
* Revision 6.1  1997/08/27 14:15:51  kans
* gene of order should not cause partial error
*
* Revision 6.0  1997/08/25 18:08:25  madden
* Revision changed to 6.0
*
* Revision 5.24  1997/08/13 18:52:51  kans
* new packaging errors set to SEV_FATAL
*
* Revision 5.23  1997/08/13 15:36:53  kans
* added NucProtNotSegSet and SegSetNotParts (Bazhin)
*
* Revision 5.22  1997/07/07 21:28:11  kans
* existing bad start codon check was being bypassed, so new one was added
*
* Revision 5.21  1997/07/07 15:00:28  kans
* signal or transit peptide do not need names
*
* Revision 5.20  1997/07/02 19:44:09  kans
* added check for et al, changed symbol names for empty gene and prot feature
*
* Revision 5.19  1997/06/24 16:39:12  kans
* fixed Digital Unix compiler complaint
*
* Revision 5.18  1997/06/19 18:39:51  vakatov
* [WIN32,MSVC++]  Adopted for the "NCBIOBJ.LIB" DLL'ization
*
* Revision 5.17  1997/05/29 17:25:16  kans
* splice check and trans check not done if excpt
*
* Revision 5.16  1997/05/28 19:10:32  kans
* added check for empty protref
*
* Revision 5.15  1997/05/20 21:11:38  kans
* warnings for delta seq not htgs1 or 2, cds orf with product, gene with no fields, cds exception gbqual without excpt
*
* Revision 5.14  1997/04/24 20:39:20  kans
* invalid splice sites are warning level unless app property forces to error
*
 * Revision 5.13  1997/03/17  21:43:28  kans
 * added closing bracket to bioseq length indication
 *
 * Revision 5.12  1997/02/20  13:50:33  ostell
 * added length check on segmented sequence back
 *
 * Revision 5.11  1996/11/22  17:23:20  kans
 * splice errors on exon imp-feats are now severity warning, since there is
 * no way of knowing which are the unspliced ends of the first and last exon
 *
 * Revision 5.10  1996/11/04  16:29:55  kans
 * app property allows splice check for exon features, and rare GC splice
 * donor has separate warning message
 *
 * Revision 5.9  1996/10/16  20:31:16  ostell
 * added length check for delta sequences
 * added CdTrnsCheck for exception and pseudo
 *
 * Revision 5.8  1996/08/21  14:08:26  ostell
 * rmoved kludge for big sequences
 *
 * Revision 5.7  1996/08/19  02:45:49  ostell
 * added check in BioseqContect for more than 30n bioseqs to control
 * feature checkes
 *
 * Revision 5.6  1996/08/06  19:56:03  kans
 * for SEQLOC_WHOLE, must call SeqIdFindBest on bsp->id
 *
 * Revision 5.5  1996/08/01  18:58:00  kans
 * on pseudo cds, suppress CdTransCheck, SpliceCheck
 *
 * Revision 5.4  1996/06/19  00:35:32  ostell
 * added check for ragged end of CdRegion
 *
 * Revision 5.1  1996/06/16  04:16:05  ostell
 * added support for delta seq
 *
 * Revision 5.0  1996/05/28  13:23:23  ostell
 * Set to revision 5.0
 *
 * Revision 4.19  1996/05/03  18:59:13  kans
 * up to 5 stops still allows mismatch report, which includes nuc position
 *
 * Revision 4.18  1996/04/01  16:31:47  ostell
 * fix to preserver eror message count between invocations
 *
 * Revision 4.17  1996/03/15  20:01:14  ostell
 * in SpliceCheck, give accession of sequence with bad junction
 *
 * Revision 4.16  1996/03/08  14:48:02  kans
 * fixed typos in ValidateSeqEntry scope memset, use as parameter
 *
 * Revision 4.15  1996/03/06  20:43:59  ostell
 * added scoping to validation
 *
 * Revision 4.14  1996/03/05  19:54:29  kans
 * added biosource to two switch statements
 *
 * Revision 4.13  1996/03/03  16:59:34  ostell
 * added SpellCheckPub() to look at more Pub types
 *
 * Revision 4.12  1996/03/02  03:41:43  ostell
 * fix to correctly identigy splice junctions on minus strand
 *
 * Revision 4.11  1996/02/26  22:06:37  ostell
 * finished gatherized version of spell check on descriptors
 *
 * Revision 4.10  1996/02/19  19:58:05  ostell
 * added support for Code-break and tRNA.anticodon
 *
 * Revision 4.9  1996/01/23  23:10:10  kans
 * implemented onlyspell and justwarnonspell code
 *
 * Revision 4.8  1995/12/07  01:55:37  ostell
 * fix to check for NULL on bioseqset parent
 *
 * Revision 4.7  1995/12/07  01:38:56  ostell
 * added Splice error flag
 *
 * Revision 4.6  1995/12/06  22:11:23  ostell
 * changed wording of SpliceCheck message
 *
 * Revision 4.5  1995/12/06  06:08:57  ostell
 * lowered warning levels on partial messages
 * added SpliceCheck()
 *
 * Revision 4.4  1995/08/16  18:21:52  epstein
 * correct declaration of static functions to be consistent with function prototypes
 *
 * Revision 4.3  1995/08/04  18:41:02  madden
 * removed "|SpellErr|" SpellCallBack.
 *
 * Revision 4.2  1995/08/03  12:45:56  madden
 * Set ValNodePtr in SpellCheckBioseqDescr; added "SpellErr" to ErrPosting.
 *
 * Revision 4.1  1995/08/02  22:21:50  madden
 * gatherized the spell functions.
 *
 * Revision 4.0  1995/07/26  13:49:01  ostell
 * force revision to 4.0
 *
 * Revision 1.14  1995/06/03  13:45:47  ostell
 * changes made in valid to use gather functions and ErrPostItem instead
 * of previous custom functions
 *
 * Revision 1.13  1995/05/15  21:46:05  ostell
 * added Log line
 *
*
*
* ==========================================================================
*/

static char *this_module = "valid";
#define THIS_MODULE this_module

static char *this_file = __FILE__;
#define THIS_FILE this_file

#include <ncbi.h>
#include <objfdef.h>
#include <valid.h>
#include <validerr.h>
#include <sqnutils.h>
#include <gbftdef.h>
#include <gbfeat.h>
#include <objsub.h>
#include <asn2ffp.h>
#include <explore.h>
#include <subutil.h>

/*****************************************************************************
*
*   NOTE: look at all the ValidErr calls with severity=0. Some should be
*   bumped up later. Look also for string "PARSER"
*
*****************************************************************************/



#ifdef VAR_ARGS
#include <varargs.h>
#else
#include <stdarg.h>
#endif

static ValidStructPtr globalvsp;  /* for spell checker */

NLM_EXTERN void CDECL  ValidErr VPROTO((ValidStructPtr vsp, int severity, int code1, int code2, const char *fmt, ...));
static void ValidateBioseqInst (GatherContextPtr gcp);
static void ValidateBioseqContext (GatherContextPtr gcp);
static void ValidateBioseqSet (GatherContextPtr gcp);
static void SpellCheckSeqDescr(GatherContextPtr gcp);
NLM_EXTERN void CdTransCheck(ValidStructPtr vsp, SeqFeatPtr sfp);
void ValidateFeatureTable (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent);
NLM_EXTERN void ValidateSeqFeat(GatherContextPtr gcp);
NLM_EXTERN void ValidateSeqLoc(ValidStructPtr vsp, SeqLocPtr slp, CharPtr prefix);
NLM_EXTERN Boolean PatchBadSequence(BioseqPtr bsp);
NLM_EXTERN CharPtr FindIDForEntry (SeqEntryPtr sep, CharPtr buf);
NLM_EXTERN void SpellCheckSeqFeat(GatherContextPtr gcp);
NLM_EXTERN void SpellCheckString (ValidStructPtr vsp, CharPtr str);
NLM_EXTERN void SpliceCheck(ValidStructPtr vsp, SeqFeatPtr sfp);
static void SpliceCheckEx(ValidStructPtr vsp, SeqFeatPtr sfp, Boolean checkAll);
static void ValidateBioSource (ValidStructPtr vsp, GatherContextPtr gcp, BioSourcePtr biop);

/*****************************************************************************
*
*   Perform Validation Checks on a SeqEntry
*
*****************************************************************************/

NLM_EXTERN void ValidStructClear (ValidStructPtr vsp)  /* 0 out a ValidStruct */
{
	CharPtr errbuf;
	Int2 cutoff;
	Boolean patch_seq;
	SpellCheckFunc spellfunc;
	SpellCallBackFunc spellcallback;
	Boolean onlyspell;
	Boolean justwarnonspell;
	Boolean useSeqMgrIndexes;
	Boolean suppressContext;

	if (vsp == NULL) return;

	errbuf = vsp->errbuf;
	cutoff = vsp->cutoff;
	patch_seq = vsp->patch_seq;
	spellfunc = vsp->spellfunc;
	spellcallback = vsp->spellcallback;
	onlyspell = vsp->onlyspell;
	justwarnonspell = vsp->justwarnonspell;
	useSeqMgrIndexes = vsp->useSeqMgrIndexes;
	suppressContext = vsp->suppressContext;
	MemSet((VoidPtr)vsp, 0, sizeof(ValidStruct));
	vsp->errbuf = errbuf;
	vsp->cutoff = cutoff;
	vsp->patch_seq = patch_seq;
	vsp->spellfunc = spellfunc;
	vsp->spellcallback = spellcallback;
	vsp->onlyspell = onlyspell;
	vsp->justwarnonspell = justwarnonspell;
	vsp->useSeqMgrIndexes = useSeqMgrIndexes;
	vsp->suppressContext = suppressContext;
	return;
}

NLM_EXTERN ValidStructPtr ValidStructNew (void)
{
	ValidStructPtr vsp;

	vsp = (ValidStructPtr)MemNew(sizeof(ValidStruct));
	return vsp;
}

NLM_EXTERN ValidStructPtr ValidStructFree (ValidStructPtr vsp)
{
	if (vsp == NULL) return vsp;

	MemFree(vsp->errbuf);
	return (ValidStructPtr)MemFree(vsp);
}

/*****************************************************************************
*
*   ValidErr()
*
*****************************************************************************/

static void ChangeSeqIdToBestID (SeqIdPtr sip)

{
  BioseqPtr  bsp;
  SeqIdPtr   id;
  Pointer    pnt;

  if (sip == NULL) return;
  bsp = BioseqFindCore (sip);
  if (bsp == NULL) return;
  id = SeqIdDup (SeqIdFindWorst (bsp->id));
  if (id == NULL) return;
  /* now remove SeqId contents to reuse SeqId valnode */
  pnt = sip->data.ptrvalue;
  switch (sip->choice) {
        case SEQID_LOCAL:      /* local */
            ObjectIdFree((ObjectIdPtr)pnt);
            break;
        case SEQID_GIBBSQ:      /* gibbseq */
        case SEQID_GIBBMT:      /* gibbmt */
            break;
        case SEQID_GIIM:      /* giimid */
            GiimFree((GiimPtr)pnt);
            break;
        case SEQID_GENBANK:      /* genbank */
        case SEQID_EMBL:      /* embl */
        case SEQID_PIR:      /* pir   */
        case SEQID_SWISSPROT:      /* swissprot */
        case SEQID_OTHER:     /* other */
        case SEQID_DDBJ:
		case SEQID_PRF:
            TextSeqIdFree((TextSeqIdPtr)pnt);
            break;
        case SEQID_PATENT:      /* patent seq id */
            PatentSeqIdFree((PatentSeqIdPtr)pnt);
            break;
        case SEQID_GENERAL:     /* general */
            DbtagFree((DbtagPtr)pnt);
            break;
        case SEQID_GI:     /* gi */
            break;
		case SEQID_PDB:
			PDBSeqIdFree((PDBSeqIdPtr)pnt);
			break;
  }
  sip->choice = id->choice;
  sip->data.ptrvalue = id->data.ptrvalue;
  SeqIdStripLocus (sip);
}

static void ChangeSeqLocToBestID (SeqLocPtr slp)

{
  SeqLocPtr      loc;
  PackSeqPntPtr  psp;
  SeqBondPtr     sbp;
  SeqIntPtr      sinp;
  SeqIdPtr       sip;
  SeqPntPtr      spp;

  while (slp != NULL) {
    switch (slp->choice) {
      case SEQLOC_NULL :
        break;
      case SEQLOC_EMPTY :
      case SEQLOC_WHOLE :
        sip = (SeqIdPtr) slp->data.ptrvalue;
        ChangeSeqIdToBestID (sip);
        break;
      case SEQLOC_INT :
        sinp = (SeqIntPtr) slp->data.ptrvalue;
        if (sinp != NULL) {
          sip = sinp->id;
          ChangeSeqIdToBestID (sip);
        }
        break;
      case SEQLOC_PNT :
        spp = (SeqPntPtr) slp->data.ptrvalue;
        if (spp != NULL) {
          sip = spp->id;
          ChangeSeqIdToBestID (sip);
        }
        break;
      case SEQLOC_PACKED_PNT :
        psp = (PackSeqPntPtr) slp->data.ptrvalue;
        if (psp != NULL) {
          sip = psp->id;
          ChangeSeqIdToBestID (sip);
        }
        break;
      case SEQLOC_PACKED_INT :
      case SEQLOC_MIX :
      case SEQLOC_EQUIV :
        loc = (SeqLocPtr) slp->data.ptrvalue;
        while (loc != NULL) {
          ChangeSeqLocToBestID (loc);
          loc = loc->next;
        }
        break;
      case SEQLOC_BOND :
        sbp = (SeqBondPtr) slp->data.ptrvalue;
        if (sbp != NULL) {
          spp = (SeqPntPtr) sbp->a;
          if (spp != NULL) {
            sip = spp->id;
            ChangeSeqIdToBestID (sip);
          }
          spp = (SeqPntPtr) sbp->b;
          if (spp != NULL) {
            sip = spp->id;
            ChangeSeqIdToBestID (sip);
          }
        }
        break;
      case SEQLOC_FEAT :
        break;
      default :
        break;
    }
    slp = slp->next;
  }
}

static Int2 WorstBioseqLabel (BioseqPtr bsp, CharPtr buffer, Int2 buflen, Uint1 content)
{
	CharPtr tmp;
	Char label[40];
	Int2 diff, len;
	SeqIdPtr sip;
	AsnModulePtr amp;
	AsnTypePtr ratp, matp;

	if ((bsp == NULL) || (buflen < 1))
		return 0;

	len = buflen;
	label[0] = '\0';

	if (content != OM_LABEL_TYPE)
	{
		sip = SeqIdStripLocus (SeqIdDup (SeqIdFindWorst (bsp->id)));
		SeqIdWrite (sip, label, PRINTID_FASTA_SHORT, 39);
		SeqIdFree (sip);
		if (content == OM_LABEL_CONTENT)
			return LabelCopy(buffer, label, buflen);

		diff = LabelCopyExtra(buffer, label, buflen, NULL, ": ");
		buflen -= diff;
		buffer += diff;
	}

	amp = AsnAllModPtr ();
	ratp = AsnTypeFind (amp, "Seq-inst.repr");
	matp = AsnTypeFind (amp, "Seq-inst.mol");

	label[0] = '\0';
	tmp = label;
	tmp = StringMove(tmp, AsnEnumTypeStr(ratp, (Int2)(bsp->repr)));
	tmp = StringMove(tmp, ", ");
	tmp = StringMove(tmp, AsnEnumTypeStr(matp, (Int2)(bsp->mol)));
	sprintf(tmp, " len= %ld", (long)(bsp->length));
	diff = LabelCopy(buffer, label, buflen);
	buflen -= diff;
	buffer += diff;

	if (content != OM_LABEL_SUMMARY)
		return (len - buflen);
	
	return (len - buflen);     /* SUMMARY not done yet */
}

#ifdef VAR_ARGS
NLM_EXTERN void CDECL  ValidErr (vsp, severity, code1, code2, fmt, va_alist)
  ValidStructPtr vsp;
  int severity;
  int code1;
  int code2;
  const char *fmt;
  va_dcl
#else
NLM_EXTERN void CDECL  ValidErr (ValidStructPtr vsp, int severity, int code1, int code2, const char *fmt, ...)
#endif
{
    va_list  args;
	GatherContextPtr gcp;
	CharPtr tmp, ctmp;
	Int2 buflen, diff;
	BioseqPtr bsp;
	SeqIdPtr sip;
	SeqLocPtr loc = NULL;

	if (vsp == NULL || severity < vsp->cutoff) return;

	if (vsp->errbuf == NULL)
	{
		vsp->errbuf = MemNew(1024);
		if (vsp->errbuf == NULL)
			AbnormalExit(1);
	}
	tmp = vsp->errbuf;

	vsp->errors[severity]++;

#ifdef VAR_ARGS
    va_start (args);
#else
    va_start (args, fmt);
#endif

	gcp = vsp->gcp;
	buflen = 1023;
    vsprintf (tmp, fmt, args);
	while (*tmp != '\0')
	{
		buflen--;
		tmp++;
	}

	va_end (args);

	if (vsp->sfp != NULL)
	{
		diff = LabelCopy(tmp, " FEATURE: ", buflen);
		buflen -= diff;
		tmp += diff;

		diff = FeatDefLabel(vsp->sfp, tmp, buflen, OM_LABEL_BOTH);
		buflen -= diff;
		tmp += diff;

		if (vsp->suppressContext) {
			loc = AsnIoMemCopy(vsp->sfp->location,
								(AsnReadFunc) SeqLocAsnRead,
								(AsnWriteFunc) SeqLocAsnWrite);
			ChangeSeqLocToBestID (loc);
			ctmp = SeqLocPrint(loc);
			SeqLocFree (loc);
		} else {
			ctmp = SeqLocPrint(vsp->sfp->location);
		}
		if (ctmp != NULL)
		{
			diff = LabelCopyExtra(tmp, ctmp, buflen, " [", "]");
			buflen -= diff;
			tmp += diff;
			MemFree(ctmp);
		}

		if (! vsp->suppressContext) {
			sip = SeqLocId(vsp->sfp->location);
			if (sip != NULL)
			{
				bsp = BioseqFind(sip);
				if (bsp != NULL)
				{
					diff = LabelCopy(tmp, " [", buflen);
					buflen -= diff;
					tmp += diff;

					diff = BioseqLabel(bsp, tmp, buflen, OM_LABEL_BOTH);
					buflen -= diff;
					tmp += diff;

					diff = LabelCopy(tmp, "]", buflen);
					buflen -= diff;
					tmp += diff;
				}
			}
		}
		if (vsp->sfp->product != NULL) {
			if (vsp->suppressContext) {
				loc = AsnIoMemCopy(vsp->sfp->product,
									(AsnReadFunc) SeqLocAsnRead,
									(AsnWriteFunc) SeqLocAsnWrite);
				ChangeSeqLocToBestID (loc);
				ctmp = SeqLocPrint(loc);
				SeqLocFree (loc);
			} else {
				ctmp = SeqLocPrint(vsp->sfp->product);
			}
			if (ctmp != NULL)
			{
				diff = LabelCopyExtra(tmp, ctmp, buflen, " -> [", "]");
				buflen -= diff;
				tmp += diff;
				MemFree(ctmp);
			}
		}
	}
	else if (vsp->descr != NULL)
	{
		diff = LabelCopy(tmp, " DESCRIPTOR: ", buflen);
		buflen -= diff;
		tmp += diff;

		diff = SeqDescLabel(vsp->descr, tmp, buflen, OM_LABEL_BOTH);
		buflen -= diff;
		tmp += diff;
	}

	/*
	if (vsp->suppressContext)
	{
	}
	else */ if (vsp->sfp == NULL)    /* sfp adds its own context */
	{
		if (vsp->bsp != NULL)
		{
			diff = LabelCopy(tmp, " BIOSEQ: ", buflen);
			buflen -= diff;
			tmp += diff;

			if (vsp->suppressContext) {
				diff = WorstBioseqLabel(vsp->bsp, tmp, buflen, OM_LABEL_CONTENT);
			} else {
				diff = BioseqLabel(vsp->bsp, tmp, buflen, OM_LABEL_BOTH);
			}
			buflen -= diff;
			tmp += diff;
		}
		else if (vsp->bssp != NULL)
		{
			diff = LabelCopy(tmp, " BIOSEQ-SET: ", buflen);
			buflen -= diff;
			tmp += diff;
		
			if (vsp->suppressContext) {
				diff = BioseqSetLabel(vsp->bssp, tmp, buflen, OM_LABEL_CONTENT);
			} else {
				diff = BioseqSetLabel(vsp->bssp, tmp, buflen, OM_LABEL_BOTH);
			}
			buflen -= diff;
			tmp += diff;
		}
	}

	ErrPostItem ((ErrSev) (severity), code1, code2, "%s", vsp->errbuf);
	vsp->errbuf[0] = '\0';

	return;
}

/*****************************************************************************
*
*   Valid1GatherProc(gcp)
*     top level gather callback
*     dispatches to other levels
*
*****************************************************************************/
static Boolean Valid1GatherProc (GatherContextPtr gcp)
{
	ValidStructPtr vsp;
	SeqFeatPtr sfp;
	ValNodePtr sdp;
	BioSourcePtr biop;

	vsp = (ValidStructPtr)(gcp->userdata);
	vsp->gcp = gcp;     /* needed for ValidErr */

	switch (gcp->thistype)
	{
		case OBJ_BIOSEQ:
			if (! vsp->onlyspell) {
				ValidateBioseqInst (gcp);
				ValidateBioseqContext (gcp);
			}
			break;
		case OBJ_BIOSEQSET:
			if (! vsp->onlyspell) {
				ValidateBioseqSet (gcp);
			}
			break;
		case OBJ_SEQFEAT:
			if (! vsp->onlyspell) {
				ValidateSeqFeat (gcp);
				sfp = (SeqFeatPtr)(gcp->thisitem);
				if (sfp != NULL && sfp->data.choice == SEQFEAT_BIOSRC) {
					biop = (BioSourcePtr) sfp->data.value.ptrvalue;
					ValidateBioSource (vsp, gcp, biop);
				}
			}
			SpellCheckSeqFeat(gcp);
			break;
		case OBJ_SEQDESC:
			SpellCheckSeqDescr(gcp);
			/**
			ValidateSeqDescr (gcp);
		    **/
		    sdp = (ValNodePtr)(gcp->thisitem);
		    if (sdp != NULL && sdp->choice == Seq_descr_source) {
		    	biop = (BioSourcePtr) sdp->data.ptrvalue;
		    	ValidateBioSource (vsp, gcp, biop);
		    }
			break;
		default:
			break;
		
	}
	return TRUE;
}

static void LookForAnyPubAndOrg (SeqEntryPtr sep, BoolPtr no_pub, BoolPtr no_biosrc)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  SeqAnnotPtr   sap = NULL;
  ValNodePtr    sdp = NULL;
  SeqFeatPtr    sfp;
  SeqEntryPtr   tmp;

  if (sep == NULL || no_pub == NULL || no_biosrc == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp == NULL) return;
    sap = bsp->annot;
    sdp = bsp->descr;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;
    for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
      LookForAnyPubAndOrg (tmp, no_pub, no_biosrc);
    }
    sap = bssp->annot;
    sdp = bssp->descr;
  } else return;
  while (sap != NULL) {
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      while (sfp != NULL) {
        if (sfp->data.choice == SEQFEAT_PUB) {
          *no_pub = FALSE;
        } else if (sfp->data.choice == SEQFEAT_BIOSRC) {
          *no_biosrc = FALSE;
        }
        sfp = sfp->next;
      }
    }
    sap = sap->next;
  }
  while (sdp != NULL) {
    if (sdp->choice == Seq_descr_pub) {
      *no_pub = FALSE;
    } else if (sdp->choice == Seq_descr_source) {
      *no_biosrc = FALSE;
    }
    sdp = sdp->next;
  }
}

NLM_EXTERN Boolean ValidateSeqEntry(SeqEntryPtr sep, ValidStructPtr vsp)
{
	Uint2 entityID;
	GatherScope gs;
	BioseqSetPtr bssp;
	Boolean do_many = FALSE;
	Int2 errors[6], i;
	Boolean suppress_no_pubs = TRUE;
	Boolean suppress_no_biosrc = TRUE;
	GatherContextPtr gcp = NULL;
	GatherContext gc;
	SeqEntryPtr fsep;
	BioseqPtr fbsp = NULL;
	ObjMgrDataPtr omdp;
	SeqEntryPtr  oldsep;
	SeqEntryPtr  topsep;

	for (i =0; i < 6; i++)  /* keep errors between clears */
		errors[i] = 0;

	/* if no pubs or biosource, only one message, not one per bioseq */

	LookForAnyPubAndOrg (sep, &suppress_no_pubs, &suppress_no_biosrc);

	if (IS_Bioseq_set(sep))
	{
		bssp = (BioseqSetPtr)(sep->data.ptrvalue);
		switch (bssp->_class)
		{
			case BioseqseqSet_class_genbank:
			case BioseqseqSet_class_pir:
			case BioseqseqSet_class_gibb:
			case BioseqseqSet_class_gi:
			case BioseqseqSet_class_swissprot:
				sep = bssp->seq_set;
				do_many = TRUE;
				break;
			default:
				break;
		}
	}
					  
	globalvsp = vsp;   /* for spell checker */

	while (sep != NULL)
	{
		MemSet(&gs, 0, sizeof(GatherScope));
		gs.scope = sep;    /* default is to scope to this set */

		ValidStructClear(vsp);
      	vsp->sep = sep;

      	MemSet ((Pointer) &gc, 0, sizeof (GatherContext));
      	gcp = &gc;
      	gc.entityID = ObjMgrGetEntityIDForChoice (sep);
      	gc.itemID = 1;
      	if (IS_Bioseq (sep)) {
      		gc.thistype = OBJ_BIOSEQ;
      	} else {
      		gc.thistype = OBJ_BIOSEQSET;
      	}
     	vsp->gcp = gcp; /* above needed for ValidErr */
      	vsp->suppress_no_pubs = suppress_no_pubs;
      	vsp->suppress_no_biosrc = suppress_no_biosrc;
 
		/* build seqmgr feature indices if not already done */

		if (vsp->useSeqMgrIndexes) {
			entityID = ObjMgrGetEntityIDForChoice (sep);

			if (SeqMgrFeaturesAreIndexed (entityID) == 0) {
				SeqMgrIndexFeatures (entityID, NULL);
			}
		}

		fsep = FindNthBioseq (sep, 1);
		if (fsep != NULL && IS_Bioseq (fsep)) {
		  fbsp = (BioseqPtr) fsep->data.ptrvalue;
		  /* report context as first bioseq */
		  vsp->bsp = fbsp;
		}
      	if (suppress_no_pubs) {
      		omdp = ObjMgrGetData (gc.entityID);
      		if (omdp == NULL || omdp->datatype != OBJ_SEQSUB) {
				ValidErr(vsp, SEV_ERROR, ERR_SEQ_DESCR_NoPubFound, "No publications anywhere on this entire record.");
      		}
      	}
      	if (suppress_no_biosrc) {
			ValidErr(vsp, SEV_ERROR, ERR_SEQ_DESCR_NoOrgFound, "No organism name anywhere on this entire record.");
      	}
      	vsp->bsp = NULL;

		topsep = GetTopSeqEntryForEntityID (gc.entityID);
		oldsep = SeqEntrySetScope (topsep);

      	GatherSeqEntry(sep, (Pointer)vsp, Valid1GatherProc, &gs);

		SeqEntrySetScope (oldsep);

		if (do_many)
		{
			for (i = 0; i < 6; i++)
				errors[i] += vsp->errors[i];
			sep = sep->next;
		}
		else
			sep = NULL;
	}

	if (do_many)
	{
		for (i =0; i < 6; i++)
			vsp->errors[i] = errors[i];
	}
      
   return TRUE;
}

static void ValidateSetContents (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
    BioseqPtr bsp;
    ValidStructPtr vsp;
      
    vsp = (ValidStructPtr)data;
      
	if (IS_Bioseq(sep))
	{
		bsp = (BioseqPtr)(sep->data.ptrvalue);
		if (ISA_aa(bsp->mol))
			vsp->protcnt++;
		else
			vsp->nuccnt++;
		if (bsp->repr == Seq_repr_seg)
			vsp->segcnt++;

	}
	return;
}


static CharPtr GetBioseqSetClass(Uint1 cl)
{
    if(cl == BioseqseqSet_class_nuc_prot)
        return("nuc-prot");
    if(cl == BioseqseqSet_class_segset)
        return("segset");
    if(cl == BioseqseqSet_class_conset)
        return("conset");
    if(cl == BioseqseqSet_class_parts)
        return("parts");
    if(cl == BioseqseqSet_class_gibb)
        return("gibb");
    if(cl == BioseqseqSet_class_gi)
        return("gi");
    if(cl == BioseqseqSet_class_genbank)
        return("genbank");
    if(cl == BioseqseqSet_class_pir)
        return("pir");
    if(cl == BioseqseqSet_class_pub_set)
        return("pub-set");
    if(cl == BioseqseqSet_class_equiv)
        return("equiv");
    if(cl == BioseqseqSet_class_swissprot)
        return("swissprot");
    if(cl == BioseqseqSet_class_pdb_entry)
        return("pdb-entry");
    if(cl == BioseqseqSet_class_mut_set)
        return("mut-set");
    if(cl == BioseqseqSet_class_pop_set)
        return("pop-set");
    if(cl == BioseqseqSet_class_phy_set)
        return("phy-set");
    if(cl == BioseqseqSet_class_other)
        return("other");
    return("not-set");
}

static void ValidateNucProtSet(BioseqSetPtr bssp, ValidStructPtr vsp)
{
    SeqEntryPtr  sep;
    BioseqSetPtr bssp1;

    if(bssp->_class != BioseqseqSet_class_nuc_prot)
        return;

    for(sep = bssp->seq_set; sep != NULL; sep = sep->next)
    {
        if(!IS_Bioseq_set(sep))
            continue;

        bssp1 = sep->data.ptrvalue;
        if(bssp1 == NULL)
            continue;

        if(bssp1->_class != BioseqseqSet_class_segset)
        {
            ValidErr(vsp, SEV_FATAL, ERR_SEQ_PKG_NucProtNotSegSet,
                     "Nuc-prot Bioseq-set contains wrong Bioseq-set, its class is \"%s\".",
                     GetBioseqSetClass(bssp1->_class));
            break;
        }
    }
}

static void ValidateSegmentedSet(BioseqSetPtr bssp, ValidStructPtr vsp)
{
    SeqEntryPtr  sep;
    BioseqSetPtr bssp1;
    BioseqPtr    bsp;
    Uint1        mol = 0;

    if(bssp->_class != BioseqseqSet_class_segset)
        return;

    for(sep = bssp->seq_set; sep != NULL; sep = sep->next)
    {
        if (IS_Bioseq (sep)) {
        	bsp = (BioseqPtr) sep->data.ptrvalue;
        	if (bsp != NULL) {
        		if (mol == 0 || mol == Seq_mol_other) {
        			mol = bsp->mol;
        		} else if (bsp->mol != Seq_mol_other) {
        			if (ISA_na (bsp->mol) != ISA_na (mol)) {
        			    ValidErr(vsp, SEV_FATAL, ERR_SEQ_PKG_SegSetMixedBioseqs,
                  			   "Segmented set contains mixture of nucleotides and proteins");
        			}
        		}
        	}
        }

        if(!IS_Bioseq_set(sep))
            continue;

        bssp1 = sep->data.ptrvalue;
        if(bssp1 == NULL)
            continue;

        if(bssp1->_class != BioseqseqSet_class_parts)
        {
            ValidErr(vsp, SEV_FATAL, ERR_SEQ_PKG_SegSetNotParts,
                     "Segmented set contains wrong Bioseq-set, its class is \"%s\".",
                     GetBioseqSetClass(bssp1->_class));
            break;
        }
    }
}

static void ValidatePartsSet(BioseqSetPtr bssp, ValidStructPtr vsp)
{
    SeqEntryPtr  sep;
    BioseqSetPtr bssp1;
    BioseqPtr    bsp;
    Uint1        mol = 0;

    if(bssp->_class != BioseqseqSet_class_parts)
        return;

    for(sep = bssp->seq_set; sep != NULL; sep = sep->next)
    {
        if (IS_Bioseq (sep)) {
        	bsp = (BioseqPtr) sep->data.ptrvalue;
        	if (bsp != NULL) {
        		if (mol == 0 || mol == Seq_mol_other) {
        			mol = bsp->mol;
        		} else if (bsp->mol != Seq_mol_other) {
        			if (ISA_na (bsp->mol) != ISA_na (mol)) {
        			    ValidErr(vsp, SEV_FATAL, ERR_SEQ_PKG_PartsSetMixedBioseqs,
                  			   "Parts set contains mixture of nucleotides and proteins");
                  		break;
        			}
        		}
        	}
        }
    }

    for(sep = bssp->seq_set; sep != NULL; sep = sep->next)
    {
        if (IS_Bioseq_set(sep)) {
    	    bssp1 = sep->data.ptrvalue;
     	   if(bssp1 == NULL)
         	   continue;

        	 ValidErr(vsp, SEV_FATAL, ERR_SEQ_PKG_PartsSetHasSets,
                     "Parts set contains unwanted Bioseq-set, its class is \"%s\".",
                     GetBioseqSetClass(bssp1->_class));
            break;
        }
    }
}

static Boolean CheckForInconsistentBiosources (SeqEntryPtr sep, ValidStructPtr vsp, OrgRefPtr PNTR orpp)

{
	BioseqPtr bsp;
	BioseqSetPtr bssp;
	SeqEntryPtr tmp;
	ValNodePtr sdp;
	SeqFeatPtr sfp;
	SeqMgrDescContext dcontext;
	SeqMgrFeatContext fcontext;
	BioSourcePtr biop;
	OrgRefPtr orp;
	OrgRefPtr firstorp;
	GatherContextPtr gcp;
	Uint2 entityID, oldEntityID;
	Uint2 itemID, oldItemID;
	Uint2 itemtype, oldItemtype;

	if (sep == NULL || vsp == NULL || orpp == NULL) return FALSE;
	gcp = vsp->gcp;

	if (IS_Bioseq_set (sep)) {
		bssp = (BioseqSetPtr) sep->data.ptrvalue;
		if (bssp == NULL) return FALSE;
		for (tmp = bssp->seq_set; tmp != NULL; tmp = tmp->next) {
			if (CheckForInconsistentBiosources (tmp, vsp, orpp)) return TRUE;
		}
		return FALSE;
	}

	if (! IS_Bioseq (sep)) return FALSE;
	bsp = (BioseqPtr) sep->data.ptrvalue;
	if (bsp == NULL) return FALSE;

	biop = NULL;
	sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
	if (sdp != NULL) {
		biop = (BioSourcePtr) sdp->data.ptrvalue;
		entityID = dcontext.entityID;
		itemID = dcontext.itemID;
		itemtype = OBJ_SEQDESC;
	} else {
		sfp = SeqMgrGetNextFeature (bsp, NULL, SEQFEAT_BIOSRC, 0, &fcontext);
		if (sfp != NULL) {
			biop = (BioSourcePtr) sfp->data.value.ptrvalue;
			entityID = fcontext.entityID;
			itemID = fcontext.itemID;
			itemtype = OBJ_SEQFEAT;
		}
	}
	if (biop == NULL) return FALSE;
	orp = biop->org;
	if (orp == NULL) return FALSE;

	firstorp = *orpp;
	if (firstorp == NULL) {
	  *orpp = orp;
	  return FALSE;
	}

	if (StringNICmp (orp->taxname, "Influenza virus ", 16) == 0 &&
		StringNICmp (firstorp->taxname, "Influenza virus ", 16) == 0 &&
		StringNICmp (orp->taxname, firstorp->taxname, 17) == 0) {
		return FALSE;
	}

	if (StringICmp (orp->taxname, firstorp->taxname) == 0) return FALSE;

	oldEntityID = gcp->entityID;
	oldItemID = gcp->itemID;
	oldItemtype = gcp->thistype;

	gcp->entityID = entityID;
	gcp->itemID = itemID;
	gcp->thistype = itemtype;

	/* only report the first one that doesn't match */

	ValidErr (vsp, SEV_ERROR, ERR_SEQ_DESCR_InconsistentBioSources,
				"Population set contains inconsistent organisms.");

	gcp->entityID = oldEntityID;
	gcp->itemID = oldItemID;
	gcp->thistype = oldItemtype;

	return TRUE;
}

static void ValidatePopSet(BioseqSetPtr bssp, ValidStructPtr vsp)

{
	OrgRefPtr orp = NULL;
	SeqEntryPtr sep;

    if (bssp->_class != BioseqseqSet_class_pop_set) return;

    for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
    	if (CheckForInconsistentBiosources (sep, vsp, &orp)) return;
    }
}

static void ValidateBioseqSet (GatherContextPtr gcp)
{
	BioseqSetPtr bssp;
	ValidStructPtr vsp;
	SeqEntryPtr sep;

	vsp = (ValidStructPtr)(gcp->userdata);
	bssp = (BioseqSetPtr)(gcp->thisitem);
	vsp->bssp = bssp;
	vsp->bsp = NULL;
	vsp->descr = NULL;
	vsp->sfp = NULL;

	if (vsp->non_ascii_chars)  /* non_ascii chars in AsnRead step */
	{
		ValidErr(vsp, SEV_ERROR, ERR_GENERIC_NonAsciiAsn, "Non-ascii chars in input ASN.1 strings");
		vsp->non_ascii_chars = FALSE;   /* only do once */
	}

	vsp->nuccnt = 0;
	vsp->segcnt = 0;
	vsp->protcnt = 0;

	sep = gcp->sep;

	SeqEntryExplore(sep, (Pointer)vsp, ValidateSetContents);

	switch (bssp->_class)
	{
		case 1:     /* nuc-prot */
		    if (vsp->nuccnt == 0)
				ValidErr(vsp, SEV_ERROR, ERR_SEQ_PKG_NucProtProblem, "No nucleotides in nuc-prot set");
			if (vsp->protcnt == 0)
				ValidErr(vsp, SEV_ERROR, ERR_SEQ_PKG_NucProtProblem, "No proteins in nuc-prot set");
			ValidateNucProtSet(bssp, vsp);
			break;
		case 2:     /* seg set */
			if (vsp->segcnt == 0)
				ValidErr(vsp, SEV_ERROR, ERR_SEQ_PKG_SegSetProblem, "No segmented Bioseq in segset");
			ValidateSegmentedSet(bssp, vsp);
			break;
		case 4:     /* seg set */
			ValidatePartsSet(bssp, vsp);
			break;
		case BioseqseqSet_class_pop_set: /* population set */
			ValidatePopSet(bssp, vsp);
			break;
		default:
			if (! ((vsp->nuccnt) || (vsp->protcnt)))
				ValidErr(vsp, SEV_WARNING, ERR_SEQ_PKG_EmptySet, "No Bioseqs in this set");
			break;
	}
	return;
}

/*****************************************************************************
*
*   ValidateBioseqInst(gcp)
*      Validate one Bioseq Seq-inst
*
*****************************************************************************/
static void ValidateBioseqInst (GatherContextPtr gcp)
{
	Boolean retval = TRUE;
	Int2 i, start_at, num;
	Boolean errors[4], check_alphabet;
	static char * repr [8] = {
		"virtual", "raw", "segmented", "constructed",
		"reference", "consensus", "map", "delta" };
	SeqPortPtr spp;
	Int2 residue, x, termination;
	Int4 len, divisor = 1, len2;
	ValNode head;
	ValNodePtr vnp, vnp2;
	BioseqContextPtr bcp;
	Boolean got_partial, is_invalid;
	int seqtype, terminations;
	ValidStructPtr vsp;
	BioseqPtr bsp;
	SeqIdPtr sip1, sip2;
	Char buf1 [41], buf2 [41];
	SeqLitPtr slitp;
	SeqCodeTablePtr sctp;
	MolInfoPtr mip;
	Boolean litHasData;
	SeqMgrDescContext context;
	SeqFeatPtr cds;
	GeneRefPtr grp;
	SeqFeatPtr gene;
	SeqMgrFeatContext genectxt;
	CharPtr genelbl;
	SeqFeatPtr prot;
	SeqMgrFeatContext protctxt;
	CharPtr protlbl;

							/* set up data structures */

	vsp = (ValidStructPtr)(gcp->userdata);
	bsp = (BioseqPtr)(gcp->thisitem);
	vsp->bsp = bsp;
	vsp->descr = NULL;
	vsp->sfp = NULL;
	vsp->bssp = (BioseqSetPtr)(gcp->parentitem);
	vsp->bsp_partial_val = 0;

	if (vsp->non_ascii_chars)  /* non_ascii chars in AsnRead step */
	{
		ValidErr(vsp, SEV_FATAL, ERR_GENERIC_NonAsciiAsn, "Non-ascii chars in input ASN.1 strings");
		vsp->non_ascii_chars = FALSE;   /* only do once */
	}

	if (bsp->id == NULL)
	{
		ValidErr(vsp, SEV_FATAL, ERR_SEQ_INST_NoIdOnBioseq, "No ids on a Bioseq");
		return;
	}

	for (sip1 = bsp->id; sip1 != NULL; sip1 = sip1->next) {
		for (sip2 = sip1->next; sip2 != NULL; sip2 = sip2->next) {
			if (SeqIdComp (sip1, sip2) != SIC_DIFF) {
				SeqIdWrite (sip1, buf1, PRINTID_FASTA_SHORT, 40);
				SeqIdWrite (sip2, buf2, PRINTID_FASTA_SHORT, 40);
				ValidErr(vsp, SEV_ERROR, ERR_SEQ_INST_ConflictingIdsOnBioseq,
					"Conflicting ids on a Bioseq: (%s - %s)", buf1, buf2);
			}
		}
	}

	for (i = 0; i < 4; i++)
		errors[i] = FALSE;

	switch (bsp->repr)
	{
		case Seq_repr_virtual:
			if ((bsp->seq_ext_type) || (bsp->seq_ext != NULL))
				errors[0] = TRUE;
			if ((bsp->seq_data_type) || (bsp->seq_data != NULL))
				errors[3] = TRUE;
			break;
		case Seq_repr_map:
			if ((bsp->seq_ext_type != 3) || (bsp->seq_ext == NULL))
				errors[1] = TRUE;
			if ((bsp->seq_data_type) || (bsp->seq_data != NULL))
				errors[3] = TRUE;
			break;
		case Seq_repr_ref:
			if ((bsp->seq_ext_type != 2) || (bsp->seq_ext == NULL))
				errors[1] = TRUE;
			if ((bsp->seq_data_type) || (bsp->seq_data != NULL))
				errors[3] = TRUE;
			break;
		case Seq_repr_seg:
			if ((bsp->seq_ext_type != 1) || (bsp->seq_ext == NULL))
				errors[1] = TRUE;
			if ((bsp->seq_data_type) || (bsp->seq_data != NULL))
				errors[3] = TRUE;
			break;
		case Seq_repr_raw:
		case Seq_repr_const:
			if ((bsp->seq_ext_type) || (bsp->seq_ext != NULL))
				errors[0] = TRUE;
			if ((bsp->seq_data_type < 1) || (bsp->seq_data_type > 11)
				|| (bsp->seq_data == NULL))
				errors[2] = TRUE;
			break;
		case Seq_repr_delta:
			if ((bsp->seq_ext_type != 4) || (bsp->seq_ext == NULL))
				errors[1] = TRUE;
			if ((bsp->seq_data_type) || (bsp->seq_data != NULL))
				errors[3] = TRUE;
			break;
		default:
			ValidErr(vsp, SEV_FATAL, ERR_SEQ_INST_ReprInvalid, "Invalid Bioseq->repr = %d", (int)(bsp->repr));
			return;
	}

	if (errors[0] == TRUE)
	{
		ValidErr(vsp, SEV_ERROR, ERR_SEQ_INST_ExtNotAllowed, "Bioseq-ext not allowed on %s Bioseq", repr[bsp->repr - 1]);
		retval = FALSE;
	}

	if (errors[1] == TRUE)
	{
		ValidErr(vsp, SEV_ERROR, ERR_SEQ_INST_ExtBadOrMissing, "Missing or incorrect Bioseq-ext on %s Bioseq", repr[bsp->repr - 1]);
		retval = FALSE;
	}

	if (errors[2] == TRUE)
	{
		ValidErr(vsp, SEV_ERROR, ERR_SEQ_INST_SeqDataNotFound, "Missing Seq-data on %s Bioseq", repr[bsp->repr - 1]);
		retval = FALSE;
	}

	if (errors[3] == TRUE)
	{
		ValidErr(vsp, SEV_ERROR, ERR_SEQ_INST_SeqDataNotAllowed, "Seq-data not allowed on %s Bioseq", repr[bsp->repr - 1]);
		retval = FALSE;
	}

	if (! retval) return;

	if (ISA_aa(bsp->mol))
	{
		if (bsp->topology > 1)   /* not linear */
		{
			ValidErr(vsp, SEV_ERROR, ERR_SEQ_INST_CircularProtein, "Non-linear topology set on protein");
		}
		if (bsp->strand > 1)
		{
			ValidErr(vsp, SEV_ERROR, ERR_SEQ_INST_DSProtein, "Protein not single stranded");
		}

	}
	else
	{
		if (! bsp->mol)
			ValidErr(vsp, SEV_ERROR, ERR_SEQ_INST_MolNotSet, "Bioseq.mol is 0");
		else if (bsp->mol == Seq_mol_other)
			ValidErr(vsp, SEV_ERROR, ERR_SEQ_INST_MolOther, "Bioseq.mol is type other");
	}
	                                          /* check sequence alphabet */
	if ((bsp->repr == Seq_repr_raw) || (bsp->repr == Seq_repr_const))
	{
		if (bsp->fuzz != NULL)
		{
			ValidErr(vsp, SEV_ERROR, ERR_SEQ_INST_FuzzyLen, "Fuzzy length on %s Bioseq", repr[bsp->repr - 1]);
		}

		if (bsp->length < 1)
		{
			ValidErr(vsp, SEV_ERROR, ERR_SEQ_INST_InvalidLen, "Invalid Bioseq length [%ld]", (long)bsp->length);
		}

		seqtype = (int)(bsp->seq_data_type);
		switch (seqtype)
		{
			case Seq_code_iupacna:
			case Seq_code_ncbi2na:
			case Seq_code_ncbi4na:
			case Seq_code_ncbi8na:
			case Seq_code_ncbipna:
				if (ISA_aa(bsp->mol))
				{
					ValidErr(vsp, SEV_FATAL, ERR_SEQ_INST_InvalidAlphabet,"Using a nucleic acid alphabet on a protein sequence");
					return;
				}
				break;
			case Seq_code_iupacaa:
			case Seq_code_ncbi8aa:
			case Seq_code_ncbieaa:
			case Seq_code_ncbipaa:
			case Seq_code_ncbistdaa:
				if (ISA_na(bsp->mol))
				{
					ValidErr(vsp, SEV_FATAL, ERR_SEQ_INST_InvalidAlphabet,"Using a protein alphabet on a nucleic acid");
					return;
				}
				break;
			default:
				ValidErr(vsp, SEV_FATAL, ERR_SEQ_INST_InvalidAlphabet, "Using illegal sequence alphabet [%d]",(int)bsp->seq_data_type);
				return;
		}

		check_alphabet = FALSE;
		switch (seqtype)
		{
			case Seq_code_iupacaa:
			case Seq_code_iupacna:
			case Seq_code_ncbieaa:
			case Seq_code_ncbistdaa:
				check_alphabet = TRUE;

			case Seq_code_ncbi8na:
			case Seq_code_ncbi8aa:
				divisor = 1;
				break;

			case Seq_code_ncbi4na:
				divisor = 2;
				break;

			case Seq_code_ncbi2na:
				divisor = 4;
				break;

			case Seq_code_ncbipna:
				divisor = 5;
				break;

			case Seq_code_ncbipaa:
				divisor = 21;
				break;
		}

		len = bsp->length;
		if (len % divisor) len += divisor;
		len /= divisor;
		len2 = BSLen(bsp->seq_data);
		if (len > len2)
		{
			ValidErr(vsp, SEV_FATAL, ERR_SEQ_INST_SeqDataLenWrong,"Bioseq.seq_data too short [%ld] for given length [%ld]",	(long)(len2 * divisor), (long)bsp->length);
			return;
		}
		else if (len < len2)
		{
			ValidErr(vsp, SEV_FATAL, ERR_SEQ_INST_SeqDataLenWrong,"Bioseq.seq_data is larger [%ld] than given length [%ld]",(long)(len2 * divisor), (long)bsp->length);
		}

		if (check_alphabet)				  /* check 1 letter alphabets */
		{
			switch (seqtype)
			{
				case Seq_code_iupacaa:
				case Seq_code_ncbieaa:
					termination = '*';
					break;
				case Seq_code_ncbistdaa:
					termination = 25;
					break;
                default:
                    termination = '\0';
                    break;
			}
			spp = SeqPortNew(bsp, 0, -1, 0, 0);
			if (spp == NULL)
			{
				ValidErr(vsp, SEV_FATAL, ERR_SEQ_INST_SeqPortFail, "Can't open SeqPort");
				return;
			}

			i = 0;
			terminations = 0;
			for (len=0; len < bsp->length; len++)
			{
				residue = SeqPortGetResidue(spp);
				if (! IS_residue(residue))
				{
					i++;
					if (i > 10)
					{
						ValidErr(vsp, SEV_FATAL,ERR_SEQ_INST_InvalidResidue,"More than 10 invalid residues. Checking stopped");
						SeqPortFree(spp);
						if (vsp->patch_seq)
							PatchBadSequence(bsp);
						return;
			 		}
					else
					{
						BSSeek(bsp->seq_data, len, SEEK_SET);
						x = BSGetByte(bsp->seq_data);
						if (bsp->seq_data_type == Seq_code_ncbistdaa)
							ValidErr(vsp, SEV_FATAL, ERR_SEQ_INST_InvalidResidue, "Invalid residue [%d] in position [%ld]",	(int)x, (long)(len+1));
						else
							ValidErr(vsp, SEV_FATAL, ERR_SEQ_INST_InvalidResidue, "Invalid residue [%c] in position [%ld]",	(char)x, (long)(len+1));
					}
				}
				else if (residue == termination)
					terminations++;
			}
			SeqPortFree(spp);
			if (terminations)
			{
				cds = SeqMgrGetCDSgivenProduct (bsp, NULL);
				grp = SeqMgrGetGeneXref (cds);
				genelbl = NULL;
				if (grp == NULL && cds != NULL) {
				  gene = SeqMgrGetOverlappingGene (cds->location, &genectxt);
				  if (gene != NULL) {
				    grp = (GeneRefPtr) gene->data.value.ptrvalue;
				  }
				}
				if (grp != NULL && (! SeqMgrGeneIsSuppressed (grp))) {
					if (grp->locus != NULL)
						genelbl = (grp->locus);
					else if (grp->desc != NULL)
						genelbl = (grp->desc);
					else if (grp->syn != NULL)
						genelbl = (CharPtr)(grp->syn->data.ptrvalue);
				}
				prot = SeqMgrGetBestProteinFeature (bsp, &protctxt);
				protlbl = protctxt.label;
				if (StringHasNoText (genelbl)) {
					genelbl = "gene?";
				}
				if (StringHasNoText (protlbl)) {
					protlbl = "prot?";
				}
				ValidErr(vsp, SEV_ERROR, ERR_SEQ_INST_StopInProtein,"[%d] termination symbols in protein sequence (%s - %s)", terminations, genelbl, protlbl);
				if (! i)
					return;
			}
			if (i)
			{
				if (vsp->patch_seq)
					PatchBadSequence(bsp);
				return;
			}

		}
	}

	if ((bsp->repr == Seq_repr_seg) ||
		(bsp->repr == Seq_repr_ref))/* check segmented sequence */
	{
		head.choice = SEQLOC_MIX;
		head.data.ptrvalue = bsp->seq_ext;
		head.next = NULL;
		ValidateSeqLoc(vsp, (SeqLocPtr)&head, "Segmented Bioseq");
		                            /* check the length */
		len = 0;
		vnp = NULL;
		while ((vnp = SeqLocFindNext(&head, vnp)) != NULL)
		{
			len2 = SeqLocLen(vnp);
			if (len2 > 0)
				len += len2;
		}
		if (bsp->length > len)
		{
			ValidErr(vsp, SEV_FATAL, ERR_SEQ_INST_SeqDataLenWrong,"Bioseq.seq_data too short [%ld] for given length [%ld]",	(long)(len), (long)bsp->length);
		}
		else if (bsp->length < len)
		{
			ValidErr(vsp, SEV_FATAL, ERR_SEQ_INST_SeqDataLenWrong,"Bioseq.seq_data is larger [%ld] than given length [%ld]",(long)(len), (long)bsp->length);
		}
		

		vsp->bsp_partial_val = SeqLocPartialCheck((SeqLocPtr)(&head));
		if ((vsp->bsp_partial_val) && (ISA_aa(bsp->mol)))
		{
			bcp = NULL;
			vnp = NULL;
			got_partial = FALSE;
			if (vsp->useSeqMgrIndexes) {
			  vnp = SeqMgrGetNextDescriptor(bsp, vnp, Seq_descr_modif, &context);
			} else {
			  bcp = BioseqContextNew(bsp);
			  vnp = BioseqContextGetSeqDescr(bcp, Seq_descr_modif, vnp, NULL);
			}
			while (vnp != NULL)
			{
				for (vnp2 = (ValNodePtr)(vnp->data.ptrvalue); vnp2 != NULL; vnp2 = vnp2->next)
				{
					switch(vnp2->data.intvalue)
					{
						case 10:   /* partial */
							got_partial = TRUE;
							break;
						case 16:   /* no-left */
							if (! (vsp->bsp_partial_val & SLP_START))
								ValidErr(vsp, SEV_ERROR, ERR_SEQ_INST_PartialInconsistent, "GIBB-mod no-left inconsistent with segmented SeqLoc");
							got_partial = TRUE;
							break;
						case 17:   /* no-right */
							if (! (vsp->bsp_partial_val & SLP_STOP))
								ValidErr(vsp, SEV_ERROR, ERR_SEQ_INST_PartialInconsistent, "GIBB-mod no-right inconsistent with segmented SeqLoc");
							got_partial = TRUE;
							break;

					}
				}
				if (vsp->useSeqMgrIndexes) {
				  vnp = SeqMgrGetNextDescriptor(bsp, vnp, Seq_descr_modif, &context);
				} else {
				  vnp = BioseqContextGetSeqDescr(bcp, Seq_descr_modif, vnp, NULL);
				}
			}
			if (! vsp->useSeqMgrIndexes) {
			  	BioseqContextFree(bcp);
			}
			if (! got_partial)
				ValidErr(vsp, SEV_ERROR, ERR_SEQ_INST_PartialInconsistent, "Partial segmented sequence without GIBB-mod");
		}
	}

	mip = NULL;

	if (bsp->repr == Seq_repr_delta)
	{
		len = 0;
		for (vnp = (ValNodePtr)(bsp->seq_ext); vnp != NULL; vnp = vnp->next)
		{
			if (vnp->data.ptrvalue == NULL)
				ValidErr(vsp, SEV_ERROR, ERR_SEQ_INST_SeqDataLenWrong, "NULL pointer in delta seq_ext valnode");
			else
			{
				switch (vnp->choice)
				{
					case 1:		 /* SeqLocPtr */
						len2 = SeqLocLen((SeqLocPtr)(vnp->data.ptrvalue));
						if (len2 < 0)
							  ValidErr(vsp, SEV_ERROR, ERR_SEQ_INST_SeqDataLenWrong, "-1 length on seq-loc of delta seq_ext");
						else
							len += len2;
						break;
					case 2:     /* SeqLitPtr */
						slitp = (SeqLitPtr)(vnp->data.ptrvalue);
						if (slitp->seq_data != NULL)
						{
							sctp = SeqCodeTableFind(slitp->seq_data_type);
							if (sctp == NULL)
							{
                                ValidErr(vsp, SEV_FATAL, ERR_SEQ_INST_InvalidAlphabet, "Using illegal sequence alphabet [%d] in SeqLitPtr",
									(int)slitp->seq_data_type);
								len += slitp->length;
								break;
							}

							start_at = (Int2)(sctp->start_at);
							num = (Int2)(sctp->num);

							switch (slitp->seq_data_type)
							{
								case Seq_code_iupacaa:
								case Seq_code_iupacna:
								case Seq_code_ncbieaa:
								case Seq_code_ncbistdaa:
									BSSeek(slitp->seq_data, 0, SEEK_SET);
									for (len2 = 1; len2 <= (slitp->length); len2++)
									{
										is_invalid = FALSE;
										residue = BSGetByte(slitp->seq_data);
										i = residue - start_at;
										if ((i < 0) || (i >= num))
											is_invalid = TRUE;
										else if (*(sctp->names[i]) == '\0')
											is_invalid = TRUE;
										if (is_invalid)
										{
											if (slitp->seq_data_type == Seq_code_ncbistdaa)
												ValidErr(vsp, SEV_FATAL, ERR_SEQ_INST_InvalidResidue, "Invalid residue [%d] in position [%ld]",	(int)residue, (long)(len+len2));
											else
												ValidErr(vsp, SEV_FATAL, ERR_SEQ_INST_InvalidResidue, "Invalid residue [%c] in position [%ld]",	(char)residue, (long)(len+len2));
										}
									}
									break;
								default:
									break;
							}
						}
						len += slitp->length;
						break;
					default:
						ValidErr(vsp, SEV_ERROR, ERR_SEQ_INST_ExtNotAllowed, "Illegal choice [%d] in delta chain",
							(int)(vnp->choice));
						break;
				}
			}
		}
		if (bsp->length > len)
		{
			ValidErr(vsp, SEV_FATAL, ERR_SEQ_INST_SeqDataLenWrong,"Bioseq.seq_data too short [%ld] for given length [%ld]",	(long)(len), (long)bsp->length);
		}
		else if (bsp->length < len)
		{
			ValidErr(vsp, SEV_FATAL, ERR_SEQ_INST_SeqDataLenWrong,"Bioseq.seq_data is larger [%ld] than given length [%ld]",(long)(len), (long)bsp->length);
		}
		vnp = NULL;
		if (vsp->useSeqMgrIndexes) {
		  vnp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &context);
		} else {
		  bcp = BioseqContextNew(bsp);
		  vnp = BioseqContextGetSeqDescr (bcp, Seq_descr_molinfo, NULL, NULL);
		  BioseqContextFree (bcp);
		}
		if (vnp != NULL) {
			mip = (MolInfoPtr) vnp->data.ptrvalue;
			if (mip != NULL) {
				if (mip->tech != MI_TECH_htgs_0 && mip->tech != MI_TECH_htgs_1 && mip->tech != MI_TECH_htgs_2) {
					ValidErr(vsp, SEV_ERROR, ERR_SEQ_INST_BadDeltaSeq,"Delta seq technique should not be [%d]",(int)(mip->tech));
				}
			}
		}
	}

	if (ISA_aa(bsp->mol))
	{
		if ((bsp->length <= 3) && (bsp->length >= 0))
		{
			ValidErr(vsp, SEV_WARNING, ERR_SEQ_INST_ShortSeq, "Sequence only %ld residues", (long)(bsp->length));
		}

	}
	else
	{
		if ((bsp->length <= 10) && (bsp->length >= 0))
		{
			ValidErr(vsp, SEV_WARNING, ERR_SEQ_INST_ShortSeq, "Sequence only %ld residues", (long)(bsp->length));
		}
	}

	if (bsp->length > 350000)
	{
		if (bsp->repr == Seq_repr_delta) {
			if (mip != NULL) {
				if (mip->tech == MI_TECH_htgs_0 || mip->tech == MI_TECH_htgs_1 || mip->tech == MI_TECH_htgs_2) {
					ValidErr(vsp, SEV_WARNING, ERR_SEQ_INST_LongHtgsSequence,"Phase 0, 1 or 2 HTGS sequence exceeds 350kbp limit");
				} else if (mip->tech == MI_TECH_htgs_3) {
					ValidErr(vsp, SEV_WARNING, ERR_SEQ_INST_SequenceExceeds350kbp,"Phase 3 HTGS sequence exceeds 350kbp limit");
				} else {
					len = 0;
					litHasData = FALSE;
					for (vnp = (ValNodePtr)(bsp->seq_ext); vnp != NULL; vnp = vnp->next) {
						if (vnp->choice == 2) {
							slitp = (SeqLitPtr)(vnp->data.ptrvalue);
							if (slitp != NULL) {
								if (slitp->seq_data != NULL) {
									litHasData = TRUE;
								}
								len += slitp->length;
							}
						}
					}
					if (len > 350000 && litHasData) {
						ValidErr(vsp, SEV_FATAL, ERR_SEQ_INST_LongLiteralSequence,"Length of sequence literals exceeds 350kbp limit");
					}
				}
			}
		} else if (bsp->repr == Seq_repr_raw) {
			vnp = NULL;
			if (vsp->useSeqMgrIndexes) {
				vnp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &context);
			} else {
				bcp = BioseqContextNew (bsp);
				vnp = BioseqContextGetSeqDescr (bcp, Seq_descr_molinfo, NULL, NULL);
				BioseqContextFree (bcp);
			}
			if (vnp != NULL) {
				mip = (MolInfoPtr) vnp->data.ptrvalue;
			}
			if (mip != NULL) {
				if (mip->tech == MI_TECH_htgs_0 || mip->tech == MI_TECH_htgs_1 || mip->tech == MI_TECH_htgs_2) {
					ValidErr(vsp, SEV_WARNING, ERR_SEQ_INST_LongHtgsSequence,"Phase 0, 1 or 2 HTGS sequence exceeds 350kbp limit");
				} else if (mip->tech == MI_TECH_htgs_3) {
					ValidErr(vsp, SEV_WARNING, ERR_SEQ_INST_SequenceExceeds350kbp,"Phase 3 HTGS sequence exceeds 350kbp limit");
				} else {
					ValidErr(vsp, SEV_WARNING, ERR_SEQ_INST_SequenceExceeds350kbp,"Length of sequence exceeds 350kbp limit");
				}
			} else {
				ValidErr(vsp, SEV_WARNING, ERR_SEQ_INST_SequenceExceeds350kbp,"Length of sequence exceeds 350kbp limit");
			}
		} else {
			/* Could be a segset header bioseq that is > 350kbp */
			/* No-op for now? Or generate a warning? */
		}
	}

	return;
}

/*****************************************************************************
*
*   ValidatePubdesc(gcp)
*      Check pubdesc for missing information
*
*****************************************************************************/
static Boolean HasNoText (CharPtr str)

{
  Char  ch;

  if (str != NULL) {
    ch = *str;
    while (ch != '\0') {
      if (ch > ' ') {
        return FALSE;
      }
      str++;
      ch = *str;
    }
  }
  return TRUE;
}

static Boolean HasNoName (ValNodePtr name)

{
	AuthorPtr  ap;
	NameStdPtr  nsp;
	PersonIdPtr  pid;

	if (name != NULL) {
		ap = name->data.ptrvalue;
		if (ap != NULL) {
			pid = ap->name;
			if (pid != NULL) {
				if (pid->choice == 2) {
					nsp = pid->data;
					if (nsp != NULL) {
						if (! HasNoText (nsp->names [0])) {
							return FALSE;
						}
					}
				}
			}
		}
	}
	return TRUE;
}

static void ValidatePubdesc (ValidStructPtr vsp, PubdescPtr pdp)

{
	AuthListPtr  alp;
	CitArtPtr  cap;
	Boolean  hasName, hasTitle;
	ValNodePtr  name;
	ValNodePtr  title;
	ValNodePtr  vnp;

	if (vsp == NULL || pdp == NULL) return;
	for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
		switch (vnp->choice) {
			case PUB_Article :
				cap = (CitArtPtr) vnp->data.ptrvalue;
				hasName = FALSE;
				hasTitle = FALSE;
				if (cap != NULL) {
					for (title = cap->title; title != NULL; title = title->next) {
						if (! HasNoText ((CharPtr) title->data.ptrvalue)) {
							hasTitle = TRUE;
						}
					}
					if (! hasTitle) {
						ValidErr (vsp, SEV_ERROR, ERR_GENERIC_MissingPubInfo,
							"Publication has no title");
					}
					alp = cap->authors;
					if (alp != NULL) {
						if (alp->choice == 1) {
							for (name = alp->names; name != NULL; name = name->next) {
								if (! HasNoName (name)) {
									hasName = TRUE;
								}
							}
						} else if (alp->choice == 2 || alp->choice == 3) {
							for (name = alp->names; name != NULL; name = name->next) {
								if (! HasNoText ((CharPtr) name->data.ptrvalue)) {
									hasName = TRUE;
								}
							}
						}
					}
					if (! hasName) {
						ValidErr (vsp, SEV_ERROR, ERR_GENERIC_MissingPubInfo,
							"Publication has no author names");
					}
				}
				break;
			default :
				break;
		}
	}
}

typedef struct bioseqvalid {
	ValidStructPtr vsp;
	Boolean is_aa;                 /* bioseq is protein? */
	Boolean is_mrna;               /* molinfo is mrna? */
	Boolean is_prerna;             /* molinfo is precursor rna? */
	Boolean got_a_pub;
	int last_na_mol,
		last_na_mod,
		last_organelle,
		last_partialness,
		last_left_right,
		last_biomol,
		last_tech,
		last_completeness,
		num_full_length_src_feat,   /* number full length src feats */
	    num_full_length_prot_ref;
	ValNodePtr last_gb,
		last_embl,
		last_prf,
		last_pir,
		last_sp,
		last_pdb,
		last_create,
		last_update,
		last_biosrc,
		last_orgref;
	OrgRefPtr last_org;
	GatherContextPtr gcp;
} BioseqValidStr, PNTR BioseqValidStrPtr;

/*****************************************************************************
*
*   ValidateSeqFeatContext(gcp)
*      Gather callback helper function for validating context on a Bioseq
*
*****************************************************************************/
static Boolean ValidateSeqFeatCommon (SeqFeatPtr sfp, BioseqValidStrPtr bvsp, ValidStructPtr vsp,
                                      Int4 left, Int4 right, Uint2 featitemid, Boolean farloc)
{
	GatherContextPtr gcp = NULL;
	ImpFeatPtr ifp;
	Uint2 olditemtype;
	Uint2 olditemid;
	RnaRefPtr rrp;

	vsp->descr = NULL;
	vsp->sfp = sfp;

	if (featitemid > 0) {
		gcp = vsp->gcp;
		if (gcp != NULL) {
			olditemid = gcp->itemID;
			olditemtype = gcp->thistype;
			gcp->itemID = featitemid;
			gcp->thistype = OBJ_SEQFEAT;
		}
	}

	if (bvsp->is_aa)
	{
		if (sfp->data.choice == SEQFEAT_PROT)
		{
			if ((left == 0) &&
				(right == ((vsp->bsp->length) -1 )))
				bvsp->num_full_length_prot_ref++;
		}

		switch (sfp->data.choice)
		{
			case SEQFEAT_CDREGION:
			case SEQFEAT_RNA:
			case SEQFEAT_RSITE:
			case SEQFEAT_TXINIT:
				ValidErr(vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidForType, "Invalid feature for a protein Bioseq.");
				break;
			default:
				break;
		}

	}
	else
	{
		switch (sfp->data.choice)
		{
			case SEQFEAT_PROT:
			case SEQFEAT_PSEC_STR:
				ValidErr(vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidForType, "Invalid feature for a nucleotide Bioseq.");
				break;
			default:
				break;
		}

	}

	if (bvsp->is_mrna) {
		switch (sfp->data.choice) {
			case SEQFEAT_RNA:
				rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
				if (rrp != NULL && rrp->type == 2) {
					ValidErr(vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidForType, "mRNA feature is invalid on an mRNA (cDNA) Bioseq.");
				}
				break;
			case SEQFEAT_IMP:
				ifp = (ImpFeatPtr) sfp->data.value.ptrvalue;
				if (ifp != NULL && ifp->key != NULL && (! HasNoText (ifp->key))) {
					if (StringCmp (ifp->key, "intron") == 0 ||
						StringCmp (ifp->key, "CAAT_signal") == 0) {
						ValidErr(vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidForType, "Invalid feature for an mRNA Bioseq.");
					}
				}
				break;
			default:
				break;
		}
	} else if (bvsp->is_prerna) {
		switch (sfp->data.choice) {
			case SEQFEAT_IMP:
				ifp = (ImpFeatPtr) sfp->data.value.ptrvalue;
				if (ifp != NULL && ifp->key != NULL && (! HasNoText (ifp->key))) {
					if (StringCmp (ifp->key, "CAAT_signal") == 0) {
						ValidErr(vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidForType, "Invalid feature for an pre-RNA Bioseq.");
					}
				}
				break;
			default:
				break;
		}
	}

	if (farloc) {
		ValidErr(vsp, SEV_WARNING, ERR_SEQ_FEAT_FarLocation, "Feature has 'far' location - accession not packaged in record");
	}

	if ((sfp->data.choice == SEQFEAT_PUB) ||
		(sfp->cit != NULL))
		bvsp->got_a_pub = TRUE;


	if (gcp != NULL) {
		gcp->itemID = olditemid;
		gcp->thistype = olditemtype;
	}

	return TRUE;
}

static Boolean LIBCALLBACK ValidateSeqFeatIndexed (SeqFeatPtr sfp, SeqMgrFeatContextPtr context)
{
	ValidStructPtr vsp;
	BioseqValidStrPtr bvsp;

	bvsp = (BioseqValidStrPtr) context->userdata;
	vsp = bvsp->vsp;

	return ValidateSeqFeatCommon (sfp, bvsp, vsp, context->left, context->right, context->itemID, context->farloc);
}

static void ValidateSeqFeatContext (GatherContextPtr gcp)
{
	ValidStructPtr vsp;
	BioseqValidStrPtr bvsp;
	SeqFeatPtr sfp;

	bvsp = (BioseqValidStrPtr)(gcp->userdata);
	vsp = bvsp->vsp;
	sfp = (SeqFeatPtr)(gcp->thisitem);

	ValidateSeqFeatCommon (sfp, bvsp, vsp, gcp->extremes.left, gcp->extremes.right, 0, FALSE);
}

/*****************************************************************************
*
*   CountryIsValid(name)
*      Validates subsource country against official country names
*
*****************************************************************************/

static CharPtr countrycodes [] = {
  "Afghanistan",
  "Albania",
  "Algeria",
  "American Samoa",
  "Andorra",
  "Angola",
  "Anguilla",
  "Antarctica",
  "Antigua and Barbuda",
  "Argentina",
  "Armenia",
  "Aruba",
  "Ashmore and Cartier Islands",
  "Australia",
  "Austria",
  "Azerbaijan",
  "Bahamas",
  "Bahrain",
  "Baker Island",
  "Bangladesh",
  "Barbados",
  "Bassas da India",
  "Belarus",
  "Belgium",
  "Belize",
  "Benin",
  "Bermuda",
  "Bhutan",
  "Bolivia",
  "Bosnia and Herzegovina",
  "Botswana",
  "Bouvet Island",
  "Brazil",
  "British Virgin Islands",
  "Brunei",
  "Bulgaria",
  "Burkina Faso",
  "Burma",
  "Burundi",
  "Cambodia",
  "Cameroon",
  "Canada",
  "Cape Verde",
  "Cayman Islands",
  "Central African Republic",
  "Chad",
  "Chile",
  "China",
  "Christmas Island",
  "Clipperton Island",
  "Cocos Islands",
  "Colombia",
  "Comoros",
  "Cook Islands",
  "Coral Sea Islands",
  "Costa Rica",
  "Cote d' Ivoire",
  "Croatia",
  "Cuba",
  "Cyprus",
  "Czech Republic",
  "Democratic Republic of the Congo",
  "Denmark",
  "Djibouti",
  "Dominica",
  "Dominican Republic",
  "Ecuador",
  "Egypt",
  "El Salvador",
  "Equatorial Guinea",
  "Eritrea",
  "Estonia",
  "Ethiopia",
  "Europa Island",
  "Falkland Islands (Islas Malvinas)",
  "Faroe Islands",
  "Fiji",
  "Finland",
  "France",
  "French Guiana",
  "French Polynesia",
  "French Southern and Antarctic Lands",
  "Gabon",
  "Galapagos Islands",
  "Gambia",
  "Gaza Strip",
  "Georgia",
  "Germany",
  "Ghana",
  "Gibraltar",
  "Glorioso Islands",
  "Greece",
  "Greenland",
  "Grenada",
  "Guadeloupe",
  "Guam",
  "Guatemala",
  "Guernsey",
  "Guinea",
  "Guinea-Bissau",
  "Guyana",
  "Haiti",
  "Heard Island and McDonald Islands",
  "Honduras",
  "Hong Kong",
  "Howland Island",
  "Hungary",
  "Iceland",
  "India",
  "Indonesia",
  "Iran",
  "Iraq",
  "Ireland",
  "Isle of Man",
  "Israel",
  "Italy",
  "Jamaica",
  "Jan Mayen",
  "Japan",
  "Jarvis Island",
  "Jersey",
  "Johnston Atoll",
  "Jordan",
  "Juan de Nova Island",
  "Kazakhstan",
  "Kenya",
  "Kingman Reef",
  "Kiribati",
  "Kuwait",
  "Kyrgyzstan",
  "Laos",
  "Latvia",
  "Lebanon",
  "Lesotho",
  "Liberia",
  "Libya",
  "Liechtenstein",
  "Lithuania",
  "Luxembourg",
  "Macau",
  "Macedonia",
  "Madagascar",
  "Malawi",
  "Malaysia",
  "Maldives",
  "Mali",
  "Malta",
  "Marshall Islands",
  "Martinique",
  "Mauritania",
  "Mauritius",
  "Mayotte",
  "Mexico",
  "Micronesia",
  "Midway Islands",
  "Moldova",
  "Monaco",
  "Mongolia",
  "Montserrat",
  "Morocco",
  "Mozambique",
  "Namibia",
  "Nauru",
  "Navassa Island",
  "Nepal",
  "Netherlands",
  "Netherlands Antilles",
  "New Caledonia",
  "New Zealand",
  "Nicaragua",
  "Niger",
  "Nigeria",
  "Niue",
  "Norfolk Island",
  "North Korea",
  "Northern Mariana Islands",
  "Norway",
  "Oman",
  "Pakistan",
  "Palau",
  "Palmyra Atoll",
  "Panama",
  "Papua New Guinea",
  "Paracel Islands",
  "Paraguay",
  "Peru",
  "Philippines",
  "Pitcairn Islands",
  "Poland",
  "Portugal",
  "Puerto Rico",
  "Qatar",
  "Republic of the Congo",
  "Reunion",
  "Romania",
  "Russia",
  "Rwanda",
  "Saint Helena",
  "Saint Kitts and Nevis",
  "Saint Lucia",
  "Saint Pierre and Miquelon",
  "Saint Vincent and the Grenadines",
  "Samoa",
  "San Marino",
  "Sao Tome and Principe",
  "Saudi Arabia",
  "Senegal",
  "Seychelles",
  "Sierra Leone",
  "Singapore",
  "Slovakia",
  "Slovenia",
  "Solomon Islands",
  "Somalia",
  "South Africa",
  "South Georgia and the South Sandwich Islands",
  "South Korea",
  "Spain",
  "Spratly Islands",
  "Sri Lanka",
  "Sudan",
  "Suriname",
  "Svalbard",
  "Swaziland",
  "Sweden",
  "Switzerland",
  "Syria",
  "Taiwan",
  "Tajikistan",
  "Tanzania",
  "Thailand",
  "Togo",
  "Tokelau",
  "Tonga",
  "Trinidad and Tobago",
  "Tromelin Island",
  "Tunisia",
  "Turkey",
  "Turkmenistan",
  "Turks and Caicos Islands",
  "Tuvalu",
  "Uganda",
  "Ukraine",
  "United Arab Emirates",
  "United Kingdom",
  "Uruguay",
  "USA",
  "Uzbekistan",
  "Vanuatu",
  "Venezuela",
  "Viet Nam",
  "Virgin Islands",
  "Wake Island",
  "Wallis and Futuna",
  "West Bank",
  "Western Sahara",
  "Yemen",
  "Yugoslavia",
  "Zambia",
  "Zimbabwe"
};

static Boolean CountryIsValid (CharPtr name)

{
  Int2     L, R, mid;
  CharPtr  ptr;
  Char     str [256];

  if (StringHasNoText (name)) return FALSE;
  StringNCpy_0 (str, name, sizeof (str));
  ptr = StringChr (str, ':');
  if (ptr != NULL) {
    *ptr = '\0';
  }

  L = 0;
  R = sizeof (countrycodes) / sizeof (countrycodes [0]);

  while (L < R) {
    mid = (L + R) / 2;
    if (StringICmp (countrycodes [mid], str) < 0) {
      L = mid + 1;
    } else {
      R = mid;
    }
  }

  if (StringICmp (countrycodes [R], str) == 0) {
    return TRUE;
  }

  return FALSE;
}

/*****************************************************************************
*
*   ValidateSeqDescrContext(gcp)
*      Gather callback helper function for validating context on a Bioseq
*
*****************************************************************************/
static void ValidateBioSource (ValidStructPtr vsp, GatherContextPtr gcp, BioSourcePtr biop)

{
	CharPtr countryname;
	ValNodePtr db;
	DbtagPtr dbt;
	OrgRefPtr orp;
	SubSourcePtr ssp;

	if (biop == NULL) return;
	ssp = biop->subtype;
	while (ssp != NULL) {
		if (ssp->subtype == SUBSRC_country) {
			if (! CountryIsValid (ssp->name)) {
				countryname = ssp->name;
				if (StringHasNoText (countryname)) {
					countryname = "?";
				}
				ValidErr(vsp, SEV_WARNING, ERR_SEQ_DESCR_BadCountryCode, "Bad country name [%s]",
					countryname);
			}
		}
		ssp = ssp->next;
	}
	orp = biop->org;
	if (orp == NULL) return;
	if (GetAppProperty ("InternalNcbiSequin") == NULL) return;
	for (db = orp->db; db != NULL; db = db->next) {
		dbt = (DbtagPtr) db->data.ptrvalue;
		if (dbt != NULL) {
			if (StringICmp (dbt->db, "taxon") == 0) return;
		}
	}
	ValidErr(vsp, SEV_WARNING, ERR_SEQ_DESCR_NoTaxonID, "BioSource is missing taxon ID");
}

static Boolean ValidateSeqDescrCommon (ValNodePtr sdp, BioseqValidStrPtr bvsp, ValidStructPtr vsp, Uint2 descitemid)
{
	ValNodePtr vnp, vnp2;
	OrgRefPtr this_org=NULL, that_org=NULL;
	int tmpval;
	Char buf1[20], buf2[20];
	PubdescPtr pdp;
	MolInfoPtr mip;
	Uint2 olditemtype;
	Uint2 olditemid;
	BioSourcePtr biop;
	GatherContextPtr gcp = NULL;
	static char * badmod = "Inconsistent GIBB-mod [%d] and [%d]";

	vsp->sfp = NULL;
	vnp = sdp;
	vsp->descr = vnp;

	if (descitemid > 0) {
		gcp = vsp->gcp;
		if (gcp != NULL) {
			olditemid = gcp->itemID;
			olditemtype = gcp->thistype;
			gcp->itemID = descitemid;
			gcp->thistype = OBJ_SEQDESC;
		}
	}

	switch (vnp->choice)
	{
		case Seq_descr_mol_type:
			tmpval = (int)(vnp->data.intvalue);
			switch (tmpval)
			{
				case 8:  /* peptide */
					if (! bvsp->is_aa)
						ValidErr(vsp, SEV_ERROR, ERR_SEQ_DESCR_InvalidForType, "Nuclic acid with GIBB-mol = peptide");
					break;
				case 0:     /* unknown */
				case 255:   /* other */
					ValidErr(vsp,SEV_ERROR, ERR_SEQ_DESCR_InvalidForType , "GIBB-mol unknown or other used");
					break;
				default:    /* the rest are nucleic acid */
					if (bvsp->is_aa)
					{
						ValidErr(vsp, SEV_ERROR, ERR_SEQ_DESCR_InvalidForType, "GIBB-mol [%d] used on protein",
							tmpval);
					}
					else
					{
						if (bvsp->last_na_mol)
						{
							if (bvsp->last_na_mol != (int)vnp->data.intvalue)
							{
								ValidErr(vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "Inconsistent GIBB-mol [%d] and [%d]",
									bvsp->last_na_mol, tmpval);
							}
						}
						else
							bvsp->last_na_mol = tmpval;
					}
					break;
			}
			break;
		case Seq_descr_modif:
			for (vnp2 = (ValNodePtr)(vnp->data.ptrvalue); vnp2 != NULL; vnp2 = vnp2->next)
			{
				tmpval = (int)(vnp2->data.intvalue);
				switch (tmpval)
				{
					case 0:   /* dna */
					case 1:	  /* rna */
						if (bvsp->is_aa)	   /* only temporarily on 0 */
						{
							ValidErr(vsp,SEV_ERROR,ERR_SEQ_DESCR_InvalidForType , "Nucleic acid GIBB-mod [%d] on protein",
								tmpval);
						}
						else if (bvsp->last_na_mod)
						{
							if (tmpval != bvsp->last_na_mod)
							{
								ValidErr(vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, badmod,
									bvsp->last_na_mod, tmpval);
							}
						}
						else
							bvsp->last_na_mod = tmpval;
						break;
					case 4:   /* mitochondria */
					case 5:   /* chloroplast */
					case 6:   /* kinetoplast */
					case 7:   /* cyanelle */
					case 18:  /* macronuclear */
						if (bvsp->last_organelle)
						{
							if (tmpval != bvsp->last_na_mod)
							{
								ValidErr(vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, badmod,
									bvsp->last_organelle, tmpval);
							}
						}
						else
							bvsp->last_organelle = tmpval;
						break;
					case 10:  /* partial */
					case 11:  /* complete */
						if (bvsp->last_partialness)
						{
							if (tmpval != bvsp->last_partialness)
							{
								ValidErr(vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, badmod,
									bvsp->last_partialness, tmpval);
							}
						}
						else
							bvsp->last_partialness = tmpval;
						if ((bvsp->last_left_right) && (tmpval == 11))
						{
							ValidErr(vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, badmod,
								bvsp->last_left_right, tmpval);
						}
						break;
					case 16:   /* no left */
					case 17:   /* no right */
						if (bvsp->last_partialness == 11)  /* complete */
						{
							ValidErr(vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, badmod,
								bvsp->last_partialness, tmpval);
						}
						bvsp->last_left_right = tmpval;
						break;
					case 255:  /* other */
						ValidErr(vsp, SEV_ERROR, ERR_SEQ_DESCR_Unknown, "GIBB-mod = other used");
						break;
					default:
						break;
	
				}
			}
			break;
		case Seq_descr_method:
			if (! bvsp->is_aa)
			{
				ValidErr(vsp, SEV_ERROR,ERR_SEQ_DESCR_InvalidForType, "Nucleic acid with protein sequence method");
			}
			break;
		case Seq_descr_genbank:
			if (bvsp->last_gb != NULL)
				ValidErr(vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "Multiple GenBank blocks");
			else
				bvsp->last_gb = vnp;
			break;
		case Seq_descr_embl:
			if (bvsp->last_embl != NULL)
				ValidErr(vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "Multiple EMBL blocks");
			else
				bvsp->last_embl = vnp;
			break;
		case Seq_descr_pir:
			if (bvsp->last_pir != NULL)
				ValidErr(vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "Multiple PIR blocks");
			else
				bvsp->last_pir = vnp;
			break;
		case Seq_descr_sp:
			if (bvsp->last_sp != NULL)
				ValidErr(vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "Multiple SWISS-PROT blocks");
			else
				bvsp->last_sp = vnp;
			break;
		case Seq_descr_pdb:
			if (bvsp->last_pdb != NULL)
				ValidErr(vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "Multiple PDB blocks");
			else
				bvsp->last_pdb = vnp;
			break;
		case Seq_descr_prf:
			if (bvsp->last_prf != NULL)
				ValidErr(vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "Multiple PRF blocks");
			else
				bvsp->last_prf = vnp;
			break;
		case Seq_descr_create_date:
			if (bvsp->last_create != NULL)
			{
				tmpval = (int)DateMatch((DatePtr)vnp->data.ptrvalue,
					(DatePtr)(bvsp->last_create->data.ptrvalue), FALSE);
				if (tmpval)
				{
					DatePrint((DatePtr)(vnp->data.ptrvalue), buf1);
					DatePrint((DatePtr)(bvsp->last_create->data.ptrvalue), buf2);
					ValidErr(vsp, SEV_WARNING, ERR_SEQ_DESCR_Inconsistent, "Inconsistent create_dates [%s] and [%s]",
						buf1, buf2);
				}
			}
			else
				bvsp->last_create = vnp;
			if (bvsp->last_update != NULL)
			{
				tmpval = (int)DateMatch((DatePtr)vnp->data.ptrvalue,
					(DatePtr)(bvsp->last_update->data.ptrvalue), FALSE);
				if (tmpval == 1)
				{
					DatePrint((DatePtr)(vnp->data.ptrvalue), buf1);
					DatePrint((DatePtr)(bvsp->last_update->data.ptrvalue), buf2);
					ValidErr(vsp, SEV_WARNING, ERR_SEQ_DESCR_Inconsistent, "Inconsistent create_date [%s] and update_date [%s]",
						buf1, buf2);
				}
			}
			break;
		case Seq_descr_update_date:
			if (bvsp->last_create != NULL)
			{
				tmpval = (int)DateMatch((DatePtr)bvsp->last_create->data.ptrvalue,
					(DatePtr)(vnp->data.ptrvalue), FALSE);
				if (tmpval == 1)
				{
					DatePrint((DatePtr)(bvsp->last_create->data.ptrvalue), buf1);
					DatePrint((DatePtr)(vnp->data.ptrvalue), buf2);
					ValidErr(vsp, SEV_WARNING, ERR_SEQ_DESCR_Inconsistent, "Inconsistent create_date [%s] and update_date [%s]",
						buf1, buf2);
				}
			}
			if (bvsp->last_update == NULL)
				bvsp->last_update = vnp;
			break;
		case Seq_descr_source:
			biop = (BioSourcePtr) vnp->data.ptrvalue;
			/* ValidateBioSource (vsp, gcp, biop); */
			this_org = biop->org;
			/* fall into Seq_descr_org */
		case Seq_descr_org:
			if (this_org == NULL)
				this_org = (OrgRefPtr)(vnp->data.ptrvalue);
			if (bvsp->last_org != NULL)
			{
				if ((this_org->taxname != NULL) && (bvsp->last_org->taxname != NULL))
				{
					if (StringCmp(this_org->taxname, bvsp->last_org->taxname))
					{
						ValidErr(vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "Inconsistent taxnames [%s] and [%s]",
							this_org->taxname, bvsp->last_org->taxname);
					}
				}
			}
			else
				bvsp->last_org = this_org;
			for (vnp2 = vnp->next; vnp2 != NULL; vnp2 = vnp2->next) {
				if (vnp2->choice == Seq_descr_source || vnp2->choice == Seq_descr_org) {
					that_org = NULL;
					if (vnp2->choice == Seq_descr_source) {
						that_org = ((BioSourcePtr)(vnp2->data.ptrvalue))->org;
					}
					if (that_org == NULL) {
						that_org = (OrgRefPtr)(vnp2->data.ptrvalue);
					}
					if (that_org != NULL) {
						if ((this_org->taxname != NULL) && (that_org->taxname != NULL) &&
							StringCmp (this_org->taxname, that_org->taxname) == 0) {
							ValidErr(vsp, SEV_ERROR, ERR_SEQ_DESCR_MultipleBioSources,
								"Undesired multiple source descriptors");
						}
					}
				}
			}
			break;
		case Seq_descr_pub:
			bvsp->got_a_pub = TRUE;
			pdp = (PubdescPtr) vnp->data.ptrvalue;
			ValidatePubdesc (vsp, pdp);
			break;
		case Seq_descr_molinfo:
			mip = (MolInfoPtr) vnp->data.ptrvalue;
			if (mip != NULL) {
				switch (mip->biomol) {
					case MOLECULE_TYPE_PEPTIDE:  /* peptide */
						if (! bvsp->is_aa) {
							ValidErr(vsp, SEV_ERROR, ERR_SEQ_DESCR_InvalidForType, "Nuclic acid with Molinfo-biomol = peptide");
						}
						break;
					case 0:     /* unknown */
					case 255:   /* other */
						ValidErr(vsp,SEV_ERROR, ERR_SEQ_DESCR_InvalidForType , "Molinfo-biomol unknown or other used");
						break;
					default:   /* the rest are nucleic acid */
						if (bvsp->is_aa) {
							ValidErr(vsp, SEV_ERROR, ERR_SEQ_DESCR_InvalidForType, "Molinfo-biomol [%d] used on protein", (int) mip->biomol);
						} else {
							if (bvsp->last_biomol) {
								if (bvsp->last_biomol != (int) mip->biomol) {
									ValidErr(vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "Inconsistent Molinfo-biomol [%d] and [%d]",
										bvsp->last_biomol, (int) mip->biomol);
								}
							} else {
								bvsp->last_biomol = (int) mip->biomol;
							}
						}
						break;
				}
				if (! bvsp->is_aa) {
					switch (mip->tech) {
						case MI_TECH_concept_trans:
						case MI_TECH_seq_pept:
						case MI_TECH_both:
						case MI_TECH_seq_pept_overlap:
						case MI_TECH_seq_pept_homol:
						case MI_TECH_concept_trans_a:
							ValidErr(vsp, SEV_ERROR,ERR_SEQ_DESCR_InvalidForType, "Nucleic acid with protein sequence method");
							break;
						default:
							break;
					}
				}
				if (bvsp->last_tech) {
					if (bvsp->last_tech != (int) mip->tech) {
						ValidErr(vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "Inconsistent Molinfo-tech [%d] and [%d]",
							bvsp->last_tech, (int) mip->tech);
								}
				} else {
					bvsp->last_tech = (int) mip->tech;
				}
				if (bvsp->last_completeness) {
					if (bvsp->last_completeness != (int) mip->completeness) {
						ValidErr(vsp, SEV_ERROR, ERR_SEQ_DESCR_Inconsistent, "Inconsistent Molinfo-completeness [%d] and [%d]",
							bvsp->last_completeness, (int) mip->completeness);
								}
				} else {
					bvsp->last_completeness = (int) mip->completeness;
				}
			}
			break;
		default:
			break;
	}


	if (gcp != NULL) {
		gcp->itemID = olditemid;
		gcp->thistype = olditemtype;
	}

	return TRUE;
}

static Boolean LIBCALLBACK ValidateSeqDescrIndexed (ValNodePtr sdp, SeqMgrDescContextPtr context)

{
	ValidStructPtr vsp;
	BioseqValidStrPtr bvsp;

	bvsp = (BioseqValidStrPtr) context->userdata;
	vsp = bvsp->vsp;

	return ValidateSeqDescrCommon (sdp, bvsp, vsp, context->itemID);
}

static void ValidateSeqDescrContext (GatherContextPtr gcp)

{
	ValidStructPtr vsp;
	BioseqValidStrPtr bvsp;
	ValNodePtr sdp;

	bvsp = (BioseqValidStrPtr)(gcp->userdata);
	vsp = bvsp->vsp;
	sdp = (ValNodePtr)(gcp->thisitem);

	ValidateSeqDescrCommon (sdp, bvsp, vsp, 0);
}

/*****************************************************************************
*
*   ValidateBioseqContextGather(gcp)
*      Gather callback for validating context on a Bioseq
*
*****************************************************************************/
static void ValidateCitSub (ValidStructPtr vsp, CitSubPtr csp)

{
	AuthListPtr alp;
	ValNodePtr name;
	Boolean hasName = FALSE;

	if (vsp == NULL || csp == NULL) return;
	alp = csp->authors;
	if (alp != NULL) {
		if (alp->choice == 1) {
			for (name = alp->names; name != NULL; name = name->next) {
				if (! HasNoName (name)) {
					hasName = TRUE;
				}
			}
		} else if (alp->choice == 2 || alp->choice == 3) {
			for (name = alp->names; name != NULL; name = name->next) {
				if (! HasNoText ((CharPtr) name->data.ptrvalue)) {
					hasName = TRUE;
				}
			}
		}
	}
	if (! hasName) {
		ValidErr (vsp, SEV_ERROR, ERR_GENERIC_MissingPubInfo,
			"Submission citation has no author names");
	}
}

static Boolean ValidateBioseqContextIndexed (BioseqPtr bsp, BioseqValidStrPtr bvsp)
{
	ValidStructPtr vsp;
	CitSubPtr csp;
	ObjMgrDataPtr omdp;
	SeqSubmitPtr ssp;
	SubmitBlockPtr sbp;
	GatherContextPtr gcp;
	SeqFeatPtr sfp;
	SeqMgrFeatContext fcontext;
	Uint2 featdeftype;
	SeqFeatPtr last = NULL;
	Boolean leave;
	CharPtr label;
	CharPtr comment;
	Int4 left;
	Int4 right;
	Uint1 strand;
	Int2 numivals;
	Int4Ptr ivals;
	Boolean ivalssame;
	SeqAnnotPtr sap;
	Uint2 olditemtype;
	Uint2 olditemid;
	Int2 i;
	Int2 j;
	int severity;

	gcp = bvsp->gcp;
	vsp = bvsp->vsp;
	vsp->descr = NULL;
	vsp->sfp = NULL;
	vsp->gcp = gcp;   /* needed for ValidErr */

	SeqMgrExploreFeatures (bsp, (Pointer) bvsp, ValidateSeqFeatIndexed, NULL, NULL, NULL);

	if (gcp != NULL) {
		olditemid = gcp->itemID;
		olditemtype = gcp->thistype;
	}
	sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
	while (sfp != NULL) {
		leave = TRUE;
		if (last != NULL) {
			if (fcontext.left == left && fcontext.right == right &&
				fcontext.featdeftype == featdeftype) {
				if (fcontext.strand == strand ||
					strand == Seq_strand_unknown ||
					fcontext.strand == Seq_strand_unknown) {
					ivalssame = TRUE;
					if (fcontext.numivals != numivals || fcontext.ivals == NULL || ivals == NULL) {
						ivalssame = FALSE;
					} else {
						for (i = 0, j = 0; i < numivals; i++, j += 2) {
							if (fcontext.ivals [j] != ivals [j]) {
								ivalssame = FALSE;
							}
							if (fcontext.ivals [j + 1] != ivals [j + 1]) {
								ivalssame = FALSE;
							}
						}
					}
					if (ivalssame && /* StringICmp (fcontext.label, label) == 0 && */ fcontext.sap == sap) {
						if (gcp != NULL) {
							gcp->itemID = fcontext.itemID;
							gcp->thistype = OBJ_SEQFEAT;
						}
						vsp->descr = NULL;
						vsp->sfp = sfp;
						severity = SEV_ERROR;
						if (StringICmp (fcontext.label, label) != 0 ||
							StringICmp (sfp->comment, comment) != 0 ||
							featdeftype == FEATDEF_PUB) {
							severity = SEV_WARNING;
						}
						ValidErr (vsp, severity,ERR_SEQ_FEAT_DuplicateFeat, "Possible duplicate feature");
						vsp->sfp = NULL;
						if (gcp != NULL) {
							gcp->itemID = olditemid;
							gcp->thistype = olditemtype;
						}
					}
				}
			}
		}
		if (leave) {
			last = sfp;
			left = fcontext.left;
			right = fcontext.right;
			label = fcontext.label;
			comment = sfp->comment;
			strand = fcontext.strand;
			featdeftype = fcontext.featdeftype;
			numivals = fcontext.numivals;
			ivals = fcontext.ivals;
			sap = fcontext.sap;
		}
		sfp = SeqMgrGetNextFeature (bsp, sfp, 0, 0, &fcontext);
	}

	SeqMgrExploreDescriptors (bsp, (Pointer) bvsp, ValidateSeqDescrIndexed, NULL);

	omdp = ObjMgrGetData (gcp->entityID);
	if (omdp != NULL && omdp->datatype == OBJ_SEQSUB) {
		ssp = (SeqSubmitPtr) omdp->dataptr;
		if (ssp != NULL) {
			sbp = ssp->sub;
			if (sbp != NULL) {
				bvsp->got_a_pub = TRUE;
				csp = sbp->cit;
				/* csp = (CitSubPtr) gcp->thisitem; */
				ValidateCitSub (vsp, csp);
			}
		}
	}

	return TRUE;
}

static Boolean ValidateBioseqContextGather (GatherContextPtr gcp)
{
	ValidStructPtr vsp;
	BioseqValidStrPtr bvsp;
	CitSubPtr csp;

	bvsp = (BioseqValidStrPtr)(gcp->userdata);
	vsp = bvsp->vsp;
	vsp->descr = NULL;
	vsp->sfp = NULL;
	vsp->gcp = gcp;   /* needed for ValidErr */

	switch (gcp->thistype)
	{
		case OBJ_SEQFEAT:
			ValidateSeqFeatContext (gcp);
			break;
		case OBJ_SEQDESC:
			ValidateSeqDescrContext (gcp);
			break;
		case OBJ_SEQSUB_CIT:
			bvsp->got_a_pub = TRUE;
			csp = (CitSubPtr) gcp->thisitem;
			ValidateCitSub (vsp, csp);
			break;
		default:
			break;
	}
	return TRUE;
}

/*****************************************************************************
*
*   ValidateBioseqContext(gcp)
*      Validate one Bioseq for descriptors, features, and context
*      This is done as a second Gather, focussed on the Bioseq in
*        question.
*
*****************************************************************************/
static void ValidateBioseqContext (GatherContextPtr gcp)
{
	ValidStructPtr vsp;
	BioseqPtr bsp;
	GatherScope gs;
	BioseqValidStr bvs;
	SeqFeatPtr sfp;
	ValNode fake_whole;
	SeqIdPtr sip;
	ValNodePtr vnp = NULL;
	MolInfoPtr mip = NULL;
	SeqMgrDescContext context;
	BioseqContextPtr bcp;
	ObjMgrDataPtr omdp;

	vsp = (ValidStructPtr)(gcp->userdata);
	bsp = (BioseqPtr)(gcp->thisitem);
	vsp->bsp = bsp;
	vsp->descr = NULL;
	vsp->sfp = NULL;
	vsp->bssp = (BioseqSetPtr)(gcp->parentitem);

	MemSet(&gs, 0, sizeof(GatherScope));
	fake_whole.choice = SEQLOC_WHOLE;
	sip = SeqIdFindBest(bsp->id, 0);

	fake_whole.data.ptrvalue = sip;

	fake_whole.next = NULL;
	gs.target = &fake_whole;
	gs.get_feats_location = TRUE;
	gs.nointervals = TRUE;
	MemSet((Pointer)(gs.ignore), (int)TRUE, (size_t)(sizeof(Boolean) * OBJ_MAX));
	gs.ignore[OBJ_SEQDESC] = FALSE;
	gs.ignore[OBJ_SEQFEAT] = FALSE;
	gs.ignore[OBJ_SEQANNOT] = FALSE;
	gs.ignore[OBJ_SUBMIT_BLOCK] = FALSE;
	gs.ignore[OBJ_SEQSUB_CIT] = FALSE;

	gs.scope = vsp->sep;

	MemSet(&bvs, 0, sizeof(BioseqValidStr));
	bvs.vsp = vsp;

	/* now looking for molinfo on every bioseq (okay on segset) */
	if (bsp != NULL) {
		vnp = NULL;
		if (vsp->useSeqMgrIndexes) {
		  vnp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &context);
		} else {
		  bcp = BioseqContextNew(bsp);
		  vnp = BioseqContextGetSeqDescr (bcp, Seq_descr_molinfo, NULL, NULL);
		  BioseqContextFree (bcp);
		}
		if (vnp != NULL) {
			mip = (MolInfoPtr) vnp->data.ptrvalue;
		}
	}

	bvs.is_mrna = FALSE;
	bvs.is_prerna = FALSE;
	if (bsp != NULL && ISA_na (bsp->mol)) {
		if (mip != NULL) {
			if (mip->biomol == MOLECULE_TYPE_MRNA) {
				bvs.is_mrna = TRUE;
			} else if (mip->biomol == MOLECULE_TYPE_PRE_MRNA) {
			  bvs.is_prerna = TRUE;
			}
		} else if (bsp->mol == Seq_mol_rna) {
			bvs.is_mrna = TRUE;  /* if no molinfo, assume rna is mrna */
		}
	}

	if (ISA_aa(bsp->mol))
	{
		bvs.is_aa = TRUE;
		    /* check proteins in nuc-prot set have a CdRegion */
		if (vsp->bssp != NULL)
		{
			if (vsp->bssp->_class == 1)   /* in a nuc-prot set */
			{
				if (vsp->useSeqMgrIndexes) {
					sfp = SeqMgrGetCDSgivenProduct (bsp, NULL);
				} else {
					sfp = SeqEntryGetSeqFeat(vsp->sep, 3, NULL, NULL, 1, bsp);
				}
				if (sfp == NULL)   /* no CdRegion points to this bsp */
					ValidErr(vsp, SEV_ERROR, ERR_SEQ_PKG_NoCdRegionPtr,	"No CdRegion in nuc-prot set points to this protein");
			}
		}
	}

	if (vsp->useSeqMgrIndexes) {
		bvs.gcp = gcp;
		ValidateBioseqContextIndexed (bsp, &bvs);
	} else {
		GatherSeqEntry (vsp->sep, &bvs, ValidateBioseqContextGather, &gs);
	}

	vsp->gcp = gcp;    /* reset the gcp pointer changed in previous gather */
	vsp->descr = NULL;
	vsp->sfp = NULL;

	if ((! bvs.got_a_pub) && (! vsp->suppress_no_pubs)) {
		omdp = NULL;
		if (gcp != NULL) {
			omdp = ObjMgrGetData (gcp->entityID);
		}
		if (omdp == NULL || omdp->datatype != OBJ_SEQSUB) {
			ValidErr(vsp, SEV_ERROR, ERR_SEQ_DESCR_NoPubFound, "No publications refer to this Bioseq.");
		}
	}

	if ((! bvs.last_org) && (! vsp->suppress_no_biosrc))
		ValidErr(vsp, SEV_ERROR, ERR_SEQ_DESCR_NoOrgFound, "No organism name has been applied to this Bioseq.");


	if ((bvs.is_aa) && (! bvs.num_full_length_prot_ref))
		ValidErr(vsp, SEV_ERROR, ERR_SEQ_FEAT_NoProtRefFound, "No full length Prot-ref feature applied to this Bioseq");

	/* for now only flag missing molinfo in Sequin */
	if (mip == NULL && GetAppProperty ("SpliceValidateAsError") != NULL) {
		ValidErr(vsp, SEV_ERROR, ERR_SEQ_DESCR_NoMolInfoFound, "No Mol-info applies to this Bioseq");
	}

	return;

}

/*****************************************************************************
*
*   ValidateSeqFeat(gcp)
*
*****************************************************************************/
static Boolean EmptyOrNullString (CharPtr str)

{
  Char  ch;

  if (str == NULL) return TRUE;
  ch = *str;
  while (ch != '\0') {
    if (ch > ' ' && ch <= '~') return FALSE;
    str++;
    ch = *str;
  }
  return TRUE;
}

static void ValidateImpFeat (ValidStructPtr vsp, GatherContextPtr gcp, SeqFeatPtr sfp, ImpFeatPtr ifp)

{
  Boolean    found;
  GBQualPtr  gbqual;
  Int2       i;
  Int2       index;
  CharPtr    key;
  Int2       qual;
  Int2       val;

  if (vsp == NULL || gcp == NULL || sfp == NULL || ifp == NULL) return;
  if (StringCmp (ifp->key, "-") == 0) {
    key = StringSave ("misc_feature");
  } else {
    key = StringSaveNoNull (ifp->key);
  }
  index = GBFeatKeyNameValid (&key, FALSE);
  if (index == -1) {
    if (key != NULL) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UnknownImpFeatKey, "Unknown feature key %s", key);
    } else {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UnknownImpFeatKey, "NULL feature key");
    }
  }
  if (StringICmp (key, "mat_peptide") == 0 ||
      StringICmp (key, "sig_peptide") == 0 ||
      StringICmp (key, "transit_peptide") == 0) {
    ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_InvalidForType,
              "Peptide processing feature should be converted to the appropriate protein feature subtype");
  }
  for (gbqual = sfp->qual; gbqual != NULL; gbqual = gbqual->next) {
    if (StringCmp (gbqual->qual, "gsdb_id") == 0) {
      continue;
    }
    val = GBQualNameValid (gbqual->qual);
    if (val == -1) {
      if (gbqual->qual != NULL) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UnknownImpFeatQual, "Unknown qualifier %s", gbqual->qual);
      } else {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_UnknownImpFeatQual, "NULL qualifier");
      }
    } else if (index != -1) {
      found = FALSE;
      for (i = 0; i < ParFlat_GBFeat [index].opt_num; i++) {
        qual = ParFlat_GBFeat [index].opt_qual [i];
        if (qual == val) {
          found = TRUE;
          break;
        }
      }
      if (! found) {
        for (i = 0; i < ParFlat_GBFeat [index].mand_num; i++) {
          qual = ParFlat_GBFeat [index].mand_qual [i];
          if (qual == val) {
            found = TRUE;
            break;
          }
        }
        if (! found) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_WrongQualOnImpFeat, "Wrong qualifier %s for feature %s", gbqual->qual, key);
        }
      }
    }
  }
  if (index != -1 && ParFlat_GBFeat [index].mand_num > 0) {
    for (i = 0; i < ParFlat_GBFeat [index].mand_num; i++) {
      found = FALSE;
      qual = ParFlat_GBFeat [index].mand_qual [i];
      for (gbqual = sfp->qual; gbqual != NULL; gbqual = gbqual->next) {
        val = GBQualNameValid (gbqual->qual);
        if (qual == val) {
          found = TRUE;
          break;
        }
      }
      if (! found) {
        if (qual == GBQUAL_citation && sfp->cit != NULL) {
          found = TRUE;
        }
      }
      if (! found) {
        ValidErr (vsp, SEV_WARNING, ERR_SEQ_FEAT_MissingQualOnImpFeat, "Missing qualifier %s for feature %s", ParFlat_GBQual_names [qual].name, key);
      }
    }
  }
  MemFree (key);
}

/* PartialAtSpliceSite uses code taken from SpliceCheckEx */
static Boolean PartialAtSpliceSite (SeqLocPtr head, Uint2 slpTag)

{
  BioseqPtr   bsp;
  Int2        residue1, residue2;
  Boolean     rsult = FALSE;
  SeqIdPtr    sip;
  SeqLocPtr   slp = NULL, first = NULL, last = NULL;
  SeqPortPtr  spp = NULL;
  Uint1       strand;
  Int4        strt, stp, donor, acceptor, len;

  if (slpTag != SLP_NOSTART && slpTag != SLP_NOSTOP) return FALSE;
  while ((slp = SeqLocFindPart (head, slp, EQUIV_IS_ONE)) != NULL) {
    if (first == NULL) {
      first = slp;
    }
    last = slp;
  }
  if (first == NULL) return FALSE;

  strand = SeqLocStrand (first);
  if (SeqLocStrand (last) != strand) return FALSE;

  if (slpTag == SLP_NOSTART) {
    slp = first;
  } else {
    slp = last;
  }
  sip = SeqLocId (slp);
  if (sip == NULL) return FALSE;
  acceptor = SeqLocStart (slp);
  donor = SeqLocStop (slp);
  bsp = BioseqLockById (sip);
  if (bsp == NULL) return FALSE;
  len = bsp->length;
  spp = SeqPortNew (bsp, 0, -1, strand, Seq_code_ncbi4na);
  BioseqUnlock (bsp);
  if (spp == NULL) return FALSE;

  if (strand != Seq_strand_minus) {
    strt = acceptor;
    stp = donor;
  } else {
    strt = donor;
    donor = acceptor;
    acceptor = strt;
    stp = len - donor - 1;
    strt = len - acceptor - 1;
  }

  if (slpTag == SLP_NOSTOP && stp < len - 2) {
    SeqPortSeek (spp, (stp + 1), SEEK_SET);
    residue1 = SeqPortGetResidue (spp);
    residue2 = SeqPortGetResidue (spp);
    if (IS_residue (residue1) && IS_residue (residue2)) {
      if ((residue1 & 4) && (residue2 & 8)) {
        rsult = TRUE;
      } else if ((residue1 & 4) && (residue2 & 2)) {
        rsult = TRUE;
      }
    }
  } else if (slpTag == SLP_NOSTART && strt > 1) {
    SeqPortSeek (spp, (strt - 2), SEEK_SET);
    residue1 = SeqPortGetResidue (spp);
    residue2 = SeqPortGetResidue (spp);
    if (IS_residue (residue1) && IS_residue (residue2)) {
      if ((residue1 & 1) && (residue2 & 4)) {
        rsult = TRUE;
      }
    }
  }

  spp = SeqPortFree (spp);
  return rsult;
}

static void CheckTrnaCodons (ValidStructPtr vsp, GatherContextPtr gcp, SeqFeatPtr sfp, tRNAPtr trp)

{
  Uint1           aa;
  BioseqPtr       bsp;
  Int2            code;
  CharPtr         codes = NULL;
  Uint1           from;
  GeneticCodePtr  gncp;
  Int2            j;
  SeqEntryPtr     sep;
  SeqMapTablePtr  smtp;
  Uint1           taa;
  ValNodePtr      vnp;

  if (vsp == NULL || gcp == NULL || sfp == NULL || trp == NULL) return;
  for (j = 0; j < 6; j++) {
    if (trp->codon [j] < 64) {
      if (codes == NULL) {
        bsp = GetBioseqGivenSeqLoc (sfp->location, gcp->entityID);
        sep = GetBestTopParentForData (gcp->entityID, bsp);
        code = SeqEntryToGeneticCode (sep, NULL, NULL, 0);
        gncp = GeneticCodeFind (code, NULL);
        if (gncp == NULL) {
          gncp = GeneticCodeFind (1, NULL);
        }
        if (gncp == NULL) return;
        for (vnp = (ValNodePtr) gncp->data.ptrvalue; vnp != NULL; vnp = vnp->next) {
          if (vnp->choice == 3) {
            codes = (CharPtr) vnp->data.ptrvalue;
          }
        }
      }
      if (codes == NULL) return;
      taa = codes [trp->codon [j]];
      aa = 0;
      if (trp->aatype == 2) {
        aa = trp->aa;
      } else {
        from = 0;
        switch (trp->aatype) {
          case 0 :
            from = 0;
            break;
          case 1 :
            from = Seq_code_iupacaa;
            break;
          case 2 :
            from = Seq_code_ncbieaa;
            break;
          case 3 :
            from = Seq_code_ncbi8aa;
            break;
          case 4 :
            from = Seq_code_ncbistdaa;
            break;
          default:
            break;
        }
        smtp = SeqMapTableFind (Seq_code_ncbieaa, from);
        if (smtp != NULL) {
          aa = SeqMapTableConvert (smtp, trp->aa);
        }
      }
      if (aa > 0 && aa != 255) {
        if (taa != aa) {
          ValidErr (vsp, SEV_ERROR, ERR_SEQ_FEAT_TrnaCodonWrong, "tRNA codon does not match genetic code");
        }
      }
    }
  }
}

NLM_EXTERN void ValidateSeqFeat( GatherContextPtr gcp)
{
	Int2 type, i, j;
	static char * errclass = "BadLoc";
	static char * parterr[2] = {"PartialProduct", "PartialLocation"};
	static char * parterrs[4] = {
		"Start does not include first/last residue of sequence",
		"Stop does not include first/last residue of sequence",
		"Internal partial intervals do not include first/last residue of sequence",
		"Improper use of partial (greater than or less than)" };
	Uint2 partials[2], errtype;
	Char buf[80];
	CharPtr tmp;
	ValidStructPtr vsp;
	SeqFeatPtr sfp;
	CdRegionPtr crp;
	CodeBreakPtr cbp;
	RnaRefPtr rrp;
	tRNAPtr trp;
	GBQualPtr gbq;
	Boolean pseudo, excpt;
	ImpFeatPtr ifp;
	GeneRefPtr grp;
	ProtRefPtr prp;
	ValNodePtr vnp;
	BioseqPtr bsp;
	BioseqContextPtr bcp;
	BioSourcePtr biop;
	OrgNamePtr onp;
	OrgRefPtr orp;
	Int2 biopgencode;
	Int2 cdsgencode;
	GeneticCodePtr gc;
	PubdescPtr pdp;
	DbtagPtr db = NULL;
	Int4 id = -1;
	SeqMgrDescContext context;
	GeneRefPtr grpx;
	SeqFeatPtr sfpx;
	Boolean redundantgenexref;
	SeqMgrFeatContext fcontext;
	CharPtr syn1, syn2, label = NULL;

	vsp = (ValidStructPtr)(gcp->userdata);
	sfp = (SeqFeatPtr)(gcp->thisitem);
	vsp->descr = NULL;
	vsp->sfp = sfp;
	type = (Int2)(sfp->data.choice);

	ValidateSeqLoc(vsp, sfp->location, "Location");

	ValidateSeqLoc(vsp, sfp->product, "Product");

	partials[0] = SeqLocPartialCheck(sfp->product);
	partials[1] = SeqLocPartialCheck(sfp->location);
	if ((partials[0] != SLP_COMPLETE) || (partials[1] != SLP_COMPLETE) || (sfp->partial))  /* partialness */
	{
		                           /* a feature on a partial sequence should be partial -- if often isn't */
		if ((! sfp->partial) && (partials[1] != SLP_COMPLETE) &&
			(sfp->location->choice == SEQLOC_WHOLE))
		{
			ValidErr(vsp, SEV_INFO, ERR_SEQ_FEAT_PartialProblem, "On partial Bioseq, SeqFeat.partial should be TRUE");
		}
								  /* a partial feature, with complete location, but partial product */
		else if ((sfp->partial) && (sfp->product != NULL) &&
			(partials[1] == SLP_COMPLETE) && (sfp->product->choice == SEQLOC_WHOLE)
			&& (partials[0] != SLP_COMPLETE))
		{
			ValidErr(vsp, SEV_WARNING, ERR_SEQ_FEAT_PartialProblem,
				"When SeqFeat.product is a partial Bioseq, SeqFeat.location should also be partial");
		}
		                         /* gene on segmented set is now 'order', should also be partial */
		else if (type == SEQFEAT_GENE && sfp->product == NULL && partials [1] == SLP_INTERNAL) {
			if (! sfp->partial) {
				ValidErr(vsp, SEV_WARNING, ERR_SEQ_FEAT_PartialProblem,
					"Gene of 'order' with otherwise complete location should have partial flag set");
			}
		}
		                         /* inconsistent combination of partial/complete product,location,partial flag */
		else if (((partials[0] == SLP_COMPLETE) && (sfp->product != NULL)) ||
			(partials[1] == SLP_COMPLETE) ||
			(! sfp->partial))
		{
			tmp = StringMove(buf, "Inconsistent: ");
			if (sfp->product != NULL)
			{
				tmp = StringMove(tmp, "Product= ");
				if (partials[0])
					tmp = StringMove(tmp, "partial, ");
				else
					tmp = StringMove(tmp, "complete, ");
			}
			tmp = StringMove(tmp, "Location= ");
			if (partials[1])
				tmp = StringMove(tmp, "partial, ");
			else
				tmp = StringMove(tmp, "complete, ");
			tmp = StringMove(tmp, "Feature.partial= ");
			if (sfp->partial)
				tmp = StringMove(tmp, "TRUE");
			else
				tmp = StringMove(tmp, "FALSE");
			ValidErr(vsp, SEV_WARNING, ERR_SEQ_FEAT_PartialProblem, buf);
		}

		                      /* may have other error bits set as well */
		for (i = 0; i < 2; i++)
		{
			errtype = SLP_NOSTART;
			for (j = 0; j < 4; j++)
			{
				if (partials[i] & errtype)
				{
					if (i == 1 && j < 2 && PartialAtSpliceSite (sfp->location, errtype)) {
						ValidErr(vsp, SEV_INFO, ERR_SEQ_FEAT_PartialProblem,
						         "%s: %s (but is at consensus splice site)", parterr[i], parterrs[j]);
					} else {
					  ValidErr(vsp, SEV_WARNING, ERR_SEQ_FEAT_PartialProblem, "%s: %s", parterr[i], parterrs[j]);
					}
				}
				errtype <<= 1;
			}
		}

	}

	for (vnp = sfp->dbxref; vnp != NULL; vnp = vnp->next) {
		id = -1;
		db = vnp->data.ptrvalue;
		if (db && db->db) {
			for (i =0; i < DBNUM; i++) {
				if (StringCmp(db->db, dbtag[i]) == 0) {
					id = i;
					break;
				}
			}
			if (id == -1 || (type != SEQFEAT_CDREGION && id < 4)) {
				ValidErr(vsp, SEV_WARNING, ERR_SEQ_FEAT_IllegalDbXref,
					     "Illegal db_xref type %s", db->db);
			}
		}
	}

	switch (type)
	{
		case 1:        /* Gene-ref */
			grp = (GeneRefPtr) (sfp->data.value.ptrvalue);
			if (grp != NULL) {
			  if (EmptyOrNullString (grp->locus) &&
			      EmptyOrNullString (grp->allele) &&
			      EmptyOrNullString (grp->desc) &&
			      EmptyOrNullString (grp->maploc) &&
			      grp->db == NULL && grp->syn == NULL) {
			  	ValidErr(vsp, SEV_WARNING, ERR_SEQ_FEAT_GeneRefHasNoData,
				        "There is a gene feature where all fields are empty");
			  }
			for (vnp = grp->db; vnp != NULL; vnp = vnp->next) {
				id = -1;
				db = vnp->data.ptrvalue;
				if (db && db->db) {
					for (i =0; i < DBNUM; i++) {
						if (StringCmp(db->db, dbtag[i]) == 0) {
							id = i;
							break;
						}
					}
					if (id == -1 || (type != SEQFEAT_CDREGION && id < 4)) {
						ValidErr(vsp, SEV_WARNING, ERR_SEQ_FEAT_IllegalDbXref,
							     "Illegal db_xref type %s", db->db);
					}
				}
			}
			}
			break;
		case 2:        /* Org-ref */
			break;
		case 3:        /* Cdregion */
			pseudo = sfp->pseudo; /* now also uses new feature pseudo flag */
			excpt = FALSE;
			gbq = sfp->qual;
			while (gbq != NULL) {
				if (StringICmp (gbq->qual, "pseudo") == 0) {
					pseudo = TRUE;
				}
				if (StringICmp (gbq->qual, "exception") == 0) {
					excpt = TRUE;
				}
				gbq = gbq->next;
			}
			if (! pseudo) {
				CdTransCheck(vsp, sfp);
				SpliceCheck(vsp, sfp);
			}
			crp = (CdRegionPtr)(sfp->data.value.ptrvalue);
			if (crp != NULL) {
			for (cbp = crp->code_break; cbp != NULL; cbp = cbp->next)
			{
				i = SeqLocCompare(cbp->loc, sfp->location);
				if ((i != SLC_A_IN_B) && (i != SLC_A_EQ_B))
					ValidErr(vsp, SEV_ERROR, ERR_SEQ_FEAT_Range, "Code-break location not in coding region");
			}
			if (excpt && (! sfp->excpt)) {
				ValidErr(vsp, SEV_WARNING, ERR_SEQ_FEAT_ExceptInconsistent,
				        "Exception flag should be set in coding region");
			}
			if (crp->orf && sfp->product != NULL) {
				ValidErr(vsp, SEV_WARNING, ERR_SEQ_FEAT_OrfCdsHasProduct,
				        "An ORF coding region should not have a product");
			}
			if (pseudo && sfp->product != NULL) {
				ValidErr(vsp, SEV_WARNING, ERR_SEQ_FEAT_PsuedoCdsHasProduct,
				        "A pseudo coding region should not have a product");
			}
			biopgencode = 0;
			cdsgencode = 0;
			bsp = GetBioseqGivenSeqLoc (sfp->location, gcp->entityID);
			if (bsp != NULL) {
				vnp = NULL;
				if (vsp->useSeqMgrIndexes) {
					vnp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &context);
				} else {
					bcp = BioseqContextNew (bsp);
					vnp = BioseqContextGetSeqDescr (bcp, Seq_descr_source, NULL, NULL);
				}
					if (vnp != NULL && vnp->data.ptrvalue != NULL) {
						biop = (BioSourcePtr) vnp->data.ptrvalue;
						orp = biop->org;
						if (orp != NULL && orp->orgname != NULL) {
							onp = orp->orgname;
							if (biop->genome == 4 || biop->genome == 5) {
								biopgencode = onp->mgcode;
							} else {
								biopgencode = onp->gcode;
							}
							gc = crp->genetic_code;
							if (gc != NULL) {
								for (vnp = gc->data.ptrvalue; vnp != NULL; vnp = vnp->next) {
									if (vnp->choice == 2) {
										cdsgencode = (Int2) vnp->data.intvalue;
									}
								}
							}
							if (biopgencode != cdsgencode) {
								ValidErr(vsp, SEV_WARNING, ERR_SEQ_FEAT_GenCodeMismatch,
								        "Genetic code conflict between CDS (code %d) and BioSource (code %d)",
								        (int) cdsgencode, (int) biopgencode);
							}
						}
					}
					if (! vsp->useSeqMgrIndexes) {
						BioseqContextFree (bcp);
					}
			}
			}
			break;
		case 4:        /* Prot-ref */
			prp = (ProtRefPtr) (sfp->data.value.ptrvalue);
			if (prp != NULL) {
			  if (prp->processed != 3 && prp->processed != 4) {
			    vnp = prp->name;
			    if ((vnp == NULL || EmptyOrNullString ((CharPtr) vnp->data.ptrvalue)) &&
			        EmptyOrNullString (prp->desc) &&
					prp->ec == NULL && prp->activity == NULL && prp->db == NULL) {
					ValidErr(vsp, SEV_WARNING, ERR_SEQ_FEAT_ProtRefHasNoData,
					"There is a protein feature where all fields are empty");
			    }
			  }
			for (vnp = prp->db; vnp != NULL; vnp = vnp->next) {
				id = -1;
				db = vnp->data.ptrvalue;
				if (db && db->db) {
					for (i =0; i < DBNUM; i++) {
						if (StringCmp(db->db, dbtag[i]) == 0) {
							id = i;
							break;
						}
					}
					if (id == -1 || (type != SEQFEAT_CDREGION && id < 4)) {
						ValidErr(vsp, SEV_WARNING, ERR_SEQ_FEAT_IllegalDbXref,
							     "Illegal db_xref type %s", db->db);
					}
				}
			}
			}
			break;
		case 5:        /* RNA-ref */
			rrp = (RnaRefPtr)(sfp->data.value.ptrvalue);
			if (rrp->type == 2) { /* mRNA */
				SpliceCheck(vsp, sfp);
			}
			if (rrp->ext.choice == 2)   /* tRNA */
			{
				trp = (tRNAPtr)(rrp->ext.value.ptrvalue);
				if (trp->anticodon != NULL)
				{
					i = SeqLocCompare(trp->anticodon, sfp->location);
					if ((i != SLC_A_IN_B) && (i != SLC_A_EQ_B))
						ValidErr(vsp, SEV_ERROR, ERR_SEQ_FEAT_Range, "Anticodon location not in tRNA");
				}
				CheckTrnaCodons (vsp, gcp, sfp, trp);
			}
			if (rrp->type == 0) {
				ValidErr(vsp, SEV_WARNING, ERR_SEQ_FEAT_RNAtype0, "RNA type 0 (unknown) not supported");
			}
			break;
		case 6:        /* Pub */
			pdp = (PubdescPtr) sfp->data.value.ptrvalue;
			ValidatePubdesc (vsp, pdp);
			break;
		case 7:        /* Seq */
			break;
		case 8:        /* Imp-feat */
			ifp = (ImpFeatPtr) sfp->data.value.ptrvalue;
			if (GetAppProperty ("ValidateExons") != NULL) {
				
				if (ifp != NULL && StringICmp (ifp->key, "exon") == 0) {
					SpliceCheckEx (vsp, sfp, TRUE);
				}
			}
			if (ifp != NULL) {
				ValidateImpFeat (vsp, gcp, sfp, ifp);
			}
			break;
		case 9:        /* Region */
			break;
		case 10:        /* Comment */
			break;
		case 11:        /* Bond */
			break;
		case 12:        /* Site */
			break;
		case 13:        /* Rsite-ref */
			break;
		case 14:        /* User-object */
			break;
		case 15:        /* TxInit */
			break;
		case 16:        /* Numbering */
			break;
		case 17:        /* Secondary Structure */
			break;
		case 18:        /* NonStdRes*/
			break;
		case 19:        /* Heterogen*/
			break;
		case 20:        /* BioSource*/
			/*
			biop = (BioSourcePtr) sfp->data.value.ptrvalue;
			ValidateBioSource (vsp, gcp, biop);
			*/
			break;
		default:
			ValidErr(vsp, SEV_ERROR, ERR_SEQ_FEAT_InvalidType, "Invalid SeqFeat type [%d]",
				(int)(type));
			break;
	}
	if (type != SEQFEAT_GENE) {
	  grp = SeqMgrGetGeneXref (sfp);
	  if (grp == NULL || SeqMgrGeneIsSuppressed (grp)) return;
	  sfpx = SeqMgrGetOverlappingGene (sfp->location, &fcontext);
	  if (sfpx == NULL || sfpx->data.choice != SEQFEAT_GENE) return;
	  grpx = (GeneRefPtr) sfpx->data.value.ptrvalue;
	  if (grpx == NULL) return;
	  redundantgenexref = FALSE;
	  label = fcontext.label;
	  if ((! StringHasNoText (grp->locus)) && (! StringHasNoText (grpx->locus))) {
	    if ((StringICmp (grp->locus, grpx->locus) == 0)) {
	      redundantgenexref = TRUE;
	      label = grp->locus;
	    }
	  } else if (grp->syn != NULL && grpx->syn != NULL) {
	    syn1 = (CharPtr) grp->syn->data.ptrvalue;
	    syn2 = (CharPtr) grpx->syn->data.ptrvalue;
	    if ((! StringHasNoText (syn1)) && (! StringHasNoText (syn2))) {
	      if ((StringICmp (syn1, syn2) == 0)) {
	        redundantgenexref = TRUE;
	        label = syn1;
	      }
	    }
	  }
	  if (redundantgenexref) {
	    if (StringHasNoText (label)) {
	      label = "?";
	    }
		ValidErr(vsp, SEV_WARNING, ERR_SEQ_FEAT_UnnecessaryGeneXref, "Unnecessary gene cross-reference %s", label);
	  }
	}
	return;
}

/*****************************************************************************
*
*   CdTransCheck(sfp)
*   	Treatment of terminal 'X'
*          If either the protein or the translation end in 'X' (usually
*          due to partial last codon) it is ignored to minimize conflicts
*          between approaches to add the X or not in this case.
*
*****************************************************************************/
static CharPtr MapToNTCoords (SeqFeatPtr sfp, SeqIdPtr protID, Int4 pos)

{
  SeqLocPtr  nslp;
  SeqLocPtr  pslp;
  CharPtr    rsult;
  SeqPntPtr  spntp;

  rsult = NULL;
  if (sfp != NULL && protID != NULL && pos >= 0) {
    spntp = SeqPntNew ();
    pslp = ValNodeNew (NULL);
    pslp->choice = SEQLOC_PNT;
    pslp->data.ptrvalue = (Pointer) spntp;
    spntp->point = pos;
    spntp->id = SeqIdDup (protID);
    nslp = aaLoc_to_dnaLoc (sfp, pslp);
    if (nslp != NULL) {
      rsult = SeqLocPrint (nslp);
    }
    SeqLocFree (pslp);
    SeqLocFree (nslp);
  }
  return rsult;
}

NLM_EXTERN void CdTransCheck(ValidStructPtr vsp, SeqFeatPtr sfp)
{
	ByteStorePtr newprot = NULL;
	BioseqPtr prot1seq=NULL, prot2seq=NULL;
	SeqLocPtr slp=NULL, curr = NULL;
	Int4 prot1len = 0, prot2len, i, len;
	CdRegionPtr crp;
	SeqIdPtr protid=NULL;
	Int2 residue1, residue2, stop_count = 0, mismatch = 0, ragged = 0;
	Boolean got_stop = FALSE, test_it = TRUE;
	SeqPortPtr spp=NULL;
	Uint2 part_loc=0, part_prod=0;
	Boolean no_end = FALSE, no_beg = FALSE, show_stop = FALSE,
		got_dash = FALSE, done;
	GBQualPtr gb;
	ValNodePtr vnp, code;
	int gccode = 0;
	Boolean transl_except = FALSE, prot_ok = TRUE;
	CharPtr nuclocstr;
	CodeBreakPtr cbp;
	Int4 pos1, pos2, pos;
	SeqLocPtr tmp;

	if (sfp == NULL) return;

	if (sfp->excpt)		 /* biological exception */
		return;

	for (gb = sfp->qual; gb != NULL; gb = gb->next)	  /* pseuogene */
	{
		if (! StringICmp("pseudo", gb->qual))
			return;
	}

	crp = (CdRegionPtr)(sfp->data.value.ptrvalue);
	if (crp->code_break == NULL)  /* check for unparsed transl_except */
	{
		for (gb = sfp->qual; gb != NULL; gb = gb->next)
		{
			if (! StringCmp(gb->qual, "transl_except"))
			{
				transl_except = TRUE;
				break;
			}
		}
	}

	if (crp->genetic_code != NULL)
	{
		for (vnp = crp->genetic_code->data.ptrvalue; ((vnp != NULL) && (! gccode)); vnp = vnp->next)
		{
			switch (vnp->choice)
			{
				case 0:
				    break;
				case 1:     /* name */
					code = GeneticCodeFind(0, (CharPtr)(vnp->data.ptrvalue));
					if (code != NULL)
					{
						for (vnp = code->data.ptrvalue; ((vnp != NULL) && (! gccode)); vnp = vnp->next)
						{
							if (vnp->choice == 2) /* id */
								gccode = (int)(vnp->data.intvalue);
						}
					}
					break;
				case 2:    /* id */
				    gccode = (int)(vnp->data.intvalue);
					break;
				default:
				    gccode = 255;
					break;
			}
		}
	}
					
	newprot = ProteinFromCdRegion(sfp, TRUE);   /* include stop codons */
	if (newprot == NULL)
	{
		ValidErr(vsp, SEV_ERROR, ERR_SEQ_FEAT_CdTransFail, "Unable to translate");
		prot_ok = FALSE;
		goto erret;
	}

	part_loc = SeqLocPartialCheck(sfp->location);
	part_prod = SeqLocPartialCheck(sfp->product);
	if ((part_loc & SLP_STOP) || (part_prod & SLP_STOP))
		no_end = TRUE;
	else    /* complete stop, so check for ragged end */
	{
		len = SeqLocLen(sfp->location);
		if (crp->frame > 1)
			len -= (Int4)(crp->frame - 1);
		ragged = (Int2)(len % (Int4)(3));
		if (ragged) {
		len = SeqLocLen(sfp->location);
		cbp = crp->code_break;
		while (cbp != NULL)
		{
			pos1 = INT4_MAX;
			pos2 = -10;
			tmp = NULL;
			while ((tmp = SeqLocFindNext(cbp->loc, tmp)) != NULL)
			{
				pos = GetOffsetInLoc(tmp, sfp->location, 
SEQLOC_START);
				if (pos < pos1)
					pos1 = pos;
				pos = GetOffsetInLoc(tmp, sfp->location, 
SEQLOC_STOP);
				if (pos > pos2)
					pos2 = pos;
			}
			pos = pos2 - pos1; /* codon length */
			if (pos >= 0 && pos <= 1 && pos2 == len - 1)   /*  a codon */
			/* allowing a partial codon at the end */
			{
				ragged = 0;
			}

			cbp = cbp->next;
		}
		}
	}

		/* check for code break not on a codon */
		len = SeqLocLen(sfp->location);
		cbp = crp->code_break;
		while (cbp != NULL)
		{
			pos1 = INT4_MAX;
			pos2 = -10;
			tmp = NULL;
			while ((tmp = SeqLocFindNext(cbp->loc, tmp)) != NULL)
			{
				pos = GetOffsetInLoc(tmp, sfp->location, 
SEQLOC_START);
				if (pos < pos1)
					pos1 = pos;
				pos = GetOffsetInLoc(tmp, sfp->location, 
SEQLOC_STOP);
				if (pos > pos2)
					pos2 = pos;
			}
			pos = pos2 - pos1; /* codon length */
			/* check for code break not on a codon */
			if (pos == 2 || (pos >= 0 && pos <= 1 && pos2 == len - 1)) {
				if (crp->frame == 2)
					pos = 1;
				else if (crp->frame == 3)
					pos = 2;
				else
					pos = 0;
				if ((pos1 % 3) != pos) {
					ValidErr(vsp,SEV_WARNING, ERR_SEQ_FEAT_TranslExceptPhase, "transl_except qual out of frame.");
				}
			}


			cbp = cbp->next;
		}
		
	if (crp->frame > 1) {
		if (! (part_loc & SLP_START)) {
			ValidErr(vsp, SEV_WARNING, ERR_SEQ_FEAT_PartialProblem, "Suspicious CDS location - frame > 1 but not 5' partial");
		} else if ((part_loc & SLP_NOSTART) && (! PartialAtSpliceSite (sfp->location, SLP_NOSTART))) {
			ValidErr(vsp, SEV_INFO, ERR_SEQ_FEAT_PartialProblem, "Suspicious CDS location - frame > 1 and not at consensus splice site");
		}
	}

	if ((part_loc & SLP_START) || (part_prod & SLP_START))
		no_beg = TRUE;

	prot2len = BSLen(newprot);
	len = prot2len;
	BSSeek(newprot, 0, SEEK_SET);
	for (i =0 ; i < len; i++)
	{
		residue1 = BSGetByte(newprot);
		if ((i == 0) && (residue1 == '-'))
			got_dash = TRUE;
		if (residue1 == '*')
		{
			if (i == (len - 1))
				got_stop = TRUE;
			else
				stop_count++;
		}
	}

	if (stop_count)
	{
		if (got_dash)
			ValidErr(vsp, SEV_ERROR, ERR_SEQ_FEAT_StartCodon,
				"Illegal start codon and %ld internal stops. Probably wrong genetic code [%d]",
				(long)stop_count, gccode);
		else
			ValidErr(vsp, SEV_ERROR, ERR_SEQ_FEAT_InternalStop, "%ld internal stops. Genetic code [%d]",
				(long)stop_count, gccode);
		prot_ok = FALSE;
		if (stop_count > 5)
			goto erret;
	} else if (got_dash) {
		ValidErr(vsp, SEV_ERROR, ERR_SEQ_FEAT_StartCodon,
			"Illegal start codon used. Wrong genetic code [%d] or protein should be partial", gccode);
	}

	show_stop = TRUE;

	protid = SeqLocId(sfp->product);
	if (protid != NULL)
	{
		prot1seq = BioseqFind(protid);
		if (prot1seq != NULL)
			prot1len = prot1seq->length;
	}

	if (prot1seq == NULL)
	{
		if (prot2len > 6)
			ValidErr(vsp, SEV_ERROR, ERR_SEQ_FEAT_NoProtein, "No protein Bioseq given");
		goto erret;
	}

	len = prot2len;

	if ((got_stop)&&(len == (prot1len + 1)))  /* ok, got stop */
	{
		len--;
	}

	spp = SeqPortNew(prot1seq, 0, -1, 0, Seq_code_ncbieaa);
	if (spp == NULL) goto erret;

	    /* ignore terminal 'X' from partial last codon if present */

	done = FALSE;
	while ((! done) && (prot1len))
	{
		SeqPortSeek(spp, (prot1len - 1), SEEK_SET);
		residue1 = SeqPortGetResidue(spp);
		if (residue1 == 'X')   /* remove terminal X */
			prot1len--;
		else
			done = TRUE;
	}
	done = FALSE;
	while ((! done) && (len))
	{
		BSSeek(newprot, (len-1), SEEK_SET);
		residue2 = BSGetByte(newprot);
		if (residue2 == 'X')
			len--;
		else
			done = TRUE;
	}

	if (len == prot1len)   /* could be identical */
	{
		SeqPortSeek(spp, 0, SEEK_SET);
	   	BSSeek(newprot, 0, SEEK_SET);
		for (i = 0; i < len; i++)
		{
			residue1 = BSGetByte(newprot);
			residue2 = SeqPortGetResidue(spp);
			if (residue1 != residue2)
			{
				prot_ok = FALSE;
				if (residue2 == INVALID_RESIDUE)
					residue2 = '?';
				if (mismatch == 10)
				{
					ValidErr(vsp, SEV_ERROR, ERR_SEQ_FEAT_MisMatchAA, "More than 10 mismatches. Genetic code [%d]", gccode);
					break;
				}
				else if (i == 0)
				{
					if ((sfp->partial) && (! no_beg) && (! no_end))  /* ok, it's partial */
						ValidErr(vsp, SEV_ERROR, ERR_SEQ_FEAT_PartialProblem, "Start of location should probably be partial");
					else if (residue1 == '-')
						ValidErr(vsp, SEV_ERROR, ERR_SEQ_FEAT_StartCodon, "Illegal start codon used. Wrong genetic code [%d] or protein should be partial", gccode);
					else
					{
						nuclocstr = MapToNTCoords (sfp, protid, i);
						if (nuclocstr != NULL) {
							ValidErr(vsp, SEV_ERROR, ERR_SEQ_FEAT_MisMatchAA,
							"Residue %ld in protein [%c] != translation [%c] at %s",
								(long)(i+1), (char)residue2, (char)residue1, nuclocstr);
						} else {
							ValidErr(vsp, SEV_ERROR, ERR_SEQ_FEAT_MisMatchAA,
							"Residue %ld in protein [%c] != translation [%c]",
								(long)(i+1), (char)residue2, (char)residue1);
						}
						MemFree (nuclocstr);
					}
				}
				else
				{
					nuclocstr = MapToNTCoords (sfp, protid, i);
					if (nuclocstr != NULL) {
						ValidErr(vsp, SEV_ERROR, ERR_SEQ_FEAT_MisMatchAA,
							"Residue %ld in protein [%c] != translation [%c] at %s",
							(long)(i+1), (char)residue2, (char)residue1, nuclocstr);
					} else {
						ValidErr(vsp, SEV_ERROR, ERR_SEQ_FEAT_MisMatchAA,
							"Residue %ld in protein [%c] != translation [%c]",
							(long)(i+1), (char)residue2, (char)residue1);
					}
					MemFree (nuclocstr);
				}
				mismatch++;
			}
		}
		spp = SeqPortFree(spp);
	}
	else
	{
		ValidErr(vsp, SEV_ERROR, ERR_SEQ_FEAT_TransLen, "Given protein length [%ld] does not match translation length [%ld]",
			prot1len, len);
	}

	if ((sfp->partial) && (! mismatch))
	{
		if ((! no_beg) && (! no_end))   /* just didn't label */
		{
			if (! got_stop)
			{
				ValidErr(vsp, SEV_ERROR, ERR_SEQ_FEAT_PartialProblem, "End of location should probably be partial");
			}
			else
			{
				ValidErr(vsp, SEV_ERROR, ERR_SEQ_FEAT_PartialProblem, "This SeqFeat should not be partial");
			}
			show_stop = FALSE;
		}
	}
		

erret:
	if (show_stop)
	{
		if ((! got_stop) && (! no_end))
		{
			ValidErr(vsp, SEV_ERROR, ERR_SEQ_FEAT_NoStop, "Missing stop codon");
		}
		else if ((got_stop) && (no_end))
		{
			ValidErr(vsp, SEV_ERROR, ERR_SEQ_FEAT_PartialProblem, "Got stop codon, but 3'end is labeled partial");
		}
		else if ((got_stop) && (! no_end) && (ragged))
		{
			ValidErr(vsp, SEV_ERROR, ERR_SEQ_FEAT_TransLen, "Coding region extends %d base(s) past stop codon", (int)ragged);
		}
	}

	if (! prot_ok)
	{
		if (transl_except)
			ValidErr(vsp,SEV_WARNING, ERR_SEQ_FEAT_TranslExcept, "Unparsed transl_except qual. Skipped");
	}

	if (prot2seq != NULL)
		BioseqFree(prot2seq);
	else
		BSFree(newprot);
	SeqPortFree(spp);
	return;
}
/*****************************************************************************
*
*   SpliceCheck(sfp)
*      checks for GT/AG rule at splice junctions
*
*****************************************************************************/
#define NOVALUE 0
#define HADGT 1
#define NOGT 2

static void SpliceCheckEx(ValidStructPtr vsp, SeqFeatPtr sfp, Boolean checkAll)
{
	SeqLocPtr slp, nxt, head;
	Uint1 strand = Seq_strand_unknown;
	SeqPortPtr spp=NULL;
	SeqIdPtr last_sip=NULL, sip;
	Int2 total, ctr;
	BioseqPtr bsp = NULL;
	Int4 strt, stp, len = 0, donor, acceptor;
	Int2 residue1, residue2;
	Char tbuf[40];
	Boolean reportAsError, first, last, firstPartial, lastPartial;
	int severity;
	Uint2 partialflag;

	if (sfp == NULL) return;

	if (sfp->excpt)		 /* biological exception */
		return;

	head = sfp->location;
	if (head == NULL) return;

	reportAsError = FALSE;
	if (GetAppProperty ("SpliceValidateAsError") != NULL) {
		reportAsError = TRUE;
	}

	slp = NULL;
	total = 0;
	while ((slp = SeqLocFindPart(head, slp, EQUIV_IS_ONE)) != NULL)
	{
		total++;
		if (slp->choice == SEQLOC_EQUIV)
			return;  /* bail on this one */
		if (total == 1)
			strand = SeqLocStrand(slp);
		else
		{
			if (strand != SeqLocStrand(slp))	 /* bail on mixed strand */
				return;
		}
	}

	if ((! checkAll) && total < 2) return;
	if (total < 1) return;

	slp = NULL;
	ctr = 0;

	first = TRUE;
	last = FALSE;
	firstPartial = FALSE;
	lastPartial = FALSE;

	slp = SeqLocFindPart(head, slp, EQUIV_IS_ONE);
	while (slp != NULL)
	{
		nxt = SeqLocFindPart(head, slp, EQUIV_IS_ONE);
		last = (Boolean) (nxt == NULL);
		partialflag = SeqLocPartialCheck (slp);
		firstPartial = (Boolean) (first && (partialflag & SLP_START));
		lastPartial = (Boolean) (last && (partialflag & SLP_STOP));
		ctr++;
		sip = SeqLocId(slp);
		if (sip == NULL) break;
		if ((ctr == 1) || (! SeqIdMatch(sip, last_sip)))
		{
			spp = SeqPortFree(spp);
			bsp = BioseqLockById(sip);
			if (bsp == NULL) break;
			len = bsp->length;
			spp = SeqPortNew(bsp, 0, -1, strand, Seq_code_ncbi4na);
			BioseqUnlock(bsp);
			if (spp == NULL) break;
			last_sip = sip;
		}
		acceptor = SeqLocStart(slp);
		donor = SeqLocStop(slp);

		if (strand != Seq_strand_minus)
		{
			strt = acceptor;
			stp = donor;
		}
		else
		{
			strt = donor;
			donor = acceptor;
			acceptor = strt;
			stp = len - donor - 1;	/* orient to reverse complement seqport */
			strt = len - acceptor - 1;
		}

		if (((checkAll && (! lastPartial)) || ctr < total) && (stp < (len - 2)))   /* check donor on all but last exon and on sequence */
		{
			SeqPortSeek(spp, (stp+1), SEEK_SET);
			residue1 = SeqPortGetResidue(spp);
			residue2 = SeqPortGetResidue(spp);
			if (IS_residue(residue1) && IS_residue(residue2))
			{
				if ((! (residue1 & 4)) ||        /* not G or */
					(! (residue2 & 8)))          /* not T */
				{
				if ((residue1 & 4) && (residue2 & 2)) { /* GC minor splice site */
					tbuf[39] = '\0';
					if (vsp->suppressContext) {
						WorstBioseqLabel(bsp, tbuf, 39, OM_LABEL_CONTENT);
					} else {
						BioseqLabel(bsp, tbuf, 39, OM_LABEL_CONTENT);
					}
					ValidErr(vsp, SEV_WARNING, ERR_SEQ_FEAT_NotSpliceConsensus,
						"Rare splice donor consensus (GC) found instead of (GT) after exon ending at position %ld of %s",
				    	  (long)(donor+1), tbuf);
				} else {
					if (checkAll) {
					  severity = SEV_WARNING;
					} else if (reportAsError) {
					  severity = SEV_ERROR;
					} else {
					  severity = SEV_WARNING;
					}
					tbuf[39] = '\0';
					if (vsp->suppressContext) {
						WorstBioseqLabel(bsp, tbuf, 39, OM_LABEL_CONTENT);
					} else {
						BioseqLabel(bsp, tbuf, 39, OM_LABEL_CONTENT);
					}
					ValidErr(vsp, severity, ERR_SEQ_FEAT_NotSpliceConsensus,
						"Splice donor consensus (GT) not found after exon ending at position %ld of %s",
				     	 (long)(donor+1), tbuf);
				}
				}
			}
		}

		if (((checkAll && (! firstPartial))  || ctr != 1) && (strt > 1))
		{
			SeqPortSeek(spp, (strt - 2), SEEK_SET);
			residue1 = SeqPortGetResidue(spp);
			residue2 = SeqPortGetResidue(spp);
			if (IS_residue(residue1) && IS_residue(residue2))
			{
				if ((! (residue1 & 1)) ||        /* not A or */
					(! (residue2 & 4)))          /* not G */
				{
					if (checkAll) {
					  severity = SEV_WARNING;
					} else if (reportAsError) {
					  severity = SEV_ERROR;
					} else {
					  severity = SEV_WARNING;
					}
					tbuf[39] = '\0';
					if (vsp->suppressContext) {
						WorstBioseqLabel(bsp, tbuf, 39, OM_LABEL_CONTENT);
					} else {
						BioseqLabel(bsp, tbuf, 39, OM_LABEL_CONTENT);
					}
					ValidErr(vsp, severity, ERR_SEQ_FEAT_NotSpliceConsensus,
					"Splice acceptor consensus (AG) not found before exon starting at position %ld of %s",
				      (long)(acceptor+1), tbuf);
				}
			}
		}
		first = FALSE;
		slp = nxt;
	}

	SeqPortFree(spp);
    return;
}

NLM_EXTERN void SpliceCheck(ValidStructPtr vsp, SeqFeatPtr sfp)

{
	SpliceCheckEx (vsp, sfp, FALSE);
}

/*****************************************************************************
*
*   ValidateSeqLoc(vsp, slp, prefix)
*
*****************************************************************************/
NLM_EXTERN void ValidateSeqLoc(ValidStructPtr vsp, SeqLocPtr slp, CharPtr prefix)
{
	SeqLocPtr tmp, prev;
	Boolean retval = TRUE, tmpval, mixed_strand = FALSE, ordered=TRUE;
	CharPtr ctmp;
	Uint1 strand2, strand1;
	SeqIntPtr sip1, sip2;
	SeqPntPtr spp;
	PackSeqPntPtr pspp;
	SeqIdPtr id1 = NULL, id2;

	if (slp == NULL) return;

	tmp = NULL;
	prev = NULL;
	sip1 = NULL;
	strand1 = Seq_strand_other;
	while ((tmp = SeqLocFindNext(slp, tmp)) != NULL)
	{
		tmpval = TRUE;
		switch (tmp->choice)
		{
			case SEQLOC_INT:
				sip2 = (SeqIntPtr)(tmp->data.ptrvalue);
				strand2 = sip2->strand;
				id2 = sip2->id;
				tmpval = SeqIntCheck (sip2);
				if ((tmpval) && (sip1 != NULL) && (ordered))
				{
					if (SeqIdForSameBioseq(sip1->id, sip2->id))
					{
						if (strand2 == Seq_strand_minus)
						{
							if (sip1->to < sip2->to)
								ordered = FALSE;
						}
						else if (sip1->to > sip2->to)
							ordered = FALSE;
					}
				}
				break;
			case SEQLOC_PNT:
				spp = (SeqPntPtr)(tmp->data.ptrvalue);
				strand2 = spp->strand;
				id2 = spp->id;
				tmpval = SeqPntCheck (spp);
				break;
			case SEQLOC_PACKED_PNT:
				pspp = (PackSeqPntPtr)(tmp->data.ptrvalue);
				strand2 = pspp->strand;
				id2 = pspp->id;
				tmpval = PackSeqPntCheck (pspp);
				break;
			default:
			    strand2 = Seq_strand_other;
				id2 = NULL;
				break;
		}
		if (! tmpval)
		{
			retval = FALSE;
			ctmp = SeqLocPrint(tmp);
			ValidErr(vsp, SEV_FATAL, ERR_SEQ_FEAT_Range, "%s: SeqLoc [%s] out of range", prefix, ctmp);
			MemFree(ctmp);

		}

		if ((strand1 != Seq_strand_other) && (strand2 != Seq_strand_other))
		{
			if (SeqIdForSameBioseq(id1, id2))
			{
				if (strand1 != strand2)
					mixed_strand = TRUE;
			}
		}

		strand1 = strand2;
		id1 = id2;
	}

	if ((mixed_strand) || (! ordered))
	{
		ctmp = SeqLocPrint(slp);
		if (mixed_strand)
			ValidErr(vsp, SEV_ERROR, ERR_SEQ_FEAT_MixedStrand, "%s: Mixed strands in SeqLoc [%s]", prefix, ctmp);
		if (! ordered)
			ValidErr(vsp, SEV_ERROR, ERR_SEQ_FEAT_SeqLocOrder, "%s: Intervals out of order in SeqLoc [%s]", prefix, ctmp);
		MemFree(ctmp);
	}

	return;
}

/*****************************************************************************
*
*   PatchBadSequence(bsp)
*
*****************************************************************************/
NLM_EXTERN Boolean PatchBadSequence(BioseqPtr bsp)
{
	ByteStorePtr newseq;
	SeqPortPtr spp;
	Boolean is_na;
	Uint1 seqcode;
	Int2 repchar, residue;
	Int4 i, len;

	if (bsp == NULL) return FALSE;
	if (! ((bsp->repr == Seq_repr_raw) || (bsp->repr == Seq_repr_const)))
		return FALSE;

	is_na = ISA_na(bsp->mol);
	if (is_na)
	{
		seqcode = Seq_code_iupacna;
		repchar = (Int2)'N';   /* N */
	}
	else
	{
		seqcode = Seq_code_iupacaa;
		repchar = (Int2)'X';
	}

	spp = SeqPortNew(bsp, 0, -1, 0, seqcode);
	if (spp == NULL) return FALSE;

	len = bsp->length;
	newseq = BSNew(len);
	if (newseq == NULL)
	{
		SeqPortFree(spp);
		return FALSE;
	}

	for (i = 0; i < len; i++)
	{
		residue = SeqPortGetResidue(spp);
		if (residue == INVALID_RESIDUE)
		{
			residue = repchar;
		}
		BSPutByte(newseq, residue);
	}

	SeqPortFree(spp);
	BSFree(bsp->seq_data);
	bsp->seq_data = newseq;
	bsp->seq_data_type = seqcode;

	BioseqRawPack(bsp);

	return TRUE;
}

static void FindABioseq(SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
	BioseqPtr PNTR bp;
	BioseqPtr bsp;

	bp = (BioseqPtr PNTR)data;
	if (*bp != NULL)   /* already got one */
		return;

	if (IS_Bioseq(sep))
	{
		bsp = (BioseqPtr)(sep->data.ptrvalue);
		*bp = bsp;
	}
	return;
}

NLM_EXTERN CharPtr FindIDForEntry (SeqEntryPtr sep, CharPtr buf)
{
	BioseqPtr bsp = NULL;
	
	if ((sep == NULL) || (buf == NULL))
		return NULL;

	*buf = '\0';
	SeqEntryExplore (sep, (Pointer)(&bsp), FindABioseq);

	if (bsp == NULL) return NULL;

	SeqIdPrint(bsp->id, buf, PRINTID_FASTA_LONG);
	return buf;
}

static CharPtr TrimSpacesOnEitherSide (CharPtr str)

{
  Uchar    ch;
  CharPtr  dst;
  CharPtr  ptr;

  if (str != NULL && str [0] != '\0') {
    dst = str;
    ptr = str;
    ch = *ptr;
    while (ch != '\0' && ch <= ' ') {
      ptr++;
      ch = *ptr;
    }
    while (ch != '\0') {
      *dst = ch;
      dst++;
      ptr++;
      ch = *ptr;
    }
    *dst = '\0';
    dst = NULL;
    ptr = str;
    ch = *ptr;
    while (ch != '\0') {
      if (ch != ' ') {
        dst = NULL;
      } else if (dst == NULL) {
        dst = ptr;
      }
      ptr++;
      ch = *ptr;
    }
    if (dst != NULL) {
      *dst = '\0';
    }
  }
  return str;
}

static void CopyLetters (CharPtr dest, CharPtr source, size_t maxsize)

{
  Char     ch;
  CharPtr  tmp;

  if (dest == NULL || maxsize < 1) return;
  *dest = '\0';
  if (source == NULL) return;
  maxsize--;
  tmp = dest;
  ch = *source;
  while (maxsize > 1 && ch != '\0') {
    if (ch != '.') {
      *dest = ch;
      dest++;
      maxsize--;
    }
    source++;
    ch = *source;
  }
  *dest = '\0';
  TrimSpacesOnEitherSide (tmp);
}

static void LookForEtAl (ValidStructPtr vsp, ValNodePtr tmp)

{
  AuthorPtr    ap;
  AuthListPtr  authors = NULL;
  CitArtPtr    cap;
  CitBookPtr   cbp;
  CitGenPtr    cgp;
  CitSubPtr    csp;
  Char         first [64];
  Char         initials [16];
  Char         last [64];
  ValNodePtr   names;
  NameStdPtr   nsp;
  PersonIdPtr  pid;

  if (vsp == NULL || tmp == NULL) return;
  switch (tmp->choice) {
		case PUB_Article:
			cap = (CitArtPtr)(tmp->data.ptrvalue);
			authors = cap->authors;
			break;
		case PUB_Man:
		case PUB_Book:
		case PUB_Proc:
			cbp = (CitBookPtr)(tmp->data.ptrvalue);
			authors = cbp->authors;
			break;
		case PUB_Gen:
			cgp = (CitGenPtr)(tmp->data.ptrvalue);
			authors = cgp->authors;
			break;
		case PUB_Sub:
			csp = (CitSubPtr)(tmp->data.ptrvalue);
			authors = csp->authors;
			break;
		default:
			break;
  }
  if (authors == NULL || authors->choice != 1) return;
  for (names = authors->names; names != NULL; names = names->next) {
    ap = names->data.ptrvalue;
    if (ap != NULL) {
      pid = ap->name;
      if (pid != NULL && pid->choice == 2) {
        nsp = pid->data;
        if (nsp != NULL && nsp->names [0] != NULL) {
          CopyLetters (last, nsp->names [0], sizeof (last));
          CopyLetters (first, nsp->names [1], sizeof (first));
          CopyLetters (initials, nsp->names [4], sizeof (initials));
          if ((StringICmp (last, "et al") == 0) ||
              (StringCmp (initials, "al") == 0 &&
               StringCmp (last, "et") == 0 &&
               first [0] == '\0')) {
            if (names->next == NULL) {
              ValidErr (vsp, SEV_WARNING, ERR_GENERIC_AuthorListHasEtAl,
                        "Author list ends in et al.");
            } else {
              ValidErr (vsp, SEV_WARNING, ERR_GENERIC_AuthorListHasEtAl,
                        "Author list contains et al.");
            }
          }
        }
      }
    }
  }
}

static void SpellCheckPub(ValidStructPtr vsp, ValNodePtr tmp)
{
	CitArtPtr cap;
	CitBookPtr cbp;
	CitGenPtr cgp;
	ValNodePtr titles = NULL;

	if ((vsp == NULL) || (tmp == NULL))
		return;

	switch (tmp->choice)
	{
		case PUB_Article:
			cap = (CitArtPtr)(tmp->data.ptrvalue);
			titles = cap->title;
			break;
		case PUB_Man:
		case PUB_Book:
		case PUB_Proc:
			cbp = (CitBookPtr)(tmp->data.ptrvalue);
			titles = cbp->title;
			break;
		case PUB_Gen:
			cgp = (CitGenPtr)(tmp->data.ptrvalue);
			if (cgp->cit != NULL)
				SpellCheckString(vsp, cgp->cit);
			if (cgp->title != NULL)
				SpellCheckString(vsp, cgp->title);
			break;
		default:
			break;
	}

	if (titles != NULL)
	{
		for (; titles != NULL; titles = titles->next)
		{
			if (titles->choice == Cit_title_name)
				SpellCheckString(vsp, (CharPtr)(titles->data.ptrvalue));
		}
	}

	return;
}

static void SpellCheckSeqDescr(GatherContextPtr gcp)

{
	PubdescPtr pdp;
	ValNodePtr tmp, vnp;
	ValidStructPtr vsp;

	vsp = (ValidStructPtr)(gcp->userdata);
	if (vsp == NULL)
		return;

	vnp = (ValNodePtr)(gcp->thisitem);
	if (vnp == NULL)
		return;

	vsp->descr = vnp;
	vsp->sfp = NULL;

	if (vnp->choice == Seq_descr_pub) {
		pdp = (PubdescPtr)(vnp->data.ptrvalue);
		for (tmp = pdp->pub; tmp != NULL; tmp = tmp->next) {
			LookForEtAl (vsp, tmp);
		}
	}

	if (vsp->spellfunc == NULL) return;

	switch (vnp->choice)
	{
		case Seq_descr_title:
		case Seq_descr_region:
		case Seq_descr_comment:
			SpellCheckString(vsp, (CharPtr)(vnp->data.ptrvalue));
			break;
		case Seq_descr_pub:
			pdp = (PubdescPtr)(vnp->data.ptrvalue);
			for (tmp = pdp->pub; tmp != NULL; tmp = tmp->next)
			{
				SpellCheckPub(vsp, tmp);
			}
			break;
		default:
			break;
	}
	return;
}

NLM_EXTERN void SpellCheckSeqFeat(GatherContextPtr gcp)
{
	PubdescPtr pdp;
	SeqFeatPtr sfp;
	ProtRefPtr prp;
	ValidStructPtr vsp;
	ValNodePtr vnp;

	vsp = (ValidStructPtr)(gcp->userdata);
	if (vsp == NULL)
		return;

	sfp = (SeqFeatPtr)(gcp->thisitem);
	if (sfp == NULL)
		return;

	vsp->descr = NULL;
	vsp->sfp = sfp;

	if (sfp->data.choice == SEQFEAT_PUB) {
		pdp = (PubdescPtr)(sfp->data.value.ptrvalue);
		for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
			LookForEtAl (vsp, vnp);
		}
	}

	if (vsp->spellfunc == NULL) return;

	SpellCheckString(vsp, sfp->comment);

	switch (sfp->data.choice)
	{
		case 1:        /* Gene-ref */
			break;
		case 2:        /* Org-ref */
			break;
		case 3:        /* Cdregion */
			break;
		case 4:        /* Prot-ref */
			prp = (ProtRefPtr)(sfp->data.value.ptrvalue);
			for (vnp = prp->name; vnp != NULL; vnp = vnp->next)
				SpellCheckString(vsp, (CharPtr)(vnp->data.ptrvalue));
			SpellCheckString(vsp, prp->desc);
			break;
		case 5:        /* RNA-ref */
			break;
		case 6:        /* Pub */
				pdp = (PubdescPtr)(sfp->data.value.ptrvalue);
				for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next)
				{
					SpellCheckPub(vsp, vnp);
				}
			break;
		case 7:        /* Seq */
			break;
		case 8:        /* Imp-feat */
			break;
		case 9:        /* Region */
			SpellCheckString(vsp, (CharPtr)(sfp->data.value.ptrvalue));
			break;
		case 10:        /* Comment */
			break;
		case 11:        /* Bond */
			break;
		case 12:        /* Site */
			break;
		case 13:        /* Rsite-ref */
			break;
		case 14:        /* User-object */
			break;
		case 15:        /* TxInit */
			break;
		case 16:        /* Numbering */
			break;
		case 17:        /* Secondary Structure */
			break;
		case 18:        /* NonStdRes*/
			break;
		case 19:        /* Heterogen*/
			break;
		case 20:        /* BioSource*/
			break;
		default:
			break;
	}

	return;
}

NLM_EXTERN void SpellCheckString (ValidStructPtr vsp, CharPtr str)
{
	if ((vsp == NULL) || (str == NULL))
		return;

	if (vsp->spellfunc == NULL) return;

	(* (vsp->spellfunc))((char *)str, (vsp->spellcallback));

	return;
}

NLM_EXTERN void SpellCallBack (char * str)
{
	ErrSev sev;

	sev = SEV_ERROR;
	if (globalvsp != NULL && globalvsp->justwarnonspell) {
		sev = SEV_WARNING;
	}
	ValidErr(globalvsp, sev, ERR_GENERIC_Spell, "[ %s ]", (CharPtr)str);
	return;
}

