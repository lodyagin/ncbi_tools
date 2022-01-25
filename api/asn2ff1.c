/*   asn2ff1.c
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*            National Center for Biotechnology Information (NCBI)
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government do not place any restriction on its use or reproduction.
*  We would, however, appreciate having the NCBI and the author cited in
*  any work or product based on this material
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
* ===========================================================================
*
* File Name:  asn2ff1.c
*
* Author:  Karl Sirotkin, Tom Madden, Tatiana Tatusov
*
* Version Creation Date:   7/15/95
*
* $Revision: 6.41 $
*
* File Description:  files that go with "asn2ff"
*
* Modifications:  
* --------------------------------------------------------------------------
* $Log: asn2ff1.c,v $
* Revision 6.41  1998/09/28 18:17:12  bazhin
* Added new function
* "CharPtr PNTR SeqEntryToStrArray(SeqEntryPtr sep, Uint1 format, Uint1 mode)",
* which works like SeqEntryToFlat(), but prints flat entry into memory
* as an array of strings instead of printing into file on disc. It returns
* that array.
*
* Revision 6.40  1998/09/23 23:47:00  tatiana
* MI_TECH_fli_cdna added to keywords
*
* Revision 6.39  1998/09/15 16:04:51  tatiana
* removed redundant MemFree in seqEntryToFlatAjp
*
* Revision 6.38  1998/09/14 16:37:57  tatiana
* SeqEntryToFlatAjp added
*
* Revision 6.37  1998/09/11 19:10:58  tatiana
* set error_msgs to TRUE
*
* Revision 6.36  1998/09/08 20:41:42  tatiana
* made non-static Asn2ffJobCreate
*
* Revision 6.35  1998/08/21 16:56:17  shavirin
* Added new function SeqEntryToGBFlatNoSeq()
*
* Revision 6.34  1998/07/27 19:58:22  kans
* SeqEntryToFlatEx was not returning at proper place when nesting pop/phy/mut sets
*
* Revision 6.33  1998/07/16 16:06:51  kans
* use ObjMgrGetEntityIDForChoice instead of ObjMgrGetEntityIDForPointer for SeqEntryPtr
*
* Revision 6.32  1998/07/16 14:44:48  kans
* handles segmented sets within pop/phy sets (Tatiana)
*
* Revision 6.31  1998/07/15 22:07:12  kans
* implemented sequence manager indexes for non-segmented nucleotides
*
* Revision 6.30  1998/06/15 14:56:38  tatiana
* UNIX compiler warnings fixed
*
* Revision 6.29  1998/05/11 21:58:15  tatiana
* some functions moved to asn2ff6.c
*
* Revision 6.28  1998/05/08 21:54:58  tatiana
* SeqEntryToPartRpt() added
*
* Revision 6.27  1998/05/05 16:56:45  kans
* MakeBaseLocAwp needed to initialize tsip
*
* Revision 6.26  1998/05/05 15:35:28  tatiana
* added SEQID_OTHER to MakeLocus()
*
* Revision 6.25  1998/04/30 21:36:48  tatiana
* *** empty log message ***
*
* Revision 6.22  1998/04/15 19:11:49  tatiana
*  GetAppProperty added
*
* Revision 6.21  1998/03/27 23:04:30  tatiana
* memory leaks cleanup
*
* Revision 6.20  1998/03/24 16:42:57  tatiana
* a bug fixed in SeqSubmitToFlat()
*
* Revision 6.19  1998/03/20 03:07:51  kans
* genpept in sequin_mode can use local ID for locus
*
* Revision 6.18  1998/03/09 21:40:40  tatiana
* accession length increased to 60
*
* Revision 6.16  1998/02/12 14:15:53  kans
* switch to StringNCpy_0 was truncating division code in locus line to two characters
*
* Revision 6.15  1998/02/11 19:29:59  tatiana
* cleaning memory leaks
*
* Revision 6.14  1998/02/11 18:27:20  tatiana
* StringNCpy changed to StringNCpy_0
*
* Revision 6.12  1998/02/03 21:28:07  tatiana
* fixed SeqLocToFlat and SeqEntryToFlatEx
*
* Revision 6.11  1998/01/16 19:00:02  tatiana
* improved the perfomance in SeqEntryToFlatEx
*
* Revision 6.10  1998/01/16 18:27:06  tatiana
* fix mol-type in LOCUS line
*
* Revision 6.8  1998/01/08 23:25:16  tatiana
* restore the usage of -h parameter
*
* Revision 6.7  1998/01/06 23:52:51  tatiana
* fixed start position in CheckSeqPort()
*
* Revision 6.6  1997/12/30 21:46:05  tatiana
* fixed scRNA in PrintLocus()
*
* Revision 6.4  1997/11/10 18:04:52  tatiana
* changes in SeqLocToFlat
*
* Revision 6.3  1997/11/03 20:47:12  shavirin
* Removed memory leak
*
* Revision 6.2  1997/10/23 16:56:15  tatiana
* changes in get_pub allow to gather CitSub directly from SeqSubmit
*
* Revision 6.1  1997/09/12 20:03:31  tatiana
* added source feature in genome_view
*
* Revision 6.0  1997/08/25 18:04:38  madden
* Revision changed to 6.0
*
* Revision 5.57  1997/08/06 22:47:35  tatiana
*  changes in CheckSeqPort() for printing a region
*
* Revision 5.55  1997/07/29 20:52:43  tatiana
* Array bounds write fixed in GetFlatRetract()
*
* Revision 5.54  1997/07/29 16:16:34  tatiana
* *** empty log message ***
*
* Revision 5.53  1997/07/29 15:50:07  tatiana
* SeqEntryToFlatEx will work for GENPEPT format, asn2gp_setup modofied
*
* Revision 5.52  1997/07/28 19:03:53  vakatov
* [WIN32,MSVC++]  Restored lost "NCBIOBJ.LIB" pro-DLL modifications
*
* Revision 5.51  1997/07/24 18:58:58  tatiana
* memory corruption fixed in asn2ff_cleanup()
*
* Revision 5.50  1997/07/24 16:50:19  tatiana
*  fixed bugs in sequence printing
*
* Revision 5.48  1997/07/16 20:52:37  tatiana
* PrintGenome() moved to wprint.c
*
* Revision 5.45  1997/06/19 18:36:52  vakatov
* [WIN32,MSVC++]  Adopted for the "NCBIOBJ.LIB" DLL'ization
*
 * Revision 5.42  1997/03/13  15:44:59  tatiana
 * added asn2hp_setup
 *
 * Revision 5.40  1997/02/03  15:23:43  tatiana
 * a bug fixed in asn2ff_print (AsnIoClose removed)
 *
 * Revision 5.39  1997/01/31  17:18:50  tatiana
 * SeqSubmitToFlat changed for EBI
 *
 * Revision 5.38  1997/01/27  19:13:39  tatiana
 * more changes to SeqSubmitToFlat()
 *
 * Revision 5.37  1997/01/27  18:34:14  tatiana
 * SeqSubmitToFlat() changed to produce EBI submissions
 *
 * Revision 5.36  1997/01/16  23:04:11  tatiana
 * a typo fixed in PrintOrganismLine()
 *
 * Revision 5.35  1997/01/13  22:34:15  tatiana
 * show_gene = TRUE
 *
 * Revision 5.34  1997/01/07  20:54:38  tatiana
 * asn2ff_setup fixed for targeted bioseq in segmented set
 *
 * Revision 5.32  1996/12/17  22:52:19  tatiana
 * added AddKeywords()
 *
 * Revision 5.31  1996/11/19  22:46:48  tatiana
 * global Boolean Template_load added
 *
 * Revision 5.29  1996/11/01  17:51:53  tatiana
 * GetDefline added to embl and GenPept formats
 *
 * Revision 5.27  1996/10/25  22:10:12  tatiana
 * HTG division is legal
 *
 * Revision 5.26  1996/09/18  20:20:18  tatiana
 * positions fixed in PrintGenome
 *
 * Revision 5.25  1996/09/17  14:58:10  tatiana
 * SeqSubmitToFlat needs show_gene argument
 *
 * Revision 5.24  1996/09/12  17:51:29  tatiana
 * a bug fixed in PrintSourceFeat
 *
 * Revision 5.23  1996/09/09  13:36:02  kans
 * moved BioseqGetGBDivCode from toasn.[ch] to asn2ff.h/asn2ff6.c
 *
 * Revision 5.22  1996/09/06  21:05:10  tatiana
 * change name GBGetDivision to BioseqGetGBDivCode
 *
 * Revision 5.21  1996/09/06  20:56:34  tatiana
 * GetDivision changed to call new function BioseqGetGBDivCode
 *
 * Revision 5.20  1996/09/03  19:49:25  tatiana
 * asn2ff_cleanup changed to free new_loc from Gather
 *
 * Revision 5.19  1996/08/27  22:51:40  tatiana
 * add ajp->only_one for gathering pubs
 *
 * Revision 5.18  1996/08/27  22:11:53  tatiana
 * change GetDivision to keep PAT and SYN from GBBlock
 *
 * Revision 5.17  1996/08/27  19:12:40  tatiana
 * calls SeqIdSelect before SeqIdWrite in GetLocusPartsAwp to get the best ID for accession number and locus name
 *
 * Revision 5.16  1996/08/22  18:46:08  tatiana
 * bug fixed
 *
 * Revision 5.15  1996/08/16  20:31:14  tatiana
 * CreateDefLine() call added
 *
 * Revision 5.14  1996/08/09  17:30:22  tatiana
 * ErrPostEx changed to ErrpostStr for const strings
 *
 * Revision 5.13  1996/08/06  20:30:46  kans
 * SeqIdFindBest called to handle local IDs and genbank IDs coexisting
 *
 * Revision 5.12  1996/07/31  19:06:13  tatiana
 * empty KEYWORD are not mapped to GBBlock in sequin
 *
 * Revision 5.11  1996/07/31  16:31:55  tatiana
 * fix in PrintGBOrganismLine
 *
 * Revision 5.10  1996/07/31  15:23:24  tatiana
 * minor change in PrintDefifnitionLine
 *
 * Revision 5.9  1996/07/30  16:36:40  tatiana
 * PrintDefinitionLine changed for htgs
 *
 * Revision 5.7  1996/07/19  21:37:42  tatiana
 * HTG keywords and deflines added
 *
 * Revision 5.4  1996/07/03  20:59:29  tatiana
 * need to free gbp->descr in GetDivision
 *
 * Revision 5.3  1996/07/02  19:42:52  tatiana
 * support for delta sequence added
 *
 * Revision 5.2  1996/06/14  18:02:56  tatiana
 * GetDivision changes
 *
 * Revision 5.1  1996/06/11  15:26:00  tatiana
 * add PrintNID
 *
 * Revision 4.43  1996/05/16  20:55:44  tatiana
 * source_info added to GBEntry structure
 *
 * Revision 4.42  1996/05/06  16:09:14  tatiana
 * a bug fixed in PrintKeyword()
 *
 * Revision 4.41  1996/05/02  20:32:51  tatiana
 * muid from PUB_Medline added
 *
 * Revision 4.40  1996/05/02  17:42:52  tatiana
 * GetSubmitDescr() added that will show CitSub.descr in REMARK
 *
 * Revision 4.39  1996/04/29  18:47:42  tatiana
 * create independent paragraph for each comment block
 *
 * Revision 4.39  1996/04/29  18:47:42  tatiana
 * create independent paragraph for each comment block
 *
 * Revision 4.38  1996/04/15  14:35:35  tatiana
 * memory leaks cleaning
 *
 * Revision 4.37  1996/04/12  03:42:23  tatiana
 *  : a bug fixed
 *
 * Revision 4.36  1996/04/10  22:50:47  kans
 * source and organism lines split, feature line has null gbp->descr,
 * comment and region, etc., in separate paragraphs (TT)
 *
 * Revision 4.34  1996/04/09  14:03:01  tatiana
 * print COMMENT in blocks
 *
 * Revision 4.33  1996/03/25  22:24:59  tatiana
 * a bug fixed in ValidateAccession
 *
 * Revision 4.32  1996/03/25  17:47:28  tatiana
 * 2+6 accession handling
 *
 * Revision 4.31  1996/03/18  23:37:33  tatiana
 * fix a bug in PrintOrganism to take correct taxin id for www hotlink
 *
 * Revision 4.30  1996/03/12  21:35:41  tatiana
 * GetRetract() added to handle Medline erratum
 *
 * Revision 4.29  1996/02/28  04:53:06  ostell
 * changes to support segmented master seeuquences
 *
 * Revision 4.28  1996/02/26  03:45:24  ostell
 * fixed usage of GatherDescrListByChoice and changed function to allocate
 * first DescrStruct instead using static storage
 *
 * Revision 4.25  1996/02/18  21:14:53  tatiana
 * memory leaks cleaned up, GetPubNum() added
 *
 * Revision 4.23  1996/02/15  15:50:07  tatiana
 * Gather for temp load items added
 *
 * Revision 4.21  1995/12/20  22:37:01  tatiana
 * NID field turned on!
 *
 * Revision 4.20  1995/12/18  20:58:54  tatiana
 * a bug fixed in CheckXrefLine (using DescrStruct)
 *
 * Revision 4.19  1995/12/15  19:38:08  kans
 * bioseq is selectable in flat file (TT)
 *
 * Revision 4.18  1995/12/14  19:13:33  kans
 * fixed GatherDescrListByChoice
 *
 * Revision 4.17  1995/12/13  16:30:04  tatiana
 * itemID ... added to descriptors
 *
 * Revision 4.15  1995/12/05  22:18:30  tatiana
 * a bug fixed in CheckXrefLine()
 *
 * Revision 4.14  1995/12/05  17:16:49  tatiana
 * bug fixed.
 *
 * Revision 4.13  1995/12/05  16:35:57  kans
 * ajp->asn2ffwep could dangle in asn2ff_setup and asn2ff_cleanup (TT)
 *
 * Revision 4.12  1995/11/29  15:47:07  tatiana
 * a big fixed in GetDivision
 *
 * Revision 4.11  1995/11/22  18:58:56  tatiana
 * memory leaks cleanup
 *
 * Revision 4.10  1995/11/21  17:51:35  tatiana
 * changes in FFPrintArray
 *
 * Revision 4.9  1995/11/17  21:28:35  kans
 * asn2ff now uses gather (Tatiana)
 *
 * Revision 4.1  1995/07/31  19:46:39  tatiana
 * accept /transl_table from taxon in RELEASE mode
 *
 * Revision 1.56  1995/07/17  19:33:20  kans
 * parameters combined into Asn2ffJobPtr structure
 *
 * Revision 1.46  1995/06/19  21:40:02  kans
 * Tatiana's first major reorganization, moving printing, adding HTML
 *
 * Revision 1.44  1995/05/22  14:51:08  tatiana
 * remove old MatchNAFeatToGene function and add ASN2FF_SHOW_ALL_PUBS
 *
 * Revision 1.43  1995/05/19  21:25:06  kans
 * gene match code moved to sequtil (ostell)
* ==========================================================================
*/

#include <seqmgr.h>
#include <gather.h>
#include <asn2ffg.h>
#include <asn2ffp.h>
#include <a2ferrdf.h>
#include <parsegb.h>
#include <gbfeat.h>
#include <ffprint.h>
#include <tofasta.h>
#include <subutil.h>
#include <explore.h>

#ifdef ENABLE_ENTREZ
#include <accentr.h>
#endif

/* The following corresponds to NUM_SEQ_LINES lines each with 60 
residues/basepairs */
#define SEQ_BLK_SIZE (60*NUM_SEQ_LINES)
#define A2F_OTHER ( (Uint1)0)
#define A2F_SOURCE_FEATURE ( (Uint1)1)
#define A2F_FEATURE ( (Uint1)2)
#define A2F_REFERENCE ( (Uint1)3)
#define A2F_FEATURE_NEW ( (Uint1)4)
#define A2F_COMMENT ( (Uint1)5)

/* static Int2 pap_total; -- UNUSED */
static Boolean Template_load = FALSE;

/* ---------------Function Prototypes ---------------*/
Int4 asn2pr_setup PROTO ((Asn2ffJobPtr ajp, FFPrintArrayPtr PNTR papp));
Int4 asn2hp_setup PROTO ((Asn2ffJobPtr ajp, FFPrintArrayPtr PNTR papp));
Int4 asn2gb_setup PROTO ((Asn2ffJobPtr ajp, FFPrintArrayPtr PNTR papp));
Int4 asn2embl_setup PROTO ((Asn2ffJobPtr ajp, FFPrintArrayPtr PNTR papp));
Int4 asn2gp_setup PROTO ((Asn2ffJobPtr ajp, FFPrintArrayPtr PNTR papp));
Int4 asn2ep_setup PROTO ((Asn2ffJobPtr ajp, FFPrintArrayPtr PNTR papp));
void LoadPap PROTO ((FFPrintArrayPtr pap, FFPapFct fct, Asn2ffJobPtr ajp, Int4 index, Uint1 last, Uint1 printxx, Int2 estimate, Uint1 element_type, GBEntryPtr gbp));

void CheckSeqPort PROTO ((Asn2ffJobPtr ajp, GBEntryPtr gbp, Int4 start));
void PrintGenome PROTO ((Asn2ffJobPtr ajp, GBEntryPtr gbp));
void GetMolInfo PROTO ((Asn2ffJobPtr ajp, CharPtr buffer, GBEntryPtr gbp));
CharPtr GetPDBSourceLine PROTO ((PdbBlockPtr pdb));
void PrintDateLines PROTO ((Asn2ffJobPtr ajp, GBEntryPtr gbp));
void PrintXrefLine PROTO ((Asn2ffJobPtr ajp, GBEntryPtr gbp));
Boolean CheckXrefLine PROTO ((Asn2ffJobPtr ajp, GBEntryPtr gbp));
void set_flags PROTO ((Asn2ffJobPtr ajp));
void PrintSeqBlk PROTO ((Asn2ffJobPtr ajp, GBEntryPtr gbp));


#define  TOTAL_ESTKW         11
#define  TOTAL_STSKW         5
#define  TOTAL_GSSKW         2

static CharPtr EST_kw_array[TOTAL_ESTKW] = {
     "EST", "EST PROTO((expressed sequence tag)", "expressed sequence tag",
     "EST (expressed sequence tag)", "EST(expressed sequence tag)",
     "partial cDNA sequence", "transcribed sequence fragment", "TSR",
     "putatively transcribed partial sequence", "UK putts"
     };

static CharPtr GSS_kw_array[TOTAL_GSSKW] = {
     "GSS", "trapped exon"
     };
static CharPtr STS_kw_array[TOTAL_STSKW] = {
     "STS", "STS(sequence tagged site)", "STS (sequence tagged site)", 
     "STS sequence", "sequence tagged site"
     };

static Int2 MatchArrayString(CharPtr array_string[], Int2 totalstr, CharPtr text)
{
   Int2 i;

   for (i = 0; i < totalstr && text != NULL; i++)
       if (StringCmp(array_string[i], text) == 0)
          return (i);

   return (-1);

} /* MatchArrayString */

/***************************************************************************
 *	 Using the chain that was spliced on, we can reconize the splice    
 *	 and break it.                     
 ****************************************************************************/
void FlatSpliceOff (SeqEntryPtr the_set, ValNodePtr desc)
{
      BioseqSetPtr bss;
      BioseqPtr bs;
      ValNodePtr PNTR desc_head=NULL;
      ValNodePtr PNTR desc_target=NULL;
      ValNodePtr scan;

      if (IS_Bioseq(the_set) ){
         bs = (BioseqPtr) the_set -> data.ptrvalue;
         desc_head = & (bs -> descr);
      }else{
         bss = (BioseqSetPtr) the_set -> data.ptrvalue;
         desc_head = & (bss -> descr);
      }  
      if (* desc_head){
         desc_target = desc_head;
         for (scan = * desc_head; scan; scan = scan -> next){
            if (scan == desc){
               * desc_target = NULL;
               break;
            }
            desc_target = & (scan -> next);
         }
      }   

}

void FlatSpliceOn (SeqEntryPtr the_set, ValNodePtr desc)
{
      BioseqSetPtr bss;
      BioseqPtr bs;

      if (IS_Bioseq(the_set) ){
         bs = (BioseqPtr) the_set -> data.ptrvalue;
         bs -> descr = tie_next(bs -> descr, desc);
      } else {
         bss = (BioseqSetPtr) the_set -> data.ptrvalue;
         bss -> descr = tie_next(bss -> descr, desc);
      }  
}

/**************************************************************************
*	Get the ValNodePtr associated with a certain reference.
**************************************************************************/

static void GetPapRefPtr (Asn2ffJobPtr ajp, GBEntryPtr gbp, Int4 ext_index, Int4 pap_index, FFPrintArrayPtr pap)

{
	PubStructPtr psp=NULL;
	ValNodePtr vnp;
	Int4 i;
	DescrStructPtr dsp;
	
	for (vnp=gbp->Pub, i=0; vnp && i < ext_index; vnp=vnp->next, i++);
	if (vnp == NULL) {
		return;
	}
	psp = vnp->data.ptrvalue;
	if (psp == NULL) {
		return;
	}
	if ((dsp = pap[pap_index].descr) == NULL) {
		dsp =  (DescrStructPtr) MemNew(sizeof(DescrStruct));
		pap[pap_index].descr = dsp;
	}
	dsp->entityID = psp->entityID;
	dsp->itemID = psp->itemID;
	dsp->itemtype = psp->itemtype;
	
	return;
}

/**************************************************************************
*	Get the Comment structure associated with a certain comment block.
**************************************************************************/

static void GetPapCommPtr (Asn2ffJobPtr ajp, GBEntryPtr gbp, Int4 ext_index, Int4 pap_index, FFPrintArrayPtr pap)

{
	ComStructPtr s=NULL;
	Int4 i;
	DescrStructPtr dsp;
	
	for (s=gbp->comm, i=0; s && i < ext_index; s=s->next, i++);
	if (s == NULL) {
		return;
	}
	if ((dsp = pap[pap_index].descr) == NULL) {
		dsp =  (DescrStructPtr) MemNew(sizeof(DescrStruct));
		pap[pap_index].descr = dsp;
	}
	dsp->entityID = s->entityID;
	dsp->itemID = s->itemID;
	dsp->itemtype = s->itemtype;
	
	return;
}

/**************************************************************************
*	Find the SeqFeatPtr that is associated with this entry in the 
*	FFPrintArrayPtr.
*************************************************************************/

static void GetPapSeqFeatPtr (GBEntryPtr gbp, Int4 ext_index, Int4 pap_index, FFPrintArrayPtr pap)

{
	Int2 feat_index, index, listsize;
	OrganizeFeatPtr ofp;
	DescrStructPtr dsp;
	
	if (gbp == NULL || gbp->feat == NULL) {
		return;
	}
	ofp = gbp->feat;
	listsize=ofp->sfpListsize;
	index = (Int2) ext_index;

	feat_index = index - listsize;
	if (feat_index < 0) {
		if ((dsp = pap[pap_index].descr) == NULL) {
			dsp =  (DescrStructPtr) MemNew(sizeof(DescrStruct));
			pap[pap_index].descr = dsp;
		}
		dsp->entityID = ofp->List[index].entityID;
		dsp->itemID = ofp->List[index].itemID;
		dsp->itemtype = ofp->List[index].itemtype;
	}
	return;
}

NLM_EXTERN Boolean asn2ff_print (Asn2ffJobPtr ajp)
{
	AsnIoPtr			aip;
	CharPtr          string;
	FFPrintArrayPtr  pap = NULL;
	Int4             index, pap_size;
	Boolean 	 result = FALSE, hold = TRUE;

	if ((ajp->sep == NULL && ajp->slp == NULL) || ajp->fp == NULL)
		return FALSE;
	if (ajp->no_hold)
		hold = FALSE;
	if (hold)
		ObjMgrSetHold();   /* hold any autoloaded records in memory */

    pap_size = asn2ff_setup (ajp, &pap);
    if (ajp->ssp && ajp->format == EMBL_FMT) {
		aip = AsnIoNew(ASNIO_TEXT_OUT, ajp->fp, NULL, NULL, NULL);
		SubmitBlockAsnWrite(ajp->ssp->sub, aip, NULL);
		AsnIoFlush(aip);
		AsnIoReset(aip);
    }
    if (pap_size > 0) {
		head_www(ajp->fp, ajp->sep);
		asn2ff_set_output (NULL, "\n");
		for (index = 0; index < pap_size; index++) {
			string = FFPrint (pap, index, pap_size);
			if (string != NULL && *string != '\0') {
				ff_print_string (ajp->fp, string, "\n");
				string = MemFree (string);
			} else if (ajp->null_str) {
				ErrPostStr(SEV_WARNING, ERR_PRINT_NullString, 
				"CAUTION: NULL String returned\n");
			}
			if (pap[index].descr) {
				pap[index].descr = MemFree(pap[index].descr);
			}
		}
		tail_www(ajp->fp);
		result = TRUE;
		MemFree(pap);
	}
	free_buff();
	asn2ff_cleanup (ajp); 
	if (hold)
		ObjMgrClearHold();
	if (ajp->free_cache)
		ObjMgrFreeCache(0);
	
	return result;
}

Asn2ffJobPtr Asn2ffJobCreate(SeqEntryPtr sep, SeqSubmitPtr ssp, SeqLocPtr slp, FILE *fp, Uint1 format, Uint1 mode, StdPrintOptionsPtr	Spop)
{
	Asn2ffJobPtr	ajp;
	Uint2			entityID, itemID=0;
	
	ajp = (Asn2ffJobPtr) MemNew(sizeof(Asn2ffJob));
	ajp->show_seq = TRUE;
	ajp->show_gi = TRUE;
	ajp->error_msgs = TRUE;
	ajp->null_str = FALSE;
	ajp->non_strict = TRUE;
	ajp->format = format;
	ajp->mode = mode;
	ajp->show_gene = TRUE;
	ajp->gb_style = TRUE;
	ajp->fp = fp;
	ajp->Spop = Spop;
	if (ssp != NULL) {
		if ((entityID = ObjMgrGetEntityIDForPointer(ssp)) == 0) {
			ErrPostStr(SEV_WARNING, 0, 0, "Couldn't get entityID");
			MemFree(ajp);
			return NULL;
		}
		ajp->ssp = ssp;
		ajp->sep = (SeqEntryPtr) ssp->data;
	} else if (slp != NULL) {
		if ((entityID = BioseqFindEntity(SeqLocId(slp), &itemID)) == 0) {
			ErrPostStr(SEV_WARNING, 0, 0, "Couldn't get entityID");
			MemFree(ajp);
			return NULL;
		}
		ajp->slp = slp;
    	ajp->sep = NULL;
	} else {
		if ((entityID = ObjMgrGetEntityIDForChoice(sep)) == 0) {
			ErrPostStr(SEV_WARNING, 0, 0, "Couldn't get entityID");
			MemFree(ajp);
			return NULL;
		}
		ajp->sep = sep;
	} 
	ajp->entityID = entityID;
	
	return ajp;
}

typedef struct _link_str {
    CharPtr line;
    struct _link_str PNTR next;
} LinkStr, PNTR LinkStrPtr;

/**********************************************************/
static LinkStrPtr asn2ff_print_to_mem(Asn2ffJobPtr ajp, LinkStrPtr lsp)
{
    AsnIoPtr        aip;
    CharPtr         string;
    FFPrintArrayPtr pap = NULL;
    Int4            index, pap_size;
    Boolean         hold = TRUE;

    if(ajp->sep == NULL && ajp->slp == NULL)
        return(lsp);

    if(ajp->no_hold)
        hold = FALSE;
    else
        ObjMgrSetHold();        /* hold any autoloaded records in memory */

    pap_size = asn2ff_setup(ajp, &pap);
    if(ajp->ssp != NULL && ajp->format == EMBL_FMT && ajp->fp != NULL)
    {
        aip = AsnIoNew(ASNIO_TEXT_OUT, ajp->fp, NULL, NULL, NULL);
        SubmitBlockAsnWrite(ajp->ssp->sub, aip, NULL);
        AsnIoFlush(aip);
        AsnIoReset(aip);
    }

    if(pap_size > 0)
    {
        for(index = 0; index < pap_size; index++)
        {
            asn2ff_set_output(NULL, "\n");
            string = FFPrint(pap, index, pap_size);
            if(string != NULL && *string != '\0')
            {
                lsp->next = (LinkStrPtr) MemNew(sizeof(LinkStr));
                lsp = lsp->next;
                lsp->next = NULL;
                lsp->line = string;
                string = NULL;
            }
            else if(ajp->null_str != FALSE)
            {
                ErrPostStr(SEV_WARNING, ERR_PRINT_NullString,
                           "CAUTION: NULL String returned\n");
            }
            if(pap[index].descr != NULL)
            {
                pap[index].descr = MemFree(pap[index].descr);
            }
        }
        MemFree(pap);
    }
    free_buff();
    asn2ff_cleanup(ajp);
    if(hold != FALSE)
        ObjMgrClearHold();
    if(ajp->free_cache)
        ObjMgrFreeCache(0);

    return(lsp);
}

/**********************************************************/
static LinkStrPtr SeqEntryToLinkStr(Asn2ffJobPtr ajp, SeqEntryPtr sep,
                                    LinkStrPtr lsp, Uint1 format, Uint1 mode)
{
    StdPrintOptionsPtr Spop = NULL;
    BioseqSetPtr       bssp;

    if(sep == NULL)
        return(lsp);

    if(format == GENPEPT_FMT)
    {
        if(AllObjLoad() && SubmitAsnLoad() && SeqCodeSetLoad())
        {
            ErrShow();
        }
        if(Template_load == FALSE)
        {
            PrintTemplateSetLoad("asn2ff.prt");
            Template_load = TRUE;
        }
        Spop = StdPrintOptionsNew(NULL);
        if(Spop != NULL)
        {
            Spop->newline = "~";
            Spop->indent = "";
        }
        else
        {
            ErrPostStr(SEV_FATAL, 0, 0, "StdPrintOptionsNew failed");
            return(lsp);
        }
    }
    if(IS_Bioseq_set(sep) != 0)
    {
        bssp = (BioseqSetPtr) sep->data.ptrvalue;
        if(bssp != NULL && (bssp->_class == 7 || bssp->_class == 13 ||
           bssp->_class == 14 || bssp->_class == 15))
        {
            for(sep = bssp->seq_set; sep != NULL; sep = sep->next)
            {
                lsp = SeqEntryToLinkStr(ajp, sep, lsp, format, mode);
            }
            return(lsp);
        }
    }

    if(ajp == NULL)
        ajp = Asn2ffJobCreate(sep, NULL, NULL, NULL, format, mode, Spop);
    else
        ajp->sep = sep;

    if(ajp == NULL)
        return(lsp);

    lsp = asn2ff_print_to_mem(ajp, lsp);

    StdPrintOptionsFree(ajp->Spop);

    return(lsp);
}

/**********************************************************/
NLM_EXTERN CharPtr PNTR SeqEntryToStrArray(SeqEntryPtr sep, Uint1 format,
                                           Uint1 mode)
{
    LinkStrPtr   lsp;
    LinkStrPtr   tlsp;
    CharPtr PNTR res;
    CharPtr PNTR tres;
    Int4         num;

    lsp = (LinkStrPtr) MemNew(sizeof(LinkStr));
    lsp->next = NULL;
    lsp->line = NULL;
    SeqEntryToLinkStr(NULL, sep, lsp, format, mode);
    tlsp = lsp;
    lsp = lsp->next;
    MemFree(tlsp);

    for(tlsp = lsp, num = 1; tlsp != NULL; tlsp = tlsp->next, num++)
        continue;

    if(num == 1)
        return(NULL);

    res = (CharPtr PNTR) MemNew(sizeof(CharPtr) * num);
    for(tres = res; lsp != NULL; lsp = tlsp, tres++)
    {
        tlsp = lsp->next;
        *tres = lsp->line;
        MemFree(lsp);
    }
    *tres = NULL;
    return(res);
}

/***********************************************************************
*
*	SeqEntryToFlat is a stand-alone function that takes a SeqEntryPtr
*	and writes a flat file to a disk file.  If the formatting is
*	successful, TRUE is returned; otherwise FALSE is returned.
*
Choices for the Uint1's format and mode are defined in asn2ff.h.

For format they are:

GENBANK_FMT 	standard GenBank flat file for nucleotides
EMBL_FMT	standard EMBL flat file  for nucleotides
GENPEPT_FMT 	standard GenBank flat file for proteins
PSEUDOEMBL_FMT  a flavor of the EMBL flat file used by the "Authorin" program

The modes are:

RELEASE_MODE	this mode assures that all the requirements (e.g., identifiers
		features, references as described in the GenBank release notes
		and the feature table) are met.
		are met 
DUMP_MODE 	dump out the ASN.1 to a flat file
SEQUIN_MODE 	mode used by sequin
CHROMO_MODE 	mode used by Chromoscope
DIRSUB_MODE 	mode used by NCBI indexers during the "dirsub" process.
DIRSUB_DEBUG_MODE 	mode used by NCBI indexers during the "dirsub" process.
REVISE_MODE 	mode used by the "revise" program at NCBI (for in-house
		editing of entries).
*
**************************************************************************/

NLM_EXTERN Boolean SeqEntryToFlat (SeqEntryPtr sep, FILE *fp, Uint1 format, Uint1 mode)

{
	Boolean				rsult = FALSE;
	Asn2ffJobPtr		ajp;
	StdPrintOptionsPtr	Spop = NULL;
	BioseqSetPtr		bssp;

	if (sep == NULL) {
		return FALSE;
	}
	if (format == GENPEPT_FMT) {
		if (AllObjLoad () && SubmitAsnLoad () && SeqCodeSetLoad ()) {
			ErrShow();
		}
		if (!Template_load) {
			PrintTemplateSetLoad ("asn2ff.prt");
			Template_load = TRUE;
		}
		Spop = StdPrintOptionsNew(NULL);
		if (Spop) {
			Spop->newline = "~";
			Spop->indent = "";
		} else {
			ErrPostStr(SEV_FATAL,0,0, "StdPrintOptionsNew failed");;
			return FALSE;
		}
	}
	if (IS_Bioseq_set (sep)) {
    	bssp = (BioseqSetPtr) sep->data.ptrvalue;
    	if (bssp != NULL && (bssp->_class == 7 || bssp->_class == 13 ||
                         bssp->_class == 14 || bssp->_class == 15)) {
      		for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        		rsult = SeqEntryToFlat (sep, fp, format, mode);
			}
			return rsult;
		}
	}

	ajp = Asn2ffJobCreate(sep, NULL, NULL, fp, format, mode, Spop);

	if (ajp == NULL) {
		return FALSE;
	}
	rsult = asn2ff_print(ajp);
	StdPrintOptionsFree(ajp->Spop);
	MemFree(ajp);
	
	return rsult;
}

NLM_EXTERN Boolean SeqEntryToFlatAjp (Asn2ffJobPtr ajp, SeqEntryPtr sep, FILE *fp, Uint1 format, Uint1 mode)

{
	Boolean				rsult = FALSE;
	StdPrintOptionsPtr	Spop = NULL;
	BioseqSetPtr		bssp;
	
	if (sep == NULL) {
		return FALSE;
	}
	if (format == GENPEPT_FMT) {
		if (AllObjLoad () && SubmitAsnLoad () && SeqCodeSetLoad ()) {
			ErrShow();
		}
		if (!Template_load) {
			PrintTemplateSetLoad ("asn2ff.prt");
			Template_load = TRUE;
		}
		Spop = StdPrintOptionsNew(NULL);
		if (Spop) {
			Spop->newline = "~";
			Spop->indent = "";
		} else {
			ErrPostStr(SEV_FATAL,0,0, "StdPrintOptionsNew failed");;
			return FALSE;
		}
	}
	if (IS_Bioseq_set (sep)) {
    	bssp = (BioseqSetPtr) sep->data.ptrvalue;
    	if (bssp != NULL && (bssp->_class == 7 || bssp->_class == 13 ||
                         bssp->_class == 14 || bssp->_class == 15)) {
      		for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        		rsult = SeqEntryToFlatAjp (ajp, sep, fp, format, mode);
			}
			return rsult;
		}
	}
	if (sep == NULL) {
		StdPrintOptionsFree(ajp->Spop);
		return rsult;
	}
	if (ajp == NULL) {
		if ((ajp = Asn2ffJobCreate(sep, NULL, NULL, fp, format, mode, Spop))
																	== NULL) {
			return FALSE;
		}
	} else {
		if ((ajp->entityID=ObjMgrGetEntityIDForPointer(sep)) == 0) {
			ErrPostStr(SEV_WARNING, 0, 0, "Couldn't get entityID");
			return rsult;
		}
		ajp->sep = sep;
	}
	rsult = asn2ff_print(ajp);
	StdPrintOptionsFree(ajp->Spop);
	
	return rsult;
}

/**************************************************************************
 *	Prints out flat file in GenBank format WITHOUT Sequence
 **************************************************************************/
NLM_EXTERN Boolean SeqEntryToGBFlatNoSeq (SeqEntryPtr sep, 
                                          FILE *fp, Uint1 mode)
{
    Boolean	  rsult;
    Asn2ffJobPtr  ajp;
    BioseqSetPtr  bssp;
    
    if (sep == NULL) {
        return FALSE;
    }

    if (IS_Bioseq_set (sep)) {
    	bssp = (BioseqSetPtr) sep->data.ptrvalue;
    	if (bssp != NULL && (bssp->_class == 7 || bssp->_class == 13 ||
                             bssp->_class == 14 || bssp->_class == 15)) {
            for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
                rsult = SeqEntryToFlat (sep, fp, GENBANK_FMT, mode);
            }
            return rsult;
        }
    }
    
    ajp = Asn2ffJobCreate(sep, NULL, NULL, fp, GENBANK_FMT, mode, NULL);
    
    if (ajp == NULL)
        return FALSE;
    
    ajp->show_seq = FALSE; /* This is the point */
    
    rsult = asn2ff_print(ajp);
    MemFree(ajp);
	
    return rsult;
}

/***********************************************************************
*
*	SeqEntryToFlatEx is a stand-alone function works as SeqEntryToFlat
*	takes SeqIdPtr and various types of the output
*
*	successful, TRUE is returned; otherwise FALSE is returned.
*
*	Choices for the Uint1's type are defined in asn2ff.h.
*	FF_REGULAR 			0
*	FF_TOP_COMPLETE		1
*	FF_TOP_CONTIG		2
*
**************************************************************************/

NLM_EXTERN Boolean SeqEntryToFlatEx (SeqEntryPtr sep, FILE *fp, Uint1 format, Uint1 mode, SeqIdPtr seqid, Uint1 type)
{
	Boolean				rsult, repr = FALSE;
	Asn2ffJobPtr		ajp;
	StdPrintOptionsPtr	Spop = NULL;
	BioseqPtr 			bsp;
	BioseqSetPtr		bssp;

	rsult = FALSE;
	if (format == GENPEPT_FMT) {
		if (AllObjLoad () && SubmitAsnLoad () && SeqCodeSetLoad ()) {
			ErrShow();
		}
		if (!Template_load) {
			PrintTemplateSetLoad ("asn2ff.prt");
			Template_load = TRUE;
		}
		Spop = StdPrintOptionsNew(NULL);
		if (Spop) {
			Spop->newline = "~";
			Spop->indent = "";
		} else {
			Message (MSG_FATAL, "StdPrintOptionsNew failed");
			return rsult;
		}
	}
	if (IS_Bioseq_set (sep)) {
    	bssp = (BioseqSetPtr) sep->data.ptrvalue;
    	if (bssp != NULL && (bssp->_class == 7 || bssp->_class == 13 ||
                         bssp->_class == 14 || bssp->_class == 15)) {
      		for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        		rsult = SeqEntryToFlatEx (sep, fp, format, mode, seqid, type);
			}
			return rsult;
		}
	}
	ajp = Asn2ffJobCreate(sep, NULL, NULL, fp, format, mode, Spop);
	if (ajp == NULL) {
		return FALSE;
	}
	if (seqid != NULL) {
		ajp->gb_style = FALSE;
		ajp->id_print = seqid;
		bsp = BioseqFind(seqid);
    	if (ISA_na (bsp->mol) && bsp->repr == Seq_repr_seg) {
    		repr = TRUE;
    	}
		if (repr) {
			if (type == FF_REGULAR) {
				ajp->gb_style = TRUE;
				ajp->id_print = NULL;
			}
			if (type == FF_TOP_COMPLETE) {
				ajp->gb_style = FALSE;
				ajp->only_one = TRUE;
				ajp->ignore_top = TRUE;
			}
			if (type == FF_TOP_CONTIG) {
				ajp->gb_style = FALSE;
				ajp->only_one = TRUE;
				ajp->ignore_top = TRUE;
				ajp->genome_view = TRUE;
			}
			ajp->sep = sep;
		} else {
		ajp->sep = SeqMgrGetSeqEntryForData((Pointer)bsp);
		}
	} else {
		ajp->sep = sep;
	}
	if ((ajp->entityID = ObjMgrGetEntityIDForChoice(ajp->sep)) == 0) {
		ErrPostStr(SEV_WARNING, 0, 0, "Couldn't get entityID");
	}
	rsult = asn2ff_print(ajp);

	StdPrintOptionsFree(ajp->Spop);
	MemFree(ajp);
	
	return rsult;
}

/**************************************************************************
*	Prints out short flat file report in GenBank format
**************************************************************************/
NLM_EXTERN Boolean SeqEntryToPartRpt (SeqEntryPtr sep, FILE *fp)
{
	Boolean			rsult;
	Asn2ffJobPtr	ajp;

	ajp = Asn2ffJobCreate(sep, NULL, NULL, fp, GENBANK_FMT, PARTIAL_MODE, NULL);
	if (ajp == NULL) {
		return FALSE;
	}
	rsult = asn2ff_print(ajp);
	MemFree(ajp);
	
	return rsult;
}

NLM_EXTERN Boolean SeqSubmitToFlat (SeqSubmitPtr ssp, FILE *fp, Uint1 mode, Boolean show_gi, Uint1 format, Boolean show_gene)
{
	Boolean			rsult = FALSE;
	Asn2ffJobPtr	ajp;

	if (ssp == NULL) {
		return rsult;
	}
	if (ssp->datatype != 1) {
		return rsult;
	}
	ajp = Asn2ffJobCreate(NULL, ssp, NULL, fp, format, mode, NULL);
	if (ajp == NULL) {
		return FALSE;
	}
	rsult = asn2ff_print(ajp);
	MemFree(ajp);
	return rsult;
}

NLM_EXTERN Boolean SeqGenomeToFlat (SeqEntryPtr sep, FILE *fp, Uint1 format, Uint1 mode)

{
	Boolean          rsult;
	Asn2ffJobPtr		ajp;
	StdPrintOptionsPtr Spop = NULL;

	rsult = FALSE;
	ajp = Asn2ffJobCreate(sep, NULL, NULL, fp, format, mode, Spop);
	if (ajp == NULL) {
		return FALSE;
	}
    ajp->only_one = TRUE;
	ajp->ignore_top = TRUE;
	ajp->genome_view = TRUE;
	
	rsult = asn2ff_print(ajp);
	MemFree(ajp);
	
	return rsult;
}

NLM_EXTERN Boolean SeqGenomeToFlatEx (SeqEntryPtr sep, FILE *fp, Uint1 format, Uint1 mode, Boolean map_view)

{
	Boolean				rsult;
	Asn2ffJobPtr		ajp;
	StdPrintOptionsPtr Spop = NULL;

	rsult = FALSE;
	ajp = Asn2ffJobCreate(sep, NULL, NULL, fp, format, mode, Spop);
	if (ajp == NULL) {
		return FALSE;
	}
    ajp->only_one = TRUE;
	ajp->ignore_top = TRUE;
	ajp->genome_view = TRUE;
	ajp->map_view = map_view;

	rsult = asn2ff_print(ajp);
	MemFree(ajp);
	return rsult;
}

NLM_EXTERN Boolean SeqLocToFlat (SeqLocPtr slp, FILE *fp, Uint1 format, Uint1 mode)
{
	Boolean          rsult;
	Asn2ffJobPtr		ajp;
	StdPrintOptionsPtr Spop = NULL;
	Uint2            itemID=0;

	rsult = FALSE;
	if (format == GENPEPT_FMT) {
		if (AllObjLoad () && SubmitAsnLoad () && SeqCodeSetLoad () &&
								PrintTemplateSetLoad ("asn2ff.prt")) {
			ErrShow();
		}
		Spop = StdPrintOptionsNew(NULL);
		if (Spop) {
			Spop->newline = "~";
			Spop->indent = "";
		} else {
			Message (MSG_FATAL, "StdPrintOptionsNew failed");
			return rsult;
		}
	}

	ajp = Asn2ffJobCreate(NULL, NULL, slp, fp, format, mode, Spop);
	if (ajp == NULL) {
		return FALSE;
	}
    ajp->only_one = TRUE;
	ajp->ignore_top = TRUE;
	ajp->id_print = SeqLocId(slp);
	
	rsult = asn2ff_print(ajp);
	return rsult;
}


/***************************************************************************
*
*	Setup the FFPrintArrayPtr to be used by "FFPrint", the number 
*	returned is the number of entries in the array.
***************************************************************************/

NLM_EXTERN Int4 asn2ff_setup (Asn2ffJobPtr ajp, FFPrintArrayPtr PNTR papp)

{
	Int4 pap_size = -1;
	Asn2ffWEPtr awp;
	SeqIdPtr sip;
	Uint2 itemID;
	GatherScope gs;
	Uint1 focus;

  	MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
	MemSet ((Pointer) (gs.ignore), (int)(TRUE), (size_t) (OBJ_MAX * sizeof(Boolean)));
	gs.ignore[OBJ_BIOSEQ] = FALSE;
	gs.ignore[OBJ_BIOSEQSET] = FALSE;

	if (ajp->format == EMBLPEPT_FMT) /* Turn off Validators for EMBLPEPT */
		ajp->mode = DUMP_MODE;

	if (ajp->format == EMBLPEPT_FMT || ajp->format == GENPEPT_FMT)
		ajp->gb_style = FALSE;
	set_flags(ajp);

	flat2asn_install_accession_user_string("SET-UP");
	flat2asn_install_locus_user_string("SET-UP");
	

	ajp->sfp_out = 	MakeSyntheticSeqFeat();
	awp = (Asn2ffWEPtr) MemNew(sizeof(Asn2ffWE));
	awp->seg = NULL;
	awp->parts = NULL;
	ajp->hup = FALSE;
	if (ajp->ssp && ajp->ssp->sub) {
		ajp->hup = ajp->ssp->sub->hup;
	}
	ajp->asn2ffwep = awp;
	if (ajp->entityID == 0) {
	  if (ajp->sep != NULL) {
	    ajp->entityID = ObjMgrGetEntityIDForChoice (ajp->sep);
	  } else if (ajp->ssp != NULL) {
	    ajp->entityID = ObjMgrGetEntityIDForPointer (ajp->ssp);
	  } else if (ajp->slp != NULL) {
	    sip = SeqLocId (ajp->slp);
	    if (sip != NULL) {
	      ajp->entityID = BioseqFindEntity (sip, &itemID);
	    }
	  }
	}
	if (ajp->entityID != 0) {
		focus = (Uint1)FocusSeqEntry(ajp->sep, &gs);
		MemSet ((Pointer) (gs.ignore), (int)(TRUE), (size_t) (OBJ_MAX * sizeof(Boolean)));
		gs.ignore[OBJ_BIOSEQ] = FALSE;
		gs.ignore[OBJ_BIOSEQSET] = FALSE;
			GatherEntity(ajp->entityID, (Pointer) ajp, SeqToAwp, &gs);
		if (focus == FOCUS_INITIALIZED) {
			SeqLocFree(gs.target);
		}
		awp = ajp->asn2ffwep;
		if (!ajp->only_one && awp->gbp == NULL) {
			if (awp) {
				ajp->asn2ffwep = MemFree(awp);
			}
			return 0;
		}
		if (awp->seg == NULL && awp->parts == NULL) {
			awp->total_seg = 0;
			if (awp->gbp) {
				awp->gbp->num_seg = 0;
			}
			if (SeqMgrFeaturesAreIndexed (ajp->entityID)) {
				ajp->useSeqMgrIndexes = TRUE;  /* initial use of new indexes */
			}
		}
		if (awp->gbp && awp->gbp->next == NULL) {
			awp->total_seg = 0;
			awp->gbp->num_seg = 0;
			if (SeqMgrFeaturesAreIndexed (ajp->entityID)) {
				ajp->useSeqMgrIndexes = TRUE;  /* initial use of new indexes */
			}
		}
		ajp->asn2ffwep = awp;
		GetGIs(ajp);
	}
	init_buff();
	ajp->pseudo = FALSE;
	if (ajp->format == SELECT_FMT) {    /* quick fix 07.17.95 change later */
		ajp->format = GENBANK_FMT;
		ajp->pseudo = TRUE;
	}
	if (ajp->format == PSEUDOEMBL_FMT) {
		ajp->pseudo = TRUE;
	}
	if (ajp->help) {
		pap_size = asn2hp_setup(ajp, papp);
		return pap_size;
	}
	if (ajp->format == GENBANK_FMT || ajp->format == SELECT_FMT) {
		if (ajp->mode == PARTIAL_MODE) {
			pap_size = asn2pr_setup(ajp, papp);
		} else {
			pap_size = asn2gb_setup(ajp, papp);
		}
	} else if (ajp->format == EMBL_FMT || ajp->format == PSEUDOEMBL_FMT) {
		pap_size = asn2embl_setup(ajp, papp);
	} else if (ajp->format == EMBLPEPT_FMT) {
			pap_size = asn2ep_setup(ajp,  papp);
	}else if (ajp->format == GENPEPT_FMT) {
		pap_size = asn2gp_setup(ajp, papp);
	}
	
	return pap_size;
}	/* asn2ff_setup */

/****************************************************************************
*void set_flags (Asn2ffJobPtr ajp)
*
*	set_flags to determine which tasks to perform.
*****************************************************************************/
void set_flags (Asn2ffJobPtr ajp)

{

/* The defines are:
ASN2FF_LOCAL_ID                 asn2ff_flags[0]
	If FALSE then entries with "local" id's are NOT formatted 
ASN2FF_LOOK_FOR_SEQ             asn2ff_flags[1]
	If TRUE BioseqFind is run in an attempt to "find" entries that
	have been loaded into memory and are referenced by an entry 
ASN2FF_VALIDATE_FEATURES        asn2ff_flags[2]
	If TRUE then validation is run on features.  If they are invalid
	they are dropped.
ASN2FF_IGNORE_PATENT_PUBS            asn2ff_flags[3]
	This flag only applies to patent pubs.  If FlatIgnoreThisPatentPub
	is true and this flag is TRUE, that pub is dropped.  ALL OTHER
	PUBS are validated all the time.
ASN2FF_DROP_SHORT_AA            asn2ff_flags[4]
	Drop amino acid sequences that are too short.  Only applies to 
	GenPept (i.e., protein) format  
ASN2FF_AVOID_LOCUS_COLL         asn2ff_flags[5]
	If TRUE Check for LOCUS collisions with Karl's algorithm
	Otherwise Use the LOCUS in the id field.
ASN2FF_DATE_ERROR_MSG		asn2ff_flags[6]
	If TRUE report a missing date.  SHould be FALSE for indexing
	work when no date for a record has been set.
ASN2FF_IUPACAA_ONLY		asn2ff_flags[7]
	Use only iupaca characters if TRUE.  Only iupacaa is the flat
	file standard.
ASN2FF_TRANSL_TABLE		asn2ff_flags[8]
	If TRUE print the transl_table qualifiers.  Set to FALSE until
	the database correctly reflects transl_tables.
ASN2FF_REPORT_LOCUS_COLL     	asn2ff_flags[9]
	If TRUE, report locus collisions via ErrPostEx
ASN2FF_SHOW_ALL_PUBS	        asn2ff_flags[10]
	if TRUE don't drop CitGen reference or replace CitGen->cit with
	"Unpublished"
ASN2FF_SHOW_ERROR_MSG	        asn2ff_flags[11]
ASN2FF_SHOW_GB_STYLE	        asn2ff_flags[12]
	show only features comleted on this bioseq or location - gb_style
	
*/

	asn2ff_flags[11] = ajp->error_msgs;
	asn2ff_flags[12] = ajp->gb_style;
	if (ajp->mode == RELEASE_MODE)
	{
		asn2ff_flags[0] = FALSE; 
		if (GetAppProperty ("InternalNcbiSequin") != NULL) {
			asn2ff_flags[0] = TRUE; 
		}
		asn2ff_flags[1] = FALSE; 
		asn2ff_flags[2] = TRUE; 
		asn2ff_flags[3] = TRUE; 
		asn2ff_flags[4] = TRUE; 
		asn2ff_flags[5] = TRUE; 
		asn2ff_flags[6] = TRUE; 
		asn2ff_flags[7] = TRUE; 
		asn2ff_flags[8] = TRUE; 
		asn2ff_flags[9] = FALSE; 
		asn2ff_flags[10] = FALSE; 
	}
	else if (ajp->mode == DIRSUB_MODE)
	{
		asn2ff_flags[0] = FALSE; 
		asn2ff_flags[1] = FALSE; 
		asn2ff_flags[2] = TRUE; 
		asn2ff_flags[3] = TRUE; 
		asn2ff_flags[4] = TRUE; 
		asn2ff_flags[5] = TRUE; 
		asn2ff_flags[6] = FALSE; 
		asn2ff_flags[7] = FALSE; 
		asn2ff_flags[8] = TRUE; 
		asn2ff_flags[9] = FALSE; 
		asn2ff_flags[10] = FALSE;
		ajp->show_gi = FALSE; 
	}
	else if (ajp->mode == DIRSUB_DEBUG_MODE)
	{
		asn2ff_flags[0] = FALSE; 
		asn2ff_flags[1] = FALSE; 
		asn2ff_flags[2] = FALSE; 
		asn2ff_flags[3] = TRUE; 
		asn2ff_flags[4] = TRUE; 
		asn2ff_flags[5] = TRUE; 
		asn2ff_flags[6] = FALSE; 
		asn2ff_flags[7] = FALSE; 
		asn2ff_flags[8] = TRUE; 
		asn2ff_flags[9] = FALSE; 
		asn2ff_flags[10] = FALSE;
		ajp->show_gi = FALSE; 
	}
	else if (ajp->mode == REVISE_MODE)
	{
		asn2ff_flags[0] = TRUE; 
		asn2ff_flags[1] = FALSE; 
		asn2ff_flags[2] = FALSE; 
		asn2ff_flags[3] = FALSE; 
		asn2ff_flags[4] = FALSE; 
		asn2ff_flags[5] = FALSE; 
		asn2ff_flags[6] = TRUE; 
		asn2ff_flags[7] = FALSE; 
		asn2ff_flags[8] = TRUE; 
		asn2ff_flags[9] = TRUE; 
		asn2ff_flags[10] = FALSE; 
	}
	else if (ajp->mode == DUMP_MODE)
	{
		asn2ff_flags[0] = TRUE; 
		asn2ff_flags[1] = FALSE; 
		asn2ff_flags[2] = FALSE; 
		asn2ff_flags[3] = FALSE; 
		asn2ff_flags[4] = FALSE; 
		asn2ff_flags[5] = FALSE; 
		asn2ff_flags[6] = TRUE; 
		asn2ff_flags[7] = FALSE; 
		asn2ff_flags[8] = TRUE; 
		asn2ff_flags[9] = FALSE; 
		asn2ff_flags[10] = TRUE; 
	}
	else if (ajp->mode == SEQUIN_MODE)
	{
		asn2ff_flags[0] = TRUE; 
		asn2ff_flags[1] = FALSE; 
		asn2ff_flags[2] = FALSE; 
		asn2ff_flags[3] = TRUE; 
		asn2ff_flags[4] = TRUE; 
		asn2ff_flags[5] = TRUE; 
		asn2ff_flags[6] = FALSE; 
		asn2ff_flags[7] = FALSE; 
		asn2ff_flags[8] = TRUE; 
		asn2ff_flags[9] = FALSE; 
		asn2ff_flags[10] = FALSE; 
	}
	else if (ajp->mode == CHROMO_MODE)
	{
		asn2ff_flags[0] = TRUE; 
		asn2ff_flags[1] = TRUE; 
		asn2ff_flags[2] = FALSE; 
		asn2ff_flags[3] = TRUE; 
		asn2ff_flags[4] = FALSE; 
		asn2ff_flags[5] = FALSE; 
		asn2ff_flags[6] = FALSE; 
		asn2ff_flags[7] = FALSE; 
		asn2ff_flags[8] = TRUE; 
		asn2ff_flags[9] = FALSE; 
		asn2ff_flags[10] = FALSE; 
	}
}

static Boolean check_whole(SeqFeatPtr f, Int4 len)
{
	Boolean whole = FALSE;
	SeqLocPtr		slp;
	SeqIntPtr		sip;
	
		slp = f->location;
		if (slp->choice == SEQLOC_WHOLE) {
			whole = TRUE;
		} else if (slp->choice == SEQLOC_INT) {
			sip = slp->data.ptrvalue;
			if (sip->from == 0 && sip->to == len-1) {
				whole = TRUE;
			}
		}
	return whole;
}
Boolean get_pubs (GatherContextPtr gcp)
{
	ValNodePtr	tmp, vnp, v;
	PubdescPtr 	pdp;
	ValNodePtr	PNTR vnpp;
	BioseqPtr	bsp;
	SeqLocPtr	slp;
	SeqFeatPtr	sfp;
	ImpFeatPtr	ifp;
	SubmitBlockPtr sbp;
	CitSubPtr 	the_cit;
	
	vnpp = gcp->userdata;
	vnp = *vnpp;
	switch (gcp->thistype)
	{
		case OBJ_SEQDESC:
			tmp = (ValNodePtr) (gcp->thisitem);
			if (gcp->parenttype == OBJ_BIOSEQ) {
				bsp = (BioseqPtr) (gcp->parentitem);
			} else {
				bsp = NULL;
			}
			if (tmp->choice == Seq_descr_pub) {
				vnp = StorePub(bsp, vnp, tmp, NULL, 1, gcp->entityID,
					gcp->itemID, gcp->thistype);
			}
			break;
		case OBJ_SEQFEAT:
			sfp = (SeqFeatPtr) (gcp->thisitem);
			if (sfp->data.choice == SEQFEAT_PUB) {
				slp = sfp->location;
				bsp = BioseqFindCore(SeqLocId(slp));
				if (check_whole(sfp, bsp->length)) {
					tmp = ValNodeNew(NULL);
					tmp->choice = Seq_descr_pub;
					tmp->data.ptrvalue = (PubdescPtr) sfp->data.value.ptrvalue;
					vnp = StorePub(bsp, vnp, tmp, NULL, 1, gcp->entityID,
					gcp->itemID, gcp->thistype);
				} else {
					vnp = StorePub(NULL, vnp, NULL, sfp, 2, gcp->entityID,
					gcp->itemID, gcp->thistype);
				}
			} 
			if (sfp->data.choice == SEQFEAT_IMP) {
				ifp = sfp->data.value.ptrvalue;
				if (StringCmp(ifp->key, "Site-ref") == 0) {
					if (sfp->cit != NULL) {
						vnp = StorePub(NULL, vnp, NULL, sfp, 3, gcp->entityID,
						gcp->itemID, gcp->thistype);
					}
				}
			}
			break;
		case OBJ_SUBMIT_BLOCK:
			sbp = (SubmitBlockPtr) (gcp->thisitem);
			the_cit = AsnIoMemCopy(sbp->cit, (AsnReadFunc) CitSubAsnRead,
			(AsnWriteFunc) CitSubAsnWrite);
			v = ValNodeNew(NULL);
			v->choice = PUB_Sub;
			v->data.ptrvalue = the_cit;
			pdp = PubdescNew();
			pdp->pub = v;
			tmp = ValNodeNew(NULL);
			tmp->choice = Seq_descr_pub;
			tmp->data.ptrvalue = pdp;
			vnp = StorePub(NULL, vnp, tmp, NULL, 1, gcp->entityID,
						gcp->itemID, gcp->thistype);
			PubdescFree(pdp);
			break;
	/*	case OBJ_SEQSUB_CIT:
			csp = (CitSubPtr) (gcp->thisitem);
			the_cit = AsnIoMemCopy(csp, (AsnReadFunc) CitSubAsnRead,
			(AsnWriteFunc) CitSubAsnWrite);
			v = ValNodeNew(NULL);
			v->choice = PUB_Sub;
			v->data.ptrvalue = the_cit;
			pdp = PubdescNew();
			pdp->pub = v;
			tmp = ValNodeNew(NULL);
			tmp->choice = Seq_descr_pub;
			tmp->data.ptrvalue = pdp;
			vnp = StorePub(NULL, vnp, tmp, NULL, 1, gcp->entityID,
						gcp->itemID, gcp->thistype);
			MemFree(csp);
			break;
	*/
		case OBJ_SEQFEAT_CIT:
/***** not used now ********/
			tmp = (ValNodePtr) (gcp->thisitem); /* PubSet->data.ptrvalue */
			break;
		default:
			break;
	}
	*vnpp = vnp;
	return TRUE;
}

static Int2 GetPubNum(GBEntryPtr gbp)
{
	ValNodePtr v;
	Int4 i;

	for (v = gbp->Pub, i=0; v != NULL; v= v->next, i++);
	
	return (Int2)i;
}
static void CheckSourceFeat(Asn2ffJobPtr ajp, GBEntryPtr gbp)
{
	OrgRefPtr orp;
	BioSourcePtr biosp;
	ValNodePtr vnp;
	DescrStructPtr ds;
	
	if (gbp && gbp->feat) {
		if (gbp->feat->sfpSourcesize != 0) 
			return;
	}
	ds = gbp->source_info;
	if ((vnp=GatherDescrByChoice(ajp, gbp, Seq_descr_source)) != NULL) {
		biosp = vnp->data.ptrvalue;
		orp = (OrgRefPtr) biosp->org;	
		if (orp) {
			if (ds == NULL) {
				ds = (DescrStructPtr) MemNew(sizeof(DescrStruct));
				gbp->source_info = ds;
			}
			ds->vnp = vnp;
			ds->entityID = gbp->descr->entityID;
			ds->itemID = gbp->descr->itemID;
			ds->itemtype = gbp->descr->itemtype;
			gbp->descr = MemFree(gbp->descr);
			return;
		}
	}
	if (gbp && gbp->descr) {
			MemFree(gbp->descr);
	}
	if ((vnp=GatherDescrByChoice(ajp, gbp, Seq_descr_org)) != NULL) {
		orp = (OrgRefPtr) vnp->data.ptrvalue;
		if (orp) {
			if (ds == NULL) {
				ds = (DescrStructPtr) MemNew(sizeof(DescrStruct));
				gbp->source_info = ds;
			}
			ds->vnp = vnp;
			ds->entityID = gbp->descr->entityID;
			ds->itemID = gbp->descr->itemID;
			ds->itemtype = gbp->descr->itemtype;
			gbp->descr = MemFree(gbp->descr);
			return;
		}
	}
	if (gbp && gbp->descr) {
			gbp->descr = MemFree(gbp->descr);
	}
	return;
}

Int4 asn2hp_setup(Asn2ffJobPtr ajp, FFPrintArrayPtr PNTR papp)
{
	FFPrintArrayPtr pap;
	Int4 index, total, pub_num;
	GBEntryPtr gbp;

	GetLocusPartsAwp(ajp); 
	total=2;
	for (gbp=ajp->asn2ffwep->gbp; gbp; gbp = gbp->next) {
		gbp->descr = NULL;	
		if (GB_GetSeqDescrComms(ajp, gbp) > 0) {
			total += gbp->comm_num;
		}
		pub_num = (Int2)GetPubsAwp(ajp, gbp);
		total += pub_num; 
		GetGBDate(ajp, gbp);
	}
	*papp = (FFPrintArrayPtr) MemNew((size_t) total*sizeof(FFPrintArray));
	pap = *papp;
	/* pap_total = total; -- NO EFFECT */
	LoadPap(NULL, NULL, ajp, 0, (Uint1)0, (Uint1)0, 0, A2F_OTHER, NULL);
	for (gbp=ajp->asn2ffwep->gbp; gbp; gbp = gbp->next)
	{
		GetDefinitionLine(ajp, gbp);
		LoadPap(pap, PrintDefinitionLine, ajp, 0, (Uint1)0, (Uint1)0, 
			line_estimate[1], A2F_OTHER, gbp);
		LoadPap(pap, PrintGBOrganismLine, ajp, 0, (Uint1)0, (Uint1)0, 
			line_estimate[4], A2F_OTHER, gbp);
		pub_num = GetPubNum(gbp);
		for (index=0; index < pub_num; index++) {
			LoadPap(pap, 
				PrintPubsByNumber, ajp, index, (Uint1)0, (Uint1)0, 
				line_estimate[5], A2F_REFERENCE, gbp);
		}
		for (index=0; index < gbp->comm_num; index++) {
			if (index == 0) {
				LoadPap(pap, 
					PrintFirstComment, ajp, index, (Uint1)0, (Uint1)0, 
						line_estimate[5], A2F_COMMENT, gbp);
			} else {
				LoadPap(pap, 
					PrintCommentByNumber, ajp, index, (Uint1)0, (Uint1)0, 
						line_estimate[5], A2F_COMMENT, gbp);
			}
		}
	}

	return total;
}

static void PrintLastLine (Asn2ffJobPtr ajp, GBEntryPtr gbp)
{
	PrintTerminator ();
}

Int4 asn2gb_setup(Asn2ffJobPtr ajp, FFPrintArrayPtr PNTR papp)
{
	FFPrintArrayPtr pap;
	Int4 index, total, pub_num, seqblks_num;
	GBEntryPtr gbp;
	SeqIdPtr sip;

	GetLocusPartsAwp(ajp);
	if (!ajp->genome_view) {
		GetSeqFeat(ajp);
	}
	total=0;
	for (gbp=ajp->asn2ffwep->gbp; gbp; gbp = gbp->next) {
		if (gbp->bsp && ajp->id_print) {
			sip = SeqIdFindBest(gbp->bsp->id, SEQID_GI);
			if (SeqIdComp(sip, ajp->id_print) != SIC_YES) {
				continue;
			}
		}
		CheckSourceFeat(ajp, gbp);
		if (gbp->map == TRUE && ajp->show_seq == FALSE) {
			total += 7;   
		} else if (ajp->genome_view) {
			total += 6;   
		} else {
			total += 8;
		}
		if (gbp->gi != -1 || gbp->ni != NULL) {
			total++;
		}
		if (ajp->asn2ffwep->total_seg > 0) {
			total++;
		}	
		gbp->descr = NULL;	
		if (GB_GetSeqDescrComms(ajp, gbp) > 0) {
			total += gbp->comm_num;
		}
		if	(gbp->feat && gbp->feat->sfpCommsize > 0) {
			total++;
		}
		if (ajp->genome_view || ajp->map_view) {
			total += 2;				/* FEATURES and 'source' feature*/
			total ++;				/* last line '//' */
			if (gbp->map) {
				gbp->feat_num = GetMapFeats(ajp, gbp);
				total += gbp->feat_num;
			} else {
				total ++;
			}
		} else {
			if (gbp->feat) {
				total += 2;				/* FEATURES and 'source' feature*/
				gbp->feat_num = gbp->feat->sfpListsize;
				total += gbp->feat_num;
			}
			if (ajp->show_seq == TRUE) {
				seqblks_num = (Int2)GetNumOfSeqBlks(ajp, gbp);
				total += seqblks_num;
			}
		}
		pub_num = (Int2)GetPubsAwp(ajp, gbp);
		total += pub_num; 
		GetGBDate(ajp, gbp);
	}
	*papp = (FFPrintArrayPtr) MemNew((size_t) total*sizeof(FFPrintArray));
	pap = *papp;
	/* pap_total = total; -- NO EFFECT */
	LoadPap(NULL, NULL, ajp, 0, (Uint1)0, (Uint1)0, 0, A2F_OTHER, NULL);
	for (gbp=ajp->asn2ffwep->gbp; gbp; gbp = gbp->next)
	{
		if (gbp->bsp && ajp->id_print) {
			sip = SeqIdFindBest(gbp->bsp->id, SEQID_GI);
			if (SeqIdComp(sip, ajp->id_print) != SIC_YES) {
				continue;
			}
		}
		LoadPap(pap, PrintLocusLine, ajp, 0, (Uint1)0, (Uint1)0,
									line_estimate[0], A2F_OTHER, gbp);
		flat2asn_delete_locus_user_string();
		flat2asn_install_locus_user_string(gbp->locus);
		GetDefinitionLine(ajp, gbp);
		LoadPap(pap, PrintDefinitionLine, ajp, 0, (Uint1)0, (Uint1)0, 
			line_estimate[1], A2F_OTHER, gbp);
		if (gbp->descr) {
			gbp->descr = MemFree(gbp->descr);
			gbp->descr = NULL;
		}
		LoadPap(pap, PrintAccessLine, ajp, 0, (Uint1)0, (Uint1)0, 
			line_estimate[2], A2F_OTHER, gbp);
		flat2asn_delete_accession_user_string();
		flat2asn_install_accession_user_string(gbp->accession);
		if (gbp->gi != -1) {
			LoadPap(pap, PrintNCBI_GI, ajp, 0, (Uint1)0, (Uint1)0, 
				line_estimate[2], A2F_OTHER, gbp);
		} else if (gbp->ni != NULL) {
			LoadPap(pap, PrintNID, ajp, 0, (Uint1)0, (Uint1)0, 
				line_estimate[2], A2F_OTHER, gbp);
		}
		LoadPap(pap, PrintKeywordLine, ajp, 0, (Uint1)0, (Uint1)0, 
			line_estimate[3], A2F_OTHER, gbp);
		if (ajp->asn2ffwep->total_seg > 0)
			LoadPap(pap, PrintSegmentLine, ajp, 0, (Uint1)0, (Uint1)0, 
				line_estimate[0], A2F_OTHER, gbp);
		LoadPap(pap, PrintGBSourceLine, ajp, 0, (Uint1)0, (Uint1)0, 
			line_estimate[4], A2F_OTHER, gbp);
		LoadPap(pap, PrintGBOrganismLine, ajp, 0, (Uint1)0, (Uint1)0, 
			line_estimate[4], A2F_OTHER, gbp);
		pub_num = GetPubNum(gbp);
		for (index=0; index < pub_num; index++) {
			LoadPap(pap, 
				PrintPubsByNumber, ajp, index, (Uint1)0, (Uint1)0, 
				line_estimate[5], A2F_REFERENCE, gbp);
		}
		for (index=0; index < gbp->comm_num; index++) {
			if (index == 0) {
				LoadPap(pap, 
					PrintFirstComment, ajp, index, (Uint1)0, (Uint1)0, 
						line_estimate[5], A2F_COMMENT, gbp);
			} else {
				LoadPap(pap, 
					PrintCommentByNumber, ajp, index, (Uint1)0, (Uint1)0, 
						line_estimate[5], A2F_COMMENT, gbp);
			}
		}
		if (gbp->feat && gbp->feat->sfpCommsize > 0) {
			LoadPap(pap, GBDescrComFeat, ajp, 0, (Uint1)0, (Uint1)0, 
					line_estimate[2], A2F_OTHER, gbp);
		}
		if (ajp->genome_view && gbp->map == FALSE) {
			LoadPap(pap, PrintFeatHeader, ajp, 0, (Uint1)0, (Uint1)0, 
					line_estimate[6], A2F_OTHER, gbp);
			LoadPap(pap, PrintSourceFeat, ajp, 0, (Uint1)0, (Uint1)0, 
					line_estimate[8], A2F_SOURCE_FEATURE, gbp);
			LoadPap(pap, PrintGenome, ajp, 0, (Uint1)0, (Uint1)0, 
				line_estimate[0], A2F_OTHER, gbp);
		} else {
			if (gbp->feat) {
				LoadPap(pap, 	PrintFeatHeader, ajp, 0, (Uint1)0, (Uint1)0, 
					line_estimate[6], A2F_OTHER, gbp);
				LoadPap(pap, PrintSourceFeat, ajp, 0, (Uint1)0, (Uint1)0, 
					line_estimate[8], A2F_SOURCE_FEATURE, gbp);
				for (index=0; index < gbp->feat_num; index++) {
					LoadPap(pap, PrintNAFeatByNumber, ajp, index, 
						(Uint1)0, (Uint1)0, line_estimate[8], A2F_FEATURE, gbp);
				}
			}
			if (gbp->map == FALSE && ajp->show_seq == TRUE) {
				LoadPap(pap, PrintBaseCount, ajp, 0, (Uint1)0, (Uint1)0, 
					line_estimate[0], A2F_OTHER, gbp);
				LoadPap(pap, PrintOriginLine, ajp, 0, (Uint1)0, (Uint1)0, 
					line_estimate[0], A2F_OTHER, gbp);
				seqblks_num = GetNumOfSeqBlks(ajp, gbp);
				for (index=0; index < seqblks_num; index++) {
					if (seqblks_num == index+1) {
						LoadPap(pap, PrintSeqBlk, ajp, index, 
						(Uint1)1, (Uint1)0, line_estimate[9], A2F_OTHER, gbp);
					} else {
						LoadPap(pap, PrintSeqBlk, ajp, index, 
						(Uint1)0, (Uint1)0, line_estimate[9], A2F_OTHER, gbp);
					}
				}
			} else {
				LoadPap(pap, PrintLastLine, ajp, 0, (Uint1)0, (Uint1)0, 
					line_estimate[0], A2F_OTHER, gbp);
			}
		}
	}

	return total;
}

/*______________________________________________________________________
**
**	This code is not currently used.
**	I do not remove this piece of code, just comment it out.
**	-- Dmitri Lukyanov
*/
#if 0

static Int2 GetCDSNumber (OrganizeFeatPtr feat)

{
	SortStructPtr p;
	SeqFeatPtr sfp;
	Int2 i,j;

	for (p = feat->List, j=0, i=0; i < feat->sfpListsize; i++, p++) {
		sfp = p->sfp;
		if (sfp->data.choice != SEQFEAT_CDREGION) {
			continue;
		}
		j++;
	}
	return j;
}

static Int2 GetPubsRptNum (ValNodePtr pubs)
{
	PubStructPtr 	psp;
	ValNodePtr 		vnp;
	Int2 			i;

	for (vnp=pubs, i=0; vnp; vnp=vnp->next) {
		psp = vnp->data.ptrvalue;
		if (psp->choice == PUB_Sub) {
			continue;
		}
		i++;
	}
	return i;
}

#endif
/*______________________________________________________________________
*/

Int4 asn2pr_setup(Asn2ffJobPtr ajp, FFPrintArrayPtr PNTR papp)
{
	FFPrintArrayPtr pap;
	Int4 index, total, feat_num, pub_num;
	GBEntryPtr gbp;

	GetLocusPartsAwp(ajp);
	GetSeqFeat(ajp);
	total=0;
	for (gbp=ajp->asn2ffwep->gbp; gbp; gbp = gbp->next) {
		total += 2;
		if (gbp->feat) {
			total++;
			feat_num = gbp->feat->sfpListsize;
			total += feat_num;
		}
		pub_num = GetPubsAwp(ajp, gbp);
		total += pub_num; 
	}
	*papp = (FFPrintArrayPtr) MemNew((size_t) total*sizeof(FFPrintArray));
	pap = *papp;
	/* pap_total = total; -- NO EFFECT */
	LoadPap(NULL, NULL, ajp, 0, (Uint1)0, (Uint1)0, 0, A2F_OTHER, NULL);
	for (gbp=ajp->asn2ffwep->gbp; gbp; gbp = gbp->next) {
		LoadPap(pap, PrintAccessLine, ajp, 0, (Uint1)0, (Uint1)0, 
			line_estimate[2], A2F_OTHER, gbp);
		pub_num = (Int2)GetPubsAwp(ajp, gbp);
		for (index=0; index < pub_num; index++) {
			LoadPap(pap, 
				PrintPubsByNumber, ajp, index, (Uint1)0, (Uint1)0, 
				line_estimate[5], A2F_REFERENCE, gbp);
		}
		if (gbp->feat) {
			LoadPap(pap, 	PrintFeatHeader, ajp, 0, (Uint1)0, (Uint1)0, 
				line_estimate[6], A2F_OTHER, gbp);
			for (index=0; index < gbp->feat->sfpListsize; index++) {
				LoadPap(pap, PrintNAFeatByNumber, ajp, index, 
					(Uint1)0, (Uint1)0, line_estimate[8], A2F_FEATURE, gbp);
			}
		}
		LoadPap(pap, PrintLastLine, ajp, 0, (Uint1)0, (Uint1)0, 
				line_estimate[0], A2F_OTHER, gbp);
	}

	return total;
}

Int4 asn2embl_setup(Asn2ffJobPtr ajp, FFPrintArrayPtr PNTR papp)
{
	FFPrintArrayPtr pap;
	Int4 index, max, total, pub_num, seqblks_num;
	GBEntryPtr gbp;

	GetLocusPartsAwp(ajp);
	GetSeqFeat(ajp);

	total=0;
	for (gbp=ajp->asn2ffwep->gbp; gbp; gbp = gbp->next) {
		CheckSourceFeat(ajp, gbp);
		total += 7;
		gbp->xref_present = FALSE;
		if (CheckXrefLine(ajp, gbp) == TRUE) {
			total ++;
			gbp->xref_present = TRUE;
		}
		if (gbp->gi != -1 || gbp->ni != NULL) {
			total++;
		}
		gbp->descr = NULL;	
		if (GB_GetSeqDescrComms(ajp, gbp) > 0) {
			total += gbp->comm_num;
		}
		if	(gbp->feat && gbp->feat->sfpCommsize > 0) {
			total++;
		}
		if (gbp->feat) {
			total += 2;				/* FEATURES and 'source' feature*/
			total += gbp->feat->sfpListsize;
		}
		seqblks_num = GetNumOfSeqBlks(ajp, gbp);
		total += seqblks_num;
		pub_num = GetPubsAwp(ajp, gbp); 
		total += pub_num; 

		GetEMBLDate(ajp, gbp);
		GetEntryVersion(gbp);
	}
	if (ajp->ssp && ajp->hup)
		total --;
	*papp = (FFPrintArrayPtr) MemNew((size_t) total*sizeof(FFPrintArray));
	pap = *papp;

	LoadPap(NULL, NULL, ajp, 0, (Uint1)0, (Uint1)0, 0, A2F_OTHER, NULL);
	for (gbp=ajp->asn2ffwep->gbp; gbp; gbp = gbp->next) {
		LoadPap(pap,
			PrintLocusLine, ajp, 0,(Uint1)0,(Uint1)1,line_estimate[0], 
																A2F_OTHER, gbp);
		flat2asn_delete_locus_user_string();
		flat2asn_install_locus_user_string(gbp->locus);
		LoadPap(pap, 
			PrintAccessLine,ajp,0,(Uint1)0,(Uint1)1,line_estimate[2], 
																A2F_OTHER, gbp);
		flat2asn_delete_accession_user_string();
		flat2asn_install_accession_user_string(gbp->accession);
		if (gbp->gi != -1) {
			LoadPap(pap, PrintNCBI_GI, ajp, 0, (Uint1)0, (Uint1)0, 
				line_estimate[2], A2F_OTHER, gbp);
		} else if (gbp->ni != NULL) {
			LoadPap(pap, PrintNID, ajp, 0, (Uint1)0, (Uint1)0, 
				line_estimate[2], A2F_OTHER, gbp);
		}
		if (ajp->ssp == NULL || ajp->hup == FALSE) {
			LoadPap(pap, 
			PrintDateLines, ajp,0,(Uint1)0,(Uint1)1,line_estimate[10],
																A2F_OTHER, gbp);
		}
		GetDefinitionLine(ajp, gbp);
		LoadPap(pap,
			 PrintDefinitionLine,ajp,0,(Uint1)0,(Uint1)1,line_estimate[1],
			 								 					A2F_OTHER, gbp);
		LoadPap(pap, PrintKeywordLine,ajp,0,(Uint1)0,(Uint1)1,line_estimate[3],
															 A2F_OTHER, gbp);
		LoadPap(pap,
		 	PrintOrganismLine,ajp,0,(Uint1)0,(Uint1)0,line_estimate[11], 
		 													A2F_OTHER, gbp);
		pub_num = GetPubNum(gbp);
		for (index=0; index < pub_num; index++) {
			LoadPap(pap, 
				PrintPubsByNumber, ajp, index, (Uint1)0, (Uint1)0, 
				line_estimate[5], A2F_REFERENCE, gbp);
		}
		for (index=0; index < gbp->comm_num; index++) {
			if (index == 0) {
				LoadPap(pap, 
					PrintFirstComment, ajp, index, (Uint1)0, (Uint1)0, 
						line_estimate[5], A2F_COMMENT, gbp);
			} else {
				LoadPap(pap, 
					PrintCommentByNumber, ajp, index, (Uint1)0, (Uint1)0, 
						line_estimate[5], A2F_COMMENT, gbp);
			}
		}
		if (gbp->feat && gbp->feat->sfpCommsize > 0) {
			LoadPap(pap, GBDescrComFeat, ajp, 0, (Uint1)0, (Uint1)0, 
				line_estimate[2], A2F_OTHER, gbp);
		}
		if (gbp->xref_present == TRUE) {
			LoadPap(pap, 
				PrintXrefLine, ajp, 0, (Uint1)0, (Uint1)0, line_estimate[0],
															 A2F_OTHER, gbp);
		}
		if (gbp->feat) {
				LoadPap(pap, 	PrintFeatHeader, ajp, 0, (Uint1)0, (Uint1)0, 
					line_estimate[6], A2F_OTHER, gbp);
		}
		if (gbp->feat) {
			max = gbp->feat->sfpListsize;
		}
		LoadPap(pap, 
				PrintSourceFeat, ajp, 0, (Uint1)0, (Uint1)1, line_estimate[8],
													 A2F_SOURCE_FEATURE, gbp);
		for (index=0; index< max; index++) {
			if (max == index+1) {
				LoadPap(pap, 
					PrintNAFeatByNumber,ajp,index,(Uint1)0,
						(Uint1)1,line_estimate[8], A2F_FEATURE, gbp);
			} else {
				LoadPap(pap, 
					PrintNAFeatByNumber, ajp, index,(Uint1)0,
						(Uint1)0,line_estimate[8], A2F_FEATURE, gbp);
			}
		}
		LoadPap(pap, 
			PrintBaseCount, ajp, 0, (Uint1)0, (Uint1)0, line_estimate[0], 
																A2F_OTHER, gbp);
		seqblks_num = GetNumOfSeqBlks(ajp, gbp);
		for (index=0; index < seqblks_num; index++) {
			if (seqblks_num == index+1) {
				LoadPap(pap, 
				PrintSeqBlk, ajp, index, (Uint1)1, (Uint1)0, line_estimate[9],
															 A2F_OTHER, gbp);
			} else {
				LoadPap(pap, 
				PrintSeqBlk, ajp, index, (Uint1)0, (Uint1)0, line_estimate[9],
															 A2F_OTHER, gbp);
			}
		}
	}

	return total;
}

/*************************************************************************
*asn2gp_setup
*
*	This code calls the routines to output a GenPept Flat File
*
**************************************************************************/
Int4 asn2gp_setup(Asn2ffJobPtr ajp, FFPrintArrayPtr PNTR papp)
{

	BioseqPtr bsp;
	FFPrintArrayPtr pap;
	Int2 feat_num, pub_num;
	Int4 index, total;
	Int4 seqblks_num;
	GBEntryPtr gbp;
	SeqIdPtr sip;

	if (ajp->mode == SEQUIN_MODE) {
		GetLocusPartsAwp(ajp);
	} else {
		UseGIforLocus(ajp);
	}
	GetSeqFeat(ajp);

	total=0;
	for (gbp=ajp->asn2ffwep->gbp; gbp; gbp = gbp->next) {
		if (gbp->bsp && ajp->id_print) {
			sip = SeqIdFindBest(gbp->bsp->id, SEQID_GI);
			if (SeqIdComp(sip, ajp->id_print) != SIC_YES) {
				continue;
			}
		}
		CheckSourceFeat(ajp, gbp);
		bsp = gbp->bsp;
		if (ASN2FF_DROP_SHORT_AA == TRUE &&
			ajp->asn2ffwep->total_seg == 0 && bsp->length < GENPEPT_MIN) {
			flat2asn_delete_accession_user_string();
			flat2asn_delete_locus_user_string();
			flat2asn_install_accession_user_string(gbp->accession);
			flat2asn_install_locus_user_string(gbp->locus);
			if (ajp->error_msgs == TRUE)
				ErrPostStr(SEV_INFO, ERR_ENTRY_Partial_peptide, 
					"Entry dropped due to length.");
			continue;
		}
		total += 8;
		if (gbp->gi != -1 || gbp->ni != NULL) {
			total++;
		}
		if (ajp->asn2ffwep->total_seg > 0) {
			total++;
		}
		gbp->descr = NULL;	
		if (GP_GetSeqDescrComms(ajp, gbp) > 0) {
			total += gbp->comm_num;
		}
		if	(gbp->feat && gbp->feat->sfpCommsize > 0) {
			total++;
		}
		if (gbp->feat) {
				total += 2;				/* FEATURES and 'source' feature*/
				feat_num = gbp->feat->sfpListsize;
				total += feat_num;
		}
		seqblks_num = GetNumOfSeqBlks(ajp, gbp);
		total += seqblks_num;
		pub_num = (Int2)GetPubsAwp(ajp, gbp);
		total += pub_num; 

		GetGPDate(ajp, gbp);
	}

	*papp = (FFPrintArrayPtr) MemNew((size_t) total*sizeof(FFPrintArray));
	pap = *papp;

	LoadPap(NULL, NULL, ajp, 0, (Uint1)0, (Uint1)0, 0, A2F_OTHER, NULL);
	for (gbp=ajp->asn2ffwep->gbp; gbp; gbp = gbp->next) {
		bsp = gbp->bsp;
		if (bsp && ajp->id_print) {
			sip = SeqIdFindBest(bsp->id, SEQID_GI);
			if (SeqIdComp(sip, ajp->id_print) != SIC_YES) {
				continue;
			}
		}
		if (ASN2FF_DROP_SHORT_AA == TRUE &&
			ajp->asn2ffwep->total_seg == 0 && bsp->length < GENPEPT_MIN) {
			continue;
		}
		LoadPap(pap, 
			PrintLocusLine, ajp, 0, (Uint1)0, (Uint1)0, line_estimate[0],
																A2F_OTHER, gbp);
		flat2asn_delete_locus_user_string();
		flat2asn_install_locus_user_string(gbp->locus);
		GetDefinitionLine(ajp, gbp);
		LoadPap(pap,
			PrintDefinitionLine, ajp, 0, (Uint1)0, (Uint1)0, line_estimate[1], 
											A2F_OTHER, gbp);
                MemFree(gbp->descr);
		gbp->descr = NULL;
		LoadPap(pap, 
			PrintAccessLine, ajp, 0, (Uint1)0, (Uint1)0, line_estimate[2],
															 A2F_OTHER, gbp);
		flat2asn_delete_accession_user_string();
		flat2asn_install_accession_user_string(gbp->accession);
		if (gbp->gi != -1) {
			LoadPap(pap, PrintNCBI_GI, ajp, 0, (Uint1)0, (Uint1)0, 
				line_estimate[2], A2F_OTHER, gbp);
		} else if (gbp->ni != NULL) {
			LoadPap(pap, PrintNID, ajp, 0, (Uint1)0, (Uint1)0, 
				line_estimate[2], A2F_OTHER, gbp);
		}
		LoadPap(pap, 
			PrintDBSourceLine, ajp, 0, (Uint1)0, (Uint1)0, line_estimate[12], 
																A2F_OTHER, gbp);
		LoadPap(pap, 
			PrintKeywordLine, ajp, 0, (Uint1)0, (Uint1)0, line_estimate[3], 
																A2F_OTHER, gbp);
		if (ajp->asn2ffwep->total_seg > 0)
			LoadPap(pap, 
				PrintSegmentLine, ajp, 0, (Uint1)0, (Uint1)0, line_estimate[0], 
																A2F_OTHER, gbp);
		LoadPap(pap, PrintGBSourceLine, ajp, 0, (Uint1)0, (Uint1)0, 
			line_estimate[4], A2F_OTHER, gbp);
		LoadPap(pap, PrintGBOrganismLine, ajp, 0, (Uint1)0, (Uint1)0, 
			line_estimate[4], A2F_OTHER, gbp);
		pub_num = GetPubNum(gbp);
		for (index=0; index < pub_num; index++) {
			LoadPap(pap, 
				PrintPubsByNumber, ajp, index, (Uint1)0, (Uint1)0, 
				line_estimate[5], A2F_REFERENCE, gbp);
		}
		for (index=0; index < gbp->comm_num; index++) {
			if (index == 0) {
				LoadPap(pap, 
					PrintFirstComment, ajp, index, (Uint1)0, (Uint1)0, 
						line_estimate[5], A2F_COMMENT, gbp);
			} else {
				LoadPap(pap, 
					PrintCommentByNumber, ajp, index, (Uint1)0, (Uint1)0, 
						line_estimate[5], A2F_COMMENT, gbp);
			}
		}
		if (gbp->feat && gbp->feat->sfpCommsize > 0) {
			LoadPap(pap, GBDescrComFeat, ajp, 0, (Uint1)0, (Uint1)0, 
					line_estimate[2], A2F_OTHER, gbp);
		}
		if (gbp->feat) {
			LoadPap(pap, 	PrintFeatHeader, ajp, 0, (Uint1)0, (Uint1)0, 
				line_estimate[6], A2F_OTHER, gbp);
			LoadPap(pap, PrintSourceFeat, ajp, 0, (Uint1)0, (Uint1)0, 
				line_estimate[8], A2F_SOURCE_FEATURE, gbp);
			if (gbp->feat) {
				feat_num = gbp->feat->sfpListsize;
			}
			for (index=0; index < feat_num; index++) {
				LoadPap(pap, PrintAAFeatByNumber, ajp, index, 
					(Uint1)0, (Uint1)0, line_estimate[8], A2F_FEATURE, gbp);
			}
		}
		LoadPap(pap, 
			PrintOriginLine, ajp, 0, (Uint1)0, (Uint1)0, line_estimate[0], 
																A2F_OTHER, gbp);
		seqblks_num = GetNumOfSeqBlks(ajp, gbp);
		for (index=0; index < seqblks_num; index++) {
			if (seqblks_num == index+1) {
				LoadPap(pap, 
					PrintSeqBlk, ajp, index, (Uint1)1, (Uint1)0, 
											line_estimate[9], A2F_OTHER, gbp);
			} else {
				LoadPap(pap, 
					PrintSeqBlk, ajp, index, (Uint1)0, (Uint1)0, 
											line_estimate[9], A2F_OTHER, gbp);
			}
		}
	}

	return total;
} 

/*************************************************************************
*asn2ep_setup
*
*	This code calls the routines to output an "EMBLPept" Flat File
*
**************************************************************************/

Int4 asn2ep_setup(Asn2ffJobPtr ajp, FFPrintArrayPtr PNTR papp)
{

	BioseqPtr bsp;
	FFPrintArrayPtr pap;
	Int4 index, total;
	Int2 feat_num, pub_num;
	Int4 seqblks_num;
	GBEntryPtr gbp;


	ajp->format = GENPEPT_FMT;
	GetLocusPartsAwp(ajp);
	GetSeqFeat(ajp);
	ajp->format = EMBLPEPT_FMT;

	total=0;
	for (gbp=ajp->asn2ffwep->gbp; gbp; gbp = gbp->next) {
		CheckSourceFeat(ajp, gbp);
		bsp = gbp->bsp;
		if (ASN2FF_DROP_SHORT_AA == TRUE && ajp->asn2ffwep->total_seg == 0 && 
				bsp->length < GENPEPT_MIN) {
			flat2asn_delete_accession_user_string();
			flat2asn_delete_locus_user_string();
			flat2asn_install_accession_user_string(gbp->accession);
			flat2asn_install_locus_user_string(gbp->locus);
			if (ajp->error_msgs == TRUE)
				ErrPostStr(SEV_INFO, ERR_ENTRY_Partial_peptide, 
					"Entry dropped due to length.");
			continue;
		}
		total += 8;
		if (ajp->asn2ffwep->total_seg > 0) {
			total++;
		}
		if (GP_GetSeqDescrComms(ajp, gbp) > 0) {
			total += gbp->comm_num;
		}
		if	(gbp->feat && gbp->feat->sfpCommsize > 0) {
			total++;
		}
		if (gbp->feat) {
				total += 2;				/* FEATURES and 'source' feature*/
				feat_num = gbp->feat->sfpListsize;
				total += feat_num;
		}
		seqblks_num = GetNumOfSeqBlks(ajp, gbp);
		total += seqblks_num;
		pub_num = (Int2)GetPubsAwp(ajp, gbp);
		total += pub_num; 

		GetEMBLDate(ajp, gbp);
		GetEntryVersion(gbp);
	}

	*papp = (FFPrintArrayPtr) MemNew((size_t) total*sizeof(FFPrintArray));
	pap = *papp;

	LoadPap(NULL, NULL, ajp, 0, (Uint1)0, (Uint1)0, 0, A2F_OTHER, NULL);
	for (gbp=ajp->asn2ffwep->gbp; gbp; gbp = gbp->next) {
		bsp = gbp->bsp;
		if (ASN2FF_DROP_SHORT_AA == TRUE && ajp->asn2ffwep->total_seg == 0 && 
				bsp->length < GENPEPT_MIN) {
			continue;
		}
		LoadPap(pap, 
			PrintEPLocusLine, ajp, 0, (Uint1)0, (Uint1)1, line_estimate[0],
															 	A2F_OTHER, gbp);
		flat2asn_delete_locus_user_string();
		flat2asn_install_locus_user_string(gbp->locus);
		LoadPap(pap, 
			PrintAccessLine, ajp, 0, (Uint1)0, (Uint1)1, line_estimate[2],
															 	A2F_OTHER, gbp);
		flat2asn_delete_accession_user_string();
		flat2asn_install_accession_user_string(gbp->accession);
		LoadPap(pap, 
			PrintDateLines, ajp,  0, (Uint1)0, (Uint1)1, line_estimate[10], 
																A2F_OTHER, gbp);
		GetDefinitionLine(ajp, gbp);
		LoadPap(pap, 
			PrintDefinitionLine,  ajp, 0, (Uint1)0, (Uint1)1, line_estimate[1], 
																A2F_OTHER, gbp);
		LoadPap(pap, 
			PrintKeywordLine, ajp, 0, (Uint1)0, (Uint1)1, line_estimate[3], 
																A2F_OTHER, gbp);
		LoadPap(pap, 
			PrintOrganismLine, ajp, 0, (Uint1)0, (Uint1)0, line_estimate[11], 
																A2F_OTHER, gbp);
		pub_num = GetPubNum(gbp);
		for (index=0; index < pub_num; index++) {
			LoadPap(pap, PrintPubsByNumber, ajp, index, (Uint1)0, (Uint1)0, 
								line_estimate[5], A2F_REFERENCE, gbp);
		} 
		for (index=0; index < gbp->comm_num; index++) {
			if (index == 0) {
				LoadPap(pap, 
					PrintFirstComment, ajp, index, (Uint1)0, (Uint1)0, 
						line_estimate[5], A2F_COMMENT, gbp);
			} else {
				LoadPap(pap, 
					PrintCommentByNumber, ajp, index, (Uint1)0, (Uint1)0, 
						line_estimate[5], A2F_COMMENT, gbp);
			}
		}
		if (gbp->feat && gbp->feat->sfpCommsize > 0) {
			LoadPap(pap, GBDescrComFeat, ajp, 0, (Uint1)0, (Uint1)0, 
					line_estimate[2], A2F_OTHER, gbp);
		}
		if (gbp->feat) {
			LoadPap(pap, 	PrintFeatHeader, ajp, 0, (Uint1)0, (Uint1)0, 
				line_estimate[6], A2F_OTHER, gbp);
			LoadPap(pap, PrintSourceFeat, ajp, 0, (Uint1)0, (Uint1)0, 
				line_estimate[8], A2F_SOURCE_FEATURE, gbp);
			for (index=0; index < gbp->feat->sfpListsize; index++) {
				LoadPap(pap, PrintAAFeatByNumber, ajp, index, 
					(Uint1)0, (Uint1)0, line_estimate[8], A2F_FEATURE, gbp);
			}
		}
		seqblks_num = GetNumOfSeqBlks(ajp, gbp);
		for (index=0; index < seqblks_num; index++) {
			if (seqblks_num == index+1) {
				LoadPap(pap, 
					PrintSeqBlk, ajp, index, (Uint1)1, (Uint1)0, 
											line_estimate[9], A2F_OTHER, gbp);
			} else {
				LoadPap(pap, 
					PrintSeqBlk, ajp, index, (Uint1)0, (Uint1)0, 
											line_estimate[9], A2F_OTHER, gbp);
			}
		}
	}
	return total;
}	

static void FreeSortStructLoc(Int2 size, SortStructPtr p)
{
	Int2 index;
	Int2 size_loc;
	
				size_loc = p->extra_loc_cnt;
				if (size_loc > 0) {
					for (index=0; index < size_loc; index++) {
						SeqLocFree(p->extra_loc[index]);
					}
					MemFree(p->extra_loc);
				}
				for (index=0; index < size; index++) {
					if (p[index].feat_free == TRUE) {
						SeqFeatFree(p[index].sfp);
					}
				}
	return;
}

NLM_EXTERN void asn2ff_cleanup (Asn2ffJobPtr ajp)
{
	GBEntryPtr gbp, next;
	Int2 size, index;
	ValNodePtr v, vnext;
	PubStructPtr psp;
	ComStructPtr s, snext;
	SortStructPtr p;
	
	if (ajp->asn2ffwep != NULL) {
		for (gbp = ajp->asn2ffwep->gbp; gbp; gbp = next) {
			next = gbp->next;
			if (gbp->spp)
				SeqPortFree(gbp->spp);
			if (gbp->base_cnt_line)
				MemFree(gbp->base_cnt_line);
			if (gbp->feat) {
				size = gbp->feat->sfpListsize;
				p = gbp->feat->List;
				for (index=0; index < size; index++, p++) {
					if (p && p->gsp)
						GeneStructFree(p->gsp);
					if (p && p->nsp)
						NoteStructFree(p->nsp);
				}
				if (gbp->feat->sfpListsize > 0) {
					FreeSortStructLoc(gbp->feat->sfpListsize, 
							gbp->feat->List);
					MemFree(gbp->feat->List);
				}
				NoteStructFree(gbp->feat->source_notes);
				if (gbp->feat->sfpCommsize > 0) {
					FreeSortStructLoc(gbp->feat->sfpCommsize, 	
						gbp->feat->Commlist);
					MemFree(gbp->feat->Commlist);
				}
				if (gbp->feat->sfpGenesize > 0) {
					size = gbp->feat->sfpGenesize;
					p = gbp->feat->Genelist;
					for (index=0; index < size; index++, p++) {
					if (p && p->nsp)
						NoteStructFree(p->nsp);
				}
					FreeSortStructLoc(gbp->feat->sfpGenesize, 
						gbp->feat->Genelist);
					MemFree(gbp->feat->Genelist);
				}
				if (gbp->feat->sfpOrgsize > 0) {
					FreeSortStructLoc(gbp->feat->sfpOrgsize, 
						gbp->feat->Orglist);
					MemFree(gbp->feat->Orglist);
				}
				if (gbp->feat->sfpSitesize > 0) {
					FreeSortStructLoc(gbp->feat->sfpSitesize, 
									gbp->feat->Siteslist);
					MemFree(gbp->feat->Siteslist);
				}
				if (gbp->feat->sfpSourcesize > 0) {
					FreeSortStructLoc(gbp->feat->sfpSourcesize,
							 gbp->feat->Sourcelist);
					MemFree(gbp->feat->Sourcelist);
				}
				if (gbp->feat->sfpXrefsize > 0) {
					FreeSortStructLoc(gbp->feat->sfpXrefsize, 
						gbp->feat->Xreflist);
					MemFree(gbp->feat->Xreflist);
				}
				MemFree(gbp->feat);
			}
			for (v=gbp->Pub; v; v = vnext) {
				vnext = v->next;
				psp = (PubStructPtr) v->data.ptrvalue;
				FreePubStruct(psp);
				MemFree(v);
			}
			for (s=gbp->comm; s; s = snext) {
				snext = s->next;
				MemFree(s->string);
				MemFree(s);
			}
			if (gbp->source_info) {
				MemFree(gbp->source_info);
			}
			if (gbp->defline) {
				MemFree(gbp->defline);
			}
			MemFree(gbp);
		}
	}
	MemFree(ajp->asn2ffwep);
	SeqFeatFree(ajp->sfp_out);
/* Delete these strings so they don't interfere with others */
    flat2asn_delete_locus_user_string();
    flat2asn_delete_accession_user_string();

	return;
}


/*****************************************************************************
*	void LoadPap(FFPrintArrayPtr pap, FFPapFct fct, Asn2ffJobPtr ajp, 
*	Int4 index, Uint1 last, Uint1 printxx, Uint1 element_type) 
*
*	This function places the parameters in the correct spaces in the
*	FFPrintArrayPtr.
*
****************************************************************************/
void LoadPap(FFPrintArrayPtr pap, FFPapFct fct, Asn2ffJobPtr ajp, Int4 index, Uint1 last, Uint1 printxx, Int2 estimate, Uint1 element_type, GBEntryPtr gbp) 
{
	static Int4 pap_index;
	DescrStructPtr dsp;
	
	if (! pap) {
		pap_index=0;
	} else {
		pap[pap_index].fct = fct;
		pap[pap_index].ajp = ajp;
		pap[pap_index].gbp = gbp;
		pap[pap_index].index = index;
		pap[pap_index].last = last;
		pap[pap_index].printxx = printxx;
		pap[pap_index].estimate = estimate;
		pap[pap_index].descr = NULL;
		if (element_type == A2F_SOURCE_FEATURE) {
			dsp =  (DescrStructPtr) MemNew(sizeof(DescrStruct));
			pap[pap_index].descr = dsp;
			if (gbp->feat && gbp->feat->Sourcelist != NULL) {
				dsp->entityID = gbp->feat->Sourcelist[0].entityID;
				dsp->itemID = gbp->feat->Sourcelist[0].itemID;
				dsp->itemtype = gbp->feat->Sourcelist[0].itemtype;
			} else if (gbp->source_info != NULL) {
				dsp->entityID = gbp->source_info->entityID;
				dsp->itemID = gbp->source_info->itemID;
				dsp->itemtype = gbp->source_info->itemtype;
			}
		} else if (element_type == A2F_FEATURE && gbp->feat) {
			GetPapSeqFeatPtr (gbp, index, pap_index, pap);
		} else if (element_type == A2F_REFERENCE) {
			GetPapRefPtr (ajp, gbp, index, pap_index, pap);
		} else if (element_type == A2F_FEATURE_NEW && gbp->feat) {
			dsp =  (DescrStructPtr) MemNew(sizeof(DescrStruct));
			pap[pap_index].descr = dsp;
			dsp->entityID = gbp->feat->List[index].entityID;
			dsp->itemID = gbp->feat->List[index].itemID;
			dsp->itemtype = gbp->feat->List[index].itemtype;
		} else if (element_type == A2F_COMMENT) {
			GetPapCommPtr (ajp, gbp, index, pap_index, pap);
		} else {
			if (gbp->descr != NULL) {
				dsp =  (DescrStructPtr) MemNew(sizeof(DescrStruct));
				pap[pap_index].descr = dsp;
				dsp->entityID = gbp->descr->entityID;
				dsp->itemID = gbp->descr->itemID;
				dsp->itemtype = gbp->descr->itemtype;
			}
		}
		pap_index++;
	}

	return;
}

/****************************************************************************
*	This function checks a SeqPortPtr, maintained on the Biotable Ptr, 
*	and compares it's BioseqPtr to that of the BioseqPtr associated
*	with segment count of the btp.  At present, used only for nucleic
*	acids (4/14/94).
****************************************************************************/

void CheckSeqPort (Asn2ffJobPtr ajp, GBEntryPtr gbp, Int4 start)
{
	BioseqPtr bsp=gbp->bsp;
	SeqPortPtr spp=gbp->spp;
	Int4 start1;

	if (spp) {
		if (ajp->slp == NULL && bsp == spp->bsp) {
			if (spp->curpos != start)
				SeqPortSeek(spp, start, SEEK_SET);
		} else {
			SeqPortFree(spp);
			if (ajp->slp) {
				spp = SeqPortNewByLoc(ajp->slp, Seq_code_iupacna);
			} else {
				spp = SeqPortNew(bsp, 0, -1, 0, Seq_code_iupacna);
			}
			start1 = start - spp->start;
			if (start1 != spp->curpos)
				SeqPortSeek(spp, start1, SEEK_SET);
		}
	} else {
		if (ajp->slp) {
			spp = SeqPortNewByLoc(ajp->slp, Seq_code_iupacna);
		} else {
			spp = SeqPortNew(bsp, 0, -1, 0, Seq_code_iupacna);
		}
		start1 = start - spp->start;
		if (start1 != spp->curpos) 
			SeqPortSeek(spp, start1, SEEK_SET);
	}

	gbp->spp = spp;

	return;
}


/***************************************************************************
*
*	"GetMolInfo" gets information about the molecule for the locus
*	line.  The formatted information is in "buffer".
*
***************************************************************************/

void GetMolInfo (Asn2ffJobPtr ajp, CharPtr buffer, GBEntryPtr gbp)
{
	static CharPtr strand [4]= { "   ", "ss-", "ds-","ms-"};
	
	static CharPtr mol [9] = {"    ", "DNA ", "RNA ", "mRNA", "rRNA", "tRNA", "uRNA", "scRNA", " AA "};

	static CharPtr embl_mol [8] = {"xxx", "DNA", "RNA", "RNA", "RNA", "RNA", "RNA", "AA "};

	BioseqPtr bsp;
	Int2 istrand, imol;
	ValNodePtr vnp = NULL;
	MolInfoPtr mfp;
	Int4 length;
	bsp = gbp->bsp;
	istrand = bsp->strand;
	if (istrand > 3) 
		istrand = 0;

	imol = bsp->mol;
	if (imol > 3)
		imol = 0;

	if (ajp->slp) {
		length = SeqLocLen(ajp->slp);
	} else {
		length = bsp->length;
	}
/*keep both old and new style, get new first*/
	if ((vnp=GatherDescrByChoice(ajp, gbp, Seq_descr_molinfo)) != NULL) {
		mfp = vnp->data.ptrvalue;
		if (mfp->biomol <= 8) {
			imol = (Int2) (mfp->biomol);
		}
	} else {
		for (vnp = bsp->descr; vnp; vnp = vnp->next) {
			if (vnp->choice == Seq_descr_mol_type) {
				if (vnp->data.intvalue <= 8) {
					imol = (Int2) (vnp->data.intvalue);
				}
				break;
			}
		}
	}
	if (imol < 2) {  /* check inst.mol if mol-type is not-set or genomic */
		imol = bsp->mol;
		if (imol == 3)
			imol = 8;
		if (imol == 4)
			imol = 0;
	}

/* if ds-DNA don't show ds */
	if (imol == 1 && istrand == 2) { 
		istrand = 0;
	} 
/* ss-any RNA don't show ss */
	if (imol > 2 && istrand == 1) { 
		istrand = 0;
	} 
	if (ajp->slp != NULL) {
		bsp->topology = 1;
	}
	
	if (ajp->format == GENBANK_FMT || ajp->format == SELECT_FMT) {
		if (bsp->topology == 2) {
			sprintf(buffer, "%7ld bp %s%-4s  circular", 
							length, strand[istrand], mol[imol]);
		} else {
			sprintf(buffer, "%7ld bp %s%-4s          ", 
				length, strand[istrand], mol[imol]);
		}
	} else if (ajp->format == GENPEPT_FMT) {
			sprintf(buffer, "%7ld aa", length);
	} else if (ajp->format == EMBL_FMT || ajp->format == PSEUDOEMBL_FMT ||
		ajp->format == EMBLPEPT_FMT) {
		if (ajp->pseudo == FALSE) { /* do authentic EMBL */
			if (imol < 8) {
				if (bsp->topology == 2)
					sprintf(buffer, "circular %s", embl_mol[imol]);
				else
					sprintf(buffer, "%s", embl_mol[imol]);
			}
		} else {	/* Use GenBank molecule names */
			if (bsp->topology == 2)
				sprintf(buffer, "circular %s", mol[imol]);
			else
				sprintf(buffer, "%s", mol[imol]);
		}
	}
	return;
}

/*************************************************************************
*	Checks if there is a Xref in EMBL format.
*	Used ONLY to make EMBL output.
*This could probably be done more efficiently???????????????????
**************************************************************************/

Boolean CheckXrefLine (Asn2ffJobPtr ajp, GBEntryPtr gbp)

{
	Boolean ret_val=FALSE;
	Char buffer[20];
	CharPtr name;
	EMBLBlockPtr eb=NULL;
	EMBLXrefPtr xref=NULL;
	ValNodePtr descr=NULL, ds_vnp, tvnp;
	DescrStructPtr dsp;

	tvnp = GatherDescrListByChoice(ajp, gbp, Seq_descr_embl);
	for (descr=	tvnp;
				descr; descr=descr->next) {
		dsp = descr->data.ptrvalue;
		ds_vnp = dsp->vnp;
		eb = (EMBLBlockPtr) ds_vnp->data.ptrvalue;
		for (xref=eb->xref; xref; xref=xref->next) {
			name=NULL;
			if (xref->_class) {
				if (xref->_class == 5)
					StringCpy(buffer, "SWISS-PROT");
				else if (xref->_class == 8)
					StringCpy(buffer, "EPD");
				else if (xref->_class == 10)
					StringCpy(buffer, "TFD");
				else if (xref->_class == 11)
					StringCpy(buffer, "FLYBASE");
				name = &(buffer[0]);
			} else if (xref->name) {
				name = xref->name;
			}
			if (name && xref->id)
				ret_val = TRUE;
			else
				ret_val = FALSE;
		}
		MemFree(ds_vnp);
	}
	ValNodeFreeData(tvnp);
	return ret_val;
}

void PrintLocusLine(Asn2ffJobPtr ajp, GBEntryPtr gbp)
{
	BioseqPtr bsp;
	Char buffer[30];

	if (gbp == NULL)
		return;
	gbp->descr = NULL;
	bsp=gbp->bsp;
	if (ajp->format == EMBL_FMT || ajp->format == PSEUDOEMBL_FMT ||
				ajp->format == EMBLPEPT_FMT)
	{
		ff_StartPrint(5, 0, ASN2FF_EMBL_MAX, "ID");
		ff_AddString(gbp->locus);
		if (ajp->hup == TRUE) {
			ff_AddString(" confidential; ");
		} else {
			ff_AddString(" standard; ");
		}
		GetMolInfo(ajp, buffer, gbp);
		ff_AddString( buffer);
		ff_AddString("; ");
		if (ajp->ssp && ajp->format == EMBL_FMT && *(gbp->div) == ' ') {
			ff_AddString("UNA");			
		} else {
			ff_AddString(gbp->div);
		}
		ff_AddString("; ");
		if (ajp->slp) {
			ff_AddInteger("%ld", (long) SeqLocLen(ajp->slp));
		} else {
			ff_AddInteger("%ld", (long) bsp->length);
		}
		ff_AddString(" BP.");
		ff_EndPrint();
	} else {
		ff_StartPrint(0, 0, ASN2FF_GB_MAX, NULL);
		ff_AddString("LOCUS");
		TabToColumn(13);
		ff_AddString( gbp->locus);
		GetMolInfo(ajp, buffer, gbp);
		ff_AddString( buffer);
		TabToColumn(53);
		ff_AddString(gbp->div);
		TabToColumn(63);
		ff_AddString(gbp->date);
		ff_EndPrint();
	}
}

void PrintEPLocusLine(Asn2ffJobPtr ajp, GBEntryPtr gbp)

{
	BioseqPtr bsp=gbp->bsp;
	Char buffer[30];

	if (gbp == NULL)
		return;
	gbp->descr = NULL;	
	ff_StartPrint(5, 0, ASN2FF_EMBL_MAX, "ID");
	ff_AddString(gbp->locus);
	ff_AddString(" standard; ");
	GetMolInfo(ajp, buffer, gbp);
	ff_AddString(buffer);
	ff_AddString("; ");
	ff_AddString(gbp->div);
	ff_AddString("; ");
	ff_AddInteger("%ld", (long) bsp->length);
	ff_AddString(" RS.");
	ff_EndPrint();
}


void PrintAccessLine(Asn2ffJobPtr ajp, GBEntryPtr gbp)

{

	if (gbp == NULL)
		return;
	gbp->descr = NULL;
	if (ajp->format == EMBL_FMT || ajp->format == PSEUDOEMBL_FMT ||
				ajp->format == EMBLPEPT_FMT)
	{
		ff_StartPrint(5, 5, ASN2FF_EMBL_MAX, "AC");
	}
	else
	{
		ff_StartPrint(0, 12, ASN2FF_GB_MAX, NULL);
		ff_AddString( "ACCESSION");
		TabToColumn(13);
	}
	if (ajp->ssp && ajp->hup) {
		ff_AddChar(';');
	} else {
		ff_AddString(gbp->accession);
	}
	AddExtraAccessions(ajp, gbp);
	ff_EndPrint();
	return;
}

void PrintNCBI_GI(Asn2ffJobPtr ajp, GBEntryPtr gbp)

{

	if (gbp == NULL)
		return;
	gbp->descr = NULL;

	if (ajp->format == EMBL_FMT || ajp->format == PSEUDOEMBL_FMT ||
				ajp->format == EMBLPEPT_FMT) {
		ff_StartPrint(5, 5, ASN2FF_EMBL_MAX, "NI");
	} else {
		ff_StartPrint(0, 12, ASN2FF_GB_MAX, NULL);
		if (ajp->format == GENBANK_FMT) {
			ff_AddString( "NID");
		} else if (ajp->format == GENPEPT_FMT) {
			ff_AddString( "PID");
		}
		TabToColumn(13);
	}
	ff_AddChar('g');
	ff_AddInteger("%ld", (long) gbp->gi);
/*		if (ajp->format == GENBANK_FMT) {
			TabToColumn(26);
			ff_AddString( "GI:");
			ff_AddInteger("%ld", (long) gbp->gi);
		}
*/
	ff_EndPrint();
	if (ajp->format == EMBL_FMT || ajp->format == PSEUDOEMBL_FMT ||
				ajp->format == EMBLPEPT_FMT) {
		PrintXX();
	}
	return;
}
void PrintNID(Asn2ffJobPtr ajp, GBEntryPtr gbp)

{

	if (gbp == NULL)
		return;
	gbp->descr = NULL;

	if (ajp->format == EMBL_FMT || ajp->format == PSEUDOEMBL_FMT ||
				ajp->format == EMBLPEPT_FMT) {
		ff_StartPrint(5, 5, ASN2FF_EMBL_MAX, "NI");
		ff_AddString(gbp->ni);
	} else {
		ff_StartPrint(0, 12, ASN2FF_GB_MAX, NULL);
		if (ajp->format == GENBANK_FMT) {
			ff_AddString( "NID");
		} else if (ajp->format == GENPEPT_FMT) {
			ff_AddString( "PID");
		}
		TabToColumn(13);
		ff_AddString(gbp->ni);
/*		if (gbp->gi != -1) {
			TabToColumn(26);
			ff_AddString( "GI:");
			ff_AddInteger("%ld", (long) gbp->gi);
		}
*/
	}
	ff_EndPrint();
	if (ajp->format == EMBL_FMT || ajp->format == PSEUDOEMBL_FMT ||
				ajp->format == EMBLPEPT_FMT) {
		PrintXX();
	}
	return;
}
void PrintDateLines (Asn2ffJobPtr ajp, GBEntryPtr gbp)

{
	if (gbp == NULL)
		return;
	gbp->descr = NULL;
	ff_StartPrint(5, 5, ASN2FF_EMBL_MAX, "DT");
	if (gbp->update_date)
	{	/* both create and update date exist.	*/
		if (ajp->pseudo == FALSE)
		{ /* In pseudo-EMBL mode only one date line */
			if (gbp->create_date) {
				ff_AddString(gbp->create_date);
				NewContLine();
			}
		}
		ff_AddString(gbp->update_date);
		if (gbp->embl_rel)
		{
			ff_AddString(" (Rel. ");
			ff_AddString(gbp->embl_rel);
			ff_AddString(", Last updated, Version ");
			ff_AddInteger("%ld", (long) gbp->embl_ver);
			ff_AddChar(')');
		}
	}
	else
	{	/* only create date exists.	*/
		ff_AddString(gbp->create_date);
		if (gbp->embl_rel)
		{
			ff_AddString(" (Rel. ");
			ff_AddString(gbp->embl_rel);
			ff_AddString(", Last updated, Version ");
			ff_AddInteger("%ld", (long) gbp->embl_ver);
			ff_AddChar(')');
		}
		if (ajp->pseudo == FALSE)
		{ /* In pseudo-EMBL only one date line. */
			NewContLine();
			ff_AddString(gbp->create_date);
			if (gbp->embl_rel)
			{
				ff_AddString(" (Rel. ");
				ff_AddString(gbp->embl_rel);
				ff_AddString(", Last updated, Version ");
				ff_AddInteger( "%ld", (long) gbp->embl_ver);
				ff_AddChar(')');
			}
		}
	}
	ff_EndPrint();
}	/* PrintDateLines */

void PrintSegmentLine(Asn2ffJobPtr ajp, GBEntryPtr gbp)
{

	if (gbp == NULL)
		return;
	gbp->descr = NULL;
	if (ajp->asn2ffwep->total_seg > 1)
	{
		ff_StartPrint(0, 0, ASN2FF_GB_MAX, NULL);
		ff_AddString("SEGMENT");
		TabToColumn(13);
		ff_AddInteger("%ld", (long) gbp->num_seg);
		ff_AddString(" of ");
		ff_AddInteger("%ld", (long) ajp->asn2ffwep->total_seg);
		ff_EndPrint();
	}
}

static ValNodePtr AddKeyword(ValNodePtr key, CharPtr add)
{
	ValNodePtr vnp;
	
	for (vnp = key; vnp; vnp = vnp->next) {
		if (StringCmp(vnp->data.ptrvalue, add) == 0) {
			return key;
		}
	}
	vnp = ValNodeNew(NULL);
	vnp->data.ptrvalue = StringSave(add);
	key = tie_next(key, vnp);
	
	return key;
}

static Boolean CheckSpecialKeyword(Boolean is_est, Boolean is_sts, Boolean is_gss, CharPtr kwd)
{
	if (is_est == FALSE && is_sts == FALSE && is_gss == FALSE) {
		return TRUE;
	}
	if (is_est) {
		if (MatchArrayString(STS_kw_array, TOTAL_STSKW, kwd) != -1) {
			return FALSE;
		}
		if (MatchArrayString(GSS_kw_array, TOTAL_GSSKW, kwd) != -1) {
			return FALSE;
		}
		return TRUE;
	}
	if (is_sts) {
		if (MatchArrayString(EST_kw_array, TOTAL_ESTKW, kwd) != -1) {
			return FALSE;
		}
		if (MatchArrayString(GSS_kw_array, TOTAL_GSSKW, kwd) != -1) {
			return FALSE;
		}
		return TRUE;
	}
	if (is_gss) {
		if (MatchArrayString(STS_kw_array, TOTAL_STSKW, kwd) != -1) {
			return FALSE;
		}
		if (MatchArrayString(EST_kw_array, TOTAL_ESTKW, kwd) != -1) {
			return FALSE;
		}
		return TRUE;
	}
	return TRUE;
}

ValNodePtr GetKeywordLine(Asn2ffJobPtr ajp, GBEntryPtr gbp)

{
	ValNodePtr block, keyword=NULL, v, vnp;
	GBBlockPtr gbblk;
	EMBLBlockPtr ebp;
	PirBlockPtr pbp;
	PrfBlockPtr prfp;
	SPBlockPtr spbp;
	MolInfoPtr mfp;
	Boolean is_est=FALSE, is_sts=FALSE, is_gss=FALSE;

	if ((block = GatherDescrByChoice(ajp, gbp, Seq_descr_molinfo)) != NULL) 
	{
		mfp = (MolInfoPtr) block->data.ptrvalue;
			switch (mfp->tech) {
				case MI_TECH_htgs_1:
					keyword = ValNodeNew(NULL);
					keyword->data.ptrvalue = StringSave("HTG");
					v = ValNodeNew(keyword);
					v->data.ptrvalue = StringSave("HTGS_PHASE1");
				break;
				case MI_TECH_htgs_2:
					keyword = ValNodeNew(NULL);
					keyword->data.ptrvalue = StringSave("HTG");
					v = ValNodeNew(keyword);
					v->data.ptrvalue = StringSave("HTGS_PHASE2");
				break;
				case MI_TECH_htgs_3:
					keyword = ValNodeNew(NULL);
					keyword->data.ptrvalue = StringSave("HTG");
				break;
				case MI_TECH_est:
					is_est = TRUE;
					keyword = ValNodeNew(NULL);
					keyword->data.ptrvalue = StringSave("EST");
				break;
				case MI_TECH_sts:
					is_sts = TRUE;
					keyword = ValNodeNew(NULL);
					keyword->data.ptrvalue = StringSave("STS");
				break;
				case MI_TECH_survey:
					is_gss = TRUE;
					keyword = ValNodeNew(NULL);
					keyword->data.ptrvalue = StringSave("GSS");
				break;
				case MI_TECH_fli_cdna:
					keyword = ValNodeNew(NULL);
					keyword->data.ptrvalue = StringSave("FLI_CDNA");
				break;
				default:
					break;
			}
	}

	if ((block = GatherDescrByChoice(ajp, gbp, Seq_descr_genbank)) != NULL)
	{
		gbblk = (GBBlockPtr) block->data.ptrvalue;
		if (gbblk->keywords != NULL) {
			for (vnp = gbblk->keywords; vnp; vnp = vnp->next) {
				if (CheckSpecialKeyword(is_est, is_sts, is_gss, vnp->data.ptrvalue) == TRUE) {
					keyword = AddKeyword(keyword, vnp->data.ptrvalue);
				}
			}
			return keyword;
		} else {
			if (gbp->descr) {
				gbp->descr = MemFree(gbp->descr);
			}
		}
	}
	if ((block = GatherDescrByChoice(ajp, gbp, Seq_descr_embl)) != NULL)
	{
		ebp = (EMBLBlockPtr) block->data.ptrvalue;
		if (ebp->keywords != NULL) {
			for (vnp = ebp->keywords; vnp; vnp = vnp->next) {
				if (CheckSpecialKeyword(is_est, is_sts, is_gss, vnp->data.ptrvalue) == TRUE) {
					keyword = AddKeyword(keyword, vnp->data.ptrvalue);
				}
			}
			return keyword;
		} else {
			if (gbp->descr) {
				gbp->descr = MemFree(gbp->descr);
			}
		}
	} 
	if ((block = GatherDescrByChoice(ajp, gbp, Seq_descr_pir)) != NULL)
	{
		pbp = (PirBlockPtr) block->data.ptrvalue;
		if (pbp->keywords != NULL) {
			for (vnp = pbp->keywords; vnp; vnp = vnp->next) {
				keyword = AddKeyword(keyword, vnp->data.ptrvalue);
			}
			return keyword;
		} else {
			if (gbp->descr) {
				gbp->descr = MemFree(gbp->descr);
			}
		}
	}
	if ((block = GatherDescrByChoice(ajp, gbp, Seq_descr_prf)) != NULL) 
	{
		prfp = (PrfBlockPtr) block->data.ptrvalue;
		if (prfp->keywords != NULL) {
			for (vnp = prfp->keywords; vnp; vnp = vnp->next) {
				keyword = AddKeyword(keyword, vnp->data.ptrvalue);
			}
			return keyword;
		} else {
			if (gbp->descr) {
				gbp->descr = MemFree(gbp->descr);
			}
		}
	}
	if ((block = GatherDescrByChoice(ajp, gbp, Seq_descr_sp)) != NULL) 
	{
		spbp = (SPBlockPtr) block->data.ptrvalue;
		if (spbp->keywords != NULL) {
			for (vnp = spbp->keywords; vnp; vnp = vnp->next) {
				keyword = AddKeyword(keyword, vnp->data.ptrvalue);
			}
			return keyword;
		} else {
			if (gbp->descr) {
				gbp->descr = MemFree(gbp->descr);
			}
		}
	}
	return keyword;

}	/* GetKeywordLine */


void PrintKeywordLine(Asn2ffJobPtr ajp, GBEntryPtr gbp)
{
	Boolean line_return = FALSE;
	Boolean first = TRUE;
	CharPtr string;
	Int2 tab_length=12;
	ValNodePtr keyword, vnp;

	gbp->descr = NULL;
	keyword = GetKeywordLine(ajp, gbp);

	if (ajp->format == EMBL_FMT || ajp->format == PSEUDOEMBL_FMT ||
			ajp->format == EMBLPEPT_FMT) {
		ff_StartPrint(5, 5, ASN2FF_EMBL_MAX, "KW");
	} else {
		ff_StartPrint(0, tab_length, ASN2FF_GB_MAX, NULL);
		ff_AddString("KEYWORDS");
		TabToColumn((Int2)(tab_length+1));
	}
	if (keyword != NULL) {	/* the next line initializes the length */
		for (vnp=keyword; vnp != NULL; vnp=vnp->next) {
			string = vnp->data.ptrvalue;
			if (first == TRUE) {
				first = FALSE;
			} else {
				if (line_return)
					NewContLine();
			}
	
			ff_AddString(string);
			if (vnp->next != NULL) {
				ff_AddChar(';');
				ff_AddChar(' ');
			}
		}
		ValNodeFreeData(keyword);
	} else if (gbp->descr) {
		MemFree(gbp->descr);
		gbp->descr = NULL;
	}
	ff_AddChar('.');
	ff_EndPrint();
	

}	/* PrintKeywordLine */

void PrintDefinitionLine(Asn2ffJobPtr ajp, GBEntryPtr gbp)
{
	if (ajp->format == EMBL_FMT || ajp->format == PSEUDOEMBL_FMT ||
				ajp->format == EMBLPEPT_FMT) {
		ff_StartPrint(5, 5, ASN2FF_EMBL_MAX, "DE");
	} else {
		ff_StartPrint(0, 12, ASN2FF_GB_MAX, NULL);
		ff_AddString("DEFINITION");
		TabToColumn(13);
	}
	ff_AddString(gbp->defline);
	ff_EndPrint();
	return;
}

void GetDefinitionLine(Asn2ffJobPtr ajp, GBEntryPtr gbp)
{
	CharPtr string, string_start, title=NULL;
	ValNodePtr vnp = NULL;
	MolInfoPtr mfp;
	CharPtr buf;
	Int2 buflen = 1001;
	ItemInfoPtr iip;
	DescrStructPtr dsp = NULL;
	Uint1 tech = 0;
	
	buf = MemNew(buflen+1);
	gbp->descr = NULL;
/*  deflines for htg sequences */
	vnp=GatherDescrByChoice(ajp, gbp, Seq_descr_molinfo);
	if (vnp != NULL) {
		mfp = (MolInfoPtr)vnp->data.ptrvalue;
		if (mfp) {
			tech = mfp->tech;
		}
	}
	if (gbp && gbp->descr) {
		gbp->descr = MemFree(gbp->descr);
	}
	
	iip = MemNew(sizeof(ItemInfo));
	CreateDefLine(iip, gbp->bsp, buf, buflen, tech, NULL, NULL);
	if (iip != NULL) {
		dsp = MemNew(sizeof(DescrStruct));
		dsp->entityID = iip->entityID;
		dsp->itemID = iip->itemID;
		dsp->itemtype = iip->itemtype;
	}
	MemFree(iip);
	gbp->descr = dsp;
	title = buf;
	string_start = string = CheckEndPunctuation(title, '.');

	while (*string != '\0')
	{
		if (*string == '\"')
			*string = '\'';
		string++;
	}

	gbp->defline = StringSave(string_start);
	MemFree(string_start);
	MemFree(buf);
}

void PrintOriginLine(Asn2ffJobPtr ajp, GBEntryPtr gbp)

{
	Char buffer[68];
	CharPtr origin=NULL;
	GBBlockPtr gb;
	Int2 length=0;
	ValNodePtr vnp=NULL;

	gbp->descr = NULL;	
	ff_StartPrint(0, 12, ASN2FF_GB_MAX, NULL);
	ff_AddString("ORIGIN");
	TabToColumn(13);
	if ((vnp=GatherDescrByChoice(ajp, gbp, Seq_descr_genbank)) != NULL)
	{
		gb = (GBBlockPtr) vnp->data.ptrvalue;
		if (gb)
		{
			if (gb->origin && (length=StringLen(gb->origin)) > 0)
			{ /*???? What if gb->origin is longer than 68 chars. */
				StringNCpy(buffer, gb->origin, 66);
				if (length < 66)
					buffer[length] = '\0';
				else
					buffer[66] = '\0';
				origin = CheckEndPunctuation(buffer, '.');
				ff_AddString(origin);
			}
			if (length > 66)
				ErrPostStr(SEV_WARNING, ERR_ENTRY_OriginTooLong, "");
		}
	}
	if (origin != NULL)
		MemFree(origin);
	ff_EndPrint();

}
static void print_source(Asn2ffJobPtr ajp, CharPtr source, OrgRefPtr orp)
{
	CharPtr		newsource, s;
	Boolean 	has_point = FALSE;
	ValNodePtr 	v;

	ff_StartPrint(0, 12, ASN2FF_GB_MAX, NULL);
	ff_AddString("SOURCE");
	TabToColumn(13);
	if (source) {
		newsource = CheckEndPunctuation(source, '.');
		ff_AddString(newsource);
		MemFree(source);
		MemFree(newsource);
	} else if (orp) {
		source = orp->common?orp->common:orp->taxname;
		ff_AddString(source);
		if (orp->mod == NULL && source != NULL) {
			if (*(source + StringLen(source) -1) == '.')
				has_point = TRUE;
		}
		for (v = orp->mod; v; v = v->next) {
			has_point = FALSE;
			s = (CharPtr) (v->data.ptrvalue);
			if (*(s + StringLen(s) -1) == '.')
				has_point = TRUE;
			ff_AddString(" ");
			ff_AddString(s);
			
		}
		if (!has_point)
			ff_AddChar('.');
	} else {
		ff_AddString("Unknown.");
		if (ajp->error_msgs == TRUE)
			ErrPostStr(SEV_WARNING, ERR_ENTRY_No_source_line, "");
	}
	ff_EndPrint();
}

/***************************************************************************
*PrintGBSourceLine
*
*	"PrintGBSourceLine" to print the source ONLY line for 
*	genbank and genpept FlatFiles. (modified from PrintSourceLine)
*
****************************************************************************/
void PrintGBSourceLine (Asn2ffJobPtr ajp, GBEntryPtr gbp)
{
	CharPtr 		source=NULL;
	GBBlockPtr 		gb=NULL;
	OrgRefPtr 		orp=NULL;
	BioSourcePtr 	biosp;
	ValNodePtr 		vnp=NULL;
	SeqFeatPtr 		sfp;
	SortStructPtr 	p;

	
	if (gbp == NULL) {
		return;
	}
	gbp->descr = NULL;	
	if ((vnp=GatherDescrByChoice(ajp, gbp, Seq_descr_genbank)) != NULL) {
		gb = (GBBlockPtr) vnp->data.ptrvalue;
		if (gb)
			source = GetGBSourceLine(gb);
	}
	if (source) {
		print_source(ajp, source, NULL);
		return;
	}
	if (gbp->descr) {
		gbp->descr = MemFree(gbp->descr);
	}
	if ((vnp=GatherDescrByChoice(ajp, gbp, Seq_descr_source)) != NULL) {
		biosp = (BioSourcePtr) vnp->data.ptrvalue;
		if (biosp != NULL) {
			orp = biosp->org;
		}
	} else if ((vnp=GatherDescrByChoice(ajp, gbp, Seq_descr_org)) != NULL) {
		orp = (OrgRefPtr) vnp->data.ptrvalue;
	} else if (gbp->feat && gbp->feat->sfpOrgsize != 0) {
		p = gbp->feat->Orglist;
		if ((sfp = p->sfp) == NULL) {
		GatherItemWithLock(p->entityID,
			p->itemID, p->itemtype, &sfp, find_item);
		}
		if (sfp != NULL) {
			orp = (OrgRefPtr) sfp->data.value.ptrvalue;
		}
	}
	print_source(ajp, NULL, orp);
	return;
}

static void print_organism(Asn2ffJobPtr ajp, GBEntryPtr gbp, OrgRefPtr orp, CharPtr lineage)
{
	DbtagPtr dbp;
	Int4	id = -1;
	CharPtr organelle, taxonomy=NULL;

	if (orp) {
		if(orp->common && !orp->taxname)
			orp->taxname = TaxNameFromCommon(orp->common);
		if (lineage == NULL && orp->orgname) {
			lineage = orp->orgname->lineage;
		}
	}
	organelle = FlatOrganelle(ajp, gbp);
	ff_StartPrint(2, 12, ASN2FF_GB_MAX, NULL);
	ff_AddString("ORGANISM");
	TabToColumn(13);
	if (orp && orp->taxname) {
		if (organelle) {
			ff_AddString(organelle);
		}
		if (orp->db != NULL) {
			dbp = (DbtagPtr) (orp->db)->data.ptrvalue;
			if (StringCmp(dbp->db, "taxon") == 0)
				id = dbp->tag->id;
		}
		www_organism(orp->taxname, id);
	} else {
		ff_AddString("Unknown.");
	}
	MemFree(organelle);
	ff_EndPrint();

	ff_StartPrint(12, 12, ASN2FF_GB_MAX, NULL);
	if (lineage) {
		taxonomy = CheckEndPunctuation(lineage, '.');
		ff_AddString(taxonomy);
		MemFree(taxonomy);
	} else {
		ff_AddString("Unclassified.");
	}
	ff_EndPrint();
}

/***************************************************************************
*PrintGBOrganismLine
*
*	"PrintGBOrganismLine" to print the ONLY organism field for 
*	genbank and genpept FlatFiles. (modified from PrintSourceLine)
*
****************************************************************************/
void PrintGBOrganismLine (Asn2ffJobPtr ajp, GBEntryPtr gbp)
{
	CharPtr 		lineage = NULL;
	GBBlockPtr 		gb=NULL;
	OrgRefPtr 		orp=NULL;
	BioSourcePtr 	biosp;
	ValNodePtr		vnp=NULL;
	SeqFeatPtr 		sfp;
	SortStructPtr 	p;
	DescrStructPtr 	dsp;

	
	if (gbp == NULL) {
		return;
	}
	gbp->descr = NULL;
/* Find Biosource with focus 
	v=GatherDescrListByChoice(ajp, gbp, Seq_descr_source);
	for (tvnp=v; tvnp; tvnp=tvnp->next) {
		dsp = (DescrStructPtr) tvnp->data.ptrvalue;
		vnp = dsp->vnp;
		biosp = (BioSourcePtr) vnp->data.ptrvalue;
		if (biosp->is_focus == TRUE) {
			orp = biosp->org;
			gbp->descr = dsp; 
			break;
		}
	}
	ValNodeFreeData(v);
	if (orp == NULL && gbp->feat && gbp->feat->biosrcsize != 0) { 
		p = gbp->feat->Biosrclist;
		for (i = 0; i < gbp->feat->biosrcsize; i++, p++) {
			if ((sfp = p->sfp) == NULL) {
				GatherItemWithLock(p->entityID,
					p->itemID, p->itemtype, &sfp, find_item);
			}
			if (sfp != NULL) {
				biosp = (BioSourcePtr) sfp->data.value.ptrvalue;
				if (biosp->is_focus == TRUE) {
					orp = biosp->org;
					dsp = MemNew(sizeof(DescrStruct));
					gbp->descr = dsp;
					dsp->entityID = p->entityID;
					dsp->itemID = p->itemID;
					dsp->itemtype = p->itemtype;
					break;
				}
			}
		}
	}
	if (orp != NULL) { 
		if (orp->orgname) {
			lineage = orp->orgname->lineage;
		}
		print_organism(ajp, gbp, orp, lineage);
		return;
	}
*/
/* BioSource descr*/
	if ((vnp=GatherDescrByChoice(ajp, gbp, Seq_descr_source)) != NULL) 
	{
		biosp = vnp->data.ptrvalue;
		orp = (OrgRefPtr) biosp->org;
		if (orp && orp->orgname) {
			lineage = orp->orgname->lineage;
		}
	} 
/* try to find lineage in GenBank block */
	if (lineage == NULL) {
		dsp = gbp->descr;
		gbp->descr = NULL;
		if ((vnp=GatherDescrByChoice(ajp, gbp, Seq_descr_genbank)) != NULL) {
			gb = (GBBlockPtr) vnp->data.ptrvalue;
			if (gb)
				lineage = gb->taxonomy;
		}
		gbp->descr = MemFree(gbp->descr);
		gbp->descr = dsp;  /* keep Seq_descr_source dsp for sequin */
	}
	if (orp) {
		print_organism(ajp, gbp, orp, lineage);
		return;
	}
/* Organism descr*/
	if ((vnp=GatherDescrByChoice(ajp, gbp, Seq_descr_org)) != NULL) {
		orp = (OrgRefPtr) vnp->data.ptrvalue;
		print_organism(ajp, gbp, orp, lineage);
		return;
	}
/* OrgRef feature */
	gbp->descr = MemFree(gbp->descr);
	if (gbp->feat && gbp->feat->sfpOrgsize != 0) {
		p = gbp->feat->Orglist;
		if ((sfp = p->sfp) == NULL) {
		GatherItemWithLock(p->entityID,
			p->itemID, p->itemtype, &sfp, find_item);
		}
		if (sfp != NULL) {
			orp = (OrgRefPtr) sfp->data.value.ptrvalue;
			dsp = MemNew(sizeof(DescrStruct));
			gbp->descr = dsp;
			dsp->entityID = p->entityID;
			dsp->itemID = p->itemID;
			dsp->itemtype = p->itemtype;
		}
	}
	print_organism(ajp, gbp, orp, lineage);
	return;
}

/***************************************************************************
*PrintOrganismLine
*
*	"PrintOrganismLine" to print the source and organism entries for 
*	EMBL FlatFiles.
*
*Rewrite for better logic???? (11/30/93 & 12/13/93)
*To take care of pdb, pir cases???? (12/13/93)
*Note: two organism lines are searched for!!
****************************************************************************/

void PrintOrganismLine (Asn2ffJobPtr ajp, GBEntryPtr gbp)

{
	ValNodePtr vnp=NULL;
	OrgRefPtr orp=NULL, orp1=NULL;
	CharPtr organelle, taxonomy=NULL, lineage = NULL;
	BioSourcePtr biosp = NULL;
	GBBlockPtr gb=NULL;
	DescrStructPtr dsp;
	ValNodePtr tvnp;
	SeqFeatPtr	sfp;
	SortStructPtr p;

	if (gbp == NULL) {
		return;
	}	
/* new first */
	gbp->descr = NULL;	
	if ((tvnp=GatherDescrListByChoice(ajp, gbp, Seq_descr_source)) != NULL) {
		dsp = (DescrStructPtr) tvnp->data.ptrvalue;
		vnp = dsp->vnp;
		biosp = (BioSourcePtr) vnp->data.ptrvalue;
		orp = (OrgRefPtr) biosp->org;
		if (tvnp->next != NULL) {
			dsp = (DescrStructPtr) tvnp->next->data.ptrvalue;
			vnp = dsp->vnp;
			biosp = (BioSourcePtr) vnp->data.ptrvalue;
			orp1 = (OrgRefPtr) biosp->org;
		}
		ValNodeFreeData(tvnp);
	}
	if (orp && orp->orgname) {
		lineage = orp->orgname->lineage;
	}
	if (orp == NULL) {
		if ((tvnp=GatherDescrListByChoice(ajp, gbp, Seq_descr_org)) != NULL) {
			dsp = (DescrStructPtr) tvnp->data.ptrvalue;
			vnp = dsp->vnp;
			orp = (OrgRefPtr) vnp->data.ptrvalue;
			if (tvnp->next != NULL) {
				dsp = (DescrStructPtr) tvnp->next->data.ptrvalue;
				vnp = dsp->vnp;
				orp1 = (OrgRefPtr) (vnp->data.ptrvalue);
			}
			ValNodeFreeData(tvnp);
		} else if (gbp->feat && gbp->feat->sfpOrgsize != 0) {
			p = gbp->feat->Orglist;	/* gbp->feat->Orglist[0] */
			if ((sfp = p->sfp) == NULL) {
				GatherItemWithLock(p->entityID,
					p->itemID, p->itemtype, &sfp, find_item);
			}
			if (sfp != NULL) {
				orp = (OrgRefPtr) sfp->data.value.ptrvalue;
			}
			p++;   /* gbp->feat->Orglist[1] */
			if ((sfp = p->sfp) == NULL) {
				GatherItemWithLock(p->entityID,
					p->itemID, p->itemtype, &sfp, find_item);
			}
			if (sfp != NULL) {
				orp1 = (OrgRefPtr) sfp->data.value.ptrvalue;
			}
		}
	}

	if (orp)
		if(orp->common && !orp->taxname)
			orp->taxname = TaxNameFromCommon(orp->common);

	ff_StartPrint(5, 5, ASN2FF_EMBL_MAX, "OS");
	if (orp && orp->taxname)
	{
		ff_AddString(orp->taxname);
		if (orp->common)
		{
			ff_AddString(" (");
			ff_AddString(orp->common);
			ff_AddChar(')');
		}
	}
	else
		ff_AddString("Unclassified.");

	ff_EndPrint();

	ff_StartPrint(5, 5, ASN2FF_EMBL_MAX, "OC");
	if (lineage == NULL) {
		if ((vnp=GatherDescrByChoice(ajp, gbp, Seq_descr_genbank)) != NULL){
			gb = (GBBlockPtr) vnp->data.ptrvalue;
			lineage = gb->taxonomy;
		}
	}
	if (lineage) {
		taxonomy = CheckEndPunctuation(lineage, '.');
		ff_AddString(taxonomy);
		MemFree(taxonomy);
	} else {
		ff_AddString("Unclassified.");
	}
	ff_EndPrint();
	
	if (orp1) {	/* second organism */
		if (orp1 && orp1->taxname) {
			PrintXX();
			ff_StartPrint(5, 5, ASN2FF_EMBL_MAX, "OS");
			ff_AddString(orp1->taxname);
			if (orp1->common) {
				ff_AddString(" (");
				ff_AddString(orp1->common);
				ff_AddChar(')');
			}
			ff_EndPrint();
		}
	}

/* What about plasmids on the OG line???????????????*/
/* Get this info from a qual of the SourceFeat that has qual "plasmid"??*/
	organelle = FlatOrganelle(ajp, gbp);
	if (organelle) {
		PrintXX();
		ff_StartPrint(5, 5, ASN2FF_EMBL_MAX, "OG");
		ff_AddString(organelle);
		ff_EndPrint();
		MemFree(organelle);
	}

}	/* PrintOrganismLine */



/****************************************************************************
*GetPDBSourceLine
*
*	Gets the source from the PDBBlock.
*
****************************************************************************/

CharPtr GetPDBSourceLine (PdbBlockPtr pdb)

{
	CharPtr source = NULL;
	ValNodePtr vnp;

	if(pdb && pdb->source)
	{
		vnp = pdb->source;
		source = StringSave(vnp->data.ptrvalue);
	}

	return source;
}

/***********************************************************************
*       This function prints out a block of the sequence (at most
*       of size SEQ_BLK_SIZE).
*       After the last sequence block, the terminator is printed also.
***********************************************************************/
 
void PrintSeqBlk (Asn2ffJobPtr ajp, GBEntryPtr gbp)
 
{
        Int4 start, stop, index=ajp->pap_index;
        Uint1 last=ajp->pap_last;
		DescrStructPtr dsp;

		dsp = MemNew(sizeof(DescrStruct));
		gbp->descr = dsp;
		dsp->entityID = gbp->entityID;
		dsp->itemID = gbp->itemID;
		dsp->itemtype = gbp->itemtype;
        if (index == 0) {
			if (ajp->slp != NULL) {
				start = SeqLocStart(ajp->slp);
			} else {
				start = 0;
            }    
        } else {
			if (ajp->slp != NULL) {
				start = index*SEQ_BLK_SIZE + SeqLocStart(ajp->slp);
			} else {
                start = index*SEQ_BLK_SIZE;
            }
        }
        if (last != LAST) {
			if (ajp->slp != NULL) {
                stop = SeqLocStart(ajp->slp) + (index+1)*SEQ_BLK_SIZE - 1;
            } else {
                stop = (index+1)*SEQ_BLK_SIZE - 1;
            }
        } else {
			if (ajp->slp != NULL) {
				stop = SeqLocStop(ajp->slp);
            } else {
                stop = -1;
            }
        }
        if (ajp->format == EMBLPEPT_FMT) {
        	PrintEPSequence(ajp, gbp, start, stop);
        } else {
        	PrintSequence(ajp, gbp, start, stop);
        }
        if (last ==  LAST)
                PrintTerminator();
}

void PrintPubsByNumber (Asn2ffJobPtr ajp, GBEntryPtr gbp)

{
	PubStructPtr psp;
	ValNodePtr vnp;
	Int2 i, index = (Int2) ajp->pap_index;

	for (vnp=gbp->Pub, i=0; vnp && i < index; vnp=vnp->next, i++);
	if (vnp) {
		psp = vnp->data.ptrvalue;
		if (ajp->format == EMBL_FMT || ajp->format == PSEUDOEMBL_FMT ||
			ajp->format == EMBLPEPT_FMT) {
			EMBL_PrintPubs(ajp, gbp, psp);
		} else {
			if (ajp->mode == PARTIAL_MODE && psp->choice == PUB_Sub) {
				return;
			} else {
				GB_PrintPubs(ajp, gbp, psp);
			}
		}
	}
}
void PrintFeatHeader (Asn2ffJobPtr ajp, GBEntryPtr gbp)

{
	gbp->descr = NULL;
	if (ajp->format == EMBL_FMT || ajp->format == PSEUDOEMBL_FMT ||
		ajp->format == EMBLPEPT_FMT) {
		PrintXX();
		ff_StartPrint( 5, 0, ASN2FF_EMBL_MAX, "FH");
		ff_AddString("Key");
		TabToColumn(22);
		ff_AddString("Location/Qualifiers");
		NewContLine();
	} else {
		ff_StartPrint(0, 0, ASN2FF_GB_MAX, NULL);
		ff_AddString("FEATURES");
		TabToColumn(22);
		ff_AddString("Location/Qualifiers");
	}
	ff_EndPrint();
}


/**************************************************************************
*void PrintTerminator ()
*
*	Prints the double slash (//) at the end of an entry.
**************************************************************************/

void PrintTerminator (void)

{
	ff_StartPrint(0, 0, 0, NULL);
	ff_AddChar( '/');
	ff_AddChar('/');
	ff_EndPrint();
}

/*************************************************************************
*	Prints out the cross-refs from the EMBL block, in the descriptor.
*	Used ONLY to make EMBL output.
**************************************************************************/

void PrintXrefLine (Asn2ffJobPtr ajp, GBEntryPtr gbp)

{
	Boolean done_once=FALSE;
	Char buffer[20], buffer1[20], buffer2[20];
	CharPtr name, string;
	EMBLBlockPtr eb=NULL;
	EMBLXrefPtr xref=NULL;
	ObjectIdPtr oip;
	ValNodePtr descr=NULL, id;

	gbp->descr = NULL;	
	if ((descr=GatherDescrByChoice(ajp, gbp, Seq_descr_embl)) != NULL) 
	{
		eb = (EMBLBlockPtr) descr->data.ptrvalue;
		for (xref=eb->xref; xref; xref=xref->next)
		{
			name=NULL;
			if (xref->_class) {
				if (xref->_class == 5)
					StringCpy(buffer, "SWISS-PROT");
				else if (xref->_class == 8)
					StringCpy(buffer, "EPD");
				else if (xref->_class == 10)
					StringCpy(buffer, "TFD");
				else if (xref->_class == 11)
					StringCpy(buffer, "FLYBASE");
				name = &(buffer[0]);
			}
			else if (xref->name)
				name = xref->name;
			if (name && xref->id)
			{
				id=xref->id;
			
				oip = id->data.ptrvalue;
				if (oip->str)
					StringCpy(buffer1, oip->str);
				else if (oip->id)
					sprintf(buffer1, "%ld", (long) (oip->id));
				id = id->next;
				if (id)
				{
					oip = id->data.ptrvalue;
					if (oip->str)
						StringCpy(buffer2, oip->str);
					else if (oip->id)
						sprintf(buffer2, "%ld", (long) (oip->id));
				}
				if (done_once == FALSE) {
					PrintXX();
					ff_StartPrint(5, 5, ASN2FF_EMBL_MAX, "DR");
					done_once=TRUE;
				} else {
					NewContLine();
				}
				ff_AddString(name);
				ff_AddString("; ");
				ff_AddString(buffer1);
				ff_AddString("; ");
				string = CheckEndPunctuation(buffer2, '.');
				ff_AddString(string);
				string = MemFree(string);
			}
		}
	}
	if (done_once)
	{
		ff_EndPrint();
	/*	PrintXX();*/
	}
}
/****************************************************************************
*
*	"PrintBaseCount" counts and prints the number of a, c, g, t, and 
*	other in a sequence.
*
****************************************************************************/

void PrintBaseCount (Asn2ffJobPtr ajp, GBEntryPtr gbp)

{
	CharPtr buffer;
	Int4 base_count[5], total=0;
	SeqPortPtr spp = NULL;
	Uint1 residue;
	DescrStructPtr dsp;
	BioseqPtr bsp = gbp->bsp;
	
	dsp = (DescrStructPtr) MemNew(sizeof(DescrStruct));
	gbp->descr = dsp;
	dsp->entityID = 0;
	dsp->itemID = 0;
	dsp->itemtype = 0;		
	if (gbp->base_cnt_line)
	{	/* Been there (at least once), done that.	*/
		buffer = gbp->base_cnt_line;
	} else {
		base_count[0]=0;
		base_count[1]=0;
		base_count[2]=0;
		base_count[3]=0;
		base_count[4]=0;

		if (ajp->slp) {
			spp = SeqPortNewByLoc(ajp->slp, Seq_code_iupacna);
		} else {
			spp = SeqPortNew(gbp->bsp, 0, -1, 0, Seq_code_iupacna);
		}
		if (bsp->repr == Seq_repr_delta) {
			SeqPortSet_do_virtual(spp, TRUE);
		}
		while ((residue=SeqPortGetResidue(spp)) != SEQPORT_EOF)
		{
			if ( !IS_residue(residue) && residue != INVALID_RESIDUE )
				continue;
	
			total++;
			switch (residue) {
				case 'A':
					base_count[0]++;
					break;
				case 'C':
					base_count[1]++;
					break;
				case 'G':
					base_count[2]++;
					break;
				case 'T':
					base_count[3]++;
					break;
				default:
						base_count[4]++;
						break;
			}
		}
		buffer = (CharPtr) MemNew(80*sizeof(Char));
		if (ajp->format == EMBL_FMT || ajp->format == PSEUDOEMBL_FMT ||
			ajp->format == EMBLPEPT_FMT)
		{
			sprintf(buffer, 
			"%ld BP; %ld A; %ld C; %ld G; %ld T; %ld other;", 
			(long) total, (long) base_count[0], (long) base_count[1], 
			(long) base_count[2], (long) base_count[3], (long) base_count[4]);
		}
		else 	/* GENBANK format */
		{
			if (base_count[4] == 0)
			{
			sprintf(buffer, 
			"%7ld a%7ld c%7ld g%7ld t", 
			(long) base_count[0], (long) base_count[1], 
				(long) base_count[2], (long) base_count[3]);
			}
			else
			{
			sprintf(buffer, 
			"%7ld a%7ld c%7ld g%7ld t%7ld others", 
			(long) base_count[0], (long) base_count[1], 
			(long) base_count[2], (long) base_count[3], (long) base_count[4]);
			}
		}
		gbp->base_cnt_line = buffer;
		if (spp) {
			SeqPortFree(spp);
		}
	}

		if (ajp->format == EMBL_FMT || ajp->format == PSEUDOEMBL_FMT ||
			ajp->format == EMBLPEPT_FMT)
	{
		ff_StartPrint(5, 5, ASN2FF_EMBL_MAX, "SQ");
		ff_AddString("Sequence ");
		ff_AddString(buffer);
	}
	else
	{
		ff_StartPrint(0, 0, ASN2FF_GB_MAX, NULL);
		ff_AddString("BASE COUNT");
		TabToColumn(13);
		ff_AddString( buffer);
	}

	ff_EndPrint();
}	/* PrintBaseCount */
/*****************************************************************************
*
*	"PrintSequence" to get the biological sequence (in iupacna or
*	iupacaa format)	and put it in a buffer suitable for Genbank 
*	or EMBL format.
*
*	The variables "start" and "stop" allow one to read from a point 
*	not at the beginning of the sequence to a point not at the end
*	of the sequence.
*
*	Rewrite to store in a buffer and print out more at once????????
*****************************************************************************/

void PrintSequence (Asn2ffJobPtr ajp, GBEntryPtr gbp, Int4 start, Int4 stop)

{
	BioseqPtr bsp=gbp->bsp;
	Char buffer[MAX_BTP_BUF], num_buffer[10];
	CharPtr ptr = &(buffer[0]), num_ptr;
	Int4 index, inner_index, inner_stop, total=start;
	Int4 loc_start;
	SeqPortPtr spp;
	Uint1 residue;


	if ((loc_start = SeqLocStart(ajp->slp)) == -1) {
		loc_start = 0;
	}
	if (ajp->format == GENBANK_FMT || ajp->format == SELECT_FMT)
	{
		ff_StartPrint(0, 0, ASN2FF_GB_MAX, NULL);
		sprintf(ptr, "%9ld ", (long) (total+1 - loc_start));
		ptr += StringLen(ptr);
		CheckSeqPort(ajp, gbp, start);
		spp = gbp->spp;
		if (bsp->repr == Seq_repr_delta) {
			SeqPortSet_do_virtual(spp, TRUE);
		}
		if (stop == -1)
			stop = spp->stop;
		for (index=start; index<=stop; index += 10) {
			if (stop < (index+10)) {
				inner_stop = stop;
			} else {
				inner_stop = index+9;
			}
			for (inner_index=index; inner_index<=inner_stop; inner_index++) {
				residue=SeqPortGetResidue(spp);
/*********/
				if (ajp->only_one) {
					if (residue == SEQPORT_VIRT) {
						*ptr = '\0';
						ff_AddString(buffer);
						NewContLine();
						MemSet(buffer, ' ', ptr - buffer);
						inner_index--;
						continue;
					}
				}
/**********/		
				if ( !IS_residue(residue) && residue != INVALID_RESIDUE) {
					continue;
				}
				if (residue == INVALID_RESIDUE) {
					residue = (Uint1) 'X';
				}
				*ptr++ = TO_LOWER(residue);
			}
			total = inner_stop+1;
		/* Put in a space every ten, unless it's the end of a row. */
			if (ROUNDUP(total-start, 60) == (total-start)) {
				if (total != (start+1) && total != (stop+1)) {
					*ptr = '\0';
					ptr = &buffer[0];
					ff_AddString(ptr);
					NewContLine();
					sprintf(ptr, "%9ld ", (long) (total+1 - loc_start));
					ptr += StringLen(ptr);
				}
			}
			else if (ROUNDUP(total-start, 10) == total-start)
			{
				*ptr = ' '; ptr++;
			}
		}
		*ptr = '\0';
		ptr = &buffer[0];
		ff_AddString( ptr);
	}
	else if (ajp->format == GENPEPT_FMT)
	{
		total++;

		ff_StartPrint(0, 0, ASN2FF_GB_MAX, NULL);
		sprintf(ptr, "%9ld ", (long) (total - loc_start));
		ptr += StringLen(ptr);
		if (ASN2FF_IUPACAA_ONLY == TRUE)
			spp = SeqPortNew(bsp, start, stop, 0, Seq_code_iupacaa);
		else
			spp = SeqPortNew(bsp, start, stop, 0, Seq_code_ncbieaa);
		while ((residue=SeqPortGetResidue(spp)) != SEQPORT_EOF)
		{
			if ( !IS_residue(residue) && residue != INVALID_RESIDUE )
				continue;
			if (residue == INVALID_RESIDUE)
				residue = (Uint1) 'X';
	
			*ptr = TO_LOWER(residue); ptr++;
			if (ROUNDUP(total - start, 10) == total - start) 
			{
				if (ROUNDUP(total - start, 60) == total - start)
				{
					if (total != (start+1) && total != (stop+1))
					{
						*ptr = '\0';
						ptr = &buffer[0];
						ff_AddString(ptr);
						NewContLine();
						num_ptr = &(num_buffer[0]);
						sprintf(num_ptr, "%9ld", (long) (total+1 - loc_start));
						while ((*ptr = *num_ptr) != '\0')
						{
							ptr++; num_ptr++;
						}
						*ptr = ' '; ptr++;
					}
				}
				else
				{
					*ptr = ' '; ptr++;
				}
			}
			total++;
		}
		*ptr = '\0';
		ptr = &buffer[0];
		ff_AddString(ptr);
		SeqPortFree(spp);
	}
	else if (ajp->format == EMBL_FMT || ajp->format == PSEUDOEMBL_FMT ||
			ajp->format == EMBLPEPT_FMT)
	{	/* numbers at far right, let line go to MAX_BTP_BUF */

		ff_StartPrint(5, 5, 0, NULL);
		CheckSeqPort(ajp, gbp, start);
		spp = gbp->spp;
		if (stop == -1)
			stop = spp->stop;
		for (index=start; index<=stop; index += 10)
		{
			if (stop < (index+10))
				inner_stop = stop;
			else
				inner_stop = index+9;
			for (inner_index=index; inner_index<=inner_stop; inner_index++)
			{
				residue=SeqPortGetResidue(spp);
				if ( !IS_residue(residue) && residue != INVALID_RESIDUE )
					continue;
				if (residue == INVALID_RESIDUE)
					residue = (Uint1) 'X';
		
				*ptr = TO_LOWER(residue); ptr++;
			}
			total = inner_index;
			if (ROUNDUP(total - start, 10) == total - start) 
			{
				if (ROUNDUP(total - start, 60) == total - start)
				{
					*ptr = '\0';
					ptr = &buffer[0];
					ff_AddString(ptr);
					TabToColumn(73);
					ff_AddInteger("%8ld", (long) (total - loc_start));
					if (ROUNDUP(total, SEQ_BLK_SIZE) != total)
						NewContLine();
				}
				else
				{
					*ptr = ' '; ptr++;
				}
			}
		}
		total = stop+1;
		if (ROUNDUP(total - start, 60) != total - start)
		{
			*ptr = '\0';
			ptr = &buffer[0];
			ff_AddString(ptr);
			TabToColumn(73);
			ff_AddInteger("%8ld", (long) (total - loc_start));
		}
	}

	ff_EndPrint();


}	/* PrintSequence */

/*****************************************************************************
*
*	"PrintEPSequence" to get the biological sequence (in iupacna or
*	iupacaa format)	and put it in a buffer suitable for Genbank 
*	or EMBL format.
*
*	The variables "start" and "stop" allow one to read from a point 
*	not at the beginning of the sequence to a point not at the end
*	of the sequence.
*
*	Rewrite to store in a buffer and print out more at once????????
*****************************************************************************/

void PrintEPSequence (Asn2ffJobPtr ajp, GBEntryPtr gbp, Int4 start, Int4 stop)

{
	BioseqPtr bsp=gbp->bsp;
	Char buffer[MAX_BTP_BUF];
	CharPtr ptr = &(buffer[0]);
	Int4 index, inner_index, inner_stop, total=start;
	SeqPortPtr spp;
	Uint1 residue;


	/* numbers at far right, let line go to MAX_BTP_BUF */

	ff_StartPrint(5, 5, 0, NULL);
	if (ASN2FF_IUPACAA_ONLY == TRUE)
		spp = SeqPortNew(bsp, start, stop, 0, Seq_code_iupacaa);
	else
		spp = SeqPortNew(bsp, start, stop, 0, Seq_code_ncbieaa);
	if (stop == -1)
		stop = spp->stop;
	for (index=start; index<=stop; index += 10)
	{
		if (stop < (index+10))
			inner_stop = stop;
		else
		inner_stop = index+9;
		for (inner_index=index; inner_index<=inner_stop; inner_index++)
		{
			residue=SeqPortGetResidue(spp);
			if ( !IS_residue(residue) && residue != INVALID_RESIDUE )
				continue;
			if (residue == INVALID_RESIDUE)
				residue = (Uint1) 'X';
	
			*ptr = TO_LOWER(residue); ptr++;
		}
		total = inner_index;
		if (ROUNDUP(total, 10) == total) 
		{
			if (ROUNDUP(total, 60) == total)
			{
				*ptr = '\0';
				ptr = &buffer[0];
				ff_AddString(ptr);
				TabToColumn(73);
				ff_AddInteger("%8ld", (long) total);
				if (ROUNDUP(total, SEQ_BLK_SIZE) != total)
					NewContLine();
			}
			else
			{
				*ptr = ' '; ptr++;
			}
		}
	}
	total = stop+1;
	if (ROUNDUP(total, 60) != total)
	{
		*ptr = '\0';
		ptr = &buffer[0];
		ff_AddString(ptr);
		TabToColumn(73);
		ff_AddInteger("%8ld", (long) total);
	}

	ff_EndPrint();

	SeqPortFree(spp);


}	/* PrintEPSequence */

void GatherItemWithLock(Uint2 entityID, Uint2 itemID, Uint2 itemtype,
                                   Pointer userdata, GatherItemProc userfunc)
{
	GatherItem(entityID, itemID, itemtype, userdata,  userfunc);
    return;
}

Boolean find_item (GatherContextPtr gcp)
{
	SeqFeatPtr sfp;
	SeqFeatPtr PNTR sfpp;


	sfpp = gcp->userdata;
	switch (gcp->thistype) {
		case OBJ_SEQFEAT:
			sfp = (SeqFeatPtr) (gcp->thisitem);
			*sfpp = sfp;
		break;
		default:
		break;
	}
	return TRUE;
}
