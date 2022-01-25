/*   accentr.h
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
* RCS $Id: accentr.h,v 6.0 1997/08/25 18:12:31 madden Exp $
*
* Author:  Ostell
*
* Version Creation Date:   4/23/92
*
* File Description: 
*       entrez index access library for Entrez
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
* 06-29-94 Schuler     Structure definitions and prototypes for EntrezInfo 
*                      AsnRead/AsnWrite functions moved to objentr.h
* 08-16-94 Brylawski   Prototypes for EntrezHierarchyGet, on-the-fly
*                      text neighboring, and enabling of medline features
* 09-08-94 Schuler     Define for TYP_ST
* 10-06-94 Schuler     Added #include <objmmdb1.h>
* 10-06-94 Schuler     Added EntrezBiostrucGet() 
* 11-20-94 Brylawski   Moved EntrezNeighborText structure to here from 
*                      objneten.h .
* 11-25-94 Schuler     Moved CdTerm structure to here from cdromlib.h 
* 04-11-95 Schuler     Defined TYP_ES, TYP_CH; changed value of TYP_SEQ
* 05-16-94 Schuler     Removed TYP_ES (decided not to use it)
* 05-16-94 Schuler     Added RCS Log directive.  Comments will henceforth
*                      be inserted automatically into the source at the
*                      time it is checked into the system. 
$Log: accentr.h,v $
Revision 6.0  1997/08/25 18:12:31  madden
Revision changed to 6.0

Revision 5.10  1997/07/28 13:30:41  ostell
Moved GetUniGeneIDForSeqId() to seqmgr.c

Revision 5.9  1997/06/26 21:55:14  vakatov
[PC] DLL'd "ncbicdr.lib", "ncbiacc.lib", "ncbinacc.lib" and "ncbicacc.lib"

Revision 5.8  1996/12/02 22:20:52  brandon
*** empty log message ***

 * Revision 5.7  1996/12/02  22:14:52  brandon
 * *** empty log message ***
 *
 * Revision 5.6  1996/12/02  22:10:46  brandon
 * *** empty log message ***
 *
 * Revision 5.5  1996/09/23  21:21:33  shavirin
 * Added definition for the function AccessionToFasta
 *
 * Revision 5.4  1996/08/14  19:56:09  epstein
 * add APIs for fetching pieces of biostruc annots
 *
 * Revision 5.3  1996/08/14  17:08:29  brandon
 * *** empty log message ***
 *
 * Revision 5.1  1996/08/14  15:15:05  brandon
 * added date parameter to tleval functions
 *
 * Revision 5.0  1996/05/28  13:55:34  ostell
 * Set to revision 5.0
 *
 * Revision 4.8  1996/05/14  21:02:10  epstein
 * adjust field descriptions to reflect current reality
 *
 * Revision 4.7  1996/03/29  18:52:51  epstein
 * add support for structure alignments
 *
 * Revision 4.6  1996/03/11  21:51:03  ostell
 * made GetUniGeneIDForSeqId() externally visible
 *
 * Revision 4.5  1995/10/11  13:39:29  epstein
 * add EntrezIsInited() function
 *
 * Revision 4.4  1995/10/02  02:35:50  epstein
 * add range-checking
 *
 * Revision 4.3  1995/08/28  23:20:47  kans
 * Biostruc_supported restored
 *
 * Revision 4.2  1995/08/21  19:42:16  epstein
 * add cluster analysis
 *
 * Revision 4.1  1995/08/11  20:25:53  epstein
 * add max-models support for biostrucs
 *
 * Revision 4.0  1995/07/26  13:50:32  ostell
 * force revision to 4.0
 *
 * Revision 2.44  1995/07/25  18:47:57  kans
 * revert to no Biostruc_supported
 *
 * Revision 2.43  1995/07/05  19:15:46  ostell
 * changed prototype for BiostrucAvail to be included only once
 *
 * Revision 2.42  1995/06/29  15:57:46  epstein
 * added Complexity argument when fetching structures
 *
 * Revision 2.41  95/06/26  18:07:13  kans
 * restored #define Biostruc_supported
 * 
 * Revision 2.39  1995/06/23  16:02:43  kans
 * support for accmmdbs.c stub to resolve symbols without MMDB link
 *
 * Revision 2.38  1995/06/22  21:21:22  kans
 * #define Biostruc_supported
 *
 * Revision 2.37  1995/05/16  14:12:46  schuler
 * Automatic comment insertion enabled
 *


* ==========================================================================
*/

#ifndef _ACCENTR_
#define _ACCENTR_

#ifndef Biostruc_supported
#define Biostruc_supported
#endif

#ifndef _NCBI_Seqset_
#include <objsset.h>
#endif

#ifndef _NCBI_Medline_
#include <objmedli.h>
#endif

#ifndef _NCBI_Access_
#include <objacces.h>
#endif

#include <objentr.h>

#ifdef _PMENTREZ_
#include <pmaccs.h>
#endif

#ifdef Biostruc_supported
#include <mmdbapi1.h>
#include <mmdbapi2.h>
#include <mmdbapi3.h>
#else
typedef struct { int bogus; } Biostruc, *BiostrucPtr;

/*****************************************************************************
*
*   BiostrucAvail is in the accmmdbs.c stub (returning FALSE) and in
*      bios2ff7.c (returning TRUE).  The stub file also has BiostrucAsnRead
*      and BiostrucAsnWrite empty functions for proper linking without MMDB.
*
*****************************************************************************/

NLM_EXTERN Boolean LIBCALL BiostrucAvail PROTO((void));
#endif


#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* --- Type Definitions --- */

typedef DocUid *DocUidPtr;
typedef Int2   DocType;
typedef Int2   DocField;

typedef struct CdTerm {
	DocType type;                 /* class of term */
	DocField field;               /* field of term */
	CharPtr term;                 /* the term */
	CharPtr highRange;            /* the top end of a "range" */
	Int4 special_count ,
		total_count;
	Int4 offset;                  /* offset into the postings file */
	Int2 page;
	struct CdTerm *next;
} CdTerm, *CdTermPtr;

typedef struct FastaSeq {
  Int4 gi;
  CharPtr label;
  CharPtr seq;
} FastaSeq, PNTR FastaSeqPtr;


/*****************************************************************************
*
*   PreDefined Entrez types and fields
*
*****************************************************************************/

/*--- doc type codes  ---
 *
 *	Although NTYPE is _NOT_ the correct number of document types, its
 *	value must not be changed or the code in cdromlib.c will break.
 *	Since the number of types is really determined at run-time, not
 *	compile-time, a more appropriate way to determine this value is to
 *	examine the field in the EntrezInfo structure.
 */

#define NTYPE       3   /* DO NOT CHANGE !! */

#define TYP_ML      0   /*  MEDLINE */
#define TYP_AA      1   /*  Amino Acid sequence */
#define TYP_NT      2   /*  Nucleotide sequence */
#define TYP_ST      3   /*  3-D Structure (not supported by cdromlib.c) */
#define TYP_CH      4   /*  Chromosome (not supported by cdromlib.c) */

#define TYP_SEQ     127   /*  either aa or na used only for uid lookups */


/*--- field codes ---*/
#define NFLD        19  /* number of fields == 19 */
#define FLD_WORD    0   /*  Words           */
#define FLD_MESH    1   /*  MeSH terms      */
#define FLD_KYWD    2   /*  Keyword         */
#define FLD_AUTH    3   /*  Authors         */
#define FLD_JOUR    4   /*  Journal title   */
#define FLD_ORGN    5   /*  Organism        */
#define FLD_ACCN    6   /*  Accession number */
#define FLD_GENE    7   /*  Gene Symbol     */
#define FLD_PROT    8   /*  Protein name    */
#define FLD_ECNO    9   /*  E.C. number     */
#define FLD_ORGN_HIER 10 /* Organism hierarchy */
#define FLD_DATE   11   /* date of entry*/
#define FLD_FKEY   12   /* Feature key  */
#define FLD_PROP   13   /* Properties */
#define FLD_SUBS   14   /* Substance */
#define FLD_DATM   15   /* Modification Date */
#define FLD_SLEN   16   /* Sequence Length */
#define FLD_AFFL   17   /* Affilication */
#define FLD_SQID   18   /* SeqId */

#define FLD_MESH_HIER 10 /* Overload for MeSH */

#define TERM__EXPLODE 1


typedef Boolean (*DocSumListCallBack) PROTO((DocSumPtr dsp, DocUid uid));


typedef struct EntrezHierarchyChild {
    CharPtr name;
    Boolean isLeafNode;
    Int4 special;
    Int4 total;
} EntrezHierarchyChild, *EntrezHierarchyChildPtr;

typedef struct EntrezHierarchy {
    Int2 numInLineage;
    Int2 numChildren;
    DocType db;
    DocField fld;
    CharPtr term;
    CharPtr PNTR lineage;
    EntrezHierarchyChildPtr children;
    CharPtr canonicalForm;
} EntrezHierarchy, *EntrezHierarchyPtr;

typedef struct EntrezNeighborText {
   Int4   fld;
   Int4   percent_terms_to_use;
   Int4   max_neighbors;
   Int4   min_score;
   CharPtr   normalText;
   CharPtr   specialText;
} EntrezNeighborText, PNTR EntrezNeighborTextPtr;


#define NULLSYM     0     /* for building booleans */
#define LPAREN      1
#define RPAREN      2
#define ANDSYMBL    3
#define ORSYMBL     4
#define BUTNOTSYMBL 5
#define SPECIALTERM 6
#define TOTALTERM   7


/**** Initialize and close session *********************/

/* Note:                                                                      */
/*   The EntrezInitWithExtras() capability is provided for backwards          */
/*   compatability until May, 1994.  This name was introduced as a temporary  */
/*   interface to EntrezInit(), to smooth the transition for EntrezInit()'s   */
/*   use of three parameters; previously it took no parameters.               */
/*                                                                            */
/*                                  - J. Epstein, 17 Feb 1994                 */

#define EntrezInitWithExtras(a,b,c) EntrezInit(a,b,c)
NLM_EXTERN Boolean LIBCALL EntrezInit PROTO((CharPtr appl_id, Boolean no_warnings, BoolPtr is_network));
NLM_EXTERN Boolean LIBCALL EntrezIsInited PROTO((void));
NLM_EXTERN void LIBCALL EntrezFini PROTO((void));

/**** Get names and numbers of fields and types ********/

NLM_EXTERN EntrezInfoPtr LIBCALL EntrezGetInfo PROTO((void));

/**** Creates a term node from the uid parameter ********/
NLM_EXTERN void LIBCALL EntrezCreateNamedUidList PROTO((CharPtr term, DocType type, DocField field, Int4 num, DocUidPtr uids));
NLM_EXTERN void LIBCALL EntrezCreateNamedUidListX PROTO((CharPtr term, DocType type, DocField field, ByteStorePtr bsp));

/**** Get detailed text information about the current status *****/
NLM_EXTERN CharPtr LIBCALL EntrezDetailedInfo PROTO((void));

/**** Get Links and Neighbors **************************/

NLM_EXTERN Int4 LIBCALL EntrezGetMaxLinks PROTO((void));
NLM_EXTERN Int4 LIBCALL EntrezSetUserMaxLinks PROTO((Int4 usermax));
NLM_EXTERN Int4 LIBCALL EntrezGetUserMaxLinks PROTO((void));

NLM_EXTERN LinkSetPtr LIBCALL EntrezUidLinks PROTO((DocType type, DocUid uid, DocType link_to_type));
NLM_EXTERN Int2 LIBCALL EntrezLinkUidList PROTO((LinkSetPtr PNTR result, DocType type, DocType link_to_type, Int2 numuid, Int4Ptr uids, Boolean mark_missing));

/**** Get Summaries ************************************/

NLM_EXTERN DocSumPtr LIBCALL EntrezDocSum PROTO((DocType type, DocUid uid));
NLM_EXTERN DocSumPtr LIBCALL DocSumFree PROTO((DocSumPtr dsp));
NLM_EXTERN Int2 LIBCALL EntrezDocSumListGet PROTO((Int2 numuid, DocType type, DocUidPtr uids, DocSumListCallBack callback));
NLM_EXTERN Int2 LIBCALL EntrezMlSumListGet PROTO((DocSumPtr PNTR result, Int2 numuid, Int4Ptr uids));
NLM_EXTERN Int2 LIBCALL EntrezSeqSumListGet PROTO((DocSumPtr PNTR result, Int2 numuid, DocType type, Int4Ptr uids));

/**** Get Term Lists ***********************************/

typedef Boolean (*TermListProc) PROTO((CharPtr term, Int4 special, Int4 total));

NLM_EXTERN Int2 LIBCALL EntrezTermListByPage PROTO((DocType type, DocField field, Int2 page, Int2 numpage, TermListProc proc));
NLM_EXTERN Int2 LIBCALL EntrezTermListByTerm PROTO((DocType type, DocField field, CharPtr term, Int2 numterms, TermListProc proc, Int2Ptr first_page));

NLM_EXTERN Boolean LIBCALL EntrezFindTerm PROTO((DocType type, DocField field, CharPtr term, Int4Ptr spcl, Int4Ptr totl));

/**** Look up terms with Boolean operations ************/

NLM_EXTERN ValNodePtr LIBCALL EntrezTLNew PROTO((DocType type));
NLM_EXTERN ValNodePtr LIBCALL EntrezTLAddTerm PROTO((ValNodePtr elst, CharPtr term, DocType type, DocField field, Boolean special));
NLM_EXTERN ValNodePtr LIBCALL EntrezTLAddTermWithRange PROTO((ValNodePtr elst, CharPtr term, DocType type, DocField field, Boolean special, CharPtr highRange));
NLM_EXTERN ValNodePtr LIBCALL EntrezTLAddLParen PROTO((ValNodePtr elst));
NLM_EXTERN ValNodePtr LIBCALL EntrezTLAddRParen PROTO((ValNodePtr elst));
NLM_EXTERN ValNodePtr LIBCALL EntrezTLAddAND PROTO((ValNodePtr elst));
NLM_EXTERN ValNodePtr LIBCALL EntrezTLAddOR PROTO((ValNodePtr elst));
NLM_EXTERN ValNodePtr LIBCALL EntrezTLAddBUTNOT PROTO((ValNodePtr elst));
NLM_EXTERN ValNodePtr LIBCALL EntrezTLFree PROTO((ValNodePtr elst));
NLM_EXTERN LinkSetPtr LIBCALL EntrezPMTLEval PROTO((ValNodePtr elst, void * edc));
NLM_EXTERN Int4 LIBCALL EntrezPMTLEvalCount PROTO((ValNodePtr elst, void * edc));
NLM_EXTERN ByteStorePtr LIBCALL EntrezPMTLEvalX PROTO((ValNodePtr elst, void * edc));
NLM_EXTERN LinkSetPtr LIBCALL EntrezTLEval PROTO((ValNodePtr elst));
NLM_EXTERN Int4 LIBCALL EntrezTLEvalCount PROTO((ValNodePtr elst));
NLM_EXTERN ByteStorePtr LIBCALL EntrezTLEvalX PROTO((ValNodePtr elst));


/**** Look Up a Uid from a SeqId using the Terms list ****/

NLM_EXTERN Int4 LIBCALL EntrezFindSeqId PROTO((SeqIdPtr sip));

/**** Look Up the source SeqId given a GI ****************/

NLM_EXTERN SeqIdPtr LIBCALL EntrezSeqIdForGI PROTO ((Int4 gi));

/**** Look Up the Fasta entry given a GI or Accession ****/

NLM_EXTERN FastaSeqPtr LIBCALL AccessionToFasta PROTO ((CharPtr string));

/**** Get Sequence or MEDLINE data **********************/

NLM_EXTERN Int2 LIBCALL EntrezSeqEntryListGet PROTO((SeqEntryPtr PNTR result, Int2 numuid, Int4Ptr uids, Int2 retcode, Boolean mark_missing));
NLM_EXTERN SeqEntryPtr LIBCALL EntrezSeqEntryGet PROTO((Int4 uid, Int2 retcode));

NLM_EXTERN Int2 LIBCALL EntrezMedlineEntryListGet PROTO((MedlineEntryPtr PNTR result, Int2 numuid, Int4Ptr uids, Boolean mark_missing));
NLM_EXTERN MedlineEntryPtr LIBCALL EntrezMedlineEntryGet PROTO((Int4 uid));

#ifdef _PMENTREZ_
Int2 LIBCALL EntrezPubmedEntryListGet PROTO((PubmedEntryPtr PNTR result, Int2 numuid, Int4Ptr uids, Boolean mark_missing));
PubmedEntryPtr LIBCALL EntrezPubmedEntryGet PROTO((Int4 uid));
#endif

NLM_EXTERN BiostrucPtr LIBCALL EntrezBiostrucGet PROTO((DocUid uid, Int4 mdlLvl, Int4 maxModels));
NLM_EXTERN BiostrucAnnotSetPtr LIBCALL EntrezBiostrucAnnotSetGet PROTO((DocUid uid));
NLM_EXTERN LinkSetPtr LIBCALL EntrezBiostrucFeatIds PROTO((DocUid mmdbid, Int2 feature_type, Int4 feature_set_id));
NLM_EXTERN BiostrucAnnotSetPtr LIBCALL EntrezBiostrucAnnotSetGetByFid PROTO((DocUid mmdbid, Int4 feature_id, Int4 feature_set_id));


NLM_EXTERN EntrezHierarchyPtr LIBCALL EntrezHierarchyGet PROTO((CharPtr term,DocType type, DocField fld));
NLM_EXTERN EntrezHierarchyPtr LIBCALL EntrezHierarchyFree PROTO ((EntrezHierarchyPtr ehp));

/***** Neighbor text on-the-fly ***********/
NLM_EXTERN EntrezNeighborTextPtr LIBCALL EntrezNeighborTextFree PROTO ((EntrezNeighborTextPtr ));
NLM_EXTERN EntrezNeighborTextPtr LIBCALL EntrezNeighborTextNew PROTO (( void ));
NLM_EXTERN LinkSetPtr LIBCALL EntrezDoNeighborText PROTO((EntrezNeighborTextPtr entp));
NLM_EXTERN Boolean LIBCALL EntrezCanNeighborText PROTO((void));


NLM_EXTERN Boolean LIBCALL EntrezExpandedMedlineFeatures PROTO((void));
NLM_EXTERN Int4 LIBCALL EntrezClusterAnalysis PROTO((DocUidPtr uids, Int4 numuids, DocField fld, Int4 minCluster, Int4 maxCluster, Int4 maxTerms, CharPtr *terms, Int4Ptr termTotals));

/***** on-the-fly BLAST *****/
NLM_EXTERN LinkSetPtr LIBCALL EntrezBlastBioseq PROTO((BioseqPtr bsp, DocType db, CharPtr program, CharPtr database, CharPtr options, Boolean usemonitor));
NLM_EXTERN Boolean LIBCALL EntrezCanBlast PROTO((void));

/*****************************************************************************
*
*   The Following two functions allow access by BioseqFetch using the
*   SeqMgr.  The application should call EntrezBioseqFetchEnable() at the start
*   of the application and EntrezBioseqFetchDisable() at the end; This
*   will make EntrezBioseqFetch() the "remote" access procedure for the
*   SeqMgr. EntrezInit() will only be called on the first fetch unless "now"
*   is true;
*
*   If you add your own fetch function after calling EntrezBioseqFetchEnable,
*     it will be called BEFORE EntrezBioseqFetchEnable. Add yours after this
*     call, and yours will be call AFTER entrez.
*   
*****************************************************************************/
NLM_EXTERN Boolean LIBCALL EntrezBioseqFetchEnable PROTO((CharPtr progname, Boolean now));
NLM_EXTERN void LIBCALL EntrezBioseqFetchDisable PROTO((void));

#ifdef __cplusplus
}
#endif

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif


