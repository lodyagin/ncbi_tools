/*   asn2gnbk.h
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
* File Name:  asn2gnbk.h
*
* Author:  Karl Sirotkin, Tom Madden, Tatiana Tatusov, Jonathan Kans
*
* Version Creation Date:   10/21/98
*
* $Revision: 6.3 $
*
* File Description:  New GenBank flatfile generator
*
* Modifications:  
* --------------------------------------------------------------------------
* ==========================================================================
*/

#ifndef _ASN2NGNBK_
#define _ASN2NGNBK_

#include <ncbi.h>
#include <objall.h>
#include <seqport.h>

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define GENBANK_FMT   ((Uint1) 0)
#define EMBL_FMT      ((Uint1) 1)
#define GENPEPT_FMT   ((Uint1) 2)
#define EMBLPEPT_FMT  ((Uint1) 5)

#define RELEASE_MODE  ((Uint1)  6)
#define SEQUIN_MODE   ((Uint1)  7)
#define DUMP_MODE     ((Uint1)  8)

#define SEGMENTED_STYLE   ((Uint1) 10)
#define MASTER_STYLE      ((Uint1) 11)

#define SEQUENCE_STYLE     ((Uint1) 15)
#define CONTIG_STYLE       ((Uint1) 16)

#define LOCUS_BLOCK       ((Uint1)  1)
#define DEFLINE_BLOCK     ((Uint1)  2)
#define ACCESSION_BLOCK   ((Uint1)  3)
#define VERSION_BLOCK     ((Uint1)  4)
#define PID_BLOCK         ((Uint1)  5)
#define DBSOURCE_BLOCK    ((Uint1)  6)
#define DATE_BLOCK        ((Uint1)  7)
#define KEYWORDS_BLOCK    ((Uint1)  8)
#define SEGMENT_BLOCK     ((Uint1)  9)
#define ORGANISM_BLOCK    ((Uint1) 10)
#define REFERENCE_BLOCK   ((Uint1) 11)
#define COMMENT_BLOCK     ((Uint1) 12)
#define FEATHEADER_BLOCK  ((Uint1) 13)
#define SOURCE_BLOCK      ((Uint1) 14)
#define FEATURE_BLOCK     ((Uint1) 15)
#define BASECOUNT_BLOCK   ((Uint1) 16)
#define ORIGIN_BLOCK      ((Uint1) 17)
#define SEQUENCE_BLOCK    ((Uint1) 18)
#define CONTIG_BLOCK      ((Uint1) 19)
#define SLASH_BLOCK       ((Uint1) 20)

#define ASN2GB_BASE_BLOCK \
  Uint2             entityID;  \
  Uint2             itemID;    \
  Uint2             itemtype;  \
  Int2              section;   \
  Int2              blocktype; \
  CharPtr           string;    \

/* base block structure for most paragraph types */

typedef struct asn2gb_base_block {
  ASN2GB_BASE_BLOCK
} BaseBlock, PNTR BaseBlockPtr;

/* version block includes VERSION and NID sections */
/* organism block includes SOURCE and ORGANISM sections */
/* source (feat) block should be the same as the organism block */

/* references are grouped by published, unpublished, sites, and cit-subs */

#define REF_CAT_PUB 1
#define REF_CAT_UNP 2
#define REF_CAT_SIT 3
#define REF_CAT_SUB 4

typedef struct reference_block {
  ASN2GB_BASE_BLOCK
  Int4              pmid;
  Int4              muid;
  CharPtr           uniquestr;
  Int2              serial;
  Int2              category;
} ReferenceBlock, PNTR ReferenceBlockPtr;

/* sequences are broken up into paragraphs and use the section's SeqPort */

typedef struct sequence_block {
  ASN2GB_BASE_BLOCK
  Int4              start;
  Int4              stop;
} SequenceBlock, PNTR SequenceBlockPtr;


/* structure for single segment or pop/phy/mut set component */

typedef struct asn2gbsection {

  /* data identifiers for individual accession report */

  BioseqPtr          parent;
  BioseqPtr          bsp;
  SeqLocPtr          slp;
  Uint2              seg;
  Int2               numsegs;
  Int4               from;
  Int4               to;

  /* SeqPort for section's sequence */

  SeqPortPtr         spp;

  /* local array pointing to all blocks in this section */

  BaseBlockPtr       PNTR blockArray;
  Int2               numBlocks;

  /* referenceks for feature citation matching, serial number assignment */

  ReferenceBlockPtr  PNTR referenceArray;
  Int2               numReferences;

} Asn2gbSection, PNTR Asn2gbSectionPtr;


/* master pointer returned to application */

typedef struct asn2gb_job {

  /* data identifiers for sequence or sequences to report */

  Uint2               entityID;
  BioseqPtr           bsp;
  BioseqSetPtr        bssp;
  SeqLocPtr           slp;

  /* flags for customizing type of report */

  Int2                format;
  Int2                mode;

  /* each accession report from LOCUS to // is a single section */

  Asn2gbSectionPtr    PNTR sectionArray;
  Int2                numSections;

  /* master array pointing to all blocks in all sections */

  BaseBlockPtr        PNTR paragraphArray;
  Int2                numParagraphs;

} Asn2gbJob, PNTR Asn2gbJobPtr;



/*
  asn2gnbk_setup creates a structure laying out the flatfile structure.

    Of the first three parameters (bsp, bssp, slp), only one should not be NULL.
    If bsp is passed in, that is the target bioseq.  (If it is segmented, then
    the report currently shows multiple segments.)  If bssp is passed in, it is
    expected to be a population/phylogenetic/mutation study, and all records are
    shown, but with no segment numbers.  slp is not yet implemented.

  asn2gnbk_format creates a string for a given paragraph.

  asn2gnbk_cleanup frees the structure and all components.
*/

NLM_EXTERN Asn2gbJobPtr asn2gnbk_setup (BioseqPtr bsp, BioseqSetPtr bssp,
                                        SeqLocPtr slp, Uint1 format, Uint1 mode,
                                        Uint1 segstyle, Uint1 seqstyle);

NLM_EXTERN CharPtr asn2gnbk_format (Asn2gbJobPtr ajp, Int2 paragraph);

NLM_EXTERN Asn2gbJobPtr asn2gnbk_cleanup (Asn2gbJobPtr ajp);

/*
   SeqEntryToGnbk calls asn2gnbk_setup, _format, and _cleanup internally.
*/

NLM_EXTERN Boolean SeqEntryToGnbk (SeqEntryPtr sep, Uint1 format, Uint1 mode,
                                   Uint1 segstyle, Uint1 seqstyle, FILE *fp);


#ifdef __cplusplus
}
#endif

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif /* ndef _ASN2NGNBK_ */

