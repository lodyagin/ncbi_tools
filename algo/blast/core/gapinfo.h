/* $Id: gapinfo.h,v 1.13 2004/07/02 19:23:38 dondosha Exp $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's offical duties as a United States Government employee and
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
 * Author: Ilya Dondoshansky
 *
 */

/** @file gapinfo.h
 * Structures definitions from gapxdrop.h in ncbitools 
 * @todo FIXME: doxygen comments
 */

#ifndef __GAPINFO__
#define __GAPINFO__

#include <algo/blast/core/blast_def.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Operation types within the edit script*/
typedef enum EGapAlignOpType { 
   eGapAlignDel = 0, /**< Deletion: a gap in query */
   eGapAlignDel2 = 1,/**< Frame shift deletion of two nucleotides */
   eGapAlignDel1 = 2,/**< Frame shift deletion of one nucleotide */
   eGapAlignSub = 3, /**< Substitution */
   eGapAlignIns1 = 4,/**< Frame shift insertion of one nucleotide */
   eGapAlignIns2 = 5,/**< Frame shift insertion of two nucleotides */
   eGapAlignIns = 6, /**< Insertion: a gap in subject */
   eGapAlignDecline = 7 /**< Non-aligned region */
} EGapAlignOpType;

/** Edit script: linked list of correspondencies between two sequences */
typedef struct GapEditScript {
   EGapAlignOpType op_type;    /**< Type of operation */
   Int4 num;                   /**< Number of operations */
   struct GapEditScript* next; /**< Pointer to next link */
} GapEditScript;

/** Editing block structure containing all information returned from a 
 * single gapped extension. */
typedef struct GapEditBlock {
   Int4 start1,  start2,       /**< starts of alignments. */
      length1, length2,       /**< total lengths of the sequences. */
      original_length1, original_length2; /**< Untranslated lengths of the 
                                             sequences. */
   Int2 frame1, frame2;	    /**< frames of the sequences. */
   Boolean translate1,	translate2; /**< are either of these be translated. */
   Boolean reverse;	/**< reverse sequence 1 and 2 when producing SeqALign? */
   Boolean is_ooframe; /**< Is this out_of_frame edit block? */
   Boolean discontinuous; /**< Is this OK to produce discontinuous SeqAlign? */
   GapEditScript* esp; /**< Editing script for the traceback. */
} GapEditBlock;

/** Structure to keep memory for state structure. */ 
typedef struct GapStateArrayStruct {
	Int4 	length,		/**< length of the state_array. */
		used;		/**< how much of length is used. */
	Uint1* state_array;	/**< array to be used. */
	struct GapStateArrayStruct* next; /**< Next link in the list. */
} GapStateArrayStruct;

/** Initialize the gap editing script structure. */
GapEditScript* 
GapEditScriptNew (GapEditScript* old);

/** Free the gap editing script structure. */
GapEditScript* GapEditScriptDelete (GapEditScript* esp);

/** Initialize the gap editing block structure. 
 * @param start1 Offset to start alignment in first sequence. [in]
 * @param start2 Offset to start alignment in second sequence. [in]
 */
GapEditBlock* GapEditBlockNew (Int4 start1, Int4 start2);

/** Free the edit block structure. */
GapEditBlock* GapEditBlockDelete (GapEditBlock* edit_block);

/** Free the gap state structure. */
GapStateArrayStruct* 
GapStateFree(GapStateArrayStruct* state_struct);

#ifdef __cplusplus
}
#endif
#endif /* !__GAPINFO__ */
