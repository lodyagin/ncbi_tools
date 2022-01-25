/*
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
*/
/*****************************************************************************

File name: lookup.h

Author: Tom Madden

Contents: defines and prototype used by lookup.c.

******************************************************************************/

/* File Name: lookup.h
*
* Author: Tom Madden
*
* Version Creation Date:   10/26/95
*
* $Revision: 6.4 $
*
* File Description: 
*       Functions that format traditional BLAST output.
*
* Modifications:  
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
* ==========================================================================
*
*
* RCS Modification History:
* $Log: lookup.h,v $
* Revision 6.4  1998/09/22 16:27:48  madden
* Added function lookup_position_aux_destruct
*
* Revision 6.3  1998/07/27 18:14:18  madden
* lookup_get_memory replaces call to MemNew
*
* Revision 6.2  1998/04/15 20:25:58  madden
* Auxillary structure added to speed-up saving words for queries
*
* Revision 6.1  1998/02/26 22:34:40  madden
* Changes for 16 bit windows
*
* Revision 6.0  1997/08/25 18:53:25  madden
* Revision changed to 6.0
*
* Revision 1.4  1996/09/25 14:17:30  madden
* removed discontiguous options.
*
 * Revision 1.3  1996/09/12  21:12:59  madden
 * Added new function to save an already computed index, lookup_add_index.
 *
 * Revision 1.2  1996/08/15  18:57:10  madden
 * Changed context from Int1 to Int2.
 *
 * Revision 1.1  1996/08/05  19:48:50  madden
 * Initial revision
 *
 * Revision 1.10  1996/07/18  22:00:02  madden
 * Changes for multiple contexts.
 *
 * Revision 1.9  1996/06/20  17:00:11  madden
 * Added "__cplusplus" define.
 *
 * Revision 1.8  1996/06/20  16:15:57  madden
 * Replaced int's with Int4's.
 *
 * Revision 1.7  1996/05/22  20:22:43  madden
 * removed unused prototypes.
 *
 * Revision 1.6  1996/05/16  19:50:15  madden
 * Added documentation block.
 *
 * Revision 1.5  1996/05/03  19:55:24  madden
 * Removed FnPtr for lookup_find.
 *
 * Revision 1.4  1996/04/04  20:47:55  madden
 * Made lookup_find into two functions, that are statics.
 *
 * Revision 1.3  1996/02/28  21:38:36  madden
 * Changed prototypes for discontiguous words.
 *
 * Revision 1.2  1995/12/26  20:28:54  madden
 * *** empty log message ***
 *
 * Revision 1.1  1995/12/08  15:48:23  madden
 * Initial revision
 *
*/

#ifndef _LOOKUP_
#define _LOOKUP_

#ifdef __cplusplus
extern "C" {
#endif

#include <ncbi.h>

/* 
	Structure that keeps a block of memory for the functions in lookup.c
	Most often used instead of a MemNew call in AddPositionToLookupTable.
*/
typedef struct lookup_memory {
	Uint1Ptr current,	/* Current position of pointer to memory. */
		start;		/* Start of memory (for deallocation). */
	size_t	remaining;	/* How much memory remains in allocated block. */
	struct 	lookup_memory *next;
} LookupMemory, PNTR LookupMemoryPtr;



typedef struct lookup_position {
	Int4	position;		/* position along the query */
	Int1 context;	/* "context" of hit, for use by BLASTContextStructPtr.*/
	struct lookup_position *next;	/* other positions */
} LookupPosition, PNTR LookupPositionPtr;

typedef struct lookup_table {
	Int4	char_size,		/* number of bits per residue/bp */
		wordsize,		/* size of "word" */
		array_size,		/* size of table's array. */
		mask;		/* Used to mask off top set of bits. */
	LookupPositionPtr PNTR position;	/* positions of the hits. */
	LookupPositionPtr PNTR position_aux;	/* auxillary structure for keeping track of 
						the last saved hit, to speed up saving of hits
						on very long sequences. */
	LookupMemoryPtr mem_struct,	/* contains memory. */
			mem_struct_start;	/* Start of LookupMemoryPtr chain. */
	size_t		memory_chunk;	/* chunk size of memory. */
} LookupTable, PNTR LookupTablePtr;


/*********************************************************************
	Structure for the BLAST_WordFinder, used to find the initial 
	word-hits.
*********************************************************************/
typedef struct _blast_wordfinder {
                Int4     wordsize;
		Int2 compression_ratio;
                LookupTablePtr  lookup;
        } BLAST_WordFinder, PNTR BLAST_WordFinderPtr;


LookupTablePtr LIBCALL lookup_new PROTO((Int2 alphabet_size, Int2 wordsize));

LookupTablePtr LIBCALL lookup_destruct PROTO((LookupTablePtr lookup));

void LIBCALL lookup_add PROTO((LookupTablePtr lookup, CharPtr string, Int4 position, Int1 context));

void LIBCALL lookup_add_index PROTO((LookupTablePtr lookup, Int4 lookup_index, Int4 position, Int1 context));

CharPtr LIBCALL lookup_find_init PROTO((LookupTablePtr lookup, Int4 PNTR lookup_index, CharPtr string));

Boolean lookup_position_aux_destruct PROTO((LookupTablePtr lookup));

BLAST_WordFinderPtr BLAST_WordFinderDestruct PROTO((BLAST_WordFinderPtr wfp));
BLAST_WordFinderPtr BLAST_WordFinderNew PROTO((Int2 alphabet_size, Int2 wordsize, Int2 compression_ratio, Boolean round_down));


#ifdef __cplusplus
}
#endif

#endif /* _LOOKUP_ */
