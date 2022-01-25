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
* $Revision: 6.15 $
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
* Revision 6.15  2000/03/03 18:00:11  dondosha
* Added prototype for MegaBlastLookupTableDup, plus cosmetic changes
*
* Revision 6.14  2000/02/17 19:00:45  shavirin
* Removed theCacheSize parameter from everywhere.
*
* Revision 6.13  2000/02/11 20:52:25  dondosha
* Added Webb Miller s lookup table and related prototypes
*
* Revision 6.12  2000/02/01 21:05:56  dondosha
* Added prototype for mb_lookup_position_aux_destruct - modification for megablast
*
* Revision 6.11  2000/01/11 17:27:02  shavirin
* Added elements theTable and theCacheSize into LookupTable
* Added parameters theCacheSize into lookup_new and BLAST_WordFinderNew
*
* Revision 6.9  1999/12/21 21:33:24  shavirin
* Added parameter mod_lookup_table_memory to the lookup structure.
*
* Revision 6.8  1999/09/29 17:20:12  madden
* Deallocate ModLookupPositionPtr memory
*
* Revision 6.7  1999/09/28 20:14:34  madden
* Joerg changes to mimize cache misses
*
* Revision 6.6  1999/09/16 16:56:39  madden
* Changes to lookup_new for long words
*
* Revision 6.5  1999/09/16 14:16:27  madden
* lookup_find_init returns Uint1Ptr instead of CharPtr
*
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



typedef struct orig_lookup_position {
	Int4	position;		/* position along the query */
	Int1 context;	/* "context" of hit, for use by BLASTContextStructPtr.*/
	struct orig_lookup_position *next;	/* other positions */
} OrigLookupPosition, PNTR OrigLookupPositionPtr;

/* some defines for the pv_array, as this changes from 32-bit to 64-bit systems. */
#if LONG_BIT==64

#define PV_ARRAY_TYPE Uint8	/* The pv_array 'native' type. */
#define PV_ARRAY_BYTES 8	/* number of BYTES in 'native' type. */
#define PV_ARRAY_BTS 6		/* bits-to-shift from lookup_index to pv_array index. */
#define PV_ARRAY_MASK 63	/* amount to mask off. */

#else

#define PV_ARRAY_TYPE Uint4	/* The pv_array 'native' type. */
#define PV_ARRAY_BYTES 4	/* number of BYTES in 'native' type. */
#define PV_ARRAY_BTS 5		/* bits-to-shift from lookup_index to pv_array index. */
#define PV_ARRAY_MASK 31	/* amount to mask off. */

#endif

#define PV_ARRAY_FACTOR 0.5	/* The fraction of sites that must have at least one hit to not use PV_ARRAY. */

/* Define LookupPosition (ie a hit_info) such that it can be loaded into a register */
typedef Uint4 ModLookupPosition;
typedef Uint4 *ModLookupPositionPtr;

#define hinfo_get_pos(h) ((h)&0xffffff)
#define hinfo_get_context(h) (((h)>>24)&0xff)

#define hinfo_set(ptr_hinfo,p,c) ( *(ptr_hinfo) = (( ((Uint4)(c)) << 24) | ((Uint4)(p))) )

#if 1                           /* To be deleted */
typedef struct mod_lt_entry {
    Int4 num_used;       /* num valid positions */
    ModLookupPosition entries[3]; /* first postion */
} ModLAEntry, PNTR ModLAEntryPtr;
#endif

#define MBSTACK_SIZE 5000

typedef struct megablast_stack {
    int diag, level, length;
} MbStack, PNTR MbStackPtr;

typedef struct megablast_lookup_table {
   Int2 width;          /* Number of bytes in hash value   */
   Int4 hashsize;      /* = 2^(8*width) - 1               */ 
   Int2 lpm;            /* Minimal length of perfect match */
   Int2 max_positions;  /* Maximal number of positions for one hash value */
   Int4Ptr hashtable;   /* Array of positions              */
   Int4Ptr next_pos;    /* Extra positions stored here     */
   Uint4 mask; 
   Int4 stack_index;    
   MbStackPtr estack; /* Most recent hits stored here */
} MbLookupTable, PNTR MbLookupTablePtr;

typedef struct lookup_table {
    Int4 char_size,		/* number of bits per residue/bp */
        wordsize,		/* size of "word" */
        reduced_wordsize,       /* size of word */
        array_size,		/* size of table's array. */
        mask;		/* Used to mask off top set of bits. */
    OrigLookupPositionPtr PNTR position;	/* positions of the hits. */
    OrigLookupPositionPtr PNTR position_aux;	/* auxillary structure for keeping track of the last saved hit, to speed up saving of hits on very long sequences. */
    LookupMemoryPtr mem_struct,	/* contains memory. */
        mem_struct_start;	/* Start of LookupMemoryPtr chain. */
    size_t memory_chunk;	/* chunk size of memory. */
    PV_ARRAY_TYPE *pv_array;                /* presence vector, used for quick presence check */
    Int4 num_pos_added;
    Int4 num_unique_pos_added;
#if 1                           /* To be deleted */
    ModLAEntry *mod_lt;	        /* The new&improved lookup table */
#endif
    Uint4Ptr theTable;          /* Lookup table with variable length */
    ModLookupPositionPtr 	mod_lookup_table_memory; /* Memory for new and improved lookup table. */
    Int4 mod_lookup_table_size;/* Size for new and improved lookup table */
   MbLookupTablePtr mb_lt; /* MegaBlast version of the lookup table */
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


LookupTablePtr LIBCALL lookup_new PROTO((Int2 alphabet_size, Int2 wordsize, Int2 reduced_wordsize));

LookupTablePtr LIBCALL lookup_destruct PROTO((LookupTablePtr lookup));

void LIBCALL lookup_add PROTO((LookupTablePtr lookup, CharPtr string, Int4 position, Int1 context));

void LIBCALL lookup_add_index PROTO((LookupTablePtr lookup, Int4 lookup_index, Int4 position, Int1 context));

Uint1Ptr LIBCALL lookup_find_init PROTO((LookupTablePtr lookup, Int4 PNTR lookup_index, Uint1Ptr string));

Boolean lookup_position_aux_destruct PROTO((LookupTablePtr lookup));
Boolean mb_lookup_position_aux_destruct PROTO((LookupTablePtr lookup));

LookupTablePtr MegaBlastLookupTableDup PROTO((LookupTablePtr lookup));
LookupTablePtr MegaBlastLookupTableDestruct PROTO((LookupTablePtr lookup));

BLAST_WordFinderPtr BLAST_WordFinderDestruct PROTO((BLAST_WordFinderPtr wfp));
BLAST_WordFinderPtr BLAST_WordFinderNew PROTO((Int2 alphabet_size, Int2 wordsize, Int2 compression_ratio, Boolean round_down));



#ifdef __cplusplus
}
#endif

#endif /* _LOOKUP_ */
