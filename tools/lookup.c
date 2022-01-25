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

File name: lookup.c

Author: Tom Madden

Contents: Functions to store "hashed" word hits and their positions and
	in arrays and look these up again.

Detailed Contents:

        - the number of bits required is given by ln(alphabet size)/ln(2).

	- the number of bits times the wordsize must be less than 32.

        - contiguous words can be of any length, discontiguous words can
	only have the patten XX0X and be of length three, or word-width of
	four.

******************************************************************************/

/*******************************************************************************
* File Name: lookup.c
*
* Author: Tom Madden
*
* Version Creation Date:   10/26/95
*
* $Revision: 6.20 $
*
* File Description: 
*       Functions to store "words" from a query and perform lookups against
*	database sequences. 
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
* $Log: lookup.c,v $
* Revision 6.20  2000/01/11 18:37:08  shavirin
* Added define DYNAMIC_CACHE to distinduish dynamic and static lookup.
*
* Revision 6.19  2000/01/11 17:58:50  shavirin
* Added new function make_mod_lt_new, which created dynamic lookup table
*
* Revision 6.18  2000/01/11 17:33:11  shavirin
* Added handling of the parameter theCacheSize.
*
* Revision 6.17  1999/12/21 21:33:23  shavirin
* Added parameter mod_lookup_table_memory to the lookup structure.
*
* Revision 6.15  1999/11/08 21:54:54  madden
* Fix typo
*
* Revision 6.14  1999/10/12 18:00:13  madden
* Add one to array lengths for purify errors.
*
* Revision 6.13  1999/10/06 14:25:20  madden
* Only allocate pv_array if needed
*
* Revision 6.12  1999/09/29 17:20:12  madden
* Deallocate ModLookupPositionPtr memory
*
* Revision 6.11  1999/09/29 13:28:44  madden
* Fill in pv_array correctly
*
* Revision 6.10  1999/09/28 20:14:33  madden
* Joerg changes to mimize cache misses
*
* Revision 6.9  1999/09/16 16:56:39  madden
* Changes to lookup_new for long words
*
* Revision 6.8  1999/09/16 14:16:27  madden
* lookup_find_init returns Uint1Ptr instead of CharPtr
*
* Revision 6.7  1999/04/27 15:40:12  madden
* Use aux lookup only for short words
*
* Revision 6.6  1998/09/22 16:27:47  madden
* Added function lookup_position_aux_destruct
*
* Revision 6.5  1998/07/27 18:14:17  madden
* lookup_get_memory replaces call to MemNew
*
* Revision 6.4  1998/04/15 20:25:58  madden
* Auxillary structure added to speed-up saving words for queries
*
* Revision 6.3  1998/02/26 22:34:39  madden
* Changes for 16 bit windows
*
* Revision 6.2  1998/01/05 20:32:38  kans
* AddPositionToLookupTable checks for NULL lookup parameter
*
* Revision 6.1  1998/01/05 17:38:05  madden
* Check that position is non-NULL in lookup_destruct
*
* Revision 6.0  1997/08/25 18:53:23  madden
* Revision changed to 6.0
*
* Revision 1.7  1997/08/19 18:19:23  madden
* Cast arg of log to Nlm_FloatHi
*
* Revision 1.6  1997/04/03 19:58:27  madden
* Added check for NULL lookup->position.
*
 * Revision 1.5  1996/09/25  14:17:30  madden
 * removed discontiguous options.
 *
 * Revision 1.4  1996/09/12  21:12:59  madden
 * Added new function to save an already computed index, lookup_add_index.
 *
 * Revision 1.3  1996/08/21  21:28:31  madden
 * Added casts to quiet NT compiler warning.s
 *
 * Revision 1.2  1996/08/15  18:33:35  madden
 * Changed Int2 to Int1.
 *
 * Revision 1.1  1996/08/05  19:48:50  madden
 * Initial revision
 *
 * Revision 1.11  1996/07/24  12:01:28  madden
 * Changes for blastx
 *
 * Revision 1.10  1996/07/18  22:00:02  madden
 * Changes for multiple contexts.
 *
 * Revision 1.9  1996/06/20  16:52:23  madden
 * Changed "pow" to "Nlm_Powi".
 *
 * Revision 1.8  1996/06/20  16:15:57  madden
 * Replaced int's with Int4's.
 *
 * Revision 1.7  1996/05/22  20:22:01  madden
 * Removed unused variable lookup_find_discontig.
 *
 * Revision 1.6  1996/05/16  19:50:15  madden
 * Added documentation block.
 *
 * Revision 1.5  1996/05/03  19:55:24  madden
 * Combined two functions.
 *
 * Revision 1.4  1996/04/04  20:47:55  madden
 * Made lookup_find into two functions, that are statics.
 *
 * Revision 1.3  1996/02/28  21:37:14  madden
 * changes for discontiguous words.
 *
 * Revision 1.2  1995/12/26  20:28:36  madden
 * removed lookup_add_word function.
 *
 * Revision 1.1  1995/12/08  15:48:23  madden
 * Initial revision
 *
*/

#include <ncbi.h>
#include <lookup.h>

/* mod_lt[] is a modified lookup table.  -cfj
 * It contains the same info that is accessible via lookup->position, but in a form that reduces the
 * number of cache misses.
 * Doing a lookup via the original lookup table cost one cache miss per access, and then one miss per 
 *   entry as we walk down the hit list.
 *
 * mod_lt replaces the linked list with an array of hits, to reduce the number of misses.
 * Also, if there are 3 or fewer hits, they are placed directly in the table.
 *       if there are more than 3 hits, the first one is placed directly in the table (followed by
 *       a pointer to a contiguous array of all the remaining hits).
 *
 * This should probably be allocated and deallocated along with the standard lookup table, but this 
 * is not currently done.
 *
 * The make_mod_lt function should be called to build the mod_lt table after all hits have been placed
 * into the orginal lookup table.
 */

#ifdef DYNAMIC_CACHE
void make_mod_lt(LookupTablePtr lookup)
#else
void make_mod_lt_new(LookupTablePtr lookup)
#endif
{
    Int4 index, j, pv_size;
    Uint4 len;
    OrigLookupPositionPtr PNTR lookup_pos;
    OrigLookupPositionPtr list;
    ModLookupPositionPtr next_free;
    PV_ARRAY_TYPE *pv_array=NULL;
    Uint4Ptr theTable;
    Int4 theCacheSize, theArraySize;

    theTable = lookup->theTable;
    theCacheSize = lookup->theCacheSize;
    theArraySize = theCacheSize + 1;
    
    if(lookup->theCacheSize < 3) {
        ErrPostEx(SEV_ERROR, 0, 0, "Invalid Cache size for lookup table %d", 
                  lookup->theCacheSize);
        return;
    }
    
    len = lookup->array_size;
    lookup_pos = lookup->position;
    
    /* Add 1 to silence purify messages. */
    /* Initial allocation - more, than necessary though */
    
    next_free = MemNew((1+lookup->num_pos_added)*sizeof(ModLookupPosition));
    lookup->mod_lookup_table_size = 0;
    
    if (next_free == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Unable to allocate lookup->position");
        return;
    }
    lookup->mod_lookup_table_memory = next_free;
    
    
    /* pv_array is an array of 'presence bits', one bit per entry 
       in the lookup table, packed 64bits per Uint8 (on a 64-bit binary); 
       otherwise a Uint4 is used. It can be used to quickly see if there 
       is an entry in the lookup table. It's used on short-to-medium 
       length queries since most lookups will be empty. -cfj
    */
    
    if (lookup->num_unique_pos_added < PV_ARRAY_FACTOR*lookup->array_size) {
	pv_size = (lookup->array_size+PV_ARRAY_MASK)/(PV_ARRAY_MASK+1);    /* Size measured in PV_ARRAY_TYPE's */
	pv_array = MemNew(pv_size*PV_ARRAY_BYTES);
    }

    /* ---------
       Elements in the Table references as:
       
       *(theTable + theArraySize*index) == .num_used
       
       if(.num_used <= theCasheSize)
       *(theTable + theArraySize*index + [1...theCasheSize]) == .entries
       else
       *(theTable + theArraySize*index + [1...theCasheSize-2]) == .entries
       *(theTable + theArraySize*index + [1...theCasheSize-1]) &
       *(theTable + theArraySize*index + [1...theCasheSize]) == address of
       extended memory chunk
       ----------- */
    
    /* Walk through table, copying info into mod_lt[] */
    for(index = 0; index < len; index++) {
        if(lookup_pos[index] == NULL) {
            *(theTable + theArraySize*index) = 0;    
        } else {
            /* count num entries */
            int count=0;
            list = lookup_pos[index];
            
            if (pv_array)
                pv_array[(index>>PV_ARRAY_BTS)] |= (((PV_ARRAY_TYPE) 1)<<(index&PV_ARRAY_MASK));
            
            while (list != NULL) {
                count++;
                list = list->next;
            }
            
            if(count <= theCacheSize) {
                list = lookup_pos[index];
                *(theTable + theArraySize*index) = count;
                for(j = 0;j < count; j++){
                    hinfo_set(theTable + theArraySize*index + j + 1, list->position, list->context);
                    list = list->next;
                }	    
            } else {
                ModLookupPositionPtr * lpp= (ModLookupPositionPtr *)(theTable + theArraySize*index + theCacheSize - 1);
                list=lookup_pos[index];
                *(theTable + theArraySize*index) = count;
                for(j = 0;j < theCacheSize-2;j++) {
                    
                    hinfo_set(theTable + theArraySize*index + j + 1, list->position, list->context);
                    list = list->next;
                    *lpp  = next_free;
                }
                for(j = theCacheSize-2; j < count; j++){
                    hinfo_set(next_free,list->position,list->context);
                    next_free++;
                    lookup->mod_lookup_table_size++;
                    list = list->next;
                }
            }
        }
    }
    
    lookup->pv_array = pv_array;

    return;
}
#ifdef DYNAMIC_CACHE
void make_mod_lt_old(LookupTablePtr lookup)
#else
void make_mod_lt(LookupTablePtr lookup)
#endif
{
    Int4 index, j, pv_size;
    Uint4 len;
    OrigLookupPositionPtr PNTR lookup_pos;
    OrigLookupPositionPtr list;
    ModLookupPositionPtr next_free;
    PV_ARRAY_TYPE *pv_array=NULL;
    
    
    ModLAEntry *mod_lt;
    mod_lt=lookup->mod_lt;
    
    len = lookup->array_size;
    lookup_pos = lookup->position;
    
    /* Add 1 to silence purify messages. */
    /* Initial allocation - more, than necessary though */

    next_free = MemNew((1+lookup->num_pos_added)*sizeof(ModLookupPosition));
    lookup->mod_lookup_table_size = 0;
    
    if (next_free == NULL){
        ErrPostEx(SEV_ERROR, 0, 0, "Unable to allocate lookup->position");
        return;
    }
    lookup->mod_lookup_table_memory = next_free;
    
    
    /* pv_array is an array of 'presence bits', one bit per entry 
       in the lookup table, packed 64bits per Uint8 (on a 64-bit binary); 
       otherwise a Uint4 is used. It can be used to quickly see if there 
       is an entry in the lookup table. It's used on short-to-medium 
       length queries since most lookups will be empty. -cfj
    */
    
    if (lookup->num_unique_pos_added < PV_ARRAY_FACTOR*lookup->array_size) {
	pv_size = (lookup->array_size+PV_ARRAY_MASK)/(PV_ARRAY_MASK+1);    /* Size measured in PV_ARRAY_TYPE's */
	pv_array = MemNew(pv_size*PV_ARRAY_BYTES);
    }
    
    /* Walk through table, copying info into mod_lt[] */
    for(index = 0; index < len; index++) {
        if(lookup_pos[index] == NULL){
            mod_lt[index].num_used = 0;
        } else {
            /* count num entries */
            int count=0;
            list = lookup_pos[index];
            
            if (pv_array)
                pv_array[(index>>PV_ARRAY_BTS)] |= (((PV_ARRAY_TYPE) 1)<<(index&PV_ARRAY_MASK));
            
            while (list != NULL){
                count++;
                list = list->next;
            }
            
            if (count<=3) {
                list = lookup_pos[index];
                mod_lt[index].num_used=count;
                for(j=0;j<count;j++){
                    hinfo_set(&mod_lt[index].entries[j],list->position,list->context);
                    list = list->next;
                }	    
            } else {
                ModLookupPositionPtr * lpp= (ModLookupPositionPtr *) &mod_lt[index].entries[1];
                list=lookup_pos[index];
                mod_lt[index].num_used=count;
                hinfo_set(&mod_lt[index].entries[0],list->position,list->context);
                list = list->next;
                *lpp  = next_free;
                for(j = 1; j < count; j++){
                    hinfo_set(next_free,list->position,list->context);
                    next_free++;
                    lookup->mod_lookup_table_size++;
                    list = list->next;
                }
            }
        }
    }

    lookup->pv_array = pv_array;
}


static void AddPositionToLookupTable PROTO((LookupTablePtr lookup, Int4 index, Int4 position, Int1 context));

#define LOOKUP_MEMORY_MIN 128*1024	/* Minimum block of memory (in bytes) for LookupMemoryPtr. */

static VoidPtr lookup_deallocate_memory (LookupTablePtr lookup)
{
    LookupMemoryPtr next=NULL, mem_struct;
    
    if (lookup == NULL)
        return NULL;
    
    mem_struct = lookup->mem_struct_start;
    
    while (mem_struct) {
        next = mem_struct->next;
        MemFree(mem_struct->start);
        MemFree(mem_struct);
        mem_struct = next;
    }

    return NULL;
}

/*
  returns the requested amount of memory.  Large blocks are obtained at
  once and then passed out. 
  
  LookupTablePtr lookup:	contains all information on Lookup.
  size_t required: how much memory has been requested.
*/

static VoidPtr lookup_get_memory (LookupTablePtr lookup, size_t required)
{
    LookupMemoryPtr last, mem_struct;
    size_t memory_requested;
    Uint1Ptr new;
    
    mem_struct = lookup->mem_struct;
    
    if (!mem_struct || required > mem_struct->remaining) {
        if (mem_struct) {
            last = mem_struct;
        } else {
            last = NULL;
        }
        
        memory_requested = MAX(lookup->memory_chunk, LOOKUP_MEMORY_MIN);
        mem_struct = (LookupMemoryPtr) MemNew(sizeof(LookupMemory));
        mem_struct->start = (Uint1Ptr) MemNew(memory_requested);
        mem_struct->current = mem_struct->start;
        mem_struct->remaining = memory_requested;
        
        /* align to 64Byte block. -cfj */
        if(1){
            Uint8 extra = ((Uint8)mem_struct->current) & ((Uint8) 63);
            if (extra>0){
                extra=64-extra;
                mem_struct->current += extra;
                mem_struct->remaining -= extra;
            }
        }
        
        if (last) {
            last->next = mem_struct;	/* used at end for deallocation. */
            lookup->mem_struct = mem_struct; /* Used as a shortcut around a long while loop. */
        } else {
            lookup->mem_struct = mem_struct;
            lookup->mem_struct_start = mem_struct;
        }
    }
    
    new = mem_struct->current;
    mem_struct->current += required;
    mem_struct->remaining -= required;
    
    return (VoidPtr) new;	/* return pointer to new memory. */
}

/*
  LookupTablePtr LIBCALL lookup_new(Int2 alphabet_size, Int2 wordsize)
  
  alphabet_size: number of symbols in alphabet
  
  wordsize: size of initial word hit (often three)
  
  Perform the initial setup of the lookup table.  "lookup" 
  expects to work with a "zero-offset" alphabet such as
  NCBIstdaa or NCBI2na.
  
  Each letter (residue or basepair) in the word goes into it's
  own set of bits (calculated as num_of_bits below).  This is
  to allow masking and shifting of bits.  As an example consider
  a three-letter word with an alphabet containing 20 residues.  
  The minimum num_of_bits is five (able to represent up to 31)
  and the total "array_size" is 15.
  "mask" is used to mask the wordsize-1 lower bits (i.e., to set the
  bits the highest num_of_bits to zero).
*/
LookupTablePtr LIBCALL
lookup_new(Int2 alphabet_size, Int2 wordsize, Int2 reduced_wordsize, Int4 theCacheSize)
     
{
    Int4 num_of_bits;
    LookupTablePtr lookup;
    
    /* How many bits are needed to hold the alphabet? */
    num_of_bits = Nlm_Nint(log((Nlm_FloatHi)alphabet_size)/NCBIMATH_LN2);

    /* 32 bits is 4 bytes */ 
    if (num_of_bits*wordsize > 32) {
        ErrPostEx(SEV_ERROR, 0, 0, "alphabet times wordsize > 32");
        return NULL;
    }
    
    lookup = (LookupTablePtr) MemNew(sizeof(LookupTable));
    if (lookup == NULL)
        return NULL;
    
    lookup->char_size = num_of_bits;
    lookup->wordsize = (Int4) wordsize;
    if (reduced_wordsize <= 0)
        reduced_wordsize = wordsize;
    lookup->reduced_wordsize = (Int4) reduced_wordsize;
    lookup->array_size = (Int4) Nlm_Powi(2.0, (num_of_bits*reduced_wordsize));
    lookup->mask = (Int4) Nlm_Powi(2.0, (num_of_bits*(reduced_wordsize-1))) - 1;
    lookup->mod_lt = NULL;    /* will be allocated when built */
    lookup->theTable = NULL;  /* will be allocated when built */
    lookup->theCacheSize = theCacheSize;
    
    lookup->position = (OrigLookupPositionPtr PNTR) MemNew((lookup->array_size)*sizeof(OrigLookupPositionPtr)); 
    if (lookup->position == NULL) {
        lookup = lookup_destruct(lookup);
        ErrPostEx(SEV_ERROR, 0, 0, "Unable to allocate lookup->position");
    }
    
    /* Only allocate the auxillary structure if it's less than 1 Meg. */
    if ((lookup->array_size)*sizeof(OrigLookupPositionPtr) <= (1<<20)) {
        lookup->position_aux =
            (OrigLookupPositionPtr PNTR) MemNew((lookup->array_size)*sizeof(OrigLookupPositionPtr));
    }
    
    return lookup;
}

LookupTablePtr LIBCALL
lookup_destruct(LookupTablePtr lookup)

{
    if (lookup == NULL)
        return lookup;
    
    if(lookup->position)
        MemFree(lookup->position);
    if(lookup->position_aux) 
        MemFree(lookup->position_aux);
    if(lookup->pv_array) 
        MemFree(lookup->pv_array);
    if(lookup->mod_lt) 
        MemFree(lookup->mod_lt);
    if(lookup->theTable)
        MemFree(lookup->theTable);
    if (lookup->mod_lookup_table_memory)
        MemFree(lookup->mod_lookup_table_memory);
    lookup_deallocate_memory(lookup);
    lookup = MemFree(lookup);
    
    return lookup;
}

/*
  Deallocated position_aux, which can be large for
  large word-sizes, and is no longer needed.
  
  This now does more than just dealloc lookup->position_aux.
  This is called afeter the original lookup table is fully built.
  So this is where we convert it to the new (mod_lt) table.
  -cfj
*/

Boolean
lookup_position_aux_destruct(LookupTablePtr lookup)
     
{
    if (lookup == NULL)
        return FALSE;
    
    if(lookup->position_aux) 
        lookup->position_aux = MemFree(lookup->position_aux);
    
    
    /* create the new lookup table, Add 1 to silence purify messages. */
    
#ifdef DYNAMIC_CACHE
    lookup->theTable = (Uint4Ptr) MemNew((1+lookup->array_size)*sizeof(Uint4)*(1 + lookup->theCacheSize));
#else
    lookup->mod_lt = (ModLAEntry PNTR) MemNew((1+lookup->array_size)*sizeof(ModLAEntry)); 
#endif
    
    make_mod_lt(lookup);
    
    /* deallocate parts of the old table no longer needed */
    if(lookup->position) lookup->position = MemFree(lookup->position);
    
    return TRUE;
}

/*
  This function adds an (already computed) index to the lookup table 
  for some set of residues.  The most likely use of this is for
  compressed residues, as in blastn.
*/
void LIBCALL
lookup_add_index(LookupTablePtr lookup, Int4 lookup_index, Int4 position, Int1 context)
     
{
    AddPositionToLookupTable(lookup, lookup_index, position, context);
    
    return;
}

/*
  This function adds an index to the lookup table for any number of 
  contiguous resdiues. 
  
*/
void LIBCALL
lookup_add(LookupTablePtr lookup, CharPtr string, Int4 position, Int1 context)

{
    Int4 char_size, lookup_index=0, wordsize;
    
    char_size = lookup->char_size;
    wordsize = lookup->wordsize;
    
    lookup_index = *string;
    wordsize--;
    while (wordsize > 0) {
        lookup_index <<= char_size;
        string++;
        lookup_index += *string;
        wordsize--;
    }
    AddPositionToLookupTable(lookup, lookup_index, position, context);
    
    return;
}

static void 
AddPositionToLookupTable(LookupTablePtr lookup, Int4 index, Int4 position, Int1 context)
     
{
    OrigLookupPositionPtr new, last, PNTR lookup_pos, PNTR lookup_pos_aux;
    
    if (lookup == NULL) return;
    lookup->num_pos_added++;
    
    lookup_pos = lookup->position;
    lookup_pos_aux = lookup->position_aux;
    new = (OrigLookupPositionPtr) lookup_get_memory(lookup, sizeof(OrigLookupPosition));
    if (*(lookup_pos+index) == NULL) {
        lookup->num_unique_pos_added++;
        lookup_pos[index] = new;
        lookup_pos[index]->position = position;
        lookup_pos[index]->context = context;
        if (lookup_pos_aux) {
            lookup_pos_aux[index] = new;
        }
    } else {	/* Go to last link, use aux struct if it exists. */
        if (lookup_pos_aux) {
            last = lookup_pos_aux[index];
            lookup_pos_aux[index] = new;
        } else {
            last = lookup_pos[index];
            while (last->next)
                last = last->next;
        }
        last->next = new;
        new->position = position;
        new->context = context;
    }
    
    return;
}

/*
  Initializes the "lookup_index" for "string".  For contiguous word
  matching the first wordsize-1 letters of string are placed into the 
  lookup_index.  The "new" starting point of the string is then returned.
*/ 

Uint1Ptr LIBCALL
lookup_find_init(LookupTablePtr lookup, Int4 PNTR lookup_index, Uint1Ptr string)
     
{
    
    Int4 char_size, wordsize;
    
    char_size = lookup->char_size;
    wordsize = lookup->wordsize;
    
    *lookup_index = *string;
    
    /* Fill in wordsize-2 spaces. */
    wordsize -= 2;
    while (wordsize > 0) {
        string++;
        *lookup_index <<= char_size;
        *lookup_index += *string;
        wordsize--;
    }
    
    return string;
}
 
/* 
   Allocates memory for the BLAST_WordFinder, and lookup table.
*/

BLAST_WordFinderPtr
BLAST_WordFinderNew (Int2 alphabet_size, Int2 wordsize, Int2 compression_ratio, Boolean round_down, Int4 theCacheSize)
     
{
    BLAST_WordFinderPtr wfp;
    Int2 reduced_wordsize;
    
    wfp = (BLAST_WordFinderPtr) MemNew(sizeof(BLAST_WordFinder));
    
    if (wfp != NULL) {
	/* If the compression_ratio is greater than one and round_down is TRUE, then we want to round 
           the wordsize down, as a multiple of compression_ratio should be three less.  round_down is only
           TRUE in blastn if we want to guarantee that a certain wordsize is used (i.e., 11). */
        if (compression_ratio > 1) {
            if (round_down) {
                reduced_wordsize = (Int2) MIN(2, ((wordsize-3)/compression_ratio));
                wfp->lookup = lookup_new(alphabet_size, (Int2) ((wordsize-3)/compression_ratio), (Int2) reduced_wordsize, theCacheSize);
            } else {
                reduced_wordsize = wordsize/compression_ratio;
                wfp->lookup = lookup_new(alphabet_size, (wordsize/compression_ratio), (Int2) reduced_wordsize, theCacheSize);
            }
        } else {
            wfp->lookup = lookup_new(alphabet_size, (Int2) (wordsize/compression_ratio), 0, theCacheSize);
        }
        wfp->compression_ratio = compression_ratio;
        wfp->wordsize = wordsize;
    }
    
    return wfp;
}

BLAST_WordFinderPtr
BLAST_WordFinderDestruct (BLAST_WordFinderPtr wfp)
     
{
    
    if (wfp != NULL) {
        wfp->lookup = lookup_destruct(wfp->lookup);
        wfp = MemFree(wfp);
    }
    
    return wfp;
}

