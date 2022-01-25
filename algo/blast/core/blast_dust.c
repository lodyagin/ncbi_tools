static char const rcsid[] = "$Id";
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
*
* File Name:  blast_dust.c
*
* Author(s): Ilya Dondoshansky
*   
* Version Creation Date: 05/28/2003
*
* $Revision: 1.16 $
*
* File Description: A utility to find low complexity NA regions.
*                   This parallels functionality of dust.c from the C toolkit,
*                   but without using the structures generated from ASN.1 spec.
* ==========================================================================
*/

#include <algo/blast/core/blast_dust.h>

/* local, file scope, structures and variables */

typedef struct DREGION { /* endpoints */
	struct	DREGION*	next;
	Int4	from, to;
} DREGION;

typedef struct DCURLOC { /* localcurrents */
	Int4	curlevel, curstart, curend;
} DCURLOC;

/* local functions */

static Int4 wo (Int4, Uint1*, Int4, DCURLOC*, Uint1* seq);
static void wo1 (Int4, Uint1*, Int4, DCURLOC*);
static Int4 dust_triplet_find (Uint1*, Int4, Int4, Uint1*);

/* entry point for dusting */

static Int4 dust_segs (Uint1* sequence, Int4 length, Int4 start,
		       DREGION* reg,
		       Int4 level, Int4 windowsize, Int4 minwin, Int4 linker)
{
   Int4    len;
   Int4	i;
   Int4 retlen;
   Uint1* seq;
   DREGION* regold = NULL;
   DCURLOC	cloc;
   Int4	nreg, windowhalf;
   
   /* defaults are more-or-less in keeping with original dust */
   if (level < 2 || level > 64) level = 20;
   if (windowsize < 8 || windowsize > 64) windowsize = 64;
   if (minwin < 4 || minwin > 128) minwin = 4;
   if (linker < 1 || linker > 32) linker = 1;
   windowhalf = windowsize / 2;
   
   nreg = 0;
   seq = (Uint1*) calloc(1, windowsize);			/* triplets */
   if (!seq) {
      return -1;
   }

   for (i = 0; i < length; i += windowhalf) {
      len = (Int4) ((length > i+windowsize) ? windowsize : length-i);
      len -= 2;
      retlen = wo (len, sequence, i, &cloc, seq);
      
      /* get rid of itsy-bitsy's, dust_triplet_find aborts - 
         move 1 triplet away */
      if ((cloc.curend - cloc.curstart + 1) < minwin) {
         if (retlen != len) {
            i += (retlen - windowhalf + 3);
         }
         continue;
      }

      if (cloc.curlevel > level) {
         if (nreg &&
             regold->to + linker >= cloc.curstart+i+start &&
             regold->from <= cloc.curstart + i + start) {
            /* overlap windows nicely if needed */
            regold->to = cloc.curend + i + start;
         } else	{
            /* new window or dusted regions do not overlap */
            reg->from = cloc.curstart + i + start;
            reg->to = cloc.curend + i + start;
            /* 5' edge effects - 3' edge effects - are best handled 
               interactively.
               It probabbly would be good to put 'linker' as an interactive 
               option.
               Interactive means wrapping stuff up in a graphics shell
            */
            regold = reg;
            reg = (DREGION*) calloc(1, sizeof(DREGION));
            if (!reg) {
               sfree(seq);
               return -1;
            }
            reg->next = NULL;
            regold->next = reg;
            nreg++;
         }
         /* kill virtually all 3' tiling anomalies */
         if (cloc.curend < windowhalf)
            i += (cloc.curend - windowhalf);
      }				/* end 'if' high score	*/
   }					/* end for */
   sfree (seq);
   return nreg;
}

static Int4 wo (Int4 len, Uint1* seq_start, Int4 iseg, DCURLOC* cloc, 
                Uint1* seq)
{
	Int4 i, flen;

	cloc->curlevel = 0;
	cloc->curstart = 0;
	cloc->curend = 0;

        /* get the chunk of sequence in triplets */

        memset (seq,0,len+2);        /* Zero the triplet buffer */
	flen = dust_triplet_find (seq_start, iseg, len, seq);

        /* dust the chunk */
	for (i = 0; i < flen; i++)
           wo1 (flen-i, seq+i, i, cloc);

	cloc->curend += cloc->curstart;

	return flen;
}

static void wo1 (Int4 len, Uint1* seq, Int4 iwo, DCURLOC* cloc)
{
   Uint4 sum;
	Int4 loop;
	Int4 newlevel;

	Int2* countsptr;
	Int2 counts[4*4*4];
	memset (counts, 0, sizeof (counts));
/* zero everything */
	sum = 0;
	newlevel = 0;

/* dust loop -- specific for triplets	*/
	for (loop = 0; loop < len; loop++)
	{
		countsptr = &counts[*seq++];
		if (*countsptr)
		{
			sum += (Uint4)(*countsptr);

			newlevel = 10 * sum / loop;

			if (cloc->curlevel < newlevel)
			{
				cloc->curlevel = newlevel;
				cloc->curstart = iwo;
				cloc->curend = loop + 2; /* triplets */
			}
		}
		(*countsptr)++;
	}
	return;
}

#define NCBI_2NA_MASK 0x03
/** Fill an array with 2-bit encoded triplets.
 * @param seq_start Pointer to the start of the sequence in blastna 
 *                  encoding [in]
 * @param icur Offset at which to start extracting triplets [in]
 * @param max Maximal length of the sequence segment to be processed [in]
 * @param s1 Array of triplets [out]
 * @return How far was the sequence processed?
 */
static Int4 
dust_triplet_find (Uint1* seq_start, Int4 icur, Int4 max, Uint1* s1)
{
   Int4 n;
   Uint1* s2,* s3;
   Int2 c;
   Boolean flagVD;
   Uint1* seq = &seq_start[icur];
   Uint1 end_byte = NCBI4NA_TO_BLASTNA[NULLB];
   
   n = 0;
   
   s2 = s1 + 1;
   s3 = s1 + 2;
   
   /* set up needs streamlining */
   /* start again at segment or virtual sequence bounderies */
   /* set up 1 */
   if ((c = *seq++) == end_byte)
      return n;
   c &= NCBI_2NA_MASK;
   *s1 |= c;
   *s1 <<= 2;
   
   /* set up 2 */
   if ((c = *seq++) == end_byte)
      return n;
   c &= NCBI_2NA_MASK;
   *s1 |= c;
   *s2 |= c;
   
   /* triplet fill loop */
   flagVD = TRUE;
   while ((c = *seq++) != end_byte && n < max) {
      c &= NCBI_2NA_MASK;
      *s1 <<= 2;
      *s2 <<= 2;
      *s1 |= c;
      *s2 |= c;
      *s3 |= c;
      s1++;
      s2++;
      s3++;
      n++;
   }
   
   return n;
}

/* look for dustable locations (as slpDust from dust.c */

static Int2 
GetDustLocations (ListNode** loc, DREGION* reg, Int4 nreg)
{
   Int4 i;
   ListNode* last_loc = NULL;
   DoubleInt* dintp;
        
   if (!loc)
      return -1;
   
   *loc = NULL;

   /* point to dusted locations */
   if (nreg > 0) {
      for (i = 0; reg && i < nreg; i++) {
         dintp = (DoubleInt*) calloc(1, sizeof(DoubleInt));
         if (!dintp) {
            return -1;
         }
         dintp->i1 = reg->from;
         dintp->i2 = reg->to;
         if (!last_loc)
            last_loc = ListNodeAddPointer (loc, 0, dintp);
         else 
            last_loc = ListNodeAddPointer (&last_loc, 0, dintp);
         reg = reg->next;
      }
   }
   return 0;
}

Int2 SeqBufferDust (Uint1* sequence, Int4 length, Int4 offset,
                    Int2 level, Int2 window, Int2 minwin, Int2 linker,
                    BlastSeqLoc** dust_loc)
{
	DREGION* reg,* regold;
	Int4 nreg;
   Int2 status = 0;

        /* place for dusted regions */
	regold = reg = (DREGION*) calloc(1, sizeof(DREGION));
	if (!reg)
           return -1;

        nreg = dust_segs (sequence, length, offset, reg, (Int4)level, 
                  (Int4)window, (Int4)minwin, (Int4)linker);

        status = GetDustLocations(dust_loc, reg, nreg);
        /* find tail - this way avoids referencing the pointer */
        while (reg->next) reg = reg->next;


        /* clean up memory */
	reg = regold;
	while (reg)
	{
		regold = reg;
		reg = reg->next;
		sfree (regold);
	}

	return status;
}
