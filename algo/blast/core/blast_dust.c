/* $Id: blast_dust.c,v 1.24 2004/05/19 14:52:02 camacho Exp $
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
 * ==========================================================================
 *
 * Authors: Richa Agarwala (based upon versions variously worked upon by Roma 
 *          Tatusov, John Kuzio, and Ilya Dondoshansky).
 *   
 * ==========================================================================
 */

/** @file blast_dust.c
 * A utility to find low complexity NA regions. This parallels functionality 
 * of dust.c from the C toolkit, but without using the structures generated 
 * from ASN.1 spec.
 */

static char const rcsid[] = 
    "$Id: blast_dust.c,v 1.24 2004/05/19 14:52:02 camacho Exp $";

#include <algo/blast/core/blast_dust.h>
#include <algo/blast/core/blast_util.h>
#include <algo/blast/core/blast_encoding.h>

/* local, file scope, structures and variables */

typedef struct DREGION { /* endpoints */
	struct	DREGION*	next;
	Int4	from, to;
} DREGION;

typedef struct DCURLOC { /* localcurrents */
	Int4	curlevel, curstart, curend;
} DCURLOC;

/* local functions */

static void wo (Int4, Uint1*, Int4, DCURLOC*, Uint1*, Boolean, Int4);
static Boolean wo1 (Int4, Uint1*, Int4, DCURLOC*);
static Int4 dust_triplet_find (Uint1*, Int4, Int4, Uint1*);

/* entry point for dusting */

static Int4 dust_segs (Uint1* sequence, Int4 length, Int4 start,
		       DREGION* reg,
		       Int4 level, Int4 windowsize, Int4 minwin, Int4 linker)
{
   Int4    len;
   Int4	i;
   Uint1* seq;
   DREGION* regold = NULL;
   DCURLOC	cloc;
   Int4	nreg;
   
   /* defaults are more-or-less in keeping with original dust */
   if (level < 2 || level > 64) level = 20;
   if (windowsize < 8 || windowsize > 64) windowsize = 64;
   if (minwin < 4 || minwin > 128) minwin = 4;
   if (linker < 1 || linker > 32) linker = 1;
   
   nreg = 0;
   seq = (Uint1*) calloc(1, windowsize);			/* triplets */
   if (!seq) {
      return -1;
   }

   len = (Int4) ((length > windowsize) ? windowsize : length);
   len -= 2;
   dust_triplet_find (sequence, 0, len-1, seq+1);

   for (i = 0; i < length-2; i++) {
      len = (Int4) ((length > i+windowsize) ? windowsize : length-i);
      len -= 2;
      if ((length >= i+windowsize) || (i==0))
          wo (len, sequence, i, &cloc, seq, TRUE, level);
      else /* remaining portion of sequence is less than windowsize */
          wo (len, sequence, i, &cloc, seq, FALSE, level);
      
      if (cloc.curlevel > level) {
         if (nreg &&
             regold->to + linker >= cloc.curstart+i+start &&
             regold->from <= cloc.curend + i + start + linker) {
            /* overlap windows nicely if needed */
            if (regold->to < cloc.curend + i + start)
                regold->to = cloc.curend + i + start;
            if (regold->from > cloc.curstart + i + start)
                regold->from = cloc.curstart + i + start;
         } else	{
            /* new window or dusted regions do not overlap */
            reg->from = cloc.curstart + i + start;
            reg->to = cloc.curend + i + start;
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
      }				/* end 'if' high score	*/
   }					/* end for */
   sfree (seq);
   return nreg;
}

static void wo (Int4 len, Uint1* seq_start, Int4 iseg, DCURLOC* cloc, 
                Uint1* seq, Boolean FIND_TRIPLET, Int4 level)
{
	Int4 smaller_window_start, mask_window_end;
        Boolean SINGLE_TRIPLET;

	cloc->curlevel = 0;
	cloc->curstart = 0;
	cloc->curend = 0;

	if (len < 1)
		return;

        /* get the chunk of sequence in triplets */
	if (FIND_TRIPLET==TRUE) /* Copy suffix as prefix and find one */
	{
		memmove(seq,seq+1,(len-1)*sizeof(Uint1));
		seq[len-1] = seq[len] = seq[len+1] = 0;
		dust_triplet_find (seq_start, iseg+len-1, 1, seq+len-1);
	}
	else /* Copy suffix */
		memmove(seq,seq+1,len*sizeof(Uint1));

        /* dust the chunk */
	SINGLE_TRIPLET = wo1 (len, seq, 0, cloc); /* dust at start of window */

        /* consider smaller windows only if anything interesting 
           found for starting position  and smaller windows have potential of
           being at higher level */
	if ((cloc->curlevel > level) && (!SINGLE_TRIPLET)) {
		mask_window_end = cloc->curend-1;
		smaller_window_start = 1;
                while ((smaller_window_start < mask_window_end) &&
                       (!SINGLE_TRIPLET)) {
			SINGLE_TRIPLET = wo1(mask_window_end-smaller_window_start,
                             seq+smaller_window_start, smaller_window_start, cloc);
                	smaller_window_start++;
	        }
	}

	cloc->curend += cloc->curstart;
}

/* returns TRUE if there is single triplet in the sequence considered */
static Boolean wo1 (Int4 len, Uint1* seq, Int4 iwo, DCURLOC* cloc)
{
   Uint4 sum;
	Int4 loop;
	Int4 newlevel;

	Int2* countsptr;
	Int2 counts[4*4*4];
	Uint1 triplet_count = 0;

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
		else
			triplet_count++;
		(*countsptr)++;
	}

	if (triplet_count > 1)
		return(FALSE);
	return(TRUE);
}

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
   Uint1* seq = &seq_start[icur];
   Uint1 end_byte = NCBI4NA_TO_BLASTNA[NULLB];
   
   n = 0;
   
   s2 = s1 + 1;
   s3 = s1 + 2;
   
   /* set up 1 */
   if ((c = *seq++) == end_byte)
      return n;
   c &= NCBI2NA_MASK;
   *s1 |= c;
   *s1 <<= 2;
   
   /* set up 2 */
   if ((c = *seq++) == end_byte)
      return n;
   c &= NCBI2NA_MASK;
   *s1 |= c;
   *s2 |= c;
   
   /* triplet fill loop */
   while (n < max && (c = *seq++) != end_byte) {
      c &= NCBI2NA_MASK;
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

/* look for dustable locations (as slpDust from dust.c) */

static Int2 
GetDustLocations (ListNode** loc, DREGION* reg, Int4 nreg)
{
   Int4 i;
   ListNode* last_loc = NULL;
   SSeqRange* dintp;
        
   if (!loc)
      return -1;
   
   *loc = NULL;

   /* point to dusted locations */
   if (nreg > 0) {
      for (i = 0; reg && i < nreg; i++) {
         dintp = (SSeqRange*) calloc(1, sizeof(SSeqRange));
         if (!dintp) {
            return -1;
         }
         dintp->left = reg->from;
         dintp->right = reg->to;
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
