/**************************************************************************
*                                                                         *
*                             COPYRIGHT NOTICE                            *
*                                                                         *
* This software/database is categorized as "United States Government      *
* Work" under the terms of the United States Copyright Act.  It was       *
* produced as part of the author's official duties as a Government        *
* employee and thus can not be copyrighted.  This software/database is    *
* freely available to the public for use without a copyright notice.      *
* Restrictions can not be placed on its present or future use.            *
*                                                                         *
* Although all reasonable efforts have been taken to ensure the accuracy  *
* and reliability of the software and data, the National Library of       *
* Medicine (NLM) and the U.S. Government do not and can not warrant the   *
* performance or results that may be obtained by using this software,     *
* data, or derivative works thereof.  The NLM and the U.S. Government     *
* disclaim any and all warranties, expressed or implied, as to the        *
* performance, merchantability or fitness for any particular purpose or   *
* use.                                                                    *
*                                                                         *
* In any work or product derived from this material, proper attribution   *
* of the author(s) as the source of the software or data would be         *
* appreciated.                                                            *
*                                                                         *
**************************************************************************/
/*****************************************************************************

File name: blastasn.c

Author: Tom Madden

Contents: Functions to produce BLAST ASN.1 (specification 1.8) from data
	produced by BLAST.

Detailed Contents:

        - Functions that produce a SeqAlign from BLAST structures.

        - Functions that produce a BLAST0Result from BLAST structures.

******************************************************************************/

/* $Log: blastasn.c,v $
 * Revision 6.3  1998/04/15 20:13:00  madden
 * MakeBLAST0Result is non-NULL if results are NULL
 *
 * Revision 6.2  1998/03/26 14:15:47  madden
 * Added second (NULL) argument to GetScoreSetFromBlastResultHsp
 *
 * Revision 6.1  1997/12/12 16:41:26  madden
 * charPtr to Uint1Ptr in GetBLAST0SeqData
 *
 * Revision 6.0  1997/08/25 18:33:57  madden
 * Revision changed to 6.0
 *
 * Revision 1.28  1997/07/07 13:58:35  madden
 * switched order of gi and textid
 *
 * Revision 1.27  1997/02/11 01:32:26  kans
 * cast second parameter of GetBLAST0Segment to CharPtr (for CodeWarrior)
 *
 * Revision 1.26  1997/01/15  13:35:06  madden
 * Changes for nucl. alphabets.
 *
 * Revision 1.25  1997/01/13  17:55:44  madden
 * Replaced tabs by spaces.
 *
 * Revision 1.24  1997/01/13  17:16:22  madden
 * Save ncbi4na matrix for blastn.
 *
 * Revision 1.23  1997/01/03  20:30:55  madden
 * Fixed formatting of results.
 *
 * Revision 1.22  1997/01/03  14:51:00  madden
 * Fixed memory leak.
 *
 * Revision 1.21  1996/12/23  14:13:01  madden
 * Changed gap stats printed out.
 *
 * Revision 1.20  1996/12/17  21:36:55  madden
 * Changes to allow deflines for inidividual entries to be retrieved.
 *
 * Revision 1.19  1996/12/17  19:41:20  madden
 * Added gapped stats.
 *
 * Revision 1.18  1996/11/27  16:42:19  madden
 * Produced DbDesc from Readdb information.
 *
 * Revision 1.17  1996/11/14  16:22:52  madden
 * Fixed problem in input to GetTranslation/
 *
 * Revision 1.16  1996/11/13  22:37:28  madden
 * Added call to readdb_get_bioseq for tblast[nx].
 *
 * Revision 1.15  1996/11/08  21:50:48  madden
 * *** empty log message ***
 *
 * Revision 1.14  1996/11/07  22:34:16  madden
 * Replaced call to readdb_get_partial_unpacked_sequence with call
 * to readdb_get_bioseq.
 *
 * Revision 1.13  1996/11/06  22:12:09  madden
 * Changed call to BlastTranslateUnambiguousSequence.
 *
 * Revision 1.12  1996/10/03  13:33:47  madden
 * Changed kbp to kbp[0].
 *
 * Revision 1.11  1996/09/30  21:57:30  madden
 * Replaced ncbi2na alphabet, for query, with blastna alphabet.
 *
 * Revision 1.10  1996/09/25  20:02:43  madden
 * Added "stats" Boolean if stats were collected.
 *
 * Revision 1.9  1996/09/25  14:18:47  madden
 * Removed discontiguous references.
 *
 * Revision 1.8  1996/09/24  18:41:34  madden
 * Fixed conversion to ncbi4na for blastn.
 *
 * Revision 1.7  1996/09/23  17:36:48  madden
 * Changed CharPtr's to Uint1Ptr's.
 *
 * Revision 1.6  1996/09/17  15:29:25  madden
 * Use readdb_get_partial_unpacked_sequence function to fetch sequences
 * for blastn.
 *
 * Revision 1.5  1996/09/17  14:29:43  madden
 * Changes for blastn formatting.
 *
 * Revision 1.4  1996/08/27  21:52:40  madden
 * Added tblastx.
 *
 * Revision 1.3  1996/08/22  21:35:30  madden
 * Corrections to calculation of offsets on negative frames.
 *
 * Revision 1.2  1996/08/14  17:20:40  madden
 * Subject sequence now translated for tblastn.
 *
 * Revision 1.1  1996/08/07  14:07:13  madden
 * Initial revision
 *
 * Revision 1.2  1996/08/06  15:23:58  madden
 * Changes in FillInStdSegInfo to produce seqalign for blastx.
 *
 * Revision 1.1  1996/08/05  19:47:16  madden
 * Initial revision
 *
 * Revision 1.34  1996/08/05  17:26:56  madden
 * Added check for query and subject frame to FillInStdSegInfo.
 *
 * Revision 1.33  1996/07/25  20:45:20  madden
 * Change to calling convention of readdb_get_sequence; protected
 * against some NULL pointers.
 *
 * Revision 1.32  1996/07/25  12:55:57  madden
 * readdb_get_sequence call changed to allow for systems w/o mmap.
 *
 * Revision 1.31  1996/07/24  13:16:28  madden
 * Removed defunct functions asn_coord0 and asn_coord1.
 *
 * Revision 1.30  1996/07/24  12:01:28  madden
 * Changes for blastx
 *
 * Revision 1.29  1996/07/18  13:36:48  madden
 * Addition of the BLASTContextStructPtr.
 *
 * Revision 1.28  1996/07/12  16:32:37  madden
 * Added function GetTheSeqAlignID, fixed errors.
 *
 * Revision 1.27  1996/07/11  16:03:58  madden
 * Corrected production of SeqAlign.
 *
 * Revision 1.26  1996/06/26  15:54:25  madden
 * Dropoff scores for both passes printed out.
 *
 * Revision 1.25  1996/06/20  16:15:57  madden
 * Replaced int's with Int4's.
 *
 * Revision 1.24  1996/06/20  14:11:24  madden
 * Removed unused parameters in FillInDenseDiagInfo and FillInStdSegInfo.
 *
 * Revision 1.23  1996/06/19  14:19:53  madden
 * Allowed SeqId to be added for a sequence comparison.
 *
 * Revision 1.22  1996/06/18  16:05:23  madden
 * fixed problem with length of SeqAlign.
 *
 * Revision 1.21  1996/06/17  21:06:45  madden
 * Added function to produce SeqAlign from a hitlist.
 *
 * Revision 1.20  1996/06/07  13:05:36  madden
 * contiguous/discontiguous printed as parameter.
 *
 * Revision 1.19  1996/06/06  17:55:00  madden
 * Cutoff for second pass now printed.
 *
 * Revision 1.18  1996/06/04  15:33:31  madden
 * Changed GetParameterStack for different programs.
 *
 * Revision 1.17  1996/05/28  14:14:14  madden
 * GetParameterStack changed to include statistics info.
 *
 * Revision 1.16  1996/05/22  20:19:50  madden
 * Changed location used in GetBLAST0Segment to correct one.
 *
 * Revision 1.15  1996/05/16  19:50:15  madden
 * Added documentation block.
 *
 * Revision 1.14  1996/05/01  14:58:51  madden
 * Functions to convert BlastResult structures to BLAST0 structs.
 *
 * Revision 1.13  1996/04/24  12:52:41  madden
 * *** empty log message ***
 *
 * Revision 1.12  1996/04/22  21:40:31  madden
 * *** empty log message ***
 *
 * Revision 1.11  1996/04/03  19:15:02  madden
 * removed DuplicateBLAST0SeqDesc.
 *
 * Revision 1.10  1996/03/29  21:27:06  madden
 * Added functions to produce a BLAST0Result from a SeqAlignPtr.
 *
 * Revision 1.9  1996/03/29  14:08:54  madden
 * GetSeqAlignForSparseHitList added.
 *
 * Revision 1.8  1996/01/23  16:31:31  madden
 * changes to print out parameters.
 *
 * Revision 1.7  1996/01/17  23:18:15  madden
 * Added GetParameterStack function.
 *
 * Revision 1.6  1996/01/08  23:23:22  madden
 * *** empty log message ***
 *
 * Revision 1.5  1996/01/06  18:57:19  madden
 * Changed hsp->n to hsp->num to agree with current defs.
 *
 * Revision 1.4  1995/12/30  19:21:24  madden
 * Added GetBLAST0KABlk.
 *
 * Revision 1.3  1995/12/26  23:03:48  madden
 * moved GetBLAST0Matrix here.
 *
 * Revision 1.2  1995/12/19  22:31:39  madden
 * *** empty log message ***
 *
 * Revision 1.1  1995/12/08  15:48:23  madden
 * Initial revision
 *
 * */

/* VAR_ARGS is used by PrintToValNode. */
#ifdef VAR_ARGS
#include <varargs.h>
#else
#include <stdarg.h>
#endif

#include <blastpri.h>
#include <blastasn.h>
#include <seqport.h>


#define BLASTASN_BUF_SIZE 128

static BLAST0SegmentPtr GetBLAST0Segment PROTO((BlastSearchBlkPtr search, Uint1Ptr sequence, Int4 total_length, Int4 from, Int4 length, Boolean get_seq, Int2 frame, Uint1 alphabet));

static ValNodePtr GetBLAST0SeqData PROTO((BlastSearchBlkPtr search, Uint1Ptr sequence, Int4 offset, Int4 length, Uint1 alphabet));

static BLAST0HitListPtr GetBLAST0HitList PROTO((BlastSearchBlkPtr search, Int4 count, Boolean get_query_seq, Boolean get_db_seq));

static BLAST0SequencePtr GetBLAST0Sequence PROTO((BlastSearchBlkPtr search, Int4 sequence_number, Boolean get_seq));

/* These are used by BLAST to print to a ValNodePtr. */
static void PrintToValNode VPROTO((ValNodePtr PNTR stk, char *format, ...));
static void PrintNewLineToValNode PROTO((ValNodePtr PNTR stp));

/*
	This function assembles a BLAST0DbDescPtr from the information
	in the ReadDBFILEPtr.
*/
BLAST0DbDescPtr LIBCALL
MakeBLAST0DbDesc(ReadDBFILEPtr rdfp)

{
	BLAST0DbDescPtr dbdesc;

	dbdesc = BLAST0DbDescNew();
	if (dbdesc != NULL)
	{
		dbdesc->name = StringSave(readdb_get_filename(rdfp));
		if (readdb_is_prot(rdfp) == TRUE)
			dbdesc->type = 1;
		else
			dbdesc->type = 2;
		dbdesc->def = StringSave(readdb_get_title(rdfp));
		dbdesc->rel_date = NULL; /* Not used or supported. */
		dbdesc->bld_date = StringSave(readdb_get_date(rdfp));
		dbdesc->count = readdb_get_num_entries(rdfp);
		dbdesc->totlen = readdb_get_dblen(rdfp);
		dbdesc->maxlen = readdb_get_maxlen(rdfp);
	}

	return dbdesc;
}

/*************************************************************************
*
*	This function assembles all the components of the BLAST0Result
*	structure.  The result is then passed back to a calling program
*	or made into ASN.1.
***************************************************************************/

BLAST0ResultPtr LIBCALL
MakeBLAST0Result(BlastSearchBlkPtr search, Boolean get_query_seq, Boolean get_db_seq)
{
	BLAST0ResultPtr		brp;
	BLAST0HitListPtr	bhlp;
	BLASTResultsStructPtr   result_struct=NULL;
	Int4			count, total;


	if ((brp = (BLAST0ResultPtr) MemNew(sizeof(BLAST0Result))) == NULL)
		return NULL;

	if (search)
	{
		result_struct = search->result_struct;
		brp->count = result_struct->hitlist_count;
	}
	total = brp->count;

	if (total > 0 && search)
	{
		bhlp = GetBLAST0HitList(search, 0, get_query_seq, get_db_seq);
		brp->hitlists = bhlp;
		for (count=1; count<total; count++)
		{
			bhlp->next = GetBLAST0HitList(search, count, get_query_seq, get_db_seq);
			bhlp = bhlp->next;
		}
	}
	else
	{
		brp->hitlists = NULL;
	}
		
	return brp;
}

static BLAST0HitListPtr
GetBLAST0HitList(BlastSearchBlkPtr search, Int4 count, Boolean get_query_seq, Boolean get_db_seq)

{
	BioseqPtr		bsp;
	BLAST0HitListPtr	bhlp;
	BLASTResultHitlistPtr	result_hitlist;
	BLASTResultHspPtr	hsp;
	BLAST0HSPPtr		blhsp, blhsp_tmp;
	Int4			index, hspcnt, length, hsp_length, prot_length;
	Int4 start, stop;
	SeqPortPtr		spp;
	Uint1Ptr		buffer_start, buffer, sequence, prot_seq;
	Uint1 alphabet;

	if ((bhlp = (BLAST0HitListPtr) MemNew(sizeof(BLAST0HitList))) == NULL)
		return NULL;
	
	result_hitlist = search->result_struct->results[count];
	hspcnt = result_hitlist->hspcnt;
	buffer_start=NULL;

	prot_seq = NULL;
	sequence = NULL;
	blhsp=NULL;
	bsp = readdb_get_bioseq(search->rdfp, result_hitlist->subject_id);

	for (index=0; index<hspcnt; index++)
	{
		hsp = &(result_hitlist->hsp_array[index]);
		blhsp_tmp = (BLAST0HSPPtr) MemNew(sizeof(BLAST0HSP));
		blhsp_tmp->scores = GetScoreSetFromBlastResultHsp(hsp, NULL);
		blhsp_tmp->len = hsp->query_length;
		blhsp_tmp->segs = GetBLAST0Segment(search, search->context[hsp->context].query->sequence+hsp->query_offset, search->context[hsp->context].query->original_length, hsp->query_offset, hsp->query_length, get_query_seq, hsp->query_frame, search->sbp->alphabet_code);
		alphabet = search->sbp->alphabet_code;
		if (get_db_seq)
		{
		    if (StringCmp("blastn", search->prog_name) == 0)
		    {
			length = readdb_get_sequence_length(search->rdfp, result_hitlist->subject_id);
			spp = SeqPortNew(bsp, hsp->subject_offset, hsp->subject_offset + hsp->subject_length - 1, Seq_strand_plus, Seq_code_ncbi4na);
			hsp_length = hsp->subject_length;
			buffer_start = buffer = MemNew(hsp_length*sizeof(Uint1));
			while (hsp_length > 0)
			{
				*buffer = SeqPortGetResidue(spp);
				buffer++;
				hsp_length--;
			}
			spp = SeqPortFree(spp);

			sequence = buffer_start;
		/* readdb_get_partial_unpacked_sequence returns ncbi4na here. */
			alphabet = Seq_code_ncbi4na;
		    }
		    else
		    {
		    	if (StringCmp("tblastn", search->prog_name) == 0 || StringCmp("tblastx", search->prog_name) == 0)
			{
				length = bsp->length;
				if (hsp->subject_frame > 0)
				{
					start = CODON_LENGTH*(hsp->subject_offset) + ABS(hsp->subject_frame) - 1;
					stop = start + CODON_LENGTH*(hsp->subject_length) - 1;
					spp = SeqPortNew(bsp, start, stop, Seq_strand_plus, Seq_code_ncbi4na);
				}
				else
				{
					start = bsp->length - CODON_LENGTH*(hsp->subject_offset + hsp->subject_length) + hsp->subject_frame + 1;
					stop = bsp->length - CODON_LENGTH*(hsp->subject_offset) + hsp->subject_frame;
					spp = SeqPortNew(bsp, start, stop, Seq_strand_minus, Seq_code_ncbi4na);
				}
				hsp_length = CODON_LENGTH*hsp->subject_length;
				buffer_start = buffer = MemNew(hsp_length*sizeof(Uint1));
				while (hsp_length > 0)
				{
					*buffer = SeqPortGetResidue(spp);
					buffer++;
					hsp_length--;
				}
				spp = SeqPortFree(spp);
				hsp_length = CODON_LENGTH*hsp->subject_length;
				prot_seq = GetTranslation(buffer_start, hsp_length, 1, &prot_length, search->db_genetic_code);

			/* The translated sequence starts with a (sentinel) NULLB */
				sequence = prot_seq+1;
		    	}
			else
			{
		    		length = readdb_get_sequence(search->rdfp, result_hitlist->subject_id, &sequence);
				sequence += hsp->subject_offset;
			}
		     }
		}
		else
		{
			length = readdb_get_sequence_length(search->rdfp, result_hitlist->subject_id);
		}
		blhsp_tmp->segs->next = 
			GetBLAST0Segment(search, sequence, length, hsp->subject_offset, hsp->subject_length, get_db_seq, hsp->subject_frame, alphabet);
		if (blhsp)
		{
			blhsp->next = blhsp_tmp;
			blhsp = blhsp->next;
		}
		else
		{
			bhlp->hsps = blhsp_tmp;
			blhsp = blhsp_tmp;
		}
		if (buffer_start)
			buffer_start = MemFree(buffer_start);
		if (prot_seq)
			prot_seq = MemFree(prot_seq);
	}

	bsp = BioseqFree(bsp);

	bhlp->seqs = GetBLAST0Sequence(search, result_hitlist->subject_id, FALSE);

	return bhlp;
}

static BLAST0SequencePtr
GetBLAST0Sequence (BlastSearchBlkPtr search, Int4 sequence_number, Boolean get_seq)

{
	BLAST0SeqDescPtr        desc;
	BLAST0SequencePtr	blsp;
	Boolean not_done;
	Char textid[100];
        CharPtr definition=NULL;
	Uint1Ptr	sequence; 
	Int4		length, gi;
	Uint4		index;
        SeqIdPtr sip=NULL;
        ValNodePtr new_id, vnp;


	if ((blsp = (BLAST0SequencePtr) MemNew(sizeof(BLAST0Sequence))) == NULL)
		return NULL;

	desc = blsp->desc = BLAST0SeqDescNew();
	not_done = TRUE;
	index = 0;
	while (not_done == TRUE)
	{
        	not_done = readdb_get_header(search->rdfp, sequence_number, &index, &sip, &definition);

		new_id = NULL;
/* Only save an id as a textid if it's not a gi! */
       		textid[0] = NULLB;
     		if ((vnp = sip) != NULL)
        	{
       		         while (vnp)
               		 {
                       		 if (vnp->choice == SEQID_GI)
                        	 {
                         		       gi = vnp->data.intvalue;
						ValNodeAddInt(&new_id, BLAST0SeqId_giid, gi);
                                		break;
                       		 }
                       		 vnp = vnp->next;
                	 }        

			 vnp = sip;
               		 while (vnp)
               		 {                   
                     		if (vnp->choice != SEQID_GI)
                     		{
                       		     SeqIdPrint(vnp, textid, PRINTID_FASTA_LONG);
				     ValNodeCopyStr(&new_id, BLAST0SeqId_textid, textid);
                       		     break;
                     		}
                     		vnp = vnp->next;
                	 }    
        	}
		sip = SeqIdSetFree(sip);
		desc->defline = definition;	
		desc->id = new_id;
		if (not_done == TRUE)
		{
			desc->next = BLAST0SeqDescNew();
			desc = desc->next;
		}
	}


	/* This needs to be set up for other codes. */
	blsp->gcode = 1;	/* "1" is the default code for the toolbox. */

	/* seq BLAST0-Seq-data OPTIONAL */
	if (get_seq)
	{
		sequence=NULL;
		length = readdb_get_sequence(search->rdfp, sequence_number, &sequence);
		blsp->seq = GetBLAST0SeqData(search, sequence, 0, length, search->sbp->alphabet_code);
	}
	else
	{
		length = readdb_get_sequence_length(search->rdfp, sequence_number);
	}
	blsp->length = length;

	return blsp;
}

static BLAST0SegmentPtr
GetBLAST0Segment (BlastSearchBlkPtr search, Uint1Ptr sequence, Int4 total_length, Int4 from, Int4 length, Boolean get_seq, Int2 frame, Uint1 alphabet_code)

{
	BLAST0SegmentPtr	blsp;
	BLAST0SeqIntervalPtr	b0sip;
	Uint2		strand;


	if ((blsp = (BLAST0SegmentPtr) MemNew(sizeof(BLAST0Segment))) == NULL)
		return NULL;

	if ((b0sip = (BLAST0SeqIntervalPtr) MemNew(sizeof(BLAST0SeqInterval))) == NULL)
		return NULL;

	strand = 0;
	if (frame > 0)
		strand = BLAST0_Seq_interval_strand_plus;
	else if (frame < 0)
		strand = BLAST0_Seq_interval_strand_minus;

	b0sip->strand = strand;
	if (alphabet_code == Seq_code_ncbistdaa)
	{
		if (frame == 0)
		{ /* blastp */
		     b0sip->from = from;
		     b0sip->to = from + length - 1;
		}
		else if (frame > 0)
		{ /* blastx, tblast[nx] */
		     b0sip->from = CODON_LENGTH*from + frame - 1;
		     b0sip->to = CODON_LENGTH*(from+length) + frame;
		}
		else
		{ /* This may seem counterintuitive, but formatter views it this way. */
		  /* blastx, tblast[nx] */
		     b0sip->from = total_length - CODON_LENGTH*(from+length) + frame + 1;
		     b0sip->to = total_length - CODON_LENGTH*(from) + frame;
		}
	}
	else
	{ /* blastn. */
		if (frame > 0)
		{
		     b0sip->from = from;
		     b0sip->to = from + length - 1;
		}
		else
		{
		     b0sip->to = total_length - from - 1;
		     b0sip->from = total_length - from - length - 2;
		}
	}

	blsp->loc = b0sip;

	/* str BLAST0-Seq-data OPTIONAL */
	if (get_seq)
	{
	/* Save sequence starting with "from" for a distance of "length",
	   where "length" is the length of the HSP. */
		blsp->str = GetBLAST0SeqData(search, sequence, from, length, alphabet_code);
	}
		
	return blsp;
}

static ValNodePtr
GetBLAST0SeqData(BlastSearchBlkPtr search, Uint1Ptr sequence, Int4 offset, Int4 length, Uint1 alphabet_code)
{
	ByteStorePtr byte_sp=NULL;
	ValNodePtr vnp=NULL;
	Int4 index, enc_index, enclen, remainder;
	Uint1Ptr buffer;
        Uint1 byte;

	if (alphabet_code == Seq_code_ncbistdaa)
	{
		byte_sp = BSNew(length);
		BSWrite(byte_sp, sequence, length);
		ValNodeAddPointer(&vnp, BLAST0SeqData_ncbistdaa, byte_sp);

	}
	else if (alphabet_code == Seq_code_ncbi4na)
	{
		enclen = length/2;
		byte_sp = BSNew(enclen + length%2);
		buffer = MemNew((enclen+1)*sizeof(Char));
		enc_index=0;
		for (index=0; index<2*enclen; index += 2)
		{
			byte = sequence[index];
			byte <<= 4;
			byte += sequence[index+1];
			buffer[enc_index] = byte;
			enc_index++;
		}

		remainder = length%2;
		if (remainder > 0)
		{
			byte = sequence[index];
			byte <<= 4;
			buffer[enc_index] = byte;
		}
		BSWrite(byte_sp, buffer, enclen+length%2);
		ValNodeAddPointer(&vnp, BLAST0SeqData_ncbi4na, byte_sp);
		buffer = MemFree(buffer);
	}
	else 
	{ /* Convert blastna to ncbi4na and save. */
		enclen = length/2;
		byte_sp = BSNew(enclen + length%2);
		buffer = MemNew((enclen+1)*sizeof(Char));
		enc_index=0;
		for (index=0; index<2*enclen; index += 2)
		{
			byte = blastna_to_ncbi4na[sequence[index]];
			byte <<= 4;
			byte += blastna_to_ncbi4na[sequence[index+1]];
			buffer[enc_index] = byte;
			enc_index++;
		}

		remainder = length%2;
		if (remainder > 0)
		{
			byte = blastna_to_ncbi4na[sequence[index]];
			byte <<= 4;
			buffer[enc_index] = byte;
		}
		
		BSWrite(byte_sp, buffer, enclen+length%2);
		ValNodeAddPointer(&vnp, BLAST0SeqData_ncbi4na, byte_sp);
		buffer = MemFree(buffer);
	}
	
	return vnp;
}

/**************************************************************************
*	
*	Translates the information in the BLAST_ScoreBlkPtr into the
*	BLAST0MatrixPtr.  If there is a problem (e.g., the BLAST0MatrixPtr
*	cannot be allocated), NULL is returned.
*	The Boolean fullreport specifies whether comments and scores 
*	should be reported.
*
**************************************************************************/

BLAST0MatrixPtr LIBCALL
GetBLAST0Matrix(BLAST_ScoreBlkPtr sbp, Boolean fullreport)
{
	BLAST0MatrixPtr		b0mp;
	BLAST_Score	score;
	Scores_scaled_intsPtr	sco;
	SeqCodeTablePtr sctp;
	Int4		index1, index2, ncbi4na_index1, ncbi4na_index2, total, start_at;
	Uint1		alphabet_code;
	ValNodePtr	vnp, vnp1;

	if (sbp == NULL)
		return NULL;

	if ((b0mp = (BLAST0MatrixPtr) MemNew(sizeof(BLAST0Matrix))) == NULL)
		return NULL;

	alphabet_code=sbp->alphabet_code;

	b0mp->matid = sbp->matid;
	if (sbp->name)
		b0mp->name = StringSave(sbp->name);
	else
		b0mp->name = StringSave("unknown matrix");

	if (fullreport && sbp->comments != NULL) 
	{
		vnp = NULL;
		for (vnp1 = sbp->comments; vnp1 != NULL; vnp1 = vnp1->next) 
			ValNodeCopyStr(&vnp, 0, vnp1->data.ptrvalue);
		b0mp->comments = vnp;
	}

/* Does this need to be fixed for blastn? */
	if (alphabet_code != BLASTNA_SEQ_CODE)
	{
		b0mp->qalpha = alphabet_code;
		b0mp->salpha = alphabet_code;
	}
	else
	{
		b0mp->qalpha = Seq_code_ncbi4na;
		b0mp->salpha = Seq_code_ncbi4na;
	}

	if (fullreport)
	{	/* If a full report is requested, add the scores */
		if (alphabet_code != BLASTNA_SEQ_CODE)
		{
			sctp = SeqCodeTableFindObj(alphabet_code);
			start_at = sctp->start_at;
			total = sctp->start_at + sctp->num;
		}
		else
		{
			start_at = 0;
			total = sbp->mat_dim1;
		}

		vnp=NULL;
		if (alphabet_code == BLASTNA_SEQ_CODE)
		{
			for (index1=start_at; index1<total; index1++) 
			{
				for (index2 = start_at; index2<total; index2++)
				{
					ncbi4na_index1 = blastna_to_ncbi4na[index1];
					ncbi4na_index2 = blastna_to_ncbi4na[index2];
					score = sbp->matrix[ncbi4na_index1][ncbi4na_index2];
					if (score < BLAST_SCORE_1MIN)
						score = INT4_MIN;
					ValNodeAddInt(&vnp, 0, score);
				}
			}
		}
		else
		{
			for (index1=start_at; index1<total; index1++) 
			{
				for (index2 = start_at; index2<total; index2++)
				{
					score = sbp->matrix[index1][index2];
					if (score < BLAST_SCORE_1MIN)
						score = INT4_MIN;
					ValNodeAddInt(&vnp, 0, score);
				}
			}
		}
		sco = (Scores_scaled_intsPtr) MemNew(sizeof(Scores_scaled_ints));
		sco->ints = vnp;
		vnp1=NULL;
		ValNodeAddPointer(&vnp1, Scores_scores_Scores_ScaledInts,sco);
		b0mp->Scores_scores = vnp1;
	}

	return b0mp;
}

/*
	Function to produce the BLAST0KABlkPtr from information
	in the BLAST_ScoreBlkPtr.  Note that most of this
	information (everything except the matid) is in the 
	BLAST_KarlinBlkPtr.
*/

BLAST0KABlkPtr LIBCALL
GetBLAST0KABlk(BLAST_ScoreBlkPtr sbp)

{
	BLAST0KABlkPtr          b0kbp;
	BLAST_KarlinBlkPtr	kbp;
	ValNodePtr		vnp=NULL;


	if (sbp == NULL || sbp->kbp == NULL)
		return NULL;

	kbp = sbp->kbp[0];

	b0kbp = (BLAST0KABlkPtr) MemNew(sizeof(BLAST0KABlk));

	if (b0kbp != NULL)
	{	
		b0kbp->matid = sbp->matid;
		ValNodeAddInt(&vnp, 0, kbp->q_frame);
		ValNodeAddInt(&vnp, 0, kbp->s_frame);
		b0kbp->frames = vnp;
		b0kbp->lambda = kbp->Lambda;
		b0kbp->k = kbp->K;
		b0kbp->h = kbp->H;
	}

	return b0kbp;
}

/*
	this function prints out the (Karlin) parameters to
	the ValNodePtr stp.

	If the Boolean "old" is TRUE, print out stats for one-pass (old)
	method.
	
*/
ValNodePtr LIBCALL
GetParameterStack(BlastSearchBlkPtr search, ValNodePtr stp, Boolean old, Boolean stats)

{
	BLAST_ParameterBlkPtr	pbp;
	BLAST_ScoreBlkPtr	sbp;
	BLAST_KarlinBlkPtr	kbp;
	BLAST0DbDescPtr		dbdesc;

	if (search == NULL || search->pbp == NULL || search->sbp == NULL || search->sbp->kbp == NULL)
		return NULL;

	sbp = search->sbp;
	kbp = sbp->kbp[0];
	pbp = search->pbp;


	PrintNewLineToValNode(&stp);
	PrintNewLineToValNode(&stp);
	PrintToValNode(&stp, "Lambda     K      H");
	PrintNewLineToValNode(&stp);

	if (kbp->Lambda > 0.)
		PrintToValNode(&stp, "%#8.3lg ", kbp->Lambda);
	else
		PrintToValNode(&stp, "     NA  ");
	if (kbp->K > 0.)
		PrintToValNode(&stp, "%#7.3lg ", kbp->K);
	else
		PrintToValNode(&stp, "    NA  ");
	if (kbp->H > 0.)
		PrintToValNode(&stp, "%#7.3lg", kbp->H);
	else
		PrintToValNode(&stp, "    NA  ");


	if (pbp->gapped_calculation)
	{
		PrintNewLineToValNode(&stp);

		PrintToValNode(&stp, "Gapped");
		PrintNewLineToValNode(&stp);
		PrintToValNode(&stp, "Lambda     K      H");
		PrintNewLineToValNode(&stp);
		kbp = sbp->kbp_gap[0];

		if (kbp->Lambda > 0.)
			PrintToValNode(&stp, "%#8.3lg ", kbp->Lambda);
		else
			PrintToValNode(&stp, "     NA  ");
		if (kbp->K > 0.)
			PrintToValNode(&stp, "%#7.3lg ", kbp->K);
		else
			PrintToValNode(&stp, "    NA  ");
		if (kbp->H > 0.)
			PrintToValNode(&stp, "%#7.3lg", kbp->H);
		else
			PrintToValNode(&stp, "    NA  ");

		PrintNewLineToValNode(&stp);
		PrintNewLineToValNode(&stp);
		PrintToValNode(&stp, "OpenGap      ExtendGap     GapX     trigger gapping");
		PrintNewLineToValNode(&stp);
		PrintToValNode(&stp, "%ld               %ld      %ld      %3.1f",
			(long) pbp->gap_open, (long) pbp->gap_extend, (long) pbp->gap_x_dropoff, pbp->gap_trigger);
	}

	PrintNewLineToValNode(&stp);
	PrintNewLineToValNode(&stp);

	if (pbp->two_pass_method == FALSE)
	{
		PrintToValNode(&stp, "E     S     T     X");
	}
	else
	{
		PrintToValNode(&stp, "Cutoff to enter 2nd pass: >= %ld (%4.1f bits)", pbp->cutoff_s_first, pbp->number_of_bits);
		PrintNewLineToValNode(&stp);
		PrintNewLineToValNode(&stp);
		PrintToValNode(&stp, "E     S     T1     T2     X1     X2     W     Gap");
	}

	PrintNewLineToValNode(&stp);

	if (pbp->two_pass_method == FALSE)
	{
		PrintToValNode(&stp, "%3.1f     %ld      %ld      %ld",
			pbp->cutoff_e, (long) pbp->cutoff_s, (long) pbp->threshold_second, (long) ABS(pbp->X));
	}
	else
	{
		PrintToValNode(&stp, "%3.1f      %ld      %ld      %ld      %ld      %ld      %ld      %ld",
			pbp->cutoff_e, (long) pbp->cutoff_s, (long) pbp->threshold_first, (long) pbp->threshold_second, (long) pbp->dropoff_1st_pass, (long) pbp->dropoff_2nd_pass, pbp->window_size, pbp->gap_size);
	}

	PrintNewLineToValNode(&stp);
	PrintNewLineToValNode(&stp);

	dbdesc = MakeBLAST0DbDesc(search->rdfp);
	if (dbdesc)
	{
		PrintToValNode(&stp, "Database:  %s", dbdesc->def);
		PrintNewLineToValNode(&stp);
		PrintToValNode(&stp, "  Posted date:  %s", dbdesc->bld_date);
		PrintNewLineToValNode(&stp);
		PrintToValNode(&stp, "# of letters in database:  %s", Ultostr(dbdesc->totlen,1));
		PrintNewLineToValNode(&stp);
		PrintToValNode(&stp, "# of sequences in database:  %s", Ultostr(dbdesc->count,1));
		PrintNewLineToValNode(&stp);
	}
	dbdesc = BLAST0DbDescFree(dbdesc);

	PrintNewLineToValNode(&stp);

	if (stats)
	{
	  if (pbp->two_pass_method == FALSE)
	  {
	    PrintToValNode(&stp, 
		"Number of Hits to DB: %ld", search->second_pass_hits);
	    PrintNewLineToValNode(&stp);

	    PrintToValNode(&stp, 
		"Number of Sequences: %ld", readdb_get_num_entries(search->rdfp));
	    PrintNewLineToValNode(&stp);

	    PrintToValNode(&stp, 
		"Number of extensions: %ld", search->second_pass_extends);
	    PrintNewLineToValNode(&stp);

	    PrintToValNode(&stp, 
		"Number of successful extensions: %ld", search->second_pass_good_extends);
	   }
	   else
	   {
	    PrintNewLineToValNode(&stp);
	    PrintNewLineToValNode(&stp);
	    
	    PrintToValNode(&stp, 
		"Number of Hits to DB: 1st pass: %ld, 2nd pass: %ld", 
			search->first_pass_hits, search->second_pass_hits);
	    PrintNewLineToValNode(&stp);

	    PrintToValNode(&stp, 
		"Number of Sequences: 1st pass: %ld, 2nd pass: %ld", 
			readdb_get_num_entries(search->rdfp), search->second_pass_trys);
	    PrintNewLineToValNode(&stp);

	    PrintToValNode(&stp, 
		"Number of extensions: 1st pass: %ld, 2nd pass: %ld", 
			search->first_pass_extends, search->second_pass_extends);
	    PrintNewLineToValNode(&stp);

	    PrintToValNode(&stp, 
		"Number of successful extensions: 1st pass: %ld, 2nd pass: %ld", 
			search->first_pass_good_extends, search->second_pass_good_extends);
	  }

	  PrintNewLineToValNode(&stp);
	  PrintToValNode(&stp, 
		"Number of sequences better than %ld: %ld", (long) pbp->cutoff_e, (long) search->number_of_seqs_better_E);

	  PrintNewLineToValNode(&stp);
	  PrintNewLineToValNode(&stp);

	  if (pbp->gapped_calculation)
	  {
	    PrintToValNode(&stp, 
		"Number of HSP's better than %ld without gapping: %ld", (long) pbp->cutoff_e, (long) search->prelim_gap_no_contest);
	    PrintNewLineToValNode(&stp);
	    PrintToValNode(&stp, 
		"Number of HSP's successfully gapped in prelim test: %ld", (long) search->prelim_gap_passed);
	    PrintNewLineToValNode(&stp);
	    PrintToValNode(&stp, 
		"Number of HSP's that attempted gapping in prelim test: %ld", (long) search->prelim_gap_attempts);
	    PrintNewLineToValNode(&stp);
	    PrintToValNode(&stp, 
		"Number of HSP's gapped (non-prelim): %ld", (long) search->real_gap_number_of_hsps);
	    PrintNewLineToValNode(&stp);

		
	  }
	}

	return stp;
}
/*
	PrintToValNode and PrintNewLineToValNode (originally from 
	Warren Gish) are used to print strings to a ValNodePtr.
	These are AsnWritten and decoded by the client.

*/ 

static void
#ifdef VAR_ARGS
PrintToValNode(stk, format, va_alist)
	ValNodePtr	PNTR stk;
	char	*format;
	va_dcl
#else
PrintToValNode(ValNodePtr PNTR stk, char * format, ...)
#endif
{
	va_list	args;
	char	buf[4096];

	if (stk == NULL)
		return;

#ifdef VAR_ARGS
	va_start(args);
#else
	va_start(args, format);
#endif
	vsprintf(buf, format, args);
	va_end(args);

	ValNodeCopyStr(stk, 0, buf);
}

/*
	This function prints a NULLB to the ValNodePtr.  This
	indicates to the formatting programs (BlastPrintValNodeStack 
	in blast2) that a new-line should be inserted.
*/
static void
PrintNewLineToValNode(ValNodePtr PNTR stp)
{
	if (stp == NULL)
		return;
	ValNodeCopyStr(stp, 0, "");
}
