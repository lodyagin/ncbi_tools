/* ===========================================================================
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
* ===========================================================================*/

/*****************************************************************************

File name: blastool.c

Author: Tom Madden

Contents: Utilities for BLAST

******************************************************************************/
/*
* $Revision: 6.74 $
* $Log: blastool.c,v $
* Revision 6.74  2000/01/21 22:24:09  madden
* Use Nlm_Int8tostr in place of Ltostr
*
* Revision 6.73  2000/01/13 18:10:43  madden
* Fix problem with incorrect stat values for blastn and missing hits
*
* Revision 6.72  2000/01/07 16:01:24  madden
* Use readdb_get_totals_ex to get db number to report
*
* Revision 6.71  1999/12/31 14:23:19  egorov
* Add support for using mixture of real and maks database with gi-list files:
* 1. Change logic of creating rdfp list.
* 2. BlastGetDbChunk gets real databases first, then masks.
* 3. Propoper calculation of database sizes using alias files.
* 4. Change to CommonIndex to support using of mask databases.
* 5. Use correct gis in formated output (BlastGetAllowedGis()).
* 6. Other small changes
*
* Revision 6.70  1999/12/21 20:03:23  egorov
* readdb_gi2seq() has new parameter.  Use NULL here.
*
* Revision 6.69  1999/12/21 16:58:00  madden
* Fix command-line parser for options in case there is no space between option and value
*
* Revision 6.68  1999/12/17 20:47:04  egorov
* Fix 'gcc -Wall' warnings
*
* Revision 6.67  1999/12/16 19:16:49  egorov
* Return a value from not-void function
*
* Revision 6.66  1999/12/14 15:35:13  madden
* Added BlastPrintFilterWarning
*
* Revision 6.65  1999/11/30 19:00:51  madden
* Added Nlm_SwapUint4 calls for the ordinal ID list
*
* Revision 6.64  1999/11/26 22:26:49  madden
* Change gap_x_dropoff value for blastn
*
* Revision 6.63  1999/11/09 14:16:53  madden
* made sum_stats default again, rolling back rev 6.61
*
* Revision 6.62  1999/11/02 15:22:18  madden
* Add BlastParceInputString and BlastGetLetterIndex
*
* Revision 6.61  1999/10/27 21:01:32  madden
* made do_sum_stats not default
*
* Revision 6.60  1999/09/29 17:14:50  shavirin
* Fixed memory leak in BLASTOptionsDelete()
*
* Revision 6.59  1999/09/22 20:59:20  egorov
* Add blast db mask stuff
*
* Revision 6.58  1999/09/16 16:54:43  madden
* Allow longer wordsizes
*
* Revision 6.57  1999/09/14 19:56:39  shavirin
* Fixed bug in PHI-Blast when number of hits to DB == 0
*
* Revision 6.56  1999/09/09 18:00:25  madden
* formatting problem for effective db length
*
* Revision 6.55  1999/08/27 20:22:51  shavirin
* Added default value for decline_align in the function BLASTOptionNew().
*
* Revision 6.53  1999/08/17 14:10:05  madden
* Validation for smith_waterman and tweak_parameters options
*
* Revision 6.52  1999/05/25 13:37:15  madden
* Call readdb_get_sequence_length only for seqs in db
*
* Revision 6.51  1999/04/27 17:22:27  madden
* Set hsp_num_max for ungapped BLAST
*
* Revision 6.50  1999/04/14 14:53:48  madden
* Correction for databases over 2 Gig
*
* Revision 6.49  1999/04/01 21:42:47  madden
* Fix memory leaks when gi list is used
*
* Revision 6.48  1999/03/22 15:19:03  beloslyu
* corrections to compile on Red Hat Linux v5.2
*
* Revision 6.47  1999/03/18 16:43:57  shavirin
* Added function Boolean HeyIAmInMemory(Int4 program)
*
* Revision 6.46  1999/03/17 13:21:07  madden
* Fix comment in comment problem
*
* Revision 6.45  1999/03/12 15:03:45  egorov
* Add proper Int4-long type casting
*
* Revision 6.44  1999/02/19 14:18:25  madden
* Added back check for negative nucl. penalty
*
* Revision 6.43  1999/02/19 14:16:20  madden
* list manipulation bug for seed
*
* Revision 6.41  1999/01/26 18:26:54  madden
* make updateLambdaK public
*
* Revision 6.40  1999/01/08 22:08:42  madden
* BlastScaleMatrix returns factor as FloatHi
*
 * Revision 6.39  1998/12/31 18:17:06  madden
 * Added strand option
 *
 * Revision 6.38  1998/12/29 17:45:07  madden
 * Add do_sum_stats flag
 *
 * Revision 6.37  1998/12/03 15:19:32  madden
 * Changes to speed up BlastFreeHeap and InsertToHeap
 *
 * Revision 6.36  1998/11/04 01:36:06  egorov
 * Add support for entrez-query and org-name to blast3
 *
 * Revision 6.35  1998/10/13 22:06:24  madden
 * Fixed AdjustDbNumbers
 *
 * Revision 6.34  1998/09/28 12:29:24  madden
 * Check for problem in rescaling code
 *
 * Revision 6.33  1998/09/22 18:46:43  egorov
 * Add BlastErrorPrintExtra()
 *
 * Revision 6.32  1998/09/17 19:53:02  madden
 * Added fillCandLambda
 *
 * Revision 6.31  1998/09/16 18:59:35  madden
 * Print subset information if entire database not searched
 *
 * Revision 6.30  1998/09/14 15:48:36  madden
 * Fixed PHI-BLAST reference
 *
 * Revision 6.29  1998/09/14 15:11:15  egorov
 * Add support for Int8 length databases; remove unused variables
 *
 * Revision 6.28  1998/09/10 22:36:09  madden
 * Added convertSeqAlignListToValNodeList and convertValNodeListToSeqAlignList
 *
 * Revision 6.27  1998/09/09 21:18:09  madden
 * Added PrintKAParametersExtra
 *
 * Revision 6.26  1998/09/04 14:45:42  madden
 * Moved code from blast.c blastool.c
 *
 * Revision 6.25  1998/08/28 21:21:29  madden
 * Changed PhiBlast ref
 *
 * Revision 6.24  1998/08/25 14:16:22  madden
 * Added BlastGetPhiReference and BlastPrintPhiReference
 *
 * Revision 6.23  1998/07/28 15:23:47  madden
 * Changed Number of sequences better printout again
 *
 * Revision 6.22  1998/07/27 21:54:02  madden
 * Fixed printing of non-integral E-values
 *
 * Revision 6.21  1998/07/21 20:58:06  madden
 * Changes to allow masking at hash only
 *
 * Revision 6.20  1998/07/17 15:40:01  madden
 * Changes for Effective search space.
 *
 * Revision 6.19  1998/07/02 22:15:31  madden
 * Fixed bug in BlastAdjustDbNumbers
 *
 * Revision 6.18  1998/06/17 18:10:24  madden
 * Validate for isPatternSearch
 *
 * Revision 6.17  1998/06/12 15:52:49  madden
 * Fixed warnings
 *
 * Revision 6.16  1998/06/05 15:50:35  madden
 * Return 1 from BLASTOptionValidateEx if wordsize incorrect.
 *
 * Revision 6.15  1998/06/03 17:40:48  madden
 * Added blastn check for wordsize in BLASTOptionValidateEx
 *
 * Revision 6.14  1998/05/28 19:59:42  madden
 * Changed hsp_range_max
 *
 * Revision 6.13  1998/05/21 19:44:45  egorov
 * Make word "Reference" be HTML link in case HTML output requested
 *
 * Revision 6.12  1998/05/03 17:21:04  madden
 * Added error message for expect_value <= 0, fix typo setting open and expect values
 *
 * Revision 6.11  1998/05/01 18:33:56  egorov
 * Add new parametes to BLASTOptionSetGapParam()
 *
 * Revision 6.10  1998/04/30 14:28:46  madden
 * Raise thresholds for blastx, tblast[nx]
 *
 * Revision 6.9  1998/04/29 14:28:07  madden
 * Fix reference formatting problem
 *
 * Revision 6.8  1998/04/27 16:47:34  madden
 * Added window and threshold to BLASTOptionSetGapParams
 *
 * Revision 6.7  1998/04/24 19:28:29  madden
 * Added BlastScaleMatrix (and other rescaling code moved from posit.c)
 *
 * Revision 6.6  1998/04/13 20:29:57  madden
 * Add one to length of array for NULLB at end
 *
 * Revision 6.5  1998/03/24 15:38:23  madden
 * Use BlastDoubleInt4Ptr to keep track of gis and ordinal_ids
 *
 * Revision 6.4  1998/03/18 14:14:18  madden
 * Support random access by gi list
 *
 * Revision 6.3  1998/02/28 17:24:30  madden
 * Default window_size zero for blastn
 *
 * Revision 6.2  1998/02/27 16:52:05  madden
 * Added BlastGetSequenceFromBioseq
 *
 * Revision 6.1  1998/02/27 14:30:30  madden
 * Tools (or utilities) for the BLAST programs
 *
*/

#include <ncbi.h>
#include <blastpri.h>
#include <objcode.h>
#include <objseq.h>
#include <sequtil.h>
#include <readdb.h>
#include <ncbithr.h>
#include <txalign.h>
#include <posit.h>
#include <seed.h>

#ifdef OS_UNIX
/* Here is function to calculate number of BLAST jobs running on spesific
   computer. This function uses array of 5 semaphores with the single key
   0xACACACAC. Each semaphore in this array represents spesific BLAST
   program, so at any given time we know exactly number of processes for
   each BLAST program. Function will increase specific semaphore by one
   and OS will reset this back whenever process will exit from memory */

#include <sys/types.h>
#if defined(__linux__)
/* that is needed because of the bug in Linux Red Hat v5.2 */
#define __USE_XOPEN
typedef __key_t key_t;
#endif
#include <sys/ipc.h>
#include <sys/sem.h>

#define BLAST_SERVER_KEY 0xACACACAC

static struct sembuf BLASTSemaCmd[5] =
{
    {blast_type_blastn - 1,  1, SEM_UNDO},  /* Increase 0 semaphore by one  - blastn  */
    {blast_type_blastp - 1,  1, SEM_UNDO},  /* Increase 1 semaphore by one  - blastp  */
    {blast_type_blastx - 1,  1, SEM_UNDO},  /* Increase 2 semaphore by one  - blastx  */
    {blast_type_tblastn - 1, 1, SEM_UNDO},  /* Increase 3 semaphore by one  - tblastn */
    {blast_type_tblastx - 1, 1, SEM_UNDO}  /* Increase 4 semaphore by one  - tblastx */
};

static Boolean HeyAlreadyCalled = FALSE;

Boolean HeyIAmInMemory(Int4 program)
{
    register int id;
    
    if(HeyAlreadyCalled)
        return TRUE;

    if(program < blast_type_blastn || program > blast_type_tblastx)
        return FALSE;
    
    if ((id = semget (BLAST_SERVER_KEY, 5, 0666 | IPC_CREAT)) < 0) {
        ErrPostEx(SEV_ERROR, 0, 0, "Cannot create semaphore\n");
        return FALSE;     /* permission problem or tables full */
    }
    
    /* Increasing "program" semaphore by one */
    
    if(semop(id, &BLASTSemaCmd[program - 1], 1) < 0) {
        ErrPostEx(SEV_ERROR, errno, 0, "semop: cannot increase by one");
        return FALSE;
    }
    
    HeyAlreadyCalled = TRUE;
    
    return TRUE;
}

#endif

#define BUFFER_LENGTH 255

/*
	adds the new string to the buffer, separating by a tilde.
	Checks the size of the buffer for FormatBlastParameters and
	allocates longer replacement if needed.
*/

static Boolean 
add_string_to_bufferEx(CharPtr buffer, CharPtr *old, Int2Ptr old_length, Boolean add_tilde)

{
	CharPtr new, ptr;
	Int2 length, new_length;

	length = (StringLen(*old));

	if((StringLen(buffer)+length+3) > *old_length)
	{
		new_length = *old_length + 255;
		new = MemNew(new_length*sizeof(Char));
		if (*old_length > 0 && *old != NULL)
		{
			MemCpy(new, *old, *old_length);
			*old = MemFree(*old);
		}
		*old = new;
		*old_length = new_length;
	}

	ptr = *old;
	ptr += length;
	if (add_tilde)
	{
		*ptr = '~';
		ptr++;
	}

	while (*buffer != NULLB)
	{
		*ptr = *buffer;
		buffer++; ptr++;
	}

	return TRUE;
}

static Boolean 
add_string_to_buffer(CharPtr buffer, CharPtr *old, Int2Ptr old_length)

{
	return add_string_to_bufferEx(buffer, old, old_length, TRUE);
}

/*
	Formats the BLAST parameters for the BLAST report.
	One CharPtr is returned, newlines are indicated by tildes ('~').
*/	


CharPtr
FormatBlastParameters(BlastSearchBlkPtr search)

{
	BLAST_ParameterBlkPtr   pbp;
	BLAST_Score cutoff;
	Char buffer[128];
	CharPtr ret_buffer;
	Int2 ret_buffer_length;
	Int4 num_entries;
	Int8 total_length;
	Nlm_FloatHi evalue;

	pbp = search->pbp;

	ret_buffer = NULL;
	ret_buffer_length = 0;

	
	sprintf(buffer, "Matrix: %s", search->sbp->name);
	add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
	
	if (pbp->gapped_calculation)
	{
		sprintf(buffer, "Gap Penalties: Existence: %ld, Extension: %ld", (long) search->pbp->gap_open, (long) search->pbp->gap_extend);
		add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
	}

	if (pbp->two_pass_method == FALSE)
	{
		sprintf(buffer, "Number of Hits to DB: %ld", (long) search->second_pass_hits);
		add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
		
		readdb_get_totals_ex(search->rdfp, &total_length, &num_entries, TRUE);
		
	    	sprintf(buffer, "Number of Sequences: %ld", (long) num_entries);
		add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
		sprintf(buffer, "Number of extensions: %ld", (long) search->second_pass_extends);
		add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
		sprintf(buffer, "Number of successful extensions: %ld", (long) search->second_pass_good_extends);
		add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
	}
	else
	{
		sprintf(buffer, "Number of Hits to DB: 1st pass: %ld, 2nd pass: %ld", 
			(long) search->first_pass_hits, (long) search->second_pass_hits);
		add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
		readdb_get_totals_ex(search->rdfp, &total_length, &num_entries, TRUE);
		sprintf(buffer, "Number of Sequences: 1st pass: %ld, 2nd pass: %ld", 
			(long) num_entries, (long) search->second_pass_trys);
		add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
		sprintf(buffer, "Number of extensions: 1st pass: %ld, 2nd pass: %ld", 
			(long) search->first_pass_extends, (long) search->second_pass_extends);
		add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
		sprintf(buffer, "Number of successful extensions: 1st pass: %ld, 2nd pass: %ld", 
			(long) search->first_pass_good_extends, (long) search->second_pass_good_extends);
		add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
	}

	if (pbp->cutoff_e > 0.1)
	{
		sprintf(buffer, "Number of sequences better than %4.1f: %ld", 
			pbp->cutoff_e, (long) search->number_of_seqs_better_E);
	}
	else
	{
		sprintf(buffer, "Number of sequences better than %3.1e: %ld", 
			pbp->cutoff_e, (long) search->number_of_seqs_better_E);
	}
	add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);

	if (pbp->gapped_calculation &&
		StringCmp(search->prog_name, "blastn") != 0)
	{
		sprintf(buffer, "Number of HSP's better than %4.1f without gapping: %ld", 
			pbp->cutoff_e, (long) search->prelim_gap_no_contest);
		add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
		sprintf(buffer, "Number of HSP's successfully gapped in prelim test: %ld", 
			(long) search->prelim_gap_passed);
		add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
		sprintf(buffer, "Number of HSP's that attempted gapping in prelim test: %ld", (long) search->prelim_gap_attempts);
		add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
		sprintf(buffer, "Number of HSP's gapped (non-prelim): %ld", (long) search->real_gap_number_of_hsps);
		add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
	}

	sprintf(buffer, "length of query: %ld", (long) search->context[search->first_context].query->length);
	add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
	sprintf(buffer, "length of database: %s", Nlm_Int8tostr ((Int8) search->dblen, 1));
	add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);

	sprintf(buffer, "effective HSP length: %ld", (long) search->length_adjustment);
	add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
	sprintf(buffer, "effective length of query: %ld", (long) search->context[search->first_context].query->effective_length);
	add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
	sprintf(buffer, "effective length of database: %s", Nlm_Int8tostr ((Int8) search->dblen_eff, 1));
	add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
	sprintf(buffer, "effective search space: %8.0f", ((Nlm_FloatHi) search->dblen_eff)*((Nlm_FloatHi) search->context[search->first_context].query->effective_length));
	add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
	sprintf(buffer, "effective search space used: %8.0f", (Nlm_FloatHi) search->searchsp_eff);
	add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);


	if (StringCmp(search->prog_name, "blastx") == 0 || StringCmp(search->prog_name, "tblastn") == 0 || StringCmp(search->prog_name, "tblastx") == 0)
	{
		sprintf(buffer, "frameshift window, decay const: %ld, %4.1f",
			(long) pbp->gap_size, pbp->gap_decay_rate);
		add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
	}

	sprintf(buffer, "T: %ld", (long) search->pbp->threshold_second);
	add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
	sprintf(buffer, "A: %ld", (long) search->pbp->window_size);
	add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
	sprintf(buffer, "X1: %ld (%4.1f bits)", (long) (-search->pbp->dropoff_1st_pass), ((-search->pbp->dropoff_1st_pass)*(search->sbp->kbp[search->first_context]->Lambda/NCBIMATH_LN2)));
	add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
	if (StringCmp(search->prog_name, "blastn") == 0 || search->pbp->gapped_calculation == FALSE)
	{
		sprintf(buffer, "X2: %ld (%4.1f bits)", (long) search->pbp->gap_x_dropoff, ((search->pbp->gap_x_dropoff)*(search->sbp->kbp[search->first_context]->Lambda/NCBIMATH_LN2)));
		add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
	}
	else
	{
		sprintf(buffer, "X2: %ld (%4.1f bits)", (long) search->pbp->gap_x_dropoff, ((search->pbp->gap_x_dropoff)*(search->sbp->kbp_gap[search->first_context]->Lambda/NCBIMATH_LN2)));
		add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
		sprintf(buffer, "X3: %ld (%4.1f bits)", (long) search->pbp->gap_x_dropoff_final, ((search->pbp->gap_x_dropoff_final)*(search->sbp->kbp_gap[search->first_context]->Lambda/NCBIMATH_LN2)));
		add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
	}
	sprintf(buffer, "S1: %ld (%4.1f bits)", (long) search->pbp->gap_trigger, ((((search->pbp->gap_trigger)*(search->sbp->kbp[search->first_context]->Lambda))-(search->sbp->kbp[search->first_context]->logK))/NCBIMATH_LN2));
	add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);
	cutoff = 0;
	evalue = pbp->cutoff_e;
	if (StringCmp(search->prog_name, "blastn") == 0 || search->pbp->gapped_calculation == FALSE)
	{
		BlastCutoffs(&cutoff, &evalue, search->sbp->kbp[search->first_context], (Nlm_FloatHi) search->context[search->first_context].query->effective_length, (Nlm_FloatHi) search->dblen_eff, FALSE);
		sprintf(buffer, "S2: %ld (%4.1f bits)", (long) cutoff, (((cutoff)*(search->sbp->kbp[search->first_context]->Lambda))-(search->sbp->kbp[search->first_context]->logK))/NCBIMATH_LN2);
	}
	else
	{
		BlastCutoffs(&cutoff, &evalue, search->sbp->kbp_gap[search->first_context], (Nlm_FloatHi) search->context[search->first_context].query->effective_length, (Nlm_FloatHi) search->dblen_eff, FALSE);
		sprintf(buffer, "S2: %ld (%4.1f bits)", (long) cutoff, (((cutoff)*(search->sbp->kbp_gap[search->first_context]->Lambda))-(search->sbp->kbp_gap[search->first_context]->logK))/NCBIMATH_LN2);
	}
	add_string_to_buffer(buffer, &ret_buffer, &ret_buffer_length);

	return ret_buffer;
}

/*
	Print the buffer, adding newlines where tildes are found.
*/

Boolean LIBCALL
PrintTildeSepLines(CharPtr buffer, Int4 line_length, FILE *outfp)

{
	if (outfp == NULL || buffer == NULL)
		return FALSE;

	asn2ff_set_output(outfp, NULL);

	ff_StartPrint(0, 0, line_length, NULL);
	while (*buffer != NULLB)
	{
		if (*buffer != '~')
			ff_AddChar(*buffer);
		else
			NewContLine();
		buffer++;
	}
	ff_EndPrint();

	return TRUE;
}

/*
	Print the Karlin-Altschul parameters.

	if gapped is TRUE, then slightly different formatting is used.
*/

Boolean LIBCALL
PrintKAParameters(Nlm_FloatHi Lambda, Nlm_FloatHi K, Nlm_FloatHi H, Int4 line_length, FILE *outfp, Boolean gapped)

{
	return PrintKAParametersExtra(Lambda, K, H, 0.0, line_length, outfp, gapped);
}

Boolean LIBCALL
PrintKAParametersExtra(Nlm_FloatHi Lambda, Nlm_FloatHi K, Nlm_FloatHi H, Nlm_FloatHi C, Int4 line_length, FILE *outfp, Boolean gapped)

{
	Char buffer[BUFFER_LENGTH];

	if (outfp == NULL)
		return FALSE;

	asn2ff_set_output(outfp, NULL);

	ff_StartPrint(0, 0, line_length, NULL);
	if (gapped)
	{
		ff_AddString("Gapped");
		NewContLine();
	}
	
	if (C == 0.0)
		ff_AddString("Lambda     K      H");
	else
		ff_AddString("Lambda     K      H      C");
	NewContLine();
	sprintf(buffer, "%#8.3g ", Lambda);
	ff_AddString(buffer);
	sprintf(buffer, "%#8.3g ", K);
	ff_AddString(buffer);
	sprintf(buffer, "%#8.3g ", H);
	ff_AddString(buffer);
	if (C != 0.0)
	{
		sprintf(buffer, "%#8.3g ", C);
		ff_AddString(buffer);
	}
	NewContLine();
	ff_EndPrint();

	return TRUE;

}

/*
	Deallocates *BlastErrorMsgPtr produced by BlastConstructErrorMessage.
*/

BlastErrorMsgPtr BlastDestroyErrorMessage(BlastErrorMsgPtr error_msg)

{
	if (error_msg == NULL)
		return NULL;

	MemFree(error_msg->msg);
	MemFree(error_msg);

	return NULL;
}

/* 
	Prepares error message and appends ValNodePtr, containing BlastErrorMsgPtr, to
	end of chain.  The beginning of the ValNodePtr chain is returned.
*/

ValNodePtr BlastConstructErrorMessage(CharPtr function, CharPtr message, Uint1 level, ValNodePtr PNTR vnpp)

{
	Char buffer[BUFFER_LENGTH];
	CharPtr ptr;
	BlastErrorMsgPtr error_msg;

	if (vnpp == NULL)
		return NULL;

	buffer[0] = NULLB;
	ptr = buffer;
	if (function != NULL)
	{
		sprintf(buffer, "%s: ", function);
		ptr = buffer + StringLen(buffer);
	}
	
	if (message != NULL)
	{
		sprintf(ptr, "%s", message);
	}

	error_msg = (BlastErrorMsgPtr) MemNew(sizeof(BlastErrorMsg));
	error_msg->msg = StringSave(buffer);
	error_msg->level = level;

	ValNodeAddPointer(vnpp, 0, error_msg);

	return *vnpp;
}

/*
	Destroys a chain of ValNodes and the BlastErrorMsgPtr data.
*/

ValNodePtr BlastErrorChainDestroy(ValNodePtr vnp)

{
	ValNodePtr start = vnp;

	while (vnp)
	{
		BlastDestroyErrorMessage(vnp->data.ptrvalue);
		vnp->data.ptrvalue = NULL;
		vnp = vnp->next;
	}

	ValNodeFree(start);

	return NULL;
}

/*
	Prints the error messages.
*/

void LIBCALL BlastErrorPrint(ValNodePtr error_return)

{
	BlastErrorMsgPtr error_msg;

	if (error_return == NULL)
		return;

	while (error_return)
	{
		error_msg = error_return->data.ptrvalue;
		switch (error_msg->level)
		{
			case 0:
				ErrPostEx(SEV_INFO, 0, 0, "%s", error_msg->msg);
				break;
			case 1:
				ErrPostEx(SEV_WARNING, 0, 0, "%s", error_msg->msg);
				break;
			case 2:
				ErrPostEx(SEV_ERROR, 0, 0, "%s", error_msg->msg);
				break;
			case 3:
				ErrPostEx(SEV_FATAL, 0, 0, "%s", error_msg->msg);
				break;
			default:
				ErrPostEx(SEV_WARNING, 0, 0, "Unknown BLAST error level");
				break;
		}
		error_return = error_return->next;
	}
	return;
	
}

void LIBCALL BlastErrorPrintExtra(ValNodePtr error_return, Boolean errpostex, FILE* fp)
{
    BlastErrorMsgPtr	error_msg;
    ErrSev		err_sev;
    CharPtr		default_msg = "Unknown BLAST error level", msg;
    CharPtr		errsevmsg,
        errsevmsg_0 = "INFO",
        errsevmsg_1 = "WARNING",
        errsevmsg_2 = "ERROR",
        errsevmsg_3 = "FATAL";
    

    
    if (error_return == NULL)
        return;
    
    while (error_return)
    {
        error_msg = error_return->data.ptrvalue;
        msg = error_msg->msg;
        
        switch (error_msg->level)
        {
            case 0:
                err_sev = SEV_INFO;
                errsevmsg = errsevmsg_0;
                break;
            case 1:
                err_sev = SEV_WARNING;
                errsevmsg = errsevmsg_1;
                break;
            case 2:
                err_sev = SEV_ERROR;
                errsevmsg = errsevmsg_2;
              break;
            case 3:
                err_sev = SEV_FATAL;
                errsevmsg = errsevmsg_3;
           break;
            default:
                err_sev = SEV_WARNING;
                msg = default_msg;
                errsevmsg = errsevmsg_1;
        break;
        }

        if (errpostex)
            ErrPostEx(err_sev, 0, 0, "%s", msg);

        if (fp)
            fprintf(fp, "\n%s: %s", errsevmsg, msg);
        
        error_return = error_return->next;
    }
    return;
}


TxDfDbInfoPtr LIBCALL 
TxDfDbInfoNew (TxDfDbInfoPtr old)

{
	TxDfDbInfoPtr dbinfo;
	dbinfo = MemNew(sizeof(TxDfDbInfo));
	if (old)
		old->next = dbinfo;
	return dbinfo;
}

TxDfDbInfoPtr LIBCALL 
TxDfDbInfoDestruct (TxDfDbInfoPtr dbinfo)

{
	TxDfDbInfoPtr next;

	if (dbinfo == NULL)
		return NULL;

	while (dbinfo)
	{
		dbinfo->name = MemFree(dbinfo->name);
		dbinfo->definition = MemFree(dbinfo->definition);
		dbinfo->date = MemFree(dbinfo->date);
		next = dbinfo->next;
		dbinfo = MemFree(dbinfo);
		dbinfo = next;
	}

	return dbinfo;
}

Boolean LIBCALL
PrintDbReport(TxDfDbInfoPtr dbinfo, Int4 line_length, FILE *outfp)

{

	if (dbinfo == NULL || outfp == NULL)
		return FALSE;

	asn2ff_set_output(outfp, NULL);

	ff_StartPrint(2, 2, line_length, NULL);

	if (dbinfo->subset == FALSE)
	{
		ff_AddString("Database: ");
		ff_AddString(dbinfo->definition);
		NewContLine();
		ff_AddString("  Posted date:  ");
		ff_AddString(dbinfo->date);
		NewContLine();
		ff_AddString("Number of letters in database: "); 
		ff_AddString(Nlm_Int8tostr((Int8) dbinfo->total_length, 1));
		NewContLine();
		ff_AddString("Number of sequences in database:  ");
		ff_AddString(Ltostr((long) dbinfo->number_seqs, 1));
		NewContLine();
	}
	else
	{
		ff_AddString("Subset of the database(s) listed below");
		NewContLine();
		ff_AddString("   Number of letters searched: "); 
		ff_AddString(Nlm_Int8tostr((Int8) dbinfo->total_length, 1));
		NewContLine();
		ff_AddString("   Number of sequences searched:  ");
		ff_AddString(Ltostr((long) dbinfo->number_seqs, 1));
		NewContLine();
	}
	ff_EndPrint();

	return TRUE;
}

/*
	Prints an acknowledgement of the Blast Query, in the standard
	BLAST format.
*/


Boolean LIBCALL
AcknowledgeBlastQuery(BioseqPtr bsp, Int4 line_length, FILE *outfp, Boolean believe_query, Boolean html)

{
	Char buffer[BUFFER_LENGTH];

	if (bsp == NULL || outfp == NULL)
		return FALSE;
	
	asn2ff_set_output(outfp, NULL);

	ff_StartPrint(0, 0, line_length, NULL);
	if (html)
		ff_AddString("<b>Query=</b> ");
	else
		ff_AddString("Query= ");
	if (bsp->id && believe_query)
	{
		SeqIdWrite(bsp->id, buffer, PRINTID_FASTA_LONG, BUFFER_LENGTH);
		ff_AddString(buffer);
		ff_AddChar(' ');
	}
	ff_AddString(BioseqGetTitle(bsp));
	NewContLine();
	TabToColumn(10);
	ff_AddChar('(');
	ff_AddString(Ltostr((long) BioseqGetLen(bsp), 1));
	ff_AddString(" letters)");
	NewContLine();
        ff_EndPrint();

        return TRUE;
}

/*
	return the version of BLAST as a char. string.
*/
CharPtr LIBCALL
BlastGetReleaseDate (void)

{
	return BLAST_RELEASE_DATE;
}


/*
	return the version of BLAST as a char. string.
*/
CharPtr LIBCALL
BlastGetVersionNumber (void)

{
	return BLAST_ENGINE_VERSION;
}

Boolean BlastPrintVersionInfo (CharPtr program, Boolean html, FILE *outfp)

{
	return BlastPrintVersionInfoEx(program, html, BlastGetVersionNumber(), BlastGetReleaseDate(), outfp);
}

Boolean BlastPrintVersionInfoEx (CharPtr program, Boolean html, CharPtr version, CharPtr date, FILE *outfp)

{
	CharPtr ret_buffer;


	if (outfp == NULL)
		return FALSE;

	ret_buffer = StringSave(program);
	Nlm_StrUpper(ret_buffer);
	if (html)
		fprintf(outfp, "<b>%s %s [%s]</b>\n", ret_buffer, version, date);
	else
		fprintf(outfp, "%s %s [%s]\n", ret_buffer, version, date);
	ret_buffer = MemFree(ret_buffer);

	return TRUE;
}

/* 
	Returns a reference for the header.
	The newlines are represented by tildes, use PrintTildeSepLines
	to print this.
*/

CharPtr LIBCALL
BlastGetReference(Boolean html)

{
	CharPtr ret_buffer;
	Int2 ret_buffer_length;

	ret_buffer = NULL;
	ret_buffer_length = 0;

	
	if (html) {
		add_string_to_bufferEx("<b><a href=\"http://www.ncbi.nlm.nih.gov/htbin-post/Entrez/query?uid=9254694&form=6&db=m&Dopt=r\">Reference</a>:</b>", &ret_buffer, &ret_buffer_length, TRUE);
		add_string_to_bufferEx("Altschul, Stephen F., Thomas L. Madden, Alejandro A. Sch&auml;ffer, ", &ret_buffer, &ret_buffer_length, TRUE);
	} else
		add_string_to_bufferEx("Reference: Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schaffer, ", &ret_buffer, &ret_buffer_length, TRUE);
	add_string_to_bufferEx("Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), ", &ret_buffer, &ret_buffer_length, TRUE);
	add_string_to_bufferEx("\"Gapped BLAST and PSI-BLAST: a new generation of protein database search", &ret_buffer, &ret_buffer_length, TRUE);
	add_string_to_bufferEx("programs\",  Nucleic Acids Res. 25:3389-3402.", &ret_buffer, &ret_buffer_length, TRUE);
	
	return ret_buffer;
}

Boolean LIBCALL
BlastPrintReference(Boolean html, Int4 line_length, FILE *outfp)

{
	CharPtr ret_buffer;
	
	if (outfp == NULL)
		return FALSE;
	
        ret_buffer = BlastGetReference(html);
        PrintTildeSepLines(ret_buffer, line_length, outfp);
        ret_buffer = MemFree(ret_buffer);

	return TRUE;
}

/* 
	Returns a reference for the header.
	The newlines are represented by tildes, use PrintTildeSepLines
	to print this.
*/

CharPtr LIBCALL
BlastGetPhiReference(Boolean html)

{
	CharPtr ret_buffer;
	Int2 ret_buffer_length;

	ret_buffer = NULL;
	ret_buffer_length = 0;

	
	if (html) {
		add_string_to_bufferEx("<b><a href=\"http://www.ncbi.nlm.nih.gov/htbin-post/Entrez/query?uid=9705509&form=6&db=m&Dopt=r\">Reference</a>:</b>", &ret_buffer, &ret_buffer_length, TRUE);
		add_string_to_bufferEx("Zhang, Zheng, Alejandro A. Sch&auml;ffer, Webb Miller, Thomas L. Madden, ", &ret_buffer, &ret_buffer_length, TRUE);
	} else
		add_string_to_bufferEx("Reference: Zhang, Zheng, Alejandro A. Schaffer, Webb Miller, Thomas L. Madden, ", &ret_buffer, &ret_buffer_length, TRUE);
	add_string_to_bufferEx("David J. Lipman, Eugene V. Koonin, and Stephen F. Altschul (1998), ", &ret_buffer, &ret_buffer_length, TRUE);
	add_string_to_bufferEx("\"Protein sequence similarity searches using patterns as seeds\", ", &ret_buffer, &ret_buffer_length, TRUE);
	add_string_to_bufferEx("Nucleic Acids Res. 26:3986-3990.", &ret_buffer, &ret_buffer_length, TRUE);
	
	return ret_buffer;
}

Boolean LIBCALL
BlastPrintPhiReference(Boolean html, Int4 line_length, FILE *outfp)

{
	CharPtr ret_buffer;
	
	if (outfp == NULL)
		return FALSE;
	
        ret_buffer = BlastGetPhiReference(html);
        PrintTildeSepLines(ret_buffer, line_length, outfp);
        ret_buffer = MemFree(ret_buffer);

	return TRUE;
}
CharPtr scan_to_break (CharPtr ptr)
{

	while (*ptr != NULLB)
	{
		if (*ptr == ';')
		{
			ptr++;
			break;
		}
		ptr++;
	}

	return ptr;
}

Boolean LIBCALL
BlastPrintFilterWarning (CharPtr filter_string, Int4 line_length, FILE *outfp, Boolean html)

{
	CharPtr ptr;

	ptr = filter_string;

	if (filter_string == NULL || outfp == NULL)
		return FALSE;

	while (*ptr != NULLB)
	{
		if (*ptr == 'S')
		{
			ptr = scan_to_break(ptr);
		}
		else if (*ptr == 'C')
		{
			ptr = scan_to_break(ptr);
		}
		else if (*ptr == 'D')
		{
			ptr = scan_to_break(ptr);
		}
		else if (*ptr == 'R')
		{
			ptr = scan_to_break(ptr);
			if (html)
				fprintf(outfp, "<B>NOTE:</B>");
			else
				fprintf(outfp, "NOTE:");
			fprintf(outfp, " This query has been filtered for human repeats.\n");
			fprintf(outfp, "This filtering is effective for 70-90%% of all repeats.\n\n");
		}
		else if (*ptr == 'L')
		{ /* do low-complexity filtering; dust for blastn, otherwise seg.*/
			ptr = scan_to_break(ptr);
		}
		else
		{
			ptr++;
		}
	}

	return TRUE;

}

/*
	Initialize the options structure.

	The fields should be set to default values, that depend on the program.
*/
BLAST_OptionsBlkPtr LIBCALL 
BLASTOptionNew(CharPtr progname, Boolean gapped)

{
	BLAST_OptionsBlkPtr options;

	options = (BLAST_OptionsBlkPtr) MemNew(sizeof(BLAST_OptionsBlk));

	options->perform_culling = FALSE;	/* Results should not be culled at all right now. */

	options->program_name = StringSave(progname);
	options->required_start = 0;
	options->required_end = -1;	/* -1 indicates the end of the query. */
	options->cutoff_s = 0;
	options->cutoff_s2 = 0;
	options->db_length = 0;		/* zero means that real size will be used. */
	options->searchsp_eff = 0;	/* zero means that real size will be used. */

	options->block_width = 20;
	options->hsp_range_max = 100;
	options->entrez_query = NULL;
	options->gifile = NULL;
	options->gilist = NULL;

	if (gapped)
	{
		options->gapped_calculation = TRUE;
/* for testing
		options->do_sum_stats = FALSE;
*/
		options->do_sum_stats = TRUE;
	}
	else
	{
		options->gapped_calculation = FALSE;
		options->do_sum_stats = TRUE;
		options->hsp_num_max = 100;
	}

	options->discontinuous = FALSE;	/* discontinuous is default. */
	if (StringICmp(progname, "blastn") == 0)
	{
		options->gap_decay_rate = 0.5;
		options->gap_prob = 0.5;
		options->gap_size = 50;
		options->window_size = 0;
		options->threshold_first = 0;
		options->threshold_second = 0;
		options->expect_value  = 10;
		options->hitlist_size = 500;
		options->two_pass_method  = FALSE;
		options->multiple_hits_only  = FALSE;
		options->number_of_bits  = 0.0;
		/* 1st pass not done for blastn. */
		options->dropoff_2nd_pass  = 20;
		options->matrix  = NULL;
		options->old_stats  = FALSE;
		options->wordsize  = 11;
		options->penalty  = -3;
		options->reward  = 1;
		options->e2 = 0.05;
		/* Used in the post-process gapping of the blastn result. */
		options->gap_open  = 5;
		options->gap_extend  = 2;
                options->decline_align = INT2_MAX;
		options->gap_x_dropoff  = 20;
		options->gap_x_dropoff_final  = 50;
		options->gap_trigger  = 25.0;
		options->strand_option  = BLAST_BOTH_STRAND;
		options->no_check_score  = FALSE;
	}
	else
	{
		options->gap_size = 50;
		options->window_size = 40;
		options->expect_value  = 10;
		options->hitlist_size = 500;
		options->two_pass_method  = TRUE;
		options->multiple_hits_only  = FALSE;
		options->number_of_bits  = 0.0;
		options->dropoff_1st_pass  = 7;
		options->dropoff_2nd_pass  = 10;
		options->matrix  = StringSave("BLOSUM62");
		options->old_stats  = FALSE;
		options->wordsize  = 3;
		options->penalty  = 0;
		options->reward  = 0;
		options->gap_decay_rate = 0.5;
		options->gap_prob = 0.5;
		options->no_check_score  = TRUE;
		if (gapped)
		{
			options->two_pass_method = FALSE;
			options->multiple_hits_only  = TRUE;
			options->dropoff_2nd_pass  = options->dropoff_1st_pass;
			options->gap_decay_rate = 0.1;
			options->gap_prob = 1.0;
		}

		options->gap_open  = 11;
		options->gap_extend  = 1;
                options->decline_align = INT2_MAX;
		options->gap_x_dropoff  = 15;
		options->gap_x_dropoff_final  = 25;
		options->gap_trigger  = 22.0;

		if (StringICmp(progname, "blastp") == 0)
		{
			options->e2 = 0.5;
			options->threshold_first = 11;
			options->threshold_second = 11;
		}
		else if (StringICmp(progname, "blastx") == 0)
		{
			options->e2 = 0.25;
			options->genetic_code = 1;
			options->threshold_first = 12;
			options->threshold_second = 12;
			options->do_sum_stats = TRUE;
			options->strand_option  = BLAST_BOTH_STRAND;
		}
		else if (StringICmp(progname, "tblastn") == 0)
		{
			options->e2 = 0.15;
			options->db_genetic_code = 1;
			options->threshold_first = 13;
			options->threshold_second = 13;
			options->do_sum_stats = TRUE;
		}
		else if (StringICmp(progname, "tblastx") == 0)
		{
			options->e2 = 0.1;
			options->genetic_code = 1;
			options->db_genetic_code = 1;
			options->threshold_first = 13;
			options->threshold_second = 13;
			options->gap_open  = 0;
			options->gap_extend  = 0;
			options->gap_x_dropoff  = 0;
			options->gap_x_dropoff_final  = 0;
			options->gapped_calculation = FALSE;
			options->do_sum_stats = TRUE;
			options->strand_option  = BLAST_BOTH_STRAND;
			options->hsp_num_max = 100;
		}
	}

	return options;
}

/*
	Delete the Options structure.
*/
BLAST_OptionsBlkPtr LIBCALL 
BLASTOptionDelete(BLAST_OptionsBlkPtr options)
{
    if (options == NULL)
        return NULL;
    
    if (options->matrix != NULL)
        MemFree(options->matrix);
    
    if (options->program_name != NULL)
        MemFree(options->program_name);
    
    MemFree(options->filter_string);
    
    options = MemFree(options);
    return options;
}


/*
	Validate the Options structure.  If an invalid option is found,
	call BLASTOptionDelete and issue an error message.
*/
BLAST_OptionsBlkPtr LIBCALL 
BLASTOptionValidate(BLAST_OptionsBlkPtr options, CharPtr progname)

{
	Int2 status;
	ValNodePtr error_return=NULL;

	status = BLASTOptionValidateEx(options, progname, &error_return);

	if (status != 0)
		options = BLASTOptionDelete(options);
	
	BlastErrorPrint(error_return);
	BlastErrorChainDestroy(error_return);

	return options;
}

/*
	Validate the Options structure.  If an invalid option is found,
	call BLASTOptionDelete and issue an error message.
*/
Int2 LIBCALL 
BLASTOptionValidateEx (BLAST_OptionsBlkPtr options, CharPtr progname, ValNodePtr PNTR error_return)

{
	Int2 status=0;

	if (options->hitlist_size < 1)
	{
		BlastConstructErrorMessage("BLASTOptionValidateEx", "No hits are being saved", 1, error_return);
		return 1;
	}

	if (options->expect_value <= 0.0 && options->cutoff_s == 0)
	{
		BlastConstructErrorMessage("BLASTOptionValidateEx", "expect value must be greater than zero", 1, error_return);
		return 1;
	}

	if (options->wordsize <= 0)
	{
		BlastConstructErrorMessage("BLASTOptionValidateEx", "wordsize must be non-zero", 1, error_return);
		return 1;
	}

	if (StringICmp(progname, "blastn") == 0)
	{
               if (options->penalty >= 0)
               {
                       BlastConstructErrorMessage("BLASTOptionValidateEx", "BLASTN penalty must be negative", 1, error_return);
                       return 1;
               }

		if (options->wordsize < 7)
		{
			BlastConstructErrorMessage("BLASTOptionValidateEx", "Wordsize must be 7 or greater", 1, error_return);
			return 1;
		}
		if (options->threshold_first != 0 || options->threshold_second != 0)
		{
			BlastConstructErrorMessage("BLASTOptionValidateEx", "non-zero threshold not permitted with blastn", 1, error_return);
			return 1;
		}

		if (options->two_pass_method == TRUE || options->number_of_bits != 0.0)
		{
			BlastConstructErrorMessage("BLASTOptionValidateEx", "two-passes not available for blastn", 1, error_return);
			return 1;
		}

		if ((options->strand_option | BLAST_BOTH_STRAND) == 0)
		{
			BlastConstructErrorMessage("BLASTOptionValidateEx", "invalid strand specified", 1, error_return);
			return 1;
		}

/*
		if (options->multiple_hits_only == TRUE)
		{
			BlastConstructErrorMessage("BLASTOptionValidateEx", "multiple hits not available for blastn", 1, error_return);
			return 1;
		}
*/
	}
	else
	{
		if (StringICmp(progname, "blastx") == 0 || StringICmp(progname, "tblastx") == 0)
		{
			if ((options->strand_option | BLAST_BOTH_STRAND) == 0)
			{
				BlastConstructErrorMessage("BLASTOptionValidateEx", "invalid strand specified", 1, error_return);
				return 1;
			}
		}
		if (options->wordsize < 2 || options->wordsize > 3)
		{
			BlastConstructErrorMessage("BLASTOptionValidateEx", "Valid wordsize range is 2 to 3", 1, error_return);
			return 1;
		}
		if (options->threshold_second == 0)
		{
			BlastConstructErrorMessage("BLASTOptionValidateEx", "non-zero threshold required", 1, error_return );
			return 1;
		}
		if (options->penalty != 0 || options->reward != 0)
		{
			BlastConstructErrorMessage("BLASTOptionValidateEx", "penalty or reward can only be non-zero for blastn", 1, error_return);
			return 1;
		}

		if (StringICmp(progname, "tblastx") == 0)
		{
			if (options->gapped_calculation == TRUE)
			{
				BlastConstructErrorMessage("BLASTOptionValidateEx", "gapped calculations not available with tblastx", 1, error_return);
				return 1;
			}
		}
		
		if (options->gapped_calculation == TRUE)
		{
			status = BlastKarlinBlkGappedCalc(NULL, options->gap_open, options->gap_extend, options->matrix, error_return);
		}
	}

        if (options->isPatternSearch && (!options->gapped_calculation))
	    {
		BlastConstructErrorMessage("BLASTOptionValidateEx", "PHI-BLAST cannot use ungapped alignments", 1, error_return);
		return 1;
	    }

        if (options->tweak_parameters && (!options->gapped_calculation))
	    {
		BlastConstructErrorMessage("BLASTOptionValidateEx", "parameter adjustment not supported with ungapped alignments", 1, error_return);
		return 1;
	    }

        if (options->smith_waterman && (!options->gapped_calculation))
	    {
		BlastConstructErrorMessage("BLASTOptionValidateEx", "locally optimal alignments not supported with ungapped alignments", 1, error_return);
		return 1;
	    }

	return status;
}

/*
	Changes the matrix value to the one given and sets the 
	default parameters for that Matrix. 
*/

Int2 LIBCALL 
BLASTOptionSetGapParams (BLAST_OptionsBlkPtr options, CharPtr matrix_name, Int4 open, Int4 extended)

{
	Boolean found_matrix=FALSE, threshold_set=FALSE;

	if (options == NULL || matrix_name == NULL)
		return -1;

	/* blastn is different. */
	if (StringICmp("blastn", options->program_name) == 0)
	{
		options->gap_open  = 5;
		options->gap_extend  = 2;
		return 0;
	}

	if (StringICmp(matrix_name, "BLOSUM62") == 0)
	{
		options->gap_open  = 11;
		options->gap_extend  = 1;
		options->window_size = 40;
		options->threshold_first = 11;
		options->threshold_second = 11;
		found_matrix = TRUE;
		threshold_set = TRUE;
	}
	else if (StringICmp(matrix_name, "BLOSUM45") == 0)
	{
		options->gap_open  = 14;
		options->gap_extend  = 2;
		options->window_size = 60;
		options->threshold_first = 14;
		options->threshold_second = 14;
		found_matrix = TRUE;
		threshold_set = TRUE;
	}
	else if (StringICmp(matrix_name, "BLOSUM50") == 0)
	{
		options->gap_open  = 13;
		options->gap_extend  = 2;
		found_matrix = TRUE;
	}
	else if (StringICmp(matrix_name, "PAM250") == 0)
	{
		options->gap_open  = 14;
		options->gap_extend  = 2;
		found_matrix = TRUE;
	}
	else if (StringICmp(matrix_name, "BLOSUM62_20") == 0)
	{
		options->gap_open  = 11;
		options->gap_extend  = 1;
		found_matrix = TRUE;
	}
	else if (StringICmp(matrix_name, "BLOSUM90") == 0)
	{
		options->gap_open  = 10;
		options->gap_extend  = 1;
		found_matrix = TRUE;
	}
	else if (StringICmp(matrix_name, "BLOSUM80") == 0)
	{
		options->gap_open  = 10;
		options->gap_extend  = 1;
		options->window_size = 25;
		options->threshold_first = 12;
		options->threshold_second = 12;
		found_matrix = TRUE;
		threshold_set = TRUE;
	}
	else if (StringICmp(matrix_name, "PAM30") == 0)
	{
		options->gap_open  = 9;
		options->gap_extend  = 1;
		options->window_size = 15;
		options->threshold_first = 16;
		options->threshold_second = 16;
		found_matrix = TRUE;
		threshold_set = TRUE;
	}
	else if (StringICmp(matrix_name, "PAM70") == 0)
	{
		options->gap_open  = 10;
		options->gap_extend  = 1;
		options->window_size = 20;
		options->threshold_first = 14;
		options->threshold_second = 14;
		found_matrix = TRUE;
		threshold_set = TRUE;
	}

	if (open)
	    options->gap_open  = open;
	if (extended)
	    options->gap_extend  = extended;

	if (matrix_name)
	{
		if (options->matrix)
			MemFree(options->matrix);
		options->matrix = StringSave(matrix_name);
	}
	
	if (!found_matrix)
		return -1;

	if (threshold_set)
	{
		if (StringICmp(options->program_name, "blastx") == 0)
		{
			options->threshold_first++;
			options->threshold_second++;
		}
		else if (StringICmp(options->program_name, "tblastn") == 0 || StringICmp(options->program_name, "tblastx") == 0)
		{
			options->threshold_first += 2;
			options->threshold_second += 2;
		}
	}
	return 0;
}

/*
	This function obtains the sequence from a BioseqPtr in ASCII alphabet.
	The return value is a Uint1Ptr containing the sequence, the Int4Ptr
	length inidcates the length of the seqeunce.
*/	

Uint1Ptr
BlastGetSequenceFromBioseq (BioseqPtr bsp, Int4Ptr length)

{
	Int4 index;
	SeqPortPtr spp;
	Uint1 residue;
	Uint1Ptr sequence;

	*length = 0;
	if (bsp == NULL)
		return NULL;

	sequence = MemNew((1+bsp->length)*sizeof(Uint1));
	if (sequence == NULL)
		return NULL;

        if (ISA_na(bsp->mol))
        	spp = SeqPortNew(bsp, 0, -1, Seq_strand_plus, Seq_code_iupacna);
	else
        	spp = SeqPortNew(bsp, 0, -1, Seq_strand_unknown, Seq_code_ncbieaa);

	index=0;
	while ((residue=SeqPortGetResidue(spp)) != SEQPORT_EOF)
	{
		if (residue == SEQPORT_VIRT)
			continue;
		sequence[index] = residue;
		index++;
	}
	sequence[index] = NULLB;
	spp = SeqPortFree(spp);

	*length = index;
	return sequence;
}


/*
	Adjusts the length and number of sequences in a database according
	to the gi_list or seqIdPtr list given.
*/

Boolean
BlastAdjustDbNumbers (ReadDBFILEPtr rdfp, Int8Ptr db_length, Int4Ptr db_number, SeqIdPtr seqid_list, BlastDoubleInt4Ptr gi_list, OIDListPtr oidlist, BlastDoubleInt4Ptr PNTR gi_list_pointers, Int4 gi_list_total)

{
	Int4 count, db_number_start, index, ordinal_id;
        Int8	db_length_private, db_length_start;
	SeqIdPtr sip;

		
	count = 0;
	db_length_private = 0;

	if (seqid_list)
	{
		sip = seqid_list;
		while (sip)
		{
			ordinal_id = SeqId2OrdinalId(rdfp, sip);
			if (ordinal_id > 0)
				count++;
			sip = sip->next;
		}
	}
	else if (oidlist) {

	    Uint4	mask, base=0, maskindex, i;
	    Uint4	total_mask = oidlist->total/MASK_WORD_SIZE + 1;

	    maskindex = 0;

	    while (maskindex < total_mask) {
		/* for each long-word mask */
		mask = Nlm_SwapUint4(oidlist->list[maskindex]);
		i = 0;
		while (mask) {
		    if (mask & (((Uint4)0x1)<<(MASK_WORD_SIZE-1))) {
			count++;
			db_length_private +=
			    readdb_get_sequence_length(rdfp, base + i);
		    }
		    mask <<= 1;
		    i++;
		}
		maskindex++;
		base += MASK_WORD_SIZE;
	    }
	}

	*db_length = db_length_private;
	*db_number = count;

	return TRUE;
}

/*
	Deletes only the BlastGiListPtr, not
	the associated arrays.
*/
BlastGiListPtr 
BlastGiListDestruct(BlastGiListPtr blast_gi_list, Boolean contents)

{
	if (blast_gi_list == NULL)
		return NULL;

	if (contents)
	{ /* On the main thread.  Deallocate contents. */
		MemFree(blast_gi_list->gi_list);
		MemFree(blast_gi_list->gi_list_pointer);
	}

	return MemFree(blast_gi_list);
}

/*
	Allocates BlastGiListPtr.  The caller still owns
	the gi_list and must delete it.
*/

BlastGiListPtr 
BlastGiListNew(BlastDoubleInt4Ptr gi_list, BlastDoubleInt4Ptr PNTR gi_list_pointers, Int4 total)

{
	BlastGiListPtr blast_gi_list;

	if (gi_list == NULL || total == 0)
		return NULL;

	blast_gi_list = MemNew(sizeof(BlastGiList));
	blast_gi_list->gi_list = gi_list;
	blast_gi_list->gi_list_pointer = gi_list_pointers;
	blast_gi_list->total = total;
	blast_gi_list->current = 0;
	
	return blast_gi_list;
}


#define POSIT_PERCENT 0.05
#define POSIT_NUM_ITERATIONS 10

#define POSIT_SCALE_FACTOR 1000
#define Xchar   21    /*character for low-complexity columns*/


/* This should be more than large enough for any alphabet. */
#define POSIT_MAX_ALPHABET_SIZE 30

/*Compute probabilities for each score in posMatrix,
also sets minScore and maxScore*/
static BLAST_ScoreFreqPtr 
BlastComputeProbs(BlastMatrixRescalePtr matrix_rescale, Boolean position_dependent)
{
   Int2 alphabet_total;
   Int4 c;  /*index on characters*/
   Int4 p;  /*index on positions*/
   Int4 s;  /*index on scores */
   Int4 dim1, dim2;
   BLAST_Score score_min, score_max;
   BLAST_ScoreFreqPtr sfp;
   Int4 numberOfScores; /* number of distinct scores*/
   Int4 score;  /*one score in the matrix*/
   Nlm_FloatHi increment;  /*Increment in probability due to one score*/
   Int4 effectiveLength;
   Int4Ptr *matrix;
   Uint1 std_alphabet[POSIT_MAX_ALPHABET_SIZE];

   alphabet_total = BlastGetStdAlphabet(Seq_code_ncbistdaa, std_alphabet, POSIT_MAX_ALPHABET_SIZE);
   if (alphabet_total <= 0)
	return NULL;
    
   matrix = matrix_rescale->matrix;
   score_min = 0;
   score_max = 0;
   if(position_dependent)
   {
	dim1 = matrix_rescale->query_length;
	dim2 = alphabet_total;
   }
   else
   {
	dim1 = alphabet_total;
	dim2 = alphabet_total;
   }

   effectiveLength = 0;
   for (p = 0; p < matrix_rescale->query_length; p++)
     if (Xchar != matrix_rescale->query[p])
       effectiveLength++;
   for (p = 0; p < dim1; p++)
     if (Xchar != matrix_rescale->query[p])
       for (c = 0; c < dim2; c++) {
	 if (matrix[p][std_alphabet[c]] <= BLAST_SCORE_MIN || matrix[p][std_alphabet[c]] >= BLAST_SCORE_MAX)
		continue;
	 if (matrix[p][std_alphabet[c]] < (score_min))
	   (score_min) = matrix[p][std_alphabet[c]];
	 if (matrix[p][std_alphabet[c]] > (score_max))
	   (score_max) = matrix[p][std_alphabet[c]];
       }
   sfp = BlastScoreFreqNew(score_min, score_max);
   sfp->obs_min = sfp->score_min;
   sfp->obs_max = sfp->score_max;
   numberOfScores = (sfp->score_max) - (sfp->score_min) + 1;
   for (p = 0; p < dim1; p++)
     if (Xchar != matrix_rescale->query[p])
       for (c = 0; c < dim2; c++) {
	 /*Increment the weight for the score in position [p][std_alphabet[c]] */
	 score = matrix[p][std_alphabet[c]];
	 increment =
	   (matrix_rescale->standardProb[std_alphabet[c]]/ effectiveLength);
	 sfp->sprob[score]+= increment;
       }

   sfp->score_avg = 0.0;
   for(s = sfp->score_min; s <= sfp->score_max; s++)
     sfp->score_avg += s*sfp->sprob[s];
   return(sfp);
}

void LIBCALL
updateLambdaK(BlastMatrixRescalePtr matrix_rescale, Boolean position_dependent)
{
  BLAST_ScoreFreqPtr sfp;
 
  sfp = BlastComputeProbs(matrix_rescale, position_dependent);
  BlastKarlinBlkCalc(matrix_rescale->kbp_psi[0], sfp);
  matrix_rescale->kbp_gap_psi[0]->K = (matrix_rescale->kbp_psi[0]->K)*(matrix_rescale->kbp_gap_std[0]->K)/matrix_rescale->K_ideal;
  matrix_rescale->kbp_gap_psi[0]->logK = log(matrix_rescale->kbp_gap_psi[0]->K);
  sfp = BlastScoreFreqDestruct(sfp);
}

Nlm_FloatHi 
BlastScaleMatrix(BlastMatrixRescalePtr matrix_rescale, Boolean position_dependent)
{
   Int4 dim1, dim2;
   Int4 a,c; /*loop indices*/
   Boolean too_high=TRUE, done, first_time;
   Nlm_FloatHi factor, factor_low=1.0, factor_high=1.0;
   Nlm_FloatHi lambda, new_lambda; /*Karlin-Altschul parameter*/
   Int4 index; /*loop index for binary search*/
   Int4Ptr *private_matrix;
   Int4Ptr *matrix;

   private_matrix = matrix_rescale->private_matrix;
   matrix = matrix_rescale->matrix;


/* Bracket the values. */
   if (position_dependent)
   {
   	dim1 = matrix_rescale->query_length;
   	dim2 = matrix_rescale->alphabet_size;
   }
   else
   {
   	dim1 = matrix_rescale->alphabet_size;
   	dim2 = matrix_rescale->alphabet_size;
   }
   lambda = matrix_rescale->lambda_ideal;

   done = FALSE;
   first_time = TRUE;
   factor = 1.0;
   while (done != TRUE)
   {
   	for(c = 0; c < dim1; c++)
   	{
       	    for(a = 0; a < dim2; a++)
	    {
		if (private_matrix[c][a] == BLAST_SCORE_MIN)
		{
			matrix[c][a] = BLAST_SCORE_MIN;
		}
		else
		{
			matrix[c][a] = (factor*private_matrix[c][a])/POSIT_SCALE_FACTOR;
		}
	    }
        }

        updateLambdaK(matrix_rescale, position_dependent);
	new_lambda = matrix_rescale->kbp_psi[0]->Lambda;
	if (new_lambda > lambda)
	{
		if (first_time)
		{
			factor_high = 1.0 + POSIT_PERCENT;
			factor = factor_high;
			factor_low = 1.0;
			too_high = TRUE;
			first_time = FALSE;
		}
		else
		{
			if (too_high == FALSE)
				break;
			factor_high += (factor_high-1.0);
			factor = factor_high;
		}
	}
	else 
	{
		if (first_time)
		{
			factor_high = 1.0;
			factor_low = 1.0 - POSIT_PERCENT;
			factor = factor_low;
			too_high = FALSE;
			first_time = FALSE;
		}
		else
		{
			if (too_high == TRUE)
				break;
			factor_low += (factor_low-1.0);
			factor = factor_low;
		}
	}
   }

/* binary search for ten times. */
   for (index=0; index<POSIT_NUM_ITERATIONS; index++)
   {
        factor = 0.5*(factor_high+factor_low);
   	for(c = 0; c < dim1; c++)
   	{
       	    for(a = 0; a < dim2; a++)
	    {
		if (private_matrix[c][a] == BLAST_SCORE_MIN)
		{
			matrix[c][a] = BLAST_SCORE_MIN;
		}
		else
		{
			matrix[c][a] = (factor*private_matrix[c][a])/POSIT_SCALE_FACTOR;
		}
	    }
   	}
        updateLambdaK(matrix_rescale, position_dependent);
	new_lambda = matrix_rescale->kbp_psi[0]->Lambda;
	if (new_lambda > lambda)
	{
		factor_low = factor;
	}
	else
	{
		factor_high = factor;
	}
    }

   for(a = 0; a < matrix_rescale->alphabet_size; a++) 
     matrix[dim1][a] = BLAST_SCORE_MIN;

	return factor;
}


/*
	Deallocates only memory for BlastMatrixRescalePtr.
*/
BlastMatrixRescalePtr
BlastMatrixRescaleDestruct (BlastMatrixRescalePtr matrix_rescale)

{
	if (matrix_rescale == NULL)
		return NULL;

	return MemFree(matrix_rescale);
}

/*
	Allocates and fills the BlastMatrixRescalePtr.
*/

BlastMatrixRescalePtr
BlastMatrixRescaleNew(Int4 alphabet_size, Int4 query_length, Uint1Ptr query,  Nlm_FloatHiPtr standardProb, Int4Ptr *matrix, Int4Ptr *private_matrix, BLAST_KarlinBlkPtr *kbp_std, BLAST_KarlinBlkPtr *kbp_psi, BLAST_KarlinBlkPtr *kbp_gap_std, BLAST_KarlinBlkPtr *kbp_gap_psi, Nlm_FloatHi lambda_ideal,  Nlm_FloatHi K_ideal)

{
	BlastMatrixRescalePtr matrix_rescale;

	matrix_rescale = (BlastMatrixRescalePtr) MemNew(sizeof(BlastMatrixRescale));

	matrix_rescale->alphabet_size = alphabet_size;
	matrix_rescale->query_length = query_length;
	matrix_rescale->query = query;
	matrix_rescale->standardProb = standardProb;
	matrix_rescale->matrix = matrix;
	matrix_rescale->private_matrix = private_matrix;
	matrix_rescale->kbp_std = kbp_std;
	matrix_rescale->kbp_psi = kbp_psi;
	matrix_rescale->kbp_gap_std = kbp_gap_std;
	matrix_rescale->kbp_gap_psi = kbp_gap_psi;
	matrix_rescale->lambda_ideal = lambda_ideal;
	matrix_rescale->K_ideal = K_ideal;

	return matrix_rescale;
}

/*
	BlastSeqLocCount counts the number of SeqLoc's 
*/
static Int4 
BlastSeqLocCount (SeqLocPtr mask_slp)

{
        Int4 index=0;
	SeqLocPtr slp;
       
	while (mask_slp)
	{
		slp=NULL;
        	while((slp = SeqLocFindNext(mask_slp, slp))!=NULL)
        	{
			index++;
        	}
		mask_slp = mask_slp->next;
	}

	return index;
}

static Boolean 
BlastIntervalSort (BlastDoubleInt4Ptr *link_ptr, Int4 link_value, Int4 total, int (LIBCALLBACK *callback )PROTO ((Nlm_VoidPtr, Int4 start, Int4 stop)), VoidPtr ptr)

{
	Boolean do_callback, start=TRUE;
	Int4 index, start_pos, stop_pos, largest_stop_pos;

	index=0;
	while (index < total)
	{
		if (start == TRUE)
		{
			start_pos = link_ptr[index]->gi;
			start = FALSE;
			largest_stop_pos = 0;
		}
		else
		{
			/* Keep track of largest stop position. */
			largest_stop_pos = MAX(largest_stop_pos, link_ptr[index]->ordinal_id);
			do_callback = FALSE;
			if (index == total-1)	/* Last one. */
			{
				stop_pos = link_ptr[index]->ordinal_id;
				start = TRUE;
				do_callback = TRUE;
			}
			else if (largest_stop_pos+link_value < link_ptr[index+1]->gi)
			{ /* Check overlap with next one. */
				stop_pos = link_ptr[index]->ordinal_id;
				start = TRUE;
				do_callback = TRUE;
			}
			
			if (do_callback)
			{
				callback(ptr, start_pos, MAX(largest_stop_pos, stop_pos));
			}
			index++;
		}
	}

	return TRUE;
}

static int LIBCALLBACK
list_ptr_compare(VoidPtr v1, VoidPtr v2)

{
	BlastDoubleInt4Ptr h1, h2;
	BlastDoubleInt4Ptr *hp1, *hp2;

	hp1 = (BlastDoubleInt4Ptr PNTR) v1;
	hp2 = (BlastDoubleInt4Ptr PNTR) v2;
	h1 = *hp1;
	h2 = *hp2;

	if (h1->gi < h2->gi)
		return -1;
	if (h1->gi > h2->gi)
		return 1;

	return 0;
}

int LIBCALLBACK
slp_callback(VoidPtr ptr, Int4 start, Int4 stop)

{
	ValNodePtr *vnpp;


	vnpp = (ValNodePtr PNTR) ptr;	

	ValNodeAddInt(vnpp, 0, start);
	ValNodeAddInt(vnpp, 1, stop);

	return 1;
}


ValNodePtr
BlastSeqLocFillDoubleInt (SeqLocPtr mask_slp, Int4 max_length, Boolean reverse)

{
	Int4 count, index, start, stop;
	BlastDoubleInt4Ptr list_pri, *list_ptr_pri;
	SeqLocPtr slp;
	ValNodePtr vnp;

	vnp = NULL;
	if (mask_slp == NULL)
	{
		ValNodeAddInt(&vnp, 1, -1);
		ValNodeAddInt(&vnp, 0, max_length);
		return vnp;
	}

	count = BlastSeqLocCount (mask_slp);
	list_pri = (BlastDoubleInt4Ptr) MemNew(count*sizeof(BlastDoubleInt4)); 
	list_ptr_pri = (BlastDoubleInt4Ptr PNTR) MemNew(count*sizeof(BlastDoubleInt4Ptr)); 

	index=0;
	while (mask_slp)
	{
		slp=NULL;
        	while((slp = SeqLocFindNext(mask_slp, slp))!=NULL)
        	{
			if (reverse)
			{
				start = max_length - 1 - SeqLocStop(slp);
				stop = max_length - 1 - SeqLocStart(slp);
			}
			else
			{
              			start = SeqLocStart(slp);
              			stop = SeqLocStop(slp);
			}

			list_pri[index].gi = start;
			list_pri[index].ordinal_id = stop;
			list_ptr_pri[index] = &(list_pri[index]);
			index++;
        	}
		mask_slp = mask_slp->next;
	}

	HeapSort(list_ptr_pri, count, sizeof(BlastHitRangePtr PNTR), list_ptr_compare);

	/* Allows the proper start. */
	ValNodeAddInt(&vnp, 1, -1);
	BlastIntervalSort(list_ptr_pri, 0, count, slp_callback, (VoidPtr) &vnp);

	return vnp;
}



static Boolean
small(BLASTResultHspPtr a, BLASTResultHspPtr b)
{
  if (a->e_value > b->e_value) return TRUE;
  if (a->e_value < b->e_value) return FALSE;
  if (a->score > b->score) return FALSE;
  if (a->score < b->score) return TRUE;
  if (a->point_back->subject_id > b->point_back->subject_id) return FALSE;
  return TRUE;
}


/*
	this is some sort of HeapSort (or it makes the heap
	as it is used in HeapSort).
*/
static void
BlastHeapify(BLASTHeapPtr which_heap, Int4 position)
{
  Int4 heap_size, index, lim, left_son, small_son;
  BLASTResultHspPtr tmp, PNTR heap;

  heap_size = which_heap->num_in_heap;
  heap = which_heap->heap;
  index = position; lim = heap_size/2;

  while (index < lim) {
    left_son = 2*index + 1;
    if (left_son == heap_size-1)
      small_son = left_son;
    else {
      if (small(heap[left_son],heap[left_son+1])) small_son = left_son;
      else small_son = left_son+1;
    }
/* If heap[small_son] is less significant than heap[index], then
   switch them.  Otherwise exit the loop. 
*/
    if (small(heap[small_son], heap[index])) {
      tmp = heap[index];
      heap[index] = heap[small_son];
      heap[small_son] = tmp;
      index = small_son;
    } else
      break;
  }
}

/*
	Insert a new element into the BLASTHeapPtr list, which
	appears to be ordered (here) from least to most significant.
	the order does not seem to be strictly from least to
	most significant, or is it?
*/
static void
BlastInsertHeap(BLASTHeapPtr which_heap, BLASTResultHspPtr hp)
{
  Int4 heap_size, index, father;
  BLASTResultHspPtr PNTR heap;

  hp->point_back->num_ref+= 1;
  /* Increase heap size by one for new entry. */
  heap_size = (which_heap->num_in_heap+=1);
  heap = which_heap->heap;
  index = heap_size-1; 
 
  while (index > 0) {
    father = (index-1)/2;
    /* If heap[father] is LESS significant than hp, then exit loop. */
    if (small(heap[father], hp)) {
      break;
    }
    /* The 1st time this is called heap[index] is NULL as it's one past
	the last filled in element.  The idea here is obviously to
	move the more significant elements to the end, but don't
	they get out of order then?  
    */
    heap[index] = heap[father];
    index = father;
  }
  /* Fill in new element. */
  heap[index] = hp;
} 

static Int4
BlastDeleteHeap(BLASTHeapPtr which_heap, Int4 position)
{
  Int4 last, return_value;
  BLASTResultHspPtr PNTR heap;

  last = (which_heap->num_in_heap-=1);
  heap = which_heap->heap;
  return_value = (heap[position]->point_back->num_ref -= 1);
  if (position != last) {
    heap[position] = heap[last];
    BlastHeapify(which_heap, position);
  }
  return return_value;
}

static void
BlastInsertWholeHeap(BlastSearchBlkPtr search, BLASTHeapPtr which_heap, Int4 cutvalue)
{
  BLASTHeapPtr hp;
  Int4 i;

  hp = (BLASTHeapPtr) MemNew(sizeof(BLASTHeapStruct));
  hp->heap = (BLASTResultHspPtr PNTR) MemNew(sizeof(BLASTResultHspPtr)*search->pbp->hsp_range_max);
  hp->next = which_heap;
  hp->num_of_ref = 0;
  hp->cutvalue = cutvalue;
  if (which_heap->prev) {
    hp->prev = which_heap->prev;
    which_heap->prev->next = hp;
  } else {
    hp->prev = NULL;
    search->result_struct->heap_ptr = hp;
  }
  which_heap->prev = hp;
  hp->num_in_heap = which_heap->num_in_heap;
  for (i = 0; i < which_heap->num_in_heap; i++) {
    hp->heap[i] = which_heap->heap[i];
    which_heap->heap[i]->point_back->num_ref +=1;
  }
}

static void
BlastDeleteWholeHeap(BlastSearchBlkPtr search, BLASTHeapPtr which_heap)
{
  Int4 i;
  if ((which_heap->num_of_ref -= 1) >0) return;
  if (which_heap->next) which_heap->next->prev = which_heap->prev;
  if (which_heap->prev) which_heap->prev->next = which_heap->next;
  else search->result_struct->heap_ptr = which_heap->next;
  for (i = 0; i < which_heap->num_in_heap; i++) {
    which_heap->heap[i]->point_back->num_ref -= 1;
  }
  which_heap->heap = MemFree(which_heap->heap);
  which_heap = MemFree(which_heap);
}

static Int2
BlastPossibleDeleteWholeHeap(BlastSearchBlkPtr search, BLASTHeapPtr PNTR hhp,  BLASTResultHspPtr heap0)
     /* if the deleted hit result a remove of the right end point, 
	return 1 to make the program rerun 
	*/
{
  BLASTHeapPtr hp = *hhp;

  if (heap0 == NULL || hp == NULL)
	return 0;

  if (heap0->back_left && hp->prev && heap0->back_left == hp->prev) {
    heap0->back_left = NULL;
    if ((hp->prev->num_of_ref -= 1) == 0) {
      BlastDeleteWholeHeap(search, hp->prev);
    }
  }
  if (heap0->back_right == hp) {
    heap0->back_right = NULL;
    if ((hp->num_of_ref -= 1)==0) {
      *hhp = hp->next;    
      BlastDeleteWholeHeap(search, hp);
      *hhp = (*hhp)->prev;
      return 1;
    }
  }
  return 0;
}



Int2
BlastInsertList2Heap(BlastSearchBlkPtr search, BLASTResultHitlistPtr result_hitlist)
{
    BLASTResultHspPtr PNTR heap, hsp;
    BLASTResultHspPtr hsp_array;
    Int4 index, hsp_range_max;
    Int4 begin, end, hspcnt;
    Boolean hsp_deleted, new_inserted;
    BLASTHeapPtr hp;

    if (search->pbp->perform_culling == FALSE)  /* Culling is turned off. */
        return 3;

    hsp_deleted = new_inserted  = FALSE;  
    hsp_range_max = search->pbp->hsp_range_max;
    hspcnt = result_hitlist->hspcnt;
    hsp_array = result_hitlist->hsp_array;

    for (index = 0; index < hspcnt; index++) {
      hsp = &hsp_array[index];    
      begin = hsp->query_offset;
      end = (hsp->query_offset+hsp->query_length-1);
      for (hp = search->result_struct->heap_ptr; hp; hp = hp->next) 
	if (hp->cutvalue >= begin) break;
      if (hp->num_in_heap < hsp_range_max || small(hp->heap[0], hsp)) {
	if (!hp->prev || hp->prev->cutvalue != begin-1) {
	  BlastInsertWholeHeap(search, hp, begin-1);
	}
	hp->prev->num_of_ref +=1;
	hsp->back_left = hp->prev;
      } else hsp->back_left = NULL;
      for (; hp; hp = hp->next) {
	if (end <= hp->cutvalue) break;
	heap = hp->heap;
	if (hp->num_in_heap >= hsp_range_max) {
	  if (small(heap[0], hsp)) {
	    if (BlastPossibleDeleteWholeHeap(search, &hp, heap[0])) continue;
	    if ((heap[0]->point_back->num_ref-=1)==0) 
	      hsp_deleted = TRUE;
	    heap = hp->heap;
	    heap[0] = hsp;
	    hsp->point_back->num_ref +=1;
	    BlastHeapify(hp, 0);
	    new_inserted = TRUE;
	  } 
	} else {
	  BlastInsertHeap(hp, hsp);
	  new_inserted = TRUE;
	} 
      }
      if (hp->num_in_heap < hsp_range_max || small(hp->heap[0], hsp)) {
	if (end != hp->cutvalue) {
	  BlastInsertWholeHeap(search, hp, end);
	  hp = hp->prev;
	} 
	hp->num_of_ref += 1;
	heap = hp->heap;
	if (hp->num_in_heap >= hsp_range_max) {
	    BlastPossibleDeleteWholeHeap(search, &hp, heap[0]);
	    if ((heap[0]->point_back->num_ref-=1)==0) 
	      hsp_deleted = TRUE;
	    heap = hp->heap;
	    heap[0] = hsp;
	    hsp->point_back->num_ref+=1;
	    BlastHeapify(hp, 0);
	    new_inserted = TRUE; 
	} else {
	  BlastInsertHeap(hp, hsp);
	  new_inserted = TRUE;
	}
	hsp->back_right = hp;
      } else hsp->back_right = NULL;
    }
    if (hsp_deleted) return 1;
    if (new_inserted) return 2;
    return 0;
}

void
BlastFreeHeap(BlastSearchBlkPtr search, BLASTResultHitlistPtr result_hitlist)
{ 
    BLASTResultHspPtr PNTR heap, hsp;
    BLASTResultHspPtr hsp_array;
    Int4 index, block_width, hsp_range_max;
    Int4 begin, end, i, hspcnt;
    BLASTHeapPtr hp;

    if (search->pbp->perform_culling == FALSE)  /* Culling is turned off. */
        return;

    block_width = search->pbp->block_width;  
    hsp_range_max = search->pbp->hsp_range_max;
    hspcnt = result_hitlist->hspcnt;
    hsp_array = result_hitlist->hsp_array;

    for (index = 0; index < hspcnt; index++) {
      hsp = &hsp_array[index];
      begin = hsp->query_offset;
      end = hsp->query_offset+hsp->query_length-1;
      if (hsp->back_left) {
	hp = hsp->back_left->next;
	if (BlastPossibleDeleteWholeHeap(search, &hp, hsp)) continue;
      } else {
	for (hp = search->result_struct->heap_ptr; hp; hp = hp->next) {
	  if (hp->cutvalue >= begin) break;
	}
      }
      for (; hp; hp = hp->next) {
	if (hp->cutvalue > end ) break;
	heap = hp->heap;
	for (i = 0; i < hp->num_in_heap; i++) {
	  if (heap[i] == hsp) {
	    BlastDeleteHeap(hp, i);
	    break;
	  }
	}
      }
      hp = hp->prev;
      BlastPossibleDeleteWholeHeap(search, &hp, hsp);
    }
}

/*converts a 1-level list of SeqAligns to a 2-level
list store in a ValNodePtr; the list lastSeqAligns
specifies where to break up seqAlignList */
ValNodePtr convertSeqAlignListToValNodeList(SeqAlignPtr seqAlignList, SeqAlignPtr * lastSeqAligns, Int4 numLastSeqAligns)
{

   ValNodePtr  returnValNodePtr, thisValNodePtr, nextValNodePtr;
   SeqAlignPtr thisSeqAlign;
   Int4 lastAlignIndex;

   returnValNodePtr = (ValNodePtr) MemNew (sizeof(ValNode));
   returnValNodePtr->data.ptrvalue = seqAlignList;
   returnValNodePtr->next = NULL;
   thisValNodePtr = returnValNodePtr;
   thisSeqAlign = seqAlignList;
   lastAlignIndex = 0;
   while ((NULL != thisSeqAlign) && (lastAlignIndex < (numLastSeqAligns -1))) {
   /*last in sublist but not last overall*/

     if ((thisSeqAlign == lastSeqAligns[lastAlignIndex]) &&
          (NULL != (thisSeqAlign->next))) {
       nextValNodePtr = (ValNodePtr) MemNew (sizeof(ValNode));
       nextValNodePtr->data.ptrvalue = thisSeqAlign->next;
       thisValNodePtr->next = nextValNodePtr;
       nextValNodePtr->next = NULL;
       thisValNodePtr = nextValNodePtr;
       thisSeqAlign = thisSeqAlign->next;
       lastSeqAligns[lastAlignIndex]->next = NULL;
       lastAlignIndex++;
     } 
     else      
       thisSeqAlign = thisSeqAlign->next;
   }
   return (returnValNodePtr);
}

/*converts a 2-level list of SeqAligns stored as a ValNodePtr */
 
SeqAlignPtr 
convertValNodeListToSeqAlignList(ValNodePtr seqAlignDoubleList, 
                                 SeqAlignPtr ** lastSeqAligns, 
                                 Int4 * numLastSeqAligns)
{
    ValNodePtr thisValNodePtr;
    SeqAlignPtr returnSeqAlign, thisSeqAlign;
    Int4 numValNodePtrs, indexValNodePtrs;
    
    thisValNodePtr = seqAlignDoubleList;
    numValNodePtrs = 0;

    while (NULL != thisValNodePtr) {
        numValNodePtrs++;
        thisValNodePtr = thisValNodePtr->next;
    }

    if (seqAlignDoubleList == NULL || 
        seqAlignDoubleList->data.ptrvalue == NULL) {
        *numLastSeqAligns = 0;
        *lastSeqAligns = NULL;
        return NULL;
    } else {
        indexValNodePtrs = 0;
        *numLastSeqAligns = numValNodePtrs;
        *lastSeqAligns = (SeqAlignPtr *) MemNew (numValNodePtrs * sizeof(SeqAlignPtr));
        returnSeqAlign = seqAlignDoubleList->data.ptrvalue;
        thisValNodePtr = seqAlignDoubleList;
        while (NULL != thisValNodePtr->next) {
            thisSeqAlign = thisValNodePtr->data.ptrvalue;
            while (thisSeqAlign->next != NULL) 
                thisSeqAlign = thisSeqAlign->next;
            (*lastSeqAligns)[indexValNodePtrs] = thisSeqAlign;
            indexValNodePtrs++;
            thisSeqAlign->next = thisValNodePtr->next->data.ptrvalue;
            thisValNodePtr = thisValNodePtr->next;
        }
        return returnSeqAlign;
    }
}

void LIBCALL
fillCandLambda(seedSearchItems * seedSearch, Char *matrixName, BLAST_OptionsBlkPtr options)
{
  if (0 == StringCmp("BLOSUM62", matrixName)) {
    seedSearch->paramC = 0.50;
    if ((11 == options->gap_open) && (1 == options->gap_extend)) {
      seedSearch->paramLambda = 0.270;
      seedSearch->paramK = 0.047;
      return;
    }
    if ((9 == options->gap_open) && (2 == options->gap_extend)) {
      seedSearch->paramLambda = 0.285;
      seedSearch->paramK = 0.075;
      return;
    }
    if ((8 == options->gap_open) && (2 == options->gap_extend)) {
      seedSearch->paramLambda = 0.265;
      seedSearch->paramK = 0.046;
      return;
    }
    if ((7 == options->gap_open) && (2 == options->gap_extend)) {
      seedSearch->paramLambda = 0.243;
      seedSearch->paramK = 0.032;
      return;
    }
    if ((12 == options->gap_open) && (1 == options->gap_extend)) {
      seedSearch->paramLambda = 0.281;
      seedSearch->paramK = 0.057;
      return;
    }
    if ((10 == options->gap_open) && (1 == options->gap_extend)) {
      seedSearch->paramLambda = 0.250;
      seedSearch->paramK = 0.033;
      return;
    }
    ErrPostEx(SEV_FATAL, 0, 0, "The combination %d for gap opening cost and %d for gap extension is not supported in PHI-BLAST with matrix %s\n", options->gap_open, options->gap_extend, matrixName);
  }
  else {
    if (0 == StringCmp("PAM30", matrixName)) { 
      seedSearch->paramC = 0.30;
      if ((9 == options->gap_open) && (1 == options->gap_extend)) {
	seedSearch->paramLambda = 0.295;
	seedSearch->paramK = 0.13;
	return;
      }
      if ((7 == options->gap_open) && (2 == options->gap_extend)) {
	seedSearch->paramLambda = 0.306;
	seedSearch->paramK = 0.15;
	return;
      }
      if ((6 == options->gap_open) && (2 == options->gap_extend)) {
	seedSearch->paramLambda = 0.292;
	seedSearch->paramK = 0.13;
	return;
      }
      if ((5 == options->gap_open) && (2 == options->gap_extend)) {
	seedSearch->paramLambda = 0.263;
	seedSearch->paramK = 0.077;
	return;
      }
      if ((10 == options->gap_open) && (1 == options->gap_extend)) {
	seedSearch->paramLambda = 0.309;
	seedSearch->paramK = 0.15;
	return;
      }
      if ((8 == options->gap_open) && (1 == options->gap_extend)) {
	seedSearch->paramLambda = 0.270;
	seedSearch->paramK = 0.070;
	return;
      }
      ErrPostEx(SEV_FATAL, 0, 0, "The combination %d for gap opening cost and %d for gap extension is not supported in PHI-BLAST with matrix %s\n", options->gap_open, options->gap_extend, matrixName);
    }
    else {
      if (0 == StringCmp("PAM70", matrixName)) { 
	seedSearch->paramC = 0.35;
	if ((10 == options->gap_open) && (1 == options->gap_extend)) {
	  seedSearch->paramLambda = 0.291;
	  seedSearch->paramK = 0.089;
	  return;
	}
	if ((8 == options->gap_open) && (2 == options->gap_extend)) {
	  seedSearch->paramLambda = 0.303;
	  seedSearch->paramK = 0.13;
	  return;
	}
	if ((7 == options->gap_open) && (2 == options->gap_extend)) {
	  seedSearch->paramLambda = 0.287;
	  seedSearch->paramK = 0.095;
	  return;
	}
	if ((6 == options->gap_open) && (2 == options->gap_extend)) {
	  seedSearch->paramLambda = 0.269;
	  seedSearch->paramK = 0.079;
	  return;
	}
	if ((11 == options->gap_open) && (1 == options->gap_extend)) {
	  seedSearch->paramLambda = 0.307;
	  seedSearch->paramK = 0.13;
	  return;
	}
	if ((9 == options->gap_open) && (1 == options->gap_extend)) {
	  seedSearch->paramLambda = 0.269;
	  seedSearch->paramK = 0.058;
	  return;
	}
	ErrPostEx(SEV_FATAL, 0, 0, "The combination %d for gap opening cost and %d for gap extension is not supported in PHI-BLAST with matrix %s\n", options->gap_open, options->gap_extend, matrixName);
      }
      else {
	if (0 == StringCmp("BLOSUM80", matrixName)) { 
	  seedSearch->paramC = 0.40;
	  if ((10 == options->gap_open) && (1 == options->gap_extend)) {
	    seedSearch->paramLambda = 0.300;
	    seedSearch->paramK = 0.072;
	    return;
	  }
	  if ((8 == options->gap_open) && (2 == options->gap_extend)) {
	    seedSearch->paramLambda = 0.308;
	    seedSearch->paramK = 0.089;
	    return;
	  }
	  if ((7 == options->gap_open) && (2 == options->gap_extend)) {
	    seedSearch->paramLambda = 0.295;
	    seedSearch->paramK = 0.077;
	    return;
	  }
	  if ((6 == options->gap_open) && (2 == options->gap_extend)) {
	    seedSearch->paramLambda = 0.271;
	    seedSearch->paramK = 0.051;
	    return;
	  }
	  if ((11 == options->gap_open) && (1 == options->gap_extend)) {
	    seedSearch->paramLambda = 0.314;
	    seedSearch->paramK = 0.096;
	    return;
	  }
	  if ((9 == options->gap_open) && (1 == options->gap_extend)) {
	    seedSearch->paramLambda = 0.277;
	    seedSearch->paramK = 0.046;
	    return;
	  }
	  ErrPostEx(SEV_FATAL, 0, 0, "The combination %d for gap opening cost and %d for gap extension is not supported in PHI-BLAST with matrix %s\n", options->gap_open, options->gap_extend, matrixName);
	}
	else {
	  if (0 == StringCmp("BLOSUM45", matrixName)) { 
	    seedSearch->paramC = 0.60;
	    if ((14 == options->gap_open) && (2 == options->gap_extend)) {
	      seedSearch->paramLambda = 0.199;
	      seedSearch->paramK = 0.040;
	      return;
	    }
	    if ((13 == options->gap_open) && (3 == options->gap_extend)) {
	      seedSearch->paramLambda = 0.209;
	      seedSearch->paramK = 0.057;
	      return;
	    }
	    if ((12 == options->gap_open) && (3 == options->gap_extend)) {
	      seedSearch->paramLambda = 0.203;
	      seedSearch->paramK = 0.049;
	      return;
	    }
	    if ((11 == options->gap_open) && (3 == options->gap_extend)) {
	      seedSearch->paramLambda = 0.193;
	      seedSearch->paramK = 0.037;
	      return;
	    }
	    if ((10 == options->gap_open) && (3 == options->gap_extend)) {
	      seedSearch->paramLambda = 0.182;
	      seedSearch->paramK = 0.029;
	      return;
	    }
	    if ((15 == options->gap_open) && (2 == options->gap_extend)) {
	      seedSearch->paramLambda = 0.206;
	      seedSearch->paramK = 0.049;
	      return;
	    }
	    if ((13 == options->gap_open) && (2 == options->gap_extend)) {
	      seedSearch->paramLambda = 0.190;
	      seedSearch->paramK = 0.032;
	      return;
	    }
	    if ((12 == options->gap_open) && (2 == options->gap_extend)) {
	      seedSearch->paramLambda = 0.177;
	      seedSearch->paramK = 0.023;
	      return;
	    }
	    if ((19 == options->gap_open) && (1 == options->gap_extend)) {
	      seedSearch->paramLambda = 0.209;
	      seedSearch->paramK = 0.049;
	      return;
	    }
	    if ((18 == options->gap_open) && (1 == options->gap_extend)) {
	      seedSearch->paramLambda = 0.202;
	      seedSearch->paramK = 0.041;
	      return;
	    }
	    if ((17 == options->gap_open) && (1 == options->gap_extend)) {
	      seedSearch->paramLambda = 0.195;
	      seedSearch->paramK = 0.034;
	      return;
	    }
	    if ((16 == options->gap_open) && (1 == options->gap_extend)) {
	      seedSearch->paramLambda = 0.183;
	      seedSearch->paramK = 0.024;
	      return;
	    }
	    ErrPostEx(SEV_FATAL, 0, 0, "The combination %d for gap opening cost and %d for gap extension is not supported in PHI-BLAST with matrix %s\n", options->gap_open, options->gap_extend, matrixName);

	  }
	  else {
	    ErrPostEx(SEV_FATAL, 0, 0, "Matrix %s not allowed in PHI-BLAST\n", matrixName);
          }
        }
      }
    }
  }
}

/* 
	Used by BlastParceInputString to get the 'index' of the
	option and to check it's validity.
*/
Int4 BlastGetLetterIndex(CharPtr letters, Char ch)
{
    Int4 index;

    for(index = 0; letters[index] != NULLB; index++) {
	if (letters[index] == ch) {
	    return index;
	}
    }
    return -1;
}

/*
	Parses a string of input options.
	For use in the Web page and filtering options. 
	Not for use in command-line programs - use GetArgs.
*/

Boolean BlastParceInputString(CharPtr string, 
	CharPtr letters, 
	CharPtr PNTR *values_in,
	CharPtr PNTR ErrorMessage)
{
    CharPtr chptr;
    Int4 i, index = 0, max_par_num;
    Char option[1024];
    CharPtr PNTR values;
    Char message[1024];

    if(string == NULL || letters == NULL || 
	    *letters == '\0' || values_in == NULL) {
	return FALSE;
    }

    max_par_num = StringLen(letters);

    values = (CharPtr PNTR)MemNew(max_par_num * sizeof(CharPtr));
    *values_in = values;

    chptr = string;

    while(1) {
	while(IS_WHITESP(*chptr)) /* Rolling spaces */
	    chptr++;

	if(*chptr == NULLB)   /* Check for NULLB */
	    break;

	if (*chptr != '-') {   /* Check for the option sign */
	    sprintf(message, "Invalid input string started from \"%s\"", 
		    chptr);
	    *ErrorMessage = StringSave(message);
	    return FALSE;
	} else {
	    chptr++;
	}

	/* checking index in options table */

	if((index = BlastGetLetterIndex(letters, *chptr)) < 0) {
	    sprintf(message, "Character \'%c\' is not a valid option", 
		    *chptr);
	    *ErrorMessage = StringSave(message);
	    return FALSE;
	}

	if(*chptr == NULLB)   /* Check for NULLB */
	    break;

	chptr++;

	while(IS_WHITESP(*chptr)) /* Rolling spaces */
	    chptr++;

	if(*chptr == NULLB)   /* Check for NULLB */
	    break;

	for(i=0; !IS_WHITESP(*chptr) && *chptr != NULLB; i++, chptr++) {
	    option[i] = *chptr;
	}

	option[i] = NULLB;

	MemFree(values[index]);
	values[index] = StringSave(option);
    }

    return TRUE;
}
