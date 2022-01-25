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
* File Name: blastcl2.c
*
* Author:  Roman L. Tatusov, Jonathan Epstein, Tom Madden
*
* Version Creation Date:   06/16/95
*
* $Revision: 6.0 $
*
* File Description: 
*       Simulates "traditional" BLAST output
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
* $Log: blastcl2.c,v $
* Revision 6.0  1997/08/25 18:34:02  madden
* Revision changed to 6.0
*
* Revision 5.3  1996/07/26 13:09:52  madden
* Added option to adjust level of dust filtering.
*
 * Revision 5.2  1996/06/12  18:06:18  madden
 * Added filtering (dust or seg) as an option ("-f").
 *
 * Revision 5.1  1996/06/04  12:19:37  madden
 * Removed define for MAX_SEQ_LEN.
 *
 * Revision 5.0  1996/05/28  14:09:11  ostell
 * Set to revision 5.0
 *
 * Revision 4.9  1996/05/22  12:52:10  madden
 * Removed unused variable "ptr".
 *
 * Revision 4.8  1996/05/06  13:15:34  madden
 * Added "PRE" messages around Motd for HTML format.
 *
 * Revision 4.7  1996/05/03  21:54:45  madden
 * Added options to produce HTML format and suppress motd.
 *
 * Revision 4.6  1996/04/04  22:43:07  madden
 * Added "queue" progress monitor.
 *
 * Revision 4.5  1996/02/24  19:02:54  madden
 * Removed second (superfluous) "main" function that did not use GetArgs.
 *
 * Revision 4.4  1995/10/24  15:57:56  madden
 * removed PrintTemplate stuff, cleaned up.
 *
 * Revision 4.3  1995/09/01  16:46:01  madden
 * Change to callback procedure so submissions can be cancelled.
 *
 * Revision 4.2  1995/08/15  13:06:12  madden
 * Removed unused SeqAnnotPtr (sap) from function blast.
 *
 * Revision 4.1  1995/08/03  21:21:10  madden
 * replaced fprintf to stderr with ErrPostEx.
 *
 * Revision 4.0  1995/07/26  13:55:34  ostell
 * force revision to 4.0
 *
 * Revision 1.6  1995/07/25  15:02:28  madden
 * Error messages returned from server printed out.
 *
 * Revision 1.5  1995/07/24  17:34:02  madden
 * Changed HitData to BLAST0Result
 *
 * Revision 1.4  1995/07/12  17:45:25  madden
 * Call BlastBioseq with new argument to perform masking.
 *
 * Revision 1.3  1995/06/23  22:14:00  madden
 * sixth argument in BlastBioseq call is now zero, rather than NULL.
 *
 * Revision 1.2  1995/06/22  17:08:16  madden
 * Added "output" argument to call to BlastBioseq.
 *
 * Revision 1.1  1995/06/16  11:26:33  epstein
 * Initial revision
 *
 * Revision 1.16  95/05/17  17:59:18  epstein
 * add RCS log revision history
 * 
*/
#define BLASTCLI_BUF_SIZE 255
#include <sequtil.h>
#include <prtutil.h>
#include <tofasta.h>
#include <netblap2.h>
#include <dust.h>


/* find the last nucleotide bioseq in the bioseqset */
static void FindNuc(SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
    BioseqPtr PNTR bp;
    BioseqPtr local_bsp;

    bp = (BioseqPtr PNTR) data;
    if (IS_Bioseq(sep))
    {
        local_bsp = (BioseqPtr) sep->data.ptrvalue;
        if (ISA_na(local_bsp->mol))
          *bp = local_bsp;
    }
}

/* find the last protein bioseq in the bioseqset */
static void FindProt(SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
    BioseqPtr PNTR bp;
    BioseqPtr local_bsp;

    bp = (BioseqPtr PNTR) data;
    if (IS_Bioseq(sep))
    {
        local_bsp = (BioseqPtr) sep->data.ptrvalue;
        if (ISA_aa(local_bsp->mol))
          *bp = local_bsp;
    }
}
/*
	Montior hook to print to stderr for UNIX clients.
*/

static int LIBCALLBACK UNIXMontiorHook(Nlm_MonitorPtr mon, MonCode code)

{
  switch (code) 
  {
#ifdef OS_UNIX
    case MonCode_Create :
	fprintf(stderr, "%s\n", (Nlm_CharPtr) mon->strTitle);
      break;
    case MonCode_StrValue :
	fprintf(stderr, "%s\n", (Nlm_CharPtr) mon->strValue);
      break;
#endif
    default :
      break;
  }
  return 0;

}

/********************************************************************
*
*	Add a string ("option" to the buffer).  "length" gives the
*	total length of the buffer.  If option can be added then
*	TRUE is returned, otherwise FALSE is returned.
********************************************************************/

static Boolean
AddOptionToBuffer(CharPtr buffer, CharPtr option, Int2 length)

{
	CharPtr ptr=buffer;
	Int2 count=0;

	if (option == NULL)
		return TRUE;

	while(*ptr != NULLB)
	{
		ptr++;
		count++;
	}

	if ((count+StringLen(option)) >= length)
		return FALSE;
	else
	{
		while (*option !=NULLB)
		{
			*ptr=*option;
			ptr++;
			option++;
		}
		*ptr=' ';
		ptr++;
		*ptr=NULLB;
	}
	return TRUE;
}

/*
	Callback to show progress in Queue, NOT in actual calculation
	of results.  
*/

static Boolean LIBCALLBACK
callback (BLAST0ResponsePtr brp, Boolean PNTR cancel)
{
	static MonitorPtr mon = NULL;
	Boolean retval;
	BLAST0QueuedPtr queue;
	Char buffer[40];

	if (brp->choice == BLAST0Response_queued)
	{
#ifdef OS_UNIX
		Nlm_SetMonitorHook(UNIXMontiorHook);
#endif
		if((queue=brp->data.ptrvalue) == NULL)
			return FALSE;

		if (mon == NULL)
		{
		   sprintf(buffer, "Waiting for %ld jobs to finish", queue->length);
		   mon=Nlm_MonitorStrNew (buffer, 40);
		}
		else
		{
		    retval = Nlm_MonitorStrValue (mon, "still waiting");
                    if (retval == FALSE)
                    { /* If cancelled, then shutdown monitor */
                        *cancel = TRUE;
                        mon = MonitorFree(mon);
                        return FALSE;
                    }
                }  
		return TRUE;
	}
	else
	{
		if (mon != NULL)
                        mon = MonitorFree(mon);
			
	}
	return FALSE;
}

static void
PrintMotd(BLAST0ResponsePtr blastResponse, FILE *fp, Boolean html_format)

{
	Char buffer[100];
	CharPtr string, ptr;

	if (blastResponse->choice != BLAST0Response_motd)
		return;

	string = blastResponse->data.ptrvalue;
	if (string == NULL)
		return;

	buffer[0] = NULLB;
	ptr = buffer;

	if (html_format)
	{
		fprintf(fp, "<PRE>\n");
	}

	while (*string != NULLB)
	{
		if (*string == '~')
		{
			*ptr = NULLB;
			fprintf(fp, "%s\n", buffer);
			buffer[0] = NULLB;
			ptr = buffer;
			string++;
			if (*string == NULLB)
				break;
		}
		else
		{
			*ptr=*string;
			ptr++;  string++;
		}
	}
	*ptr = NULLB;
	fprintf(fp, "%s\n", buffer);

	if (html_format)
	{
		fprintf(fp, "</PRE>\n");
	}

	fflush(fp);
}

/***********************************************************************
*	This function gets the Message-Of-The-Day, submits a request 
*	using BlastBioseq, formats the output, and reports any error
*	messages.
***********************************************************************/
static Int2 blast(BioseqPtr bsp, CharPtr blast_program, CharPtr blast_database, CharPtr cmd_options, Boolean HTML_option, Boolean print_motd, SeqLocPtr dust_slp, FILE *outfp)

{
	BLAST0PrefacePtr preface;
	BLAST0RequestPtr blreqp;
	BLAST0ResponsePtr brp, brp1, blastResponse;
	BLAST0ResultPtr result;
	BLAST0StatusPtr status;
	ValNodePtr vnp;

	if (bsp == NULL) {
		ErrPostEx(SEV_FATAL, 0, 0, "Couldn't read sequences");
		return (5);
	}
#ifdef MAX_SEQ_LEN
	if (bsp != NULL && bsp->length > MAX_SEQ_LEN) {
		ErrPostEx(SEV_FATAL, 0, 0, "Cannot process sequences > %d base pairs (this sequence is %d bp long)", MAX_SEQ_LEN, bsp->length);
		return (4);
	}
#endif
	
/* Get the MOTD from the server. */
	if (print_motd)
	{
		blreqp = ValNodeNew(NULL);
		blreqp->choice = BLAST0Request_motd;
		blreqp->data.ptrvalue = StringSave(blast_program);
		blastResponse = SubmitInfoRequest(blreqp);
		blreqp->data.ptrvalue = MemFree(blreqp->data.ptrvalue);
		blreqp = BLAST0RequestFree(blreqp);
		if (blastResponse == NULL)
		{
			ErrPostEx(SEV_FATAL, 0, 0, "BLAST connection to server failed");
			return (3);
		}
		PrintMotd(blastResponse, outfp, HTML_option);
		blastResponse = BLAST0ResponseFree(blastResponse);
	}
	brp = NULL;
	result = BlastBioseq(bsp, blast_program, blast_database, cmd_options, &brp, dust_slp, 0, callback);
	if (result != NULL) 
	{
#ifdef DEBUG
		AsnIoPtr aip;

		aip = AsnIoOpen("blastout.asn", "w");
		BLAST0ResultAsnWrite(result, aip, NULL);
		AsnIoClose(aip);
#endif /* DEBUG */
		if (outfp != NULL) 
		{
			if (HTML_option)
			{
				TraditionalBlastOutputHTML(result, brp, blast_program, outfp);
			}
			else
			{
				TraditionalBlastOutput(result, brp, blast_program, outfp);
			}
			result = BLAST0ResultFree(result);
			while (brp)
			{
				brp1 = brp;
				brp = brp->next;
				brp1->next = NULL;
				BLAST0ResponseFree(brp1);
			}
		} 
		else 
		{
			ErrPostEx(SEV_FATAL, 0, 0, "FileOpen failed");
			return (3);
		}
	}
	else if ((CheckIfBlastJobCancelled()) == TRUE)
	{
		return (7);
	}
	else
	{
		if (brp)
		{
			if (brp->choice == BLAST0Response_status)
			{
				blreqp = ValNodeNew(NULL);
				blreqp->choice = BLAST0Response_usage_info;
				blreqp->data.ptrvalue = StringSave(blast_program);
				blastResponse = SubmitInfoRequest(blreqp);
				blreqp = BLAST0RequestFree(blreqp);
				preface = blastResponse->data.ptrvalue;
				vnp = preface->prog_usage;
				print_usage(stderr, vnp);
				status = brp->data.ptrvalue;
				ErrPostEx(SEV_FATAL, 0, 0,
					"%s\nEXIT CODE %ld", status->reason, status->code);
			}
			else
			{
				brp1 = brp->next;
				while (brp1)
				{
				     if (brp1->choice == BLAST0Response_status)
				     {
					status = brp1->data.ptrvalue;
					ErrPostEx(SEV_FATAL, 0, 0,
						"%s\nEXIT CODE %ld", status->reason, status->code);
				     }
				     brp1 = brp1->next;
				}
			}
		}
		else
		{
			ErrPostEx(SEV_FATAL, 0, 0, "BLAST search failed");
			return (3);
		}
	}

	return 0;
}

#define NUMARG 9

static Args myargs [NUMARG] = {
  { "Program Name (blastn,blastp,blastx,tblastn,tblastx)", "", NULL, NULL, FALSE, 'p', ARG_STRING, 0.0, 0, NULL},
  { "Database", "nr", NULL, NULL, FALSE, 'd', ARG_STRING, 0.0, 0, NULL},
  { "Query File", "stdin", NULL, NULL, FALSE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
  { "Output File", "stdout", NULL, NULL, FALSE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
  { "HTML format", "F", NULL, NULL, FALSE, 'h', ARG_BOOLEAN, 0.0, 0, NULL},
  { "print MOTD", "T", NULL, NULL, FALSE, 'm', ARG_BOOLEAN, 0.0, 0, NULL},
  { "filter query sequence (nucleotides with dust, proteins with seg)", "F", NULL, NULL, FALSE, 'f', ARG_BOOLEAN, 0.0, 0, NULL},
  { "level of dust filtering", "-1", NULL, NULL, FALSE, 'l', ARG_INT, 0.0, 0, NULL},
  { "BLAST options (enclosed in double quotes, separated by spaces)", NULL, NULL, NULL, TRUE, 'b', ARG_STRING, 0.0, 0, NULL} 
};

/*********************************************************************
*	"main" function to call blast for the client.  
*
*	This function checks the command-line arguments, opens the
*	connection to the server, processes all the entries in
*	the FASTA file (obtained using FastaToSeqEntry), and
*	closes the connection.
*********************************************************************/
Int2 Main (void)

{
        CharPtr blast_program;
        CharPtr blast_database;
        CharPtr blast_inputfile;
        CharPtr blast_outputfile;
        CharPtr blast_params;
	Char buffer[BLASTCLI_BUF_SIZE+1];
        Boolean isprot = FALSE;
	BioseqPtr bsp;
	Int2 num_of_queries, retval;
	SeqEntryPtr sep;
	SeqLocPtr dust_slp=NULL;
	FILE *fp, *outfp;

	if (! GetArgs ("blastcli", NUMARG, myargs)) 
	{
		exit(1);
	}

	blast_program = myargs [0].strvalue;

	if (blast_program == NULL)
	{
		ErrPostEx(SEV_FATAL, 0, 0, 
			"blast: No name for BLAST program specified:\n");
		return (2);
	}

	if (StringCmp(blast_program, "blastp") == 0 ) {
		isprot = TRUE;
	} else if (StringCmp(blast_program, "blastx") == 0) {
		isprot = FALSE;
	} else if (StringCmp(blast_program, "blastn") == 0) {
		isprot = FALSE;
	} else if (StringCmp(blast_program, "tblastn") == 0) {
		isprot = TRUE;
	} else if (StringCmp(blast_program, "tblastx") == 0) {
		isprot = FALSE;
	} else {
		ErrPostEx(SEV_FATAL, 0, 0, "blast: Bad name for BLAST program: \"%s\"\n", blast_program);
		exit(1);
	}
	blast_database = myargs [1].strvalue;
	blast_inputfile = myargs [2].strvalue;
	blast_outputfile = myargs [3].strvalue;

	if ((fp = FileOpen(blast_inputfile, "r")) == NULL)
	{
		ErrPostEx(SEV_FATAL, 0, 0, "blast: Unable to open input file %s\n", blast_inputfile);
		return (1);
	}
	outfp = FileOpen (blast_outputfile, "w");

	if (! BlastInit("blastcl2", FALSE)) {
	      ErrPostEx(SEV_FATAL, 0, 0, "Unable to initialize BLAST service");
		return (1);
	}

	buffer[0] = NULLB;
	AddOptionToBuffer(buffer, myargs[8].strvalue, BLASTCLI_BUF_SIZE);
	if (myargs[6].intvalue)
	{
		if (StringCmp("blastn", blast_program) != 0)
			AddOptionToBuffer(buffer, "-filter seg", BLASTCLI_BUF_SIZE);
	}
	if (buffer[0] != NULLB)
	{
		blast_params = StringSave(buffer);
	}
	else
	{
		blast_params = NULL;
	}

	num_of_queries=0;
	retval=0;
	while ((sep = FastaToSeqEntry(fp, !isprot)) != NULL)
	{
		bsp = NULL;
		SeqEntryExplore(sep, &bsp, isprot? FindProt : FindNuc);
		if (bsp == NULL)
		{
	  	 	ErrPostEx(SEV_FATAL, 0, 0, "Unable to obtain bioseq\n");
	   		retval = 2;
			break;
		}

		if (myargs[6].intvalue)
		{
			if (StringCmp("blastn", blast_program) == 0)
				dust_slp = BioseqDust(bsp, 0, -1, myargs[7].intvalue, -1, -1, -1);
		}

		/* Put a space between the reports. */
		if (num_of_queries != 0)
			fprintf(outfp, "\n");
		retval = blast(bsp, blast_program, blast_database, blast_params, myargs [4].intvalue, myargs [5].intvalue, dust_slp, outfp);
		if (retval != 0)
			break;

		if (dust_slp)
			dust_slp = SeqLocSetFree(dust_slp);

		num_of_queries++;
	}
	FileClose(fp);
	FileClose (outfp);
	BlastFini();
	return retval;
}
