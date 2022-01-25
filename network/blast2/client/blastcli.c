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
* File Name: blastcli.c
*
* Author:  Tom Madden
*
* Version Creation Date:   08/22/95
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
* $Log: blastcli.c,v $
* Revision 6.0  1997/08/25 18:34:05  madden
* Revision changed to 6.0
*
* Revision 5.1  1997/04/11 15:18:38  madden
* Removed static function callbackWithMon, used Blast2callbackWithMon instead.
*
 * Revision 5.0  1996/05/28  14:09:11  ostell
 * Set to revision 5.0
 *
 * Revision 1.39  1996/04/25  16:51:56  kans
 * added file and edit menus to main window
 *
 * Revision 1.38  1996/04/23  21:26:41  madden
 * Added version number, hello message; cleaned out unused variable.
 *
 * Revision 1.37  1996/04/18  13:04:17  madden
 * Moved sequence input window into main window.
 *
 * Revision 1.36  1996/04/17  14:37:41  madden
 * Replaced many calls to "TraditionalBlast..." functions with call to
 * TraditionalBlastOutput.  Replaced Histogram check-box with Show GI
 * checkbox.
 *
 * Revision 1.35  1996/04/17  13:59:07  madden
 * Added "Dismiss" and "Clear" buttons to sequence input window.
 *
 * Revision 1.34  1996/04/16  20:56:14  madden
 * Sequence can now be entered from a window.
 *
 * Revision 1.33  1996/04/16  19:16:51  madden
 * Removed old functions used only with command-line interface.
 *
 * Revision 1.32  1996/04/16  19:10:04  madden
 * Bioseq now found in PerformBlastSearch rather than Calculate.
 *
 * Revision 1.31  1996/04/12  19:21:27  madden
 * Added "queue" monitor that shows when job is waiting to start.
 *
 * Revision 1.30  1996/03/25  17:34:09  shavirin
 * Changes due to new function TraditionalBlastOutputHTML(),
 * that produce HTML output from the BLAST server
 *
 * Revision 1.29  1996/03/19  21:24:27  madden
 * Removed opening of non-existent file under windows.
 *
 * Revision 1.28  1996/03/18  19:23:25  madden
 * Removed suggested output file for DOS and windows.
 *
 * Revision 1.27  1996/03/15  17:00:14  madden
 * Popup menu added for database list.
 *
 * Revision 1.26  1996/02/23  22:39:10  kans
 * removed extraneous fileclose, moved fileremove
 *
 * Revision 1.25  1995/12/15  15:03:21  epstein
 * add small-screen support
 *
 * Revision 1.24  1995/10/18  21:33:34  madden
 * Added check that blastResponse is not NULL before dereferencing it.
 *
 * Revision 1.23  1995/10/17  15:04:18  madden
 * Changes to accommodate CheckIfBlastJobCancelled.
 *
 * Revision 1.22  1995/10/16  21:43:28  madden
 * Changes to support command-line version of blastcli.
 *
 * Revision 1.21  1995/09/06  12:47:30  madden
 * Added message to help facility.
 *
 * Revision 1.20  1995/09/05  18:09:46  madden
 * missing BlastFini added in.
 *
 * Revision 1.19  1995/09/05  18:06:56  madden
 * Removed commented out callback function.
 *
 * Revision 1.18  1995/09/05  18:05:56  madden
 * BlastInits moved to avoid a problem.
 *
 * Revision 1.17  1995/09/05  17:20:38  madden
 * FileClose called before DisplayFile.
 *
 * Revision 1.16  1995/09/05  15:21:32  madden
 * Changed BLASTCLI_SMALL_BUF_SIZE to PATH_MAX on TmpNam call.
 *
 * Revision 1.15  1995/09/01  22:07:02  madden
 * converted some variable to "static" if they are only referenced in this file.
 *
 * Revision 1.14  1995/09/01  22:00:25  madden
 * Input file now read in with GetInputFileName.
 *
 * Revision 1.13  1995/09/01  16:43:11  madden
 * Changes to allow cancellation of jobs.
 *
 * Revision 1.12  1995/09/01  15:05:41  madden
 * enable "dismiss" button when another window is open.
 *
 * Revision 1.11  1995/09/01  13:36:40  madden
 * Default output file now named after input file.
 *
 * Revision 1.10  1995/08/31  22:02:30  madden
 * Added a second Main function that uses GetArgs.
 *
 * Revision 1.9  1995/08/31  20:53:23  madden
 * WIN_DUMB added as an ifndef.
 *
 * Revision 1.8  1995/08/31  15:02:27  madden
 * Moved num_of_descrp, num_of_align, num_of_expected to global variables.
 *
 * Revision 1.7  1995/08/31  14:21:33  madden
 * Added (number) extension to output file name.
 *
 * Revision 1.6  1995/08/31  12:04:41  madden
 * Changed from desktop GetArgs to customized Vibrant interface.
 *
 * Revision 1.5  1995/08/29  21:51:04  madden
 * removed unused variables that were lint complaints.
 *
 * Revision 1.4  1995/08/29  13:39:30  madden
 * Request no matrix or sequence data if alignments will not be shown.
 *
 * Revision 1.3  1995/08/28  21:24:46  madden
 * corrected index on myargs for call to AddOptionToBuffer.
 *
 * Revision 1.2  1995/08/28  17:34:18  madden
 * replaced custom functions with call to GetOutputFileName.
 *
 * Revision 1.1  1995/08/25  22:35:04  madden
 * Initial revision
 * */

#include <ncbi.h>
#include <sequtil.h>
#include <prtutil.h>
#include <tofasta.h>
#include <netblap2.h>
#include <dust.h>
#include <vibrant.h>
#include <document.h>
#include <blast2.h>

/* If WIN_DUMB is defined, don't use vibrant. */
#ifndef WIN_DUMB
#define BLASTCLI_WITH_VIBRANT
#endif

#define X_VIEW 750
#define Y_VIEW 500

#define BLASTCLI_SMALL_BUF_SIZE 25
#define BLASTCLI_BUF_SIZE 100

static DoC document;
static Nlm_WindoW results_window=NULL;
static Nlm_TexT input_text=NULL, sequence=NULL;
static ButtoN options_button, dismiss_button, submit_button, help_button;

CharPtr blast_output=NULL, new_blast_output=NULL, blast_program=NULL;

/* The following (global) variables are used as input for the BLAST
programs.  The large number of variables is to accomodate a command-line
GetArgs interface and a customized vibrant interface. */
static Boolean complex_filter, /* Mask low-complexity (dust on blastn, others seg). */
	get_gi;	/* Should gi's be displayed? */

/* The vibrant version stores the info in these buffers (e.g., expect_buffer), 
the information is then read into the Int4 below (e.g., num_of_expected). */
static Char
	expect_buffer[BLASTCLI_SMALL_BUF_SIZE], 
	descrp_buffer[BLASTCLI_SMALL_BUF_SIZE], 
	align_buffer[BLASTCLI_SMALL_BUF_SIZE], 
	expert_options[BLASTCLI_BUF_SIZE],
	score_buffer[BLASTCLI_SMALL_BUF_SIZE], 
	score2_buffer[BLASTCLI_SMALL_BUF_SIZE],
	blast_inputfile[PATH_MAX],
	blast_database[BLASTCLI_SMALL_BUF_SIZE];

static Int4 score, score2, num_of_descrp, num_of_align, num_of_expected;

/* contains list of databases, filled in GetListOfDatabases */
static ValNodePtr database_vnp=NULL;

static void PrintMotd PROTO((BLAST0ResponsePtr blastResponse, FILE *fp));

static void PrintBlastOptions PROTO((FILE *fp, BLAST0ResponsePtr brp));

static Int2 PerformBlastSearch PROTO((BioseqPtr bsp));

static Boolean AddOptionToBuffer PROTO((CharPtr buffer, CharPtr option, Int2 length));

static void PresentReportPanel PROTO((FILE *outfp));


/* 
	Clears the window users paste sequence into.
*/
static void ClearSeqInputWindow (ButtoN b)
{
	SetTitle(sequence, "");

}


static void HideResultsWindow (IteM i)
{
	if (results_window)
		Hide(results_window);
/* These were disabled in PresentReportPanel */
	Enable(help_button);
	Enable(options_button);
	Enable(submit_button);
/* Select this to keep the dismiss button from becoming the default. */
	Select(input_text);
}

static void ButtonQuitProc (ButtoN b)
{
	QuitProgram();
}

static void ItemQuitProc (IteM i)
{
	QuitProgram();
}

/*
	A function to check that the database and the blast program 
	are compatible.

	CharPtr program: name of program, e.g., blastn, blastp, etc.
	CharPtr database: database name, e.g., "nr", "est", etc.
	ValNodePtr vnp: ValNodePtr containing list of DB's (created in
	function GetListOfDatabases).	
*/

static Boolean
VerifyProgramDbCompatibility(CharPtr program, CharPtr database, ValNodePtr dblist)

{ 
	Uint2 type=0;

	if (program == NULL || *program == NULLB || database == NULL 
		|| *database == NULLB)
		return FALSE;

	for (dblist; dblist; dblist=dblist->next)
	{
		if (StringICmp(dblist->data.ptrvalue, database) == 0)
		{
			type = dblist->choice;
			break;
		}
	}

	if (type == 0)	/* Not found or not set. */
		return FALSE;

/* Check for compatability. */
	if ((type == BLAST0_Alphatype_nucleic_acid
	  && StringICmp("blastn", blast_program) != 0
            && StringICmp("tblastn", blast_program) != 0
              && StringICmp("tblastx", blast_program) != 0)
		|| (type == BLAST0_Alphatype_amino_acid
	  		&& StringICmp("blastp", blast_program) != 0
            			&& StringICmp("blastx", blast_program) != 0))
	{
		ErrPostEx(SEV_WARNING, 0, 0, "Incorrect Database type for %s search", blast_program);
		return FALSE;
	}

	return TRUE;
	
}

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

static void Calculate (ButtoN b)
{
	Int2 retval;
	BioseqPtr bsp;
	Boolean isprot;
	CharPtr input_sequence, last_char;
	SeqEntryPtr sep;
	FILE *infp;
	

	if (blast_program == NULL || *blast_program == NULLB)
	{
		ErrPostEx(SEV_WARNING, 0, 0, 
				"BLAST program name required for search");
		return;
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
		ErrPostEx(SEV_WARNING, 0, 0, "blast: Bad name for BLAST program: \"%s\"\n", blast_program);
		return;
	}

	if (blast_inputfile[0] != NULLB)
	{
		if ((infp = FileOpen(blast_inputfile, "r")) == NULL)
		{
			ErrPostEx(SEV_WARNING, 0, 0, "blast: Unable to open input file %s\n", blast_inputfile);
			return ;
		}

		sep = FastaToSeqEntry(infp, !isprot);
		FileClose(infp);
	}
	else if ((input_sequence=SaveStringFromText(sequence)) != NULL)
	{
		last_char = NULL;
		sep = FastaToSeqBuff(input_sequence, &last_char, !isprot);
		input_sequence = MemFree(input_sequence);
	}
	else
	{
		ErrPostEx(SEV_WARNING, 0, 0, "No input file or sequence specified\n");
		return;
	}
	

	bsp = NULL;
	SeqEntryExplore(sep, &bsp, isprot? FindProt : FindNuc);
	if (bsp == NULL)
	{
	   ErrPostEx(SEV_WARNING, 0, 0, "Unable to obtain sequence\n");
	   return;
	}


	if (blast_database == NULLB || *blast_database == NULLB)
	{
		ErrPostEx(SEV_WARNING, 0, 0, "Database name required");
		return;
	}

	if (VerifyProgramDbCompatibility(blast_program, blast_database, database_vnp) == FALSE)
		return;

	input_sequence = SaveStringFromText(sequence);

	WatchCursor();
/* These must be disabled to prevent conflicts. */
	Disable(help_button);
	Disable(options_button);
	Disable(submit_button);
/* Select this to keep the dismiss button from becoming the default. */
	Select(input_text);

	sscanf(descrp_buffer, "%ld", &num_of_descrp);
	sscanf(align_buffer, "%ld", &num_of_align);
	sscanf(expect_buffer, "%ld", &num_of_expected);
	score=0;
	if (score_buffer[0] != NULLB && StringCmp(score_buffer, "default") != 0)
		sscanf(score_buffer, "%ld", &score);
	score2=0;
	if (score2_buffer[0] != NULLB && 
		StringCmp(score2_buffer, "default") != 0)
		sscanf(score2_buffer, "%ld", &score2);
	retval = PerformBlastSearch(bsp);

	if (retval != 0)
	{	/* Panel never displayed */
		Enable(help_button);
		Enable(options_button);
		Enable(submit_button);
/* Select this to keep the dismiss button from becoming the default. */
		Select(input_text);
	}
		
	ArrowCursor();

	if (sep)
		sep = SeqEntryFree(sep);

		return;
}

static void GetOptionHelp (ButtoN b)
{
	FILE *outfp;
	
	if (blast_program == NULL || *blast_program == NULLB)
	{
		ErrPostEx(SEV_WARNING, 0, 0, 
				"BLAST program name required for option list");
		return;
	}

	blast_output = MemNew(PATH_MAX*sizeof(Char));
	blast_output = TmpNam(blast_output);

	outfp = FileOpen (blast_output, "w");
	PrintBlastOptions(outfp, NULL);
	FileRemove(blast_output);
	blast_output = MemFree(blast_output);

	return;
}


static void GetBlastHelp (ButtoN b)
{

	FILE *outfp;
	
	blast_output = MemNew(PATH_MAX*sizeof(Char));
	blast_output = TmpNam(blast_output);

	outfp = FileOpen (blast_output, "w");
	fprintf(outfp, "This help facility briefly describes the BLAST client interface.  To obtain more BLAST documentation\n");
	fprintf(outfp, "send a message consisting of just the word ``help'' (without quotes) to blast@ncbi.nlm.nih.gov\n\n");
	fprintf(outfp, "BLAST Help\n\n");
	fprintf(outfp, "Input FASTA file: file containing query sequence.\n\n");
	fprintf(outfp, "BLAST program: select one from the popup menu, choices are:\n");
	fprintf(outfp, "       blastn: compares nt query with nt database.\n");
	fprintf(outfp, "       blastp: compares aa query with aa database.\n");
	fprintf(outfp, "       blastx: compares translated nt query with aa database.\n");
	fprintf(outfp, "       tblastn: compares aa query with translated nt database.\n");
	fprintf(outfp, "       tblastx: compares translated nt query with translated nt database.\n");
	fprintf(outfp, "       (tblastx is only available for smaller databases such as dBEST or dbSTS)\n\n");
	fprintf(outfp, "Database: the default is non-redundant, or nr.\n\n");
	fprintf(outfp, "Number of one-line descriptions: number of high-scoring database sequences to be summarized.\n\n");
	fprintf(outfp, "Number of alignments: number of high-scoring database-query alignments to be shown.\n\n");
	fprintf(outfp, "Expectation number: number of matches expected by chance alone.\n\n");
	fprintf(outfp, "Minimum Score to report a database sequence:\n");
	fprintf(outfp, "	the default value is calculated from the Expectation number (above) and the database.\n\n");
	fprintf(outfp, "Minimum Score to report an HSP:\n");
	fprintf(outfp, "	a database sequence may have multiple High-Scoring Sequence Pairs (HSP).  This option\n");
	fprintf(outfp, "	sets a minimum score for each HSP, the default value depends on the program and database.\n\n");
	fprintf(outfp, "Show GI's: if checked, GI's are reported\n\n");
	fprintf(outfp, "Low-complexity filtering: for blastn dust is used, other programs use seg.\n\n");
	fprintf(outfp, "Other (expert) BLAST options may be entered.\n\n");
	fprintf(outfp, "A list of valid BLAST options may be requested from the server.\n\n");

/* These must be disabled to prevent conflicts. */
	Disable(help_button);
	Disable(options_button);
	Disable(submit_button);
/* Select this to keep the dismiss button from becoming the default. */
	Select(input_text);

	if (new_blast_output)
		new_blast_output = MemFree(new_blast_output);
	new_blast_output = StringSave("blhelp");

	PresentReportPanel(outfp);
	FileRemove(blast_output);

	blast_output = MemFree(blast_output);

	return;
}


static void SaveAsProc (IteM i)
{
	Char buffer[PATH_MAX];
	FILE *fp;

	if(GetOutputFileName(buffer, PATH_MAX, new_blast_output))
	{
		if (buffer[0] != NULLB)
		{
			fp = FileOpen(buffer, "w");
			SaveDocument(document, fp);
			FileClose(fp);
			HideResultsWindow(i);
		}
	}
}

static void
GetProgramType(PopuP p)

{
	Int2 program_number;
	
	if (blast_program != NULL)
		blast_program = MemFree(blast_program);

	program_number = GetValue (p);
	if (program_number == 1)
		blast_program = StringSave("blastn");
	else if (program_number == 2)
		blast_program = StringSave("blastp");
	else if (program_number == 3)
		blast_program = StringSave("blastx");
	else if (program_number == 4)
		blast_program = StringSave("tblastn");
	else if (program_number == 5)
		blast_program = StringSave("tblastx");
	return;
}

/*
	A Callback to get the database name selected by the user.
*/  
static void
GetDatabaseName(PopuP list)

{
	Int2 database_number, index;
	ValNodePtr vnp;

	database_number = GetValue (list);

	for (vnp=database_vnp, index=0; vnp; vnp=vnp->next)
	{
		index++;
		if (index == database_number)
			break;
	}
	VerifyProgramDbCompatibility(blast_program, vnp->data.ptrvalue, database_vnp);

	/* Set the string blast_database to the new db as the popup list 
	shows the new db as THE database. */
	StringNCpy(blast_database, vnp->data.ptrvalue, BLASTCLI_SMALL_BUF_SIZE);
	return;
}

/* 
	Gets the list of databases that apply to a give program.
	May be called multiple times if the program used changes.

	Two different databases with the same name (e.g., nr for 
	proteins and nucleic acids) are merged.  The two "types"
	(BLAST0_Alphatype_nucleic_acid and BLAST0_Alphatype_amino_acid)
	are added together to show that the database applies to 
	both molecules.
	
*/
static void
GetListOfDatabases (BLAST0ResponsePtr blresp)

{
	BLAST0DbDescPtr dbdesc, dbdesc1, last;
	Uint2 type;

	if (database_vnp)
		database_vnp = ValNodeFreeData(database_vnp); 

	/* Collect list of database names. */
	for (dbdesc=blresp->data.ptrvalue; dbdesc; dbdesc=dbdesc->next)
	{
		type = dbdesc->type;
		/* Check for redundant names (e.g., "nr"). */
		last = dbdesc;
		dbdesc1 = dbdesc->next;
		while (dbdesc1)
		{
			if (StringICmp(dbdesc->name, dbdesc1->name) == 0)
			{
				if (dbdesc1->type != type)
					type += dbdesc1->type;
				last->next = dbdesc1->next;
				dbdesc1->next = NULL;
				dbdesc1 = BLAST0DbDescFree(dbdesc1);
				dbdesc1 = last->next;
			}
			else
			{
				last = dbdesc1;
				dbdesc1 = dbdesc1->next;
			}
				
		}
		ValNodeCopyStr(&database_vnp, type, dbdesc->name);
	}
}

static void FileInActnProc (TexT text)

{
	GetTitle (text, blast_inputfile, PATH_MAX);
	return;
}

static void GetInputFile (ButtoN b)
{
	Char buffer[PATH_MAX];

	buffer[0] = NULLB;
	if(GetInputFileName(buffer, PATH_MAX, "*", "TEXT"))
	{
		
		if (buffer[0] != NULLB)
		{
			StringNCpy(blast_inputfile, buffer, PATH_MAX);
			SetTitle (input_text, blast_inputfile);
		}
	}
	return;
}

static void ExpertOptionActnProc (TexT text)

{
	GetTitle (text, expert_options, BLASTCLI_BUF_SIZE);
	return;
}

static void num_of_descrpActnProc (TexT text)

{
	GetTitle (text, descrp_buffer, BLASTCLI_SMALL_BUF_SIZE);
	return;
}

static void num_of_alignActnProc (TexT text)

{
	GetTitle (text, align_buffer, BLASTCLI_SMALL_BUF_SIZE);
	return;
}

static void num_of_expectActnProc (TexT text)

{
	GetTitle (text, expect_buffer, BLASTCLI_SMALL_BUF_SIZE);
	return;
}

static void ScoreActnProc (TexT text)

{
	GetTitle (text, score_buffer, BLASTCLI_SMALL_BUF_SIZE);
	return;
}

static void Score2ActnProc (TexT text)

{
	GetTitle (text, score2_buffer, BLASTCLI_SMALL_BUF_SIZE);
	return;
}

static void GetGi(ButtoN b)

{
	get_gi = GetStatus(b);
	return;
}

static void SaveFilter(ButtoN b)

{
	complex_filter = GetStatus(b);
	return;
}


/* The Boolean PNTR "cancel" is here only for PROTOTYPE "consistency" */
static Boolean LIBCALLBACK
callback (BLAST0ResponsePtr brp, Boolean PNTR cancel)
{
	return FALSE;
}

/*********************************************************************
*	"main" function to call blast for the client.
*	This function sets up the customized vibrant interface
*	for BLAST.  The actual call to BlastBioseq is with
*	PerformBlastSearch.
*********************************************************************/

Int2 Main (void)

{
        BLAST0RequestPtr blreqp;
        BLAST0ResponsePtr blresp;
	Char buffer[40];
	CharPtr blastcli_version_number="2.00";
	WindoW       wdialog;
	MenU     	 menu;
	GrouP        checkbox_group, decision_group, scroll_group,
			dialog_group=NULL, input_group=NULL, scroll_decision_group;
	PopuP        popd, db_popd;
	ButtoN       input_button, get_gi_button, complex_filter_button;
	ValNodePtr vnp;


	/* BlastInit needed for request of db-list. */
	if (! BlastInit("blastcli", FALSE)) 
	{
		ErrPostEx(SEV_FATAL, 0, 0, "Unable to initialize BLAST service");
		return (3);
	}

	sprintf(buffer, "BLAST Search %s", blastcli_version_number);
	wdialog = FixedWindow (-20, -60, -10, -10, buffer, NULL);

	menu = PulldownMenu (wdialog, "File");
	CommandItem (menu, "Quit", ItemQuitProc);
	menu = PulldownMenu (wdialog, "Edit");
	CommandItem(menu, "Cut", StdCutTextProc);
	CommandItem(menu, "Copy", StdCopyTextProc);
	CommandItem(menu, "Paste", StdPasteTextProc);
	CommandItem(menu, "Clear", StdDeleteTextProc);

	input_group = HiddenGroup (wdialog, 3, 0, NULL);

	StaticPrompt (input_group, "FASTA file input:", 0, dialogTextHeight, systemFont, '1');

	input_button = PushButton (input_group, "File", GetInputFile);

	input_text = DialogText (input_group, blast_inputfile, 20, FileInActnProc);

	scroll_group = HiddenGroup (wdialog, 1, 0, NULL);

	StaticPrompt (scroll_group, "Or cut and paste FASTA query into sequence window:", 0, 
		dialogTextHeight, systemFont, '1');
	blast_inputfile[0] = NULLB;

	sequence = ScrollText (scroll_group, 33, 6, programFont, TRUE, NULL);

	scroll_decision_group = HiddenGroup(scroll_group, 4, 0, NULL);

	PushButton (scroll_decision_group, "Clear sequence window", ClearSeqInputWindow);

	Break(wdialog);

	dialog_group = HiddenGroup (wdialog, 2, 100, NULL);

	StaticPrompt (dialog_group, "BLAST program", 0, popupMenuHeight, systemFont, '1');
	popd = PopupList (dialog_group, TRUE, GetProgramType);
	PopupItem (popd, "blastn");
	PopupItem (popd, "blastp");
	PopupItem (popd, "blastx");
	PopupItem (popd, "tblastn");
	PopupItem (popd, "tblastx");
/* Set the popup to blastn to start. */
	SetValue(popd, 1);
	if (blast_program != NULL)
		blast_program = MemFree(blast_program);
	blast_program = StringSave("blastn");


/* Use the same ValNode for two requests. */
/* Send Hello to server. */
	blreqp = ValNodeNew(NULL);
	blreqp->choice = BLAST0Request_hello;
	blreqp->data.ptrvalue = buffer;
	blresp = SubmitInfoRequest(blreqp);
	blresp = BLAST0ResponseFree(blresp);
	blreqp->data.ptrvalue = NULL;
/* Get the list of valid databases. */
	blreqp->choice = BLAST0Request_db_info;
	blresp = SubmitInfoRequest(blreqp);
	GetListOfDatabases(blresp);
	blreqp = BLAST0RequestFree(blreqp);
	blresp = BLAST0ResponseFree(blresp);


	StaticPrompt (dialog_group, "Database", 0, popupMenuHeight, systemFont, '1');
	db_popd = PopupList (dialog_group, TRUE, GetDatabaseName);
	for (vnp=database_vnp; vnp; vnp=vnp->next)
	{
		PopupItem (db_popd, (CharPtr) vnp->data.ptrvalue);
	}
/*Set the popup to the first database (nr?) */
	SetValue(db_popd, 1);
	StringNCpy(blast_database, database_vnp->data.ptrvalue, BLASTCLI_SMALL_BUF_SIZE);

	StaticPrompt (dialog_group, "Number of one-line descriptions to be reported (V)", 0, dialogTextHeight, systemFont, '1');
	StringNCpy(descrp_buffer, "250", BLASTCLI_SMALL_BUF_SIZE);
	DialogText (dialog_group, descrp_buffer, 5, num_of_descrpActnProc);

	StaticPrompt (dialog_group, "Number of alignments to be reported (B)", 0, dialogTextHeight, systemFont, '1');
	StringNCpy(align_buffer, "250", BLASTCLI_SMALL_BUF_SIZE);
	DialogText (dialog_group, align_buffer, 5, num_of_alignActnProc);

	StaticPrompt (dialog_group, "Expectation number (E)", 0, dialogTextHeight, systemFont, '1');
	StringNCpy(expect_buffer, "10", BLASTCLI_SMALL_BUF_SIZE);
	DialogText (dialog_group, expect_buffer, 5, num_of_expectActnProc);

	StaticPrompt (dialog_group, "Minimum Score to report a database sequence (S)", 0, dialogTextHeight, systemFont, '1');
	StringNCpy(score_buffer, "default", BLASTCLI_SMALL_BUF_SIZE);
	DialogText (dialog_group, score_buffer, 5, ScoreActnProc);

	StaticPrompt (dialog_group, "Minimum Score to report an HSP (S2)", 0, dialogTextHeight, systemFont, '1');
	StringNCpy(score2_buffer, "default", BLASTCLI_SMALL_BUF_SIZE);
	DialogText (dialog_group, score2_buffer, 5, Score2ActnProc);

	checkbox_group = HiddenGroup(dialog_group, 4, 0, NULL);
	get_gi_button = CheckBox(checkbox_group, "Show GI's", GetGi);
	SetStatus(get_gi_button, FALSE);
	complex_filter_button = CheckBox(checkbox_group, "Low-complexity filtering", SaveFilter);
	SetStatus(complex_filter_button, TRUE);

/* Set these in case their CheckBoxes are never touched. */
	complex_filter=TRUE;
	get_gi=FALSE;

	Break(dialog_group);
	StaticPrompt (dialog_group, "Other BLAST Options", 0, dialogTextHeight, systemFont, '1');
	DialogText (dialog_group, expert_options, 5, ExpertOptionActnProc);
	options_button = PushButton (dialog_group, "Option List", GetOptionHelp);
	Break(dialog_group);


	Break(dialog_group);
	decision_group = HiddenGroup(dialog_group, 4, 0, NULL);
	submit_button = DefaultButton (decision_group, "Submit", Calculate );
	dismiss_button = PushButton (decision_group, "Dismiss", ButtonQuitProc);
	help_button = PushButton (decision_group, "Help", GetBlastHelp);


	Show (wdialog); 
	ProcessEvents();
	BlastFini();

	return 0;
}


#define BLASTCLI_FILE_NAME_SIZE 8 /* choose this size for DOS */
static Int2
PerformBlastSearch (BioseqPtr bsp)

{
	BLAST0RequestPtr blreqp;
	BLAST0StatusPtr status;
	CharPtr filename, ptr, ptr1;
	Char buf_tmp[BLASTCLI_SMALL_BUF_SIZE+1], buffer[BLASTCLI_BUF_SIZE];
	BLAST0ResponsePtr brp, brp1, blastResponse;
	BLAST0ResultPtr result;
	Int2 index;
	static Int2 file_index=0;
	SeqLocPtr dust_slp=NULL;
	FILE *outfp;
	Uint4 output=0;

/* Get the MOTD from the server. */
	blreqp = ValNodeNew(NULL);
	blreqp->choice = BLAST0Request_motd;
	blreqp->data.ptrvalue = StringSave(blast_program);
	blastResponse = SubmitInfoRequest(blreqp);
	if (blastResponse == NULL)
	{
		ErrPostEx(SEV_FATAL, 0, 0, "BLAST connection to server failed");
		return (3);
	}
	blast_output = MemNew(PATH_MAX*sizeof(Char));
	blast_output = TmpNam(blast_output);
/* Extract the filename from the entire path. */ 
	filename = FileNameFind(blast_inputfile);
	StringNCpy(buf_tmp, filename, BLASTCLI_FILE_NAME_SIZE);
	index = StringLen(filename);
/* Both of these must be done to assure no garbage in string. */ 
	buf_tmp[index] = NULLB;
	buf_tmp[BLASTCLI_FILE_NAME_SIZE] = NULLB;
	ptr = buffer;
	for (index=0; 
		index<BLASTCLI_FILE_NAME_SIZE, buf_tmp[index] != NULLB; 
			index++)
	{
		*ptr = buf_tmp[index];
		ptr++;
	}
		
	*ptr = '.';
	ptr++;
	file_index++; /* incremented every time this function is called. */
	ptr1 = Nlm_Ltostr((long) file_index, 1);
	while (*ptr1 != NULLB)
	{
		*ptr = *ptr1;
		ptr++; ptr1++;
	}
	
	*ptr = NULLB;
	if (StringCmp(buffer, filename) == 0)
	{
		*ptr = 'X';
		ptr++;
	}
	*ptr = NULLB;

/* IF WIN_MSWIN is set, then DO NOT suggest an output file, as it will probably
collide with the input file because of the 8+3 naming convention of DOS. */
#ifdef WIN_MSWIN
	new_blast_output = NULL;
#else
/* There is no input filename if the sequence is entered through a window. */
	if (filename != NULL && filename[0] != NULLB)
	{
		new_blast_output = StringSave(buffer);
	}
	else
	{
		new_blast_output = NULL;
		
	}
#endif

	outfp = FileOpen (blast_output, "w");
	PrintMotd(blastResponse, outfp);

/* Compose the "options" line. */

	buffer[0]=NULLB;
	/* "V" and "B" options that limit the number of deflines and alignment*/
	sprintf(buf_tmp, "V=%ld B=%ld", (long) num_of_descrp, (long) num_of_align);
	AddOptionToBuffer(buffer, buf_tmp, BLASTCLI_BUF_SIZE);
	/* E (expected number of hits). */
	if (num_of_expected)
	{
		sprintf(buf_tmp, "E=%ld", (long) num_of_expected);
		AddOptionToBuffer(buffer, buf_tmp, BLASTCLI_BUF_SIZE);
	}

	/* S (minimum score to report a db sequence) */
	if (score)
	{
		sprintf(buf_tmp, "S=%ld", (long) score);
		AddOptionToBuffer(buffer, buf_tmp, BLASTCLI_BUF_SIZE);
	}

	/* S2 (minimum score to report an HSP) */
	if (score2)
	{
		sprintf(buf_tmp, "S2=%ld", (long) score2);
		AddOptionToBuffer(buffer, buf_tmp, BLASTCLI_BUF_SIZE);
	}
		
	/* SHould GI's be shown? */
	if (get_gi)
		AddOptionToBuffer(buffer, "-gi", BLASTCLI_BUF_SIZE);
	
	/* Low comlexity filtering? */
	dust_slp=NULL;
	if (complex_filter)
	{
		if (StringCmp("blastn", blast_program) != 0)
			AddOptionToBuffer(buffer, "-filter seg", BLASTCLI_BUF_SIZE);
		else
			dust_slp = BioseqDust(bsp, 0, -1, -1, -1, -1, -1);
	}

	/* Other (expert) BLAST options. */
	if (expert_options != NULL && *expert_options != NULLB)
		AddOptionToBuffer(buffer, expert_options, BLASTCLI_BUF_SIZE);

/* If no alignments are requested, then the matrix and sequence data is not
needed. */
	if (num_of_align == 0)
	{
		output = BLAST_SERVER_OMIT_MATRIX;
		output += BLAST_SERVER_OMIT_QUERY_SEQ_IN_SEG;
		output += BLAST_SERVER_OMIT_DB_SEQ_IN_SEG;
		ptr=expert_options;
		/* Only a very crude test is done for "V" as an option.*/
		if (ptr != NULL)
		{
			while (*ptr != NULLB)
			{
				if (*ptr == 'V' || *ptr == 'v')
				{
					if (*(ptr+1) == '=')
					{
						output=0;
						break;
					}
				}
				ptr++;
			}
		}
	}

/* Submit the request. */
	brp = NULL;

	result = BlastBioseq(bsp, blast_program, blast_database, buffer, &brp, dust_slp, output, Blast2callbackWithMon);

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
			TraditionalBlastOutput(result, brp, blast_program, outfp);

			BLAST0ResultFree(result);

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
			ErrPostEx(SEV_WARNING, 0, 0, "FileOpen failed");
			return (3);
		}
	}
	else if ((CheckIfBlastJobCancelled()) == TRUE)
	{
		return (7);	/*non-zero return value enables buttons */
	}
	else
	{
		if (brp)
		{
			if (brp->choice == BLAST0Response_status)
			{
				PrintBlastOptions(outfp, brp);
				return (5);
			}
			else
			{
				brp1 = brp->next;
				while (brp1)
				{
				     if (brp1->choice == BLAST0Response_status)
				     {
					status = brp1->data.ptrvalue;
					ErrPostEx(SEV_WARNING, 0, 0,
						"%s\nEXIT CODE %ld", status->reason, status->code);
				     }
				     brp1 = brp1->next;
				}
				return (6);
			}
		}
		
		else
		{
			ErrPostEx(SEV_FATAL, 0, 0, "BLAST search failed");
			return (3);
		}
	}
	if (dust_slp)
		dust_slp = SeqLocSetFree(dust_slp);


	PresentReportPanel(outfp);
/* Don't remove the output file if command-line version is compiled. */
#ifdef BLASTCLI_WITH_VIBRANT
	FileRemove(blast_output);
#endif

	return 0;
}

static void
PrintBlastOptions(FILE *fp, BLAST0ResponsePtr brp)
{
	BLAST0RequestPtr blreqp;
	BLAST0ResponsePtr blastResponse;
	BLAST0PrefacePtr preface;
	BLAST0StatusPtr status=NULL;
	ValNodePtr vnp;

	blreqp = ValNodeNew(NULL);
	blreqp->choice = BLAST0Request_usage_info;
	blreqp->data.ptrvalue = StringSave(blast_program);
	blastResponse = SubmitInfoRequest(blreqp);
	preface = blastResponse->data.ptrvalue;
	vnp = preface->prog_usage;

	if (brp != NULL)
	{
		status = brp->data.ptrvalue;
#ifdef BLASTCLI_WITH_VIBRANT
		if (fp != NULL) 
			fprintf(fp, "%s\nEXIT CODE %ld\n\n", status->reason, status->code);
#endif
	}

	if (fp != NULL) 
		print_usage(fp, vnp);

#ifdef BLASTCLI_WITH_VIBRANT
/* These must be disabled to prevent conflicts. */
	Disable(help_button);
	Disable(options_button);
	Disable(submit_button);
/* Select this to keep the dismiss button from becoming the default. */
	Select(input_text);

	if (new_blast_output)
		new_blast_output = MemFree(new_blast_output);
	new_blast_output = StringSave("blopt");

	PresentReportPanel(fp);
#else
	if (status != NULL)
		ErrPostEx(SEV_WARNING, 0, 0, "%s\nEXIT CODE %ld\n Valid options have been written to the file %s", status->reason, status->code, blast_output);
#endif
	blreqp = ValNodeFreeData(blreqp); 
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

static void
PrintMotd(BLAST0ResponsePtr blastResponse, FILE *fp)

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

	fflush(fp);
}

static void
PresentReportPanel(FILE *outfp)

{
	MenU      menu;
	Int2      doc_x;
	Int2      doc_y;

#ifdef BLASTCLI_WITH_VIBRANT
	WatchCursor();
	if (results_window == NULL)
	{
		results_window = FixedWindow(-50, -33, -10, -10, "blastcli", NULL);
		menu = PulldownMenu (results_window, "File");
		CommandItem(menu, "Save As", SaveAsProc);
		CommandItem(menu, "Quit", HideResultsWindow);
		doc_x = MIN(X_VIEW, (screenRect.right - 50));
		doc_y = MIN(Y_VIEW, (screenRect.bottom - 50));
		document = DocumentPanel (results_window, doc_x, doc_y);
	}
	FileClose (outfp);
	DisplayFile (document, blast_output, programFont);
	ArrowCursor();
	Show (results_window);
/* blast_output was created by TmpNam() */
#else
	FileClose (outfp);
#endif
}

