

/**********************************************************
*
*	provide a Graphic User Interface for PowerBlast
*
***********************************************************/
#include <vibrant.h>
#include <document.h>
#include <pblproc.h>
#include <lsqfetch.h>
#include <accentr.h>
#include <sqnutils.h>


/* for the input file name or sequences from cut and paste */
static Nlm_TexT input_text=NULL, sequence=NULL;
static Nlm_TexT repeat_text = NULL;
static ButtoN	search_button;


/* 
	Clears the window users paste sequence into.
*/
static void ClearSeqInputWindow (ButtoN b)
{
	WindoW w;
	PBlastOptionPtr pbop;

	SetTitle(sequence, "");

	w = ParentWindow(b);
	pbop = (PBlastOptionPtr)GetObjectExtra(w);
	if(pbop != NULL)
	{
		if(pbop->seq_data != NULL)
			pbop->seq_data = MemFree(pbop->seq_data);
		if(pbop->file_name[0] == '\0')
			Disable(search_button);
	}
}

static void FreePBlastOptionData(GraphiC g, VoidPtr data)
{
	PBlastOptionPtr pbop;

	pbop = (PBlastOptionPtr)data;
	if(pbop->errfp != NULL)
		FileClose(pbop->errfp);
	FreePBlastOption(pbop);
}


static void ButtonQuitProc (ButtoN b)
{
	WindoW w;

	w = ParentWindow(b);
	Remove(w);

	QuitProgram();
}

static void ItemQuitProc (IteM i)
{

	WindoW w;

	w = ParentWindow(i);
	Remove(w);
	QuitProgram();
}

static void CloseMainWinProc (WindoW w)
{
	Remove(w);
	BioseqFetchDisable();
	EntrezBioseqFetchDisable();
	BlastFini();
	QuitProgram();
}


static void FileInActnProc (TexT text, Boolean is_repeat)
{
	WindoW w;
	PBlastOptionPtr pbop;

	w = ParentWindow(text);
	pbop = (PBlastOptionPtr)GetObjectExtra(w);
	if(pbop == NULL)
		return;
	
	if(is_repeat)
	{
		pbop->repeat_library[0] = '\0';
		GetTitle(text, pbop->repeat_library, PATH_MAX);
	}
	else
	{
		pbop->file_name[0] = '\0';
		GetTitle (text, pbop->file_name, PATH_MAX);
		if(pbop->file_name[0] != '\0')
			Enable(search_button);
		else if(pbop->seq_data == NULL)
			Disable(search_button);
	}
}

static void FileInForInput(TexT text)
{
	FileInActnProc(text, FALSE);
}

static void FileInForRepeat(TexT text)
{
	FileInActnProc(text, TRUE);
}

static void GetInputFile (ButtoN b, Boolean is_repeat)
{
	Char buffer[PATH_MAX];
	WindoW w;
	PBlastOptionPtr pbop;

	w = ParentWindow(b);
	pbop = (PBlastOptionPtr)GetObjectExtra(w);
	if(pbop == NULL)
		return;

	buffer[0] = NULLB;
	if(GetInputFileName(buffer, PATH_MAX, NULL, NULL))
	{
		
		if (buffer[0] != NULLB)
		{
			if(is_repeat)
			{
				StringNCpy(pbop->repeat_library, buffer, PATH_MAX);
				SetTitle(repeat_text, pbop->repeat_library);
			}
			else
			{
				StringNCpy(pbop->file_name, buffer, PATH_MAX);
				SetTitle (input_text, pbop->file_name);
				if(pbop->file_name[0] != '\0')
					Enable(search_button);
			}
		}
	}
}

static void GetInputForBlast (ButtoN b)
{
	GetInputFile(b, FALSE);
}

static void GetInputForRepeat (ButtoN b)
{
	GetInputFile(b, TRUE);
}

static void CancelSettingProc(ButtoN b)
{
	WindoW w;
	PBlastOptionPtr pbop;

	w = ParentWindow(b);
	pbop = (PBlastOptionPtr)GetObjectExtra(w);
	if(pbop->temp_bpp != NULL)
	{
		free_blast_param(pbop->temp_bpp);
		pbop->temp_bpp = NULL;
	}
	Remove(w);
}

static void AcceptSettingProc(ButtoN b)
{
	WindoW w;
	PBlastOptionPtr pbop;

	w = ParentWindow(b);
	pbop = (PBlastOptionPtr)GetObjectExtra(w);
	if(pbop->bpp != NULL)
	{
		free_blast_param(pbop->bpp);
		pbop->bpp = pbop->temp_bpp;
		pbop->temp_bpp = NULL;
	}
	Remove(w);
}

static BlastParamPtr dup_blast_param (BlastParamPtr bpp)
{
	BlastParamPtr new_bpp;

	new_bpp = MemNew(sizeof(BlastParameter));
	new_bpp->blast_program = bpp->blast_program;
	new_bpp->dna_db= bpp->dna_db;
	new_bpp->prot_db = bpp->prot_db;
	new_bpp->max_blast_n = bpp->max_blast_n;
	new_bpp->max_blast_p = bpp->max_blast_p;
	new_bpp->max_blast_x = bpp->max_blast_x;
	new_bpp->max_blast_t = bpp->max_blast_t;

	if(bpp->n_param != NULL)
		new_bpp->n_param = StringSave(bpp->n_param);
	if(bpp->p_param != NULL)
		new_bpp->p_param = StringSave(bpp->p_param);
	if(bpp->t_param != NULL)
		new_bpp->t_param = StringSave(bpp->t_param);
	if(bpp->x_param != NULL)
		new_bpp->x_param = StringSave(bpp->x_param);
	if(bpp->other_dna!= NULL)
		new_bpp->other_dna= StringSave(bpp->other_dna);
	if(bpp->other_protein!= NULL)
		new_bpp->other_protein= StringSave(bpp->other_protein);
	return new_bpp;
}


/************************************************************
*
*	functions related to the set up of the database
*
*************************************************************/
static void SetBlastDB (ButtoN b, Uint4 value, Boolean is_protein)
{
	WindoW w;
	PBlastOptionPtr pbop;

	w = ParentWindow(b);
	pbop = (PBlastOptionPtr)GetObjectExtra(w);
	if(pbop != NULL && pbop->temp_bpp != NULL)
	{
		if(is_protein)
		{
			if(pbop->temp_bpp->prot_db & value)
				pbop->temp_bpp->prot_db -= value;
			else
				pbop->temp_bpp->prot_db |= value;
		}
		else
		{
			if(pbop->temp_bpp->dna_db & value)
				pbop->temp_bpp->dna_db -= value;
			else
				pbop->temp_bpp->dna_db |= value;
		}
	}
}

static void ResetNRForDNA(ButtoN b)
{
	SetBlastDB(b, BLAST_NR, FALSE);
}
				
static void ResetESTForDNA(ButtoN b)
{
	SetBlastDB(b, BLAST_EST, FALSE);
}
				
static void ResetSTSForDNA(ButtoN b)
{
	SetBlastDB(b, BLAST_STS, FALSE);
}
				
static void ResetMonthForDNA(ButtoN b)
{
	SetBlastDB(b, BLAST_MONTH, FALSE);
}
				
static void ResetVectorForDNA(ButtoN b)
{
	SetBlastDB(b, BLAST_VECTOR, FALSE);
}
				
static void ResetMitoForDNA(ButtoN b)
{
	SetBlastDB(b, BLAST_MITO, FALSE);
}

static void ResetEPDForDNA(ButtoN b)
{
	SetBlastDB(b, BLAST_EPD, FALSE);
}

static void ResetAluForDNA(ButtoN b)
{
	SetBlastDB(b, BLAST_ALU, FALSE);
}

static void ResetYeastForDNA(ButtoN b)
{
	SetBlastDB(b, BLAST_YEAST, FALSE);
}

static void ResetGSSForDNA(ButtoN b)
{
	SetBlastDB(b, BLAST_GSS, FALSE);
}

static void ResetPDBForDNA(ButtoN b)
{
	SetBlastDB(b, BLAST_PDB, FALSE);
}

static void ResetHTGSForDNA(ButtoN b)
{
	SetBlastDB(b, BLAST_HTGS, FALSE);
}


static void ResetNRForProt(ButtoN b)
{
	SetBlastDB(b, BLAST_NR, TRUE);
}
				
static void ResetPDBForProt(ButtoN b)
{
	SetBlastDB(b, BLAST_PDB, TRUE);
}
				
static void ResetSwissProtForProt(ButtoN b)
{
	SetBlastDB(b, BLAST_SWISSPROT, TRUE);
}
				
static void ResetMonthForProt(ButtoN b)
{
	SetBlastDB(b, BLAST_MONTH, TRUE);
}
				
static void ResetKabatForDNA(ButtoN b)
{
	SetBlastDB(b, BLAST_KABAT, TRUE);
}

static void ResetKabatForProt(ButtoN b)
{
	SetBlastDB(b, BLAST_KABAT, TRUE);
}

static void ResetYeastForProt(ButtoN b)
{
	SetBlastDB(b, BLAST_YEAST, TRUE);
}

static void ResetAluForProt(ButtoN b)
{
	SetBlastDB(b, BLAST_ALU, TRUE);
}

				

static void ResetOtherDatabase(TexT t, Boolean is_protein)
{
	WindoW w;
	PBlastOptionPtr pbop;
	Char buf[101];

	w = ParentWindow(t);
	pbop = (PBlastOptionPtr)GetObjectExtra(w);
	if(pbop != NULL && pbop->temp_bpp != NULL)
	{
		GetTitle(t, buf, 100);
		if(is_protein)
		{
			MemFree(pbop->temp_bpp->other_protein);
			pbop->temp_bpp->other_protein = StringSave(buf);
		}
		else
		{
			MemFree(pbop->temp_bpp->other_dna);
			pbop->temp_bpp->other_dna= StringSave(buf);
		}
	}
}

static void ResetBlastParameter(TexT t, Uint1 program)
{
	WindoW w;
	PBlastOptionPtr pbop;
	Char buf[201];

	w = ParentWindow(t);
	pbop = (PBlastOptionPtr)GetObjectExtra(w);
	if(pbop != NULL && pbop->temp_bpp != NULL)
	{
		GetTitle(t, buf, 200);
		switch(program)
		{
			case BLASTN_PROGRAM:
				MemFree(pbop->temp_bpp->n_param);
				pbop->temp_bpp->n_param = StringSave(buf);
				break;

			case BLASTP_PROGRAM:
				MemFree(pbop->temp_bpp->p_param);
				pbop->temp_bpp->p_param = StringSave(buf);
				break;

			case TBLASTN_PROGRAM:
				MemFree(pbop->temp_bpp->t_param);
				pbop->temp_bpp->t_param = StringSave(buf);
				break;

			case BLASTX_PROGRAM:
				MemFree(pbop->temp_bpp->x_param);
				pbop->temp_bpp->x_param = StringSave(buf);
				break;

			default:
				break;
		}
	}
}

static void ResetBlastNParameter(TexT t)
{
	ResetBlastParameter(t, BLASTN_PROGRAM);
}

static void ResetBlastPParameter(TexT t)
{
	ResetBlastParameter(t, BLASTP_PROGRAM);
}

static void ResetTBlastNParameter(TexT t)
{
	ResetBlastParameter(t, TBLASTN_PROGRAM);
}

static void ResetBlastXParameter(TexT t)
{
	ResetBlastParameter(t, BLASTX_PROGRAM);
}

static void ResetOtherNucDB(TexT t)
{
	ResetOtherDatabase(t, FALSE);
}

static void ResetOtherProtDB(TexT t)
{
	ResetOtherDatabase(t, TRUE);
}

static void HelpSettingProc(ButtoN b)
{
	Char buf[101];

	sprintf(buf, "Send a message consisting of just the word \n\"help\" (without quotes) to blast@ncbi.nlm.nih.gov");
	Message(MSG_OK, "%s", buf);
}

static void ResetBlastProgram (ButtoN b, Uint1 choice)
{
	WindoW w;
	PBlastOptionPtr pbop;
	BlastParamPtr bpp;
	Uint1 value;

	w = ParentWindow(b);
	pbop = (PBlastOptionPtr)GetObjectExtra(w);
	if(pbop == NULL || pbop->temp_bpp == NULL)
		return;
	switch(choice)
	{
		case 1:
			value = BLASTN_PROGRAM;
			break;
		case 2:
			value = TBLASTN_PROGRAM;
			break;
		case 3:
			value = BLASTP_PROGRAM;
			break;
		case 4:
			value = BLASTX_PROGRAM;
			break;

		default:
			return;
	}
	bpp = pbop->temp_bpp;
	if(bpp->blast_program & value)
		bpp->blast_program -= value;
	else
		bpp->blast_program |= value;
}

static void ResetBLASTN(ButtoN b)
{
	ResetBlastProgram (b, 1);
}

static void ResetTBLASTN(ButtoN b)
{
	ResetBlastProgram (b, 2);
}

static void ResetBLASTP(ButtoN b)
{
	ResetBlastProgram (b, 3);
}

static void ResetBLASTX(ButtoN b)
{
	ResetBlastProgram (b, 4);
}

		
static void CloseParamWinProc (WindoW w)
{
	PBlastOptionPtr pbop;

	pbop = (PBlastOptionPtr)GetObjectExtra(w);
	if(pbop->temp_bpp != NULL)
	{
		free_blast_param(pbop->temp_bpp);
		pbop->temp_bpp = NULL;
	}
	Remove(w);

}

static void SetBlastParamProc(ButtoN b)
{
	WindoW w;
	PBlastOptionPtr pbop;
	BlastParamPtr bpp;
	WindoW pw;
	GrouP g, gg, tg;
	ButtoN bb;

	w = ParentWindow(b);
	pbop = (PBlastOptionPtr)GetObjectExtra(w);
	if(pbop == NULL || pbop->bpp == NULL)
		return;
	bpp = pbop->bpp;
	if(pbop->temp_bpp != NULL)
		free_blast_param(pbop->temp_bpp);
	pbop->temp_bpp = dup_blast_param(pbop->bpp);

	pw = MovableModalWindow(-50, -33, -10, -10, 
		"Parameter and DataBase for Blast Search", CloseParamWinProc);


	/*set up to search the nucleotide database */
	tg = NormalGroup(pw, 0, 2, "Search Nucleotide DataBase", systemFont, NULL);
	g = HiddenGroup(tg, 4, 0, NULL);
	SetGroupSpacing(g, 5, 0);
	bb = CheckBox(g, "nr", ResetNRForDNA);
	if(bpp->dna_db & BLAST_NR)
		SetStatus(bb, TRUE);
	bb = CheckBox(g, "est", ResetESTForDNA);
	if(bpp->dna_db & BLAST_EST)
		SetStatus(bb, TRUE);
	bb = CheckBox(g, "sts", ResetSTSForDNA);
	if(bpp->dna_db & BLAST_STS)
		SetStatus(bb, TRUE);
	bb = CheckBox(g, "month", ResetMonthForDNA);
	if(bpp->dna_db & BLAST_MONTH)
		SetStatus(bb, TRUE);
	bb = CheckBox(g, "htgs", ResetHTGSForDNA);
	if(bpp->dna_db & BLAST_HTGS)
		SetStatus(bb, TRUE);
	bb = CheckBox(g, "vector", ResetVectorForDNA);
	if(bpp->dna_db & BLAST_VECTOR)
		SetStatus(bb, TRUE);
	bb = CheckBox(g, "mito", ResetMitoForDNA);
	if(bpp->dna_db & BLAST_MITO)
		SetStatus(bb, TRUE);
	bb = CheckBox(g, "kabat", ResetKabatForDNA);
	if(bpp->dna_db & BLAST_KABAT)
		SetStatus(bb, TRUE);
	bb = CheckBox(g, "epd", ResetEPDForDNA);
	if(bpp->dna_db & BLAST_EPD)
		SetStatus(bb, TRUE);
	bb = CheckBox(g, "pdb", ResetPDBForDNA);
	if(bpp->dna_db & BLAST_PDB)
		SetStatus(bb, TRUE);
	bb = CheckBox(g, "yeast", ResetYeastForDNA);
	if(bpp->dna_db & BLAST_YEAST)
		SetStatus(bb, TRUE);
	bb = CheckBox(g, "gss", ResetGSSForDNA);
	if(bpp->dna_db & BLAST_GSS)
		SetStatus(bb, TRUE);
	bb = CheckBox(g, "alu", ResetAluForDNA);
	if(bpp->dna_db & BLAST_ALU)
		SetStatus(bb, TRUE);
	StaticPrompt (g, "other", 0, 
		dialogTextHeight, programFont, '1');
	DialogText (g, bpp->other_dna, 5, ResetOtherNucDB);

	gg = HiddenGroup(tg, 2, 0, NULL);
	SetGroupSpacing(gg, 5, 0);
	g = HiddenGroup(gg, 0, 3, NULL);
	SetGroupSpacing(g, 0, 10);
	StaticPrompt (g, "Program", 0, 
		dialogTextHeight, programFont, '1');
	bb = CheckBox(g, "BLASTN", ResetBLASTN);
	if(bpp->blast_program & BLASTN_PROGRAM)
		SetStatus(bb, TRUE);
	bb = CheckBox(g, "TBLASTN", ResetTBLASTN);
	if(bpp->blast_program & TBLASTN_PROGRAM)
		SetStatus(bb, TRUE);

	g = HiddenGroup(gg, 0, 3, NULL);
	StaticPrompt (g, "Parameter", 0, 
		dialogTextHeight, programFont, '1');
	DialogText (g, bpp->n_param, 15, ResetBlastNParameter);
	DialogText (g, bpp->t_param, 15, ResetTBlastNParameter); 
	

	tg = NormalGroup(pw, 0, 2, "Search Protein DataBase", systemFont, NULL);
	g = HiddenGroup(tg, 3, 0, NULL);
	SetGroupSpacing(g, 5, 0);
	bb = CheckBox(g, "nr", ResetNRForProt);
	if(bpp->prot_db & BLAST_NR)
		SetStatus(bb, TRUE);
	bb = CheckBox(g, "pdb", ResetPDBForProt);
	if(bpp->prot_db & BLAST_PDB)
		SetStatus(bb, TRUE);
	bb = CheckBox(g, "swissprot", ResetSwissProtForProt);
	if(bpp->prot_db & BLAST_SWISSPROT)
		SetStatus(bb, TRUE);
	bb = CheckBox(g, "month", ResetMonthForProt);
	if(bpp->prot_db & BLAST_MONTH)
		SetStatus(bb, TRUE);
	bb = CheckBox(g, "kabat", ResetKabatForProt);
	if(bpp->prot_db & BLAST_KABAT)
		SetStatus(bb, TRUE);
	bb = CheckBox(g, "yeast", ResetYeastForProt);
	if(bpp->prot_db & BLAST_YEAST)
		SetStatus(bb, TRUE);
	bb = CheckBox(g, "alu", ResetAluForProt);
	if(bpp->prot_db & BLAST_ALU)
		SetStatus(bb, TRUE);
	StaticPrompt (g, "other", 0, 
		dialogTextHeight, programFont, '1');
	DialogText(g, bpp->other_protein, 5, ResetOtherProtDB);


	gg = HiddenGroup(tg, 2, 0, NULL);
	SetGroupSpacing(gg, 5, 0);
	g = HiddenGroup(gg, 0, 3, NULL);
	SetGroupSpacing(g, 0, 10);
	StaticPrompt (g, "Program", 0, 
		dialogTextHeight, programFont, '1');
	bb = CheckBox(g, "BLASTX", ResetBLASTX);
	if(bpp->blast_program & BLASTX_PROGRAM)
		SetStatus(bb, TRUE);
	bb = CheckBox(g, "BLASTP", ResetBLASTP);
	if(bpp->blast_program & BLASTP_PROGRAM)
		SetStatus(bb, TRUE);

	g = HiddenGroup(gg, 0, 3, NULL);
	StaticPrompt (g, "Parameter", 0, 
		dialogTextHeight, programFont, '1');
	DialogText (g, bpp->x_param, 15, ResetBlastXParameter);
	DialogText (g, bpp->p_param, 15, ResetBlastPParameter); 


	/*set up t search the protein database */
	/* g = NormalGroup(pw, 4, 0, "Search Protein DataBase", systemFont, NULL);
	SetGroupSpacing(g, 5, 0);
	bb = CheckBox(g, "nr", ResetNRForProt);
	if(bpp->prot_db & BLAST_NR)
		SetStatus(bb, TRUE);
	bb = CheckBox(g, "pdb", ResetPDBForProt);
	if(bpp->prot_db & BLAST_PDB)
		SetStatus(bb, TRUE);
	bb = CheckBox(g, "swissprot", ResetSwissProtForProt);
	if(bpp->prot_db & BLAST_SWISSPROT)
		SetStatus(bb, TRUE);
	bb = CheckBox(g, "month", ResetMonthForProt);
	if(bpp->prot_db & BLAST_MONTH)
		SetStatus(bb, TRUE);
	bb = CheckBox(g, "kabat", ResetKabatForProt);
	if(bpp->prot_db & BLAST_KABAT)
		SetStatus(bb, TRUE);
	StaticPrompt (g, "other", 0, 
		dialogTextHeight, systemFont, '1');
	DialogText(g, bpp->other_protein, 10, ResetOtherProtDB);
	Break(pw);

	gg = NormalGroup(pw, 2, 0, "BLAST Program and Parameter", systemFont, NULL);
	SetGroupSpacing(gg, 5, 0);
	g = HiddenGroup(gg, 0, 2, NULL);
	SetGroupSpacing(g, 0, 15);
	bb = CheckBox(g, "BLASTP", ResetBLASTP);
	if(bpp->blast_program & BLASTP_PROGRAM)
		SetStatus(bb, TRUE);
	bb = CheckBox(g, "BLASTX", ResetBLASTX);
	if(bpp->blast_program & BLASTX_PROGRAM)
		SetStatus(bb, TRUE);
	g = HiddenGroup(gg, 0, 2, NULL);
	DialogText (g, bpp->p_param, 20, ResetTBlastNParameter);
	DialogText (g, bpp->x_param, 20, ResetBlastXParameter); */

	g = HiddenGroup(pw, 3, 0, NULL);
	SetGroupSpacing(g, 60, 0);
	PushButton (g, "Cancel", CancelSettingProc);
	PushButton(g, "Help", HelpSettingProc);
	DefaultButton(g, "Accept", AcceptSettingProc);

	SetObjectExtra(pw, pbop, NULL);
	Show(pw);
}




static void SetGapAlignment(GrouP g)
{
	PBlastOptionPtr pbop;
	WindoW w;
	Int2 value;

	w = ParentWindow(g);
	pbop = (PBlastOptionPtr)GetObjectExtra(w);
	if(pbop == NULL)
		return;
	value = GetValue(g);
	if(value > 0)
		pbop->gap_alignment = (Uint1)(value -1);
}

static void SetOrgFilter (GrouP g)
{
	PBlastOptionPtr pbop;
	WindoW w;
	Int2 value;

	value = GetValue(g);
	if(value >0 && value < 4)
	{
		w = ParentWindow(g);
		pbop = (PBlastOptionPtr)GetObjectExtra(w);
		if(pbop != NULL)
			pbop->filter_org = (Uint1)(value -1);
	}
}


static void CheckDustProc(ButtoN b)
{
	PBlastOptionPtr pbop;
	WindoW w;

	w = ParentWindow(b);
	pbop = (PBlastOptionPtr)GetObjectExtra(w);

	if(pbop != NULL)
		pbop->dust = 1 - pbop->dust;
}

static void CheckQueryProc (ButtoN b)
{
	PBlastOptionPtr pbop;
	WindoW w;

	w = ParentWindow(b);
	pbop = (PBlastOptionPtr)GetObjectExtra(w);

	if(pbop != NULL)
		pbop->filter_self = 1 - pbop->filter_self;
}

static void SetOrgName (TexT t)
{
	
	PBlastOptionPtr pbop;
	WindoW w;
	Char buf[201];

	w = ParentWindow(t);
	pbop = (PBlastOptionPtr)GetObjectExtra(w);

	if(pbop != NULL)
	{
		buf[0] = '\0';
		GetTitle(t, buf, 200);
		if(pbop->organism != NULL)
			pbop->organism = MemFree(pbop->organism);
		if(buf[0] != '\0')
			pbop->organism = StringSave(buf);
	}
}

static void ResetOutputFormat (ButtoN b, Uint1 choice)
{
	PBlastOptionPtr pbop;
	WindoW w;
	Uint1 value;

	w = ParentWindow(b);
	pbop = (PBlastOptionPtr)GetObjectExtra(w);
	if(pbop != NULL)
	{
		switch (choice)
		{
			case 1:
				value = OUTPUT_TEXT;
				break;
			case 2:
				value = OUTPUT_SEQALIGN;
				break;
			case 3:
				value = OUTPUT_SEQENTRY;
				break;
			case 4:
				value = OUTPUT_HTML;
				break;
			default:
				value = 0;
				break;
		}
		if(pbop->output_format & value)
			pbop->output_format -= value;
		else
			pbop->output_format |= value;
	}
}

static void ResetOutputText (ButtoN b)
{
	ResetOutputFormat(b, 1);
}
	
static void ResetOutputSeqAlign(ButtoN b)
{
	ResetOutputFormat(b, 2);
}
	
static void ResetOutputSeqEntry(ButtoN b)
{
	ResetOutputFormat(b, 3);
}

static void ResetOutputHTML (ButtoN b)
{
	ResetOutputFormat(b, 4);
}
	
	
static void ResetOutputPath (TexT t)
{
	PBlastOptionPtr pbop;
	WindoW w;

	w = ParentWindow(t);
	pbop = (PBlastOptionPtr)GetObjectExtra(w);

	if(pbop != NULL)
	{
		pbop->output_path[0] = '\0';
		GetTitle(t, pbop->output_path, PATH_MAX -1);
	}
}

static void SaveSettingToConfig (ButtoN b)
{
	PBlastOptionPtr pbop;
	WindoW w;

	w = ParentWindow(b);
	pbop = (PBlastOptionPtr)GetObjectExtra(w);
	if(pbop != NULL)
		load_pboption_to_configure(pbop);
}

static void GetSeqFormat (PopuP p)
{
	PBlastOptionPtr pbop;
	WindoW w;

	w = ParentWindow(p);
	pbop = (PBlastOptionPtr)GetObjectExtra(w);
	if(pbop != NULL)
		pbop->seq_format = (Uint1)GetValue(p);
}

static void LoadSeqData ( TexT t)
{
	PBlastOptionPtr pbop;
	WindoW w;

	w = ParentWindow(t);
	pbop = (PBlastOptionPtr)GetObjectExtra(w);
	if(pbop != NULL)
	{
		if(pbop->seq_data == NULL)
			MemFree(pbop->seq_data);
		pbop->seq_data = SaveStringFromText(t);
		if(pbop->seq_data != NULL)
			Enable(search_button);
		else if(pbop->file_name[0] == '\0')
			Disable(search_button);
	}
}


static void SearchFunc(ButtoN b)
{
	PBlastOptionPtr pbop;
	WindoW w;
	Uint1 retval;
	

	w = ParentWindow(b);
	pbop = (PBlastOptionPtr)GetObjectExtra(w);
	if(pbop != NULL)
	{
		Hide(w);
		/* retval = RunPowerBlast(pbop, TRUE); */
		retval = RunPowerBlast(pbop); 
		if(retval == POWBLAST_FATAL)
		{
			Remove(w);
			BioseqFetchDisable();
			EntrezBioseqFetchDisable();
			BlastFini();
			QuitProgram();
		}
		else
			Show(w);
	}
}


static void CheckMonitorProc (ButtoN b)
{
	PBlastOptionPtr pbop;
	WindoW w;

	w = ParentWindow(b);
	pbop = (PBlastOptionPtr)GetObjectExtra(w);
	if(pbop != NULL)
		pbop->monitor = 1- pbop->monitor;
}

	
static void ChangeFeatureProc (IteM i)
{
	PBlastOptionPtr pbop;
	WindoW w;

	w = ParentWindow(i);
	pbop = (PBlastOptionPtr)GetObjectExtra(w);
	if(pbop != NULL)
		pbop->hide_alignment_feature = 1- pbop->hide_alignment_feature;
}

static void OpenErrorLogProc (IteM i)
{
	PBlastOptionPtr pbop;
	WindoW w;
	Char f_name[101];

	w = ParentWindow(i);
	pbop = (PBlastOptionPtr)GetObjectExtra(w);
	if(pbop->errfp != NULL)
		Disable(i);
	else
	{
		if(GetOutputFileName(f_name, sizeof(f_name), ""))	
		{
			pbop->errfp = FileOpen(f_name, "w");
			if(pbop->errfp == NULL)
				ErrPostEx(SEV_ERROR, 0, 0, "File to Open Error Log File %s", pbop->errfp);
			else
				Disable(i);
		}
	}
}
				
	

/*********************************************************************
*	"main" function to call blast for the client.
*	This function sets up the customized vibrant interface
*	for BLAST.  The actual call to BlastBioseq is with
*	PerformBlastSearch.
*********************************************************************/

Int2 Main (void)

{
	WindoW               w;
	MenU     	         m;
	GrouP                g, gg;
	ButtoN               b;
	PBlastOptionPtr      pbop;
	PopuP			p;
	ObjMgrPtr            omp;
	IteM			i;

	if(!UseLocalAsnloadDataAndErrMsg ())
	{
		if(!SeqEntryLoad())
			return 1;
	}

	if(!EntrezBioseqFetchEnable ("powblast", TRUE))
	{
		ErrPostEx(SEV_ERROR, 0, 0, "Unable to initialize Entrez service");
		return 1;
	}

	BioseqFetchInit(TRUE);

	if (! BlastInit("PowerBlast", FALSE))
	{
		ErrPostEx(SEV_ERROR, 0, 0, "Unable to initialize BLAST service");
		EntrezBioseqFetchDisable();
		return 1;
	}

	omp = ObjMgrGet();
	omp->maxtemp = 60;

	pbop = MemNew(sizeof(PBlastOption));
	pbop->file_name[0] = '\0';
	pbop->seq_data = NULL;
	load_default_blast_param(pbop);

	w = FixedWindow(-50, -33, -10, -10, "PowerBlast Search 1.01", CloseMainWinProc);
	m = PulldownMenu (w, "File");
	CommandItem(m, "Save Error Log", OpenErrorLogProc);
	CommandItem (m, "Quit", ItemQuitProc);
	m = PulldownMenu (w, "Edit");
	CommandItem(m, "Cut", StdCutTextProc);
	CommandItem(m, "Copy", StdCopyTextProc);
	CommandItem(m, "Paste", StdPasteTextProc);
	CommandItem(m, "Clear", StdDeleteTextProc);

	m = PulldownMenu(w, "Option");
	i = StatusItem(m, "ShowFeature", ChangeFeatureProc);
	SetStatus(i, pbop->hide_alignment_feature? FALSE : TRUE);
	


	g = HiddenGroup (w, 3, 0, NULL);
	b = PushButton (g, "Read Input File", GetInputForBlast);
	SetObjectExtra(b, pbop, NULL);
	input_text = DialogText (g, pbop->file_name, 20, FileInForInput);
	PushButton (g, "Blast Program", SetBlastParamProc);

	g = HiddenGroup (w, 4, 0, NULL);
	SetGroupSpacing(g, 5, 0);
	StaticPrompt (g, "Or Paste Query Formated As:", 0, 
		dialogTextHeight, systemFont, '1');
	p = PopupList (g, TRUE, GetSeqFormat);
	PopupItem(p, "FASTA");
	PopupItem(p, "Accession");
	PopupItem(p, "GI");
	SetValue(p, pbop->seq_format);
	PushButton (g, "Clear Window", ClearSeqInputWindow);
	b = CheckBox(g, "Use Monitor", CheckMonitorProc);
	SetStatus(b, pbop->monitor);
	
	sequence = ScrollText (w, 35, 6, programFont, TRUE, LoadSeqData);

	g = HiddenGroup(w, 2, 0, NULL);
	PushButton(g, "Mask Repeats", GetInputForRepeat);
	if(pbop->repeat_library[0] != '\0')
		repeat_text = DialogText (g, pbop->repeat_library, 20, FileInForRepeat);
	else
		repeat_text = DialogText (g, " ", 20, FileInForRepeat);

	gg = HiddenGroup(w, 2, 0, NULL);
	SetGroupSpacing(gg, 20, 0);
	g = NormalGroup(gg, 4, 0, "Gapped Alignment Algorithm", systemFont, SetGapAlignment);
	RadioButton(g, "None");
	RadioButton(g, "SIM");
	RadioButton(g, "SIM2");
	RadioButton(g, "SIM3");
	SetValue(g, (Int2)(pbop->gap_alignment +1));

	g = NormalGroup(gg, 2, 0, "Filter", systemFont, NULL);
	b = CheckBox(g, "Low Complexity", CheckDustProc);
	SetStatus(b, pbop->dust); 
	b = CheckBox(g, "Self Hit", CheckQueryProc);
	SetStatus(b, pbop->filter_self);


	g = NormalGroup(w, 4, 0, "Organism Specific Search", systemFont, SetOrgFilter);
	RadioButton(g, "None");
	RadioButton(g, "Include");
	RadioButton(g, "Exclude");
	SetValue(g, (Int2)(pbop->filter_org +1));
	DialogText (g, pbop->organism, 20, SetOrgName);

	/* gg = HiddenGroup(w, 2, 0, NULL); */
	gg = NormalGroup(w, 0, 2, "Export Output File ", systemFont, NULL);
	/*directory for the output file*/
	g = HiddenGroup(gg, 2, 0, NULL);
	StaticPrompt(g, "Path", 0, dialogTextHeight, systemFont, '1');
	if(pbop->output_path[0] != '\0')
		DialogText (g, pbop->output_path, 12, ResetOutputPath);
	else
		DialogText (g, " ", 12, ResetOutputPath);

	/*format for the output file */
	g = HiddenGroup(gg, 4, 0, NULL);
	SetGroupSpacing(g, 10, 0);
	/* StaticPrompt(g, "Format for Export File", 0, dialogTextHeight, systemFont, '1'); */
	b = CheckBox(g, "Text(*.ali)", ResetOutputText);
	SetStatus(b, (pbop->output_format & OUTPUT_TEXT));
	b = CheckBox(g, "HTML(*.html)", ResetOutputHTML);
	SetStatus(b, (pbop->output_format & OUTPUT_HTML));
	b = CheckBox(g, "Seq-align(*.sat)", ResetOutputSeqAlign);
	SetStatus(b, (pbop->output_format & OUTPUT_SEQALIGN));
	b = CheckBox(g, "Seq-entry(*.ent)", ResetOutputSeqEntry);
	SetStatus(b, (pbop->output_format & OUTPUT_SEQENTRY));
	
	
	g = HiddenGroup(w, 3, 0, NULL);
	SetGroupSpacing(g, 120, 0);
	PushButton(g, "Save Setting", SaveSettingToConfig);
	search_button = DefaultButton(g, "Search", SearchFunc);
	Disable(search_button);
	PushButton(g, "Quit", ButtonQuitProc);

	SetObjectExtra(w, pbop, FreePBlastOptionData);
	Show (w); 
	ProcessEvents();

	return 0;
}

