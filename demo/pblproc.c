#include <pobutil.h>
#include <txalign.h>
#include <pblproc.h>
#include <netblap2.h>
#include <accutils.h>
/* #include <getfa.h> */

#define MAX_GI_NUM 1000
#define MIN_GAP_LEN 10	/*for converting the Seq-align into feature interval*/



			

/********************************************************
*
*	free the settings in the blast parameters
*
********************************************************/
void free_blast_param(BlastParamPtr bpp)
{
	if(bpp->n_param != NULL)
		MemFree(bpp->n_param);
	if(bpp->p_param != NULL)
		MemFree(bpp->p_param);
	if(bpp->t_param != NULL)
		MemFree(bpp->t_param);
	if(bpp->x_param != NULL)
		MemFree(bpp->x_param);
	if(bpp->other_dna != NULL)
		MemFree(bpp->other_dna);
	if(bpp->other_protein != NULL)
		MemFree(bpp->other_protein);

	MemFree(bpp);
}

void FreePBlastOption(PBlastOptionPtr pbop)
{
	if(pbop == NULL)
		return;
	if(pbop->seq_data != NULL)
		MemFree(pbop->seq_data);
	if(pbop->bpp != NULL)
		free_blast_param (pbop->bpp);
	if(pbop->temp_bpp != NULL)
		free_blast_param(pbop->temp_bpp);
	if(pbop->organism)
		MemFree(pbop->organism);
	if(pbop->aa_gi_bsp)
		BSFree(pbop->aa_gi_bsp);
	if(pbop->na_gi_bsp)
		BSFree(pbop->na_gi_bsp);
	if(pbop->alu_sep_list || pbop->alu_slp_list)
		free_alu_list(pbop->alu_sep_list, pbop->alu_slp_list);
	if(pbop->output_file != NULL)
		ValNodeFreeData(pbop->output_file);

	MemFree(pbop);
}


void load_default_blast_param(PBlastOptionPtr pbop)
{
	Char buf[101];
	BlastParamPtr bpp;

	pbop->bpp = MemNew(sizeof(BlastParameter));
	bpp = pbop->bpp;
	/*load the pre-configured blast program parameter and database*/
	

	/*		blastn		*/
	if(GetAppParam ("POWBLAST", "BLASTN", "SEARCH", "FALSE", buf, sizeof (buf)))
	{
		if(StringCmp(buf, "TRUE") == 0)
			bpp->blast_program |= BLASTN_PROGRAM;
	}
	buf[0] = '\0';
	if(GetAppParam ("POWBLAST", "BLASTN", "PARAMETER", "", buf, sizeof (buf)))
	{
		if(buf[0] != '\0')
			bpp->n_param = StringSave(buf);
	}

	buf[0] = '\0';
	if(GetAppParam ("POWBLAST", "BLASTN", "MAXLEN", "", buf, sizeof (buf)))
	{
		if(buf[0] != '\0')
			bpp->max_blast_n = MIN(MAX_BLASTN_LEN, atol(buf));
		else
			bpp->max_blast_n = MAX_BLASTN_LEN;
	}
	else
		bpp->max_blast_n = MAX_BLASTN_LEN;


	/*		blastp		*/
	if(GetAppParam ("POWBLAST", "BLASTP", "SEARCH", "FALSE", buf, sizeof (buf)))
	{
		if(StringCmp(buf, "TRUE") == 0)
			bpp->blast_program |= BLASTP_PROGRAM;
	}
	buf[0] = '\0';
	if(GetAppParam ("POWBLAST", "BLASTP", "PARAMETER", "", buf, sizeof (buf)))
	{
		if(buf[0] != '\0')
			bpp->p_param = StringSave(buf);
	}

	buf[0] = '\0';
	if(GetAppParam ("POWBLAST", "BLASTP", "MAXLEN", "", buf, sizeof (buf)))
	{
		if(buf[0] != '\0')
			bpp->max_blast_p = MIN(MAX_BLASTP_LEN, atol(buf));
		else
			bpp->max_blast_p = MAX_BLASTN_LEN;
	}
	else
		bpp->max_blast_p = MAX_BLASTN_LEN;

	/*		tblastn			*/
	if(GetAppParam ("POWBLAST", "TBLASTN", "SEARCH", "FALSE", buf, sizeof (buf)))
	{
		if(StringCmp(buf, "TRUE") == 0)
			bpp->blast_program |= TBLASTN_PROGRAM;
	}
	buf[0] = '\0';
	if(GetAppParam ("POWBLAST", "TBLASTN", "PARAMETER", "", buf, sizeof (buf)))
	{
		if(buf[0] != '\0')
			bpp->t_param = StringSave(buf);
	}

	buf[0] = '\0';
	if(GetAppParam ("POWBLAST", "TBLASTN", "MAXLEN", "", buf, sizeof (buf)))
	{
		if(buf[0] != '\0')
			bpp->max_blast_t = MIN(MAX_TBLASTN_LEN, atol(buf));
		else
			bpp->max_blast_t = MAX_TBLASTN_LEN;
	}
	else
		bpp->max_blast_t = MAX_TBLASTN_LEN;


	/*		blastx			*/
	if(GetAppParam ("POWBLAST", "BLASTX", "SEARCH", "FALSE", buf, sizeof (buf)))
	{
		if(StringCmp(buf, "TRUE") == 0)
			bpp->blast_program |= BLASTX_PROGRAM;
	}
	buf[0] = '\0';
	if(GetAppParam ("POWBLAST", "BLASTX", "PARAMETER", "", buf, sizeof (buf)))
	{
		if(buf[0] != '\0')
			bpp->x_param = StringSave(buf);
	}

	buf[0] = '\0';
	if(GetAppParam ("POWBLAST", "BLASTX", "MAXLEN", "", buf, sizeof (buf)))
	{
		if(buf[0] != '\0')
			bpp->max_blast_x = MIN(MAX_BLASTX_LEN, atol(buf));
		else
			bpp->max_blast_x = MAX_BLASTX_LEN;
	}
	else
		bpp->max_blast_x = MAX_BLASTX_LEN;


	/*		the database 		*/
	if(GetAppParam ("POWBLAST", "DATABASE", "DNA_DB", "1", buf, sizeof (buf)))
	{
		bpp->dna_db= atol(buf);
	}
	else
		bpp->dna_db= 1;
	if(GetAppParam ("POWBLAST", "DATABASE", "PROT_DB", "1", buf, sizeof (buf)))
		bpp->prot_db= atol(buf);
	else
		bpp->prot_db= 1;
	
	buf[0] = '\0';
	if(GetAppParam ("POWBLAST", "DATABASE", "OTHER_DNA", "", buf, sizeof (buf)))
	{
		if(buf[0] != '\0')
			bpp->other_dna = StringSave(buf);
	}

	buf[0] = '\0';
	if(GetAppParam ("POWBLAST", "DATABASE", "OTHER_PROTEIN", "", buf, sizeof (buf)))
	{
		if(buf[0] != '\0')
			bpp->other_protein = StringSave(buf);
	}

	buf[0] = '\0';
	if(GetAppParam ("POWBLAST", "FILTER", "REPEAT", "", buf, sizeof (buf)))
	{
		if(buf[0] != '\0')
			StringCpy(pbop->repeat_library, buf);
		else
			pbop->repeat_library[0] = '\0';
	}
	else
		pbop->repeat_library[0] = '\0';

	if(GetAppParam ("POWBLAST", "FILTER", "DUST", "TRUE", buf, sizeof (buf)))
	{
		if(StringCmp(buf, "TRUE") == 0)
			pbop->dust = TRUE;
	}
	else
		pbop->dust = TRUE;

	if(GetAppParam ("POWBLAST", "FILTER", "QUERY", "TRUE", buf, sizeof (buf)))
	{
		if(StringCmp(buf, "TRUE") == 0)
			pbop->filter_self = TRUE;
	}

	if(GetAppParam ("POWBLAST", "FILTER", "ORGANISM", "0", buf, sizeof (buf)))
		pbop->filter_org = (Uint1)atoi(buf);

	buf[0] = '\0';
	if(GetAppParam ("POWBLAST", "FILTER", "ORG_NAME", "", buf, sizeof (buf)))
	{
		if(buf[0] != '\0')
			pbop->organism = StringSave(buf);
	}

	if(GetAppParam ("POWBLAST", "ALIGN", "SIM", "0", buf, sizeof (buf)))
		pbop->gap_alignment = (Uint1)atoi(buf);

	if(GetAppParam ("POWBLAST", "OUTPUT", "FORMAT", "1", buf, sizeof (buf)))
		pbop->output_format= (Uint1)atoi(buf);
	else
		pbop->output_format = 1;

	pbop->output_path[0] = '\0';
	GetAppParam ("POWBLAST", "OUTPUT", "PATH", "", pbop->output_path, sizeof (pbop->output_path));

	pbop->seq_format = SEQFMT_FASTA;
	if(GetAppParam ("POWBLAST", "SEQWIN", "FORMAT", "1", buf, sizeof (buf)))
		pbop->seq_format = MAX(1, (Uint1)atoi (buf));
	

	if(GetAppParam ("POWBLAST", "RUNTIME", "MONITOR", "TRUE", buf, sizeof (buf)))
	{
		if(StringCmp(buf, "TRUE") == 0)
			pbop->monitor = TRUE;
		else if(StringCmp(buf, "FALSE") == 0)
			pbop->monitor = FALSE;
	}
	else
		pbop->monitor = TRUE;
		
}

/****************************************************************
*
*	write the current parameter to the configure file 
*	.powblastrc
*
*
****************************************************************/
Boolean load_pboption_to_configure(PBlastOptionPtr pbop)
{
	

	Uint1 blast_program;	/*0 for default */
	Char section[201];
	CharPtr param;
	Char value[201];
	BlastParamPtr bpp;
	Int4 max_len;

	if(pbop == NULL)
		return FALSE;
	if(pbop->bpp != NULL)
	{
		bpp = pbop->bpp;
		for(blast_program = BLASTN_PROGRAM; 
			blast_program <= TBLASTN_PROGRAM; blast_program *=2)
		{
			param = NULL;
			switch(blast_program)
			{
				case BLASTN_PROGRAM:
					StringCpy(section, "BLASTN");
					param = pbop->bpp->n_param;
					max_len = pbop->bpp->max_blast_n;
					break;

				case TBLASTN_PROGRAM:
					StringCpy(section, "TBLASTN");
					param = pbop->bpp->t_param;
					max_len = pbop->bpp->max_blast_t;
					break;

				case BLASTP_PROGRAM:
					StringCpy(section, "BLASTP");
					param = pbop->bpp->p_param;
					max_len = pbop->bpp->max_blast_p;
					break;

				case BLASTX_PROGRAM:
					StringCpy(section, "BLASTX");
					param = pbop->bpp->x_param;
					max_len = pbop->bpp->max_blast_x;
					break;

				default:
					break;
			}
			if(pbop->bpp->blast_program & blast_program)
				SetAppParam ("POWBLAST", section, "SEARCH", "TRUE");
			else
				SetAppParam ("POWBLAST", section, "SEARCH", "FALSE");
			SetAppParam ("POWBLAST", section, "PARAMETER", param);
			sprintf(value, "%ld", max_len);
			SetAppParam ("POWBLAST", section, "MAXLEN", value);
		}

		/* Setting up the databases */
		sprintf(value, "%ld", bpp->dna_db);
		SetAppParam("POWBLAST", "DATABASE", "DNA_DB", value);
		SetAppParam("POWBLAST", "DATABASE", "OTHER_DNA", bpp->other_dna);

		sprintf(value, "%ld", bpp->prot_db);
		SetAppParam("POWBLAST", "DATABASE", "PROT_DB", value);
		SetAppParam("POWBLAST", "DATABASE", "OTHER_PROTEIN", bpp->other_dna);

	}

	/*Setting up the filter functions */
	SetAppParam ("POWBLAST", "FILTER", "REPEAT", pbop->repeat_library);
	if(pbop->dust)
		StringCpy(value, "TRUE");
	else
		StringCpy(value, "FALSE");
	SetAppParam ("POWBLAST", "FILTER", "DUST", value);

	if(pbop->filter_self)
		StringCpy(value, "TRUE");
	else
		StringCpy(value, "FALSE");
	SetAppParam ("POWBLAST", "FILTER", "QUERY", value);

	sprintf(value, "%d", pbop->filter_org);
	SetAppParam ("POWBLAST", "FILTER", "ORGANISM", value);
	SetAppParam ("POWBLAST", "FILTER", "ORG_NAME", pbop->organism);

	sprintf(value, "%d", pbop->gap_alignment);
	SetAppParam ("POWBLAST", "ALIGN", "SIM", value);

	sprintf(value, "%d", pbop->output_format);
	SetAppParam ("POWBLAST", "OUTPUT", "FORMAT", value);
	SetAppParam ("POWBLAST", "OUTPUT", "PATH", pbop->output_path);

	sprintf(value, "%d", pbop->seq_format);
	SetAppParam ("POWBLAST", "SEQWIN", "FORMAT", value);

	if(pbop->monitor)
		StringCpy(value, "TRUE");
	else
		StringCpy(value, "FALSE");
	SetAppParam ("POWBLAST", "RUNTIME", "MONITOR", value);


	return TRUE;
		
}



/***********************************************************
*
*	check_FASTA_valid() check to see if the FASTA sequence 
*	is valid or not and figure out the type of the FASTA 
*	file
*
************************************************************/
static Boolean check_FASTA_valid (CharPtr str, Uint1Ptr p_program, BoolPtr is_aa)
{
	Int2 ck_len;
	Int4 max_len, j;
	Int2 non_DNA;
	Int2 bad_residue;
	Int2 empty_residue;
	Uint1 program;


	max_len = StringLen(str);
	non_DNA = 0;
	bad_residue = 0;
	empty_residue = 0;
	j = 0;
	ck_len = 0;
	while(ck_len <100 && j < max_len )
	{
		++j;
		if(!IS_ALPHA(*str))
		{
			if(*str == '-' || *str == '*')
			{
				++non_DNA;
				++ck_len;
			}
			else
			{
				 if(!IS_WHITESP(*str) && ! IS_DIGIT(*str))
				{
					++bad_residue;
				}
			}
		}
		else
		{
			++ck_len;
			if(StrChr("ACGTNacgtn", *str)==NULL)
			{
				++non_DNA;
			}
		}
		++str;
	} 

	if(ck_len == 0)
	{
		ErrPostEx(SEV_ERROR, 0, 0, "Empty sequence data");
		return FALSE;
	}

	if(bad_residue > 0)
	{
		if(ck_len / bad_residue <=2)
		{
			ErrPostEx(SEV_ERROR, 0, 0, "Two many bad residues");
			return FALSE;
		}
	}

	*is_aa = (non_DNA >=ck_len/4);
	program = *p_program;
	if(program == BLAST_DEFAULT)
	{
		if(*is_aa)
			*p_program = BLASTP_PROGRAM;
		else
			*p_program = BLASTN_PROGRAM;
	}
	
	if(program & BLASTN_PROGRAM || program & BLASTX_PROGRAM)
	{
		if(*is_aa == FALSE || 
			(!(program & BLASTP_PROGRAM) &&!(program & TBLASTN_PROGRAM)))
		{
			if(ck_len < MIN_BLASTN_LEN)
			{
				ErrPostEx(SEV_ERROR, 0, 0, "The search length %d is too short", ck_len);
				return FALSE;
			}
			*is_aa = FALSE;
			return TRUE;
		}
	}

	if(program & BLASTP_PROGRAM || program & program & TBLASTN_PROGRAM)
	{
		*is_aa = TRUE;
		if(ck_len < MIN_BLASTP_LEN)
		{
			ErrPostEx(SEV_ERROR, 0, 0, "The search length %ld is too short", ck_len);
			return FALSE;
		}
		return TRUE;
	}

	return FALSE;
}

static Boolean check_fasta_seq(CharPtr buf, BoolPtr is_aa, Uint1Ptr program)
{
	Int2 length, i;
	CharPtr str;
	Boolean stop;

	str = buf;
	length = 2000;
	if(*str == '>')
	{
		i = 0;
		stop = FALSE;
		while(!stop && str && i<length)
		{
			if(*str == '\n' || *str == '\r')
			{
				stop = TRUE;
				break;
			}
			++str;
			++i;
		}
		if(!stop)
		{
			ErrPostEx(SEV_ERROR, 0, 0, "Empty sequence data");
			return FALSE;
		}
		++str;
	}


	return check_FASTA_valid (str, program, is_aa);
}

/* return 0 for non-FASTA files, 1 for good FASTA files, 2 for bad FASTA files */
static Uint1 is_fasta_file(FILE *fp, BoolPtr is_aa, Uint1Ptr program)
{
	Char buf[1001];

	FileGets(buf, 1000, fp);
	if(buf[0] == '>')
	{
		FileGets(buf, 1000, fp);	/*read the first line*/
		if(check_FASTA_valid(buf, program, is_aa))
			return 1;
		else
			return 2;
	}
	else
		return 0;
}

static void clean_up_this_session (MonitorPtr mp)
{
	if(mp != NULL)
		MonitorFree(mp);
	/* BlastFini(); */
	SeqMgrFreeCache();
}
	
static Uint1 powblast_one_entry PROTO((BioseqPtr bsp, SeqLocPtr slp, PBlastOptionPtr pbop, MonitorPtr mp));

static Uint1 process_one_bioseq (BioseqPtr bsp, SeqLocPtr slp, PBlastOptionPtr pbop, MonitorPtr mp)
{
	Char temp[201];
	Char seq_name[101];
	Uint1 val;

	if(bsp == NULL || slp == NULL)
		return POWBLAST_NONE;

	seqid_name(bsp->id, seq_name, FALSE, FALSE);
	if(mp != NULL)
	{
		sprintf(temp, "Running PowerBlast for %s", seq_name);
		MonitorStrValue(mp, temp);
	}
	val = powblast_one_entry(bsp, slp,  pbop, mp);
	if(val == POWBLAST_FATAL)
		clean_up_this_session(mp);
	else if(val == POWBLAST_HIT)
	{
		if(mp != NULL)
		{
			sprintf(temp, "Find Hits for %s", seq_name);
			MonitorStrValue(mp, temp);
		}
	}
	else
	{
		if(mp != NULL)
		{
			sprintf(temp, "No Hit for  %s", seq_name);
			MonitorStrValue(mp, temp);
		}
	}

	return val;
}

#define SEPCHARS " \t\r\n"
/* process a data string which can be a FASTA sequence name, gi or acc 
w/o specified locations */
static Uint1 powblast_one_string (CharPtr data, PBlastOptionPtr pbop, MonitorPtr mp)
{
	CharPtr str;
	SeqLocPtr slp;
	BioseqPtr bsp;
	Uint1 retval, val;

	retval = POWBLAST_NONE;
	str = StrTok(data, SEPCHARS);
	while(str)
	{
		slp = prepare_align_data(str, &bsp);
		if(slp != NULL && bsp != NULL)
		{
			val = process_one_bioseq (bsp, slp, pbop, mp);
			BioseqUnlock(bsp);
			SeqLocFree(slp);
			if(val == POWBLAST_FATAL)
				return val;
			else if(val == POWBLAST_HIT)
				retval = POWBLAST_HIT;
		}
		else
			ErrPostEx(SEV_ERROR, 0, 0, "Fail to get data for %s", str);
		str = StrTok(NULL, SEPCHARS);
	}

	return retval;
}

static ByteStorePtr get_gi_list(CharPtr tax_name, Uint1 mol, Int4Ptr count)
{
	Char str[200];
	DocType db;

	db = (mol == Seq_mol_aa) ? TYP_AA : TYP_NT;
	sprintf(str, "\"%s\"[ORGN]", tax_name);
	*count = EntrezTLEvalCountString(str, db, -1, NULL, NULL);
	if(*count > 0 && *count < MAX_GI_NUM)
		return EntrezTLEvalXString(str, db, -1, NULL, NULL);
	else
		return NULL;
}


/***************************************************************
*
*	some powblast operations, like the repeat filtering 
*	and organism filtering will be operated on all the 
*	sequences. InitPowerBlastOption will load the repeat 
*	library and find if there is any errors in organism 
*	field and decide weather the organism filtering is done 
*	with Eval or the other way 
*
***************************************************************/
Boolean InitPowerBlastOption (PBlastOptionPtr pbop)
{
	Int4 count;
	Int2 answer;
	
	
	/*chekc the organism name first */
	/*clean up */
	if(pbop->na_gi_bsp != NULL)
		pbop->na_gi_bsp = BSFree(pbop->na_gi_bsp);
	if(pbop->aa_gi_bsp != NULL)
		pbop->aa_gi_bsp = BSFree(pbop->aa_gi_bsp);

	if(pbop->filter_org && pbop->organism != NULL)
	{

		if(pbop->bpp->blast_program & BLASTN_PROGRAM || 
			pbop->bpp->blast_program & TBLASTN_PROGRAM)
		{
			pbop->na_gi_bsp = get_gi_list(pbop->organism, Seq_mol_dna, &count);
			if(count == 0)
			{
				answer = Message(MSG_YN, "No DNA Record for Organism %s. Proceed?", pbop->organism);
				if(answer == ANS_NO)
					return FALSE;
			}
		}

		if(pbop->bpp->blast_program & BLASTP_PROGRAM || 
				pbop->bpp->blast_program & BLASTX_PROGRAM)
		{
			pbop->aa_gi_bsp = get_gi_list(pbop->organism, Seq_mol_aa, &count);
			if(count == 0)
			{
				answer = Message(MSG_YN, "No Protein Record for %s. Proceed?", pbop->organism);
				if(answer == ANS_NO)
					return FALSE;
			}
		}
	}
	if(pbop->repeat_library[0] != '\0')
	{
		/*clean up the left overs */
		if(pbop->alu_sep_list != NULL)
		{
			free_alu_list(pbop->alu_sep_list, pbop->alu_slp_list);
			pbop->alu_sep_list = NULL;
			pbop->alu_slp_list = NULL;
		}
			
		pbop->alu_sep_list = make_repeat_lib(pbop->repeat_library, 
                                   &(pbop->alu_slp_list), FILE_IO);
		if(pbop->alu_slp_list == NULL)
		{
			answer = Message(MSG_YN, "Fail to Load Repeat Sequence from  %s. Proceed?",  pbop->repeat_library);
			if(answer == ANS_NO)
				return FALSE;
		}
	}

	return TRUE;
}

static void print_seq_entry (SeqEntryPtr sep)
{
	AsnIoPtr aip;
	Char name[101];

	StringCpy(name, "temp.out");
	aip = AsnIoOpen(name, "w");
	SeqEntryAsnWrite(sep, aip, NULL);
	AsnIoClose(aip);

}
/***********************************************************************
*
*	Run one session of powblast. return POWBLAST_NONE for no hit
*	return POWBLAST_HIT to indicate there is a hit
*	return POWBLAST_FATAL to stop the process and indicate there is 
*	fatal error 
*
***********************************************************************/
Uint1 RunPowerBlast(PBlastOptionPtr pbop)
{
	FILE *fp;
	CharPtr StartFasta, NextFasta;
	Boolean is_aa, is_na;
	Uint1 is_fasta;
	Char temp[201], buf[201];
	SeqEntryPtr sep;
	BioseqPtr bsp;
	SeqLocPtr slp;
	Uint1 val, retval;
	MonitorPtr mp;
	Boolean use_monitor;
	Boolean reset_blast;

	if(pbop == NULL)
		return POWBLAST_FATAL;
	use_monitor = pbop->monitor;

	fp = NULL;
	if(pbop->file_name[0] != '\0')
	{
		fp = FileOpen(pbop->file_name, "r");
		if(fp == NULL)
			ErrPostEx(SEV_WARNING, 0, 0, "Fail to Open File %s", pbop->file_name);

	}
	if(fp == NULL && pbop->seq_data == NULL)
	{
		ErrPostEx(SEV_FATAL, 0, 0, "No Input Sequence Data");
		return POWBLAST_NONE;
	}

	/* if (! BlastInit("PowerBlast", FALSE)) {
		ErrPostEx(SEV_ERROR, 0, 0, "Unable to initialize BLAST service");
		if(fp != NULL)
			FileClose(fp);
   		return POWBLAST_FATAL;
  	} */


	/* if(!EntrezBioseqFetchEnable ("powblast", TRUE))
	{
		ErrPostEx(SEV_ERROR, 0, 0, "Unable to initialize Entrez service");
		if(fp != NULL)
			FileClose(fp);
		BlastFini();
   		return POWBLAST_FATAL;
  	} */

	if(!InitPowerBlastOption(pbop))
	{
		/* BlastFini(); */
		return POWBLAST_FATAL;
	}

	/* BioseqFetchInit(TRUE); */

	mp = NULL;
	if(use_monitor)
		mp = MonitorStrNew("Power Blast", 40);
	reset_blast = (pbop->bpp->blast_program == BLAST_DEFAULT);
	if(fp == NULL)
	{	/*read from the buffer */
		switch(pbop->seq_format)
		{
		case SEQFMT_FASTA:
			for(StartFasta = pbop->seq_data; StartFasta!= NULL; )
			{
				if(reset_blast)
					pbop->bpp->blast_program = BLAST_DEFAULT;
				if(mp != NULL)
					MonitorStrValue(mp, "Load Query Sequence");
				sep = NULL;
				if(check_fasta_seq(StartFasta, &is_aa, &(pbop->bpp->blast_program)))
				{
					is_na = (is_aa == FALSE);
					if((sep = FastaToSeqBuff(StartFasta, &NextFasta, is_na)) == NULL)
						ErrPostEx(SEV_ERROR, 0, 0, "Fail In FastaToSeqBuff. Skip");
					else
					{
						bsp = sep->data.ptrvalue;
						slp = SeqLocIntNew(0, bsp->length-1, Seq_strand_both, bsp->id);
						val = process_one_bioseq (bsp, slp, pbop, mp);
						SeqEntryFree(sep);
						SeqLocFree(slp);
						if(val == POWBLAST_FATAL)
							return POWBLAST_NONE;
						else if(val == POWBLAST_HIT)
							retval = POWBLAST_HIT;
					}
				}
				else
				{
					is_na = (is_aa == FALSE);
					sep = FastaToSeqBuff(StartFasta, &NextFasta, is_na);
					ErrPostEx(SEV_ERROR, 0, 0, "Bad FASTA Sequence. Skip");
					if(sep != NULL)
						SeqEntryFree(sep);
				}
					
				StartFasta = NextFasta;
			}
			break;

		default:
			val = powblast_one_string (pbop->seq_data, pbop, mp);
			if(val == POWBLAST_FATAL)
				return POWBLAST_NONE;
			else if(val == POWBLAST_HIT)
				retval = POWBLAST_HIT;
			break;
		}
	}
	else
	{	/*read it from file */
		is_fasta = is_fasta_file(fp, &is_aa, &(pbop->bpp->blast_program));	/*check if it is a fasta file?*/
		is_na = (is_aa == FALSE);
		rewind (fp);
		if(is_fasta == 1)	/*it is a good FASTA file */
		{
			while((sep = FastaToSeqEntry(fp, is_na)) != NULL)
			{
				bsp = sep->data.ptrvalue;
				slp = SeqLocIntNew(0, bsp->length-1, Seq_strand_both, bsp->id);
				val = process_one_bioseq (bsp, slp, pbop, mp);
				SeqEntryFree(sep);
				SeqLocFree(slp);
				if(val == POWBLAST_FATAL)
					return POWBLAST_NONE;
				else if(val == POWBLAST_HIT)
					retval = POWBLAST_HIT;
			}
		}
		else if (is_fasta == 0)
		{
			while(FileGets(temp, 200, fp))
			{
				sscanf(temp, "%s\n", buf);
				val = powblast_one_string (buf, pbop, mp);
				if(val == POWBLAST_FATAL)
				{
					FileClose(fp);
					return POWBLAST_NONE;
				}
				else if(val == POWBLAST_HIT)
					retval = POWBLAST_HIT;
			}
		}
		FileClose(fp);
		if(reset_blast)
			pbop->bpp->blast_program = BLAST_DEFAULT;
	}

	clean_up_this_session (mp);
	return retval;
}


static void link_two_annot(SeqAnnotPtr PNTR head, SeqAnnotPtr new_sap)
{
	SeqAnnotPtr curr;

	if(new_sap == NULL)
		return;
	if(*head == NULL)
		*head = new_sap;
	else
	{
		curr = *head;
		curr->data = link_align((SeqAlignPtr)(new_sap->data), (SeqAlignPtr)(curr->data));
		new_sap->data = NULL;
		SeqAnnotFree(new_sap);
	}
}

static Boolean print_blast_error(SeqAnnotPtr result, BLAST0ResponsePtr brp, Boolean prt_msg)
{
	BLAST0StatusPtr status;
	Char buf[501];

	if(result == NULL)
	{
		if(CheckIfBlastJobCancelled())
			return TRUE;
		while(brp)
		{
			if(brp->choice == BLAST0Response_status)
			{
				status = brp->data.ptrvalue;
				if(status->reason != NULL)
				{
					sprintf(buf, "%s\nEXIT CODE %ld", status->reason, status->code);
					if(prt_msg)
						ErrPostEx(SEV_ERROR, 0, 0, buf);
					if((status->code == 16 ) && StrStr(
						status->reason, "valid contexts"))
					{	/*clean up and re-establish the connection*/
						BlastFini();
						BlastInit("PowerBlast", TRUE);
						return FALSE;
					}
					if((status->code == 5) && StrStr(status->reason, "cpu time"))
					{	/*clean up and re-establish the connection*/
						BlastFini();
						BlastInit("PowerBlast", TRUE);
						return FALSE;
					}
					ErrPostEx(SEV_FATAL, 0, 0,
						"%s\nEXIT CODE %ld", status->reason, status->code);
					return TRUE;
				}
			}
			brp = brp->next;
		}
	}
	return FALSE;

}
			
				
static Uint1 get_blast_type (SeqAnnotPtr annot)
{
	ValNodePtr desc;
	ObjectIdPtr oip;
	UserFieldPtr ufp;
	UserObjectPtr uop;

	if(annot->desc == NULL)
		return 0;
	for(desc = annot->desc; desc != NULL; desc = desc->next)
	{
		if(desc->choice == Annot_descr_user)
		{
			uop = desc->data.ptrvalue;
			oip = uop->type;
			if(oip && oip->str && 
				StringCmp(oip->str, "Blast Type") == 0)
			{
				ufp = uop->data;
				if(ufp && ufp->choice == 2)
					return (Uint1)(ufp->data.intvalue);
			}
		}
	}

	return 0;
}
				
			
static void link_annot_to_end(SeqAnnotPtr PNTR head, SeqAnnotPtr new_sap)
{
	SeqAnnotPtr sap;

	if(*head == NULL)
		*head = new_sap;
	else
	{
		sap = *head;
		while(sap->next != NULL)
			sap = sap->next;
		sap->next = new_sap;
	}
}


/*
*	functions related to filtering the organisms
*
*/

static void print_gi_list(ByteStorePtr bsp)
{
	FILE *fp;
	Int4 gi;
	Int4 i;

	fp = FileOpen("temp.out", "w");
	BSSeek(bsp, 0, SEEK_SET);
	 i = 0;
	while(i < bsp->totlen)
	{
		BSRead(bsp, &gi, 4);
		fprintf(fp, "%ld\n", gi);
		i+=4;
	}
	FileClose(fp);
}


static ByteStorePtr get_filtered_gi_list(SeqAnnotPtr annot, CharPtr organism, BioseqPtr bsp, Uint1 blast_type, Int4 seq_count)
{
	Int4 i, n;
	Int4Ptr gis;
	ValNodePtr gi_list, curr;
	SeqAlignPtr align;
	DenseDiagPtr ddp;
	DenseSegPtr dsp;
	StdSegPtr ssp;
	SeqIdPtr sip;
	Char str[201], name[201];
	ByteStorePtr gi_bsp;
	DocType db;


	n = 0;
	gi_list = NULL;
	align = annot->data;
	while(align)
	{
		sip = NULL;
		switch(align->segtype)
		{
			case 1:
				ddp = align->segs;
				sip = ddp->id->next;
				break;
			case 2:
				dsp = align->segs;
				sip = dsp->ids->next;
				break;
			case 3:
				ssp = align->segs;
				sip = SeqLocId(ssp->loc->next);
				break;
			default:
				break;
		}
		if(sip != NULL && sip->choice == SEQID_GI)
		{
			++n;
			ValNodeAddInt(&gi_list, 0, sip->data.intvalue);
		}
		align = align->next;
	}
	gis = MemNew((size_t)n * sizeof(Int4));
	for(curr = gi_list, i=0; curr != NULL; curr = curr->next)
		gis[i++] = curr->data.intvalue;
	ValNodeFree(gi_list);
	bb_sort(gis, (Int2)n);
	switch(blast_type)
	{
		case DO_BLAST_P:
			db = TYP_AA;
			break;
		case DO_BLAST_N:
			db = TYP_NT;
			break;
		case DO_BLAST_X:
			db = TYP_AA;
			break;
		case DO_T_BLAST_N:
			db = TYP_NT;
			break;
		default:
			db = TYP_NT;
			break;
	}
	sprintf(name, "*pblast%ld", seq_count);
	EntrezCreateNamedUidList(name, db, FLD_WORD, n, gis);

	sprintf(str, "%s [WORD] & \"%s\" [ORGN]", name, organism);
	gi_bsp = EntrezTLEvalXString(str, db, -1, NULL, NULL);
	MemFree(gis);

	if(gi_bsp != NULL && gi_bsp->totlen == 0)
	{	/*there is no record in the ByteStore */
		BSFree(gi_bsp);
		return NULL;
	}
	return gi_bsp;
}
	
	
static Boolean bin_search_gi(Int4 low, Int4 high, Int4 gi, ByteStorePtr gi_bsp)
{
	Int4 mid;
	Int4 value;

	if(low > high)
		return FALSE;
	mid = (low + high)/2;
	BSSeek(gi_bsp, mid *4, SEEK_SET);
	BSRead(gi_bsp, &value, 4);
	if(gi == value)
		return TRUE;
	if(gi < value)
		return bin_search_gi(low, mid-1, gi, gi_bsp);
	else
		return bin_search_gi(mid+1, high, gi, gi_bsp);
}
		
static Boolean is_gi_in_list(Int4 gi, ByteStorePtr gi_bsp)
{
	Int4 size;

	if(gi_bsp == NULL)
		return TRUE;
	size = gi_bsp->totlen;
	if(size == 0)
		return FALSE;
	if(size%4 != 0)
	{
		Message(MSG_ERROR, "Incorrect size for ByteStore");
		return FALSE;
	}

	return bin_search_gi(0, (size/4-1), gi, gi_bsp);
}

/***************************************************************
*
*	filter_organism(sap, bsp, organism)
*	filter the output of Seq-align by specific organisms
*
***************************************************************/
static void filter_organism(SeqAnnotPtr sap, BioseqPtr bsp, ByteStorePtr gi_bsp, CharPtr organism, Uint1 filter_org, Uint1 blast_type, Int4 seq_count)
{
	SeqAlignPtr align, prev, next;
	DenseDiagPtr ddp;
	StdSegPtr ssp;
	DenseSegPtr dsp;
	ByteStorePtr ngi_bsp;
	SeqIdPtr sip;
	Boolean keep;

	if(filter_org == 0 || (gi_bsp == NULL && organism == NULL))
		return;
	if(gi_bsp == NULL)
	{
		/*use the second strategy*/
		ngi_bsp = get_filtered_gi_list(sap, organism, bsp, blast_type, seq_count);
		if(ngi_bsp == NULL)
		{
			if(filter_org == 1)
			{
				align = sap->data;
				sap->data = NULL;
				while(align)
				{
					next = align->next;
					SeqAlignFree(align);
					align = next;
				}
			}
			return;
		}
	}
	else
		ngi_bsp = gi_bsp;
	

	align = sap->data;
	prev = NULL;
	while(align)
	{
		next = align->next;
		keep = TRUE;
		sip = NULL;
		switch(align->segtype)
		{
			case 1:
				ddp = align->segs;
				sip = ddp->id->next;
				break;
			case 2:
				dsp = align->segs;
				sip = dsp->ids->next;
				break;
			case 3:
				ssp = align->segs;
				sip = SeqLocId(ssp->loc->next);
				break;
			default:
				break;
		}
		if(sip != NULL && sip->choice == SEQID_GI)
		{
			if(!BioseqMatch(bsp, sip))
			{
				keep = is_gi_in_list(sip->data.intvalue, ngi_bsp);
				if(filter_org == 2)
					keep = 1 - keep;
			}
		}
		if(keep)
			prev = align;
		else
		{
			if(prev)
				prev->next = align->next;
			else
				sap->data = align->next;
			align->next = NULL;
			SeqAlignFree(align);
		}
		align = next;
	}
	if(gi_bsp == NULL)
		BSFree(ngi_bsp);
}


static Boolean load_blast_score_for_sim (SeqAlignPtr blast_sap, SeqAlignPtr sim_sap)
{
	Int4 b_score, s_score;
	DenseSegPtr dsp;
	DenseDiagPtr ddp;

	if(blast_sap->segtype == 1)
	{
		ddp = blast_sap->segs;
		dsp = sim_sap->segs;
		b_score = get_score_value(ddp->scores);
		s_score = get_score_value(dsp->scores);

		if(s_score == -1 || b_score > s_score)
		{
			if(dsp->scores != NULL)
				ScoreSetFree(dsp->scores);
			dsp->scores = ddp->scores;
			ddp->scores = NULL;
			return TRUE;
		}
	}
	return FALSE;
}

static Boolean replace_with_sim_sap(SeqAnnotPtr sap, SeqAnnotPtr sim_sap)
{
	SeqAlignPtr sim_align, align, prev, next;
	DenseSegPtr dsp;
	DenseDiagPtr ddp;
	Boolean found;

	if(sap == NULL || sim_sap == NULL)
		return FALSE;

	if(sim_sap->data == NULL)
	{
		SeqAnnotFree(sim_sap);
		return FALSE;
	}

	align = (SeqAlignPtr)(sap->data);
	prev = NULL;
	while(align)
	{
		next = align->next;
		found = FALSE;
		if(align->segtype == 1)
		{
			ddp = align->segs;
			for(sim_align = (SeqAlignPtr)(sim_sap->data); 
				sim_align != NULL; sim_align = sim_align->next)
			{
				dsp = sim_align->segs;
				if(SeqIdMatch(dsp->ids->next, ddp->id->next))
				{
					/* load_blast_score_for_sim (align, sim_align); */
					found = TRUE;
					break;
				}
			}
		}
		if(found)
		{
			if(prev == NULL)
				sap->data = next;
			else
				prev->next = next;
			SeqAlignFree(align);
		}
		else
			prev = align;
		align = next;
	}

	sap->data = link_align((SeqAlignPtr)(sap->data), (SeqAlignPtr)(sim_sap->data));
	sim_sap->data = NULL;
	SeqAnnotFree(sim_sap);
	return TRUE;
}

static Boolean write_text_output(SeqAnnotPtr annot, CharPtr name, Boolean io_type, Boolean is_html, Boolean hide_feature)
{
	FILE *fp;
	Char out_name[100];
	Uint1 f_order[FEATDEF_ANY], g_order[FEATDEF_ANY];
	Uint4 option;

	if(annot->data == NULL)
		return FALSE;

	if(is_html)
		sprintf(out_name, "%s.html", name);
	else
		sprintf(out_name, "%s.ali", name);
	
	option = 0;
	option |= TXALIGN_MASTER;
	option |= TXALIGN_MISMATCH;
	if(is_html)
		option |= TXALIGN_HTML;
	option |= TXALIGN_SHOW_RULER;
	
	if(io_type == FILE_IO) {
		fp = FileOpen(out_name, "w");
		if(fp == NULL) {
			Message(MSG_ERROR, "Fail to open file %s", out_name);
			return FALSE;
		}
	} else {
		fp = stdout;
	}
	if(hide_feature)
	{
		MemSet((Pointer)(g_order), 0, (size_t)(FEATDEF_ANY* sizeof(Uint1)));
		MemSet((Pointer)(f_order), 0, (size_t)(FEATDEF_ANY* sizeof(Uint1)));
		ShowTextAlignFromAnnot(annot, LINE, fp, f_order, g_order, option, NULL, NULL, NULL);
	}
	else
		ShowTextAlignFromAnnot(annot, LINE, fp,  NULL, NULL, option, NULL, NULL, NULL);
	if(io_type == FILE_IO) 
		FileClose(fp);
	return TRUE;
}

/*
*	clean up any empty Seq-annot
*
*/
static void clean_up_empty_sap(SeqAnnotPtr PNTR blast_sap)
{
	SeqAnnotPtr sap, next, prev;

	prev = NULL;
	sap = *blast_sap;

	while(sap)
	{
		next = sap->next;
		if(sap->data == NULL)
		{
			if(prev != NULL)
				prev->next = next;
			else
				*blast_sap = next;
			sap->next = NULL;
			SeqAnnotFree(sap);
		}
		else
			prev = sap;
		sap = next;
	}
}
				


/*
*	search for one record
*
*/

static Uint1 powblast_one_entry(BioseqPtr bsp, SeqLocPtr slp, 
				PBlastOptionPtr pbop, MonitorPtr mp)
{
	Uint1 blast_program;
	ValNodePtr ends_list;
	SeqLocPtr blastloc, cloc;
	SeqLocPtr aluloc, dustloc, filterloc;
	SeqEntryPtr sep;
	BioseqPtr dust_bsp;
	CharPtr bpam;	/*parameters for blast search */
	Char blast_prog[10];
	BLAST0ResponsePtr brp = NULL;
	Uint4 db_option;
	Int2 k, total;
	Uint4 i;
	Char db[21];
	CharPtr other_db;
	SeqAnnotPtr sap, n_sap, t_sap;
	SeqAnnotPtr blast_sap;
	SeqAnnotPtr sim_sap;
	Char buf[101];
	Uint1 blast_type;
	Char seq_name[101];
	Char out_name[201];
	Int4 overlap_len;
	Int4 max_len;
	Boolean use_monitor;
	Boolean skip;	/*skip the search type that mismatches with the molecule */
	

	if(bsp == NULL || pbop == NULL || slp == NULL)
		return POWBLAST_NONE;
	

	if(mp == NULL)
		use_monitor = FALSE;
	else
		use_monitor = TRUE;
	/* handle the repeats and low complexity masking */
	ends_list = NULL;
	filterloc = NULL;
	if(bsp->mol != Seq_mol_aa && (pbop->dust || pbop->alu_slp_list != NULL))
	{
		blastloc = break_blast_job(slp, bsp->id, MAX_BLASTN_LEN, 100);
		total = get_vnp_num(blastloc);
		k = 0;
		for(cloc = blastloc; cloc != NULL; cloc = cloc->next)
		{
			++k;
			if(pbop->alu_slp_list != NULL)
			{
				if(mp != NULL)
				{
					sprintf(buf, "Filter Repeats for %d out of %d Pieces", k, total);
					MonitorStrValue(mp, buf);
				}
				aluloc = filter_repeats(cloc, pbop->alu_slp_list, &ends_list);
				ValNodeLink(&filterloc, aluloc);
			}

			if(pbop->dust)
			{
				if(mp != NULL)
				{
					sprintf(buf, "Dust Low Complexity for %d out of %d Pieces", k, total);
					MonitorStrValue(mp, buf);
				}
				dustloc = SeqLocDust(cloc, -1, 2, -1, -1);
				ValNodeLink(&filterloc, dustloc);
			}
		}
		if(ends_list != NULL)
		{
			ends_list = CleanNewList(ends_list);
			SaveRepeats(bsp, ends_list);
			ValNodeFreeData(ends_list);
		}
		SeqLocSetFree(blastloc);
	}
	dust_bsp = make_dust_bsp(bsp, 0, (bsp->length-1), filterloc);
	/* prt_FASTA_file(dust_bsp, NULL, FILE_IO); */
	if(filterloc != NULL)
		SeqLocSetFree(filterloc);

	/*masking the existing repeat features in there is any*/
	if(pbop->alu_slp_list == NULL || bsp->mol == Seq_mol_aa)
	{
		sep = SeqEntryFind(SeqLocId(slp));
		mask_with_repeat_feature(bsp, sep, dust_bsp);
	}

	/********************************************************
	*
	*	start to run the blast search
	*
	********************************************************/

	overlap_len = OVERLAP_SPACE;
	blast_sap = NULL;
	for(blast_program = BLASTN_PROGRAM; blast_program <= TBLASTN_PROGRAM; blast_program *=2)
	{
		if(pbop->bpp->blast_program & blast_program)
		{
			skip = FALSE;
			switch(blast_program)
			{
			case BLASTN_PROGRAM:
				max_len = pbop->bpp->max_blast_n;
				bpam = pbop->bpp->n_param;
				StringCpy(blast_prog, "blastn");
				db_option = pbop->bpp->dna_db;
				other_db = pbop->bpp->other_dna;
				blast_type = 1;
				skip = (bsp->mol == Seq_mol_aa);
				
				break;
			case BLASTP_PROGRAM:
				max_len = pbop->bpp->max_blast_p;
				bpam = pbop->bpp->p_param;
				StringCpy(blast_prog, "blastp");
				db_option = pbop->bpp->prot_db;
				other_db = pbop->bpp->other_protein;
				blast_type = 2;
				skip = (bsp->mol != Seq_mol_aa);
				break;
			case BLASTX_PROGRAM:
				max_len = pbop->bpp->max_blast_x;
				bpam = pbop->bpp->x_param;
				StringCpy(blast_prog, "blastx");
				db_option = pbop->bpp->prot_db;
				other_db = pbop->bpp->other_protein;
				blast_type = 3;
				skip = (bsp->mol == Seq_mol_aa);
				break;
			case TBLASTN_PROGRAM:
				max_len = pbop->bpp->max_blast_t;
				bpam = pbop->bpp->t_param;
				StringCpy(blast_prog, "tblastn");
				db_option = pbop->bpp->dna_db;
				other_db = pbop->bpp->other_dna;
				blast_type = 4;
				skip = (bsp->mol != Seq_mol_aa);

				break;
			default:
				skip = TRUE;
				break;
			}

			if(skip == FALSE)
			{
				blastloc = break_blast_job(slp, dust_bsp->id, max_len, overlap_len);
				total = get_vnp_num(blastloc);
				sap = NULL;
				for(cloc = blastloc, k = 1; cloc != NULL; cloc = cloc->next, ++k)
				{
					n_sap = NULL;
					for(i = (other_db != NULL) ? 0 : BLAST_NR; i <= BLAST_THC;)
					{
						if(i & db_option)
						{
							switch(i)
							{
							case 0:
								StringCpy(db, other_db);
								break;
							case BLAST_NR:
								StringCpy(db, "nr");
								break;
							case BLAST_EST:
								StringCpy(db, "est");
								break;
							case BLAST_MONTH:
								StringCpy(db, "month");
								break;
							case BLAST_STS:
								StringCpy(db, "sts");
								break;
							case BLAST_HTGS:
								StringCpy(db, "htgs");
								break;
							case BLAST_VECTOR:
								StringCpy(db, "vector");
								break;
							case BLAST_MITO:
								StringCpy(db, "mito");
								break;
							case BLAST_KABAT:
								StringCpy(db, "kabat");
								break;
							case BLAST_SWISSPROT:
								StringCpy(db, "swissprot");
								break;
							case BLAST_PDB:
								StringCpy(db, "pdb");
								break;
							case BLAST_ALU:
								StringCpy(db, "alu");
								break;
							case BLAST_EPD:
								StringCpy(db, "epd");
								break;
							case BLAST_YEAST:
								StringCpy(db, "yeast");
								break;
							case BLAST_GSS:
								StringCpy(db, "gss");
								break;
							case BLAST_THC:
								StringCpy(db, "thc");
								break;
							default:
								break;
							}

							if(mp != NULL)
							{
								sprintf(buf, "%s  %s for %ld of %ld Pieces", blast_prog, db, k, total); 
								MonitorStrValue(mp, buf);
							}
							brp = NULL;
							t_sap = BlastSeqLoc2(cloc, blast_prog, db, bpam, &brp, NULL, use_monitor);
							if(print_blast_error(t_sap, brp, use_monitor))
							{
								if(sap != NULL)	/*sap for one program */
									SeqAnnotFree(sap);
								if(n_sap != NULL)	/*sap for one piece */
									SeqAnnotFree(n_sap);
								if(blast_sap != NULL)	/*sap for multiple program */
									SeqAnnotFree(blast_sap);
								return POWBLAST_FATAL;
							}
							else if(t_sap != NULL)
							{
								link_two_annot(&n_sap, t_sap);
								if(mp != NULL)
								{
									sprintf(buf, "Find Hits in %s for %d of %d Pieces", db, k, total); 
									MonitorStrValue(mp, buf);
								}
							}
							if(brp != NULL)
								BLAST0ResponseFree(brp);
	
						}
						if(i == 0)
							++i;
						else i*=2;
					}
					if(n_sap != NULL)
					{
						restore_blast_id(n_sap, SeqLocId(slp));
						SortAlignByLocation(n_sap);
						clean_all_internal_repeats((SeqAlignPtr PNTR )(&(n_sap->data)));
						merge_blast_result(&sap, n_sap);
						clean_all_internal_repeats((SeqAlignPtr PNTR )(&(sap->data)));
						clean_empty_seqalign((SeqAlignPtr PNTR )(&(sap->data)));
					}
				}
				SeqLocSetFree(blastloc);
	
				if(sap != NULL)
				{
					AddAlignInfoToSeqAnnot(sap, blast_type);
					link_annot_to_end(&blast_sap, sap);
				}
			}
		}
	}
	BioseqFree(dust_bsp);

	/*Postprocess blast search */
	for(sap = blast_sap; sap != NULL; sap = sap->next)
	{
		/*filter self */
		if(pbop->filter_self)
			filter_blast_query(sap, bsp);

		blast_type = get_blast_type (sap);

		/* organism filtering */
		if(pbop->filter_org && pbop->organism != NULL)
		{
			if(mp != NULL)
			{
				sprintf(buf, "Filter Hits for Organism %s", pbop->organism);
				MonitorStrValue(mp, buf);
			}
				
			++(pbop->seq_count);
			if(blast_type == 1 || blast_type == 4)
				filter_organism(sap, bsp, pbop->na_gi_bsp, 
					pbop->organism, pbop->filter_org, blast_type, pbop->seq_count);
			else
				filter_organism(sap, bsp, pbop->aa_gi_bsp, 
					pbop->organism, pbop->filter_org, blast_type, pbop->seq_count);
		}

		/*run gapped alignment */
		if(sap != NULL && sap->data != NULL && pbop->gap_alignment && (blast_type == 1 || blast_type == 2))
		{
			if(mp != NULL)
				MonitorStrValue(mp, "Running Gapped Alignment. Takes Time");
			sim_sap = NULL;
			if(blast_type == 1)
				sim_sap = sim_for_blast(sap,  SeqLocId(slp), pbop->gap_alignment);
			else
				sim_sap = sim_for_blast(sap, SeqLocId(slp), RUN_SIM_1);
			if(sim_sap != NULL)
			{
				replace_with_sim_sap(sap, sim_sap);
			}
		}

		if(sap != NULL && sap->data != NULL && pbop->map_align_to_feature)
			load_align_to_interval(sap, MIN_GAP_LEN, bsp);
	}

	clean_up_empty_sap(&blast_sap);

	seqid_name(bsp->id, seq_name, FALSE, TRUE);
	if(pbop->output_path)
		sprintf(out_name, "%s%s", pbop->output_path, seq_name);
	else
		sprintf(out_name, "%s", seq_name);
	if(blast_sap == NULL)
	{
		if(use_monitor)
			ErrPostEx(SEV_WARNING, 0, 0, "No Hit was found for %s", seq_name);
		if(pbop->errfp != NULL)
			fprintf(pbop->errfp, "No Hit was found for %s\n", seq_name);
		if(pbop->output_format & OUTPUT_SEQENTRY)
		{
			write_asn_output(out_name, bsp, NULL, 2, FILE_IO);
			ValNodeCopyStr(&(pbop->output_file), OUTPUT_SEQENTRY, out_name);
		}
		return POWBLAST_NONE;
	}
	link_annot_to_end (&(bsp->annot), blast_sap);
	if(pbop->output_format & OUTPUT_TEXT)
	{
		write_text_output(blast_sap, out_name, FILE_IO, FALSE, pbop->hide_alignment_feature);
		ValNodeCopyStr(&(pbop->output_file), OUTPUT_TEXT, out_name);
	}

	if(pbop->output_format & OUTPUT_HTML)
	{
		write_text_output(blast_sap, out_name, FILE_IO, TRUE, pbop->hide_alignment_feature);
		ValNodeCopyStr(&(pbop->output_file), OUTPUT_TEXT, out_name);
	}

	if(pbop->output_format & OUTPUT_SEQALIGN)
	{
		save_SeqAnnot(blast_sap, out_name, FILE_IO);
		ValNodeCopyStr(&(pbop->output_file), OUTPUT_SEQALIGN, out_name);
	}
	if(pbop->output_format & OUTPUT_SEQENTRY)
	{
		write_asn_output(out_name, bsp, NULL, 2, FILE_IO);
		ValNodeCopyStr(&(pbop->output_file), OUTPUT_SEQENTRY, out_name);
	}
	return POWBLAST_HIT;
}


