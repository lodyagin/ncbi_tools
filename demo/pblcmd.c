#include <pblproc.h>
#include <netblap2.h>
#include <accentr.h>
#include <tofasta.h>
#include <accutils.h>
#include <lsqfetch.h>
#include <edutil.h>
/* #include <getfa.h> */



#define NUMARGS 17
			

Args myargs[NUMARGS] = { 
	{ "The file name for power blast job", NULL, NULL, NULL, FALSE, 'i', ARG_STRING, 0.0, 0, NULL},
	{ "Reset the options 0=No 1=Reset 2=Reset+Save 3=Modify 4=Modify+Save\n", "0", "0", "4", TRUE, 'c', ARG_INT, 0.0, 0, NULL}, 
	{ "The repeat FASTA library file for filtering", NULL, NULL, NULL, TRUE, 'l', ARG_STRING, 0.0, 0, NULL},
	{ "dust the sequence before blast", "TRUE", NULL, NULL, TRUE, 'd', ARG_BOOLEAN, 0.0, 0, NULL},
	{ "filter the blast output with the organism? 0=NO 1=Keep 2=Filter\n", "0", "0", "2", TRUE, 'f', ARG_INT, 0.0, 0, NULL},
	{ "the name for organism for filtering", NULL, NULL, NULL, TRUE, 'o', ARG_STRING, 0.0, 0, NULL},
	{ "compute gapped alignment 0=No 1=sim1 2=sim2 3=sim3", "0", "0", "3", TRUE, 's', ARG_INT, 0.0, 0, NULL},
	{ "export the results as 1=text(*.ali) 8=HTML(*.html) \n"
	"     2=Seq-align(*.sat) 4=Seq-entry(*.ent)\n", "0", "0", "15", TRUE, 'a', ARG_INT, 0.0, 0, NULL},
	{ "type of blast 0=default 1=blastn 2=blastp 4=blastx 8=tblastn\n", "0", "0", "15", TRUE, 'b', ARG_INT, 0.0, 0, NULL}, 
	{ "Search Nucleotide databases: 1=nr 2=est 4=sts 8=month 16=htgs 32=vector \n"
          "     64=mito 128=kabat 512=pDB epd=1024 yeast=2048 gss=4096 alu=8192\n", "1", NULL, NULL, TRUE, 'N', ARG_STRING, 0.0, 0, NULL},
	{ "Search Protein databases: 1=nr 8=month 128=kabat 256=swissprot 512=pdb \n"
	"      yeast=2048 alu=8192\n", "1", NULL, NULL, TRUE, 'A', ARG_STRING, 0.0, 0, NULL},
	{ "Parameters for BLASTN, use quote", NULL, NULL, NULL, TRUE, 'n', ARG_STRING, 0.0, 0, NULL},
	{ "Parameters for BLASTX, use quote", NULL, NULL, NULL, TRUE, 'x', ARG_STRING, 0.0, 0, NULL},
	{ "Parameters for BLASTP, use quote", NULL, NULL, NULL, TRUE, 'p', ARG_STRING, 0.0, 0, NULL},
	{ "Parameters for TBLASTN, use quote", NULL, NULL, NULL, TRUE, 't', ARG_STRING, 0.0, 0, NULL},
	{ "filter out the GenBank query itself", "FALSE", NULL, NULL, TRUE, 'q', ARG_BOOLEAN, 0.0, 0, NULL},
	{ "Enable the Monitor", "TRUE", NULL, NULL, TRUE, 'm', ARG_BOOLEAN, 0.0, 0, NULL}
	/* { "Map alignments to feature intervals", "FALSE", NULL, NULL, TRUE, 'm', ARG_BOOLEAN, 0.0, 0, NULL} */
};

static Boolean FileExists (CharPtr dirname, CharPtr subname, CharPtr filename)

{
  Char  path [PATH_MAX];

  StringNCpy (path, dirname, sizeof (path));
  FileBuildPath (path, subname, NULL);
  FileBuildPath (path, NULL, filename);
  return (Boolean) (FileLength (path) > 0);
}

static Boolean CheckAsnloadPath (CharPtr dirname, CharPtr subdir)

{
#ifdef ASNLOAD_NEEDED
  Char  fname [16];
  int   i;

  for (i = 60; i <= 69; ++i) {
    sprintf (fname, "asnmedli.l%02d", (int) i);
    if (FileExists (dirname, subdir, fname)) {
      return TRUE;
    }
  }
  return FALSE;
#else
  return TRUE;
#endif
}


static Boolean CheckDataPath (CharPtr dirname, CharPtr subdir)

{
  return (Boolean) (FileExists (dirname, subdir, "seqcode.val"));
}


static void SetTransientPath (CharPtr dirname, CharPtr subname, CharPtr file,
                              CharPtr section, CharPtr type)

{
  Char  path [PATH_MAX];

  StringNCpy (path, dirname, sizeof (path));
  FileBuildPath (path, subname, NULL);
  TransientSetAppParam (file, section, type, path);
}

static Boolean UseLocalAsnloadAndData (void)

{
  Boolean  asnFound;
  Boolean  dataFound;
  Char     path [PATH_MAX];
  CharPtr  ptr;

  ProgramPath (path, sizeof (path));
  ptr = StringRChr (path, DIRDELIMCHR);
  if (ptr != NULL) {
    *ptr = '\0';
  }
  asnFound = CheckAsnloadPath (path, "asnload");
  dataFound = CheckDataPath (path, "data");
  if (asnFound && dataFound) {
    SetTransientPath (path, "asnload", "NCBI", "NCBI", "ASNLOAD");
    SetTransientPath (path, "data", "NCBI", "NCBI", "DATA");
	return TRUE;
  }
  return FALSE;
}



static Boolean check_input_file(CharPtr f_name)
{
  Char     path [PATH_MAX];
  CharPtr  ptr;

  ProgramPath (path, sizeof (path));
  ptr = StringRChr (path, DIRDELIMCHR);
  if (ptr != NULL) {
    *ptr = '\0';
  }
  return FileExists (path, NULL, f_name);
}


	
Int2 Main(void)
{
	PBlastOptionPtr pbop;
	BlastParamPtr bpp;
	Int4 value;
	Int4 setting;
	Boolean reset;

	if(!GetArgs("power blast", NUMARGS, myargs))
		return 1;

	if(!UseLocalAsnloadAndData ())
	{
		if(!SeqEntryLoad())
			return 1;
	}

	if(!check_input_file(myargs[0].strvalue))
	{
		Message(MSG_ERROR, "Fail to open file %s", myargs[0].strvalue);
		return 1;
	}

  	if (! BlastInit("power blast", FALSE)) {
    		ErrPostEx(SEV_FATAL, 0, 0, "Unable to initialize BLAST service");
    		return (1);
  	}

	if(!EntrezBioseqFetchEnable ("powblast", TRUE))
	{
		Message(MSG_ERROR, "Fail to Connect to the Entrez Server");
		BlastFini();
		return 1;
	}

	BioseqFetchInit(TRUE);
	/* kludge so that it doesn't try to use brokered Entrez */
	TransientSetAppParam("NCBI", "NET_SERV", "SOME_BROKERED", "FALSE");

	pbop = MemNew(sizeof(PBlastOption));
	StringCpy(pbop->file_name, myargs[0].strvalue);
	pbop->seq_data = NULL;
	load_default_blast_param(pbop);
	bpp = pbop->bpp;

	setting = myargs[1].intvalue;
	if(setting != 0)
	{
		reset = (setting == 1 || setting == 2);
		if(myargs[2].strvalue)
			StringCpy(pbop->repeat_library, myargs[2].strvalue);
		else if(reset)
			pbop->repeat_library[0] = '\0'; 
		pbop->dust = (Boolean)(myargs[3].intvalue);
		if(reset || pbop->filter_org == 0)
			pbop->filter_org = myargs[4].intvalue;
		if(myargs[5].strvalue || reset)
		{
			if(pbop->organism != NULL)
				pbop->organism = MemFree(pbop->organism);
			if(myargs[5].strvalue)
				pbop->organism = StringSave(myargs[5].strvalue);
		}
		else if(reset && pbop->organism)
			pbop->organism = MemFree(pbop->organism);
		if(reset || pbop->gap_alignment == 0)
			pbop->gap_alignment = myargs[6].intvalue;
		if(myargs[7].intvalue != 1 || reset)
			pbop->output_format = myargs[7].intvalue;
		if(myargs[8].intvalue != 0 || reset)
			bpp->blast_program = myargs[8].intvalue;
		value = atol(myargs[9].strvalue);
		if(value == 0)
		{
			if(bpp->other_dna != NULL)
				MemFree(bpp->other_dna);
			bpp->other_dna = StringSave(myargs[9].strvalue);
		}
		else if(value != 1 || reset)
			bpp->dna_db = (Uint4)value;

		value = atol(myargs[10].strvalue);
		if(value == 0)
		{
			if(bpp->other_protein != NULL)
				MemFree(bpp->other_protein);
			bpp->other_protein = StringSave(myargs[10].strvalue);
		}
		else if(value != 1 || reset)
			bpp->prot_db = (Uint4)value;


		if(myargs[11].strvalue || reset)
		{
			if(bpp->n_param)
				bpp->n_param = MemFree(bpp->n_param);
			if(myargs[11].strvalue)
				bpp->n_param = StringSave(myargs[11].strvalue);
		}

		if(myargs[12].strvalue || reset)
		{
			if(bpp->x_param)
				bpp->x_param = MemFree(bpp->x_param);
			if(myargs[12].strvalue)
				bpp->x_param = StringSave(myargs[12].strvalue);
		}

		if(myargs[13].strvalue || reset)
		{
			if(bpp->p_param)
				bpp->p_param = MemFree(bpp->p_param);
			if(myargs[13].strvalue)
				bpp->p_param = StringSave(myargs[13].strvalue);
		}

		if(myargs[14].strvalue || reset)
		{
			if(bpp->t_param)
				bpp->t_param = MemFree(bpp->t_param);
			if(myargs[14].strvalue)
				bpp->t_param = StringSave(myargs[14].strvalue);
		}

		if(myargs[15].intvalue == 0|| reset)
			pbop->filter_self = (Boolean)(myargs[15].intvalue != 0);
		if(myargs[16].intvalue == 0 || reset)
			pbop->monitor = (Boolean)(myargs[16].intvalue != 0);

		if(myargs[1].intvalue ==2 || myargs[1].intvalue == 4)
			/*write the current session to the configuration file */
			load_pboption_to_configure(pbop);
	}

	RunPowerBlast(pbop);

	FreePBlastOption(pbop);
	BioseqFetchDisable();
	EntrezBioseqFetchDisable();
	BlastFini();

	return 0;
}
