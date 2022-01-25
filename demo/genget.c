/*****************************************************************************
*
*   seqget.c
*     entrez version
*
*   "Fasta style" SeqIds include a string indicating the class of SeqId,
*      vertical bar, then fields from the SeqId separated by vertical bar
*      If an (OPTIONAL) field is missing, the vertical bar must still be
*      there.
*
* local = lcl|integer or string
* gibbsq = bbs|integer
* gibbmt = bbm|integer
* giim = gim|integer
* genbank = gb|accession|locus
* embl = emb|accession|locus
* pir = pir|accession|name
* swissprot = sp|accession|name
* patent = pat|country|patent number (string)|seq number (integer)
* other = oth|accession|name|release
* general = gnl|database(string)|id (string or number)
* gi = gi|integer
* ddbj = dbj|accession|locus
* prf = prf|accession|name
* pdb = pdb|entry name (string)|chain id (char)
*
*****************************************************************************/
#include <accentr.h>
#include <sequtil.h>
#include <tofasta.h>
#include <asn2ff.h>
#include <asn2ffg.h>

#define NUMARGS 8
Args myargs[NUMARGS] = {
	{"GI id for single Bioseq to dump" ,"0","1","99999999",
				TRUE,'g',ARG_INT,0.0,0,NULL},
	{"Input filename ","stdin", NULL,NULL,TRUE,'i',ARG_FILE_IN, 0.0,0,NULL},
	{"Input asnfile in binary mode",
			"F",NULL,NULL,TRUE,'b',ARG_BOOLEAN,0.0,0,NULL},
	{"Output complete flat file?","T", NULL ,NULL ,TRUE,'c',
			ARG_BOOLEAN,0.0,0,NULL},
	{"Output short flat file?","T", NULL ,NULL ,TRUE,'s',
				ARG_BOOLEAN,0.0,0,NULL},
	{"Output FASTA file?","T", NULL ,NULL ,TRUE,'f',ARG_BOOLEAN,0.0,0,NULL},
	{"New output filename", NULL, NULL ,NULL ,TRUE,'o',ARG_STRING,0.0,0,NULL},
	{"Output compress file?","F", NULL ,NULL ,TRUE,'z',ARG_BOOLEAN,0.0,0,NULL},
	};

Int2 Main(void)
{
	Int2         retcode;  /* Default is -1 (genomes)     */
	Int4         gi = 0;
	SeqEntryPtr  sep;
	FILE *gbk=NULL; 
	FILE *gbs=NULL; 
	FILE *fsa=NULL;
	AsnIoPtr aip = NULL, aop = NULL;
	Boolean      is_network;
	Char tbuf[40];
	CharPtr outmode;
	Asn2ffJobPtr ajp;
	Uint2 entityID;
	CharPtr	infile, infile0, s, tarfile = NULL;
	static Char command[255];
	static CharPtr afile, bfile, cfile, dfile, efile;
	Int2 f;

	if ( !GetArgs("GenGet 1.0", NUMARGS, myargs) ) return 1;

	if ( !EntrezInit("GenGet", FALSE, &is_network) ) {
			ErrPostEx(SEV_FATAL,0,0, "Can't initialize Entrez");
			return 1;
	}
	gi = myargs[0].intvalue;

	if (gi != 0) {
		retcode = -1;
		sep = EntrezSeqEntryGet(gi, retcode);
		if (sep == NULL) {
			ErrPostEx(SEV_FATAL,0,0,
				"Could not retrieve entry for GI %ld", 	(long)gi);
			return 1;
		}
		
	} else {
		if ((aip = 
		AsnIoOpen (myargs[1].strvalue, myargs[2].intvalue?"rb":"r")) == NULL) {
			ErrPostEx(SEV_FATAL,0,0,"Could not open file %s", 	myargs[1].strvalue);
			return 1;
		}
		sep = SeqEntryAsnRead(aip, NULL);
		if (sep == NULL) {
			ErrPostEx(SEV_FATAL,0,0,"SeqEntryAsnRead failed");
			return 1;
		}
	}
	
	/* create output filenames */
	if (myargs[6].strvalue != NULL) {
		infile0 = (CharPtr) MemNew(StringLen(myargs[6].strvalue) + 1);
		infile = infile0;
		StringCpy(infile, myargs[6].strvalue);
	} else if (gi) {
		f = (Int2) log10((double) gi) + 1;
		infile0 = (CharPtr) MemNew(f + 2);
		infile = infile0;
		sprintf(infile, "%ld", gi);
		
	} else {
		infile0 = (CharPtr) MemNew(StringLen(myargs[1].strvalue) + 1);
		infile = infile0;
		StringCpy(infile, myargs[1].strvalue);
		if ((s = strrchr(infile, '/')) != NULL) {
			infile = s + 1;
		}
		if ((s = strchr(infile, '.')) != NULL) {
			*s = '\0';
		}
	}
	if (myargs[7].intvalue) {
		tarfile = (CharPtr) MemNew(StringLen(infile) + 5);
		StringCpy(tarfile, infile);
		StringCat(tarfile, ".tar");
		dfile = (CharPtr) MemNew(StringLen(infile) + 5);
		StringCpy(dfile, infile);
		StringCat(dfile, ".prt");
		efile = (CharPtr) MemNew(StringLen(infile) + 5);
		StringCpy(efile, infile);
		StringCat(efile, ".val");

	}
	if (myargs[3].intvalue) {
		afile = (CharPtr) MemNew(StringLen(infile) + 5);
		StringCpy(afile, infile);
		StringCat(afile, ".gbk");

		if ( (gbk = FileOpen (afile, "w")) == NULL) {
			ErrPostEx(SEV_ERROR,0,0, "Can't open %s", afile);
			exit (1);
		}
	}
	if (myargs[4].intvalue) {
		bfile = (CharPtr) MemNew(StringLen(infile) + 5);
		StringCpy(bfile, infile);
		StringCat(bfile, ".gbs");
		if ( (gbs = FileOpen (bfile, "w")) == NULL) {
			ErrPostEx(SEV_ERROR,0,0, "Can't open %s", bfile);
			exit (1);
		}
	}
	if (myargs[5].intvalue) {
		cfile = (CharPtr) MemNew(StringLen(infile) + 5);
		StringCpy(cfile, infile);
		StringCat(cfile, ".fsa");
		if ( (fsa = FileOpen (cfile, "w")) == NULL) {
			ErrPostEx(SEV_ERROR,0,0, "Can't open %s", cfile);
			exit (1);
		}
	}
	EntrezBioseqFetchEnable ("genget", FALSE); 
	ajp = (Asn2ffJobPtr) MemNew(sizeof(Asn2ffJob));
	ajp->show_seq = TRUE;
	ajp->show_gi = TRUE;
	ajp->error_msgs = FALSE;
	ajp->null_str = FALSE;
	ajp->non_strict = TRUE;
	ajp->format = GENBANK_FMT;
	ajp->mode = RELEASE_MODE;
	ajp->gb_style = FALSE;
	ajp->only_one = TRUE;
	ajp->ignore_top = TRUE;
	ajp->sep = sep;
	entityID = ObjMgrGetEntityIDForPointer(sep);
	if (entityID == 0)
		ErrPostStr(SEV_WARNING, 0, 0, "Couldn't get entityID");
				ajp->entityID = entityID;

	if (gbk != NULL) {
		ajp->only_one = TRUE;
		ajp->fp = gbk;
		asn2ff_print(ajp);
		FileClose(gbk);
	}
	if (gbs != NULL) {
		ajp->only_one = TRUE;
		ajp->genome_view = TRUE;
		ajp->fp = gbs;
		asn2ff_print(ajp);
		FileClose(gbs);
	}
	if (fsa != NULL) {
		SeqEntrysToFasta(sep, fsa, TRUE, 1); /* for genomes only */
		FileClose(fsa);
	}

	if (tarfile && gbk && gbs && fsa) {
		if ((aop = AsnIoOpen (dfile, "w")) == NULL) {
			ErrPostEx(SEV_FATAL,0,0,"Could not open file %s", dfile);
		} else {
			SeqEntryAsnWrite(sep, aop, NULL);
			AsnIoClose(aop);
		}
		if ((aop = AsnIoOpen (efile, "wb")) == NULL) {
			ErrPostEx(SEV_FATAL,0,0,"Could not open file %s", efile);
		} else {
			SeqEntryAsnWrite(sep, aop, NULL);
			AsnIoClose(aop);
		}
		sprintf(command, "tar cvf %s %s %s %s", tarfile, afile, bfile, cfile, dfile, efile);
		system(command);
		sprintf(command, "compress -v %s", tarfile);
		system(command);
	}
	MemFree(infile0);
	MemFree(afile);
	MemFree(bfile);
	MemFree(cfile);
	MemFree(dfile);
	MemFree(efile);
	MemFree(ajp);
	EntrezBioseqFetchDisable ();
	EntrezFini();
	if (aip) 
		AsnIoClose(aip);
	SeqEntryFree(sep);

	return 0;
}


