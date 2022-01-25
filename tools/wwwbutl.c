/* $Id: wwwbutl.c,v 6.1 2000/05/17 15:53:40 shavirin Exp $
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
* File Name:  $RCSfile: wwwbutl.c,v $
*
* Author:  Sergei Shavirin
*
* Initial Version Creation Date: 04/21/2000
*
* $Revision: 6.1 $
*
* File Description:
*         WWW BLAST/PSI/PHI utilities
*
* $Log: wwwbutl.c,v $
* Revision 6.1  2000/05/17 15:53:40  shavirin
* Initial revision.
*
*
* ==========================================================================
*/

#include <wwwblast.h>

static Int4 GlobalAlignNumber = 0; /* For PSI Blast printing ONLY !*/

void WWWBlastInfoFree(WWWBlastInfoPtr theInfo)
{
    Int4 i;
    WWWInfoFree(theInfo->info);
    BLASTOptionDelete(theInfo->options);
    MemFree(theInfo->database);
    MemFree(theInfo->program);

    for(i = 0; i < MAX_DB_NUM; i++) {
        MemFree(theInfo->blast_config->allow_db[i]);
    }
    MemFree(theInfo->blast_config);

    /* if(!theInfo->believe_query)
       fake_bsp = BlastDeleteFakeBioseq(fake_bsp); */

    MemFree(theInfo->blast_type);
    MemFree(theInfo->ConfigFile);
    MemFree(theInfo);

    return;
}
void WWWBlastErrMessage(BLASTErrCode error_code, CharPtr error_msg)
{
    CharPtr delim = "<BR>";

    if(error_code == BLASTNoError)
	return;

    printf("<HTML>\n");
    printf("<TITLE>BLAST Search Error</TITLE>\n"); 
    fflush(stdout);
    
    printf("<BODY BGCOLOR=\"#FFFFFF\" LINK=\"#0000FF\" "
           "VLINK=\"#660099\" ALINK=\"#660099\">\n");
    printf("<A HREF=\"blast_form.map\"> \r"
           "<IMG SRC=\"images/blast_results.gif\" "
           "BORDER=0 ISMAP>\r</A><P>\n");
    
    fprintf(stdout, "<FONT color=red><h3>");
    fprintf(stdout, "Error %d in submitting BLAST query", labs(error_code));

    fprintf(stdout, "</h3></FONT><HR>\n<b>");
    fprintf(stdout, "Short error description:");

    fprintf(stdout, "</b><BR><BR>\n");
    switch(error_code) {
        
    case BLASTEntrez:
        
        fprintf(stdout,
                "Your input sequence may not be found in Entrez %s"
                "or Entrez access interface currently unavailable. %s"
                "Please send message to blast_help@ncbi.nlm.nih.gov %s"
                "with description of your query", 
                delim, delim, delim);
        break;
        
    case BLASTFastaToSE:
        
        fprintf(stdout,
                "Your input sequence formatted incorrectly. %s"
                "Please read blast help if you have problems with formatting %s"
                "or send request to blast help account.", 
                delim, delim);
        break;
        
    case BLASTErrNoSequence:
        
        fprintf(stdout,
                "Input sequence for the BLAST search, probably missing. %s"
                "Please see the blast help for a description %s"
                "of the FASTA sequence format.", 
                delim, delim);
        break;
        
    case BLASTErrCombination:
        
        fprintf(stdout, 
                "The combination of database and program, that you provided in your %s"
                "message is invalid or not acceptable by BLAST search system. %s"
                "Please look at current possible combinations in BLAST help. ",
                delim, delim);
        break;
        
    case BLASTErrAccType:
        
        fprintf(stdout, 
                "You specified a protein (or nucleotide) sequence identifier, %s"
                "but a nucleotide (or protein) sequence is required for your search.",
                delim);
        
        break;
        
    case BLASTErrDatalib:
        
        fprintf(stdout, "No database was specified.  ");
        break;
        
    case BLASTErrNoQueue:
        fprintf(stdout, 
                "Unable to accept more BLAST jobs right now, %s"
                "Queue overloaded. Please try again later.", 
                delim);
        break;

    case BLASTOptionStr:
        
        if(error_msg != NULL) {
            fprintf(stdout, "%s", error_msg);
        }
        break;
        
    case BLASTMiscError:
    default:
        
        if(error_msg != NULL) {
            fprintf(stdout, "%s %s", error_msg, delim);
        } else {
            fprintf(stdout, 
                    "There were some internal software problems while processing %s"
                    "your request. Please contact blast support with a full %s"
                    "description of your query to BLAST as soon as possible.",
                    delim, delim); 
        }
        break;
    }
    
    fprintf(stdout, "\n<HR>\n");

    printf("</BODY>\n");
    printf("</HTML>\n");
    fflush(stdout);
    return;
}

Boolean BLAST_Time(CharPtr string, Int4 len, time_t seconds)
{
    CharPtr chptr;
    
    if(string == NULL || len < 25)
	return FALSE;
    
    if(!seconds) {
	seconds = GetSecs();
    }
    
    if((chptr = ctime(&seconds)) != NULL) 
        StringCpy(string, chptr);
    
    string[24] = NULLB;
    return TRUE;
}

WWWBlastInfoPtr WWWBlastReadArgs(CharPtr type)
{
    WWWBlastInfoPtr theInfo;
    WWWErrorCode error = WWWErrOk; 
    Int4 index;
    CharPtr blast_type, hostname, chptr;
    Char tmp[256];
    FILE *log_fd;

    theInfo = MemNew(sizeof(WWWBlastInfo));
    
    if((error = WWWGetArgs(&theInfo->info)) != WWWErrOk) {
        WWWBlastErrMessage(BLASTMiscError, NULL);
        return NULL;    
    }

    if((chptr = WWWGetQuery(theInfo->info)) == NULL || *chptr == NULLB) {
        fprintf(stdout, "<META HTTP-EQUIV=\"Refresh\" "
                "CONTENT=\"2; URL=%s.html\">", type ? type : "blast");
        return NULL;
    }
    
#ifdef PRINT_ALL_INPUT  /* Printing out all coming data for debugging */
    for(index= 0; index < WWWGetNumEntries(theInfo->info); index ++) {
        printf("%s : %s<BR>", 
               WWWGetNameByIndex(theInfo->info, index), 
               WWWGetValueByIndex(theInfo->info, index));
    }
#endif
    
    if(getenv("DEBUG_COMMAND_LINE") != NULL) {
        FILE *fd;
        fd = FileOpen("/tmp/__web.in", "w");
        fprintf(fd, "%s", ((WWWInfoDataPtr)theInfo->info)->query);
        FileClose(fd);
    }

    /* Root path for PSI/PHI Blast images */ 
    theInfo->www_root_path = getenv("WWW_ROOT_PATH"); /* May be NULL */
    
    if ( !ErrSetLogfile ("/dev/null", ELOG_APPEND) ) {
        fprintf(stdout, "Cannot set logfile exiting....\n");
        return FALSE;
    } else {
        ErrSetOpts (ERR_CONTINUE, ERR_LOG_ON);
    }
    
    /* Config file with program/database relationsship */
    
    blast_type = WWWGetValueByName(theInfo->info, "WWW_BLAST_TYPE");
    
    if(blast_type == NULL || *blast_type == NULLB) {
        theInfo->blast_type = StringSave(type ? type : "blast");
        sprintf(tmp, "%s.rc", theInfo->blast_type); 
        theInfo->ConfigFile = StringSave(tmp);
    } else {
        sprintf(tmp, "%s.rc", blast_type);
        theInfo->blast_type = StringSave(blast_type);
        theInfo->ConfigFile = StringSave(tmp);
    }
    
    sprintf(tmp, "%s.log", blast_type == NULL? "wwwblast" : blast_type);

    log_fd = FileOpen(tmp, "a"); 

    if(log_fd == NULL) /* If log_fd == NULL - no problem */
        return theInfo;
    
    BLAST_Time(tmp, sizeof(tmp), 0);

    if((hostname = getenv("PROXIED_IP")) == NULL)
        hostname = WWWGetAddress(theInfo->info);
    
    fprintf(log_fd, "\n%d|%s|%s|%s",
            getpid(), tmp, hostname == NULL? "host_not_set" : hostname,
            WWWGetAgent(theInfo->info));
    
    FileClose(log_fd);
    
    return theInfo;
}
static Int4 GetLetterIndex(CharPtr letters, Char ch)
{
    Int4 index;

    for(index = 0; letters[index] != NULLB; index++) {
	if (letters[index] == ch) {
	    return index;
	}
    }
    return -1;
}
static Boolean ParseInputString(CharPtr string, 
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

    max_par_num = strlen(letters);

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

	if((index = GetLetterIndex(letters, *chptr)) < 0) {
	    sprintf(message, "Character \'%c\' is not a valid option", 
		    *chptr);
	    *ErrorMessage = StringSave(message);
	    return FALSE;
	}

	if(*chptr == NULLB)   /* Check for NULLB */
	    break;

	while(!IS_WHITESP(*chptr)) { /* Finding first space */
	    if(*chptr == NULLB)   /* Check for NULLB */
		break;
	    chptr++;
	}

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
/* Set of functions to handle BLAST custom configuration file */
static BLASTConfigPtr BLASTConfigNew(void)
{
    BLASTConfigPtr config;
    
    if((config = (BLASTConfigPtr) MemNew(sizeof(BLASTConfig))) == NULL)
	return NULL;
    
    config->run_max = DEFAULT_RUN_MAX;
    config->queue_max = DEFAULT_QUEUE_MAX;
    config->num_cpu = NUM_CPU_TO_USE;
    MemSet(config->allow_db, 0, sizeof(CharPtr)*MAX_DB_NUM);

    return config;
}
static Int4 BLASTEatWs (FILE* fp)
{
    Int4 ch;

    while ((ch = fgetc (fp)) != EOF) {
	if (ch != ' ' && ch != '\t')
	    return ch;
    }
    return ch;
}

static void BLASTConfigGetWord(CharPtr word, CharPtr line) 
{
    Int4 x = 0, y = 0;
    
    for(x=0; line[x] && IS_WHITESP(line[x]); x++);

    while(TRUE) {
	if(!(word[y] = line[x]))
	    break;
	if(IS_WHITESP(line[x]))
	    if((!x) || (line[x-1] != '\\'))
		break;
	if(line[x] != '\\') ++y;
	++x;
    }
    word[y] = '\0';

    while(line[x] && IS_WHITESP(line[x])) ++x;

    for(y=0;(line[y] = line[x]);++x,++y);
}
static Int4 BLASTConfigGetLine (CharPtr s, Int4 n, FILE* fp)
{
    int   len = 0, ch;

    ch = BLASTEatWs(fp);

    while (TRUE) {
	if (ch == EOF || ch == '\n' || (len == n-1)) {
	    if (len && s[len - 1] == ' ') s[len - 1] = '\0';
	    else s[len] = '\0';
	    return feof(fp) ? 1 : 0;
	}
	s[len++] = ch;
	ch = fgetc (fp);

	if (ch == '\t' || ch == ' ') {
	    s[len++] = ch;
	    ch = BLASTEatWs(fp);
	}
    }
}

#define MAX_LINE_SIZE 2048
static BLASTConfigPtr BLASTReadConfigFile(CharPtr filename, CharPtr program)
{
    FILE *fd;
    BLASTConfigPtr config;
    Char line[MAX_LINE_SIZE], word[MAX_LINE_SIZE];
    Int4 value, i;
    
    if(filename == NULL)
	return NULL;
    
    if((config = BLASTConfigNew()) == NULL)
	return NULL;
    
    if((fd = FileOpen(filename, "r")) == NULL)
	return NULL;

    while(!(BLASTConfigGetLine(line, MAX_LINE_SIZE, fd))) {
	if((line[0] != '#') && (line[0] != '\0')) {
	    BLASTConfigGetWord(word, line);

	    if(!StringICmp(word, "RunMaxProcesses") && 
		    (value = atoi(line)) != 0) {
		config->run_max = value;
	    } else if(!StringICmp(word, "QueueMaxJobs") && 
		    (value = atoi(line)) != 0) {
		config->queue_max = value;
	    } else if(!StringICmp(word, "NumCpuToUse") && 
		    (value = atoi(line)) != 0) {
		config->num_cpu = value;
	    } else if(!StringICmp(word, "NiceValue") && 
		    (value = atoi(line)) != 0) {
		config->niceval = value;
	    } else if(!StringICmp(word, program)) {
		for(i = 0 ; line[0] != NULLB && i < MAX_DB_NUM; i++) {
		    BLASTConfigGetWord(word, line);
		    config->allow_db[i] = StringSave(word);
		}
	    }
	}
    }

    FileClose(fd);
    return config;
}

static Boolean ValidateCombinationsEx(WWWBlastInfoPtr theInfo, 
                                      CharPtr database)
{
    Int4 i;
    
    for(i = 0; theInfo->blast_config->allow_db[i] != NULL; i++) {
	if(!StringICmp(database, theInfo->blast_config->allow_db[i]))
	    return TRUE;
    }
    return FALSE;
}

/* This will work if search require to choose few databases */
static Boolean WWWParseDatabases(WWWBlastInfoPtr theInfo)
{
    Int4  count, index;
    Boolean done, datalib_found;
    Char buffer[4096], buffer1[4096]; /* is 4096 always long enough? XXX */
    CharPtr ptr, chptr;

    count = WWWGetNumEntries(theInfo->info);
    datalib_found = FALSE;
    ptr = buffer;
    
    for (index=0; index<count; index++) {
	chptr = WWWGetNameByIndex(theInfo->info, index);
	if (StringCmp(chptr, "DATALIB") == 0) {
	    datalib_found = TRUE;
	    chptr = WWWGetValueByIndex(theInfo->info, index);
	    done = FALSE;
            
            /* Parse string if multiple database names. */
	    while (done == FALSE) { 
		done = readdb_parse_db_names(&chptr, buffer1);
		if (ValidateCombinationsEx(theInfo, buffer1) == TRUE) {
                    
		    CharPtr prefix = WWWGetValueByName(theInfo->info, "DB_DIR_PREFIX");
		    Char tmpbuf[1024];
                    
		    if (prefix) {
			sprintf(tmpbuf, "%s%c%s", prefix, DIRDELIMCHR, buffer1);
		    } else {
			sprintf(tmpbuf, "%s", buffer1);
		    }
		    
		    StringCpy(ptr, tmpbuf);
		    ptr += StringLen(tmpbuf);
		    *ptr = ' '; ptr++;
		} else {
                    WWWBlastErrMessage(BLASTErrCombination, NULL);
                    return FALSE;
		}
	    }
	}
    }
    
    if (datalib_found) {
	ptr--;
	*ptr = NULLB;
	theInfo->database = StringSave(buffer);
    } else {
        WWWBlastErrMessage(BLASTErrDatalib, NULL);
	return FALSE;
    }
    
    /* Processing database aliases  */
    
    if(StringStr(theInfo->database, "E.coli") != NULL) {
	MemFree(theInfo->database);
	theInfo->database = StringSave("ecoli");
    }

    return TRUE;
}
#if defined(NCBI_CLIENT_SERVER) || defined (NCBI_ENTREZ_CLIENT)

static Int4 AccessionToGi (CharPtr string) 
{
    Char buffer[32];
    CharPtr chptr;
    Int2 version;
    Int4 gi, index;
    SeqIdPtr sip;
    TextSeqIdPtr tsip;
    PDBSeqIdPtr  psip;
    long tmplong;
    Boolean digit;

    for(chptr = string, digit = TRUE; *chptr != NULLB; chptr++) {
        if(!IS_DIGIT(*chptr)) {
            digit = FALSE;
            break;
        }
    }
        
    if(digit) {
        if((gi = atol(string)) > 0)
            return gi;
    }

    /* all letters in accesion should be upper */
    string = Nlm_StringUpper(string);
    
    gi = 0;

    if((sip = ValNodeNew (NULL)) == NULL)
        return -1;
    
    index = 0; version = 0;
    while (*string != '\0' && index < 16) {
        if (*string == '.')
            break;
        buffer[index] = *string;
        string++;
        index++;
    }

    buffer[index] = '\0';
    if (*string == '.' && *(string+1) != '\0') {
        sscanf((string+1), "%ld", &tmplong);
        version = (Int2) tmplong;
    }
    
    if((tsip = TextSeqIdNew ()) == NULL)
        return -1;
    
    tsip->accession = StringSave(buffer);
    tsip->version = version;
    
    /* GenBank, EMBL, and DDBJ. */
    sip->choice = SEQID_GENBANK;
    sip->data.ptrvalue = (Pointer) tsip;
    gi = ID1FindSeqId (sip);
    
    if (gi == 0) {
        /* SwissProt. */
        sip->choice = SEQID_SWISSPROT;
        gi = ID1FindSeqId (sip);
    } else {
        goto retpoint;
    }

    if (gi == 0) {
        /* PIR */
        sip->choice = SEQID_PIR;
        gi = ID1FindSeqId (sip);
    } else {
        goto retpoint;
    }

    if (gi == 0) {
        /* PRF */
        sip->choice = SEQID_PRF;
        gi = ID1FindSeqId (sip);
    } else {
        goto retpoint;
    }

    if (gi == 0) {
        /* OTHER, probably 'ref' */
        sip->choice = SEQID_OTHER;
        gi = ID1FindSeqId (sip);
    }

    if(gi != 0)
        goto retpoint;

    /* OK. We failed to find gi using string as TextSeqId. Now trying
       last time - with PDBSeqIdPtr */

    if((psip = PDBSeqIdNew()) == NULL)
        return -1;
    
    sip->choice = SEQID_PDB;
    sip->data.ptrvalue = psip;
    
    psip->mol = StringSave(buffer);
    psip->chain = version;

    gi = ID1FindSeqId (sip);

    SeqIdFree(sip);

 retpoint:
    TextSeqIdFree(tsip);
    return gi;
}
#endif

Boolean WWWCreateSearchOptions(WWWBlastInfoPtr theInfo)
{
    CharPtr chptr, ptr, sequence, outptr;
    SeqEntryPtr sep;
    Boolean gapped_set;
    CharPtr opt_str = "GErqeWvbKLY";
    BLAST_OptionsBlkPtr options;
    Char tmp[128];
    Int4 value;

    /* PROGRAM */
    
    if((chptr = WWWGetValueByName(theInfo->info, "PROGRAM")) != NULL) {
	theInfo->program = StringSave(chptr);
    } else {
        WWWBlastErrMessage(BLASTErrProgram, NULL);
	return FALSE;
    }

    /* Configuration file set program/database relations */

    if((theInfo->blast_config = 
        BLASTReadConfigFile(theInfo->ConfigFile, theInfo->program)) == NULL) {
        WWWBlastErrMessage(BLASTConfigFile, NULL);
	return FALSE;
    }

    /* DATALIB */
    if(!WWWParseDatabases(theInfo))
        return FALSE;
    
    /* SEQUENCE or SEQFILE */
    
    if((sequence = WWWGetValueByName(theInfo->info, "SEQUENCE")) == NULL ||
       sequence[0] == NULLB) {
        if((sequence = WWWGetValueByName(theInfo->info, "SEQFILE")) == NULL ||
           sequence[0] == NULLB) {
            WWWBlastErrMessage(BLASTErrNoSequence, NULL);
            return FALSE;
        }
    }
    
    theInfo->align_type = BlastGetTypes(theInfo->program, &theInfo->query_is_na, &theInfo->db_is_na);
#if defined(NCBI_CLIENT_SERVER) || defined (NCBI_ENTREZ_CLIENT)
    
    if((chptr = WWWGetValueByName(theInfo->info, "INPUT_TYPE")) != NULL && 
       !StringNICmp(chptr, "Accession", 9)) {

        Int4 gi, number, title_length, id_length;
        CharPtr accession, new_defline;
        BioseqPtr bsp_tmp;
        SeqIdPtr sip;
        ObjectIdPtr  oid;
        SeqPortPtr spp;
        Int2 retval, buf_length=512;
        Uint1 buf[512];
        Char tmp[255];

	/* This is request by Accession/GI - asking ID */
        
        if (!ID1BioseqFetchEnable("web-blasts", TRUE)) {
            WWWBlastErrMessage(BLASTEntrez, NULL);
	    return FALSE;
        }
        
	chptr = sequence;
        
	/* Strip off blanks at beginning. */
        while (IS_WHITESP(*chptr) && *chptr != NULLB)
            chptr++;
        outptr = chptr; /* Beginning of the valuable info */
        
	/* Strip off non-alphanumerics, except for '_' (used in SP) and '.' (soon to be used by  the collaboration. */
        while (IS_ALPHANUM(*chptr) || *chptr == '_' || *chptr == '.')
            chptr++;
        *chptr = NULLB; 
        
	/* accession = sequence; */
        accession = outptr;

	sip = NULL;
	gi = AccessionToGi(accession);
        
	if (gi > 0) {
	    ValNodeAddInt(&sip, SEQID_GI, gi);
	} else {
            WWWBlastErrMessage(BLASTEntrez, NULL);
	    return FALSE;
	}	

        /* If is is not found - it is not found - ID1 is down? */
        
	if((bsp_tmp = BioseqLockById(sip)) == NULL) {
            WWWBlastErrMessage(BLASTEntrez, NULL);
	    return FALSE;
        }

	if (ISA_na(bsp_tmp->mol) != theInfo->query_is_na) {
            WWWBlastErrMessage(BLASTErrAccType, NULL);
            return FALSE;
	}
        
	theInfo->query_bsp = BioseqNew();
	theInfo->query_bsp->length = bsp_tmp->length;
	theInfo->query_bsp->mol = bsp_tmp->mol;
	theInfo->query_bsp->repr = Seq_repr_raw;
	theInfo->query_bsp->seq_data = BSNew(theInfo->query_bsp->length);

	if (ISA_na(theInfo->query_bsp->mol)) {
            spp = SeqPortNew(bsp_tmp, 0, -1, Seq_strand_plus, 
                             Seq_code_iupacna);
            theInfo->query_bsp->seq_data_type = Seq_code_iupacna;
	} else {
            spp = SeqPortNew(bsp_tmp, 0, -1, Seq_strand_unknown, 
                             Seq_code_ncbieaa);
            theInfo->query_bsp->seq_data_type = Seq_code_ncbieaa;
	}

	SeqPortSet_do_virtual(spp, TRUE);
	number = 0;
	while (number < theInfo->query_bsp->length) {
            retval = SeqPortRead(spp, buf, buf_length);
            if (retval <= 0)
                break;
            BSWrite(theInfo->query_bsp->seq_data, buf, retval);
            number += retval;
	}

        SeqPortFree(spp);
        
	title_length = StringLen(BioseqGetTitle(bsp_tmp));
	SeqIdWrite(bsp_tmp->id, tmp, PRINTID_FASTA_LONG, 255);
	id_length = StringLen(tmp);
	title_length += id_length;
	title_length +=3;
	new_defline = (CharPtr) MemNew(title_length*sizeof(Char));
	StringCpy(new_defline, tmp);
	*(new_defline+id_length) = ' ';
	StringCpy(new_defline+id_length+1, BioseqGetTitle(bsp_tmp)); 
	*(new_defline+title_length-1) = NULLB;
	theInfo->query_bsp->descr = ValNodeAddStr(NULL, Seq_descr_title, 
                                                  new_defline);
	theInfo->query_bsp->id = ValNodeNew(NULL);
	oid = ObjectIdNew();
	oid->str = StringSave("blast_tmp");
	theInfo->query_bsp->id->choice = SEQID_LOCAL;
	theInfo->query_bsp->id->data.ptrvalue = (Pointer) oid;
        
	SeqMgrDeleteFromBioseqIndex(bsp_tmp);
        
	BioseqUnlock(bsp_tmp);
        
	BioseqPack(theInfo->query_bsp);
        ID1BioseqFetchDisable();
    }
#endif
    
    /* Creating Bioseq */
        
    if(theInfo->query_bsp == NULL) {
        if((sep = FastaToSeqBuffEx(sequence, &outptr, theInfo->query_is_na, 
                                   NULL, FALSE)) == NULL) {
            WWWBlastErrMessage(BLASTFastaToSE, NULL);
            return FALSE;
        }

        theInfo->query_bsp = (BioseqPtr) sep->data.ptrvalue;
    }

    /* The last check of Bioseq - if length of sequence too small ? */
    
    if(theInfo->query_bsp == NULL || 
       theInfo->query_bsp->length <= 0) {
        WWWBlastErrMessage(BLASTFastaToSE, NULL);
        return FALSE;
    }
    
    /* This will prevent from incorrect formating in case when input
       sequence has valig SeqId, but in fact this SeqId do not correspond
       to the real sequence  - XXX */

    if(!theInfo->believe_query)
        theInfo->fake_bsp = BlastMakeFakeBioseq(theInfo->query_bsp, NULL);
    else
        theInfo->fake_bsp = theInfo->query_bsp;
    
    /* OVERVIEW */
    
    if (WWWGetValueByName(theInfo->info, "OVERVIEW") != NULL)
        theInfo->show_overview = TRUE;

    /* UNGAPPED_ALIGNMENT */
    gapped_set = TRUE;
    if(WWWGetValueByName(theInfo->info, "UNGAPPED_ALIGNMENT") != NULL)
	gapped_set = FALSE;

    if((options = BLASTOptionNew(theInfo->program, gapped_set)) == NULL) {
        WWWBlastErrMessage(BLASTErrOptions, NULL);
	return FALSE; 
    }
    
    theInfo->options = options;

    /* Set default values fot matrix and gap parameters */
    BLASTOptionSetGapParams (options, NULL, 0, 0);
    
    /* Read MAT_PARAM if set */
    
    if(StringICmp("blastn", theInfo->program) && 
       (chptr = WWWGetValueByName(theInfo->info, "MAT_PARAM")) != NULL) {
        Char    matrixname[64];
        Int4    opencost, extendedcost;
        /* Get matrix name and gap costs */
        if (chptr[1] != '-' || chptr[2] != '-') {
            sscanf(chptr, "%s\t %ld\t %ld", matrixname, &opencost, 
                   &extendedcost);
            /* set the parameters */
            options->gap_open  = opencost;
            options->gap_extend  = extendedcost;
            if (options->matrix)
                MemFree(options->matrix);
            options->matrix = StringSave(matrixname);
        }
    } 

    if((chptr = WWWGetValueByName(theInfo->info, "GAP_OPEN")) != NULL &&
	    StringStr(chptr, "default") == NULL) {
	options->gap_open = atoi(chptr);
    }

    if((chptr = WWWGetValueByName(theInfo->info, "GAP_EXTEND")) != NULL &&
	    StringStr(chptr, "default") == NULL) {
	options->gap_extend = atoi(chptr);
    }

    if((chptr = WWWGetValueByName(theInfo->info, "GAP_VALUES")) != NULL &&
	    StringStr(chptr, "default") == NULL) {
	sscanf(chptr, "%d,%d", &options->gap_open, &options->gap_extend);
    }

    if((chptr = WWWGetValueByName(theInfo->info, "X_DROPOFF")) != NULL &&
	    StringStr(chptr, "default") == NULL) {
	options->gap_x_dropoff = atoi(chptr);
    }

    if (!StringICmp(theInfo->program, "blastn")) {
	options->penalty = -3;
	options->reward = 1;
    }

    if((chptr = WWWGetValueByName(theInfo->info, "GAP_SIZE")) != NULL &&
	    StringStr(chptr, "default") == NULL) {
	options->gap_size = atoi(chptr);
    }

    if((chptr = WWWGetValueByName(theInfo->info, 
	    "WINDOW_SIZE")) != NULL &&
	    StringStr(chptr, "default") == NULL) {
	options->window_size = atoi(chptr);
    }

    /* For BLASTX we could set genetic code */

    if (!StringICmp(theInfo->program, "blastx")) {
        BioSourcePtr source;

	options->genetic_code = 1;

	if((chptr = WWWGetValueByName(theInfo->info, 
		"GENETIC_CODE")) != NULL &&
		(StringStr(chptr, "default") == NULL)) {
	    chptr = StringChr(chptr, '(');
	    sscanf(chptr, "(%d", &value);
	    if(value != 0) {
		options->genetic_code = value; 

		source = BioSourceNew();
		source->org = OrgRefNew();
		source->org->orgname = OrgNameNew();
		source->org->orgname->gcode = options->genetic_code;
		ValNodeAddPointer(&theInfo->query_bsp->descr, 
			Seq_descr_source, source);
	    }
	}
    }

    /* For TBLASTN and TBLASTX we could set DB_GENETIC_CODE */

    if (!StringICmp(theInfo->program, "tblastn") || 
        !(StringICmp(theInfo->program, "tblastx"))) {
	options->db_genetic_code = 1;
	if((chptr = WWWGetValueByName(theInfo->info, 
                                      "DB_GENETIC_CODE")) != NULL &&
           (StringStr(chptr, "default") == NULL)) {
	    chptr = StringChr(chptr, '(');
	    sscanf(chptr, "(%d", &value);
	    if(value != 0) {
		options->genetic_code = value; 
	    }      
	}
    }

    if((chptr = WWWGetValueByName(theInfo->info, 
	    "THRESHOLD_1")) != NULL &&
	    (StringStr(chptr, "default") == NULL)) {
	options->threshold_first = atoi(chptr);
    } 

    if((chptr = WWWGetValueByName(theInfo->info, 
	    "THRESHOLD_2")) != NULL &&
	    (StringStr(chptr, "default") == NULL)) {
	options->threshold_second = atoi(chptr);
    }

    if((chptr = WWWGetValueByName(theInfo->info, 
	    "REQUIRED_START")) != NULL &&
	    StringStr(chptr, "default") != NULL) {
	options->required_start = atoi(chptr);
    }

    if((chptr = WWWGetValueByName(theInfo->info, 
	    "REQUIRED_END")) != NULL &&
	    StringStr(chptr, "default") != NULL) {
	options->required_end = atoi(chptr);
    }

    if((chptr = WWWGetValueByName(theInfo->info, 
	    "DROPOFF_1")) != NULL &&
	    (StringStr(chptr, "default") == NULL)) {
	options->dropoff_1st_pass = atof(chptr);
    }

    if((chptr = WWWGetValueByName(theInfo->info, 
	    "CUTOFF")) != NULL &&
	    (StringStr(chptr, "default") == NULL)) {
	options->cutoff_s = atof(chptr);
    }


    if((chptr = WWWGetValueByName(theInfo->info, 
	    "DROPOFF_2")) != NULL &&
	    (StringStr(chptr, "default") == NULL)) {
	options->dropoff_2nd_pass = atof(chptr);
    }

    /* MATRIX: */

    if((chptr = WWWGetValueByName(theInfo->info, 
	    "MATRIX")) != NULL &&
	    (StringStr(chptr, "default") == NULL)) {
	options->matrix = StringSave(chptr);
    }

    /* EXPECT */

    options->expect_value  = DEFAULT_EXPECT;

    if((chptr = WWWGetValueByName(theInfo->info, 
	    "EXPECT")) != NULL &&
	    StringStr(chptr, "default") == NULL) {
	options->expect_value = atof(chptr);
    }

    /* NUMBER OF BITS: */

    if((chptr = WWWGetValueByName(theInfo->info, 
	    "NUM_OF_BITS")) != NULL &&
	    (StringStr(chptr, "default") == NULL)) {
	options->number_of_bits = atof(chptr);
    }

    /* Number of CPU to use in BLAST Search: */

    if(theInfo->blast_config->num_cpu != 0)
        options->number_of_cpus = theInfo->blast_config->num_cpu;
    else
        options->number_of_cpus = NUM_CPU_TO_USE;

    /* CPU time limit. */
    
    options->cpu_limit = DEFAULT_CPU_LIMIT;

    /* FILTER: */


    options->filter = FILTER_NONE;

    if(WWWGetMethod(theInfo->info) == WWW_GET ||
	    (chptr = WWWGetValueByName(theInfo->info, "FSET")) != NULL) {
	if (!StringICmp(theInfo->program, "blastn")) {
	    options->filter = FILTER_DUST;
	} else {
	    options->filter = FILTER_SEG;
	}
    }
    
    {
        Char	tmpbuf[4096];
        /* Filter settings */
	int	i, num_entries = WWWGetNumEntries(theInfo->info);
        
	StringCpy(tmpbuf, "");

	for(i = 0; i < num_entries; i++) {
	    if((chptr = WWWGetNameByIndex(theInfo->info, i)) != NULL && 
		    !StringICmp(chptr, "FILTER")) {

		chptr = WWWGetValueByIndex(theInfo->info, i);
		/* add the filter */
		StringCat(tmpbuf, chptr);
		StringCat(tmpbuf, ";");
	    }
	}
	options->filter_string = StringSave(tmpbuf);
    }

    /* DESCRIPTIONS: */

    theInfo->number_of_descriptions = DEFAULT_DESCRIPTIONS;

    if((chptr = WWWGetValueByName(theInfo->info, 
                                  "DESCRIPTIONS")) != NULL && 
       StringStr(chptr, "default") == NULL) {
	theInfo->number_of_descriptions = atoi(chptr);
    }
    
    /* ALIGNMENTS */
    theInfo->number_of_alignments = DEFAULT_ALIGNMENTS;

    if((chptr = WWWGetValueByName(theInfo->info, "ALIGNMENTS")) != NULL &&
       StringStr(chptr, "default") == NULL) {
	theInfo->number_of_alignments = atoi(chptr);  
    }
    
    /* Now processing OTHER_ADVANCED_OPTIONS */
    
    if((chptr = WWWGetValueByName(theInfo->info, 
                                  "OTHER_ADVANCED")) != NULL) {
        Int4 index;
        CharPtr ErrorMessage = NULL;
        CharPtr PNTR values = NULL;

	if(!ParseInputString(chptr, opt_str, 
                             &values, &ErrorMessage)) {
            
            WWWBlastErrMessage(BLASTOptionStr, ErrorMessage);
	    return FALSE;
	}
        
	/* -G  gap open cost */
        
	index = GetLetterIndex(opt_str, 'G');
	if(values[index] != NULL) {
	    options->gap_open = atoi(values[index]);
	}
        
	/* -E gap extend cost */
        
	index = GetLetterIndex(opt_str, 'E');
	if(values[index] != NULL) {
	    options->gap_extend = atoi(values[index]);
	}
        
	/* -q penalty for nucleotide mismatch. */

	index = GetLetterIndex(opt_str, 'q');
	if(values[index] != NULL) {
	    options->penalty = atoi(values[index]);
	}

	/* -r reward for nucleotide match. */

	index = GetLetterIndex(opt_str, 'r');
	if(values[index] != NULL) {
	    options->reward = atoi(values[index]);
	}

	/* -e expect value. */

	index = GetLetterIndex(opt_str, 'e');
	if(values[index] != NULL) {
	    options->expect_value = atof(values[index]);
	}

	/* -W wordsize. */

	index = GetLetterIndex(opt_str, 'W');
	if(values[index] != NULL) {
	    options->wordsize = atoi(values[index]);
	}

	/* -v Number of descriptions to print. */

	index = GetLetterIndex(opt_str, 'v');
	if(values[index] != NULL) {
	    theInfo->number_of_descriptions = atoi(values[index]);
	}

	/* -b Number of alignments to show. */

	index = GetLetterIndex(opt_str, 'b');
	if(values[index] != NULL) {
	    theInfo->number_of_alignments = atoi(values[index]);
	}

        /* -K Number of best hits from a region to keep. */

        index = GetLetterIndex(opt_str, 'K');
        if(values[index] != NULL) {
            options->hsp_range_max = atoi(values[index]);
            if (options->hsp_range_max != 0)
                   options->perform_culling = TRUE;
        }

        /* -L Number of best hits from a region to keep. */

        index = GetLetterIndex(opt_str, 'L');
        if(values[index] != NULL) {
            options->block_width = atoi(values[index]);
        }

        /* -Y effective search space. */

        index = GetLetterIndex(opt_str, 'Y');
        if(values[index] != NULL) {
            options->searchsp_eff = atof(values[index]);
        }
        
        MemFree(values);
    }
        /*
          options->hsp_range_max = 100*options->hitlist_size;
          options->block_width = theInfo->bsp->length;
          */

    options->perform_culling = FALSE;


    /* Some values for PSI-Blast */
    value = 0;
    if((chptr = WWWGetValueByName(theInfo->info, "STEP_NUMBER")) != NULL)
	value = atoi(chptr);
    
    sprintf(tmp, "PSI_BEGIN%d", value-1);
    
    if((chptr = WWWGetValueByName(theInfo->info, tmp)) != NULL)
	options->required_start = atoi(chptr) - 1;
    
    sprintf(tmp, "PSI_END%d", value-1);
    if((chptr = WWWGetValueByName(theInfo->info, tmp)) != NULL)
	options->required_end = atoi(chptr) - 1;    
    
    if((chptr = WWWGetValueByName(theInfo->info, "E_THRESH")) != NULL)
	options->ethresh = atof(chptr);
    
    if((chptr = WWWGetValueByName(theInfo->info, "PHI_BLAST")) != NULL) {
	theInfo->is_phi_blast = TRUE;
        options->number_of_cpus = 1;
    }
    
    /* ------------------------ */

    /* HITLIST_SIZE: */
    
    theInfo->options->hitlist_size = MAX(theInfo->number_of_descriptions, 
                                         theInfo->number_of_alignments);
    
    /* ALIGNMENT VIEWS */
    
    if((chptr = WWWGetValueByName(theInfo->info, 
                                  "ALIGNMENT_VIEW")) != NULL &&
       StringStr(chptr, "default") == NULL) {
	theInfo->align_view = atoi(chptr);  
    }
    
    if (WWWGetValueByName(theInfo->info, "NCBI_GI") != NULL)
	theInfo->show_gi = TRUE;
    
    if (WWWGetValueByName(theInfo->info, "OVERVIEW") != NULL)
	theInfo->show_overview = TRUE;

    if (WWWGetValueByName(theInfo->info, "TAX_BLAST") != NULL)
        theInfo->show_tax_blast = TRUE;
    
    /* COLOR_SCHEMA */
    if((chptr = WWWGetValueByName(theInfo->info, "COLOR_SCHEMA")) != NULL &&
       StringStr(chptr, "No color schema") == NULL) {
	theInfo->color_schema = atoi(chptr);  
    }
    /* Now seting print and align options for BLAST output printing */

    theInfo->print_options = 0;
    theInfo->align_options = 0;
    theInfo->align_options += TXALIGN_COMPRESS;
    theInfo->align_options += TXALIGN_END_NUM;

    if (StringICmp("blastx", theInfo->program) == 0) {
        theInfo->align_options += TXALIGN_BLASTX_SPECIAL;
    }

    if (theInfo->show_gi) {
        theInfo->align_options += TXALIGN_SHOW_GI;
        theInfo->print_options += TXALIGN_SHOW_GI;
    }

    if (!gapped_set)
        theInfo->print_options += TXALIGN_SHOW_NO_OF_SEGS;
    
    if (theInfo->align_view) {
        theInfo->align_options += TXALIGN_MASTER;
        if (theInfo->align_view == 1 || theInfo->align_view == 3)
            theInfo->align_options += TXALIGN_MISMATCH;
        if (theInfo->align_view == 3 || theInfo->align_view == 4 || 
            theInfo->align_view == 6)
            theInfo->align_options += TXALIGN_FLAT_INS;
        if (theInfo->align_view == 5 || theInfo->align_view == 6)
            theInfo->align_options += TXALIGN_BLUNT_END;
    } else {
        theInfo->align_options += TXALIGN_MATRIX_VAL;
        theInfo->align_options += TXALIGN_SHOW_QS;
    }
    
    /* Always HTML in WWW case */
    
    theInfo->align_options += TXALIGN_HTML;
    theInfo->print_options += TXALIGN_HTML;
    
    return TRUE;
}

Boolean WWWValidateOptions(WWWBlastInfoPtr theInfo)
{
    ValNodePtr error_return=NULL;
    BlastErrorMsgPtr error_msg;
    CharPtr msg = NULL;

    if(theInfo == NULL || theInfo->options == NULL)
	return FALSE;
    
    if(BLASTOptionValidateEx(theInfo->options, theInfo->program, 
                             &error_return) != 0) {
	if (error_return) {
	    error_msg = (BlastErrorMsgPtr) error_return->data.ptrvalue;
            msg = error_msg->msg;
	}
        WWWBlastErrMessage(BLASTErrOptions, msg);

	return FALSE;
    }
    return TRUE;
}

/* Used for PSI/PHI Blast searches */

static Int1 S62ToInt(Uint1 ch)
{
    if(isupper(ch)) /* negative value */
	return(64 - ch);
    else if (isdigit(ch)) /* positive less than 10 */
	return(ch - 48); 
    else if (!isupper(ch)) /* positive more or eq  10 */
	return(ch - 87);
    return 0;
}

static Uint1 IntTo62S(Int1 value)
{
    if(value < 0)
	return(64-value);
    else if (value < 10)
	return(value + 48); 
    else if (value < 36)
	return(value + 87);
    return 0;
}
static Int4Ptr PNTR Decode62Matrix(CharPtr Matrix62, Int4 length, Int4 size)
{
    Int4Ptr PNTR posMatrix;
    register Int4 i, j, k = 0;

    if(Matrix62 == NULL || Matrix62[0] == NULLB)
	return NULL;

    posMatrix = (Int4Ptr PNTR) MemNew(sizeof(Int4)*(length+1));

    for(i = 0; i <= length; i++)
	posMatrix[i] = (Int4Ptr) MemNew(sizeof(Int4Ptr)*size);

    for(i = 0; i <= length; i++) {
	for(j =0; j < size; j++) {
	    if(Matrix62[k] == 'z')
		posMatrix[i][j] = BLAST_SCORE_MIN;
	    else if (Matrix62[k] == 'Z')
		posMatrix[i][j] = BLAST_SCORE_MAX; 
	    else
		posMatrix[i][j]= S62ToInt(Matrix62[k]);
	    k++;
	}
    }
    return posMatrix;
}

static CharPtr Encode62Matrix(Int4Ptr PNTR posMatrix, Int4 length, Int4 size)
{
    register Int4 i, j, k=0;
    CharPtr Matrix62;

    Matrix62 = (CharPtr) MemNew((length+1)*size+1);

    for(i = 0; i <= length; i++) {
	for(j =0; j < size; j++) {

	    if(posMatrix[i][j] < -26)
		Matrix62[k] = 'z';
	    else if (posMatrix[i][j] > 35)
		Matrix62[k] = 'Z';
	    else
		Matrix62[k] = IntTo62S(posMatrix[i][j]);

	    k++;
	}
    }
    return Matrix62;
}

void BLASTPrintDataFree(BLASTPrintDataPtr data)
{
    GIListPtr glp, glp_next;
    ValNodePtr vnp;

    if(data == NULL)
	return;

    TxDfDbInfoDestruct(data->dbinfo);
    MemFree(data->ka_params);
    MemFree(data->ka_params_gap);
    MemFree(data->buffer);
    ValNodeFreeData(data->info_vnp);
    
    if(data->psidata != NULL) {
        MemFree(data->psidata->matrix62);

        for(glp = data->psidata->PrevCheckedGIs; glp != NULL; glp = glp_next) {
            glp_next = glp->next;
            MemFree(glp);
        }

        for(glp = data->psidata->PrevGoodGIs; glp != NULL; glp = glp_next) {
            glp_next = glp->next;
            MemFree(glp);
        }
        MemFree(data->psidata);
    }
    
    if(data->seqalign != NULL)
        SeqAlignSetFree(data->seqalign);
    
    SeqLocFree(data->seqloc);
    
    for(vnp = data->vnp; vnp != NULL; vnp=vnp->next) {
        SeqAlignSetFree((SeqAlignPtr) vnp->data.ptrvalue);
    }
    
    ValNodeFree(data->vnp);
    
    MemFree(data);
    return;
}

/* We got seqalign list, now divide it into two lists:
   the first one will contain alignments with Evalue
   better than threshold, other worse than threshold
*/

Boolean 
SplitSeqAlign(SeqAlignPtr seqalign, SeqAlignPtr *GoodSeqAlignment_ptr, 
              SeqAlignPtr *BadSeqAlignment_ptr, SeqAlignPtr *lastGood_ptr, 
              Int2Ptr *marks_ptr, Int2Ptr countBad_ptr, 
              Int2Ptr countGood_ptr, Nlm_FloatHi ethresh_old)
{
    
    Boolean first_time;
    SeqIdPtr last_id, subject_id;
    SeqAlignPtr         gsl, seqalign_var, last_seqalign;
    SeqAlignPtr		BadSeqAlignments, GoodSeqAlignments, lastGood = NULL;
    Nlm_FloatHi         bit_score, evalue;
    Int2Ptr		marks;
    Int2		countGood, countBad;
    Int4                number, score;

    last_id = NULL;
    first_time = TRUE;
    GoodSeqAlignments = seqalign_var = seqalign;
    
    BadSeqAlignments = NULL;
    while (seqalign_var) {
	subject_id = SeqIdDup(TxGetSubjectIdFromSeqAlign(seqalign_var));
	if (last_id == NULL || SeqIdComp(subject_id, last_id) == SIC_NO) {
            SeqIdSetFree(last_id);
            last_id = subject_id;
            GetScoreAndEvalue(seqalign_var, &score, &bit_score, &evalue, &number);
            if (evalue > ethresh_old) {
                if (first_time == TRUE) {
                    GoodSeqAlignments = NULL;
                    lastGood = NULL;
                    last_seqalign = NULL;	/* split the good and bad lists. */
                } else {
                    lastGood = last_seqalign;
                    last_seqalign->next = NULL;	/* split the good and bad lists. */
                }
                BadSeqAlignments = seqalign_var;
                break;	
            }
	} else {
            SeqIdSetFree(subject_id);
        }

	first_time = FALSE;
	last_seqalign = seqalign_var;
	seqalign_var = seqalign_var->next;
    }
    
    /* count number of good and bad alignments */
    for (gsl = GoodSeqAlignments, countGood = 0; gsl; gsl = gsl->next, 
             countGood++);
    for (gsl = BadSeqAlignments, countBad = 0; gsl; gsl = gsl->next, 
             countBad++);
    
    if (countGood + countBad)
	/* allocate memo for marks array */
	marks = (Int2Ptr) MemNew(sizeof(Int2) * (countGood + countBad));
    else
	marks = NULL;
    
    *GoodSeqAlignment_ptr = GoodSeqAlignments;
    *BadSeqAlignment_ptr = BadSeqAlignments;
    *lastGood_ptr = lastGood;
    *marks_ptr = marks;
    *countBad_ptr = countBad;
    *countGood_ptr = countGood;
    
    return TRUE;
}

static Boolean TestSTDOut(void)
{
    if(write(1, "", 0) < 0) {
	return FALSE;
    }
    return TRUE;
}

static int LIBCALLBACK WWWTickCallback(Int4 sequence_number, 
                                       Int4 number_of_positive_hits)
{
    if(!TestSTDOut()) {
	return -1;
    }
    
    fprintf(stdout, ".");
    fflush(stdout);
    
    return 1;
}

static void PrintRequestHeader(WWWBlastInfoPtr theInfo)
{
    Char TimeNowStr[30], TimeModStr[30];
    struct stat buf;

    printf("<HTML>\n");
    printf("<HEAD>\n");
    printf("<TITLE>BLAST Search Results </TITLE>\n");
    printf("</HEAD>\n");
    printf("<BODY BGCOLOR=\"#FFFFFF\" LINK=\"#0000FF\" "
           "VLINK=\"#660099\" ALINK=\"#660099\">\n");
    printf("<A HREF=\"%s/blast_form.map\">"
           "<IMG SRC=\"%s/images/psi_blast.gif\" "
           "BORDER=0 ISMAP></A>\n", 
           theInfo->www_root_path == NULL? "/blast" : theInfo->www_root_path, 
           theInfo->www_root_path == NULL? "/blast" : theInfo->www_root_path);
    printf("<BR><BR><PRE>\n");
   
    init_buff_ex(90);

    BlastPrintVersionInfo(theInfo->program, TRUE, stdout);
    fprintf(stdout, "\n");
    
    if (theInfo->is_phi_blast) {
	BlastPrintPhiReference(TRUE, 90, stdout);
    } else {
    	BlastPrintReference(TRUE, 90, stdout);
    }
    
    fprintf(stdout, "\n");

    AcknowledgeBlastQuery(theInfo->query_bsp, 70, 
                          stdout, theInfo->believe_query, TRUE);
    
    /*
       dbinfo = GetDbInfoFromReadDb(search_data->database, 
       !search_data->db_is_na);
       PrintDbInformation(dbinfo, 70, search_data->outfp, 
       search_data->html);
     */

    PrintDbInformation(theInfo->database, !theInfo->db_is_na, 70, 
                       stdout, TRUE);
    free_buff();

    fprintf(stdout, "Searching");
    
    fflush(stdout);
    
    return;
}

static ValNodePtr seed_core_private (BlastSearchBlkPtr search, CharPtr program_name, BLAST_OptionsBlkPtr options, SeqLocPtr *seqloc_ptr, CharPtr patfile, CharPtr pattern, Uint1Ptr query, Uint1Ptr unfilter_query, Int4 queryLength, Boolean show_diagnostics, Nlm_FloatHi *paramC, ValNodePtr PNTR info_vnp)
     
{
    Boolean tmp_file_made = FALSE;
    Char buffer[PATH_MAX];
    ImpFeatPtr ifp;
    Int4 hitlist_count, hspcnt, index, index1;
    Int4 number_of_descriptions, number_of_alignments;
    ObjectIdPtr obidp;
    posSearchItems *posSearch;
    SeqFeatPtr sfp;
    SeqIdPtr subject_id;
    SeqLocPtr next, seqloc;

    FILE *patfp; 
    seedSearchItems *seedSearch;
    Int4 program_flag;
    Boolean is_dna = FALSE; /*cannot use DNA queries in blastpgp*/
    Int4 i; /*index over characters*/
    ValNodePtr vnp;
    
    
    if (search == NULL)
        return NULL;
    
    program_flag = convertProgramToFlag(program_name, &is_dna);
    
    if (options->isPatternSearch) {
        if (pattern) {
            /* open and fill a temporary file if there's a pattern. XXX */
            TmpNam(buffer);
            patfp = FileOpen(buffer, "w");
            fprintf(patfp, "ID \n");
            fprintf(patfp, "PA   %s\n", pattern);
            fflush(patfp);
            FileClose(patfp);
            patfp = NULL;
            tmp_file_made = TRUE;
        }

        if (patfile)	/* If a file was give, use it. */
            StringCpy(buffer, patfile);
        if ((patfp = FileOpen(buffer, "r")) == NULL) {
            ErrPostEx(SEV_FATAL, 0, 0, "blast: Unable to open pattern file %s\n", buffer);
            return NULL;
        }
        seedSearch = (seedSearchItems *) MemNew(sizeof(seedSearchItems));
        fillCandLambda(seedSearch, options->matrix, options);
    } else {
        ErrPostEx(SEV_FATAL, 0, 0, "Must be a pattern search");
        return NULL;
    }
    
    if (paramC)
        *paramC = seedSearch->paramC;
    
    search->gap_align = GapAlignBlkNew(1,1);
    search->gap_align->gap_open = options->gap_open;
    search->gap_align->gap_extend = options->gap_extend;
    search->gap_align->decline_align = (-(BLAST_SCORE_MIN));
    search->gap_align->x_parameter = options->gap_x_dropoff
        *NCBIMATH_LN2/seedSearch->paramLambda;
    search->gap_align->matrix = search->sbp->matrix;
    initProbs(seedSearch);
    init_order(search->gap_align->matrix,program_flag,is_dna,seedSearch);

    for(i = 0; i < queryLength; i++)
        query[i] = seedSearch->order[query[i]];
    if (unfilter_query) {
        for(i = 0; i < queryLength; i++)
            unfilter_query[i] = seedSearch->order[unfilter_query[i]];
    }
    
    seqloc = NULL;
    posSearch = (posSearchItems *) MemNew(sizeof(posSearchItems));
    
    vnp = seedEngineCore(search, options, query, unfilter_query, readdb_get_filename(search->rdfp), patfile, program_flag, patfp, is_dna, FALSE, seedSearch, options->ethresh, 0.0, posSearch, &seqloc, show_diagnostics, info_vnp);
    
    if (tmp_file_made)	/* Remove temporary pattern file if it exists. */
        FileRemove(buffer);
    
    MemFree(seedSearch);

    MemFree(posSearch->posResultSequences);
    MemFree(posSearch);

    *seqloc_ptr = seqloc;
    
    return vnp;
}

static  BLAST_ScorePtr PNTR GetPSIMatrix(PSIDataPtr psidata, WWWBlastInfoPtr theInfo, Nlm_FloatHi *karlinK_out, Nlm_FloatHi *karlinK_nogap) {
    
    CharPtr     chptr, Matrix62_last=NULL, pattern;
    Int4        i, j;
    Nlm_FloatHi karlinK;
    FILE *fd;
    Int2        num_entries;
    Int4 queryLength; /*length of query sequence*/
    Int4 numSeqAligns;
    SeqAlignPtr	seqalign;
    SeqAlignPtr *lastSeqAligns=NULL;
    SeqLocPtr seg_slp; /*pointer to structure for seg filtering*/
    Uint1Ptr query = NULL; /*query sequence read in*/
    Uint1Ptr unfilter_query = NULL; /*needed if seg will filter query*/
    ValNodePtr info_vnp, phi_vnp = NULL;
    SeqLocPtr phi_seqloc = NULL;
    Uint4      gi;
    SeqIdPtr   sip, all_sip = NULL;

    BlastSearchBlkPtr search;
    BLAST_ScorePtr PNTR posMatrix;
    compactSearchItems *compactSearch = NULL;
    posSearchItems *posSearch;

    /* first step; return NULL, means to use default matrix */
    if (psidata->StepNumber == 0) 
	return NULL;
        
    /* The second step;  means, that there is list of GI's to 
       recalculate matrix 
       using default or read matrix and list of seqaligns 
       from limited search;  limited search -
       search the query not in whole database but just
       in subset specified byt the list of GI's
     */

    /* Create list of SeqIdPtr for the limited search */

    num_entries = WWWGetNumEntries(theInfo->info);

    for(i = 0; i < num_entries; i++) {
	GIListPtr	good_gil, checked_gil;
        
	if((chptr = WWWGetNameByIndex(theInfo->info, i)) != NULL && 
           !StringICmp(chptr, "checked_GI")) {
            
	    if((chptr = WWWGetValueByIndex(theInfo->info, i)) != NULL) {
                gi = atoi(chptr);                
		sip = ValNodeAddInt(NULL, SEQID_GI, gi);
                ValNodeLink(&all_sip, sip);

		/* Create list of GI's, it will be used when we
		   test convergence */
                
		if (!psidata->PrevCheckedGIs) {
		    /* first one */
		    checked_gil = (GIListPtr) MemNew(sizeof(GIList));
		    checked_gil->gi = gi;
		    checked_gil->next = NULL;
		    psidata->PrevCheckedGIs = checked_gil;
		} else {
		    checked_gil->next = (GIListPtr) MemNew(sizeof(GIList));
		    checked_gil = checked_gil->next;
		    checked_gil->gi = gi;
		    checked_gil->next = NULL;
		}
	    }
	} else if((chptr = WWWGetNameByIndex(theInfo->info, i)) != NULL && 
                  !StringICmp(chptr, "good_GI")) {
            
	    if((chptr = WWWGetValueByIndex(theInfo->info, i)) != NULL) {
                gi = atoi(chptr);
                
		if (!psidata->PrevGoodGIs) {
		    /* first one */
		    good_gil = (GIListPtr) MemNew(sizeof(GIList));
		    good_gil->gi = gi;
		    good_gil->next = NULL;
		    psidata->PrevGoodGIs = good_gil;
		} else {
		    good_gil->next = (GIListPtr) MemNew(sizeof(GIList));
		    good_gil = good_gil->next;
		    good_gil->gi = gi;
		    good_gil->next = NULL;
		}
	    }
	}
    }

    /* So this search will be limited to the list of gis */
    theInfo->options->gilist = all_sip;
    
    /* Some additional parameters required for the matrix recalculation */
    
    theInfo->options->use_best_align = TRUE;
    theInfo->options->use_real_db_size = TRUE;
    
    /* the search */

    if((search = BLASTSetUpSearchWithReadDb(theInfo->fake_bsp, theInfo->program, theInfo->query_bsp->length, theInfo->database, theInfo->options, WWWTickCallback ))  == NULL) {
	return NULL;
    } 
    
    search->positionBased = FALSE;
    
    if (psidata->StepNumber > 1) {
	/* The third and rest of steps;  means, that
	   before we are recalculating matrix, we 
	   should read old matrix from the previous step;
        */
        
	/* Read matrix from the previous step */
        
        Matrix62_last = WWWGetValueByName(theInfo->info, "PSI_MATRIX");
        
        if((chptr = WWWGetValueByName(theInfo->info, 
                                      "PSI_KARLIN_K")) != NULL) {
            karlinK = atof(chptr);
            search->sbp->kbp_gap_psi[0]->K = karlinK;
            search->sbp->kbp_gap_psi[0]->logK = log(karlinK);
        }
	
	/* Decode read matrix */
        
	if(Matrix62_last != NULL && Matrix62_last[0] != NULLB) {
	    search->positionBased = TRUE;
	    search->sbp->posMatrix = 
		Decode62Matrix (Matrix62_last, 
                                search->context[0].query->length, 
                                search->sbp->alphabet_size);
            
	}
    } /* end reread the matrix */    

    search->thr_info->tick_callback = NULL;
    
    pattern = WWWGetValueByName(theInfo->info, "PHI_PATTERN");
    
    /* If pattern is non-NULL, then it is a PHI-BLAST search. */
    if (pattern) {
        query = BlastGetSequenceFromBioseq(theInfo->fake_bsp, &queryLength);
        seg_slp = BlastBioseqFilter(theInfo->fake_bsp, theInfo->options->filter_string);
        unfilter_query = NULL;
        if (seg_slp) {
            unfilter_query = MemNew((queryLength + 1) * sizeof(Uint1));
            for (i = 0; i < queryLength; i++)
                unfilter_query[i] = query[i];
            BlastMaskTheResidues(query,queryLength,21,seg_slp,FALSE, 0);
        }
        
	theInfo->options->isPatternSearch = TRUE;
	phi_vnp = seed_core_private(search, "patseedp", theInfo->options, &phi_seqloc, NULL, pattern, query, unfilter_query, queryLength, FALSE, NULL, &info_vnp);
        ValNodeFreeData(info_vnp);

        MemFree(query);
        MemFree(unfilter_query);
        
	seqalign = convertValNodeListToSeqAlignList(phi_vnp, &lastSeqAligns, &numSeqAligns);
        ValNodeFree(phi_vnp);
    } else {
    	seqalign = BioseqBlastEngineCore(search, theInfo->options, 
                                         search->sbp->posMatrix);
    }
    
    if(search->sbp->posMatrix != NULL) {
        for(i = 0; i <= theInfo->fake_bsp->length; i++) {
            MemFree(search->sbp->posMatrix[i]);
        }    
        MemFree(search->sbp->posMatrix);
        search->sbp->posMatrix = NULL;
    }

    /* Now finaly calculating matrix that will be used at this step */
    
    if(seqalign) {

        ReadDBBioseqFetchEnable("psiblast", theInfo->database, 
                                FALSE, TRUE);
        compactSearch = compactSearchNew(compactSearch);
        copySearchItems(compactSearch, search);

        compactSearch->pseudoCountConst = 10;
        
        posMatrix = WposComputation(compactSearch, seqalign);

        MemFree(compactSearch->standardProb);
        MemFree(compactSearch);

        ReadDBBioseqFetchDisable();
        
	/* Encode matrix for the use in the next step*/
        
	psidata->matrix62 = Encode62Matrix(posMatrix, 
                                           search->context[0].query->length, 
                                           search->sbp->alphabet_size);
    }

    *karlinK_out = search->sbp->kbp_gap_psi[0]->K;
    *karlinK_nogap = search->sbp->kbp[0]->K;

    SeqAlignSetFree(seqalign);
    SeqLocFree(phi_seqloc);
    
    search = BlastSearchBlkDestruct(search);

    theInfo->options->use_best_align = FALSE;
    theInfo->options->use_real_db_size = FALSE;

    SeqIdSetFree(theInfo->options->gilist);
    theInfo->options->gilist = NULL;
    theInfo->options->isPatternSearch = FALSE;
    
    return posMatrix;
} /* end of GetPSIMatrix() */

BLASTPrintDataPtr PSIBlastSearch(WWWBlastInfoPtr theInfo)
{
    BLASTPrintDataPtr print_data;
    ValNodePtr vnp, other_returns= NULL;
    Int4 i, num_entries;
    PSIDataPtr psidata;
    CharPtr    chptr;
    Char       matrixname[64] = "BLOSUM62";
    Int4       opencost = 0, extendedcost = 0;
    BlastSearchBlkPtr search;
    BLAST_ScorePtr PNTR posMatrix;
    Nlm_FloatHi karlinK;
    Nlm_FloatHi karlinK_nogap;
    BLAST_MatrixPtr  matrix = NULL;

    if(theInfo == NULL)
	return NULL;

    PrintRequestHeader(theInfo); 
    print_data = (BLASTPrintDataPtr) MemNew(sizeof(BLASTPrintData));
    
    psidata = MemNew(sizeof(PSIData));
    psidata->PrevGoodGIs = NULL;
    psidata->PrevCheckedGIs = NULL;
    
    /* initialize the search */
    theInfo->options->pseudoCountConst = 10;

    if((search = BLASTSetUpSearchWithReadDb(theInfo->fake_bsp, theInfo->program, theInfo->query_bsp->length, theInfo->database, theInfo->options, WWWTickCallback ))  == NULL) {
	return NULL;
    }

    /* Matrix and StepNumber for PSI-Blast */
    
    if((chptr = WWWGetValueByName(theInfo->info, "STEP_NUMBER")) != NULL)
	psidata->StepNumber = atoi(chptr);
    
    if((posMatrix = GetPSIMatrix(psidata, theInfo, 
                                 &karlinK, &karlinK_nogap)) != NULL) {
        search->positionBased = TRUE;
        search->sbp->kbp_gap_psi[0]->K = karlinK;
        search->sbp->kbp_gap_psi[0]->logK = log(karlinK);
        search->sbp->kbp[0]->K = karlinK_nogap;
        search->sbp->kbp[0]->logK = log(karlinK_nogap);
    }
    
    search->sbp->posMatrix = posMatrix;
    
    search->thr_info->tick_callback = WWWTickCallback;
    
    print_data->seqalign = BioseqBlastEngineCore(search, theInfo->options, posMatrix);

    if(posMatrix != NULL) {
        for(i = 0; i <= theInfo->fake_bsp->length; i++) {
            MemFree(posMatrix[i]);
        }
        
        MemFree(posMatrix);
        search->sbp->posMatrix = NULL;
    }

    print_data->psidata = psidata;
    
    /*  Blast search */
    
    other_returns = BlastOtherReturnsPrepare(search);

    print_data->mask_loc = NULL;

    for (vnp=other_returns; vnp; vnp = vnp->next) {
	switch (vnp->choice) {
        case TXDBINFO:
            print_data->dbinfo = vnp->data.ptrvalue;
            break;
        case TXKABLK_NOGAP:
            print_data->ka_params = 
                (BLAST_KarlinBlkPtr) vnp->data.ptrvalue;
            break;
        case TXKABLK_GAP:
            print_data->ka_params_gap = 
                (BLAST_KarlinBlkPtr) vnp->data.ptrvalue;
            psidata->karlinK = print_data->ka_params_gap->K;
            break;
        case TXPARAMETERS:
            print_data->buffer = vnp->data.ptrvalue;
            break;
        case TXMATRIX:
            matrix = (BLAST_MatrixPtr) vnp->data.ptrvalue;
            BLAST_MatrixDestruct(matrix);
            vnp->data.ptrvalue = NULL;
            break;
        case SEQLOC_MASKING_NOTSET:
        case SEQLOC_MASKING_PLUS1:
        case SEQLOC_MASKING_PLUS2:
        case SEQLOC_MASKING_PLUS3:
        case SEQLOC_MASKING_MINUS1:
        case SEQLOC_MASKING_MINUS2:
        case SEQLOC_MASKING_MINUS3:
            ValNodeAddPointer(&(print_data->mask_loc), vnp->choice, vnp->data.ptrvalue);
            break;
        default:
            break;
	}
    }   

    ValNodeFree(other_returns);

    search = BlastSearchBlkDestruct(search);

    return print_data;
}

BLASTPrintDataPtr PHIBlastSearch(WWWBlastInfoPtr theInfo)
{
    BLASTPrintDataPtr print_data;
    ValNodePtr vnp, other_returns= NULL;
    PSIDataPtr psidata;
    Nlm_FloatHi paramC;
    Int4 i, num_entries;
    Int4 queryLength; /*length of query sequence*/
    CharPtr chptr;
    Char	matrixname[64] = "BLOSUM62";
    Int4	opencost = 0, extendedcost = 0;
    SeqLocPtr seqloc=NULL;
    SeqLocPtr seg_slp; /*pointer to structure for seg filtering*/
    Uint1Ptr query = NULL; /*query sequence read in*/
    Uint1Ptr unfilter_query = NULL; /*needed if seg will filter query*/
    ValNodePtr info_vnp;
    BlastSearchBlkPtr search;
    BLAST_MatrixPtr  matrix = NULL;

    if(theInfo == NULL)
	return NULL;
    
    PrintRequestHeader(theInfo); 
    print_data = (BLASTPrintDataPtr) MemNew(sizeof(BLASTPrintData));
    
    psidata = MemNew(sizeof(PSIData));
    psidata->PrevGoodGIs = NULL;
    psidata->PrevCheckedGIs = NULL;
    
    print_data->psidata = psidata;
    
    /* Get matrix name and gap costs */
    
    if((chptr = WWWGetValueByName(theInfo->info, "MAT_PARAM")) != NULL) {
	if (chptr[1] != '-' || chptr[2] != '-')
	    sscanf(chptr, "%s\t %ld\t %ld", matrixname, &opencost, &extendedcost);
    }

    /* Change matrix parameters */
    
    BLASTOptionSetGapParams (theInfo->options, matrixname, 
                             opencost, extendedcost);
    
    chptr = WWWGetValueByName(theInfo->info, "PHI_PATTERN");
    
    /* Reguilar PHI-Blast search */
    theInfo->options->isPatternSearch = TRUE;
    
    if((search = BLASTSetUpSearchWithReadDb(theInfo->fake_bsp, "blastp", theInfo->query_bsp->length, theInfo->database, theInfo->options, WWWTickCallback))  == NULL) {
	return NULL;
    }
    
    query = BlastGetSequenceFromBioseq(theInfo->fake_bsp, &queryLength);
    seg_slp = BlastBioseqFilter(theInfo->fake_bsp, 
                                theInfo->options->filter_string);
    
    unfilter_query = NULL;
    if (seg_slp) {
        unfilter_query = MemNew((queryLength + 1) * sizeof(Uint1));
        for (i = 0; i < queryLength; i++)
            unfilter_query[i] = query[i];
        BlastMaskTheResidues(query,queryLength,21,seg_slp,FALSE, 0);
    }
    
    print_data->vnp = seed_core_private(search, "patseedp", theInfo->options, &(print_data->seqloc), NULL, chptr, query, unfilter_query, queryLength, TRUE, &paramC, &info_vnp);
    
    print_data->info_vnp = info_vnp;
    
    MemFree(query);
    MemFree(unfilter_query);
    
    /*  Blast search */
    
    other_returns = BlastOtherReturnsPrepare(search);
    print_data->mask_loc = NULL;
    for (vnp=other_returns; vnp; vnp = vnp->next) {
	switch (vnp->choice) {
        case TXDBINFO:
            print_data->dbinfo = vnp->data.ptrvalue;
            break;
        case TXKABLK_GAP:
            print_data->ka_params_gap = vnp->data.ptrvalue;
            /* print_data->ka_params_gap->paramC = paramC; ?? */
            break;
        case TXKABLK_NOGAP:
            print_data->ka_params = vnp->data.ptrvalue;
            break;
        case TXPARAMETERS:
            print_data->buffer = vnp->data.ptrvalue;
            break;
        case TXMATRIX:
            matrix = (BLAST_MatrixPtr) vnp->data.ptrvalue;
            BLAST_MatrixDestruct(matrix);
            vnp->data.ptrvalue = NULL;
            break;
        case SEQLOC_MASKING_NOTSET:
        case SEQLOC_MASKING_PLUS1:
        case SEQLOC_MASKING_PLUS2:
        case SEQLOC_MASKING_PLUS3:
        case SEQLOC_MASKING_MINUS1:
        case SEQLOC_MASKING_MINUS2:
        case SEQLOC_MASKING_MINUS3:
            ValNodeAddPointer(&(print_data->mask_loc), vnp->choice, vnp->data.ptrvalue);
            break;
        default:
            break;
        }
    }
    
    search = BlastSearchBlkDestruct(search);
    
    fflush(stdout);
    return print_data;
}

static	void	printSubmitButton(FILE* fp, Int4 step)
{
    fprintf(fp, "<INPUT TYPE=\"submit\" NAME=\"NEXT_I\" "
            "VALUE=\"Run PSI-Blast iteration %d\">\n", 
            step);
}

static Int4 get_number_alignment(SeqAlignPtr align)
{
    Int4 num = 0;
    
    while(align) {
	++num;
	align = align->next;
    }
    
    return num;
}

Boolean PHIPrintOutput(WWWBlastInfoPtr theInfo,
	BLASTPrintDataPtr print_data, 
	ValNodePtr vnp, Nlm_FloatHi ethresh_old)
{
    Uint4 align_options, print_options;
    SeqAnnotPtr seqannot;
    BlastTimeKeeper time_keeper;
    BlastPruneSapStructPtr prune;
    Uint1 f_order[FEATDEF_ANY], g_order[FEATDEF_ANY];
    Char	hostname[30], buffer[32];
    Char	href[1024];
    CharPtr	chptr;
    Char	f_name[64], title[1024];
    Int4	align_num;
    Int2	count, countBad, countGood;
    Int2Ptr 	marks;
    Int4 numSeqAligns;
    SeqAlignPtr	seqalign;
    SeqAlignPtr *lastSeqAligns=NULL;
    SeqAlignPtr	lastGood, BadSeqAlignments, GoodSeqAlignments;
    SeqLocPtr	seqloc;
    ValNodePtr	vnp_var;

    MemSet((Pointer)(g_order), 0, (size_t)(FEATDEF_ANY* sizeof(Uint1)));
    MemSet((Pointer)(f_order), 0, (size_t)(FEATDEF_ANY* sizeof(Uint1)));
    
    if(print_data == NULL) {
        WWWBlastErrMessage(BLASTMiscError, NULL);                
	return FALSE;
    }

    print_options = 0;
    align_options = 0;

    align_options += TXALIGN_COMPRESS;
    align_options += TXALIGN_END_NUM;

    if (theInfo->show_gi) {
	align_options += TXALIGN_SHOW_GI;
	print_options += TXALIGN_SHOW_GI;
    }

    if (theInfo->options->gapped_calculation == FALSE)
	print_options += TXALIGN_SHOW_NO_OF_SEGS;


    if (theInfo->align_view) {
	align_options += TXALIGN_MASTER;

	if (theInfo->align_view == 1 || theInfo->align_view == 3)
	    align_options += TXALIGN_MISMATCH;
        
	if (theInfo->align_view == 3 || theInfo->align_view == 4 || 
            theInfo->align_view == 6)
	    align_options += TXALIGN_FLAT_INS;

	if (theInfo->align_view == 5 || theInfo->align_view == 6)
	    align_options += TXALIGN_BLUNT_END;
    } else {
	align_options += TXALIGN_MATRIX_VAL;
	align_options += TXALIGN_SHOW_QS;
    }

    /* align_options += TXALIGN_MATRIX_VAL;
       align_options += TXALIGN_SHOW_QS; */

    align_options += TXALIGN_HTML;
    print_options += TXALIGN_HTML; 

    ReadDBBioseqFetchEnable ("phiblast", 
                             theInfo->database,  theInfo->db_is_na, TRUE);
    
    seqannot = SeqAnnotNew();
    seqannot->type = 2;
    AddAlignInfoToSeqAnnot(seqannot, theInfo->align_type);

    init_buff();

    /*    gethostname(hostname, sizeof(hostname)); */

    sprintf(href, "%s/nph-viewgif.cgi?", theInfo->www_root_path == NULL? "/blast" : theInfo->www_root_path);
    
    seqalign = convertValNodeListToSeqAlignList(print_data->vnp, &lastSeqAligns, &numSeqAligns);
    seqannot->data = seqalign;
    if (theInfo->show_overview) {
        sprintf(f_name, "%ld%ld.gif", (long)random(), (long)getpid());
        align_num = get_number_alignment((SeqAlignPtr)(seqannot->data)); 
        sprintf(title, "<H3><a href=\"%s/docs/newoptions.html#graphical-overview\"> "
                "Distribution of %ld Blast Hits on the Query Sequence</a></H3>\n", theInfo->www_root_path == NULL? "/blast" : theInfo->www_root_path, (long)align_num);  
        
        PrintAlignmentOverview(seqannot, stdout, "PSI_BLAST", href, f_name, title); 
    }

    seqannot->data = NULL;
    seqannot = SeqAnnotFree(seqannot);
    
    print_data->vnp = convertSeqAlignListToValNodeList(seqalign, lastSeqAligns, numSeqAligns);
    
    print_options += TXALIGN_DO_NOT_PRINT_TITLE; 
    print_options += TXALIGN_CHECK_BOX;
    if (print_data->psidata->StepNumber)
        print_options += TXALIGN_NEW_GIF;
    
    /* submit button */
    printSubmitButton(stdout, 
                      print_data->psidata->StepNumber+1);
    
    if (print_data->psidata->StepNumber && theInfo->number_of_descriptions) {
        printf("<HR><p><b>Legend:</b><p>\
<IMG SRC=\"%s/images/new.gif\" WIDTH=25 HEIGHT=15> - means that \
the alignment score was below the threshold on the previous iteration \
<p>\
<IMG SRC=\"%s/images/checked.gif\" WIDTH=15 HEIGHT=15> - means that \
the alignment was checked on the previous iteration \
</p>", theInfo->www_root_path == NULL? "/blast" : theInfo->www_root_path,
               theInfo->www_root_path == NULL? "/blast" : theInfo->www_root_path);
    }
    
    /*
      if (theInfo->number_of_descriptions) {
      if (print_data->psidata->StepNumber)
      printf("\n<IMG SRC=\"/BLAST/bg.gif\" WIDTH=65 HEIGHT=15>");
      printf("                                                                     Score    E");
      printf("\nSequences producing significant alignments:");
      if (print_data->psidata->StepNumber)
      printf("<IMG SRC=\"/BLAST/bg.gif\" WIDTH=65 HEIGHT=15>");
      printf("                          (bits) Value\n\n");
      }
    */
    
    vnp_var = vnp;
    seqloc = print_data->seqloc;
    marks = NULL;
    while (vnp_var) {
        SplitSeqAlign((SeqAlignPtr) vnp_var->data.ptrvalue, &GoodSeqAlignments, &BadSeqAlignments, &lastGood, &marks, 
                      &countBad, &countGood, ethresh_old);
	
        printf("<HR><CENTER><b><FONT color=\"green\">"
               "Sequences with pattern at position %d and E-value BETTER than threshold</FONT></b></CENTER>\n",
               SeqLocStart(seqloc)+1);
        
        if (print_data->psidata->StepNumber)
            printf("\n<IMG SRC=\"%s/images/bg.gif\" WIDTH=65 HEIGHT=15>", theInfo->www_root_path == NULL? "/blast" : theInfo->www_root_path);
        printf("                                                                     Score    E");
        printf("\nSequences producing significant alignments:");
        if (print_data->psidata->StepNumber)
            printf("<IMG SRC=\"%s/images/bg.gif\" WIDTH=65 HEIGHT=15>",
                   theInfo->www_root_path == NULL? "/blast" : theInfo->www_root_path);
        printf("                          (bits) Value\n\n");
        
        fflush(stdout);
        print_options += TXALIGN_CHECK_BOX_CHECKED;
        PrintDefLinesFromSeqAlignEx(GoodSeqAlignments, 80, stdout, print_options, FIRST_PASS, marks, theInfo->number_of_descriptions);
        
        print_options -= TXALIGN_CHECK_BOX_CHECKED;
        
        if (print_data->psidata->StepNumber == 0)
            printf("<a name = Evalue> </a>");
	
        if (theInfo->number_of_descriptions > countGood && BadSeqAlignments) {
            /* submit button */
            printSubmitButton(stdout, print_data->psidata->StepNumber+1);
            
            printf("<HR><CENTER><b><FONT color=\"green\">"
                   "Sequences with pattern at position %d and E-value WORSE than threshold</FONT></b></CENTER>\n", 
                   SeqLocStart(seqloc)+1);
            
            PrintDefLinesFromSeqAlignEx(BadSeqAlignments, 80, stdout, print_options, FIRST_PASS, &marks[countGood], theInfo->number_of_descriptions - countGood);
        }

        marks = MemFree(marks);
        
        /* merge lists */
        if (lastGood)
            lastGood->next = BadSeqAlignments;
        
        vnp_var = vnp_var->next;
        seqloc = seqloc->next;
    }
    
    if (theInfo->number_of_descriptions) {
        /* submit button */
        printSubmitButton(stdout, 
                          print_data->psidata->StepNumber+1);
    }
    
    free_buff();
    fflush(stdout);

    GlobalAlignNumber = 0;

    fprintf(stdout, "<HR>");

    if (theInfo->number_of_alignments) {
	fprintf(stdout, "<CENTER><b><FONT color=\"green\">"
		"Alignments</FONT></b></CENTER>\n");
        
        f_order[FEATDEF_REGION] = 1;
        g_order[FEATDEF_REGION] = 1;
        if(theInfo->align_view == 0) {
            ShowTextAlignFromAnnotExtra(theInfo->query_bsp, 
                                        print_data->vnp, 
                                        print_data->seqloc, 60, 
                                        stdout, 
                                        f_order, g_order, align_options, 
                                        NULL, print_data->mask_loc, 
                                        FormatScoreFunc);
        } else {
            ShowTextAlignFromAnnotExtra(theInfo->query_bsp, 
                                        print_data->vnp, 
                                        print_data->seqloc, 60, 
                                        stdout, 
                                        f_order, g_order, align_options, 
                                        NULL, print_data->mask_loc, 
                                        NULL);
            printf("<P>\n");
        }
    }

    fflush(stdout);
    ObjMgrClearHold(); 

    printf("<PRE>\n");

    BlastTimeFillStructure(&time_keeper);

    fprintf(stdout, "CPU time: %8.2f user secs.\t%8.2f sys. "
	    "secs\t%8.2f total secs.\n\n", 
	    time_keeper.user, time_keeper.system, time_keeper.total);    
    
    init_buff();
    PrintDbReport(print_data->dbinfo, 70, stdout);
    
    fflush(stdout);
    if (print_data->ka_params_gap) {
        PrintKAParameters(print_data->ka_params_gap->Lambda, 
                          print_data->ka_params_gap->K, 
                          print_data->ka_params_gap->H, 
                          70, stdout, TRUE);
    }
    fflush(stdout);
    
    PGPOutTextMessages(print_data->info_vnp, stdout);

    PrintTildeSepLines(print_data->buffer, 70, stdout);
    free_buff();

    fflush(stdout);

    ReadDBBioseqFetchDisable();

    return TRUE;
}

Boolean PSIPrintOutput(WWWBlastInfoPtr theInfo,
	BLASTPrintDataPtr print_data, 
	SeqAlignPtr BadSeqAlignments, SeqAlignPtr GoodSeqAlignments,
	SeqAlignPtr lastGood,
	Int2Ptr marks, Int2 countBad, Int2 countGood,
	Nlm_FloatHi ethresh_old)
{
    Uint4 align_options, print_options;
    SeqAnnotPtr seqannot;
    BlastTimeKeeper time_keeper;
    BlastPruneSapStructPtr prune;
    Uint1 f_order[FEATDEF_ANY], g_order[FEATDEF_ANY];
    Char	hostname[30], buffer[32];
    Char	href[1024];
    CharPtr	chptr;
    Char	f_name[64], title[1024];
    Int4	align_num;
    Int2	count;

    MemSet((Pointer)(g_order), 0, (size_t)(FEATDEF_ANY* sizeof(Uint1)));
    MemSet((Pointer)(f_order), 0, (size_t)(FEATDEF_ANY* sizeof(Uint1)));
    
    if(print_data == NULL) {
        WWWBlastErrMessage(BLASTMiscError, NULL);                
	return FALSE;
    }
    
    print_options = 0;
    align_options = 0;
    
    align_options += TXALIGN_COMPRESS;
    align_options += TXALIGN_END_NUM;
    
    if (theInfo->show_gi) {
	align_options += TXALIGN_SHOW_GI;
	print_options += TXALIGN_SHOW_GI;
    }
    
    if (theInfo->options->gapped_calculation == FALSE)
	print_options += TXALIGN_SHOW_NO_OF_SEGS;
    
    if (theInfo->align_view) {
	align_options += TXALIGN_MASTER;
        
	if (theInfo->align_view == 1 || theInfo->align_view == 3)
	    align_options += TXALIGN_MISMATCH;
        
	if (theInfo->align_view == 3 || theInfo->align_view == 4 || 
            theInfo->align_view == 6)
	    align_options += TXALIGN_FLAT_INS;

	if (theInfo->align_view == 5 || theInfo->align_view == 6)
	    align_options += TXALIGN_BLUNT_END;
    } else {
	align_options += TXALIGN_MATRIX_VAL;
	align_options += TXALIGN_SHOW_QS;
    }
    
    /* align_options += TXALIGN_MATRIX_VAL;
       align_options += TXALIGN_SHOW_QS; */

    align_options += TXALIGN_HTML;
    print_options += TXALIGN_HTML; 

    ReadDBBioseqFetchEnable ("psiblast", 
	    theInfo->database, 
	    theInfo->db_is_na, 
	    TRUE);

    seqannot = SeqAnnotNew();
    seqannot->type = 2;
    AddAlignInfoToSeqAnnot(seqannot, theInfo->align_type);
    seqannot->data = print_data->seqalign;

    init_buff();

    /* merge lists */
    if (lastGood)
        lastGood->next = BadSeqAlignments;
    
    /* gethostname(hostname, sizeof(hostname)); */
    
    sprintf(href, "%s/nph-viewgif.cgi?",
            theInfo->www_root_path == NULL? "/blast" : theInfo->www_root_path);
    
    if (theInfo->show_overview) {
        sprintf(f_name, "%ld%ld.gif", (long)random(), (long)getpid());
        align_num = get_number_alignment((SeqAlignPtr)(seqannot->data)); 
        sprintf(title, "<H3><a href=\"%s/docs/newoptions.html#graphical-overview\"> "
                "Distribution of %ld Blast Hits on the Query Sequence</a></H3>\n", theInfo->www_root_path == NULL? "/blast" : theInfo->www_root_path, (long)align_num);  
        
        PrintAlignmentOverview(seqannot, stdout, "PSI_BLAST", href, f_name, title); 
    }
    
    /* separate lists */
    if (lastGood)
        lastGood->next = NULL;
    
    print_options += TXALIGN_DO_NOT_PRINT_TITLE; 
    print_options += TXALIGN_CHECK_BOX;
    print_options += TXALIGN_CHECK_BOX_CHECKED;
    if (print_data->psidata->StepNumber)
        print_options += TXALIGN_NEW_GIF;
    
    /* submit button */
    printSubmitButton(stdout, 
                      print_data->psidata->StepNumber+1);
    
    if (print_data->psidata->StepNumber && theInfo->number_of_descriptions) {
        printf("<HR><p><b>Legend:</b><p>\
<IMG SRC=\"%s/images/new.gif\" WIDTH=25 HEIGHT=15> - means that \
the alignment score was below the threshold on the previous iteration \
<p>\
<IMG SRC=\"%s/images/checked.gif\" WIDTH=15 HEIGHT=15> - means that \
the alignment was checked on the previous iteration \
</p>", theInfo->www_root_path == NULL? "/blast" : theInfo->www_root_path,
               theInfo->www_root_path == NULL? "/blast" : theInfo->www_root_path);
    }
    
    if (theInfo->number_of_descriptions) {
        printf("<HR><CENTER><b><FONT color=\"green\">"
               "Sequences with E-value BETTER than threshold </FONT></b></CENTER>\n");
        if (print_data->psidata->StepNumber)
            printf("\n<IMG SRC=\"%s/images/bg.gif\" WIDTH=65 HEIGHT=15>", theInfo->www_root_path == NULL? "/blast" : theInfo->www_root_path);
        printf("                                                                     Score    E");
        printf("\nSequences producing significant alignments:");
        if (print_data->psidata->StepNumber)
            printf("<IMG SRC=\"%s/images/bg.gif\" WIDTH=65 HEIGHT=15>", theInfo->www_root_path == NULL? "/blast" : theInfo->www_root_path);
        printf("                          (bits) Value\n\n");
    }
    
    PrintDefLinesFromSeqAlignEx(GoodSeqAlignments, 80, stdout, 
                                print_options, FIRST_PASS, marks, theInfo->number_of_descriptions);
    
    print_options -= TXALIGN_CHECK_BOX_CHECKED;
    
    if (print_data->psidata->StepNumber == 0)
        printf("<a name = Evalue> </a>");
    
    if (theInfo->number_of_descriptions > countGood && BadSeqAlignments) {

        printSubmitButton(stdout, 
                          print_data->psidata->StepNumber+1);
        
        printf("<HR><CENTER><b><FONT color=\"green\">"
               "Sequences with E-value WORSE than threshold </FONT></b></CENTER>\n");
    
        PrintDefLinesFromSeqAlignEx(BadSeqAlignments, 80, stdout, print_options, FIRST_PASS, &marks[countGood], theInfo->number_of_descriptions - countGood);
    }

    if (theInfo->number_of_descriptions) {
	printSubmitButton(stdout, 
                          print_data->psidata->StepNumber+1);
    }

    free_buff();
    fflush(stdout);

    GlobalAlignNumber = 0;

    /* merge lists */
    if (lastGood)
	lastGood->next = BadSeqAlignments;

    prune = BlastPruneHitsFromSeqAlign((SeqAlignPtr) seqannot->data, theInfo->number_of_alignments, NULL);
    seqannot->data = prune->sap;

    fprintf(stdout, "<HR>");

    if (theInfo->number_of_alignments) {
	fprintf(stdout, "<CENTER><b><FONT color=\"green\">"
		"Alignments</FONT></b></CENTER>\n");

        /* New DDV formating requested */
        if(theInfo->color_schema != 0) {
            if(!DDV_DisplayBlastPairList(prune->sap, print_data->mask_loc, 
                                         stdout, 
                                         theInfo->query_is_na, align_options, 
                                         theInfo->color_schema)) { 
                fprintf(stdout, 
                        "\n\n!!!\n   "
                        "    --------  Failure to print alignment...  --------"
                        "\n!!!\n\n");
                fflush(stdout);
            }
        } else {   /* Old type formating */
            if (theInfo->align_view == 0) {
                ShowTextAlignFromAnnot2(seqannot, 60, stdout, f_order,
                                        g_order, align_options, NULL, 
                                        print_data->mask_loc, 
                                        FormatScoreFunc, theInfo->database, 
                                        "psiblast");
            } else { 
                ShowTextAlignFromAnnot(seqannot, 60, stdout, f_order,
                                       g_order, align_options, NULL, 
                                       print_data->mask_loc, NULL);
                printf("<P>\n");
            }
        }
    }

    /* separate lists */
    if (lastGood)
	lastGood->next = NULL;

    fflush(stdout);
    ObjMgrClearHold(); 

    prune = BlastPruneSapStructDestruct(prune);

    seqannot->data = NULL;
    seqannot = SeqAnnotFree(seqannot);

    printf("<PRE>\n");

    BlastTimeFillStructure(&time_keeper);

    fprintf(stdout, "CPU time: %8.2f user secs.\t%8.2f sys. "
	    "secs\t%8.2f total secs.\n\n", 
	    time_keeper.user, time_keeper.system, time_keeper.total);    

    init_buff();
    PrintDbReport(print_data->dbinfo, 70, stdout);
    
    if (print_data->ka_params) {
	PrintKAParameters(print_data->ka_params->Lambda, 
                          print_data->ka_params->K, 
                          print_data->ka_params->H, 70, 
                          stdout, FALSE);
    }
    
    if (print_data->ka_params_gap) {
	PrintKAParameters(print_data->ka_params_gap->Lambda, 
                          print_data->ka_params_gap->K, 
                          print_data->ka_params_gap->H, 
                          70, stdout, TRUE);
    }

    PrintTildeSepLines(print_data->buffer, 70, stdout);
    free_buff();
    
    fflush(stdout);

    ReadDBBioseqFetchDisable();

    return TRUE;
}