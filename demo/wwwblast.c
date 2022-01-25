/* $Id: wwwblast.c,v 6.4 2000/04/21 18:10:59 shavirin Exp $
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
* File Name:  $RCSfile: wwwblast.c,v $
*
* Author:  Sergei Shavirin
*
* Initial Creation Date: 03/15/2000
*
* $Revision: 6.4 $
*
* File Description:
*        Standalone WWW Blast CGI program.
*
* $Log: wwwblast.c,v $
* Revision 6.4  2000/04/21 18:10:59  shavirin
* Added possibility to print Patrick's alignment.
*
* Revision 6.3  2000/03/28 14:44:20  shavirin
* Changed function ctime_r to ctime() for compatibility.
*
* Revision 6.2  2000/03/24 16:05:37  shavirin
* Added option to be used as NCBI client/server.
*
* Revision 6.1  2000/03/20 19:01:00  shavirin
* Initial revision.
*
*
* ==========================================================================
*/

#include <ncbi.h>
#include <blastdef.h>
#include <blast.h>
#include <blastpri.h>
#include <sequtil.h>
#include <txalign.h>
#include <salogif.h>
#include <ddvcreate.h>

#ifdef NCBI_CLIENT_SERVER
#include <objblst3.h>
#include <netblap3.h>
#endif

typedef enum {
    BLASTNoError         =   0,   /* no error             */
    BLASTNetwork         =  -1,
    BLASTNoSpace         =  -2,
    BLASTBadFileName     =  -3,
    BLASTNotImplemented  =  -4,
    BLASTErrProgram      =  -5,  /* program missing from posting data    */
    BLASTErrDatalib      =  -6,  /* datalib missing from posting data    */
    BLASTErrProgName     =  -7,  /* program name is not supported        */
    BLASTErrNoSequence   =  -8,  /* NULL sequence passed to the engine   */
    BLASTErrCombination  =  -9,  /* Invalid program/database combination */
    BLASTNoMemory        =  -10, /* too bad ...                          */
    BLASTNCBI_DATA       =  -11, /* DATA files missing on any path       */
    BLASTFastaToSE       =  -12, /* FastaToSeqEntry() returned NULL      */
    BLASTErrOptions      =  -13, /* Badly formatted options              */
    BLASTErrNoQueue      =  -14, /* Queue overloaded */
    BLASTConfigFile      =  -15, /* Error reading config file */
    BLASTEntrez          =  -16, /* Cannot connect to Entrez */  
    BLASTAccesssion      =  -17, /* Invalid of unavailable accession */  
    BLASTSendmail        =  -18, /* Cannot start sendmail process */ 
    BLASTAddress         =  -19, /* Invalid return address */
    BLASTOptionStr       =  -20, /* Invalidly formatted advanced string */
    BLASTErrAccType      =  -21, /* wrong type of sequence identifier */
    BLASTErrClient       =  -22, /* cannot connect to the Blast service */
    BLASTErrServer       =  -23, /* Error from the server side */
    BLASTMiscError       =  -99  /* undefined internal error             */
} BLASTErrCode;

#define MAX_DB_NUM 256
#define DEFAULT_CONFIG_FILE "blast.rc"

/* Max. total number of concurrent running requests */
#define DEFAULT_RUN_MAX            2 

/* Max. total number of waiting requests */
#define DEFAULT_QUEUE_MAX          100 
#define NUM_CPU_TO_USE             4
#define DEFAULT_DESCRIPTIONS    100
#define DEFAULT_ALIGNMENTS       50 
#define DEFAULT_EXPECT           10

/* CPU time limit. */
#define DEFAULT_CPU_LIMIT 3600

typedef struct BLASTConfig {
    Int4 run_max;
    Int4 queue_max;
    Int4 num_cpu;
    Int4 niceval;
    CharPtr allow_db[MAX_DB_NUM];
} BLASTConfig, PNTR BLASTConfigPtr;

typedef struct _www_blast_info {
    BLAST_OptionsBlkPtr options;
    WWWInfoPtr info;
    BLASTErrCode error_code;
    CharPtr ConfigFile;
    CharPtr program, database, blast_type;
    BioseqPtr query_bsp;
    BioseqPtr fake_bsp;
    Int4 number_of_descriptions, number_of_alignments;
    Boolean query_is_na, db_is_na, align_type, show_gi, show_overview;
    Boolean believe_query;
    Uint4 align_options, print_options;
    Int4 align_view, input_type, color_schema;
    BLASTConfigPtr blast_config;
} WWWBlastInfo, PNTR WWWBlastInfoPtr;

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
void BLASTConfigGetWord(CharPtr word, CharPtr line) 
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
BLASTConfigPtr BLASTReadConfigFile(CharPtr filename, CharPtr program)
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

static Boolean BLAST_Time(CharPtr string, Int4 len, time_t seconds)
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

WWWBlastInfoPtr WWWBlastReadArgs(void)
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
                "CONTENT=\"2; URL=blast.html\">");
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

    if ( !ErrSetLogfile ("/dev/null", ELOG_APPEND) ) {
        fprintf(stdout, "Cannot set logfile exiting....\n");
        return FALSE;
    } else {
        ErrSetOpts (ERR_CONTINUE, ERR_LOG_ON);
    }

    /* Config file with program/database relationsship */
    
    blast_type = WWWGetValueByName(theInfo->info, "WWW_BLAST_TYPE");
    
    if(blast_type == NULL || *blast_type == NULLB) {
        theInfo->ConfigFile = StringSave(DEFAULT_CONFIG_FILE);
        theInfo->blast_type = StringSave("blast");
    } else {
        sprintf(tmp, "%s.rc", blast_type);
        theInfo->ConfigFile = StringSave(tmp);
        theInfo->blast_type = StringSave(blast_type);
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
Boolean ValidateCombinationsEx(WWWBlastInfoPtr theInfo, CharPtr database)
{
    Int4 i;
    
    for(i = 0; theInfo->blast_config->allow_db[i] != NULL; i++) {
	if(!StringICmp(database, theInfo->blast_config->allow_db[i]))
	    return TRUE;
    }
    return FALSE;
}

/* This will work if search require to choose few databases */
Boolean WWWParseDatabases(WWWBlastInfoPtr theInfo)
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
        Int4 value;
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
        Int4 value;
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
    
    if((chptr = WWWGetValueByName(theInfo->info, 
                                  "DESCRIPTIONS")) != NULL && 
       StringStr(chptr, "default") == NULL) {
	theInfo->number_of_descriptions = atoi(chptr);
        /* if(theInfo->NumDescr == 0)
           theInfo->NumDescr = DEFAULT_DESCRIPTIONS; */
    }
    
    /* ALIGNMENTS */
    
    if((chptr = WWWGetValueByName(theInfo->info, "ALIGNMENTS")) != NULL &&
       StringStr(chptr, "default") == NULL) {
	theInfo->number_of_alignments = atoi(chptr);  
        /* if(theInfo->NumAlign == 0)
           theInfo->NumAlign = DEFAULT_ALIGNMENTS; */
    }
    
    /* Now processing OTHER_ADVANCED_OPTIONS */
    
    if((chptr = WWWGetValueByName(theInfo->info, 
                                  "OTHER_ADVANCED")) != NULL) {
        Int4 index;
        CharPtr ErrorMessage;
        CharPtr PNTR values;

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


    }
        /*
          options->hsp_range_max = 100*options->hitlist_size;
          options->block_width = theInfo->bsp->length;
          */

    options->perform_culling = FALSE;

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

    if ((chptr = WWWGetValueByName(theInfo->info, "COLOR_SCHEMA")) != NULL)
	theInfo->color_schema = atoi(chptr);    

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
static int LIBCALLBACK
tick_callback(Int4 sequence_number, Int4 number_of_positive_hits)

{

    /* do nothing */
    return 0;
}

static Int4 get_number_alignment(SeqAlignPtr align)
{
    Int4 num = 0;

    while(align)
    {
	++num;
	align = align->next;
    }

    return num;
}
static void
PrintMotd(CharPtr string, FILE *fp, Boolean html_format)

{
    Char buffer[100];
    CharPtr ptr;
    
    if (string == NULL)
        return;
    
    buffer[0] = NULLB;
    ptr = buffer;
    
    if (html_format) {
        fprintf(fp, "<PRE>\n");
    }
    
    while (*string != NULLB) {
        if (*string == '~') {
            *ptr = NULLB;
            fprintf(fp, "%s\n", buffer);
            buffer[0] = NULLB;
            ptr = buffer;
            string++;
            if (*string == NULLB)
                break;
        } else {
            *ptr=*string;
            ptr++;  string++;
        }
    }
    *ptr = NULLB;
    fprintf(fp, "%s\n", buffer);
    
    if (html_format) {
        fprintf(fp, "</PRE>\n");
    }
    
    fflush(fp);

    return;
}
#ifdef NCBI_CLIENT_SERVER
static Boolean
TraditionalBlastReportWithImage(BioseqPtr bsp, BLAST_OptionsBlkPtr options, BlastNet3Hptr bl3hp, CharPtr program, CharPtr database, Uint4 print_options, Uint4 align_options, Int4 number_of_descriptions, Int4 number_of_alignments, Boolean show_overview, CharPtr blast_type, Int4 color_schema)

{
    BlastDbinfoPtr dbinfo;
    BlastKABlkPtr ka_params=NULL, ka_params_gap=NULL;
    BlastPruneSapStructPtr prune;
    BLAST_MatrixPtr matrix;
    Boolean query_is_na, db_is_na;
    Boolean status;
    CharPtr params_buffer=NULL;
    Int4 number_of_hits_private=0, length;
    SeqAlignPtr seqalign;
    SeqAnnotPtr seqannot=NULL;
    TxDfDbInfoPtr tx_dbinfo=NULL, tx_dbinfo_head;
    ValNodePtr mask_loc, mask_loc_start, other_returns, error_returns, vnp, vnp1=NULL;
    Uint1 align_type;
    Uint1 f_order[FEATDEF_ANY], g_order[FEATDEF_ANY];
    
    MemSet((Pointer)(g_order), 0, (size_t)(FEATDEF_ANY* sizeof(Uint1)));
    MemSet((Pointer)(f_order), 0, (size_t)(FEATDEF_ANY* sizeof(Uint1)));
    
    if (bsp == NULL)
        return FALSE;
    
    if (bl3hp == NULL || program == NULL || database == NULL)
        return FALSE;
    
    align_type = BlastGetTypes(program, &query_is_na, &db_is_na);

    init_buff_ex(85);
    dbinfo = BlastRequestDbInfo(bl3hp, database, !db_is_na);

    if (dbinfo)
        PrintDbInformationBasic(database, !db_is_na, 70, dbinfo->definition, dbinfo->number_seqs, dbinfo->total_length, stdout, TRUE);
    dbinfo = BlastDbinfoFree(dbinfo);
    free_buff();

    if (bsp) {
        seqalign = BlastBioseqNetCore(bl3hp, bsp, program, database, options, &other_returns, &error_returns, NULL, NULL, &status);
    }
    
    BlastErrorPrintExtra(error_returns, TRUE, stdout);
    
    mask_loc = NULL;
    for (vnp=other_returns; vnp; vnp = vnp->next) {
        switch (vnp->choice) {
        case TXDBINFO:
            tx_dbinfo = (TxDfDbInfoPtr) vnp->data.ptrvalue;
            break;
        case TXKABLK_NOGAP:
            ka_params = (BlastKABlkPtr) vnp->data.ptrvalue;
            break;
        case TXKABLK_GAP:
            ka_params_gap = (BlastKABlkPtr) vnp->data.ptrvalue;
            break;
        case TXPARAMETERS:
            params_buffer = (CharPtr) vnp->data.ptrvalue;
            break;
        case TXMATRIX:
            /* This function should not use matrix */
            matrix = (BLAST_MatrixPtr) vnp->data.ptrvalue;
            BLAST_MatrixDestruct(matrix);
            break;
        case SEQLOC_MASKING_NOTSET:
        case SEQLOC_MASKING_PLUS1:
        case SEQLOC_MASKING_PLUS2:
        case SEQLOC_MASKING_PLUS3:
        case SEQLOC_MASKING_MINUS1:
        case SEQLOC_MASKING_MINUS2:
        case SEQLOC_MASKING_MINUS3:
            ValNodeAddPointer(&mask_loc, vnp->choice, vnp->data.ptrvalue);
            break;
        default:
            break;
        }
    }	
    
    if (seqalign) {
        seqannot = SeqAnnotNew();
        seqannot->type = 2;
        AddAlignInfoToSeqAnnot(seqannot, align_type);
        seqannot->data = seqalign;
        
        if(show_overview) {
            Char f_name[64], title[1024], href[64];
            Int4 align_num;       

            sprintf(f_name, "%ld%ld.gif", (long)random(), (long)getpid());
            sprintf(href, "nph-viewgif.cgi?");
        
            align_num = get_number_alignment(seqalign); 
            sprintf(title, 
                    "<H3><a href=\"docs/newoptions.html#graphical-overview\"> "
                    "Distribution of %ld Blast Hits on the Query Sequence</a> "
                    "</H3>\n", (long)align_num);  
            
            
            /* Open HTML form */
            fprintf(stdout, "<FORM NAME=\"BLASTFORM\">\n");
            fflush(stdout);
            
            PrintAlignmentOverview(seqannot, stdout, 
                                   "BLASTFORM", href, f_name, title); 
        }

        prune = BlastPruneHitsFromSeqAlign(seqalign, number_of_descriptions, NULL);
        ObjMgrSetHold();
        init_buff_ex(85);

        PrintDefLinesFromSeqAlignEx2(prune->sap, 80, stdout, print_options, FIRST_PASS, NULL, number_of_descriptions, database, blast_type);
        free_buff();
        
        prune = BlastPruneHitsFromSeqAlign(seqalign, number_of_alignments, prune);
        seqannot->data = prune->sap;

#ifdef NEW_ALIGNMENT

        if(!DDV_DisplayBlastPairList(prune->sap, mask_loc, stdout, 
                                     query_is_na, align_options, 
                                     color_schema)) { 
            fprintf(stdout, 
                    "\n\n!!!\n   "
                    "    --------  Failure to print alignment...  --------"
                    "\n!!!\n\n");
            fflush(stdout);
        }
#else

        if (align_options & TXALIGN_MASTER)
            ShowTextAlignFromAnnot(seqannot, 60, stdout, f_order, g_order, align_options, NULL, mask_loc, NULL);
        else
            ShowTextAlignFromAnnot(seqannot, 60, stdout, f_order, g_order, align_options, NULL, mask_loc, FormatScoreFunc);
#endif
        printf("<PRE>\n");

        seqannot->data = seqalign;
        number_of_hits_private = prune->original_number; 
        prune = BlastPruneSapStructDestruct(prune);
        ObjMgrClearHold();
    }
    
    mask_loc_start = mask_loc;
    while (mask_loc) {
        SeqLocSetFree(mask_loc->data.ptrvalue);
        mask_loc = mask_loc->next;
    }
    ValNodeFree(mask_loc_start);
    
    init_buff_ex(85);
    tx_dbinfo_head = tx_dbinfo;
    while (tx_dbinfo) {
        PrintDbReport(tx_dbinfo, 70, stdout);
        tx_dbinfo = tx_dbinfo->next;
    }
    tx_dbinfo_head = TxDfDbInfoDestruct(tx_dbinfo_head);
    
    if (ka_params) {
        PrintKAParameters(ka_params->lambda, ka_params->k, ka_params->h, 
                          70, stdout, FALSE);
        MemFree(ka_params);
    }
    
    if (ka_params_gap) {
        PrintKAParameters(ka_params_gap->lambda, ka_params_gap->k, ka_params_gap->h, 70, stdout, TRUE);
        MemFree(ka_params_gap);
    }
    
    PrintTildeSepLines(params_buffer, 70, stdout);
    MemFree(params_buffer);
    free_buff();
    
    other_returns = ValNodeFree(other_returns);
    
    if (seqannot)
        seqannot = SeqAnnotFree(seqannot);
    
    return status;
}

Boolean WWWBlastDoClientSearch(WWWBlastInfoPtr theInfo)
{
    BlastNet3Hptr	bl3hp;
    BlastResponsePtr	response;
    BioseqPtr fake_bsp; /* Needed for correct formating */
    BlastVersionPtr	blast_version;
    CharPtr		date, motd, version;
    Boolean status;
    
    if(theInfo == NULL)
	return FALSE;
    
    /* This will prevent from incorrect formating in case when input
       sequence has valig SeqId, but in fact this SeqId do not correspond
       to the real sequence  - XXX */
    
    if(!theInfo->believe_query)
        fake_bsp = BlastMakeFakeBioseq(theInfo->query_bsp, NULL);
    else
        fake_bsp = theInfo->query_bsp;
    
    
    if (!BlastInit("blastcl3", &bl3hp, &response)) {
        WWWBlastErrMessage(BLASTErrClient, NULL);                
        return FALSE;
    }
    
    if (response && response->choice == BlastResponse_init) {
        blast_version = response->data.ptrvalue;
        version = blast_version->version;
        date = blast_version->date;
    } else {
        WWWBlastErrMessage(BLASTErrClient, NULL);                
        return FALSE;
    }
    
    BlastNetBioseqFetchEnable(bl3hp, theInfo->database, 
                              theInfo->db_is_na, TRUE);
    
#ifdef BLAST_PRINT_MOTD
    motd = Blast3GetMotd(bl3hp);
    PrintMotd(motd, stdout, TRUE);
    motd = MemFree(motd);
#endif
    
    status = TraditionalBlastReportWithImage(fake_bsp, theInfo->options, 
                                    bl3hp, theInfo->program, 
                                    theInfo->database,
                                    theInfo->print_options, 
                                    theInfo->align_options, 
                                    theInfo->number_of_descriptions, 
                                    theInfo->number_of_alignments, 
                                    theInfo->show_overview,
                                    theInfo->blast_type,
                                    theInfo->color_schema);
    if (status == FALSE) {
        WWWBlastErrMessage(BLASTErrServer, NULL);
        return FALSE;
    }
    
    return TRUE;
}
#endif

Boolean WWWBlastDoSearch(WWWBlastInfoPtr theInfo)
{
    BioseqPtr fake_bsp; /* Needed for correct formating */
    SeqAlignPtr  seqalign;
    ValNodePtr  mask_loc, mask_loc_start, vnp, other_returns, error_returns;
    TxDfDbInfoPtr dbinfo=NULL, dbinfo_head;
    BLAST_KarlinBlkPtr ka_params=NULL, ka_params_gap=NULL;
    BLAST_MatrixPtr matrix;
    CharPtr params_buffer=NULL;
    SeqAnnotPtr seqannot;
    BlastPruneSapStructPtr prune;
    
    if(theInfo == NULL)
	return FALSE;
    
    PrintDbInformation(theInfo->database, !theInfo->db_is_na, 
                       70, stdout, TRUE);

    /* This will prevent from incorrect formating in case when input
       sequence has valig SeqId, but in fact this SeqId do not correspond
       to the real sequence  - XXX */
    
    if(!theInfo->believe_query)
        fake_bsp = BlastMakeFakeBioseq(theInfo->query_bsp, NULL);
    else
        fake_bsp = theInfo->query_bsp;
    
    other_returns = NULL;
    error_returns = NULL;
    
    seqalign = BioseqBlastEngine(fake_bsp, theInfo->program, 
                                 theInfo->database, 
                                 theInfo->options, 
                                 &other_returns, 
                                 &error_returns, 
                                 tick_callback);
    
    BlastErrorPrint(error_returns);
    
    dbinfo = NULL;
    ka_params = NULL;
    ka_params_gap = NULL;
    params_buffer = NULL;
    mask_loc = NULL;
    matrix = NULL;
    for (vnp=other_returns; vnp; vnp = vnp->next) {
        switch (vnp->choice) {
        case TXDBINFO:
            dbinfo = vnp->data.ptrvalue;
            break;
        case TXKABLK_NOGAP:
            ka_params = vnp->data.ptrvalue;
            break;
        case TXKABLK_GAP:
            ka_params_gap = vnp->data.ptrvalue;
            break;
        case TXPARAMETERS:
            params_buffer = vnp->data.ptrvalue;
            break;
        case TXMATRIX:
            matrix = vnp->data.ptrvalue;
            break;
        case SEQLOC_MASKING_NOTSET:
        case SEQLOC_MASKING_PLUS1:
        case SEQLOC_MASKING_PLUS2:
        case SEQLOC_MASKING_PLUS3:
        case SEQLOC_MASKING_MINUS1:
        case SEQLOC_MASKING_MINUS2:
        case SEQLOC_MASKING_MINUS3:
            ValNodeAddPointer(&mask_loc, vnp->choice, vnp->data.ptrvalue);
            break;
        default:
            break;
        }
    }	
    
    
    fflush(stdout);
    
    ReadDBBioseqFetchEnable ("blastall", theInfo->database, 
                             theInfo->db_is_na, TRUE);
    
    if (seqalign) {
        seqannot = SeqAnnotNew();
        seqannot->type = 2;
        AddAlignInfoToSeqAnnot(seqannot, theInfo->align_type);
        seqannot->data = seqalign;

        /* Now printing nice gif with alignment overview */

        if(theInfo->show_overview) {
            Char f_name[64], title[1024], href[64];
            Int4 align_num;
            
            sprintf(f_name, "%ld%ld.gif", (long)random(), (long)getpid());
            sprintf(href, "nph-viewgif.cgi?");
            
            align_num = get_number_alignment(seqalign); 
            sprintf(title, 
                    "<H3><a href=\"docs/newoptions.html#graphical-overview\"> "
                    "Distribution of %ld Blast Hits on the Query Sequence</a> "
                    "</H3>\n", (long)align_num);  
	

            /* Open HTML form */
            fprintf(stdout, "<FORM NAME=\"BLASTFORM\">\n");
            fflush(stdout);
            
            PrintAlignmentOverview(seqannot, stdout, 
                                   "BLASTFORM", href, f_name, title); 
        }

        prune = BlastPruneHitsFromSeqAlign(seqalign, theInfo->number_of_descriptions, NULL);
        ObjMgrSetHold();
        init_buff_ex(85);
        PrintDefLinesFromSeqAlign(prune->sap, 80, stdout, 
                                  theInfo->print_options, FIRST_PASS, NULL);
        free_buff();
        
        prune = BlastPruneHitsFromSeqAlign(seqalign, theInfo->number_of_alignments, prune);
        seqannot->data = prune->sap;

#ifdef NEW_ALIGNMENT

        if(!DDV_DisplayBlastPairList(prune->sap, mask_loc, stdout, 
                                     theInfo->query_is_na, 
                                     theInfo->align_options, 
                                     theInfo->color_schema)) { 
            fprintf(stdout, 
                    "\n\n!!!\n   "
                    "    --------  Failure to print alignment...  --------"
                    "\n!!!\n\n");
            fflush(stdout);
        }
        
        printf("<PRE>\n");
#else
        if (theInfo->align_view != 0)
            ShowTextAlignFromAnnot(seqannot, 60, stdout, NULL, NULL, 
                                   theInfo->align_options, NULL, 
                                   mask_loc, NULL);
        else
            ShowTextAlignFromAnnot(seqannot, 60, stdout, NULL, NULL, 
                                   theInfo->align_options, NULL, mask_loc, 
                                   FormatScoreFunc);
#endif

        seqannot->data = seqalign;
        prune = BlastPruneSapStructDestruct(prune);
        ObjMgrClearHold();
        ObjMgrFreeCache(0);
        
    } else {
        fprintf(stdout, "\n\n ***** No hits found ******\n\n");
    }
    
    matrix = BLAST_MatrixDestruct(matrix);
    
    fprintf(stdout, "<PRE>");
    init_buff_ex(85);
    dbinfo_head = dbinfo;
    while (dbinfo) {
        PrintDbReport(dbinfo, 70, stdout);
        dbinfo = dbinfo->next;
    }
    dbinfo_head = TxDfDbInfoDestruct(dbinfo_head);
    
    if (ka_params) {
        PrintKAParameters(ka_params->Lambda, ka_params->K, ka_params->H, 70, stdout, FALSE);
        MemFree(ka_params);
    }
    
    if (ka_params_gap) {
        PrintKAParameters(ka_params_gap->Lambda, ka_params_gap->K, ka_params_gap->H, 70, stdout, TRUE);
        MemFree(ka_params_gap);
    }
    
    PrintTildeSepLines(params_buffer, 70, stdout);

    fprintf(stdout, "</PRE>\n</BODY>\n</HTML>\n");

    seqannot = SeqAnnotFree(seqannot);

    MemFree(params_buffer);
    free_buff();
    
    mask_loc_start = mask_loc;
    while (mask_loc) {
        SeqLocSetFree(mask_loc->data.ptrvalue);
        mask_loc = mask_loc->next;
    }
    ValNodeFree(mask_loc_start);
    
    if(!theInfo->believe_query)
        fake_bsp = BlastDeleteFakeBioseq(fake_bsp);
    
    ReadDBBioseqFetchDisable();
    other_returns = ValNodeFree(other_returns);

    return TRUE;
}

void WWWBlastPrintHeader(WWWBlastInfoPtr theInfo)
{

    fprintf(stdout, "<HTML>\n");
    fprintf(stdout, "<TITLE>BLAST Search Results</TITLE>\n"); 
    fflush(stdout);
    
    fprintf(stdout, "<BODY BGCOLOR=\"#FFFFFF\" LINK=\"#0000FF\" "
            "VLINK=\"#660099\" ALINK=\"#660099\">\n");
    fprintf(stdout, "<A HREF=\"blast_form.map\"> \r"
            "<IMG SRC=\"images/blast_results.gif\" "
            "BORDER=0 ISMAP>\r</A><P>\n");    
    fprintf(stdout, "<PRE>\n");
    init_buff_ex(90);
    BlastPrintVersionInfo(theInfo->program, TRUE, stdout);
    fprintf(stdout, "\n");
    BlastPrintReference(TRUE, 90, stdout);
    fprintf(stdout, "\n");
    AcknowledgeBlastQuery(theInfo->query_bsp, 70, stdout, 
                          theInfo->believe_query, TRUE);
    free_buff();

    return;
}

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

    MemFree(theInfo->ConfigFile);
    MemFree(theInfo);

    return;
}
Int2 Main(void)
{
    WWWBlastInfoPtr theInfo;
    
    UseLocalAsnloadDataAndErrMsg ();

    /* This function will read posting data, set-up config file and
       write small message into logfile (if it exists) */
    
    if((theInfo = WWWBlastReadArgs()) == NULL)
        return 1;
    
    /* Read options into structure */
    if(!WWWCreateSearchOptions(theInfo)) {
        return 1;
    }
    
    /* validate them */
    if(!WWWValidateOptions(theInfo)) {
        return 1;
    }

    /* Print BLAST Header */

    WWWBlastPrintHeader(theInfo);

    /* Do the search and Format output */

#ifdef NCBI_CLIENT_SERVER
    WWWBlastDoClientSearch(theInfo);
#else
    WWWBlastDoSearch(theInfo);
#endif

    WWWBlastInfoFree(theInfo);

    return 0;
}
