/*****************************************************************************

  
                          PUBLIC DOMAIN NOTICE
              National Center for Biotechnology Information

    This software/database is a "United States Government Work" under the
    terms of the United States Copyright Act.  It was written as part of    
    the author's official duties as a United States Government employee
    and thus cannot be copyrighted.  This software/database is freely
    available to the public for use. The National Library of Medicine and
    the U.S. Government have not placed any restriction on its use or
    reproduction.

    Although all reasonable efforts have been taken to ensure the accuracy
    and reliability of the software and data, the NLM and the U.S.
    Government do not and cannot warrant the performance or results that
    may be obtained by using this software or data. The NLM and the U.S.
    Government disclaim all warranties, express or implied, including
    warranties of performance, merchantability or fitness for any
    particular purpose.

    Please cite the author in any work or product based on this material.

   ***************************************************************************

   File Name:  formatdb.c

   Author:  Sergei B. Shavirin
   
   Version Creation Date: 10/01/96

   $Revision: 6.13 $

   File Description:  formats FASTA databases for use by BLAST

   $Log: formatdb.c,v $
   Revision 6.13  1998/11/16 18:34:42  madden
   Add return-value checks

   Revision 6.12  1998/07/13 15:32:17  egorov
   make error message more understandable

   Revision 6.10  1998/06/19 21:05:46  egorov
   Fix MemFree() bug

   Revision 6.9  1998/05/05 13:57:37  madden
   Print version number to log file

   Revision 6.8  1998/04/20 19:14:05  egorov
   Fix just one, but huge MLK

   Revision 6.7  1998/02/23 16:49:14  egorov
   Changes to make the tofasta.c independent on readdb.h

   Revision 6.6  1998/02/18 15:29:31  madden
   Added const to prototype for FormatdbCreateStringIndex

   Revision 6.5  1998/02/11 18:05:32  madden
   Changed program to take ASN.1 as input

   Revision 6.3  1997/12/08 21:55:00  madden
   Parse naked (no bars) as IDs

   Revision 6.2  1997/11/06 18:11:17  madden
   Added indices for naked gnl|PID and backbone entries

   Revision 6.1  1997/10/30 18:15:08  madden
   Changes to SeqIdE2Index to allow lookups by accession strings

   Revision 6.0  1997/08/25 18:20:04  madden
   Revision changed to 6.0

   Revision 1.20  1997/07/28 18:36:55  madden
   Replaced printf with ErrPostEx and fprintf

   Revision 1.19  1997/07/28 14:35:37  vakatov
   Added LIBCALLBACK to the ID_Compare() proto

   Revision 1.18  1997/06/10 18:44:11  shavirin
   Fixed return value from UpdateLookupInfo()

   Revision 1.17  1997/05/19 21:16:30  shavirin
   Changed content of string index file due to E2Iindex API logic

   Revision 1.16  1997/05/12 19:57:38  shavirin
   Added additional dump of Accessions/Locuses into string indexes

   Revision 1.15  1997/05/07 21:08:15  madden
   flipped parse argument default

   Revision 1.14  1997/05/05 17:01:42  shavirin
   Added ability to format "non-parced" seqid-deflines
   Removed not-used d if#defs  with FASTA_ASN

 * Revision 1.13  1997/05/01  17:31:32  shavirin
 * Added dumping of 2 more files: String ISAM SeqId index
 *
 * Revision 1.12  1997/02/25  22:20:39  shavirin
 * Changes in accordance to ISAM API changes
 *
 * Revision 1.11  1997/02/24  21:22:57  shavirin
 * Added dump of numeric ISAM information.
 *
 * Revision 1.10  1996/12/20  00:31:19  madden
 * Protected ambiguity data against big/little endian changes.
 *
 * Revision 1.9  1996/12/19  16:30:36  madden
 * Changes to eliminate ".nac" file for nucl.
 *
 * Revision 1.8  1996/11/27  16:40:19  madden
 * Save build date, Make "o" argument FALSE by default.
 *
 * Revision 1.7  1996/11/26  20:08:08  madden
 * BioseqRawConvert(bsp, Seq_code_ncbistdaa); only called for protein alphabets.
 *
 * Revision 1.6  1996/11/26  19:52:10  madden
 * Removed FORMATDB_VER and added readdb.h (which contains same);
 * Changed phd or nhd to phr or nhr
 *
 * Revision 1.5  1996/11/18  20:53:58  shavirin
 * Forced output protein code to Seq_code_ncbistdaa.
 *
 * Revision 1.4  1996/11/06  23:15:34  shavirin
 * Removed bug with reallocation of index tables
 *

*****************************************************************************/
#include <ncbi.h>
#include <tofasta.h>
#include <sequtil.h>
#include <readdb.h>
#include <ncbisam.h>
#include <ncbisort.h>
#include <blast.h>

#define STRLENGTH     4096
#define INDEX_ARRAY_CHUNKS 100000

#define LOOKUP_CHUNK   5
#define LOOKUP_SIZE    12
#define LOOKUP_ID_SIZE 8

#define FORMATDB_SIZE 4
#define ID_MAX_SIZE   64

#define LOOKUP_NO_ERROR  0
#define ERR_GI_FAILED    1
#define ERR_SEQID_FAILED 2

#define NON_SEQID_PREFIX "gnl|BL_ORD_ID|"
#define CREATE_DEFLINE_INDEX 1

#define SEQID_FIELD   1
#define ACCN_FIELD    2
#define DEFLINE_FIELD 4



    /* static functions */

static Boolean		FormatDbUint4Write(Uint4 number, FILE *fp);
static Int4		UpdateLookupInfo(CharPtr defline, FASTALookupPtr lookup, 
                                         Int4 num_of_seqs, FILE *fd_stmp,
                                         Boolean ParseSeqid);
static int LIBCALLBACK	ID_Compare(VoidPtr i, VoidPtr j);
static FASTALookupPtr	FASTALookupNew(void);
static void		FASTALookupFree(FASTALookupPtr lookup);
static Boolean		FormatdbCreateStringIndex(const CharPtr FileName, 
                                         Boolean ProteinDump);
static Boolean	SeqIdE2Index (SeqIdPtr anp, FILE *fd, Int4 seq_num);

static	Int2	process_sep (SeqEntryPtr sep, FormatDBPtr fdbp);
static	Int2	finish_formatDB (FormatDBPtr fdbp);


    /* program's arguments */

#define NUMARG 8

Args dump_args[NUMARG] = {
    { "Title for database file", 
      NULL, NULL, NULL, TRUE, 't', ARG_STRING, 0.0, 0, NULL},
    {"Input file for formatting (this parameter must be set)",
     NULL, NULL,NULL,FALSE,'i',ARG_FILE_IN, 0.0,0,NULL},
    {"Logfile name:",
     "formatdb.log", NULL,NULL,TRUE,'l',ARG_FILE_OUT, 0.0,0,NULL},
    {"Type of file\n"
     "         T - protein   \n"
     "         F - nucleotide", 
     "T", NULL,NULL,TRUE,'p',ARG_BOOLEAN,0.0,0,NULL},
    {"Parse options\n"
     "         T - True: Parse SeqId and create indexes.\n"
     "         F - False: Do not parse SeqId. Do not create indexes.\n",
     "F", NULL,NULL,TRUE,'o',ARG_BOOLEAN,0.0,0,NULL},
    {"Input file is database in ASN.1 format (otherwise FASTA is expected)\n"
     "         T - True, \n"
     "         F - False.\n",
     "F", NULL,NULL,TRUE,'a',ARG_BOOLEAN,0.0,0,NULL},
    {"ASN.1 database in binary mode\n"
     "         T - binary, \n"
     "         F - text mode.\n",
     "F", NULL,NULL,TRUE,'b',ARG_BOOLEAN,0.0,0,NULL},
    {"Input is a Seq-entry","F", NULL ,NULL ,TRUE,'e',ARG_BOOLEAN,0.0,0,NULL}
};

#define db_title	(const CharPtr) dump_args[0].strvalue 
#define db_file		(const CharPtr) dump_args[1].strvalue 
#define LogFileName	(const CharPtr) dump_args[2].strvalue  
#define is_protein	dump_args[3].intvalue
#define parse_mode	dump_args[4].intvalue
#define isASN		dump_args[5].intvalue
#define asnbin		dump_args[6].intvalue
#define is_seqentry	dump_args[7].intvalue


Int4 OffsetAllocated = INDEX_ARRAY_CHUNKS;

    /* main() */

Int2 Main(void) 
{
    SeqEntryPtr sep;
    FormatDBPtr	fdbp;
    Uint1 group_segs = 0;

        /* get arguments */

    if ( !GetArgs ("formatdb", NUMARG, dump_args) ) {
        return -1;
    }
    if ( !ErrSetLog (LogFileName) ) {
        ErrShow();
    } else {
        ErrSetOpts (ERR_CONTINUE, ERR_LOG_ON);
    }

        /* Initialize formatdb structure */
    if ((fdbp =
         FormatDBInit(db_file, db_title, is_protein, parse_mode, isASN, asnbin))
        == NULL)
        return -1;
     
    
        /* Input database file maybe either in ASN.1 or in FASTA format */
    
    if (!isASN) {
            /* FASTA format of input database */

            /* Get sequences */
        while ((sep = FastaToSeqEntryEx(fdbp->fd, (Boolean)!is_protein,
                                        NULL, parse_mode)) != NULL) {

            if(!IS_Bioseq(sep)) { /* Not Bioseq - failure */
                ErrLogPrintf("Error in readind Bioseq Formating failed.\n");
                return -1;
            }
            
            if (process_sep(sep, fdbp))
                return -1;
                
            SeqEntryFree(sep);
            fdbp->num_of_seqs++;  
        }
    }
    else {
            /* ASN.1 format of input database */
        AsnTypePtr atp, atp2;
        AsnModulePtr amp;

        if (! SeqEntryLoad())
            ErrShow();

            /* get pointer to all loaded ASN.1 modules */
        
            /* get pointer to all loaded ASN.1 modules */
        amp = AsnAllModPtr();
        if (amp == NULL)
        {
            ErrLogPrintf("Could not load ASN.1 modules.\n");
            return -1;
        }
        
            /* get the initial type pointers */
        
        atp = AsnFind("Bioseq-set");
        if (atp == NULL)
        {
            ErrLogPrintf("Could not get type pointer for Bioseq-set.\n");
            return -1;
        }

        atp2 = AsnFind("Bioseq-set.seq-set.E");
        if (atp2 == NULL)
        {
            ErrLogPrintf("Could not get type pointer for Bioseq-set.seq-set.E\n");
            return -1;
        }

        if ((fdbp->aip = AsnIoOpen (db_file, asnbin ? "rb":"r")) == NULL)
        {
            ErrLogPrintf("Cannot open input database file. Formating failed...\n");
            return -1;
        }

        if (is_seqentry) {
                /* Seq entry */
            sep = SeqEntryAsnRead(fdbp->aip, NULL);
            SeqEntrysToBLAST(sep, fdbp, !(is_protein), group_segs);
            SeqEntryFree(sep);
        }
        else {
            /* Bioseq-set */
        
            while ((atp = AsnReadId(fdbp->aip, amp, atp)) != NULL)
            {
                if (atp == atp2)    /* top level Seq-entry */
                {
                    sep = SeqEntryAsnRead(fdbp->aip, atp);
                    SeqEntrysToBLAST(sep, fdbp, !(is_protein), group_segs);
                    SeqEntryFree(sep);
                }
                else
                {
                    AsnReadVal(fdbp->aip, atp, NULL);
                }
            }
        } /* end "if Bioseq or Bioseq-set */
        
    } /* end "if FASTA or ASN.1" */
    
        /* write files */
    if (finish_formatDB(fdbp))
        return -1;
    
        /* Deallocate structure, arrays, etc. */
    fdbp = FormatDBDestruct(fdbp);

    return 0;
    
} /* main()*/


NLM_EXTERN Boolean SeqEntrysToBLAST (SeqEntryPtr sep, FormatDBPtr fdbp,
                                     Boolean is_na, Uint1 group_segs)
{
    FastaDat tfa;
    MyFsa mfa;
    Char buf[255];
    
    if ((sep == NULL) || (fdbp == NULL))
        return FALSE;
    
    mfa.buf	= buf;
    mfa.buflen	= 254;
    mfa.seqlen	= 70;
    mfa.mydata	= (Pointer)fdbp;
    mfa.myfunc	= BLASTFileFunc;
    mfa.bad_asn1	= FALSE;
    mfa.order		= 0;
    mfa.accession	= NULL;
    mfa.organism	= NULL;
    mfa.do_virtual	= FALSE;
    mfa.tech		= 0;
    mfa.no_sequence	= FALSE;
    mfa.formatdb	= TRUE;

    if (is_na)
            /* in case of "formatdb" we wont use this parameter */
        mfa.code = Seq_code_ncbi2na;
    else
        mfa.code = Seq_code_ncbistdaa;
    
    tfa.mfp = &mfa;
    tfa.is_na = is_na;
    if (group_segs == 3)  /* do 2 things */
    {
        mfa.do_virtual = TRUE;
        group_segs = 1;
    }
    tfa.group_segs = group_segs;
    tfa.last_indent = -1;
    tfa.parts = -1;
    tfa.seg = -1;
    tfa.got_one = FALSE;
    SeqEntryExplore(sep, (Pointer)&tfa, SeqEntryFasta);
    return tfa.got_one;
}


/*****************************************************************************
 *
 *   FastaFileFunc(key, buf, data)
 *       standard "write to file" callback
 *
 *****************************************************************************/
NLM_EXTERN Boolean BLASTFileFunc (BioseqPtr bsp, Int2 key, CharPtr buf, Uint4 buflen,
                                  Pointer data)
{
    FormatDBPtr	fdbp = (FormatDBPtr) data;
    Int4		SequenceLen;
    Uint4		i, total, index;
    
    switch (key)
    {
        case FASTA_ID:

            SequenceLen = bsp->length;
            fdbp->TotalLen += SequenceLen;
            
            if (fdbp->MaxSeqLen < SequenceLen)
                fdbp->MaxSeqLen = SequenceLen;
            
            if(OffsetAllocated <= fdbp->num_of_seqs) {
                OffsetAllocated += INDEX_ARRAY_CHUNKS;
                
                fdbp->DefOffsetTable = (Int4Ptr)Realloc(fdbp->DefOffsetTable, 
                                                        OffsetAllocated*sizeof(Uint4));
                fdbp->SeqOffsetTable = (Int4Ptr)Realloc(fdbp->SeqOffsetTable, 
                                                        OffsetAllocated*sizeof(Uint4));
                if(!fdbp->isProtein) {
                    fdbp->AmbOffsetTable = (Int4Ptr)Realloc(fdbp->AmbOffsetTable, 
                                                            OffsetAllocated*sizeof(Uint4));
                }
            }
            
            fdbp->DefOffsetTable[fdbp->num_of_seqs] = ftell(fdbp->fd_def); 
            fdbp->SeqOffsetTable[fdbp->num_of_seqs] = ftell(fdbp->fd_seq);
            
            if (FileWrite(buf, buflen, 1, fdbp->fd_def) != (Uint4) 1)
		return FALSE;
            if (FileWrite(" ", 1, 1, fdbp->fd_def) != (Uint4) 1)
		return FALSE;

                    /* Now adding new entried into lookup hash table */
            
            if((UpdateLookupInfo(buf, fdbp->lookup, 
                                 fdbp->num_of_seqs, fdbp->fd_stmp, fdbp->ParseMode)) 
               != LOOKUP_NO_ERROR) {
                return -1;
            } 
            break;
        case FASTA_DEFLINE:
	    if (FileWrite(buf, buflen, 1, fdbp->fd_def) != (Uint4) 1)
		return FALSE;
	    break;
        case FASTA_SEQLINE:
            if (FileWrite(buf, buflen, 1, fdbp->fd_seq) != (Uint4) 1)
		return FALSE;
            break;
        case FASTA_EOS:   /* end of sequence */
            if(fdbp->isProtein) {
                i=0;
                if (FileWrite(&i, 1, 1, fdbp->fd_seq) != (Uint4) 1)
			return FALSE;
            } else {
                    /* dump ambiguity characters. */
                fdbp->AmbOffsetTable[fdbp->num_of_seqs] = ftell(fdbp->fd_seq); /* Anyway... */
                
                    /* if AmbCharPtr is not NULL, then there was ambiguity. */
                if(fdbp->AmbCharPtr != NULL) 
                {
                        /* The first Uint4 holds the total number of ambig. bp. */
                    total = (*(fdbp->AmbCharPtr))+1;
                    for (index=0; index<total; index++) {
                        if (!FormatDbUint4Write(fdbp->AmbCharPtr[index], fdbp->fd_seq))
				return FALSE;
                    }
                    MemFree(fdbp->AmbCharPtr);
                    fdbp->AmbCharPtr = NULL;
                }
            }
            fdbp->num_of_seqs++;
            break;
        case FASTA_FORMATDB_AMB:
        {
            Int4 len;
            Char tmpbuff[1024];
                /* In case of "formatdb" nucleotides have to be compressed */
            
            fdbp->AmbCharPtr = NULL;
            
            if (bsp->seq_data_type == Seq_code_ncbi4na && bsp->seq_data != NULL){
                
                    /* ncbi4na require compression into ncbi2na */
                
                if((bsp->seq_data = 
                    BSCompressDNA(bsp->seq_data, bsp->length, 
                                  &(fdbp->AmbCharPtr))) == NULL) {
                    ErrLogPrintf("Error converting ncbi4na to ncbi2na. " 
                                 "Formating failed.\n");
                    return -1;
                }
                
                bsp->seq_data_type = Seq_code_ncbi2na; /* just for information */
            } else {
                    /* if sequence already in ncbi2na format we have to update last byte */
                Uint1 ch, remainder; 
                
                if((remainder = (bsp->length%4)) == 0) {
                    BSSeek(bsp->seq_data, bsp->length/4+1, SEEK_SET);
                    BSPutByte(bsp->seq_data, NULLB);
                } else {
                    BSSeek(bsp->seq_data, bsp->length/4, SEEK_SET);
                    ch = remainder + BSGetByte(bsp->seq_data);
                    BSSeek(bsp->seq_data, bsp->length/4, SEEK_SET);
                    BSPutByte(bsp->seq_data, ch);
                }
            }
                /* Now dumping sequence */
            
            BSSeek(bsp->seq_data, 0, SEEK_SET);
            while((len = BSRead(bsp->seq_data, tmpbuff, sizeof(tmpbuff))) != 0) {
                BLASTFileFunc(bsp, FASTA_SEQLINE, tmpbuff, len, data);
            }
            
            BLASTFileFunc(bsp, FASTA_EOS, NULL, 0, data);
        }
            break;
        default:
            break;
    }
    return TRUE;
}
    

/*******************************************************************************
 * Initializing FormatDB structure (see formatdb.h),
 ******************************************************************************* 
 * Parameters:
 *	dbname		- name of the input file
 *	isProtein	- true, if file with protein seqs
 *
 * Returns pointer to allocated FormatDB structure (FormatDBPtr)
 *	
 ******************************************************************************/

FormatDBPtr	FormatDBInit(const CharPtr dbname, const CharPtr dbtitle,
                             Boolean isProtein, Boolean ParseMode,
                             Boolean is_asn, Boolean asn_binmode)
{

    FormatDBPtr		fdbp;
    Char		filenamebuf[FILENAME_MAX];
    Uint4		i = 0;
    
    fdbp = (FormatDBPtr) MemNew (sizeof(*fdbp));
    
    fdbp->dbname = dbname;
    fdbp->DbTitle = dbtitle;
    fdbp->num_of_seqs = 0;
    fdbp->isProtein = isProtein;
    fdbp->ParseMode = ParseMode;
    fdbp->TotalLen=0, fdbp->MaxSeqLen=0;

        /* open input database */
    if (is_asn) {
#if 0        
        if ((fdbp->aip = AsnIoOpen (dbname, asn_binmode ? "rb":"r")) == NULL)
        {
            ErrLogPrintf("Cannot open input database file. Formating failed...\n");
            return NULL;
        }
        fdbp->fd = NULL;
#endif        
    }
    else {
        
        if((fdbp->fd = FileOpen(dbname, "r")) == NULL) {
            ErrLogPrintf("Cannot open input database file. Formating failed...\n");
            return NULL;
        }
        fdbp->aip = NULL;
    }
    
        /* open output BLAST files */
    
        /* Defline file */
    
    sprintf(filenamebuf, "%s.%chr", 
            dbname, fdbp->isProtein ? 'p' : 'n'); 
    fdbp->fd_def = FileOpen(filenamebuf, "wb");        
    
        /* Sequence file */
    
    sprintf(filenamebuf, "%s.%csq",
            dbname, fdbp->isProtein ? 'p' : 'n'); 
    fdbp->fd_seq = FileOpen(filenamebuf, "wb");        
    
    if (FileWrite(&i, 1, 1, fdbp->fd_seq) != (Uint4) 1) /* Sequence file started from NULLB */
	return NULL;
    
        /* Index file */
    
    sprintf(filenamebuf, "%s.%cin",
            dbname, fdbp->isProtein ? 'p' : 'n'); 
    fdbp->fd_ind = FileOpen(filenamebuf, "wb");      
    
        /* String (accession) index temporary file */

    fdbp->fd_stmp = NULL;
    if(ParseMode) {
        sprintf(filenamebuf, "%s.%ctm",
                dbname, fdbp->isProtein ? 'p' : 'n'); 
        fdbp->fd_stmp = FileOpen(filenamebuf, "wb");      
    }
    
    ErrLogPrintf("Version %s [%s]\n", BlastGetVersionNumber(), BlastGetReleaseDate()); 
    ErrLogPrintf("Started database file \"%s\"\n", dbname);

        /* Allocating space for offset tables */
    fdbp->DefOffsetTable = (Int4Ptr)MemNew(OffsetAllocated*sizeof(Uint4));
    fdbp->SeqOffsetTable = (Int4Ptr)MemNew(OffsetAllocated*sizeof(Uint4));
    if(!isProtein) 
        fdbp->AmbOffsetTable = (Int4Ptr)MemNew(OffsetAllocated*sizeof(Uint4));
    else
        fdbp->AmbOffsetTable = NULL;

        /* Allocating space for lookup table */
    
    if((fdbp->lookup = FASTALookupNew()) == NULL) {
        ErrLogPrintf("Error initializing Lookup structure. Formatting failed.\n");
        return NULL;
    }

    return fdbp;
}


/*******************************************************************************
 * Free memory allocated for given variable of FormatDB
 ******************************************************************************* 
 * Parameters:
 *	fdbp	- pointer to memory to be freed
 *
 * Returns NULL
 ******************************************************************************/


FormatDBPtr	FormatDBDestruct(FormatDBPtr fdbp)
{

    MemFree(fdbp->DefOffsetTable);
    MemFree(fdbp->SeqOffsetTable);
    
    if(!fdbp->isProtein) {
        MemFree(fdbp->AmbOffsetTable);
    }
    
    FASTALookupFree(fdbp->lookup);

    if (fdbp->fd)
        FileClose(fdbp->fd);
    FileClose(fdbp->fd_def);
    FileClose(fdbp->fd_ind);
    FileClose(fdbp->fd_seq);
    
    if (fdbp->aip)
        AsnIoClose(fdbp->aip);
    
    MemFree (fdbp);
    
    return NULL;
}


/*******************************************************************************
 * Pass thru each bioseq into given SeqEntry and write corresponding information
 * into "def", "index", ...., files
 ******************************************************************************* 
 * Parameters:
 *	fdbp	- pointer to memory to be freed
 *
 * Returns NULL
 ******************************************************************************/
static	Int2	process_sep (SeqEntryPtr sep, FormatDBPtr fdbp)
{

    Int4		SequenceLen;
    BioseqPtr		bsp = NULL;
    CharPtr		defline;
    Char		tmpbuff[1024];
    Int4		buffer_size=0, defline_len=0;
    CharPtr		buffer=NULL;
    Int4		len, id_length;
    Uint4Ptr		AmbCharPtr = NULL;
    Uint1		ch, remainder;
    Uint4		i, total, index;

    if (IS_Bioseq(sep))
        bsp = (BioseqPtr) sep->data.ptrvalue;
    else
            /* This is Bioseq-set.  Exit */
        return 0;

        /* Make a convertion to stadard form */

    if (fdbp->isProtein)
        BioseqRawConvert(bsp, Seq_code_ncbistdaa);
    
    SequenceLen = bsp->length;
    fdbp->TotalLen += SequenceLen;
    
    if (fdbp->MaxSeqLen < SequenceLen)
        fdbp->MaxSeqLen = SequenceLen;
        
    if(OffsetAllocated <= fdbp->num_of_seqs) {
        OffsetAllocated += INDEX_ARRAY_CHUNKS;
        
        fdbp->DefOffsetTable = (Int4Ptr)Realloc(fdbp->DefOffsetTable, 
                                          OffsetAllocated*sizeof(Uint4));
        fdbp->SeqOffsetTable = (Int4Ptr)Realloc(fdbp->SeqOffsetTable, 
                                          OffsetAllocated*sizeof(Uint4));
        if(!fdbp->isProtein) {
            fdbp->AmbOffsetTable = (Int4Ptr)Realloc(fdbp->AmbOffsetTable, 
                                              OffsetAllocated*sizeof(Uint4));
        }
    }
    
    fdbp->DefOffsetTable[fdbp->num_of_seqs] = ftell(fdbp->fd_def); 
    fdbp->SeqOffsetTable[fdbp->num_of_seqs] = ftell(fdbp->fd_seq);
    
        /* ---------------------- */

    if(fdbp->ParseMode == FALSE) 
    {
        sprintf(tmpbuff, "%s%d ", NON_SEQID_PREFIX, fdbp->num_of_seqs);
        if (FileWrite(tmpbuff, StringLen(tmpbuff), 1, fdbp->fd_def) != (Uint4) 1)
		return 1;
        defline = (CharPtr)bsp->descr->data.ptrvalue;
    }
    else
    {
        if (bsp->descr)
            defline_len = StringLen(BioseqGetTitle(bsp));
        else
            defline_len = 0;
        defline_len += 255;	/* Sufficient for an ID. */
        if (buffer_size < defline_len)
        {
            if (buffer)
                buffer = MemFree(buffer);
            buffer = MemNew((defline_len+1)*sizeof(Char));
            buffer_size = defline_len;
        }
        SeqIdWrite(bsp->id, buffer, PRINTID_FASTA_LONG, STRLENGTH);
        id_length = StringLen(buffer);
        buffer[id_length] = ' ';
        id_length++;
        StringCpy(&buffer[id_length], BioseqGetTitle(bsp));
        defline = buffer;
    }
    if (FileWrite(defline, StringLen(defline), 1, fdbp->fd_def) != (Uint4) 1)
	return 1;
    
        /* -------- Now adding new entried into lookup hash table */
    
    if((UpdateLookupInfo(defline, fdbp->lookup, 
                         fdbp->num_of_seqs, fdbp->fd_stmp, fdbp->ParseMode)) 
       != LOOKUP_NO_ERROR) {
        return -1;
    }

    defline = NULL;
    if (buffer)
	MemFree(buffer);
    
    if(!fdbp->isProtein)  {
        AmbCharPtr = NULL;
        if (bsp->seq_data_type == Seq_code_ncbi4na && bsp->seq_data != NULL){
            
                /* ncbi4na require compression into ncbi2na */
            
            if((bsp->seq_data = 
                BSCompressDNA(bsp->seq_data, bsp->length, 
                              &AmbCharPtr)) == NULL) {
                ErrLogPrintf("Error converting ncbi4na to ncbi2na. " 
                             "Formating failed.\n");
                return -1;
            }
            
            bsp->seq_data_type = Seq_code_ncbi2na; /* just for information */
        } else {
                /* if sequence already in ncbi2na format we have to update last byte */
            
            if((remainder = (bsp->length%4)) == 0) {
                BSSeek(bsp->seq_data, bsp->length/4+1, SEEK_SET);
                BSPutByte(bsp->seq_data, NULLB);
            } else {
                BSSeek(bsp->seq_data, bsp->length/4, SEEK_SET);
                ch = remainder + BSGetByte(bsp->seq_data);
                BSSeek(bsp->seq_data, bsp->length/4, SEEK_SET);
                BSPutByte(bsp->seq_data, ch);
            }
        }
    }
        /* Now dumping sequence */
    
    BSSeek(bsp->seq_data, 0, SEEK_SET);

    while((len = BSRead(bsp->seq_data, tmpbuff, sizeof(tmpbuff))) != 0) 
    {
        if (FileWrite(tmpbuff, len, 1, fdbp->fd_seq) != (Uint4) 1)
		return 1;
    }
    
    if(fdbp->isProtein) {
        i=0;
        if (FileWrite(&i, 1, 1, fdbp->fd_seq) != (Uint4) 1)
		return 1;
    } else {
            /* dump ambiguity characters. */
        fdbp->AmbOffsetTable[fdbp->num_of_seqs] = ftell(fdbp->fd_seq); /* Anyway... */
        
            /* if AmbCharPtr is not NULL, then there was ambiguity. */
        if(AmbCharPtr != NULL) 
        { /* The first Uint4 holds the total number of ambig. bp. */
            total = (*AmbCharPtr)+1;
            for (index=0; index<total; index++) {
                if (!FormatDbUint4Write(AmbCharPtr[index], fdbp->fd_seq))
			return 1;
            }
            MemFree(AmbCharPtr);
            AmbCharPtr = NULL;
        }
    }

    return 0;
}    





/*******************************************************************************
 * Finish stage - out offset tables, etc, into files.  Is to be called before
 * FormatDBDestruct()
 ******************************************************************************* 
 * Parameters:
 *	
 *
 * Returns  void
 ******************************************************************************/

static	Int2	finish_formatDB (FormatDBPtr fdbp) 
{
    Char	DBName[FILENAME_MAX];
    Int4	title_len;
    Char	dateTime[30];
    ISAMObjectPtr object;
    ISAMErrorCode error;
    Uint4	i;
    Char	filenamebuf[FILENAME_MAX];
    
   
    fdbp->DefOffsetTable[fdbp->num_of_seqs] = ftell(fdbp->fd_def); 
    
    if(!fdbp->isProtein) {
        fdbp->AmbOffsetTable[fdbp->num_of_seqs] = ftell(fdbp->fd_seq);
        fdbp->SeqOffsetTable[fdbp->num_of_seqs] =
            fdbp->AmbOffsetTable[fdbp->num_of_seqs];
    } else {
        fdbp->SeqOffsetTable[fdbp->num_of_seqs] = ftell(fdbp->fd_seq);
    }
    
        /* Parsing finished - now dumping index file */
    
    if(fdbp->ParseMode)
        FileClose(fdbp->fd_stmp);
    
        /* Information */
    
    if (!FormatDbUint4Write(FORMATDB_VER, fdbp->fd_ind))
	return 1;
    if (!FormatDbUint4Write(fdbp->isProtein, fdbp->fd_ind))
	return 1;
    
    if(fdbp->DbTitle != NULL)
        title_len = StringLen(fdbp->DbTitle);
    else
        title_len = 0;
    
    if (!FormatDbUint4Write(title_len, fdbp->fd_ind))
	return 1;
    
    if (title_len != 0)
        if (FileWrite(fdbp->DbTitle, title_len, 1, fdbp->fd_ind) != (Uint4) 1)
		return 1;
    
    Nlm_DayTimeStr(dateTime, TRUE, TRUE);
    if (!FormatDbUint4Write(StringLen(dateTime), fdbp->fd_ind))
	return 1;
    if (FileWrite(dateTime, StringLen(dateTime), 1, fdbp->fd_ind) != 1)
	return 1;
    
    if (!FormatDbUint4Write(fdbp->num_of_seqs, fdbp->fd_ind))
	return 1;
    if (!FormatDbUint4Write(fdbp->TotalLen, fdbp->fd_ind))
	return 1;
    if (!FormatDbUint4Write(fdbp->MaxSeqLen, fdbp->fd_ind))
	return 1;
    
        /* Offset tables */
    
    for(i=0; i <= fdbp->num_of_seqs; i++) {
        if (!FormatDbUint4Write(fdbp->DefOffsetTable[i], fdbp->fd_ind))
		return 1;
    }
    
    for(i=0; i <= fdbp->num_of_seqs; i++) {
        if (!FormatDbUint4Write(fdbp->SeqOffsetTable[i], fdbp->fd_ind))
		return 1;
    }
    if(!fdbp->isProtein) {
        for(i=0; i <= fdbp->num_of_seqs; i++) {
            if (!FormatDbUint4Write(fdbp->AmbOffsetTable[i], fdbp->fd_ind))
		return 1;
        }
    }
    
        /* Numeric lookup table sort & dump */
    
    if(fdbp->ParseMode && fdbp->lookup->used > 0) {
        FILE	*fd_lookup;
        
        sprintf(DBName, "%s.%cnd",
                fdbp->dbname, fdbp->isProtein ? 'p' : 'n'); 
        fd_lookup = FileOpen(DBName, "wb");          
        
        HeapSort(fdbp->lookup->table, fdbp->lookup->used/2,
                 sizeof(Uint4)*2, ID_Compare); 
        
        for(i=0; i < fdbp->lookup->used; i++) {
            if (!FormatDbUint4Write(fdbp->lookup->table[i], fd_lookup))
		return 1;
        }
        
        FileClose(fd_lookup);
        
            /* Now creating numeric ISAM index */
        
        sprintf(filenamebuf, "%s.%cni", 
                fdbp->dbname, fdbp->isProtein ? 'p' : 'n'); 
        
        if((object = ISAMObjectNew(ISAMNumeric, 
                                   DBName, filenamebuf)) == NULL) {
            ErrPostEx(SEV_ERROR, 0, 0, "Failed to create ISAM object.\n");
            return 1;
        }
        
        if((error = ISAMMakeIndex(object, 0)) != ISAMNoError) {
            if (error == ISAMNoOrder) {
                ErrPostEx(SEV_ERROR, 0, 0, "Failed to create index."
                          "  Possibly a gi included more than once in the database.\n", (long) error);
            } else {
                ErrPostEx(SEV_ERROR, 0, 0, "Failed to create index.\n", (long) error);
            }
            return 1;
        }
        ISAMObjectFree(object);
    }
    
        /* String file sorting */
    
    if(fdbp->ParseMode)
    {
        if (!FormatdbCreateStringIndex(fdbp->dbname, fdbp->isProtein))
		return 1;
    }
    
    ErrLogPrintf("Formated %d sequences\n", fdbp->num_of_seqs);

    return 0;
    
} /* end finish_formatDB() */


/*******************************************************/
static Boolean FormatdbCreateStringIndex(const CharPtr FileName, 
                                         Boolean ProteinType)
{
    SORTObjectPtr sop;
    Char filenamebuf[FILENAME_MAX], DBName[FILENAME_MAX];
    FILE *fd_out;
    CharPtr files;
    ISAMErrorCode error;
    ISAMObjectPtr isamp;

    /*  object for unique sorting */
    
    if((sop = SORTObjectNew(NULL, '\0', 0,
                            FALSE, TRUE)) == NULL) { 
        ErrPostEx(SEV_ERROR, 0, 0, "Failed to create SORT Object");
        return FALSE;
    }

    sprintf(filenamebuf, "%s.%ctm",
            FileName, ProteinType ? 'p' : 'n'); 
    
    sprintf(DBName, "%s.%csd",
            FileName, ProteinType ? 'p' : 'n'); 
    
    if((fd_out = FileOpen(DBName, "w")) == NULL)
    {
        return FALSE;
    }
    files = filenamebuf;
    
    if (SORTFiles(&files, 1, fd_out, sop) != SORTNoError)
    {
        ErrPostEx(SEV_ERROR, 0, 0, "SORTFiles failed");
	return FALSE;
    }
    SORTObjectFree(sop);

    FileClose(fd_out);
    FileRemove(filenamebuf);

    sprintf(filenamebuf, "%s.%csi",
            FileName, ProteinType ? 'p' : 'n'); 
    
    if((isamp = ISAMObjectNew(ISAMString, DBName, 
                              filenamebuf)) == NULL) {
        ErrPostEx(SEV_ERROR, 0, 0, "Creating of ISAM object failed");
        return FALSE;
    }
    
    if((error = ISAMMakeIndex(isamp, 0)) != ISAMNoError) {
        ErrPostEx(SEV_ERROR, 0, 0, "Creating of index failed with error code %ld\n", (long) error);
        ISAMObjectFree(isamp);
        return FALSE;
    } 
    
    ISAMObjectFree(isamp);
    
    return TRUE;
}

static Int4 UpdateLookupInfo(CharPtr defline, 
                             FASTALookupPtr lookup, 
                             Int4 num_of_seqs,
                             FILE *fd_stmp,
                             Boolean ParseSeqid
                             )
{
    CharPtr p, d = defline;
    Int4 i, gi = 0;
    Char TextId[ID_MAX_SIZE+1];
    SeqIdPtr sip, sip_tmp;
    
    if(defline == NULL)
        return LOOKUP_NO_ERROR;
    
    if(!ParseSeqid)
        return LOOKUP_NO_ERROR;
    
    for(p = d = defline; ;d = p + StringLen(TextId)) {
        
        MemSet(TextId, 0, sizeof(TextId));
        
        for(i=0; !isspace(*p) && i < ID_MAX_SIZE; p++,i++)
            TextId[i]=*p;
        
        if((sip = SeqIdParse(TextId)) == NULL) {/* Bad SeqId string */
            ErrLogPrintf("Sequence id \"%s\" is not parseable. "
                         "Formating failed at %s\n", TextId, defline);
            return ERR_SEQID_FAILED;
        }
     
        for(sip_tmp = sip; sip_tmp != NULL; sip_tmp = sip_tmp->next) {
            if(sip_tmp->choice == SEQID_GI) {
                gi = sip_tmp->data.intvalue;
                break;
            }
        }

        if(gi != 0) { /* GI not found */
            
            if((lookup->used + 2) >= lookup->allocated) {
                lookup->allocated += LOOKUP_CHUNK;
                lookup->table = (Int4Ptr)Realloc(lookup->table, 
                                                 lookup->allocated*(sizeof(Int4))); 
            }
            
            lookup->table[lookup->used] = gi;
            lookup->table[lookup->used+1] = num_of_seqs;
            lookup->used += 2;    
        }
        
        if(!SeqIdSetE2Index (sip, fd_stmp, num_of_seqs)) {
            ErrLogPrintf("SeIdSetE2Index failed. Exiting..\n");
            return FALSE;
        }
        
	sip = SeqIdSetFree(sip);
        
        if((p = StringChr(d, READDB_DEF_SEPARATOR)) == NULL)
            break;
        else
            p++;
    }
    return LOOKUP_NO_ERROR;
}

/* Size of variable that is manipulated, and swapped 
   for big/little endian stuff. */

static Boolean
FormatDbUint4Write(Uint4 number, FILE *fp)
  
{
  Uint4 value;

  /* If FORMATDB_SIZE changes, this must be changed. */
  value = Nlm_SwapUint4(number);	
  if (FileWrite(&(value), FORMATDB_SIZE, 1, fp) != (Uint4) 1)
	return FALSE;

  return TRUE;
}

static FASTALookupPtr FASTALookupNew(void) {
  FASTALookupPtr lookup;
  
  if((lookup = (FASTALookupPtr)MemNew(sizeof(FASTALookup))) == NULL)
    return NULL;
  if((lookup->table = (Int4Ptr)MemNew(LOOKUP_CHUNK*4)) == NULL)
    return NULL;
  
  lookup->allocated = LOOKUP_CHUNK;
  lookup->used = 0;
  return lookup;
}
static void FASTALookupFree(FASTALookupPtr lookup)
{
  MemFree(lookup->table);
  MemFree(lookup);
}
/* ------------------------------------------------------------------
                This is handler for HeapSort function
   ------------------------------------------------------------------*/
static int LIBCALLBACK ID_Compare(VoidPtr i, VoidPtr j)
{
  if (*(Int4Ptr)i > *(Int4Ptr)j)
    return (1);
  if (*(Int4Ptr)i < *(Int4Ptr)j)
    return (-1);
  return (0);
}

/*****************************************************************************
*
*   SeqIdE2Index(anp)
*   	atp is the current type (if identifier of a parent struct)
*       if atp == NULL, then assumes it stands alone (SeqId ::=)
*
*****************************************************************************/
static Boolean SeqIdE2Index (SeqIdPtr anp, FILE *fd, Int4 seq_num)
{
    Boolean retval = FALSE;
    TextSeqIdPtr tsip = NULL;
    ObjectIdPtr oid;
    PDBSeqIdPtr psip;
    Boolean do_gb = FALSE;
    Uint1 tmptype;
    CharPtr tmp, ptr=NULL;
    Char buf[81];
    Int4 length, i;
    DbtagPtr dbt;
    Uint1 chain = 0;

    if (anp == NULL)
        return FALSE;
    
    switch (anp->choice) {

    case SEQID_LOCAL:     /* local */
	oid = (ObjectIdPtr)(anp->data.ptrvalue);
	ptr = oid->str;
        break;
    case SEQID_GIBBSQ:    /* gibbseq */
        sprintf(buf, "%ld", (long)(anp->data.intvalue));
        ptr = buf;
        break;
    case SEQID_GIBBMT:    /* gibbmt */
        break;
    case SEQID_GIIM:      /* giimid */
        return TRUE;      /* not indexed */
    case SEQID_EMBL:      /* embl */
    case SEQID_DDBJ:      /* ddbj */
        do_gb = TRUE;     /* also index embl, ddbj as genbank */
    case SEQID_GENBANK:   /* genbank */
    case SEQID_PIR:       /* pir   */
    case SEQID_SWISSPROT: /* swissprot */
    case SEQID_OTHER:     /* other */
    case SEQID_PRF:       /* prf   */
        tsip = (TextSeqIdPtr)(anp->data.ptrvalue);
        break;
    case SEQID_PATENT:    /* patent seq id */
        break;
    case SEQID_GENERAL:   /* general */
        dbt = (DbtagPtr)(anp->data.ptrvalue);
        ptr = dbt->tag->str;
        break;
    case SEQID_GI:        /* gi */
        break;
    case SEQID_PDB:       /* pdb   */
	psip = (PDBSeqIdPtr)(anp->data.ptrvalue);
	ptr = psip->mol;
        chain = psip->chain;
        break;
    }
    
    SeqIdWrite(anp, buf, PRINTID_FASTA_SHORT, 80);

    length = StringLen(buf);
    for(i = 0; i < length; i++)
        buf[i] = TO_LOWER(buf[i]);
    fprintf(fd, "%s%c%d\n", buf, ISAM_DATA_CHAR, seq_num);

    if (ptr != NULL)   /* write a single string */
    {
	StringMove(buf, ptr);
        length = StringLen(buf);
        for(i = 0; i < length; i++)
            buf[i] = TO_LOWER(buf[i]);
        fprintf(fd, "%s%c%ld\n", buf, ISAM_DATA_CHAR, (long) seq_num);
 
        if (chain != 0) /* PDB only. */
        {
            fprintf(fd, "%s|%c%c%ld\n", buf, chain, ISAM_DATA_CHAR, (long) seq_num);
            fprintf(fd, "%s %c%c%ld\n", buf, chain, ISAM_DATA_CHAR, (long) seq_num);
        }
    }

    if (tsip != NULL) {   /* separately index accession and locus */
        if ((tsip->accession != NULL) && (tsip->name != NULL)) {
            tmp = tsip->accession;
            tsip->accession = NULL;
            SeqIdWrite(anp, buf, PRINTID_FASTA_SHORT, 80);
            length = StringLen(buf);
            for(i = 0; i < length; i++)
                buf[i] = TO_LOWER(buf[i]);
            fprintf(fd, "%s%c%d\n", buf, ISAM_DATA_CHAR, seq_num);
            tsip->accession = tmp;
            tmp = tsip->name;
            tsip->name = NULL;
            SeqIdWrite(anp, buf, PRINTID_FASTA_SHORT, 80);
            length = StringLen(buf);
            for(i = 0; i < length; i++)
                buf[i] = TO_LOWER(buf[i]);
            fprintf(fd, "%s%c%d\n", buf, ISAM_DATA_CHAR, seq_num);
            tsip->name = tmp;
        }

               /* now index as separate strings */
	if (tsip->name != NULL)
	{
		StringMove(buf, tsip->name);
		length = StringLen(buf);
            for(i = 0; i < length; i++)
                buf[i] = TO_LOWER(buf[i]);
            fprintf(fd, "%s%c%d\n", buf, ISAM_DATA_CHAR, seq_num);
	}
        if (tsip->accession != NULL)
        {
                StringMove(buf, tsip->accession);
                length = StringLen(buf);
            for(i = 0; i < length; i++)
                buf[i] = TO_LOWER(buf[i]);
            fprintf(fd, "%s%c%d\n", buf, ISAM_DATA_CHAR, seq_num);
	}
 
    }
    
    if (do_gb) {   /* index embl and ddbj as genbank */
        tmptype = anp->choice;
        anp->choice = SEQID_GENBANK;
        SeqIdE2Index(anp, fd, seq_num);
        anp->choice = tmptype;
    }

    retval = TRUE;
    return retval;
}

/*****************************************************************************
*
*   SeqIdSetE2Index(anp, e2p, settype, elementtype)
*
*****************************************************************************/
Boolean SeqIdSetE2Index (SeqIdPtr anp, FILE *fd, Int4 seq_num)
{
    SeqIdPtr oldanp;
    Boolean retval = FALSE;
    
    if (anp == NULL)
        return FALSE;
    
    oldanp = anp;

    while (anp != NULL) {
        if (!SeqIdE2Index(anp, fd, seq_num))
            goto erret;
        anp = anp->next;
    }
    
    retval = TRUE;
erret:
    return retval;
}
