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

   $Revision: 6.30 $

   File Description:  formats FASTA databases for use by BLAST

   $Log: formatdb.c,v $
   Revision 6.30  1999/09/10 16:30:35  shavirin
   Fixed problems with formating proteins by formatdb

   Revision 6.29  1999/09/09 18:25:51  shavirin
   Changed way to parse ASN.1. Added possibility to parse
   delta sequences.

   Revision 6.28  1999/08/25 20:20:27  shavirin
   Added -s option to create sparse indexes.

   Revision 6.27  1999/08/18 15:00:11  shavirin
   If title missing from args *.pal file will have basename as title.

   Revision 6.26  1999/08/03 16:38:56  shavirin
   Added function FD_CreateAliasFile() for multivolume formating.

   Revision 6.24  1999/07/23 18:59:01  shavirin
   Added support for creation of multivolume databases.

   Revision 6.23  1999/05/13 19:34:19  shavirin
   More changes towards dump from ID.

   Revision 6.21  1999/05/12 15:46:52  shavirin
   Changed parameter in function FDBAddSequence().

   Revision 6.20  1999/04/26 21:06:19  shavirin
   Fixed minor bug.

   Revision 6.19  1999/04/26 19:37:45  shavirin
   Dumping info set to FALSE.

   Revision 6.18  1999/04/26 14:53:16  shavirin
   Fixed memory leaks in FDBAddSequence() function.

   Revision 6.17  1999/04/21 21:44:34  shavirin
   Many functions were moved to "readdb.c" file.

   Revision 6.16  1999/03/21 19:16:59  madden
   Fix problem on round numbers

   Revision 6.15  1999/03/05 21:34:48  madden
   Changes for accession.version

   Revision 6.14  1999/02/04 18:01:48  madden
   Add -n option for basename

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

/* program's arguments */

#define NUMARG 11

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
    {"Input is a Seq-entry","F", NULL ,NULL ,TRUE,'e',ARG_BOOLEAN,0.0,0,NULL},
    { "Base name for BLAST files", 
      NULL, NULL, NULL, TRUE, 'n', ARG_STRING, 0.0, 0, NULL},
    { "Number of sequence bases to be created in the volume", 
      "0", NULL, NULL, TRUE, 'v', ARG_INT, 0.0, 0, NULL},
    { "Create indexes limited only to accessions - sparse", 
      "F", NULL, NULL, TRUE, 's', ARG_BOOLEAN, 0.0, 0, NULL}
};

/*#define db_title	(const CharPtr) dump_args[0].strvalue 
  #define db_file		(const CharPtr) dump_args[1].strvalue 
  #define LogFileName	(const CharPtr) dump_args[2].strvalue  
  #define is_protein	dump_args[3].intvalue
  #define parse_mode	dump_args[4].intvalue
  #define isASN		dump_args[5].intvalue
  #define asnbin		dump_args[6].intvalue
  #define is_seqentry	dump_args[7].intvalue
  #define base_name	(CharPtr) dump_args[8].strvalue */

#define Bases_In_Volume dump_args[9].intvalue

FDB_optionsPtr FDB_CreateCLOptions(void)
{
    FDB_optionsPtr options;
    
    options = MemNew(sizeof(FDB_options));
    
    if ( !GetArgs ("formatdb", NUMARG, dump_args) )
        return NULL;
    
    if ( !ErrSetLog (dump_args[2].strvalue) ) { /* Logfile */
        ErrShow();
    } else {
        ErrSetOpts (ERR_CONTINUE, ERR_LOG_ON);
    }
    
    options->db_title = StringSave(dump_args[0].strvalue);
    options->db_file = StringSave(dump_args[1].strvalue);
    options->LogFileName = StringSave(dump_args[2].strvalue);
    options->is_protein = dump_args[3].intvalue; 
    options->parse_mode = dump_args[4].intvalue; 
    options->isASN = dump_args[5].intvalue; 
    options->asnbin = dump_args[6].intvalue; 
    options->is_seqentry = dump_args[7].intvalue;
    options->base_name = StringSave(dump_args[8].strvalue);
    options->dump_info = FALSE;
    options->sparse_idx = dump_args[10].intvalue;

    return options;
}


Boolean FD_CreateAliasFile(CharPtr title, CharPtr basename, 
                           Int4 volumes, Boolean is_protein)
{
    Char filenamebuf[128], buffer[64];
    time_t tnow;
    Int4 i;
    FILE *fd;

    sprintf(filenamebuf, "%s.%cal", basename, is_protein? 'p' : 'n'); 

    if((fd = FileOpen(filenamebuf, "wb")) == NULL)
        return FALSE;
    
    tnow = time(NULL);
    fprintf(fd, "#\n# Alias file created %s#\n#\n", ctime(&tnow));
    
    if(title != NULL)
        fprintf(fd, "TITLE %s\n#\n", title);
    else if (basename != NULL)
        fprintf(fd, "TITLE %s\n#\n", basename);
    else
        fprintf(fd, "#TITLE\n#\n");
    
    /* Now printing volume databases */
    fprintf(fd, "DBLIST ", title);
    
    for(i = 0; i < volumes; i++) {
        fprintf(fd, "%s.%02d ", basename, i);
    }
    fprintf(fd, "\n#\n");
    
    fprintf(fd, "#GILIST\n#\n");
    fprintf(fd, "#OIDLIST\n#\n");
    
    FileClose(fd);
    
    return TRUE;
}
/* main() */

Int2 Main(void) 
{
    SeqEntryPtr sep;
    FormatDBPtr	fdbp;
    Uint1 group_segs = 0;
    FDB_optionsPtr options;
    BioseqPtr bsp;
    BioseqSetPtr bssp;
    ByteStorePtr seq_data_new;
    Int4 count = 0, volume = 0;
    Char basename[128], filenamebuf[128];
    FILE *fd;
    /* get arguments */

    if((options = FDB_CreateCLOptions()) == NULL)
        return 1;
    
    if(options->base_name != NULL)
        StringCpy(basename, options->base_name);
    else
        StringCpy(basename, options->db_file);
    
    if(Bases_In_Volume > 1) {
        sprintf(filenamebuf, "%s.%02d", basename, volume); 
        MemFree(options->base_name);
        options->base_name = StringSave(filenamebuf);
        volume++;
    }
    
    
    /* Initialize formatdb structure */
    if ((fdbp = FormatDBInit(options)) == NULL)
        return 2;        
    
    /* Input database file maybe either in ASN.1 or in FASTA format */
    
    if (!options->isASN) {
        /* FASTA format of input database */

        if((fd = FileOpen(options->db_file, "r")) == NULL)
            return 3;
        
        /* Get sequences */
        while ((sep = FastaToSeqEntryEx(fd, 
                                        (Boolean)!options->is_protein,
                                        NULL, options->parse_mode)) != NULL) {
            
            if(!IS_Bioseq(sep)) { /* Not Bioseq - failure */
                ErrLogPrintf("Error in readind Bioseq Formating failed.\n");
                return 4;
            }
            
            SeqEntrySetScope(sep);
            bsp = (BioseqPtr) sep->data.ptrvalue;            
            
            if(Bases_In_Volume >= 1) {
                if(count > Bases_In_Volume) { 
                    /* starting new volume ? */
                    count = 0;
                    if(FormatDBClose(fdbp))
                        return 9;
                    
                    if(Bases_In_Volume > 1) {
                        sprintf(filenamebuf, "%s.%02d", basename, volume); 
                        MemFree(options->base_name);
                        options->base_name = StringSave(filenamebuf);
                        volume++;
                    }
                    
                    if ((fdbp = FormatDBInit(options)) == NULL)
                        return 2;
                }
                count += bsp->length;
            }
            
            FDBAddBioseq(fdbp, bsp);            

            SeqEntryFree(sep);
        }
        FileClose(fd);

        /* Writting multi-volume pointer file */
        if(Bases_In_Volume > 1) {
            
            FD_CreateAliasFile(options->db_title, basename, volume, 
                               options->is_protein);
        }

    } else {
        /* ASN.1 format of input database */
        AsnTypePtr atp, atp2;
        AsnModulePtr amp;
        
        if (! SeqEntryLoad())
            ErrShow();
        
        /* get pointer to all loaded ASN.1 modules */
        
        /* get pointer to all loaded ASN.1 modules */
        amp = AsnAllModPtr();

        if (amp == NULL) {
            ErrLogPrintf("Could not load ASN.1 modules.\n");
            return 5;
        }
        
        /* get the initial type pointers */
        
        atp = AsnFind("Bioseq-set");
        if (atp == NULL) {
            ErrLogPrintf("Could not get type pointer for Bioseq-set.\n");
            return 6;
        }
        
        atp2 = AsnFind("Bioseq-set.seq-set.E");
        if (atp2 == NULL) {
            ErrLogPrintf("Could not get type pointer for Bioseq-set.seq-set.E\n");
            return 7;
        }
        
        if ((fdbp->aip = AsnIoOpen (options->db_file, 
                                    options->asnbin ? "rb":"r")) == NULL) {
            ErrLogPrintf("Cannot open input database file. Formating failed...\n");
            return 8;
        }
        
        if (options->is_seqentry) {
            /* Seq entry */
            sep = SeqEntryAsnRead(fdbp->aip, NULL);
            FDBAddSeqEntry(fdbp, sep); 
            SeqEntryFree(sep);
        } else {
            /* Bioseq-set */
            
            while ((atp = AsnReadId(fdbp->aip, amp, atp)) != NULL) {
                if (atp == atp2) {   /* top level Seq-entry */
                    sep = SeqEntryAsnRead(fdbp->aip, atp);

                    FDBAddSeqEntry(fdbp, sep);
                    SeqEntryFree(sep);
                } else {
                    AsnReadVal(fdbp->aip, atp, NULL);
                }
            }
        } /* end "if Bioseq or Bioseq-set */
        

    } /* end "if FASTA or ASN.1" */
    
    /* Dump indexes, deallocate structure, arrays, etc. */

    if(FormatDBClose(fdbp))
        return 9;
    
    return 0;
    
} /* main()*/



