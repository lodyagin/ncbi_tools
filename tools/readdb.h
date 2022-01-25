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
*/

/*****************************************************************************

File name: readdb.h

Author: Tom Madden

Contents: defines and prototypes used by readdb.c and formatdb.c.

******************************************************************************/

/*
* File Name: readdb.h
*
* Author: Tom Madden
*
* Version Creation Date:   3/21/95
*
* $Revision: 6.28 $
*
* File Description: 
*       Functions to rapidly read databases from files produced by formatdb.
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
* $Log: readdb.h,v $
* Revision 6.28  1999/05/06 15:25:27  egorov
* Remove static function declaration
*
* Revision 6.27  1999/04/26 14:36:29  shavirin
* Added ability to dump statistics.
*
* Revision 6.26  1999/04/21 22:55:39  kans
* was not checked in
*
* Revision 6.25  1999/02/22 21:48:03  egorov
* Optimize GIs2OIDs not reinitializing ISAM indicies for non-exclisive databases, but use already initialized rdfp's field for that.
*
* Revision 6.24  1999/02/05 13:47:05  madden
* Add basename for formatdb
*
* Revision 6.23  1998/12/14 21:49:23  egorov
* new max gi number memeber in CommonIndexHead structure and therefore no need for COMMON_INDEX_TABLE_SIZE
*
* Revision 6.22  1998/12/14 16:05:36  egorov
* *** empty log message ***
*
* Revision 6.21  1998/09/14 15:11:19  egorov
* Add support for Int8 length databases; remove unused variables
*
* Revision 6.20  1998/08/27 15:02:37  madden
* Added LIBCALL for readdb_get_sequence_ex
*
* Revision 6.19  1998/08/24 14:59:57  madden
* readdb_get_sequence_ex function
*
* Revision 6.18  1998/08/11 17:49:48  madden
* is_na becomes is_aa
*
* Revision 6.17  1998/07/01 14:03:07  egorov
* Fix bug with a thread freeing CommonIndex: add new flag to rdfp
*
* Revision 6.16  1998/06/26 16:51:15  egorov
* Fix CommonIndex bugs
*
* Revision 6.15  1998/06/24 21:03:40  egorov
* Remove memory leaks
*
* Revision 6.12  1998/05/22 20:19:54  madden
* Changes to fix multi-db search bug
*
* Revision 6.11  1998/02/26 22:34:24  madden
* Changes for 16 bit windows
*
* Revision 6.10  1998/02/11 17:49:38  madden
* Added structures and prototypes for formatdb to take ASN.1 as input
*
* Revision 6.9  1998/01/16 22:03:00  madden
* Added init_indices Boolean
*
* Revision 6.8  1997/11/26 22:48:38  madden
* Added readdb_parse_db_names for multiple db searches
*
* Revision 6.7  1997/11/07 16:16:36  shavirin
* Added definition of new function readdb_acc2fastaEx()
*
* Revision 6.6  1997/10/24 19:08:16  madden
* Added ReadDBGetDb and ReadDBGetDbId
*
* Revision 6.5  1997/09/24 22:37:06  madden
* Added readdb_destruct_element
*
* Revision 6.4  1997/09/16 16:31:40  madden
* More changes for multiple db runs
*
* Revision 6.3  1997/09/12 19:55:38  madden
* Added readdb_compare
*
* Revision 6.2  1997/09/11 18:49:40  madden
* Changes to enable searches against multiple databases.
*
* Revision 6.1  1997/08/27 14:46:59  madden
* Changes to enable multiple DB searches
*
* Revision 6.0  1997/08/25 18:53:59  madden
* Revision changed to 6.0
*
* Revision 1.26  1997/05/12 21:34:05  madden
* readdb_new allows indeterminate database type
*
* Revision 1.25  1997/05/12 21:11:42  shavirin
* Added definition for function readdb_acc2fasta()
*
* Revision 1.23  1997/05/07 21:04:02  madden
* Added prototype for SeqId2OrdinalId and changed FORMATDB_VER
*
* Revision 1.22  1997/05/01 17:26:58  shavirin
* Added definition for the function readdb_seqid2fasta()
*
 * Revision 1.21  1997/02/25  22:16:32  shavirin
 * Changes in accordance to ISAM API changes
 *
 * Revision 1.20  1997/02/25  16:28:38  shavirin
 * Added new entries in ReadDBFILEPtr structure to do search by gi
 * number.
 *
 * Revision 1.19  1996/12/19  16:29:56  madden
 * Changes to eliminate ".nac" file for nucl.
 *
 * Revision 1.18  1996/12/17  21:34:46  madden
 * Changes to allow deflines for inidividual entries to be retrieved.
 *
 * Revision 1.17  1996/12/11  18:42:36  madden
 * Added prototypes for BioseqFetch functions.
 *
 * Revision 1.16  1996/11/27  16:39:11  madden
 * Added functions to return filename and date.
 *
 * Revision 1.15  1996/11/26  19:54:27  madden
 * Added check for database in standard places.
 *
 * Revision 1.14  1996/11/22  19:05:48  madden
 * removed ifdef for OLD_BIT_ORDER.
 *
 * Revision 1.13  1996/11/08  21:45:03  madden
 * Removed function readdb_get_partial_unpacked_sequence.
 *
 * Revision 1.12  1996/11/07  22:33:00  madden
 * Added prototype for readdb_ambchar_present.
 *
 * Revision 1.11  1996/11/04  18:50:20  shavirin
 * Added definitions for ambiguity information pointers
 *
 * Revision 1.10  1996/10/31  16:29:55  shavirin
 * Changed definitions due to reverce of residues in BLAST database
 * for nucleotide sequences from (4321) to (1234)
 * New dumper now required to create BLAST databases.
 *
 * Revision 1.9  1996/09/27  19:12:17  madden
 * Added function readdb_get_bioseq to obtain a BioseqPtr from the BLAST databases.
 *
 * Revision 1.8  1996/09/26  15:09:21  madden
 * Corrected misplaced comment.
 *
 * Revision 1.7  1996/09/23  14:37:35  madden
 * Replaced CharPtr (for sequence) with Uint1Ptr.
 *
 * Revision 1.6  1996/09/20  21:59:16  madden
 * *** empty log message ***
 *
 * Revision 1.5  1996/09/13  20:01:52  madden
 * defined READDB_COMPRESSION_RATIO
 *
 * Revision 1.4  1996/09/13  18:55:04  madden
 * Added function readdb_get_partial_unpacked_sequence.
 *
 * Revision 1.3  1996/08/29  20:42:01  madden
 * memory mapping moved to the corelib (in ncbimem.[ch]).
 *
 * Revision 1.2  1996/08/07  18:32:05  madden
 * Moved define of MMAP_AVAIL from readdb.h to readdb.c
 *
 * Revision 1.1  1996/08/05  19:48:21  madden
 * Initial revision
 *
 * Revision 1.12  1996/08/02  14:20:06  madden
 * Added readdb_attach function.
 *
 * Revision 1.11  1996/07/31  13:09:17  madden
 * Changes for partial copy of ReadDB structure.
 *
 * Revision 1.10  1996/07/25  20:45:20  madden
 * Change to arguments of readdb_get_sequence.
 *
 * Revision 1.9  1996/07/25  12:56:15  madden
 * readdb_get_sequence changed to allow for systems w/o mmap.
 *
 * Revision 1.8  1996/06/20  17:00:11  madden
 * Added "__cplusplus" define.
 *
 * Revision 1.7  1996/06/20  16:16:36  madden
 * Replaced int's with Int4's.
 *
 * Revision 1.6  1996/05/16  19:50:15  madden
 * Added documentation block.
 *
 * Revision 1.5  1996/04/22  21:42:07  madden
 * New prototype for readdb_get_sequence
 *
 * Revision 1.4  1996/04/11  14:30:06  madden
 * Memory-mapping added.
 *
 * Revision 1.3  1996/03/29  21:28:30  madden
 * Added function readdb_get_sequence_length.
 *
 * Revision 1.2  1996/03/28  20:42:36  madden
 * Added functions readdb_get_title, readdb_is_prot and
 * readdb_get_formatdb_version.
 *
 * Revision 1.1  1996/03/26  19:38:08  madden
 * Initial revision
 *
 *
*/

#ifndef _READDB_
#define _READDB_

#ifdef __cplusplus
extern "C" {
#endif

/****************************************************************************/
/* INCLUDES */
/****************************************************************************/

#include <ncbi.h>
#include <objloc.h>
#include <sequtil.h>
#include <ncbisam.h>

/****************************************************************************/
/* DEFINES */
/****************************************************************************/

/* Defines used to retrieve a base out of a packed byte. */
/* x should be unsigned (Uint1) to avoid sign extension problems. */

#define READDB_UNPACK_BASE_1(x) ((x)>>6)
#define READDB_UNPACK_BASE_2(x) (((x)>>4) & 0x03)
#define READDB_UNPACK_BASE_3(x) (((x)>>2) & 0x03)
#define READDB_UNPACK_BASE_4(x) ((x) & 0x03)

/* Default location of the databases. */
#ifdef OS_UNIX    
#define BLASTDB_DIR "/usr/ncbi/db/blast"
#else
#define BLASTDB_DIR ""
#endif    

/* Compress 4 bytes to one. */
#define READDB_COMPRESSION_RATIO 4

/* Character used to separate deflines from different entries that all
belong to the same sequence. */
#define READDB_DEF_SEPARATOR '\001'

/* Choices for whether it's a protein db or not. */
#define READDB_DB_IS_NUC 0
#define READDB_DB_IS_PROT 1
#define READDB_DB_UNKNOWN 2

/* The following variables are shared by formatdb and readdb. */
/* version of formatdb. */
#define FORMATDB_VER 3

/****************************************************************************/
/* TYPEDEFS */
/****************************************************************************/

typedef struct nlm_mfile {
	Nlm_MemMapPtr mem_mapp;	/* structure containing mem-map info, 
				produced by Nlm_MemMapInit. */
	FILE PNTR fp;		/* FILE pointer. */
	Uint1Ptr  mmp_begin,	/* beginning of mmap'ed are. */
		  mmp,		/* present position of mmap'ed pointer. */
		  mmp_end;	/* end of mmap'ed area. */
	Int4	  file_size;	/* size of file that is mmap'ed. */
	Boolean   mfile_true;	/* If TRUE then mmap succeeded. */
	Boolean   contents_allocated; /* If TRUE, the contents have been allocated
					and are not merely a copy. */
} NlmMFILE, PNTR NlmMFILEPtr;

/*
Open the file and initialze the memory mapping.
*/
NlmMFILEPtr LIBCALL NlmOpenMFILE PROTO((CharPtr name));

/*
Undo the memory mapping.
*/
NlmMFILEPtr LIBCALL NlmCloseMFILE PROTO((NlmMFILEPtr mfp));

/*
Read "nitems" of size "size" from a memory mapped file into "buffer"
usig the memory-mapped file given by "mfp".
*/
Int4 LIBCALL NlmReadMFILE PROTO((Uint1Ptr buffer, size_t size, Int4 nitems, NlmMFILEPtr mfp));

/*
	"fseek" to a point in the memory mapped file.
*/
Int4 LIBCALL NlmSeekInMFILE PROTO((NlmMFILEPtr mfp, long offset, Int4 ptrname));

/*
        What is the offset (in bytes) to the beginning of the file.
        Analog to ftell.
*/
Int4 LIBCALL NlmTellMFILE PROTO((NlmMFILEPtr mfp));


/*
	Common index structures
 */

/*
#define INDEXFN		"/net/cruncher/usr/ncbi/db/disk.blast/blast2/comindex/mmfile.idx"
#define RECORDFN	"/net/cruncher/usr/ncbi/db/disk.blast/blast2/comindex/mmfile.rec"
*/
#define INDEX_FN	"comindex.mm"
#define DB_CONFIG_FN	"dblist.txt"

typedef struct  CommonIndex{
    Int4        dbmask; 	/* mask to define which db contains the GI */
    Int4        oftenOID;       /* ordinal ID for the GI in most often DB */
} CommonIndex, *CommonIndexPtr;
 
typedef	struct	CommonIndexResult {
    Int4	gi;	/* GI */
    Int4	oid;	/* OID */
    Int2	dbid;	/* database ID */
    struct CommonIndexResult *next;	/* make a list */
} CommonIndexResult, *CommonIndexResultPtr;

/* Data bases */
 
typedef struct	DataBaseID {
    CharPtr	name;	/* database name like gss, nr, etc */
    Char	id;	/* integer ID, value from 0, to 32, used for bitmasks */
    Boolean	isprot;	/* says TRUE if database contains proteins, FALSE otherwise */
} DataBaseID, *DataBaseIDPtr;

typedef struct	CommonIndexHead {
    CommonIndexPtr	ci;
    Nlm_MemMapPtr	memmap;
    Int2		num_of_DBs;
    DataBaseIDPtr	dbids;
    Int4		maxgi; /* maximum GI number permited */
} CommonIndexHead, *CommonIndexHeadPtr;

typedef struct read_db_file {
	struct read_db_file PNTR next;
/* the contents of this struct. were allocated, or not.  Does NOT include
the actual structure and buffer, below. */
	Boolean contents_allocated,
		indices_initialized;	/* The indices were initialized. */
	CharPtr filename;	/* name of the input (w/o extensions). */
/* The files pointers for "file" (above), the index file, the file 
containing the headers, and the sequence file. */
        NlmMFILEPtr indexfp, headerfp, sequencefp;
	Int4	header_index_offset;	/* offset to beginning of header index in indexfp. */
	Boolean is_prot; /* If TRUE, sequence is protein, otherwise dna. */
	CharPtr title,	/* Database Title. */
		date;	/* Date and time database was prepared. */
	Int4 num_seqs, /* Number of sequences in the database. */
	      formatdb_ver;	/* Version of formatdb used. */
	Int4 	start,	/* 1st ordinal id in this file. */
		stop;	/* last ordinal id in this file. */
	Uint4 totlen,	/* Total length of database. */
	      maxlen;	/* Length of longest sequence in database. */
/* The "index" arrays specify the offsets (in files) of the header and 
sequence information. */
	Uint4Ptr header_index,	sequence_index, ambchar_index;	
	Uint4Ptr header_index_start,	sequence_index_start, ambchar_index_start;	
/* Buffer and allocated amount of this buffer.  These should always be
NULL (i.e., NOT USED) if mem-mapping is used; only used to store sequence
if there is no mem-mapping or it failed. */
        ISAMObjectPtr nisam_opt;  /* Object for numeric search */
        ISAMObjectPtr sisam_opt;  /* Object for string search */
	Uint1Ptr buffer;
	Int4 allocated_length;
	Boolean			handle_common_index; /* TRUE only for a initial thread;  needs for proper freeing of the CommonIndex */
	CommonIndexHeadPtr	cih;	/* head of the common index */
	Int2			filebit;/* bit corresponding to the DB file */
} ReadDBFILE, PNTR ReadDBFILEPtr;

/* Function prototypes */
Int4    GI2OID(CommonIndexHeadPtr cih, Int4 gi, Int4 dbmask, Int2 *dbid, ReadDBFILEPtr rdfp);
Int2	DBShift(Int2 num_of_DBs, DataBaseIDPtr dbids, CharPtr dbname, Boolean is_prot);
CharPtr	DBName(Int2 num_of_DBs, DataBaseIDPtr dbids, Int2 shift);
Boolean	DBisProt(Int2 num_of_DBs, DataBaseIDPtr dbids, Int2 shift);
CommonIndexResultPtr	GIs2OIDs(CommonIndexHeadPtr cih,
			Int4Ptr gis, Int4 number_of_gis, Int4 dbshift, ReadDBFILEPtr rdfp);
Int2	SeniorBit(Int4	bitmask);
CommonIndexHeadPtr	CommonIndexInit(CharPtr indexfilename);
void	CommonIndexDestruct(CommonIndexHeadPtr cihp);
Int2	bit_engine_firstbit (Int4 word);
Int2Ptr	bit_engine_arr(Int4 word);
Int2	bit_engine_numofbits(Int4 word);
Int2	ParseDBConfigFile(DataBaseIDPtr *dbidsp, CharPtr path);

/* mmap's */
 
NLM_EXTERN Nlm_MemMapPtr EA_MemMapInit(const Nlm_Char PNTR name, Boolean readonly);

/****************************************************************************/
/* FINCTION DEFINITIONS */
/****************************************************************************/
/* 
Intitialize the readdb structure using the database "filename".
If no database is used, set filename to NULL. 
*/
ReadDBFILEPtr LIBCALL readdb_new PROTO((CharPtr filename, Uint1 is_prot));

/*
	init_indices should be TRUE if entire database is to be searched, otherwise
	it can be FALSE.
*/
ReadDBFILEPtr LIBCALL readdb_new_ex PROTO((CharPtr filename, Uint1 is_prot, Boolean init_indices));


/* 
Deallocate the ReadDBFILEPtr. 
*/
ReadDBFILEPtr LIBCALL readdb_destruct PROTO((ReadDBFILEPtr readdb));

ReadDBFILEPtr LIBCALL readdb_destruct_element PROTO((ReadDBFILEPtr rdfp));


/*
        Attach to an already open ReadDBFILEPtr.  Duplicate the
        indexfp, sequencefp, and headerfp structures as the pointers
        there (i.e., mmp) will need to be manipulated.  Do not
        change the FILE PNTR fp.
*/
ReadDBFILEPtr LIBCALL readdb_attach PROTO((ReadDBFILEPtr rdfp));

/*
        Checks whether a ReadDBFILEPtr is the original, or just attaced.
        It does this by checking the rdfp->contents_allocated flag.
*/
Boolean LIBCALL readdb_copy PROTO((ReadDBFILEPtr rdfp));

/*
	Checks two ReadDBFILEPtr to see if they refer to the same
	database.
*/
Boolean LIBCALL readdb_compare PROTO((ReadDBFILEPtr rdfp1, ReadDBFILEPtr rdfp2));

/*
        Get total length and number of sequences in multiple databases.
*/

Boolean LIBCALL readdb_get_totals PROTO((ReadDBFILEPtr rdfp_list, Int8Ptr total_len, Int4Ptr total_num));


/* 
Get the sequence with sequence_number and put it in buffer.  No memory
is allocated for this if memory-mapped files are used, otherwise it is.  
Return the length of the sequence.
*/
Int4 LIBCALL readdb_get_sequence PROTO((ReadDBFILEPtr rdfp, Int4 sequence_number, Uint1Ptr PNTR buffer));

/* 
	Gets the sequence number "sequence_number".  The sequence returned includes
	all ambiguity information.  THis funciton should only be used for nucleic
	acid sequences, for proteins use readdb_get_sequence.

	buffer contains the sequence and is reallocated if *buffer_length is not long enough.

	The length of the sequence requested is the return value.
	protein sequences are always returned as Seq_code_ncbistdaa,
	nucleotide sequences as Seq_code_ncbi4na.
*/
Int4 LIBCALL readdb_get_sequence_ex PROTO((ReadDBFILEPtr rdfp, Int4 sequence_number, Uint1Ptr PNTR buffer, Int4 *buffer_length));

/* Gets sequence number by gi number. Returnes -1 if gi not found or
   other negative value if NISAM library faults. Non-negative value
   means success. Use numeric ISAM indexes.
*/
Int4 LIBCALL readdb_gi2seq(ReadDBFILEPtr rdfp, Int4 gi);

/* Gets sequence number by SeqId number. Returnes -1 if gi not found or
   other negative value if SISAM library faults. Non-negative value
   means success. Use string ISAM indexes.
*/
Int4 LIBCALL readdb_seqid2fasta(ReadDBFILEPtr rdfp, SeqIdPtr sip);

/* Gets sequence number by Accession/Locus string. Returnes -1 
   if accession not found or
   other negative value if SISAM library faults. Non-negative value
   means success. Use string ISAM indexes.
*/
Int4 LIBCALL readdb_acc2fasta(ReadDBFILEPtr rdfp, CharPtr string);

/* Gets array of sequence numbers by Accession/Locus string. Returnes -1 
   if accession not found or
   other negative value if SISAM library faults. Non-negative value
   means success. Use string ISAM indexes.
*/
Int4 LIBCALL readdb_acc2fastaEx(ReadDBFILEPtr rdfp, CharPtr string,
                                Int4Ptr PNTR ids, Int4Ptr count);

/*
Gets a BioseqPtr containing the sequence in sequence_number.
*/
BioseqPtr LIBCALL readdb_get_bioseq PROTO((ReadDBFILEPtr rdfp, Int4 sequence_number));

/*
Get the length of the sequence.
*/
Int4 LIBCALL readdb_get_sequence_length PROTO((ReadDBFILEPtr rdfp, Int4 sequence_number));

/* 
Get the ID and definition for the sequence with sequence_number.
It is the caller's RESPONSIBILITY to DEALLOCATE "id" and "description". 
*/
Boolean LIBCALL readdb_get_descriptor PROTO((ReadDBFILEPtr rdfp, Int4 sequence_number, SeqIdPtr PNTR id, CharPtr PNTR description));

/*
Get the ID's and headers for a sequence. 
*/
Boolean LIBCALL
readdb_get_header PROTO((ReadDBFILEPtr rdfp, Int4 sequence_number, Uint4Ptr header_index , SeqIdPtr PNTR id, CharPtr PNTR description));

/* 
 Get the Int4Ptr to ambiguity buffer
*/
Boolean  LIBCALL readdb_get_ambchar PROTO((ReadDBFILEPtr rdfp, Int4 sequence_number, Uint4Ptr PNTR ambchar_return));

/*
	Check whether ambiguity characters are present in the sequence. 
*/
Boolean LIBCALL readdb_ambchar_present PROTO((ReadDBFILEPtr rdfp, Int4 sequence_number));

/* 
Get the total length (in bp or residues) of the database. 
*/
Int8 LIBCALL readdb_get_dblen PROTO((ReadDBFILEPtr rdfp));

/* 
Get the number of entries in the database. 
*/
Int4 LIBCALL readdb_get_num_entries PROTO((ReadDBFILEPtr rdfp));

/* 
Get the total number of entries in all the files.
*/
Int4 LIBCALL readdb_get_num_entries_total PROTO((ReadDBFILEPtr rdfp));


/* 
Get the length of the longest sequence in the database. 
*/
Int4 LIBCALL readdb_get_maxlen PROTO((ReadDBFILEPtr rdfp));

/* 
Get the title (i.e., name) of the database. 
NOTE: the CharPtr returned is not owned by the caller! 
*/
CharPtr LIBCALL readdb_get_title PROTO((ReadDBFILEPtr rdfp));

/* 
Get the name of the file used for formatting.
NOTE: the CharPtr returned is not owned by the caller! 
*/
CharPtr LIBCALL readdb_get_filename PROTO((ReadDBFILEPtr rdfp));

/* 
Get the date the database was formatted. 
NOTE: the CharPtr returned is not owned by the caller! 
*/
CharPtr LIBCALL readdb_get_date PROTO((ReadDBFILEPtr rdfp));

/* 
Is this a protein database? 
*/
Boolean LIBCALL readdb_is_prot PROTO((ReadDBFILEPtr rdfp));

/*
        Parses the databases names (if more than one) from
        'filenames' into buffer.  buffer should already be
        long enough and allocated.  The funciton should be
        repeatedly called until TRUE is returned.
*/
Boolean LIBCALL readdb_parse_db_names PROTO((CharPtr PNTR filenames, CharPtr buffer));

/* 
Get the version of formatdb used on this database. 
*/
Int4 LIBCALL readdb_get_formatdb_version PROTO((ReadDBFILEPtr rdfp));

/* For the BioseqFetch functions. */

Boolean LIBCALL ReadDBBioseqFetchEnable PROTO((CharPtr program, CharPtr dbname, Boolean is_na, Boolean now));

void LIBCALL ReadDBBioseqFetchDisable PROTO((void));

/* Converts a SeqIdPtr to an ordinal_id, which readdb can use to look
up sequences etc.  Negative numbers are returned if the SeqIdPtr
cannot be converted. */
Int4 SeqId2OrdinalId PROTO((ReadDBFILEPtr rdfp, SeqIdPtr sip));

/*
	Returns the ReadDBFILEPtr by the database ID.
*/
ReadDBFILEPtr ReadDBGetDb PROTO((ReadDBFILEPtr list, Int2 db_id));

/*
	Returns the Database ID.
*/
Int2 ReadDBGetDbId PROTO((ReadDBFILEPtr list, ReadDBFILEPtr target));


/********************/
/*     formatdb     */
    
    /* Type definitions */

typedef struct FASTALookup {
    Int4Ptr table;          /* Main buffer for gi/fasta_id pairs */
    Int4    allocated;      /* Nunber of Uint4 allocated */
    Int4    used;           /* Number of Uint4 used      */
} FASTALookup, PNTR FASTALookupPtr;

typedef struct _FDB_options {
    CharPtr db_title;
    CharPtr db_file;
    CharPtr LogFileName;
    Int4 is_protein;
    Int4 parse_mode;
    Int4 isASN;
    Int4 asnbin;
    Int4 is_seqentry;
    CharPtr base_name;
    Int4  dump_info;
} FDB_options, PNTR FDB_optionsPtr;

typedef struct formatdb 
{
    /* CharPtr	dbname;	(db_file)  name of input database */
    /* CharPtr	basename; (base_name)  base-name to be used for 
       BLAST databases. */
    /* CharPtr	DbTitle; (db_title) database title */
    
    /* file handlers */

    FILE *fd,
        *fd_ind,
        *fd_seq,
        *fd_def,
        *fd_sdi,  /* This is file for misc. info data */
        *fd_stmp;
    
    /* ASN.1 input, if the "-a" specified */
    AsnIoPtr	aip;
    
    /* Boolean	isProtein;  is_protein true, if input database 
       contains proteins */

    /* Boolean ParseMode; parse_mode  true, if ISAM needed */

    Int4 num_of_seqs;  /* number of parsed sequences */
    Int4 TotalLen, MaxSeqLen;
    
    /* offset tables */
    Int4Ptr	DefOffsetTable,	/* definitions */
        	SeqOffsetTable,	/* sequences */
        	AmbOffsetTable;	/* ambiguities */

    /* lookup table */

    FASTALookupPtr	lookup;

    /* General formatdb options */

    FDB_optionsPtr options;

    Uint4Ptr	AmbCharPtr;	/* ambiguity characters while
                                 * convert from ncbi2na->ncbi4na */

    Int4 OffsetAllocated; /* storage for allocation size */

    
} FormatDB, *FormatDBPtr;

/* Function prototypes for formatdb library*/

FormatDBPtr	FormatDBInit(FDB_optionsPtr options);
Int2 FDBAddSequence (FormatDBPtr fdbp, Int4 gi, ValNodePtr seq_id,
                     CharPtr title, Int4 tax_id, CharPtr div, 
                     Int4 owner, Uint1 seq_data_type, 
                     ByteStorePtr *seq_data, Int4 SequenceLen, Int4 date);
Int2 FormatDBClose(FormatDBPtr fdbp);


Int2 process_sep (SeqEntryPtr sep, FormatDBPtr fdbp);

NLM_EXTERN Boolean SeqEntrysToBLAST (SeqEntryPtr sep, FormatDBPtr fdbp,
                                     Boolean is_na, Uint1 group_segs);

NLM_EXTERN Boolean BLASTFileFunc (BioseqPtr bsp, Int2 key, CharPtr buf,
                                  Uint4 buflen, Pointer data);

/*
Print a summary of the database used.
*/
Boolean LIBCALL PrintDbInformation PROTO((CharPtr database, Boolean is_aa, Int4 line_length, FILE *outfp, Boolean html));

Boolean LIBCALL PrintDbInformationBasic PROTO((CharPtr database, Boolean is_aa, Int4 line_length, CharPtr definition, Int4 number_seqs, Int8 total_length, FILE *outfp, Boolean html));


#ifdef __cplusplus
}
#endif

#endif /* _READDB_ */
