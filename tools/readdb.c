
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

File name: readdb.c

Author: Tom Madden

Contents: Reads Databases formatted by formatdb.

Detailed Contents:

        - memory maps files. 

	- database sequences are identified (by these routines) by their
	order in the files.  this is based on a zero-offset.


******************************************************************************/

/* File Name: readdb.c
*
* Author: Tom Madden
*
* Version Creation Date:   3/22/95
*
* $Revision: 6.167 $
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
* $Log: readdb.c,v $
* Revision 6.167  2000/05/25 20:31:24  madden
* Do not change aliasfilebit unless it is zero
*
* Revision 6.166  2000/05/23 21:22:37  dondosha
* Do not open sequence files in shared_info when flag is set to not do it - correction to previous change
*
* Revision 6.165  2000/05/22 18:46:43  dondosha
* Merged all Boolean members in ReadDBFILE structure into a single Int4
*
* Revision 6.164  2000/05/09 15:54:19  shavirin
* Added function ReadDBBioseqSetDbGeneticCode().
*
* Revision 6.163  2000/05/03 17:41:21  madden
* Fix for readdb_get_descriptor problem when searching subset database
*
* Revision 6.162  2000/05/03 16:19:01  dondosha
* Added function FastaToBlastDB
*
* Revision 6.161  2000/05/03 12:49:45  madden
* Do not add > if not first definition
*
* Revision 6.160  2000/05/01 20:01:11  madden
* Protection against too large gis
*
* Revision 6.159  2000/04/19 17:59:23  madden
* Move setting of start and stop, adjust of indices to end, do every time in case of recursive calls or multiple databases
*
* Revision 6.158  2000/04/14 21:16:36  madden
* Fix for non-NULL aliasfilename
*
* Revision 6.157  2000/04/11 19:56:48  madden
* Set aliasfilename even if oidlist does not exist
*
* Revision 6.156  2000/04/10 18:01:46  dondosha
* Fixed FindBlastDBFile when file exists in current directory
*
* Revision 6.155  2000/04/05 19:25:09  madden
* Check for NULL searchstr in Fastacmd_Search, allow line break to be a valid delimiter for a file that is read in
*
* Revision 6.154  2000/04/03 21:17:57  dondosha
* readdb_MakeGiFileBinary will sort gis in increasing order
*
* Revision 6.153  2000/04/03 17:34:27  shavirin
* Fixed case when indexed and regular databases are mixed in multiple
* database set.
*
* Revision 6.152  2000/03/28 04:38:20  egorov
* Bug seen on Malaria page fixed
*
* Revision 6.151  2000/03/24 14:36:43  egorov
* Allow NULL alias_dbid on input of GI2OID
*
* Revision 6.150  2000/03/24 14:34:33  egorov
* Add support for month.sts, month.pataa, month.patnt month subsets
*
* Revision 6.149  2000/03/20 22:03:34  egorov
* bug with multiple alias databases mask is fixed
*
* Revision 6.148  2000/03/20 17:03:19  dondosha
* Return NULL from readdb_get_link and readdb_get_bioseq if cannot mem-map files
*
* Revision 6.147  2000/03/20 14:36:54  egorov
* Add protection from the reading out of the ISAM index file boundary when update CommonIndex.
*
* Revision 6.146  2000/03/16 19:47:18  egorov
* Db mask should be Uint4, not Int2.  Also previous change about FreeOIDList is rolled back.
*
* Revision 6.145  2000/03/16 18:09:50  dondosha
* Fixes memory leak in OIDListFree; corrects ReadDBCloseMHdrAndSeqFiles
*
* Revision 6.144  2000/03/15 21:34:30  egorov
* 1. Fix bug with using alias databases.
* 2. 2. Initialize new_defline variable.
*
* Revision 6.143  2000/03/13 18:36:37  madden
* Added insert_ctrlA Boolean to readdb_get_bioseq_ex
*
* Revision 6.142  2000/03/13 13:53:50  madden
* Check for non-NULL rdfp before dereference
*
* Revision 6.141  2000/03/10 19:16:30  shavirin
* Added multi-thread support for the function ReadDBBioseqFetchEnable().
*
* Revision 6.140  2000/03/10 18:51:33  madden
* Add prototype for readdb_get_filebits
*
* Revision 6.139  2000/03/08 22:03:32  madden
* added readdb_get_filebits
*
* Revision 6.138  2000/03/08 20:52:37  madden
* readdb_get_bioseq_ex only returns gis for subset database
*
* Revision 6.137  2000/02/28 21:50:13  egorov
* All month subsets use same criteria.
*
* Revision 6.136  2000/02/24 19:02:37  egorov
* Add support for PDB subset of nr
*
* Revision 6.135  2000/02/16 18:39:59  madden
* Fix check for nucl. alias file
*
* Revision 6.134  2000/02/11 19:59:29  shavirin
* Increased nthreads when attaching to the rdfp structure.
*
* Revision 6.133  2000/02/09 19:35:51  madden
* Added readdb_MakeGiFileBinary
*
* Revision 6.132  2000/02/07 21:15:15  madden
* Issue warning before stripping zero gi
*
* Revision 6.131  2000/02/07 20:56:08  madden
* Strip off gi|0 identifiers for formatdb
*
* Revision 6.130  2000/01/26 15:38:34  madden
* Fix for fastacmd and alias files
*
* Revision 6.129  2000/01/26 15:19:59  madden
* Return aliasfilename if present in readdb_get_filename
*
* Revision 6.128  2000/01/20 20:26:05  egorov
* Use "est_" prefix for subsets 'human', 'mouse', and 'others'
*
* Revision 6.127  2000/01/20 18:57:24  madden
* Check whether rdfp is NULL before dereference
*
* Revision 6.126  2000/01/12 21:51:50  madden
* Check for oidlist before setting aliasfilename
*
* Revision 6.125  2000/01/12 21:46:35  dondosha
* Fixed memory leak (rdfp->aliasfilename)
*
* Revision 6.124  2000/01/12 21:03:52  egorov
* 1. Introduce Fastacmd API function - Fastacmd_Search
* 2. Rearrange order of functions to have Fastacmd, ID1, and CommonIndex stuff separate.
*
* Revision 6.123  2000/01/12 20:28:31  dondosha
* Fixed readdb_new_ex2 behavior with multiple volume database
*
* Revision 6.122  2000/01/12 18:06:03  egorov
* Fix memory leak.  Remove debug stuff.
*
* Revision 6.121  2000/01/12 17:39:31  madden
* Fix readdb_parse_db_names so done is TRUE on last db
*
* Revision 6.120  2000/01/11 15:32:46  dondosha
* Fixed memory leaks in opening shared header and sequence file memory maps
*
* Revision 6.119  2000/01/07 16:00:25  madden
* Alias db length is Int8 instead of Uint4
*
* Revision 6.118  1999/12/31 14:23:20  egorov
* Add support for using mixture of real and maks database with gi-list files:
* 1. Change logic of creating rdfp list.
* 2. BlastGetDbChunk gets real databases first, then masks.
* 3. Propoper calculation of database sizes using alias files.
* 4. Change to CommonIndex to support using of mask databases.
* 5. Use correct gis in formated output (BlastGetAllowedGis()).
* 6. Other small changes
*
* Revision 6.117  1999/12/29 13:46:42  madden
* Fix for moving virtual rdfp to end, remove bad fix for infinite recursion
*
* Revision 6.116  1999/12/23 18:15:37  madden
* Move mask databases to end of all databases
*
* Revision 6.115  1999/12/22 21:54:41  dondosha
* Open header and sequence files consecutively as needed, close them when all threads have finished working with the database
*
* Revision 6.114  1999/12/21 20:02:16  egorov
* Set proper 'start' and 'stop' values for mask's rdfp.
* Add 'start' parameter into readdb_gi2seq.  This is return
* value which is set to rdfp->start where given gi was found.
*
* Revision 6.113  1999/12/17 21:33:01  egorov
* Add support for the 'month' subset.
*
* Revision 6.112  1999/12/17 20:47:05  egorov
* Fix 'gcc -Wall' warnings
*
* Revision 6.111  1999/12/15 21:57:58  egorov
* Initialize extra_bytes variable
*
* Revision 6.110  1999/12/15 17:40:07  egorov
* 1. Fix but with path to CommonIndexFile.
* 2. Add ScanDIFile() function for scanning DI index file and perform
*    callback-specified action for each record which meets database
*    subset criteria.
* 3. Change UpdateCommonIndexFile() function to use ScanDIFile.
* 4. Criteria for est_others, est_human, est_mouse, swissprot added.
*
* Revision 6.109  1999/12/14 19:27:09  dondosha
* Test against infinite recursion in IndexFileExists
*
* Revision 6.108  1999/11/30 17:07:15  egorov
* Fix problem with parsing database names when file_path is not NULL.
* Add prefix path prefix to OIDLIST and GILIST values, if any.
*
* Revision 6.107  1999/11/29 14:45:47  egorov
* Bug fixed.
*
* Revision 6.106  1999/11/24 21:43:34  madden
* Added Nlm_SwapUint4 call to make database masks work with both big and small endian systems
*
* Revision 6.105  1999/11/24 18:42:25  egorov
* It was reported by Andrei Shkeda and observed by us in neighboring software
* that using ReadDbFile structure was not MT-safe, and, as a result,
* it was impossible to format seqalign from different threads.
* Now it is fixed.  Mutecies shared same ISAM structures, so I had to put additional mutex.
*
* Revision 6.104  1999/11/24 18:01:38  egorov
* Bug fixed:  it truncated full database name to just file name if BLASTDB was specified.
* So it was impossible to have BLASTDB=/blast/db/blast and filename = "subdir/database".
* Now it works and makes it possible to use subdirectories for organism-specific databases.
*
* Revision 6.103  1999/11/23 22:02:26  madden
* Added readdb_get_totals_ex that may use alias file values
*
* Revision 6.102  1999/11/23 21:30:10  madden
* Deallocate OID list
*
* Revision 6.101  1999/11/22 16:15:36  egorov
* Remove correct return code in readdb_get_header function
*
* Revision 6.100  1999/11/15 17:42:48  egorov
* Fix bug when CommonIndex finds wrong Gi if database is not the first
* or the second in the CommonIndex list
*
* Revision 6.99  1999/11/12 14:15:54  madden
* Allow NlmOpenMFILE to simply open a file if it cannot be memory-mapped, alow other initialization states in readdb_new_ex2
*
* Revision 6.98  1999/10/07 20:40:48  madden
* Remove calls and function readdb_get_index
*
* Revision 6.97  1999/10/07 13:40:37  madden
* remove extra call to Nlm_SwapUint4
*
* Revision 6.96  1999/10/06 21:08:36  shavirin
* Cleared last bits in last byte written in function FDBAddSequence()
* These bits may be dirty in case of ASN.1 coming directly from ID.
*
* Revision 6.95  1999/10/01 18:25:07  shavirin
* Fixed bug in the function FDBAddSequence
*
* Revision 6.94  1999/09/30 20:48:24  madden
* Change static buffer to dynamically allocated
*
* Revision 6.93  1999/09/29 17:20:34  shavirin
* Fixed minor memory leak.
*
* Revision 6.92  1999/09/29 13:30:51  shavirin
* Changed sequence of allocating/deleting of oidlist structure.
*
* Revision 6.91  1999/09/28 20:45:07  shavirin
* Passed oidlist info when cloning rdfp in readdb_attach() function.
*
* Revision 6.90  1999/09/28 13:41:57  shavirin
* Freed memory of OID list in readdb_destruct().
*
* Revision 6.89  1999/09/24 16:30:25  egorov
* Remove Mac incompatible stuff.  Add two more functions for CommonIndex API.
*
* Revision 6.88  1999/09/23 18:22:30  egorov
* Do not keep private copy of index arrays (sequence_index, header_index, ambchar_index),
* but just use as it is in memory mapped file.  Big and small endian stuff is not forgot.
*
* Revision 6.87  1999/09/23 15:17:24  egorov
* Add CommonIndex API function - UpdateCommonIndexFile
*
* Revision 6.86  1999/09/23 15:10:52  egorov
* Add new fields into OIDList structure.
* Add new keywords into alias file: NSEQ and LENGTH.
* Use Nlm_Malloc instead of MemNew where MemSet is not needed.
* Create ReadOIDList function.
*
* Revision 6.85  1999/09/23 15:03:43  egorov
* Close alias file;  change name of index file;  add comments
*
* Revision 6.84  1999/09/22 21:58:07  egorov
* fix compilation bug
*
* Revision 6.83  1999/09/13 16:18:37  shavirin
* Added function readdb_get_bioseq_ex, which has possibility
* to bypass ObjMgr registration.
*
* Revision 6.82  1999/09/10 16:30:17  shavirin
* Fixed problems with formating proteins by formatdb
*
* Revision 6.81  1999/09/09 18:25:04  shavirin
* Added functions to parse ASN.1 with formatdb
*
* Revision 6.80  1999/09/02 18:02:33  madden
* No spaces after date
*
* Revision 6.79  1999/09/02 12:56:52  egorov
* Change format of the BLAST index file to set proper alignment
* for memory map.
*
* Revision 6.78  1999/08/30 18:21:29  shavirin
* Temporary return of full dumping set in SeqidE2Index() function.
*
* Revision 6.77  1999/08/26 20:55:55  shavirin
* Changed way to look for seqids.
*
* Revision 6.76  1999/08/26 14:12:50  shavirin
* Redused amount of information dumped for string indexes in regular case.
*
* Revision 6.75  1999/08/25 20:17:38  shavirin
* Added option to create and retrieve from sparse indexes.
*
* Revision 6.74  1999/08/04 18:26:41  madden
* Change databases in alias file for file path
*
* Revision 6.72  1999/08/03 19:21:44  shavirin
* Changed to dynamically allocated memory in function readdb_read_alias_file()
*
* Revision 6.71  1999/08/02 13:36:01  shavirin
* Rolled back last changes.
*
* Revision 6.69  1999/06/29 19:26:59  madden
* Took SeqIdWrite out of loop for efficiency
*
* Revision 6.68  1999/06/10 20:53:22  egorov
* Few changes to make it possible to perform multiple searches against different db's.
*
* Revision 6.67  1999/05/28 14:30:38  yaschenk
* rolling back fixes of 6.63 by shavirin, since they lead to coredump
*
* Revision 6.66  1999/05/27 21:47:12  yaschenk
* fix to the previous change
*
* Revision 6.65  1999/05/27 21:41:44  yaschenk
* dump_info file should be created fro nucleotides
*
* Revision 6.64  1999/05/27 15:51:29  shavirin
* Added function readdb_get_defline ()
*
* Revision 6.63  1999/05/27 14:40:17  shavirin
* Fixed some memory leaks.
*
* Revision 6.62  1999/05/21 17:36:52  madden
* Minor efficiencies
*
* Revision 6.61  1999/05/18 20:35:30  madden
* Changes to read an alias file for multiple db searches and ordinal ID lists
*
* Revision 6.60  1999/05/17 15:28:30  egorov
* First check that gi belongs to correct database and only then do all CommonIndex stuff
*
* Revision 6.59  1999/05/13 19:31:13  shavirin
* More changes toward dump from ID.
*
* Revision 6.58  1999/05/12 15:48:33  shavirin
* Many changes to fit new dump from ID.
*
* Revision 6.57  1999/05/10 13:47:44  madden
* NULL database not a fatal error
*
* Revision 6.56  1999/05/04 13:12:19  egorov
* Declare parse* functions as static and remove unused argument
*
* Revision 6.55  1999/05/03 21:44:33  chappey
* getline is now static function
*
* Revision 6.54  1999/04/27 17:28:17  shavirin
* Fixed few problems in the function FDBAddSequence().
*
* Revision 6.53  1999/04/26 14:55:23  shavirin
* Checked variable for not NULL.
*
* Revision 6.52  1999/04/26 14:36:04  shavirin
* Added ability to dump statistics.
*
* Revision 6.50  1999/04/21 22:59:41  kans
* added includes
*
* Revision 6.49  1999/04/21 21:43:28  shavirin
* Added set of functions, which used in "formatdb".
*
* Revision 6.48  1999/04/14 14:53:49  madden
* Correction for databases over 2 Gig
*
* Revision 6.47  1999/03/23 14:38:28  egorov
* Destruct CommonIndex structures only by thread it belongs to.
*
* Revision 6.46  1999/03/19 19:29:47  egorov
* Bug fixed.  Initialize cih.
*
* Revision 6.45  1999/03/18 16:55:22  egorov
* Previous fix was incompleete.
*
* Revision 6.44  1999/03/18 16:36:16  egorov
* Check if rdfp is not NULL before dereferencing it.
*
* Revision 6.43  1999/03/17 16:57:21  egorov
* Previously each element in rdfp list had his own CommonIndexHeadPtr
* initialized with MemMap.  But when we do search agains many databases,
* like in case of unfinished genomes, we meet limit for doing MemMap
* on SGI machines.  So now we initialize 'rdfp->cih' only for the first
* element in the list and reuse it for the others.
* Also the change contains proper freeing memory after the above change.
*
* Revision 6.42  1999/03/12 23:02:49  madden
* initialize memory in buffer_2na first
*
* Revision 6.41  1999/03/12 18:36:16  madden
* formatting fix
*
* Revision 6.40  1999/02/22 21:49:08  egorov
* Optimize GIs2OIDs using already initialized ISAM indecies from rdfp.  Use SwapUint4 function to use common index file
* on Solaris/Intel machines
*
* Revision 6.39  1999/02/18 21:19:12  madden
* ignore GIs not in common index
*
* Revision 6.38  1999/02/17 13:23:40  madden
* use MapNa2ByteToNa4String
*
* Revision 6.37  1999/01/07 14:35:01  madden
* Fix for readdb_acc2fasta for multiple databases
*
* Revision 6.36  1998/12/14 21:50:15  egorov
* new max gi number memeber in CommonIndexHead structure and therefore no need for COMMON_INDEX_TABLE_SIZE
*
* Revision 6.35  1998/09/24 15:26:41  egorov
* Fix lint complaints
*
* Revision 6.34  1998/09/14 15:11:20  egorov
* Add support for Int8 length databases; remove unused variables
*
* Revision 6.33  1998/09/03 18:43:09  egorov
* Close db config file
*
* Revision 6.32  1998/08/29 20:05:47  madden
* Fixed MemCpy length problem
*
* Revision 6.31  1998/08/24 14:59:56  madden
* readdb_get_sequence_ex function
*
* Revision 6.30  1998/07/31 19:30:11  egorov
* Fix bug when OID=0 treated as bad in common index
*
* Revision 6.29  1998/07/09 13:35:16  egorov
* remove platform dependent statement
*
* Revision 6.28  1998/07/08 14:10:53  madden
* Fix for multiple db search, use of more efficient readdb_new_ex
*
* Revision 6.27  1998/07/01 16:45:25  egorov
* Remove debug mesages
*
* Revision 6.26  1998/07/01 14:14:49  egorov
* Move FilePathFind function into ncbitoolkit remove its definition here
*
* Revision 6.25  1998/07/01 14:03:04  egorov
* Fix bug with a thread freeing CommonIndex: add new flag to rdfp
*
* Revision 6.24  1998/06/26 16:51:13  egorov
* Fix CommonIndex bugs
*
* Revision 6.23  1998/06/24 21:03:35  egorov
* Remove memory leaks
*
* Revision 6.20  1998/05/22 20:19:53  madden
* Changes to fix multi-db search bug
*
* Revision 6.19  1998/02/26 22:49:23  kans
* needed to include ffprint.h
*
* Revision 6.18  1998/02/26 22:34:21  madden
* Changes for 16 bit windows
*
* Revision 6.17  1998/01/16 22:02:03  madden
* Added readdb_new_ex with init_indices Boolean to allow faster retrieval of one sequence
*
* Revision 6.16  1997/12/12 20:39:25  madden
* Added parens for if
*
* Revision 6.15  1997/12/11 22:21:05  madden
* Removed unused variables
*
* Revision 6.14  1997/12/03 21:48:01  madden
* Check for duplicate database names
*
* Revision 6.13  1997/12/02 22:18:09  madden
* Fixed UMR
*
* Revision 6.12  1997/11/26 22:48:35  madden
* Added readdb_parse_db_names for multiple db searches
*
* Revision 6.11  1997/11/07 16:16:14  shavirin
* Added new function readdb_acc2fastaEx(), that retrieve array of hits
*
* Revision 6.10  1997/11/07 14:44:53  madden
* Sped up start up
*
* Revision 6.9  1997/11/06 21:27:19  madden
* Speeded up initialization
*
* Revision 6.8  1997/10/30 18:16:12  madden
* Change to readdb_acc2fasta to allow lookups by accession strings
*
* Revision 6.7  1997/10/24 19:08:13  madden
* Added ReadDBGetDb and ReadDBGetDbId
*
* Revision 6.6  1997/10/24 14:10:30  madden
* Changed Fetch function to speed up retrieval of cached sequences
*
* Revision 6.5  1997/09/24 22:37:03  madden
* Added readdb_destruct_element
*
* Revision 6.4  1997/09/16 16:31:36  madden
* More changes for multiple db runs
*
* Revision 6.3  1997/09/12 19:55:35  madden
* Added readdb_compare
*
* Revision 6.2  1997/09/11 18:49:37  madden
* Changes to enable searches against multiple databases.
*
* Revision 6.1  1997/08/27 14:46:56  madden
* Changes to enable multiple DB searches
*
* Revision 6.0  1997/08/25 18:53:55  madden
* Revision changed to 6.0
*
* Revision 1.52  1997/07/14 20:11:21  madden
* Removed unused variables
*
* Revision 1.51  1997/06/26 20:32:55  madden
* Only convert sequence if ambig. chars
*
* Revision 1.50  1997/05/20 14:33:32  shavirin
* Fixed retrievel by LOCUS in function readdb_acc2fasta()
*
* Revision 1.49  1997/05/19 21:14:56  shavirin
* Changed function readdb_acc2fasta() as required by E2Index() functions
* family
*
* Revision 1.48  1997/05/16 13:50:42  madden
* Fixed bug, wrong type of database opened
*
* Revision 1.47  1997/05/12 21:33:57  madden
* readdb_new allows indeterminate database type
*
* Revision 1.46  1997/05/12 21:10:31  shavirin
* Added new function readdb_acc2fasta()
*
* Revision 1.44  1997/05/07 21:03:11  madden
* Added function SeqId2OrdinalId
*
* Revision 1.43  1997/05/01 17:27:31  shavirin
* Added new function readdb_seqid2fasta()
*
 * Revision 1.42  1997/03/31  17:06:40  shavirin
 * Changed function readdb_get_bioseq to use BSRebuildDNA_4na()
 * function.
 *
 * Revision 1.41  1997/03/26  14:01:34  madden
 * Changes to Fetch function to allow cached-out structures to be read back in.
 *
 * Revision 1.40  1997/03/05  18:24:17  madden
 * Fixed MT problem introduced with use of ISAM code.
 *
 * Revision 1.39  1997/02/26  23:39:54  madden
 * Removed unused variables.
 *
 * Revision 1.38  1997/02/26  20:37:31  madden
 * Added protection against MT use to fetch function.
 *
 * Revision 1.37  1997/02/25  23:52:05  madden
 * Added readdb_gi2seq call to ReadDBBioseqFetchFunc.
 *
 * Revision 1.36  1997/02/25  22:15:33  shavirin
 * Changes in accordance to ISAM API changes
 *
 * Revision 1.35  1997/02/25  16:28:05  shavirin
 * Added function readdb_gi2seq() - returnes sequence number from gi
 *
 * Revision 1.34  1997/02/14  17:17:59  madden
 * Checked for NULL return from MemNew.
 *
 * Revision 1.33  1997/02/07  22:32:40  madden
 * Fixed bug.
 *
 * Revision 1.32  1997/01/14  23:11:27  madden
 * Cleaned ctrl-A's out of defline in readdb_get_bioseq.
 *
 * Revision 1.31  1996/12/20  00:30:20  madden
 * Protected ambiguity data against big/little endian changes.
 *
 * Revision 1.30  1996/12/19  16:29:56  madden
 * Changes to eliminate ".nac" file for nucl.
 *
 * Revision 1.29  1996/12/17  21:34:46  madden
 * Changes to allow deflines for inidividual entries to be retrieved.
 *
 * Revision 1.28  1996/12/11  18:42:36  madden
 * Added BioseqFetch functions.
 *
 * Revision 1.27  1996/12/11  17:59:42  madden
 * Fixed purify leaks.
 *
 * Revision 1.26  1996/12/08  15:19:59  madden
 * Checked for NULL pointer.
 *
 * Revision 1.25  1996/11/27  16:39:11  madden
 * Added functions to return filename and date.
 *
 * Revision 1.24  1996/11/26  19:54:27  madden
 * Added check for database in standard places.
 *
 * Revision 1.23  1996/11/22  19:05:48  madden
 * removed ifdef for OLD_BIT_ORDER.
 *
 * Revision 1.22  1996/11/18  17:28:13  madden
 * properly set contents_allocated flag for ambig. char. in readdb_attach.
 *
 * Revision 1.21  1996/11/08  21:45:03  madden
 * Removed function readdb_get_partial_unpacked_sequence.
 *
 * Revision 1.20  1996/11/07  22:31:15  madden
 * Added function readdb_ambchar_present to check for the presence
 * of ambig. characters in a db sequence.
 *
 * Revision 1.19  1996/11/04  18:48:53  shavirin
 * Added possibility to reconstruct Nucleotide sequence using function
 * readdb_get_bioseq. Added new function readdb_get_ambchar() to retrieve
 * ambiguity information.
 *
 * Revision 1.18  1996/10/31  16:29:18  shavirin
 * Multiple changes due to reverce of residues in BLAST database
 * for nucleotide sequences from (4321) to (1234)
 * New dumper now required to create BLAST databases.
 *
 * Revision 1.17  1996/09/27  19:12:17  madden
 * Added function readdb_get_bioseq to obtain a BioseqPtr from the BLAST databases.
 *
 * Revision 1.16  1996/09/26  20:18:43  madden
 * Saved filename.
 *
 * Revision 1.15  1996/09/23  17:36:20  madden
 * Removed unused variable.
 *
 * Revision 1.14  1996/09/23  14:37:35  madden
 * Replaced CharPtr (for sequence) with Uint1Ptr.
 *
 * Revision 1.13  1996/09/20  21:58:14  madden
 * Changed CharPtr's to Uint1Ptr, got remainder length out of top order bits.
 *
 * Revision 1.12  1996/09/16  13:48:51  madden
 * Removed extra increment of counter in readdb_get_partial_unpacked_sequence.
 *
 * Revision 1.11  1996/09/15  17:35:48  madden
 * readdb_get_partial_unpacked_sequence now packages ncbi4na properly.
 *
 * Revision 1.10  1996/09/13  18:55:04  madden
 * Added function readdb_get_partial_unpacked_sequence.
 *
 * Revision 1.9  1996/09/11  21:31:11  shavirin
 * Added check for NULL from function Nlm_MemMapInit(name)
 *
 * Revision 1.8  1996/08/29  20:42:01  madden
 * memory mapping moved to the corelib (in ncbimem.[ch]).
 *
 * Revision 1.7  1996/08/23  15:32:02  shavirin
 * Fixed a lot of NT compiler warnings about type mismatch
 *
 * Revision 1.6  1996/08/21  21:25:25  madden
 * Changes for reading nt. db's.
 *
 * Revision 1.5  1996/08/14  14:31:28  madden
 * Added efficiencies in readdb_get_sequence_length.
 *
 * Revision 1.4  1996/08/13  22:04:36  madden
 * Changed readdb_get_sequence to report the uncompressed length of
 * a nucl. sequence.
 *
 * Revision 1.3  1996/08/08  21:39:48  madden
 * Added code to read in nucleotide databases.
 *
 * Revision 1.2  1996/08/07  18:32:05  madden
 * Moved define of MMAP_AVAIL from readdb.h to readdb.c
 *
 * Revision 1.1  1996/08/05  19:48:21  madden
 * Initial revision
 *
 * Revision 1.14  1996/08/02  14:20:06  madden
 * Added readdb_attach function.
 *
 * Revision 1.13  1996/07/31  13:09:17  madden
 * Changes for partial copy of ReadDB structure.
 *
 * Revision 1.12  1996/07/29  19:43:35  madden
 * Changes to make BLAST big/little endian independent.
 *
 * Revision 1.11  1996/07/25  20:45:20  madden
 * Change to arguments of readdb_get_sequence.
 *
 * Revision 1.10  1996/07/25  12:56:15  madden
 * readdb_get_sequence changed to allow for systems w/o mmap.
 *
 * Revision 1.9  1996/06/20  16:16:36  madden
 * Replaced int's with Int4's.
 *
 * Revision 1.8  1996/06/07  15:05:21  madden
 * MemCpy used instead of a while loop.
 *
 * Revision 1.7  1996/05/16  21:07:33  madden
 * Added protections against missing input files.
 *
 * Revision 1.6  1996/05/16  19:50:15  madden
 * Added documentation block.
 *
 * Revision 1.5  1996/04/22  21:41:13  madden
 * memory mapping added.
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


/* Description of conception:

 * BLAST uses the concept of a virtual database, meaning that
 * we may have a few databases searched together as one.  For
 * a virtual database BLAST numbers the sequence from zero to
 * the total number in all databases minus one (the numbers
 * are called ordinal ID's or OID's):
 * 
 * 0 <= OID < total in all db's
 * 
 * Readdb is aware of these virtual OID's and handles them properly.
 * The situation has grown rather confused as we also allow
 * the specification of a gilist (that determines a subset of the
 * sequences in the database to be searched) as well as the recent
 * addition of the 'mask' database (which specifies that a subset of
 * a real database is to be searched, e.g., est_human is a mask and est is
 * the real database).  To clarify the situation Alexey and I have written
 * down some rules that describe how virtual databases should be used.
 * 
 * 1.) Ordinal ID (OID) numbering is from zero to total-1, where total
 * is the total number of sequences in all databases in the virtual
 * database.
 * 
 * 2.) OID's of 'mask' databases refer to the real (i.e., underlying)
 * database - not the mask.
 * 
 * 3.) If a gilist is used, then one virtual mask is used for all
 * databases - regardless of whether any database being searched is real.
 * 
 * 4.) If there is a mixture of real and mask databases and no gilist
 * is being used, then the mask databases should go to the end of the
 * virtual database and one virtual mask will be created for this subsection
 * of the virtual database.  (readdb_new_ex2 will be changed to move
 * the 'mask' databases to the end).
*/ 

#include <readdb.h>
#include <ncbithr.h>
#include <ffprint.h>
#include <ncbisami.h>
#include <blast.h>
#include <ncbisort.h>
#include <tofasta.h>
#include <errno.h>

/* Used by fetch functions. */
#define READDB_BUF_SIZE 255
#define READDBBF_INIT 0
#define READDBBF_READY 1
#define READDBBF_DISABLE 2

typedef struct readdbbioseqfetch {
    struct readdbbioseqfetch PNTR next;
    Uint1 ReadDBFetchState;
    CharPtr dbname;	/* Name of the database. */
    Uint2 ctr;
    Boolean is_prot; /* Is it a protein or not. */
    ReadDBFILEPtr rdfp;
    Int4 db_genetic_code;
    TNlmThread	thread_id;
} ReadDBFetchStruct, PNTR ReadDBFetchStructPtr;

typedef struct readdbfetchuserdata {
    Int4 ordinal_number;	/* ordinal number of db sequence. */
    Int2 db_id;		/* database ID, for multiple databases. */
} ReadDBFetchUserData, PNTR ReadDBFetchUserDataPtr;

static Int2 LIBCALLBACK ReadDBBioseqFetchFunc PROTO((Pointer data));
static ReadDBFILEPtr ReadDBFILENew(void);

static TNlmMutex isamsearch_mutex;	/* Mutex to regulate using ISAM;
					 rdfp->isam is common for all threads */
static TNlmMutex hdrseq_mutex;

/* Common index global variables */
Boolean	isCommonIndex = TRUE;

/**************************************************************************
*
*	Functions to perform memory mapping.
*
*	If memory mapping is not available, then these functions should
*	default to normal FILE pointers.
*
*	This is allowed with "read-only files right now.
*
**************************************************************************/

/*
	Initialize the memory-mapping.
*/
NlmMFILEPtr LIBCALL 
NlmOpenMFILE (CharPtr name)

{
	NlmMFILEPtr mfp;

	if ((mfp=(NlmMFILEPtr) MemNew(sizeof(NlmMFILE))) == NULL)
		return NULL;

        /* Default is FALSE. */
        mfp->mfile_true = FALSE;

	mfp->mmp_begin = NULL;

	if (Nlm_MemMapAvailable() == TRUE)
	{	 /* IF mem-map fails, open as a regular file. */
                if((mfp->mem_mapp = Nlm_MemMapInit(name)) != NULL)
		{ /* copy this pointer to where it's convenient. */
			mfp->mmp_begin = mfp->mmp = (Uint1Ptr) mfp->mem_mapp->mmp_begin;
			if (mfp->mmp_begin != NULL) 
			{
				mfp->mfile_true = TRUE;
				mfp->mmp_end = mfp->mmp_begin + mfp->mem_mapp->file_size;
			}
		}
	}

	if (mfp->mmp_begin == NULL)
	{
		mfp->fp = FileOpen(name, "rb");
		if (mfp->fp == NULL)
		{
			mfp = (NlmMFILEPtr) MemFree(mfp);
			return NULL;
		}
	}

	/* contents have been allocated. */
	mfp->contents_allocated = TRUE;

	return mfp;

}	/* NlmOpenMFILE */

/*
  Open the shared sequence and header files for memory mapping, if this hasn't
  already been done; duplicate this in headerfp and sequencefp 
*/
Boolean ReadDBOpenMHdrAndSeqFiles(ReadDBFILEPtr rdfp)
{
   Char buffer[PATH_MAX];
   Boolean is_prot = (Boolean) (rdfp->parameters & READDB_IS_PROT);

   if (rdfp->shared_info == NULL) {
      if(!((Boolean)(rdfp->parameters & READDB_NO_SEQ_FILE))) {
	 sprintf(buffer, "%s.%csq", rdfp->full_filename, is_prot? 'p':'n');
	 if ((rdfp->sequencefp = NlmOpenMFILE(buffer)) == NULL) {
	    ErrPostEx(SEV_WARNING, 0, 0, "Unable to open %s", buffer);
	    rdfp = readdb_destruct(rdfp);
	    return FALSE;
	 } 
      }
      sprintf(buffer, "%s.%chr", rdfp->full_filename, is_prot? 'p':'n');
      if((rdfp->headerfp = NlmOpenMFILE(buffer)) == NULL) {
	 ErrPostEx(SEV_WARNING, 0, 0, "Unable to open %s", buffer);
	 rdfp = readdb_destruct(rdfp);
	 return FALSE;
      } 
      return TRUE;
   }

   if (rdfp->shared_info->sequencefp == NULL || 
       rdfp->shared_info->headerfp == NULL) 
      rdfp->shared_info->nthreads = 1;
   else /* Just attaching another thread, not opening a new memory map */
      rdfp->shared_info->nthreads++;

   if (!((Boolean)(rdfp->parameters & READDB_NO_SEQ_FILE)) &&  
       rdfp->shared_info->sequencefp == NULL) {
      sprintf(buffer, "%s.%csq", rdfp->full_filename, is_prot? 'p':'n');
      if((rdfp->shared_info->sequencefp = NlmOpenMFILE(buffer)) == NULL) {
	 ErrPostEx(SEV_WARNING, 0, 0, "Unable to open %s", buffer);
	 rdfp = readdb_destruct(rdfp);
	 return FALSE;
      } 
   }
   rdfp->sequencefp = 
      (NlmMFILEPtr) MemDup(rdfp->shared_info->sequencefp, 
			   sizeof(NlmMFILE));
   rdfp->sequencefp->contents_allocated = FALSE;

   if (rdfp->shared_info->headerfp == NULL) {
      sprintf(buffer, "%s.%chr", rdfp->full_filename, is_prot? 'p':'n');
      if((rdfp->shared_info->headerfp = NlmOpenMFILE(buffer)) == NULL) {
	 ErrPostEx(SEV_WARNING, 0, 0, "Unable to open %s", buffer);
	 rdfp = readdb_destruct(rdfp);
	 return FALSE;
      } 
   }
   rdfp->headerfp = (NlmMFILEPtr) MemDup(rdfp->shared_info->headerfp, 
					 sizeof(NlmMFILE));
   rdfp->headerfp->contents_allocated = FALSE;
   
   return TRUE;
}

ReadDBFILEPtr ReadDBCloseMHdrAndSeqFiles(ReadDBFILEPtr rdfp)
{
   ReadDBFILEPtr start = rdfp;

   while (rdfp != NULL) {
      rdfp->sequencefp = NlmCloseMFILE(rdfp->sequencefp);  
      rdfp->headerfp = NlmCloseMFILE(rdfp->headerfp);  
      if (rdfp->shared_info) {
	 rdfp->shared_info->sequencefp = 
	    NlmCloseMFILE(rdfp->shared_info->sequencefp);
	 rdfp->shared_info->headerfp = 
	    NlmCloseMFILE(rdfp->shared_info->headerfp);
	 rdfp->shared_info->nthreads = 0;
	 /*rdfp->shared_info = 
	   (ReadDBSharedInfoPtr) MemFree(rdfp->shared_info);*/
      }
      rdfp = rdfp->next; 
   }
   return start;
}

ReadDBFILEPtr ReadDBFreeSharedInfo(ReadDBFILEPtr rdfp)
{
   ReadDBFILEPtr start = rdfp; 

   while (rdfp != NULL) {
      if (rdfp->shared_info) 
	 rdfp->shared_info = 
	    (ReadDBSharedInfoPtr) MemFree(rdfp->shared_info);
      rdfp = rdfp->next;
   }
   return start;
}

/****************************************************************************
*
*	Undo the memory-mapping.
*
*****************************************************************************/
NlmMFILEPtr LIBCALL 
NlmCloseMFILE (NlmMFILEPtr mfp)

{
	if (mfp == NULL)
		return NULL;

	/* Have the contents been allocated, or is this just an attachemnt? */
	if (mfp->contents_allocated)
	{

	        if (mfp->mfile_true == TRUE)
       	 	{
			Nlm_MemMapFini(mfp->mem_mapp);
		}

		if (mfp->fp)
		{
			FileClose(mfp->fp);
			mfp->fp = NULL;
		}
	}

	mfp = (NlmMFILEPtr) MemFree(mfp);
	return mfp;

}	/* NlmCloseMFILE */

/***********************************************************************
*
*	Analogous to ANSI-C fread.
*
************************************************************************/
Int4 LIBCALL 
NlmReadMFILE (Uint1Ptr buffer, size_t size, Int4 nitems, NlmMFILEPtr mfp)

{
	register size_t	diff, len;

	if (mfp == NULL)
		return 0;
	
	if (mfp->mfile_true == TRUE)
	{
		len = size * nitems;
		diff = mfp->mmp_end - mfp->mmp;
		if (len > diff) 
		{
			nitems = diff / size;
			len = nitems * size;
		}
		MemCpy((VoidPtr) buffer, (VoidPtr) mfp->mmp, len);
		mfp->mmp += len;
		return nitems;
	}
	
	return FileRead(buffer, size, nitems, mfp->fp);

}	/* NlmReadMFILE */

/*
	Seeks to a point in the file, analogous to fseek.
*/
Int4 LIBCALL 
NlmSeekInMFILE (NlmMFILEPtr mfp, long offset, Int4 ptrname)

{
	Uint1Ptr cp;

	if (mfp->mfile_true == TRUE)
	{
		switch (ptrname) {
			case SEEK_SET: /* relative to beginning */
				cp = mfp->mmp_begin + offset;
				if (offset < 0 || cp >= mfp->mmp_end)
					return -1;
				mfp->mmp = cp;
				break;
			case SEEK_CUR: /* relative to current position */
				cp = mfp->mmp + offset;
				if (cp >= mfp->mmp_end || cp < mfp->mmp_begin)
					return -1;
				mfp->mmp = cp;
				break;
			case SEEK_END: /* relative to end of file */
				if (offset > 0 || mfp->mem_mapp->file_size < -offset)
					return -1;
				mfp->mmp = mfp->mmp_begin + (mfp->mem_mapp->file_size + offset);
				break;
			default:
				return -1;
		}
		return 0;
	}

	return (Int4) fseek(mfp->fp, offset, ptrname);

}	/* NlmSeekInMFILE */

/*
	What is the offset (in bytes) to the beginning of the file.
	Analog to ftell.
*/
Int4 LIBCALL 
NlmTellMFILE (NlmMFILEPtr mfp)

{
	if (mfp->mfile_true == TRUE)
	{
		return (mfp->mmp - mfp->mmp_begin);
	}
	else
	{
		return (Int4) ftell(mfp->fp);
	}

}	/* NlmTellMFILE */

static ReadDBFILEPtr ReadDBFILENew(void)
{
  ReadDBFILEPtr new_t; 
  
  new_t = (ReadDBFILEPtr) MemNew(sizeof(ReadDBFILE));
  return new_t;
}

/*
	Parses the databases names (if more than one) from
	'filenames' into buffer.  buffer should already be
	long enough and allocated.  The funciton should be
	repeatedly called until TRUE is returned.
*/
Boolean LIBCALL
readdb_parse_db_names (CharPtr PNTR filenames, CharPtr buffer)

{
	Boolean done = FALSE;

	while (**filenames == ' ')
	{
		(*filenames)++;
	}

	while (**filenames != NULLB)
	{
		if (**filenames == ' ')
		{
			*buffer = NULLB;
			break;
		}
		*buffer = **filenames;
		buffer++;
		(*filenames)++;
	}

	while (**filenames == ' ')
	{
		(*filenames)++;
	}

	if (**filenames == NULLB)
	{
		*buffer = NULLB;
		done = TRUE;
	}

	return done;
}


typedef struct _readdb_alias_file {
	CharPtr title,		/* title of the database. */
		dblist,		/* list of databases. */
		gilist,		/* a gilist to be used with the database. */
		oidlist;	/* an ordinal id list to be used with this database. */
	Int8	len;		/* length of the database */
	Uint4	nseq;		/* number of seqs of the database */
} ReadDBAlias, PNTR ReadDBAliasPtr;
/*
	This function reads the 'alias' file for the BLAST databases.
*/

static void ReadDBAliasFree(ReadDBAliasPtr rdbap)
{
    MemFree(rdbap->title);
    MemFree(rdbap->dblist);
    MemFree(rdbap->gilist);
    MemFree(rdbap->oidlist);
    MemFree(rdbap);
    return;
}

static ReadDBAliasPtr
readdb_read_alias_file(CharPtr filename)

{
    Boolean done=FALSE, first=TRUE;
    CharPtr buffer;
    CharPtr file_path, ptr;
    Char file_buffer[PATH_MAX], full_buffer[PATH_MAX];
    ReadDBAliasPtr rdbap;
    FILE *fp;
    Int4 buflen, buffer_length, total_length=PATH_MAX, length;

    if (filename == NULL || (buflen = FileLength(filename)) <= 0)
        return NULL;
    
    fp = FileOpen(filename, "r");
    if (fp == NULL)
        return NULL;

    file_path = Nlm_FilePathFind(filename);
    
    buffer = MemNew(buflen + 1);
    
    rdbap = (ReadDBAliasPtr) MemNew(sizeof(ReadDBAlias));

    while (Nlm_FileGets(buffer, buflen + 1, fp) != NULL) {
        if (buffer[0] == '#')  /* ignore comments. */
            continue;
        
        if (StringNCmp(buffer, "TITLE", 5) == 0) {
            ptr = buffer;
            ptr += 5;
            while (*ptr == ' ')
                ptr++;
            
            *Nlm_StrChr(ptr, '\n') = NULLB;
            if (*ptr != NULLB)
                rdbap->title = StringSave(ptr);
            
            continue;
        }
        
        if (StringNCmp(buffer, "DBLIST", 6) == 0) {
            ptr = buffer;
            ptr += 6;
            while (*ptr == ' ')
                ptr++;
            
            *Nlm_StrChr(ptr, '\n') = NULLB;
            if (*ptr != NULLB)
	    {
		if (file_path && *file_path != NULLB)
		{ /* Prepend file_path if it exists. */
			rdbap->dblist = MemNew(total_length*sizeof(Char));
			length=0;
			while (!done)
			{
				done = readdb_parse_db_names(&ptr, file_buffer);
				sprintf(full_buffer, "%s%c%s", file_path, DIRDELIMCHR, file_buffer);
				buffer_length = StringLen(full_buffer);
				if (buffer_length+length > total_length)
				{
					rdbap->dblist = Realloc(rdbap->dblist, 2*total_length);
					total_length *= 2;
				}		
				if (!first)
				{
					StringCpy(rdbap->dblist+length, " ");
					length++;
				}
				else
					first = FALSE;
				StringCpy(rdbap->dblist+length, full_buffer);
					length += buffer_length;
			}
			
		}
		else
                	rdbap->dblist = StringSave(ptr);
	    }
            
            continue;
        }
        
        if (StringNCmp(buffer, "GILIST", 6) == 0) {
            ptr = buffer;
            ptr += 6;
            while (*ptr == ' ')
                ptr++;
            
            *Nlm_StrChr(ptr, '\n') = NULLB;
            if (*ptr != NULLB) {
		if (file_path && StrCmp(file_path,"")) {
		    /* add directory prefix, if any */
		    sprintf(full_buffer, "%s%c%s", file_path, DIRDELIMCHR, ptr);
		    rdbap->gilist = StringSave(full_buffer);
		} else {
		    rdbap->gilist = StringSave(ptr);
		}
	    }
            
            continue;
        }
        
        if (StringNCmp(buffer, "OIDLIST", 7) == 0) {
            ptr = buffer;
            ptr += 7;
            while (*ptr == ' ')
                ptr++;
            
            *Nlm_StrChr(ptr, '\n') = NULLB;
            if (*ptr != NULLB) {
		if (file_path && StrCmp(file_path, "")) {
		    /* add directory prefix, if any */
		    sprintf(full_buffer, "%s%c%s", file_path, DIRDELIMCHR, ptr);
		    rdbap->oidlist = StringSave(full_buffer);
		} else {
		    rdbap->oidlist = StringSave(ptr);
		}
	    }
            
            continue;
        }
        if (StringNCmp(buffer, "LENGTH", 6) == 0) {
            ptr = buffer;
            ptr += 6;
            while (*ptr == ' ')
                ptr++;
            
            *Nlm_StrChr(ptr, '\n') = NULLB;
            if (*ptr != NULLB)
                rdbap->len = atol(ptr);
            
            continue;
        }
        if (StringNCmp(buffer, "NSEQ", 4) == 0) {
            ptr = buffer;
            ptr += 4;
            while (*ptr == ' ')
                ptr++;
            
            *Nlm_StrChr(ptr, '\n') = NULLB;
            if (*ptr != NULLB)
                rdbap->nseq = atol(ptr);
            
            continue;
        }
    }
    
    MemFree(file_path);
    MemFree(buffer);
    FileClose(fp);
    return rdbap;
}
/* Check if .?in file exists for specified database
   and assign proper is_prot to rdfp->is_prot */

static	Int2	IndexFileExists(CharPtr full_filename, ReadDBFILEPtr PNTR rdfpp, Boolean PNTR is_prot, Uint1 init_state) 
{
    Char	buffer[PATH_MAX];
    Int4	length = 0;
    ReadDBAliasPtr rdbap;
    ReadDBFILEPtr rdfp=NULL;
    
    /* Check for protein and nucl. alias files first. */
    if (*is_prot == READDB_DB_UNKNOWN || *is_prot == READDB_DB_IS_PROT) {
        
        sprintf(buffer, "%s.pal", full_filename);
	rdbap = readdb_read_alias_file(buffer);
	if (rdbap && StringCmp(full_filename, rdbap->dblist)) {
            rdfp = readdb_new_ex2(rdbap->dblist, TRUE, init_state, rdbap->oidlist);
	    rdfp->aliasfilename = StringSave(Nlm_FileNameFind(full_filename));
	    if (rdfp->cih)
	    	rdfp->aliasfilebit = DBShift(rdfp->cih->num_of_DBs, rdfp->cih->dbids, rdfp->aliasfilename, TRUE);
            *rdfpp = rdfp;
            /* replace standard title with new one. */
            if (rdbap->title) {
                if (rdfp->title) {
                    MemFree(rdfp->title);
                }
                rdfp->title = rdbap->title;
                rdbap->title = NULL;
		/* Length of the database is already calculated in alias file */
		if (rdbap->len) {
		    rdfp->aliaslen = rdbap->len;
		}
		if (rdbap->nseq) {
		    rdfp->aliasnseq = rdbap->nseq;
		}
                rdfp = rdfp->next;
                while (rdfp) {
                    rdfp->title = MemFree(rdfp->title);
                    rdfp = rdfp->next;
                }
            }
            ReadDBAliasFree(rdbap);
            return 1;
	}
    } 
    if ((*is_prot == READDB_DB_IS_NUC) || (rdfp == NULL && *is_prot == READDB_DB_UNKNOWN)) {
       sprintf(buffer, "%s.nal", full_filename);
       rdbap = readdb_read_alias_file(buffer);
	if (rdbap && StringCmp(full_filename, rdbap->dblist)) {
            rdfp = readdb_new_ex2(rdbap->dblist, FALSE, init_state, rdbap->oidlist);
	    rdfp->aliasfilename = StringSave(Nlm_FileNameFind(full_filename));
	    if (rdfp->cih && rdfp->aliasfilebit == 0)
	    	rdfp->aliasfilebit = DBShift(rdfp->cih->num_of_DBs, rdfp->cih->dbids, rdfp->aliasfilename, FALSE);
            *rdfpp = rdfp;
            /* replace standard title with new one. */
            if (rdbap->title) {
                if (rdfp->title) {
                    MemFree(rdfp->title);
                }
                rdfp->title = rdbap->title;
                rdbap->title = NULL;
		/* Length of the database is already calculated in alias file */
		if (rdbap->len) {
		    rdfp->aliaslen = rdbap->len;
		}
		if (rdbap->nseq) {
		    rdfp->aliasnseq = rdbap->nseq;
		}
                rdfp = rdfp->next;
                while (rdfp) {
                    rdfp->title = MemFree(rdfp->title);
                    rdfp = rdfp->next;
                }
            }
            ReadDBAliasFree(rdbap);
            return 1;
	}
    }
    
    /* Now check for standard index files. */
    if (*is_prot == READDB_DB_UNKNOWN || *is_prot == READDB_DB_IS_PROT) {
        
        sprintf(buffer, "%s.pin", full_filename);
        length = FileLength(buffer);
        if (length > 0) {
            *is_prot = TRUE;
        }
    }
    
    if (length <= 0 && *is_prot != READDB_DB_IS_PROT) {
        
        sprintf(buffer, "%s.nin", full_filename);
        length = FileLength(buffer);
        if (length > 0) {
            *is_prot = FALSE;
        }
    }    
    
    if (length > 0)
        return 0;
    else
        return -1;
}

CharPtr	FindBlastDBFile (CharPtr filename)
{

    CharPtr	buffer, buffer1;
    Int4	len;
    
    /* first check current directory */

    len = FileLength(filename);
    if (len) 
       return StringSave(filename);

    buffer  = MemNew(PATH_MAX);
    buffer1 = MemNew(PATH_MAX);

#ifdef OS_UNIX
    if (getenv("BLASTDB"))
	Nlm_GetAppParam ("NCBI", "BLAST", "BLASTDB", getenv("BLASTDB"), buffer, PATH_MAX);
    else
#endif
	Nlm_GetAppParam ("NCBI", "BLAST", "BLASTDB", NULL, buffer, PATH_MAX);

    sprintf(buffer1, "%s%s%s", buffer, DIRDELIMSTR, filename);

    /* see if the file is not empty */

    len = FileLength(buffer1);

    if (len) {
	MemFree(buffer);
	return buffer1;
    } else {
    
	/* the only location where we now can find database file is a 
	   standard place #define'd as BLASTDB_DIR */

	sprintf(buffer, "%s%s%s", BLASTDB_DIR, DIRDELIMSTR, filename);
	len = FileLength(buffer);
	if (len) {
	    MemFree(buffer1);
	    return buffer;
	} else {
	    /* we cannot find directory :( */
 	    ErrPostEx(SEV_WARNING, 0, 0, "Could not find %s", filename);
	    return NULL;
	}
    }
}

/* 
	filename: name of the file to be openend.
	is_prot: three choices: protein, nucleotide, or either one.
	init_state: how much should be initialized.
		READDB_NEW_DO_ALL : initialize everything possible
		READDB_NEW_DO_REPORT : init enough for a report on db size etc.
	cih: common index
*/

static ReadDBFILEPtr
readdb_new_internal(CharPtr filename, Uint1 is_prot, Uint1 init_state, CommonIndexHeadPtr cih)
{
  ReadDBFILEPtr rdfp=NULL;
  Char buffer[PATH_MAX], buffer1[PATH_MAX];
  Char commonindex_full_filename[PATH_MAX];
  Char	database_dir[PATH_MAX] = "";
  Uint4 seq_type, formatdb_ver, date_length, title_length, value;
  Int2 status;
  Int4 length, num_seqs;
  CharPtr	charptr;
  ISAMErrorCode error;
  Boolean	localdb = FALSE;

  if (filename == NULL)
    return NULL;
  

      /* We need to find out what directory to use and which index system will
         be used for searching OID by give GI.  The algorithm is:
         Define blast database directory by present database searching
         in the following order:  current directory, .ncbirc, getenv("BLASTDB"),
         standard place;
         Then defind which index system to use.  If "CommonIndex" then
         we need all CommonIndex and ISAM files to be present in database directory,
         if ISAM, them the only ISAM index files should be present in the
         current directory
         */ 
  
 
      /* first see in the current directory */
  
  if ((status=IndexFileExists(filename, &rdfp, &is_prot, init_state)) >= 0)
  {

	if (status > 0)
		return rdfp;

          /* use current directory */
      charptr = Nlm_FilePathFind(filename);
      StringCpy(database_dir, charptr);
      MemFree(charptr);
      localdb = TRUE;
      rdfp = readdb_destruct(rdfp);
  }
  else 
  {
          /* Try to read this from .ncbirc file and if it is not specified in the file,
           then in BLASTDB environment variable (UNIX only) */
      
#ifdef OS_UNIX
      if (getenv("BLASTDB"))
	  Nlm_GetAppParam ("NCBI", "BLAST", "BLASTDB", getenv("BLASTDB"), buffer, PATH_MAX);
      else
#endif
      Nlm_GetAppParam ("NCBI", "BLAST", "BLASTDB", NULL, buffer, PATH_MAX);

      sprintf(buffer1, "%s%s%s", buffer, DIRDELIMSTR, filename);

      if ((status=IndexFileExists(buffer1, &rdfp, &is_prot, init_state)) >= 0) {
	  if (status > 0)
		return rdfp;
              /* database file is in directory 'buffer' */
          StringCpy(database_dir, buffer);
	  rdfp = readdb_destruct(rdfp);
      }
      else {
              /* the only location where we now can find database file is standard place
                 #define'd as BLASTDB_DIR */
          sprintf(buffer, "%s%s%s", BLASTDB_DIR, DIRDELIMSTR, filename);
          if ((status=IndexFileExists(buffer, &rdfp, &is_prot, init_state)) >= 0) {
	      if (status > 0)
		return rdfp;
              StringCpy(database_dir, BLASTDB_DIR);
	      rdfp = readdb_destruct(rdfp);
          }
          else {
                  /* we cannot find directory :( */
              ErrPostEx(SEV_WARNING, 0, 0, "Could not find index files for database %s", filename);
	      rdfp = readdb_destruct(rdfp); 
	      return rdfp;
          }
      }
  }

  rdfp = readdb_destruct(rdfp);
  rdfp = ReadDBFILENew();
  if (rdfp == NULL)
	return NULL;

  rdfp->filename = StringSave(filename);

  /*rdfp->is_prot = is_prot;*/
  if (is_prot)
     rdfp->parameters |= READDB_IS_PROT;
  

      /* Here we know that database is in database_dir directory */

      /* constract full file name */
  if (!StringCmp(database_dir, "")) {
      sprintf(rdfp->full_filename, "%s", Nlm_FileNameFind(filename));
  }
  else if (!localdb) {
      sprintf(rdfp->full_filename, "%s%s%s", database_dir, DIRDELIMSTR, filename);
  } else  {
      sprintf(rdfp->full_filename, "%s%s%s", database_dir, DIRDELIMSTR, Nlm_FileNameFind(filename));
  }

      /* Now let's find out which index system to use */

      /* First see if user has preferences */
  StringCpy(buffer1, "CommonIndex");
  
#ifdef OS_UNIX
  if (getenv("INDEX_SYSTEM") && StringCmp(getenv("INDEX_SYSTEM"), "CommonIndex"))
      StringCpy(buffer1, "ISAM");
#endif
  
  Nlm_GetAppParam ("NCBI", "BLAST", "INDEX_SYSTEM", buffer1, buffer, PATH_MAX);
  
  isCommonIndex = !StrCmp("CommonIndex", buffer);

      /* now we know that if isCommonIndex == TRUE, than it is prefered to use CommonIndex */

      /* test if there exist common index file */
  if (isCommonIndex) {
      if (!StringCmp(database_dir, "")) {
          sprintf(commonindex_full_filename, "%s", COMMONINDEX_FN);
      }
      else {
          sprintf(commonindex_full_filename, "%s%s%s", database_dir, DIRDELIMSTR, COMMONINDEX_FN);
      }

      if (!(length = FileLength(commonindex_full_filename))) {
              /* no CommonIndex files in this directory, try to use ISAM only */
          isCommonIndex = FALSE;
      }
  }


      /* check if present main three files: index, sequences, headers */

  sprintf(buffer, "%s.%cin", rdfp->full_filename, is_prot? 'p':'n');
  if((rdfp->indexfp = NlmOpenMFILE(buffer)) == NULL) {
    ErrPostEx(SEV_WARNING, 0, 0, "Unable to open %s", buffer);
    rdfp = readdb_destruct(rdfp);
    return rdfp;
  }
  
  if (init_state & READDB_NEW_DO_ALL)
     if (ReadDBOpenMHdrAndSeqFiles(rdfp) == FALSE)
	ErrPostEx(SEV_ERROR, 0, 0, 
		  "ReadDBOpenMHdrAndSeqFiles: failed to map files\n");

      /* fill in other fields of rdfp-> */
  NlmReadMFILE((Uint1Ptr) &value, 4, 1, rdfp->indexfp);
  formatdb_ver = Nlm_SwapUint4(value);
  if (formatdb_ver != FORMATDB_VER) {
    ErrPostEx(SEV_WARNING, 0, 0, "readdb: wrong version of formatdb "
              "was used to make database.");
    rdfp = readdb_destruct(rdfp);
    return NULL;
  }
  rdfp->formatdb_ver = formatdb_ver;
  
  NlmReadMFILE((Uint1Ptr) &value, 4, 1, rdfp->indexfp);
  seq_type = Nlm_SwapUint4(value);
  if ((is_prot && seq_type == 0) || (!is_prot && seq_type == 1)) {
    rdfp = readdb_destruct(rdfp);
    return rdfp;
  }
  NlmReadMFILE((Uint1Ptr) &value, 4, 1, rdfp->indexfp);
  title_length = Nlm_SwapUint4(value);

  if (title_length) {
    rdfp->title = (CharPtr)Nlm_Malloc((title_length+1)*sizeof(Char));
    NlmReadMFILE((Uint1Ptr) rdfp->title, title_length, 1, rdfp->indexfp);
    rdfp->title[title_length] = NULLB;
  } else {	/* Use the filename, if there is no title. */
    rdfp->title = StringSave(rdfp->filename);;
  }

  NlmReadMFILE((Uint1Ptr) &value, 4, 1, rdfp->indexfp);
  date_length = Nlm_SwapUint4(value);

  rdfp->date = (CharPtr)Nlm_Malloc((date_length+1)*sizeof(Char));
  NlmReadMFILE((Uint1Ptr) rdfp->date, date_length, 1, rdfp->indexfp);
  rdfp->date[date_length] = NULLB;

  NlmReadMFILE((Uint1Ptr) &(value), 4, 1, rdfp->indexfp);
  num_seqs = rdfp->num_seqs = Nlm_SwapUint4(value);
  NlmReadMFILE((Uint1Ptr) &(value), 4, 1, rdfp->indexfp);
  rdfp->totlen = Nlm_SwapUint4(value);
  NlmReadMFILE((Uint1Ptr) &(value), 4, 1, rdfp->indexfp);
  rdfp->maxlen = Nlm_SwapUint4(value);

  if (init_state & READDB_NEW_DO_REPORT) {
     rdfp->parameters |= READDB_CONTENTS_ALLOCATED;
     /*rdfp->contents_allocated = TRUE; */ /* Some was allocated, but index pointers are NULLs - that's OK */
      return rdfp;
  }

  if (!((title_length + date_length)%4) && rdfp->indexfp->mfile_true) {
      rdfp->header_index = (Uint4Ptr) rdfp->indexfp->mmp;
      rdfp->indexfp->mmp += 4 * (num_seqs+1);

      rdfp->sequence_index = (Uint4Ptr) rdfp->indexfp->mmp;
      rdfp->indexfp->mmp += 4 * (num_seqs+1);

      rdfp->ambchar_index = (Uint4Ptr) rdfp->indexfp->mmp;
      rdfp->indexfp->mmp += 4 * (num_seqs+1);
  } else {
      /* Use old stuff */
   
      if((rdfp->header_index = 
          (Uint4Ptr) Nlm_Malloc((num_seqs+1)*sizeof(Uint4))) == NULL) {
          rdfp = readdb_destruct(rdfp);
          return rdfp;
      }
      
      rdfp->header_index_start = rdfp->header_index;
      rdfp->header_index_offset = NlmTellMFILE(rdfp->indexfp);
      NlmReadMFILE((Uint1Ptr) rdfp->header_index, 4, num_seqs+1, rdfp->indexfp);
      
      if((rdfp->sequence_index = 
          (Uint4Ptr)Nlm_Malloc((num_seqs+1)*sizeof(Uint4))) == NULL) {
          rdfp = readdb_destruct(rdfp);
          return rdfp;
      }
      rdfp->sequence_index_start = rdfp->sequence_index;
      NlmReadMFILE((Uint1Ptr) rdfp->sequence_index, 4, num_seqs+1, rdfp->indexfp);
      
          /* For nucleotide sequence we will process ambiguity file */
      if(!is_prot) {
          if((rdfp->ambchar_index = (Uint4Ptr)Nlm_Malloc((num_seqs+1)*sizeof(Uint4))) == NULL) {
              rdfp = readdb_destruct(rdfp);
              return rdfp;
          }
          rdfp->ambchar_index_start = rdfp->ambchar_index;
          NlmReadMFILE((Uint1Ptr) rdfp->ambchar_index, 4, num_seqs+1, rdfp->indexfp);
      }
  }
  

  /* Contents were allocated above. */
  /*rdfp->contents_allocated = TRUE;*/
  rdfp->parameters |= READDB_CONTENTS_ALLOCATED;
  
  /* mmap is not being used, allocate a buffer 2 longer (for sentinel bytes)
     than the longest subject length. */
  if (rdfp->sequencefp && rdfp->sequencefp->mfile_true == FALSE) {
    rdfp->buffer = (UcharPtr)Nlm_Malloc((2+rdfp->maxlen)*sizeof(Uint1));
    if (rdfp->buffer == NULL)
    {
	rdfp = readdb_destruct(rdfp);
    	return rdfp;
    }
    rdfp->allocated_length = 2 + rdfp->maxlen;
  }

  /* Now initializing Numeric ISAM indexes */ 

  sprintf(buffer,  "%s.%cnd", rdfp->full_filename, is_prot? 'p':'n');  
  sprintf(buffer1, "%s.%cni", rdfp->full_filename, is_prot? 'p':'n');

  if(FileLength(buffer) != 0 && FileLength(buffer1) != 0) {
      if((rdfp->nisam_opt = ISAMObjectNew(ISAMNumeric, 
                                          buffer, buffer1)) == NULL) {
          ErrPostEx(SEV_WARNING, 0, 0, "Failed to create NISAM object");
          rdfp = readdb_destruct(rdfp);
          return rdfp;
      }
  }

  /* Now initializing Numeric ISAM indexes */ 

  sprintf(buffer,  "%s.%csd", rdfp->full_filename, is_prot? 'p':'n');  
  sprintf(buffer1, "%s.%csi", rdfp->full_filename, is_prot? 'p':'n');
  
  if(FileLength(buffer) != 0 && FileLength(buffer1) != 0) {

      if((rdfp->sisam_opt = ISAMObjectNew(ISAMString, 
                                          buffer, buffer1)) == NULL) {
          ErrPostEx(SEV_WARNING, 0, 0, "Failed to create SISAM object");
          rdfp = readdb_destruct(rdfp);
          return rdfp;
      }
      
#if 0
      if((error = ISAMGetIdxOption(rdfp->sisam_opt, &rdfp->sparse_idx)) < 0) {
          ErrPostEx(SEV_WARNING, 0, 0, "Failed to access string index "
                    "ISAM Error code is %d\n", error);
          rdfp = readdb_destruct(rdfp);
          return rdfp;
      }
#else
      rdfp->sparse_idx = ((ISAMDataPtr) rdfp->sisam_opt)->idx_option;
#endif
  }
  
  /* Now initializing Common index files */

  if (isCommonIndex) {
      if (cih) {
	  rdfp->cih = cih;
	  /*rdfp->handle_common_index = FALSE;*/
	  rdfp->parameters &= ~READDB_HANDLE_COMMON_INDEX;
      } else {
	  rdfp->cih = CommonIndexInit(commonindex_full_filename);
	  /*rdfp->handle_common_index = TRUE;*/
	  rdfp->parameters |= READDB_HANDLE_COMMON_INDEX;
      }

      if (!(rdfp->cih)) {
          isCommonIndex = FALSE;
          /*rdfp->handle_common_index = FALSE;*/
	  rdfp->parameters &= ~READDB_HANDLE_COMMON_INDEX;
      } else {
          rdfp->filebit = DBShift(rdfp->cih->num_of_DBs, rdfp->cih->dbids,
                                  Nlm_FileNameFind(filename), is_prot);
      }
  }

  /* Initialize shared information structure */
  rdfp->shared_info = (ReadDBSharedInfoPtr) 
     Nlm_Malloc(sizeof(ReadDBSharedInfo));
  rdfp->shared_info->nthreads = 0;
  rdfp->shared_info->headerfp = NULL;
  rdfp->shared_info->sequencefp = NULL;

  return rdfp;
}

OIDListPtr OIDListFree (OIDListPtr oidlist)

{
	if (oidlist == NULL)
		return NULL;

	if (oidlist->memory)
		MemFree(oidlist->memory);
	else
		NlmCloseMFILE(oidlist->mfp);
	
	MemFree(oidlist);

	return NULL;
}


Boolean	ReadOIDList (OIDListPtr oidlist)
{
    NlmMFILEPtr       mmfile;
    Int4 length;
	


    /* format of the file:
       <number of OID's - N>
       <oid1>
       <oid2>
       ...
       <oidN>
     */

    /* use memmap */
    mmfile = NlmOpenMFILE(oidlist->filename);
    if (!mmfile) {
	ErrPostEx(SEV_ERROR, 0, 0, "Could not open OID mask.\n");
	return FALSE;
    }

    if (mmfile->mfile_true == FALSE)
    {
	length = FileLength(oidlist->filename);
	oidlist->memory = MemNew(length);
	if (oidlist->memory == NULL)
		return FALSE;
	FileRead(oidlist->memory, length, 1, mmfile->fp);
	oidlist->list = oidlist->memory + 1;
	NlmCloseMFILE(mmfile);
    }
    else
    {
    	oidlist->list = (Uint4Ptr) mmfile->mmp_begin + 1;
	oidlist->mfp = mmfile;
    }
    oidlist->total = Nlm_SwapUint4(*((Int4Ptr) mmfile->mmp_begin));

    return TRUE;
}


ReadDBFILEPtr LIBCALL 
readdb_new_ex (CharPtr filename, Uint1 is_prot, Boolean init_indices)

{
	return readdb_new_ex2(filename, is_prot, READDB_NEW_INDEX, NULL);
}

ReadDBFILEPtr LIBCALL 
readdb_new_ex2 (CharPtr filename, Uint1 is_prot, Uint1 init_state, CharPtr oidlist)

{
    Boolean done = FALSE, duplicate_db;
    Char buffer[PATH_MAX], buffer_oidlist[PATH_MAX];
    Int4 start=0, old_start = 0;
    ReadDBFILEPtr new, tmp, var, var1, rdfp_w_oidlist, var2;
    CommonIndexHeadPtr	cih = NULL;
    
    new = NULL;
    rdfp_w_oidlist = NULL;
    buffer_oidlist[0] = NULLB;
    
    while (!done) {
        done = readdb_parse_db_names(&filename, buffer);
        if (*buffer == NULLB)
            break;
        if (oidlist) { /* NOTE: no account taken of duplicate databases?? */
            readdb_parse_db_names(&oidlist, buffer_oidlist);
            if (*buffer_oidlist == NULLB)
                break;
        }
        /* Look for duplicates of the database names. */
        duplicate_db = FALSE;
        var1 = new;
        while (var1) {
            if (StringCmp(readdb_get_filename(var1), buffer) == 0) {
                duplicate_db = TRUE;
                break;
            }
            var1 = var1->next;
        }
        if (duplicate_db)
            continue;
        
        /* 'continue' if return is NULL in case only one of many databases can't 
           be found.  Warning issued by readdb_new_internal. */
        if(!(tmp = readdb_new_internal(buffer, is_prot, init_state, cih)))
            
            continue;

        if (tmp->cih) {
            cih = tmp->cih;
        }
        
        if (tmp->oidlist) {
            /* Save these separately. */
            if (rdfp_w_oidlist == NULL) {
                rdfp_w_oidlist = tmp;
            } else {
                var = rdfp_w_oidlist;
                while (var->next)
                    var = var->next;
                var->next = tmp;
            }
            continue; /* do not execute rest of while loop. */
        } else if (/*tmp->not_first_time == FALSE*/!(tmp->parameters & READDB_NOT_FIRST_TIME)) {
	   /*tmp->not_first_time = TRUE;*/
	   tmp->parameters |= READDB_NOT_FIRST_TIME;
            if (buffer_oidlist[0] != NULLB) {
                
				/* read this OID list */
                tmp->oidlist = (OIDListPtr) MemNew (sizeof(OIDList));
                tmp->oidlist->filename = StringSave(buffer_oidlist);
                if (!ReadOIDList(tmp->oidlist)) {
                    return NULL;
                }
            }
            
/*
            tmp->start += start;
            tmp->stop = tmp->num_seqs-1+tmp->start;
            tmp->ambchar_index -= start;
            tmp->header_index -= start;
            tmp->sequence_index -= start;
            
            for(var2 = tmp; var2 != NULL; var2=var2->next) {
                start += var2->num_seqs;
            }
*/

            var = new;
            if (new == NULL) {
                new = tmp;
            } else {
                var = new;
                while(var->next)
                    var = var->next;
                var->next = tmp;
            }
        }
    }
    
    /* Attach the RDFP's with an OID. */
    if (rdfp_w_oidlist) {
        if (new == NULL) {
            new = rdfp_w_oidlist;
        } else {
            var = new;
            while(var->next)
                var = var->next;
            var->next = rdfp_w_oidlist;
        }
    }

    /* adjust all the RDFP's. */
    tmp = new;
    start = 0;
    while (tmp) {
	/* this may have been adjusted on a previous call to this function, readjust for indices. */
	old_start = tmp->start;
        tmp->start = start;
        tmp->stop = tmp->num_seqs-1+start;
        tmp->ambchar_index -= (start-old_start);
        tmp->header_index -= (start-old_start);
        tmp->sequence_index -= (start-old_start);
        
        start = tmp->stop+1; 
        tmp = tmp->next;
    }
    
    
    if (new)
       /*new->not_first_time = FALSE;*/
       new->parameters &= ~READDB_NOT_FIRST_TIME;
    return new;
}

ReadDBFILEPtr LIBCALL 
readdb_new (CharPtr filename, Uint1 is_prot)

{

	return readdb_new_ex(filename, is_prot, TRUE);
}


/*
	Get total length and number of sequences in multiple databases.
	if 'use_alias' is TRUE, values from the alias file will be used
	if non-zero.
*/

Boolean LIBCALL
readdb_get_totals_ex(ReadDBFILEPtr rdfp_list, Int8Ptr total_len, Int4Ptr total_num, Boolean use_alias)

{
	Int8 aliaslen=0;
	Int4 aliasnseq=0;

	*total_len = 0;
	*total_num = 0;

	if (rdfp_list == NULL || total_len == NULL || total_num == NULL)
		return FALSE;
	
	while (rdfp_list) {
	    if (use_alias && rdfp_list->aliasfilename) {
		if (rdfp_list->aliaslen)
			*total_len += rdfp_list->aliaslen;
		else
			*total_len += readdb_get_dblen(rdfp_list);

		if (rdfp_list->aliasnseq)
			*total_num += rdfp_list->aliasnseq;
		else
			*total_num += readdb_get_num_entries(rdfp_list);
	    } else {
		*total_len += readdb_get_dblen(rdfp_list);
		*total_num += readdb_get_num_entries(rdfp_list);
	    }
	    rdfp_list = rdfp_list->next;
	}

	return TRUE;
}

/*
	Get total length and number of sequences in multiple databases.
*/

Boolean LIBCALL
readdb_get_totals(ReadDBFILEPtr rdfp_list, Int8Ptr total_len, Int4Ptr total_num)

{
	return readdb_get_totals_ex(rdfp_list, total_len, total_num, FALSE);
}

/*
	Checks whether a ReadDBFILEPtr is the original, or just attaced.
	It does this by checking the rdfp->contents_allocated flag.
*/
Boolean LIBCALL 
readdb_copy (ReadDBFILEPtr rdfp)

{
	if (rdfp == NULL)
		return FALSE;

	/* if allocated, this is not a copy. */
	/*if (rdfp->contents_allocated)*/
	if (rdfp->parameters & READDB_CONTENTS_ALLOCATED)
		return FALSE;

	return TRUE;
}

/*
	Check whether two different ReadDBFILEPtr refer to the
	same database.

	If they are, then TRUE is returned.
*/
Boolean LIBCALL
readdb_compare(ReadDBFILEPtr rdfp1, ReadDBFILEPtr rdfp2)

{
	if (rdfp1 == NULL || rdfp2 == NULL)
		return FALSE;

	if (rdfp1 == rdfp2)
		return TRUE;

	/*if (rdfp1->is_prot != rdfp2->is_prot)*/
	if ((rdfp1->parameters & READDB_IS_PROT) != 
	    (rdfp2->parameters & READDB_IS_PROT))
		return FALSE;

	if (rdfp1->totlen != rdfp2->totlen)
		return FALSE;

	if (rdfp1->maxlen != rdfp2->maxlen)
		return FALSE;

	if (StringCmp(rdfp1->filename, rdfp2->filename) != 0)
		return FALSE;

	if (StringCmp(rdfp1->title, rdfp2->title) != 0)
		return FALSE;

	if (StringCmp(rdfp1->date, rdfp2->date) != 0)
		return FALSE;

	return TRUE;
}

/*
	Attach to an already open ReadDBFILEPtr.  Duplicate the 
	indexfp, sequencefp, and headerfp structures as the pointers
	there (i.e., mmp) will need to be manipulated.  Do not 
	change the FILE PNTR fp.
*/

ReadDBFILEPtr LIBCALL 
readdb_attach (ReadDBFILEPtr rdfp)

{
        ReadDBFILEPtr head, last, new_t;

	if (rdfp == NULL)
		return NULL;

	head = NULL;
	last = NULL;
	while (rdfp)
	{
		new_t = (ReadDBFILEPtr) MemDup(rdfp, sizeof(ReadDBFILE));

		/* 
		The contents_allocated flag DOES NOT apply to the actual
		structures indexfp, headerfp, or sequencefp.  These must always 
		be duplicated, as their pointers need to be independently 
		manipulated by threads.  They have their own allocation flags.
		*/
       		/*new_t->contents_allocated = FALSE;*/
		new_t->parameters &= ~READDB_CONTENTS_ALLOCATED;
       		new_t->indexfp = (NlmMFILEPtr) MemDup(rdfp->indexfp, 
                                              sizeof(NlmMFILE));
		new_t->indexfp->contents_allocated = FALSE;
		if (rdfp->headerfp != NULL) {
		   new_t->headerfp = (NlmMFILEPtr) MemDup(rdfp->headerfp, 
							  sizeof(NlmMFILE));
		   new_t->headerfp->contents_allocated = FALSE;
		}
		if (rdfp->sequencefp != NULL) {
		   new_t->sequencefp = (NlmMFILEPtr) MemDup(rdfp->sequencefp, 
							    sizeof(NlmMFILE));
		   new_t->sequencefp->contents_allocated = FALSE;
		}
                /*new_t->handle_common_index = FALSE;*/
		new_t->parameters &= ~READDB_HANDLE_COMMON_INDEX;

                new_t->oidlist = rdfp->oidlist;

		/* Copy address of shared information */
		new_t->shared_info = rdfp->shared_info;
                if(new_t->shared_info != NULL) 
                    rdfp->shared_info->nthreads++;

		/* Contents_allocated also does not apply to buffer, this is
		determined by allocated_length. */
		if (new_t->allocated_length > 0)
        	{
               		 new_t->buffer = (UcharPtr) MemNew((new_t->allocated_length)*sizeof(Uint1));
        	}

		if (head == NULL)
		{
			head = new_t;
		}
		else
		{
			last->next = new_t;
		}
                
		last = new_t;
		rdfp = rdfp->next;
	}

	return head;
}

ReadDBFILEPtr LIBCALL 
readdb_destruct (ReadDBFILEPtr rdfp)

{
    ReadDBFILEPtr next;
    
    if (!rdfp)
        return NULL;
    
    while (rdfp) {
        next = rdfp->next;
        rdfp = readdb_destruct_element(rdfp);
        rdfp = next;
    }
    
    return NULL;
}

/*
	Destroys a single element.
*/
ReadDBFILEPtr LIBCALL 
readdb_destruct_element (ReadDBFILEPtr rdfp)

{

    if (rdfp == NULL)
        return NULL;
    
    /* Deallocate if contents were allocated. */
    /*if (rdfp->contents_allocated) {*/
    if (rdfp->parameters & READDB_CONTENTS_ALLOCATED) {
        rdfp->filename = (CharPtr)MemFree(rdfp->filename);
        rdfp->aliasfilename = (CharPtr)MemFree(rdfp->aliasfilename);
        rdfp->title = (CharPtr)MemFree(rdfp->title);
        rdfp->date = (CharPtr)MemFree(rdfp->date);
        /* free array if they were allocated, ie no memmap */
        if (rdfp->header_index_start)
            rdfp->header_index_start = (Uint4Ptr)MemFree(rdfp->header_index_start);
        if (rdfp->sequence_index_start)
            rdfp->sequence_index_start = (Uint4Ptr)MemFree(rdfp->sequence_index_start);
        if (rdfp->ambchar_index_start)
            rdfp->ambchar_index_start  =(Uint4Ptr) MemFree(rdfp->ambchar_index_start);
        /* is it completely safe to have one rdfp->nisam_opt for all threads. */
        ISAMObjectFree(rdfp->nisam_opt); /* Terminating NISAM */
        ISAMObjectFree(rdfp->sisam_opt); /* Terminating NISAM */
        OIDListFree(rdfp->oidlist);
    }
    rdfp->indexfp = NlmCloseMFILE(rdfp->indexfp);
    if (rdfp->shared_info) {  
       if (--(rdfp->shared_info->nthreads) == 0) {
	  rdfp->shared_info->sequencefp = 
	     NlmCloseMFILE(rdfp->shared_info->sequencefp); 
	  rdfp->shared_info->headerfp = 
	     NlmCloseMFILE(rdfp->shared_info->headerfp);  
       } else if (rdfp->shared_info->nthreads == -1) {
	  rdfp->shared_info->nthreads = 0;
	  rdfp->shared_info = NULL;
       }
    } 
    rdfp->shared_info = NULL;
    rdfp->sequencefp = NlmCloseMFILE(rdfp->sequencefp);
    rdfp->headerfp = NlmCloseMFILE(rdfp->headerfp);
    
    if (rdfp->allocated_length > 0) {
        rdfp->buffer = (UcharPtr)MemFree(rdfp->buffer);
    }
    
    
    /* destruct common index only if it is permited to do it for this thread */
        
if (rdfp->cih && /*rdfp->handle_common_index*/
    (rdfp->parameters & READDB_HANDLE_COMMON_INDEX))
        CommonIndexDestruct(rdfp->cih);
    
    rdfp = (ReadDBFILEPtr) MemFree(rdfp);
    
    return NULL;
}

/*
	Goes through a chain of ReadDBfILEPtr's, looking for the one
	that contains the specified ordinal ID.
*/

static ReadDBFILEPtr
readdb_get_link(ReadDBFILEPtr rdfp, Int4 ordinal_id)

{
   ReadDBFILEPtr last;

   last = rdfp;

   while (rdfp) {
      if (rdfp->start <=ordinal_id && rdfp->stop >= ordinal_id)
	 break;
      rdfp = rdfp->next;
   }
   if (!(last->parameters & READDB_KEEP_HDR_AND_SEQ)) {
      while (rdfp != last) {
	 if (last->sequencefp != NULL || last->headerfp != NULL) {
	    if (last->shared_info) {
	       NlmMutexLockEx(&hdrseq_mutex);
	       if (--(last->shared_info->nthreads) == 0) {
		  last->shared_info->sequencefp = 
		     NlmCloseMFILE(last->shared_info->sequencefp);
		  last->shared_info->headerfp = 
		     NlmCloseMFILE(last->shared_info->headerfp);
	       } else if (last->shared_info->nthreads < 0)
		  last->shared_info->nthreads = 0;
	       NlmMutexUnlock(hdrseq_mutex);
	    } 
	    last->sequencefp = NlmCloseMFILE(last->sequencefp);
	    last->headerfp = NlmCloseMFILE(last->headerfp);
	 }
	 last = last->next;
      }
   }
   if (rdfp->sequencefp==NULL && rdfp->headerfp==NULL) {
      NlmMutexLockEx(&hdrseq_mutex);
      if (ReadDBOpenMHdrAndSeqFiles(rdfp) == FALSE) {
	 ErrPostEx(SEV_ERROR, 0, 0, 
		   "ReadDBOpenMHdrAndSeqFiles: failed to map files\n");
	 rdfp = NULL;
      }
      NlmMutexUnlock(hdrseq_mutex);
   }
   return rdfp;
}
/*
  Returnes Int4 sequence_number by gi using NISAM indexes:
  
  ReadDBFILEPtr rdfp: the main ReadDB reference,
  Int4 gi - input gi number to find
  Int4 sequence_number: which number is this sequence,
  Returned 0 indicates, that gi was found
  Returned -1 indicates, that gi was not found
  Returned negative value mean fault of NISAM library
*/

Int4 LIBCALL
readdb_gi2seq(ReadDBFILEPtr rdfp, Int4 gi, Int4Ptr start)
{

    Boolean	thereis_unknown_database = FALSE;
    ReadDBFILEPtr	rdfp_start = rdfp;
    Boolean is_prot = (Boolean) (rdfp->parameters & READDB_IS_PROT);

    if (start)
	*start = 0;

    while(rdfp) {
	if (!rdfp->filebit) {
	    thereis_unknown_database = TRUE;
	    break;
	}
	rdfp = rdfp->next;
    }
	    
    rdfp = rdfp_start;

    if (thereis_unknown_database || (!isCommonIndex)) {
	ISAMErrorCode error;
	Uint4 Value;

	while (rdfp)
	{
	    if(rdfp->nisam_opt == NULL)
		return -1;

	    if((error = NISAMSearch(rdfp->nisam_opt, gi, 
		    &Value, NULL)) < 0) {
		ErrPostEx(SEV_WARNING, 0, 0, "Failed to initialize search. "
			"ISAM Error code is %d\n", error);
		return error;
	    }

	    if(error != ISAMNotFound) {
#if 0		
	       printf("\n%d\t-> %d\t in %s (%d)", gi, Value+rdfp->start,
		      rdfp->filename, is_prot);
#endif		
		if (start)
		    *start = rdfp->start;
		return (Int4) (Value+rdfp->start);
	    }

	    rdfp = rdfp->next;
	}
	return -1;  
    } else {	    
	Int4		retval = 0;
	Int4		mask = 0, alias_mask = 0;
	CommonIndexHeadPtr	cih = rdfp->cih;
	Int2		dbid=0, alias_dbid=0;

	/* create common mask for all databases */
	while (rdfp) {
	    if (rdfp->aliasfilebit) {
		alias_mask |= (0x1 << rdfp->aliasfilebit);
	    };
	    mask |= (0x1 << rdfp->filebit);
	    rdfp = rdfp->next;
	}

	/* get OID and database id (dbid) of this OID */
	if (cih)
	    retval = GI2OID(cih, gi, mask, alias_mask, &dbid, &alias_dbid, rdfp_start);

	if (retval >= 0) {
	    /* find correct rdfp in the list */
	    rdfp = rdfp_start;
	    while (rdfp) {
		/* if the oid found in mask database */
		if (alias_mask && rdfp->aliasfilebit == alias_dbid)
		    break;
		/* if the oid found in real database */
		if (!alias_mask && (rdfp->filebit == dbid))
		    break;
		rdfp = rdfp->next;
	    }

	    if (!rdfp) {
		/* we did not find the gi where we were trying */
		return -1;
	    }

	    if (start)
		*start = rdfp->start;

	    return retval+rdfp->start;
	}
	else
	    return -1;  
    }
}

Boolean readdb_find_best_id(SeqIdPtr sip, Int4Ptr gi, CharPtr tmpbuf)
{
    TextSeqIdPtr tsip = NULL;
    ObjectIdPtr oid;
    PDBSeqIdPtr psip;
    DbtagPtr dbt;
    SeqIdPtr sip_tmp;
    
    if (sip == NULL)
        return FALSE;

    for(sip_tmp = sip; sip_tmp != NULL; sip_tmp = sip_tmp->next) {
        if(sip_tmp->choice == SEQID_GI) {
            *gi = sip_tmp->data.intvalue;
            break;
        }
    }
    
    if(*gi != 0) return TRUE;
    
    for(sip_tmp = sip; sip_tmp != NULL; sip_tmp = sip_tmp->next) {
        
        switch (sip_tmp->choice) {       
        case SEQID_LOCAL:     /* local */
            oid = (ObjectIdPtr)(sip_tmp->data.ptrvalue);
            StringCpy(tmpbuf, oid->str);
            break;
        case SEQID_GIBBSQ:    /* gibbseq */
            sprintf(tmpbuf, "%ld", (long)(sip_tmp->data.intvalue));
            break;
        case SEQID_EMBL:      /* embl */
        case SEQID_DDBJ:      /* ddbj */
        case SEQID_GENBANK:   /* genbank */
        case SEQID_OTHER:     /* other */
        case SEQID_PIR:       /* pir   */
        case SEQID_SWISSPROT: /* swissprot */
        case SEQID_PRF:       /* prf   */
            tsip = (TextSeqIdPtr)(sip_tmp->data.ptrvalue);
            break;
        case SEQID_GENERAL:   /* general */
            dbt = (DbtagPtr)(sip_tmp->data.ptrvalue);
            StringCpy(tmpbuf, dbt->tag->str);
            break;
        case SEQID_PDB:       /* pdb   */
            psip = (PDBSeqIdPtr)(sip_tmp->data.ptrvalue);
            StringCpy(tmpbuf, psip->mol);
            break;
        }
    }
    
    if(tsip != NULL) {
        if(tsip->accession != NULL)
            StringCpy(tmpbuf, tsip->accession);
        else
            StringCpy(tmpbuf, tsip->name);
    }

    return TRUE;
}

/*
  Returnes Int4 sequence_number by SeqIdPtr using SISAM indexes:
  
  ReadDBFILEPtr rdfp: the main ReadDB reference,
  SeqIdPtr sip - input SeqIdPtr to find
  Int4 sequence_number: which number is this sequence,
  Returned 0 indicates, that gi was found
  Returned -1 indicates, that gi was not found
  Returned negative value mean fault of NISAM library
*/
Int4 LIBCALL
readdb_seqid2fasta(ReadDBFILEPtr rdfp, SeqIdPtr sip)
{
    ISAMErrorCode error;
    Int4 Value;
    Char tmpbuff[81];
    CharPtr key_out = NULL, data = NULL;
    Uint4 index;
    Int4 gi = 0;
    CharPtr chptr = NULL;
    TextSeqIdPtr tsip = NULL;

    if(rdfp->sisam_opt == NULL || sip == NULL)
        return -1;
    
    if(rdfp->sparse_idx) {
        readdb_find_best_id(sip, &gi, tmpbuff);
        if(gi != 0) {
            return readdb_gi2seq(rdfp, gi, NULL);
        }
    } else {

        switch (sip->choice) {
        case SEQID_EMBL:      /* embl */
        case SEQID_DDBJ:      /* ddbj */
        case SEQID_GENBANK:   /* genbank */
        case SEQID_OTHER:     /* other */
        case SEQID_PIR:       /* pir   */
        case SEQID_SWISSPROT: /* swissprot */
        case SEQID_PRF:       /* prf   */
            tsip = (TextSeqIdPtr)(sip->data.ptrvalue);
            break;
        default:
            break;
        }

        /* We have to clear name if both accession and name exists */
        if(tsip != NULL && tsip->accession != NULL && tsip->name != NULL) {
            chptr = tsip->name;
            tsip->name = NULL;
        }
        
        if((SeqIdWrite(sip, tmpbuff, 
                       PRINTID_FASTA_SHORT, sizeof(tmpbuff))) == NULL) {
            return -1;
        }
        
        /* Returning back name */
        if(chptr != NULL) {
            tsip->name = chptr;
        }
    }
    
    while (rdfp) {
	NlmMutexLockEx(&isamsearch_mutex);
        if((error = SISAMSearch(rdfp->sisam_opt, tmpbuff, 0, &key_out,
                                &data, &index)) < 0) {
            ErrPostEx(SEV_WARNING, 0, 0, "Failed to search string index "
                      "ISAM Error code is %d\n", error);
            return error;
        }
	NlmMutexUnlock(isamsearch_mutex);
        
        MemFree(key_out); /* We need no this for now */
        
        if(data && error != ISAMNotFound) {
            Value = atol(data);
            MemFree(data);
            return Value + rdfp->start;
        }
        rdfp = rdfp->next;
    }
    return -1;
}

/*
  Returnes array of sequence numbers by accession using SISAM indexes:

  ReadDBFILEPtr rdfp: the main ReadDB reference,
  CharPtr string - input accession to find
  Int4Ptr PNTR ids - array of sequence numbers
  Int4Ptr count - number of hits
  Returned  non-negative value indicates, that hits were found
  Returned -1 indicates, that hit(s) were not found
  Returned negative value mean fault of ISAM library
*/

Int4 LIBCALL
readdb_acc2fastaEx(ReadDBFILEPtr rdfp, CharPtr string, Int4Ptr PNTR ids,
                   Int4Ptr count)
{
    ISAMErrorCode error;
    SeqIdPtr sip;
    
    if(rdfp->sisam_opt == NULL || string == NULL)
        return -1;
    
    if (StringChr(string, '|') != NULL) {
        
        sip = SeqIdParse(string);        
        **ids = readdb_seqid2fasta(rdfp, sip);
        SeqIdFree(sip);
        
        if(**ids >= 0) {
            *count = 1;
            return 1;
        } else {
            return -1;
        }
    }
                                  
    while (rdfp) {
    	error =  SISAMFindAllData(rdfp->sisam_opt, string, ids, count);
        
    	if(error != ISAMNotFound) {
            return 1;
    	} 
	rdfp = rdfp->next;
    }
    return -1;
}
/*
  Returnes Int4 sequence_number by accession/locus using SISAM indexes:
  
  ReadDBFILEPtr rdfp: the main ReadDB reference,
  CharPtr string - input accession to find
  Int4 sequence_number: which number is this sequence,
  Returned 0 indicates, that gi was found
  Returned -1 indicates, that gi was not found
  Returned negative value mean fault of ISAM library
*/

Int4 LIBCALL
readdb_acc2fasta(ReadDBFILEPtr rdfp, CharPtr string)
{
    ISAMErrorCode error;
    Int4 Value;
    CharPtr key_out = NULL, data = NULL;
    Uint4 index;
    Char tmp_str[64];
    SeqIdPtr sip;

    if(rdfp->sisam_opt == NULL || string == NULL)
        return -1;
    
    if(rdfp->sparse_idx) {
        
        Int4 seq_num, count;
        Int4Ptr ids;
        
        readdb_acc2fastaEx(rdfp, string, &ids, &count);
        if(count > 0) {
            seq_num = *ids;
            MemFree(ids);
            return seq_num;
        } else {
            return -1;
        }
    }
    
    if (StringChr(string, '|') != NULL) {
        
        sip = SeqIdParse(string);        
        Value = readdb_seqid2fasta(rdfp, sip);
        SeqIdFree(sip);
        return Value;
    }
    
    while (rdfp) {
    	/* Trying accession first */
        
        sprintf(tmp_str, "gb|%s|", string);
        
        if((error = SISAMSearch(rdfp->sisam_opt, tmp_str, 0, &key_out,
                                &data, &index)) < 0) {
            ErrPostEx(SEV_WARNING, 0, 0, "Failed to search string index "
                      "ISAM Error code is %d\n", error);
            return error;
    	}
        
        MemFree(key_out); /* We need no this for now */
        
        if(error != ISAMNotFound) {
            Value = atol(data);
            MemFree(data);
            return Value + rdfp->start;
    	}
        
        /* Now trying LOCUS */
        
    	sprintf(tmp_str, "gb||%s", string);
        
    	if((error = SISAMSearch(rdfp->sisam_opt, tmp_str, 0, &key_out,
                                &data, &index)) < 0) {
            ErrPostEx(SEV_WARNING, 0, 0, "Failed to search string index "
                      "ISAM Error code is %d\n", error);
            return error;
    	}
        
    	MemFree(key_out); /* We need no this for now */
        
        if(error != ISAMNotFound) {
            Value = atol(data);
            MemFree(data);
            return Value + rdfp->start;
    	}
        
        /* Now trying string */
        
        
    	if((error = SISAMSearch(rdfp->sisam_opt, string, 0, &key_out,
                                &data, &index)) < 0) {
            ErrPostEx(SEV_WARNING, 0, 0, "Failed to search string index "
                      "ISAM Error code is %d\n", error);
            return error;
    	}
        
    	MemFree(key_out); /* We need no this for now */
        
    	if(error != ISAMNotFound) {
            Value = atol(data);
            MemFree(data);
            return Value + rdfp->start;
    	} else {
            MemFree(data);
    	}
	rdfp = rdfp->next;
    }
    
    return -1;
}

/*
	Obtains a BioseqPtr from readdb:

	ReadDBFILEPtr rdfp: the main ReadDB reference,
	Int4 sequence_number: which number is this sequence,
*/
BioseqPtr LIBCALL
readdb_get_bioseq(ReadDBFILEPtr rdfp, Int4 sequence_number)
{
    return readdb_get_bioseq_ex(rdfp, sequence_number, TRUE, FALSE);
}

BioseqPtr LIBCALL
readdb_get_bioseq_ex(ReadDBFILEPtr rdfp, Int4 sequence_number, 
                     Boolean use_objmgr, Boolean insert_ctrlA)

{
    BioseqPtr bsp;
    Boolean not_done;
    ByteStorePtr byte_store;
    CharPtr defline, new_defline = NULL, defline_ptr, new_defline_ptr, tmp_defline;
    Char seqid_buffer[50];
    CommonIndexPtr      cigi;
    Int2 byte_value, count;
    Int4 length, compressed_length, gi, new_defline_length, defline_length;
    SeqIdPtr sip, seqid, bestid;
    Uint1Ptr buffer, buffer_4na, buffer_2na; 
    Uint4Ptr ambchar = NULL;
    Uint4 header_index;
    Boolean is_prot = (Boolean) (rdfp->parameters & READDB_IS_PROT);
    
    if ((rdfp = readdb_get_link(rdfp, sequence_number)) == NULL)
       return NULL;

    defline = NULL;

    	readdb_get_descriptor(rdfp, sequence_number, &sip, &defline);
    
	if (insert_ctrlA == FALSE)
	{
    		count = 0;
    		new_defline = NULL;
    		if (defline != NULL) {
       		 defline_ptr = defline;
        
       		 while (*defline_ptr != NULLB) {
       		     count++;
       		     if (*defline_ptr == READDB_DEF_SEPARATOR) {
       	         /* Two spaces for every ctrl-A as it will be replaced by 2. */
       		         count++;
       		     }
       		     defline_ptr++;
       		 }
        
       	    	if (count != 0) {
       	     		new_defline = (CharPtr)Nlm_Malloc((count+1)*sizeof(Char));
       	     		new_defline_ptr = new_defline;
       	     		defline_ptr = defline;
       	     		while (*defline_ptr != NULLB) {
                		if (*defline_ptr == READDB_DEF_SEPARATOR) { 	
                    			*new_defline_ptr = ' ';
                    			new_defline_ptr++;
                    			*new_defline_ptr = '>';
                    			new_defline_ptr++;
               			} else {
                    			*new_defline_ptr = *defline_ptr;
                    			new_defline_ptr++;
                		}
                	defline_ptr++;
       	     		}
       	     		*new_defline_ptr = NULLB;
       	     		defline = (CharPtr)MemFree(defline);
       	     	}
    	    }
    	}
    	else
		new_defline = defline;
    
    if((length = readdb_get_sequence(rdfp, sequence_number, &buffer)) < 1)
        return NULL;
    
    if(use_objmgr) {
        if((bsp = BioseqNew()) == NULL)
            return NULL;
    } else {
        bsp = (BioseqPtr)MemNew(sizeof(Bioseq));
        if (bsp == NULL) return bsp;
        bsp->length = -1;    /* not set */
        bsp->topology = 1;   /* DEFAULT = linear */
    }

    byte_store = BSNew(0);
    if (is_prot) {
        bsp->mol = Seq_mol_aa;
        bsp->seq_data_type = Seq_code_ncbistdaa;
        BSWrite(byte_store, (VoidPtr) buffer, length);
    } else {
        /* Nucleotide sequence require more attention */
/*
        BSWrite(byte_store, (VoidPtr) buffer, length/4);
*/
        
/*
        if (length%4 != 0) {
            byte_value = *(buffer+length/4);
            byte_value &= 252; 
	    *(buffer_2na+length/4) = byte_value;
            BSPutByte(byte_store, byte_value);
        }
*/
        
        if(!readdb_get_ambchar(rdfp, sequence_number, &ambchar)) {
            ErrPostEx(SEV_WARNING, 0, 0, 
                      "Failure to read ambiguity information");
            return NULL;
        }
	/* Convert sequence if ambiguities. */
        if(ambchar != NULL) {/* are there any ambiguity ? */
		compressed_length = (length+3)/4; /* enough bytes for all bases. */
		buffer_2na = Nlm_Malloc(compressed_length*sizeof(Uint1));
		MemCopy(buffer_2na, buffer, compressed_length);
		if (length%4 != 0)
		{
            		byte_value = *(buffer_2na+length/4);
            		byte_value &= 252; 
			*(buffer_2na+length/4) = byte_value;
		}
	    	buffer_4na = MemNew((2*compressed_length)*sizeof(Uint1));
	    	MapNa2ByteToNa4String(buffer_2na, (Uint2Ptr) buffer_4na, compressed_length);
            	BSWrite(byte_store, (VoidPtr) buffer_4na, compressed_length*2);
/*
            byte_store =  BSConvertSeq(byte_store, Seq_code_ncbi4na, Seq_code_ncbi2na, length);
*/
            	byte_store = BSRebuildDNA_4na(byte_store, ambchar);
		MemFree(buffer_2na);
		MemFree(buffer_4na);
            	MemFree(ambchar);
            	bsp->seq_data_type = Seq_code_ncbi4na;
        }
	else
	{
        	BSWrite(byte_store, (VoidPtr) buffer, length/4);
        	if (length%4 != 0) {
            		byte_value = *(buffer+length/4);
            		byte_value &= 252; 
            		BSPutByte(byte_store, byte_value);
        	}
            	bsp->seq_data_type = Seq_code_ncbi2na;
	}
        
        bsp->mol = Seq_mol_na;
    }
    
    bsp->seq_data = byte_store;
    
    bsp->length = length;
    bsp->id = sip;
    bsp->repr = Seq_repr_raw;
    if (new_defline != NULL)  {
        bsp->descr = ValNodeAddStr(NULL, Seq_descr_title, new_defline);
    }
    return bsp;
}

/*
	returns the 'filebits' associated with a certain ordinal number.
	This is done by going to the rdfp for that ordinal id and
	gathering the filebits.
*/
Boolean LIBCALL 
readdb_get_filebits (ReadDBFILEPtr rdfp, Int4 ordinal_id, Uint2Ptr filebit, Uint2Ptr aliasfilebit)

{
        rdfp = readdb_get_link(rdfp, ordinal_id);

	if (rdfp == NULL)
		return FALSE;

	if (filebit)
		*filebit = rdfp->filebit;

	if (aliasfilebit)
		*aliasfilebit = rdfp->aliasfilebit;

	return TRUE;
}

	
/* 
	Gets the sequence number "sequence_number".  If memory-mapped
	files are enabled, then *buffer points to the appropriate place
	in the memory-mapped file.  If memory-mapped files are not enabled,
	then sufficient space in *buffer is allocated (if this is not already
	the case) and this length is stored in *buffer_length. 

	The length of the sequence requested is the return value; for memory-
	mapped files this is different than *buffer_length, which is always
	zero.
*/

Int4 LIBCALL 
readdb_get_sequence (ReadDBFILEPtr rdfp, Int4 sequence_number, Uint1Ptr PNTR buffer)

{
	Int4 diff, length, nitems=0;
	Uint1 remainder;
	Boolean is_prot = (Boolean) (rdfp->parameters & READDB_IS_PROT);

        rdfp = readdb_get_link(rdfp, sequence_number);

	if (rdfp == NULL || rdfp->sequencefp == NULL)
		return 0;

	if (is_prot == FALSE)
	{
		nitems = Nlm_SwapUint4(rdfp->ambchar_index[sequence_number]) - 
		    Nlm_SwapUint4(rdfp->sequence_index[sequence_number]);
	}
	else
	{
		nitems = Nlm_SwapUint4(rdfp->sequence_index[sequence_number+1]) - 
		    Nlm_SwapUint4(rdfp->sequence_index[sequence_number]) - 1;
	}

	NlmSeekInMFILE(rdfp->sequencefp,
		Nlm_SwapUint4(rdfp->sequence_index[sequence_number]),
		SEEK_SET);

        length = sizeof(Uint1) * nitems;
	/* Use memory-mapped file, don't allocate buffer. */
	if (rdfp->sequencefp->mfile_true == TRUE)
	{
        	diff = rdfp->sequencefp->mmp_end - rdfp->sequencefp->mmp;

        	if (length > diff)
        	{
               		 nitems = diff / sizeof(Uint1);
               		 length = nitems * sizeof(Uint1);
        	}
        	*buffer = rdfp->sequencefp->mmp;
	}
	else
	{
	/* No mem-mapping, allocate a buffer for the subject sequence. */
		if (length > rdfp->allocated_length)
		{
			if (rdfp->buffer != NULL)
				rdfp->buffer = (UcharPtr)MemFree(rdfp->buffer);
			rdfp->allocated_length = length+2;
			rdfp->buffer = (UcharPtr)MemNew((rdfp->allocated_length)*sizeof(Uint1));
		}
/* For protein db's the first and last byte is the NULLB, which is a sentinel byte 
used by the extension functions. For nucl. db's there are no sentinel bytes. */
		if (is_prot)
		{
			rdfp->buffer[0] = NULLB;
			*buffer = rdfp->buffer+1;
			FileRead(*buffer, sizeof(Uint1), nitems+1, rdfp->sequencefp->fp);
		}
		else
		{
			*buffer = rdfp->buffer;
			FileRead(*buffer, sizeof(Uint1), nitems, rdfp->sequencefp->fp);
		}
	}

	/* For nucl. return "unpacked" length and get the remainder out
	of the last byte. */
	if (is_prot == FALSE)
	{
/* The first six bits in the byte holds the "remainder" (not a multiple of 4) 
and the last two bits of the byte holds the size of the remainder (0-3). */
		remainder = *(*buffer+length-1);
		remainder &= 0x3;
		length--;
/* 4 bases per byte. */
		length *= 4;
		length += remainder;
	}

	return length;
}
	
/* 
	Gets the sequence number "sequence_number".  The sequence returned includes
	all ambiguity information.  THis funciton should only be used for nucleic
	acid sequences, for proteins use readdb_get_sequence.

	buffer contains the sequence and is reallocated if *buffer_length is not long enough.

	The length of the sequence requested is the return value.
	protein sequences are always returned as Seq_code_ncbistdaa,
	nucleotide sequences as Seq_code_ncbi4na.
*/

Int4 LIBCALL 
readdb_get_sequence_ex (ReadDBFILEPtr rdfp, Int4 sequence_number, Uint1Ptr PNTR buffer, Int4 *buffer_length)

{
    	ByteStorePtr byte_store;
    	Int2 byte_value;
	Int4 index, index2, length, copy_length;
	Uint1Ptr private_buffer, buffer_4na;
    	Uint4Ptr ambchar = NULL;
	Boolean is_prot = (Boolean) (rdfp->parameters & READDB_IS_PROT);

	length = readdb_get_sequence(rdfp, sequence_number, &private_buffer);

	/* Check the length, make it one longer for ALIGN. */
	if ((length+1) > *buffer_length || *buffer == NULL)
	{
		if (*buffer)
			MemFree(*buffer);

		*buffer = MemNew((length+1)*sizeof(Uint1));
		*buffer_length = length+1;
	}

	/* Copy sequence into allocated buffer. */
	if (!is_prot)
	{
		copy_length = length/4;
		if (length%4 != 0) 
			copy_length++;
	}
	else
	{
		copy_length = length;
	}
	MemCpy((VoidPtr) *buffer, private_buffer, copy_length);

    	if (!is_prot) 
	{
    		byte_store = BSNew(0);
       	 	/* Nucleotide sequence require more attention */
        
		if (length%4 != 0)
		{
            		byte_value = *(*buffer+length/4);
            		byte_value &= 252; 
			*(*buffer+length/4) = byte_value;
		}
	    	buffer_4na = MemNew((2*copy_length)*sizeof(Uint1));
	    	MapNa2ByteToNa4String(*buffer, (Uint2Ptr) buffer_4na, copy_length);
            	BSWrite(byte_store, (VoidPtr) buffer_4na, copy_length*2);
        
        	if(!readdb_get_ambchar(rdfp, sequence_number, &ambchar)) {
            		ErrPostEx(SEV_WARNING, 0, 0, 
               	       "Failure to read ambiguity information");
            		return -1;
        	}
		/* Convert sequence if ambiguities. */
       		if(ambchar != NULL) {/* are there any ambiguity ? */
       		     byte_store = BSRebuildDNA_4na(byte_store, ambchar);
       		     MemFree(ambchar);
    		}
		MemFree(buffer_4na);

		/* Sequence is copied back to *buffer. */
		BSSeek(byte_store, 0, SEEK_SET);
		BSMerge(byte_store, *buffer);
		BSFree(byte_store);

		private_buffer = *buffer;
		index = length/2 - 1;
		index2 = length-1;
		if (length%2 != 0)
		{
			private_buffer[index2] = (private_buffer[index+1] >> 4);
			index2--;
		}
		while (index2 > 0)
		{
			private_buffer[index2] = (private_buffer[index] & 15);
			index2--; 
			private_buffer[index2] = (private_buffer[index] >> 4);
			index2--; index--;
		}
		
	}

	return length;
}

/* 
	Gets the length of sequence number "sequence_number". 
*/

Int4 LIBCALL 
readdb_get_sequence_length (ReadDBFILEPtr rdfp, Int4 sequence_number)

{

	Int4 length;
	Uint1 remainder;
	Boolean is_prot = (Boolean) (rdfp->parameters & READDB_IS_PROT);

        rdfp = readdb_get_link(rdfp, sequence_number);

	if (rdfp == NULL)
		return 0;

	if (is_prot == FALSE)
	{
		length = Nlm_SwapUint4(rdfp->ambchar_index[sequence_number]) -
		    Nlm_SwapUint4(rdfp->sequence_index[sequence_number]);
	}
	else
	{
		length = Nlm_SwapUint4(rdfp->sequence_index[sequence_number+1]) -
		    Nlm_SwapUint4(rdfp->sequence_index[sequence_number]) - 1;
	}

	/* For nucl. return "unpacked" length and get the remainder out
	of the last byte. */
	if (is_prot == FALSE)
	{
		if (rdfp->sequencefp->mfile_true == TRUE)
		{
			NlmSeekInMFILE(rdfp->sequencefp, 
				Nlm_SwapUint4(rdfp->ambchar_index[sequence_number])-1, SEEK_SET);
			remainder = *(rdfp->sequencefp->mmp);
		}
		else
		{
			NlmSeekInMFILE(rdfp->sequencefp, 
				Nlm_SwapUint4(rdfp->ambchar_index[sequence_number])-1, SEEK_SET);
			NlmReadMFILE((Uint1Ptr) &remainder, 1, 1, rdfp->sequencefp);
		}
/* The first six bits in the byte holds the "remainder" (not a multiple of 4) 
and the last two bits of the byte holds the size of the remainder (0-3). */
		remainder &= 3;
		length--;
/* 4 bases per byte. */
		length *= 4;
		length += remainder;
	}

	return length;
}
#ifdef FASTA_ASN
/*
Get the FasfaPtr (ASN.1) for the sequence with sequence_number.
It is the caller's RESPONSIBILITY to DEALLOCATE Fasta ASN.1".
*/
FdbFastaPtr LIBCALL readdb_get_fastaid PROTO((ReadDBFILEPtr rdfp,
                                           Int4 sequence_number))
{
  FdbFastaPtr fasta;
  AsnIoPtr aip;
  AsnIoMemPtr aimp;
  Int4 size;

  rdfp = readdb_get_link(rdfp, sequence_number);

  if (rdfp == NULL)
    return FALSE;

  size = Nlm_SwapUint4(rdfp->header_index[sequence_number+1]) - 
      Nlm_SwapUint4(rdfp->header_index[sequence_number]);
  
  if (rdfp->headerfp->mfile_true == TRUE) {
    NlmSeekInMFILE(rdfp->headerfp,
                   Nlm_SwapUint4(rdfp->header_index[sequence_number]),
		   SEEK_SET);
    aimp = AsnIoMemOpen("rb", rdfp->headerfp->mmp, size);    
    fasta = FdbFastaAsnRead(aimp->aip, NULL);
    AsnIoMemClose(aimp);
  } else {
    aip = AsnIoNew(ASNIO_BIN_IN, rdfp->headerfp->fp, NULL, NULL, NULL);  
    NlmSeekInMFILE(rdfp->headerfp,
                   Nlm_SwapUint4(rdfp->header_index[sequence_number]),
		   SEEK_SET);
    fasta = FdbFastaAsnRead(aimp->aip, NULL);   
    AsnIoFree(aip, FALSE);
  }
  return fasta;
}
#endif
Boolean  LIBCALL
readdb_get_ambchar (ReadDBFILEPtr rdfp, Int4 sequence_number, Uint4Ptr PNTR ambchar_return)
{
  Uint4Ptr ambchar;
  Int4 length, index, total;

  rdfp = readdb_get_link(rdfp, sequence_number);

  if((length = Nlm_SwapUint4(rdfp->sequence_index[sequence_number+1]) -
	      Nlm_SwapUint4(rdfp->ambchar_index[sequence_number])) == 0) {
    *ambchar_return = NULL;
    return TRUE;    /* no ambiguous characters available */
  }
    
    /* Each ambig. residue is represented by a Uint4, 
       but length is in bytes. */

    total = length/4;
    if((ambchar = (Uint4Ptr)MemNew(total*sizeof(Uint4))) == NULL)
      return FALSE;

    NlmSeekInMFILE(rdfp->sequencefp, 
                   Nlm_SwapUint4(rdfp->ambchar_index[sequence_number]), SEEK_SET);
    
    NlmReadMFILE((Uint1Ptr) ambchar, 4, total, rdfp->sequencefp);
    for (index=0; index<total; index++) {
      ambchar[index] = Nlm_SwapUint4(ambchar[index]);
    }
    
  *ambchar_return = ambchar;
  return TRUE;  
}

/*
	Check if ambiguity characters are present in the sequence. 
*/

Boolean LIBCALL
readdb_ambchar_present (ReadDBFILEPtr rdfp, Int4 sequence_number)

{
  	rdfp = readdb_get_link(rdfp, sequence_number);
	if (rdfp == NULL)
		return FALSE;

	if (rdfp->ambchar_index == NULL)
		return FALSE;

	if((Nlm_SwapUint4(rdfp->sequence_index[sequence_number+1]) -
		    Nlm_SwapUint4(rdfp->ambchar_index[sequence_number])) == 0)
	{
		return FALSE;
	}

	return TRUE;
}

static Boolean
readdb_adjust_local_id(ReadDBFILEPtr rdfp, SeqIdPtr sip)

{
	DbtagPtr dbtag;
	ObjectIdPtr oid;

	if (sip == NULL || sip->choice != SEQID_GENERAL)
		return FALSE;

	if (rdfp->start == 0)
		return TRUE;

	dbtag = sip->data.ptrvalue;
	oid = dbtag->tag;
	oid->id += rdfp->start;

	return TRUE;

	

}

Boolean LIBCALL
readdb_get_descriptor (ReadDBFILEPtr rdfp, Int4 sequence_number, SeqIdPtr PNTR id, CharPtr PNTR description)

{
	Boolean not_done;
        Char buffer[READDB_BUF_SIZE], id_buf[READDB_BUF_SIZE];
        CharPtr buf_ptr, defline, new_defline=NULL, tmp_defline;
	CommonIndexPtr      cigi;
	Int4 index, new_size, alias_mask=0, gi;
	Int4 defline_length, new_defline_length;
	SeqIdPtr bestid, seqid;
	Uint2 aliasfilebit=0;
	Uint4 header_index;

  	rdfp = readdb_get_link(rdfp, sequence_number);
	if (rdfp == NULL)
		return FALSE;

    if (rdfp->oidlist)
    {
	readdb_get_filebits(rdfp, sequence_number, NULL, &aliasfilebit);	
    }

    if (aliasfilebit != 0)
    {
	alias_mask |= (0x1 << aliasfilebit);

	*id = NULL;
	not_done = TRUE;
	header_index = 0;
	while (not_done)
	{
		not_done = readdb_get_header(rdfp, sequence_number, &header_index, &seqid, &defline);
		if (not_done == FALSE)
			break;
		bestid = SeqIdFindBest(seqid, SEQID_GI);
		gi = bestid->data.intvalue;
		cigi = rdfp->cih->ci + gi;
		if (alias_mask & SwapUint4(cigi->dbmask))
		{
			if (*id == NULL)
			{
				*id = seqid;
				seqid = NULL;
				new_defline = defline;
				new_defline_length = StringLen(new_defline);
				defline = NULL;
			}
			else
			{
				SeqIdWrite(seqid, id_buf, PRINTID_FASTA_LONG, READDB_BUF_SIZE);
				seqid = SeqIdSetFree(seqid);
				defline_length = new_defline_length;
				new_defline_length += StringLen(defline) + StringLen(id_buf);
				new_defline_length += 3;
				tmp_defline = MemNew(new_defline_length+1);
				MemCpy(tmp_defline, new_defline, defline_length);
				sprintf(tmp_defline+defline_length, "%c%s %s", READDB_DEF_SEPARATOR, id_buf, defline);
				defline = MemFree(defline);
				new_defline = MemFree(new_defline);
				new_defline = tmp_defline;
			}
		}
		else
		{
			seqid = SeqIdSetFree(seqid);
			defline = MemFree(defline);
		}
	}
	if (description != NULL)
		*description = new_defline;
	else
		new_defline = MemFree(new_defline);
    }
    else
    {
	new_size = Nlm_SwapUint4(rdfp->header_index[sequence_number+1]) - 
	    Nlm_SwapUint4(rdfp->header_index[sequence_number]);

	if (new_size > READDB_BUF_SIZE)
	{
		buf_ptr = (CharPtr)Nlm_Malloc(new_size*sizeof(Char) + 1);
	}
	else
	{
		buf_ptr = &buffer[0];
	}

	NlmSeekInMFILE(rdfp->headerfp, Nlm_SwapUint4(rdfp->header_index[sequence_number]),
		SEEK_SET);
	if (NlmReadMFILE((Uint1Ptr) buf_ptr, sizeof(Char), new_size, rdfp->headerfp) != new_size)
		return FALSE;
	buf_ptr[new_size] = NULLB;	/* defline saved w/o NULLB. */

	for (index=0; index<READDB_BUF_SIZE; index++)
	{
		if (buf_ptr[index] == ' ' || buf_ptr[index] == NULLB)
		{
			id_buf[index] = NULLB;
			index++;
			break;
		}
		id_buf[index] = buf_ptr[index];
	}

	*id = SeqIdParse(id_buf);
	readdb_adjust_local_id(rdfp, *id);
	if (description != NULL)
		*description = StringSave(&buf_ptr[index]);

	if (buf_ptr != &buffer[0])
		buf_ptr = (CharPtr)MemFree(buf_ptr);
   }
	
	return TRUE;
}
Boolean
readdb_get_defline (ReadDBFILEPtr rdfp, Int4 sequence_number, CharPtr PNTR description)

{
    Char buffer[READDB_BUF_SIZE];
    CharPtr buf_ptr;
    Int4 new_size;
    
    rdfp = readdb_get_link(rdfp, sequence_number);
    if (rdfp == NULL)
        return FALSE;
    
    new_size = Nlm_SwapUint4(rdfp->header_index[sequence_number+1]) -
	Nlm_SwapUint4(rdfp->header_index[sequence_number]);
    
    if (new_size > READDB_BUF_SIZE){
        buf_ptr = (CharPtr)Nlm_Malloc(new_size*sizeof(Char) + 1);
    } else {
        buf_ptr = &buffer[0];
    }
    
    NlmSeekInMFILE(rdfp->headerfp, Nlm_SwapUint4(rdfp->header_index[sequence_number]), 
                   SEEK_SET);
    if (NlmReadMFILE((Uint1Ptr) buf_ptr, sizeof(Char), new_size, 
                     rdfp->headerfp) != new_size)
        return FALSE;
    
    buf_ptr[new_size] = NULLB;	/* defline saved w/o NULLB. */
    
    *description = StringSave(buf_ptr);

    if (buf_ptr != &buffer[0])
        buf_ptr = (CharPtr)MemFree(buf_ptr);
    
    return TRUE;
}


/*
	A single sequence may be attched to several entries (as they all
	have the same sequence).  This function gets the ID and deflines for 
	each entry attched to one sequence.  On the first call the Uint4
	(*header_index) should be zero; it will be filled in by readdb_get_header.  
	Subsequent calls will use this information to know which ID and
	defline to retrieve next.  When all are retrieved, FALSE will be returned.
*/

Boolean LIBCALL
readdb_get_header (ReadDBFILEPtr rdfp, Int4 sequence_number, Uint4Ptr header_index, SeqIdPtr PNTR id, CharPtr PNTR description)

{
        CharPtr id_buf;
        CharPtr buf_ptr, buf_defline_start;
	Int4 index, size;
	Uint4 header_index_end;

  	rdfp = readdb_get_link(rdfp, sequence_number);
	if (rdfp == NULL)
		return FALSE;

	if (*header_index == 0)
		*header_index = Nlm_SwapUint4(rdfp->header_index[sequence_number]);

	header_index_end = Nlm_SwapUint4(rdfp->header_index[sequence_number+1]);

	if (*header_index >= header_index_end)
	{
		*header_index = 0;
		return FALSE;
	}

	size = header_index_end-(*header_index);
	buf_ptr = MemNew((size+1)*sizeof(Char));
	id_buf = MemNew((size+1)*sizeof(Char));

	NlmSeekInMFILE(rdfp->headerfp, (long) *header_index, SEEK_SET);
	if (NlmReadMFILE((Uint1Ptr) buf_ptr, sizeof(Char), size, rdfp->headerfp) != size)
		return FALSE;

	for (index=0; index<size; index++)
	{
		if (buf_ptr[index] == ' ')
		{
			id_buf[index] = NULLB;
			index++;
			break;
		}
		id_buf[index] = buf_ptr[index];
	}

	*id = SeqIdParse(id_buf);

	buf_defline_start = &buf_ptr[index];
	while (index < size)
	{
		if (buf_ptr[index] == READDB_DEF_SEPARATOR)
		{
			break;
		}	
		index++;
	}
	buf_ptr[index] = NULLB;
	index++;
	if (description != NULL)
	{
		*description = StringSave(buf_defline_start);
	}

	buf_ptr = MemFree(buf_ptr);
	id_buf = MemFree(id_buf);
	*header_index += index;
	
	return TRUE;
}

/*
	Obtains the total database length from the ReadDBFILE structure. 
*/
Int8 LIBCALL
readdb_get_dblen (ReadDBFILEPtr rdfp)

{
	if (rdfp == NULL)
		return 0;

	return rdfp->totlen;
}

/*
Obtains the total number of database sequences from all the ReadDBFILE structures. 
*/
Int4 LIBCALL
readdb_get_num_entries_total (ReadDBFILEPtr rdfp)

{
	Int4 total=0;
	if (rdfp == NULL)
		return 0;

	while (rdfp)
	{
		total += rdfp->num_seqs;
		rdfp = rdfp->next;
	}
	return total;
}

/*
Obtains the total number of real database sequences from all the ReadDBFILE structures. 
*/
Int4 LIBCALL
readdb_get_num_entries_total_real (ReadDBFILEPtr rdfp)

{
	Int4 total=0;
	if (rdfp == NULL)
		return 0;

	/* if first database contains no sequences, that means, that
	   the only mask should be used, so size of real databases is zero */
	if (rdfp->num_seqs == 0)
	    return 0;

	while (rdfp && !rdfp->oidlist)
	{
		total += rdfp->num_seqs;
		rdfp = rdfp->next;
	}
	return total;
}

/*
Obtains the number of database sequences from the ReadDBFILE structure. 
*/
Int4 LIBCALL
readdb_get_num_entries (ReadDBFILEPtr rdfp)

{
	if (rdfp == NULL)
		return 0;

	return rdfp->num_seqs;
}

/*
Obtains the length of the longest database seq from the ReadDBFILE structure. 
*/
Int4 LIBCALL
readdb_get_maxlen (ReadDBFILEPtr rdfp)

{
	if (rdfp == NULL)
		return 0;

	return rdfp->maxlen;
}

/*
Obtains the title of the database.  Note that the return CharPtr is not
owned by the caller.  It should be copied if the user wishes to modify it.
*/
CharPtr LIBCALL
readdb_get_filename (ReadDBFILEPtr rdfp)

{
	if (rdfp == NULL)
		return NULL;

        if (rdfp->aliasfilename)
                return rdfp->aliasfilename;

	return rdfp->filename;
}

/*
Obtains the title of the database.  Note that the return CharPtr is not
owned by the caller.  It should be copied if the user wishes to modify it.
*/
CharPtr LIBCALL
readdb_get_title (ReadDBFILEPtr rdfp)

{
	if (rdfp == NULL)
		return NULL;

	return rdfp->title;
}

/*
Obtains the date and time the database was formatted with formatdb.
Note that the return CharPtr is not owned by the caller.  It should 
be copied if the user wishes to modify it.
*/
CharPtr LIBCALL
readdb_get_date (ReadDBFILEPtr rdfp)

{
	if (rdfp == NULL)
		return NULL;

	return rdfp->date;
}

/*
Queries readdb whether the sequence is protein.
*/
Boolean LIBCALL
readdb_is_prot (ReadDBFILEPtr rdfp)

{
	if (rdfp == NULL)
		return FALSE;

	return /*rdfp->is_prot*/(Boolean) (rdfp->parameters & READDB_IS_PROT);
}

/*
Obtains the formatdb version used to format the database.
*/
Int4 LIBCALL
readdb_get_formatdb_version (ReadDBFILEPtr rdfp)

{
	if (rdfp == NULL)
		return 0;

	return rdfp->formatdb_ver;
}

/* 
	Translates a SeqIdPtr to an ordinal ID, used by the BLAST database.
	If the SeqIdPtr cannot be translated, a negative number is returned. 
	All valid ordinal numbers are >= 0.
*/

Int4 SeqId2OrdinalId(ReadDBFILEPtr rdfp, SeqIdPtr sip)

{
	DbtagPtr	dbtagptr;
	Int4 ordinal_id;

	if (rdfp == NULL || sip == NULL)
		return -2;

	switch (sip->choice)
	{
		case SEQID_GI:
			ordinal_id = readdb_gi2seq(rdfp, sip->data.intvalue, NULL);
			break;

		case SEQID_GENERAL:
			dbtagptr = (DbtagPtr) sip->data.ptrvalue;
			if (dbtagptr == NULL)
				return OM_MSG_RET_OK;
			if (StringCmp(dbtagptr->db, "BL_ORD_ID") == 0)
			{
				ordinal_id = dbtagptr->tag->id;
				break;
			}
			/* Fall through to default if not "BL_ORD_ID" */
		default:
			ordinal_id = readdb_seqid2fasta(rdfp, sip);
			break;
	}

	return ordinal_id;
}

/*************************************************************************

	Inits the ReadDBFILEPtr for the BioseqFetch functions.

**************************************************************************/

static Boolean
ReadDBInit(ReadDBFetchStructPtr rdfsp)
{

	rdfsp->rdfp = readdb_new_ex(rdfsp->dbname, rdfsp->is_prot, TRUE);

	if (rdfsp->rdfp != NULL)
		return TRUE;
	else
		return FALSE;
}

/*
	Checks the chain of ReadDBFetchStructPtr's for one
	which belongs to the calling thread. If none is found,
	NULL isreturned; otherwise the ReadDBFetchStructPtr is
	returned.
*/
static ReadDBFetchStructPtr
ReadDBFindFetchStruct(ReadDBFetchStructPtr rdfp)

{

	if (rdfp == NULL)
		return NULL;

	while (rdfp)
	{
		if (NlmThreadCompare(rdfp->thread_id, NlmThreadSelf()) == TRUE)
			break;
		rdfp = rdfp->next;
	}
	return rdfp;
}

/*
	Initializes the ReadDBFetchStructPtr and adds onto end of
	chain of ReadDBFetchStructPtr (head).  The new ReadDBFetchStructPtr
	is returned.
*/
static ReadDBFetchStructPtr
ReadDBFetchStructNew(ReadDBFetchStructPtr head, CharPtr dbname, Boolean is_na)

{
	ReadDBFetchStructPtr rdfsp, rdfsp_var;

	
	rdfsp = (ReadDBFetchStructPtr) MemNew(sizeof(ReadDBFetchStruct));
	rdfsp->dbname = StringSave(dbname);
	rdfsp->is_prot = (is_na == TRUE) ? FALSE : TRUE;
	rdfsp->thread_id = NlmThreadSelf();

	if (head != NULL)
	{
		rdfsp_var = head;
		while (rdfsp_var->next)
			rdfsp_var = rdfsp_var->next;
		rdfsp_var->next = rdfsp;
	}
	
	return rdfsp;
}

/****************************************************************
*
*	ReadDBFetchFreeFunc
*	Frees ReadDBFetchUserData.
*
****************************************************************/
	
static Pointer LIBCALLBACK ReadDBFetchFreeFunc (Pointer ptr)
{
	ReadDBFetchUserDataPtr userdata;

	userdata = (ReadDBFetchUserDataPtr) ptr;
	return MemFree(userdata);
}



/**********************************************************************

	Fetches the Bioseq, based on the ordinal number of the
	sequence in the database.

************************************************************************/

static Int2 LIBCALLBACK ReadDBBioseqFetchFunc(Pointer data)
{
	BioseqPtr bsp, core_bsp;
	Boolean status;
	Int4 ordinal_id;
	OMProcControlPtr ompcp;
        ObjMgrProcPtr ompp;
	OMUserDataPtr omdp;
	ReadDBFetchStructPtr rdfsp;
	ReadDBFILEPtr rdfp=NULL;
	ReadDBFetchUserDataPtr userdata;
	SeqIdPtr sip, best_id;
	SeqEntryPtr sep;

	ordinal_id = -1;
	
	ompcp = (OMProcControlPtr)data;
        ompp = ompcp->proc;

	rdfsp = ReadDBFindFetchStruct((ReadDBFetchStructPtr)(ompp->procdata));

	if (rdfsp->ReadDBFetchState == READDBBF_DISABLE)
	{
		return OM_MSG_RET_OK;
	}

	if (rdfsp->ReadDBFetchState == READDBBF_INIT)
	{
		status = ReadDBInit(rdfsp);
		if (status == FALSE)
			return OM_MSG_RET_OK;
		rdfsp->ReadDBFetchState = READDBBF_READY;
	}

#if 0
        if (ompcp->input_entityID)
        {
                omdp = ObjMgrGetUserData(ompcp->input_entityID, ompp->procid, OMPROC_FETCH, 0);
                if (omdp != NULL)
		{
			userdata = (ReadDBFetchUserDataPtr) (omdp->userdata.ptrvalue);	
			if (userdata != NULL)
			{
				ordinal_id = userdata->ordinal_number;
				rdfp = ReadDBGetDb(rdfsp->rdfp, userdata->db_id);
			}
		}
	}
#endif

	if (ordinal_id < 0 || rdfp == NULL)
	{
		sip = (SeqIdPtr) (ompcp->input_data);

		best_id = SeqIdFindBest(sip, SEQID_GI);

		if (best_id == NULL)
		{
			core_bsp = BioseqFindCore(sip);
			if (core_bsp)
				best_id = SeqIdFindBest(core_bsp->id, SEQID_GI);
		}

		if (best_id == NULL)
			return OM_MSG_RET_OK;

		rdfp = rdfsp->rdfp;
		while (rdfp)
		{	/* look in all databases for the proper ID. */
			ordinal_id = SeqId2OrdinalId(rdfp, best_id);
			if (ordinal_id >= 0)
				break;
			rdfp = rdfp->next;
		}

	}

	/* ordinal_id's start at zero. */
	if (ordinal_id < 0)
		return OM_MSG_RET_OK;

	/* A BioseqPtr is returned by this function. */
	bsp = readdb_get_bioseq(rdfp, ordinal_id);

        /* We have to add information about genetic code to
           the Bioseq */

        if(rdfsp->db_genetic_code > 1) {
            BioSourcePtr source;
            source = BioSourceNew();
            source->org = OrgRefNew();
            source->org->orgname = OrgNameNew();
            source->org->orgname->gcode = rdfsp->db_genetic_code;
            ValNodeAddPointer(&(bsp->descr), Seq_descr_source, source);
        }

	sep = SeqEntryNew();
	sep->choice = 1;
	sep->data.ptrvalue = bsp;
	SeqMgrSeqEntry(SM_BIOSEQ, (Pointer)bsp, sep);
	ompcp->output_data = (Pointer)bsp;
	ompcp->output_entityID = ObjMgrGetEntityIDForChoice(sep);
	omdp = ObjMgrAddUserData(ompcp->output_entityID, ompp->procid, OMPROC_FETCH, 0);
	userdata = (ReadDBFetchUserDataPtr) MemNew(sizeof(ReadDBFetchUserData));
	omdp->userdata.ptrvalue = userdata;
	userdata->ordinal_number = ordinal_id;
	userdata->db_id = ReadDBGetDbId(rdfsp->rdfp, rdfp);
	omdp->freefunc = ReadDBFetchFreeFunc;

	return OM_MSG_RET_DONE;
}

Boolean LIBCALL ReadDBBioseqSetDbGeneticCode(Int4 db_genetic_code)
{
    ReadDBFetchStructPtr rdfsp;
    ObjMgrPtr omp;
    ObjMgrProcPtr ompp;
    
    omp = ObjMgrGet();
    ompp = ObjMgrProcFind(omp, 0, "ReadDBBioseqFetch", OMPROC_FETCH);
    if (ompp != NULL) {   /* already initialized */
        rdfsp = ReadDBFindFetchStruct((ReadDBFetchStructPtr)(ompp->procdata));
        rdfsp->db_genetic_code = db_genetic_code;
        return FALSE;
    }
    return TRUE;
}

/*********************************************************************

	Enables the fetching.  Initializes needed structures and calls
	ReadDBInit.

**********************************************************************/
Boolean LIBCALL 
ReadDBBioseqFetchEnable(CharPtr program, CharPtr dbname, Boolean is_na, Boolean now)

{
    Boolean result;
    ReadDBFetchStructPtr rdfsp;
    ObjMgrPtr omp;
    ObjMgrProcPtr ompp;
    static TNlmMutex enable_lock = NULL;
    /* check if already enabled ***/

    NlmMutexInit(&enable_lock);
    NlmMutexLock(enable_lock);

    omp = ObjMgrGet();
    ompp = ObjMgrProcFind(omp, 0, "ReadDBBioseqFetch", OMPROC_FETCH);
    if (ompp != NULL) {   /* already initialized */
        rdfsp = ReadDBFindFetchStruct((ReadDBFetchStructPtr)(ompp->procdata));

        if(rdfsp == NULL) { /* Another thread */
            rdfsp = ReadDBFetchStructNew((ReadDBFetchStructPtr)(ompp->procdata), dbname, is_na);
        } else {
            if (rdfsp->is_prot == is_na || StringCmp(rdfsp->dbname, dbname)) {
                rdfsp->is_prot = (is_na == TRUE) ? FALSE : TRUE;
                rdfsp->dbname = MemFree(rdfsp->dbname);
                rdfsp->dbname = StringSave(dbname);
            }
        }
    } else { /* New element is not registered with ObjMgr */
        rdfsp = ReadDBFetchStructNew(NULL, dbname, is_na);
        ObjMgrProcLoad(OMPROC_FETCH, "ReadDBBioseqFetch", "ReadDBBioseqFetch", OBJ_SEQID, 0,OBJ_BIOSEQ,0,
                       (Pointer)rdfsp, ReadDBBioseqFetchFunc, PROC_PRIORITY_DEFAULT);
        rdfsp->ReadDBFetchState = READDBBF_INIT;
    }
    
    rdfsp->ctr++;    /* count number of enables */
    
    NlmMutexUnlock(enable_lock);
    
    if (rdfsp->ReadDBFetchState == READDBBF_READY) {
        return TRUE;
    }
    
    if (now) {
        result = ReadDBInit(rdfsp);
        if (! result) {
            return result;
        }
        rdfsp->ReadDBFetchState = READDBBF_READY;
    } else {
        rdfsp->ReadDBFetchState = READDBBF_INIT;
    }
    
    return TRUE;
}

/*****************************************************************************
*
*		ReadDBBioseqFetchDisable()
*
*	Calls readdb_destruct if necessary to deallocate resources.
*
*****************************************************************************/
void LIBCALL ReadDBBioseqFetchDisable(void)
{
        ObjMgrPtr omp;
        ObjMgrProcPtr ompp;
        ReadDBFetchStructPtr rdfsp;

        omp = ObjMgrGet();
        ompp = ObjMgrProcFind(omp, 0, "ReadDBBioseqFetch", OMPROC_FETCH);
        if (ompp == NULL)   /* not initialized */
                return;

	rdfsp = ReadDBFindFetchStruct((ReadDBFetchStructPtr)(ompp->procdata));
	if (! rdfsp->ctr)   /* no enables active */
		return;

	rdfsp->ctr--;
	if (rdfsp->ctr)   /* connection still pending */
			  return;

        if (rdfsp->ReadDBFetchState == READDBBF_READY)
	{
		rdfsp->ReadDBFetchState = READDBBF_DISABLE;  /* not active */
		rdfsp->rdfp = ReadDBCloseMHdrAndSeqFiles(rdfsp->rdfp);
		rdfsp->rdfp = ReadDBFreeSharedInfo(rdfsp->rdfp);
		rdfsp->rdfp = readdb_destruct(rdfsp->rdfp);
	}

        return;
}

/*
	Returns the ReadDBFILEPtr by the database ID.
	NULL is returned on error.
*/

ReadDBFILEPtr 
ReadDBGetDb (ReadDBFILEPtr rdfp_list, Int2 db_id)

{
	Int2 index=0;

	while (rdfp_list)
	{
		if (index == db_id)
		{
			return rdfp_list;
		}
		rdfp_list = rdfp_list->next;
		index++;
	}
	return NULL;
}

/*
	Returns the Database ID.
	-1 is returned on error.
*/

Int2 
ReadDBGetDbId (ReadDBFILEPtr list, ReadDBFILEPtr target)

{
	Int2 index=0;

	while (list)
	{
		if (readdb_compare(list, target) == TRUE)
			return index;
		list = list->next;
		index++;
	}
	return -1;
}

/*
	Formatting functions for databases formatted by formatdb.
*/
Boolean LIBCALL
PrintDbInformationBasic (CharPtr database, Boolean is_aa, Int4 line_length, CharPtr definition, Int4 number_seqs, Int8 total_length, FILE *outfp, Boolean html)

{

	asn2ff_set_output(outfp, NULL);

	ff_StartPrint(0, 0, line_length, NULL);
	if (html)
		ff_AddString("<b>Database:</b> ");
	else
		ff_AddString("Database: ");
	ff_AddString(definition);
	NewContLine();
	TabToColumn(12);
	ff_AddString(Ltostr((long) number_seqs, 1));
	ff_AddString(" sequences; ");
	ff_AddString(Nlm_Int8tostr(total_length, 1));
	ff_AddString(" total letters");
	NewContLine();
	ff_EndPrint();

	return TRUE;

}

/*
	Print a summary of the database(s) used.
*/

Boolean LIBCALL
PrintDbInformation(CharPtr database, Boolean is_aa, Int4 line_length, FILE *outfp, Boolean html)

{
	CharPtr		definition, ptr, chptr;
	Int8		total_length;
        Int4		number_seqs, length;
	ReadDBFILEPtr	rdfp, rdfp_var;
        Boolean         first_title;

	if (database == NULL || outfp == NULL)
		return FALSE;

	if (is_aa == TRUE)
		rdfp = readdb_new_ex(database, READDB_DB_IS_PROT, FALSE);
	else
		rdfp = readdb_new_ex(database, READDB_DB_IS_NUC, FALSE);

	if (rdfp == FALSE)
		return FALSE;

	rdfp_var = rdfp;
	length = 0;
	while (rdfp_var)
	{
		length += StringLen(readdb_get_title(rdfp_var));
		length += 3;
		rdfp_var = rdfp_var->next;
	}
	definition = MemNew(length*sizeof(Char));
	ptr = definition;
	rdfp_var = rdfp;
        
        first_title = TRUE;
	while (rdfp_var) {
            chptr = readdb_get_title(rdfp_var);
            
            if(chptr == NULL) {
                rdfp_var = rdfp_var->next;
                continue;
            }
            
            if(!first_title) {
                *ptr = ';';
                ptr++;
                *ptr = ' ';
                ptr++;
            }
	    else
                first_title = FALSE;
		
            
            StringCpy(ptr, readdb_get_title(rdfp_var)); 	

            length = StringLen(ptr);
            ptr += length;

            rdfp_var = rdfp_var->next;
	}

	*ptr = NULLB;
	readdb_get_totals_ex(rdfp, &(total_length), &(number_seqs), TRUE);

	rdfp = ReadDBFreeSharedInfo(rdfp);
	rdfp = readdb_destruct(rdfp);

	PrintDbInformationBasic (database, is_aa, line_length, definition, number_seqs, total_length, outfp, html);

	definition = MemFree(definition);

	return TRUE;
}

/** Common Index Stuff **/

/* Parse DB configuration file */

#define	MAX_LINE_LENGTH	1024

typedef enum {
    lexIGNORE,
    lexINT,
    lexSTRING,
    lexBOOL,
    lexEOF
} LexTokens;

static CharPtr	getline (FILE *fp, CharPtr buf)
{
    buf[0] = '\0';
    while (!buf || (buf[0] == '#') || (buf[0] == '\0') || (buf[0] == '\n')) {
	FileGets(buf, MAX_LINE_LENGTH, fp);
    }
    return buf;
}
    	    
static Int4	parseInt(CharPtr buf)
{
    Int4	retval;
    long	my_long;

    sscanf(buf, "%ld", &my_long);
    retval = my_long;

    return retval;
}

static CharPtr	parseString(CharPtr buf)
{
    CharPtr	retval = MemNew(sizeof(Char) * MAX_LINE_LENGTH);

    sscanf(buf, "%s", retval);
    return retval;
}

static Boolean	parseBool(CharPtr buf)
{
    Boolean	retval;
    CharPtr	str = parseString(buf);

    if ((!StrCmp(str, "true")) || (!StrCmp(str, "True")) || (!StrCmp(str, "TRUE")) ||
	    (!StrCmp(str, "t")) || (!StrCmp(str, "T")) ||
	    (!StrCmp(str, "1")) ||
	    (!StrCmp(str, "y")) || (!StrCmp(str, "Y")))
	retval = TRUE;
    else
	retval = FALSE;

    str = MemFree(str);
    return retval;
}

Int2	ParseDBConfigFile(DataBaseIDPtr *dbidsp, CharPtr path)
{
    Int2		number_of_DBs = 0, i;
    FILE		*fp;
    DataBaseIDPtr	retval;
    Char		buf[MAX_LINE_LENGTH], name[MAX_LINE_LENGTH];
    Char		dbid[MAX_LINE_LENGTH], isprot[MAX_LINE_LENGTH];
    Char		full_filename[PATH_MAX];

    /* open config file */
    if (path && StrCmp(path, "")) {
	sprintf(full_filename, "%s%s%s", path, DIRDELIMSTR, DB_CONFIG_FN);
    } else {
	sprintf(full_filename, "%s", DB_CONFIG_FN);
    }
    
    if (!(fp = FileOpen(full_filename, "r")))
	return 0;
    
    getline(fp, buf);

    /* first line is number of databases */
    number_of_DBs = parseInt(buf);
    
    /* allocate that much memory */
    retval = (DataBaseIDPtr) MemNew(sizeof(DataBaseID) * number_of_DBs);
    
    /* each next line is contains name, id and type of a DB */
    for (i=0; i < number_of_DBs; i++) {
	getline(fp, buf);
	sscanf(buf, "%s%s%s", name, dbid, isprot);
	(retval+i)->name   = parseString(name);
	(retval+i)->id     = parseInt(dbid);
	(retval+i)->isprot = parseBool(isprot);
    }

    FileClose(fp);
    *dbidsp = retval;
    return number_of_DBs;
}

/* ---------------------------------------------------------------------*/
/* --------- Here is set of functions, that uses in formatdb ---------- */
/* ---------------------------------------------------------------------*/

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

FormatDBPtr FormatDBInit(FDB_optionsPtr options)
{
    
    FormatDBPtr		fdbp;
    Char		filenamebuf[FILENAME_MAX];
    Uint4		i = 0;
    
    fdbp = (FormatDBPtr) MemNew (sizeof(*fdbp));
    
    fdbp->num_of_seqs = 0;
    fdbp->TotalLen=0, fdbp->MaxSeqLen=0;

    fdbp->options = options;
    
    /* If basename is NULL, use dbname. */
    if (options->base_name == NULL)
	options->base_name = options->db_file;
    
    fdbp->fd = NULL;
    fdbp->aip = NULL;
    
    /* open output BLAST files */
    
    /* Defline file */
    
    sprintf(filenamebuf, "%s.%chr", 
            options->base_name, fdbp->options->is_protein ? 'p' : 'n'); 
    fdbp->fd_def = FileOpen(filenamebuf, "wb");        
    
    /* Sequence file */
    
    sprintf(filenamebuf, "%s.%csq",
            options->base_name, fdbp->options->is_protein ? 'p' : 'n'); 
    fdbp->fd_seq = FileOpen(filenamebuf, "wb");        
    
    if (FileWrite(&i, 1, 1, fdbp->fd_seq) != (Uint4) 1) /* Sequence file started from NULLB */
	return NULL;
    
    /* Index file */
    
    sprintf(filenamebuf, "%s.%cin",
            options->base_name, fdbp->options->is_protein ? 'p' : 'n'); 
    fdbp->fd_ind = FileOpen(filenamebuf, "wb");      

    /* Misc. info dump file */

    if(options->dump_info) {
        sprintf(filenamebuf, "%s.%cdi",
                options->base_name, fdbp->options->is_protein ? 'p' : 'n'); 
        fdbp->fd_sdi = FileOpen(filenamebuf, "wb"); 
    }

    /* String (accession) index temporary file */
    
    fdbp->fd_stmp = NULL;

    if(options->parse_mode) {
        sprintf(filenamebuf, "%s.%ctm",
                options->base_name, fdbp->options->is_protein ? 'p' : 'n'); 
        fdbp->fd_stmp = FileOpen(filenamebuf, "wb");      
    }
    
    ErrLogPrintf("Version %s [%s]\n", BlastGetVersionNumber(), BlastGetReleaseDate()); 
    ErrLogPrintf("Started database file \"%s\"\n", options->db_file);
    
    /* Allocating space for offset tables */
    fdbp->DefOffsetTable = (Int4Ptr)MemNew(fdbp->OffsetAllocated*sizeof(Uint4));
    fdbp->SeqOffsetTable = (Int4Ptr)MemNew(fdbp->OffsetAllocated*sizeof(Uint4));
    if(!options->is_protein) 
        fdbp->AmbOffsetTable = (Int4Ptr)MemNew(fdbp->OffsetAllocated*sizeof(Uint4));
    else
        fdbp->AmbOffsetTable = NULL;
    
    /* Allocating space for lookup table */
    
    if((fdbp->lookup = FASTALookupNew()) == NULL) {
        ErrLogPrintf("Error initializing Lookup structure. Formatting failed.\n");
        return NULL;
    }
    
    return fdbp;
}

#ifdef REDUCED_E2INDEX_SET
/*****************************************************************************
*
*   SeqIdE2Index(anp)
*   	atp is the current type (if identifier of a parent struct)
*       if atp == NULL, then assumes it stands alone (SeqId ::=)
*
*****************************************************************************/
static Boolean SeqIdE2Index (SeqIdPtr anp, FILE *fd, Int4 seq_num, 
                             Boolean sparse)
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
    Int2 version = 0;

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
    case SEQID_OTHER:     /* other */
        tsip = (TextSeqIdPtr)(anp->data.ptrvalue);
	if ((tsip->version > 0) && (tsip->release == NULL))
		version = tsip->version;
	break;
    case SEQID_PIR:       /* pir   */
    case SEQID_SWISSPROT: /* swissprot */
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
    
    if(tsip == NULL && !sparse) {
        SeqIdWrite(anp, buf, PRINTID_FASTA_SHORT, 80);
        StringLower(buf);
        fprintf(fd, "%s%c%d\n", buf, ISAM_DATA_CHAR, seq_num);
    }

    if (ptr != NULL) {   /* write a single string */
	StringMove(buf, ptr);
        StringLower(buf);
        fprintf(fd, "%s%c%ld\n", buf, ISAM_DATA_CHAR, (long) seq_num);
 
        if (chain != 0) { /* PDB only. */
            fprintf(fd, "%s|%c%c%ld\n", buf, chain, ISAM_DATA_CHAR, 
                    (long) seq_num);
            fprintf(fd, "%s %c%c%ld\n", buf, chain, ISAM_DATA_CHAR, 
                    (long) seq_num);
        }
    }

    if (tsip != NULL) {   /* separately index accession and locus */
        if ((tsip->accession != NULL) && (tsip->name != NULL) && !sparse) {

            /* Here is we indexing accession as part of SeqId */
            tmp = tsip->name;
            tsip->name = NULL;
            SeqIdWrite(anp, buf, PRINTID_FASTA_SHORT, 80);
            StringLower(buf);
            fprintf(fd, "%s%c%d\n", buf, ISAM_DATA_CHAR, seq_num);
            tsip->name = tmp;

            /* If accession and locus are different we also print locus */
            if(StringICmp(tsip->accession, tsip->name)) {
                
                tmp = tsip->accession;
                tsip->accession = NULL;
                SeqIdWrite(anp, buf, PRINTID_FASTA_SHORT, 80);
                StringLower(buf);
                fprintf(fd, "%s%c%d\n", buf, ISAM_DATA_CHAR, seq_num);
                tsip->accession = tmp;
            }

	    if (version) { /* Index accession without verison. */
		tsip->version = 0;
            	tmp = tsip->name;
            	tsip->name = NULL;
            	SeqIdWrite(anp, buf, PRINTID_FASTA_SHORT, 80);
                StringLower(buf);
            	fprintf(fd, "%s%c%d\n", buf, ISAM_DATA_CHAR, seq_num);
            	tsip->name = tmp;
		tsip->version = version;
	    }
        }

               /* now index as separate strings */
        if (tsip->accession != NULL) {
            StringMove(buf, tsip->accession);
            StringLower(buf);
            fprintf(fd, "%s%c%d\n", buf, ISAM_DATA_CHAR, seq_num);
	    if (version && !sparse) {
            	fprintf(fd, "%s%.d%c%d\n", buf, version, 
                        ISAM_DATA_CHAR, seq_num);
            }
	}

        if(tsip->name != NULL && 
           StringICmp(tsip->accession, tsip->name) && !sparse) {
            StringMove(buf, tsip->name);
            StringLower(buf);
            fprintf(fd, "%s%c%d\n", buf, ISAM_DATA_CHAR, seq_num);
	}
    }
    
    if (do_gb && !sparse) {   /* index embl and ddbj as genbank */
        tmptype = anp->choice;
        anp->choice = SEQID_GENBANK;
        SeqIdE2Index(anp, fd, seq_num, sparse);
        anp->choice = tmptype;
    }

    retval = TRUE;
    return retval;
}
#else
/*****************************************************************************
*
*   SeqIdE2Index(anp)
*   	atp is the current type (if identifier of a parent struct)
*       if atp == NULL, then assumes it stands alone (SeqId ::=)
*
*****************************************************************************/
static Boolean SeqIdE2Index (SeqIdPtr anp, FILE *fd, Int4 seq_num, 
                             Boolean sparse)
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
    Int2 version = 0;

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
    case SEQID_OTHER:     /* other */
        tsip = (TextSeqIdPtr)(anp->data.ptrvalue);
	if ((tsip->version > 0) && (tsip->release == NULL))
		version = tsip->version;
	break;
    case SEQID_PIR:       /* pir   */
    case SEQID_SWISSPROT: /* swissprot */
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
    
    if(!sparse) {
        SeqIdWrite(anp, buf, PRINTID_FASTA_SHORT, 80);
        
        length = StringLen(buf);
        for(i = 0; i < length; i++)
            buf[i] = TO_LOWER(buf[i]);

        fprintf(fd, "%s%c%ld\n", buf, ISAM_DATA_CHAR, seq_num);
        
        /* Index without version. */
        if (version) {
            tsip->version = 0;
            SeqIdWrite(anp, buf, PRINTID_FASTA_SHORT, 80);
            
            length = StringLen(buf);
            for(i = 0; i < length; i++)
       		buf[i] = TO_LOWER(buf[i]);
            fprintf(fd, "%s%c%ld\n", buf, ISAM_DATA_CHAR, seq_num);
            tsip->version = version;
        }
    } /* if(!sparse) */

    if (ptr != NULL) {   /* write a single string */
	StringMove(buf, ptr);
        length = StringLen(buf);
        for(i = 0; i < length; i++)
            buf[i] = TO_LOWER(buf[i]);
        fprintf(fd, "%s%c%ld\n", buf, ISAM_DATA_CHAR, (long) seq_num);
 
        if (chain != 0) { /* PDB only. */
            fprintf(fd, "%s|%c%c%ld\n", buf, chain, ISAM_DATA_CHAR, 
                    (long) seq_num);
            fprintf(fd, "%s %c%c%ld\n", buf, chain, ISAM_DATA_CHAR, 
                    (long) seq_num);
        }
    }

    if (tsip != NULL) {   /* separately index accession and locus */
        if ((tsip->accession != NULL) && (tsip->name != NULL) && !sparse) {
            tmp = tsip->accession;
            tsip->accession = NULL;
            SeqIdWrite(anp, buf, PRINTID_FASTA_SHORT, 80);
            length = StringLen(buf);
            for(i = 0; i < length; i++)
                buf[i] = TO_LOWER(buf[i]);
            fprintf(fd, "%s%c%ld\n", buf, ISAM_DATA_CHAR, seq_num);
            tsip->accession = tmp;
            tmp = tsip->name;
            tsip->name = NULL;
            SeqIdWrite(anp, buf, PRINTID_FASTA_SHORT, 80);
            length = StringLen(buf);
            for(i = 0; i < length; i++)
                buf[i] = TO_LOWER(buf[i]);
            fprintf(fd, "%s%c%ld\n", buf, ISAM_DATA_CHAR, seq_num);
            tsip->name = tmp;
	    if (version)
	    { /* Index accession without verison. */
		tsip->version = 0;
            	tmp = tsip->name;
            	tsip->name = NULL;
            	SeqIdWrite(anp, buf, PRINTID_FASTA_SHORT, 80);
            	length = StringLen(buf);
            	for(i = 0; i < length; i++)
               		 buf[i] = TO_LOWER(buf[i]);
            	fprintf(fd, "%s%c%ld\n", buf, ISAM_DATA_CHAR, seq_num);
            	tsip->name = tmp;
		tsip->version = version;
	    }
        }

               /* now index as separate strings */
	if (tsip->name != NULL && !sparse) {
            StringMove(buf, tsip->name);
            length = StringLen(buf);
            for(i = 0; i < length; i++)
                buf[i] = TO_LOWER(buf[i]);
            fprintf(fd, "%s%c%ld\n", buf, ISAM_DATA_CHAR, seq_num);
	}
        if (tsip->accession != NULL) {
            StringMove(buf, tsip->accession);
            length = StringLen(buf);
            for(i = 0; i < length; i++)
                buf[i] = TO_LOWER(buf[i]);
            fprintf(fd, "%s%c%ld\n", buf, ISAM_DATA_CHAR, seq_num);
	    if (version && !sparse)
            	fprintf(fd, "%s.%d%c%ld\n", buf, version, ISAM_DATA_CHAR, seq_num);
	}
        
    }
    
    if (do_gb && !sparse) {   /* index embl and ddbj as genbank */
        tmptype = anp->choice;
        anp->choice = SEQID_GENBANK;
        SeqIdE2Index(anp, fd, seq_num, sparse);
        anp->choice = tmptype;
    }

    retval = TRUE;
    return retval;
}
#endif

/*****************************************************************************
 *
 *   SeqIdSetE2Index(anp, e2p, settype, elementtype)
 *
 *****************************************************************************/
static Boolean SeqIdSetE2Index (SeqIdPtr anp, FILE *fd, Int4 seq_num, 
                                Boolean sparse)
{
    SeqIdPtr oldanp;
    Boolean retval = FALSE;
    
    if (anp == NULL)
        return FALSE;
    
    oldanp = anp;
    
    while (anp != NULL) {
        if (!SeqIdE2Index(anp, fd, seq_num, sparse))
            goto erret;
        anp = anp->next;
    }
    
    retval = TRUE;
erret:
    return retval;
}
static SeqIdPtr SeqIdSetFree_NO_OBJ_MGR(SeqIdPtr sip)
{
	SeqIdPtr	next;
	
	while(sip != NULL){
		next=sip->next;
		switch(sip->choice) {
		 case SEQID_LOCAL:      /* local */
			ObjectIdFree(sip->data.ptrvalue);
			break;
		 case SEQID_GIBBSQ:      /* gibbseq */
		 case SEQID_GIBBMT:      /* gibbmt */
			break;
		 case SEQID_GIIM:      /* giimid */
			GiimFree(sip->data.ptrvalue);
			break;
		 case SEQID_GENBANK:      /* genbank */
		 case SEQID_EMBL:      /* embl */
		 case SEQID_PIR:      /* pir   */
		 case SEQID_SWISSPROT:      /* swissprot */
		 case SEQID_OTHER:     /* other */
		 case SEQID_DDBJ:
		 case SEQID_PRF:
			TextSeqIdFree(sip->data.ptrvalue);
			break;
		 case SEQID_PATENT:      /* patent seq id */
			PatentSeqIdFree(sip->data.ptrvalue);
			break;
		 case SEQID_GENERAL:     /* general */
			DbtagFree(sip->data.ptrvalue);
			break;
		 case SEQID_GI:     /* gi */
			break;
		 case SEQID_PDB:
                        PDBSeqIdFree(sip->data.ptrvalue);
                        break;
		}
		MemFree(sip);
		sip=next;
	}
	return NULL;
}

static Int4 UpdateLookupInfo(CharPtr defline, 
                             FASTALookupPtr lookup, 
                             Int4 num_of_seqs,
                             FILE *fd_stmp,
                             Boolean ParseSeqid,
                             Boolean sparse
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
        
        for(i=0; !isspace((int)*p) && i < ID_MAX_SIZE; p++,i++)
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
        
        if(!SeqIdSetE2Index (sip, fd_stmp, num_of_seqs, sparse)) {
            ErrLogPrintf("SeIdSetE2Index failed. Exiting..\n");
            return FALSE;
        }
        
	sip = SeqIdSetFree_NO_OBJ_MGR(sip);
        
        if((p = StringChr(d, READDB_DEF_SEPARATOR)) == NULL)
            break;
        else
            p++;
    }
    return LOOKUP_NO_ERROR;
}
static Boolean FormatdbCreateStringIndex(const CharPtr FileName, 
                                         Boolean ProteinType,
                                         Int4 sparse_idx)
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
    
    if((error = ISAMMakeIndex(isamp, 0, sparse_idx)) != ISAMNoError) {
        ErrPostEx(SEV_ERROR, 0, 0, "Creating of index failed with error code %ld\n", (long) error);
        ISAMObjectFree(isamp);
        return FALSE;
    } 
    
    ISAMObjectFree(isamp);
    
    return TRUE;
}

Int2 FDBAddSequence (FormatDBPtr fdbp, Int4 gi, CharPtr seq_id,
                     CharPtr title, Int4 tax_id, CharPtr div, 
                     Int4 owner, Uint1 seq_data_type, 
                     ByteStorePtr *seq_data, Int4 SequenceLen, Int4 date)
{
    
    CharPtr		defline;
    Char		tmpbuff[1025];
    Int4		buffer_size=0, defline_len=0;
    CharPtr		buffer=NULL;
    Int4		len, id_length;
    Uint4Ptr		AmbCharPtr = NULL;
    Uint1		ch, remainder, new_code;
    Uint4		i, total, index, hash;
    ByteStorePtr        new_data;
    
    if (fdbp->options->is_protein && seq_data_type != Seq_code_ncbistdaa) {
        new_data =  BSConvertSeq (*seq_data, Seq_code_ncbistdaa, 
                                  seq_data_type, SequenceLen);
        *seq_data = new_data;
    }

    /* BioseqRawConvert(bsp, Seq_code_ncbistdaa); */
    
    fdbp->TotalLen += SequenceLen;
    
    if (fdbp->MaxSeqLen < SequenceLen)
        fdbp->MaxSeqLen = SequenceLen;
    
    if(fdbp->OffsetAllocated <= (fdbp->num_of_seqs+1)) {
        fdbp->OffsetAllocated += INDEX_ARRAY_CHUNKS;
        
        fdbp->DefOffsetTable = (Int4Ptr)Realloc(fdbp->DefOffsetTable, 
                                                fdbp->OffsetAllocated*sizeof(Uint4));
        fdbp->SeqOffsetTable = (Int4Ptr)Realloc(fdbp->SeqOffsetTable, 
                                                fdbp->OffsetAllocated*sizeof(Uint4));
        if(!fdbp->options->is_protein) {
            fdbp->AmbOffsetTable = (Int4Ptr)Realloc(fdbp->AmbOffsetTable, 
                                                    fdbp->OffsetAllocated*sizeof(Uint4));
        }
    }
    
    fdbp->DefOffsetTable[fdbp->num_of_seqs] = ftell(fdbp->fd_def); 
    fdbp->SeqOffsetTable[fdbp->num_of_seqs] = ftell(fdbp->fd_seq);
    
    /* ---------------------- */
    
    if(fdbp->options->parse_mode == FALSE)  {
        sprintf(tmpbuff, "%s%ld ", NON_SEQID_PREFIX, fdbp->num_of_seqs);
        if (FileWrite(tmpbuff, StringLen(tmpbuff), 1, fdbp->fd_def) != (Uint4) 1)
            return 1;
        defline = title;
    } else {
        if (title != NULL)
            defline_len = StringLen(title);
        else
            defline_len = 0;
        
        defline_len += 255;	/* Sufficient for an ID. */

        if (buffer_size < defline_len) {
            if (buffer)
                buffer = MemFree(buffer);
            buffer = MemNew((defline_len+1)*sizeof(Char));
            buffer_size = defline_len;
        }

	/* IF the gi is zero and there is another ID, then do not print it. */
	if (StringNCmp(seq_id, "gi|0|", 5) == 0)
	{
		StringCpy(buffer, seq_id+5);	
            	ErrPostEx(SEV_WARNING, 0, 0, "%s: zero gi stripped", seq_id);
	}
	else
		StringCpy(buffer, seq_id);	
		
        id_length = StringLen(buffer);
        buffer[id_length] = ' ';
        id_length++;
        StringCpy(&buffer[id_length], title);
        defline = buffer;
    }
    if (FileWrite(defline, StringLen(defline), 1, fdbp->fd_def) != (Uint4) 1)
	return 1;
    
    /* -------- Now adding new entried into lookup hash table */
    
    if((UpdateLookupInfo(defline, fdbp->lookup, fdbp->num_of_seqs, 
                         fdbp->fd_stmp, fdbp->options->parse_mode, 
                         fdbp->options->sparse_idx)) != LOOKUP_NO_ERROR) {
        return -1;
    }
    
    defline = NULL;
    if (buffer)
	MemFree(buffer);
    
    if(!fdbp->options->is_protein)  {
        AmbCharPtr = NULL;
        if(seq_data_type != Seq_code_ncbi2na && seq_data_type != Seq_code_ncbi4na){
            new_data=BSPack(*seq_data,seq_data_type,SequenceLen,&new_code);
            if(new_data != NULL){
                *seq_data=new_data;
                seq_data_type=new_code;
            }
        }

        if (seq_data_type == Seq_code_ncbi4na && seq_data != NULL){
            
            /* ncbi4na require compression into ncbi2na */
            
            if((new_data = BSCompressDNA(*seq_data, SequenceLen, 
                                         &AmbCharPtr)) == NULL) {
                ErrLogPrintf("Error converting ncbi4na to ncbi2na. " 
                             "Formating failed.\n");
                return -1;
            }
            *seq_data = new_data;
            
            seq_data_type = Seq_code_ncbi2na; /* just for information */
            
        } else {
            /* if sequence already in ncbi2na format we have to update last byte */
            BSSeek(*seq_data, SequenceLen/4, SEEK_SET);
            
            if((remainder = (SequenceLen%4)) == 0) {
                BSPutByte(*seq_data, NULLB);
            } else {
                ch = remainder + (BSGetByte(*seq_data) & 0xfc);
                BSSeek(*seq_data, SequenceLen/4, SEEK_SET);
                BSPutByte(*seq_data, ch);
            }
        }
    }

    /* Now dumping sequence */
    
    BSSeek(*seq_data, 0, SEEK_SET);

    hash = 0;
    
    while((len = BSRead(*seq_data, tmpbuff, sizeof(tmpbuff)-1)) != 0) {
        if (FileWrite(tmpbuff, len, 1, fdbp->fd_seq) != (Uint4) 1)
            return 1;
        if(fdbp->options->dump_info){
            int i;
            for(i=0;i<len;i++){
                hash *= 1103515245;
                hash += (unsigned long)(tmpbuff[i]) + 12345;
            }
        }
    }
    
    if(fdbp->options->is_protein) {
        i=0;
        if (FileWrite(&i, 1, 1, fdbp->fd_seq) != (Uint4) 1)
            return 1;
    } else {
        /* dump ambiguity characters. */
        fdbp->AmbOffsetTable[fdbp->num_of_seqs] = ftell(fdbp->fd_seq); /* Anyway... */
        
        /* if AmbCharPtr is not NULL, then there was ambiguity. */
        if(AmbCharPtr != NULL) { /* The first Uint4 holds the total number of ambig. bp. */
            total = (*AmbCharPtr)+1;
            for (index=0; index<total; index++) {
                if (!FormatDbUint4Write(AmbCharPtr[index], fdbp->fd_seq))
                    return 1;
            }
            MemFree(AmbCharPtr);
            AmbCharPtr = NULL;
        }
    }


    /* Dumping misc info file now */
    if(fdbp->options->dump_info) {               
        fprintf(fdbp->fd_sdi, "%ld %ld %ld %ld %s %ld %ld %ld\n", 
                fdbp->num_of_seqs, gi, tax_id, owner, div?div:"N/A", 
                SequenceLen, hash, date);
    }

    fdbp->num_of_seqs++;  /* Finshed ... */
    
    return 0;
}    
Int2 FDBAddBioseq(FormatDBPtr fdbp, BioseqPtr bsp)
{
    Char tmpbuf[128];
    SeqIdWrite(bsp->id, tmpbuf, PRINTID_FASTA_LONG, sizeof(tmpbuf));

    return FDBAddSequence (fdbp, 0, tmpbuf, BioseqGetTitle(bsp), 
                           0, 0, 0, bsp->seq_data_type, &bsp->seq_data, 
                           bsp->length, 0);
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
Int2 process_sep (SeqEntryPtr sep, FormatDBPtr fdbp)
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
    
    if (fdbp->options->is_protein)
        BioseqRawConvert(bsp, Seq_code_ncbistdaa);
    
    SequenceLen = bsp->length;
    fdbp->TotalLen += SequenceLen;
    
    if (fdbp->MaxSeqLen < SequenceLen)
        fdbp->MaxSeqLen = SequenceLen;
    
    if(fdbp->OffsetAllocated <= (fdbp->num_of_seqs+1)) {
        fdbp->OffsetAllocated += INDEX_ARRAY_CHUNKS;
        
        fdbp->DefOffsetTable = (Int4Ptr)Realloc(fdbp->DefOffsetTable, 
                                                fdbp->OffsetAllocated*sizeof(Uint4));
        fdbp->SeqOffsetTable = (Int4Ptr)Realloc(fdbp->SeqOffsetTable, 
                                                fdbp->OffsetAllocated*sizeof(Uint4));
        if(!fdbp->options->is_protein) {
            fdbp->AmbOffsetTable = (Int4Ptr)Realloc(fdbp->AmbOffsetTable, 
                                                    fdbp->OffsetAllocated*sizeof(Uint4));
        }
    }
    
    fdbp->DefOffsetTable[fdbp->num_of_seqs] = ftell(fdbp->fd_def); 
    fdbp->SeqOffsetTable[fdbp->num_of_seqs] = ftell(fdbp->fd_seq);
    
    /* ---------------------- */
    
    if(fdbp->options->parse_mode == FALSE)  {
        sprintf(tmpbuff, "%s%ld ", NON_SEQID_PREFIX, fdbp->num_of_seqs);
        if (FileWrite(tmpbuff, StringLen(tmpbuff), 1, fdbp->fd_def) != (Uint4) 1)
            return 1;
        defline = (CharPtr)bsp->descr->data.ptrvalue;
    } else {
        if (bsp->descr)
            defline_len = StringLen(BioseqGetTitle(bsp));
        else
            defline_len = 0;
        defline_len += 255;	/* Sufficient for an ID. */
        if (buffer_size < defline_len) {
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
    
    if((UpdateLookupInfo(defline, fdbp->lookup, fdbp->num_of_seqs, 
                         fdbp->fd_stmp, fdbp->options->parse_mode, 
                         fdbp->options->sparse_idx)) != LOOKUP_NO_ERROR) {
        return -1;
    }
    
    defline = NULL;
    if (buffer)
	MemFree(buffer);
    
    if(!fdbp->options->is_protein)  {
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
    
    while((len = BSRead(bsp->seq_data, tmpbuff, sizeof(tmpbuff))) != 0) {
        if (FileWrite(tmpbuff, len, 1, fdbp->fd_seq) != (Uint4) 1)
            return 1;
    }
 
    
    if(fdbp->options->is_protein) {
        i=0;
        if (FileWrite(&i, 1, 1, fdbp->fd_seq) != (Uint4) 1)
            return 1;
    } else {
        /* dump ambiguity characters. */
        fdbp->AmbOffsetTable[fdbp->num_of_seqs] = ftell(fdbp->fd_seq); /* Anyway... */
        
        /* if AmbCharPtr is not NULL, then there was ambiguity. */
        if(AmbCharPtr != NULL) { /* The first Uint4 holds the total number of ambig. bp. */
            total = (*AmbCharPtr)+1;
            for (index=0; index<total; index++) {
                if (!FormatDbUint4Write(AmbCharPtr[index], fdbp->fd_seq))
                    return 1;
            }
            MemFree(AmbCharPtr);
            AmbCharPtr = NULL;
        }
    }
    
    fdbp->num_of_seqs++;  /* Finshed ... */

    return 0;
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

/*******************************************************************************
 * Finish stage - out offset tables, etc, into files.  Is to be called before
 * FormatDBClose()
 ******************************************************************************* 
 * Parameters:
 *	
 *
 * Returns  void
 ******************************************************************************/
#define DATETIME_LENGTH 64

static	Int2	FDBFinish (FormatDBPtr fdbp) 
{
    Char	DBName[FILENAME_MAX];
    Int4	title_len;
    Char	dateTime[DATETIME_LENGTH];
    ISAMObjectPtr object;
    ISAMErrorCode error;
    Uint4	i;
    Char	filenamebuf[FILENAME_MAX];
    Int2	tmp, extra_bytes = 0;

    if(fdbp->num_of_seqs==0){
	ErrLogPrintf("FDBFinish: Empty %s database...\n",
                     fdbp->options->is_protein?"protein":"nucleotide");
	return 0;
    }
    
   
    fdbp->DefOffsetTable[fdbp->num_of_seqs] = ftell(fdbp->fd_def); 
    
    if(!fdbp->options->is_protein) {
        fdbp->AmbOffsetTable[fdbp->num_of_seqs] = ftell(fdbp->fd_seq);
        fdbp->SeqOffsetTable[fdbp->num_of_seqs] =
            fdbp->AmbOffsetTable[fdbp->num_of_seqs];
    } else {
        fdbp->SeqOffsetTable[fdbp->num_of_seqs] = ftell(fdbp->fd_seq);
    }
    
        /* Parsing finished - now dumping index file */
    
    if(fdbp->options->parse_mode)
        FileClose(fdbp->fd_stmp);
    
        /* Information */
    
    if (!FormatDbUint4Write(FORMATDB_VER, fdbp->fd_ind))
	return 1;
    if (!FormatDbUint4Write(fdbp->options->is_protein, fdbp->fd_ind))
	return 1;
    
    if(fdbp->options->db_title != NULL) {
        title_len = StringLen(fdbp->options->db_title);
    }
    else
        title_len = 0;
    
    if (!FormatDbUint4Write(title_len, fdbp->fd_ind))
	return 1;
    
    if (title_len != 0)
        if (FileWrite(fdbp->options->db_title, title_len, 1, fdbp->fd_ind) != (Uint4) 1)
		return 1;
    
    MemSet(dateTime, 0, DATETIME_LENGTH);
    Nlm_DayTimeStr(dateTime, TRUE, TRUE);

    /* write db_title and date-time stamp eigth bytes aligned */
    tmp = title_len + StringLen(dateTime);
    if (tmp%8) {
	extra_bytes = 8 - tmp%8;
    }
    if (!FormatDbUint4Write(StringLen(dateTime) + extra_bytes, fdbp->fd_ind))
	return 1;
    if (FileWrite(dateTime, StringLen(dateTime) + extra_bytes, 1, fdbp->fd_ind) != 1)
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
    if(!fdbp->options->is_protein) {
        for(i=0; i <= fdbp->num_of_seqs; i++) {
            if (!FormatDbUint4Write(fdbp->AmbOffsetTable[i], fdbp->fd_ind))
		return 1;
        }
    }
    
        /* Numeric lookup table sort & dump */
    
    if(fdbp->options->parse_mode && fdbp->lookup->used > 0) {
        FILE	*fd_lookup;
        
        sprintf(DBName, "%s.%cnd", fdbp->options->base_name, 
                fdbp->options->is_protein ? 'p' : 'n'); 
        
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
                fdbp->options->base_name, fdbp->options->is_protein ? 'p' : 'n'); 
        
        if((object = ISAMObjectNew(ISAMNumeric, 
                                   DBName, filenamebuf)) == NULL) {
            ErrPostEx(SEV_ERROR, 0, 0, "Failed to create ISAM object.\n");
            return 1;
        }
        
        if((error = ISAMMakeIndex(object, 0, 0)) != ISAMNoError) {
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
    
    if(fdbp->options->parse_mode)
    {
        if (!FormatdbCreateStringIndex(fdbp->options->base_name, 
                                       fdbp->options->is_protein,
                                       fdbp->options->sparse_idx))
            return 1;
    }
    
    ErrLogPrintf("Formated %d sequences\n", fdbp->num_of_seqs);

    return 0;
    
} /* end FDBFinish() */

/*******************************************************************************
 * Free memory allocated for given variable of FormatDB
 ******************************************************************************* 
 * Parameters:
 *	fdbp	- pointer to memory to be freed
 *
 * Returns NULL
 ******************************************************************************/

Int2 FormatDBClose(FormatDBPtr fdbp)
{

    /* Now dumping all data to disk */
    
    if(FDBFinish (fdbp))
        return 1;
    
    /* ... and MemFree all stuff */

    MemFree(fdbp->DefOffsetTable);
    MemFree(fdbp->SeqOffsetTable);
    
    if(!fdbp->options->is_protein) {
        MemFree(fdbp->AmbOffsetTable);
    }
    
    FASTALookupFree(fdbp->lookup);
    
    if (fdbp->fd)
        FileClose(fdbp->fd);

    FileClose(fdbp->fd_def);
    FileClose(fdbp->fd_ind);
    FileClose(fdbp->fd_seq);
    FileClose(fdbp->fd_sdi);
    
    if (fdbp->aip)
        AsnIoClose(fdbp->aip);
    
    MemFree (fdbp);
    
    return 0;
}
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
    if (group_segs == 3) { /* do 2 things */
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
Boolean BLASTFileFunc (BioseqPtr bsp, Int2 key, CharPtr buf, Uint4 buflen,
                       Pointer data)
{
    FormatDBPtr	fdbp = (FormatDBPtr) data;
    Int4		SequenceLen;
    Uint4		i, total, index;
    
    switch (key) {
    case FASTA_ID:
        
        SequenceLen = bsp->length;
        fdbp->TotalLen += SequenceLen;
        
        if (fdbp->MaxSeqLen < SequenceLen)
            fdbp->MaxSeqLen = SequenceLen;
        
        if(fdbp->OffsetAllocated <= fdbp->num_of_seqs) {
            fdbp->OffsetAllocated += INDEX_ARRAY_CHUNKS;
            
            fdbp->DefOffsetTable = (Int4Ptr)Realloc(fdbp->DefOffsetTable, 
                                                    fdbp->OffsetAllocated*sizeof(Uint4));
            fdbp->SeqOffsetTable = (Int4Ptr)Realloc(fdbp->SeqOffsetTable, 
                                                    fdbp->OffsetAllocated*sizeof(Uint4));
            if(!fdbp->options->is_protein) {
                fdbp->AmbOffsetTable = (Int4Ptr)Realloc(fdbp->AmbOffsetTable, 
                                                        fdbp->OffsetAllocated*sizeof(Uint4));
            }
        }
        
        fdbp->DefOffsetTable[fdbp->num_of_seqs] = ftell(fdbp->fd_def); 
        fdbp->SeqOffsetTable[fdbp->num_of_seqs] = ftell(fdbp->fd_seq);
        
        if (FileWrite(buf, buflen, 1, fdbp->fd_def) != (Uint4) 1)
            return FALSE;
        if (FileWrite(" ", 1, 1, fdbp->fd_def) != (Uint4) 1)
            return FALSE;
        
        /* Now adding new entried into lookup hash table */
        
        if((UpdateLookupInfo(buf, fdbp->lookup, fdbp->num_of_seqs, 
                             fdbp->fd_stmp, fdbp->options->parse_mode,
                             fdbp->options->sparse_idx)) != LOOKUP_NO_ERROR) {
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
        if(fdbp->options->is_protein) {
            i=0;
            if (FileWrite(&i, 1, 1, fdbp->fd_seq) != (Uint4) 1)
                return FALSE;
        } else {
            /* dump ambiguity characters. */
            fdbp->AmbOffsetTable[fdbp->num_of_seqs] = ftell(fdbp->fd_seq); /* Anyway... */
            
            /* if AmbCharPtr is not NULL, then there was ambiguity. */
            if(fdbp->AmbCharPtr != NULL) {
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
    case FASTA_FORMATDB_AMB: {
        Int4 len;
        Char tmpbuff[1024];
        /* In case of "formatdb" nucleotides have to be compressed */
        
        fdbp->AmbCharPtr = NULL;
        
        if (bsp->seq_data_type == Seq_code_ncbi4na && bsp->seq_data != NULL){
                
            /* ncbi4na require compression into ncbi2na */
            
            if((bsp->seq_data = BSCompressDNA(bsp->seq_data, bsp->length, 
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

/* ----------------- Proccessing ASN.1 with formatdb ----------------- */

typedef struct _FDB_SEDataInfo {
    CharPtr         seqid;
    ByteStorePtr    bsp;
    Int4            length;
    Uint1           seq_data_type;
    CharPtr         defline;
    FastaDat PNTR   tfp;
    FormatDBPtr     fdbp;
} FDB_SEDataInfo, PNTR FDB_SEDataInfoPtr;

static Boolean FDB_FastaFileFunc(BioseqPtr bsp, Int2 key, CharPtr buf, 
                                 Uint4 buflen, Pointer data)
{
    FDB_SEDataInfoPtr fsedip;

    if((fsedip = data) == NULL)
        return TRUE;
    
    switch (key) {
    case FASTA_DEFLINE:
        MemCpy(fsedip->defline, buf, buflen);
        fsedip->defline[buflen] = NULLB;
        break;
    case FASTA_SEQLINE:
        BSWrite(fsedip->bsp, buf, buflen);
        fsedip->length += buflen;
        break;
    case FASTA_ID:
        MemCpy(fsedip->seqid, buf, buflen);
        fsedip->seqid[buflen] = NULLB;
        break;
    case FASTA_EOS:   /* end of sequence */
        /* Here we should add new entry to FD database and reset 
           all spaces */

        FDBAddSequence(fsedip->fdbp, 0 ,fsedip->seqid, fsedip->defline,
                       0, 0, 0, fsedip->seq_data_type, &fsedip->bsp, 
                       fsedip->length, 0);
        
        BSSeek(fsedip->bsp, 0, SEEK_SET);
        BSDelete(fsedip->bsp, BSLen(fsedip->bsp));
        fsedip->length = 0;

        break;
    }

    return TRUE;
}

FDB_SEDataInfoPtr FDB_SEDataInfoNew(void)
{
    FDB_SEDataInfoPtr fsedip;
    MyFsa PNTR mfp;

    fsedip = MemNew(sizeof(FDB_SEDataInfo));

    fsedip->tfp  = MemNew(sizeof (FastaDat));
    mfp  = MemNew(sizeof (MyFsa));
    fsedip->tfp->mfp = mfp;

    mfp->buf     = MemNew(255);
    mfp->buflen  = 254;
    mfp->seqlen  = 254;
    mfp->myfunc  = FDB_FastaFileFunc;
    mfp->bad_asn1        = FALSE;
    mfp->order           = 0;
    mfp->accession       = NULL;
    mfp->organism        = NULL;
    mfp->do_virtual      = TRUE;
    mfp->tech            = 0;
    mfp->no_sequence     = FALSE;
    mfp->formatdb        = FALSE;
    mfp->mydata          = fsedip; /* ... */

    fsedip->tfp->group_segs = 1; /*** to trigger delta's and maps ***/
    fsedip->tfp->last_indent = -1;
    fsedip->tfp->parts = -1;
    fsedip->tfp->seg = -1;
    fsedip->tfp->got_one = FALSE;


    if(fsedip->seqid == NULL)
        fsedip->seqid = MemNew(fsedip->tfp->mfp->buflen+1);
    
    fsedip->bsp = BSNew(2048); 

    if(fsedip->defline == NULL){
        fsedip->defline = MemNew(fsedip->tfp->mfp->buflen+1);
    }    
    return fsedip;
}

void FDB_SEDataInfoFree(FDB_SEDataInfoPtr fsedip)
{
    MemFree(fsedip->tfp->mfp->buf);
    MemFree(fsedip->tfp->mfp);
    MemFree(fsedip->tfp);
    BSFree(fsedip->bsp);
    MemFree(fsedip->defline);
    MemFree(fsedip->seqid);

    MemFree(fsedip);
}

static void FDBSeqEntry_callback (SeqEntryPtr sep, Pointer data, 
                                  Int4 index, Int2 indent)
{
    FDB_SEDataInfoPtr fsedip;
    BioseqPtr    bsp=NULL;
    Boolean      is_na;
    
    if((fsedip = (FDB_SEDataInfoPtr) data) == NULL)
        return;
    
    if(!IS_Bioseq(sep)) {
        SeqEntryFasta(sep, fsedip->tfp, index, indent);
        return;
    }
    
    bsp = sep->data.ptrvalue;
    is_na = ISA_na(bsp->mol);
    
    /* We will format only sequences of one kind */
    if(fsedip->fdbp->options->is_protein != !is_na) {
        fsedip->tfp->mfp->no_sequence = TRUE;
        SeqEntryFasta(sep, fsedip->tfp, index, indent);
        return;
    }

    /* Segmented and virtual sequences are not indexed */
    if(bsp->repr == Seq_repr_seg || bsp->repr == Seq_repr_virtual) {
        fsedip->tfp->mfp->no_sequence = TRUE;
        SeqEntryFasta(sep, fsedip->tfp, index, indent);
        return;
    }

    fsedip->tfp->last_indent = -1;
    
    if(bsp->repr == Seq_repr_raw || bsp->repr == Seq_repr_const){
        
        
        /* This will collect defline and seqid */

        fsedip->tfp->mfp->no_sequence = TRUE;
        SeqEntryFasta(sep, fsedip->tfp, index, indent);
        
        FDBAddSequence(fsedip->fdbp, 0 ,fsedip->seqid, fsedip->defline,
                       0, 0, 0, bsp->seq_data_type, &bsp->seq_data, 
                       bsp->length, 0);
        /* Reseting mfp structure */
        /* fsedip->tfp->mfp->accession = NULL;
           fsedip->tfp->mfp->organism  = NULL;
           *fsedip->defline = NULLB;
           *fsedip->seqid = NULLB; */

    } else {  /* This will work for example for delta seqs */
        fsedip->seq_data_type = fsedip->tfp->mfp->code;
        fsedip->length = 0;
        fsedip->tfp->mfp->no_sequence = FALSE;
        BSSeek(fsedip->bsp, 0, SEEK_SET);
        SeqEntryFasta(sep, fsedip->tfp, index, indent);
    } 

    return;
}

Boolean FDBAddSeqEntry(FormatDBPtr fdbp, SeqEntryPtr sep)
{
    FDB_SEDataInfoPtr fsedip;
    
    fsedip = FDB_SEDataInfoNew();
    fsedip->fdbp = fdbp;
    
    fsedip->tfp->is_na = !fsedip->fdbp->options->is_protein;
    
    if (fsedip->tfp->is_na){
        fsedip->tfp->mfp->code = Seq_code_iupacna;
    } else {
        fsedip->tfp->mfp->code = Seq_code_ncbistdaa;
    }
    
    SeqEntryExplore(sep, fsedip, FDBSeqEntry_callback);
    
    FDB_SEDataInfoFree(fsedip);
       
    return TRUE;
}


/* ---------------------------------------------------------------------*/
/* ------------- End of functions, that uses in formatdb -------------- */
/* ---------------------------------------------------------------------*/

/************************************************************************/
/*	The CommonIndex stuff						*/
/************************************************************************/

/* The function initializes CommonIndexPtr with give filename */

CommonIndexHeadPtr	CommonIndexInit(CharPtr indexfilename)
{

    Nlm_MemMapPtr	mmpindx;
    CommonIndexHeadPtr	cihp = (CommonIndexHeadPtr) MemNew(sizeof(CommonIndexHead));
    CharPtr		charptr = NULL;

    if (!(mmpindx = Nlm_MemMapInit(indexfilename))) {
	ErrPostEx(SEV_ERROR, 0, 0, "Could not open Common Index file.  Probably wrong path specified\n");
	return NULL;
    }

    cihp->maxgi = FileLength(indexfilename) / sizeof(CommonIndex);
    cihp->memmap = mmpindx;
    cihp->ci = (CommonIndexPtr) mmpindx->mmp_begin;

    /* read list of databases from the configuration file */

    charptr = Nlm_FilePathFind(indexfilename);
    if (!(cihp->num_of_DBs = ParseDBConfigFile(&(cihp->dbids), charptr))) {
#if 0        
	ErrPostEx(SEV_ERROR, 0, 0, "Could not find or parse DB configuration file \n");
#endif        
	return NULL;
    }
    if (charptr)
	MemFree(charptr);
    
    if (!(cihp->ci)) {
#if 0        
	ErrPostEx(SEV_ERROR, 0, 0, "Could not initialize Common Index file %s\n", indexfilename);
#endif        
	return NULL;
    } else
	return cihp;
}
  
void	CommonIndexDestruct(CommonIndexHeadPtr cihp) {

    Int2	i;

    if (cihp &&  cihp->memmap)
	Nlm_MemMapFini(cihp->memmap);

    for (i=0; i < cihp->num_of_DBs; i++) {
	if (cihp && cihp->dbids && ((cihp->dbids + i)->name))
	    MemFree((cihp->dbids + i)->name);
    }
    if (cihp && cihp->dbids)
	MemFree(cihp->dbids);
    
    MemFree(cihp);
}
/* returns shift of bit for specified DB name */

Int2	DBShift(Int2 num_of_DBs, DataBaseIDPtr dbids, CharPtr dbname, Boolean is_prot)
{
    Int2	i;

    if (!dbname) {
	ErrPostEx(SEV_ERROR, 0, 0, "Specified database name is NULL\n");
	return 0;
    }

    for(i=0; i < num_of_DBs; i++) {
	if(!StrCmp(dbname, (dbids+i)->name) && ((dbids+i)->isprot == is_prot)) {
	    return (dbids+i)->id;
	}
    }

#if 0    
    ErrPostEx(SEV_ERROR, 0, 0, "Specified database name %s is not known\n", dbname);
#endif    
    return 0;
}

/* returns name of the database by given bit shift */

CharPtr	DBName(Int2 num_of_DBs, DataBaseIDPtr dbids, Int2 shift)
{
    Int2      i;

    if (!shift) {
	ErrPostEx(SEV_ERROR, 0, 0, "Specified bit shift is zero\n");
	return NULL;
    }

    for(i=0; i < num_of_DBs; i++) {
	if((dbids+i)->id == shift) {
	    return (dbids+i)->name;
	}
    }
    ErrPostEx(SEV_ERROR, 0, 0, "Specified bit shift %d is not known\n", shift);
    return NULL;
}
 
/* say if the database contains proteins */

Boolean	DBisProt(Int2 num_of_DBs, DataBaseIDPtr dbids, Int2 shift)
{
    Int2      i;

    if (!shift) {
	ErrPostEx(SEV_ERROR, 0, 0, "Specified bit shift is zero\n");
	return FALSE;
    }

    for(i=0; i < num_of_DBs; i++) {
	if((dbids+i)->id == shift) {
	    return (dbids+i)->isprot;
	}
    }
    ErrPostEx(SEV_ERROR, 0, 0, "Specified bit shift %d is not known\n", shift);
    return FALSE;
}

void	CommonIndexResultDestruct(CommonIndexResultPtr cir)
{
    if (!cir)
	return;
    if (cir->next)
        CommonIndexResultDestruct(cir->next);
    if (cir)
        MemFree(cir);
}

/* returns OID by given GI */
Int4    GI2OID(CommonIndexHeadPtr cih, Int4 gi, Int4 dbmask, Int4 alias_dbmask,
	Int2Ptr dbid, Int2Ptr alias_dbid, ReadDBFILEPtr rdfp)
{
    CommonIndexResultPtr	cir, cir_start;
    Int4		retval=-1;
    Uint4		dbmask_tmp;
 
    /* gi is not in the database (or even in the common index).
       The most probable reason for this is that the gi was released
       after the database was built.  Return -1 to indicate that it
       is not in the database. */
    if (gi >= cih->maxgi) {
        return -1;
    }
    
    cir_start = GIs2OIDs(cih, &gi, 1, dbmask | alias_dbmask, rdfp);


    /* Get oid in a real database */
    cir = cir_start;
    while (cir && (retval==-1)) {
	if (dbmask & (0x1<<cir->dbid)) {
	    *dbid = cir->dbid;
	    retval = cir->oid;
	}
	cir = cir->next;
    }


    /* now set dbid to correct alias database, if alias database mask specified */
    dbmask_tmp = SwapUint4(cih->ci[gi].dbmask);
    if (dbmask_tmp - dbmask > 0 && alias_dbid)
    	*alias_dbid = bit_engine_firstbit(dbmask_tmp & alias_dbmask);

    CommonIndexResultDestruct(cir_start);

    return retval;
}

/*
   gets list of GI's and returns all OID for each database from the mask
   the GI belongs to.  dbmask == 0 means all databases.
   The list of OID is constructed as list of the CommonIndexResult items
   (see readdb.h for definition)
   noids - number of found oid on return
 */

CommonIndexResultPtr	GIs2OIDs(CommonIndexHeadPtr cih, Int4Ptr gis,
	Int4 number_of_gis, Int4 dbmask, ReadDBFILEPtr startrdfp)
{
    Int4		i, gi, numDB, mask;
    Int2		firstpos, curfirstpos;
    CommonIndexPtr	cigi;
    CommonIndexResultPtr cir = NULL, cirfirst = NULL;
    Boolean		first = TRUE;
    ISAMObjectPtr	nisam_opt = NULL;
    ISAMErrorCode	error;
    Uint4		value;
    ReadDBFILEPtr	rdfp;

    /* for each given GI we need to check if this gi is in list */

    for(i=0; i < number_of_gis; i++) {
	gi = gis[i];
	cigi = cih->ci + gi;

	/* mask says what DBs the GI belongs to */
        mask = SwapUint4(cigi->dbmask);

	if (dbmask && !(dbmask & mask)) {
	    /* skip the gi if it is not in dbmask databases */
	    continue;
	}

	numDB = bit_engine_numofbits(mask);

	if (numDB) {
	    /* Okay, there is at least one database which contains such GI */

	    /* Check if this is the "often" database for the GI */
	    firstpos = bit_engine_firstbit(mask);

	    /* dbmask == 0 means that we search for ALL DBs */
	    if (!dbmask || (dbmask & (0x1 << firstpos))) {
		if (first) {
		    /* create first if needed */
		    cirfirst = (CommonIndexResultPtr) MemNew(sizeof(CommonIndexResult));
		    first = FALSE;
		    cir = cirfirst;
		} else {
		    cir->next = (CommonIndexResultPtr) MemNew(sizeof(CommonIndexResult));
		    cir = cir->next;
		}
		cir->gi = gi;

		/* we know that for the first database the often field is used */
		cir->oid = SwapUint4(cigi->oftenOID);

		cir->dbid = firstpos;
		cir->next = NULL;
	    }
            curfirstpos = firstpos;
            
	    /* do for the rest of databases */
            while (--numDB) {
		/* shift mask to get next database bit shift */
                mask >>= (curfirstpos + 1);
                curfirstpos = bit_engine_firstbit(mask);
		/* update absolute bit shift */
                firstpos += curfirstpos + 1;

		if (!dbmask || (dbmask & (0x1 << firstpos))) {

		    /* find OID using ISAM old index */

		    rdfp = startrdfp;
		    while (rdfp) {
			if (rdfp->filebit == firstpos) {
			    nisam_opt = rdfp->nisam_opt;
			    break;
			}
			rdfp = rdfp->next;
		    }

		    if (!nisam_opt) {
			/* that means that the database specified by 'firstpos' is mask */
			/* skip the database */
			continue;
		    }
		    if (first) {
			cirfirst = (CommonIndexResultPtr) MemNew(sizeof(CommonIndexResult));
			first = FALSE;
			cir = cirfirst;
		    } else {
			cir->next = (CommonIndexResultPtr) MemNew(sizeof(CommonIndexResult));
			cir = cir->next;
		    }

		    cir->gi = gi;

		    /* Initialize and perform the ISAM search */
		    if((error = NISAMSearch(nisam_opt, gi, &value, NULL)) < 0) {
			ErrPostEx(SEV_ERROR, 0, 0, "Failed to initialize ISAM search");
			return NULL;
		    }

		    if(error == ISAMNotFound) {
			ErrPostEx(SEV_ERROR, 0, 0, "Internal error inside GIs2OIDs(), we expected to find this GI into the database\n");
		    }

		    cir->oid = (Int4) value;

		    cir->dbid = firstpos;
		    cir->next = NULL;
		}
            }
	}
    }
    /* return first item of the list */
    return cirfirst;
}


CharPtr	FindDBbyGI(CommonIndexHeadPtr cih, Int4 gi, Uint1 *is_prot)
{
    Int4		numDB, mask;
    Int2		firstpos;
    CommonIndexPtr	cigi;

   if (gi > cih->maxgi)
	return NULL;

    cigi = cih->ci + gi;
    mask = SwapUint4(cigi->dbmask);

    numDB = bit_engine_numofbits(mask);

    if (numDB) {
	firstpos = bit_engine_firstbit(mask);
	*is_prot = DBisProt(cih->num_of_DBs, cih->dbids, firstpos);
	return DBName(cih->num_of_DBs, cih->dbids, firstpos);
    } else {
	return NULL;
    }

}

/* returns senior (first) bit in the word */
Int2	bit_engine_firstbit (Int4 word)
{
    Int2	i;
    Int4	senior_bit = 0x1;

    for (i=0; i < 8*sizeof(Int4); i++) {
	if (word & senior_bit)
	    return i;
	senior_bit <<= 1;
    }
    return -1;
}

/* return number of bits which are ON in the give "word" */
Int2	bit_engine_numofbits(Int4 word)
{
    Int2	i;
    Int4	tmpbit = 0x1;
    Int2	count = 0;

    if (!word) {
	return 0;
    }

    for (i=0; i < 8*sizeof(Int4); i++, tmpbit <<= 1) {
	if (word & tmpbit) {
	    count++;
	}
    }
    return count;
}
/* returns:
   1. list of dbid shifts
   2. number of dbs
 */

Int2Ptr	bit_engine_arr(Int4 word)
{
    Int2	i;
    Int4	tmpbit = 0x1;
    Int2Ptr	retval;
    Int2	count = 0;

    retval = (Int2Ptr) MemNew(sizeof(Int2)*8*sizeof(Int4));

    if (!word) {
	retval[0] = 0;
	return retval;
    }

    for (i=0; i < 8*sizeof(Int4); i++, tmpbit <<= 1) {
	if (word & tmpbit) {
	    retval[count+1] = i;
	    count++;
	}
    }
    retval[0] = count;

    return retval;
}

/* This callback is used in UpdateCommonIndexFile() when DI file is given */
Boolean	DI_updateindex_callback(DI_RecordPtr direc, VoidPtr data)
{
    UpdateIndexStructPtr	uisp = (UpdateIndexStructPtr) data;
    CommonIndex			cigi;
    Int4			dbmask, tmpint;
    Int2Ptr			bitarray;	/* contains positions of DB bits */
    Int2			numDB;
    Int4			retstat;

    cigi.dbmask = 0;
    cigi.oftenOID = 0;

    retstat = fseek(uisp->cifile, direc->gi * sizeof(CommonIndex), SEEK_SET);
    retstat = FileRead(&cigi, sizeof(CommonIndex), 1, uisp->cifile);

    /* get db information for this GI */
    dbmask = SwapUint4(cigi.dbmask); 
    bitarray = bit_engine_arr(dbmask); 
    numDB = bitarray[0];

    if (!numDB) {
	/* the GI does not belong to a database */
	cigi.oftenOID = SwapUint4(direc->oid);
    } else {

	/* already exist one database with this GI */
	/* try to figure out what DB is more often */

	if (uisp->shift <= bitarray[1]) {
	    cigi.oftenOID = SwapUint4(direc->oid);
	}
    }

    /* add mask of given database */
    tmpint = SwapUint4(cigi.dbmask);
    cigi.dbmask = SwapUint4(tmpint | ((0x1 << uisp->shift)));

    fprintf(uisp->fout, "\n%ld \t%ld", direc->gi, direc->oid);

    retstat = fseek(uisp->cifile, direc->gi * sizeof(CommonIndex), SEEK_SET);
    retstat = FileWrite(&cigi, sizeof(CommonIndex), 1, uisp->cifile);
    /* printf("%s", strerror(errno)); */

    return TRUE;
}

Int4	UpdateCommonIndexFile (CharPtr dbfilename, Boolean proteins,
	FILE *fout, CharPtr difile, Int4 gi_threshold)
{
    Char		full_filename[PATH_MAX];
    Int4		num_of_gis;
    CommonIndexPtr	cigi = MemNew(sizeof(CommonIndex));
    FILE		*fildes;
    Nlm_MemMapPtr	mmpisam;
    NISAMKeyDataPtr	kdi;
    Uint4		gi, oid;
    Char		buf[1024], p_or_n;
    Int2		shift;
    Int2Ptr		bitarray;	/* contains positions of DB bits */
    Int2		numDB;
    Int4		dbmask, tmpint;
    Int2		num_of_DBs;
    DataBaseIDPtr	dbidsp = MemNew(sizeof(DataBaseID));


    /* read and parse DB config file */
    if (!(num_of_DBs = ParseDBConfigFile(&dbidsp, getenv("BLASTDB")))) {
	ErrPostEx(SEV_ERROR, 0, 0, "DB configuration file is not present or corrupted");
	return -1;
    }
    
    /* find shift of the database */
    shift = DBShift(num_of_DBs, dbidsp, dbfilename, proteins);
    if (!shift) {
	ErrPostEx(SEV_ERROR, 0, 0, "Could not find db ID for %s", dbfilename);
    }

    fprintf(fout, "The database ID is %d\n", shift);

    /* create full name of the common index file commindex.mm */
    if (getenv("BLASTDB")) {
	sprintf(full_filename, "%s%s%s", getenv("BLASTDB"), DIRDELIMSTR, COMMONINDEX_FN);
    } else {
	sprintf(full_filename, "%s", COMMONINDEX_FN);
    }


    /* open common index file */
    if (!(fildes = FileOpen(full_filename, "rb+"))) {
	ErrPostEx (SEV_ERROR, 0, 0,
		"Cannot open/create file %s\n", full_filename);
	return -1;
    };

    /* To update Common index we use old ISAM files produced by formatdb */
    /* Read ISAM index file */
    if (proteins)
	p_or_n = 'p';
    else
	p_or_n = 'n';

    if (getenv("BLASTDB")) {
	sprintf(buf, "%s/%s.%cnd", getenv("BLASTDB"), dbfilename, p_or_n);
    } else {
	sprintf(buf, "%s.%cnd", dbfilename, p_or_n);
    }

    if (difile) {
	char	subset[4096];
	UpdateIndexStruct	uis;
	uis.shift = shift;
	uis.cifile = fildes;
	uis.fout = fout;

	if (!StringCmp(dbfilename, "month.nt") ||
		!StringCmp(dbfilename, "month.est") ||
		!StringCmp(dbfilename, "month.gss") ||
		!StringCmp(dbfilename, "month.sts") ||
		!StringCmp(dbfilename, "month.pataa") ||
		!StringCmp(dbfilename, "month.patnt") ||
		!StringCmp(dbfilename, "month.htgs")) {
	    sprintf(subset, "month");
	} else {
	    sprintf(subset, "%s", dbfilename);
	}

	ScanDIFile(difile, subset, DI_updateindex_callback, (VoidPtr) &uis, 
			fout, gi_threshold);
	/* Update CommonIndex using .?di file */
    } else {
	/* Update CommonIndex using ISAM files */

	/* create memory map for ISAM index file */
	if ((mmpisam = Nlm_MemMapInit(buf)) == NULL) {
	    ErrPostEx(SEV_ERROR, 0, 0,
		    "\nCould not init mmap for file %s\n", buf);
	    return -1;
	}

	if (!(kdi = (NISAMKeyDataPtr)mmpisam->mmp_begin)) {
	    ErrPostEx(SEV_ERROR, 0, 0, "Could not initialize ISAM index files");
	    return -1;
	}

	/* update common index */

	for (num_of_gis = 0;
		num_of_gis*8 < mmpisam->file_size && kdi->key;
		kdi++, num_of_gis++) {
	    /* get next record in KeyData */
	    gi = SwapUint4(kdi->key);
	    oid = SwapUint4(kdi->data);

	    fseek(fildes, gi * sizeof(CommonIndex), SEEK_SET);
	    FileRead(cigi, sizeof(CommonIndex), 1, fildes);

	    /* get db information for this GI */
	    dbmask = SwapUint4(cigi->dbmask); 
	    bitarray = bit_engine_arr(dbmask); 
	    numDB = bitarray[0];

	    if (!numDB) {
		/* the GI does not belong to a database */
		cigi->oftenOID = SwapUint4(oid);
	    } else {

		/* already exist one database with this GI */
		/* try to figure out what DB is more often */

		if (shift <= bitarray[1]) {
		    cigi->oftenOID = SwapUint4(oid);
		}
	    }

	    /* add mask of given database */
	    tmpint = SwapUint4(cigi->dbmask);
	    cigi->dbmask = SwapUint4(tmpint | ((0x1 << shift)));

	    fprintf(fout, "\n%ld \t%ld", gi, oid);

	    fseek(fildes, gi * sizeof(CommonIndex), SEEK_SET);
	    FileWrite(cigi, sizeof(CommonIndex), 1, fildes);

	}

	Nlm_MemMapFini(mmpisam);
    }

    FileClose(fildes);
    return num_of_gis;
}



/************************************************************************/
/* END	The CommonIndex stuff						*/
/************************************************************************/

/************************************************************************/
/*	The functions used with ID1 dump stuff				*/
/************************************************************************/

Boolean	DB_Subset (Boolean (*f)(DI_Record direc), DI_Record direc)
{
    return (*f)(direc);
}

Boolean	is_EST_HUMAN (DI_Record direc)
{
    return (direc.taxid == 9606);
}
Boolean	is_EST_MOUSE (DI_Record direc)
{
    return (direc.taxid == 10090 ||
	    direc.taxid == 10091 ||
	    direc.taxid == 10092 ||
	    direc.taxid == 35531 ||
	    direc.taxid == 80274 ||
	    direc.taxid == 57486
	    );
}
Boolean	is_EST_OTHERS (DI_Record direc)
{
    return (!is_EST_HUMAN(direc) && !is_EST_MOUSE(direc));
}

Boolean	is_SWISSPROT (DI_Record direc)
{
#if 0
    SeqIdPtr	sip;

    sip = ID1SeqIdForGI(direc.gi);

    if (sip) {
	return (sip->choice == SEQID_SWISSPROT);
    } else {
	printf("\nWARNING: could not get SeqID for GI: %ld\n", direc.gi);
    }
#endif
    return (direc.owner == 6);
}

Boolean	is_MONTH (DI_Record direc)
{
    return (direc.gi > direc.gi_threshold);
}

Boolean	is_PDB (DI_Record direc)
{
    return (direc.owner == 10);
}

/* Function scans .pdi or .ndi file and do callback() for this gi and oid */

Boolean	ScanDIFile(CharPtr difilename, CharPtr subset,
	Boolean (*callback)(DI_RecordPtr direc, VoidPtr data), VoidPtr data,
	FILE *out, Int4 gi_threshold)
{

    FILE		*fdi;
    DI_Record	direc;
    Char		skipstr[128];
    int			readstat, total=0, progress_count=0;
    Boolean		(*criteria)(DI_Record);
    Int4		progress_chunk = 100;

    /* open index file */
    fdi = FileOpen(difilename, "r");

    /* check what subset was wanted */
	if (!strcmp(subset, "est_human")) {
		criteria = &is_EST_HUMAN;
	}	else if (!strcmp(subset, "est_mouse")) {
		criteria = &is_EST_MOUSE;
	}	else if (!strcmp(subset, "est_others")) {
		criteria = &is_EST_OTHERS;
	}	else if (!strcmp(subset, "swissprot")) {
		criteria = &is_SWISSPROT;
	}	else if (!strcmp(subset, "pdb")) {
		criteria = &is_PDB;
	}	else if (!strcmp(subset, "month")) {
		criteria = &is_MONTH;
	}

	/* set gi threshold for month subset */
	direc.gi_threshold = gi_threshold;

    /* each line in index file looks like: */
    /* 2933800 5769963 9606 8 EST 427 -2021038615 38990825 */

    fprintf(out, "\nScan index files");
    while ((readstat = fscanf(fdi, "%ld %ld %ld %ld %s %ld %ld %ld",
		&direc.oid, &direc.gi, &direc.taxid, &direc.owner, 
		skipstr, &direc.len, &direc.hash, &direc.div)) > 0) {

	if (DB_Subset(criteria, direc)) {
	    callback(&direc, data);
	    progress_count++;
	}
	if (!(total % progress_chunk)) {
	    if (progress_count < progress_chunk/3) { 
		printf (".");
	    } else if (progress_count > 2*progress_chunk/3) {
		printf ("X");
	    } else {
		printf ("x");
	    }
	    progress_count = 0;
	    fflush(out);
	}

	total++;
    }

    if (readstat != EOF) {
	fprintf(out, "Error occured");
    }
    FileClose(fdi);

	return TRUE;
}
/************************************************************************/
/* END	The functions used with ID1 dump stuff				*/
/************************************************************************/

/************************************************************************/
/*    	Fastacmd API                           				*/
/************************************************************************/

typedef struct FCMDAccList {
    CharPtr acc;
    Int4 gi;
    struct FCMDAccList *next;
} FCMDAccList, PNTR FCMDAccListPtr;

static FCMDAccListPtr GetAccList(CharPtr file, 
                                 Int4Ptr TotalItems)
{
    Char TmpBuff[128];
    Int4 i, j, k;
    Int4 FileLen = 0;
    FCMDAccListPtr AccList = NULL;
    FCMDAccListPtr AccListTmp, AccListLast;
    Int4 NumNotValid = 0;
    Int4 gi = 0;
  
  if(file == NULL || file[0] == NULLB) {
    *TotalItems = 0;
    return NULL;
  }
  
  FileLen = StringLen(file);
  
  for(i = 0; i < FileLen; i++) {
      
      if(isspace((int)file[i]) || file[i] == ',') /* Rolling spaces */
          continue;
      
      /* This is defence from badly formatted requests */
      
      if(NumNotValid > 10) {
          ErrPostEx(SEV_ERROR, 0, 0, "**** ERROR: Too many invalid Gis/Accessions, "
                 "parsing aborted\n");
          *TotalItems = 0;
          return NULL;
      }
      
      /* Rolling spaces */
      
      j= 0;
      while (j < 128  && i < FileLen) { 
          TmpBuff[j] = file[i];
          j++; i++;
          if(file[i] == ',' || file[i] == '\n')  /* Comma is valid delimiter */
              break;
      }
      TmpBuff[j] = NULLB;
    
      /* Is gi/accession too long ??? */

      if(j == 128) {
          ErrPostEx(SEV_WARNING, 0, 0, "Gi/Accession \"%s\" is too long\r\n", 
                 TmpBuff);
          NumNotValid++;
          
          while(!isspace((int)file[i]) || 
                file[i] == ',' || 
                file[i] == NULLB) /* Rolling until spaces */
              i++;
          continue;  /* Next may be valid ... who knows...?? */   
      }
      
      /* Now validating accession/gi */
      
      for(k =0; k < j; k++) {
          if(!IS_DIGIT(TmpBuff[k])) {
              break;
          }
      }

      if(k == j)
          gi = atol(TmpBuff);
      
      /* If this is valid Accession check and tranfer it to gi */
      
      /* It we come here - we got valid text ID */
      
      if(AccList == NULL) { /* first element */
          AccList = (FCMDAccListPtr) MemNew(sizeof(FCMDAccList));
          AccListTmp = AccList;
          AccListTmp->acc = StringSave(TmpBuff);
          AccListTmp->gi = gi;
          AccListTmp->next = NULL;
          AccListLast=AccListTmp;
          *TotalItems = *TotalItems +1; 
      } else {
          AccListTmp = (FCMDAccListPtr) MemNew(sizeof(FCMDAccList));
          AccListLast->next = AccListTmp;
          AccListTmp->acc = StringSave(TmpBuff);
          AccListTmp->gi = gi;
          AccListTmp->next = NULL;
          AccListLast = AccListTmp;
          *TotalItems = *TotalItems +1;
      }
  }
  if(NumNotValid) {
      ErrPostEx(SEV_ERROR, 0, 0, "**** %d invalid Gi%s/Accession%s present in fastacmd "
             "request\r\n", 
             NumNotValid,
             NumNotValid == 1 ? "" : "s",
             NumNotValid == 1 ? "" : "s"
             );
  }
  return AccList;
}

static void FCMDAccListFree(FCMDAccListPtr falp)
{
    FCMDAccListPtr falp_tmp, falp_next;

    if(falp == NULL)
        return;

    for(falp_tmp = falp; falp_tmp != NULL; falp_tmp=falp_next) {
        falp_next = falp_tmp->next;
        MemFree(falp_tmp->acc);
        MemFree(falp_tmp);
    }
}

Int2 Fastacmd_Search (CharPtr searchstr, CharPtr database,
	CharPtr batchfile, Boolean dupl, Int4 linelen, FILE *out)
{
    BioseqPtr		bsp;
    ReadDBFILEPtr	rdfp;
    Int4		i, fid, TotalItems=0, count;
    FCMDAccListPtr	falp, falp_tmp;
    CharPtr		buffer = NULL;
    FILE		*fd;
    Int4Ptr		ids = NULL;
    Int4		guess_gi = -1;
    Uint1		is_prot = READDB_DB_UNKNOWN;

    if (searchstr)
    	guess_gi = atol(searchstr);

    /* if database name is not specified then use CommonIndex */
    if (!database && guess_gi > 0) {
	CommonIndexHeadPtr	cih = NULL;
	CharPtr			commonindex_file;
	CharPtr			dbname;

	commonindex_file = FindBlastDBFile(COMMONINDEX_FN);

	cih = CommonIndexInit(commonindex_file);
	commonindex_file = MemFree(commonindex_file);

	if (cih) {
	    dbname = FindDBbyGI(cih, guess_gi, &is_prot);

	    if (!dbname) {
		ErrPostEx(SEV_ERROR, 0, 0, "The gi \"%ld\" not found\n", guess_gi);
		return (1);
	    }
	    database = MemNew(1024);
	    StrCpy(database, dbname);
	    CommonIndexDestruct(cih);
	}
    }

    if (!database) {
	database = MemNew(1024);
	StrCpy(database, "nr");
    }

    if((rdfp = readdb_new_ex(database, is_prot, TRUE)) == NULL) {
	ErrPostEx(SEV_ERROR, 0, 0, "ERROR: Cannot initialize readdb engine. Exiting...\n");
	return(1);
    }
    database = MemFree(database);

    if(searchstr != NULL) {
	if((falp =  GetAccList(searchstr, &TotalItems)) == NULL) {
	    ErrPostEx(SEV_ERROR, 0, 0, "ERROR: No valid Gis/Accessions found. Exiting...\n");
	    return(1);
	}
    } else if(batchfile != NULL){
	if((fd = FileOpen(batchfile, "r")) == NULL)
	    return (1);

	buffer = WWWReadFileInMemory(fd, 0, TRUE);

	if((falp =  GetAccList(buffer, &TotalItems)) == NULL) {
	    ErrPostEx(SEV_ERROR, 0, 0, "ERROR: No valid Gis/Accessions found. Exiting...\n");
	    return(1);
	}
    }

    for (falp_tmp = falp; falp_tmp != NULL; falp_tmp = falp_tmp->next) {  

	if(falp_tmp->gi != 0) {
	    fid = readdb_gi2seq(rdfp, falp_tmp->gi, NULL);
	} else {
	    if(!dupl) {
		fid = readdb_acc2fasta(rdfp, falp_tmp->acc);
	    } else {
		count = 0;
		fid = readdb_acc2fastaEx(rdfp, falp_tmp->acc, &ids, &count);
	    }
	}

	if(fid < 0 && fid != -1) { 
	    ErrPostEx(SEV_ERROR, 0, 0, "Accesion search failed for \"%s\" "
		    "with error code %d\n", falp_tmp->acc, fid);
	    return 1;
	} else if (fid == -1) {
	    ErrPostEx(SEV_ERROR, 0, 0, "Entry \"%s\" not found\n", falp_tmp->acc);
	} else if (ids == NULL) { /* gi or SeqId */
	    bsp = readdb_get_bioseq(rdfp, fid);
	    BioseqRawToFastaExtra(bsp, out, linelen);
	    BioseqFree(bsp);  
	} else {
	    for(i = 0; i < count; i++) {
		bsp = readdb_get_bioseq(rdfp, ids[i]);
		BioseqRawToFastaExtra(bsp, out, linelen);
		BioseqFree(bsp);
	    }
	}   
    }
    readdb_destruct(rdfp);
    MemFree(buffer);
    FCMDAccListFree(falp);
    return 0;
} 


/************************************************************************/
/* END	Fastacmd API                           				*/
/************************************************************************/


/*************************************************************************
	This function reads in a list of gi's from a text file
and make a binary gilist file.

The binary gilist format has the following construction:

1.) 1st 4 bytes: a 'magic' number: UINT4_MAX
2.) 2nd 4 bytes: total number of gi's in the file (call this value 'number').
3.) 'number' set of 4 bytes, allowing 4 bytes for each gi.

The function GetGisFromFile first checks what the first 4 bytes
of a file are, if they are the 'magic' number, then it proceeds
to read values assuming a binary format.  If they are not the
'magic' number, then a text format is assumed.

*************************************************************************/

static int LIBCALLBACK
compare_gis(VoidPtr v1, VoidPtr v2)
{
   Uint4 gi1 = *(Uint4Ptr) v1;
   Uint4 gi2 = *(Uint4Ptr) v2;

   return ((gi1<gi2) ? -1 : ((gi1>gi2) ? 1 : 0));
}


#define	GIFILE_LINE_LEN	1024
Int4 LIBCALL
readdb_MakeGiFileBinary (CharPtr input_file, CharPtr output_file)
{
	FILE		*infp=NULL, *outfp=NULL;
  	Int4		index = 0, value, chunk_size = 24, gilist_size;
        Int2		status;
        Char		line[GIFILE_LINE_LEN];
	long		tmplong;
  	Uint4Ptr	gi_list;

	if (!(infp = FileOpen(input_file, "r"))) {
            ErrPostEx(SEV_ERROR, 0, 0, "Unable to open file %s", input_file);
            return -1;
        }
        
	if (!(outfp = FileOpen(output_file, "wb"))) {
            ErrPostEx(SEV_ERROR, 0, 0, "Unable to open file %s", output_file);
            return -1;
        }

	gi_list = MemNew(chunk_size * sizeof(Uint4));

	while (FileGets(line, GIFILE_LINE_LEN, infp)) 
	{
	    /* do correct casting */
	    status = sscanf(line, "%ld", &tmplong);
	    value = tmplong;

	    /* skip non-valid lines */
	    if (status > 0) {
		/* do we have enough space in gi_list ? */
		if (chunk_size < index + 1) {
		    chunk_size *= 2;
		    gi_list = Realloc(gi_list, chunk_size * sizeof(Uint4));
		}

		gi_list[index++] = value;
	    }
	}

	FormatDbUint4Write(READDB_MAGIC_NUMBER, outfp);
	FormatDbUint4Write(index, outfp);
	
	gilist_size = index;
	HeapSort(gi_list, gilist_size, sizeof(Uint4), compare_gis);

	for (index=0; index<gilist_size; index++)
	{
		FormatDbUint4Write(gi_list[index], outfp);
	}

	gi_list = MemFree(gi_list);

	FileClose(infp);
	FileClose(outfp);

	return gilist_size;
}

Int4 FastaToBlastDB(FDB_optionsPtr options, CharPtr basename, 
		    Int4 Bases_In_Volume)
{
   FILE *fd;
   FormatDBPtr fdbp;
   SeqEntryPtr sep;
   BioseqPtr bsp;
   Char filenamebuf[128];
   Int4 count, volume;

   if ((fdbp = FormatDBInit(options)) == NULL)
      return 2;
   if((fd = FileOpen(options->db_file, "r")) == NULL)
      return 3;
   
   /* Get sequences */
   while ((sep = FastaToSeqEntryEx(fd, (Boolean)!options->is_protein, 
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
	       sprintf(filenamebuf, "%s.%02ld", basename, volume);
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
   
   if(FormatDBClose(fdbp))
      return 9;
}
