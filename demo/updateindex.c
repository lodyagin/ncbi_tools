/* $Id: updateindex.c,v 6.5 1999/12/17 21:34:38 egorov Exp $ 
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
* File Name:  $RCSfile: updateindex.c,v $
*
* Author:  Alexey Egorov
*
* Version Creation Date: June 23, 1998
*
* $Revision: 6.5 $
*
* File Description:
*         Update Common Index files
*
* $Log: updateindex.c,v $
* Revision 6.5  1999/12/17 21:34:38  egorov
* Add support for the 'month' subset
*
* Revision 6.4  1999/12/17 20:48:55  egorov
* Fix 'gcc -Wall' warnings and remove old stuff.
*
* Revision 6.3  1999/12/14 19:31:55  egorov
* Add possibility to use DI file for updating CommonIndex
*
* Revision 6.2  1999/11/29 14:45:52  egorov
* Bug fixed.
*
* Revision 6.1  1999/09/24 19:07:49  egorov
* The program updates CommonIndex with new information about a database
*
*
* ==========================================================================
*/

#include <sys/mman.h>
#include <fcntl.h>
#include <readdb.h>
#include <ncbi.h>
#include <math.h>
#include <blast.h>
#include <blastdef.h>
#include <ncbisami.h>

#define NUMARG 4

Args dump_args[NUMARG] = {
    {"Input file for formatting",
	"", NULL,NULL,FALSE,'i',ARG_FILE_IN, 0.0,0,NULL},
    { "Data base contains proteins",
	"T", NULL, NULL, FALSE, 'p', ARG_BOOLEAN, 0.0, 0, NULL},
	{"Full path to DI index file",
	NULL, NULL,NULL,TRUE,'d',ARG_FILE_IN, 0.0,0,NULL},
	{ "GI threshold for 'month' subset", 
	0, NULL, NULL, TRUE, 'g', ARG_INT, 0.0, 0, NULL},
};

Int2 Main(void) 
{
    Uint4		num_of_gis;
    CharPtr		dbfilename, difile;
    Boolean		proteins;
	Int4		gi_threshold;

    if ( !GetArgs ("updateindex", NUMARG, dump_args) ) {
	return -1;
    }

    /* Update Common Index with specified database */
    puts("Updating current index files");

    dbfilename = (const CharPtr) dump_args[0].strvalue;
    proteins = dump_args[1].intvalue;
    difile = dump_args[2].strvalue;
    gi_threshold = dump_args[3].intvalue;

    printf("Database:  %s ", dbfilename);
    if (proteins)
	printf("[prot]\n");
    else
	printf("[nucl]\n");


    num_of_gis = UpdateCommonIndexFile (dbfilename, proteins,
			stdout, difile, gi_threshold);

    printf("\nUpdated %ld GIs", num_of_gis);

    return 0;
}
