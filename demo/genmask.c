/*  $RCSfile: genmask.c,v $  $Revision: 6.4 $  $Date: 1999/12/17 21:34:37 $
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
* Author:  Alexey Egorov
*
* File Description:
*   The program generates mask files for subsets of databases
*
* ---------------------------------------------------------------------------
* $Log: genmask.c,v $
* Revision 6.4  1999/12/17 21:34:37  egorov
* Add support for the 'month' subset
*
* Revision 6.3  1999/12/17 20:48:55  egorov
* Fix 'gcc -Wall' warnings and remove old stuff.
*
* Revision 6.2  1999/12/14 19:33:32  egorov
* Move body of the main() to ScanDIFile() function to readdb.c,
* and use this function with proper callback here, in genmask.
* Few other enhancments and bug fixes.
*
* Revision 6.1  1999/09/24 19:07:14  egorov
* Program generates mask file for virtual subsets of databases
*
* ===========================================================================
*/

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

#include <ncbi.h>
#include <objseq.h>
#include <objsset.h>
#include <sequtil.h>
#include <seqport.h>
#include <tofasta.h>
#include <blast.h>
#include <blastpri.h>
#include <simutil.h>
#include <txalign.h>
#include <gapxdrop.h>
#include <sqnutils.h>
#include <readdb.h>

#define NUMARG 6

static Args myargs [NUMARG] = {
    { "DI file name",
	NULL, NULL, NULL, FALSE, 'i', ARG_STRING, 0.0, 0, NULL},
    { "output alias file", 
	NULL, NULL, NULL, FALSE, 'a', ARG_STRING, 0.0, 0, NULL},
    { "output mask file", 
	NULL, NULL, NULL, FALSE, 'm', ARG_STRING, 0.0, 0, NULL},
    { "output subset", 
	NULL, NULL, NULL, FALSE, 'd', ARG_STRING, 0.0, 0, NULL},
    { "output subset type", 
	FALSE, NULL, NULL, FALSE, 'p', ARG_BOOLEAN, 0.0, 0, NULL},
    { "GI threshold for 'month' subset", 
	0, NULL, NULL, TRUE, 'g', ARG_INT, 0.0, 0, NULL},
};

/* maximum number of sequences in a database */
#define	MAX_DB_SZ	9000000

typedef	struct gms_struct {
    Int4	*memmask;
    Int4	count, dblen, maxoid;
} GenMaskStruct, *GenMaskStructPtr;

Boolean	DI_genmask_callback(DI_RecordPtr direc, VoidPtr data)
{
    GenMaskStructPtr	gmsp = (GenMaskStructPtr) data;

    gmsp->memmask[(int)direc->oid/MASK_WORD_SIZE] |= 
	0x1 << (MASK_WORD_SIZE -1 - direc->oid % MASK_WORD_SIZE);
    gmsp->count++;
    gmsp->dblen += direc->len;

    gmsp->maxoid = MAX (direc->oid, gmsp->maxoid); 

    return TRUE;
}

Int2 Main (void)
{
    CharPtr		db_file, mask_file, alias_file;
    FILE		*fmask, *falias;
    Int4		*memmask = MemNew(MAX_DB_SZ/8 + 1);
    Boolean		is_prot;
    CharPtr		subset;
    Char		TmpBuf[1024];
    GenMaskStruct	gms;
	Int4		gi_threshold;

    printf("\nRead input parameters");

    if (! GetArgs ("genmask", NUMARG, myargs)) {
	return (1);
    }

    db_file = myargs [0].strvalue;
    alias_file = myargs [1].strvalue;
    mask_file = myargs [2].strvalue;
    subset = myargs [3].strvalue;
    is_prot = myargs [4].intvalue;
    gi_threshold = myargs [5].intvalue;

    if (getenv("BLASTDB")) {
	strcpy(TmpBuf, getenv("BLASTDB"));
	strcat(TmpBuf, "/");
	strcat(TmpBuf, db_file);
    } else {
	strcpy(TmpBuf, db_file);
    }

    if (is_prot) {
	strcat(TmpBuf, ".pdi");
    } else {
	strcat(TmpBuf, ".ndi");
    }

    /* initialize */
    gms.memmask = memmask;
    gms.count   = 0;
    gms.dblen   = 0;
    gms.maxoid  = 0;

    /* scan di file with our callback */
    ScanDIFile(TmpBuf, subset, DI_genmask_callback,
			(VoidPtr) &gms, (stdout), gi_threshold);

    /* Generate Alias file */
    printf("\nGenerate alias file");
    falias = FileOpen(alias_file, "w");
    fprintf(falias, "# This is program generated alias file\n#\n");
    fprintf(falias, "TITLE    This is a subset of %s database for %s\n", db_file, subset);
    fprintf(falias, "DBLIST   %s\n", db_file);
    fprintf(falias, "OIDLIST  %s\n", mask_file);
    fprintf(falias, "LENGTH   %ld\n", gms.dblen);
    fprintf(falias, "NSEQ     %ld\n", gms.count);
    fprintf(falias, "MAXOID   %ld\n", gms.maxoid);
    fprintf(falias, "# end of the file");
    FileClose(falias);

    /* Generate Mask file */
    printf("\nWrite mask file to disk");
    fmask = FileOpen(mask_file, "w");
    FileWrite (&gms.maxoid, sizeof(Int4), 1, fmask);
    FileWrite (memmask, gms.maxoid/8 + 1, 1, fmask);
    FileClose(fmask);
    printf("\n%ld sequences found\nTotal length is %ld\n", gms.count, gms.dblen);
    printf("\nDone\n");
	return 0;
}
