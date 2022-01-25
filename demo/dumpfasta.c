/*  $RCSfile: dumpfasta.c,v $  $Revision: 6.6 $  $Date: 2001/03/26 16:49:43 $
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
*   The program dumps FASTA file for specified database
*
* ---------------------------------------------------------------------------
* $Log: dumpfasta.c,v $
* Revision 6.6  2001/03/26 16:49:43  dondosha
* Allow databases without common index
*
* Revision 6.5  2001/03/05 13:36:09  egorov
* Clean up
*
* Revision 6.4  2001/02/27 21:51:10  madden
* Call BioseqToFastaDump instead of BioseqToFasta
*
* Revision 6.3  2000/03/15 21:35:45  egorov
* Use readdb_get_bioseq_ex()
*
* Revision 6.2  2000/03/02 21:27:11  egorov
* Bug fixed
*
* Revision 6.1  2000/03/01 20:07:23  egorov
* New program added.  It dumps databases in FASTA format.
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

#define NUMARG 4

static Args myargs [NUMARG] = {
    { "Database", 
	NULL, NULL, NULL, FALSE, 's', ARG_STRING, 0.0, 0, NULL},
    { "Database type", 
	FALSE, NULL, NULL, FALSE, 'p', ARG_BOOLEAN, 0.0, 0, NULL},
    { "Output file name", 
	"stdout", NULL, NULL, TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
    { "GI threshold for month subset", 
	0, NULL, NULL, TRUE, 'm', ARG_INT, 0.0, 0, NULL},
};

Int2 Main (void)
{
    Boolean		is_prot;
    CharPtr		database;
    BioseqPtr		bsp;
    CharPtr		filename;
    FILE		*fp;
    ReadDBFILEPtr	rdfp;
    Int4		gi, oid;
    Int4		total = 0, countall = 0, countdumped = 0;
    Int2		progress_chunk = 100;
    Int4		start;
    Uint4		*oidmask;
    Uint4		maskindex, bit;
    Uint4		dumped_already;


    printf("\nRead input parameters");

    if (! GetArgs ("genmask", NUMARG, myargs)) {
	return (1);
    }

    database = myargs [0].strvalue;
    is_prot = myargs [1].intvalue;
    filename = myargs [2].strvalue;

    if (!(fp = FileOpen(filename, "w+"))) {
	ErrPostEx(SEV_ERROR, 0, 0, "Could not open file %s", filename);
    }

    /* creat ReadDBFile */

    if (!(rdfp = readdb_new(database, is_prot))) {
	ErrPostEx(SEV_ERROR, 0, 0, "Could not initialize ReadDBFILEPtr");
    }

    printf("\nStart dumping...\n");
    if (rdfp->cih)
       total = rdfp->cih->maxgi + 1;
    else {
       Int8 tot_len;
       readdb_get_totals(rdfp, &tot_len, &total);
    }

    /* create oid mask to mark those oid which are already being dumped */
    oidmask = (Uint4Ptr) MemNew(total/(sizeof(Uint4)) + 4);

    for (gi=0; gi<total; gi++) {
	if (!(gi%progress_chunk)) {
	    printf("\b\b\b\b%3d%%", (int)((100*gi)/total));
	    fflush(stdout);
	}
        if (rdfp->cih)
           oid = readdb_gi2seq(rdfp, gi, &start);
        else 
           oid = gi;

	if (oid >= 0) {
	    countall++;
	    maskindex = oid/(8*sizeof(Uint4));
	    bit = oid - 8*sizeof(Uint4)*maskindex;
	    dumped_already = oidmask[maskindex] & (0x1<<bit);

	    if (!dumped_already) {
		bsp = readdb_get_bioseq_ex(rdfp, oid, TRUE, TRUE);
		if (!BioseqToFastaDump (bsp, fp, !is_prot)) {
		    ErrPostEx(SEV_ERROR, 0, 0, "Could not convert Bioseq to FASTA");
		}
		BioseqFree(bsp);
		countdumped++;
		/* mark this oid as 'dumped' */
		oidmask[maskindex] |= (0x1<<bit);
	    }
	}
    }

    FileClose(fp);

    printf("\nDatabase: %s, found: %d, dumped: %d.  Total scanned: %d\n",
	    database, countall, countdumped, total);

    return 0;
}
