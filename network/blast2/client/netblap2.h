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
*
* File Name: netblap2.h
*
* Author:  Jonathan Epstein, Tom Madden
*
* Version Creation Date:   06/16/95
*
* $Revision: 6.1 $
*
* File Description: 
*       Application Programming Interface (API) for BLAST network server
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
* $Log: netblap2.h,v $
* Revision 6.1  1998/02/19 21:13:49  shavirin
* Added parameter to MT save function BlastBioseqMT()
*
* Revision 6.0  1997/08/25 18:34:17  madden
* Revision changed to 6.0
*
* Revision 5.2  1997/04/11 15:18:38  madden
* renamed callbackWithMon Blast2callbackWithMon, made non-static.
*
 * Revision 5.1  1996/07/01  21:20:00  shavirin
 * Added new function BlastSeqLocMon() allows to user set callback.
 *
 * Revision 5.0  1996/05/28  14:09:11  ostell
 * Set to revision 5.0
 *
 * Revision 4.16  1996/05/10  20:44:49  madden
 * MT safe function calls added.
 *
 * Revision 4.15  1996/03/22  17:28:49  madden
 * Removed formatting prototypes to blast2.h
 *
 * Revision 4.14  1995/11/07  17:55:24  madden
 * removed prototype for PrintTraditionalBlastPreface.
 *
 * Revision 4.13  1995/11/01  15:41:25  madden
 * removed prototypes for TraditionalBlastReportSetUp and TraditionalBlastReportCleanUp
 *
 * Revision 4.12  1995/10/27  21:15:57  madden
 * Removed blast_report_struct (moved to blast2.h).
 *
 * Revision 4.11  1995/10/17  15:03:12  madden
 * added prototype for CheckIfBlastJobCancelled.
 *
 * Revision 4.10  1995/09/20  20:57:29  madden
 * Comments added that describe how to call BlastBioseq.
 *
 * Revision 4.9  1995/09/01  16:43:46  madden
 * changes to allow cancellation of jobs through "Cancel" button of prog. bar.
 *
 * Revision 4.8  1995/08/24  17:23:29  madden
 * Added "get_gi" element to blast_report_struct.
 *
 * Revision 4.7  1995/08/23  13:32:14  madden
 * Added num_of_defline and num_of_align to blast_report_struct.
 *
 * Revision 4.6  1995/08/02  16:20:00  madden
 * Added "qoffset" to BlastReportStruct.
 *
 * Revision 4.5  1995/08/01  21:30:45  madden
 * removed define for BLASTCLI_MATRIX_SIZE.
 *
 * Revision 4.4  1995/08/01  20:44:33  madden
 * Added typdef for BlastReportStruct, changed "traditional" function
 * prototypes to include BlastReportStruct.
 *
 * Revision 4.3  1995/07/28  16:16:37  madden
 * Added prototypes for formatting individual parts of the traditional
 * BLAST report.
 *
 * Revision 4.2  1995/07/28  12:29:36  madden
 * Added prototype for BlastAnIUPACString.
 *
 * Revision 4.1  1995/07/26  22:29:03  kans
 * added prototype for print_usage
 *
 * Revision 4.0  1995/07/26  13:55:34  ostell
 * force revision to 4.0
 *
 * Revision 1.9  1995/07/24  16:34:49  madden
 * removed define for NETBLAP1_BUFFER_SIZE.
 *
 * Revision 1.8  1995/07/18  20:05:36  madden
 * Added "start" (for start of SeqLoc) to PrepRequestInfoPtr.
 *
 * Revision 1.7  1995/07/12  17:44:32  madden
 * Changed prototypes to allow SeqLocPtr for masking of sequence.
 *
 * Revision 1.6  1995/06/28  18:28:07  madden
 * Changed BlastBioseq2 and BlastSeqLoc2 to take an additional argument:
 * a "BLAST0ResponsePtr PNTR".
 *
 * Revision 1.5  1995/06/23  21:39:51  madden
 * Add defines BLAST_SERVER_OMIT_MATRIX, BLAST_SERVER_OMIT_QUERY etc.
 *
 * Revision 1.4  1995/06/22  16:08:08  madden
 * Added prototypes for SubmitResultBlastRequest and SubmitSeqAlignBlastRequest.
 *
 * Revision 1.3  1995/06/22  13:19:01  madden
 * Added five Booleans to PrepRequestInfoPtr that control the amount of
 * output, these correspond to the five Uint1's on BLAST0SearchPtr that
 * determine whether a matrix, query seq, db seq, etc is returned.
 *
 * Revision 1.2  1995/06/22  12:54:45  madden
 * Changed HitDataPtr to BLAST0ResultPtr.
 *
 * Revision 1.1  1995/06/16  11:27:08  epstein
 * Initial revision
 *
 * Revision 1.15  95/05/17  17:59:30  epstein
 * add RCS log revision history
 * 
*/

#ifndef _NETBLAP2_
#define _NETBLAP2_

#include <objblst2.h>
#include <objseq.h>
#include <seqport.h>
#include <prtutil.h>
#include <blast2.h>
#include <ni_types.h>

#define BLAST_SERVER_OMIT_MATRIX 1
#define BLAST_SERVER_OMIT_QUERY 2
#define BLAST_SERVER_OMIT_QUERY_SEQ_IN_SEG 4
#define BLAST_SERVER_OMIT_DB_SEQ_IN_SEG 8

/*------------the structure for the function PrepareRequest (netblap2.c)---*/

typedef struct preprequestinfo {
        Boolean is_na;	/* Is this a nucleic acid? */
	CharPtr defline, /* The fasta definition line. */
		options, /* BLAST command line options */
		dbname, /* name of the database */
		textid, /* textid of the entry */
		progname; /* BLAST program name (e.g., blastn, blastp...) */
	Int4 	gi,	/* gi number */
		length, /* number of residues in entry */
		start;	/* BLAST run starts at this residue */ 
	SeqPortPtr spp; /* SeqPort to get sequence data from. */
/* These Booleans determine the amount of output */
	Boolean return_matrix,		/* Should matrix be returned? */
		return_query, 		/* Should query be sent back? */
		return_BLAST0result,	/* Should Blast0ResultPtr be returned;
					if not, then a Seq-align is returned */
		return_query_seq_in_seg, /* Should the query sequence for a hit
					be returned? */ 
		return_db_seq_in_seg;	/* Should the db sequence for a hit be
					returned? */
	SeqLocPtr mask;	/* which positions should be masked */
} PrepRequestInfo, PNTR PrepRequestInfoPtr;


/* the following are for backwards-compatability */
#define HitDataPtr BLAST0ResultPtr
#define HitDataFree BLAST0ResultFree

typedef Boolean (LIBCALLBACK *ProgressCallback) PROTO((BLAST0ResponsePtr, Boolean PNTR cancel));

/******************************************************************
 *             Multi-Thread safe Blast handle
 *****************************************************************/ 

typedef struct BlastMTHandle {
  NI_HandPtr svcp;
  AsnIoPtr asnin;
  AsnIoPtr asnout;
  Int4 socket;
  CharPtr error;
} BlastMTHandle, PNTR BlastMTHandlePtr;

/******************************************************************
*
*	These initialize and close the connection to the BLAST server.
*	BlastInit should be called at the beginning of a session,
*	BlastFini at the end.
*****************************************************************/ 
Boolean LIBCALL BlastInit PROTO((CharPtr clientId, Boolean ignoreErrs));
Boolean LIBCALL BlastFini PROTO((void));

Boolean LIBCALL BlastInitMT PROTO((BlastMTHandlePtr BlastMThp, 
                                   CharPtr clientId, 
                                   Boolean ignoreErrs));

Boolean LIBCALL BlastFiniMT PROTO((BlastMTHandlePtr BlastMThp));

/******************************************************************
*
*	BlastBioseq submit the sequence data in a Bioseq ("bsp") to the
*	server.  The arguments are:
*		bsp: BioseqPtr containing sequence data,
*		progname: CharPtr with name of program (blastn, blastp, blastx,
			tblastn, tblastx),
		dbname: CharPtr with name of the database (most requests
			use "nr"),
		cmdLineOptions: CharPtr contains "expert" BLAST options that
			control output, score cutoffs, etc.
		blrespPtr: BLAST0ResponsePtr that is filled in by the server.
			This contains the matrix, progress messages, 
			error messages.
		mask_seqloc: SeqLocPtr specifying sequence to be masked (with
			X's for proteins, N's for nucleotides).
		output: Uint4 that specifies amount of output, this corresponds
			to four of the five Uint1's in BLAST0Search:

	1st bit, if set omit matrix (return_matrix is FALSE);
	2nd bit, if set omit query (return_query is FALSE);
	3rd bit, if set omit query seq (return_query_seq_in_seg is FALSE);
	4th bit, if set omit db seq (return_db_seq_in_seg is FALSE);

		ProgressCallback: function that produces progress message.

**********************************************************************/
BLAST0ResultPtr LIBCALL BlastBioseq PROTO((BioseqPtr bsp, CharPtr progname, CharPtr dbname, CharPtr cmdLineOptions, BLAST0ResponsePtr PNTR blrespPtr, SeqLocPtr mask_seqloc, Uint4 output, ProgressCallback progCallback));

BLAST0ResultPtr LIBCALL BlastBioseqMT PROTO((BioseqPtr bsp, CharPtr progname, CharPtr dbname, CharPtr cmdLineOptions, BLAST0ResponsePtr PNTR blrespPtr, SeqLocPtr mask_seqloc, Uint4 output, ProgressCallback progCallback, BlastMTHandlePtr BlastMThp));

BLAST0ResultPtr LIBCALL SimpleBlastBioseq PROTO((BioseqPtr bsp, CharPtr progname, CharPtr dbname, CharPtr cmdLineOptions, Boolean useMonitors));
SeqAnnotPtr LIBCALL BlastBioseq2 PROTO ((BioseqPtr bsp, CharPtr progname, CharPtr dbname, CharPtr blast_params, BLAST0ResponsePtr PNTR blrespPtr, SeqLocPtr mask_seqloc, Boolean useMonitors));
SeqAnnotPtr LIBCALL BlastSeqLoc2 PROTO ((SeqLocPtr slp, CharPtr progname, CharPtr dbname, CharPtr blast_params, BLAST0ResponsePtr PNTR blrespPtr, SeqLocPtr mask_seqloc, Boolean useMonitors));
SeqAnnotPtr LIBCALL BlastSeqLocMon PROTO ((SeqLocPtr slp, CharPtr progname, CharPtr dbname, CharPtr blast_params, BLAST0ResponsePtr PNTR blrespPtr, SeqLocPtr mask_seqloc, ProgressCallback userCallback));

SeqAnnotPtr LIBCALL HitDataToSeqAnnot PROTO((BLAST0ResultPtr, SeqIdPtr));
SeqAnnotPtr LIBCALL HitDataToSeqAnnotAlignment PROTO((BLAST0ResultPtr, SeqIdPtr));

CharPtr FormatResultWithTemplate PROTO ((BLAST0ResultPtr blresp, StdPrintOptionsPtr Spop));

/**************************************************************************
 *        This defines MT Safe function SubmitInfoRequestMT() and 
 *              not MT - safe SubmitInfoRequest() analog
 *************************************************************************/

BLAST0ResponsePtr SubmitInfoRequest PROTO ((BLAST0RequestPtr blreqp));
BLAST0ResponsePtr SubmitInfoRequestMT PROTO ((BLAST0RequestPtr blreqp, 
                                            BlastMTHandlePtr BlastMThp));

BLAST0ResultPtr SubmitResultBlastRequest PROTO ((BLAST0RequestPtr blreqp, BLAST0ResponsePtr PNTR blrespPtr, ProgressCallback progCallback));

SeqAlignPtr SubmitSeqAlignBlastRequest PROTO ((BLAST0RequestPtr blreqp, BLAST0ResponsePtr PNTR blrespPtr, ProgressCallback progCallback));

int print_usage(FILE *fp, ValNodePtr vnp);

BLAST0ResultPtr LIBCALL BlastAnIUPACString PROTO((CharPtr sequence, Int4 length, CharPtr id, CharPtr defline, CharPtr progname, CharPtr dbname, CharPtr cmdLineOptions, BLAST0ResponsePtr PNTR blrespPtr, SeqLocPtr mask_seqloc, Uint4 output, ProgressCallback progCallback));


/*
        CheckIfBlastJobCancelled returns the state of the Boolean 
	"job_cancelled".  If it is "TRUE", then the BLAST job has 
	been cancelled through the monitor (e.g., callbackWithMon).  
	This function need only be called if a "suspicious" return 
	value exists (i.e., NULL was returned by the API function)
	and the client wishes to distinguish a "cancellation" from
	an interruption of service or other error.    
	
*/
Boolean LIBCALL  CheckIfBlastJobCancelled PROTO((void));

Boolean LIBCALLBACK Blast2callbackWithMon PROTO((BLAST0ResponsePtr brp, Boolean PNTR cancel));


#endif /* _NETBLAP2_ */
