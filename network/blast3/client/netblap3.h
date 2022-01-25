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
* File Name: netblap3.h
*
* Author: Tom Madden
*
* Version Creation Date:   05/12/97
*
* File Description: 
*       Application Programming Interface (API) for BLAST network server
*
* RCS Modification History:
* $Log: netblap3.h,v $
* Revision 1.9  1998/09/22 16:14:07  egorov
* Add prototype for parametersToOptions()
*
* Revision 1.8  1998/09/01 20:17:04  madden
*  Fixed uninitialzed problem in BlastNetBioseqFetchDisable, changed prototype
*
* Revision 1.7  1998/05/08 21:41:47  vakatov
* use NetProgressCallback in the function prototypes
*
* Revision 1.6  1998/05/08 20:56:23  madden
* Fix PC compile warnings, rename callback
*
* Revision 1.5  1998/04/23 14:18:42  egorov
* Add number_of_hits parameter to TraditionalBlastReportLoc
*
* Revision 1.4  1998/04/22 18:10:07  egorov
* Add support for SeqLoc to blastcl3
*
* Revision 1.3  1998/04/16 19:35:32  madden
* Added Int4Ptr arg to TraditionalBlastReport specifying the numbers of hits
*
* Revision 1.2  1997/11/25 14:43:36  madden
* Changes to allow iterative searches
*
* Revision 1.1  1997/10/08 19:27:22  madden
* Network support for gapped blast
*
*/

#include <ncbi.h>
#include <sequtil.h>
#include <blstspc.h>
#include <objblst3.h>
#include <ni_types.h>
#include <blastdef.h>
#include <blastpri.h>

typedef Boolean (LIBCALLBACK *NetProgressCallback) PROTO((BlastResponsePtr response, Boolean PNTR cancel));


typedef struct _blastnet3h {
	NI_HandPtr svcp;
} BlastNet3H, PNTR BlastNet3Hptr;

typedef struct _blastnet3_block {
	BioseqPtr bsp;		/* The query. */
	Uint2 prog_type;	/* blast[npx], tblast[nx] */
	CharPtr dbname;		/* Name of database. */
	BlastParametersPtr parameters;
	BlastResponsePtr response;	/* Response from server. */
	SeqLocPtr mask;	/* Sequence to be masked */
	BlastNet3Hptr bl3hptr;	/* BlastNet 3 handler returned from BlastInit. */
	NetProgressCallback callback;
	BLAST_MatrixPtr blast_matrix;	/* Matrix to be sent to server. */
} BlastNet3Block, PNTR BlastNet3BlockPtr;

/*
BlastNet3Hptr PNTR bl3hpp must not be NULL when this is called,
it must be saved and used to call BlastBioseq, as well as BlastFini.
BlastResponsePtr PNTR resp contains a response pointer of type BlastResponse_init.
This function is MT-safe.
*/

Boolean LIBCALL BlastInit PROTO((CharPtr program_name, BlastNet3Hptr PNTR bl3hpp, BlastResponsePtr PNTR resp));

/*
The BlastNet3Hptr returned by BlastInitMt must be passed in.
*/

Boolean LIBCALL BlastFini PROTO((BlastNet3Hptr bl3hptr));

/*
	Function to get a Bioseq for a give SeqIdPtr.  Should be used 
	for BioseqFetch function.
*/

BioseqPtr LIBCALL BlastGetBioseq PROTO((BlastNet3BlockPtr blnet3blkptr, SeqIdPtr sip));

/*
	Function to submit blast queries.  

*/

SeqAlignPtr LIBCALL BlastBioseq PROTO ((BlastNet3BlockPtr blnet3blkptr));

BlastNet3BlockPtr LIBCALL BlastNet3BlockNew PROTO((CharPtr program, CharPtr dbname));

BlastNet3BlockPtr LIBCALL BlastNet3BlockDestruct PROTO((BlastNet3BlockPtr blnet));

CharPtr LIBCALL Blast3GetMotd PROTO((BlastNet3Hptr bl3hptr));

SeqLocPtr LIBCALL BlastGetMaskedLoc PROTO((BlastNet3BlockPtr blnet3blkptr));

Boolean LIBCALL BlastNetBioseqFetchEnable PROTO((BlastNet3Hptr bl3hp, CharPtr dbname, Boolean is_na, Boolean now));

void LIBCALL BlastNetBioseqFetchDisable PROTO((BlastNet3Hptr bl3hp, CharPtr dbname, Boolean is_na));

BlastParametersPtr LIBCALL BlastOptionsToParameters PROTO((BLAST_OptionsBlkPtr options));


CharPtr LIBCALL BlastGetParameterBuffer PROTO((BlastNet3BlockPtr blnet3blkptr));

BlastKABlkPtr LIBCALL BlastGetKaParams PROTO((BlastNet3BlockPtr blnet3blkptr, Boolean gapped));

BlastDbinfoPtr LIBCALL BlastRequestDbInfo PROTO((BlastNet3Hptr bl3hp, CharPtr database, Boolean is_prot));

BlastDbinfoPtr LIBCALL BlastGetDbInfo PROTO((BlastNet3BlockPtr blnet3blkptr));

TxDfDbInfoPtr LIBCALL NetDbinfo2TxDbinfo PROTO((BlastDbinfoPtr net_dbinfo));

Boolean LIBCALL TraditionalBlastReport PROTO((BioseqPtr bsp, BLAST_OptionsBlkPtr options, BlastNet3Hptr bl3hp, CharPtr program, CharPtr database, Boolean html, FILE *outfp, Boolean verbose, Uint4 print_options, Uint4 align_options, Int4 number_of_descriptions, Int4 number_of_alignments, Int4Ptr number_of_hits));

Boolean LIBCALL TraditionalBlastReportLoc PROTO((SeqLocPtr slp, BLAST_OptionsBlkPtr options, BlastNet3Hptr bl3hp, CharPtr program, CharPtr database, Boolean html, FILE *outfp, Boolean verbose, Uint4 print_options, Uint4 align_options, Int4 number_of_descriptions, Int4 number_of_alignments, Int4Ptr number_of_hits));

SeqAlignPtr LIBCALL BlastBioseqNet PROTO((BlastNet3Hptr bl3hp, BioseqPtr bsp, CharPtr program, CharPtr database, BLAST_OptionsBlkPtr options, ValNodePtr *other_returns, ValNodePtr *error_returns, NetProgressCallback callback));

SeqAlignPtr LIBCALL BlastSeqLocNet PROTO((BlastNet3Hptr bl3hp, SeqLocPtr slp, CharPtr program, CharPtr database, BLAST_OptionsBlkPtr options, ValNodePtr *other_returns, ValNodePtr *error_returns, NetProgressCallback callback));

SeqAlignPtr LIBCALL BlastBioseqNetCore PROTO((BlastNet3Hptr bl3hp, BioseqPtr bsp, CharPtr program, CharPtr database, BLAST_OptionsBlkPtr options, ValNodePtr *other_returns, ValNodePtr *error_returns, NetProgressCallback callback, BLAST_MatrixPtr blast_matrix));

SeqAlignPtr LIBCALL BlastSeqLocNetCore PROTO((BlastNet3Hptr bl3hp, SeqLocPtr slp, CharPtr program, CharPtr database, BLAST_OptionsBlkPtr options, ValNodePtr *other_returns, ValNodePtr *error_returns, NetProgressCallback callback, BLAST_MatrixPtr blast_matrix));

BLAST_MatrixPtr LIBCALL BlastNetMatrixToBlastMatrix PROTO((BlastMatrixPtr net_matrix));

BlastMatrixPtr LIBCALL BlastMatrixToBlastNetMatrix PROTO((BLAST_MatrixPtr matrix));

BLAST_OptionsBlkPtr parametersToOptions (BlastParametersPtr parameters, CharPtr program,
	ValNodePtr PNTR error_returns);
