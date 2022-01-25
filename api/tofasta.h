/*  tofasta.h
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
* File Name:  tofasta.h
*
* Author:  James Ostell
*   
* Version Creation Date: 7/12/91
*
* $Revision: 6.4 $
*
* File Description:  various sequence objects to fasta output
*
* Modifications:  
* --------------------------------------------------------------------------
* Date	   Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
* $Log: tofasta.h,v $
* Revision 6.4  1998/02/23 16:51:27  egorov
* Changes to make the tofasta.c independent on readdb.h
*
* Revision 6.2  1998/01/27 20:28:10  madden
* Added BioseqRawToFastaExtra with line_length arg
*
* Revision 6.1  1997/10/22 16:44:07  shavirin
* Added definitions for functions: FastaReadSequence() and FastaReadSequenceMem()
*
* Revision 6.0  1997/08/25 18:07:48  madden
* Revision changed to 6.0
*
* Revision 5.6  1997/06/19 18:39:21  vakatov
* [WIN32,MSVC++]  Adopted for the "NCBIOBJ.LIB" DLL'ization
*
* Revision 5.5  1996/10/22 16:00:50  shavirin
* Added new Boolean no_sequence in MyFsa structure to disable
* sequence printing
*
 * Revision 5.4  1996/10/21  21:37:24  shavirin
 * Added definition for function SeqEntrysToDefline()
 *
 * Revision 5.3  1996/10/08  22:27:05  shavirin
 * Moved definition of functions FastaToSeqEntryEx and FastaToSeqBuffEx
 * into include file
 *
 * Revision 5.2  1996/08/15  18:15:23  tatiana
 * CreateDefLine() added
 *
 * Revision 5.1  1996/06/15  17:29:44  ostell
 * changed MyFsa structure by adding do_virtual and tech fields
 * added value of 3 for group_segs
 * addes support of tech to FastaDefLine()
 *
 * Revision 5.0  1996/05/28  13:23:23  ostell
 * Set to revision 5.0
 *
 * Revision 4.1  1996/03/13  19:50:23  shavirin
 * Added definition for new external function FastaToSeqBuff()
 *
 * Revision 4.0  1995/07/26  13:49:01  ostell
 * force revision to 4.0
 *
 * Revision 2.7  1995/05/09  18:43:09  ostell
 * added support for (accession) on GenPept deflines and [organism] on
 * GenPept and PRF deflines
 *
*
* ==========================================================================
*/

#ifndef _NCBI_Tofasta_
#define _NCBI_Tofasta_

#include <seqport.h>

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

#ifdef __cplusplus
extern "C" {
#endif

	                            /* keys returned by FastaWriteFunc */

#define FASTA_ID ((Int2)1)       /* the SeqId string */
#define FASTA_DEFLINE ((Int2)2)	 /* the definition line */
#define FASTA_SEQLINE ((Int2)3)	 /* a line of sequence */
#define FASTA_EOS ((Int2)4)		 /* all the sequence has been printed */
#define FASTA_FORMATDB_AMB ((Int2)5) /* make conversion to standard form for nucleotides */

typedef Boolean (* FastaWriteFunc) PROTO ((BioseqPtr bsp, Int2 key,
                                           CharPtr buf, Uint4 buflen, Pointer data));

typedef struct iteminfo {   /* struct used for defline */
	Uint2	entityID,	   
		    itemID,			   
			itemtype;
} ItemInfo, PNTR ItemInfoPtr;

typedef struct descrinfo {   /* struct used for defline */
	BioseqPtr bsp;
	ValNodePtr vnp;
	ItemInfoPtr iip;
	Uint1 choice;
} DescrInfo, PNTR DescrInfoPtr;

typedef struct myfsa {   /* struct used for fasta searches */
    CharPtr buf;         /* buffer for strings (suggest 255 minimum) */
    Int2 buflen,         /* length of buf */
        seqlen;          /* length of sequence lines (must be buflen-1 or less) */
    Pointer mydata;      /* pointer to your own data */
    FastaWriteFunc myfunc;  /* callback to process parts of fasta file */
    BioseqPtr bsp;       /* Bioseq being processed */
    Boolean bad_asn1; /* set if error in input object like mol not set */
    CharPtr accession;   /* used internally for GenPept def lines */
    CharPtr organism;    /* used internally for GenPept/PRF def lines */
    Uint1 order;          /* used to order def lines for BLAST */
    Boolean do_virtual;   /* if TRUE, instantiate virtual sequences */
    Uint1 tech;           /* for MolInfo.tech */
    Boolean no_sequence;  /* used to disable sequence printing */
    Uint1	code;	/* coding of sequence */
    Boolean	formatdb; /* TRUE, if is used in formatdb */
} MyFsa, PNTR MyFsaPtr;

typedef struct tofasta {
    Boolean is_na;
    Boolean got_one;
    MyFsaPtr mfp;
	Uint1 group_segs;
	Int2 last_indent,
		parts, seg;
} FastaDat, PNTR FastaPtr;
    

NLM_EXTERN Boolean BioseqRawToFastaExtra PROTO((BioseqPtr bsp, FILE *fp, Int2 line_length));
NLM_EXTERN Boolean BioseqRawToFasta PROTO((BioseqPtr bsp, FILE * fp, Boolean is_na));
NLM_EXTERN Boolean SeqEntryToFasta PROTO((SeqEntryPtr sep, FILE * fp, Boolean is_na));
NLM_EXTERN Boolean BioseqToFasta PROTO((BioseqPtr bsp, FILE *fp, Boolean is_na));

NLM_EXTERN void	SeqEntryFasta PROTO ((SeqEntryPtr sep, Pointer data,
                                     Int4 index, Int2 indent));
   

/*****************************************************************************
*
*   SeqEntrysToFasta(sep, fp, is_na, group_segs)
*
*   	group_segs = 0 ... take only raw Bioseqs
*       group_segs = 1 ... group segmented seqs into single entry.. no parts
*       group_segs = 2 ... show only parts of segmented seqs
*       group_segs = 3 ... like 1, but also instantiate virtual sequences
*   
*****************************************************************************/
NLM_EXTERN Boolean SeqEntrysToFasta PROTO((SeqEntryPtr sep, FILE *fp, Boolean is_na, Uint1 group_segs));

NLM_EXTERN Boolean SeqEntrysToFastaX PROTO((SeqEntryPtr sep, MyFsaPtr mfp, Boolean is_na, Uint1 group_segs));
NLM_EXTERN Boolean SeqEntrysToDefline PROTO((SeqEntryPtr sep, 
                                  FILE *fp, Boolean is_na, Uint1 group_segs));
NLM_EXTERN Boolean BioseqRawToFastaX PROTO((BioseqPtr bsp, MyFsaPtr mfp, Boolean is_na));
NLM_EXTERN Boolean BioseqToFastaX PROTO((BioseqPtr bsp, MyFsaPtr mfp, Boolean is_na));

/*****************************************************************************
*
*   FastaFileFunc(key, buf, data)
*   	standard "write to file" callback
*
*****************************************************************************/
NLM_EXTERN Boolean FastaFileFunc PROTO((BioseqPtr bsp, Int2 key,
                                        CharPtr buf, Uint4 buflen, Pointer data));


/*****************************************************************************
*
*   Reads a Fasta File into a SeqEntry structure
*   Conventions:
*   >name def
*   agaggagagagagag
*   agagagagagagagag
*
*   "name" = string is considered SEQID_LOCAL until first white space
*   "def"  = after first white space, and before first newline will be "title"
*   "agaga.." = sequence follows until EOF. can be upper or lower case IUPAC
*
*****************************************************************************/
NLM_EXTERN SeqEntryPtr FastaToSeqEntry PROTO((FILE *fp, Boolean is_na));

/*****************************************************************************
*
*   Reads a Fasta Buffer into a SeqEntry structure
*   Conventions:
*   >name def
*   agaggagagagagag
*   agagagagagagagag
*
*   "name" = string is considered SEQID_LOCAL until first white space
*   "def"  = after first white space, and before first newline will be "title"
*   "agaga.." = sequence follows until EOF. can be upper or lower case IUPAC
*
*****************************************************************************/
NLM_EXTERN SeqEntryPtr FastaToSeqBuff PROTO((CharPtr buffer, 
                                  CharPtr PNTR last_char,Boolean is_na));

NLM_EXTERN SeqEntryPtr FastaToSeqBuffEx(CharPtr buffer, CharPtr PNTR last_char, 
                             Boolean is_na, CharPtr PNTR errormsg,
                             Boolean parseSeqId);
NLM_EXTERN SeqEntryPtr FastaToSeqEntryEx (FILE *fp, Boolean is_na, 
                               CharPtr PNTR errormsg,
                               Boolean parseSeqId); 

/*****************************************************************************
*
*   Boolean FastaReadSequence() - read sequence from file
*
*****************************************************************************/

Boolean FastaReadSequence
(
 FILE *fd,                 /* input file pointer) */
 Boolean is_na,            /* type of sequence */
 Int4Ptr seq_length,       /* Returned length of sequence in residues */
 ByteStorePtr PNTR bs_out, /* Returned pointer to sequence ByteStore */
 CharPtr PNTR errormsg     /* error message for debugging */
 );

/*****************************************************************************
*
*   Boolean FastaReadSequenceMem() - read sequence from buffer
*
*****************************************************************************/

Boolean FastaReadSequenceMem
(
 CharPtr buffer,           /* input buffer with sequence */
 CharPtr PNTR next_char,   /* returned pointer to next FASTA sequence */
 Boolean is_na,            /* type of sequence */
 Int4Ptr seq_length,       /* Returned length of sequence in residues */
 ByteStorePtr PNTR bs_out, /* Returned pointer to sequence ByteStore */
 CharPtr PNTR errormsg     /* error message for debugging */
);

/*****************************************************************************
*
*   This routines get parts needed to make FASTA format from ASN.1
*
*****************************************************************************/
/*****************************************************************************
*
*   FastaId(bsp, buf, buflen)
*      Makes the string for the id part of fasta format.
*      buf should be at least 40 bytes
*
*****************************************************************************/
NLM_EXTERN Boolean FastaId PROTO((BioseqPtr bsp, CharPtr buf, Int2 buflen));

/*****************************************************************************
*
*   FastaDefLine(bsp, buf, buflen, accession, organism)
*   	Finds or makes a FASTA format defline (just locates the string)
*       buf should be very long if possible
*       function truncates if buf not long enough
*       a few deflines are longer than 255
*       if (accession != NULL) prefixes defline with (accession)
*          used for translated GenBank records
*       if (organism != NULL) adds [organism] to end
*       if (tech == MI_TECH_phase1 or phase2, adds order comment to defline)
*
*****************************************************************************/
NLM_EXTERN Boolean FastaDefLine PROTO((BioseqPtr bsp, CharPtr buf, Int2 buflen, CharPtr accession, CharPtr organism, Uint1 tech));

NLM_EXTERN Boolean CreateDefLine PROTO((ItemInfoPtr dip, BioseqPtr bsp, CharPtr buf, Int2 buflen, Uint1 tech, CharPtr accession, CharPtr organism));
/*****************************************************************************
*
*   FastaSeqPort(bsp, is_na, do_virtual)
*   	opens a SeqPort for a fasta output of bsp
*
*****************************************************************************/
NLM_EXTERN SeqPortPtr FastaSeqPort PROTO((BioseqPtr bsp, Boolean is_na,
                                          Boolean do_virtual, Uint1 code));

/*****************************************************************************
*
*   FastaSeqLine(spp, buf, linelen)
*     an open seqport is passed in.
*     fills buf with linelen bases
*     assumes buf[linelen] = '\0'
*     returns FALSE when no more residues to print
*
*****************************************************************************/
NLM_EXTERN Boolean FastaSeqLine PROTO((SeqPortPtr spp, CharPtr buf, Int2 linelen, Boolean is_na));

#ifdef __cplusplus
}
#endif

#undef NLM_EXTERN
#ifdef NLM_EXPORT
#define NLM_EXTERN NLM_EXPORT
#else
#define NLM_EXTERN
#endif

#endif
