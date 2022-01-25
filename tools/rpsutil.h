/* $Id: rpsutil.h,v 6.12 2000/10/27 15:39:33 kans Exp $
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
* File Name:  $RCSfile: rpsutil.h,v $
*
* Author:  Sergei Shavirin
*
* Initial Version Creation Date: 12/14/1999
*
* $Revision: 6.12 $
*
* File Description:
*         Reversed PSI BLAST utilities file
*
* $Log: rpsutil.h,v $
* Revision 6.12  2000/10/27 15:39:33  kans
* added AnnotateRegionsFromCDD and FreeCDDRegions for common use by ripen, Sequin, and RefSeq processor
*
* Revision 6.11  2000/10/16 19:35:18  shavirin
* Function createFakeProtein() become external.
*
* Revision 6.10  2000/09/28 18:50:11  shavirin
* Added parameter BioseqPtr query_bsp to print results callback.
*
* Revision 6.9  2000/09/27 19:09:41  shavirin
* Significantly redesigned external interface to RPS Blast.
*
* Revision 6.8  2000/09/21 13:49:55  madden
* Rename CddNew and CddDestruct to CddHitNew and CddHitDestruct
*
* Revision 6.7  2000/09/20 22:13:26  madden
* Added code to get CddHit structure
*
* Revision 6.6  2000/09/13 21:26:11  lewisg
* add RPSUpdateDbSize declaration
*
* Revision 6.5  2000/05/02 17:57:19  shavirin
* Corrected path to RPS Databases changed definition of RPSInit() function.
*
* Revision 6.4  2000/03/28 20:32:27  shavirin
* Added functions RPSInfoDetach and RPSInfoAttach.
*
* Revision 6.3  2000/02/25 17:01:14  madden
* Added copyMatrix field
*
* Revision 6.2  2000/02/11 20:44:08  shavirin
* Added possibility to search PSSM database against DNA sequence.
*
* Revision 6.1  1999/12/29 19:39:47  shavirin
* Initial revision.
*
*
* ==========================================================================
*/
#ifndef _RPSUTIL_H_
#define _RPSUTIL_H_

#include <ncbi.h>
#include <blastdef.h>
#include <blastkar.h>
#include <blast.h>
#include <blastpri.h>
#include <readdb.h>
#include <profiles.h>
#include <salpacc.h>

#ifdef __cplusplus
extern "C" {
#endif

/****************************************************************************/
/* DEFINES */
/****************************************************************************/
#define RPS_ALPHABET_SIZE 26    
#define RPS_WORD_SIZE     3
#define RPS_ARRAY_SIZE 32768    /* Size of the lookup table in ModLAEntry */

/****************************************************************************/
/* TYPEDEFS */
/****************************************************************************/

typedef BLAST_Score RPScoreRow[RPS_ALPHABET_SIZE];

typedef struct _RPSLookupHeader {
    Int4 info[8];               /* Information to be determined */
} RPSLookupHeader, PNTR RPSLookupHeaderPtr;

typedef struct _RPSLookup {
    Nlm_MemMapPtr mmLookup;     /* Memory map pointer for the lookup file */
    Int4Ptr  looktbl;           /* Main pointer to the set of lookups */
    Int4Ptr offsets;            /* Offsets to specific lookup tables */
    Int4Ptr header;             /* Lookup header information */
    Int4    entries;            /* Number of entries in the table */
} RPSLookup, PNTR RPSLookupPtr;

typedef struct _RPSequence {   
    UcharPtr sequence;          /* Pointer to the mem-mapped file */
    Int4 number;                /* Number in the database */
    Int4 seqlen;                /* Length of the sequence */
    CharPtr description;        /* Actually this is defline */
    Int4Ptr PNTR posMatrix;     /* PSI Matrix of the sequence */
    Int4Ptr PNTR copyMatrix;     /* copy of Matrix for rescaling. */
    ModLAEntry  *mod_lt;        /* Lookup table for the sequence */
    ModLookupPositionPtr mod_lookup_table_memory; /* Memory tail */
    Int4 num_pos_added;         /* Elements of lookup structure */
    Int4 num_unique_pos_added;  /* Elements of lookup structure */
    Int4 mod_lookup_table_size; /* Elements of lookup structure */
    SeqIdPtr  seqid;            /* Sequence ID in BLAST Database */
} RPSequence, PNTR RPSequencePtr;    

typedef struct _RPSInfo {
    ReadDBFILEPtr rdfp;         /* Handle of the sequence database */
    Nlm_MemMapPtr mmMatrix;     /* Memory map pointer for the matrix file */
    Int4 matrixCount;           /* Total number of PSSM matrixes */
    Int4Ptr offsets;            /* Offsets of matrixes in the file */
    RPScoreRow PNTR bigMatrix;  /* PSI Matrixes for all sequences */
    RPSLookupPtr lookup;        /* Precalculated lookup tables */
    Boolean query_is_prot;      /* Do we need translate query sequence ? */
} RPSInfo, PNTR RPSInfoPtr;
    
/* Aron Bauer's structure used to save CDD hits. */
typedef struct _cdd_hit {
  CharPtr             CDDid;
  CharPtr             ShortName;
  CharPtr             Definition;
  Int4                start;
  Int4                stop;
  Int4                score;
  Nlm_FloatHi         evalue;
  Nlm_FloatHi         bit_score;
  struct _cdd_hit PNTR next;
} CddHit, PNTR CddHitPtr;

typedef struct _rps_blast_options {
    Boolean query_is_protein;
    CharPtr rps_database;
    BLAST_OptionsBlkPtr options;
    Int4 number_of_descriptions, number_of_alignments;
    Boolean html;
    Boolean believe_query;
    Uint4 align_options, print_options;
    Boolean is_xml_output;
    Int4 num_threads;

    /* These parameters are for foprmating convinience only */

    ReadDBFILEPtr rdfp;         /* Handle of the sequence database */
    FILE *outfp;                /* Output file opened descriptor */
    CharPtr out_filename;       /* Output filename */

} RPSBlastOptions, PNTR RPSBlastOptionsPtr;

/* Definitions of multi-threaded batch RPS Blast search */    
typedef SeqEntryPtr (LIBCALLBACK *RPSReadBSPCallback)(SeqLocPtr PNTR slp, 
                                                      VoidPtr user_data);
typedef Boolean   (LIBCALLBACK *RPSHandleResultsCallback)(BioseqPtr bsp, RPSBlastOptionsPtr rpsbop, SeqAlignPtr sap, ValNodePtr other_returns, ValNodePtr error_returns, VoidPtr user_data);
Boolean RPSBlastSearchMT(RPSBlastOptionsPtr rpsbop, 
                         RPSReadBSPCallback bsp_callback, 
                         VoidPtr bsp_user_data,
                         RPSHandleResultsCallback print_callback, 
                         VoidPtr print_user_data);

/****************************************************************************/
/* FINCTION DEFINITIONS */
/****************************************************************************/

/* ----------------------  RPSInit --------------------------
   Purpose:     Initialize main structures of the RPS Search
                
   Parameters:  database - BLAST database of the sequence set corresponding 
                           to PSI matrix set
   Returns:     Poiner to created RPSInfoPtr
  ------------------------------------------------------------------*/
RPSInfoPtr RPSInit(CharPtr database, Int4 query_is_prot);

/* ----------------------  RPSClose --------------------------
   Purpose:     De-initialize main structures of the RPS Search
                
   Parameters:  rpsinfo -  Poiner to created RPSInfoPtr

   Returns:     none
  ------------------------------------------------------------------*/
void RPSClose(RPSInfoPtr rpsinfo);

/* ----------------------  RPSGetSequence --------------------------
   Purpose:     To get sequence information by sequence number
                
   Parameters:  rpsinfo -  Poiner to created RPSInfoPtr
                seqnum - Number of the sequence to be retrieved
   Returns:     Poiner to created RPSequencePtr
  ------------------------------------------------------------------*/
RPSequencePtr RPSGetSequence(RPSInfoPtr rpsinfo, Int4 seqnum);

/* ----------------------  RPSequenceFree --------------------------
   Purpose:     De-initialize RPSequencePtr structure
                
   Parameters:  rpseqp -  Poiner to RPSequencePtr structure

   Returns:     none
  ------------------------------------------------------------------*/
void RPSequenceFree(RPSequencePtr rpseqp);

/* ----------------------  RPSequenceFree --------------------------
   Purpose:     Basic function to do search of the query sequence
                against RPS Blast database

   Parameters:  search - main Blast search structure
                query_bsp - input query sequence Bioseq
                rpsinfo - main object of RPS database

   Returns:     SeqAlign pointer.
  ------------------------------------------------------------------*/
SeqAlignPtr RPSBlastSearch (BlastSearchBlkPtr search,
                            BioseqPtr query_bsp, RPSInfoPtr rpsinfo);


RPSInfoPtr RPSInfoAttach(RPSInfoPtr rpsinfo);
void RPSInfoDetach(RPSInfoPtr rpsinfo);
void RPSUpdateDbSize(BLAST_OptionsBlkPtr options, RPSInfoPtr rpsinfo, Int4 query_length);

BioseqPtr createFakeProtein(void);

/*
	Functions to retrieve CDD hit information.
*/

CddHitPtr CddHitDestruct(CddHitPtr cdd); /* frees memory for every CddHitPtr in linked list. */
CddHitPtr CddHitNew(void); /* produces only one CddHitPtr. */
CddHitPtr RPSBgetCddHits(SeqAlignPtr sap);

/* functions to take BlastBioseqNet result and annotate region features on proteins */

NLM_EXTERN void AnnotateRegionsFromCDD (BioseqPtr bsp, SeqAlignPtr salp, FloatHi expectValue);
NLM_EXTERN void FreeCDDRegions (SeqEntryPtr topsep);


#ifdef __cplusplus
}
#endif

#endif /* _RPSUTIL_H_ */
