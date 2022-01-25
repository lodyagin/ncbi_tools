/* $Id: rpsutil.h,v 6.5 2000/05/02 17:57:19 shavirin Exp $
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
* $Revision: 6.5 $
*
* File Description:
*         Reversed PSI BLAST utilities file
*
* $Log: rpsutil.h,v $
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


#ifdef __cplusplus
}
#endif

#endif /* _RPSUTIL_H_ */
