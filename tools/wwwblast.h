/* $Id: wwwblast.h,v 6.1 2000/05/17 15:52:42 shavirin Exp $
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
* File Name:  $RCSfile: wwwblast.h,v $
*
* Author:  Sergei Shavirin
*
* Initial Creation Date: 03/15/2000
*
* $Revision: 6.1 $
*
* File Description:
*        Definitions for standalone WWW Blast CGI program.
*
* $Log: wwwblast.h,v $
* Revision 6.1  2000/05/17 15:52:42  shavirin
* Initial revision in new location.
*
*
* ==========================================================================
*/

#include <ncbi.h>
#include <blastdef.h>
#include <blast.h>
#include <blastpri.h>
#include <sequtil.h>
#include <txalign.h>
#include <salogif.h>
#include <ddvcreate.h>
#include <objblst3.h>

#ifdef NCBI_CLIENT_SERVER
#include <netblap3.h>
#endif

typedef enum {
    BLASTNoError         =   0,   /* no error             */
    BLASTNetwork         =  -1,
    BLASTNoSpace         =  -2,
    BLASTBadFileName     =  -3,
    BLASTNotImplemented  =  -4,
    BLASTErrProgram      =  -5,  /* program missing from posting data    */
    BLASTErrDatalib      =  -6,  /* datalib missing from posting data    */
    BLASTErrProgName     =  -7,  /* program name is not supported        */
    BLASTErrNoSequence   =  -8,  /* NULL sequence passed to the engine   */
    BLASTErrCombination  =  -9,  /* Invalid program/database combination */
    BLASTNoMemory        =  -10, /* too bad ...                          */
    BLASTNCBI_DATA       =  -11, /* DATA files missing on any path       */
    BLASTFastaToSE       =  -12, /* FastaToSeqEntry() returned NULL      */
    BLASTErrOptions      =  -13, /* Badly formatted options              */
    BLASTErrNoQueue      =  -14, /* Queue overloaded */
    BLASTConfigFile      =  -15, /* Error reading config file */
    BLASTEntrez          =  -16, /* Cannot connect to Entrez */  
    BLASTAccesssion      =  -17, /* Invalid of unavailable accession */  
    BLASTSendmail        =  -18, /* Cannot start sendmail process */ 
    BLASTAddress         =  -19, /* Invalid return address */
    BLASTOptionStr       =  -20, /* Invalidly formatted advanced string */
    BLASTErrAccType      =  -21, /* wrong type of sequence identifier */
    BLASTErrClient       =  -22, /* cannot connect to the Blast service */
    BLASTErrServer       =  -23, /* Error from the server side */
    BLASTMiscError       =  -99  /* undefined internal error             */
} BLASTErrCode;

#define MAX_DB_NUM 256

/* Max. total number of concurrent running requests */
#define DEFAULT_RUN_MAX            2 

/* Max. total number of waiting requests */
#define DEFAULT_QUEUE_MAX          100 
#define NUM_CPU_TO_USE             4
#define DEFAULT_DESCRIPTIONS    100
#define DEFAULT_ALIGNMENTS       50 
#define DEFAULT_EXPECT           10

/* CPU time limit. */
#define DEFAULT_CPU_LIMIT 3600

typedef struct BLASTConfig {
    Int4 run_max;
    Int4 queue_max;
    Int4 num_cpu;
    Int4 niceval;
    CharPtr allow_db[MAX_DB_NUM];
} BLASTConfig, PNTR BLASTConfigPtr;

typedef struct _www_blast_info {
    BLAST_OptionsBlkPtr options;
    WWWInfoPtr info;
    BLASTErrCode error_code;
    CharPtr ConfigFile;
    CharPtr program, database, blast_type;
    BioseqPtr query_bsp;
    BioseqPtr fake_bsp;
    Int4 number_of_descriptions, number_of_alignments;
    Boolean query_is_na, db_is_na, align_type, show_gi, show_overview;
    Boolean believe_query;
    Uint4 align_options, print_options;
    Int4 align_view, input_type, color_schema;
    Boolean is_phi_blast;
    Boolean show_tax_blast;
    BLASTConfigPtr blast_config;
    CharPtr www_root_path;
} WWWBlastInfo, PNTR WWWBlastInfoPtr;


/* Structures used in PHI/PSI Blast searches */

typedef struct GIList {
    Int4  gi;
    struct GIList  *next;
} GIList, PNTR GIListPtr;

typedef struct PSIData {
    CharPtr  matrix62;
    /* BlastSearchBlkPtr search; */
    GIListPtr	PrevCheckedGIs;
    GIListPtr	PrevGoodGIs;
    Int4 StepNumber;
    Nlm_FloatHi karlinK;
} PSIData, PNTR PSIDataPtr;

typedef struct BLASTPrintData {
    SeqAlignPtr        seqalign;
    BLAST_KarlinBlkPtr ka_params, ka_params_gap;
    TxDfDbInfoPtr      dbinfo;
    CharPtr            buffer;  
    ValNodePtr         mask_loc;	  
    PSIDataPtr         psidata;
    SeqLocPtr          seqloc;
    ValNodePtr	       vnp;	/* PHI-BLAST output. */
    ValNodePtr         info_vnp; /* PHI-Blast info strings */
} BLASTPrintData, PNTR BLASTPrintDataPtr;

/* ------------------------------------------- */
void WWWBlastInfoFree(WWWBlastInfoPtr theInfo);

void WWWBlastErrMessage(BLASTErrCode error_code, CharPtr error_msg);

WWWBlastInfoPtr WWWBlastReadArgs(CharPtr type);

Boolean BLAST_Time(CharPtr string, Int4 len, time_t seconds);

Boolean WWWValidateOptions(WWWBlastInfoPtr theInfo);

Boolean WWWCreateSearchOptions(WWWBlastInfoPtr theInfo);

/* PSI/PHI Blast-related function */

void BLASTPrintDataFree(BLASTPrintDataPtr data);

Boolean SplitSeqAlign(SeqAlignPtr seqalign, SeqAlignPtr *GoodSeqAlignment_ptr, SeqAlignPtr *BadSeqAlignment_ptr, SeqAlignPtr *lastGood_ptr, Int2Ptr *marks_ptr, Int2Ptr countBad_ptr, Int2Ptr countGood_ptr, Nlm_FloatHi ethresh_old);

BLASTPrintDataPtr PSIBlastSearch(WWWBlastInfoPtr theInfo);

BLASTPrintDataPtr PHIBlastSearch(WWWBlastInfoPtr theInfo);

