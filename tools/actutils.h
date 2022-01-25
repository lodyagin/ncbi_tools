/* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*            National Center for Biotechnology Information (NCBI)
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government do not place any restriction on its use or reproduction.
*  We would, however, appreciate having the NCBI and the author cited in
*  any work or product based on this material.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
* ===========================================================================
*
* File Name:  actutils.h
*
* Author:  Sarah Wheelan
*
* Version Creation Date:   2/00
*
* $Revision: 6.4 $
*
* File Description: utility functions for alignments
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: actutils.h,v $
* Revision 6.4  2000/04/18 13:57:14  wheelan
* added AlnMgrForcePairwiseContinuousEx
*
* Revision 6.3  2000/03/14 11:25:48  wheelan
* added ACT_ProfileFree functions
*
* Revision 6.2  2000/03/02 21:11:06  lewisg
* use bandalign for import sequence, make standalone ddv use viewmgr, make dialogs modal, send color update
*
* Revision 6.1  2000/02/11 17:31:45  kans
* initial checkin of functions depending upon blast/bandalign (SW)
*
* ==========================================================================
*/

#ifndef __ACTUTILS__
#define __ACTUTILS__

#undef NLM_EXTERN
#ifdef NLM_IMPORT
#define NLM_EXTERN NLM_IMPORT
#else
#define NLM_EXTERN extern
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include <ncbi.h>
#include <alignmgr.h>
#include <sqnutils.h>
#include <bandalgn.h>

#define CG_WINDOWSIZE 200

#define MAX_LEN 12000


/* defines for profile builder -- how many options are allowed for
*  protein or nucleotide sequences */
#define ACT_NUCLEN 5
#define ACT_PROTLEN 27

typedef struct act_sitelist {
   Int4     start;
   CharPtr  name;
   struct act_sitelist PNTR next;
} ACT_sitelist, PNTR ACT_sitelistPtr;

typedef struct act_topscore {
   FloatHi  score;
   Int4     pos;
   struct act_topscore PNTR next;
} ACT_TopScore, PNTR ACT_TopScorePtr;

typedef struct act_position {
   Int4Ptr  posarray;
   FloatHi  score;
} ACT_Position, PNTR ACT_PositionPtr;

typedef struct act_profile
{
   FloatHiPtr PNTR freq;
   Boolean    nuc;
   Int4       len;
   Int4       numseq;
   FloatHi    confidence;
   struct act_profile PNTR next;
} ACTProfile, PNTR ACTProfilePtr;

typedef struct act_placebounds {
   ACT_TopScorePtr PNTR ats;
   Int4            len;
   ACTProfilePtr   app;
   ACT_PositionPtr  apos;
   Int4Ptr         currpos;
   Int4            currprof;
   Int4            nprof;
   Int4Ptr         boundarray;
   Int4Ptr         numats;
   ACT_TopScorePtr PNTR currats;
} ACT_PlaceBounds, PNTR ACT_PlaceBoundsPtr;


typedef struct act_cginfo {
   Int4             from;
   Int4             to;
   Int4             length;
   Int4             cg;
   Int4             a;
   Int4             c;
   Int4             t;
   Int4             g;
   Int4             n;
   ACT_sitelistPtr  asp;
   CharPtr          sequence;
   struct act_cginfo PNTR next;
} ACT_CGInfo, PNTR ACT_CGInfoPtr;

typedef struct act_statetable
{
   Int4Ptr PNTR  state;
} ACT_Statetab, PNTR ACT_StatetabPtr;

NLM_EXTERN ACT_CGInfoPtr ACT_FindCpG(BioseqPtr bsp);
NLM_EXTERN Uint1 ACT_GetResidue(Int4 pos, Uint1Ptr buf, Int4Ptr offset, BioseqPtr bsp);
NLM_EXTERN ACTProfilePtr ACT_ProfileNew(Boolean nuc);

/***************************************************************************
*
*  ACT_ProfileFree frees a single profile; ACT_ProfileSetFree frees an
*  entire linked list of profiles.
*
***************************************************************************/
NLM_EXTERN ACTProfilePtr ACT_ProfileFree(ACTProfilePtr app);
NLM_EXTERN ACTProfilePtr ACT_ProfileSetFree(ACTProfilePtr app);

NLM_EXTERN void ACT_BuildProfile(SeqLocPtr slp, ACTProfilePtr app);
NLM_EXTERN FloatHi ACT_ScoreProfile(BioseqPtr bsp, Int4 pos, Uint1 strand, ACTProfilePtr app);
NLM_EXTERN void ACT_EstimateConfidence(ACTProfilePtr app);
NLM_EXTERN ACTProfilePtr ACT_SortProfilesByConfidence(ACTProfilePtr app);
NLM_EXTERN int LIBCALLBACK ACT_CompareProfileConfidence(VoidPtr base, VoidPtr large_son);
NLM_EXTERN ACTProfilePtr ACT_MakeProfileFromSA(SeqAlignPtr sap);
NLM_EXTERN Boolean ACT_AddBioseqToSAByProfile(SeqAlignPtr sap, BioseqPtr bsp);
NLM_EXTERN ACT_TopScorePtr PNTR ACT_SortAndTruncate(ACT_TopScorePtr PNTR ats);
NLM_EXTERN ACT_TopScorePtr ACT_FindPeakScores(FloatHiPtr scorearray, Int4 len);
NLM_EXTERN ACT_PositionPtr ACT_PlaceByScore(ACT_PlaceBoundsPtr abp);
NLM_EXTERN FloatHi ACT_CalcScore(ACT_PlaceBoundsPtr abp);
NLM_EXTERN SeqAlignPtr ACT_GlobalAlignSimple(BioseqPtr bsp1, BioseqPtr bsp2,
                                             Boolean Default);
NLM_EXTERN SeqAlignPtr ACT_GlobalAlignTwoSeq(BioseqPtr bsp1, BioseqPtr bsp2);
NLM_EXTERN SeqAlignPtr AlnMgrForcePairwiseContinuous(SeqAlignPtr sap);
NLM_EXTERN SeqAlignPtr AlnMgrForcePairwiseContinuousEx(SeqAlignPtr sap, Int4 start_1, Int4 stop_1, Int4 start_2, Int4 stop_2);

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