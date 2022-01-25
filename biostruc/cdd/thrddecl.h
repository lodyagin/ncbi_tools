/* $Id: thrddecl.h,v 1.5 2000/09/22 22:31:33 hurwitz Exp $
*===========================================================================
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
* File Name:  thrddecl.h
*
* Author:  Stephen Bryant
*
* Initial Version Creation Date: 08/16/2000
*
* $Revision: 1.5 $
*
* File Description: threader
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: thrddecl.h,v $
* Revision 1.5  2000/09/22 22:31:33  hurwitz
* added memory management of ThdTbl (results structure)
*
* Revision 1.4  2000/09/14 22:23:39  hurwitz
* adding memory management for contact-potential, annealing-schedule, and folding-motif
*
* Revision 1.3  2000/08/30 21:33:55  hurwitz
* added new and free functions for Seq_Mtf and Qry_Seq
*
* Revision 1.2  2000/08/23 20:10:42  hurwitz
* added memory management functions for Cor_Def and extern C in .h files
*
* Revision 1.1  2000/08/16 20:45:21  hurwitz
* initial check in of threading routines
*
* ==========================================================================
*/

#if !defined(THRDDECL_H)
#define THRDDECL_H

#ifdef __cplusplus
extern "C" {
#endif

int atd(Fld_Mtf* mtf, Cor_Def* cdf, Qry_Seq* qsq, Rcx_Ptl* pmf,
        Gib_Scd* gsp, Thd_Tbl* ttb, Seq_Mtf* psm, float* trg);

int cprl(Fld_Mtf* mtf, Cor_Def* cdf, Rcx_Ptl* pmf, Cxl_Los** cpr, int ct);

int ttb0(Thd_Tbl* ttb);

int sgoi(int is, int it, Rnd_Smp* pvl, Seg_Ord* sgo);

int slo0(Fld_Mtf* mtf, Cor_Def* cdf, Qry_Seq* qsq, Cur_Loc* sli,
         int cs, int ct, int* mn, int* mx);

int rsmp(Rnd_Smp* pvl);

int sal0(Cor_Def* cdf, Qry_Seq* qsq, Cur_Loc* sli, Cur_Aln* sai,
         int cs, int* mn, int* mx);

int spci(Cor_Def* cdf, Qry_Seq* qsq, Cur_Loc* sli, Cur_Aln* sai,
         int n, Seg_Cmp* spc);

int spni(Cxl_Los** cpl, Cur_Loc* sli, int n, Seg_Nsm* spn);

int cxei(Seg_Nsm* spn, Seg_Cmp* spc, Rcx_Ptl* pmf, Thd_Cxe* cxe);

int cpll(Cor_Def* cdf, Rcx_Ptl* pmf, Qry_Seq* qsq, Cxl_Los** cpr,
         Cur_Aln* sai, Cxl_Los** cpl);

int spel(Cxl_Los** cpl, Cur_Aln* sai, Cur_Loc* sli, int cs, Seg_Gsm* spe,
         Cor_Def* cdf, Seq_Mtf* psm, Seg_Cmp* spc);

int ttbi(Cur_Aln* sai, Cur_Loc* sli, Thd_Gsm* tdg, Thd_Tbl* ttb,
         int nrs, int o);

int cpal(Rcx_Ptl* pmf, Cxl_Los** cpr, Cur_Loc* sli, Cxl_Als** cpa);

int spea(Cor_Def* cdf, Cxl_Als** cpa, Cur_Aln* sai, Cur_Loc* sli,
         int n, int h, Seg_Gsm* spe, Seq_Mtf* psm, Seg_Cmp* spc);

int salr(Cor_Def* cdf, Qry_Seq* qsq, Cur_Loc* sli, Cur_Aln* sai,
         int cs, int* mn, int* mx);

int salu(Cor_Def* cdf, Qry_Seq* qsq, Cur_Loc* sli,
         int cs, int al, Cur_Aln* sai);

int algs(Rnd_Smp* pvl, int tm);

int slor(Fld_Mtf* mtf, Cor_Def* cdf, Qry_Seq* qsq, Cur_Loc* sli, Cur_Aln* sai,
         int cs, int ct, int* mn, int*mx);

int slou(Fld_Mtf* mtf, Cor_Def* cdf, int cs, int ct, int of, Cur_Loc* sli);

float bwfi(Thd_Tbl* ttb, Gib_Scd* gsp, Thd_Tst* tts);

float zsc(Thd_Tbl* ttb, Seq_Mtf* psm, Qry_Seq* qsq, Cxl_Los** cpr,
          Cor_Def* cdf, Rcx_Ptl* pmf, Seg_Gsm* spe, Cur_Aln* sai, Rnd_Smp* pvl);

int g(Seg_Gsm* spe, Seq_Mtf* psm, Thd_Gsm* tdg, Seg_Cmp* spc);

float g0(Seg_Nsm* spn, Thd_Cxe* cxe);

float dgri(Seg_Gsm* spe, Seg_Nsm* spn, Thd_Cxe* cxe, Thd_Gsm* tdg,
           Seq_Mtf* psm, Seg_Cmp* spc);

Cor_Def*  NewCorDef(int NumBlocks);
Cor_Def*  FreeCorDef(Cor_Def* cdf);

Seq_Mtf*  NewSeqMtf(int NumResidues, int AlphabetSize);
Seq_Mtf*  FreeSeqMtf(Seq_Mtf* psm);

Qry_Seq*  NewQrySeq(int NumResidues, int NumBlocks);
Qry_Seq*  FreeQrySeq(Qry_Seq* qsq);

Rcx_Ptl*  NewRcxPtl(int NumResTypes, int NumDistances, int PeptideIndex);
Rcx_Ptl*  FreeRcxPtl(Rcx_Ptl* pmf);

Gib_Scd*  NewGibScd(int NumTempSteps);
Gib_Scd*  FreeGibScd(Gib_Scd* gsp);

Fld_Mtf*  NewFldMtf(int NumResidues, int NumResResContacts, int NumResPepContacts);
Fld_Mtf*  FreeFldMtf(Fld_Mtf* mtf);

Thd_Tbl*  NewThdTbl(int NumResults, int NumCoreElements);
Thd_Tbl*  FreeThdTbl(Thd_Tbl* ttb);

#ifdef __cplusplus
}
#endif

#endif
