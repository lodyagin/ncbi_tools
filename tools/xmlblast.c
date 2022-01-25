/* $Id: xmlblast.c,v 6.7 2000/10/23 22:12:46 shavirin Exp $ */
/**************************************************************************
*                                                                         *
*                             COPYRIGHT NOTICE                            *
*                                                                         *
* This software/database is categorized as "United States Government      *
* Work" under the terms of the United States Copyright Act.  It was       *
* produced as part of the author's official duties as a Government        *
* employee and thus can not be copyrighted.  This software/database is    *
* freely available to the public for use without a copyright notice.      *
* Restrictions can not be placed on its present or future use.            *
*                                                                         *
* Although all reasonable efforts have been taken to ensure the accuracy  *
* and reliability of the software and data, the National Library of       *
* Medicine (NLM) and the U.S. Government do not and can not warrant the   *
* performance or results that may be obtained by using this software,     *
* data, or derivative works thereof.  The NLM and the U.S. Government     *
* disclaim any and all warranties, expressed or implied, as to the        *
* performance, merchantability or fitness for any particular purpose or   *
* use.                                                                    *
*                                                                         *
* In any work or product derived from this material, proper attribution   *
* of the author(s) as the source of the software or data would be         *
* appreciated.                                                            *
*                                                                         *
************************************************************************** 
* File Name:  xmlblast.c
*
* Author:  Sergei B. Shavirin
*   
* Version Creation Date: 05/17/2000
*
* $Revision: 6.7 $
*
* File Description:  Functions to print simplified BLAST output (XML)
*
* 
* $Log: xmlblast.c,v $
* Revision 6.7  2000/10/23 22:12:46  shavirin
* Added possibility to create XML with error message in case of failure or
* no hits.
*
* Revision 6.6  2000/10/23 19:56:06  dondosha
* AsnIo should be opened and closed outside function BXMLPrintOutput
*
* Revision 6.5  2000/10/12 21:35:31  shavirin
* Added support for OOF alignment.
*
* Revision 6.4  2000/08/22 14:53:00  shavirin
* Changed variable from enthropy to entropy.
*
* Revision 6.3  2000/08/11 16:54:17  kans
* return FALSE instead of NULL for Mac compiler
*
* Revision 6.2  2000/08/10 14:42:33  shavirin
* Added missing comment.
*
* Revision 6.1  2000/08/09 20:40:04  shavirin
* Initial revision.
* *
*
*/

#include <xmlblast.h>

static Int4Ptr PNTR glb_matrix;

Boolean BXMLGetSeqLineForDenseDiag(DenseDiagPtr ddp, HspPtr hsp, Int4 length, Boolean is_aa, Int4Ptr PNTR matrix)
{
    SeqInt si;
    SeqLoc sl;
    Int4 i;
    Uint1 m_res, t_res;
    SeqPortPtr m_spp, t_spp;
    
    if(ddp == NULL || hsp == NULL)
        return FALSE;
    
    hsp->qseq = MemNew(length+1);
    hsp->hseq = MemNew(length+1);
    hsp->midline = MemNew(length+1);

    sl.choice = SEQLOC_INT;
    sl.data.ptrvalue = &si;

    /* !! Here we assume, that this is 2-dimensional DenseDiag where
       first element is query and second element is database sequence 
       This is the case for ALL non-gapped blastp and blastn SeqAligns */
    
    /* SeqLoc for query sequence */
    
    si.id = ddp->id;
    si.from = hsp->query_from;
    si.to = hsp->query_to;
    si.strand = (ddp->strands == NULL) ? 0 : ddp->strands[0];

    
    m_spp = SeqPortNewByLoc(&sl, is_aa ? Seq_code_ncbieaa : Seq_code_iupacna);
    
    /* SeqLoc for the subject */
    
    si.id = ddp->id->next;
    si.from = hsp->hit_from;
    si.to = hsp->hit_to;
    si.strand = (ddp->strands == NULL) ? 0 : ddp->strands[1];
    
    t_spp = SeqPortNewByLoc(&sl, is_aa ? Seq_code_ncbieaa : Seq_code_iupacna);
    
    if(m_spp == NULL || t_spp == NULL) {
        if(m_spp == NULL)
            SeqPortFree(m_spp);
        if(t_spp != NULL)
            SeqPortFree(t_spp);
        return FALSE;
    }
    
    for(i = 0; i < length; ++i) {
        m_res = SeqPortGetResidue(m_spp);
        t_res = SeqPortGetResidue(t_spp);
        
        hsp->qseq[i] = m_res;
        hsp->hseq[i] = t_res;
        
        if(m_res == t_res) {
            hsp->midline[i] = is_aa ? m_res : '|';
        } else if(matrix != NULL && is_aa) {
            if (IS_residue(m_res) && IS_residue(t_res) && 
                matrix[m_res][t_res] >0) {
                hsp->midline[i] = '+';
            } else {
                hsp->midline[i] = ' ';
            }
        } else {
            hsp->midline[i] = ' ';
        }
    }

    SeqPortFree(m_spp);
    SeqPortFree(t_spp);
    
    return TRUE;
}
Boolean BXMLGetSeqLineForDenseSeg(DenseSegPtr dsp, HspPtr hsp, Int4 length,
                                  Boolean is_aa, Int4Ptr PNTR matrix)
{

    SeqInt msi, tsi;
    SeqIntPtr sint;
    SeqLoc sl;
    Int2 i, k;
    Int2 dim;
    Int2 m_order, t_order;	/*order of the master and the target sequence*/
    Int4 index, abs_index;
    Int4 j, val, t_val;
    Uint1 m_res, t_res, stdaa_res;
    SeqIdPtr sip;
    SeqPortPtr m_spp, t_spp;
    SeqMapTablePtr smtp;
    
    
    if(dsp == NULL || hsp == NULL)
        return FALSE;
    
    hsp->qseq = MemNew(length+1);
    hsp->hseq = MemNew(length+1);
    hsp->midline = MemNew(length+1);

    m_order = 0;
    t_order = 1;
    dim = 2;
    
    msi.id = dsp->ids;
    msi.from = hsp->query_from;
    msi.to = hsp->query_to;
    msi.strand = (dsp->strands == NULL) ? 0 : dsp->strands[m_order];
    
    tsi.id = dsp->ids->next;
    tsi.from = hsp->hit_from;
    tsi.to = hsp->hit_to;
    tsi.strand = (dsp->strands == NULL) ? 0 : dsp->strands[t_order];
        
    sl.choice = SEQLOC_INT;
    sl.data.ptrvalue = &msi;
    m_spp = SeqPortNewByLoc(&sl, 
                            (is_aa) ? Seq_code_ncbieaa : Seq_code_iupacna);
    
    sl.choice = SEQLOC_INT;
    sl.data.ptrvalue = &tsi;
    t_spp = SeqPortNewByLoc(&sl, 
                            (is_aa) ? Seq_code_ncbieaa : Seq_code_iupacna);
    
    abs_index = 0;
    for(i = 0; i<dsp->numseg; ++i) {
        val = dsp->starts[i*dim + m_order];
        t_val = dsp->starts[i*dim + t_order];
        if(val == -1 || t_val == -1) {
            
            if(val == -1 && t_val == -1) /* never should happend */
                continue;
            
            if(val != -1) {
                index = dsp->lens[i];
                while (index > 0) {
                    index--;
                    m_res = SeqPortGetResidue(m_spp);
                    hsp->qseq[abs_index] = m_res;
                    hsp->hseq[abs_index] = '-';
                    hsp->midline[abs_index] = ' ';
                    abs_index++;
                }
            }
            if(t_val != -1) {
                index = dsp->lens[i];
                while (index > 0) {
                    index--;
                    t_res = SeqPortGetResidue(t_spp);
                    hsp->qseq[abs_index] = '-';
                    hsp->hseq[abs_index] = t_res;
                    hsp->midline[abs_index] = ' ';
                    abs_index++;
                }
            }
        } else {
            for(j = 0; j<dsp->lens[i]; ++j) {
                m_res = SeqPortGetResidue(m_spp);
                t_res = SeqPortGetResidue(t_spp);

                hsp->qseq[abs_index] = m_res;
                hsp->hseq[abs_index] = t_res;

                if(m_res == t_res) {
                    if(is_aa)
                        hsp->midline[abs_index] = m_res;
                    else
                        hsp->midline[abs_index] = '|';
                    
                } else if(matrix != NULL && is_aa && 
                          IS_residue(m_res) && IS_residue(t_res)) {
                    
                    if(matrix[m_res][t_res] >0)
                        hsp->midline[abs_index] = '+';
                    else
                        hsp->midline[abs_index] = ' ';
                    
                } else {
                    hsp->midline[abs_index] = ' ';
                }
                
                abs_index++;
            }
        }
    }
    
    SeqPortFree(m_spp);
    SeqPortFree(t_spp);
    return TRUE;
}
Boolean BXMLGetSeqLineForStdSeg(StdSegPtr ssp, HspPtr hsp, Int4 length,
                                Int4Ptr PNTR matrix)
{
    Boolean master_is_translated=FALSE, both_translated=FALSE;
    Boolean target_is_translated = FALSE;
    CharPtr genetic_code1, genetic_code2;
    SeqPortPtr spp1, spp2;
    Uint1 codon[4], residue1, residue2;
    Int4 abs_index;
    Boolean  ungapped_align = FALSE;
    
    if(ssp == NULL || hsp == NULL)
        return FALSE;
    
    hsp->qseq = MemNew(length+1);
    hsp->hseq = MemNew(length+1);
    hsp->midline = MemNew(length+1);
    
    /* Check for valid sequence. */
    if (SeqLocLen(ssp->loc) == 3*SeqLocLen(ssp->loc->next))
        master_is_translated = TRUE;
    else if (3*SeqLocLen(ssp->loc) == SeqLocLen(ssp->loc->next))
        target_is_translated = TRUE;	
    else if (SeqLocLen(ssp->loc) == SeqLocLen(ssp->loc->next))
        both_translated = TRUE;
    else
        return FALSE;
    
    if (master_is_translated) {
        genetic_code1 = GetGeneticCodeFromSeqId(ssp->ids);
    } else if (both_translated) {
        genetic_code1 = GetGeneticCodeFromSeqId(ssp->ids);
        genetic_code2 = GetGeneticCodeFromSeqId(ssp->ids->next);
    } else {
        genetic_code1 = GetGeneticCodeFromSeqId(ssp->ids->next);
    }
    
    
    for(abs_index = 0; ssp != 0 && abs_index < length; ssp = ssp->next) {

        if (ssp->loc->choice != SEQLOC_EMPTY || ssp->loc->next->choice != SEQLOC_EMPTY) {
            /* Non - gap element */
            if (ssp->loc->choice != SEQLOC_EMPTY && ssp->loc->next->choice != SEQLOC_EMPTY) {
                if (both_translated) { /* TBLASTX */
                    spp1 = SeqPortNewByLoc(ssp->loc, Seq_code_ncbi4na);
                    spp2 = SeqPortNewByLoc(ssp->loc->next, Seq_code_ncbi4na);
                    while ((codon[0]=SeqPortGetResidue(spp2)) != SEQPORT_EOF) {
                        codon[1] = SeqPortGetResidue(spp2);
                        codon[2] = SeqPortGetResidue(spp2);
                        residue1 = AAForCodon(codon, genetic_code1);
                        codon[0] = SeqPortGetResidue(spp1);
                        codon[1] = SeqPortGetResidue(spp1);
                        codon[2] = SeqPortGetResidue(spp1);
                        residue2 = AAForCodon(codon, genetic_code2);
                        
                        hsp->qseq[abs_index] = residue1;
                        hsp->hseq[abs_index] = residue2;
                        
                        if (residue1 == residue2)
                            hsp->midline[abs_index] = residue1;
                        else if (matrix != NULL && matrix[residue1][residue2] >0)
                            hsp->midline[abs_index] = '+';
                        else
                            hsp->midline[abs_index] = ' ';
                        
                        abs_index++;
                        
                    }
                } else {     /* TBLASTN or BLASTX */
                    if (master_is_translated) { /* BLASTX */
                        spp1 = SeqPortNewByLoc(ssp->loc, Seq_code_ncbi4na);
                        spp2 = SeqPortNewByLoc(ssp->loc->next, Seq_code_ncbieaa);
                    } else {                    /* TBLASTN */
                        spp2 = SeqPortNewByLoc(ssp->loc, Seq_code_ncbieaa);
                        spp1 = SeqPortNewByLoc(ssp->loc->next, Seq_code_ncbi4na);
                    }
                    
                    while ((residue1=SeqPortGetResidue(spp2)) != SEQPORT_EOF) {
                        codon[0] = SeqPortGetResidue(spp1);
                        codon[1] = SeqPortGetResidue(spp1);
                        codon[2] = SeqPortGetResidue(spp1);
                        residue2 = AAForCodon(codon, genetic_code1);
                        
                        if (master_is_translated) { /* BLASTX */
                            hsp->qseq[abs_index] = residue2;
                            hsp->hseq[abs_index] = residue1;
                        } else {                    /* TBLASTN */
                            hsp->qseq[abs_index] = residue1;
                            hsp->hseq[abs_index] = residue2;
                        }
                        
                        if (residue1 == residue2)
                            hsp->midline[abs_index] = residue1;
                        else if (matrix != NULL && matrix[residue1][residue2] >0)
                            hsp->midline[abs_index] = '+';
                        else
                            hsp->midline[abs_index] = ' ';
                        
                        abs_index++;
                    }
                }
                SeqPortFree(spp1);
                SeqPortFree(spp2);
                /* Check if this is an ungapped alignment;
                   in this case do not go to next link */
                if (ssp->next && 
                    ssp->next->loc->choice != SEQLOC_EMPTY && 
                    ssp->next->loc->next->choice != SEQLOC_EMPTY)
                    ungapped_align = TRUE;

            } else if (ssp->loc->choice == SEQLOC_EMPTY) { /* gap in query */

                if (target_is_translated) { /* TBLASTN */
                    spp1 = SeqPortNewByLoc(ssp->loc->next, Seq_code_ncbi4na);
                    
                    while(TRUE) {
                        if((codon[0] = SeqPortGetResidue(spp1)) == SEQPORT_EOF)
                            break;
                        if((codon[1] = SeqPortGetResidue(spp1)) == SEQPORT_EOF)
                            break;
                        if((codon[2] = SeqPortGetResidue(spp1)) == SEQPORT_EOF)
                            break;
                        
                        residue1 = AAForCodon(codon, master_is_translated ? genetic_code2 : genetic_code1);
                        
                        hsp->qseq[abs_index]  = '-'; /* gap character */
                        hsp->hseq[abs_index] = residue1;
                        hsp->midline[abs_index] = ' ';
                        abs_index++;
                    }
                    SeqPortFree(spp1);
                } else {        /* BLASTX */
                    spp1 = SeqPortNewByLoc(ssp->loc->next, Seq_code_ncbieaa);

                    while ((residue1=SeqPortGetResidue(spp1)) != SEQPORT_EOF) {
                        hsp->qseq[abs_index] = '-'; /* gap character */
                        hsp->hseq[abs_index] = residue1;
                        hsp->midline[abs_index] = ' ';
                        abs_index++;
                    }
                    
                    SeqPortFree(spp1);
                }
            } else if (ssp->loc->next->choice == SEQLOC_EMPTY) {
                if (master_is_translated) { /* BLASTX */
                    spp1 = SeqPortNewByLoc(ssp->loc, Seq_code_ncbi4na);

                    while(TRUE) {
                        if((codon[0] = SeqPortGetResidue(spp1)) == SEQPORT_EOF)
                            break;
                        if((codon[1] = SeqPortGetResidue(spp1)) == SEQPORT_EOF)
                            break;
                        if((codon[2] = SeqPortGetResidue(spp1)) == SEQPORT_EOF)
                            break;

                        residue1 = AAForCodon(codon, genetic_code1);
                        
                        hsp->qseq[abs_index] = residue1;
                        hsp->hseq[abs_index] = '-'; /* gap character */
                        hsp->midline[abs_index] = ' ';
                        abs_index++;
                    }
                    SeqPortFree(spp1);
                } else {                    /* TBLASTN */
                    
                    spp1 = SeqPortNewByLoc(ssp->loc, Seq_code_ncbieaa);
                    
                    while ((residue1=SeqPortGetResidue(spp1)) != SEQPORT_EOF) {
                        hsp->qseq[abs_index] = residue1;
                        hsp->hseq[abs_index] = '-'; /* gap character */
                        hsp->midline[abs_index] = ' ';
                        abs_index++;
                    }

                    SeqPortFree(spp1);
                }
            }
            
            if (both_translated || ungapped_align)	
                /* for tblastx perform only one StdSegPtr. so far...*/
                break;
        } else {
            continue;    /* Both EMPTY - never should happened */
        }
    }
    
    return TRUE;
}

Boolean BXMLGetSeqLines(SeqAlignPtr align, HspPtr hsp, Int4 length,
                        Boolean is_aa, Int4 chain, Int4Ptr PNTR matrix)
{
    DenseDiagPtr ddp;
    DenseSegPtr dsp;
    StdSegPtr ssp;
    Uint2 order = 0;
    SeqAlignPtr sap;
    ScorePtr    sp;

    if(align == NULL)
        return FALSE;
    
    switch (align->segtype) {
    case 1: /*Dense-diag*/
        ddp = (DenseDiagPtr) align->segs;

        while(ddp) {
            ++order;
            if(order == chain) {
                BXMLGetSeqLineForDenseDiag(ddp, hsp, length, is_aa, matrix);
                break;
            }
            ddp = ddp->next;
        }

        break;
    case 2:
        dsp = (DenseSegPtr) align->segs;
        BXMLGetSeqLineForDenseSeg(dsp, hsp, length, is_aa, matrix);
        break;
    case 3:
        ssp = (StdSegPtr) align->segs;
        while(ssp) {
            ++order;
            if(order == chain) {
                BXMLGetSeqLineForStdSeg(ssp, hsp, length, matrix);
                break;
            }
            ssp = ssp->next;
        }
        break;
    case 5: /* Discontinuous alignment */
        sap = (SeqAlignPtr) align->segs;
        return BXMLGetSeqLines(sap, hsp, length, is_aa, chain, matrix);
    default:
        break;
    }

    return TRUE;
}

HspPtr BXMLGetHspFromSeqAlign(SeqAlignPtr sap, Boolean is_aa, Int4 chain,
                              Boolean is_ooframe)
{
    HspPtr hsp;
    AlignSum as;
    ScorePtr score, sp;
    Boolean matrix_allocated = FALSE;

    if((hsp = HspNew()) == NULL)
        return NULL;

    MemSet(&as, NULLB, sizeof(AlignSum));
    
    as.master_sip = TxGetQueryIdFromSeqAlign(sap);
    as.target_sip = TxGetSubjectIdFromSeqAlign(sap);

    if(glb_matrix != NULL)
        as.matrix = glb_matrix;
    else {
        as.matrix = load_default_matrix ();
        matrix_allocated = TRUE;
    }
    as.is_aa = is_aa;
    as.ooframe = is_ooframe;

    if(chain == 0)
        chain = 1;
    
    if((score = find_score_in_align(sap, chain, &as)) == NULL) {
        if(matrix_allocated)
            free_default_matrix(as.matrix);
        HspFree(hsp);
        return NULL;
    }

    for(sp = score; sp != NULL; sp = sp->next) {
        if(!(StringICmp(sp->id->str, "e_value")) || 
           !(StringICmp(sp->id->str, "sum_e")))
            hsp->evalue = sp->value.realvalue;
        if(!StringICmp(sp->id->str, "bit_score")) {
            hsp->score = sp->value.realvalue;
        }
    }

    hsp->query_frame = as.m_frame;
    hsp->hit_frame = as.t_frame;
    hsp->identity = as.identical;
    hsp->positive = as.positive + as.identical;
    hsp->gaps = as.gaps;

    hsp->query_from = as.master_from;
    hsp->query_to = as.master_to;
    hsp->hit_from = as.target_from;
    hsp->hit_to = as.target_to;
    
    hsp->density = 0;             /* ???? */

    BXMLGetSeqLines(sap, hsp, as.totlen, is_aa, chain, as.matrix);
    
    hsp->query_from++;
    hsp->query_to++;
    hsp->hit_from++;
    hsp->hit_to++;
    
    if(matrix_allocated)
        free_default_matrix(as.matrix);
    
    return hsp;
}

HitPtr BXMLSeqAlignToHits(SeqAlignPtr seqalign, Boolean ungapped, 
                          Boolean is_ooframe)
{
    HitPtr hitp, hitp_head;
    HspPtr hspp;
    SeqAlignPtr sap, sap2;
    SeqIdPtr subject_sip, new_sip;
    BioseqPtr bsp;
    Char buffer[526];
    HspPtr hsp;
    Int4 hit_count, hsp_count, chain;
    Boolean is_aa;
    
    /* We assume, that  seqalignes in this set ordered by 
       database sequences and by e-values for a given database sequence */

    hitp_head = NULL;
    hit_count = 1;              /* Hits starting with 1 */
    
    for(sap = seqalign; sap != NULL;) {
        
        if(hitp_head == NULL) { /* first element */
            hitp_head = hitp = HitNew();
        } else {
            hitp->next = HitNew();
            hitp = hitp->next;
        }
        hitp->num = hit_count;
        hit_count++;
        
        subject_sip = TxGetSubjectIdFromSeqAlign(sap);
        bsp = BioseqLockById(subject_sip);
        is_aa = (bsp->mol == Seq_mol_aa);

        hitp->def = StringSave(BioseqGetTitle(bsp));

        SeqIdWrite(bsp->id, buffer, PRINTID_FASTA_LONG, sizeof(buffer));
        hitp->id = StringSave(buffer);

        if(bsp->id->choice == SEQID_GI && bsp->id->next != NULL)
            SeqIdWrite(bsp->id->next, buffer, PRINTID_TEXTID_ACCESSION, sizeof(buffer));
        else
            SeqIdWrite(bsp->id, buffer, PRINTID_TEXTID_ACCESSION, sizeof(buffer));
        
        hitp->accession = StringSave(buffer);

        hitp->len = bsp->length;

        hsp_count = 1;          /* Hsps starting with 1 */
        chain = 1;

        for(sap2 = sap; sap2 != NULL;) {
            
            /* Filling info about specific alignments */

            while(TRUE) {
                if(hitp->hsps == NULL) {
                    hspp = hitp->hsps = BXMLGetHspFromSeqAlign(sap2, is_aa, chain, is_ooframe);
                } else {
                    if((hspp->next = BXMLGetHspFromSeqAlign(sap2, is_aa, chain, is_ooframe)) == NULL)
                        break;
                    else
                        hspp = hspp->next;
                }
                
                if(hspp == NULL)
                    break;
                
                if(!ungapped) {  /* Only one chain for gapped */
                    hspp->num = hsp_count;
                    hsp_count++;
                    break;
                }
                chain++;
            }

            sap2 = sap2->next;
            new_sip = TxGetSubjectIdFromSeqAlign(sap2);
            
            if(SeqIdMatch(subject_sip, new_sip)) {
                continue;
            } else {
                sap = sap2;
                break;
            }
        }
    }

    return hitp_head;
}

Boolean BXMLPrintOutput(AsnIoPtr aip, SeqAlignPtr seqalign, 
                        BLAST_OptionsBlkPtr options,CharPtr program,
                        CharPtr database, BioseqPtr query, 
                        ValNodePtr other_returns, Int4 option, CharPtr message)
{
    BlastOutputPtr boutp;
    Char buffer[1024];
    HitPtr hitp;
    TxDfDbInfoPtr dbinfo=NULL;
    ValNodePtr vnp;
    Boolean ungapped = FALSE, is_aa;
    SeqPortPtr spp;
    Int4 i;
    BLAST_KarlinBlkPtr ka_params=NULL, ka_params_gap=NULL;

    if((boutp = BlastOutputNew()) == NULL)
        return FALSE;

    /* For optimization BLOSUM62 may be loaded ones */
    glb_matrix = load_default_matrix ();
    
    SeqIdWrite(query->id, buffer, PRINTID_FASTA_LONG, sizeof(buffer));
    boutp->query_ID = StringSave(buffer);
    boutp->query_def = StringSave(BioseqGetTitle(query));
    boutp->query_len = query->length;

    if(option & BXML_INCLUDE_QUERY) {
        boutp->query_seq = MemNew(query->length+1);
        is_aa = (query->mol == Seq_mol_aa);
        spp = SeqPortNew(query, 0, -1, Seq_strand_plus, 
                         (is_aa) ? Seq_code_ncbieaa : Seq_code_iupacna);
        
        for (i = 0; i < query->length; i++) {
            boutp->query_seq[i] = SeqPortGetResidue(spp);
        }
        spp = SeqPortFree(spp);
    } else {
        boutp->query_seq = NULL;    /* Do we need sequence here??? */
    }
    if(options->gapped_calculation == FALSE || !StringICmp(program, "tblastx"))
        ungapped = TRUE;

    if(seqalign != NULL) {
        boutp->hits = BXMLSeqAlignToHits(seqalign, ungapped, 
                                         options->is_ooframe);
    }

    boutp->program = StringSave(program);
    boutp->db = StringSave(database);
    
    sprintf(buffer, "%s %s [%s]", program, BlastGetVersionNumber(), 
            BlastGetReleaseDate());
    boutp->version = StringSave(buffer);

    /* Reference */

    boutp->reference = BlastGetReference(FALSE);

    /* Filling parameters */
    
    boutp->param = ParametersNew();
    boutp->stat = StatisticsNew();
    
    boutp->param->expect = options->expect_value;
    boutp->param->matrix = StringSave(options->matrix);
    
    boutp->param->sc_match = options->reward;
    boutp->param->sc_mismatch = options->penalty;
    boutp->param->gap_open = options->gap_open;
    boutp->param->gap_extend = options->gap_extend;
    boutp->param->include = options->ethresh;
    
    if(options->filter_string != NULL)
        boutp->param->filter = StringSave(options->filter_string);

    if(other_returns != NULL) {    
        for (dbinfo = NULL, vnp=other_returns; vnp; vnp = vnp->next) {
            switch (vnp->choice) {
            case TXDBINFO:
                dbinfo = vnp->data.ptrvalue;
                break;
            case TXKABLK_GAP:
                ka_params_gap = vnp->data.ptrvalue;
                break;
                
            case EFF_SEARCH_SPACE:
                boutp->stat->eff_space = vnp->data.realvalue;
                break;
            default:
                break;
            }
        }
        
        if(dbinfo != NULL) {
            boutp->stat->db_num= dbinfo->number_seqs;
            boutp->stat->db_len = dbinfo->total_length;
        }
        
        if(ka_params_gap != NULL) {
            boutp->stat->lambda = ka_params_gap->Lambda;
            boutp->stat->kappa = ka_params_gap->K;
            boutp->stat->entropy = ka_params_gap->H;
        }
    
    }
    
    boutp->message = StringSave(message);
    if (aip != NULL)       
        BlastOutputAsnWrite(boutp, aip, NULL);
    
    free_default_matrix(glb_matrix);
    
    BlastOutputFree(boutp);
    
    ObjMgrFreeCache(0);
    
    return TRUE;
}

