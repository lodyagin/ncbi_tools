/* ===========================================================================
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
 Collection of SeqAlign Accession utilities.
 Maintainer: Hugues Sicotte
 Authors of the original routines: Hugues Sicotte, Colombe Chappey, Tom Madden, Jinghui Zhang

*/
#include <sequtil.h>
#include <salsap.h>
#include <salpacc.h>
/* 
   Function to return a pointers to a linked list of SeqAligns
   that are in the SeqAligns. The user may NOT free these
   SeqIds, as they are in the SeqAlign data-structure
*/

NLM_EXTERN SeqIdPtr LIBCALL SeqIdPtrFromSeqAlign(SeqAlignPtr salp) {

    SeqIdPtr sip=NULL;
    if(salp==NULL) 
        return NULL;
    if(salp->segtype==1) { /* DenseDiag */
        DenseDiagPtr ddp;
        ddp=(DenseDiagPtr)salp->segs;
        if(ddp!=NULL) 
            sip=ddp->id;
    } else if (salp->segtype==2) { /* DenseSeg */
        DenseSegPtr dsp;
        dsp = (DenseSegPtr) salp->segs;
        if(dsp!=NULL)
            sip=dsp->ids;
    } else if (salp->segtype==3) { /* StdSeg */
        StdSegPtr ssp;
        ssp = (StdSegPtr)salp->segs;
        if(ssp!=NULL) {
            Int4 numids=0;
            sip=ssp->ids;
            while(sip) {/* SeqIds are optional at this level */
                numids++;
                sip=sip->next;
            }
            if(numids==ssp->dim) { /* salp->dim is optional, while ssp->dim is not */
                sip = ssp->ids;
            } else {/* Look in the SeqLocs */
                SeqLocPtr slp;
                SeqIdPtr sip_tmp,sip_last=NULL;
                if (ssp->ids!=NULL) {
                    ssp->ids = SeqIdSetFree(ssp->ids);
                    numids=0;
                }
                slp = ssp->loc;
                while(slp!=NULL) {
                    sip_tmp = SeqIdDup(SeqLocId(slp));
                    if(ssp->ids==NULL) /* "Cache" SeqIds in ids slot */
                        ssp->ids = sip_tmp;
                    else if(sip_last!=NULL)
                        sip_last->next=sip_tmp;
                    if(sip_tmp!=NULL)
                        sip_last=sip_tmp;
                    else {
                        ErrPostEx(SEV_WARNING,0,0,"SeqLoc in StdSeg has NULL SeqLocId()\n");
                    }
                    numids++;
                    slp=slp->next;
                }
                sip = ssp->ids;
                if(numids!=ssp->dim) {
                    ErrPostEx(SEV_WARNING,0,0,"Number of SeqLocs differs from SeqAlign ssp->dim\n");
                    ssp->dim = numids;
                }
            }
        }
    } else if(salp->segtype==4) { /* Packed Seg. Optimal for editing alignments */
        PackSegPtr psp;
        psp = (PackSegPtr) salp->segs;
        if (psp!=NULL)
            sip = psp->ids;
    } else if (salp->segtype ==5) {/* Assume that Linked list of SeqAlign set
                                      always refers to the same set of sequences */
        SeqAlignPtr salp_tmp;
        salp_tmp = (SeqAlignPtr) salp->segs;
        sip = SeqIdPtrFromSeqAlign(salp_tmp);
    } else if (salp->segtype == 6) { /* Temporarely Supported SeqAlign Type */
        CompSegPtr csp;
        csp = (CompSegPtr) salp->segs;
        sip = csp->ids;
    } else {
        ErrPostEx(SEV_WARNING,0,0,"Undefined Seqalign segtype=%d\n",(int)salp->segtype);
    }
    return sip;
}
