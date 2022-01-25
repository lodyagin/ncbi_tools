#include <sequtil.h> /* SeqIdDupList */
#include <salpedit.h>
#include <salptool.h>
#include <salpacc.h>
#include <blast.h>
#include <simutil.h>
#include <alignval.h>
#include <blastpri.h>
#include <txalign.h>
#include <sqnutils.h>
#include <satutil.h>

#define MIN_SEG_SIZE 10 /* used by MergeTwoDspBySIM4 , check_align_match_diagnol<-ModFilterDenseSegAlign<-ModifyAlignList<-SeqAlignConsistentDiagFilter */


/* Blast can now (2/1999) handle longer queries more efficiently.. though
 there is still a hardcoded limit to the number of HSP's */

#define SPLIT_LEN       10000  /* used in compute_alignment */
#define OVERLAP_LEN     500 /* used in compute_alignment */

#define PERCENT_A_RES	90 /* used by get_polyA_index<-check_polyA_tail<-SeqLocTrimPolyAtail<-compute_alignment  */


#define MAX_HANG 100 /* used in need_recompute<-compute_alignment */

#define MAX_EXTEND_OVERHANG	50 /* used by compute_alignment and SeqAlignSetFlobalFromLocal */

static Boolean is_bad_align (BioseqPtr m_bsp, BioseqPtr s_bsp, SeqAlignPtr align, Int4 min_align_len,Int4 loc_len);

static void OrderInt4(Int4Ptr x, Int4Ptr y) /* replace jzmisc: swap */
{
  Int4 temp;

	if((*x) > (*y)){
	  temp = *x;
	  *x = *y;
	  *y = temp;
	}
}

static void SeqAnnotWrite(SeqAlignPtr align)
{
        SeqAnnotPtr annot;
        AsnIoPtr aip;

        annot = SeqAnnotNew();
        annot->type = 2;
        annot->data = align;

        aip = AsnIoOpen("temp2.sat", "w");
        SeqAnnotAsnWrite(annot, aip, NULL);
        AsnIoClose(aip);

        annot->data = NULL;
        SeqAnnotFree(annot);
}

static Boolean check_align_match_diagnol(SeqAlignPtr align, ValNodePtr ddp_list, Uint1 strand)
{
	DenseSegPtr dsp;
	Int2 i;
	Int4 m_start, s_start;
	Boolean minus;
	Int4 total_len, match_len;
	ValNodePtr curr;
	DenseDiagPtr ddp;

	/*check for orientations first*/
	dsp = (DenseSegPtr) align->segs;
	minus = FALSE;
	if(dsp->strands !=NULL)
	{
		if(dsp->strands[0] != dsp->strands[1])
		{
			if(dsp->strands[1] == Seq_strand_minus  || 
				dsp->strands[0] == Seq_strand_minus)
				minus = TRUE;
		}
	}

	if(strand == Seq_strand_minus && !minus)
		return FALSE;
	if(minus && strand != Seq_strand_minus)
		return FALSE;

	total_len = 0;
	match_len = 0;
	for(i = 0; i<dsp->numseg; ++i)
	{
		m_start = dsp->starts[2*i];
		s_start = dsp->starts[2*i+1];
		if(m_start != -1 && s_start != -1 && dsp->lens[i] > MIN_SEG_SIZE)
		{
			total_len += dsp->lens[i];
			for(curr = ddp_list; curr != NULL; curr = curr->next)
			{
				ddp = (DenseDiagPtr) curr->data.ptrvalue;
				if(ddp->starts[0] == m_start 
					&& ddp->starts[1] == s_start)
				{
					match_len += dsp->lens[i];
					break;
				}
			}
		}
	}

	if(match_len < 2*MIN_SEG_SIZE || total_len < 2*MIN_SEG_SIZE)
		return FALSE;

	if(match_len == total_len)
		return TRUE;

	/*at least 40% need to be matching the fragment*/
	if(match_len > 4 * MIN_SEG_SIZE)
	{
		if((match_len *100 )/40 >= total_len)
			return TRUE;
	}
	else if(total_len < 4* MIN_SEG_SIZE)
		if((match_len *100 )/80 >= total_len)
			return TRUE;

	return FALSE;
}

				
			
/*
*
*	Filter the Dense-seg alignment. Get rid of the alignments that were not 
*	in the main diagnol
*	all the alignments MUST be for the SAME sequences
*	sip is the id of the sequence that is used for the primary ordering
*/
static Boolean ModFilterDenseSegAlign(SeqAlignPtr PNTR align, SeqIdPtr sip)
{
	SeqAlignPtr ddp_align;
	ValNodePtr head = NULL;
	SeqAlignPtr curr, prev, next;
	Uint1 strand;
        /* 	DenseSegPtr dsp; */
	BioseqPtr bsp;

	if(align == NULL || *align == NULL || sip == NULL)
		return FALSE;
	/*check if it has more than two alignments*/
	curr = *align;
	if(curr->segtype != 2)
		return FALSE;
	/* dsp = curr->segs;
	if(dsp->ids->next->choice == SEQID_LOCAL)
	{
           	SeqIdPtr t_sip;
          	ObjectIdPtr oip;
		t_sip = dsp->ids->next;
		oip = t_sip->data.ptrvalue;
		printf("%s\n", oip->str);
	} */
	

	ddp_align = SeqAlignConvertDspToDdpList(*align);
	if(ddp_align == NULL)
		return FALSE;

	bsp = BioseqLockById(sip);
	if(bsp == NULL)
		return FALSE;
	head = FilterSeqAlign(ddp_align, sip, &strand);
	BioseqUnlock(bsp);
	if(head == NULL)
	{
		SeqAlignFree(ddp_align);
		return FALSE;
	}

	curr = *align;
	prev = NULL;
	while(curr)
	{
		next = curr->next;
		if(!check_align_match_diagnol(curr, head, strand))
		{
			if(prev == NULL)
				*align = next;
			else
				prev->next = next;
			curr->next = NULL;
			SeqAlignFree(curr);
		}
		else
			prev = curr;
		curr = next;
	}

	ValNodeFree(head);
	SeqAlignFree(ddp_align);

	return TRUE;
}

static Boolean reload_seq_loc(SeqLocPtr slp, Int4 from, Int4 to, Uint1 strand)
{
	SeqIntPtr sint;

	if(slp == NULL || slp->choice != SEQLOC_INT)
		return FALSE;
	sint = (SeqIntPtr) slp->data.ptrvalue;
	sint->from = from;
	sint->to = to;
	sint->strand = strand;

	return TRUE;
}
/*
*
*	modify the locations on loc_1 and loc_2 to record the from, to for the 
*	end segment. if (head), return the first segment >= min_seg_len, else 
*	return the last segment
*/
static Boolean get_end_seg(DenseSegPtr dsp, Int4 min_seg_len, SeqLocPtr loc_1, SeqLocPtr loc_2, Boolean head)
{
	Int4 i;
	Int4 start_1, start_2;

	if(dsp == NULL || loc_1 == NULL || loc_2 == NULL)
		return FALSE;
	if(head)
		i = 0;
	else
		i = dsp->numseg -1;

	for(; head? i<dsp->numseg : i>=0; head? ++i: --i)
	{
		start_1 = dsp->starts[2*i];
		start_2 = dsp->starts[2*i+1];
		if(start_1 != -1 && start_2 != -1)
		{
			if(dsp->lens[i] >= min_seg_len)
			{
				reload_seq_loc(loc_1, start_1, 
					start_1+dsp->lens[i]-1, dsp->strands[0]);
				reload_seq_loc(loc_2, start_2, 
					start_2+dsp->lens[i]-1, dsp->strands[1]);

				return TRUE;
			}
		}
	}

	return FALSE;
}


/*
*
*	align_1 and align_2 will be  re-sorted by the ascending order in the first 
*	sequence in the alignment. The first sequence needs to aligned 
*	in the plus strand. They have to be alignment on the same 
*	diagnol
*	return NULL if any merge fails
*
*/
static SeqAlignPtr MergeTwoDspBySIM4(SeqAlignPtr align_1, SeqAlignPtr align_2,Int4 MaxGap)
{
	SeqLocPtr loc_1_1, loc_1_2, loc_2_1, loc_2_2;
	Int4 start, stop;
	Int4 match_start, match_stop;
	Boolean found = FALSE;
	BioseqPtr bsp;
	SeqAlignPtr align, t_align;
	DenseSegPtr dsp_1, dsp_2;
	BioseqPtr bsp_1, bsp_2;
	Int4 max_align_size;
	Uint1 strand;
	Int4 gap=0,gap_1, gap_2;
	Int4 start_1, stop_1, start_2, stop_2;
	SeqIdPtr sip;
	Boolean reverse = FALSE;
	Boolean has_overlap;
	Int2 order;

	if(align_1 == NULL || align_2 == NULL)
		return NULL;

	if(align_1->segtype != 2 || align_2->segtype != 2)
		return NULL;
        if(SeqAlignStart(align_1,0)>SeqAlignStart(align_2,0)) {
            VoidPtr segs;
            segs = align_1->segs;
            align_1->segs = align_2->segs;
            align_2->segs = segs;
        }

	dsp_1 = (DenseSegPtr) align_1->segs;
	if(dsp_1 && dsp_1->strands &&  dsp_1->strands[0] == Seq_strand_minus)
	{
		SeqAlignReverse(align_1, 0);
		reverse = TRUE;
	}
	dsp_2 = (DenseSegPtr) align_2->segs;
	if(dsp_2 && dsp_2->strands &&  dsp_2->strands[0] == Seq_strand_minus)
	{
		SeqAlignReverse(align_2, 0);
		reverse = TRUE;
	}

	max_align_size = MIN(SeqAlignLength(align_1), 
		SeqAlignLength(align_2));

	if(dsp_1 == NULL || dsp_2 == NULL)
		return NULL;

	SeqAlignStartStopById(align_1, dsp_1->ids, &start_1, &stop_1, &strand); 
	SeqAlignStartStopById(align_2, dsp_1->ids, &start_2, &stop_2, &strand); 
        /* Redundant alignment covering same region of query */
	if(start_2 >= start_1 && stop_2 <= stop_1)
		return NULL;

	SeqAlignStartStopById(align_1, dsp_1->ids->next, &start_1, &stop_1, &strand); 
	SeqAlignStartStopById(align_2, dsp_1->ids->next, &start_2, &stop_2, &strand); 
        /* Redundant alignment covering same region of subject */
	if(start_2 >= start_1 && stop_2 <= stop_1)
		return NULL;

	SeqAlignStartStopById(align_1, dsp_1->ids, &start, &stop, &strand); 
	loc_1_1 = SeqLocIntNew(start, stop, strand, dsp_1->ids);
	SeqAlignStartStopById(align_1, dsp_1->ids->next, &start, &stop, &strand); 
	loc_1_2 = SeqLocIntNew(start, stop, strand, dsp_1->ids->next);

	/* loc_1_1 = SeqLocIntNew(0, -1, 0, dsp_1->ids);
	loc_1_2 = SeqLocIntNew(0, -1, 0, dsp_1->ids->next);
	if(!get_end_seg(dsp_1, MIN_SEG_SIZE*2, loc_1_1, loc_1_2, FALSE))
	{
		SeqLocFree(loc_1_1);
		SeqLocFree(loc_1_2);
		return NULL;
	} */

	loc_2_1 = SeqLocIntNew(0, -1, 0, dsp_1->ids);
	loc_2_2 = SeqLocIntNew(0, -1, 0, dsp_1->ids->next);
	if(!get_end_seg(dsp_2, MIN_SEG_SIZE*2, loc_2_1, loc_2_2, TRUE))
	{
		SeqLocFree(loc_1_1);
		SeqLocFree(loc_1_2);
		SeqLocFree(loc_2_1);
		SeqLocFree(loc_2_2);
		return NULL;
	}

	found = FALSE;
	has_overlap = FALSE;
	if(select_overlap_loc(loc_1_1, loc_2_1, loc_1_2, loc_2_2, &order))
	{
		if(order == 0)
			is_loc_overlap(loc_1_1, loc_2_1, &start, &stop, 0, max_align_size, NULL);
		else
			is_loc_overlap(loc_1_2, loc_2_2, &start, &stop, 0, max_align_size, NULL);
		found = TRUE;
		has_overlap = TRUE;
	}

	if(!found)
	{
		SeqLocFree(loc_1_1);
		SeqLocFree(loc_1_2);
		loc_1_1 = SeqLocIntNew(0, -1, 0, dsp_1->ids);
		loc_1_2 = SeqLocIntNew(0, -1, 0, dsp_1->ids->next);
		if(!get_end_seg(dsp_1, MIN_SEG_SIZE*2, loc_1_1, loc_1_2, FALSE))
		{
			SeqLocFree(loc_1_1);
			SeqLocFree(loc_1_2);
			SeqLocFree(loc_2_1);
			SeqLocFree(loc_2_2);
			return NULL;
		} 
		/*check for the overlap between the two segments*/
		if(is_loc_overlap(loc_1_1, loc_2_1, &start, &stop, 0, max_align_size, NULL))
		{
			found = TRUE;
			order = 0;
		}
		else if(is_loc_overlap(loc_1_2, loc_2_2, &start, &stop, 0, max_align_size, NULL))
		{
			found = TRUE;
			order = 1;
		}
		else
		{
			start_1 = start;
			stop_1 = stop;
			start_2 = start;
			stop_2 = stop;

			if(is_loc_overlap(loc_1_1, loc_2_1, &start_1, &stop_1, MIN_SEG_SIZE*10, max_align_size, &gap_1))
			{
				found = TRUE;
				order = 0;
                                gap = gap_1;
			}
			if(is_loc_overlap(loc_1_2, loc_2_2, &start_2, &stop_2, MIN_SEG_SIZE*10, max_align_size, &gap_2))
			{
				if(!found || gap_2 < gap_1)
				{
					order = 1;
					found = TRUE;
                                        gap = gap_2;
				}
			}
			if(found)
			{
				if(order == 0)
				{
					start = start_1;
					stop = stop_1;
				}
				else
				{
					start = start_2;
					stop = stop_2;
				}
			}
		}
	}

	SeqLocFree(loc_1_1);
	SeqLocFree(loc_1_2);
	SeqLocFree(loc_2_1);
	SeqLocFree(loc_2_2);

	if(!found)
	{
		return NULL;
	}

	/*
		===================>
			====================>
			^^^^^^^^^^^^^^^^^
			           |----------
                                    need to have some offset into the 
				    second sequence
	*/
	if(order == 0)
		sip = dsp_1->ids;
	else
		sip = dsp_1->ids->next;
	SeqAlignStartStopById(align_1, sip, &start_1, &stop_1, &strand); 
	SeqAlignStartStopById(align_2, sip, &start_2, &stop_2, &strand); 

        start = MAX(start_2, start_1);
        stop = MIN(stop_1, stop_2);
        OrderInt4(&start,&stop); /* Need reordering for gaps */

	if(dsp_1->strands[order] == Seq_strand_minus)
	{
		match_start = find_matching_position(dsp_1, stop, order);
		match_stop = find_matching_position(dsp_2, start, order);
	}
	else
	{
		match_start = find_matching_position(dsp_1, start, order);
		match_stop = find_matching_position(dsp_2, stop, order);
	}

	if(match_start == -1 || match_stop == -1)
		return NULL;
	bsp = NULL;
	if(order == 0)
		bsp = BioseqLockById(dsp_1->ids->next);
	else
		bsp = BioseqLockById(dsp_1->ids);

	if(bsp == NULL)
		return NULL;
        { 
            Int4 old_start=start,old_stop = stop;
            start = MAX(start-100,MIN(start_1,start_2));
            stop = MIN(stop+100,MAX(stop_1,stop_2));
            if(match_start<=match_stop) {
                match_start -=(old_start-start);
                match_stop +=(stop-old_stop);
                if(match_start<0) {
                    start+=match_start;
                    match_start=0;
                } 
                if(match_stop>bsp->length-1) {
                    stop-=(match_stop-(bsp->length-1));
                    match_stop = bsp->length-1;
                }
            }
            else {
                match_start +=(old_start-start);
                match_stop -=(stop-old_stop);
                if(match_stop<0) {
                    stop+=match_stop;
                    match_stop=0;
                } 
                if(match_start>bsp->length-1) {
                    start-=(match_start-(bsp->length-1));
                    match_start = bsp->length-1;
                }

            }

        }
	BioseqUnlock(bsp);


	OrderInt4(&match_start, &match_stop);
	OrderInt4(&start, &stop);
        /* HS Fix to Disallow Very Long Gaps */
        if(abs(abs(stop-start) -abs(match_stop-match_start))>MaxGap+3000)
            return NULL;
	if(order == 0)
	{
		loc_1_1 = SeqLocIntNew(start, stop, dsp_1->strands[order], dsp_1->ids);
		loc_1_2 = SeqLocIntNew(match_start, match_stop, 
			dsp_1->strands[1-order], dsp_1->ids->next);
	}
	else
	{
		loc_1_1 = SeqLocIntNew(match_start, match_stop, 
			dsp_1->strands[1-order], dsp_1->ids);
		loc_1_2 = SeqLocIntNew(start, stop, 
			dsp_1->strands[order], dsp_1->ids->next);
	}

	bsp_1 = BioseqLockById(SeqLocId(loc_1_1));
	if(bsp_1 == NULL)
	{
		SeqLocFree(loc_1_1);
		SeqLocFree(loc_1_2);
		return NULL;
	}
	bsp_2 = BioseqLockById(SeqLocId(loc_1_2));
	if(bsp_2 == NULL)
	{
		BioseqUnlock(bsp_1);
		SeqLocFree(loc_1_1);
		SeqLocFree(loc_1_2);
		return NULL;
	}

	align = SIM4ALN_choice(loc_1_1, loc_1_2, 200, 8);
	SeqLocFree(loc_1_1);
	SeqLocFree(loc_1_2);

	if(align != NULL)
	{
		if(order == 1)
		{
			start = match_start;
			stop = match_stop;
			order = 0;
		}
		t_align = (SeqAlignPtr) AsnIoMemCopy((Pointer)align_1, (AsnReadFunc)SeqAlignAsnRead, (AsnWriteFunc)SeqAlignAsnWrite);

		/* for debugging purpose*/
		/* t_align->next = align;
		align->next = (SeqAlignPtr) AsnIoMemCopy((Pointer)align_2, (AsnReadFunc)SeqAlignAsnRead, (AsnWriteFunc)SeqAlignAsnWrite);
		return t_align;    */

		if(MergeTwoAlignList(t_align, &align, start, stop, order))
		{
			if(align != NULL)
				SeqAlignFree(align);
			align = t_align;
			align->next = NULL;

			t_align = (SeqAlignPtr) AsnIoMemCopy((Pointer)align_2, (AsnReadFunc)SeqAlignAsnRead, (AsnWriteFunc)SeqAlignAsnWrite);
			t_align->next = NULL;
			MergeTwoAlignList(align, &t_align, start, stop, order);
			if(t_align != NULL)
				SeqAlignFree(t_align);
		}
		else
		{
			if(align != NULL)
				align = SeqAlignSetFree(align);
			SeqAlignFree(t_align);
		}
	}
	/*for debugging */
	/* else
	{
		align_1->next = align_2;
		align_2->next = NULL;
		SeqAnnotWrite(align_1);
		exit(1);
	} */
	BioseqUnlock(bsp_1);
	BioseqUnlock(bsp_2);

	if(reverse && align != NULL)
		SeqAlignReverse(align, 1);
	return align;
}

		

static SeqAlignPtr MergeToOneAlignment(SeqAlignPtr PNTR palign,Int4 MaxGap) 
{
	SeqAlignPtr next, prev, align;
	SeqAlignPtr merge_align;
	SeqAlignPtr curr;

	if(palign == NULL || *palign == NULL)
		return NULL;
	align = *palign;
	if(align->segtype != 2)
		return NULL;
	*palign = SeqAlignSortByRegion(align, 0);

	prev = NULL;
	align = *palign;

	while(align)
	{
		next = align->next;
		if(next != NULL)
		{
			merge_align = NULL;
			curr = *palign;
			/* if(next->next == NULL)
				printf("stop here\n"); */
			merge_align = MergeTwoDspBySIM4(align, next,MaxGap);
			if(merge_align != NULL)
			{
				if(prev == NULL)
					*palign = merge_align;
				else
					prev->next = merge_align;
				
				curr = merge_align;
				while(curr->next != NULL)
				{
					prev = curr;
					curr = curr->next;
				}
				curr->next = next->next;
				next->next = NULL;
				SeqAlignSetFree(align);
				align = curr;
			}
			else
			{
				prev = align;
				align = align->next;
			}
		}
		
		else
			align = align->next;
	}

	if(*palign != NULL && (*palign)->next!=NULL)
		*palign = SeqAlignSortByLength(*palign);
	return (*palign);
}

/*
*
*	Filter out any alignments that were not in the main 
*	diagnol
*
*/
static Boolean ModifyAlignList(SeqAlignPtr PNTR palign)
{
	SeqAlignPtr align, next;
	DenseSegPtr dsp;
	SeqAlignPtr h_align = NULL;

	if(palign == NULL || *palign == NULL)
		return FALSE;

	align = *palign;
	while(align)
	{
		next = align->next;
		if(align->segtype == 2)
		{
			dsp = (DenseSegPtr) align->segs;
			align->next = SeqAlignExtractByIds(&next,dsp->ids, dsp->ids->next);
			ModFilterDenseSegAlign(&align, dsp->ids);
			if(align != NULL)
				h_align = SeqAlignLink(align, h_align);
		}
		else
		{
			align->next = NULL;
			SeqAlignFree(align);
		}
		align = next;
	}

	*palign = h_align;
	return TRUE;
}

				

			


static Int2 get_polyA_index (Uint1Ptr res_buf, Int2 len)
{
	Int2 i, total;
	Uint1 res;
	Int2 count_A;

	count_A = 0;
	total = 0;

	for( i= len -1; i>=0; --i)
	{
		res = res_buf[i];
		if(res == 'a' || res == 'A' || res == 'n' || res == 'N')
		{
			++count_A;
			++total;
		}
		else
		{
			if(total > 5 && count_A * 100 < (total + 1) * PERCENT_A_RES)
				break;
			else
				++total;
		}
	}
	/*trim the non-A and non-Ns*/
	for(i = MAX(i, 0); i<len; )
	{
		res = res_buf[i];
		if(res != 'n' && res != 'N' && res != 'a' && res != 'A')
			++i;
		else	/*N/A residue*/
			break;
	}

	if(len-1 - i > 2)
		return i;

	return -1;
}


static SeqLocPtr check_polyA_tail (BioseqPtr bsp)
{
	Int2 len = 100;
	SeqPortPtr spp;
	Uint1 res_buf[101];
	Uint1 res;
	Int2 i;
	SeqLocPtr end_loc, begin_loc;
	Int4 start, stop;

	begin_loc = NULL;
	end_loc = NULL;

	len = MIN(bsp->length-1, 100);
	/*check for the end of the sequence for polyA tail*/
	spp = SeqPortNew(bsp, bsp->length-1-(len-1), bsp->length-1, 
		Seq_strand_plus, Seq_code_iupacna);
	i = 0;
	while((res = SeqPortGetResidue(spp)) != SEQPORT_EOF)
	{
		if(IS_ALPHA(res))
			res_buf[i++] = res;
	}
	SeqPortFree(spp);
	if(i > 10)
	{
		len = i;
		i = get_polyA_index (res_buf, len);
		if(i != -1)
		{
			stop = bsp->length -1;
			start = bsp->length - 1 - (len -1) + i;

			end_loc = SeqLocIntNew(start, stop, 
				Seq_strand_plus, SeqIdFindBest(bsp->id, 0));
		}

	}

	/*check for the beginning on the reverse complement of the sequence 
          for polyA tail. Usually the 3' clone of the EST*/
	spp = SeqPortNew(bsp, 0, len-1, 
		Seq_strand_minus, Seq_code_iupacna);
	i = 0;
	while((res = SeqPortGetResidue(spp)) != SEQPORT_EOF)
	{
		if(IS_ALPHA(res))
			res_buf[i++] = res;
	}
	SeqPortFree(spp);
	if(i > 10)
	{
		len = i;
		i = get_polyA_index (res_buf, len);
		if(i != -1)
		{
			stop = len -1 - i;
			start = 0;

			begin_loc = SeqLocIntNew(start, stop, 
				Seq_strand_minus, SeqIdFindBest(bsp->id, 0));
		}

	}
	if(begin_loc == NULL)
		return end_loc;
	else
	{
		if(end_loc != NULL)
		{
			if(SeqLocLen(end_loc) > SeqLocLen(begin_loc))
			{
				SeqLocFree(begin_loc);
				return end_loc;
			}
			else
			{
				SeqLocFree(end_loc);
				return begin_loc;
			}
		}
		return begin_loc;
	}
}

static void SeqLocTrimPolyATail(SeqLocPtr loc_2, BioseqPtr bsp_2) {
	SeqLocPtr poly_loc;
	SeqIntPtr sint;

	poly_loc = check_polyA_tail (bsp_2);
	if(poly_loc != NULL)
	{
		sint = (SeqIntPtr) loc_2->data.ptrvalue;
		if(SeqLocStart(poly_loc) > (bsp_2->length)/2)
		    sint->to = SeqLocStart(poly_loc) -1;
		else
			sint->from = SeqLocStop(poly_loc) + 1;
		SeqLocFree(poly_loc);
	}
}

static Boolean is_bad_blast_alignment(SeqAlignPtr align, SeqLocPtr slp1, SeqLocPtr slp2, Int4 min_align_len)
{
	SeqAlignPtr best_align;
	BioseqPtr bsp1, bsp2;
	Int4 max_score;

	best_align = find_best_align(align, &max_score);
	if(best_align == NULL)
		return TRUE;
	bsp1 = BioseqFind(SeqLocId(slp1));
	bsp2 = BioseqFind(SeqLocId(slp2));

	if(bsp1 == NULL || bsp2 == NULL)
		return (max_score > 50);
	return is_bad_align(bsp1, bsp2, best_align,
		min_align_len, MIN(SeqLocLen(slp1), SeqLocLen(slp2)));
}

NLM_EXTERN SeqAlignPtr SeqAlignSplitBlastTwoSeq(SeqLocPtr slp1, SeqLocPtr slp2, 
		Int4 split_len, Int4 overlap_len, BLAST_OptionsBlkPtr options)
{
	Boolean split = FALSE;
	Int4 len1, len2;
	SeqAlignPtr align = NULL, t_align;
	Int4 num, i;
	SeqLocPtr slp;
	SeqIntPtr sint;
	Int4 pstop;


	if(slp1 == NULL || slp2 == NULL)
		return NULL;
	len1 = SeqLocLen(slp1);
	len2 = SeqLocLen(slp2);
	if(len1 == 0 || len2 == 0)
		return NULL;


	if(MAX(len1, len2) > 100000 && MIN(len1, len2) > split_len)
		split = TRUE;

	align = NULL;
	if(split)
	{
		num = (len1 - overlap_len)/(split_len - overlap_len);
		if((len1 - overlap_len)%(split_len - overlap_len) > split_len/2)
			++num;
		slp = SeqLocIntNew(SeqLocStart(slp1), SeqLocStop(slp1), 
			SeqLocStrand(slp1), SeqLocId(slp1));
		sint = (SeqIntPtr) slp->data.ptrvalue;
		pstop = -1;
		for(i = 0; i<num; ++i)
		{
			if(i == 0)
			{
				sint->from = SeqLocStart(slp1);
				sint->to = sint->from + split_len -1;
			}
			else
			{
				sint->from = sint->to -overlap_len;
				sint->to = MIN(sint->from + split_len -1, SeqLocStop(slp1));
			}
			if(i == num -1)
				sint->to = SeqLocStop(slp1);
			sint->strand = SeqLocStrand(slp1);
			t_align = BlastTwoSequencesByLoc(slp, slp2, "blastn", options);
			if(t_align != NULL)
			{
				if(align == NULL)
					align = t_align;
				else if(pstop != -1)
				{
					MergeTwoAlignList(align, &t_align, 
						sint->from, pstop, 0);
					if(t_align != NULL)
						align = SeqAlignLink(t_align, align);
				}
				pstop = sint->to;
			}
		}
		SeqLocFree(slp);
	} else
	    align = BlastTwoSequencesByLoc(slp1, slp2, "blastn", options);

	return align;
}

NLM_EXTERN PSeqAlignInfoPtr SeqAlignToPSeqAlignInfo (SeqAlignPtr sap)
{
        PSeqAlignInfoPtr alip, alip_head, newalip;
        Int2 found;
        SeqAlignPtr alip_sap;
        SeqIdPtr sip;

        if (!sap)
                return NULL;
        alip = (PSeqAlignInfoPtr)MemNew(sizeof(PSeqAlignInfo));
        alip_head = alip;
        sip = SeqIdPtrFromSeqAlign(sap);
        alip->sap = sap;
        alip->sip = SeqIdDupList(sip);
        alip->next = NULL;
        sap = sap->next;
        alip->sap->next=NULL;
        while (sap)
        {
                sip = SeqIdPtrFromSeqAlign(sap);
                found = 0;
                alip = alip_head;
                while (alip && !found)
                {
                        if (SeqIdComp(sip->next, alip->sip->next))
                        {
                                alip_sap = alip->sap;
                                while(alip_sap->next !=NULL)
                                        alip_sap = alip_sap->next;
                                alip_sap->next = sap;
                                sap = sap->next;
                                alip_sap->next->next = NULL;
                                found = 1;
                        }
                        alip = alip->next;
                }
                if (found == 0)
                {
                        alip = alip_head;
                        while(alip->next != NULL)
                                alip = alip->next;
                        newalip = (PSeqAlignInfoPtr)MemNew(sizeof(PSeqAlignInfo));
                        newalip->sip = SeqIdDupList(sip);
                        newalip->sap = sap;
                        newalip->next = NULL;
                        sap = sap->next;
                        newalip->sap->next = NULL;
                        alip->next = newalip;
                }
        }
        return alip_head;
}

NLM_EXTERN SeqAlignPtr ReassembleSeqAlignFromPSeqAlignInfo(PSeqAlignInfoPtr alip)
{
        SeqAlignPtr sap, head_sap, sap_curr;
        PSeqAlignInfoPtr prevalip;

        head_sap = NULL;
        while(alip)
        {
                if (!head_sap)
                {
                        head_sap = alip->sap;
                        sap_curr = head_sap;
                } else
                {
                        sap = alip->sap;
                        sap_curr->next = sap;
                        while (sap_curr->next)
                                sap_curr = sap_curr->next;
                }
                prevalip = alip;
                alip = alip->next;
                prevalip->sip = SeqIdSetFree(prevalip->sip);
                MemFree(prevalip);
        }
        return head_sap;
}

NLM_EXTERN SeqAlignPtr SeqAlignSplitGappedBlast(SeqLocPtr slp1, CharPtr progname, CharPtr database, ValNodePtr *other_returns, ValNodePtr *error_returns, Int4 split_len, Int4 overlap_len, BLAST_OptionsBlkPtr options)
{
        Boolean split = FALSE;
        Int4 len1;
        SeqAlignPtr t_align;
        Int4 num, i;
        SeqLocPtr slp;
        SeqIntPtr sint;
        Int4 pstop;
        PSeqAlignInfoPtr alip, t_alip, head_alip=NULL, head_t_alip=NULL;
        Int2 beg=0;


        if(slp1 == NULL)
                return NULL;
        len1 = SeqLocLen(slp1);
        if(len1 == 0)
                return NULL;


        if(len1 > 10000 && len1 > split_len)
                split = TRUE;

        alip = NULL;
        if(split)
        {
                num = (len1 - overlap_len)/(split_len - overlap_len);
                if((len1 - overlap_len)%(split_len - overlap_len) > split_len/2)
                        ++num;
                slp = SeqLocIntNew(SeqLocStart(slp1), SeqLocStop(slp1),
                        SeqLocStrand(slp1), SeqLocId(slp1));
                sint = (SeqIntPtr) slp->data.ptrvalue;
                pstop = -1;
                for(i = 0; i<num; ++i)
                {

                        if(i == 0)
                        {
                                sint->from = SeqLocStart(slp1);
                                sint->to = sint->from + split_len -1;
                        }
                        else
                        {
                                sint->from = sint->to -overlap_len;
                                sint->to = MIN(sint->from + split_len -1, SeqLocStop(slp1));
                        }
                        if(i == num -1)
                                sint->to = SeqLocStop(slp1);
                        sint->strand = SeqLocStrand(slp1);
                        t_align = BioseqBlastEngineByLoc(slp, progname, database, options, other_returns, error_returns, NULL);
                        if(t_align != NULL)
                        {
                                t_alip = SeqAlignToPSeqAlignInfo(t_align);
                                if(alip == NULL)
                                {
                                        alip = t_alip;
                                        head_alip = alip;
                                }
                                else if(pstop != -1)
                                {
                                        alip = head_alip;
                                        head_t_alip = t_alip;
                                        while (alip != NULL)
                                        {
                                                t_alip = head_t_alip;
                                                while (t_alip != NULL)
                                                {
                                                        if (SeqIdComp(t_alip->sip->next, alip->sip->next)

                                                                && t_alip->used == FALSE)
                                                        {
                                                                MergeTwoAlignList(alip->sap, &t_alip->sap,
                                                                        sint->from, pstop, 0);
                                                                if (t_alip->sap !=NULL)
                                                                {
                                                                        alip->sap = SeqAlignLink(t_alip->sap,
                                                                                alip->sap);
                                                                }
                                                                t_alip->used = TRUE;
                                                        }
                                                        t_alip = t_alip->next;
                                                }
                                                alip = alip->next;
                                        }
                                        alip = head_alip;
                                        while (alip->next != NULL)
                                                alip = alip->next;
                                        t_alip = head_t_alip;
                                        while (t_alip)
                                        {
                                                if (t_alip->used == FALSE)
                                                {
                                                        alip->next = t_alip;
                                                        t_alip = t_alip->next;
                                                        alip = alip->next;
                                                        alip->next = NULL;
                                                } else
                                                        t_alip = t_alip->next;
                                        }
                                        t_alip = head_t_alip = NULL;
                                }
                                pstop = sint->to;
                        }
                }
                SeqLocFree(slp);
                t_align = ReassembleSeqAlignFromPSeqAlignInfo(head_alip);
        } else
        {
            t_align = BioseqBlastEngineByLoc(slp1, progname, database, options, other_returns, error_returns, NULL);
        }

        return t_align;
}
/* Filter SeqAlign by removing bad hits and hits off the main diagonal (as defined by the
   first hit .. which is the highest scoring in blast */
static SeqAlignPtr SeqAlignConsistentDiagFilter(SeqAlignPtr align,SeqLocPtr slp1, SeqLocPtr slp2, 
					 FILE* err_fp,Int4 MaxGap) {
    Int4 old_num,curr_num,len,len1,len2;

    SeqAlignPtr prev,t_align;

    if(slp1 == NULL || slp2 == NULL || align == NULL)
	return NULL;
    len1 = SeqLocLen(slp1);
    len2 = SeqLocLen(slp2);
    if(len1 == 0 || len2 == 0)
	return NULL;

    if(align != NULL) {
#ifdef DEBUG
	    save_output (SeqLocId(slp1), "blast.out", align);
#endif
		/*
		  Checks for harcoded conditions on the FIRST alignment. 
		  Alignment at least
		  50 bases or 90% of length if sequences are smaller than 50.
                  Also do not allow more than 10% mismatches.
		  */
	if(is_bad_blast_alignment(align, slp1, slp2, 50)) {
	    SeqAlignSetFree(align);
	    return NULL;
	}
	    /* Remove hits that are off the main diagonal of the
		   longest hit */
	ModifyAlignList(&align);
#ifdef DEBUG
	    save_output (SeqLocId(slp1), "merge.out", align);
#endif
	if(is_bad_blast_alignment(align, slp1, slp2, 50))
	    {
		SeqAlignSetFree(align);
		return NULL;
	    }
	
	if(align != NULL)
	    {
		old_num = SeqAlignCount(align);
		MergeToOneAlignment(&align,MaxGap);
		curr_num = SeqAlignCount(align);
		
		if(align == NULL)
		    {
			fprintf(stderr, "Fail in MergeToOneAlignment\n");
		    } 
		while(curr_num > 1 && curr_num < old_num)
		    {
			old_num = curr_num;
			MergeToOneAlignment(&align,MaxGap);
			curr_num = SeqAlignCount(align);
		    }
		
		if(curr_num > 1)
		    {
			prev = align;
			t_align = align->next;
			while(t_align)
			    {
				len = SeqAlignLength(t_align);
				if(len < 50 || SeqAlignLength(prev)/len > 5)
				    {
					prev->next = NULL;
					SeqAlignSetFree(t_align);
						break;
				    }
				else
				    {
					prev = t_align;
						t_align = t_align->next;
				    }
			    }
		    }
						
		
	    }
#ifdef DEBUG
	    save_output (SeqLocId(slp1), "sim4.out", align);
#endif
	if(is_bad_blast_alignment(align, slp1, slp2, MIN(50, MIN(len1, len2)/5)))
	    {
		SeqAlignSetFree(align);
		return NULL;
	    }
	
    }
    return align;
}

static SeqAlignPtr SplitBlastTwoSeq(SeqLocPtr slp1, SeqLocPtr slp2, 
		Int4 split_len, Int4 overlap_len, BLAST_OptionsBlkPtr options, FILE* err_fp)
{
	Int4 len1, len2;
        BioseqPtr bsp,t_bsp;
	SeqAlignPtr align = NULL;
        Int4 MaxGap;
        MaxGap = (Int4)(options->gap_x_dropoff*0.6);

	if(slp1 == NULL || slp2 == NULL)
		return NULL;
	len1 = SeqLocLen(slp1);
	len2 = SeqLocLen(slp2);
	if(len1 == 0 || len2 == 0)
		return NULL;
	/* Align the sequence "a la Powerblast" */
	align = SeqAlignSplitBlastTwoSeq(slp1, slp2,split_len, overlap_len, options);
        bsp = BioseqLockById(SeqLocId(slp1));
        t_bsp = BioseqLockById(SeqLocId(slp2));
        if(bsp && t_bsp) {
            align = SeqAlignSetGlobalFromLocal(align,slp1, slp2, bsp, t_bsp, err_fp,  (Int4)(options->gap_x_dropoff*0.6));
            /* Filter out anything that is off the main diagonal */
            align = SeqAlignConsistentDiagFilter(align,slp1,slp2,err_fp,MaxGap);
        }
        if(bsp)
            BioseqUnlock(bsp);
        if(t_bsp)
            BioseqUnlock(t_bsp);
	return align;
}




static Boolean is_bad_align (BioseqPtr m_bsp, BioseqPtr s_bsp, SeqAlignPtr align, Int4 min_align_len,Int4 loc_len)
{
	DenseSegPtr dsp;
	Int4 m_start, m_stop, s_start, s_stop;
	SeqPortPtr m_spp, s_spp;
	Int2  i;
	Uint1 m_res, s_res;
	Int4 mismatch = 0;
	Uint1 code;
	Int4 j;
	Uint1 strand;
	Int4 align_len;

	code = Seq_code_iupacna;
	dsp = (DenseSegPtr) align->segs;
	SeqAlignStartStopById(align, dsp->ids->next, &s_start, &s_stop, &strand);
	if(s_start<0 || s_stop>= s_bsp->length) {
	    ErrPostEx(SEV_WARNING,0,0,"find_mismatch_residue: corrupted alignment \n");
	    return TRUE;
	} else {
	    s_spp = SeqPortNew(s_bsp, s_start, s_stop, dsp->strands[1], code);
	}
	/* if(s_stop - s_start + 1 < MIN(min_align_len, loc_len-MAX(10, loc_len/10)))
		return TRUE; */

	SeqAlignStartStopById(align, dsp->ids, &m_start, &m_stop, &strand);
	if(m_start<0 || m_stop>= m_bsp->length) {
	    ErrPostEx(SEV_WARNING,0,0,"find_mismatch_residue: corrupted alignment \n");
	    return TRUE;
	} else {

	    m_spp = SeqPortNew(m_bsp, m_start, m_stop, strand, code);
	}


	align_len = 0;
	for(i = 0; i<dsp->numseg; ++i)
	{
		m_start = dsp->starts[2*i];
		s_start = dsp->starts[2*i + 1];
		if(m_start == -1)
		{
			if(dsp->lens[i] <=5)
			{
				mismatch += 1;
				++align_len;
			}
			if(s_start != -1)
				SeqPortSeek(s_spp, dsp->lens[i], SEEK_CUR);
		}
		else if(s_start == -1)
		{
			if(dsp->lens[i] <=5)
			{
				mismatch += 1;
				++align_len;
			}
			if(m_start != -1)
				SeqPortSeek(m_spp, dsp->lens[i], SEEK_CUR);
		}
		else
		{
			for(j = 0; j<dsp->lens[i]; ++j)
			{
				m_res = SeqPortGetResidue(m_spp);
				while(!IS_ALPHA(m_res) && m_res != SEQPORT_EOF)
					m_res = SeqPortGetResidue(m_spp);
				s_res = SeqPortGetResidue(s_spp);
				while(!IS_ALPHA(s_res) && s_res != SEQPORT_EOF)
					s_res = SeqPortGetResidue(s_spp);
				if(m_res != 'n' && m_res != 'N' && s_res != 'n' && s_res != 'N')
				{
					if(m_res != s_res)
						++mismatch;
				}
			}
			align_len += dsp->lens[i];
		}
		
	}

	SeqPortFree(m_spp);
	SeqPortFree(s_spp);
	if(align_len < MIN(min_align_len, loc_len-MAX(10, loc_len/10)))
		return TRUE;
	if(mismatch > 0)
		return ((mismatch*100)/align_len >10);
	else
		return FALSE;
}

static Boolean need_recompute(SeqAlignPtr align, Int4 m_len, Int4 s_len)
{
        DenseSegPtr dsp;
        Int4 m_start, m_stop;
        Int4 s_start, s_stop;
        Int2 end, i;
        Int2 small_seg_num;

        if(m_len > 10000 || s_len > 10000)
                return FALSE;
        dsp = (DenseSegPtr) align->segs;
        end = dsp->numseg-1;
        m_start = dsp->starts[0];
        m_stop = dsp->starts[2*end] + dsp->lens[end] -1;

        if(dsp->strands[1] == Seq_strand_minus)
        {
                s_start = dsp->starts[2*end +1];
                s_stop = dsp->starts[1] + dsp->lens[0] -1;
                if(s_start > MAX_HANG && (m_len - 1 - m_stop) > MAX_HANG)
                        return TRUE;
                if((s_len-1-s_stop)>MAX_HANG && m_start > MAX_HANG)
                        return TRUE;
        }
        else
        {
                s_start = dsp->starts[1];
                s_stop = dsp->starts[2*end+1] + dsp->lens[end] -1;
                if(s_start > MAX_HANG && m_start > MAX_HANG)
                        return TRUE;
                if((s_len-1-s_stop) > MAX_HANG && (m_len-1-m_stop)> MAX_HANG)
                        return TRUE;
        }

        small_seg_num = 0;
        for(i = 0; i<dsp->numseg; ++i)
        {
                if(dsp->lens[i] < 10)
                        ++small_seg_num;
        }

        return (small_seg_num >= 20);

}

static void extend_dsp_ends(Int4 first_len, Int4 second_len, SeqAlignPtr align, Int4 max_overhang)
{
	Int4 first, last;
	DenseSegPtr dsp;
	Int2 endseg;
	Int4 extend_len;
	Int4 offset;

	dsp = (DenseSegPtr) align->segs;
	endseg = dsp->numseg-1;
	if(dsp->strands[1] == Seq_strand_minus)
	{
		first = dsp->starts[2*endseg+1];
		if(first > 0 && first < max_overhang)
		{
                    if(dsp->strands[0] != Seq_strand_minus)
                        extend_len = MIN(first, first_len - dsp->starts[2*endseg] - dsp->lens[endseg]);
                    else
                        extend_len = MIN(first, dsp->starts[2*endseg]);
			if(extend_len > 0)
			{
				dsp->lens[endseg] += extend_len;
                                dsp->starts[2*endseg+1]-=extend_len;
                                if(dsp->strands[0]==Seq_strand_minus)
                                    dsp->starts[2*endseg] -= extend_len;
			}
		}

		last = dsp->starts[1] + dsp->lens[0];
		offset = second_len - last;
		if(offset > 0 && offset < max_overhang)
		{
                    if(dsp->strands[0]!=Seq_strand_minus) 
			extend_len = MIN(offset,dsp->starts[0]);
                    else
                        extend_len = MIN(offset,first_len - dsp->starts[0]-dsp->lens[0]);
			if(extend_len > 0)
			{
				dsp->lens[0] += extend_len;
                                if(dsp->strands[0]==Seq_strand_plus)
                                    dsp->starts[0]-=extend_len;
			}
		}
	}
	else
	{
		first = dsp->starts[1];
		if(first> 0 && first < max_overhang)
		{
                    if(dsp->strands[0]!=Seq_strand_minus) 
			extend_len = MIN(first, dsp->starts[0]);
                    else
                        extend_len = MIN(first, first_len -dsp->starts[0]-dsp->lens[0]);
			if(extend_len > 0)
			{
				dsp->lens[0] += extend_len;
				dsp->starts[1] -= extend_len;
                                if(dsp->strands[0]!=Seq_strand_minus)
                                    dsp->starts[0] -= extend_len;
			}
		}
		last = dsp->starts[2*endseg+ 1] + dsp->lens[endseg];
		offset = second_len  - last;

		if(offset > 0 && offset < max_overhang)
		{
                    if(dsp->strands[0]!=Seq_strand_minus)
			extend_len = MIN(offset, first_len - dsp->starts[2*endseg] - dsp->lens[endseg]);
                    else
			extend_len = MIN(offset, dsp->starts[2*endseg]);
			if(extend_len > 0)
			{
				dsp->lens[endseg] += extend_len;
                                if(dsp->strands[0]==Seq_strand_minus)
                                    dsp->starts[2*endseg]-=extend_len;
			}
		}
	}

}
	


/* Tries to Make a global alignment from a Set of local SeqAligns. 
 */
NLM_EXTERN SeqAlignPtr SeqAlignSetGlobalFromLocal(SeqAlignPtr align,SeqLocPtr loc_1, SeqLocPtr loc_2, BioseqPtr bsp_1, BioseqPtr bsp_2, FILE *err_fp,Int4 MaxGap)
{
	Int4 len_1, len_2;
	Char label[101];
        if(align==NULL || loc_1 == NULL || loc_2 == NULL || bsp_1 == NULL || bsp_2 == NULL)
            return NULL;
	len_1 = SeqLocLen(loc_1);
	if(len_1 <=0)
	{
            if(err_fp) {
		MuskSeqIdWrite(SeqLocId(loc_1), label, 100, PRINTID_TEXTID_ACCESSION, TRUE, TRUE);
		fprintf(err_fp, "bad length %ld in sequence %s\n", (long)len_1, label);
            }
            return NULL;
	}
	len_2 = SeqLocLen(loc_2);
	if(len_2 <=0)
	{
            if(err_fp) {
		MuskSeqIdWrite(SeqLocId(loc_2), label, 100, PRINTID_TEXTID_ACCESSION, TRUE, TRUE);
		fprintf(err_fp, "bad length %ld in sequence %s\n", (long) len_1, label);
            }
            return NULL;
	}

	if(align != NULL)
	{

	    align = SeqAlignConsistentDiagFilter(align,loc_1, loc_2,err_fp,MaxGap);
	}

	if(align != NULL)
	{
		extend_dsp_ends(bsp_1->length, bsp_2->length, align, MAX_EXTEND_OVERHANG);
		SeqAlignSwapOrder(align);
		extend_dsp_ends(bsp_2->length, bsp_1->length, align, MAX_EXTEND_OVERHANG);
		SeqAlignSwapOrder(align);
	}
	return align;
}

/* Attempt to Make a Global Alignment Using  local Alignments tools 
 */
NLM_EXTERN SeqAlignPtr compute_alignment(SeqLocPtr loc_1, SeqLocPtr loc_2, BioseqPtr bsp_1, BioseqPtr bsp_2, FILE *err_fp, BLAST_OptionsBlkPtr options, Boolean ck_polyA)
{
	SeqAlignPtr align = NULL;
        Char label[101];
	Int4 len_1, len_2;

	if(ck_polyA)
	    SeqLocTrimPolyATail(loc_2, bsp_2);

	len_1 = SeqLocLen(loc_1);
	len_2 = SeqLocLen(loc_2);

	if(len_1 <=0)
	{
            if(err_fp) {
		MuskSeqIdWrite(SeqLocId(loc_1), label, 100, PRINTID_TEXTID_ACCESSION, TRUE, TRUE);
		fprintf(err_fp, "bad length %ld in sequence %s\n", (long)len_1, label);
            }
            return NULL;
	}

	if(len_2 <=0)
	{
            if(err_fp) {
		MuskSeqIdWrite(SeqLocId(loc_2), label, 100, PRINTID_TEXTID_ACCESSION, TRUE, TRUE);
		fprintf(err_fp, "bad length %ld in sequence %s\n", (long)len_1, label);
            }
            return NULL;
	}

	if(len_1 <10000 && len_2 < 10000) 
	    align = SIM4ALN_choice(loc_1, loc_2, 1000, 8);

	if(align != NULL)
	{
		if(is_bad_align (bsp_1, bsp_2, align, 1000, MIN(len_1, len_2)) || 
			need_recompute(align, len_1, len_2))
			align = SeqAlignFree(align);
	}
	
	if(align == NULL)
	{
		/*align by splitting the FIRST sequence if it is too small*/
		align = SplitBlastTwoSeq(loc_2, loc_1,  SPLIT_LEN, 
			OVERLAP_LEN, options, err_fp);
		if(align != NULL)
			SeqAlignSwapOrder(align);
	}

	if(align != NULL)
	{
		extend_dsp_ends(bsp_1->length, bsp_2->length, align, MAX_EXTEND_OVERHANG);
		SeqAlignSwapOrder(align);
		extend_dsp_ends(bsp_2->length, bsp_1->length, align, MAX_EXTEND_OVERHANG);
		SeqAlignSwapOrder(align);
	}
	return align;
}

/***************************************************************************************
***
***  ValidateSeqAlignandACC
***	calls ValidateSeqAlign (in api directory)
***	and tests for occurrence of ACC string in sequence ID.
***	ACC|ACC# will be compared with the corresponding sequence (ACC#)
***	in the database and replaced by a far pointer if the sequences
***	are identical.
***
***************************************************************************************/
typedef struct saval {
  Boolean     message;
  Boolean     msg_success;
  Boolean     find_remote_bsp;
  Boolean     find_acc_bsp;
  Boolean     delete_salp;
  Boolean     delete_bsp;
  Boolean     retdel;
  ValNodePtr  ids;
  Uint2       entityID;
  Boolean     dirty;
} SaVal, PNTR SaValPtr;

static ValNodePtr errorp = NULL;

/******************************************************************
Output error message according to code defined in alignval.h.  
id refers to seqid of the sequence that causes the error 
and idcontext refers to other sequences in the same segment.  
Intvalue is used to indicate 1) the segment where the sequence 
with error is, or 2) the segtype in case of segtype error.  
Please note that not all errors report all three 
parameters(id, idcontext, Intvalue)
******************************************************************/ 
static void ValMessage (Int1 MessageCode, ErrSev errlevel, SeqIdPtr id, SeqIdPtr idcontext , Int4 Intvalue) 
{
  
  Char     buf[256], 
           buf3[64],
           string1[64],
           string2[252];

  string1[0] = '\0';
  string2[0] = '\0';
  SeqIdWrite(id, buf, PRINTID_FASTA_LONG, sizeof(buf)-1);
  switch(MessageCode)
  {
    case Err_SeqId:
      sprintf(string1, "SeqId");
      sprintf(string2, "Invalid Seq_id: %s\n", buf);
      break;

    case Err_Strand_Rev:      
      SeqIdWrite (idcontext, buf3, PRINTID_REPORT, sizeof (buf3));
      sprintf(string1, "Strand");
      sprintf(string2, "Alignment strand is reversed in segment %d for Seq ID: %s in the context of%s\n", Intvalue, buf, buf3);
      break;

    case Err_Denseg_Len_Start:
      SeqIdWrite (idcontext, buf3, PRINTID_REPORT, sizeof (buf3));
      sprintf(string1, "Start/Length");
      sprintf(string2, "Error in length and/or starts in segment %d for sequence ID: %s in the context of %s\n", Intvalue, buf, buf3);
      break;

    case  Err_Start_Less_Than_Zero:
      SeqIdWrite (idcontext, buf3, PRINTID_REPORT, sizeof (buf3));
      sprintf(string1, "Start");
      sprintf(string2, "Start point is less than zero in segment %d for sequence ID: %s in the context of %s\n", Intvalue, buf, buf3);
      break;

    case Err_Start_More_Than_Biolen:      
      SeqIdWrite (idcontext, buf3, PRINTID_REPORT, sizeof (buf3));
      sprintf(string1, "Start");
      sprintf(string2, "Start point is greater than total bioseq length in segment %d for sequence ID: %s in the context of%s\n", Intvalue, buf, buf3);
      break;

    case Err_End_Less_Than_Zero:
      SeqIdWrite (idcontext, buf3, PRINTID_REPORT, sizeof (buf3));
      sprintf(string1, "Length");
      sprintf(string2, "End point is less than zero in segment %d for sequence ID: %s in the context of %s\n", Intvalue, buf, buf3);
      break;

    case Err_End_More_Than_Biolen:
      SeqIdWrite (idcontext, buf3, PRINTID_REPORT, sizeof (buf3));
      sprintf(string1, "Length");
      sprintf(string2, "End point is greater than total bioseq length in segment %d for sequence ID: %s in the context of%s\n", Intvalue, buf, buf3);
      break;

    case Err_Len_Less_Than_Zero:
      SeqIdWrite (idcontext, buf3, PRINTID_REPORT, sizeof (buf3));
      sprintf(string1, "Length");
      sprintf(string2, "Segment length is less than zero in segment %d for sequence ID: %s in the context of %s\n", Intvalue, buf, buf3); 
      break;

    case Err_Len_More_Than_Biolen:
      SeqIdWrite (idcontext, buf3, PRINTID_REPORT, sizeof (buf3));
      sprintf(string1, "Length");
      sprintf(string2, "Segment length is greater than total bioseq length in segment %d for sequence ID: %s in the context of %s\n", Intvalue, buf, buf3);
      break; 
 
    case Err_Sum_Len_Start:
      SeqIdWrite (idcontext, buf3, PRINTID_REPORT, sizeof (buf3));
      sprintf(string1, "Start");
      sprintf(string2, "Sum of start point and segment is greater than total bioseq length in segment %d  for sequence ID: %s in the context of %s\n",  Intvalue, buf, buf3); 
      break;

    case Err_SeqAlign_DimSeqId_Not_Match:
      SeqIdWrite (idcontext, buf3, PRINTID_REPORT, sizeof (buf3));
      sprintf(string1, "SeqId");
      sprintf(string2, "The number of SeqId does not match the dimensions for sequence ID's %s\n", buf3); 
      break;

    case Err_Segs_DimSeqId_Not_Match:
      SeqIdWrite (idcontext, buf3, PRINTID_REPORT, sizeof (buf3));
      sprintf(string1, "SeqId");
      sprintf(string2, "The number of SeqId does not match the dimensions in segment %d for  sequence ID's %s\n", Intvalue, buf3); 
      break;

    case Err_Fastalike:
      SeqIdWrite (idcontext, buf3, PRINTID_REPORT, sizeof (buf3));
      sprintf(string1, "Fasta");
      sprintf(string2, "This may be a fasta-like alignment for SeqId: %s in the context of %s\n", buf, buf3); 
      break;

    case Err_Null_Segs:
      sprintf(string1, "Segs");
      sprintf(string2, "This alignment contains a null segs\n");
      break;

    case Err_Segment_Gap:
      SeqIdWrite (idcontext, buf3, PRINTID_REPORT, sizeof (buf3));
      sprintf(string1, "Segs");
      sprintf(string2, "Segment %d is a gap for all sequence with the following ID's: %s\n", Intvalue, buf3); 
      break;

    case Err_Segs_Dim_One:
      SeqIdWrite (idcontext, buf3, PRINTID_REPORT, sizeof (buf3));
      sprintf(string1, "Segs");
      sprintf(string2, "There is only one dimension in segment %d for  sequence ID's %s\n", Intvalue, buf3); 
      break;

    case Err_SeqAlign_Dim_One:
      SeqIdWrite (idcontext, buf3, PRINTID_REPORT, sizeof (buf3));
      sprintf(string1, "Dim");
      sprintf(string2, "There is only one dimension for sequence ID's %s\n", buf3); 
      break;

    case Err_Segtype :
      sprintf(string1, "Segs");
      sprintf(string2, "This alignment has a undefined or unsupported Seqalign segtype %d\n", Intvalue);
      break;

    defaulf:
      break;
  }
  if (StringLen(string1) > 0)
     errorp = BlastConstructErrorMessage (string1, string2, errlevel, &errorp);
}

 
/********************************************************
***
*** NormalizeSeqAlignId
***   Checks local seqid . if a seqid string contains "acc"
***   seqid has a correspondant sequence in db.
***   This local seqid is replaced by its gi number. 
***
***   The local sequence is compared to the sequence from 
***   the db. If the local sequence is a region of the db sequence
***   the positions in the seqalign are updaded with the offset.
***
***          Thanks to Mark for this useful function! 
**********************************************************/
static Int4 getlengthforid (SeqIdPtr sip)
{
  BioseqPtr        bsp;
  Int4             lens=0;
 
  if (sip==NULL)
     return 0;
  bsp = BioseqLockById (sip);
  if (bsp != NULL) {
     lens = bsp->length;
     BioseqUnlock (bsp);
  }
  return lens;
}

static ValNodePtr nrSeqIdAdd (ValNodePtr vnp, SeqIdPtr sip)
{
  ValNodePtr vnptmp=NULL;
  SeqIdPtr   siptmp;

  if (vnp!=NULL) {
     for (vnptmp=vnp; vnptmp!=NULL; vnptmp=vnptmp->next) {
        siptmp=(SeqIdPtr)vnptmp->data.ptrvalue;
           if (SeqIdForSameBioseq(sip, siptmp))
              break;
     }
  }
  if (vnptmp==NULL)
     ValNodeAddPointer(&vnp, 0, sip);
  return vnp;
}

static SeqIdPtr SeqIdReplaceID (SeqIdPtr head, SeqIdPtr pre, SeqIdPtr sip, SeqIdPtr next)
{
  SeqIdPtr tmp;

  if (pre == NULL)
  {
     head = SeqIdDup(sip);
     head->next = next;
     return head;
  }
  tmp = pre->next;
  pre->next = NULL;
  tmp->next = NULL;
  SeqIdFree (tmp);
  pre->next = SeqIdDup(sip);
  pre->next->next = next;
  return head;
}

static SeqAlignPtr LIBCALL SeqAlignBestHit (SeqAlignPtr salp, Int4 length, Int4 threshold)
{
  SeqAlignPtr  tmp, 
               ret=NULL;
  Int4         len;
  DenseDiagPtr ddp, 
               ddptmp;


  if (salp->segtype==2) {
     for (tmp=salp; tmp!=NULL; tmp=tmp->next)
     {
        len = SeqAlignLength(salp);
        if (100*((float)len/(float)length) >= threshold)
           break;
     }
     if (tmp!=NULL) {
        ret = (SeqAlignPtr) AsnIoMemCopy (tmp, (AsnReadFunc) SeqAlignAsnRead, (AsnWriteFunc) SeqAlignAsnWrite);  
     }
  }
  else if (salp->segtype ==1) {
     for (tmp=salp; tmp!=NULL; tmp=tmp->next)
     {
        for (ddp=salp->segs; ddp!=NULL; ddp=ddp->next) {
           if (100*((float)ddp->len/(float)length) >= threshold)
              break; 
        }
        if (ddp) {
           ddptmp = (DenseDiagPtr) AsnIoMemCopy (ddp,  (AsnReadFunc) DenseDiagAsnRead, (AsnWriteFunc) DenseDiagAsnWrite);
           if (ddptmp) {
              ret = SeqAlignNew();
              ret->segtype=1;
              ret->segs=ddptmp;
           }
           break;
        }
     }
  }
  return ret;
}

static void SeqAlignStartUpdate (SeqAlignPtr salp, SeqIdPtr target_sip, Int4 offset, Uint1 strand)
{
  SeqAlignPtr salptmp;
  DenseSegPtr dsp;
  SeqIdPtr    pre, sip, next;
  Int4Ptr     lenp,
              startp;
  Uint1Ptr    strandp;
  Int4        len_sum,
              start;
  Int2        index, k, j;

  if (salp==NULL || offset<=0)
     return;
  for (salptmp=salp; salptmp!=NULL; salptmp=salptmp->next)
  {
     if (salptmp->segtype == 2)
     {
        dsp = (DenseSegPtr) salptmp->segs;
        pre = NULL;
        index=0;
        sip=dsp->ids;
        while (sip)
        {
           next=sip->next;
           if (SeqIdForSameBioseq(target_sip, sip))
           {
              if (strand == Seq_strand_minus)
              {
                 strandp=dsp->strands;
                 strandp+=index;
                 for (j=0; j < dsp->numseg && strandp!=NULL; j++)
                 {
                    if (*strandp == Seq_strand_minus)
                       *strandp = Seq_strand_plus;
                    else if (*strandp == Seq_strand_plus)
                       *strandp = Seq_strand_minus;
                    strandp+=dsp->dim;
                 }
                 lenp=dsp->lens;
                 startp=dsp->starts;
                 start=startp[index];
                 len_sum=0;
                 j=dsp->dim*dsp->numseg-dsp->dim;
                 k=dsp->numseg-1;
                 for (; j>=0; j-=dsp->dim) {
                    startp[j+index]=start+len_sum;
                    len_sum+=lenp[k];
                    k--;
                 }
              }
              for (j=0; j<dsp->numseg; j++) {
                 if (dsp->starts[dsp->dim*j+index] != -1)
                    dsp->starts[dsp->dim*j+index] += offset;
              }
           }
           pre=sip;
           sip=next;
           index++;
        }
     }
  }
}
static void showtextalign_fromalign (SeqAlignPtr salp, CharPtr path, FILE *fp)
{
  SeqAnnotPtr sap;
  Int4        line = 80;     
  Uint4       option = 0;
  Boolean     do_close = TRUE;

  if (salp == NULL || (path==NULL && fp==NULL))
     return;
  if (path != NULL && fp == NULL) {
     fp = FileOpen (path, "a");
  }
  else
     do_close = FALSE;
  if (fp != NULL) {
/********
{{
     option = RULER_TOP;
     option|=DISPE_SHOWBLOCK;
     option|=VIEW_VARIA;
     option|=DISP_FULL_TXT;
     DDV_DisplayDefaultAlign (salp, 0, -1,-1, option, NULL, fp); 
}}
*******/
     sap = SeqAnnotNew ();
     sap->type = 2;
     sap->data = (Pointer) salp;
     option += TXALIGN_MASTER;
     option += TXALIGN_MISMATCH;
     ShowTextAlignFromAnnot (sap, line, fp, NULL, NULL, option, NULL, NULL, NULL);
     sap->data = NULL;
     SeqAnnotFree (sap);
     if (do_close) {
        FileClose(fp);
     }
  } 
}

static ValNodePtr CCNormalizeSeqAlignId (SeqAlignPtr salp, ValNodePtr vnp)
{
  BLAST_OptionsBlkPtr options;
  SeqLocPtr           slp1, slp2;
  MsgAnswer           ans;
  DenseSegPtr         dsp;
  SeqIdPtr            sip,
                      dbsip = NULL,
                      lclsip,
                      presip, 
                      next;
  SeqAlignPtr         seqalign = NULL;
  SeqAlignPtr         bestsalp;
  CharPtr             TmpBuff, tmp;
  Char                str [52];
  Int4                gi = 0,
                      offset,
                      totlenlcl, totlendb;
  Int4                j, k;
  Int2                index;
  Uint1               strand;
  Boolean             ok, 
                      found;
  
  Char                strLog[50];

  if (salp!=NULL) {
     if (salp->segtype == 2) {
        dsp = (DenseSegPtr) salp->segs;
        presip = NULL;
        sip = dsp->ids;
        index = 0;
        found = FALSE;
        while (sip != NULL) 
        {
           next = sip->next;
           lclsip = SeqIdDup (sip);
           SeqIdWrite (lclsip, str, PRINTID_FASTA_LONG, 50);
           tmp = StringStr (str, "acc");
           if (tmp==NULL) 
           {
              tmp = StringStr (str, "ACC");
           }
           if (tmp!=NULL) {
              tmp++; tmp++; tmp++;
              if (*tmp == '|')
                 tmp++;   
              TmpBuff = tmp;
              while (*tmp!='\0' && *tmp != '|' && *tmp!='\n')
                 tmp++;
              *tmp = '\0';

              ok = FALSE;
              j = StringLen (TmpBuff);
              for(k =0; k < j; k++) {
                 if(!isdigit(TmpBuff[k])) {
                    break;
                 }
              }
              dbsip=NULL;
              if(k != j) {
                 ok=(IS_ntdb_accession(TmpBuff) || IS_protdb_accession(TmpBuff));
                 if (ok) {
                    dbsip = SeqIdFromAccession (TmpBuff, 0, NULL);
                 }
              }
              else {
                 gi = (Int4)atol(TmpBuff);
                 if (gi>0) {
                    dbsip = ValNodeNew (NULL);
                    if (dbsip) {
                       dbsip->choice = SEQID_GI;
                       dbsip->data.intvalue = (Int4)gi;
                    }
                 }
              }
              if (dbsip!=NULL) {
                 totlendb = getlengthforid(dbsip); 
                 totlenlcl = getlengthforid(lclsip);
                  
                 slp1 = SeqLocIntNew (0, totlenlcl-1, Seq_strand_both, lclsip);
                 slp2 = SeqLocIntNew (0, totlendb-1, Seq_strand_both, dbsip);
                 options = BLASTOptionNew("blastn", FALSE);
/*
                 options->penalty = -5; 
                 options->cutoff_s = totlenlcl; 
*/
                 seqalign = BlastTwoSequencesByLoc (slp1, slp2, NULL, options);
                 
                 bestsalp = SeqAlignBestHit (seqalign, totlenlcl, 100);
                 if (bestsalp) 
                 {
                    SeqIdWrite (dbsip, strLog, PRINTID_TEXTID_ACCESSION, 50);
                    ans = Message (MSG_OKC, "This alignment contains \"%s\" that is already in GenBank. \n Do you wish to replace it?", strLog);
                    if (ans != ANS_CANCEL) 
                    {
                       offset = SeqAlignStart(bestsalp, 1)-SeqAlignStart(bestsalp, 0);
                       if (SeqAlignStrand(bestsalp, 0)==Seq_strand_minus || SeqAlignStrand(bestsalp, 1)==Seq_strand_minus)
                          strand=Seq_strand_minus;
                       else
                          strand=Seq_strand_plus;
                       SeqAlignStartUpdate (salp, lclsip, offset, strand);
                       dsp->ids = SeqIdReplaceID(dsp->ids, presip, dbsip, next); 
                       if (presip)
                          sip = presip->next;
                       else
                          sip = dsp->ids;
                       SeqAlignReplaceId (lclsip, dbsip, salp);
                       vnp = nrSeqIdAdd (vnp, lclsip);
                       found = TRUE;
                    }
                    SeqAlignFree (bestsalp);
                 }
                 else {
                    SeqIdWrite (dbsip, strLog, PRINTID_TEXTID_ACCESSION, 50);
                    ans = Message (MSG_OKC, "This alignment contains \"%s\" that is already in GenBank.\n However, the local version is not identical to the database version.\n Do you wish to replace it anyway ?\n If you cancel, the alignment of the local and the database versions \nof \"%s\" will be saved in the error file \"error.log\"", strLog, strLog);
                    if (ans != ANS_CANCEL) 
                    {
                       bestsalp = seqalign;
                       offset = SeqAlignStart(bestsalp, 1)-SeqAlignStart(bestsalp, 0);
                       if (SeqAlignStrand(bestsalp, 0)==Seq_strand_minus || SeqAlignStrand(bestsalp, 1)==Seq_strand_minus)
                          strand=Seq_strand_minus;
                       else
                          strand=Seq_strand_plus;
                       SeqAlignStartUpdate (salp, lclsip, offset, strand);
                       dsp->ids = SeqIdReplaceID(dsp->ids, presip, dbsip, next); 
                       if (presip)
                          sip = presip->next;
                       else
                          sip = dsp->ids;
                       SeqAlignReplaceId (lclsip, dbsip, salp);
                       vnp = nrSeqIdAdd (vnp, lclsip);
                       found = TRUE;
                    }
                    else {
                       showtextalign_fromalign (seqalign, "error.log", NULL);
                    }
                    sip->next = next;
                 }
                 if (seqalign)
                    seqalign = SeqAlignFree (seqalign);
              } 
              else {
                 SeqIdWrite (sip, strLog, PRINTID_TEXTID_ACCESSION, 50);
                 ans = Message (MSG_OK, "This alignment contains \"%s\" that can not be found in GenBank.\nPlease check the accession number.\n", strLog);
                 sip->next = next;
              }
           }
           presip = sip;
           sip = next;
           index++;
           found = FALSE;
        }
     }
  }
  return vnp;  
}


static Boolean check_dbid_seqalign (SeqAlignPtr salp)
{
  DenseSegPtr dsp;
  SeqIdPtr    sip, next;
  Char        str [52];
  CharPtr     TmpBuff, tmp;
  Int4        j, k;
  Boolean     found = FALSE;

  if (salp!=NULL) 
  {
     if (salp->segtype == 2) 
     {
        dsp = (DenseSegPtr) salp->segs;
        sip = dsp->ids;
        while (!found && sip != NULL) 
        {
           next = sip->next;
           SeqIdWrite (sip, str, PRINTID_FASTA_LONG, 50);
           tmp = StringStr (str, "acc");
           if (tmp!=NULL) {
              tmp++; tmp++; tmp++;
              if (*tmp == '|')
                 tmp++;
              TmpBuff = tmp;
              while (*tmp!='\0' && *tmp != '|' && *tmp!='\n')
                 tmp++;
              *tmp = '\0';

              j = StringLen (TmpBuff);
              for(k =0; k < j; k++) {
                 if(!isdigit(TmpBuff[k])) {
                    break;
                 }
              }
              if(k != j) {
                found=(IS_ntdb_accession(TmpBuff) || IS_protdb_accession(TmpBuff));
              }
           }  
           sip = next;
        }     
     }
  }     
  return found;
}


/******************************************************************
validate each alignment sequentially.  
This function will subject the seqalign to all validation functions
******************************************************************/ 
/*********************************************************/
static void delete_bioseqs (ValNodePtr ids, Uint2 entityID)
{
  SeqEntryPtr  sep_top;
  SeqEntryPtr  sep_del;
  ValNodePtr   vnp;
  SeqIdPtr     sip;
  SeqLocPtr    slp;
  BioseqPtr    bsp;
  ObjMgrDataPtr  omdptop;
  ObjMgrData     omdata;
  Uint2          parenttype;
  Pointer        parentptr;

  if (ids == NULL)
     return;
  sep_top = GetTopSeqEntryForEntityID (entityID);
  SaveSeqEntryObjMgrData (sep_top, &omdptop, &omdata);
  GetSeqEntryParent (sep_top, &parentptr, &parenttype);

  vnp=ids;
  while (vnp!=NULL)
  {
     sip = (SeqIdPtr) vnp->data.ptrvalue;
     if (sip!=NULL) {
        slp = (SeqLocPtr)ValNodeNew (NULL);
        slp->choice = SEQLOC_WHOLE;
        slp->data.ptrvalue = sip;
        bsp = GetBioseqGivenSeqLoc (slp, entityID);
        if (bsp!=NULL) {
           sep_del=GetBestTopParentForData (entityID, bsp);
           RemoveSeqEntryFromSeqEntry (sep_top, sep_del, FALSE);
        }
        slp->data.ptrvalue = NULL;
        SeqLocFree (slp);
     }
     vnp=vnp->next;
  }
  SeqMgrLinkSeqEntry (sep_top, parenttype, parentptr);
  RestoreSeqEntryObjMgrData (sep_top, omdptop, &omdata);
  RenormalizeNucProtSets (sep_top, TRUE);

  for (vnp=ids; vnp!=NULL; vnp=vnp->next) {
     SeqIdFree ((SeqIdPtr) vnp->data.ptrvalue);
     vnp->data.ptrvalue = NULL;
  }
  ValNodeFree (vnp);
  return;
}

NLM_EXTERN Boolean ValidateSeqAlignandACC (SeqAlignPtr salp, Uint2 entityID, Boolean message,
                         Boolean msg_success, Boolean find_remote_bsp,Boolean find_acc_bsp,
                         Boolean delete_bsp, Boolean delete_salp, BoolPtr dirty)
{  
  SeqAlignPtr  pre,
               salptmp;
  SaVal        sv;
  SaValPtr     svp;
  MsgAnswer    ans;
  Int2         err_count=0,
               salp_count=0;
  Boolean      ok; 

  if(salp!=NULL)
  {
        sv.message = message;
        sv.msg_success = msg_success;
        sv.find_remote_bsp = find_remote_bsp;
        sv.find_acc_bsp = find_acc_bsp;
        sv.delete_salp = delete_salp;
        sv.delete_bsp = delete_bsp;
        sv.retdel = TRUE;
        sv.ids = NULL;
        sv.entityID = entityID; 
        sv.dirty = FALSE;   
        svp = &sv;   
     pre=NULL;
     salptmp=salp; 
     while (salptmp)
     {
        salp_count++;
        if(salp->segtype==5)
        {
           ValidateSeqAlignandACC ((SeqAlignPtr) (salptmp->segs), entityID, message, msg_success, find_remote_bsp, find_acc_bsp, delete_bsp, delete_salp, &svp->dirty);
        } 
        else if (salp->segtype<1 || salp->segtype>4)
        {
           ValMessage (Err_Segtype, SEV_ERROR, NULL, NULL, salptmp->segtype);
        }
        else {
           ValidateSeqAlign (salptmp, svp->entityID, svp->message, svp->msg_success, svp->find_remote_bsp, svp->delete_bsp, svp->delete_salp, &svp->dirty);
           if (svp->find_acc_bsp) {
	      ok = check_dbid_seqalign (salptmp);
	      if (ok) {
                 svp->ids = CCNormalizeSeqAlignId (salptmp, svp->ids);
                 if (svp->ids!=NULL && svp->entityID > 0) {
                    if (svp->delete_bsp)
                       delete_bioseqs (svp->ids, svp->entityID); 
                    svp->dirty = TRUE;
                 }
              }       	
           }
        }     	
       	if (errorp)
       	{
       	   if(svp->message)
       	   {
              BlastErrorPrint (errorp);
       	      errorp = BlastErrorChainDestroy (errorp);
       	   }
       	   if (svp->delete_salp)
       	   {
            if (pre==NULL) {
              salp=salptmp->next;
              salptmp->next = NULL;
              SeqAlignFree (salptmp);
              salptmp = salp;
            }
            else {
              pre->next = salptmp->next;
              salptmp->next = NULL;
              SeqAlignFree (salptmp);
              salptmp = pre->next;
            }
           }
       	   else {
       	      salptmp = salptmp->next;
       	   }
       	   err_count++;
           svp->retdel=FALSE;
        }
       	else {
       	   salptmp = salptmp->next;
       	}
     }
     if (err_count==0 && svp->msg_success) {
        if (salp_count>1)
           ans = Message (MSG_OK, "Validation test of %d alignments succeded", salp_count);
        else
           ans = Message (MSG_OK, "Validation test of the alignment succeded");
     }
     if (dirty)
        *dirty = svp->dirty;
  }   
  return svp->retdel;
} 


/******************************************************************
call back function for REGISTER_ALIGNVALIDATION defined in sequin4.c.  
Starting point for seqalign validation if user clicked on 
SeqalignValidation under menu Filer/Alignment.  
Either individual alignment or alignment block 
should be highlighted for this validation to work
******************************************************************/ 

static SeqAlignPtr LIBCALL is_salp_in_sap (SeqAnnotPtr sap, Uint1 choice)
{
  SeqAlignPtr      salp = NULL;

  if (sap != NULL) {
     for (; sap!= NULL; sap=sap->next) {
        if (sap->type == choice) {
           salp = (SeqAlignPtr) sap->data;
           return salp;
        }
     }   
  }
  return NULL;
}

static Pointer LIBCALL sap_empty (SeqAnnotPtr sap, Uint1 type, Pointer PNTR ptr)
{
  SeqAlignPtr      salp = NULL;

  if (sap != NULL) {
     for (; sap!= NULL; sap=sap->next) {
        if (sap->type == type) {
           salp = (SeqAlignPtr) sap->data;
           if (ptr!=NULL)
              *ptr = (Pointer) sap;
           break;
        }
     }
  }
  return salp;
}

NLM_EXTERN Int2 LIBCALLBACK ValidateSeqAlignandACCFromData (Pointer data)
{ 
 
  OMProcControlPtr  ompcp;
  SeqAlignPtr       salp=NULL;
  SeqAnnotPtr       sap=NULL;
  SeqEntryPtr       sep=NULL;
  
  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL || ompcp->proc == NULL) return OM_MSG_RET_ERROR;
  
  if (ompcp->input_data == NULL) return OM_MSG_RET_ERROR;
  
  switch(ompcp->input_itemtype)
    {
    case OBJ_BIOSEQ :
      sep = SeqMgrGetSeqEntryForData (ompcp->input_data);
      break;
    case OBJ_BIOSEQSET :
      sep = SeqMgrGetSeqEntryForData (ompcp->input_data);
      break;
      /*if clicked on alignment block*/
    case OBJ_SEQANNOT:
      sap=(SeqAnnotPtr) (ompcp->input_data);
      break;
      /*if clicked on individual alignment*/
    case OBJ_SEQALIGN:
      salp=(SeqAlignPtr) (ompcp->input_data);
      break;
    case 0 :
      return OM_MSG_RET_ERROR;
    default :
      return OM_MSG_RET_ERROR;
  }
  
  ErrSetMessageLevel(SEV_ERROR);
  if(sap!=NULL)
  {
     salp=is_salp_in_sap(sap, 2);
     ValidateSeqAlignandACC (salp, 0, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, NULL);
  }
  if (salp!=NULL) {
     ValidateSeqAlignandACC (salp, 0, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, NULL);
  }
  if (sep!=NULL) {
     ValidateSeqAlignandACCInSeqEntry (sep, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE);
  }
  return OM_MSG_RET_DONE;
}

static void ValidateSeqAlignandACCCallback (SeqEntryPtr sep, Pointer mydata,
                                          Int4 index, Int2 indent)
{
  BioseqPtr          bsp;
  BioseqSetPtr       bssp;
  SeqAlignPtr        salp;
  SaValPtr           svp;

  if (sep != NULL && sep->data.ptrvalue && mydata != NULL) {
     svp = (SaValPtr)mydata;
     if (IS_Bioseq(sep)) {
        bsp = (BioseqPtr) sep->data.ptrvalue;
        if (bsp!=NULL) {
           salp=sap_empty(bsp->annot, 2, NULL);
           if (salp!=NULL) {
              ValidateSeqAlignandACC (salp, svp->entityID, svp->message, svp->msg_success, svp->find_remote_bsp, svp->find_acc_bsp, svp->delete_bsp, svp->delete_salp, &svp->dirty);
           }
        }
     }   
     else if(IS_Bioseq_set(sep)) {
        bssp = (BioseqSetPtr)sep->data.ptrvalue;
        if (bssp!=NULL) {
           salp=sap_empty(bssp->annot, 2, NULL);
           if (salp!=NULL) {
              ValidateSeqAlignandACC (salp, svp->entityID, svp->message, svp->msg_success, svp->find_remote_bsp, svp->find_acc_bsp, svp->delete_bsp, svp->delete_salp, &svp->dirty);
           }
        }
     }
  }
}



NLM_EXTERN Boolean ValidateSeqAlignandACCInSeqEntry (SeqEntryPtr sep, Boolean message, 
                                 Boolean msg_success, Boolean find_remote_bsp, Boolean find_acc_bsp,
                                 Boolean delete_bsp, Boolean delete_salp)
{
  SeqEntryPtr      sep_head;
  Uint2            entityID;
  SaVal            sv;
  Boolean          success=TRUE;

  entityID = ObjMgrGetEntityIDForChoice (sep);
  if (entityID > 0) {
     sep_head = GetTopSeqEntryForEntityID (entityID);
     if (sep_head != NULL) {
        sv.message = message;
        sv.msg_success = msg_success;
        sv.find_remote_bsp = find_remote_bsp;
        sv.find_acc_bsp = find_acc_bsp;
        sv.delete_salp = delete_salp;
        sv.delete_bsp = delete_bsp;
        sv.retdel = TRUE;
        sv.ids = NULL;
        sv.entityID = entityID; 
        sv.dirty = FALSE;
        SeqEntryExplore (sep_head, (Pointer)&sv, ValidateSeqAlignandACCCallback);
        if (sv.dirty) {
           ObjMgrSetDirtyFlag (entityID, TRUE);
           ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
        }
        success = sv.retdel;
     }
  }
  return success;
}