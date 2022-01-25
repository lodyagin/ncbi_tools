#include <blast.h>
#include <jzcoll.h>
#include <jzmisc.h>
#include <simutil.h>
#include <salpedit.h>

#define MIN_SEG_SIZE 10
#define MIN_DIAG_LEN	20	/*minimal amount of overlap required*/

#define SIM2_CUTOFF 20.0

/* BLast can handle longer queries more efficiently.. though
 there is still a hardcoded limit to the number of HSP's */
#define SPLIT_LEN       10000 
#define OVERLAP_LEN     500

#define PERCENT_A_RES	90


#define MAX_HANG 100

#define MAX_EXTEND_OVERHANG	50
static SeqAlignPtr SeqAlignInsert(SeqIdPtr from_id, Int4 from, Int4 to, Uint1 strand, SeqIdPtr to_id, Int4 pos, SeqAlignPtr head, Boolean merge, Boolean add_insertion);
static Int2 get_t_order(SeqIdPtr ids, SeqIdPtr t_id);
static Boolean get_seg(ValNodePtr v_starts, Int2 dim, Int4Ptr starts);
static ValNodePtr Free_node(ValNodePtr v_starts, ValNodePtr pv_starts, Int2 dim);
static Boolean is_gap_segment(Int4Ptr c_starts, Int2 dim);
static void mod_prev_start(ValNodePtr p_starts, Int4Ptr c_starts, Int4Ptr c_strands, Int2 dim);
static Boolean mod_alignment(ValNodePtr v_starts, ValNodePtr v_strands, ValNodePtr v_lens, Int2 dim, Int2 m_order, Int2Ptr n_segs);
static void load_new_data(Int4Ptr starts, Uint1Ptr strands, Int4Ptr lens, ValNodePtr v_starts, ValNodePtr v_strands, ValNodePtr v_lens, Int2 dim, Int2 cur_seg);
static SeqIdPtr get_nth_id(SeqIdPtr ids, Int2 order);
static void count_seg_num(ValNodePtr vnp, Int2 dim);
static DenseSegPtr make_dsp(ValNodePtr v_starts, ValNodePtr v_lens, ValNodePtr v_strands, Int2 m_order, Int2 dim, DenseSegPtr prev, Int2 cur_seg, Boolean merge, SeqIdPtr m_sip, SeqIdPtr s_sip);
static Boolean DenseSegInsert(SeqIdPtr from_id, Int4 from, Int4 to, Uint1 strand, SeqIdPtr to_id, Int4 pos, DenseSegPtr dsp, Boolean merge, BoolPtr seq_insert);
static ValNodePtr store_data(Int4 val, ValNodePtr head);
static ValNodePtr store_reverse_data(Int4 val, ValNodePtr head);
static ValNodePtr find_last_node(ValNodePtr head);
static ValNodePtr link_to_end(ValNodePtr new_vnp, ValNodePtr head);
static ValNodePtr find_nth_node(ValNodePtr vnp, Int2 order);
static Int4 get_last_pos(ValNodePtr v_starts, ValNodePtr v_lens);
static Boolean merge_two_seg(ValNodePtr PNTR v_starts, ValNodePtr vs_starts, ValNodePtr PNTR v_lens, ValNodePtr vs_lens, ValNodePtr PNTR v_strands, ValNodePtr vs_strands, Int4 pos, Uint1 s_strand, Int2 numseg);
static Boolean dsp_process(DenseSegPtr f_dsp, Int4 from, Int4 to, Uint1 strand, SeqIdPtr from_id, SeqIdPtr to_id, Int4 pos, ValNodePtr PNTR vs_starts, ValNodePtr PNTR vs_strands, ValNodePtr PNTR vs_lens);
static Int2 get_cur_seg(ValNodePtr lens);
static Boolean DenseSegMerge(DenseSegPtr f_dsp, SeqIdPtr from_id, Int4 from, Int4 to, Uint1 strand, SeqIdPtr to_id, Int4 pos, DenseSegPtr t_dsp);
static SeqAlignPtr SeqAlignLink(SeqAlignPtr new_salp, SeqAlignPtr align);
static SeqAlignPtr make_new_align(DenseSegPtr f_dsp, Int4 from, Int4 to, Uint1 strand, SeqIdPtr from_id, SeqIdPtr to_id, Int4 pos, SeqAlignPtr head);
static SeqAlignPtr JZSeqAlignMerge(SeqAlignPtr f_align, SeqIdPtr from_id, Int4 from, Int4 to, Uint1 strand, SeqIdPtr to_id, Int4 pos, SeqAlignPtr t_align);
static SeqAlignPtr AttachSeqInAlign(SeqIdPtr from_id, Int4 from, Int4 to, Uint1 strand, SeqIdPtr to_id, Int4 pos, SeqAlignPtr head);
static void mod_old_align(SeqAlignPtr PNTR old, SeqIdPtr sip);
static void attach_new_align(SeqAlignPtr new_salp, SeqAlignPtr PNTR old);
static void delete_bad_node(ValNodePtr PNTR head, Uint1 choice);
static Int4 count_coverage (BoolPtr used_list, Int4 size);

static int LIBCALLBACK CompareChainProc(VoidPtr ptr1, VoidPtr ptr2);
static Uint1 has_too_much_overlap (Int4 start_1, Int4 stop_1, Int4 start_2, Int4 stop_2);
static Boolean ck_delete_next_chain (DenseDiagPtr ddp, ValNodePtr PNTR ddp_list, Int2 s_order);
static void delete_redundant_node (ValNodePtr PNTR ddp_list, Int2 s_order);
static Int4Ptr convert_ddp_to_array (ValNodePtr ddp_list, Uint1 strand, Int2 s_order, Int4Ptr p_num);
static Boolean clean_up_discrepancies (Int4Ptr val, Int4 size, Boolean descend, Int4Ptr score_order);
static void clean_genomic_order (ValNodePtr PNTR head, Int2 s_order, Uint1 strand, Int4Ptr score_order);
static Int4 get_score_from_list (ScorePtr scores);
static Int4 get_score_from_DenseDiag (DenseDiagPtr ddp);
static Boolean filter_this_seg (DenseDiagPtr curr, DenseDiagPtr PNTR h_list, Uint1 order);
static Boolean is_same_orientation(Uint1 strand_1, Uint1 strand_2);
static void filter_wrong_orientation (SeqAlignPtr align, Boolean same_orient);
static FloatHi get_overlap(BoolPtr status, Int4 from, Int4 len, Int4Ptr poverlap);
static Boolean is_overlap(BoolPtr status, Int4 from, Int4 len);
static void set_ddp_status(ValNodePtr ddp_list, Int4 first_len, Int4 second_len);
static Boolean free_diag_in_list(ValNodePtr ddp_list, DenseDiagPtr ddp);
static void filter_repeats (ValNodePtr PNTR ddp_list, Int4 first_len, Int4 second_len);
static Int4Ptr get_process_list (ValNodePtr ddp_list);
static ValNodePtr get_nth_node(ValNodePtr vnp, Int4 order);
static Uint1 get_alignment_orientation(ValNodePtr PNTR pddp_list, BioseqPtr bsp, Int2 s_order);
static Boolean load_ddp_in_order(ValNodePtr PNTR h_list, DenseDiagPtr ddp, Boolean reverse, Boolean add);
static Int4Ptr load_ddp_scores(ValNodePtr ddp_list, Int4Ptr num);
static Boolean is_ddp_overlap(DenseDiagPtr ddp, BoolPtr first_status, BoolPtr second_status);
static Int4 lns(Int4Ptr S, Int4 len_s, Int4Ptr LNS);
static void sort_list_order(ValNodePtr PNTR ph_list, Uint1 strand);
static ValNodePtr check_ddp_order(ValNodePtr ddp_list, Uint1 strand, Int4 first_len, Int4 second_len);
static ValNodePtr FilterSeqAlign(SeqAlignPtr align, SeqIdPtr s_sip, Uint1Ptr p_strand);
static Int4 get_dsp_align_len(DenseSegPtr dsp);
static ScorePtr make_seg_score(Int4 align_len, Int4 seg_len, Int4 total_score);
static DenseDiagPtr covert_dsp_to_ddp(SeqAlignPtr align, DenseDiagPtr PNTR h_ddp);
static SeqAlignPtr ConvertDspAlignList(SeqAlignPtr align);
static SeqAlignPtr clean_up_align_list(SeqAlignPtr PNTR h_align, ValNodePtr head);
static Boolean check_align_match_diagnol(SeqAlignPtr align, ValNodePtr ddp_list, Uint1 strand);
static Boolean ModFilterDenseSegAlign(SeqAlignPtr PNTR align, SeqIdPtr sip);
static SeqAlignPtr FilterDenseSegAlign(SeqAlignPtr align, SeqIdPtr sip);
static Boolean FigureAlignRange(SeqAlignPtr align, SeqLocPtr m_loc, SeqLocPtr s_loc);
static Boolean reload_seq_loc(SeqLocPtr slp, Int4 from, Int4 to, Uint1 strand);
static Boolean get_end_seg(DenseSegPtr dsp, Int4 min_seg_len, SeqLocPtr loc_1, SeqLocPtr loc_2, Boolean head);
static Boolean is_loc_overlap(SeqLocPtr loc_1, SeqLocPtr loc_2, Int4Ptr min_val, Int4Ptr max_val, Int4 gap_len, Int4 max_align_size, Int4Ptr p_gap_val);
static Int4 find_matching_position(DenseSegPtr dsp, Int4 pos, Uint1 pos_order);
static void SeqAnnotWrite(SeqAlignPtr align);
static Boolean select_overlap_loc(SeqLocPtr loc_1_1, SeqLocPtr loc_2_1, SeqLocPtr loc_1_2, SeqLocPtr loc_2_2, Uint1Ptr p_order);
static SeqAlignPtr MergeTwoDspBySIM4(SeqAlignPtr align_1, SeqAlignPtr align_2);
static Int4 get_align_len_in_overlap (SeqAlignPtr align, Int4 from, Int4 to, Int2 order);
static int LIBCALLBACK AlignOverlapCompProc (VoidPtr ptr1, VoidPtr ptr2);
static ValNodePtr build_overlap_list(SeqAlignPtr align, Int4 from, Int4 to, Int2 order);
static SeqAnnotPtr open_annot(CharPtr f_name);
static int LIBCALLBACK SegOrderCompProc (VoidPtr ptr1, VoidPtr ptr2);
static ValNodePtr find_overlap_diagnol(SeqAlignPtr align, Int4 from, Int4 to, Int2 order);
static Boolean is_dsp_same(DenseSegPtr dsp_1, DenseSegPtr dsp_2);
static Int4 get_score_value(ScorePtr sp);
static void load_big_score (SeqAlignPtr align_1, SeqAlignPtr align_2);
static Boolean merge_two_align(SeqAlignPtr align_1, SeqAlignPtr align_2, Int4 from, Int4 to, Int2 order);
static SeqAlignPtr extract_align_from_list(SeqAlignPtr PNTR list, SeqAlignPtr align);
static void reverse_alignment(SeqAlignPtr align, Uint1 order);
static void link_this_align(SeqAlignPtr PNTR h_align, SeqAlignPtr align);
static Int4 get_align_length(SeqAlignPtr align);
static int LIBCALLBACK CompareAlignInfoProc(VoidPtr ptr1, VoidPtr ptr2);
static SeqAlignPtr sort_align_by_length(SeqAlignPtr align);
static Int4 get_align_mid_point(SeqAlignPtr align, Uint1 order);
static SeqAlignPtr sort_align_by_region(SeqAlignPtr align, Uint1 order);
static Int4 get_align_list_len(SeqAlignPtr align);
static SeqAlignPtr MergeToOneAlignment(SeqAlignPtr PNTR palign);
static SeqAlignPtr extract_align_with_sameID(SeqAlignPtr PNTR align, SeqIdPtr sip_1, SeqIdPtr sip_2);
static Boolean ModifyAlignList(SeqAlignPtr PNTR palign);
static Int2 get_polyA_index (Uint1Ptr res_buf, Int2 len);
static SeqLocPtr check_polyA_tail (BioseqPtr bsp);
static Int4 get_align_num(SeqAlignPtr align);
static void save_output (SeqIdPtr sip, CharPtr f_name, SeqAlignPtr align);
static SeqAlignPtr find_best_align(SeqAlignPtr align, Int4Ptr pmax_score);
static Boolean is_bad_align PROTO((BioseqPtr m_bsp, BioseqPtr s_bsp, SeqAlignPtr align, Int4 min_align_len,Int4 loc_len));
static Boolean is_bad_blast_alignment(SeqAlignPtr align, SeqLocPtr slp1, SeqLocPtr slp2, Int4 min_align_len);
static void re_order_alignment(SeqAlignPtr align);
static Boolean need_recompute(SeqAlignPtr align, Int4 m_len, Int4 s_len);
static void extend_dsp_ends(Int4 first_len, Int4 second_len, SeqAlignPtr align, Int4 max_overhang);
                             static void SeqLocTrimPolyATail(SeqLocPtr loc_2, BioseqPtr bsp_2);
static SeqAlignPtr SeqAlignSetGlobalFromLocal(SeqAlignPtr align,SeqLocPtr loc_1, SeqLocPtr loc_2, BioseqPtr bsp_1, BioseqPtr bsp_2, FILE *fp, Boolean record_err);
static SeqAlignPtr compute_alignment(SeqLocPtr loc_1, SeqLocPtr loc_2, BioseqPtr bsp_1, BioseqPtr bsp_2, FILE *fp, BLAST_OptionsBlkPtr options, Boolean ck_polyA, Boolean record_err);
static void reverse_seqloc (SeqLocPtr s_loc);
static BioseqPtr create_fake_bioseq(SeqLocPtr seq_loc, CharPtr fake_name);
static SeqLocPtr dup_mul_locs(SeqLocPtr slp);
static Boolean second_seq_has_gap(BioseqPtr s_bsp);
static SeqIdPtr get_best_id(BioseqPtr bsp);
static SeqAlignPtr process_multi_cds(BioseqPtr fake_gcontig_bsp, SeqLocPtr gcontig_loc, BioseqPtr s_bsp, BioseqPtr contig_bsp, Boolean is_unigene, FILE *err_fp, BLAST_OptionsBlkPtr options);
static SeqIdPtr find_gcontig_id(BioseqPtr bsp);
static Boolean is_virtual_segment(SeqIdPtr sip);
static SeqLocPtr build_slp_for_gcontig (BioseqPtr gcontig_bsp, SeqIdPtr gcontig_sip);
static SeqAlignPtr SplitBlastTwoSeq(SeqLocPtr slp1, SeqLocPtr slp2, 
		Int4 split_len, Int4 overlap_len, BLAST_OptionsBlkPtr options, Boolean store_err);

static SeqAlignPtr SeqAlignConsistentDiagFilter(SeqAlignPtr align,SeqLocPtr slp1, SeqLocPtr slp2, 
					 Boolean store_err);


/**************************************************************************
***
*	Functions for editing a Seq-align object
*
***************************************************************************
***/
static SeqAlignPtr SeqAlignInsert(SeqIdPtr from_id, Int4 from, Int4 to, Uint1 strand, SeqIdPtr to_id, Int4 pos, SeqAlignPtr head, Boolean merge, Boolean add_insertion)
{
	DenseSegPtr dsp;
        DenseDiagPtr ddp;
	StdSegPtr ssp;
	Boolean seq_insert;
	BioseqPtr bsp;
	SeqAlignPtr align;

	align = head;
	seq_insert = FALSE;
	while(align)
	{
		switch(align->segtype)
		{
			case 0:
				break;
			case 1:
				ddp = (DenseDiagPtr)(align->segs);
				break;
			case 2:
				dsp = (DenseSegPtr)(align->segs);
				DenseSegInsert(from_id, from, to, strand, to_id, pos, dsp, merge, &seq_insert);
				break;
		        case 3:
				ssp = (StdSegPtr)(align->segs);
				break;
			default:
				break;
		}
		align = align->next;
	}

	/*sequence is not inserted*/
	if(!seq_insert && add_insertion)
	{
		bsp = BioseqFind(from_id);
		if(bsp->repr != Seq_repr_seg)
			head= AttachSeqInAlign(from_id, from, to, strand, to_id, pos, head);
	}
	return head;

}

static Int2 get_t_order(SeqIdPtr ids, SeqIdPtr t_id)
{
	Int2 t_order = 0;

	if(t_id == NULL)
		return -1;
	while(ids)
	{
		if(SeqIdForSameBioseq(t_id, ids))
			return t_order;
		++t_order;
		ids = ids->next;
	}
	return -1;

}



static Boolean get_seg(ValNodePtr v_starts, Int2 dim, Int4Ptr starts)
{
	Int2 i;
	
	for(i=0; i<dim && v_starts != NULL; ++i, v_starts = v_starts->next)
		starts[i] = v_starts->data.intvalue;
	return (i == dim);
}


static ValNodePtr Free_node(ValNodePtr v_starts, ValNodePtr pv_starts, Int2 dim)
{
	Int2 i;
	ValNodePtr next = NULL;

	for(i =0; i<dim; ++i)
	{
		next= v_starts->next;
		v_starts->next = NULL;
		ValNodeFree(v_starts);
		v_starts = next;
	 }

	while(dim > 1)
	{
	    pv_starts = pv_starts->next;
	    --dim;
	}
	 pv_starts->next = next;
	 return next;

}

static Boolean is_gap_segment(Int4Ptr c_starts, Int2 dim)
{
	Int2 i;
	Int2 numgap = 0;

	for(i =0; i<dim; ++i)
	{
		if(c_starts[i] == -1)
			++numgap;
	}
	return (numgap >= (dim -1));
}

static void mod_prev_start(ValNodePtr p_starts, Int4Ptr c_starts, Int4Ptr c_strands, Int2 dim)
{
	Int2 i;

	for(i =0; i<dim; ++i)
	{
		if(c_strands[i] == (Int4)Seq_strand_minus)
			p_starts->data.intvalue = c_starts[i];
		p_starts = p_starts->next;
	}
}

/****************************************************************************
***
*	return TRUE if needs another round of modification
*
*****************************************************************************
***/
static Boolean mod_alignment(ValNodePtr v_starts, ValNodePtr v_strands, ValNodePtr v_lens, Int2 dim, Int2 m_order, Int2Ptr n_segs)
{
    ValNodePtr pv_starts, pv_lens, pv_strands;
    Boolean has_prev;
    Int4 c_len, p_len=0;
    Int4 p_start, c_start;
    Int2 i;
    Boolean retval;
    Int4Ptr cur_starts, prev_starts;
    Int4Ptr cur_strands, prev_strands;
    Boolean merge;

	retval = FALSE;
	cur_starts = (Int4Ptr) MemNew((size_t)dim * sizeof(Int4));
	prev_starts = (Int4Ptr) MemNew((size_t)dim * sizeof(Int4));
	cur_strands = (Int4Ptr) MemNew((size_t)dim * sizeof(Int4));
	prev_strands = (Int4Ptr) MemNew((size_t)dim * sizeof(Int4));
	
	has_prev = FALSE;
	pv_starts = NULL;
	pv_strands = NULL;
	pv_lens = NULL;
	while(v_starts && v_strands && v_lens)
	{
		get_seg(v_starts, dim, cur_starts);
		get_seg(v_strands, dim, cur_strands);
		c_len = v_lens->data.intvalue;
		merge = FALSE;
		
		if(has_prev) /*check if the two can be merged*/
		{
			merge = TRUE;
			for(i =0; i<dim; ++i)
			{
				p_start = prev_starts[i];
				c_start = cur_starts[i];
				if(cur_strands[i] != prev_strands[i])
				{
					merge = FALSE;
					break;
				}
				if(p_start == -1 && c_start != -1)
				{
					merge = FALSE;
					break;
				}
				if(p_start != -1 && c_start == -1)
				{
					merge = FALSE;
					break;
				}
				if(p_start != -1 && c_start != -1)
				{
					if(cur_strands[i] == (Int4)Seq_strand_minus)
					{
						if(c_start + c_len != p_start)
						{
							merge = FALSE;
							break;
						}
					}
					else
					{
                                                if(p_start + p_len != c_start)
						{
							merge = FALSE;
							break;
						}
					}
				}
			}
		}
		if(merge)
		{
			p_len += c_len;
			pv_lens->data.intvalue = p_len;
			mod_prev_start(pv_starts, cur_starts, cur_strands, dim);
			v_starts = Free_node(v_starts, pv_starts, dim);
			v_strands = Free_node(v_strands, pv_strands, dim);
			v_lens = Free_node(v_lens, pv_lens, 1);
			-- (*n_segs);
			retval = TRUE;
		}
		else
		{
			pv_starts = v_starts;
			pv_strands = v_strands;
			pv_lens = v_lens;
			MemCopy(prev_starts, cur_starts, (size_t)dim * sizeof(Int4));
			MemCopy(prev_strands, cur_strands, (size_t)dim * sizeof(Int4));
			p_len = c_len;
			for(i =0; i<dim; ++i)
			{
				v_starts = v_starts->next;
	  			v_strands = v_strands->next;
			}
			v_lens = v_lens->next;
			has_prev = TRUE;
 
		}
		
	}
	
	MemFree(prev_starts);
	MemFree(cur_starts);
	MemFree(prev_strands);
	MemFree(cur_strands);
	return retval;

}

static void load_new_data(Int4Ptr starts, Uint1Ptr strands, Int4Ptr lens, ValNodePtr v_starts, ValNodePtr v_strands, ValNodePtr v_lens, Int2 dim, Int2 cur_seg)
{
  Int2 i, j;

	for (i =0; i<cur_seg; ++i)
	{
	   for(j =0; j<dim; ++j)
	   {
	      starts[i*dim+j] = v_starts->data.intvalue;
	      strands[i*dim+j] = (Uint1)(v_strands->data.intvalue);
	      v_starts = v_starts->next;
	      v_strands = v_strands->next;
	    }
	    lens[i] = v_lens->data.intvalue;
	    v_lens = v_lens->next;
	}

}


static SeqIdPtr get_nth_id(SeqIdPtr ids, Int2 order)
{
   Int2 i;

	i =0;
	while(ids){
		if(order == i)
			return ids;
		++i;
		ids = ids->next;
	}

	return NULL;
}

static void count_seg_num(ValNodePtr vnp, Int2 dim)
{
	Int2 num = 0;

	while(vnp)
	{
		++num;
		vnp = vnp->next;
	}

	printf("num is %d, dim is %d\n", num, dim);
}


static DenseSegPtr make_dsp(ValNodePtr v_starts, ValNodePtr v_lens, ValNodePtr v_strands, Int2 m_order, Int2 dim, DenseSegPtr prev, Int2 cur_seg, Boolean merge, SeqIdPtr m_sip, SeqIdPtr s_sip)
{
   Int4Ptr starts, lens;
   Uint1Ptr strands;
   DenseSegPtr dsp;

    /*while(merge)*/
   if(merge)
         merge = mod_alignment(v_starts, v_strands, v_lens, dim, m_order, &cur_seg);
    starts = (Int4Ptr) MemNew((size_t)(dim* cur_seg) * sizeof(Int4));
    lens= (Int4Ptr) MemNew((size_t)(cur_seg) * sizeof(Int4));
    strands= (Uint1Ptr) MemNew((size_t)(dim * cur_seg) * sizeof(Uint1));

	load_new_data(starts, strands, lens, v_starts, v_strands, v_lens, dim, cur_seg);
	ValNodeFree(v_starts);
 	ValNodeFree(v_strands);
	ValNodeFree(v_lens);
	if(prev)
    {
		MemFree(prev->starts);
		MemFree(prev->lens);
		MemFree(prev->strands);
		prev->numseg = cur_seg;
		prev->starts = starts;
		prev->strands = strands;
		prev->lens = lens;
		return prev;
	}
	else
	{
		if(m_sip == NULL || s_sip == NULL)
		{
			Message(MSG_ERROR, "Fail to get the new DenseSeg");
			MemFree(starts);
			MemFree(lens);
			MemFree(strands);
			return NULL;
		}
		dsp = DenseSegNew();
		dsp->dim = dim;
		dsp->ids = SeqIdDup(m_sip);
		dsp->ids->next = SeqIdDup(s_sip);
		dsp->numseg = cur_seg;
		dsp->starts = starts;
		dsp->lens = lens;
		dsp->strands = strands;
		return dsp;
	}

}


static Boolean DenseSegInsert(SeqIdPtr from_id, Int4 from, Int4 to, Uint1 strand, SeqIdPtr to_id, Int4 pos, DenseSegPtr dsp, Boolean merge, BoolPtr seq_insert)
{
   Int2 cur_seg, i, j;
   Int2 t_order;      /*the order of target sequence*/
   Int4 t_start, t_stop, t_max;
   Int4 cur_start;
   Int2 dim;
   Int4Ptr starts, lens;
   Uint1 cur_strand;
   Uint1Ptr strands;
   ValNodePtr v_starts, v_lens, v_strands;
   SeqIdPtr c_id;
   Int4 ins_len;
   Boolean add_insert_seg;
   Boolean has_source_id;
   Int4 gap_size, s_gap_size;
   Int4 s_last, t_last, s_first, t_first;
   Int2 s_order = -1;
 


	if(pos <0)
		return FALSE;
	t_max = 0;
	t_order = get_t_order(dsp->ids, to_id);
	if(t_order == -1)
		return FALSE;
	dim = dsp->dim;
	starts = dsp->starts;
	lens = dsp->lens;
	strands = dsp->strands;

	cur_seg = 0;
	v_starts = NULL;
	v_lens = NULL;
	v_strands = NULL;
	ins_len = to - from +1;
	has_source_id = FALSE;
	for(c_id = dsp->ids, i =0; c_id != NULL; c_id = c_id->next, ++i)
	{
		if(SeqIdMatch(c_id, from_id))
		{
			s_order = i;
			has_source_id = TRUE;
			break;
		}
	}
			
	
	/*check if it can be inserted before the start of the alignment*/
	gap_size = -1;
	s_gap_size = -1;
	if(*seq_insert == FALSE && has_source_id)
	{
		t_first = dsp->starts[t_order];
		cur_strand = dsp->strands[t_order];
		if(t_first != -1)
		{
			if(cur_strand == Seq_strand_minus)
			{
				t_first += dsp->lens[0];
				if(t_first <= pos)
					gap_size = pos - t_first;
			}
			else	/*plus strand*/
			{
				if(t_first >=pos )
					gap_size = t_first - pos;
			}
		}

		s_first = dsp->starts[s_order];
		cur_strand = dsp->strands[s_order];
		if(s_first != -1)
		{

			if(cur_strand == Seq_strand_minus)
			{
				s_first += dsp->lens[0];
				if(s_first <= from)
					s_gap_size = from - s_first;
			}
			else
			{
				if(s_first > 0)
				{
					--s_first;
					if(s_first >=to)
						s_gap_size = s_first - to;
				}
			}
		}
		if(gap_size >= 0 && s_gap_size >= 0  && gap_size == s_gap_size)
		{
			for(j = 0; j<dim; ++j)
			{
				cur_start = -1;
				cur_strand = dsp->strands[dim*(dsp->numseg-1) + j];
				if(j == t_order)
				{
					cur_start = pos;
				}
				if(j == s_order)
				{
					if(cur_strand == Seq_strand_minus)
						cur_start = from - gap_size;
					else
						cur_start = from;
				}
				ValNodeAddInt(&v_starts, 0, cur_start);
				ValNodeAddInt(&v_strands, 0, cur_strand);
			}
			ValNodeAddInt(&v_lens, 0, (ins_len + gap_size));
			*seq_insert = TRUE;
			++cur_seg;
		}
	}



	for(i =0; i<dsp->numseg; ++i)
	{
		t_start = starts[dim*i +t_order];
		if(t_start != -1)
		{
			t_stop = t_start + lens[i] -1;
			t_max = MAX(t_stop, t_max);

			/*insertion is downstrems of the current segment*/
			/*copy everything*/
			if(pos  > t_stop)
			{	/*downstream insertions*/
				for(j=0; j<dim; ++j)
				{
					cur_start = starts[dim*i + j];
					cur_strand = (strands[dim*i + j]);
					ValNodeAddInt(&v_starts, 0, cur_start);
					ValNodeAddInt(&v_strands, 0, (Int4)cur_strand);
				}
				ValNodeAddInt(&v_lens, 0, lens[i]);
				++cur_seg;
			}
			
			/*insertion is upstream, add offset*/
			if(pos < t_start)
			{
				for(j=0; j<dim; ++j)
				{
					cur_start = starts[dim*i + j];
					cur_strand = (strands[dim*i + j]);
					if(j == t_order)
						cur_start += ins_len;
					ValNodeAddInt(&v_starts, 0, cur_start);
					ValNodeAddInt(&v_strands, 0, (Int4)cur_strand);
				}
				ValNodeAddInt(&v_lens, 0, lens[i]);
				++cur_seg;
			}

			/*insert within the current segment*/
			if(pos >= t_start && pos <= t_stop) /*insert within*/
			{				
				/*if the segment is broken, add the first half*/
				if(pos > t_start)
				{
					for(j = 0; j<dim; ++j)
					{
						cur_start = starts[dim*i+j];
						cur_strand = strands[dim*i + j];
						ValNodeAddInt(&v_starts, 0, cur_start);
						ValNodeAddInt(&v_strands, 0, (Int4)cur_strand);
					}
					ValNodeAddInt(&v_lens, 0, (pos - t_start));
					++cur_seg;
				}

				/*add the inserted segment*/
				if(i == 0)
					add_insert_seg = (has_source_id && *seq_insert == FALSE);
				else
					add_insert_seg = TRUE;
				
				if(add_insert_seg)
				{
					for(j = 0; j<dim ; ++j)
					{
						cur_start = -1;
						cur_strand = strands[dim*i + j];
						if(has_source_id && j == s_order)
						{
							*seq_insert = TRUE;
							cur_start = from;
							cur_strand = strand;
						}
						if(j == t_order)
							cur_start = pos;
						ValNodeAddInt(&v_starts, 0, cur_start);
						ValNodeAddInt(&v_strands, 0, (Int4)cur_strand);
					}
					ValNodeAddInt(&v_lens, 0, ins_len);
					++cur_seg;
				}
					
				/*add the leftover segment*/
				for(j=0; j<dim; ++j)
				{
					cur_start = starts[dim*i+j];
					cur_strand = strands[dim*i + j];
					if(j == t_order)
						cur_start = pos + ins_len;
					ValNodeAddInt(&v_starts, 0, cur_start);
					ValNodeAddInt(&v_strands, 0, (Int4)cur_strand);
				}
				ValNodeAddInt(&v_lens, 0, (lens[i] - (pos - t_start)));
				++cur_seg;
			}
		}
		
		else	/*gaps in the target sequence*/
		{
			for(j=0; j<dim; ++j)
			{
				cur_start = starts[dim*i + j];
				cur_strand = (strands[dim*i + j]);
				ValNodeAddInt(&v_starts, 0, cur_start);
				ValNodeAddInt(&v_strands, 0, (Int4)cur_strand);
			}
			ValNodeAddInt(&v_lens, 0, lens[i]);
			++cur_seg;
		}
	}


	/*attach to the end of the sequence*/
	gap_size = -1;
	s_gap_size = -1;
	if(*seq_insert == FALSE && has_source_id)
	{
		t_last = dsp->starts[dim*(dsp->numseg-1) +t_order];
		cur_strand = dsp->strands[dim*(dsp->numseg-1) + t_order];
		if(t_last != -1)
		{
			if(cur_strand == Seq_strand_minus)
			{
				t_last -= 1;
				if(pos <= t_last)
					gap_size = t_last - pos;
			}
			else	/*plus strand*/
			{
				t_last += dsp->lens[dsp->numseg-1];
				if(pos >= t_last)
					gap_size = pos - t_last;
			}
		}

		s_last = dsp->starts[dim*(dsp->numseg-1) +s_order];
		cur_strand = dsp->strands[dim*(dsp->numseg-1) + s_order];
		if(s_last != -1)
		{

			if(cur_strand == Seq_strand_minus)
			{
				s_last -= 1;
				s_gap_size = s_last - to;
			}
			else
			{
				s_last += dsp->lens[dsp->numseg-1];
				s_gap_size = from - s_last;
			}
		}
		if(gap_size >=0 && s_gap_size >=0 && gap_size == s_gap_size)
		{
			for(j = 0; j<dim; ++j)
			{
				cur_start = -1;
				cur_strand = dsp->strands[dim*(dsp->numseg-1) + j];
				if(j == t_order)
				{
					if(cur_strand == Seq_strand_minus)
						cur_start = t_last - (ins_len + gap_size -1);
					else
						cur_start = t_last;
				}
				if(j == s_order)
				{
					if(cur_strand == Seq_strand_minus)
						cur_start = from;
					else
						cur_start = from - gap_size;
				}
				ValNodeAddInt(&v_starts, 0, cur_start);
				ValNodeAddInt(&v_strands, 0, cur_strand);
			}
			ValNodeAddInt(&v_lens, 0, (ins_len + gap_size));
			*seq_insert = TRUE;
			++cur_seg;
		}
	}
			
		
	dsp = make_dsp(v_starts, v_lens, v_strands, t_order, dim, dsp, cur_seg, merge, NULL, NULL);
	return TRUE;
}


static ValNodePtr store_data(Int4 val, ValNodePtr head)
{
	ValNodeAddInt(&head, 0, val);
	return head;
}

static ValNodePtr store_reverse_data(Int4 val, ValNodePtr head)
{
   ValNodePtr new_vnp;

	new_vnp = ValNodeNew(NULL);
	new_vnp->data.intvalue = val;

	new_vnp->next = head;
	return new_vnp;
}

static ValNodePtr find_last_node(ValNodePtr head)
{
	if(head != NULL)
	{
	   while(head->next != NULL)
	      head = head->next;
	   return head;
	}
	return NULL;

}

static ValNodePtr link_to_end(ValNodePtr new_vnp, ValNodePtr head)
{
     ValNodePtr l_node;

     if(head == NULL)
	return new_vnp;
      l_node = find_last_node(head);
      l_node->next = new_vnp;
	return head;
}

static ValNodePtr find_nth_node(ValNodePtr vnp, Int2 order)
{
   Int2 i;

	i =0;
	while(vnp ){
	   if(i == order)
	      return vnp;
	   vnp = vnp->next;
	   ++i;
	}
	return NULL;
}

static Int4 get_last_pos(ValNodePtr v_starts, ValNodePtr v_lens)
{
	ValNodePtr prev;
	Int4 start, pos;

	prev = NULL;
	while(v_starts->next != NULL)
	{
		prev = v_starts;
		v_starts = v_starts->next;
	}
	if(prev == NULL)
		return -1;
	start = prev->data.intvalue;
	while(v_lens->next != NULL)
		v_lens = v_lens->next;
	pos = start + v_lens->data.intvalue -1;
	return pos;
}




/*************************************************************************
***
*	only consider the master sequence is on plus strand
*	and m_order is always == 0
*
**************************************************************************
***/
static Boolean merge_two_seg(ValNodePtr PNTR v_starts, ValNodePtr vs_starts, ValNodePtr PNTR v_lens, ValNodePtr vs_lens, ValNodePtr PNTR v_strands, ValNodePtr vs_strands, Int4 pos, Uint1 s_strand, Int2 numseg)
{
      Int4 gap, last_pos;
      ValNodePtr l_start, l_len;

      gap = (*v_starts)->data.intvalue - (pos);
      if(gap >=0)
      {
	  /*if(gap > 0)
	  {
	     *v_starts = store_reverse_data(-1, *v_starts);
	     *v_starts = store_reverse_data(pos, *v_starts);
	     *v_strands = store_reverse_data(s_strand, *v_strands);
	     *v_strands = store_reverse_data(1, *v_strands);
	     *v_lens = store_reverse_data(gap+1, *v_lens);
	  }*/
	  last_pos = get_last_pos(vs_starts, vs_lens);
	  if(last_pos == -1)
		return FALSE;

	  gap = (*v_starts)->data.intvalue - last_pos -1;
	  if(gap >0)
		return FALSE;
	  *v_starts = link_to_end(*v_starts, vs_starts);
	  *v_strands = link_to_end(*v_strands, vs_strands);
	  *v_lens = link_to_end(*v_lens, vs_lens);
	  return TRUE;
      }

      l_start = find_nth_node(*v_starts, 2*(numseg-1));
      l_len= find_nth_node(*v_lens, (numseg-1));
      last_pos = l_start->data.intvalue + l_len->data.intvalue -1;
      gap = pos - last_pos;
      if(gap >=0)
      {
	  /*if(gap >0)
	  {
	     *v_starts = store_data(pos, *v_starts);
	     *v_starts = store_data(-1, *v_starts);
	     *v_strands = store_data(1, *v_strands);
	     *v_strands = store_data(s_strand, *v_strands);
	     *v_lens = store_data(gap+1, *v_lens);
	   }*/
	  gap = vs_starts->data.intvalue - last_pos -1;
	  if(gap >0)
		return FALSE;
	  *v_starts = link_to_end(vs_starts, *v_starts);
	  *v_strands = link_to_end(vs_strands, *v_strands);
	  *v_lens = link_to_end(vs_lens, *v_lens);
	  return TRUE;
	}

	return FALSE;
}


/************************************************************************
***
*	load the portion of alignment (from, to, strand) in f_dsp to
*	the interval where to_id points to
*	from_id has already been inserted into to_id, so the region
*	can be mapped to the interval on to_id
*
*************************************************************************
***/
static Boolean dsp_process(DenseSegPtr f_dsp, Int4 from, Int4 to, Uint1 strand, SeqIdPtr from_id, SeqIdPtr to_id, Int4 pos, ValNodePtr PNTR vs_starts, ValNodePtr PNTR vs_strands, ValNodePtr PNTR vs_lens)
{
   Int4Ptr starts, lens;
   Uint1Ptr strands;
   Uint1 i_strand;
   SeqLocPtr m_loc;
   Int4 m_start, m_stop, m_pos;
   Int4 off_start, off_stop;
   Int4 t_start, t_stop;
   Int4 s_start, s_len;
   Int2 m_order, i;
   Boolean is_found;
   BioseqPtr tbsp;


   if(f_dsp == NULL)
	return FALSE;
   tbsp = BioseqFind(to_id);
   if(tbsp == NULL)
	return FALSE;
   starts = f_dsp->starts;
   lens = f_dsp->lens;
   strands = f_dsp->strands;
   m_loc = SeqLocIntNew(0, 0, 1, from_id);
   m_pos = -1;
   m_order = get_t_order(f_dsp->ids, from_id);
   is_found = FALSE;
   *vs_starts = NULL;
   *vs_strands = NULL;
   *vs_lens = NULL;

   for(i = 0; i<f_dsp->numseg; ++i)
   {
	m_start = starts[2*i+m_order];
	if(m_start != -1){
	    m_stop = m_start + lens[i] -1;
	    m_pos = m_stop +1;
	}
	else
	{
	    m_start = m_pos;
	    m_stop = m_pos;
	}
	if(m_start != -1 && m_stop != -1)
	{
	    if(!(m_start > to) || (m_stop < from))
	    {
	       off_start = MAX(0, from-m_start);
	       off_stop = MAX(0, m_stop - to);
	       /*if(tbsp->repr == Seq_repr_seg)
	       {
	         m_start = MAX(m_start, from);
	         m_stop = MIN(m_stop, to);
	         update_seq_loc(m_start, m_stop, strand, m_loc);
	         t_start = GetOffsetInBioseq(m_loc, tbsp, SEQLOC_LEFT_END);
	         t_stop = GetOffsetInBioseq(m_loc, tbsp, SEQLOC_RIGHT_END);
	       }
	       else
	       {*/
		 if(strand == Seq_strand_minus)
		 {
		    t_start = (pos-1) -(to - from) +  MAX((to-m_stop), 0);
		    t_stop = pos -1 - MAX((m_start-from), 0);
		 }
		 else
		 {
		    t_start = (pos-1) - (to - from)  + MAX((m_start-from), 0);
		    t_stop = pos-1 - MAX((to - m_stop), 0);
		 }
		/*}*/
	       if(t_start != -1 && t_stop != -1)
	       {
	       swap(&t_start, &t_stop);
	       s_start = starts[2*i +1-m_order];
	       if(s_start != -1)
	       {
	         if(strands[2*i] == strands[2*i+1])
		   s_start += off_start;
		 else
		   s_start += off_stop;
		}
		s_len = lens[i] - (off_start + off_stop);
		if(s_len >0)
		{
		 is_found = TRUE;
		 i_strand = strands[2*i+1-m_order];
		 if(starts[2*i+m_order] == -1)
		    t_start = -1;
		 if(strand != strands[2*i+m_order]){ /*reverse the seg order*/
		    i_strand = 3 - i_strand;
		    *vs_starts = store_reverse_data(s_start, *vs_starts);
		    *vs_starts = store_reverse_data(t_start, *vs_starts);
		    *vs_strands = store_reverse_data((Int4)i_strand, *vs_strands);
		    *vs_strands = store_reverse_data(1, *vs_strands);
		    *vs_lens = store_reverse_data(s_len, *vs_lens);
		  }
		  else
		  {
		    *vs_starts = store_data(t_start, *vs_starts);
		    *vs_starts = store_data(s_start, *vs_starts);
		    *vs_strands = store_data(1, *vs_strands);
		    *vs_strands = store_data((Int4)i_strand, *vs_strands);
		    *vs_lens = store_data(s_len, *vs_lens);
		  }
		 }  /*end of if(s_len>0)*/
	       }
	     }
	}
    }
    SeqLocFree(m_loc);
    return is_found;

}


static Int2 get_cur_seg(ValNodePtr lens)
{
  Int2 i=0;

	while(lens)
	{
	  ++i;
	  lens = lens->next;
	}
	return i;
}

/************************************************************************
***
*	DenseSegMerge(only for the two-dimentional pairwise alignment
*	pos is the position on the merged sequence
*
*************************************************************************
***/
static Boolean DenseSegMerge(DenseSegPtr f_dsp, SeqIdPtr from_id, Int4 from, Int4 to, Uint1 strand, SeqIdPtr to_id, Int4 pos, DenseSegPtr t_dsp)
{
   Int2 f_order, t_order;
   SeqIdPtr fs_id, ts_id;
   Int2 i, j;
   ValNodePtr v_starts, v_lens, v_strands, vs_starts, vs_lens, vs_strands;
   Int4Ptr starts, lens;
   Uint1Ptr strands;
   Int2 cur_seg;
   Uint1 i_strand;
   BioseqPtr t_bsp;

   t_order = get_t_order(t_dsp->ids, to_id);
   if(t_order == -1)
	return FALSE;
   f_order = get_t_order(f_dsp->ids, from_id);
   if(f_order == -1)
	return FALSE;

   fs_id = get_nth_id(f_dsp->ids, 1-f_order);
   ts_id = get_nth_id(t_dsp->ids, 1-t_order);
   if(!SeqIdForSameBioseq(fs_id, ts_id))
	return FALSE;

   vs_starts =NULL;
   vs_strands = NULL;
   vs_lens = NULL;
   if(!dsp_process(f_dsp, from, to, strand, from_id, to_id, pos, &vs_starts, &vs_strands, &vs_lens))		/*the  alignment is not within [from, to]*/
   	return TRUE;
   i_strand = (Uint1)(vs_strands->next->data.intvalue);
   t_bsp = BioseqFind(to_id);
   if(pos == BioseqGetLen(t_bsp))		/*insertion at the end*/
	pos = pos - (to - from +1) -1;

   v_starts =NULL;
   v_strands = NULL;
   v_lens = NULL;
   starts = t_dsp->starts;
   lens = t_dsp->lens;
   strands = t_dsp->strands;
   for(i =0; i<t_dsp->numseg; ++i)
   {
       for(j =0; j<2; ++j)
       {
          v_starts = store_data(starts[2*i +j], v_starts);
	  v_strands = store_data(strands[2*i+j], v_strands);
	}
	v_lens = store_data(lens[i], v_lens);
    }

    if(!merge_two_seg(&v_starts, vs_starts, &v_lens, vs_lens, &v_strands, vs_strands, pos, i_strand, t_dsp->numseg))	/*fail to brought togather*/
    {
	ValNodeFree(v_starts);
	ValNodeFree(vs_starts);
	ValNodeFree(v_lens);
	ValNodeFree(vs_lens);
	ValNodeFree(v_strands);
	ValNodeFree(vs_strands);
	return FALSE;
    }

    cur_seg = get_cur_seg(v_lens);
    t_dsp = make_dsp(v_starts, v_lens, v_strands, t_order, 2, t_dsp, cur_seg, TRUE, to_id, fs_id);
    return TRUE;
}

static SeqAlignPtr SeqAlignLink(SeqAlignPtr new_salp, SeqAlignPtr align)
{
   if(align == NULL)
	return new_salp;
   while(align->next != NULL)
	align = align->next;
   align->next = new_salp;
   return new_salp;

}

static SeqAlignPtr make_new_align(DenseSegPtr f_dsp, Int4 from, Int4 to, Uint1 strand, SeqIdPtr from_id, SeqIdPtr to_id, Int4 pos, SeqAlignPtr head)
{
   ValNodePtr vs_starts, vs_strands, vs_lens;
   DenseSegPtr dsp;
   Int2 cur_seg;
   SeqAlignPtr align;
   SeqIdPtr fs_id;

   if(!dsp_process(f_dsp, from, to, strand, from_id, to_id, pos, &vs_starts, &vs_strands, &vs_lens))		/*the  alignment is not within [from, to]*/
	return head;

    fs_id = f_dsp->ids->next;
    cur_seg = get_cur_seg(vs_lens);
    dsp = make_dsp(vs_starts, vs_lens, vs_strands, 0, 2, NULL, cur_seg, TRUE, to_id, fs_id);
    if(dsp)
    {
       align = SeqAlignNew();
       align->segtype = 2;
       align->segs = dsp;
       align->type =3;
       if(head == NULL)
	  head = align;
	else
	  SeqAlignLink(align, head);
    }
    return head;
}



static SeqAlignPtr JZSeqAlignMerge(SeqAlignPtr f_align, SeqIdPtr from_id, Int4 from, Int4 to, Uint1 strand, SeqIdPtr to_id, Int4 pos, SeqAlignPtr t_align)
{
   SeqAlignPtr new_salp, curr;
   DenseSegPtr f_dsp, t_dsp;
   Boolean is_found;

	if(f_align == NULL)
		return t_align;

	new_salp = NULL;
	while(f_align)
	{
	   if(f_align->segtype == 2)
	   {
	      f_dsp = (DenseSegPtr) f_align->segs;
	      is_found = FALSE;
	      curr = t_align;
	      while(curr && !is_found)
	      {
		 if(curr->segtype ==2)
		 {
		    t_dsp = (DenseSegPtr) curr->segs;
		    is_found = DenseSegMerge(f_dsp, from_id, from, to, strand, to_id, pos, t_dsp);
		 }
		 curr = curr->next;
	       }

	       if(!is_found)
	          new_salp = make_new_align(f_dsp, from, to, strand, from_id, to_id, pos, new_salp);


		f_align = f_align->next;
	     }
	}

	if(new_salp != NULL){
	       if(t_align == NULL)
		   t_align = new_salp;
	       else
	           SeqAlignLink(new_salp, t_align);
	}
	return t_align;


}



/***************************************************************************
***
*	AttachSeqInAlign(): for a sequence not found in the previous alignment
*	data and it is used in producing the alignment, it would be attached
*	as a part of the alignment node for the seq-hist of the new sequence
*
******************************************************************************
***/
static SeqAlignPtr AttachSeqInAlign(SeqIdPtr from_id, Int4 from, Int4 to, Uint1 strand, SeqIdPtr to_id, Int4 pos, SeqAlignPtr head)
{
   DenseSegPtr dsp;
   ValNodePtr v_starts, v_lens, v_strands;
   SeqAlignPtr align;

   v_starts = NULL;
   v_lens = NULL;
   v_strands = NULL;


   ValNodeAddInt(&v_starts, 0, pos);
   ValNodeAddInt(&v_starts, 0, from);
   ValNodeAddInt(&v_strands, 0, (Int4)Seq_strand_plus);
   ValNodeAddInt(&v_strands, 0, (Int4)strand);
   ValNodeAddInt(&v_lens, 0, (to - from +1));
   dsp = make_dsp(v_starts, v_lens, v_strands, 0, 2, NULL, 1, TRUE, to_id, from_id);
   align = SeqAlignNew();
   align->segtype = 2;
   align->segs = dsp;
   align->type =3;
   if(head == NULL)
        head = align;
   else
	SeqAlignLink(align, head);
   return head;

}


static void mod_old_align(SeqAlignPtr PNTR old, SeqIdPtr sip)
{
	SeqAlignPtr curr, prev;
	SeqIdPtr c_sip;

	prev = NULL;
	curr = *old;

	while(curr)
	{
		c_sip = get_align_id(curr, 1);
		if(SeqIdForSameBioseq(c_sip, sip))
		{
			if(prev == NULL)
				*old= curr->next;
			else
				prev->next = curr->next;
			curr->next = NULL;
			SeqAlignFree(curr);
			return;
		}
		prev = curr;
		curr = curr->next;
	}
}

static void attach_new_align(SeqAlignPtr new_salp, SeqAlignPtr PNTR old) {
	SeqIdPtr sip;
	SeqAlignPtr align;

	if(new_salp == NULL)
		return;
	align = new_salp;
	while(align)
	{
		sip = get_align_id(align, 1);
		if(sip)
			mod_old_align(old, sip);
		align= align->next;
	}

	if(*old== NULL)
		*old = new_salp;
	else
		SeqAlignLink(new_salp, (*old));

}

static void delete_bad_node(ValNodePtr PNTR head, Uint1 choice)
{
	ValNodePtr prev, next, curr;

	curr = *head;
	prev = NULL;

	while(curr)
	{
		next = curr->next;
		if(curr->choice != choice)
		{
			if(prev == NULL)
				*head = next;
			else
				prev->next = next;
			curr->next = NULL;
			ValNodeFree(curr);
		}
		else
			prev = curr;

		curr = next;
	}
}
static Int4 count_coverage (BoolPtr used_list, Int4 size)
{
	Int4 i;
	Int4 count;

	count = 0;
	for(i = 0; i<size; ++i)
	{
		if(used_list[i] != 0)
			++count;
	}
	return count;
}


static Int2 filter_seqid_order;
static int LIBCALLBACK CompareChainProc(VoidPtr ptr1, VoidPtr ptr2)
{
  ValNodePtr   vnp1;
  ValNodePtr   vnp2;
  DenseDiagPtr ddp1, ddp2;
  Int4 start1, stop1, start2, stop2;

  start1 = -1;
  start2 = -1;
  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      ddp1 = (DenseDiagPtr) vnp1->data.ptrvalue;
      ddp2 = (DenseDiagPtr) vnp2->data.ptrvalue;
	  start1 = ddp1->starts[filter_seqid_order];
	  stop1 = start1 + ddp1->len -1;
	  start2 = ddp2->starts[filter_seqid_order];
	  stop2 = start2 + ddp2->len -1;
      if(stop1> stop2)
		return 1;
      else if(stop1 < stop2)
		return -1;
      else if (stop1 == stop2)
      {
		if(ddp1->len > ddp2->len)
			return -1;
		else if(ddp1->len < ddp2->len)
			return 1;
		else
			return 0;
      }
    }
  }
  return 0;
}

		
/*
*	return 0 for no overlap
*	return 1 for _2 overlap too much with _1
*	return 2 for _1 overlap too much with _1
*/
static Uint1 has_too_much_overlap (Int4 start_1, Int4 stop_1, Int4 start_2, Int4 stop_2)
{
	Int4 m_start, m_stop;
	FloatHi value;
	Int4 len_1, len_2;

	if(start_2 > stop_1 || stop_2 < start_1) /*no overlap */
		return 0;

	if(start_2 >= start_1 && stop_2 <= stop_1)	/*completely contained within*/
		return 1;

	if(start_1 >= start_2 && stop_1 <= start_2)
		return 2;

	m_start = MAX(start_1, start_2);
	m_stop = MIN(stop_1, stop_2);

	len_1 = stop_1 - start_1 +1 ;
	len_2 = stop_2 - start_2 + 1;

	value = (FloatHi)(m_stop - m_start + 1)/(FloatHi)(MIN(len_1, len_2));
	if(value >= 0.7)
	{
		if(len_1 >= len_2)
			return 1;
		else
			return 2;
	}
	else
		return 0;
}
	


	
/*
*	return TRUE if the ddp needs to be deleted. return FALSE if ddp doesn't need to 
*	be deleted
*/
static Boolean ck_delete_next_chain (DenseDiagPtr ddp, ValNodePtr PNTR ddp_list, Int2 s_order)
{
	ValNodePtr curr, prev, next;
	DenseDiagPtr n_ddp;
	Uint1 overlap_val;
	Int4 start_1, stop_1, start_2, stop_2;

	curr = *ddp_list;
	prev = NULL;
	
	start_1 = ddp->starts[s_order];
	stop_1 = ddp->starts[s_order] + ddp->len -1;


	while(curr)
	{
		next = curr->next;
		n_ddp = (DenseDiagPtr) curr->data.ptrvalue;
		start_2 = n_ddp->starts[s_order];
		stop_2 = n_ddp->starts[s_order] + n_ddp->len -1;
		overlap_val = has_too_much_overlap (start_1, stop_1, start_2, stop_2);
		if(overlap_val == 2)
		{
			return TRUE;
		}
		else if(overlap_val == 1)
		{
			if(prev == NULL)
				*ddp_list = next;
			else
				prev->next = next;
			curr->next = NULL;
			ValNodeFree(curr);
		}
		else if(overlap_val == 0)
			prev = curr;
		curr = next;
	}

	return FALSE;
}


static void delete_redundant_node (ValNodePtr PNTR ddp_list, Int2 s_order)
{
	ValNodePtr curr, prev, next;
	DenseDiagPtr ddp;

	prev = NULL;
	curr = *ddp_list;
	while(curr)
	{
		next = curr->next;
		ddp = (DenseDiagPtr) curr->data.ptrvalue;
		if(ck_delete_next_chain (ddp, &next, s_order))
		{
			if(prev == NULL)
				*ddp_list = next;
			else
				prev->next = next;
		}
		else
		{
			prev = curr;
			curr->next = next;
		}

		curr = next;
	}
}

static Int4Ptr convert_ddp_to_array (ValNodePtr ddp_list, Uint1 strand, Int2 s_order, Int4Ptr p_num)
{
	Int4Ptr val; 
	Int4 num;
	ValNodePtr curr;
	DenseDiagPtr ddp;
	Int4 i;

	num = 0;
	for(curr = ddp_list; curr != NULL; curr = curr->next)
	{
		++num;
	}

	val = (Int4Ptr) MemNew((size_t)num * sizeof(Int4));

	/*if strand == minus, check the descending order of the stop. 
	  otherwise, check the ascending order of the start */

	i = 0;
	for(curr = ddp_list; curr != NULL; curr = curr->next)
	{
		ddp = (DenseDiagPtr) curr->data.ptrvalue;
		if(strand == Seq_strand_minus)
			val[i++] = ddp->starts[1-s_order] + ddp->len -1;
		else
			val[i++] = ddp->starts[1-s_order];
	}

	*p_num = num;
	return val;
}

static Boolean clean_up_discrepancies (Int4Ptr val, Int4 size, Boolean descend, Int4Ptr score_order)
{
	Int4 i, j;
	Int4 num_pos, num_neg;
	Int4Ptr temp;
	Int4 order;

	if(size == 1)
		return TRUE;
	temp = (Int4Ptr) MemNew((size_t)size * sizeof(Int4));
	MemSet(temp, -1, (size_t)size * sizeof(Int4));

	for(i = 0; i<size; ++i)
	{
		order = score_order[i];
		if(i == 0)
			temp[order] = val[order];
		else
		{
			num_pos = 0;
			num_neg = 0;
			for(j = 0; j<size; ++j)
			{
				if(temp[j] != -1 && j != order)
				{
					if(descend)
					{
						if(order > j)
						{
							if(temp[j] > val[order])
								++num_pos;
							else
								++num_neg;
						}
						else
						{
							if(temp[j] > val[order])
								++num_neg;
							else
								++num_pos;
						}
					}
					else
					{
						if(order > j)
						{
							if(temp[j] > val[order])
								++num_neg;
							else
								++num_pos;
						}
						else
						{
							if(temp[j] > val[order])
								++num_pos;
							else
								++num_neg;
						}
					}
				}
			}
			if(num_neg > 0)
				temp[order] = -1;
			else
				temp[order] = val[order];
		}
	}

	MemCopy(val, temp, (size_t)size * sizeof(Int4));
	MemFree(temp);
	return TRUE;

}
					
				
/*
*
*	clean up the stuffs in the chain of ddp_list that contains out of order 
*	genomic seqeunce list
*
*/
static void clean_genomic_order (ValNodePtr PNTR head, Int2 s_order, Uint1 strand, Int4Ptr score_order)
{

	Int4Ptr val;
	Int4 size, i;
	ValNodePtr ddp_list, prev, next, curr;

	ddp_list = *head;

	/*convert the ddp_list to an arrya of the genomic sequences */
	val = convert_ddp_to_array (ddp_list, strand, s_order, &size);

	/*identify the out of bound entries, and set the value to -1 */
	clean_up_discrepancies (val, size, (Boolean)(strand == Seq_strand_minus), score_order);

	/*do the real clean up work */
	i = 0;
	curr = *head; 
	prev = NULL;
	while(curr)
            {		
                /* DenseSegPtr ddp = curr->data.ptrvalue;
		 printf("%ld, %ld, len = %ld\n", ddp->starts[s_order], ddp->starts[1-s_order], ddp->len);  */
		next = curr->next; 
		if(val[i] == -1)
		{
			/* printf("clean\n");  */
			if(prev == NULL)
				*head = next;
			else
				prev->next = next;
			curr->next = NULL;
			ValNodeFree(curr);
		}
		else
			prev = curr;
		curr = next;
		++i;
	}
	MemFree(val);
}
		

static Int4 get_score_from_list (ScorePtr scores)
{
	ObjectIdPtr oip;

	if(scores == NULL)
		return 0;
	while(scores)
	{
		oip = scores->id;
		if(oip && StringCmp(oip->str, "score") == 0)
			return scores->value.intvalue;
		scores = scores->next;
	}

	return 0;
}

	
static Int4 get_score_from_DenseDiag (DenseDiagPtr ddp)
{
	return get_score_from_list (ddp->scores);
}
	
static Boolean filter_this_seg (DenseDiagPtr curr, DenseDiagPtr PNTR h_list, Uint1 order)
{
	Int4 start, stop, t_start, t_stop;
	Int4 a_start, a_stop;
	Int4 score, t_score;
	DenseDiagPtr prev, next, list;
	Boolean found;
	Int4 len;

	start = curr->starts[order];
	stop = start + curr->len -1;
	score = get_score_from_DenseDiag (curr);
	prev = NULL;
	list = *h_list;
	while(list)
	{
		next = list->next;
		t_start = list->starts[order];
		t_stop = t_start + list->len -1;
		t_score = get_score_from_DenseDiag (list);
		
		found = FALSE;
		if(!(t_start > stop || t_stop < start))
		{
			a_start = MAX(t_start, start);
			a_stop = MIN(t_stop, stop);
			len = MAX(t_stop, stop) - MIN(t_start, start);
			if((FloatHi)(a_stop - a_start)/(FloatHi)len > 0.80)
			{
				if(score > t_score)
				{
					if(prev == NULL)
						*h_list = next;
					else
						prev->next = next;
					list->next = NULL;
					DenseDiagFree(list);
					found = TRUE;
				}
			}
		}

		if(!found)
			prev = list;
		list = next;
	}
		
	list = *h_list;
	while(list)
	{
		t_start = list->starts[order];
		t_stop = t_start+ list->len -1;
		t_score = get_score_from_DenseDiag (list);
		
		if(!(t_start > stop || t_stop < start))
		{
			a_start = MAX(t_start, start);
			a_stop = MIN(t_stop, stop);
			if((FloatHi)(a_stop - a_start)/(FloatHi)(stop - start) > 0.50)
			{
				if(score < t_score)
					return TRUE;
			}
		}
		list = list->next;
	}

	return FALSE;
}


static Boolean is_same_orientation(Uint1 strand_1, Uint1 strand_2)
{
	if(strand_1 == strand_2)
		return TRUE;
	if(strand_1 == Seq_strand_minus)
		return (strand_2 == Seq_strand_minus);
	else
		return (strand_2 != Seq_strand_minus);
}

static void filter_wrong_orientation (SeqAlignPtr align, Boolean same_orient)
{
	DenseDiagPtr ddp, prev, next;
	Boolean del;

	ddp = (DenseDiagPtr) align->segs;
	prev = NULL;

	while(ddp)
	{
		next = ddp->next;
		del = (is_same_orientation(ddp->strands[0], ddp->strands[1]) != same_orient);
		if(del)
		{
			if(prev == NULL)
				align->segs = next;
			else
				prev->next = next;
			ddp->next = NULL;
			DenseDiagFree(ddp);
		}
		else
			prev = ddp;
		ddp = next;
	}
}


static FloatHi get_overlap(BoolPtr status, Int4 from, Int4 len, Int4Ptr poverlap)
{
	Int4 i;
	Int4 overlap_len;

	overlap_len = 0;
	for(i = 0; i<len; ++i)
	{
		if(status[i+from])
			++overlap_len;
	}

	*poverlap = overlap_len;
	return (FloatHi)overlap_len/(FloatHi)len;
}


static Boolean is_overlap(BoolPtr status, Int4 from, Int4 len)
{
	Int4 i;
	Int4 overlap_len;

	overlap_len = 0;
	for(i = 0; i<len; ++i)
	{
		if(status[i+from])
			++overlap_len;
	}

	if((FloatHi)overlap_len/(FloatHi)len> 0.7)
		return TRUE;
	else
	{
		MemSet(status+from, 1, (size_t)(len) * sizeof(Boolean));
		return FALSE;
	}
}


static void set_ddp_status(ValNodePtr ddp_list, Int4 first_len, Int4 second_len)
{
	BoolPtr first_status, second_status;
	DenseDiagPtr ddp;
	FloatHi val_1, val_2;
	Int4 overlap_1, overlap_2;

	first_status = (BoolPtr) MemNew((size_t)first_len * sizeof(Boolean));
	second_status = (BoolPtr) MemNew((size_t)second_len * sizeof(Boolean));

	while(ddp_list)
	{
		ddp = (DenseDiagPtr) ddp_list->data.ptrvalue;
		val_1 = get_overlap(first_status, ddp->starts[0], ddp->len, &overlap_1);
		val_2 = get_overlap(second_status, ddp->starts[1], ddp->len, &overlap_2);

		if(val_1 >= 0.95 || val_2 >=0.95  
			|| (ddp->len - MAX(overlap_1, overlap_2) < MIN(50, ddp->len-10)) 
			|| (val_1 > 0.7 && val_2 >0.7))
		{
			/* if(ddp->starts[0] == 35)
				printf("val_1 = %lf, val_2 = %lf, overlap_1 = %ld, overlap_2 = %ld, len = %ld\n", val_1, val_2, overlap_1, overlap_2, ddp->len);
			*/
			ddp_list->choice = 0;
		}
		else
		{
			MemSet(first_status+ddp->starts[0], 1, (size_t)(ddp->len) * sizeof(Boolean));
			MemSet(second_status+ddp->starts[1], 1, (size_t)(ddp->len) * sizeof(Boolean));
			ddp_list->choice = 1;
		}
		ddp_list = ddp_list->next;
	}

	MemFree(first_status);
	MemFree(second_status);
}

static Boolean free_diag_in_list(ValNodePtr ddp_list, DenseDiagPtr ddp)
{
	DenseDiagPtr c_ddp;

	while(ddp_list)
	{
		c_ddp = (DenseDiagPtr) ddp_list->data.ptrvalue;
		if(c_ddp == ddp)
			return (ddp_list->choice == 1);
		ddp_list = ddp_list->next;
	}

	return FALSE;
}
		
static void filter_repeats (ValNodePtr PNTR ddp_list, Int4 first_len, Int4 second_len)
{

	set_ddp_status(*ddp_list, first_len, second_len);
	delete_bad_node(ddp_list, 1);
}

static Int4Ptr get_process_list (ValNodePtr ddp_list)
{
	Int4 num, i, j;
	ValNodePtr curr;
	Int4Ptr order, score_list;
	DenseDiagPtr ddp;
	Boolean change = TRUE;
	Int4 temp;

	num = 0;
	for(curr = ddp_list; curr != NULL; curr = curr->next)
		++num;

	score_list = (Int4Ptr) MemNew((size_t)num * sizeof(ValNodePtr));
	order = (Int4Ptr) MemNew((size_t)num * sizeof(Int4));

	i = 0;
	for(curr = ddp_list; curr != NULL; curr = curr->next)
	{
		ddp = (DenseDiagPtr) curr->data.ptrvalue;
		score_list[i] = get_score_from_DenseDiag (ddp);
		order[i] = i;
		++i;
	}

	for(i = 0; i<num-1 && change; ++i)
	{
		change = FALSE;
		for(j = 0; j<num-i-1; ++j)
		{
			if(score_list[j+1] > score_list[j])
			{

				change =TRUE;
				temp = score_list[j];
				score_list[j] = score_list[j+1];
				score_list[j+1] = temp;

				temp = order[j];
				order[j] = order[j+1];
				order[j+1] = temp;
			}
		}
	}

	MemFree(score_list);
	return order;
}


static ValNodePtr get_nth_node(ValNodePtr vnp, Int4 order)
{
	Int4 i;

	i = 0;
	while(vnp)
	{
		if(i == order)
			return vnp;
		++i;
		vnp = vnp->next;
	}

	return NULL;
}

static Uint1 get_alignment_orientation(ValNodePtr PNTR pddp_list, BioseqPtr bsp, Int2 s_order)
{
	Int4 num;
	Int4Ptr score_order;
	ValNodePtr curr;
	DenseDiagPtr ddp;
	BoolPtr plus_used_list, minus_used_list;
	Int4 plus_count, minus_count;
	Int4 i;
	ValNodePtr list, ddp_list;
	Uint1 strand, t_strand;

	ddp_list = *pddp_list;
	num = 0;
	for(curr = ddp_list; curr != NULL; curr = curr->next)
		++num;

	score_order = get_process_list (ddp_list);
	plus_used_list = (BoolPtr) MemNew((size_t)bsp->length * sizeof(Boolean));
	minus_used_list = (BoolPtr) MemNew((size_t)bsp->length * sizeof(Boolean));

	

	for(i = 0; i<MIN(num, 10); ++i)
	{
		curr = get_nth_node(ddp_list, score_order[i]);
		if(curr != NULL)
		{
			ddp = (DenseDiagPtr) curr->data.ptrvalue;
			if(ddp->strands[0] != ddp->strands[1] && 
				(ddp->strands[0] == Seq_strand_minus || 
					ddp->strands[1] == Seq_strand_minus))
				MemSet(minus_used_list+ddp->starts[s_order], 1,
                        	(size_t)(ddp->len) * sizeof(Boolean));
			else
				MemSet(plus_used_list+ddp->starts[s_order], 1,
                        	(size_t)(ddp->len) * sizeof(Boolean));
		}
	}

	plus_count = count_coverage (plus_used_list, bsp->length);
	minus_count = count_coverage (minus_used_list, bsp->length);
	MemFree(plus_used_list);
	MemFree(minus_used_list);

	if(plus_count == 0 && minus_count == 0)
	{
		MemFree(score_order);
		return 0;
	}
	else if(plus_count > minus_count)
		strand = Seq_strand_plus;
	else
		strand = Seq_strand_minus;

	/*only load the ddp_list with the right orientation*/
	list = NULL;
	for(i = 0; i<num; ++i)
	{
		curr = get_nth_node(ddp_list, score_order[i]);
		ddp = (DenseDiagPtr) curr->data.ptrvalue;
		if(ddp->strands[0] != ddp->strands[1] && 
			(ddp->strands[0] == Seq_strand_minus || 
				ddp->strands[1] == Seq_strand_minus))
			t_strand = Seq_strand_minus;
		else
			t_strand = Seq_strand_plus;
		if(t_strand == strand)
			ValNodeAddPointer(&list, 0, ddp);
	}
	ValNodeFree(*pddp_list);
	*pddp_list = list;
	MemFree(score_order);
	return strand;


}



/*
*	check to see if ddp falls in the order of the existing list
*
*/
static Boolean load_ddp_in_order(ValNodePtr PNTR h_list, DenseDiagPtr ddp, Boolean reverse, Boolean add)
{
	Int4 m_stop, s_stop;
	Int4 mt_stop, st_stop;
	ValNodePtr curr, prev;
	DenseDiagPtr t_ddp;
	Boolean load;
	ValNodePtr vnp;
	Int4 i;

	m_stop = ddp->starts[0] + ddp->len -1;
	s_stop = ddp->starts[1] + ddp->len -1;

	for(curr = *h_list; curr != NULL; curr = curr->next)
	{
		t_ddp = (DenseDiagPtr) curr->data.ptrvalue;
		if(ddp->len <= t_ddp->len)
		{
			for(i = 0; i<2; ++i)
			{
				if(ddp->starts[i] >= t_ddp->starts[i])
				{
					if(t_ddp->starts[i] + t_ddp->len >= ddp->starts[i] + ddp->len)
						return FALSE;
				}
			}
		}
	}
	prev = NULL;
	curr = *h_list;
	while(curr)
	{
		t_ddp = (DenseDiagPtr) curr->data.ptrvalue;
		mt_stop = t_ddp->starts[0] + t_ddp->len -1;
		st_stop = t_ddp->starts[1] + t_ddp->len -1;
		if(m_stop <= mt_stop)
		{
			if(reverse)
				load = (s_stop >= st_stop);
			else
				load = (s_stop <= st_stop);
			/*check with the previous order*/
			if(load)
			{
				if(prev != NULL)
				{
					t_ddp = (DenseDiagPtr) prev->data.ptrvalue;
					mt_stop = t_ddp->starts[0] + t_ddp->len -1;
					st_stop = t_ddp->starts[1] + t_ddp->len -1;

					if(m_stop < mt_stop)
						load = FALSE;
					else
					{
						if(reverse)
							load = (s_stop <= st_stop);
						else
							load = (s_stop >= st_stop);
					}
				}
			}
			if(load)
			{
				if(add)
				{
					vnp = ValNodeNew(NULL);	
					vnp->choice = 1;
					vnp->data.ptrvalue = ddp;
					if(prev == NULL)
						*h_list = vnp;
					else
						prev->next = vnp;
					vnp->next = curr;
				}
				return TRUE;
			}
			else
				return FALSE;
		}
		prev = curr;
		curr = curr->next;
	}

	if(prev != NULL)
	{
		if(reverse)
			load = (s_stop <= st_stop);
		else
			load = (s_stop >= st_stop);

		if(load)
		{
			vnp = ValNodeNew(NULL);	
			vnp->choice = 1;
			vnp->data.ptrvalue = ddp;
			prev->next = vnp;
			return TRUE;
		}
	}

	return FALSE;
}

	

static Int4Ptr load_ddp_scores(ValNodePtr ddp_list, Int4Ptr num)
{
	DenseDiagPtr ddp;
	Int4 i;
	Int4Ptr score_list;
	ValNodePtr curr;

	i = 0;
	for(curr = ddp_list; curr != NULL; curr = curr->next)
		++i;
	if(i == 0)
		return NULL;
	*num = i;
	score_list = (Int4Ptr) MemNew((size_t)i * sizeof(Int4));
	i = 0;
	for(curr = ddp_list; curr != NULL; curr = curr->next)
	{
		ddp = (DenseDiagPtr) curr->data.ptrvalue;
		score_list[i] = get_score_from_list (ddp->scores);
		++i;
	}
	return score_list;
}
	

static Boolean is_ddp_overlap(DenseDiagPtr ddp, BoolPtr first_status, BoolPtr second_status)
{
	FloatHi val_1, val_2;
	Int4 overlap_1, overlap_2;

	val_1 = get_overlap(first_status, ddp->starts[0], ddp->len, &overlap_1);
	val_2 = get_overlap(second_status, ddp->starts[1], ddp->len, &overlap_2);

	if(val_1 >=0.95 || val_2 >=0.95 || (val_1 > 0.7 && val_2 >0.7))
		return TRUE;

	/*maximum non-overlap region is less than 50-bp*/
	if(ddp->len - MAX(overlap_1, overlap_2) < MIN(50, ddp->len-10)) 
		return TRUE;
	return FALSE;
}


/*
* lns.c - Given a sequence of numbers, this program computes
*	  the Longest Nondecreasing Subsequence.
*
* Feb. 20, 1997
*
* Program by Jinghui Zhang & Kun-Mao Chao
* Toolkit version : H. Sicotte, 1/04/1999
*/

/* Input :
    call lns()
	S : a sequence of numbers;
	len_s : the length of the sequences;
   Output:
	LNS : the positions of the Longest Nondecreasing Subsequence;
	return value : the length of the LNS;
   Note: "-1" entries are allowed.
*/

/* Input :
        SS : a sequence of numbers; (without "-1" entries)
        len_ss : the length of the sequences;
   Output:
        LNS : the positions of the Longest Nondecreasing Subsequence;
        return value : the length of the LNS;
*/
static Int4 pure_lns(Int4Ptr SS, Int4 len_ss, Int4Ptr LNS)
{
	Int4 len_lns;
        Int4Ptr PREV=NULL,BEST=NULL,BEST_POS=NULL;
	Int4 i, low, high, mid;
	Int4 t;

	BEST[0] = SS[0];
	BEST_POS[0] = 0;
	PREV[0] = -1;
	len_lns = 1;
        if(len_ss) {
            PREV = (Int4Ptr) Malloc(sizeof(Int4)*len_ss);
            BEST = (Int4Ptr) Malloc(sizeof(Int4)*len_ss);
            BEST_POS = (Int4Ptr) Malloc(sizeof(Int4)*len_ss);
        }
	for (i=1; i<len_ss; ++i) {
	   t = SS[i];
	   low = 0;
	   high = len_lns-1;
	   while (low < high) {
		mid = (low + high) / 2;
		if (t >= BEST[mid]) low = mid + 1;
		else high = mid;
	   }
	   if (t < BEST[low]) {
		if (low == 0) PREV[i] = -1;
		else PREV[i] = BEST_POS[low-1];
		BEST[low] = t;
		BEST_POS[low] = i;
	   } else {
		BEST[len_lns] = t;
		BEST_POS[len_lns] = i;
		PREV[i] = BEST_POS[len_lns-1];
		++len_lns;
	   }
	}

	/* trace back */
	
	t = BEST_POS[len_lns-1];
	for (i=len_lns-1; i>=0; --i) {
		LNS[i] = t;
		t = PREV[t];
	}

        if(PREV)
            Free(PREV);
        if(BEST)
            Free(BEST);
        if(BEST_POS)
            Free(BEST_POS);
	return len_lns;
}
static Int4 lns(Int4Ptr S, Int4 len_s, Int4Ptr LNS)
{
	Int4Ptr SS=NULL;
	Int4 len_ss, len_lns;
	Int4 i, j;
        SS=(Int4Ptr)Malloc(len_s*sizeof(Int4));
	/* screen out those "-1" entries */
	len_ss = 0;
	for (i=0; i<len_s; ++i)
	   if (S[i] >= 0)
		SS[len_ss++] = S[i];

	len_lns = pure_lns(SS, len_ss, LNS);

	/* add the offset */
	for (i=0, j=0; j<len_lns; ++i)
	   if (S[i] == SS[LNS[j]]) LNS[j++] = i;

        if(SS)
            Free(SS);
	return len_lns;
}




static void sort_list_order(ValNodePtr PNTR ph_list, Uint1 strand)
{
	ValNodePtr h_list;
	ValNodePtr curr;
	DenseDiagPtr ddp;
	Int4 val;
	ValNodePtr best_list = NULL;
	

	h_list = *ph_list;
	if(h_list == NULL || h_list->next == NULL)
		return;
	ddp = (DenseDiagPtr) h_list->data.ptrvalue;
	ValNodeAddPointer(&best_list, 0, ddp);
	for(curr = h_list->next; curr != NULL; curr = curr->next)
	{
		ddp = (DenseDiagPtr) curr->data.ptrvalue;
		load_ddp_in_order(&best_list, ddp, strand==Seq_strand_minus, TRUE);
	}
	

	filter_seqid_order = 0;
	h_list = SortValNode(h_list, CompareChainProc);
        {
            Int4 len_S=0,len_lns=0,rlen_lns=0;
            Int4Ptr S=NULL,rS=NULL,LNS=NULL,rLNS=NULL,this_LNS;
            Int4 i,j;
            len_S = (Int4)get_vnp_num(h_list);
            if(len_S>0) {
                S=(Int4Ptr)Malloc(len_S*sizeof(Int4));
                rS=(Int4Ptr)Malloc(len_S*sizeof(Int4));
                LNS=(Int4Ptr)Malloc(len_S*sizeof(Int4));
                rLNS=(Int4Ptr)Malloc(len_S*sizeof(Int4));
                for(curr = h_list,i=0; curr != NULL; curr = curr->next,i++)
                    {
                        ddp = (DenseDiagPtr) curr->data.ptrvalue;
                        if(strand == Seq_strand_minus)
                            S[i++] = ddp->starts[1];
                        else
                            S[i++] = ddp->starts[1] + ddp->len-1;
                    }
                for(i=len_S-1,j=0;i>=0;i--,j++)
                    rS[j]=S[i];
                len_lns = lns(S, len_S, LNS);
                rlen_lns = lns(rS, len_S, rLNS);
            }
            if (len_lns <= 0 || rlen_lns <=0) {
                val=-1;
                Message(MSG_ERROR, "Fail in get_lns");
                if(len_S) {
                    Free(S);
                    Free(rS);
                    Free(LNS);
                    Free(rLNS);
                }
                exit(2);
            } else if(rlen_lns > len_lns) {
                val = 1; /* usr rLNS */
                len_lns=rlen_lns;
                this_LNS = rLNS;
            } else {
                val =0; /* use LNS */
                this_LNS = LNS;
                printf("0\n");
            }
                /* inconsistent orientation from LNS */
            if((strand == Seq_strand_minus && val == 0) || 
               (strand != Seq_strand_minus && val== 1)) 
                {
                    ValNodeFree(h_list);
                    *ph_list = best_list;
                    if(len_S) {
                        Free(S);
                        Free(rS);
                        Free(LNS);
                        Free(rLNS);
                    }
                    return;
                }            
            
            ValNodeFree(best_list);
            for(i=0;i<len_lns;i++) {
                val = this_LNS[i];
                if(strand == Seq_strand_minus)
                    val =len_S -1 - val;
                curr = get_nth_node(h_list, val);
                if(curr != NULL)
                    curr->choice = 1;
            }
            if(len_S) {
                Free(S);
                Free(rS);
                Free(LNS);
                Free(rLNS);
            }

        }
        /* 
        {
            FILE *tmp_fp;
            tmp_fp = FileOpen("Ins.test", "w");
            for(curr = h_list; curr != NULL; curr = curr->next)
                {
                    ddp = curr->data.ptrvalue;
                    if(strand == Seq_strand_minus)
                        fprintf(tmp_fp, "%ld\t", ddp->starts[1]);
                    else
                        fprintf(tmp_fp, "%ld\t", ddp->starts[1]+ddp->len-1);
                }
            fprintf(tmp_fp, "\n");
            FileClose(tmp_fp);
            system("get_lns < Ins.test >out");
            
            tmp_fp = FileOpen("out", "r");
            fscanf(tmp_fp, "%ld\n", &val);
            if(val == -1)
                {
                    Message(MSG_ERROR, "Fail in get_lns");
                    exit(1);
                }
                /* inconsistent orientation from LNS //
            if((strand == Seq_strand_minus && val == 0) || 
               (strand != Seq_strand_minus && val== 1)) 
                {
                    FileClose(tmp_fp);
                    ValNodeFree(h_list);
                    *ph_list = best_list;
                    return;
                }
            
            
            ValNodeFree(best_list);
            num = (Int4)get_vnp_num(h_list);
            while(fscanf(tmp_fp, "%ld\n", &val) != EOF)
                {
                    if(strand == Seq_strand_minus)
			val = num -1 - val;
                    curr = get_nth_node(h_list, val);
                    if(curr != NULL)
			curr->choice = 1;
                }
            FileClose(tmp_fp);
        }
*/
	delete_bad_node(&h_list, 1);

	if(h_list == NULL)
	{
		Message(MSG_ERROR, "No good node in sort_list_order");
		exit(1);
	}

	*ph_list = h_list;
}



static ValNodePtr check_ddp_order(ValNodePtr ddp_list, Uint1 strand, Int4 first_len, Int4 second_len)
{
	DenseDiagPtr ddp;
	ValNodePtr h_list = NULL, t_list;
	ValNodePtr curr;
	Int4 num, i;
	Boolean reverse;
	BoolPtr first_status, second_status;
	Int4 p_score;
	Int4Ptr score_list;
	Boolean is_end;
	

	if(ddp_list == NULL || ddp_list->next == NULL)
		return ddp_list;
	score_list = load_ddp_scores(ddp_list, &num);
	if(score_list == NULL)
		return ddp_list;

	for(curr = ddp_list; curr != NULL; curr = curr->next)
		curr->choice = 0;
	p_score = score_list[0];
	first_status = (BoolPtr) MemNew((size_t)first_len * sizeof(Boolean));
	second_status = (BoolPtr) MemNew((size_t)second_len * sizeof(Boolean));

	reverse = (strand == Seq_strand_minus);
	h_list = NULL;
	is_end = FALSE;
	while(!is_end)
	{
		t_list = NULL;
		for(i= 0, curr = ddp_list; curr != NULL; curr = curr->next)
		{
 			if(p_score/(score_list[i]) <=2)
			{
				if(curr->choice == 0)
				{
					curr->choice = 1;
					ddp = (DenseDiagPtr) curr->data.ptrvalue;
					if(h_list == NULL || 
						load_ddp_in_order(&h_list, ddp, reverse, FALSE))
					{
						if(!is_ddp_overlap(ddp, first_status, second_status))
							ValNodeAddPointer(&t_list, 0, ddp);
					}
				}
				if(curr->next == NULL)
				{
					p_score = score_list[i];
					is_end = TRUE;
					break;
				}
			}
			else
			{
				p_score = score_list[i];
				break;
			}
			++i;
		}
		if(t_list != NULL)
		{
			sort_list_order(&t_list, strand);
			if(h_list == NULL)
			{
				h_list = t_list;
				for(curr = t_list; curr != NULL; curr = curr->next)
				{
					ddp = (DenseDiagPtr) curr->data.ptrvalue;
					MemSet(first_status+ddp->starts[0], 1, (size_t)(ddp->len) * sizeof(Boolean));
					MemSet(second_status+ddp->starts[1], 1, (size_t)(ddp->len) * sizeof(Boolean));
				}
			}
			else
			{
				for(curr = t_list; curr != NULL; curr = curr->next)
				{
					ddp = (DenseDiagPtr) curr->data.ptrvalue;
					if(load_ddp_in_order(&h_list, ddp, reverse, TRUE))
					{
						MemSet(first_status+ddp->starts[0], 1, (size_t)(ddp->len) * sizeof(Boolean));
						MemSet(second_status+ddp->starts[1], 1, (size_t)(ddp->len) * sizeof(Boolean));
					}
				}
				ValNodeFree(t_list);
			}
		}
		if(is_end)
			break;
	}
			
	if(h_list == NULL)
	{
		Message(MSG_ERROR, "No good node in check_ddp_order");
		exit(1);
	}

	MemFree(score_list);
	MemFree(first_status);
	MemFree(second_status);
	ValNodeFree(ddp_list);
	return h_list;


}

	

/*
*
*	assume the align is the alignment computed from BLAST. Get rid of 
*	all the smaller segments that are out of the main diagnol. the ValNodePtr 
*	returns a list of dense-diags that are on the main diagnol. s_sip specifies 
*	the sequence that will have non-redundant coverage
*/
static ValNodePtr FilterSeqAlign(SeqAlignPtr align, SeqIdPtr s_sip, Uint1Ptr p_strand)
{
	Int2 order = 0, s_order;
	DenseDiagPtr ddp;
	SeqIdPtr sip;
	Uint1 strand;
	BioseqPtr bsp;
	ValNodePtr ddp_list;
	BioseqPtr first_bsp, second_bsp;

	if(s_sip == NULL || align == NULL)
		return NULL;
	if(align->segtype != 1)
		return NULL;
	bsp = BioseqFind(s_sip);
	if(bsp == NULL)
	{
		ErrPostEx(SEV_WARNING, 0, 0, "Fail to find Bioseq for s_sip");
		return NULL;
	}

	ddp_list = NULL;
	s_order = -1;
	for(ddp = (DenseDiagPtr) align->segs; ddp != NULL; ddp = ddp->next)
	{
		ValNodeAddPointer(&ddp_list, 0, ddp);
		if(s_order == -1)
		{
			for(sip =  ddp->id, order = 0; sip != NULL; sip = sip->next, ++order)
			{
				if(SeqIdMatch(s_sip, sip))
				{
					s_order = order;
				}
			}
		}
	}
	if(s_order == -1)
	{
		printf("error in finding the right order\n");
		exit(1);
	}

	/*ddp_list only contains the alignment with the right orientation*/
	/* it also has been sorted according to the order of the scores*/
	strand = get_alignment_orientation(&ddp_list, bsp, s_order);

	if(strand == 0)
		return NULL;

	ddp = (DenseDiagPtr) align->segs;
	first_bsp = BioseqFind(ddp->id);
	second_bsp = BioseqFind(ddp->id->next);
	/* filter_wrong_orientation (align, strand == Seq_strand_plus); */
	filter_repeats(&ddp_list, first_bsp->length, second_bsp->length);

	if(ddp_list == NULL)
	{
		Message(MSG_ERROR, "Nothing leftover after filter_repeats");
		exit(1);
	}
	*p_strand = strand;
	return check_ddp_order(ddp_list, strand, first_bsp->length, second_bsp->length);

}

static Int4 get_dsp_align_len(DenseSegPtr dsp)
{
	Int2 i;
	Int4 len = 0;

	for(i = 0; i<dsp->numseg; ++i)
	{
		if(dsp->starts[2*i] != -1 && dsp->starts[2*i+1] != -1)
			len += dsp->lens[i];
	}

	return len;
}

/*
*
*	convert the Dense-seg to Dense-diag
*	take the end point of aligned segment and make an alignment of 
*	Dense-diag
*
*/

static ScorePtr make_seg_score(Int4 align_len, Int4 seg_len, Int4 total_score)
{
	Int4 score;
	ScorePtr sp;
	ObjectIdPtr oip;

	/* score = (seg_len * total_score)/align_len; */
	score = total_score;
	sp = ScoreNew();
	oip = ObjectIdNew();
	oip->str = StringSave("score");
	sp->id = oip;
	sp->choice = 1;
	sp->value.intvalue = score;
	return sp;
}


static DenseDiagPtr covert_dsp_to_ddp(SeqAlignPtr align, DenseDiagPtr PNTR h_ddp)
{
	DenseSegPtr dsp;
	DenseDiagPtr ddp, prev;
	Int2 i;
	Int4 score;
	Int4 align_len;

	if(align == NULL || align->segtype != 2)
		return NULL;
	dsp = (DenseSegPtr) align->segs;
	if(dsp == NULL || dsp->dim != 2)
		return NULL;
	align_len = get_dsp_align_len(dsp);
	if(align_len < MIN_SEG_SIZE)
		return NULL;

	if(*h_ddp != NULL)
	{
		prev = *h_ddp;
		while(prev->next != NULL)
			prev = prev->next;
	}
	else
		prev = NULL;
	score = get_score_from_list (align->score);
	for(i = 0; i<dsp->numseg; ++i)
	{
		if(dsp->starts[2*i] != -1 && 
			dsp->starts[2*i+1] != -1 && dsp->lens[i] > MIN_SEG_SIZE)
		{
			ddp = DenseDiagNew();
			ddp->dim = 2;
			ddp->id = SeqIdDup(dsp->ids);
			ddp->id->next = SeqIdDup(dsp->ids->next);
			ddp->starts = (Int4Ptr) MemNew((size_t)2 * sizeof(Int4));
			ddp->starts[0] = dsp->starts[2*i];
			ddp->starts[1] = dsp->starts[2*i + 1];
			ddp->len = dsp->lens[i];
			ddp->strands = (Uint1Ptr) MemNew((size_t)2 * sizeof(Uint1));
			if(dsp->strands != NULL)
			{
				ddp->strands[0] = dsp->strands[2*i];
				ddp->strands[1] = dsp->strands[2*i+1];
			}
			ddp->scores = make_seg_score(align_len, dsp->lens[i], score);
			if(*h_ddp == NULL)
				*h_ddp= ddp;
			else
				prev->next = ddp;
			prev = ddp;
		}
	}

	return (*h_ddp);

}

/*
* convert a list of DenseSeg to DenseDiag
*
*/
static SeqAlignPtr ConvertDspAlignList(SeqAlignPtr align)
{
	SeqAlignPtr n_align;
	DenseDiagPtr h_ddp = NULL;


	while(align)
	{
		covert_dsp_to_ddp(align, &h_ddp);
		align = align->next;
	}
	if(h_ddp != NULL)
	{
		n_align = SeqAlignNew();
		n_align->segtype = 1;
		n_align->segs = h_ddp;
		return n_align;
	}
	else
		return NULL;

}

static SeqAlignPtr clean_up_align_list(SeqAlignPtr PNTR h_align, ValNodePtr head)
{
	DenseDiagPtr ddp, n_ddp, p_ddp, c_ddp;
	SeqAlignPtr prev, next;
	Boolean found;
	ValNodePtr curr;
	SeqAlignPtr align;


	prev = NULL;
	align = *h_align;
	while(align)
	{
		next = align->next;

		ddp = (DenseDiagPtr) align->segs;
		p_ddp = NULL;
		while(ddp)
		{
			found = FALSE;
			n_ddp = ddp->next;
			for(curr = head; curr != NULL; curr = curr->next)
			{
				c_ddp = (DenseDiagPtr) curr->data.ptrvalue;
				if(ddp == c_ddp)
				{
					found = TRUE;
					break;
				}
			}
			if(!found)
			{
				if(p_ddp == NULL)
					align->segs = n_ddp;
				else
					p_ddp->next = n_ddp;
				ddp->next = NULL;
				DenseDiagFree(ddp);
			}
			else
				p_ddp = ddp;
			ddp = n_ddp;
		}
		if(align->segs == NULL)
		{
			if(prev == NULL)
				*h_align = next;
			else
				prev->next = next;
			align->next = NULL;
			SeqAlignFree(align);
		}
		else
			prev = align;
		align = next;
	}

	return (*h_align);
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
	

	ddp_align = ConvertDspAlignList(*align);
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



/*
*	return a Seq-align that was made up from the DenseDiag in the 
*	right order. The Seq-align is freshly created
*	the chain of the Seq-align should be for the same Seq-id
*
*/
static SeqAlignPtr FilterDenseSegAlign(SeqAlignPtr align, SeqIdPtr sip)
{
	SeqAlignPtr h_align;
	Uint1 strand;
	ValNodePtr head = NULL;

	h_align = ConvertDspAlignList(align);
	if(h_align == NULL)
		return NULL;
	head = FilterSeqAlign(h_align, sip, &strand);

	if(head == NULL)
	{
		SeqAlignFree(h_align);
		return NULL;
	}
	clean_up_align_list(&h_align, head);
	ValNodeFree(head);
	return h_align;
}



/*
*
*	from the align, figure out the maximum and minimum range on the two sequences. 
*   the two Seq-locs can then be sent for the seond path algorithm to compute 
*	sequence alignment
*	align is normally computed from a heuristic method, such as BLAST
*	return TRUE for success and FLASE for failure
*/
static Boolean FigureAlignRange(SeqAlignPtr align, SeqLocPtr m_loc, SeqLocPtr s_loc)
{
	DenseDiagPtr ddp;
	ValNodePtr head, vnp;
	SeqIntPtr sint;
	Int4 s_max =-1 , s_min = -1;
	Int4 m_max = -1, m_min = -1;
	Uint1 strand;
	Int2 s_order, order;
	SeqIdPtr sip;
	Uint1 m_strand, s_strand;



	m_strand = 0;
	s_strand = 0;
	head = FilterSeqAlign(align, SeqLocId(s_loc), &strand);
	if(head == NULL)
		return FALSE;

	ddp = (DenseDiagPtr) head->data.ptrvalue;
	s_order = -1;
	order = 0;
	for(sip = ddp->id; sip != NULL; sip = sip->next, ++order)
	{
		if(SeqIdForSameBioseq(sip, SeqLocId(s_loc)))
		{
			s_order = order;
			break;
		}
	}
	if(s_order == -1)
	{
		ValNodeFree(head);
		return FALSE;
	}
	s_min = ddp->starts[s_order];
	if(strand == Seq_strand_minus)
	{
		m_max = ddp->starts[1-s_order] + ddp->len -1;
	}
	else
		m_min = ddp->starts[1-s_order];
	vnp = head;
	while(vnp->next != NULL)
	{
		ddp = (DenseDiagPtr) vnp->data.ptrvalue;
		/* printf("%ld, %ld, len = %ld\n", ddp->starts[0], ddp->starts[1], ddp->len); */
		vnp = vnp->next;
	}
	ddp = (DenseDiagPtr) vnp->data.ptrvalue;
	/* printf("%ld, %ld, len = %ld\n", ddp->starts[0], ddp->starts[1], ddp->len); */ 
	s_max = ddp->starts[s_order] + ddp->len -1;
	if(strand == Seq_strand_minus)
	{
		m_min = ddp->starts[1-s_order];
	}
	else
		m_max = ddp->starts[1-s_order] + ddp->len -1;

	if(m_min > m_max)
	{

		Message(MSG_ERROR, "m_min = %ld, m_max = %ld", m_min, m_max);
		ValNodeFree(head);
		return FALSE;
	}
	ValNodeFree(head);
	/* printf("m_min = %ld m_max = %ld, s_min = %ld, s_max = %ld \n", m_min, m_max, s_min, s_max); */
	

	if(s_min == -1 || s_max == -1 || m_min == -1 || m_max == -1)
		return FALSE;

	m_strand = Seq_strand_plus;
	s_strand = strand;

	/* printf("s_min = %ld, s_max = %ld, s_strand = %d\n", s_min, s_max, s_strand);
	printf("m_min = %ld, m_max = %ld, m_strand = %d\n", m_min, m_max, m_strand); */


	sint = (SeqIntPtr) m_loc->data.ptrvalue;
	sint->from = m_min;
	sint->to = m_max;
	sint->strand = m_strand;

	sint = (SeqIntPtr) s_loc->data.ptrvalue;
	sint->from = s_min;
	sint->to = s_max;
	sint->strand = s_strand;

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
					start_2+dsp->lens[i]-1, dsp->strands[0]);

				return TRUE;
			}
		}
	}

	return FALSE;
}

static Boolean is_loc_overlap(SeqLocPtr loc_1, SeqLocPtr loc_2, Int4Ptr min_val, Int4Ptr max_val, Int4 gap_len, Int4 max_align_size, Int4Ptr p_gap_val)
{
	Int4 start_1, stop_1, start_2, stop_2;

	start_1 = SeqLocStart(loc_1);
	stop_1 = SeqLocStop(loc_1);
	start_2 = SeqLocStart(loc_2);
	stop_2 = SeqLocStop(loc_2);

	/*strict overlap*/
	if(!(stop_1 < start_2 || start_1 > stop_2))
	{
		*min_val = MAX(start_1, start_2);
		*max_val = MIN(stop_1, stop_2);
		return TRUE;
	}

	if(gap_len == 0)
		return FALSE;
	if(stop_2 > stop_1)
	{
		/* if(start_2 - stop_1 <= gap_len)  */
		if(start_2 - stop_1 <= max_align_size)
		{
			*min_val = stop_1;
			*max_val = start_2;
			if(p_gap_val != NULL)
				*p_gap_val = start_2 - stop_1;
			return TRUE;
		}
		else
			return FALSE;
	}
	else
	{
		/* if(start_1 - stop_2 <= gap_len) */
		if(start_1 - stop_2 <= max_align_size)
		{
			*min_val = stop_2;
			*max_val = start_1;
			if(p_gap_val != NULL)
				*p_gap_val = start_1 - stop_2;
			return TRUE;
		}
		else
			return FALSE;
	}

}

	
/*
*
*	ppos is the position of the known sequence, trying to find the matching sequence 
*	position
*/
static Int4 find_matching_position(DenseSegPtr dsp, Int4 pos, Uint1 pos_order)
{
	Int2 i;
	Int4 start_1, start_2;
	Int4 offset;


	for(i = 0; i<dsp->numseg; ++i)
	{
		start_1 = dsp->starts[2*i + pos_order];
		start_2 = dsp->starts[2*i + 1 - pos_order];

		if(start_1 != -1 && start_2 != -1)
		{
			if(pos >= start_1 && pos <= start_1 + dsp->lens[i] -1)
			{
				if(dsp->strands[pos_order] == Seq_strand_minus)
				{
					offset = start_1 + dsp->lens[i] - 1 - pos;
				}
				else
					offset = pos - start_1;

				if(dsp->strands[1 - pos_order] == Seq_strand_minus)
					return (start_2 + dsp->lens[i] -1 - offset);
				else
					return (start_2 + offset);
			}
		}
	}


	return -1;
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

static Boolean select_overlap_loc(SeqLocPtr loc_1_1, SeqLocPtr loc_2_1, SeqLocPtr loc_1_2, SeqLocPtr loc_2_2, Uint1Ptr p_order)
{
	Int4 start_1, stop_1, start_2, stop_2;
	Int4 overlap_1, overlap_2;

	overlap_1 = 0;
	overlap_2 = 0;

	start_1 = SeqLocStart(loc_1_1);
	stop_1 = SeqLocStop(loc_1_1);

	start_2 = SeqLocStart(loc_2_1);
	stop_2 = SeqLocStop(loc_2_1);

	if(!(start_2 > stop_1 || stop_2 < start_1))
	{
		if(stop_2 >= stop_1)
			overlap_1 = MAX(0, stop_1 - start_2 + 1);
		else
			overlap_1 = MAX(0, stop_2 - start_1 + 1);
	}

	start_1 = SeqLocStart(loc_1_2);
	stop_1 = SeqLocStop(loc_1_2);

	start_2 = SeqLocStart(loc_2_2);
	stop_2 = SeqLocStop(loc_2_2);

	if(!(start_2 > stop_1 || stop_2 < start_1))
	{
		if(stop_2 >= stop_1)
			overlap_2 = MAX(0, stop_1 - start_2 + 1);
		else
			overlap_2 = MAX(0, stop_2 - start_1 + 1);
	}

	if(overlap_1 == 0 && overlap_2 == 0)
		return FALSE;
	else
	{
		if(overlap_1 >= overlap_2)
			*p_order = 0;
		else
			*p_order = 1;
		return TRUE;
	}
}
					

/*
*
*	align_1 and align_2 are sorted by the ascending order in the first 
*	sequence in the alignment. The first sequence needs to aligned 
*	in the plus strand. They have to be alignment on the same 
*	diagnol
*	return NULL if any merge fails
*
*/
static SeqAlignPtr MergeTwoDspBySIM4(SeqAlignPtr align_1, SeqAlignPtr align_2)
{
	SeqLocPtr loc_1_1, loc_1_2, loc_2_1, loc_2_2;
	Int4 start, stop;
	Int4 match_start, match_stop;
	Boolean found = FALSE;
	Uint1 order;
	BioseqPtr bsp;
	SeqAlignPtr align, t_align;
	DenseSegPtr dsp_1, dsp_2;
	BioseqPtr bsp_1, bsp_2;
	Int4 max_align_size;
	Uint1 strand;
	Int4 gap_1, gap_2;
	Int4 start_1, stop_1, start_2, stop_2;
	SeqIdPtr sip;
	Boolean reverse = FALSE;
	Boolean has_overlap;

	if(align_1 == NULL || align_2 == NULL)
		return NULL;

	if(align_1->segtype != 2 || align_2->segtype != 2)
		return NULL;

	dsp_1 = (DenseSegPtr) align_1->segs;
	if(dsp_1 && dsp_1->strands &&  dsp_1->strands[0] == Seq_strand_minus)
	{
		reverse_alignment(align_1, 0);
		reverse = TRUE;
	}
	dsp_2 = (DenseSegPtr) align_2->segs;
	if(dsp_2 && dsp_2->strands &&  dsp_2->strands[0] == Seq_strand_minus)
	{
		reverse_alignment(align_2, 0);
		reverse = TRUE;
	}

	max_align_size = MIN(get_align_length(align_1), 
		get_align_length(align_2));

	if(dsp_1 == NULL || dsp_2 == NULL)
		return NULL;
	/* if(dsp_2->starts[1] == 40737)
		printf("stop here\n"); */

	get_align_ends(align_1, dsp_1->ids, &start_1, &stop_1, &strand); 
	get_align_ends(align_2, dsp_1->ids, &start_2, &stop_2, &strand); 
	if(start_2 >= start_1 && stop_2 <= stop_1)
		return NULL;


	get_align_ends(align_1, dsp_1->ids->next, &start_1, &stop_1, &strand); 
	get_align_ends(align_2, dsp_1->ids->next, &start_2, &stop_2, &strand); 
	if(start_2 >= start_1 && stop_2 <= stop_1)
		return NULL;

	get_align_ends(align_1, dsp_1->ids, &start, &stop, &strand); 
	loc_1_1 = SeqLocIntNew(start, stop, strand, dsp_1->ids);
	get_align_ends(align_1, dsp_1->ids->next, &start, &stop, &strand); 
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
			}
			if(is_loc_overlap(loc_1_2, loc_2_2, &start_2, &stop_2, MIN_SEG_SIZE*10, max_align_size, &gap_2))
			{
				if(!found || gap_2 < gap_1)
				{
					order = 1;
					found = TRUE;
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

	bsp = NULL;
	if(order == 0)
		bsp = BioseqLockById(dsp_1->ids);
	else
		bsp = BioseqLockById(dsp_1->ids->next);
	if(bsp == NULL)
		return NULL;

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
	get_align_ends(align_1, sip, &start_1, &stop_1, &strand); 
	get_align_ends(align_2, sip, &start_2, &stop_2, &strand); 

	if(dsp_1->strands[order] != Seq_strand_minus)
	{
		if(has_overlap)
		{
			stop = MIN(stop_1 + 100, stop_2);
			start = MAX(start_2 -100, start_1);
		}
		else
		{
			start -= 100;
			start = MAX(start, start_1);
			stop += 100;
			stop = MIN(stop, stop_2);
		}
			
	}
	else
	{
		if(has_overlap)
		{
			start = MAX(start_1 -100, start_2);
			stop = MIN(stop_2+100, stop_1);
		}
		else
		{
			start -= 100;
			start = MAX(start, start_2);
			stop += 100;
			stop = MIN(stop, stop_1);
		}
	}
	BioseqUnlock(bsp);


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

	swap(&match_start, &match_stop);
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
		reverse_alignment(align, 1);
	return align;
}

		
				
/*###############################################################
*
*	functions for merge the two alignments
*
################################################################*/


static Int4 get_align_len_in_overlap (SeqAlignPtr align, Int4 from, Int4 to, Int2 order)
{
	DenseSegPtr dsp;
	Int2 i;
	Int4 start, stop;
	Int4 align_len = 0;

	dsp = (DenseSegPtr) align->segs;
	for(i = 0; i<dsp->numseg; ++i)
	{
		if(dsp->starts[2*i] != -1 && dsp->starts[2*i+1] != -1)
		{
			start = dsp->starts[2*i + order];
			stop = dsp->starts[2*i + order] + dsp->lens[i] -1;
			if(!(start > to || stop < from))
			{
				start = MAX(start, from);
				stop = MIN(stop, to);

				if(stop >= start)
					align_len += (stop - start + 1);
			}
		}
	}

	return align_len;
}


typedef struct align_overlap {
	SeqAlignPtr align;
	Int4 overlap_len;
	Int4 align_len;
}AlignOverlap, PNTR AlignOverlapPtr;

static int LIBCALLBACK AlignOverlapCompProc (VoidPtr ptr1, VoidPtr ptr2)
{
  ValNodePtr vnp1, vnp2;
  AlignOverlapPtr sop1, sop2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
	sop1 = (AlignOverlapPtr) vnp1->data.ptrvalue;
	sop2 = (AlignOverlapPtr) vnp2->data.ptrvalue;
	if(sop1 != NULL && sop2 != NULL)
	{
		if(sop1->overlap_len > sop2->overlap_len)
			return -1;
		else if(sop1->overlap_len < sop2->overlap_len)
			return 1;
		else if(sop1->align_len > sop2->align_len)
			return -1;
		else if(sop1->align_len < sop2->align_len)
			return 1;
		else
			return 0;
	}
    }
  }

  return 0;
}
 

static ValNodePtr build_overlap_list(SeqAlignPtr align, Int4 from, Int4 to, Int2 order)
{
	Int4 t_len;
	ValNodePtr list = NULL;
	AlignOverlapPtr aop;

	/*sort the alignment by score*/
	while(align)
	{
		/*check if the alignment is significant*/
		t_len = get_align_len_in_overlap (align, from, to, order);
		if(t_len > MIN_DIAG_LEN)
		{
			aop = (AlignOverlapPtr) MemNew(sizeof(AlignOverlap));
			aop->overlap_len = t_len;
			aop->align = align;
			aop->align_len = get_align_length(align);
			ValNodeAddPointer(&list, 0, aop);
		}
		align = align->next;
	}

	if(list != NULL)
		return SortValNode(list, AlignOverlapCompProc);
	else
		return NULL;

}

static SeqAnnotPtr open_annot(CharPtr f_name)
{
	SeqAnnotPtr annot;
	AsnIoPtr aip;

	aip = AsnIoOpen(f_name, "r");
	annot = SeqAnnotAsnRead(aip, NULL);
	AsnIoClose(aip);
	return annot;
}

typedef struct segdata {
	Int2 seg_order;
	Int4 overlap_len;
}SegOrder, PNTR SegOrderPtr;

/*
*	find the longest diagnol in the overlapping region
*
*/
static int LIBCALLBACK SegOrderCompProc (VoidPtr ptr1, VoidPtr ptr2)
{
  ValNodePtr vnp1, vnp2;
  SegOrderPtr sop1, sop2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
	sop1 = (SegOrderPtr) vnp1->data.ptrvalue;
	sop2 = (SegOrderPtr) vnp2->data.ptrvalue;
	if(sop1 != NULL && sop2 != NULL)
	{
		if(sop1->overlap_len > sop2->overlap_len)
			return -1;
		else if(sop1->overlap_len < sop2->overlap_len)
			return 1;
		else
			return 0;
	}
    }
  }

  return 0;
}
 
static ValNodePtr find_overlap_diagnol(SeqAlignPtr align, Int4 from, Int4 to, Int2 order)
{
	Int2 i;
	DenseSegPtr dsp;
	Int4 start, stop;
	SegOrderPtr sop;
	ValNodePtr seg_list;

	seg_list = NULL;
	dsp = (DenseSegPtr) align->segs;
	for(i = 0; i<dsp->numseg; ++i)
	{
		if(dsp->starts[2*i] != -1 && dsp->starts[2*i+1] != -1)
		{
			start = dsp->starts[2*i+order];
			stop = start + dsp->lens[i] -1;
			if(!(start > to || stop < from))
			{
				
				start = MAX(start, from);
				stop = MIN(stop, to);
				sop = (SegOrderPtr) MemNew(sizeof(SegOrder));
				sop->seg_order = i;
				sop->overlap_len = stop - start + 1;
				ValNodeAddPointer(&seg_list, 0, sop);
			}
		}
	}

	if(seg_list != NULL)
		return SortValNode(seg_list, SegOrderCompProc);
	else
		return NULL;
}

static Boolean is_dsp_same(DenseSegPtr dsp_1, DenseSegPtr dsp_2)
{
	if(dsp_1 == NULL || dsp_2 == NULL)
		return FALSE;
	if(dsp_1->ids == NULL || dsp_2->ids == NULL)
		return FALSE;
	if(!SeqIdMatch(dsp_1->ids, dsp_2->ids))
		return FALSE;

	if(dsp_1->ids->next == NULL || dsp_2->ids->next == NULL)
		return FALSE;
	if(!SeqIdMatch(dsp_1->ids->next, dsp_2->ids->next))
		return FALSE;

	if(dsp_1->strands == NULL || dsp_2->strands == NULL)
		return FALSE;

	if(!is_same_orientation(dsp_1->strands[0], dsp_2->strands[0]))
		return FALSE;

	if(!is_same_orientation(dsp_1->strands[1], dsp_2->strands[1]))
		return FALSE;

	return TRUE;
}

static Int4 get_score_value(ScorePtr sp)
{
        ObjectIdPtr oip;

        while(sp)
        {
                if(sp->id)
                {
                        oip = sp->id;
                        if(oip && oip->str && StringCmp(oip->str, "score") == 0)
                        {
                                if(sp->choice == 1)
                                        return sp->value.intvalue;
                        }
                }
                sp = sp->next;
        }

        return -1;
}

static void load_big_score (SeqAlignPtr align_1, SeqAlignPtr align_2)
{
	Int4 val_1, val_2;

	val_1 = get_score_value(align_1->score);
	val_2 = get_score_value(align_2->score);

	if(val_1 < val_2)
	{
		ScoreSetFree(align_1->score);
		align_1->score = align_2->score;
		align_2->score  = NULL;
	}
}



static Boolean merge_two_align(SeqAlignPtr align_1, SeqAlignPtr align_2, Int4 from, Int4 to, Int2 order)
{
	ValNodePtr seg_list_1, seg_list_2;
	SegOrderPtr sop_1, sop_2;
	ValNodePtr curr_1, curr_2;
	Boolean found;
	Int2 seg_order_1, seg_order_2;
	DenseSegPtr dsp_1, dsp_2;
	Int4 start_1, stop_1, start_2, stop_2;
	Int4 overlap_1, overlap_2;
	Int4 extend_len;
	Int2 n_seg_num, i, j;
	Int4Ptr starts, lens;
	Uint1Ptr strands;
	Int4 m_start_1, s_start_1, len_1;
	Int4 m_start_2, s_start_2, len_2;
	Int4 offset_1;

	dsp_1 = (DenseSegPtr) align_1->segs;
	dsp_2 = (DenseSegPtr) align_2->segs;

	if(!is_dsp_same(dsp_1, dsp_2))
		return FALSE;
	seg_list_1 = find_overlap_diagnol(align_1, from, to, order);
	if(seg_list_1 == NULL)
		return FALSE;
	seg_list_2 = find_overlap_diagnol(align_2, from, to, order);
	if(seg_list_2 == NULL)
	{
		ValNodeFreeData(seg_list_1);
		return FALSE;
	}

	found = FALSE;
	for(curr_1 = seg_list_1; curr_1 != NULL && !found; curr_1 = curr_1->next)
	{
		sop_1 = (SegOrderPtr) curr_1->data.ptrvalue;
		seg_order_1 = sop_1->seg_order;
		m_start_1 = dsp_1->starts[2*seg_order_1 + order];
		s_start_1 = dsp_1->starts[2*seg_order_1 + 1 - order];
		len_1 = dsp_1->lens[seg_order_1];
		for(curr_2 = seg_list_2; curr_2 != NULL && !found; curr_2 = curr_2->next)
		{
			sop_2 = (SegOrderPtr) curr_2->data.ptrvalue;
			seg_order_2 = sop_2->seg_order;

			m_start_2 = dsp_2->starts[2*seg_order_2 + order];
			s_start_2 = dsp_2->starts[2*seg_order_2 + 1 - order];
			len_2 = dsp_2->lens[seg_order_2];
		

			found = TRUE;
			start_1 = dsp_1->starts[2*seg_order_1 + order];
			stop_1 = start_1 + dsp_1->lens[seg_order_1] -1;
	
			start_2 = dsp_2->starts[2*seg_order_2 + order];
			stop_2 = start_2 + dsp_2->lens[seg_order_2] -1;

			/*check if the overlapping region on the sequence overlaps??*/
			if(start_1 > stop_2 || stop_1 < start_2)
				found = FALSE;

			offset_1 = 0;
			/*trim the tail of the first alignment and the head of the second*/
			if(found)
			{
				if(dsp_1->strands[order] == Seq_strand_minus)
				{
					/*trim the head of the second alignment*/
					/*    <====-------===================== first
					      <================================  second
					trim  ^^^^^^^^^^^^XXXXXXXXXXXXXXXXXXXXX  overlap
					*/
					if(stop_2 > stop_1)
					{
						offset_1 = stop_2 - stop_1;
						dsp_2->lens[seg_order_2] -= offset_1;
						if(dsp_2->strands[1-order] != Seq_strand_minus)
							dsp_2->starts[2*seg_order_2 + 1 - order] += offset_1;
					
						stop_2 = stop_1;
					}
					/*trim the tail of the first alignment*/
					/*
						                           ^^^^^^^^^^ trim
					         <=================================== first
					            <=======================------=========  second
						     XXXXXXXXXXXXXXXXXXXXXXX	overlap
					*/
					if(start_1 < start_2)
					{
						offset_1 = start_2 - start_1;
						dsp_1->starts[2*seg_order_1 + order] += offset_1;
						dsp_1->lens[seg_order_1] -= offset_1;
						if(dsp_1->strands[1-order] == Seq_strand_minus)
							dsp_1->starts[2*seg_order_1 + 1 - order] += offset_1;
					
						start_1 = start_2;
					}
					overlap_1 = start_1 - stop_2;
				}
				else
				{
					/*trim the head of the second alignment*/
					/*
						====----====================> first
						================================ second
						^^^^^^^^XXXXXXXXXXXXXXXXXXXX 
						trim     overlap
					*/
					if(start_2 < start_1)
					{
						offset_1 = start_1 - start_2;
						dsp_2->starts[2*seg_order_2 + order] += offset_1;
						if(dsp_2->strands[1-order] != Seq_strand_minus)
					
							dsp_2->starts[2*seg_order_2 + 1 - order] += offset_1;
						dsp_2->lens[seg_order_2] -= offset_1;
						start_2 = start_1;
					}
					/*trim the tail of the first alignment*/
					/*
						  overlap               trim
						XXXXXXXXXXXXXXXXXXXXXXX^^^^^^^^^^
						=================================> first
						========================------============= second
					*/
					
					if(stop_1 > stop_2)
					{
						offset_1 = stop_1 - stop_2;
						dsp_1->lens[seg_order_1] -= offset_1;
						if(dsp_1->strands[1-order] == Seq_strand_minus)
					
							dsp_1->starts[2*seg_order_1 + 1 - order] += offset_1;
						stop_1 = stop_2;
					}
					overlap_1 = stop_1 - start_2;
				}
				if(found && overlap_1 <=0)
					found = FALSE;
			}


			/*check the second sequence*/
			if(found)
			{
				start_1 = dsp_1->starts[2*seg_order_1 + 1 - order];
				stop_1 = start_1 + dsp_1->lens[seg_order_1] -1;
	
				start_2 = dsp_2->starts[2*seg_order_2 + 1 - order];
				stop_2 = start_2 + dsp_2->lens[seg_order_2] -1;

				/*check if the overlapping region on the sequence overlaps??*/
				if(start_1 > stop_2 || stop_1 < start_2)
					found = FALSE;
			}
			if(found)
			{
				if(dsp_1->strands[ 1 - order] == Seq_strand_minus)
				{
					overlap_2 = stop_2 - start_1;
				}
				else
				{
					overlap_2 = stop_1 - start_2;
				}
				if(found && overlap_1 <=0)
					found = FALSE;
			}

			if(found && overlap_2 != overlap_1)
				found = FALSE;

			if(found)
			{
				if(dsp_1->strands[0] == Seq_strand_minus)
				{
					
					stop_1 = dsp_1->starts[2*seg_order_1 + 1 - order] + 
						dsp_1->lens[seg_order_1] -1;
					stop_2 = dsp_2->starts[2*seg_order_2 + 1 - order] + 
						dsp_2->lens[seg_order_2] -1;
					extend_len = stop_1 - stop_2;
				}
			
				else
				{
					start_1 = dsp_1->starts[2*seg_order_1 + order];
					start_2 = dsp_2->starts[2*seg_order_2 + order];
					extend_len = start_2 - start_1;
				}
			
				/* if(dsp_1->strands[order] == Seq_strand_minus)
				{
					
					start_1 = dsp_1->starts[2*seg_order_1 + order];
					start_2 = dsp_2->starts[2*seg_order_2 + order];
					extend_len = MAX(0, start_1 - start_2);
				}
				else
				{
					stop_1 = dsp_1->starts[2*seg_order_1 + order] + 
						dsp_1->lens[seg_order_1] -1;
					stop_2 = dsp_2->starts[2*seg_order_2 + order] + 
						dsp_2->lens[seg_order_2] -1;
					extend_len = MAX(0, stop_2 - stop_1);
				} */
				if(dsp_2->strands[order] != Seq_strand_minus)
				{
					dsp_2->starts[2*seg_order_2 + order] -= extend_len;
				}
				if(dsp_2->strands[1-order] != Seq_strand_minus)
				{
					dsp_2->starts[2*seg_order_2 + 1 - order] -= extend_len;
				}

				dsp_2->lens[seg_order_2] += extend_len;

				n_seg_num = seg_order_1 + (dsp_2->numseg - seg_order_2);
				if(n_seg_num < 1)
					found = FALSE;
			}
			if(found)
			{

				starts = (Int4Ptr) MemNew((size_t)(2*n_seg_num) * sizeof(Int4));
				strands = (Uint1Ptr) MemNew((size_t)(2*n_seg_num) * sizeof(Uint1));
				lens = (Int4Ptr) MemNew((size_t) n_seg_num * sizeof(Int4));

				j = 0;
				if(seg_order_1 > 0)
				{
					for(i = 0; i<seg_order_1; ++i)
					{
						starts[2*j] = dsp_1->starts[2*i];
						starts[2*j+1] = dsp_1->starts[2*i+1];
						strands[2*j] = dsp_1->strands[2*i];
						strands[2*j+1] = dsp_1->strands[2*i+1];
			
						lens[j] = dsp_1->lens[i];
						++j;
					}
				}
			
				for(i = seg_order_2; i<dsp_2->numseg; ++i)
				{
					starts[2*j] = dsp_2->starts[2*i];
					starts[2*j+1] = dsp_2->starts[2*i+1];
					strands[2*j] = dsp_2->strands[2*i];
					strands[2*j+1] = dsp_2->strands[2*i+1];
					lens[j] = dsp_2->lens[i];
					++j;
				}
			
				dsp_1->numseg = n_seg_num;
				MemFree(dsp_1->starts);
				MemFree(dsp_1->strands);
				MemFree(dsp_1->lens);
			
				dsp_1->starts = starts;
				dsp_1->strands = strands;
				dsp_1->lens = lens;
				load_big_score (align_1, align_2);
				break;
			}
			else
			{
				dsp_1->starts[2*seg_order_1 + order] = m_start_1;
				dsp_1->starts[2*seg_order_1 + 1 - order] = s_start_1;
				dsp_1->lens[seg_order_1] = len_1;

				dsp_2->starts[2*seg_order_2 + order] = m_start_2;
				dsp_2->starts[2*seg_order_2 + 1 - order] = s_start_2;
				dsp_2->lens[seg_order_2] = len_2;
			}
		}
	}

	ValNodeFreeData(seg_list_1);
	ValNodeFreeData(seg_list_2);
	return found;
			
}

static SeqAlignPtr extract_align_from_list(SeqAlignPtr PNTR list, SeqAlignPtr align)
{
	SeqAlignPtr curr, prev;

	prev = NULL;
	curr = *list;
	while(curr)
	{
		if(curr == align)
		{
			if(prev == NULL)
				*list = curr->next;
			else
				prev->next = curr->next;
			curr->next = NULL;
			return curr;
		}
		prev = curr;
		curr = curr->next;
	}
	return NULL;
}

static void reverse_alignment(SeqAlignPtr align, Uint1 order)
{
	DenseSegPtr dsp;
	Int4Ptr starts, lens;
	Int2 numseg, i, j;

	if(align->segtype != 2)
		return;
	dsp = (DenseSegPtr) align->segs;
	if(dsp == NULL || dsp->strands == NULL || dsp->strands[order] != Seq_strand_minus)
		return;

	numseg = dsp->numseg;
	starts = (Int4Ptr) MemNew((size_t)2*numseg * sizeof(Int4));
	lens = (Int4Ptr) MemNew((size_t)numseg * sizeof(Int4));

	for(i = numseg-1, j = 0; i>=0; --i, ++j)
	{
		starts[2*j] = dsp->starts[2*i];
		starts[2*j+1] = dsp->starts[2*i+1];
		lens[j] = dsp->lens[i];
		dsp->strands[2*i] = 3 - dsp->strands[2*i];
		dsp->strands[2*i + 1] = 3 - dsp->strands[2*i + 1];
	}

	MemFree(dsp->starts);
	dsp->starts = starts;
	MemFree(dsp->lens);
	dsp->lens = lens;
}

NLM_EXTERN Boolean MergeTwoAlignList (SeqAlignPtr h_align_1, SeqAlignPtr PNTR p_align_2, Int4 from, Int4 to, Int2 order)
{

	ValNodePtr list_1, list_2;
	AlignOverlapPtr aop_1, aop_2;
	ValNodePtr curr_1, curr_2;
	SeqAlignPtr align_1, align_2;
	Boolean found;
	SeqAlignPtr h_align_2;
	Boolean retval = FALSE;

	h_align_2 = *p_align_2;
	if(p_align_2 == NULL || h_align_1 == NULL || (h_align_2 = *p_align_2)== NULL)
		return FALSE;

	list_1 = build_overlap_list(h_align_1, from, to, order);
	if(list_1 == NULL)
		return FALSE;
	list_2 = build_overlap_list(h_align_2, from, to, order);
	if(list_2 == NULL)
	{
		ValNodeFreeData(list_1);
		return FALSE;
	}

	for(curr_1 = list_1; curr_1 != NULL; curr_1 = curr_1->next)
	{
		aop_1 = (AlignOverlapPtr) curr_1->data.ptrvalue;
		align_1 = aop_1->align;
		reverse_alignment(align_1, order);

		found = FALSE;
		for(curr_2 = list_2; !found && curr_2 != NULL; curr_2 = curr_2->next)
		{
			aop_2 = (AlignOverlapPtr) curr_2->data.ptrvalue;
			if(aop_2->align != NULL)
			{
				align_2 = aop_2->align;
				reverse_alignment(align_2, order);
				if(merge_two_align(align_1, align_2, from, to, 0))
				{
					retval = TRUE;
					aop_2->align = NULL;
					if(extract_align_from_list(p_align_2, align_2) != NULL)
						SeqAlignFree(align_2);
					found = TRUE;
				}
			}
		}
	}

	ValNodeFreeData(list_1);
	ValNodeFreeData(list_2);

	return retval;
}

static void link_this_align(SeqAlignPtr PNTR h_align, SeqAlignPtr align)
{
	SeqAlignPtr curr;

	if(*h_align == NULL)
		*h_align = align;
	else
	{
		curr = *h_align;
		while(curr->next != NULL)
			curr = curr->next;
		curr->next = align;
	}
}

					

/*
*
*	functions related to reduce the redundant level of the view
*/
static Int4 get_align_length(SeqAlignPtr align)
{
	Int4 len = 0;
	DenseSegPtr dsp;
	DenseDiagPtr ddp;
	/* StdSegPtr ssp; */
	Int2 i;

	len = 0;
	switch(align->segtype)
	{
	case 1:
		ddp = (DenseDiagPtr) align->segs;
		while(ddp)
		{
			len += ddp->len;
			ddp = ddp->next;
		}
		break;
	case 2:
		dsp = (DenseSegPtr) align->segs;
		for(i = 0; i<dsp->numseg; ++i)
		{
			if(dsp->starts[2*i] != -1 && dsp->starts[2*i+1] != -1)
				len += dsp->lens[i];
		}
		break;
	default:
		break;
	}

	return len;
}

typedef struct align_info {
	SeqAlignPtr align;
	Int4 align_len;
}AlignInfo, PNTR AlignInfoPtr;


static int LIBCALLBACK CompareAlignInfoProc(VoidPtr ptr1, VoidPtr ptr2)
{
  ValNodePtr   vnp1;
  ValNodePtr   vnp2;
  AlignInfoPtr aip_1, aip_2;

  if (ptr1 != NULL && ptr2 != NULL) {
    vnp1 = *((ValNodePtr PNTR) ptr1);
    vnp2 = *((ValNodePtr PNTR) ptr2);
    if (vnp1 != NULL && vnp2 != NULL) {
      aip_1 = (AlignInfoPtr) vnp1->data.ptrvalue;
      aip_2 = (AlignInfoPtr) vnp2->data.ptrvalue;
	  if(aip_1->align_len > aip_2->align_len)
		  return -1;
	  else if(aip_1->align_len < aip_2->align_len)
		  return 1;
	  else return 0;
	}
  }
  return 0;
}

static SeqAlignPtr sort_align_by_length(SeqAlignPtr align)
{
	AlignInfoPtr aip;
	ValNodePtr vnp = NULL, curr;
	SeqAlignPtr h_align, prev;

	while(align)
	{
		aip = (AlignInfoPtr) MemNew(sizeof(AlignInfo));
		aip->align_len = get_align_length(align);
		aip->align = align;
		ValNodeAddPointer(&vnp, 0, aip);
		align = align->next;
	}
	vnp = SortValNode (vnp, CompareAlignInfoProc);
	h_align = NULL;
	prev = NULL;
	for(curr = vnp; curr != NULL; curr = curr->next)
	{
		aip = (AlignInfoPtr) curr->data.ptrvalue;
		if(prev == NULL)
			h_align = aip->align;
		else
			prev->next = aip->align;
		prev = aip->align;
	}
	if(prev != NULL)
		prev->next = NULL;
	ValNodeFreeData(vnp);

	return h_align;
}

	
static Int4 get_align_mid_point(SeqAlignPtr align, Uint1 order)
{
	DenseSegPtr dsp;
	Int4 start, stop;
	Uint1 strand, val;
	SeqIdPtr sip;

	if(align->segtype == 2)
	{
		dsp = (DenseSegPtr) align->segs;
		val = 0;
		for(sip = dsp->ids; sip != NULL; sip = sip->next)
		{
			if(val == order)
			{
				get_align_ends(align, sip, &start, &stop, &strand); 
				return (start + stop)/2;
			}
			++val;
		}
	}

	return -1;
}
		
static SeqAlignPtr sort_align_by_region(SeqAlignPtr align, Uint1 order)
{
	AlignInfoPtr aip;
	ValNodePtr vnp = NULL, curr;
	SeqAlignPtr h_align;

	while(align)
	{
		aip = (AlignInfoPtr) MemNew(sizeof(AlignInfo));
		aip->align_len = get_align_mid_point(align, order);
		aip->align = align;
		ValNodeAddPointer(&vnp, 0, aip);
		align = align->next;
	}
	vnp = SortValNode (vnp, CompareAlignInfoProc);
	h_align = NULL;
	for(curr = vnp; curr != NULL; curr = curr->next)
	{
		aip = (AlignInfoPtr) curr->data.ptrvalue;
		aip->align->next = NULL;
		if(h_align == NULL)
			h_align = aip->align;
		else
		{
			aip->align->next = h_align;
			h_align = aip->align;
		}
	}
	ValNodeFreeData(vnp);

	return h_align;
}


static Int4 get_align_list_len(SeqAlignPtr align)
{
	Int4 len = 0;
	
	while(align)
	{
		len += get_align_length(align);
		align = align->next;
	}

	return len;
}


static SeqAlignPtr MergeToOneAlignment(SeqAlignPtr PNTR palign)
{
	SeqAlignPtr next, prev, align;
	SeqAlignPtr merge_align;
	SeqAlignPtr curr;

	if(palign == NULL || *palign == NULL)
		return NULL;
	align = *palign;
	if(align->segtype != 2)
		return NULL;
	*palign = sort_align_by_region(align, 0);

	prev = NULL;
	align = *palign;
	/* align->next->next->next->next->next->next->next->next->next = NULL; */
	while(align)
	{
		next = align->next;
		if(next != NULL)
		{
			merge_align = NULL;
			curr = *palign;
			/* if(next->next == NULL)
				printf("stop here\n"); */
			merge_align = MergeTwoDspBySIM4(align, next);
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

	if(*palign != NULL)
		*palign = sort_align_by_length(*palign);
	return (*palign);
}

static SeqAlignPtr extract_align_with_sameID(SeqAlignPtr PNTR align, SeqIdPtr sip_1, SeqIdPtr sip_2)
{
	DenseSegPtr dsp;
	SeqAlignPtr prev, next, curr;
	Boolean match;
	SeqAlignPtr h_align, p_align;

	prev = NULL;
	curr = *align;
	p_align = NULL;
	h_align = NULL;
	while(curr)
	{
		match = FALSE;
		next = curr->next;
		if(curr->segtype == 2)
		{
			dsp = (DenseSegPtr) curr->segs;
			if(SeqIdMatch(dsp->ids, sip_1) 
				&& SeqIdMatch(dsp->ids->next, sip_2))
				match = TRUE;
		}

		if(match)
		{
			if(prev == NULL)
				*align = next;
			else
				prev->next = next;
			curr->next = NULL;

			if(p_align == NULL)
				h_align = curr;
			else
				p_align->next = curr;
			p_align = curr;
		}
		else
			prev = curr;
		curr = next;
	}

	return h_align;
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
			align->next = extract_align_with_sameID(&next, 
				dsp->ids, dsp->ids->next);
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

static Int4 get_align_num(SeqAlignPtr align)
{
	Int4 num = 0;

	while(align)
	{
		++num;
		align = align->next;
	}

	return num;
}

static void save_output (SeqIdPtr sip, CharPtr f_name, SeqAlignPtr align)
{
	AsnIoPtr aip;
	BioseqPtr bsp;
	SeqEntryPtr sep;

	bsp = BioseqFind(sip);
	sep = SeqEntryFind(sip);

	if(bsp == NULL || sep == NULL)
		return;

	bsp->hist = SeqHistNew();
	bsp->hist->assembly = align;

	aip = AsnIoOpen(f_name, "w");
	SeqEntryAsnWrite(sep, aip, NULL);
	AsnIoClose(aip);
	
	bsp->hist->assembly = NULL;
	SeqHistFree(bsp->hist);
	bsp->hist = NULL;
}


static SeqAlignPtr find_best_align(SeqAlignPtr align, Int4Ptr pmax_score)
{
        SeqAlignPtr curr, best_align;
        Int4 score, max_score = 0, number;
        Nlm_FloatHi bit_score, evalue;

        best_align = NULL;
        for(curr = align; curr != NULL; curr = curr->next)
        {
                GetScoreAndEvalue(curr, &score, &bit_score, &evalue , &number);
                if(score > max_score)
                {
                        max_score = score;
                        best_align = curr;
                }
        }


	*pmax_score = max_score;
        return best_align;
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


/* Filter SeqAlign by removing bad hits and hits off the main diagonal (as defined by the
   first hit .. which is the highest scoring in blast */
static SeqAlignPtr SeqAlignConsistentDiagFilter(SeqAlignPtr align,SeqLocPtr slp1, SeqLocPtr slp2, 
					 Boolean store_err) {
    Int4 old_num,curr_num,len,len1,len2;

    SeqAlignPtr prev,t_align;

    if(slp1 == NULL || slp2 == NULL)
	return NULL;
    len1 = SeqLocLen(slp1);
    len2 = SeqLocLen(slp2);
    if(len1 == 0 || len2 == 0)
	return NULL;

    if(align != NULL) {
	if(store_err)
	    save_output (SeqLocId(slp1), "blast.out", align);
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
	if(store_err)
	    save_output (SeqLocId(slp1), "merge.out", align);
	if(is_bad_blast_alignment(align, slp1, slp2, 50))
	    {
		SeqAlignSetFree(align);
		return NULL;
	    }
	
	if(align != NULL)
	    {
		old_num = get_align_num(align);
		MergeToOneAlignment(&align); 
		curr_num = get_align_num(align);
		
		if(align == NULL)
		    {
			fprintf(stderr, "Fail in MergeToOneAlignment\n");
		    } 
		while(curr_num > 1 && curr_num < old_num)
		    {
			old_num = curr_num;
			MergeToOneAlignment(&align);
			curr_num = get_align_num(align);
		    }
		
		if(curr_num > 1)
		    {
			prev = align;
			t_align = align->next;
			while(t_align)
			    {
				len = get_align_length(t_align);
				if(len < 50 || get_align_length(prev)/len > 5)
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
	if(store_err)
	    save_output (SeqLocId(slp1), "sim4.out", align);
	if(is_bad_blast_alignment(align, slp1, slp2, MIN(50, MIN(len1, len2)/5)))
	    {
		SeqAlignSetFree(align);
		return NULL;
	    }
	
    }
    return align;
}

static SeqAlignPtr SplitBlastTwoSeq(SeqLocPtr slp1, SeqLocPtr slp2, 
		Int4 split_len, Int4 overlap_len, BLAST_OptionsBlkPtr options, Boolean store_err)
{
	Int4 len1, len2;
	SeqAlignPtr align = NULL;

	if(slp1 == NULL || slp2 == NULL)
		return NULL;
	len1 = SeqLocLen(slp1);
	len2 = SeqLocLen(slp2);
	if(len1 == 0 || len2 == 0)
		return NULL;
	/* Align the sequence "a la Powerblast" */
	align = SeqAlignSplitBlastTwoSeq(slp1, slp2,split_len, overlap_len, options);
	/* Filter out anything that is off the main diagonal */
	if(align)
	    align = SeqAlignConsistentDiagFilter(align,slp1,slp2,store_err);
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
	get_align_ends(align, dsp->ids->next, &s_start, &s_stop, &strand);
	if(s_start<0 || s_stop>= s_bsp->length) {
	    ErrPostEx(SEV_WARNING,0,0,"find_mismatch_residue: corrupted alignment \n");
	    return TRUE;
	} else {
	    s_spp = SeqPortNew(s_bsp, s_start, s_stop, dsp->strands[1], code);
	}
	/* if(s_stop - s_start + 1 < MIN(min_align_len, loc_len-MAX(10, loc_len/10)))
		return TRUE; */

	get_align_ends(align, dsp->ids, &m_start, &m_stop, &strand);
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

/*flip the order of the first and the second sequences*/
static void re_order_alignment(SeqAlignPtr align)
{
	SeqAlignPtr curr;
	DenseSegPtr dsp;
	SeqIdPtr sip;
	Int2 i;
	Int4 temp;
	Uint1 strand;

	for(curr = align; curr != NULL; curr = curr->next)
	{
		reverse_alignment(curr, 1);
		dsp = (DenseSegPtr) curr->segs;

		/*switching the Seq-ids*/
		sip = dsp->ids->next;
		dsp->ids->next = NULL;
		sip->next = dsp->ids;
		dsp->ids = sip;

		for(i = 0; i<dsp->numseg; ++i)
		{
			temp = dsp->starts[2*i];
			dsp->starts[2*i] = dsp->starts[2*i+1];
			dsp->starts[2*i+1] = temp;

			strand = dsp->strands[2*i];
			dsp->strands[2*i] = dsp->strands[2*i+1];
			dsp->strands[2*i+1] = strand;
		}

	}
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

/* Tries to MAke a global alignment from a Set of local SeqAligns. 
 */
static SeqAlignPtr SeqAlignSetGlobalFromLocal(SeqAlignPtr align,SeqLocPtr loc_1, SeqLocPtr loc_2, BioseqPtr bsp_1, BioseqPtr bsp_2, FILE *fp, Boolean record_err)
{
	Int4 len_1, len_2;
	Char label[101];

	len_1 = SeqLocLen(loc_1);
	if(len_1 <=0)
	{
		MuskSeqIdWrite(SeqLocId(loc_1), label, 100, PRINTID_TEXTID_ACCESSION, TRUE, TRUE);
		fprintf(fp, "bad length %ld in sequence %s\n", len_1, label);
		return NULL;
	}
	len_2 = SeqLocLen(loc_2);
	if(len_2 <=0)
	{
		MuskSeqIdWrite(SeqLocId(loc_2), label, 100, PRINTID_TEXTID_ACCESSION, TRUE, TRUE);
		fprintf(fp, "bad length %ld in sequence %s\n", len_1, label);
		return NULL;
	}

	if(align != NULL)
	{

	    align = SeqAlignConsistentDiagFilter(align,loc_1, loc_2,record_err);
	}

	if(align != NULL)
	{
		extend_dsp_ends(bsp_1->length, bsp_2->length, align, MAX_EXTEND_OVERHANG);
		re_order_alignment(align);
		extend_dsp_ends(bsp_2->length, bsp_1->length, align, MAX_EXTEND_OVERHANG);
		re_order_alignment(align);
	}
	return align;
}

/* Attempt to Make a Global Alignment Using  local Alignments tools */
static SeqAlignPtr compute_alignment(SeqLocPtr loc_1, SeqLocPtr loc_2, BioseqPtr bsp_1, BioseqPtr bsp_2, FILE *fp, BLAST_OptionsBlkPtr options, Boolean ck_polyA, Boolean record_err)
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
		MuskSeqIdWrite(SeqLocId(loc_1), label, 100, PRINTID_TEXTID_ACCESSION, TRUE, TRUE);
		fprintf(fp, "bad length %ld in sequence %s\n", len_1, label);
		return NULL;
	}

	if(len_2 <=0)
	{
		MuskSeqIdWrite(SeqLocId(loc_2), label, 100, PRINTID_TEXTID_ACCESSION, TRUE, TRUE);
		fprintf(fp, "bad length %ld in sequence %s\n", len_1, label);
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
			OVERLAP_LEN, options, record_err);
		if(align != NULL)
			re_order_alignment(align);
	}

	if(align != NULL)
	{
		extend_dsp_ends(bsp_1->length, bsp_2->length, align, MAX_EXTEND_OVERHANG);
		re_order_alignment(align);
		extend_dsp_ends(bsp_2->length, bsp_1->length, align, MAX_EXTEND_OVERHANG);
		re_order_alignment(align);
	}
	return align;
}


static void reverse_seqloc (SeqLocPtr s_loc)
{
	SeqLocPtr slp, t_slp, n_slp;
	SeqIntPtr sint;
	Uint1 strand;

	if(s_loc == NULL || s_loc->choice != SEQLOC_MIX)
		return;
	slp = NULL;
	t_slp = NULL;
	while((slp = SeqLocFindNext(s_loc, slp)) != NULL)
	{
		if(slp->choice == SEQLOC_INT)
		{
			sint = (SeqIntPtr) slp->data.ptrvalue;
			strand = sint->strand;
			n_slp = SeqLocIntNew(sint->from, sint->to, 3-strand, sint->id);
			if(t_slp == NULL)
				t_slp = n_slp;
			else
			{
				n_slp->next = t_slp;
				t_slp = n_slp;
			}
		}
	}
	SeqLocSetFree((SeqLocPtr) s_loc->data.ptrvalue);
	s_loc->data.ptrvalue = t_slp;
}
	

static BioseqPtr create_fake_bioseq(SeqLocPtr seq_loc, CharPtr fake_name)
{

	BioseqPtr bsp;
	SeqPortPtr spp;
	ByteStorePtr b_store;
	Uint1 residue;
	Uint1 code;
	SeqLocPtr slp;
	Int4 length;

	code = Seq_code_iupacna;
	length = SeqLocLen(seq_loc);
	b_store = BSNew(length + 2);
	BSSeek(b_store, 0, SEEK_SET);

	slp = NULL;
	length = 0;
	while((slp = SeqLocFindNext(seq_loc, slp)) != NULL)
	{
		if(slp->choice == SEQLOC_INT || slp->choice == SEQLOC_WHOLE)
		{
			spp = SeqPortNewByLoc(slp, code);
			residue = SeqPortGetResidue(spp);
			while(residue != SEQPORT_EOF)
			{
				if(IS_ALPHA(residue))
					BSPutByte(b_store, (Int2)residue);
				residue = SeqPortGetResidue(spp);
			}
			SeqPortFree(spp);
			length += SeqLocLen(slp);
		}
	}

	bsp = BioseqNew();
	bsp->id = local_id_make(fake_name);
	bsp->seq_data = b_store;
	bsp->repr = Seq_repr_raw;
	bsp->mol = Seq_mol_dna;
	bsp->length = length;
	bsp->topology = 1;
	bsp->seq_data_type = code;
	return bsp;


}

static SeqLocPtr dup_mul_locs(SeqLocPtr slp)
{
	Int4 start, stop;
	Uint1 strand;
	SeqLocPtr t_slp, h_slp;
	SeqLocPtr curr;

	curr = NULL;
	h_slp = NULL;
	while((curr = SeqLocFindNext(slp, curr)) != NULL)
	{
		start = SeqLocStart(curr);
		stop = SeqLocStop(curr);
		strand = SeqLocStrand(curr);
		t_slp = SeqLocIntNew(start, stop, strand, SeqLocId(curr));
		ValNodeLink(&h_slp, t_slp);
	}

	if(h_slp->next == NULL)
		return h_slp;
	else
	{
		t_slp = ValNodeNew(NULL);
		t_slp->choice = SEQLOC_MIX;
		t_slp->data.ptrvalue = h_slp;
		return t_slp;
	}
}


static Boolean second_seq_has_gap(BioseqPtr s_bsp)
{
	SeqLocPtr slp;
	SeqIdPtr sip;
	BioseqPtr t_bsp;

	for(slp = (SeqLocPtr) s_bsp->seq_ext; slp != NULL; slp = slp->next)
	{
		if(slp->choice != SEQLOC_NULL)
		{
			sip = SeqLocId(slp);
			if(sip != NULL)
			{
				t_bsp = BioseqFind(sip);
				if(t_bsp != NULL)
				{
					if(t_bsp->repr == Seq_repr_virtual)
						return TRUE;
				}
			}
		}
	}

	return FALSE;
}


static SeqIdPtr get_best_id(BioseqPtr bsp)
{
	SeqIdPtr sip;

	sip = SeqIdFindBest(bsp->id, 0);
	if(sip == NULL || sip->choice != SEQID_GI)
	{
		for(sip = bsp->id; sip != NULL; sip = sip->next)
		{
			if(sip->choice == SEQID_LOCAL)
				return sip;
		}
		return bsp->id;
	}
	else
		return sip;
}
	
static SeqAlignPtr process_multi_cds(BioseqPtr fake_gcontig_bsp, SeqLocPtr gcontig_loc, BioseqPtr s_bsp, BioseqPtr contig_bsp, Boolean is_unigene, FILE *err_fp, BLAST_OptionsBlkPtr options)
{
	SeqLocPtr s_loc, m_loc;
	Int4 start, stop;
	Int4 seg_start, seg_stop;
	SeqLocPtr curr,align_loc;
	Uint1 strand;
	SeqAlignPtr align, c_align;
	DenseSegPtr dsp;
	Int2 i;
	Int4 gap_len;
	Int4 p_right, p_left, left, right;
	ValNodePtr v_starts = NULL, v_lens = NULL;
	ValNodePtr v_strands = NULL;
	Int2 curseg;
	Int4 offset;
	Int4 length;
        Int4 m_start;
	BioseqPtr fake_second_bsp = NULL;
	ValNode vn;
	Char label[101], label_1[101];
	Boolean record_err = FALSE;


	MuskSeqIdWrite(s_bsp->id, label_1, 100, 
		PRINTID_FASTA_LONG, TRUE, TRUE);
	MuskSeqIdWrite(contig_bsp->id, label, 100, 
		PRINTID_FASTA_LONG, TRUE, TRUE);
	/* if(StrStr(label, "992"))
		record_err = TRUE; */
	if(s_bsp->repr == Seq_repr_seg)
	{
		if(second_seq_has_gap(s_bsp))
		{
			fprintf(err_fp, "%s: Has Gap\n", label);
			return NULL;
		}
		vn.choice = SEQLOC_MIX;
		vn.data.ptrvalue = s_bsp->seq_ext;
		vn.next = NULL;
		fake_second_bsp = create_fake_bioseq((SeqLocPtr)&vn, "temp2");
		if(fake_second_bsp->length == 0)
		{
			printf("%s\n", label);
			BioseqFree(fake_second_bsp);
			return NULL;
		}
		s_loc = SeqLocIntNew(0, fake_second_bsp->length-1, 
			Seq_strand_both, fake_second_bsp->id);
	}
	else
		s_loc = SeqLocIntNew(0, s_bsp->length-1, 
			Seq_strand_both, SeqIdFindBest(s_bsp->id, 0));
	m_loc = SeqLocIntNew(0, fake_gcontig_bsp->length-1, Seq_strand_plus, 
		fake_gcontig_bsp->id);

	if(fake_second_bsp != NULL)
		align = compute_alignment(m_loc, s_loc, fake_gcontig_bsp, 
			fake_second_bsp, err_fp, options, is_unigene, record_err);
	else
		align = compute_alignment(m_loc, s_loc, fake_gcontig_bsp, 
			s_bsp, err_fp, options, is_unigene, record_err);
	SeqLocFree(s_loc);
	SeqLocFree(m_loc);


	if(align != NULL)
	{
		for(c_align = align; c_align != NULL; c_align = c_align->next)
		{
			dsp = (DenseSegPtr) c_align->segs;
			SeqIdSetFree(dsp->ids);
			dsp->ids = SeqIdDup(get_best_id(contig_bsp));
			dsp->ids->next = SeqIdDup(get_best_id(s_bsp));
		}
	}
	if(fake_second_bsp != NULL)
		BioseqFree(fake_second_bsp);

	if(align == NULL)
	{
		return NULL;
	}

	if(gcontig_loc->choice != SEQLOC_MIX)
		return align;

	/*the original code has put the MIX loc as the second Bioseq*/
	re_order_alignment(align);
	dsp = (DenseSegPtr) align->segs;
	s_loc = dup_mul_locs(gcontig_loc);
	if(dsp->strands[1] == Seq_strand_minus)
		reverse_seqloc(s_loc);
	strand = dsp->strands[1];
	align_loc = (SeqLocPtr) s_loc->data.ptrvalue;



	curseg = 0;
	for(i = 0; i<dsp->numseg; ++i)
	{
		start = dsp->starts[2*i +1];
		m_start = dsp->starts[2*i];
		if(start != -1)
		{
			stop = start + dsp->lens[i] -1;
			p_right = -1;
			p_left = -1;
			gap_len = 0;

			/*keep track of the coordinates on the fake-bioseq*/
			if(strand == Seq_strand_minus)
			{
				seg_stop = SeqLocLen(s_loc) -1;
				seg_start = seg_stop + 1;
			}
			else
			{
				seg_start = 0;
				seg_stop = -1;
			}
			for(curr = align_loc; curr != NULL; curr = curr->next)
			{
				left = GetOffsetInBioseq(curr, contig_bsp, SEQLOC_LEFT_END);
				right = GetOffsetInBioseq(curr, contig_bsp, SEQLOC_RIGHT_END);
				swap (&left, &right);

				if(strand == Seq_strand_minus)
				{
					if(p_left != -1)
						gap_len = p_left - (right +1);
					else
						gap_len = 0;
				}
				else
				{
					if(p_right != -1)
						gap_len = left - (p_right+1);
					else
						gap_len = 0;
				}
				p_left = left;
				p_right = right;
					
				if(strand == Seq_strand_minus)
					seg_start -= SeqLocLen(curr);
				else
					seg_stop += SeqLocLen(curr);
				/*current alignment overlaps with THIS segment*/
				if(!(start > seg_stop || stop < seg_start))
				{
					if(strand == Seq_strand_minus)
						offset = MAX(0, seg_stop - stop);
					else
						offset = MAX(0, start - seg_start);

					/*check for the gaps. ignore the end gaps*/
					if(curseg > 0 && gap_len > 0)
					{
						if(offset == 0) /*add the gaps*/
						{
							ValNodeAddInt(&v_starts, 0, -1);
							if(strand == Seq_strand_minus)
								ValNodeAddInt(&v_starts, 0, right+1);
							else
								ValNodeAddInt(&v_starts, 0, (left - gap_len));
							ValNodeAddInt(&v_lens, 0, gap_len);
							ValNodeAddInt(&v_strands, 0, 1);
							ValNodeAddInt(&v_strands, 0, (Int4)strand);
							++curseg;
						}
					}
					length  = MIN(stop, seg_stop) - MAX(start, seg_start) +1;
					ValNodeAddInt(&v_starts, 0, m_start);
					if(strand == Seq_strand_minus)
						ValNodeAddInt(&v_starts, 0, (right - offset - (length - 1)));
					else
						ValNodeAddInt(&v_starts, 0, left + offset);
					ValNodeAddInt(&v_lens, 0, length);
					if(m_start != -1)
						m_start += length;
					if(strand == Seq_strand_minus)
						stop -= length;
					else
						start += length;
					ValNodeAddInt(&v_strands, 0, 1);
					ValNodeAddInt(&v_strands, 0, (Int4)strand);
					++curseg;
				}
				if(strand == Seq_strand_minus)
					seg_stop = seg_start -1;
				else
					seg_start = seg_stop +1;
			}
		}
		else	/*gap on the mixed Seq-loc*/
		{
			ValNodeAddInt(&v_starts, 0, m_start);
			ValNodeAddInt(&v_starts, 0, -1);
			ValNodeAddInt(&v_lens, 0, dsp->lens[i]);
			ValNodeAddInt(&v_strands, 0, 1);
			ValNodeAddInt(&v_strands, 0, (Int4)strand);
			++curseg;
		}
	}
	dsp = make_dsp(v_starts, v_lens, v_strands, 1, 2, dsp, curseg, TRUE, NULL, NULL);
	align->segs = dsp;
	re_order_alignment(align);
	SeqLocFree(s_loc);
	return align;
}
								

static SeqIdPtr find_gcontig_id(BioseqPtr bsp)
{
	ObjectIdPtr oip;

	SeqIdPtr sip;

	
	sip = SeqIdFindBest(bsp->id, SEQID_GI);
	if(sip && sip->choice == SEQID_GI)
		return sip;

	for(sip = bsp->id; sip != NULL; sip = sip->next)
	{
		if(sip->choice == SEQID_GENERAL)
			return sip;
	}

	for(sip = bsp->id; sip != NULL; sip = sip->next)
	{
		if(sip->choice == SEQID_LOCAL)
		{
			oip = (ObjectIdPtr) sip->data.ptrvalue;
			if(oip->str && StringNCmp(oip->str, "HGM", 3) == 0)
				return sip;
		}

	}

	return (bsp->id);
}

/*gap = unknown length is separated as separate segments. This 
  is to make sure that create_fake_bioseq won't include those 
  segments
 */

static Boolean is_virtual_segment(SeqIdPtr sip)
{
	BioseqPtr bsp;
	ObjectIdPtr oip;

	bsp = BioseqFind(sip);
	if(bsp != NULL)
	{
		if(bsp->repr == Seq_repr_virtual || bsp->repr == Seq_repr_map)
			return TRUE;
		else
			return FALSE;
	}
	if(sip->choice == SEQID_LOCAL)
	{
		oip = (ObjectIdPtr) sip->data.ptrvalue;
		return (oip->str && StringCmp(oip->str, "virtual") == 0);
	}

	return FALSE;
}
	

static SeqLocPtr build_slp_for_gcontig (BioseqPtr gcontig_bsp, SeqIdPtr gcontig_sip)
{
	SeqLocPtr slp, h_slp = NULL, t_slp;
	Int4 start = 0, stop = -1;
	SeqIdPtr sip;
	Boolean update;


	for(slp = (SeqLocPtr) gcontig_bsp->seq_ext; slp != NULL; slp = slp->next)
	{
		if(slp->choice != SEQLOC_NULL)
		{
			sip = SeqLocId(slp);
			/* if(sip->choice != SEQID_GI) */
			if(sip != NULL)
			{
				update = FALSE;
				if(is_virtual_segment(sip))
				{
					t_slp = SeqLocIntNew(start, stop, Seq_strand_plus, gcontig_sip);
					ValNodeLink(&h_slp, t_slp);
					update = TRUE;
				}
				stop += SeqLocLen(slp);
				if(update)
					start = stop + 1;
				
			}
		}
		else
		{
			if(start < stop)
			{
				t_slp = SeqLocIntNew(start, stop, Seq_strand_plus, gcontig_sip);
				ValNodeLink(&h_slp, t_slp);
				start = stop + 1;
			}
		}
	}

	if(start < stop)
	{
		t_slp = SeqLocIntNew(start, stop, Seq_strand_plus, gcontig_sip);
		ValNodeLink(&h_slp, t_slp);
	}

	if(h_slp == NULL)
		return NULL;
	
	if(h_slp->next != NULL)
	{
		t_slp = ValNodeNew(NULL);
		t_slp->choice = SEQLOC_MIX;
		t_slp->data.ptrvalue = h_slp;
		return t_slp;
	}
	else
		return h_slp;
}
	

