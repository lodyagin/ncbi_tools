/* $Id: txalign.c,v 6.148 2000/10/27 17:54:17 madden Exp $
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
* ===========================================================================
*
* File Name:  txalign.c
*
* $Revision: 6.148 $
* 
* File Description:  Formating of text alignment for the BLAST output
*
* Modifications:
* --------------------------------------------------------------------------
* $Log: txalign.c,v $
* Revision 6.148  2000/10/27 17:54:17  madden
* Changes to make_dumpgnl_links for new hs genome page
*
* Revision 6.147  2000/10/12 21:37:32  shavirin
* Adjusted calculation of ends of alignment in minus strand.
*
* Revision 6.146  2000/10/06 19:30:51  shavirin
* Added printing of initial frame number in OOF case. Fixed some spacing
* to be the same as in regular case.
*
* Revision 6.145  2000/10/06 17:55:44  shavirin
* Added usage of correct matrix in OOF case.
*
* Revision 6.144  2000/10/06 17:23:21  shavirin
* Added BioseqUnlock in printing OOF alignment.
*
* Revision 6.142  2000/10/02 22:03:26  shavirin
* Changed function OOFShowSingleAlignment to have correct spacing.
*
* Revision 6.141  2000/09/28 15:51:44  dondosha
* Open <PRE> block in PrintDefLinesFromSeqAlignEx2 - needed for PSI BLAST
*
* Revision 6.140  2000/09/28 15:03:15  dondosha
* Added boolean splice_junction score type
*
* Revision 6.139  2000/09/27 20:57:55  shavirin
* Fixed bug with printing DNA line ends in case of minus strand.
*
* Revision 6.138  2000/09/25 19:22:12  shavirin
* Fixed start of protein sequencer in OOF alignment. Added check for NULL
* return from BioseqLocById in FormatScoreFromSeqAlignEx() function.
*
* Revision 6.137  2000/09/13 22:24:30  dondosha
* Corrected the printing of </PRE> at the end of PrintDefLinesFromSeqAlignEx2
*
* Revision 6.136  2000/09/13 21:15:39  dondosha
* Removed opening <PRE> in PrintDefLinesFromSeqAlignEx2
*
* Revision 6.135  2000/09/01 18:43:38  shavirin
* Adjusted start and stop of every line for OOF alignment printout.
*
* Revision 6.134  2000/08/31 16:53:10  shavirin
* Fixed memory leak in OOFShowSingleAlignment().
*
* Revision 6.133  2000/08/30 14:18:42  shavirin
* Fixed case for printing OOF alignment.
*
* Revision 6.132  2000/08/25 19:02:37  shavirin
* Corrected calculation of number of mismatches, gaps,positives etc.
* for discontinuous alignments.
*
* Revision 6.131  2000/08/24 18:14:54  shavirin
* Added "<a name=..." links to every HSP Score for Greg's BLAST page.
* This easyly may be changed for general Blast output case.
*
* Revision 6.130  2000/07/25 16:47:13  shavirin
* Changed function to print OOF alignment.
*
* Revision 6.129  2000/07/18 22:37:22  shavirin
* Adjusted end_of_line values in the function OOFShowSingleAlignment()
*
* Revision 6.128  2000/07/17 14:11:47  shavirin
* Adjusted function OOFShowSingleAlignment()
*
* Revision 6.127  2000/07/14 16:02:41  shavirin
* Initialixed variable ooframe to FALSE.
*
* Revision 6.126  2000/07/11 20:51:04  shavirin
* Added major functions for displaying Out-Of-Frame alignments.
*
* Revision 6.125  2000/07/10 20:45:53  shavirin
* Added parameter ooframe for Out-Of-frame alignment and corresponding changes
* to accomodate this parameter.
*
* Revision 6.124  2000/06/22 18:56:55  egorov
* Add a protection against empty deflines.
*
* Revision 6.123  2000/06/19 12:53:18  madden
* Do SeqIdWrite for both HTML and text
*
* Revision 6.122  2000/06/16 18:25:43  shavirin
* Fixed problem with removing full path of db in make_dumpgnl_links()
*
* Revision 6.121  2000/06/16 16:18:46  madden
* Roll back change from rev. 6.114
*
* Revision 6.120  2000/06/15 15:35:47  shavirin
* Fixed Uninitialized memory read error in make_dumpgnl_links() function
*
* Revision 6.119  2000/06/13 18:58:51  shavirin
* Adjusted region of database sequence in the function
* load_align_sum_for_DenseDiag()
*
* Revision 6.118  2000/06/12 16:50:03  shavirin
* Fixed bug with calculation total length of the StdSeg alignment and
* adjusted calculation of translation frame.
*
* Revision 6.117  2000/06/09 19:00:05  shavirin
* Function GetGeneticCodeFromSeqId() made external and added to header file.
*
* Revision 6.116  2000/06/08 20:44:49  shavirin
* Added calculation of start/stop values in the function find_score_in_align().
*
* Revision 6.115  2000/06/08 17:13:53  dondosha
* Fixed bug with wrong scores reported for ungapped blastx alignments
*
* Revision 6.114  2000/06/06 16:40:20  shavirin
* Use plain database name if TXALIGN_NO_ENTREZ option is set.
*
* Revision 6.113  2000/05/16 16:32:17  shavirin
* Added check for WWW_ROOT_PATH environment for PSI Blast.
*
* Revision 6.112  2000/05/05 20:23:08  shavirin
* Do not make gnl-link if passwd or tool_url == NULL.
*
* Revision 6.111  2000/05/05 20:03:48  shavirin
* Rolled back revision 6.110.
*
* Revision 6.109  2000/05/01 19:09:58  shavirin
* Removed function SeqIdSetDup()
*
* Revision 6.108  2000/05/01 16:26:14  shavirin
* Added multiple-highligted deflines in the BLAST output.
*
* Revision 6.107  2000/04/25 18:13:43  shavirin
* Do not link to anything ids with BL_ORD_ID.
*
* Revision 6.106  2000/04/04 21:52:50  madden
* Roll-back last change
*
* Revision 6.105  2000/04/03 17:45:42  shavirin
* Changed way to print multiple deflines. Removed define to print
* old Entrez links.
*
* Revision 6.104  2000/03/24 16:04:24  shavirin
* Added hack for Drosophila BLAST page.
*
* Revision 6.103  2000/03/23 15:02:12  shavirin
* Added possibility to use environment variables in the function
* make_dumpgnl_links()
*
* Revision 6.102  2000/03/14 17:16:11  shavirin
* Cleared AlignSum buffer in the function FormatScoreFromSeqAlign
*
* Revision 6.101  2000/03/07 21:58:40  shavirin
* Now will use PSSM Matrix to show positives in PSI Blast
*
* Revision 6.100  2000/03/02 16:25:09  shavirin
* Fixed bug with very long deflines in FilterTheDefline() function.
*
* Revision 6.99  2000/01/19 21:54:31  madden
* Moved vecscreen stuff to vecscrn.[ch]
*
* Revision 6.98  2000/01/07 21:08:28  shavirin
* Fixed minor memory leak in ShowTextAlignFromAnnot2().
*
* Revision 6.97  1999/12/02 19:36:10  shavirin
* Fixed hundreds of errors detected by C++ compiler. Fixed formating
* of deflines in PHI/PSI Blast search.
*
* Revision 6.96  1999/11/24 21:24:31  vakatov
* Fixed for the C++ and/or MSVC DLL compilation
*
* Revision 6.95  1999/11/16 20:44:07  egorov
* Close all <PRE>'s
*
* Revision 6.94  1999/11/12 19:01:40  madden
* Fix memory leaks
*
* Revision 6.93  1999/11/10 18:27:11  madden
* Fix problem with determining if same sequence as last is being examined
*
* Revision 6.92  1999/11/09 22:15:07  shavirin
* Added parameter follower to the Blast score printing function
*
* Revision 6.91  1999/11/01 16:50:25  shavirin
* Fixed typo in the function get_seqid_for_textbuf()
*
* Revision 6.90  1999/11/01 15:38:13  shavirin
* Turned on producing of new Entrez links in the BLAST output.
*
* Revision 6.89  1999/10/19 18:25:43  shavirin
* Added prefix "images" to psi_blast gifs.
*
* Revision 6.88  1999/10/07 16:08:04  shavirin
* Passed matrix to the function FormatScoreFromSeqAlign().
*
* Revision 6.87  1999/10/07 13:13:44  egorov
* Fix bug when preprocessor direrective and C statement were at the same line.
*
* Revision 6.86  1999/09/30 20:43:34  madden
* change ray links for VecScreen
*
* Revision 6.85  1999/09/29 17:15:37  shavirin
* Added new funtion FormatScoreFromSeqAlign()
*
* Revision 6.84  1999/09/29 14:07:56  shavirin
* Fixed typo in webb_blossum62 matrix.
*
* Revision 6.83  1999/09/28 20:11:34  shavirin
* Added Id, Revision and Log information.
*
*
* ==========================================================================
*/
#include <txalign.h>
#include <codon.h>
#include <ncbimisc.h>
#include <salpacc.h>
#include <salpstat.h>

#define BUFFER_LENGTH 2048
#define MIN_INS_SPACE 50
#define MAX_GI_NUM    10
#define MAX_DB_NUM    10

#define TXALIGN_HREF "http://www.ncbi.nlm.nih.gov"

#define NEW_ENTREZ_HREF "http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi"

/* Used in make_dumpgnl_links, set in getreq.cpp or getreqcmd.cpp */
const char *RID_glb;

/*
	Used by the functions that format the one-line descriptions.
*/
typedef struct _txdfline_struct {
        struct _txdfline_struct *next;
        SeqAlignPtr seqalign;
        SeqIdPtr id;
        Char *buffer_id;
        Char *title;
        Nlm_FloatHi bit_score;
        Nlm_FloatHi evalue;
        Int4 score;
        Int4 number;
	Boolean is_na;
	Boolean found_score;
	Boolean isnew;		/* used to print mark "new.gif" near defline */
	Boolean waschecked;	/* used to print some another .gif near such defline */
	CharPtr segs_str;	/* Used to print segs for dumpgnl program. */
	size_t  segs_buflen,
		segs_used;
} TxDfLineStruct, *TxDfLineStructPtr;


static Int4 get_num_empty_space(Boolean compress)
{
	return (compress ? (8+5 +1) : B_SPACE+POS_SPACE+STRAND_SPACE +1);
}

static Boolean SeqAlignSegsStr(SeqAlignPtr salp, Int2 index, CharPtr *dst, size_t *size, size_t *used);
static CharPtr StringAppend(CharPtr *dst, size_t *size, CharPtr src, size_t *used);


static ValNodePtr ProcessTextInsertion PROTO((AlignNodePtr anp, Int4 m_left, Int4 m_right, BioseqPtr bsp, Int4 line_len, Int1 frame));

static ValNodePtr ExtractCurrentAlignNode(ValNodePtr PNTR anp_list)
{
	ValNodePtr head, curr, prev = NULL;
	AlignNodePtr anp;
	Uint2 itemID, chain;

	head = *anp_list;
	while(head && head->choice == OBJ_SEQANNOT)
		head = head->next;
	if(head == NULL)
		return NULL;

	anp = (AlignNodePtr) head->data.ptrvalue;
	itemID = anp->itemID;
	chain = anp->chain;

	head = *anp_list;
	curr = *anp_list;
	while(curr)
	{
		if(curr->choice != OBJ_SEQANNOT)
		{
			anp = (AlignNodePtr) curr->data.ptrvalue;
			if(anp->itemID != itemID || anp->chain != chain)
			{
				*anp_list = curr;
				if(prev !=NULL)
					prev->next = NULL;
				return head;
			}
		}
		else
		{
			*anp_list = curr;
			if(prev != NULL)
				prev->next = NULL;
			return head;
		}
		prev = curr;
		curr = curr->next;
	}
	*anp_list = NULL;
	return head;
}
	

static void modify_kludge_itemID (ValNodePtr anp_list, Uint2 itemID)
{
	AlignNodePtr anp;

	while(anp_list)
	{
		if(anp_list->choice != OBJ_SEQANNOT)
		{
			anp = (AlignNodePtr) anp_list->data.ptrvalue;
			anp->itemID = itemID;
		}
		anp_list = anp_list->next;
	}
}

/******************************************************************
*
*	LoadFollowerForSameId(anp_list)
*	if the same sequence appears multiple times in the anp_list, 
*	it will be moved to the sequence that are the head of this 
*	list. The field anp->follower is set as the order of the 
*	repeats in this list
*
******************************************************************/
static void LoadFollowerForSameId(ValNodePtr anp_list)
{
	ValNodePtr curr, n_curr;
	AlignNodePtr anp, n_anp;

	curr = anp_list;
	while(curr)
	{
		if(curr->choice != OBJ_SEQANNOT)
		{
			anp = (AlignNodePtr) curr->data.ptrvalue;
			if(anp->is_master == FALSE && anp->follower == FALSE)
			{
				for(n_curr = curr->next; n_curr != NULL; n_curr = n_curr->next)
				{
					if(n_curr->choice != OBJ_SEQANNOT)
					{
						n_anp = (AlignNodePtr) n_curr->data.ptrvalue;
						if(n_anp->is_master == FALSE && n_anp->follower == FALSE)
						{
							if(SeqIdMatch(n_anp->sip, anp->sip))
								n_anp->follower = TRUE;
						}
					}
				}
			}
		}
		curr = curr->next;
	}
}


static void MaskWithLowComplexity(ByteStorePtr bsp, SeqLocPtr maskloc, Uint1 mol)
{
	SeqLocPtr slp = NULL;
	Int4 start, stop;
	Uint1 res = 'N';
	

	if(mol == Seq_mol_aa)
		res = 'X';

	while(maskloc)
	{
		slp = NULL;
 		while((slp = SeqLocFindNext(maskloc, slp))!=NULL)
		{
			start = SeqLocStart(slp);
			stop = SeqLocStop(slp);
			BSSeek(bsp, start, SEEK_SET);
			for(; start <=stop; ++start)
                BSPutByte(bsp, (Int2)res);
		}
		maskloc = maskloc->next;
	}
}

static ByteStorePtr create_byte_store_from_bsp (BioseqPtr bsp)
{
	SeqPortPtr spp;
	Uint1 code;
	ByteStorePtr b_store;
	Uint1 residue;

	if(bsp == NULL)
		return NULL;
	if(bsp->mol == Seq_mol_aa)
		code = Seq_code_iupacaa;
	else
		code = Seq_code_iupacna;

	spp = SeqPortNew(bsp, 0, bsp->length-1, Seq_strand_plus, code);
	b_store = BSNew(bsp->length +1);
	BSSeek(b_store, 0, SEEK_SET);
	while((residue = SeqPortGetResidue(spp)) != SEQPORT_EOF)
		BSPutByte(b_store, (Int2)residue);
	SeqPortFree(spp);
	return b_store;
}

static ValNodePtr CreateMaskByteStore (ValNodePtr mask_list)
{
	BioseqPtr bsp;
	SeqLocPtr slp;
	SeqIdPtr sip;
	ValNodePtr list, curr;
	ByteStorePtr b_store, c_store;
	Uint1 mol;

	list = NULL;
	b_store = NULL;
	while(mask_list)
	{
		curr = ValNodeNew(list);
		curr->choice = mask_list->choice;
		if(list == NULL)
			list = curr;
		slp = (ValNodePtr) mask_list->data.ptrvalue;
		if(slp != NULL)
		{
			if(b_store == NULL)
			{
				sip = SeqLocId(slp);
				if(sip != NULL)
				{
					bsp = BioseqLockById(sip);
					if(bsp != NULL)
					{
						b_store = create_byte_store_from_bsp (bsp);
						mol = bsp->mol;
						BioseqUnlock(bsp);
					}
				}
			}
			if(b_store != NULL)
			{
				if(mask_list->next == NULL)
				{
					c_store = b_store;
					b_store = NULL;
				}
				else
					c_store = BSDup(b_store);
				MaskWithLowComplexity(c_store, slp, mol);
				curr->data.ptrvalue = c_store;
			}
		}

		mask_list = mask_list->next;
	}

	if(b_store != NULL)
		BSFree(b_store);
	return list;
}

static void FreeByteStoreList (ValNodePtr bs_list)
{
	ByteStorePtr bsp;
	ValNodePtr curr;

	for(curr = bs_list; curr != NULL; curr = curr->next)
	{
		bsp = (ByteStorePtr) curr->data.ptrvalue;
		if(bsp != NULL)
			BSFree(bsp);
	}
	ValNodeFree(bs_list);
}

static Boolean replace_bytestore_data (BioseqPtr bsp, ValNodePtr bs_list, Uint1 frame)
{
	ByteStorePtr b_store;
	Uint1 code;

	if(bsp == NULL)
		return FALSE;

	if(bsp->mol == Seq_mol_aa)
		code = Seq_code_iupacaa;
	else
		code = Seq_code_iupacna;

	while(bs_list)
	{
		if(bs_list->choice == frame)
		{
			b_store = (ByteStorePtr) bs_list->data.ptrvalue;
			if(b_store != NULL)
			{
				bsp->repr = Seq_repr_raw;
				bsp->seq_data = b_store;
				bsp->seq_data_type = code;
				return TRUE;
			}
		}
		bs_list = bs_list->next;
	}

	return FALSE;
}


/*can the current alignnode be printed for text view*/
NLM_EXTERN Boolean PrintAlignForText(AnnotInfoPtr info, AlignNodePtr anp)
{
        if(anp == NULL || anp->segs == NULL)
                return FALSE;
        if(anp->segs->type == STD_SEG)
        {
                if(info == NULL)
                        return FALSE;
                if(info->annot_type != ANNOT_BLAST)
                        return FALSE;
                if(info->blast_type != ALIGN_BLASTX && 
                        info->blast_type != ALIGN_TBLASTN && 
                                info->blast_type != ALIGN_TBLASTX)
                        return FALSE;
        }

        return TRUE;
}


/*
*	for tblastn and blastx, return the frame of the non-master
*	for tblastx, return the frame for the master sequence
*/

static Boolean is_master_alignment (AlignNodePtr anp, BioseqPtr m_bsp)
{
	return (anp->is_master || BioseqMatch(m_bsp, anp->sip));
}

static Int1 get_alignment_frame(ValNodePtr anp_list, BioseqPtr m_bsp)
{
	Uint1 c_type = 0;
	AlignNodePtr anp;
	AnnotInfoPtr annot_info;
	ValNodePtr  curr;
	Boolean found;

	if(anp_list == NULL)
		return -1;

	annot_info = NULL;
	found = FALSE;
	for(curr = anp_list; curr != NULL; curr = curr->next)
	{
		if(curr->choice == OBJ_SEQANNOT)
		{
			annot_info = (AnnotInfoPtr) anp_list->data.ptrvalue;
			c_type = get_alignment_type(annot_info);
		}
		else
		{
			anp = (AlignNodePtr) curr->data.ptrvalue;
			if(!is_master_alignment(anp, m_bsp))
				found = (c_type != ALIGN_TDNA_TO_TDNA);
			else
				found = (c_type == ALIGN_TDNA_TO_TDNA);
			if(found)
			{
				if(!PrintAlignForText(annot_info, anp))
					return -1;
				if(c_type == ALIGN_NORMAL || c_type == ALIGN_PROT_TO_DNA)
					return 0;
				if(c_type == ALIGN_DNA_TO_PROT || c_type == ALIGN_TDNA_TO_TDNA)
					return anp->m_frame;
			}
		}
	}

	return -1;
}

/**********************************************************************************
*
*       Given a chain of annots (ValNodePtrs) they are all printed out, one pattern
*       at a time.
*	
*	For a give annot all alignments from one database sequence are assumed to be grouped together.
*	
*	The Alignments from one databases sequence are currently ranked by expect value.
*	It has been suggested that this be changed and should not be relied on indefinitely.
*
*************************************************************************************/

NLM_EXTERN Boolean LIBCALL
ShowTextAlignFromAnnotExtra(BioseqPtr bsp, ValNodePtr vnp, SeqLocPtr seqloc,
        Int4 line_len, FILE *fp,
        Uint1Ptr featureOrder, Uint1Ptr groupOrder, Uint4 option, Int4Ptr PNTR matrix,
        ValNodePtr mask_loc, int (LIBCALLBACK *fmt_score_func)PROTO((AlignStatOptionPtr)))
{
        Int4 index=0;
        SeqAnnotPtr seqannot;
        SeqAnnotPtr annot;
        SeqFeatPtr sfp;
        SeqLocPtr next;

        seqannot = SeqAnnotNew();
        seqannot->type = 2;
        AddAlignInfoToSeqAnnot(seqannot, 2);

        while (vnp && seqloc)
        {
                index++;
                seqannot->data = vnp->data.ptrvalue;
                if (bsp->annot)
                        bsp->annot = SeqAnnotFree(bsp->annot);
                annot = bsp->annot = SeqAnnotNew();
                bsp->annot->type = 1;   /* ftable. */
                next = seqloc->next;
                sfp = SeqFeatNew();
                seqloc->next = NULL;
                sfp->location = seqloc;
                sfp->data.choice = SEQFEAT_REGION;
                sfp->data.value.ptrvalue = StringSave("pattern");
                annot->data = sfp;
                fprintf(fp, "\nSignificant alignments for pattern occurrence %ld at position %ld\n\n",
                        (long) index, (long) (SeqLocStart(seqloc)+1));
                ShowTextAlignFromAnnot(seqannot, line_len, fp, featureOrder, groupOrder, option, matrix, mask_loc, fmt_score_func);
                seqloc->next = next;
                seqloc = seqloc->next;
                vnp = vnp->next;
        }
        seqannot->data = NULL;
	seqannot = SeqAnnotFree(seqannot);
	
        return TRUE;
}
static Boolean load_master_translate_frame PROTO((ValNodePtr anp_list, Int4 m_len, BioseqPtr m_bsp));
static AlignNodePtr get_master_align_node PROTO((ValNodePtr anp_list));


/*****************************************************************************
*
*	ShowTextAlignFromAnnot(annot, locus, line_len, fp, master, f_order)
*	display the alignment stored in a Seq-annot in a text file
*	annot: the Seq-annot pointer
*	locus: if TRUE, show the locus name as the sequence label, otherwise, 
*		use the accession
*	line_len: the number of sequence char per line
*	fp: The file pointer to store the text output
*	master: if TRUE, show the result as a master-slave type multiple pair
*	wise alignment. if FALSE, display one alignment after the other
*	f_order: the user selected feature type and order to be shown together
*	with the alignment
*	is_html: print out the format as an HTML page?
*	return TRUE for success, FALSE for fail
*
*****************************************************************************/
/* This modification of the function to pass position-specific matrix */
NLM_EXTERN Boolean ShowTextAlignFromAnnot3(SeqAnnotPtr hannot, Int4 line_len, FILE *fp, Uint1Ptr featureOrder, Uint1Ptr groupOrder, Uint4 option, Int4Ptr PNTR matrix, ValNodePtr mask_loc, int (LIBCALLBACK *fmt_score_func)PROTO((AlignStatOptionPtr)), CharPtr db_name, CharPtr www_blast_type, Int4Ptr PNTR posMatrix)
{
    SeqAlignPtr align, h_align, n_align, prev;
    SeqLocPtr m_loc;
    Boolean retval, matrix_loaded=FALSE;
    SeqIdPtr m_sip;
    Int4 m_start, m_stop, t_start, t_stop;
    ValNodePtr anp_node, curr_list;
    ValNodePtr annot_head;
    Uint1 style;
    SeqAnnotPtr annot;
    Boolean master;
    Boolean flat_insert;
    Uint1 m_strand;
    Uint1 annot_type;
    Char annotDB[101];
    Uint2 order;
    Uint1   blast_type;
    Boolean load_matrix;
    ValNodePtr bs_list;	/*store the ByteStores that were masked master seqences*/
    BioseqPtr m_bsp;
    Int1 frame;
    Uint1 code;
    Uint1 repr;
    ByteStorePtr seq_data;
    
    
    annot = hannot;
    flat_insert = (Boolean)(option & TXALIGN_FLAT_INS); 
    /* flat_insert = TRUE; */
    master = (Boolean)(option & TXALIGN_MASTER);
    if(annot->type != 2)
        return FALSE;
    m_start = -1;
    m_stop = -1;
    m_sip = NULL;
    annot = hannot;
    while(annot) {
        if(annot->type == 2) {
            align = (SeqAlignPtr) annot->data;
            if(m_sip == NULL)
                m_sip = make_master(align);
            if(m_sip != NULL) {
                get_boundary(m_sip, &t_start, &t_stop, align);
                if(m_start == -1 || m_start > t_start)
                    m_start = t_start;
                if(m_stop == -1 || m_stop < t_stop)
                    m_stop = t_stop;
            }
        }
        annot = annot->next;
    }
    if(m_sip == NULL || m_start == -1 || m_stop == -1)
        return FALSE;
    
    if(master)
        style = COLLECT_MP;
    else
        style = COLLECT_MD;
    anp_node = NULL;
    m_loc = SeqLocIntNew(m_start, m_stop, Seq_strand_plus, m_sip);
    annot = hannot;
    load_matrix = FALSE;	/*if there is any protein sequence, set the load_matrix to TRUE*/
    while(annot) {
        if(annot->type == 2) {
            annotDB[0] = '\0';
            blast_type = get_align_annot_qual(annot, annotDB, 100, &annot_type);
            if(blast_type == ALIGN_BLASTX || blast_type == ALIGN_TBLASTN 
               || blast_type == ALIGN_TBLASTX)
                load_matrix = TRUE;
            if(blast_type == ALIGN_TBLASTX || 
               (blast_type == ALIGN_BLASTX && annot_type == ANNOT_BLAST
                && (option & TXALIGN_BLASTX_SPECIAL))) { /*!!!!!!!this messes up all the itemIDs and entityIDs !!!!!!*/
                align = (SeqAlignPtr) annot->data;
                prev = NULL;
                h_align = NULL;
                order = 0;
                while(align) {
                    ++order;
                    n_align = align->next;
                    align->next = NULL;
                    if(get_align_ends(align, m_sip, &m_start, &m_stop, &m_strand)) {
                        /* sint = m_loc->data.ptrvalue;
                           sint->strand = m_strand; */
                        update_seq_loc(m_start, m_stop, m_strand, m_loc); 
                        style = COLLECT_MD;
                        master = FALSE;
                        
                        annot->data = align;
                        curr_list  = CollAlignFromSeqAnnot(annot, m_loc, 
                                                           featureOrder, groupOrder, style, FALSE, master, flat_insert);
                        if(curr_list != NULL) {
                            modify_kludge_itemID (curr_list, order);
                            ValNodeLink(&anp_node, curr_list);
                        }
                    }
                    if(prev == NULL)
                        h_align = align;
                    else
                        prev->next = align;
                    prev = align;
                    align = n_align;
                }
                annot->data = h_align;
            } else {
                curr_list  = CollAlignFromSeqAnnot(annot, m_loc, featureOrder, groupOrder, style, FALSE, master, flat_insert);
                if(curr_list != NULL)
                    ValNodeLink(&anp_node, curr_list);
            }
        }
        annot = annot->next;
    }
    SeqLocFree(m_loc);
    if(anp_node == NULL)
        return FALSE;
    
    m_bsp = BioseqLockById(m_sip);
    if(m_bsp == NULL) {
        FreeAlignNode(anp_node);
        return FALSE;
    }
    
    if(mask_loc != NULL)
        bs_list = CreateMaskByteStore (mask_loc);
    else
        bs_list = NULL;
    
    repr = m_bsp->repr;
    seq_data = m_bsp->seq_data;
    code = m_bsp->seq_data_type;
    
    if(matrix == NULL && (option & TXALIGN_MATRIX_VAL || load_matrix)) {
        matrix = load_default_matrix();
        matrix_loaded = TRUE;
    }
    
    if(fmt_score_func != NULL) {
        free_buff();
        init_buff_ex(MAX(80, line_len + 23 + 12));
    }
    if(master) {
        frame = get_alignment_frame(anp_node, m_bsp);
        if(frame != -1 && bs_list != NULL) {
            load_master_translate_frame(anp_node, m_bsp->length, m_bsp);
            
            if(!replace_bytestore_data (m_bsp, bs_list, (Uint1)frame)) {
                m_bsp->repr = repr;
                m_bsp->seq_data = seq_data;
                m_bsp->seq_data_type = code;
            }
        }
        retval = ShowAlignNodeText2(anp_node, -1, line_len, fp, -1, -1, option, matrix, fmt_score_func, db_name, www_blast_type, posMatrix);
        FreeAlignNode(anp_node);
    } else {
        annot_head = NULL;
        if(fmt_score_func != NULL)
            LoadFollowerForSameId(anp_node);
        while(anp_node) {
            if(anp_node->choice == OBJ_SEQANNOT) {
                if(annot_head != NULL) {
                    annot_head->next = NULL;
                    FreeAlignNode(annot_head);
                }
                annot_head = anp_node;
                anp_node = anp_node->next;
            } else {
                curr_list = ExtractCurrentAlignNode(&anp_node);
                if(curr_list) {
                    if(annot_head != NULL) {
                        annot_head->next = curr_list;
                        load_master_translate_frame(annot_head, m_bsp->length, m_bsp);
                        frame = get_alignment_frame(annot_head, m_bsp);
                    } else
                        frame = 0;
                    if(frame != -1 && bs_list != NULL) {
                        if(!replace_bytestore_data (m_bsp, bs_list, (Uint1)frame)) {
                            m_bsp->repr = repr;
                            m_bsp->seq_data = seq_data;
                            m_bsp->seq_data_type = code;
                        }
                    }
                    
                    if(annot_head != NULL) {
                        retval = ShowAlignNodeText2(annot_head, -1, line_len, fp, -1, -1, option, matrix, fmt_score_func, db_name, www_blast_type, posMatrix);
                        annot_head->next = NULL;
                    } else
                        retval = ShowAlignNodeText2(curr_list, -1, line_len, fp, -1, -1, option, matrix, fmt_score_func, db_name, www_blast_type, posMatrix);
                    if(retval == TRUE)
                        fprintf(fp, "\n\n");
                    FreeAlignNode(curr_list);
                }
            }
        }
        if(annot_head != NULL) {
            annot_head->next = NULL;
            FreeAlignNode(annot_head);
        }
    }
    
    m_bsp->repr = repr;
    m_bsp->seq_data = seq_data;
    m_bsp->seq_data_type = code;
    
    if (matrix_loaded)
        free_default_matrix(matrix);
    
    if(fmt_score_func != NULL)
        free_buff();
    if(bs_list != NULL)
        FreeByteStoreList (bs_list);
    BioseqUnlock(m_bsp);
    return retval;
}

NLM_EXTERN Boolean ShowTextAlignFromAnnot(SeqAnnotPtr hannot, Int4 line_len, FILE *fp, Uint1Ptr featureOrder, Uint1Ptr groupOrder, Uint4 option, Int4Ptr PNTR matrix, ValNodePtr mask_loc, int (LIBCALLBACK *fmt_score_func)PROTO((AlignStatOptionPtr)))
{
    return ShowTextAlignFromAnnot2(hannot, line_len, fp, featureOrder, groupOrder, option, matrix, mask_loc, fmt_score_func, NULL, NULL);

}
NLM_EXTERN Boolean ShowTextAlignFromAnnot2(SeqAnnotPtr hannot, Int4 line_len, FILE *fp, Uint1Ptr featureOrder, Uint1Ptr groupOrder, Uint4 option, Int4Ptr PNTR matrix, ValNodePtr mask_loc, int (LIBCALLBACK *fmt_score_func)PROTO((AlignStatOptionPtr)), CharPtr db_name, CharPtr www_blast_type)
{
    return ShowTextAlignFromAnnot3(hannot, line_len, fp, featureOrder,
                                   groupOrder, option, matrix, mask_loc, 
                                   fmt_score_func, db_name, www_blast_type,
                                   NULL); 
}		

/*************************************************************************
*
*	functions and structure related to create a text buffer for the 
*	alignment
*
*************************************************************************/

NLM_EXTERN ValNodePtr FreeTextAlignList(ValNodePtr tdp_list)
{
	TextAlignBufPtr tdp;
	ValNodePtr next;
	Int2 i;

	while(tdp_list)
	{
		next = tdp_list->next;
		tdp_list->next = NULL;
		tdp = (TextAlignBufPtr) tdp_list->data.ptrvalue;
		if(tdp->label)
			MemFree(tdp->label);
		if(tdp->buf)
			MemFree(tdp->buf);
		if(tdp->matrix_val)
			MemFree(tdp->matrix_val);
		if(tdp->exonCount > 0)
		{
			for(i =0; i<3; ++i)
				MemFree(tdp->codon[i]);
		}
		MemFree(tdp);
		MemFree(tdp_list);
		tdp_list = next;
	}

	return NULL;
}

static ValNodePtr load_tdp_data(ValNodePtr PNTR head, CharPtr label, CharPtr text, Uint2 itemID, Uint2 entityID, Uint2 seqEntityID, Uint2 bsp_itemID)
{
	TextAlignBufPtr tdp;

	tdp = (TextAlignBufPtr) MemNew(sizeof(TextAlignBuf));
	tdp->pos = -1;
	tdp->label = label;
	tdp->buf= text;
	tdp->itemID = itemID;
	tdp->entityID = entityID;
	tdp->seqEntityID = seqEntityID;
	tdp->bsp_itemID = bsp_itemID;

	return ValNodeAddPointer(head, 0, (Pointer)tdp);
}

/*********************************************************************
*
*	functions used for producing the Web browser
*
*********************************************************************/
static Boolean find_seqid_for_bioseq(GatherContextPtr gcp)
{
	BioseqPtr bsp;
	ValNodePtr vnp;

	if(gcp->thistype != OBJ_BIOSEQ)
		return FALSE;
	bsp = (BioseqPtr)(gcp->thisitem);
	if(bsp == NULL)
		return FALSE;
	vnp  = (ValNodePtr)(gcp->userdata);
	vnp->choice = bsp->mol;
	vnp->data.ptrvalue = bsp->id;
	return TRUE;
}
	
static Boolean find_seqid_for_seqfeat(GatherContextPtr gcp)
{
	SeqFeatPtr sfp;
	ValNodePtr vnp;
	BioseqPtr bsp;

	if(gcp->thistype != OBJ_SEQFEAT)
		return FALSE;
	sfp = (SeqFeatPtr)(gcp->thisitem);
	if(sfp == NULL || sfp->product == NULL)
		return FALSE;
	vnp  = (ValNodePtr)(gcp->userdata);
	vnp->choice = Seq_mol_aa;
	bsp = BioseqFindCore(SeqLocId(sfp->product));
	if(bsp != NULL)
		vnp->data.ptrvalue = bsp->id;
	else
		vnp->data.ptrvalue = SeqLocId(sfp->product);
	return TRUE;
}
	
static Boolean add_html_label(TextAlignBufPtr tdp)
{
	if(tdp == NULL)
		return FALSE;
	if(tdp->label == NULL)
		return FALSE;
	if(tdp->seqEntityID == 0)
		return FALSE;
	/*only for the coding region features and Bioseqs */
	return ((tdp->feattype == SEQFEAT_CDREGION) ||(tdp->feattype == 0 
		&& tdp->bsp_itemID != 0)); 
	/* return (tdp->feattype == 0); */
}

static SeqIdPtr get_seqid_for_textbuf(TextAlignBufPtr tdp, CharPtr HTML_db,
                                      CharPtr HTML_dopt)
{
	ValNode vn;
	SeqIdPtr sip;


	if(!add_html_label(tdp))
		return 0;

	vn.choice = 0;
	vn.data.ptrvalue = NULL;
	if(tdp->feattype == SEQFEAT_CDREGION)
		GatherItem(tdp->seqEntityID, tdp->itemID, OBJ_SEQFEAT, (Pointer)(&vn), find_seqid_for_seqfeat);
	else
		GatherItem(tdp->seqEntityID, tdp->bsp_itemID, OBJ_BIOSEQ, (Pointer)(&vn), find_seqid_for_bioseq);
	sip = (SeqIdPtr)(vn.data.ptrvalue);
	if(sip != NULL) {
            if(vn.choice == Seq_mol_aa) {
                StringCpy(HTML_dopt, "GenPept");
                StringCpy(HTML_db,   "Protein");
            } else {
                StringCpy(HTML_dopt, "GenBank");
                StringCpy(HTML_db, "Nucleotide");
            }
	}
	return sip;

}


/******************************************************************
*
*	DrawTextToBuffer(tdp_list, m_buf)
*	write the text into a buffer instead of a FILE
*	return the buffer
*
******************************************************************/
static CharPtr DrawTextToBuffer(ValNodePtr tdp_list, CharPtr PNTR m_buf, Boolean is_html, Int4 label_size, Int4 num_size, Boolean compress, Int4Ptr PNTR matrix, Int4 stop_val, Int4 line_len, Boolean show_strand, Boolean strip_semicolon, SeqIdPtr *already_linked)
{
	Boolean already_done;
	TextAlignBufPtr tdp;
	CharPtr docbuf = NULL;
	Int2 i;
	Int4 pos;
	ValNodePtr curr;
	Int4 max_len;	/*maximum length for each line*/
	Int4 size;
	CharPtr HTML_buffer;
	CharPtr matrix_buf;
	Char HTML_db[32], HTML_dopt[16];
	Int4 html_len;
	SeqIdPtr sip;
	DbtagPtr db_tag;
	ObjectIdPtr oip;
	Boolean load;
	Int4 num_empty, max_pos_val;
	Uint1 res;
	Char temp[21];
	Boolean is_first;
	SeqIdPtr seqid_var;

	if(tdp_list==NULL)
		return NULL;
	tdp = (TextAlignBufPtr) tdp_list->data.ptrvalue;
	if(tdp->buf == NULL)
		return NULL;
	if(compress)
	{
		num_empty = label_size + 1 + num_size + 1;
		if(show_strand)
			num_empty += STRAND_SPACE;
	}
	else
	{
		label_size = B_SPACE;
		num_size = POS_SPACE;
		num_empty = B_SPACE + STRAND_SPACE + POS_SPACE + 2;
	}
	/* max_len = 150; */
	max_len = line_len + num_empty + 20;
	if(is_html) {
            Char buffer[1024];
            sprintf(buffer, "%s?cmd=Retrieve&db=&list_uids=&"
                    "dopt=GenPept", NEW_ENTREZ_HREF);
            html_len = StringLen(buffer);
            html_len += 1 + MAX_GI_NUM + 10 + 20 + MAX_DB_NUM;
            HTML_buffer = (CharPtr) MemNew((size_t)html_len * sizeof(Char));
	}
	size = 0;
	max_pos_val = 12;
	for(curr = tdp_list; curr !=NULL; curr = curr->next)
	{
	   tdp = (TextAlignBufPtr) curr->data.ptrvalue;
	   if(tdp->exonCount > 0 && tdp->buf == NULL)	/*it is a codon*/
		size += (3 * max_len);
	   else
	   {
		if(is_html && add_html_label(tdp))
			size += (max_len + html_len);
		else
			size += max_len;
	   }
	   if(tdp->matrix_val)
		   size += max_len;
	}
	if(size == 0)
	{
		if(is_html)
			MemFree(HTML_buffer);
		return NULL;
	}
	size += max_pos_val;
	docbuf = (CharPtr) MemNew((size_t)(size) * sizeof(Char));
	matrix_buf = (CharPtr) MemNew((size_t)max_len * sizeof (Char));

	pos = 0;
	is_first = TRUE;
	while(tdp_list)
	{
	   tdp = (TextAlignBufPtr) tdp_list->data.ptrvalue;
	   if(tdp->exonCount > 0 && tdp->buf == NULL)	/*it is a codon*/
	   {
		for(i =0; i<3; ++i)
		{
			if(i == tdp->frame)
				pos+= print_label_to_buffer_all_ex(docbuf+pos, tdp->label, tdp->pos, 
					tdp->strand, tdp->extra_space, FALSE, label_size, num_size, show_strand, strip_semicolon);
			else
				pos+= print_label_to_buffer_all_ex(docbuf+pos, NULL, -1, 
					0, tdp->extra_space, FALSE, label_size, num_size, show_strand, strip_semicolon);
			sprintf(docbuf+pos, "%s\n", tdp->codon[i]);
			pos += (StringLen(tdp->codon[i]) +1);
		}
	    }
	    if(tdp->exonCount == 0 && tdp->buf !=NULL)
	    {
		if(tdp->matrix_val)	/*print the matrix of the alignment*/
		{
			MemSet(matrix_buf, ' ', (size_t)max_len* sizeof(Char));
			matrix_buf[max_len-1] = '\0';
			size = StringLen(tdp->buf);
			if(matrix) /*protein alignment*/
			{
				for(i = 0; i<size; ++i)
				{
					res = tdp->buf[i];
					if(tdp->matrix_val[i] > 0)
					{
                                            if(tdp->matrix_val[i] == matrix[res][res] || tdp->matrix_val[i] == INT2_MAX)
                                                matrix_buf[i+num_empty] = res;
                                            else
                                                matrix_buf[i+num_empty] = '+';
					}	
				}
			}
			else /*DNA alignment*/
			{
				for(i = 0; i<size; ++i)
					if(tdp->matrix_val[i] != 0)
						matrix_buf[i+num_empty] = (Uint1)(tdp->matrix_val[i]);
			}
			matrix_buf[i+num_empty] = '\0';
			sprintf(docbuf+pos, "%s\n", matrix_buf);
			pos += (StringLen(matrix_buf) +1);
		}
		load = FALSE;
		if(is_html)
		{
			sip = get_seqid_for_textbuf(tdp, HTML_db, HTML_dopt);
			while(!load && sip)
			{
				if(sip->choice == SEQID_GI && sip->data.intvalue != 0)
				{
					seqid_var = *already_linked;
					already_done = FALSE;
					while (seqid_var)
					{
						if (SeqIdMatch(sip, seqid_var) == TRUE)
						{
							already_done = TRUE;
							break;
						}
						seqid_var = seqid_var->next;
					}
					if (already_done) {
       sprintf(HTML_buffer, 
               "<a href=%s?cmd=Retrieve&db=%s&list_uids=%08ld&dopt=%s>", 
               NEW_ENTREZ_HREF, HTML_db, (long)sip->data.intvalue, HTML_dopt);
					} else {

       sprintf(HTML_buffer, "<a name = %ld></a>"
               "<a href=%s?cmd=Retrieve&db=%s&list_uids=%08ld&dopt=%s>", 
               (long) sip->data.intvalue, NEW_ENTREZ_HREF, HTML_db, 
               (long)sip->data.intvalue, HTML_dopt);

       ValNodeAddInt(already_linked, SEQID_GI, sip->data.intvalue);
					}
					html_len = StringLen(HTML_buffer);
					sprintf(docbuf+pos, HTML_buffer);
					pos += html_len;
					pos += print_label_to_buffer_all_ex(docbuf+pos, tdp->label, tdp->pos, 
						tdp->strand, FALSE, TRUE, label_size, num_size, show_strand, strip_semicolon);
					load = TRUE;
				}
				else if(sip->choice == SEQID_GENERAL)
				{
					db_tag = (DbtagPtr) sip->data.ptrvalue;
					if(db_tag->db && StringCmp(db_tag->db, "THC") == 0)
					{
						oip = db_tag->tag;
						if(oip->id != 0)
						{
							sprintf(HTML_buffer, "<a name = THC%ld></a><a href=\"http://www.tigr.org/docs/tigr-scripts/hgi_scripts/thc_report.spl?est=THC%ld&report_type=n\">", (long) oip->id, (long) oip->id);

							html_len = StringLen(HTML_buffer);
							sprintf(docbuf+pos, HTML_buffer);
							pos += html_len;
							pos += print_label_to_buffer_all_ex(docbuf+pos, tdp->label, tdp->pos, 
								tdp->strand, FALSE, TRUE, label_size, num_size, show_strand, strip_semicolon);
							load = TRUE;
						}
					}
				}
				sip = sip->next;
			}
		}

		if(!load)
			pos += print_label_to_buffer_all_ex(docbuf+pos, tdp->label, tdp->pos, tdp->strand, 
				FALSE, FALSE, label_size, num_size, show_strand, strip_semicolon);
		sprintf(docbuf+pos, "%s", tdp->buf);
		pos += StringLen(tdp->buf);
		if(stop_val >=0 && is_first)
		{
			sprintf(temp, " %ld\n", (long) (stop_val+1));
			sprintf(docbuf+pos, "%s", temp);
			pos += StringLen(temp);
		}
		else
		{
			sprintf(docbuf+pos, "\n");
			pos += 1;
		}

		if(m_buf != NULL && *m_buf == NULL)
			*m_buf = StringSave(tdp->buf);

	    }
	    tdp_list = tdp_list->next;
	    is_first = FALSE;
	}

	if(is_html)
		MemFree(HTML_buffer);
	MemFree(matrix_buf);
	docbuf[pos] = '\0';
	return docbuf;
}



/*************************************************************************
*
*	DrawTextList(tdp_list, fp)
*	returns the tdbp->text of the first node. It is used as the master
*	sequence to compare the mismatches
*
*************************************************************************/
/* static CharPtr DrawTextList(ValNodePtr tdp_list, FILE *fp)
{
	TextAlignBufPtr tdp;
	CharPtr m_buf = NULL;
	Int2 i;

	if(tdp_list==NULL)
		return NULL;

	while(tdp_list)
	{
	   tdp = tdp_list->data.ptrvalue;
	   if(tdp->exonCount > 0 && tdp->buf == NULL)	
	   {
		for(i =0; i<3; ++i)
		{
			if(i == tdp->frame)
				print_label(fp, tdp->label, tdp->pos, tdp->strand, tdp->extra_space);
			else
				print_label(fp, NULL, -1, 0, tdp->extra_space);
			fprintf(fp, "%s\n", tdp->codon[i]);
		}
	    }
	    if(tdp->exonCount == 0 && tdp->buf !=NULL)
	    {
		print_label(fp, tdp->label, tdp->pos, tdp->strand, FALSE);
		fprintf(fp, "%s\n", tdp->buf);
		if(m_buf == NULL)
			m_buf = tdp->buf;
	    }
	    tdp_list = tdp_list->next;
	}

	return m_buf;
} */

NLM_EXTERN Boolean make_scale_bar_str(CharPtr PNTR bar, CharPtr PNTR num_str, 
		Int4 num_empty, Int4 line_len)
{
	Int4 i, j;
	CharPtr curr;
	Char temp[100];
	Int4 len;

	if(bar == NULL || num_str == NULL)
		return FALSE;
	*bar = (CharPtr) MemNew((size_t)(line_len+num_empty+2) * sizeof(Char));
	*num_str = (CharPtr) MemNew((size_t)(line_len+num_empty+2) * sizeof(Char));
	make_empty(*bar, (Int2)(line_len+num_empty));
	make_empty(*num_str, (Int2)(line_len+num_empty));
	for(i =0; i<line_len; ++i)
	{
		curr = *bar;
		if((i+1)%5 ==0)
			curr[i+num_empty]= '|';
		curr = *num_str;
		if((i+1)%10==0)
		{
			sprintf(temp, "%ld", (long) (i+1));
			len = StringLen(temp);
			for(j = 0; j<len; ++j)
				curr[i+num_empty-(len-1-j)] = temp[j];
		}
		
	}
	return TRUE;
}

static AlignNodePtr get_master_align_node(ValNodePtr anp_list)
{
	AlignNodePtr anp, first_anp = NULL;

	while(anp_list)
	{
		if(anp_list->choice != OBJ_SEQANNOT)
		{
			anp = (AlignNodePtr) anp_list->data.ptrvalue;
			if(anp->is_master)
				return anp;
			else if(first_anp == NULL)
				first_anp = anp;
		}
		anp_list = anp_list->next;
	}
	return first_anp;
}
	

/**********************************************************************
*
*	figure out the DNA-protein alignment, the frame of translation 
*	in the DNA sequence compared with the protein sequence
*	m_len is used to figure out the reading frame of the minus 
*	strand translation
*
**********************************************************************/
static Boolean load_master_translate_frame(ValNodePtr anp_list, Int4 m_len, BioseqPtr m_bsp)
{
	AlignNodePtr anp;
	AlignSegPtr asp;
	Int4 g_left, g_right;
	Int4 start_pos;
	Uint1 strand;
	Int4 val;
	AnnotInfoPtr annot_info;
	Uint1 align_type;
	AlignNodePtr master_anp;
	Int4 offset;
	Boolean found;


	master_anp = get_master_align_node(anp_list);
	if(master_anp == NULL)
		return FALSE;

	while(anp_list)
	{
		align_type = 0;
		annot_info = NULL;
		while(anp_list != NULL)
		{
			if(anp_list->choice == OBJ_SEQANNOT)
			{
				annot_info = (AnnotInfoPtr) anp_list->data.ptrvalue;
				align_type = get_alignment_type(annot_info);
			}
			else
			{
				if(annot_info == NULL)
					break;
				else 
				{
					if(align_type != 0 || annot_info->annot_type == ANNOT_BLAST)
						break;
				}
			}
			anp_list = anp_list->next;
		}
	
		while(anp_list != NULL && anp_list->choice != OBJ_SEQANNOT)
		{
			if(align_type == ALIGN_DNA_TO_PROT || align_type == ALIGN_TDNA_TO_TDNA)
			{
				anp = (AlignNodePtr) anp_list->data.ptrvalue;
				/*for tblastx, need to figure out the translation frame of the master, 
				   for blastx (reverted to tblastn) and tblastn, need to figure out
				   the frame for the query
				*/
				if((align_type == ALIGN_TDNA_TO_TDNA && BioseqMatch(m_bsp, anp->sip))|| 
					(!(anp->is_master) && !BioseqMatch(m_bsp, anp->sip)))
				{
					g_left = anp->extremes.left;
					g_right = anp->extremes.right;
					if(align_type == ALIGN_TDNA_TO_TDNA)
					{
						strand = anp->extremes.strand;
						offset = 0;
						found = TRUE;
					}
					else
					{
						strand = Seq_strand_plus;
						if(anp->extremes.strand != master_anp->extremes.strand)
						{
							if(anp->extremes.strand == Seq_strand_minus || 
								master_anp->extremes.strand == Seq_strand_minus)
									strand = Seq_strand_minus;
						}
						if(anp->extremes.strand == Seq_strand_minus)
							g_left = g_right;

						offset = 0;
						found = FALSE;
						for(asp = master_anp->segs; asp != NULL; asp = asp->next)
						{
							if(asp->type != GAP_SEG)
							{
								if(asp->type == INS_SEG)
								{
									if(offset > 0)
										offset += asp->gr.right;
								}
								else
								{
									if(asp->gr.right < g_left)
										offset += (asp->gr.right - asp->gr.left + 1);
									else
									{
										if(g_left >= asp->gr.left && g_left <= asp->gr.right)
										{
											offset += MAX(0, g_left - asp->gr.left);
											found = TRUE;
											break;
										}
									}
								}
							}
						}
					}
	
					if(found)
					{
						start_pos = ABS(master_anp->seqpos + offset);
						if(strand == Seq_strand_minus)
						{
							val = (m_len -1 - start_pos)%3L;
						}
						else
						{
							val = start_pos%3L;
						}

						switch(val)
						{
							case 0:
								if(strand == Seq_strand_minus)
									anp->m_frame = 4;
								else
									anp->m_frame = 1;
								break;
							case 1:
								if(strand == Seq_strand_minus)
									anp->m_frame = 5;
								else
									anp->m_frame = 2;
								break;
							case 2:
								if(strand == Seq_strand_minus)
									anp->m_frame = 6;
								else
									anp->m_frame = 3;
								break;
							default:
								break;
						}
					}
				}
			}
			anp_list = anp_list->next;
		}
	}
	return TRUE;
}


/**********************************************************************
*
*	figure out in the current range (m_left, m_right), the total number
*	of reading frames that the hit proteins have
*	return the AlignNode for the master sequence (master is always the 
*	DNA sequence
*
***********************************************************************/											
static Boolean get_current_master_frame(ValNodePtr list, Int4 m_left, Int4 m_right, Uint1Ptr all_frame)
{
	ValNodePtr anp_list;
	AlignNodePtr anp;
	Int4 g_left, g_right;
	Uint1 i;
	Boolean retval;

	MemSet((Pointer)all_frame, 0, (size_t)6 * sizeof(Uint1));
	while(list)
	{
		anp_list = (ValNodePtr) list->data.ptrvalue;
		anp = (AlignNodePtr) anp_list->data.ptrvalue;
		if(!(anp->is_master))
		{
			g_left = anp->extremes.left;
			g_right = anp->extremes.right;
			if(!(g_left > m_right || g_right < m_left))
			{
				if(anp->m_frame != 0)
				{
					for(i = 0; i<6; ++i)
					{
						if(all_frame[i] == anp->m_frame)
							break;
						else if(all_frame[i] == 0)
						{
							all_frame[i] = anp->m_frame;
							break;
						}
					}
					retval = TRUE;
				}
			}
		}

		list = list->next;
	}

	return retval;
}

NLM_EXTERN SeqFeatPtr make_fake_cds(BioseqPtr m_bsp, Int4 start, Int4 stop, Uint1 strand)
{
	SeqFeatPtr sfp;
	CdRegionPtr crp;
	IntFuzzPtr ifp_from, ifp_to;
	Uint1 g_code = 0;
	SeqDescrPtr descr;
	SeqIntPtr seq_int;
	SeqLocPtr slp;
	BioSourcePtr source;
	OrgRefPtr org;
	OrgNamePtr orgname;
	ValNodePtr vnp;

	descr = m_bsp->descr;
	while(descr)
	{
		/*look into BioSource to get the genetic code*/
		if(descr->choice == Seq_descr_source)
		{
			source = (BioSourcePtr) descr->data.ptrvalue;
			if(source != NULL)
			{
				org = source->org;
				if(org != NULL)
				{
					orgname = org->orgname;
					if(orgname != NULL)
					{
						g_code = orgname->gcode;
						break;
					}
				}
			}
		}
		descr = descr->next;
	}

	crp = CdRegionNew();
	if(g_code != 0)
	{
		vnp = ValNodeNew(NULL);
		vnp->choice = 2;
		vnp->data.intvalue = (Int4)g_code;
		ValNodeAddPointer(&(crp->genetic_code), 254, (Pointer)vnp);
	}

	sfp = SeqFeatNew();
	sfp->data.choice = 3;
	sfp->data.value.ptrvalue = crp;
	sfp->partial = TRUE;
	sfp->product = NULL;
	slp = SeqLocIntNew(start, stop, strand, m_bsp->id);
	seq_int = (SeqIntPtr) slp->data.ptrvalue;
	ifp_from = IntFuzzNew();
	ifp_from->choice = 4;
	ifp_from->a = 2;
	seq_int->if_from = ifp_from;
	ifp_to = IntFuzzNew();
	ifp_to->choice = 4;
	ifp_to->a = 1;
	seq_int->if_to = ifp_to;
	sfp->location = slp;
	

	return sfp;
}


static CharPtr translate_faked_cds(SeqFeatPtr fake_cds, Uint1 frame, Int4 c_start, Int4 c_stop, Int4 master_len, AlignNodePtr anp)
{
	Uint1 c_frame;
	SeqLocPtr slp, t_slp;
	SeqIntPtr sint;
	Uint1 strand;
	CdRegionPtr crp;
	CharPtr buf;
	Int4 from, to;
	Int4 n;
	AlignSegPtr asp;
	Int4 m_start, m_stop;
	Int4 t_start, t_stop;
	Int4 l_pos = 0;
	Int4 offset;


	buf = (CharPtr) MemNew((size_t)(c_stop - c_start+2) * sizeof(Char));
	buf[0] = '\0';
	offset = 0;
	for(asp = anp->segs; asp != NULL; asp = asp->next)
	{
		if(!((asp->gr.left > c_stop) || ( asp->gr.right < c_start)))
		{
			t_start = MAX(asp->gr.left, c_start);
			t_stop = MIN(asp->gr.right, c_stop);
			m_start = ABS(anp->seqpos + t_start - offset - anp->extremes.left);
			m_stop = ABS(anp->seqpos + t_stop - offset - anp->extremes.left);
			if(asp->type == GAP_SEG)
			{
				MemSet(buf+l_pos, ' ', (size_t)(t_stop - t_start + 1) * sizeof (Char));
				l_pos += t_stop - t_start + 1;
				buf[l_pos] = '\0';
			}
			else if(asp->type != INS_SEG)
			{
				if(frame > 3)
				{
					strand = Seq_strand_minus;
					n = MAX(0, (master_len-1 - m_stop)/3 -1);
					from = m_start;
					to = master_len -1 - n*3;
					c_frame = frame -3;
					from = MAX(0, from -3);
				}
				else
				{
					strand = Seq_strand_plus;
					n = MAX(0, (m_start/3-1));
					from = n * 3;
					to = m_stop;
					c_frame = frame;
					to = MIN(master_len-1, to + 3);
				}
				slp = fake_cds->location;
				sint = (SeqIntPtr) slp->data.ptrvalue;
				sint->from = from;
				sint->to = to;
				sint->strand = strand;
				crp = (CdRegionPtr) fake_cds->data.value.ptrvalue;
				crp->frame = c_frame;
			
				t_slp = SeqLocIntNew(m_start, m_stop, Seq_strand_plus, SeqLocId(slp));
				print_protein_for_cds(fake_cds, buf+l_pos, t_slp, TRUE);
				l_pos += t_stop - t_start + 1;
				SeqLocFree(t_slp);
			}
		}
		else if(asp->gr.left > c_stop)
			break;
		if(asp->type == GAP_SEG)
			offset += (asp->gr.right - asp->gr.left +1 );
	}
	return buf;
}

static ValNodePtr load_fake_protein_buf(CharPtr buf, Uint1 frame, AlignNodePtr master_anp)
{
	Char temp[20];
	ValNodePtr head = NULL;
	TextAlignBufPtr tdp;

	tdp = (TextAlignBufPtr) MemNew(sizeof(TextAlignBuf));
	tdp->pos = -1;
	if(frame >3)
	{
		tdp->strand = Seq_strand_minus;
		sprintf(temp, "frame=-%d", frame-3);
	}
	else
	{
		tdp->strand = Seq_strand_plus;
		sprintf(temp, "frame=+%d", frame);
	}
	tdp->label = StringSave(temp);
	tdp->buf = StringSave(buf);
	tdp->itemID = 0;
	tdp->feattype = 0;
	tdp->subtype = 0;
	tdp->entityID = master_anp->entityID;
	tdp->seqEntityID = master_anp->seq_entityID;
	tdp->bsp_itemID = master_anp->bsp_itemID;
	ValNodeAddPointer(&head, 0, tdp);
	return head;
}

static Boolean has_hit_in_region(Uint1Ptr all_frame)
{
	Uint1 i;

	for(i = 0; i<6; ++i)
		if(all_frame[i] != 0)
			return TRUE;
	return FALSE;
}

	
static Boolean check_bsp_id(ValNodePtr PNTR id_list, SeqIdPtr sip)
{
	ValNodePtr curr;

	if(*id_list == NULL)
	{
		ValNodeAddPointer(id_list, 0, sip);
		return FALSE;
	}
	curr = *id_list;
	while(curr)
	{
		if(curr->data.ptrvalue == sip || 
			SeqIdMatch((SeqIdPtr)(curr->data.ptrvalue), sip))
			return TRUE;
		curr = curr->next;
	}

	
	ValNodeAddPointer(id_list, 0, sip);
	return FALSE;
}

static Boolean find_align_proc(GatherContextPtr gcp)
{
	SeqAlignPtr PNTR p_align;

	p_align = (SeqAlignPtr PNTR)(gcp->userdata);
	*p_align = (SeqAlignPtr)gcp->thisitem;
	return TRUE;
}


/*functions related to load the alignment summary, such as the number of 
identical, positive residues, # of gaps, to the structure
*/ 

static Boolean load_align_sum_for_DenseDiag(DenseDiagPtr ddp, AlignSumPtr asp)
{
    SeqInt si;
    SeqLoc sl;
    Int4 i;
    Int2 m_order, t_order;	/*order of the master and the target sequence*/
    Uint1 m_res, t_res;
    SeqIdPtr sip;
    SeqPortPtr m_spp, t_spp;
    Int2 dim;
    SeqPortPtr spp;
    
    if(ddp == NULL || asp == NULL)
        return FALSE;
    m_order = -1;
    t_order = -1;
    dim = 0;
    for(i = 0, sip = ddp->id; sip != NULL; sip = sip->next, ++i) {
        if(SeqIdMatch(sip, asp->master_sip) && m_order == -1)
            m_order = i;
        else if(SeqIdMatch(sip, asp->target_sip) && t_order == -1)
            t_order = i;
        ++dim;
    }
    
    if(m_order == -1 || t_order == -1)
        return FALSE;
    
    asp->m_frame_set = FALSE;
    asp->t_frame_set = FALSE;
    
    for(i = 0; i<dim; ++i) {
        if(i == m_order || i == t_order) {
            
            if(i == m_order)
                si.id = asp->master_sip;
            else
                si.id = asp->target_sip;
            si.from = ddp->starts[i];
            si.to = si.from + ddp->len -1;
            if(ddp->strands != NULL)
                si.strand = ddp->strands[i];
            else
                si.strand = 0;
            
            
            if (asp->is_aa) {
                asp->m_strand = Seq_strand_unknown;
                asp->t_strand = Seq_strand_unknown;
            } else {
                if(i == m_order) {
                    asp->m_strand = si.strand;
                } else {
                    asp->t_strand = si.strand;
                }
            }
            
            sl.choice = SEQLOC_INT;
            sl.data.ptrvalue = &si;
            
            spp = SeqPortNewByLoc(&sl, (asp->is_aa) ? Seq_code_ncbieaa : Seq_code_iupacna);
            if(i == m_order) {
                asp->master_from = si.from;
                asp->master_to = si.to;
                m_spp = spp;
            } else {
                asp->target_from = si.from;
                asp->target_to = si.to;
                t_spp = spp;
            }
        }
    }
    
    if(m_spp == NULL || t_spp == NULL) {
        if(m_spp == NULL)
            SeqPortFree(m_spp);
        if(t_spp != NULL)
            SeqPortFree(t_spp);
        return FALSE;
    }
    
    for(i = 0; i<ddp->len; ++i) {
        m_res = SeqPortGetResidue(m_spp);
        t_res = SeqPortGetResidue(t_spp);
        if(m_res == t_res)
            ++(asp->identical);
        else if(asp->matrix != NULL && asp->is_aa) {
            if (IS_residue(m_res) && IS_residue(t_res)) {
                if(asp->matrix[m_res][t_res] >0)
                    ++(asp->positive);
            }
        }
    }
    asp->totlen = ddp->len;
    
    SeqPortFree(m_spp);
    SeqPortFree(t_spp);
    return TRUE;
}


static Boolean load_align_sum_for_DenseSeg(DenseSegPtr dsp, AlignSumPtr asp)
{
    SeqInt msi, tsi;
    SeqIntPtr sint;
    SeqLoc sl;
    Int2 i, k;
    Int2 dim;
    Int2 m_order, t_order;	/*order of the master and the target sequence*/
    Int4 index;
    Int4 j, val, t_val;
    Uint1 m_res, t_res, stdaa_res;
    SeqIdPtr sip;
    SeqPortPtr m_spp, t_spp;
    SeqMapTablePtr smtp;
    
    
    if(dsp == NULL || asp == NULL)
        return FALSE;
    
    if(asp->posMatrix != NULL) { 
        if((smtp = SeqMapTableFindObj(Seq_code_ncbistdaa, 
                                      Seq_code_ncbieaa)) == NULL)
            return FALSE;
    }
    
    m_order = -1;
    t_order = -1;
    dim = 0;
    for(i = 0, sip = dsp->ids; sip != NULL; sip = sip->next, ++i) {
        if(SeqIdMatch(sip, asp->master_sip) && m_order == -1)
            m_order = i;
        else if(SeqIdMatch(sip, asp->target_sip) && t_order == -1)
            t_order = i;
        ++dim;
    }
    
    if(m_order == -1 || t_order == -1)
        return FALSE;
    
    msi.id = asp->master_sip;
    msi.from = -1;
    msi.to = -1;
    msi.strand = (dsp->strands == NULL) ? 0 : dsp->strands[m_order];
    
    tsi.id = asp->target_sip;
    tsi.from = -1;
    tsi.to = -1;
    tsi.strand = (dsp->strands == NULL) ? 0 : dsp->strands[t_order];
    
    for(i = 0; i<dsp->numseg; ++i) {
        for(k = 0; k<dim; ++k) {
            val = dsp->starts[i*dim + k];
            if(val != -1 && (k == m_order || k == t_order)) {
                sint = (k == m_order) ? (&msi) : (&tsi);
                if(sint->from == -1 || sint->from > val)
                    sint->from = val;
                if(sint->to == -1 || sint->to < (val + dsp->lens[i] -1))
                    sint->to = val + dsp->lens[i] -1;
            }
        }
    }

    asp->master_from = msi.from;
    asp->master_to = msi.to;
    asp->target_from = tsi.from;
    asp->target_to = tsi.to;
    
    if (asp->is_aa) {
        asp->m_strand = Seq_strand_unknown;
        asp->t_strand = Seq_strand_unknown;
    } else {
        asp->m_strand = dsp->strands[m_order];
        asp->t_strand = dsp->strands[t_order];
    }
    asp->m_frame_set = FALSE;
    asp->t_frame_set = FALSE;
    
    sl.choice = SEQLOC_INT;
    sl.data.ptrvalue = &msi;
    m_spp = SeqPortNewByLoc(&sl, (asp->is_aa) ? Seq_code_ncbieaa : Seq_code_iupacna);
    
    sl.choice = SEQLOC_INT;
    sl.data.ptrvalue = &tsi;
    t_spp = SeqPortNewByLoc(&sl, (asp->is_aa) ? Seq_code_ncbieaa : Seq_code_iupacna);
    
    for(i = 0; i<dsp->numseg; ++i) {
        val = dsp->starts[i*dim + m_order];
        t_val = dsp->starts[i*dim + t_order];
        if(val == -1 || t_val == -1) {
            asp->gaps += dsp->lens[i];
            if(val != -1) {
                index = dsp->lens[i];
                while (index > 0) {
                    index--;
                    m_res = SeqPortGetResidue(m_spp);
                }
            }
            if(t_val != -1) {
                index = dsp->lens[i];
                while (index > 0) {
                    index--;
                    t_res = SeqPortGetResidue(t_spp);
                }
            }
        } else {
            for(j = 0; j<dsp->lens[i]; ++j) {
                m_res = SeqPortGetResidue(m_spp);
                t_res = SeqPortGetResidue(t_spp);
                if(m_res == t_res)
                    ++(asp->identical);
                else if(asp->matrix != NULL && asp->is_aa) {
                    if (IS_residue(m_res) && IS_residue(t_res)) {

                        if(asp->posMatrix != NULL) {
                            stdaa_res = SeqMapTableConvert(smtp, t_res);
                            if(asp->posMatrix[val+j][stdaa_res] > 0)
                                ++(asp->positive);
                        } else {
                            if(asp->matrix[m_res][t_res] >0)
                                ++(asp->positive);
                        }
                    }
                }
            }
            
        }
        asp->totlen += dsp->lens[i];
    }
    SeqPortFree(m_spp);
    SeqPortFree(t_spp);
    return TRUE;
}

/*
	Obtains the genetic code from a BioseqPtr, assuming that a fetch function
	has been enabled.
*/
NLM_EXTERN CharPtr
GetGeneticCodeFromSeqId (SeqIdPtr sip)

{
	BioseqPtr bsp;
	BioSourcePtr source;
	CharPtr genetic_code=NULL;
	GeneticCodePtr gcp;
	Int4 gen_code_val=1;	/* std genetic code if nothing found. */
	ValNodePtr vnp;


	bsp = BioseqLockById(sip);

	if (bsp)
	{
		vnp = BioseqGetSeqDescr(bsp, Seq_descr_source, NULL);
		if (vnp)
		{
			source = (BioSourcePtr) vnp->data.ptrvalue;
			gen_code_val = source->org->orgname->gcode;
		}
		BioseqUnlock(bsp);
	}

	gcp = GeneticCodeFind(gen_code_val, NULL);
	for (vnp = (ValNodePtr)gcp->data.ptrvalue; vnp != NULL; vnp = vnp->next)
	{
		if (vnp->choice == 3)   /* ncbieaa */
		{
			genetic_code = (CharPtr)vnp->data.ptrvalue;
			break;
		}
	}

	return genetic_code;
}
NLM_EXTERN CharPtr OOFTranslateDNAInAllFrames(Uint1Ptr dna, Int4 length, 
                                              SeqIdPtr query_id)
{
    CharPtr dnap;
    CharPtr codes;
    Int4 i;
    Uint1 codon[3];

    if(dna == NULL || length == 0)
        return NULL;
    
    dnap = (CharPtr) MemNew(length+1);
    codes = GetGeneticCodeFromSeqId(query_id);

    dnap[0] = dnap[1] = dnap[2] = 0;
    
    for (i = 2; i < length; i++) {
        codon[0] = dna[i-2];
        codon[1] = dna[i-1];
        codon[2] = dna[i];        
        dnap[i+1] = AAForCodon(codon, codes);
    }
    return dnap;
}

NLM_EXTERN Uint1 AAForCodon (Uint1Ptr codon, CharPtr codes);

static Boolean load_align_sum_for_StdSeg(StdSegPtr ssp, AlignSumPtr asp)
{
    Boolean master_is_translated=FALSE, both_translated=FALSE;
    Boolean target_is_translated = FALSE;
    BioseqPtr bsp;
    CharPtr genetic_code1, genetic_code2;
    SeqPortPtr spp1, spp2;
    Uint1 codon[4], residue1, residue2;
    Boolean ungapped_align = FALSE;
    StdSegPtr ssp_last;
    
    if(ssp == NULL || asp == NULL)
        return FALSE;

    if(asp->ooframe) {
        master_is_translated = TRUE;
        target_is_translated = FALSE;	
    } else {
        /* Check for valid sequence. */
        if (SeqLocLen(ssp->loc) == 3*SeqLocLen(ssp->loc->next))
            master_is_translated = TRUE;
        else if (3*SeqLocLen(ssp->loc) == SeqLocLen(ssp->loc->next))
            target_is_translated = TRUE;	
        else if (SeqLocLen(ssp->loc) == SeqLocLen(ssp->loc->next))
            both_translated = TRUE;
        else
            return FALSE;
    }
    
    if (master_is_translated) {
        genetic_code1 = GetGeneticCodeFromSeqId(ssp->ids);
    } else if (both_translated) {
        genetic_code1 = GetGeneticCodeFromSeqId(ssp->ids);
        genetic_code2 = GetGeneticCodeFromSeqId(ssp->ids->next);
    } else {
        genetic_code1 = GetGeneticCodeFromSeqId(ssp->ids->next);
    }
    
    asp->m_frame_set = FALSE;
    asp->t_frame_set = FALSE;
    
    if (master_is_translated || both_translated) {
        asp->m_strand = SeqLocStrand(ssp->loc);
        asp->m_frame = SeqLocStart(ssp->loc);
        if (SeqLocStrand(ssp->loc) == Seq_strand_minus) {
            bsp = BioseqLockById(SeqLocId(ssp->loc));
            asp->m_frame += SeqLocLen(ssp->loc);
            asp->m_frame = -(1+(bsp->length - asp->m_frame)%3);
            asp->m_frame_set = TRUE;
            BioseqUnlock(bsp);
        } else {
            asp->m_frame = (1+(asp->m_frame)%3);
            asp->m_frame_set = TRUE;
        }
    }
    
    if (!master_is_translated || both_translated) {
        asp->t_strand = SeqLocStrand(ssp->loc->next);
        asp->t_frame = SeqLocStart(ssp->loc->next);
        if (SeqLocStrand(ssp->loc->next) == Seq_strand_minus) {
            if (bsp = BioseqLockById(SeqLocId(ssp->loc->next))) {
                asp->t_frame += SeqLocLen(ssp->loc->next);
                asp->t_frame = -(1+(bsp->length - asp->t_frame)%3);
                asp->t_frame_set = TRUE;
                BioseqUnlock(bsp);
            } else {
                return FALSE;
            }
        } else {
            asp->t_frame = (1+(asp->t_frame)%3);
            asp->t_frame_set = TRUE;
        }
    }
    
    if (SeqLocStrand(ssp->loc) == Seq_strand_minus) {
        asp->master_from = SeqLocStop(ssp->loc);
    } else {
        asp->master_from = SeqLocStart(ssp->loc);
    }
    
    if (SeqLocStrand(ssp->loc->next) == Seq_strand_minus) {
        asp->target_from = SeqLocStop(ssp->loc->next);
    } else {
        asp->target_from = SeqLocStart(ssp->loc->next);
    }
    
    
    while (ssp) {
        if (ssp->loc->choice != SEQLOC_EMPTY && ssp->loc->next->choice != SEQLOC_EMPTY) {
            if (both_translated) {
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
                    if (residue1 == residue2)
                        ++(asp->identical);
                    else if (asp->matrix != NULL && asp->matrix[residue1][residue2] >0)
                        ++(asp->positive);
                }
            } else {
                if (master_is_translated) {
                    spp1 = SeqPortNewByLoc(ssp->loc, Seq_code_ncbi4na);
                    spp2 = SeqPortNewByLoc(ssp->loc->next, Seq_code_ncbieaa);
                } else {
                    spp2 = SeqPortNewByLoc(ssp->loc, Seq_code_ncbieaa);
                    spp1 = SeqPortNewByLoc(ssp->loc->next, Seq_code_ncbi4na);
                }
                
                while ((residue1=SeqPortGetResidue(spp2)) != SEQPORT_EOF) {
                    codon[0] = SeqPortGetResidue(spp1);
                    codon[1] = SeqPortGetResidue(spp1);
                    codon[2] = SeqPortGetResidue(spp1);
                    residue2 = AAForCodon(codon, genetic_code1);
                    if (residue1 == residue2)
                        ++(asp->identical);
                    else if (asp->matrix != NULL && asp->matrix[residue1][residue2] >0)
                        ++(asp->positive);
                }
            }
            SeqPortFree(spp1);
            SeqPortFree(spp2);
            /* Check if this is an ungapped alignment;
               in this case do not go to next link */
            if (!asp->ooframe && ssp->next && 
                ssp->next->loc->choice != SEQLOC_EMPTY && 
                ssp->next->loc->next->choice != SEQLOC_EMPTY)
                ungapped_align = TRUE;
        } else {	/* Count only gaps in the top (master) strand. */
            if (ssp->loc->choice == SEQLOC_EMPTY)
                {
                    if (!master_is_translated || both_translated)
                        asp->gaps += SeqLocLen(ssp->loc->next)/3;
                    else
                        asp->gaps += SeqLocLen(ssp->loc->next);
                }
        }

        if(asp->ooframe) {
            if(ssp->loc->next->choice != SEQLOC_EMPTY)
                asp->totlen += SeqLocLen(ssp->loc->next);
            else
                asp->totlen += SeqLocLen(ssp->loc)/3;
        } else {
        
            if (ssp->loc->choice != SEQLOC_EMPTY) {
                if (master_is_translated || both_translated)
                    asp->totlen += SeqLocLen(ssp->loc)/3;
                else
                    asp->totlen += SeqLocLen(ssp->loc);
            } else {
                if (target_is_translated || both_translated)
                    asp->totlen += SeqLocLen(ssp->loc->next)/3;
                else
                    asp->totlen += SeqLocLen(ssp->loc->next);
            }
        }

        ssp_last = ssp;
        
        if (both_translated || ungapped_align)	
            /* for tblastx perform only one StdSegPtr. */
            break;
        
        ssp = ssp->next;
    }

    if (SeqLocStrand(ssp_last->loc) == Seq_strand_minus) {
        asp->master_to = SeqLocStart(ssp_last->loc);
    } else {
        asp->master_to = SeqLocStop(ssp_last->loc);
    }
    
    if (SeqLocStrand(ssp_last->loc->next) == Seq_strand_minus) {
        asp->target_to = SeqLocStart(ssp_last->loc->next);
    } else {
        asp->target_to = SeqLocStop(ssp_last->loc->next);
    }
    
    return TRUE;
}


/*****************************************************************
*
*	find_score_in_align(align, chain, asp)
*	align: the Seq-align point
*	chain: for multiple segment Seq-aligns, such as DenseDiag and 
*	StdSeg, the order within the Seq-align
*	asp:   the structure that records and stores the positive, 
*			identical residues
*
*****************************************************************/
NLM_EXTERN ScorePtr find_score_in_align(SeqAlignPtr align, Uint2 chain, 
                                        AlignSumPtr asp)
{
    DenseDiagPtr ddp;
    DenseSegPtr dsp;
    StdSegPtr ssp;
    Uint2 order = 0;
    SeqAlignPtr sap;
    ScorePtr    sp;
    
    if(align == NULL)
        return NULL;
    
    if(asp != NULL) {
        asp->totlen = 0;
        asp->positive = 0;
        asp->identical = 0;
        asp->gaps = 0;
    }
    switch (align->segtype) {
    case 1: /*Dense-diag*/
        ddp = (DenseDiagPtr) align->segs;
        while(ddp) {
            ++order;
            if(order == chain) {
                if(asp != NULL) 
                    load_align_sum_for_DenseDiag(ddp, asp);
                return ddp->scores;
            }
            ddp = ddp->next;
        }
        break;
    case 2:
        dsp = (DenseSegPtr) align->segs;
        if(asp != NULL)
            load_align_sum_for_DenseSeg(dsp, asp);
        if (dsp->scores)
            return dsp->scores;
        else
            return align->score;
    case 3:
        ssp = (StdSegPtr) align->segs;
        while(ssp) {
            ++order;
            if(order == chain) {
                if(asp != NULL) 
                    load_align_sum_for_StdSeg(ssp, asp);
                return ssp->scores;
            }
            ssp = ssp->next;
        }
        break;
    case 5: /* Discontinuous alignment */
        sap = (SeqAlignPtr) align->segs;
        
        if((sp = find_score_in_align(sap, chain, asp)) == NULL)
            return align->score;
        else
            return sp;
    default:
        break;
    }
    return NULL;
}

/*try to decide if this fit the prototype of reversing the BLASTX 
  result to make a TBLASTN output
*/
static Boolean reverse_blastx_order (ValNodePtr anp_list)
{
	Int2 num = 0;
	ValNodePtr c_list;
	AnnotInfoPtr annot_info;
	Uint1 align_type = 0;

	for(c_list = anp_list; c_list != NULL; c_list = c_list->next)
	{
		if(c_list->choice == OBJ_SEQANNOT)
		{
			annot_info = (AnnotInfoPtr) c_list->data.ptrvalue;
			align_type = get_alignment_type(annot_info);
			if(align_type != ALIGN_DNA_TO_PROT)
				return FALSE;
		}
		else
		{
			++num;
			if(num > 2)	/*only for pairwise alignment */
				return FALSE;
		}
	}

	/* return (align_type == ALIGN_DNA_TO_PROT && num == 2); */
	return (align_type == ALIGN_DNA_TO_PROT);
}

/*
* change the alignnode of blastx to a pseudo tblastn and switch the 
* master sequence so that the blastx display will be the same as the 
* traditional blastx
*	
*/
static void modify_gather_range (GatherRangePtr grp, Boolean expand)
{
	Int4 len;

	len = grp->right - grp->left + 1;
	if(expand)
	{
		grp->left *= 3;
		grp->right = grp->left + len * 3 -1;
	}
	else
	{
		grp->left /=3;
		grp->right = grp->left + len/3 -1;
	}
}

static Boolean change_blastx_master(ValNodePtr anp_list, AlignNodePtr PNTR master_anp)
{
	
	ValNodePtr c_list;
	AnnotInfoPtr annot_info = NULL;
	AlignNodePtr anp, m_anp, t_anp;
	AlignSegPtr asp;

	anp = NULL;
	m_anp = NULL;
	annot_info = NULL;
	for(c_list = anp_list; c_list != NULL; c_list = c_list->next)
	{
		if(c_list->choice == OBJ_SEQANNOT)
		{
			annot_info = (AnnotInfoPtr) c_list->data.ptrvalue;
			if(annot_info != NULL && 
				get_alignment_type(annot_info) != ALIGN_DNA_TO_PROT)
				annot_info = NULL;
		}
		else
		{
			t_anp = (AlignNodePtr) c_list->data.ptrvalue;
			if(t_anp->is_master || t_anp == *master_anp)
				m_anp = t_anp;
			else
				anp = t_anp;
		}
	}

	if(m_anp == NULL || anp == NULL || annot_info == NULL)
		return FALSE;
				
	/*shrink the interval */
	for(c_list = anp_list; c_list != NULL; c_list = c_list->next)
	{
		if(c_list->choice != OBJ_SEQANNOT)
		{
			t_anp = (AlignNodePtr) c_list->data.ptrvalue;
			modify_gather_range (&(t_anp->extremes), FALSE);
			for(asp = t_anp->segs; asp != NULL; asp = asp->next)
			{
				if(asp->type != INS_SEG)
					modify_gather_range (&(asp->gr), FALSE);
			}
		}
	}
			
	annot_info->blast_type = ALIGN_TBLASTN;
	*master_anp = anp;
	anp->is_master = TRUE;
	m_anp->is_master = FALSE;
	return TRUE;
}

static Int4 get_max_feature_label (AlignNodePtr anp)
{
	AlignSegPtr asp;
	Int4 len = 0, f_len;
	FeatNodePtr fnp;
	ValNodePtr vnp;

	for(asp = anp->segs; asp != NULL; asp = asp->next)
	{
		if(asp->type != GAP_SEG && asp->type != INS_SEG)
		{
			for(vnp = asp->cnp; vnp != NULL; vnp = vnp->next)
			{
				fnp = (FeatNodePtr) vnp->data.ptrvalue;
				if(fnp !=NULL && fnp->label != NULL)
				{
					f_len = StringLen(fnp->label);
					if(f_len > len)
						len = f_len;
				}
			}
		}
	}

	return len;
}


/*
*
*	look through the list of alignnode to figure out the maximum 
*	length required to print the coordinates of the sequence in 
*	alignment
*
*/ 
static Int4 get_max_coordinates_len (ValNodePtr anp_list, Int4Ptr max_label_size)
{
	AlignNodePtr anp;
	AlignSegPtr asp;
	Int4 max_num, seqpos;
	Char buf[101];
	Int4 flabel_len;

	max_num = 0;
	*max_label_size = 0;
	while(anp_list)
	{
		if(anp_list->choice != OBJ_SEQANNOT)
		{
			anp = (AlignNodePtr) anp_list->data.ptrvalue;
			if(anp->seqpos < 0)
				max_num = MAX(ABS(anp->seqpos), max_num);
			else
			{
				seqpos = anp->seqpos - 1;
				for(asp = anp->segs; asp != NULL; asp = asp->next)
				{
					if(asp->type != GAP_SEG)
					{
						if(asp->type == INS_SEG)
							seqpos += (asp->gr.right -1);
						else
							seqpos += (asp->gr.right - asp->gr.left + 1);
					}
				}
				max_num = MAX(seqpos, max_num);
			}
			if(anp->label != NULL)
				*max_label_size = MAX(*max_label_size, (Int4)StringLen(anp->label));
			flabel_len = get_max_feature_label (anp);
			if(flabel_len > (*max_label_size))
				*max_label_size = flabel_len; 
		}
		anp_list = anp_list->next;
	}
	buf[0] = '\0';
	sprintf(buf, "%ld", (long) (max_num+1));
	return StringLen(buf);
}



/*
*	do the converse of change blastx_master
*	revert it from tblastn to blastx
*/
static Boolean revert_blastx_alignment (ValNodePtr anp_list, AlignNodePtr master_anp)
{
	ValNodePtr c_list;
	AnnotInfoPtr annot_info;
	AlignNodePtr t_anp;
	AlignSegPtr asp;

	for(c_list = anp_list; c_list != NULL; c_list = c_list->next)
	{
		if(c_list->choice == OBJ_SEQANNOT)
		{
			annot_info = (AnnotInfoPtr) c_list->data.ptrvalue;
			if(annot_info != NULL && 
				get_alignment_type(annot_info) == ALIGN_PROT_TO_DNA)
				annot_info->blast_type = ALIGN_BLASTX;
		}
		else
		{
			t_anp = (AlignNodePtr) c_list->data.ptrvalue;
			if(t_anp == master_anp)
				t_anp->is_master = TRUE;
			else
				t_anp->is_master = FALSE;
			/*expand to the original interval */
			modify_gather_range (&(t_anp->extremes), TRUE);
			for(asp = t_anp->segs; asp != NULL; asp = asp->next)
			{
				if(asp->type != INS_SEG)
					modify_gather_range (&(asp->gr), TRUE);
			}
		}
	}
			
	return TRUE;
}

static void reverse_print(FILE *fp, CharPtr doc)
{
	Int4 i;
	CharPtr str;

	i = 0;
	for(str = doc; *str != '\n' && *str != '\0'; ++str)
		++i;

	if(*str == '\n')
	{
		++i;
		fprintf(fp, "%s", doc+i);
		*str = '\0';
		fprintf(fp, "%s\n", doc);
		*str = '\n';
	}
	else
		fprintf(fp, "%s", doc);
}

static ValNodePtr get_anp_list_for_aligntype(ValNodePtr anp_list, Uint1 align_type, 
	Int4 left, Int4 right)
{
	Uint1 c_type = 0;
	ValNodePtr list, prev;
	AlignNodePtr anp;
	AnnotInfoPtr annot_info;
	ValNodePtr  curr;
	Boolean first;
	Boolean extract;

	if(anp_list == NULL)
		return NULL;


	list = NULL;
	extract = (align_type == 0);
	prev = NULL;
	first = TRUE;
	annot_info = NULL;
	while(anp_list != NULL)
	{
		if(anp_list->choice == OBJ_SEQANNOT)
		{
			annot_info = (AnnotInfoPtr) anp_list->data.ptrvalue;
			c_type = get_alignment_type(annot_info);
			extract = (c_type == align_type);
			first = TRUE;
		}
		else if(extract)
		{
			anp = (AlignNodePtr) anp_list->data.ptrvalue;
			if(!anp->is_master)
			{
				/*check the first alignnode to see if it is legal*/
				if(first)
				{
					if(!PrintAlignForText(annot_info, anp))
						extract = FALSE;
					first = FALSE;
				}
				if(extract && 
				!(anp->extremes.left > right || anp->extremes.right < left))
				{	
					curr = ValNodeNew(NULL);
					curr->data.ptrvalue = anp_list;
					if(prev)
						prev->next = curr;
					else
						list = curr;
					prev = curr;
				}
			}
		}
		anp_list = anp_list->next;
	}

	return list;
}

static Boolean modify_separation_bar(CharPtr buf, Int4 size, Int1 frame)
{
	Char temp[21];
	Int4 len, start;
	Int4 i;

	sprintf(temp, "BLASTX: frame = %d", frame);
	len = StringLen(temp);
	start = (size - len)/2;
	if(start < 0)
		return FALSE;
	else
	{
		for(i = 0; i<len; ++i)
		{
			buf[start+i] = temp[i];
		}
		return TRUE;
	}
}


	
/*
*
*	for hardline old blast users who prefer to see the label as Sbjct/Query
*
*/
static void convert_label_to_query_subject (ValNodePtr anp_list)
{
	AlignNodePtr anp;
	Boolean first = TRUE;

	while(anp_list)
	{
		if(anp_list->choice != OBJ_SEQANNOT)
		{
			anp = (AlignNodePtr) anp_list->data.ptrvalue;
			if(anp->is_master)
			{
				MemFree(anp->label);
				anp->label = StringSave("Query:");
				first = FALSE;
			}
			else
			{
				MemFree(anp->label);
				if(first)
					anp->label = StringSave("Query:");
				else
					anp->label = StringSave("Sbjct:");
				first = FALSE;
			}
		}
		anp_list = anp_list->next;
	}
}

/*
*	request from Detlef: convert the sequence label to gi
*
*/
static void convert_label_to_gi(ValNodePtr anp_list)
{
	AlignNodePtr anp;
	Char temp[101];

	while(anp_list)
	{
		if(anp_list->choice != OBJ_SEQANNOT)
		{
			anp = (AlignNodePtr) anp_list->data.ptrvalue;
			if(anp && anp->sip && anp->sip->choice == SEQID_GI &&
				!anp->keep_label)
			{
				sprintf(temp, "%ld", (long) anp->sip->data.intvalue);
				if(anp->label != NULL)
					MemFree(anp->label);
				anp->label = StringSave(temp);
			}
		}
		anp_list = anp_list->next;
	}
}

/*
*
*	for the display of tblastx, only one Seq-annot is allowed at 
*	any given time
*
*/
static Boolean illegal_tblastx_anp (ValNodePtr anp_list, BoolPtr has_tblastx)
{
	AnnotInfoPtr info;
	Int2 info_num = 0;

	*has_tblastx = FALSE;

	while(anp_list)
	{
		if(anp_list->choice == OBJ_SEQANNOT)
		{
			info = (AnnotInfoPtr) anp_list->data.ptrvalue;
			if(info->blast_type == ALIGN_TBLASTX)
				*has_tblastx = TRUE;
			++info_num;
		}
		anp_list = anp_list->next;
	}

	if(*has_tblastx && info_num > 1)
		return TRUE;
	else
		return FALSE;
}

static Int4 expand_position(Int4 pos, Int4 exp_val, Boolean inverse)
{
	return inverse ? (pos/exp_val) : pos*exp_val;
}


static void modify_tblastx_value (ValNodePtr anp_list, Int4 val, Boolean inverse)
{
	AlignNodePtr anp;
	AlignSegPtr asp;
	Int4 left;

	while(anp_list)
	{
		if(anp_list->choice != OBJ_SEQANNOT)
		{
			anp = (AlignNodePtr) anp_list->data.ptrvalue;
			left = anp->extremes.left;
			/* anp->extremes.left = expand_position(anp->extremes.left, val, inverse); */
			anp->extremes.right = left + expand_position((anp->extremes.right -left), val, inverse);
			for(asp = anp->segs; asp != NULL; asp = asp->next)
			{
				if(asp->type == INS_SEG)
				{
					asp->ins_pos = left + expand_position((asp->ins_pos -left), val, inverse);
					asp->gr.right = expand_position(asp->gr.right, val, inverse);
				}
				else
				{
					asp->gr.left = left + expand_position(asp->gr.left - left , val, inverse);
					asp->gr.right = left + expand_position(asp->gr.right - left, val, inverse);
				}
			}
		}

		anp_list = anp_list->next;
	}
}

static
SeqIdPtr
get_seq_id(SeqAlignPtr sap, Int2 index)
{
	SeqIdPtr sip = NULL;
	if (sap->segtype == 1) {
		DenseDiagPtr ddp = (DenseDiagPtr)sap->segs;
		sip = ddp->id;
	}
	else if (sap->segtype == 2) {
		DenseSegPtr dsp = (DenseSegPtr)sap->segs;
		sip = dsp->ids;
	}
	for (; sip != NULL && --index >= 0; sip = sip->next);
	return sip;
}

/***********************************************************************
*
*	ShowAlignNodeText(anp_list, num_node, line_len, locus,
*	fp)
*	convert the alignment data in the list of AlignNode into text written
*	to a file
*	anp_list: a list (ValNodePtr) of AlignNode processed from Seq-aligns
*	num_node: the number of AlignNode to be processed currently. It can
*	be used in the cases where only the top num_node in the anp_list is 
*	going to be processed. This can be useful to make vertically cashed
*	buffer
*	line_len: the length of sequence char per line
*	locus: if TRUE, show the locus name
*	fp: the file Pointer
*	left: the leftmost position for display
*	right: the rightmost position for display
*	align_type:	the type of alignment. DNA-protein alignment?
*
*	return TRUE for success, FALSE for fail
*
************************************************************************/
NLM_EXTERN Boolean ShowAlignNodeText(ValNodePtr anp_list, Int2 num_node, Int4 line_len, FILE *fp, Int4 left, Int4 right, Uint4 option, Int4Ptr PNTR matrix, int (LIBCALLBACK *fmt_score_func)PROTO((AlignStatOptionPtr)))
{
    return ShowAlignNodeText2(anp_list, num_node, line_len, fp, left, right,
                              option, matrix, fmt_score_func, NULL, NULL, NULL);
}

NLM_EXTERN Boolean ShowAlignNodeText2(ValNodePtr anp_list, Int2 num_node, Int4 line_len, FILE *fp, Int4 left, Int4 right, Uint4 option, Int4Ptr PNTR matrix, int (LIBCALLBACK *fmt_score_func)PROTO((AlignStatOptionPtr)), CharPtr db_name, CharPtr blast_type, Int4Ptr PNTR posMatrix)
{
    CharPtr bar, sep_bar;
    CharPtr num_str;
    Int4 i, j;
    Int4 num;
    
    Int4 c_start, c_stop;
    CharPtr m_buf, cm_buf;		/*text for the master sequence*/
    BioseqPtr m_bsp;
    Int4 m_len;			/*length of the master sequence*/
    
    ValNodePtr list;	/*list of DrawText*/
    AlignNodePtr anp, master_anp;
    Int4Ptr p_stop;
    Boolean is_end, strip_semicolon=TRUE;
    CharPtr docbuf, master_docbuf;
    Uint1 all_frame[6];
    SeqFeatPtr fake_cds;
    Int1 frame;
    Uint1 align_type;
    ValNodePtr curr;
    ValNodePtr c_list, PNTR pc_list;
    ValNodePtr a_list;
    Boolean is_html;
    Int4Ptr PNTR t_matrix;
    Boolean compress;
    Int4 last_pos;
    Boolean load_last_pos;
    ValNodePtr id_list;
    BioseqPtr bsp;
    SeqAlignPtr align;
    ScorePtr sp;
    AlignStatOption aso;
    AlignSum as;
    AlignSumPtr asp;
    Boolean reverse_display;
    Boolean has_data;
    Boolean show_score;
    Int4 max_label_size;
    Int4 max_num_size;
    Int4 empty_space;
    Boolean show_strand;
    Boolean has_tblastx;
    SeqIdPtr already_linked=NULL;
    
    
    if(anp_list == NULL)
        return FALSE;
    
    /*for tblastx, only one Seq-annot at time*/
    if(illegal_tblastx_anp (anp_list, &has_tblastx))
        return FALSE;
    
    /*for alignment that is not a same-molecule, needs to have a master*/
    master_anp = get_master_align_node(anp_list);
    if(master_anp == NULL) {
        Message(MSG_ERROR, "Fail to the master AlignNode");
        return FALSE;
    }
    is_html = (Boolean)(option & TXALIGN_HTML);
    load_last_pos = (Boolean)(option & TXALIGN_END_NUM);
    
    
    /* for hard line old blast user only!!!!! */
    if(option & TXALIGN_SHOW_QS)
        convert_label_to_query_subject (anp_list);
    else if(option & TXALIGN_SHOW_GI)
        convert_label_to_gi(anp_list);
    
    compress = (Boolean)(option & TXALIGN_COMPRESS);
    if(compress)
        max_num_size = get_max_coordinates_len (anp_list, &max_label_size);
    else {
        max_num_size = POS_SPACE;
        max_label_size = B_SPACE;
    }

    /* for display of the traditional regular blastX output */
    reverse_display = FALSE;
    if(option & TXALIGN_BLASTX_SPECIAL) {
        reverse_display = reverse_blastx_order (anp_list);
        if(reverse_display)
            change_blastx_master(anp_list, &master_anp);
    } 
    
    if(has_tblastx)
        modify_tblastx_value (anp_list, 3, TRUE);
    if(left == -1)
        left = master_anp->extremes.left;
    if(right == -1)
        right = master_anp->extremes.right;
    if(left > master_anp->extremes.right || right < master_anp->extremes.left) {
        if(reverse_display)
            revert_blastx_alignment (anp_list, master_anp);
        return FALSE;
    }
    
    /*check for the molecule type of not-normal DNA-protein alignment*/
    fake_cds = NULL;
    frame = 0;
    m_bsp = BioseqLockById(master_anp->sip);
    if(m_bsp == NULL) {
        if(reverse_display)
            revert_blastx_alignment (anp_list, master_anp);
        return FALSE;
    }
    
    if(m_bsp->mol!= Seq_mol_aa) { /*a nucleotide sequence*/
        m_len = m_bsp->length;
        fake_cds = make_fake_cds(m_bsp, 0, m_bsp->length-1, Seq_strand_plus);
        load_master_translate_frame(anp_list, m_len, m_bsp);
    }
    
    
    ObjMgrSetHold();
    left = MAX(left, master_anp->extremes.left);
    right = MIN(right, master_anp->extremes.right);
    
    
    for(curr=anp_list, i=0; curr!=NULL; curr= curr->next) {	/*initiate the position*/
	
        if(curr->choice != OBJ_SEQANNOT) {
            anp = (AlignNodePtr) curr->data.ptrvalue;
            anp->align_num = i;
            ++i;
        }
    }
    p_stop = (Int4Ptr) MemNew((size_t)(i) * sizeof(Int4));
    for(j=0; j<i; ++j)
        p_stop[j] = -1;
    if(compress) {
        empty_space = max_num_size + 1 + max_label_size + 1;
        if(option & TXALIGN_SHOW_STRAND) {
            empty_space += STRAND_SPACE;
            show_strand = TRUE;
        } else
            show_strand = FALSE;
    } else {
        empty_space = get_num_empty_space(compress);
        show_strand = TRUE;
    }
    
    make_scale_bar_str(&bar, &num_str, compress? empty_space : empty_space+1, line_len);
    num = line_len + empty_space;
    sep_bar = (CharPtr) MemGet((size_t)(num+1)*sizeof(Char), MGET_ERRPOST);
    MemSet((Pointer)sep_bar, '-',(size_t)num* sizeof(Char));
    sep_bar[num] = '\0';
    
    if(is_html)
        fprintf(fp, "<PRE>\n");
    c_start = left;
    
    /*format the summary for the score */
    if(fmt_score_func != NULL) {
        aso.indent_len = (Int2)empty_space;
        aso.line_len = (Int2)(line_len + empty_space);
        aso.html_hot_link_relative = FALSE;
        if (option & TXALIGN_NO_ENTREZ)
            aso.no_entrez = TRUE;
        else
            aso.no_entrez = FALSE;
        
        if (option & TXALIGN_NO_DUMPGNL)
            aso.no_dumpgnl = TRUE;
        else
            aso.no_dumpgnl = FALSE;
        
        if (option & TXALIGN_HTML) {
            aso.html_hot_link = TRUE;
            if (option & TXALIGN_HTML_RELATIVE)
                aso.html_hot_link_relative = TRUE;
        } else {
            aso.html_hot_link = FALSE;
        }
        if (option & TXALIGN_SHOW_GI)
            aso.show_gi = TRUE;
        else
            aso.show_gi = FALSE;
        aso.fp = fp;
        aso.buf = NULL;
        id_list = NULL;
        aso.segs = NULL;
        if (blast_type == NULL) {
            aso.blast_type = StringSave("UNFIN_GEN");
        } else {
            aso.blast_type = StringSave(blast_type);
            StringUpper(aso.blast_type);
        }
        for(curr = anp_list; curr != NULL; curr = curr->next) {
            if(curr->choice != OBJ_SEQANNOT) {
                anp = (AlignNodePtr) curr->data.ptrvalue;
                show_score = FALSE;
                if((reverse_display && anp == master_anp) || (!reverse_display && anp != master_anp)) {
                    if(!check_bsp_id(&id_list, anp->sip))
                        show_score = TRUE;
                }
                if(show_score) {
                    /*the first time it sees the Bioseq*/
                    bsp = BioseqLockById(anp->sip);
                    align = NULL;
                    GatherItem(anp->entityID, anp->itemID, (Uint2)(curr->choice), (Pointer)(&align), find_align_proc);
                    if(align != NULL) {
                        if(align->segtype == 1 || align->segtype == 2 || align->segtype == 3 || align->segtype == 5) {
                            as.matrix = matrix;
                            as.posMatrix = posMatrix;
                            as.master_sip = master_anp->sip;
                            as.target_sip = anp->sip;
                            as.is_aa = (m_bsp->mol == Seq_mol_aa);
                            as.ooframe = FALSE; /* Not supported */
                            asp = &as;
                        } else
                            asp = NULL;
                        sp = find_score_in_align(align, anp->chain, asp);
                        if(sp != NULL) {
                            aso.follower = anp->follower;
                            aso.bsp = bsp;
                            aso.sp = sp;
                            aso.db_name = db_name;
                            if(asp != NULL) {
                                aso.gaps = asp->gaps;
                                aso.positive = asp->positive;
                                aso.identical = asp->identical;
                                aso.align_len = asp->totlen;

                                if (asp->m_frame_set) {
                                        aso.m_frame = asp->m_frame;
                                } else {
                                    aso.m_frame = 255;
                                }

                                if (asp->t_frame_set) {
                                    aso.t_frame = asp->t_frame;
                                } else {
                                    aso.t_frame = 255;
                                }

                                aso.m_strand = asp->m_strand;
                                aso.t_strand = asp->t_strand;
                            } else {
                                aso.align_len = 0;
                            }

                            {{
                                SeqAlignPtr sap;
                                size_t size = 0;
                                size_t used = 0;
                                
                                aso.segs = NULL;
                                for (sap = align; sap != NULL; sap = sap->next) {
                                    if (SeqIdMatch(TxGetSubjectIdFromSeqAlign(align), TxGetSubjectIdFromSeqAlign(sap))) {
                                        if (aso.segs != NULL) {
                                            StringAppend(&aso.segs, &size, ",", &used);
                                        }
                                        SeqAlignSegsStr(sap, 1, &aso.segs, &size, &used);
                                    } else
                                        break;
                                }
                                if (aso.segs == NULL) {
                                    /**
                                     * Something is really wrong if we're here
                                     */
                                    aso.segs = StringSave("");
                                }
                            }}
                            fmt_score_func(&aso);
                        }
                        aso.segs = (CharPtr) MemFree(aso.segs);
                    }
                    if(bsp != NULL)
                        BioseqUnlock(bsp);
                }
            }
        }
        aso.blast_type = (CharPtr) MemFree(aso.blast_type);
        ValNodeFree(id_list);
    }
    
    pc_list = (ValNodePtr *) MemNew((size_t)(ALIGN_MAX_TYPE +1) * sizeof(ValNodePtr));
    has_data = FALSE;

    for(i = 0; i<=ALIGN_MAX_TYPE; ++i) {
        pc_list[i]= get_anp_list_for_aligntype(anp_list, (Uint1)i, left, right);
        if(pc_list[i] != NULL)
            has_data = TRUE;
    }
    
    if(option & TXALIGN_SHOW_QS) {
        is_html = FALSE;
        strip_semicolon = FALSE;
    }
    master_docbuf = NULL;
    if(has_data) {
	while(c_start <= right)	{ /*process line by line*/
            
            c_stop = MIN(right, (c_start+line_len -1));
            m_buf = NULL;
            is_end = FALSE;
            docbuf = NULL;
            
            /*process the master sequence*/
            if(has_tblastx)
                list = ProcessTextAlignNode(master_anp, c_start, c_stop, &(p_stop[master_anp->align_num]), NULL, line_len, -1, option, matrix);
            else
                list = ProcessTextAlignNode(master_anp, c_start, c_stop, &(p_stop[master_anp->align_num]), NULL, line_len, 0, option, NULL);
            if(list != NULL) {
                if(option & TXALIGN_SHOW_RULER) {
                    fprintf(fp, "%s\n", num_str);	/*show scale*/
                    fprintf(fp, "%s\n", bar);
                }
                last_pos = load_last_pos ? p_stop[master_anp->align_num] : -1;
                docbuf = DrawTextToBuffer(list, &m_buf, is_html, max_label_size, max_num_size, compress, matrix, last_pos, line_len, show_strand, strip_semicolon, &already_linked);
                if(docbuf !=NULL) {
                    if(reverse_display)
                        master_docbuf = docbuf;
                    else {
                        fprintf(fp, "%s", docbuf);
                        MemFree(docbuf);
                    }
                }
                list = FreeTextAlignList(list);
            }
            
            for(align_type = 0; align_type <= ALIGN_MAX_TYPE; ++align_type) {
                c_list = pc_list[align_type];
                if(c_list != NULL) {
                    if(align_type == ALIGN_DNA_TO_PROT) {
                        /*process the hit protein sequence*/
                        if(get_current_master_frame(c_list, c_start, c_stop, all_frame)) {
                            list = NULL;
                            for(j = 0; j<6; ++j) {
                                frame = all_frame[j];
                                if(frame > 0) {
                                    /*translate the master sequence in the specified frame*/
                                    cm_buf = translate_faked_cds(fake_cds, frame, c_start, c_stop, m_len, master_anp); 
                                    list = load_fake_protein_buf(cm_buf, frame, master_anp);
                                    docbuf = DrawTextToBuffer(list, NULL, is_html, max_label_size, max_num_size, compress, matrix, -1, line_len, show_strand, strip_semicolon, &already_linked);
                                    if(docbuf != NULL) {
                                        modify_separation_bar(sep_bar, num, frame);
                                        fprintf(fp, "%s\n", sep_bar);
                                        fprintf(fp, "%s", docbuf);
                                        MemFree(docbuf);
                                    }
                                    FreeTextAlignList(list);
                                    
                                    for(curr = c_list; curr != NULL; curr = curr->next) {
                                        a_list = (ValNodePtr) curr->data.ptrvalue;
                                        anp = (AlignNodePtr) a_list->data.ptrvalue;
                                        if(anp != master_anp) {
                                            list = ProcessTextAlignNode(anp, c_start, c_stop, &(p_stop[anp->align_num]), cm_buf, line_len, frame, option, matrix);
                                            if(list != NULL) {
                                                docbuf = DrawTextToBuffer(list, NULL, is_html, max_label_size, max_num_size, compress, matrix, -1, line_len, show_strand, strip_semicolon, &already_linked);
                                                if(docbuf != NULL) {
                                                    fprintf(fp, "%s", docbuf);
                                                    MemFree(docbuf);
                                                }
                                                FreeTextAlignList(list);
                                            }
                                        }
                                    }
                                    MemFree(cm_buf);
                                } /*end of frame > 0 */
                            }
                        }
                    } else {
                        if(align_type == ALIGN_PROT_TO_DNA || align_type == ALIGN_TDNA_TO_TDNA)
                            frame = -1;
                        else
                            frame = 0;
                        if(frame == 0 && m_bsp->mol != Seq_mol_aa)
                            t_matrix = NULL;
                        else
                            t_matrix = matrix;
                        is_end = FALSE;
                        for(curr = c_list; curr !=NULL; curr = curr->next) {
                            a_list = (ValNodePtr) curr->data.ptrvalue;
                            anp = (AlignNodePtr) a_list->data.ptrvalue;
                            if(anp != master_anp) {
                                /*generate the DrawText buffer*/
                                if(align_type == ALIGN_NORMAL && m_bsp->mol != Seq_mol_aa)
                                    /*DNA alignment */
                                    list = ProcessTextAlignNode(anp, c_start, c_stop, &(p_stop[anp->align_num]), m_buf, line_len, frame, option, NULL);
                                else
                                    list = ProcessTextAlignNode2(anp, c_start, c_stop, &(p_stop[anp->align_num]), m_buf, line_len, frame, option, t_matrix, posMatrix, master_anp->seqpos);
                                
                                last_pos = load_last_pos ? p_stop[anp->align_num] : -1;
                                if(list != NULL) {
                                    /*DrawTextList(list, fp);*/
                                    docbuf = DrawTextToBuffer(list, NULL, is_html, max_label_size, max_num_size, compress, t_matrix, last_pos, line_len, show_strand, strip_semicolon, &already_linked);
                                    if(docbuf !=NULL) {
                                        if(reverse_display) {
                                            reverse_print(fp, docbuf);
                                            fprintf(fp, "%s", master_docbuf);
                                            MemFree(master_docbuf);
                                        } else
                                            fprintf(fp, "%s", docbuf);
                                        MemFree(docbuf);
                                    }
                                    list = FreeTextAlignList(list);
                                }
                            }
                        }
                    }	/*end of else*/
                }
            }
            
            if(m_buf != NULL)
                MemFree(m_buf);
            if(c_stop < right)
                fprintf(fp, "\n"); 
            c_start = c_stop+1;
	}
    }
    for(i = 0; i<=ALIGN_MAX_TYPE; ++i) {
        if(pc_list[i] != NULL)
            ValNodeFree(pc_list[i]);
    }
    MemFree(pc_list);
    
    if(option & TXALIGN_HTML)
        fprintf(fp, "</PRE>\n");
    
    already_linked = ValNodeFree(already_linked);
    
    if(fake_cds != NULL)
        SeqFeatFree(fake_cds);
    BioseqUnlock(m_bsp);
    if(has_tblastx)
        modify_tblastx_value (anp_list, 3, FALSE);
    MemFree(num_str);
    MemFree(sep_bar);
    MemFree(bar);
    MemFree(p_stop);
    ObjMgrClearHold();
    if(reverse_display)
        revert_blastx_alignment (anp_list, master_anp);
    return has_data;
}
				
	


/*######################################################################
#
#	functions related to ProcessTextAlignNode
#
#######################################################################*/


/******************************************************************************
*
*	load_text(bsp, pos1, pos2, l_seq, l_pos, mbuf, maxlen)
*	load the sequence into text
*	bsp: the Bioseq
*	pos1: the first position on the sequence. 
*	pos2: the second position on the sequence. 
*		if(pos1 and pos2 are negative val, indicate the region in on the*	minus strand
*	l_seq: the buffer for loading the sequence
*	l_pos: the current position in l_seq. Will be updated after the sequence
*		is loaded
*	mbuf: buffer from the master sequence. For checking mismatches and positive scores
*	maxlen: the maximum length per line. Used to determine the special 
*	format used for long insertions
*	spacing is the space between the two adjacent residues
*	mismatch: if TRUE, show the identical residue with 
*
*****************************************************************************/

static Boolean load_text(BioseqPtr bsp, Int4 pos1, Int4 pos2, CharPtr l_seq, Int4Ptr l_pos, CharPtr mbuf, Int2 maxlen, Int2 spacing, Boolean translate, Boolean mismatch, Int2Ptr matrix_val, Int4Ptr PNTR matrix, Uint1 strand, Int4Ptr PNTR posMatrix, Int4 q_start)
{
    SeqPortPtr spp = NULL;
    ByteStorePtr b_store = NULL;
    Uint1 code;
    Int4 start, stop;
    Uint1 m_res, t_res, stdaa_res;
    Int2 i;
    Int2 val;
    Int4 length, s_len;
    Int2 c_pos;
    Char temp[100];
    Boolean protein;
    Boolean overflow;
    Boolean reverse;
    Boolean is_real;
    SeqFeatPtr fake_cds;
    Boolean check_neg;	/*if aa is negative, load it as lower case char*/
    SeqMapTablePtr smtp;
    
    if(*l_pos >= maxlen )
        return FALSE;

    /* posMatrix uses NCBIstdaa encoding */
    
    if(posMatrix != NULL) { 
        if((smtp = SeqMapTableFindObj(Seq_code_ncbistdaa, 
                                      Seq_code_ncbieaa)) == NULL)
            return FALSE;
    }

    protein = (bsp->mol == Seq_mol_aa);
    reverse = FALSE;
    if(protein)
        code = Seq_code_ncbieaa;
    else
        code = Seq_code_iupacna;
    check_neg = (matrix_val == NULL && matrix != NULL);
    if(strand == Seq_strand_minus) {	/*on the minus strand*/

        start = -pos2;
        stop = -pos1;
        
        if(protein) {
            strand = Seq_strand_plus;
            reverse = TRUE;
        }
        
    } else {
        start = pos1;
        stop = pos2;
    }
    if(translate) {
        fake_cds = make_fake_cds(bsp, start, stop, strand);
        b_store = ProteinFromCdRegionEx(fake_cds, TRUE, FALSE);
        SeqFeatFree(fake_cds);
        if(b_store == NULL)
            return FALSE;
        length = (stop - start +1)/3;
        BSSeek(b_store, 0, SEEK_SET);
    } else {
        spp = SeqPortNew(bsp, start, stop, strand, code);
        length = stop - start +1;
    }
    c_pos = (Int2)(*l_pos);
    overflow = (c_pos >= maxlen);
    if(maxlen>0 && (length > maxlen)) {	/*large insertions*/
	
        for(i =0; i<5 && !overflow; ++i) {
            if(translate)
                l_seq[c_pos++] = (Uint1)BSGetByte(b_store);
            else {
                if(reverse)
                    SeqPortSeek(spp, length-1 -i, SEEK_SET);
                l_seq[c_pos++] = SeqPortGetResidue(spp);
            }
            overflow = (c_pos >= maxlen);
        }
        for(i =0; i<3 && !overflow; ++i) {
            l_seq[c_pos++] = '.';
            overflow = (c_pos >= maxlen);
        }
        if(!overflow) {
            if(translate)
                BSSeek(b_store, length-1, SEEK_SET);
            else if(!reverse)
                SeqPortSeek(spp, length-5, SEEK_SET);
            for(i =0; i<5 && !overflow; ++i) {
                if(translate)
                    l_seq[c_pos++] = (Uint1)BSGetByte(b_store);
                else {
                    if(reverse)
                        SeqPortSeek(spp, 4-i, SEEK_SET);
                    l_seq[c_pos++] = SeqPortGetResidue(spp);
                }
                overflow = (c_pos >= maxlen);
            }
        }
        if(overflow)
            l_seq[maxlen-1] = '\0';
        else
            l_seq[c_pos] = '\0';
        sprintf(temp, "(length=%ld)", (long) length);
        s_len = StringLen(temp);
        StringCat(l_seq, temp);
        *l_pos = c_pos+s_len;
    } else {
        if(translate) {
            while((val = BSGetByte(b_store)) != EOF) {
                t_res = (Uint1)val;
                l_seq[c_pos]= t_res;
                if(mbuf != NULL) {
                    m_res = mbuf[c_pos];
                    if(matrix_val && matrix)
                        matrix_val[c_pos] = (Int2)matrix[m_res][t_res];
                    if(mismatch && t_res == m_res)
                        l_seq[c_pos] = '.';
                    else if(check_neg && matrix[t_res][m_res] < 0)
                        l_seq[c_pos] = TO_LOWER(t_res);
                }
                c_pos += spacing;
                if(c_pos >= maxlen) {
                    c_pos = maxlen;
                    break;
                }
            }
        } else {
            if(reverse)
                SeqPortSeek(spp, length-1, SEEK_SET);
            s_len = 0;
            while((t_res = SeqPortGetResidue(spp)) != SEQPORT_EOF) {
                is_real = IS_ALPHA(t_res);
                if(is_real || t_res == '*' || t_res == '-') {
                    if(is_real && !protein)
                        t_res = TO_LOWER(t_res);
                    l_seq[c_pos] = t_res;
                    if(mbuf != NULL) {
                        m_res = mbuf[c_pos];
                        if(matrix_val) {
                            if(matrix) {
                                if(posMatrix != NULL) {
                                    if(t_res == m_res) /* complete match */
                                        matrix_val[c_pos] = INT2_MAX;
                                    else {
                                        stdaa_res = SeqMapTableConvert(smtp, t_res);
                                        matrix_val[c_pos] = (Int2)posMatrix[c_pos + q_start][stdaa_res];

                                        /* 
                                     if(posMatrix[c_pos + q_start][t_res] == 
                                     matrix[t_res][t_res]) {
                                     printf("Got it!");
                                     } */

                                    }
                                } else {
                                    matrix_val[c_pos] = (Int2)matrix[m_res][t_res];
                                }
                                
                            } else if(t_res == m_res)
                                matrix_val[c_pos] = '|';
                        }
                        
                        if(mismatch && t_res == m_res)
                            l_seq[c_pos] = '.';
                        else if(posMatrix != NULL) {
                            stdaa_res = SeqMapTableConvert(smtp, m_res);
                            if(check_neg && posMatrix[c_pos + q_start][stdaa_res] < 0)
                                l_seq[c_pos] = TO_LOWER(t_res);
                        } else { /*regular BLOSSUM62*/
                            if(check_neg && matrix[t_res][m_res] < 0)
                                l_seq[c_pos] = TO_LOWER(t_res);
                        }
                    }
                    c_pos += spacing;
                    if(c_pos >= maxlen) {
                        c_pos = maxlen;
                        break;
                    }
                    ++s_len;
                }
                if(reverse) {
                    if(s_len == length)
                        break;
                    else
                        SeqPortSeek(spp, length -1 - s_len, SEEK_SET);
                }
            }
        }
        *l_pos = c_pos;
    }
    
    if(translate)
        BSFree(b_store);
    else
        SeqPortFree(spp);
    return TRUE;
}

/*##########################################################################
#
#	functions related to add the features to the alignment
#
###########################################################################*/


typedef struct protbuf{	/*for loading the translation of a CDs*/
	CharPtr buf;	/*load the protein sequence*/
	Int4 start;	/*start position in CDs*/
	Int4 stop;	/*stop position in CDs*/
	Int4 pos;	/*position for the feature*/
	Boolean load_codon;	/*load the codon data for aa sequence*/
	ValNodePtr cvp_list;	/*list for loading the codon of an aa*/
}ProtBuf, PNTR ProtBufPtr;

		
			
/************************************************************************
*
*	check the protein sequence from CDs feature into the buffer
*
*************************************************************************/
static Boolean load_prot_seq(GatherContextPtr gcp)
{
	SeqFeatPtr sfp;
	ProtBufPtr pbp;
	SeqLocPtr loc;

	if(gcp->thistype != OBJ_SEQFEAT)
		return FALSE;
	sfp = (SeqFeatPtr)(gcp->thisitem);
	if(sfp->data.choice !=3)
		return FALSE;

	pbp = (ProtBufPtr)(gcp->userdata);
	if(pbp->load_codon)	/*looking for codon in aa sequence*/
	{
		pbp->cvp_list = aa_to_codon(sfp, pbp->start, pbp->stop);
		return (pbp->cvp_list !=NULL);
	}
		

	if(pbp->start <0)/*minus strand*/
		loc = SeqLocIntNew((-pbp->stop), (-pbp->start), Seq_strand_minus, SeqLocId(sfp->location));
	else
		loc = SeqLocIntNew(pbp->start, pbp->stop, Seq_strand_plus, SeqLocId(sfp->location));

	pbp->pos = print_protein_for_cds(sfp, pbp->buf, loc, TRUE);
	SeqLocFree(loc);
	return (pbp->pos != -1);
}



static Boolean buffer_for_feature(Int4 c_left, Int4 c_right, Int4 seq_start, Int4 seq_stop, ValNodePtr fnp_node, Boolean load_codon, ProtBufPtr pbp)
{
	FeatNodePtr fnp;
	Uint2 itemtype;
	CharPtr buf = NULL;
	Int2 i=0;
	Char symbol;
	ValNodePtr curr;
	IvalNodePtr inp;
	Int4 i_left, i_right;
	Int4 f_len;


	itemtype = (Uint2)(fnp_node->choice);

	if(itemtype!= OBJ_SEQFEAT)
		return FALSE;
	fnp = (FeatNodePtr) fnp_node->data.ptrvalue;
	f_len = seq_stop - seq_start +1; 
	if(load_codon)
		pbp->buf = NULL;
	else
		pbp->buf = (CharPtr) MemNew((size_t)(f_len +1)*sizeof(Char));
	pbp->start = seq_start;
	pbp->stop = seq_stop;
	pbp->pos = -1;
	pbp->load_codon= load_codon;
	pbp->cvp_list = NULL;

	buf = pbp->buf;
	if(buf !=NULL)
		MemSet((Pointer)buf,  '~', (size_t)(f_len) * sizeof(Char));
	switch(fnp->feattype)/*check symbol for different features*/
	{
		case FEATDEF_GENE:
			symbol = '+';
			break;
		case FEATDEF_mRNA:
			symbol = '^';
			break;
		case FEATDEF_CDS:
			symbol = '$';
			break;
		default:
			symbol = '*';
			break;
	}
	if(fnp->feattype ==FEATDEF_CDS)

		GatherItem(fnp->entityID, fnp->itemID, itemtype, (Pointer)(pbp), load_prot_seq);
	else
	{
		if(fnp->interval !=NULL)
		{
			for(curr = fnp->interval; curr !=NULL; curr = curr->next)
			{
				inp = (IvalNodePtr) curr->data.ptrvalue;
				i_left = inp->gr.left;
				i_right = inp->gr.right;
				if(!(i_left > c_right || i_right < c_left))
				{
					i_left = MAX(i_left, c_left);
					i_right = MIN(i_right, c_right);
					i_left -= c_left;
					i_right -=c_left;
					for(; i_left<=i_right; ++i_left)
						buf[i_left] = symbol;
				}
			}
		}
		else
		{
			i_left = fnp->extremes.left;
			i_right = fnp->extremes.right;
			if(!(i_left > c_right || i_right < c_left))
			{
				i_left = MAX(i_left, c_left);
				i_right = MIN(i_right, c_right);
				i_left -= c_left;
				i_right -=c_left;
				for(; i_left<=i_right; ++i_left)
					buf[i_left] = symbol;
			}
		}
			
	}
	if(buf!=NULL)
		buf[f_len]= '\0';
	if(pbp->pos == -1)
		pbp->pos = ABS(seq_start);

	if(pbp->buf != NULL || pbp->cvp_list !=NULL)
		return TRUE;
	else
		return FALSE;
}


	
static Boolean load_feature_data(ProtBufPtr pbp, FeatNodePtr fnp, Int4 pos, Int4 maxlen, ValNodePtr PNTR fbp_head)
{
	Boolean found;
	TextAlignBufPtr fbp;
	ValNodePtr curr, pcvp;
	CodonVectorPtr cvp;
	Boolean load_codon;
	CharPtr PNTR codon;
	Int2 i;
	Int4 f_len;
	Char label[100];
	CharPtr buf;
	Boolean locus = FALSE;
	
	if(pbp == NULL)
		return FALSE;
	if(pbp->buf == NULL && pbp->cvp_list == NULL)
		return FALSE;
	load_codon = (pbp->cvp_list !=NULL);
	f_len = pbp->stop - pbp->start +1;
	
	found = FALSE;
	for(curr = *fbp_head; curr !=NULL; curr = curr->next)
	{
	   fbp = (TextAlignBufPtr) curr->data.ptrvalue;
	   if(fbp->itemID == fnp->itemID)
	   {
	     if(load_codon)
	     {
		for(pcvp = pbp->cvp_list; pcvp!=NULL; pcvp= pcvp->next)
		{
		   cvp = (CodonVectorPtr) pcvp->data.ptrvalue;
		   if(cvp->exonCount == fbp->exonCount)
		   {
			codon = fbp->codon;
			for(i =0; i<3; ++i)
			{
			   if(pos > fbp->f_pos)
				make_empty(codon[i] + fbp->f_pos, (Int2)(pos - fbp->f_pos));
			   StringCat(codon[i], (cvp->buf[i]+cvp->aa_index));
			}
			cvp->exonCount = 0;
			fbp->f_pos = pos + f_len;
		   }
		}/*end of for*/
	     }
	     else
	     {
		if(fbp->pos == -1)
		   fbp->pos = pbp->pos+1;
		if(pos > fbp->f_pos)
		   make_empty(fbp->buf+fbp->f_pos, (Int2)(pos - fbp->f_pos));	
		StringCat(fbp->buf, pbp->buf);
		fbp->f_pos = pos + f_len;
		found = TRUE;
	     }
	   }
	}


	if(load_codon)
	{
	   for(pcvp = pbp->cvp_list; pcvp!=NULL; pcvp= pcvp->next)
	   {
		cvp = (CodonVectorPtr) pcvp->data.ptrvalue;
		if(cvp->exonCount !=0)
		{
		   fbp = (TextAlignBufPtr) MemNew(sizeof(TextAlignBuf));
		   fbp->seqEntityID = fnp->entityID;
		   fbp->pos = cvp->dna_pos +1;
		   fbp->strand = cvp->strand;
		   seqid_name(cvp->sip, label, locus, FALSE);
		   fbp->label = StringSave(label);
		   fbp->buf = NULL;
		   for(i =0; i<3; ++i)
		   {
		   	buf = (CharPtr) MemNew((size_t)(maxlen+1+1+1) * sizeof(Char));
			/*1 for partial start, 1 for partial stop*/
			if(pos > 0)
				make_empty(buf, (Int2)pos);
			StringCat(buf, cvp->buf[i]+cvp->aa_index);
			fbp->codon[i] = buf;
		   }
		   fbp->frame = cvp->frame;
		   fbp->f_pos = pos+f_len;
		   fbp->exonCount = cvp->exonCount;
		   fbp->itemID = fnp->itemID;
		   fbp->itemID = fnp->itemID;
		   fbp->feattype = fnp->feattype;
		   fbp->subtype = fnp->subtype;
		   fbp->entityID = fnp->entityID;
		   fbp->extra_space = (cvp->aa_index == 0);
		   ValNodeAddPointer(fbp_head, 0, fbp);
		
		}
	     }
	}
	else
	{
	   	if(!found)	/*create a new node*/
		{
		   fbp = (TextAlignBufPtr) MemNew(sizeof(TextAlignBuf));
		   buf = (CharPtr) MemNew((size_t)(maxlen+1) * sizeof(Char));
		   if(pos > 0)
			make_empty(buf, (Int2)pos);
		   StringCat(buf, pbp->buf);
		   fbp->seqEntityID = fnp->entityID;
		   fbp->f_pos = pos + f_len;
		   fbp->itemID = fnp->itemID;
		   fbp->buf = buf;
		   fbp->pos = pbp->pos+1;
		   if(fnp->label !=NULL)
			fbp->label = StringSave(fnp->label);
		   fbp->strand = fnp->extremes.strand;
		   fbp->itemID = fnp->itemID;
		   fbp->feattype = fnp->feattype;
		   fbp->subtype = fnp->subtype;
		   fbp->entityID = fnp->entityID;
		   fbp->exonCount = 0;
		   ValNodeAddPointer(fbp_head, 0, fbp);
		}
	}
	if(pbp->buf)
		MemFree(pbp->buf);
	if(pbp->cvp_list)
		free_cvp_list(pbp->cvp_list);
	return TRUE;
}


	
/**************************************************************************
*
*	collect_feature_buf(fnp_list, g_left, g_right, seq_start, l_pos, 
*	fbp_head, max_len)
*	collect the features to be shown together with the alignment
*	fnp_list: a list of FeatNode associated with the current segment
*	g_left: the left position 
*
***************************************************************************/
static ValNodePtr collect_feature_buf(ValNodePtr fnp_list, Int4 g_left, Int4 g_right, Int4 seq_start, Int4 l_pos, ValNodePtr fbp_head, Int4 maxlen, Boolean is_aa)
{
	ProtBuf pb;
	FeatNodePtr fnp;
	Int4 c_left, c_right;
	Int4 pos;
	Int4 fseq_start, fseq_stop;	/*map sequence start stop to the feature*/
	Int4 f_len;		/*length of the feature*/
	Boolean load_codon;

	if(fnp_list == NULL)
		return fbp_head;

	
	while(fnp_list)
	{
	   fnp = (FeatNodePtr) fnp_list->data.ptrvalue;
	   c_left = fnp->extremes.left;
	   c_right = fnp->extremes.right;
	   load_codon = (is_aa && fnp->feattype == FEATDEF_CDS);
	   if(!(c_left > g_right || c_right < g_left))
	   {
		if(c_left > g_left)	/*map the seq pos from the graphic pos*/
			fseq_start = seq_start + (c_left-g_left);
		else
			fseq_start = seq_start;
		c_left = MAX(c_left, g_left);
		c_right = MIN(c_right, g_right);
		f_len = c_right - c_left+1;
		fseq_stop = fseq_start+f_len-1;

		if(c_left > g_left)
		   pos = l_pos + (c_left - g_left);
		else
		   pos = l_pos;

		if(buffer_for_feature(c_left, c_right, fseq_start, fseq_stop, fnp_list, load_codon, &pb))
		
			load_feature_data(&pb, fnp, pos, maxlen, &fbp_head);
	   }
	   fnp_list = fnp_list->next;
	}

	return fbp_head;
}
	
static Int4 map_position_by_spacing(Int4 distance, Int4 spacing, Boolean is_head)
{
	Int4 pos, left_over;

	if(spacing == 1)
		return distance;

	pos = distance/spacing;
	left_over = distance%spacing;

	if(left_over == 0 && !is_head)
		pos = MAX(pos-1, 0);
	else if(left_over == 2 && is_head)
		++pos;
	return pos;
}
			 
	

/***********************************************************************
*
*	ProcessTextAlignNode(anp, m_left, m_right, p_stop, m_buf, locus)
*	Process the AlignNode to make a list of text buffer on the 
*	current region
*	anp: AlignNodePtr
*	m_left, m_right: the region on the alignment. Mapped in response
*	to anp->extremes.left, and anp->extremes.right
*	p_stop: the stop position of the previous segment. Used to label
*	the position of a line composed entirely of gaps
*	m_buf: buffer for the master sequence. Used to compare mismatches
*	locus: if TRUE, show the locus name of the alignment
*
*	frame: frame >0, those are the hits from blastx. So, the 
*	protein need to be displayed to the proper frame
*	frame 1-3: match to the plus strand of the master
*	frame 4-6: match to the minus strand of the master
*	frame 0:   no tranlsation, no frame match to the master
*	frame -1:  translate the DNA sequence
*	option:    option for display the alignments
*   matrix:	   the protein alignment matrix
*
*
************************************************************************/

NLM_EXTERN ValNodePtr ProcessTextAlignNode2(AlignNodePtr anp, Int4 m_left, Int4 m_right, Int4Ptr p_stop, CharPtr m_buf, Int4 line_len, Int1 m_frame, Uint4 option, Int4Ptr PNTR matrix, Int4Ptr PNTR posMatrix, Int4 q_start)
{
    Int4 maxlen;
    Int4 g_left, g_right;
    Int4 len;		/*length of the segment*/
    CharPtr l_seq;	/*the buffer for the sequence*/
    Int2Ptr matrix_val;	/*value of each residue in alignment matrix*/
    Int4 l_pos;		/*the start position on the line*/	
    Int4 offset;
    BioseqPtr bsp;
    SeqEntryPtr sep;
    
    AlignSegPtr asp;
    Int4 seq_offset, off_len;
    Int4 seq_start, seq_stop;
    Int4 s_start, s_stop;	/*for marking the position on one line*/
    CharPtr str;
    
    ValNodePtr head = NULL, ins_node;
    ValNodePtr fbuf_list = NULL;
    TextAlignBufPtr tdp;
    Boolean is_aa;
    Int4 spacing;
    Boolean translate;
    Int4 seq_expand;
    Boolean show_mismatch;
    Boolean set_matrix;
    Uint1 strand;
    
    
    if(m_frame > 6 || m_frame < -1)	/*check the m_frame. -1 for translate the hits*/
        return NULL;
    
    
    g_left = anp->extremes.left;
    g_right = anp->extremes.right;
    if(m_left > g_right || m_right < g_left)/*no overlap*/ {
        if(m_frame > 0) {
            if(anp->m_frame != m_frame)
                return NULL;
            if(m_buf == NULL)
                return NULL;
        }
        if(option & TXALIGN_BLUNT_END) {
            maxlen = m_right - m_left +1;
            l_seq = (CharPtr) MemGet((size_t)(maxlen+1)*sizeof(Char), MGET_ERRPOST);
            MemSet((Pointer)l_seq, '-',(size_t)(maxlen) * sizeof(Char));
            l_seq[maxlen] = '\0';
            tdp = (TextAlignBufPtr) MemNew(sizeof(TextAlignBuf));
            tdp->pos = *p_stop;
            tdp->strand = anp->extremes.strand;
            tdp->label = StringSave(anp->label);
            tdp->buf = l_seq;
            tdp->matrix_val = NULL;
            tdp->itemID = anp->itemID;
            tdp->feattype = 0;
            tdp->subtype = 0;
            tdp->entityID = anp->entityID;
            tdp->seqEntityID = anp->seq_entityID;
            tdp->bsp_itemID = anp->bsp_itemID;
            ValNodeAddPointer(&head, 0, tdp);
            return head;
        }
        else
            return NULL;
    }
    
    strand = Seq_strand_plus;
    if(anp->seqpos < 0)
        strand = Seq_strand_minus;
    else if(anp->seqpos == 0 && anp->extremes.strand == Seq_strand_minus)
        strand = Seq_strand_minus;
    
    l_pos = 0;
    spacing = 1;
    offset = 0;
    if(m_frame > 0) {
        if(anp->m_frame != m_frame)
            return NULL;
        if(m_buf == NULL)
            return NULL;
        /*add the empty space to reflect the reading frame*/
        for(str = m_buf; *str != '\n' && *str != '\0'; ++str) {
            if(IS_WHITESP(*str))
                ++offset;
            else
                break;
        }
        spacing = 3;
    }
    if(m_left < g_left) {
        l_pos += (g_left - m_left);
        if(m_frame > 0)
            ++l_pos;
    } else
        l_pos += offset;
    
    bsp = BioseqLockById(anp->sip);
    if(bsp == NULL)
        return NULL;
    is_aa = (bsp->mol == Seq_mol_aa);
    if((m_frame > 0 && !is_aa) || (m_frame == -1 && is_aa)) {
        BioseqUnlock(bsp);
        return NULL;
    }
    if(anp->seq_entityID == 0) {
        sep = SeqEntryFind(bsp->id);
        anp->seq_entityID = SeqMgrGetEntityIDForSeqEntry(sep);
    }
    if(anp->bsp_itemID == 0)
        anp->bsp_itemID = get_bioseq_itemID(bsp, anp->seq_entityID);
    
    if(m_frame == -1) {
        translate = TRUE;
        seq_expand = 3;
    } else {
        translate = FALSE;
        seq_expand = 1;
    }
    
    maxlen = m_right - m_left +1;
    l_seq = (CharPtr) MemGet((size_t)(maxlen+1)*sizeof(Char), 
                             MGET_ERRPOST);
    if(option & TXALIGN_BLUNT_END)
        MemSet((Pointer)l_seq, '-',(size_t)maxlen * sizeof(Char));
    else
        MemSet((Pointer)l_seq, ' ',(size_t)maxlen * sizeof(Char));
    l_seq[maxlen] = '\0';
    
    
    set_matrix = FALSE;
    if(m_frame == 0 && bsp->mol != Seq_mol_aa) { /*DNA-DNA alignment*/
        if(option & TXALIGN_MATRIX_VAL)
            set_matrix = TRUE;
    } else {
        if(matrix != NULL && (option & TXALIGN_MATRIX_VAL))
            set_matrix = TRUE;
    }
    if(set_matrix) {
        matrix_val = (Int2Ptr) MemGet((size_t)(maxlen+1)*sizeof(Int2), MGET_ERRPOST);
        MemSet((Pointer)matrix_val, 0,(size_t)maxlen * sizeof(Int2));
    } else
        matrix_val = NULL;
    show_mismatch = (Boolean)(option & TXALIGN_MISMATCH);
    
    
    /*process  the GAPs and the DIAGs segs*/
    s_start = -1;
    s_stop = -1;
    off_len = 0;
    for(asp = anp->segs; asp !=NULL; asp = asp->next) { 
        g_left = asp->gr.left;
        g_right = asp->gr.right;
        if(!(g_left > m_right || g_right < m_left)) {
            switch(asp->type) {
            case GAP_SEG:
                g_left = MAX(m_left, g_left);
                g_right = MIN(m_right, g_right);
                len = g_right - g_left +1;
                MemSet((Pointer)(l_seq +l_pos), '-',(size_t)len * sizeof(Char));
                l_pos += len;
                break;
                
            case REG_SEG:
            case DIAG_SEG:
            case STD_SEG:	/* Std-seg only works if the m_frame != 0 */
                if(m_left > g_left)
                    len = off_len + m_left - g_left;
                else
                    len = off_len;
                seq_offset = map_position_by_spacing(len, spacing, TRUE) * seq_expand;
                seq_start = anp->seqpos + seq_offset;
                g_left = MAX(m_left, g_left);
                g_right = MIN(m_right, g_right);
                len += (g_right - g_left);
                seq_stop = anp->seqpos + map_position_by_spacing(len, spacing, FALSE) * seq_expand + seq_expand -1;
                
                if(seq_start <= seq_stop) {	/*the order of start and stop is reversed*/
                    if(s_start == -1)	/*record the end point*/
                        s_start = ABS(seq_start);
                    s_stop = ABS(seq_stop);
                    
                    if(m_frame == 0)
                        fbuf_list = collect_feature_buf(asp->cnp, g_left, g_right, seq_start, l_pos, fbuf_list, maxlen, is_aa);	/*check the features first*/
                    load_text(bsp, seq_start, seq_stop, l_seq, &l_pos, m_buf, (Int2)maxlen, 
                              (Int2)spacing, translate, show_mismatch, matrix_val, matrix, strand, posMatrix, q_start);
                    
                }
                break;
                
            default:
                break;
            }
        }
        if(asp->type == INS_SEG)
            off_len += (asp->gr.right * spacing);
        if(asp->type == REG_SEG || asp->type == DIAG_SEG || asp->type == STD_SEG)
            off_len+=(asp->gr.right - asp->gr.left +1);
    }
    
    
    /*the first segment in the layout is a gap segment*/
    if(s_start == -1)
        s_start = *p_stop;
    if(s_stop == -1)	/*gap across the entire region*/
        s_stop = *p_stop;
    *p_stop = s_stop	/*update the stop value*/;
    tdp = (TextAlignBufPtr) MemNew(sizeof(TextAlignBuf));
    tdp->pos = s_start+1;
    tdp->strand = anp->extremes.strand;
    tdp->label = StringSave(anp->label);
    tdp->buf = l_seq;
    tdp->matrix_val = matrix_val;
    tdp->itemID = anp->itemID;
    tdp->feattype = 0;
    tdp->subtype = 0;
    tdp->entityID = anp->entityID;
    tdp->seqEntityID = anp->seq_entityID;
    tdp->bsp_itemID = anp->bsp_itemID;
    ValNodeAddPointer(&head, 0, tdp);
    ValNodeLink(&head, fbuf_list);
    
    ins_node = ProcessTextInsertion(anp, m_left, m_right, bsp, line_len, m_frame);
    ValNodeLink(&head, ins_node);
    BioseqUnlock(bsp);
    return head;
}

NLM_EXTERN ValNodePtr ProcessTextAlignNode(AlignNodePtr anp, Int4 m_left, Int4 m_right, Int4Ptr p_stop, CharPtr m_buf, Int4 line_len, Int1 m_frame, Uint4 option, Int4Ptr PNTR matrix)
{
    return ProcessTextAlignNode2(anp, m_left, m_right, p_stop, m_buf, line_len, m_frame, option, matrix, NULL, 0);
}
	
static void add_empty_space(CharPtr buf, Int4 maxlen)
{
	Int4 buf_len;

	buf_len = StringLen(buf);	
	if(buf_len < maxlen)
		make_empty(buf+buf_len, (Int2)(maxlen-buf_len));
}

static Int4 get_long_insert_len(Int4 length)
{
	Char temp[50];

	sprintf(temp, "(length=%ld)", (long) length);
	return (StringLen(temp) + 13);
}

static void copy_insertion_bar(CharPtr buf, CharPtr ins_2, Int2 sym_pos, Int4 len)
{
	Int2 k;

	if(buf == NULL || ins_2 == NULL)
		return;
	add_empty_space(buf, len);
	for(k = 0; k<sym_pos; ++k)
		if(ins_2[k] == '|' && buf[k] == ' ')
			buf[k] = '|';
}


/******************************************************************************
*
*	ProcessTextInsertion(anp, m_left, m_right, bsp)
*	convert the insertions that are located within [m_left, m_right] into
*	text buffer (a list of TextDrawPtr)
*	anp: AlignNodePtr
*	m_left, m_right: the current region for selection
*	bsp: the BioseqPtr for this anp
*
*	return a list of TextDrawPtr
*
******************************************************************************/
static ValNodePtr ProcessTextInsertion(AlignNodePtr anp, Int4 m_left, Int4 m_right, BioseqPtr bsp, Int4 line_len, Int1 m_frame)
{
	AlignSegPtr asp;
	Int4 inslen;		/*length of insertion*/
	Int2 insnum;		/*the number of insertions*/
	Int2 i, j;
	Int4Ptr inslevel;	/*for layout the level of insertions*/
	Int4 level;
	Int4 inspos;		/*position for insertion*/
	Int4 left;
	Int4 len;
	Int4 last_ins;

	CharPtr ins_1;	/* \ symbols for insertions*/
	CharPtr ins_2;	/*| symbols for insertion*/
	CharPtr ins_seq;
	Int4 sym_pos;
	Int4 l_pos;
	Int4 seq_offset, seq_start, seq_stop;
	ValNodePtr head = NULL;
	ValNodePtr fbuf_list = NULL, curr;
	TextAlignBufPtr fbp;
	Int4 g_left, g_right;
	Boolean is_aa;
	Int4 seq_expand;
	Int4 spacing;
	Boolean translate;
	Uint1 strand;

	strand = Seq_strand_plus;
	if(anp->seqpos < 0)
		strand = Seq_strand_minus;
	else if(anp->seqpos == 0 && anp->extremes.strand == Seq_strand_minus)
		strand = Seq_strand_minus;
	spacing = 1;
	if(m_frame > 0)
		spacing = 3;
	if(m_frame  == -1)
	{
		translate = TRUE;
		seq_expand = 3;
	}
	else
	{
		seq_expand = 1;
		translate = FALSE;
	}
	is_aa = (bsp->mol == Seq_mol_aa);
	insnum = 0;
	for(asp = anp->segs; asp !=NULL; asp = asp->next)
	/*checking the insertion numbers*/
	{ 
	   if(asp->type == INS_SEG)
	   {
	   	inspos = asp->ins_pos;
		if (inspos >= m_left && inspos<=m_right)
		{
			++insnum;
			asp->line = 0;
		}
		else
			asp->line = -1;
	   }
	}
	if(insnum == 0)
		return head;

	/*layout the insertions*/
	inslevel = (Int4Ptr) MemNew((size_t)(2*insnum) * sizeof(Int4));	/*layout insert*/
	level = 0;
	len = MAX(m_right - m_left +1, line_len);
	for(asp = anp->segs; asp !=NULL; asp = asp->next)
	{
	   if(asp->type == INS_SEG && asp->line == 0)
	   {
	   	inspos = asp->ins_pos;
		inslen = asp->gr.right/seq_expand;
		/* if(inslen > (m_right-m_left+1)) */
		if(inslen > len)
			inslen = get_long_insert_len(inslen);
		inspos -= m_left;
		asp->line = find_insert_ypos(&left, inslen, inspos, 0, len-1, inslevel, 2, insnum);
		asp->gr.left = left;
		level = MAX(asp->line, level);
	   }	
	}
	MemFree(inslevel);
	

	/*comput the insertion text*/
	for(j = 0; j< (level+1); ++j)
	{
	   l_pos = 0;
	   sym_pos = 0;
	   fbuf_list = NULL;
	   ins_seq = (CharPtr) MemNew((size_t)(len+1) * sizeof(Char));
	   ins_2 = (CharPtr) MemNew((size_t)(len+1) * sizeof(Char));
	   if(j == 0)
		ins_1 = (CharPtr) MemNew((size_t)(len+1) * sizeof(Char));
	   seq_offset = 0;
	   for(asp = anp->segs; asp !=NULL; asp = asp->next)
	   {
	   	if(asp->type == INS_SEG && asp->line >=j)
		{
	
			inspos = asp->ins_pos - m_left;
			if(inspos > sym_pos)
			{
			   if(j == 0)	/*the first level*/
				make_empty(ins_1+sym_pos, (Int2)(inspos-sym_pos));
			   make_empty(ins_2+sym_pos, (Int2)(inspos-sym_pos));
			   sym_pos = inspos;
			}
			if(j == 0)
				ins_1[sym_pos] = '\\';
			ins_2[sym_pos] = '|';
			if(asp->line == j)
				last_ins = inspos+1;
			++sym_pos;
			
			if(asp->line == j)
			{
			   seq_start = anp->seqpos + seq_offset;
			   seq_stop = seq_start + asp->gr.right -1;
			   /* seq_stop = seq_start + map_position_by_spacing(asp->gr.right, spacing, FALSE) * seq_expand + seq_expand -1; */
			   if(asp->gr.left > l_pos)
			   {
				make_empty(ins_seq+l_pos, (Int2)(asp->gr.left-l_pos));
				l_pos = asp->gr.left;
			   }

			   g_left = asp->ins_pos;
			   g_right = asp->ins_pos + asp->gr.right -1;
			   /* g_left = asp->gr.left;
			   g_right = g_left + asp->gr.right -1;*/

			   if((seq_stop - seq_start+1)>len)/*long insertions*/
			   {
		   	      fbuf_list = collect_feature_buf(asp->cnp, g_left, (g_left+4), seq_start, l_pos, fbuf_list, len, is_aa);	/*check the features first*/
		   	      fbuf_list = collect_feature_buf(asp->cnp, g_left, (g_left+4), seq_stop-4, l_pos+8, fbuf_list, len, is_aa);	/*check the features ffirst. 3 is the 3 dots*/
			   }
			   else
			      fbuf_list = collect_feature_buf(asp->cnp, g_left, g_right, seq_start, l_pos, fbuf_list, len, is_aa);

			   load_text(bsp, seq_start, seq_stop, ins_seq, &l_pos, NULL, (Int2)len, 1, translate, FALSE, NULL, NULL, strand, NULL, 0);
			}

		}
		if(asp->type == INS_SEG)
			seq_offset += asp->gr.right;
		if(asp->type == DIAG_SEG || asp->type == REG_SEG || asp->type == STD_SEG)
			seq_offset += map_position_by_spacing(asp->gr.right - asp->gr.left +1, 
				spacing, TRUE) * seq_expand; 
			/* seq_offset += (asp->gr.right - asp->gr.left +1) * seq_expand; */
	   }

	   ins_2[sym_pos] = '\0';
	   ins_seq[l_pos] = '\0';
	   if(j == 0)
	   {
		ins_1[sym_pos] = '\0';	
		load_tdp_data(&head, NULL, ins_1, 0, 0, 0, 0);
	   }

	   for(curr = head; curr !=NULL; curr = curr->next)
	   /*for(curr = fbuf_list; curr !=NULL; curr = curr->next)*/
	   {
		fbp = (TextAlignBufPtr) curr->data.ptrvalue;
		if(fbp->buf != NULL)
			copy_insertion_bar(fbp->buf, ins_2, (Int2)sym_pos, len);
		else
		{
			for(i =0; i<3; ++i)
				copy_insertion_bar(fbp->codon[i], ins_2, (Int2)sym_pos, len);
		}
	   }

	   copy_insertion_bar(ins_seq, ins_2, (Int2)sym_pos, len);
	   load_tdp_data(&head, NULL, ins_2, 0, 0, 0, 0);
	   load_tdp_data(&head, NULL, ins_seq, anp->itemID, anp->entityID, anp->seq_entityID, anp->bsp_itemID);
	   ValNodeLink(&head, fbuf_list);
	   fbuf_list = head;
	}

	return head;

}


static Char pchars[] = "ARNDCQEGHILKMFPSTWYVBZX";	/* amino acid names */
static Int4 webb_blosum62[WEBB_asize][WEBB_asize] = {
   { 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,-2,-1, 0 },
   {-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,-1, 0,-1 },
   {-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3, 3, 0,-1 },
   {-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3, 4, 1,-1 },
   { 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2 },
   {-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2, 0, 3,-1 },
   {-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1 },
   { 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,-1,-2,-1 },
   {-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3, 0, 0,-1 },
   {-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-3,-3,-1 },
   {-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-4,-3,-1 },
   {-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2, 0, 1,-1 },
   {-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-3,-1,-1 },
   {-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-3,-3,-1 },
   {-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-2,-1,-2 },
   { 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2, 0, 0, 0 },
   { 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,-1,-1, 0 },
   {-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,-4,-3,-2 },
   {-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-3,-2,-1 },
   { 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-3,-2,-1 },
   {-2,-1, 3, 4,-3, 0, 1,-1, 0,-3,-4, 0,-3,-3,-2, 0,-1,-4,-3,-3, 4, 1,-1 },
   {-1, 0, 0, 1,-3, 3, 4,-2, 0,-3,-3, 1,-1,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1 },
   { 0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-1,-1,-1 },
 };

NLM_EXTERN Int4Ptr PNTR load_default_matrix (void)
{
	Int4Ptr PNTR ss;
	Int2 i, j;

	ss = (Int4Ptr PNTR) MemNew((size_t)TX_MATRIX_SIZE * sizeof (Int4Ptr));
	for(i = 0; i<TX_MATRIX_SIZE; ++i)
		ss[i] = (Int4Ptr) MemNew((size_t)TX_MATRIX_SIZE * sizeof (Int4));

	for(i = 0; i < TX_MATRIX_SIZE; i++) 
		for(j = 0; j < TX_MATRIX_SIZE;j++) 
			ss[i][j] = -1000;
	for(i = 0; i < WEBB_asize; ++i)
		for(j = 0; j < WEBB_asize; ++j)
			ss[pchars[i]][pchars[j]] = webb_blosum62[i][j];
	for(i = 0; i < WEBB_asize; ++i)
		ss[pchars[i]]['*'] = ss['*'][pchars[i]] = -4;
	ss['*']['*'] = 1;
	return ss;
}

NLM_EXTERN void free_default_matrix (Int4Ptr PNTR matrix)
{
	Int2 i;

	for(i = 0; i<TX_MATRIX_SIZE; ++i)
		MemFree(matrix[i]);
	MemFree(matrix);
}



/*
	Adds tdsp to the end of a chain of TxDfLineStructPtr's.
	Returns the new TxDfLineStructPtr.
*/

static TxDfLineStructPtr
TxDfLineStructAdd(TxDfLineStructPtr PNTR head, TxDfLineStructPtr tdsp)

{
	TxDfLineStructPtr var;

	if (*head == NULL)
	{
		*head = tdsp;
	}
	else
	{
		var = *head;
		while (var->next)
		{
			var = var->next;
		}
		var->next = tdsp;
	}

	return tdsp;
}


static SeqIdPtr 
ScorePtrUseThisGi (ScorePtr sp)

{
    ObjectIdPtr obid;
    ScorePtr scrp;
    SeqIdPtr gilist=NULL;
    
    for (scrp=sp; scrp; scrp = scrp->next) {
        obid = scrp->id;
        if(obid && obid->str) {
            if (StringICmp(obid->str, "use_this_gi") == 0) {
                ValNodeAddInt(&gilist, SEQID_GI, scrp->value.intvalue);
            }
        }
    }
    
    return gilist;
}

/*
	GetUseThisGi(SeqAlignPtr) looks for the "use_this_gi" flag in the ScorePtr.
*/

NLM_EXTERN SeqIdPtr LIBCALL
GetUseThisGi(SeqAlignPtr seqalign)

{
	Boolean retval=FALSE;
        DenseDiagPtr ddp;
        DenseSegPtr dsp;
        ScorePtr sp;
	SeqIdPtr gilist=NULL;
        StdSegPtr ssp;
	
	sp = seqalign->score;
	if (sp == NULL)
	{
		switch (seqalign->segtype)
		{
			case 1: /*Dense-diag*/
				ddp = (DenseDiagPtr) seqalign->segs;
				while (ddp)
				{
					sp = ddp->scores;
					if (sp)
						break;
					ddp = ddp->next;
				}
				break;
			case 2:
				dsp = ( DenseSegPtr) seqalign->segs;
				if (dsp)
				{
					sp = dsp->scores;
				}
				break;
			case 3:
				ssp = (StdSegPtr) seqalign->segs;
				while (ssp)
				{
					sp = ssp->scores;
					if (sp)
						break;
					ssp = ssp->next;
				}
				break;
			default:
				break;
		}
	}


	gilist = ScorePtrUseThisGi(sp);
	return gilist;
}


/*
  Filters the FASTA definition lines based on the SeqIdPtr's given 
  as input. The gi_list is used to pull sequences out of the BioseqPtr, 
  The buffer_id contains FASTA formatted ID, title contains the rest 
  of the title.
*/

NLM_EXTERN Boolean LIBCALL
FilterTheDefline (BioseqPtr bsp, SeqIdPtr gi_list_head, CharPtr buffer_id, Int4 buffer_id_length, CharPtr PNTR titlepp)

{
    Boolean first_time, found_gi, found_first_gi, not_done;
    CharPtr bsp_title, bsp_title_ptr, title, title_ptr;
    Char buffer[BUFFER_LENGTH], id_buf[255];
    Int4 index;
    SeqIdPtr gi_list, sip;
    
    if (bsp == NULL || gi_list_head == NULL)
        return FALSE;
    
    bsp_title = BioseqGetTitle(bsp);
    bsp_title_ptr = bsp_title;
    /* This is the longest it could be, this could be done more efficiently. */
    title = (CharPtr) MemNew((256+StringLen(bsp_title))*sizeof(Char));
    title_ptr = title;
    *titlepp = title;
    
    /*
      if (bsp_title_ptr == NULL)
      return FALSE;
    */
    
    first_time = TRUE;
    found_first_gi = TRUE;
    not_done = TRUE;
    while (not_done) {
        if (!first_time) {
            index=0;
            id_buf[0] = NULLB;
            if (bsp_title_ptr) {
                while (*bsp_title_ptr != NULLB) {
                    if (*bsp_title_ptr == ' ') {
                        id_buf[index] = NULLB;
                        break;
                    }
                    id_buf[index] = *bsp_title_ptr;
                    bsp_title_ptr++;
                    index++;
                }
            }
            if (id_buf[0] == NULLB)
                break;
            sip = SeqIdParse(id_buf);
        } else {
            sip = bsp->id;
        }
	
        found_gi = FALSE;
        gi_list = gi_list_head;
        while (gi_list) {
            if(SeqIdIn(gi_list, sip) == TRUE) {
                found_gi = TRUE;
                break;
            }
            gi_list = gi_list->next;
        }
        
        if (found_gi) {
            if (!found_first_gi) {
                *title_ptr = '>';
                title_ptr++;
                SeqIdWrite(sip, buffer, PRINTID_FASTA_LONG, BUFFER_LENGTH);
                StringCpy(title_ptr, buffer);
                title_ptr += StringLen(buffer);
            } else {
                SeqIdWrite(sip, buffer_id, PRINTID_FASTA_LONG, buffer_id_length);
                found_first_gi = FALSE;
            }
            
            if (bsp_title_ptr) {
                while (*bsp_title_ptr != '>' && *bsp_title_ptr != NULLB) {
                    *title_ptr = *bsp_title_ptr;
                    bsp_title_ptr++;
                    title_ptr++;
                }
            }
        } else {
            if (bsp_title_ptr) {
                while (*bsp_title_ptr != '>' && *bsp_title_ptr != NULLB)
                    bsp_title_ptr++;
            }
        }
        
        if (first_time) {
            first_time = FALSE;
        } else {
            sip = SeqIdSetFree(sip);
        }
        
        if (bsp_title_ptr) {
            if (*bsp_title_ptr == '>')
                bsp_title_ptr++;
            
            if (*bsp_title_ptr == NULLB) {
                *title_ptr = NULLB;
                break;
            }
        }
    }
    return TRUE;
}


#define STATS_LENGTH 14
NLM_EXTERN Boolean LIBCALL
PrintDefLinesFromAnnot(SeqAnnotPtr seqannot, Int4 line_length, FILE *outfp, Uint4 options, Int4 mode, Int2Ptr marks)

{
	Boolean retval;

	if (seqannot == NULL || seqannot->type != 2)
	{
		return FALSE;
	}

	retval = PrintDefLinesFromSeqAlign((SeqAlignPtr) seqannot->data, line_length, outfp, options, mode, marks);

	return retval;
}

NLM_EXTERN Boolean LIBCALL
PrintDefLinesExtra(ValNodePtr vnp, Int4 line_length, FILE *outfp, Uint4 options, Int4 mode, Int2Ptr marks, SeqLocPtr seqloc)

{
	Boolean retval;
	Char buffer[128];
	Int2 titleIdAllocated;
	Int4 index=0;
	SeqAlignPtr seqalign;

	if (vnp == NULL || seqloc == NULL)
	{
		return FALSE;
	}

	/* Disable printing of title. */
	if (!(options & TXALIGN_DO_NOT_PRINT_TITLE)) 
		options += TXALIGN_DO_NOT_PRINT_TITLE;


	asn2ff_set_output(outfp, NULL);

	ff_StartPrint(0, 0, (Int2)(line_length+2), NULL);

	titleIdAllocated = line_length - STATS_LENGTH;

	NewContLine();
	TabToColumn((Int2)(titleIdAllocated+2));
	ff_AddString("Score     E");
	NewContLine();
	TabToColumn((Int2)(titleIdAllocated+2));
	if (options & TXALIGN_SHOW_NO_OF_SEGS) {
	    ff_AddString("(bits)  Value  N");
	} else {
	   ff_AddString("(bits)  Value");
	}
	ff_EndPrint(); 

        while (vnp && seqloc)
        {
		ff_StartPrint(0, 0, (Int2)(line_length+2), NULL);
		index++;
                seqalign = (SeqAlignPtr) vnp->data.ptrvalue;
		sprintf(buffer, "\nSignificant matches for pattern occurrence %ld at position %ld\n\n",
                        (long) index, (long) (SeqLocStart(seqloc)+1));
                ff_AddString(buffer);
		ff_EndPrint(); 
		retval = PrintDefLinesFromSeqAlign(seqalign, line_length, outfp, options, mode, marks);
		vnp = vnp->next;
		seqloc = seqloc->next;
	}

	return retval;
}

NLM_EXTERN Boolean LIBCALL
PrintDefLinesFromSeqAlignEx(SeqAlignPtr seqalign, Int4 line_length, FILE *outfp, Uint4 options,
		Int4 mode, Int2Ptr marks, Int4 number_of_descriptions)
{
	return PrintDefLinesFromSeqAlignEx2(seqalign, line_length, outfp, options,
			mode, marks, number_of_descriptions, (CharPtr)NULL, (CharPtr)NULL);
}

/**
 * transforms a string so that it becomes safe to be used as part of URL
 * the function converts characters with special meaning (such as
 * semicolon -- protocol separator) to escaped hexadecimal (%xx)
 */
static
CharPtr
MakeURLSafe(CharPtr src)
{
	static Char HEXDIGS[] = "0123456789ABCDEF";
	CharPtr buf;
	size_t len;
	CharPtr p;
	Char c;

	if (src == NULL) {
		return NULL;
	}
	/* first pass to calculate required buffer size */
	for (p = src, len = 0; (c = *(p++)) != '\0'; ) {
		switch (c) {
		default:
			if (c < '0' || (c > '9' && c < 'A') ||
					(c > 'Z' && c < 'a') || c > 'z') {
				len += 3;
				break;
			}
		case '-': case '_': case '.': case '!': case '~':
		case '*': case '\'': case '(': case ')':
			++len;
		}
	}
	buf = (CharPtr)MemNew(len + 1);
	/* second pass -- conversion */
	for (p = buf; (c = *(src++)) != '\0'; ) {
		switch (c) {
		default:
			if (c < '0' || (c > '9' && c < 'A') ||
					(c > 'Z' && c < 'a') || c > 'z') {
				*(p++) = '%';
				*(p++) = HEXDIGS[(c >> 4) & 0xf];
				*(p++) = HEXDIGS[c & 0xf];
				break;
			}
		case '-': case '_': case '.': case '!': case '~':
		case '*': case '\'': case '(': case ')':
			*(p++) = c;
		}
	}
	*p = '\0';
	return buf;
}

static
CharPtr
StringAppend(CharPtr *dst, size_t *size, CharPtr src, size_t *used)
{
	size_t pos, len;

	if (*dst == NULL) {
		*size = 1;
		pos = 0;
	}
	else {
		pos = *used;
	}
	if (src == NULL) {
		return *dst;
	}
	len = StringLen(src);
	*used += len;
	if (*dst == NULL || pos + len + 1 > *size) {
		/**
		 * extending destination buffer
		 */
		CharPtr old = *dst;
		for (; pos + len + 1 > *size; *size *= 2);
		*dst = (CharPtr)MemNew(*size);
		**dst = '\0';
		if (old != NULL) {
			StringCpy(*dst, old);
			MemFree(old);
		}
	}
	StringCpy((*dst) + pos, src);
	return *dst;
}


static
Boolean
SeqAlignSegsStr(SeqAlignPtr sap, Int2 index, CharPtr *dst, size_t *size, size_t *used)
{
	Char buf[128];
	Int4 start, stop;

	start = SeqAlignStart(sap, 1);
	stop = SeqAlignStop(sap, 1);

	if (sap == NULL) {
		return FALSE;
	}

	sprintf(buf, "%ld-%ld", (long)(start), (long)(stop));
	StringAppend(dst, size, buf, used);

	return TRUE;
}

/**
	* links to incomplete genomes
**/
static void
make_dumpgnl_links(SeqIdPtr sip, CharPtr blast_type, CharPtr segs, CharPtr dbname, Boolean is_na, FILE *fp, CharPtr sip_buffer, Boolean nodb_path)
{
    BioseqPtr bsp;
    Char gnl[256];
    CharPtr str, chptr, dbtmp;
    Uchar buf[32];
    Int4 i, j, length, gi;
    MD5Context context;
    Char passwd[128], tool_url[128], tmpbuff[256];
    SeqIdPtr bestid;
    
    /* We do need to make security protected link to BLAST gnl */
    if (StringStr(sip_buffer, "gnl|BL_ORD_ID") != NULL)
        return;
    
    *passwd = NULLB;
    *tool_url = NULLB;
    
    str = NULL;
#ifdef OS_UNIX
    str = getenv("DUMPGNL_PASSWD");
#endif
    if(str != NULL) {
        StringCpy(passwd, str);
    } else {
        GetAppParam("NCBI", blast_type, "PASSWD", "", passwd, 
                    sizeof(passwd));
    }
    
    str = NULL;
#ifdef OS_UNIX
    str = getenv("DUMPGNL_TOOL_URL");
#endif
    if(str != NULL) {
        StringCpy(tool_url, str);
    } else {
        GetAppParam("NCBI", blast_type, "TOOL_URL", "", tool_url, 
                    sizeof(tool_url));
    }

    if(*passwd == NULLB || *tool_url == NULLB)
        return;
    
    if(nodb_path) {

        length = StringLen(dbname);
        dbtmp = MemNew(sizeof(Char)*length + 2); /* aditional space and NULLB */
        
        for(i = 0; i < length; i++) {
            
            if(isspace(dbname[i]) || dbname[i] == ',') /* Rolling spaces */
                continue;

            j = 0;
            while (!isspace(dbname[i]) && j < 256  && i < length) { 
                tmpbuff[j] = dbname[i];
                j++; i++;
                if(dbname[i] == ',')  /* Comma is valid delimiter */
                    break;
            }
            tmpbuff[j] = NULLB;

            if((chptr = strrchr(tmpbuff, '/')) != NULL) { 
                StringCat(dbtmp, chptr+1);
            } else {
                StringCat(dbtmp, tmpbuff);
            }

            StringCat(dbtmp, " ");            
        }
    } else {
        dbtmp = dbname;
    }
    
    bsp = BioseqLockById(sip);
    if (bsp)
	sip = bsp->id;

    bestid = SeqIdFindBest(sip, SEQID_GENERAL);
    if (bestid && bestid->choice != SEQID_GENERAL)
    {
    	bestid = SeqIdFindBest(sip, SEQID_OTHER);
    }
    /*
     * Need to protect start and stop positions
     * to avoid web users sending us hand-made URLs
     * to retrive full sequences
     */
    if (bestid && (bestid->choice == SEQID_GENERAL || bestid->choice == SEQID_OTHER))
    {
    	MD5Init(&context);
    	length = StringLen(passwd);
	MD5Update(&context, (UcharPtr)passwd, (Uint4)length);
	SeqIdWrite(bestid, gnl, PRINTID_FASTA_SHORT, sizeof(gnl));
	MD5Update(&context, (UcharPtr)gnl, (Uint4)StringLen(gnl));
	MD5Update(&context, (UcharPtr)segs, (Uint4)StringLen(segs));
	MD5Update(&context, (UcharPtr)passwd, (Uint4)length);
	MD5Final(&context, (UcharPtr)buf);
    }
    else
    {
	gnl[0] = NULLB;
    }

    gi = -1;
    bestid = SeqIdFindBest(sip, SEQID_GI);
    if (bestid && bestid->choice == SEQID_GI)
    {
	gi = bestid->data.intvalue;
    }
    
    str = MakeURLSafe(dbtmp == NULL ? "nr" : dbtmp);
    fprintf(fp, "<a href=\"%s?db=%s&na=%d&", tool_url, str, is_na);
    str = (CharPtr) MemFree(str);
    if (gnl[0] != NULLB)
    {
    	str = MakeURLSafe(gnl);
    	fprintf(fp, "gnl=%s&", str);
    	str = (CharPtr) MemFree(str);
    }
    if (gi != -1)
    {
    	fprintf(fp, "gi=%ld&", (long) gi);
    }
    if (RID_glb)
    {
    	fprintf(fp, "RID=%s&", RID_glb);
    }

    fprintf(fp,
            "segs=%s&seal=%02X%02X%02X%02X"
            "%02X%02X%02X%02X%02X%02X%02X%02X%02X%02X%02X%02X\">",
            segs,
            buf[0], buf[1], buf[2], buf[3],
            buf[4], buf[5], buf[6], buf[7],
            buf[8], buf[9], buf[10], buf[11],
            buf[12], buf[13], buf[14], buf[15]);
    
    BioseqUnlock(bsp);
    if(nodb_path)
        MemFree(dbtmp);
    
    return;
}

NLM_EXTERN Boolean LIBCALL
PrintDefLinesFromSeqAlignEx2(SeqAlignPtr seqalign, Int4 line_length, FILE *outfp, Uint4 options,
		Int4 mode, Int2Ptr marks, Int4 number_of_descriptions,
		CharPtr db_name, CharPtr blast_type)

{
    BioseqPtr bsp;
    Boolean found_next_one, found_gnl_id, same_id, found_score=FALSE, make_link=FALSE;
    Char buffer[BUFFER_LENGTH+1], buffer1[BUFFER_LENGTH+1], eval_buff[10], bit_score_buff[10]; 
    Char HTML_buffer[BUFFER_LENGTH+1], HTML_database[32], HTML_dopt[16], id_buffer[BUFFER_LENGTH+1];
    Char *ptr, *ptr_start, *eval_buff_ptr, *bit_score_buff_ptr;
    Int2 pos, title_length, title_allocated, titleIdAllocated;
    Nlm_FloatHi bit_score, evalue;
    Int4 gi = 0, number, score;
    SeqIdPtr bestid, gi_list, subject_id, sip_list=NULL, last_id;
    TxDfLineStructPtr txsp = NULL, txsp_head, txsp_var;
    Boolean retval = FALSE;
    Boolean firstnew = TRUE, nodb_path = FALSE;
    Int4 countdescr = number_of_descriptions;
    Int4 numalign;
    Char passwd[128];
    Char tool_url[128];
    DbtagPtr db_tag;
    ObjectIdPtr oip;
    CharPtr www_root_path = NULL;

    if (outfp == NULL) {
        return FALSE;
    }

    if (seqalign == NULL) {	/* Two line returns so that the alignments or db report is not all bunched up. */
        NewContLine();
        NewContLine();
        return FALSE;
    }

#ifdef OS_UNIX
    www_root_path = getenv("WWW_ROOT_PATH");
#endif
    
    if(!StringICmp(blast_type, "fruitfly")) {
        fprintf(stdout, "<IMG SRC=\"/BLAST/images/map_mark.gif\" BORDER=0> - please follow this image for the map location of the sequence<P>\n");
    }
    
    asn2ff_set_output(outfp, NULL);
    
    ff_StartPrint(0, 0, (Int2)(line_length+2), NULL);
    
    titleIdAllocated = line_length - STATS_LENGTH;
    
    if (options & TXALIGN_SHOW_NO_OF_SEGS) {
        titleIdAllocated -= 4;
    }
    
    if (options & TXALIGN_CHECK_BOX) {
        titleIdAllocated += 2;
    }

    if(options & TXALIGN_NEW_GIF)
        titleIdAllocated += 3;
    
    /* <PRE> block should be already opened outside of this function, 
       but open it here just in case */
    if (options & TXALIGN_HTML) {
       ff_AddString("<PRE>");
       NewContLine();
    }

    /*AAS*/
    if (!(options & TXALIGN_DO_NOT_PRINT_TITLE)) {
        if ((mode == FIRST_PASS) || (mode == NOT_FIRST_PASS_REPEATS)) {
            NewContLine();
            NewContLine();
            TabToColumn((Int2)(titleIdAllocated+2));
            
            ff_AddString("Score     E");
            NewContLine();
            ff_AddString("Sequences producing significant alignments:");
            TabToColumn((Int2)(titleIdAllocated+2));
            if (options & TXALIGN_SHOW_NO_OF_SEGS) {
                ff_AddString("(bits)  Value  N");
            } else {
                ff_AddString("(bits)  Value");
            }
            NewContLine();
        }
        
        if (mode == NOT_FIRST_PASS_REPEATS) {
            ff_AddString("Sequences used in model and found again:");
            NewContLine();
        }
        if (mode == NOT_FIRST_PASS_NEW) {
            ff_AddString("Sequences not found previously or not previously below threshold:");
            NewContLine();
        }
        ff_EndPrint(); 
    }
    
    numalign = 0;
    last_id = NULL;
    txsp_head = NULL;
    while (seqalign) {
        if ((mode == FIRST_PASS) ||
            ((mode == NOT_FIRST_PASS_REPEATS) && marks && marks[numalign] & SEQ_ALIGN_MARK_REPEAT) ||
            ((mode == NOT_FIRST_PASS_NEW) && marks && (!(marks[numalign] & SEQ_ALIGN_MARK_REPEAT)))) {
            
            subject_id = SeqIdDup(TxGetSubjectIdFromSeqAlign(seqalign));
            same_id = FALSE;
            if(last_id && SeqIdComp(subject_id, last_id) == SIC_YES) {
                same_id = TRUE;
            }
            
            last_id = SeqIdFree(last_id);
            last_id = SeqIdDup(subject_id);
            
            found_score = GetScoreAndEvalue(seqalign, &score, &bit_score, &evalue, &number);
            /* if the ID has been seen before, check that proper values are saved. */
            if (same_id == TRUE) {
                if (score > txsp->score)
                    txsp->score = score;
                if (bit_score > txsp->bit_score)
                    txsp->bit_score = bit_score;
                if (evalue < txsp->evalue)
                    txsp->evalue = evalue;
                if (number < txsp->number)
                    txsp->number = number;
                StringAppend(&txsp->segs_str, &txsp->segs_buflen, ",", &txsp->segs_used);
                SeqAlignSegsStr(seqalign, 1, &txsp->segs_str, &txsp->segs_buflen, &txsp->segs_used);
                subject_id = SeqIdFree(subject_id);
            } else {
                bsp = BioseqLockById(subject_id);	
                txsp = (TxDfLineStructPtr) MemNew(sizeof(TxDfLineStruct));
                txsp->segs_str = NULL;
                txsp->segs_buflen = 0;
                if(bsp != NULL) {
                    gi_list = GetUseThisGi(seqalign);
                    if (gi_list) {
                        FilterTheDefline(bsp, gi_list, buffer, BUFFER_LENGTH, &(txsp->title));
                        gi_list = SeqIdSetFree(gi_list);
                        subject_id = SeqIdFree(subject_id);
                        txsp->id = SeqIdParse(buffer);
                    } else {
                        SeqIdWrite(bsp->id, buffer, PRINTID_FASTA_LONG, BUFFER_LENGTH);
                        txsp->title = StringSave(BioseqGetTitle(bsp));
                        txsp->id = subject_id;
                    }
                    txsp->is_na = (bsp->mol != Seq_mol_aa);
                } else {
                    SeqIdWrite(subject_id, buffer, PRINTID_FASTA_LONG, BUFFER_LENGTH);
                    txsp->title = StringSave("Unknown");
                    txsp->is_na = FALSE;
                    txsp->id = subject_id;
                }
                txsp->seqalign = seqalign;
                txsp->buffer_id = StringSave(buffer);
                txsp->score = score;
                txsp->bit_score = bit_score;
                txsp->evalue = evalue;
                txsp->number = number;
                txsp->found_score = found_score;
                SeqAlignSegsStr(seqalign, 1, &txsp->segs_str, &txsp->segs_buflen, &txsp->segs_used);
                if (marks) {
                    /* seq is new if it was not Good on previous iteration */
                    txsp->isnew = (Boolean) !(marks[numalign] & SEQ_ALIGN_MARK_PREVGOOD);
                    txsp->waschecked = (Boolean) marks[numalign] & SEQ_ALIGN_MARK_PREVCHECKED;
                } else {
                    txsp->isnew = FALSE;
                    txsp->waschecked = FALSE;
                }
                
                txsp = TxDfLineStructAdd(&txsp_head, txsp);
                if(bsp != NULL)
                    BioseqUnlock(bsp);
                retval = TRUE;
            }
        }
        seqalign = seqalign->next;
        numalign++;
    }
    last_id = SeqIdFree(last_id);
    
    if(retval == FALSE)
        return FALSE;
    
    
    /* Used for dumpgnl reports if GNL id's. (overwrite parameter!) */
    if (blast_type == NULL)  {
        blast_type = StringSave("UNFIN_GEN");
    } else {
        blast_type = StringSave(blast_type);
        StringUpper(blast_type);
    }

    /* If option TXALIGN_NO_ENTREZ set full database name will be stripped
       to the database fileneme */
    
    if(options & TXALIGN_NO_ENTREZ)
        nodb_path = TRUE;

    txsp = txsp_head;
    while (txsp) {
        found_next_one = FALSE;
        if (options & TXALIGN_HTML) {

            if (txsp->is_na) {
                StringCpy(HTML_dopt, "GenBank");
                StringCpy(HTML_database, "Nucleotide");
            } else {
                StringCpy(HTML_dopt, "GenPept");
                StringCpy(HTML_database, "Protein");
            }
            gi = 0;
            make_link = FALSE;
            bestid = SeqIdFindBest(txsp->id, SEQID_GI);
            if (bestid != NULL && bestid->choice == SEQID_GI && !(options & TXALIGN_NO_ENTREZ)) {
                gi = bestid->data.intvalue;
                if (options & TXALIGN_CHECK_BOX && options & TXALIGN_CHECK_BOX_CHECKED) {
                    if (countdescr == -1 || countdescr > 0) {

   sprintf(HTML_buffer, 
           "<INPUT TYPE=\"checkbox\" NAME=\"checked_GI\" "
           "VALUE=\"%d\" CHECKED> "
           "<a href=%s?cmd=Retrieve&db=%s&list_uids=%08ld&dopt=%s>"
           "<INPUT TYPE=\"hidden\" NAME =\"good_GI\" VALUE = \"%d\">",
           gi, NEW_ENTREZ_HREF, HTML_database, (long) gi, HTML_dopt, gi);

                                    } else {
   sprintf(HTML_buffer, 
           "<INPUT TYPE=\"hidden\" NAME=\"checked_GI\" VALUE=\"%d\">"
           "<INPUT TYPE=\"hidden\" NAME =\"good_GI\" VALUE = \"%d\">",
           gi, gi);
                                    }
                } else if (options & TXALIGN_CHECK_BOX) {
                    if (countdescr == -1 || countdescr > 0) 

   sprintf(HTML_buffer, 
           "<INPUT TYPE=\"checkbox\" NAME=\"checked_GI\" VALUE=\"%d\"> "
           "<a href=\"%s?cmd=Retrieve&db=%s&list_uids=%08ld&dopt=%s\">",
           gi, NEW_ENTREZ_HREF, HTML_database, (long) gi, HTML_dopt);

                } else {
                    if (countdescr == -1 || countdescr > 0)
   if(!StringICmp(blast_type, "fruitfly")) {
       sprintf(HTML_buffer, 
               "<a href=\"http://www.ncbi.nlm.nih.gov\">"
               "<IMG SRC=\"/BLAST/images/map_mark.gif\" BORDER=0></a>"
               "&nbsp;&nbsp;<a href=\"%s?"
               "cmd=Retrieve&db=%s&list_uids=%08ld&dopt=%s\">", 
               NEW_ENTREZ_HREF, HTML_database, (long) gi, HTML_dopt);
   } else {

   sprintf(HTML_buffer, 
           "<a href=\"%s?cmd=Retrieve&db=%s&list_uids=%08ld&dopt=%s\">", 
           NEW_ENTREZ_HREF, HTML_database, (long) gi, HTML_dopt);
   }
                }

                if (options & TXALIGN_NEW_GIF && (countdescr == -1 || countdescr > 0)) {
                    if (txsp->isnew) {
                        if (firstnew) {
                            firstnew = FALSE;
                            fprintf(outfp, "<a name = Evalue></a>");
                        }
                        fprintf(outfp, "<br><IMG SRC=\"%s/images/new.gif\" WIDTH=25 HEIGHT=15>", www_root_path == NULL? "/blast" : www_root_path);
                    } else {
                        fprintf(outfp, "<br><IMG SRC=\"%s/images/bg.gif\" WIDTH=25 HEIGHT=15>",  www_root_path == NULL? "/blast" : www_root_path);
                    }
                    if (txsp->waschecked) {
                        fprintf(outfp, "<IMG SRC=\"%s/images/checked.gif\" WIDTH=15 HEIGHT=15>",  www_root_path == NULL? "/blast" : www_root_path);
                    } else {
                        fprintf(outfp, "<IMG SRC=\"%s/images/bg.gif\" WIDTH=15 HEIGHT=15>",  www_root_path == NULL? "/blast" : www_root_path);
                    }
                }
                fprintf(outfp, "%s", HTML_buffer);
                make_link = TRUE;
                /* If not SEQID_GI */
            } else if (bestid != NULL && (bestid->choice == SEQID_GENERAL && !(options & TXALIGN_NO_DUMPGNL)) || (options & TXALIGN_NO_ENTREZ))
                {
                    if (bestid->choice != SEQID_GENERAL && bestid->choice != SEQID_OTHER)
                        { /* HACK, HACK, use SEQID_GENERAL for Greg's page, even though GI is present. */
                            bsp = BioseqLockById(bestid);	
                            bestid = SeqIdFindBest(bsp->id, SEQID_OTHER);
                            BioseqUnlock(bsp);
                        }
                    if (bestid->choice == SEQID_GENERAL) {
                        db_tag = (DbtagPtr) bestid->data.ptrvalue;
                        if(db_tag->db && StringCmp(db_tag->db, "THC") == 0) {
                            oip = db_tag->tag;
                            if(oip->id != 0) {
                                fprintf(outfp, "<a href=\"http://www.tigr.org/docs/tigr-scripts/hgi_scripts/thc_report.spl?est=THC%ld&report_type=n\">", (long) oip->id);
                                
                            }
                        } else {
                            make_dumpgnl_links(txsp->id, blast_type, txsp->segs_str, db_name, ISA_na(bsp->mol), outfp, txsp->buffer_id, nodb_path);
                        }
                    } else
                        make_dumpgnl_links(txsp->id, blast_type, txsp->segs_str, db_name, ISA_na(bsp->mol), outfp, txsp->buffer_id, nodb_path);
                    make_link = TRUE;
                }
        }
        
        sprintf(buffer, "%s", txsp->buffer_id);
        if (!(options & TXALIGN_SHOW_GI)) {
            if (StringNCmp(buffer, "gi|", 3) == 0) {
                ptr = &buffer[3];
                while (*ptr != NULLB && *ptr != ' ') { /* Is there another ID beside the GI? */
                    if (*ptr == '|') {
                        ptr++;
                        found_next_one = TRUE;
                        break;
                    }
                    ptr++;
                }
            }
        }
        if (found_next_one == FALSE) {
            ptr = buffer;
            ptr_start = buffer;
        } else {
            ptr_start = ptr;
        }
        
        found_gnl_id = FALSE;
        /* Check for an ID of type general from BLAST */

        if (StringNCmp(buffer, "gnl|BL_ORD_ID", 13) == 0) {
            ptr = buffer;
            /* look for end of gnl ID. */
            while (*ptr != NULLB && *ptr != ' ')
                ptr++;
            /* Clear out all spaces. */	
            while (*ptr != NULLB && *ptr == ' ')
                ptr++;
            ptr_start = ptr;
            found_gnl_id = TRUE;
            make_link = FALSE;
        }

        if (StringNCmp(ptr, "lcl|", 4) == 0) {
            ptr_start += 4;
        }
        
        pos = StringLen(ptr);
        
        if ((options & TXALIGN_HTML) && make_link) {
            StringCpy(ptr+pos, "</a> ");
            pos++;		/* One for the space after "</a>" */
            pos += 4;	/* for "</a>" */
            title_allocated = titleIdAllocated - pos;
        }
        
        title_allocated = titleIdAllocated - pos;
        
        if (pos >= titleIdAllocated) {
            pos = titleIdAllocated+1; /* no space to definition. */
            sprintf(ptr+pos-3, "...");
            *(ptr+pos) = ' ';
            pos++;
            *(ptr+pos) = NULLB; /* in case no scores are printed. */
        } else  {
            if (found_gnl_id == FALSE) {
                *(ptr + pos) = ' ';
                pos++;
            }
            
            title_length = StringLen(txsp->title);
            if (title_length > title_allocated) {
                title_length = title_allocated;
                title_length -= 3;	/* For "..." */
                if (txsp->title) {
                    StringNCpy((ptr+pos), txsp->title, title_length);
                    pos += title_length;
                }
                sprintf((ptr+pos), "...");
                pos += 3;
            } else {
                if (txsp->title) {
                    StringNCpy((ptr+pos), txsp->title, title_length);
                    pos += title_length;
                } 
                while (title_length < title_allocated) {
                    *(ptr + pos) = ' ';
                    title_length++;
                    pos++;
                }
            }
            *(ptr + pos) = ' ';
            pos++;
            
            /* set to NULLB in case no scores have been found. */
            *(ptr + pos) = NULLB;
        }
        
#ifdef OS_MAC
        if (txsp->found_score) {
            evalue = txsp->evalue;
            eval_buff_ptr = eval_buff;
            if (evalue < 1.0e-180) {
                sprintf(eval_buff, "0.0");
            } else if (evalue < 1.0e-99) {
                sprintf(eval_buff, "%2.0Le", evalue);
                eval_buff_ptr++; 	/* Knock off digit. */
            } else if (evalue < 0.0009) {
                sprintf(eval_buff, "%3.0Le", evalue);
            } else if (evalue < 0.1) {
                sprintf(eval_buff, "%4.3Lf", evalue);
            } else if (evalue < 1.0) {
                sprintf(eval_buff, "%3.2Lf", evalue);
            } else if (evalue < 10.0) {
                sprintf(eval_buff, "%2.1Lf", evalue);
            } else {
                sprintf(eval_buff, "%5.0Lf", evalue);
            }
	
            bit_score = txsp->bit_score;
            if (bit_score > 9999) {
                sprintf(bit_score_buff, "%4.3Le", bit_score);
            } else if (bit_score > 99.9) {
                sprintf(bit_score_buff, "%4.0ld", (long) bit_score);
            } else {
                sprintf(bit_score_buff, "%4.0ld", (long) bit_score); /* %4.0Lf is bad on 68K Mac, so cast to long */
            }
            
            if (options & TXALIGN_HTML) {
                if (gi != 0)
                    sprintf(id_buffer, "%ld", (long) gi);
                else
                    sprintf(id_buffer, "%s", txsp->buffer_id);
                bit_score_buff_ptr = bit_score_buff; 
                if (*bit_score_buff_ptr == ' ') {
                    bit_score_buff_ptr++;
                    sprintf(buffer1, " <a href = #%s>%s</a>", id_buffer, bit_score_buff_ptr);
                } else {
                    sprintf(buffer1, "<a href = #%s>%s</a>", id_buffer, bit_score_buff_ptr);
                }
            } else {
                sprintf(buffer1, "%s", bit_score_buff);
            }
            
            if (options & TXALIGN_SHOW_NO_OF_SEGS) 
                sprintf((ptr+pos), " %s  %s  %ld", buffer1, eval_buff_ptr, (long) txsp->number);
            else
                sprintf((ptr+pos), " %s  %s", buffer1, eval_buff_ptr);
        }
#else
        if (txsp->found_score) {
            evalue = txsp->evalue;
            eval_buff_ptr = eval_buff;
            if (evalue < 1.0e-180) {
                sprintf(eval_buff, "0.0");
            } else if (evalue < 1.0e-99) {
                sprintf(eval_buff, "%2.0le", evalue);
                eval_buff_ptr++; 	/* Knock off digit. */
            } else if (evalue < 0.0009) {
                sprintf(eval_buff, "%3.0le", evalue);
            } else if (evalue < 0.1) {
                sprintf(eval_buff, "%4.3lf", evalue);
            } else if (evalue < 1.0) {
                sprintf(eval_buff, "%3.2lf", evalue);
            } else if (evalue < 10.0) {
                sprintf(eval_buff, "%2.1lf", evalue);
            } else {
                sprintf(eval_buff, "%5.0lf", evalue);
            }
            
            bit_score = txsp->bit_score;
            if (bit_score > 9999) {
                sprintf(bit_score_buff, "%4.3le", bit_score);
            } else if (bit_score > 99.9) {
                sprintf(bit_score_buff, "%4.0ld", (long) bit_score);
            } else {
                sprintf(bit_score_buff, "%4.0lf", bit_score);
            }
            
            if (options & TXALIGN_HTML) {
                if (gi != 0)
                    sprintf(id_buffer, "%ld", (long) gi);
                else {
                    /*
                      sprintf(id_buffer, "%s", txsp->buffer_id);
                    */
                    MuskSeqIdWrite(txsp->id, id_buffer, BUFFER_LENGTH, PRINTID_TEXTID_ACCESSION, FALSE, FALSE);
                }
                bit_score_buff_ptr = bit_score_buff; 
                if (*bit_score_buff_ptr == ' ') {
                    bit_score_buff_ptr++;
                    sprintf(buffer1, " <a href = #%s>%s</a>", id_buffer, bit_score_buff_ptr);
                } else {
                    sprintf(buffer1, "<a href = #%s>%s</a>", id_buffer, bit_score_buff_ptr);
                }
            } else {
                sprintf(buffer1, "%s", bit_score_buff);
            }
            
            if (options & TXALIGN_SHOW_NO_OF_SEGS) 
                sprintf((ptr+pos), " %s  %s  %ld", buffer1, eval_buff_ptr, (long) txsp->number);
            else
                sprintf((ptr+pos), " %s  %s", buffer1, eval_buff_ptr);
        }
#endif
        if (countdescr == -1 || countdescr > 0)
            fprintf(outfp, "%s\n", ptr_start);
        
        txsp = txsp->next;
        if (countdescr >= 0)
            countdescr--;
    }
    
    if (options & TXALIGN_HTML) {
        ff_AddString("</PRE>");
        NewContLine();
    } else
       fprintf(outfp, "\n");
    
    /* blast_type (overwriting parameter) allocated before last while loop. */
    blast_type = (CharPtr) MemFree(blast_type);
    
    txsp = txsp_head;
    while (txsp) {
        txsp->title = (CharPtr) MemFree(txsp->title);
        txsp->buffer_id = (CharPtr) MemFree(txsp->buffer_id);
        txsp->id = SeqIdSetFree(txsp->id);
        txsp->segs_str = (CharPtr) MemFree(txsp->segs_str);
        txsp_var = txsp;
        txsp = txsp->next;
        MemFree(txsp_var);
    }
    
    return TRUE;
}

NLM_EXTERN Boolean LIBCALL
PrintDefLinesFromSeqAlign(SeqAlignPtr seqalign, Int4 line_length, FILE *outfp, Uint4 options, Int4 mode, Int2Ptr marks)
{
    Boolean	retval;

    retval = PrintDefLinesFromSeqAlignEx(seqalign, line_length, outfp, options, mode, marks, -1);

    return retval;
}

/* 
	Converts a number into a frame.
*/
static CharPtr 
NumToFrame(Int2 frame, CharPtr buffer)
{
	if (buffer)
	{
		if (frame > 0)
		{
			sprintf(buffer, "+%d", frame);
		}
		else
		{
			sprintf(buffer, "%d", frame);
		}
	}

	return buffer;
}

/* This function transfer SeqAlignPtr into AlignStatOptionPtr */

NLM_EXTERN Boolean FormatScoreFromSeqAlignEx(SeqAlignPtr sap, Uint4 option, FILE *fp, Int4Ptr PNTR matrix, Boolean follower, Boolean ooframe)
{
    AlignStatOptionPtr asop;
    Int4 empty_space, line_len;
    AlignSum as;
    SeqAlignPtr sap_tmp;

    asop = (AlignStatOptionPtr) MemNew(sizeof(AlignStatOption));
    MemSet(&as, 0, sizeof(AlignSum));

    empty_space = 12; line_len = 60; /* TO BE DETERMINED !!!! */

    asop->indent_len = (Int2) empty_space;
    asop->line_len = (Int2) (line_len + empty_space);
    asop->html_hot_link_relative = FALSE;
    
    if (option & TXALIGN_NO_ENTREZ)
        asop->no_entrez = TRUE;
    else
        asop->no_entrez = FALSE;
    
    if (option & TXALIGN_NO_DUMPGNL)
        asop->no_dumpgnl = TRUE;
    else
        asop->no_dumpgnl = FALSE;
    
    if (option & TXALIGN_HTML) {
        asop->html_hot_link = TRUE;
        if (option & TXALIGN_HTML_RELATIVE)
            asop->html_hot_link_relative = TRUE;
    } else {
        asop->html_hot_link = FALSE;
    }
    if (option & TXALIGN_SHOW_GI)
        asop->show_gi = TRUE;
    else
        asop->show_gi = FALSE;

    asop->fp = fp;
    asop->buf = NULL;
    asop->segs = NULL;
    as.matrix = matrix;

    as.master_sip = TxGetQueryIdFromSeqAlign(sap);
    as.target_sip = TxGetSubjectIdFromSeqAlign(sap);

    if((asop->bsp = BioseqLockById(as.target_sip)) == NULL) {
        Char tmp[128];
        SeqIdWrite(as.target_sip, tmp, PRINTID_FASTA_LONG, sizeof(tmp));
        ErrPostEx(SEV_ERROR, 0, 0, "Failure to get Bioseq for %s\n", tmp);
        return FALSE;
    }

    as.is_aa = (asop->bsp->mol == Seq_mol_aa);
    as.ooframe = ooframe;
    
    asop->sp = NULL;
    if(sap->segtype == SAS_DISC) {

        Int4 last_m_to = 0, last_t_to = 0;
        Int4 m_adj = 0, t_adj = 0;
        for(sap_tmp = (SeqAlignPtr)sap->segs; sap_tmp != NULL; 
            sap_tmp = sap_tmp->next) {
            
            /* We cannot find score this way .. :-) this fuction just
               calculates number of positives,identities etc. */
            
            find_score_in_align(sap_tmp, 1, &as);
            
            asop->gaps += as.gaps;
            asop->positive += as.positive;
            asop->identical += as.identical;
            asop->align_len += as.totlen;
            
            /* Adjustment for unaligned regions not counted in the
               function above */
            
            if(last_m_to != 0) {
                m_adj = as.master_from - last_m_to - 1;
            }
            
            asop->align_len += m_adj;
            asop->gaps += m_adj;
            
            last_m_to = as.master_to;            
        }
        asop->sp = sap->score;
    } else {
        asop->sp = find_score_in_align(sap, 1, &as);
        asop->gaps = as.gaps;
        asop->positive = as.positive;
        asop->identical = as.identical;
        asop->align_len = as.totlen;
    }

    asop->db_name = NULL;

    if (as.m_frame_set) {
        asop->m_frame = as.m_frame;
    } else {
        asop->m_frame = 255;
    }
    
    if (as.t_frame_set) {
        asop->t_frame = as.t_frame;
    } else {
        asop->t_frame = 255;
    }
    
    asop->m_strand = as.m_strand;
    asop->t_strand = as.t_strand;    

    /*    if(!ooframe) {
          asop->m_frame = 255;
          asop->t_frame = 255;
          } else {
          asop->m_frame = as.m_frame;
          asop->t_frame = as.t_frame;
          } */

    /* asop->m_strand = Seq_strand_unknown;
       asop->t_strand = Seq_strand_unknown; */

    asop->follower = follower;
    
    init_buff_ex(255);
    FormatScoreFunc(asop);
    free_buff();

    BioseqUnlock(asop->bsp);
    
    MemFree(asop);

    return TRUE;
}
/* This function transfer SeqAlignPtr into AlignStatOptionPtr */

NLM_EXTERN Boolean FormatScoreFromSeqAlign
(SeqAlignPtr sap, Uint4 option, FILE *fp,
 Int4Ptr PNTR matrix, Boolean follower)
{
    return FormatScoreFromSeqAlignEx(sap, option, fp, matrix, follower, FALSE);
}

static CharPtr FSFPrintOneDefline(AlignStatOptionPtr asop, Boolean is_na, 
                                  SeqIdPtr sip, CharPtr defline)
{
    Char HTML_database[32], HTML_dopt[16], id_buffer[BUFFER_LENGTH+1];
    Char buffer[BUFFER_LENGTH+1];
    SeqIdPtr bestid;
    Boolean make_link = FALSE, found_next_one, found_gnl_id;
    DbtagPtr db_tag;
    ObjectIdPtr oip;
    CharPtr ptr;
    Int4 gi, seqid_len = 0;

    /* Printing full label to the buffer */
     SeqIdWrite(sip, buffer, PRINTID_FASTA_LONG, BUFFER_LENGTH);

    if (asop->html_hot_link == TRUE)  {
        /* if (ISA_na(bsp->seq_data_type)) */
        
        if (is_na) {
            StringCpy(HTML_dopt, "GenBank");
            StringCpy(HTML_database, "Nucleotide");
        } else {
            StringCpy(HTML_dopt, "GenPept");
            StringCpy(HTML_database, "Protein");
        }

        bestid = SeqIdFindBest(sip, SEQID_GI);
        make_link = FALSE;
        gi = 0;
        if (bestid != NULL) {
            if (bestid->choice == SEQID_GI && asop->no_entrez == FALSE) {
                gi = bestid->data.intvalue;
                make_link = TRUE;
                sprintf(id_buffer, "%ld", (long) gi);
            } else {
                MuskSeqIdWrite(bestid, id_buffer, BUFFER_LENGTH, PRINTID_TEXTID_ACCESSION, FALSE, FALSE);
            }
            
            fprintf(asop->fp, "<a name = %s></a>", id_buffer);
            if (make_link) {
                fprintf(asop->fp, 
                        "<a href=\"%s?cmd=Retrieve&db=%s&list_uids=%08ld&dopt=%s\">", 
                        NEW_ENTREZ_HREF, HTML_database, (long) gi, HTML_dopt);
                
            } else if ((bestid->choice == SEQID_GENERAL && asop->no_dumpgnl == FALSE) || asop->no_entrez == TRUE) {
                if (bestid->choice != SEQID_GENERAL && bestid->choice != SEQID_OTHER)
                    { /* HACK, HACK, use SEQID_OTHER for Greg's page, even though GI is present. */
                        /* bsp is already present. */
                        bestid = SeqIdFindBest(sip, SEQID_OTHER);
                    }
                if (bestid->choice == SEQID_GENERAL) {
                    db_tag = (DbtagPtr) bestid->data.ptrvalue;
                    if(db_tag->db && StringCmp(db_tag->db, "THC") == 0) {
                        oip = db_tag->tag;
                        if(oip->id != 0) {
                            fprintf(asop->fp, "<a href=\"http://www.tigr.org/docs/tigr-scripts/hgi_scripts/thc_report.spl?est=THC%ld&report_type=n\">", (long) oip->id);
                            
                        }
                    } else {
                        /** * links to incomplete genomes */
                        make_dumpgnl_links(sip, asop->blast_type, asop->segs, asop->db_name, is_na, asop->fp, buffer, asop->no_entrez);
                    }
                } else {
                    make_dumpgnl_links(sip, asop->blast_type, asop->segs, asop->db_name, is_na, asop->fp, buffer, asop->no_entrez);
                }
                make_link = TRUE;
            }
        }
    }

    /* else {
       fprintf(asop->fp, ">");
       } */
        
    found_next_one = FALSE;
    if (asop->show_gi == FALSE) {
        if (StringNCmp(buffer, "gi|", 3) == 0) {
            ptr = &buffer[3];
            while (*ptr != NULLB && *ptr != ' ') { 
                /* Is there another ID beside the GI? */
                if (*ptr == '|') {
                    ptr++;
                    found_next_one = TRUE;
                    break;
                }
                ptr++;
            }
        }
    }
    if (found_next_one == FALSE) /* If TRUE, then ptr set above. */
        ptr = buffer;
    
    /* Remove local ID's. */
    if (StringNCmp(ptr, "lcl|", 4) == 0) {
        ptr += 4;
    }
    
    found_gnl_id = TRUE;
    /* Check for an ID of type general from BLAST */
    if (StringNCmp(buffer, "gnl|BL_ORD_ID", 13) != 0) {
        fprintf(asop->fp, "%s", ptr);
        seqid_len = StringLen(ptr);
        found_gnl_id = FALSE;
    } else {
        make_link = FALSE;
    }
    
    if (asop->html_hot_link == TRUE && make_link == TRUE) {
        fprintf(asop->fp, "</a> ");
    } else if (found_gnl_id == FALSE) {
        fprintf(asop->fp, " ");
    }
    
    /* Subtract 10 off the lines length as the ID is not printed
       with ffprint functions. */

    ff_StartPrint(0, asop->indent_len, 
                  (Int2)(asop->line_len+asop->indent_len-15), NULL); 
    ff_AddString(defline);
    ff_EndPrint();
    
    return NULL;
}

NLM_EXTERN int LIBCALLBACK FormatScoreFunc(AlignStatOptionPtr asop)

{
    BioseqPtr bsp;
    Boolean allocated, first;
    CharPtr defline, ptr, eval_buff_ptr, dline_buf, chptr, sptr;
    CharPtr new_defline;
    Char buf1[5], buf2[5];
    Char buffer[BUFFER_LENGTH+1], eval_buff[10], bit_score_buff[10];
    Char HTML_buffer[BUFFER_LENGTH+1], seqid_buf[128];
    Char id_buffer[BUFFER_LENGTH+1];
    Nlm_FloatHi bit_score, evalue; 
    Int4 percent_identical, percent_positive;
    Int4 number, score, gi, len, i;
    ObjectIdPtr obid;
    SeqIdPtr gilist, sip, new_sip, sip_tmp;
    ScorePtr	scrp, sp;
    Boolean splice_junction = FALSE;
    
    sp = asop->sp;
    bsp = asop->bsp;
    
    asn2ff_set_output(asop->fp, NULL);
    
    bit_score = 0.0;
    score = 0;
    evalue = 0.0;
    defline = NULL;
    *id_buffer = NULLB;

    if (bsp && asop->follower == FALSE) {
        /* Is the defline and sip allocated? */
        allocated = FALSE;
        gilist = ScorePtrUseThisGi(sp);
        if (gilist) {
            FilterTheDefline(bsp, gilist, buffer, BUFFER_LENGTH, &(defline));
            gilist = SeqIdSetFree(gilist);
            sip = SeqIdParse(buffer);
            allocated = TRUE;
        } else {
            SeqIdWrite(bsp->id, buffer, PRINTID_FASTA_LONG, BUFFER_LENGTH);
            sip = SeqIdSetDup(bsp->id);
            defline = StringSave(BioseqGetTitle(bsp));
        }
        
        /* Here we print all defline one by one */

        dline_buf = defline;
        chptr = defline;
        new_sip = NULL;
        new_defline = NULL;
        first = TRUE;
        while (TRUE) {
            if((chptr = StringChr(chptr, '>')) != NULL) {
                
                /* If ">" character exists in the defline - we have to check,
                   that this is start of new SeqId string */
                
                for (ptr = chptr+1, sptr = seqid_buf; 
                     *ptr != ' ' && *ptr != NULLB; 
                     ptr++, sptr++) {
                    *sptr = *ptr;
                }
                *sptr = NULLB;
                
                if((new_sip = SeqIdParse(seqid_buf)) == NULL) {
                    chptr++;
                    continue;
                }
                
                *chptr = NULLB;
                if(*ptr == ' ')
                    new_defline = ptr + 1;
                else 
                    new_defline = NULL;
            } else {
                new_sip = NULL;
                new_defline = NULL;
            }

            if(sip != NULL) {
                
                if(first) {
                    fprintf(asop->fp, ">");
                    first = FALSE;
                } else {
                    fprintf(asop->fp, " ");
                }
                len = StringLen(defline);

                /* Trimming tail white spaces if any */
                for(i = len; i > 0 && IS_WHITESP(defline[i-1]); i++)
                    defline[i-1] = NULLB;

                FSFPrintOneDefline(asop, ISA_na(bsp->mol), sip, defline);
                sip =SeqIdSetFree(sip);
            }

            if(new_sip != NULL && new_defline != NULL) {
                chptr = defline = new_defline;
                sip = new_sip;
            } else {
                break;
            }
        }

        fprintf(asop->fp, "          Length = %ld\n\n", 
                (long) BioseqGetLen(bsp));
        
        dline_buf = (CharPtr) MemFree(dline_buf);
    }

    if (asop->no_entrez == TRUE && 
        asop->html_hot_link == TRUE && bsp != NULL) {
        
        /* For Gregs and Human Genome stuff we will add links to every
           HSP */
        SeqIdPtr bestid;
        
        gi = 0;
        bestid = SeqIdFindBest(bsp->id, SEQID_GI);
        
        if (bestid != NULL) {
            if (bestid->choice == SEQID_GI) {
                gi = bestid->data.intvalue;
                sprintf(id_buffer, "%ld", (long) gi);
            }            
        }
    }
    
    number=1;
    for (scrp=sp; scrp; scrp = scrp->next) {
        obid = scrp->id;
        if(obid != NULL) {
            if (StringICmp(obid->str, "score") == 0) {
                score = scrp->value.intvalue;
                continue;
            }
            else if (StringICmp(obid->str, "e_value") == 0 || StringICmp(obid->str, "sum_e") == 0) {
                evalue = scrp->value.realvalue;
                continue;
            } else if (StringICmp(obid->str, "sum_n") == 0) {
                number = scrp->value.intvalue;
                continue;
            } else if (StringICmp(obid->str, "bit_score") == 0) {
                bit_score = scrp->value.realvalue;
                continue;
            } else if (StringICmp(obid->str, "splice_junction") == 0) {
               splice_junction = TRUE;
            }
        } else {
            if(scrp->choice == 1) {
                score = scrp->value.intvalue;
                continue;
            } else if(scrp->choice == 2) {
                bit_score = scrp->value.realvalue;
                continue;
            }
        }
    }
    
    ff_StartPrint(0, 0, (Int2)(asop->line_len+asop->indent_len), NULL);
    
    eval_buff_ptr = eval_buff;
#ifdef OS_MAC
    if (evalue < 1.0e-180) {
        sprintf(eval_buff, "0.0");
    } else if (evalue < 1.0e-99) {
        sprintf(eval_buff, "%2.0Le", evalue);
        eval_buff_ptr++;	/* Knock off digit. */
    } else if (evalue < 0.0009) {
        sprintf(eval_buff, "%3.0Le", evalue);
    } else if (evalue < 0.1) {
        sprintf(eval_buff, "%4.3Lf", evalue);
    } else if (evalue < 1.0) {
        sprintf(eval_buff, "%3.2Lf", evalue);
    } else if (evalue < 10.0) {
        sprintf(eval_buff, "%2.1Lf", evalue);
    } else {
        sprintf(eval_buff, "%5.0Lf", evalue);
    }
    
    if (bit_score > 9999) {
        sprintf(bit_score_buff, "%4.3Le", bit_score);
    } else if (bit_score > 99.9) {
        sprintf(bit_score_buff, "%4.0ld", (long) bit_score);
    } else {
        sprintf(bit_score_buff, "%4.0ld", (long) bit_score); /* %4.1Lf is bad on 68K Mac, so cast to long */
    }
#else
    if (evalue < 1.0e-180) {
        sprintf(eval_buff, "0.0");
    } else if (evalue < 1.0e-99) {
        sprintf(eval_buff, "%2.0le", evalue);
        eval_buff_ptr++;	/* Knock off digit. */
    } else if (evalue < 0.0009) {
        sprintf(eval_buff, "%3.0le", evalue);
    } else if (evalue < 0.1) {
        sprintf(eval_buff, "%4.3lf", evalue);
    } else if (evalue < 1.0) {
        sprintf(eval_buff, "%3.2lf", evalue);
    } else if (evalue < 10.0) {
        sprintf(eval_buff, "%2.1lf", evalue);
    } else {
        sprintf(eval_buff, "%5.0lf", evalue);
    }

    if (bit_score > 9999) {
        sprintf(bit_score_buff, "%4.3le", bit_score);
    } else if (bit_score > 99.9) {
        sprintf(bit_score_buff, "%4.0ld", (long) bit_score);
    } else {
        sprintf(bit_score_buff, "%4.1lf", bit_score);
    }
#endif

    if(asop->html_hot_link == TRUE && *id_buffer != NULLB) {
        sprintf(buffer, " <a name = %s_%ld></a>Score = %s bits (%ld), ", 
                id_buffer, (long) score, bit_score_buff, (long) score);
    } else {
        sprintf(buffer, " Score = %s bits (%ld), ", 
                bit_score_buff, (long) score);
    }
    ff_AddString(buffer);

    if (number == 1)
        sprintf(buffer, "Expect = %s", eval_buff_ptr);
    else if (!splice_junction)
        sprintf(buffer, "Expect(%ld) = %s", (long) number, eval_buff_ptr);
    else
        sprintf(buffer, "Expect(%ld+) = %s", (long) number, eval_buff_ptr);

    ff_AddString(buffer);
    NewContLine();
    if (asop->align_len > 0) {
        percent_identical = (100*asop->identical)/ (asop->align_len);
        percent_positive = (100*asop->positive)/ (asop->align_len);
        /* Don't show positives for blastn, which has these set to 255. */
        if (asop->m_frame == 255 && asop->t_frame == 255 &&
            asop->m_strand != Seq_strand_unknown && asop->t_strand != Seq_strand_unknown)
            sprintf(buffer, " Identities = %ld/%ld (%ld%%)", (long) asop->identical, (long) asop->align_len, (long) percent_identical);
        else
            sprintf(buffer, " Identities = %ld/%ld (%ld%%), Positives = %ld/%ld (%ld%%)", (long) asop->identical, (long) asop->align_len, (long) percent_identical, (long) (asop->positive+asop->identical), (long) asop->align_len, (long) (percent_identical+percent_positive));
        ff_AddString(buffer);
        if (asop->gaps > 0) {
            sprintf(buffer, ", Gaps = %ld/%ld (%ld%%)", (long) asop->gaps, 
                    (long) asop->align_len, 
                    (long) (100*asop->gaps)/(asop->align_len));
            ff_AddString(buffer);
        }
        NewContLine();

        /* for testing. */
        if (asop->m_frame != 255 || asop->t_frame != 255) {
            if (asop->m_frame != 255 && asop->t_frame != 255) {
                sprintf(buffer, " Frame = %s / %s", NumToFrame(asop->m_frame, buf1), NumToFrame(asop->t_frame, buf2));
            } else if (asop->m_frame != 255) {
                sprintf(buffer, " Frame = %s", NumToFrame(asop->m_frame, buf2));
            }
            else if (asop->t_frame != 255) {
                sprintf(buffer, " Frame = %s", NumToFrame(asop->t_frame, buf2));
            }
            ff_AddString(buffer);
            NewContLine();
        } else if (asop->m_strand != Seq_strand_unknown && asop->t_strand != Seq_strand_unknown) {
            if (asop->m_strand == Seq_strand_minus)
                sprintf(buffer, " Strand = Plus / Minus");
            else
                sprintf(buffer, " Strand = Plus / Plus");
            ff_AddString(buffer);
            NewContLine();
        }
        /* for testing. */
    }
    ff_EndPrint();
    
    return 0;
}

/*
*
*	determine the option for alignment based on the named tx_option
*
*/
NLM_EXTERN Uint4 GetTxAlignOptionValue (Uint1 tx_option, BoolPtr hide_feature, 
	BoolPtr print_score, BoolPtr split_display)
{
	Uint4 option;

	option = 0;
	*print_score = FALSE;
	*split_display = FALSE;
	switch (tx_option)
	{
		/*multiple pairwise alignment */
		case TEXT_MP:
		case TEXT_MP_MISMATCH:
			option |= TXALIGN_MASTER;
			option |= TXALIGN_SHOW_RULER;
			option |= TXALIGN_SHOW_STRAND;
			if(tx_option == TEXT_MP_MISMATCH)
				option |= TXALIGN_MISMATCH;
			break;
		/*FLAT multiple pairwise alignment*/
		case TEXT_MPFLAT:
		case TEXT_MPFLAT_MISMATCH:
			option |= TXALIGN_MASTER;
			option |= TXALIGN_FLAT_INS;
			option |= TXALIGN_END_NUM;
			option |= TXALIGN_COMPRESS;
			if(tx_option == TEXT_MPFLAT_MISMATCH)
				option |= TXALIGN_MISMATCH;
			*split_display = TRUE;
			break;
		case TEXT_BLAST:
			option |= TXALIGN_END_NUM;
			option |= TXALIGN_BLASTX_SPECIAL;
			option |= TXALIGN_MATRIX_VAL;
			option |= TXALIGN_SHOW_QS;
			*hide_feature = TRUE;
			*print_score = TRUE;
			*split_display = TRUE;
			break;
		default:
			option |= TXALIGN_MASTER;
			option |= TXALIGN_MISMATCH;
			option |= TXALIGN_SHOW_RULER;
			option |= TXALIGN_SHOW_STRAND;
			break;
	}
	if(*hide_feature)
		option |= TXALIGN_COMPRESS;
	return option;
}

Int4 OOFGetDNAStrand(StdSegPtr sseg)
{
    Int4 dna_strand;
    SeqIntPtr seq_int1;
    SeqLocPtr slp1;

    for(; sseg != NULL; sseg= sseg->next) {
        slp1 = sseg->loc;
        
        if(slp1->choice == SEQLOC_INT) {
            seq_int1 = (SeqIntPtr) slp1->data.ptrvalue;
            return seq_int1->strand;
        }
    }
    return Seq_strand_unknown;
}
static SetDNALineEnd(Int4 dna_index, Int4 dna_strand)
{
    Int4 dna_line_end;
    
    if(dna_strand != Seq_strand_minus)
        dna_line_end = dna_index == 0 ? 0 : dna_index -3;
    else
        dna_line_end = dna_index < 1 ? dna_index : dna_index -1;
    
    return dna_line_end;
}

static Int4 GetDigitsInINT(Int4 number)
{
    Int4 count;

    for(count = 1; number > 9; count++) 
        number = number/10;

    return count;
}

static Int4 GetMaxFROMDigits(StdSegPtr sseg)
{
    StdSegPtr ssp, ssp_last;
    Int4 master_from, target_from, master_to, target_to;
    Int4 max_number, count;

    master_from = SeqLocStart(sseg->loc);
    target_from = SeqLocStart(sseg->loc->next);

    for(ssp_last = ssp = sseg; ssp != NULL; ssp = ssp->next)
        ssp_last = ssp;

    master_to = SeqLocStop(ssp_last->loc);
    target_to = SeqLocStop(ssp_last->loc->next);

    max_number = MAX(MAX(master_from, master_to), 
                     MAX(target_from, target_to));

    count = GetDigitsInINT(max_number);

    return count;
}

#define WIDTH 60
static Boolean OOFShowSingleAlignment(SeqAlignPtr sap, ValNodePtr mask,
                                      Int4Ptr PNTR matrix, FILE *fp)
{
    StdSegPtr sseg;
    SeqIntPtr seq_int1, seq_int2;
    SeqLocPtr slp, slp1, slp2;
    SeqIdPtr sip1, sip2;
    SeqFeatPtr fake_cds;
    ByteStorePtr b_store = NULL;
    Char line1[128], line2[128], line3[128];
    Char tmpbuf[256];
    Int4 line_index, length_dna, length_pro, length;
    Int4 dna_index, pro_index, dna_line_start, pro_line_start;
    Int4 dna_line_end, pro_line_end, dna_to, dna_from;
    BioseqPtr bsp;
    SeqPortPtr spp;
    Int4 i, lines, k, shift_info = 0;
    Char  c1, c2, c3;
    Int4 dna_strand, max_digits, num_pad;

    if(sap == NULL || sap->segtype != 3) /* Should be StdSeg here! */
        return FALSE;
    
    line_index = 0;
    lines = 0;
    dna_index =0;
    pro_index = 0;
    pro_line_end = 0;
    dna_line_end = 0;
    
    dna_strand = OOFGetDNAStrand((StdSegPtr) sap->segs);

    /* Needed for printing nice alignment with normal spacing */
    max_digits = GetMaxFROMDigits((StdSegPtr) sap->segs);

    for(sseg = (StdSegPtr) sap->segs; sseg != NULL; sseg= sseg->next) {
        
        /* Now starting new alignment region */
        
        length_dna = 0;
        length_pro = 0;
        b_store = NULL;
        
        slp1 = sseg->loc;
        
        if(slp1->choice == SEQLOC_INT) 
            seq_int1 = (SeqIntPtr) slp1->data.ptrvalue;
        else if (slp1->choice == SEQLOC_EMPTY)
            seq_int1 = NULL;
        else
            return FALSE;       /* Invalid SeqLoc */
        
        slp2 = sseg->loc->next;
        
        if(slp2->choice == SEQLOC_INT)
            seq_int2 = (SeqIntPtr) slp2->data.ptrvalue;
        else if (slp2->choice == SEQLOC_EMPTY)
            seq_int2 = NULL;
        else
            return FALSE;       /* Invalid SeqLoc */

        /* Ignore double gap */
        if(seq_int1 == NULL && seq_int2 == NULL)
            continue;
        
        sip1 = sseg->ids;       /* DNA */
        sip2 = sseg->ids->next; /* Protein */

        /* printf("shift_info = %d\n", shift_info); */

        if(shift_info%3)
            dna_index -= (3 - shift_info); /* adjustment for frameshift */

        switch(shift_info) {
        case 1:
            line1[line_index] = '\\';
            line2[line_index] = ' ';
            line3[line_index] = ' ';
            line_index++;

            if(line_index == WIDTH) {
                dna_line_end = SetDNALineEnd(dna_index, dna_strand);
                pro_line_end = pro_index;
            }
        case 2:
            line1[line_index] = '\\';
            line2[line_index] = ' ';
            line3[line_index] = ' ';
            line_index++;
            
            if(line_index == WIDTH) {
                dna_line_end = SetDNALineEnd(dna_index, dna_strand);
                pro_line_end = pro_index;
            }
           break;
        case 5:
            line1[line_index] = '/';
            line2[line_index] = ' ';
            line3[line_index] = ' ';
            line_index++;

            if(line_index == WIDTH) {
                dna_line_end = SetDNALineEnd(dna_index, dna_strand);
                pro_line_end = pro_index;
            }

        case 4:
            line1[line_index] = '/';
            line2[line_index] = ' ';
            line3[line_index] = ' ';
            line_index++;

            if(line_index == WIDTH) {
                dna_line_end = SetDNALineEnd(dna_index, dna_strand);
                pro_line_end = pro_index;
            }

            break;
        case 0:
        default:
            break;
        }
  
        /* Looking if any frame shift is followed next */
        if(seq_int1 != NULL && seq_int2 != NULL) {            
            shift_info = (seq_int1->to - seq_int1->from + 1) -
                (seq_int2->to - seq_int2->from)*3;
        } else if(seq_int1 != NULL) {
            shift_info = (seq_int1->to - seq_int1->from + 1)%3 + 3;
        } else {
            shift_info = 0;
        }

        if(seq_int1 != NULL) {

            if(dna_strand != Seq_strand_minus)
                dna_index = seq_int1->from;
            else 
                dna_index = seq_int1->to;
            
            length_dna = (seq_int1->to - seq_int1->from + 1)/3; 
        }        
        
        if(seq_int2 != NULL) {
            pro_index = seq_int2->from;
            length_pro = seq_int2->to - seq_int2->from + 1;
        }

        if(line_index == 0) {
            dna_line_start = dna_index; 
            pro_line_start = pro_index + 1;
        }

        if (dna_line_start == 0)
            dna_line_start = dna_index; 
        
        if(pro_line_start == 0)
            pro_line_start = pro_index + 1;

        if(seq_int1 != NULL) {

            /* if(length_dna == 0) insertion
               continue; */
            
            /* Byte store for DNA */
            bsp = BioseqLockById(sip1);


            dna_from = seq_int1->from;
            dna_to = seq_int1->to;
            
            if(0 < shift_info && shift_info < 3) {                
                if(dna_strand != Seq_strand_minus)
                    dna_to = seq_int1->to + 3 - shift_info;
                else 
                    dna_from = seq_int1->from - 3 + shift_info;
            }

            if(dna_from >= dna_to) {
                BioseqUnlock(bsp);
                continue;
            }

            fake_cds = make_fake_cds(bsp, dna_from, dna_to, 
                                     seq_int1->strand);
            BioseqUnlock(bsp);
            
            b_store = ProteinFromCdRegionEx(fake_cds, TRUE, FALSE);
            SeqFeatFree(fake_cds);

            if(b_store == NULL) {
                return FALSE;
            }

            BSSeek(b_store, 0, SEEK_SET);

            /* length_dna = BSLen(b_store); */
        }

        if(seq_int2 != NULL) {
            /* Seq port for protein */
            bsp = BioseqLockById(sip2);
            spp = SeqPortNew(bsp, seq_int2->from, 
                             seq_int2->to, 0, Seq_code_ncbieaa);
            BioseqUnlock(bsp);
        } else {
            spp = NULL;
        }
        
        if(length_dna == 0) length_dna = length_pro;
        if(length_pro == 0) length_pro = length_dna;

        length = MAX(length_pro, length_dna);
        /* length = MIN(length_pro, length_dna); */
        
        /* printf("length = %d\n", length); */
        for(i = 0; i < length; i++) {

            if(seq_int1 != NULL) {
                
                /* This line should be checked for correctness */
                if((line1[line_index] = BSGetByte(b_store)) == EOF)
                    line1[line_index] = '?';
                
                if(dna_strand != Seq_strand_minus)
                    dna_index += 3;
                else
                    dna_index -= 3;
                
            } else {
                line1[line_index] = '-';
            }
            
            if(seq_int2 != NULL) {
                line2[line_index] = SeqPortGetResidue(spp);
                pro_index++;
            } else {
                line2[line_index] = '-';
            }

            if(line1[line_index] == line2[line_index])
                line3[line_index] = line1[line_index];
            else if(matrix[line1[line_index]][line2[line_index]] > 0)
                line3[line_index] = '+';
            else
                line3[line_index] = ' ';

            line_index++;

            if(line_index == WIDTH) {
                dna_line_end = SetDNALineEnd(dna_index, dna_strand);
                pro_line_end = pro_index;
            }
            
            if(line_index > WIDTH) { /* Printout */
                
                line1[line_index] = line2[line_index] = line3[line_index] = '\0';
#ifdef SHOW_RULER
                fprintf(fp, "%5d",
 WIDTH*lines++);
                
                for (k = 10; k <= WIDTH; k+=10) 
                    fprintf(fp, "    .    :");
                if (k-5 < WIDTH) fprintf(fp, "    ."); 
#endif
                
                c1 = line1[WIDTH]; c2 = line2[WIDTH]; c3 = line3[WIDTH];
                line1[WIDTH] = line2[WIDTH] = line3[WIDTH] = '\0';

                /* ------- Printout of the alignment ------------- */ 

                fprintf(fp, "Query: %d", dna_line_start+1);
                
                num_pad = max_digits - GetDigitsInINT(dna_line_start+1) + 1;

                for(k=0; k < num_pad; k++)
                    fprintf(fp, " ");

                fprintf(fp, "%s %d\n", line1, dna_line_end+3);
                
                num_pad = 8 + max_digits;
                
                for(k=0; k < num_pad; k++)
                    fprintf(fp, " ");

                fprintf(fp, "%s\nSbjct: %d", line3, pro_line_start);

                num_pad = max_digits - GetDigitsInINT(pro_line_start) + 1;

                for(k=0; k < num_pad; k++)
                    fprintf(fp, " ");

                fprintf(fp, "%s %d\n\n", line2, pro_line_end);

                /* --------------------------------------------------- */

                if(dna_line_end != 0) {

                    if(dna_strand != Seq_strand_minus)
                        dna_line_start = dna_line_end+3; /*takes 3 bases*/
                    else
                        dna_line_start = dna_line_end+1; /*takes 3 bases*/
                }
                if(pro_line_end != 0) 
                    pro_line_start = pro_line_end+1;
                
                line1[WIDTH] = c1; line2[WIDTH] = c2; line3[WIDTH] = c3;
                strcpy(line1, &line1[WIDTH]);
                strcpy(line2, &line2[WIDTH]);
                strcpy(line3, &line3[WIDTH]);
                line_index = line_index - WIDTH;
            }
        }

        SeqPortFree(spp);       /* Protein SeqPort */
        BSFree(b_store);        /* DNA Byte store  */
    }
    
    /* Printing out remaining tail ... if any */
    line1[line_index] = line2[line_index] = line3[line_index] = '\0';

#ifdef SHOW_RULER
    fprintf(fp, "%5d", WIDTH*lines);

    for (k = 10; k < line_index; k+=10) 
        fprintf(fp, "    .    :");
    
    if (k-5 < line_index) fprintf(fp, "    .");
#endif

    dna_line_end = SetDNALineEnd(dna_index, dna_strand);
    pro_line_end = pro_index;


    /* ------- Printout of the alignment remainder ------- */ 
    
    fprintf(fp, "Query: %d", dna_line_start+1);
    
    num_pad = max_digits - GetDigitsInINT(dna_line_start+1) + 1;
    
    for(k=0; k < num_pad; k++)
        fprintf(fp, " ");
    
    fprintf(fp, "%s %d\n", line1, dna_line_end+3);
    
    num_pad = 8 + max_digits;
    
    for(k=0; k < num_pad; k++)
        fprintf(fp, " ");
    
    fprintf(fp, "%s\nSbjct: %d", line3, pro_line_start);
    
    num_pad = max_digits - GetDigitsInINT(pro_line_start) + 1;
    
    for(k=0; k < num_pad; k++)
        fprintf(fp, " ");
    
    fprintf(fp, "%s %d\n\n\n", line2, pro_line_end);
    
    /* --------------------------------------------------- */
    
    /*    fprintf(fp, "\nQuery: %-5d %s %-5d\n             "
          "%s\nSbjct: %-5d %s %-5d\n\n", 
          dna_line_start+1, line1, dna_line_end+3, line3, 
          pro_line_start, line2, pro_line_end); */
    
    return TRUE;
}

/*******************************************************************************

  Function : OOFShowBlastAlignment();
  
  Purpose : function to display a BLAST output with Out-of-Frame
            information
  
  Parameters : 	sap; seqalign
                mask; list of masked regions in the query
                fp; output file;
                tx_option; some display options
				
  Return value : FALSE if failure

*******************************************************************************/
NLM_EXTERN Boolean OOFShowBlastAlignment(SeqAlignPtr sap, ValNodePtr mask,
                                         FILE *fp, Uint4 tx_option, 
                                         Int4Ptr PNTR matrix)
{
    SeqAlignPtr     sap4;
    SeqIdPtr        new_id = NULL, old_id = NULL;    
    Uint4           option,i;
    Boolean         bRet, follower= FALSE, matrix_loaded = FALSE;
    
    if(sap == NULL || fp == NULL) 
        return FALSE;
    
    bRet = TRUE;

    /* get the matrix */
    
    if(matrix == NULL) { 
        if((matrix = load_default_matrix()) == NULL)
            return FALSE;
        matrix_loaded = TRUE;
    }

    for(sap4 = sap; sap4 != NULL; sap4 = sap4->next) {
        
        /* Attempt to print score for the alignment */
        new_id = TxGetSubjectIdFromSeqAlign(sap4);
        if(old_id != NULL) {
            if(SeqIdMatch(new_id, old_id))
                follower = TRUE;
        }
        
        old_id = new_id;
        if(!FormatScoreFromSeqAlignEx(sap4, tx_option, fp, matrix, 
                                      follower, TRUE)){
            bRet=FALSE;
            break;
        }
        
        follower = FALSE;
        
        /*display a SeqAlign*/
        if (!OOFShowSingleAlignment(sap4, mask, matrix, fp)) {
            bRet=FALSE;
            break;
        }
    }
    
    if (matrix_loaded){
        for(i = 0; i<TX_MATRIX_SIZE; ++i)
            MemFree(matrix[i]);
        MemFree(matrix);
    } 
    
    return(bRet);
	
}

NLM_EXTERN void OOFDisplayTraceBack1(Int4Ptr a, CharPtr dna, 
                                     CharPtr pro, Int4 ld, Int4 lp, 
                                     Int4 q_start, Int4 p_start)
{
    int len = 0, i, j, x, y, lines, k;
    static char line1[100], line2[100], line3[100],
        tmp[10] = "         ", *st;
    char *dna1, c1, c2, c3;

    dna1 = Malloc(ld+2);
    MemCpy(dna1+1, dna, ld);
    dna1[0] = ' '; dna1[1] = ' ';
    
    line1[0] = line2[0] = line3[0] = '\0'; x= q_start; y = p_start;
    printf("dna=%d pro=%d\n", y, x);
    
    for (len = 0, j = 0, lines = 0; x < lp && y < ld; j++) {
        i = a[j];
        switch(i) {
        case 0: 
            line1[len] = '-';
            line3[len] = ' ';
            line2[len++] = pro[x++];
            break;
        case 1:
        case 5:
            if (i == 1) line1[len]  = '\\';
            else line1[len] = '/';
            line2[len] = line3[len] = ' ';
            len++;
        case 2:
        case 4:	  
            if (i < 3) line1[len]  = '\\';
            else line1[len] = '/';
            line2[len] = line3[len] = ' ';
            len++;
        case 3:
            line1[len] = dna1[y+i-2]; y+= i;
            line2[len] = pro[x++];
            if (line1[len] == line2[len]) line3[len++] = '|';
            else line3[len++] = ' ';
            break;
        case 6:
            line1[len] = dna1[y+1]; y+= 3;
            line2[len] = '-';
            line3[len++] = ' ';
        }
        if (len >= WIDTH) {
            line1[len] = line2[len] = line3[len] = '\0';
            printf("\n%5d", WIDTH*lines++);
            for (k = 10; k <= WIDTH; k+=10) 
                printf("    .    :");
            if (k-5 < WIDTH) printf("    .");
            c1 = line1[WIDTH]; c2 = line2[WIDTH]; c3 = line3[WIDTH];
            line1[WIDTH] = line2[WIDTH] = line3[WIDTH] = '\0';
            printf("\n     %s\n     %s\n     %s\n", line1, line3, line2);
            line1[WIDTH] = c1; line2[WIDTH] = c2; line3[WIDTH] = c3;
            strcpy(line1, &line1[WIDTH]);
            strcpy(line2, &line2[WIDTH]);
            strcpy(line3, &line3[WIDTH]);
            len = len - WIDTH;
        }
    }
    printf("\n%5d", WIDTH*lines);
    line1[len] = line2[len] = line3[len] = '\0';
    for (k = 10; k < len; k+=10) 
        printf("    .    :");
    if (k-5 < len) printf("    .");
    printf("\n     %s\n     %s\n     %s\n", line1, line3, line2);

    MemFree(dna1);

    return;
}
NLM_EXTERN void OOFDisplayTraceBack2(Int4Ptr a, CharPtr dna, CharPtr pro, 
                                     Int4 ld, Int4 lp, 
                                     Int4 q_start, Int4 p_start)
{
    int len = 0, i, j, x, y, lines, k;
    static char line1[100], line2[100], line3[100],
        tmp[10] = "         ", *st;
    char *dna1, c1, c2, c3;
    
    dna1 = Malloc(ld+2);
    printf("%d %d\n", q_start, p_start);
    
    dna1 = Malloc(ld+2);
    MemCpy(dna1+1, dna, ld);
    dna1[0] = ' '; dna1[1] = ' ';
    
    line1[0] = line2[0] = line3[0] = '\0'; 
    x= q_start; 
    y = p_start;
    
    for (len = 0, j = 0, lines = 0; x < lp && y < ld; j++) {
        i = a[j];
        /*printf("%d %d %d\n", i, len, b->j);*/
        if (i > 0 && i < 6) {
            if (i == 1) {
                tmp[0] = pro[x++];
                len--;
                y--;
                i++;
            } else tmp[i-2] = pro[x++];
        }
        if (i == 6) {
            i = 3; tmp[0] = tmp[1] = tmp[2] = '-';
            if (a[j+1] == 2) tmp[2] = ' ';
        }
        if (i > 0) {
            strncpy(&line1[len], &dna1[y], i); y+=i;
        } else {line1[len] = '-'; i = 1; tmp[0] = pro[x++];}
        strncpy(&line2[len], tmp, i);
        for (k = 0; k < i; k++) {
            if (tmp[k] != ' ' && tmp[k] != '-') {
                if (k >= 2) tmp[k] = '\\';
                else if (k == 1) tmp[k] = '|';
                else tmp[k] = '/';
            } else tmp[k] = ' ';
        }
        if (i == 1) tmp[0] = ' ';
        strncpy(&line3[len], tmp, i); 
        tmp[0] = tmp[1] =  tmp[2] = ' ';
        len += i;
        line1[len] = line2[len] =line3[len]  = '\0'; 
        if (len >= WIDTH) {
            printf("\n%5d", WIDTH*lines++);
            for (k = 10; k <= WIDTH; k+=10) 
                printf("    .    :");
            if (k-5 < WIDTH) printf("    .");
            c1 = line1[WIDTH]; c2 = line2[WIDTH]; c3 = line3[WIDTH];
            line1[WIDTH] = line2[WIDTH] = line3[WIDTH] = '\0';
            printf("\n     %s\n     %s\n     %s\n", line1, line3, line2);
            line1[WIDTH] = c1; line2[WIDTH] = c2; line3[WIDTH] = c3;
            strcpy(line1, &line1[WIDTH]);
            strcpy(line2, &line2[WIDTH]);
            strcpy(line3, &line3[WIDTH]);
            len = len - WIDTH;
        }
    }
    printf("\n%5d", WIDTH*lines);
    for (k = 10; k < len; k+=10) 
        printf("    .    :");
    if (k-5 < len) printf("    .");
    printf("\n     %s\n     %s\n     %s\n", line1, line3, line2);
}
