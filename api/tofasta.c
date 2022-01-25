/*  tofasta.c
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
* File Name:  tofasta.c
*
* Author:  James Ostell
*   
* Version Creation Date: 7/12/91
*
* $Revision: 6.25 $
*
* File Description:  various sequence objects to fasta output
*
* Modifications:  
* --------------------------------------------------------------------------
* Date	   Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
* $Log: tofasta.c,v $
* Revision 6.25  1998/11/03 21:43:57  kans
* call SeqMapTableConvert to map asterisk to appropriate alphabet - BLAST is using ncbistdaa, not ncbieaa or iupacaa
*
* Revision 6.24  1998/10/26 22:35:53  kans
* filter out terminal asterisk in protein sequence
*
* Revision 6.23  1998/10/06 21:04:43  tatiana
* a bug fixed in CreateDefline()
*
* Revision 6.22  1998/09/04 13:09:22  volodya
* fix the memory bug and the memory leaks in CreateDefLine()
*
* Revision 6.21  1998/07/10 12:47:56  egorov
* Sequence code bug fixed
*
* Revision 6.20  1998/07/07 18:43:06  egorov
* Revert the BioseqToFastaX() to use FASTA_ID as it was before
*
* Revision 6.19  1998/03/19 21:14:29  kans
* defline says map instead of from for map source qualifier
*
* Revision 6.18  1998/03/12 17:40:52  kans
* added GSS to defline handling
*
* Revision 6.17  1998/02/25 20:06:30  kans
* default deflines for est, sts records
*
* Revision 6.16  1998/02/23 16:51:24  egorov
* Changes to make the tofasta.c independent on readdb.h
*
* Revision 6.15  1998/02/12 21:44:57  madden
* Set code before callto BioseqToFastaX
*
* Revision 6.14  1998/02/11 20:36:40  tatiana
* ABW error fixed
*
* Revision 6.13  1998/02/11 19:05:35  madden
* Changes to allow formatdb to take ASN.1 as input
*
* Revision 6.12  1998/02/09 18:54:07  kans
* MakeSeqID only in sqnutils, makes oid->id if numeric
*
* Revision 6.11  1998/01/27 20:28:08  madden
* Added BioseqRawToFastaExtra with line_length arg
*
* Revision 6.10  1997/12/08 18:15:19  ostell
* Hsiu-chuan sets tech value back to zero for repeated loading in HSP
*
* Revision 6.9  1997/12/02 19:04:45  shavirin
* Fixed typecast warnings and removed unused functions
*
* Revision 6.8  1997/11/28 15:55:47  shavirin
* Fixed some warnings. Removed redundant "*" and added function return type
*
* Revision 6.7  1997/11/14 15:49:47  tatiana
* added pdbip protection in CreateDefline()
*
* Revision 6.6  1997/11/13 19:39:59  tatiana
* *** empty log message ***
*
* Revision 6.5  1997/11/12 21:24:07  tatiana
* PDB special case in CreateDefline
*
* Revision 6.4  1997/11/05 19:25:23  kans
* special_symbol returned for &{}[]
*
* Revision 6.3  1997/11/03 20:46:54  shavirin
* Removed memory leak
*
* Revision 6.2  1997/10/22 16:43:06  shavirin
* Added new external functions FastaReadSequence() and FastaReadSequenceMem()
*
* Revision 6.1  1997/08/27 14:35:02  kans
* fixed memory leak in FastaReadSequenceInternal where an error message buffer was allocated but no bad characters were found
*
* Revision 6.0  1997/08/25 18:07:44  madden
* Revision changed to 6.0
*
* Revision 5.46  1997/06/20 19:13:34  kans
* treat - as illegal character for errmsg
*
* Revision 5.45  1997/06/19 18:39:17  vakatov
* [WIN32,MSVC++]  Adopted for the "NCBIOBJ.LIB" DLL'ization
*
* Revision 5.44  1997/05/28 18:38:17  tatiana
* a typo fixed in CreateDefLine
*
* Revision 5.43  1997/05/21 17:42:09  tatiana
* fixed uninitialized iip in UseOrgMods
*
* Revision 5.42  1997/05/21 17:22:00  tatiana
* CreateDefLine changed to take info from BioSource.mod
*
* Revision 5.41  1997/04/17 17:19:06  shavirin
* Removed truncation of one character in sequence without '>' in
* the beginning.
* FastaToSeqEntry()
*
 * Revision 5.40  1997/04/10  20:03:47  shavirin
 * Fixed bug in function FastaToSeqBuff()
 *
 *
 * Revision 5.39  1997/04/01  21:36:49  shavirin
 * Function FastaToSeqEntryInternal() was disintegrated into 2
 * functions: new function FastaReadSequenceInternal() was added
 * for this reason. This function may use used to create sequence
 * reding functions without creation SeqEntryPtr.
 *
 * Revision 5.38  1997/03/28  21:41:40  tatiana
 * turn off HTG phase 3 in CreateDefLine()
 *
 * Revision 5.37  1997/03/13  15:21:05  tatiana
 * CreateFeLine change to protect from 'HTG phase' duplication
 *
 * Revision 5.36  1997/02/20  15:07:54  kans
 * added call to SeqMgrAddToBioseqIndex, removed obsolete FastaIntSec
 * and GetLineFromBuff functions, added kludge to convert U to T for
 * nucleotides prior to SeqMapTableConvert lookup
 *
 * Revision 5.35  1997/02/11  01:14:10  kans
 * CountGapsInDeltaSeq takes 0, not NULL, as last parameter
 *
 * Revision 5.34  1997/02/10  22:08:20  tatiana
 * a bug fixed in CreateDefLine()
 *
 * Revision 5.32  1997/01/21  17:03:47  tatiana
 * added a check for SEQUENCING IN PROGRESS already existing in title
 * in CreateDefLine()
 *
 * Revision 5.31  1997/01/15  17:35:06  tatiana
 * a bug fixed in CreateDefLine()
 *
 * Revision 5.30  1997/01/13  23:20:54  tatiana
 * added SEQUENCING IN PROGRESS for htg phase 1 and 2 in CreateDefLine()
 *
 * Revision 5.28  1997/01/08  22:49:50  tatiana
 * buf and buflen arguments added to CountGapsInDeltaSeq()
 *
 * Revision 5.27  1996/11/25  18:59:52  shavirin
 * Reading from file in FastaToSeqEntryInternal() changed from FileRead
 * to getc back (from revision 5.16 due to Windows specific problems
 *
 * Revision 5.26  1996/11/19  22:45:09  tatiana
 * CreateDefLine changed for HTG
 *
 * Revision 5.25  1996/11/18  20:36:39  shavirin
 * Default coding for protein is set to ncbieaa in
 * FastaToSeqEntryInternal()
 *
 * Revision 5.24  1996/11/15  17:44:51  shavirin
 * Removed difference in upper and lower bad characters
 *
 * Revision 5.23  1996/11/04  17:59:33  tatiana
 * memory leaks cleaning in CreateDefline
 *
 * Revision 5.22  1996/10/30  16:52:07  kans
 * SeqIdDup now calls SeqIdFindBest first.  Also, removed unused variable.
 *
 * Revision 5.21  1996/10/28  21:51:47  tatiana
 * a bug fixed in FastaToSeqEntryInternal()
 *
 * Revision 5.20  1996/10/25  21:48:14  shavirin
 * Windows NT specific change - in FastaToSeqEntryInternal()
 * changed \n \r check to opposite id defline processing.
 *
 * Revision 5.19  1996/10/22  16:01:25  shavirin
 * Added support to disable sequence printing with new
 * parameter no_sequence and removed unused variables and functions
 *
 * Revision 5.18  1996/10/21  21:36:28  shavirin
 * Added new functions : SeqEntrysToDefline() and static DeflineFunc()
 *
 * Revision 5.17  1996/10/18  18:54:30  shavirin
 * Changed getc() macros to FileRead() function,
 * exit(1) to return NULL end ErrLogPrintf to ErrPostEx()
 *
 * Revision 5.16  1996/10/16  17:41:26  zjing
 * replace ch in getc from Char to Int1
 *
 * Revision 5.15  1996/10/11  20:30:09  shavirin
 * *** empty log message ***
 *
 * Revision 5.14  1996/10/11  18:47:01  shavirin
 * In FastaToSeqEntryInternal removed bad characters processing if
 * we do not require it
 *
 * Revision 5.13  1996/10/09  00:06:39  shavirin
 * Removed few unused variables and function definitions
 *
 * Revision 5.12  1996/10/08  22:25:52  shavirin
 * Rewritten engine for FastaToSeqEntry functions set.
 * Removed hardcoded limit of 255 characters for defline.
 *
 * Revision 5.11  1996/09/27  20:00:32  ostell
 * made FastaToSeqEntry call BioseqPack instead of BioseqConvert(ncbi4na)
 *
 * Revision 5.10  1996/08/30  14:57:11  tatiana
 * memory leak cleaned up in GatherDescrOnBioseq()
 *
 * Revision 5.9  1996/08/26  20:45:43  kans
 * FindProtDefLine length was uninitialized, also now protects against
 * prp->name with NULL ->data.ptrvalue.  (The latter is fixed in cdrgn.c.)
 *
 * Revision 5.8  1996/08/19  18:02:27  tatiana
 * a bug fixed in FindProtDefLine
 *
 * Revision 5.6  1996/08/15  18:15:52  tatiana
 * CreateDefLine() added
 *
 * Revision 5.5  1996/07/29  15:20:12  kans
 * converted errpost to errpostex
 *
 * Revision 5.4  1996/07/25  02:32:26  ostell
 * changed fastadefline() for elta seqs
 *
 * Revision 5.3  1996/07/11  16:25:31  kans
 * FastaToSeqEntry now calls MakeNewProteinSeqId to make unique ID
 *
 * Revision 5.2  1996/07/05  15:57:07  shavirin
 * Added support to parce buffers with multiple fasta entries in FastaToSeqBuff()
 *
 * Revision 5.1  1996/06/15  17:29:44  ostell
 * changed MyFsa structure by adding do_virtual and tech fields
 * added value of 3 for group_segs
 * addes support of tech to FastaDefLine()
 *
 * Revision 5.0  1996/05/28  13:23:23  ostell
 * Set to revision 5.0
 *
 * Revision 4.14  1996/05/07  21:25:17  madden
 * GetLineFromBuff has new check for length of buffer.
 *
 * Revision 4.13  1996/05/01  21:41:47  vakatov
 * GetLineFromBuff() function rewritten so that it can process DOS end-of-line
 * sequence
 *
 * Revision 4.12  1996/04/30  21:23:41  shavirin
 * Fixed bug for Macintosh in FastaToSeqBuff() function
 *
 * Revision 4.10  1996/04/03  15:15:06  kans
 * added parseSeqId to "Ex" functions
 *
 * Revision 4.9  1996/03/27  04:10:18  shavirin
 * Fixed problem with one-line sequence in
 * the function FastaToSeqBuff()
 *
 * Revision 4.8  1996/03/13  19:42:50  shavirin
 * Added new functions FastaToSeqBuff() and FastaToSeqBuffEx() that allow to transfer
 * FASTA entries into SeqEntry structure.
 * Redefined functions FastaToSeqEntry() and FastaToSeqEntryEx() using common engine -
 * internal function FastaIntSeq.
 * Added auxiliary function CharPtr GetLineFromBuff()
 *
 * Revision 4.7  1996/02/23  15:44:51  kans
 * FastaToSeqEntry now allows selenocysteine (U)
 *
 * Revision 4.6  1996/02/19  20:29:08  kans
 * FastaToSeqEntryEx now reports counts of individual illegal characters
 *
 * Revision 4.5  1996/02/17  00:51:54  kans
 * new FastaToSeqEntryEx hidden function
 *
 * Revision 4.4  1996/01/17  18:32:34  kans
 * added TrimSpacesAroundString, called by FastaToSeqEntry for each line
 *
 * Revision 4.3  1995/10/13  21:36:59  madden
 * ID's that contain a "|", but are not legal, are put in a local
 * ID in MakeSeqID.
 *
 * Revision 4.2  1995/09/29  15:03:06  ostell
 * changed FastaToSeqentry() so that if seqid fails to parse, makes a type local
 *
 * Revision 4.1  1995/09/26  12:52:57  ostell
 * removed extra newline from PATENT def line
 *
 * Revision 4.0  1995/07/26  13:49:01  ostell
 * force revision to 4.0
 *
 * Revision 2.26  1995/07/20  19:42:53  tatiana
 * change SeqIdPrint to SeqIdWrite
 *
 * Revision 2.25  1995/06/07  18:52:39  kans
 * kludge to allow - in nucleotide sequence, convert to Seq_code_ncbi4na)
 *
 * Revision 2.24  1995/06/07  16:59:41  kans
 * allow - and * in proteins, use Seq_code_ncbieaa
 *
 * Revision 2.23  1995/05/26  20:34:57  kans
 * gap of unknown length now makes virtual sequence of length -1
 *
 * Revision 2.22  1995/05/09  18:43:09  ostell
 * added support for (accession) on GenPept deflines and [organism] on
 * GenPept and PRF deflines
 *
*
* ==========================================================================
*/

#include <tofasta.h>
#include <gather.h>
#include <sqnutils.h>  /* MakeSeqID */

static Uint1 na_order[NUM_SEQID] = {   /* order of nucleic acid deflines */
 	255, /* 0 = not set */
	230, /* 1 = local Object-id */
	30,  /* 2 = gibbsq */
	30,  /* 3 = gibbmt */
	255, /* 4 = giim Giimport-id */
	20, /* 5 = genbank */
	20, /* 6 = embl */
	255, /* 7 = pir */
	255, /* 8 = swissprot */
	40,  /* 9 = patent */
	50, /* 10 = other TextSeqId */
	50, /* 11 = general Dbtag */
	120,  /* 12 = gi */
	20, /* 13 = ddbj */
	255, /* 14 = prf */
	30  /* 15 = pdb */
    };

static Uint1 aa_order[NUM_SEQID] = {   /* order of nucleic acid deflines */
 	255, /* 0 = not set */
	230, /* 1 = local Object-id */
	40,  /* 2 = gibbsq */
	40,  /* 3 = gibbmt */
	255, /* 4 = giim Giimport-id */
	60, /* 5 = genbank */
	60, /* 6 = embl */
	30, /* 7 = pir */
	20, /* 8 = swissprot */
	80,  /* 9 = patent */
	90, /* 10 = other TextSeqId */
	90, /* 11 = general Dbtag */
	120,  /* 12 = gi */
	60, /* 13 = ddbj */
	70, /* 14 = prf */
	50  /* 15 = pdb */
    };

#define PATENT_ORDER 110         /* order for any patent */

/*****************************************************************************
*
*   The above sets the ordering to be, lowest to highest
*

Nucleic Acids:
    GenBank/EMBL/DDBJ
    PDB
    Patents
    Anything else

Proteins:
    SWISSPROT
    PIR
    NCBI BackBone (but not in GenBank)
    PDB
    GenBank/EMBL/DDBJ translations
    PRF
    Patents
    Anything else

*
*****************************************************************************/



/*****************************************************************************
*
*   Traversal routine for SeqEntryToFasta
*
*****************************************************************************/
void	SeqEntryFasta (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)

{
    FastaPtr tfa;
    BioseqPtr bsp;
	BioseqSetPtr bssp;
    MyFsaPtr mfp;
    Boolean is_na;
	SeqIdPtr sip;
	TextSeqIdPtr tsip;
	ValNodePtr vnp;
	OrgRefPtr orp;
	BioSourcePtr bsrcp;
	MolInfoPtr mip;

    tfa = (FastaPtr) data;
    mfp = tfa->mfp;
	if (tfa->group_segs == 2)       /* put out only segments */
	{
		if (tfa->parts != -1)      /* in parts set */
		{
			if (indent <= tfa->parts)   /* out of parts set */
			{
				tfa->parts = -1;
				tfa->seg = -1;
			}
		}
	}

	if (IS_Bioseq(sep))
	{
		bsp = (BioseqPtr)(sep->data.ptrvalue);
		vnp = bsp->descr;
	}
	else
	{
		bssp = (BioseqSetPtr)(sep->data.ptrvalue);
		vnp = bssp->descr;
	}

	orp = NULL;
	mip = NULL;
	while (vnp != NULL)   /* check for organism info */
	{
		switch (vnp->choice)
		{
			case Seq_descr_source:
				bsrcp = (BioSourcePtr)(vnp->data.ptrvalue);
				orp = bsrcp->org;
				break;
			case Seq_descr_org:
				orp = (OrgRefPtr)(vnp->data.ptrvalue);
				break;
			case Seq_descr_molinfo:
				mip = (MolInfoPtr)(vnp->data.ptrvalue);
				break;
			default:
				break;
		}
		vnp = vnp->next;
	}

	if (orp != NULL)
	{
		if (orp->taxname != NULL)
			mfp->organism = orp->taxname;
		else if (orp->common != NULL)
			mfp->organism = orp->common;
	}

	if (mip != NULL)
            mfp->tech = mip->tech;
        else
            mfp->tech = 0 ;

   if (! IS_Bioseq(sep))	/* check for taking only parts of seg seqs */
	{
		if (tfa->group_segs == 2)    /* put out only segments */
		{
			if (bssp->_class == 2)   /* segset */
				tfa->seg = indent;
			else if (bssp->_class == 4)   /* parts */
			{
				if ((tfa->seg >= 0) && (tfa->seg < indent))
				{
					tfa->parts = indent;   /* in parts set */
				}
			}
		}
		return;
	}

   is_na = tfa->is_na;

	if ((! is_na) && (! ISA_aa(bsp->mol))) /* check for translations */
	{
		for (sip = bsp->id; sip != NULL; sip = sip->next)
		{
		 	switch (sip->choice)
			{
				case SEQID_GENBANK:
				case SEQID_EMBL:
				case SEQID_DDBJ:
					tsip = (TextSeqIdPtr)(sip->data.ptrvalue);
					if (tsip->accession != NULL)
						mfp->accession = tsip->accession;
					break;
				default:
					break;
			}
		}
	}

	if (tfa->last_indent != -1)   /* putting out segments together */
	{
		if (indent > tfa->last_indent)
			return;
		tfa->last_indent = -1;
	}
									 /* do raw bioseqs only */
	if (! tfa->group_segs)
	{
		if (BioseqRawToFastaX(bsp, mfp, is_na))
			tfa->got_one = TRUE;
	}
	else if (tfa->group_segs == 1)    /* do segmented sets */
	{
		if (BioseqToFastaX(bsp, mfp, is_na))
		{
			tfa->got_one = TRUE;
			if (bsp->repr == Seq_repr_seg)
				tfa->last_indent = indent;
		}
	}
	else if (tfa->group_segs == 2)    /* take only the parts */
	{
		if (tfa->parts >= 0)    /* in segmented parts set */
		{
			if (BioseqRawToFastaX(bsp, mfp, is_na))
				tfa->got_one = TRUE;
		}
	}
    return;
}

/*****************************************************************************
*
*   SeqEntryToFasta(sep, fp, is_na)
*
*****************************************************************************/
NLM_EXTERN Boolean SeqEntryToFasta (SeqEntryPtr sep, FILE *fp, Boolean is_na)
{
	if (IS_Bioseq(sep))
		return SeqEntrysToFasta(sep, fp, is_na, 3);
	else
		return SeqEntrysToFasta(sep, fp, is_na, 0);
}

/*****************************************************************************
*
*   FastaFileFunc(key, buf, data)
*   	standard "write to file" callback
*
*****************************************************************************/
NLM_EXTERN Boolean FastaFileFunc (BioseqPtr bsp, Int2 key, CharPtr buf,
                                  Uint4 buflen, Pointer data)
{
	FILE * fp;

	fp = (FILE *)data;

	switch (key)
	{
		case FASTA_ID:
			fprintf(fp, ">%s ", buf);
			break;
		case FASTA_DEFLINE:
			fprintf(fp, "%s\n", buf);
			break;
		case FASTA_SEQLINE:
			fprintf(fp, "%s\n", buf);
			break;
		case FASTA_EOS:   /* end of sequence */
			break;
		default:
			break;
	}
	return TRUE;
}

/*****************************************************************************
*
*   SeqEntrysToFasta(sep, fp, is_na, group_segs)
*
*   	group_segs = 0 ... take only raw Bioseqs
*       group_segs = 1 ... group segmented seqs into single entry.. no parts
*       group_segs = 2 ... show only parts of segmented seqs
*       group_segs = 3 ... like 1, but instantiate virtual Bioseqs
*   
*****************************************************************************/
NLM_EXTERN Boolean SeqEntrysToFasta (SeqEntryPtr sep, FILE *fp, Boolean is_na, Uint1 group_segs)

{
    FastaDat tfa;
    MyFsa mfa;
    Char buf[255];
    
    if ((sep == NULL) || (fp == NULL))
        return FALSE;
    
    mfa.buf = buf;
    mfa.buflen = 254;
    mfa.seqlen = 70;
    mfa.mydata = (Pointer)fp;
    mfa.myfunc = FastaFileFunc;
    mfa.bad_asn1 = FALSE;
    mfa.order = 0;
    mfa.accession = NULL;
    mfa.organism = NULL;
    mfa.do_virtual = FALSE;
    mfa.tech = 0;
    mfa.no_sequence = FALSE;
    mfa.formatdb	= FALSE;
    
    tfa.mfp = &mfa;
    tfa.is_na = is_na;

    if (is_na)
        mfa.code = Seq_code_iupacna;
    else
        mfa.code = Seq_code_ncbieaa;
    
    if (group_segs == 3)  /* do 2 things */
    {
        mfa.do_virtual = TRUE;
        group_segs = 1;
    }
    tfa.group_segs = group_segs;
    tfa.last_indent = -1;
    tfa.parts = -1;
    tfa.seg = -1;
    tfa.got_one = FALSE;
    SeqEntryExplore(sep, (Pointer)&tfa, SeqEntryFasta);
    return tfa.got_one;
}

/*****************************************************************************
*
*   SeqEntrysToFastaX(sep, mfa, is_na, group_segs)
*
*****************************************************************************/
NLM_EXTERN Boolean SeqEntrysToFastaX (SeqEntryPtr sep, MyFsaPtr mfp, Boolean is_na, Uint1 group_segs)

{
    FastaDat tfa;

	if ((sep == NULL) || (mfp == NULL))
		return FALSE;

    tfa.mfp = mfp;
    tfa.is_na = is_na;
	if (group_segs == 3)  /* do 2 things */
	{
		mfp->do_virtual = TRUE;
		group_segs = 1;
	}
	tfa.group_segs = group_segs;
	tfa.last_indent = -1;
	tfa.parts = -1;
	tfa.seg = -1;
	tfa.got_one = FALSE;
        SeqEntryExplore(sep, (Pointer)&tfa, SeqEntryFasta);
	return tfa.got_one;
}
/*****************************************************************************
*
*   SeqEntrysToDefline(sep, mfa, is_na, group_segs)
*
*****************************************************************************/
#define DEFLINE_MAX_LEN 2048

NLM_EXTERN Boolean SeqEntrysToDefline(SeqEntryPtr sep, 
                           FILE *fp, Boolean is_na, Uint1 group_segs)

{
  FastaDat tfa;
  MyFsa mfa;
  
  if ((sep == NULL) || (fp == NULL))
    return FALSE;
  
  mfa.buf = (CharPtr) MemNew(DEFLINE_MAX_LEN);
  mfa.buflen = DEFLINE_MAX_LEN-1;
  mfa.seqlen = DEFLINE_MAX_LEN;
  mfa.mydata = (Pointer)fp;
  mfa.myfunc = FastaFileFunc;
  mfa.no_sequence = TRUE;
  mfa.bad_asn1 = FALSE;
  mfa.order = 0;
  mfa.accession = NULL;
  mfa.organism = NULL;
  mfa.do_virtual = FALSE;
  mfa.formatdb = FALSE;
  mfa.tech = 0;

  tfa.mfp = &mfa;
  tfa.is_na = is_na;
  if (group_segs == 3)  /* do 2 things */
    {
      mfa.do_virtual = TRUE;
      group_segs = 1;
    }
  tfa.group_segs = group_segs;
  tfa.last_indent = -1;
  tfa.parts = -1;
  tfa.seg = -1;
  tfa.got_one = FALSE;
  SeqEntryExplore(sep, (Pointer)&tfa, SeqEntryFasta);
  MemFree(mfa.buf);
  return tfa.got_one;
}

/*****************************************************************************
*
*   Boolean BioseqRawToFasta(bsp, fp, is_na)
*
*****************************************************************************/
NLM_EXTERN Boolean BioseqRawToFasta (BioseqPtr bsp, FILE *fp, Boolean is_na)

{
	return BioseqRawToFastaExtra(bsp, fp, 80);
}

NLM_EXTERN Boolean BioseqRawToFastaExtra (BioseqPtr bsp, FILE *fp, Int2 line_length)

{
    MyFsa mfa;
    Char buf[255];

    if ((bsp == NULL) || (fp == NULL))
        return FALSE;

	mfa.buf = buf;
	mfa.buflen = 254;
	mfa.seqlen = line_length;
	mfa.mydata = (Pointer)fp;
	mfa.myfunc = FastaFileFunc;
	mfa.bad_asn1 = FALSE;
	mfa.order = 0;
	mfa.accession = NULL;
	mfa.organism = NULL;
	mfa.do_virtual = FALSE;
	mfa.tech = 0;
        mfa.no_sequence = FALSE;
    	mfa.formatdb	= FALSE;
    

 	return BioseqRawToFastaX(bsp, &mfa, ISA_na(bsp->mol));
}

/*****************************************************************************
*
*   Boolean BioseqRawToFastaX(bsp, mfp, is_na)
*
*****************************************************************************/
NLM_EXTERN Boolean BioseqRawToFastaX (BioseqPtr bsp, MyFsaPtr mfp, Boolean is_na)

{
	Uint1 repr;

    if ((bsp == NULL) || (mfp == NULL))
        return FALSE;

    repr = Bioseq_repr(bsp);
    if (! ((repr == Seq_repr_raw) || (repr == Seq_repr_const)))
        return FALSE;

	return BioseqToFastaX(bsp, mfp, is_na);
}

/*****************************************************************************
*
*   Boolean BioseqToFasta(bsp, fp, is_na)
*
*****************************************************************************/
NLM_EXTERN Boolean BioseqToFasta (BioseqPtr bsp, FILE *fp, Boolean is_na)

{
	MyFsa mfa;
	Char buf[255];

    if ((bsp == NULL) || (fp == NULL))
        return FALSE;

	mfa.buf = buf;
	mfa.buflen = 254;
	mfa.seqlen = 80;
	mfa.mydata = (Pointer)fp;
	mfa.myfunc = FastaFileFunc;
	mfa.bad_asn1 = FALSE;
	mfa.order = 0;
	mfa.accession = NULL;
	mfa.organism = NULL;
	mfa.do_virtual = FALSE;
	mfa.tech = 0;
        mfa.no_sequence = FALSE;
    	mfa.formatdb	= FALSE;

	return BioseqToFastaX(bsp, &mfa, is_na);
}

/*****************************************************************************
*
*   Boolean BioseqToFastaX(bsp, mfp, is_na)
*
*****************************************************************************/
NLM_EXTERN Boolean BioseqToFastaX (BioseqPtr bsp, MyFsaPtr mfp, Boolean is_na)

{
    SeqPortPtr spp;
    Uint1 repr, code;
    Char buf[41];
    SeqIdPtr sip;
    Uint1 order = 255;
    Boolean is_patent = FALSE, is_genbank = FALSE;
    Uint1Ptr order_array;
    int i;
    CharPtr organism = NULL;
    

    if ((bsp == NULL) || (mfp == NULL))
        return FALSE;

    repr = Bioseq_repr(bsp);

    if (ISA_na(bsp->mol))
    {
        if (! is_na)
            return FALSE;
		order_array = na_order;
    }
    else if (ISA_aa(bsp->mol))
    {
		if (is_na)
			return FALSE;
		order_array = aa_order;
		if (mfp->accession != NULL)		   /* translated genbank */
		{
			order = order_array[SEQID_GENBANK];
			is_genbank = TRUE;
			organism = mfp->organism;
		}
    }
    else
    {
	   buf[0] = '\0';
		SeqIdWrite(bsp->id, buf, PRINTID_FASTA_LONG, 40);
	ErrPostEx(SEV_ERROR,0,0,"ToFasta: [%s] Unrecognized bsp->mol = %d",
		  buf, (int)(bsp->mol));
	mfp->bad_asn1 = TRUE;
	return FALSE;
    }

	mfp->bsp = bsp;

	for (sip = bsp->id; sip != NULL; sip = sip->next)
	{
		i=(int)(sip->choice);
		if (! is_genbank)    /* don't change order for translated genbank */
		{
			if (order_array[i] < order)
				order = order_array[i];
		}
		if (i == (int)SEQID_PATENT)
			is_patent = TRUE;
		else if (i == (int)SEQID_PRF)
			organism = mfp->organism;
	}

	if (is_patent)
		order = PATENT_ORDER;

	mfp->order = order;

	switch (mfp->tech)
	{
		case MI_TECH_est:
		case MI_TECH_sts:
		case MI_TECH_survey:
		case MI_TECH_htgs_1:
		case MI_TECH_htgs_2:
		case MI_TECH_htgs_3:
			organism = mfp->organism;
			break;
		default:
			break;
	}

	if (! FastaId(bsp, mfp->buf, mfp->buflen))
		return FALSE;

	(*(mfp->myfunc))(bsp, FASTA_ID, mfp->buf, StringLen(mfp->buf), mfp->mydata);

   	if (! CreateDefLine(NULL, bsp, mfp->buf, mfp->buflen, mfp->tech, mfp->accession, organism))
   		return FALSE;

	(*(mfp->myfunc))(bsp, FASTA_DEFLINE, mfp->buf, StringLen(mfp->buf), mfp->mydata);

        if (mfp->formatdb && is_na) {
            (*(mfp->myfunc))(bsp, FASTA_FORMATDB_AMB, mfp->buf, StringLen(mfp->buf), mfp->mydata);
        }
        else if(!mfp->no_sequence) {
	    if (!mfp->formatdb) {
		if (is_na)
		    code = Seq_code_iupacna;
		else
		    code = Seq_code_ncbieaa;
	    } else {
		code = mfp->code;
	    }

            spp = FastaSeqPort(bsp, is_na, mfp->do_virtual, code);
            if (spp == NULL) return FALSE;
            
            while (FastaSeqLine(spp, mfp->buf, mfp->seqlen, is_na))
                (*(mfp->myfunc))(bsp, FASTA_SEQLINE, mfp->buf, StringLen(mfp->buf),
                                 mfp->mydata);
            
            SeqPortFree(spp);
            
            (*(mfp->myfunc))(bsp, FASTA_EOS, mfp->buf, StringLen(mfp->buf),
                             mfp->mydata);
        }
        return TRUE;
}


/*****************************************************************************
*
*   Here are functions that convert FASTA format from file or from memory
*
*****************************************************************************/

/********* DEFINES *********/

#define FASTA_MEM_IO  1   /* type of reading from buffer in memory */
#define FASTA_FILE_IO 2   /* type of reading from file */
#define FTSE_BUFF_CHUNK 4096
#define BIOSEQ 1

/********* INTERNAL FUNCTIONS *********/

NLM_EXTERN SeqEntryPtr FastaToSeqEntryInternal
(
 VoidPtr input,          /* input pointer (file or memory) */
 Int4 type,              /* type of inquiry FASTA_MEM_IO or FASTA_FILE_IO */
 CharPtr PNTR last_char, /* returned pointer to next FASTA sequence */
 Boolean is_na,          /* type of sequence */
 CharPtr PNTR errormsg,  /* error messge for debugging */
 Boolean parseSeqId,     /* Parse SeqID from def line */
 CharPtr special_symbol  /* Returns special symbol if no SeqEntry */
 );
static Boolean FastaReadSequenceInternal
(
 VoidPtr input,          /* input pointer (file or memory) */
 Int4 type,              /* type of inquiry FASTA_MEM_IO or FASTA_FILE_IO */
 CharPtr PNTR last_char, /* returned pointer to next FASTA sequence */
 Boolean is_na,          /* type of sequence */
 Int4Ptr seq_length,     /* Returned length of sequence in residues */
 ByteStorePtr PNTR,      /* Returned pointer to sequence ByteStore */
 CharPtr PNTR errormsg,  /* error messge for debugging */
 CharPtr special_symbol  /* Returns special symbol if no SeqEntry */
 );
static Int4 FastaReadSequenceChunk
( 
 VoidPtr input,          /* input pointer (file or memory) */
 Int4    type,           /* type of inquiry FASTA_MEM_IO or FASTA_FILE_IO */
 CharPtr PNTR next_char, /* returned pointer to next FASTA sequence */
 Uint1Ptr sequence,      /* buffer to read sequence to */
 Int4     length,        /* size of buffer */
 CharPtr special_symbol  /* Returns special symbol if no SeqEntry */
 );

/********* FINCTIONS *********/

/*****************************************************************************
*
*   SeqEntryPtr FastaToSeqBuffEx() - function to return SeqEntryPtr from
*                                    buffer with error handling
*
*****************************************************************************/
NLM_EXTERN SeqEntryPtr FastaToSeqBuffEx
  (
    CharPtr buffer,         /* buffer in memory with FASTA sequence */
    CharPtr PNTR last_char, /* here returned pointer to next FASTA if any */
    Boolean is_na,          /* type of sequence */
    CharPtr PNTR errormsg,  /* error message for debugging */
    Boolean parseSeqId      /* Parse SeqID from def line */
  )
{
  return FastaToSeqEntryInternal((void *)buffer, FASTA_MEM_IO , 
                                 last_char, is_na, errormsg, parseSeqId, NULL);
}

/*****************************************************************************
*
*   SeqEntryPtr FastaToSeqEntryEx() - function to return SeqEntryPtr from
*                                     file with error handling
*
*****************************************************************************/
NLM_EXTERN SeqEntryPtr FastaToSeqEntryEx 
  (
    FILE *fp,               /* file to get sequence from */ 
    Boolean is_na,          /* type of sequence */
    CharPtr PNTR errormsg,  /* error message for debugginq */
    Boolean parseSeqId      /* Parse SeqID from def line */
  )
{
  return FastaToSeqEntryInternal((void *)fp, FASTA_FILE_IO,
                                 NULL,is_na, errormsg, parseSeqId, NULL);	 
}

/*****************************************************************************
*
*   SeqEntryPtr FastaToSeqEntry() - function to return SeqEntryPtr from
*                                   file without error handling
*
*****************************************************************************/
NLM_EXTERN SeqEntryPtr FastaToSeqEntry (FILE *fp, Boolean is_na)
{
  return FastaToSeqEntryEx (fp, is_na, NULL, TRUE);
}

/*****************************************************************************
*
*   SeqEntryPtr FastaToSeqBuff() - function to return SeqEntryPtr from
*                                   buffer without error handling
*
*****************************************************************************/
NLM_EXTERN SeqEntryPtr FastaToSeqBuff (CharPtr buffer, CharPtr PNTR last_char, 
                                   Boolean is_na)
{
  return FastaToSeqBuffEx (buffer, last_char, is_na, NULL, TRUE);
}

/*****************************************************************************
*
*   Boolean FastaReadSequenceChunk() - read sequence chunkfrom
*                                      file or buffer for use in
*                                      FastaReadsequenceInternal()
*****************************************************************************/
static Int4 FastaReadSequenceChunk
( 
 VoidPtr input,          /* input pointer (file or memory) */
 Int4    type,           /* type of inquiry FASTA_MEM_IO or FASTA_FILE_IO */
 CharPtr PNTR next_char, /* returned pointer to next FASTA sequence */
 Uint1Ptr sequence,      /* buffer to read sequence to */
 Int4     length,        /* size of buffer */
 CharPtr special_symbol  /* Returns special symbol if no SeqEntry */
 )

{
    const Char PNTR firstchar;
    FILE *fd;
    register Int4 i;
    Int2 ch;
    
    /* Type of input depends upon calling function */
    
    if(type == FASTA_FILE_IO) {
        fd = (FILE *) input;
    } else {   /* type == FASTA_MEM_IO */
        if((firstchar = (const Char PNTR) input) == NULL)
            return 0;
    }
    
    if(type == FASTA_FILE_IO) {
        for(i=0; i < length && (ch = getc(fd)) != EOF; i++) {      
            if((sequence[i] = (Char) ch) == '>' || (Char) ch == '&' || (Char) ch == '{' ||
                (Char) ch == '}' || (Char) ch == '[' || (Char) ch == ']') {
                if((i == 0) || 
                   (i > 0 && (sequence[i-1] == '\n' || 
                              sequence[i-1] == '\r'))) {
                    ungetc(ch, fd);
                    if (special_symbol != NULL) {
                        *special_symbol = (Char) ch;
                    }
                    break;
                }
            }
        }
    } else { /* type = FASTA_MEM_IO */
        for(i =0; i < length && (ch = *firstchar) != NULLB; firstchar++, i++) {
            if((sequence[i] = (Char) ch) == '>' || (Char) ch == '&' || (Char) ch == '{' ||
                (Char) ch == '}' || (Char) ch == '[' || (Char) ch == ']') {
                if((i == 0) || 
                   (i > 0 && (sequence[i-1] == '\n' || 
                              sequence[i-1] == '\r'))) {
                    
                    if (special_symbol != NULL) {
                        *special_symbol = (Char) ch;
                    }
                    break;
                }
            }
        }   
        
        if(ch == NULLB) /* the end of buffer */
            *next_char = NULL; 
        else 
            *next_char = (CharPtr) firstchar; 
    }
    return i;
}
/*****************************************************************************
*
*   Boolean FastaReadSequence() - read sequence from file
*
*****************************************************************************/

Boolean FastaReadSequence
(
 FILE *fd,            /* input pointer (file or memory) */
 Boolean is_na,            /* type of sequence */
 Int4Ptr seq_length,       /* Returned length of sequence in residues */
 ByteStorePtr PNTR bs_out, /* Returned pointer to sequence ByteStore */
 CharPtr PNTR errormsg     /* error message for debugging */
 )
{
    return  FastaReadSequenceInternal((VoidPtr) fd, FASTA_FILE_IO, NULL,
                                      is_na, seq_length, bs_out, errormsg, NULL);
}
/*****************************************************************************
*
*   Boolean FastaReadSequenceMem() - read sequence from buffer
*
*****************************************************************************/

Boolean FastaReadSequenceMem
(
 CharPtr buffer,           /* input buffer with sequence */
 CharPtr PNTR next_char,   /* returned pointer to next FASTA sequence */
 Boolean is_na,            /* type of sequence */
 Int4Ptr seq_length,       /* Returned length of sequence in residues */
 ByteStorePtr PNTR bs_out, /* Returned pointer to sequence ByteStore */
 CharPtr PNTR errormsg     /* error message for debugging */
)
{
    return  FastaReadSequenceInternal((VoidPtr) buffer, FASTA_MEM_IO, 
                                      next_char, is_na, seq_length, bs_out, 
                                      errormsg, NULL);
}
/*****************************************************************************
*
*   Boolean FastaReadSequenceInternal() - read sequence from
*                                         file or buffer for internal use
*
*****************************************************************************/

static Boolean FastaReadSequenceInternal
(
 VoidPtr input,            /* input pointer (file or memory) */
 Int4 type,                /* type of inquiry FASTA_MEM_IO or FASTA_FILE_IO */
 CharPtr PNTR next_char,   /* returned pointer to next FASTA sequence */
 Boolean is_na,            /* type of sequence */
 Int4Ptr seq_length,       /* Returned length of sequence in residues */
 ByteStorePtr PNTR bs_out, /* Returned pointer to sequence ByteStore */
 CharPtr PNTR errormsg,    /* error message for debugging */
 CharPtr special_symbol  /* Returns special symbol if no SeqEntry */
)
{
    SeqMapTablePtr smtp;
    Uint1Ptr       in_buff, out_buff;
    CharPtr        ptr, chptr;
    Int2           ch;
    Uint1          byte_from, uch;
    register Int4  i;
    CharPtr        badchar = NULL;
    Int4           in_index, out_index, total_read, badchars = 0;
    Int4           total_length = 0;
    Char           tmp[32];
    Boolean        Second, skip_to_eol, last_was_star;
    
    if (input == NULL)     /* empty input */
        return FALSE;
     
    /* Initializing conversion tables */
    
    if(is_na) {
        if((smtp = SeqMapTableFind(Seq_code_ncbi4na, 
                                   Seq_code_iupacna)) == NULL) {
            return FALSE;
        }
    } else {
        if((smtp = SeqMapTableFind(Seq_code_ncbistdaa, 
                                   Seq_code_ncbieaa)) == NULL) {
            return FALSE;
        }
    }

    /* Allocationg error message buffers if required */
    
    if (errormsg != NULL) {
        *errormsg = NULL;        
        if((badchar = (CharPtr) MemNew(256)) == NULL)
            return FALSE;
    }

    if((in_buff = (Uint1Ptr) MemNew(FTSE_BUFF_CHUNK)) == NULL)
        return FALSE;

    if((out_buff = (Uint1Ptr) MemNew(FTSE_BUFF_CHUNK)) == NULL)
        return FALSE;

    if((*bs_out = BSNew(FTSE_BUFF_CHUNK)) == NULL)
        return FALSE;

    Second = FALSE;
    skip_to_eol = FALSE;
    last_was_star = FALSE;
    in_index = out_index = total_read = 0;

    while(TRUE) {
        if (in_index == total_read) {
            if((total_read = FastaReadSequenceChunk(input, type, 
                                                    next_char, in_buff, 
                                                    FTSE_BUFF_CHUNK, special_symbol)) == 0)
                break; /* Here is exit from the loop */
	 
            if(type == FASTA_MEM_IO)
                input = (VoidPtr) *next_char;
	  
            in_index = 0;
	}
      	byte_from = in_buff[in_index];
      	in_index++;
        if ((! is_na) && (! last_was_star) && byte_from == '*') {
        	last_was_star = TRUE;
        } else if(byte_from != ';' && !skip_to_eol) {
            byte_from = TO_UPPER (byte_from);
                
            if (is_na && byte_from == 'U') byte_from = 'T';
                
            if((uch = SeqMapTableConvert(smtp, byte_from)) != 
               INVALID_RESIDUE) {
                if (last_was_star) {
                	total_length++;
                	out_buff[out_index] = SeqMapTableConvert(smtp, '*');
                	out_index++;
                	if(out_index == FTSE_BUFF_CHUNK) {
                		if(BSWrite(*bs_out, out_buff, out_index) != out_index) {
                			MemFree (badchar);
                			MemFree(in_buff);
                			MemFree(out_buff);                
                			return FALSE;                
                		}
                		out_index = 0;
                	}
                	last_was_star = FALSE;
                }
                total_length++;
                if(is_na) {
                    if(!Second) {
                        uch <<= 4;   
                        out_buff[out_index] = uch;
                    } else {
                        out_buff[out_index] += uch;
                        out_index++;
                    } 
                    Second = !Second;
                } else {
                    out_buff[out_index] = uch;
                    out_index++;
                }
            } else if (errormsg != NULL){ 
                if(IS_ALPHA(byte_from) || byte_from == '-') {
                    (badchar [(int) (byte_from)])++;
                    badchars++;
                }
            }
        } else {    /* ch == ';' */

            /* We have to ignore rest of the line */
            skip_to_eol = TRUE;
            
            while(in_index < total_read  && 
                  (byte_from = in_buff[in_index]) != '\n' && 
                  byte_from != '\r')
                in_index++;
            if(in_index < total_read)
                skip_to_eol = FALSE;                    
        }
        
        if(out_index == FTSE_BUFF_CHUNK) {
            if(BSWrite(*bs_out, out_buff, out_index) != out_index) {
            	MemFree (badchar);
                MemFree(in_buff);
                MemFree(out_buff);                
                return FALSE;                
            }
            out_index = 0;
        }
    }  /* while (TRUE) */ 

    /* We have to write remaining stuff in out_buff */

    if(is_na && Second) out_index++; /* Partial byte for DNA */

    if(BSWrite(*bs_out, out_buff, out_index) != out_index) {
    	MemFree (badchar);
        MemFree(in_buff);
        MemFree(out_buff);                
        return FALSE;
    }
    
    *seq_length = total_length;
    
    /* If required bad characters statistics */
    
    if (errormsg != NULL && badchars > 0) {
        if((ptr = (CharPtr) MemNew (sizeof(Char)*512)) == NULL)
            return FALSE;
        
        chptr = "";
        sprintf (ptr, "%ld illegal %s %s removed:\n", badchars, 
                 badchars == 1 ? "character" : "characters",
                 badchars == 1 ? "was" : "were"
                 );
        for (ch = 'A', i =0; ch <= 'Z'; ch++, i++) {
            if ((badchar[ch]) > 0) {
                sprintf (tmp, "%s%d %c%s", 
                         chptr, badchar[ch], ch,
                         badchar[ch] == 1 ? "" : "s");
                StringCat (ptr, tmp);
                chptr = ", ";
            }
        }
        ch = '-';
        if ((badchar[ch]) > 0) {
            sprintf (tmp, "%s%d %c%s", 
                     chptr, badchar[ch], ch,
                     badchar[ch] == 1 ? "" : "s");
            StringCat (ptr, tmp);
            chptr = ", ";
        }
        *errormsg = StringSave (ptr);
        MemFree (ptr);
    }
    
    MemFree (badchar);
    MemFree(in_buff);
    MemFree(out_buff);
    return TRUE;
}

NLM_EXTERN SeqEntryPtr FastaToSeqEntryInternal
(
 VoidPtr input,          /* input pointer (file or memory) */
 Int4 type,              /* type of inquiry FASTA_MEM_IO or FASTA_FILE_IO */
 CharPtr PNTR next_char, /* returned pointer to next FASTA sequence */
 Boolean is_na,          /* type of sequence */
 CharPtr PNTR errormsg,  /* error messge for debugging */
 Boolean parseSeqId,     /* Parse SeqID from def line */
 CharPtr special_symbol  /* Returns special symbol if no SeqEntry */
 )

{
    SeqEntryPtr    sep = NULL;
    BioseqPtr      bsp = NULL;
    ValNodePtr     vnp = NULL;
    Int2           ch;
    CharPtr        chptr = NULL, ptr = NULL;
    register Int4  i;
    CharPtr        defline, buffer= NULL;   /* Working buffers */    
    Int4           len = 0, BuffSize = FTSE_BUFF_CHUNK; 
    FILE           *fd;    
    const Char     PNTR firstchar;
  

    if (special_symbol != NULL) {
        *special_symbol = '\0';
    }
    if (input == NULL)     /* empty input */
        return NULL;

    /* Type of input depends upon calling function */
    
    if(type == FASTA_FILE_IO)
        fd = (FILE *) input;
    else    /* type == FASTA_MEM_IO */
        firstchar = (const Char PNTR) input;
    
    /* Rolling spaces to check first non-space character */
    
    if(type == FASTA_FILE_IO) {
        do {
            ch = getc(fd);
        } while (IS_WHITESP(ch));
    } else {   /* if(type == FASTA_MEM_IO*/
        while (IS_WHITESP(ch = *firstchar)) /* Rolling spaces */
            firstchar++;  
    }
    
    if(ch == EOF || ch == NULLB || ch == '&' || ch == '{' ||
       ch == '}' || ch == '[' || ch == ']') {
        /* This is empty FILE or buffer or special symbol detected */
        if (special_symbol != NULL) {
            *special_symbol = ch;
        }
        return NULL;
    }

    /* First character is valid: initializing main structures */
        
    /* Initializing Seq-entry structure */
    
    if((sep = SeqEntryNew()) == NULL) {
        MemFree(buffer);
        return NULL;
    }
    sep->choice = BIOSEQ;  /* == 1 */

    if((bsp = BioseqNew()) == NULL) {
        MemFree(buffer);
        return NULL;
    }
    
    sep->data.ptrvalue = bsp;
    SeqMgrSeqEntry (SM_BIOSEQ, (Pointer)bsp, sep); 
    if (is_na) {
        bsp->mol = Seq_mol_na;
        bsp->seq_data_type = Seq_code_ncbi4na;
    } else {
        bsp->mol = Seq_mol_aa;
        bsp->seq_data_type = Seq_code_ncbistdaa;
    }
    bsp->repr = Seq_repr_raw;

    SeqMgrAddToBioseqIndex (bsp);
    
    /*  ------------- */
    
    /* Now reading defline into memory */
    
  /* DEFLINE PROCCESSING*/
    
    if(ch == '>') {     /* Defline is present - processing */
        if((buffer = (CharPtr) MemNew(BuffSize)) == NULL)
            return NULL;        
        
        if(type == FASTA_FILE_IO) {    /* File */
            buffer[0] = (Char) ch;
            
            for(i=1; (ch = getc(fd)) != EOF; i++) {      
                
                if (i >= BuffSize) {        /* Reallocating defline buffer */
                    BuffSize = i + FTSE_BUFF_CHUNK;
                    if((buffer = (CharPtr)Realloc(buffer, 
                                                  BuffSize)) == NULL) {
                        ErrLogPrintf("Error re-allocating memory in "
                                     "FastaToSeqEntry");
                        MemFree(buffer);
                        return NULL;
                    }
                }
                if((buffer[i] = (Char) ch) == '\n' || ch == '\r') {
                    break;
                }
            }
            buffer[i] = NULLB;
            
        } else {  /* type = FASTA_MEM_IO */
            for(i =0; (ch = *firstchar) != NULLB; firstchar++, i++) {
                if (i >= BuffSize) {
                    BuffSize = i + FTSE_BUFF_CHUNK;
                    buffer = (CharPtr) Realloc(buffer, BuffSize);
                }
                if((buffer[i] = (Char) ch) == '\n' || ch == '\r') {
                    break;
                }
            }
            buffer[i] = NULLB;
        
            if(ch == NULLB) {/* the end of buffer */
                *next_char = NULL;
                input =  (VoidPtr) "\0";
            } else { 
                *next_char = (CharPtr) firstchar; 
                input = (VoidPtr) firstchar;
            }
        }
        
        defline = buffer+1;   /* Character after '>' */
        
        if(defline[0] != '?') { 
            /* Creating standard Seq-id */        
            ptr = defline;
            
            while (IS_WHITESP(*ptr))
                ptr++;
            
            if (parseSeqId) {
                if (*ptr == '"') {
                    ptr++;
                    chptr = StringChr (ptr, '"');
                } else {
                    chptr = StringChr (ptr, ' ');
                }
            }
            if (!parseSeqId) {
                chptr = ptr;
            } else if (chptr != NULL) {
                *chptr = NULLB;
                chptr++;
                bsp->id = MakeSeqID (ptr);
            } else if (*ptr != NULLB) {
                bsp->id = MakeSeqID (ptr);
            }
            
            if (bsp->id == NULL) 
                bsp->id = MakeNewProteinSeqId (NULL, NULL);
            if (chptr != NULL) {
                if((vnp = ValNodeNew(NULL)) != NULL) {
                    vnp->choice = Seq_descr_title;
                    while (IS_WHITESP(*chptr))
                        chptr++;
                    vnp->data.ptrvalue = StringSave (chptr);
                }
                bsp->descr = vnp;
            }
        } else {
            /* Unknown Seq-id */ 
            bsp->id = MakeSeqID ("lcl|gap");
            bsp->repr = Seq_repr_virtual;
            
            if(sscanf(defline + 2, "%ld", &len) == 1 && len > 0) {
                bsp->length =  len;
            } else {
                bsp->length = -1;
            }
        }

        MemFree(buffer);
        
    } else {  /* if ch == '>' EMPTY DEFLINE */ 
        /* Defline is upsent - creating default defline */
        bsp->id = MakeNewProteinSeqId (NULL, NULL);

        if(type == FASTA_FILE_IO)
            ungetc(ch, fd);
    }
        
    /* OK, now processing sequence */
        
    if(!FastaReadSequenceInternal(input, type, next_char, is_na,
                                  &bsp->length, &bsp->seq_data,
                                  errormsg, special_symbol)) {
        ErrPostEx(SEV_FATAL, 0, 0, "Failure to read sequence. "
                  "FastaToSeqEntry() failed.\n");
        return NULL;
    }
    
    BioseqPack(bsp);     /* Trying to pack Bioseq more */
    return sep;
}

/*****************************************************************************
*
*   FastaId(bsp, buf, buflen)
*      Makes the string for the id part of fasta format.
*      buf should be at least 40 bytes
*
*****************************************************************************/
NLM_EXTERN Boolean FastaId(BioseqPtr bsp, CharPtr buf, Int2 buflen)
{
	if ((bsp == NULL) || (buf == NULL)) return FALSE;

    SeqIdWrite(bsp->id, buf, PRINTID_FASTA_LONG, buflen);
	return TRUE;
}

/*****************************************************************************
*
*   FastaDefLine(bsp, buf, buflen, accession, organism)
*   	Finds or makes a FASTA format defline (just locates the string)
*       buf should be very long if possible
*       function truncates if buf not long enough
*       a few deflines are longer than 255
*
*****************************************************************************/
NLM_EXTERN Boolean FastaDefLine (BioseqPtr bsp, CharPtr buf, Int2 buflen,
                                          CharPtr accession, CharPtr organism, Uint1 tech)
{
	BioseqContextPtr bcp;
	ValNodePtr vnp;
	CharPtr tmp;
	PdbBlockPtr pbp;
	PatentSeqIdPtr psip;
	Int2 diff, phase;
	Int4 num_segs, num_gaps;
	Char tbuf[80];
	static CharPtr htgs[2] = {
		"unordered", "ordered" };

	if ((bsp == NULL) || (buf == NULL)) return FALSE;

	buflen--;
	buf[buflen] = '\0';

	if (accession != NULL)
	{
		diff = LabelCopyExtra(buf, accession, buflen, "(", ") ");
		buflen -= diff;
		buf += diff;
	}

	bcp = BioseqContextNew(bsp);
	diff = 0;

	if ((tmp = BioseqContextGetTitle(bcp)) != NULL) {
		diff = LabelCopy(buf, tmp, buflen);
		                        /* remove trailing blanks and periods */
		tmp = buf + diff - 1;   /* point at last character */
		while (((*tmp <= ' ') || (*tmp == '.')) && (diff))
		{
			*tmp = '\0';
			tmp--; diff--;
		}
	}
	else
		if ((vnp = BioseqContextGetSeqDescr(bcp, Seq_descr_pdb, NULL, NULL)) != NULL)
	{
		pbp = (PdbBlockPtr)(vnp->data.ptrvalue);
		diff = LabelCopy(buf, (CharPtr)(pbp->compound->data.ptrvalue), buflen);
	}
	else
	{
		for (vnp = bsp->id; vnp != NULL; vnp = vnp->next)
		{
			if (vnp->choice == SEQID_PATENT)
			{
				psip = (PatentSeqIdPtr)(vnp->data.ptrvalue);
				sprintf(tbuf, "Sequence %d from Patent %s %s",
					(int)psip->seqid, psip->cit->country, psip->cit->number);
				diff = LabelCopy(buf, tbuf, buflen);
				break;
			}
		}
		if (vnp == NULL)
			diff = LabelCopy(buf, "No definition line found", buflen);
	}
	buflen -= diff;
	buf += diff;
	
	BioseqContextFree(bcp);

	if ((tech >= MI_TECH_htgs_1) && (tech <= MI_TECH_htgs_3))
	{
		phase = (Int2)(tech - MI_TECH_htgs_1 + 1);
		sprintf(tbuf, "; HTGS phase %d", (int)phase);
		diff = LabelCopy(buf, tbuf, buflen);
		buflen -= diff;
		buf += diff;
		
		if (phase == 3)
		{
			if (tmp && StringStr(tmp, "complete sequence") == NULL) {
				diff = LabelCopy(buf, ", complete sequence", buflen);
				buflen -= diff;
				buf += diff;
			}
		}
		else if (bsp->repr == Seq_repr_delta)
		{
			if (CountGapsInDeltaSeq(bsp, &num_segs, &num_gaps, NULL, NULL, NULL, 0))
			{
				sprintf(tbuf, ", %ld %s pieces", (long)(num_segs - num_gaps),
					htgs[phase - 1]);
				diff = LabelCopy(buf, tbuf, buflen);
				buflen -= diff;
				buf += diff;
			}
		}

	}

	if (organism != NULL)
	{
		LabelCopyExtra(buf, organism, buflen, " [", "]");
	}
	return TRUE;
}

static Boolean is_pdb(BioseqPtr bsp)
{
	SeqIdPtr id;
	
	if (bsp ==NULL)
		return FALSE;
	for (id = bsp->id; id; id=id->next)
	{
		if (id->choice == SEQID_PDB)
			return TRUE;
	}
	return FALSE;
}

static ValNodePtr tie_next(ValNodePtr head, ValNodePtr next)
{
   ValNodePtr v;

   if (head == NULL) {
      return next;
   }
   for (v = head; v->next != NULL; v = v->next)
   		continue;    
   v->next = next;
   return head;
}

static Boolean get_descr_on_top (GatherContextPtr gcp)
{
    ValNodePtr	tmp;
    DescrInfoPtr	PNTR dspp;
    DescrInfoPtr	dsp;
    ItemInfoPtr 	iip;
    
    dspp = (DescrInfoPtr PNTR) gcp->userdata;
    dsp = *dspp;
    switch (gcp->thistype) {
    case OBJ_SEQDESC:
        tmp = (ValNodePtr) (gcp->thisitem);
        if (tmp->choice == dsp->choice) {
            if (tmp->data.ptrvalue != NULL) {
                dsp->vnp = tmp;
                iip = (ItemInfoPtr) MemNew(sizeof(ItemInfo));
                dsp->iip = iip;
                iip->entityID = gcp->entityID;
                iip->itemID = gcp->itemID;
                iip->itemtype = gcp->thistype;
            }
        } 
        break;
    default:
        break;
    }
    return TRUE;
}

static Boolean get_descr (GatherContextPtr gcp)
{
	ValNodePtr	tmp;
	DescrInfoPtr	PNTR dspp;
	DescrInfoPtr	dsp;
	ItemInfoPtr 	iip;
	BioseqPtr 		bsp;
	
	dspp = (DescrInfoPtr PNTR) gcp->userdata;
	dsp = *dspp;
	switch (gcp->thistype)
	{
		case OBJ_SEQDESC:
			tmp = (ValNodePtr) (gcp->thisitem);
			if (tmp->choice == dsp->choice) {
				bsp = (BioseqPtr) (gcp->parentitem);
				if (dsp->bsp != bsp) {
					break;
				}
				if (tmp->data.ptrvalue != NULL) {
					dsp->vnp = tmp;
					iip = (ItemInfoPtr) MemNew(sizeof(ItemInfo));
					dsp->iip = iip;
					iip->entityID = gcp->entityID;
					iip->itemID = gcp->itemID;
					iip->itemtype = gcp->thistype;
				}
			} 
			break;
		default:
			break;
	}
	return TRUE;
}

static Boolean GetFeatProt (GatherContextPtr gcp)
{
	ValNodePtr	PNTR vnpp;
	ValNodePtr tmp;
	SeqFeatPtr	sfp;
	
	vnpp = (ValNodePtr PNTR) gcp->userdata;
	switch (gcp->thistype)
	{
		case OBJ_SEQFEAT:
			sfp = (SeqFeatPtr) (gcp->thisitem);
			if (sfp->data.choice == SEQFEAT_PROT) {
				tmp = ValNodeNew(NULL);
				tmp->data.ptrvalue = sfp;
				*vnpp = tie_next(*vnpp, tmp);
			} 
			break;
		default:
			break;
	}
	return TRUE;
}

static Boolean GetFeatCDS (GatherContextPtr gcp)
{
	SeqFeatPtr	PNTR sfpp;
	SeqFeatPtr	sfp;
	
	sfpp = (SeqFeatPtr PNTR) gcp->userdata;
	switch (gcp->thistype)
	{
		case OBJ_SEQFEAT:
			sfp = (SeqFeatPtr) (gcp->thisitem);
			if (sfp->data.choice == SEQFEAT_CDREGION) {
				*sfpp = sfp;
				return FALSE;
			} 
			break;
		default:
			break;
	}
	*sfpp = NULL;
	return TRUE;
}

static Boolean GetFeatGenes (GatherContextPtr gcp)
{
	ValNodePtr	PNTR vnpp;
	ValNodePtr tmp;
	SeqFeatPtr	sfp;
	
	vnpp = (ValNodePtr PNTR) gcp->userdata;
	switch (gcp->thistype)
	{
		case OBJ_SEQFEAT:
			sfp = (SeqFeatPtr) (gcp->thisitem);
			if (sfp->data.choice == SEQFEAT_GENE) {
				tmp = ValNodeNew(NULL);
				tmp->data.ptrvalue = sfp;
				*vnpp = tie_next(*vnpp, tmp);
			} 
			break;
		default:
			break;
	}
	return TRUE;
}

static ValNodePtr GatherDescrOnBioseq(ItemInfoPtr iip, BioseqPtr bsp, Uint1 choice)
{
	ValNodePtr   vnp;

	GatherScope  gsc;
	SeqLocPtr    slp = NULL;
	Uint2        bspID;
	DescrInfoPtr dsp;

	dsp = (DescrInfoPtr) MemNew(sizeof(DescrInfo));
	dsp->choice = choice;
	dsp->bsp = bsp;
  	MemSet ((Pointer) (&gsc), 0, sizeof (GatherScope));
	MemSet ((Pointer) (gsc.ignore), (int)(TRUE), (size_t) (OBJ_MAX * sizeof(Boolean)));
	gsc.ignore[OBJ_SEQDESC] = FALSE;
	bspID = ObjMgrGetEntityIDForPointer(bsp);
	slp = ValNodeNew(NULL);
	slp->choice = SEQLOC_WHOLE;
	slp->data.ptrvalue = (SeqIdPtr) SeqIdDup (SeqIdFindBest (bsp->id, 0));
	gsc.target = slp;
	if (ISA_aa(bsp->mol) && !is_pdb(bsp)) {
		GatherEntity(bspID, &dsp, get_descr, &gsc);
	} else {
		GatherEntity(bspID, &dsp, get_descr_on_top, &gsc);
	} 
	SeqLocFree(slp);
	vnp = dsp->vnp;
	if (vnp && vnp->data.ptrvalue) {
		if (iip != NULL) {
			iip->entityID = dsp->iip->entityID;
			iip->itemID = dsp->iip->itemID;
			iip->itemtype = dsp->iip->itemtype;
		}
		MemFree(dsp->iip);
		MemFree(dsp);
		return vnp;
	}
	MemFree(dsp->iip);
	MemFree(dsp);
	return NULL;
}

static SeqFeatPtr GatherSeqFeatProt(BioseqPtr bsp)
{
	GatherScope gsc;
	SeqLocPtr slp = NULL;
	Uint2 bspID;
	SeqFeatPtr sfp = NULL;
	SeqFeatPtr f;
	ValNodePtr prot, v;
	Int4 length, longest_length=0;

  	MemSet ((Pointer) (&gsc), 0, sizeof (GatherScope));
	gsc.get_feats_location = TRUE;
	bspID = ObjMgrGetEntityIDForPointer(bsp);
	slp = ValNodeNew(NULL);
	slp->choice = SEQLOC_WHOLE;
	slp->data.ptrvalue = (SeqIdPtr) SeqIdDup (SeqIdFindBest (bsp->id, 0));
	gsc.target = slp;
	prot = NULL;
	GatherEntity(bspID, &prot, GetFeatProt, &gsc);
	for (v=prot; v; v=v->next) {
		f = (SeqFeatPtr) v->data.ptrvalue;
		if ((length=SeqLocLen(f->location)) == -1)
			continue;
		if (length > longest_length) {
			sfp = f;
			longest_length = length;
		}
	}
	ValNodeFree(prot);
	SeqLocFree(slp);
	
	return sfp;
}

static SeqFeatPtr GatherProtCDS(BioseqPtr bsp)
{
	GatherScope gsc;
	SeqLocPtr slp = NULL;
	Uint2 bspID;
	SeqFeatPtr sfp;

  	MemSet ((Pointer) (&gsc), 0, sizeof (GatherScope));
	gsc.get_feats_product = TRUE;
	bspID = ObjMgrGetEntityIDForPointer(bsp);
	slp = ValNodeNew(NULL);
	slp->choice = SEQLOC_WHOLE;
	slp->data.ptrvalue = (SeqIdPtr) SeqIdDup (SeqIdFindBest (bsp->id, 0));
	gsc.target = slp;
	sfp = NULL;
	GatherEntity(bspID, &sfp, GetFeatCDS, &gsc); 
	SeqLocFree(slp);
	
	return sfp;
}
static ValNodePtr GatherGenesForCDS(SeqLocPtr slp)
{
	GatherScope gsc;
	Uint2 bspID;
	ValNodePtr vnp;
	BioseqPtr bsp;

	bsp = BioseqFindCore(SeqLocId(slp));
	if (bsp == NULL)
		return NULL;
	bspID = ObjMgrGetEntityIDForPointer(bsp);
  	MemSet ((Pointer) (&gsc), 0, sizeof (GatherScope));
	gsc.get_feats_location = TRUE;
	gsc.target = slp;
	vnp = NULL;
	GatherEntity(bspID, &vnp, GetFeatGenes, &gsc); 
	
	return vnp;
}

static CharPtr FindProtDefLine(BioseqPtr bsp)
{
	SeqFeatPtr sfp = NULL, f;
	ProtRefPtr prp;
	SeqFeatXrefPtr xref;
	GeneRefPtr grp=NULL;
	ValNodePtr vnp, v;
	SeqLocPtr loc;
	CharPtr title = NULL, s;
	Int4 diff_lowest = 0, diff_current;
	Int2 length = 0;
	SeqFeatPtr best_gene = NULL;
		
	if (bsp == NULL) {
		return NULL;
	}
	if ((sfp = GatherSeqFeatProt(bsp)) != NULL) {
		prp = (ProtRefPtr) sfp->data.value.ptrvalue;
		if (prp && prp->name) {
			for (vnp=prp->name; vnp; vnp=vnp->next) {
				length += StringLen((CharPtr)vnp->data.ptrvalue) + 2;
			}
			s = title = (CharPtr) MemNew(length + 1);
			if (prp->name->data.ptrvalue) {
				sprintf(title, "%s", 
                                        (CharPtr) prp->name->data.ptrvalue);
			}
			s += StringLen(title);
			for (vnp=prp->name->next; vnp; vnp=vnp->next) {
				sprintf(s, "; %s", 
                                        (CharPtr) vnp->data.ptrvalue);
				s += StringLen((CharPtr)vnp->data.ptrvalue) + 2;
			}
		} else if (prp->desc) {
			title = StringSave(prp->desc);
		}
	} 
	if (title == NULL) {
		if ((sfp = GatherProtCDS(bsp)) != NULL) {
			loc = sfp->location;
			for (xref = sfp->xref; xref; xref=xref->next) {
				if (xref->data.choice == SEQFEAT_GENE) {
					grp = (GeneRefPtr) xref->data.value.ptrvalue;
				}
			}
			if (grp && grp->locus) {
				title = StringSave(grp->locus);
			}
			if (title == NULL) {
				vnp = GatherGenesForCDS(loc);
				for (v=vnp; v; v=v->next) {
					f = (SeqFeatPtr) v->data.ptrvalue;
					diff_current = SeqLocAinB(loc, f->location);
					if (! diff_current) {
						best_gene = f;
						break;
					} else if (diff_current > 0) {
						if ((diff_lowest == -1) || (diff_current<diff_lowest)) {
							diff_lowest = diff_current;
							best_gene = f;
						}
					}
				}
				ValNodeFree(vnp);
				if (best_gene != NULL) {
					grp = (GeneRefPtr) best_gene->data.value.ptrvalue;
					if (grp && grp->locus) {
						s = (CharPtr) MemNew(StringLen(grp->locus) + 15);
						sprintf(s, "%s gene product", grp->locus);
						title = s;
					}
				}
			}
		}
	}
	return title;
}

static CharPtr UseOrgMods(BioseqPtr bsp)
{
	ItemInfoPtr 		iip = NULL;
	ValNodePtr 			vnp;
	BioSourcePtr 		biop;
	OrgModPtr          	mod;
	OrgNamePtr         	onp;
	OrgRefPtr          	orp;
	SubSourcePtr       	ssp;
	CharPtr				name = NULL, chr = NULL, str = NULL,
						cln = NULL, map = NULL, def=NULL;
	Int2 				deflen = 0;
		
	if (bsp == NULL) {
		return NULL;
	}
	if ((vnp=GatherDescrOnBioseq(iip, bsp, Seq_descr_source)) == NULL) {
		return NULL;
	}
	biop = (BioSourcePtr) vnp->data.ptrvalue;
    orp = biop->org;
    if (orp && orp->taxname) {
	    name = StringSave(orp->taxname);
	    deflen += StringLen(orp->taxname);
    }
	for (ssp = biop->subtype; ssp; ssp=ssp->next) {
		if (ssp->subtype == 1) { /* chromosome */
			if (ssp->name != NULL) {
				chr = (CharPtr) MemNew(StringLen(ssp->name) + 13);
				deflen += StringLen(ssp->name) + 13;
				sprintf(chr, " chromosome %s", ssp->name);		
			}
		}
		if (ssp->subtype == 3) { /* clone */
			if (ssp->name != NULL) {
				cln = (CharPtr) MemNew(StringLen(ssp->name) + 8);
				deflen += StringLen(ssp->name) + 8;
				sprintf(cln, " clone %s", ssp->name);		
			}
		}
		if (ssp->subtype == 2) { /* map */
			if (ssp->name != NULL) {
				map = (CharPtr) MemNew(StringLen(ssp->name) + 7);
				deflen += StringLen(ssp->name) + 7;
				sprintf(map, " map %s", ssp->name);		
			}
		}
	}
    if (orp != NULL) {
		onp = orp->orgname;
		if (onp != NULL) {
			mod = onp->mod;
			if (mod && mod->subtype == 2) { /* strain */
				if (mod->subname != NULL) {
					str = (CharPtr) MemNew(StringLen(mod->subname) + 9);
					deflen += StringLen(mod->subname) + 9;
					sprintf(str, " strain %s", mod->subname);		
				}
			}
		}
	}
	def = (CharPtr) MemNew(deflen+1);
	if (name) {
		def = StringCat(def, name);
		MemFree(name);
	}
	if (chr) {
		def = StringCat(def, chr);
		MemFree(chr);
	}
	if (cln) {
		def = StringCat(def, cln);
		MemFree(cln);
	}
	if (map) {
		def = StringCat(def, map);
		MemFree(map);
	}
	if (str) {
		def = StringCat(def, str);
		MemFree(str);
	}
	return def;
}

/*****************************************************************************
*
*   CreateDefLine(iip, bsp, buf, buflen, tech)
*   	Finds or makes a FASTA format defline using Gather functions
*       buf should be very long if possible
*       function truncates if buf not long enough
*       a few deflines are longer than 255
*
*		ItemInfoPtr iip is used in flat file generator to keep entityId, itemId
*		and itemtype
*****************************************************************************/
NLM_EXTERN Boolean CreateDefLine (ItemInfoPtr iip, BioseqPtr bsp, CharPtr buf, Int2 buflen, Uint1 tech, CharPtr accession, CharPtr organism)
{
	ValNodePtr vnp;
	CharPtr tmp = NULL, title = NULL;
	PdbBlockPtr pbp;
	PatentSeqIdPtr psip;
	PDBSeqIdPtr	pdbip;	
	Int2 diff, phase;
	Int4 num_segs, num_gaps;
	static Char tbuf[80];
	static CharPtr htgs[2] = {
		"unordered", "ordered" };
	Boolean htg_tech = FALSE;

	if ((bsp == NULL) || (buf == NULL) || buflen == 0) return FALSE;

	buflen--;
	buf[buflen] = '\0';
	tbuf[0] = '\0';
	if ((tech >= MI_TECH_htgs_1) && (tech <= MI_TECH_htgs_3)) {
		htg_tech = TRUE;
	}
	if (iip == NULL && accession != NULL) {
		diff = LabelCopyExtra(buf, accession, buflen, "(", ") ");
		buflen -= diff;
		buf += diff;
	}
	diff = 0;
	vnp=GatherDescrOnBioseq(iip, bsp, Seq_descr_title);
	if (vnp != NULL) {
		title = StringSave((CharPtr)vnp->data.ptrvalue);
	}
	if (tech == MI_TECH_htgs_1 || tech == MI_TECH_htgs_2) {
		if (title == NULL || *title == '\0') {
			title = UseOrgMods(bsp);
		}
		if (StringStr(title, "*** SEQUENCING IN PROGRESS *** ") == NULL) {
			diff = LabelCopy(buf, 
							"*** SEQUENCING IN PROGRESS *** ", buflen);
			buflen -= diff;
			buf += diff;
		}
	} else if (tech == MI_TECH_est || tech == MI_TECH_sts || tech == MI_TECH_survey) {
		if (title == NULL || *title == '\0') {
			title = UseOrgMods(bsp);
		}
	}
/* some titles may have zero length */
	if (title != NULL && *title != '\0') {
		diff = LabelCopy(buf, title, buflen);
		                        /* remove trailing blanks and periods */
		tmp = buf + diff - 1;   /* point at last character */
		while (tmp >= buf && ((*tmp <= ' ') || (*tmp == '.'))) {
			*tmp = '\0';
			tmp--;
			diff--;
		}
	} else if ((vnp = GatherDescrOnBioseq(iip, bsp, Seq_descr_pdb)) != NULL) {
		pbp = (PdbBlockPtr)(vnp->data.ptrvalue);
		for (vnp = bsp->id; vnp != NULL; vnp = vnp->next) {
			if (vnp->choice == SEQID_PDB) {
				pdbip = (PDBSeqIdPtr)(vnp->data.ptrvalue);
				if (pdbip && pdbip->chain > 32) {
					sprintf(tbuf, "Chain %c, ", pdbip->chain);
					diff = LabelCopy(buf, tbuf, buflen);
					buflen -= diff;
					buf += diff;
					break;
				}
			}
		}
		if (pbp && pbp->compound) {
			diff = LabelCopy(buf, (CharPtr)(pbp->compound->data.ptrvalue), buflen);
		}
	} else {
		for (vnp = bsp->id; vnp != NULL; vnp = vnp->next) {
			if (vnp->choice == SEQID_PATENT)
			{
				psip = (PatentSeqIdPtr)(vnp->data.ptrvalue);
				if (psip) {
					sprintf(tbuf, "Sequence %d from Patent %s %s",
					(int)psip->seqid, psip->cit->country, psip->cit->number);
					diff = LabelCopy(buf, tbuf, buflen);
					break;
				}
			}
		}
		if (vnp == NULL) {
			if (ISA_aa(bsp->mol)) {
				title = FindProtDefLine(bsp);
			}
			if (title != NULL) {
				diff = LabelCopy(buf, title, buflen);
			} else if (!htg_tech) {
				title = UseOrgMods(bsp);
				if (title != NULL) {
					diff = LabelCopy(buf, title, buflen);
				} else {
					diff = LabelCopy(buf, "No definition line found", buflen);
				}
			}
		}
	}
	buflen -= diff;
	buf += diff;
	
	if (htg_tech) {
		phase = (Int2)(tech - MI_TECH_htgs_1 + 1);
		if (title == NULL|| *title == '\0') {
			title = UseOrgMods(bsp);
			if (title != NULL) {
				diff = LabelCopy(buf, title, buflen);
				buflen -= diff;
				buf += diff;
			}
		}
		if (phase == 3)
		{
			if (title) {
				if (title && StringStr(title, "complete sequence") == NULL) {
					diff = LabelCopy(buf, ", complete sequence", buflen);
					buflen -= diff;
					buf += diff;
				}
			}
		} else {
			if (StringStr(title, "HTGS phase") == NULL) {
				if (diff != 0) {
					sprintf(tbuf, "; HTGS phase %d", (int)phase);
				} else {
					sprintf(tbuf, "HTGS phase %d", (int)phase);
				}
				diff = LabelCopy(buf, tbuf, buflen);
				buflen -= diff;
				buf += diff;
			}
			if (bsp->repr == Seq_repr_delta) {
				if (CountGapsInDeltaSeq(bsp, 
						&num_segs, &num_gaps, NULL, NULL, NULL, 0))
				{
					sprintf(tbuf, ", %ld %s pieces", 
							(long)(num_segs - num_gaps), htgs[phase - 1]);
					diff = LabelCopy(buf, tbuf, buflen);
					buflen -= diff;
					buf += diff;
				}
			}
		}

	} else if (tech == MI_TECH_est || tech == MI_TECH_sts || tech == MI_TECH_survey) {
		if (title == NULL|| *title == '\0') {
			title = UseOrgMods(bsp);
			if (title != NULL) {
				diff = LabelCopy(buf, title, buflen);
				buflen -= diff;
				buf += diff;
			}
		}
		if (tech == MI_TECH_est) {
			if (title) {
				if (title && StringStr(title, "mRNA sequence") == NULL) {
					diff = LabelCopy(buf, ", mRNA sequence", buflen);
					buflen -= diff;
					buf += diff;
				}
			}
		} else if (tech == MI_TECH_sts) {
			if (title) {
				if (title && StringStr(title, "sequence tagged site") == NULL) {
					diff = LabelCopy(buf, ", sequence tagged site", buflen);
					buflen -= diff;
					buf += diff;
				}
			}
		} else if (tech == MI_TECH_survey) {
			if (title) {
				if (title && StringStr(title, "genomic survey sequence") == NULL) {
					diff = LabelCopy(buf, ", genomic survey sequence", buflen);
					buflen -= diff;
					buf += diff;
				}
			}
		}
	}

	if (iip == NULL && organism != NULL) {
		LabelCopyExtra(buf, organism, buflen, " [", "]");
	}
        MemFree(title);
	return TRUE;
}

/*****************************************************************************
*
*   FastaSeqPort(bsp, is_na, do_virtual)
*   	opens a SeqPort for a fasta output of bsp
*
*****************************************************************************/
NLM_EXTERN SeqPortPtr FastaSeqPort(BioseqPtr bsp, Boolean is_na, Boolean do_virtual,
                                   Uint1 code)
{
	SeqPortPtr spp = NULL;

	if (bsp == NULL) return spp;

    spp = SeqPortNew(bsp, 0, -1, 0, code);
	if (do_virtual)
		SeqPortSet_do_virtual(spp, TRUE);
    SeqPortSeek(spp, 0, SEEK_SET);
	return spp;
}

/*****************************************************************************
*
*   FastaSeqLine(spp, buf, linelen)
*     an open seqport is passed in.
*     fills buf with linelen bases
*     assumes buf[linelen] = '\0'
*     returns FALSE when no more residues to print
*
*****************************************************************************/
NLM_EXTERN Boolean FastaSeqLine(SeqPortPtr spp, CharPtr buf, Int2 linelen, Boolean is_na)
{
	Int2 ctr = 0;
	Uint1 residue;
	Int4 pos;
	Char idbuf[40];

	if ((spp == NULL) || (buf == NULL)) return FALSE;

    while ((residue = SeqPortGetResidue(spp)) != SEQPORT_EOF)
    {
		if (! IS_residue(residue))
		{
			if (residue == INVALID_RESIDUE)
			{
				if (is_na)
					residue = 'N';
				else
					residue = 'X';
				FastaId(spp->bsp, idbuf, 39);
				pos = SeqPortTell(spp);
				ErrPostEx(SEV_ERROR,0,0, "ToFastA: Invalid residue at position %ld in %s",
					pos, idbuf);
			}
			else
			{
				if (residue == SEQPORT_VIRT)  /* gap */
				{
					if (ctr)  /* got some residues already */
					{
		    	        buf[ctr] = '\0';
						SeqPortSeek(spp, -1, SEEK_CUR);  /* back up one */
                                     /* can only seek to a real residue, so go past it */
                                                residue = SeqPortGetResidue(spp);
                                                if (residue == SEQPORT_VIRT)
                                                   SeqPortSeek(spp, -1, SEEK_CUR);
						return TRUE;
					}
					else			   /* first one */
					{
						buf[ctr] = '-';
						buf[ctr + 1] = '\0';
						return TRUE;
					}
				}

				residue = '\0';
			}

		}

		if (residue != '\0')
		{
			buf[ctr] = residue;
			ctr++;
			if (ctr == linelen)
			{
				buf[ctr] = '\0';
				return TRUE;
			}
		}
    }

    buf[ctr] = '\0';
    if (ctr)
		return TRUE;
	else
		return FALSE;
}
