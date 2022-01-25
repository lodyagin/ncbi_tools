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
*
* File Name:  medutil.c
*
* Author:  James Ostell
*   
* Version Creation Date: 8/31/93
*
* $Revision: 6.1 $
*
* File Description:  Medline Utilities for MedArch
*   Assumes user calls MedArchInit and Fini
*
* Modifications:  
* --------------------------------------------------------------------------
* Date	   Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
*
* ==========================================================================
*
*
* RCS Modification History:
* $Log: medutil.c,v $
* Revision 6.1  1998/09/14 20:54:51  bazhin
* Changes in "ten_authors()" : if successful medline lookup has returned
* no authors, then will use input ones.
*
* Revision 6.0  1997/08/25 18:35:48  madden
* Revision changed to 6.0
*
* Revision 5.0  1996/05/28 14:11:11  ostell
* Set to revision 5.0
*
 * Revision 4.4  1996/05/22  17:30:09  kans
 * changed name of error file, initialize muid in lookup procedure
 *
 * Revision 4.3  1996/03/12  22:21:07  tatiana
 * add bullet proof to MedlineToISO function
 *
 * Revision 4.2  1995/11/02  23:39:28  tatiana
 * do not accept Medline article with 0 authors
 *
 * Revision 4.1  1995/08/16  14:28:21  tatiana
 * modified get_std_auth to handle spaces between initials
 *
 * Revision 4.0  1995/07/26  13:55:12  ostell
 * force revision to 4.0
 *
 * Revision 1.32  1995/06/06  14:34:39  tatiana
 * bug fixed in print_pub()
 *
 * Revision 1.30  1995/05/25  22:17:14  kans
 * fixed FixPubEquiv
 *
 * Revision 1.29  1995/05/25  18:14:35  kans
 * additional long and int casts in ErrPostEx
 *
 * Revision 1.28  1995/05/24  16:13:47  tatiana
 * add in_press() to avoid lookup for in press articles
 *
 * Revision 1.24  1995/05/17  17:54:05  epstein
 * add RCS log revision history
 *
*/

/** for ErrPostEx() ****/

static char *this_module = "medarch";
#ifdef THIS_MODULE
#undef THIS_MODULE
#endif

#define THIS_MODULE this_module
static char *this_file = __FILE__;
#define THIS_FILE this_file

/**********************/
#include <medutil.h>		   /* the interface */
#include <medarch.h>           /* the medarch interface */
#include <mdrcherr.h>           /* the errors interface */

/*--------------------------- print_pub() ---------------------------*/

/****************************************************************************
*    print_pub   
* 											
****************************************************************************/
void print_pub(ValNodePtr pub, Boolean found, Boolean auth, Int4 muid)
{
	ValNodePtr		v, title = NULL;
	CitArtPtr		art;
	CitJourPtr		jour;
	CitBookPtr		book;
    NameStdPtr  	namestd;
    AuthorPtr  		aup;
    CharPtr			s_title=NULL, page=NULL, vol=NULL, last=NULL, first=NULL;
    Int2			year = 0;
	ImprintPtr		imp = NULL;
	DatePtr			dp;
	
	if (pub == NULL) {
		ErrPostEx(SEV_WARNING, ERR_PRINT_Failed, "Citation NULL");
    	return;
	}
	if (pub->data.ptrvalue == NULL) {
		ErrPostEx(SEV_WARNING, ERR_PRINT_Failed, "Citation NULL");
    	return;
	}
	if (pub->choice != PUB_Article)
		return;

	art = pub->data.ptrvalue;
	first = "";
	last = "";
	if (art->authors == NULL) {
		ErrPostEx(SEV_WARNING, ERR_PRINT_Failed, "Authors NULL");
	} else {
	   v = art->authors->names;
	   if (v == NULL) {
	      ErrPostEx(SEV_WARNING, ERR_PRINT_Failed, "Authors NULL");
	   } else {
	      
	      /* warning! processing only the first of the author name valnodes */
	      /* should loop through them all in case first is uninformative */
	      
	      if (art->authors->choice == 1) {
		 	aup = v->data.ptrvalue;
		 	if (aup != NULL) {
		    	namestd = aup->name->data;
		    	if (namestd->names[0]) {
		    		last = namestd->names[0];
		    	}
		    	if (namestd->names[4]) {
		    		first = namestd->names[4];
		    	}
		 	}
	      } else  {
		 	first = "";
		 	last = v->data.ptrvalue;
	      }
	   }
	}
	if (art->from == 1) {
		jour = art->fromptr;
		if (jour != NULL) {
			title = jour->title;
			imp = jour->imp;
		}
	} else if (art->from == 2) {
		book = art->fromptr;
		if (book != NULL) {
			title = book->title;
			imp = book->imp;
		}
	}
	if (title != NULL) {
		s_title = title->data.ptrvalue;
	}
	if (s_title == NULL) {
		s_title = "journal unknown";
	}
	if (imp != NULL) {
		vol = imp->volume;
		page = imp->pages;
		if (imp->date != NULL) {
			year = imp->date->data[1]+1900;
		}
	}
	if (page == NULL) {
		page = "no page number";
	}
	if (vol == NULL) {
		vol = "no volume number";
	}
	if (auth) {
		ErrPostEx(SEV_WARNING, ERR_REFERENCE_MedlineMatchIgnored,  
    	"Too many author name differences: %ld|%s %s|%s|(%d)|%s|%s", 
    				(long) muid, last, first, s_title, (int) year, vol, page);
    	return;
	}
	if (imp != NULL && imp->prepub == 2) { /* in-press */ 
		dp = DateCurr();
		if (year && (Int2) (dp->data[1]) + 1900 - year > 2) {
			ErrPostEx(SEV_WARNING, ERR_REFERENCE_OldInPress, 
    "encountered in-press article more than 2 years old: %s %s|%s|(%d)|%s|%s",
            last, first, s_title, (int) year, vol, page);
		}
		DateFree(dp);
	} else {
		if (found) {
			ErrPostEx(SEV_INFO, ERR_REFERENCE_SuccessfulMuidLookup, 
    	"%ld|%s %s|%s|(%d)|%s|%s", (long) muid, last, first, s_title, (int) year, vol, page);
    	} else if (muid == 0) {
			ErrPostEx(SEV_WARNING, ERR_REFERENCE_MuidNotFound,  
    		"%s %s|%s|(%d)|%s|%s", last, first, s_title, (int) year, vol, page);
    	}
    	else {
			ErrPostEx(SEV_WARNING, ERR_REFERENCE_MuidNotFound,  
    		">>%ld<<|%s %s|%s|(%d)|%s|%s", (long) muid, last, first, s_title, (int) year, vol, page);
    	}
    }
    return;
}

/*--------------------------- ten_authors() ---------------------------*/

/****************************************************************************
*    ten_authors   
* 											
****************************************************************************/
static Boolean ten_authors(CitArtPtr art, CitArtPtr art_tmp) 
{
	Int2 num, numnew, i, match, n;
	ValNodePtr v, pub;
	CharPtr ptr, mu[10];
	AuthorPtr	aup;
    NameStdPtr  namestd;
    Boolean		no_compare = TRUE, ret = TRUE;
		
		if (art == NULL || art_tmp == NULL) {
			return TRUE;
		}
		numnew = 0;
		num = 0;
		match = 0;
		if (art->authors != NULL)
		{
			for (v = art->authors->names; v != NULL;
									v = v->next, num++);
		}
		if (art_tmp->authors != NULL)
		{
			for (v = art_tmp->authors->names; v != NULL && numnew < 10;
								v = v->next, numnew++) {
				aup = v->data.ptrvalue;
				namestd	= aup->name->data;
				mu[numnew] = namestd->names[0];	
			}
		}

		if (art->authors != NULL && art_tmp->authors != NULL) {
			if (art->authors->choice == 1 && art_tmp->authors->choice == 1) {
				no_compare = FALSE;
				for (v = art->authors->names; v != NULL; v = v->next) {
					aup = v->data.ptrvalue;
					namestd	= aup->name->data;
					ptr = namestd->names[0];
					for (i = 0; i < numnew; i++) {
						if (StringICmp(ptr, mu[i]) == 0) {
							match++;
							break;
						}
					}	
				}
			} 
		}
		if (no_compare) {
			ret = TRUE;  /* replace */
			if (num == 0 && numnew == 0) {
					ret = TRUE;			/* replace */
			} else if (num == 0) {   /* no original authors - replace*/
				ret = TRUE;
			} else if (numnew == 0) { /* no Medline authors -  don't replace*/
				ret = FALSE;
			}
		} else {
				n = num < numnew?num:numnew;
				if (n > 3*match) { 
					ret = FALSE;
				}
		}
		if (ret && (num > 10 || numnew == 0))     /* original article has > 10 authors */
		{
			AuthListFree(art_tmp->authors);  /* keep original authors */
			art_tmp->authors = art->authors;
			art->authors = NULL;
		}
		return ret;
}
/*****************************************************************************
*
*   FindPub
*   	SeqEntryExplore to process all pubs
*
*****************************************************************************/


void FindPub(SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
	BioseqPtr bsp;
	BioseqSetPtr bssp;
	ValNodePtr vnp, pub, newpub, prev, next, head;
	SeqAnnotPtr sap;
	SeqFeatPtr sfp;
	PubdescPtr pdp;
	FindPubOptionPtr fpop;

	fpop = (FindPubOptionPtr)data;

	if (IS_Bioseq(sep))
	{
		bsp = (BioseqPtr)(sep->data.ptrvalue);
		vnp = bsp->descr;
		sap = bsp->annot;
	}
	else
	{
		bssp = (BioseqSetPtr)(sep->data.ptrvalue);
		vnp = bssp->descr;
		sap = bssp->annot;
	}

	for (vnp; vnp != NULL; vnp = vnp->next)
	{
		if (vnp->choice == Seq_descr_pub)
		{
			pdp = (PubdescPtr)(vnp->data.ptrvalue);
		 	newpub = FixPubEquiv(pdp->pub, fpop);
			pdp->pub = newpub;
		}
	}

	for (sap; sap != NULL; sap = sap->next)
	{
		if (sap->type == 1)   /* feature table */
		{
			for (sfp = (SeqFeatPtr)(sap->data); sfp != NULL; sfp = sfp->next)
			{
				if (sfp->data.choice == 6)   /* pub feature */
				{
		 			pdp = (PubdescPtr)(sfp->data.value.ptrvalue);
					newpub = FixPubEquiv(pdp->pub, fpop);
					pdp->pub = newpub;
				}
				
				if (sfp->cit != NULL)
				{
					head = sfp->cit;
					prev = NULL;
					for (pub = (ValNodePtr)(head->data.ptrvalue);
						pub != NULL; pub = next)
					{
						next = pub->next;
						newpub = FixPub(pub, fpop);
						if (prev == NULL)
							head->data.ptrvalue = newpub;
						else
							prev->next = newpub;
						prev = newpub;
					}
					head->choice = 1;   /* set to Pub-set.pub */
				}
			}
		}
	}

	return;
}

/*****************************************************************************
*
*   FixPub(pub, fpop)
*   	Tries to make any Pub into muid/cit-art
*
*****************************************************************************/
ValNodePtr FixPub (ValNodePtr pub, FindPubOptionPtr fpop)
{
	ValNodePtr newpub, tmp, v;
	Int4 		muid;
	Int2 		num, numnew;
	CitArtPtr	cit;
	

	if (pub == NULL) return NULL;

	pub->next = NULL;   /* ready for return */
	newpub = pub;
	switch (pub->choice)
	{
		case PUB_Medline:           /* medline, just split to pubequiv */
			tmp = SplitMedlineEntry((MedlineEntryPtr)(pub->data.ptrvalue));
			newpub->choice = PUB_Equiv;
			newpub->data.ptrvalue = tmp;
			break;
		case PUB_Article:
			cit = pub->data.ptrvalue;
			if (cit->from == 2 || in_press(cit)) { 
									/* article from book or in press*/
				return pub;
			}
			fpop->lookups_attempted++;
			muid = MedArchCitMatch(pub);
			if (muid)              /* matched it */
			{
				print_pub(pub, TRUE, FALSE, muid);
				fpop->lookups_succeeded++;
				if (fpop->replace_cit)
				{
					fpop->fetches_attempted++;
					tmp = FetchPub(muid);
				
					if (tmp != NULL)
					{
						if (ten_authors(pub->data.ptrvalue,
									tmp->data.ptrvalue)) 
						{

							fpop->fetches_succeeded++;
							PubFree(pub);
							pub = ValNodeNew(tmp);
							pub->choice = PUB_Muid;
							pub->data.intvalue = muid;
							newpub = ValNodeNew(NULL);
							newpub->choice = PUB_Equiv;
							newpub->data.ptrvalue = tmp;
						}
						else 
						{
							print_pub(pub, FALSE, TRUE, muid);
							newpub = pub;
							MedlineToISO(pub);
						}
					}
				}
				else
				{
					print_pub(pub, FALSE, FALSE, muid);
					newpub = pub;
					MedlineToISO(pub);
				}
			}
			break;
		case PUB_Equiv:
			tmp = FixPubEquiv((ValNodePtr)(pub->data.ptrvalue), fpop);
			pub->data.ptrvalue = tmp;
			newpub = pub;
			break;
		default:
			break;
	}
	return newpub;
}

/*****************************************************************************
*
*   SplitMedlineEntry(mep)
*      splits a medline entry into 2 pubs (1 muid, 1 Cit-art)
*      converts Cit-art to ISO/GenBank style
*      returns pointer to first (chained) pub
*      deletes original medline entry
*
*****************************************************************************/
ValNodePtr SplitMedlineEntry(MedlineEntryPtr mep)
{
	ValNodePtr uid, cit;

	uid = ValNodeNew(NULL);
	uid->choice = PUB_Muid;
	uid->data.intvalue = mep->uid;

	cit = ValNodeNew(uid);
	cit->choice = PUB_Article;
	cit->data.ptrvalue = mep->cit;
	mep->cit = NULL;
	MedlineToISO(cit);

	MedlineEntryFree(mep);

	return uid;
}

/*****************************************************************************
*
*   FixPubEquiv()
*
*****************************************************************************/
ValNodePtr FixPubEquiv(ValNodePtr pube, FindPubOptionPtr fpop)
{
	ValNodePtr tmp, newpub, muidptr=NULL, citartptr=NULL, mepptr = NULL,
		otherptr = NULL, tmp2, next, new;
	Int2 muidctr = 0, citartctr = 0, mepctr = 0, otherctr = 0;
	Int4 muid = 0, oldmuid=0;
	Int2 num, numnew;
	ValNodePtr v;
	CitArtPtr	cit;

	if (pube == NULL) return NULL;

	for (tmp = pube; tmp != NULL; tmp = next)
	{
		next = tmp->next;
		tmp->next = NULL;
		switch (tmp->choice)
		{
			case PUB_Muid:
				muidctr++;
				if (muidptr == NULL)
					muidptr = tmp;
				else
				{
					for (tmp2 = muidptr; tmp2->next != NULL; tmp2 = tmp2->next)
						continue;
					tmp2->next = tmp;
				}
				break;
			case PUB_Article:
				cit = tmp->data.ptrvalue;
				if (cit->from == 2 || in_press(cit)) {
					otherctr++;
					if (otherptr == NULL)
						otherptr = tmp;
					else
					{
						for (tmp2 = otherptr; tmp2->next != NULL; tmp2 = tmp2->next)
							continue;
						tmp2->next = tmp;
					}
				} else {
					citartctr++;
					if (citartptr == NULL)
						citartptr = tmp;
					else
					{
						for (tmp2 = citartptr; tmp2->next != NULL;
													tmp2 = tmp2->next)
							continue;
						tmp2->next = tmp;
					}
				}
				break;
			case PUB_Medline:
				mepctr++;
				if (mepptr == NULL)
					mepptr = tmp;
				else
				{
					for (tmp2 = mepptr; tmp2->next != NULL; tmp2 = tmp2->next)
						continue;
					tmp2->next = tmp;
				}
				break;
			default:
				otherctr++;
				if (otherptr == NULL)
					otherptr = tmp;
				else
				{
					for (tmp2 = otherptr; tmp2->next != NULL; tmp2 = tmp2->next)
						continue;
					tmp2->next = tmp;
				}
				break;
		}
	}

	pube = otherptr;   /* put back others */
	tmp2 = otherptr;
	if (tmp2 != NULL)
	{
		for (tmp2 = otherptr; tmp2->next != NULL; tmp2 = tmp2->next)
			continue;
	}
	
	if ((muidctr != 0) && (! fpop->always_look))   /* got an muid */
	{
		if (mepptr != NULL)
		{
			if (tmp2 == NULL)
			{
				tmp2 = mepptr;
				pube = tmp2;
			}
			else
				tmp2->next = mepptr;
			
			while (tmp2->next != NULL)
				tmp2 = tmp2->next;
		}
		if (muidptr != NULL)
		{
			if (tmp2 == NULL)
			{
				tmp2 = muidptr;
				pube = tmp2;
			}
			else
				tmp2->next = muidptr;
			
			while (tmp2->next != NULL)
				tmp2 = tmp2->next;
		}
		if (citartptr != NULL)
		{
			if (tmp2 == NULL)
			{
				tmp2 = citartptr;
				pube = tmp2;
			}
			else
				tmp2->next = citartptr;
		}
		return pube;   /* no changes */
	}

	if (mepctr != 0)    /* have a medline entry, take it first */
	{
		if (mepctr > 1)
		{
			ErrPostEx(SEV_WARNING, ERR_REFERENCE_Multiple_ref, "More than one Medline entry in Pub-equiv");
			mepptr->next = PubEquivFree(mepptr->next);
		}
		PubEquivFree(muidptr);   /* ditch others */
		PubEquivFree(citartptr);
		tmp = SplitMedlineEntry((MedlineEntryPtr)(mepptr->data.ptrvalue));
		ValNodeFree(mepptr);
		if (tmp2 == NULL)
			pube = tmp;
		else
			tmp2->next = tmp;
		return pube;
	}

	if (muidctr != 0)    /* have an muid */
	{
		oldmuid = muidptr->data.intvalue;
		if (muidctr > 1)     /* more than one, just take the first */
		{
			for (tmp = muidptr->next; tmp != NULL; tmp = tmp->next)
			{
				if (tmp->data.intvalue != muid)
				{
					ErrPostEx(SEV_WARNING, ERR_REFERENCE_Multiple_muid, "Two different muids in Pub-equiv [%ld] [%ld]\n",
						(long) oldmuid, (long) (tmp->data.intvalue));
				}
			}
		 	muidptr->next = PubEquivFree(muidptr->next);
		}
	}

	if (citartctr > 0)
	{
		if (citartctr > 1)  /* ditch extras */
		{
			ErrPostEx(SEV_WARNING, ERR_REFERENCE_Multiple_ref, "More than one Cit-art in Pub-equiv");
			citartptr->next = PubEquivFree(citartptr->next);
		}

		fpop->lookups_attempted++;
		muid = MedArchCitMatch(citartptr);
		if (muid != 0)   /* success */
		{
			print_pub(citartptr, TRUE, FALSE, muid);
			fpop->lookups_succeeded++;
			if (oldmuid)   /* already had an muid */
			{
				if (oldmuid != muid)
				{
					ErrPostEx(SEV_ERROR, ERR_REFERENCE_MuidMissmatch, "OldMUID=%ld doesn't match lookup (%ld). Keeping lookup.",
						(long)oldmuid, (long)muid);
				}
			}
			if (fpop->replace_cit)
			{
				fpop->fetches_attempted++;
				new = FetchPub(muid);
				if (new != NULL) 
				{
					fpop->fetches_succeeded++;

					if (ten_authors(citartptr->data.ptrvalue,
											new->data.ptrvalue))
					{
						if (muidptr != NULL)
							tmp = muidptr;
						else
						{
							tmp = ValNodeNew(tmp2);
							tmp->choice = PUB_Muid;
						}
						if (tmp2 == NULL) {
							pube = tmp;
						} else {
							tmp2->next = tmp;
						}
						tmp->data.intvalue = muid;
						tmp->next = new;
						
						PubFree(citartptr);
					} else {
						print_pub(citartptr, FALSE, TRUE, muid);
						PubFree(new);
						if (tmp2 == NULL) {
							pube = citartptr;
						} else {
							tmp2->next = citartptr;
						}
						if (muidptr != NULL)
							citartptr->next = muidptr;
					}
				}
			}
			else
			{
				if (muidptr != NULL)
					tmp = muidptr;
				else
				{
					tmp = ValNodeNew(tmp2);
					tmp->choice = PUB_Muid;
				}
				if (tmp2 == NULL)
					pube = tmp;
				tmp->data.intvalue = muid;
				tmp2 = tmp;
				tmp->next = citartptr;
				MedlineToISO(citartptr);
			}
			return pube;
		}
		print_pub(citartptr, FALSE, FALSE, oldmuid);
		if (tmp2 == NULL)
			pube = citartptr;
		else
			tmp2->next = citartptr;
		if (muidptr != NULL)   /* ditch the mismatched muid */
			PubFree(muidptr);
		return pube;

	}
		
	if (oldmuid != 0)    /* have an muid but no cit-art */
	{
		fpop->fetches_attempted++;
		tmp = MedArchGetPub(oldmuid);
		if (tmp == NULL)
		{
			ErrPostEx(SEV_WARNING, ERR_REFERENCE_No_reference, "Cant find article for muid [%ld]", (long) oldmuid);
		}
		else
		{
			fpop->fetches_succeeded++;
			if (fpop->replace_cit)
			{
				tmp = MedlineToISO(tmp);
				if (tmp2 == NULL)
				 	pube = tmp;
				else
					tmp2->next = tmp;
				tmp2 = tmp;
			} else {
				PubFree(tmp);
			}
		}
		if (tmp2 == NULL)
		 	pube = muidptr;
		else
			tmp2->next = muidptr;
	}

	return pube;
}


/*****************************************************************************
*
*   FetchPub(muid)
*   	get a pub given a uid
*   	convert to ISO/GenBank Style
*
*****************************************************************************/
ValNodePtr FetchPub (Int4 muid)
{
	ValNodePtr tmp;

	tmp = MedArchGetPub(muid);
	if (tmp == NULL) return NULL;
	tmp = MedlineToISO(tmp);
	return tmp;
}

/*****************************************************************************
*
*   MedlineToISO(tmp)
*   	converts a MEDLINE citation to ISO/GenBank style
*
*****************************************************************************/
ValNodePtr MedlineToISO (ValNodePtr tmp)
{
	ValNodePtr titlenode, oldnames, curr, tmp2, v;
	AuthListPtr alp;
	AuthorPtr ap;
	PersonIdPtr pip;
	NameStdPtr nsp;
	CitArtPtr cap;
	CitJourPtr cjp;
	Int4 titlenum;
	CharPtr titles[1], title;
	Int1 types[1];
	ImprintPtr ip;
	Boolean is_iso;

	if (tmp == NULL) return NULL;

	if (tmp->choice != PUB_Article) return tmp;

	cap = (CitArtPtr)(tmp->data.ptrvalue);
	alp = cap->authors;
	if (alp && alp->choice == 2)   /* ml style names */
	{
		oldnames = alp->names;
		alp->names = NULL;
		curr = NULL;
		alp->choice = 1;     /* make std names */
		for (tmp2 = oldnames; tmp2 != NULL; tmp2 = tmp2->next)
		{
			curr = get_std_auth((CharPtr)(tmp2->data.ptrvalue), ML_REF);
			if (alp->names == NULL) {
				alp->names = curr;
			} else {
				for (v = alp->names; v->next != NULL; v = v->next);
				v->next = curr;
			}
		}
		ValNodeFreeData(oldnames);
	}
	if (cap->from == 1)   /* from a journal - get iso_jta */
	{
		if ((cjp = (CitJourPtr)(cap->fromptr)) != NULL)
		{
			is_iso = FALSE;
			for (titlenode = cjp->title; titlenode != NULL; titlenode = titlenode->next)
			{
				if (titlenode->choice == Cit_title_iso_jta)  /* have it */
				{
					is_iso = TRUE;
					break;
				}
			}
			if (! is_iso)
			{
				titlenode = cjp->title;
				title = (CharPtr)(titlenode->data.ptrvalue);
				titlenum = MedArchGetTitles(titles, types, title, (Int1)titlenode->choice,
					Cit_title_iso_jta, 1);
				if (titlenum)
				{
					MemFree(title);
					titlenode->choice = types[0];
					titlenode->data.ptrvalue = titles[0];
				}
			}
			ip = cjp->imp;     /* remove Eng language */
			if (! StringCmp(ip->language, "Eng"))
			{
				ip->language = MemFree(ip->language);
			}
		}
	}
	return tmp;
}

Boolean AllUpperCase(CharPtr p )
{
	if (p == NULL) return FALSE;
    while( *p != '\0' )
	{
		if( ! IS_UPPER( *p ) )
			return FALSE;
		p++;
	}
    return TRUE;
}


void SplitMlAuthorName( CharPtr name, CharPtr last, CharPtr initials, CharPtr suffix )
{
    CharPtr	p, p2;
    Int2 i;
	Char sbuf[20], ibuf[20];

    /* Clear the ibuf field and transfer the entire name to 'last',
	excluding leading and trailing spaces */

	if (name == NULL) return;

    ibuf[0] = '\0';
    sbuf[0] = '\0';
	last[0] = '\0';
	initials[0] = '\0';
	suffix[0] = '\0';
    while(*name <= ' ')
	{
		name++;
		if (*name == '\0')
			return;
	}
    StringCpy( last, name );

    for(i=StringLen(last)-1;((i >= 0) && (last[i] <= ' ')); i--)
		last[i] = '\0';

    /* Strip off the last token (initials or name suffix (Jr, Sr, suffix.) */

	p = StringRChr(last, (int)' ');
	if (p != NULL)   /* more than just last name */
	{
	    /* Separate the token from the last name */
		
	    p2 = p + 1;
    	while((p > last) && (*p == ' '))
		{
			*p = '\0';
			*p--;
		}

	    /* If the last token is not all upper case, and there are more than
		two tokens, see if the next to the last are initials (upper case) */

    	if ( ! AllUpperCase(p2) && (p = StringRChr( last, (int)' ' )) != NULL )
	    {
		/* We have at least three tokens, is the next to last initials? */

			if( AllUpperCase( p + 1 ) )
			{
			    /* Yes - concatenate the last two tokens as initials */

	    		StringCpy( ibuf, p + 1 );
			    StringCpy( sbuf, p2 );
			    while( p > last && (*p == ' ') )
				{
					*p = '\0';
					p--;
				}
			}
	    }
		
		if (ibuf[0] == '\0') /* Only the last token goes in ibuf */
	    	StringCpy( ibuf, p2 );
	}

	     /* now add periods to ibuf and convert suffix */

	for (p = initials, p2 = ibuf; *p2 != '\0'; p2++, p++)
	{
		*p = *p2;
		if (! IS_LOWER(*(p2 + 1)))   /* watch out for foreign names */
		{
			p++;
			*p = '.';
		}
	}
	*p = '\0';

	if (sbuf[0])
	{
		if ( StringCmp(sbuf, "1d")==0)
			p = StringMove (suffix, "I.");
		else if ( StringCmp(sbuf, "2d")==0)
			p = StringMove (suffix, "II.");
		else if ( StringCmp(sbuf, "3d")==0)
			p = StringMove (suffix, "III.");
		else if ( StringCmp(sbuf, "4th")==0)
			p = StringMove (suffix, "IV.");
		else if ( StringCmp(sbuf, "5th")==0)
			p = StringMove (suffix, "V.");
		else if ( StringCmp(sbuf, "6th")==0)
			p = StringMove (suffix, "VI.");
		else if ( StringCmp(sbuf, "Sr")==0)
			p = StringMove (suffix, "Sr.");
		else if ( StringCmp(sbuf, "Jr")==0)
			p = StringMove (suffix, "Jr.");
		else
			p = StringMove (suffix, sbuf);
	}

	return;
}

ValNodePtr get_std_auth(CharPtr token, Uint1 format)
{
   static CharPtr      auth, eptr;
	Char last[80], initials[20], suffix[20];
   ValNodePtr  tmp;
   NameStdPtr  namestd = NULL;
   AuthorPtr  aup;
   PersonIdPtr pid;
   
	if (token == NULL) {
		return NULL;
	}
	eptr = token + StringLen(token) - 1;
	for (; eptr > token && *eptr == ' '; eptr--);
	namestd = NameStdNew();
	if (format == PIR_REF || format == GB_REF) {
		for(auth = token; *auth != ',' && *auth != '\0'; auth++);
		if (*auth == ',') {
			*auth = '\0'; 
			if (*(auth+1) != '\0') {
				namestd->names[4] = StringSave(auth + 1);
			}
		}
		namestd->names[0] = StringSave(token);
	} else if (format == PDB_REF) {
		for (auth = eptr; auth > token && *auth != '.'; auth--);
		if (*auth == '.') {
			if (*(auth+1) != '\0' && *(auth+1) != '.') {
				namestd->names[0] = StringSave(auth+1);
			}
			*(auth+1) = '\0';
			namestd->names[4] = StringSave(token);
		} else {
			namestd->names[0] = StringSave(token);
		}
	} else if (format == EMBL_REF || format == SP_REF) {
		for(auth = eptr; *auth != ' ' && auth > token; auth--);
		if (*auth == ' ') {
			if (*(auth-1) == '.') {
				for(auth--; *auth != ' ' && auth > token; auth--);
			}
			if (*auth == ' ') {
				*auth = '\0'; 
				namestd->names[0] = StringSave(token);
				if (*(auth+1) != '\0') {
					namestd->names[4] = StringSave(auth+1);
				}
			} else {
				namestd->names[0] = StringSave(token);
			}
		} else {
			namestd->names[0] = StringSave(token);
		}
	} else if (format == ML_REF) {
		SplitMlAuthorName(token, last, initials, suffix);
		namestd->names[0] = StringSave(last);
		if (initials[0] != '\0')
			namestd->names[4] = StringSave(initials);
		if (suffix[0] != '\0')
			namestd->names[5] = StringSave(suffix);
	}
	if (namestd != NULL && namestd->names[0] != NULL) {
		pid = PersonIdNew();
		pid->choice = 2;    /*name*/
		pid->data = namestd;
		aup = AuthorNew();
		aup->name = pid;
		tmp = ValNodeNew(NULL);
		tmp->data.ptrvalue = aup;
	} 
	
	return tmp;
}

Boolean in_press(CitArtPtr cit)
{
	ImprintPtr	imp;
	CitJourPtr	jour;
	
	if (cit == NULL) {
		return TRUE;
	}
	if ((jour = cit->fromptr) == NULL || jour->title == NULL) {
		return TRUE;
	}
	if ((imp = jour->imp) == NULL) {
		return TRUE;
	}
	if (imp->volume == NULL || imp->pages == NULL || imp->date == NULL) {
		return TRUE;
	}
	
	return FALSE;
}