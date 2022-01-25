#include <stdio.h>
#include <ncbi.h>
#include <sequtil.h>
#include <toasn3.h>
#include <toporg.h>
#include <tfuns.h>
#include <utilpub.h>

/* includes for new cleanup functions from Sequin */
#include <sqnutils.h>
#include <gather.h>

SeqDescrPtr remove_descr PROTO((SeqDescrPtr head, SeqDescrPtr x));
/****************************************************************************
* 	move org-ref, modif, mol_type, date, title and pubs
*	to seg-set level in segmented set, if not there
* 	move org-ref, modif, date, title to nuc-prot level, if not there
*	notice that mol_type and pub are not moved to nucprot level
*****************************************************************************/
void toporg(SeqEntryPtr sep)
{
	
	SeqEntryExplore(sep, NULL, ChkSegset);
	SeqEntryExplore(sep, NULL, ChkNucProt);
	
	return;
}

/****************************************************************************
* 	check for backbone entry and expand (to whole) OrgRef feature 
*	if it's not whole and the only one OrgRef in the entry
* 	whole fetures would be converted to descr later in FindOrg function
*
* 	check if org-ref, modif, mol_type, date and title
* 	are the same for all segments
* 	move them to seg-set level in segmented set, if not already there
*****************************************************************************/
void ChkSegset (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{

	BioseqSetPtr 	bssp, tmp;
	BioseqPtr		bsp;
	SeqEntryPtr	 	segsep, parts;
	ValNodePtr	 	vnp = NULL, set_vnp = NULL;
    ValNodePtr    	org, modif, mol, date, v /*, title */;
    SeqAnnotPtr		sap = NULL;
    SeqFeatPtr		tmp_sfp, sfp0, sfp;
    SeqIdPtr		sidp;
    Pointer			pnt;
    SeqLocPtr		slp;
    SeqIntPtr		sip;
    OrgRefPtr		orp;
    Boolean		 	is_org = FALSE, is_modif = FALSE, is_title = FALSE;
    Boolean		 	is_date = FALSE, is_mol = FALSE;
	Boolean 		is_na = FALSE, is_bb = FALSE, whole = FALSE;
	Int2			count = 0;
	Int4			len;
	static Char			msg[51];

	if (IS_Bioseq(sep)) {
		bsp = (BioseqPtr)(sep->data.ptrvalue);
		if ((bsp->repr != Seq_repr_raw) && (bsp->repr != Seq_repr_const))
			return;
		vnp = bsp->descr;
		sap = bsp->annot;
		len = bsp->length;
		if (bsp->mol != Seq_mol_aa) {
			is_na = TRUE;
		}
		if (bsp == NULL || !is_na) {
			return;
		}
		for (sidp = bsp->id; sidp != NULL; sidp = sidp->next) {
    		pnt = sidp->data.ptrvalue;
    		switch (sidp->choice) 
    		{
        		case SEQID_LOCAL:      /* local */
        		case SEQID_GIIM:      /* giimid */
        		case SEQID_PATENT:      /* patent seq id */
        		case SEQID_GENERAL:     /* general */
        		case SEQID_GI:     /* gi */
				case SEQID_PDB:
        		case SEQID_GENBANK:      /* genbank */
        		case SEQID_EMBL:      /* embl */
        		case SEQID_PIR:      /* pir   */
        		case SEQID_SWISSPROT:      /* swissprot */
        		case SEQID_OTHER:     /* other */
        		case SEQID_DDBJ:
				case SEQID_PRF:
            		continue;
        		case SEQID_GIBBSQ:      /* gibbseq */
        		case SEQID_GIBBMT:      /* gibbmt */
					is_bb = TRUE;
        			break;
       		 	default:
        			continue;
    		}
    	}
		if (!is_bb) {
			return;
		} 
		if (sap != NULL && sap->type == 1)   /* feature table */
		{
			tmp_sfp = (SeqFeatPtr) (sap->data);
			sfp0 = SeqFeatExtractList(&(tmp_sfp), SEQFEAT_ORG);
			for (sfp = sfp0; sfp != NULL; sfp = sfp->next) {
				orp = (OrgRefPtr)(sfp->data.value.ptrvalue);
				count++;
				if ((whole = check_whole(sfp, bsp->length)) == TRUE) {
					break;
				}
			}
			if (!whole && count == 1) {
				StringNCpy(msg, SeqLocPrint(sfp0->location), 50);
				ErrPostEx(SEV_WARNING, 0, 2,
					 "Backbone entry source with bad OrgRef feature: %s", msg);
				slp = sfp0->location;
				sip = slp->data.ptrvalue;
				sip->from = 0;
				sip->to = len-1;
			}
			tmp_sfp = tie_feat(tmp_sfp, sfp0);
			sap->data = tmp_sfp;
			bsp->annot = sap;
		}
		return;	
	}
	
	bssp = (BioseqSetPtr)(sep->data.ptrvalue);
	if (bssp->_class != 2) {    /*  do the rest for segset only */
		return;
	}
	segsep = bssp->seq_set;
	set_vnp = bssp->descr;
	if (segsep->next == NULL) {
		return;
	}
	if (!IS_Bioseq(segsep->next)) {
		tmp = (BioseqSetPtr) (segsep->next->data.ptrvalue); /*segsep->next=parts*/
		parts = tmp->seq_set;
		if (parts == NULL) {
			return;
		}
	}
/*	title = SrchSegChoice(parts, Seq_descr_title);*/
	org = SrchSegChoice(parts, Seq_descr_org);
	mol = SrchSegChoice(parts, Seq_descr_mol_type);
	modif = SrchSegChoice(parts, Seq_descr_modif);
	date = SrchSegChoice(parts, Seq_descr_update_date);
	for(v=set_vnp; v != NULL; v = v->next) {
		if (v->choice == Seq_descr_org) {
			is_org = TRUE;
		}
		if (v->choice == Seq_descr_org) {
			is_title = TRUE;
		}
		if (v->choice == Seq_descr_mol_type) {
			if (mol && mol->data.intvalue == v->data.intvalue) {
				is_mol = TRUE;
				ValNodeFree(mol);
			}
		}
		if (v->choice == Seq_descr_modif) {
			is_modif = TRUE;
		}
		if (v->choice == Seq_descr_update_date) {
			is_date = TRUE;
		}
	}
/*	
	if (!is_title) {
		set_vnp = tie_next(set_vnp, title);
	}
*/
	if (!is_modif) {
		set_vnp = tie_next(set_vnp, modif);
	}
	if (!is_org) {
		set_vnp = tie_next(set_vnp, org);
	}
	if (!is_mol) {
		set_vnp = tie_next(set_vnp, mol);
	}
	if (!is_date) {
		set_vnp = tie_next(set_vnp, date);
	}
	SrchSegSeqMol(parts);
	
	return;
}

static void RemovePubFromParts(SeqEntryPtr sep, ValNodePtr pub)
{
	BioseqPtr b;
	ValNodePtr v, vnp, next;
	SeqEntryPtr s;
	PubdescPtr pdp, p;
	
	for (vnp = pub; vnp; vnp = vnp->next) {
		pdp = vnp->data.ptrvalue;
		for (s= sep; s; s=s->next) {
			b = (BioseqPtr)(s->data.ptrvalue);
			for (v=b->descr; v; v=next) {
				next = v->next;
				if (v->choice != Seq_descr_pub)
					continue;
				p = v->data.ptrvalue;
				if (PubLabelMatch(pdp->pub, p->pub) == 0) {
					if (pdp->name != NULL || pdp->fig != NULL 
						|| pdp->num != NULL || pdp->maploc != NULL
							|| pdp->comment != NULL) {
						continue;
					} else {
						b->descr = remove_descr(b->descr, v);
					}
				}
			}
		}
	}
	return;
}
/***************************************************************************
*	0  match
*	1 no match
*	2 type unknown or not implemented
***************************************************************************/
static Int2 NumberingMatch(ValNodePtr num1, ValNodePtr num2)
{
	NumContPtr nc1, nc2;
	NumEnumPtr ne1, ne2;
	NumRealPtr nr1, nr2;
	
	if (num1 == NULL || num2 == NULL)
		return 0;
	if (num1->choice != num2->choice)
		return 1;
	switch(num1->choice)
	{
		case Numbering_cont:
			nc1 = (NumContPtr) num1->data.ptrvalue;
			nc2 = (NumContPtr) num2->data.ptrvalue;
			if (nc1->refnum == nc2->refnum) {
				return 0;
			}
			break;
		case Numbering_enum:
			ne1 = (NumEnumPtr) num1->data.ptrvalue;
			ne2 = (NumEnumPtr) num2->data.ptrvalue;
			if (ne1->num != ne2->num)
				return 1;
			return (StringCmp(ne1->buf, ne2->buf));
		case Numbering_ref_source:
			break;
		case Numbering_ref_align:
			break;
		case Numbering_real:
			nr1 = (NumRealPtr) num1->data.ptrvalue;
			nr2 = (NumRealPtr) num2->data.ptrvalue;
			return (StringCmp(nr1->units, nr2->units));
		default:
			break;
	}
	return 2;
}

static Boolean PubdescMatch(PubdescPtr p1, PubdescPtr p2)
{
	if (p1 == NULL || p2 == NULL)
		return TRUE;
	if (p1->name && p2->name) {
		if (StringCmp(p1->name, p2->name) != 0)
			return FALSE;
	}
	if (p1->fig && p2->fig) {
		if (StringCmp(p1->fig, p2->fig) != 0)
			return FALSE;
	}
	if (p1->maploc && p2->maploc) {
		if (StringCmp(p1->maploc, p2->maploc) != 0)
			return FALSE;
	}
	if (p1->comment && p2->comment) {
		if (StringCmp(p1->comment, p2->comment) != 0)
			return FALSE;
	}
	if (p1->num && p2->num) {
		if (NumberingMatch(p1->num, p2->num) != 0) 
			return FALSE;
	}
	return TRUE;
}

/* return list of pubs that are the same in all segmets */
static ValNodePtr CheckSegsForPub(SeqEntryPtr sep)
{
	BioseqPtr	bsp, b;
	ValNodePtr	vnp, v, list = NULL, vnpnext, new, next;
	PubdescPtr	pdp, p, new_p;
	Boolean		same;
	SeqEntryPtr	s;
	
	if (sep == NULL) {
		return NULL;
	}
	if (!IS_Bioseq(sep)) {
		return NULL;
	}
	bsp = (BioseqPtr)(sep->data.ptrvalue);
				       /* first bioseq from parts */
	for (vnp=bsp->descr; vnp; vnp=vnpnext) {
		vnpnext = vnp->next;
		if (vnp->choice != Seq_descr_pub)
			continue;
		pdp = vnp->data.ptrvalue;
		for (s= sep->next, same = FALSE; s; s=s->next) {
			b = (BioseqPtr)(s->data.ptrvalue);
			for (v=b->descr; v; v=next) {
				next = v->next;
				if (v->choice != Seq_descr_pub)
					continue;
				p = v->data.ptrvalue;
				if (PubLabelMatch(pdp->pub, p->pub) == 0) {
					if (PubdescMatch(pdp, p) == TRUE) {
						same = TRUE;
						break;
					}
				}
			}
			if (v == NULL) {
				same = FALSE;
				break;
			}
		}
		if (same == TRUE) {
			new = ValNodeNew(NULL);
			new->choice = Seq_descr_pub;
			new_p = AsnIoMemCopy(pdp, (AsnReadFunc) PubdescAsnRead,
											(AsnWriteFunc) PubdescAsnWrite);
			new->data.ptrvalue = new_p;
			list = tie_next(list, new);
		}
	}
		return list;
}
/* move identical pubs in segmented set to the set level */
void MoveSegmPubs (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{

	BioseqSetPtr 	bssp, tmp;
	SeqEntryPtr	 	segsep, parts;
	ValNodePtr	 	vnp = NULL;
    ValNodePtr    	v, pub, vv, next;
	PubdescPtr		pdp, pdpv;

	if (IS_Bioseq(sep)) {
		return;
	}
	bssp = (BioseqSetPtr)(sep->data.ptrvalue);
	if (bssp->_class != 2) {    /*  do the rest for segset only */
		return;
	}
	segsep = bssp->seq_set;
	if (segsep->next == NULL) {
		return;
	}
	if (!IS_Bioseq(segsep->next)) {
		tmp = (BioseqSetPtr) (segsep->next->data.ptrvalue); /*segsep->next=parts*/
		parts = tmp->seq_set;
		if (parts == NULL) {
			return;
		}
	}
/* find identical pubs in all segments */
	pub = CheckSegsForPub(parts);
	if (pub) {
		RemovePubFromParts(parts, pub);
	}
/* check if pub is already on the set descr */
	for(v=bssp->descr; v != NULL; v = v->next) {
		if (v->choice != Seq_descr_pub)
			continue;
		for (vv = pub; vv; vv = next) {
			next = vv->next;
			pdp = vv->data.ptrvalue;
			pdpv = v->data.ptrvalue;
			if (PubLabelMatch(pdp->pub, pdpv->pub) == 0) {
				PubdescFree(pdp);
				pub = remove_node(pub, vv);
			}
		}
	}
	
	bssp->descr = tie_next(bssp->descr, pub);
	
	return;
}

void ChkNucProt (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{

	BioseqSetPtr bssp, tmp;
	BioseqPtr	 bsp;
	SeqEntryPtr	 seqsep;
	ValNodePtr	 descr = NULL, vnp;
	Boolean 	 is_org = FALSE, is_modif = FALSE, is_title = FALSE;
	Boolean 	 is_date = FALSE, is_pub = FALSE;

	if (IS_Bioseq(sep)) {
		return;
	}
	bssp = (BioseqSetPtr)(sep->data.ptrvalue);
	if (bssp->_class != 1) {    /*  do the rest for nuc-prot only */
		return;
	}
	seqsep = bssp->seq_set;
	if (seqsep->choice == 1) {
		bsp = seqsep->data.ptrvalue;
		descr = bsp->descr;
	}
	if (seqsep->choice == 2) {
	      tmp = seqsep->data.ptrvalue;
	      descr = tmp->descr;
	}
	if (bssp->descr == NULL) {
		bssp->descr = GetDescr(&descr);
	} else {
		for (vnp = bssp->descr; vnp!= NULL; vnp= vnp->next) {
			if (vnp->choice == Seq_descr_title) {
				is_title = TRUE;
			}
			if (vnp->choice == Seq_descr_org) {
				is_org = TRUE;
			}
			if (vnp->choice == Seq_descr_modif) {
				is_modif = TRUE;
			}
			if (vnp->choice == Seq_descr_update_date) {
				is_date = TRUE;
			}
		/*	if (vnp->choice == Seq_descr_pub) {
				is_pub = TRUE;
			}
		*/
		}
		if (!is_title) {
			   vnp = ValNodeExtractList(&descr, Seq_descr_title);
			bssp->descr = ValNodeLink(&(bssp->descr), vnp);
		}
		if (!is_modif && check_GIBB(descr)) {
			   vnp = ValNodeExtractList(&descr, Seq_descr_modif);
			bssp->descr = ValNodeLink(&(bssp->descr), vnp);
		}
		if (!is_org) {
			   vnp = ValNodeExtractList(&descr, Seq_descr_org);
			bssp->descr = ValNodeLink(&(bssp->descr), vnp);
		}
		if (!is_date) {
			   vnp = ValNodeExtractList(&descr, Seq_descr_update_date);
			bssp->descr = ValNodeLink(&(bssp->descr), vnp);
		}
	/*	vnp = ValNodeExtractList(&descr, Seq_descr_pub);
		if (!is_pub)
			bssp->descr = ValNodeLink(&(bssp->descr), vnp);
	*/
	}
	if (seqsep->choice == 1) {
		bsp = seqsep->data.ptrvalue;
		bsp->descr = descr;
	}
	if (seqsep->choice == 2) {
	      tmp = seqsep->data.ptrvalue;
	      tmp->descr = descr;
	}
	return;
}
void MoveNPPubs (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{

	BioseqSetPtr 	bssp, tmp;
	BioseqPtr	 	bsp;
	SeqEntryPtr	 	seqsep;
	ValNodePtr		descr = NULL, vnp = NULL, vnext, v, v_copy;
	PubdescPtr		pdp, pdp_copy;

	if (IS_Bioseq(sep)) {
		return;
	}
	bssp = (BioseqSetPtr)(sep->data.ptrvalue);
	if (bssp->_class != 1) {    /*  do the rest for nuc-prot only */
		return;
	}
	seqsep = bssp->seq_set;
	if (seqsep->choice == 1) {
		bsp = seqsep->data.ptrvalue;
		descr = bsp->descr;
	}
	if (seqsep->choice == 2) {
	      tmp = seqsep->data.ptrvalue;
	      descr = tmp->descr;
	}
	if (bssp->descr == NULL) {
		bssp->descr = GetDescr(&descr);
	} else {
/* move pubs to nuc-prot level */	
		vnp = ValNodeExtractList(&descr, Seq_descr_pub);
		for (v=vnp; v; v=vnext) {
			vnext = v->next;
			pdp = (PubdescPtr) v->data.ptrvalue;
			if (pdp->num != NULL || pdp->name != NULL || pdp->fig != NULL
				|| pdp->comment !=NULL) {
				pdp_copy = AsnIoMemCopy(pdp, (AsnReadFunc) PubdescAsnRead,
						(AsnWriteFunc) PubdescAsnWrite);
				v_copy = ValNodeNew(NULL);
				v_copy->choice = Seq_descr_pub;
				v_copy->data.ptrvalue = pdp_copy;
				descr = ValNodeLink(&(descr), v_copy);
				vnp = remove_node(vnp, v);
			}
		}
		bssp->descr = ValNodeLink(&(bssp->descr), vnp);
	}
	if (seqsep->choice == 1) {
		bsp = seqsep->data.ptrvalue;
		bsp->descr = descr;
	}
	if (seqsep->choice == 2) {
	      tmp = seqsep->data.ptrvalue;
	      tmp->descr = descr;
	}
	return;
}

Boolean CmpOrgById(BioSourcePtr b1, BioSourcePtr b2)
{
	DbtagPtr d1 = NULL, d2 = NULL;
	ValNodePtr vnp;
	
	if (b1 == NULL || b2 == NULL) {
		return FALSE;
	}
	if (b1->org ==  NULL || b2->org == NULL) {
		return FALSE;
	}
	for (vnp = b1->org->db; vnp; vnp = vnp->next) {
		d1 = (DbtagPtr) vnp->data.ptrvalue;
		if (StringCmp(d1->db, "taxon") == 0) {
			break;
		}
	}
	for (vnp = b2->org->db; vnp; vnp = vnp->next) {
		d2 = (DbtagPtr) vnp->data.ptrvalue;
		if (StringCmp(d2->db, "taxon") == 0) {
			break;
		}
	}
	if (d1 && d2) {
		if (d1->tag->id == d2->tag->id) {
			return TRUE;
		} else {
		}
	} else if (StringICmp(b1->org->taxname, b2->org->taxname) == 0) {
			return TRUE;
	}
	return FALSE;
}

BioSourcePtr BioSourceMerge(BioSourcePtr host, BioSourcePtr guest)
{
	SubSourcePtr ssp, sp;
	OrgModPtr omp, homp;
	OrgNamePtr	onp;
	
	if (host == NULL && guest == NULL) {
		return NULL;
	}
	if (host == NULL && guest != NULL) {
		host = AsnIoMemCopy(guest, (AsnReadFunc) BioSourceAsnRead, 
		   						(AsnWriteFunc) BioSourceAsnWrite);
		return host;
	}
	if (host != NULL && guest == NULL) {
		return host;
	}
	if (host->genome == 0 && guest->genome != 0) {
		host->genome = guest->genome;
	}
	if (host->origin == 0 && guest->origin != 0) {
		host->origin = guest->origin;
	}
	for (ssp = guest->subtype; ssp; ssp = ssp->next) {
			sp = AsnIoMemCopy(ssp, (AsnReadFunc) SubSourceAsnRead, 
		   						(AsnWriteFunc) SubSourceAsnWrite);
		host->subtype = tie_next_subtype(host->subtype, sp);
	}
	if (guest->org->orgname) {
		for (omp = guest->org->orgname->mod; omp; omp = omp->next) {
			homp = AsnIoMemCopy(omp, (AsnReadFunc) OrgModAsnRead, 
		   						(AsnWriteFunc) OrgModAsnWrite);
		   	if ((onp = host->org->orgname)	== NULL) {
		   		onp = OrgNameNew();
		   		host->org->orgname = onp;
		   	}				
			onp->mod = tie_next_OrgMod(onp->mod, homp);
		}
	}
	return host;
}

BioSourcePtr BioSourceCommon(BioSourcePtr host, BioSourcePtr guest)
{
	SubSourcePtr ssp, sp, spnext;
	OrgModPtr omp, om, ompnext;
	
	if (host->genome != guest->genome) {
		host->genome = 0;
	}
	if (host->origin != guest->origin) {
		host->origin = 0;
	}
	for (sp = host->subtype; sp; sp = spnext) {
		spnext = sp->next;
		for (ssp = guest->subtype; ssp; ssp = ssp->next) {
			if (sp->subtype == ssp->subtype && 
					StringCmp(sp->name, ssp->name) == 0) {
				break;
			}
		}
		if (ssp == NULL) {
			host->subtype = remove_subtype(host->subtype, sp);
		}
	}
	if (CmpOrgById(host, guest) == FALSE) {
		OrgRefFree(host->org);
		host->org = NULL;
		return host;
	}
	if (StringExtCmp(host->org->common, guest->org->common) != 0) {
		MemFree(host->org->common);
		host->org->common = NULL;
	}
	if (guest->org->orgname == NULL) {
		MemFree(host->org->orgname);
		host->org->orgname = NULL;
	} else {
		if (host->org->orgname) {
			for (omp = host->org->orgname->mod; omp; omp = ompnext) {
				ompnext = omp->next;
				for (om = guest->org->orgname->mod; om; om = om->next) {
					if (om->subtype == omp->subtype && 
							StringCmp(om->subname, omp->subname) == 0) {
						break;
					}
				}
			}
			if (omp == NULL) {
				host->org->orgname->mod = 
					remove_OrgMod(host->org->orgname->mod, omp);
			}
		}
	}
	return host;
}

static Boolean EmptyBioSource(BioSourcePtr bio)
{
	OrgRefPtr org;
	
	if (bio->genome == 0 && bio->origin == 0 && bio->org == NULL)
		return TRUE;
	if ((org = bio->org) != NULL) {
		if (org->taxname == NULL && org->common == NULL && org->db == NULL)
			return TRUE;
	}
	return FALSE;
}

void StripBSfromTop (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
	BioseqSetPtr bssp, tmp;
	ValNodePtr	 vnp, bio_vnp = NULL;
	SeqEntryPtr	 	segsep, parts, cur;
	BioseqPtr bsp;
	BioSourcePtr biotop = NULL, bio = NULL;
	Boolean first = TRUE;

	if (IS_Bioseq(sep)) {
		return;
	}
	bssp = (BioseqSetPtr)(sep->data.ptrvalue);
	if (bssp->_class != 2) {    /*  do the rest for segset only */
		return;
	}
	bio_vnp = ValNodeExtractList(&(bssp->descr), Seq_descr_source);
	if (bio_vnp) {
        biotop = ((BioSourcePtr)bio_vnp->data.ptrvalue);
	} else {
		bio_vnp = ValNodeNew(NULL);
		bio_vnp->choice = Seq_descr_source;
	}
	segsep = bssp->seq_set;
	if (segsep->next == NULL) {
		return;
	}
	if (!IS_Bioseq(segsep->next)) {
		tmp = (BioseqSetPtr) (segsep->next->data.ptrvalue); /*segsep->next=parts*/
		parts = tmp->seq_set;
		if (parts == NULL) {
			return;
		}
		for (cur = parts; cur; cur = cur->next) {
			bsp = cur->data.ptrvalue;
			for (vnp = bsp->descr; vnp; vnp=vnp->next) {
				if (vnp->choice == Seq_descr_source) {
        			bio = ((BioSourcePtr)vnp->data.ptrvalue);
					break;
				}
			}
			if (bio) {
				if (biotop == NULL && first == TRUE) {
					biotop = AsnIoMemCopy(bio, (AsnReadFunc) BioSourceAsnRead, 
		   						(AsnWriteFunc) BioSourceAsnWrite);
		   			bio_vnp->data.ptrvalue = biotop;
				 	first = FALSE;
				} else {
					biotop = BioSourceCommon(biotop, bio);
				}
			}
		}
		if (biotop != NULL) {
			if (EmptyBioSource(biotop)) {
				BioSourceFree(biotop);
				ValNodeFree(bio_vnp);
			} else {
				bssp->descr = tie_next(bssp->descr, bio_vnp);
			}
		}
	}
	return;
}

void StripBSfromParts (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{

	BioseqSetPtr bssp, tmp;
	BioseqPtr	 bsp;
	SeqEntryPtr	 segsep, parts = NULL, cur;
	ValNodePtr	 descr, vnp, set_vnp;
	BSMapPtr 	bsmp;
	Boolean 	empty;
	BioSourcePtr biosp = NULL;

	bsmp = (BSMapPtr) data;
	if (IS_Bioseq(sep)) {
		return;
	}
	bssp = (BioseqSetPtr)(sep->data.ptrvalue);
	if (bssp->_class != 2) {    /*  do the rest for segset only */
		return;
	}
	vnp = ValNodeExtractList(&(bssp->descr), Seq_descr_source);
	if (vnp) {
		biosp = (BioSourcePtr) vnp->data.ptrvalue;
        ValNodeFree(vnp);
        vnp=NULL;
	}	  
	segsep = bssp->seq_set;
	if (segsep->next == NULL) {
		return;
	}
	set_vnp = bssp->descr;
	if (!IS_Bioseq(segsep->next)) {
		tmp = (BioseqSetPtr) (segsep->next->data.ptrvalue); /*segsep->next=parts*/
		parts = tmp->seq_set;
	}
	for (cur = parts, empty = TRUE; cur; cur = cur->next) {
		bsp = cur->data.ptrvalue;
		vnp = ValNodeExtractList(&(bsp->descr), Seq_descr_source);
		if (vnp) {
			biosp = BioSourceMerge(biosp, vnp->data.ptrvalue);
			BioSourceFree((BioSourcePtr) (vnp->data.ptrvalue));
			ValNodeFree(vnp);
		}
	}
	if (biosp) {
		descr = ValNodeNew(NULL);
		descr->choice = Seq_descr_source;
		descr->data.ptrvalue = biosp;
		bssp->descr = ValNodeLink(&(bssp->descr), descr);
	}

	return;
}

/*------------------------ GetDescr() --------------------------*/
/*****************************************************************************
*  GetDescr: 
*                                                                    8-12-93
******************************************************************************/
ValNodePtr GetDescr(ValNodePtr PNTR descr)
{
   ValNodePtr    vnp, hvnp = NULL;



   vnp = ValNodeExtractList(descr, Seq_descr_title);
   hvnp = ValNodeLink(&hvnp, vnp);

   vnp = ValNodeExtractList(descr, Seq_descr_org);
   hvnp = ValNodeLink(&hvnp, vnp);

   if ( check_GIBB(*descr)) {
       vnp = ValNodeExtractList(descr, Seq_descr_modif);
       hvnp = ValNodeLink(&hvnp, vnp); 
   }
   vnp = ValNodeExtractList(descr, Seq_descr_comment);
   hvnp = ValNodeLink(&hvnp, vnp);

   vnp = ValNodeExtractList(descr, Seq_descr_pub);
   hvnp = ValNodeLink(&hvnp, vnp);

   vnp = ValNodeExtractList(descr, Seq_descr_update_date);
   hvnp = ValNodeLink(&hvnp, vnp);

   return (hvnp);

} /* GetDescr */

/*------------------------ check_GIBB() --------------------------*/
/*****************************************************************************
*  check_GIBB: 
*                                                                    8-12-93
******************************************************************************/
Boolean check_GIBB(ValNodePtr descr)
{
	ValNodePtr	vnp, modif;
	Int4		gmod;
	
	if (descr == NULL) {
		return FALSE;
	}
	for (vnp = descr; vnp && vnp->choice != Seq_descr_modif; vnp = vnp->next)
		continue;
	if (vnp == NULL) {
		return FALSE;
	}
	modif = (ValNodePtr) vnp->data.ptrvalue;
	if (modif == NULL) {
		return FALSE;
	}
	gmod = modif->data.intvalue;
	if (gmod == Seq_descr_GIBB_mod_dna || gmod == Seq_descr_GIBB_mod_rna ||
    	gmod == Seq_descr_GIBB_mod_est || gmod == Seq_descr_GIBB_mod_complete
                    	|| gmod == Seq_descr_GIBB_mod_partial) {
		return FALSE;
	}
	return TRUE;
}

/*----------------------------- SrchSegChoice() --------------------------*/
/*****************************************************************************
*  SrchSegChoice:
******************************************************************************/
ValNodePtr SrchSegChoice(SeqEntryPtr sep, Uint1 choice)
{
   BioseqPtr     bsp;
   ValNodePtr    hvnp = NULL;

	if (sep == NULL) {
		return NULL;
	}
	if (IS_Bioseq(sep)) {
		bsp = (BioseqPtr)(sep->data.ptrvalue);
					       /* first bioseq from parts */
		if (CheckSegDescrChoice(sep, choice)) {     /*identical */
			hvnp = ValNodeExtractList(&(bsp->descr), choice);
			CleanUpSeqDescrChoice(sep->next, choice);
		}
	}
	return (hvnp);

} /* SrchSegChoice */

/*---------------------------- SrchSegSeqMol() --------------------------*/
/*************************************************************************
*  SrchSegSeqMol:
*                                                             5-14-93
**************************************************************************/
void SrchSegSeqMol(SeqEntryPtr sep)
{
   BioseqPtr     bsp = NULL;
   SeqEntryPtr   cursep;
   Uint1         mol;
   /*
   CharPtr       str1, str2;
   */

	if (sep == NULL || sep->next) {
		return;
	}
	if (IS_Bioseq(sep)) {
   		bsp = sep->data.ptrvalue;
	}
	if (bsp == NULL) {
		return;
	}
	mol = bsp->mol;

	for (cursep = sep->next; cursep != NULL; cursep = cursep->next) {
		if (IS_Bioseq(sep)) {
			bsp = cursep->data.ptrvalue;
		} else {
			continue;
		}
		if (mol != bsp->mol) {
          break;
		}
   }

   return;

} /* SrchSegSeqMol */

/*------------------------ CheckSegDescrChoice() -------------------------*/
/*****************************************************************************
*  CheckSegDescrChoice:
*                                                                  5-18-93
******************************************************************************/
Boolean CheckSegDescrChoice(SeqEntryPtr sep, Uint1 choice)
{
   BioseqPtr     bsp;
   SeqEntryPtr   cursep;
   ValNodePtr    vnp, mvnp;
   Boolean       same;
   Boolean       no_choice = TRUE;
   BioSourcePtr	 biosp = NULL, biosp_first = NULL;
   OrgRefPtr     orp;
   CharPtr		 title = NULL;
   DatePtr       dp = NULL;
   CharPtr       org = NULL;
   Int4          modif = -1, mol = -1;
   PubdescPtr	pdp = NULL;
 
   for (cursep = sep, same = TRUE;
                              cursep != NULL && same; cursep = cursep->next) {
       bsp = cursep->data.ptrvalue;
       for (vnp = bsp->descr; vnp != NULL && vnp->choice != choice;
                                                            vnp = vnp->next)
           continue;

       if (vnp == NULL) {
          same = FALSE;
          no_choice = TRUE;
       } else if (choice == Seq_descr_org) {
          no_choice = FALSE;
          orp = vnp->data.ptrvalue;

          if (org == NULL)
             org = orp->taxname;
          else if (StringCmp(org, orp->taxname) != 0)
             same = FALSE;
       } else if (choice == Seq_descr_source) {
          no_choice = FALSE;
          if (biosp == NULL) {
			biosp = vnp->data.ptrvalue;
          } else if (BSComparison(biosp, 
          				(BioSourcePtr) vnp->data.ptrvalue) != 0) {
          	 same = FALSE;
          }
       } else if (choice == Seq_descr_mol_type) {
          no_choice = FALSE;
          if (mol == -1)
             mol = vnp->data.intvalue;
          else if (mol != vnp->data.intvalue)
          	 same = FALSE;
       } else if (choice == Seq_descr_modif) {
          no_choice = FALSE;
          mvnp = vnp->data.ptrvalue;

          if (modif == -1)
             modif = mvnp->data.intvalue;
          else if (modif != mvnp->data.intvalue)
             same = FALSE;
       }else if (choice == Seq_descr_update_date) { 
          no_choice = FALSE;
          if (dp == NULL)
             dp = vnp->data.ptrvalue;
          else if (DateMatch(dp, vnp->data.ptrvalue, TRUE) != 0)
             same = FALSE;
       } else if (choice == Seq_descr_pub) { 
          no_choice = FALSE;
          if (pdp == NULL)
             pdp = vnp->data.ptrvalue;
          else if (PubMatch(pdp->pub, 
          	((PubdescPtr)(vnp->data.ptrvalue))->pub) != 0)
             same = FALSE;
       } else if (choice == Seq_descr_title) { 
          no_choice = FALSE;
          if (title == NULL)
             title = vnp->data.ptrvalue;
          else if (StringCmp(title, (CharPtr) vnp->data.ptrvalue) != 0)
             same = FALSE;
       } else {
          	no_choice = FALSE;
				ErrPostEx(SEV_WARNING, 0, 2,
					 "Unrecognized choice: %d", choice);
       		same = FALSE;
       }
   }
	if (same == FALSE && no_choice == TRUE) {
		same = TRUE;
	}
   return (same);

} /* CheckSegDescrChoice */


void StripProtXref (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
	BioseqPtr       bsp, prot = NULL;
	BioseqSetPtr    bssp;
	ValNodePtr		vnp = NULL;
	SeqAnnotPtr     sap = NULL, ap, pap = NULL;
	SeqFeatPtr		sfp, psfp = NULL, head;
	ProtRefPtr		prp, pprp = NULL;
	SeqFeatXrefPtr	xrp, xrpnext;
	SeqIdPtr		sid;
	SeqLocPtr		slp;

	if (IS_Bioseq(sep)) {
		bsp = (BioseqPtr)(sep->data.ptrvalue);
		sap = bsp->annot;
	}
	else {
		bssp = (BioseqSetPtr)(sep->data.ptrvalue);
		sap = bssp->annot;
	}
	for (ap = sap; ap != NULL; ap = ap->next) {
		if (ap->type != 1) {
			continue;
		}
		head = (SeqFeatPtr)(ap->data);
		for (sfp = head; sfp; sfp = sfp->next) {
     		if (sfp->data.choice != SEQFEAT_CDREGION) {
     			return;
     		}
     		if ((vnp = sfp->product) != NULL) {
     			if (vnp->choice == SEQLOC_WHOLE) {
     				sid = vnp->data.ptrvalue;
     				prot = BioseqFind(sid);
     			}
     		}
     		if (prot != NULL) {
				pap = prot->annot;
				if (pap != NULL) {
					for (psfp = pap->data; psfp; psfp=psfp->next) {
						if (psfp->data.choice == SEQFEAT_PROT) {
     						pprp = psfp->data.value.ptrvalue;
     						break;
     					}
     				}
				}
			}
			if (vnp) {     	/* sfp->product != NULL */	
     			for (xrp = sfp->xref; xrp != NULL; xrp = xrpnext) {
     				xrpnext = xrp->next;
     				if (xrp->data.choice == SEQFEAT_PROT) {
     					prp = xrp->data.value.ptrvalue;
     					if (pap != NULL && pprp == NULL) {
     						if (psfp == NULL) {
     							psfp = SeqFeatNew();
     							psfp->data.choice = SEQFEAT_PROT;
     							slp = ValNodeNew(NULL);
     							slp->choice = SEQLOC_WHOLE;
     							slp->data.ptrvalue = SeqIdDup(sid);
     							psfp->location = slp;
     							pap->data = tie_feat(pap->data, psfp);
     						}
     						psfp->data.value.ptrvalue = AsnIoMemCopy(prp, 
								(AsnReadFunc) ProtRefAsnRead, 
		   						(AsnWriteFunc) ProtRefAsnWrite);
     					}
     					sfp->xref = remove_xref(sfp->xref, xrp);
     				}
				}
     		}
     	}
     	ap->data = head;
    }
	if (IS_Bioseq(sep)) {
		bsp = (BioseqPtr)(sep->data.ptrvalue);
		bsp->annot = sap;
	}
	else {
		bssp = (BioseqSetPtr)(sep->data.ptrvalue);
		bssp->annot = sap;
	}
	return;
}

static SeqLocPtr GetAnticodonFromObject(SeqFeatPtr sfp)
{
	UserObjectPtr usop;
	UserFieldPtr ufp;
	Int4Ptr ints;
	SeqLocPtr 	slp;
	SeqIntPtr 	sip;
	Int4		from, to;

	if (sfp == NULL) {
		return NULL;
	}
	if ((usop = sfp->ext) == NULL) {
		return NULL;
	}
	ufp = usop->data;
	if (ufp && ufp->choice == 8) {  /* ints */
		ints = (Int4Ptr) ufp->data.ptrvalue;
		from = ints[0];
		to = ints[1];
	}
	sip = SeqIntNew();
	sip->from = from;
	sip->to = to;
	sip->id = SeqIdDup(SeqLocId(sfp->location));
	slp = ValNodeNew(NULL);
	slp->choice = SEQLOC_INT;
	slp->data.ptrvalue = sip;
	MemFree(usop);
	return slp;

}

/*--------------------------- CheckMaps() --------------------------*/
/***************************************************************************
*   CheckMaps:
*   -- find all /map and Gene-ref
*		if all maps are the same put it to Biosource and remove quals
*   -- change User-Object anticodon in tRNA to SeqLoc
****************************************************************************/
void CheckMaps (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
	BioseqPtr       bsp;
	BioseqSetPtr    bssp;
	ValNodePtr		descr;
	SeqAnnotPtr     sap = NULL, ap;
	CharPtr			qval;
	SeqFeatPtr		sfp;
	GeneRefPtr		grp;
	QualMapPtr		qmp;
   	RnaRefPtr  		rrp;
	tRNAPtr			trna;
	GBQualPtr		q, qnext;
		
	qmp = data;
	if (qmp->same == FALSE) {
		return;
	}
	if (IS_Bioseq(sep)) {
		bsp = (BioseqPtr)(sep->data.ptrvalue);
		if ((bsp->repr != Seq_repr_raw) && (bsp->repr != Seq_repr_const))
			return;
		descr = bsp->descr;
		sap = bsp->annot;
	}
	else {
		bssp = (BioseqSetPtr)(sep->data.ptrvalue);
		descr = bssp->descr;
		sap = bssp->annot;
	}
	
/* look for all the same maploc and place it to SubSource*/
	for (ap = sap; ap != NULL; ap = ap->next) {
		if (ap->type != 1) {
			continue;
		}
		for (sfp = (SeqFeatPtr)(ap->data); sfp; sfp = sfp->next) {
     		if ((qval = get_qvalue(sfp->qual, "map")) != NULL) {
			 	if (qmp->name == NULL) {
				 	qmp->name = StringSave(qval);
			 	} else if (StringCmp(qmp->name, qval) != 0) {
			 		qmp->same = FALSE;
				 	break;
			 	}
			}
			if (sfp->data.choice == SEQFEAT_GENE) {
				grp = sfp->data.value.ptrvalue;
				if ((qval = grp->maploc) != NULL) {
					if (qmp->name == NULL) {
					 	qmp->name = StringSave(qval);
				 	} else if (StringCmp(qmp->name, qval) != 0) {
				 		qmp->same = FALSE;
					 	break;
				 	}
			 	}
     		}
			if (sfp->data.choice == SEQFEAT_RNA) {
				rrp = sfp->data.value.ptrvalue;	
   				if (rrp->type == 3 && rrp->ext.choice == 2) {
   					trna = rrp->ext.value.ptrvalue;
   					if (sfp->ext != NULL && trna->anticodon == NULL) {
   						trna->anticodon = GetAnticodonFromObject(sfp);
   						sfp->ext = NULL;
   						for (q = sfp->qual; q; q = qnext) {
   							qnext = q->next;
   							if (StringCmp(q->qual, "anticodon") == 0) {
   								sfp->qual = remove_qual(sfp->qual, q);
   							}
   						}
   					}
   				}
			}
		}
	}
	return;
}

void StripMaps (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
	BioseqPtr       bsp;
	BioseqSetPtr    bssp;
	BioSourcePtr	biosp = NULL;
	SubSourcePtr	ssp;
	ValNodePtr		descr = NULL, vnp;
	SeqAnnotPtr     sap = NULL, ap;
	CharPtr			qval = NULL;
	SeqFeatPtr		sfp;
	GeneRefPtr		grp;
	QualMapPtr		qmp;
	SeqFeatXrefPtr	xrp;

	qmp = data;
	if (qmp->same == FALSE || qmp->name == NULL) {
		return;
	}
	if (IS_Bioseq(sep)) {
		bsp = (BioseqPtr)(sep->data.ptrvalue);
		if ((bsp->repr != Seq_repr_raw) && (bsp->repr != Seq_repr_const))
			return;
		descr = bsp->descr;
		sap = bsp->annot;
	}
	else {
		bssp = (BioseqSetPtr)(sep->data.ptrvalue);
		descr = bssp->descr;
		sap = bssp->annot;
	}
	for ( vnp = descr; vnp != NULL; vnp = vnp->next) {
		if (vnp->choice == Seq_descr_source) {
			biosp = vnp->data.ptrvalue;
			break;
		}
	}
	for (ap = sap; ap != NULL; ap = ap->next) {
	 if (ap->type == 1) {
		for (sfp = (SeqFeatPtr)(ap->data); sfp; sfp = sfp->next) {
     		qval = qvalue_extract(&(sfp->qual), "map");
		if(qval != NULL)
		{
			MemFree(qval);
			qval = NULL;
		}
     		if (sfp->data.choice == SEQFEAT_GENE) {
     			grp = sfp->data.value.ptrvalue;
			if(grp->maploc != NULL)
			{
				MemFree(grp->maploc);
     				grp->maploc = NULL;
			}
     		}
     		for (xrp = sfp->xref; xrp != NULL; xrp = xrp->next) {
     			if (xrp->data.choice == SEQFEAT_GENE) {
     				grp = xrp->data.value.ptrvalue;
				if(grp->maploc != NULL)
				{
					MemFree(grp->maploc);
	     				grp->maploc = NULL;
				}
     			}
     		}
     	}
      }
    }
	if (biosp != NULL) {  /* has biosource */
			ssp = SubSourceNew();
			ssp->subtype = 2;    /*map */
			ssp->name = StringSave(qmp->name);
			biosp->subtype = tie_next_subtype(biosp->subtype, ssp);	
    }
/*    if (qval) {
		MemFree(qval);
	}*/
	return;
}

static Boolean GBQualPresent(CharPtr ptr, GBQualPtr gbqual)

{
	Boolean present=FALSE;
	GBQualPtr qual;

	for (qual=gbqual; qual; qual=qual->next)
		if (StringCmp(ptr, qual->qual) == 0)
		{
			present = TRUE;
			break;
		}

	return present;
}

void MapsToGenref (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
	BioseqPtr       bsp;
	BioseqSetPtr    bssp;
	BioSourcePtr	biosp = NULL;
	ValNodePtr		descr = NULL;
	SeqAnnotPtr     sap = NULL, ap;
	CharPtr			qval= NULL, name;
	SeqFeatPtr		sfp, cur;
	SeqLocPtr		loc;
	GeneRefPtr		grp;
	
	if (IS_Bioseq(sep)) {
		bsp = (BioseqPtr)(sep->data.ptrvalue);
		if ((bsp->repr != Seq_repr_raw) && (bsp->repr != Seq_repr_const))
			return;
		descr = bsp->descr;
		sap = bsp->annot;
	}
	else {
		bssp = (BioseqSetPtr)(sep->data.ptrvalue);
		descr = bssp->descr;
		sap = bssp->annot;
	}

/* maps are different */
  for (ap = sap; ap != NULL; ap = ap->next) {
	 if (ap->type == 1) {
       for (sfp = (SeqFeatPtr)(ap->data); sfp; sfp = sfp->next) {
    	if (sfp->data.choice == SEQFEAT_GENE) {
    		grp = sfp->data.value.ptrvalue;
    		name = NULL;
    		loc = NULL;
    		for (cur = (SeqFeatPtr)(sap->data); cur; cur = cur->next) {
    			if ((GBQualPresent("map", cur->qual)) == FALSE) {
    				continue;
    			}
    			if (SeqLocAinB(cur->location, sfp->location) < 0 ) {
    				continue;
    			}
    			qval = qvalue_extract(&(cur->qual), "map");
    			if (name == NULL) {
    				name = qval;
    				loc = cur->location;
    			} else {
    				if (SeqLocLen(cur->location) >= SeqLocLen(loc)) {
    					name = qval;
    					loc = cur->location;
    				}
    			}
    		} 
    		if (grp->maploc == NULL && name != NULL) {
    			grp->maploc = StringSave(name);
    			MemFree(qval);
    		}
    	} /* if SEQFEAT_GENE */
	  }
  	} /* if ftable */
  }
    return;
}

static Boolean CheckMinPub(ValNodePtr pub)
{
	CitGenPtr	gen;

	if (pub == NULL) {
		return TRUE;
	}
	if (pub->choice == PUB_Muid) {
		if (pub->next == NULL) {
			return TRUE;
		} else {
			return (CheckMinPub(pub->next));
		}
	}
	if (pub->choice == PUB_Gen) {
		gen = pub->data.ptrvalue;
		if (gen->cit != NULL && gen->journal == NULL && gen->authors == NULL
			&& gen->volume == NULL && gen->pages == NULL)
			return TRUE;
	}
	return FALSE;
}

ValNodePtr AddToList(ValNodePtr list, ValNodePtr check, PubdescPtr pdp)
{
	ValNodePtr	v, vnext;
	PubdescPtr	vpdp;
	PubStructPtr psp;
	ValNodePtr pubequ1 = NULL, pubequ2 = NULL;
	Boolean is_1;
		
	if (pdp == NULL) {
		return NULL;
	}
	for (v = check; v != NULL; v = v->next) {
		psp = v->data.ptrvalue;
		if (psp->start != 2) {
			continue;
		} 
		if (PubLabelMatch(psp->pub, pdp->pub) == 0) {
			return list;
		}
	}
	if (pdp->name == NULL && pdp->fig == NULL && pdp->fig == NULL) {
		if (CheckMinPub(pdp->pub) == TRUE) {   /* do not add minimum pub */
			return list;
		}
	}
	for (v = list; v != NULL; v = vnext) {
		vnext = v->next;
		vpdp = v->data.ptrvalue;
		if (vpdp->pub->next != NULL) {
			pubequ1 = ValNodeNew(NULL);
			is_1 = TRUE;
			pubequ1->choice = PUB_Equiv;
			pubequ1->data.ptrvalue = vpdp->pub;
		} else {
			is_1 = FALSE;
			pubequ1 = vpdp->pub;
		}
		if (pdp->pub->next != NULL) {
			pubequ2 = ValNodeNew(NULL);
			pubequ2->choice = PUB_Equiv;
			pubequ2->data.ptrvalue = pdp->pub;
		} else {
			pubequ2 = pdp->pub;
		}
		if (PubLabelMatch(pubequ1, pubequ2) == 0) {
			if (pdp->reftype == 2 && vpdp->reftype == 1) {
				vpdp->reftype = 2;
			}
			if (pdp->reftype == 1 && vpdp->reftype == 2) {
				pdp->reftype = 2;
			}
			if (SelectBestPub(pubequ1, pubequ2) >= 0) {
				if (is_1) {
					ValNodeFree(pubequ1);
				}
				if (pdp->pub->next != NULL) {
					ValNodeFree(pubequ2);
				}
				return list;
			}
			PubdescFree((PubdescPtr) (v->data.ptrvalue));
			list = remove_node(list, v);
		}
		if (is_1) {
			ValNodeFree(pubequ1);
		}
		if (pdp->pub->next != NULL) {
			ValNodeFree(pubequ2);
		}
	}
	if (pdp && pdp->pub) {
		v = ValNodeNew(NULL);
		v->choice = Seq_descr_pub;
		v->data.ptrvalue = AsnIoMemCopy(pdp, (AsnReadFunc) PubdescAsnRead, 
									   (AsnWriteFunc) PubdescAsnWrite);
	} 
	list = ValNodeLink(&list, v);
/*  may be sort ???? */
	return list;
}
void CheckCitSubNew(ValNodePtr vnp)
{
	CitSubPtr csp;
	AuthListPtr  alp;
	ImprintPtr   imp;
	
	if (vnp == NULL)
		return;
	if (vnp->choice != PUB_Sub)
		return;
	csp = (CitSubPtr) vnp->data.ptrvalue;
    if (csp != NULL) {
		alp = csp->authors;
		imp = csp->imp;
		if (alp != NULL && alp->affil == NULL &&
				imp != NULL && imp->pub != NULL) {
			alp->affil = imp->pub;
			imp->pub = NULL;
		}
		if (csp->date == NULL && imp != NULL && imp->date != NULL) {
			csp->date = imp->date;
			imp->date = NULL;
		}
		if (imp != NULL && imp->pub == NULL) {
			csp->imp = ImprintFree (csp->imp);
		}
	}
	return;
}

void ChangeCitSub (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
	BioseqPtr       bsp = NULL;
	BioseqSetPtr    bssp = NULL;
	ValNodePtr		descr = NULL, vnp, v;
	PubdescPtr 		pdp;
	
	if (IS_Bioseq(sep)) {
		bsp = (BioseqPtr)(sep->data.ptrvalue);
		if ((bsp->repr != Seq_repr_raw) && (bsp->repr != Seq_repr_const))
			return;
		descr = bsp->descr;
	}
	else {
		bssp = (BioseqSetPtr)(sep->data.ptrvalue);
		descr = bssp->descr;
	}
	for (vnp = descr; vnp; vnp=vnp->next) {
		if (vnp->choice != Seq_descr_pub) {
			continue;
		}
		if ((pdp = vnp->data.ptrvalue) == NULL) {
			continue;
		}
		for (v = pdp->pub; v; v=v->next) {
			if (v->choice == PUB_Sub) {
				CheckCitSubNew(v);
			}
		}
	}
}

/***************************************************************************
*   NewPubs:
*   -- find all ImpFeat "sites"
*		change to pubdesc with reftype 'sites'
*	-- find all other sfp->cit
*		change to pubdesc with reftype 'feats'
*	-- pubs are moved from SeqAnnot to Seqdescr on the same level
****************************************************************************/
void NewPubs (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
	BioseqPtr       bsp = NULL;
	BioseqSetPtr    bssp = NULL;
	ValNodePtr		descr = NULL, pubset, tmp, pubequ;
	ValNodePtr		next_pubequ, pub, min_pub;
	SeqAnnotPtr     sap = NULL, ap, apnext;
	SeqFeatPtr		sfp, cur, curnext;
	ImpFeatPtr		ifp;
	PubdescPtr		pubdesc;
	ValNodePtr		publist = NULL, check = NULL, np_list = NULL;
	
	if (IS_Bioseq(sep)) {
		bsp = (BioseqPtr)(sep->data.ptrvalue);
		if ((bsp->repr != Seq_repr_raw) && (bsp->repr != Seq_repr_const))
			return;
		descr = bsp->descr;
		sap = bsp->annot;
	}
	else {
		bssp = (BioseqSetPtr)(sep->data.ptrvalue);
		descr = bssp->descr;
		sap = bssp->annot;
	}
   	tmp = ValNodeExtractList(&descr, Seq_descr_pub);
	np_list = ValNodeLink(&np_list, tmp);
	for (ap = sap; ap != NULL; ap = apnext) {
		apnext = ap->next;
		if (ap->type != 1) {
			continue;
		}
		sfp = (SeqFeatPtr)(ap->data);
	 	for (cur = sfp; cur; cur = curnext) {
       		curnext = cur->next;
       		if (cur->cit == NULL) {
       			continue;
       		}
			pubset = cur->cit;
			pub = NULL;
			min_pub = NULL;
			pubequ = pubset->data.ptrvalue; 
			while (pubequ) {
				next_pubequ = pubequ->next; 
				pubdesc = PubdescNew();
				if (pubequ->choice == PUB_Equiv) {
					pubdesc->pub = pubequ->data.ptrvalue;
				} else {
					pubdesc->pub = pubequ;
				}
				if (cur->data.choice == SEQFEAT_IMP) {
				    ifp = cur->data.value.ptrvalue;
				    if (StringCmp(ifp->key, "Site-ref") == 0) {
					    pubdesc->reftype = 1; /* sites */
					    np_list = AddToList(np_list, check, pubdesc);
						min_pub = MinimizePub(pubequ);
						pub = tie_next(pub, min_pub);
						MemFree(pubdesc);
				     } else {
						pubdesc->reftype = 2;
						np_list = AddToList(np_list, check, pubdesc);
						min_pub = MinimizePub(pubequ);
						pub = tie_next(pub, min_pub);
						MemFree(pubdesc);
				     }
				} else {
					pubdesc->reftype = 2;
					np_list = AddToList(np_list, check, pubdesc);
					min_pub = MinimizePub(pubequ);
					pub = tie_next(pub, min_pub);
					MemFree(pubdesc);
				 }
				 PubFree(pubequ);
				 pubequ = next_pubequ;
			 }
			 if (pub && pubset) {
				 pubset->data.ptrvalue = pub;
			 }
		}
		ap->data = sfp;
	}
	
	for (ap = sap; ap != NULL; ap = apnext) {
		apnext = ap->next;
		if (ap->data == NULL) {
			sap = remove_annot(sap, ap);
		}
	}
	
	if (bssp != NULL) {
		descr = tie_next(descr, np_list);
		bssp->descr = descr;
		bssp->annot = sap;
	} else {
		descr = tie_next(descr, np_list);
		bsp->descr = descr;
		bsp->annot = sap;
	}
	return;
}

void CmpPub (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
	BioseqPtr bsp;
	PubdescPtr pdp, vpdp;
	ValNodePtr	pubequ1, pubequ2, v;
	PubdescPtr PNTR ppdp;
	Boolean is_1;
	
	ppdp = data;
	pdp = *ppdp;
	if (pdp == NULL) {
		return;
	}
	if (!IS_Bioseq(sep)) {
		return;
	}
	bsp = (BioseqPtr)(sep->data.ptrvalue);
	if (bsp->mol == Seq_mol_aa) {
		return;
	}
	for(v = bsp->descr; v; v=v->next) {
		if (v->choice != Seq_descr_pub) {
			continue;
		}
		vpdp = v->data.ptrvalue;
		if (vpdp->pub->next != NULL) {
			pubequ1 = ValNodeNew(NULL);
			is_1 = TRUE;
			pubequ1->choice = PUB_Equiv;
			pubequ1->data.ptrvalue = vpdp->pub;
		} else {
			is_1 = FALSE;
			pubequ1 = vpdp->pub;
		}
		if (pdp->pub->next != NULL) {
			pubequ2 = ValNodeNew(NULL);
			pubequ2->choice = PUB_Equiv;
			pubequ2->data.ptrvalue = pdp->pub;
		} else {
			pubequ2 = pdp->pub;
		}
		if (PubMatch(pubequ1, pubequ2) == 0) {
			if (is_1) {
				ValNodeFree(pubequ1);
			}
			if (pdp->pub->next != NULL) {
				ValNodeFree(pubequ2);
			}
			break;
		}
		if (is_1) {
			ValNodeFree(pubequ1);
		}
		if (pdp->pub->next != NULL) {
			ValNodeFree(pubequ2);
		}
	}
	if (v == NULL) {
		PubdescFree(*ppdp);
		*ppdp = NULL;
	}

	return;
}

/***********************************************************************
*	delete pubs from Bioseqs if they are already on the top level
*	don't delete Pubdesc if additional info (name fig, num etc) is present
************************************************************************/
void DeletePubs (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
	BioseqPtr bsp;
	PubdescPtr pdp, vpdp;
	ValNodePtr	pubequ1, pubequ2, v, vnext, descr = NULL;
	Boolean is_1, is_2;
	
	
	pdp = data;
	if (pdp == NULL) {
		return;
	}
	if (!IS_Bioseq(sep)) {
		return;
	}
	bsp = (BioseqPtr)(sep->data.ptrvalue);
	if ((bsp->repr != Seq_repr_raw) && (bsp->repr != Seq_repr_const))
		return;
	/*
	if (bsp->mol == Seq_mol_aa) {
		return;
	}
	*/
	descr = bsp->descr;
	for(v = descr; v; v = vnext) {
		vnext = v->next;
		if (v->choice != Seq_descr_pub) {
			continue;
		}
		vpdp = v->data.ptrvalue;   /* from the Bioseq */
		if (vpdp->name != NULL || vpdp->fig != NULL 
						|| vpdp->num != NULL || vpdp->maploc != NULL
							|| vpdp->comment != NULL) {
			continue;
		}
		if (vpdp->pub->next != NULL) {
			is_1 = TRUE;
			pubequ1 = ValNodeNew(NULL);
			pubequ1->choice = PUB_Equiv;
			pubequ1->data.ptrvalue = vpdp->pub;
		} else {
			is_1 = FALSE;
			pubequ1 = vpdp->pub;
		}
		if (pdp->pub->next != NULL) {  /* from the set */
			is_2 = TRUE;
			pubequ2 = ValNodeNew(NULL);
			pubequ2->choice = PUB_Equiv;
			pubequ2->data.ptrvalue = pdp->pub;
		} else {
			is_2 = FALSE;
			pubequ2 = pdp->pub;
		}
		if (PubMatch(pubequ1, pubequ2) == 0) {
			PubdescFree((PubdescPtr) (v->data.ptrvalue));
			descr = remove_node(descr, v);
		}
		if (is_1) {
			ValNodeFree(pubequ1);
		}
		if (is_2) {
			ValNodeFree(pubequ2);
		}
	}
	bsp->descr = descr;
	return;
}

void MoveSetPubs (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
	BioseqPtr 		bsp = NULL;
	BioseqSetPtr	bssp;
	Boolean			first;
	PubSetListPtr	psp;
	ValNodePtr 		descr = NULL,v, vnext, tmp, set_list;
	PubdescPtr		tmp_pdp, pdp;
	
	psp = data;
	set_list = psp->list;
	first = psp->first;
	if (IS_Bioseq(sep) && (first == TRUE)) {
		bsp = sep->data.ptrvalue;
		if ((bsp->repr != Seq_repr_raw) && (bsp->repr != Seq_repr_const))
			return;
		descr = bsp->descr;
	} else {
		bssp = sep->data.ptrvalue;
		if (bssp->_class == 4 && bssp->seq_set != NULL) {
			bsp = (BioseqPtr) bssp->seq_set->data.ptrvalue;
			descr = bsp->descr;
		}
	}
	if (bsp == NULL) {
		return;
	}
	if (first == FALSE) {
		return;
	}
	for (v = descr; v; v=vnext) {
		vnext = v->next;
		if (v->choice == Seq_descr_pub) {
			pdp = v->data.ptrvalue;
			tmp_pdp = AsnIoMemCopy(pdp, 
				(AsnReadFunc) PubdescAsnRead, (AsnWriteFunc) PubdescAsnWrite);
			SeqEntryExplore(sep, &tmp_pdp, CmpPub);
			if (tmp_pdp != NULL) {
				tmp = ValNodeNew(NULL);
				tmp->choice = Seq_descr_pub;
				tmp->data.ptrvalue = tmp_pdp;
				set_list = tie_next(set_list, tmp);
				SeqEntryExplore(sep, tmp_pdp, DeletePubs);
				first = FALSE;
			}
		}
	}
	psp->list = set_list;
	psp->first = first;
	data = psp;
	return;
}

static void AddFeat (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
	BioseqPtr       bsp ;
	PubdescPtr		pdp, tmp_pdp;
	SeqAnnotPtr     ap;
	ValNodePtr		pubset, tmp;
	SeqFeatPtr		sfp;
	ImpFeatPtr		ifp;
	SeqIdPtr		sip;
	
	pdp = data;
	if (!IS_Bioseq(sep)) {
		return;
	}
	bsp = (BioseqPtr)(sep->data.ptrvalue);
	if (bsp->mol == Seq_mol_aa) {
		return;
	}
	for (ap = bsp->annot; ap != NULL; ap = ap->next) {
		if (ap->type != 1) {
			continue;
		}
		sfp = SeqFeatNew();
		sfp->data.choice = SEQFEAT_IMP;
	    ifp = ImpFeatNew();

	    ifp->key = StringSave("Site-ref");
	    ifp->loc = StringSave("sites");
	    sfp->data.value.ptrvalue = ifp;

	    sfp->location = ValNodeNew(NULL);
	    sfp->location->choice = SEQLOC_WHOLE;
		sip = SeqIdDup(bsp->id);
	    sfp->location->data.ptrvalue = sip ;
	    pubset = ValNodeNew(NULL);
	    pubset->choice = 1;       
		tmp = ValNodeNew(NULL);
		tmp->choice = PUB_Equiv;
		tmp_pdp = AsnIoMemCopy(pdp, 
		(AsnReadFunc) PubdescAsnRead, (AsnWriteFunc) PubdescAsnWrite);
		tmp->data.ptrvalue = tmp_pdp->pub;
		pubset->data.ptrvalue = tmp;
		sfp->cit = pubset;
		ap->data = tie_feat(ap->data, sfp);
		
	}
	return;
}

void FindOldLineage (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
	BioseqPtr       bsp;
	BioseqSetPtr    bssp;
	ValNodePtr		descr = NULL, vnp;
	GBBlockPtr		gb;
	CharPtr	PNTR	linp;
	CharPtr			lineage;
	
	linp = (CharPtr PNTR) data;
	lineage = *linp;
	if (IS_Bioseq(sep)) {
		bsp = (BioseqPtr)(sep->data.ptrvalue);
		if ((bsp->repr != Seq_repr_raw) && (bsp->repr != Seq_repr_const))
			return;
		descr = bsp->descr;
	} else {
		bssp = (BioseqSetPtr)(sep->data.ptrvalue);
		descr = bssp->descr;
	}
	for (vnp = descr; vnp; vnp= vnp->next) {
		if (vnp->choice == Seq_descr_genbank) {
			gb = vnp->data.ptrvalue;
			if (gb->taxonomy) {
				if (*linp) {
					MemFree(*linp);
				}
				*linp = gb->taxonomy;
				gb->taxonomy = NULL;
			}
			break;
		}
	}
}

void FindNewLineage (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
	BioseqPtr       bsp;
	BioseqSetPtr    bssp;
	ValNodePtr		descr = NULL, vnp;
	BioSourcePtr	biosp;
	OrgNamePtr		onp;
	OrgRefPtr		orp;
	CharPtr	PNTR	linp;
	CharPtr			lineage;
	
	linp = (CharPtr PNTR) data;
	lineage = *linp;
	if (lineage != NULL && *lineage != '\0') {
		return;
	}
	
	if (IS_Bioseq(sep)) {
		bsp = (BioseqPtr)(sep->data.ptrvalue);
		if ((bsp->repr != Seq_repr_raw) && (bsp->repr != Seq_repr_const))
			return;
		descr = bsp->descr;
	}
	else {
		bssp = (BioseqSetPtr)(sep->data.ptrvalue);
		descr = bssp->descr;
	}
	for (vnp = descr; vnp; vnp= vnp->next) {
		if (vnp->choice == Seq_descr_source) {
			biosp = vnp->data.ptrvalue;
			orp = (OrgRefPtr) biosp->org;
			if (orp && orp->orgname) {
				onp = orp->orgname;
				*linp = onp->lineage;
				orp->orgname->lineage = NULL;
				orp->orgname->gcode = 0;
				orp->orgname->mgcode = 0;
				if (onp->lineage == NULL && onp->mod == NULL &&
					onp->gcode == 0 && onp->mgcode == 0 &&
					onp->data == NULL && onp->attrib == NULL) {
					orp->orgname = OrgNameFree (orp->orgname);
				}
			}
			break;
		}
	}
}

void NewLineage (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
	BioseqPtr       bsp;
	BioseqSetPtr    bssp;
	ValNodePtr		descr = NULL, vnp;
	BioSourcePtr	biosp;
	OrgRefPtr		orp = NULL;
	OrgNamePtr		omp;
	CharPtr	PNTR	linp;
	CharPtr			lineage;

	linp = (CharPtr PNTR) data;
	lineage = *linp;
	if (IS_Bioseq(sep)) {
		bsp = (BioseqPtr)(sep->data.ptrvalue);
		if ((bsp->repr != Seq_repr_raw) && (bsp->repr != Seq_repr_const))
			return;
		descr = bsp->descr;
	}
	else {
		bssp = (BioseqSetPtr)(sep->data.ptrvalue);
		descr = bssp->descr;
	}
	for (vnp = descr; vnp; vnp= vnp->next) {
		if (vnp->choice == Seq_descr_source) {
			biosp = vnp->data.ptrvalue;
			orp = (OrgRefPtr) biosp->org;
			break;
		}
	}
	if (orp && lineage) {
			if (orp->orgname == NULL) {
				omp = OrgNameNew();
				orp->orgname = omp;
			}
			if(orp->orgname->lineage != NULL)
			{
				MemFree(orp->orgname->lineage);
			}
			orp->orgname->lineage = StringSave(lineage);
	}
	
}

void OldLineage (SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent)
{
	BioseqPtr       bsp;
	BioseqSetPtr    bssp;
	ValNodePtr		descr = NULL, vnp;
	GBBlockPtr		gb = NULL;
	CharPtr	PNTR	linp;
	CharPtr			lineage;

	linp = (CharPtr PNTR) data;
	lineage = *linp;
	if (IS_Bioseq(sep)) {
		bsp = (BioseqPtr)(sep->data.ptrvalue);
		if ((bsp->repr != Seq_repr_raw) && (bsp->repr != Seq_repr_const))
			return;
		descr = bsp->descr;
	}
	else {
		bssp = (BioseqSetPtr)(sep->data.ptrvalue);
		descr = bssp->descr;
	}
	for (vnp = descr; vnp; vnp= vnp->next) {
		if (vnp->choice == Seq_descr_genbank) {
			gb = vnp->data.ptrvalue;
			break;
		}
	}
	if (gb && lineage) {
			gb->taxonomy = StringSave(lineage);
			MemFree(lineage);
	}
}

void OldPubs (SeqEntryPtr sep)
{
	BioseqPtr       bsp = NULL;
	BioseqSetPtr    bssp;
	ValNodePtr		descr = NULL, cur, first, pre= NULL, next_f;
	PubdescPtr		pdp;
		
	if (IS_Bioseq(sep)) {
		bsp = (BioseqPtr)(sep->data.ptrvalue);
		if ((bsp->repr != Seq_repr_raw) && (bsp->repr != Seq_repr_const))
			return;
		descr = bsp->descr;
	}
	else {
		bssp = (BioseqSetPtr)(sep->data.ptrvalue);
		descr = bssp->descr;
	}
	 for (first = cur = descr; cur; cur = next_f) {
       	next_f = cur->next;
		if (cur->choice != Seq_descr_pub) {
			continue;
		}
		pdp = cur->data.ptrvalue;
		if (pdp->reftype == 1 || pdp->reftype == 2) {
/* create ImpFeat "Sites-ref" for every nuc bioseqs in the SeqEntry */
			SeqEntryExplore(sep, pdp, AddFeat);
					/* delete valnode from descr */
			DelNodeFromList(&descr, cur, pre);
		} 
		if (pre) {
			if (pre->next != next_f) {
				pre = cur;
			}
		} else {
			if (first == descr) {
				pre = cur;
			} else {
				first = descr;
			}
		}
		if (descr == NULL) {
			break;
		}
	}
	if (bsp) {
		bsp->descr = descr;
	} else {
		bssp->descr = descr;
	}
}	

/****************************************************************************
*  delete_valnode:
*****************************************************************************/
static ValNodePtr delete_valnode(ValNodePtr host, Uint1 choice)
{
   Boolean       first;
   ValNodePtr    curvnp, prevnp;
   
   for (curvnp = host, first = TRUE; curvnp != NULL
	     && curvnp->choice != choice; curvnp = curvnp->next) {

       if (first) {
	  prevnp = curvnp;
	  first = FALSE;
       }
       else
	  prevnp = prevnp->next;
   }

   if (curvnp == NULL) {
		return host;
   }
   if (first) {
   		host = curvnp->next;
   } else {
      prevnp->next = curvnp->next;
   }
  curvnp->next = NULL;

  switch (choice) {
  case Seq_descr_org:
     OrgRefFree(curvnp->data.ptrvalue);
     break;
  case Seq_descr_modif:
     ValNodeFree(curvnp->data.ptrvalue);
     break;
  case Seq_descr_update_date:
     DateFree(curvnp->data.ptrvalue);
     break;
  case Seq_descr_mol_type:
  	 break;
  default:
  	 break;
  }

  ValNodeFree(curvnp);
  
  return host;
}

/*------------------- CleanUpSeqDescrChoice() -------------------------*/
/****************************************************************************
*  CleanUpSeqDescrChoice:
*                                                                    5-21-93
*****************************************************************************/
void CleanUpSeqDescrChoice(SeqEntryPtr sep, Uint1 choice)
{
   BioseqPtr     bsp;
   SeqEntryPtr   cursep;
  
   for (cursep = sep; cursep != NULL; cursep = cursep->next) {
       bsp = cursep->data.ptrvalue;

	   bsp->descr = delete_valnode(bsp->descr, choice);
   }

} /* CleanUpSeqDescrChoice */

/**********************************************************/
SeqDescrPtr remove_descr(SeqDescrPtr head, SeqDescrPtr x)
{
    SeqDescrPtr v;
    SeqDescrPtr p;
	
    if(head == NULL)
        return(NULL);

    if(x == head)
    {
        head = x->next;
        x->next = NULL;
        SeqDescFree(x);
        return(head);
    }
    for(v = head; v != NULL && v != x; v = v->next)
        p = v;

    if(v != NULL)
    {
        p->next = x->next;
        x->next = NULL;
        SeqDescFree(x);
    }
    return(head);
}


/* Cleanup functions originally from Sequin */

static void FindConsistentMolInfo (SeqEntryPtr sep, MolInfoPtr PNTR mipp, BoolPtr consist)

{
  BioseqPtr     bsp = NULL;
  BioseqSetPtr  bssp = NULL;
  MolInfoPtr    mip;
  ValNodePtr    sdp = NULL;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sdp = bsp->descr;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sdp = bssp->descr;
  } else return;
  while (sdp != NULL) {
    if (sdp->choice == Seq_descr_molinfo) {
      mip = (MolInfoPtr) sdp->data.ptrvalue;
      if (mip != NULL) {
        if (*mipp == NULL) {
          *mipp = mip;
        } else {
          if ((*mipp)->biomol != mip->biomol ||
              (*mipp)->tech != mip->tech ||
              (*mipp)->completeness != mip->completeness ||
              StringICmp ((*mipp)->techexp, mip->techexp) != 0) {
            *consist = FALSE;
          }
        }
      }
    }
    sdp = sdp->next;
  }
  if (bssp == NULL) return;
  for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
    FindConsistentMolInfo (sep, mipp, consist);
  }
}

static void RemoveMolInfoCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  ValNodePtr    nextsdp;
  Pointer PNTR  prevsdp;
  ValNodePtr    sdp;

  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sdp = bsp->descr;
    prevsdp = (Pointer PNTR) &(bsp->descr);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sdp = bssp->descr;
    prevsdp = (Pointer PNTR) &(bssp->descr);
  } else return;
  while (sdp != NULL) {
    nextsdp = sdp->next;
    if (sdp->choice == Seq_descr_molinfo) {
      *(prevsdp) = sdp->next;
      sdp->next = NULL;
      SeqDescFree (sdp);
    } else {
      prevsdp = (Pointer PNTR) &(sdp->next);
    }
    sdp = nextsdp;
  }
}

extern void NormalizeSegSeqMolInfo (SeqEntryPtr sep)

{
  BioseqSetPtr  bssp;
  Boolean       consistent;
  MolInfoPtr    master;
  MolInfoPtr    mip;
  ValNodePtr    sdp;

  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;
    if (bssp != NULL && bssp->_class != BioseqseqSet_class_segset) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        NormalizeSegSeqMolInfo (sep);
      }
      return;
    }
    if (bssp != NULL && bssp->_class == BioseqseqSet_class_segset) {
      mip = NULL;
      consistent = TRUE;
      FindConsistentMolInfo (sep, &mip, &consistent);
      if (mip != NULL && consistent) {
        master = MolInfoNew ();
        if (master == NULL) return;
        master->biomol = mip->biomol;
        master->tech = mip->tech;
        master->completeness = mip->completeness;
        master->techexp = StringSaveNoNull (mip->techexp);
        SeqEntryExplore (sep, NULL, RemoveMolInfoCallback);
        sdp = CreateNewDescriptor (sep, Seq_descr_molinfo);
        if (sdp != NULL) {
          sdp->data.ptrvalue = (Pointer) master;
        }
      }
    }
  }
}

static void CollectPseudoCdsProducts (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr        bsp;
  BioseqSetPtr     bssp;
  BioseqPtr        product;
  Boolean          pseudo;
  GBQualPtr        gbqual;
  SeqAnnotPtr      sap;
  SeqFeatPtr       sfp;
  ValNodePtr PNTR  vnpp;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  vnpp = (ValNodePtr PNTR) mydata;
  if (vnpp == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
  } else return;
  while (sap != NULL) {
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      while (sfp != NULL) {
        if (sfp->data.choice == SEQFEAT_CDREGION) {
          pseudo = sfp->pseudo;
          if (! pseudo) {
            gbqual = sfp->qual;
            while (gbqual != NULL) {
              if (StringICmp (gbqual->qual, "pseudo") == 0) {
                pseudo = TRUE;
              }
              gbqual = gbqual->next;
            }
          }
          if (pseudo) {
            product = BioseqFind (SeqLocId (sfp->product));
            if (product != NULL) {
              ValNodeAddPointer (vnpp, 0, (Pointer) product);
              sfp->product = SeqLocFree (sfp->product);
            }
          }
        }
        sfp = sfp->next;
      }
    }
    sap = sap->next;
  }
}

extern void CleanUpPseudoProducts (Uint2 entityID, SeqEntryPtr sep)

{
  BioseqPtr      bsp;
  Uint2          itemID;
  ValNodePtr     list;
  OMProcControl  ompc;
  ValNodePtr     vnp;

  if (entityID == 0 || sep == NULL) return;
  list = NULL;
  SeqEntryExplore (sep, &list, CollectPseudoCdsProducts);
  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    bsp = (BioseqPtr) vnp->data.ptrvalue;
    itemID = GetItemIDGivenPointer (entityID, OBJ_BIOSEQ, (Pointer) bsp);
    if (itemID > 0) {
      MemSet ((Pointer) (&ompc), 0, sizeof (OMProcControl));
      ompc.do_not_reload_from_cache = TRUE;
      ompc.input_entityID = entityID;
      ompc.input_itemID = itemID;
      ompc.input_itemtype = OBJ_BIOSEQ;
      if (! DetachDataForProc (&ompc, FALSE)) {
        Message (MSG_POSTERR, "DetachDataForProc failed");
      }
    }
  }
  ValNodeFree (list);
}

extern void CleanupGenbankCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  Boolean        empty;
  GBBlockPtr     gbp;
  ValNodePtr     nextsdp;
  Pointer PNTR   prevsdp;
  ValNodePtr     sdp;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sdp = bsp->descr;
    prevsdp = (Pointer PNTR) &(bsp->descr);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sdp = bssp->descr;
    prevsdp = (Pointer PNTR) &(bssp->descr);
  } else return;
  while (sdp != NULL) {
    nextsdp = sdp->next;
    empty = FALSE;
    if (sdp->choice == Seq_descr_genbank && sdp->data.ptrvalue != NULL) {
      gbp = (GBBlockPtr) sdp->data.ptrvalue;
      gbp->source = MemFree (gbp->source);
      gbp->origin = MemFree (gbp->origin);
      gbp->taxonomy = MemFree (gbp->taxonomy);
      if (gbp->extra_accessions == NULL && gbp->source == NULL &&
          gbp->keywords == NULL && gbp->origin == NULL &&
          gbp->date == NULL && gbp->entry_date == NULL &&
          gbp->div == NULL && gbp->taxonomy == NULL) {
        empty = TRUE;
      }
    }
    if (empty) {
      *(prevsdp) = sdp->next;
      sdp->next = NULL;
      SeqDescFree (sdp);
    } else {
      prevsdp = (Pointer PNTR) &(sdp->next);
    }
    sdp = nextsdp;
  }
}

static Boolean EmptyOrNullString (CharPtr str)

{
  Char  ch;

  if (str == NULL) return TRUE;
  ch = *str;
  while (ch != '\0') {
    if (ch > ' ' && ch <= '~') return FALSE;
    str++;
    ch = *str;
  }
  return TRUE;
}

extern void MergeAdjacentAnnotsCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  SeqAnnotPtr   nextsap;
  SeqAnnotPtr   sap;
  SeqFeatPtr    sfp;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
  } else return;
  while (sap != NULL) {
    nextsap = sap->next;
    if (sap->type == 1 && nextsap != NULL && nextsap->type == 1) {
      if (sap->id == NULL && nextsap->id == NULL &&
          sap->name == NULL && nextsap->name == NULL &&
          sap->db == 0 && nextsap->db == 0 &&
          sap->desc == NULL && nextsap->desc == NULL &&
          sap->data != NULL && nextsap->data != NULL) {
        sfp = (SeqFeatPtr) sap->data;
        while (sfp->next != NULL) {
          sfp = sfp->next;
        }
        sfp->next = (SeqFeatPtr) nextsap->data;
        nextsap->data = NULL;
        sap->next = nextsap->next;
        SeqAnnotFree (nextsap);
        nextsap = sap->next;
      }
    }
    sap = nextsap;
  }
}

extern void CleanupEmptyFeatCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  Boolean        empty;
  GeneRefPtr     grp;
  SeqAnnotPtr    nextsap;
  SeqFeatPtr     nextsfp;
  Pointer PNTR   prevsap;
  Pointer PNTR   prevsfp;
  ProtRefPtr     prp;
  SeqAnnotPtr    sap;
  SeqFeatPtr     sfp;
  ValNodePtr     vnp;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
    prevsap = (Pointer PNTR) &(bsp->annot);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
    prevsap = (Pointer PNTR) &(bssp->annot);
  } else return;
  while (sap != NULL) {
    nextsap = sap->next;
    if (sap->type == 1) {
      sfp = (SeqFeatPtr) sap->data;
      prevsfp = (Pointer PNTR) &(sap->data);
      while (sfp != NULL) {
        nextsfp = sfp->next;
        empty = FALSE;
        if (sfp->data.choice == SEQFEAT_GENE && sfp->data.value.ptrvalue != NULL) {
          grp = (GeneRefPtr) sfp->data.value.ptrvalue;
          if (EmptyOrNullString (grp->locus)) {
            grp->locus = MemFree (grp->locus);
          }
          if (EmptyOrNullString (grp->allele)) {
            grp->allele = MemFree (grp->allele);
          }
          if (EmptyOrNullString (grp->desc)) {
            grp->desc = MemFree (grp->desc);
          }
          if (EmptyOrNullString (grp->maploc)) {
            grp->maploc = MemFree (grp->maploc);
          }
          if (EmptyOrNullString (grp->locus) &&
			  EmptyOrNullString (grp->allele) &&
			  EmptyOrNullString (grp->desc) &&
			  EmptyOrNullString (grp->maploc) &&
			  grp->db == NULL && grp->syn == NULL) {
            empty = TRUE;
          }
        } else if (sfp->data.choice == SEQFEAT_PROT && sfp->data.value.ptrvalue != NULL) {
          prp = (ProtRefPtr) sfp->data.value.ptrvalue;
          if (prp->processed != 3 && prp->processed != 4) {
            vnp = prp->name;
            if ((vnp == NULL || EmptyOrNullString ((CharPtr) vnp->data.ptrvalue)) &&
                EmptyOrNullString (prp->desc) &&
                prp->ec == NULL && prp->activity == NULL && prp->db == NULL) {
              empty = TRUE;
			}
          }
        } else if (sfp->data.choice == SEQFEAT_COMMENT && EmptyOrNullString (sfp->comment)) {
          empty = TRUE;
        }
        if (empty) {
          *(prevsfp) = sfp->next;
          sfp->next = NULL;
          SeqFeatFree (sfp);
        } else {
          prevsfp = (Pointer PNTR) &(sfp->next);
        }
        sfp = nextsfp;
      }
    }
    if (sap->data == NULL) {
      *(prevsap) = sap->next;
      sap->next = NULL;
      SeqAnnotFree (sap);
    } else {
      prevsap = (Pointer PNTR) &(sap->next);
    }
    sap = nextsap;
  }
}

extern void RemoveBioSourceOnPopSet (SeqEntryPtr sep, OrgRefPtr master)

{
  BioSourcePtr  biop;
  BioseqSetPtr  bssp;
  OrgRefPtr     orp;
  ValNodePtr    sdp;

  if (sep == NULL) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;
    if (bssp->_class == 7 || bssp->_class == 13 ||
        bssp->_class == 14 || bssp->_class == 15) { /* now on phy and mut sets */
      sdp = SeqEntryGetSeqDescr (sep, Seq_descr_source, NULL);
      if (sdp == NULL) return;
      biop = (BioSourcePtr) sdp->data.ptrvalue;
      if (biop == NULL) return;
      orp = biop->org;
      if (orp == NULL) return;
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        RemoveBioSourceOnPopSet (sep, orp);
      }
      sdp = ValNodeExtract (&(bssp->descr), Seq_descr_source);
      SeqDescrFree (sdp);
      return;
    }
    /* if (bssp->_class == 7 || bssp->_class == 13 || bssp->_class == 15) return; */
    if (bssp->_class == 7) { /* also handle genbank supersets */
      orp = NULL;
      sdp = SeqEntryGetSeqDescr (sep, Seq_descr_source, NULL);
      if (sdp != NULL) {
        biop = (BioSourcePtr) sdp->data.ptrvalue;
        if (biop != NULL) {
          orp = biop->org;
        }
      }
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        RemoveBioSourceOnPopSet (sep, orp);
      }
      sdp = ValNodeExtract (&(bssp->descr), Seq_descr_source);
      SeqDescrFree (sdp);
      return;
    }
  }
  if (master == NULL) return;
  sdp = SeqEntryGetSeqDescr (sep, Seq_descr_source, NULL);
  if (sdp != NULL) return;
  biop = BioSourceNew ();
  if (biop == NULL) return;
  orp = OrgRefNew ();
  if (orp == NULL) return;
  biop->org = orp;
  orp->taxname = StringSave (master->taxname);
  orp->common = StringSave (master->common);
  sdp = CreateNewDescriptor (sep, Seq_descr_source);
  if (sdp == NULL) return;
  sdp->data.ptrvalue = (Pointer) biop;
}

extern Boolean NoBiosourceOrTaxonId (SeqEntryPtr sep)

{
  BioSourcePtr  biop;
  BioseqSetPtr  bssp;
  DbtagPtr      dbt;
  ObjectIdPtr   oid;
  OrgRefPtr     orp;
  ValNodePtr    sdp;
  ValNodePtr    vnp;

  if (sep == NULL) return TRUE;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && (bssp->_class == 7 || bssp->_class == 13 ||
                         bssp->_class == 14 || bssp->_class == 15)) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        if (NoBiosourceOrTaxonId (sep)) return TRUE;
      }
      return FALSE;
    }
  }
  sdp = SeqEntryGetSeqDescr (sep, Seq_descr_source, NULL);
  if (sdp == NULL) return TRUE;
  biop = (BioSourcePtr) sdp->data.ptrvalue;
  if (biop == NULL) return TRUE;
  orp = biop->org;
  if (orp == NULL) return TRUE;
  vnp = orp->db;
  if (vnp == NULL) return TRUE;
  while (vnp != NULL) {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt != NULL) {
      if (StringCmp (dbt->db, "taxon") == 0) {
        oid = dbt->tag;
        if (oid != NULL) {
          if (oid->str == NULL && oid->id > 0) {
            return FALSE;
          }
        }
      }
    }
    vnp = vnp->next;
  }
  return TRUE;
}

static void CollectGeneFeatures (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr        bsp;
  BioseqSetPtr     bssp;
  SeqAnnotPtr      sap;
  SeqFeatPtr       sfp;
  ValNodePtr PNTR  vnpp;

  if (sep == NULL || sep->data.ptrvalue == NULL || mydata == NULL) return;
  vnpp = (ValNodePtr PNTR) mydata;
  sap = NULL;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
  } else return;
  while (sap != NULL) {
    if (sap->type == 1 && sap->data != NULL) {
      sfp = (SeqFeatPtr) sap->data;
      while (sfp != NULL) {
        if (sfp->data.choice == SEQFEAT_GENE) {
          ValNodeAddPointer (vnpp, 0, (Pointer) sfp);
        }
        sfp = sfp->next;
      }
    }
    sap = sap->next;
  }
}

static void ExtendGeneWithinNucProt (SeqEntryPtr sep)

{
  BioseqSetPtr  bssp;
  ValNodePtr    vnp;

  if (sep == NULL) return;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp == NULL) return;
    if (bssp->_class == 7 || bssp->_class == 13 ||
        bssp->_class == 14 || bssp->_class == 15) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        ExtendGeneWithinNucProt (sep);
      }
      return;
    }
  }
  vnp = NULL;
  SeqEntryExplore (sep, (Pointer) (&vnp), CollectGeneFeatures);
  if (vnp != NULL && vnp->next == NULL) {
    SeqEntryExplore (sep, NULL, CorrectGeneFeatLocation);
  }
  ValNodeFree (vnp);
}

extern void ExtendGeneFeatIfOnMRNA (Uint2 entityID, SeqEntryPtr sep)

{
  if (entityID < 1 && sep == NULL) return;
  if (entityID > 0 && sep == NULL) {
    sep = GetTopSeqEntryForEntityID (entityID);
  }
  if (sep == NULL) return;
  ExtendGeneWithinNucProt (sep);
}

static Int4 LoopSeqEntryToAsn3 (SeqEntryPtr sep, Boolean strip, Boolean correct, SeqEntryFunc taxfun, SeqEntryFunc taxmerge)

{
  BioseqSetPtr  bssp;
  SeqEntryPtr   oldscope;
  Int4          rsult;
  Boolean       taxserver;

  rsult = 0;
  if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    if (bssp != NULL && (bssp->_class == 7 || bssp->_class == 13 ||
                         bssp->_class == 14 || bssp->_class == 15)) {
      for (sep = bssp->seq_set; sep != NULL; sep = sep->next) {
        rsult += LoopSeqEntryToAsn3 (sep, strip, correct, taxfun, taxmerge);
      }
      return rsult;
    }
  }
  oldscope = SeqEntrySetScope (sep);
  taxserver = (Boolean) (taxfun != NULL || taxmerge != NULL);
  rsult = SeqEntryToAsn3Ex (sep, strip, correct, taxserver, taxfun, taxmerge);
  SeqEntrySetScope (oldscope);
  return rsult;
}

extern void SeriousSeqEntryCleanup (SeqEntryPtr sep, SeqEntryFunc taxfun, SeqEntryFunc taxmerge)

{
  Uint2   entityID;
  ErrSev  lsev;
  ErrSev  msev;

  if (sep == NULL) return;
  msev = ErrSetMessageLevel (SEV_MAX);
  lsev = ErrSetLogLevel (SEV_MAX);
  entityID = SeqMgrGetEntityIDForSeqEntry (sep);
  BasicSeqEntryCleanup (sep);
  SeqEntryExplore (sep, NULL, CleanupGenbankCallback);
  SeqEntryExplore (sep, NULL, CleanupEmptyFeatCallback);
  SeqEntryExplore (sep, NULL, MergeAdjacentAnnotsCallback);
  EntryChangeImpFeat(sep);     /* change any CDS ImpFeat to real CdRegion */
  MoveRnaGBQualProductToName (sep); /* move rna gbqual product to rna-ref.ext.name */
  MoveProtGBQualProductToName (sep); /* move prot gbqual product to prot-ref.name */
  MoveCdsGBQualProductToName (sep); /* move cds gbqual product to prot-ref.name */
  MoveFeatGBQualsToFields (sep); /* move feature partial, exception to fields */
  ExtendGeneFeatIfOnMRNA (0, sep); /* gene on mRNA is full length */
  RemoveBioSourceOnPopSet (sep, NULL);
  SeqEntryExplore (sep, NULL, DeleteMultipleTitles);
  LoopSeqEntryToAsn3 (sep, TRUE, FALSE, taxfun, taxmerge);
  EntryStripSerialNumber(sep); /* strip citation serial numbers */
  EntryChangeGBSource (sep);   /* at least remove redundant information in GBBlocks */
  EntryCheckGBBlock (sep);
  SeqEntryMoveDbxrefs (sep); /* db_xref gbqual to sfp->dbxref */
  EntryMergeDupBioSources (sep);
  SeqEntryExplore (sep, NULL, GetRidOfEmptyFeatsDescCallback);
  CleanUpPseudoProducts (entityID, sep);
  RenormalizeNucProtSets (sep, TRUE);
  StripTitleFromProtsInNucProts (sep);
  move_cds (sep);
  CdCheck (sep, NULL);
  ErrSetMessageLevel (msev);
  ErrSetLogLevel (lsev);
}

