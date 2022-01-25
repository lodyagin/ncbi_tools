/*----------------*/

#include <stdlib.h>
#include <ncbi.h>
#include <taxinc.h>
#include <txclient.h>
#define REALTAXsyb

#define MAX_ORG_LIST	10

#define BUFF_SIZE 16
#define TAX_READ 0
#define TAX_WRITE 1


static TreePtr tax_tree= NULL;

static Int2 VRL_div= 0;
static Int2 PHG_div= 0;

static Int2 SpeciesRank= 26;
static Int2 SubspeciesRank= 27;
static Int2 GenusRank= 22;
static Int2 FamilyRank= 0;
static Int2 OrderRank= 0;
static Int2 ClassRank= 0;
static Int2 SYNONYM= 0;
static Int2 COMMON_NAME= 0;
static Int2 PREF_COMMON= 0;

static int my_timer= 0;

static struct t_or_buff {
  Int4 tax_id;
  OrgRefPtr p_org_ref;
  int timer;
  char div[16];
  char embl[4];
  int is_species;
} or_buff[BUFF_SIZE];

static CharPtr DB_PATH= "TAXCLIENT_OS";

static Boolean we_want_synonyms= 0;

static OrgRefPtr getFromBuff(Int4 id, int* is_sp, CharPtr div, CharPtr embl);
static void loadInBuff(Int4 id);
static void bldOrgRefOut(OrgRefPtr dst, OrgRefPtr src, Int4 tax_id);

Boolean tax1_setSynonyms(Boolean on_off)
{
    Boolean ret;

    ret= we_want_synonyms;
    we_want_synonyms= on_off;
    return ret;
}

  
static Boolean tc2_toNode(TreeCursorPtr cursor, Int4 tax_id)
{
    if(!tax_ptree_toTaxId(cursor, tax_id, FALSE)) {
	/* this node is not in our tree */
	return (tax_ptree_addNode(tax_tree, tax_id)) ? 
	    tax_ptree_toTaxId(cursor, tax_id, FALSE) : FALSE;
    }
    return TRUE;
}

static void lockBuff(int mode)
{
    mode= mode;
}

static void unlockBuff(void)
{
}

static void initBuff(void)
{
    int i;

    my_timer= 0;
  
    for(i= 0; i < BUFF_SIZE; i++) {
	or_buff[i].tax_id= 0;
	or_buff[i].p_org_ref= NULL;
    }
}

static Int4 getLiveId(Int4 id)
{
    TreeCursorPtr cursor= tree_openCursor(tax_tree, NULL, NULL);
    Uint2 s;
    TXC_TreeNodePtr tnp;
    
    if((cursor == NULL) || (!tc2_toNode(cursor, id))) return 0;

    tnp= tree_getNodeData(cursor, &s);
    tree_closeCursor(cursor);
    return (tnp != NULL)? tnp->tax_id : 0;
}


/**************************************************************************
 *
 *	InitTaxDB
 *
 **************************************************************************/

int InitTaxDB(void)
{
    CharPtr tmp;

    if((tmp=getenv("TAXDBPATH")) != NULL) DB_PATH= tmp;

    if(!txc_connect2Server(DB_PATH, "soussov", "vladimir", "tax2cl")) {

	return 0;
    }
       
    if((!txc_loadNameClasses()) || (!txc_loadRanks()) || 
       (!txc_loadDivisions()) || (!txc_loadGCs())) {

	return 0;
    }


    SpeciesRank=    tax_getRankId("species");
    SubspeciesRank= tax_getRankId("subspecies");
    GenusRank=      tax_getRankId("genus");
    FamilyRank=     tax_getRankId("family");
    OrderRank=      tax_getRankId("order");
    ClassRank=      tax_getRankId("class");


    VRL_div= tax_getDivisionId("VRL", NULL);
    PHG_div= tax_getDivisionId("PHG", NULL);

    SYNONYM=     tax_getClass_cde("synonym");
    COMMON_NAME= tax_getClass_cde("common name");
    PREF_COMMON= tax_getClass_cde("preferred common name");

    initBuff();
    tax_tree= tax_ptree_new();
  
    return (tax_tree == NULL)? 0 : 1;
}

/**************************************************************************
 *
 *	CloseTaxDB
 *
 **************************************************************************/

int CloseTaxDB(void)
{
    int i;

    if(tax_tree != NULL) {
	tree_delete(tax_tree);

	txc_close();

	for(i= 0; i < BUFF_SIZE; i++) {
	    if(or_buff[i].p_org_ref != NULL) {
		OrgRefFree(or_buff[i].p_org_ref);
	    }
	}
    }
    return 1;
}

Int4 tax1_getParent(Int4 id_tax)
{
    TreeCursorPtr cursor;
    Int4 ret_id= -1;

    if(tax_tree == NULL) return -1;

    if(id_tax == 1) return 0;

    cursor= tree_openCursor(tax_tree, NULL, NULL);
    if(cursor == NULL) return -1;

    if(tc2_toNode(cursor, id_tax)) {
	TXC_TreeNodePtr tnp;
	Uint2 s;

	tree_parent(cursor);
	tnp= tree_getNodeData(cursor, &s);
	if(tnp != NULL) ret_id= tnp->tax_id;
    }
    
    tree_closeCursor(cursor);
    return ret_id;
}

int tax1_getChildren(Int4 id_tax, Int4** ids_out)
{
    int n= 0;
    Int4* ids;
    TreeCursorPtr cursor= tree_openCursor(tax_tree, NULL, NULL);

    *ids_out= NULL;

    if(cursor == NULL) return 0;

    if(tc2_toNode(cursor, id_tax)) {
	if(tax_ptree_addChildren(cursor) && tree_child(cursor)) {
	    TXC_TreeNodePtr tnp;
	    Uint2 s;

	    do {
		n++;
	    }
	    while(tree_sibling(cursor));

	    *ids_out= ids= MemNew(n*sizeof(Int4));

	    tree_parent(cursor);
	    tree_child(cursor);
	    n= 0;
	    do {
		tnp= tree_getNodeData(cursor, &s);
		if(tnp != NULL) ids[n++]= tnp->tax_id;
	    }
	    while(tree_sibling(cursor));
	}
    }

    tree_closeCursor(cursor);
    return n;
}
    

/* find last common ancestor for two nodes */
Int4 tax1_join(Int4 taxid1, Int4 taxid2)
{
    TreeNodeId nid;
    Int4 aid= 0;
    TreeCursorPtr cursor1= tree_openCursor(tax_tree, NULL, NULL);
    TreeCursorPtr cursor2= tree_openCursor(tax_tree, NULL, NULL);

    if((cursor1 == NULL) || (cursor2 == NULL) || 
       (!tc2_toNode(cursor1, taxid1)) || (!tc2_toNode(cursor2, taxid2))) {
	if(cursor1 != NULL) tree_closeCursor(cursor1);
	if(cursor2 != NULL) tree_closeCursor(cursor2);
	return -1;
    }

    nid= tree_getAncestor(cursor1, cursor2);

    if(tree_toNode(cursor1, nid)) {
	TXC_TreeNodePtr tnp;
	Uint2 s;

	tnp= tree_getNodeData(cursor1, &s);
	if(tnp != NULL) aid= tnp->tax_id;
    }

    tree_closeCursor(cursor1);
    tree_closeCursor(cursor2);
    return aid;
}
	

/***************************************************
 * Get tax_id by organism name
 * returns:
 *       tax_id if one node found
 *       0      no organism found
 *       -tax_id if more than one node found
 */
Int4 tax1_getTaxIdByName(CharPtr orgname)
{
    return tax_getIdByName(orgname, NULL, 0);
}

/***************************************************
 * Get all tax_id by organism name
 * returns:
 *       Number of tax ids found
 */
Int4 tax1_getAllTaxIdByName(CharPtr orgname, Int4 **Ids_out)
{
    int i;
    TaxNamePtr nameList;
    Int4 *Ids;
    int n= tax_findByName(orgname, TAX_NAME_SEARCH, &nameList);

    if(n < 1) return 0;

    *Ids_out= Ids= MemNew(n*sizeof(Int4));
    if(Ids != NULL) {
	for(i= 0; i < n; i++) {
	    Ids[i]= nameList[i].tax_id;
	    if(nameList[i].name_txt != NULL) MemFree(nameList[i].name_txt);
	    if(nameList[i].unique_name != NULL) MemFree(nameList[i].unique_name);
	}
    }
    else {
	for(i= 0; i < n; i++) {
	    if(nameList[i].name_txt != NULL) MemFree(nameList[i].name_txt);
	    if(nameList[i].unique_name != NULL) MemFree(nameList[i].unique_name);
	}
	n= 0;
    }
    MemFree(nameList);
    return n;
}

static Boolean goodOrgMode(Uint1 t)
{
    return (t != 254) && (t != 20);
}
	
static Int4 getIdByOrgRef(CharPtr sname, OrgNamePtr orNm)
{
    if(orNm != NULL) {
	OrgModPtr o_mod= orNm->mod;
	Boolean om_flag= FALSE;
	Int4 tax_id, id;
	CharPtr altname= NULL;
	int nof_mods= 0;
	_subspecPtr ss;
	_subspec src;

	/* first try to search using search name */
	for(;o_mod != NULL; o_mod= o_mod->next) {
	    if(o_mod->subtype == 254) {
		/* search name */
		altname= o_mod->subname;
	    }
	    else if(o_mod->subtype != 20) nof_mods++;
	}

	if(nof_mods == 0) {
	    /* we have no modifiers */
	    if(altname != NULL) {
		if((tax_id= tax_getIdByName(altname, NULL, 0)) > 0) return tax_id;
		return tax_getIdByName(sname, altname, 254);
	    }
	    return tax_getIdByName(sname, NULL, 0);
	}
		
	if(nof_mods == 1) {
	    /* we have one valuable modifier */
	    for(o_mod= orNm->mod; o_mod != NULL; o_mod= o_mod->next) {
		if(goodOrgMode(o_mod->subtype)) {
		    if(altname != NULL) {
			if((tax_id= tax_getIdByName(altname, o_mod->subname, o_mod->subtype)) > 0) 
			    return tax_id; /* find by old name and modifier */
			if((tax_id= tax_getIdByName(sname, o_mod->subname, o_mod->subtype)) > 0) 
			    return tax_id; /* find by new name and modifier */
			return tax_getIdByName(sname, altname, 254); /* find by new name and old name */
		    }
		    return tax_getIdByName(sname, o_mod->subname, o_mod->subtype);
		}
	    }
	    return 0;
	}

	/* we have more than one modifier */
	/* first try to find organism using just names */
	if(altname != NULL) {
	    if((tax_id= tax_getIdByName(altname, NULL, 0)) == 0) {
		tax_id= tax_getIdByName(sname, altname, 254);
	    }
	}
	else {
	    tax_id= tax_getIdByName(sname, NULL, 0);
	}

	if(tax_id == 0) return 0; /* we have no such names */
	if(tax_id > 0) {
	    /* we have found just one node, check it against modifiers */
	    id= 0;
	    for(o_mod= orNm->mod; o_mod != NULL; o_mod= o_mod->next) {
		if(o_mod->subtype != 20) {
		    src.stype= o_mod->subtype;
		    src.sname= o_mod->subname;
		    src.rname= NULL;
		    if((ss= tax_SSget(tax_id, &src)) != NULL) {
			if(ss->rname != NULL) MemFree(ss->rname);
			if((ss->r_id != 0) && (ss->r_id != tax_id)) {
			    if(id == 0) id= ss->r_id;
			    else if(id != ss->r_id) {
				id= -id; /* conflict in mapping */
				break;
			    }
			}
		    }
		}
	    }
	    if(id == 0) return tax_id;
	    if(id < 0) return -tax_id; /* we have a mapping conflict */
	    
	    /* we have a mapping without conflict, we try to make the best assumption */
	    return id;
	}

	if(tax_id < 0) {
	    /* more than one tax_id was found */
	    Int4Ptr ids;
	    Int4 n;

	    if(altname != NULL) {
		n= tax1_getAllTaxIdByName(altname, &ids);
		if(n < 1) n= tax1_getAllTaxIdByName(sname, &ids);
	    }
	    else n= tax1_getAllTaxIdByName(sname, &ids);

	    id= 0;
	    while(n-- > 0) {
		for(o_mod= orNm->mod; o_mod != NULL; o_mod= o_mod->next) {
		    if(goodOrgMode(o_mod->subtype)) {
			src.stype= o_mod->subtype;
			src.sname= o_mod->subname;
			src.rname= NULL;
			if((ss= tax_SSget(ids[n], &src)) != NULL) {
			    if(ss->rname != NULL) MemFree(ss->rname);
			    if(ss->r_id != 0) {
				if(id == 0) id= ss->r_id;
				else if(id != ss->r_id) id= -id;
			    }
			}
		    }
		}
	    }
	    if(ids != NULL) MemFree(ids);
	    if(id > 0) return id;
	}
	return tax_id;
    }
    else {
	/* we have no modifiers */
	return tax_getIdByName(sname, NULL, 0);
    }

}

Int4 tax1_getTaxIdByOrgRef(OrgRefPtr orgRef)
{
#ifdef TAXSERVICE
    return txc_getTaxIdByOrgRef(orgRef);
#else
    Int4 tax_id, id;

    if(orgRef == NULL) return 0;

    tax_id= 0;

    if((orgRef->taxname != NULL) &&
       ((tax_id= getIdByOrgRef(orgRef->taxname, orgRef->orgname)) > 0)) return tax_id;
	
    if((orgRef->common != NULL) &&
       ((tax_id= getIdByOrgRef(orgRef->common, orgRef->orgname)) > 0)) return tax_id;

    if(orgRef->syn != NULL) {
	ValNodePtr synonym;
	
	id= 0;

	for(synonym= orgRef->syn; (synonym != NULL) && (id < 1); synonym= synonym->next) {
	    id= getIdByOrgRef(synonym->data.ptrvalue, orgRef->orgname);
	}
    }

    return (id > 0)? id : tax_id;
#endif
}

/*******************************************************************
 * Get tax_id by organism name (it could be "unique" variant of name)
 * returns:
 *       tax_id if one node found
 *       0      no organism found
 *       -tax_id if more than one node found
 */
Int4 tax1_findTaxIdByName(CharPtr orgname)
{
    Int4 id= tax_getIdByName(orgname, NULL, 0);

    if(id < 1) {
	Int4 idu= tax_uniqueName(orgname, 0);

	if(idu > 0) id= idu;
    }
    return id;
}


/*************************************************************************
 * Get all tax_id by organism name (it could be "unique" variant of name)
 * returns:
 *       Number of tax ids found
 */
Int4 tax1_findAllTaxIdByName(CharPtr orgname, Int4 **Ids_out)
{
    Int4 id= tax1_findTaxIdByName(orgname);
    
    if(id > 0) {
	*Ids_out= MemNew(sizeof(Int4));
	if(*Ids_out != NULL) {
	    **Ids_out= id;
	    return 1;
	}
	else return 0;
    }
    
    if(id < 0) {
	return tax1_getAllTaxIdByName(orgname, Ids_out);
    }
   
    return 0;
}

Int2 tax1_getAllNames(Int4 tax_id, CharPtr **out_names, Boolean unique)
{
    TaxNamePtr nameList;
    Int2 n= tax_getOrgNames(tax_id, &nameList);
    Int2 i;
    CharPtr* names;

    if(n < 1) return 0;

    *out_names= names= MemNew(n*sizeof(CharPtr));
    if(names != NULL) {
        for(i= 0; i < n; i++) {
	    if(unique && (nameList[i].unique_name != NULL) && (nameList[i].unique_name[0] != '\0')) {
		names[i]= nameList[i].unique_name;
		nameList[i].unique_name= NULL;
	    }
	    else {
		names[i]= nameList[i].name_txt;
		nameList[i].name_txt= NULL;
	    }
	}
    }
	
    for(i= 0; i < n; i++) {
	if(nameList[i].name_txt != NULL) MemFree(nameList[i].name_txt);
	if(nameList[i].unique_name != NULL) MemFree(nameList[i].unique_name);
    }
    
    MemFree(nameList);
    return n;
}


CharPtr tax1_getGCName(Int2 gc_id)
{
    return tax_getGCName(gc_id);
}

OrgRefPtr tax1_getOrgRef(Int4 tax_id, int* is_species, CharPtr div, CharPtr embl_cde)
{
    OrgRefPtr orp;

    tax_id= getLiveId(tax_id);
    if(tax_id == 0) return NULL;

    if((orp= getFromBuff(tax_id, is_species, div, embl_cde)) != NULL) {
	/* OrgRef is already in buffer */
	return orp;
    }

    lockBuff(TAX_WRITE);
    loadInBuff(tax_id);
    unlockBuff();

    return getFromBuff(tax_id, is_species, div, embl_cde);
}

Taxon1DataPtr tax1_getbyid(Int4 tax_id)
{
    Taxon1DataPtr res;
    OrgRefPtr db_orgRef;
    int is_species;

    if(tax_id <= 0) return NULL;
    res= Taxon1DataNew();
    res->div= MemNew(16);
    res->embl_code= MemNew(4);
    db_orgRef= tax1_getOrgRef(tax_id, &is_species, res->div, NULL/*res->embl_code*/);
    res->embl_code[0]= '\0';
    if(db_orgRef == NULL) {
	Taxon1DataFree(res);
	return NULL;
    }

    /* make new orgref */
    res->org= OrgRefNew();
    res->org->db= NULL;
    res->org->orgname= NULL;
    res->is_species_level= is_species;

    /* fill-up orgref based on db_orgRef */
    bldOrgRefOut(res->org, db_orgRef, getLiveId(tax_id));
    return res;
}

/*************************************************************************/
/* return pointer to first non-blank character in str1 after prefix str2 */
/* if str2 is not prefix for str1 then return str1                       */
/*************************************************************************/
static CharPtr strTail(CharPtr str1, CharPtr str2)
{
    CharPtr c;

    if(StringStr(str1, str2) != str1) return str1;

    c= str1 + StringLen(str2);

    while((*c != '\0') && IS_WHITESP(*c)) c++;

    return c;
}



static OrgModPtr bldOrgMod(TreeCursorPtr cursor)
{
    TXC_TreeNodePtr parent= NULL;
    Uint2 s;
    TXC_TreeNodePtr me= tree_getNodeData(cursor, &s);
    TXC_TreeNodePtr tnp;
    TreeNodeId nid= tree_getId(cursor);
    Int2 rank, prank;
    OrgModPtr orgMdf= OrgModNew();

    while(tree_parent(cursor)) {
	if((tnp= tree_getNodeData(cursor, &s)) == NULL) continue;
	prank= tnp->flags & 0xFF;
	--prank;
	if((prank == SubspeciesRank) || 
	   (prank == SpeciesRank) ||
	   (prank == GenusRank)) {
	    parent= tnp;
	    break;
	}
    }
    tree_toNode(cursor, nid);

    if(parent != NULL) {
	orgMdf->subname= StringSave(strTail(me->node_label, parent->node_label));
    }
    else {
	orgMdf->subname= StringSave(me->node_label);
    }

    rank= me->flags & 0xFF;

    if(--rank == SubspeciesRank) {
	orgMdf->subtype= 22; /* subspecies */
    }
    else if(rank == tax_getRankId("varietas")) {
	orgMdf->subtype= 6; /* variety */
    }
    else if(rank == tax_getRankId("forma")) {
	orgMdf->subtype= 2; /* strain */
    }
    else if((parent != NULL) && (prank == SubspeciesRank)) {
	orgMdf->subtype= 2; /* strain */
    }
    else {
	orgMdf->subtype= 255; /* other */
    }

    orgMdf->attrib= NULL;

    return orgMdf;
}


/***************************************************************/  
/* if name is binomial name build the correspondent structures */
/* otherwise return 0                                          */
/***************************************************************/
static int binomialName(TreeCursorPtr cursor, OrgNamePtr onp)
{
    TXC_TreeNodePtr spec= NULL;
    TXC_TreeNodePtr subspec= NULL;
    TXC_TreeNodePtr genus= NULL;
    TXC_TreeNodePtr tnp;
    Uint2 s;
    TreeNodeId nid= tree_getId(cursor);
    Int2 rank;
    BinomialOrgNamePtr bName;

    do {
	tnp= tree_getNodeData(cursor, &s);
	if(tnp == NULL) continue;
	rank= tnp->flags & 0xFF;
	if(--rank == SubspeciesRank) subspec= tnp;
	else if(rank == SpeciesRank) spec= tnp;
	else if(rank == GenusRank) {
	    genus= tnp;
	    break;
	}
    }
    while(tree_parent(cursor));

    tree_toNode(cursor, nid);

    if(genus == NULL) {
	/* try to find subgenus */
	do {
	    tnp= tree_getNodeData(cursor, &s);
	    if(tnp == NULL) continue;
	    rank= tnp->flags & 0xFF;
	    if(--rank == (GenusRank + 1)) {
		genus= tnp;
		break;
	    }
	}
	while(tree_parent(cursor));
	tree_toNode(cursor, nid);
    }

    if(genus == NULL) return 0; /* no genus - no binomial */

    onp->choice= 1; /*binomial*/

    onp->data= bName= BinomialOrgNameNew();
  
    bName->genus= StringSave(genus->node_label);
    

    if(spec != NULL) {
	/* we have a species in lineage */
	bName->species= StringSave(strTail(spec->node_label, genus->node_label));

	if(subspec != NULL) {
	    /* we also have a subspecies in lineage */
	    bName->subspecies= StringSave(strTail(subspec->node_label, spec->node_label));
	}
	else {
	    bName->subspecies= NULL;
	}
	tnp= tree_getNodeData(cursor, &s);

	onp->mod= (tnp == spec)? NULL : bldOrgMod(cursor);    
	return 1;
    }
  
    /* no species in lineage */

    if(subspec != NULL) {
	/* we have no species but we have subspecies */
	bName->species= NULL;
	bName->subspecies= StringSave(strTail(subspec->node_label, genus->node_label));
	onp->mod= bldOrgMod(cursor);
	return 1;
    }
  
    /* we have no species, no subspecies but we are under species level (varietas or forma) */

    bName->species= NULL;
    bName->subspecies= NULL;
    onp->mod= bldOrgMod(cursor);
    return 1;
}

      
static void partialName(TreeCursorPtr cursor, OrgNamePtr onp)
{
    TaxElementPtr taxElem;
    Uint2 s;
    TXC_TreeNodePtr tnp= tree_getNodeData(cursor, &s);
    Int2 rank_id= tnp->flags & 0xFF;

    onp->choice= 5; /* partial */
    onp->data= taxElem= TaxElementNew();
  
    if(--rank_id == FamilyRank) {
	taxElem->fixed_level= 1; /* family */
	taxElem->level= NULL;
    }
    else if(rank_id == OrderRank) {
	taxElem->fixed_level= 2;
	taxElem->level= NULL;
    }
    else if(rank_id == ClassRank) {
	taxElem->fixed_level= 3;
	taxElem->level= NULL;
    }
    else {
	taxElem->fixed_level= 0;
	taxElem->level= StringSave(tax_getRank(rank_id));
    }

    taxElem->name= StringSave(tnp->node_label);
    taxElem->next= NULL;
}

  
/*****************************************************************
 * build synonyms valnodes
 * this routine include in valnodes synonyms and common synonyms
 */
static ValNodePtr bldSynValNodes(TaxNamePtr syn, Int2 n)
{
    ValNodePtr list= NULL;
    ValNodePtr header= NULL;
    Int2 i;

    for(i= 1; i < n; i++) {
	if((syn[i].class_cde == SYNONYM) || (syn[i].class_cde == COMMON_NAME)) {
	    list= ValNodeNew(list);
	    list->choice= (syn[i].class_cde == SYNONYM)? 1 : 0;
	    list->data.ptrvalue= syn[i].name_txt;
	    syn[i].name_txt= NULL;
	    if(header == NULL) header= list;
	}
    }
    return header;
}
	   

/**************************************************************************
 *
 *	FsTaxGetLineage
 *
 **************************************************************************/

static CharPtr bldLineage(TreeCursorPtr cursor)
{
    TXC_TreeNodePtr tnp;
    Uint2 s;
    TreeNodeId nid= tree_getId(cursor);
    CharPtr lineage= NULL;
    CharPtr tmp, t;
    Int2 rank;
    

    while(tree_parent(cursor)) {
	if((tnp= tree_getNodeData(cursor, &s)) != NULL) {
	    rank= tnp->flags & 0xFF;
	    if(rank > SpeciesRank) {
		if(lineage != NULL) {
		  lineage = MemFree(lineage);
		}
		continue;
	    }
	    if(tnp->tax_id < 2) break;

	    if((tnp->flags & TXC_GBHIDE) == 0) {
		s= StringLen(tnp->node_label);
		if(lineage != NULL) {
		    s+= StringLen(lineage) + 2;
		}

		tmp= MemNew(s+2);
		if(tmp == NULL) continue;
		t= StringMove(tmp, tnp->node_label);
		if(lineage != NULL) {
		    t= StringMove(t, "; ");
		    t= StringMove(t, lineage);
		    MemFree(lineage);
		}
		lineage= tmp;
	    }
	}
    }

    tree_toNode(cursor, nid);

    return lineage;
}

static ValNodePtr bldDBId(Int4 id)
{
    ValNodePtr dbnode;
    DbtagPtr dbtag;
    ObjectIdPtr object_id;

    /* populate tax_id */
    dbnode= ValNodeNew(NULL);
    dbnode->data.ptrvalue= dbtag= DbtagNew();
    dbtag->db = StringSave("taxon");
    dbtag->tag= object_id= ObjectIdNew();
    object_id->str= NULL;
    object_id->id = id;
    return dbnode;
}

static OrgNamePtr bldOrgName(TreeCursorPtr cursor, int* is_species_out, CharPtr div, CharPtr embl)
{
    OrgNamePtr onp;
    Uint2 s;
    TXC_TreeNodePtr tnp= tree_getNodeData(cursor, &s);
    CharPtr div_abbr;
    Int2 rank_id;
    int is_species;
    Int2 div_id;
    TreeNodeId nid= tree_getId(cursor);
    Int2 rank;

    /*Int4 p_id, s_id;*/

    if(tnp == NULL) return NULL;

    onp= OrgNameNew();

    rank_id= tnp->flags & 0xFF;
    rank_id--;
      
    div_id= (tnp->flags >> 8) & 0x3F;
    onp->gcode= (tnp->flags >> (8+6)) & 0x3F;
    onp->mgcode= (tnp->flags >> (8+6+6)) & 0x3F;
    onp->lineage= bldLineage(cursor);
    if(embl != NULL) *embl= '\0';
	   
    is_species= (rank_id >= SpeciesRank)? 1 : 0;
    /* correct level by lineage if node has no rank */
    if(rank_id < 0) {

	while(tree_parent(cursor)) {
	    tnp= tree_getNodeData(cursor, &s);
	    if(tnp != NULL) {
		rank= tnp->flags & 0xFF;
		if(rank != 0) {
		    is_species= (rank >= SpeciesRank) ? 1 : 0;
		    break;
		}
	    }
	}
	tree_toNode(cursor, nid);
	tnp= tree_getNodeData(cursor, &s);
    }

    if(tax_getDivision(div_id, &div_abbr, NULL)) {
	onp->div= StringSave(div_abbr);
	StringCpy(div, div_abbr);
    }
    *is_species_out= is_species;

    if(is_species) {
	/* we are on species level or below */
	     
	/* check for viruses */
	if((div_id == VRL_div) || (div_id == PHG_div)) {
	    /* this is a virus */
	    onp->choice= 2; /* virus */
	    if(rank_id == SpeciesRank) {
		/* we are on species level */
		onp->data= StringSave(tnp->node_label);
		onp->mod= NULL;
	    }
	    else {
		/* we are below species */
		/* first try to find species or min rank which below species but above us */
		TreeNodeId s_id;
		
		s_id.idi= 0;
		while(tree_parent(cursor)) {
		    tnp= tree_getNodeData(cursor, &s);
		    if(tnp != NULL) {
			rank= tnp->flags & 0xFF;
			if(--rank >= SpeciesRank) {
			    s_id= tree_getId(cursor);
			    if(rank == SpeciesRank) break;
			}
			else if(rank >= 0) break;
		    }
		}
		if(s_id.idi != 0) {
		    /* we have species or something above us */
		    tree_toNode(cursor, s_id);
		}
		else {
		    /* no species above */
		    tree_toNode(cursor, nid);
		}
		
		tnp= tree_getNodeData(cursor, &s);
		onp->data= StringSave(tnp->node_label);
		if(s_id.idi != 0) {
		    tree_toNode(cursor, nid);
		}

		onp->mod= bldOrgMod(cursor);
	    }
	}		
	else if(!binomialName(cursor, onp)) {
	    /* name is not binomial: set partial */
	    partialName(cursor, onp);
	}
    }
    else {
	/* above species */
	partialName(cursor, onp);
    }

  return onp;
}


static Boolean bldOrgRef(Int4 id, OrgRefPtr orp, int* is_species, CharPtr div, CharPtr embl)
{
    TreeCursorPtr cursor= tree_openCursor(tax_tree, NULL, NULL);
    TaxNamePtr nameList;
    Int2 n, i; 
    
    if((cursor == NULL) || (!tc2_toNode(cursor, id))) return FALSE;

    *is_species= 0;
    *div= *embl= '\0';

    n= tax_getOrgNames(id, &nameList);

    if(n < 1) {
	tree_closeCursor(cursor);
	return FALSE;
    }
    
    orp->taxname= nameList[0].name_txt;
    nameList[0].name_txt= NULL;
    

    /* fill-up preferred common name */
    orp->common= NULL;
    for(i= 1; i < n; i++) {
	if(nameList[i].class_cde == PREF_COMMON) {
	    orp->common= nameList[i].name_txt;
	    nameList[i].name_txt= NULL;
	    break;
	}
    }

    /* fill-up synonyms */
    orp->syn= bldSynValNodes(nameList, n);
    for(i= 0; i < n; i++) {
	if(nameList[i].name_txt != NULL) MemFree(nameList[i].name_txt);
	if(nameList[i].unique_name != NULL) MemFree(nameList[i].unique_name);
    }
    
    MemFree(nameList);

    orp->mod= NULL;
    orp->db= bldDBId(id);
    orp->orgname= bldOrgName(cursor, is_species, div, embl);
    

    tree_closeCursor(cursor);
    return TRUE;	
}
  

static void loadInBuff(Int4 id)
{
    int i, k= -1;
    Int4 t= my_timer + 1;
    /*Int4 bt;*/
  
    for(i= 0; i < BUFF_SIZE; i++) {
	if(or_buff[i].tax_id == 0) {
	    k= i;
	    break;
	}
	if(or_buff[i].timer < t) {
	    t= or_buff[i].timer;
	    k= i;
	}
    }

    if(k >= 0) {
	if(or_buff[k].p_org_ref != NULL) OrgRefFree(or_buff[k].p_org_ref);
	or_buff[k].tax_id= id;
	or_buff[k].p_org_ref= OrgRefNew();
	or_buff[k].timer= ++my_timer;
	if(!bldOrgRef(id, or_buff[k].p_org_ref, &or_buff[k].is_species, or_buff[k].div, or_buff[k].embl)) {
	    OrgRefFree(or_buff[k].p_org_ref);
	    or_buff[k].tax_id= 0;
	}
    }
}

static OrgRefPtr getFromBuff(Int4 id, int* is_sp, CharPtr div, CharPtr embl)
{
    int i;
    OrgRefPtr orp= NULL;

    lockBuff(TAX_READ);

    for(i= 0; i < BUFF_SIZE; i++) {
	if(or_buff[i].tax_id == id) {
	    or_buff[i].timer= ++my_timer;
	    orp= or_buff[i].p_org_ref;
	    if(is_sp != NULL) *is_sp= or_buff[i].is_species;
	    if(div != NULL) StringCpy(div, or_buff[i].div);
	    if(embl != NULL) StringCpy(embl, ""/*or_buff[i].embl*/);
	    break;
	}
    }
    unlockBuff();
    return orp;
}

static OrgModPtr fixModifier(Int4 tax_id, OrgModPtr omp)
{
    _subspec src_ss;
    _subspecPtr ss= NULL;

    memset(&src_ss, 0, sizeof(_subspec));

    if(omp->subtype < 2) {
	omp->next= NULL;
	OrgModFree(omp);
	return NULL;
    }
      

    if((omp->subname != NULL) && (omp->subtype != 0)) {
	src_ss.stype= omp->subtype;
	src_ss.sname= omp->subname;
	src_ss.rname= NULL;

	ss= tax_SSget(tax_id, &src_ss);
    }

    if((ss != NULL) && (ss->r_id == tax_id) && (ss->stype == 0)) {
	/* remove it */
	if(ss->rname != NULL) MemFree(ss->rname);
	omp->next= NULL;
	OrgModFree(omp);
	return NULL;
    }

    if((ss != NULL) && (ss->r_id == tax_id) && (ss->stype != 0)) {		
	MemFree(omp->subname);
	omp->subname= src_ss.rname;
	omp->subtype= src_ss.rtype;
	return omp;
    }

    if(src_ss.rname != NULL) MemFree(src_ss.rname);
    return omp;
}

static void CleanOrgMod(Int4 tax_id, OrgNamePtr onp)
{
    OrgModPtr omp, omp_p, omp_n;

    for(omp_p= NULL, omp= onp->mod; omp != NULL; omp= omp_n) {
	omp_n= omp->next;
	if((omp= fixModifier(tax_id, omp)) == NULL) {
	    /* exclude this modifier */
	    if(omp_p == NULL) {
		onp->mod= omp_n;
	    }
	    else omp_p->next= omp_n;
	}
	else omp_p= omp;
    }
}
	    
static void cleanOrgName(Int4 tax_id, OrgNamePtr onp)
{
    if(onp->lineage != NULL) MemFree(onp->lineage);
    if(onp->div != NULL) MemFree(onp->div);
#if 1
    /* #if 0 means that we will trust to initial modifier */
    if(onp->mod != NULL) CleanOrgMod(tax_id, onp);
#endif
    if(onp->next != NULL) OrgNameSetFree(onp->next);
    if(onp->data != NULL) {
	switch(onp->choice) {
	case 1 : /* binomial name */
	    BinomialOrgNameFree(onp->data);
	    break;
	case 2 : /* virus name */
	    MemFree(onp->data);
	    break;
	case 5 : /* partial name */
	    TaxElementSetFree(onp->data);
	    break;
	}
    }
}

  
static BinomialOrgNamePtr copyBinomial(BinomialOrgNamePtr src)
{
    BinomialOrgNamePtr dst;

    if(src == NULL) return NULL;

    dst= BinomialOrgNameNew();
    dst->genus= (src->genus != NULL)? StringSave(src->genus) : NULL;
    dst->species= (src->species != NULL)? StringSave(src->species) : NULL;
    dst->subspecies= (src->subspecies != NULL)? StringSave(src->subspecies) : NULL;

    return dst;
}

static TaxElementPtr copyPartial(TaxElementPtr src)
{
    TaxElementPtr dst;

    if(src == NULL) return NULL;

    dst= TaxElementNew();
    dst->fixed_level= src->fixed_level;
    dst->level= (src->level != NULL)? StringSave(src->level) : NULL;
    dst->name= (src->name != NULL)? StringSave(src->name) : NULL;
    dst->next= (src->next != NULL)? copyPartial(src->next) : NULL;
    return dst;
}

static OrgModPtr copyOrgMod(OrgModPtr src)
{
#if 0
    OrgModPtr dst;

    if(src == NULL) return NULL;

    dst= OrgModNew();
    dst->subtype= src->subtype;
    dst->subname= (src->subname != NULL)? StringSave(src->subname) : NULL;
    dst->attrib= (src->attrib != NULL)? StringSave(src->attrib) : NULL;
    dst->next= (src->next != NULL)? copyOrgMod(src->next) : NULL;

    return dst;
#endif
    return NULL;
}
  
static ValNodePtr removeDbtag(ValNodePtr vnp)
{
    ValNodePtr vnn, vnf, vnl= NULL;
    DbtagPtr dbtag;

    for(vnf= vnp; vnp != NULL; vnp= vnn) {
	dbtag= vnp->data.ptrvalue;
	vnn= vnp->next;
	if(dbtag == NULL) return NULL;
	if(StringCmp(dbtag->db, "taxon") == 0) {
	    /* taxon tag, remove it */
	    if(vnl == NULL) {
		vnf= vnn;
	    }
	    else {
		vnl->next= vnn;
	    }
	    DbtagFree(dbtag);
	    MemFree(vnp);
	}
	else {
	    vnl= vnp;
	}
    }
    return vnf;
}
	

static void bldOrgRefOut(OrgRefPtr dst, OrgRefPtr src, Int4 tax_id)
{
    ValNodePtr vnp, vnl;
    DbtagPtr dbtag;
    ObjectIdPtr object_id;
    OrgNamePtr onp;

    dst->taxname= StringSave(src->taxname);
    dst->common= (src->common != NULL)? StringSave(src->common) : NULL;

    /* populate tax_id */
    vnp= ValNodeNew(NULL);
    if (dst->db != NULL) {
	dst->db= removeDbtag(dst->db);
    }
    vnp->next= dst->db;
    dst->db= vnp;
    vnp->data.ptrvalue= dbtag= DbtagNew();
    dbtag->db = StringSave("taxon");
    dbtag->tag= object_id= ObjectIdNew();
    object_id->str= NULL;
    object_id->id = getLiveId(tax_id);
    
    /* copy the synonym list */
    dst->syn= NULL; vnl= NULL;
    if(we_want_synonyms) {
	for(vnp= src->syn; vnp != NULL; vnp= vnp->next) {
	    vnl= ValNodeNew(vnl);
	    vnl->choice= vnp->choice;
	    vnl->data.ptrvalue= StringSave(vnp->data.ptrvalue);
	    if(dst->syn == NULL) dst->syn= vnl;
	}
    }
  
    /* copy orgname */
    if(dst->orgname == NULL) dst->orgname= onp= OrgNameNew();
    else onp= dst->orgname;

    onp->choice= src->orgname->choice;

    switch(src->orgname->choice) {
    case 1 : /*binomial*/
	onp->data= copyBinomial(src->orgname->data);
	break;
    case 2 : /* virus */
	onp->data= (src->orgname->data != NULL)? StringSave(src->orgname->data) : NULL;
	break;
    case 5 : /* partial */
	onp->data= copyPartial(src->orgname->data);
	break;
    default : /* can't handle */
	onp->data= NULL;
    }
  
    if(onp->mod == NULL) onp->mod= copyOrgMod(src->orgname->mod);
    onp->lineage= (src->orgname->lineage != NULL)? StringSave(src->orgname->lineage) : NULL;
    onp->gcode= src->orgname->gcode;
    onp->mgcode= src->orgname->mgcode;
    onp->div= StringSave(src->orgname->div);
}

static void populateReplaced(OrgRefPtr orp, CharPtr oldName)
{
    OrgNamePtr onp;
    OrgModPtr omp;

    if((orp->taxname != NULL) && (StringICmp(orp->taxname, oldName) == 0)) {
	MemFree(oldName);
	return;
    }

    if((orp->common != NULL) && (StringICmp(orp->common, oldName) == 0)) {
	MemFree(oldName);
	return;
    }

    /* organism name was changed */
    onp= orp->orgname;
    if(onp != NULL) {
	omp= OrgModNew();
	omp->next= onp->mod;
	omp->subtype= 254;
	omp->subname= oldName;
	onp->mod= omp;
    }
    else {
	MemFree(oldName);
    }
}

Taxon1DataPtr tax1_lookup(OrgRefPtr inp_orgRef, int merge)
{
    Taxon1DataPtr res;
    Int4 tax_id;
    OrgRefPtr db_orgRef;
    int is_species;
    Boolean need_search_name= TRUE;
    CharPtr hit_name;

    tax_id= tax1_getTaxIdByOrgRef(inp_orgRef);
    if(tax_id <= 0) return NULL;
    res= Taxon1DataNew();
    res->div= MemNew(16);
    res->embl_code= MemNew(4);
    db_orgRef= tax1_getOrgRef(tax_id, &is_species, res->div, NULL /*res->embl_code*/);
    res->embl_code[0]= '\0';
    if(db_orgRef == NULL) {
	Taxon1DataFree(res);
	return NULL;
    }

    res->is_species_level= is_species;

    /* populate search name if necessary */
    if(inp_orgRef->taxname != NULL) {
	if((db_orgRef->taxname != NULL) && (StringICmp(inp_orgRef->taxname, db_orgRef->taxname) == 0)) {
	    need_search_name= FALSE;
	}
	else if((db_orgRef->common != NULL) && (StringICmp(inp_orgRef->taxname, db_orgRef->common) == 0)) {
	    need_search_name= FALSE;
	}
    }

    if(need_search_name && (inp_orgRef->common != NULL)) {
	if((db_orgRef->taxname != NULL) && (StringICmp(inp_orgRef->common, db_orgRef->taxname) == 0)) {
	    need_search_name= FALSE;
	}
	else if((db_orgRef->common != NULL) && (StringICmp(inp_orgRef->common, db_orgRef->common) == 0)) {
	    need_search_name= FALSE;
	}
    }

    if(need_search_name && (inp_orgRef->orgname != NULL)) {
	/* check if search name already exists */
	OrgModPtr omp;

	for(omp= inp_orgRef->orgname->mod; omp != NULL; omp= omp->next) {
	    if(omp->subtype == 254) {
		need_search_name= FALSE;
		break;
	    }
	}
    }

    hit_name= NULL;
    if(need_search_name) {
	if((inp_orgRef->taxname != NULL) && (inp_orgRef->taxname[0] != '\0')) {
	    hit_name= StringSave(inp_orgRef->taxname);
	}
	else if((inp_orgRef->common != NULL) && (inp_orgRef->common[0] != '\0')) {
	    hit_name= StringSave(inp_orgRef->common);
	}
    }

    if(merge) {
	/* we have to merge old orgref with the new one */
	res->org= inp_orgRef;
	/* clean-up old information */
	if(inp_orgRef->taxname != NULL) MemFree(inp_orgRef->taxname);
	if(inp_orgRef->common != NULL) MemFree(inp_orgRef->common);
	if(inp_orgRef->syn != NULL) ValNodeFreeData(inp_orgRef->syn);
	if(inp_orgRef->orgname != NULL) cleanOrgName(tax_id, inp_orgRef->orgname);
    }
    else {
	/* make new orgref */
	res->org= OrgRefNew();
	res->org->db= NULL;
	res->org->orgname= NULL;
    }
    /* fill-up orgref based on db_orgRef */
    bldOrgRefOut(res->org, db_orgRef, tax_id);
    if(need_search_name && (hit_name != NULL)) populateReplaced(res->org, hit_name);
    return res;
}
  
  
Boolean tax1_init(void)
{
    return InitTaxDB();
}

void tax1_fini(void)
{
    CloseTaxDB();
}

TreePtr tax1e_getTaxTreePtr(void)
{
    return tax_tree;
}

Boolean tax1e_invokeNode(Int4 tax_id)
{
    return tax_ptree_addNode(tax_tree, tax_id);
}

Boolean tax1e_invokeChildren(Int4 tax_id)
{
    Boolean res= FALSE;

    TreeCursorPtr cursor= tree_openCursor(tax_tree, NULL, NULL);
    if(tax_id < 0) {
	res= TRUE;
	tax_id= -tax_id;
    }

    if((cursor == NULL) || (!tc2_toNode(cursor, tax_id))) return FALSE;

    res= res? tax_ptree_addSubtree(cursor) : tax_ptree_addChildren(cursor);

    tree_closeCursor(cursor);
    return res;
}

Boolean tax1e_toNode(TreeCursorPtr cursor, Int4 tax_id)
{
    return tc2_toNode(cursor, tax_id);
}

Int4 tax1e_getTaxId(TreeCursorPtr cursor)
{
    Uint2 s;
    TXC_TreeNodePtr tnp;
    
    tnp= tree_getNodeData(cursor, &s);

    return (tnp != NULL)? tnp->tax_id : 0;
}

CharPtr tax1e_getTaxName(TreeCursorPtr cursor)
{
    Uint2 s;
    TXC_TreeNodePtr tnp;
    
    tnp= tree_getNodeData(cursor, &s);

    return (tnp != NULL)? StringSave(tnp->node_label) : NULL;
}
    
Int4 tax1_getTaxId4Str(CharPtr str, CharPtr* substring, Int4Ptr *Ids_out)
{
    CharPtr b, e;
    Int4 tax_id;
    /*Int4Ptr Ids;*/
    int k;
    char c;
    
    *substring= NULL;
    tax_id= tax1_getTaxIdByName(str);

    if(tax_id > 1) {
	*Ids_out= MemNew(sizeof(Int4));
	**Ids_out= tax_id;
	*substring= StringSave(str);
	return 1;
    }
    else if(tax_id < 0) {
	*substring= StringSave(str);
	return tax1_getAllTaxIdByName(str, Ids_out);
    }

    /* whole string matches nothing */
    /* try the whole string inside first set of parenthesis */
    for(b= str; *b != '\0'; b++) {
	if(*b == '(') {
	    k= 0;
	    for(e= b+1; *e != '\0'; e++) {
		if(*e == '(') {
		    k++;
		    continue;
		}
		if(*e == ')') {
		    if(k > 0) {
			k--; 
			continue;
		    }

		    *e= '\0';
		    tax_id= tax1_getTaxIdByName(b+1);

		    if(tax_id > 1) {
			*substring= StringSave(b+1);
			*e= ')';
			*Ids_out= MemNew(sizeof(Int4));
			**Ids_out= tax_id;
			return 1;
		    }
		    else if(tax_id < 0) {
			*substring= StringSave(b+1);
			*e= ')';
			return tax1_getAllTaxIdByName(*substring, Ids_out);
		    }

		    /* whole string won't help lets try truncate it at first comma*/
		    *e= ')';
		    for(e= b+1; *e != '\0'; e++) {
			if(*e == ',') {
			    *e= '\0';
			    tax_id= tax1_getTaxIdByName(b+1);

			    if(tax_id > 1) {
				*substring= StringSave(b+1);
				*e= ',';
				*Ids_out= MemNew(sizeof(Int4));
				**Ids_out= tax_id;
				return 1;
			    }
			    else if(tax_id < 0) {
				*substring= StringSave(b+1);
				*e= ',';
				return tax1_getAllTaxIdByName(*substring, Ids_out);
			    }
			    *e= ',';
			    break;
			}
		    }
		    break;
		}
	    }
	    break;
	}
    }
			
    /* we still have got nothing */
    /* try the substring before first '(' */
    
    if(*b == '(') {
	/* we are staying on the first '(' */
	*b= '\0';
	
	tax_id= tax1_getTaxIdByName(str);

	if(tax_id > 1) {
	    *Ids_out= MemNew(sizeof(Int4));
	    **Ids_out= tax_id;
	    *substring= StringSave(str);
	    *b= '(';
	    return 1;
	}
	else if(tax_id < 0) {
	    *substring= StringSave(str);
	    *b= '(';
	    return tax1_getAllTaxIdByName(*substring, Ids_out);
	}
	*b= '(';
    }

    b= StringStr(str, "Organism");
    if(b == NULL) b= StringStr(str, "organism");
    if(b == NULL) b= StringStr(str, "ORGANISM");

    if(b != NULL) {
	e= StringChr(b, ':');
	if(e != NULL) {
	    b= e+1;
	    tax_id= tax1_getTaxIdByName(b);

	    if(tax_id > 1) {
		*Ids_out= MemNew(sizeof(Int4));
		**Ids_out= tax_id;
		*substring= StringSave(b);
		return 1;
	    }
	    else if(tax_id < 0) {
		*substring= StringSave(b);
		return tax1_getAllTaxIdByName(*substring, Ids_out);
	    }
	
	    /* if multiple lines  or ; , ( */
	    for(++e; *e != '\0'; e++) {
		if((*e == '\n') || (*e == ';') || (*e == ',') || (*e == '(')) {
		    c= *e;
		    *e= '\0';
		    tax_id= tax1_getTaxIdByName(b);

		    if(tax_id > 1) {
			*substring= StringSave(b);
			*e= c;
			*Ids_out= MemNew(sizeof(Int4));
			**Ids_out= tax_id;
			return 1;
		    }
		    else if(tax_id < 0) {
			*substring= StringSave(b);
			*e= c;
			return tax1_getAllTaxIdByName(*substring, Ids_out);
		    }
		    *e= c;
		    break;
		}
	    }
	}
    }
    return 0;
}    

static int nameCmp(CharPtr s1, CharPtr s2)
{
    if((s1 == NULL) && (s2 == NULL)) return 0;
    if((s1 == NULL) || (s2 == NULL)) return 1;

    return strcmp(s1, s2);
}

static Int4 storedTaxId(OrgRefPtr orp)
{
    ValNodePtr vnp;
    DbtagPtr dbtag;
    ObjectIdPtr object_id;

    for(vnp= orp->db; vnp != NULL; vnp= vnp->next) {
	dbtag= vnp->data.ptrvalue;
	if((dbtag != NULL) && (StringCmp(dbtag->db, "taxon") == 0)) {
	    ObjectIdPtr object_id= dbtag->tag;
	    return object_id->id;
	}
    }

    return 0;
}

static Int4 OrgModCmp(OrgModPtr omp1, OrgModPtr omp2)
{
    OrgModPtr omp;
    int found;

    if(omp2 == NULL) return 0;
    if(omp1 == NULL) return 100;
    
    for(;omp2 != NULL; omp2= omp2->next) {
	found= 0;
	for(omp= omp1; omp != NULL; omp= omp->next) {
	    if((omp2->subtype == omp->subtype) &&
	       (nameCmp(omp2->subname, omp->subname) == 0)) {
		found= 1;
		break;
	    }
	}
	if(!found) return 100;
    }
    return 0;
}

static Int4 OrgRefCmp(OrgRefPtr orp1, OrgRefPtr orp2)
{
    OrgNamePtr onp1= orp1->orgname;
    OrgNamePtr onp2= orp2->orgname;

    if(onp1 == NULL) return 4;
    if(onp2 == NULL) return -2;

    if(onp1->gcode != onp2->gcode) return 2;
    if(onp1->mgcode != onp2->mgcode) return 3;

    if(nameCmp(orp1->taxname, orp2->taxname) != 0) return 10;

    if(nameCmp(onp1->lineage, onp2->lineage) != 0) return 20;

    if(nameCmp(orp1->common, orp2->common) != 0) return 30;

    if(onp1->choice != onp2->choice) return 40;

    if(nameCmp(onp1->div, onp2->div) != 0) return 50;

    return OrgModCmp(onp1->mod, onp2->mod);
    
}
 
Int4 tax1e_needUpdate(OrgRefPtr inp_orgRef)
{
    Taxon1DataPtr res;
    Int4 tax_id;
    OrgRefPtr db_orgRef;
    int is_species;
    Boolean need_search_name= TRUE;
    CharPtr hit_name;

    tax_id= tax1_getTaxIdByOrgRef(inp_orgRef);
    if(tax_id <= 0) return -1;

    if(tax_id != storedTaxId(inp_orgRef)) return 1;
    
    db_orgRef= tax1_getOrgRef(tax_id, NULL, NULL, NULL /*res->embl_code*/);
    if(db_orgRef == NULL) {
	return -2;
    }

    return OrgRefCmp(inp_orgRef, db_orgRef);
}

Boolean tax1_isAlive(void)
{
    return (tax1_getTaxIdByName("dog") > 1)? TRUE : FALSE;
}

