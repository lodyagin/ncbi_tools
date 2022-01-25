/*----------------*/

#include <stdlib.h>
#include <ncbi.h>
#include <taxinc.h>
#include <tax0.h>
#define REALTAXsyb

#define MAX_ORG_LIST	10

#define BUFF_SIZE 16
#define TAX_READ 0
#define TAX_WRITE 1

static int my_timer= -1;

static char hit_name[256];

static struct t_or_buff {
  Int4 tax_id;
  OrgRefPtr p_org_ref;
  int timer;
  char div[16];
  char embl[4];
  int is_species;
} or_buff[BUFF_SIZE];

static Boolean we_want_synonyms= 0;

Boolean tax1_setSynonyms(Boolean on_off)
{
  Boolean ret;

  ret= we_want_synonyms;
  we_want_synonyms= on_off;
  return ret;
}

static void lockBuff(int mode)
{
  mode= mode;
}

static void unlockBuff()
{
}

static void initBuff()
{
  int i;

  my_timer= 0;
  
  for(i= 0; i < BUFF_SIZE; i++) {
    or_buff[i].tax_id= 0;
    or_buff[i].p_org_ref= NULL;
  }
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

/* mark synonyms (name begins with '~') */
static void markSynonym(ValNodePtr syn)
{
  CharPtr name;
  Int2 i;

  for(;syn != NULL; syn= syn->next) {
    name= syn->data.ptrvalue;
    if(*name != '~') {
      syn->choice= 0;
    }
    else {
      syn->choice= 1;
      for(i= 0; name[i] != '\0'; i++) {
	name[i]= name[i+1];
      }
    }
  }
}

static void loadInBuff(Int4 id)
{
  int i, k= -1;
  Int4 t= my_timer + 1;
  Int4 bt;

  if(my_timer < 0) {
    initBuff();
    t= my_timer+1;
  }
  
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
    or_buff[k].p_org_ref= Tax0GetRef(-id);
    or_buff[k].timer= ++my_timer;
    or_buff[k].is_species= -1;
    if((or_buff[k].p_org_ref != NULL) && (or_buff[k].p_org_ref->syn != NULL)) {
      markSynonym(or_buff[k].p_org_ref->syn);
    }
  }
}

static Boolean getDivSpecFromOrgName(OrgNamePtr onp, int* is_spec, CharPtr div)
{
    Boolean has_div= FALSE;
    Boolean has_is_spec= FALSE;

    if(onp == NULL) return FALSE;

    if(onp->div != NULL) {
	div= StringNCpy(div, onp->div, 16);
	has_div= TRUE;
    }
    switch(onp->choice) {
    case 1: /* binomial name */
    case 2: /* virus name */
	*is_spec= 1;
	has_is_spec= TRUE;
	break;
    case 5: /* partial name */
	*is_spec= 0;
	has_is_spec= FALSE; /*TRUE;*/
	break;
    default:
	*is_spec= -1;
    }
    return has_div & has_is_spec;
}

static void loadComplete(Int4 id)
{
  int i;
  TaxCompleteListPtr tclp;
  TaxonIdName tinp;

  if(my_timer < 0) initBuff();

  lockBuff(TAX_WRITE);

  for(i= 0; i < BUFF_SIZE; i++) {
    if(or_buff[i].tax_id == id) {
      or_buff[i].timer= ++my_timer;
      if(!getDivSpecFromOrgName(or_buff[i].p_org_ref->orgname, &or_buff[i].is_species, or_buff[i].div)) {
	 
	  tinp.choice= TaxonIdName_id;
	  tinp.data.intvalue= id;
	  tclp= Tax0GetComplete(&tinp);
	  if((tclp != NULL) && (tclp->num > 0) && (tclp->info != NULL)) {
	      or_buff[i].is_species= tclp->info->is_species_level;
	      if(tclp->info->gb_div != NULL) StringNCpy(or_buff[i].div, tclp->info->gb_div, 16);
	      else or_buff[i].div[0]= '\0';
	      /*if(tclp->info->embl_code != NULL) StringNCpy(or_buff[i].embl, tclp->info->embl_code, 4);
		else*/ or_buff[i].embl[0]= '\0';
		TaxCompleteListFree(tclp);
	  }
      }
      break;
    }
  }
  unlockBuff();
}

static OrgRefPtr getFromBuff(Int4 id, int* is_sp, CharPtr div, CharPtr embl)
{
  int i;
  OrgRefPtr orp= NULL;

  if(my_timer < 0) initBuff();

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

OrgRefPtr tax1_getOrgRef(Int4 tax_id, int* is_species, CharPtr div, CharPtr embl_cde)
{
  OrgRefPtr orp;


  if((orp= getFromBuff(tax_id, is_species, div, embl_cde)) == NULL) {
    /* organism is not in buffer already */
    lockBuff(TAX_WRITE);
    loadInBuff(tax_id);
    unlockBuff();

    orp= getFromBuff(tax_id, is_species, div, embl_cde);
  }
  if((orp != NULL) && (is_species != NULL) && (*is_species == -1)) {
    /* we need, but we don't have is_species, div and embl_cde in buffer */
    loadComplete(tax_id);
    orp= getFromBuff(tax_id, is_species, div, embl_cde);
  }
  return orp;
    
}

static Int4 findIdInBuff(CharPtr orgname)
{
  int i;
  OrgRefPtr orp;
  ValNodePtr vnp;
  Int4 id= 0;

  lockBuff(TAX_READ);

  for(i= 0; i < BUFF_SIZE; i++) {
    if((or_buff[i].tax_id > 0) && (or_buff[i].p_org_ref != NULL)) {
      orp= or_buff[i].p_org_ref;
      if((orp->taxname != NULL) && (StringICmp(orp->taxname, orgname) == 0)) {
	id= or_buff[i].tax_id;
	break;
      }
      if((orp->common != NULL) && (StringICmp(orp->common, orgname) == 0)) {
	id= or_buff[i].tax_id;
	break;
      }
      for(vnp= orp->syn; vnp != NULL; vnp= vnp->next) {
	if((vnp->data.ptrvalue != NULL) && (StringICmp(vnp->data.ptrvalue, orgname) == 0)) {
	  id= or_buff[i].tax_id;
	  break;
	}
      }
      if(id != 0) break;
    }
  }
  unlockBuff();
  return id;
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
  TaxonIdListPtr tilp;
  TaxonName tn;
  Int4 tax_id;

  if((tax_id= findIdInBuff(orgname)) > 0) {
    StringNCpy(hit_name, orgname, 256);
    return tax_id;
  }

  tn.data.ptrvalue= orgname;
  tn.choice= TaxonName_taxname;
  tilp= Tax0GetTaxId(&tn);
  if(tilp == NULL) return 0;
  if(tilp->_num_ids == 0) {
    tax_id= 0;
  }
  else if(tilp->_num_ids == 1) {
    tax_id= tilp->ids[0];
  }
  else {
    tax_id= -tilp->ids[0];
  }
  TaxonIdListFree(tilp);
  StringNCpy(hit_name, orgname, 256);
  return tax_id;
}


/***************************************************
 * Get all tax_id by organism name
 * returns:
 *       Number of tax ids found
 */
Int4 tax1_getAllTaxIdByName(CharPtr orgname, Int4 **Ids_out)
{
  TaxonName tn;
  TaxonIdListPtr id_list;
  Int4 nof_ids= 0;

  tn.data.ptrvalue= orgname;
  tn.choice= TaxonName_taxname;

  id_list= Tax0GetTaxId(&tn);
  if(id_list == NULL) return 0;
  
  if(id_list->_num_ids > 0) {
    *Ids_out= MemNew(id_list->_num_ids * sizeof(Int4));
    MemCopy(*Ids_out, id_list->ids, id_list->_num_ids * sizeof(Int4));
    nof_ids= id_list->_num_ids;
  }

  TaxonIdListFree(id_list);
  return nof_ids;
}


static CharPtr getSearchName(OrgModPtr omp)
{
    while(omp != NULL) {
	if(omp->subtype == 254) {
	    return omp->subname;
	}
	omp= omp->next;
    }
    return NULL;
}

Int4 tax1_getTaxIdByOrgRef(OrgRefPtr orgRef)
{
  Int4 tax_id, id;
  ValNodePtr synonym;
  CharPtr search_name;

  tax_id= id= 0;

  if((orgRef->orgname != NULL) && (orgRef->orgname->mod != NULL) &&
     ((search_name= getSearchName(orgRef->orgname->mod)) != NULL)) {
      tax_id= tax1_getTaxIdByName(search_name);
      if(tax_id > 0) return tax_id;
      else tax_id = 0;
  }

  if(orgRef->taxname != NULL) tax_id= tax1_getTaxIdByName(orgRef->taxname);
  if(tax_id > 0) return tax_id;
  if(orgRef->common != NULL) id= tax1_getTaxIdByName(orgRef->common);
  if(id > 0) return id;

  if(tax_id == 0) tax_id= id;
  
  if(orgRef->syn != NULL) {
    id= 0;

    for(synonym= orgRef->syn; (synonym != NULL) && (id < 1); synonym= synonym->next) {
      id= tax1_getTaxIdByName(synonym->data.ptrvalue);
    }
  }

  return (id > 0)? id : tax_id;

}

static void cleanOrgName(OrgNamePtr onp)
{
  if(onp->lineage != NULL) MemFree(onp->lineage);
  /*  if(onp->mod != NULL) OrgModSetFree(onp->mod); */
  if(onp->div != NULL) MemFree(onp->div);
  if(onp->next != NULL) OrgNameSetFree(onp->next);
  if(onp->data != NULL) {
    switch(onp->choice) {
    case 1: /* binomial name */
      BinomialOrgNameFree(onp->data);
      break;
    case 2: /* virus name */
      MemFree(onp->data);
      break;
    case 5: /* partial name */
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

static Int4 extractId(DbtagPtr dbtag)
{
  ObjectIdPtr object_id;

  if((dbtag == NULL) || ((object_id= dbtag->tag) == NULL)) return 0;
  return object_id->id;
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
  object_id->id = extractId(src->db->data.ptrvalue);
    
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
  case 1: /*binomial*/
    onp->data= copyBinomial(src->orgname->data);
    break;
  case 2: /* virus */
    onp->data= (src->orgname->data != NULL)? StringSave(src->orgname->data) : NULL;
    break;
  case 5: /* partial */
    onp->data= copyPartial(src->orgname->data);
    break;
  default: /* can't handle */
    onp->data= NULL;
  }
  
  if(onp->mod == NULL) onp->mod= copyOrgMod(src->orgname->mod);
  onp->lineage= (src->orgname->lineage != NULL)? StringSave(src->orgname->lineage) : NULL;
  onp->gcode= src->orgname->gcode;
  onp->mgcode= src->orgname->mgcode;
  onp->div= (src->orgname->div != NULL) ? StringSave(src->orgname->div) : NULL;
}

static void populateReplaced(OrgRefPtr orp, CharPtr oldName)
{
  OrgNamePtr onp;
  OrgModPtr omp;

  if((orp->taxname != NULL) && (StringICmp(orp->taxname, oldName) == 0)) return;
  if((orp->common != NULL) && (StringICmp(orp->common, oldName) == 0)) return;

  /* organism name was changed */
  onp= orp->orgname;
  if((onp != NULL) && (getSearchName(onp->mod) == NULL)) {
    omp= OrgModNew();
    omp->next= onp->mod;
    omp->subtype= 254;
    omp->subname= StringSave(oldName);
    onp->mod= omp;
  }
}

Taxon1DataPtr tax1_lookup(OrgRefPtr inp_orgRef, int merge)
{
  Taxon1DataPtr res;
  Int4 tax_id;
  OrgRefPtr db_orgRef;
  int is_species;

  tax_id= tax1_getTaxIdByOrgRef(inp_orgRef);
  if(tax_id <= 0) return NULL;
  res= Taxon1DataNew();
  res->div= MemNew(16);
  res->embl_code= MemNew(4);
  db_orgRef= tax1_getOrgRef(tax_id, &is_species, res->div, res->embl_code);
  res->embl_code[0]= '\0';
  if(db_orgRef == NULL) {
    Taxon1DataFree(res);
    return NULL;
  }

  res->is_species_level= is_species;
  if(merge) {
    /* we have to merge old orgref with new one */
    res->org= inp_orgRef;
    /* clean-up old information */
    if(inp_orgRef->taxname != NULL) MemFree(inp_orgRef->taxname);
    if(inp_orgRef->common != NULL) MemFree(inp_orgRef->common);
    if(inp_orgRef->syn != NULL) ValNodeFreeData(inp_orgRef->syn);
    if(inp_orgRef->orgname != NULL) cleanOrgName(inp_orgRef->orgname);
  }
  else {
    /* make new orgref */
    res->org= OrgRefNew();
    res->org->db= NULL;
    res->org->orgname= NULL;
  }
  /* fill-up orgref based on db_orgRef */
  bldOrgRefOut(res->org, db_orgRef, tax_id);
  populateReplaced(res->org, hit_name);
  return res;
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
  db_orgRef= tax1_getOrgRef(tax_id, &is_species, res->div, res->embl_code);
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
  bldOrgRefOut(res->org, db_orgRef, tax_id);
  return res;
}
  
Boolean tax1_init()
{
  initBuff();
  return Tax0Init();
}

void tax1_fini()
{
  Tax0Fini();
}
  
Int4 tax1_getParent(Int4 id_tax)
{
    TaxonIdListPtr tilp;
    Int4 res= 0;

    if((tilp= Tax0GetParents(id_tax)) != NULL) {
	if(tilp->_num_ids > 0) res= tilp->ids[0];
	TaxonIdListFree(tilp);
    }

    return res;
}

int tax1_getChildren(Int4 id_tax, Int4** ids_out)
{
    int n= 0;
    TaxonIdListPtr tilp;
    
    *ids_out= NULL;
    if((tilp= Tax0GetChildren(id_tax)) != NULL) {
	if(tilp->_num_ids > 0) *ids_out= tilp->ids;
	n= tilp->_num_ids;
	tilp->_num_ids= 0;
	tilp->ids= NULL;
	TaxonIdListFree(tilp);
    }
    return n;
}

/* find the nearest ancestor for two nodes */
Int4 tax1_join(Int4 taxid1, Int4 taxid2)
{
    TaxonIdListPtr tilp1, tilp2;
    Int4 res= 0;
    Int2 i, j;
    
    if(taxid1 == taxid2) return taxid1;

    tilp1= Tax0GetParents(-taxid1);
    if(tilp1 == NULL) return 1;
    if(tilp1->_num_ids < 1) {
	TaxonIdListFree(tilp1);
	return 1;
    }

    for(i= 0; i < tilp1->_num_ids; i++) {
	if(tilp1->ids[i] == taxid2) {
	    TaxonIdListFree(tilp1);
	    return taxid2;
	}
    }

    tilp2= Tax0GetParents(-taxid2);
    if(tilp2 == NULL) {
	TaxonIdListFree(tilp1);
	return 1;
    }
    if(tilp2->_num_ids >= 1) {
	for(i= 0; i < tilp2->_num_ids; i++) {
	    if(tilp2->ids[i] == taxid1) {
		res= taxid1;
		break;
	    }
	}

	if(res == 0) {
	    for(i= 0; i < tilp1->_num_ids; i++) {
		taxid1= tilp1->ids[i];
		for(j= 0; j < tilp2->_num_ids; j++) {
		    if(tilp2->ids[j] == taxid1) {
			res= taxid1;
			break;
		    }
		}
		if(res) break;
	    }
	}

    }

    TaxonIdListFree(tilp2);
    TaxonIdListFree(tilp1);
    return ((res > 0)? res : 1);
}
		    
		    


