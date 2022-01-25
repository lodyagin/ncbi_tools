/*----------------*/

#include <stdlib.h>
#include <ncbi.h>
#include "taxinc.h"
#include "tax_cmmn.h"
#define REALTAXsyb

#define MAX_ORG_LIST	10

#define BUFF_SIZE 16
#define TAX_READ 0
#define TAX_WRITE 1

static Int2 VRL_div= 0;
static Int2 PHG_div= 0;

static _taxDBCtlPtr DBctl= NULL;
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

static char hit_name[256];

static struct t_or_buff {
  Int4 tax_id;
  OrgRefPtr p_org_ref;
  int timer;
  char div[16];
  char embl[4];
  int is_species;
} or_buff[BUFF_SIZE];

static CharPtr DB_PATH= "/netopt/genbank/subtool/taxdb/";

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

static Int4 getLiveId(Int4 id)
{
  if((DBctl->tree == NULL) || (id < 0) || (id > tax_getMaxId(DBctl))) return 0;

  if(DBctl->tree[id].parent < 1) {
    /* deleted node */
    return 0;
  }

  while(DBctl->tree[id].child == -1) {
    id=DBctl->tree[id].parent;
  }

  return id;
}

static void msgCallBack(int msg_type, int msg_code, CharPtr msg_txt)
{
  if(msg_type == TAX_ERROR) {
    ErrPostEx(SEV_ERROR, 100, msg_code, msg_txt);
  }
  else {
    ErrPostEx(SEV_WARNING, 100, msg_code, msg_txt);
  }
}

/**************************************************************************
 *
 *	InitTaxDB
 *
 **************************************************************************/

InitTaxDB()
{
  CharPtr tmp;

  if((tmp=getenv("TAXDBPATH")) != NULL) DB_PATH= tmp;

  tax_setMsgHndl(msgCallBack);
  DBctl= tax_loadDBCtl(DB_PATH, "dbctl.tax");
  if(DBctl == NULL) return 0;
  if(tax_loadClasses(DBctl) || tax_loadRanks(DBctl) || tax_loadDivs(DBctl) ||
     tax_loadGCs(DBctl) || tax_loadTree(DBctl) || tax_loadNames(DBctl)) {
    return 0;
  }

  tax_getNode(DBctl, -1);

  SpeciesRank= tax_getRankId(DBctl, "species");
  SubspeciesRank= tax_getRankId(DBctl, "subspecies");
  GenusRank= tax_getRankId(DBctl, "genus");
  FamilyRank= tax_getRankId(DBctl, "family");
  OrderRank= tax_getRankId(DBctl, "order");
  ClassRank= tax_getRankId(DBctl, "class");


  VRL_div= tax_getDivId(DBctl, "VRL", TAX_DIV_CDE);
  PHG_div= tax_getDivId(DBctl, "PHG", TAX_DIV_CDE);

  SYNONYM= tax_getClassId(DBctl, "synonym");
  COMMON_NAME= tax_getClassId(DBctl, "common name");
  PREF_COMMON= tax_getClassId(DBctl, "preferred common name");
  initBuff();
  return 1;
}

/**************************************************************************
 *
 *	CloseTaxDB
 *
 **************************************************************************/

void CloseTaxDB()
{
  int i;

  if(DBctl != NULL) {
    tax_freeNames(DBctl);
    tax_freeTree(DBctl);
    tax_freeGCs(DBctl);
    tax_freeDivs(DBctl);
    tax_freeClasses(DBctl);
    tax_freeRanks(DBctl);
   
    for(i= 0; i < BUFF_SIZE; i++) {
      if(or_buff[i].p_org_ref != NULL) {
	OrgRefFree(or_buff[i].p_org_ref);
      }
    }
  }
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

static OrgModPtr bldOrgMod(Int4 id, _taxNamePtr o_name)
{
  Int4 p_id;
  Int2 rank;
  _taxNamePtr parent;
  OrgModPtr orgMdf= OrgModNew();
  CharPtr name= o_name->name_txt;

  for(p_id= DBctl->tree[id].parent; p_id > 1; p_id= DBctl->tree[p_id].parent) {
    rank= DBctl->tree[p_id].rank;
    if((rank == SubspeciesRank) || 
       (rank == SpeciesRank) ||
       (rank == GenusRank)) break;
  }

  parent= (p_id > 1)? tax_findOrgName(DBctl, p_id) : NULL;

  if(parent != NULL) {
    orgMdf->subname= StringSave(strTail(name, parent->name_txt));
  }
  else {
    orgMdf->subname= StringSave(name);
  }

  rank= DBctl->tree[id].rank;

  if(rank == SubspeciesRank) {
    orgMdf->subtype= 22; /* subspecies */
  }
  else if(rank == tax_getRankId(DBctl, "varietas")) {
    orgMdf->subtype= 6; /* variety */
  }
  else if(rank == tax_getRankId(DBctl, "forma")) {
    orgMdf->subtype= 2; /* strain */
  }
  else if((parent != NULL) && (DBctl->tree[p_id].rank == SubspeciesRank)) {
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
static int binomialName(Int4 id, _taxNamePtr o_name, OrgNamePtr onp)
{
  int i;
  Int4 genusId= 0;
  Int4 speciesId= 0;
  Int4 subspeciesId= 0;
  Int4 p_id;
  Int2 rank;
  BinomialOrgNamePtr bName;
  _taxNamePtr genus, species, subspecies;


  for(p_id= id; p_id > 1; p_id= DBctl->tree[p_id].parent) {
    rank= DBctl->tree[p_id].rank;
    if(rank == SpeciesRank) speciesId= p_id;
    else if(rank == SubspeciesRank) subspeciesId= p_id;
    else if(rank == GenusRank) genusId= p_id;
  }

  if(genusId == 0) {
    /* try to find subgenus if genus anavalable */
    for(p_id= id; p_id > 1; p_id= DBctl->tree[p_id].parent) {
      if(DBctl->tree[p_id].rank == (GenusRank + 1)) genusId= p_id;
    }
  }

  genus= (genusId == 0)? NULL : tax_findOrgName(DBctl, genusId);
  species= (speciesId == 0)? NULL : tax_findOrgName(DBctl, speciesId);
  subspecies= (subspeciesId == 0)? NULL : tax_findOrgName(DBctl, subspeciesId);

  if(genus == NULL) return 0; /* no genus - no binomial */

  onp->choice= 1; /*binomial*/

  onp->data= bName= BinomialOrgNameNew();
  
  bName->genus= StringSave(genus->name_txt);

  if(species != NULL) {
    /* we have a species in lineage */
    bName->species= StringSave(strTail(species->name_txt, genus->name_txt));
    if(subspecies != NULL) {
      /* we also have a subspecies in lineage */
      bName->subspecies= StringSave(strTail(subspecies->name_txt, species->name_txt));
    }
    else {
      bName->subspecies= NULL;
    }
    onp->mod= (id == speciesId)? NULL : bldOrgMod(id, o_name);    
    return 1;
  }
  
  /* no species in lineage */

  if(subspecies != NULL) {
    /* we have no species but we have subspecies */
    bName->species= NULL;
    bName->subspecies= StringSave(strTail(subspecies->name_txt, genus->name_txt));
    onp->mod= bldOrgMod(id, o_name);
    return 1;
  }
  
  /* we have no species, no subspecies but we are under species level (varietas or forma) */

  bName->species= NULL;
  bName->subspecies= NULL;
  onp->mod= bldOrgMod(id, o_name);
  return 1;
}

      
static void partialName(Int4 id, _taxNamePtr o_name, OrgNamePtr onp)
{

  TaxElementPtr taxElem;
  Int2 rank_id= DBctl->tree[id].rank;

  onp->choice= 5; /* partial */
  onp->data= taxElem= TaxElementNew();
  
  if(rank_id == FamilyRank) {
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
    taxElem->level= StringSave(tax_getRankById(DBctl, rank_id));
  }

  taxElem->name= StringSave(o_name->name_txt);
  taxElem->next= NULL;
}

  
/*****************************************************************
 * build synonyms valnodes
 * this routine include in valnodes synonyms and common synonyms
 */
static ValNodePtr bldSynValNodes(Int4 id, _taxNamePtr syn)
{
  ValNodePtr list= NULL;
  ValNodePtr header= NULL;
  int i;

  for(i= 0; syn[i].tax_id == id; i++) {
    if((syn[i].class_cde == SYNONYM) || (syn[i].class_cde == COMMON_NAME)) {
      list= ValNodeNew(list);
      list->choice= (syn[i].class_cde == SYNONYM)? 1 : 0;
      list->data.ptrvalue= StringSave(syn[i].name_txt);
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

static CharPtr bldLineage(Int4 id)
{
  ValNodePtr head=NULL, this=NULL;
  int len = 2;
  CharPtr temp, retval = NULL;
  Int4 lineage[64];
  int n, i;
  _taxNamePtr orgName;
  
  n= tax_getLin(DBctl, id, lineage);

  for(i= 1; i < n; i++) {
    if((DBctl->tree[lineage[i]].flags & GB_HIDDEN) != 0) continue;

    if((DBctl->tree[lineage[i]].rank >= SpeciesRank) && head) break;

    this= ValNodeNew(this);
    if (head == NULL){
      head = this;
    }

    if((orgName= tax_findOrgName(DBctl, lineage[i])) != NULL) {
	this->data.ptrvalue= StringSave(orgName->name_txt);
	len+= StringLen(orgName->name_txt) + 3;
    }
    else {
	ErrPostEx(SEV_ERROR, 100, lineage[i], "Error in taxonomy lineage");
    }	
  }
  
  /* all data into linked list, now format it */
  if (len > 2){
    temp = retval = MemNew(len);
    for (this = head; this; this = this -> next){
      if (this != head){
	temp = StringMove(temp, "; " );
      }
      temp = StringMove(temp, this -> data.ptrvalue);
    }
    ValNodeFreeData(head);
  }

  if(retval == NULL) {
    /* empty lineage */
    if((orgName= tax_findOrgName(DBctl, id)) != NULL) 
      retval= StringSave(orgName->name_txt);
  }

  return retval;
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

static OrgNamePtr bldOrgName(Int4 id, _taxNamePtr o_name, int* is_species_out, CharPtr div, CharPtr embl)
{
  OrgNamePtr onp;
  _taxNodePtr node;
  Int2 rank_id= DBctl->tree[id].rank;
  int is_species;
  Int4 p_id, s_id;
  CharPtr div_abbr;
  Int2 div_id;
  _taxNamePtr s_name;

  onp= OrgNameNew();
  node= tax_getNode(DBctl, id);
  if(node == NULL) return NULL;
  
  onp->gcode= node->gc;
  onp->mgcode= node->mgc;
  onp->lineage= bldLineage(id);
  if(embl != NULL) StringCpy(embl, ""/*node->embl_cde*/);
	   
  is_species= (rank_id >= SpeciesRank)? 1 : 0;
  /* correct level by lineage if node has no rank */
  if(rank_id < 0) {
    p_id= id;
    while(p_id > 1) {
      p_id= DBctl->tree[p_id].parent;
      if(DBctl->tree[p_id].rank >= SpeciesRank) {
	is_species= 1;
	break;
      }
    }
  }

  div_id= node->division;
  tax_freeNode(node);
  onp->div= StringSave(tax_getDivById(DBctl, div_id, TAX_DIV_CDE));
  StringCpy(div, tax_getDivById(DBctl, div_id, TAX_DIV_CDE));
  *is_species_out= is_species;

  if(is_species) {
    /* we are on species level or below */
	     
    /* check for viruses */
    if((div_id == VRL_div) || (div_id == PHG_div)) {
      /* this is a virus */
      onp->choice= 2; /* virus */
      if(rank_id == SpeciesRank) {
	/* we are on species level */
	onp->data= StringSave(o_name->name_txt);
	onp->mod= NULL;
      }
      else {
	/* we are below species */
	/* first try to find species or min rank which below species but above us */
	p_id= id; s_id= 0;
	while(p_id > 1) {
	  p_id= DBctl->tree[p_id].parent;
	  if(DBctl->tree[p_id].rank == SpeciesRank) break;
	  else if(DBctl->tree[p_id].rank > SpeciesRank) s_id= p_id;
	}
	if(p_id <= 1) p_id= s_id;

	if(p_id > 1) {
	  /* we below species but we have species above us */
	  s_name= tax_findOrgName(DBctl, p_id);
	  if(s_name != NULL) {
	    onp->data= StringSave(s_name->name_txt);
	    onp->mod= bldOrgMod(id, o_name);
	  }
	}
	else {
	  /* we probably below species but no species above us */
	  onp->data= StringSave(o_name->name_txt);
	  onp->mod= bldOrgMod(id, o_name);
	}
      }
    }
    else if(!binomialName(id, o_name, onp)) {
      /* name is not binomial: set partial */
      partialName(id, o_name, onp);
    }
  }
  else {
    /* above species */
    partialName(id, o_name, onp);
  }

  return onp;
}


static void bldOrgRef(Int4 id, OrgRefPtr orp, int* is_species, CharPtr div, CharPtr embl)
{
  _taxNamePtr o_name;
  int i;

  *is_species= 0;
  *div= *embl= '\0';
  o_name= tax_findOrgName(DBctl, id);
  if(o_name == NULL) return;
  orp->taxname= StringSave(o_name->name_txt);

  /* fill-up preferred common name */
  orp->common= NULL;
  for(i= 1; o_name[i].tax_id == id; i++) {
    if(o_name[i].class_cde == PREF_COMMON) {
      orp->common= StringSave(o_name[i].name_txt);
      break;
    }
  }

  /* fill-up synonyms */
  orp->syn= bldSynValNodes(id, o_name);
  orp->mod= NULL;
  orp->db= bldDBId(id);
  orp->orgname= bldOrgName(id, o_name, is_species, div, embl);
}
  

static void loadInBuff(Int4 id)
{
  int i, k= -1;
  Int4 t= my_timer + 1;
  Int4 bt;
  
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
    bldOrgRef(id, or_buff[k].p_org_ref, &or_buff[k].is_species, or_buff[k].div, or_buff[k].embl);
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

/***************************************************
 * Get tax_id by organism name
 * returns:
 *       tax_id if one node found
 *       0      no organism found
 *       -tax_id if more than one node found
 */
Int4 tax1_getTaxIdByName(CharPtr orgname)
{
  int n, i, j, f;
  Int4 tax_id, spec_id;
  _taxNamePtr* nameList;

#if 0
  if((n= tax_findByName(DBctl, orgname, TAX_NAME_SEARCH, &nameList)) == 0) {
    /* do not find exact name */
    if((n= tax_findByName(DBctl, orgname, TAX_RE_SEARCH, &nameList)) == 0) {
      /* do not find again */
      if((n= tax_findByName(DBctl, orgname, TAX_TOKEN_SEARCH, &nameList)) == 0) {
	/* no such name */
	return 0;
      }
    }
  }
#else
  if((n= tax_findByName(DBctl, orgname, TAX_NAME_SEARCH, &nameList)) == 0) {
    return 0;
  }
#endif

  tax_id= nameList[0]->tax_id;
  for(i= 1; i < n; i++) {
    if(tax_id != nameList[i]->tax_id) {
      tax_id= tax_getDesignator(orgname);
      break;
    }
  }
  
  if(tax_id == 0) {
    /* try to find species */
    tax_id= nameList[0]->tax_id;
    spec_id= (DBctl->tree[tax_id].rank == SpeciesRank)? tax_id : 0;
    for(i= 1; i < n; i++) {
      if(tax_id != nameList[i]->tax_id) {
	tax_id= nameList[i]->tax_id;
	if(DBctl->tree[tax_id].rank == SpeciesRank) {
	  if(spec_id == 0) spec_id= tax_id;
	  else if(spec_id != tax_id) {
	    spec_id= -spec_id;
	    break;
	  }
	}
      }
    }
    tax_id= (spec_id == 0)? -nameList[0]->tax_id : spec_id;
  }
      

  MemFree(nameList);
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
  int n, i, j, f;
  _taxNamePtr* nameList;
  Int4 *Ids;

#if 0
  if((n= tax_findByName(DBctl, orgname, TAX_NAME_SEARCH, &nameList)) == 0) {
    /* do not find exact name */
    if((n= tax_findByName(DBctl, orgname, TAX_RE_SEARCH, &nameList)) == 0) {
      /* do not find again */
      if((n= tax_findByName(DBctl, orgname, TAX_TOKEN_SEARCH, &nameList)) == 0) {
	/* no such name */
	*Ids_out= NULL;
	return 0;
      }
    }
  }
#else
  if((n= tax_findByName(DBctl, orgname, TAX_NAME_SEARCH, &nameList)) == 0) {
    *Ids_out= NULL;
    return 0;
  }
#endif

  *Ids_out= Ids= MemNew(n*sizeof(Int4));
  
  Ids[0]= nameList[0]->tax_id;
  j= 0;
  for(i= 1; i < n; i++) {
    if(Ids[j] != nameList[i]->tax_id) {
      Ids[++j]= nameList[i]->tax_id;
    }
  }

  MemFree(nameList);
  return ++j;
}

/***************************************************
 * Get tax_id by organism name
 * returns:
 *       tax_id if one node found
 *       0      no organism found
 *       -tax_id if more than one node found
 */
Int4 tax1_findTaxIdByName(CharPtr orgname)
{
  int n, i, j, f;
  Int4 tax_id, spec_id;
  _taxNamePtr* nameList;

#if 1
  if((n= tax_findByName(DBctl, orgname, TAX_ALPHANUM_SEARCH, &nameList)) == 0) {
      CharPtr tail;
      if((tail= StringChr(orgname, '<')) != NULL) {
	  *tail= '\0';
	  if((n= tax_findByName(DBctl, orgname, TAX_ALPHANUM_SEARCH, &nameList)) != NULL) {
	      *tail= '<';
	      for(i= 0; i < n; i++) {
		  if((nameList[i]->unique_name != NULL) && 
		     (StringICmp(nameList[i]->unique_name, orgname) == 0)) {
		      tax_id= nameList[i]->tax_id;
		      StringNCpy(hit_name, orgname, 256);
		      return tax_id;
		  }
	      }
	  }
	  *tail= '<';
	  return 0;
      }
      /* do not find exact name */
      if((n= tax_findByName(DBctl, orgname, TAX_RE_SEARCH, &nameList)) == 0) {
	  /* do not find again */
	  if((n= tax_findByName(DBctl, orgname, TAX_TOKEN_SEARCH, &nameList)) == 0) {
	      /* no such name */
	      return 0;
	  }
      }
  }
#else
  if((n= tax_findByName(DBctl, orgname, TAX_NAME_SEARCH, &nameList)) == 0) {
    return 0;
  }
#endif

  tax_id= nameList[0]->tax_id;
  for(i= 1; i < n; i++) {
    if(tax_id != nameList[i]->tax_id) {
      tax_id= tax_getDesignator(orgname);
      break;
    }
  }
  
  if(tax_id == 0) {
    /* try to find species */
    tax_id= nameList[0]->tax_id;
    spec_id= (DBctl->tree[tax_id].rank == SpeciesRank)? tax_id : 0;
    for(i= 1; i < n; i++) {
      if(tax_id != nameList[i]->tax_id) {
	tax_id= nameList[i]->tax_id;
	if(DBctl->tree[tax_id].rank == SpeciesRank) {
	  if(spec_id == 0) spec_id= tax_id;
	  else if(spec_id != tax_id) {
	    spec_id= -spec_id;
	    break;
	  }
	}
      }
    }
    tax_id= (spec_id == 0)? -nameList[0]->tax_id : spec_id;
  }
      

  MemFree(nameList);
  StringNCpy(hit_name, orgname, 256);
  return tax_id;
}


/***************************************************
 * Get all tax_id by organism name
 * returns:
 *       Number of tax ids found
 */
Int4 tax1_findAllTaxIdByName(CharPtr orgname, Int4 **Ids_out)
{
  int n, i, j, f;
  _taxNamePtr* nameList;
  Int4 *Ids;

#if 1
  if((n= tax_findByName(DBctl, orgname, TAX_ALPHANUM_SEARCH, &nameList)) == 0) {
      CharPtr tail;
      if((tail= StringChr(orgname, '<')) != NULL) {
	  *tail= '\0';
	  if((n= tax_findByName(DBctl, orgname, TAX_ALPHANUM_SEARCH, &nameList)) != NULL) {
	      *tail= '<';
	      for(i= 0; i < n; i++) {
		  if((nameList[i]->unique_name != NULL) && 
		     (StringICmp(nameList[i]->unique_name, orgname) == 0)) {
		      *Ids_out= MemNew(sizeof(Int4));
		      **Ids_out= nameList[i]->tax_id;
		      return 1;
		  }
	      }
	  }
	  *tail= '<';
	  *Ids_out= NULL;
	  return 0;
      }
      /* do not find exact name */
      if((n= tax_findByName(DBctl, orgname, TAX_RE_SEARCH, &nameList)) == 0) {
	  /* do not find again */
	  if((n= tax_findByName(DBctl, orgname, TAX_TOKEN_SEARCH, &nameList)) == 0) {
	      /* no such name */
	      *Ids_out= NULL;
	      return 0;
	  }
      }
  }
#else
  if((n= tax_findByName(DBctl, orgname, TAX_NAME_SEARCH, &nameList)) == 0) {
    *Ids_out= NULL;
    return 0;
  }
#endif

  *Ids_out= Ids= MemNew(n*sizeof(Int4));
  
  Ids[0]= nameList[0]->tax_id;
  j= 0;
  for(i= 1; i < n; i++) {
    if(Ids[j] != nameList[i]->tax_id) {
      Ids[++j]= nameList[i]->tax_id;
    }
  }

  MemFree(nameList);
  return ++j;
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
  if(onp->div != NULL) MemFree(onp->div);
#if 0
  /* #if 0 means that we will trust to initial modifier */
  if(onp->mod != NULL) OrgModSetFree(onp->mod);
#endif
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
  onp->div= StringSave(src->orgname->div);
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
  db_orgRef= tax1_getOrgRef(tax_id, &is_species, res->div, NULL /*res->embl_code*/);
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
  
Boolean tax1_init()
{
  return InitTaxDB();
}

void tax1_fini()
{
  CloseTaxDB();
}
  
Int4 tax1_getParent(Int4 id_tax)
{
  return (id_tax == 1)? 0 : tax_getParent(DBctl, id_tax);
}

int tax1_getChildren(Int4 id_tax, Int4** ids_out)
{
  int n, i;
  Int4* ids;
  Int4 id;

  for(n= 0, id= tax_getChild(DBctl, id_tax); id > 0; id= tax_getSibling(DBctl, id)) {
    n++;
  }

  if(n > 0) {
    ids= MemNew(n*sizeof(Int4));
    for(i= 0, id= tax_getChild(DBctl, id_tax); id > 0; id= tax_getSibling(DBctl, id)) {
      ids[i++]= id;
    }
    *ids_out= ids;
  }
  return n;
}

CharPtr tax1_getGCName(Int2 gc_id)
{
  return tax_getGCById(DBctl, gc_id, TAX_GC_NM);
}

/* find last common ancestor for two nodes */
Int4 tax1_join(Int4 taxid1, Int4 taxid2)
{
    Int4 lin1[64];
    Int2 i= 1, j;

    lin1[0]= taxid1;
    if(taxid1 > 1) {
	for(i= 1; i < 64; i++) {
	    lin1[i]= tax_getParent(DBctl, lin1[i-1]);
	    if(lin1[i] <= 1) break;
	    if(lin1[i] == taxid2) return taxid2;
	}
    }

    while(taxid2 > 1) {
	for(j= 0; j < i; j++) {
	    if(taxid2 == lin1[j]) return taxid2;
	}
	taxid2= tax_getParent(DBctl, taxid2);
    }

    return 1;
}

Int2 tax1_getAllNames(Int4 tax_id, CharPtr **out_names, Boolean unique)
{
    CharPtr* names;
    _taxNamePtr org_n; 
    Int2 k, i;

    *out_names= NULL;

    org_n= tax_findOrgName(DBctl, tax_id);
    if(org_n == NULL) return 0;

    for(k= 0; (k <= tax_getNofNames(DBctl)) && (org_n[k].tax_id == tax_id); k++);

    *out_names= names= MemNew(k*sizeof(CharPtr));

    for(i= 0; i < k; i++) {
	if(unique) {
	    names[i]= StringSave((org_n[i].unique_name == NULL)?
				 org_n[i].name_txt : org_n[i].unique_name);
	}
	else {
	    names[i]= StringSave(org_n[i].name_txt);
	}
    }

    return k;
}
	    
Int4 tax1_getTaxId4Str(CharPtr str, CharPtr* substring, Int4Ptr *Ids_out)
{
    CharPtr b, e;
    Int4 n, tax_id;
    Int4Ptr Ids;
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
    

