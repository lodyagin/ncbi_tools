/*****   taxprocs.joel.c   *****/
/*------------

old name                    Joel name
-------------              ---------------
tn_get_taxid                  get_tax_id 
   SELECT   tax_id, (int)
                name_txt, (varchar(96))
                class_cde (int)
   FROM  #ls
   HAVING   pri = min(pri)




ts_get_children               get_children 
   SELECT   DISTINCT name_txt,
               tax_id
   FROM     #tax_get
   WHERE tax_id != @search_tax_id
   ORDER BY name_txt


ts_get_parents                get_parent_tax_id
   SELECT   parent_tax_id
   FROM  tax_node
   WHERE tax_id = @search_tax_id


tn_get_orgref                 get_names
   SELECT   synonym.name_txt,
		    class.class_cde,
            class.class_txt, (varchar(30))
      class.priority_no (int)
   FROM     synonym, class
   WHERE    synonym.tax_id = @search_tax_id AND
      synonym.class_cde = class.class_cde
   ORDER BY class.priority_no


tn_get_lineage_by_taxid       get_lin
 SELECT name_txt,
   tax_id
 FROM   #tax
 ORDER BY tax_level_no DESC


ts_get_gc                     get_tax_info ?
tn_get_taxinfo                get_tax_info

   SELECT   gc_id, (smallint)
      mgc_id, (smallint)
      gc_nm, (varchar(128))
      mgc_nm, (varchar(128))
      div_id,(smallint)
      div_cde, char(3)
      div_txt, varchar(64)
      rank_id, smallint
      rank_txt, (varchar(32))
      embl_cde, char(4)
      inherit_gc_ind, tinyint
      inherit_mgc_ind, tinyint
      inherit_div_ind, tinyint
      tax_id,
      parent_tax_id,
      tl_ind, smallint
      name_txt varchar(96)
   FROM  #info

the new procedures are in getorg.sp
in /export/home/plotkin/tax/stored_procs
on peony...

----------------*/

#include <stdlib.h>
#include <ncbi.h>
#include <ncbinet.h>
#include <accentr.h>
#include <objfeat.h>
#include <objtaxon.h>
#define REALTAXsyb
#include <tax_cmmn.h>

#ifdef  SYSV
/* This is a quicky BSD -> SYSV conversion pack */
#define bcopy(s,d,n)    memmove(d,s,n)
#define bzero(d,n)      memset(d,0,n)
#define bcmp(s1,s2,n)   (memcmp(s1,s2,n) ? 1 : 0)
#define index(s,c)      strchr(s,c)
#define rindex(s,c)     strrchr(s,c)
#endif

#define MAX_ORG_LIST	10

#define AssignValue( s ) ( s[0]=='\0' ? NULL : s )

static _taxDBCtlPtr DBctl;
static int IsSpeciesRank= 26;
static Int2 SpeciesRank= 26;
static Int2 SubspeciesRank= 27;
static Int2 GenusRank= 22;
static Int2 FamilyRank= 0;
static Int2 OrderRank= 0;
static Int2 ClassRank= 0;
static Int2 SYNONYM= 0;
static Int2 COMMON_NAME= 0;
static Int2 PREF_COMMON= 0;

static Int2 VRL_div= 0;
static Int2 PHG_div= 0;

static CharPtr DB_PATH= "/netopt/genbank/subtool/taxdb/";


/**************************************************************************
 *
 *	InitTaxDB
 *
 **************************************************************************/

InitTaxDB()
{
  CharPtr tmp;

  if((tmp=getenv("TAXDBPATH")) != NULL) DB_PATH= tmp;

  DBctl= tax_loadTaxDB(DB_PATH, "dbctl.tax");
  if(DBctl == NULL) return 0;
  if(tax_loadTree(DBctl) || tax_loadNames(DBctl)) {
    return 0;
  }

  IsSpeciesRank= tax_getRankId(DBctl, "species");

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
  return 1;
}

/**************************************************************************
 *
 *	CloseTaxDB
 *
 **************************************************************************/

CloseTaxDB()
{
   tax_freeNames(DBctl);
   tax_freeTree(DBctl);
   tax_freeGCs(DBctl);
   tax_freeDivs(DBctl);
   tax_freeClasses(DBctl);
   tax_freeRanks(DBctl);
}

/***************************************************************************
 *
 *  function:	CopyString
 *
 ***************************************************************************/
static
CopyString( CharPtr new, CharPtr text )
{
   CharPtr p, np;
   Char quote = '"';

   for( p = text, np = new; *p; *np++ = *p++  )
	if( *p == quote ) *np++ = *p;

   *np = '\0';
}


/**************************************************************************
 *
 *	FsTaxGetTaxId
 *
 **************************************************************************/
TaxonIdListPtr FsTaxGetTaxId(TaxonNamePtr tnp)
{
  int n, i, j, f;
  Int4 tax_id;
  _taxNamePtr* nameList;
  CharPtr orgname;
  TaxonIdListPtr tilp = NULL;  

  orgname= tnp->data.ptrvalue;

#if 0
  if((n= tax_findByName(DBctl, orgname, TAX_NAME_SEARCH, &nameList)) == 0) {
    /* do not find exact name */
    if((n= tax_findByName(DBctl, orgname, TAX_RE_SEARCH, &nameList)) == 0) {
      /* do not find again */
      if((n= tax_findByName(DBctl, orgname, TAX_TOKEN_SEARCH, &nameList)) == 0) {
	/* no such name */
#ifdef WARN_ON_FAILURE
	ErrPostEx(SEV_WARNING,1,2,"FsTaxGetTaxId: Not found: <%s>", orgname);
#endif
	return NULL;
      }
    }
  }
#else
  if((n= tax_findByName(DBctl, orgname, TAX_NAME_SEARCH, &nameList)) == 0) {
    /* no such name */
#ifdef WARN_ON_FAILURE
    ErrPostEx(SEV_WARNING,1,2,"FsTaxGetTaxId: Not found: <%s>", orgname);
#endif
    return NULL;
  }
#endif

  tilp= TaxonIdListNew();
  tilp->_num_ids= 0;
  tilp->ids= MemNew(n*sizeof(Int4));

  for(i= 0; i != n; i++) {
    tax_id= nameList[i]->tax_id;
    /* check for duplicate */
    f= 0;
    for(j= 0; j < tilp->_num_ids; j++) {
      if(tilp->ids[j] == tax_id) {
	f= 1;
	break;
      }
    }
    if(!f) {
      tilp->ids[tilp->_num_ids++]= tax_id;
    }
  }
  if(tilp->_num_ids > 1) {
    tax_id= tax_getDesignator(orgname);
    if(tax_id > 0) {
      tilp->_num_ids= 1;
      tilp->ids[0]= tax_id;
    }
    else {
      tax_id= 0;
      for(i= 0; i < tilp->_num_ids; i++) {
	if(DBctl->tree[tilp->ids[i]].rank == SpeciesRank) {
	  if(tax_id == 0) tax_id= tilp->ids[i];
	  else if(tax_id != tilp->ids[i]) {
	    tax_id= 0;
	    break;
	  }
	}
      }
      if(tax_id != 0) {
	tilp->_num_ids= 1;
	tilp->ids[0]= tax_id;
      }
    }
  }
      
  MemFree(nameList);
  return tilp;
}

static Int4 getIdById(Int4 id)
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
      
/**************************************************************************
 *
 *	FsTaxGetChildren
 *
 **************************************************************************/

TaxonIdListPtr FsTaxGetChildren( int id_tax )
{
  TaxonIdListPtr tilp = NULL;
  int n;
  Int4 id;


  if((id= getIdById(id_tax)) == 0) {
#ifdef WARN_ON_FAILURE
    ErrPostEx(SEV_WARNING,1,2,"FsTaxGetChildren: Not found: <%d>", id_tax);
#endif
    return NULL;
  }

  id_tax= id;
  
  for(n= 0, id= tax_getChild(DBctl, id); id > 0; id= tax_getSibling(DBctl, id)) {
    if((DBctl->tree[id].flags & TAX_HIDDEN) == 0) n++;
  }

  if(n > 0) {
    tilp= TaxonIdListNew();
    tilp->_num_ids= 0;
    tilp->ids= MemNew(n*sizeof(Int4));
    for(id= tax_getChild(DBctl, id_tax); id > 0; id= tax_getSibling(DBctl, id)) {
      if((DBctl->tree[id].flags & TAX_HIDDEN) == 0) {
	tilp->ids[tilp->_num_ids++]= id;
      }
    }
  }
  return tilp;
}


/**************************************************************************
 *
 *	FsTaxGetParents
 *
 **************************************************************************/


TaxonIdListPtr FsTaxGetParents( int id_tax )
{
  TaxonIdListPtr tilp = NULL;
  Int4 id;
  Int2 n= 1;  


  if(id_tax < 0) {
      id_tax= -id_tax;
      n= 0;
  }

  if((id= getIdById(id_tax)) == 0) {
#ifdef WARN_ON_FAILURE
    ErrPostEx(SEV_WARNING,1,2,"FsTaxGetParents: Not found: <%d>", id_tax);
#endif
    return NULL;
  }

  tilp= TaxonIdListNew();
  if(n == 0) {
      Int2 i;

      for(id_tax= id; id_tax > 1; id_tax= tax_getParent(DBctl, id_tax)) n++;
      
      if(n == 0) n= 1;
      tilp->_num_ids= n;
      tilp->ids= MemNew(n*sizeof(Int4));
      for(i= 0; i < n; i++) {
	  tilp->ids[i]= id= tax_getParent(DBctl, id);
      }
  }
  else {
      tilp->_num_ids= 1;
      tilp->ids= MemNew(sizeof(Int4));
      tilp->ids[0]= tax_getParent(DBctl, id);
  }
  return tilp;
}

#define CLASS_TXT_SCI 0
#define CLASS_TXT_COM 7
#define CLASS_TXT_SYN 3
#define CLASS_TXT_SYNCOM 8

/**************************************************************************
 *
 *	FsTaxGetRef
 *
 **************************************************************************/


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
  static char buff[256];

  buff[0]= '~';
  for(i= 0; syn[i].tax_id == id; i++) {
    if((syn[i].class_cde == SYNONYM) || (syn[i].class_cde == COMMON_NAME)) {
      list= ValNodeNew(list);
      if(syn[i].class_cde == SYNONYM) {
	StringCpy(&buff[1], syn[i].name_txt);
	list->data.ptrvalue= StringSave(buff);
      }
      else {
	list->data.ptrvalue= StringSave(syn[i].name_txt);
      }
      if(header == NULL) header= list;
    }
  }
  return header;
}
	   

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

    orgName= tax_findOrgName(DBctl, lineage[i]);
    this->data.ptrvalue= StringSave(orgName->name_txt);
    len+= StringLen(orgName->name_txt) + 3;
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

static OrgNamePtr bldOrgName(Int4 id, _taxNamePtr o_name)
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
  onp->div= StringSave(tax_getDivById(DBctl, node->division, TAX_DIV_CDE));
  onp->lineage= bldLineage(id);
	   
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

static void bldOrgRef(Int4 id, OrgRefPtr orp)
{
  _taxNamePtr o_name;
  int i;

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
  orp->orgname= bldOrgName(id, o_name);
}
 

static
ValNodePtr InsertSynonymList( ValNodePtr vnp, CharPtr s )
{

   if ( vnp == NULL ) {
	vnp = ValNodeNew(NULL);
   }
   else {
	vnp->next = ValNodeNew(NULL);
	vnp = vnp->next;
   }

   vnp->data.ptrvalue = StringSave( s );
   return( vnp );
}

OrgRefPtr FsTaxGetRef( int id_tax )
{
  OrgRefPtr	orp = NULL;
  ValNodePtr	vnp=NULL;
  CharPtr	p, s;
  DbtagPtr dbtag=NULL;
  ObjectIdPtr object_id;
  _taxNamePtr orgName;
  Int4 id;
  int new_stile= 0;

  if(id_tax < 0) {
    id_tax= -id_tax;
    new_stile= 1;
  }
  if(((id= getIdById(id_tax)) == 0) || ((orgName= tax_findOrgName(DBctl, id)) == NULL)) {
#ifdef WARN_ON_FAILURE
    ErrPostEx(SEV_WARNING,1,2,"FsTaxGetRef: Not found: <%d>", id_tax);
#endif
    return NULL;
  }
  
  orp = OrgRefNew();
  if(new_stile) {
    bldOrgRef(id, orp);
    return orp;
  }
#if 1
  orp->mod = NULL;
  orp->db = ValNodeNew(NULL);
  orp -> db -> data.ptrvalue = dbtag = DbtagNew();
  dbtag -> db = StringSave("taxon");
  dbtag -> tag = object_id = ObjectIdNew();
  object_id -> str = NULL;
  object_id -> id = id;
  orp->taxname = NULL;
  orp->common = NULL;
  orp->syn = NULL;
  
  while(orgName->tax_id == id) {
    switch(orgName->class_cde) {
    case CLASS_TXT_SCI :
      orp->taxname = StringSave(orgName->name_txt);
      break;
      
    case CLASS_TXT_COM:
      orp->common = StringSave(orgName->name_txt);
      break;

    default:
      break;
    }
    orgName++;
  }
#endif

  return orp;
}

/**************************************************************************
 *
 *	FsTaxGetLineage
 *
 **************************************************************************/

static
CharPtr s_FsTaxGetLineage(int id_tax)
{
  ValNodePtr head=NULL, this=NULL;
  int len = 2;
  CharPtr temp, retval = NULL;
  static Int4 lineage[128];
  int n, i;
  Int4 id;
  _taxNamePtr orgName;
  
  if((id= getIdById(id_tax)) == 0) {
#ifdef WARN_ON_FAILURE
    ErrPostEx(SEV_WARNING,1,2,"s_FsTaxGetLineage: Not found: <%d>", id_tax);
#endif
    return NULL;
  }

  n= tax_getLin(DBctl, id, lineage);

  for(i= 1; i < n; i++) {
    if((DBctl->tree[lineage[i]].flags & GB_HIDDEN) != 0) continue;

    if((DBctl->tree[lineage[i]].rank >= IsSpeciesRank) && head) break;

    this= ValNodeNew(this);
    if (head == NULL){
      head = this;
    }

    orgName= tax_findOrgName(DBctl, lineage[i]);
    this->data.ptrvalue= StringSave(orgName->name_txt);
    len+= StringLen(orgName->name_txt) + 3;
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


/**************************************************************************
 *
 *	FsTaxGetTaxonline
 *
 **************************************************************************/

ValNodePtr FsTaxGetTaxonline( int id_tax )
{
   ValNodePtr vnp = NULL;
   char s[255];
   CharPtr lineage = s_FsTaxGetLineage(id_tax);

   if (lineage != NULL) {
     if ( * lineage != '\0' ) {
       vnp = ValNodeNew(NULL);
       vnp->data.ptrvalue = lineage;
     }else{
       MemFree(lineage);
     }
   }

   return( vnp );
}


/**************************************************************************
 *
 *	FsTaxGetWithinDiv
 *
 **************************************************************************/

ValNodePtr FsTaxGetWithinDiv( int id_tax )
{
  ValNodePtr vnp = NULL;
  Int4 id;
  _taxNodePtr node;

  if((id= getIdById(id_tax)) == 0) {
#ifdef WARN_ON_FAILURE
    ErrPostEx(SEV_WARNING,1,2,"FsTaxGetWithinDiv: Not found: <%d>", id_tax);
#endif
    return NULL;
  }
  
  if((node= tax_getNode(DBctl, id)) != NULL) {
    vnp= ValNodeNew(NULL);
    vnp->data.ptrvalue= StringSave(tax_getDivById(DBctl, node->division, TAX_DIV_CDE));
  }
  return vnp;
}

/**************************************************************************
 *
 *	FsTaxGetGeneticCode
 *
 **************************************************************************/

GeneticCodeListPtr FsTaxGetGeneticCode( int id_tax )
{
  GeneticCodeListPtr gclp = NULL;
  Int4 id;
  _taxNodePtr node;

  if((id= getIdById(id_tax)) == 0) {
#ifdef WARN_ON_FAILURE
    ErrPostEx(SEV_WARNING,1,2,"FsTaxGetGeneticCode: Not found: <%d>", id_tax);
#endif
    return NULL;
  }
  
  if((node= tax_getNode(DBctl, id)) != NULL) {
    gclp= GeneticCodeListNew();
    gclp->genomic= node->gc;
    gclp->mitochondrial = node->mgc;
  }

  return gclp;

}

/**************************************************************************
 *
 *	SybaseTaxGetComplete
 *

As Per Jim Ostell, there is a paradigm shift here.
If a single tax_id can be found for the name, when a
name is given, then that tax id is used for the lookup,
else, if zero or multiple tax ids found, it fails.

-Karl  12/8/94
 **************************************************************************/
TaxCompleteListPtr FsTaxGetComplete( TaxonIdNamePtr tinp )
{
   TaxCompletePtr tcp;
   TaxCompleteListPtr tclp = NULL;
   Char name[255];
   TaxonNamePtr tnp = ValNodeNew(NULL)	;
   TaxonIdListPtr tilp=NULL;
   OrgRefPtr orp=NULL;
   Int4 tax_id = 0;
  _taxNodePtr node;


   if (tinp->choice == TaxonIdName_id) {
     tax_id= getIdById(tinp->data.intvalue);
   }
   else if(tinp->choice == TaxonIdName_name) {

     CopyString( name, (CharPtr)tinp->data.ptrvalue );
     tnp -> data.ptrvalue = name;
     tilp= FsTaxGetTaxId( tnp );

     if (tilp)
       if (tilp->_num_ids == 1) {
	 tax_id= tilp->ids[0];
       }else{
	 TaxonIdListFree(tilp);
	 return tclp;
       }
     TaxonIdListFree(tilp);
   }else{
     return tclp;
   }

   if((node= tax_getNode(DBctl, tax_id)) != NULL) {
     
     if ( tclp == NULL ) {
       tclp = TaxCompleteListNew();
       tclp->num = 0;
       tclp->info = TaxCompleteNew();
       tcp = tclp->info;
     } 
     tclp->num++;
     tcp->next = NULL;

     tcp->id_gc = node->gc;
     tcp->id_mgc = node->mgc;
     tcp->is_species_level= (DBctl->tree[tax_id].rank >= IsSpeciesRank)? 1 : 0;

     tcp->gb_div= StringSave(tax_getDivById(DBctl, node->division, TAX_DIV_CDE));
     tcp->embl_code= StringSave(node->embl_cde);
     tcp->name_gc= StringSave(tax_getGCById(DBctl, node->gc, TAX_GC_NM));
     tcp->name_mgc= StringSave(tax_getGCById(DBctl, node->mgc, TAX_GC_NM));

     orp=  FsTaxGetRef( tax_id );
     if (orp){
       tcp->sciname   = StringSave( orp -> taxname?orp->taxname:"" );
       tcp->comname   = StringSave( orp -> common?orp->common :"");
       OrgRefFree(orp);
     }else{
     }
     tcp->synonyms  = StringSave( "Call Jim Ostell if you need these synonyms" );

     tcp->lineage   = s_FsTaxGetLineage(tax_id);
   }
   else{
#ifdef WARN_ON_FAILURE
     ErrPostEx(SEV_WARNING,1,2,"FsTaxGetComplete: Not found: <%d>", tax_id);
#endif
   }
   return( tclp );
}


