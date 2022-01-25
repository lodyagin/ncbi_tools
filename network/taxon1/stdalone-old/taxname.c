/****************************
 * File: taxname.c
 * Description: taxonomy name/search functions
 */

#include <stdlib.h>
#include <rex_util.h>
#include "tax_cmmn.h"

static Int4 nof_des= 0;

typedef struct t_NameDesign{
  Int4 tax_id;
  CharPtr name_txt;
} _NameDesign;

static _NameDesign* design= NULL;

static int get_hashValue(CharPtr str)
{
  int i;
  int v;

  if(str == NULL) return 0;

  for(v= 0,i= 0; *str != '\0'; str++, i++) {
    if(isalnum(*str)) {
      v+= toupper(*str) - 0x20;
      v&= 0xFFF;
    }
  }
  return v;
}

static void insert_name(_hashItemPtr hTbl, _taxNamePtr name)
{
  _hashItemPtr item;
  
  item= &hTbl[get_hashValue(name->name_txt)];
  if(item->tax_name == NULL) {
    /* this hash value still free */
    item->next= NULL;
    item->tax_name= name;
  }
  else {
    /* already occupied */
    _hashItemPtr new_item= MemNew(sizeof(_hashItem));
    
    if(new_item != NULL) {
      new_item->next= item->next;
      item->next= new_item;
      new_item->tax_name= name;
    }
  }
}

/* Compare strings ignore non alphanum symbols and case */
static int strWCmp(CharPtr s1, CharPtr s2)
{
  int i;

  while((*s1 != '\0') && (*s2 != '\0')) {
    if(!isalnum(*s1)) {
      ++s1;
      continue;
    }
    if(!isalnum(*s2)) {
      ++s2;
      continue;
    }
    
    i= toupper(*s1) - toupper(*s2);
    if(i != 0) return i;
    
    ++s1; ++s2;
  }

  while(*s1 != '\0') {
    if(isalnum(*s1)) return 1;
    ++s1;
  }

  while(*s2 != '\0') {
    if(isalnum(*s2)) return -1;
    ++s2;
  }

  return 0;
}


static int find_exactName(_hashItemPtr hashTbl, CharPtr sname, _taxNamePtr** nameList)
{
  _hashItemPtr item;
  _taxNamePtr* list= NULL;
  _hashItemPtr tmp;
  int n= 0;

  item= &hashTbl[get_hashValue(sname)];

  if(item->tax_name == NULL) return 0; /* no name found */

  /* on the first pass we'll search for exact hit */

  for(tmp= item; tmp != NULL; tmp= tmp->next) {
    if(strcmp(tmp->tax_name->name_txt, sname) == 0) {
      /* we've found it! */
      if(list == NULL) {
	list= MemNew(sizeof(_taxNamePtr));
	if(list != NULL) {
	  *list= tmp->tax_name;
	  n= 1;
	}
	else {
	  tax_outMsg(TAX_ERROR, TAX_ERR_NO_MEMORY, "find_exactName: not enaugh memory");
	  return -1;
	}
      }
      else {
	++n;
	list= MemMore(list, n*sizeof(_taxNamePtr));
	if(list != NULL) {
	  list[n-1]= tmp->tax_name;
	}
	else {
	  tax_outMsg(TAX_ERROR, TAX_ERR_NO_MEMORY, "find_exactName: not enaugh memory");
	  return -1;
	}
      }
    }
  }

  if(n > 0) {
    /* we've found something */
    *nameList= list;
    return n;
  }

  /* First pass was unsuccessfull, let's try the case ignore search */

  for(tmp= item; tmp != NULL; tmp= tmp->next) {
    if(StringICmp(tmp->tax_name->name_txt, sname) == 0) {
      /* we've found it! */
      if(list == NULL) {
	list= MemNew(sizeof(_taxNamePtr));
	if(list != NULL) {
	  *list= tmp->tax_name;
	  n= 1;
	}
	else {
	  tax_outMsg(TAX_ERROR, TAX_ERR_NO_MEMORY, "find_exactName: not enaugh memory");
	  return -1;
	}
      }
      else {
	++n;
	list= MemMore(list, n*sizeof(_taxNamePtr));
	if(list != NULL) {
	  list[n-1]= tmp->tax_name;
	}
	else {
	  tax_outMsg(TAX_ERROR, TAX_ERR_NO_MEMORY, "find_exactName: not enaugh memory");
	  return -1;
	}
      }
    }
  }

  if(n > 0) {
    /* we've found something */
    *nameList= list;
  }
  
  return n;
}


static int find_alphanumName(_hashItemPtr hashTbl, CharPtr sname, _taxNamePtr** nameList)
{
  _hashItemPtr item;
  _taxNamePtr* list= NULL;
  _hashItemPtr tmp;
  int n= 0;

  item= &hashTbl[get_hashValue(sname)];

  if(item->tax_name == NULL) return 0; /* no name found */

  /* on the first pass we'll search for exact hit */

  for(tmp= item; tmp != NULL; tmp= tmp->next) {
    if(strcmp(tmp->tax_name->name_txt, sname) == 0) {
      /* we've found it! */
      if(list == NULL) {
	list= MemNew(sizeof(_taxNamePtr));
	if(list != NULL) {
	  *list= tmp->tax_name;
	  n= 1;
	}
	else {
	  tax_outMsg(TAX_ERROR, TAX_ERR_NO_MEMORY, "find_exactName: not enaugh memory");
	  return -1;
	}
      }
      else {
	++n;
	list= MemMore(list, n*sizeof(_taxNamePtr));
	if(list != NULL) {
	  list[n-1]= tmp->tax_name;
	}
	else {
	  tax_outMsg(TAX_ERROR, TAX_ERR_NO_MEMORY, "find_exactName: not enaugh memory");
	  return -1;
	}
      }
    }
  }

  if(n > 0) {
    /* we've found something */
    *nameList= list;
    return n;
  }

  /* First pass was unsuccessfull, let's try the case ignore search */

  for(tmp= item; tmp != NULL; tmp= tmp->next) {
    if(StringICmp(tmp->tax_name->name_txt, sname) == 0) {
      /* we've found it! */
      if(list == NULL) {
	list= MemNew(sizeof(_taxNamePtr));
	if(list != NULL) {
	  *list= tmp->tax_name;
	  n= 1;
	}
	else {
	  tax_outMsg(TAX_ERROR, TAX_ERR_NO_MEMORY, "find_exactName: not enaugh memory");
	  return -1;
	}
      }
      else {
	++n;
	list= MemMore(list, n*sizeof(_taxNamePtr));
	if(list != NULL) {
	  list[n-1]= tmp->tax_name;
	}
	else {
	  tax_outMsg(TAX_ERROR, TAX_ERR_NO_MEMORY, "find_exactName: not enaugh memory");
	  return -1;
	}
      }
    }
  }

  if(n > 0) {
    /* we've found something */
    *nameList= list;
    return n;
  }

  /* Nothing found on second pass too, let's try to exclude all non alphanum symbols 
     and go through again 
  */

  for(tmp= item; tmp != NULL; tmp= tmp->next) {
    if(strWCmp(tmp->tax_name->name_txt, sname) == 0) {
      /* we've found it! */
      if(list == NULL) {
	list= MemNew(sizeof(_taxNamePtr));
	if(list != NULL) {
	  *list= tmp->tax_name;
	  n= 1;
	}
	else {
	  tax_outMsg(TAX_ERROR, TAX_ERR_NO_MEMORY, "find_exactName: not enaugh memory");
	  return -1;
	}
      }
      else {
	++n;
	list= MemMore(list, n*sizeof(_taxNamePtr));
	if(list != NULL) {
	  list[n-1]= tmp->tax_name;
	}
	else {
	  tax_outMsg(TAX_ERROR, TAX_ERR_NO_MEMORY, "find_exactName: not enaugh memory");
	  return -1;
	}
      }
    }
  }

  if(n > 0) {
    /* we've found something */
    *nameList= list;
  }
  
  return n;
}

static void freeHashTbl(_taxDBCtlPtr ctl)
{
  _hashItemPtr item;
  int i;
  _hashItemPtr nxt;

  if(ctl->hashTbl == NULL) return;

  /* free all chain in hash structure */
  for(i= 0; i <= 0xFFF; i++) {
    for(item= ctl->hashTbl[i].next; item != NULL; item= nxt) {
      nxt= item->next;
      MemFree(item);
    }
  }


  /* free table itself */
  MemFree(ctl->hashTbl);
  ctl->hashTbl= NULL;
}


static int iniHashTbl(_taxDBCtlPtr ctl)
{
  _hashItemPtr tbl;
  int i;

  if(ctl->hashTbl != NULL) freeHashTbl(ctl);

  if((tbl= MemNew(sizeof(_hashItem) * (0xFFF + 1))) == NULL) {
    tax_outMsg(TAX_ERROR, TAX_ERR_NO_MEMORY, "iniHashTbl: not enaugh memory");
    return -1;
  }

  for(i= 0; i <= 0xFFF; i++) {
    tbl[i].tax_name= NULL;
    tbl[i].next= NULL;
  }

  ctl->hashTbl= tbl;
  return 0;
}


/*------------------------------------------------
 * load_bin - load tree from binary file
 */
static int load_bin(_taxDBCtlPtr ctl)
{
  FILE* f;
  _taxNamePtr nTbl;
  Int4 i;
  Int4 nof_names= tax_getNofNames(ctl);

  if(nof_names < 1) {
    tax_outMsg(TAX_WARN, TAX_WARN_NO_NAMES, "load_bin: no names in taxonomy");
    return 1;
  }
	
  if((f= fopen(tax_getFileName(ctl, TAX_NAME_FILE), "rb")) == NULL) {
    tax_outMsg(TAX_ERROR, TAX_ERR_NO_FILE, "load_bin: Can't open taxonomy names file");
    return -1;
  }

  if((nTbl= malloc(sizeof(_taxName) * nof_names)) == NULL) {
    tax_outMsg(TAX_ERROR, TAX_ERR_NO_MEMORY, "load_bin: Not enaugh memory for taxonomy names");
    return -1;
  }

  /* read file and fill-up array */
  for(i= 0; i < nof_names; i++) {
    nTbl[i].tax_id= fsb_getInt4(f);
    nTbl[i].class_cde= fsb_getInt1(f);
    nTbl[i].name_txt= fsb_getString(f);
    nTbl[i].unique_name= fsb_getString(f);
    insert_name(ctl->hashTbl, &nTbl[i]);   /* insert name into hash table */
  }

  ctl->names= nTbl;

  nof_des= fsb_getInt4(f);
  if(nof_des > 0) {
    design= MemNew(nof_des*sizeof(_NameDesign));
    for(i= 0; i < nof_des; i++) {
      design[i].tax_id= fsb_getInt4(f);
      design[i].name_txt= fsb_getString(f);
    }
  }
  fclose(f);
  return 0;
}

/*---------------------------------------------------
 * load_txt - load taxonomy tree from text file
 */
static int load_txt(_taxDBCtlPtr ctl)
{
  FILE* f;
  Int4 i;
  _taxNamePtr nTbl;
  Int4 nof_names= tax_getNofNames(ctl);


  if(nof_names < 1) {
    tax_outMsg(TAX_WARN, TAX_WARN_NO_NAMES, "load_txt: no taxonomy names");
    return 1;
  }

  if((f= fopen(tax_getFileName(ctl, TAX_NAME_FILE), "r")) == NULL) {
    tax_outMsg(TAX_ERROR, TAX_ERR_NO_FILE, "load_txt: Can't open taxonomy names file");
    return -1;
  }

  if((nTbl= malloc(sizeof(_taxName) * nof_names)) == NULL) {
    tax_outMsg(TAX_ERROR, TAX_ERR_NO_MEMORY, "load_txt: Not enaugh memory for taxonomy names");
    return -1;
  }

  /* read file and fill-up array */
  for(i= 0; i < nof_names; i++) {
    nTbl[i].tax_id=   fs_getInt4(f, tax_getMarker(ctl));
    nTbl[i].class_cde= fs_getInt2(f, tax_getMarker(ctl));
    nTbl[i].name_txt=  fs_getString(f, tax_getMarker(ctl));
    nTbl[i].unique_name=  fs_getString(f, tax_getMarker(ctl));
    insert_name(ctl->hashTbl, &nTbl[i]);   /* insert name into hash table */
  }

  ctl->names= nTbl;
  fclose(f);
  return 0;
}

static int load_asn(_taxDBCtlPtr ctl)
{
    return -10; /* not done yet */
}

/*==============================================================
 * tax_freeNames - free memory
 */
void tax_freeNames(_taxDBCtlPtr ctl)
{
  _taxNamePtr nTbl= ctl->names;
  Int4 nof_names= tax_getNofNames(ctl);
  Int4 i;
  
  if((nof_names > 0) && (nTbl != NULL)) {
    for(i= 0; i != nof_names; i++) {
      if(nTbl[i].name_txt != NULL) MemFree(nTbl[i].name_txt);
      if(nTbl[i].unique_name != NULL) MemFree(nTbl[i].unique_name);
    }
    MemFree(nTbl);
  }

  freeHashTbl(ctl);
  ctl->names= NULL;
  if(design != NULL) {
    for(i= 0; i < nof_des; i++) {
      if(design[i].name_txt != NULL) MemFree(design[i].name_txt);
    }
    MemFree(design);
    design= NULL;
    nof_des= 0;
  }
}

int tax_loadNames(_taxDBCtlPtr ctl)
{

  if(ctl->names != NULL) tax_freeNames(ctl);

  if(iniHashTbl(ctl)) return -1;

  switch(tax_getDbType(ctl)) {
  case TAX_BIN: return load_bin(ctl);
  case TAX_TXT: return load_txt(ctl);
  case TAX_ASN: return load_asn(ctl);
  default: break;
  }

  tax_outMsg(TAX_ERROR, TAX_ERR_FILE_TYPE, "tax_loadNames: Wrong file type for taxonomy names");
  return -1;
}

static int rex_search(_taxDBCtlPtr ctl, CharPtr sname, _taxNamePtr** nameList)
{
  char nBuff[256];
  Int4 nof_names= tax_getNofNames(ctl);
  _taxNamePtr nTbl= ctl->names;
  _taxNamePtr* list= NULL;
  Int4 i, k, n= 0;
  rex_handler rh;
  CharPtr tmp;

  strncpy(&nBuff[1], sname, 250);
  nBuff[0]= '@';
  k= strlen(nBuff);
  nBuff[k]= '@';
  nBuff[k+1]= '\0';
  rh= rex_setExpr(nBuff, REX_NO_CASE | REX_NO_SPACE);
  if(rh == NULL) return 0;

  nBuff[0]= '@';
  for(i= 0; i < nof_names; i++) {
    tmp= nTbl[i].name_txt;
    for(k= 0; (k < 250) && (tmp[k] != '\0'); k++) {
      nBuff[k+1]= toupper(tmp[k]);
    }
    k++;
    nBuff[k]= '@';
    nBuff[k+1]= '\0';

    if(rex_cmp(rh, nBuff)) {
      /* name found */
      if(list == NULL) {
	if((list= MemNew(sizeof(_taxNamePtr))) == NULL) {
	  tax_outMsg(TAX_ERROR, TAX_ERR_NO_MEMORY, "rex_search: not enaugh memory for name list");
	  return 0;
	}
	*list= &nTbl[i];
	n= 1;
      }
      else {
	++n;
	if((list= MemMore(list, sizeof(_taxNamePtr)*n)) == NULL) {
	  tax_outMsg(TAX_ERROR, TAX_ERR_NO_MEMORY, "rex_search: not enaugh memory for name list");
	  return 0;
	}
	list[n-1]= &nTbl[i];
      }
    }
  }

  if(n > 0) *nameList= list;
  free(rh);
  return n;
}

static CharPtr get_token(CharPtr str, CharPtr token)
{
  int i;

  token[2]= '\0';

  while(IS_WHITESP(*str)) {
    if(*str == '\0') return NULL;
    ++str;
  }

  if(*str == '\0') return NULL;
  
  for(i= 1; i < 250; i++) {
    if(IS_WHITESP(*str)) {
      token[i]= ' ';
      token[i+1]= '\0';
      return str;
    }

    if(*str == '\0') {
      token[i]= ' ';
      token[i+1]= '\0';
      return NULL;
    }

    token[i]= TO_UPPER(*str);
    str++;
  }

  token[i]= ' ';
  token[i+1]= '*';
  token[i+2]= '\0';
  return str;
}
  
static int token_search(_taxDBCtlPtr ctl, CharPtr sname, _taxNamePtr** nameList)
{
  char nBuff[256];
  Int4 nof_names= tax_getNofNames(ctl);
  _taxNamePtr nTbl= ctl->names;
  _taxNamePtr* list= NULL;
  Int4 i, k, n= 0;
  rex_handler rh[16];
  Int2 nof_rh, j, res;
  CharPtr tail= sname;

  nBuff[0]= '*';
  nBuff[1]= ' ';
  for(nof_rh= 0; (nof_rh != 16) && (tail != NULL); nof_rh++) {
    tail= get_token(tail, nBuff);
    rh[nof_rh]= rex_setExpr(nBuff, REX_NO_CASE | REX_NO_SPACE);
    if(rh[nof_rh] == NULL) nof_rh--;
  }
  if(nof_rh < 1) return 0;

  nBuff[0]= ' ';
  for(i= 0; i < nof_names; i++) {
    tail= nTbl[i].name_txt;
    for(k= 0; (k < 250) && (tail[k] != '\0'); k++) {
      nBuff[k+1]= toupper(tail[k]);
    }
    k++;
    nBuff[k]= ' ';
    nBuff[k+1]= '\0';

    for(res= 1, j= 0; (j < nof_rh) && (res > 0); j++) {
      res= rex_cmp(rh[j], nBuff);
    }

    if(res) {
      /* name found */
      if(list == NULL) {
	if((list= MemNew(sizeof(_taxNamePtr))) == NULL) {
	  tax_outMsg(TAX_ERROR, TAX_ERR_NO_MEMORY, "rex_search: not enaugh memory for name list");
	  return 0;
	}
	*list= &nTbl[i];
	n= 1;
      }
      else {
	++n;
	if((list= MemMore(list, sizeof(_taxNamePtr)*n)) == NULL) {
	  tax_outMsg(TAX_ERROR, TAX_ERR_NO_MEMORY, "rex_search: not enaugh memory for name list");
	  return 0;
	}
	list[n-1]= &nTbl[i];
      }
    }
  }

  if(n > 0) *nameList= list;
  for(j= 0; j < nof_rh; j++) {
    free(rh[j]);
  }
  return n;
}

int tax_findByName(_taxDBCtlPtr ctl, CharPtr inp_sname, int mode, _taxNamePtr** nameList)
{
    char sname[128];
    int i, k;

    for(i= 0; inp_sname[i] != '\0'; i++) {
	if(!isspace(inp_sname[i])) break;
    }

    sname[0]= inp_sname[i++];

    for(k= 0; k < 120; i++) {
	if(inp_sname[i] == '\0') {
	    break;
	}
	if(isspace(inp_sname[i])) {
	    if(isspace(sname[k])) continue;
	    else {
		sname[++k]= ' ';
	    }
	}
	else {
	    sname[++k]= inp_sname[i];
	}
    }

    if(sname[k] != ' ') k++;

    sname[k]= '\0';

    switch(mode) {
    case TAX_RE_SEARCH:    return rex_search(ctl, sname, nameList);
    case TAX_TOKEN_SEARCH: return token_search(ctl, sname, nameList);
    case TAX_ALPHANUM_SEARCH: return find_alphanumName(ctl->hashTbl, sname, nameList);
    default:               return find_exactName(ctl->hashTbl, sname, nameList);

    }
}

_taxNamePtr tax_findOrgName(_taxDBCtlPtr ctl, Int4 id)
{
  _taxNamePtr nameList;
  Int4 m, l= 0, r= ctl->NofNames - 1;
  

  nameList= ctl->names;
  if(nameList == NULL) return NULL;

  /* binary search */
  while((r - l) > 1) {
    m= (r + l)/2;
    if(nameList[m].tax_id >= id) r= m;
    else l= m;
  }

  if((nameList[r].tax_id < id) || (nameList[l].tax_id  > id)) return NULL;

  if(nameList[l].tax_id == id) {
    while(l > 0) {
      if(nameList[l-1].tax_id == id) l--;
      else break;
    }
    return &nameList[l];
  }

  if(nameList[r].tax_id == id) return &nameList[r];
  
  return NULL;
}
  
Int4 tax_getDesignator(CharPtr name)
{
  Int4 i;

  for(i= 0; i < nof_des; i++) {
    if(strcmp(name, design[i].name_txt) == 0) {
      return design[i].tax_id;
    }
  }

  /* try to find ignore case */
  for(i= 0; i < nof_des; i++) {
    if(StringICmp(name, design[i].name_txt) == 0) {
      return design[i].tax_id;
    }
  }
  
  /* try to find ignore case and non alphanum */
  for(i= 0; i < nof_des; i++) {
    if(strWCmp(name, design[i].name_txt) == 0) {
      return design[i].tax_id;
    }
  }
  return 0;
}


