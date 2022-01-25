/****************************
 * File: tree.c
 * Description: taxonomy tree functions
 */

#include <stdlib.h>
#include "tax_cmmn.h"

static _taxNodePtr nodeSet= NULL;

/*------------------------------------------------
 * load_bin - load tree from binary file
 */
static int load_bin(_taxDBCtlPtr ctl)
{
  FILE* f;
  _treePtr node;
  long nof_ids= tax_getMaxId(ctl);

  if(nof_ids < 1) {
    tax_outMsg(TAX_WARN, TAX_WARN_EMPTY_TREE, "load_bin: taxonomy tree is empty");
    return 1;
  }
	
  if((f= fopen(tax_getFileName(ctl, TAX_TREE_FILE), "rb")) == NULL) {
    tax_outMsg(TAX_ERROR, TAX_ERR_NO_FILE, "load_bin: Can't open taxonomy tree file");
    return -1;
  }

  if((node= malloc(sizeof(_tree) * (nof_ids+1))) == NULL) {
    tax_outMsg(TAX_ERROR, TAX_ERR_NO_MEMORY, "load_bin: Not enaugh memory for taxonomy tree");
    return -1;
  }
  memset(node, 0, sizeof(_tree) * (nof_ids+1));

  /* read file and fill-up array */
  if(fread(node, sizeof(_tree), nof_ids, f) < nof_ids) {
    free(node);
    tax_outMsg(TAX_ERROR, TAX_ERR_WRONG_FILE, "load_bin: taxonomy tree currupted");
    return -1;
  }

  ctl->tree= node;
  fclose(f);
  return 0;
}

/*---------------------------------------------------
 * load_txt - load taxonomy tree from text file
 */
static int load_txt(_taxDBCtlPtr ctl)
{
  FILE* f;
  long i;
  _treePtr nmc;
  long nof_ids= tax_getMaxId(ctl);
  _treePtr node;


  if(nof_ids < 1) {
    tax_outMsg(TAX_WARN, TAX_WARN_EMPTY_TREE, "load_txt: taxonomy tree is empty");
    return 1;
  }

  if((f= fopen(tax_getFileName(ctl, TAX_TREE_FILE), "r")) == NULL) {
    tax_outMsg(TAX_ERROR, TAX_ERR_NO_FILE, "load_txt: Can't open taxonomy tree file");
    return -1;
  }

  if((node= malloc(sizeof(_tree) * nof_ids)) == NULL) {
    tax_outMsg(TAX_ERROR, TAX_ERR_NO_MEMORY, "load_txt: Not enaugh memory for taxonomy tree");
    return -1;
  }

  /* read file and fill-up array */
  for(i= 0; i < nof_ids; i++) {
    node[i].parent=   fs_getLong(f, tax_getMarker(ctl));
    node[i].child=    fs_getLong(f, tax_getMarker(ctl));
    node[i].sibling=  fs_getLong(f, tax_getMarker(ctl));
    node[i].dataKey=  fs_getLong(f, tax_getMarker(ctl));
    node[i].rank=     fs_getInt2(f, tax_getMarker(ctl));
    node[i].flags=    fs_getInt2(f, tax_getMarker(ctl));
  }

  ctl->tree= node;
  fclose(f);
  return 0;
}

static int load_asn(_taxDBCtlPtr ctl)
{
    return -10; /* not done yet */
}

/*==============================================================
 * tax_freeClasses - free memory
 */
void tax_freeTree(_taxDBCtlPtr ctl)
{
  Int4 i;

  if(nodeSet != NULL) {
    for(i= 0; i < tax_getNofNodes(ctl); i++) {
      if(nodeSet[i].comment != NULL) MemFree(nodeSet[i].comment);
    }
    MemFree(nodeSet);
    nodeSet= NULL;
  }

  if((tax_getMaxId(ctl) > 0) && (ctl->tree != NULL)) {
    free(ctl->tree);
  }
  ctl->tree= NULL;
}

int tax_loadTree(_taxDBCtlPtr ctl)
{

  if(ctl->tree != NULL) tax_freeTree(ctl);

  switch(tax_getDbType(ctl)) {
  case TAX_BIN: return load_bin(ctl);
  case TAX_TXT: return load_txt(ctl);
  case TAX_ASN: return load_asn(ctl);
  default: break;
  }

  tax_outMsg(TAX_ERROR, TAX_ERR_FILE_TYPE, "tax_loadTree: Wrong file type for taxonomy tree");
  return -1;
}

static CharPtr emblcode(CharPtr name)
{
  static char embl[4];
  int i;

  embl[0]= embl[1]= embl[2]= embl[3]= '\0';
  for(i= 0; name[i] != '\0'; i++) {
    if(IS_ALPHA(name[i])) {
      embl[0]= TO_UPPER(name[i]);
      break;
    }
  }

  if(name[i] != '\0') {
    for(;name[i] != '\0'; i++) {
      if(IS_WHITESP(name[i])) break;
    }
    while(name[i++] != '\0') {      
      if(IS_ALPHANUM(name[i])) {
	embl[1]= TO_UPPER(name[i]);
	break;
      }
    }
  }

  if(embl[0] == '\0') {
    embl[0]=embl[1]= 'X';
  }
  else if(embl[1] == '\0') {
    embl[1]= 'X';
  }

  return embl;
}


static int loadNodesBin(_taxDBCtlPtr ctl)
{
  _taxNodePtr node;
  Int4 nofNodes= tax_getNofNodes(ctl);
  Int4 i;

  if(ctl->nodefile == NULL) {
    /* open node file */
    if((ctl->nodefile= fopen(tax_getFileName(ctl, TAX_NODE_FILE), "rb")) == NULL) {
      tax_outMsg(TAX_ERROR, TAX_ERR_NO_FILE, "loadNodesBin: Can't open node file");
      return NULL;
    }
  }
  else {
    fseek(ctl->nodefile, 0L, SEEK_SET);
  }
    

  nodeSet= MemNew(sizeof(_taxNode)*nofNodes);
  if(nodeSet == NULL) {
    tax_outMsg(TAX_ERROR, TAX_ERR_NO_MEMORY, "loadNodesBin: Not enaugh memory");
    return -1;
  }

  for(i= 0; i < nofNodes; i++) {
    if(fread(&nodeSet[i], sizeof(_taxNode), 1, ctl->nodefile) != 1) {
      tax_outMsg(TAX_ERROR, TAX_ERR_WRONG_FILE, "loadNodesBin: File currapted");
      return -2;
    }
    if(nodeSet[i].comment_size > 0) nodeSet[i].comment= fsb_getString(ctl->nodefile);
    else nodeSet[i].comment= NULL;
  }
  fclose(ctl->nodefile);
  ctl->nodefile= NULL;
  return 0;
}


static _taxNodePtr getNodeMem(_taxDBCtlPtr ctl, Int4 id)
{
  Int4 l, u, m;
  _taxNodePtr node;
  _taxNamePtr orgName;

  l= 0;
  u= tax_getNofNodes(ctl);

  while((u-l) > 1) {
    m= (u + l)/2;
    if(nodeSet[m].tax_id > id) u= m;
    else l= m;
  }

  if(nodeSet[l].tax_id == id) m= l;
  else m= u;

  if((node= MemNew(sizeof(_taxNode))) == NULL) {
    tax_outMsg(TAX_ERROR, TAX_ERR_NO_MEMORY, "getNodeMem: Not enaugh memory for tax node");
    return NULL;
  }
  
  *node= nodeSet[m];
  if(node->comment_size > 0) node->comment= StringSave(nodeSet[m].comment);
  if((node->embl_cde[0] == '\0') && (ctl->tree[id].rank > 25)) {
    /* we need to calculate the default EMBL code */
    if(ctl->names == NULL) tax_loadNames(ctl);
    orgName= tax_findOrgName(ctl, id);
    strcpy(node->embl_cde, emblcode(orgName->name_txt));
  }
   return node;
}
  

static _taxNodePtr getNodeBin(_taxDBCtlPtr ctl, Int4 id)
{
  _taxNodePtr node;
  _taxNamePtr orgName;

  if(ctl->nodefile == NULL) {
    /* open node file */
    if((ctl->nodefile= fopen(tax_getFileName(ctl, TAX_NODE_FILE), "rb")) == NULL) {
      tax_outMsg(TAX_ERROR, TAX_ERR_NO_FILE, "getNodeBin: Can't open node file");
      return NULL;
    }
  }

  if((node= MemNew(sizeof(_taxNode))) == NULL) {
    tax_outMsg(TAX_ERROR, TAX_ERR_NO_MEMORY, "getNodeBin: Not enaugh memory for tax node");
    return NULL;
  }

  fseek(ctl->nodefile, ctl->tree[id].dataKey, SEEK_SET);
  fread(node, sizeof(_taxNode), 1, ctl->nodefile);
  if(node->comment_size > 0) node->comment= fsb_getString(ctl->nodefile);
  if((node->embl_cde[0] == '\0') && (ctl->tree[id].rank > 25)) {
    /* we need to calculate the default EMBL code */
    if(ctl->names == NULL) tax_loadNames(ctl);
    orgName= tax_findOrgName(ctl, id);
    strcpy(node->embl_cde, emblcode(orgName->name_txt));
  }
   return node;
}



static _taxNodePtr getNodeTxt(_taxDBCtlPtr ctl, Int4 id)
{
  _taxNodePtr node;
  FILE* f;
  CharPtr marker;
  CharPtr t;
  _taxNamePtr orgName;
  
  if(ctl->nodefile == NULL) {
    /* open node file */
    if((ctl->nodefile= fopen(tax_getFileName(ctl, TAX_NODE_FILE), "r")) == NULL) {
      tax_outMsg(TAX_ERROR, TAX_ERR_NO_FILE, "getNodeTxt: Can't open node file");
      return NULL;
    }
  }

  if((node= MemNew(sizeof(_taxNode))) == NULL) {
    tax_outMsg(TAX_ERROR, TAX_ERR_NO_MEMORY, "getNodeTxt: Not enaugh memory for tax node");
    return NULL;
  }

  f= ctl->nodefile;
  marker= tax_getMarker(ctl);

  fseek(f, ctl->tree[id].dataKey, SEEK_SET);
  node->tax_id= fs_getInt4(f, marker);
  t= fs_getString(f, marker);
  strncpy(node->embl_cde, t, 4);
  MemFree(t);
  if((node->embl_cde[0] == '\0') && (ctl->tree[id].rank > 25)) {
    /* we need to calculate the default EMBL code */
    if(ctl->names == NULL) tax_loadNames(ctl);
    orgName= tax_findOrgName(ctl, id);
    strcpy(node->embl_cde, emblcode(orgName->name_txt));
  }
    
  node->division= fs_getInt2(f, marker);
  node->gc= fs_getInt2(f, marker);
  node->mgc= fs_getInt2(f, marker);
  node->comment_size= fs_getInt2(f, marker);
  if(node->comment_size > 0) node->comment= fs_getString(f, marker);
  return node;
}
  
static _taxNodePtr getNodeAsn(_taxDBCtlPtr ctl, Int4 id)
{
  return NULL;
}

void tax_freeNode(_taxNodePtr node)
{
  if(node != NULL) {
    if(node->comment_size > 0) MemFree(node->comment);
    MemFree(node);
  }
}
     

_taxNodePtr tax_getNode(_taxDBCtlPtr ctl, Int4 id)
{

  if((id < 0) && (nodeSet == NULL)) {
    loadNodesBin(ctl);
    return NULL;
  }

  if((ctl->tree == NULL) || (id < 0) || (id > tax_getMaxId(ctl))) return NULL;

  if(ctl->tree[id].parent < 1) {
    /* deleted node */
    return NULL;
  }

  while(ctl->tree[id].child == -1) {
    id=ctl->tree[id].parent;
  }

  if(nodeSet != NULL) return getNodeMem(ctl, id);

  switch(tax_getDbType(ctl)) {
  case TAX_BIN: return getNodeBin(ctl, id);
  case TAX_TXT: return getNodeTxt(ctl, id);
  case TAX_ASN: return getNodeAsn(ctl, id);
  default: break;
  }

  tax_outMsg(TAX_ERROR, TAX_ERR_FILE_TYPE, "tax_getNode: Wrong file type for taxonomy node");
  return NULL;
}
  
Int2 tax_getLin(_taxDBCtlPtr ctl, Int4 id, Int4Ptr lin)
{
  int n, k;
  Int4 i;

  if((ctl->tree == NULL) || (id < 0) || (id > tax_getMaxId(ctl))) return 0;

  if(ctl->tree[id].parent < 1) {
    /* deleted node */
    return 0;
  }

  while(ctl->tree[id].child == -1) {
    id=ctl->tree[id].parent;
  }

  i= ctl->tree[id].parent;
  if(i == id) return 0; /* No lineage */
  
  /* calculate the number of nodes in lineage */
  
  for(k= 1; i != ctl->tree[i].parent; i= ctl->tree[i].parent) {
    k++;
  }
  
  n= k;

  for(i= ctl->tree[id].parent; k-- > 0; i= ctl->tree[i].parent) {
    lin[k]= i;
  }

  return n;
}
