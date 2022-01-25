/**************************
 * File: tax_cmmn.h
 * Description: common types, constants and variables
 */

#ifndef TAX_CMMN_H_DONE
#define TAX_CMMN_H_DONE

#include <ncbi.h>
#include "taxobj.h"
#include "fsutils.h"

#define NOF_FILES 7

/*----------------------------------
 * File types
 */

#define TAX_BIN 'b'
#define TAX_TXT 't'
#define TAX_ASN '1'

/*----------------------------------
 * message types
 */
#define TAX_MSG 0
#define TAX_WARN 1
#define TAX_ERROR -1

/*----------------------------------
 * message codes
 */
/* warnings */
#define TAX_WARN_NO_RANKS       100
#define TAX_WARN_NO_CLASSES     101
#define TAX_WARN_NO_DIVS        102
#define TAX_WARN_NO_GC          103
#define TAX_WARN_EMPTY_TREE     104
#define TAX_WARN_NO_NAMES       105

/* errors */
#define TAX_ERR_NO_MEMORY      1000
#define TAX_ERR_NO_FILE        1001
#define TAX_ERR_FILE_TYPE      1002
#define TAX_ERR_WRONG_FILE     1003

/*----------------------------------
 * file indexes
 */
#define TAX_TREE_FILE  0
#define TAX_NAME_FILE  1
#define TAX_NODE_FILE  2
#define TAX_GC_FILE    3
#define TAX_RANK_FILE  4
#define TAX_CLASS_FILE 5
#define TAX_DIVS_FILE  6
/*----------------------------------
 * various constants
 */
#define CTL_FILE_MARKER "\n"

#define TAX_DIV_TXT 0
#define TAX_DIV_CDE 1
#define TAX_DIV_COM 2

#define TAX_GC_NM     0
#define TAX_GC_ABBREV 1
#define TAX_GC_CDE    2
#define TAX_GC_STARTS 3

#define TAX_NAME_SEARCH  0
#define TAX_RE_SEARCH    1
#define TAX_TOKEN_SEARCH 2
#define TAX_ALPHANUM_SEARCH 3

/*---------------------------------
 * private types
 */

typedef struct t_rank {
  Int2 rank_id;
  CharPtr rank_txt;
} _rank, PNTR _rankPtr;

typedef struct t_nmclass {
  Int2 class_cde;
  CharPtr class_txt;
} _nmclass, PNTR _nmclassPtr;

typedef struct t_division {
  Int2 div_id;
  char div_cde[4];
  CharPtr div_txt;
  CharPtr div_comments_txt;
} _division, PNTR _divisionPtr;

typedef struct t_gencode {
  Int2 gc_id;
  char gc_abbrev[6];
  CharPtr gc_nm;
  CharPtr gc_cde;
  CharPtr gc_starts;
} _gencode, PNTR _gencodePtr;

typedef struct t_hashItem {
  _taxNamePtr tax_name;
  struct t_hashItem* next;
} _hashItem, PNTR _hashItemPtr;
  
#if 0
typedef struct t_taxName {
  int tax_id;
  _TINY class_cde;
  CharPtr name_txt;
} _taxName, PNTR _taxNamePtr;
#endif

typedef _taxCtlNodePtr _treePtr;
typedef _taxCtlNode _tree;

/*=================================
 * public types
 */

typedef struct t_taxDBCtl {
  _taxCtlNodePtr tree;
  _hashItemPtr hashTbl;
  _taxNamePtr names;
  _rankPtr rank;
  _nmclassPtr nmclass;
  _divisionPtr division;
  _gencodePtr gc;
  FILE* nodefile;
  Int4 MaxId;
  Int4 NofNodes;
  Int4 NofNames;
  Int2 NofRanks;
  Int2 NofDivs;
  Int2 NofClasses;
  Int2 NofGCs;
  CharPtr path;
  CharPtr FileName[NOF_FILES];
  char DbType;
  char DbFieldMarker[7];
} _taxDBCtl;

typedef _taxDBCtl PNTR _taxDBCtlPtr;

/*===================================
 * prototypes and macros
 */

typedef void (*tax_MessageHandler)(int, int, CharPtr);

#define tax_getNodeRank(c, i) ((c)->tree[i].rank)
#define tax_nodeGBHide(c, i) ((c)->tree[i].flags & GB_HIDDEN)

_taxDBCtlPtr tax_loadDBCtl(CharPtr DBPath, CharPtr ctlFileName);
_taxDBCtlPtr tax_loadTaxDB(CharPtr DBPath, CharPtr ctlFileName);
Int4 tax_getMaxId(_taxDBCtlPtr ctl);
Int4 tax_getNofNodes(_taxDBCtlPtr ctl);
Int4 tax_getNofNames(_taxDBCtlPtr ctl);
Int2 tax_getNofRanks(_taxDBCtlPtr ctl);
Int2 tax_getNofDivs(_taxDBCtlPtr ctl);
Int2 tax_getNofClasses(_taxDBCtlPtr ctl);
Int2 tax_getNofGCs(_taxDBCtlPtr ctl);
void tax_closeDbCtl(_taxDBCtlPtr ctl);
void tax_outMsg(int msg_type, int msg_code, CharPtr msg_txt);
void tax_setMsgHndl(tax_MessageHandler f);

#define tax_getDbType(x) ((x)->DbType)
#define tax_getMarker(x) ((x)->DbFieldMarker)
#define tax_getFileName(x, n) fs_fileName((x)->path, (x)->FileName[n])

/*-------------------------------------
 * rank.c
 */
void tax_freeRanks(_taxDBCtlPtr ctl);
int tax_loadRanks(_taxDBCtlPtr ctl);
CharPtr tax_getRankById(_taxDBCtlPtr ctl, Int2 id);
Int2 tax_getRankId(_taxDBCtlPtr ctl, CharPtr rank_name);

/*-------------------------------------
 * nmclass.c
 */
void tax_freeClasses(_taxDBCtlPtr ctl);
int tax_loadClasses(_taxDBCtlPtr ctl);
CharPtr tax_getClassById(_taxDBCtlPtr ctl, Int2 id);
Int2 tax_getClassId(_taxDBCtlPtr ctl, CharPtr class_name);


/*-------------------------------------
 * division.c
 */
void tax_freeDivs(_taxDBCtlPtr ctl);
int tax_loadDivs(_taxDBCtlPtr ctl);
CharPtr tax_getDivById(_taxDBCtlPtr ctl, Int2 id, Int2 mode);
Int2 tax_getDivId(_taxDBCtlPtr ctl, CharPtr txt, Int2 mode);

/*-------------------------------------
 * gc.c
 */
void tax_freeGCs(_taxDBCtlPtr ctl);
int tax_loadGCs(_taxDBCtlPtr ctl);
CharPtr tax_getGCById(_taxDBCtlPtr ctl, Int2 id, Int2 mode);
Int2 tax_getGCId(_taxDBCtlPtr ctl, CharPtr txt, Int2 mode);

/*-------------------------------------
 * tree.c
 */
void tax_freeTree(_taxDBCtlPtr ctl);
int tax_loadTree(_taxDBCtlPtr ctl);
void tax_freeNode(_taxNodePtr node);
_taxNodePtr tax_getNode(_taxDBCtlPtr ctl, Int4 id);
Int2 tax_getLin(_taxDBCtlPtr ctl, Int4 id, Int4Ptr lin);
#define tax_getParent(c, i) ((c)->tree[i].parent)
#define tax_getChild(c, i) ((c)->tree[i].child)
#define tax_getSibling(c, i) ((c)->tree[i].sibling)
#define tax_getNodeKey(c, i) ((c)->tree[i].dataKey)
#define tax_getNamesKey(c, i) ((c)->tree[i].namesKey)


/*-------------------------------------
 * taxname.c
 */
void tax_freeNames(_taxDBCtlPtr ctl);
int tax_loadNames(_taxDBCtlPtr ctl);
int tax_findByName(_taxDBCtlPtr ctl, CharPtr sname, int mode, _taxNamePtr** nameList);
_taxNamePtr tax_findOrgName(_taxDBCtlPtr ctl, Int4 id);
Int4 tax_getDesignator(CharPtr name);

#endif
