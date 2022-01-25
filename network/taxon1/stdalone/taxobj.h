/****************************
 * File: taxobj.h
 * Description: header file for taxonomy objects
 */


#ifndef TAXOBJ_H_DONE
#define TAXOBJ_H_DONE

#include <stdio.h>
#include <ncbi.h>

typedef unsigned char _TINY;
typedef Int4 _DBKEY;

#define INHERIT_DIV 0x1
#define INHERIT_GC 0x2
#define INHERIT_MGC 0x4
#define GB_HIDDEN 0x8
#define TAX_HIDDEN 0x10


typedef struct t_taxCursor {
  FILE* f;
  _DBKEY  position;
} _taxCursor;

_taxCursor* tq_openCursor(char* fileName);
void tq_closeCursor(_taxCursor* cursor);
int tq_setCursor(_taxCursor* cursor, _DBKEY pos);

typedef struct t_taxCtlNode {
  Int4 parent, child, sibling;
  _DBKEY dataKey;
  Int2 rank;
  Int2 flags;
} _taxCtlNode, PNTR _taxCtlNodePtr;

typedef struct t_taxName {
  int tax_id;
  _TINY class_cde;
  CharPtr name_txt;
  CharPtr unique_name;
} _taxName, PNTR _taxNamePtr;

typedef struct t_taxNode {
  int tax_id;
  char* comment;
  char embl_cde[4];
  _TINY division;
  _TINY gc;
  _TINY mgc;
  _TINY comment_size; /* number of 32 bytes pices in comment */
} _taxNode, PNTR _taxNodePtr;
  

int tq_loadNode(_taxCursor* cursor, _taxNode* node);
int tq_loadNodeNames(_taxCursor* cursor, _taxName* name);


#endif
