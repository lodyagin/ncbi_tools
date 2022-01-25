/****************************
 * File: nmclass.c
 * Description: taxonomy name class functions
 */

#include <stdlib.h>
#include "tax_cmmn.h"

/*------------------------------------------------
 * load_bin - load name classes from binary file
 */
static int load_bin(_taxDBCtlPtr ctl)
{
  FILE* f;
  int i;
  Int2 t;
  _nmclassPtr nmc;
  int nof_nmcs= tax_getNofClasses(ctl);

  if(nof_nmcs < 1) {
    tax_outMsg(TAX_WARN, TAX_WARN_NO_CLASSES, "load_bin: The number of name classes <= 1");
    return 1;
  }
	
  if((f= fopen(tax_getFileName(ctl, TAX_CLASS_FILE), "rb")) == NULL) {
    tax_outMsg(TAX_ERROR, TAX_ERR_NO_FILE, "load_bin: Can't open name classes file");
    return -1;
  }

  if((nmc= malloc(sizeof(_nmclass) * nof_nmcs)) == NULL) {
    tax_outMsg(TAX_ERROR, TAX_ERR_NO_MEMORY, "load_bin: Not enaugh memory for name classes");
    return -1;
  }

  /* read file and fill-up array */
  for(i= 0; i < nof_nmcs; i++) {
    nmc[i].class_cde= fsb_getInt2(f);
    nmc[i].class_txt= fsb_getString(f);
  }

  ctl->nmclass= nmc;
  fclose(f);
  return 0;
}

/*---------------------------------------------------
 * load_txt - load name classes from text file
 */
static int load_txt(_taxDBCtlPtr ctl)
{
  FILE* f;
  int i;
  _nmclassPtr nmc;
  int nof_nmcs= tax_getNofClasses(ctl);


  if(nof_nmcs < 1) {
    tax_outMsg(TAX_WARN, TAX_WARN_NO_CLASSES, "load_txt: The number of name classes <= 1");
    return 1;
  }

  if((f= fopen(tax_getFileName(ctl, TAX_CLASS_FILE), "r")) == NULL) {
    tax_outMsg(TAX_ERROR, TAX_ERR_NO_FILE, "load_txt: Can't open name classes file");
    return -1;
  }

  if((nmc= malloc(sizeof(_nmclass) * nof_nmcs)) == NULL) {
    tax_outMsg(TAX_ERROR, TAX_ERR_NO_MEMORY, "load_txt: Not enaugh memory for name classes");
    return -1;
  }

  /* read file and fill-up array */
  for(i= 0; i < nof_nmcs; i++) {
    nmc[i].class_cde= fs_getInt(f, tax_getMarker(ctl));
    nmc[i].class_txt= fs_getString(f, tax_getMarker(ctl));
  }

  ctl->nmclass= nmc;
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
void tax_freeClasses(_taxDBCtlPtr ctl)
{
  int nof_nmcs= tax_getNofClasses(ctl);
  int i;
  _nmclassPtr r= ctl->nmclass;

  if((r == NULL) || (nof_nmcs < 1)) return;

  for(i= 0; i != nof_nmcs; i++) {
    if(r[i].class_txt != NULL) free(r[i].class_txt);
  }
  free(r);
  ctl->nmclass= NULL;
}

int tax_loadClasses(_taxDBCtlPtr ctl)
{

  if(ctl->nmclass != NULL) tax_freeClasses(ctl);

  switch(tax_getDbType(ctl)) {
  case TAX_BIN: return load_bin(ctl);
  case TAX_TXT: return load_txt(ctl);
  case TAX_ASN: return load_asn(ctl);
  default: break;
  }

  tax_outMsg(TAX_ERROR, TAX_ERR_FILE_TYPE, "tax_loadClasses: Wrong file type for name classes");
  return -1;
}

char* tax_getClassById(_taxDBCtlPtr ctl, Int2 id)
{
  int i;
  _nmclassPtr r= ctl->nmclass;
  int nof_nmcs= tax_getNofClasses(ctl);

  for(i= 0; i < nof_nmcs; i++) {
    if(r[i].class_cde == id) return r[i].class_txt;
  }

  return NULL;
}

Int2 tax_getClassId(_taxDBCtlPtr ctl, CharPtr class_name)
{
  int i;
  _nmclassPtr r= ctl->nmclass;
  int nof_nmcs= tax_getNofClasses(ctl);

  for(i= 0; i < nof_nmcs; i++) {
    if(strcmp(r[i].class_txt, class_name) == 0) {
      return r[i].class_cde;
    }
  }
  return 99;
}
  



