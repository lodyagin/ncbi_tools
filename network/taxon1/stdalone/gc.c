/****************************
 * File: gc.c
 * Description: taxonomy GC functions
 */

#include <stdlib.h>
#include "tax_cmmn.h"

/*------------------------------------------------
 * load_bin - load GCs from binary file
 */
static int load_bin(_taxDBCtlPtr ctl)
{
  FILE* f;
  int i;
  _gencodePtr GC;
  int nof_GCs= tax_getNofGCs(ctl);

  if(nof_GCs < 1) {
    tax_outMsg(TAX_WARN, TAX_WARN_NO_GC, "load_bin: The number of GCs <= 1");
    return 1;
  }
	
  if((f= fopen(tax_getFileName(ctl, TAX_GC_FILE), "rb")) == NULL) {
    tax_outMsg(TAX_ERROR, TAX_ERR_NO_FILE, "load_bin: Can't open GC file");
    return -1;
  }

  if((GC= malloc(sizeof(_gencode) * nof_GCs)) == NULL) {
    tax_outMsg(TAX_ERROR, TAX_ERR_NO_MEMORY, "load_bin: Not enaugh memory for GCs");
    return -1;
  }

  /* read file and fill-up array */
  for(i= 0; i < nof_GCs; i++) {
    GC[i].gc_id= fsb_getInt2(f);
    fsb_getBytes(f, 6, GC[i].gc_abbrev);
    GC[i].gc_nm= fsb_getString(f);
    GC[i].gc_cde= fsb_getString(f);
    GC[i].gc_starts= fsb_getString(f);
  }

  ctl->gc= GC;
  fclose(f);
  return 0;
}

/*---------------------------------------------------
 * load_txt - load GC from text file
 */
static int load_txt(_taxDBCtlPtr ctl)
{
  FILE* f;
  int i;
  _gencodePtr GC;
  int nof_GCs= tax_getNofGCs(ctl);
  CharPtr tmp;


  if(nof_GCs < 1) {
    tax_outMsg(TAX_WARN, TAX_WARN_NO_RANKS, "load_txt: The number of GCs <= 1");
    return 1;
  }
	
  if((f= fopen(tax_getFileName(ctl, TAX_GC_FILE), "r")) == NULL) {
    tax_outMsg(TAX_ERROR, TAX_ERR_NO_FILE, "load_txt: Can't open GC file");
    return -1;
  }

  if((GC= malloc(sizeof(_gencode) * nof_GCs)) == NULL) {
    tax_outMsg(TAX_ERROR, TAX_ERR_NO_MEMORY, "load_txt: Not enaugh memory for GCs");
    return -1;
  }

  /* read file and fill-up array */
  for(i= 0; i < nof_GCs; i++) {
    GC[i].gc_id= fs_getInt(f, tax_getMarker(ctl));
    if((tmp= fs_getString(f, tax_getMarker(ctl))) != NULL) {
      strncpy(GC[i].gc_abbrev, tmp, 6);
      GC[i].gc_abbrev[5]= '\0';
      free(tmp);
    }
    GC[i].gc_nm= fs_getString(f, tax_getMarker(ctl));
    GC[i].gc_cde= fs_getString(f, tax_getMarker(ctl)); 
    GC[i].gc_starts= fs_getString(f, tax_getMarker(ctl)); 
  }

  ctl->gc= GC;
  fclose(f);
  return 0;
}

static int load_asn(_taxDBCtlPtr ctl)
{
    return -10; /* not done yet */
}

/*==============================================================
 * tax_freeDivs - free division's memory
 */
void tax_freeGCs(_taxDBCtlPtr ctl)
{
  int nof_GCs= tax_getNofGCs(ctl);
  int i;
  _gencodePtr GC= ctl->gc;

  if((GC == NULL) || (nof_GCs < 1)) return;

  for(i= 0; i != nof_GCs; i++) {
    if(GC[i].gc_nm != NULL) free(GC[i].gc_nm);
    if(GC[i].gc_cde != NULL) free(GC[i].gc_cde);
    if(GC[i].gc_starts != NULL) free(GC[i].gc_starts);
  }
  free(GC);
  ctl->gc= NULL;
}

int tax_loadGCs(_taxDBCtlPtr ctl)
{

  if(ctl->gc != NULL) tax_freeGCs(ctl);

  switch(tax_getDbType(ctl)) {
  case TAX_BIN: return load_bin(ctl);
  case TAX_TXT: return load_txt(ctl);
  case TAX_ASN: return load_asn(ctl);
  default: break;
  }

  tax_outMsg(TAX_ERROR, TAX_ERR_FILE_TYPE, "tax_loadDivs: Wrong file type for GCs");
  return -1;
}

CharPtr tax_getGCById(_taxDBCtlPtr ctl, Int2 id, Int2 mode)
{
  int i;
  _gencodePtr GC= ctl->gc;
  int nof_GCs= tax_getNofGCs(ctl);

  for(i= 0; i < nof_GCs; i++) {
    if(GC[i].gc_id == id) {
      switch(mode) {
      case TAX_GC_CDE:    return GC[i].gc_cde;
      case TAX_GC_ABBREV: return GC[i].gc_abbrev;
      case TAX_GC_STARTS: return GC[i].gc_starts;
      default:            return GC[i].gc_nm;
      }
    }
  }

  return NULL;
}

Int2 tax_getGCId(_taxDBCtlPtr ctl, CharPtr txt, Int2 mode)
{
  int i;
  _gencodePtr GC= ctl->gc;
  int nof_GCs= tax_getNofGCs(ctl);

  for(i= 0; i < nof_GCs; i++) {
    if(mode == TAX_GC_ABBREV) {
      if(StringNICmp(GC[i].gc_abbrev, txt, 6) == 0) return GC[i].gc_id;
    }
    else {
      if(StringICmp(GC[i].gc_nm, txt) == 0) return GC[i].gc_id;
    }
  }
  return -1;
}
  



