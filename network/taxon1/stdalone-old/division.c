/****************************
 * File: division.c
 * Description: taxonomy division functions
 */

#include <stdlib.h>
#include "tax_cmmn.h"

/*------------------------------------------------
 * load_bin - load divisions from binary file
 */
static int load_bin(_taxDBCtlPtr ctl)
{
  FILE* f;
  int i;
  _divisionPtr div;
  int nof_divs= tax_getNofDivs(ctl);

  if(nof_divs < 1) {
    tax_outMsg(TAX_WARN, TAX_WARN_NO_DIVS, "load_bin: The number of divisions <= 1");
    return 1;
  }
	
  if((f= fopen(tax_getFileName(ctl, TAX_DIVS_FILE), "rb")) == NULL) {
    tax_outMsg(TAX_ERROR, TAX_ERR_NO_FILE, "load_bin: Can't open divisions file");
    return -1;
  }

  if((div= malloc(sizeof(_division) * nof_divs)) == NULL) {
    tax_outMsg(TAX_ERROR, TAX_ERR_NO_MEMORY, "load_bin: Not enaugh memory for divisions");
    return -1;
  }

  /* read file and fill-up array */
  for(i= 0; i < nof_divs; i++) {
    div[i].div_id= fsb_getInt2(f);
    fsb_getBytes(f, 4, div[i].div_cde);
    div[i].div_txt= fsb_getString(f);
    div[i].div_comments_txt= fsb_getString(f);
  }

  ctl->division= div;
  fclose(f);
  return 0;
}

/*---------------------------------------------------
 * load_txt - load divisions from text file
 */
static int load_txt(_taxDBCtlPtr ctl)
{
  FILE* f;
  int i;
  _divisionPtr div;
  int nof_divs= tax_getNofDivs(ctl);
  CharPtr tmp;


  if(nof_divs < 1) {
    tax_outMsg(TAX_WARN, TAX_WARN_NO_RANKS, "load_txt: The number of divisions <= 1");
    return 1;
  }
	
  if((f= fopen(tax_getFileName(ctl, TAX_DIVS_FILE), "r")) == NULL) {
    tax_outMsg(TAX_ERROR, TAX_ERR_NO_FILE, "load_txt: Can't open divisions file");
    return -1;
  }

  if((div= malloc(sizeof(_division) * nof_divs)) == NULL) {
    tax_outMsg(TAX_ERROR, TAX_ERR_NO_MEMORY, "load_txt: Not enaugh memory for divisions");
    return -1;
  }

  /* read file and fill-up array */
  for(i= 0; i < nof_divs; i++) {
    div[i].div_id= fs_getInt(f, tax_getMarker(ctl));
    if((tmp= fs_getString(f, tax_getMarker(ctl))) != NULL) {
      strncpy(div[i].div_cde, tmp, 4);
      div[i].div_cde[3]= '\0';
      free(tmp);
    }
    div[i].div_txt= fs_getString(f, tax_getMarker(ctl));
    div[i].div_comments_txt= fs_getString(f, tax_getMarker(ctl)); 
  }

  ctl->division= div;
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
void tax_freeDivs(_taxDBCtlPtr ctl)
{
  int nof_divs= tax_getNofDivs(ctl);
  int i;
  _divisionPtr div= ctl->division;

  if((div == NULL) || (nof_divs < 1)) return;

  for(i= 0; i != nof_divs; i++) {
    if(div[i].div_txt != NULL) free(div[i].div_txt);
    if(div[i].div_comments_txt != NULL) free(div[i].div_comments_txt);
  }
  free(div);
  ctl->division= NULL;
}

int tax_loadDivs(_taxDBCtlPtr ctl)
{

  if(ctl->division != NULL) tax_freeDivs(ctl);

  switch(tax_getDbType(ctl)) {
  case TAX_BIN: return load_bin(ctl);
  case TAX_TXT: return load_txt(ctl);
  case TAX_ASN: return load_asn(ctl);
  default: break;
  }

  tax_outMsg(TAX_ERROR, TAX_ERR_FILE_TYPE, "tax_loadDivs: Wrong file type for divisions");
  return -1;
}

CharPtr tax_getDivById(_taxDBCtlPtr ctl, Int2 id, Int2 mode)
{
  int i;
  _divisionPtr div= ctl->division;
  int nof_divs= tax_getNofDivs(ctl);

  for(i= 0; i < nof_divs; i++) {
    if(div[i].div_id == id) {
      switch(mode) {
      case TAX_DIV_CDE: return div[i].div_cde;
      case TAX_DIV_COM: return div[i].div_comments_txt;
      default:          return div[i].div_txt;
      }
    }
  }

  return NULL;
}

Int2 tax_getDivId(_taxDBCtlPtr ctl, CharPtr txt, Int2 mode)
{
  int i;
  _divisionPtr div= ctl->division;
  int nof_divs= tax_getNofDivs(ctl);

  for(i= 0; i < nof_divs; i++) {
    if(mode == TAX_DIV_CDE) {
      if(StringNICmp(div[i].div_cde, txt, 3) == 0) return div[i].div_id;
    }
    else {
      if(StringICmp(div[i].div_txt, txt) == 0) return div[i].div_id;
    }
  }
  return -1;
}
  



