/****************************
 * File: rank.c
 * Description: taxonomy rank functions
 */

#include <stdlib.h>
#include "tax_cmmn.h"

/*------------------------------------------------
 * load_bin - load ranks from binary file
 */
static int load_bin(_taxDBCtlPtr ctl)
{
  FILE* f;
  int i;
  Int2 t;
  _rankPtr rank;
  int nof_ranks= tax_getNofRanks(ctl);

  if(nof_ranks < 1) {
    tax_outMsg(TAX_WARN, TAX_WARN_NO_RANKS, "load_bin: The number of ranks <= 1");
    return 1;
  }
	
  if((f= fopen(tax_getFileName(ctl, TAX_RANK_FILE), "rb")) == NULL) {
    tax_outMsg(TAX_ERROR, TAX_ERR_NO_FILE, "load_bin: Can't open rank file");
    return -1;
  }

  if((rank= malloc(sizeof(_rank) * nof_ranks)) == NULL) {
    tax_outMsg(TAX_ERROR, TAX_ERR_NO_MEMORY, "load_bin: Not enaugh memory for ranks");
    return -1;
  }

  /* read file and fill-up array */
  for(i= 0; i < nof_ranks; i++) {
    rank[i].rank_id= fsb_getInt2(f);
    rank[i].rank_txt= fsb_getString(f);
  }

  ctl->rank= rank;
  fclose(f);
  return 0;
}

/*---------------------------------------------------
 * load_txt - load ranks from text file
 */
static int load_txt(_taxDBCtlPtr ctl)
{
  FILE* f;
  int i;
  _rankPtr rank;
  int nof_ranks= tax_getNofRanks(ctl);


  if(nof_ranks < 1) {
    tax_outMsg(TAX_WARN, TAX_WARN_NO_RANKS, "load_txt: The number of ranks <= 1");
    return 1;
  }
	
  if((f= fopen(tax_getFileName(ctl, TAX_RANK_FILE), "r")) == NULL) {
    tax_outMsg(TAX_ERROR, TAX_ERR_NO_FILE, "load_txt: Can't open rank file");
    return -1;
  }

  if((rank= malloc(sizeof(_rank) * nof_ranks)) == NULL) {
    tax_outMsg(TAX_ERROR, TAX_ERR_NO_MEMORY, "load_txt: Not enaugh memory for ranks");
    return -1;
  }

  /* read file and fill-up array */
  for(i= 0; i < nof_ranks; i++) {
    rank[i].rank_id= fs_getInt(f, tax_getMarker(ctl));
    rank[i].rank_txt= fs_getString(f, tax_getMarker(ctl));
  }

  ctl->rank= rank;
  fclose(f);
  return 0;
}

static int load_asn(_taxDBCtlPtr ctl)
{
    return -10; /* not done yet */
}

/*==============================================================
 * tax_freeRanks - free rank's memory
 */
void tax_freeRanks(_taxDBCtlPtr ctl)
{
  int nof_ranks= tax_getNofRanks(ctl);
  int i;
  _rankPtr r= ctl->rank;

  if((r == NULL) || (nof_ranks < 1)) return;

  for(i= 0; i != nof_ranks; i++) {
    if(r[i].rank_txt != NULL) free(r[i].rank_txt);
  }
  free(r);
  ctl->rank= NULL;
}

int tax_loadRanks(_taxDBCtlPtr ctl)
{

  if(ctl->rank != NULL) tax_freeRanks(ctl);

  switch(tax_getDbType(ctl)) {
  case TAX_BIN: return load_bin(ctl);
  case TAX_TXT: return load_txt(ctl);
  case TAX_ASN: return load_asn(ctl);
  default: break;
  }

  tax_outMsg(TAX_ERROR, TAX_ERR_FILE_TYPE, "tax_loadRanks: Wrong file type for rank");
  return -1;
}

char* tax_getRankById(_taxDBCtlPtr ctl, Int2 id)
{
  int i;
  _rankPtr r= ctl->rank;
  int nof_ranks= tax_getNofRanks(ctl);

  for(i= 0; i < nof_ranks; i++) {
    if(r[i].rank_id == id) return r[i].rank_txt;
  }

  return NULL;
}

Int2 tax_getRankId(_taxDBCtlPtr ctl, CharPtr rank_name)
{
  int i;
  _rankPtr r= ctl->rank;
  int nof_ranks= tax_getNofRanks(ctl);

  for(i= 0; i < nof_ranks; i++) {
    if(strcmp(r[i].rank_txt, rank_name) == 0) {
      return r[i].rank_id;
    }
  }
  return -1;
}
  



