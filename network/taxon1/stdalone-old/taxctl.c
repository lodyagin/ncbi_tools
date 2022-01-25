/**************************
 * File: taxctl.c
 * Description: catalog file utility
 */

#include "tax_cmmn.h"
#include "fsutils.h"

#define DEFAULT_MARKER "\t|\n"

static tax_MessageHandler msg_hndl= NULL;

void tax_setMsgHndl(tax_MessageHandler f)
{
  msg_hndl= f;
}

void tax_outMsg(int msg_type, int msg_code, CharPtr msg_txt)
{
  if(msg_hndl != NULL) {
    (*msg_hndl)(msg_type, msg_code, msg_txt);
  }
  else {
    char* et;

    et= (msg_type == TAX_MSG)? "NOTE" : ((msg_type == TAX_WARN)? "WARNING" : "ERROR");
      
    fprintf(stderr, "%s: File based taxonomy: %s\n", et, msg_txt);
    if(msg_type == TAX_ERROR) {
      fprintf(stderr, "The taxonomy files are corrupted\n");
      exit(1);
    }
  }
}

static char default_DB[]= "dbctl.tax";

_taxDBCtlPtr tax_loadDBCtl(CharPtr DBPath, CharPtr ctlFileName)
{
  FILE* f;
  _taxDBCtlPtr pctl;
  CharPtr resp;
  int i;

  if(ctlFileName == NULL) ctlFileName= default_DB;

  if((f= fopen(fs_fileName(DBPath, ctlFileName), "r")) == NULL) return NULL;

  if((pctl= malloc(sizeof(_taxDBCtl))) == NULL) return NULL;
  memset(pctl, 0, sizeof(_taxDBCtl));

  pctl->path= StringSave(DBPath);

  /* Get DB type (bin, txt or ASN */
  resp= fs_getString(f, CTL_FILE_MARKER);
  if(resp == NULL) {
    free(pctl);
    return NULL;
  }

  pctl->DbType= resp[0];
  free(resp);
  if(pctl->DbType != TAX_ASN) {
    /* not ASN.1 database */
    pctl->MaxId= fs_getLong(f, CTL_FILE_MARKER);
    pctl->NofNodes= fs_getLong(f, CTL_FILE_MARKER);
    pctl->NofNames= fs_getLong(f, CTL_FILE_MARKER);
    pctl->NofRanks= fs_getInt(f, CTL_FILE_MARKER);
    pctl->NofDivs= fs_getInt(f, CTL_FILE_MARKER);
    pctl->NofClasses= fs_getInt(f, CTL_FILE_MARKER);
    pctl->NofGCs= fs_getInt(f, CTL_FILE_MARKER);

    for(i= 0; i != NOF_FILES; i++) {
      pctl->FileName[i]= fs_getString(f, CTL_FILE_MARKER);
    }

    resp= fs_getString(f, CTL_FILE_MARKER);
    if((resp != NULL) && (resp[0] != '\0')) {
      strncpy(pctl->DbFieldMarker, resp, 6);
      pctl->DbFieldMarker[6]= '\0';
    }
    else {
      strcpy(pctl->DbFieldMarker, DEFAULT_MARKER);
    }

    if(resp != NULL) free(resp);
  }
  else {
    /* ASN.1 files */
    /* ... */
  }
  pctl->nodefile= NULL;

  return pctl;
}
    
Int4 tax_getMaxId(_taxDBCtlPtr ctl)
{
  if(ctl->MaxId > 0) return ctl->MaxId;

  return (ctl->DbType == TAX_ASN)? (-1) : 0;
}
    
Int4 tax_getNofNodes(_taxDBCtlPtr ctl)
{
  if(ctl->NofNodes > 0) return ctl->NofNodes;

  return (ctl->DbType == TAX_ASN)? (-1) : 0;
}
    
Int4 tax_getNofNames(_taxDBCtlPtr ctl)
{
  if(ctl->NofNames > 0) return ctl->NofNames;

  return (ctl->DbType == TAX_ASN)? (-1) : 0;
}
    
Int2 tax_getNofRanks(_taxDBCtlPtr ctl)
{
  if(ctl->NofRanks > 0) return ctl->NofRanks;

  return (ctl->DbType == TAX_ASN)? (-1) : 0;
}

Int2 tax_getNofDivs(_taxDBCtlPtr ctl)
{
  if(ctl->NofDivs > 0) return ctl->NofDivs;

  return (ctl->DbType == TAX_ASN)? (-1) : 0;
}

Int2 tax_getNofClasses(_taxDBCtlPtr ctl)
{
  if(ctl->NofClasses > 0) return ctl->NofClasses;

  return (ctl->DbType == TAX_ASN)? (-1) : 0;
}

Int2 tax_getNofGCs(_taxDBCtlPtr ctl)
{
  if(ctl->NofGCs > 0) return ctl->NofGCs;

  return (ctl->DbType == TAX_ASN)? (-1) : 0;
}

void tax_closeDbCtl(_taxDBCtlPtr ctl)
{
  int i;

  for(i= 0; i != NOF_FILES; i++) {
    if(ctl->FileName[i] != NULL) free(ctl->FileName[i]);
  }

  free(ctl);
}

/*------------------------------------------------
 * load_div - load divisions from binary file
 */
static int load_div(_taxDBCtlPtr ctl, FILE* f)
{
  int i;
  _divisionPtr div;
  int nof_divs= tax_getNofDivs(ctl);

  if(nof_divs < 1) {
    tax_outMsg(TAX_WARN, TAX_WARN_NO_DIVS, "load_bin: The number of divisions <= 1");
    return 1;
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
  return 0;
}

/*------------------------------------------------
 * load_gc - load GCs from binary file
 */
static int load_gc(_taxDBCtlPtr ctl, FILE* f)
{
  int i;
  _gencodePtr GC;
  int nof_GCs= tax_getNofGCs(ctl);

  if(nof_GCs < 1) {
    tax_outMsg(TAX_WARN, TAX_WARN_NO_GC, "load_bin: The number of GCs <= 1");
    return 1;
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
  return 0;
}

/*------------------------------------------------
 * load_classes - load name classes from binary file
 */
static int load_classes(_taxDBCtlPtr ctl, FILE* f)
{
  int i;
  Int2 t;
  _nmclassPtr nmc;
  int nof_nmcs= tax_getNofClasses(ctl);

  if(nof_nmcs < 1) {
    tax_outMsg(TAX_WARN, TAX_WARN_NO_CLASSES, "load_bin: The number of name classes <= 1");
    return 1;
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
  return 0;
}

/*------------------------------------------------
 * load_rank - load ranks from binary file
 */
static int load_rank(_taxDBCtlPtr ctl, FILE* f)
{
  int i;
  Int2 t;
  _rankPtr rank;
  int nof_ranks= tax_getNofRanks(ctl);

  if(nof_ranks < 1) {
    tax_outMsg(TAX_WARN, TAX_WARN_NO_RANKS, "load_bin: The number of ranks <= 1");
    return 1;
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
  return 0;
}


static char data_file[]= "dbdata.dat";

_taxDBCtlPtr tax_loadTaxDB(CharPtr DBPath, CharPtr ctlFileName)
{
  FILE* f;
  _taxDBCtlPtr pctl;
  CharPtr resp;
  int i;
  FILE* df;

  if(ctlFileName == NULL) ctlFileName= default_DB;

  if((f= fopen(fs_fileName(DBPath, ctlFileName), "r")) == NULL) return NULL;

  if((pctl= malloc(sizeof(_taxDBCtl))) == NULL) return NULL;
  memset(pctl, 0, sizeof(_taxDBCtl));

  pctl->path= StringSave(DBPath);

  /* Get DB type (bin, txt or ASN */
  resp= fs_getString(f, CTL_FILE_MARKER);
  if(resp == NULL) {
    free(pctl);
    return NULL;
  }

  pctl->DbType= resp[0];
  free(resp);
  if(pctl->DbType != TAX_ASN) {
    /* not ASN.1 database */
    pctl->MaxId= fs_getLong(f, CTL_FILE_MARKER);
    pctl->NofNodes= fs_getLong(f, CTL_FILE_MARKER);
    pctl->NofNames= fs_getLong(f, CTL_FILE_MARKER);
    pctl->NofRanks= fs_getInt(f, CTL_FILE_MARKER);
    pctl->NofDivs= fs_getInt(f, CTL_FILE_MARKER);
    pctl->NofClasses= fs_getInt(f, CTL_FILE_MARKER);
    pctl->NofGCs= fs_getInt(f, CTL_FILE_MARKER);

    for(i= 0; i != NOF_FILES; i++) {
      pctl->FileName[i]= fs_getString(f, CTL_FILE_MARKER);
    }

    resp= fs_getString(f, CTL_FILE_MARKER);
    if((resp != NULL) && (resp[0] != '\0')) {
      strncpy(pctl->DbFieldMarker, resp, 6);
      pctl->DbFieldMarker[6]= '\0';
    }
    else {
      strcpy(pctl->DbFieldMarker, DEFAULT_MARKER);
    }

    if(resp != NULL) free(resp);
  }
  else {
    /* ASN.1 files */
    /* ... */
  }
  pctl->nodefile= NULL;

  if((df= fopen(fs_fileName(DBPath, data_file), "rb")) != NULL) {
      load_div(pctl, df);
      load_gc(pctl, df);
      load_classes(pctl, df);
      load_rank(pctl, df);
      fclose(df);
  }
 
  return pctl;
}
