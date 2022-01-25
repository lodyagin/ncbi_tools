/*************************
 * File: tax1proc.c
 * Description: old stile client functions
 */

#include <taxinc.h>
#include <objsset.h>
#include <taxarch.h>

#if 0
Boolean TaxArchInit(void)
{
  return tax1_init();
}

Boolean TaxArchFini(void)
{
  tax1_fini();
  return TRUE;
}
#endif

TaxonIdListPtr TaxArchGetTaxId(TaxonNamePtr tnp)
{
  TaxonIdListPtr res= TaxonIdListNew();
  
  res->_num_ids= tax1_getAllTaxIdByName(tnp->data.ptrvalue, &res->ids);
  return res;
}

TaxonIdListPtr TaxArchGetChildren(Int4 id_tax)
{
  TaxonIdListPtr res= TaxonIdListNew();

  res->_num_ids= tax1_getChildren(id_tax, &res->ids);
  return res;
}

TaxonIdListPtr TaxArchGetParents(Int4 id_tax)
{
  TaxonIdListPtr res= TaxonIdListNew();
  Int4 pid;
  
  pid= tax1_getParent(id_tax);
  if(pid > 0) {
    res->_num_ids= 1;
    res->ids= MemNew(sizeof(Int4));
    *(res->ids)= pid;
  }
  else {
    res->_num_ids= 0;
    res->ids= NULL;
  }
  return res;
}

OrgRefPtr TaxArchGetRef(Int4 id_tax)
{
  Taxon1DataPtr res;
  OrgRefPtr orp;

  res= tax1_getbyid(id_tax);
  if(res == NULL) return NULL;
  orp= res->org;
  res->org= NULL;
  Taxon1DataFree(res);
  return orp;
}

CharPtr TaxArchGetTaxonLine(Int4 id_tax)
{
  OrgRefPtr orp;
  CharPtr lin;

  orp= tax1_getOrgRef(id_tax, NULL, NULL, NULL);
  if((orp == NULL) || (orp->orgname == NULL)) {
    return NULL;
  }
  lin= StringSave(orp->orgname->lineage);
  return lin;
}

CharPtr TaxArchGetWithinDiv(Int4 id_tax)
{
  CharPtr division= MemNew(8);
  OrgRefPtr orp= tax1_getOrgRef(id_tax, NULL, division, NULL);

  if(orp == NULL) {
    MemFree(division);
    return NULL;
  }

  return division;
}

GeneticCodeListPtr TaxArchGetGeneticCode(Int4 id_tax)
{
  OrgRefPtr orp;
  GeneticCodeListPtr gclp;

  orp= tax1_getOrgRef(id_tax, NULL, NULL, NULL);
  if((orp == NULL) || (orp->orgname == NULL)) {
    return NULL;
  }
  gclp= GeneticCodeListNew();
  gclp->genomic= orp->orgname->gcode;
  gclp->mitochondrial= orp->orgname->mgcode;
  return gclp;
}
  
TaxCompleteListPtr TaxArchGetComplete(TaxonIdNamePtr tinp)
{
  OrgRefPtr orp;
  TaxCompleteListPtr tclp;
  TaxCompletePtr tcp;
  Int4 id_tax;
  char embl[8], divis[8];
  int is_spec;

  if(tinp->choice == TaxonIdName_id) {
    id_tax= tinp->data.intvalue;
  }
  else if(tinp->choice == TaxonIdName_name) {
    id_tax= tax1_getTaxIdByName(tinp->data.ptrvalue);
  }
  else {
    id_tax= -1;
  }

  if(id_tax <= 0) return NULL;

  *embl= *divis= '\0';
  orp= tax1_getOrgRef(id_tax, &is_spec, divis, embl);
  if(orp == NULL) return NULL;

  tclp= TaxCompleteListNew();
  tclp->num= 1;
  tclp->info= tcp= TaxCompleteNew();

  tcp->next= NULL;
  tcp->sciname= (orp->taxname == NULL)? NULL : StringSave(orp->taxname);
  tcp->comname= (orp->common == NULL)? NULL : StringSave(orp->common);
  tcp->synonyms= NULL;
  tcp->gb_div= (*divis == '\0')? NULL : StringSave(divis);
  tcp->embl_code= (*embl == '\0')? NULL : StringSave(embl);
  tcp->is_species_level= is_spec;
  if(orp->orgname != NULL) {
    tcp->id_gc= orp->orgname->gcode;
    tcp->id_mgc= orp->orgname->mgcode;
    tcp->lineage= (orp->orgname->lineage == NULL)? NULL : StringSave(orp->orgname->lineage);
  }
  tcp->name_gc= tax1_getGCName(tcp->id_gc);
  tcp->name_mgc= tax1_getGCName(tcp->id_mgc);

  return tclp;
}
    


