
#ifndef _TAXUTIL_H_
#define _TAXUTIL_H_ taxutil

#include <ncbi.h>
#include <asn.h>
#include <sequtil.h>
#include <taxinc.h>

#define MAXIDLIST 50  

typedef struct taxonomy_block {
	Int4  			hits;
	Taxon1DataPtr	tax;
} TaxBlk, PNTR TaxBlkPtr;

typedef struct struct_GeneticCodeList {
	Int4   genomic;
	Int4   mitochondrial;
} GeneticCodeList, PNTR GeneticCodeListPtr;

void tax_init PROTO((void));
Int4 taxname_replace PROTO((CharPtr iname, Taxon1DataPtr new));
Int4 taxname_match PROTO((CharPtr orgname, Boolean err));
OrgRefPtr check_org_ref PROTO((OrgRefPtr orp, Boolean replace));
OrgRefPtr get_tax_org PROTO((CharPtr name));
CharPtr get_lineage PROTO((CharPtr name));
GeneticCodeListPtr get_gcode PROTO((CharPtr name));
GeneticCodeListPtr get_gcode_from_lineage PROTO((CharPtr name));
OrgRefPtr replace_org PROTO((OrgRefPtr orp, Boolean replace));
OrgRefPtr replace_org_err PROTO((OrgRefPtr orp, Boolean replace));
CharPtr get_tax_division PROTO((OrgRefPtr orp));
CharPtr get_embl_code PROTO((OrgRefPtr orp));
OrgRefPtr check_org_ref PROTO((OrgRefPtr orp, Boolean replace));
void GetTaxserverOrg PROTO((SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent));
Boolean CheckTaxId PROTO((SeqEntryPtr sep));

/* temporary simulant functions for old-new Taxon switch period */
Boolean TaxArchInit PROTO((void));
Boolean TaxArchFini PROTO((void));
void TaxMergeBSinDescr PROTO((SeqEntryPtr sep, Pointer data, Int4 index, Int2 indent));

#endif /* _TAXUTIL_H_ */
