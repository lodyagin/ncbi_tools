/*****************************************
 * File: taxinc.h
 * Description: taxon1 API
 */

#ifndef TAXINC_H_DONE
#define TAXINC_H_DONE

#include <ncbi.h>
#include <objfeat.h>

#define struct_Org_ref orgref

#include "objtax1.h"

/* API functions prototypes */

/*---------------------------------------------
 * Taxon1 server init
 * Returns: TRUE - OK
 *          FALSE - Can't open taxonomy database
 */
Boolean tax1_init();

/*---------------------------------------------
 * Taxon1 server fini (frees memory)
 */

void tax1_fini();

/*---------------------------------------------
 * Get organism by tax_id
 * Returns: pointer to Taxon1Data if organism exists
 *          NULL - if tax_id wrong
 *
 * typedef struct struct_Taxon1_data {
 *   OrgRefPtr   org;
 *   CharPtr   div;                 genbank division ("MAM", "VIR", "PHG", ...)
 *   CharPtr   embl_code;           in common case first letters from genus-species pair
 *   Uint1   is_species_level;      1 if node on species level or below
 * } Taxon1Data, PNTR Taxon1DataPtr;
 * 
 * NOTE:
 * Caller gets his own copy of OrgRef structure.
 */

Taxon1DataPtr tax1_getbyid(Int4 tax_id);

/*----------------------------------------------
 * Get organism by OrgRef
 * Returns: pointer to Taxon1Data if organism exists
 *          NULL - if no such organism in taxonomy database
 *
 * NOTE:
 * 1. This functions uses the following data from inp_orgRef to find organism
 *    in taxonomy database. It uses taxname first. If no organism was found (or multiple nodes found)
 *    then it tryes to find organism using common name. If nothing found, then it tryes to find
 *    organism using synonyms. tax1_lookup never uses tax_id to find organism.
 * 2. If merge == 0 this function makes the new copy of OrgRef structure and puts into it data
 *    from taxonomy database.
 *    If merge != 0 this function updates inp_orgRef structure.
 */
Taxon1DataPtr tax1_lookup(OrgRefPtr inp_orgRef, int merge);

/*-----------------------------------------------
 * Get tax_id by OrgRef
 * Returns: tax_id - if organism found
 *          0 - no organism found
 *          -tax_id - if multiple nodes found (where tax_id is id of one of the nodes)
 * NOTE:
 * This function uses the same information from inp_orgRef as a tax1_lookup
 */ 
Int4 tax1_getTaxIdByOrgRef(OrgRefPtr inp_orgRef);

/*----------------------------------------------
 * Get tax_id by organism name
 * Returns: tax_id - if organism found
 *          0 - no organism found
 *          -tax_id - if multiple nodes found (where tax_id is id of one of the nodes)
 * NOTE:
 * orgname can be a regular expression.
 */
Int4 tax1_getTaxIdByName(CharPtr orgname);


/*----------------------------------------------
 * Get ALL tax_id by organism name
 * Returns: number of organism found
 * NOTE:
 * 1. orgname can be a regular expression.
 * 2. Ids consists of tax ids. Caller is responsible to free this memory
 */
Int4 tax1_getAllTaxIdByName(CharPtr orgname, Int4 **Ids);

/*----------------------------------------------
 * Get organism by tax_id
 * Returns: pointer to OrgRef structure if organism found
 *          NULL - if no such organism in taxonomy database
 * NOTE:
 * This function does not make a copy of OrgRef structure, so, caller can not use
 * OrgRefFree function for returned pointer.
 */
OrgRefPtr tax1_getOrgRef(Int4 tax_id, int* is_species, CharPtr div, CharPtr embl_cde);

/*---------------------------------------------
 * Set mode for synonyms in OrgRef
 * Returns: previous mode
 * NOTE:
 * Default mode: do not include synonyms in OrgRef
 */
Boolean tax1_setSynonyms(Boolean on_off);

/*---------------------------------------------
 * Get parent tax_id
 */
Int4 tax1_getParent(Int4 id_tax);

/*---------------------------------------------
 * Get taxids for all children.
 * Returns: number of children
 */
int tax1_getChildren(Int4 id_tax, Int4** children_ids);

/*---------------------------------------------
 * Get genetic code name by genetic code id
 */
CharPtr tax1_getGCName(Int2 gc_id);


#endif