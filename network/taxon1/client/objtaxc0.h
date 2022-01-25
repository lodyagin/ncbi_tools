/*
*
*
* RCS Modification History:
* $Log: objtaxc0.h,v $
* Revision 6.0  1997/08/25 18:41:16  madden
* Revision changed to 6.0
*
* Revision 5.0  1996/05/28 14:15:13  ostell
* Set to revision 5.0
*
 * Revision 1.1  1996/03/06  17:05:19  soussov
 * Initial revision
 *
 * Revision 4.0  1995/07/26  13:55:46  ostell
 * force revision to 4.0
 *
 * Revision 1.4  1995/05/17  17:58:27  epstein
 * add RCS log revision history
 *
*/



/**************************************************
*
*    Generated objects for Module NCBI-Taxon0
*
**************************************************/


/**************************************************
*
*    Taxon0Req
*
**************************************************/
typedef ValNode  Taxon0Req;
typedef ValNodePtr Taxon0ReqPtr;
#define Taxon0Req_init 1 	/*  	NULL  	*/
#define Taxon0Req_getid 2 	/*  	Taxon-name  	*/
#define Taxon0Req_getref 3 	/*  	INTEGER  	*/
#define Taxon0Req_getchildren 4 	/*  	INTEGER  	*/
#define Taxon0Req_getparents 5 	/*  	INTEGER  	*/
#define Taxon0Req_getgeneticcode 6 	/*  	INTEGER  	*/
#define Taxon0Req_gettaxonline 7 	/*  	INTEGER  	*/
#define Taxon0Req_getdivision 8 	/*  	INTEGER  	*/
#define Taxon0Req_getcomplete 9 	/*  	Taxon-id-name  	*/
#define Taxon0Req_fini 10 	/*  	NULL  	*/


Taxon0ReqPtr LIBCALL Taxon0ReqFree PROTO ((Taxon0ReqPtr ));
Taxon0ReqPtr LIBCALL Taxon0ReqAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL Taxon0ReqAsnWrite PROTO (( Taxon0ReqPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    TaxonName
*
**************************************************/
typedef ValNode  TaxonName;
typedef ValNodePtr TaxonNamePtr;
#define TaxonName_taxname 1 	/*  	VisibleString  	*/
#define TaxonName_common 2 	/*  	VisibleString  	*/
#define TaxonName_tax_synonym 3 	/*  	VisibleString  	*/
#define TaxonName_com_synonym 4 	/*  	VisibleString  	*/


TaxonNamePtr LIBCALL TaxonNameFree PROTO ((TaxonNamePtr ));
TaxonNamePtr LIBCALL TaxonNameAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL TaxonNameAsnWrite PROTO (( TaxonNamePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    TaxonIdName
*
**************************************************/
typedef ValNode  TaxonIdName;
typedef ValNodePtr TaxonIdNamePtr;
#define TaxonIdName_id 1 	/*  	INTEGER  	*/
#define TaxonIdName_name 2 	/*  	VisibleString  	*/


TaxonIdNamePtr LIBCALL TaxonIdNameFree PROTO ((TaxonIdNamePtr ));
TaxonIdNamePtr LIBCALL TaxonIdNameAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL TaxonIdNameAsnWrite PROTO (( TaxonIdNamePtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    Taxon0Resp
*
**************************************************/
typedef ValNode  Taxon0Resp;
typedef ValNodePtr Taxon0RespPtr;
#define Taxon0Resp_error 1 	/*  	INTEGER  	*/
#define Taxon0Resp_init 2 	/*  	NULL  	*/
#define Taxon0Resp_getid 3 	/*  	Taxon-id-list  	*/
#define Taxon0Resp_getref 4 	/*  	Org-ref  	*/
#define Taxon0Resp_gettaxon 5 	/*  	Taxon-id-list  	*/
#define Taxon0Resp_getgeneticcode 6 	/*  	GeneticCodeList  	*/
#define Taxon0Resp_gettaxonline 7 	/*  	VisibleString  	*/
#define Taxon0Resp_getdivision 8 	/*  	VisibleString  	*/
#define Taxon0Resp_getcomplete 9 	/*  	Tax-complete-list  	*/
#define Taxon0Resp_fini 10 	/*  	NULL  	*/


Taxon0RespPtr LIBCALL Taxon0RespFree PROTO ((Taxon0RespPtr ));
Taxon0RespPtr LIBCALL Taxon0RespAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL Taxon0RespAsnWrite PROTO (( Taxon0RespPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    TaxonIdList
*
**************************************************/
typedef struct struct_Taxon_id_list {
	Int4   _num_ids;
	Int4 PNTR   ids;
} TaxonIdList, PNTR TaxonIdListPtr;


TaxonIdListPtr LIBCALL TaxonIdListFree PROTO ((TaxonIdListPtr ));
TaxonIdListPtr LIBCALL TaxonIdListNew PROTO (( void ));
TaxonIdListPtr LIBCALL TaxonIdListAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL TaxonIdListAsnWrite PROTO (( TaxonIdListPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    GeneticCodeList
*
**************************************************/
typedef struct struct_GeneticCodeList {
	Int4   genomic;
	Int4   mitochondrial;
} GeneticCodeList, PNTR GeneticCodeListPtr;


GeneticCodeListPtr LIBCALL GeneticCodeListFree PROTO ((GeneticCodeListPtr ));
GeneticCodeListPtr LIBCALL GeneticCodeListNew PROTO (( void ));
GeneticCodeListPtr LIBCALL GeneticCodeListAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL GeneticCodeListAsnWrite PROTO (( GeneticCodeListPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    TaxonId
*
**************************************************/
typedef struct struct_Taxon_id {
	Int4   id;
} TaxonId, PNTR TaxonIdPtr;


TaxonIdPtr LIBCALL TaxonIdFree PROTO ((TaxonIdPtr ));
TaxonIdPtr LIBCALL TaxonIdNew PROTO (( void ));
TaxonIdPtr LIBCALL TaxonIdAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL TaxonIdAsnWrite PROTO (( TaxonIdPtr , AsnIoPtr, AsnTypePtr));



/**************************************************
*
*    TaxComplete
*
**************************************************/
typedef struct struct_Tax_complete {
	 struct struct_Tax_complete PNTR next;
	CharPtr   sciname;
	CharPtr   comname;
	CharPtr   synonyms;
	Int4   id_gc;
	CharPtr   name_gc;
	Int4   id_mgc;
	CharPtr   name_mgc;
	CharPtr   gb_div;
	CharPtr   embl_code;
	CharPtr   lineage;
	Uint1   is_species_level;
} TaxComplete, PNTR TaxCompletePtr;


TaxCompletePtr LIBCALL TaxCompleteFree PROTO ((TaxCompletePtr ));
TaxCompletePtr LIBCALL TaxCompleteNew PROTO (( void ));
TaxCompletePtr LIBCALL TaxCompleteAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL TaxCompleteAsnWrite PROTO (( TaxCompletePtr , AsnIoPtr, AsnTypePtr));


/**************************************************
*
*    TaxCompleteList
*
**************************************************/
typedef struct struct_Tax_complete_list {
	Int4   num;
	TaxCompletePtr   info;
} TaxCompleteList, PNTR TaxCompleteListPtr;


TaxCompleteListPtr LIBCALL TaxCompleteListFree PROTO ((TaxCompleteListPtr ));
TaxCompleteListPtr LIBCALL TaxCompleteListNew PROTO (( void ));
TaxCompleteListPtr LIBCALL TaxCompleteListAsnRead PROTO (( AsnIoPtr, AsnTypePtr));
Boolean LIBCALL TaxCompleteListAsnWrite PROTO (( TaxCompleteListPtr , AsnIoPtr, AsnTypePtr));

