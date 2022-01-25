/*
*
*
* RCS Modification History:
* $Log: objtax0.c,v $
* Revision 6.0  1997/08/25 18:42:20  madden
* Revision changed to 6.0
*
* Revision 1.1  1996/02/05 14:51:02  soussov
* Initial revision
*
*/

#include <asn.h>
#include <ncbi.h>
#include <objfeat.h>

#ifdef NLM_OBJ_INCL
#include NLM_OBJ_INCL
#endif
#include "objtax0.h"

static Boolean loaded = FALSE;

#include "asntaxon.h"

#ifndef NLM_EXTERN_LOADS
#define NLM_EXTERN_LOADS {}
#endif

static Boolean _AsnLoad(void)
{

	if ( ! loaded) {
		NLM_EXTERN_LOADS

		if ( ! AsnLoad ())
			return FALSE;
		loaded = TRUE;
	}

	return TRUE;
}



/**************************************************
*
*    Generated object loaders for Module NCBI-Taxon0
*
**************************************************/


/**************************************************
*
*    Taxon0ReqFree()
*
**************************************************/
Taxon0ReqPtr LIBCALL Taxon0ReqFree ( Taxon0ReqPtr ptr)
{
	if (ptr == NULL) return NULL;

	{Pointer pt = ptr -> data.ptrvalue;
	switch (ptr ->  choice){
		case  Taxon0Req_getid :
		TaxonNameFree(pt);
			break;
		case  Taxon0Req_getcomplete :
		TaxonIdNameFree(pt);
			break;
	}}
	MemFree(ptr);
	return NULL;
}


/**************************************************
*
*    TaxonNameFree()
*
**************************************************/
TaxonNamePtr LIBCALL TaxonNameFree ( TaxonNamePtr ptr)
{
	if (ptr == NULL) return NULL;

	{Pointer pt = ptr -> data.ptrvalue;
	switch (ptr ->  choice){
		case  TaxonName_taxname :
		case  TaxonName_common :
		case  TaxonName_tax_synonym :
		case  TaxonName_com_synonym :
		MemFree(pt);
			break;
	}}
	MemFree(ptr);
	return NULL;
}


/**************************************************
*
*    TaxonIdNameFree()
*
**************************************************/
TaxonIdNamePtr LIBCALL TaxonIdNameFree ( TaxonIdNamePtr ptr)
{
	if (ptr == NULL) return NULL;

	{Pointer pt = ptr -> data.ptrvalue;
	switch (ptr ->  choice){
		case  TaxonIdName_name :
		MemFree(pt);
			break;
	}}
	MemFree(ptr);
	return NULL;
}


/**************************************************
*
*    Taxon0RespFree()
*
**************************************************/
Taxon0RespPtr LIBCALL Taxon0RespFree ( Taxon0RespPtr ptr)
{
	if (ptr == NULL) return NULL;

	{Pointer pt = ptr -> data.ptrvalue;
	switch (ptr ->  choice){
		case  Taxon0Resp_getid :
		TaxonIdListFree(pt);
			break;
		case  Taxon0Resp_getref :
		OrgRefFree(pt);
			break;
		case  Taxon0Resp_gettaxon :
		TaxonIdListFree(pt);
			break;
		case  Taxon0Resp_getgeneticcode :
		GeneticCodeListFree(pt);
			break;
		case  Taxon0Resp_gettaxonline :
		case  Taxon0Resp_getdivision :
		MemFree(pt);
			break;
		case  Taxon0Resp_getcomplete :
		TaxCompleteListFree(pt);
			break;
	}}
	MemFree(ptr);
	return NULL;
}


/**************************************************
*
*    TaxonIdListFree()
*
**************************************************/
TaxonIdListPtr LIBCALL TaxonIdListFree ( TaxonIdListPtr ptr)
{
	if (ptr == NULL) return NULL;

	MemFree(ptr -> ids);
	MemFree(ptr);
	return NULL;
}


/**************************************************
*
*    GeneticCodeListFree()
*
**************************************************/
GeneticCodeListPtr LIBCALL GeneticCodeListFree ( GeneticCodeListPtr ptr)
{
	if (ptr == NULL) return NULL;

	MemFree(ptr);
	return NULL;
}


/**************************************************
*
*    TaxCompleteListFree()
*
**************************************************/
TaxCompleteListPtr LIBCALL TaxCompleteListFree ( TaxCompleteListPtr ptr)
{
	if (ptr == NULL) return NULL;

	TaxCompleteFree(  ptr -> info);
	MemFree(ptr);
	return NULL;
}


/**************************************************
*
*    TaxonIdFree()
*
**************************************************/
TaxonIdPtr LIBCALL TaxonIdFree ( TaxonIdPtr ptr)
{
	if (ptr == NULL) return NULL;

	MemFree(ptr);
	return NULL;
}


/**************************************************
*
*    TaxCompleteFree()
*
**************************************************/
TaxCompletePtr LIBCALL TaxCompleteFree ( TaxCompletePtr head_ptr)
{
	 TaxCompletePtr ptr, hold_ptr;

	if (head_ptr == NULL) return NULL;

 	for (ptr = head_ptr; ptr; ptr = hold_ptr){
		hold_ptr = ptr -> next;
		MemFree(ptr -> sciname);
		MemFree(ptr -> comname);
		MemFree(ptr -> synonyms);
		MemFree(ptr -> name_gc);
		MemFree(ptr -> name_mgc);
		MemFree(ptr -> gb_div);
		MemFree(ptr -> embl_code);
		MemFree(ptr -> lineage);
		MemFree(ptr);
	 }
	return NULL;
}


/**************************************************
*
*    TaxonIdListNew()
*
**************************************************/

TaxonIdListPtr LIBCALL
TaxonIdListNew()
{

	return (TaxonIdListPtr) MemNew(sizeof(TaxonIdList));

}



/**************************************************
*
*    GeneticCodeListNew()
*
**************************************************/

GeneticCodeListPtr LIBCALL
GeneticCodeListNew()
{

	return (GeneticCodeListPtr) MemNew(sizeof(GeneticCodeList));

}



/**************************************************
*
*    TaxCompleteListNew()
*
**************************************************/

TaxCompleteListPtr LIBCALL
TaxCompleteListNew()
{

	return (TaxCompleteListPtr) MemNew(sizeof(TaxCompleteList));

}



/**************************************************
*
*    TaxonIdNew()
*
**************************************************/

TaxonIdPtr LIBCALL
TaxonIdNew()
{

	return (TaxonIdPtr) MemNew(sizeof(TaxonId));

}



/**************************************************
*
*    TaxCompleteNew()
*
**************************************************/

TaxCompletePtr LIBCALL
TaxCompleteNew()
{

	return (TaxCompletePtr) MemNew(sizeof(TaxComplete));

}



/**************************************************
*
*    Taxon0ReqAsnRead()
*
**************************************************/
Taxon0ReqPtr LIBCALL Taxon0ReqAsnRead ( AsnIoPtr aip, AsnTypePtr orig)
{
	DataVal av;
	AsnTypePtr atp;
	ValNodePtr vnp;
	Uint1 choice ;
	AsnReadFunc func;

	if (aip == NULL) return NULL;


	if ( ! loaded){
		if ( ! _AsnLoad ())
			goto erret;
		loaded = TRUE;
	}
	if (orig == NULL)
		atp = AsnReadId(aip, amp, TAXON0_REQ);
	else
		atp = AsnLinkType(orig, TAXON0_REQ);
	if(atp == NULL) return NULL;

 	vnp = ValNodeNew(NULL);
 	if (vnp == NULL) goto erret;
 	if ( AsnReadVal(aip, atp, &av) <= 0) goto erret;
 	if ( (atp = AsnReadId(aip, amp, atp)) == NULL) goto erret;
 	func = NULL;
		if (atp == TAXON0_REQ_init){
			choice = Taxon0Req_init;
			AsnReadVal(aip,atp,&av);
			vnp -> data.ptrvalue = av.ptrvalue;
		}
	  		else
		if (atp == TAXON0_REQ_getid){
			choice = Taxon0Req_getid;
			func = (AsnReadFunc) TaxonNameAsnRead;
		}
	  		else
		if (atp == TAXON0_REQ_getref){
			choice = Taxon0Req_getref;
			AsnReadVal(aip,atp,&av);
			vnp -> data.ptrvalue = av.ptrvalue;
		}
	  		else
		if (atp == TAXON0_REQ_getchildren){
			choice = Taxon0Req_getchildren;
			AsnReadVal(aip,atp,&av);
			vnp -> data.ptrvalue = av.ptrvalue;
		}
	  		else
		if (atp == TAXON0_REQ_getparents){
			choice = Taxon0Req_getparents;
			AsnReadVal(aip,atp,&av);
			vnp -> data.ptrvalue = av.ptrvalue;
		}
	  		else
		if (atp == TAXON0_REQ_getgeneticcode){
			choice = Taxon0Req_getgeneticcode;
			AsnReadVal(aip,atp,&av);
			vnp -> data.ptrvalue = av.ptrvalue;
		}
	  		else
		if (atp == TAXON0_REQ_gettaxonline){
			choice = Taxon0Req_gettaxonline;
			AsnReadVal(aip,atp,&av);
			vnp -> data.ptrvalue = av.ptrvalue;
		}
	  		else
		if (atp == TAXON0_REQ_getdivision){
			choice = Taxon0Req_getdivision;
			AsnReadVal(aip,atp,&av);
			vnp -> data.ptrvalue = av.ptrvalue;
		}
	  		else
		if (atp == TAXON0_REQ_getcomplete){
			choice = Taxon0Req_getcomplete;
			func = (AsnReadFunc) TaxonIdNameAsnRead;
		}
	  		else
		if (atp == TAXON0_REQ_fini){
			choice = Taxon0Req_fini;
			AsnReadVal(aip,atp,&av);
			vnp -> data.ptrvalue = av.ptrvalue;
		}
	  	else
 			goto erret;
	  	vnp -> choice = choice;
 		if (func != NULL)
 			vnp -> data.ptrvalue = (*func) (aip,atp);
ret:
	AsnUnlinkType(orig);
	return vnp;
erret:
	vnp=Taxon0ReqFree(vnp);
	goto ret;
}


/**************************************************
*
*    TaxonNameAsnRead()
*
**************************************************/
TaxonNamePtr LIBCALL TaxonNameAsnRead ( AsnIoPtr aip, AsnTypePtr orig)
{
	DataVal av;
	AsnTypePtr atp;
	ValNodePtr vnp;
	Uint1 choice ;
	AsnReadFunc func;

	if (aip == NULL) return NULL;


	if ( ! loaded){
		if ( ! _AsnLoad ())
			goto erret;
		loaded = TRUE;
	}
	if (orig == NULL)
		atp = AsnReadId(aip, amp, TAXON_NAME);
	else
		atp = AsnLinkType(orig, TAXON_NAME);
	if(atp == NULL) return NULL;

 	vnp = ValNodeNew(NULL);
 	if (vnp == NULL) goto erret;
 	if ( AsnReadVal(aip, atp, &av) <= 0) goto erret;
 	if ( (atp = AsnReadId(aip, amp, atp)) == NULL) goto erret;
 	func = NULL;
		if (atp == TAXON_NAME_taxname){
			choice = TaxonName_taxname;
			AsnReadVal(aip,atp,&av);
			vnp -> data.ptrvalue = av.ptrvalue;
		}
	  		else
		if (atp == TAXON_NAME_common){
			choice = TaxonName_common;
			AsnReadVal(aip,atp,&av);
			vnp -> data.ptrvalue = av.ptrvalue;
		}
	  		else
		if (atp == TAXON_NAME_tax_synonym){
			choice = TaxonName_tax_synonym;
			AsnReadVal(aip,atp,&av);
			vnp -> data.ptrvalue = av.ptrvalue;
		}
	  		else
		if (atp == TAXON_NAME_com_synonym){
			choice = TaxonName_com_synonym;
			AsnReadVal(aip,atp,&av);
			vnp -> data.ptrvalue = av.ptrvalue;
		}
	  	else
 			goto erret;
	  	vnp -> choice = choice;
 		if (func != NULL)
 			vnp -> data.ptrvalue = (*func) (aip,atp);
ret:
	AsnUnlinkType(orig);
	return vnp;
erret:
	vnp=TaxonNameFree(vnp);
	goto ret;
}


/**************************************************
*
*    TaxonIdNameAsnRead()
*
**************************************************/
TaxonIdNamePtr LIBCALL TaxonIdNameAsnRead ( AsnIoPtr aip, AsnTypePtr orig)
{
	DataVal av;
	AsnTypePtr atp;
	ValNodePtr vnp;
	Uint1 choice ;
	AsnReadFunc func;

	if (aip == NULL) return NULL;


	if ( ! loaded){
		if ( ! _AsnLoad ())
			goto erret;
		loaded = TRUE;
	}
	if (orig == NULL)
		atp = AsnReadId(aip, amp, TAXON_ID_NAME);
	else
		atp = AsnLinkType(orig, TAXON_ID_NAME);
	if(atp == NULL) return NULL;

 	vnp = ValNodeNew(NULL);
 	if (vnp == NULL) goto erret;
 	if ( AsnReadVal(aip, atp, &av) <= 0) goto erret;
 	if ( (atp = AsnReadId(aip, amp, atp)) == NULL) goto erret;
 	func = NULL;
		if (atp == TAXON_ID_NAME_id){
			choice = TaxonIdName_id;
			AsnReadVal(aip,atp,&av);
			vnp -> data.ptrvalue = av.ptrvalue;
		}
	  		else
		if (atp == TAXON_ID_NAME_name){
			choice = TaxonIdName_name;
			AsnReadVal(aip,atp,&av);
			vnp -> data.ptrvalue = av.ptrvalue;
		}
	  	else
 			goto erret;
	  	vnp -> choice = choice;
 		if (func != NULL)
 			vnp -> data.ptrvalue = (*func) (aip,atp);
ret:
	AsnUnlinkType(orig);
	return vnp;
erret:
	vnp=TaxonIdNameFree(vnp);
	goto ret;
}


/**************************************************
*
*    Taxon0RespAsnRead()
*
**************************************************/
Taxon0RespPtr LIBCALL Taxon0RespAsnRead ( AsnIoPtr aip, AsnTypePtr orig)
{
	DataVal av;
	AsnTypePtr atp;
	ValNodePtr vnp;
	Uint1 choice ;
	AsnReadFunc func;

	if (aip == NULL) return NULL;


	if ( ! loaded){
		if ( ! _AsnLoad ())
			goto erret;
		loaded = TRUE;
	}
	if (orig == NULL)
		atp = AsnReadId(aip, amp, TAXON0_RESP);
	else
		atp = AsnLinkType(orig, TAXON0_RESP);
	if(atp == NULL) return NULL;

 	vnp = ValNodeNew(NULL);
 	if (vnp == NULL) goto erret;
 	if ( AsnReadVal(aip, atp, &av) <= 0) goto erret;
 	if ( (atp = AsnReadId(aip, amp, atp)) == NULL) goto erret;
 	func = NULL;
		if (atp == TAXON0_RESP_error){
			choice = Taxon0Resp_error;
			AsnReadVal(aip,atp,&av);
			vnp -> data.ptrvalue = av.ptrvalue;
		}
	  		else
		if (atp == TAXON0_RESP_init){
			choice = Taxon0Resp_init;
			AsnReadVal(aip,atp,&av);
			vnp -> data.ptrvalue = av.ptrvalue;
		}
	  		else
		if (atp == TAXON0_RESP_getid){
			choice = Taxon0Resp_getid;
			func = (AsnReadFunc) TaxonIdListAsnRead;
		}
	  		else
		if (atp == TAXON0_RESP_getref){
			choice = Taxon0Resp_getref;
			func = (AsnReadFunc) OrgRefAsnRead;
		}
	  		else
		if (atp == TAXON0_RESP_gettaxon){
			choice = Taxon0Resp_gettaxon;
			func = (AsnReadFunc) TaxonIdListAsnRead;
		}
	  		else
		if (atp == TAXON0_RESP_getgeneticcode){
			choice = Taxon0Resp_getgeneticcode;
			func = (AsnReadFunc) GeneticCodeListAsnRead;
		}
	  		else
		if (atp == TAXON0_RESP_gettaxonline){
			choice = Taxon0Resp_gettaxonline;
			AsnReadVal(aip,atp,&av);
			vnp -> data.ptrvalue = av.ptrvalue;
		}
	  		else
		if (atp == TAXON0_RESP_getdivision){
			choice = Taxon0Resp_getdivision;
			AsnReadVal(aip,atp,&av);
			vnp -> data.ptrvalue = av.ptrvalue;
		}
	  		else
		if (atp == TAXON0_RESP_getcomplete){
			choice = Taxon0Resp_getcomplete;
			func = (AsnReadFunc) TaxCompleteListAsnRead;
		}
	  		else
		if (atp == TAXON0_RESP_fini){
			choice = Taxon0Resp_fini;
			AsnReadVal(aip,atp,&av);
			vnp -> data.ptrvalue = av.ptrvalue;
		}
	  	else
 			goto erret;
	  	vnp -> choice = choice;
 		if (func != NULL)
 			vnp -> data.ptrvalue = (*func) (aip,atp);
ret:
	AsnUnlinkType(orig);
	return vnp;
erret:
	vnp=Taxon0RespFree(vnp);
	goto ret;
}


/**************************************************
*
*    TaxonIdListAsnRead()
*
**************************************************/
TaxonIdListPtr LIBCALL TaxonIdListAsnRead ( AsnIoPtr aip, AsnTypePtr orig)
{
	DataVal av;
	AsnTypePtr atp;
	TaxonIdListPtr ptr;

	if (aip == NULL) return NULL;


	if ( ! loaded){
		if ( ! _AsnLoad ())
			goto erret;
		loaded = TRUE;
	}
	if (orig == NULL)
		atp = AsnReadId(aip, amp, TAXON_ID_LIST);
	else
		atp = AsnLinkType(orig, TAXON_ID_LIST);
	if(atp == NULL) return NULL;

	ptr = TaxonIdListNew();
	if (ptr == NULL) goto erret;

	if (AsnReadVal(aip, atp, &av) <= 0) goto erret;

	atp = AsnReadId(aip,amp, atp);
	if (atp == TAXON_ID_LIST_ids)
	{
		AsnTypePtr now_atp;
		ValNodePtr current;
		ValNodePtr head = NULL;
		ValNodePtr prev = NULL;
		Int4 _num_ = 0;
		if (AsnReadVal(aip, atp, &av) <= 0) /* START_STUCT */
			goto erret;

		now_atp = atp;

		while ((now_atp = AsnReadId(aip, amp, now_atp)) != atp){
			if (now_atp == NULL)  goto erret;

			current = ValNodeNew(prev);
			AsnReadVal(aip, now_atp, & current -> data);
			_num_ ++;
			if (current == NULL)  goto erret;

			if (head == NULL)
				head = current;
			else
				prev -> next = current;
			prev=current;
		}
		if (AsnReadVal(aip, atp, &av) <=0) /* END_STRUCT */				goto erret;
		ptr -> _num_ids = _num_;
		ptr -> ids = MemNew( _num_ * sizeof (Pointer) );
		for (_num_ = 0, current = head; current; current = current -> next, _num_ ++){
			(ptr -> ids) [_num_] = current -> data.intvalue;
		}
		ValNodeFree(head);
		atp = AsnReadId(aip,amp, atp);
	}
	if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
ret:
	AsnUnlinkType(orig);
	return ptr;
erret:
	ptr=TaxonIdListFree(ptr);
	goto ret;
}


/**************************************************
*
*    GeneticCodeListAsnRead()
*
**************************************************/
GeneticCodeListPtr LIBCALL GeneticCodeListAsnRead ( AsnIoPtr aip, AsnTypePtr orig)
{
	DataVal av;
	AsnTypePtr atp;
	GeneticCodeListPtr ptr;

	if (aip == NULL) return NULL;


	if ( ! loaded){
		if ( ! _AsnLoad ())
			goto erret;
		loaded = TRUE;
	}
	if (orig == NULL)
		atp = AsnReadId(aip, amp, GENETICCODELIST);
	else
		atp = AsnLinkType(orig, GENETICCODELIST);
	if(atp == NULL) return NULL;

	ptr = GeneticCodeListNew();
	if (ptr == NULL) goto erret;

	if (AsnReadVal(aip, atp, &av) <= 0) goto erret;

	atp = AsnReadId(aip,amp, atp);
	if (atp == GENETICCODELIST_genomic){
		if (AsnReadVal(aip, atp, &av) <= 0)
			goto erret;
		ptr->genomic = av.intvalue;
		atp = AsnReadId(aip,amp, atp);
	}
	if (atp == GENETICCODELIST_mitochondrial){
		if (AsnReadVal(aip, atp, &av) <= 0)
			goto erret;
		ptr->mitochondrial = av.intvalue;
		atp = AsnReadId(aip,amp, atp);
	}
	if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
ret:
	AsnUnlinkType(orig);
	return ptr;
erret:
	ptr=GeneticCodeListFree(ptr);
	goto ret;
}


/**************************************************
*
*    TaxCompleteListAsnRead()
*
**************************************************/
TaxCompleteListPtr LIBCALL TaxCompleteListAsnRead ( AsnIoPtr aip, AsnTypePtr orig)
{
	DataVal av;
	AsnTypePtr atp;
	TaxCompleteListPtr ptr;

	if (aip == NULL) return NULL;


	if ( ! loaded){
		if ( ! _AsnLoad ())
			goto erret;
		loaded = TRUE;
	}
	if (orig == NULL)
		atp = AsnReadId(aip, amp, TAX_COMPLETE_LIST);
	else
		atp = AsnLinkType(orig, TAX_COMPLETE_LIST);
	if(atp == NULL) return NULL;

	ptr = TaxCompleteListNew();
	if (ptr == NULL) goto erret;

	if (AsnReadVal(aip, atp, &av) <= 0) goto erret;

	atp = AsnReadId(aip,amp, atp);
	if (atp == TAX_COMPLETE_LIST_num){
		if (AsnReadVal(aip, atp, &av) <= 0)
			goto erret;
		ptr->num = av.intvalue;
		atp = AsnReadId(aip,amp, atp);
	}
	if (atp == TAX_COMPLETE_LIST_info)
	{
		AsnTypePtr now_atp;
		TaxCompletePtr current;
		TaxCompletePtr head = NULL;
		TaxCompletePtr prev = NULL;
		if (AsnReadVal(aip, atp, &av) <= 0) /* START_STUCT */
			goto erret;

		now_atp = atp;

		while ((now_atp = AsnReadId(aip, amp, now_atp)) != atp){
			if (now_atp == NULL)  goto erret;

			current= TaxCompleteAsnRead(aip, now_atp);
			if (current == NULL)  goto erret;

			if (head == NULL)
				head = current;
			else
				prev -> next = current;
			prev=current;
		}
		if (AsnReadVal(aip, atp, &av) <=0) /* END_STRUCT */				goto erret;
		ptr -> info = head;
		atp = AsnReadId(aip,amp, atp);
	}
	if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
ret:
	AsnUnlinkType(orig);
	return ptr;
erret:
	ptr=TaxCompleteListFree(ptr);
	goto ret;
}


/**************************************************
*
*    TaxonIdAsnRead()
*
**************************************************/
TaxonIdPtr LIBCALL TaxonIdAsnRead ( AsnIoPtr aip, AsnTypePtr orig)
{
	DataVal av;
	AsnTypePtr atp;
	TaxonIdPtr ptr;

	if (aip == NULL) return NULL;


	if ( ! loaded){
		if ( ! _AsnLoad ())
			goto erret;
		loaded = TRUE;
	}
	if (orig == NULL)
		atp = AsnReadId(aip, amp, TAXON_ID);
	else
		atp = AsnLinkType(orig, TAXON_ID);
	if(atp == NULL) return NULL;

	ptr = TaxonIdNew();
	if (ptr == NULL) goto erret;

	if (AsnReadVal(aip, atp, &av) <= 0) goto erret;

	atp = AsnReadId(aip,amp, atp);
	if (atp == TAXON_ID_id){
		if (AsnReadVal(aip, atp, &av) <= 0)
			goto erret;
		ptr->id = av.intvalue;
		atp = AsnReadId(aip,amp, atp);
	}
	if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
ret:
	AsnUnlinkType(orig);
	return ptr;
erret:
	ptr=TaxonIdFree(ptr);
	goto ret;
}


/**************************************************
*
*    TaxCompleteAsnRead()
*
**************************************************/
TaxCompletePtr LIBCALL TaxCompleteAsnRead ( AsnIoPtr aip, AsnTypePtr orig)
{
	DataVal av;
	AsnTypePtr atp;
	TaxCompletePtr ptr;

	if (aip == NULL) return NULL;


	if ( ! loaded){
		if ( ! _AsnLoad ())
			goto erret;
		loaded = TRUE;
	}
	if (orig == NULL)
		atp = AsnReadId(aip, amp, TAX_COMPLETE);
	else
		atp = AsnLinkType(orig, TAX_COMPLETE);
	if(atp == NULL) return NULL;

	ptr = TaxCompleteNew();
	if (ptr == NULL) goto erret;

	if (AsnReadVal(aip, atp, &av) <= 0) goto erret;

	atp = AsnReadId(aip,amp, atp);
	if (atp == TAX_COMPLETE_sciname){
		if (AsnReadVal(aip, atp, &av) <= 0)
			goto erret;
		ptr->sciname = (CharPtr)av.ptrvalue;
		atp = AsnReadId(aip,amp, atp);
	}
	if (atp == TAX_COMPLETE_comname){
		if (AsnReadVal(aip, atp, &av) <= 0)
			goto erret;
		ptr->comname = (CharPtr)av.ptrvalue;
		atp = AsnReadId(aip,amp, atp);
	}
	if (atp == TAX_COMPLETE_synonyms){
		if (AsnReadVal(aip, atp, &av) <= 0)
			goto erret;
		ptr->synonyms = (CharPtr)av.ptrvalue;
		atp = AsnReadId(aip,amp, atp);
	}
	if (atp == TAX_COMPLETE_id_gc){
		if (AsnReadVal(aip, atp, &av) <= 0)
			goto erret;
		ptr->id_gc = av.intvalue;
		atp = AsnReadId(aip,amp, atp);
	}
	if (atp == TAX_COMPLETE_name_gc){
		if (AsnReadVal(aip, atp, &av) <= 0)
			goto erret;
		ptr->name_gc = (CharPtr)av.ptrvalue;
		atp = AsnReadId(aip,amp, atp);
	}
	if (atp == TAX_COMPLETE_id_mgc){
		if (AsnReadVal(aip, atp, &av) <= 0)
			goto erret;
		ptr->id_mgc = av.intvalue;
		atp = AsnReadId(aip,amp, atp);
	}
	if (atp == TAX_COMPLETE_name_mgc){
		if (AsnReadVal(aip, atp, &av) <= 0)
			goto erret;
		ptr->name_mgc = (CharPtr)av.ptrvalue;
		atp = AsnReadId(aip,amp, atp);
	}
	if (atp == TAX_COMPLETE_gb_div){
		if (AsnReadVal(aip, atp, &av) <= 0)
			goto erret;
		ptr->gb_div = (CharPtr)av.ptrvalue;
		atp = AsnReadId(aip,amp, atp);
	}
	if (atp == TAX_COMPLETE_embl_code){
		if (AsnReadVal(aip, atp, &av) <= 0)
			goto erret;
		ptr->embl_code = (CharPtr)av.ptrvalue;
		atp = AsnReadId(aip,amp, atp);
	}
	if (atp == TAX_COMPLETE_lineage){
		if (AsnReadVal(aip, atp, &av) <= 0)
			goto erret;
		ptr->lineage = (CharPtr)av.ptrvalue;
		atp = AsnReadId(aip,amp, atp);
	}
	if (atp == TAX_COMPLETE_is_species_level){
		if (AsnReadVal(aip, atp, &av) <= 0)
			goto erret;
		ptr->is_species_level = (Uint1)av.boolvalue;
		atp = AsnReadId(aip,amp, atp);
	}
	if (AsnReadVal(aip, atp, &av) <= 0) goto erret;
ret:
	AsnUnlinkType(orig);
	return ptr;
erret:
	ptr=TaxCompleteFree(ptr);
	goto ret;
}


/**************************************************
*
*    Taxon0ReqAsnWrite()
*
**************************************************/
Boolean LIBCALL Taxon0ReqAsnWrite (Taxon0ReqPtr ptr,  AsnIoPtr aip, AsnTypePtr orig)
{
	Boolean retval = FALSE;
	DataVal av;
	AsnTypePtr atp;
	AsnWriteFunc func;
	AsnTypePtr writetype = NULL;
	Pointer pnt;

	if (aip == NULL) return FALSE;


	if ( ! loaded){
		if ( ! _AsnLoad ())
			goto erret;
		amp = AsnAllModPtr();
		loaded = TRUE;
	}
	atp = AsnLinkType(orig, TAXON0_REQ);
	if(atp == NULL) {AsnNullValueMsg(aip,atp); goto erret;}
 	av.ptrvalue = (Pointer) ptr;
 	if (! AsnWriteChoice(aip, atp, ptr -> choice, & av)) goto erret;
		pnt = ptr -> data.ptrvalue;
		av.intvalue = ptr -> data.intvalue;
 	switch (ptr -> choice){
		case  Taxon0Req_init :
		retval = AsnWrite(aip, TAXON0_REQ_init, &av); break;
		case  Taxon0Req_getid :
		writetype = TAXON0_REQ_getid;
		 func = (AsnWriteFunc) TaxonNameAsnWrite;
			break;
		case  Taxon0Req_getref :
		retval = AsnWrite(aip, TAXON0_REQ_getref, &av); break;
		case  Taxon0Req_getchildren :
		retval = AsnWrite(aip, TAXON0_REQ_getchildren, &av); break;
		case  Taxon0Req_getparents :
		retval = AsnWrite(aip, TAXON0_REQ_getparents, &av); break;
		case  Taxon0Req_getgeneticcode :
		retval = AsnWrite(aip, TAXON0_REQ_getgeneticcode, &av); break;
		case  Taxon0Req_gettaxonline :
		retval = AsnWrite(aip, TAXON0_REQ_gettaxonline, &av); break;
		case  Taxon0Req_getdivision :
		retval = AsnWrite(aip, TAXON0_REQ_getdivision, &av); break;
		case  Taxon0Req_getcomplete :
		writetype = TAXON0_REQ_getcomplete;
		 func = (AsnWriteFunc) TaxonIdNameAsnWrite;
			break;
		case  Taxon0Req_fini :
		retval = AsnWrite(aip, TAXON0_REQ_fini, &av); break;
		}
 	if (writetype != NULL)
 			retval = (*func) (pnt, aip,writetype);
ret:
	AsnUnlinkType(orig);
	return retval;
erret:
	ptr=Taxon0ReqFree(ptr);
	goto ret;
}


/**************************************************
*
*    TaxonNameAsnWrite()
*
**************************************************/
Boolean LIBCALL TaxonNameAsnWrite (TaxonNamePtr ptr,  AsnIoPtr aip, AsnTypePtr orig)
{
	Boolean retval = FALSE;
	DataVal av;
	AsnTypePtr atp;
	AsnWriteFunc func;
	AsnTypePtr writetype = NULL;
	Pointer pnt;

	if (aip == NULL) return FALSE;


	if ( ! loaded){
		if ( ! _AsnLoad ())
			goto erret;
		amp = AsnAllModPtr();
		loaded = TRUE;
	}
	atp = AsnLinkType(orig, TAXON_NAME);
	if(atp == NULL) {AsnNullValueMsg(aip,atp); goto erret;}
 	av.ptrvalue = (Pointer) ptr;
 	if (! AsnWriteChoice(aip, atp, ptr -> choice, & av)) goto erret;
		pnt = ptr -> data.ptrvalue;
		av.intvalue = ptr -> data.intvalue;
 	switch (ptr -> choice){
		case  TaxonName_taxname :
		retval = AsnWrite(aip, TAXON_NAME_taxname, &av); break;
		case  TaxonName_common :
		retval = AsnWrite(aip, TAXON_NAME_common, &av); break;
		case  TaxonName_tax_synonym :
		retval = AsnWrite(aip, TAXON_NAME_tax_synonym, &av); break;
		case  TaxonName_com_synonym :
		retval = AsnWrite(aip, TAXON_NAME_com_synonym, &av); break;
		}
 	if (writetype != NULL)
 			retval = (*func) (pnt, aip,writetype);
ret:
	AsnUnlinkType(orig);
	return retval;
erret:
	ptr=TaxonNameFree(ptr);
	goto ret;
}


/**************************************************
*
*    TaxonIdNameAsnWrite()
*
**************************************************/
Boolean LIBCALL TaxonIdNameAsnWrite (TaxonIdNamePtr ptr,  AsnIoPtr aip, AsnTypePtr orig)
{
	Boolean retval = FALSE;
	DataVal av;
	AsnTypePtr atp;
	AsnWriteFunc func;
	AsnTypePtr writetype = NULL;
	Pointer pnt;

	if (aip == NULL) return FALSE;


	if ( ! loaded){
		if ( ! _AsnLoad ())
			goto erret;
		amp = AsnAllModPtr();
		loaded = TRUE;
	}
	atp = AsnLinkType(orig, TAXON_ID_NAME);
	if(atp == NULL) {AsnNullValueMsg(aip,atp); goto erret;}
 	av.ptrvalue = (Pointer) ptr;
 	if (! AsnWriteChoice(aip, atp, ptr -> choice, & av)) goto erret;
		pnt = ptr -> data.ptrvalue;
		av.intvalue = ptr -> data.intvalue;
 	switch (ptr -> choice){
		case  TaxonIdName_id :
		retval = AsnWrite(aip, TAXON_ID_NAME_id, &av); break;
		case  TaxonIdName_name :
		retval = AsnWrite(aip, TAXON_ID_NAME_name, &av); break;
		}
 	if (writetype != NULL)
 			retval = (*func) (pnt, aip,writetype);
ret:
	AsnUnlinkType(orig);
	return retval;
erret:
	ptr=TaxonIdNameFree(ptr);
	goto ret;
}


/**************************************************
*
*    Taxon0RespAsnWrite()
*
**************************************************/
Boolean LIBCALL Taxon0RespAsnWrite (Taxon0RespPtr ptr,  AsnIoPtr aip, AsnTypePtr orig)
{
	Boolean retval = FALSE;
	DataVal av;
	AsnTypePtr atp;
	AsnWriteFunc func;
	AsnTypePtr writetype = NULL;
	Pointer pnt;

	if (aip == NULL) return FALSE;


	if ( ! loaded){
		if ( ! _AsnLoad ())
			goto erret;
		amp = AsnAllModPtr();
		loaded = TRUE;
	}
	atp = AsnLinkType(orig, TAXON0_RESP);
	if(atp == NULL) {AsnNullValueMsg(aip,atp); goto erret;}
 	av.ptrvalue = (Pointer) ptr;
 	if (! AsnWriteChoice(aip, atp, ptr -> choice, & av)) goto erret;
		pnt = ptr -> data.ptrvalue;
		av.intvalue = ptr -> data.intvalue;
 	switch (ptr -> choice){
		case  Taxon0Resp_error :
		retval = AsnWrite(aip, TAXON0_RESP_error, &av); break;
		case  Taxon0Resp_init :
		retval = AsnWrite(aip, TAXON0_RESP_init, &av); break;
		case  Taxon0Resp_getid :
		writetype = TAXON0_RESP_getid;
		 func = (AsnWriteFunc) TaxonIdListAsnWrite;
			break;
		case  Taxon0Resp_getref :
		writetype = TAXON0_RESP_getref;
		 func = (AsnWriteFunc) OrgRefAsnWrite;
			break;
		case  Taxon0Resp_gettaxon :
		writetype = TAXON0_RESP_gettaxon;
		 func = (AsnWriteFunc) TaxonIdListAsnWrite;
			break;
		case  Taxon0Resp_getgeneticcode :
		writetype = TAXON0_RESP_getgeneticcode;
		 func = (AsnWriteFunc) GeneticCodeListAsnWrite;
			break;
		case  Taxon0Resp_gettaxonline :
		retval = AsnWrite(aip, TAXON0_RESP_gettaxonline, &av); break;
		case  Taxon0Resp_getdivision :
		retval = AsnWrite(aip, TAXON0_RESP_getdivision, &av); break;
		case  Taxon0Resp_getcomplete :
		writetype = TAXON0_RESP_getcomplete;
		 func = (AsnWriteFunc) TaxCompleteListAsnWrite;
			break;
		case  Taxon0Resp_fini :
		retval = AsnWrite(aip, TAXON0_RESP_fini, &av); break;
		}
 	if (writetype != NULL)
 			retval = (*func) (pnt, aip,writetype);
ret:
	AsnUnlinkType(orig);
	return retval;
erret:
	ptr=Taxon0RespFree(ptr);
	goto ret;
}


/**************************************************
*
*    TaxonIdListAsnWrite()
*
**************************************************/
Boolean LIBCALL TaxonIdListAsnWrite (TaxonIdListPtr ptr,  AsnIoPtr aip, AsnTypePtr orig)
{
	Boolean retval = FALSE;
	DataVal av;
	AsnTypePtr atp;

	if (aip == NULL) return FALSE;


	if ( ! loaded){
		if ( ! _AsnLoad ())
			goto erret;
		amp = AsnAllModPtr();
		loaded = TRUE;
	}
	atp = AsnLinkType(orig, TAXON_ID_LIST);
	if(atp == NULL) {AsnNullValueMsg(aip,atp); goto erret;}
	if (! AsnOpenStruct (aip, atp, (Pointer) ptr)) goto erret;

	if (ptr -> ids != NULL){
	if (! AsnOpenStruct (aip, TAXON_ID_LIST_ids, (Pointer) ptr -> ids))
			goto erret;

	{Int4   _num_ ;
		for( _num_ = 0; _num_ < ptr -> _num_ids; _num_ ++){
			av.intvalue = (ptr -> ids) [_num_];
			if ( ! AsnWrite(aip, TAXON_ID_LIST_ids_E, & av)) goto erret;
		}}

	if (! AsnCloseStruct (aip, TAXON_ID_LIST_ids, (Pointer) ptr -> ids))
			goto erret;
	}
	if (! AsnCloseStruct (aip, atp, (Pointer) ptr)) goto erret;
	retval = TRUE;
erret:AsnUnlinkType(orig);
	return retval;
}


/**************************************************
*
*    GeneticCodeListAsnWrite()
*
**************************************************/
Boolean LIBCALL GeneticCodeListAsnWrite (GeneticCodeListPtr ptr,  AsnIoPtr aip, AsnTypePtr orig)
{
	Boolean retval = FALSE;
	DataVal av;
	AsnTypePtr atp;

	if (aip == NULL) return FALSE;


	if ( ! loaded){
		if ( ! _AsnLoad ())
			goto erret;
		amp = AsnAllModPtr();
		loaded = TRUE;
	}
	atp = AsnLinkType(orig, GENETICCODELIST);
	if(atp == NULL) {AsnNullValueMsg(aip,atp); goto erret;}
	if (! AsnOpenStruct (aip, atp, (Pointer) ptr)) goto erret;

		av.intvalue = (int) ptr -> genomic;
		if ( ! AsnWrite(aip,GENETICCODELIST_genomic, &av)) goto erret;
		av.intvalue = (int) ptr -> mitochondrial;
		if ( ! AsnWrite(aip,GENETICCODELIST_mitochondrial, &av)) goto erret;
	if (! AsnCloseStruct (aip, atp, (Pointer) ptr)) goto erret;
	retval = TRUE;
erret:AsnUnlinkType(orig);
	return retval;
}


/**************************************************
*
*    TaxCompleteListAsnWrite()
*
**************************************************/
Boolean LIBCALL TaxCompleteListAsnWrite (TaxCompleteListPtr ptr,  AsnIoPtr aip, AsnTypePtr orig)
{
	Boolean retval = FALSE;
	DataVal av;
	AsnTypePtr atp;

	if (aip == NULL) return FALSE;


	if ( ! loaded){
		if ( ! _AsnLoad ())
			goto erret;
		amp = AsnAllModPtr();
		loaded = TRUE;
	}
	atp = AsnLinkType(orig, TAX_COMPLETE_LIST);
	if(atp == NULL) {AsnNullValueMsg(aip,atp); goto erret;}
	if (! AsnOpenStruct (aip, atp, (Pointer) ptr)) goto erret;

		av.intvalue = (int) ptr -> num;
		if ( ! AsnWrite(aip,TAX_COMPLETE_LIST_num, &av)) goto erret;
	if (ptr -> info != NULL){
	if (! AsnOpenStruct (aip, TAX_COMPLETE_LIST_info, (Pointer) ptr -> info))
			goto erret;

		{TaxCompletePtr current;
		for(current=ptr -> info; current; current = current -> next){
		if ( ! TaxCompleteAsnWrite(current,aip,TAX_COMPLETE_LIST_info_E)) goto erret;
		}}

	if (! AsnCloseStruct (aip, TAX_COMPLETE_LIST_info, (Pointer) ptr -> info))
			goto erret;
	}
	if (! AsnCloseStruct (aip, atp, (Pointer) ptr)) goto erret;
	retval = TRUE;
erret:AsnUnlinkType(orig);
	return retval;
}


/**************************************************
*
*    TaxonIdAsnWrite()
*
**************************************************/
Boolean LIBCALL TaxonIdAsnWrite (TaxonIdPtr ptr,  AsnIoPtr aip, AsnTypePtr orig)
{
	Boolean retval = FALSE;
	DataVal av;
	AsnTypePtr atp;

	if (aip == NULL) return FALSE;


	if ( ! loaded){
		if ( ! _AsnLoad ())
			goto erret;
		amp = AsnAllModPtr();
		loaded = TRUE;
	}
	atp = AsnLinkType(orig, TAXON_ID);
	if(atp == NULL) {AsnNullValueMsg(aip,atp); goto erret;}
	if (! AsnOpenStruct (aip, atp, (Pointer) ptr)) goto erret;

		av.intvalue = (int) ptr -> id;
		if ( ! AsnWrite(aip,TAXON_ID_id, &av)) goto erret;
	if (! AsnCloseStruct (aip, atp, (Pointer) ptr)) goto erret;
	retval = TRUE;
erret:AsnUnlinkType(orig);
	return retval;
}


/**************************************************
*
*    TaxCompleteAsnWrite()
*
**************************************************/
Boolean LIBCALL TaxCompleteAsnWrite (TaxCompletePtr ptr,  AsnIoPtr aip, AsnTypePtr orig)
{
	Boolean retval = FALSE;
	DataVal av;
	AsnTypePtr atp;

	if (aip == NULL) return FALSE;


	if ( ! loaded){
		if ( ! _AsnLoad ())
			goto erret;
		amp = AsnAllModPtr();
		loaded = TRUE;
	}
	atp = AsnLinkType(orig, TAX_COMPLETE);
	if(atp == NULL) {AsnNullValueMsg(aip,atp); goto erret;}
	if (! AsnOpenStruct (aip, atp, (Pointer) ptr)) goto erret;

	if (ptr -> sciname != NULL){
		av.ptrvalue = (Pointer) ptr -> sciname;
		if ( ! AsnWrite(aip,TAX_COMPLETE_sciname, &av)) goto erret;
	}
	if (ptr -> comname != NULL){
		av.ptrvalue = (Pointer) ptr -> comname;
		if ( ! AsnWrite(aip,TAX_COMPLETE_comname, &av)) goto erret;
	}
	if (ptr -> synonyms != NULL){
		av.ptrvalue = (Pointer) ptr -> synonyms;
		if ( ! AsnWrite(aip,TAX_COMPLETE_synonyms, &av)) goto erret;
	}
		av.intvalue = (int) ptr -> id_gc;
		if ( ! AsnWrite(aip,TAX_COMPLETE_id_gc, &av)) goto erret;
	if (ptr -> name_gc != NULL){
		av.ptrvalue = (Pointer) ptr -> name_gc;
		if ( ! AsnWrite(aip,TAX_COMPLETE_name_gc, &av)) goto erret;
	}
		av.intvalue = (int) ptr -> id_mgc;
		if ( ! AsnWrite(aip,TAX_COMPLETE_id_mgc, &av)) goto erret;
	if (ptr -> name_mgc != NULL){
		av.ptrvalue = (Pointer) ptr -> name_mgc;
		if ( ! AsnWrite(aip,TAX_COMPLETE_name_mgc, &av)) goto erret;
	}
	if (ptr -> gb_div != NULL){
		av.ptrvalue = (Pointer) ptr -> gb_div;
		if ( ! AsnWrite(aip,TAX_COMPLETE_gb_div, &av)) goto erret;
	}
	if (ptr -> embl_code != NULL){
		av.ptrvalue = (Pointer) ptr -> embl_code;
		if ( ! AsnWrite(aip,TAX_COMPLETE_embl_code, &av)) goto erret;
	}
	if (ptr -> lineage != NULL){
		av.ptrvalue = (Pointer) ptr -> lineage;
		if ( ! AsnWrite(aip,TAX_COMPLETE_lineage, &av)) goto erret;
	}
		av.boolvalue = (int) ptr -> is_species_level;
		if ( ! AsnWrite(aip,TAX_COMPLETE_is_species_level, &av)) goto erret;
	if (! AsnCloseStruct (aip, atp, (Pointer) ptr)) goto erret;
	retval = TRUE;
erret:AsnUnlinkType(orig);
	return retval;
}
